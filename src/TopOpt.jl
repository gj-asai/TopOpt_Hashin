module TopOpt

using Ferrite
using Ferrite: nnodes_per_cell
using Tensors, SparseArrays, LinearAlgebra
using ForwardDiff

using Nonconvex, NonconvexMMA
Nonconvex.@load MMA

export StressLimits, Isotropic2D, Orthotropic2D
export NodalLoad, LinearLoad
export FEModel, OptimOpts
export Compliance, MatrixFailure, FibreFailure
export get_centers
export topopt

struct StressLimits
    Xt
    Xc
    Yt
    Yc
    Sln
    Stn
end
StressLimits() = StressLimits(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

struct Material{dim}
    C::SymmetricTensor{4,dim}
    stress_limits::StressLimits
end

function Isotropic2D(; E, ν, limits=StressLimits())
    S = [1/E -ν/E 0;
        -ν/E 1/E 0;
        0 0 2(1+ν)/E]
    C = fromvoigt(SymmetricTensor{4,2}, inv(S))
    return Material(C, limits)
end

function Orthotropic2D(; El, Et, νlt, Glt, limits=StressLimits())
    S = [1/El -νlt/El 0;
        -νlt/El 1/Et 0;
        0 0 1/Glt]
    C = fromvoigt(SymmetricTensor{4,2}, inv(S))
    return Material(C, limits)
end

abstract type Load{dim} end
struct NodalLoad{dim} <: Load{dim}
    nodeset_name::AbstractString
    F::NTuple{dim}
end
struct LinearLoad{dim} <: Load{dim}
    faceset_name::AbstractString
    F::NTuple{dim}
end

struct FEModel{dim}
    mat::Material{dim}
    grid::Grid
    elemvol::Vector{<:Real}

    constraints::Vector{Dirichlet}
    loads::Vector{<:Load{dim}}

    cellvalues::CellValues
    facevalues::FaceValues
    dh::DofHandler
    ch::ConstraintHandler
end

function FEModel(; mat::Material{dim}, grid::Grid, ip::Interpolation, qr::QuadratureRule{dim,shape}, constraints::Vector{Dirichlet}, loads::Vector{<:Load}) where {dim,shape}
    # element type and quadrature rule
    cellvalues = CellVectorValues(qr, ip)

    qr_order = length(getpoints(qr))
    face_qr = QuadratureRule{dim - 1,shape}(qr_order)
    facevalues = FaceVectorValues(face_qr, ip)

    # degrees of freedom
    dh = DofHandler(grid)
    add!(dh, :u, dim)
    close!(dh)

    # constraints
    ch = ConstraintHandler(dh)
    for cc in constraints
        add!(ch, cc)
    end
    close!(ch)

    # elemental volume vector
    elemvol = zeros(getncells(grid))
    for cell in CellIterator(dh)
        reinit!(cellvalues, cell)
        for q_point in 1:getnquadpoints(cellvalues)
            elemvol[cellid(cell)] += getdetJdV(cellvalues, q_point)
        end
    end

    return FEModel(mat, grid, elemvol, constraints, loads, cellvalues, facevalues, dh, ch)
end

abstract type Objective end
struct Compliance <: Objective end
struct MatrixFailure <: Objective end
struct FibreFailure <: Objective end

@kwdef struct OptimOpts
    maxiter::Integer
    volfrac::Real
    rρ::Real
    rθ::Real
    reltol::Real
    objective::Objective
end

function topopt(model::FEModel{dim}, opts::OptimOpts) where {dim}
    Hρ = convolution_filter(opts.rρ, model)
    Hθ = convolution_filter(opts.rθ, model)
    f = global_force(model)
    penal = 3

    # preallocate stiffness and sensitivities
    n_basefuncs = getnbasefunctions(model.cellvalues)
    K = create_sparsity_pattern(model.dh)
    ∂Ke∂x = fill(zeros(n_basefuncs, n_basefuncs), 2 * getncells(model.grid))
    Keinv = fill(zeros(n_basefuncs, n_basefuncs), getncells(model.grid))

    # preallocate solution
    u = zeros(dim * getnnodes(model.grid))

    # preallocate failure indices and compliance
    c = [0.0] # vector only to be used in local scopes, should always be called as c[1]
    IFm = zeros(getnnodes(model.grid))
    IFf = zeros(getnnodes(model.grid))

    # preallocate sensitivities
    dcdx = zeros(2 * getncells(model.grid))
    dIFmdx = zeros(2 * getncells(model.grid))
    dIFfdx = zeros(2 * getncells(model.grid))

    # optimisation history
    c_hist = Float64[]
    IF_hist = Float64[]

    function fea(x::AbstractVector)
        # filter orientations
        x[2:2:end] .= Hθ * x[2:2:end]

        # assemble linear system
        global_stiffness!(K, ∂Ke∂x, x, penal, model)
        apply!(K, f, model.ch)

        # get stiffness matrix slices corresponding to each element and invert
        for cell in CellIterator(model.dh)
            Keinv[cellid(cell)] = K[celldofs(cell), celldofs(cell)] |> collect |> inv
        end

        # solve linear system
        u .= K \ f
        c[1] = u' * K * u

        # compliance sensitivity
        for cell in CellIterator(model.dh)
            ue = u[celldofs(cell)]
            dcdx[2*cellid(cell)-1] = -ue' * ∂Ke∂x[2*cellid(cell)-1] * ue
            dcdx[2*cellid(cell)] = -ue' * ∂Ke∂x[2*cellid(cell)] * ue
        end

        # filter density
        dcdx[1:2:end] .= Hρ * (x[1:2:end] .* dcdx[1:2:end]) ./ x[1:2:end]

        ρ = x[1:2:end]
        θ = x[2:2:end]

        # compute failure indices and sensitivities
        dIFmdx = zeros(2 * getncells(model.grid))
        dIFfdx = zeros(2 * getncells(model.grid))
        for cell in CellIterator(model.dh)
            e = cellid(cell)
            ue = u[celldofs(cell)]
            xe_ue = [ρ[e], θ[e], ue...]

            # compute stresses
            σ = stress(xe_ue, penal, e, model)
            σ = Tensor{2,dim}(σ)

            # using ForwardDiff's AD because Tensors's AD does not support using Vec's of arbitrary dimension
            jac = ForwardDiff.jacobian(xe_ue -> stress(xe_ue, penal, e, model), xe_ue)
            ∂σ∂ρ = reshape(jac[:, 1], (dim, dim)) |> Tensor{2,dim}
            ∂σ∂θ = reshape(jac[:, 2], (dim, dim)) |> Tensor{2,dim}
            ∂σ∂u = [reshape(d, (dim, dim)) |> Tensor{2,dim} for d in eachcol(jac[:, 3:end])]

            ∂u∂ρ = -Keinv[e] * ∂Ke∂x[2*e-1] * ue
            ∂u∂θ = -Keinv[e] * ∂Ke∂x[2*e] * ue

            # chain rule
            dσdρ = ∂σ∂ρ + sum(∂σ∂u .* ∂u∂ρ)
            dσdθ = ∂σ∂θ + sum(∂σ∂u .* ∂u∂θ)

            # compute failure indices
            ∂IFm∂σ, IFm[e] = gradient(σ -> hashin_matrix(σ, model.mat), σ, :all)
            ∂IFf∂σ, IFf[e] = gradient(σ -> hashin_fibre(σ, model.mat), σ, :all)

            # chain rule again
            dIFmdx[2*e-1] = ∂IFm∂σ ⊡ dσdρ
            dIFmdx[2*e] = ∂IFm∂σ ⊡ dσdθ
            dIFfdx[2*e-1] = ∂IFf∂σ ⊡ dσdρ
            dIFfdx[2*e] = ∂IFf∂σ ⊡ dσdθ
        end
    end

    function obj(x::AbstractVector)
        fea(x)
        opts.objective isa Compliance && return c[1]
        opts.objective isa MatrixFailure && return -IFm
        opts.objective isa FibreFailure && return -IFf
    end
    function dobj(x::AbstractVector)
        opts.objective isa Compliance && return dcdx
        opts.objective isa MatrixFailure && return -dIFmdx
        opts.objective isa FibreFailure && return -dIFfdx
    end

    # volume fraction constraint function/sensitivities
    function constraint(x)
        return x[1:2:end] ⋅ model.elemvol / sum(model.elemvol) - opts.volfrac
    end
    function dconstraint(x)
        ∂g∂x = zero(x)
        ∂g∂x[1:2:end] .= model.elemvol / sum(model.elemvol)
        return ∂g∂x
    end

    function post(solution; update=false)
        push!(c_hist, c[1])
        push!(IF_hist, max(IFm..., IFf...))
    end

    # create optimisation model
    optim = Model(CustomGradFunction(obj, dobj))
    for i = 1:getncells(model.grid)
        addvar!(optim, 1e-3, 1) # ρ
        addvar!(optim, -pi, pi) # θ
    end
    add_ineq_constraint!(optim, CustomGradFunction(constraint, dconstraint))

    # initial values
    x = zeros(2 * getncells(model.grid))
    x[1:2:end] .= opts.volfrac
    x[2:2:end] .= deg2rad(0)

    while penal <= 3
        r = optimize(
            optim,
            MMA(),
            x,
            options=MMAOptions(
                maxiter=opts.maxiter,
                convcriteria=GenericCriteria(),
                tol=Tolerance(x=0.0, fabs=0.0, frel=opts.reltol),
            ),
            callback=post,
        )

        x .= r.minimizer
        penal += 1
    end

    ρ = x[1:2:end]
    θ = x[2:2:end]

    return ρ, θ, c_hist, IF_hist, IFm, IFf
end

function global_stiffness!(K, ∂Ke∂x, x::AbstractVector, p::Real, model::FEModel)
    assembler = start_assemble(K)

    # preallocate elemental stiffness
    n_basefuncs = getnbasefunctions(model.cellvalues)
    Ke = zeros(n_basefuncs, n_basefuncs)

    for cell in CellIterator(model.dh)
        reinit!(model.cellvalues, cell)

        e = cellid(cell)
        xe = [x[2*e-1], x[2*e]]

        # compute Ke and derivatives
        jac = ForwardDiff.jacobian((Ke, xe) -> element_stiffness!(Ke, xe, p, model), Ke, xe)

        ∂Ke∂x[2*e-1] = reshape(jac[:, 1], (n_basefuncs, n_basefuncs))
        ∂Ke∂x[2*e] = reshape(jac[:, 2], (n_basefuncs, n_basefuncs))

        assemble!(assembler, celldofs(cell), Ke)
    end
end

function element_stiffness!(Ke::Matrix{T}, xe::Vector{T}, p::Real, model::FEModel) where {T<:Real}
    fill!(Ke, zero(T))

    ρe = xe[1]
    θe = xe[2]

    cellvalues = model.cellvalues
    @inbounds for q_point in 1:getnquadpoints(cellvalues)
        dΩ = getdetJdV(cellvalues, q_point)
        for i in 1:getnbasefunctions(cellvalues)
            δεi = shape_symmetric_gradient(cellvalues, q_point, i)
            for j in 1:i
                δεj = shape_symmetric_gradient(cellvalues, q_point, j)
                Ke[i, j] += ρe^p * (δεi ⊡ rotate(model.mat.C, θe) ⊡ δεj) * dΩ
            end
        end
    end

    Ke .= Symmetric(Ke, :L)
end

function global_force(model::FEModel{dim}) where {dim}
    f = zeros(ndofs(model.dh))

    # nodal forces
    for cell in CellIterator(model.dh)
        dofs = celldofs(cell)
        for (i, node) in enumerate(getnodes(cell))
            for force in model.loads
                force isa NodalLoad || continue
                node in getnodeset(model.grid, force.nodeset_name) || continue
                f[dofs[dim*(i-1)+1:dim*i]] .+= force.F
            end
        end
    end

    # linear forces
    for cell in CellIterator(model.dh)
        dofs = celldofs(cell)
        for face in 1:nfaces(cell)
            for force in model.loads
                force isa LinearLoad || continue
                (cellid(cell), face) in getfaceset(model.grid, force.faceset_name) || continue

                reinit!(model.facevalues, cell, face)
                for q_point in 1:getnquadpoints(model.facevalues)
                    dΓ = getdetJdV(model.facevalues, q_point)
                    for i in 1:getnbasefunctions(model.facevalues)
                        δu = shape_value(model.facevalues, q_point, i)
                        f[dofs[i]] += (δu ⋅ force.F) * dΓ
                    end
                end
            end
        end
    end

    return f
end

function stress(xe_ue::Vector{T}, p::Real, el::Integer, model::FEModel{dim}) where {T<:Real,dim}
    σe = zero(Tensor{2,dim})
    ρe = xe_ue[1]
    θe = xe_ue[2]
    ue = xe_ue[3:end]

    cellvalues = model.cellvalues
    for q_point in 1:getnquadpoints(cellvalues)
        dΩ = getdetJdV(cellvalues, q_point)
        ϵ = function_symmetric_gradient(cellvalues, q_point, ue)
        σe += ρe^p * rotate(model.mat.C, θe) ⊡ ϵ * dΩ
    end
    σe = rotate(σe, -θe) # rotate to material coordinates
    σe /= model.elemvol[el]

    return reshape(reinterpret(T, σe), (dim, dim))
end

function hashin_matrix(σ, mat::Material)
    sl, st, slt = tovoigt(σ)
    lim = mat.stress_limits
    CYc = (lim.Yc / 2lim.Stn)^2 - 1

    if st > 0
        modem = "MT"
        IFm = (st / lim.Yt)^2 + (slt / lim.Sln)^2
        FSm = 1 / sqrt(IFm)
    else
        modem = "MC"

        # IF = (st/(2*Stn))^2 + Cyc*(st/Yc) + (slt/Sln)^2
        # F^2 * ((st/(2*Stn))^2 + (slt/Sln)^2) + F * CYc*(st/Yc) - 1 = 0
        a = (st / 2lim.Stn)^2 + (slt / lim.Sln)^2
        b = -CYc * st / lim.Yc

        IFm = a + b
        FSm = (-b + sqrt(b^2 + 4a)) / 2a
    end

    return IFm
end
function hashin_fibre(σ, mat::Material)
    sl, st, slt = tovoigt(σ)
    lim = mat.stress_limits

    if sl > 0
        modef = "FT"
        IFf = (sl / lim.Xt)^2 + (slt / lim.Sln)^2
        FSf = 1 / sqrt(IFf)
    else
        modef = "FC"
        IFf = -sl / lim.Xc
        FSf = 1 / IFf
    end

    return IFf
end

function get_centers(model::FEModel{dim}) where {dim}
    centers = zeros(getncells(model.grid), dim)
    for cell in CellIterator(model.dh)
        id = cellid(cell)
        for node in getcoordinates(model.grid, id)
            centers[id, :] .+= node
        end
        centers[id, :] ./= nnodes_per_cell(model.grid, id)
    end
    return centers
end

function convolution_filter(rmin::Real, model::FEModel)
    rmin == 0 && return I

    centers = get_centers(model)

    # convolution weights
    iH, jH = [], []
    sH = Vector{Float64}()
    for celli in CellIterator(model.dh)
        for cellj in CellIterator(model.dh)
            dist = sqrt(sum((centers[cellid(cellj), :] - centers[cellid(celli), :]) .^ 2))
            dist > rmin && continue
            append!(iH, cellid(celli))
            append!(jH, cellid(cellj))
            append!(sH, rmin - dist)
        end
    end

    # build sparse matrix
    H = sparse(iH, jH, sH)
    H ./= sum(H, dims=2)
    dropzeros!(H)

    return H
end

end # module
