module TopOpt

using Ferrite
using Ferrite: nnodes_per_cell
using Tensors, SparseArrays, LinearAlgebra
using ForwardDiff

export StressLimits, Isotropic2D, Orthotropic2D
export NodalLoad, LinearLoad
export FEModel, OptimOpts
export fea, topopt
export get_centers

export global_stiffness, global_force

struct StressLimits
    Xt
    Xc
    Yt
    Yc
    Sln
    Stn
end
StressLimits() = StressLimits(0, 0, 0, 0, 0, 0)

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

function FEModel(mat::Material{dim}, grid::Grid, ip::Interpolation, qr::QuadratureRule{dim,shape}, constraints::Vector{Dirichlet}, loads::Vector{<:Load}) where {dim,shape}
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

struct OptimOpts
    maxiter::Integer
    volfrac::Real
    p::Real
    rmin::Real
    move::Real
end

function topopt(model::FEModel, opts::OptimOpts)
    x = Dict(
        :ρ => fill(opts.volfrac, getncells(model.grid)),
        :θ => zeros(getncells(model.grid))
    )

    H = convolution_filter(opts.rmin, model)

    # history
    c_hist = []
    FS_hist = []
    mode_hist = []

    # final result only
    u, σ, IFm, IFf = 0, 0, 0, 0
    for _ in 1:opts.maxiter
        (c, ∂c∂x), (u, σ), (FS, mode, IFm, IFf) = fea(x, opts.p, model)
        push!(c_hist, c)
        push!(FS_hist, FS)
        push!(mode_hist, mode)

        # apply filter and update design variables
        ∂c∂x[:ρ] .= H * (x[:ρ] .* ∂c∂x[:ρ]) ./ x[:ρ]
        update_x!(x, ∂c∂x, σ, model, opts)
    end
    return c_hist, x, u, σ, FS_hist, mode_hist, IFm, IFf
end

function fea(x::Dict{Symbol,<:Vector{<:Real}}, p::Real, model::FEModel{dim}) where {dim}
    # assemble linear system
    K, ∂Ke∂x = global_stiffness(x, p, model)
    f = global_force(model) # TODO not necessary to rebuild in every iteration
    apply!(K, f, model.ch)

    # solve
    u = K \ f

    # compliance
    c = u' * K * u

    # sensitivities
    ∂c∂x = Dict(var => zero(xi) for (var, xi) in x)
    for cell in CellIterator(model.dh)
        ue = u[celldofs(cell)]
        for var in keys(∂c∂x)
            ∂c∂x[var][cellid(cell)] = -ue' * ∂Ke∂x[var][cellid(cell)] * ue
        end
    end

    # stress - averaged inside element
    σ = zeros(Tensor{2,dim}, getncells(model.grid))
    for cell in CellIterator(model.dh)
        reinit!(model.cellvalues, cell)
        e = cellid(cell)
        for q_point in 1:getnquadpoints(model.cellvalues)
            dΩ = getdetJdV(model.cellvalues, q_point)
            ϵ = function_symmetric_gradient(model.cellvalues, q_point, u[celldofs(cell)])
            σ[e] += x[:ρ][e]^p * rotate(model.mat.C, x[:θ][e]) ⊡ ϵ * dΩ
        end
        σ[e] /= model.elemvol[e]
    end

    # failure
    FS = zeros(getncells(model.grid))
    mode = ["" for _ in 1:getncells(model.grid)]
    IFm = zeros(getncells(model.grid))
    IFf = zeros(getncells(model.grid))
    for cell in CellIterator(model.dh)
        e = cellid(cell)
        sl, st, slt = tovoigt(rotate(σ[e], -x[:θ][e]))
        FS[e], mode[e], (IFm[e], IFf[e]) = hashin(σ[e], x[:θ][e], model.mat)
    end

    crit = argmin(FS)
    FS = FS[crit]
    mode = mode[crit]

    return (c, ∂c∂x), (u, σ), (FS, mode, IFm, IFf)
end

function update_x!(x::Dict{Symbol,<:Vector{<:Real}}, ∂c∂x::Dict{Symbol,<:Vector{<:Real}}, σ, model::FEModel, opts::OptimOpts)
    # update density
    ρ = x[:ρ]
    ∂c∂ρ = ∂c∂x[:ρ]
    ρnew = zero(x[:ρ])

    λ1, λ2 = 0, 100000
    while λ2 - λ1 > 1e-4
        λmid = (λ1 + λ2) / 2
        @. ρnew = max(1e-3, max(ρ - opts.move, min(ρ * sqrt(-∂c∂ρ / λmid), min(ρ + opts.move, 1))))
        if ρnew ⋅ model.elemvol - opts.volfrac * sum(model.elemvol) > 0
            λ1 = λmid
        else
            λ2 = λmid
        end
    end

    x[:ρ] .= ρnew

    ## update orientation
    # lim = model.mat.stress_limits
    # CYc = (lim.Yc / 2lim.Stn)^2 - 1
    # for cell in CellIterator(model.dh)
    #     e = cellid(cell)
    #     principal = eigen(σ[e])
    #     s1, s2 = sort(abs.(principal.values), rev=true)
    #     sl, st, slt = tovoigt(rotate(σ[e], -x[:θ][e]))

    #     if st > 0
    #         Qy = lim.Sln^2 * (s1 + s2) / (lim.Sln^2 - lim.Yt^2) / (s1 - s2) # matrix tensile
    #     else
    #         Qy = lim.Sln^2 * (4 * CYc * lim.Stn^2 + lim.Yc * (s1 + s2)) / (lim.Yc * (lim.Sln^2 - 4 * lim.Stn^2) * (s1 - s2)) # matrix compressive
    #     end
    #     Qy = min(-1, max(Qy, 1))

    #     idx_max = argmax(abs.(principal.values))
    #     direction = principal.vectors[:, idx_max]
    #     x[:θ][e] = atan(direction[2] / direction[1]) + 0.5 * acos(Qy)
    # end
end

# function global_stiffness(x::Dict{Symbol,<:Vector{<:Real}}, p::Real, model::FEModel)
function global_stiffness(x, p::Real, model::FEModel)
    # preallocate global stiffness
    K = create_sparsity_pattern(model.dh)
    assembler = start_assemble(K)

    # preallocate elemental stiffness
    n_basefuncs = getnbasefunctions(model.cellvalues)
    Ke = zeros(n_basefuncs, n_basefuncs)

    # Dict{Symbol, Vector{Matrix}}
    # ∂Ke∂x[design variable][element] is a square matrix of order n_basefuncs with the derivatives of Ke
    ∂Ke∂x = Dict(var => fill(zeros(n_basefuncs, n_basefuncs), getncells(model.grid)) for var in keys(x))

    for cell in CellIterator(model.dh)
        reinit!(model.cellvalues, cell)

        e = cellid(cell)
        xe = [var[e] for var in values(x)]

        # compute Ke and derivatives
        jac = ForwardDiff.jacobian((Ke, xe) -> element_stiffness!(Ke, xe, p, model), Ke, xe)

        for (i, var) in enumerate(keys(∂Ke∂x))
            ∂Ke∂x[var][e] = reshape(jac[:, i], (n_basefuncs, n_basefuncs))
        end

        assemble!(assembler, celldofs(cell), Ke)
    end

    return K, ∂Ke∂x
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

function hashin(σ, θ, mat::Material)
    sl, st, slt = tovoigt(rotate(σ, -θ))
    lim = mat.stress_limits
    CYc = (lim.Yc / 2lim.Stn)^2 - 1

    # matrix failure
    if st > 0
        modem = "MT"
        IFm = (st / lim.Yt)^2 + (slt / lim.Sln)^2
        FSm = 1 / sqrt(IFm)
    else
        modem = "MC"

        # IF = (st/(2*Stn))^2 + Cyc*(st/Yc) + (slt/Sln)^2
        # F^2 * ((st/(2*Stn))^2 + (slt/Sln)^2) + F * CYc*(st/Yc) - 1 = 0
        a = (st / 2lim.Stn)^2 + (slt / lim.Sln)^2
        b = CYc * st / lim.Yc
        c = -1

        IFm = a + b
        FSm = (-b + sqrt(b^2 - 4 * a * c)) / 2a
    end

    # fiber failure
    if sl > 0
        modef = "FT"
        IFf = (sl / lim.Xt)^2 + (slt / lim.Sln)^2
        FSf = 1 / sqrt(IFf)
    else
        modef = "FC"
        IFf = -sl / lim.Xc
        FSf = 1 / IFf
    end

    (FS, mode) = FSm < FSf ? (FSm, modem) : (FSf, modef)
    return FS, mode, (IFm, IFf)
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
