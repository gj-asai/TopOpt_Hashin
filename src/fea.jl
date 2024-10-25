mutable struct FEData
    K::AbstractArray
    f::AbstractVector
    ∂Ke∂x::AbstractVector{AbstractArray}
    Keinv::AbstractVector{AbstractArray}

    u::AbstractVector
    c::Real
    IFm::AbstractVector
    IFf::AbstractVector

    dcdx::AbstractVector
    dIFmdx::AbstractVector
    dIFfdx::AbstractVector
end
function FEData(model::FEModel{dim}) where {dim}
    f = global_force(model)

    # preallocate stiffness and sensitivities
    n_basefuncs = getnbasefunctions(model.cellvalues)
    K = create_sparsity_pattern(model.dh)
    ∂Ke∂x = fill(zeros(n_basefuncs, n_basefuncs), 2 * getncells(model.grid))
    Keinv = fill(zeros(n_basefuncs, n_basefuncs), getncells(model.grid))

    # preallocate solution
    u = zeros(dim * getnnodes(model.grid))
    c = 0

    # preallocate failure indices
    IFm = zeros(getncells(model.grid))
    IFf = zeros(getncells(model.grid))

    # preallocate sensitivities
    dcdx = zeros(2 * getncells(model.grid))
    dIFmdx = zeros(2 * getncells(model.grid))
    dIFfdx = zeros(2 * getncells(model.grid))

    return FEData(K, f, ∂Ke∂x, Keinv, u, c, IFm, IFf, dcdx, dIFmdx, dIFfdx)
end

function fea(x::AbstractVector, penal::Real, model::FEModel{dim}, data::FEData) where {dim}
    K = data.K
    f = data.f
    ∂Ke∂x = data.∂Ke∂x
    Keinv = data.Keinv
    u = data.u
    IFm = data.IFm
    IFf = data.IFf
    dcdx = data.dcdx
    dIFmdx = data.dIFmdx
    dIFfdx = data.dIFfdx

    # assemble linear system
    global_stiffness!(K, ∂Ke∂x, x, penal, model)
    apply!(K, f, model.ch)

    # get stiffness matrix slices corresponding to each element and invert
    for cell in CellIterator(model.dh)
        Keinv[cellid(cell)] = K[celldofs(cell), celldofs(cell)] |> collect |> inv
    end

    # solve linear system
    u .= K \ f

    # compliance and sensitivity
    data.c = u' * K * u
    for cell in CellIterator(model.dh)
        ue = u[celldofs(cell)]
        dcdx[2*cellid(cell)-1] = -ue' * ∂Ke∂x[2*cellid(cell)-1] * ue
        dcdx[2*cellid(cell)] = -ue' * ∂Ke∂x[2*cellid(cell)] * ue
    end

    ρ = x[1:2:end]
    θ = x[2:2:end]

    # failure indices and sensitivities
    for cell in CellIterator(model.dh)
        e = cellid(cell)
        ue = u[celldofs(cell)]
        xe_ue = [ρ[e], θ[e], ue...]

        # compute stresses and sensitivities
        # using ForwardDiff's AD because Tensors's AD does not support using Vec's of arbitrary dimension
        # requires stress function to work with Arrays instead of Tensors, need to convert back the results
        σ = stress(xe_ue, penal, e, model) |> Tensor{2,dim}
        jac = ForwardDiff.jacobian(xe_ue -> stress(xe_ue, penal, e, model), xe_ue)
        ∂σ∂ρ = reshape(jac[:, 1], (dim, dim)) |> Tensor{2,dim}
        ∂σ∂θ = reshape(jac[:, 2], (dim, dim)) |> Tensor{2,dim}
        ∂σ∂u = [reshape(d, (dim, dim)) for d in eachcol(jac[:, 3:end])] .|> Tensor{2,dim}

        # derivatives of u from the equilibrium Ku = f
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

# combines xe and ue in the same vector to get all derivatives via AD in the same run
# ForwardDiff's AD requires output to be an Array - Tensor's docs says it is ok to just convert the output
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
        IFm = (st / lim.Yt)^2 + (slt / lim.Sln)^2
        FSm = 1 / sqrt(IFm)
    else
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
        IFf = (sl / lim.Xt)^2 + (slt / lim.Sln)^2
        FSf = 1 / sqrt(IFf)
    else
        IFf = -sl / lim.Xc
        FSf = 1 / IFf
    end

    return IFf
end
