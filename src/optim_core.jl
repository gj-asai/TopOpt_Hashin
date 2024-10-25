@kwdef struct OptimOpts
    maxiter::Integer
    volfrac::Real
    rρ::Real
    rθ::Real
    reltol_comp::Real
    reltol_fail::Real
end

mutable struct OptimData
    Hρ::AbstractArray
    Hθ::AbstractArray
    penal::Real
    model::FEModel
    data::FEData
    opts::OptimOpts
    c_hist::AbstractVector
    IFm_hist::AbstractVector
    IFf_hist::AbstractVector
end
function OptimData(penal::Real, opts::OptimOpts, model::FEModel)
    Hρ = convolution_filter(opts.rρ, model)
    Hθ = convolution_filter(opts.rθ, model)
    data = FEData(model)

    c_hist = Float64[]
    IFm_hist = Float64[]
    IFf_hist = Float64[]

    return OptimData(Hρ, Hθ, penal, model, data, opts, c_hist, IFm_hist, IFf_hist)
end

compose(ρ, θ) = zip(ρ, θ) |> Iterators.flatten |> collect

function compliance(x::AbstractVector, optim_data::OptimData)
    x[2:2:end] .= optim_data.Hθ * x[2:2:end]
    fea(x, optim_data.penal, optim_data.model, optim_data.data)
    return optim_data.data.c
end
function dcompliance(x::AbstractVector, optim_data::OptimData)
    data = optim_data.data
    data.dcdx[1:2:end] .= optim_data.Hρ * (x[1:2:end] .* data.dcdx[1:2:end]) ./ x[1:2:end]
    return data.dcdx
end

function failure(x::AbstractVector, optim_data::OptimData)
    data = optim_data.data
    x[2:2:end] .= optim_data.Hθ * x[2:2:end]
    fea(x, optim_data.penal, optim_data.model, data)
    return sum([data.IFm; data.IFf] .^ 2)^(1 / 2)
end
function dfailure(x::AbstractVector, optim_data::OptimData)
    data = optim_data.data
    data.dIFmdx[1:2:end] .= optim_data.Hρ * (x[1:2:end] .* data.dIFmdx[1:2:end]) ./ x[1:2:end]
    data.dIFfdx[1:2:end] .= optim_data.Hρ * (x[1:2:end] .* data.dIFfdx[1:2:end]) ./ x[1:2:end]
    return sum([data.IFm; data.IFf] .^ 2) .^ (-1 / 2) .* (compose(data.IFm, data.IFm) .* data.dIFmdx + compose(data.IFf, data.IFf) .* data.dIFfdx)
end

function volume(x::AbstractVector, optim_data::OptimData)
    model = optim_data.model
    opts = optim_data.opts
    return x[1:2:end] ⋅ model.elemvol / sum(model.elemvol) - opts.volfrac
end
function dvolume(x::AbstractVector, optim_data::OptimData)
    model = optim_data.model
    ∂g∂x = zero(x)
    ∂g∂x[1:2:end] .= model.elemvol / sum(model.elemvol)
    return ∂g∂x
end

function min_compliance(x0::AbstractVector, optim_data::OptimData)
    data = optim_data.data
    model = optim_data.model
    opts = optim_data.opts

    optim = Model(CustomGradFunction(x -> compliance(x, optim_data), x -> dcompliance(x, optim_data)))
    for i = 1:getncells(model.grid)
        addvar!(optim, 1e-3, 1) # ρ
        addvar!(optim, -π, π) # θ
    end
    add_ineq_constraint!(optim, CustomGradFunction(x -> volume(x, optim_data), x -> dvolume(x, optim_data)))

    function post(solution; update=false)
        push!(optim_data.c_hist, data.c)
        push!(optim_data.IFm_hist, max(data.IFm...))
        push!(optim_data.IFf_hist, max(data.IFf...))
    end

    r = optimize(
        optim,
        MMA(),
        x0,
        options=MMAOptions(
            maxiter=opts.maxiter,
            convcriteria=GenericCriteria(),
            tol=Tolerance(x=0.0, fabs=0.0, frel=opts.reltol_comp),
        ),
        callback=post,
    )
    return r.minimizer
end

function min_failure(x, optim_data)
    data = optim_data.data
    model = optim_data.model
    opts = optim_data.opts

    ρ = x[1:2:end]

    optim = Model(CustomGradFunction(θ -> failure(compose(ρ, θ), optim_data), θ -> dfailure(compose(ρ, θ), optim_data)[2:2:end]))

    for i = 1:getncells(model.grid)
        addvar!(optim, -π, π) # θ
    end

    function post(solution; update=false)
        push!(optim_data.c_hist, data.c)
        push!(optim_data.IFm_hist, max(data.IFm...))
        push!(optim_data.IFf_hist, max(data.IFf...))
    end

    r = optimize(
        optim,
        MMA(),
        x[2:2:end],
        options=MMAOptions(
            maxiter=opts.maxiter,
            convcriteria=GenericCriteria(),
            tol=Tolerance(x=0.0, fabs=0.0, frel=opts.reltol_fail),
        ),
        callback=post,
    )

    return compose(ρ, r.minimizer)
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
