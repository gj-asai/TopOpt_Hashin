module TopOpt

using Ferrite, Tensors
using SparseArrays, LinearAlgebra
using ForwardDiff
using Nonconvex, NonconvexMMA
Nonconvex.@load MMA

export StressLimits, Isotropic2D, Orthotropic2D
export NodalLoad, LinearLoad
export FEModel, get_centers
export OptimOpts, topopt

include("material.jl")
include("femodel.jl")
include("fea.jl")
include("optim_core.jl")

function topopt(model::FEModel, opts::OptimOpts)
    # initial values
    x = zeros(2 * getncells(model.grid))
    x[1:2:end] .= opts.volfrac
    x[2:2:end] .= deg2rad(0)
    # x[2:2:end] .= 2π * rand(Float64, getncells(model.grid)) .- π

    penal = 3
    data = OptimData(3, opts, model)

    x .= min_compliance(x, data)
    x .= min_failure(x, data)

    ρ = x[1:2:end]
    θ = x[2:2:end]
    θ = data.Hθ * θ # smooth out last design

    return ρ, θ, data.c_hist, data.IFm_hist, data.IFf_hist, data.data.IFm, data.data.IFf
end

end # module
