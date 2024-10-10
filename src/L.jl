using Ferrite, FerriteGmsh

using FerriteViz, GLMakie

include("./TopOpt.jl")
import .TopOpt

# mesh
grid = togrid("models/L.msh")
addfaceset!(grid, "fix", x -> x[2] ≈ 200.0) # top edge
addnodeset!(grid, "force", x -> x[1] ≈ 200.0 && x[2] ≈ 50.0) # bottom right corner
# addfaceset!(grid, "force", x -> x[1] ≈ 200.0) # right edge

# create and solve FE model
model = TopOpt.FEModel(
    mat=TopOpt.Orthotropic2D(El=138e3, Et=11e3, νlt=0.280, Glt=5.50e3, limits=TopOpt.StressLimits(1500, 900, 27, 200, 80, 42.426)), # cf
    # mat=Orthotropic2D(El=53.48e3, Et=17.7e3, νlt=0.278, Glt=5.83e3, limits=StressLimits(1140, 570, 35, 114, 72, 36.469)), # gf
    grid=grid,
    ip=Lagrange{2,RefCube,1}(), # linear elements
    qr=QuadratureRule{2,RefCube}(2), # 2 point quadrature
    constraints=[
        Dirichlet(:u, getfaceset(grid, "fix"), x -> 0 * x),
    ],
    loads=[
        # TopOpt.LinearLoad("force", (0.0, -1.0)),
        TopOpt.NodalLoad("force", (0.0, -1.0)),
    ],
)

opts = TopOpt.OptimOpts(maxiter=150, volfrac=0.4, rρ=10.0, rθ=10.0, reltol=1e-4, objective=TopOpt.Compliance())
@time ρ, θ, c_hist, IF_hist, IFm, IFf = TopOpt.topopt(model, opts)

# plot convergence
f1 = Figure()

ax1 = Axis(f1[1, 1], ylabel="Compliance (N.mm)",
    xtickalign=1, ytickalign=1, xticksmirrored=true, yticksmirrored=true)
lines!(ax1, 0:length(c_hist)-1, c_hist)

ax2 = Axis(f1[2, 1], ylabel="Failure index", xlabel="Iteration",
    xtickalign=1, ytickalign=1, xticksmirrored=true, yticksmirrored=true)
lines!(ax2, 0:length(IF_hist)-1, IF_hist)

wait(display(f1))

# plot failure
f2 = Figure()
plotter = MakiePlotter(model.dh, zeros(getnnodes(model.grid)))

# structure
ax = Axis(f2[1, 1], aspect=DataAspect())
hidedecorations!(ax, ticks=false, ticklabels=false)
hidespines!(ax)

centers = TopOpt.get_centers(model)
GLMakie.arrows!(centers[:, 1], centers[:, 2], cos.(θ), sin.(θ),
    arrowsize=0, lengthscale=3.0, align=:center, color=ρ, colormap=:binary)

# matrix
axm = Axis(f2[2, 1][1, 1][1, 1], aspect=DataAspect(), title="Matrix failure index")
hidedecorations!(axm)
hidespines!(axm)

pm = cellplot!(plotter, IFm, colormap=:viridis)
f2[2, 1][1, 1][2, 1] = GLMakie.Colorbar(f2[2, 1][1, 1][1, 1], pm, vertical=false)

# fibre
axf = Axis(f2[2, 1][1, 2][1, 1], aspect=DataAspect(), title="Fibre failure index")
hidedecorations!(axf)
hidespines!(axf)

pf = cellplot!(plotter, IFf, colormap=:viridis)
f2[2, 1][1, 2][2, 1] = GLMakie.Colorbar(f2[2, 1][1, 2][1, 1], pf, vertical=false)

wait(display(f2))

nothing
