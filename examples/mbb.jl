using Ferrite, FerriteGmsh
using FerriteViz, GLMakie

include("../src/TopOpt.jl")
import .TopOpt

# mesh
grid = togrid("models/mbb.msh")
addfaceset!(grid, "symmetry", x -> x[1] ≈ 0.0) # left edge
addnodeset!(grid, "support", x -> x[1] ≈ 168.0 && x[2] ≈ 0.0) # bottom right corner
addnodeset!(grid, "force", x -> x[1] ≈ 0.0 && x[2] ≈ 80.0) # top left corner

# create and solve FE model
model = TopOpt.FEModel(
    mat=TopOpt.Orthotropic2D(El=138e3, Et=11e3, νlt=0.280, Glt=5.50e3, limits=TopOpt.StressLimits(1500, 900, 27, 200, 80, 42.426)), # cf
    # mat=Orthotropic2D(El=53.48e3, Et=17.7e3, νlt=0.278, Glt=5.83e3, limits=StressLimits(1140, 570, 35, 114, 72, 36.469)), # gf
    grid=grid,
    ip=Lagrange{2,RefCube,1}(), # linear elements
    qr=QuadratureRule{2,RefCube}(2), # 2 point quadrature
    constraints=[
        Dirichlet(:u, getfaceset(grid, "symmetry"), (x, t) -> 0.0, [1]), # block x displacement
        Dirichlet(:u, getnodeset(grid, "support"), (x, t) -> 0.0, [2]), # block y displacement
    ],
    loads=[
        TopOpt.NodalLoad("force", (0.0, -100.0)),
    ],
)

opts = TopOpt.OptimOpts(maxiter=300, volfrac=0.4, rρ=5.0, rθ=5.0, reltol_comp=5e-4, reltol_fail=5e-4)
@time ρ, θ, c_hist, IFm_hist, IFf_hist, IFm, IFf = TopOpt.topopt(model, opts)
IF = max(IFm_hist[end], IFf_hist[end])

## plot convergence
f1 = Figure(resolution=(700, 500))

ax1 = Axis(f1[1, 1][1, 1], ylabel="Compliance (N.mm)",
    xtickalign=1, ytickalign=1, xticksmirrored=true, yticksmirrored=true)
lines!(ax1, 0:length(c_hist)-1, c_hist)

ax2 = Axis(f1[1, 1][2, 1], ylabel="Failure index", xlabel="Iteration",
    xtickalign=1, ytickalign=1, xticksmirrored=true, yticksmirrored=true)
lines!(ax2, 0:length(IFm_hist)-1, IFm_hist, label="Matrix")
lines!(ax2, 0:length(IFf_hist)-1, IFf_hist, label="Fibre")
axislegend(ax2)

wait(display(f1))
# save("mbb_fail_convergence.png", f1)

## plot structure
f2 = Figure(resolution=(900, 450))
plotter = MakiePlotter(model.dh, zeros(getnnodes(model.grid)))

ax = Axis(f2[1, 1], aspect=DataAspect(),
    title="Compliance = $(round(c_hist[end], digits=2)) N.mm, IF = $(round(IF, digits=3))")
hidedecorations!(ax, ticks=false, ticklabels=false)
hidespines!(ax)

centers = TopOpt.get_centers(model)
GLMakie.arrows!(centers[:, 1], centers[:, 2], cos.(θ), sin.(θ),
    arrowsize=0, lengthscale=1.5, align=:center, color=ρ, colormap=:binary)

wait(display(f2))
# save("mbb_fail.png", f2)

## plot failure
f3 = Figure(resolution=(600, 600))

# matrix
axm = Axis(f3[1, 1][1, 1], aspect=DataAspect(), title="Matrix failure index")
hidedecorations!(axm)
hidespines!(axm)

pm = cellplot!(plotter, IFm, colormap=:viridis, colorrange=(0, IF))

# fibre
axf = Axis(f3[1, 1][2, 1], aspect=DataAspect(), title="Fibre failure index")
hidedecorations!(axf)
hidespines!(axf)

pf = cellplot!(plotter, IFf, colormap=:viridis, colorrange=(0, IF))
f3[2, 1] = GLMakie.Colorbar(f3[1, 1][2, 1], limits=(0, IF), colormap=:viridis, vertical=false)

wait(display(f3))
# save("mbb_fail_IF.png", f3)

nothing
