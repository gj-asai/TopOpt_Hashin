using Ferrite, FerriteGmsh
using FerriteViz, GLMakie

include("TopOpt.jl")
using .TopOpt

# mesh
grid = togrid("models/mbb.msh")
addfaceset!(grid, "symmetry", x -> x[1] ≈ 0.0) # left edge
addnodeset!(grid, "support", x -> x[1] ≈ 60.0 && x[2] ≈ 0.0) # bottom right corner
addnodeset!(grid, "force", x -> x[1] ≈ 0.0 && x[2] ≈ 20.0) # top left corner

# create and solve FE model
model = FEModel(
    mat=Orthotropic2D(El=138e3, Et=11e3, νlt=0.280, Glt=5.50e3, limits=StressLimits(1500, 900, 27, 200, 80, 42.426)), # cf
    # mat=Orthotropic2D(El=53.48e3, Et=17.7e3, νlt=0.278, Glt=5.83e3, limits=StressLimits(1140, 570, 35, 114, 72, 36.469)), # gf
    grid=grid,
    ip=Lagrange{2,RefCube,1}(), # linear elements
    qr=QuadratureRule{2,RefCube}(2), # 2 point quadrature
    constraints=[
        Dirichlet(:u, getfaceset(grid, "symmetry"), (x, t) -> 0.0, [1]), # block x displacement
        Dirichlet(:u, getnodeset(grid, "support"), (x, t) -> 0.0, [2]), # block y displacement
    ],
    loads=[
        NodalLoad("force", (0.0, -1.0)),
    ],
)
@time ρ, θ, c_hist, FS_hist, mode_hist, IFm, IFf = topopt(model, OptimOpts(200, 0.4, 3, 1.5, 2.0, 3e-5))

# plot convergence
f1 = Figure()

ax1 = Axis(f1[1, 1], ylabel="Compliance (N.mm)",
    xtickalign=1, ytickalign=1, xticksmirrored=true, yticksmirrored=true)
lines!(ax1, 0:length(c_hist)-1, c_hist)

ax2 = Axis(f1[2, 1], ylabel="Failure load (N)", xlabel="Iteration",
    xtickalign=1, ytickalign=1, xticksmirrored=true, yticksmirrored=true)
lines!(ax2, 0:length(FS_hist)-1, FS_hist)

wait(display(f1))

# plot failure
f2 = Figure()
plotter = MakiePlotter(model.dh, zeros(getnnodes(model.grid)))

# structure
ax = Axis(f2[1, 1], aspect=DataAspect())
hidedecorations!(ax, ticks=false, ticklabels=false)
hidespines!(ax)

centers = get_centers(model)
GLMakie.arrows!(centers[:, 1], centers[:, 2], cos.(θ), sin.(θ),
    arrowsize=0, lengthscale=0.9, align=:center, color=ρ, colormap=:binary)

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

display(f2)

nothing
