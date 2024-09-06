using Ferrite, FerriteGmsh
using FerriteViz, GLMakie
using TopOpt

# mesh
grid = togrid("models/cantilever_plate.msh")
ip = Lagrange{2,RefTetrahedron,1}()
qr = QuadratureRule{2,RefTetrahedron}(1)

# loads
addfaceset!(grid, "force", x -> x[1] ≈ 70.0)
forces = [
    # LinearLoad("force", (0.0, 15.3)),
    LinearLoad("force", (0.0, 1.0)),
]

# constraints
addnodeset!(grid, "fix", x -> x[1] ≈ 0.0 && x[2] ≈ 0.0)
addfaceset!(grid, "support", x -> x[1] ≈ 0.0)
constraints = [
    Dirichlet(:u, getfaceset(grid, "support"), (x, t) -> 0, [1]),
    Dirichlet(:u, getnodeset(grid, "fix"), (x, t) -> 0, [2]),
]

# material
limits_cf = StressLimits(1500, 900, 27, 200, 80, 42.426)
limits_gf = StressLimits(1140, 570, 35, 114, 72, 36.469)
limits_cf3d = StressLimits(493.9, 323.9, 13.5, 20.25, 35, 8.482)
cf = Orthotropic2D(El=138e3, Et=11e3, νlt=0.280, Glt=5.50e3, limits=limits_cf)
gf = Orthotropic2D(El=53.48e3, Et=17.7e3, νlt=0.278, Glt=5.83e3, limits=limits_gf)
cf3d = Orthotropic2D(El=50e3, Et=2.322e3, νlt=0.333, Glt=0.624e3, limits=limits_cf3d)

# FE model
model = FEModel(cf, grid, ip, qr, constraints, forces)
opts = OptimOpts(100, 1.0, 3, 1.5, 0.2)
@time c_hist, x, u, σ, FS_hist, mode_hist, IFm, IFf = topopt(model, opts)

display(c_hist[end] / 0.51)
display(FS_hist[end] * 0.51)
println(mode_hist[end])

f = Figure()
ax1 = Axis(f[1, 1])
lines!(ax1, 1:length(c_hist), c_hist ./ 0.51)
ax2 = Axis(f[2, 1])
lines!(ax2, 1:length(FS_hist), FS_hist .* 0.51)
display(f)

# plotter = MakiePlotter(model.dh, u)

# f = Figure()
# ax = Axis(f[1, 1], aspect=DataAspect())
# hidedecorations!(ax, ticks=false, ticklabels=false)
# hidespines!(ax)
# p = cellplot!(plotter, IF, colormap=:viridis)
# FerriteViz.wireframe!(grid, plotnodes=false)
# f[1, 2] = GLMakie.Colorbar(f[1, 1], p)
# # display(f)

# f = Figure()
# ax = Axis(f[1, 1], aspect=DataAspect())
# hidedecorations!(ax, ticks=false, ticklabels=false)
# hidespines!(ax)
# FerriteViz.wireframe!(grid, plotnodes=false)

# centers = get_centers(model)
# GLMakie.arrows!(centers[:, 1], centers[:, 2], cos.(x[:θ]), sin.(x[:θ]), arrowsize=0, lengthscale=4.0, align=:center)
# display(f)
