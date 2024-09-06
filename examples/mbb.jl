using Ferrite, FerriteGmsh
using FerriteViz, GLMakie
using TopOpt

# mesh
grid = togrid("models/mbb.msh")
ip = Lagrange{2,RefCube,1}()
qr = QuadratureRule{2,RefCube}(2)

# loads
addnodeset!(grid, "force", x -> x[1] ≈ 0.0 && x[2] ≈ 20.0)
forces = [
    NodalLoad("force", (0.0, -1.0)),
]

# constraints
addfaceset!(grid, "symmetry", x -> x[1] ≈ 0.0)
addnodeset!(grid, "support", x -> x[1] ≈ 60.0 && x[2] ≈ 0.0)
constraints = [
    Dirichlet(:u, getfaceset(grid, "symmetry"), (x, t) -> 0, [1]),
    Dirichlet(:u, getnodeset(grid, "support"), (x, t) -> 0, [2])
]

# material
mat = Isotropic2D(E=1.0, ν=0.3)

# FE model
model = FEModel(mat, grid, ip, qr, constraints, forces)
opts = OptimOpts(100, 0.5, 3, 1.5, 0.2)
@time c, x, u, σ = topopt(model, opts)

println(c)

plotter = MakiePlotter(model.dh, u)

f = Figure()
ax = Axis(f[1, 1], aspect=DataAspect())
hidedecorations!(ax, ticks=false, ticklabels=false)
hidespines!(ax)
cellplot!(plotter, x[:ρ], colormap=:binary)
display(f)
