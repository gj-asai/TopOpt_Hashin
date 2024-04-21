import os, sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from optim import TopOpt, Plotter
from utils import plot_config

plot_config.config()

# top(60,20,0.5,3.0,1.5)
solver = TopOpt(inputfile='models/mbb.db', res_dir='results/simp/', jobname='simp', topology=True, orientation=None)
solver.set_filters(r_rho=1.5)
solver.set_optim_options(volfrac=0.5, max_iter=100, tol=1e-4)

solver.run()
solver.save()

# solver = TopOpt.load('results/simp/topopt.json')
p = Plotter(solver, dpi=150)
p.plot_convergence()
p.plot_density()