import os, sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from optim import TopOpt, Plotter
from utils import plot_config

plot_config.config()

solver = TopOpt(inputfile='models/plate_hole.db', res_dir='results/principal/', jobname='principal', topology=False, orientation='principal')
solver.set_material(Ex=53.48e3, Ey=17.7e3, nuxy=0.278, Gxy=5.83e3, Xt=1140, Xc=570, Yt=35, Yc=114, S12=72, S23=36.469)
solver.set_initial_conditions(angle_type='fix', theta0=0)
solver.set_optim_options(max_iter=30)

solver.run()
solver.save()

# solver = TopOpt.load('results/principal/topopt.json')
p = Plotter(solver, dpi=150)
p.plot_convergence()
p.plot_orientation()
p.plot_failure()