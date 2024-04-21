import numpy as np

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import MaxNLocator
from matplotlib.patches import Polygon
from matplotlib.animation import FuncAnimation
from functools import partial

class Plotter():
    def __init__(self, solver, dpi=250):
        self.solver = solver
        self.quiver_scale = 1/(0.7*np.sqrt(self.solver.elemvol))
        self.dpi = dpi
        
    def plot_convergence(self, end_iter=-1, start_iter=0, compliance=True, failure_load=True, save=True, fig=None, ax1=None, ax2=None):
        if not self.solver.fail_eval:
            failure_load = False
        
        if not compliance and not failure_load:
            raise ValueError('No data to plot. Set compliance or failure_load to True')

        if fig is None:
            if compliance and failure_load:
                fig, (ax1,ax2) = plt.subplots(2, 1, dpi=self.dpi)
            elif compliance:
                fig, ax1 = plt.subplots(dpi=self.dpi)
            elif failure_load:
                fig, ax2 = plt.subplots(dpi=self.dpi)
        else:
            if (compliance and ax1 is None) or (failure_load and ax2 is None):
                raise ValueError('Missing axis: ax1 must not be None for a compliance plot and ax2 must not be None for a failure load plot')

        if compliance:
            comp = self.solver.comp_hist[start_iter:end_iter]
            ax1.cla()
            ax1.plot(np.arange(len(comp))+start_iter, comp)
            ax1.set_ylabel('Compliance')
            ax1.set_xlabel('Iteration')
            ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
        
        if failure_load:
            fail = self.solver.fail_load_hist[start_iter:end_iter]
            ax2.cla()
            ax2.plot(np.arange(len(fail))+start_iter, fail)
            ax2.set_ylabel('Failure load')
            ax2.set_xlabel('Iteration')
            ax2.xaxis.set_major_locator(MaxNLocator(integer=True))
        
        if save: plt.savefig(self.solver.res_dir / 'convergence.png')

    def animate_convergence(self, compliance=True, failure_load=True):
        if not self.solver.fail_eval:
            failure_load = False
        
        if compliance and failure_load:
            fig, (ax1,ax2) = plt.subplots(2, 1, dpi=self.dpi)
        elif compliance:
            fig, ax1 = plt.subplots(dpi=self.dpi)
            ax2 = None
        elif failure_load:
            fig, ax2 = plt.subplots(dpi=self.dpi)
            ax1 = None

        anim = FuncAnimation(fig, partial(self.plot_convergence, compliance=compliance, failure_load=failure_load, save=False, fig=fig, ax1=ax1, ax2=ax2), frames=len(self.solver.rho_hist))
        anim.save(self.solver.res_dir / 'convergence.gif')
        fig.savefig(self.solver.res_dir / 'convergence.png')

    def plot_orientation(self, iteration=-1, principal=True, save=True, fig=None, ax=None):
        if fig is None: fig, ax = plt.subplots(dpi=self.dpi)
        ax.cla()
        ax.set_aspect('equal')
        plt.xlim(np.min(self.solver.node_coord[:,0]),np.max(self.solver.node_coord[:,0]))
        plt.ylim(np.min(self.solver.node_coord[:,1]),np.max(self.solver.node_coord[:,1]))

        x = self.solver.centers[:,0]
        y = self.solver.centers[:,1]
        rho   = self.solver.rho_hist[iteration]
        theta = self.solver.theta_hist[iteration]

        for i in range(self.solver.num_elem):
            xy = self.solver.node_coord[self.solver.elmnodes[i],:2]
            ax.add_patch(Polygon(xy, facecolor='grey', edgecolor='k', alpha=0.2*rho[i]))
        
        u = np.cos(theta)
        v = np.sin(theta)
        ax.quiver(x, y, u, v, color='black', alpha=rho, pivot='mid',
            headwidth=0, headlength=0, headaxislength=0, angles='xy', scale_units='xy', width=3e-3, scale=self.quiver_scale)

        if principal:
            s1    = self.solver.s1_hist[iteration]
            ang_p = self.solver.ang_p_hist[iteration]
            norm  = np.abs(s1)/np.max(np.abs(s1))
            u = norm * np.cos(ang_p)
            v = norm * np.sin(ang_p)

            positive = np.where(s1 > 0)[0]
            ax.quiver(x[positive], y[positive], -u[positive], -v[positive], color='blue', alpha=rho[positive],
                headwidth=3, headlength=3, headaxislength=3, angles='xy', scale_units='xy', width=1.5e-3, scale=self.quiver_scale[positive]/2)
            ax.quiver(x[positive], y[positive], u[positive], v[positive], color='blue', alpha=rho[positive],
                headwidth=3, headlength=3, headaxislength=3, angles='xy', scale_units='xy', width=1.5e-3, scale=self.quiver_scale[positive]/2)

            negative = np.where(s1 < 0)[0]
            ax.quiver(x[negative], y[negative], -u[negative], -v[negative], color='blue', alpha=rho[negative],
                headwidth=-3, headlength=-3, headaxislength=-3, angles='xy', scale_units='xy', width=1.5e-3, scale=self.quiver_scale[negative]/2)
            ax.quiver(x[negative], y[negative], u[negative], v[negative], color='blue', alpha=rho[negative],
                headwidth=-3, headlength=-3, headaxislength=-3, angles='xy', scale_units='xy', width=1.5e-3, scale=self.quiver_scale[negative]/2)

        if save: plt.savefig(self.solver.res_dir / 'orientation.png')

    def animate_orientation(self, principal=True):
        fig, ax = plt.subplots(dpi=self.dpi)
        anim = FuncAnimation(fig, partial(self.plot_orientation, principal=principal, save=False, fig=fig, ax=ax), frames=len(self.solver.rho_hist))
        anim.save(self.solver.res_dir / 'orientation.gif')
        fig.savefig(self.solver.res_dir / 'orientation.png')

    def plot_density(self, iteration=-1, save=True, fig=None, ax=None):
        if fig is None: fig, ax = plt.subplots(dpi=self.dpi)
        ax.cla()
        ax.set_aspect('equal')
        plt.xlim(np.min(self.solver.node_coord[:,0]),np.max(self.solver.node_coord[:,0]))
        plt.ylim(np.min(self.solver.node_coord[:,1]),np.max(self.solver.node_coord[:,1]))
        
        rho = self.solver.rho_hist[iteration]
        for i in range(self.solver.num_elem):
            xy = self.solver.node_coord[self.solver.elmnodes[i],:2]
            ax.add_patch(Polygon(xy, facecolor='k', edgecolor=None, alpha=rho[i]))

        if save: plt.savefig(self.solver.res_dir / 'density.png')

    def animate_density(self):
        fig, ax = plt.subplots(dpi=self.dpi)
        anim = FuncAnimation(fig, partial(self.plot_density, save=False, fig=fig, ax=ax), frames=len(self.solver.rho_hist))
        anim.save(self.solver.res_dir / 'density.gif')
        fig.savefig(self.solver.res_dir / 'density.png')

    def plot_failure(self, iteration=-1, save=True, fig=None, ax=None):
        if fig is None:
            fig, ax = plt.subplots(2, 2, dpi=self.dpi)
            cb = fig.colorbar(cm.ScalarMappable(cmap=cm.jet), ax=ax.flatten().tolist(), shrink=0.8, location='left')
            cb.ax.set_xlabel('Failure\nIndex')
        ax = ax.flatten()

        titles = ['Matrix Tensile', 'Matrix Compressive', 'Fiber Tensile', 'Fiber Compressive']
        for i in range(4):
            ax[i].cla()
            ax[i].set_title(titles[i])
            ax[i].set_aspect('equal')
            ax[i].set_xlim(np.min(self.solver.node_coord[:,0]),np.max(self.solver.node_coord[:,0]))
            ax[i].set_ylim(np.min(self.solver.node_coord[:,1]),np.max(self.solver.node_coord[:,1]))
        
        rho = self.solver.rho_hist[iteration]
        fail_MT_hist = self.solver.fail_MT_hist[iteration]
        fail_MC_hist = self.solver.fail_MC_hist[iteration]
        fail_FT_hist = self.solver.fail_FT_hist[iteration]
        fail_FC_hist = self.solver.fail_FC_hist[iteration]
        for i in range(self.solver.num_elem):
            xy = self.solver.node_coord[self.solver.elmnodes[i],:2]
            ax[0].add_patch(Polygon(xy, facecolor=cm.jet(fail_MT_hist[i]), edgecolor=None, alpha=rho[i]))
            ax[1].add_patch(Polygon(xy, facecolor=cm.jet(fail_MC_hist[i]), edgecolor=None, alpha=rho[i]))
            ax[2].add_patch(Polygon(xy, facecolor=cm.jet(fail_FT_hist[i]), edgecolor=None, alpha=rho[i]))
            ax[3].add_patch(Polygon(xy, facecolor=cm.jet(fail_FC_hist[i]), edgecolor=None, alpha=rho[i]))

        if save: plt.savefig(self.solver.res_dir / 'failure.png')
        
    def animate_failure(self):
        fig, ax = plt.subplots(2, 2, dpi=self.dpi)
        cb = fig.colorbar(cm.ScalarMappable(cmap=cm.jet), ax=ax.flatten().tolist(), shrink=0.8, location='left')
        cb.ax.set_xlabel('Failure\nIndex')
        anim = FuncAnimation(fig, partial(self.plot_failure, save=False, fig=fig, ax=ax), frames=len(self.solver.rho_hist))
        anim.save(self.solver.res_dir / 'failure.gif')
        fig.savefig(self.solver.res_dir / 'failure.png')