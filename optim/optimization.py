import time
import os, shutil, glob
from pathlib import Path
import numpy as np
import json, jsonpickle

from ansys.mapdl import core as pymapdl

from .filters import DensityFilter, OrientationFilter

class TopOpt():
    def __init__(self, *, inputfile, res_dir, topology=True, orientation='hashin', jobname='file', echo=True, **kwargs):
        '''
        orientation: None, 'principal', 'hashin'
        '''

        if not topology and orientation is None:
            raise ValueError('topology is False and orientation is None. No optmization to be performed')

        self.ansys_timeout = 180
        
        self.inputfile   = Path(inputfile)
        self.topology    = topology
        self.orientation = orientation
        self.jobname     = jobname
        self.echo        = echo

        self.res_dir     = Path(res_dir)
        self.res_dir.mkdir(parents=True, exist_ok=True)
        
        self.num_elem, self.num_node, self.centers, self.elemvol, self.elmnodes, self.node_coord = [
            kwargs.pop(x,None) for x in ['num_elem','num_node', 'centers', 'elemvol', 'elmnodes', 'node_coord']
        ]
        if self.num_elem is None:
            self.num_elem, self.num_node, self.centers, self.elemvol, self.elmnodes, self.node_coord = self.__get_mesh_data()
        
        # initial setup with default parameters
        self.set_material()
        self.set_optim_options()
        self.set_filters()
        self.set_initial_conditions()
        
    def set_material(self, *, Ex=1, Ey=1, nuxy=0.3, Gxy=1/(2*(1+0.3)), **kwargs):
        self.Ex   = Ex
        self.Ey   = Ey
        self.nuxy = nuxy
        self.Gxy  = Gxy

        self.Xt  = kwargs.pop('Xt', None)
        self.Xc  = kwargs.pop('Xc', None)
        self.Yt  = kwargs.pop('Yt', None)
        self.Yc  = kwargs.pop('Yc', None)
        self.S12 = kwargs.pop('S12', None)
        self.S23 = kwargs.pop('S23', None)

        self.CYc = (self.Yc/(2*self.S23))**2 - 1 if self.Yc is not None else None

        self.fail_eval = self.Xt is not None
        
    def set_filters(self, *, r_rho=0, r_theta=0):
        self.r_rho          = r_rho
        self.density_filter = DensityFilter(r_rho, self.centers)
        
        self.r_theta            = r_theta
        self.orientation_filter = OrientationFilter(r_theta, self.centers)
        
    def set_initial_conditions(self, angle_type='random', **kwargs):
        self.iter = 0
        self.rho  = self.volfrac * np.ones(self.num_elem)

        if not self.topology:
            self.rho = np.ones(self.num_elem)

        if self.orientation is None:
            self.theta = np.zeros(self.num_elem)

        # initial angles are given
        elif angle_type == 'fix':
            theta0 = kwargs.pop('theta0')
            self.theta = np.deg2rad(theta0) * np.ones(self.num_elem)
        
        # random numbers between -pi/2 and pi/2
        elif angle_type == 'random':
            self.theta = np.random.default_rng().uniform(-np.pi/2, np.pi/2, self.num_elem)

        else:
            raise ValueError('Unsupported value for angle_type')
            
    def set_optim_options(self, *, volfrac=0.3, max_iter=200, tol=0):
        self.rho_min = 1e-3
        self.move    = 0.2
        self.volfrac = volfrac

        self.max_iter = max_iter
        self.tol      = tol
        self.penal    = 3

        # Elemental history
        self.rho_hist       = []
        self.theta_hist     = []
        self.sx_hist        = []
        self.sy_hist        = []
        self.sxy_hist       = []
        self.sl_hist        = []
        self.st_hist        = []
        self.slt_hist       = []
        self.s1_hist        = []
        self.s2_hist        = []
        self.fail_MT_hist   = []
        self.fail_MC_hist   = []
        self.fail_FT_hist   = []
        self.fail_FC_hist   = []

        # Scalar history
        self.comp_hist      = []
        self.ang_p_hist     = []
        self.fail_load_hist = []
        self.fail_mode_hist = []
        
        self.time        = 0
        self.fea_time    = 0
        self.update_time = 0
        
    def run(self):
        t0 = time.time()

        if self.topology: self.update = 'density'
        else:             self.update = 'orientation'

        mapdl = pymapdl.launch_mapdl(jobname=self.jobname, run_location=self.res_dir.absolute(), override=True, start_timeout=self.ansys_timeout)
        mapdl.resume(fname=self.inputfile.absolute())

        for _ in range(self.max_iter):
            if self.echo: print('Iteration {:3d}... '.format(self.iter), end=' ')

            # Handle instability on MAPDL connection: if it is broken, creates a new instance and retries
            for attempt in range(5):
                try:
                    self.__fea(mapdl)
                except pymapdl.errors.MapdlExitedError:
                    mapdl = pymapdl.launch_mapdl(jobname=self.jobname, run_location=self.res_dir.absolute(), override=True, start_timeout=self.ansys_timeout)
                    mapdl.resume(fname=self.inputfile.absolute())
                else:
                    break
            else:
                raise ConnectionError('Too many failed attempts to reconnect to MAPDL')
            
            # stopping criterion
            if self.iter >= 1 and np.abs(self.comp_hist[-1]-self.comp_hist[-2])/self.comp_hist[-2] < self.tol:
                if self.orientation != 'hashin': break
                if np.abs(self.fail_load_hist[-1]-self.fail_load_hist[-2])/self.fail_load_hist[-2] < self.tol: break

            if self.update == 'density':    
                self.rho = self.__compliance_oc()
            elif self.update == 'orientation':
                self.theta = self.__orientation_oc()

            self.iter += 1
        else:        
            # Evaluate result from last iteration
            if self.echo: print('Iteration {:3d}... '.format(self.iter), end=' ')
            self.__fea(mapdl)

        mapdl.save()
        mapdl.exit()
        self.__clear_files()
        
        self.time += time.time() - t0
        if self.echo: print('\nTotal elapsed time     {:7.2f}s'.format(self.time))
            
    def print_timing(self):
        print('Total elapsed time     {:7.2f}s'.format(self.time))
        print('FEA time               {:7.2f}s'.format(self.fea_time))
        print('Derivation time        {:7.2f}s'.format(self.deriv_time))
        print('Variable updating time {:7.2f}s'.format(self.mma.update_time))

    def save(self, filename=None):
        if filename is None: filename = self.res_dir / 'topopt.json'
        
        json_str = json.dumps(jsonpickle.encode(self.__dict__))
        with open(filename, 'w') as f:
            f.write(json_str)
    
    @staticmethod
    def load(filename):
        with open(filename, 'r') as f:
            data = json.load(f)
        
        dict_rebuilt = jsonpickle.decode(data)
        solver_rebuilt = TopOpt(**dict_rebuilt)
        solver_rebuilt.__dict__ = dict_rebuilt

        return solver_rebuilt

# -------------------------------------------- Internal functions -------------------------------------------- 
    def __get_mesh_data(self):
        mapdl = pymapdl.launch_mapdl(jobname=self.jobname, run_location=self.res_dir.absolute(), override=True, start_timeout=self.ansys_timeout)
        mapdl.resume(fname=self.inputfile.absolute())
        
        num_elem   = mapdl.mesh.n_elem
        num_node   = mapdl.mesh.n_node
        node_coord = mapdl.mesh.nodes
        elemvol    = mapdl.get_array(entity='elem', item1='geom')
        elmnodes   = np.array(mapdl.mesh.elem)[:,10:] - 1
        
        centers    = np.vstack((
                        mapdl.get_array('elem', item1='cent', it1num='x'),
                        mapdl.get_array('elem', item1='cent', it1num='y'),
                        mapdl.get_array('elem', item1='cent', it1num='z'),
                    )).T

        mapdl.exit()
        self.__clear_files()
        
        return num_elem, num_node, centers, elemvol, elmnodes, node_coord

    def __fea(self, mapdl):
        rho   = self.rho
        theta = self.theta
        
        # Generate 1000 discrete materials properties
        NUM_MAT = 1000
        rho_disc = np.linspace(0.001, 1, NUM_MAT)
        Ex   = rho_disc**self.penal * self.Ex
        Ey   = rho_disc**self.penal * self.Ey
        nuxy = self.nuxy * np.ones(NUM_MAT)
        Gxy  = rho_disc**self.penal * self.Gxy

        mapdl.prep7()
        with mapdl.non_interactive:
            for i in range(NUM_MAT):
                mapdl.mp('ex',i+1,Ex[i])
                mapdl.mp('ey',i+1,Ey[i])
                mapdl.mp('prxy',i+1,nuxy[i])
                mapdl.mp('gxy',i+1,Gxy[i])

            for i in range(self.num_elem):
                mapdl.emodif(i+1,'mat',int(NUM_MAT*rho[i]))
                mapdl.clocal(i+100,0,*self.centers[i,:],np.rad2deg(theta[i]),0,0)
                mapdl.emodif(i+1,'esys',i+100)
                mapdl.csys(0)
        
        mapdl.slashsolu()
        mapdl.solve()

        mapdl.post1()
        mapdl.set('last')

        mapdl.etable('energy','sene')

        mapdl.rsys(0)
        mapdl.etable('sx','s','x')
        mapdl.etable('sy','s','y')
        mapdl.etable('sxy','s','xy')

        mapdl.rsys('solu')
        mapdl.etable('sl','s','x')
        mapdl.etable('st','s','y')
        mapdl.etable('slt','s','xy')

        self.energy = mapdl.get_array(entity='elem', item1='etable', it1num='energy')
        self.c = 2*self.energy.sum()

        sl  = mapdl.get_array(entity='elem', item1='etable', it1num='sl')
        st  = mapdl.get_array(entity='elem', item1='etable', it1num='st')
        slt = mapdl.get_array(entity='elem', item1='etable', it1num='slt')

        sx  = mapdl.get_array(entity='elem', item1='etable', it1num='sx')
        sy  = mapdl.get_array(entity='elem', item1='etable', it1num='sy')
        sxy = mapdl.get_array(entity='elem', item1='etable', it1num='sxy')

        s1 = 0.5*(sx+sy) + np.sqrt((0.5*(sx-sy))**2 + sxy**2)
        s2 = 0.5*(sx+sy) - np.sqrt((0.5*(sx-sy))**2 + sxy**2)
        angle_p = 0.5*np.arctan2(2*sxy, sx-sy)

        self.angle_p = np.where(np.abs(s1) > np.abs(s2), angle_p, angle_p + np.pi/2)
        self.s1 = np.where(np.abs(s1) > np.abs(s2), s1, s2)
        self.s2 = np.where(np.abs(s1) < np.abs(s2), s1, s2)
        self.sl, self.st, self.slt = sl, st, slt

        if self.echo: print('c = {:10.4f}'.format(self.c), end=' ')

        self.rho_hist   += [rho]
        self.theta_hist += [theta]
        self.comp_hist  += [self.c]
        self.ang_p_hist += [self.angle_p]

        self.sx_hist  += [sx]
        self.sy_hist  += [sy]
        self.sxy_hist += [sxy]
        self.sl_hist  += [sl]
        self.st_hist  += [st]
        self.slt_hist += [slt]
        self.s1_hist  += [s1]
        self.s2_hist  += [s2]

        # Failure criteria
        if self.fail_eval:
            FS, mode, IFyt, IFyc, IFxt, IFxc = self.__hashin(sl, st, slt)

            self.fail_load_hist += [FS]
            self.fail_mode_hist += [mode]
            self.fail_MT_hist   += [IFyt]
            self.fail_MC_hist   += [IFyc]
            self.fail_FT_hist   += [IFxt]
            self.fail_FC_hist   += [IFxc]

            if self.echo: print(' FS = {:5.4f}'.format(self.fail_load_hist[-1]), end=' ')

        if self.echo: print()

    def __hashin(self, sl, st, slt):
        # matrix tensile
        IFyt = (st/self.Yt)**2 + (slt/self.S12)**2
        FSmt = 1/np.sqrt(IFyt)

        # matrix compressive
        # IFyc = (st/(2*S23))**2 + ((Yc/(2*S23))**2-1)*(st/Yc) + (slt/S12)**2
        # F**2 * ((st/(2*S23))**2 + (slt/S12)**2) + F * ((Yc/(2*S23))**2-1)*(st/Yc) - 1 = 0
        IFyc = (st/(2*self.S23))**2 + ((self.Yc/(2*self.S23))**2-1)*(st/self.Yc) + (slt/self.S12)**2
        a = (st/(2*self.S23))**2 + (slt/self.S12)**2
        b = ((self.Yc/(2*self.S23))**2-1)*(st/self.Yc)
        c = -1
        FSmc = (-b + np.sqrt(b**2-4*a*c))/(2*a)

        # fiber tensile
        IFxt = (sl/self.Xt)**2 + (slt/self.S12)**2
        FSft = 1/np.sqrt(IFxt)

        # fiber compressive
        IFxc = -sl/self.Xt
        FSfc = 1/IFxc

        FSm  = np.where(st > 0, FSmt, FSmc)
        FSf  = np.where(sl > 0, FSft, FSfc)
        FS   = np.where(FSm < FSf, FSm, FSf)
        mode = np.where(FSm < FSf, np.where(st > 0, 'MT', 'MC'), np.where(sl > 0, 'FT', 'FC'))

        return np.min(FS), mode[np.argmin(FS)], IFyt, IFyc, IFxt, IFxc

    def __compliance_oc(self):
        rho = self.rho
        dcdrho = self.__dcdrho()
        l1 = 0
        l2 = 100000
        while l2-l1 > 1e-4:
            lmid = (l2+l1)/2
            rhonew = np.maximum(self.rho_min, np.maximum(rho-self.move, np.minimum(1, np.minimum(rho+self.move,rho*np.sqrt(-dcdrho/lmid)))))
            if np.dot(rhonew,self.elemvol) - self.volfrac*np.sum(self.elemvol) > 0:
                l1 = lmid
            else:
                l2 = lmid

        return rhonew

    def __orientation_oc(self):
        s1, s2 = self.s1, self.s2
        if self.orientation == 'hashin':
            Qyt  = self.S12**2*(s1+s2)/(self.S12**2-self.Yt**2)/(s1-s2) # matrix tensile failure
            Qyc  = self.S12**2*(4*self.CYc*self.S23**2+self.Yc*(s1+s2))/(self.Yc*(self.S12**2-4*self.S23**2)*(s1-s2)) # matrix compressive failure
            Qy   = np.where(self.st > 0, Qyt, Qyc)
            Qy   = np.maximum(-1, np.minimum(Qy, 1))
            beta = 0.5*np.arccos(Qy)
        elif self.orientation == 'principal':
            beta = 0

        theta = self.angle_p + beta
        theta = np.arctan(np.tan(theta)) # wrap from -90 to 90
        theta = self.orientation_filter.filter(self.rho, theta, self.penal)
        return theta

    def __dcdrho(self):
        dcdrho = -self.penal/self.rho * 2*self.energy
        dcdrho = self.density_filter.filter(self.rho, dcdrho)
        return dcdrho

    def __clear_files(self):
        # clear Ansys temporary files
        for filename in list(set(glob.glob(f'{self.res_dir/self.jobname}.*')) - set(glob.glob(f'{self.res_dir/self.jobname}.db')) - set(glob.glob(f'{self.res_dir/self.jobname}.rst'))): os.remove(filename)
        for filename in glob.glob(f"{self.res_dir/'.__tmp__.*'}"): os.remove(filename)
        for filename in glob.glob(f"{self.res_dir/'*_tmp_*'}"): os.remove(filename)