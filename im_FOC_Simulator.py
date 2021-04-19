# encoding:utf-8
# h.c.
# 2015.6.20
 
from pylab import * # numpy + matplotlib
import copy # copy the old value of im.st
from scipy.integrate import ode
# from scipy.io import loadmat
# import pdb
from pprint import pprint
import h5py 
import os


import random
def prbs():
    while True:
        yield choice([-1., 1.])
## SysId needs
# random.seed(20)
# sgen1 = prbs()
# sgen2 = prbs()

def cmd(t, delay_time, **kwarg):
    '''
    # not ready for negative rotate (the limit is positive)
    # try add a LPF for rpm cmd here
    '''
    ret_rpm, ret_Tl = 0., 0.
    for (instant, rpm) in kwarg['rpm']:
        if t-delay_time < instant:
            break
        ret_rpm = rpm
    for (instant, Tload) in kwarg['Tload']:
        if t-delay_time < instant:
            break
        ret_Tl = Tload
    return ret_rpm, ret_Tl

def hpf(x, x_last, y_last, dt, RC):
    alp = RC / (RC + dt)
    state = y_last - x_last
    # return alp * (y_last + x - x_last)
    return alp * (state + x)

from dopri5 import dopri54
class Simulator(object):
    ''' Package all numerate integrate into one, which is the actually valid practice.
    '''
    def __init__(self):
        self.backend = 'dopri5'
        # self.backend = 'dop853'
        # self.backend = 'vode'

        global initials
        # initials = zeros(7+4+3); initials[-3]=1.1; initials[-2]=8.
        # initials = zeros(7+4+3); initials[-3]=1.; initials[-2]=7.
        initials = zeros(7+4+3); initials[-3]=0.9; initials[-2]=7.
        # initials = zeros(7+4+3); initials[-3]=1.1; initials[-2]=6.

        # initials = \
        # [ -9.58128996e+00,  -2.48649405e+00,   3.36390916e+00,   5.51131720e+00,
        #    7.94680277e-03,  -4.08831145e-02,   0.00000000e+00,  -8.91673011e-01,
        #    5.22575389e-01,  -1.07786943e+00,   4.74254630e-01,   1.,
        #    7.,               0.] # 保证磁链已经达到稳态，转矩电流也已经补偿了初始的负载转矩

        # initials = \
        # [ -9.58128996e+00,  -2.48649405e+00,   3.36390916e+00,   5.51131720e+00,
        #    7.94680277e-03,  -4.08831145e-02,   0.00000000e+00,  -8.91673011e-01,
        #    5.22575389e-01,  -1.07786943e+00,   4.74254630e-01,   1.1,
        #    8.,               0.] # 保证磁链已经达到稳态，转矩电流也已经补偿了初始的负载转矩

        # self.solver = ode(self.dynamics).set_integrator(self.backend).set_initial_value(initials, 0.)
        # self.solver = ode(self.dynamics).set_integrator(self.backend, atol=1e2, rtol=0, nsteps=1, first_step=1e-4, max_step=1e-4).set_initial_value(initials, 0.)
        self.solver = ode(self.dynamics).set_integrator(self.backend, atol=1e2, rtol=0, nsteps=1).set_initial_value(initials, 0.)
    
    def dynamics(self, t, x):
        # print '%4f|'%(t), 
        # # for el in x:
        # #     print '%g, '%(el),
        # print ', '.join('%6.3f'%(el) for el in x)

        global im, ualbe_prev
        if True:
            U_isNotTVInput = im.u
        else:
            U_isNotTVInput = ualbe_prev + (t-timer.t)/im.Ts*(im.u-ualbe_prev)
        ## STEP ONE: collect all the states
        rcur = x[0:4]
        rpsi = dot(im.L,rcur)
        omega_r = x[4]
        # 其实MT系下仿真也是可以考虑的，见Ned2014书
        if rotor_frame == True:
            omega_g = omega_r
        else:
            omega_g = 0.
        theta_r = x[5]
        theta_g = x[6]
        
        ## STEP TWO: electromagnetic model
        # turn albe frame into rotor frame
        tmp_abc = im.fr.dqn2abc(np.append(U_isNotTVInput,0), 0.)
        ru = im.fr.abc2dqn(tmp_abc, theta_g)[0:2] # im.imt = im.fr.abc2dqn(tmp_abc, theta_psi_r)
        # get derivative of current
        em_derivative = dot(im.L_inv, np.append(ru,[0.,0.]) \
                             - dot(im.R + omega_g*im.G1 + (omega_g-omega_r)*im.G2, rcur) )
        
        ## STEP THREE: mechanical model
        im.Tem = im.mech_compensation * im.P/2 * (rpsi[0]*rcur[1]-rpsi[1]*rcur[0])
        mech_derivative = im.P/2./im.J * (im.Tem - im.Tl)

        # return np.append(em_derivative,[mech_derivative,omega_r,omega_g])

        ##########################################################################
        global ctrl
        est_CM  = x[7:9]
        tmp_psi = x[9:11]
        rs, alpha, omg = tuple(x[11:14])
        J = array([[0,-1],[1,0]],dtype='float64')
        U = U_isNotTVInput
        # turn rotor frame into albe frame
        tmp_abc = im.fr.dqn2abc(np.append(rcur[0:2],0), theta_g) # key is here!?
        I = im.fr.abc2dqn(tmp_abc, 0.)[0:2]
        ## use this to align the ui.dat 
        # print timer.t, U,  I
        # ERROR
        est_VM = im.Lr/im.Lm * (tmp_psi - im.sigma*im.Ls * I)
        est_err = est_VM - est_CM

        # omega_elec, or you can try ctrl.omega_syn
        if sqrt(tmp_psi[0]**2+tmp_psi[1]**2)>1e-4:
            E = (U - rs*I)
            # omega_elec = - (E[0]*-tmp_psi[1]+E[1]*tmp_psi[0])/sqrt(tmp_psi[0]**2+tmp_psi[1]**2) # 愚蠢的加了sqrt
            omega_elec = - (E[0]*-tmp_psi[1]+E[1]*tmp_psi[0])/(tmp_psi[0]**2+tmp_psi[1]**2)
            # derivative_CM = - alpha*(est_CM-im.Lm*I) + omg*dot(J,est_CM) + ctrl.mras.k_CM*est_err - 1.0*ctrl.omega_syn*dot(J,est_err)
        else:
            omega_elec = 0.0

        # VM
        derivative_VM = (U - rs*I) - ctrl.mras.k_VM*est_err
        # derivative_VM = (U - rs*I) - ctrl.mras.k_VM*est_err + 0.0*omega_elec*dot(J,est_err)
        
        # CM
        # derivative_CM = - alpha*(est_CM-im.Lm*I) + omg*dot(J,est_CM) + ctrl.mras.k_CM*est_err
        derivative_CM = - alpha*(est_CM-im.Lm*I) + omg*dot(J,est_CM) + ctrl.mras.k_CM*est_err - omega_elec*dot(J,est_err)

        global rampOrflat, turnOffRsAdapt
        # TODO: 这里应该直接修改Gamma，而不是修改Phi
        Phi = array([  [im.Lr/im.Lm*I[0], -(est_CM[0]-im.Lm*I[0]), -est_CM[1]],
                       [im.Lr/im.Lm*I[1], -(est_CM[1]-im.Lm*I[1]), est_CM[0]]], dtype='float64')
        Gamma = copy.deepcopy(ctrl.mras.Gamma_inv)
        if rampOrflat == True:
            Gamma[0][0] = 0.0
        else: # flat
            Gamma[1][1] = 0.0
            if turnOffRsAdapt == True:
                Gamma[0][0] = 0.0
        # update in alpha-beta frame
        ## No decouple of adaption
        thet_derivative = dot(Gamma, dot(Phi.T, est_err)) ###
        # thet_derivative[0] = 0.
        # thet_derivative[1] = 0.
        ##########################################
        ## update rs and alpha in MT frame
        est_err_MT = im.fr.albe2MT(est_err, est_CM[0], est_CM[1])
        IMT = im.fr.albe2MT(I, est_CM[0], est_CM[1])
        thet_derivative[0] = Gamma[0][0]*IMT[0]*est_err_MT[0]
        LmirMT = im.fr.albe2MT(-(est_CM-Lm*I), est_CM[0], est_CM[1])
        thet_derivative[1] = Gamma[1][1]*LmirMT[0]*est_err_MT[0]
        ##########################################

        return np.append( np.append(em_derivative,[mech_derivative,omega_r,omega_g]), \
                          np.append(np.append(derivative_CM, derivative_VM), \
                                    thet_derivative))
        # return np.append( np.append(em_derivative,[mech_derivative,omega_r,omega_g]), \
        #                   np.append(np.append(array([0.,0.]), array([0.,0.])), \
        #                             array([0,0,0],dtype='float64')))

    def update(s, t):
        s.solver.integrate(t, step=True)
        # s.solver.y[9:11] - im.Lr/im.Lm * im.sigma*im.Ls * s.solver.y[0:2] is only valid for stator frame

        return s.solver.y[0:4], s.solver.y[4], s.solver.y[5], s.solver.y[6], s.solver.y[7:9], \
                im.Lr/im.Lm * (s.solver.y[9:11] - im.sigma*im.Ls * s.solver.y[0:2]), \
                s.solver.y[11:14]

    # Techinique 1, with arrays:
    # n = len(x)
    def rKN(self, x, n, hs):
        k1, k2, k3, k4, xk = [], [], [], [], []

        fx = self.dynamics(timer.t, x)
        for i in range(n):
            k1.append(fx[i]*hs)
        for i in range(n):
            xk.append(x[i] + k1[i]*0.5)

        fx = self.dynamics(timer.t+hs/2., xk)
        for i in range(n):
            k2.append(fx[i]*hs)
        for i in range(n):
            xk[i] = x[i] + k2[i]*0.5

        fx = self.dynamics(timer.t+hs/2., xk)
        for i in range(n):
            k3.append(fx[i]*hs)
        for i in range(n):
            xk[i] = x[i] + k3[i]
    
        fx = self.dynamics(timer.t+hs, xk)
        for i in range(n):
            k4.append(fx[i]*hs)
        for i in range(n):
            x[i] = x[i] + (k1[i] + 2*(k2[i] + k3[i]) + k4[i])/6.
        return x

    def dP5(self, y0, t0, hs, maxiter=1): # maxiter=1 means fixed step
        hmax = hs
        hmin = hs/100. # 1e-6
        y, flag, _ = dopri54(self.dynamics, t0, y0, t0+hs, hmax, hmin, maxiter, 1e-5) 
        if flag != 0:
            print('@%6.4g sec, ode flat=%d.' %(t0,flag))
        return y 

class Timer(object):
    def __init__(self, t_ini, Ts):
        ''' discret fixed time interval
        '''
        self.t = t_ini
        self.Ts = Ts

    def ticktock(self):
        self.t += self.Ts
        if int((self.t/self.Ts)%(10000))==0:
            print('t=%g'%(t))

class Scope(object):
    ''' record data as a list for plotting
    '''
    id = 0
    def __init__(self, key_list, title='', printKey=False):
        self.time = []
        self.data = {}
        self.key_list = key_list
        for key in self.key_list:
            self.data[key] = []
        self.title = title

        # plot style
        plt.style.use('ggplot') 
        # plot setting
        mpl.rcParams['legend.fontsize'] = 12

        # print key or not
        self.printKey = printKey

        Scope.id += 1
        self.id = Scope.id

    def recordTime(self, t):
        self.time.append(t)
        
    def record(self, key, next_data):
        self.data[key].append(next_data)

    def plot(self, label_list):

        self.label_index = -1

        fig, axes = plt.subplots(ncols=1, nrows=len(self.data));
        fig.suptitle(self.title)
        fig.subplots_adjust(right=0.85, bottom=0.1, top=0.95, hspace=0.25)
        ax_list = axes.ravel()

        ax_index = -1
        for key in self.key_list: # don't use "for key in self.data:". it returns key in random order
            ax_index += 1
            ax = ax_list[ax_index]
            ax.set_xlabel(r'$time / sec$', fontsize=15)
            #ax.set_ylabel('y')

            # how much lines need to be plot for the key?
            to_plot_list = []
            for i in range(len(self.data[key][0])):
                to_plot_list.append([el[i] for el in self.data[key]])

            for i in range(len(to_plot_list)):
                self.label_index += 1
                ax.plot(self.time, to_plot_list[i], label=label_list[self.label_index])
                if self.printKey == True:
                    print(key)
            #[left, bottom, width, height]
            ax.legend(bbox_to_anchor=(1.08,0.5), borderaxespad=0., loc='center', shadow=True) #http://matplotlib.org/1.3.1/users/legend_guide.html
            
            #ax.set_frame_on(True) # make it transparent

    def locus(self, key, target, xlabel, ylabel):
        ''' for instance, to plot locus of id & iq.
            key = 'idqns'            
            target = [0,1], 0 for x, 1 for y
        '''
        x = [el[target[0]] for el in self.data[key]]
        y = [el[target[1]] for el in self.data[key]]
        
        fig = figure('locus');
        fig.suptitle('locus')
        ax = fig.gca()
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.plot(x,y)

    def getPlotData(self, key, which_one):
        return self.time, [el[which_one] for el in self.data[key]]

    def getData(self, key, which_one):
        if which_one == 'all':
            lTmp = []
            print('\t', len(self.data[key][0]))
            for i in range(len(self.data[key][0])):
                lTmp.append([el[i] for el in self.data[key]])
            return lTmp            
        else:
            return [el[which_one] for el in self.data[key]]

    def saveAsHDF5(self, where=None):
        if where == None:
            f_name = './h5/%d_%s.h5'%(self.id,self.title[0:4])
        else:
            f_name = './h5/%s/%d_%s.h5'%(where,self.id,self.title[0:4])
            if not os.path.exists('./h5/%s/'%(where)):
                os.makedirs('./h5/%s/'%(where))
        print('save', self.title)
        open(f_name,'w').close() # create file
        fh5 = h5py.File(f_name,'w') # h5py won't create directory for you. you have to create yourself.
        for key in self.key_list:
            print('\t', key)
            fh5.create_dataset(key, data=self.getData(key, 'all'))
        fh5.close()

class Frame2(object):
    def __init__(self, trans_type, theta_g_ini):
        # transfer type
        self.trans_type = trans_type

        # the position of the gyrate frame.
        # theta_g === <d,a> according to matrix MT
        self.theta_g = theta_g_ini
        # self.omega_g = 0. 

        # constant of transformation
        if self.trans_type == 'constant_power':
            # aka power invariant
            self.CT = sqrt(2./3.) # (2-49~51)
        elif self.trans_type == 'constant_magnitude':
            # aka not power invariant
            self.CT = 2./3. # (2-58~60)
        else:
            print('WRONG TRANSFER TYPE')

    def dqn2abc(s, dqn, theta_g):
        ''' Let theta_g be <d,a>, that is when theta_g = 0, d-q frame <=> alpha-beta frame.
            The inverse transformation matrix is the transpose of the forward one in power invariant.
        '''
        # matrix of transformation by Lipo@p88 (2-49~51) (2-58~60)
        # but rotage positive 90 degree
        MT = array([[cos(theta_g),        -sin(theta_g),        1./sqrt(2)],
                    [cos(theta_g-2*pi/3), -sin(theta_g-2*pi/3), 1./sqrt(2)],
                    [cos(theta_g-4*pi/3), -sin(theta_g-4*pi/3), 1./sqrt(2)]], dtype='float64')

        # get abc variables
        if s.trans_type == 'constant_power':
            abc = s.CT * dot(MT, dqn)
        else:
            abc = 1 * dot(MT, dqn)
        return abc

    def abc2dqn(s, abc, theta_g):
        ''' refer to dqn2abc()
        '''
        # matrix of transformation by Lipo@p88 (2-49~51) (2-58~60)
        # but rotage positive 90 degree
        MT = array([[ cos(theta_g),  cos(theta_g-2./3.*pi),  cos(theta_g-4./3.*pi)],
                    [-sin(theta_g), -sin(theta_g-2./3.*pi), -sin(theta_g-4./3.*pi)],
                    [    1./sqrt(2),              1./sqrt(2),              1./sqrt(2)]], dtype='float64')

        # get dqn variables
        dqn = s.CT * dot(MT, abc)
        return dqn

    def albe2MT(s, albe, psiral, psirbe):

        amplitute = sqrt(psiral**2 + psirbe**2)
        if amplitute <1e-4:
            cosThetaM = 1
            sinThetaM = 0
        else:
            cosThetaM = psiral / amplitute
            sinThetaM = psirbe / amplitute

        # matrix of transformation
        MT = array([[ cosThetaM, sinThetaM],
                    [-sinThetaM, cosThetaM]], dtype='float64')

        # get Mag/Torque variables
        magtorq = dot(MT, albe)
        return magtorq

    def MT2albe(s, magtorq, psiral, psirbe):

        amplitute = sqrt(psiral**2 + psirbe**2)
        if amplitute == 0.:
            amplitute = 1e-3
            psiral = amplitute
        cosThetaM = psiral / amplitute
        sinThetaM = psirbe / amplitute

        # matrix of transformation
        MT = array([[ cosThetaM, -sinThetaM],
                    [ sinThetaM, cosThetaM]], dtype='float64')

        # get Mag/Torque variables
        albe = dot(MT, magtorq)
        return albe

    def albe2rotor(s, albe, psiral, psirbe, slip):
        ''' 变换到转子坐标系下的旋转变换，变换的角度采用转子磁链的估计减去给定的转差频率决定。
            slip的单位是rad/s
        '''
        amplitute = sqrt(psiral**2 + psirbe**2)
        if amplitute == 0.:
            amplitute = 1e-3
            psiral = amplitute
        cosThetaM = psiral / amplitute
        sinThetaM = psirbe / amplitute
        cos_angleSynMinusSlip = cosThetaM*cos(slip) + sinThetaM*sin(slip)
        sin_angleSynMinusSlip = sinThetaM*cos(slip) - cosThetaM*sin(slip)
        cos_angleSynMinusSlip = cos(im.omega_r)
        sin_angleSynMinusSlip = sin(im.omega_r)
        # matrix of transformation
        MT = array([[ cos_angleSynMinusSlip, sin_angleSynMinusSlip],
                    [-sin_angleSynMinusSlip, cos_angleSynMinusSlip]], dtype='float64')

        # get dq variables
        dq = dot(MT, albe)
        return dq

class ThreePhaseInductionMachine(object):
    """ This simulation model uses ids, iqs, psi_dr, psi_qr as states.
        And the fr must be a stationary one, i.e. the simulation equation is only valid for stationary frame im.

        Future Modify:
        self.p = {’m’: 1.0, ’b’: 0.7, ’c’: 5.0, ’func’: ’y’,
                ’A’: 5.0, ’w’: 2*pi, ’y0’: 0.2,
                ’tstop’: 30.0, ’dt’: 0.05}
        self.p.update(kwargs)
    """ 
    def __init__(self, Ts, \
                 rs,  rr, \
                 Lm, Lls, Llr, \
                 fr, J, P, initials):
        # sample time
        self.Ts = Ts

        # frame
        self.fr = fr

        ## states are [ids, iqs, psi_dr, psi_qr] in alpha-beta frame, and omega_r.
        self.st = array([0,0,0,0], dtype='float64')
        # self.st = array([1e-3,1e-3,1e-5,1e-5], dtype='float64')
        # omega_r = 2 * P * omega_rm
        self.omega_r = 0.
        self.rpm = 0.
        # ds axis coincides with as axis
        self.theta_r = 0.
        # just flux linkage. not hybrid flux linkage, nor unit volt.
        # self.psi = array([0,0,0,0,0,0], dtype='float64')
        
        ## params.
        # poles
        self.P = P
        # inertia of rotating
        self.J = J
        # stator & rotor resistance
        self.rs = rs
        self.rr = rr
        # stator & rotor leakage
        self.Lls = Lls
        self.Llr = Llr
        # direct axis & quadrature axis magnetizing inductance
        self.Lm = Lm
        # stator * rotor self inductance
        self.Ls = Lls + Lm
        self.Lr = Llr + Lm
        # viscous friction is neglected for now
        self.viscous = 0.

        # be an im
        self.sigma = 1. - self.Lm**2./(self.Ls*self.Lr)
        self.Tr = self.Lr/self.rr

        # just for convenience
        self.k1 = float(self.Lm/self.Tr)
        self.k2 = float(self.Lm/(self.sigma*self.Ls*self.Lr))
        self.k3 = 1./(self.sigma*self.Ls)

        # equivalent inductance
        self.Leq = self.Lm**2/self.Lr
        self.alpha = 1.0/self.Tr
        # equivalent rotor resistance
        self.rreq = self.Leq*self.alpha
        # leakage inductance
        self.Lsigma = self.Ls * self.sigma

        # generate matrix form
        # self.genMatrices()

        # mechanical compensation due to transformation
        if self.fr.trans_type == 'constant_power':
            self.mech_compensation = 1.
        elif self.fr.trans_type == 'constant_magnitude':
            self.mech_compensation = 1.5
        else:
            print('WRONG TRANSFER TYPE')

        # using is, ir as states
        self.rst = array([0,0,0,0],dtype='float64')
        self.rcur = array([0,0,0,0],dtype='float64')
        self.rpsi = array([0,0,0,0],dtype='float64')
        self.genMatrices2()

        self.psi_oneStepAhead = array([0,0],dtype='float64')
        self.tmp = copy.deepcopy(self.st)
        self.tmp_derivative = array([0,0],dtype='float64')

        # all in one
        self.x = initials

    def genMatrices2(s):
        # Reactance
        s.L = array([[s.Ls,          0,      s.Lm,           0],
                     [0,             s.Ls,      0,        s.Lm],
                     [s.Lm,          0,      s.Lr,           0],
                     [0,             s.Lm,      0,        s.Lr]], dtype='float64')
        s.L_inv = np.linalg.inv(s.L)

        # Resistance
        s.R = diag([s.rs, s.rs, s.rr, s.rr])

        # Gyrate matrix without omega_r which will be added back in em_derivative
        # G = omega_r J L

        # for rotor frame
        s.G1 = array([[0,     -s.Ls,       0,      -s.Lm],
                      [s.Ls,     0,        s.Lm,      0],
                      [0,        0,        0,         0],
                      [0,        0,        0,         0]], dtype='float64')
        # for albe frame
        s.G2 = array([[0,     0,        0,      0],
                      [0,     0,        0,      0],
                      [0,    -s.Lm,     0,  -s.Lr],
                      [s.Lm,  0,     s.Lr,      0]], dtype='float64')
        # s.G2 = array([[0,     0,        0,      0],
        #              [0,     0,        0,      0],
        #              [0,    s.Lm,      0,   s.Lr],
        #              [-s.Lm, 0,     -s.Lr,     0]], dtype='float64')

    def updateMatrices2(s, t):
        # # if timer.t>1.8 and timer.t<1.8+0.5:
        # if isclose(t, 3.5):            
        #     s.rs = s.rs + 0.3
        #     s.rr = s.rr + 0.2 # 为什么改成0.1反而发散了？
        #     s.Tr = s.Lr/s.rr
        #     s.k1 = float(s.Lm/s.Tr)
        #     global mras_which, ctrl
        #     if mras_which == 'separateModel':
        #         ctrl.mras.thet_real[1] = 1./s.Tr
        #         ctrl.mras.thet_real[0] = s.rs
        #     # Resistance
        #     s.R = diag([s.rs, s.rs, s.rr, s.rr])

        pass
        # ## change in simulated Lm
        # if isclose(t, 0.1):
        #     s.Lm += 0.005
        #     # stator * rotor self inductance
        #     s.Ls = s.Lls + s.Lm + 0.01*0
        #     s.Lr = s.Llr + s.Lm + 0.01*0

        #     # be an im
        #     s.sigma = 1. - s.Lm**2./(s.Ls*s.Lr)
        #     s.Tr = s.Lr/s.rr

        #     # just for convenience
        #     s.k1 = float(s.Lm/s.Tr)
        #     s.k2 = float(s.Lm/(s.sigma*s.Ls*s.Lr))
        #     s.k3 = 1./(s.sigma*s.Ls)

        #     s.genMatrices2() # 不调用这个的效果和实验中出现的rs在不同磁链幅值下不同收敛值的现象一致。

        #     ## use incorrect value in tdda observer
        #     # s.Lm *= 0.90
        #     # s.Lm *= 0.99
        #     # s.Lm *= 1.01
        #     # s.Lm *= 1.1 # 大对小，1.8~2.6
        #     print s.Lm
        #     s.Lm *= 1.01
        #     print s.Lm
        #     # stator * rotor self inductance
        #     s.Ls = s.Lls + s.Lm + 0.01*0
        #     s.Lr = s.Llr + s.Lm + 0.01*0

        #     # be an im
        #     s.sigma = 1. - s.Lm**2./(s.Ls*s.Lr)
        #     s.Tr = s.Lr/s.rr

        #     # just for convenience
        #     s.k1 = float(s.Lm/s.Tr)
        #     s.k2 = float(s.Lm/(s.sigma*s.Ls*s.Lr))
        #     s.k3 = 1./(s.sigma*s.Ls)

        pass
        # # change in Lm used in observer
        # if isclose(t, 0.001):
        #     s.Lm -= 0.02
        #     # stator * rotor self inductance
        #     s.Ls = s.Lls + s.Lm + 0.01*0
        #     s.Lr = s.Llr + s.Lm + 0.01*0

        #     # be an im
        #     s.sigma = 1. - s.Lm**2./(s.Ls*s.Lr)
        #     s.Tr = s.Lr/s.rr

    def update_accurately(self, u, Tl):
        '''dqn frame with rotating speed omega_g'''

        # to be called in IMDynamic
        self.u = u
        self.Tl = Tl

        ## STEP ONE: update time-varing params of IMs
        global timer, sim
        self.updateMatrices2(timer.t)
        if False:
            # if not sim.solver.successful():
            #     print 'Fail @ %g' % timer.t
            # print sim.solver.t, timer.t
            self.rcur, self.omega_r, self.theta_r, self.fr.theta_g, ctrl.mras.est_CM, ctrl.mras.est_VM, ctrl.mras.thet = sim.update(timer.t)
        else:
            # RK4
            self.x = array(sim.rKN(self.x, len(self.x), im.Ts))
            x = self.x
            self.rcur=x[0:4]; self.omega_r=x[4]; self.theta_r=x[5]; self.fr.theta_g=x[6]; 
            ctrl.mras.est_CM=x[7:9]; ctrl.mras.est_VM=im.Lr/im.Lm*(x[9:11] - im.sigma*im.Ls * x[0:2]); ctrl.mras.thet=x[11:14]

            # # DP5
            # if timer.t == Ts:
            #     print timer.t, 'initialize!'
            #     x=initials
            # else:
            #     x = self.rcur.tolist()+[self.omega_r,self.theta_r,self.fr.theta_g]+ctrl.mras.est_CM.tolist()+ctrl.mras.est_VM.tolist()+ctrl.mras.thet.tolist()
            # x = sim.dP5(x, timer.t, im.Ts)
            # self.rcur=x[0:4]; self.omega_r=x[4]; self.theta_r=x[5]; self.fr.theta_g=x[6]; ctrl.mras.est_CM=x[7:9]; ctrl.mras.est_VM=x[9:11]; ctrl.mras.thet=x[11:14]

        # get flux linakage
        self.rpsi = dot(self.L,self.rcur)
        # get electric states
        self.rst = np.append(self.rcur[0:2],self.rpsi[2:4])
        # moderate mech quantities
        self.rpm = self.omega_r / 2. / pi * 60 / self.P * 2. 
        self.theta_r = self.theta_r - self.theta_r//(2*pi)*(2*pi)

        ######################
        # turn rotor frame into albe frame before update frame theta
        tmp_abc = self.fr.dqn2abc(np.append(self.rst[0:2],0), self.fr.theta_g)
        self.st[0:2] = self.fr.abc2dqn(tmp_abc, 0)[0:2] 
        tmp_abc = self.fr.dqn2abc(np.append(self.rst[2:4],0), self.fr.theta_g)
        self.st[2:4] = self.fr.abc2dqn(tmp_abc, 0)[0:2]

        self.is_MT = self.fr.albe2MT(self.rst[0:2], self.rst[2], self.rst[3])
        self.psi_MT = self.fr.albe2MT(self.rst[2:4], self.rst[2], self.rst[3])
        self.ir_MT = (self.psi_MT - self.Lm*self.is_MT)/self.Lr
        # self.ir_MT = self.fr.albe2MT(self.cur[2:4], self.rst[2], self.rst[3])

class Controller2(object):

    def __init__(self, acm, timer):
        ''' The argument acm & timer are reference of the variable sm & timer.
            By reference it means if sm or timer is modified outside the Controller class,
            the value of the ctrl.acm or ctrl.timer will coincide with the outside one(s).
        '''
        # target to control
        self.acm = acm
        self.timer = timer

        # vvvf relevent
        self.f_synchro = 0.
        self.f_increse_step = 0.05
        # linear vvvf coefficient = rated phase magnitude value / freq corresponding to rated rpm
        self.cvf0 = acm.rated_volt_mag / acm.rated_freq

        # FOC rotor flux estimator wi VCM
        self.tmp_alpha = 0.
        self.tmp_beta = 0.
        self.psi_r_alpha_est = 0.
        self.psi_r_beta_est = 0.
        self.theta_psi_r_est = 0.000
        self.mag_psi_r_est = 1e-3
        # FOC rotor flux estimator wi CSM
        self.psi_CM_r_alpha_est = 0.
        self.psi_CM_r_beta_est = 0.
        self.theta_psi_CM_r_est = 0.000
        self.mag_psi_CM_r_est = 1e-3
        
        # IFOC relevent
        self.theta_M = 0. # had to be defined
        self.omega_syn = 0.

        # for Feedback Control
        self.umtn = zeros(3)
        # self.hc = hysteresisComparator()
        self.i_mag = 0.

        # for PI Control the Integration Register
        self.tmp_IR = {}
        self.tmp_IR['torque'] = 0.
        self.tmp_IR['field'] = 0.
        self.tmp_IR['ud'] = 0.
        self.tmp_IR['uq'] = 0.

        # for MRAS
        global mras_which
        if mras_which == 'separateModel':
            ## Tr & Omg
            self.mras = MRAS_RotorFluxCurrentSpeedModel_LmKnown(self.acm.Ts)
            ## rs
            # self.identifier4rs = TimeScaleSeparationIdentifierOfStatorResistance(self.acm.Ts)
            ## STA
            # self.sta = SuperTwistingAlgorithm(self.acm.Ts)

        # cnt for TDDA
        self.cntVS = int(0.15/1e-4) + 1
        self.cntFlat = 0 

    def PI_returnIPart(self,error,P,I,key,I_lim):
        # self.tmp_IR[key] += I*error # wrong
        self.tmp_IR[key] += self.acm.Ts*I*error
        # 积分限幅是必需的
        if self.tmp_IR[key] > I_lim:
            self.tmp_IR[key] = I_lim
        elif self.tmp_IR[key] < -I_lim:
            self.tmp_IR[key] = -I_lim

        ## for tuning I coefficient
        # global timer
        # if int((timer.t-self.acm.Ts)/0.1) != int(timer.t/0.1):
        #     if key == 'torque':
        #         print timer.t, key, 'I_part', self.tmp_IR[key]

        if key == 'torque':
            return P*error + self.tmp_IR[key], self.tmp_IR[key]
        else:
            return P*error + self.tmp_IR[key]
    
    def PI(self,error,P,I,key,I_lim):
        # self.tmp_IR[key] += I*error # wrong
        self.tmp_IR[key] += self.acm.Ts*I*error
        # 积分限幅是必需的
        if self.tmp_IR[key] > I_lim:
            self.tmp_IR[key] = I_lim
        elif self.tmp_IR[key] < -I_lim:
            self.tmp_IR[key] = -I_lim

        return P*error + self.tmp_IR[key]

    def P(self,error,P):
        return P*error

    def givenOmegarReturnPsir(self, rpm_cmd):
        ''' What the argument should be, omega_r or omega_r_cmd?
        '''
        global rampOrflat, turnOffRsAdapt
        global lstInstant, lstAmpl

        self.cntVS += 1 # cnt of variable structure
                        # the value of this can determine which structure comes first
        noBreak = True
        while True:
            for i in range(len(lstAmpl)):
                if self.cntVS*im.Ts<lstInstant[i]:
                    noBreak = False
                    deriv_ExtraModulus = (lstAmpl[i]-lstAmpl[i-1])/.15
                    if isclose(deriv_ExtraModulus,0.):
                        rampOrflat = False 
                        self.cntFlat += 1
                        if self.cntFlat*im.Ts<0.6: # 1.2 /2
                            turnOffRsAdapt = True 
                        else:
                            turnOffRsAdapt = False
                    else: 
                        rampOrflat = True
                        self.cntFlat = 0
                    extraModulus = lstAmpl[i] - (lstInstant[i] - self.cntVS*im.Ts)*deriv_ExtraModulus # i=0时，循环到末尾lstAmpl[-1]
                    break
            if noBreak==False:
                break
            print('cntVS clean up at', self.cntVS)
            self.cntVS = 1

        global timer, scp_tdda
        scp_tdda.recordTime(timer.t)
        scp_tdda.record('A', array([extraModulus]))
        scp_tdda.record('B', array([deriv_ExtraModulus]))
        scp_tdda.record('C', array([rampOrflat, turnOffRsAdapt]))

        alpha = 1./im.Tr
        return 1.0 + extraModulus + deriv_ExtraModulus/alpha
        # return 1.0 + extraModulus

    def feedback_control(self, rpm_cmd, Omg):

        acm = self.acm

        # omega_r_cmd = rpm_cmd / 60. * acm.P/2 * 2*pi
        omega_r_cmd = rpm_cmd / 60. * acm.P * pi
        
        # 转速环、磁链查表
        # 根据转速给定查表磁链给定，除以电感，确定M轴电流给定。
        # 根据转速给定和转速反馈之间的误差，经过PI，得到转矩给定，除以转矩因数，得到T轴电流给定。
        self.imt_cmd = zeros(2)

        # #0
        # magnitude of psi_r command
        self.magPsi_SteadyStateCmd = self.givenOmegarReturnPsir(rpm_cmd)
        # global t
        # if t>0.2: # t=0.7 转速突变
        #     # self.magPsi_SteadyStateCmd += 0.2 * cos(2.*pi*1. * t)
        #     self.magPsi_SteadyStateCmd += 0.1 * cos(2.*pi*5. * t)

        if sensorless == True:
            # magPsi_fdbk = sqrt(self.mras.est_VM[0]**2+self.mras.est_VM[1]**2)
            # self.magPsi_PICmd = self.PI(self.magPsi_SteadyStateCmd-magPsi_fdbk, 10., 0.05*1e4, 'field', 20)
            # no need of doing this, tune PI for u_M is enough
            self.magPsi_PICmd = self.magPsi_SteadyStateCmd
        else:
            self.magPsi_PICmd = self.magPsi_SteadyStateCmd
            # magPsi_fdbk = sqrt(im.st[2]**2+im.st[3]**2)
            # self.magPsi_PICmd = self.PI(self.magPsi_SteadyStateCmd-magPsi_fdbk, 10., 0.05*1e4, 'field', 20)
        # Excitation for satisfyi PE condition
        # if self.hc.comp(omega_r_cmd-Omg, 1. / 60. * acm.P * pi, 10. / 60. * acm.P * pi):
            # self.magPsi_cmd = self.magPsi_cmd + 0.5 * cos(2.*pi*10. * t)
            # self.magPsi_cmd = self.magPsi_cmd + 0.2 * cos(2.*pi*.5 * t)
            # pass

        # ''' 逻辑是这样的，为了让实际磁链跟踪给定的磁链，加入PI环，实际磁链能跟踪给定的磁链。
        #     但是由于给定iMs时是根据iMs*=PI(psi_err)/Lm执行的，当iMs跟踪iMs*时，满足psi-Lm*iMs=0的磁链值是magPsi_PICmd，
        #     而不是我们期望的magPsi_SteadyStateCmd（=magPsi_fdbk），所以转子电流在M轴上有值。
        # '''
        # Lr*iMr = psi_Mr - Lm*iMs
        # self.imt_cmd[0] = self.magPsi_SteadyStateCmd/acm.Lm
        self.imt_cmd[0] = self.magPsi_PICmd/acm.Lm
        # self.imt_cmd[0] = (self.magPsi_PICmd - acm.Lr*acm.ir_MT[0])/acm.Lm
        
        # self.i_mag += acm.Ts / acm.Tr * (-self.i_mag + self.imt[0])
        # self.imt_cmd[0] = self.PI(self.magPsi_SteadyStateCmd-acm.Lm*self.i_mag, 10., 0.05*1e4, 'field', 20)
        # self.magPsi_PICmd = self.imt_cmd[0] * acm.Lm

        # self.torque_cmd, self.torque_cmd_Ipart = self.PI_returnIPart(omega_r_cmd-Omg, 2., 0.05*1e4, 'torque', 90.)
        self.torque_cmd, self.torque_cmd_Ipart = self.PI_returnIPart(omega_r_cmd-Omg, 10., 0.05*1e4, 'torque', 90.)
        # self.torque_cmd = self.PI(omega_r_cmd-Omg, 6, 0., 'torque', 90)
        # self.torque_cmd = self.PI(omega_r_cmd-Omg, 1, 0., 'torque', 90)

        if self.torque_cmd > 90:
            self.torque_cmd = 90
        elif self.torque_cmd < -90:
            self.torque_cmd = -90

        # note self.magPsi_SteadyStateCmd is the command value, see 《电机控制》@p146
        self.imt_cmd[1] = self.torque_cmd / self.magPsi_SteadyStateCmd / (acm.P/2) * acm.Lr/acm.Lm 

        if self.imt_cmd[1] > 50:
            self.imt_cmd[1] = 50
        elif self.imt_cmd[1] < -50:
            self.imt_cmd[1] = -50

        # 电流环，根据电流误差确定MT系电压给定，经过观测得的MT系位置，得到abc系下电压给定。
        self.umtn[0] = self.PI(self.imt_cmd[0]-self.imt[0], 50, 0.05*1e4, 'ud', 500)
        self.umtn[1] = self.PI(self.imt_cmd[1]-self.imt[1], 50, 0.05*1e4, 'uq', 500)
        # just limit to the mag
        if self.umtn[0]>500:
            self.umtn[0]=500
        elif self.umtn[0]<-500:
            self.umtn[0]=-500
        if self.umtn[1]>800:
            self.umtn[1]=800
        elif self.umtn[1]<-800:
            self.umtn[1]=-800
        # # if saturation then shrink by propotional
        # mag = sqrt(self.umtn[0]**2 + self.umtn[1]**2)
        # lim = 500
        # if mag > lim:
        #     self.umtn[0] = self.umtn[0] / mag * lim
        #     self.umtn[1] = self.umtn[1] / mag * lim

        # if(Vs_req>Vs_max)
        # {
        #     Vd1=(1.0*Vs_max/Vs_req)*Vd1;
        #     Vq1=(1.0*Vs_max/Vs_req)*Vq1;
        # }

    def ifoc(self, rpm_cmd, ialbe, Omg, Tr):
        ''' [INTRO]     Indirect Field Oriented Control
                        A.k.a slip freqency control
            [FEATURE]   Using PI controller for speed & current loop both, which means it needs feedback of the imt.
                        Using no decoule of the voltage & current.
            [Arguments] ialbe is the measured currents
                        Omg & Tr can be identified values or from im
        '''
        acm = self.acm

        global im
        
        self.theta_M = arctan2(im.st[3],im.st[2]) #//////////

        # how to get imt??? ZLH: just use theta_M (before be calc'd)
        ialbe = array(ialbe.tolist()+[0.], dtype='float64')
        tmp_abc = acm.fr.dqn2abc(ialbe, acm.fr.theta_g)
        self.imt = acm.fr.abc2dqn(tmp_abc, self.theta_M)

        # the feedback control for speed & mt current
        self.feedback_control(rpm_cmd, Omg)

        # # decouple of voltage & current 自己推导看看
        # # note that we use measured mt current val not cmd val
        # self.umtn[0] -= self.omega_syn*acm.sigma*acm.Ls*self.imt[1]
        # self.umtn[1] += self.omega_syn*acm.sigma*acm.Ls*self.imt[0]

        # mt to dqn frame
        self.umtn[2] = 0 # just for alignment in acm.fr.dqn2abc()
        tmp_abc = acm.fr.dqn2abc(self.umtn, self.theta_M)
        udqns = acm.fr.abc2dqn(tmp_abc, acm.fr.theta_g)

        # to consistent with the model acm
        ualbe = udqns[0:2]

        # omega_sl_cmd == its_cmd/acm.Tr/ims_cmd
        # NOTE THIS EXPRESSION ONLY VALID WHEN PSI_R IS STEADY
        # AND I AM CONFUSED HOW COME iM* = psi_r/Lm EVEN WHEN STEADY NO LONGER
        # Update theta_M before calc. iMT or after???
        self.omega_sl = self.imt_cmd[1]/Tr/self.imt_cmd[0]
        self.omega_syn = Omg + self.omega_sl
        # theta_psi_r is approximated to theta_M
        self.theta_M += acm.Ts * self.omega_syn
        if self.theta_M > pi:
            self.theta_M -= 2*pi    

        return ualbe

    def foc_VM(self, rpm_cmd, ialbe, Omg, psial, psibe):
        ''' [INTRO]     Direct Field Oriented Control
                        A.k.a. vector control
                        The key of FOC is to get the position of the rotor flux vector, 
                        and the left things are no more than PI control of speed & current loop.
            [FEATURE]   Using PI controller for speed & current loop both, which means it needs feedback of the imt.
                        Using no decoule of the voltage & current.
            [Arguments] ialbe is the measured currents
                        psial, psibe can be observed values or from im
        '''
        # global theta_psi_r
        acm = self.acm

        # 根据测得的电流ialbe，和观测得的转子磁链，将电流变换到MT系下。
        self.imt = acm.fr.albe2MT(ialbe, psial, psibe)

        # the feedback control for speed & mt current
        self.feedback_control(rpm_cmd, Omg)

        # mt to dqn frame
        ualbe = acm.fr.MT2albe(self.umtn[0:2], psial, psibe)
        return ualbe

    def certainFreqExcite(self, Vs, f_e):
        omega_e = 2.*pi*f_e
        ea = Vs * cos(omega_e * t)
        eb = Vs * cos(omega_e * t - 2./3.*pi)
        ec = Vs * cos(omega_e * t - 4./3.*pi)
        u_line2line = array([ea-ec, eb-ec, 0], dtype = 'float64')
        ualbe = im.fr.abc2dqn(u_line2line, im.fr.theta_g)
        # neutral 2 ground voltage is actually known with Y connection assumption
        u_n2g = 1./3.*(ea+eb+ec) # Lipo@(2.122)
        # ucs is actually ucn here. Lipo use 's' for neutral point as well as 'stator'.
        ucs = ec - u_n2g # Lipo@(2.118)
        ualbe[2] += im.fr.CT * 3./sqrt(2.) * ucs
        return ualbe[0:2]

class MRAS_RotorFluxCurrentSpeedModel_LmKnown(object):
    ''' 
    '''
    def __init__(self, Ts):
        self.Ts = Ts
        global initials
        self.est_CM = array([initials[-7],initials[-6]])
        self.tmp_psi = array([initials[-5],initials[-4]])
        self.est_VM = self.tmp_psi - im.Lr/im.Lm * im.sigma*im.Ls * array([initials[0],initials[1]])
        self.est_err = zeros(2)
        
        # self.Lr = im.Lr
        # self.Lm = im.Lm
        # self.thet = array([1, 5, 0], dtype='float64')
        self.thet = array([1, 7, 0], dtype='float64')
        self.thet_real = array([im.rs, 1./im.Tr, im.omega_r], dtype='float64')
        print('thet_real:', self.thet_real, 'the third is a time-varing one.')

        if sensorless:
            self.Gamma_inv = diag([500,1e5,2e5]) # for sensorless case, alpha 2e5 for alpha leads to instability
            self.Gamma_inv = diag([200,5e4,2e5]) # for sensorless case, 2e5 for alpha leads to instability
            self.Gamma_inv = diag([50,5e3,3e5]) # for sensorless case,
            self.Gamma_inv = diag([50,2.5e5,3e5]) # for sensorless case,
        else:
            self.Gamma_inv = diag([5000,1e6,2e6]) # for sensed case
            self.Gamma_inv = diag([100,5e4,2e6]) # for sensed case and a lot of speed changes

        self.k_VM = 200
        self.k_CM = 200
        if sensorless:
            self.k_VM = 500 # gamma_alpha = 1e5 won't diverge
            self.k_CM = 500
            self.k_VM = 800 # gamma_alpha = 5e5 won't diverge
            self.k_CM = 800

        # self.k_VM=0.
        # self.k_CM=0.

    def update(self, U, I):
        ''' U is from ualbe at t=k
            I is from the tmp_st, at t=k,
            and PSI is from tmp_psi_ole at t=k.
            Now t=k+1.
        '''
        ## time-varing parameters updates its real value here
        # global tmp_omega_r
        # self.thet_real[2] = tmp_omega_r
        # self.thet = copy.deepcopy(self.thet_real)
        # global im, timer
        # self.U=U
        # self.I=I

        ## Step 1: compute regressor Phi at t=k
        

        ## Step 2: observer output est at t=k+1 using quantites at t=k
        # self.est_CM, self.est_VM = EstDynamic(timer.t)
        self.est_err = self.est_VM - self.est_CM

        self.mag_psi_ole = sqrt(dot(self.est_VM, self.est_VM))
        if isclose(self.mag_psi_ole, 0.):
            self.phase_psi_ole = 0.
        else:
            self.phase_psi_ole = self.est_VM[0] / self.mag_psi_ole

        ## Step 3: We use quantities (including est_err) at t=k+1 to get thet at t=k+1

        # self.thet_real[2] = im.omega_r
        # self.thet = copy.deepcopy(self.thet_real)
        # self.thet[0] = copy.deepcopy(self.thet_real[0])
        # self.thet[1] = copy.deepcopy(self.thet_real[1])
        # self.thet[2] = copy.deepcopy(self.thet_real[2])
        # PI adpation law for speed. 
        # self.thet[1] += 20 * thet_derivative[1] / self.Gamma_inv[1][1]

        ## Step 4: error computing
        self.thet_err = self.thet_real - self.thet

        ## move this to other place
        ## See the observability between Omg and Tr
        # self.mag = sqrt(self.est_VM[0]**2 + self.est_VM[1]**2)
        self.mag = sqrt(self.est_CM[0]**2 + self.est_CM[1]**2)

if __name__ == '__main__':

    controlStrategy_which = 'focVM_woDecouple'
    controlStrategy_which = 'ifoc_woDecouple' # Note ifoc needs Tr. I.e., an initial value of Tr is needed
    print('controller:', controlStrategy_which)
    mras_which = 'separateModel'
    rotor_locked = False
    rotor_frame = False
    sensorless = True
    extra = False # True only valid for im.update()
    if sensorless:
        print('Speed sensorless')
    else:
        print('Speed sensed')
            

    print('Independent resistances adaptation is already on.')

    ## Simulation arguments.
    Tctrl = 1e-4
    Ts = 5e-6
    Ts = 1e-4
    cntSim = 0
    cntCtrl = 0
    delay_time = 0. #
    t_max = 120 + delay_time
    t_max = 12*2 + delay_time
    # t_max = 9 + delay_time
    # t_max = 0.001 + delay_time
    t_max = 8 + delay_time
    print('t_max =', t_max)
    timer = Timer(0., Ts)

    sim = Simulator()

    if True:
        ## initialize dqn frame arguments to determine dqn frame's position @ t=0s
        theta_g_ini = 0. 
        frame_albe = Frame2('constant_power', theta_g_ini)

        ## im: non-salient pole & symmetrical damper winding
        print("XiaoBo's machine")
        rs = 1.1; rr = 1.3
        Lm = 0.15
        Lls = 0.0064; Llr = 0.0125
        J = 0.06
        P = 3 * 2
        # print 1.6/(Lm+Llr)
        # quit()
        # print "Marino's machine"
        # rs = 5.3; rr = 3.3
        # Lm = 0.34
        # Lls = 0.365-Lm; Llr = 0.375-Lm
        # J = 0.075
        # P = 3 * 2
        # print "KongBo's machine"
        # rs = 11e-3; rr = 8.1e-3
        # Lm = 1.63e-3;
        # Lls = 0.0064e-2; Llr = 0.0125e-2
        # J = 0.06
        # P = 2 * 2 
        # print "JinHai's machine"
        # rs = 1.4; rr = 0.8
        # Lm = 0.135
        # Lls = 0.005; Llr = 0.005
        # J = 0.078
        # P = 2 * 2
        im = ThreePhaseInductionMachine( Ts, \
                                         rs, rr, \
                                         Lm, Lls, Llr, \
                                         frame_albe, J, P, initials)
        ## Controller 
        im.rated_volt_mag = 220.*sqrt(2) # for vvvf
        im.rated_freq = 50. #64.7 # for vvvf
        ctrl = Controller2(im, timer)
        ualbe = array([0,0], dtype='float64')
        ualbe_prev = array([0,0], dtype='float64')


        # print im.k1, im.k2, im.k3, im.sigma, im.Lm/im.sigma/im.Ls/im.Lr, 1/im.sigma/im.Ls*(im.rs+im.Lm**2/im.Lr/im.Tr)
        # quit()

        ## Scope plotting arguments.
        # for control
        if controlStrategy_which.find('ifoc')!=-1:
            key_list = ['ualbe','ialbe','psi_r','mech','rpm']
            scp1 = Scope(key_list, 'Basic Electrical & mechanical Quantities')
            key_list = ['umt','im','it','irotor','fo_err']
            scp2 = Scope(key_list, 'IFOC Control Relevent')
            # key_list = ['pos_of_rotor_flux','omega_sl']
            # scp3 = Scope(key_list, 'pos_of_rotor_flux')
        elif controlStrategy_which == 'focVM_woDecouple':
            key_list = ['ualbe','ialbe','psi_r','mech','rpm']
            scp1 = Scope(key_list, 'Basic Electrical & mechanical Quantities')
            key_list = ['umt','im','it','irotor']
            scp2 = Scope(key_list, 'FOC Control Relevent')
        # for mras
        if mras_which == 'separateModel':
            key_list = ['est','est_err','rr','omg']
            scp_mras2 = Scope(key_list, r'MRAS for $\tau_r$ & $\omega_r$')
            key_list = ['est_err','rs','mag','phase']
            scp_mrasRs = Scope(key_list, r'MRAS for $r_s$')

        scp_tdda = Scope(['A','B','C'], 
                        'TEMP')

        # for prediction use
        tmp_st = array([0,0,0,0],dtype='float64')
        tmp_omega_r = 0.
        # tmp_sta_psi = array([0,0],dtype='float64')
        # tmp_CMFlux = array([0,0],dtype='float64')
        # tmp_ole_psi = array([0,0],dtype='float64')

        ## Evaluation for different control stradegy.
        evalCtrl = 0.
        evalCtrl_wiPenaltyOnLargeInput = 0.

    global rampOrflat
    rampOrflat = True
    turnOffRsAdapt = True
    DoF = 3.2
    # DoF = 0.0
    DoR = 0.15
    ampl = 0.05
    slope = ampl/DoR
    dutyRsAdaptOff = 0.5
    lstInstant, lstAmpl = zeros(6), zeros(6)
    lstInstant[0] = DoR; #array([0.15, 1.35, 1.5, 1.65, 2.85, 3])
    lstInstant[1] = lstInstant[0] + DoF; # 1.35
    lstInstant[2] = lstInstant[1] + DoR; # 1.50
    lstInstant[3] = lstInstant[2] + DoR; # 1.65
    lstInstant[4] = lstInstant[3] + DoF; # 2.85
    lstInstant[5] = lstInstant[4] + DoR; # 3.00
    lstAmpl[0] =  -1.0*ampl; #array([   1,    1,   0,   -1,   -1, 0])*ampl
    lstAmpl[1] =  -1.0*ampl;
    lstAmpl[2] =  0.0*ampl;
    lstAmpl[3] =  1.0*ampl;
    lstAmpl[4] =  1.0*ampl;
    lstAmpl[5] =  0.0*ampl;


    # ＭＡＩＮ　ＬＯＯＰ
    # the conclusion is that:
    # 1. i can be either 0 or 1, for omega_e*0 == omega_g*0. (omega_e is the excitation by vvvf)
    # 2. omega_g should update after omega_e, since omega_g == omega_r which is dependent on omega_e.
    for i in range(1, int(t_max / Ts)): # i should start from 1 for omega_e*t matching omega_g*t in vvvf
        #####
        ## Time evolves only when the states of machine are updated.
        timer.ticktock()
        ## reference to timer
        t = timer.t

        #####
        ## Exogenous input
        # rpm command & laod torque
        # rpm_cmd, Tload = cmd(t, delay_time, rpm=[(0,500),(1,500),(6.,10)], Tload=[(0,0),(0.01,20.)])
        # rpm_cmd, Tload = cmd(t, delay_time, rpm=[(0,0),(1,900),(6.,10)], Tload=[(0,0),(0.01,20.)])
        # rpm_cmd, Tload = cmd(t, delay_time, rpm=[(0,0),(0.3,500),(6.,900)], Tload=[(0,20.),(3,20.)])
        # rpm_cmd, Tload = cmd(t, delay_time, rpm=[(0,0),(0.3,500),(6.,500)], Tload=[(0,20.),(3,20.)])

        # rpm_cmd, Tload = cmd(t, delay_time, rpm=[(0,0),(0.3,20),(6.,20)], Tload=[(0,20.),(9,0.0),(15,20.)])
        # rpm_cmd, Tload = cmd(t, delay_time, rpm=[(0,0),(0.3,20),(6.,20)], Tload=[(0,0.),(9,0.0),(15,0.)])
        
        rpm_cmd, Tload = cmd(t, delay_time, 
                            rpm=[(0,0),(0.3,20),(1.,300),(2.,600),(3.,900),(4.,700),(5.,500),(6.,200),(7.,20)],
                            Tload=[(0,20.),(3,20.)])
        # Tload = 0.

        # ualbe = array([100., 0.])

        #####
        ## Machine of simulation
        # im.update2(ualbe, Tload)
        im.update_accurately(ualbe, Tload)
        ualbe_prev = copy.deepcopy(ualbe)
        # record ualbe cannot be located here

        #####
        ## MRAS
        if mras_which == 'separateModel':
            ctrl.mras.update(ualbe, tmp_st[0:2])
            # ctrl.mras.update(ualbe, im.st[0:2])

        #####
        ## Controller's operating period is Tctrl.

        # # ISE is applied on tracking error.
        # evalCtrl += Ts*(rpm_cmd - im.rpm)**2
        # # and penalty on large input should be added 
        # evalCtrl_wiPenaltyOnLargeInput += Ts*(ualbe[0]**2 + ualbe[1]**2)

        # cntSim += 1
        ## Separate the controllor and simulation of acm with different simulation frequency (Do this for higher accuracy to real machine)
        # if int((timer.t-Ts)/Tctrl) != int(timer.t/Tctrl) or timer.t == 0.:
            # cntCtrl += 1
        if controlStrategy_which == 'ifoc_woDecouple':
            # IFOC
            if sensorless == False:
                ualbe = ctrl.ifoc(rpm_cmd, im.st[0:2], im.omega_r, im.Tr)
            else:
                ualbe = ctrl.ifoc(rpm_cmd, im.st[0:2], ctrl.mras.thet[2], 1./ctrl.mras.thet[1])
                # ualbe = ctrl.ifoc(rpm_cmd, im.st[0:2], ctrl.mras.thet[2], im.Tr)
        elif controlStrategy_which == 'focVM_woDecouple':
            # FOC
            if sensorless == False:
                ualbe = ctrl.foc_VM(rpm_cmd, im.st[0:2], im.omega_r, im.st[2], im.st[3])
            else:
                ualbe = ctrl.foc_VM(rpm_cmd, im.st[0:2], ctrl.mras.thet[2], ctrl.mras.est_VM[0], ctrl.mras.est_VM[1])
        
        ## this should be put after the all estimation when new st should be recorded for next prediction
        tmp_st = copy.deepcopy(im.st) # Copy the value! Don't use the reference of an array (im.st is array) like tmp_st = im.st
        tmp_omega_r = copy.deepcopy(im.omega_r)
        # tmp_CMFlux = ctrl.mras.est
        # if mras_which == 'separateModel':
            # tmp_sta_psi = copy.deepcopy(ctrl.sta.itg_z)
            # tmp_ole_psi = copy.deepcopy(ctrl.mras.psi_ole)

        # if True:
        #     ualbe += ctrl.certainFreqExcite(100, 100) + ctrl.certainFreqExcite(100, 110)

        #####
        ## Scope
        ## Control
        # record for foc or ifoc with rotor flux estimating
        if controlStrategy_which.find('ifoc') != -1:

            # ualbe to record should be the evaluated output of controller with respect to measured currents
            # 逻辑很简单，其实就是新的电压激励到旧的电流状态，更新得到新的电流状态。
            scp1.recordTime(t)
            scp1.record('ualbe', array([ualbe[0], ualbe[1]]))
            scp1.record('ialbe', array([im.st[0], im.st[1]]))
            scp1.record('psi_r', array([ctrl.magPsi_SteadyStateCmd+0.2, sqrt(im.st[2]**2+im.st[3]**2), ctrl.magPsi_PICmd]))
            scp1.record('mech', array([im.theta_r, im.Tem,Tload]))
            scp1.record('rpm', array([im.rpm, rpm_cmd, rpm_cmd-im.rpm]))
            scp2.recordTime(t)
            scp2.record('umt', array([ctrl.umtn[0], ctrl.umtn[1]]))
            scp2.record('im', array([ctrl.imt_cmd[0], ctrl.imt[0]]))
            scp2.record('it', array([ctrl.imt_cmd[1], ctrl.imt[1]]))
            scp2.record('irotor', array([im.ir_MT[0], im.ir_MT[1]]))
            scp2.record('fo_err', array([arctan2(im.st[3],im.st[2]) - ctrl.theta_M, 0.]))
            # scp3.recordTime(t)
            # scp3.record('pos_of_rotor_flux', array([ctrl.theta_M, theta_psi_r, ctrl.theta_psi_r_est]))
            # scp3.record('omega_sl', array([ctrl.omega_sl, ctrl.omega_syn]))

            # omega_sl_cmd & omega_sl shall be observed.
            # omega_sl_cmd & omega_sl shall be observed.
            # omega_sl_cmd & omega_sl shall be observed.
        elif controlStrategy_which == 'focVM_woDecouple':
            # record for foc or ifoc with rotor flux estimating
            # key, data
            scp1.recordTime(t)
            scp1.record('ualbe', array([ualbe[0], ualbe[1]]))
            scp1.record('ialbe', array([im.st[0], im.st[1]]))
            scp1.record('psi_r', array([ctrl.magPsi_SteadyStateCmd, sqrt(im.st[2]**2+im.st[3]**2), ctrl.magPsi_PICmd]))
            scp1.record('mech', array([im.theta_r, im.Tem, Tload]))
            scp1.record('rpm', array([im.rpm, rpm_cmd, rpm_cmd-im.rpm]))
            scp2.recordTime(t)
            scp2.record('umt', array([ctrl.umtn[0], ctrl.umtn[1]]))
            scp2.record('im', array([ctrl.imt_cmd[0], ctrl.imt[0]]))
            scp2.record('it', array([ctrl.imt_cmd[1], ctrl.imt[1]]))
            scp2.record('irotor', array([im.ir_MT[0], im.ir_MT[1]]))
            
        ## Observer
        if mras_which == 'separateModel':
            # MRAS CM
            scp_mras2.recordTime(t)
            # scp_mras2.record('est', array([im.st[2],im.st[3],ctrl.mras.est_CM[0],ctrl.mras.est_CM[1],ctrl.mras.mag]))
            scp_mras2.record('est', array([ctrl.mras.est_VM[0],ctrl.mras.est_VM[1],ctrl.mras.est_CM[0],ctrl.mras.est_CM[1],ctrl.mras.mag]))
            scp_mras2.record('est_err', array([ctrl.mras.est_err[0],ctrl.mras.est_err[1]]))
            scp_mras2.record('rr', array([ctrl.mras.thet_real[1],ctrl.mras.thet[1]]))
            scp_mras2.record('omg', array([im.omega_r, ctrl.mras.thet[2]]))
            # scp_mras2.record('omg', array([im.omega_r/50., im.omega_r - ctrl.mras.thet[1]]))

            # MRAS VM
            PSIA = im.st[2]
            PSIB = im.st[3]
            if isclose(sqrt(PSIA**2+PSIB**2), 0.):
                PHASE = 0.
            else:
                PHASE = PSIA/sqrt(PSIA**2+PSIB**2)
            phase_err = PHASE-ctrl.mras.phase_psi_ole
            if isclose(phase_err, 1.):
                phase_err = 0.
            scp_mrasRs.recordTime(t)
            # scp_mrasRs.record('est_err', array([ctrl.mras.est_CM[0],ctrl.mras.est_VM[0]]))
            scp_mrasRs.record('est_err', array([im.st[2]-ctrl.mras.est_CM[0],im.st[2]-ctrl.mras.est_VM[0]]))
            scp_mrasRs.record('rs', array([im.rs, ctrl.mras.thet[0]]))
            scp_mrasRs.record('mag', array([sqrt(PSIA**2+PSIB**2),ctrl.mras.mag_psi_ole]))
            scp_mrasRs.record('phase', array([PHASE,ctrl.mras.phase_psi_ole]))            
    # end of ＭＡＩＮ　ＬＯＯＰ


    print('Begin to plot')

    if mras_which == 'separateModel':
        print('Rs', ctrl.mras.thet[0], im.rs)
        print('Tr', 1./ctrl.mras.thet[1], im.Tr)
        print('omega_r', ctrl.mras.thet[2], im.omega_r)

    print('ISE of tracking for %s is:'%(controlStrategy_which), evalCtrl)
    print('ISE of penalty for %s is:'%(controlStrategy_which), evalCtrl_wiPenaltyOnLargeInput, 'out of %g sec'%(t_max))
    
    if True:
        if controlStrategy_which.find('ifoc')!=-1:
            ## for foc or ifoc
            label_list = [r'$u_\alpha$',r'$u_\beta$',
                          r'$i_\alpha$',r'$i_\beta$',
                          r'$|\psi_r|^{*}$',r'$|\psi_r|$',r'$PI|\psi_r|^{*}$',
                          r'$\theta_r$',r'$T_{em}$',r'$T_{load}$',
                          r'$rpm$',r'$rpm_{cmd}$',r'$err$']
            scp1.plot(label_list)
            label_list = [r'$u_m$',r'$u_t$',
                          r'$i_m^{cmd}$',r'$i_m$',
                          r'$i_t^{cmd}$',r'$i_t$',
                          r'$i_Mr$',r'$i_Tr$',
                          r'$real\theta_M$',r'$given\theta_M$']
            scp2.plot(label_list)
            # label_list = ['theta_M','theta_psi_r','theta_psi_r_est',
            #               'omega_sl_cmd','omega_syn_cmd']
            # scp3.plot(label_list)
        elif controlStrategy_which == 'focVM_woDecouple':
            ## for foc or ifoc
            label_list = [r'$u_\alpha$',r'$u_\beta$',
                          r'$i_\alpha$',r'$i_\beta$',
                          r'$|\psi_r|^{*}$',r'$|\psi_r|$',r'$PI|\psi_r|^{*}$',
                          r'$\theta_r$',r'$T_{em}$',r'$T_{load}$',
                          r'$rpm$',r'$rpm_{cmd}$',r'$err$']
            scp1.plot(label_list)
            label_list = [r'$u_m$',r'$u_t$',
                          r'$i_m^{cmd}$',r'$i_m$',
                          r'$i_t^{cmd}$',r'$i_t$',
                          r'$i_Mr$',r'$i_Tr$']
            scp2.plot(label_list)
            # label_list = ['theta_M','theta_psi_r','theta_psi_r_est',
            #               'omega_sl_cmd','omega_syn_cmd']
            # scp3.plot(label_list)
        
    if mras_which == 'separateModel':
        # for mras
        label_list = [r'$\psi_{\alpha r}^{V\!M}$',r'$\psi_{\beta r}^{V\!M}$',r'$\hat\psi_{\alpha r}^{C\!M}$',r'$\hat\psi_{\beta r}^{C\!M}$',r'$|\psi_r|$',
                      r'$\varepsilon_\alpha$',r'$\varepsilon_\beta$',
                      r'$\alpha$',r'$\hat\alpha$',
                      r'$\omega_r$',r'$\hat\omega_r$']
        scp_mras2.plot(label_list)
        label_list = [r'$im-CM$',r'$im-VM$',
                      r'$r_s$',r'$\hat r_s$',
                      r'$real mag$',r'$est mag$',
                      r'$real phase$',r'$est phase$']
        scp_mrasRs.plot(label_list)

    scp_tdda.plot(['extraModulus','deriv','rampOrflat','offRsAdapt'])

    ## for plot in paper TDDA
    where = 'LmUncertainty'
    scp1.saveAsHDF5()
    scp2.saveAsHDF5(where)
    scp_mras2.saveAsHDF5(where)
    scp_mrasRs.saveAsHDF5(where)
    scp_mras2.saveAsHDF5(where)


    # ## for GCC TDDA.c
    # f_name = 'C:/MinGW/TDDA_SIM/ui.dat'
    # def WriteToFiles(f_name):
    #     lual = scp1.getData('ualbe',0)
    #     lube = scp1.getData('ualbe',1)
    #     lial = scp1.getData('ialbe',0)
    #     libe = scp1.getData('ialbe',1)
    #     lrampOrflat = scp_tdda.getData('C',0)
    #     lturnOffRsAdapt = scp_tdda.getData('C',1)
    #     lomg = scp_mras2.getData('omg',0)
    #     s = ''
    #     for ual, ube, ial, ibe, rampOrflat, turnOffRsAdapt, omg in zip(lual,lube,lial,libe,lrampOrflat,lturnOffRsAdapt,lomg):
    #         s += '%f, %f, %f, %f, %d, %d, %f\n'%(ual, ube, ial, ibe, rampOrflat, turnOffRsAdapt, omg)
    #         # s += '%f, %f, %f, %f, %d, %d\n'%(ual, ube, ial, ibe, rampOrflat, turnOffRsAdapt)
    #     with open(f_name, mode='w') as f:
    #         f.write(s)
    # WriteToFiles(f_name)

    ## for MATLAB TDDA.m
    # fh5 = h5py.File('./h5/uiout.h5','w') # h5py won't create directory for you. you have to create yourself.
    # fh5.create_dataset('ualbe', data=scp1.getData('ualbe','all')) #<- print 2??
    # fh5.create_dataset('ialbe', data=scp1.getData('ialbe','all'))
    # fh5.create_dataset('lrampOrflat', data=scp_tdda.getData('C',0))
    # fh5.create_dataset('lturnOffRsAdapt', data=scp_tdda.getData('C',1))
    # fh5.create_dataset('omg', data=scp_mras2.getData('omg',0))
    # fh5.close()

    show()


# print globals()
# print locals()

