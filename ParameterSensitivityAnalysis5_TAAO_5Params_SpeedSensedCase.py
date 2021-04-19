# coding:u8
# we encourage individuality, expressionism and unity through diversity.

from im_FOC_Simulator import Frame2, ThreePhaseInductionMachine
from pylab import *
from scipy.interpolate import griddata
def plotContour(coordinate, potential, tit, range_x, range_y):
    plt.style.use('ggplot') 
    mpl.rcParams['legend.fontsize'] = 8
    figure(num=None, figsize=(8, 8), dpi=80) #facecolor='w', edgecolor='k'

    X_mesh, Y_mesh = meshgrid(range_x, range_y)

    A = griddata(coordinate, potential, (X_mesh,Y_mesh))

    # CS = contour(X_mesh, Y_mesh, A, 20)
    CS = plt.contour(X_mesh, Y_mesh, A, 3, linewidths=0.5, colors='k')
    CS = plt.contourf(X_mesh, Y_mesh, A, 3, cmap=plt.cm.summer, alpha=0.3)    
    # CS = plt.contour(X_mesh, Y_mesh, A*scale, 15, linewidths=0.5, colors='k')
    # CS = plt.contourf(X_mesh, Y_mesh, A*scale, 15, cmap=plt.cm.rainbow)
    # CB = colorbar(CS, shrink=0.8, extend='both', ticks=linspace(0,1,5,endpoint=True)) # 2 ways: http://stackoverflow.com/questions/21952100/setting-the-limits-on-a-colorbar-in-matplotlib 
    CB = colorbar(CS, shrink=0.8, extend='both') # 2 ways: http://stackoverflow.com/questions/21952100/setting-the-limits-on-a-colorbar-in-matplotlib 
    # plt.clabel(CS, inline=1, fontsize=10) # 在线上显示云图的深度数值，过密的线会导致色块缺失。

    # xlabel(r'$\omega_r (rad/s)$')
    # ylabel(r'$Load Torque (Nm)$')
    xlabel(r'$\omega_r (%\omega_{rN})$')
    ylabel(r'$Load Torque (%\T_{LN})$')
    title('Sensitivity of ' + tit)

    # xlim([range_x[0],range_x[-1]])
    # ylim([range_y[0],range_y[-1]])

def mag(vector):
    return sqrt(vector[0]**2 + vector[1]**2)

def square(vector):
    return (vector[0]**2 + vector[1]**2)

initials = zeros(7+4+3); initials[-3]=0.9; initials[-2]=7. # not good practice
frame_albe = Frame2(trans_type='constant_power', theta_g_ini=0.)
im = ThreePhaseInductionMachine( Ts=1e-4, \
                                 rs=3.04, rr=1.26, \
                                 Lm=0.47, Lls=0.0126, Llr=0.0126, \
                                 # Lm=0.4144, Lls=0.0126, Llr=0.0126, \
                                 fr=frame_albe, J=0.06, P=2*2, initials=initials)

psiMm = 3.0
psiMm = 1.5
psiMm = 1.2
# psiMm = 0.5

ppsiMm = 5 # ramp
ppsiMm = 0 # flat


# max_Tload = 40 
# max_omega = (3000/60.*(im.P*pi))
# range_Tload = arange(-max_Tload, max_Tload+40/200., 40/200.)
# range_omega = arange(-max_omega, max_omega+10, 10)

max_Tload = 4000.0/(1440.0/60*2*pi) # 26.5 Nm
max_omega = (1440/60.*(2*pi)) # 此处应取整
range_Tload = arange(-max_Tload, max_Tload+.2, .2)
range_omega = arange(-max_omega, max_omega+1, 1)
# range_omega = arange(-200/60*(im.P*pi), 200/60*(im.P*pi)+1,1)


Jadot2009 = False

S_rs, S_Lsigma, S_rreq, S_Leq, S_omegar = [], [], [], [], []
S_alpha = []
if Jadot2009:
    S_rreq, S_Leq = [], []
coordinate = []
for rind, Tload in enumerate(range_Tload):
    for cind, omega_cmd in enumerate(range_omega):

        iMs_for_alpha =  ppsiMm/im.rreq + psiMm/im.Leq 
        piMs = 0 # steady state of ramp struct., the derivative of iMs is zero, but ppsiMm is nonzero.
        piTs = 0
        iMs = iMs_for_alpha 

        iTs = Tload / (im.P*0.5*psiMm)
        omega_sl = im.rreq*iTs/psiMm
        omega_e = omega_cmd + omega_sl

        # coordinate.append((omega_cmd, Tload))
        # Normalized
        coordinate.append((omega_cmd/max_omega, Tload/max_Tload*2.0))

        if not Jadot2009:
            # 1. Sensivity Check by Jadot2009. Do note that the \Delta\theta vector in Jadot2009 is normalised already!
            # check out different norms: max or sum
            a1 = mag((             iMs, iTs))
            a2 = mag((piMs-iTs*omega_e, piTs+iMs*omega_e))
            a3 = mag((iMs-psiMm/im.Leq, iTs)) # is = imu - ireq
            a4 = mag((    psiMm/im.Leq, 0))
            # a5 = mag((              -0, psiMm))
            den = sqrt(sum([a1**2      *(im.rs)**2,
                            a2**2      *(im.Lsigma)**2,
                            a3**2      *(im.rreq)**2,
                            a4**2      *(im.rreq)**2]))
                            # a5**2      *(omega_cmd)**2]))
            
            S_rs.append(        a1 *(im.rs)       /den)
            S_Lsigma.append(    a2 *(im.Lsigma)   /den)
            S_rreq.append(      a3 *(im.rreq)     /den)
            S_Leq.append(       a4 *(im.rreq)     /den)
            # S_omegar.append(    a5 *abs(omega_cmd)/den)

        else:
            # Jadot2009
            # use max or sum to get deffirent results
            den = sqrt(sum([square((im.alpha, -omega_sl))   *(omega_e*omega_sl*psiMm/(im.alpha**2+omega_sl**2))**2,
                            square((omega_sl, im.alpha))    *(-omega_e*im.alpha*psiMm/(im.alpha**2+omega_sl**2))**2,
                            square((iMs, iTs))              *(im.rs)**2,
                            square((iTs, iMs))              *(omega_e*im.Lsigma)**2]))
            # den = sqrt(max([square((im.alpha, -omega_sl))     *(omega_e*omega_sl*psiMm/(im.alpha**2+omega_sl**2))**2,
            #                 square((omega_sl, im.alpha))      *(-omega_e*im.alpha*psiMm/(im.alpha**2+omega_sl**2))**2,
            #                 square((iMs, iTs))                *(im.rs)**2,
            #                 square((iTs, iMs))                *(omega_e*im.Lsigma)**2]))

            S_rreq.append(   abs(mag((im.alpha, -omega_sl))  *(omega_e*omega_sl*psiMm/(im.alpha**2+omega_sl**2))     /den))
            S_Leq.append(    abs(mag((omega_sl, im.alpha))   *(-omega_e*im.alpha*psiMm/(im.alpha**2+omega_sl**2))    /den))
            S_rs.append(     abs(mag((iMs, iTs))             *(im.rs)             /den))
            S_Lsigma.append( abs(mag((iTs, iMs))             *(omega_e*im.Lsigma) /den))


def limitToLessThanOne(S):
    for ind,el in enumerate(S):
        if el>1.:
            S[ind] = 1.


if Jadot2009:
    range_x = array(range_omega)/max_omega * 1
    range_y = array(range_Tload)/max_Tload * 2

    print('''只从磁链幅值中提取rs，在不带载的情况下无法辨识。注意到如果要辨识转速，必须牺牲磁链的相位供之使用（用Hory2015_decoupleAdaption说明），
    既然如此，我们只需要分析各个参数在磁链幅值中的敏感性即可，反正磁链的相位要用于转速辨识。''')
    plotContour(coordinate, S_rreq, 'rreq', range_x, range_y)
    plotContour(coordinate, S_Leq, 'Leq', range_x, range_y)
    plotContour(coordinate, S_rs, 'rs', range_x, range_y)
    plotContour(coordinate, S_Lsigma, 'Lsigma', range_x, range_y)
    # plotContour(coordinate, S_omegar, 'omega_r', range_x, range_y)
    print('能不能在接近原点的地方剖分得密一些，google一下，不同精度的mesh？')
else:

    limitToLessThanOne(S_rs)
    limitToLessThanOne(S_Lsigma)
    limitToLessThanOne(S_rreq)
    limitToLessThanOne(S_Leq)
    # limitToLessThanOne(S_omegar)

    range_x = array(range_omega)/max_omega * 1
    range_y = array(range_Tload)/max_Tload * 2

    print('''只从磁链幅值中提取rs，在不带载的情况下无法辨识。注意到如果要辨识转速，必须牺牲磁链的相位供之使用（用Hory2015_decoupleAdaption说明），
    既然如此，我们只需要分析各个参数在磁链幅值中的敏感性即可，反正磁链的相位要用于转速辨识。''')
    # plotContour(coordinate, S_rs, 'rs', range_x, range_y)
    # plotContour(coordinate, S_Lsigma, 'Lsigma', range_x, range_y)
    # plotContour(coordinate, S_alpha, 'alpha', range_x, range_y)
    # plotContour(coordinate, S_omegar, 'omega_r', range_x, range_y)

    print('能不能在接近原点的地方剖分得密一些，google一下，不同精度的mesh？')

    plt.style.use('ggplot') 
    # plt.style.use('grayscale') # print plt.style.available # get [u'dark_background', u'bmh', u'grayscale', u'ggplot', u'fivethirtyeight']# 
    mpl.rcParams['legend.fontsize'] = 14
    # figure(num=None, figsize=(8, 8), dpi=75) #facecolor='w', edgecolor='k'
    fig, axes = plt.subplots(ncols=2, nrows=3, dpi=150, sharex=True, sharey=True); # modifying dpi would change the size of fonts
    # ax_omg = subplot(3, 1, 3) # the 3rd of 3-row-1-col subplot
    fig.subplots_adjust(right=0.95, bottom=0.1, top=0.95, left=0.05, hspace=0.2, wspace=0.01)
    ax_list = axes.ravel()
    X_mesh, Y_mesh = meshgrid(range_x, range_y)

    # print ax_list
    # quit()
    def ct(A,ax):
        # ax = fig.add_subplot(2,2,i)
        CS = ax.contour(X_mesh*100, Y_mesh*100, A, 3, linewidths=0.5, colors='k', vmin=0, vmax=1)
        CS = ax.contourf(X_mesh*100, Y_mesh*100, A, 3, cmap=plt.cm.gray, alpha=1, vmin=0, vmax=1)
        CB = colorbar(CS, shrink=1.0, extend='both', ax=ax) # 2 ways: http://stackoverflow.com/questions/21952100/setting-the-limits-on-a-colorbar-in-matplotlib 
    
        ax.set_xlim(-1*100,1*100)
        ax.set_xticks(arange(-100, 101, 10), minor=True) # 标出点，但是没有数字
        ax.set_ylim(-2*100,2*100)
        ax.set_yticks(arange(-200, 201, 100))                                                       
        ax.set_yticks(arange(-200, 201, 50), minor=True)         
        return CS

    for S, tit, ax in zip([S_rs, S_Lsigma, S_rreq, S_Leq], [r'$\tilde r_s$', r'$D\!\!\;\tilde L_\sigma$', r'$\tilde r_{req}$', r'$\tilde L_{\mu}$'], ax_list):
        if S!=None:
            A = griddata(coordinate, S, (X_mesh,Y_mesh))
            if tit == r'$\tilde \omega$':
                ax = ax_omg
                CS = ct(A, ax)
            else:
                ct(A, ax)    

            ax.set_title('Normalized Sensitivity of '+tit, fontsize=10, **{'fontname':'Times New Roman'})
            # ax.locator_params(nbins=4)

    # fig.subplots_adjust(right=0.8)
    # cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    # fig.colorbar(CS, cax=cbar_ax)


    ax.set_ylabel(r'$T_L [\%T_{LN}]$')
    ax.set_xlabel(r'$\omega_r [\%\omega_{rN}]$ with $\omega_{rN}=1440$ rpm')
    for i in [0,2]:
        ax_list[i].set_ylabel(r'$T_L [\%T_{LN}]$')
    


    fig.tight_layout() # in order to avoid overlapping labels.
    #savefig(r'C:\Dr.H\[03]Ref\p\IEMDC\TotallyAdaptOb_TEX\pic\fig01', dpi=300)
    
    # if isclose(ppsiMm, 0.):
    #     savefig(r'./SensitivityAnalysisLeqKnown', dpi=300)
    # else:
    #     savefig(r'./SensitivityAnalysisLeqKnown_Ramp', dpi=300)    

show()


