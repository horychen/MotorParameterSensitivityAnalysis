# coding:u8
# we encourage individuality, expressionism and unity through diversity.

from pylab import *
from im_FOC_Simulator import Frame2, ThreePhaseInductionMachine

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
    CB = colorbar(CS, shrink=0.8, extend='both')
    # plt.clabel(CS, inline=1, fontsize=10) # 在线上显示云图的深度数值，过密的线会导致色块缺失。

    # xlabel(r'$\omega_r (rad/s)$')
    # ylabel(r'$Load Torque (Nm)$')
    xlabel(r'$\omega_r (%\omega_{rN})$')
    ylabel(r'$Load Torque (%\T_{LN})$')
    title('Sensitivity w.r.t. ' + tit)

    # xlim([range_x[0],range_x[-1]])
    # ylim([range_y[0],range_y[-1]])

def mag(vector):
    return sqrt(vector[0]**2 + vector[1]**2)

def square(vector):
    return (vector[0]**2 + vector[1]**2)

initials = zeros(7+4+3); initials[-3]=0.9; initials[-2]=7. # not good practice
frame_albe = Frame2(trans_type='constant_power', 
                    theta_g_ini=0.)
im = ThreePhaseInductionMachine( Ts=1e-4, \
                                 rs=1.1, rr=1.3, \
                                 Lm=0.15, Lls=0.0064, Llr=0.00125, \
                                 fr=frame_albe, J=0.06, P=3*2, initials=initials)

psiMm = 3.0
psiMm = 1.0
# psiMm = 0.5
ppsiMm = 0.0

range_Tload = arange(-200,200,1)
range_omega = arange(-1000/60*(im.P*pi), 1000/60*(im.P*pi)+5,5)
# range_omega = arange(-200/60*(im.P*pi), 200/60*(im.P*pi)+1,1)

S_rs, S_Lsigma, S_omegar, S_Leq, S_rreq = [], [], [], [], []
coordinate = []
for rind, Tload in enumerate(range_Tload):
    for cind, omega_cmd in enumerate(range_omega):

        iMs =  ppsiMm/im.rreq + psiMm/im.Leq
        iTs = Tload / (im.P*0.5*psiMm)
        omega_sl = im.rreq*iTs/psiMm
        omega_e = omega_cmd + omega_sl

        # coordinate.append((omega_cmd, Tload))
        # Normalized
        coordinate.append((omega_cmd*0.00318309886184, Tload*0.01))


        if False:
            den = sqrt(max([square((iMs, iTs))      *(im.rs)**2,
                            square((-iTs, iMs))     *(omega_e*im.Lsigma)**2,
                            square((0, psiMm))      *(omega_cmd)**2,
                            square((psiMm, 0))      *(im.alpha)**2, # old Leq
                                # square((psiMm, 0))      *(omega_sl)**2, # new Leq
                            square((-im.alpha*psiMm+iMs, iTs*im.rreq))]))

            S_rs.append(        mag((iMs, iTs))                 *(im.rs)                    /den)
            S_Lsigma.append(    mag((-iTs, iMs))                *abs((omega_e*im.Lsigma))   /den)
            S_omegar.append(    mag((0, psiMm))                 *abs(omega_cmd)             /den)
            S_Leq.append(       mag((psiMm, 0))                 *(im.alpha)                 /den) # old Leq
                # S_Leq.append(       mag((psiMm, 0))                 *abs(omega_sl)                 /den) # new Leq
            S_rreq.append(      mag((-im.alpha/im.rreq**2*psiMm*0+iMs, iTs/im.rreq))          /den)
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

# limitToLessThanOne(S_rs)
# limitToLessThanOne(S_Lsigma)
# limitToLessThanOne(S_omegar)
# limitToLessThanOne(S_Leq)
# limitToLessThanOne(S_rreq)


if False:
    plotContour(coordinate, S_rs, 'rs', range_omega, range_Tload)
    plotContour(coordinate, S_Lsigma, 'Lsigma', range_omega, range_Tload)
    plotContour(coordinate, S_Leq, 'Leq', range_omega, range_Tload)
    plotContour(coordinate, S_rreq, 'rreq', range_omega, range_Tload)
else:
    range_x = array(range_omega)*0.00318309886184
    range_y = array(range_Tload)*0.01

    print('''只从磁链幅值中提取rs，在不带载的情况下无法辨识。注意到如果要辨识转速，必须牺牲磁链的相位供之使用（用Hory2015_decoupleAdaption说明），
    既然如此，我们只需要分析各个参数在磁链幅值中的敏感性即可，反正磁链的相位要用于转速辨识。''')
    plotContour(coordinate, S_rreq, 'rreq', range_x, range_y)
    plotContour(coordinate, S_Leq, 'Leq', range_x, range_y)
    plotContour(coordinate, S_rs, 'rs', range_x, range_y)
    plotContour(coordinate, S_Lsigma, 'Lsigma', range_x, range_y)
    # plotContour(coordinate, S_omegar, 'omega_r', range_x, range_y)
    print('能不能在接近原点的地方剖分得密一些，google一下，不同精度的mesh？')
show()


