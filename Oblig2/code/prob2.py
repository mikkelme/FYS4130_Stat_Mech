from analyze import *
from scipy import interpolate


def prob2b():
    filename = '1D_chain_2b.txt'
    system = analyze(filename)

    plt.figure(num=0, dpi=80, facecolor='w', edgecolor='k')

    r_num = np.arange(system.L + 1)
    r_ana = np.linspace(r_num[0], r_num[-1], int(1e4))

    alpha_std = 0.50
    alpha_ana = 0.75
    linewidth_ana = 1.5
    error_dis = 10

    for i in range(len(system.C_r)):
        C_extended = np.insert(system.C_r[i], system.C_r.shape[1] , system.C_r[i,0])
        C_extended_std = np.insert(system.std_C_r[i], system.std_C_r.shape[1] , system.std_C_r[i,0])


        plt.plot(r_num, C_extended, 'o',color = color_cycle(i), label = f'T/J = {system.T[i]}')


        if i < len(system.C_r) -1:
            plt.fill_between(r_num, C_extended - error_dis*C_extended_std, C_extended + error_dis*C_extended_std, alpha = 0.25, color = color_cycle(i))
            plt.plot(r_ana, C(r_ana, system.T[i], system.L), '--', linewidth = linewidth_ana, color = color_cycle(0), alpha = alpha_ana)

        else: 
            plt.fill_between(r_num, C_extended - error_dis*C_extended_std, C_extended + error_dis*C_extended_std, alpha = 0.25, color = color_cycle(i),  label =  "$\pm$ " + str(error_dis) + "s")
            analytic, = plt.plot(r_ana, C(r_ana, system.T[i], system.L), '--', linewidth = linewidth_ana, color = color_cycle(1), alpha = alpha_ana, label = "Analytic")


    plt.xlabel(r"$r$", fontsize=14)
    plt.ylabel(r"$C(r)$", fontsize=14)
    plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
    leg = plt.legend(fontsize = 13)
    leg.legendHandles[2].set_color('grey')
    leg.legendHandles[3].set_color('grey')

    plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
    plt.savefig("../article/figures/fig2b.pdf", bbox_inches="tight")
    plt.show()



def prob2c():
    system = analyze(filename = '2D_chain_2c.txt')
    error_dis = 1
    alpha_std = 0.25
  
    plt.figure(num=0, dpi=80, facecolor='w', edgecolor='k')
    plt.plot(system.T, system.m, 'o-', color = color_cycle(3), label = "Data points")
    plt.fill_between(system.T, system.m - error_dis*system.std_m, system.m + error_dis*system.std_m, alpha = alpha_std, color = color_cycle(0), label =  "$\pm$ " + str(error_dis) + "s" )
    
    plt.xlabel(r"$T/J$", fontsize=14)
    plt.ylabel(r"$\langle m \rangle$", fontsize=14)
    plt.legend(fontsize = 13)
    plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
    plt.savefig("../article/figures/fig2c.pdf", bbox_inches="tight")
    plt.show()


def prob2d():
    sys = analyze(filename = '2D_chain_2c.txt')
    sys_details = analyze(filename = '2D_chain_2d.txt')
    
    alpha_std = 0.50
    error_dis = 100


    T1 = sys_details.T[0]
    T2 = sys_details.T[-1]
    T3 = 4
    sec1 = np.argwhere(sys.T <= T1).ravel()
    sec2 = np.argwhere(np.logical_and(T2 <= sys.T, sys.T <= T3)).ravel()
    T_c_theory = 1 / np.log(1 + np.sqrt(3))
  
    plt.figure(num=0, dpi=80, facecolor='w', edgecolor='k')

    # sec1
    plt.plot(sys.T[sec1], sys.m2[sec1], 'o', color = color_cycle(4))
    plt.fill_between(sys.T[sec1], sys.m2[sec1] - error_dis*sys.std_m2[sec1], sys.m2[sec1] + error_dis*sys.std_m2[sec1], alpha = alpha_std, color = color_cycle(0) )
    
    # sec2
    plt.plot(sys.T[sec2], sys.m2[sec2], 'o-', color = color_cycle(4), label = r"$\Delta T/J = 0.1$")
    plt.fill_between(sys.T[sec2], sys.m2[sec2] - error_dis*sys.std_m2[sec2], sys.m2[sec2] + error_dis*sys.std_m2[sec2], alpha = alpha_std, color = color_cycle(0) )

    # details     
    plt.plot(sys_details.T, sys_details.m2, 'o-', color = color_cycle(1),  label = r"$\Delta T/J = 0.01$")
    plt.fill_between(sys_details.T, sys_details.m2 - error_dis*sys_details.std_m2, sys_details.m2 + error_dis*sys_details.std_m2, alpha = alpha_std, color = color_cycle(0), label =  "$\pm$ " + str(error_dis) + "s" )

    plt.vlines(T_c_theory, 0, 1, linestyle = '--', linewidth = 1, color = 'k', label = '$T_c/J$ (exact)')

    plt.xlabel(r"$T/J$", fontsize=14)
    plt.ylabel(r"$\langle |m|^2 \rangle$", fontsize=14)
    plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
    plt.legend(fontsize = 13)

    plt.savefig("../article/figures/fig2d.pdf", bbox_inches="tight")
    plt.show()


def prob2f():
    sysL8 = analyze('2D_chain_2f_L8.txt')
    sysL16 = analyze('2D_chain_2f_L16.txt')
    sysL32 = analyze('2D_chain_2f_L32.txt')
    system = [sysL8, sysL16, sysL32]

    error_dis = 1
    alpha_std = 0.3

    plt.figure(num=0, dpi=80, facecolor='w', edgecolor='k')
    for i, sys in enumerate(system):
       
        plt.plot(sys.T, sys.gamma, '-o', markersize = 4, label = f"L = {sys.L}", color = color_cycle(i))

        if i < len(system) - 1:
            plt.fill_between(sys.T, sys.gamma - error_dis*sys.std_gamma, sys.gamma + error_dis*sys.std_gamma, alpha = alpha_std, color = color_cycle(i) )
        else: 
            plt.fill_between(sys.T, sys.gamma - error_dis*sys.std_gamma, sys.gamma + error_dis*sys.std_gamma, alpha = alpha_std, color = color_cycle(i), label =  "$\pm$ " + str(error_dis) + "s" )


    # Find intersections
    f8 = interpolate.interp1d(sysL8.T, sysL8.gamma)
    f16 = interpolate.interp1d(sysL16.T, sysL16.gamma)
    f32 = interpolate.interp1d(sysL32.T, sysL32.gamma)
    systems = [sysL8, sysL16, sysL32]
    intsec = []
    intsec_err = []
    
    target_temp = 0.995
    temp_range = 0.005
    search_points = int(1e5)
    search_temp = np.linspace(target_temp - temp_range/2, target_temp + temp_range/2, search_points)
    search_temp = np.linspace(sysL8.T[0], sysL8.T[-1], search_points)

    for i in range(len(systems) - 1):
        f1_low = interpolate.interp1d(systems[i].T, systems[i].gamma - systems[i].std_gamma)
        f1 = interpolate.interp1d(systems[i].T, systems[i].gamma)
        f1_high = interpolate.interp1d(systems[i].T, systems[i].gamma + systems[i].std_gamma)

        f2_low = interpolate.interp1d(systems[i+1].T, systems[i+1].gamma - systems[i+1].std_gamma)
        f2 = interpolate.interp1d(systems[i+1].T, systems[i+1].gamma)
        f2_high = interpolate.interp1d(systems[i+1].T, systems[i+1].gamma + systems[i+1].std_gamma)

        lower_intsec = search_temp[np.argmin(np.abs(f1_low(search_temp) - f2_high(search_temp)))]
        intsec.append(search_temp[np.argmin(np.abs(f1(search_temp) - f2(search_temp)))])
        high_intsec = search_temp[np.argmin(np.abs(f1_high(search_temp) - f2_low(search_temp)))]    
        intsec_err.append(np.abs(lower_intsec - high_intsec)/2)


    print("#---Intersections---#")
    print(f'L8  >< L16, T_c = {intsec[0]:.{decimals(intsec_err[0])}f} +- {intsec_err[0]:1.0e}')
    print(f'L16 >< L32, T_c = {intsec[1]:.{decimals(intsec_err[0])}f} +- {intsec_err[1]:1.0e}')
    print("#---Averaged intersection---#")
    avg_intsec = np.mean(intsec)
    avg_intsec_err = np.mean(intsec_err)
    print(f'Avg intersection, T_c = {avg_intsec:.{decimals(avg_intsec_err)}f} +- {avg_intsec_err:1.0e}')


    plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
    plt.xlabel(r"$T/J$", fontsize=14)
    plt.ylabel(r"$\Gamma$", fontsize=14)
    leg = plt.legend(fontsize = 13)
    leg.legendHandles[3].set_color('grey')
    plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
    plt.savefig("../article/figures/fig2f.pdf", bbox_inches="tight")
    plt.show()




if __name__ == '__main__':
    prob2b()
    prob2c()
    prob2d()
    prob2f()



