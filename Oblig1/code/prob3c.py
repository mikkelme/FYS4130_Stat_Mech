from prob3b import *




def equil_states_V(V_start, V_end, N, T):
    """ Gather equlibrium F, G, P and configurations for
        increasing V for number of particles N and temperature T """
    V = np.linspace(V_start, V_end, V_end - V_start + 1)
    F = np.zeros(len(V))
    G = np.zeros(len(V))
    P = np.zeros(len(V))
    states = np.zeros((len(V), 3))

    for i in range(len(V)):
        F[i], states[i] = minimum_F(N, V[i], T)
        P[i] = cal_P(states[i], V[i], T)
        G[i] = F[i] + P[i]*V[i]
        print(f'\r{i}/{len(V)}, V = {V[i]:.2f}, ngamma = {gamma*N/V[i]:.2f}, state = {states[i]}', end='')
    print()

    dPdV = single_derivative(V, P)
    ddFddV = -dPdV
    phase_trans = np.argwhere(ddFddV < 0).ravel()


    save = True
    axis_label_size = 16

    plt.figure(num=0, dpi=80, facecolor='w', edgecolor='k')
    plt.plot(P[phase_trans], V[phase_trans], linewidth = 6, alpha = 0.3, color = color_cycle(1))
    plt.plot(P, V, '-o', markersize = 5)
    plot_area(P[phase_trans], plt.gca(), color = color_cycle(4))
    plt.gca().axvline(0.474811, linestyle = '--', color = color_cycle(1), alpha = 0.5)

    plt.xlabel(r"$P_{eq}$", fontsize = axis_label_size)
    plt.ylabel(r"$\tilde{V}$", fontsize = axis_label_size)
    plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
    if save:
        plt.savefig("../article/figures/VeqP.pdf", bbox_inches="tight")

    plt.figure(num=1, dpi=80, facecolor='w', edgecolor='k')
    plt.plot(P[phase_trans], G[phase_trans], linewidth = 6, alpha = 0.3, color = color_cycle(1))
    plt.plot(P, G, '-o', markersize = 5)
    plot_area(P[phase_trans], plt.gca(), color = color_cycle(4))

    plt.gca().axvline(0.474811, linestyle = '--', color = color_cycle(1), alpha = 0.5)
    plt.xlabel(r"$P_{eq}$", fontsize = axis_label_size)
    plt.ylabel(r"$G_{eq}$", fontsize = axis_label_size)
    plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
    if save:
        plt.savefig("../article/figures/Geq.pdf", bbox_inches="tight")
    plt.xlim([0.34, 0.66])
    plt.ylim(-2.5, 5)
    if save:
        plt.savefig("../article/figures/Geq_zoom.pdf", bbox_inches="tight")


    plt.figure(num=2, dpi=80, facecolor='w', edgecolor='k')
    plt.plot(P[phase_trans], -dPdV[phase_trans], linewidth = 6, alpha = 0.3, color = color_cycle(1))
    plt.plot(P, ddFddV, '-o', markersize = 5)
    plt.gca().axvline(0.474811, linestyle = '--', color = color_cycle(1), alpha = 0.5)
    plot_area(P[phase_trans], plt.gca(), color = color_cycle(4))

    plt.xlabel(r"$P_{eq}$", fontsize = axis_label_size)
    plt.ylabel(r"$\frac{\partial^2 F_{eq}}{\partial V^2} = - \frac{\partial P_{eq}}{\partial V}$", fontsize = 20)
    plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
    if save:
        plt.savefig("../article/figures/ddFddV.pdf", bbox_inches="tight")

    plt.show()


def single_derivative(x, y):
    """ First order central finite difference """
    h = 1
    d = np.zeros(len(x))
    d[:] = np.nan
    for i in range(h, len(y)-h):
        d[i] = (y[i+h] -  y[i-h])/(2*h)

    return d


if __name__ == "__main__":
    N = 5
    V_start = N
    V_end = 100
    T = 1
    equil_states_V(V_start, V_end, N, T)
