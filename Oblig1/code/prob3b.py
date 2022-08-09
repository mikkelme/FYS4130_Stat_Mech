import numpy as np
import matplotlib.pyplot as plt
from plot_set import *


alpha = 1
gamma = 10
def cal_F(N_arr, V, T):
    """ return F with repsect to configuration N_rr,
        temperature T and dimensionless volume V """
    sum = 0
    for N_ in N_arr[N_arr > 0]:
        sum += N_*np.log(alpha*N_/V)
    return T*(sum + gamma*(N_arr[0]*N_arr[1] + N_arr[1]*N_arr[2] + N_arr[2]*N_arr[0])/V)

def cal_P(N_arr, V, T):
    """ Calculate pressure with respect configuration N_rr
        dimensionless volume V and temperature T """
    Nx, Ny, Nz = N_arr
    P = T*(Nx/V + Ny/V + Nz/V + gamma*(Nx*Ny + Ny*Nz + Nz*Nx)/V**2)
    return P


def minimum_F(N, V, T):
    """ Go through all possible configurations and
        return minimum F with respect to number of particles N,
        dimensionless volume and temperature T """
    minF = 1e10
    N = int(N)
    state = np.zeros(3)
    trial_state = np.zeros(3)
    for Nx in range(N+1):
        for Ny in range(N-Nx+1):
            Nz = N - Nx - Ny
            trial_state[:] = (Nx, Ny, Nz)

            F = cal_F(trial_state, V, T)

            if F < minF:
                minF = F
                state[:] = trial_state[:]
    state.sort()
    return minF, state


def equil_states_N(N_start, N_end, V, T):
    """ Gather equlibrium F, P and configurations for
        increasing N for dimensionsless volume V and temperature T """
    N = np.linspace(N_start, N_end, N_end - N_start + 1)
    F = np.zeros(len(N))
    P = np.zeros(len(N))
    states = np.zeros((len(N), 3))
    for i in range(len(N)):
        F[i], states[i] = minimum_F(N[i], V, T)
        P[i] = cal_P(states[i], V, T)
        print(f'\rN = {int(N[i])}/{int(N[-1])}, state = {states[i]}', end='')
    print()

    ddFddN = double_derivative(N, F)
    dd_neg = np.argwhere(ddFddN < 0).ravel()
    phase_trans = [dd_neg[0]]
    for i in range(1, len(dd_neg)):
        if dd_neg[i] - dd_neg[i-1] > 2:
            break
        phase_trans.append(dd_neg[i])
    phase_trans = np.array(phase_trans)

    save = True
    n = N/V
    print(f'Estimated phase trans, n: {n[phase_trans][0]}, {n[phase_trans][-1]}')

    axis_label_size = 16
    plt.figure(num=0, dpi=80, facecolor='w', edgecolor='k')
    plt.plot(n, F)
    plot_area(n[phase_trans], plt.gca())
    plt.xlabel(r"$\tilde{n}$", fontsize = axis_label_size)
    plt.ylabel(r"$F_{eq}$", fontsize = axis_label_size)
    plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
    if save:
        plt.savefig("../article/figures/Feq.pdf", bbox_inches="tight")

    plt.figure(num=1, dpi=80, facecolor='w', edgecolor='k')
    plt.plot(n, P)
    plot_area(n[phase_trans], plt.gca())

    plt.xlabel(r"$\tilde{n}$", fontsize = axis_label_size)
    plt.ylabel(r"$P_{eq}$", fontsize = axis_label_size)
    plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
    if save:
        plt.savefig("../article/figures/Peq.pdf", bbox_inches="tight")

    plt.figure(num=2, dpi=80, facecolor='w', edgecolor='k')
    plt.plot(n, ddFddN)
    plot_area(n[phase_trans], plt.gca())
    plt.xlabel(r"$\tilde{n}$", fontsize = axis_label_size)
    plt.ylabel(r"$\frac{\partial^2 F_{eq}}{\partial N^2}$", fontsize = 20)
    plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
    if save:
        plt.savefig("../article/figures/ddFddN.pdf", bbox_inches="tight")

    plt.figure(num=3, dpi=80, facecolor='w', edgecolor='k')
    plt.plot(n, states[:,0], color = color_cycle(0), label = "# $N_x$")
    plt.plot(n, states[:,1], color = color_cycle(2), label = "# $N_y$")
    plt.plot(n, states[:,2], color = color_cycle(3), label = "# $N_z$")
    plot_area(n[phase_trans], plt.gca())
    plt.xlabel(r"$\tilde{n}$", fontsize= axis_label_size)
    plt.ylabel(r"Rod configuration", fontsize=14)
    plt.legend(fontsize = 13)

    plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
    if save:
        plt.savefig("../article/figures/Rod_conf.pdf", bbox_inches="tight")


    plt.show()


def plot_area(x, ax, color = color_cycle(1)):
    """ Shade x domain area with red color """
    ax = plt.gca()
    ylim = ax.get_ylim()
    plt.fill_between(x, ylim[0], ylim[1], alpha = 0.1, color = color)
    ax.set_ylim(ylim)

def double_derivative(x, y):
    """ Second order central finite difference """
    h = 1
    dd = np.zeros(len(x))
    dd[:] = np.nan
    for i in range(h, len(y)-h):
        dd[i] = (y[i+h] - 2*y[i] + y[i-h])/h**2

    return dd


if __name__ == "__main__":
    V = 400
    T = 1
    N_start = 3
    N_end = V
    equil_states_N(N_start, N_end, V, T)
