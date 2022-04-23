from analyze import *


def plot_Cr():
    T_list = [0.25, 0.5, 0.75, 10, 100] 
    L = 10
    r = np.linspace(0, L, int(1e4))

    plt.figure(num=0, dpi=80, facecolor='w', edgecolor='k')
    for T in T_list:   
        plt.plot(r, C(r, T, L), label = f'T = {T} J')
        plt.xlabel(r"$r$", fontsize=14)
        plt.ylabel(r"$C(r)$", fontsize=14)
        plt.legend(fontsize = 13)
    plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
    plt.savefig("../article/figures/Cr_L10.pdf", bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    plot_Cr()








