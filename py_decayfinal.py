"""This version is used to plot final decay-profiles 
      in the current working directory.
"""
import os
import pathlib
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def plot_profile_final(snap, rho0, alpha0):
    filebase = os.path.splitext(snap)[0]
    length = len(filebase)
    frame = int(filebase[5:length - 2])
    df = np.array(pd.read_csv(snap, delimiter=",", header=None))

    fig, axs = plt.subplots(2, 1)
    axs[0].plot(df[:, 0], df[:, 1], "b-", lw=2)
    axs[0].set_xlim(0, np.max(df[:, 0]))

    axs[1].plot(df[:, 0], df[:, 2], "r-", lw=2, label=r"$f_1^R$")
    # axs[1].plot(df[:, 0], df[:, 3], "g-", lw=2, label=r"$f_1^I$")
    plt.legend(loc=0, fontsize=14)
    axs[1].set_xlim(0, np.max(df[:, 0]))
    axs[1].set_xlabel(r"$X$", fontsize=18)
    # plt.suptitle(r"$\alpha_0={0:.2f}$".format(alpha0), fontsize=18)

    fig.savefig(
        r"../../pk2={0:.2f}_alpha0={1:.2f}.png".format(rho0, alpha0),
        dpi=300,
        bbox_inches='tight',
        pad_inches=0)
    plt.close()


def profile_final():
    path = os.path.join(os.getcwd(), "WithKernelAlpha0")
    rho_start = 0.6
    rho_end = 1.0
    alpha0_start = 1
    alpha0_end = 20
    alpha_inter = 1
    Rho0 = np.linspace(rho_start, rho_end, num=5, endpoint=True)
    Alpha0 = np.linspace(
        alpha0_start,
        alpha0_end,
        num=(alpha0_end - alpha0_start) // alpha_inter + 1,
        endpoint=True)
    for rho0 in Rho0:
        cwd0 = os.path.join(path, r"Rho0=1.0")
        # cwd0 = os.path.join(cwd0, r"Pk21_Pk22_alpha0=10_rho0=1.0")
        for alpha0 in Alpha0:
            cwd1 = os.path.join(
                cwd0, r"pk2={0:.2f}_alpha0={1:.2f}".format(rho0, alpha0))
            cwd2 = os.path.join(cwd1, "TKS")
            os.chdir(cwd2)
            print(os.path.basename(cwd1))
            snaps = glob.glob(r"decay*sA.dat")
            snaps.sort(key=lambda x: int(x[5:len(os.path.basename(x)) - 6]))
            plot_profile_final(snaps[-1], rho0, alpha0)


def main():
    profile_final()


if __name__ == '__main__':
    main()
