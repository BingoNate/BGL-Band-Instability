"""This version is used to plot last decay-profiles 
      in the directory containing the file.
"""
import os
import pathlib
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def plot_profile(snap, p12, p22):
    filebase = os.path.splitext(snap)[0]
    length = len(filebase)
    frame = int(filebase[5:length - 2])
    df = np.array(pd.read_csv(snap, delimiter=",", header=None))

    fig, axs = plt.subplots(2, 1)
    axs[0].plot(df[:, 0], df[:, 1], "b-", lw=2, label=r"$\rho$")
    axs[0].set_xlim(0, np.max(df[:, 0]))
    axs[0].legend(loc=0, fontsize=14)

    axs[1].plot(df[:, 0], df[:, 2], "r-", lw=2, label=r"$f_1^R$")
    # axs[1].plot(df[:, 0], df[:, 3], "g-", lw=2, label=r"$f_1^I$")
    axs[1].legend(loc=0, fontsize=14)
    axs[1].set_xlim(0, np.max(df[:, 0]))
    axs[1].set_xlabel(r"$X$", fontsize=18)
    # plt.suptitle(r"$\alpha_0={0:.2f}$".format(alpha0), fontsize=18)
    fig.savefig(
        r"../pk12={0:.2f}_pk22={1:.6f}.png".format(p12, p22),
        dpi=300,
        bbox_inches='tight',
        pad_inches=0)
    plt.pause(1)
    plt.close()


def profiles():
    pathname = r"WithKernelAlpha0"
    file_dir = os.path.join(os.getcwd(), pathname)
    file_dir = os.path.join(file_dir, "TKS")
    # file_dir = pathlib.Path.joinpath(pathlib.Path.cwd(), "WithoutKernel", "Rho0=1.0", "TKS")
    os.chdir(file_dir)

    snaps = glob.glob(r"decay*sA.dat")
    snaps.sort(key=lambda x: int(x[5:len(os.path.basename(x)) - 6]))

    # pk12 = pathname.split("/")[-1].split("_")[0].split("=")[-1]
    # pk22 = pathname.split("/")[-1].split("_")[-1].split("=")[-1]
    pk12 = 0.8
    pk22 = 1.0
    plot_profile(snaps[-1], float(pk12), float(pk22))


def main():
    profiles()


if __name__ == '__main__':
    main()
