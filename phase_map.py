"""This script is for plotting the phase map 
   depending on the band profile.
"""
import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def phase():
    path = os.path.join(os.getcwd(), "WithoutKernel")
    pk_start = 0.6
    pk_end = 1.0
    Pk = np.linspace(
        pk_start, pk_end, num=(pk_end - pk_start) // 0.10 + 1, endpoint=True)
    alpha1 = np.linspace(20, 200, num=10, endpoint=True)
    alpha2 = np.linspace(300, 1000, num=8, endpoint=True)
    alpha3 = np.linspace(25, 75, num=11, endpoint=True)
    alpha4 = np.linspace(2, 18, num=9, endpoint=True)
    Alpha0 = np.hstack((alpha1, alpha2))
    Alpha0 = np.hstack((Alpha0, alpha3))
    Alpha0 = np.hstack((Alpha0, alpha4))
    for pk in Pk:
        cwd0 = os.path.join(path, r"Rho0=1.0")
        cwd0 = os.path.join(cwd0, r"same_distribution")
        for alpha0 in alpha4:
            cwd1 = os.path.join(
                cwd0, r"pk22={0:.2f}_alpha0={1:.2f}".format(pk, alpha0))
            cwd2 = os.path.join(cwd1, "TKS")
            os.chdir(cwd2)
            print(os.path.basename(cwd1))
            snaps = glob.glob(r"decay*sA.dat")
            snaps.sort(key=lambda x: int(x[5:len(os.path.basename(x)) - 6]))
            df = np.array(pd.read_csv(snaps[-1], delimiter=",", header=None))
            band = np.max(df[:, 1]) - np.min(df[:, 1])
            if band <= 0.0005:
                plt.scatter(pk, alpha0, color="red")
            else:
                plt.scatter(pk, alpha0, color="blue")
    plt.xlim(0.5, 1.1)
    plt.ylim(0, 20)
    plt.xlabel(r"$P_{2}^{(2)}$", fontsize=16)
    plt.ylabel(r"$\alpha_0$", fontsize=16)
    plt.tight_layout()
    plt.savefig("../../../phasemap.png", dpi=300)
    plt.show()


def phase_pk12_pk22():
    path = os.path.join(os.getcwd(), "WithKernelAlpha0")
    Pk = np.linspace(0.6, 1.0, num=5, endpoint=True)
    Alpha0 = np.linspace(1, 20.0, num=20, endpoint=True)
    for pk in Pk:
        cwd0 = os.path.join(path, r"Rho0=1.0")
        # cwd0 = os.path.join(cwd0, r"Pk21_Pk22_alpha0=10_rho0=1.0")
        for alpha0 in Alpha0:
            cwd1 = os.path.join(
                cwd0, r"pk2={0:.2f}_alpha0={1:.2f}".format(pk, alpha0))
            cwd2 = os.path.join(cwd1, "TKS")
            os.chdir(cwd2)
            print(os.path.basename(cwd1))
            snaps = glob.glob(r"decay*sA.dat")
            snaps.sort(key=lambda x: int(x[5:len(os.path.basename(x)) - 6]))
            df = np.array(pd.read_csv(snaps[-1], delimiter=",", header=None))
            band = np.max(df[:, 1]) - np.min(df[:, 1])
            if band <= 0.0005:
                plt.scatter(pk, alpha0, color="red")
            else:
                plt.scatter(pk, alpha0,  color="blue")
    plt.xlim(0.58, 1.02)
    plt.ylim(0.98, 20.02)
    plt.xlabel(r"$\hat{P}_{2}$", fontsize=16)
    plt.ylabel(r"$\alpha_0$", fontsize=16)
    plt.tight_layout()
    plt.savefig("../../../phasemap.png", dpi=300)
    plt.show()


def main():
    # phase()
    phase_pk12_pk22()


if __name__ == '__main__':
    main()
