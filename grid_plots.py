import matplotlib.pyplot as plt
import numpy as np
import pickle
from data_generation.variables import *

pk_dd_interp = pickle.load(open("interps/pks_dd_interp.p", "rb"))


def scientific_format(x):
    s = "%.1e" % x
    mantissa, exponent = s.split("e")
    return r"${} \times 10^{{{}}}$".format(mantissa, int(exponent))


def dmu_idm_dr(
    T_rec,
    A_rec,
    z,
    case="recoupling",
    a_idm_dr=a_idm_dr,
    nindex_idm_dr=nindex_idm_dr,
    omega0_cdm=omega0_cdm,
    f_idm_dr=f_idm_dr,
    h=h,
    xi=xi_idr,
):
    """
    Calculate the comoving scattering rate for a given T_rec, A_rec, z.
    """
    base_rate = (
        a_idm_dr * ((1 + z) / (1e7)) ** nindex_idm_dr * omega0_cdm * f_idm_dr * h ** 2
    )
    T_idr = 2.7255 * xi

    # This is for a phase transition like scenario, not a delta function
    if case == "recoupling":
        if T_rec > T_idr * (1 + z):
            return base_rate * (1 + A_rec)

    return base_rate


def plot_varied_recoupling_grid_collapsed(Tr0, Ar0, Trs, Ars, save=True):
    """
    Plot a grid of recoupling scenarios with varied location and recoupling strength.
    Normalized to the first scenario (0), presumably no recoupling.
    """

    grid_size = max(len(Trs), len(Ars))

    fig, axes = plt.subplots(
        nrows=len(Trs),
        ncols=1,
        sharex=True,
        sharey=True,
        figsize=(grid_size * 2, grid_size * 3),
    )
    for i, Tr in enumerate(Trs):
        for j, Ar in enumerate(Ars):
            axes[i].plot(
                kk,
                (1 - (pk_dd_interp((Tr, Ar, kk)) / pk_dd_interp((Tr0, Ar0, kk))))
                / (1 - (pk_dd_interp((Tr, Ars[0], kk)) / pk_dd_interp((Tr0, Ar0, kk)))),
                label="$A_\mathrm{rec}$=%s" % (scientific_format(Ar)),
            )
            # axes[i].plot(
            #     kk,
            #     pk_dd_interp((Tr, Ar, kk)) / pk_dd_interp((Tr, Ars[0], kk)),
            #     label="$A_\mathrm{rec}$=%s" % (scientific_format(Ar)),
            # )
        axes[i].set_xscale("log")
        axes[i].set_yscale("log")
        axes[i].set_xlim(5e0, 1e2)
        axes[i].set_ylim(7e-1, 5e2)
        axes[i].set_xticks([1e1, 1e2])
        axes[i].set_title(
            "$T_\mathrm{rec}$=%s eV " % (scientific_format(Tr * ktoev)), fontsize=12,
        )
        axes[i].legend(loc="upper left")
    fig.supxlabel("$k [Mpc^{-1}$]")
    fig.supylabel("$(1 - P(k)/P(k)_0)/(1 - P(k)_\mathrm{low}/P(k)_0)$")

    plt.tight_layout()

    if save:
        plot_dir = "Figures/"
        filename = "varying_recoupling_grid_collapsed.pdf"
        plt.savefig(plot_dir + filename)
        plt.clf()
    else:
        plt.show()

    fig, axes = plt.subplots(
        nrows=len(Trs),
        ncols=1,
        sharex=True,
        sharey=True,
        figsize=(grid_size * 2, grid_size * 3),
    )

    zs = np.logspace(5, 8, 1000)
    for i, Tr in enumerate(Trs):
        for j, Ar in enumerate(Ars):
            dmus = [dmu_idm_dr(Tr, Ar, z) for z in zs]
            axes[i].plot(
                zs, dmus, label="$A_\mathrm{rec}$=%s" % (scientific_format(Ar)),
            )
            axes[i].plot(np.logspace(5, 8, 100), np.ones(100), "--k")
            axes[i].plot(np.logspace(5, 8, 100), 1e-3 * np.ones(100), "--k")

        axes[i].set_xscale("log")
        axes[i].set_yscale("log")

        axes[i].set_xlim(1e5, 1e8)
        axes[i].set_xticks([1e6, 1e7, 1e8])
        axes[i].set_ylim(1e-6, 1e3)

        axes[i].legend(loc="upper left")

        axes[i].set_title(
            "$T_\mathrm{rec}$=%s eV" % (scientific_format(Tr * ktoev)), fontsize=12,
        )
    fig.supxlabel("$z$")
    fig.supylabel("$\Gamma_{\mathrm{DM-DR}} / \mathcal{H}$")

    plt.tight_layout()

    if save:
        plot_dir = "Figures/"
        filename = "varying_recoupling_rate_grid_collapsed.pdf"
        plt.savefig(plot_dir + filename)
        plt.clf()
    else:
        plt.show()


N_points = 100
# Values over which dd interpolation is defined (This is hardcoded)
T_rec_arr = np.logspace(5, 7, N_points)
A_rec_arr = np.logspace(-1, 3, N_points)

idx = 65
Trs = [
    T_rec_arr[idx],
    T_rec_arr[idx + 15],
    T_rec_arr[idx + 30],
]
Ars = [
    A_rec_arr[idx - 10],
    A_rec_arr[idx],
    A_rec_arr[idx + 10],
    A_rec_arr[idx + 20],
    A_rec_arr[idx + 30],
]

plot_varied_recoupling_grid_collapsed(T_rec_arr[0], A_rec_arr[0], Trs, Ars)
