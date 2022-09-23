import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import UnivariateSpline
import h5py
from scipy.integrate import quadrature
import pickle

# Load interpolations
pk_sd_interp = pickle.load(open("interps/pks_sd_interp.p", "rb"))
pk_dd_interp = pickle.load(open("interps/pks_dd_interp.p", "rb"))

plt.style.use("/Users/nolansmyth/Dropbox/kinetic_recoupling/figures/style.mplstyle")


def dmu_idm_dr(
    T_rec,
    A_rec,
    z,
    case="recoupling",
    a_idm_dr=1,
    nindex_idm_dr=4,
    omega0_cdm=0.12038,
    f_idm_dr=1.0,
    h=0.67556,
    xi=0.3,
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


def scientific_format(x):
    s = "%.1e" % x
    mantissa, exponent = s.split("e")
    return r"${} \times 10^{{{}}}$".format(mantissa, int(exponent))


def plot_varied_recoupling(Tr0, Ar0, Tr1, Tr2, Ar1, Ar2, num_interps=7):
    """
    Plot the interpolation between two scenarios with varied 
    recoupling temperature and/or strength (1,2).
    Normalized to the first scenario (0), presumably no recoupling.
    """
    T_interps = np.logspace(np.log10(Tr1), np.log10(Tr2), num_interps)
    A_interps = np.logspace(np.log10(Ar1), np.log10(Ar2), num_interps)

    zs = np.logspace(5, 8, 1000)
    dmus1 = [dmu_idm_dr(Tr1, Ar1, z) for z in zs]
    dmus2 = [dmu_idm_dr(Tr2, Ar2, z) for z in zs]

    pk_max = 1e2
    kk = np.logspace(-4, np.log10(pk_max), 500)

    fig = plt.figure(1, figsize=(8, 8))
    plt.subplot(211)
    plt.plot(
        kk,
        pk_dd_interp((Tr1, Ar1, kk)) / pk_dd_interp((Tr0, Ar0, kk)),
        "b",
        label="$T_\mathrm{rec}$=%s, $A_\mathrm{rec}$=%s"
        % (scientific_format(Tr1), scientific_format(Ar1)),
    )
    plt.plot(
        kk,
        pk_dd_interp((Tr2, Ar2, kk)) / pk_dd_interp((Tr0, Ar0, kk)),
        "r",
        label="$T_\mathrm{rec}$=%s, $A_\mathrm{rec}$=%s"
        % (scientific_format(Tr2), scientific_format(Ar2)),
    )
    for num_interp in range(num_interps):
        plt.plot(
            kk,
            pk_dd_interp((T_interps[num_interp], A_interps[num_interp], kk))
            / pk_dd_interp((Tr0, Ar0, kk)),
            "--k",
            alpha=1 / (num_interp + 1),
        )
    plt.xscale("log")
    plt.yscale("log")
    plt.xlim(1e1, 1e2)
    plt.ylim(0.97, 1)
    plt.xlabel("k [h/Mpc]")
    plt.ylabel("$P(k)/P(k)_0$")
    plt.legend()
    plt.title("Double Decoupling, Varying $A_\mathrm{rec}$")
    plt.subplot(212)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("z")
    plt.ylabel("$\Gamma_{\mathrm{DR-DM}}$")
    plt.plot(
        zs,
        dmus1,
        "b",
        label="$T_\mathrm{rec}$=%s, $A_\mathrm{rec}$=%s"
        % (scientific_format(Tr1), scientific_format(Ar1)),
    )
    plt.plot(
        zs,
        dmus2,
        "r--",
        label="$T_\mathrm{rec}$=%s, $A_\mathrm{rec}$=%s"
        % (scientific_format(Tr2), scientific_format(Ar2)),
    )
    plt.plot(np.logspace(5, 8, 100), 1e-3 * np.ones(100), "--k")
    plt.plot(np.logspace(5, 8, 100), np.ones(100), "--k")

    plt.xlim(zs[0], zs[-1])

    plt.legend()
    fig.tight_layout(h_pad=2)

    plot_dir = "Figures/"
    filename = "varying_recoupling_{:.1e}{:.1e}{:.1e}{:.1e}.pdf".format(
        Tr1, Ar1, Tr2, Ar2
    )

    plt.savefig(plot_dir + filename)
    plt.clf()
