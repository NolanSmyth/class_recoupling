import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import UnivariateSpline
import h5py
from scipy.integrate import quadrature
import pickle
import pandas as pd

# Load interpolations
pk_sd_interp = pickle.load(open("interps/pks_sd_interp.p", "rb"))
pk_dd_interp = pickle.load(open("interps/pks_dd_interp.p", "rb"))

plt.style.use("/Users/nolansmyth/Dropbox/kinetic_recoupling/figures/style.mplstyle")

# Constants for plotting
pk_max = 1e2
kk = np.logspace(-4, np.log10(pk_max), 500)


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
    plt.title("Double Decoupling")
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


def plot_observations(Tr, Ar, best_fit_a, mwarm):
    path = "observation_data/"

    mwarm_path = "output/warm" + mwarm + ".dat"

    try:
        dfWarm = pd.read_csv(
            mwarm_path, header=None, names=["k", "P(k)"], skiprows=4, delimiter="\s+",
        )
    except FileNotFoundError:
        print("ERROR: No file found for that warm dm mass")
        return

    dfDES = pd.read_csv(path + "DESY1.csv")
    dfDES = dfDES.assign(
        ylow=dfDES["Y"] - dfDES["-DeltaY"], yhigh=dfDES["+DeltaY"] - dfDES["Y"]
    )
    dfDES = dfDES.assign(
        xlow=dfDES["X"] - dfDES["-DeltaX"], xhigh=dfDES["+DeltaX"] - dfDES["X"]
    )
    yerrDES = np.array([dfDES["ylow"], dfDES["yhigh"]])
    xerrDES = np.array([dfDES["xlow"], dfDES["xhigh"]])

    yerrDESdimless = np.array(
        [
            dfDES["ylow"] * (dfDES["X"] ** 3) / (2 * (np.pi ** 2)),
            dfDES["yhigh"] * (dfDES["X"] ** 3) / (2 * (np.pi ** 2)),
        ]
    )

    columns = ["k", "P(k)", "delta P(k)+", "delta P(k)-"]
    dfBOSS = pd.read_csv(path + "BOSS.csv", header=None)
    dfBOSS.columns = columns

    dfBOSS = dfBOSS.assign(
        ylow=dfBOSS["P(k)"] - dfBOSS["delta P(k)-"],
        yhigh=dfBOSS["delta P(k)+"] - dfBOSS["P(k)"],
    )
    yerrBOSS = np.array([dfBOSS["ylow"], dfBOSS["yhigh"]])

    yerrBOSSdimless = np.array(
        [
            dfBOSS["ylow"] * (dfBOSS["k"] ** 3) / (2 * (np.pi ** 2)),
            dfBOSS["yhigh"] * (dfBOSS["k"] ** 3) / (2 * (np.pi ** 2)),
        ]
    )

    dfHERA = pd.read_csv(path + "HeraProjected.csv")
    dfEDGES = pd.read_csv(path + "EDGESProjected.csv")

    dfHERA = dfHERA.assign(
        ylow=dfHERA["Y"] - dfHERA["-DeltaY"], yhigh=dfHERA["+DeltaY"] - dfHERA["Y"]
    )
    dfHERA = dfHERA.assign(
        xlow=dfHERA["X"] - dfHERA["-DeltaX"], xhigh=dfHERA["+DeltaX"] - dfHERA["X"]
    )
    yerrHERA = np.array([dfHERA["ylow"], dfHERA["yhigh"]])
    xerrHERA = np.array([dfHERA["xlow"], dfHERA["xhigh"]])

    dfEDGES = dfEDGES.assign(
        ylow=dfEDGES["Y"] - dfEDGES["-DeltaY"], yhigh=dfEDGES["+DeltaY"] - dfEDGES["Y"]
    )
    dfEDGES = dfEDGES.assign(
        xlow=dfEDGES["X"] - dfEDGES["-DeltaX"], xhigh=dfEDGES["+DeltaX"] - dfEDGES["X"]
    )
    yerrEDGES = np.array([dfEDGES["ylow"], dfEDGES["yhigh"]])
    xerrEDGES = np.array([dfEDGES["xlow"], dfEDGES["xhigh"]])

    # Adjust projections to agree with lcdm line
    dflcdm = pd.read_csv(
        "output/lambdacdm00_pk.dat",
        header=None,
        names=["k", "P(k)"],
        skiprows=4,
        delimiter="\s+",
    )

    dfEDGES["Y"] = [
        1.05 * dflcdm["P(k)"][121] * (dflcdm["k"][121] ** 3) / (2 * np.pi ** 2),
        dflcdm["P(k)"][129] * (dflcdm["k"][129] ** 3) / (2 * np.pi ** 2),
    ]
    dfHERA["Y"] = [
        1.05 * dflcdm["P(k)"][121] * (dflcdm["k"][121] ** 3) / (2 * np.pi ** 2),
        dflcdm["P(k)"][128] * (dflcdm["k"][128] ** 3) / (2 * np.pi ** 2),
        dflcdm["P(k)"][130] * (dflcdm["k"][130] ** 3) / (2 * np.pi ** 2),
    ]

    plt.plot(
        dflcdm["k"],
        dflcdm["P(k)"] * (dflcdm["k"] ** 3) / (2 * np.pi ** 2),
        "k",
        label=r"$\Lambda_{CDM}$",
    )

    plt.errorbar(
        dfBOSS["k"],
        dfBOSS["P(k)"] * (dfBOSS["k"]) ** 3 / (2 * np.pi ** 2),
        yerr=yerrBOSSdimless,
        marker="o",
        ms=5,
        color="mediumpurple",
        ls="none",
        label=r"BOSS DR9 Ly-$\alpha$ forest",
    )
    plt.errorbar(
        dfDES["X"],
        dfDES["Y"] * (dfDES["X"] ** 3) / (2 * np.pi ** 2),
        yerr=yerrDESdimless,
        xerr=xerrDES,
        marker="o",
        ms=5,
        color="goldenrod",
        ls="none",
        capsize=6,
        label="DES Y1 cosmic Shear",
    )

    # No conversion because these are already dimensionless
    plt.errorbar(
        dfEDGES["X"],
        dfEDGES["Y"],
        xerr=xerrEDGES,
        yerr=yerrEDGES,
        marker="^",
        ms=10,
        color="r",
        ls="none",
        label=r"EDGES 21-cm proj.",
    )
    plt.errorbar(
        dfHERA["X"],
        dfHERA["Y"],
        xerr=xerrHERA,
        yerr=yerrHERA,
        marker="^",
        ms=10,
        color="teal",
        ls="none",
        label=r"HERA 21-cm proj.",
    )

    plt.plot(
        kk,
        pk_dd_interp((Tr, Ar, kk)) * (kk ** 3) / (2 * np.pi ** 2),
        "b",
        label="Double Decoupling",
    )
    plt.plot(
        kk[-200:],
        pk_sd_interp((best_fit_a, kk[-200:])) * (kk[-200:] ** 3) / (2 * np.pi ** 2),
        "g-.",
        label="Single Decoupling",
    )

    # plt.plot(dfWarmDig['k'], dfWarmDig['delta^2(k)'], 'r--', label='Warm DM')
    plt.plot(
        dfWarm["k"],
        dfWarm["P(k)"] * (dfWarm["k"] ** 3) / (2 * np.pi ** 2),
        "r--",
        label="Warm DM, m=80 keV",
    )

    plt.xscale("log")
    plt.yscale("log")
    plt.xlim(1, 100)
    plt.ylim(2, 6e1)
    plt.xlabel("k [h/Mpc]")
    plt.ylabel("$\Delta^2_m(k)$")
    plt.title("Dimensionless Power Spectrum")

    if mwarm == "30kev":
        plt.legend(loc="upper left")
    else:
        plt.legend(loc="lower right")

    plot_dir = "Figures/"
    filename = "Power_spectrum{:.1e}{:.1e}.pdf".format(Tr, Ar)
    plt.savefig(plot_dir + filename)
    plt.clf()

