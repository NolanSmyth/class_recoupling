import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import UnivariateSpline
import h5py
from scipy.integrate import quad
import pickle
import pandas as pd
import warnings
from data_generation.variables import *
import matplotlib.ticker as plticker

warnings.filterwarnings("ignore")

# Load interpolations
pk_sd_interp = pickle.load(open("interps/pks_sd_interp.p", "rb"))
pk_dd_interp = pickle.load(open("interps/pks_dd_interp.p", "rb"))

# paths to styles and data
plt.style.use("Figures/style.mplstyle")
h5pydir = "h5py_dat/"

# Which A_recs to use for delta recoupling rate
A_recs = [1e11, 1e12, 1e14]


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


def scientific_format(x):
    s = "%.1e" % x
    mantissa, exponent = s.split("e")
    return r"${} \times 10^{{{}}}$".format(mantissa, int(exponent))


def plot_varied_recoupling(Tr0, Ar0, Tr1, Tr2, Ar1, Ar2, num_interps=7, save=True):
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

    if save:
        plot_dir = "Figures/"
        filename = "varying_recoupling_{:.1e}{:.1e}{:.1e}{:.1e}.pdf".format(
            Tr1, Ar1, Tr2, Ar2
        )

        plt.savefig(plot_dir + filename)
        plt.clf()
    else:
        plt.show()


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
        "output/lambdacdm.dat",
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

    plt.plot(
        dfWarm["k"],
        dfWarm["P(k)"] * (dfWarm["k"] ** 3) / (2 * np.pi ** 2),
        "r--",
        label="Warm DM, m={} keV".format(mwarm.split("k")[0]),
    )

    plt.xscale("log")
    plt.yscale("log")
    plt.xlim(1, 100)
    plt.ylim(2, 6e1)
    plt.xlabel("k [h/Mpc]")
    plt.ylabel("$\Delta^2_m(k)$")
    plt.title("Dimensionless Power Spectrum")

    plt.legend(loc="upper left")

    plot_dir = "Figures/"
    filename = "Power_spectrum{:.1e}{:.1e}.pdf".format(Tr, Ar)
    plt.savefig(plot_dir + filename)
    plt.clf()


def plot_delta_recoupling_rate():
    for i, A_rec in reversed(list(enumerate(A_recs))):
        plt.plot(
            tau_arr[i],
            kappa_dot_taus_arr[i](tau_arr[i]),
            label="A_rec = " + scientific_format(A_rec),
        )

    plt.plot([1e-3, 1e10], [1, 1], "k--")
    plt.plot([1e-3, 1e10], [1e-3, 1e-3], "k:")
    plt.plot([1e-3, 1e10], [1e3, 1e3], "k:")

    plt.xlim(2e-1, 2e-0)
    plt.ylim(1e-8, 1e8)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("$\\tau$ [Mpc]", fontsize=16)

    plt.ylabel("$\Gamma_{\mathrm{DR-DM}}$", fontsize=16)
    plt.title("Comoving Scattering Rate")
    plt.legend()
    plot_dir = "Figures/"
    filename = "Scattering_rate_delta.pdf"
    plt.savefig(plot_dir + filename)
    plt.clf()


# Get data for no recoupling
data_file = h5pydir + "class_model_data_no_rec.hdf5"
with h5py.File(data_file, "r") as f:
    tau_data = np.array(f["scalar"]["k=" + str(k)]["tau [Mpc]"])
    delta_chi_data = np.array(f["scalar"]["k=" + str(k)]["delta_idm_dr"])
    phi_data = np.array(f["scalar"]["k=" + str(k)]["phi"])
    theta_dr_data = np.array(f["scalar"]["k=" + str(k)]["theta_idr"])
    theta_data = np.array(f["scalar"]["k=" + str(k)]["theta_idm_dr"])
    Pks_no_rec = np.array(f["power_spectrum"]["Pk"])
    kks_no_rec = np.array(f["power_spectrum"]["kk"])
    kappa_dot_data_no_rec = np.flip(np.array(f["thermodynamics"]["dmu_idm_dr"]))
    z_data_no_rec = np.flip(np.array(f["thermodynamics"]["z"]))
    rho_idr_no_rec = np.array(f["background"]["(.)rho_idr"])
    rho_idm_no_rec = np.array(f["background"]["(.)rho_idm_dr"])
    bkg_zs_no_rec = np.array(f["background"]["z"])

delta_chi_no_rec = UnivariateSpline(tau_data, delta_chi_data, **spline_pars)
delta_chi_dot_no_rec = delta_chi_no_rec.derivative()
phi_no_rec = UnivariateSpline(tau_data, phi_data, **spline_pars)
phi_dot_no_rec = phi_no_rec.derivative()
theta_dr_no_rec = UnivariateSpline(tau_data, theta_dr_data, **spline_pars)
theta_chi_no_rec = UnivariateSpline(tau_data, theta_data, **spline_pars)
theta_chi_dot_no_rec = theta_chi_no_rec.derivative()

# Get data for recouplings
delta_chi_arr = []
delta_chi_dot_arr = []
delta_chi_ddot_arr = []
phi_dot_arr = []
phi_ddot_arr = []
Pk_arr = []
kappa_dot_taus_arr = []
kappa_dot_zs_arr = []
z_arr = []
tau_arr = []
thermo_tau_arr = []
theta_dr_arr = []
rho_idr_arr = []
rho_idm_arr = []
theta_chi_dot_arr = []
theta_chi_arr = []
cx2_chi_arr = []
a_prime_arr = []
a_arr = []
psi_arr = []

for A_rec in A_recs:
    data_file = h5pydir + "class_model_data_" + "%.2e" % A_rec + ".hdf5"

    with h5py.File(data_file, "r") as f:
        tau_data = np.array(f["scalar"]["k=" + str(k)]["tau [Mpc]"])
        delta_chi_data = np.array(f["scalar"]["k=" + str(k)]["delta_idm_dr"])
        phi_data = np.array(f["scalar"]["k=" + str(k)]["phi"])
        kappa_dot_data = np.flip(np.array(f["thermodynamics"]["dmu_idm_dr"]))
        z_data = np.flip(np.array(f["thermodynamics"]["z"]))
        thermo_tau_data = np.flip(np.array(f["thermodynamics"]["conf. time [Mpc]"]))
        cx2_data = np.flip(np.array(f["thermodynamics"]["c_idm_dr^2"]))
        theta_data = np.array(f["scalar"]["k=" + str(k)]["theta_idm_dr"])
        theta_dr_data = np.array(f["scalar"]["k=" + str(k)]["theta_idr"])
        psi_data = np.array(f["scalar"]["k=" + str(k)]["psi"])
        rho_idr_data = np.array(f["background"]["(.)rho_idr"])
        rho_idm_data = np.array(f["background"]["(.)rho_idm_dr"])
        bkg_tau_data = np.array(f["background"]["conf. time [Mpc]"])
        Pks = np.array(f["power_spectrum"]["Pk"])
        kks = np.array(f["power_spectrum"]["kk"])

    rho_idr = UnivariateSpline(bkg_tau_data, rho_idr_data, **spline_pars)
    rho_idm = UnivariateSpline(bkg_tau_data, rho_idm_data, **spline_pars)

    delta_chi = UnivariateSpline(tau_data, delta_chi_data, **spline_pars)
    delta_chi_dot = delta_chi.derivative()
    delta_chi_ddot = delta_chi_dot.derivative()

    phi = UnivariateSpline(tau_data, phi_data, **spline_pars)
    phi_dot = phi.derivative()
    phi_ddot = phi_dot.derivative()

    theta_chi = UnivariateSpline(tau_data, theta_data, **spline_pars)
    theta_chi_dot = theta_chi.derivative()
    theta_dr = UnivariateSpline(tau_data, theta_dr_data, **spline_pars)

    psi = UnivariateSpline(tau_data, psi_data, **spline_pars)
    cx2_chi = UnivariateSpline(thermo_tau_data, cx2_data, **spline_pars)

    kappa_dot_taus = UnivariateSpline(thermo_tau_data, kappa_dot_data, **spline_pars)
    kappa_dot_zs = UnivariateSpline(
        np.flip(z_data), np.flip(kappa_dot_data), **spline_pars
    )
    thermo_taus = UnivariateSpline(
        np.flip(z_data), np.flip(thermo_tau_data), **spline_pars
    )

    a_data = 1 / (1 + z_data)
    a = UnivariateSpline(thermo_tau_data, a_data, **spline_pars)
    aprime = a.derivative()

    delta_chi_arr.append(delta_chi)
    delta_chi_dot_arr.append(delta_chi_dot)
    delta_chi_ddot_arr.append(delta_chi_ddot)
    Pk_arr.append(Pks)
    kappa_dot_taus_arr.append(kappa_dot_taus)
    kappa_dot_zs_arr.append(kappa_dot_zs)
    z_arr.append(z_data)
    tau_arr.append(tau_data)
    phi_dot_arr.append(phi_dot)
    phi_ddot_arr.append(phi_ddot)
    thermo_tau_arr.append(thermo_taus)
    theta_dr_arr.append(theta_dr)
    rho_idr_arr.append(rho_idr)
    rho_idm_arr.append(rho_idm)
    theta_chi_arr.append(theta_chi)
    theta_chi_dot_arr.append(theta_chi_dot)
    cx2_chi_arr.append(cx2_chi)
    a_prime_arr.append(aprime)
    a_arr.append(a)
    psi_arr.append(psi)


def lhs221NonInt(tau, idx):
    return (
        3 * phi_ddot_arr[idx](tau)
        - delta_chi_ddot_arr[idx](tau)
        - cx2_chi_arr[idx](tau) * k ** 2 * delta_chi_arr[idx](tau)
        + a_prime_arr[idx](tau) / a_arr[idx](tau) * theta_chi_arr[idx](tau)
        - k ** 2 * psi_arr[idx](tau)
    ) / (
        3 * phi_dot_arr[idx](tau) - delta_chi_dot_arr[idx](tau) - theta_dr_arr[idx](tau)
    )


def rhs221NonInt(tau, idx):
    return (
        4
        / 3
        * rho_idr_arr[idx](tau)
        / rho_idm_arr[idx](tau)
        * -1
        * kappa_dot_taus_arr[idx](tau)
    )


def lhs221(tau, idx):
    return quad(lhs221NonInt, 6e-1, tau, args=(idx), epsrel=1e-8, limit=20)[0]


def rhs221(tau, idx):
    return quad(rhs221NonInt, 6e-1, tau, args=(idx), epsrel=1e-8, limit=20)[0]


def plot_delta_effect():
    taus = np.linspace(5.6e-1, 7e-1, 50)

    lhsarr0 = np.array([lhs221(t, 0) for t in taus])
    rhsarr0 = np.array([rhs221(t, 0) for t in taus])

    lhsarr1 = np.array([lhs221(t, 1) for t in taus])
    rhsarr1 = np.array([rhs221(t, 1) for t in taus])

    lhsarr2 = np.array([lhs221(t, 2) for t in taus])
    rhsarr2 = np.array([rhs221(t, 2) for t in taus])

    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    plt.loglog(
        taus,
        -1 * lhsarr0,
        # label="LHS, A_rec = " + scientific_format(A_recs[0]),
        color="k",
    )
    plt.loglog(
        taus,
        -1 * rhsarr0,
        "--",
        label="RHS, A_rec = " + scientific_format(A_recs[0]),
        color=colors[2],
    )

    plt.loglog(
        taus,
        -1 * lhsarr1,
        # label="LHS, A_rec = " + scientific_format(A_recs[1]),
        color="k",
    )
    plt.loglog(
        taus,
        -1 * rhsarr1,
        "--",
        label="RHS, A_rec = " + scientific_format(A_recs[1]),
        color=colors[1],
    )

    plt.loglog(
        taus,
        -1 * lhsarr2,
        # label="LHS, A_rec = " + scientific_format(A_recs[2]),
        color="k",
    )
    plt.loglog(
        taus,
        -1 * rhsarr2,
        "--",
        label="RHS, A_rec = " + scientific_format(A_recs[2]),
        color=colors[0],
    )

    plt.title("Eq 2.21 from draft")
    plt.xlabel("$\\tau$ [Mpc]")
    plt.xlim(5.95e-1, 7e-1)
    plt.ylim(1e-2, 1e4)

    plt.legend(loc="upper left")
    plot_dir = "Figures/"
    filename = "Delta_effect.pdf"
    plt.savefig(plot_dir + filename)
    plt.clf()


def plot_delta_power_spectrum():
    lines = ["-", "--", "-."]
    for i, A_rec in reversed(list(enumerate(A_recs))):
        plt.plot(
            kk,
            Pk_arr[i] / Pks_no_rec,
            ls=lines[i % len(lines)],
            label="A_rec = " + scientific_format(A_rec),
        )

    plt.plot(kk, Pks_no_rec / Pks_no_rec, "--", label="No Rec")
    plt.yscale("log")
    plt.xscale("log")
    plt.xlabel("k")
    plt.ylabel("$P(k)/P(k)_0$")
    plt.title("Matter Power Spectrum Ratio")
    plt.xlim(1, 1e2)
    plt.ylim(1e-3, 2)
    plt.legend()
    plot_dir = "Figures/"
    filename = "delta_power_spectrum.pdf"
    plt.savefig(plot_dir + filename)
    plt.clf()


def plot_delta_power_spectrum_dimless():

    lines = ["-", "--", "-."]
    for i, A_rec in reversed(list(enumerate(A_recs))):
        plt.plot(
            kk,
            Pk_arr[i] * kk ** 3 / (2 * np.pi ** 2),
            ls=lines[i % len(lines)],
            label="A_rec = " + scientific_format(A_rec),
        )

    plt.plot(kk, Pks_no_rec * kk ** 3 / (2 * np.pi ** 2), "--", label="No Rec")
    plt.yscale("log")
    plt.xscale("log")
    plt.xlabel("k")
    plt.ylabel("$\Delta^2_m(k)$")
    plt.title("Dimensionless Matter Power Spectrum - $\delta$ recoupling")
    plt.xlim(1, 1e2)
    plt.ylim(1e-3, 3e1)
    plt.legend()
    plot_dir = "Figures/"
    filename = "delta_power_spectrum_dimless.pdf"
    plt.savefig(plot_dir + filename)
    plt.clf()


def plot_delta_power_spectra_both():
    lines = ["-", "--", "-."]

    plt.figure(1, figsize=(8, 8))
    plt.subplot(211)

    for i, A_rec in reversed(list(enumerate(A_recs))):
        plt.plot(
            kk,
            Pk_arr[i] * kk ** 3 / (2 * np.pi ** 2),
            ls=lines[i % len(lines)],
            label="A_rec = " + scientific_format(A_rec),
        )

    plt.plot(kk, Pks_no_rec * kk ** 3 / (2 * np.pi ** 2), "--", label="No Rec")
    plt.yscale("log")
    plt.xscale("log")
    plt.ylabel("$\Delta^2_m(k)$")
    plt.title("Dimensionless Matter Power Spectrum - $\delta$ recoupling")
    plt.xlim(1, 1e2)
    plt.ylim(1e-3, 3e1)
    plt.legend()

    plt.subplot(212)

    for i, A_rec in reversed(list(enumerate(A_recs))):
        plt.plot(
            kk,
            Pk_arr[i] / Pks_no_rec,
            ls=lines[i % len(lines)],
            label="A_rec = " + scientific_format(A_rec),
        )

    plt.plot(kk, Pks_no_rec / Pks_no_rec, "--", label="No Rec")
    plt.yscale("log")
    plt.xscale("log")
    plt.xlabel("k [h/Mpc]")
    plt.ylabel("$P(k)/P(k)_0$")
    plt.title("Ratio Compared to No Recoupling")
    plt.xlim(1, 1e2)
    plt.ylim(1e-3, 2)

    # plt.legend()
    plot_dir = "Figures/"
    filename = "delta_power_spectra.pdf"
    plt.savefig(plot_dir + filename)
    plt.clf()


def plot_delta_effect_both():

    taus = np.linspace(5.6e-1, 7.05e-1, 50)

    lhsarr0 = np.array([lhs221(t, 0) for t in taus])
    rhsarr0 = np.array([rhs221(t, 0) for t in taus])

    lhsarr1 = np.array([lhs221(t, 1) for t in taus])
    rhsarr1 = np.array([rhs221(t, 1) for t in taus])

    lhsarr2 = np.array([lhs221(t, 2) for t in taus])
    rhsarr2 = np.array([rhs221(t, 2) for t in taus])

    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    fig = plt.figure(1, figsize=(8, 8))

    for i, A_rec in reversed(list(enumerate(A_recs))):
        plt.plot(
            tau_arr[i],
            kappa_dot_taus_arr[i](tau_arr[i]),
            label="A_rec = " + scientific_format(A_rec),
        )

    plt.plot([1e-3, 1e10], [1, 1], "k--")
    plt.plot([1e-3, 1e10], [1e-3, 1e-3], "k:")
    plt.plot([1e-3, 1e10], [1e3, 1e3], "k:")

    plt.xlim(1e-1, 1e-0)
    plt.ylim(1e-8, 1e8)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("$\\tau$ [Mpc]", fontsize=16)
    plt.legend(loc="lower left")

    plt.ylabel("$\Gamma_{\mathrm{DR-DM}}$", fontsize=16)
    plt.title("Comoving Scattering Rate")

    left, bottom, width, height = [0.25, 0.58, 0.25, 0.25]
    ax2 = fig.add_axes([left, bottom, width, height])

    plt.loglog(
        taus,
        -1 * lhsarr0,
        # label="LHS, A_rec = " + scientific_format(A_recs[0]),
        color="k",
    )
    plt.loglog(
        taus,
        -1 * rhsarr0,
        "--",
        label="RHS, A_rec = " + scientific_format(A_recs[0]),
        color=colors[2],
    )

    plt.loglog(
        taus,
        -1 * lhsarr1,
        # label="LHS, A_rec = " + scientific_format(A_recs[1]),
        color="k",
    )
    plt.loglog(
        taus,
        -1 * rhsarr1,
        "--",
        label="RHS, A_rec = " + scientific_format(A_recs[1]),
        color=colors[1],
    )

    plt.loglog(
        taus,
        -1 * lhsarr2,
        # label="LHS, A_rec = " + scientific_format(A_recs[2]),
        color="k",
    )
    plt.loglog(
        taus,
        -1 * rhsarr2,
        "--",
        label="RHS, A_rec = " + scientific_format(A_recs[2]),
        color=colors[0],
    )

    # plt.title("Eq 2.21 from draft")
    plt.xlabel("$\\tau$", fontsize=16)
    plt.xlim(5.95e-1, 7.05e-1)
    ax2.set_xticks([6e-1, 7e-1])
    ax2.get_xaxis().set_major_formatter(plticker.ScalarFormatter())
    ax2.get_xaxis().set_minor_formatter(plticker.NullFormatter())

    plt.ylim(1e-2, 1e4)

    plot_dir = "Figures/"
    filename = "Delta_effect_both.pdf"
    plt.savefig(plot_dir + filename)
    plt.clf()
