import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import UnivariateSpline
import h5py
from scipy.integrate import quadrature
import pickle
import pandas as pd
import warnings
from matplotlib.ticker import MultipleLocator, NullFormatter
import os

dirname = os.path.dirname(__file__)
h5pydir = os.path.join(dirname, "../h5py_dat/")


warnings.filterwarnings("ignore")

# Load interpolations
pk_sd_interp = pickle.load(open("interps/pks_sd_interp.p", "rb"))
pk_dd_interp = pickle.load(open("interps/pks_dd_interp.p", "rb"))

# paths to styles and data
plt.style.use("/Users/nolansmyth/Dropbox/kinetic_recoupling/figures/style.mplstyle")

# Constants for plotting
pk_max = 1e2
kk = np.logspace(-4, np.log10(pk_max), 500)

# Which A_recs to use for delta recoupling rate
A_recs = [1e3]
k = 10  # use one k mode for now. This is the only one I generated data for


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


def plot_delta_recoupling_rate():

    # zs used for plotting
    z_plot = np.geomspace(1e4, 2e7, int(1e5))

    for i, A_rec in reversed(list(enumerate(A_recs))):
        plt.plot(
            z_plot,
            kappa_dot_zs_arr[i](z_plot),
            label="A_rec = " + scientific_format(A_rec),
        )

    plt.plot(z_data_no_rec, kappa_dot_data_no_rec, label="no rec ")

    plt.plot([1e-3, 1e10], [1, 1], "k--")
    plt.plot([1e-3, 1e10], [1e-3, 1e-3], "k:")
    plt.plot([1e-3, 1e10], [1e3, 1e3], "k:")

    plt.xlim(1e6, 6e7)
    plt.ylim(1e-4, 1e5)

    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("z", fontsize=16)
    plt.ylabel("$\Gamma_{\mathrm{DR-DM}}$", fontsize=16)
    plt.title("Comoving Scattering Rate")
    plt.legend()

    plot_dir = "Figures/"
    filename = "Scattering_rate_delta_test.pdf"
    plt.savefig(plot_dir + filename)
    plt.clf()


# Get data for no recoupling
spline_pars = {"k": 3, "s": 0.0}
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

plot_delta_recoupling_rate()
