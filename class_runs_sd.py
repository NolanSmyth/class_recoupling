from classy import Class
import matplotlib.pyplot as plt
import numpy as np
import sys


def getPk(classObj):
    kk = np.logspace(-4, np.log10(pk_max), 500)  # k in h/Mpc
    Pk = []  # P(k) in (Mpc/h)**3
    h = classObj.h()  # get reduced Hubble for conversions to 1/Mpc
    for k in kk:
        Pk.append(classObj.pk(k * h, z_pk) * h ** 3)  # function .pk(k,z)
    return Pk


def getdmu_idm_dr(classObj):
    return classObj.get_thermodynamics()["dmu_idm_dr"]


def getzs(classObj):
    return classObj.get_thermodynamics()["z"]


z_pk = 0.0
# redshift at which Pk is determined
pk_max = 1.0e2
# maximum k for Pk
kk = np.logspace(-4, np.log10(pk_max), 500)
f_idm_dr = 1.0
omega0_cdm = 0.12038

N_points = 50  # Number of grid points

a_idm_dr_arr = np.logspace(-5, 5, N_points)
a_idm_dr_idx = int(sys.argv[1])


commonset = {
    "omega_b": 0.022032,
    "omega_cdm": omega0_cdm,
    "h": 0.67556,
    "A_s": 2.215e-9,
    "n_s": 0.9619,
    "tau_reio": 0.0925,
    "output": "tCl,pCl,lCl,mPk",
    "lensing": "yes",
    "P_k_max_1/Mpc": pk_max,
    "z_pk": 0.0,
}

idrset = {
    # Dark matter/radiation parameters
    "f_idm_dr": f_idm_dr,  # Amount of dm that is interacting
    "xi_idr": 0.3,
    "stat_f_idr": 0.875,  # fermionic
    "nindex_idm_dr": 4.0,
    "m_idm": 1.0e3,
    # Scattering rate parameters
    "a_idm_dr": a_idm_dr_arr[a_idm_dr_idx],
    "rec_case": 3,  # 1 = power, 2 = Theta, 3 = no recoupling
}

idrNoRec = Class()
idrNoRec.set(commonset)
idrNoRec.set(idrset)

idrNoRec.compute()

zs = getzs(idrNoRec)
Pk = getPk(idrNoRec)
dmu_idm_dr = getdmu_idm_dr(idrNoRec)

np.savez(
    "./data_sd/run_sd_" + str(sys.argv[1]) + ".npz",
    # zs=zs,
    Pk=Pk,
    # dmu_idm_dr=dmu_idm_dr,
    a_idm_dr=a_idm_dr_arr[a_idm_dr_idx],
)
