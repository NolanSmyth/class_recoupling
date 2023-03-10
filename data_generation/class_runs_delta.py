from classy import Class
import matplotlib.pyplot as plt
import numpy as np
import sys

# This dense version runs several instances of class sqeuentially to help with the maximum allowed job array size of 1001 on the cluster
# Could be parallelized but runs should be short enough to not be a problem


def getPk(classObj):
    kk = np.logspace(-4, np.log10(pk_max), 500)  # k in h/Mpc
    Pk = []  # P(k) in (Mpc/h)**3
    h = classObj.h()  # get reduced Hubble for conversions to 1/Mpc
    for k in kk:
        Pk.append(classObj.pk(k * h, z_pk) * h**3)  # function .pk(k,z)
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
sigma_fac = 0.02

N_points = 100  # Number of grid points in each dimension

A_rec_arr = np.logspace(12, 21, N_points)
T_rec_arr = np.logspace(np.log10(6e6), np.log10(6e3), N_points)

idx = int(sys.argv[1])

T_rec = T_rec_arr[idx]
A_rec = A_rec_arr[idx]

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
    # "a_idm_dr": 1.0e0,
    "a_idm_dr": 1.0e-3,
    "rec_case": 4,  # 1 = power, 2 = Theta, 3 = delta function, 4 = no recoupling
}


cos = Class()
cos.set(commonset)
cos.set(idrset)
cos.set(
    {
        "rec_case": 3,
        "A_rec": A_rec,
        "T_rec": T_rec,
        "sigma": sigma_fac * T_rec,  # Gaussian width
    }
)
cos.compute()

zs = getzs(cos)
Pk = getPk(cos)
dmu_idm_dr = getdmu_idm_dr(cos)

np.savez(
    "./data_delta/run_" + str(sys.argv[1]) + ".npz",
    # zs=zs,
    Pk=Pk,
    # dmu_idm_dr=dmu_idm_dr,
    T_rec=T_rec,
    A_rec=A_rec,
)

# clear content of cos (to reuse it for another model)
cos.struct_cleanup()
# reset parameters to default
cos.empty()
