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

N_points = 100  # Number of grid points in each dimension

T_rec_arr = np.logspace(4, 7, N_points)
A_rec_arr = np.logspace(-1, 7, N_points)

T_rec_idx = int(sys.argv[1])
# A_rec_idx = int(sys.argv[1]) // len(T_rec_arr)

# print("T_rec_idx: ", T_rec_idx, "A_rec_idx: ", A_rec_idx)

T_rec = T_rec_arr[T_rec_idx]

# print("T_rec: ", T_rec, "A_rec: ", A_rec)

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
    "nindex_idm_dr": 2.0,
    "m_idm": 1.0e3,
    # Scattering rate parameters
    # "a_idm_dr": 1.0e0,
    "a_idm_dr": 1.0e1,
    "rec_case": 3,  # 1 = power, 2 = Theta, 3 = no recoupling
}

for A_rec_idx in range(len(A_rec_arr)):

    recset = {
        "T_rec": T_rec,
        "A_rec": A_rec_arr[A_rec_idx],
        "rec_case": 2,
    }

    idrRec = Class()
    idrRec.set(commonset)
    idrRec.set(idrset)
    idrRec.set(recset)

    idrRec.compute()

    zs = getzs(idrRec)
    Pk = getPk(idrRec)
    dmu_idm_dr = getdmu_idm_dr(idrRec)

    np.savez(
        "./data_dense_n2/run_" + str(sys.argv[1]) + "_" + str(A_rec_idx) + ".npz",
        # zs=zs,
        Pk=Pk,
        # dmu_idm_dr=dmu_idm_dr,
        T_rec=T_rec,
        A_rec=A_rec_arr[A_rec_idx],
    )
    idrRec.empty()
