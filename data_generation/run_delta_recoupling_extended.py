from classy import Class
import numpy as np
import h5py
import os
from variables import *

dirname = os.path.dirname(__file__)
h5pydir = os.path.join(dirname, "../h5py_dat/")

A_recs = [1e12, 10 ** (13.5), 1e15, 10 ** (16.5), 1e18]
T_recs = [6e6, 6e5 * np.sqrt(10), 6e5, 6e4 * np.sqrt(10), 6e4]
A_recs = A_recs[::-1]
T_recs = T_recs[::-1]
BM_KS = [str(round(n, 3)) for n in np.geomspace(0.1, 10, 100)]


def getPk(classObj):
    Pk = []  # P(k) in (Mpc/h)**3
    h = classObj.h()  # get reduced Hubble for conversions to 1/Mpc
    for k in kk:
        Pk.append(classObj.pk(k * h, z_pk) * h ** 3)  # function .pk(k,z)
    return Pk


def save_class_obj(class_obj, A_rec=""):
    model = class_obj
    data_file = h5pydir + "delta_data_" + "%.2e" % A_rec + ".hdf5"
    with h5py.File(data_file, "w") as f:
        # Scalar group
        data = model.get_perturbations()["scalar"]
        sub_group = f.create_group("scalar")
        for i, k in enumerate(BM_KS):
            sub_sub_group = sub_group.create_group(f"k={k}")
            d = data[i]
            for key, val in d.items():
                sub_sub_group.create_dataset(key, data=val)

        # Background group
        data = model.get_background()
        sub_group = f.create_group("background")
        for key, val in data.items():
            sub_group.create_dataset(key, data=val)

        # Thermo group
        data = model.get_thermodynamics()
        sub_group = f.create_group("thermodynamics")
        for key, val in data.items():
            sub_group.create_dataset(key, data=val)

        # Power Spectrum group
        data = getPk(model)
        sub_group = f.create_group("power_spectrum")
        sub_group.create_dataset("kk", data=kk)
        sub_group.create_dataset("Pk", data=data)


commonset = {
    "omega_b": omega_b,
    "omega_cdm": omega0_cdm,
    "h": h,
    "A_s": A_s,
    "n_s": n_s,
    "tau_reio": tau_reio,
    "output": "tCl,pCl,lCl,mPk",
    "lensing": "yes",
    "P_k_max_1/Mpc": pk_max,
    "z_pk": z_pk,
    # "k_output_values": k_output_values,
    "k_output_values": ",".join(BM_KS),
}

idrset = {
    # Dark matter/radiation parameters
    "f_idm_dr": f_idm_dr,  # Amount of dm that is interacting
    "xi_idr": xi_idr,
    "stat_f_idr": stat_f_idr,  # fermionic
    "nindex_idm_dr": nindex_idm_dr,
    "m_idm": m_idm,
    # Scattering rate parameters
    "a_idm_dr": 1e-3,
    "rec_case": 4,  # 1 = power, 2 = Theta, 3 = delta function, 4 = no recoupling
}

for i, A_rec in enumerate(A_recs):
    cos = Class()
    cos.set(commonset)
    cos.set(idrset)
    cos.set(
        {
            "rec_case": 3,
            "A_rec": A_recs[i],
            "T_rec": T_recs[i],
            "sigma": sigma_fac * T_recs[i],  # Gaussian width
        }
    )
    cos.compute()

    save_class_obj(cos, A_rec)

    # clear content of cos (to reuse it for another model)
    cos.struct_cleanup()
    # reset parameters to default
    cos.empty()
