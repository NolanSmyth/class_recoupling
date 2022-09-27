from classy import Class
import numpy as np
import h5py
import os

dirname = os.path.dirname(__file__)
h5pydir = os.path.join(dirname, "../h5py_dat/")

A_recs = [1e8]

pk_max = 1.0e2
# maximum k for Pk
kk = np.logspace(-4, np.log10(pk_max), 500)
# redshift at which Pk is determined
z_pk = 0.0
BM_KS = ["10"]
omega0_cdm = 0.12038
f_idm_dr = 1.0
xi = 0.3
T_rec = 6.0e5
k = 10  # use one k mode for now


def getPk(classObj):
    kk = np.logspace(-4, np.log10(pk_max), 500)  # k in h/Mpc
    Pk = []  # P(k) in (Mpc/h)**3
    h = classObj.h()  # get reduced Hubble for conversions to 1/Mpc
    for k in kk:
        Pk.append(classObj.pk(k * h, z_pk) * h ** 3)  # function .pk(k,z)
    return Pk


def save_class_obj(class_obj, A_rec=""):
    model = class_obj
    data_file = h5pydir + "class_model_data_saturation" + "%.2e" % A_rec + ".hdf5"
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
    "k_output_values": str(k),
}

idrset = {
    # Dark matter/radiation parameters
    "f_idm_dr": f_idm_dr,  # Amount of dm that is interacting
    "xi_idr": xi,
    "stat_f_idr": 0.875,  # fermionic
    "nindex_idm_dr": 4.0,
    "m_idm": 1.0e3,
    # Scattering rate parameters
    "a_idm_dr": 1.0e0,
    "rec_case": 4,  # 1 = power, 2 = Theta, 3 = delta function, 4 = no recoupling
}

sigma_fac = 0.01

for A_rec in A_recs:
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

    save_class_obj(cos, A_rec)

    # clear content of cos (to reuse it for another model)
    cos.struct_cleanup()
    # reset parameters to default
    cos.empty()
