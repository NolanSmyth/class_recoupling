import numpy as np
import glob
from scipy.interpolate import RegularGridInterpolator
import pickle

# load sd data
data_sd_dir = "./data_sd_n2/"
npz_sd_files = glob.glob(data_sd_dir + "*.npz")
npz_sd_files.sort()

pk_max = 1e2
kk = np.logspace(-4, np.log10(pk_max), 500)

# get sd data
pks_sd = []
a_idm_drs = []
for npzfile in npz_sd_files:
    data = np.load(npzfile)
    pk = data["Pk"]
    pks_sd.append(pk)
    a_idm_drs.append(data["a_idm_dr"])

a_idm_drs_unique = np.unique(a_idm_drs)

# populate pk grid such that pk is a function of T_rec and A_rec
pk_sd_dat = np.zeros((len(a_idm_drs), len(kk)))
for i, a_idm_dr in enumerate(a_idm_drs_unique):
    pk_sd_dat[i, :] = pks_sd[a_idm_drs.index(a_idm_dr)]

# Create sd interpolation
pks_sd_interp = RegularGridInterpolator((a_idm_drs_unique, kk), pk_sd_dat)

pickle.dump(pks_sd_interp, open("interps/pks_sd_interp_n2.p", "wb"))
