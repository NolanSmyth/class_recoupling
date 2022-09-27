import numpy as np
import glob
from scipy.interpolate import RegularGridInterpolator
import pickle

# load dd data
data_dd_dir = "./data_dense/"
npz_dd_files = glob.glob(data_dd_dir + "*.npz")
npz_dd_files.sort()

# load sd data
data_sd_dir = "./data_sd/"
npz_sd_files = glob.glob(data_sd_dir + "*.npz")
npz_sd_files.sort()


pk_max = 1e2
kk = np.logspace(-4, np.log10(pk_max), 500)

pks_dd = []
T_recs = []
A_recs = []
T_A_recs = []

# get dd data
for npzfile in npz_dd_files:
    data = np.load(npzfile)
    pk = data["Pk"]
    pks_dd.append(pk)
    T_recs.append(data["T_rec"].item())
    A_recs.append(data["A_rec"].item())
    T_A_recs.append((data["T_rec"].item(), data["A_rec"].item()))

# Returns sorted unique list of each
T_recs_unique = np.unique(T_recs)
A_recs_unique = np.unique(A_recs)

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

pk_dd_dat = np.zeros((len(T_recs_unique), len(A_recs_unique), len(kk)))
for i, T in enumerate(T_recs_unique):
    for j, A in enumerate(A_recs_unique):
        pk_dd_dat[i, j, :] = pks_dd[T_A_recs.index((T, A))]

# Create dd interpolation
pks_dd_interp = RegularGridInterpolator((T_recs_unique, A_recs_unique, kk), pk_dd_dat)

# Create sd interpolation
pks_sd_interp = RegularGridInterpolator((a_idm_drs_unique, kk), pk_sd_dat)

pickle.dump(pks_dd_interp, open("interps/pks_dd_interp.p", "wb"))

pickle.dump(pks_sd_interp, open("interps/pks_sd_interp.p", "wb"))
