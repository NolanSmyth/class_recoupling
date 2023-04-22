import numpy as np
import glob
from scipy.interpolate import RegularGridInterpolator
import pickle

# load dd data
data_delta_dir = "./data_delta/"
npz_delta_files = glob.glob(data_delta_dir + "*.npz")
npz_delta_files.sort()

pk_max = 1e2
kk = np.logspace(-4, np.log10(pk_max), 500)

pks_dd = []
T_recs = []
A_recs = []
T_A_recs = []

# get delta data
for npzfile in npz_delta_files:
    data = np.load(npzfile)
    pk = data["Pk"]
    pks_dd.append(pk)
    T_recs.append(data["T_rec"].item())
    A_recs.append(data["A_rec"].item())

# Returns sorted unique list of each
T_recs_unique = np.unique(T_recs)
A_recs_unique = np.unique(A_recs)

T_recs_sorted = np.flip(sorted(T_recs))
# T_recs_sorted = sorted(T_recs)

pks_dd_sorted = [x for _, x in sorted(zip(T_recs, pks_dd))]

# print(np.log10(np.min(T_recs_unique)), np.log10(np.max(T_recs_unique)))
# print(np.log10(np.min(A_recs_unique)), np.log10(np.max(A_recs_unique)))

# populate pk grid such that pk is a function of T_rec

pk_delta_dat = np.zeros((len(T_recs_sorted), len(kk)))
for i, T in enumerate(T_recs_unique):
    pk_delta_dat[i, :] = pks_dd_sorted[np.where(T_recs_sorted == T)[0][0]]

# Create dd interpolation
pks_delta_interp = RegularGridInterpolator((T_recs_unique, kk), pk_delta_dat[::-1])

# Create dd interpolation

pickle.dump(pks_delta_interp, open("interps/pks_delta_interp_new.p", "wb"))
