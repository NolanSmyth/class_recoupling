import numpy as np
import glob
from scipy.interpolate import RegularGridInterpolator
import pickle

# load dd data
data_delta_dir = "./data_delta/"
npz_delta_files = glob.glob(data_delta_dir + "*.npz")
npz_delta_files.sort()

T_recs = []
A_recs = []
zs = []
dmus = []

# get delta data
for npzfile in npz_delta_files:
    data = np.load(npzfile)
    T_recs.append(data["T_rec"].item())
    A_recs.append(data["A_rec"].item())
    zs.append(data["zs"])
    dmus.append(data["dmu_idm_dr"])

# Returns sorted unique list of each
T_recs_unique = np.unique(T_recs)
A_recs_unique = np.unique(A_recs)

T_recs_sorted = sorted(T_recs)

zs_sorted = [x for _, x in sorted(zip(T_recs, zs))]
dmus_sorted = [x for _, x in sorted(zip(T_recs, dmus))]

# print(T_recs_sorted, zs_sorted[0], dmus_sorted[0])

print(np.log10(np.min(T_recs_unique)), np.log10(np.max(T_recs_unique)))
pks_delta_interp = RegularGridInterpolator((T_recs_sorted, zs_sorted[0]), dmus_sorted)

# Create interpolation

pickle.dump(pks_delta_interp, open("interps/pks_delta_interp_rate.p", "wb"))
