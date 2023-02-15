import numpy as np
import glob
from scipy.interpolate import RegularGridInterpolator
import pickle

# load dd data
data_strong_dir = "./data_dense_strong/"
npz_strong_files = glob.glob(data_strong_dir + "*.npz")
npz_strong_files.sort()

pk_max = 1e2
kk = np.logspace(-4, np.log10(pk_max), 500)

pks_dd = []
T_recs = []
A_recs = []
T_A_recs = []

# get dd data
for npzfile in npz_strong_files:
    data = np.load(npzfile)
    pk = data["Pk"]
    pks_dd.append(pk)
    T_recs.append(data["T_rec"].item())
    A_recs.append(data["A_rec"].item())
    T_A_recs.append((data["T_rec"].item(), data["A_rec"].item()))

# Returns sorted unique list of each
T_recs_unique = np.unique(T_recs)
A_recs_unique = np.unique(A_recs)

print(np.log10(np.min(T_recs_unique)), np.log10(np.max(T_recs_unique)))
print(np.log10(np.min(A_recs_unique)), np.log10(np.max(A_recs_unique)))


# populate pk grid such that pk is a function of T_rec and A_rec

pk_strong_dat = np.zeros((len(T_recs_unique), len(A_recs_unique), len(kk)))
for i, T in enumerate(T_recs_unique):
    for j, A in enumerate(A_recs_unique):
        try:
            pk_strong_dat[i, j, :] = pks_dd[T_A_recs.index((T, A))]
        except ValueError:
            continue

        # pk_strong_dat[i, j, :] = pks_dd[T_A_recs.index((T, A))]


# Create dd interpolation
pks_strong_interp = RegularGridInterpolator(
    (T_recs_unique, A_recs_unique, kk), pk_strong_dat
)

# Create dd interpolation

pickle.dump(pks_strong_interp, open("interps/pks_strong_interp_late.p", "wb"))
