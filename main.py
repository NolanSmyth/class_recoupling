import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import UnivariateSpline
import h5py
from scipy.integrate import quadrature
import pickle
import plots as p

plt.style.use("/Users/nolansmyth/Dropbox/kinetic_recoupling/figures/style.mplstyle")

N_points = 100
# Values over which dd interpolation is defined (This is hardcoded)
T_rec_arr = np.logspace(5, 7, N_points)
A_rec_arr = np.logspace(-1, 3, N_points)

# plot interpolation between two scenarios varying recoupling strength
# idx = 65
# idx_add = idx + 8

# print(
#     "Showing interpolation between points with A_rec = %.2e, T_rec = %.2e and A_rec = %.2e, T_rec = %.2e"
#     % (A_rec_arr[idx], T_rec_arr[idx], A_rec_arr[idx_add], T_rec_arr[idx])
# )

# # Varied Recoupling Strength Plot
# p.plot_varied_recoupling(
#     T_rec_arr[0],
#     A_rec_arr[0],
#     T_rec_arr[idx],
#     T_rec_arr[idx],
#     A_rec_arr[idx],
#     A_rec_arr[idx_add],
# )

# idx = 65
# idx_add = idx + 4

# print(
#     "Showing interpolation between points with A_rec = %.2e, T_rec = %.2e and A_rec = %.2e, T_rec = %.2e"
#     % (A_rec_arr[idx], T_rec_arr[idx], A_rec_arr[idx], T_rec_arr[idx_add])
# )

# # Varied Recoupling Temperature Plot
# p.plot_varied_recoupling(
#     T_rec_arr[0],
#     A_rec_arr[0],
#     T_rec_arr[idx],
#     T_rec_arr[idx_add],
#     A_rec_arr[idx],
#     A_rec_arr[idx],
# )

# idx = 68
# idx_add = idx + 4
# idx_sub = idx - 8

# print(
#     "Showing interpolation between points with A_rec = %.2e, T_rec = %.2e and A_rec = %.2e, T_rec = %.2e"
#     % (A_rec_arr[idx], T_rec_arr[idx], A_rec_arr[idx_sub], T_rec_arr[idx_add])
# )

# # Same Peak $\Gamma_{\mathrm{DR-DM, peak}}$
# ## Todo Need to add k cutoff scale to this plot
# p.plot_varied_recoupling(
#     T_rec_arr[0],
#     A_rec_arr[0],
#     T_rec_arr[idx],
#     T_rec_arr[idx_add],
#     A_rec_arr[idx],
#     A_rec_arr[idx_sub],
# )

# # Current and future observations
# best_fits_saved = np.load("./interps/best_fit_sds.npz")
# best_a_dimless_arr = best_fits_saved["best_a_dimless_arr"]
# plot_idx = 77

# p.plot_observations(
#     T_rec_arr[plot_idx],
#     A_rec_arr[plot_idx],
#     best_a_dimless_arr[plot_idx - 50, plot_idx - 50],
#     mwarm="80kev",
# )

# plot_idx = 98

# p.plot_observations(
#     T_rec_arr[plot_idx],
#     A_rec_arr[plot_idx],
#     best_a_dimless_arr[plot_idx - 50, plot_idx - 50],
#     mwarm="30kev",
# )

# Plot delta function recoupling rate
p.plot_delta_recoupling_rate()

# Plot equation 2.21 in paper
# p.plot_delta_recoupling()
