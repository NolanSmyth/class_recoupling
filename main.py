import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import UnivariateSpline
import h5py
from scipy.integrate import quadrature
import pickle
import plots as p

N_points = 100
# Values over which dd interpolation is defined (This is hardcoded)
T_rec_arr = np.logspace(5, 7, N_points)
A_rec_arr = np.logspace(-1, 3, N_points)

# # # plot interpolation between two scenarios varying recoupling strength
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

# print("{:.2e}, {:.2e}, {:.2e}".format( T_rec_arr[plot_idx],
#     A_rec_arr[plot_idx],
#     best_a_dimless_arr[plot_idx - 50, plot_idx - 50]))

# p.plot_observations(
#     T_rec_arr[plot_idx],
#     A_rec_arr[plot_idx],
#     best_a_dimless_arr[plot_idx - 50, plot_idx - 50],
#     mwarm="80kev",
# )

# plot_idx = 81

# print("{:.2e}, {:.2e}, {:.2e}".format( T_rec_arr[plot_idx],
#     A_rec_arr[plot_idx],
#     best_a_dimless_arr[plot_idx - 50, plot_idx - 50]))

# p.plot_observations(
#     T_rec_arr[plot_idx],
#     A_rec_arr[plot_idx],
#     best_a_dimless_arr[plot_idx - 50, plot_idx - 50],
#     # mwarm="30kev",
#     mwarm="50kev",
# )

# Plot delta function recoupling rate
# p.plot_delta_recoupling_rate()

# Plot equation 2.21 in paper
# p.plot_delta_effect()

# p.plot_delta_power_spectrum()

# p.plot_delta_power_spectrum_dimless()

# p.plot_delta_power_spectra_both()

# p.plot_delta_effect_both()


idx = 65
Trs = [
    # T_rec_arr[idx - 10],
    T_rec_arr[idx],
    T_rec_arr[idx + 15],
    T_rec_arr[idx + 30],
    # T_rec_arr[99],
]
Ars = [
    A_rec_arr[idx - 10],
    A_rec_arr[idx],
    A_rec_arr[idx + 10],
    A_rec_arr[idx + 20],
    A_rec_arr[idx + 30],
]

# idx = 65

# Trs = [T_rec_arr[idx], T_rec_arr[idx + 15], T_rec_arr[idx + 30]]
# Ars = [A_rec_arr[idx], A_rec_arr[idx - 4], A_rec_arr[idx - 8]]

# p.plot_varied_recoupling_grid(T_rec_arr[0], A_rec_arr[0], Trs, Ars)
p.plot_varied_recoupling_grid_collapsed(T_rec_arr[0], A_rec_arr[0], Trs, Ars)

