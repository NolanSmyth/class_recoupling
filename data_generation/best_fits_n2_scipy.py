import numpy as np
import pickle as pickle
from scipy.optimize import minimize


# range of k's for which pk interpolations are valid
pk_max = 1e2
kk = np.logspace(-4, np.log10(pk_max), 500)

N_points = 100
# Values over which dd interpolation is defined

T_rec_strong_arr = np.logspace(4, 7, N_points)
A_rec_strong_arr = np.logspace(-1, 7, N_points)

T_rec_arr_test = T_rec_strong_arr
A_rec_arr_test = A_rec_strong_arr

N_points_a = 100  # Number of grid points
a_idm_dr_arr = np.logspace(-1, 9, N_points_a)

pk_sd_interp = pickle.load(open("interps/pks_sd_interp_n2.p", "rb"))
pk_dd_interp = pickle.load(open("interps/pks_strong_interp_n2.p", "rb"))


def objective_function(a, T_rec, A_rec, kk):
    return np.sum(
        [
            (
                pk_dd_interp((T_rec, A_rec, k)) * (k**3) / (2 * np.pi**2)
                - pk_sd_interp((a, k)) * (k**3) / (2 * np.pi**2)
            )
            ** 2
            for k in kk[-200:]
        ]
    )


def best_sd_fit_dimless(T_rec, A_rec, a_idm_dr_arr):
    """
    For a given double decoupling case, finds the best fit a_idm_dr for a single decoupling case using an l2 metric for the dimensionless power spectrum
    """
    # initial guess for a
    a0 = a_idm_dr_arr[N_points_a // 2]

    result = minimize(
        objective_function,
        a0,
        args=(T_rec, A_rec, kk),
        bounds=[(a_idm_dr_arr.min(), a_idm_dr_arr.max())],
    )

    if result.success:
        return result.x.item()
    else:
        return -1


def create_best_sd_fit_dimless(Tr_arr, Ar_arr, a_idm_dr_arr):
    # Scan over params
    # l2_best_dimless_arr = np.zeros((len(Tr_arr), len(Ar_arr)))
    best_a_dimless_arr = np.zeros((len(Tr_arr), len(Ar_arr)))

    for i, T_rec in enumerate(Tr_arr):
        for j, A_rec in enumerate(Ar_arr):
            best_a_dimless_arr[i, j] = best_sd_fit_dimless(T_rec, A_rec, a_idm_dr_arr)

    np.savez(
        "./interps/best_fit_sds_n2.npz",
        # l2_best_dimless_arr=l2_best_dimless_arr,
        best_a_dimless_arr=best_a_dimless_arr,
        T_rec_arr_test=Tr_arr,
        A_rec_arr_test=Ar_arr,
    )


create_best_sd_fit_dimless(T_rec_arr_test, A_rec_arr_test, a_idm_dr_arr)
