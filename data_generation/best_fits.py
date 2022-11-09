import numpy as np
import pickle as pickle

# range of k's for which pk interpolations are valid
pk_max = 1e2
kk = np.logspace(-4, np.log10(pk_max), 500)

N_points = 100
# Values over which dd interpolation is defined
T_rec_arr = np.logspace(5, 7, N_points)
A_rec_arr = np.logspace(-1, 3, N_points)
T_rec_arr_test = T_rec_arr[50:]
A_rec_arr_test = A_rec_arr[50:]

N_points_a = 50
# Values over which sd interpolation is defined
a_idm_dr_arr = np.logspace(-5, 5, N_points_a)

pk_sd_interp = pickle.load(open("interps/pks_sd_interp.p", "rb"))
pk_dd_interp = pickle.load(open("interps/pks_dd_interp.p", "rb"))


def best_sd_fit_dimless(T_rec, A_rec, a_idm_dr_arr):
    """
    For a given double decoupling case, finds the best fit a_idm_dr for a single decoupling case using an l2 metric for the dimensionless power spectrum
    """
    l2best = np.inf
    l2best_a = 0
    for (
        a
    ) in (
        a_idm_dr_arr
    ):  # only searching over a_idm_dr values that were defined in the sd interpolation
        l2 = np.sum(
            [
                (
                    pk_dd_interp((T_rec, A_rec, k)) * (k ** 3) / (2 * np.pi ** 2)
                    - pk_sd_interp((a, k)) * (k ** 3) / (2 * np.pi ** 2)
                )
                ** 2
                for k in kk[-200:]
            ]
        )  # Only look at highest k
        if l2 < l2best:
            l2best = l2
            l2best_a = a
    return l2best, l2best_a


def create_best_sd_fit_dimless(Tr_arr, Ar_arr, a_idm_dr_arr):
    # Scan over params in dimless space
    l2_best_dimless_arr = np.zeros((len(Tr_arr), len(Ar_arr)))
    best_a_dimless_arr = np.zeros((len(Tr_arr), len(Ar_arr)))

    for i, T_rec in enumerate(Tr_arr):
        for j, A_rec in enumerate(Ar_arr):
            l2_best_dimless_arr[i, j], best_a_dimless_arr[i, j] = best_sd_fit_dimless(
                T_rec, A_rec, a_idm_dr_arr
            )

    np.savez(
        "./interps/best_fit_sds_2.npz",
        l2_best_dimless_arr=l2_best_dimless_arr,
        best_a_dimless_arr=best_a_dimless_arr,
        T_rec_arr_test=Tr_arr,
        A_rec_arr_test=Ar_arr,
    )


create_best_sd_fit_dimless(T_rec_arr_test, A_rec_arr_test, a_idm_dr_arr)
