import matplotlib.pyplot as plt
import numpy as np
import pickle
import multiprocessing as mp
import time

def best_sd_fit(T_rec, A_rec, a_idm_dr_arr):
    '''
    For a given double decoupling case, finds the best fit a_idm_dr for a single decoupling case using an l2 metric
    #currently ~4.5 secs. Can improve?
    '''
    l2best = np.inf
    l2best_a = 0
    for a in a_idm_dr_arr: #only searching over a_idm_dr values that were defined in the sd interpolation
        l2 = np.sum([(pk_dd_interp((T_rec, A_rec, k)) - pk_sd_interp((a, k)))**2 for k in kk[-200:]]) #Only look at highest k
        if l2 < l2best:
            l2best = l2
            l2best_a = a
    return l2best, l2best_a

def best_sd_fit_parallel(index):
    T_rec = T_rec_arr_test[index//len(A_rec_arr_test)]
    A_rec = A_rec_arr_test[index%len(A_rec_arr_test)]
    best_sd_fit(T_rec, A_rec, a_idm_dr_arr)

def run_best_sd_fit(operation, index, pool):
    # pool.map(operation, (T_rec_arr_test[T_rec_idx], A_rec_arr_test[A_rec_idx], a_idm_dr_arr))
    pool.map(operation, index)

#Load pk interpolations from pickle file
pk_sd_interp = pickle.load(open('interps/pks_sd_interp.p','rb'))
pk_dd_interp = pickle.load(open('interps/pks_dd_interp.p','rb'))

pk_max = 1e2
kk = np.logspace(-4, np.log10(pk_max), 500)

N_points = 100
#Values over which dd interpolation is defined
T_rec_arr = np.logspace(5, 7, N_points)
A_rec_arr = np.logspace(-1, 3, N_points)

N_points_a = 50
#Values over which sd interpolation is defined
a_idm_dr_arr = np.logspace(-5, 5, N_points_a)

#Not parallelized
T_rec_arr_test = T_rec_arr[50:52]
A_rec_arr_test = A_rec_arr[50:52]

#Parallelized
if __name__ == '__main__':

    ti_un = time.perf_counter()
    for i, T_rec in enumerate(T_rec_arr_test):
        for j, A_rec in enumerate(A_rec_arr_test):
            best_sd_fit(T_rec, A_rec, a_idm_dr_arr)
    tf_un = time.perf_counter()
    print('Time for unparallelized:', tf_un - ti_un)

    ti_par = time.perf_counter()
    with mp.Pool(mp.cpu_count()) as p:
        # run_best_sd_fit(best_sd_fit_parallel, range(len(T_rec_arr_test)), range(len(A_rec_arr_test)), p)
        run_best_sd_fit(best_sd_fit_parallel, range(len(T_rec_arr_test)*len(A_rec_arr_test)), p)
    tf_par = time.perf_counter()
    print('Time for parallelized:', tf_par - ti_par)    


