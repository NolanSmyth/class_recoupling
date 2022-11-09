import numpy as np

# CLASS variables
pk_max = 1.0e2
# maximum k for Pk
kk = np.logspace(-4, np.log10(pk_max), 500)
# redshift at which Pk is determined
z_pk = 0.0
BM_KS = ["10"]
omega_b = 0.022032
omega0_cdm = 0.12038
stat_f_idr = 0.875  # Fermionic
sigma_fac = 0.01
f_idm_dr = 1.0
xi_idr = 0.3
nindex_idm_dr = 4.0
m_idm = 1.0e3
a_idm_dr = 1.0e0
T_rec = 6.0e5
k = 10  # use one k mode for now
h = 0.67556
A_s = 2.215e-9
n_s = 0.9667
tau_reio = 0.0925
z_pk = 0.0
k_output_values = str(k)

# Other variables
spline_pars = {"k": 3, "s": 0.0}
