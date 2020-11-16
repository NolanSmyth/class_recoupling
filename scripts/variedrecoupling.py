from classy import Class
import matplotlib.pyplot as plt
from math import pi
import numpy as np

z_pk = 0.0
pk_max = 1.e2
m_idm = 1.0e3;

commonset = {'omega_b':0.022032,
           'omega_cdm':0.12038,
           'h':0.67556,
           'A_s':2.215e-9,
           'n_s':0.9619,
           'tau_reio':0.0925,
           'output':'tCl,pCl,lCl,mPk',
           'lensing':'yes',
           'P_k_max_1/Mpc':pk_max,
           'z_pk':z_pk}

idmset = {'f_idm_dr':1., 
           'xi_idr':0.5, 
           'stat_f_idr': 1.0, #bosonic
           'a_idm_dr':1.e3, 
           'nindex_idm_dr':4., 
           'm_idm':m_idm}

#Define (z, Gamma) points for multiple runs
z1s = [10**7,10**7]
z2s = [10**5.2,10**5.2,]
z3s = [10**5.15,10**5.15]
z4s = [10**5.0,10**5.0]

g1s = [1, 1]
g2s = [10**(-6), 10**(-6)]
g3s = [10**(2), 10**(3)]
g4s = [10**(-4), 10**(-4)]

############################################
#Run for different Gamma/H points
#############################################
                
#
# arrays for output
#
kvec = np.logspace(-4,np.log10(pk_max),500)
legarray = []
twopi = 2.*pi

Pkarr = []

# loop over varying parameter values
#
for i,var in enumerate(g3s):
    #
    print(' * Compute with %s=%e'%('g3',var))
    #
    # deal with colors and legends

    #    
    # call CLASS
    #
    
    z1 = z1s[i]; z2 = z2s[i]; z3 = z3s[i]; z4 = z1s[i]
    g1 = g1s[i]; g2 = g2s[i]; g3 = g3s[i]; g4 = g1s[i]
    
    M = Class()
    M.set(commonset)
    M.set(idmset)
    M.set({'z1':z1,
          'z2':z2,
          'z3':z3,
          'z4':z4,
          'g1':g1,
          'g2':g2,
          'g3':g3,
          'g4':g4,})
    M.compute()
   
    #
    # get P(k) for common k values
    #

    pkM = []
    Mh = M.h()
    for k in kvec:
        pkM.append(M.pk(k*Mh,z_pk)*Mh**3)
    
    #Save results in array so you don't have to re-run
    Pkarr.append(pkM)
    
    #
    # reset CLASS
    #
    M.struct_cleanup()
    M.empty()    
