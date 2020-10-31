from classy import Class
import matplotlib.pyplot as plt
from math import pi
import numpy as np

z_pk = 0.0
pk_max = 1.e2
m_idm = 1.0e3;

############################################
#
# Varying parameter (others fixed to default)
#
var_name = 'z_recouple'
var_array = [1.0e5, 2.0e5, 3.0e5, 4.0e5, 5.0e5]
var_num = len(var_array)
var_legend = 'z_recouple'
var_figname = 'z_recouple'
#
#############################################
#
# Fixed settings
#
common_settings = {'output':'mPk',
#                    'output':'tCl,pCl,lCl,mPk', #Try not calculating these to speed up computation
#                    'lensing':'yes',
                   # LambdaCDM parameters
                   'h':0.67556,
                   'omega_b':0.022032,
                   'omega_cdm':0.12038,
                   'A_s':2.215e-9,
                   'n_s':0.9619,
                   'tau_reio':0.0925,
                   # other output and precision parameters
                   'P_k_max_1/Mpc':pk_max,
                   # idm_dr parameters
                   'z_pk':z_pk,
                   'f_idm_dr':1.,
                   'xi_idr':0.5, 
                   'a_idm_dr':1.e3,
                   'nindex_idm_dr':4.,
                   'm_idm':1.0e3,
                   'zd1':1.0e7,
                   'zd2':1.0e5}

                
#
# arrays for output
#
kvec = np.logspace(-4,np.log10(pk_max),500)
legarray = []
twopi = 2.*pi

Pkarr = []
# llarr = []
# clTTarr = []

#create figures
fig_Pk, ax_Pk = plt.subplots()

# loop over varying parameter values
#
for i,var in enumerate(var_array):
    #
    print(' * Compute with %s=%e'%(var_name,var))
    #
    # deal with colors and legends
    #
#     if i == 0:
#         var_color = 'k'
#         var_alpha = 1.
#         legarray.append(r'ref. $\Lambda CDM$')
#     else:
    var_color = 'r'
    var_alpha = 1.*(i+1)/(var_num)
    if i == var_num-1:
        legarray.append(var_legend)  
    #    
    # call CLASS
    #
    M = Class()
    M.set(common_settings)
    M.set({var_name:var})
    M.compute()
    #
    # get Cls
    #
    
#     clM = M.lensed_cl(2500)
#     ll = clM['ell'][2:]
#     clTT = clM['tt'][2:]
    
    '''
    clEE = clM['ee'][2:]
    clPP = clM['pp'][2:]
    '''

    #
    # get P(k) for common k values
    #

    pkM = []
    Mh = M.h()
    for k in kvec:
        pkM.append(M.pk(k*Mh,z_pk)*Mh**3)
    #    
    # plot P(k)
    #
    ax_Pk.plot(kvec,np.array(pkM),color=var_color,alpha=var_alpha,linestyle='-',label=f'z={var:.2E}')
    #
    # plot C_l^TT
    #
    '''
    ax_TT.semilogx(ll,clTT*ll*(ll+1)/twopi,color=var_color,alpha=var_alpha,linestyle='-')
    #
    # plot Cl EE 
    #
    ax_EE.loglog(ll,clEE*ll*(ll+1)/twopi,color=var_color,alpha=var_alpha,linestyle='-')
    #
    # plot Cl phiphi
    #
    ax_PP.loglog(ll,clPP*ll*(ll+1)*ll*(ll+1)/twopi,color=var_color,alpha=var_alpha,linestyle='-')
    '''
    
    #Save results in array so you don't have to re-run
    Pkarr.append(pkM)
#     llarr.append(ll)
#     clTTarr.append(clTT)
    
    #
    # reset CLASS
    #
    M.struct_cleanup()
    M.empty()    

ax_Pk.set_xlim([1.e-4,pk_max])
ax_Pk.set_xlabel(r'$k \,\,\,\, [h/\mathrm{Mpc}]$')
ax_Pk.set_ylabel(r'$P(k) \,\,\,\, [\mathrm{Mpc}/h]^3$')
ax_Pk.loglog()
ax_Pk.legend()
fig_Pk.tight_layout()
# fig_Pk.savefig('../figures/pk_varied_decoupling.pdf')