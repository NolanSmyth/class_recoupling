
# coding: utf-8

# In[ ]:

# import classy module
from classy import Class
import matplotlib.pyplot as plt
from math import pi
import numpy as np


z_pk = 0.0
pk_max = 110.

# create instance of the class "Class"
LambdaCDM = Class()
# pass input parameters
LambdaCDM.set({'omega_b':0.022032,'omega_cdm':0.12038,'h':0.67556,'A_s':2.215e-9,'n_s':0.9619,'tau_reio':0.0925})
LambdaCDM.set({'output':'tCl,pCl,lCl,mPk','lensing':'yes','P_k_max_1/Mpc':pk_max})
LambdaCDM.set({'z_pk':z_pk})
# run class
LambdaCDM.compute()


# # get all C_l output
# cls = LambdaCDM.lensed_cl(2500)
# # To check the format of cls
# # cls.viewkeys()


# ll = cls['ell'][2:]
# clTT = cls['tt'][2:]
# clEE = cls['ee'][2:]
# clPP = cls['pp'][2:]



# # get P(k) at redhsift z=0

# kk = np.logspace(-4,np.log10(3),1000) # k in h/Mpc
kk = np.logspace(np.log10(2),np.log10(pk_max),500) # k in h/Mpc
Pk = [] # P(k) in (Mpc/h)**3
h = LambdaCDM.h() # get reduced Hubble for conversions to 1/Mpc
for k in kk:
    Pk.append(LambdaCDM.pk(k*h,z_pk)*h**3) # function .pk(k,z)


# optional: clear content of LambdaCDM (to reuse it for another model)
LambdaCDM.struct_cleanup()
# optional: reset parameters to default
LambdaCDM.empty()

# create instance of the class "Class"
mycos = Class()
# pass input parameters
mycos.set({'omega_b':0.022032,'omega_cdm':0.12038,'h':0.67556,'A_s':2.215e-9,'n_s':0.9619,'tau_reio':0.0925})
mycos.set({'output':'tCl,pCl,lCl,mPk','lensing':'yes','P_k_max_1/Mpc':pk_max})

mycos.set({'f_idm_dr':1., 'xi_idr':0.5, 'a_idm_dr':1.e3, 'nindex_idm_dr':4., 'm_idm':1.0e11})
mycos.set({'z_pk':z_pk})
mycos.set({'z_scale':1.e-8})

# run class
mycos.compute()


# # get all C_l output
# mycls = mycos.lensed_cl(2500)

# myll = mycls['ell'][2:]
# myclTT = mycls['tt'][2:]
# myclEE = mycls['ee'][2:]
# myclPP = mycls['pp'][2:]

mykk = np.logspace(np.log10(2),np.log10(pk_max),500) # k in h/Mpc
myPk = [] # P(k) in (Mpc/h)**3
myh = mycos.h() # get reduced Hubble for conversions to 1/Mpc
for k in mykk:
    myPk.append(mycos.pk(k*myh,z_pk)*myh**3) # function .pk(k,z)

mycos.struct_cleanup()
mycos.empty()

# mycos = Class()
# # pass input parameters
# mycos.set({'omega_b':0.022032,'omega_cdm':0.12038,'h':0.67556,'A_s':2.215e-9,'n_s':0.9619,'tau_reio':0.0925})
# mycos.set({'output':'tCl,pCl,lCl,mPk','lensing':'yes','P_k_max_1/Mpc':pk_max})

# mycos.set({'f_idm_dr':1., 'xi_idr':0.5, 'a_idm_dr':1e2, 'nindex_idm_dr':4})
# mycos.set({'z_pk':z_pk})
# # run class
# mycos.compute()


# # get all C_l output
# mycls = mycos.lensed_cl(2500)

# myll = mycls['ell'][2:]
# myclTT = mycls['tt'][2:]
# myclEE = mycls['ee'][2:]
# myclPP = mycls['pp'][2:]

# myPk2 = [] # P(k) in (Mpc/h)**3
# myh = mycos.h() # get reduced Hubble for conversions to 1/Mpc
# for k in mykk:
#     myPk2.append(mycos.pk(k*myh,z_pk)*myh**3) # function .pk(k,z)

# plot Matter Power Spectrum
# plt.figure(2)
# plt.xscale('log'); plt.yscale('log');
# plt.xlabel(r'$k \,\,\,\, [h/\mathrm{Mpc}]$')
# plt.ylabel(r'$P(k) \,\,\,\, [\mathrm{Mpc}/h]^3$')
# plt.title('P(k)')
# plt.xlim(kk[0],kk[-1])
# plt.plot(kk,Pk,'b-', label=r'$\Lambda CDM$')

# plt.savefig('figures/lambda_cdm_all.pdf')

# plt.figure(2)
# plt.xscale('log');plt.yscale('log');plt.xlim(mykk[0],mykk[-1])
# plt.xlabel(r'$k \,\,\,\, [h/\mathrm{Mpc}]$')
# plt.ylabel(r'$P(k) \,\,\,\, [\mathrm{Mpc}/h]^3$')
# plt.plot(mykk,myPk,'r-', label='idm_dr')

# CMBratio = [(myclTT[i]*myll[i]*(myll[i]+1)/2./pi)/(clTT[i]*ll[i]*(ll[i]+1)/2./pi) for i in range(len(ll))]


# plt.figure(1)
# plt.xscale('log');plt.yscale('linear');plt.xlim(2,2500)
# plt.xlabel(r'$\ell$')
# plt.ylabel(r'$[\ell(\ell+1)/2\pi]  C_\ell^\mathrm{TT}$')
# plt.plot(ll,CMBratio,'k-')
# plt.plot(ll,clTT*ll*(ll+1)/2./pi,'k-')
# plt.plot(myll,myclTT*myll*(myll+1)/2./pi,'r-')

# plt.savefig('figures/cmb_ang_ratio')


# plt.legend()

pkratio = [myPk[i]/Pk[i] for i in range(len(Pk))]
# pkratio2 = [myPk2[i]/Pk[i] for i in range(len(Pk))]

plt.figure(2)
plt.xscale('log');plt.yscale('linear');plt.xlim(mykk[0],mykk[-1])
plt.ylim(0,1)
plt.xlabel(r'$k \,\,\,\, [h/\mathrm{Mpc}]$')
plt.ylabel(r'$P(k)/P(k)_{\Lambda CDM} \,\,\,\,$')
plt.title('z = 0, a scaled')
plt.plot(mykk,pkratio,'k-', label=r'$\xi = 0.5, a = 1e3$')
# plt.plot(mykk,pkratio2,'m-', label=r'$\xi = 0.5, a = 1e2$')



plt.savefig('figures/pk_ratio_z_scale_1e-8_test.pdf')
