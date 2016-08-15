# -*- coding: utf-8 -*-
"""
Created on Sat Jan 30 18:19:12 2016

@author: admin-bellei
"""
import pickle
import numpy as np
import matplotlib.pylab as plt
import sys
from TDMAsolver import TDMAsolver
from pandas import DataFrame

import csv

qe_C = 1.6e-19

nn = 316
LDIAG = list( csv.reader(open('LDIAG.csv')) )
LDIAG = np.asarray(LDIAG).astype(float)
LDIAG = LDIAG[nn:]
DIAG = list( csv.reader(open('DIAG.csv')) )
DIAG = np.asarray(DIAG).astype(float)
DIAG = DIAG[nn:]
UDIAG = list( csv.reader(open('UDIAG.csv')) )
UDIAG = np.asarray(UDIAG).astype(float)
UDIAG = UDIAG[nn:]
BB = list( csv.reader(open('BB.csv')) )
BB = np.asarray(BB).astype(float)
BB = BB[nn:]
T = list( csv.reader(open('T.csv')) )
T = np.asarray(T).astype(float).ravel()
try:
    Tnew_code = list( csv.reader(open('Tnew.csv')) )
    Tnew_code = np.asarray(Tnew_code).astype(float).ravel()
except:
    pass
try:
    T_before = list( csv.reader(open('T_before.csv')) )
    T_before = np.asarray(T_before).astype(float).ravel()
except:
    pass
mu_ion = list( csv.reader(open('mu_ion.csv')) )
mu_ion = np.asarray(mu_ion).astype(float).ravel()
try:
    muold = list( csv.reader(open('muold.csv')) )
    muold = np.asarray(muold).astype(float).ravel()
except:
    pass
#U1Dnew = list( csv.reader(open('U1Dnew.csv')) )
#U1Dnew = np.asarray(U1Dnew).astype(float).ravel()
#U1D = list( csv.reader(open('U1D.csv')) )
#U1D = np.asarray(U1D).astype(float).ravel()
rho = list( csv.reader(open('rho.csv')) )
rho = np.asarray(rho).astype(float).ravel()
u = list( csv.reader(open('u.csv')) )
u = np.asarray(u).astype(float).ravel()
u_old = list( csv.reader(open('u_old.csv')) )
u_old = np.asarray(u_old).astype(float).ravel()

plt.figure(figsize=(8,8))
plt.subplot(2,2,1)
plt.plot(r,u)
plt.plot(r,u_old)
plt.title("u,uold")
plt.subplot(2,2,2)
plt.plot(rho)
plt.title("rho")
plt.subplot(2,2,3)
plt.plot(mu_ion)
try:
    plt.plot(muold)
except:
    pass
plt.title("mu_ion,muold")
plt.subplot(2,2,4)
plt.plot(r,T/qe_C)
plt.plot(r,Tnew_code/qe_C)
plt.title("T, Tnew [eV]")


plt.figure()
uold = np.copy(u)
unew = np.copy(u)
X2_1 = TDMAsolver(LDIAG,DIAG,UDIAG,BB)
unew[nn+1:-1] = u[nn+1:-1] - X2_1.ravel()
plt.plot(unew)
plt.plot(uold)


Told = np.copy(Tnew_code)
Tnew = np.copy(Tnew_code)
X2_1 = TDMAsolver(LDIAG,DIAG,UDIAG,BB)

Tnew[nn+1:-1] = Tnew_code[nn+1:-1] - X2_1.ravel()

plt.figure()
plt.plot(Tnew/qe_C)
plt.plot(Told/qe_C)
plt.plot(Tnew_code/qe_C)

res2 = np.max(abs(Tnew-Tnew_code)/Tnew_code)

sys.exit(0)


J = np.zeros((len(DIAG),len(DIAG)))
for i in range(len(DIAG)):
    J[i,i] = DIAG[i]
    if i>0:
        J[i,i-1] = LDIAG[i-1]
    if i<len(DIAG)-1:
        J[i,i+1] = UDIAG[i]
    
Told = np.copy(T)
Tnew = np.copy(T)
Jinverse = np.linalg.inv
X2_2 = np.dot(Jinverse,BB)

Tnew[nn+1:-1] = T[nn+1:-1] - X2_2.ravel()
residue2 = abs(Tnew[nn+1:-1]-Told[nn+1:-1])/Told[nn+1:-1]

plt.figure()
plt.semilogy(residue1)
plt.semilogy(residue2)


sys.exit(0)


T = list( csv.reader(open('T.csv')) )
T = np.asarray(T).astype(float).ravel()
mu_ion = list( csv.reader(open('mu_ion.csv')) )
mu_ion = np.asarray(mu_ion).astype(float).ravel()
r = list( csv.reader(open('r.csv')) )
r = np.asarray(r).astype(float).ravel()
rho = list( csv.reader(open('rho.csv')) )
rho = np.asarray(rho).astype(float).ravel()
u = list( csv.reader(open('u.csv')) )
u = np.asarray(u).astype(float).ravel()
u_old = list( csv.reader(open('u_old.csv')) )
u_old = np.asarray(u_old).astype(float).ravel()
U1D = list( csv.reader(open('U1D.csv')) )
U1D = np.asarray(U1D).astype(float).ravel()
L_ab = list( csv.reader(open('L_ab.csv')) )
L_ab = np.asarray(L_ab).astype(float).ravel()

rr = list( csv.reader(open('rr.csv')) )
rr = np.asarray(rr).astype(float).ravel()



#f = open("variables", "r")
#T = np.asarray(pickle.load(f))
#mi = pickle.load(f)
#mu_ion = np.asarray(pickle.load(f))
#r = np.asarray(pickle.load(f))
#rho = np.asarray(pickle.load(f))
#rr = np.asarray(pickle.load(f))
#u = np.asarray(pickle.load(f))
#u_old = np.asarray(pickle.load(f))
#U1D = np.asarray(pickle.load(f))
#f.close()

g = 5./3
qe_C = 1.6e-19
Ai = 6.5
C = 5.016e22
mi = 1.0871874140331813E-026 

mu_ip12 = 0.5 * ( mu_ion[2:] + mu_ion[1:-1] )
#mu_ip12 = np.insert(mu_ip12,0,0.)
#mu_ip12 = np.insert(mu_ip12,len(mu_ip12),0.)

mu_im12 = 0.5 * ( mu_ion[1:-1] + mu_ion[0:-2] )
#mu_im12 = np.insert(mu_im12,0,0.)
#mu_im12 = np.insert(mu_im12,len(mu_im12),0.)

u_ip12  = 0.5 * ( u[2:] + u[1:-1] )
#u_ip12 = np.insert(u_ip12,0,0.)
#u_ip12 = np.insert(u_ip12,len(u_ip12),0.)

u_im12  = 0.5 * ( u[1:-1] + u[0:-2] )
#u_im12 = np.insert(u_im12,0,0.)
#u_im12 = np.insert(u_im12,len(u_im12),0.)

u_oldip12  = 0.5 * ( u_old[2:] + u_old[1:-1] )
#u_oldip12 = np.insert(u_oldip12,0,0.)
#u_oldip12 = np.insert(u_oldip12,len(u_oldip12),0.)

u_oldim12  = 0.5 * ( u_old[1:-1] + u_old[0:-2] )
#u_oldim12 = np.insert(u_oldim12,0,0.)
#u_oldim12 = np.insert(u_oldim12,len(u_oldim12),0.)

r_ip12 = 0.5 * ( r[2:] + r[1:-1] )
#r_ip12 = np.insert(r_ip12,0,0.)
#r_ip12 = np.insert(r_ip12,len(r_ip12),0.)

r_im12 = 0.5 * ( r[1:-1] + r[0:-2] )
#r_im12 = np.insert(r_im12,0,0.)
#r_im12 = np.insert(r_im12,len(r_im12),0.)


Tnew = np.copy(T)

    
#diff_diag = np.zeros(len(F))
#diff_ldiag = np.zeros(len(F)-1)
#diff_udiag = np.zeros(len(F)-1)
#
#for i in xrange(0,len(F)):
#    diff_diag[i] = abs(dFdt[i,i]-DIAG[i])/abs(dFdt[i,i])
#    if (i>0):    
#        diff_ldiag[i-1] = abs(dFdt[i,i-1]-LDIAG[i-1])/abs(dFdt[i,i-1])
#    if (i<len(F)-1):
#        diff_udiag[i] = abs(dFdt[i,i+1]-UDIAG[i])/abs(dFdt[i,i+1])
#
#plt.figure()
#plt.semilogy(diff_diag)
#plt.semilogy(diff_ldiag)
#plt.semilogy(diff_udiag)



