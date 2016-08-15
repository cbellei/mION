# -*- coding: utf-8 -*-
"""
Created on Sat Feb 27 11:43:09 2016

@author: admin-bellei
"""
import pickle
import numpy as np
import matplotlib.pylab as plt
import sys
from TDMAsolver import TDMAsolver
from pandas import DataFrame

import csv
nn = 316
rr = list( csv.reader(open('rr.csv')) )
rr = np.float(np.asarray(rr).astype(float))
mi = list( csv.reader(open('mi.csv')) )
mi = np.float(np.asarray(mi).astype(float))
LDIAG = list( csv.reader(open('LDIAG.csv')) )
LDIAG = np.asarray(LDIAG).astype(float)
DIAG = list( csv.reader(open('DIAG.csv')) )
DIAG = np.asarray(DIAG).astype(float)
UDIAG = list( csv.reader(open('UDIAG.csv')) )
UDIAG = np.asarray(UDIAG).astype(float)
BB = list( csv.reader(open('BB.csv')) )
BB = np.asarray(BB).astype(float)
try:
    BB_dopo = list( csv.reader(open('BB2_dopo.csv')) )
    BB_dopo = np.asarray(BB_dopo).astype(float)
except:
    pass
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
U1D = list( csv.reader(open('U1D.csv')) )
U1D = np.asarray(U1D).astype(float).ravel()
rho = list( csv.reader(open('rho.csv')) )
rho = np.asarray(rho).astype(float).ravel()
u = list( csv.reader(open('u.csv')) )
u = np.asarray(u).astype(float).ravel()
u_old = list( csv.reader(open('u_old.csv')) )
u_old = np.asarray(u_old).astype(float).ravel()
L_ab = list( csv.reader(open('L_ab.csv')) )
L_ab = np.asarray(L_ab).astype(float).ravel()
r = list( csv.reader(open('radius.csv')) )
r = np.asarray(r).astype(float).ravel()

g = 5./3
qe_C = 1.6e-19
Ai = 6.5
C = 5.016e22
CC = 5.016e22
nz = len(T)

Tnew = np.copy(T)
U1Dnew = np.copy(U1D)
res2 = 1.e5
count2 = 1
BB2 = np.zeros(nz-2)
LDIAG2 = np.zeros(nz-3)
UDIAG2 = np.zeros(nz-3)
DIAG2  = np.zeros(nz-2)
AA2  = np.zeros((nz-2,nz-2))


for k in np.arange(0,20):
    is_diag_dominant = True
    for i in range(nn-1,nz-3): #Because AA has size (nz-2)x(nz-2)

        mu_ion_ip12 = 0.5 * (  mu_ion[i+1] + mu_ion[i+2] )
        mu_ion_im12 = 0.5 * (  mu_ion[i] + mu_ion[i+1] )
        u_ip12  = 0.5 * ( u[i+1] + u[i+2] )
        u_im12  = 0.5 * ( u[i] + u[i+1] )
        u_oldip12  = 0.5 * ( u_old[i+1] + u_old[i+2] )
        u_oldim12  = 0.5 * ( u_old[i] + u_old[i+1] )
        r_ip12 = 0.5 * ( r[i+1] + r[i+2] )
        r_im12 = 0.5 * ( r[i] + r[i+1] )
        
        BB2[i] = U1Dnew[i+1] - U1D[i+1] - 0.5 * rr/r[i+1]**2 *  \
                    (  \
                        r_ip12**2*mu_ion_ip12*u_ip12*u[i+2]  \
                        - (r_ip12**2*mu_ion_ip12*u_ip12+r_im12**2*mu_ion_im12*u_im12)*u[i+1]  \
                        + r_im12**2*mu_ion_im12*u_im12*u[i] \
                        + r_ip12**2*mu_ion_ip12*u_oldip12*u_old[i+2]   \
                        - (r_ip12**2*mu_ion_ip12*u_oldip12+r_im12**2*mu_ion_im12*u_oldim12)*u_old[i+1] \
                        + r_im12**2*mu_ion_im12*u_oldim12*u_old[i]   \
                        )
#        print U1Dnew[i+1], U1D[i+1]
#        print r_ip12**2, mu_ion_ip12, u[i+2]
#
        #diagonal
        DIAG2[i] = 3./2*rho[i+1]/mi - 0.25*rr*  \
                    (  \
                        r_ip12**2*u_ip12*u[i+2]      \
                        - (r_ip12**2*u_ip12+r_im12**2*u_im12)*u[i+1]  \
                        + r_im12**2*u_im12*u[i] \
                        + r_ip12**2*u_oldip12*u_old[i+2]   \
                        - (r_ip12**2*u_oldip12+r_im12**2*u_oldim12)*u_old[i+1] \
                        + r_im12**2*u_oldim12*u_old[i]   \
                    ) * 3./2 * CC / L_ab[i+1] * np.sqrt(Tnew[i+1])
        AA2[i,i] = DIAG2[i]
        
#        print DIAG2[i]
#        print rho[i+1], rr, r_ip12, u_ip12,u[i+2] 
#        print r[i+1], r[i+2]
#        print CC, L_ab[i+1], Tnew[i+1] ###OK!!!!
#        
#        print "rho = ", rho[i+1]
#        print 3./2*rho[i+1]/mi
#        print 0.25*rr
#        print r_ip12**2*u_ip12*u[i+2]
#        print (r_ip12**2*u_ip12+r_im12**2*u_im12)*u[i+1]
#        print r_im12**2*u_im12*u[i]
#        print r_ip12**2*u_oldip12*u_old[i+2]
#        print (r_ip12**2*u_oldip12+r_im12**2*u_oldim12)*u_old[i+1]
#        print r_im12**2*u_oldim12*u_old[i]
#        print 3./2 * CC / L_ab[i+1] * np.sqrt(Tnew[i+1])
#        sys.exit(0)
                            
        if (i<nz-3):  #take care of upper diagonal       
            UDIAG2[i] = -0.25*rr* \
                (  \
                    r_ip12**2*u_ip12*u[i+2]  \
                    - r_ip12**2*u_ip12*u[i+1]  \
                    + r_ip12**2*u_oldip12*u_old[i+2]   \
                    - r_ip12**2*u_oldip12*u_old[i+1] \
                ) * 3./2 * CC / L_ab[i+2] * np.sqrt(Tnew[i+2])
            AA2[i,i+1] = UDIAG2[i]
     
        if (i>nn-1): #take care of lower diagonal
             LDIAG2[i-1] = -0.25*rr* \
                 (  \
                     - r_im12**2*u_im12*u[i+1]  \
                     + r_im12**2*u_im12*u[i] \
                     - r_im12**2*u_oldim12*u_old[i+1] \
                     + r_im12**2*u_oldim12*u_old[i]   \
                     ) * 3./2 * CC / L_ab[i] * np.sqrt(Tnew[i])
             AA2[i,i-1] = LDIAG2[i-1]
        if(is_diag_dominant):
         	if ( (i>nn-1) and (i<nz-3) and \
			 ( abs(DIAG2[i]) < abs(UDIAG2[i])+abs(LDIAG2[i-1]) ) ):
				is_diag_dominant = False
		elif ( (i==nn-1) and ( abs(DIAG2[i]) < abs( UDIAG2[i])) ):
				is_diag_dominant = False
		elif ( (i==nz-3) and ( abs(DIAG2[i]) < abs( LDIAG2[i-1])) ):
				is_diag_dominant = False

    #sys.exit(0)
    plt.plot(BB2/qe_C)
#    print "prima: ", BB2[nn-1]
    if(is_diag_dominant):
        BB2[nn-1:] = TDMAsolver(LDIAG[nn-1:],DIAG[nn-1:],UDIAG[nn-1:],BB2[nn-1:])
    else:
        Jinverse = np.linalg.inv(AA2[nn-1:,nn-1:])
        BB2[nn-1:] = np.dot(Jinverse,BB2[nn-1:])
#    print "dopo: ", BB2[nn-1]    

    plt.plot(BB2/qe_C)
        
    if (count2>1):
        T_before = np.copy(Tnew)         

    
    Tnew[1:-1] = Tnew[1:-1] - BB2
    #update U1Dnew and viscosity
    TeV = Tnew / qe_C
    p = rho * Tnew / mi
    U1Dnew = p / (g-1) + 0.5 * rho * u**2
    #NRL page 38 - excluding density because it cancels out anyways
    taui =  2.09e7 * TeV**1.5 / L_ab * np.sqrt(Ai)
    mu_ion = 0.96 * ( 1.6e-12 * TeV )  * taui #all cgs
    mu_ion = 1.e-1 * mu_ion # [mu_ion] = [M*L^-1*T^-1]
    #Redefine ion viscosity to include the factor 1./3 in front of derivative
    mu_ion = 1./3 * mu_ion
    
    plt.plot(r,Tnew)
				
    if (count2 > 1):
        res2 = np.max(np.abs(Tnew-T_before)/Tnew)
        #plt.plot(np.abs(1.-T_before/Tnew))
        print "res2 = ", res2        
    count2 += 1

sys.exit(0)
plt.figure()
#T1 = list( csv.reader(open('T1.csv')) )
#T1 = np.asarray(T1).astype(float).ravel()
T2 = list( csv.reader(open('T2.csv')) )
T2 = np.asarray(T2).astype(float).ravel()
T3 = list( csv.reader(open('T3.csv')) )
T3 = np.asarray(T3).astype(float).ravel()

#plt.plot(r,T1)
plt.plot(r,T2)
plt.plot(r,T3)
res2 = np.max(np.abs(T3-T2)/T3)
print res2

