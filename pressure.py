import csv
import os
import matplotlib.pyplot as plt
import sys
import matplotlib
import numpy as np
from scipy.ndimage.filters import gaussian_filter1d

#formatting for figure
matplotlib.rc('font', size=22)
matplotlib.rc('font', family='Arial')

param = list( csv.reader(open('parameters.csv')) )[0]

dtm = np.asarray(param[0]).astype(float) 
nz  = np.asarray(param[1]).astype(int) 
dt_print = np.asarray(param[2]).astype(float)
maxind = np.asarray(param[3]).astype(int)
tm_quiet = np.asarray(param[4]).astype(float)
L = 1.e2 * np.asarray(param[5]).astype(float)  # in cm

x = list( csv.reader(open('r.csv')) )
x = np.asarray(x).astype(float).flatten()
print "done with x"


den1 = list( csv.reader(open('den1.csv')) )
den1 = np.asarray(den1).astype(float).flatten()
print "done with den"
den2 = list( csv.reader(open('den2.csv')) )
den2 = np.asarray(den2).astype(float).flatten()
den3 = list( csv.reader(open('den3.csv')) )
den3 = np.asarray(den3).astype(float).flatten()

temp1 = list( csv.reader(open('temp1.csv')) )
temp1 = np.asarray(temp1).astype(float).flatten()
print "done with temp"
temp2 = list( csv.reader(open('temp2.csv')) )
temp2 = np.asarray(temp2).astype(float).flatten()
print "done with temp"
temp3 = list( csv.reader(open('temp3.csv')) )
temp3 = np.asarray(temp3).astype(float).flatten()
print "done with temp"

xrange = [0.,L];
nfiles = np.int( x.size / nz )

n_total = int(np.ceil(len(x)/float(nz)))

print nz, n_total


print "rmin = ", min(x)



fusion = "DT"
print "fusion = ", fusion

if fusion=="D3He":
	#calculate D3He fusion reactions
	C0 = 10.572
	C1 = 151.16e-16
	C2 = 6.4192e-3
	C3 = -2.029e-3
	C4 = -0.019108e-3
	C5 = 0.13578e-3
	C6 = 0
	C7 = 0
	mfus1 = 2. * 1836 * 9.11e-31 #Deuterium
	mfus2 = 3. * 1836 * 9.11e-31 #3He
elif fusion=="DD":
	#calculate DD ->  p + T
	C0a = 6.2696
	C1a = 3.7212e-16
	C2a = 3.4127e-3
	C3a = 1.9917e-3
	C4a = 0.
	C5a = 0.010506e-3
	C6a = 0
	C7a = 0
	#calculate DD ->  n + 3He
	C0b = 6.2696
	C1b = 3.5741e-16
	C2b = 5.8577e-3
	C3b = 7.6822e-3
	C4b = 0.
	C5b = -0.002964e-3
	C6b = 0
	C7b = 0
	mfus1 = 2. * 1836 * 9.11e-31 #Deuterium
	mfus2 = 2. * 1836 * 9.11e-31 #Deuterium
elif fusion=="DT":
	#calculate DT fusion reactions
	C0 = 6.6610
	C1 = 643.41e-16
	C2 = 15.136e-3
	C3 = 75.189e-3
	C4 = 4.6064e-3
	C5 = 13.500e-3
	C6 = -0.10675e-3
	C7 = -0.01366e-3
	mfus1 = 2. * 1836 * 9.11e-31 #Deuterium
	mfus2 = 3. * 1836 * 9.11e-31 #3He

#mp_g = 1*1836*9.11e-28;
# m_Ti = 47.87 * mp_g #mass of Titanium
#mfus1 = 2. * 1836 * 9.11e-31 #Deuterium
#mfus2 = mfus1 #Deuterium


time = np.array([])
R = np.array([])

for n in range(1,n_total):

	print n
	
	time = np.append(time,dt_print * (n-1))

	start = (n - 1) * nz; finish = start + nz
	
	mu = x[start:finish] #here in m
	n1 = den1[start:finish]
	n2 = den2[start:finish]
	T1 = temp1[start:finish]
	T2 = temp2[start:finish]
	T3 = temp3[start:finish]
	
	mu = 1.e2 * mu #in cm
	mu = np.flipud(mu)
		
	nfus1 = 0.5*(np.flipud(n1)+np.flipud(n2)) #D ions
	nfus2 = nfus1 #T ions
	Tfus1 = 0.5*(np.flipud(T1)+np.flipud(T2)) #D ions
	Tfus2 = Tfus1 #T ions
		
	Teff = 1.e-3 * ( mfus2 * Tfus1 + mfus1 * Tfus2 ) / ( mfus1 + mfus2 )  #in keV
	 
	if fusion=="D3He" or fusion=="DT":
		zeta = 1. - (C2*Teff+C4*Teff**2+C6*Teff**3) / ( 1 + C3*Teff + C5*Teff**2 + C7*Teff**3 )
		xi = C0 / Teff**(1./3)
		sigma_v = C1 * zeta**(-5./6) * xi**2 * np.exp(-3*zeta**(1./3)*xi)
		sigma_v[Teff<1.] = 0.
		sigma_v[mu>0.05] = 0.
	elif fusion=="DD":
		#reaction #1
		zeta = 1. - (C2a*Teff+C4a*Teff**2+C6a*Teff**3) / ( 1 + C3a*Teff + C5a*Teff**2 + C7a*Teff**3 )
		xi = C0a / Teff**(1./3)
		sigma_v_A = C1a * zeta**(-5./6) * xi**2 * np.exp(-3*zeta**(1./3)*xi)

		#reaction #2
		zeta = 1. - (C2b*Teff+C4b*Teff**2+C6b*Teff**3) / ( 1 + C3b*Teff + C5b*Teff**2 + C7b*Teff**3 )
		xi = C0b / Teff**(1./3)
		sigma_v_B = C1b * zeta**(-5./6) * xi**2 * np.exp(-3*zeta**(1./3)*xi)
		
		sigma_v = sigma_v_A + sigma_v_B

	if fusion == "D3He" or fusion=="DT":
		r = np.trapz( 4 * np.pi * mu**2 * nfus1 * nfus2 * sigma_v, mu)
	elif fusion == "DD":
		r = 0.5 * np.trapz( 4 * np.pi * mu**2 * nfus1 * nfus2 * sigma_v, mu)		

 	R = np.append(R,r)

tmax = time[np.where(R==np.max(R))]
print R

n = int(np.round(tmax / dt_print) + 1)
print "n = ", n, tmax/1.e-9
start = (n - 1) * nz; finish = start + nz

mu = x[start:finish]
n1 = den1[start:finish]
n2 = den2[start:finish]
ne = den3[start:finish]
T1 = temp1[start:finish]
T2 = temp2[start:finish]
T3 = temp3[start:finish]

mu = 1.e2 * mu #in cm
mu = np.flipud(mu)	
nfus1 = 0.5*(np.flipud(n1)+np.flipud(n2)) #D ions
nfus2 = nfus1 #T ions
Tfus1 = 0.5*(np.flipud(T1)+np.flipud(T2)) #D ions
Tfus2 = Tfus1 #T ions
	
Teff = 1.e-3 * ( mfus2 * Tfus1 + mfus1 * Tfus2 ) / ( mfus1 + mfus2 )  #in keV 
if fusion=="D3He" or fusion=="DT":
	zeta = 1. - (C2*Teff+C4*Teff**2+C6*Teff**3) / ( 1 + C3*Teff + C5*Teff**2 + C7*Teff**3 )
	xi = C0 / Teff**(1./3)
	sigma_v = C1 * zeta**(-5./6) * xi**2 * np.exp(-3*zeta**(1./3)*xi)
	sigma_v[Teff<1.] = 0.
	sigma_v[mu>0.05] = 0.
elif fusion=="DD":
	#reaction #1
	zeta = 1. - (C2a*Teff+C4a*Teff**2+C6a*Teff**3) / ( 1 + C3a*Teff + C5a*Teff**2 + C7a*Teff**3 )
	xi = C0a / Teff**(1./3)
	sigma_v_A = C1a * zeta**(-5./6) * xi**2 * np.exp(-3*zeta**(1./3)*xi)

	#reaction #2
	zeta = 1. - (C2b*Teff+C4b*Teff**2+C6b*Teff**3) / ( 1 + C3b*Teff + C5b*Teff**2 + C7b*Teff**3 )
	xi = C0b / Teff**(1./3)
	sigma_v_B = C1b * zeta**(-5./6) * xi**2 * np.exp(-3*zeta**(1./3)*xi)
	
	sigma_v = sigma_v_A + sigma_v_B

if fusion == "D3He" or fusion=="DT":
	r = np.trapz( 4 * np.pi * mu**2 * nfus1 * nfus2 * sigma_v, mu)
elif fusion == "DD":
	r = 0.5 * np.trapz( 4 * np.pi * mu**2 * nfus1 * nfus2 * sigma_v, mu)		

j_dndrmax = np.where(np.diff(nfus1)==np.max(np.diff(nfus1))) 
p_total = 1.6e-13 * (n1*T1 + n2*T2 + ne*T3)
p_total = p_total / 1.e14 #Gbar
p_total = np.flipud(p_total)
rhoR = 2.5*1836*9.11e-28*(nfus1+nfus2)*mu
radius = 1.e-2*mu[mu<=mu[j_dndrmax]]
pressure = p_total[mu<=mu[j_dndrmax]]
temp = Teff[mu<=mu[j_dndrmax]]
rhoR = rhoR[mu<=mu[j_dndrmax]]
print "rmax = ", radius[-1]
Rhs = radius[-1]

p_average = np.trapz(4*np.pi*radius**2*pressure,radius)
p_average = p_average / (4./3*np.pi*radius[-1]**3.)

T_average = np.trapz(4*np.pi*radius**2*temp,radius)
T_average = T_average / (4./3*np.pi*radius[-1]**3.)

rhoR_average = np.trapz(4*np.pi*radius**2*rhoR,radius)
rhoR_average = rhoR_average / (4./3*np.pi*radius[-1]**3.)

print "p_average = ", p_average
print "T_average = ", T_average
print "rhoR_average = ", rhoR_average

plt.figure()
plt.subplot(2,2,1)
plt.plot(1.e4*mu,nfus1,mu,nfus2)
plt.plot(1.e4*mu[j_dndrmax],nfus1[j_dndrmax],'o',markersize=8)
plt.xlim([0.,100.])
plt.ylabel("nfus [cm-3]")
plt.subplot(2,2,2)
plt.plot(1.e4*mu,Teff)
plt.xlim([0.,100.])
plt.ylim([0.,20.])
plt.ylabel("Temp [keV]")
plt.subplot(2,2,3)
plt.plot(1.e4*mu,p_total)
plt.xlim([0.,100.])
plt.ylabel("p_total [Gbar]")
plt.subplot(2,2,4)
plt.plot(1.e4*mu,4 * np.pi * mu**2 * nfus1 * nfus2 * sigma_v)
plt.xlim([0.,100.])
plt.ylabel("R_{fusion}")
plt.show()

np.savez("p_average.npz",p_average,T_average,rhoR_average,mu,nfus1+nfus2,p_total,Teff,Rhs)
