import csv
import os
import matplotlib.pyplot as plt
import sys
import matplotlib
import numpy as np

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
x = np.asarray(x).astype(float)
print "done with x"

vel1 = list( csv.reader(open('vel1.csv')) )
vel1 = np.asarray(vel1).astype(float)
print "done with vel"
vel2 = list( csv.reader(open('vel2.csv')) )
vel2 = np.asarray(vel2).astype(float)
print "done with vel"
vel3 = list( csv.reader(open('vel3.csv')) )
vel3 = np.asarray(vel3).astype(float)
print "done with vel"
vel4 = list( csv.reader(open('vel4.csv')) )
vel4 = np.asarray(vel4).astype(float)
print "done with vel"

den1 = list( csv.reader(open('den1.csv')) )
den1 = np.asarray(den1).astype(float)
print "done with den"
den2 = list( csv.reader(open('den2.csv')) )
den2 = np.asarray(den2).astype(float)
print "done with den"
den3 = list( csv.reader(open('den3.csv')) )
den3 = np.asarray(den3).astype(float)
print "done with den"
den4 = list( csv.reader(open('den4.csv')) )
den4 = np.asarray(den4).astype(float)
print "done with den"

temp1 = list( csv.reader(open('temp1.csv')) )
temp1 = np.asarray(temp1).astype(float)
print "done with temp"
temp2 = list( csv.reader(open('temp2.csv')) )
temp2 = np.asarray(temp2).astype(float)
print "done with temp"
temp3 = list( csv.reader(open('temp3.csv')) )
temp3 = np.asarray(temp3).astype(float)
print "done with temp"
temp4 = list( csv.reader(open('temp4.csv')) )
temp4 = np.asarray(temp4).astype(float)
print "done with temp"

efield = list( csv.reader(open('efield.csv')) )
efield = np.asarray(efield).astype(float)
print "done with efield"

# friction = list( csv.reader(open('friction.csv')) )
# friction = np.asarray(friction).astype(float)
# print "done with friction"

xrange = [0.,L];
nfiles = np.int( x.size / nz )

n_total = int(np.ceil(len(x)/float(nz)))

print nz, n_total


print "rmin = ", min(x)

time = np.array([])
mix = np.array([])

mp_g = 1*1836*9.11e-28;
m_Ti = 47.87 * mp_g #mass of Titanium

for n in range(1,n_total):

	print n
	
	time = np.append(time,dt_print * (n-1))

	start = (n - 1) * nz; finish = start + nz
	
	mu = x[start:finish]
	u1 = vel1[start:finish]
	u2 = vel2[start:finish]
	u3 = vel3[start:finish]
	ue = vel4[start:finish]
	n1 = den1[start:finish]
	n2 = den2[start:finish]
	n3 = den3[start:finish]
	ne = den4[start:finish]
	T1 = temp1[start:finish]
	T2 = temp2[start:finish]
	T3 = temp3[start:finish]
	Te = temp4[start:finish]
	Efield = efield[start:finish]
	p1 = n1 * 1.e6 * T1 * 1.6e-19 
	pe = ne * 1.e6 * Te * 1.6e-19 

	mu = 1.e2 * mu
	mu = np.flipud(mu)
	n2 = np.flipud(n2)
 
	mx = np.trapz(4*np.pi*mu[mu<=55.e-4]**2*n2[mu<=55.e-4]*m_Ti,mu[mu<=55.e-4])

	mix = np.append(mix,mx)
 


plt.semilogy(time/1.e-9,mix/1.e-9)
plt.plot([0,10],[1.3,1.3])	
plt.xlabel('Time [ns]')
plt.ylabel('Mass of Ti in core [ng]')
plt.ylim([0,5])
plt.xlim([0,1.5])

plt.show()

if mix[-1]/1.e-9>1.2:
	f = np.where(abs(mix-1.e-9*1.2)==np.min(abs(mix-1.e-9*1.2)))

start = (int(f[0]) - 1) * nz; finish = start + nz
	
mu = x[start:finish]
n1 = den1[start:finish]
n2 = den2[start:finish]
n3 = den3[start:finish]
ne = den4[start:finish]

plt.figure()
plt.plot(1.e6*mu,1.e-19*n2,'g',linewidth=2)
plt.xlim([0,50])
plt.ylim([0,2.])
plt.xlabel('Radius [$\mu$m]')
plt.ylabel('n$_{Ti}$ [$10^{19}$ cm-3]')
plt.savefig('Ti_mix.pdf',format='pdf')	

plt.show()

