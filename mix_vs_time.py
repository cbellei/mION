from pylab import *
import os
import matplotlib.pyplot as pyplot
import csv
import os
import sys

rc('font', size=28)
rc('text', usetex=False)
matplotlib.rc('font', family='Arial')

#============================================
matplotlib.rcParams['xtick.major.pad']='14'
matplotlib.rcParams['ytick.major.pad']='15'
figsize1 = 10
figsize2 = 8

fig = figure(figsize=(figsize1,figsize2))
#============================================

param = list( csv.reader(open('parameters.csv')) )[0]

dtm = np.asarray(param[0]).astype(float)
nz  = np.asarray(param[1]).astype(int)
nprint = np.asarray(param[2]).astype(int)
maxind = np.asarray(param[3]).astype(int)
tm_quiet = np.asarray(param[4]).astype(float)
L = 1.e2 * np.asarray(param[5]).astype(float)  # in cm

x = list( csv.reader(open('x.csv')) )
x = np.asarray(x).astype(float)
print "done with x"

den = list( csv.reader(open('den.csv')) )
den = np.asarray(den).astype(float)
print "done with den"

temp = list( csv.reader(open('temp.csv')) )
temp = np.asarray(temp).astype(float)
print "done with temp"

nfiles = x.size / nz

xrange = [0.,L];

logic = 0;

print nfiles

mu12 = np.zeros(nfiles-1)
tt   = np.zeros(nfiles-1)

cont = 0

for i in range(nfiles-1):

    tt[i] = nprint * dtm * i
    
    start = (i + 1) * nz; finish = start + nz
    mu = x[start:finish,0]
    n1 = den[start:finish,0]
    n2 = den[start:finish,1]
    T2 = temp[start:finish,1]

    mu = 1.e2 * mu #cm

    if (i==0):
        n2back = 7.19e19
        dmu = mu[1] - mu[0]
        L_extra = 870e-4 - L/2;
        print "L_extra = ", L_extra
        N_extra = int(L_extra / dmu)

    mu = np.append(mu,np.linspace(L,L+L_extra,N_extra))
    n1 = np.append(n1,np.zeros(N_extra))
    n2 = np.append(n2,n2back * np.ones(N_extra))

    I1 = np.trapz(n1 * n2 / ( n1 + n2), mu)
    I2 = np.trapz(n2 * n2 / ( n1 + n2), mu) #870.e-4 * 7.19e19
#    I2 = 870.e-4 * 7.19e19

    mu12[i] = I1 / I2

#semilogy(mu,n2,'--')
#semilogy(mu,n1 * n2 / ( n1 + n2),linewidth=3)
#semilogy(mu,n1,'--')
#semilogy(mu,n2 * n2 / ( n1 + n2),linewidth=3)

plot(mu,n2,'--')
plot(mu,n1 * n2 / ( n1 + n2),linewidth=3)
plot(mu,n1,'--')
plot(mu,n2 * n2 / ( n1 + n2),linewidth=3)
ylim([1.e18, 1.e21])

#============================================
matplotlib.rcParams['xtick.major.pad']='14'
matplotlib.rcParams['ytick.major.pad']='15'
figsize1 = 10
figsize2 = 8

fig = figure(figsize=(figsize1,figsize2))
#============================================
plot(1.e9 * tt, mu12,'g',linewidth=2)
ylabel('$\mu_{12}(t)$ ')
xlabel('Time [ns]')

plt.tight_layout()

show()




