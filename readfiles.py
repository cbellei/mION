import numpy as np
import matplotlib.pyplot as plt
import sys

try:
	nz = int(sys.argv[1])
	minden = float(sys.argv[2])
except:
	pass

def readfiles(nz=1000,minden=1.e19*1.e6):
	"""nz: number of grid points
	minden: minimum ion number density [m-3])"""

	rho1 = np.array([]) #C ions
	rho2 = np.array([]) #D ions
	rho3 = np.array([]) #3He ions
	rho4 = np.array([]) #electrons
	r = np.array([])
	temp1 = np.array([])
	temp2 = np.array([])
	temp3 = np.array([])
	temp4 = np.array([]) #electrons
	vel1 = np.array([])
	vel2 = np.array([])
	vel3 = np.array([])
	vel4 = np.array([]) #electrons

	Z1 = 2 #Carbon
	Z2 = 1 #D
	Z3 = 2 #3He
	A1 = 12
	A2 = 2
	A3 = 3

	m1 = A1*1836*9.11e-31 #kg
	m2 = A2*1836*9.11e-31
	m3 = A3*1836*9.11e-31 #kg
	me = 9.11e-31 #kg

	qe = 1.60218e-19 #Coulomb
	
	with open('eden.dat') as f:
		for i,line in enumerate(f):
			data = np.asarray(line.split(),dtype=float)
			data[1] *= 1.e6 #now it is in m-3
			r = np.append(r,1.e-2*data[0])
			if i<=100:
				rho1 = np.append(rho1,m1*minden)
				rho2 = np.append(rho2,m2*minden)
				rho3 = np.append(rho3,m3*data[1]/Z3)
			else:
				rho1 = np.append(rho1,m1*data[1]/(Z1+Z2))
				rho2 = np.append(rho2,m2*data[1]/(Z1+Z2))
				rho3 = np.append(rho3,m3*minden)

	#electrons			
	rho4 = me * ( Z1*rho1/m1 + Z2*rho2/m2 + Z3*rho3/m3 )
			
	with open('etemp.dat') as f:
		for i,line in enumerate(f):
			data = np.asarray(line.split(),dtype=float)
			temp4 = np.append(temp4,1.e3*qe*data[1])
			
	with open('itemp.dat') as f:
		for i,line in enumerate(f):
			data = np.asarray(line.split(),dtype=float)
			temp1 = np.append(temp1,1.e3*qe*data[1])
			temp2 = np.append(temp2,1.e3*qe*data[1])
			temp3 = np.append(temp3,1.e3*qe*data[1])
		
	with open('vel.dat') as f:
		for i,line in enumerate(f):
			data = np.asarray(line.split(),dtype=float)
			vel1 = np.append(vel1,1.e3*data[1])
			vel2 = np.append(vel2,vel1[i])
			vel3 = np.append(vel3,vel1[i])

	#electrons
	vel4 = ( Z1*rho1*vel1/m1 + Z2*rho2*vel2/m2 + Z3*rho3*vel3/m3 ) / (rho4/me)

	L = r[-1]
	rmin = 2.5e-6

	#now interpolates to nz
	R = np.linspace(rmin,L,nz)
	Rho1 = np.interp(R,r,rho1)
	Rho2 = np.interp(R,r,rho2)
	Rho3 = np.interp(R,r,rho3)
	Rho4 = np.interp(R,r,rho4)
	Temp1 = np.interp(R,r,temp1)
	Temp2 = np.interp(R,r,temp2)
	Temp3 = np.interp(R,r,temp3)
	Temp4 = np.interp(R,r,temp4)
	Vel1 = np.interp(R,r,vel1)
	Vel2 = np.interp(R,r,vel2)
	Vel3 = np.interp(R,r,vel3)
	Vel4 = np.interp(R,r,vel4)

	#flip orientation
	R = np.flipud(R)
	Rho1 = np.flipud(Rho1)
	Rho2 = np.flipud(Rho2)
	Rho3 = np.flipud(Rho3)
	Rho4 = np.flipud(Rho4)
	Temp1 = np.flipud(Temp1)
	Temp2 = np.flipud(Temp2)
	Temp3 = np.flipud(Temp3)
	Temp4 = np.flipud(Temp4)
	Vel1 = np.flipud(Vel1)
	Vel2 = np.flipud(Vel2)
	Vel3 = np.flipud(Vel3)
	Vel4 = np.flipud(Vel4)

	#write out files
	f1 = 'r.dat'
	f2 = 'rho1.dat'
	f3 = 'rho2.dat'
	f4 = 'rho3.dat'
	f5 = 'rho4.dat'
	f6 = 'temp1.dat'
	f7 = 'temp2.dat'
	f8 = 'temp3.dat'
	f9 = 'temp4.dat'
	f10 = 'vel1.dat'
	f11 = 'vel2.dat'
	f12 = 'vel3.dat'
	f13 = 'vel4.dat'

	with open(f1,'w') as f:
		for i in range(nz):
			f.write('%.5e' % (R[i])+'\n')		
	with open(f2,'w') as f:
		for i in range(nz):
			f.write('%.5e' % (Rho1[i])+'\n')		
	with open(f3,'w') as f:
		for i in range(nz):
			f.write('%.5e' % (Rho2[i])+'\n')		
	with open(f4,'w') as f:
		for i in range(nz):
			f.write('%.5e' % (Rho3[i])+'\n')		
	with open(f5,'w') as f:
		for i in range(nz):
			f.write('%.5e' % (Rho4[i])+'\n')		
	with open(f6,'w') as f:
		for i in range(nz):
			f.write('%.5e' % (Temp1[i])+'\n')		
	with open(f7,'w') as f:
		for i in range(nz):
			f.write('%.5e' % (Temp2[i])+'\n')		
	with open(f8,'w') as f:
		for i in range(nz):
			f.write('%.5e' % (Temp3[i])+'\n')		
	with open(f9,'w') as f:
		for i in range(nz):
			f.write('%.5e' % (Temp4[i])+'\n')		
	with open(f10,'w') as f:
		for i in range(nz):
			f.write('%.5e' % (Vel1[i])+'\n')		
	with open(f11,'w') as f:
		for i in range(nz):
			f.write('%.5e' % (Vel2[i])+'\n')		
	with open(f12,'w') as f:
		for i in range(nz):
			f.write('%.5e' % (Vel3[i])+'\n')		
	with open(f13,'w') as f:
		for i in range(nz):
			f.write('%.5e' % (Vel4[i])+'\n')		


	plt.figure()
	plt.subplot(2,2,1)
	plt.semilogy(1.e6*R,Z1*Rho1/m1+Z2*Rho2/m2+Z3*Rho3/m3,1.e6*R,Rho4/me,'r--')
	plt.subplot(2,2,2)
	plt.semilogy(1.e6*r,rho1,1.e6*R,Rho1)
	plt.semilogy(1.e6*r,rho2,1.e6*R,Rho2)
	plt.semilogy(1.e6*r,rho3,1.e6*R,Rho3)
	plt.semilogy(1.e6*r,rho4,1.e6*R,Rho4,'k')
	plt.subplot(2,2,3)
	plt.plot(1.e6*r,temp1,1.e6*R,Temp1)
	plt.plot(1.e6*r,temp2,1.e6*R,Temp2)
	plt.plot(1.e6*r,temp3,1.e6*R,Temp3)
	plt.plot(1.e6*r,temp3,1.e6*R,Temp4,'k')
	plt.subplot(2,2,4)
	plt.plot(1.e6*r,vel1,1.e6*R,Vel1)
	plt.plot(1.e6*r,vel2,1.e6*R,Vel2)
	plt.plot(1.e6*r,vel3,1.e6*R,Vel3)
	plt.plot(1.e6*r,vel4,1.e6*R,Vel4,'k')
	plt.show()

if __name__ == "__main__":
	print "------------"
	print "nz = " + str(nz)
	print "minden = " + str(minden)
	print "------------"
	readfiles(nz,minden)

		