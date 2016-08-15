import csv
from pylab import *
import os
import matplotlib.pyplot as pyplot
import sys

matplotlib.rc('font', size=24)
matplotlib.rc('font', family='Arial')

nspec = int(sys.argv[1])
maxn = sys.argv[2]
dframe = sys.argv[3]

def movie(nspec=4,maxn=10,dframe=1):

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
	print "done with vel1"
	vel2 = list( csv.reader(open('vel2.csv')) )
	vel2 = np.asarray(vel2).astype(float)
	print "done with vel2"
	vel3 = list( csv.reader(open('vel3.csv')) )
	vel3 = np.asarray(vel3).astype(float)
	print "done with vel3"
	vel4 = list( csv.reader(open('vel4.csv')) )
	vel4 = np.asarray(vel4).astype(float)
	print "done with vel4"
	if (nspec==4):
		vel5 = list( csv.reader(open('vel5.csv')) )
		vel5 = np.asarray(vel5).astype(float)
		print "done with vel5"

	den1 = list( csv.reader(open('den1.csv')) )
	den1 = np.asarray(den1).astype(float)
	print "done with den1"
	den2 = list( csv.reader(open('den2.csv')) )
	den2 = np.asarray(den2).astype(float)
	print "done with den2"
	den3 = list( csv.reader(open('den3.csv')) )
	den3 = np.asarray(den3).astype(float)
	print "done with den3"
	den4 = list( csv.reader(open('den4.csv')) )
	den4 = np.asarray(den4).astype(float)
	print "done with den4"
	if (nspec==4):
		den5 = list( csv.reader(open('den5.csv')) )
		den5 = np.asarray(den5).astype(float)
		print "done with den5"

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
	if (nspec==4):
		temp5 = list( csv.reader(open('temp5.csv')) )
		temp5 = np.asarray(temp4).astype(float)
		print "done with temp"

	efield = list( csv.reader(open('efield.csv')) )
	efield = np.asarray(efield).astype(float)
	print "done with efield"

# 	mu_ion = list( csv.reader(open('mu_ion.csv')) )
# 	mu_ion = np.asarray(mu_ion).astype(float)
# 	print "done with mu_ion"


	# friction = list( csv.reader(open('friction.csv')) )
	# friction = np.asarray(friction).astype(float)
	# print "done with friction"

	xrange = [0.,0.04];#[0.,L];#[0.015,0.025]# #[0.015,0.025]; #
	Trange = [0., 100.]#1.]
	
	nfiles = np.int( x.size / nz )

	n_total = int(ceil(len(x)/float(nz)))

	print nz, n_total

	f1 = figure(figsize=(10,10))
	f2 = figure(figsize=(10,5))


	print "rmin = ", min(x)

	count = 0

	for n in range(1,min(int(maxn),n_total),int(dframe)):

		print n
		count += 1
			
		time = dt_print * (n-1)

		start = (n - 1) * nz; finish = start + nz
		mu = x[start:finish]
		u1 = vel1[start:finish]
		u2 = vel2[start:finish]
		u3 = vel3[start:finish]
		if (nspec==4):
			u4 = vel4[start:finish]
			ue = vel5[start:finish]
		else:
			ue = vel4[start:finish]
		n1 = den1[start:finish]
		n2 = den2[start:finish]
		n3 = den3[start:finish]
		if (nspec==4):
			n4 = den4[start:finish]
			ne = den5[start:finish]
		else:
			ne = den4[start:finish]
		
		T1 = temp1[start:finish]
		T2 = temp2[start:finish]
		T3 = temp3[start:finish]
		if (nspec==4):
			T4 = temp4[start:finish]
			Te = temp5[start:finish]
		else:
			Te = temp4[start:finish]
		
		Efield = efield[start:finish]
		p1 = n1 * 1.e6 * T1 * 1.6e-19 
		pe = ne * 1.e6 * Te * 1.6e-19 
#  		Mu_ion = mu_ion[start:finish]
		
		mu = 1.e2 * mu

		matplotlib.pyplot.figure(1) #select figure 1

		if (nspec==4):
			semilogy(mu,n1/1.e23,linewidth=3,label='CH')
			semilogy(mu,n3/1.e23,linewidth=3,color='r',label='D')
			semilogy(mu,n2/1.e23,linewidth=3,linestyle='--',color='b',label='C')
			semilogy(mu,n4/1.e23,linewidth=3,color='g',label='3He')
		else:
			semilogy(mu,n1/1.e23,linewidth=2,label='CH')
			semilogy(mu,n2/1.e23,linewidth=2,color='r',label='D')
			semilogy(mu,n3/1.e23,linewidth=2,color='g',label='3He')
		
		# semilogy(1.e9*mu,n1/1.e23,linewidth=2)
		xlabel('R [cm]')
		#ylabel('Density [1.e22 cm-3]')
		ylabel('Density [1.e23 cm-3]')
		xlim(xrange)
		#ylim([1.e-8, 20])
		ylim([1.e-6, 5.])
		#title('t =' + str(time))
		plt.legend(loc='lower right')
		plt.grid()
		title(time)
		pyplot.locator_params(axis='x', nbins=4)
			
		filename = 'foo' + str(count).zfill(3) + '.jpg';
		
		print filename

		savefig(filename)
	
		clf()
	
	os.system("ffmpeg -y -i 'foo%03d.jpg' output.m4v")
	#os.system("ffmpeg -y -i 'friction%03d.jpg' friction.m4v")
	#os.system("avconv -i 'foo%03d.jpg' -r 10 movie.avi")
	os.system("rm -f *.jpg")

	print "done"


if __name__ == "__main__":
	print "------------"
	print "nspec = " + str(nspec)
	print "maxn = " + maxn
	print "dframe = " + dframe
	print "------------"
	movie(nspec,maxn,dframe)