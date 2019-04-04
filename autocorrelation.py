import glob
import numpy as np

euclraw = glob.glob('./risultati-imfreq-free-64/IMFREQmean.*')
euclraw.sort()
eucl = euclraw[:4000]
data = np.loadtxt('./risultati-imfreq-free-64/IMFREQmean')
autocorrel = np.ndarray(shape=(len(data[:,0]),len(eucl)),dtype=complex)
coefficients = np.ndarray(shape=(len(data[:,0])),dtype=int)

for i,file in enumerate(eucl):
	corr = np.loadtxt('%s'%file, dtype=complex)
	autocorrel[:,i] = corr[:,2]-data[:,2]

for i in xrange(len(data[:,0])):
	autocorrel[i,:] = np.fft.ifft(autocorrel[i,:])
	autocorrel[i,:] = autocorrel[i,:]*np.conj(autocorrel[i,:])
	autocorrel[i,:] = np.fft.fft(autocorrel[i,:])

for j in xrange(len(data[:,0])):
	for i in xrange(1,len(eucl)):
		if autocorrel[0,i] < autocorrel[0,i]*0.08:
			coefficients[j]=i
		break

np.savetxt("nope.txt",coefficients)
