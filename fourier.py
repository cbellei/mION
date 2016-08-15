# -*- coding: utf-8 -*-
"""
Created on Sun Feb 28 23:46:51 2016

@author: admin-bellei
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft, ifft, fftfreq

x = np.linspace(0,20)
noise = 5.*np.random.normal(-1,3,len(x))
y = x**2 + noise


#now filter out noise
freq_max = 0.5 /(x[1]-x[0])
f = fft(y)
freq = fftfreq(len(y),d=x[1]-x[0])

#filter out 
f[np.abs(freq)>0.8*np.max(freq)] = 0

y_filtered = ifft(f)

plt.plot(x,y)
plt.plot(x,y_filtered)
