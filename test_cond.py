# -*- coding: utf-8 -*-
"""
Created on Sat Mar 19 11:05:03 2016

@author: admin-bellei
"""
import numpy as np

D  = np.array([3.0,2.3,-5.0,-0.9,7.1])
DU = np.array([2.1,-1.0,1.9,8.0])
DL = np.array([3.4,3.6,7.0,-6.0])

A = np.zeros((5,5))
for i in range(5):
    A[i,i] = D[i]
    if(i<4): A[i,i+1] = DU[i]
    if(i>=1): A[i,i-1] = DL[i-1]
    
norm = np.linalg.norm(A,ord=1)
print norm