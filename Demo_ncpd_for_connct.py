# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 16:00:32 2020

===============================================================
Simulation of Non-negative CP decomposition (ncpd) for time-freq connectivity
===============================================================

@author: Yongjie Zhu
"""

import tensortools as tt
import matplotlib.pyplot as plt
import scipy.io as scio
import numpy as np
from mne.viz import  plot_connectivity_circle


Sim = scio.loadmat("simulatedData.mat")
X = Sim['Sim'] # X[i,j,k] represents the connectivity value in time points k, freq bin j and paires (nodes) i
freq_index=np.linspace(3,45,169)
# Fit CP tensor decomposition .
U = tt.ncp_hals(X, rank=3, verbose=True) # _bcd or _als for options
fig, ax, po = tt.plot_factors(U.factors)
plt.show()
# visulazation
Nroi = 64 # number of rois
R = 3 # number of components
for index in range(R):
    # temporal courses
    plt.subplot(R,R,index*R+1)
    plt.plot(U.factors.factors[2][:,index])
    if not index:
        plt.title('Temporal courses')
        plt.xlabel('Time/s')
    # spectra
    plt.subplot(R,R,index*R+2)
    plt.plot(freq_index, U.factors.factors[1][:,index])
    if not index:
        plt.title('Spectra')
        plt.xlabel('Frequence/Hz')
    # connectivity
    # convert the vector of connections to adjacent matrix for visulazation
    Corrm = np.zeros([Nroi,Nroi])
    flag = 0
    for num in range(Nroi-1):
        Corrm[num:,num] = U.factors.factors[0][flag+1:flag+Nroi-num+1,index]
        flag = flag+Nroi-num
    plt.subplot(R,R,index*R+3)
    plt.imshow(Corrm)
#    plot_connectivity_circle(Corrm, node_Name, n_lines=50,
#                         title='All-to-All Connectivity left-Auditory '
#                               'Condition (DTI)')
    if not index:
        plt.title('Adjacent matrix')
        plt.xlabel('Nodes')
        plt.ylabel('Nodes')
    
