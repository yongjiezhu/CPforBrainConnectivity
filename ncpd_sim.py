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

Sim = scio.loadmat("simulatedData.mat")
X = Sim['Sim']
# Fit CP tensor decomposition .
U = tt.ncp_hals(X, rank=3, verbose=True) # _bcd or _als for options
fig, ax, po = tt.plot_factors(U.factors)
plt.show()
# convert the vector of connections to adjacent matrix for visulazation

