# -*- coding: utf-8 -*-
"""
Created on Wed Feb 9 01:13:00 2022

Runfile of DP4GINO software -- Linux Version

@author: Xuyu Chen 
@Insititution: School of Mathematical Sciences, Fudan University
@Email: chenxy18@fudan.edu.cn
"""
import lib_DP4GNIO_Linux as gnio
import numpy as np

# An example on how to use the software
n = 100000
data = np.random.uniform(-100, 100, n)
w = 0.5 * np.ones(n)
lbd = np.random.uniform(0, 1000, n-1)
mu = np.random.uniform(0, 1000, n-1)

results_l1 = gnio.l1gnio(data, w, lbd, mu)
results_l2 = gnio.l2gnio(data, w, lbd, mu)

