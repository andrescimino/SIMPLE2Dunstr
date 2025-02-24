#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 19:20:47 2024

@author: sebastian
"""

import numpy as np
import matplotlib.pyplot as plt
x = np.random.rand(100)
y = np.random.rand(100)
z = np.sin(x)+np.cos(y)
f, ax = plt.subplots(1,2, sharex=True, sharey=True)
ax[0].tripcolor(x,y,z)
ax[1].tricontourf(x,y,z, 20) # choose 20 contour levels, just to show how good its interpolation is
ax[1].plot(x,y, 'ko ')
ax[0].plot(x,y, 'ko ')