#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 12:38:38 2020

@author: amelielaurens
"""


import numpy
import matplotlib.pyplot as plt
from math import *
from models import *


def Reynolds(Fiber_diameter, ua, fiber_velocity, va):
    
    """Melt Blown Reynolds number
    
    Re=D*abs(ua-u)/va
    
     :Input:
     - *Fiber_diameter* (numpy.ndarray) - fiber diameter (m)
     - *ua* (float) - x-component of air velocity (m/s)
     - *fiber_velocity* (numpy.ndarray) - fiber velocity (m/s)
     - *va* (float) - air kinematic viscosity (m^2/s)
     
     :Returns:
     (numpy.ndarray) - Reynolds number
     
     """
    return Fiber_diameter*abs(ua-fiber_velocity)/va

N=20
fiber_velocity = numpy.linspace(1, 40, N)
rho = 900.
G = 0.000012
Fiber_diameter=[]
for k in range(N):
    Fiber_diameter.append(continuity_equation(G, fiber_velocity[k], rho))
Fiber_diameter = numpy.array(Fiber_diameter)   
ua = 25.7
va = 1.169*10**(-5)
x_hat = numpy.linspace(0, 0.2, N)

fig = plt.figure(figsize=(16,6))
axes = fig.add_subplot(1, 2, 1)
for k in range(N):
    axes.plot(x_hat[k], Reynolds(Fiber_diameter[k], ua, fiber_velocity[k], va), 'bo')          
axes.set_xlabel("Axial position x (m)", fontsize=16)
axes.set_ylabel("Reynolds number", fontsize=16)
axes.set_title("Reynolds number as a function of the axial position", fontsize=18)
axes.grid()
axes.ticklabel_format(axis='y', style='sci', scilimits=(-1,1))