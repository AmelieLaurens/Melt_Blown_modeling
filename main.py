#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 10:55:08 2020

@author: amelielaurens

Reference : Modeling Polymer Air Drawing in the Melt Blowing Nonwoven Process
by Chen et Huang, 2003

"""

import numpy
import matplotlib.pyplot as plt
from math import *
from models import *

# Discretisation
N = 20   # The higher the discretisation number N is, the finer the discretisation will be,
        # there will be more points on the graphic.
x_position = numpy.linspace(0, 20, N)
delta_x = x_position[1] - x_position[0]

# Constants initialization
rho = 900.
pol_flow_rate = 0.000012
gravitational_acceleration = 9.81
rho_a = 1.20
ua = 25.7
eta = 3000
cst_m = 0.78
ua = 0.257
va = 1.169*10**(-5)
beta=0.78
cst_n = 0.61
coef_j = numpy.linspace(-1, 1, N)

# Initialization : enter the initial conditions
Fiber_diameter = [0.000407]
Cxx = [366757.]
Cyy = [-183378.]
Fr = [0.007335137]
Re = [275.3]
Cf = [0.0253403]
fiber_velocity = [0.257]

# Compute datas
for k in range(N):
    fiber_velocity.append(4*pol_flow_rate/(pi*Fiber_diameter[k]**2*rho))
    Cxx.append(constitutive_equation_Cxx(eta, cst_m, fiber_velocity[k], fiber_velocity[k+1], N))
    Cyy.append(constitutive_equation_Cyy(eta, cst_m, fiber_velocity[k], fiber_velocity[k+1], N))
    Fr.append(rheological_force(Fiber_diameter[k], Cxx[k], Cyy[k]))
    Re.append(Reynolds(Fiber_diameter[k], ua, fiber_velocity[k], va))
    Cf.append(Drag_coefficient(Re[k], beta, cst_n))
    Fr1 = Fr[k]
    Fr2 = Fr[k+1]
    Fiber_diameter.append(momentum_equation(Fiber_diameter[k], delta_x, rho, pol_flow_rate, Fr1, Fr2, gravitational_acceleration, Cf[k], rho_a, ua, coef_j[k], fiber_velocity[k]))
    
Fiber_diameter = numpy.array(Fiber_diameter)
# import pdb; pdb.set_trace()

# Draw the figure 
fig = plt.figure()
axes = fig.add_subplot(1, 1, 1)

for i in range(N):
    axes.plot(x_position[i], Fiber_diameter[i], 'ro')
axes.grid()
axes.set_xlabel("Axial position x (m)", fontsize=16)
axes.set_ylabel("Fiber diameter (m)", fontsize=16)
axes.set_title("Fiber diameter as a function of the axial position ", fontsize=16)
plt.show()
