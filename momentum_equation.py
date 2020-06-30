#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 10:55:08 2020

@author: amelielaurens
"""

import numpy
import matplotlib.pyplot as plt
from math import *
from models import *

def momentum_equation(Fiber_diameter, delta_x, rho, pol_flow_rate, Fr, gravitational_acceleration, Cf, rho_a, ua, coef_j, k):
    return (1/Fiber_diameter[k]**2+delta_x*pi*rho/(4*pol_flow_rate**2)*((Fr[k+1]-Fr[k])/delta_x+pi*Fiber_diameter[k]**2*rho*gravitational_acceleration/4-coef_j*pi*Fiber_diameter[k]*Cf[k]*rho_a/2*(ua-4*pol_flow_rate/(pi*Fiber_diameter[k]**2*rho))**2))**(-1/2)


# Discrétisation
N = 20 #Plus N est grand plus la discrétisation est fine
x_hat = numpy.linspace(0, 0.2, N)
delta_x = x_hat[1] - x_hat[0]

# Définition des constantes
rho = 900.
pol_flow_rate = 0.000012
gravitational_acceleration = 9.81
rho_a = 1.20
ua = 25.7
eta = 3000
cst_m = 0.78
ua = 25.7
va = 1.169*10**(-5)
beta=0.78
cst_n = 0.61

# Initialisation : on rentre les conditions initiales
Fiber_diameter = [0.000407]
Cxx = [366757]
Cyy = [-183378]
Fr = [0.007335137]
Re = [275.3]
Cf = [0.0253403]
fiber_velocity = [ua]

for k in range(N//2):
    coef_j=1
    fiber_velocity.append(4*pol_flow_rate/(pi*Fiber_diameter[k]**2*rho))
    Cxx, Cyy = constitutive_equation(eta, cst_m, fiber_velocity, N)
    Fr.append(rheological_force(Fiber_diameter[k], Cxx[k], Cyy[k]))
    Re.append(Reynolds(Fiber_diameter[k], ua, fiber_velocity[k], va))
    #import pdb; pdb.set_trace()
    Cf.append(Drag_coefficient(Re[k], beta, cst_n))
    Fiber_diameter.append(momentum_equation(Fiber_diameter, delta_x, rho, pol_flow_rate, Fr, gravitational_acceleration, Cf, rho_a, ua, coef_j, k))
for k in range(N//2+1,N):
    coef_j=-1
    fiber_velocity.append(4*pol_flow_rate/(pi*Fiber_diameter[k]**2*rho))
    Cxx, Cyy = constitutive_equation(eta, cst_m, fiber_velocity, N)    
    Fr.append(rheological_force(Fiber_diameter[k], Cxx[k], Cyy[k]))
    Re.append(Reynolds(Fiber_diameter[k], ua, fiber_velocity[k+1], va))
    Cf.append(Drag_coefficient(Re[k], beta, cst_n))
    Fiber_diameter.append(momentum_equation(Fiber_diameter, delta_x, rho, pol_flow_rate, Fr, gravitational_acceleration, Cf, rho_a, ua, coef_j, k))

Fiber_diameter=numpy.array(Fiber_diameter)


fig = plt.figure()
axes = fig.add_subplot(1, 1, 1)

for i in range(N):
    axes.plot(x_hat[i], Fiber_diameter[i], 'ro')
axes.grid()
axes.set_xlabel("Axial position x (m)", fontsize=16)
axes.set_ylabel("Fiber diameter (m)", fontsize=16)
axes.set_title("Fiber diameter as a function of the axial position ", fontsize=16)

