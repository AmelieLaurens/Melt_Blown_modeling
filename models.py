#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 16:03:49 2020

@author: amelielaurens
"""

from math import *
import numpy

def continuity_equation(G, u, rho):
    """Melt Blown continuity equation
        
    D=sqrt(4*G/(pi*u*rho))
        
    :Input:
    - *G* (float) - polymer flow rate (kg/s)
    - *u* (float) - fiber velocity (m/s)
    - *rho* (float) - polymer density (kg/m^3)
         
    :Returns:
    (float) - fiber diameter (m)
         
    """
    return sqrt(4*G/(pi*u*rho))

def rheological_force(D, Cxx, Cyy):
    return D**2*pi*(Cxx-Cyy)/4
    """Melt Blown Rheological force
    
    Fr=D**2*pi*(Cxx-Cyy)/4
    
     :Input:
     - *D* (float) - fiber diameter (m)
     - *Cxx* (float) - axial tensile stress of polymer
     - *Cyy* (float) - transversal tensile stress of polymer
     
     :Returns:
     (float) - rheological force (N)
     
     """

def Reynolds(D, ua, u, va):
    
    """Melt Blown Reynolds number
    
    Re=D*abs(ua-u)/va
    
     :Input:
     - *D* (float) - fiber diameter (m)
     - *ua* (float) - x-component of air velocity (m/s)
     - *u* (float) - fiber velocity (m/s)
     - *va* (float) - air kinematic viscosity (m^2/s)
     
     :Returns:
     (float) - Reynolds number
     
     """
    return D*abs(ua-u)/va

def Nusselt(Re, gamma, k1):
    return gamma*Re**k1
    
    """Melt Blown Nusselt number
    
    Nu = gamma*Re**k1
    
     :Input:
     - *Re* (float) - Reynolds number
     - *gamma* (float) - constant
     - *k1* (float) - constant
     
     :Returns:
     (float) - Nusselt number
     
     """
     
def Drag_coefficient(Rey, beta, n):
    return beta*Rey**(-n)
    """Melt Blown Drag_coefficient
    
    Cf=beta*Re**(-n)
    
     :Input:
     - *Re* (float) - Reynolds number (m)
     - *beta* (float) - constant of Matsui’s relation
     - *n* (float) - constant of Matsui’s relation
     
     :Returns:
     (float) - Drag coefficient
     
     """
     
def constitutive_equation(eta, cst_m, fiber_velocity, N):
    Cxx=[]
    Cyy=[]
    x_hat = numpy.linspace(0, 0.2, N)
    delta_x = x_hat[1] - x_hat[0]
    u_prime=0.
    for i in range(1,len(fiber_velocity)):
        u_prime=(fiber_velocity[i]-fiber_velocity[i-1])/delta_x
        Cxx.append(2*eta*(u_prime)**cst_m)
        Cyy.append(-eta*(u_prime)**cst_m)
    Cxx = numpy.array(Cxx)
    Cyy = numpy.array(Cyy)
    return Cxx, Cyy