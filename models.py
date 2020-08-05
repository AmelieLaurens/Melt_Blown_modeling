#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 16:03:49 2020

@author: amelielaurens

Reference : Modeling Polymer Air Drawing in the Melt Blowing Nonwoven Process
by Chen et Huang, 2003

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
    """Melt Blown Rheological force

    Fr=D**2*pi*(Cxx-Cyy)/4

    :Input:
    - *D* (float) - fiber diameter (m)
    - *Cxx* (float) - axial tensile stress of polymer
    - *Cyy* (float) - transversal tensile stress of polymer

    :Returns:
    (float) - rheological force (N)

    """
    return D**2*pi*(Cxx-Cyy)/4

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
    """Melt Blown Nusselt number

    Nu = gamma*Re**k1

    :Input:
    - *Re* (float) - Reynolds number
    - *gamma* (float) - constant
    - *k1* (float) - constant

    :Returns:
    (float) - Nusselt number

    """
    return gamma*Re**k1

def Drag_coefficient(Rey, beta, n):
    """Melt Blown Drag coefficient

    Cf=beta*Re**(-n)

    :Input:
    - *Re* (float) - Reynolds number (m)
    - *beta* (float) - constant of Matsui’s relation
    - *n* (float) - constant of Matsui’s relation

    :Returns:
    (float) - Drag coefficient

    """
    if Rey == 0 :
        return 10000
    else:
        return beta*Rey**(-n)


def momentum_equation(fiber_diameter, delta_x, rho, pol_flow_rate, Fr1, Fr2, gravitational_acceleration, Cf, rho_a, ua, coef_j, fiber_velocity):
    """Melt Blown Momentum equation

    :Input:
    - *fiber_diameter* (float) - fiber diameter (m)
    - *delta_x* (float) - dicretisation on the axis x : delta_x = x_position[1] - x_position[0]
    - *rho* (float) - polymer density (kg/m^3)
    - *pol_flow_rate* (float) - polymer flow rate (kg/s)
    - *Fr1* (float) - Rheological force (N) at the position x(n)
    - *Fr2* (float) - Rheological force (N) at the position x(n+1) = x(n) + delta_x
    - *gravitational_acceleration* (float) - gravitational acceleration (m/s^2)
    - *Cf* (float) - Drag coefficient / air drawing coefficient
    - *rho_a* (float) - air density (kg/m^3)
    - *ua* (float) - x-component of air velocity (m/s)
    - *coef_j* (float) - factor j : this factor deals with the fact that near the spinneret, the airflow
    acts as a positive (downward) force on the polymer melt, but on the far side of the spinneret, 
    the force is negative.
    - *fiber_velocity* (float) -  Fiber velocity (m/s)

    :Returns:
    (float) - Fiber diameter (m) (at a position n+1)

    """
    return fiber_diameter+delta_x*pi*rho*fiber_diameter**3/(8*pol_flow_rate**2)*((Fr1-Fr2)/delta_x+coef_j*pi*fiber_diameter*Cf*rho_a/2*(ua-fiber_velocity)**2-pi*fiber_diameter**2*rho*gravitational_acceleration/4)


def constitutive_equation_Cxx(eta, cst_m, fiber_velocity1, fiber_velocity2, delta_x):
    """Melt Blown constitutive equation Cxx

    :Input:
    - *eta* (float) - Shear viscosity (Pa.s)
    - *cst_m* (float) - Power-law exponent
    - *fiber_velocity1* (float) - Fiber velocity (m/s) at the position x(n)
    - *fiber_velocity2* (float) - Fiber velocity (m/s) at the position x(n+1) = x(n) + delta_x
    - *delta_x* (float) - dicretisation on the axis x : delta_x = x_position[1] - x_position[0]
   
    :Returns:
    (float) - Axial tensile stress of polymer Cxx (Pa)

    """
    if fiber_velocity1==fiber_velocity2:
        return 0
    else:
        u_prime=(fiber_velocity2-fiber_velocity1)/delta_x
        return 2*eta*(u_prime)**cst_m


def constitutive_equation_Cyy(eta, cst_m, fiber_velocity1, fiber_velocity2, delta_x):
    """Melt Blown constitutive equation Cyy

    :Input:
    - *eta* (float) - Shear viscosity (Pa.s)
    - *cst_m* (float) - Power-law exponent
    - *fiber_velocity1* (float) - Fiber velocity (m/s) at the position x(n)
    - *fiber_velocity2* (float) - Fiber velocity (m/s) at the position x(n+1) = x(n) + delta_x
    - *delta_x* (float) - dicretisation on the axis x : delta_x = x_position[1] - x_position[0]
   
    :Returns:
    (float) - Transversal tensile stress of polymer Cyy (Pa)

    """
    if fiber_velocity1==fiber_velocity2:
        return 0
    else:
        u_prime=(fiber_velocity2-fiber_velocity1)/delta_x
        return -eta*(u_prime)**cst_m
