# -*- coding: utf-8 -*-
"""
Created on Sun Feb 28 15:55:52 2021

@author: aaftab
"""

import sympy as sp
import numpy as np
import math
from scipy.integrate import quad
import math
# air density (kg.m-3) and dynamic viscosity (Pa.s) at cruise altitude
rho_air, mu = 0.089, 1.4216e-5
# air speed (m.s-1) at cruise altitude
v = 20
#interference factor
interference_factor=2
l=300.1
a,b,c,d=7.447,2.072,9.010,7.981
x=sp.symbols("x")
def f(x):
    return a*(l-x)*(b*x-l*c**0.5+(c*l**2-d*l*x)**0.5)/64
j,err=quad(f,0,300.1)
volume=j*np.pi

import scipy.integrate as spi


def surface_area(f,a0,b0,h=0.00001,N=10000):
    '''Approximate the arc length of y=f(x) from x=a to x=b.

    Parameters
    ----------
    f : (vectorized) function of one variable
    a,b : numbers defining the interval [a,b]
    h : step size to use in difference formulas
    N : number of subintervals in trapezoid method

    Returns
    -------
    Approximation of the integral \int_a^b \sqrt{1 + (f'(x))^2} dx
    representing the arc length of y=f(x) from x=a to x=b.
    '''
    l=300.1
    p,q,c,d=7.447,2.072,9.010,7.981
    x = np.linspace(a0,b0,N+1)
   
    i=(p*(l-x)*(q*x-l*c**0.5+(c*l**2-d*l*x)**0.5))**0.5/8
    y = f(x)
    # Compute central difference formula for x_k for 1 &lt;= k &lt;= N-1
    h = np.min([h,(b0-a0)/N]) # Make sure that h is smaller than the size of the subintervals
    x_interior = x[1:-1]
    df_interior = (f(x_interior + h) - f(x_interior - h))/(2*h)

    # Use forward/backward difference formula at the endpoints
    df_a0 = (f(a0 + h) - f(a0))/h
    df_b0 = (f(b0) - f(b0 - h))/h
    df = np.hstack([[df_a0],df_interior,[df_b0]])

    # Compute values of the integrand in arc length formula
    y = i*np.sqrt(1 + df**2)

    # Compute the integral
    Surface_area = spi.trapz(y,x)
    

    return Surface_area
surface_area=2*np.pi*surface_area(lambda x:np.sqrt(7.447*(300.1-x)*(2.072*x-300.1*9.010**0.5+(9.010*300.1**2-7.981*300.1*x)**0.5))/8,0,300.1)
surface_area
def surface_area_of_solar_cell():
    cov_ratio=0.45
    return cov_ratio*surface_area
# helium density(kg.m-3) t cruise altitude
rho_he=0.012287
def buo():
    return (rho_air-rho_he)*volume
def CDV():
    """ Calculate the drag coefficient. """
    Re = rho_air * v * l / mu       # Reynold's number
    r = 3.87                  # "Fineness" ratio
    return (0.172 * r**(1/3) + 0.252 / r**1.2 + 1.032 / r**2.7)/Re**(1/6)
CDV()

# Drag calculation
def D():                   
    """ Return the total drag on the airship envelope. """
    return 0.5 * interference_factor * rho_air * v**2 * volume**(2/3) * CDV()
# power required by payload
power_by_payload=10000
# power required by control system is equial to 14% of payload power
control_system_power=1000

pro_eff=0.85   #propulsion efficiency
gear_eff=0.8   #gear efficiency
def thrust_power():
    return D()*v/(pro_eff*gear_eff)
def total_req_power():
    return power_by_payload+thrust_power()+control_system_power

# weight of solar cell configuration
def w_sol_array():                   
    cov_ratio=0.45
    mass_density_of_sol_cells=0.25
    
    return 1.3*surface_area*mass_density_of_sol_cells*cov_ratio

# propulsion system calculation
def wt_pro_sys():         
    
    power_density=75
    return thrust_power()/power_density

def wt_rfc():             # weight of regenertive fuel cell
    rfc_eff=0.6
    rfc_density=250
    irradiance=480
    discharge_time=14
    
    
    return total_req_power()*discharge_time/(rfc_eff*rfc_density)
volume_hull=volume
surface_hull=surface_area
mass_energy=w_sol_array()+wt_rfc()
mass_thrust=wt_pro_sys()
def structural_weight():
    
    mol_wt_of_air=28.97  # molecular weght of air in kg/kmol
    mol_wt_of_He=4.003    # molecular weght of helium in kg/kmol
    rho_fabric =0.2      # density of fabric used in volume
    surface_fin=0.0121*volume_hull
    mass_of_gas=rho_air*(mol_wt_of_He/mol_wt_of_air)*volume_hull
    mass_of_hull=1.2*rho_fabric*surface_hull
    mass_of_fin=1.2*rho_fabric*surface_fin
    mass_of_other_component=0.25*(mass_of_hull+mass_of_fin+mass_thrust+mass_energy)
    return mass_of_gas+mass_of_hull+ mass_of_fin+mass_of_other_component
def total_weight():
    payload_weight=1000
    return structural_weight()+wt_pro_sys()+mass_energy+payload_weight


print('volume:', volume,
      '\nsurface area of envelope',surface_area,
      '\nsurface area of solar array:',surface_area_of_solar_cell(),
      '\nthrust required(drag):',D() ,'volumetric drag coefficient:',CDV(),
      '\nthrust power required:',thrust_power(),
      '\ntotal required power:',total_req_power(),
      '\nbuoyancy:',buo(),
      '\nweight of solar array:',w_sol_array(),
      '\nweight of propulsion system:',wt_pro_sys(),
      '\nweight of rfc system:',wt_rfc(),
      '\nweight of sturcture',structural_weight(),
      '\ntotal weight:',total_weight(),
      )

