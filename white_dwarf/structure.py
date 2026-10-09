########################################################################
# Team <your team name>: <names>
# AST304, Fall 2020
# Michigan State University
########################################################################

"""
This module solves the structure equations for white dwarf stars, including
the mass-radius relationship, using the equation of state (EOS) for degenerate
matter and numerical methods such as the Runge-Kutta integration method (rk4).
The central pressure is guessed using the virial theorem.
"""

import numpy as np
from eos import density, pressure   # fill this in
from ode import rk4    # fill this in
from astro_const import G, Ke  # fill this in

def stellar_derivatives(m,z,mue):
    """
    RHS of Lagrangian differential equations for radius and pressure
    
    Arguments
        m
            current value of the mass
        z (array)
            current values of (radius, pressure)
        mue
            ratio, nucleons to electrons.  For a carbon-oxygen white dwarf, 
            mue = 2.
        
    Returns
        dzdm (array)
            Lagrangian derivatives dr/dm, dP/dm
    """
    
    r, p = z  # Get radius and pressure
    rho = density(p, mue)  # Calculate density using EOS
    
    dzdm = np.zeros_like(z)
    dzdm[0] = 1 / (4 * np.pi * r**2 * rho)  # dr/dm
    dzdm[1] = -G * m / (4 * np.pi * r**4)  # dP/dm
        
    return dzdm

def central_values(Pc,delta_m,mue):
    """
    Constructs the boundary conditions at the edge of a small, constant density 
    core of mass delta_m with central pressure P_c
    
    Arguments
        Pc
            central pressure (units = Pascal)
        delta_m
            infinitesimal core mass (units = kg)
        mue
            nucleon/electron ratio
    
    Returns
        z = array([ r, p ])
            central values of radius and pressure (units = m, Pascal)
    """
    rho_c = density(Pc, mue)
    r = (3 * delta_m / (4 * np.pi * rho_c))**(1/3)
    z = np.array([r, Pc])
    return z
    
def lengthscales(m,z,mue):
    """
    Computes the radial length scale H_r and the pressure length H_P
    
    Arguments
        m
            current mass coordinate (units = kg)
        z (array)
           [ r, p ] (units = m, Pascal)
        mue
            mean electron weight
    
    Returns
        z/|dzdm| (units = m)
    """

    # fill this in
    dzdm = stellar_derivatives(m, z, mue)    
    return z/np.abs(dzdm)
    
def integrate(Pc,delta_m,eta,xi,mue,max_steps=10000):
    """
    Integrates the scaled stellar structure equations

    Arguments
        Pc
            central pressure (units = Pascal)
        delta_m
            initial offset from center (units = kg)
        eta
            The integration stops when P < eta * Pc
        xi
            The stepsize is set to be xi*min(p/|dp/dm|, r/|dr/dm|)
        mue
            mean electron mass
        max_steps
            solver will quit and throw error if this more than max_steps are 
            required (default is 10000)
                        
    Returns
        m_step, r_step, p_step
            arrays containing mass coordinates, radii and pressures during 
            integration (units = kg, m, Pascal)
    """
        
    m_step = np.zeros(max_steps)
    r_step = np.zeros(max_steps)
    p_step = np.zeros(max_steps)
    
    # set starting conditions using central values
    z = central_values(Pc, delta_m, mue)  # Starting values for radius and pressure
    m = delta_m  # Initial mass coordinate

    Nsteps = 0
    for step in range(max_steps):
        radius = z[0]
        pressure = z[1]
        # are we at the surface?
        if (pressure < eta * Pc):
            break
        # store the step
        m_step[Nsteps] = m
        r_step[Nsteps] = radius
        p_step[Nsteps] = pressure
        
        # set the stepsize
        H_r, H_p = lengthscales(m, z, mue)
        h = xi * min(H_r, H_p)
        
        # take a step
        z = rk4(stellar_derivatives, m, z, h, mue)
        m += h  # Update mass coordinate
        
        # increment the counter
        Nsteps += 1
    # if the loop runs to max_steps, then signal an error
    else:
        raise Exception('too many iterations')
        
    return m_step[0:Nsteps], r_step[0:Nsteps], p_step[0:Nsteps]

def pressure_guess(m,mue):
    """
    Computes a guess for the central pressure based on virial theorem and mass-
    radius relation. 
    
    Arguments
        m
            mass of white dwarf (units are kg)
        mue
            mean electron mass
    
    Returns
        P
            guess for pressure
    """
    # fill this in
    Pguess = (G**5 / Ke**4) * (m * mue**2)**(10/3)
    return Pguess