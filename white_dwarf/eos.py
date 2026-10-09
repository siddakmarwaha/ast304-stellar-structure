########################################################################
# Team <your team name>: <names>
# AST304, Fall 2020
# Michigan State University
########################################################################

"""
This module provides functions to compute the electron degeneracy 
pressure and mass density for white dwarf stars, based on the equation 
of state for relativistic electron gases. 
"""

import astro_const as ac
import numpy as np

def pressure(rho, mue):
    """
    Arguments
        rho
            mass density (kg/m**3)
        mue
            baryon/electron ratio
    
    Returns
        electron degeneracy pressure (Pascal)
    """
    
    # replace following lines with body of routine
    p = 1/5 * (3/(8*np.pi))**(2/3) * ac.h**2/ac.m_e * (rho/(mue*ac.m_u))**(5/3)
    return p

def density(p, mue):
    """
    Arguments
        p
            electron degeneracy pressure (Pascal)
        mue
            baryon/electron ratio
        
    Returns
        mass density (kg/m**3)
    """
    
    # replace following lines with body of routine
    rho = (5*p * (8*np.pi/3)**(2/3) * ac.m_e/ac.h**2 * (mue*ac.m_u)**(5/3))**(3/5)
    return rho