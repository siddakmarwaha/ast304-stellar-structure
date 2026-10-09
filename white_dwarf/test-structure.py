import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import brentq
from structure import integrate, pressure_guess
from astro_const import Msun, Rsun, G
from observations import MassRadiusObservations
from eos import density, pressure

def compute_mass_radius_relation(Mwant, mue=2, eta=0.01, xi=0.01):
    """
    Compute the mass-radius relation for a white dwarf with a desired mass Mwant.

    Parameters:
    - Mwant: Desired mass in kg
    - mue, eta, xi: EOS and numerical parameters

    Returns:
    - radius (in solar radii), central pressure (in Pascal)
    """
    def f(Pc, Mwant):
    # Calculate the total mass for a given Pc
        m_step, r_step, p_step = integrate(Pc, delta_m=1e-8 * Mwant, eta=eta, xi=xi, mue=mue)
        total_mass = m_step[-1]
        return total_mass - Mwant
        
    # Central pressure bounds for brentq
    P_guess = pressure_guess(Mwant, mue)
    Pc_min = 1e-2 * P_guess
    Pc_max = 1e2 * P_guess

    # Find the correct central pressure using brentq (Finds a root of a function in a bracketing interval using Brent’s method)
    Pc_solution = brentq(f, Pc_min, Pc_max, args=(Mwant), xtol=1e-5)

    # Integrate to find radius and mass for Pc_solution
    m_step, r_step, p_step = integrate(Pc_solution, delta_m=1e-8 * Mwant, eta=eta, xi=xi, mue=mue)
    total_radius = r_step[-1]

    return (total_radius / Rsun), Pc_solution

def generate_mass_radius_table():
    """
    Generate a table of mass-radius relations for white dwarf stars.
    """
    masses = np.linspace(0.1, 1.0, 10)  # Masses from 0.1 to 1.0 M_sun
    mass_radius_data = []

    for M in masses:
        M_kg = M * Msun  # Convert mass to kg
        radius, Pc = compute_mass_radius_relation(M_kg)  # Compute radius and central pressure
        Pc_dimensionless = Pc / (G * M_kg**2 * radius**(-4))  # Pc / (GM^2R^(-4))
        rho_c = density(Pc, 2)  # in kg/m^3
        rho_c_dimensionless = rho_c / (3 * M_kg / (4 * np.pi * radius**3))  # rho_c / [3M/(4πR^3)]
        mass_radius_data.append([M, radius, Pc, Pc_dimensionless, rho_c, rho_c_dimensionless])
    return np.array(mass_radius_data)

def plot_mass_radius_relation():
    """
    Plot the mass-radius relation for model white dwarfs and compare it with observational data.
    """
    # Get the mass-radius data for model white dwarfs
    mass_radius_data = generate_mass_radius_table()
    # Load observational data
    obs = MassRadiusObservations()
    # Plot
    plt.figure(figsize=(8,6))
    plt.plot(mass_radius_data[:, 0], mass_radius_data[:, 1], label="Model", color='b', marker='o', linestyle='--')
    plt.scatter(obs.masses, obs.radii*1e-2, color='r', label="Observations")
    plt.xlabel('Mass (M_sun)')
    plt.ylabel('Radius (R_sun)')
    plt.legend()
    plt.grid()
    plt.show()
    # Convert data into a dataframe to display the required table
    mass_radius_data = pd.DataFrame(
        mass_radius_data,
        columns=[
            "M/Msun",
            "R/Rsun",
            "Pc (MKS)",
            "Pc/(GM^2R^(-4))",
            "rho_c (MKS)",
            "rho_c/[3M/(4*pi*R3)]",
        ]
    )
    print("The table required is: \n", mass_radius_data)
plot_mass_radius_relation()