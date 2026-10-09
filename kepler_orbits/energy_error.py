import numpy as np
import matplotlib.pyplot as plt
from kepler_template import set_initial_conditions, integrate_orbit

# Q4

# Set parameters
a = 1.0
m = 1.0
e = 0.5

# Get initial conditions, energy, and period
z0, eps0, Tperiod = set_initial_conditions(a, m, e)

# Initial step size
h0 = 0.1 * Tperiod

# Number of periods to integrate
num_periods = 3
tend = num_periods * Tperiod

# Step size factors
h_factors = [1, 1/2, 1/4, 1/8, 1/16, 1/32, 1/64, 1/128, 1/256, 1/512, 1/1024]

# Integration methods to compare
methods = ['Euler', 'RK2', 'RK4']

# Prepare to store results
errors = {method: [] for method in methods}

# Loop through each integration method
for method in methods:
    for factor in h_factors:
        h = h0 * factor
        
        # Integrate the orbit
        ts, Xs, Ys, KEs, PEs, TEs = integrate_orbit(z0, m, tend, h, method=method)
        
        # Compute relative error in energy at the end
        E_initial = TEs[0]
        E_final = TEs[-1]
        relative_error = np.abs(E_final - E_initial) / np.abs(E_initial)
        
        # Store the result
        errors[method].append(relative_error)

# Plot the error vs step size (log plot)
plt.figure()
for method in methods:
    plt.loglog(h0 * np.array(h_factors), errors[method], label=method)

plt.xlabel('Step size (h)')
plt.ylabel('Relative Energy Error')
plt.legend()
plt.title('Energy Error vs Step Size for Different Methods')
plt.grid()
plt.show()


# It does scale as expected and it's better to use the logarithmic scale instead of linear.
# That's because the logarithmic scale is separating the three clearly in lower step sizes, while the linear scale was not.

# For Q5

h_large = h0  # Largest
h_small = h0 / 1024  # Smallest

# Plotting for each method
for method in methods:
    for h in [h_large, h_small]:
        ts, Xs, Ys, KEs, PEs, TEs = integrate_orbit(z0, m, tend, h, method=method)

        # Plot trajectory (x vs y)
        plt.figure()
        plt.plot(Xs, Ys, label=f'{method} with h = {h}')
        plt.xlabel('x (AU)')
        plt.ylabel('y (AU)')
        plt.title(f'Trajectory for {method} (h = {h})')
        plt.grid(True)
        plt.axis('equal')
        plt.legend()
        plt.show()

        # Plot energies (KE, PE, TE over time)
        plt.figure()
        plt.plot(ts, KEs, label='Kinetic Energy')
        plt.plot(ts, PEs, label='Potential Energy')
        plt.plot(ts, TEs, label='Total Energy')
        plt.xlabel('Time (years)')
        plt.ylabel('Energy')
        plt.title(f'Energies for {method} (h = {h})')
        plt.grid(True)
        plt.legend()
        plt.show()
