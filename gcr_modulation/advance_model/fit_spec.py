# Import modules
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Broken power law
def broken_pow_law(T, J, T_b, a1, a2, d):
   return np.log(J * (T / 1000)**a1 / (1.0 + (T / T_b)**(a2/d))**d)

# Import data
voyager = np.loadtxt("Voyager_2012_2015.dat")
pamela = np.loadtxt("PAMELA_2010_2011_2012_2014.dat")
energy_V = voyager[:,0]
flux_V = voyager[:,1]
energy_P = pamela[62:,0] * 1000
flux_P = pamela[62:,1] / 1000
energy_F = np.exp(np.linspace(np.log(energy_V[0]), np.log(energy_P[-1]), num=100))

# Fit data
energy = np.concatenate((energy_V, energy_P))
flux = np.concatenate((flux_V, flux_P))
opt_params, cov = curve_fit(broken_pow_law, energy, np.log(flux),
                            p0=[1000.0, 600, 0.5, 3.0, 5.0], maxfev=10000,
                            bounds=([10.0, 10, 0.0, 2.0, 1.0], [2000.0, 10000, 1.0, 4.0, 10.0]))
print(opt_params)

# Plot
fig, ax = plt.subplots()
ax.loglog(energy_V, flux_V)
ax.loglog(energy_P, flux_P)
ax.loglog(energy_F, np.exp(broken_pow_law(energy_F, *opt_params)))

ax.set(xlabel='Kinetic Energy (MeV)', ylabel='Proton Flux (some units)')

plt.show()