from zwindstroom import *

# The neutrino species
M_nu = [0.05, 0.07] # eV
deg_nu = [2.0, 1.0] # degeneracies
N_nu = len(M_nu)

# Initialise a unit system (default uses Mpc lengths and km/s velocities)
unit_system, physical_consts = units.init_units()

# We want to integrate the cosmological tables starting at this scale factor
a_start = 1e-3

# Set up a cosmological model
params = {"h": 0.67,
          "Omega_b": 0.048,
          "Omega_c": 0.0,
          "N_nu": N_nu,
          "M_nu": M_nu,
          "deg_nu": deg_nu,
          "T_nu_0": 1.95,
          "T_CMB_0": 2.728,
          "w0": -1.0,
          "Omega_dcdmdr_0": 0.19,
          "Gamma_dcdm": 5.}
model = cosmology.MODEL()
model.set(params)
model.compute(unit_system, physical_consts, a_start)

Omega_dcdm = model.get_Omega_dcdm_of_a(1.0)
Omega_dr = model.get_Omega_dr_of_a(1.0)
print("The present density fractions in decaying dark matter and dark radiation are:")
print("Omega_dcdm:", Omega_dcdm)
print("Omega_dr:", Omega_dr)

print("Total:", Omega_dcdm + Omega_dr)

