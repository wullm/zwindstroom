from zwindstroom import *

# The neutrino species
M_nu = [0.05, 0.07] # eV
deg_nu = [2.0, 1.0] # degeneracies
N_nu = len(M_nu)

# Initialise a non-standard unit system. Let's use the CLASS system,
# which is (U_L, U_T) = (Mpc, c/Mpc) for length and time units
unit_system, physical_consts = units.init_units(
                LengthMetres = units.__Mpc,
                TimeSeconds = units.__kms / units.__c_kms,
                MassKilogram = units.__MsolE10,
                TemperatureKelvin = 1.0,
                CurrentAmpere = 1.0)

# We want to integrate the cosmological tables starting at this scale factor
a_start = 1e-3

# Set up a cosmological model
params = {"h": 0.67,
          "Omega_b": 0.048,
          "Omega_c": 0.242,
          "N_nu": N_nu,
          "M_nu": M_nu,
          "deg_nu": deg_nu,
          "T_nu_0": 1.95,
          "T_CMB_0": 2.728,
          "w0": -1.0}
model = cosmology.MODEL()
model.set(params)
model.compute(unit_system, physical_consts, a_start)

H_0 = model.get_H_of_a(1.0)
print("The Hubble rate today is H_0 =", H_0, "1/U_T")
