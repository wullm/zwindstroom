from zwindstroom import *
import numpy as np

# The neutrino species
M_nu = [0.06] # eV
deg_nu = [1.0] # degeneracies
N_nu = len(M_nu)

# Initialise a unit system (default uses Mpc lengths and km/s velocities)
unit_system, physical_consts = units.init_units()

# We want to integrate the cosmological tables starting at this scale factor
a_start = 1e-3

# Set up a cosmological model
params = {"h": 0.681,
          "Omega_b": 0.0486,
          "Omega_c": 0.2560110606,
          "N_nu": N_nu,
          "M_nu": M_nu,
          "deg_nu": deg_nu,
          "N_ur": 2.0308,
          "T_nu_0": 2.7255 * 0.71611,
          "T_CMB_0": 2.7255,
          "w0": -1.0}
model = cosmology.MODEL()
model.set(params)
model.compute(unit_system, physical_consts, a_start)

# Prepare a table of scale factors and Hubble rates
a_tab_start = 1e-3
a_tab_final = 1.0
na = 1000
avec = np.exp(np.linspace(np.log(a_tab_start), np.log(a_tab_final), na))
Hvec = np.zeros_like(avec)

# Extract Hubble rate as a function of time
for i in range(na):
    Hvec[i] = model.get_H_of_a(avec[i])

# Write the output to a text file
H_table = np.array([avec, Hvec])
header = "Hubble rates between a_begin = %g and a_end = %g" % (a_tab_start, a_tab_final) + "\n"
header += "H(a) in inverse time unit 1/U_T = %g s^-1 = %g km/s/Mpc\n" % (1.0 / unit_system.UnitTimeSeconds, units.__kms / unit_system.UnitTimeSeconds)
header += "Neutrino masses: [" + "".join("%g, " % m for m in M_nu)[:-2] + "] eV\n"
header += "Degeneracies: [" + "".join("%g, " % m for m in deg_nu)[:-2] + "]\n"
header += "H_0: %g 1/U_T = %g km/s/Mpc\n\n" % ((model.h * 100) / (units.__kms / unit_system.UnitTimeSeconds), (model.h * 100))
header += "a H(a)"
fname = "hubble.txt"
np.savetxt(fname, H_table.T, header = header)

print("Table with Hubble rates written to", fname)
    
