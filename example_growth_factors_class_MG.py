import zwindstroom.units as units
import zwindstroom.cosmology as cosmology
import zwindstroom.backend as backend
import zwindstroom.fluid as fluid
import zwindstroom.class_link as class_link
import ctypes
import numpy as np

# The neutrino species
M_nu = [0.05, 0.07] # eV
deg_nu = [2.0, 1.0] # degeneracies
N_nu = len(M_nu)

# Paths to external files
path_E = list(b"example_Ea_table.txt")       # Hubble rate
path_G = list(b"example_Geff_table.txt")     # Gravitational constant

# Initialise a unit system (default uses Mpc lengths and km/s velocities)
unit_system, physical_consts = units.init_units()

# We want to integrate the cosmological tables starting at this scale factor
a_start = 1e-3

# Set up a cosmological model
params = {"h": 0.67,
          "Omega_b": 0.048,
          "Omega_c": 0.242,
          "N_ur": 0.0,
          "N_nu": N_nu,
          "M_nu": M_nu,
          "deg_nu": deg_nu,
          "T_nu_0": 1.95,
          "T_CMB_0": 2.728,
          "w0": -1.0,
          "wa": 0.0,
          "has_external_E": 1,
          "path_external_E": path_E,
          "has_external_G_eff": 1,
          "path_external_G_eff": path_G}
model = cosmology.MODEL()
model.set(params)

# Let zwindstroom know how neutrinos are treated in the cosmological simulation
# The options are:
# + 0 (fully relativistic)
# + 1 (relativistic Hubble rate, but constant particle masses)
# + 2 (fully non-relativistic)
model.set_sim_type(1)

print("Integrating cosmological tables.")

# Integrate tables of cosmological background quantities
model.compute(unit_system, physical_consts, a_start)

# Tolerance and step parameter for the fluid integration
tol = 1e-10
hstart = 1e-10

# We want to integrate the fluid equations between these two points
# For rescaling, this would be the starting and final times of the simulation
a_start_fl = 1.0 / 32.0
a_final_fl = 1.0

# Prepare the fluid integrator
fluid.prepare_fluid_integrator(model, unit_system, physical_consts, model.tables, tol, hstart)

print("Running CLASS.")

# Maximum wavenumber for the growth factors (default units = 1 / Mpc)
k_max = 10.0

# Run CLASS on this model
cosmo = class_link.run_class(model, unit_system, a_start_fl, k_max)

# Extract growth factors, growth rates, and wavenumbers from CLASS
k = class_link.get_wavenumbers(model, cosmo, unit_system)
g = class_link.get_growth_rates(model, cosmo, a_start_fl)
D = class_link.get_growth_factors(model, cosmo, a_start_fl)
nk = len(k)

print("Integrating fluid equations.")

# Compute Newtonian growth factors with the fluid approximation
D_cdm = np.zeros(nk)
D_b = np.zeros(nk)
D_nu = np.zeros((nk, N_nu))
for i in range(nk):
    delta_n = [D["d_ncdm[%d]" % j][i] for j in range(N_nu)]
    gn = [g["d_ncdm[%d]" % j][i] for j in range(N_nu)]
    Dn = [0.0 for j in range(N_nu)]

    growth_factors = fluid.GROWTH_FACTORS()
    growth_factors.k = k[i]
    growth_factors.delta_c = D["d_cdm"][i]
    growth_factors.delta_b = D["d_b"][i]
    growth_factors.delta_n = (ctypes.c_double * N_nu)(*delta_n)
    growth_factors.gc = g["d_cdm"][i]
    growth_factors.gb = g["d_b"][i]
    growth_factors.gn = (ctypes.c_double * N_nu)(*gn)
    growth_factors.Dc = 0.
    growth_factors.Db = 0.
    growth_factors.Dn = (ctypes.c_double * N_nu)(*Dn)

    fluid.integrate_fluid_equations(model, unit_system, physical_consts, model.tables, growth_factors, a_start_fl, a_final_fl)

    D_cdm[i] = growth_factors.Dc
    D_b[i] = growth_factors.Db
    D_nu[i,:] = np.array([growth_factors.Dn[i] for i in range(N_nu)])

# Clean up the integrator
fluid.free_fluid_integrator()
cosmology.free_cosmology_tables(model.tables)

print("Done with integrations.")

# Write the output to a text file
D_table = np.append(np.array([k, D_cdm, D_b]), D_nu.T).reshape((3+N_nu,nk))
header = "Growth factors between a_begin = %g and a_end = %g" % (a_start_fl, a_final_fl) + "\n"
header += "Wavenumbers in inverse length unit 1/U_L = %g m^-1 = %g Mpc^-1 (no h)\n" % (1.0 / unit_system.UnitLengthMetres, units.__Mpc / unit_system.UnitLengthMetres)
header += "Neutrino masses: [" + "".join("%g, " % m for m in M_nu)[:-2] + "] eV\n"
header += "Degeneracies: [" + "".join("%g, " % m for m in deg_nu)[:-2] + "]\n\n"
header += "k D_cdm D_b " + "".join("D_nu[%g] " % d for d in range(N_nu))[:-1]
fname = "table_MG.txt"
np.savetxt(fname, D_table.T, header = header)

print("Table with growth factors written to", fname)
