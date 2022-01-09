"""
Perform first-order fluid calculations for a coupled CDM-Baryon-Neutrino model.
"""

import ctypes
import trifa.backend as backend
import trifa.units as units
import trifa.cosmology as cosmology

class GROWTH_FACTORS(ctypes.Structure):
    _fields_ = [("k", ctypes.c_double),
                ("delta_c", ctypes.c_double),
                ("delta_b", ctypes.c_double),
                ("delta_n", ctypes.POINTER(ctypes.c_double)),
                ("gc", ctypes.c_double),
                ("gb", ctypes.c_double),
                ("gn", ctypes.POINTER(ctypes.c_double)),
                ("Dc", ctypes.c_double),
                ("Db", ctypes.c_double),
                ("Dn", ctypes.POINTER(ctypes.c_double))]

backend.engine.prepare_fluid_integrator.argtypes = [ctypes.POINTER(cosmology.MODEL), ctypes.POINTER(units.UNITS), ctypes.POINTER(units.PHYSICAL_CONSTS), ctypes.POINTER(cosmology.TABLES), ctypes.c_double, ctypes.c_double]

backend.engine.integrate_fluid_equations.argtypes = [ctypes.POINTER(cosmology.MODEL), ctypes.POINTER(units.UNITS), ctypes.POINTER(units.PHYSICAL_CONSTS), ctypes.POINTER(cosmology.TABLES), ctypes.POINTER(GROWTH_FACTORS), ctypes.c_double, ctypes.c_double]

backend.engine.free_fluid_integrator.argtypes = []

def prepare_fluid_integrator(model, units, physical_consts, tables, tol, hstart):
    """
    Prepare the system of ODEs to be integrated.
    Parameters
    ----------
    model : MODEL
        a cosmological model
    units: UNITS
        a unit system
    physical_consts: PHYSICAL_CONSTS
        physical constants corresponding to the unit system
    tables : TABLES
        the cosmological tables to be calculated
    tol : double
        tolerance for the integrator
    hstart : double
        initial step size
    Return
    ------
    bool
        True on success
    """

    backend.engine.prepare_fluid_integrator(ctypes.byref(model), ctypes.byref(units), ctypes.byref(physical_consts), ctypes.byref(tables), tol, hstart)

    return True

def integrate_fluid_equations(model, units, physical_consts, tables, growth_factors, a_start, a_final):
    """
    Integrate the fluid equations.
    Parameters
    ----------
    model : MODEL
        a cosmological model
    units: UNITS
        a unit system
    physical_consts: PHYSICAL_CONSTS
        physical constants corresponding to the unit system
    tables : TABLES
        the cosmological tables to be calculated
    growth_factors: GROWTH_FACTORS
        structure with input and output growth factors and rates (state variables)
    a_start : double
        starting scale factor time of the integration
    a_final : double
        final scale factor time of the integration
    Return
    ------
    bool
        True on success
    """

    backend.engine.integrate_fluid_equations(ctypes.byref(model), ctypes.byref(units), ctypes.byref(physical_consts), ctypes.byref(tables), ctypes.byref(growth_factors), a_start, a_final)

    return True

def free_fluid_integrator():
    """
    Clean up the fluid integrator.
    Return
    ------
    bool
        True on success
    """

    backend.engine.free_fluid_integrator()

    return True
