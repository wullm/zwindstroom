"""
Set up a unit system for the cosmological and fluid calculations. By default,
we use GADGET units without reference to the h parameter:
+ Unit Length = Mpc
+ Unit Velocity = km/s
+ Unit Mass = 10^10 M_sol
+ Unit Temperature = K
+ Unit Current = A
"""

import ctypes
import zwindstroom.backend as backend

class UNITS(ctypes.Structure):
    _fields_ = [("UnitLengthMetres", ctypes.c_double),
                ("UnitTimeSeconds", ctypes.c_double),
                ("UnitMassKilogram", ctypes.c_double),
                ("UnitTemperatureKelvin", ctypes.c_double),
                ("UnitCurrentAmpere", ctypes.c_double)]

class PHYSICAL_CONSTS(ctypes.Structure):
    _fields_ = [("SpeedOfLight", ctypes.c_double),
                ("GravityG", ctypes.c_double),
                ("hPlanck", ctypes.c_double),
                ("kBoltzmann", ctypes.c_double),
                ("ElectronVolt", ctypes.c_double),
                ("SoundSpeedNeutrinos", ctypes.c_double)]

# The default unit system
__Mpc = 3.085677581282e22
__kms = 3.08567758148957E+019
__MsolE10 = 1.988435e40
default_units = UNITS()
default_units.UnitLengthMetres = __Mpc # Mpc (no h)
default_units.UnitTimeSeconds = __kms # km/s speeds
default_units.UnitMassKilogram = __MsolE10 # 10^10 M_sol (no h)
default_units.UnitTemperatureKelvin = 1.0 # Kelvin
default_units.UnitCurrentAmpere = 1.0 # Ampere

backend.engine.set_physical_constants.argtypes = [ctypes.POINTER(UNITS), ctypes.POINTER(PHYSICAL_CONSTS)]

def init_units(LengthMetres = default_units.UnitLengthMetres,
               TimeSeconds = default_units.UnitTimeSeconds,
               MassKilogram = default_units.UnitMassKilogram,
               TemperatureKelvin = default_units.UnitTemperatureKelvin,
               CurrentAmpere = default_units.UnitCurrentAmpere):
    """
    Initialise a unit system and compute physical constants.
    Parameters
    ----------
    LengthMetres : double
        unit length in metres (default 1 Mpc)
    TimeSeconds : double
        unit time in seconds (default v = 1 km/s)
    MassKilogram : double
        unit mass in kilogrammes (default 10^10 M_sol)
    TemperatureKelvin : double
        unit temperature in Kelvin (default 1 K)
    CurrentAmpere : double
        unit current in amperes (default 1 A)
    Return
    ------
    (unit_system, physical_consts) : (UNITS, PHYSICAL_CONSTS)
        The unit system and physical constants
    """

    unit_system = UNITS(LengthMetres, TimeSeconds, MassKilogram, TemperatureKelvin, CurrentAmpere)
    physical_consts = PHYSICAL_CONSTS()
    backend.engine.set_physical_constants(ctypes.byref(unit_system), ctypes.byref(physical_consts))

    return (unit_system, physical_consts)
