"""
Set up the cosmological model and compute tables of cosmological quanties,
such as the Hubble rate and neutrino density per species as a function of
time.
"""

import ctypes
import zwindstroom.backend as backend
import zwindstroom.units as units

class TABLES(ctypes.Structure):
    _fields_ = [("avec", ctypes.POINTER(ctypes.c_double)),
                ("Avec", ctypes.POINTER(ctypes.c_double)),
                ("Bvec", ctypes.POINTER(ctypes.c_double)),
                ("Hvec", ctypes.POINTER(ctypes.c_double)),
                ("f_nu_nr", ctypes.POINTER(ctypes.c_double)),
                ("f_nu_nr_tot", ctypes.POINTER(ctypes.c_double)),
                ("size", ctypes.c_int)]

class MODEL(ctypes.Structure):
    _fields_ = [("h", ctypes.c_double),
                ("Omega_b", ctypes.c_double),
                ("Omega_c", ctypes.c_double),
                ("Omega_k", ctypes.c_double),
                ("N_ur", ctypes.c_double),
                ("N_nu", ctypes.c_int),
                ("M_nu", ctypes.POINTER(ctypes.c_double)),
                ("deg_nu", ctypes.POINTER(ctypes.c_double)),
                ("c_s_nu", ctypes.POINTER(ctypes.c_double)),
                ("T_nu_0", ctypes.c_double),
                ("T_CMB_0", ctypes.c_double),
                ("w0", ctypes.c_double),
                ("wa", ctypes.c_double),
                ("sim_neutrino_nonrel_masses", ctypes.c_int),
                ("sim_neutrino_nonrel_Hubble", ctypes.c_int)]

    # Cosmological tables with functions of time
    tables = TABLES()

    def default_neutrino_arrays(self):
        """
        Initialise arrays with default neutrino parameters if we know the
        number of species, but have not set all the parameters.
        Defaults:
        + M_nu = [0.] * N
        + deg_nu = [1.0] * N
        + c_s_nu = [0.0] * N
        """
        if self.N_nu == 0:
            return

        if not bool(self.M_nu):
            M_nu = [0.] * self.N_nu
            self.M_nu = (ctypes.c_double * self.N_nu)(*M_nu)
        if not bool(self.deg_nu):
            deg_nu = [1.] * self.N_nu
            self.deg_nu = (ctypes.c_double * self.N_nu)(*deg_nu)
        if not bool(self.c_s_nu):
            c_s_nu = [0.] * self.N_nu
            self.c_s_nu = (ctypes.c_double * self.N_nu)(*c_s_nu)

    def set_parameter(self, key, value):
        # Ensure that we know the length when storing neutrino arrays
        array_params = {"M_nu", "deg_nu", "c_s_nu"}
        if key in array_params:
            if self.N_nu == 0 and len(value) > 0:
                raise ValueError("Specify the number of neutrino species (N_nu) before specifying arrays with values per species")
            # Store the address of the array
            store = (ctypes.c_double * self.N_nu)(*value)
        else:
            # Store the value itself
            store = value

        # Store the parameter in the corresponding field
        if key == "h":
            self.h = store
        elif key == "Omega_b":
            self.Omega_b = store
        elif key == "Omega_c":
            self.Omega_c = store
        elif key == "Omega_k":
            self.Omega_k = store
        elif key == "N_ur":
            self.N_ur = store
        elif key == "N_nu":
            self.N_nu = store
        elif key == "M_nu":
            self.M_nu = store
        elif key == "deg_nu":
            self.deg_nu = store
        elif key == "c_s_nu":
            self.c_s_nu = store
        elif key == "T_nu_0":
            self.T_nu_0 = store
        elif key == "T_CMB_0":
            self.T_CMB_0 = store
        elif key == "w0":
            self.w0 = store
        elif key == "wa":
            self.wa = store
        elif key == "sim_neutrino_nonrel_masses":
            self.sim_neutrino_nonrel_masses = store
        elif key == "sim_neutrino_nonrel_Hubble":
            self.sim_neutrino_nonrel_Hubble = store
        else:
            raise KeyError("Unknown parameter")

    def set(self, params):
        for key in params:
            self.set_parameter(key, params[key])
        self.default_neutrino_arrays()
        return True

    def set_sim_type(self, type):
        """
        Set the type of cosmological simulation in terms of how massive
        neutrinos are handled:
        + 0 (fully relativistic)
        + 1 (relativistic Hubble rate, but constant particle masses)
        + 2 (fully non-relativistic)
        Parameters
        ----------
        model : MODEL
            a cosmological model
        type: int
            the type of simulation (0,1,2)
        """
        if not (type == 0 or type == 1 or type == 2):
            raise ValueError("Invalid simulation type")
        self.sim_neutrino_nonrel_masses = 1 if (type == 1 or type == 2) else 0
        self.sim_neutrino_nonrel_Hubble = 1 if type == 2 else 0

    def compute(self, units, physical_consts, a_start, a_final = 1.05, size = 1000):
        integrate_cosmology_tables(self, units, physical_consts, self.tables, a_start, a_final, size)

    def get_H_of_a(self, a):
        return get_H_of_a(self.tables, a)

    def get_f_nu_nr_tot_of_a(self, a):
        return get_f_nu_nr_tot_of_a(self.tables, a)


backend.engine.integrate_cosmology_tables.argtypes = [ctypes.POINTER(MODEL), ctypes.POINTER(units.UNITS), ctypes.POINTER(units.PHYSICAL_CONSTS), ctypes.POINTER(TABLES), ctypes.c_double, ctypes.c_double, ctypes.c_int]

backend.engine.free_cosmology_tables.argtypes = [ctypes.POINTER(TABLES)]

backend.engine.set_neutrino_sound_speeds.argtypes = [ctypes.POINTER(MODEL), ctypes.POINTER(units.UNITS), ctypes.POINTER(units.PHYSICAL_CONSTS)]

backend.engine.get_H_of_a.argtypes = [ctypes.POINTER(TABLES), ctypes.c_double]
backend.engine.get_H_of_a.restype = ctypes.c_double

backend.engine.get_f_nu_nr_tot_of_a.argtypes = [ctypes.POINTER(TABLES), ctypes.c_double]
backend.engine.get_f_nu_nr_tot_of_a.restype = ctypes.c_double

def integrate_cosmology_tables(model, units, physical_consts, tables, a_start, a_final, size):
    """
    Integrate cosmological tables for a given model.
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
    a_start: double
        starting scale factor time of the tables
    a_final: double
        final scale factor time of the tables
    size: int
        length of the tables
    Return
    ------
    bool
        True on success
    """

    backend.engine.integrate_cosmology_tables(ctypes.byref(model), ctypes.byref(units), ctypes.byref(physical_consts), ctypes.byref(tables), a_start, a_final, size)

    return True

def free_cosmology_tables(tables):
    """
    Free the cosmological tables.
    Parameters
    ----------
    tables : TABLES
        the cosmological tables to be freed
    Return
    ------
    bool
        True on success
    """

    backend.engine.free_cosmology_tables(ctypes.byref(tables))

    return True

def set_neutrino_sound_speeds(model, units, physical_consts):
    """
    Set neutrino sound speeds to default values based on the masses.
    Parameters
    ----------
    model : MODEL
        a cosmological model
    units: UNITS
        a unit system
    physical_consts: PHYSICAL_CONSTS
        physical constants corresponding to the unit system
    Return
    ------
    bool
        True on success
    """

    backend.engine.set_neutrino_sound_speeds(ctypes.byref(model), ctypes.byref(units), ctypes.byref(physical_consts))

    return True


def get_H_of_a(tables, a):
    """
    Get the Hubble rate as a function of scale factor a
    Parameters
    ----------
    tables : TABLES
        the cosmological tables to be calculated
    a: double
        the scale factor time
    Return
    ------
    double
        The Hubble rate in internal units
    """

    return backend.engine.get_H_of_a(ctypes.byref(tables), a)

def get_f_nu_nr_tot_of_a(tables, a):
    """
    Get the mass fraction of non-relativistic neutrinos summed over all
    neutrino species, as a function of scale factor a
    Parameters
    ----------
    tables : TABLES
        the cosmological tables to be calculated
    a: double
        the scale factor time
    Return
    ------
    double
        The non-relativistic neutrino fraction
    """

    return backend.engine.get_f_nu_nr_tot_of_a(ctypes.byref(tables), a)
