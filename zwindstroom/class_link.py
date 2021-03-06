"""
Contains methods to run CLASS and extract relevant quantities.
"""

from classy import Class
import zwindstroom.units as units
from math import log, exp

def get_class_name(name):
    dictionary = {"Omega_c": "Omega_cdm",
                  "N_nu": "N_ncdm",
                  "M_nu": "m_ncdm",
                  "deg_nu": "deg_ncdm",
                  "T_nu_0": "T_ncdm",
                  "T_CMB_0": "T_cmb",
                  "w0": "w0_fld",
                  "wa": "wa_fld"}

    return dictionary[name] if name in dictionary else name

def double_array_to_string(array):
    return "".join("%g, " % d for d in array)[:-2]

def run_class(model, unit_system, a_first, k_max, extra_params = {}):
    """
    Run CLASS on a cosmological model.
    Parameters
    ----------
    model : MODEL
        a cosmological model
    units: UNITS
        a unit system
    a_first: double
        starting scale factor time of the CLASS tables
    k_max: double
        maximum wavenumber in internal units to compute transfer functions
    extra_params: dictionary
        optionally pass along additional parameters to CLASS
    Return
    ------
    cosmo
        A CLASS object
    """
    model_params = {}
    model_params["h"] = model.h
    model_params["Omega_b"] = model.Omega_b
    model_params["Omega_c"] = model.Omega_c
    model_params["Omega_k"] = model.Omega_k
    model_params["N_ur"] = model.N_ur
    model_params["N_nu"] = model.N_nu
    model_params["M_nu"] = double_array_to_string([model.M_nu[i] for i in range(model.N_nu)])
    model_params["deg_nu"] = double_array_to_string([model.deg_nu[i] for i in range(model.N_nu)])
    model_params["T_nu_0"] = double_array_to_string([model.T_nu_0 / model.T_CMB_0 for i in range(model.N_nu)])
    model_params["T_CMB_0"] = model.T_CMB_0
    model_params["w0"] = model.w0
    model_params["wa"] = model.wa

    class_params = {}
    for key in model_params:
        class_params[get_class_name(key)] = model_params[key]

    # Use fld parametrisation instead
    class_params["Omega_Lambda"] = 0.0

    # Make sure that CLASS outputs power spectra and transfer functions
    class_params["output"] = "mPk dTk vTk"

    # Specify the maximum redshift
    z_max = 1.0 / a_first - 1.0
    class_params["z_max_pk"] = z_max

    # Specify the maximum wavenumber
    k_max_Mpc = k_max * units.__Mpc / unit_system.UnitLengthMetres
    class_params["P_k_max_1/Mpc"] = k_max_Mpc

    # Add the additional parameters
    class_params = {**class_params, **extra_params}

    # Run CLASS
    cosmo = Class()
    cosmo.set(class_params)
    cosmo.compute()
    return cosmo

def get_growth_rates(model, cosmo, a):
    """
    Extract logarithmic growth rates from CLASS at a given scale factor a
    Parameters
    ----------
    model : MODEL
        a cosmological model
    cosmo: CLASS
        a CLASS object
    a: double
        the scale factor
    Return
    ------
    growth_rates
        Dictionary with growth rates for the fluid species
    """
    # Prepare the output
    growth_rates = {}
    keys = ["d_cdm", "d_b"]
    keys += ["d_ncdm[%d]" % i for i in range(model.N_nu)]

    # We use a central difference approximation
    loga = log(a)
    dloga = 0.001
    a_prev = exp(loga - 0.5 * dloga)
    a_next = min(exp(loga + 0.5 * dloga), 1)
    z_prev = 1.0 / a_prev - 1.0
    z_next = 1.0 / a_next - 1.0
    z = 1.0 / a - 1.0

    transfer = cosmo.get_transfer(z = z)
    transfer_prev = cosmo.get_transfer(z = z_prev)
    transfer_next = cosmo.get_transfer(z = z_next)

    for key in keys:
        d_next = transfer_next[key]
        d_prev = transfer_prev[key]
        d = transfer[key]
        dlogd_dloga = (d_next - d_prev) / (d * dloga)
        growth_rates[key] = dlogd_dloga

    return growth_rates

def get_growth_factors(model, cosmo, a):
    """
    Extract growth factors (density transfer functions) from CLASS at a given
    scale factor a
    Parameters
    ----------
    model : MODEL
        a cosmological model
    cosmo: CLASS
        a CLASS object
    a: double
        the scale factor
    Return
    ------
    growth_factors
        Dictionary with growth factors for the fluid species
    """
    # Prepare the output
    growth_factors = {}
    keys = ["d_cdm", "d_b"]
    keys += ["d_ncdm[%d]" % i for i in range(model.N_nu)]


    # Extract the transfer functions as redshift z
    z = 1.0 / a - 1.0
    transfer = cosmo.get_transfer(z = z)

    # Normalise by the cdm density transfer function
    d_cdm = transfer["d_cdm"]

    for key in keys:
        d = transfer[key] / d_cdm
        growth_factors[key] = d

    return growth_factors

def get_wavenumbers(model, cosmo, unit_system):
    """
    Extract vector of wavenumbers from CLASS in internal units
    Parameters
    ----------
    model : MODEL
        a cosmological model
    cosmo: CLASS
        a CLASS object
    unit_system: UNITS
        the unit system
    Return
    ------
    wavenumbers: vector of doubles
        The wavenumbers
    """
    transfer = cosmo.get_transfer()
    return transfer["k (h/Mpc)"] * model.h * units.__Mpc / unit_system.UnitLengthMetres
