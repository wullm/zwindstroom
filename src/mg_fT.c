/*******************************************************************************
 * This file is part of Zwindstroom.
 * Copyright (c) 2024 Willem Elbers (whe@willemelbers.com)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf.h>

#include "../include/mg_fT.h"
#include "../include/strooklat.h"

/* Right-hand size of dH/da for fT gravity */
int fT_dfunc(double a, const double y[], double f[], void *params_ptr) {
    struct fT_df_params *lparams = (struct fT_df_params*) params_ptr;
    struct strooklat *spline = lparams->spline;
    struct model *m = lparams->m;
    struct units *us = lparams->us;
    struct physical_consts *pcs = lparams->pcs;
    const int size = lparams->size;

    /* Pull down some constants */
    const double Omega_CMB = lparams->Omega_CMB;
    const double Omega_ur = lparams->Omega_ur;
    const double Omega_c = lparams->Omega_c;
    const double Omega_b = lparams->Omega_b;
    const double h = m->h;
    const double b = m->b;

    /* Vector with neutrino densities (per species) */
    double *Omega_nu = lparams->Omega_nu;
    /* Vector with neutrino equations of state (per species) */
    double *w_nu = lparams->w_nu;

    /* Define further constants */
    const double H_0 = h * 100.0 * KM_METRES / MPC_METRES * us->UnitTimeSeconds;
    const double H_0_2 = H_0 * H_0;
    const double inv_a = 1. / a;
    const double inv_a2 = inv_a * inv_a;
    const double inv_a4 = inv_a2 * inv_a2;
    
    /* Density and pressure from photons and ultra-relativistic particles */
    const double rho_r = 3.0 * H_0_2 * (Omega_CMB + Omega_ur) * inv_a4;
    const double p_r = rho_r / 3.0;
    /* Density and pressure from CDM and baryons (no neutrinos here) */
    const double rho_m = 3.0 * H_0_2 * (Omega_c + Omega_b) * inv_a2 * inv_a;
    const double p_m = 0.;
    /* Density and pressure from all massive neutrino species */
    double rho_nu = 0.;
    double p_nu = 0.;
    double Omega_nu_tot = 0.;
    for (int j=0; j<m->N_nu; j++) {
        const double O_nu = strooklat_interp(spline, Omega_nu + j * size, a);
        const double w = strooklat_interp(spline, w_nu + j * size, a);
        const double rho = 3.0 * H_0_2 * O_nu * inv_a4;
        rho_nu += rho;
        p_nu += w * rho;
        Omega_nu_tot += O_nu;
    }

    /* Hence, total pressure and ratiation */
    const double p_tot = p_r + p_m + p_nu;
    const double rho_tot = rho_r + rho_m + rho_nu;
    
    /* Density parameters */
    const double Omega_m = Omega_c + Omega_b;
    const double Omega_r = Omega_CMB + Omega_ur + Omega_nu_tot;

    /* Intermediate steps */
    const double T0 = -6.0 * H_0_2;

    /* Intermediate steps */
    const double half_on_b = 0.5 / b;
    const double alpha = gsl_sf_lambert_W0(-(Omega_m + Omega_r) * 0.5 * exp(-0.5 / b) / b) + half_on_b;
    const double T = -6.0 * y[0] * y[0];
    const double T0_on_T_pow_b = pow(T0 / T, b);
    const double exp_T0_on_T_pow_b_alpha = exp(T0_on_T_pow_b * alpha);
    const double FT = -1.0 + exp_T0_on_T_pow_b_alpha * (1.0 - b * T0_on_T_pow_b * alpha);
    const double FTT = (b * exp_T0_on_T_pow_b_alpha * T0_on_T_pow_b * alpha * (-1.0 + b + b * T0_on_T_pow_b * alpha)) / T;

    /* Final answer, dH/da */
    f[0] = (-0.5 * (rho_tot + p_tot) / (1.0 + FT + 2.0 * T * FTT)) / (a * y[0]);

    return GSL_SUCCESS;
}
