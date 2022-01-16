/*******************************************************************************
 * This file is part of Zwindstroom.
 * Copyright (c) 2021 Willem Elbers (whe@willemelbers.com)
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <sys/time.h>

#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>

#include "../include/fluid_equations.h"
#include "../include/strooklat.h"

struct ode_params {
    struct strooklat *spline;
    struct cosmology_tables *tab;
    double k;
    double f_b;
    double *c_s_nu;
    int N_nu;
};

int func (double log_a, const double y[], double f[], void *params) {
    struct ode_params *p = (struct ode_params *) params;
    struct strooklat *spline = p->spline;
    struct cosmology_tables *tab = p->tab;

    /* Pull down some constants */
    const int N_nu = p->N_nu;
    const double k = p->k;
    const double f_b = p->f_b;

    /* Cosmological functions of time */
    const double a = exp(log_a);
    const double A = strooklat_interp(spline, tab->Avec, a);
    const double B = strooklat_interp(spline, tab->Bvec, a);
    const double H = strooklat_interp(spline, tab->Hvec, a);

    /* Non-relativistic neutrino density fractions */
    double f_nu_nr_tot = 0.;
    double f_nu_nr[N_nu];
    for (int i = 0; i < N_nu; i++) {
        f_nu_nr[i] = strooklat_interp(spline, tab->f_nu_nr + i * tab->size, a);
        f_nu_nr_tot += f_nu_nr[i];
    }

    /* Neutrino free-streaming scales (squared) */
    double k_fs2[N_nu];
    for (int i = 0; i < N_nu; i++) {
        double c_s = p->c_s_nu[i] / a;
        k_fs2[i] = -B * H * H / (c_s * c_s) * (a * a);
    }

    /* The weighted cdm and baryon density perturbation */
    double D_cb = (1.0 - f_b) * y[0] + f_b * y[2];

    /* The weighted total neutrino density perturbation */
    double D_nu_tot = 0.;
    for (int i = 0; i < N_nu; i++) {
        D_nu_tot += y[4 + 2 * i] * f_nu_nr[i] / f_nu_nr_tot;
    }

    /* Now we can write down the ODE */

    /* CDM */
    f[0] = -y[1];
    f[1] = A * y[1] + B * ((1.0 - f_nu_nr_tot) * D_cb + f_nu_nr_tot * D_nu_tot);
    /* Baryons */
    f[2] = -y[3];
    f[3] =  A * y[3] + B * ((1.0 - f_nu_nr_tot) * D_cb + f_nu_nr_tot * D_nu_tot);
    /* Neutrinos */
    for (int i = 0; i < N_nu; i++) {
        f[4 + 2 * i] = -y[5 + 2 * i];
        f[5 + 2 * i] = A * y[5 + 2 * i] + B * ((1.0 - f_nu_nr_tot) * D_cb + f_nu_nr_tot * D_nu_tot - (k * k) / k_fs2[i] * y[4 + 2 * i]);
    }

    return GSL_SUCCESS;
}

/* GSL ODE integrator */
struct ode_params odep;
gsl_odeiv2_system sys = {func, NULL, 0, &odep};
gsl_odeiv2_driver *d;

struct strooklat spline_cosmo;

void prepare_fluid_integrator(struct model *m, struct units *us,
                              struct physical_consts *pcs,
                              struct cosmology_tables *tab, double tol,
                              double hstart) {

    /* Set the dimension of the system (2nd order ODE for cdm, baryons, and
     * N neutrino species makes for 4 + 2N degrees of greedom) */
    sys.dimension = 4 + 2 * m->N_nu;

    /* Prepare the parameters for the fluid ODEs */
    odep.spline = &spline_cosmo;
    odep.tab = tab;
    odep.f_b = m->Omega_b / (m->Omega_c + m->Omega_b);
    odep.N_nu = m->N_nu;
    odep.c_s_nu = m->c_s_nu;

    /* Check if the sounds speeds have been set */
    int need_to_set_speeds = 0;
    for (int i = 0; i < m->N_nu; i++) {
        need_to_set_speeds += m->c_s_nu[i] == 0.;
    }

    /* Set sounds speeds equal to default values based on the masses */
    if (need_to_set_speeds) {
        set_neutrino_sound_speeds(m, us, pcs);
    }

    /* Allocate GSL ODE driver */
    d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, hstart, tol, tol);

    /* Prepare a spline for the cosmological tables */
    spline_cosmo.x = tab->avec;
    spline_cosmo.size = tab->size;
    init_strooklat_spline(&spline_cosmo, 100);
}

void integrate_fluid_equations(struct model *m, struct units *us,
                               struct physical_consts *pcs,
                               struct cosmology_tables *tab,
                               struct growth_factors *gfac,
                               double a_start, double a_final) {

    /* The wavenumber of interest */
    odep.k = gfac->k;

    /* The number of neutrino species */
    const int N_nu = m->N_nu;

    /* Dimension of the problem (2nd order ODE for cdm, baryons, and N
     * neutrino species makes for 4 + 2N degrees of greedom) */
    const int dim = 4 + 2 * N_nu;

    /* The initial conditions */
    double Dc_ini = gfac->delta_c;
    double Db_ini = gfac->delta_b;
    double *Dn_ini = gfac->delta_n;

    /* Growth rates at a_start */
    double gc = gfac->gc;
    double gb = gfac->gb;
    double *gn = gfac->gn;

    /* Prepare the state variables */
    double *y = malloc(dim * sizeof(double));

    /* CDM */
    y[0] = Dc_ini;
    y[1] = -gc * Dc_ini;
    /* Baryons */
    y[2] = Db_ini;
    y[3] = -gb * Db_ini;
    /* Neutrinos */
    for (int i = 0; i < m->N_nu; i++) {
        y[4 + 2 * i] = Dn_ini[i];
        y[5 + 2 * i] = -gn[i] * Dn_ini[i];
    }

    /* We will integrate from log(a) = log(a_start) to log(a_final) */
    double log_a = log(a_start);
    double log_a_final = log(a_final);

    /* Integrate */
    gsl_odeiv2_driver_apply(d, &log_a, log_a_final, y);

    /* Extract the final densities */
    double Dc_final = y[0];
    double Db_final = y[2];
    double Dn_final[N_nu];
    for (int i = 0; i < N_nu; i++) {
        Dn_final[i] = y[4 + 2 * i];
    }

    /* Store the relative growth factors between a_start and a_final */
    gfac->Dc = Dc_ini / Dc_final;
    gfac->Db = Db_ini / Db_final;
    for (int i = 0; i < N_nu; i++) {
        gfac->Dn[i] = Dn_ini[i] / Dn_final[i];
    }

    /* Free the state variables */
    free(y);
}

void free_fluid_integrator() {
    /* Free the GSL ODE drive */
    gsl_odeiv2_driver_free(d);

    /* Free the spline */
    free_strooklat_spline(&spline_cosmo);
}
