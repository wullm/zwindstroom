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

#include <gsl/gsl_integration.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

#include "../include/cosmology_tables.h"
#include "../include/units.h"
#include "../include/strooklat.h"


double F_integrand(double x, void *params) {
    double y = *((double*) params);
    double x2 = x * x;
    double y2 = y * y;
    return x2 * sqrt(x2 + y2) / (1.0 + exp(x));
}

double G_integrand(double x, void *params) {
    double y = *((double*) params);
    double x2 = x * x;
    double y2 = y * y;
    return y * x2 / (sqrt(x2 + y2) * (1.0 + exp(x)));
}

double w_tilde(double a, double w0, double wa) {
    return (a - 1.0) * wa - (1.0 + w0 + wa) * log(a);
}

double E2(double a, double Omega_CMB, double Omega_ur, double Omega_nu,
          double Omega_c, double Omega_b, double Omega_lambda, double Omega_k,
          double Omega_dcdm, double Omega_dr, double w0, double wa) {

    const double a_inv = 1.0 / a;
    const double a_inv2 = a_inv * a_inv;
    const double E2 = (Omega_CMB + Omega_ur + Omega_nu + Omega_dr) * a_inv2 * a_inv2 +
                      (Omega_c + Omega_b + Omega_dcdm) * a_inv2 * a_inv +
                      Omega_k * a_inv2 +
                      Omega_lambda * exp(3. * w_tilde(a, w0, wa));
    return E2;
}

struct dcdm_ode_params {
    struct strooklat *spline;
    struct model *m;
    double *Omega_nu_tot;
    double Omega_CMB;
    double Omega_lambda;
    double Omega_ur;
    double H_0;
};

/* The right-hand side of the ODE */
int dcdm_func (double log_a, const double y[], double f[], void *params) {
    struct dcdm_ode_params *p = (struct dcdm_ode_params *) params;
    struct strooklat *spline = p->spline;
    struct model *m = p->m;
    double *Omega_nu_tot = p->Omega_nu_tot;

    /* Pull down constants */
    const double Omega_c = m->Omega_c;
    const double Omega_b = m->Omega_b;
    const double Omega_k = m->Omega_k;
    const double Gamma_dcdm = m->Gamma_dcdm;
    const double w0 = m->w0;
    const double wa = m->wa;
    const double Omega_lambda = p->Omega_lambda;
    const double Omega_CMB = p->Omega_CMB;
    const double Omega_ur = p->Omega_ur;
    const double H_0 = p->H_0;

    /* Cosmological functions of time */
    const double a = exp(log_a);

    double Omega_nu_a = strooklat_interp(spline, Omega_nu_tot, a);
    double Omega_dcdm = y[0] * (a * a * a);
    double Omega_dr = y[1] * (a * a * a * a);

    double E2i = E2(a, Omega_CMB, Omega_ur, Omega_nu_a, Omega_c,
                    Omega_b, Omega_lambda, Omega_k, Omega_dcdm, Omega_dr,
                    w0, wa);
    if (E2i < 0) {
        f[0] = 0;
        f[1] = 0;
        return GSL_FAILURE;
    }
    double H = sqrt(E2i) * H_0;

    /* Now we can write down the ODE */

    /* Decaying dark matter */
    f[0] = -3 * y[0] - Gamma_dcdm / H * y[0];
    /* Dark radiation */
    f[1] = -4 * y[1] + Gamma_dcdm / H * y[0];

    return GSL_SUCCESS;
}

struct dcdm_root_params {
    gsl_odeiv2_driver *d;
    double a_start;
    double Omega_dcdmdr_0;
};

/* Compute the value of Omega_dcdm + Omega_dr at z = 0, for a given value of
 * Omega_dcdm at some initial time. */
double shoot_Omega_dcdm(double Omega_dcdm_ini, void *params) {
    struct dcdm_root_params *p = (struct dcdm_root_params *) params;
    gsl_odeiv2_driver *d = p->d;

    /* Pull down constants */
    const double a_start = p->a_start;
    const double Omega_dcdmdr_0_target = p->Omega_dcdmdr_0;

    /* Allocate the state vector */
    double x[2];

    /* Reset the integration to the start */
    x[0] = Omega_dcdm_ini / (a_start * a_start * a_start);
    x[1] = 0;

    /* Integrate the ODE from the start to this row of the table */
    double log_a = log(a_start);
    double log_a_final = 0; // present day (a=1)

    /* Integrate */
    int status = gsl_odeiv2_driver_apply(d, &log_a, log_a_final, x);

    /* Catch errors */
    if (status != GSL_SUCCESS) {
        printf("ERROR: GSL ODE driver returned error code %d.\n", status);
        return status;
    }

    double Omega_dcdmdr_0 = x[0] + x[1];
    return Omega_dcdmdr_0 - Omega_dcdmdr_0_target;
}

void integrate_cosmology_tables(struct model *m, struct units *us,
                                struct physical_consts *pcs,
                                struct cosmology_tables *tab, double a_start,
                                double a_final, int size) {


    /* Prepare interpolation tables of F(y) and G(y) with y > 0 */
    const int table_size = 500;
    const double y_min = 1e-4;
    const double y_max = 1e6;
    const double log_y_min = log(y_min);
    const double log_y_max = log(y_max);
    const double delta_log_y = (log_y_max - log_y_min) / table_size;

    /* Allocate the tables */
    double *y = malloc(table_size * sizeof(double));
    double *Fy = malloc(table_size * sizeof(double));
    double *Gy = malloc(table_size * sizeof(double));

    /* Prepare GSL integration workspace */
    const int gsl_workspace_size = 1000;
    const double abs_tol = 1e-10;
    const double rel_tol = 1e-10;

    /* Allocate the workspace */
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(gsl_workspace_size);
    gsl_function func_F = {F_integrand};
    gsl_function func_G = {G_integrand};

    /* Perform the numerical integration */
    for (int i=0; i<table_size; i++) {
        y[i] = exp(log_y_min + i * delta_log_y);

        /* Integration result and absolute error */
        double res, abs_err;

        /* Evaluate F(y) by integrating on [0, infinity) */
        func_F.params = &y[i];
        gsl_integration_qagiu(&func_F, 0.0, abs_tol, rel_tol, gsl_workspace_size, workspace, &res, &abs_err);
        Fy[i] = res;

        /* Evaluate G(y) by integrating on [0, infinity) and dividing by F(y) */
        func_G.params = &y[i];
        gsl_integration_qagiu(&func_G, 0.0, abs_tol, rel_tol, gsl_workspace_size, workspace, &res, &abs_err);
        Gy[i] = y[i] * res / Fy[i];
    }

    /* Free the workspace */
    gsl_integration_workspace_free(workspace);

    /* Prepare an interpolation spline for the y argument of F and G */
    struct strooklat spline_y = {y, table_size};
    init_strooklat_spline(&spline_y, 100);

    /* We want to interpolate the scale factor */
    const double a_min = a_start;
    const double a_max = fmax(a_final, 1.01);
    const double log_a_min = log(a_min);
    const double log_a_max = log(a_max);
    const double delta_log_a = (log_a_max - log_a_min) / size;

    tab->size = size;
    tab->avec = malloc(size * sizeof(double));
    tab->Hvec = malloc(size * sizeof(double));
    double *Ga = malloc(size * sizeof(double));
    double *E2a = malloc(size * sizeof(double));
    double *dHdloga = malloc(size * sizeof(double));

    for (int i=0; i<size; i++) {
        tab->avec[i] = exp(log_a_min + i * delta_log_a);
        Ga[i] = sqrt(tab->avec[i] + 1);
    }

    /* Prepare a spline for the scale factor */
    struct strooklat spline = {tab->avec, size};
    init_strooklat_spline(&spline, 100);

    /* The critical density */
    const double h = m->h;
    const double H_0 = h * 100 * KM_METRES / MPC_METRES * us->UnitTimeSeconds;
    const double G_grav = pcs->GravityG;
    const double rho_crit_0 = 3.0 * H_0 * H_0 / (8.0 * M_PI * G_grav);

    /* First, calculate the present-day CMB density from the temperature */
    const double h_bar = pcs->hPlanck / (2.0 * M_PI);
    const double kT = m->T_CMB_0 * pcs->kBoltzmann;
    const double hc = h_bar * pcs->SpeedOfLight;
    const double kT4 = kT * kT * kT * kT;
    const double hc3 = hc * hc * hc;
    const double c2 = pcs->SpeedOfLight * pcs->SpeedOfLight;
    const double Omega_CMB = M_PI * M_PI / 15.0 * (kT4 / hc3) / (rho_crit_0 * c2);

    /* Other density components */
    const double Omega_c = m->Omega_c;
    const double Omega_b = m->Omega_b;
    const double Omega_cb = Omega_c + Omega_b;
    const double Omega_k = m->Omega_k;
    const double Omega_dcdmdr_0 = m->Omega_dcdmdr_0;

    /* Next, calculate the ultra-relativistic density */
    const double ratio = 4. / 11.;
    const double ratio4 = ratio * ratio * ratio * ratio;
    const double Omega_ur = m->N_ur * (7. / 8.) * cbrt(ratio4) * Omega_CMB;

    /* Now, we want to evaluate the neutrino density and equation of state */
    const int N_nu = m->N_nu;
    const double kT_nu_eV_0 = m->T_nu_0 * pcs->kBoltzmann / pcs->ElectronVolt;
    const double T_on_pi = m->T_nu_0 / m->T_CMB_0 / M_PI;
    const double pre_factor = Omega_CMB * 15.0 * T_on_pi * T_on_pi * T_on_pi * T_on_pi;
    double *Omega_nu = malloc(N_nu * size * sizeof(double));
    double *w_nu = malloc(N_nu * size * sizeof(double));

    /* For each neutrino species */
    for (int j=0; j<N_nu; j++) {
        const double M_nu = m->M_nu[j];
        const double deg_nu = m->deg_nu[j];

        /* For each time step, interpolate the distribution function */
        for (int i=0; i<size; i++) {
            /* Compute the density */
            const double arg = tab->avec[i] * M_nu / kT_nu_eV_0;
            const double Farg = strooklat_interp(&spline_y, Fy, arg);
            const double Onu_ij = deg_nu * pre_factor * Farg;
            Omega_nu[j * size + i] = Onu_ij;

            /* Also compute the equation of state */
            const double Garg = strooklat_interp(&spline_y, Gy, arg);
            w_nu[j * size + i] = (1.0 - Garg) / 3.0;
        }

    }

    /* Compute the total neutrino density at each step */
    double *Omega_nu_tot = malloc(size * sizeof(double));

    for (int i=0; i<size; i++) {
        Omega_nu_tot[i] = 0.0;
        for (int j=0; j<N_nu; j++) {
            const double O_nu = Omega_nu[j * size + i];
            Omega_nu_tot[i] += O_nu;
        }
    }

    /* The total neutrino density at z = 0 */
    const double Omega_nu_tot_0 = strooklat_interp(&spline, Omega_nu_tot, 1.0);

    /* The neutrino density per species at z = 0 */
    double *Omega_nu_0 = malloc(N_nu * sizeof(double));
    for (int i = 0; i < N_nu; i++) {
        Omega_nu_0[i] = strooklat_interp(&spline, Omega_nu + i * size, 1.0);
    }

    /* Close the universe */
    const double Omega_lambda = 1.0 - Omega_nu_tot_0 - Omega_k - Omega_ur - Omega_CMB - Omega_c - Omega_b - Omega_dcdmdr_0;
    const double w0 = m->w0;
    const double wa = m->wa;

    /* Solve for the decaying dark matter and dark radiation fractions */
    tab->Omega_dcdm = malloc(size * sizeof(double));
    tab->Omega_dr = malloc(size * sizeof(double));

    /* Do we have decaying dark matter and dark radiation? */
    if (Omega_dcdmdr_0 > 0.) {

        /* GSL ODE integrator */
        struct dcdm_ode_params dcdm_odep;
        gsl_odeiv2_system sys = {dcdm_func, NULL, 2, &dcdm_odep};
        gsl_odeiv2_driver *d;

        dcdm_odep.spline = &spline;
        dcdm_odep.m = m;
        dcdm_odep.Omega_nu_tot = Omega_nu_tot;
        dcdm_odep.Omega_lambda = Omega_lambda;
        dcdm_odep.Omega_CMB = Omega_CMB;
        dcdm_odep.Omega_ur = Omega_ur;
        dcdm_odep.H_0 = H_0;

        double tol = 1e-10;
        double hstart = 1e-10;

        /* Allocate GSL ODE driver */
        d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, hstart, tol, tol);

        /* Find the initial decaying dark matter density */
        double Omega_dcdm_ini = m->Omega_dcdmdr_0;
        /* The value will be in this range */
        double Omega_dcdm_ini_min = 0.0;
        double Omega_dcdm_ini_max = 1.0;

        /* Parameters for root finding */
        struct dcdm_root_params dcdm_root_pars;
        dcdm_root_pars.d = d;
        dcdm_root_pars.a_start = a_start;
        dcdm_root_pars.Omega_dcdmdr_0 = m->Omega_dcdmdr_0;

        /* The function to be solved */
        gsl_function F = {shoot_Omega_dcdm, &dcdm_root_pars};

        /* Set up the GSL root finder */
        const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
        gsl_root_fsolver *s = gsl_root_fsolver_alloc (T);
        gsl_root_fsolver_set (s, &F, Omega_dcdm_ini_min, Omega_dcdm_ini_max);

        int status;
        int iter = 0, max_iter = 100;
        while (iter == 0 || (status == GSL_CONTINUE && iter < max_iter)) {
            iter++;
            status = gsl_root_fsolver_iterate(s);
            Omega_dcdm_ini = gsl_root_fsolver_root(s);
            Omega_dcdm_ini_min = gsl_root_fsolver_x_lower(s);
            Omega_dcdm_ini_max = gsl_root_fsolver_x_upper(s);
            status = gsl_root_test_interval(Omega_dcdm_ini_min, Omega_dcdm_ini_max, 0, tol);
        }

        if (status != GSL_SUCCESS) {
            printf("ERROR: GSL root finder returned error code %d.\n", status);
            return;
        }

        /* Clean up the root finder */
        gsl_root_fsolver_free (s);

        /* Allocate the state vector */
        double x[2];

        for (int i=0; i<size; i++) {
            /* Reset the integration to the start */
            x[0] = Omega_dcdm_ini / (a_start * a_start * a_start);
            x[1] = 0;

            /* Integrate the ODE from the start to this row of the table */
            double log_a = log(a_start);
            double log_a_final = log(tab->avec[i]);

            /* Integrate */
            int status = gsl_odeiv2_driver_apply(d, &log_a, log_a_final, x);

            /* Catch errors */
            if (status != GSL_SUCCESS) {
                printf("ERROR: GSL ODE driver returned error code %d.\n", status);
                return;
            }

            /* Store the decaying dark matter and dark radiation densities */
            tab->Omega_dcdm[i] = x[0] * (tab->avec[i] * tab->avec[i] * tab->avec[i]);
            tab->Omega_dr[i] = x[1] * (tab->avec[i] * tab->avec[i] * tab->avec[i] * tab->avec[i]);
        }

        /* Free the GSL ODE drive */
        gsl_odeiv2_driver_free(d);
    } else {
        /* Otherwise, the densities are always zero */
        for (int i=0; i<size; i++) {
            tab->Omega_dcdm[i] = 0.;
            tab->Omega_dr[i] = 0.;
        }
    }

    /* Split the neutrino densities into relativistic and non-relativistic parts */
    double *Omega_nu_nr = malloc(size * sizeof(double));
    double *Omega_r = malloc(size * sizeof(double));
    double *Omega_m = malloc(size * sizeof(double));
    tab->f_nu_nr = malloc(size * N_nu * sizeof(double));
    tab->f_nu_nr_tot = malloc(size * sizeof(double));

    for (int i=0; i<size; i++) {

        /* Start with non-neutrino contributions to radiation & matter */
        Omega_r[i] = Omega_CMB + Omega_ur + tab->Omega_dr[i];
        Omega_m[i] = Omega_c + Omega_b + tab->Omega_dcdm[i];
        Omega_nu_nr[i] = 0.0;

        /* Add the massive neutrino species */
        for (int j=0; j<N_nu; j++) {
            const double O_nu = Omega_nu[j * size + i];
            const double w = w_nu[j * size + i];
            Omega_nu_nr[i] += (1.0 - 3.0 * w) * O_nu;
            Omega_r[i] += 3.0 * w * O_nu;
            /* We rescale by 1/a, since this is in fact Omega_m * E^2 * a^3 and
             * Omega_nu is in fact Omega_nu * E^2 * a^4 */
            Omega_m[i] += (1.0 - 3.0 * w) * O_nu / tab->avec[i];
        }

        /* Fraction of non-relativistic neutrinos in matter */
        tab->f_nu_nr_tot[i] = Omega_nu_nr[i] / Omega_m[i] / tab->avec[i];

        /* Fraction per species */
        for (int j=0; j<N_nu; j++) {
            const double O_nu = Omega_nu[j * size + i];
            const double w = w_nu[j * size + i];
            const double O_nu_nr = (1.0 - 3.0 * w) * O_nu;
            tab->f_nu_nr[j * size + i] = O_nu_nr / Omega_m[i] / tab->avec[i];
        }
    }

    /* If neutrinos are to be treated as non-relativistic for the purpose of
     * calculating the Hubble rate, we need to replace the previous calculation
     * of Omega_nu(a) with Omega_nu_0 now (i.e. before calculating Hvec). */
    if (m->sim_neutrino_nonrel_Hubble) {
        for (int i=0; i<size; i++) {
            Omega_m[i] = Omega_cb + Omega_nu_tot_0 + tab->Omega_dcdm[i];
            Omega_nu_tot[i] = Omega_nu_tot_0;
            Omega_nu_nr[i] = Omega_nu_tot_0;
            tab->f_nu_nr_tot[i] = Omega_nu_tot_0 / (Omega_cb + Omega_nu_tot_0 + tab->Omega_dcdm[i]);
            for (int j=0; j<N_nu; j++) {
                tab->f_nu_nr[j * size + i] = Omega_nu_0[j] / (Omega_cb + Omega_nu_tot_0 + tab->Omega_dcdm[i]);
            }
        }
    }

    /* Now, create a table with the Hubble rate */
    for (int i=0; i<size; i++) {
        double Omega_nu_a = strooklat_interp(&spline, Omega_nu_tot, tab->avec[i]);
        double Omega_dcdm = tab->Omega_dcdm[i];
        double Omega_dr = tab->Omega_dr[i];
        E2a[i] = E2(tab->avec[i], Omega_CMB, Omega_ur, Omega_nu_a, Omega_c,
                       Omega_b, Omega_lambda, Omega_k, Omega_dcdm, Omega_dr,
                       w0, wa);
        tab->Hvec[i] = sqrt(E2a[i]) * H_0;
    }

    /* Now, differentiate the Hubble rate */
    for (int i=0; i<size; i++) {
        /* Forward at the start, five-point in the middle, backward at the end */
        if (i < 2) {
            dHdloga[i] = (tab->Hvec[i+1] - tab->Hvec[i]) / delta_log_a;
        } else if (i < size - 2) {
            dHdloga[i]  = tab->Hvec[i-2];
            dHdloga[i] -= tab->Hvec[i-1] * 8.0;
            dHdloga[i] += tab->Hvec[i+1] * 8.0;
            dHdloga[i] -= tab->Hvec[i+2];
            dHdloga[i] /= 12.0 * delta_log_a;
        } else {
            dHdloga[i] = (tab->Hvec[i] - tab->Hvec[i-1]) / delta_log_a;
        }
    }

    /* If neutrino particle masses are not varied in the N-body simulation to
     * account for the relativistic energy density, we need to replace the
     * previous calculation of Omega_nu(a) with Omega_nu_0. However, this
     * must be done after calculating the Hubble rate, if we do take the
     * relativistic contribution into account there. */
    if (m->sim_neutrino_nonrel_masses && !m->sim_neutrino_nonrel_Hubble) {
        for (int i=0; i<size; i++) {
            Omega_m[i] = Omega_cb + Omega_nu_tot_0 + tab->Omega_dcdm[i];
            Omega_nu_tot[i] = Omega_nu_tot_0;
            Omega_nu_nr[i] = Omega_nu_tot_0;
            tab->f_nu_nr_tot[i] = Omega_nu_tot_0 / (Omega_cb + Omega_nu_tot_0 + tab->Omega_dcdm[i]);
            for (int j=0; j<N_nu; j++) {
                tab->f_nu_nr[j * size + i] = Omega_nu_0[j] / (Omega_cb + Omega_nu_tot_0 + tab->Omega_dcdm[i]);
            }
        }
    }

    /* Now, create the A and B functions */
    tab->Avec = malloc(size * sizeof(double));
    tab->Bvec = malloc(size * sizeof(double));

    for (int i=0; i<size; i++) {
        double a = tab->avec[i];
        tab->Avec[i] = -(2.0 + dHdloga[i] / tab->Hvec[i]);
        tab->Bvec[i] = -1.5 * Omega_m[i] / (a * a * a) / E2a[i];
    }

    free(Omega_nu_nr);
    free(Omega_r);
    free(Omega_m);
    free(Omega_nu);
    free(Omega_nu_tot);
    free(Omega_nu_0);
    free(w_nu);
    free(dHdloga);
    free(E2a);
    free(Ga);

    /* Free the interpolation tables */
    free(y);
    free(Fy);
    free(Gy);

    free_strooklat_spline(&spline);
    free_strooklat_spline(&spline_y);
}

double get_H_of_a(struct cosmology_tables *tab, double a) {
    /* Prepare a spline for the scale factor */
    struct strooklat spline = {tab->avec, tab->size};
    init_strooklat_spline(&spline, 100);

    /* Interpolate */
    double Ha = strooklat_interp(&spline, tab->Hvec, a);

    /* Free the spline */
    free_strooklat_spline(&spline);

    return Ha;
}

double get_f_nu_nr_tot_of_a(struct cosmology_tables *tab, double a) {
    /* Prepare a spline for the scale factor */
    struct strooklat spline = {tab->avec, tab->size};
    init_strooklat_spline(&spline, 100);

    /* Interpolate */
    double f_nu_nr_tot = strooklat_interp(&spline, tab->f_nu_nr_tot, a);

    /* Free the spline */
    free_strooklat_spline(&spline);

    return f_nu_nr_tot;
}

double get_Omega_dcdm_of_a(struct cosmology_tables *tab, double a) {
    /* Prepare a spline for the scale factor */
    struct strooklat spline = {tab->avec, tab->size};
    init_strooklat_spline(&spline, 100);

    /* Interpolate */
    double Omega_dcdm = strooklat_interp(&spline, tab->Omega_dcdm, a);

    /* Free the spline */
    free_strooklat_spline(&spline);

    return Omega_dcdm;
}

double get_Omega_dr_of_a(struct cosmology_tables *tab, double a) {
    /* Prepare a spline for the scale factor */
    struct strooklat spline = {tab->avec, tab->size};
    init_strooklat_spline(&spline, 100);

    /* Interpolate */
    double Omega_dr = strooklat_interp(&spline, tab->Omega_dr, a);

    /* Free the spline */
    free_strooklat_spline(&spline);

    return Omega_dr;
}

void set_neutrino_sound_speeds(struct model *m, struct units *us,
                               struct physical_consts *pcs) {

    /* Use the estimate from Blas+14 as default */
    for (int i = 0; i < m->N_nu; i++) {
        m->c_s_nu[i] = pcs->SoundSpeedNeutrinos / m->M_nu[i];
    }

}

void free_cosmology_tables(struct cosmology_tables *tab) {
    free(tab->avec);
    free(tab->Avec);
    free(tab->Bvec);
    free(tab->Hvec);
    free(tab->f_nu_nr);
    free(tab->f_nu_nr_tot);
    free(tab->Omega_dcdm);
    free(tab->Omega_dr);
}
