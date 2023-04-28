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
          double w0, double wa) {

    const double a_inv = 1.0 / a;
    const double a_inv2 = a_inv * a_inv;
    const double E2 = (Omega_CMB + Omega_ur + Omega_nu) * a_inv2 * a_inv2 +
                      (Omega_c + Omega_b) * a_inv2 * a_inv +
                      Omega_k * a_inv2 +
                      Omega_lambda * exp(3. * w_tilde(a, w0, wa));
    return E2;
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
    tab->Gvec = malloc(size * sizeof(double));
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

    /* Split the neutrino densities into relativistic and non-relativistic parts */
    double *Omega_nu_nr = malloc(size * sizeof(double));
    double *Omega_nu_tot = malloc(size * sizeof(double));
    double *Omega_r = malloc(size * sizeof(double));
    double *Omega_m = malloc(size * sizeof(double));
    tab->f_nu_nr = malloc(size * N_nu * sizeof(double));
    tab->f_nu_nr_tot = malloc(size * sizeof(double));

    for (int i=0; i<size; i++) {

        /* Start with constant contributions to radiation & matter */
        Omega_r[i] = Omega_CMB + Omega_ur;
        Omega_m[i] = Omega_c + Omega_b;
        Omega_nu_nr[i] = 0.0;
        Omega_nu_tot[i] = 0.0;

        /* Add the massive neutrino species */
        for (int j=0; j<N_nu; j++) {
            const double O_nu = Omega_nu[j * size + i];
            const double w = w_nu[j * size + i];
            Omega_nu_tot[i] += O_nu;
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

    /* The total neutrino density at z = 0 */
    const double Omega_nu_tot_0 = strooklat_interp(&spline, Omega_nu_tot, 1.0);

    /* The neutrino density per species at z = 0 */
    double *Omega_nu_0 = malloc(N_nu * sizeof(double));
    for (int i = 0; i < N_nu; i++) {
        Omega_nu_0[i] = strooklat_interp(&spline, Omega_nu + i * size, 1.0);
    }

    /* Close the universe */
    const double Omega_lambda = 1.0 - Omega_nu_tot_0 - Omega_k - Omega_ur - Omega_CMB - Omega_c - Omega_b;
    const double w0 = m->w0;
    const double wa = m->wa;

    /* If neutrinos are to be treated as non-relativistic for the purpose of
     * calculating the Hubble rate, we need to replace the previous calculation
     * of Omega_nu(a) with Omega_nu_0 now (i.e. before calculating Hvec). */
    if (m->sim_neutrino_nonrel_Hubble) {
        for (int i=0; i<size; i++) {
            Omega_m[i] = Omega_cb + Omega_nu_tot_0;
            Omega_nu_tot[i] = Omega_nu_tot_0;
            Omega_nu_nr[i] = Omega_nu_tot_0;
            tab->f_nu_nr_tot[i] = Omega_nu_tot_0 / (Omega_cb + Omega_nu_tot_0);
            for (int j=0; j<N_nu; j++) {
                tab->f_nu_nr[j * size + i] = Omega_nu_0[j] / (Omega_cb + Omega_nu_tot_0);
            }
        }
    }

    /* Now, create a table with the Hubble rate */
    if (m->has_external_E) {
        /* Read the external Hubble rate file */
        char *line = NULL;
        size_t len = 0;
        ssize_t read;

        FILE *f = fopen(m->path_external_E, "r");
        if (f == NULL) exit(EXIT_FAILURE);

        /* Count the number of non-comment lines */
        int lines = 0;
        while ((read = getline(&line, &len, f)) != -1) {
            /* Skip comments */
            if (line[0] != '#') lines++;
        }
        rewind(f);

        /* Allocate memory for the input data */
        double *ext_avec = malloc(lines * sizeof(double));
        double *ext_Evec = malloc(lines * sizeof(double));

        int ctr = 0;
        while ((read = getline(&line, &len, f)) != -1) {
            /* Skip comments */
            if (line[0] == '#') continue;
            /* Parse scale factor and dimensionless Hubble rate (= E(a)) */
            float a, E;
            sscanf(line, "%f %f", &a, &E);
            ext_avec[ctr] = a;
            ext_Evec[ctr] = E;
            ctr++;
        }

        /* Prepare an interpolation spline */
        struct strooklat ext_spline = {ext_avec, lines};
        init_strooklat_spline(&ext_spline, 100);

        /* Interpolate and store the result in the internal Hubble vector */
        for (int i=0; i<size; i++) {
            double E = strooklat_interp(&ext_spline, ext_Evec, tab->avec[i]);
            E2a[i] = E * E;
            tab->Hvec[i] = E * H_0;
        }

        /* Free the interpolation spline */
        free_strooklat_spline(&ext_spline);

        /* Free the memory and close the file */
        free(ext_avec);
        free(ext_Evec);
        fclose(f);
    } else {
        /* Compute the Hubble rate */
        for (int i=0; i<size; i++) {
            double Omega_nu_a = strooklat_interp(&spline, Omega_nu_tot, tab->avec[i]);
            E2a[i] = E2(tab->avec[i], Omega_CMB, Omega_ur, Omega_nu_a, Omega_c,
                           Omega_b, Omega_lambda, Omega_k, w0, wa);
            tab->Hvec[i] = sqrt(E2a[i]) * H_0;
        }
    }

    /* Now, create a table with the effective Gravitational constant */
    if (m->has_external_G_eff) {
        /* Read the external Gravitational constant file */
        char *line = NULL;
        size_t len = 0;
        ssize_t read;

        FILE *f = fopen(m->path_external_G_eff, "r");
        if (f == NULL) exit(EXIT_FAILURE);

        /* Count the number of non-comment lines */
        int lines = 0;
        while ((read = getline(&line, &len, f)) != -1) {
            /* Skip comments */
            if (line[0] != '#') lines++;
        }
        rewind(f);

        /* Allocate memory for the input data */
        double *ext_avec = malloc(lines * sizeof(double));
        double *ext_Gvec = malloc(lines * sizeof(double));

        int ctr = 0;
        while ((read = getline(&line, &len, f)) != -1) {
            /* Skip comments */
            if (line[0] == '#') continue;
            /* Parse scale factor and G_eff */
            float a, G;
            sscanf(line, "%f %f", &a, &G);
            ext_avec[ctr] = a;
            ext_Gvec[ctr] = G;
            ctr++;
        }

        /* Prepare an interpolation spline */
        struct strooklat ext_spline = {ext_avec, lines};
        init_strooklat_spline(&ext_spline, 100);

        /* Interpolate and store the result in the internal G_eff vector */
        for (int i=0; i<size; i++) {
            double G_factor = strooklat_interp(&ext_spline, ext_Gvec, tab->avec[i]);
            tab->Gvec[i] = G_factor * G_grav;
        }

        /* Free the interpolation spline */
        free_strooklat_spline(&ext_spline);

        /* Free the memory and close the file */
        free(ext_avec);
        free(ext_Gvec);
        fclose(f);
    } else {
        /* The Newtonian gravitational constant is just constant */
        for (int i=0; i<size; i++) {
            tab->Gvec[i] = G_grav;
        }
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
            Omega_m[i] = Omega_cb + Omega_nu_tot_0;
            Omega_nu_tot[i] = Omega_nu_tot_0;
            Omega_nu_nr[i] = Omega_nu_tot_0;
            tab->f_nu_nr_tot[i] = Omega_nu_tot_0 / (Omega_cb + Omega_nu_tot_0);
            for (int j=0; j<N_nu; j++) {
                tab->f_nu_nr[j * size + i] = Omega_nu_0[j] / (Omega_cb + Omega_nu_tot_0);
            }
        }
    }

    /* Now, create the A and B functions */
    tab->Avec = malloc(size * sizeof(double));
    tab->Bvec = malloc(size * sizeof(double));

    for (int i=0; i<size; i++) {
        double a = tab->avec[i];
        tab->Avec[i] = -(2.0 + dHdloga[i] / tab->Hvec[i]);
        tab->Bvec[i] = -1.5 * Omega_m[i] / (a * a * a) / E2a[i] * (tab->Gvec[i] / G_grav);
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
    free(tab->Gvec);
    free(tab->f_nu_nr);
    free(tab->f_nu_nr_tot);
}
