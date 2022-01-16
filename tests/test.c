/*******************************************************************************
 * This file is part of Zwindstroom.
 * Copyright (c) 2020 Willem Elbers (whe@willemelbers.com)
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <zwindstroom.h>

int main() {

    /* 3FA structs */
    struct model m;
    struct units us;
    struct physical_consts pcs;
    struct cosmology_tables tab;
    struct growth_factors gfac;

    /* Choice of neutrino masses */
    double M_nu[1] = {0.067666666};
    double deg_nu[1] = {3.0};
    double c_s_nu[1];

    /* Specify the cosmological model */
    m.h = 0.6771;
    m.Omega_b = 0.0495;
    m.Omega_c = 0.2491464152;
    m.Omega_k = 0.;
    m.N_ur = 0.00441;
    m.N_nu = 1;
    m.M_nu = M_nu;
    m.deg_nu = deg_nu;
    m.c_s_nu = c_s_nu;
    m.T_nu_0 = 1.951757805;
    m.T_CMB_0 = 2.7255;
    m.w0 = -1.0;
    m.wa = 0.0;
    m.sim_neutrino_nonrel_masses = 1;

    /* Set up 3FA unit system */
    us.UnitLengthMetres = MPC_METRES;
    us.UnitTimeSeconds = 1.0;
    us.UnitMassKilogram = 1.0;
    us.UnitTemperatureKelvin = 1.0;
    us.UnitCurrentAmpere = 1.0;
    set_physical_constants(&us, &pcs);

    /* Set default neutrino sound speed */
    set_neutrino_sound_speeds(&m, &us, &pcs);

    /* Starting and final redshifts */
    const double a_start = 1.0 / 128.0;
    const double a_final = 1.0;

    printf("Integrating cosmological tables.\n");

    /* Integrate the cosmological tables */
    integrate_cosmology_tables(&m, &us, &pcs, &tab, a_start, a_final, 1000);

    /* Get the Hubble rate at a_start */
    double H_start = get_H_of_a(&tab, a_start);
    double H0 = get_H_of_a(&tab, 1.0);
    double f_nu_nr_tot_0 = get_f_nu_nr_tot_of_a(&tab, 1.0);
    printf("H(a_start) = %g 1/U_t\n", H_start);
    printf("H(a_today) = %.10g km/s/Mpc\n", H_start / H0 * 100 * m.h);
    printf("f_nu_nr_tot(a_today) = %.10g\n", f_nu_nr_tot_0);

    printf("Integrating fluid equations.\n");

    /* Prepare integrating the fluid equations */
    const double tol = 1e-12;
    const double hstart = 1e-12;
    prepare_fluid_integrator(&m, &us, &pcs, &tab, tol, hstart);

    /* Integrate */
    int k_size = 10;
    double log_k_min = log(1e-5);
    double log_k_max = log(1e2);
    for (int i=0; i<k_size; i++) {
        double k = exp(log_k_min + i * (log_k_max - log_k_min) / k_size);

        double delta_n[1] = {exp(-k*k)};
        double gn[1] = {0.6};
        double Dn[1];

        gfac.k = k;
        gfac.delta_b = 1.0;
        gfac.delta_c = 1.0;
        gfac.delta_n = delta_n;
        gfac.gb = 1.0;
        gfac.gc = 1.0;
        gfac.gn = gn;
        gfac.Dn = Dn;

        integrate_fluid_equations(&m, &us, &pcs, &tab, &gfac, a_start, a_final);

        printf("%g %g %g %g\n", gfac.k, gfac.Dc, gfac.Db, gfac.Dn[0]);
    }

    /* Done with integration */
    free_fluid_integrator();

    /* Free the cosmological spline and tables */
    free_cosmology_tables(&tab);


    printf("All done.\n");

    return 0;
}
