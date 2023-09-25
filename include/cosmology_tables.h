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

#ifndef COSMOLOGY_TABLES_H
#define COSMOLOGY_TABLES_H

#include "units.h"

struct model {
    double h;
    double Omega_b;
    double Omega_c;
    double Omega_k;
    double N_ur;
    int N_nu;
    double *M_nu;
    double *deg_nu;
    double *c_s_nu;
    double T_nu_0;
    double T_CMB_0;
    double w0;
    double wa;

    double Omega_dcdmdr_0;
    double Gamma_dcdm;

    /* Do the simulation particles not have masses that vary with w_nu(a)? */
    int sim_neutrino_nonrel_masses;
    /* Does the simulation have a Hubble rate with Omega_nu(a) = constant? */
    int sim_neutrino_nonrel_Hubble;
};

struct cosmology_tables {
    double *avec;
    double *Avec;
    double *Bvec;
    double *Hvec;
    double *f_nu_nr;
    double *f_nu_nr_tot;
    double *Omega_dcdm;
    double *Omega_dr;
    int size;
};

void integrate_cosmology_tables(struct model *m, struct units *us,
                                struct physical_consts *pcs,
                                struct cosmology_tables *tab, double a_start,
                                double a_final, int size);
void free_cosmology_tables(struct cosmology_tables *tab);
void set_neutrino_sound_speeds(struct model *m, struct units *us,
                               struct physical_consts *pcs);

double get_H_of_a(struct cosmology_tables *tab, double a);
double get_f_nu_nr_tot_of_a(struct cosmology_tables *tab, double a);
double get_Omega_dcdm_of_a(struct cosmology_tables *tab, double a);
double get_Omega_dr_of_a(struct cosmology_tables *tab, double a);

#endif
