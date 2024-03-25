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

#ifndef MG_FT_H
#define MG_FT_H

#include "../include/cosmology_tables.h"
#include "../include/units.h"
#include "../include/strooklat.h"

struct fT_df_params {
  double Omega_CMB;
  double Omega_ur;
  double Omega_c;
  double Omega_b;
  double Omega_lambda;
  double *Omega_nu;
  double *w_nu;
  struct strooklat *spline;
  struct model *m;
  struct units *us;
  struct physical_consts *pcs;
  int size;
};

int fT_dfunc(double a, const double y[], double f[], void *params_ptr);

#endif
