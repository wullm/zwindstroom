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
#include <gsl/gsl_roots.h>

#include "../include/mg_fR.h"

double fR_Q(double R, void *p) {
    struct fR_shoot_params *params = (struct fR_shoot_params *) p;
    double b = params->b;
    double Lambda = params->Lambda;
    double R_on_b = R / b;
    return - exp(-R_on_b / Lambda) * (2.0 * R_on_b + 4.0 * Lambda) + 4.0 * Lambda - R;
}

double fR_shoot_func(double R, void *p) {
    struct fR_shoot_params *params = (struct fR_shoot_params *) p;
    double b = params->b;
    double Lambda = params->Lambda;
    double Q = params->Q;
    double R_on_b = R / b;
    
    double f = R - 2.0 * Lambda * (1.0 - exp(-R_on_b / Lambda));
    double fR = 1.0 - 2.0 / b * exp(-R_on_b / Lambda);
    
    return fR * R - 2 * f - Q;
}

double fR_shoot_func_deriv(double R, void *p) {
    struct fR_shoot_params *params = (struct fR_shoot_params *) p;
    double b = params->b;
    double Lambda = params->Lambda;
    double Q = params->Q;
    double R_on_b = R / b;
    
    double fR = 1.0 - 2.0 / b * exp(-R_on_b / Lambda);
    double fRR = 2.0 / (Lambda * b * b) *  exp(-R_on_b / Lambda);
    
    return fRR * R - fR;
}

void fR_shoot_func_fdf(double R, void *p, double *f, double *df) {
    *f = fR_shoot_func(R, p);
    *df = fR_shoot_func_deriv(R, p);
}

double fR_solve_for_R(double b, double Lambda, double Q) {
    /* Set parameters for the root finding algorithm */
    struct fR_shoot_params fR_pars;
    fR_pars.b = b;
    fR_pars.Lambda = Lambda;
    fR_pars.Q = Q;

    /* The initial guess for R */
    const double R_ini = 1e5;
    double R;

    /* Set up the root finder */
    const gsl_root_fdfsolver_type *T = gsl_root_fdfsolver_newton;
    gsl_root_fdfsolver *s = gsl_root_fdfsolver_alloc(T);
    gsl_function_fdf F;
    F.f = &fR_shoot_func;
    F.df = &fR_shoot_func_deriv;
    F.fdf = &fR_shoot_func_fdf;
    F.params = &fR_pars;
    gsl_root_fdfsolver_set(s, &F, R_ini);

    /* Error handling */
    int status = GSL_CONTINUE;
    int iter = 0;
    const int max_iter = 1000;
    double residual;
            
    while (status == GSL_CONTINUE && iter < max_iter) {
        iter++;
        status = gsl_root_fdfsolver_iterate(s);
        R = gsl_root_fdfsolver_root(s);
        residual = fR_shoot_func(R, &fR_pars);
        status = gsl_root_test_residual(residual, 1e-6);
    }


    /* Free the root finder */
    gsl_root_fdfsolver_free (s);
    
    return R;
}