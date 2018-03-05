/*
 * Laplacian.cpp
 *
 *  Created on: 28 feb. 2018
 *      Author: dirac89
 */

#ifndef LAPLACIAN_H_
#define LAPLACIAN_H_

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include <gsl/gsl_blas.h>            // for blas usage (gsl)
#include <gsl/gsl_math.h>            // for mathematical functions (gsl)
#include <gsl/gsl_vector.h>          // for vectors (gsl)
#include <gsl/gsl_matrix.h>          // for matrices (gsl)
#include <gsl/gsl_eigen.h>           // for eigensystems (gsl)
#include <gsl/gsl_sf_coulomb.h>      // for Coulomb wavefunctions (gsl)
#include <gsl/gsl_sf_legendre.h>     // for Legendre polynomials (gsl)
#include <gsl/gsl_sf_laguerre.h>     // for Laguerre polynomials (gsl)
#include <gsl/gsl_sf_gamma.h>        // for gamma function and factorial (gsl)
#include <gsl/gsl_errno.h>           // for errors (gsl)
#include <gsl/gsl_interp.h>          // for general interpolation (gsl)
#include <gsl/gsl_spline.h>          // for spline interpolation (gsl)
#include <gsl/gsl_sort.h>            // for simple sorting (gsl)
#include <gsl/gsl_sort_vector.h>     // for vector elements sorting (gsl)
#include <gsl/gsl_sf_log.h>          // for vector elements sorting (gsl)
#include <gsl/gsl_sf_exp.h>

#define pi  4.*atan(1.)

double DOlm(int l, int m, double theta);
double DPhilm(int l, int m, double theta);
double DRnl(int n, int l, double r);
double Spherical(int l, int m, double x1);

#endif
