/*
 * Laplacian.cpp
 *
 *  Created on: 28 feb. 2018
 *      Author: dirac89
 */



#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>


	
#include <gsl/gsl_blas.h>                                                        // for blas usage (gsl)
#include <gsl/gsl_math.h>                                                      // for mathematical functions (gsl)
#include <gsl/gsl_vector.h>                                                     // for vectors (gsl)
#include <gsl/gsl_matrix.h>                                                     // for matrices (gsl)
#include <gsl/gsl_eigen.h>                                                      // for eigensystems (gsl)
#include <gsl/gsl_sf_coulomb.h>                                            // for Coulomb wavefunctions (gsl)
#include <gsl/gsl_sf_legendre.h>                                            // for Legendre polynomials (gsl)
#include <gsl/gsl_sf_laguerre.h>                                            // for Laguerre polynomials (gsl)
#include <gsl/gsl_sf_gamma.h>                                             // for gamma function and factorial (gsl)
#include <gsl/gsl_errno.h>                                                      // for errors (gsl)
#include <gsl/gsl_interp.h>                                                     // for general interpolation (gsl)
#include <gsl/gsl_spline.h>                                                     // for spline interpolation (gsl)
#include <gsl/gsl_sort.h>                                                        // for simple sorting (gsl)
#include <gsl/gsl_sort_vector.h>                                             // for vector elements sorting (gsl)
#include <gsl/gsl_sf_log.h>                                             // for vector elements sorting (gsl)
#include <gsl/gsl_sf_exp.h>

#define pi  4.*atan(1.)

double DOlm(int l, int m, double theta);
double DPhilm(int l, int m, double theta);
double DRnl(int n, int l, double r);
double Spherical(int l, int m, double x1);




double DOlm(int l, int m, double theta){
	
	
	double x = cos(theta) ; // NOTE : neq +/-1.
	  double term2=0. ;
	  double term1=0. ;
	  if(l>0){
	    double C1 = sqrt((2.*l+1.)*gsl_sf_fact(l-m)) ;
	    double C2 = sqrt((4.*pi)*gsl_sf_fact(l+m)) ;
	    double Clm = C1/C2*0.5 ;

	    if (m+1>=0) {
	      if (m+1<=l) {
		term1 = gsl_sf_legendre_Plm(l,m+1,x) ; }
	      if (m+1>l) {
		term1 = 0. ; }
	    }

	    if (m+1<0) {
	      if (-m-1<=l) {
		term1 = pow(-1.,1+m)*gsl_sf_fact(l+1+m)/gsl_sf_fact(l-1-m)*gsl_sf_legendre_Plm(l,-m-1,x) ; }
	      if (-m-1>l) {
		term1 = 0. ; }
	    }

	    if ((m-1<0)) {
	      if (-m+1<=l){
	      term2 = (l+m)*(l-m+1)*
		gsl_sf_legendre_Plm(l,1-m,x)*pow(-1.,1-m)*gsl_sf_fact(l-1+m)/gsl_sf_fact(l+1-m); }}
	    if (-m+1>l){ term2=0.;}

	    if ((m-1>=0)) {
	      if (m-1<=l){
		term2 = (l+m)*(l-m+1)*gsl_sf_legendre_Plm(l,m-1,x) ; }
	      if (m-1>l){term2=0.;}

	    }

	    return Clm*(term1-term2);}

	  else {return 0.;}
}


double DPhilm(int l, int m, double theta){

  double x = cos(theta) ; // NOTE : neq +/-1.
  double term2=0. ;
  double term1=0.;

  if(l>0){
    double C1 = sqrt((2.*l+1.)*gsl_sf_fact(l-m)) ;
    double C2 = sqrt((4.*pi)*gsl_sf_fact(l+m)) ;
    double Clm = -C1/C2*0.5 ;
    if ((m+1>=0)) {
      if (m+1<= l+1){
       term1 = gsl_sf_legendre_Plm(l+1,m+1,x);
      }
      if (m+1> l+1){
	term1 = 0.;
      }
    }

    if ((m+1<0)) { // if m+1<0 we have to take care of the Associated Legendre Function
       if (-m-1<= l+1){
       term1 = gsl_sf_legendre_Plm(l+1,-m-1,x)*pow(-1.,1+m)*gsl_sf_fact(l+1+1+m)/gsl_sf_fact(l+1-1-m) ;
       }
      if (-m-1> l+1){
	 term1 = 0.;
      }
    }

    if ((m-1<0)) { // if m-1<0 we have to take care of the Associated Legendre Function
      if (-m+1<=l+1){
	term2 = (l-m+1)*(l+2-m)*
	  gsl_sf_legendre_Plm(l+1,1-m,x)*pow(-1.,1-m)*gsl_sf_fact(l+1-1+m)/gsl_sf_fact(l+1+1-m); }
      if (-m+1>l+1){
	term2 =0.;  }
    }
    if ((m-1>=0)) {
      if (m-1<=l+1){
	term2 = (l-m+1)*(l-m+2)*gsl_sf_legendre_Plm(l+1,m-1,x) ; }
      if (m-1>l+1){
	term2 = 0.;}
    }
    return Clm*(term1+term2);}

  else {return 0.;}

}


double DRnl(int n, int l, double r){

	if(l<n-1){double xn = 2.*r/n ;
	                 double N1 = 2.* sqrt(gsl_sf_fact(n-l-1)) ;
	                 double N2 = pow(n,2.)* sqrt(gsl_sf_fact(n+l)) ;
	                 double Nnl = N1/N2 ;
	                 double term1 = ((l/r)-(1./n))*gsl_sf_hydrogenicR(n, l, 1., r) ;
	                 double term2 = pow(xn,l)*exp(-r/n)*gsl_sf_laguerre_n(n-l-2, 2*l+2, xn) ;
	                 return  term1 - Nnl* term2*2./n ;
	      }
	      else { return 0.; }

	/*
	COMMENTS :  Derivative of the radial hydrogenic wave functions. Note that r is not 0 !
	*/
}

double Spherical(int l, int m, double x1){

     if(m<0){
			return  gsl_sf_legendre_sphPlm(l,-m,x1)*pow(-1,-m)*gsl_sf_fact(l+m)/gsl_sf_fact(l-m);
			}  
	else { 
		
		return gsl_sf_legendre_sphPlm(l, m,x1);
				}

}
 

