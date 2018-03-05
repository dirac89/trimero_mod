/*
 * FermiPotentials.cpp
 *
 *  Created on: 28 feb. 2018
 *      Author: dirac89
 */

#include "FermiPotentials.h"
#include "Laplacian.h"


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

FermiPotentials::FermiPotentials(int s1, int n1, int n21, int li1, int lj1, int mi1, int mj1, double r11, double theta11, double As11, double wave11,double wave21, double Ap11, double Dwave11, double Dwave21) {
	// TODO Auto-generated constructor stub
s=s1; 
n=n1;
n2=n21; 
li=li1; 
lj=lj1;
mi=mi1;
mj=mj1;
r1=r11; 
theta1=theta11; 
As1=As11;
wave1=wave11;
wave2=wave21;
Ap1=Ap11; 
Dwave1=Dwave11; 
Dwave2=Dwave21;
}

FermiPotentials::~FermiPotentials() {
	// TODO Auto-generated destructor stub
}



double FermiPotentials::Vs(){
	
	theta1=this->theta1;
	li=this->li;
	lj=this->lj;
	As1=this->As1;
	mi=this->mi;
	mj=this->mj;
	n=this->n;
	r1=this->r1;
	n2=this->n2;
	wave1=this->wave1;
	wave2=this->wave2;
	
	
	double x1 = cos(theta1) ; 
     if (( li>2) && (lj>2 )) {  //""" Degenerated Part       

          double Vs1 = 2.*pi*As1*Spherical(li, mi, x1) *Spherical(lj, mj, x1) 
		 *gsl_sf_hydrogenicR (n, li, 1., r1) *gsl_sf_hydrogenicR (n2, lj, 1., r1) ;

               return Vs1 ;
         }

        if (( li<=2)&& (lj>2 )) {  //"Degenerated Part       

	   double Vs1 = 2.*pi*As1*Spherical(li, mi, x1) *Spherical(lj, mj, x1) 
		 *gsl_sf_hydrogenicR (n2, lj, 1., r1) *wave1;
		 return Vs1; 
         } 
         
	 if (( li>2)&& (lj<=2 )) {  //"Degenerated Part       

	   double Vs1 = 2.*pi*As1*Spherical(li, mi, x1) *Spherical(lj, mj, x1) 
		 *gsl_sf_hydrogenicR (n, li, 1., r1) *wave2;
		 return Vs1; 
         }
                  
         else { //"Non-Degenerated Part
              double Vs1 = 2.*pi* As1* Spherical(li, mi, x1) *Spherical(lj, mj, x1) * wave1* wave2 ;
				return Vs1; 
         }   
         
	}


double FermiPotentials::Vp(){
	
	theta1=this->theta1;
	li=this->li;
	lj=this->lj;
	As1=this->As1;
	mi=this->mi;
	mj=this->mj;
	n=this->n;
	r1=this->r1;
	n2=this->n2;
	wave1=this->wave1;
	wave2=this->wave2;
	Dwave1=this->Dwave1;
	Dwave2=this->Dwave2;
	
	double x1 = cos(theta1) ;  
  
  if (( li>2)&& (lj>2 )) {  //"""""""""""" Degenerated Part       
    
    double VpA1 = 6.*pi*Ap1*Spherical(li, mi, x1)*Spherical(lj, mj, x1)
      *DRnl(n, li, r1)* DRnl(n2, lj, r1) ;
    
    double VpB1 = 6.*pi*Ap1* DOlm(li, mi, theta1)* DOlm(lj, mj, theta1)* gsl_sf_hydrogenicR (n2, lj, 1., r1)
      * gsl_sf_hydrogenicR (n, li, 1., r1) *pow(r1, -2.) ;
    
    double VpC1 =6.*pi*Ap1* DPhilm(li, mi, theta1)* DPhilm(lj, mj, theta1)
      *gsl_sf_hydrogenicR (n, li, 1., r1)* gsl_sf_hydrogenicR (n2, lj, 1., r1)*pow(r1, -2.) ;
     
    return  (VpA1+VpB1+VpC1);
  }
  
  if (( li<=2 )&&(lj>2 )) {  //"""""""""""" Degenerated Part       
    
    double VpA1 = 6.*pi*Ap1* Spherical(li, mi, x1)*Spherical(lj, mj, x1)*DRnl(n2, lj, r1)*Dwave1 ;
    
    double VpB1 = 6.*pi*Ap1* DOlm(li, mi, theta1)* DOlm(lj, mj, theta1)* gsl_sf_hydrogenicR (n2, lj, 1., r1)
      * wave1 *pow(r1, -2.) ;
    
    double VpC1 =6.*pi*Ap1* DPhilm(li, mi, theta1)* DPhilm(lj, mj, theta1)*
      gsl_sf_hydrogenicR (n2, lj, 1., r1)*wave1*pow(r1, -2.) ; 
    
    return  (VpA1+VpB1+VpC1) ; 
  }
  
  if (( li>2 )&&(lj<=2 )) {  //"""""""""""" Degenerated Part       
    
    double VpA1 = 6.*pi*Ap1* Spherical(li, mi, x1)*Spherical(lj, mj, x1)*DRnl(n, li, r1)*Dwave2 ;
    
    double VpB1 = 6.*pi*Ap1* DOlm(li, mi, theta1)* DOlm(lj, mj, theta1)* gsl_sf_hydrogenicR (n, li, 1., r1)
      * wave2 *pow(r1, -2.) ;
    
    double VpC1 =6.*pi*Ap1* DPhilm(li, mi, theta1)* DPhilm(lj, mj, theta1)*
      gsl_sf_hydrogenicR (n, li, 1., r1)*wave2*pow(r1, -2.) ; 
    
    return  (VpA1+VpB1+VpC1);    
  }
  
  else { //"""""""""" Non-Degenerated Part
    
    double VpA1 = 6.*pi*Ap1*Spherical(li, mi, x1)*Spherical(lj, mj, x1)* Dwave1* Dwave2 ;
    
    double VpB1 = 6.*pi*Ap1* DOlm(li, mi, theta1)* DOlm(lj, mj, theta1)* wave1* wave2* pow(r1, -2.) ; 
    
    double VpC1 = 6.*pi*Ap1* DPhilm(li, mi,theta1)* DPhilm(lj, mj, theta1)*wave1*wave2*pow(r1, -2.);
    
    return  (VpA1+VpB1+VpC1); 
  } 
	}   
	
	
double FermiPotentials::Vsp(){
	
	s=this->s;

		return FermiPotentials::Vs() + s*FermiPotentials::Vp() ;
		
		}
		
/*
void FermiPotentials::setWave2(double givenwave2){
	
	wave2=givenwave2;
	
	}


void FermiPotentials::setDWave2(double givendwave2){
	
	Dwave2=givendwave2;
	
	}


void FermiPotentials::setWave1(double givenwave1){
	
	wave1=givenwave1;
	
	}


void FermiPotentials::setDWave1(double givendwave1){
	
	Dwave1=givendwave1;
	
	}
	
*/
