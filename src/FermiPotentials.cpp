/*
 * FermiPotentials.cpp
 *
 *  Created on: 28 feb. 2018
 *      Author: dirac89
 */

#include "FermiPotentials.h"

#define pi  4.*atan(1.)

FermiPotentials::FermiPotentials(int s1, int n1, int n21, int li1, int lj1, int mi1, int mj1, double r11, double theta11, double As11, double wave11,double wave21, double Ap11, double Dwave11, double Dwave21) {
    s_=s1;
    n_=n1;
    n2_=n21;
    li_=li1;
    lj_=lj1;
    mi_=mi1;
    mj_=mj1;
    r1_=r11;
    theta1_=theta11;
    as1_=As11;
    wave1_=wave11;
    wave2_=wave21;
    ap1_=Ap11;
    dWave1_=Dwave11;
    dWave2_=Dwave21;
}

double FermiPotentials::Vs() const {

    double x1 = cos(theta1_) ;
    if (( li_>2) && (lj_>2 )) {  //""" Degenerated Part

        double Vs1 = 2.*pi*as1_*Spherical(li_, mi_, x1) *Spherical(lj_, mj_, x1)
		         *gsl_sf_hydrogenicR (n_, li_, 1., r1_) *gsl_sf_hydrogenicR (n2_, lj_, 1., r1_) ;

        return Vs1 ;
    }

    if (( li_<=2)&& (lj_>2 )) {  //"Degenerated Part

        double Vs1 = 2.*pi*as1_*Spherical(li_, mi_, x1) *Spherical(lj_, mj_, x1)
		         *gsl_sf_hydrogenicR (n2_, lj_, 1., r1_) *wave1_;
        return Vs1;
    }

    if (( li_>2)&& (lj_<=2 )) {  //"Degenerated Part

        double Vs1 = 2.*pi*as1_*Spherical(li_, mi_, x1) *Spherical(lj_, mj_, x1)
		         *gsl_sf_hydrogenicR (n_, li_, 1., r1_) *wave2_;
        return Vs1;
    }

    else { //"Non-Degenerated Part
        double Vs1 = 2.*pi* as1_* Spherical(li_, mi_, x1) *Spherical(lj_, mj_, x1) * wave1_* wave2_ ;
        return Vs1;
    }

}


double FermiPotentials::Vp() const{

    double x1 = cos(theta1_) ;

    if (( li_>2)&& (lj_>2 )) {  //"""""""""""" Degenerated Part

        double VpA1 = 6.*pi*ap1_*Spherical(li_, mi_, x1)*Spherical(lj_, mj_, x1)
              *DRnl(n_, li_, r1_)* DRnl(n2_, lj_, r1_) ;

        double VpB1 = 6.*pi*ap1_* DOlm(li_, mi_, theta1_)* DOlm(lj_, mj_, theta1_)* gsl_sf_hydrogenicR (n2_, lj_, 1., r1_)
              * gsl_sf_hydrogenicR (n_, li_, 1., r1_) *pow(r1_, -2.) ;

        double VpC1 =6.*pi*ap1_* DPhilm(li_, mi_, theta1_)* DPhilm(lj_, mj_, theta1_)
              *gsl_sf_hydrogenicR (n_, li_, 1., r1_)* gsl_sf_hydrogenicR (n2_, lj_, 1., r1_)*pow(r1_, -2.) ;

        return  (VpA1+VpB1+VpC1);
    }

    if (( li_<=2 )&&(lj_>2 )) {  //"""""""""""" Degenerated Part

        double VpA1 = 6.*pi*ap1_* Spherical(li_, mi_, x1)*Spherical(lj_, mj_, x1)*DRnl(n2_, lj_, r1_)*dWave1_ ;

        double VpB1 = 6.*pi*ap1_* DOlm(li_, mi_, theta1_)* DOlm(lj_, mj_, theta1_)* gsl_sf_hydrogenicR (n2_, lj_, 1., r1_)
              * wave1_ *pow(r1_, -2.) ;

        double VpC1 =6.*pi*ap1_* DPhilm(li_, mi_, theta1_)* DPhilm(lj_, mj_, theta1_)*
                gsl_sf_hydrogenicR (n2_, lj_, 1., r1_)*wave1_*pow(r1_, -2.) ;

        return  (VpA1+VpB1+VpC1) ;
    }

    if (( li_>2 )&&(lj_<=2 )) {  //"""""""""""" Degenerated Part

        double VpA1 = 6.*pi*ap1_* Spherical(li_, mi_, x1)*Spherical(lj_, mj_, x1)*DRnl(n_, li_, r1_)*dWave2_ ;

        double VpB1 = 6.*pi*ap1_* DOlm(li_, mi_, theta1_)* DOlm(lj_, mj_, theta1_)* gsl_sf_hydrogenicR (n_, li_, 1., r1_)
              * wave2_ *pow(r1_, -2.) ;

        double VpC1 =6.*pi*ap1_* DPhilm(li_, mi_, theta1_)* DPhilm(lj_, mj_, theta1_)*
                gsl_sf_hydrogenicR (n_, li_, 1., r1_)*wave2_*pow(r1_, -2.) ;

        return  (VpA1+VpB1+VpC1);
    }

    else { //"""""""""" Non-Degenerated Part

        double VpA1 = 6.*pi*ap1_*Spherical(li_, mi_, x1)*Spherical(lj_, mj_, x1)* dWave1_* dWave2_ ;

        double VpB1 = 6.*pi*ap1_* DOlm(li_, mi_, theta1_)* DOlm(lj_, mj_, theta1_)* wave1_* wave2_* pow(r1_, -2.) ;

        double VpC1 = 6.*pi*ap1_* DPhilm(li_, mi_,theta1_)* DPhilm(lj_, mj_, theta1_)*wave1_*wave2_*pow(r1_, -2.);

        return  (VpA1+VpB1+VpC1);
    }
}


double FermiPotentials::Vsp() const {
    return Vs() + s_ * Vp() ;
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
