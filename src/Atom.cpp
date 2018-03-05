/*
 * Atom.cpp
 *
 *  Created on: 28 feb. 2018
 *      Author: dirac89
 */

#include "Atom.h"
#include <math.h>

Atom::Atom(int n2, int l2) {
    // TODO Auto-generated constructor stub
    n=n2;
    l=l2;
}

Atom::~Atom() {
    // TODO Auto-generated destructor stub
}

double Atom::E_Rb(){
    n=this->n;
    l=this->l;
    double muns_Rb = 3.1311804 + 0.1745312* pow( n - 3.1311804 , -2. ) ;
    double Ens_Rb = -0.5* pow( n - muns_Rb , -2. ) ;     // s non-Degenerated energy

    double munp_Rb = 2.6482793 + 0.2925324* pow( n - 2.6482793 , -2. ) ;
    double Enp_Rb = -0.5* pow( n - munp_Rb , -2. ) ;    // p non-Degenerated energy

    double mund_Rb = 1.3472787  - 0.5994376* pow( n - 1.3472787 , -2. ) ;
    double End_Rb = -0.5* pow( n - mund_Rb , -2. ) ;    // d non-Degenerated energy

    double En = -0.5* pow( n , -2. ) ;             // l=l_high Degenerated energy (Hydrogen)

    if       ((l==0)) { return Ens_Rb ; }
    else if  ((l==1)) { return Enp_Rb ; }
    else if  ((l==2)) { return End_Rb ; }
    else              { return En ;}

}

double Atom::Vfield (int li, int lj, int mi, int mj, double radial, double strength){


    return strength*Angular_dc_field(li, lj, mi,  mj)*radial;
}

double Atom::Angular_dc_field(int l, int l1, int m, int m1){
    double term1=0.;

    if (m==m1){
        if (l==l1-1){
            term1=sqrt((l1+m1)*(l1-m1)/(4.*l1*l1-1));
        }
        if (l==l1+1){
            term1=sqrt((l1+m1+1)*(l1-m1+1)/(2.*l1+1)/(2.*l1+3));
        }
        return term1;
    }

    else{ return 0.;} ;
}
