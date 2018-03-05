/*
 * parameter.h
 *
 *  Created on: 28 feb. 2018
 *      Author: dirac89
 */

//#ifndef PARAMETER_H_
//#define PARAMETER_H_

#define pi  4.*atan(1.)
#define ex  exp(1.)
#define EhtoHz  (6.579683920)*pow(10,15)                           // [Hartree(Eh)] = (6.579683) 10^15 [Hz]

#define EhtoGHz (6.579683920)*pow(10,6)                           // [Hartree(Eh)] = (6.579683) 10^15 [Hz]
#define eVtoHz  (1.509190450)* pow(10,33)                           // [eV] = (1.509190450) 10^33 [Hz]
#define eVtoAu  (0.003674932379)                                               // [eV] = (0.003674932379) [au]

#define autovperm (5.14220652e11) //*pow(10,11) // From au to V/m
#define vpermtoau (1.9446929e-12)


int pqn=35;
int oqn=3;


double dc_field=300.; // Aqui se define el campo electrico en V/m
double dc_field_au=dc_field*vpermtoau;// Aqui esta pasando a unidades atomicas.

//#endif /* PARAMETER_H_ */
