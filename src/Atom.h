/*
 * Atom.h
 *
 *  Created on: 28 feb. 2018
 *      Author: dirac89
 */

#ifndef ATOM_H_
#define ATOM_H_

class Atom {
private:
	int n,l;
public:
	Atom(int n2, int l2);

	double E_Rb();
	double Vfield(int li, int lj, int mi, int mj, double radial, double strength);
	double Angular_dc_field(int l, int l1, int m, int m1);
    
    virtual ~Atom();
};

#endif /* ATOM_H_ */
