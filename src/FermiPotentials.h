/*
 * FermiPotentials.h
 *
 *  Created on: 28 feb. 2018
 *      Author: dirac89
 */

#ifndef FERMIPOTENTIALS_H_
#define FERMIPOTENTIALS_H_

#include "Laplacian.h"

class FermiPotentials {
public:
	FermiPotentials(int s, int n, int n2,
	                int li, int lj,
	                int mi, int mj,
	                double r1, double theta1, double As1,
	                double wave1, double wave2,
	                double Ap1, double Dwave1, double Dwave2);
	double Vs() const;
	double Vp() const;
	double Vsp() const;
	/*void setWave2(double givenwave2);
	void setDWave2(double givendwave2);
	void setWave1(double givenwave1);
	void setDWave1(double givendwave1);

	*/
private:
    int s_, n_, n2_, li_, lj_, mi_, mj_;
    double r1_, theta1_, as1_, wave1_, wave2_;
    double ap1_, dWave1_, dWave2_;
};

#endif /* THETADERIVATE_H_ */

