/*
 * FermiPotentials.h
 *
 *  Created on: 28 feb. 2018
 *      Author: dirac89
 */

#ifndef FermiPotentials_H_
#define FermiPotentials_H_

class FermiPotentials {
private:
	int s, n, n2, li, lj, mi, mj;
	double r1, theta1, As1, wave1, wave2;
	double Ap1, Dwave1, Dwave2;
public:
	FermiPotentials(int s1, int n1, int n21, int li1, int lj1, int mi1, int mj1, double r11, double theta11, double As11, double wave11,double wave21, double Ap11, double Dwave11, double Dwave21);
	virtual ~FermiPotentials();
	double Vs ();
	double Vp ();
	double Vsp ();
	/*void setWave2(double givenwave2);
	void setDWave2(double givendwave2);
	void setWave1(double givenwave1);
	void setDWave1(double givendwave1);

	*/
	
};

#endif /* THETADERIVATE_H_ */

