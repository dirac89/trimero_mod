//============================================================================
// Name        : trimero.cpp
// Author      : javieraguilerafernandez
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

//...................................... G++ LIBS
	#include <iostream>
	#include <iomanip>
	#include <fstream>
	#include <string>
	#include <ctime>
	#include <math.h>
	#include <stdlib.h>
	#include <stdio.h>

//...................................................... GSL LIBS
//
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
#include <gsl/gsl_sf_exp.h>                                             // for vector elements sorting (gsl)


//************************** Classes **************************
#include "Atom.h"
#include "FermiPotentials.h"


#include "parameter.h"

void Trimer_energies_field (int n1,double dc_field_au);

using namespace std;

//**************************************************************
int main() {
	
	      // Dimer_energies_one( 35, dc_field_au);
	      //  Dimer_energies_field(35,dc_field_au);
	      // Trimer_energies_one( 35, dc_field_au);
	      // Trimer_energies_one_field_free(35,dc_field_au);
	       Trimer_energies_field (pqn, dc_field_au);
	    
	     // Atom rubidio(n, l);
	     
	      //cout << "Energia Rubidio: " << rubidio.E_Rb()<< endl;
	     
	     
	      
	      
 
	
	return 0 ;
}
//**************************************************************


void Trimer_energies_field (int n1,double dc_field_au){
   
  Atom rubidium(n1, oqn); 

  
  
  const int lc = 2 ;//lc=2
  const int cols = 2 ;                                              
  const int rows =776 ; 
  const double theta=0.;

  const double theta1 = pi;
  int k1,k2;
  int row,row1 ;
  int s;
  double delta_ij,field_contri;
  double R, AS, AP, R1, AS1, AP1, En, VShift ;
  double wave_R, wave1_R;
  double Dwave_R, Dwave1_R;
  double wave_R1, wave1_R1;
  double Dwave_R1, Dwave1_R1;
  double R38s_R, R36d_R, R37p_R, DR38s_R, DR36d_R, DR37p_R ;
  double radial;
    
   double wave_1, wave1_2, Dwave_1, Dwave1_2;
	      
   double wave_2, wave2_1, Dwave_2, Dwave2_1;
   
   int n11, n12, n21,n22;
  
  
  gsl_matrix * As = gsl_matrix_alloc (rows,cols) ; FILE* fAs = fopen("rvsAS.dat","r")  ;  
  gsl_matrix_fscanf(fAs, As) ;  fclose(fAs) ;   
  
  gsl_matrix * Ap = gsl_matrix_alloc (rows,cols) ; FILE* fAp = fopen("rvsAP.dat","r")  ;  
  gsl_matrix_fscanf(fAp, Ap) ;  fclose(fAp) ; 
  
  gsl_matrix * As1 = gsl_matrix_alloc (rows,cols) ; FILE* fAs1 = fopen("rvsAS.dat","r")  ;  
  gsl_matrix_fscanf(fAs1, As1) ;  fclose(fAs1) ;   

  gsl_matrix * Ap1 = gsl_matrix_alloc (rows,cols) ; FILE* fAp1 = fopen("rvsAP.dat","r")  ;  
  gsl_matrix_fscanf(fAp1, Ap1) ;  fclose(fAp1) ;      
  
  cout << "Manifold n1 = "<<n1<<" plots of E(R) type are included !\n" ;


  
  gsl_matrix * R38s = gsl_matrix_alloc (rows,cols) ; FILE* fR38s  = fopen("Wavefunction/rvsR38s.dat","r")  ;  
  gsl_matrix_fscanf(fR38s, R38s)    ; fclose(fR38s) ;   
    
  gsl_matrix * R36d = gsl_matrix_alloc (rows,cols) ; FILE* fR36d  = fopen("Wavefunction/rvsR36d.dat","r") ;  
  gsl_matrix_fscanf(fR36d, R36d)    ; fclose(fR36d) ;  
    
  gsl_matrix * R37p = gsl_matrix_alloc (rows,cols) ; FILE* fR37p  = fopen("Wavefunction/rvsR37p.dat","r") ;  
  gsl_matrix_fscanf(fR37p, R37p)    ; fclose(fR37p) ; 
    
  gsl_matrix * DR38s = gsl_matrix_alloc (rows,cols) ; FILE* fDR38s  = fopen("Wavefunction/rvsDR38s.dat","r")  ;  
  gsl_matrix_fscanf(fDR38s, DR38s)    ; fclose(fDR38s) ;   
    
  gsl_matrix * DR36d = gsl_matrix_alloc (rows,cols) ; FILE* fDR36d  = fopen("Wavefunction/rvsDR36d.dat","r") ;  
  gsl_matrix_fscanf(fDR36d, DR36d)    ; fclose(fDR36d) ;  
    
  gsl_matrix * DR37p = gsl_matrix_alloc (rows,cols) ; FILE* fDR37p  = fopen("Wavefunction/rvsDR37p.dat","r") ;  
  gsl_matrix_fscanf(fDR37p, DR37p)    ; fclose(fDR37p) ; 
    
      

  ofstream store2DShift("Trimer_R_sp_wave_N35_R_300_GHz.dat") ; 
  ofstream store2DShift_au("Trimer_R_sp_wave_N35_R_300_au.dat") ;     
      
  const int max = n1*n1;
  
  //cout<<max<<endl;

  gsl_matrix * spV = gsl_matrix_calloc(max, max) ;                                                               
  gsl_vector * evalues = gsl_vector_calloc (max) ;                                                            
  gsl_matrix * evectors = gsl_matrix_calloc (max, max) ;                                                   
  gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (max) ;  
  
  int m1, m2;
  m1=0;
  m2=0;
  s=1; //s=0 only s-wave, s=1 s-wave + p-wave
 
  int li ; 
  int i = li- lc -1 ;  
  int lj ; 
  int j = lj- lc -1 ;  

  gsl_matrix * field = gsl_matrix_calloc (max, max) ;  
                                                  
  double value,r_exp;
  
  gsl_matrix * Radial = gsl_matrix_calloc(max, max) ;                                                               
  
  ifstream infile;
  infile.open ("exp_val_r.txt");
  
 
  for(i = 0; i <= n1-1-1 ; i++){ 
    infile >> value;
    gsl_matrix_set (Radial, i, i+1,value);
  }
  cout << "done reading radial field "<< endl;
  
  
  
  for(i = 0; i <= n1-1 ; i=i+2){ 
    k1=0;

    for(k1 = 1 ;k1 <= 2*i+1 ; k1++){ m1=k1-i-1;		
 
      for(j = 0 ; j <= n1-1; j=j+2){ 
	k2=0;
	
	for(k2 = 1 ;k2 <= 2*j+1 ; k2++){ m2=k2-j-1;		
	  
	  if ((i==j+1)&&(m1==m2)){r_exp = gsl_matrix_get(Radial,j,i) ;
		  
	    gsl_matrix_set (field, i*i+k1-1,j*j+k2-1, rubidium.Vfield (i, j, m1, m2, r_exp,dc_field_au) );
	    
	    
	  }
	  if ((i==j-1)&&(m1==m2)){r_exp = gsl_matrix_get(Radial,i,j) ;
		  
	    gsl_matrix_set (field, i*i+k1-1,j*j+k2-1, rubidium.Vfield (i, j, m1, m2,r_exp, dc_field_au) );
	  }
	}}
	}}
  cout << "done field matrix trimer full nlm"<< endl;
 // cout<<"Metod \t"<<rubidium.Vfield (i, j, m1, m2, r_exp,dc_field_au)<<endl;


  for(row=297;row<298;row++){
      
    row1=row;
     
    // getting the radial possition and the s an p-wave scattering lengths
     R = gsl_matrix_get(As,row,0) ; 
    AS = gsl_matrix_get(As,row,1) ;  
    AP = gsl_matrix_get(Ap,row,1) ;      
   
    R1 = gsl_matrix_get(As1,row1,0) ; 
   AS1 = gsl_matrix_get(As1,row1,1) ;  
   AP1 = gsl_matrix_get(Ap1,row1,1) ;    
        
	cout<<"Building the Hamiltonian matrix"<<endl;	
	     
	
    for(i = 0; i <= n1-1 ; i++){ 
      k1=0;	
      for(k1 = 1 ;k1 <= 2*i+1 ; k1++){ m1=k1-i-1;
	 	
	for(j = 0 ; j <= n1-1 ; j++){ 

	  // constructing the Hamiltonian  matrix 
	  //  cout<<m1<<"\t"<<m2<<"\t"<<gsl_matrix_get(spV,i*i+k1-1,j*j+k2-1)<<endl;
	  
	//  FermiPotentials fermi_1(s, n11, n12,i, j, m1, m2, R, theta, AS, AP,wave_1,wave1_2, Dwave_1,Dwave1_2);
	//  FermiPotentials fermi_2(s, n21, n22,i, j, m1, m2, R1, theta1, AS1, AP1,wave_2,wave2_1, Dwave_2,Dwave2_1); 
   

	  if  ((i<3)&&(j<3)) {
		  
	    if ((i==1)&& (j==1)){
	      wave_R = gsl_matrix_get(R37p,row,1) ; 	
	      
	      Dwave_R = gsl_matrix_get(DR37p,row,1) ;
	      
	      wave_R1 = gsl_matrix_get(R37p,row1,1) ; 	
	      
	      Dwave_R1 = gsl_matrix_get(DR37p,row1,1) ;
	      }
	    
	    if ((i==2)&&(j==2)){
	      
	      wave_R = gsl_matrix_get(R36d,row,1) ; 	
	      
	      Dwave_R = gsl_matrix_get(DR36d,row,1) ;
	      
	      wave_R1 = gsl_matrix_get(R36d,row1,1) ; 	
	      
	      Dwave_R1 = gsl_matrix_get(DR36d,row1,1) ;
	      }
	    
	    if ((i==0)&&(j==0)){
			
	      wave_R = gsl_matrix_get(R38s,row,1) ; 	
	      
	      Dwave_R = gsl_matrix_get(DR38s,row,1) ;
	      
	      wave_R1 = gsl_matrix_get(R38s,row1,1) ; 	
	      
	      Dwave_R1 = gsl_matrix_get(DR38s,row1,1) ;
	      }
	      
	    k2=0;
              
	      for(k2 = 1 ;k2 <= 2*j+1 ; k2++){ m2=k2-j-1;

	      delta_ij=0;
	      
	      if ((i==j)&&(m1==m2)){delta_ij=1.;}
	      field_contri=0.;
	      if ((i==j+1)&&(m1==m2)){field_contri=gsl_matrix_get(field,i*i+k1-1,j*j+k2-1); }
	      if ((i==j-1)&&(m1==m2)){field_contri=gsl_matrix_get(field,i*i+k1-1,j*j+k2-1); }
	      
	      
	      n11=n1+3-i;
	      n21=n1+3-i;
	      
	      n12=n1+3-j;
	      n22=n1+3-j;
	      
	      wave_1=wave_R;
	      wave1_2=wave1_R;
	      
	      Dwave_1=Dwave_R;
	      Dwave1_2=Dwave1_R;
	      
	      wave_2=wave_R1;
	      wave2_1=wave1_R1;
	      
	      Dwave_2=Dwave_R1;
	      Dwave2_1=Dwave1_R1;
	      
	     //FermiPotentials fermi_1(s, n1+3-i, n1+3-j,i, j, m1, m2, R, theta, AS, AP,wave_R,wave1_R, Dwave_R,Dwave1_R);
	     // FermiPotentials fermi_2(s, n1+3-i, n1+3-j,i, j, m1, m2, R1, theta1, AS1, AP1,wave_R1,wave1_R1, Dwave_R1,Dwave1_R1);
	      
	     Atom rubidium_mod(n11, i); 
	     FermiPotentials fermi_1(s, n11, n12,i, j, m1, m2, R, theta, AS, AP,wave_1,wave1_2, Dwave_1,Dwave1_2);
	     FermiPotentials fermi_2(s, n21, n22,i, j, m1, m2, R1, theta1, AS1, AP1,wave_2,wave2_1, Dwave_2,Dwave2_1); 
		
	      gsl_matrix_set(spV, i*i+k1-1, j*j+k2-1, rubidium_mod.E_Rb()*delta_ij + 
			     fermi_1.Vsp() +
			     fermi_2.Vsp()+field_contri);
			     
			  //   cout<<"Energia\t"<<rubidium_mod.E_Rb()<<"i"<<i<<"n11"<<n11<<endl;
	    }} 

/*

	  if  ((i<3)&&(j>2)) {
	    if (i==1){
	      wave_R = gsl_matrix_get(R37p,row,1) ; 	Dwave_R = gsl_matrix_get(DR37p,row,1) ;
	      wave_R1 = gsl_matrix_get(R37p,row1,1) ; 	Dwave_R1 = gsl_matrix_get(DR37p,row1,1) ;}
	    if (i==2){
	      wave_R = gsl_matrix_get(R36d,row,1) ; 	Dwave_R = gsl_matrix_get(DR36d,row,1) ;
	      wave_R1 = gsl_matrix_get(R36d,row1,1) ; 	Dwave_R1 = gsl_matrix_get(DR36d,row1,1) ;}
	    if (i==0){
	      wave_R = gsl_matrix_get(R38s,row,1) ; 	Dwave_R = gsl_matrix_get(DR38s,row,1) ;
	      wave_R1 = gsl_matrix_get(R38s,row1,1) ; 	Dwave_R1 = gsl_matrix_get(DR38s,row1,1) ; }

	    k2=0;
	   
	    for(k2 = 1 ;k2 <= 2*j+1 ; k2++){ m2=k2-j-1;
			
			 delta_ij=0;
	      if ((i==j)&&(m1==m2)){delta_ij=1.;}
	      field_contri=0.;
	      if ((i==j+1)&&(m1==m2)){field_contri=gsl_matrix_get(field,i*i+k1-1,j*j+k2-1); }
	      if ((i==j-1)&&(m1==m2)){field_contri=gsl_matrix_get(field,i*i+k1-1,j*j+k2-1); }
	      
	      FermiPotentials fermi_11(s, n1+3-i, n1,i, j, m1, m2, R, theta, AS, AP,wave_R,0., Dwave_R,0.);
	      FermiPotentials fermi_21(s, n1+3-i, n1,i, j, m1, m2, R1, theta1, AS1, AP1,wave_R1,0., Dwave_R1,0.);

	      gsl_matrix_set(spV, i*i+k1-1, j*j+k2-1,
			     fermi_11.Vsp()+
			     fermi_21.Vsp()+field_contri );
	    }}
	
	
	  if  ((i>2)&&(j<3)) {
	    if (j==0){
	      wave_R = gsl_matrix_get(R38s,row,1) ; 	Dwave_R = gsl_matrix_get(DR38s,row,1) ;
	      wave_R1 = gsl_matrix_get(R38s,row1,1) ; 	Dwave_R1 = gsl_matrix_get(DR38s,row1,1) ;}
	    if (j==1){
	      wave_R = gsl_matrix_get(R37p,row,1) ; 	Dwave_R = gsl_matrix_get(DR37p,row,1) ;
	      wave_R1 = gsl_matrix_get(R37p,row1,1) ; 	Dwave_R1 = gsl_matrix_get(DR37p,row1,1) ;}
	    if (j==2){
	      wave_R = gsl_matrix_get(R36d,row,1) ; 	Dwave_R = gsl_matrix_get(DR36d,row,1) ;
	      wave_R1 = gsl_matrix_get(R36d,row1,1) ; 	Dwave_R1 = gsl_matrix_get(DR36d,row1,1) ;
	    }

	    
	    k2=0;
	  
	      
	    for(k2 = 1 ;k2 <= 2*j+1 ; k2++){ m2=k2-j-1;
		  
		  delta_ij=0;
	      if ((i==j)&&(m1==m2)){delta_ij=1.;}
	      field_contri=0.;
	      if ((i==j+1)&&(m1==m2)){field_contri=gsl_matrix_get(field,i*i+k1-1,j*j+k2-1); }
	      if ((i==j-1)&&(m1==m2)){field_contri=gsl_matrix_get(field,i*i+k1-1,j*j+k2-1); }
	      
	      FermiPotentials fermi_12(s, n1, n1+3-j,i, j, m1, m2, R, theta, AS, AP, 0.,wave_R,0., Dwave_R);
	      FermiPotentials fermi_22(s, n1, n1+3-j,i, j, m1, m2, R1, theta1, AS1, AP1,0.,wave_R1,0., Dwave_R1);


	      gsl_matrix_set(spV, i*i+k1-1, j*j+k2-1,
			     fermi_12.Vsp() +
			     fermi_22.Vsp()+field_contri);
	    }}


	  if ( (i>lc) && (j>lc)){ 
	    k2=0;
	    
	  	    for(k2 = 1 ;k2 <= 2*j+1 ; k2++){ m2=k2-j-1;
	      delta_ij=0;
	      if ((i==j)&&(m1==m2)){delta_ij=1.; }

	      field_contri=0.;
	      if ((i==j+1)&&(m1==m2)){field_contri=gsl_matrix_get(field,i*i+k1-1,j*j+k2-1); }
	      if ((i==j-1)&&(m1==m2)){field_contri=gsl_matrix_get(field,i*i+k1-1,j*j+k2-1); }
	      
	      Atom rubidium35(n1,lc+1);
	      
	      FermiPotentials fermi_13(s, n1,n1, i, j, m1, m2, R,  theta, AS, AP,  0., 0., 0., 0.);
	      FermiPotentials fermi_23(s, n1,n1, i, j, m1, m2, R1, theta1, AS1, AP1, 0., 0., 0., 0.);

	      gsl_matrix_set(spV, i*i+k1-1, j*j+k2-1, rubidium35.E_Rb()*delta_ij + 
			     fermi_13.Vsp()+
			     fermi_23.Vsp()+field_contri) ;
	    }}
	*/}
	}}


	cout<<"Before Diagonalization"<<endl;

    gsl_eigen_symmv (spV, evalues, evectors, w) ;                                                          
	gsl_eigen_symmv_sort (evalues, evectors, GSL_EIGEN_SORT_ABS_ASC) ;   
	
	cout<<"After Diagonalization"<<endl;	
	
	//cout<<rubidium.E_Rb()<<endl;	
	
	   

    VShift = gsl_vector_get (evalues, max-1) ;

    //---------------- 2D PLOTS -----------

    store2DShift.precision(15) ; store2DShift << R << "\t" ; 
    store2DShift_au.precision(15) ; store2DShift_au << R << "\t" ;
     

    for (i=0; i<=max-1;i++){store2DShift << (gsl_vector_get (evalues, max-1-i )-rubidium.E_Rb())*EhtoGHz << "\t";} 
    store2DShift<< endl ;
    
    for (i=0; i<=max-1;i++){store2DShift_au << (gsl_vector_get (evalues, max-1-i )-rubidium.E_Rb()) << "\t";} 
    store2DShift_au<< endl ;

    } 
 
      
  gsl_eigen_symmv_free (w) ;    
  gsl_vector_free (evalues);      
  gsl_matrix_free (evectors);     
  gsl_matrix_free (spV) ;                            
  gsl_matrix_free (field) ;
  gsl_matrix_free (Radial) ;

  store2DShift.close();  
   
 
  return ;
}

