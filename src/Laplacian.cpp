
#include "Laplacian.h"

double DOlm(int l, int m, double theta)
{
	double x = cos(theta) ; // NOTE : neq +/-1.
    double term2=0. ;
    double term1=0. ;

    if(l>0)
    {
        double C1 = sqrt((2.*l+1.)*gsl_sf_fact(l-m)) ;
        double C2 = sqrt((4.*pi)*gsl_sf_fact(l+m)) ;
        double Clm = C1/C2*0.5 ;

        if (m+1>=0)
        {
            if (m+1<=l)
            {
                term1 = gsl_sf_legendre_Plm(l,m+1,x) ;
            }
            if (m+1>l)
            {
                term1 = 0. ;
            }
        }

        if (m+1<0)
        {
            if (-m-1<=l)
            {
                term1 = pow(-1.,1+m)*gsl_sf_fact(l+1+m)/gsl_sf_fact(l-1-m)*gsl_sf_legendre_Plm(l,-m-1,x) ;
            }
            if (-m-1>l)
            {
                term1 = 0. ;
            }
        }

        if ((m-1<0))
        {
            if (-m+1<=l)
            {
                term2 = (l+m)*(l-m+1)*
                        gsl_sf_legendre_Plm(l,1-m,x)*pow(-1.,1-m)*gsl_sf_fact(l-1+m)/gsl_sf_fact(l+1-m);
            }
        }
        if (-m+1>l)
        {
        	term2=0.;
        }

        if ((m-1>=0))
        {
            if (m-1<=l)
            {
                term2 = (l+m)*(l-m+1)*gsl_sf_legendre_Plm(l,m-1,x) ;
            }
            if (m-1>l)
            {
            	term2=0.;
            }

        }

        return Clm*(term1-term2);
    }

    else
    {
    	return 0.;
    }
}

double DPhilm(int l, int m, double theta)
{
    double x = cos(theta) ; // NOTE : neq +/-1.
    double term2=0. ;
    double term1=0.;

    if(l>0)
    {
        double C1 = sqrt((2.*l+1.)*gsl_sf_fact(l-m)) ;
        double C2 = sqrt((4.*pi)*gsl_sf_fact(l+m)) ;
        double Clm = -C1/C2*0.5 ;

        if ((m+1>=0))
        {
            if (m+1<= l+1)
            {
                term1 = gsl_sf_legendre_Plm(l+1,m+1,x);
            }
            if (m+1> l+1)
            {
                term1 = 0.;
            }
        }

        if ((m+1<0)) // if m+1<0 we have to take care of the Associated Legendre Function
        {
            if (-m-1<= l+1)
            {
                term1 = gsl_sf_legendre_Plm(l+1,-m-1,x)*pow(-1.,1+m)*gsl_sf_fact(l+1+1+m)/gsl_sf_fact(l+1-1-m) ;
            }
            if (-m-1> l+1)
            {
                term1 = 0.;
            }
        }

        if ((m-1<0)) // if m-1<0 we have to take care of the Associated Legendre Function
        {
            if (-m+1<=l+1)
            {
                term2 = (l-m+1)*(l+2-m)
                        *gsl_sf_legendre_Plm(l+1,1-m,x)*pow(-1.,1-m)*gsl_sf_fact(l+1-1+m)/gsl_sf_fact(l+1+1-m);
            }
            if (-m+1>l+1)
            {
                term2 =0.;
            }
        }
        if ((m-1>=0))
        {
            if (m-1<=l+1)
            {
                term2 = (l-m+1)*(l-m+2)*gsl_sf_legendre_Plm(l+1,m-1,x) ;
            }
            if (m-1>l+1)
            {
                term2 = 0.;
            }
        }
        return Clm*(term1+term2);
    }

    else
    {
    	return 0.;
    }

}

double DRnl(int n, int l, double r)
{

    if(l<n-1)
    {
    	double xn = 2.*r/n ;
    	double N1 = 2.* sqrt(gsl_sf_fact(n-l-1)) ;
    	double N2 = pow(n,2.)* sqrt(gsl_sf_fact(n+l)) ;
    	double Nnl = N1/N2 ;
    	double term1 = ((l/r)-(1./n))*gsl_sf_hydrogenicR(n, l, 1., r) ;
    	double term2 = pow(xn,l)*exp(-r/n)*gsl_sf_laguerre_n(n-l-2, 2*l+2, xn) ;
    	return  term1 - Nnl* term2*2./n ;
    }
    else
    {
    	return 0.;
    }

    /*
    COMMENTS :  Derivative of the radial hydrogenic wave functions. Note that r is not 0 !
     */
}

double Spherical(int l, int m, double x1)
{

    if(m<0)
    {
        return  gsl_sf_legendre_sphPlm(l,-m,x1)
        			*pow(-1,-m)*gsl_sf_fact(l+m)/gsl_sf_fact(l-m);
    }
    else {

        return gsl_sf_legendre_sphPlm(l, m,x1);
    }

}

