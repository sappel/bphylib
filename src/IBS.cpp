#include<string>
#include<complex>
#include<stdio.h>
#include<iostream>

using namespace std;

#include "PhyConstants.h"


/*!
Intrabeam scattering rates from Jie Wei
*/

double ibs_weifun(double x)
{
   double chifun;	
   if(x < 1.0) 
     chifun=atan(sqrt((1.0-x)/x))/sqrt(x*(1.0-x));
   else  
     chifun=atanh(sqrt((x-1.0)/x))/sqrt(x*(x-1.0));
   return ((1.0+2.0*x)*chifun-3.0)/(1.0-x);
}

void ibs_rates_wei(double rates[], double rms_emittance_h, double rms_emittance_v, double rms_mom_spread,
                       double beta_h, double beta_v, double D, double current, double gamma0,
		       double beta0, double Z, double A)
  {     
   double rc=qe*qe/(4.0*PI*eps0*mp*clight*clight);
   double dconst=D*rms_mom_spread/sqrt(beta_h*rms_emittance_h+pow(D*rms_mom_spread,2));
   double aconst=beta_h*dconst/(D*gamma0);
   double bconst=aconst*sqrt(beta_v*rms_emittance_h/(rms_emittance_v*beta_h));
   double A0=2.0*sqrt(PI)/(32.0*PI*PI)*rc*rc*current*pow(Z,3)/
    (qe*A*A*rms_emittance_h*rms_emittance_v*rms_mom_spread*pow(gamma0*beta0,4));
   double coul_log=10.0;
   double nbeam=2.0;   // 2=unbunched 1=bunched
   double xx=(aconst*aconst+bconst*bconst)/2.0;
   double Ffunc=ibs_weifun(xx);

   double rate_l= 4.0*PI*A0*coul_log*Ffunc*nbeam*(1.0-dconst*dconst);
   double rate_h= 4.0*PI*A0*coul_log*Ffunc*(dconst*dconst-0.5*aconst*aconst);
   double rate_v= 4.0*PI*A0*coul_log*Ffunc*(-0.5*bconst*bconst);
   rates[0]=rate_l;
   rates[1]=rate_h;
   rates[2]=rate_v;
  }


/*!
Plasma IBS model
*/

double plasma_long_diffusion(double rms_emittance, double rms_radius, double current, 
                             double gamma0,double beta0, double Z, double A)
{
 double ri=qe*qe*Z*Z/(4.0*PI*eps0*A*mp*clight*clight);
 double coul_log=20.0;
 return sqrt(PI)*ri*ri*current*coul_log/(4.0*pow(beta0*gamma0,3)*beta0*qe*Z
        *rms_emittance*rms_radius);
}


double cooling_betaf(double dp, double rate0, double beta0, double veff)
{
  double dp_eff=veff/(beta0*clight);
  return rate0; //*pow(dp_eff,3)/pow(dp_eff*dp_eff+dp*dp,1.5);
}

