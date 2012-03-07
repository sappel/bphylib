#include "bphylib.h"
#include "mathtools.h"

/*!
Initialzation of the target parameters.
Target input parameters: charge state, mass number, area density [m^-2]
Beam input parameters: charge state, relativistic gamma.
*/

target::target(double Zt, double At, double nt, double EI, double Zp, double Ap, double gamma0, double circum)
{
 It=EI;	
 beta0  = sqrt((gamma0*gamma0-1.0)/(gamma0*gamma0)) ;
 xi=1.0e6*0.1535*pow(Zp,2)*Zt/pow(beta0,2)*mp*1.0e3*nt*1.0e-4; // in eV
 Emax=1.0/qe*2.0*me*pow(clight*beta0*gamma0,2)/(1.0+2.0*gamma0*me/(Ap*mp)+pow(me/(Ap*mp),2));  // in eV
 kappa=xi/Emax;

 double f0 = beta0*clight/circum;
 E_rms_2=xi*Emax*(1.0-0.5*pow(beta0,2)); // in eV
 double Ekin=Ap*(gamma0-1.0)*mp*pow(clight,2)/qe;  // in eV
 delta_rms_2=pow(gamma0/(gamma0+1.0),2)*E_rms_2/pow(Ekin,2); 
 E0=gamma0*mp*Ap*pow(beta0*clight,2)/qe; // in eV
}

/*!
Mean energy loss from Bethe-Bloch.
Target input parameter: Ionization energy [eV].
Beam input: gamma
Return value: mean energy loss [eV]
*/

double target::energyloss(double f0, double dt) 
{
 //return 2.0*xi*f0*dt*log(Emax/It-pow(beta0,2));
 return -xi*f0*dt*log(Emax/It);
}


double target::straggling_gauss(double f0, double dt, long* d)
{
 double Rz=gasdev(d); 
 return Rz*sqrt(delta_rms_2*f0*dt);  
}


double target::straggling_landau(double f0, double dt, long *d)
{
 float rantmp;
 double lambda_mean=-0.422-pow(beta0,2)-log(f0*dt*kappa);
 double lambda_max=0.60715+1.1934*lambda_mean+(0.67794+0.052382*lambda_mean)*exp(0.94753+0.74442*lambda_mean);
 double ranl; 
 do {
     rantmp=ran1(d);
     ranl=ranlan_(&rantmp);
    } while (ranl > lambda_max );  
 return xi*f0*dt*(ranl+0.422+pow(beta0,2)+log(kappa*f0*dt))/(pow(beta0,2)*E0);

}


double target::straggling_vavilov(double f0, double dt, long *d)
{
 float tmp=ran1(d);
 float rka=dt*f0*kappa;
 if(rka < 0.01 || rka > 10.0)
	 { cout << "Vavilov: kappa=" << rka << endl; exit(0); }
 float bet2=pow(beta0,2);
 const int j=-1;
 coedin_(&rka,&bet2,&j);
 return xi*dt*f0*dinvav_(&tmp)/E0;	
}


double target::straggling_single(double f0, double dt, long *d)
{
 float rantmp=ran1(d);
 if (xi*dt*f0/It > 1.0 ) { cout << "target: xi/I >1 !" << endl; exit(0); }
 if (rantmp > xi*dt*f0/It) return 0.0; 
 double denergy=1.0/(1.0/It-rantmp/(xi*dt*f0));
 if (denergy > Emax) return 0.0; //Emax/(pow(beta0,2)*E0);
 return denergy/(pow(beta0,2)*E0);
}



