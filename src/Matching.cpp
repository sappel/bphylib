#include "bphylib.h"
#include "mathtools.h"

//! distribution function 
Grid2D distf;

/*!
Single rf: Sets the rf vector rf and the rf potential vector Yrf.
*/
void set_rf_single(Grid1D& rf, Grid1D& Yrf, double harmonic)
{
 int j;
 int dimx=rf.get_size();
 double dx=rf.get_dz();
 double circum=dimx*dx;
 for(j=0; j<dimx; j++)
  {
   rf[j]=-sin(2.0*PI*harmonic*rf.z[j]/circum);
   Yrf[j]=(circum/(2.0*PI*harmonic))*cos(2.0*PI*harmonic*rf.z[j]/circum)-(circum/(2.0*PI*harmonic));
  }
}


/*!
Dual rf: Sets the rf vector rf and the rf potential vector Yrf.
*/
void set_rf_dual(Grid1D& rf, Grid1D& Yrf, double harmonic)
{
 int j;
 int dimx=rf.get_size();
 double dx=rf.get_dz();
 double circum=dimx*dx;
 for(j=0; j<dimx; j++)
  {
   rf[j]=-sin(2.0*PI*harmonic*rf.z[j]/circum)+0.5*sin(4.0*PI*harmonic*rf.z[j]/circum);
   Yrf[j]=(circum/(2.0*PI*harmonic))*cos(2.0*PI*harmonic*rf.z[j]/circum)-(circum/(2.0*PI*harmonic))
        -0.5*(circum/(4.0*PI*harmonic))*cos(4.0*PI*harmonic*rf.z[j]/circum)-0.5*(circum/(4.0*PI*harmonic));
  }
}

void set_rfamplitude(double voltage0, Grid1D& rf)
{
	int dimx=rf.get_size();
	for(int j=0; j<dimx; j++)
		rf[j]=voltage0*rf[j];
}


/*!
Returns a matched elliptic bunch distribution distf for input 
rf voltage form rf and rf potential Yrf. 
*/
double match_elliptic(Grid1D& rf, Grid1D& Yrf, double zm2, double mom_spread, 
	double Zs, double Np, double beta0, double eta0, double charge, double mass)
{
  int j,k; 
  double density0=0.0, f0=0.0;
  double dx=rf.get_dz();
  int dimx=rf.get_size();
  double delta_v=-eta0*beta0*clight*mom_spread; 
  double circum=dimx*dx;
  double gamma0=1.0/sqrt(1.0-beta0*beta0);
  double Y0=Yrf.get_grid_lin(zm2);
  double tem=0.0, um=0.0;
  
   // init distf
  int dimv=160;
  Grid2D distf_tmp(dimx,dimv,dx,4.0*delta_v/dimv);
  distf=distf_tmp;
   
  // normalize bunch profile 
  for(j=0; j<dimx; j++)
    {
     tem=Yrf.get_grid_lin(Yrf.z[j])-Y0;
     if( tem*Y0 < 0.0 ) um+=tem*dx;
    }
   
  // max. line density and peak dist. function 
  density0=-Np*Y0/um;
  f0=2.0*density0/(PI*delta_v);
 
  // space charge contribution 
  double voltage_s=charge*beta0*clight*circum/(2.0*PI)*Zs*Np/um;

  // total voltage amplitude:
  double voltage0=(mass*gamma0*circum*pow(delta_v,2))/
    (2.0*charge*eta0*Y0)+voltage_s;

  /*
  cout << "Voltage_sc:" << voltage_s << endl;
  cout << "Voltage:" << voltage0 << endl;
  cout << "f_s:" << sqrt(fabs(eta0)*voltage0*charge
			/(2.0*PI*pow(circum,2)*mass*gamma0)) << endl;
  cout << "Sigma:" << voltage_s/voltage0/(1.0-voltage_s/voltage0) << endl;	
  */

  // init distribution function:
  for(j=0; j<dimx; j++)
      for(k=0; k<dimv; k++)
        {	
         tem=1.0-Yrf[j]/Y0-pow(distf.y[k]/delta_v,2);
         if ( tem > 0.0 ) distf[j*dimv+k]=f0*sqrt(tem);
         else tem=0.0;
        } 

  return voltage0;  
}


/*!
  nonlinear matched bunch with elliptic distribution 
  plus matched hole.
  Input parameter
  zh : hole half-length
  Output parameter
  voltage0: matched voltage amplitude
*/


double match_elliptic_with_hole(Grid1D& rf, Grid1D& Yrf, double zm2, double mom_spread, 
	double Zs, double Np, double beta0, double eta0, double charge, double mass, double zh)
{
  int j,k;  
  double density0, f0, fh, delta_v_h;
  double dx=rf.get_dz();
  int dimx=rf.get_size();
  double delta_v=-eta0*beta0*clight*mom_spread;
  double circum=dimx*dx;
  double gamma0=1.0/sqrt(1.0-beta0*beta0);
  double YRFm=Yrf.get_grid_lin(zm2);
  double YRFh=Yrf.get_grid_lin(zh);
  double Yfunc, Ym, Yh;
  double tem=0.0, um=0.0, uh=0.0;

  // init distf
  int dimv=160;
  Grid2D distf_tmp(dimx,dimv,dx,4.0*delta_v/dimv);
  distf=distf_tmp;

  for(j=0; j<dimx; j++)
    {
     tem=Yrf.get_grid_lin(distf.x[j])-YRFm;
     if( tem*YRFm < 0.0 && distf.x[j] <= zm2 )
      um+=tem*dx;
    }
   
   for(j=0; j<dimx; j++)
    {
     tem=Yrf.get_grid_lin(distf.x[j])-YRFh;
     if( tem*YRFh < 0.0 && distf.x[j] <= zh )
      uh+=tem*dx;
    }

   // flat line density:
   density0=Np/(um-uh)*(YRFh-YRFm);

  // space charge contribution without hole
  double Npnew=Np/(1.0-uh/um);
  double voltage_s=charge*beta0*clight*circum/(2.0*PI)*Zs*Npnew/um;

  // space charge contribution with hole
  voltage_s*=(1.0-YRFh/YRFm);

  // total voltage amplitude:
  double voltage0=(mass*gamma0*circum*pow(delta_v,2))/
    (2.0*charge*eta0*YRFm)+voltage_s;

  cout << "Voltage_sc:" << voltage_s << endl;
  cout << "Voltage:" << voltage0 << endl;
  cout << "f_s:" << sqrt(fabs(eta0)*voltage0*charge
			/(2.0*PI*pow(circum,2)*mass*gamma0)) << endl;

  double YmdYh=(1.0-voltage_s/voltage0)*YRFm/YRFh+voltage_s/voltage0;

  // init distribution without hole:

  f0=2.0/PI*density0/delta_v*1.0/(1.0-1.0/YmdYh);
  Ym=(1.0-voltage_s/voltage0)*YRFm+voltage_s/voltage0*YRFh;
  Yh=YRFh;

  for(j=0; j<dimx; j++)
      for(k=0; k<dimv; k++)
        {
	 if( distf.x[j] > zh || distf.x[j] < -zh )
	   Yfunc=(1.0-voltage_s/voltage0)*Yrf.get_grid_lin(distf.x[j])
	   +voltage_s/voltage0*Yrf.get_grid_lin(zh);
         if( distf.x[j]<= zh && distf.x[j] >= -zh )
	   Yfunc=Yrf.get_grid_lin(distf.x[j]);
         tem=1.0-Yfunc/Ym-pow(distf.y[k]/delta_v,2);
         if ( tem > 0.0 && distf.x[j] <= zm2 ) distf[j*dimv+k]=f0*sqrt(tem);
         else tem=0.0;
        } 

  
  // now the hole:

  // determine delta_v_h:
  tem=2.0*charge*eta0*voltage0/(circum*mass*gamma0)*YRFh;
  //if ( tem >= 0.0)
  delta_v_h=sqrt(abs(tem));
  //else { cout << "delta_v_hole < 0.0" << endl; exit(0); }
  if (delta_v_h >0.0)
    fh=2.0/PI*density0/delta_v_h*1.0/(YmdYh-1.0);
  

  // substract hole
   for(j=0; j<dimx; j++)
      for(k=0; k<dimv; k++)
        {
         tem=1.0-Yrf.get_grid_lin(distf.x[j])/Yh-pow(distf.y[k]/delta_v_h,2);
         if ( tem > 0.0 && distf.x[j] <= zh ) distf[j*dimv+k]-=fh*sqrt(tem);
         else tem=0.0;
        }
 
  cout << "Sigma:" << voltage_s/voltage0/(1.0-voltage_s/voltage0) << endl;	
	
  return voltage0;  

}

//! Matched Gaussian bunch for equivalent elliptic distribution
/*!
Returns a matched Gaussian bunch distribution distf for input 
rf voltage form rf and rf potential Yrf. 
*/
double match_gauss(Grid1D& rf, Grid1D& Yrf, double zm2, double mom_spread, 
	double Zs, double Np, double beta0, double eta0, double charge, double mass)
{
	
	// init elliptic distf 
	double V0=match_elliptic(rf,Yrf,zm2,mom_spread,Zs,Np,beta0,eta0,charge,mass);
	
	int j,k;
	double dx=rf.get_dz();
	double dv=distf.get_dy();
	int dimx=distf.get_NX();
	int dimv=distf.get_NY();
	double circum=dimx*dx;
        double delta_v=-eta0*beta0*clight*mom_spread;
	double gamma0=1.0/sqrt(1.0-beta0*beta0);
	
	vektor ldy_0(dimx); 
	vektor ldy_1(dimx);
	vektor ldy_2(dimx);
	double tem=0.0;
	double mred=mass*gamma0/fabs(eta0);
	
	// calculate rms velocity of elliptic distribution
	double v2rms = 0.0;
	for(int k=0; k<dimv; k++)
		for(int j=0; j<dimx; j++)
			v2rms+=pow(distf.y[k],2)*distf[j*dimv+k]*dx*dv/Np;
	double vrms=sqrt(v2rms); // vm/sqrt(5.0); 
	
	// line density from elliptic f0:

	for(j=0; j<dimx; j++)
	{
		ldy_0[j]=0.0;
		for(k=0; k<dimv; k++)
		   ldy_0[j]+=distf[j*dimv+k]*dv;
		ldy_1[j]=ldy_0[j];
	}
	
	// line density iteration:
	
	for (int it=0; it<30; it++)
	{	
	
	// new line density:
	
	for(j=0; j<dimx; j++)
		{
		 ldy_2[j]=exp(charge*V0*Yrf.get_grid_lin(distf.x[j])/(mred*circum*pow(vrms,2)))
			*exp(-pow(charge,2)*beta0*clight/(mred*2.0*PI*pow(vrms,2))*Zs*((1.0/3.0)*ldy_0[j]+(2.0/3.0)*ldy_1[j]));
		}
	
	// normalize new line density:
	
	tem = 0.0;
	for(int j=0; j<dimx; j++) 
			tem+=ldy_2[j]*dx;
	for(int j=0; j<dimx; j++)
			ldy_2[j]*=Np/tem;
	
	// update old line densities
	
	for(int j=0; j<dimx; j++)
	  {	
			ldy_0[j]=ldy_1[j];
			ldy_1[j]=ldy_2[j];
	  }
	 
        }
	
       // new f with space charge:
       
       for(j=0; j<dimx; j++)
		for(k=0; k<dimv; k++)
			distf[j*dimv+k]=ldy_1[j]/(sqrt(2.0*PI)*vrms)*exp(-0.5*pow((distf.y[k])/vrms,2));
			
      return V0;  
}


void cut_distribution(Grid1D& Yrf, double voltage0, double charge, double meff, double zcut)
{
 int i,j;
 int dimx=distf.get_NX();
 int dimv=distf.get_NY();
 double dx=Yrf.get_dz();
 double circum=dimx*dx;
 double Ycut=Yrf.get_grid_lin(zcut);
 
 for(i=0; i<dimx; i++)
   for(j=0; j<dimv; j++)
     	if( (meff*circum/(2.0*charge*voltage0))*pow(distf.y[j],2)-Yrf.get_grid_lin(distf.x[i]) > -Ycut )  
           distf(i,j)=0.0;
 
}


//! Interface to PIC codes. Random number generator.

void matchPIC(vektor& x, vektor& dp, double beta0, double eta0, long* d)
{	
  double xran, vran, fran;
  double fmax=0.0, ftem;
  int dimx=distf.get_NX();
  int dimv=distf.get_NY();
  double xmax=-distf.get_beginx();
  double vmax=-distf.get_beginy();
  
  // find fmax:
  for(int j=0; j<dimx; j++)
    for(int l=0; l<dimv; l++)
      {
       ftem=distf(j,l);
       if (ftem > fmax) fmax=ftem;
      } 

   for (int j=0 ; j < x.size() ; j++)
	 {
         do { xran = (2.0 * ran1(d) - 1.0)*xmax ;  
	     vran = (2.0 * ran1(d) - 1.0)*vmax ;
         fran = ran1(d)*fmax;              
         } while (fran > distf.Grid2PIC(xran,vran)) ;
               
          dp[j]   =  -vran/(eta0*beta0*clight);
          x[j] = xran;
      }
}



