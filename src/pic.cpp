#include "bphylib.h"
#include "pic.h"
#include "mathtools.h"


// T0 in eV

void PIC::set_const_density_xy(double boundary, double T0, long *d1)
{ 
 double u,v,w;	
 double ut=clight*sqrt(pow(fabs(charge_p)*T0/(mass_p*clight*clight)+1.0,2)-1.0);
 for (int j=0 ; j < part.size(); j++)
  {
   	do {
        u = boundary*(2.0*ran1(d1)-1.0); 
     	v = boundary*(2.0*ran1(d1)-1.0);
	   } while (sqrt(pow(u,2)+pow(v,2)) > boundary);
    		
	part[j].x=u;  	
	part[j].y=v;
	part[j].z=length * ran1(d1);
	
	do { u = 2.0 * ran1(d1) - 1.0 ;  
	     v = 2.0 * ran1(d1) - 1.0 ;              
	     w = u*u + v*v ; }
	     while (w >= 1.0) ;

		part[j].ux = ut * u * sqrt(-2.0*log(w)/w) ;
	    part[j].uy = ut * v * sqrt(-2.0*log(w)/w) ;
		part[j].uz = 0.0; 
	 
  }	
}

/*!
inject beam
*/


void PIC::inject_beam(double boundary, double vel, double angle, double ndensity, double Tempxy, double Tempz, double dt, long *d1)
{
 double Ninj=vel*ndensity*PI*pow(boundary,2)*dt/macroN; 
 double gamma=1.0/sqrt(1.0-pow(vel/clight,2));
 double utherm=clight*sqrt(pow(charge_p*Tempxy/(mass_p*clight*clight)+1.0,2)-1.0);
 double r_amp, phase;
 double u, v, w;
 for(int j=0; j<Ninj; j++) 
   { 
	Particle part_inj;
	r_amp=boundary*sqrt(ran1(d1));
	phase=2.0*PI*ran1(d1);	
	part_inj.x=r_amp*cos(phase);
	part_inj.y=r_amp*sin(phase);
	part_inj.z=-0.5*length;
	
	do { u = 2.0 * ran1(d1) - 1.0 ;  
	     v = 2.0 * ran1(d1) - 1.0 ;              
	     w = u*u + v*v ; }
	     while (w >= 1.0) ;
	
	part_inj.ux = part_inj.x/r_amp*vel*angle + utherm * u * sqrt(-2.0*log(w)/w) ;
	part_inj.uy = part_inj.y/r_amp*vel*angle + utherm * v * sqrt(-2.0*log(w)/w) ;
	part_inj.uz = gamma*vel;
    part.push_back(part_inj);
   }	
}



void PIC::densityXY(Grid2D& target)
{
  long j;
  target.reset();
  for(j=0; j<part.size(); j++)
	  target.Pic2Grid(macroN*charge_p/length,part[j].x,part[j].y);	
}


void PIC::densityRZ(Grid2D_rz& target)
{
  long j;
  target.reset();
  for(j=0; j<part.size(); j++)
	    target.Pic2Grid(macroN*charge_p,sqrt(pow(part[j].x,2)+pow(part[j].y,2)),part[j].z);
}

void PIC::get_potential_xy(Grid2D& pote, Grid2D& rhoe, Greenfb& gf)
{
	pote=rhoe;
	poisson_xy(pote,gf);		
}


void PIC::get_potential_rz(Grid2D_rz& pote, Grid2D_rz& rhoe)
{
	pote=rhoe;
	poisson_rz(pote);
}

/*!
coordinate shift with absorbing boundary
*/

void PIC::shift_xyz(double dt, double boundary)
{
  long j;
  double gamma, u2;
  for(j=0; j<part.size(); j++)
   {
	u2=pow(part[j].uy,2)+pow(part[j].ux,2)+pow(part[j].uz,2);
	gamma=sqrt(1.0+u2/pow(clight,2));
	part[j].x+=part[j].ux*dt/gamma;
	part[j].y+=part[j].uy*dt/gamma;	
	part[j].z+=part[j].uz*dt/gamma;
	if ( sqrt(pow(part[j].x,2)+pow(part[j].y,2)) > boundary || part[j].z > 0.5*length || part[j].z < -0.5*length )
	  part.erase(part.begin()+j);
	  
   }
}

/*!
coordinate shift with SEY at boundary
*/

double PIC::shift_xyz(double dt, double boundary, double delta_max, double Emax, double E0, long *d1)
{
  long j;
  int ref=0, sec=0;
  double gamma, u2, r2;
  double T0=2.0, phi_emit=0.0, energy=0.0;
  double ut=clight*sqrt(pow(fabs(charge_p)*T0/(mass_p*clight*clight)+1.0,2)-1.0);
  vector<Particle> part_new; 
  for(j=0; j<part.size(); j++)
   {
	u2=pow(part[j].uy,2)+pow(part[j].ux,2)+pow(part[j].uz,2);
	gamma=sqrt(1.0+u2/pow(clight,2));
	part[j].x+=part[j].ux*dt/gamma;
	part[j].y+=part[j].uy*dt/gamma;
    r2=pow(part[j].y,2)+pow(part[j].x,2);
	if ( sqrt(r2) > boundary ) 
	  {
	   ref=reflection(get_energy(j),E0,d1);
	   if (delta_max > 0.0) sec=secondary(get_energy(j),delta_max,Emax,d1);	
	   if ( ref == 1)
	    	{ 
		     //phi_emit=0.5*PI*(2.0*ran1(d1)-1.0);
			 //part[j].ux=-ut*(cos(phi_emit)*part[j].x/sqrt(r2)-sin(phi_emit)*part[j].y/sqrt(r2));
			 //part[j].uy=-ut*(sin(phi_emit)*part[j].y/sqrt(r2)+cos(phi_emit)*part[j].y/sqrt(r2));
		     part[j].ux=-part[j].ux;
		     part[j].uy=-part[j].uy;
		     part[j].x=0.99*boundary*part[j].x/sqrt(r2);  
			 part[j].y=0.99*boundary*part[j].y/sqrt(r2);
		    } 
	   else { 
		     if (sec == 0) 
		       { 	
				energy+=get_energy(j);
		        part.erase(part.begin()+j);
		       }
	         else if (sec == 1)
	          { 
		       energy+=get_energy(j);
			   phi_emit=0.1*PI*(2.0*ran1(d1)-1.0);
			   part[j].ux=-ut*(cos(phi_emit)*part[j].x/sqrt(r2)-sin(phi_emit)*part[j].y/sqrt(r2));
			   part[j].uy=-ut*(sin(phi_emit)*part[j].y/sqrt(r2)+cos(phi_emit)*part[j].y/sqrt(r2));
		       part[j].x=0.99*boundary*part[j].x/sqrt(r2);  
			   part[j].y=0.99*boundary*part[j].y/sqrt(r2);  
		      }
	         else if (sec == 2)
	          { 
		       energy+=get_energy(j);
		       Particle part_inj;
		       phi_emit=0.1*PI*(2.0*ran1(d1)-1.0);
		       part[j].ux=-ut*(cos(phi_emit)*part[j].x/sqrt(r2)-sin(phi_emit)*part[j].y/sqrt(r2));
		       part[j].uy=-ut*(sin(phi_emit)*part[j].y/sqrt(r2)+cos(phi_emit)*part[j].y/sqrt(r2));
		       phi_emit=0.1*PI*(2.0*ran1(d1)-1.0);
		       part_inj.ux = -ut*(cos(phi_emit)*part[j].x/sqrt(r2)-sin(phi_emit)*part[j].y/sqrt(r2));
		   	   part_inj.uy = -ut*(sin(phi_emit)*part[j].y/sqrt(r2)+cos(phi_emit)*part[j].y/sqrt(r2));
		       part_inj.uz = 0.0;
		       part[j].x=0.99*boundary*part[j].x/sqrt(r2);  
			   part[j].y=0.99*boundary*part[j].y/sqrt(r2);  
		       part_inj.x = part[j].x;
			   part_inj.y = part[j].y; 
			   part_inj.z = part[j].z;
			   part_new.push_back(part_inj);
		     }
		    else if (sec == 3)
	         { 
		       energy+=get_energy(j);
			   Particle part_inj, part_inj_2;
			   phi_emit=0.1*PI*(2.0*ran1(d1)-1.0);
		       part[j].ux=-ut*(cos(phi_emit)*part[j].x/sqrt(r2)-sin(phi_emit)*part[j].y/sqrt(r2));
		       part[j].uy=-ut*(sin(phi_emit)*part[j].y/sqrt(r2)+cos(phi_emit)*part[j].y/sqrt(r2));
		       phi_emit=0.1*PI*(2.0*ran1(d1)-1.0);
		       part_inj.ux = -ut*(cos(phi_emit)*part[j].x/sqrt(r2)-sin(phi_emit)*part[j].y/sqrt(r2));
		   	   part_inj.uy = -ut*(sin(phi_emit)*part[j].y/sqrt(r2)+cos(phi_emit)*part[j].y/sqrt(r2));
		       part_inj.uz = 0.0;
		       phi_emit=0.1*PI*(2.0*ran1(d1)-1.0);
			   part_inj_2.ux = -ut*(cos(phi_emit)*part[j].x/sqrt(r2)-sin(phi_emit)*part[j].y/sqrt(r2));
		   	   part_inj_2.uy = -ut*(sin(phi_emit)*part[j].y/sqrt(r2)+cos(phi_emit)*part[j].y/sqrt(r2));
		       part_inj_2.uz = 0.0;
		       part[j].x=0.99*boundary*part[j].x/sqrt(r2);  
			   part[j].y=0.99*boundary*part[j].y/sqrt(r2);  
		       part_inj.x = part[j].x;
			   part_inj.y = part[j].y; 
			   part_inj.z = part[j].z;
			   part_inj_2.x = part[j].x;
			   part_inj_2.y = part[j].y; 
			   part_inj_2.z = part[j].z;
			   part_new.push_back(part_inj);
			   part_new.push_back(part_inj_2);
		     }
			else {cout << "shift_xyz: n>4" << endl; exit(0);}
		    }
	  }	
   }
   for(j=0; j<part_new.size(); j++)
	    part.push_back(part_new[j]);
	return macroN*energy;
}



void PIC::shift_uxy(Grid2D& pot, double dt)
{
  long j;
  double Ex, Ey;
  for(j=0; j<part.size(); j++)
   {
	Ex=-pot.Grid2dx(part[j].x,part[j].y);
	Ey=-pot.Grid2dy(part[j].x,part[j].y);
	
	part[j].ux+=charge_p*Ex*dt/mass_p;
	part[j].uy+=charge_p*Ey*dt/mass_p;
   }
}


void PIC::shift_urz(Grid2D_rz& pot, double dt)
{
  long j;
  double Er, Ez, rpart;
  for(j=0; j<part.size(); j++)
   {
	Er=-pot.Grid2dr(sqrt(pow(part[j].x,2)+pow(part[j].y,2)),part[j].z);
	Ez=-pot.Grid2dz(sqrt(pow(part[j].x,2)+pow(part[j].y,2)),part[j].z);
	rpart=sqrt(pow(part[j].x,2)+pow(part[j].y,2));
	
	part[j].ux+=part[j].x/rpart*charge_p*Er*dt/mass_p;
	part[j].uy+=part[j].y/rpart*charge_p*Er*dt/mass_p;
	part[j].uz+=charge_p*Ez*dt/mass_p;
   }
}

/*!
Boris particle updater for a B-field in x,y.
The push is done in 3D cartesian coordinates.
*/


void PIC::Boris_uxy(Grid2D& Bx, Grid2D& By, double dt)
{
 long j;	
 double tx, ty, tabs2, Bxj, Byj, Babs, u2, gamma, ux0, uy0, uz0;	
 for(j=0; j<part.size(); j++)
   { 
    Bxj=Bx.Grid2PIC(part[j].x,part[j].y);
	Byj=By.Grid2PIC(part[j].x,part[j].y); 
	Babs=sqrt(pow(Bxj,2)+pow(Byj,2));
	u2=pow(part[j].uy,2)+pow(part[j].ux,2)+pow(part[j].uz,2);
	gamma=sqrt(1.0+u2/pow(clight,2));
	ux0=part[j].ux;uy0=part[j].uy;uz0=part[j].uz;
	if (Babs > 0.0) {
		tx=Bxj/Babs*tan(charge_p*dt*Babs/(2.0*gamma*mass_p));
		ty=Byj/Babs*tan(charge_p*dt*Babs/(2.0*gamma*mass_p)); }
	else {tx=0.0; ty=0.0;}
	
	tabs2=pow(tx,2)+pow(ty,2); 
	
	part[j].ux+=-part[j].uz*ty;
	part[j].uy+=part[j].uz*tx;
    part[j].uz+=part[j].ux*ty-part[j].uy*tx;
	
	part[j].ux=ux0-part[j].uz*ty*2.0/(1.0+tabs2);
	part[j].uy=uy0+part[j].uz*tx*2.0/(1.0+tabs2);
    part[j].uz=uz0+part[j].ux*ty*2.0/(1.0+tabs2)-part[j].uy*tx*2.0/(1.0+tabs2);
   }
}

/*!
Boris particle updater for a pure dipole B-field in y-direction.
The push is done in 3D cartesian coordinates.
*/


void PIC::Boris_uxy(double By, double dt)
{
 long j;	
 double ty, tabs2, u2, gamma, ux0, uy0, uz0;	
 for(j=0; j<part.size(); j++)
   { 
	u2=pow(part[j].uy,2)+pow(part[j].ux,2)+pow(part[j].uz,2);
	gamma=sqrt(1.0+u2/pow(clight,2));
	ux0=part[j].ux;uy0=part[j].uy;uz0=part[j].uz;
	ty=By/fabs(By)*tan(charge_p*dt*fabs(By)/(2.0*gamma*mass_p));
	
	tabs2=pow(ty,2); 
	
	part[j].ux+=-part[j].uz*ty;
    part[j].uz+=part[j].ux*ty;
	
	part[j].ux=ux0-part[j].uz*ty*2.0/(1.0+tabs2);
    part[j].uz=uz0+part[j].ux*ty*2.0/(1.0+tabs2);
   }
}

/*!
Boris particle updater for a B-field given in cylindrical coordinates.
The push is done in 3D cartesian coordinates.
*/

void PIC::Boris_urz(Grid2D_rz& Br, Grid2D_rz& Bz, double dt)
{
 long j;	
 double tx, ty, tz, tabs2, Bxj, Byj, Bzj, Babs, u2, gamma, ux0, uy0, uz0, rpart;	
 for(j=0; j<part.size(); j++)
   { 
	rpart=sqrt(pow(part[j].x,2)+pow(part[j].y,2));
    Bxj=part[j].x/rpart*Br.Grid2PIC(rpart,part[j].z);
	Byj=part[j].y/rpart*Br.Grid2PIC(rpart,part[j].z);
	Bzj=Bz.Grid2PIC(rpart,part[j].z);
	Babs=sqrt(pow(Bxj,2)+pow(Byj,2)+pow(Bzj,2));
	u2=pow(part[j].uy,2)+pow(part[j].ux,2)+pow(part[j].uz,2);
	gamma=sqrt(1.0+u2/pow(clight,2));
	ux0=part[j].ux;uy0=part[j].uy;uz0=part[j].uz;
	if (Babs > 0.0) {
		tx=Bxj/Babs*tan(charge_p*dt*Babs/(2.0*gamma*mass_p));
		ty=Byj/Babs*tan(charge_p*dt*Babs/(2.0*gamma*mass_p)); 
		tz=Bzj/Babs*tan(charge_p*dt*Babs/(2.0*gamma*mass_p)); 
		}
	else {tx=0.0; ty=0.0; tz=0.0;}
	
	tabs2=pow(tx,2)+pow(ty,2); 
	
	part[j].ux+=part[j].uy*tz-part[j].uz*ty;
	part[j].uy+=part[j].uz*tx-part[j].ux*tz;
    part[j].uz+=part[j].ux*ty-part[j].uy*tx;
	
	part[j].ux=ux0+(part[j].uy*tz-part[j].uz*ty)*2.0/(1.0+tabs2);
	part[j].uy=uy0+(part[j].uz*tx-part[j].ux*tz)*2.0/(1.0+tabs2);
    part[j].uz=uz0+(part[j].ux*ty-part[j].uy*tx)*2.0/(1.0+tabs2);
   }
}


double PIC::Tx_rms()
{
 long j, n=part.size();
 double tem1=0.0, tem2=0.0;
 for(j=0; j<n; j++)
   {
    tem1+=pow(get_vx(j),2);
    tem2+=get_vx(j);
   }
 return mass_p*pow(sqrt(tem1/part.size()-pow(tem2/part.size(),2)),2)/qe; 
}


double PIC::vth_xy()
{
 long j, n=part.size();
 double tem1=0.0, tem2=0.0;
 for(j=0; j<n; j++)
   {
    tem1+=pow(0.5*(get_vx(j)+get_vy(j)),2);
    tem2+=0.5*(get_vx(j)+get_vy(j));
   }
 return sqrt(tem1/part.size()-pow(tem2/part.size(),2)); 
}


double PIC::x_rms()
{
 long j, n=part.size();
 double tem1=0.0, tem2=0.0;
 for(j=0; j<n; j++)
   {
    tem1+=pow(part[j].x,2);
    tem2+=part[j].x;
   }
 return sqrt(tem1/part.size()-pow(tem2/part.size(),2)); 
}


double PIC::offset_x()
{
 long j, n=part.size();
 double tem1=0.0;
 for(j=0; j<n; j++)
    tem1+=part[j].x;  
 return tem1/part.size(); 
}

double PIC::offset_y()
{
 long j, n=part.size();
 double tem1=0.0;
 for(j=0; j<n; j++)
    tem1+=part[j].y;  
 return tem1/part.size(); 
}

double PIC::total_charge()
{
	return charge_p*macroN*part.size();   // charge/qe
}

// energy/length in eV/m

double PIC::total_kinetic_energy()
{
 long j, n=part.size();
 double tem1=0.0;
 for(j=0; j<n; j++)
	     tem1+=macroN*get_energy(j);
 return tem1;
}


double PIC::total_energy(Grid2D& pot)
{
	long j, n=part.size();
	double ekin=total_kinetic_energy();	
	double epot=0.0;
	for(j=0; j<n; j++)
		  epot+=macroN*pot.Grid2PIC(part[j].x,part[j].y);
	return ekin+epot;
}


int PIC::reflection(double energy, double E0, long *d)
{
	//double E0=150.0;  //  eV
	double delta=pow(sqrt(energy)-sqrt(energy+E0),2)/pow(sqrt(energy)+sqrt(energy+E0),2);
	if ( ran1(d) <= delta ) return 1;
	else return 0;
}


int PIC::secondary(double energy, double delta_max, double Emax, long *d)
{
	double s=1.35;
	double x=energy/Emax;
	double nmax=3.0;
	double deltas=delta_max*s*x/(s-1.0+pow(x,s));  // true secondaries
	double P0=1.0*pow(1.0-deltas/nmax,3.0); 
	double P1=3.0*deltas/nmax*pow(1.0-deltas/nmax,2.0);	
	double P2=3.0*pow(deltas/nmax,2.0)*(1.0-deltas/nmax);
	double P3=1.0*pow(deltas/nmax,3.0);
	double rannum=ran1(d);
	if (rannum < P0) return 0;
	else if (rannum < P0+P1) return 1;
	else if (rannum < P0+P1+P2) return 2;
	else if (rannum < P0+P1+P2+P3) return 3;
	else return 0;
}



