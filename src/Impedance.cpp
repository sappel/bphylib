#include "bphylib.h"
#include "mathtools.h"
#include <algorithm>


BBImpedance::BBImpedance(int NZ)
{
	Z=vector<komplex>(NZ/2);
	for(int j=0; j<NZ/2; j++) Z[j]=0.0;
}   


NBImpedance::NBImpedance(int h, int N_t)
: time_window(2*N_t)
{
	harmonic=h;
	Z=vector<komplex>(N_t);
	for(int j=0; j<N_t; j++) Z[j]=0.0;
}   


void BBImpedance::cavity(double Rs, double Q, double harmonic, double omega0, double off_res)
{
	int j;
	komplex i(0.0,1.0);
	double omega_h=harmonic*omega0+off_res;
	for(j=1; j<Z.size(); j++)  
		Z[j]+=Rs/(1.0+i*Q*((omega0*j)/omega_h-omega_h/(omega0*j)));
}


void BBImpedance::MA_cavity(double Zabs,double phase0, double phase1, double fmax, double f0)
{
 double fj;		
 for(int j=1; j<Z.size(); j++)  
   { 
    fj=j*f0;	 
    Z[j]+=std::polar(Zabs,PI*(phase0+(phase1-phase0)*fj/fmax)/180.0);	
   }
}


void NBImpedance::cavity(double Rs, double Q, double omega0, double off_res, double dt)
{
	int N_t=time_window.size()/2;
	if( N_t==0 ) return;
	int j;
	komplex i(0.0,1.0);
	double omega_res=harmonic*omega0+off_res;
	for(j=0; j<Z.size(); j++)
	{
		double N_t=time_window.size()/2.0;
		double omega=omega0*harmonic+j*2.0*PI/(N_t*dt)-PI/dt;
		if( omega != 0.0 )
			Z[j]=Rs/(1.0+i*Q*(omega/omega_res-omega_res/omega));
	}
}


void BBImpedance::spacecharge(double Z0, double h_c)
{
	int j;
	komplex i(0,1.0);
	for(j=1; j<Z.size(); j++)  
		Z[j]+=-i*Z0*(double)(j)/(1.0+pow((double)j/h_c,2));
}


// `beam loading' field ( returns amplitude, efield, potential)
//  input current, dt, dz, circum, e_kin:


void BBImpedance::InducedField(Grid1D& efield, Grid1D& line_current_density)
{
	double bl_real, bl_imag; 
	double dz=efield.get_dz();
	int n=efield.get_size();
	double circum=dz*n;
	vektor current(n);
	for(int j=0; j<n; j++) current[j]=line_current_density[j]; 

	realft(current,1);

	for(int h=0; h < n/2-1; h++)
	{

		komplex temc(current[2*h],current[2*h+1]);
		komplex temc2=-temc*Z[h];

		bl_real = temc2.real()/circum; 
		bl_imag = temc2.imag()/circum;         

		current[2*h]=bl_real;
		current[2*h+1]=bl_imag;
	}

	realft(current,-1);

	for(int j=0; j<n; j++)
		efield[j]=current[j];    

	return;
}


void BBImpedance::InducedPotential(Grid1D& potential, Grid1D& line_current_density)
{
	komplex i(0.0,1.0);
	double bl_real=0.0, bl_imag=0.0; 
	double dz=potential.get_dz();
	int n=potential.get_size();
	double circum=dz*n;
	double R=circum/(2.0*PI);
	vektor current(n);
	for(int j=0; j<n; j++) current[j]=line_current_density[j]; 

	realft(current,1);

	for(int h=1; h < n/2-1; h++)
	{

		komplex temc(current[2*h],current[2*h+1]);
		komplex temc2=-i*temc*Z[h]*R/((double)h);

		bl_real = temc2.real(); 
		bl_imag = temc2.imag();         

		current[2*h]=bl_real;
		current[2*h+1]=bl_imag;
	}

	realft(current,-1);

	for(int j=0; j<n; j++)
		potential[j]=current[j];    

	return;
}


void NBImpedance::InducedField(Grid1D& efield, Grid1D& line_current_density)
{
	int n=efield.get_size();
	int N_t=time_window.size()/2;
	
	if ( N_t == 0) {
		for(int j=0;j<n;j++) efield[j]=0.0;
		return;
	}
	
	static int t=0;
	double bl_amp, bl_real, bl_imag, bl_phase; 
	double dz=efield.get_dz();
	
	double circum=dz*n;
	double bucket_h = circum/(double)harmonic;
	vektor current(n);
	for(int j=0; j<n; j++) current[j]=line_current_density[j]; 

	realft(current,1);

	if(t==0)
	{
		for(int j=0;j<2*N_t;j+=2)
		{
			time_window[j]=current[2*(int)harmonic];
			time_window[j+1]=current[2*(int)harmonic+1];
		}
		t=1;
	}

	time_window.pop_front();
	time_window.pop_front();
	time_window.push_back(current[2*(int)harmonic]);
	time_window.push_back(current[2*(int)harmonic+1]);

	vektor tem(2*N_t);
	for(int j=0;j<2*N_t;j++)
		tem[j]=time_window[j]; 
	four1(tem,1);

	for(int j=0;j<=N_t/2;j++)
	{ 
		komplex temc(tem[2*j],tem[2*j+1]);
		komplex temc2(0.0,0.0);
		temc2=-temc*Z[N_t/2-j];
		tem[2*j]=temc2.real(); 
		tem[2*j+1]=temc2.imag(); 
	}

	for(int j=1;j<N_t/2;j++)
	{ 
		komplex temc(tem[2*N_t-2*j-3],tem[2*N_t-2*j-2]);
		komplex temc2(0.0,0.0);
		temc2=-temc*Z[j+N_t/2];
		tem[2*N_t-2*j-3]=temc2.real(); 
		tem[2*N_t-2*j-2]=temc2.imag(); 
	}

	four1(tem,-1);

	bl_real=tem[2*N_t-2]/circum;
	bl_imag=tem[2*N_t-1]/circum;
	bl_amp   = sqrt(bl_real*bl_real+bl_imag*bl_imag) ;
	bl_phase = atan2(bl_imag,bl_real) ;

	for(int j=0;j<n;j++)
		efield[j]=bl_amp*cos(2.0*PI/bucket_h*(j+0.5)*dz-bl_phase);

	return;

}


