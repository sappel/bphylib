// Energy loss (straggeling) due to an internal target

extern "C" 
{
	 extern float ranlan_(float*);
	 extern float dinlan_(float*);
	 extern float dinvav_(float*);
	 extern float disvav_(float*,int*);
	 extern void coedin_(float*,float*,const int*);
	 extern void coedis_(float*,float*,int*,const int*); 
}

class target {
	
    double xi, kappa, Emax, It, E_rms_2, delta_rms_2, E0, beta0;
    
public:

 target(double Zt, double At, double nt, double EI, double Zp, double Ap, double gamma0, double circum);
 double get_xi() { return xi; }
 double get_kappa() { return kappa; }
 double get_Emax() { return Emax; }
 double get_delta_rms_2() { return delta_rms_2; }
 double get_Erms() { return sqrt(E_rms_2); }
 double energyloss(double f0,double dt);
 double straggling_gauss(double f0, double dt, long* d);
 double straggling_landau(double f0, double dt, long *d);
 double straggling_vavilov(double f0, double dt, long *d);
 double straggling_single(double f0, double dt, long *d);
};
