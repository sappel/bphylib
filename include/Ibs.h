// intrabeam scattering diffusion/rates

void ibs_rates_wei(double rates[], double rms_emittance_h,double rms_emittance_v,double rms_mom_spread,
                       double beta_h,double beta_v,double D,double current,double gamma0,
		       double beta0, double Z, double A);
double plasma_long_diffusion(double rms_emittance, double rms_radius, double current, 
                             double gamma0,double beta0, double Z, double A);

			     
// cooling force

double cooling_betaf(double dp, double rate0, double beta0, double veff);


