void set_rf_single(Grid1D& rf, Grid1D& Yrf,double harmonic);
void set_rf_dual(Grid1D& rf, Grid1D& Yrf,double harmonic);
void set_rfamplitude(double voltage0, Grid1D& rf);
double match_elliptic(Grid1D& rf, Grid1D& Yrf, double zm2, double mom_spread, 
	double Zs, double Np, double beta0, double eta0, double charge, double mass);
double match_gauss(Grid1D& rf, Grid1D& Yrf, double zm2, double mom_spread, 
	double Zs, double Np, double beta0, double eta0, double charge, double mass);
double match_elliptic_with_hole(Grid1D& rf, Grid1D& Yrf, double zm2, double mom_spread, 
	double Zs, double Np, double beta0, double eta0, double charge, double mass, double zh);
void cut_distribution(Grid1D& Yrf, double voltage0, double charge, double meff, double zcut);
void matchPIC(vektor& x, vektor& dp, double beta0, double eta0, long* d);

