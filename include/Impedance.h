// calculates the induced field 
// resulting from a complex impedance Z  

//! Broadband impedance:

class BBImpedance {

	vector<komplex> Z; 

public:  

	BBImpedance() {}
	BBImpedance(int NZ);
	~BBImpedance() {}

	double get_ZR(int j) { return Z[j].real(); }
	double get_ZI(int j) { return Z[j].imag(); }

	void cavity(double Rs, double Q, double harmonic, double omega0, double off_res);
	void MA_cavity(double Zabs, double phase0, double phase1, double fmax, double f0);
	void spacecharge(double Z0, double h_c);

	void InducedField(Grid1D& efield, Grid1D& line_current_density); 
	void InducedPotential(Grid1D& potential, Grid1D& line_current_density);

};

//! Narrowband impedance:

class NBImpedance {

	double harmonic; 
	vector<komplex> Z;
	deque<double> time_window;    // time window 

public:  

	NBImpedance() {}
	NBImpedance(int h, int N_t);
	~NBImpedance() {}

	double get_ZR(int j) { return Z[j].real(); }
	double get_ZI(int j) { return Z[j].imag(); }

	void cavity(double Rs, double Q, double omega0, double off_res, double dt);

	void InducedField(Grid1D& efield, Grid1D& line_current_density); 

};


