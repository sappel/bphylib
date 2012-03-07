//! 3D particle structure.
struct Particle {
 //!  x, y, z positions. x,y: Transverse coordinates. z: longitudinal 
 double x, y, z; 
 //! x, y, z velocities.
 double ux, uy, uz;   
};


//! Particle container

class PIC {
 
protected:

 vector<Particle> part;
 double macroN;   //!< number of physical particles per macro particle
 double charge_p, mass_p;   //!< particle charge and charge numbers
 double length;  //!< length of the missing dimension;

public:
 
 PIC() {}
 PIC(double macro, double q, double m, double len, int n) : part(n)
	{ macroN=macro; charge_p=q; mass_p=m; length=len;}
 virtual ~PIC() {}

 int get_size() { return part.size(); }
 double get_macro() { return macroN; }
 double set_macro(double macro) {macroN=macro;}
 double get_m() { return mass_p; }
 double get_q() { return charge_p; }
 double get_gamma(int j) {return sqrt(1.0+(pow(part[j].ux,2)+pow(part[j].uy,2)+pow(part[j].uz,2))/pow(clight,2));}
 double get_vx(int j) {return part[j].ux/get_gamma(j);}
 double get_vy(int j) {return part[j].uy/get_gamma(j);}
 double get_vz(int j) {return part[j].uz/get_gamma(j);}
 double get_energy(int j) { return clight*clight*mass_p*(get_gamma(j)-1.0)/qe;}   //!< kinetic energy in eV
 double get_x(int j) {return part[j].x;}
 double get_y(int j) {return part[j].y;}
 double get_z(int j) {return part[j].z;}

 // Initial distributions

 void set_const_density_xy(double boundary, double T0, long *d1);
 void inject_beam(double boundary, double vel, double angle, double ndensity, double Tempxy, double Tempz, double dt, long *d1); 

 // relativistic particle pusher

 void shift_xyz(double dt, double boundary);
 double shift_xyz(double dt, double boundary, double delta_max, double Emax, double E0, long *d1);
 void shift_uxy(Grid2D& pot, double dt);
 void shift_urz(Grid2D_rz& pot, double dt);
 void Boris_uxy(double By, double dt);
 void Boris_uxy(Grid2D& Bx, Grid2D& By, double dt);
 void Boris_urz(Grid2D_rz& Br, Grid2D_rz& Bz, double dt);

 // ES potential

 void get_potential_xy(Grid2D& pote, Grid2D& rhoe, Greenfb& gf);
 void get_potential_rz(Grid2D_rz& pote, Grid2D_rz& rhoe);

 // Moments

 void densityXY(Grid2D& target);
 void densityRZ(Grid2D_rz& target);
 
 double Tx_rms();
 double vth_xy();
 double x_rms();
 double offset_x();
 double offset_y();

 double total_charge();
 double total_kinetic_energy();
 double total_energy(Grid2D& pot);

 // SEY

 int reflection(double energy, double E0, long *d);
 int secondary(double energy, double delta_max, double Emax, long *d);

}; 

