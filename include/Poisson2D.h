// Green function

class Greenfb {

	vector<Grid2D> greenf; 
	vektor xb, yb;
	double a0, b0;
	
public:

  	Greenfb() {}
	~Greenfb() {}
	Greenfb(Grid2D& g, double a, double b);
	Grid2D& operator[](int j) { return greenf[j]; } 
	double open_2D(double x, double x0, double y, double y0);	
	double circle_2D(double x, double x0, double y, double y0);
	double get_potential_b(int j, Grid2D& rho);
	void add_boundary_potential(Grid2D& rho);
};




// 2D Poisson solver

void poisson_xy(Grid2D& rho, Greenfb& gf);
void poisson_xy(Grid2D& Ex, Grid2D& Ey, Grid2D& rho, Greenfb& gf);
void poisson_rz(Grid2D& rho);

// 3D sliced Poisson solver

void poisson_xyz(Grid3D& Ex, Grid3D& Ey, Grid3D& rho, Greenfb& gf);



