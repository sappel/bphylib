/*! 
One dimensional grid with periodic boundary conditions
and interpolation (linear, spline) routines for PIC codes
*/

class Grid1D {

	vektor grid;
	double dz, begin; 
	ofstream outfile;

public:

	vektor z;

	Grid1D() {}
	Grid1D(int n, double diff, double beg); 
	Grid1D(int n, double diff, double beg, string filename); 
	~Grid1D() { if(outfile.is_open()) outfile.close(); }

	Grid1D& operator=(const Grid1D& g);
	double& operator[](int j); 
        Grid1D& operator*=(double a);
        Grid1D& operator+=(Grid1D& g);
        void set_grid(int j, double a) { grid[j]=a; }
        double get_grid(int j) { return grid[j]; }
	double* get_grid();
        double get_grid_lin(double z0);
        double get_grid_spline(double z0); 
	void regrid(double diff, double beg);
	int get_size() { return grid.size(); } 
	double get_dz() { return dz; } 
	double get_begin() { return begin; }
	double get_max_abs();
	double Field2Pic(double pic_pos); 
	double Field2PicTSC(double pic_pos); 
	void Pic2Field(double pic_quantity ,double pic_pos);
	void Pic2FieldTSC(double pic_quantity ,double pic_pos);
	void smooth();
	void print();
	void reset();

};


/*! 
Two dimensional grid with periodic boundary conditions
and interpolation routines (linear) for PIC codes
*/

class Grid2D {

private:	
	
	double *grid;
    ofstream outfile;

protected:
	
	int NX, NY;
	double dx, dy;
    double beginx, beginy;
	double boundary;
	
public:

	vektor x, y; 

	Grid2D() { grid=0; NX=0; NY=0;}  
	Grid2D(int dimx, int dimy, double diffx, double diffy);  
	Grid2D(int dimx, int dimy, double diffx, double diffy, double begx, double begy);  
	Grid2D(int dimx, int dimy, double diffx, double diffy, string filename); 
	~Grid2D() { delete [] grid; if(outfile.is_open()) outfile.close();}
	Grid2D(const Grid2D& g);
	Grid2D& operator=(const Grid2D& g);
	double* get_grid() {return grid;}
	int get_size() { return NX*NY; }
	int get_NX() {return NX;}
	int get_NY() {return NY;}
	double get_dx() {return dx;}
	double get_dy() {return dy;}
	double get_beginx() { return beginx; }
	double get_beginy() { return beginy; }
	double& operator()(int j, int l);     
	double& operator[](int j);
	Grid2D& operator+=(const Grid2D& g);
	Grid2D& operator-=(const Grid2D& g);
    Grid2D& operator-=(double a);
	Grid2D& operator*=(double a);
	void reset();
	void filtering();
	double Grid2dx(double pic_pos_x, double pic_pos_y);
	double Grid2dy(double pic_pos_x, double pic_pos_y);
	void Pic2Grid(double pic_quantity, double pic_pos_x, double pic_pos_y); 
	double Grid2PIC(double pic_pos_x, double pic_pos_y);
	void regrid(double diffy, double begy); 
	void print();
};


class Grid2D_rz : public Grid2D {

 public:
 
	Grid2D_rz(int dimr, int dimz, double diffr, double diffz) : Grid2D(dimr,dimz,diffr,diffz) 
		{ beginx=0.0; for(int j=0; j<dimr; j++) x[j]=j*dx; };
	double& operator()(int j, int l); 	    
	double Grid2dr(double pic_pos_r, double pic_pos_z);
	double Grid2dz(double pic_pos_r, double pic_pos_z);
	void Pic2Grid(double pic_quantity, double pic_pos_r, double pic_pos_z); 
	double Grid2PIC(double pic_pos_r, double pic_pos_z);
};


class Grid3D {

    vector<Grid2D> grid;
	// longitudinal boundaries of the bunch
	double zleft, zright;
	Grid2D ghostl, ghostr;

public:

	Grid3D() {}
	Grid3D(int NZ, double zl, double zr, Grid2D& g);
	~Grid3D() {}
	Grid2D& operator[](int j) { return grid[j]; } 
	void reset();
	// left and right bunch ends can be updated
	double& get_zleft() { return zleft; }
	double& get_zright() { return zright; }
	// This dz can be different from the dz in Grid1D. 
	double get_dz() { return (zright-zleft)/grid.size(); }
	// dx,dy taken from Grid2D
	double get_dx() { return grid[0].get_dx(); }
	double get_dy() { return grid[0].get_dy(); }
	int get_NX() { return grid[0].get_NX(); }
	int get_NY() { return grid[0].get_NY(); }
	int get_NZ() { return grid.size(); } 
	void Pic2Grid(double pic_quantity, double pic_pos_x, double pic_pos_y, double pic_pos_z); 
	double Grid2PIC(double pic_pos_x, double pic_pos_y, double pic_pos_z);
	// ghost cells
	double* get_ghostl() {return ghostl.get_grid();}
    double* get_ghostr() {return ghostr.get_grid();}
	
	
};



/*!
Interpolation routines for advancing the distribution function (Vlasov)
*/

namespace Interpolate
{
	double twopoint(vektor& v, int j, double shift);
	double threepoint(vektor& v, int j, double shift);
	double fourpoint(vektor& v, int j, double shift);
	double flux_balance(vektor& v, int j, double shift);
}	
