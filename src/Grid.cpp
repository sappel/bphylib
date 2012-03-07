#include "bphylib.h"
#include "mathtools.h"

Grid1D::Grid1D(int n, double diff, double beg) : grid(n), z(n)
{ 
	dz=diff; 
	begin=beg; 
	for(int j=0; j<n; j++) 
		z[j]=begin+(j+0.5)*dz;
}

Grid1D::Grid1D(int n, double diff, double beg, string filename) : grid(n), z(n)
{
  dz=diff; 
  begin=beg; 
  for(int j=0; j<n; j++) 
		z[j]=begin+(j+0.5)*dz;
  outfile.open(filename.c_str(),ios::out);	 
}

double* Grid1D::get_grid()
{
 return (double*)&grid[0];
}


double Grid1D::get_grid_lin(double z0)
{
  int j1, j2;
  double f1, f2;
  double dist0=z0-begin;
  j1=(int)floor(dist0/dz-0.5);
  j2=j1+1;
  f1=((j2+0.5)*dz-dist0)/dz;
  f2=(dist0-(j1+0.5)*dz)/dz;
  return this->operator[](j1)*f1
    + this->operator[](j2)*f2;
}


double Grid1D::get_grid_spline(double z0)
{
  int j1, j2, j3;
  double f1, f2, f3;
  double dist0=z0-begin;
  j2=(int)floor(dist0/dz);
  j3=j2+1;
  j1=j2-1;
  f1=0.5*pow(1.5-fabs(dist0/dz-(j1+0.5)),2);
  f2=3.0/4.0-pow(dist0/dz-(j2+0.5),2);
  f3=0.5*pow(1.5-fabs(dist0/dz-(j3+0.5)),2);
  return this->operator[](j1)*f1
    + this->operator[](j2)*f2
    + this->operator[](j3)*f3;
}


void Grid1D::regrid(double diff, double beg)
{
	int n=grid.size();
	dz=diff; 
	begin=beg; 
	for(int j=0; j<n; j++) 
	      z[j]=begin+(j+0.5)*dz;
}


// returns field with periodic boundary conditions:


double& Grid1D::operator[](int j)
{
	int n = grid.size();
	if(j >= n)
		return grid[j-n];
	else if( j < 0)
		return grid[n+j];
	else return grid[j];
} 

// assignment constructor

Grid1D& Grid1D::operator=(const Grid1D& g)
{
	if ( this != &g )
	{ 
		grid=g.grid;  
		dz=g.dz;
                z=g.z;
		begin=g.begin;
	}
	return *this;
}


//! multiplication with double:


Grid1D& Grid1D::operator*=(double a)
{
  for(int j=0; j < grid.size(); j++)
    {
      grid[j]*=a;
    } 

 return *this; 
}


Grid1D& Grid1D::operator+=(Grid1D& g)
{
  for(int j=0; j < grid.size(); j++)
    {
      grid[j]+=g.grid[j];
    } 
 return *this; 
}


//! Interpolates PIC particle to grid:

void Grid1D::Pic2Field(double pic_quantity, double pic_pos)
{
	int j1, j2, n = grid.size();
	double f1, f2;
	double dist0=pic_pos-begin;
	j1=(int)floor(dist0/dz-0.5);
	j2=j1+1;
	f1=((j2+0.5)*dz-dist0)/dz;
	f2=(dist0-(j1+0.5)*dz)/dz;

	this->operator[](j1) += pic_quantity*f1;
	this->operator[](j2) += pic_quantity*f2;

}


//! Triangular shaped density cloud (TSC), see Hockney


void Grid1D::Pic2FieldTSC(double pic_quantity, double pic_pos)
{
	int j1, j2, j3, n = grid.size();
	double f1, f2, f3;
	double dist0=pic_pos-begin;
	j2=(int)floor(dist0/dz);
	j1=j2-1;
	j3=j2+1;
	double dist2=dist0-(j2+0.5)*dz;
	f1=0.5*pow(0.5-dist2/dz,2);
	f2=3.0/4.0-pow(dist2/dz,2);
	f3=0.5*pow(0.5+dist2/dz,2);

	this->operator[](j1) += pic_quantity*f1;
	this->operator[](j2) += pic_quantity*f2;
	this->operator[](j3) += pic_quantity*f3;

}


// Interpolates field quanity to PIC particle:


double Grid1D::Field2Pic(double pic_pos)
{
  return get_grid_lin(pic_pos);
}


double Grid1D::Field2PicTSC(double pic_pos)
{
  return get_grid_spline(pic_pos);
}




// smooth field:

void Grid1D::smooth()
{
	int j;
	vektor tem(grid);
	double alpha;

	alpha=0.5;
	for(j=1; j<grid.size()-1; j++)
		grid[j]=(alpha*tem[j-1]+tem[j]+alpha*tem[j+1])/(1.0+2.0*alpha);
	j=0;
	grid[j]=(alpha*tem[grid.size()-1]+tem[j]+alpha*tem[j+1])/(1.0+2.0*alpha);
	j=grid.size()-1;
	grid[j]=(alpha*tem[j-1]+tem[j]+alpha*tem[0])/(1.0+2.0*alpha);
}  


void Grid1D::print()
{
 if (outfile.is_open())	
 {
  int n=grid.size(); 
  float *tem=new float[n];
  for(int j=0; j<n; j++)
     tem[j]=grid[j];
  outfile.write((char*)tem,sizeof(float)*n); 
  outfile.flush();
  delete tem;  
 }
}

void Grid1D::reset()
{
 int n=grid.size();
 for(int j=0; j<n; j++)
	 grid[j]=0.0;
}

//------------------min/max--------------------------------


double Grid1D::get_max_abs()
{
	double tem=0.0;
	for(int j=0; j<grid.size(); j++)
		if( tem < fabs(grid[j]) ) tem=grid[j];
	return tem;
}


//-----------------------Grid2D----------------------------------------


Grid2D::Grid2D(int dimx, int dimy, double diffx, double diffy, double begx, double begy) //: x(dimx), y(dimy)
{
	x=vektor(dimx);
	y=vektor(dimy);
	NX=dimx;  
	NY=dimy;
	dx=diffx;
	dy=diffy;
	beginx=begx;
	beginy=begy;
	grid=new double[NX*NY];
	boundary=0.0;   
	for(int j=0; j<dimx; j++) x[j]=beginx+(j+0.5)*dx;
	for(int j=0; j<dimy; j++) y[j]=beginy+(j+0.5)*dy;
}



Grid2D::Grid2D(int dimx, int dimy, double diffx, double diffy) : x(dimx), y(dimy)
{
	NX=dimx;  
	NY=dimy;
	dx=diffx;
	dy=diffy;
	beginx=-0.5*dimx*dx;
	beginy=-0.5*dimy*dy;
	grid=new double[NX*NY];
	boundary=0.0;   
	for(int j=0; j<dimx; j++) x[j]=beginx+(j+0.5)*dx;
	for(int j=0; j<dimy; j++) y[j]=beginy+(j+0.5)*dy;
}


Grid2D::Grid2D(int dimx, int dimy, double diffx, double diffy, string filename) : x(dimx), y(dimy) 
{
    NX=dimx;  
    NY=dimy;
	dx=diffx;
	dy=diffy;
	beginx=-0.5*dimx*dx;
	beginy=-0.5*dimy*dy;
	grid=new double[NX*NY];
	boundary=0.0;   
    for(int j=0; j<dimx; j++) x[j]=beginx+(j+0.5)*dx;
	for(int j=0; j<dimy; j++) y[j]=beginy+(j+0.5)*dy;
	outfile.open(filename.c_str(),ios::out); 
}  



Grid2D::Grid2D(const Grid2D& g)
{
if ( this != &g )
	{ 
		grid=g.grid;  
		dx=g.dx;
		dy=g.dy;
		NX=g.NX;
		NY=g.NY;
		beginx=g.beginx;
		beginy=g.beginy;
		x=g.x;
		y=g.y;
                //if(g.outfile.is_open()) g.outfile.close();
		//if(outfile.is_open()) outfile.close();
	}
}

Grid2D& Grid2D::operator=(const Grid2D& g)
{
	if ( this != &g )
	{ 
		if ( NX==0 && NY==0 )
		  grid=new double[g.NX*g.NY];
		for(int l=0; l<g.NY; l++)
		   for(int j=0; j<g.NX; j++)
			grid[j*g.NY+l]=g.grid[j*g.NY+l];
		dx=g.dx;
		dy=g.dy;
		NX=g.NX;
		NY=g.NY;
		beginx=g.beginx;
		beginy=g.beginy;
		x=g.x;
		y=g.y;
	        //if(g.outfile.is_open()) g.outfile.close();
		//if(outfile.is_open()) outfile.close();
	}
	return *this;
}


// periodic bc in x and y : 

/*
double& Grid2D::operator()(int j, int l) 
{ 
	//if(l >= NY || l < 0 )
	//{ boundary=0.0;
	//  return boundary;}
	if (l < 0)
        return grid[l+NY+j*NY];
	else if (l >= NY)
	    return grid[l-NY+j*NY];
	else if(j < 0)  
		return grid[l+(NX+j)*NY];
	else if(j >= NX)
		return grid[l+(j-NX)*NY];
	else
		return grid[l+j*NY];
}
*/

double& Grid2D::operator()(int j, int l)
{
       if (l < 0)
               return grid[l+NY+j*NY];
       else if (l >= NY)
             {
               if(j>0) return grid[l-NY+j*NY];
               else return grid[l-NY+(NX+j)*NY];
               }
       else if(j < 0)
               return grid[l+(NX+j)*NY];
       else if(j >= NX)
               return grid[l+(j-NX)*NY];
       else
               return grid[l+j*NY];
}



double& Grid2D::operator[](int j)
{
	return grid[j];
} 


void Grid2D::reset()
{
	int j,l;
	for(l=0; l<NY; l++)
		for(j=0; j<NX; j++)
			grid[j*NY+l]=0.0;	
}


void Grid2D::filtering()
{
	int j,l;
	double a=0.5;
	vektor tem(NX); 

	for(l=0; l<NY; l++)
	{
		for(j=0;j<NX;j++)
			tem[j]=grid[j*NY+l];
		for(j=1;j<NX-1;j++)
			grid[j*NY+l]=(a*tem[j-1]+tem[j]+a*tem[j+1])/(1.0+2.0*a);
		j=0;
		grid[j*NY+l]=(a*tem[NX-1]+tem[j]+a*tem[j+1])/(1.0+2.0*a);
		j=NX-1;
		grid[j*NY+l]=(a*tem[j-1]+tem[j]+a*tem[0])/(1.0+2.0*a);     
	} 

}

double Grid2D::Grid2dx(double pic_pos_x, double pic_pos_y)
{
  double g1, g2;
  g2=Grid2PIC(pic_pos_x+dx,pic_pos_y);
  g1=Grid2PIC(pic_pos_x-dx,pic_pos_y);
  return (g2-g1)/(2.0*dx);	
}


double Grid2D::Grid2dy(double pic_pos_x, double pic_pos_y)
{
  double g1, g2;
  g2=Grid2PIC(pic_pos_x,pic_pos_y+dy);
  g1=Grid2PIC(pic_pos_x,pic_pos_y-dy);
  return (g2-g1)/(2.0*dy);	
}


void Grid2D::Pic2Grid(double pic_charge, double pic_pos_x, double pic_pos_y)
{
	int j1,j2,l1,l2;
	double a1, a2, a3, a4;

	j1=(int)floor((pic_pos_x+0.5*(NX-1.0)*dx)/dx);
	j2=j1+1;
	l1=(int)floor((pic_pos_y+0.5*(NY-1.0)*dy)/dy); 
	l2=l1+1;

	a1=(-0.5*(NX-1.0)*dx+j2*dx-pic_pos_x)*(pic_pos_y+0.5*(NY-1.0)*dy-l1*dy);
	a2=(-0.5*(NX-1.0)*dx+j2*dx-pic_pos_x)*(-0.5*(NY-1.0)*dy+l2*dy-pic_pos_y);
	a3=(pic_pos_x+0.5*(NX-1.0)*dx-j1*dx)*(pic_pos_y+0.5*(NY-1.0)*dy-l1*dy);
	a4=(pic_pos_x+0.5*(NX-1.0)*dx-j1*dx)*(-0.5*(NY-1.0)*dy+l2*dy-pic_pos_y);
	
	this->operator()(j1,l1) += a2*pic_charge/pow(dx*dy,2);
	this->operator()(j1,l2) += a1*pic_charge/pow(dx*dy,2);
	this->operator()(j2,l1) += a4*pic_charge/pow(dx*dy,2);
	this->operator()(j2,l2) += a3*pic_charge/pow(dx*dy,2);

}


double Grid2D::Grid2PIC(double pic_pos_x, double pic_pos_y)
{
	int j1,j2,l1,l2;
	double a1, a2, a3, a4;

	j1=(int)floor((pic_pos_x+0.5*(NX-1.0)*dx)/dx);
	j2=j1+1;
	l1=(int)floor((pic_pos_y+0.5*(NY-1.0)*dy)/dy); 
	l2=l1+1;
	
	a1=(-0.5*(NX-1.0)*dx+j2*dx-pic_pos_x)*(pic_pos_y+0.5*(NY-1.0)*dy-l1*dy);
	a2=(-0.5*(NX-1.0)*dx+j2*dx-pic_pos_x)*(-0.5*(NY-1.0)*dy+l2*dy-pic_pos_y);
	a3=(pic_pos_x+0.5*(NX-1.0)*dx-j1*dx)*(pic_pos_y+0.5*(NY-1.0)*dy-l1*dy);
	a4=(pic_pos_x+0.5*(NX-1.0)*dx-j1*dx)*(-0.5*(NY-1.0)*dy+l2*dy-pic_pos_y);

	return this->operator()(j1,l1)*a2/(dx*dy)+this->operator()(j1,l2)*a1/(dx*dy)+
	this->operator()(j2,l1)*a4/(dx*dy)+this->operator()(j2,l2)*a3/(dx*dy);

}

Grid2D& Grid2D::operator+=(const Grid2D& g)
{
 	for(int l=0; l<g.NY; l++)
	   for(int j=0; j<g.NX; j++)
		grid[j*g.NY+l]+=g.grid[j*g.NY+l];	
	return *this;
}

Grid2D& Grid2D::operator-=(const Grid2D& g)
{
 	for(int l=0; l<g.NY; l++)
	   for(int j=0; j<g.NX; j++)
		grid[j*g.NY+l]-=g.grid[j*g.NY+l];	
	return *this;
}

Grid2D& Grid2D::operator-=(double a)
{
 	for(int l=0; l<NY; l++)
	   for(int j=0; j<NX; j++)
		grid[j*NY+l]-=a;	
	return *this;
}

Grid2D& Grid2D::operator*=(double a)
{
 	for(int l=0; l<NY; l++)
	   for(int j=0; j<NX; j++)
		grid[j*NY+l]*=a;	
	return *this;
}


void Grid2D::regrid(double diffy, double begy)
{
	vektor y_new(NY), ys(NY), tem(NY);
	for(int l=0; l<NY; l++) y_new[l]=begy+(0.5+l)*diffy;
	for(int j=0; j<NX; j++)
	{
		for(int l=0; l<NY; l++) tem[l]=grid[l+j*NY];
		spline(&ys[0],&tem[0],NY);
		for(int l=0; l<NY; l++)
			grid[j*NY+l]=splint(&tem[0],&ys[0],&y[0],y_new[l],NY);

	}
	for(int l=0; l<NY; l++) y[l]=y_new[l];
	beginy=begy;
	dy=diffy;
}	


void Grid2D::print()
{
 if (outfile.is_open())	
 {
  float *tem=new float[NY];
  for(int j=0; j<NX; j++)
   {
    for(int l=0; l<NY; l++)
	tem[l]=grid[l+NY*j];
    outfile.write((char*)tem,sizeof(float)*NY); 
   }
  outfile.flush();
  delete tem;
 }
}

//----------------------------r-z------------------------------

double& Grid2D_rz::operator()(int j, int l) // corners not correct !!!!!
{ 
	// four corners
	if ( l < 0 && j < 0)
      return this->operator[](l+NY+abs(j)*NY);  	
	if ( l < 0 && j >= NX)
      return boundary;
	if (l >= NY && j < 0)
	    return this->operator[](l-NY+abs(j)*NY);
	if (l >= NY && j >= NX)
	    return boundary;
	
	// boundaries
	if (l < 0)
        return this->operator[](l+NY+j*NY);
	else if (l >= NY)
	    return this->operator[](l-NY+j*NY);
	else if(j < 0)  
		return this->operator[](l+abs(j)*NY);
	else if(j >= NX)
		return boundary;
	else
		return this->operator[](l+j*NY);   // default
}


void Grid2D_rz::Pic2Grid(double pic_charge, double pic_pos_x, double pic_pos_y)
{
	int j1,j2,l1,l2;
	double a1, a2, a3, a4;

	j1=(int)floor(pic_pos_x/dx);
	j2=j1+1;
	l1=(int)floor((pic_pos_y+0.5*(NY-1.0)*dy)/dy); 
	l2=l1+1;
  
	a1=(j2*dx-pic_pos_x)*(pic_pos_y+0.5*(NY-1.0)*dy-l1*dy);
	a2=(j2*dx-pic_pos_x)*(-0.5*(NY-1.0)*dy+l2*dy-pic_pos_y);
	a3=(pic_pos_x-j1*dx)*(pic_pos_y+0.5*(NY-1.0)*dy-l1*dy);
	a4=(pic_pos_x-j1*dx)*(-0.5*(NY-1.0)*dy+l2*dy-pic_pos_y);
	
	if(j1 == 0)
	 {
      		this->operator()(j1,l1) += a2*pic_charge/(PI/3.0*pow(dx,3)*pow(dy,2));	
			this->operator()(j1,l2) += a1*pic_charge/(PI/3.0*pow(dx,3)*pow(dy,2));	
			this->operator()(j2,l1) += a4*pic_charge/(2.0*PI*j2*pow(dx,3)*pow(dy,2));
			this->operator()(j2,l2) += a3*pic_charge/(2.0*PI*j2*pow(dx,3)*pow(dy,2));
 	 }
    else { 
	   this->operator()(j1,l1) += a2*pic_charge/(2.0*PI*j1*pow(dx,3)*pow(dy,2));
	   this->operator()(j1,l2) += a1*pic_charge/(2.0*PI*j1*pow(dx,3)*pow(dy,2));
	   this->operator()(j2,l1) += a4*pic_charge/(2.0*PI*j2*pow(dx,3)*pow(dy,2));
	   this->operator()(j2,l2) += a3*pic_charge/(2.0*PI*j2*pow(dx,3)*pow(dy,2));
     }	
}


double Grid2D_rz::Grid2PIC(double pic_pos_x, double pic_pos_y)
{
	int j1,j2,l1,l2;
	double a1, a2, a3, a4;

	j1=(int)floor(pic_pos_x/dx);
	j2=j1+1;
	l1=(int)floor((pic_pos_y+0.5*(NY-1.0)*dy)/dy); 
	l2=l1+1;
		
	a1=(j2*dx-pic_pos_x)*(pic_pos_y+0.5*(NY-1.0)*dy-l1*dy);
	a2=(j2*dx-pic_pos_x)*(-0.5*(NY-1.0)*dy+l2*dy-pic_pos_y);
	a3=(pic_pos_x-j1*dx)*(pic_pos_y+0.5*(NY-1.0)*dy-l1*dy);
	a4=(pic_pos_x-j1*dx)*(-0.5*(NY-1.0)*dy+l2*dy-pic_pos_y);

	return this->operator()(j1,l1)*a2/(dx*dy)+this->operator()(j1,l2)*a1/(dx*dy)+
	this->operator()(j2,l1)*a4/(dx*dy)+this->operator()(j2,l2)*a3/(dx*dy);

}


double Grid2D_rz::Grid2dr(double pic_pos_x, double pic_pos_y)
{
  double g1, g2;
  g2=this->Grid2PIC(pic_pos_x+dx,pic_pos_y);
  g1=this->Grid2PIC(pic_pos_x-dx,pic_pos_y);
  return (g2-g1)/(2.0*dx);	
}


double Grid2D_rz::Grid2dz(double pic_pos_x, double pic_pos_y)
{
  double g1, g2;
  g2=this->Grid2PIC(pic_pos_x,pic_pos_y+dy);
  g1=this->Grid2PIC(pic_pos_x,pic_pos_y-dy);
  return (g2-g1)/(2.0*dy);	
}



//---------------------------------3D----------------------------------


Grid3D::Grid3D(int NZ, double zl, double zr, Grid2D& g) : grid(NZ)
{
 zleft=zl; zright=zr;	
 for(int j=0; j<NZ; j++)
 {
   grid[j]=g;
   grid[j].reset();
 }
 ghostl=g;
 ghostl.reset();
 ghostr=g;
 ghostr.reset();	 	 
}


void Grid3D::reset()
{
 for(int j=0; j<grid.size(); j++)
     grid[j].reset();
 ghostl.reset();
 ghostr.reset();	
}

/*
void Grid3D::Pic2Grid(double pic_quantity, double pic_pos_x, double pic_pos_y, double pic_pos_z)
{
 double dist0=pic_pos_z-zleft;
 double dz=get_dz(); 
 int j0=(int)floor(dist0/dz);
 if (j0 < 0) return; // was j0=0;
 if (j0 >= grid.size() ) return; // j0=grid.size()-1; 
 grid[j0].Pic2Grid(pic_quantity,pic_pos_x,pic_pos_y);	
}
*/

void Grid3D::Pic2Grid(double pic_quantity, double pic_pos_x, double pic_pos_y, double pic_pos_z)
{
 double dist0=pic_pos_z-zleft;
 double dz=get_dz(); 
 int j1=(int)floor(dist0/dz-0.5);
 int j2=j1+1;
 double f1=((j2+0.5)*dz-dist0)/dz;
 double f2=(dist0-(j1+0.5)*dz)/dz; 

 if (j1 < 0) return; // was j0=0;
 if (j2 >= grid.size() ) return; // j0=grid.size()-1; 
 
 if (j1==1)  {
    ghostl.Pic2Grid(f1*pic_quantity,pic_pos_x,pic_pos_y);	
	grid[j1].Pic2Grid(f2*pic_quantity,pic_pos_x,pic_pos_y);	
   }
 else if (j2 == grid.size()-1)  {
    ghostr.Pic2Grid(f2*pic_quantity,pic_pos_x,pic_pos_y);	
	grid[j2].Pic2Grid(f1*pic_quantity,pic_pos_x,pic_pos_y);	
   }
 else {  
  grid[j1].Pic2Grid(f1*pic_quantity,pic_pos_x,pic_pos_y);	
  grid[j2].Pic2Grid(f2*pic_quantity,pic_pos_x,pic_pos_y);
 }	
}

/*
double Grid3D::Grid2PIC(double pic_pos_x, double pic_pos_y, double pic_pos_z)
{
 double dist0=pic_pos_z-zleft;
 double dz=get_dz(); 
 int j0=(int)floor(dist0/dz);
 if (j0 < 0) return 0.0; // j0=0;
 if (j0 >= grid.size() ) return 0.0; //j0=grid.size()-1; 
 return grid[j0].Grid2PIC(pic_pos_x,pic_pos_y);	
}
*/


double Grid3D::Grid2PIC(double pic_pos_x, double pic_pos_y, double pic_pos_z)
{
 double dist0=pic_pos_z-zleft;
 double dz=get_dz(); 
 int j1=(int)floor(dist0/dz-0.5);
 int j2=j1+1;
 double f1=((j2+0.5)*dz-dist0)/dz;
 double f2=(dist0-(j1+0.5)*dz)/dz; 

 if (j1 < 0) return 0.0; // was j0=0;
 if (j2 >= grid.size() ) return 0.0; // j0=grid.size()-1; 
 if (j1==1)  return f1*ghostl.Grid2PIC(pic_pos_x,pic_pos_y)+f2*grid[j1].Grid2PIC(pic_pos_x,pic_pos_y);
 if (j2 == grid.size()-1) return f2*ghostr.Grid2PIC(pic_pos_x,pic_pos_y)+f1*grid[j2].Grid2PIC(pic_pos_x,pic_pos_y);
 return f1*grid[j1].Grid2PIC(pic_pos_x,pic_pos_y)+f2*grid[j2].Grid2PIC(pic_pos_x,pic_pos_y);	 
}


//-----------------Interpolation schemes--------------------------


double Interpolate::twopoint(vektor& v, int j, double shift)
{
	int bc=0; // periodic boundary conditions
	int n = v.size();

	if(fabs(shift) > 1.0) 
	{ 
		printf("shift > 1 !! j=%d  bc=%d shift=%g\n",j,bc,shift); 
		exit(0);
	}

	if(shift < 0.0 && j < n-1 )
		return v[j]-shift*(v[j+1]-v[j]);
	if(shift >= 0.0 && j > 0 )
		return v[j]-shift*(v[j]-v[j-1]);

	if( bc == 0)  // periodic
	{
		if( j == 0 ) 
			return v[j]-shift*(v[j]-v[n-1]);   
		if( j == n-1 )
			return v[j]-shift*(v[0]-v[j]);
	}

	if( bc == 1)  // v=0 at boundary
	{
		if( j == 0 ) 
			return v[j]-shift*v[j];   
		if( j == n-1 )
			return v[j]-shift*(-v[j]);
	}

	printf("something wrong in interpolate \n");
	return 0.0;
}


double Interpolate::threepoint(vektor& v, int j, double shift)
{
	int bc=0; // periodic boundary conditions
	int n = v.size();

	if(fabs(shift) > 1.0) 
	{ 
		printf("shift > 1 !! j=%d  bc=%d shift=%g\n",j,bc,shift); 
		exit(0);
	}

	if(shift < 0.0 && j < n-2 )   
		return v[j]-shift*(2.0*v[j+1]-1.5*v[j]-0.5*v[j+2])
			+0.5*pow(shift,2)*(v[j+2]-2.0*v[j+1]+v[j]);
	if(shift >= 0.0 && j > 1 )   
		return v[j]-shift*(1.5*v[j]-2.0*v[j-1]+0.5*v[j-2])
			+0.5*pow(shift,2)*(v[j]-2.0*v[j-1]+v[j-2]); 

	if(bc == 0) // periodic 
	{
		if( j == n-2 )
			return v[j]-shift*(2.0*v[j+1]-1.5*v[j]-0.5*v[0])
				+0.5*pow(shift,2)*(v[0]-2.0*v[j+1]+v[j]);
		if( j == n-1 )
			return v[j]-shift*(2.0*v[0]-1.5*v[j]-0.5*v[1])
				+0.5*pow(shift,2)*(v[1]-2.0*v[0]+v[j]); 
		if( j == 1 )
			return v[j]-shift*(1.5*v[j]-2.0*v[j-1]+0.5*v[n-1])
				+0.5*pow(shift,2)*(v[j]-2.0*v[j-1]+v[n-1]);  
		if( j == 0 )
			return v[j]-shift*(1.5*v[j]-2.0*v[n-1]+0.5*v[n-2])
				+0.5*pow(shift,2)*(v[j]-2.0*v[n-1]+v[n-2]);   
	}

	if(bc == 1) // v=0 at boundary
	{
		if( j == n-2 )
			return v[j]-shift*(2.0*v[j+1]-1.5*v[j])
				+0.5*pow(shift,2)*(-2.0*v[j+1]+v[j]);
		if( j == n-1 )
			return v[j]-shift*(-1.5*v[j])
				+0.5*pow(shift,2)*(-2.0*v[0]); 
		if( j == 1 )
			return v[j]-shift*(1.5*v[j]-2.0*v[j-1])
				+0.5*pow(shift,2)*(v[j]-2.0*v[j-1]);  
		if( j == 0 )
			return v[j]-shift*(1.5*v[j])
				+0.5*pow(shift,2)*(v[j]);   
	}

	printf("something wrong in interpolate \n");
	return 0.0;
}


double Interpolate::fourpoint(vektor& v, int j0, double shift0)
{
	int bc=0; // periodic boundary conditions
	int n = v.size(), j=j0;
	double shift=shift0;

	// reduce shift < 1, if > 1:
	if(fabs(shift0) >= 1.0) { 
		if(shift0 > 0.0)
		{
			shift=shift0-floor(fabs(shift0));
			j=j0-(int)floor(fabs(shift0));
			if(j < 0) return 0.0;
		}
		else if(shift0 < 0.0)
		{
			shift=shift0+floor(fabs(shift0));
			j=j0+(int)floor(fabs(shift0));
			if(j >= n) return 0.0;
		}
	}	

	if(shift < 0.0 && j < n-2 && j > 0)   
		return v[j]-shift*(v[j+1]-0.5*v[j]-1.0/3.0*v[j-1]-1.0/6.0*v[j+2])
			+0.5*pow(shift,2)*(v[j+1]-2.0*v[j]+v[j-1])
		-pow(shift,3)/6.0*(v[j+2]-3.0*v[j+1]+3.0*v[j]-v[j-1]);    
	if(shift >= 0.0 && j > 1 && j < n-1 ) 
		return v[j]-shift*(1.0/3.0*v[j+1]+0.5*v[j]-v[j-1]+1.0/6.0*v[j-2])
			+0.5*pow(shift,2)*(v[j+1]-2.0*v[j]+v[j-1])
		-pow(shift,3)/6.0*(v[j+1]-3.0*v[j]+3.0*v[j-1]-v[j-2]);


	if(bc == 0) // periodic 
	{
		if(shift < 0.0 && j == n-2 )
			return v[j]-shift*(v[j+1]-0.5*v[j]-1.0/3.0*v[j-1]-1.0/6.0*v[0])
				+0.5*pow(shift,2)*(v[j+1]-2.0*v[j]+v[j-1])
			-pow(shift,3)/6.0*(v[0]-3.0*v[j+1]+3.0*v[j]-v[j-1]);  
		if(shift < 0.0 && j == n-1 )
			return v[j]-shift*(v[0]-0.5*v[j]-1.0/3.0*v[j-1]-1.0/6.0*v[1])
				+0.5*pow(shift,2)*(v[0]-2.0*v[j]+v[j-1])
			-pow(shift,3)/6.0*(v[1]-3.0*v[0]+3.0*v[j]-v[j-1]);
		if(shift < 0.0 && j == 0 )
			return v[j]-shift*(v[j+1]-0.5*v[j]-1.0/3.0*v[n-1]-1.0/6.0*v[j+2])
				+0.5*pow(shift,2)*(v[j+1]-2.0*v[j]+v[n-1])
			-pow(shift,3)/6.0*(v[j+2]-3.0*v[j+1]+3.0*v[j]-v[n-1]);    
		if(shift >= 0.0 && j == 0 )        
			return v[j]-shift*(1.0/3.0*v[j+1]+0.5*v[j]-v[n-1]+1.0/6.0*v[n-2])
				+0.5*pow(shift,2)*(v[j+1]-2.0*v[j]+v[n-1])
			-pow(shift,3)/6.0*(v[j+1]-3.0*v[j]+3.0*v[n-1]-v[n-2]); 
		if(shift >= 0.0 && j == 1 )
			return v[j]-shift*(1.0/3.0*v[j+1]+0.5*v[j]-v[j-1]+1.0/6.0*v[n-1])
				+0.5*pow(shift,2)*(v[j+1]-2.0*v[j]+v[j-1])
			-pow(shift,3)/6.0*(v[j+1]-3.0*v[j]+3.0*v[j-1]-v[n-1]); 
		if(shift >= 0.0 && j == n-1 ) 
			return v[j]-shift*(1.0/3.0*v[0]+0.5*v[j]-v[j-1]+1.0/6.0*v[j-2])
				+0.5*pow(shift,2)*(v[0]-2.0*v[j]+v[j-1])
			-pow(shift,3)/6.0*(v[0]-3.0*v[j]+3.0*v[j-1]-v[j-2]); 
	}

	if(bc == 1) // v=0 at boundary
	{
		if(shift < 0.0 && j == n-2 )
			return v[j]-shift*(v[j+1]-0.5*v[j]-1.0/3.0*v[j-1])
				+0.5*pow(shift,2)*(v[j+1]-2.0*v[j]+v[j-1])
			-pow(shift,3)/6.0*(-3.0*v[j+1]+3.0*v[j]-v[j-1]);  
		if(shift < 0.0 && j == n-1 )
			return v[j]-shift*(-0.5*v[j]-1.0/3.0*v[j-1])
				+0.5*pow(shift,2)*(-2.0*v[j]+v[j-1])
			-pow(shift,3)/6.0*(3.0*v[j]-v[j-1]);
		if(shift < 0.0 && j == 0 )
			return v[j]-shift*(v[j+1]-0.5*v[j]-1.0/6.0*v[j+2])
				+0.5*pow(shift,2)*(v[j+1]-2.0*v[j])
			-pow(shift,3)/6.0*(v[j+2]-3.0*v[j+1]+3.0*v[j]);    
		if(shift >= 0.0 && j == 0 )        
			return v[j]-shift*(1.0/3.0*v[j+1]+0.5*v[j])
				+0.5*pow(shift,2)*(v[j+1]-2.0*v[j])
			-pow(shift,3)/6.0*(v[j+1]-3.0*v[j]); 
		if(shift >= 0.0 && j == 1 )
			return v[j]-shift*(1.0/3.0*v[j+1]+0.5*v[j]-v[j-1])
				+0.5*pow(shift,2)*(v[j+1]-2.0*v[j]+v[j-1])
			-pow(shift,3)/6.0*(v[j+1]-3.0*v[j]+3.0*v[j-1]); 
		if(shift >= 0.0 && j == n-1 ) 
			return v[j]-shift*(0.5*v[j]-v[j-1]+1.0/6.0*v[j-2])
				+0.5*pow(shift,2)*(-2.0*v[j]+v[j-1])
			-pow(shift,3)/6.0*(-3.0*v[j]+3.0*v[j-1]-v[j-2]); 
	}  

	printf("something wrong in interpolate \n");
	return 0.0;

}


double Interpolate::flux_balance(vektor& v, int j, double shift)
{
	double tem1, tem2;

	int n = v.size();

	if(fabs(shift) >= 1.0) 
	{ 
		printf("flux_balance: shift > 1 !! j=%d  shift=%g\n",j,shift); 
		exit(0);
	}

	if( shift >= 0.0 )
	{
		if( j>1 && j<n-1)
		{
			tem2=shift*(v[j]+(v[j+1]-v[j-1])*(1.0-shift)/4.0);
			tem1=shift*(v[j-1]+(v[j]-v[j-2])*(1.0-shift)/4.0);
		}
		else if( j==0 )
		{
			tem2=shift*(v[j]+(v[j+1]-v[n-1])*(1.0-shift)/4.0);
			tem1=shift*(v[n-1]+(v[j]-v[n-2])*(1.0-shift)/4.0);
		}
		else if( j==1 )
		{
			tem2=shift*(v[j]+(v[j+1]-v[j-1])*(1.0-shift)/4.0);
			tem1=shift*(v[j-1]+(v[j]-v[n-1])*(1.0-shift)/4.0);
		}
		else if( j==n-1 )
		{
			tem2=shift*(v[j]+(v[0]-v[j-1])*(1.0-shift)/4.0);
			tem1=shift*(v[j-1]+(v[j]-v[j-2])*(1.0-shift)/4.0);
		}
		else
		{  
			printf("flux_balance: wrong j !!"); 
			exit(0);
		}
	}
	else if( shift < 0.0 )
	{
		if( j>0 && j<n-2 )
		{
			tem2=fabs(shift)*(v[j]-(v[j+1]-v[j-1])*(1.0-fabs(shift))/4.0);
			tem1=fabs(shift)*(v[j+1]-(v[j+2]-v[j])*(1.0-fabs(shift))/4.0);          
		}
		else if( j==0 )
		{
			tem2=fabs(shift)*(v[j]-(v[j+1]-v[n-1])*(1.0-fabs(shift))/4.0);
			tem1=fabs(shift)*(v[j+1]-(v[j+2]-v[j])*(1.0-fabs(shift))/4.0);          
		} 
		else if( j==n-1 ) 
		{
			tem2=fabs(shift)*(v[j]-(v[0]-v[j-1])*(1.0-fabs(shift))/4.0);
			tem1=fabs(shift)*(v[0]-(v[1]-v[j])*(1.0-fabs(shift))/4.0);          
		}
		else if( j==n-2 )
		{
			tem2=fabs(shift)*(v[j]-(v[j+1]-v[j-1])*(1.0-fabs(shift))/4.0);
			tem1=fabs(shift)*(v[j+1]-(v[0]-v[j])*(1.0-fabs(shift))/4.0);          
		}
		else
		{  
			printf("flux_balance: wrong j !!"); 
			exit(0);
		}

	}

	return v[j]+tem1-tem2;

}
