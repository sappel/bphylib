#include "PhyConstants.h"
#include "Definitions.h"
#include "Grid.h"
#include "mathtools.h"
#include "Poisson2D.h"
//#include "fftw3.h"

Greenfb::Greenfb(Grid2D& g, double a, double b) : greenf(2*g.get_NX()+2*g.get_NY()-8)
{
  int i, j, v, u=0;
  int Nb=greenf.size(); 	
  int NX=g.get_NX();
  int NY=g.get_NY();
  xb=vektor(Nb);
  yb=vektor(Nb);
  a0=a;
  b0=b;

  for(j=0; j<Nb; j++)
	 {
	   greenf[j]=g;
	   greenf[j].reset();
	 }

  for(i=1; i<NX; i++)
      {
		xb[u]=g.x[i];
		yb[u]=g.y[NY-1];
		u=u+1;
	  }	
  for(j=1; j<NY-1; j++)	
      {
		xb[u]=g.x[NX-1];
		yb[u]=g.y[j];
		u=u+1;
      }   
  for(i=1; i<NX-1; i++)
      {
		xb[u]=g.x[i];
		yb[u]=g.y[1];
	    u=u+1;
	  }
  for(j=2; j<NY-1; j++)	  
	 {
		xb[u]=g.x[1];
		yb[u]=g.y[j];
		u=u+1;
	  }
     	
  for(int j=0; j<Nb; j++)
    for (u=0; u<NX; u++)
		for(v=0; v<NY; v++)
		{
		 if (a0 == 0.0) greenf[j](u,v)=open_2D( g.x[u], xb[j], g.y[v], yb[j]);
		 else if (a0 > 0.0 && b0 == 0.0) greenf[j](u,v)=circle_2D( g.x[u], xb[j], g.y[v], yb[j]);
		 else greenf[j](u,v)=0.0;
		}

}


double Greenfb::open_2D(double x, double x0, double y, double y0)
{
	double R=sqrt(pow(x-x0,2)+pow(y-y0,2));
	if ( R == 0.0 ) return 0.0;
	return -1.0/(2.0*PI*eps0)*log(R);
}


double Greenfb::circle_2D(double x, double x0, double y, double y0)
{
 double R=sqrt(pow(x-x0,2)+pow(y-y0,2));
 if ( R == 0.0 ) return 0.0;
 double r0=sqrt(pow(x0,2)+pow(y0,2));	
 double Rd=sqrt(pow(pow(r0,2)*x-pow(a0,2)*x0,2)+pow(pow(r0,2)*y-pow(a0,2)*y0,2))/pow(r0,2);
 return 1.0/(2.0*PI*eps0)*log(Rd*r0/(R*a0));	
}


double Greenfb::get_potential_b(int j, Grid2D& rho)
{
	double tem=0.0;
	int NX=rho.get_NX();
	int NY=rho.get_NY();
	int u,v;
	for (u=0; u<NX; u++)
		for(v=0; v<NY; v++)
		    tem+=greenf[j](u,v)*rho(u,v);
	return tem;
	
}


void Greenfb::add_boundary_potential(Grid2D& rho)
{
	if (a0 > 0.0 && b0 > 0.0) return;
	int i,j, u=0;
	int NX=rho.get_NX();
	int NY=rho.get_NY();
	double dx=rho.get_dx();
	double dy=rho.get_dy();
	Grid2D rho_tmp(NX,NY,dx,dy);
	for(i=0; i<NX; i++)
		for(j=0; j<NY; j++)	 
		  rho_tmp(i,j)=rho(i,j);	
	for(i=1; i<NX; i++)
	   {
		rho(i,NY-1)+=eps0*get_potential_b(u,rho_tmp);  
		u=u+1;
	   }
	for(j=1; j<NY-1; j++)	
	   {  
        rho(NX-1,j)+=eps0*get_potential_b(u,rho_tmp); 
		u=u+1;
       }
	for(i=1; i<NX-1; i++)
	   {
		rho(i,1)+=eps0*get_potential_b(u,rho_tmp);
		u=u+1;
	   }
	for(j=2; j<NY-1; j++)
	   {	  
        rho(1,j)+=eps0*get_potential_b(u,rho_tmp);
		u=u+1;
       }
}


void poisson_xy(Grid2D& rho, Greenfb& gf)
{
 int j, l,NX=rho.get_NX(),NY=rho.get_NX();
 double dx=rho.get_dx(), dy=rho.get_dy();
 double *temx, *temy;
 temx=new double[NX];
 temy=new double[NY]; 

 gf.add_boundary_potential(rho);

 for(j=0; j<NX; j++)
   {  
    for(l=0;l<NY;l++)
     temy[l]=rho[l+NY*j];
    sinft(temy-1,NY);
    for(l=0;l<NY;l++)
     rho[l+NY*j]=temy[l];
   }
 
 for(l=0; l<NY; l++)
   {
    for(j=0;j<NX;j++)
     temx[j]=rho[l+j*NY];
    sinft(temx-1,NX);
    for(j=0;j<NX;j++)
     rho[l+NY*j]=temx[j];
   }
 
  for(j=1;j<NX;j++)
   for(l=1;l<NY;l++)
     rho[l+j*NY]*=-1.0/eps0             
                  *0.5*dx*dy/(cos(PI*j/NX)+cos(PI*l/NY)-2.0); 
		  
 for(j=0; j<NX; j++)
  {
   for(l=0;l<NY;l++)
    temy[l]=rho[l+NY*j];   
   sinft(temy-1,NY);
   for(l=0;l<NY;l++)
     rho[l+NY*j]=2.0*temy[l]/NY;   
  }
 
 for(l=0; l<NY; l++)
   {
    for(j=0;j<NX;j++)
     temx[j]=rho[l+j*NY];
    sinft(temx-1,NX);
    for(j=0;j<NX;j++)
     rho[l+NY*j]=2.0*temx[j]/NX;
   }

 delete temx;
 delete temy;
 
}


/*
void poisson_xy(float* cdensity, int NX, int NY, float dx, float dy)
{
 int j, l;
 float *temx, *temy;
 temx=new float[NX];
 temy=new float[NY]; 

 for(j=0; j<NX; j++)
   {  
    for(l=0;l<NY;l++)
     temy[l]=cdensity[l+NY*j];
    sinft(temy-1,NY);
    for(l=0;l<NY;l++)
     cdensity[l+NY*j]=temy[l];
   }
 
 for(l=0; l<NY; l++)
   {
    for(j=0;j<NX;j++)
     temx[j]=cdensity[l+j*NY];
    sinft(temx-1,NX);
    for(j=0;j<NX;j++)
     cdensity[l+NY*j]=temx[j];
   }
 
  for(j=1;j<NX;j++)
   for(l=1;l<NY;l++)
     cdensity[l+j*NY]*=-1.0/eps0             
                  *0.5*dx*dy/(cos(PI*j/NX)+cos(PI*l/NY)-2.0); 
		  
 for(j=0; j<NX; j++)
  {
   for(l=0;l<NY;l++)
    temy[l]=cdensity[l+NY*j];   
   sinft(temy-1,NY);
   for(l=0;l<NY;l++)
     cdensity[l+NY*j]=2.0*temy[l]/NY;   
  }
 
 for(l=0; l<NY; l++)
   {
    for(j=0;j<NX;j++)
     temx[j]=cdensity[l+j*NY];
    sinft(temx-1,NX);
    for(j=0;j<NX;j++)
     cdensity[l+NY*j]=2.0*temx[j]/NX;
   }

 delete temx;
 delete temy;
 
}

*/

//! calculates electric field

void poisson_xy(Grid2D& Ex, Grid2D& Ey, Grid2D& rho, Greenfb& gf)
{
 int NX=rho.get_NX(),NY=rho.get_NX();
 double dx=rho.get_dx(), dy=rho.get_dy();
 poisson_xy(rho,gf);  
 for(int j=3;j<NX-2;j++)
  for(int l=3;l<NY-2;l++)
   {
    Ex(j,l)=-(rho(j+1,l)-rho(j-1,l))/(2.0*dx);
    Ey(j,l)=-(rho(j,l+1)-rho(j,l-1))/(2.0*dy);
   }
}


void poisson_xyz(Grid3D& Ex, Grid3D& Ey, Grid3D& rho, Greenfb& gf)
{
 int NX=rho.get_NX(),NY=rho.get_NX(), NZ=rho.get_NZ();
 double dx=rho.get_dx(), dy=rho.get_dy();
 for(int i=0; i<NZ; i++)
 {
  poisson_xy(rho[i],gf);  
  for(int j=1;j<NX-1;j++)
   for(int l=1;l<NY-1;l++)
    {
     Ex[i](j,l)=-(rho[i](j+1,l)-rho[i](j-1,l))/(2.0*dx);
     Ey[i](j,l)=-(rho[i](j,l+1)-rho[i](j,l-1))/(2.0*dy);
    }
 }
}


// needed for poisson_rz

void poisson_r(double *q, int m, int NR, int NZ, double dr, double dz)
 {
  int j;
  double *a,*b,*c,*r,alpha,beta;
  alpha=2.0*(cos(2.0*PI*m/(NZ))-1.0)/pow(dz,2);
  beta=4.0/pow(dr,2);
  a = new double[NR]; 
  b = new double[NR]; 
  c = new double[NR]; 
  r = new double[NR]; 
  for(j=1; j<NR; j++)
    {
     a[j]=(1.0-0.5/((double)j))/pow(dr,2);
     c[j]=(1.0+0.5/((double)j))/pow(dr,2);
     b[j]=2.0*((cos(2.0*PI*m/(NZ))-1.0)/pow(dz,2)-1.0/pow(dr,2));
     r[j]=-q[j]/eps0;
    }
  b[0]=alpha-beta;
  c[0]=beta; 
  r[0]=-q[0]/eps0; 
  tridiag(q,a,b,c,r,NR);
  delete a;
  delete b;
  delete c;
  delete r;
 }


// r,z poisson solver with FFT in z and tridiag in r.

void poisson_rz(Grid2D& rho)
{
 int j,m;
 int NR=rho.get_NX();
 int NZ=rho.get_NY();
 double dr=rho.get_dx();
 double dz=rho.get_dy();

 double *q1=new double[NR];
 double *q2=new double[NR];
 double *tem=new double[NZ];

 for(j=0; j<NR; j++)
   {
    for(m=0;m<NZ;m++)
     tem[m]=rho[m+NZ*j];
    realft(tem-1,NZ,1);
    //cosft1(tem-1,NZ);
    //sinft(tem-1,NZ);
    for(m=0;m<NZ;m++)
     rho[m+NZ*j]=tem[m];
   }

 for(m=2;m<NZ;m+=2)
    {
     for(j=0;j<NR;j++)
      {
       q1[j]=rho[m+NZ*j];
       q2[j]=rho[m+1+NZ*j];
      }       
     poisson_r(q1,(int)(0.5*m),NR,NZ,dr,dz);
     poisson_r(q2,(int)(0.5*m),NR,NZ,dr,dz);
     for(j=0;j<NR;j++)
      {
       rho[m+NZ*j]=q1[j];
       rho[m+1+NZ*j]=q2[j];
      }    
    }

  for(j=0;j<NR;j++)
       q1[j]=rho[1+NZ*j];
  poisson_r(q1,(int)(0.5*NZ),NR,NZ,dr,dz);
  for(j=0;j<NR;j++)
     rho[1+NZ*j]=q1[j];


  for(j=0;j<NR;j++)               
    q1[j]=rho[j*NZ];          
  poisson_r(q1,0,NR,NZ,dr,dz);      


  for(j=0;j<NR;j++)
    rho[0+NZ*j]=q1[j]; 


  for(j=0; j<NR; j++)
    { 
     for(m=0;m<NZ;m++)
      tem[m]=rho[m+NZ*j];
     realft(tem-1,NZ,-1);
     //cosft1(tem-1,NZ);
     //sinft(tem-1,NZ);
     for(m=0;m<NZ;m++)
      rho[m+NZ*j]=2.0*tem[m]/NZ;     
    }

  delete q1;
  delete q2;
  delete tem;
}





