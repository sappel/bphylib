#include "Definitions.h"
#include "mathtools.h"

// 
// Spline interpolation routines
// from Num. Rec. p. 115-116 
//

void tridiag_special(double *x, double a, double b, double c, double *r, int n)
 {
  int j;
  double bet,*gam;
  gam=new double[n];
  x[0]=r[0]/(bet=b);
  for(j=1;j<n;j++)
    {
      gam[j]=c/bet;
      bet=b-a*gam[j];
      x[j]=(r[j]-a*x[j-1])/bet;
    }
  for(j=(n-2);j>=0;j--)
    x[j]-=gam[j+1]*x[j+1];
  delete gam;
 }


//
// calculate second derivatives:
//

void mathtools::spline(double *ys, double *y, int n)
 {
  int j;
  double *r=new double[n];

  for(j=1; j<n-1; j++)
    r[j]=6.0*(y[j+1]-2.0*y[j]+y[j-1]);    
  j=0;
  r[j]=6.0*(y[j+1]-2.0*y[j]); 
  j=n-1;
  r[j]=6.0*(-2.0*y[j]+y[j-1]); 

  double a=1.0, b=4.0;
  tridiag_special(ys,a,b,a,r,n);

  delete r;

 }


//
// search in an ordered table:
//


double mathtools::splint(double *ya, double *ys, double *xa, double x, int n)
{
  int klo, khi; 
  double a, b, dx=xa[1]-xa[0]; 

  klo=(int)floor((x-xa[0])/dx); 
  khi=(int)ceil((x-xa[0])/dx);
  if ( klo < 0 || khi > n-1 ) return 0.0;
  if (klo==khi) return ya[klo];
  
  a=(xa[khi]-x)/dx; 
  b=(x-xa[klo])/dx;
  return a*ya[klo]+b*ya[khi]+((a*a*a-a)*ys[klo]+(b*b*b-b)*ys[khi])/6.0;
}


