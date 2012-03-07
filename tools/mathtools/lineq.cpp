#include "Definitions.h"
#include "mathtools.h"

// solve tridiagonal system (from Num. Rec. p. )


void mathtools::tridiag(float* x,float* a,float* b,float* c,float* r,int n)
 {
  int j;
  double bet,*gam;
  gam=new double[n];
  x[0]=r[0]/(bet=b[0]);
  for(j=1;j<n;j++)
    {
      gam[j]=c[j-1]/bet;
      bet=b[j]-a[j]*gam[j];
      x[j]=(r[j]-a[j]*x[j-1])/bet;
    }
  for(j=(n-2);j>=0;j--)
    x[j]-=gam[j+1]*x[j+1];
  delete gam;
 }


void mathtools::tridiag(double* x,double* a,double* b,double* c,double* r,int n)
 {
  int j;
  double bet,*gam;
  gam=new double[n];
  x[0]=r[0]/(bet=b[0]);
  for(j=1;j<n;j++)
    {
      gam[j]=c[j-1]/bet;
      bet=b[j]-a[j]*gam[j];
      x[j]=(r[j]-a[j]*x[j-1])/bet;
    }
  for(j=(n-2);j>=0;j--)
    x[j]-=gam[j+1]*x[j+1];
  delete gam;
 }
