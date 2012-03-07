#include "Definitions.h"
#include "mathtools.h"
#include <stdlib.h>

void mathtools::rk4(double *y, int n, double x ,double h, double *yout,
         void (*derivs)(double,double *,double *))
 {
  int i;
  double  xh,hh,h6,*dym,*dyt,*yt,*dydx;

  dym=(double*)malloc(sizeof(double)*n); 
  dyt=(double*)malloc(sizeof(double)*n);
  yt=(double*)malloc(sizeof(double)*n);
  dydx=(double*)malloc(sizeof(double)*n);
  
  (*derivs)(x,y,dydx); 

  hh=h*0.5;
  h6=h/6.0;
  xh=x+hh;
  for(i=0;i<n;i++) yt[i]=y[i]+hh*dydx[i];
   (*derivs)(xh,yt,dyt);
  for(i=0;i<n;i++) yt[i]=y[i]+hh*dyt[i];
   (*derivs)(xh,yt,dym);
  for(i=0;i<n;i++) 
   {
    yt[i]=y[i]+h*dym[i];
    dym[i]+=dyt[i];
   }
  (*derivs)(x+h,yt,dyt);
  for(i=0;i<n;i++)
   yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);

  free(dym);
  free(dyt);
  free(yt);
  free(dydx);
 }


