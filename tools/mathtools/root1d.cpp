#include "Definitions.h"
#include "mathtools.h"

#define MAXIT 100
#define ITMAX 100
#define EPS 3.0E-8
#define SIGN(a,b) ((b)>=0.0 ? fabs(a) : -fabs(b))

/* Root Searching, Numerical Recipes p. 366  */


double mathtools::rtsafe(void (*funcd)(double,double*,double*),
                         double x1,double x2,double xacc)
 {
  int j;
  double df,dx,dxold,f,fh,fl;
  double temp,xh,xl,rts;
  (*funcd)(x1,&fl,&df);
  (*funcd)(x2,&fh,&df);
  if((fl>0.0 && fh > 0.0) || (fl<0.0 && fh<0.0))
    printf("error\n");
  if(fl==0.0) return x1;
  if(fh==0.0) return x2;
  if(fl<0.0) {
    xl=x1;
    xh=x2;
  } else {
    xh=x1;
    xl=x2;
  }
  rts=0.5*(x1+x2);
  dxold=fabs(x2-x1);
  dx=dxold;
  (*funcd)(rts,&f,&df);
  for(j=1;j<=MAXIT;j++) {
   if((((rts-xh)*df-f)*((rts-xl)*df-f) >= 0.0) 
      || (fabs(2.0*f) > fabs(dxold*df))) {
      dxold=dx;
      dx=0.5*(xh-xl);
      rts=xl+dx;
      if(xl==rts) return rts;
   } else {
      dxold=dx;
      dx=f/df;
      temp=rts;
      rts-=dx;
      if(temp==rts) return rts;
   }
   if(fabs(dx)<xacc) return rts;
   (*funcd)(rts,&f,&df);
   if(f<0.0)
    xl=rts;
   else
    xh=rts;
  }
  printf("maxit!\n");
  return(0.0);
 }      
  
/* Root Searching, Numerical Recipes p. 361  */

/*

double zbrent(double (*func)(double),double x1,double x2,double tol)
 {
  int iter;
  double a=x1,b=x2,c=x2,e,min1,min2;
  double fa=(*func)(a),fb=(*func)(b),fc,p,q,r,s,tol1,xm;
  if((fa>0.0 && fb>0.0) || (fa<0.0 && fb<0.0))
    printf("error\n");
  fc=fb;
  for(iter=1;iter<=ITMAX;iter++) {
    if((fb>0.0 && fc>0.0) || (fb<0.0 && fc<0.0)) {
      c=a;
      fc=fa;
      e=d=b-a;
    }
    if(fabs(fc)<fabs(fb)) {
      a=b;
      b=c;
      c=a;
      fa=fb;
      fb=fc;
      fc=fa;
    }
    tol1=2.0*EPS*fabs(b)+0.5*tol;
    xm=0.5*(c-b);
    if(fabs(xm) <= tol1 || fb==0.0) return b;
    if(fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
      s=fb/fa;
      if(a==c) {
         p=2.0*xm*s;
         q=1.0-s;
      } else {
         q=fa/fc;
         r=fb/fc;
         p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
         q=(q-1.0)*(r-1.0)*(s-1.0);
      }
      if(p>0.0) q=-q;
      p=fabs(p);
      min1=3.0*xm*q-fabs(tol1*q);
      min2=fabs(e*q);
      if(2.0*p < (min1 < min2 ? min1 : min2)) {
         e=d;
         d=p/q;
      } else {
         d=xm;
         e=d;
      }
   } else {
      d=xm;
      e=d;
     }   
     a=b;
     fa=fb;
     if(fabs(d) > tol1)
       b+=d;
     else
       b+=SIGN(tol1,xm);
     fb=(*func)(b);
  }
  printf("error2\n");
  return(0.0);
 }

 */




