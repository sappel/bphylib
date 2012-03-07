/* Newton-Raphson Method  S.379 Num. Rec.  */

#include<math.h>
#include<stdlib.h>
#include<stdio.h>

#define TINY 1.0e-20

extern "C" void mnewt(int ntrial,double *x,int n,float tolx,float tolf);

void usrfun(double *x,int n,double *fvec,double **fjac);

void ludcmp(double **a,int n,int *indx,double *d)
 {
  int i,imax,j,k;
  double big,dum,sum,temp;
  double *vv;
  vv=(double*)malloc(sizeof(double)*n);
  *d=1.0;
  for(i=0;i<n;i++) {
   big=0.0;
   for(j=0;j<n;j++)
    if ((temp=fabs(a[i][j])) > big) big=temp;
   if(big==0.0) printf("error\n");
   vv[i]=1.0/big;
  }
  for(j=0;j<n;j++) {
   for(i=0;i<j;i++) {
    sum=a[i][j];
    for(k=0;k<i;k++) sum-=a[i][k]*a[k][j];
    a[i][j]=sum;
   }
   big=0.0;
   for(i=j;i<n;i++) {
     sum=a[i][j];
     for(k=1;k<j;k++)
       sum-=a[i][k]*a[k][j];
     a[i][j]=sum;
     if((dum=vv[i]*fabs(sum)) >= big) {
      big=dum;
      imax=i;
     }
   }
   if(j != imax) {
    for(k=0;k<n;k++) {
      dum=a[imax][k];
      a[imax][k]=a[j][k];
      a[j][k]=dum;
    }
    *d=-(*d);
    vv[imax]=vv[j];
   }   
   indx[j]=imax;
   if(a[j][j]==0.0) a[j][j]=TINY;
   if(j != n-1) {
    dum=1.0/(a[j][j]);
    for(i=j+1;i<n;i++) a[i][j]*=dum;
   } 
  }
  free(vv);
 }


void lubksb(double **a,int n,int *indx,double *b)
 {
  int i,ii=-1,ip,j;
  double sum;
  for(i=0;i<n;i++) {
   ip=indx[i];
   sum=b[ip];
   b[ip]=b[i];
   if(ii+1)
      for(j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
   else if(sum) ii=i;
   b[i]=sum;
  }
  for(i=n-1;i>=0;i--) {
     sum=b[i];
     for(j=i+1;j<n;j++) sum -= a[i][j]*b[j];
     b[i]=sum/a[i][i];
   }
 }
  

void mnewt(int ntrial,double *x,int n,float tolx,float tolf)
 {
  int k,i,*indx;
  double errx,errf,d,*fvec,**fjac,*p;
  indx=(int*)malloc(sizeof(int)*n);
  fvec=(double*)malloc(sizeof(double)*n);
  p=(double*)malloc(sizeof(double)*n);
  fjac=(double**)malloc(sizeof(double*)*n);
  for(i=0;i<n;i++)
   fjac[i]=(double*)malloc(sizeof(double)*n);
  for(k=0;k<ntrial;k++) {
    usrfun(x,n,fvec,fjac);
    errf=0.0;
    for(i=0;i<n;i++) errf += fabs(fvec[i]);
    if(errf <= tolf) {
      free(indx);free(fvec);free(fjac);free(p);printf("1\n");return;}
    for(i=0;i<n;i++) p[i]=-fvec[i];
    ludcmp(fjac,n,indx,&d);
    lubksb(fjac,n,indx,p);
    errx=0.0;
    for(i=0;i<n;i++) {  
      errx+=fabs(p[i]);
      x[i] += p[i];
    }
    if(errx <= tolx) {
     free(indx);free(fvec);free(fjac);free(p);printf("2\n");return;}
   }
   free(indx);free(fvec);free(fjac);free(p);return;
 }








