#include "Definitions.h"
#include "mathtools.h"


/*-------------------------C++ FFT Interfaces----------------------------*/


void mathtools::four1(vector<double>& data,int isign)
 {
  int n=data.size();
  four1(&data[0]-1,n/2,isign);
  if(isign==-1)
   for(int j=0;j<n;j++)
     data[j]*=2.0/(double)n;
 }  


void mathtools::realft(vector<double>& data,int isign)
 {
  int n=data.size();
  realft(&data[0]-1,n,isign);
  if(isign==1)
   for(int j=0;j<n;j++)
     data[j]*=2.0/n;
 }   


/*-------------------------FFT aus Num. Rec.----------------------------*/

#define SWAP(a,b) tempr=(a) ; (a)=(b) ; (b)=tempr
#define PI 3.14159265359 

void mathtools::four1(double *data,unsigned long nn,int isign)
 {
  unsigned long n,mmax,m,j,istep,i;
  double wtemp,wr,wpr,wpi,wi,theta;
  double tempr,tempi;
  n=nn << 1;
  j=1;
  for(i=1;i<n;i+=2) {
     if(j>i) {
         SWAP(data[j],data[i]);
         SWAP(data[j+1],data[i+1]);
        }
  m=n >> 1;
  while(m>=2 && j>m) {
     j-=m;
     m>>=1;
  }
  j+=m;
 }

mmax=2;
while(n>mmax){
  istep=mmax<<1;
  theta=isign*(6.28318530717959/mmax);
  wtemp=sin(0.5*theta);
  wpr=-2.0*wtemp*wtemp;
  wpi=sin(theta);
  wr=1.0;
  wi=0.0;
  for(m=1;m<mmax;m+=2) {
   for(i=m;i<=n;i+=istep) {
    j=i+mmax;
    tempr=wr*data[j]-wi*data[j+1];
    tempi=wr*data[j+1]+wi*data[j];
    data[j]=data[i]-tempr;
    data[j+1]=data[i+1]-tempi;
    data[i]+=tempr;
    data[i+1]+=tempi;
   }
   wr=(wtemp=wr)*wpr-wi*wpi+wr;
   wi=wi*wpr+wtemp*wpi+wi;
  }
  mmax=istep;
 }
}

void mathtools::four1(float *data,unsigned long nn,int isign)
 {
  unsigned long n,mmax,m,j,istep,i;
  float wtemp,wr,wpr,wpi,wi,theta;
  float tempr,tempi;
  n=nn << 1;
  j=1;
  for(i=1;i<n;i+=2) {
     if(j>i) {
         SWAP(data[j],data[i]);
         SWAP(data[j+1],data[i+1]);
        }
  m=n >> 1;
  while(m>=2 && j>m) {
     j-=m;
     m>>=1;
  }
  j+=m;
 }

mmax=2;
while(n>mmax){
  istep=mmax<<1;
  theta=isign*(6.28318530717959/mmax);
  wtemp=sin(0.5*theta);
  wpr=-2.0*wtemp*wtemp;
  wpi=sin(theta);
  wr=1.0;
  wi=0.0;
  for(m=1;m<mmax;m+=2) {
   for(i=m;i<=n;i+=istep) {
    j=i+mmax;
    tempr=wr*data[j]-wi*data[j+1];
    tempi=wr*data[j+1]+wi*data[j];
    data[j]=data[i]-tempr;
    data[j+1]=data[i+1]-tempi;
    data[i]+=tempr;
    data[i+1]+=tempi;
   }
   wr=(wtemp=wr)*wpr-wi*wpi+wr;
   wi=wi*wpr+wtemp*wpi+wi;
  }
  mmax=istep;
 }
}


void mathtools::realft(double *data,unsigned long n,int isign)
 {
  unsigned long i,i1,i2,i3,i4,np3;
  double c1=0.5,c2,h1r,h1i,h2r,h2i;
  double wr,wi,wpr,wpi,wtemp,theta;
  theta=PI/(double) (n>>1);
  if(isign == 1) {
    c2=-0.5;
    four1(data,n>>1,1);
 } else {
    c2=0.5;
    theta=-theta;
   }
  wtemp=sin(0.5*theta);
  wpr=-2.0*wtemp*wtemp;
  wpi=sin(theta);
  wr=1.0+wpr;
  wi=wpi;
  np3=n+3;
  for(i=2;i<=(n>>2);i++) {
     i4=1+(i3=np3-(i2=1+(i1=i+i-1)));
     h1r=c1*(data[i1]+data[i3]);
     h1i=c1*(data[i2]-data[i4]);
     h2r=-c2*(data[i2]+data[i4]);
     h2i=c2*(data[i1]-data[i3]);
     data[i1]=h1r+wr*h2r-wi*h2i;
     data[i2]=h1i+wr*h2i+wi*h2r;
     data[i3]=h1r-wr*h2r+wi*h2i;
     data[i4]=-h1i+wr*h2i+wi*h2r;
     wr=(wtemp=wr)*wpr-wi*wpi+wr;
     wi=wi*wpr+wtemp*wpi+wi;
  }
  if(isign==1) {
    data[1]=(h1r=data[1])+data[2];
    data[2]=h1r-data[2];
  } else {
    data[1]=c1*((h1r=data[1])+data[2]);
    data[2]=c1*(h1r-data[2]);
    four1(data,n>>1,-1);
  }
 }


void mathtools::realft(float *data, unsigned long n,int isign)
 {
  unsigned long i,i1,i2,i3,i4,np3;
  float c1=0.5,c2,h1r,h1i,h2r,h2i;
  float wr,wi,wpr,wpi,wtemp,theta;
  theta=PI/(float) (n>>1);
  if(isign == 1) {
    c2=-0.5;
    four1(data,n>>1,1);
 } else {
    c2=0.5;
    theta=-theta;
   }
  wtemp=sin(0.5*theta);
  wpr=-2.0*wtemp*wtemp;
  wpi=sin(theta);
  wr=1.0+wpr;
  wi=wpi;
  np3=n+3;
  for(i=2;i<=(n>>2);i++) {
     i4=1+(i3=np3-(i2=1+(i1=i+i-1)));
     h1r=c1*(data[i1]+data[i3]);
     h1i=c1*(data[i2]-data[i4]);
     h2r=-c2*(data[i2]+data[i4]);
     h2i=c2*(data[i1]-data[i3]);
     data[i1]=h1r+wr*h2r-wi*h2i;
     data[i2]=h1i+wr*h2i+wi*h2r;
     data[i3]=h1r-wr*h2r+wi*h2i;
     data[i4]=-h1i+wr*h2i+wi*h2r;
     wr=(wtemp=wr)*wpr-wi*wpi+wr;
     wi=wi*wpr+wtemp*wpi+wi;
  }
  if(isign==1) {
    data[1]=(h1r=data[1])+data[2];
    data[2]=h1r-data[2];
  } else {
    data[1]=c1*((h1r=data[1])+data[2]);
    data[2]=c1*(h1r-data[2]);
    four1(data,n>>1,-1);
  }
 }

 
void mathtools::sinft(double *y, unsigned long n)
 {
  int j,n2=n+2;
  double sum,y1,y2;
  double theta,wi=0.0,wpi,wpr,wr=1.0,wtemp;
  theta=PI/((double)n);
  wtemp=sin(0.5*theta);
  wpr=-2.0*wtemp*wtemp;
  wpi=sin(theta);
  y[1]=0.0;
  for(j=2;j<=(n>>1)+1;j++) {
     wr=(wtemp=wr)*wpr-wi*wpi+wr;
     wi=wi*wpr+wtemp*wpi+wi;
     y1=wi*(y[j]+y[n2-j]);
     y2=0.5*(y[j]-y[n2-j]);
     y[j]=y1+y2;
     y[n2-j]=y1-y2;
   }
  realft(y,n,1);
  y[1]*=0.5;
  sum=y[2]=0.0;
  for(j=1;j<=n-1;j+=2) {
    sum+=y[j];
    y[j]=y[j+1];
    y[j+1]=sum;
  }
 }


void mathtools::sinft(float *y,int n)
 {
  int j,n2=n+2;
  float sum,y1,y2;
  float theta,wi=0.0,wpi,wpr,wr=1.0,wtemp;
  theta=PI/((float)n);
  wtemp=sin(0.5*theta);
  wpr=-2.0*wtemp*wtemp;
  wpi=sin(theta);
  y[1]=0.0;
  for(j=2;j<=(n>>1)+1;j++) {
     wr=(wtemp=wr)*wpr-wi*wpi+wr;
     wi=wi*wpr+wtemp*wpi+wi;
     y1=wi*(y[j]+y[n2-j]);
     y2=0.5*(y[j]-y[n2-j]);
     y[j]=y1+y2;
     y[n2-j]=y1-y2;
   }
  realft(y,n,1);
  y[1]*=0.5;
  sum=y[2]=0.0;
  for(j=1;j<=n-1;j+=2) {
    sum+=y[j];
    y[j]=y[j+1];
    y[j+1]=sum;
  }
 }


void mathtools::cosft1(float *y,int n)
 {
  int j,n2;
  float sum,y1,y2;
  float theta,wi=0.0,wpi,wpr,wr=1.0,wtemp;
  theta=PI/n;
  wtemp=sin(0.5*theta);
  wpr=-2.0*wtemp*wtemp;
  wpi=sin(theta);
  sum=0.5*(y[1]-y[n]);
  y[1]=0.5*(y[1]+y[n]);
  n2=n+2;
  for(j=2;j<=(n>>1);j++) {
     wr=(wtemp=wr)*wpr-wi*wpi+wr;
     wi=wi*wpr+wtemp*wpi+wi;
     y1=0.5*(y[j]+y[n2-j]);
     y2=(y[j]-y[n2-j]);
     y[j]=y1-wi*y2;
     y[n2-j]=y1+wi*y2;
     sum+=wr*y2;
   }
  realft(y,n,1);
  y[2]=sum;
  for(j=4;j<=n;j+=2) {
    sum+=y[j];
    y[j]=sum;
  }
 }

 void mathtools::cosft1(double *y,int n)
 {
  int j,n2;
  double sum,y1,y2;
  double theta,wi=0.0,wpi,wpr,wr=1.0,wtemp;
  theta=PI/n;
  wtemp=sin(0.5*theta);
  wpr=-2.0*wtemp*wtemp;
  wpi=sin(theta);
  sum=0.5*(y[1]-y[n]);
  y[1]=0.5*(y[1]+y[n]);
  n2=n+2;
  for(j=2;j<=(n>>1);j++) {
     wr=(wtemp=wr)*wpr-wi*wpi+wr;
     wi=wi*wpr+wtemp*wpi+wi;
     y1=0.5*(y[j]+y[n2-j]);
     y2=(y[j]-y[n2-j]);
     y[j]=y1-wi*y2;
     y[n2-j]=y1+wi*y2;
     sum+=wr*y2;
   }
  realft(y,n,1);
  y[2]=sum;
  for(j=4;j<=n;j+=2) {
    sum+=y[j];
    y[j]=sum;
  }
 }

 
