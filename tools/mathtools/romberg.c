#include<math.h>
#include<stdlib.h>
#include<stdio.h>

#define EPS 1.0e-6
#define JMAX 20

/*
This routine computes the nth stage of refinement of an extended trapezoidal rule. func is input
as a pointer to the function to be integrated between limits a and b, also input. When called with
n=1, the routine returns the crudest estimate of  b
a f(x)dx. Subsequent calls with n=2,3,...
(in that sequential order) will improve the accuracy by adding 2n-2 additional interior points.
*/

double trapzd(double (*func)(double), double a, double b, int n)
{
double x,tnm,sum,del;
static double s;
int it,j;
if (n == 1) {
return (s=0.5*(b-a)*(((*func)(a))+((*func)(b))));
} else {
for (it=1,j=1;j<n-1;j++) it <<= 1;
tnm=it;
del=(b-a)/tnm; // This is the spacing of the points to be added.
x=a+0.5*del;
for (sum=0.0,j=1;j<=it;j++,x+=del) sum += ((*func)(x));
s=0.5*(s+(b-a)*sum/tnm); // This replaces s by its refined value.
return s;
}
}

/*
Returns the integral of the function func from a to b. The parameters EPS can be set to the
desired fractional accuracy and JMAX so that 2 to the power JMAX-1 is the maximum allowed
number of steps. Integration is performed by SimpsonÕs rule.
*/

double qsimp(double (*func)(double), double a, double b)
{
int j;
double s,st,ost=0.0,os=0.0;
for (j=1;j<=JMAX;j++) {
st=trapzd(func,a,b,j);
s=(4.0*st-ost)/3.0; // Compare equation (4.2.4), above.
if (j > 5) // Avoid spurious early convergence.
if (fabs(s-os) < EPS*fabs(os) ||
(s == 0.0 && os == 0.0)) return s;
os=s;
ost=st;
}
return 0.0;
}


