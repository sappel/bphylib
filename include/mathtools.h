namespace mathtools {

// FFT:

void realft(float *data,unsigned long n,int isign);
void realft(double *data,unsigned long n,int isign);
void realft(vector<double>& data,int isign);
void four1(float *data,unsigned long nn,int isign);
void four1(double *data,unsigned long nn,int isign);
void four1(vector<double>& data,int isign);
void sinft(double *y, unsigned long n);
void sinft(float *y,int n);
void cosft1(float *y,int n);
void cosft1(double *y,int n);

// Runge-Kutta:

void rk4(double *y, int n, double x ,double h, double *yout,
         void (*derivs)(double,double *,double *));

// random number generators:

float ran1(long *idum);
float gasdev(long *idum); 

// spline interpolation:

void spline(double *ys, double *y, int n);
double splint(double *ya, double *ys, double *xa, double x, int n);

// root searching:

double rtsafe(void (*funcd)(double,double*,double*),
              double x1,double x2,double xacc);
extern "C" void mnewt(int ntrial,double *x,int n,float tolx,float tolf);

// linear equations:

void tridiag(float* x,float* a,float* b,float* c,float* r,int n);
void tridiag(double* x,double* a,double* b,double* c,double* r,int n);

// integration:

extern "C" double qsimp(double (*func)(double), double a, double b);

}

using namespace std;
using namespace mathtools;
