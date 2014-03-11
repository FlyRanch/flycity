/* Subroutine to integrate the tangent vector equations of a planar geometric curve
whos bend angle is parameterized by a bspline curve.  
	
	Name: modelInt.c
	Language: C
	Subroutines called: none
	
    b1        = first term of the basis function recursion relation
    b2        = second term of the basis function recursion relation
    c         = order of the B-spline basis function
    d1[]      = array containing the derivative of the basis functions
                d1[1]) contains the derivative of the basis function associated with B1 etc.
    d2[]      = array containing the derivative of the basis functions
                d2[1] contains the derivative of the basis function associated with B1 etc.
    f1        = first term of the first derivative of the basis function recursion relation
    f2        = second term of the first derivative of the basis function recursion relation
    f3        = third term of the first derivative of the basis function recursion relation
    f4        = fourth term of the first derivative of the basis function recursion relation
    npts      = number of defining polygon vertices
    n[]       = array containing the basis functions
                n[1]) contains the basis function associated with B1 etc.
    nplusc    = constant -- npts + c -- maximum knot value
    s1        = first term of the second derivative of the basis function recursion relation
    s2        = second term of the second derivative of the basis function recursion relation
    s3        = third term of the second derivative of the basis function recursion relation
    s4        = fourth term of the second derivative of the basis function recursion relation
    t         = parameter value
    temp[]    = temporary array
    x[]       = knot vector
*/	

#include	<stdio.h>
#include        <stdlib.h>
#include	<math.h>
#include        "mex.h"
#include        "matrix.h"


/* Input Arguments */
#define P_IN    prhs[0]
#define	T_IN	prhs[1]
#define	C_IN	prhs[2]
#define	X_IN	prhs[3]

/* Output Arguments */
#define	X_OUT	plhs[0]



#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

/*----------------------------------------------------*/
void dbasisP(double N[],
	     double *c,
	     double T,
	     int npts,
	     double x[])

{
  int nplusc,numofpts;
  int j,k;
  double t;
  double b1,b2;
  double *temp;
    
  nplusc = npts + *c;
  numofpts = npts;
  
  /* Don't forget that C index starts from 0! */

  t = T;
  
  /* allows for number of defining polygon vertices = numofpts */
  temp = mxCalloc((nplusc-1),sizeof(double));		
  
  /* calculate the first order basis functions n[j] */
  for (j = 0; j <= nplusc-1-1; j++){
    if (( t >= x[j]) && (t < x[j+1]))
      temp[j] = 1;
    else
      temp[j] = 0;
  }
  
  if (t == (double)x[numofpts-1]){		/*    pick up last point	*/
    temp[numofpts-1] = 1;
    temp[numofpts] = 0;
  }
  
  /* calculate the higher order basis functions */
    
  for (k = 2; k <= *c; k++){
    for (j = 0; j <= nplusc-k-1; j++){
      if (temp[j] != 0)    /* if the lower order basis function is zero skip the calculation */
	b1 = ((t-x[j]) * (temp[j])) / (x[j+k-1] - x[j]);
      else
	b1 = 0;
      
      if (temp[j+1] != 0)     /* if the lower order basis function is zero skip the calculation */
	b2 = ((x[j+k]-t) * (temp[j+1])) / (x[j+k] - x[j+1]);
      else
	b2 = 0;
      
      temp[j] = b1 + b2;
    }
  }
  
  /* put in i,j array	*/
  for (j = 0; j <= numofpts-1; j++) {
    N[j] = temp[j];
  }
  
  mxFree(temp);
  
  return;
} 

/*----------------------------------------------------*/
double tangfunc1(double T,
		 double p[],
		 double *c,
		 int npts,
		 double x[])
{
  void dbasisP(double N[], double *c, double T, int npts, double x[]);
  double a=0.0,V;
  double *N;
  int n;
  N = mxCalloc((npts-1),sizeof(double));
  dbasisP(N,c,T,npts,x);
  
  /* Perform inner product between spline control points and basis function*/
  for (n=0; n <= npts-1; n++){
    a += N[n] * p[n];
  }
  mxFree(N);
  V = cos(a);
  return V;
}

/*----------------------------------------------------*/
double tangfunc2(double T,
		 double p[],
		 double *c,
		 int npts,
		 double x[])
{
  void dbasisP(double N[], double *c, double T, int npts, double x[]);
  double a=0.0,V;
  double *N;
  int n;
  N = mxCalloc((npts-1),sizeof(double));
  dbasisP(N,c,T,npts,x);
  
  /* Perform inner product between spline control points and basis function*/
  for (n=0; n <= npts-1; n++){
    a += N[n] * p[n];
  }
  mxFree(N);
  V = sin(a);
  return V;
}
/*----------------------------------------------------*/
static double *ss = NULL;

void cleanup(void) {
  mxFree(ss);
}
/*----------------------------------------------------*/
/* #define FUNC(x,a1,a2,a3,a4) ((*func)(x,a1,a2,a3,a4))*/

void polint(double xa[], double ya[], int n, double x, double *y, double *dy)
{
    int i,m,ns = 1;
    double den,dif,dift,ho,hp,w;
    double *c,*d;
    
    dif = fabs(x - xa[1]);
    c = mxCalloc(n,sizeof(double));
    d = mxCalloc(n,sizeof(double));
    /*c = vector(1,n);
    d = vector(1,n);*/
    
    for (i = 1; i <= n; i++){
        if ( (dift = fabs(x-xa[i])) < dif) {
            ns = i;
            dif = dift;
        }
        c[i] = ya[i];
        d[i] = ya[i];
    }
    
    *y = ya[ns--];
    for (m = 1; m < n; m++){
        for (i = 1; i < n-m; i++){
            ho = xa[i]-x;
            hp = xa[i+m] - x;
            w = c[i+1] - d[i];
            if ( (den = ho-hp) == 0.0) mexErrMsgTxt("Error in routine polint");
            den = w/den;
            d[i] = hp*den;
            c[i] = ho*den;
        }
        *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
    }
    /*free_vector(d,1,n);
    free_vector(c,1,n);*/
    mxFree(c);
    mxFree(d);
}
    

/*----------------------------------------------------*/
void trapzd(double (*func)(double, double *, double *, int, double *), 
	      double a, 
	      double b, 
	      int n,
	      double *ss,
	      double p[],
	      double *c,
	      int npts,
	      double xx[])
	      

/*This routine computes the nth stage of refinement of an extended 
trapezoidal rule. func is input as a pointer to the function to be 
integrated between limits a and b, also input. When called with n=1,
the routine returns the crudest estimate of b a f(x)dx. Subsequent calls 
with n=2,3,... (in that sequential order) will improve the accuracy by 
adding 2n-2 additional interior points. */

{
  double x,tnm,sum,del;
  /*static double ss;*/
  int it,j;
  if (n == 1) {
    (*ss) = 0.5*(b-a)*(((*func)(a,p,c,npts,xx)) + ((*func)(b,p,c,npts,xx)) );
    return;
  } 
  else {
    for (it=1,j=1;j<n-1;j++) it <<= 1;
    tnm = it;
    del = (b-a)/tnm; /*This is the spacing of the points to be added.*/
    x = a + 0.5*del;
    for (sum=0.0,j=1;j<=it;j++,x+=del) sum += ((*func)(x,p,c,npts,xx));
    (*ss) = 0.5*((*ss) + (b-a)*sum/tnm); /*This replaces s by its refined value.*/
    return; /*ss*/
  }
}
/*----------------------------------------------------*/

#define EPS 1.0e-5
#define JMAX 20
#define JMAXP (JMAX+1)
#define K 5
double qromb(double (*func)(double, double *, double *, int, double *), 
	     double a,
	     double b,
	     double *ss,
	     double p[],
	     double *c,
	     int npts,
	     double xx[])
/*Returns the integral of the function func from a to b. The parameters 
EPS can be set to the desired fractional accuracy and JMAX so that 2 to 
the power JMAX-1 is the maximum allowed number of steps. Integration is 
performed by Simpson's rule. */
{
    void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
    
    void trapzd(double (*func)(double, double *, double *, int, double *),
    double a,
    double b,
    int n,
    double *ss,
    double p[],
    double *c,
    int npts,
    double xx[]);
    
  /*void nrerror(char error_text[]);*/
    int j;
    double sss, dsss;
    double *s, *h;
    
    s = mxCalloc(JMAXP,sizeof(double));
    h = mxCalloc(JMAXP+1,sizeof(double));
    
    h[1] = 1.0;
    for (j=1; j<=JMAX; j++) {
        trapzd(func,a,b,j,ss,p,c,npts,xx);
        s[j] = *ss;
        if (j >= K) /*Avoid spurious early convergence.*/{
            polint(&h[j-K], &s[j-K], K, 0.0, &sss, &dsss);
            if (fabs(dsss) <= EPS * fabs(sss)) {
                /* This is a HACK.  I divide by 2 because for some reason, 
                 the 'polint' function seems to scale the answer by 2 */
                return sss/2;
            }
        }
        h[j+1] = 0.25*h[j];
    }
    mxFree(s);
    mxFree(h);
    /*nrerror("Too many steps in routine qsimp");*/
    return sss; /*Never get here.*/

}
/*----------------------------------------------------*/

void modelRomIntP(double X[],
	      double p[],
	      double *c,
	      double T[],
	      int npts,
	      double x[],
	      int sizeofT,
	      double *ss)
		   
{
  double tangfunc(double T, double p[], double *c, int npts, double x[]);
  double qromb(double (*func)(double, double *, double *, int, double *), 
	       double a,
	       double b,
	       double *ss,
	       double p[],
	       double *c,
	       int npts,
	       double xx[]);
  double t;
  int i;
  
  for (i = 0; i <= sizeofT - 1; i++){
    t = T[i];
    
    X[i] = qromb(tangfunc1,0,t,ss,p,c,npts,x);
    X[(i + sizeofT)] = qromb(tangfunc2,0,t,ss,p,c,npts,x);
  }
  return;  
}

/*----------------------------------------------------*/
/* mexFunction is the gateway routine for the MEX-file. */ 
void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
     
{ 
  double *X; 
  double *c,*T,*x,*p;
  int sizeofT, npts; 
  static double *ss;
  
  /* Check for proper number of arguments */
    
  if (nrhs != 4) { 
    mexErrMsgTxt("Five input arguments required."); 
  } else if (nlhs > 1) {
    mexErrMsgTxt("Too many output arguments."); 
  } 
  
  
  /* Assign pointers to the input parameters */ 
  c = mxGetPr(C_IN); 
  T = mxGetPr(T_IN);
  x = mxGetPr(X_IN);
  p = mxGetPr(P_IN);

  /* Basis function matrix is DIM = [length(t) npts]. t,npts are row
   vectors*/
  sizeofT = mxGetN(T_IN);
  npts = mxGetN(P_IN);
  

  /* Create a matrix for the return argument */ 
  X_OUT = mxCreateDoubleMatrix(sizeofT,2, mxREAL); 
    
   /* Assign pointers to the output parameters */ 
  X  = mxGetPr(X_OUT);
  
  ss = mxGetPr(mxCreateDoubleScalar(0));
  mexMakeMemoryPersistent(ss);

  /* Do the actual computations in a subroutine */
  modelRomIntP(X,p,c,T,npts,x,sizeofT,ss); 
  mexAtExit(cleanup);
  return;
  
}
