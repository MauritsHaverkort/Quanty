/* Cubic splines */

int InitCubicSpline(CubicSplineType *S)
{
  S->x  = (double*) malloc(S->N*sizeof(double));
  S->y  = (double*) malloc(S->N*sizeof(double));
  S->y2 = (double*) malloc(S->N*sizeof(double));
  return 0;
}

void FreeCubicSpline(CubicSplineType *S)
{
  free(S->x);
  free(S->y);
  free(S->y2);
  S->N=0;
}

void CreateCubicSpline(CubicSplineType *S, double yp1, double ypn)
{
  CreateCubicSpline(S->x, S->y, S->N, yp1, ypn, S->y2);
}

void CreateCubicSpline(double *x, double *y, uint32_t n, double yp1, double ypn, double *y2)
{
  int i,k;
  double p,qn,sig,un,u[n];
  if (isnan(yp1))
  {
    y2[0]=u[0]=0.0;
  }
  else
  {
    y2[0] = -0.5;
    u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
  }
  for (i=1;i<n-1;i++)
  {
    sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p=sig*y2[i-1]+2.0;
    y2[i]=(sig-1.0)/p;
    u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]); u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }
  if (isnan(ypn))
  {
    qn=un=0.0;
  }
  else
  {
    qn=0.5;
    un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
  }
  y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
  for (k=n-2;k>=0;k--)
  {
    y2[k]=y2[k]*y2[k+1]+u[k];
  }
}

double EvaluateCubicSpline(CubicSplineType *S, double x)
{
  return EvaluateCubicSpline(S->x, S->y, S->y2, S->N, x);
}

double EvaluateCubicSpline(double *xa, double *ya, double *y2a, uint32_t n, double x)
{
  int klo,khi,k;
  double h,b,a;
  klo=0;
  khi=n-1;
  while (khi-klo > 1)
  {
    k=(khi+klo) >> 1;
    if (xa[k] > x) khi=k; else klo=k;
  }
  h=xa[khi]-xa[klo];
  if (h == 0.0)
  {
    return ((ya[khi]+ya[klo])/2);
  }
  else
  {
    a=(xa[khi]-x)/h;
    b=(x-xa[klo])/h;
    return a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
  }
}
