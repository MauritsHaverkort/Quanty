#define DistanceSqr(a,b) ((a[0]-b[0])*(a[0]-b[0]))+((a[1]-b[1])*(a[1]-b[1]))+((a[2]-b[2])*(a[2]-b[2]))
#define Distance(a,b)    (sqrt((a[0]-b[0])*(a[0]-b[0]))+((a[1]-b[1])*(a[1]-b[1]))+((a[2]-b[2])*(a[2]-b[2])))
#define Min(a,b) ((a) < (b) ? (a) : (b))
#define Max(a,b) ((a) < (b) ? (b) : (a))
#define Abs(a)   ((a) < (0) ? (-(a)) : (a))
#define Sqr(a) ((a)*(a))
#define Cubed(a) ((a)*(a)*(a))
#define Dot(a,b) ((a[0])*(b[0])+(a[1])*(b[1])+(a[2])*(b[2]))
#define Crossx(a,b) ((a[1])*(b[2])-(b[1])*(a[2]))
#define Crossy(a,b) ((a[2])*(b[0])-(b[2])*(a[0]))
#define Crossz(a,b) ((a[0])*(b[1])-(b[0])*(a[1]))
#define ComplexTimes2Re(x1,y1,x2,y2) (x1*x2-y1*y2)
#define ComplexTimes2Im(x1,y1,x2,y2) (x1*y2+y1*x2)
#define ComplexTimes3Re(x1,y1,x2,y2,x3,y3) (x1*x2*x3-x1*y2*y3-y1*x2*y3-y1*y2*x3)
#define ComplexTimes3Im(x1,y1,x2,y2,x3,y3) (y1*x2*x3+x1*y2*x3+x1*x2*y3-y1*y2*y3)
#define ComplexTimes4Re(x1,y1,x2,y2,x3,y3,x4,y4) (ComplexTimes2Re(x1,y1,x2,y2)*ComplexTimes2Re(x3,y3,x4,y4)-ComplexTimes2Im(x1,y1,x2,y2)*ComplexTimes2Im(x3,y3,x4,y4))
#define ComplexTimes4Im(x1,y1,x2,y2,x3,y3,x4,y4) (ComplexTimes2Re(x1,y1,x2,y2)*ComplexTimes2Im(x3,y3,x4,y4)+ComplexTimes2Im(x1,y1,x2,y2)*ComplexTimes2Re(x3,y3,x4,y4))
#define ComplexTimes5Re(x1,y1,x2,y2,x3,y3,x4,y4,x5,y5) (ComplexTimes2Re(x1,y1,x2,y2)*ComplexTimes3Re(x3,y3,x4,y4,x5,y5)-ComplexTimes2Im(x1,y1,x2,y2)*ComplexTimes3Im(x3,y3,x4,y4,x5,y5))
#define ComplexTimes5Im(x1,y1,x2,y2,x3,y3,x4,y4,x5,y5) (ComplexTimes2Re(x1,y1,x2,y2)*ComplexTimes3Im(x3,y3,x4,y4,x5,y5)+ComplexTimes2Im(x1,y1,x2,y2)*ComplexTimes3Re(x3,y3,x4,y4,x5,y5))
#define ComplexTimes6Re(x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6) (ComplexTimes2Re(x1,y1,x2,y2)*ComplexTimes4Re(x3,y3,x4,y4,x5,y5,x6,y6)-ComplexTimes2Im(x1,y1,x2,y2)*ComplexTimes4Im(x3,y3,x4,y4,x5,y5,x6,y6))
#define ComplexTimes6Im(x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6) (ComplexTimes2Im(x1,y1,x2,y2)*ComplexTimes4Re(x3,y3,x4,y4,x5,y5,x6,y6)+ComplexTimes2Im(x1,y1,x2,y2)*ComplexTimes4Re(x3,y3,x4,y4,x5,y5,x6,y6))
#define ComplexTimes7Re(x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,x7,y7) (ComplexTimes2Re(x1,y1,x2,y2)*ComplexTimes5Re(x3,y3,x4,y4,x5,y5,x6,y6,x7,y7)-ComplexTimes2Im(x1,y1,x2,y2)*ComplexTimes5Im(x3,y3,x4,y4,x5,y5,x6,y6,x7,y7))
#define ComplexTimes7Im(x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,x7,y7) (ComplexTimes2Im(x1,y1,x2,y2)*ComplexTimes5Re(x3,y3,x4,y4,x5,y5,x6,y6,x7,y7)+ComplexTimes2Im(x1,y1,x2,y2)*ComplexTimes5Re(x3,y3,x4,y4,x5,y5,x6,y6,x7,y7))
#define ParallelalPipedVolumeMatrix(a) Abs(a[0][0]*Crossx(a[1],a[2])+a[0][1]*Crossy(a[1],a[2])+a[0][2]*Crossz(a[1],a[2]))

#include "BasicMath/BasicMath_Constants.c"
#include "BasicMath/BasicMath_StandardFunctions.c"
#include "BasicMath/BasicMath_Combinatorial.c"
#include "BasicMath/BasicMath_VectorLinearAlgabra.c"
#include "BasicMath/BasicMath_MatrixVector.c"
#include "BasicMath/BasicMath_SparseMatrixVector.c"
#include "BasicMath/BasicMath_TriDiagonalMatrixVector.c"
#include "BasicMath/BasicMath_DoubleTriDiagonalMatrixVector.c"
#include "BasicMath/BasicMath_BADoubleTriDiagonalMatrixVector.c"
#include "BasicMath/BasicMath_BlokTriDiagonalMatrixVector.c"
#include "BasicMath/BasicMath_BlokDoubleTriDiagonalMatrixVector.c"
#include "BasicMath/BasicMath_AndersonMatrixVector.c"
#include "BasicMath/BasicMath_BlokAndersonMatrixVector.c"
#include "BasicMath/BasicMath_ListOfPoles.c"
#include "BasicMath/BasicMath_MatrixTypeConversion.c"
#include "BasicMath/BasicMath_GreenFunctionArithmetic.c"

/*
 *   x^a with both x and a integers
 */
uint32_t PowerInt(uint32_t x, uint32_t a)
{
  uint32_t i;
  uint32_t retval = 1;
  for (i = 0; i < a; i++)
  {
    retval *= x;
  }
  return retval;
}

/*
*  z = x cross y for vectors of lenght 3
*/
void Cross(const double x[], const double y[], double z[])
{
  z[0] = x[1] * y[2] - y[1] * x[2];
  z[1] = x[2] * y[0] - y[2] * x[0];
  z[2] = x[0] * y[1] - y[0] * x[1];
}

/*
 *
 */
double ParallelalPipedVolume(const double x[], const double y[], const double z[])
{
  double t[3];
  Cross(x, y, t);
  return Dot(z, t);
}



void Swap(double *const a, double *const b)
{
  double c;
  c = *a;
  *a = *b;
  *b = c;
}

double CoVariance(double x[], double y[], int N)
{
  double AvverageX = 0, AvverageY = 0, Result = 0;
  int    i;
  for (i = 0; i <= N; i++)
  {
    AvverageX += x[i];
    AvverageY += y[i];
  }
  AvverageX /= ((double)(N + 1));
  AvverageY /= ((double)(N + 1));
  for (i = 0; i <= N; i++)
  {
    Result += (x[i] - AvverageX) * (y[i] - AvverageY);
  }
  return Result;
}

double Variance(double x[], int N)
{
  double AvverageX = 0, Result = 0;
  int    i;
  for (i = 0; i <= N; i++)
  {
    AvverageX += x[i];
  }
  AvverageX /= ((double)(N + 1));
  for (i = 0; i <= N; i++)
  {
    Result += Sqr(x[i] - AvverageX);
  }
  return Result;
}

double Avverage(double x[], int N)
{
  int    i;
  double Result = 0;
  for (i = 0; i <= N; i++)
  {
    Result += x[i];
  }
  Result /= ((double)(N + 1));
  return Result;
}


/*
*   given a set of points find a linear function y=ax+b that minimizes
*   \sum_i Sqr(y_i-a x_i - b)
*
*   Error     =\sum_i Sqr(y_i-a x_i - b)
*   dError/da =\sum_i -2 x_i (y_i - a x_i - b) = 0
*   dError/db =\sum_i -2 (y_i - a x_i -b)      = 0 
*
*/
void LinearInterpolate(double x[], double y[], int N, double *a, double *b)
{
  if (Variance(x, N) > (10E-10) * Avverage(x, N))
  {
    *a = CoVariance(x, y, N) / Variance(x, N);
    *b = Avverage(y, N) - (*a) * Avverage(x, N);
  }
  else
  {
    printf("Linear Interpolate on a set of dependent variables return a constant\n");
    *a = 0;
    *b = Avverage(y, N);
  }
}

/*
*  given a set of points find a quadratic polinomial y=ax^2+bx+c that minimizes
*  \sum_i Sqr(y_i-a x_i^2 - b x_i -x)
*/
void QuadraticInterpolate(double x[], double y[], int N, double *a, double *b, double *c)
{
  double SumX = 0, SumXX = 0, SumXXX = 0, SumXXXX = 0, SumY = 0, SumXY = 0, SumXXY = 0, denom;
  int    i;
  for (i = 0; i <= N; i++)
  {
    SumX += x[i];
    SumXX += x[i] * x[i];
    SumXXX += x[i] * x[i] * x[i];
    SumXXXX += x[i] * x[i] * x[i] * x[i];
    SumY += y[i];
    SumXY += x[i] * y[i];
    SumXXY += x[i] * x[i] * y[i];
  }
  denom =
    SumXX * SumXX * SumXX - 2 * SumX * SumXX * SumXXX + (N + 1) * SumXXX * SumXXX + SumX * SumX * SumXXXX - (N +
                                                                                                             1) *
    SumXX * SumXXXX;
  if (fabs(denom) > 10E-10)
  {
    *a =
      (SumX * SumX * SumXXY - (N + 1) * SumXX * SumXXY - SumX * SumXX * SumXY + (N + 1) * SumXXX * SumXY +
       SumXX * SumXX * SumY - SumX * SumXXX * SumY) / denom;
    *b =
      -(SumX * SumXX * SumXXY - (N + 1) * SumXXY * SumXXX - SumXX * SumXX * SumXY + (N + 1) * SumXXXX * SumXY +
        SumXX * SumXXX * SumY - SumX * SumXXXX * SumY) / denom;
    *c =
      -(-SumXX * SumXX * SumXXY + SumX * SumXXY * SumXXX + SumXX * SumXXX * SumXY - SumX * SumXXXX * SumXY -
        SumXXX * SumXXX * SumY + SumXX * SumXXXX * SumY) / denom;
  }
  else
  {
    printf("Quadratic Interpolate on a set of dependent variables return linear interpolation\n");
    *a = 0;
    LinearInterpolate(x, y, N, b, c);
  }
}

/*
 *
 */
double InterpolateRoot(double const *x, double const *y, uint32_t N)
{
  double xsort[N], ysort[N], xret;
  uint32_t i, j, k, firstpos;
  if (N == 1)
    return x[0] + y[0];
  xsort[0] = x[0];
  ysort[0] = y[0];
  for (i = 1; i < N; i++)
  {
    ysort[i] = y[i];
    xsort[i] = x[i];
    for (j = 0; j < i; j++)
    {
      if (y[i] < ysort[j])
      {
        for (k = i; k > j; k--)
        {
          ysort[k] = ysort[k - 1];
          xsort[k] = xsort[k - 1];
        }
        ysort[j] = y[i];
        xsort[j] = x[i];
        j = i;
      }
    }
  }
  firstpos = 0;
  for (i = 0; i < N; i++)
  {
    if (y[i] < 0)
      firstpos++;
  }
  if (firstpos == 0)
  {
    if (ysort[1] == ysort[0])
      return xsort[0];
    return (xsort[0] - ysort[0] * (xsort[1] - xsort[0]) / (ysort[1] - ysort[0]));
  }
  else if (firstpos == N)
  {
    if (ysort[N - 1] == ysort[N - 2])
      return xsort[N - 1];
    return (xsort[N - 1] - ysort[N - 1] * (xsort[N - 1] - xsort[N - 2]) / (ysort[N - 1] - ysort[N - 2]));
  }
  else
  {
    xret =
      (xsort[firstpos - 1] * ysort[firstpos] - xsort[firstpos] * ysort[firstpos - 1]) / (ysort[firstpos] -
                                                                                         ysort[firstpos - 1]);
    if ((xret < xsort[firstpos - 1]) || (xret > xsort[firstpos]))
    {
      xret = (xsort[firstpos] + xsort[firstpos - 1]) / 2.0;
    }
    if (xret == x[N - 1])
    {
      xret = (xsort[firstpos] + xsort[firstpos - 1]) / 2.0;
    }
    return xret;
  }
}

double InterpolateRootOfset(double const *x, double const *y, uint32_t N, double const ofset)
{
  double yoff[N];
  uint32_t i;
  for (i = 0; i < N; i++)
  {
    yoff[i] = y[i] - ofset;
  }
  return InterpolateRoot(x, yoff, N);
}

void InterpolateRoot2D(double const *x1, double const *x2, double const *y1, double const *y2, uint32_t N,
                       double *const x1next, double *const x2next)
{
  uint32_t const maxit1 = 10;
  uint32_t i;
  double x2step[N], y2step[N];
  if (N == 1)
  {
    *x1next = x1[0] + y1[0];
    *x2next = x2[0];
    return;
  }
  if ((N % maxit1) == 0)
  {
    /* step in second direction */
    printf("\n\n\nNew loop\n\n\n");
    for (i = maxit1 - 1; i < N; i += maxit1)
    {
      x2step[((i + 1) / maxit1) - 1] = x2[i];
      y2step[((i + 1) / maxit1) - 1] = y2[i];
      printf("%22.15E %22.15E %i\n", x2[i], y2[i], ((i + 1) / maxit1) - 1);
    }
    *x1next = x1[N - 1];
    *x2next = InterpolateRoot(x2step, y2step, (N / maxit1));
    printf("%22.15E \n\n\n\n", *x2next);
  }
  else
  {
    *x2next = x2[N - 1];
    *x1next = InterpolateRoot(&(x1[N - (N % maxit1)]), &(y1[N - (N % maxit1)]), N % maxit1);
  }
  return;
}


void InterpolateRoot2DOfset(double const *x1, double const *x2, double const *y1, double const *y2, uint32_t N,
                            double const ofsety1, double const ofsety2, double *const x1next, double *const x2next)
{
  double y1off[N], y2off[N];
  uint32_t i;
  for (i = 0; i < N; i++)
  {
    y1off[i] = y1[i] - ofsety1;
    y2off[i] = y2[i] - ofsety2;
  }
  InterpolateRoot2D(x1, x2, y1off, y2off, N, x1next, x2next);
  return;
}


/*
 *
 */
double FindSelfconsistentnumber(double const *valin, double const *valout, uint32_t N)
{
  double valdif[N];
  uint32_t i;
  for (i = 0; i < N; i++)
  {
    valdif[i] = valin[i] - valout[i];
  }
  return InterpolateRoot(valin, valdif, N);
}



void ResiprocalLatticeVector(double r[3][3], double g[3][3])
{
  double V;
  Cross(r[1], r[2], g[0]);
  Cross(r[2], r[0], g[1]);
  Cross(r[0], r[1], g[2]);
  V = Dot(r[0], g[0]);
  g[0][0] *= 2 * Pi / V;
  g[0][1] *= 2 * Pi / V;
  g[0][2] *= 2 * Pi / V;
  g[1][0] *= 2 * Pi / V;
  g[1][1] *= 2 * Pi / V;
  g[1][2] *= 2 * Pi / V;
  g[2][0] *= 2 * Pi / V;
  g[2][1] *= 2 * Pi / V;
  g[2][2] *= 2 * Pi / V;
}
