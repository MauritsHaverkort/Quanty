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
#include "BasicMath/BasicMath_BlokTriDiagonalMatrixVector.c"
#include "BasicMath/BasicMath_BlokDoubleTriDiagonalMatrixVector.c"
#include "BasicMath/BasicMath_AndersonMatrixVector.c"
#include "BasicMath/BasicMath_BlokAndersonMatrixVector.c"
#include "BasicMath/BasicMath_ListOfPoles.c"
#include "BasicMath/BasicMath_MatrixTypeConversion.c"
#include "BasicMath/BasicMath_GreenFunctionArithmetic.c"
#include "BasicMath/BasicMath_GreenFunctionMatrixArithmetic.c"

/*
 *   x^a with both x and a integers
 */
uint32_t PowerInt(uint32_t x,uint32_t a)
{
  uint32_t i;
  uint32_t retval=1;
  for(i=0;i<a;i++)
  {
    retval*=x;
  }
  return retval;
}

/*
*  z = x cross y for vectors of lenght 3
*/
void Cross(const double x[], const double y[], double z[])
{
  z[0]=x[1]*y[2]-y[1]*x[2];
  z[1]=x[2]*y[0]-y[2]*x[0];
  z[2]=x[0]*y[1]-y[0]*x[1];
}

/*
 *
 */
double ParallelalPipedVolume(const double x[], const double y[], const double z[])
{
  double t[3];
  Cross(x,y,t);
  return Dot(z,t);
}



void Swap(double *const a, double *const b)
{
  double c;
  c=*a;
  *a=*b;
  *b=c;
}

double CoVariance(double x[], double y[], int N)
{
  double AvverageX=0, AvverageY=0,Result=0;
  int i;
  for(i=0;i<=N;i++)
  {
    AvverageX+=x[i];
    AvverageY+=y[i];
  }
  AvverageX/=((double)(N+1));
  AvverageY/=((double)(N+1));
  for(i=0;i<=N;i++)
  {
    Result+=(x[i]-AvverageX)*(y[i]-AvverageY);
  }
  return Result;
}

double Variance(double x[], int N)
{
  double AvverageX=0,Result=0;
  int i;
  for(i=0;i<=N;i++)
  {
    AvverageX+=x[i];
  }
  AvverageX/=((double)(N+1));
  for(i=0;i<=N;i++)
  {
    Result+=Sqr(x[i]-AvverageX);
  }
  return Result;
}

double Avverage(double x[], int N)
{
  int i;
  double Result=0;
  for(i=0;i<=N;i++)
  {
    Result+=x[i];
  }
  Result/=((double)(N+1));
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
  if(Variance(x,N)>(10E-10)*Avverage(x,N))
  {
    *a=CoVariance(x,y,N)/Variance(x,N);
    *b=Avverage(y,N)-(*a)*Avverage(x,N);
  }
  else
  {
    printf("Linear Interpolate on a set of dependent variables return a constant\n");
    *a=0;
    *b=Avverage(y,N);
  }
}
/*
*  given a set of points find a quadratic polinomial y=ax^2+bx+c that minimizes
*  \sum_i Sqr(y_i-a x_i^2 - b x_i -x)
*/
void QuadraticInterpolate(double x[], double y[], int N, double *a, double *b, double *c)
{
  double SumX=0,SumXX=0,SumXXX=0,SumXXXX=0,SumY=0,SumXY=0,SumXXY=0,denom;
  int i;
  for(i=0;i<=N;i++)
  {
    SumX+=x[i];
    SumXX+=x[i]*x[i];
    SumXXX+=x[i]*x[i]*x[i];
    SumXXXX+=x[i]*x[i]*x[i]*x[i];
    SumY+=y[i];
    SumXY+=x[i]*y[i];
    SumXXY+=x[i]*x[i]*y[i];
  }
  denom=SumXX*SumXX*SumXX-2*SumX*SumXX*SumXXX+(N+1)*SumXXX*SumXXX+SumX*SumX*SumXXXX-(N+1)*SumXX*SumXXXX;
  if(fabs(denom)>10E-10)
  {
    *a=(SumX*SumX*SumXXY-(N+1)*SumXX*SumXXY-SumX*SumXX*SumXY+(N+1)*SumXXX*SumXY+SumXX*SumXX*SumY-SumX*SumXXX*SumY)/denom;
    *b=-(SumX*SumXX*SumXXY-(N+1)*SumXXY*SumXXX-SumXX*SumXX*SumXY+(N+1)*SumXXXX*SumXY+SumXX*SumXXX*SumY-SumX*SumXXXX*SumY)/denom;
    *c=-(-SumXX*SumXX*SumXXY+SumX*SumXXY*SumXXX+SumXX*SumXXX*SumXY-SumX*SumXXXX*SumXY-SumXXX*SumXXX*SumY+SumXX*SumXXXX*SumY)/denom;
  }
  else
  {
    printf("Quadratic Interpolate on a set of dependent variables return linear interpolation\n");
    *a=0;
    LinearInterpolate(x, y, N, b, c);
  }
}

void ResiprocalLatticeVector(double r[3][3], double g[3][3])
{
  double V;
  Cross(r[1],r[2],g[0]);
  Cross(r[2],r[0],g[1]);
  Cross(r[0],r[1],g[2]);
  V=Dot(r[0],g[0]);
  g[0][0]*=2*Pi/V;
  g[0][1]*=2*Pi/V;
  g[0][2]*=2*Pi/V;
  g[1][0]*=2*Pi/V;
  g[1][1]*=2*Pi/V;
  g[1][2]*=2*Pi/V;
  g[2][0]*=2*Pi/V;
  g[2][1]*=2*Pi/V;
  g[2][2]*=2*Pi/V;
}
