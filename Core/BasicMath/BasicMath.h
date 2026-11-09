// #include "BasicMath_Constants.h" no header file for constants, only defines
#include "BasicMath_Complex.h"
#include "BasicMath_StandardFunctions.h"
#include "Interpolation/GaussLegendre.h"
//#include "BasicMath_InterpolatingFunction.h"
//#include "BasicMath_Spline.h"
#include "BasicMath_Combinatorial.h"
#include "BasicMath_VectorLinearAlgabra.h"
#include "BasicMath_MatrixVector.h"
#include "BasicMath_CompactMatrixVector.h"
#include "BasicMath_CompactTensor.h"
#include "BasicMath_SparseMatrixVector.h"
#include "BasicMath_TriDiagonalMatrixVector.h"
#include "BasicMath_DoubleTriDiagonalMatrixVector.h"
#include "BasicMath_BADoubleTriDiagonalMatrixVector.h"
#include "BasicMath_BlokTriDiagonalMatrixVector.h"
#include "BasicMath_BlokDoubleTriDiagonalMatrixVector.h"
#include "BasicMath_AndersonMatrixVector.h"
#include "BasicMath_BlokAndersonMatrixVector.h"
#include "BasicMath_ListOfPoles.h"
#include "BasicMath_MatrixTypeConversion.h"
#include "BasicMath_GreenFunctionArithmetic.h"

#include "Interpolation/InterpolatingFunction.h"
#include "Interpolation/CubicSpline.h"
#include "Interpolation/CubicHermiteSpline.h"
#include "Interpolation/BSpline.h"
#include "Interpolation/TBSpline.h"




uint32_t PowerInt(uint32_t x, uint32_t a);
double PowerDoubleInt(double x, uint32_t a);
int pow2roundup (int x);
void Cross(const double x[], const double y[], double z[]);
int NormalizeVector(double v[3], double epsilon);
double AngleBetweenVectors(const double a[3], const double b[3]);
double ParallelalPipedVolume(const double x[], const double y[], const double z[]);
int Invert3x3(const double M[3][3], double Minv[3][3]);
void Swap(double *const a, double *const b);
double CoVariance(double x[], double y[], int N);
double Variance(double x[], int N);
double Avverage(double x[], int N);
void LinearInterpolate(double x[], double y[], int N, double *a, double *b);
void QuadraticInterpolate(double x[], double y[], int N, double *a, double *b, double *c);
double InterpolateRoot(double const *x, double const *y, uint32_t N);
double InterpolateRootOfset(double const *x, double const *y, uint32_t N, double const ofset);
void InterpolateRoot2D(double const *x1, double const *x2, double const *y1, double const *y2, uint32_t N, double *const x1next, double *const x2next);
void InterpolateRoot2DOfset(double const *x1, double const *x2, double const *y1, double const *y2, uint32_t N, double const ofsety1, double const ofsety2, double *const x1next, double *const x2next);
double FindSelfconsistentnumber(double const *valin, double const *valout, uint32_t N);
void ReciprocalLatticeVector(double r[3][3], double g[3][3], uint32_t const Units);
int findPower(double num, int power);
int findPower(double num);
double round(double num, double precision);
bool matrixSorterByFirstValue(vector<double> v1, vector<double> v2);
bool Test_uint8_t_array_in_range(const uint8_t min, const uint8_t max, uint8_t *l, const uint32_t N);
bool Test_uint16_t_array_in_range(const uint16_t min, const uint16_t max, uint16_t *l, const uint32_t N);
bool Test_uint32_t_array_in_range(const uint32_t min, const uint32_t max, uint32_t *l, const uint32_t N);
template <typename T> bool isElementOf(const T element, const vector<T> &set);
