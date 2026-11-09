/*
 * BSpline.h
 *
 * B-splines in plain C. See BSpline.cpp and Doc/Claude/BSPLINE_DESIGN.md.
 *
 * Conventions:
 *  - degree, not order (order = degree + 1), as Lua "Degree";
 *  - knot[0 .. nKnot-1] non-decreasing; there are nKnot - degree - 1 basis functions,
 *    B_i has support [knot[i], knot[i+degree+1]];
 *  - a point on an interior knot belongs to the interval to its right, the last knot to
 *    the last non-empty interval;
 *  - outside [knot[0], knot[nKnot-1]] every basis function, and every derivative, is 0.
 *
 * Martin Brass (2019) wrote the first version in C++ (namespace BSpline). Rewritten in
 * 2026 with plain arrays and de Boor's algorithm; the class TBSpline stays in TBSpline.h.
 */

#ifndef CORE_BASICMATH_INTERPOLATION_BSPLINE_H_
#define CORE_BASICMATH_INTERPOLATION_BSPLINE_H_

#define BSplineMaxDegree 30

uint32_t BSplineFindInterval(double t, const double *knot, uint32_t nKnot);

int BSplineValues(double t, uint32_t interval, const double *knot, uint32_t nKnot, uint32_t degree,
                  uint32_t nDeriv, double *val);

double BSplineEvaluate(uint32_t i, uint32_t degree, const double *knot, uint32_t nKnot, double t);

double BSplineEvaluateDerivative(uint32_t i, uint32_t degree, const double *knot, uint32_t nKnot, double t, uint32_t k);

double BSplineDerivative(double t, const double *knot, uint32_t nKnot, uint32_t degree, const double *c, uint32_t k);

double BSplineDeBoor(uint32_t idx, double t, const double *knot, uint32_t degree, const double *c);

int BSplineKnotsAveraged(const double *t, uint32_t nt, uint32_t degree, double *knot);

int BSplineCollocationMatrix(const double *knot, uint32_t nKnot, uint32_t degree, const double *t, uint32_t nt, double *B);

int BSplineKnotsFischer(int m, double hMax, double rMax, double Z, uint32_t degree, double **knot, uint32_t *nKnot);

#endif /* CORE_BASICMATH_INTERPOLATION_BSPLINE_H_ */
