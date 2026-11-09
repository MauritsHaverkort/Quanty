/*
 * BSplineMatrix.h
 *
 * Integrals over B-splines in plain C: quadrature points per knot interval, banded
 * matrices <B_i^(a)| f |B_j^(b)>, and the cell tables and assembly of radial Slater
 * integrals of functions given by B-spline coefficients. See BSplineMatrix.cpp.
 *
 * The matrix builders and the Slater cell method are Martin Brass' (2019, namespace
 * HartreeFock in Martin/Hartree-Fock). Rewritten in 2026 on the B-spline core, with the
 * same quadrature and the same table layouts; Martin's functions are now wrappers.
 */

#ifndef CORE_BASICMATH_INTERPOLATION_BSPLINEMATRIX_H_
#define CORE_BASICMATH_INTERPOLATION_BSPLINEMATRIX_H_

int BSplineQuadratureSetup(const double *knot, uint32_t nKnot, uint32_t nQ, double **rq, double **wq, uint32_t *nPoint);

int BSplineMatrixBanded(const double *knot, uint32_t nKnot, uint32_t degree,
                        const double *rq, const double *wq, uint32_t nPoint, const double *fq,
                        uint32_t a, uint32_t b, uint32_t first, uint32_t n, double *band, uint32_t ldBand);

int BSplineCellValues(const double *knot, uint32_t nKnot, uint32_t degree, uint32_t nQ, double *BSval);

int BSplineSlaterCellTables(const double *knot, uint32_t nKnot, uint32_t degree, uint32_t nQ, const double *BSval,
                            int k, double *rk1, double *rk2, double *rkD);

double BSplineSlaterIntegral(const double *ci1, const double *cj1, const double *ci2, const double *cj2,
                             const double *rk1, const double *rk2, const double *rkD, uint32_t nKnot, uint32_t degree);

#endif /* CORE_BASICMATH_INTERPOLATION_BSPLINEMATRIX_H_ */
