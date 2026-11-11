/*
 * DKBSpline.h
 *
 * One component of a radial Dirac spinor in the dual-kinetic-balance (DKB) B-spline
 * basis of the DFT/HF add-on, stage A (Doc/Claude/QuantyDFT_StageA, section 5 and 10).
 *
 * With the working B-splines B_i (all but the first and the last of the knot set) and
 * h = 1/(2c), the 2n DKB functions are
 *     u_i^(1) = ( B_i , h (B_i' + kappa B_i/r) ),   u_i^(2) = ( h (B_i' - kappa B_i/r) , B_i ),
 * and a spinor sum_i a_i u_i^(1) + b_i u_i^(2) has the components
 *     P = S_a + h (S_b' - kappa S_b/r)        (large, "G")
 *     Q = h (S_a' + kappa S_a/r) + S_b        (small, "F")
 * with the splines S_a = sum_i a_i B_i and S_b = sum_i b_i B_i. P and Q are no pure
 * B-spline expansions (they contain B' and B/r) and are always evaluated through these
 * formulas. Like the other B-splines they are 0 outside the knots, value and derivatives.
 *
 * The class keeps its own copy of the knots, kappa, c and the coefficients, so that every
 * Lua object owns its data. x holds the quadrature points of the knot set, y the values
 * of the function on them, so that code which samples a function on GetX() keeps working.
 */

#ifndef CORE_BASICMATH_INTERPOLATION_DKBSPLINE_H_
#define CORE_BASICMATH_INTERPOLATION_DKBSPLINE_H_

#include "InterpolatingFunction.h"

namespace Interpolation {

class DKBSpline: public InterpolatingFunction {
public:
	DKBSpline();
	virtual ~DKBSpline();

	std::vector<double> knots;
	uint32_t degree;
	int kappa;
	uint32_t component;          // 0 the large component P (G), 1 the small component Q (F)
	double speedOfLight;         // c of the balance h = 1/(2c)
	std::vector<double> a, b;    // coefficients of S_a and S_b over all nKnot-degree-1 B-splines,
	                             // the first and the last are 0
	uint32_t nQ;                 // Gauss points per knot interval of x and of Integrate

	double operator()(double t);
	double operator()(double t, unsigned int k);
	std::vector<double> operator()(const std::vector<double> &t);
	std::vector<double> operator()(const std::vector<double> &t, unsigned int k);
};

} /* namespace Interpolation */

#endif /* CORE_BASICMATH_INTERPOLATION_DKBSPLINE_H_ */
