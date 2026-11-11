/*
 * DKBSpline.cpp
 *
 * Radial Dirac spinor components in the dual-kinetic-balance B-spline basis, see
 * DKBSpline.h. Skeleton: the evaluation follows with the radial Dirac solver (DFT stage
 * A1); until then every value is NaN.
 */

#include "DKBSpline.h"

namespace Interpolation {

DKBSpline::DKBSpline() {
	type = "DKBSpline";
	degree = 0;
	kappa = 0;
	component = 0;
	speedOfLight = 0.;
	nQ = 0;
}

DKBSpline::~DKBSpline() {
}

double DKBSpline::operator()(double t)
{
	return NAN;
}

double DKBSpline::operator()(double t, unsigned int k)
{
	return NAN;
}

std::vector<double> DKBSpline::operator()(const std::vector<double> &t)
{
	return std::vector<double>(t.size(), NAN);
}

std::vector<double> DKBSpline::operator()(const std::vector<double> &t, unsigned int k)
{
	return std::vector<double>(t.size(), NAN);
}

} /* namespace Interpolation */
