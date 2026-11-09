/*
 * TBSpline.h
 *
 *  Created on: Sep 10, 2019
 *      Author: brass
 */

#ifndef MARTIN_MATH_TBSPLINE_H_
#define MARTIN_MATH_TBSPLINE_H_

#include "BSpline.h"
#include "InterpolatingFunction.h"

namespace BSpline {

class TBSpline: public Interpolation::InterpolatingFunction {
public:
	TBSpline();
	TBSpline(const std::vector<double> &u,const std::vector<double> &k,const std::vector<double> &v, uint32_t n);
	virtual ~TBSpline();

	std::vector<double> knots;
	uint32_t degree;
	std::vector<double> points;  // the points the spline interpolates; arithmetic samples it there

	void Create(const std::vector<double> &u, const std::vector<double> &v,uint32_t n);
	double operator()(double t);
	double operator()(double t, unsigned int k);
	std::vector<double> operator()(const std::vector<double> &t);
	std::vector<double> operator()(const std::vector<double> &t, unsigned int k);
};

void UnaryMinus(TBSpline &f, TBSpline &g);
void Add(TBSpline &f, Interpolation::InterpolatingFunction &g, TBSpline &r);
void Sub(TBSpline &f, Interpolation::InterpolatingFunction &g, TBSpline &r);
void Mul(TBSpline &f, Interpolation::InterpolatingFunction &g, TBSpline &r);
void Div(TBSpline &f, Interpolation::InterpolatingFunction &g, TBSpline &r);
void Pow(TBSpline &f, double n, TBSpline &r);
void Concat(TBSpline &f, Interpolation::InterpolatingFunction &g, TBSpline &r);
bool operator==(TBSpline &f, TBSpline &g);

} /* namespace BSpline */

#endif /* MARTIN_MATH_TBSPLINE_H_ */
