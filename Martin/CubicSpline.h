/*
 * CubicSpline.h
 *
 *  Created on: Sep 8, 2018
 *      Author: brass
 */

#ifndef CUBICSPLINE_H_
#define CUBICSPLINE_H_

#include "InterpolatingFunction.h"

namespace Interpolation {

class CubicSpline: public InterpolatingFunction {
private:
	std::vector<double> y2;
public:
	CubicSpline();
	CubicSpline(std::vector<double> u, std::vector<double> v);
	virtual ~CubicSpline();

	double operator()(double t);
	double operator()(double t, unsigned int k);
};

} /* namespace Interpolation */

#endif /* CUBICSPLINE_H_ */
