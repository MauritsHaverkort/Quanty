/*
 * CubicHermiteSpline.h
 *
 *  Created on: Oct 23, 2018
 *      Author: brass
 */

#ifndef CORE_BASICMATH_INTERPOLATION_CUBICHERMITESPLINE_H_
#define CORE_BASICMATH_INTERPOLATION_CUBICHERMITESPLINE_H_

#include "InterpolatingFunction.h"

namespace Interpolation {

class CubicHermiteSpline: public InterpolatingFunction {
private:
	std::vector<double> m;
public:
	CubicHermiteSpline();
	CubicHermiteSpline(const std::vector<double> &u, const std::vector<double> &v);
	CubicHermiteSpline(const std::vector<double> &u, const std::vector<double> &v, const std::vector<double> &w);
	virtual ~CubicHermiteSpline();

	void Create(const std::vector<double> &u, const std::vector<double> &v);
	void Create(const std::vector<double> &u, const std::vector<double> &v, const std::vector<double> &w);

	double operator()(double t);
	double operator()(double t, unsigned int k);
	std::vector<double> operator()(const std::vector<double> &t);
	std::vector<double> operator()(const std::vector<double> &t, unsigned int k);
};

} /* namespace Interpolation */

#endif /* CORE_BASICMATH_INTERPOLATION_CUBICHERMITESPLINE_H_ */
