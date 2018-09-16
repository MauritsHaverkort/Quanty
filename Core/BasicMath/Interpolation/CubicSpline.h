/*
 * CubicSpline.h
 *
 *  Created on: Sep 8, 2018
 *      Author: brass
 */

#ifndef CUBICSPLINE_H_
#define CUBICSPLINE_H_

namespace Interpolation {

class CubicSpline: public InterpolatingFunction {
private:
	std::vector<double> y2;
public:

	CubicSpline();
	CubicSpline(std::vector<double> u, std::vector<double> v);
	virtual ~CubicSpline();

	void Create(std::vector<double>, std::vector<double>);
	std::vector<double> GetY2();
	double GetY2(uint32_t i);

	double operator()(double t);
	double operator()(double t, unsigned int k);
	std::vector<double> operator()(std::vector<double> t);
};

} /* namespace Interpolation */

#endif /* CUBICSPLINE_H_ */
