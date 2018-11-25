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
	CubicSpline(const std::vector<double> &u, const std::vector<double> &v);
	virtual ~CubicSpline();

	void Create(const std::vector<double> &u, const std::vector<double> &v);
	std::vector<double> GetY2();
	double GetY2(uint32_t i);

	double operator()(double t);
	double operator()(double t, unsigned int k);
	std::vector<double> operator()(const std::vector<double> &t);
	std::vector<double> operator()(const std::vector<double> &t, unsigned int k);
};

} /* namespace Interpolation */

#endif /* CUBICSPLINE_H_ */
