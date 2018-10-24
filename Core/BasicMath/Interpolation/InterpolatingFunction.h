/*
 * InterpolatingFunction.h
 *
 *  Created on: Sep 8, 2018
 *      Author: brass
 */

#ifndef INTERPOLATINGFUNCTION_H_
#define INTERPOLATINGFUNCTION_H_

#include <vector>

namespace Interpolation {

class InterpolatingFunction {
protected:
	std::vector<double> x;
	std::vector<double> y;
public:

	InterpolatingFunction();
	virtual ~InterpolatingFunction();

	std::vector<double> GetX();
	std::vector<double> GetY();
	double GetX(uint32_t i);
	double GetY(uint32_t i);

	virtual void Create(const std::vector<double> &u, const std::vector<double> &v);

	virtual double operator()(double t);
	virtual double operator()(double t, unsigned int k);
	virtual std::vector<double> operator()(const std::vector<double> &t);
	virtual std::vector<double> operator()(const std::vector<double> &t, unsigned int k);

};

} /* namespace Interpolation */

#endif /* INTERPOLATINGFUNCTION_H_ */
