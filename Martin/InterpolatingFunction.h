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
	InterpolatingFunction(std::vector<double>, std::vector<double>);

	virtual double operator()(double t);
	virtual double operator()(double t, unsigned int k);

};

} /* namespace Interpolation */

#endif /* INTERPOLATINGFUNCTION_H_ */
