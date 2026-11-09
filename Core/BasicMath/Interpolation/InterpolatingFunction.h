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

public:

	InterpolatingFunction();
	virtual ~InterpolatingFunction();

	std::string type;
	std::vector<double> x;
	std::vector<double> y;
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

uint32_t GridInterval(const std::vector<double> &x, double t);
void MergeGrids(const std::vector<double> &x1,const std::vector<double> &x2, std::vector<double> &xout );
void Add(InterpolatingFunction &f, InterpolatingFunction &g, InterpolatingFunction &fplusg);
void Mul(InterpolatingFunction &f, InterpolatingFunction &g, InterpolatingFunction &fplusg);
void Sub(InterpolatingFunction &f, InterpolatingFunction &g, InterpolatingFunction &fplusg);
void Div(InterpolatingFunction &f, InterpolatingFunction &g, InterpolatingFunction &fplusg);
void Pow(InterpolatingFunction &f, double n, InterpolatingFunction &g);
void Concat(InterpolatingFunction &f, InterpolatingFunction &g, InterpolatingFunction &fg);
bool operator==(InterpolatingFunction &f,InterpolatingFunction &g);
double Integrate(InterpolatingFunction &f);
double Integrate(InterpolatingFunction &f,double a,double b);
double IntegrateInterval(InterpolatingFunction &f, uint32_t idx);

} /* namespace Interpolation */

#endif /* INTERPOLATINGFUNCTION_H_ */
