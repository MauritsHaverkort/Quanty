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

public:

	std::vector<double> y2;

	CubicSpline();
	CubicSpline(const std::vector<double> &u, const std::vector<double> &v);
	virtual ~CubicSpline();

	void Create(const std::vector<double> &u, const std::vector<double> &v);
    void Create(const std::vector<double> &u, const std::vector<double> &v, double yp1, double ypn);
	std::vector<double> GetY2();
	double GetY2(uint32_t i);

	double operator()(double t);
	double operator()(double t, unsigned int k);
	std::vector<double> operator()(const std::vector<double> &t);
	std::vector<double> operator()(const std::vector<double> &t, unsigned int k);
};

bool operator==(CubicSpline &f,CubicSpline &g);
void UnaryMinus(CubicSpline &f,CubicSpline &g);

} /* namespace Interpolation */

void ThomasAlgorithm(const std::vector<double> &a,const std::vector<double> &b,const std::vector<double> &c,const std::vector<double> &d, std::vector<double> &x);

#endif /* CUBICSPLINE_H_ */
