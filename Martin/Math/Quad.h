/*
 * Quad.h
 *
 *  Created on: May 2, 2019
 *      Author: brass
 */

#ifndef MARTIN_MATH_QUAD_H_
#define MARTIN_MATH_QUAD_H_

namespace Quad
{

void GaussianQuadratureWeights(uint32_t N, vector<double> &x, vector<double> &w, double a=-1, double b=1);

void GaussianQuadratureWeights(uint32_t N, const vector<double> &x0, vector<double> &x, vector<double> &w);

double LegendreP(double x, uint32_t n);

template <typename T>
double Integrate(T &f, double a, double b, uint32_t N);

}
#endif /* MARTIN_MATH_QUAD_H_ */
