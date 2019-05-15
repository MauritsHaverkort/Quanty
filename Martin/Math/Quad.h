/*
 * Quad.h
 *
 *  Created on: May 2, 2019
 *      Author: brass
 */

#ifndef MARTIN_MATH_QUAD_H_
#define MARTIN_MATH_QUAD_H_

void GaussianQuadratureWeights(uint32_t N, vector<double> &x, vector<double> &w, double a=-1, double b=1);

double LegendreP(double x, uint32_t n);

#endif /* MARTIN_MATH_QUAD_H_ */
