/*
 * BSpline.h
 *
 *  Created on: Sep 10, 2019
 *      Author: brass
 */

#ifndef MARTIN_MATH_BSPLINE_H_
#define MARTIN_MATH_BSPLINE_H_

/*
 * BSpline.cpp
 *
 *  Created on: Apr 24, 2019
 *      Author: brass
 */



namespace BSpline
{

double Evaluate(uint32_t i, uint32_t n, const vector<double> &x, double t);

double EvalDerivative(uint32_t i, uint32_t n, const vector<double> &x,  double t, uint32_t k);

double Derivative(uint32_t idx, double t, const vector<double> &x, uint32_t n, const vector<double> &c, uint32_t k);

double deBoor(uint32_t idx, double t, const vector<double> &x, uint32_t n, const vector<double> &c);

void GenerateGrid(const std::vector<double> &t,std::vector<double> &x, uint32_t n);

void GetKnotMatrix(const std::vector<double> &x, const uint32_t n, const std::vector<double> &t, std::vector<double> &B);

void CalculateCellIntegrals(const vector<double> &x, uint32_t ks,uint32_t k, vector<double> &OffDiagInts1, vector<double> &OffDiagInts2, vector<double> &DiagInts);

double Overlap(uint32_t i, uint32_t j, uint32_t ks, const vector<double> &x);

void OverlapMatrix(uint32_t ks, vector<double> x, vector<double> &S);

void Components(const vector<Interpolation::InterpolatingFunction*> &funs, uint32_t ks, const vector<double> &x, vector<double> &S, vector<double> &c);

} // end of namespace BSpline






#endif /* MARTIN_MATH_BSPLINE_H_ */
