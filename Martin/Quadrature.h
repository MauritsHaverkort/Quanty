/*
 * Quadrature.h
 *
 *  Created on: Jun 28, 2018
 *      Author: brass
 */

#ifndef QUADRATURE_H_
#define QUADRATURE_H_

#include <vector>
using std::vector;

#define NKnotsMax 7

class Quadrature
{
private:
	vector<vector<double> > weights;

	vector<vector<double> > knots;

public:
	Quadrature();
	virtual ~Quadrature();

	vector<double> GetKnots(int n);
	vector<double> GetKnots(int n, double a, double b);
	double Integrate(vector<double> f);
	double Integrate(vector<double> f,double a, double b);

};

#endif /* QUADRATURE_H_ */
