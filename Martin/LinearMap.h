/*
 * LinearMap.h
 *
 *  Created on: Nov 20, 2018
 *      Author: brass
 */

#ifndef MARTIN_LINEARMAP_H_
#define MARTIN_LINEARMAP_H_

#include <stdio.h>

namespace LinAlg {

/*
 * let F denote a field, V and W denote vectorspaces defined over the field f
 * LinearMap M: V -> W
 * 			    x -> y = M(x) by a call to M.Evaluate(x,y)
 *
 * typical types for F are double or Complex
 * typical types for V and W are std::vector<F> or WaveFunctionType
 * 		   hence V and W denote the type of the vectors that can be mapped
 *
 * the member function M.AdjointEvaluate(x,y) calculates y = adj(M)(x)
 * 			where adj(M) is the canonical hermitian adjoint of M
*/
template<typename F, typename V , typename W>
class LinearMap {
public:
	virtual void Evaluate(const V &x, W &y) = 0;
	virtual void AdjointEvaluate(const W &x, V &y) = 0;
};

} /* namespace LinAlg */

#endif /* MARTIN_LINEARMAP_H_ */
