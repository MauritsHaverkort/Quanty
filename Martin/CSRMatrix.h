/*
 * CSRMatrix.h
 *
 *  Created on: Nov 19, 2018
 *      Author: brass
 */

#ifndef MARTIN_CSRMATRIX_H_
#define MARTIN_CSRMATRIX_H_

#include <vector>
#include <map>
#include "LinearMap.h"

namespace LinAlg {

template <typename F>
struct SparseTriple_t{
	F val;
	uint32_t row, col;
};

template <typename F> // F stands for Field upon which the matrix is defined (double or Complex)
class CSRMatrix: public LinearMap<F,std::vector<F>,std::vector<F> > {
public:
	CSRMatrix();
	CSRMatrix(std::vector<SparseTriple_t<F> > &M);
	virtual ~CSRMatrix();

	std::vector<F> A;
	std::vector<uint32_t> AI, AJ;
	uint32_t Ncol = 0;

	void Print();
	void AddRow(std::map<uint32_t,F> rowEntries);
	void Evaluate(const std::vector<F> &x,std::vector<F> &y);
	void AdjointEvaluate(const std::vector<F> &x,std::vector<F> &y);
	void BlockEvaluate(const std::vector<F> &x,std::vector<F> &y, const uint32_t blocksize);

};

} /* namespace LinAlg */

#endif /* MARTIN_CSRMATRIX_H_ */
