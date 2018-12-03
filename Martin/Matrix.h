/*
 * Matrix.h
 *
 *  Created on: Nov 28, 2018
 *      Author: brass
 */

#ifndef MARTIN_MATRIX_H_
#define MARTIN_MATRIX_H_

//#include "LinearMap.h"

namespace LinAlg {

template<typename F>
class Matrix //: public LinearMap<F,std::vector<F>,std::vector<F> > {
{
public:
	Matrix();
	virtual ~Matrix();
	Matrix(uint32_t Mrows, uint32_t Ncols,bool RowMajor);
	Matrix(const Matrix<F> &A);
	uint32_t M = 0;
	uint32_t N = 0;
	bool RowMajor = false;
	F *Entries = 0;

	void Evaluate(const std::vector<F> &x, std::vector<F> &y);
	void AdjointEvaluate(const std::vector<F> &x, std::vector<F> &y);
	void BlockEvaluate(const std::vector<F> &x, std::vector<F> &y, const uint32_t blocksize);
	void BlockEvaluate(const Matrix<F> &x, Matrix<F> &y);
	F* &operator[](uint32_t i);
	void resize(uint32_t Mnew, uint32_t Nnew);
	void RemoveCol(uint32_t i);
	void print();


};

} /* namespace LinAlg */

#endif /* MARTIN_MATRIX_H_ */
