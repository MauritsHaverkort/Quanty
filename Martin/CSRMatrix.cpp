/*
 * CSRMatrix.cpp
 *
 *  Created on: Nov 19, 2018
 *      Author: brass
 */

#include "CSRMatrix.h"

namespace LinAlg {

template<typename F>
CSRMatrix<F>::CSRMatrix() {
	// TODO Auto-generated constructor stub

}

template<typename F>
CSRMatrix<F>::~CSRMatrix() {
	// TODO Auto-generated destructor stub
}

template<typename F>
int qsortcompareSparseTriple(const void *a, const void *b)
{
	SparseTriple_t<F> *A = (SparseTriple_t<F>*) a;
	SparseTriple_t<F> *B = (SparseTriple_t<F>*) b;
	if (A->row == B->row)
	{
		return A->col - B->col;
	}
	else
	{
		return A->row - B->row;
	}
}

template<typename F>
CSRMatrix<F>::CSRMatrix(std::vector<SparseTriple_t<F> > &M)
{
	uint32_t N = M.size();
	qsort(&M[0],N,sizeof(SparseTriple_t<F>),qsortcompareSparseTriple<F>);
	A.reserve(N);
	AJ.reserve(N);
	AI.push_back(0);
	AI.push_back(0);
	uint32_t currRow = 0;
	for (uint32_t i = 0 ; i < N; i++)
	{
		A.push_back(M[i].val);
		AJ.push_back(M[i].col);
		Ncol = std::max(Ncol,M[i].col+1);
		if (M[i].row > currRow)
		{
			for (uint32_t j = currRow; j < M[i].row; j++)
			{
				AI.push_back(AI[currRow+1]);
			}
			currRow = M[i].row;
		}
		AI[currRow+1]++;
	}
}
template<typename F>
void CSRMatrix<F>::AddRow(std::map<uint32_t,F> rowEntries)
{
	uint32_t M = AI.size();
	if (M == 0)
	{
		AI.push_back(0);
		M++;
	}
	AI.push_back(AI[M-1]);
	for(typename std::map<uint32_t, F>::iterator it = rowEntries.begin(); it != rowEntries.end(); ++it)
	{
		AJ.push_back(it->first);
		Ncol = std::max(Ncol,it->first+1);
		A.push_back(it->second);
		AI[M]++;
	}
}


template<typename F>
void CSRMatrix<F>::Evaluate(const std::vector<F> &x,std::vector<F> &y)
{
	uint32_t N = x.size();
	uint32_t M = AI.size()-1;

	y.reserve(M);

	for (uint32_t i = 0; i < M; i++)
	{
		y.push_back(0);
		for (uint32_t j = AI[i]; j < AI[i+1]; j++)
		{
			y[i] += A[j] * x[AJ[j]];
		}
	}
}
template<typename F>
void CSRMatrix<F>::AdjointEvaluate(const std::vector<F> &x,std::vector<F> &y)
{
	uint32_t M = AI.size()-1;

	y = std::vector<F>(Ncol);

	for (uint32_t i = 0; i < M; i++)
	{
		for (uint32_t j = AI[i]; j < AI[i+1]; j++)
		{
			y[AJ[j]] += A[j] * x[i];
		}
	}
}

template<typename F>
void CSRMatrix<F>::Print()
{
	printf("A \t");
	for (uint32_t i = 0; i < A.size();i++)
	{
		printf("%.1f\t",A[i]);
	}
	printf("\n");
	printf("AJ\t");
	for (uint32_t i = 0; i < A.size();i++)
	{
		printf("%d\t",AJ[i]);
	}
	printf("\n");
	printf("AI\t");
	for (uint32_t i = 0; i < AI.size();i++)
	{
		printf("%d\t",AI[i]);
	}
	printf("\n");
}

} /* namespace LinAlg */
