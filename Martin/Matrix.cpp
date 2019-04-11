/*
 * Matrix.cpp
 *
 *  Created on: Nov 28, 2018
 *      Author: brass
 */

#include "Matrix.h"

namespace LinAlg {
template<typename F>
Matrix<F>::Matrix() {
	//printf("let there be math\n");
}
template<typename F>
Matrix<F>::~Matrix() {
	free(Entries);
	//printf("Ceterum censeo Carthaginem esse delendam.\n");
}
template<typename F>
Matrix<F>::Matrix(uint32_t Mrows,uint32_t Ncols,bool RowMajor)
{
	M = Mrows;
	N = Ncols;
	this->RowMajor = RowMajor;
	Entries = (F*) malloc( sizeof(F)*M*N );
	if(!Entries)
	{
		printf("Failed to allocate memory for %u by %u Matrix\n",M,N);
	}
	//printf("alloc constr\n");

}

template<typename F>
void Matrix<F>::copy(const Matrix<F> &A)
{
	M = A.M;
	N = A.N;
	RowMajor = A.RowMajor;
	free(Entries);
	Entries = (F*) malloc( sizeof(F)*M*N );
	if(!Entries)
	{
		printf("Failed to allocate memory for copying a %u by %u Matrix\n",M,N);
	}
	for (uint32_t i = 0; i < M*N; i++)
	{
		Entries[i] = A.Entries[i];
	}
	//printf("zu Guttenberg-Function\n");
}
template<typename F>
void Matrix<F>::Evaluate(const std::vector<F> &x, std::vector<F> &y)
{

}
template<typename F>
void Matrix<F>::AdjointEvaluate(const std::vector<F> &x, std::vector<F> &y)
{

}
template<typename F>
void Matrix<F>::BlockEvaluate(const std::vector<F> &x, std::vector<F> &y, const uint32_t blocksize)
{

}
template<typename F>
void Matrix<F>::BlockEvaluate(const Matrix<F> &x, Matrix<F> &y)
{

}
template<typename F>
void Matrix<F>::print()
{
	if (RowMajor)
	{
		for (uint32_t i = 0; i < M; i++)
		{
			for (uint32_t j = 0; j < N; j++)
			{
				printf("%.4f\t",Entries[i*N+j]);
			}
			printf("\n");
		}
	}
	else
	{
		for (uint32_t i = 0; i < M; i++)
		{
			for (uint32_t j = 0; j < N; j++)
			{
				printf("%.16E\t",Entries[i+j*M]);
			}
			printf("\n");
		}
	}
}

template<typename F>
void Matrix<F>::RemoveCol(uint32_t i)
{
	if (RowMajor)
	{
		printf("Removing column of RowMajor Matrix? I'm too lazy for that\n");
	}
	else
	{
		memmove(&(Entries[i*M]),&(Entries[(i+1)*M]),(N-i-1)*sizeof(F)*M);
		N--;
		Entries = (F*) realloc(Entries,M*N*sizeof(F));

	}
}

template<typename F>
F* &Matrix<F>::operator[](uint32_t i)
{
	if (RowMajor)
	{
		if (i < M)
		{
			return &(Entries[i*M]);
		}
		else
		{
			printf("Error: Matrix has %u rows, but row number %d was requested\n",M,i);
			return 0;
		}
	}
	else
	{
		printf("Error: Matrix is stored in column major format, hence Matrix[i][j] is not yet defined\n");
		return 0;
	}
}
template<typename F>
void Matrix<F>::reserve(uint32_t Mnew, uint32_t Nnew)
{
	M = Mnew;
	N = Nnew;
	if(Entries)
	{
		Entries = (F*) realloc(Entries,M*N*sizeof(F));
	}
	else
	{
		Entries = (F*) malloc (M*N*sizeof(F));

	}
	if(!Entries)
	{
		printf("Failed to reserve memory for %u by %u Matrix\n",M,N);
	}
}
template<typename F>
void Matrix<F>::resize(uint32_t Mnew, uint32_t Nnew)
{
	M = Mnew;
	N = Nnew;
	Entries = (F*) realloc(Entries,M*N*sizeof(F));
}

} /* namespace LinAlg */
