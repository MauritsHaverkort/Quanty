/*
 * PlotSpectrumWithSelfEnergy.cpp
 *
 *  Created on: Oct 25, 2018
 *      Author: brass
 */

#include "CSRMatrix.h"

double ImExpectationValue(const CompactMatrixType &M, const std::vector<double> &v)
{
	double result = 0;
	for(uint32_t i = 0; i < v.size(); i++)
	{
		for(uint32_t j = 0; j < v.size(); j++)
		{
			result += v[i] * M.V[i*M.Ni+j].Im * v[j];
		}
	}
	return result;
}
double ImExpectationValue(const CompactMatrixType *M, const std::vector<double> &v)
{
	double result = 0;
	for(uint32_t i = 0; i < v.size(); i++)
	{
		for(uint32_t j = 0; j < v.size(); j++)
		{
			result += v[i] * M->V[i*M->Ni+j].Im * v[j];
		}
	}
	return result;
}

double EvaluateResponseFunction(double w, double gH, double gS,const std::vector<double> &H, const std::vector<ResponsefunctionType *> &SigmaRes, const std::vector<double> psi)
{
	CompactMatrixType sigma, temp;
	EvaluateBlockResponsefunctionType(SigmaRes[0],&sigma,w,gS);
	for (uint32_t i = 1; i < SigmaRes.size(); i++)
	{
		EvaluateBlockResponsefunctionType(SigmaRes[i],&temp,w,gS);
		CompactMatrixAdd(&sigma,&temp);
		FreeCompactMatrix(&temp);
	}
	uint32_t N = sigma.Nj;
	for (uint32_t i = 0; i < N; i++)
	{
		sigma.V[i*N +i].Re += H[i] - w;
		sigma.V[i*N +i].Im += -gH;
	}

	CompactMatrixInverse(&sigma);
	//printf("%e %e I\n",sigma.V[0].Re,sigma.V[0].Im);
	return ImExpectationValue(sigma,psi);
}

double EvaluateResponseFunction(double w, double gH, const std::vector<double> &H, CompactMatrixType * Sigma, const std::vector<double> psi)
{
	uint32_t N = Sigma->Nj;
	for (uint32_t i = 0; i < N; i++)
	{
		Sigma->V[i*N +i].Re += H[i] - w;
		Sigma->V[i*N +i].Im += -gH/2.0;
	}

	CompactMatrixInverse(Sigma);
	//printf("%e %e I\n",sigma.V[0].Re,sigma.V[0].Im);
	return ImExpectationValue(Sigma,psi);
}

CompactMatrixType *EvaluateSelfEnergy(double w, double gS, const std::vector<ResponsefunctionType *> &SigmaRes)
{
	CompactMatrixType *sigma, temp;
	sigma = (CompactMatrixType *)malloc(sizeof(CompactMatrixType));
	EvaluateBlockResponsefunctionType(SigmaRes[0],sigma,w,gS);
	for (uint32_t i = 1; i < SigmaRes.size(); i++)
	{
		EvaluateBlockResponsefunctionType(SigmaRes[i],&temp,w,gS);
		CompactMatrixAdd(sigma,&temp);
		FreeCompactMatrix(&temp);
	}

	return sigma;
}

CompactMatrixType *EvaluateSelfEnergy(double w, double gS, const std::vector<ResponsefunctionType *> &SigmaRes, const std::vector<double> &energyshift)
{
	CompactMatrixType *sigma, temp;
	sigma = (CompactMatrixType *)malloc(sizeof(CompactMatrixType));
	EvaluateBlockResponsefunctionType(SigmaRes[0],sigma,w-energyshift[0],gS);
	for (uint32_t i = 1; i < SigmaRes.size(); i++)
	{
		EvaluateBlockResponsefunctionType(SigmaRes[i],&temp,w-energyshift[i],gS);
		CompactMatrixAdd(sigma,&temp);
		FreeCompactMatrix(&temp);
	}

	return sigma;
}

void WriteDataToFile(CompactMatrixType * M, const char *filename)
{
	FILE *f = fopen(filename,"wb");

	fwrite(&(M->QComplex),sizeof(uint32_t),1,f);
	fwrite(&(M->Ni),sizeof(uint32_t),1,f);
	fwrite(&(M->Nj),sizeof(uint32_t),1,f);
	if (M->QComplex)
	{
		fwrite(M->V,sizeof(Complex16),M->Ni*M->Nj,f);
	}
	else
	{
		fwrite(M->VRe,sizeof(double),M->Ni*M->Nj,f);
	}
	fclose(f);

}

template<typename F>
void WriteDataToFile(LinAlg::CSRMatrix<F> *M, const char *filename)
{
	FILE *f = fopen(filename,"wb");

	fwrite(&(M->Ncol), sizeof(uint32_t), 1, f);
	size_t N = M->A.size();
	fwrite(&N, sizeof(size_t), 1, f);
	fwrite(&(M->A[0]), sizeof(F), N, f);
	N = M->AI.size();
	fwrite(&N, sizeof(size_t), 1, f);
	fwrite(&(M->AI[0]), sizeof(uint32_t), N, f);
	N = M->AJ.size();
	fwrite(&N, sizeof(size_t), 1, f);
	fwrite(&(M->AJ[0]), sizeof(uint32_t), N, f);
	fclose(f);

}
// M has to be allocated at input as M = (LinAlg::CSRMatrix<F> *) malloc(sizeof(LinAlg::CSRMatrix<F>));
template<typename F>
void ReadDataFromFile(LinAlg::CSRMatrix<F> *M, const char *filename)
{
	if (!M)
	{
		printf("Error in ReadDataFromFile: CSRMatrix has to be allocated at input\n");
		fflush(stdout);
		return;
	}

	FILE *f = fopen(filename,"rb");

	M = new(M) LinAlg::CSRMatrix<F> ;

	fread(&(M->Ncol), sizeof(uint32_t), 1, f);
	size_t N = 0;
	fread(&N, sizeof(size_t), 1, f);
	M->A.reserve(N);
	F value;
	for (size_t i = 0; i < N; i++)
	{
		fread(&value, sizeof(F), 1, f);
		M->A.push_back(value);
	}
	fread(&N, sizeof(size_t), 1, f);
	M->AI.reserve(N);
	uint32_t idx;
	for (size_t i = 0; i < N; i++)
	{
		fread(&idx, sizeof(uint32_t), 1, f);
		M->AI.push_back(idx);
	}
	fread(&N, sizeof(size_t), 1, f);
	M->AJ.reserve(N);
	for (size_t i = 0; i < N; i++)
	{
		fread(&idx, sizeof(uint32_t), 1, f);
		M->AJ.push_back(idx);
	}
	fclose(f);

}

void ReadDataFromFile(CompactMatrixType * M, const char *filename)
{
	FILE *f = fopen(filename,"rb");

	M = (CompactMatrixType *) malloc(sizeof(CompactMatrixType));
	fread(&(M->QComplex),sizeof(uint32_t),1,f);
	fread(&(M->Ni),sizeof(uint32_t),1,f);
	fread(&(M->Nj),sizeof(uint32_t),1,f);
	if (M->QComplex)
	{
		M->V = (Complex16 *) malloc(M->Ni*M->Nj*sizeof(Complex16));
		fread(M->V,sizeof(Complex16),M->Ni*M->Nj,f);
	}
	else
	{
		M->VRe = (double *) malloc(M->Ni*M->Nj*sizeof(double));
		fread(M->VRe,sizeof(double),M->Ni*M->Nj,f);
	}
	fclose(f);
}

void WriteDataToFile(std::vector<CompactMatrixType *> Ms, const char *filename)
{
	FILE *f = fopen(filename,"wb");
	size_t N = Ms.size();
	fwrite(&N,sizeof(size_t),1,f);
	for (uint32_t i = 0; i < N; i++)
	{
		fwrite(&(Ms[i]->QComplex),sizeof(uint32_t),1,f);
		fwrite(&(Ms[i]->Ni),sizeof(uint32_t),1,f);
		fwrite(&(Ms[i]->Nj),sizeof(uint32_t),1,f);
		if (Ms[i]->QComplex)
		{
			fwrite(Ms[i]->V,sizeof(Complex16),Ms[i]->Ni*Ms[i]->Nj,f);
		}
		else
		{
			fwrite(Ms[i]->VRe,sizeof(double),Ms[i]->Ni*Ms[i]->Nj,f);
		}
	}
	fclose(f);
}

void ReadDataFromFile(std::vector<CompactMatrixType *> &Ms, const char *filename)
{
	FILE *f = fopen(filename,"rb");

	size_t N;
	fread(&N,sizeof(size_t),1,f);
	if (Ms.size()!=0)
	{
		printf("Error: ReadDataFromFile: vector were data are stored is not empty\n");
	}
	else
	{
		Ms.reserve(N);

		for (uint32_t i = 0; i < N; i++)
		{
			Ms.push_back((CompactMatrixType *) malloc(sizeof(CompactMatrixType)));
			fread(&(Ms[i]->QComplex),sizeof(uint32_t),1,f);
			fread(&(Ms[i]->Ni),sizeof(uint32_t),1,f);
			fread(&(Ms[i]->Nj),sizeof(uint32_t),1,f);
			if (Ms[i]->QComplex)
			{
				Ms[i]->V = (Complex16 *) malloc(Ms[i]->Ni*Ms[i]->Nj*sizeof(Complex16));
				fread(Ms[i]->V,sizeof(Complex16),Ms[i]->Ni*Ms[i]->Nj,f);
			}
			else
			{
				Ms[i]->VRe = (double *) malloc(Ms[i]->Ni*Ms[i]->Nj*sizeof(double));
				fread(Ms[i]->VRe,sizeof(double),Ms[i]->Ni*Ms[i]->Nj,f);
			}
		}
	}
	fclose(f);
}
