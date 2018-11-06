/*
 * PlotSpectrumWithSelfEnergy.cpp
 *
 *  Created on: Oct 25, 2018
 *      Author: brass
 */

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

CompactMatrixType *EvaluateSelfEnergy(double w, double gS, const std::vector<ResponsefunctionType *> &SigmaRes)
{
	CompactMatrixType *sigma, temp;
	EvaluateBlockResponsefunctionType(SigmaRes[0],sigma,w,gS);
	for (uint32_t i = 1; i < SigmaRes.size(); i++)
	{
		EvaluateBlockResponsefunctionType(SigmaRes[i],&temp,w,gS);
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
		fread(M->V,sizeof(Complex16),M->Ni*M->Nj,f);
	}
	else
	{
		fread(M->VRe,sizeof(double),M->Ni*M->Nj,f);
	}
}
