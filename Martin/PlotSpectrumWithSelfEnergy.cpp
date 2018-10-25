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


