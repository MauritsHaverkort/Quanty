/*
 * Continuum-Diagonal-Method.cpp
 *
 *  Created on: Mar 12, 2019
 *      Author: brass
 */

typedef struct mypoint
{
	double w;
	Complex16 y;
}mypoint;

void PsiOppListOppListPsi(WaveFunctionType *psi, vector<OperatorType*> U, CompactMatrixType *M)
{
	OperatorType UU;
	double Re, Im;
	for (uint32_t i = 0; i < M->Ni; i++)
	{
		for (uint32_t j = i; j < M->Nj; j++)
		{
			InitOperatorForMult(&UU,2,U[i],U[j]);
			OperatorMultiply2NormalOrder(&UU, *U[j], *U[i], 'D', 'N', &QDetNotOrderedFcFaBcBa);
			PsiOperatorPsi(UU,*psi,&Re,&Im);
			M->VRe[i*M->Nj + j] = Re;
			M->VRe[j*M->Nj + i] = Re;
			FreeOperator(&UU);
		}
	}
}

void PsiOppListOperatorOppListPsi(WaveFunctionType *psi, vector<OperatorType*> U, OperatorType *H, CompactMatrixType *M)
{
	OperatorType UHU, HU;
	double Re, Im;
	if(M->QComplex)
	{
		for (uint32_t i = 0; i < M->Ni; i++)
		{
			InitOperatorForMult(&HU,2,U[i],H);
			OperatorMultiply2NormalOrder(&HU,*H,*U[i],'N','N',&QDetNotOrderedFcFaBcBa);
			for (uint32_t j = i; j < M->Nj; j++)
			{
				//InitOperatorForMult(&UHU,3,U[i],U[j],H);
				//OperatorMultiply3NormalOrder(&UHU, *U[j], *H, *U[i], 'D', 'N', 'N', &QDetNotOrderedFcFaBcBa);
				InitOperatorForMult(&UHU,2,U[j],HU);
				OperatorMultiply2NormalOrder(&UHU,*U[j],HU,'D','N',&QDetNotOrderedFcFaBcBa);
				PsiOperatorPsi(UHU,*psi,&Re,&Im);
				M->V[i*M->Nj + j].Re = Re;
				M->V[i*M->Nj + j].Im = Im;
				M->V[j*M->Nj + i] = M->V[i*M->Nj + j];
				FreeOperator(&UHU);
			}
			FreeOperator(&HU);
		}
	}
	else
	{
#pragma omp parallel for
		for (uint32_t i = 0; i < M->Ni; i++)
		{
			OperatorType hu;
			InitOperatorForMult(&hu,2,H,U[i]);
			OperatorMultiply2NormalOrder(&hu,*H,*U[i],'N','N',&QDetNotOrderedFcFaBcBa);
			for (uint32_t j = i; j < M->Nj; j++)
			{
				//InitOperatorForMult(&UHU,3,U[i],U[j],H);
				//OperatorMultiply3NormalOrder(&UHU, *U[j], *H, *U[i], 'D', 'N', 'N', &QDetNotOrderedFcFaBcBa);
				OperatorType uhu;
				InitOperatorForMult(&uhu,2,U[j],&hu);
				OperatorMultiply2NormalOrder(&uhu,*U[j],hu,'D','N',&QDetNotOrderedFcFaBcBa);
				double re, im;
				PsiOperatorPsi(uhu,*psi,&re,&im);
				M->VRe[i*M->Nj + j] = re;
				M->VRe[j*M->Nj + i] = re;
				FreeOperator(&uhu);
			}
			FreeOperator(&hu);
		}
	}

}

uint32_t GetNumberOfPolesBelowCutoff(const std::vector<double> &poles, double cutoff)
{
	uint32_t i = 0;
	while (i < (uint32_t) poles.size() && poles[i] < cutoff)
	{
		i++;
	}
	return i;
}

void CalculateResidues(const std::vector<double> &v, CompactMatrixType *M, uint32_t Nresidues, std::vector<double> &residues)
{
	if(M->Ni < Nresidues || M->Nj < (uint32_t) v.size())
	{
		printf("error in CalculateResidues: CompactMatrix has wrong dimensions: Ni = %d Nj = %d\n",M->Ni,M->Nj);
		fflush(stdout);
		return;
	}
	if(M->QComplex)
	{
		printf("error in CalculateResidues: CompactMatrix is complex. I should implement that\n");
		fflush(stdout);
		return;
	}
	residues = std::vector<double>(Nresidues);
	for (uint32_t i = 0; i < Nresidues; i++)
	{
		residues[i] = 0;
		for (uint32_t j = 0; j < (uint32_t) v.size(); j++)
		{
			residues[i] += M->VRe[i*M->Nj + j] * v[j];
		}
		residues[i] *= residues[i];
	}
}

void SelectDominantResiduesAndPoles(const std::vector<double> &poles, const std::vector<double> &residues, double eps, std::vector<double> &selectedPoles, std::vector<double> &selectedResidues)
{
	if (poles.size() < residues.size())
	{
		printf("error in SelectDominantResiduesAndPoles: poles.size() = %d < %d = residues.size()\n",poles.size(),residues.size());
		fflush(stdout);
		return;
	}
	for (uint32_t i = 0; i < (uint32_t) residues.size(); i++)
	{
		if (residues[i] > eps)
		{
			selectedPoles.push_back(poles[i]);
			selectedResidues.push_back(residues[i]);
		}
	}
}

Complex16 Resolvent(Complex16 z, const std::vector<double> &poles, const std::vector<double> &residues)
{
	Complex16 result;
	result.Re = 0;
	result.Im = 0;
	if (poles.size() < residues.size())
	{
		printf("error in Resolvent: poles.size() = %d < %d = residues.size()",poles.size(),residues.size());
		fflush(stdout);
		return result;
	}
	for (uint32_t i = 0; i < (uint32_t) residues.size(); i++)
	{
		result += residues[i]/(z - poles[i]);
	}
	return result;
}

Complex16 Resolvent(Complex16 z, const std::vector<double> &poles, const std::vector<double> &residues,const vector<Complex16> &selfenergy)
{
	Complex16 result;
	result.Re = 0;
	result.Im = 0;
	if (poles.size() < residues.size() || poles.size() != selfenergy.size())
	{
		printf("error in Resolvent: poles.size() = %d < %d = residues.size()",poles.size(),residues.size());
		fflush(stdout);
		return result;
	}
	for (uint32_t i = 0; i < (uint32_t) residues.size(); i++)
	{
		result += residues[i]/(z - poles[i] - selfenergy[i]);
	}
	return result;
}

void PlotSpectrumFromPoles(double eps,double gamma, double shift, const std::vector<double> &poles, const std::vector<double> &v, CompactMatrixType *M, std::vector<mypoint> &spectrum)
{

	double wMax = spectrum[spectrum.size()-1].w;
	uint32_t Npoles = GetNumberOfPolesBelowCutoff(poles, wMax);
	std::vector<double> residues, selectedPoles, selectedResidues;
	CalculateResidues(v,M,Npoles,residues);
	SelectDominantResiduesAndPoles(poles,residues,eps,selectedPoles,selectedResidues);
#pragma omp parallel for
	for (uint32_t i = 0; i < (uint32_t) spectrum.size(); i++)
	{
		Complex16 z;
		z.Re = spectrum[i].w - shift;
		z.Im = gamma/2;
		spectrum[i].y = Resolvent(z,selectedPoles,selectedResidues);
	}

}

void PlotSpectrumFromPoles(double eps,double gamma, double shift, const std::vector<double> &poles, const std::vector<double> &v, CompactMatrixType *M, std::vector<mypoint> &spectrum,const vector<vector<Complex16> > &selfenergy)
{

	double wMax = spectrum[spectrum.size()-1].w;
	uint32_t Npoles = GetNumberOfPolesBelowCutoff(poles, wMax);
	std::vector<double> residues, selectedPoles, selectedResidues;
	CalculateResidues(v,M,Npoles,residues);
	SelectDominantResiduesAndPoles(poles,residues,eps,selectedPoles,selectedResidues);

#pragma omp parallel for
	for (uint32_t i = 0; i < (uint32_t) spectrum.size(); i++)
	{
		Complex16 z;
		z.Re = spectrum[i].w - shift;
		z.Im = gamma/2;
		spectrum[i].y = Resolvent(z,selectedPoles,selectedResidues,selfenergy[i]);
	}

}

void lua_push(lua_State * L, const mypoint p)
{
	lua_createtable(L,0,0);
	lua_push(L,p.w);
	lua_rawseti(L,-2,1);
	Complex16 *z = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
	luaL_setmetatable(L, "Complex_Type");
	*z = p.y;
	lua_rawseti(L,-2,2);
}
