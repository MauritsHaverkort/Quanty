/*
 * Modules.cpp
 *
 *  Created on: May 14, 2019
 *      Author: brass
 */

#include "Modules.h"

void SlaterIntegralsToList(const SSlaterIntegral *integrals, const uint32_t N, vector<double> &Fk)
{
	for (uint32_t i = 0; i < N; i++)
	{
		for (uint32_t j = 0; j < integrals[i].Nval; j++)
		{
			Fk.push_back(integrals[i].val[j]);
		}
	}
}
void SlaterIntegralsToList(const SSlaterIntegral *integrals, const uint32_t N, vector<vector<double> >&Fk)
{
	Fk = vector<vector<double> >(N);
	for (uint32_t i = 0; i < N; i++)
	{
		Fk[i] = vector<double>(integrals[i].Nval);
		for (uint32_t j = 0; j < integrals[i].Nval; j++)
		{
			Fk[i][j]=integrals[i].val[j];
		}
	}
}

void CreateContinuumSlaterIntegrals(vector<Interpolation::InterpolatingFunction *> &WavefunsBoundG, vector<Interpolation::InterpolatingFunction *> &WavefunsBoundF,
							vector<Interpolation::InterpolatingFunction *> &WavefunsUnboundG, vector<Interpolation::InterpolatingFunction *> &WavefunsUnboundF,
							vector<int> &indToK, const int kappa, vector<vector<double> > &Fks)
{
	SSlaterIntegral *si;
	uint32_t kMax, n, NIntegrals, NBound, NUnbound;
	double **kernels;
	vector<vector<double> > dualfunsBound, dualfunsUnbound;

	vector<double> x = WavefunsBoundG[0]->GetX();
	vector<vector<double> > WavefunsBound = PutRelativisticSplineOnGrid(WavefunsBoundG,WavefunsBoundF,x,0);
	vector<vector<double> > WavefunsUnbound = PutRelativisticSplineOnGrid(WavefunsUnboundG,WavefunsUnboundF,x,0);
	NUnbound = (uint32_t) WavefunsUnbound.size();
	NBound = (uint32_t) WavefunsBound.size();

	unsigned int NGridPoints = (unsigned int)x.size();
	n = 2* NGridPoints;
	CreateDualFunctions(WavefunsBound,x,dualfunsBound);
	CreateDualFunctions(WavefunsUnbound,x,dualfunsUnbound);
	indToK.push_back(kappa);
	si = InitializeIntegrals(indToK, NIntegrals,&kMax,true);
	kernels = CreateAllKernels(kMax, NGridPoints, &x[0]);
	WavefunsBound.push_back(WavefunsUnbound[0]);
	dualfunsBound.push_back(dualfunsUnbound[0]);

	Fks = vector<vector<double> >(NUnbound);
	for (uint32_t i = 0; i < NUnbound; i++)
	{
		WavefunsBound[NBound] = WavefunsUnbound[i];
		dualfunsBound[NBound] = dualfunsUnbound[i];
		CalculateSlaterIntegrals(NIntegrals, si,n,kernels,WavefunsBound,dualfunsBound);
		SlaterIntegralsToList(si,NIntegrals,Fks[i]);
	}
	indToK.pop_back();
//	WavefunsBoundG.pop_back();	 //TODO: delte if all test done - it is nonsense
//	WavefunsBoundF.pop_back();

	for (uint32_t i = 0; i <= kMax; i++)
	{
		free(kernels[i]);
	}
	free(kernels);
	for (uint32_t i = 0; i < NIntegrals; i++)
	{
		free(si[i].val);
		free(si[i].k);
	}
	free(si);

}

void CreateContinuumSlaterIntegrals(vector<Interpolation::InterpolatingFunction *> &WavefunsBoundG, vector<Interpolation::InterpolatingFunction *> &WavefunsBoundF,
							vector<Interpolation::InterpolatingFunction *> &WavefunsUnboundG, vector<Interpolation::InterpolatingFunction *> &WavefunsUnboundF,
							vector<int> &indToK, const int kappa, vector<vector<vector<double> > > &Fks)
{
	SSlaterIntegral *si;
	uint32_t kMax, n, NIntegrals, NBound, NUnbound;
	double **kernels;
	vector<vector<double> > dualfunsBound, dualfunsUnbound;

	vector<double> x = WavefunsBoundG[0]->GetX();
	vector<vector<double> > WavefunsBound = PutRelativisticSplineOnGrid(WavefunsBoundG,WavefunsBoundF,x,0);
	vector<vector<double> > WavefunsUnbound = PutRelativisticSplineOnGrid(WavefunsUnboundG,WavefunsUnboundF,x,0);
	NUnbound = (uint32_t) WavefunsUnbound.size();
	NBound = (uint32_t) WavefunsBound.size();

	unsigned int NGridPoints = (unsigned int)x.size();
	n = 2* NGridPoints;
	CreateDualFunctions(WavefunsBound,x,dualfunsBound);
	CreateDualFunctions(WavefunsUnbound,x,dualfunsUnbound);
	indToK.push_back(kappa);
	si = InitializeIntegrals(indToK, NIntegrals,&kMax,true);
	kernels = CreateAllKernels(kMax, NGridPoints, &x[0]);
	WavefunsBound.push_back(WavefunsUnbound[0]);
	dualfunsBound.push_back(dualfunsUnbound[0]);

	Fks = vector<vector<vector<double> > >(NUnbound);
	for (uint32_t i = 0; i < NUnbound; i++)
	{
		WavefunsBound[NBound] = WavefunsUnbound[i];
		dualfunsBound[NBound] = dualfunsUnbound[i];
		CalculateSlaterIntegrals(NIntegrals, si,n,kernels,WavefunsBound,dualfunsBound);
		SlaterIntegralsToList(si,NIntegrals,Fks[i]);
	}
	indToK.pop_back();
//	WavefunsBoundG.pop_back(); //TODO: delte if all test done - it is nonsense
//	WavefunsBoundF.pop_back();

	for (uint32_t i = 0; i <= kMax; i++)
	{
		free(kernels[i]);
	}
	free(kernels);
	for (uint32_t i = 0; i < NIntegrals; i++)
	{
		free(si[i].val);
		free(si[i].k);
	}
	free(si);

}

void CreateContinuumSlaterIntegrals(vector<Interpolation::InterpolatingFunction *> &WavefunsBoundG,	vector<Interpolation::InterpolatingFunction *> &WavefunsUnboundG,
									vector<int> &indToK, const int kappa, vector<vector<vector<double> > > &Fks)
{
	SSlaterIntegral *si;
	uint32_t kMax, NIntegrals, NBound, NUnbound;
	double **kernels;
	vector<vector<double> > dualfunsBound, dualfunsUnbound;

	NUnbound = (uint32_t) WavefunsUnboundG.size();
	NBound = (uint32_t) WavefunsBoundG.size();

	vector<double> x = WavefunsBoundG[0]->GetX();
	vector<vector<double> > WavefunsBound = vector<vector<double> >(NBound);
	for (uint32_t i = 0; i < NBound; i++)
	{
		WavefunsBound[i] = (*WavefunsBoundG[i])(x);
	}
	vector<vector<double> > WavefunsUnbound = vector<vector<double> >(NUnbound);
	for (uint32_t i = 0; i < NUnbound; i++)
	{
		WavefunsUnbound[i] = (*WavefunsUnboundG[i])(x);
	}


	unsigned int NGridPoints = (unsigned int)x.size();
	CreateDualFunctions(WavefunsBound,x,dualfunsBound,false);
	CreateDualFunctions(WavefunsUnbound,x,dualfunsUnbound,false);
	indToK.push_back(kappa);								//TODO: make this consistent
	si = InitializeIntegrals(indToK, NIntegrals,&kMax,true);
	kernels = CreateAllKernels(kMax, NGridPoints, &x[0]);
	WavefunsBound.push_back(WavefunsUnbound[0]);
	dualfunsBound.push_back(dualfunsUnbound[0]);

	Fks = vector<vector<vector<double> > >(NUnbound);
	for (uint32_t i = 0; i < NUnbound; i++)
	{
		WavefunsBound[NBound] = WavefunsUnbound[i];
		dualfunsBound[NBound] = dualfunsUnbound[i];
		CalculateSlaterIntegrals(NIntegrals, si,NGridPoints,kernels,WavefunsBound,dualfunsBound);
		SlaterIntegralsToList(si,NIntegrals,Fks[i]);
	}
	indToK.pop_back();

	for (uint32_t i = 0; i <= kMax; i++)
	{
		free(kernels[i]);
	}
	free(kernels);
	for (uint32_t i = 0; i < NIntegrals; i++)
	{
		free(si[i].val);
		free(si[i].k);
	}
	free(si);

}

void CreateAugerMeitnerU(vector<int> &indToK, uint32_t NF, uint16_t **indices, vector<OperatorType> &UList, bool QRelativistic)
{
	uint32_t kMax, NIntegrals;
	SSlaterIntegral *integrals = InitializeIntegrals(indToK, NIntegrals,&kMax,true);

	for(uint32_t i = 0; i < NIntegrals; i++)
	{
		int *k = integrals[i].tau.kappa;
		uint32_t *n = integrals[i].tau.n;
		double *Fk = integrals[i].val;
		for (uint32_t j = 0; j < integrals[i].Nval; j++)
		{
			Fk[j] = 0;
		}
		for (uint32_t j = 0; j < integrals[i].Nval; j++)
		{
			OperatorType optemp;

			Fk[j] = 1;
			if (QRelativistic)
			{
				CreateOperatorAtomicURelativisticOneIntegral(NF,k[0],k[1],k[2],k[3],indices[n[0]],indices[n[1]],indices[n[2]],indices[n[3]],Fk,n[0],n[1],n[2],n[3],&optemp);
			}
			else
			{
			    int l[4];
			    l[0] = -k[0]-1;
			    l[1] = -k[1]-1;
			    l[2] = -k[2]-1;
			    l[3] = -k[3]-1;
				CreateOperatorAtomicUNonRelOneIntegral(NF,l[0],l[1],l[2],l[3],indices[2*n[0]],indices[2*n[0]+1],indices[2*n[1]],indices[2*n[1]+1],indices[2*n[2]],indices[2*n[2]+1],indices[2*n[3]],indices[2*n[3]+1],Fk,n[0],n[1],n[2],n[3],&optemp);
			}
			UList.push_back(optemp);
			Fk[j] = 0;
		}

	}
	for (uint32_t i = 0; i < NIntegrals; i++)
	{
		free(integrals[i].val);
		free(integrals[i].k);
	}
	free(integrals);
}

static int LuaCreateContinuumSlaterIntegrals (lua_State *L)
{
	vector<Interpolation::InterpolatingFunction *> WavefunsBoundG, WavefunsBoundF, WavefunsUnboundG,WavefunsUnboundF;
	vector<int> indToK;
	vector<vector<vector<double> > > Fks; //vector<vector<double> > Fks;
	int kappa;

	uint32_t len = (uint32_t)lua_gettop(L);
	if (len == 6)
	{
		lua_check(L,1,WavefunsBoundG);
		lua_check(L,2,WavefunsBoundF);
		lua_check(L,3,indToK);
		lua_check(L,4,WavefunsUnboundG);
		lua_check(L,5,WavefunsUnboundF);
		lua_check(L,6,kappa);
		CreateContinuumSlaterIntegrals(WavefunsBoundG,WavefunsBoundF,WavefunsUnboundG,WavefunsUnboundF,indToK,kappa,Fks);
	}
	else
	{
		if (len == 4)
		{
			lua_check(L,1,WavefunsBoundG);
			lua_check(L,2,indToK);										// read l of bound states
			for (uint32_t i = 0; i < (uint32_t) indToK.size(); i++)
			{
				indToK[i] = -(indToK[i] + 1);							// transform l to kappa
			}
			lua_check(L,3,WavefunsUnboundG);
			lua_check(L,4,kappa);			// read l of continuum
			kappa = -(kappa +1);			// transform l to kappa
			CreateContinuumSlaterIntegrals(WavefunsBoundG,WavefunsUnboundG,indToK,kappa,Fks);
		}
		else
		{
			return luaL_error(L,"Error in LuaCreateContinuumSlaterIntegrals: 4 or 6 arguments expected; got %d",len);
		}
	}
	lua_push(L,Fks);
	return 1;
}

static int LuaCoolNewSelfEnergy(lua_State *L)
{
	vector<WaveFunctionType> psiCont, phiBound;
	vector<OperatorType> Uk;

	lua_check(L,1,psiCont);
	uint32_t npsiCont = (uint32_t)psiCont.size();
	lua_check(L,2,phiBound);
	uint32_t nphiBound = (uint32_t)phiBound.size();
	lua_check(L,3,Uk);
	uint32_t nUk = (uint32_t)Uk.size();
	vector<MatrixType> UkMatrix = vector<MatrixType>(nUk);

#pragma omp parallel for
	for (uint32_t k = 0; k < nUk; k++)
	{
		UkMatrix[k].Ni = npsiCont;
		UkMatrix[k].Nj = nphiBound;
		UkMatrix[k].QComplex = psiCont[0].QComplex + phiBound[0].QComplex + Uk[k].QComplex;
		InitMatrix(&UkMatrix[k]);
		PsiPrimeOperatorPsiMatrix(Uk[k],&psiCont[0],&phiBound[0],&UkMatrix[k]);
	}

	double wMin, wMax, dw,gammaSE;
	lua_check(L,4,wMin);
	lua_check(L,5,wMax);
	lua_check(L,6,dw);
	lua_check(L,7,gammaSE);
	uint32_t Nw = floor((wMax-wMin)/dw) + 1;
	vector<vector<double> > Fk;
	lua_check(L,8,Fk);
	uint32_t nContFuns = (uint32_t) Fk.size();
	vector<double> Econt, EAuger;
	lua_check(L,9,Econt);
	lua_check(L,10,EAuger);

	vector<CompactMatrixType> SelfEnergy = vector<CompactMatrixType>(Nw);
    vector<double> numerators = vector<double>(nphiBound*nContFuns*npsiCont);
    printf("Uk calculated, begin with Nums\n");
#pragma omp parallel for
    for (uint32_t a = 0; a < nphiBound; a++)
    {      
        for (uint32_t iEps = 0; iEps < nContFuns; iEps++)
        {
            for (uint32_t ipsi = 0; ipsi < npsiCont; ipsi++)
            {
                for (uint32_t k = 0; k < nUk; k++)
                {
                    numerators[a*nContFuns*npsiCont + iEps*npsiCont + ipsi] += Fk[iEps][k] * UkMatrix[k].VRe[ipsi][a]; //TODO: complex case coden                        
                }
            }
        }
    }
	printf("Nums calculated, begin with SE\n");
	fflush(stdout);
#pragma omp parallel for
	for (uint32_t i = 0; i < Nw; i++)
	{
        SelfEnergy[i].Ni = nphiBound;
        SelfEnergy[i].Nj = nphiBound;
        SelfEnergy[i].QComplex = 1;
        InitCompactMatrix(&SelfEnergy[i]);
    }
#pragma omp parallel for
    for (uint32_t a = 0; a < nphiBound; a++)
    {
        double w;
        Complex16 denom;
        denom.Im = gammaSE/2.;
        double num,num1,num2, Eac;
        
        for (uint32_t iEps = 0; iEps < nContFuns; iEps++)
        {
            for (uint32_t ipsi = 0; ipsi < npsiCont; ipsi++)
            {
                Eac = EAuger[ipsi] + Econt[iEps];
                num = numerators[a*nContFuns*npsiCont + iEps*npsiCont + ipsi];
                num *= num;
                for (uint32_t i = 0; i < Nw; i++)
                {
                    w = wMin + i*dw;
                    denom.Re = w - Eac;
                    SelfEnergy[i].V[a*nphiBound + a] += num / denom;
                }
                for (uint32_t b = a+1; b < nphiBound; b++)
                {
                    num = numerators[a*nContFuns*npsiCont + iEps*npsiCont + ipsi] * numerators[b*nContFuns*npsiCont + iEps*npsiCont + ipsi];
                    for (uint32_t i = 0; i < Nw; i++)
                    {
                        w = wMin + i*dw;
                        denom.Re = w - Eac;
                        SelfEnergy[i].V[a*nphiBound + b] += num / denom;
                        SelfEnergy[i].V[b*nphiBound + a] += num / denom;
                    }
                }

            }
        }
        
    }
	for (uint32_t k = 0; k < nUk; k++)
	{
		FreeMatrix(&UkMatrix[k]);
	}
	lua_push(L,SelfEnergy);
	return 1;
}

static int LuaCoolNewSelfEnergyDiagonal(lua_State *L)
{
	vector<WaveFunctionType> psiCont, phiBound;
	vector<OperatorType> Uk;

	lua_check(L,1,psiCont);
	uint32_t npsiCont = (uint32_t)psiCont.size();
	lua_check(L,2,phiBound);
	uint32_t nphiBound = (uint32_t)phiBound.size();
	lua_check(L,3,Uk);
	uint32_t nUk = (uint32_t)Uk.size();

	double wMin, wMax, dw,gammaSE;
	lua_check(L,4,wMin);
	lua_check(L,5,wMax);
	lua_check(L,6,dw);
	lua_check(L,7,gammaSE);
	vector<vector<double> > Fk;
	lua_check(L,8,Fk);
	vector<double> Econt, EAuger;
	lua_check(L,9,Econt);
	lua_check(L,10,EAuger);

    vector<vector<Complex16> > SelfEnergy;
    CalcCNSE(psiCont, phiBound, Uk, wMin, wMax, dw, gammaSE, Fk, Econt, EAuger, SelfEnergy);
    lua_push(L,SelfEnergy);
	return 1;
}

template <class T>
uint32_t CheckQComplex(const vector<T> &list)
{
    for (uint32_t i = 0; i < (uint32_t) list.size(); i++)
        if ( list[i].QComplex ) return 1;
    return 0;
}

void CalcCNSE(const vector<WaveFunctionType> &psiCont, const vector<WaveFunctionType> &phiBound, const vector<OperatorType> &Uk, double wMin, double wMax, double dw,
              double gammaSE, const vector<vector<double> > &Fk, const vector<double> &Econt, const vector<double> &EAuger, vector<vector<Complex16> > &SelfEnergy)
{
	uint32_t npsiCont = (uint32_t)psiCont.size();
	uint32_t nphiBound = (uint32_t)phiBound.size();
	uint32_t nUk = (uint32_t)Uk.size();
    uint32_t Nw = floor((wMax-wMin)/dw) + 1;
	uint32_t nContFuns = (uint32_t) Fk.size();
    
    vector<MatrixType> UkMatrix = vector<MatrixType>(nUk);
    SelfEnergy = vector<vector<Complex16> >(nphiBound);
    uint32_t QComplex = CheckQComplex(psiCont);
    if (!QComplex ) QComplex = CheckQComplex(phiBound);
    if (!QComplex ) QComplex = CheckQComplex(Uk);

#pragma omp parallel for
	for (uint32_t k = 0; k < nUk; k++)
	{
		UkMatrix[k].Ni = npsiCont;
		UkMatrix[k].Nj = nphiBound;
		UkMatrix[k].QComplex = QComplex;
		InitMatrix(&UkMatrix[k]);
		PsiPrimeOperatorPsiMatrix(Uk[k],&psiCont[0],&phiBound[0],&UkMatrix[k]);
	}

    if (!QComplex)
    {
    #pragma omp parallel for
        for (uint32_t a = 0; a < nphiBound; a++)
        {
            Complex16 denom;
            double num;
            double w;
            double Eac;
            SelfEnergy[a] = vector<Complex16>(Nw);
            denom.Im = gammaSE/2.;

            for (uint32_t iEps = 0; iEps < nContFuns; iEps++)
            {
                for (uint32_t ipsi = 0; ipsi < npsiCont; ipsi++)
                {
                    Eac = EAuger[ipsi] + Econt[iEps];
                    num = 0;
                    for (uint32_t k = 0; k < nUk; k++)
                    {
                        
                        num += Fk[iEps][k] * UkMatrix[k].VRe[ipsi][a];
                            
                        
                    }
                    num *= num;
                    for (uint32_t i = 0; i < Nw; i++)
                    {
                        w = wMin + i*dw;
                        denom.Re = w - Eac;
                        SelfEnergy[a][i] += num / denom;
                    }
                }
            }

        }
    }
    else
    {   
    #pragma omp parallel for
        for (uint32_t a = 0; a < nphiBound; a++)
        {
            Complex16 denom, amp;
            double num;
            double w;
            double Eac;
            SelfEnergy[a] = vector<Complex16>(Nw);
            denom.Im = gammaSE/2.;

            for (uint32_t iEps = 0; iEps < nContFuns; iEps++)
            {
                for (uint32_t ipsi = 0; ipsi < npsiCont; ipsi++)
                {
                    Eac = EAuger[ipsi] + Econt[iEps];
                    amp.Re = 0.;
                    amp.Im = 0.;
                    for (uint32_t k = 0; k < nUk; k++)
                    {
                        
                        amp.Re += Fk[iEps][k] * UkMatrix[k].VRe[ipsi][a];
                        amp.Im += Fk[iEps][k] * UkMatrix[k].VIm[ipsi][a];
                            
                        
                    }
                    num = amp.Re * amp.Re + amp.Im * amp.Im;
                    for (uint32_t i = 0; i < Nw; i++)
                    {
                        w = wMin + i*dw;
                        denom.Re = w - Eac;
                        SelfEnergy[a][i] += num / denom;
                    }
                }
            }

        }
    }
	for (uint32_t k = 0; k < nUk; k++)
	{
		FreeMatrix(&UkMatrix[k]);
	}
}

static int LuaCreateAugerMeitnerU(lua_State *L)
{
	vector<int> indToK;
	int kappa;
	uint32_t NF, Nindices, temp, numberofArgs;
	uint16_t **indices;
	vector<OperatorType> UList;
	bool QRelativistic;

	numberofArgs = (uint32_t)lua_gettop(L);
	if (numberofArgs != 4)
	{
		return luaL_error(L,"error in LuaCreateAugerMeitnerU: expected 4, got %u",numberofArgs);
	}
	lua_check(L,1,NF);
	if (!lua_istable(L,4)) // if last argument is not a table, then we are in the relavtivistic case
	{
		QRelativistic = true;
		lua_getfield(L,2,"kappas");
		if (lua_isnil(L,-1))
		{
			return luaL_error(L,"error in LuaCreateAugerMeitnerU: table of indices does not contain field kappas\n");
		}
		else
		{
			lua_check(L,-1,indToK);
			lua_pop(L,1);

			Nindices = (uint32_t)lua_objlen(L,2) + 1;

			indices = (uint16_t **) malloc(Nindices * sizeof(uint16_t *));
			if(!indices) printf("malloc fail at indices in LuaOperatorNewFullCoulomb. size = %u\n",Nindices);

			for(uint32_t i = 0; i < Nindices-1; i++ )
			{
				lua_rawgeti(L, 2, i+1);
				indices[i] = luaL_checkarray_unit16_t(L,-1,&temp);
				lua_pop(L,1);
			}
			indices[Nindices-1] = luaL_checkarray_unit16_t(L,3,&temp);
			lua_check(L,4,kappa);
			indToK.push_back(kappa);

		}
	}
	else // if last argument is a table, then we are in the non-relavtivistic case
	{
		QRelativistic = false;
		Nindices = (uint32_t)lua_objlen(L,2) + 2;

		indices = (uint16_t **) malloc(Nindices * sizeof(uint16_t *));
		if(!indices) printf("malloc fail at indices in LuaOperatorNewFullCoulomb. size = %u\n",Nindices);

		for(uint32_t i = 0; i < Nindices-2; i++ )
		{
			lua_rawgeti(L, 2, i+1);
			indices[i] = luaL_checkarray_unit16_t(L,-1,&temp);
			lua_pop(L,1);
			if (i%2 == 0)
			{
				indToK.push_back(-(((int)temp-1)/2 + 1));
			}
		}
		indices[Nindices-2] = luaL_checkarray_unit16_t(L,3,&temp);
		indices[Nindices-1] = luaL_checkarray_unit16_t(L,4,&temp);
		indToK.push_back(-(((int)temp-1)/2 + 1));


	}

	CreateAugerMeitnerU(indToK,NF,indices,UList,QRelativistic);
	lua_push(L,UList);

	for(uint32_t i = 0; i < Nindices; i++ )
	{
		free(indices[i]);
	}
	free(indices);
	return 1;

}

Complex16 Resolvent(Complex16 z, const std::vector<double> &poles, const std::vector<double> &residues,vector<vector<mypoint> > &selfenergy, uint32_t widx)
{
	Complex16 result;
	result.Re = 0;
	result.Im = 0;
	if (poles.size() < residues.size() || poles.size() != selfenergy.size())
	{
		printf("error in Resolvent: poles.size() = %u < %u = residues.size()\n",(uint32_t)poles.size(),(uint32_t)residues.size());
		fflush(stdout);
		return result;
	}
	for (uint32_t i = 0; i < (uint32_t) residues.size(); i++)
	{
		selfenergy[i][widx].y.Re = 0;	//TODO: Also add Re ??????????????????????????
		//printf("SE = %.2E\n",selfenergy[i][widx].y.Im);
		result += residues[i]/(z - poles[i] - selfenergy[i][widx].y);
	}
	return result;
}

void ConvertSpectrum( SpectraType * const G, const vector<mypoint> &spectrum)
{

  for (uint32_t i = 0; i <= G->N; i++)
  {
    G->Re[i] = spectrum[i].y.Re;
    G->Im[i] = spectrum[i].y.Im;
  }
  G->ZerothMoment = NAN;
  G->FirstMoment = NAN;
  for (uint32_t i = 0; i < 5; i++)
  {
    G->a[i] = NAN;
    G->b[i] = NAN;
  }
}

void PrintMySpectrum(const vector<mypoint> &spec)
{
	printf("Spectrum : \n");
	for (uint32_t i = 0; i < (uint32_t) spec.size(); i++)
	{
		printf("%.12E\t%.12E\t%.12E\n",spec[i].w,spec[i].y.Re,spec[i].y.Im);
	}
}
