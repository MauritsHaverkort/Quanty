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

void CreateContinuumSlaterIntegrals(vector<Interpolation::InterpolatingFunction *> &WavefunsBoundG, vector<Interpolation::InterpolatingFunction *> &WavefunsBoundF,
							vector<Interpolation::InterpolatingFunction *> &WavefunsUnboundG, vector<Interpolation::InterpolatingFunction *> &WavefunsUnboundF,
							vector<int> &indToK, const int kappa, vector<vector<double> > &Fks)
{
	SSlaterIntegral *si;
	uint32_t kMax, n, NIntegrals, NBound, NUnbound;
	double **kernels;
	vector<vector<double> > dualfunsBound, dualfunsUnbound;
	vector<OperatorType> UList;
	vector<double> poles;
	CompactMatrixType T;


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
	WavefunsBoundG.pop_back();
	WavefunsBoundF.pop_back();

}

static int LuaCreateContinuumSlaterIntegrals (lua_State *L)
{
	vector<Interpolation::InterpolatingFunction *> WavefunsBoundG, WavefunsBoundF, WavefunsUnboundG,WavefunsUnboundF;
	vector<int> indToK;
	vector<vector<double> > Fks;
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
			lua_check(L,2,indToK);
			lua_check(L,3,WavefunsUnboundG);
			lua_check(L,4,kappa);
		}
		else
		{
			return luaL_error(L,"Error in LuaCreateContinuumSlaterIntegrals: 4 or 6 arguments expected; got %d",len);
		}
	}
	lua_push(L,Fks);
	return 1;
}


