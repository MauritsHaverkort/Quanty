/*
 * AuxFunctions.cpp
 *
 *  Created on: Feb 11, 2019
 *      Author: brass
 */

void RotateWaveFunctions(const WaveFunctionType* psi,const uint32_t N,const MatrixType &U, WaveFunctionType* &out)
{
	if (U.Ni != N || U.Nj != N)
	{
		printf("error in RotateWaveFunctions: Dimensions of U and psi do not match\n");
		return;
	}
	out = (WaveFunctionType*) malloc(N*sizeof(WaveFunctionType));
	for (uint32_t i = 0; i < N; i++)
	{
		out[i].NFermion = psi[i].NFermion;
		out[i].NBoson   = psi[i].NBoson;
		out[i].NBitsKey = psi[i].NBitsKey;
		if (psi[i].QComplex || U.QComplex)
		{
			out[i].QComplex = 1;
			if(WaveFunctionInit(&(out[i])))
			{
				printf("WaveFunctionInit failed in ResonantSpectra 00");
				return;
			}
			for (uint32_t j = 0; j < N; j++)
			{
				WaveFunctionAddComplexMultipliedUnionBasis(&(out[i]),psi[j],U.VRe[i][j],U.VIm[i][j]);
			}
			WaveFunctionTruncateBasis(&(out[i]),DBL_EPSILON);

		}
		else
		{
			out[i].QComplex = 0;
			if(WaveFunctionInit(&(out[i])))
			{
				printf("WaveFunctionInit failed in ResonantSpectra 00");
				return;
			}
			for (uint32_t j = 0; j < N; j++)
			{
				WaveFunctionAddMultipliedUnionBasis(&(out[i]),psi[j],U.VRe[i][j]);
			}
			WaveFunctionTruncateBasis(&(out[i]),DBL_EPSILON);
		}
	}
}

void CreateSpectrum(const OperatorType &H,const OperatorType &T, WaveFunctionType * const psi, const uint32_t Npsi, double epsilon, BlockTriDiagonalMatrixType &GTri, const uint32_t NTri)
{
	WaveFunctionType * Tpsi = (WaveFunctionType *) malloc(sizeof(WaveFunctionType) * Npsi) ;
	WaveFunctionType * TpsiSameBasis = (WaveFunctionType *) malloc(sizeof(WaveFunctionType) * Npsi) ;
	for (uint32_t i = 0; i < Npsi; i++)
	{
		OperatorPsi(T,psi[i],&Tpsi[i]);
	}
    WaveFunctionCopy(Tpsi[0], &TpsiSameBasis[0]);
    WaveFunctionFree(&(Tpsi[0]));
    for (uint32_t i = 1; i < Npsi; i++)
    {
    	if(WaveFunctionAddMultipliedUnionBasis(&TpsiSameBasis[0], Tpsi[i], 0))
    	{
    		printf("malloc failed in CreateSpectrum\n");
    		fflush(stdout);
    		return;
    	}
    }
    for(uint32_t i = 1; i < Npsi; i++)
    {
		if(WaveFunctionInitCopyBasis(&(TpsiSameBasis[i]), TpsiSameBasis[0]))
		{
    		printf("malloc failed in CreateSpectrum\n");
    		fflush(stdout);
    		return;
		}
		WaveFunctionAddConserveBasis(&(TpsiSameBasis[i]), Tpsi[i]);
		WaveFunctionFree(&(Tpsi[i]));
    }
    free(Tpsi);
	GTri.BlockSize = Npsi;
	GTri.N = NTri;
	GTri.QComplex=0;
	if(H.QComplex) GTri.QComplex=1;
	for(uint32_t i = 0; i < GTri.BlockSize; i++)
	{
	if(TpsiSameBasis[i].QComplex) GTri.QComplex=1;
	}
	InitBlockTridiagonalMatrix(&GTri);
	GTri.mu = 0; //psiEn[0];
	LanczosBlockTriDiagonalize(H, TpsiSameBasis, epsilon, &GTri, NAN);


}

void GetIndices(lua_State *L,int arg, vector<int> &indToK, uint16_t** &indices)
{
	lua_getfield(L,arg,"kappas");
	lua_check(L,-1,indToK);
	lua_pop(L,1);
	uint32_t Nindices = (uint32_t)lua_objlen(L,arg);

	indices = (uint16_t **) malloc(Nindices * sizeof(uint16_t *));

	if(!indices) printf("malloc fail at indices in GetIndices, size = %u\n",Nindices);

	for(uint32_t i = 0; i < Nindices; i++ )
	{
		lua_rawgeti(L, arg, i+1);
		uint32_t temp;
		indices[i] = luaL_checkarray_unit16_t(L,-1,&temp);
		lua_pop(L,1);
	}
}
