/*
 * AuxFunctions.cpp
 *
 *  Created on: Feb 11, 2019
 *      Author: brass
 */

void RotateWaveFunctions(const WaveFunctionType* psi,const uint32_t N,const MatrixType &U, WaveFunctionType* out)
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

void RotateWaveFunctions(const WaveFunctionType* psi,const uint32_t N,const CompactMatrixType &U, WaveFunctionType* out)
{
	if (U.Ni != N || U.Nj != N)
	{
		printf("error in RotateWaveFunctions: Dimensions of U and psi do not match\n");
		return;
	}
	//out = (WaveFunctionType*) malloc(N*sizeof(WaveFunctionType));
	for (uint32_t i = 0; i < N; i++)
	{
		out[i].NFermion = psi[i].NFermion;
		out[i].NBoson   = psi[i].NBoson;
		out[i].NBitsKey = psi[i].NBitsKey;
		if (U.QComplex)
		{
			out[i].QComplex = 1;
			if(WaveFunctionInit(&(out[i])))
			{
				printf("WaveFunctionInit failed in ResonantSpectra 00\n");
				return;
			}
			for (uint32_t j = 0; j < N; j++)
			{
				WaveFunctionAddComplexMultipliedUnionBasis(&(out[i]),psi[j],U.V[i*U.Nj + j].Re,U.V[i*U.Nj + j].Im);
			}
			WaveFunctionTruncateBasis(&(out[i]),DBL_EPSILON);

		}
		else
		{
			out[i].QComplex = psi[i].QComplex;
			if(WaveFunctionInit(&(out[i])))
			{
				printf("WaveFunctionInit failed in ResonantSpectra 00");
				return;
			}
			for (uint32_t j = 0; j < N; j++)
			{
				WaveFunctionAddMultipliedUnionBasis(&(out[i]),psi[j],U.VRe[i*U.Nj + j]);
			}
			WaveFunctionTruncateBasis(&(out[i]),DBL_EPSILON);
		}
	}
}

void RotateWaveFunctions(WaveFunctionType** psi,const uint32_t N, const vector<uint32_t> &BlockSize, const CompactMatrixType &U, vector<WaveFunctionType> &out)
{
	uint32_t dim = 0;
    for (uint32_t i = 0; i < N; i++)
    {
        dim += BlockSize[i];
        printf("BlockSize = %u\n",BlockSize[i]);
    }
    printf("dim = %u\nN=%u\n",dim,N);
    if (U.Ni < dim || U.Nj < dim)
	{
		printf("error in RotateWaveFunctions: Dimensions of U and psi do not match\n");
		return;
	}
	out = vector<WaveFunctionType>(U.Ni);
	for (uint32_t i = 0; i < U.Ni; i++)
	{
		out[i].NFermion = psi[0][0].NFermion;
		out[i].NBoson   = psi[0][0].NBoson;
		out[i].NBitsKey = psi[0][0].NBitsKey;
		if (U.QComplex)
		{
			out[i].QComplex = 1;
			if(WaveFunctionInit(&(out[i])))
			{
				printf("WaveFunctionInit failed in ResonantSpectra 00\n");
				return;
			}
			uint32_t offset = 0;
			for (uint32_t j = 0; j < N; j++)
			{
                for (uint32_t k = 0; k < BlockSize[j]; k++ )
                {
                    WaveFunctionAddComplexMultipliedUnionBasis(&(out[i]),psi[j][k],U.V[i*U.Nj + offset + k].Re,U.V[i*U.Nj + offset + k].Im);
                }
                offset+=BlockSize[j];
			}
			WaveFunctionTruncateBasis(&(out[i]),DBL_EPSILON);

		}
		else
		{
			out[i].QComplex = psi[0][0].QComplex;
			if(WaveFunctionInit(&(out[i])))
			{
				printf("WaveFunctionInit failed in ResonantSpectra 00");
				return;
			}
			uint32_t offset = 0;
			for (uint32_t j = 0; j < N; j++)
			{
				for (uint32_t k = 0; k < BlockSize[j]; k++ )
                {
                    WaveFunctionAddMultipliedUnionBasis(&(out[i]),psi[j][k],U.VRe[i*U.Nj + offset + k]);
                }
                offset+=BlockSize[j];
			}
			WaveFunctionTruncateBasis(&(out[i]),DBL_EPSILON);
		}
	}
}

void RotateWaveFunctions(WaveFunctionType** psi,const uint32_t N, const uint32_t BlockSize, const vector<uint32_t> &bs, const CompactMatrixType &U, vector<WaveFunctionType> &out)
{
	uint32_t dim = N*BlockSize;

    printf("dim = %u\nN=%u\n",dim,N);
    if (U.Ni < dim || U.Nj < dim)
	{
		printf("error in RotateWaveFunctions: Dimensions of U and psi do not match\n");
		return;
	}
	out = vector<WaveFunctionType>(U.Ni);
	for (uint32_t i = 0; i < U.Ni; i++)
	{
		out[i].NFermion = psi[0][0].NFermion;
		out[i].NBoson   = psi[0][0].NBoson;
		out[i].NBitsKey = psi[0][0].NBitsKey;
		if (U.QComplex)
		{
			out[i].QComplex = 1;
			if(WaveFunctionInit(&(out[i])))
			{
				printf("WaveFunctionInit failed in ResonantSpectra 00\n");
				return;
			}
			for (uint32_t j = 0; j < N; j++)
			{
                for (uint32_t k = 0; k < bs[j]; k++ )
                {
                    WaveFunctionAddComplexMultipliedUnionBasis(&(out[i]),psi[j][k],U.V[i*U.Nj + j*BlockSize + k].Re,U.V[i*U.Nj + j*BlockSize + k].Im);
                }
                
			}
			WaveFunctionTruncateBasis(&(out[i]),DBL_EPSILON);

		}
		else
		{
			out[i].QComplex = psi[0][0].QComplex;
			if(WaveFunctionInit(&(out[i])))
			{
				printf("WaveFunctionInit failed in ResonantSpectra 00");
				return;
			}
			for (uint32_t j = 0; j < N; j++)
			{
				for (uint32_t k = 0; k < bs[j]; k++ )
                {
                    WaveFunctionAddMultipliedUnionBasis(&(out[i]),psi[j][k],U.VRe[i*U.Nj + j*BlockSize + k]);
                }
			}
			WaveFunctionTruncateBasis(&(out[i]),DBL_EPSILON);
		}
	}
}

void RotateWaveFunctions(const WaveFunctionType* psi,const uint32_t N,const MatrixType &U, vector<WaveFunctionType> &out)
{
    if (U.Ni != N || U.Nj != N)
	{
		printf("error in RotateWaveFunctions: Dimensions of U and psi do not match\n");
		return;
	}
	out = vector<WaveFunctionType>(N);
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

void BlockTriMatrixToCompact(const BlockTriDiagonalMatrixType &BTM, const vector<uint32_t> &BlockSize, CompactMatrixType &CM)
{
    uint32_t dim = 0;
    for (uint32_t i = 0; i < (uint32_t) BlockSize.size(); i++)
    {
        dim += BlockSize[i];
        printf("BlockSize[%u] = %u\n",i,BlockSize[i]);
    }
    CM.QComplex = BTM.QComplex;
    CM.Ni = dim;
    CM.Nj = dim;
    InitCompactMatrix(&CM);
    if (BTM.QComplex)
    {
        for (uint32_t k = 0; k < BTM.N; k++)
        {
            uint32_t offset = 0;
            for (uint32_t i = 0; i < BlockSize[k]; i++)
            {
                for (uint32_t j = 0; j < BlockSize[k]; j++)
                {
                    CM.V[(i+offset)*dim + (j+offset)] = BTM.a[k*BTM.BlockSize*BTM.BlockSize+i*BTM.BlockSize+j];
                }
            }
            offset += BlockSize[k];
        }
    }
    else
    {
        uint32_t offset = 0;
        for (uint32_t k = 1; k < BTM.N; k++)
        {
            for (uint32_t i = 0; i < BlockSize[k-1]; i++)
            {
                for (uint32_t j = 0; j < BlockSize[k-1]; j++)
                {
                    CM.VRe[(i+offset)*dim + (j+offset)] = BTM.aRe[k*BTM.BlockSize*BTM.BlockSize+i*BTM.BlockSize+j];
                }
            }
            offset += BlockSize[k-1];
        }
        uint32_t offseti = BlockSize[0];
        uint32_t offsetj = 0;
        for (uint32_t k = 1; k < BTM.N-1; k++)
        {
            for (uint32_t i = 0; i < BlockSize[k]; i++)
            {
                for (uint32_t j = 0; j < BlockSize[k-1]; j++)
                {
                    CM.VRe[(i+offseti)*dim + (j+offsetj)] = BTM.bRe[k*BTM.BlockSize*BTM.BlockSize+j*BTM.BlockSize+i];
                    CM.VRe[(j+offsetj)*dim + (i+offseti)] = BTM.bRe[k*BTM.BlockSize*BTM.BlockSize+j*BTM.BlockSize+i];
                }
            }
            offseti += BlockSize[k];
            offsetj += BlockSize[k-1];
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

template <typename T>
void AddVectors(vector<T> &u, const vector<T> &v)
{
    if (u.size() != v.size())
    {
        printf("Error adding Vectors: different lengths: %u != %u\n",(uint32_t) u.size(),(uint32_t) v.size());
    }
    else
    {
        for (uint32_t i = 0; i < (uint32_t) u.size(); i++)
        {
            u[i] += v[i];
        }
    }
}

uint32_t GetIndices(lua_State *L,int arg, vector<int> &indToK, uint16_t** &indices)
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
	return Nindices;
}

static int LuaLinCombine(lua_State *L)
{
    vector<WaveFunctionType> basis;
    WaveFunctionType psi;
    vector<double> vec;
    lua_check(L,1,basis);
    lua_check(L,2,vec);
    psi.NFermion = basis[0].NFermion;
    psi.NBoson = basis[0].NBoson;
    psi.QComplex = basis[0].QComplex;
    psi.NBitsKey = basis[0].NBitsKey;
    WaveFunctionInit(&psi);
    for (uint32_t i = 0; i < (uint32_t) vec.size(); i++)
    {
        RealWaveFunctionAddMultipliedUnionBasis(&psi,basis[i],vec[i]);
    }
    WaveFunctionNormalize(&psi);
    lua_push(L,&psi);
    return 1;
}

static int LuaBasisChange(lua_State *L)
{
    vector<WaveFunctionType> basis, psi;
    vector<vector<double> >U;
    lua_check(L,1,basis);
    lua_check(L,2,U);
    psi = vector<WaveFunctionType>( U.size());
#pragma omp parallel for
    for (uint32_t j = 0; j < (uint32_t) U.size(); j++)
    {
        psi[j].NFermion = basis[0].NFermion;
        psi[j].NBoson = basis[0].NBoson;
        psi[j].QComplex = basis[0].QComplex;
        psi[j].NBitsKey = basis[0].NBitsKey;
        WaveFunctionInit(&(psi[j]));
        for (uint32_t i = 0; i < (uint32_t) U[j].size(); i++)
        {
            RealWaveFunctionAddMultipliedUnionBasis(&(psi[j]),basis[i],U[j][i]);
        }
        WaveFunctionNormalize(&(psi[j]));
		WaveFunctionTruncateBasis(&(psi[j]),DBL_EPSILON);

    }
    lua_push(L,psi);
    return 1;
}
























