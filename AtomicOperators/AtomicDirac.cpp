/*
 * AtomicDirac.cpp
 *
 *  Created on: Sep 10, 2018
 *      Author: brass
 */

#include <iostream>
#include <vector>
#include <map>

using std::vector;
using std::string;
using std::map;

void CreateOperatorDirac(uint32_t const NFermion,vector<vector<double> > basisfuns, vector<vector<double> > derivfuns, vector<double> x,
		SKappaGroup *groups,uint32_t Ngroups,uint32_t **indices, int Z, double c, OperatorType *op)
{
	vector<vector<double> > diracElements, dualfuns;

	//unsigned int NBasisFuns  = (unsigned int) basisfuns.size();
	unsigned int NGridpoints = (unsigned int) x.size();
	unsigned int *mj1,*mj2,n,k,j,i;
	uint16_t oadd[2] = { 0, 0 };

	CreateDualFunctions(basisfuns,x,dualfuns);

	InitOperator(op);
	InitOperatorFixedLength(op,2,8);
	op->NFermion = NFermion;
	op->NWord = Min(op->NFermion, 1 + (op->NFermion - 1) / 8);

	for ( n = 0; n < Ngroups; n++)
	{
		diracElements = CalculateDiracMatrixElements(basisfuns,derivfuns,dualfuns,groups[n].indices,groups[n].n,x,NGridpoints,groups[n].kappa,Z,c);
		//PrintMartix2(diracElements,groups[n].n);

		for ( i = 0; i < groups[n].n; i++)
		{
			for ( j = i; j < groups[n].n; j++)
			{

				mj1 = indices[groups[n].indices[i]];
				mj2 = indices[groups[n].indices[j]];
				for ( k = 0; k < (unsigned int) 2*abs(groups[n].kappa); k++)
				{

					oadd[1] = CrAnMask +  (uint16_t)mj2[k];
					oadd[0] = (uint16_t) mj1[k];
					RealOperatorAddLaderNormalOrder(op, oadd, 2, diracElements[i][j], &QDetNotOrderedFcFaBcBa);

					if ( i != j)
					{
						oadd[1] = CrAnMask +  (uint16_t) mj1[k];
						oadd[0] = (uint16_t) mj2[k];
						RealOperatorAddLaderNormalOrder(op, oadd, 2, diracElements[i][j], &QDetNotOrderedFcFaBcBa);
					}
				}
			}
		}
	}

}

void AddIndex(const uint32_t idx, const int kappa, vector<SKappaGroup> &groups)
{
	for (uint32_t i = 0; i < (uint32_t) groups.size(); i++){
		if (kappa == groups[i].kappa)
		{
			groups[i].n++;
			groups[i].indices = (unsigned int*) realloc(groups[i].indices,sizeof(unsigned int) * groups[i].n);
			groups[i].indices[groups[i].n - 1] = idx;
			return;
		}
	}
	SKappaGroup newgroup;
	newgroup.kappa = kappa;
	newgroup.n = 1;
	newgroup.indices = (unsigned int*) malloc(sizeof(unsigned int));
	newgroup.indices[0] = idx;
	groups.push_back(newgroup);
}

void GroupIndicesAccordingToKappa(const vector<int> &ind, vector<SKappaGroup> &groups)
{
	for (uint32_t i = 0; i < (uint32_t) ind.size(); i++)
	{
		AddIndex(i, ind[i], groups);
	}
}

void Print(const vector<SKappaGroup> &group)
{
	for (uint32_t i = 0; i < (uint32_t) group.size(); i++)
	{
		printf("kappa = %d\t",group[i].kappa);
		for (uint32_t j = 0; j < group[i].n; j++)
		{
			printf("%u\t",group[i].indices[j]);
		}
		printf("\n");
	}
}

static int LuaOperatorNewAtomicH0(lua_State *L)
{
	OperatorType *op;
	uint16_t **indices;
	uint32_t NF, Nindices, temp;
	uint32_t lMax = 0;
	int Z,l;
	vector<int> indToK;
	vector<BSpline::TBSpline*> G;
	lua_check(L,1,NF);
	lua_check(L,2,Z);
	Nindices = (uint32_t)lua_objlen(L,3);
	indices = (uint16_t **) malloc(Nindices * sizeof(uint16_t *));
	if(!indices) printf("malloc fail at indices in LuaOperatorNewAtomicH0. size = %u\n",Nindices);

	for(uint32_t i = 0; i < Nindices; i++ )
	{
		lua_rawgeti(L, 3, i+1);
		indices[i] = luaL_checkarray_unit16_t(L,-1,&temp);
		if (i%2 == 0)
		{
			l = ((int)temp-1)/2;
			indToK.push_back(-(l + 1));
			if (l > lMax) lMax = l;

		}
		lua_pop(L,1);
	}
	vector<vector<uint32_t> > groupedIdx = vector<vector<uint32_t> >(lMax+1);
	for (uint32_t i = 0; i < (uint32_t) indToK.size(); i++)
	{
		l = -indToK[i] - 1;
		groupedIdx[l].push_back(i);
	}
	lua_check(L,4,G);

	vector<double> Lap;
	vector<vector<double> > H0 = vector<vector<double> >(lMax+1);
	vector<vector<double> > V = vector<vector<double> >(lMax+1);
	HartreeFock::CreateLaplacian(G[0]->knots,G[0]->degree,Lap);
	uint32_t dim = (uint32_t) G[0]->y.size()-2;
	for (l = 0; l <= lMax; l++)
		HartreeFock::CreatePotentialMatrix(G[0]->knots,G[0]->degree,V[l],Z,l);

	op = (OperatorType *) lua_newuserdata(L,sizeof(OperatorType));
	luaL_setmetatable(L,"Operator_Type");
	InitOperator(op);
	InitOperatorFixedLength(op,2,8);
	op->NFermion = NF;
	op->NWord = Min(op->NFermion, 1 + (op->NFermion - 1) / 8);
	uint16_t oadd[2] = { 0, 0 };
	for (l = 0; l <= lMax; l++)
	{
		H0[l] = vector<double>(dim*dim);
		for (uint32_t a = 0; a < (uint32_t) groupedIdx[l].size(); a++)
		{
			uint32_t idx1 = groupedIdx[l][a];
			for (uint32_t b = 0; b < (uint32_t) groupedIdx[l].size(); b++)
			{
				double h0 = 0;
				uint32_t idx2 = groupedIdx[l][b];
				for (int j = 0; j < dim; j++)
				{
					for ( int i = j; i <= std::min(dim-1,j+G[0]->degree); i++)  // n = #sub-diagonals
					{
						h0 += G[idx1]->y[i+1] * (Lap[(i-j) + j*dim]/2. + V[l][(i-j) + j*dim]) * G[idx2]->y[j+1];
						if (i != j)
							h0 += G[idx1]->y[j+1] * (Lap[(i-j) + j*dim]/2. + V[l][(i-j) + j*dim]) * G[idx2]->y[i+1];
					}
				}
				for (uint32_t k = 0; k < 2*l+1; k++)
				{
					oadd[1] = CrAnMask +  (uint16_t)indices[2*idx1][k];
					oadd[0] = (uint16_t) indices[2*idx2][k];
					RealOperatorAddLaderNormalOrder(op, oadd, 2, h0, &QDetNotOrderedFcFaBcBa);

					oadd[1] = CrAnMask +  (uint16_t)indices[2*idx1+1][k];
					oadd[0] = (uint16_t) indices[2*idx2+1][k];
					RealOperatorAddLaderNormalOrder(op, oadd, 2, h0, &QDetNotOrderedFcFaBcBa);
				}

			}
		}

	}

	return 1;
}

static int LuaOperatorNewAtomicDirac(lua_State *L)
{
	OperatorType *op;
	uint32_t **indices;
	uint32_t NF, NIndices, NfunsG, NfunsF, temp;
	//SKappaGroup *groups;
	int Z;

	NF = (uint32_t)luaL_checkunsigned(L, 1);
	Z = (int) luaL_checkunsigned(L, 2);
	NIndices = (uint32_t)lua_objlen(L,3);
	NfunsG = (uint32_t)lua_objlen(L,4);
	NfunsF = (uint32_t)lua_objlen(L,5);

	if (NfunsG > NIndices)
	{
		printf("Warning: in AtomicDirac.cpp number of basis funs (%d) != number of indices (%d)\n",NfunsG , NIndices);
		fflush(stdout);
	}
	else
	{
		if (NfunsG < NIndices)
		{
			return luaL_error(L, "Errorr: in AtomicDirac.cpp number of basis funs (%d) < number of indices (%d)\n",NfunsG , NIndices);
		}
	}
	if (NfunsF != NfunsG)
	{
		return luaL_error(L,"Error in AtomicDirac.cpp: number of functions small part (%d) != large part (%d)\n",NfunsF,NfunsG);
	}
	vector<Interpolation::InterpolatingFunction *> funsG = vector<Interpolation::InterpolatingFunction *>(NfunsG);
	for(uint32_t i = 0; i < NfunsG; i++)
	{
		lua_rawgeti(L, 4, i+1);
		funsG[i] = (Interpolation::InterpolatingFunction *) luaL_checkudata(L, -1, "InterpolatingFunction_Type");
		lua_pop(L,1);
	}
	vector<Interpolation::InterpolatingFunction *> funsF = vector<Interpolation::InterpolatingFunction *>(NfunsF);
	for(uint32_t i = 0; i < NfunsF; i++)
	{
		lua_rawgeti(L, 5, i+1);
		funsF[i] = (Interpolation::InterpolatingFunction *) luaL_checkudata(L, -1, "InterpolatingFunction_Type");
		lua_pop(L,1);
	}

	/*	luaL_checktype(L, 6, LUA_TTABLE);
	NGroups = (uint32_t)lua_objlen(L,6);

	groups = (SKappaGroup *) malloc(NGroups * sizeof(SKappaGroup));
	if (!groups) printf("malloc fail at groups in LuaOperatorNewDirac. size = %u\n",NGroups);

	for(uint32_t i = 0; i < NGroups; i++ )
	{
		lua_rawgeti(L, 6, i+1);
	   	lua_rawgeti(L, -1, 1);
	   	groups[i].kappa = (int) luaL_checkunsigned(L, -1);
	   	lua_pop(L,1);
	   	lua_rawgeti(L, -1, 2);
		groups[i].indices = luaL_checkarray_unit32_t(L,-1,&(groups[i].n));
	   	lua_pop(L,1);
	   	lua_pop(L,1);
	}
	*/
	luaL_checktype(L, 3, LUA_TTABLE);
	indices = (uint32_t **) malloc(NIndices * sizeof(uint32_t *));
	if(!indices) printf("malloc fail at indices in LuaOperatorNewDirac. size = %u\n",NIndices);
	for(uint32_t i = 0; i < NIndices; i++ )
	{
	   	lua_rawgeti(L, 3, i+1);
	   	indices[i] = luaL_checkarray_unit32_t(L,-1,&temp);
	   	lua_pop(L,1);
	}
	lua_getfield(L,3,"kappas");
	if (lua_isnil(L,-1))
	{
		return luaL_error(L,"Error: NewOperator Dirac: indices have noc field kappa\n");
	}
	vector<int> indToK;
	lua_check(L,-1,indToK);
	lua_pop(L,1);

	vector<SKappaGroup> groups;
	GroupIndicesAccordingToKappa(indToK,groups);
	//Print(groups2);

	vector<double> x = funsG[0]->GetX();
	vector<vector<double> > basisfuns = PutRelativisticSplineOnGrid(funsG,funsF,x,0);
	vector<vector<double> > derivfuns = PutRelativisticSplineOnGrid(funsG,funsF,x,1);

	double c = LuaGetPhysicalConstant(L, "cAU");
	op = (OperatorType *) lua_newuserdata(L,sizeof(OperatorType));

	//CreateOperatorDirac(NF, basisfuns,derivfuns,x,groups,NGroups,indices,Z,op);
	CreateOperatorDirac(NF, basisfuns,derivfuns,x,&groups[0],(uint32_t)groups.size(),indices,Z,c,op);

	luaL_setmetatable(L,"Operator_Type");

	for (uint32_t i = 0;i < (uint32_t) groups.size(); i++ )
	{
		free(groups[i].indices);
	}
	//free(groups);
	for(uint32_t i = 0; i < NIndices; i++ )
	{
		free(indices[i]);
	}
	free(indices);

	return 1;
}


