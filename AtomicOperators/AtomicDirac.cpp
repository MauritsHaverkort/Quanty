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
		SKappaGroup *groups,uint32_t Ngroups,uint32_t **indices, int Z, OperatorType *op)
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
		diracElements = CalculateDiracMatrixElements(basisfuns,derivfuns,dualfuns,groups[n].indices,groups[n].n,x,NGridpoints,groups[n].kappa,Z);
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

static int LuaOperatorNewAtomicDirac(lua_State *L)
{
	OperatorType *op;
	uint32_t **indices;
	uint32_t NF, NGroups, NIndices, temp;
	//SKappaGroup *groups;
	int Z;

	NF = (uint32_t)luaL_checkunsigned(L, 1);
	Z = (int) luaL_checkunsigned(L, 2);
	NIndices = (uint32_t)lua_objlen(L,3);
	uint32_t NfunsG = (uint32_t)lua_objlen(L,4);
	uint32_t NfunsF = (uint32_t)lua_objlen(L,5);

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

	op = (OperatorType *) lua_newuserdata(L,sizeof(OperatorType));

	//CreateOperatorDirac(NF, basisfuns,derivfuns,x,groups,NGroups,indices,Z,op);
	CreateOperatorDirac(NF, basisfuns,derivfuns,x,&groups[0],(uint32_t)groups.size(),indices,Z,op);

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


