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

	unsigned int NBasisFuns  = basisfuns.size();
	unsigned int NGridpoints = x.size();
	unsigned int *mj1,*mj2,n,k,j,i;
	uint16_t oadd[2] = { 0, 0 };

	dualfuns = CreateDualFunctions(NBasisFuns, NGridpoints,basisfuns,x);

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
				for ( k = 0; k < 2*abs(groups[n].kappa); k++)
				{

					oadd[1] = CrAnMask +  mj2[k];
					oadd[0] = mj1[k];
					RealOperatorAddLaderNormalOrder(op, oadd, 2, diracElements[i][j], &QDetNotOrderedFcFaBcBa);

					if ( i != j)
					{
						oadd[1] = CrAnMask +  mj1[k];
						oadd[0] = mj2[k];
						RealOperatorAddLaderNormalOrder(op, oadd, 2, diracElements[i][j], &QDetNotOrderedFcFaBcBa);
					}
				}
			}
		}
	}

}

static int LuaOperatorNewAtomicDirac(lua_State *L)
{
	OperatorType *op;
	uint32_t **indices;
	uint32_t NF, NGroups, NGridpoints, NIndices, temp;
	SKappaGroup *groups;
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
		printf("Error in AtomicDirac.cpp: number of functions small part (%d) != large part (%d)\n",NfunsF,NfunsG);
		return -1;
	}
	vector<Interpolation::InterpolatingFunction *> funsG = vector<Interpolation::InterpolatingFunction *>(NfunsG);
	for(int i = 0; i < NfunsG; i++)
	{
		lua_rawgeti(L, 4, i+1);
		funsG[i] = (Interpolation::InterpolatingFunction *) luaL_checkudata(L, -1, "InterpolatingFunction_Type");
		lua_pop(L,1);
	}
	vector<Interpolation::InterpolatingFunction *> funsF = vector<Interpolation::InterpolatingFunction *>(NfunsF);
	for(int i = 0; i < NfunsF; i++)
	{
		lua_rawgeti(L, 5, i+1);
		funsF[i] = (Interpolation::InterpolatingFunction *) luaL_checkudata(L, -1, "InterpolatingFunction_Type");
		lua_pop(L,1);
	}

	luaL_checktype(L, 6, LUA_TTABLE);
	NGroups = (uint32_t)lua_objlen(L,6);

	groups = (SKappaGroup *) malloc(NGroups * sizeof(SKappaGroup));
	if (!groups) printf("malloc fail at groups in LuaOperatorNewDirac. size = %d\n",NGroups);

	for(int i = 0; i < NGroups; i++ )
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

	luaL_checktype(L, 3, LUA_TTABLE);
	indices = (uint32_t **) malloc(NIndices * sizeof(uint32_t *));
	if(!indices) printf("malloc fail at indices in LuaOperatorNewDirac. size = %d\n",NIndices);
	for(int i = 0; i < NIndices; i++ )
	{
	   	lua_rawgeti(L, 3, i+1);
	   	indices[i] = luaL_checkarray_unit32_t(L,-1,&temp);
	   	lua_pop(L,1);
	}

	vector<double> x = funsG[0]->GetX();
	vector<vector<double> > basisfuns = PutRelativisticSplineOnGrid(funsG,funsF,x,0);
	vector<vector<double> > derivfuns = PutRelativisticSplineOnGrid(funsG,funsF,x,1);
/*
	printf("%d\n",derivfuns.size());
	printf("{");
	for (int i = 0; i < derivfuns.size(); i++)
	{
		printf("{");
		for(int j = 0; j < derivfuns[i].size()-1; j++)
		{
			printf("%.16E, ",derivfuns[i][j]);
		}
		printf("%.16E},\n",derivfuns[i][derivfuns[i].size()-1]);
	}
	printf("}\n");
*/
	op = (OperatorType *) lua_newuserdata(L,sizeof(OperatorType));

	//CreateOperatorDirac(NF, basisfuns,derivfuns,x,groups,NGroups,indices,Z,op);
	CreateOperatorDirac(NF, basisfuns,derivfuns,x,groups,NGroups,indices,Z,op);

	luaL_setmetatable(L,"Operator_Type");

	for (int i = 0;i < NGroups; i++ )
	{
		free(groups[i].indices);
	}
	free(groups);
	for(int i = 0; i < NIndices; i++ )
	{
		free(indices[i]);
	}
	free(indices);

	return 1;
}
