/*
 * AtomicU.cpp
 *
 *  Created on: Sep 7, 2018
 *      Author: brass
 */

#include "../AtomicU/AtomicU.h"
#include "SlaterIntegrals.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <map>

#include "SlaterTuple.cpp"
#include "SlaterIntegrals.cpp"


using std::vector;
using std::string;
using std::ifstream;
using std::map;

static int LuaOperatorNewAtomicU (lua_State *L)
{
    OperatorType *op;
    uint32_t **tuple, NIntegrals, NGridpoints, Nindices, temp, len, NF;
    uint16_t **indices;
    vector<double> x;

    len = (uint32_t)lua_gettop(L);

    if ( len < 3 ||  len > 4)
    {
    	printf("Error in Martin.NewCoulomb: 3 or 4 arguments expected; got %d",len);
    	return -1;
    }

    NF = (uint32_t)luaL_checkunsigned(L, 1);

    vector<int> indToK;
	lua_getfield(L,2,"kappas");
	if (lua_isnil(L,-1))
	{
		lua_pop(L,1);
		Nindices = (uint32_t)lua_objlen(L,2);
		indices = (uint16_t **) malloc(Nindices * sizeof(uint16_t *));
		if(!indices) printf("malloc fail at indices in LuaOperatorNewFullCoulomb. size = %d\n",Nindices);
		for(int i = 0; i < Nindices; i++ )
		{
			lua_rawgeti(L, 2, i+1);
			lua_rawgeti(L,-1,1);
			indToK.push_back((int) luaL_checkunsigned(L, -1));
			lua_pop(L,1);
			lua_rawgeti(L,-1,2);
			indices[i] = luaL_checkarray_unit16_t(L,-1,&temp);
			lua_pop(L,1);
			lua_pop(L,1);
		}
	}
	else
	{
		uint16_t *help = luaL_checkarray_unit16_t(L,-1,&temp);
		indToK = ConvertArray(help,temp);
		free(help);

		lua_pop(L,1);
		Nindices = (uint32_t)lua_objlen(L,2);

		indices = (uint16_t **) malloc(Nindices * sizeof(uint16_t *));
		if(!indices) printf("malloc fail at indices in LuaOperatorNewFullCoulomb. size = %d\n",Nindices);

		for(int i = 0; i < Nindices; i++ )
		{
			lua_rawgeti(L, 2, i+1);
			indices[i] = luaL_checkarray_unit16_t(L,-1,&temp);
			lua_pop(L,1);
		}
	}

	uint32_t NfunsG = (uint32_t)lua_objlen(L,3);

	if (NfunsG > Nindices)
	{
		printf("Warning: Martin.Coulomb() number of basis funs (%d) != number of indices (%d)\n",NfunsG , Nindices);
		fflush(stdout);
	}
	else
	{
		if (NfunsG < Nindices)
		{
			return luaL_error(L, "Errorr: Martin.Coulomb() number of basis funs (%d) < number of indices (%d)\n",NfunsG , Nindices);
		}
	}
	if (indToK.size() != Nindices) printf("Warning: Martin.Coulomb() number of kappas (%d) != number of indices (%d)\n",indToK.size() , Nindices);

	vector<InterpolatingFunctionType *> funsG = vector<InterpolatingFunctionType *>(NfunsG);
	for(int i = 0; i < NfunsG; i++)
	{
		lua_rawgeti(L, 3, i+1);
		funsG[i] = (InterpolatingFunctionType *) luaL_checkudata(L, -1, "InterpolatingFunction_Type");
		lua_pop(L,1);
	}

    if (len == 4)
    {

    	uint32_t NfunsF = (uint32_t)lua_objlen(L,4);
    	if (NfunsF != NfunsG)
    	{
    		printf("Error in LuaOperatorNewCoulombFull: number of functions small part (%d) != large part (%d)\n",NfunsF,NfunsG);
    		return -1;
    	}
    	vector<InterpolatingFunctionType *> funsF = vector<InterpolatingFunctionType *>(NfunsF);
    	for(int i = 0; i < NfunsF; i++)
    	{
    		lua_rawgeti(L, 4, i+1);
    		funsF[i] = (InterpolatingFunctionType *) luaL_checkudata(L, -1, "InterpolatingFunction_Type");
    		lua_pop(L,1);
    	}

    	NGridpoints = ((CubicSplineType*)funsG[0]->f)->N;
    	x = ConvertArray(((CubicSplineType*)funsG[0]->f)->x,NGridpoints);

    	vector<vector<double> > wavefuns = PutRelativisticSplineOnGrid(funsG,funsF,NGridpoints);

    	op = (OperatorType *) lua_newuserdata(L,sizeof(OperatorType));
    	CreateOperatorAtomicURelativistic(NF,wavefuns,x,indToK,indices, op);
    }
    else
    {
    	printf("error: Non-relativistic Slater Integrals not yet implemented\n");
    	return -1;
    }

    luaL_setmetatable(L,"Operator_Type");


	for(int i = 0; i < Nindices; i++ )
	{
		free(indices[i]);
	}
	free(indices);

    return 1;
}

uint32_t CreateOperatorAtomicURelativistic(uint32_t const NFermion,vector<vector<double> >basisfuns, vector<double> grid,vector<int> indToK
                                               ,uint16_t **indices, OperatorType *op)
{

    OperatorType optemp;
    SSlaterIntegral *integrals;
    uint32_t NIntegrals;

    printf("%e\n",basisfuns[0][1]);

    integrals = GetSlaterIntegrals(basisfuns,grid, indToK, NIntegrals);

    int *k, i;
    unsigned int *n;
    double *Fk;

    k = integrals[0].tau.kappa;
    n = integrals[0].tau.n;
    Fk = integrals[0].val;

    CreateOperatorAtomicURelativisticOneIntegral(NFermion,k[0],k[1],k[2],k[3],indices[n[0]],indices[n[1]],indices[n[2]],indices[n[3]],Fk,n[0],n[1],n[2],n[3],op);

    for( i = 1; i < NIntegrals; i++)
    {
        k = integrals[i].tau.kappa;
        n = integrals[i].tau.n;
        Fk = integrals[i].val;

        CreateOperatorAtomicURelativisticOneIntegral(NFermion,k[0],k[1],k[2],k[3],indices[n[0]],indices[n[1]],indices[n[2]],indices[n[3]],Fk,n[0],n[1],n[2],n[3],&optemp);

        OperatorAdd(op, optemp, 1., 0.);

    }
    for ( i = 0; i < NIntegrals; i++)
    {
    	free(integrals[i].val);
    	free(integrals[i].k);
    }
    free(integrals);

    return 0;
}

void CreateOperatorAtomicURelativisticOneIntegral(uint32_t const NFermion,int const kappa1,int const kappa2, int const kappa3,int const kappa4,uint16_t const *Jindex1,
                                                  uint16_t const *Jindex2,uint16_t const *Jindex3, uint16_t const *Jindex4,double const *Fk, const unsigned int n1,
												  const unsigned int n2, const unsigned int n3, const unsigned int n4, OperatorType *const op)
{
    OperatorType optemp;

    CreateOperatorCoulombFourRelativisticShells(NFermion,kappa1,kappa2,kappa3,kappa4, Jindex1,Jindex2,Jindex3,Jindex4,Fk,op);

    if (n1 != n3)
    {
        CreateOperatorCoulombFourRelativisticShells(NFermion,kappa3,kappa2,kappa1,kappa4, Jindex3,Jindex2,Jindex1,Jindex4,Fk,&optemp);
        OperatorAdd(op, optemp, 1., 0.);
    }

    if (n2 != n4)
    {
        CreateOperatorCoulombFourRelativisticShells(NFermion,kappa1,kappa4,kappa3,kappa2, Jindex1,Jindex4,Jindex3,Jindex2,Fk,&optemp);
        OperatorAdd(op, optemp, 1., 0.);

        if (n1 != n3)
        {
            CreateOperatorCoulombFourRelativisticShells(NFermion,kappa3,kappa4,kappa1,kappa2, Jindex3,Jindex4,Jindex1,Jindex2,Fk,&optemp);
            OperatorAdd(op, optemp, 1., 0.);
        }
    }

    if ( (n1 != n2) || (n3 != n4) )
    {
        CreateOperatorCoulombFourRelativisticShells(NFermion,kappa2,kappa1,kappa4,kappa3, Jindex2,Jindex1,Jindex4,Jindex3,Fk,&optemp);
        OperatorAdd(op, optemp, 1., 0.);

        if (n1 != n3)
        {
            CreateOperatorCoulombFourRelativisticShells(NFermion,kappa2,kappa3,kappa4,kappa1, Jindex2,Jindex3,Jindex4,Jindex1,Fk,&optemp);
            OperatorAdd(op, optemp, 1., 0.);
        }

        if (n2 != n4)
        {
            CreateOperatorCoulombFourRelativisticShells(NFermion,kappa4,kappa1,kappa2,kappa3, Jindex4,Jindex1,Jindex2,Jindex3,Fk,&optemp);
            OperatorAdd(op, optemp, 1., 0.);
            if (n1 != n3)
            {
                CreateOperatorCoulombFourRelativisticShells(NFermion,kappa4,kappa3,kappa2,kappa1, Jindex4,Jindex3,Jindex2,Jindex1,Fk,&optemp);
                OperatorAdd(op, optemp, 1., 0.);
            }
        }
    }
    //FreeOperator(&optemp);

}

vector<double> ConvertArray(double *a, const int n)
{
	vector<double> res = vector<double>(n);
	for (int i = 0; i < n; i++)
	{
		res[i] = a[i];
	}
	return res;
}

vector<vector<double> > PutRelativisticSplineOnGrid(vector<InterpolatingFunctionType *> funsG, vector<InterpolatingFunctionType *> funsF, uint32_t NGridPoints)
{
	vector<vector<double> > funsCombined = vector<vector<double> >(funsG.size());
	for (int i = 0; i < funsG.size(); i++)
	{
		if (((CubicSplineType *)funsG[i]->f)->N != NGridPoints)
		{
			printf("error: Interpolating functions are on different Grids, which is not allowed, but will be allowed in the future\n");
			fflush(stdout);
			//TODO return NULL;
		}
		funsCombined[i] = vector<double>(2*NGridPoints);
		for (int j = 0; j < NGridPoints; j++)
		{
			funsCombined[i][j            ] = ((CubicSplineType *)funsG[i]->f)->y[j];
			funsCombined[i][j+NGridPoints] = ((CubicSplineType *)funsF[i]->f)->y[j];
		}
	}
	return funsCombined;

}

vector<int> ConvertArray(uint16_t *a, const int n)
{
	vector<int> res = vector<int>(n);
	for (int i = 0; i < n; i++)
	{
		res[i] = (short) a[i];
	}
	return res;
}

vector<double> ReadBasisfunction(ifstream &file, vector<double> &x)
{
	double val;
	vector<double> y;
	string word;

	file >> word;
	while (word != "#")
	{
		val = atof(word.c_str());
		if(x.size() == y.size())
			x.push_back(val);
		file >> word;
		val = atof(word.c_str());
		y.push_back(val);
		file >> word;
	}
	while ((word != "<start") && !(file.eof()))
		file >> word;
	return y;
}

string GetOrbitalName(ifstream &file)
{
	string word, name;
	while ( file >> word)
	{
		if (word == "function")
		{
			file >> word;
			file >> word;
			file >> word;
			file >> name;
			file >> word;
			name.append(word);
			name.erase(name.length()-1);
			return name;
		}
	}
	return "";
}

static int LuaReadFPLOBasisFunctions(lua_State *L)
{
	ifstream file1,file2;

	const char *name1 = luaL_checkstring(L,1);
	const char *name2 = luaL_checkstring(L,2);
	unsigned int Norbitals;
	char** orbitalNames = luaL_checkarray_string(L, 3, &Norbitals); //TODO free
	file1.open(name1);
	file2.open(name2);
	string line, word, name;
	vector<string> orbitals;
	vector<double> x, y1,y2;
	std::map<string,vector<double> > funsG, funsF;
	vector<InterpolatingFunctionType*> splinesG, splinesF;

	// Read raw data
	while (!file1.eof())
	{
		name = GetOrbitalName(file1);
		GetOrbitalName(file2);
		y1 = ReadBasisfunction(file1,x);
		y2 = ReadBasisfunction(file2,x);

		funsG[name] = y1;
		funsF[name] = y2;
	}
	file1.close();
	file2.close();

	splinesG = vector<InterpolatingFunctionType*>(Norbitals);
	splinesF = vector<InterpolatingFunctionType*>(Norbitals);

	// Push large part on stack
	lua_createtable(L,0,Norbitals);
	for (int i = 0; i < Norbitals; i++)
	{
		string charTostring(orbitalNames[i]);
		if (funsG.find(charTostring) != funsG.end())
		{

			splinesG[i] = (InterpolatingFunctionType*) lua_newuserdata(L, sizeof(InterpolatingFunctionType));
			luaL_setmetatable(L, "InterpolatingFunction_Type");
			lua_rawseti(L,-2,i+1);

			InitInterpolatingFunctionAsCubicSpline(splinesG[i], (uint32_t)x.size());
			for(int j=0;j<((CubicSplineType*)(splinesG[i]->f))->N;j++)
			{
				((CubicSplineType*)(splinesG[i]->f))->x[j]=x[j];
				((CubicSplineType*)(splinesG[i]->f))->y[j]=funsG[charTostring][j];
			}
			CreateCubicSpline((CubicSplineType*)(splinesG[i]->f), NAN, NAN);

		}
		else
		{
			printf("LuaReadFPLOBasis: Couldn't find ");
			printf("%s\n",orbitalNames[i]);
			fflush(stdout);
		}

	}
	// Push small part on stack
	lua_createtable(L,0,Norbitals);
	for (int i = 0; i < Norbitals; i++)
	{
		string charTostring(orbitalNames[i]);
		if (funsF.find(charTostring) != funsF.end())
		{

			splinesF[i] = (InterpolatingFunctionType*) lua_newuserdata(L, sizeof(InterpolatingFunctionType));
			luaL_setmetatable(L, "InterpolatingFunction_Type");
			lua_rawseti(L,-2,i+1);

			InitInterpolatingFunctionAsCubicSpline(splinesF[i], (uint32_t)x.size());
			for(int j=0;j<((CubicSplineType*)(splinesF[i]->f))->N;j++)
			{
				((CubicSplineType*)(splinesF[i]->f))->x[j]=x[j];
				((CubicSplineType*)(splinesF[i]->f))->y[j]=funsF[charTostring][j];
			}
			CreateCubicSpline((CubicSplineType*)(splinesF[i]->f), NAN, NAN);

		}
		else
		{
      printf("LuaReadFPLOBasis: Couldn't find ");
      printf("%s\n",orbitalNames[i]);
			fflush(stdout);
		}

	}


	return 2;
}
