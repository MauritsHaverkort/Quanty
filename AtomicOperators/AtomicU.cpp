/*
 * AtomicU.cpp
 *
 *  Created on: Sep 7, 2018
 *      Author: brass
 */

#include "../AtomicOperators/AtomicU.h"

using std::vector;
using std::string;
using std::ifstream;
using std::map;

static int LuaOperatorNewAtomicU (lua_State *L)
{
    OperatorType *op;
    uint32_t Nindices, temp, len, NF;
    uint16_t **indices;
    vector<double> x;
    len = (uint32_t)lua_gettop(L);
    
    struct
    {
        bool conserve;
        bool QMF;
    } Options;
    LuaOptionsType LuaEDOptions[]={
        {"conserve","bool",&(Options.conserve)},
        {"1P-Operator","bool",&(Options.QMF)},
        {NULL, NULL, NULL}};
    Options.conserve = false;
    Options.QMF = true;
    
    bool QOptions = LuaTestOptionsList(L,-1);
    if (QOptions)
    {
         LuaCheckOptions(L,-1,LuaEDOptions);
    }

    NF = (uint32_t)luaL_checkunsigned(L, 1);
    vector<int> indToK;
    if ((len == 4 && !QOptions) || (len == 5 && QOptions) )
	{
		lua_getfield(L,2,"kappas");
		if (lua_isnil(L,-1))
		{
			lua_pop(L,1);
			Nindices = (uint32_t)lua_objlen(L,2);
			indices = (uint16_t **) malloc(Nindices * sizeof(uint16_t *));
			if(!indices) printf("malloc fail at indices in LuaOperatorNewFullCoulomb. size = %u\n",Nindices);
			for(uint32_t i = 0; i < Nindices; i++ )
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
			if(!indices) printf("malloc fail at indices in LuaOperatorNewFullCoulomb. size = %u\n",Nindices);

			for(uint32_t i = 0; i < Nindices; i++ )
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
		if (indToK.size() != Nindices) printf("Warning: Martin.Coulomb() number of kappas (%d) != number of indices (%d)\n",(uint32_t)indToK.size() , Nindices);

		vector<Interpolation::InterpolatingFunction *> funsG = vector<Interpolation::InterpolatingFunction *>(NfunsG);
		for(uint32_t i = 0; i < NfunsG; i++)
		{
			lua_rawgeti(L, 3, i+1);
			funsG[i] = (Interpolation::InterpolatingFunction *) luaL_checkudata(L, -1, "InterpolatingFunction_Type");
			lua_pop(L,1);
		}

		uint32_t NfunsF = (uint32_t)lua_objlen(L,4);
		if (NfunsF != NfunsG)
		{
			return luaL_error(L,"Error in LuaOperatorNewCoulombFull: number of functions small part (%d) != large part (%d)\n",NfunsF,NfunsG);
		}
		vector<Interpolation::InterpolatingFunction *> funsF = vector<Interpolation::InterpolatingFunction *>(NfunsF);
		for(uint32_t i = 0; i < NfunsF; i++)
		{
			lua_rawgeti(L, 4, i+1);
			funsF[i] = (Interpolation::InterpolatingFunction *) luaL_checkudata(L, -1, "InterpolatingFunction_Type");
			lua_pop(L,1);
		}

		x = funsG[0]->GetX();

		vector<vector<double> > wavefuns = PutRelativisticSplineOnGrid(funsG,funsF,x,0);
		op = (OperatorType *) lua_newuserdata(L,sizeof(OperatorType));
        if (Options.conserve)
        {
            CreateOperatorAtomicURelativisticOccupationConserving(NF,wavefuns,x,indToK,indices, op);
        }
        else
        {
            CreateOperatorAtomicURelativistic(NF,wavefuns,x,indToK,indices, op);
        }
    }
    else if ((len == 3 && !QOptions )||(len == 4 && QOptions))
    {
    	Nindices = (uint32_t)lua_objlen(L,2);

		indices = (uint16_t **) malloc(Nindices * sizeof(uint16_t *));
		if(!indices) printf("malloc fail at indices in LuaOperatorNewFullCoulomb. size = %u\n",Nindices);

		for(uint32_t i = 0; i < Nindices; i++ )
		{
			lua_rawgeti(L, 2, i+1);
			indices[i] = luaL_checkarray_unit16_t(L,-1,&temp);
			if (i%2 == 0)
			{
				indToK.push_back(-(((int)temp-1)/2 + 1));
			}
			lua_pop(L,1);
		}

		uint32_t NfunsG = (uint32_t)lua_objlen(L,3);

		if (NfunsG > Nindices/2)
		{
			printf("Warning: Martin.Coulomb() number of basis funs (%d) != number of indices (%d)\n",NfunsG , Nindices/2);
			fflush(stdout);
		}
		else
		{
			if (NfunsG < Nindices/2)
			{
				return luaL_error(L, "Errorr: Martin.Coulomb() number of basis funs (%d) < number of indices (%d)\n",NfunsG , Nindices/2);
			}
		}

		vector<Interpolation::InterpolatingFunction *> funsG = vector<Interpolation::InterpolatingFunction *>(NfunsG);
		for(uint32_t i = 0; i < NfunsG; i++)
		{
			lua_rawgeti(L, 3, i+1);
			funsG[i] = (Interpolation::InterpolatingFunction *) luaL_checkudata(L, -1, "InterpolatingFunction_Type");
			lua_pop(L,1);
		}

		op = (OperatorType *) lua_newuserdata(L,sizeof(OperatorType));
		CreateOperatorAtomicUNonRel(NF,funsG,indToK,indices,op);

    }
    else if (len == 7 && QOptions)
    {
        vector<int> indToK1, indToK2;
        lua_check(L,5,indToK1);
        lua_check(L,6,indToK2);

        Nindices = (uint32_t)lua_objlen(L,2);

        indices = (uint16_t **) malloc(Nindices * sizeof(uint16_t *));
        if(!indices) printf("malloc fail at indices in LuaOperatorNewFullCoulomb. size = %u\n",Nindices);

        for(uint32_t i = 0; i < Nindices; i++ )
        {
            lua_rawgeti(L, 2, i+1);
            indices[i] = luaL_checkarray_unit16_t(L,-1,&temp);
            lua_pop(L,1);
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
		
		vector<Interpolation::InterpolatingFunction *> funsG = vector<Interpolation::InterpolatingFunction *>(NfunsG);
		for(uint32_t i = 0; i < NfunsG; i++)
		{
			lua_rawgeti(L, 3, i+1);
			funsG[i] = (Interpolation::InterpolatingFunction *) luaL_checkudata(L, -1, "InterpolatingFunction_Type");
			lua_pop(L,1);
		}

		uint32_t NfunsF = (uint32_t)lua_objlen(L,4);
		if (NfunsF != NfunsG)
		{
			return luaL_error(L,"Error in LuaOperatorNewCoulombFull: number of functions small part (%d) != large part (%d)\n",NfunsF,NfunsG);
		}
		vector<Interpolation::InterpolatingFunction *> funsF = vector<Interpolation::InterpolatingFunction *>(NfunsF);
		for(uint32_t i = 0; i < NfunsF; i++)
		{
			lua_rawgeti(L, 4, i+1);
			funsF[i] = (Interpolation::InterpolatingFunction *) luaL_checkudata(L, -1, "InterpolatingFunction_Type");
			lua_pop(L,1);
		}

		x = funsG[0]->GetX();

		vector<vector<double> > wavefuns = PutRelativisticSplineOnGrid(funsG,funsF,x,0);
		op = (OperatorType *) lua_newuserdata(L,sizeof(OperatorType));
        if (!Options.QMF)
        {
            CreateOperatorAtomicURelativisticOccupationConserving(NF,wavefuns,x,indToK1,indToK2,indices, op);
        }
        else
        {
            OperatorType U;
            WaveFunctionType psi;
            psi.NFermion = NF;
            psi.NBoson = 0;
            psi.QComplex = 0;
            psi.NBitsKey = 16;
            WaveFunctionInit(&psi);
            std::string detstring;
            for (uint32_t i = 0; i < (uint32_t) indToK1.size(); i++) 
                for (uint32_t j = 0; j < 2*abs(indToK1[i]); j++) detstring.append("1");
            for (uint32_t i = 0; i < (uint32_t) indToK2.size(); i++)
                for (uint32_t j = 0; j < 2*abs(indToK2[i]); j++) detstring.append("0");
            uint8_t * det = StringToDet(psi.NFermion, psi.NBoson,(char*) detstring.c_str());
            uint32_t detleng = psi.NBoson + Min(psi.NFermion, 1 + (psi.NFermion - 1) / 8);
            RealWaveFunctionAddElement(&psi, det, 1);
            MatrixType Rho;
            CalculateDensityMatrix(&Rho, psi);
            CompactMatrixType rho;
            MatrixToCompactMatrix(&Rho,&rho);
            CreateOperatorAtomicURelativisticOccupationConserving(NF,wavefuns,x,indToK1,indToK2,indices, &U);
            MeanFieldOperatorOld(&U,op,rho);
            WaveFunctionFree(&psi);
            FreeOperator(&U);
            //FreeMatrix(&Rho);
            FreeCompactMatrix(&rho);
        }
    }
    else
    {
        return luaL_error(L, "Errorr: NewOperator AtomicU: wrong input parameters\n");
    }
        

    luaL_setmetatable(L,"Operator_Type");


	for(uint32_t i = 0; i < Nindices; i++ )
	{
		free(indices[i]);
	}
	free(indices);

    return 1;
}

static int LuaOperatorNewContinuumU (lua_State *L)
{
    OperatorType *op;
    uint32_t Nindices, temp, len, NF;
    uint16_t **indices;
    vector<double> x;
    len = (uint32_t)lua_gettop(L);

    if ( len < 3 ||  len > 4)
    {
    	return luaL_error(L,"Error in Martin.NewCoulomb: 3 or 4 arguments expected; got %d",len);
    }

    NF = (uint32_t)luaL_checkunsigned(L, 1);
    vector<int> indToK;
    if (len == 4)
	{
		lua_getfield(L,2,"kappas");
		if (lua_isnil(L,-1))
		{
			lua_pop(L,1);
			Nindices = (uint32_t)lua_objlen(L,2);
			indices = (uint16_t **) malloc(Nindices * sizeof(uint16_t *));
			if(!indices) printf("malloc fail at indices in LuaOperatorNewFullCoulomb. size = %u\n",Nindices);
			for(uint32_t i = 0; i < Nindices; i++ )
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
			if(!indices) printf("malloc fail at indices in LuaOperatorNewFullCoulomb. size = %u\n",Nindices);

			for(uint32_t i = 0; i < Nindices; i++ )
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
		if (indToK.size() != Nindices) printf("Warning: Martin.Coulomb() number of kappas (%d) != number of indices (%d)\n",(uint32_t)indToK.size() , Nindices);

		vector<Interpolation::InterpolatingFunction *> funsG = vector<Interpolation::InterpolatingFunction *>(NfunsG);
		for(uint32_t i = 0; i < NfunsG; i++)
		{
			lua_rawgeti(L, 3, i+1);
			funsG[i] = (Interpolation::InterpolatingFunction *) luaL_checkudata(L, -1, "InterpolatingFunction_Type");
			lua_pop(L,1);
		}

		uint32_t NfunsF = (uint32_t)lua_objlen(L,4);
		if (NfunsF != NfunsG)
		{
			return luaL_error(L,"Error in LuaOperatorNewCoulombFull: number of functions small part (%d) != large part (%d)\n",NfunsF,NfunsG);
		}
		vector<Interpolation::InterpolatingFunction *> funsF = vector<Interpolation::InterpolatingFunction *>(NfunsF);
		for(uint32_t i = 0; i < NfunsF; i++)
		{
			lua_rawgeti(L, 4, i+1);
			funsF[i] = (Interpolation::InterpolatingFunction *) luaL_checkudata(L, -1, "InterpolatingFunction_Type");
			lua_pop(L,1);
		}

		x = funsG[0]->GetX();

		vector<vector<double> > wavefuns = PutRelativisticSplineOnGrid(funsG,funsF,x,0);
		op = (OperatorType *) lua_newuserdata(L,sizeof(OperatorType));
		CreateOperatorAtomicURelativistic(NF,wavefuns,x,indToK,indices, op, true);
    }
    else
    {
    	Nindices = (uint32_t)lua_objlen(L,2);

		indices = (uint16_t **) malloc(Nindices * sizeof(uint16_t *));
		if(!indices) printf("malloc fail at indices in LuaOperatorNewFullCoulomb. size = %u\n",Nindices);

		for(uint32_t i = 0; i < Nindices; i++ )
		{
			lua_rawgeti(L, 2, i+1);
			indices[i] = luaL_checkarray_unit16_t(L,-1,&temp);
			if (i%2 == 0)
			{
				indToK.push_back(-(((int)temp-1)/2 + 1));
			}
			lua_pop(L,1);
		}

		uint32_t NfunsG = (uint32_t)lua_objlen(L,3);

		if (NfunsG > Nindices/2)
		{
			printf("Warning: Martin.Coulomb() number of basis funs (%d) != number of indices (%d)\n",NfunsG , Nindices/2);
			fflush(stdout);
		}
		else
		{
			if (NfunsG < Nindices/2)
			{
				return luaL_error(L, "Errorr: Martin.Coulomb() number of basis funs (%d) < number of indices (%d)\n",NfunsG , Nindices/2);
			}
		}

		vector<Interpolation::InterpolatingFunction *> funsG = vector<Interpolation::InterpolatingFunction *>(NfunsG);
		for(uint32_t i = 0; i < NfunsG; i++)
		{
			lua_rawgeti(L, 3, i+1);
			funsG[i] = (Interpolation::InterpolatingFunction *) luaL_checkudata(L, -1, "InterpolatingFunction_Type");
			lua_pop(L,1);
		}

		op = (OperatorType *) lua_newuserdata(L,sizeof(OperatorType));
		CreateOperatorAtomicUNonRel(NF,funsG,indToK,indices,op,true);

    }

    luaL_setmetatable(L,"Operator_Type");


	for(uint32_t i = 0; i < Nindices; i++ )
	{
		free(indices[i]);
	}
	free(indices);

    return 1;
}

static int LuaGetSlaterIntegrals (lua_State *L)
{
    uint32_t len, NIntegrals;
    vector<double> x;
    len = (uint32_t)lua_gettop(L);
    SSlaterIntegral *integrals;
    vector<string> orbitalNames;
    vector<int> indToK;

    if ( len < 2 ||  len > 3)
    {
    	return luaL_error(L,"Error in LuaGetSlaterIntegrals: 2 or 3 arguments expected; got %d",len);
    }

    orbitalNames = luaL_checkarray_string(L, 1);
    indToK = NamesToKappa(orbitalNames);

    if (len == 3)
	{
		uint32_t NfunsG = (uint32_t)lua_objlen(L,2);

		if (NfunsG > (uint32_t) indToK.size())
		{
			printf("Warning: LuaGetSlaterIntegrals number of basis funs (%d) != number of indices (%d)\n",NfunsG , (uint32_t) indToK.size());
			fflush(stdout);
		}
		else
		{
			if (NfunsG < (uint32_t) indToK.size())
			{
				return luaL_error(L, "Errorr: LuaGetSlaterIntegrals number of basis funs (%d) < number of indices (%d)\n",NfunsG , (uint32_t) indToK.size());
			}
		}

		vector<Interpolation::InterpolatingFunction *> funsG = vector<Interpolation::InterpolatingFunction *>(NfunsG);
		for(uint32_t i = 0; i < NfunsG; i++)
		{
			lua_rawgeti(L, 2, i+1);
			funsG[i] = (Interpolation::InterpolatingFunction *) luaL_checkudata(L, -1, "InterpolatingFunction_Type");
			lua_pop(L,1);
		}

		uint32_t NfunsF = (uint32_t)lua_objlen(L,3);
		if (NfunsF != NfunsG)
		{
			return luaL_error(L,"Error in LuaGetSlaterIntegrals: number of functions small part (%d) != large part (%d)\n",NfunsF,NfunsG);
		}
		vector<Interpolation::InterpolatingFunction *> funsF = vector<Interpolation::InterpolatingFunction *>(NfunsF);
		for(uint32_t i = 0; i < NfunsF; i++)
		{
			lua_rawgeti(L, 3, i+1);
			funsF[i] = (Interpolation::InterpolatingFunction *) luaL_checkudata(L, -1, "InterpolatingFunction_Type");
			lua_pop(L,1);
		}

		x = funsG[0]->GetX();

		vector<vector<double> > wavefuns = PutRelativisticSplineOnGrid(funsG,funsF,x,0);

		integrals = GetSlaterIntegrals(wavefuns,x, indToK, NIntegrals);
    }
    else
    {

		uint32_t NfunsG = (uint32_t)lua_objlen(L,2);

		if (NfunsG > (uint32_t) indToK.size())
		{
			printf("Warning: LuaGetSlaterIntegrals number of basis funs (%d) != number of indices (%d)\n",NfunsG , (uint32_t) indToK.size());
			fflush(stdout);
		}
		else
		{
			if (NfunsG < (uint32_t) indToK.size())
			{
				return luaL_error(L, "Errorr: LuaGetSlaterIntegrals number of basis funs (%d) < number of indices (%d)\n",NfunsG , (uint32_t) indToK.size());
			}
		}

		vector<Interpolation::InterpolatingFunction *> funsG = vector<Interpolation::InterpolatingFunction *>(NfunsG);
		for(uint32_t i = 0; i < NfunsG; i++)
		{
			lua_rawgeti(L, 2, i+1);
			funsG[i] = (Interpolation::InterpolatingFunction *) luaL_checkudata(L, -1, "InterpolatingFunction_Type");
			lua_pop(L,1);
		}
		integrals = GetSlaterIntegrals(funsG, indToK, NIntegrals);
    }

    lua_createtable(L,NIntegrals,0);
    for(uint32_t i = 0; i < NIntegrals; i++)
    {

    	uint32_t *n = integrals[i].tau.n;
    	double *Rk = integrals[i].val;
    	uint32_t N = integrals[i].Nval;
    	uint32_t *k = integrals[i].k;

    	string orb1 = orbitalNames[n[0]];
    	string orb2 = orbitalNames[n[1]];
    	string orb3 = orbitalNames[n[2]];
    	string orb4 = orbitalNames[n[3]];

    	Lua_Push_SlaterIntegral(L,Rk,k,N,orb1,orb2,orb3,orb4);
    	Lua_Push_SlaterIntegral(L,Rk,k,N,orb3,orb2,orb1,orb4);
    	Lua_Push_SlaterIntegral(L,Rk,k,N,orb1,orb4,orb3,orb2);
    	Lua_Push_SlaterIntegral(L,Rk,k,N,orb3,orb4,orb1,orb2);

    	Lua_Push_SlaterIntegral(L,Rk,k,N,orb2,orb1,orb4,orb3);
		Lua_Push_SlaterIntegral(L,Rk,k,N,orb2,orb3,orb4,orb1);
		Lua_Push_SlaterIntegral(L,Rk,k,N,orb4,orb1,orb2,orb3);
		Lua_Push_SlaterIntegral(L,Rk,k,N,orb4,orb3,orb2,orb1);


    }


    return 1;
}

void Lua_Push_SlaterIntegral(lua_State *L, double *Rk, uint32_t *k, uint32_t N,string orb1,string orb2,string orb3,string orb4)
{
	string label = orb1;
	label.insert(label.size()," ");
	label.insert(label.size(),orb2);
	label.insert(label.size()," ");
	label.insert(label.size(),orb3);
	label.insert(label.size()," ");
	label.insert(label.size(),orb4);
	lua_pushstring(L,label);
	lua_createtable(L,0,0);
	for(uint32_t j = 0; j < N; j++)
	{
		lua_pushnumber(L,Rk[j]);
		lua_rawseti(L,-2,k[j]);
	}
	lua_settable(L,-3);
}

void Lua_Push_SlaterIntegral(lua_State *L, double *Fk, uint32_t N)
{
	lua_createtable(L,0,0);
	for (uint32_t i = 0; i < N; i++)
	{
		lua_push(L,Fk[i]);
		lua_rawseti(L,-2,i+1);
	}
}

vector<int> NamesToKappa(const vector<string> &orbitalNames)
{
	vector<int> kappas = vector<int>(orbitalNames.size());
	for(size_t i = 0; i < orbitalNames.size(); i++)
	{
		size_t idx = orbitalNames[i].find_first_of("spdfg");
		if (idx == string::npos)
		{
			printf("Error: No valid name for an orbital\n");
		}
		else
		{
			if (idx == (orbitalNames[i].size()-1))
			{
				switch (orbitalNames[i][idx]) {
					case 's':
						kappas[i] = -1;
						break;
					case 'p':
						kappas[i] = -2;
						break;
					case 'd':
						kappas[i] = -3;
						break;
					case 'f':
						kappas[i] = -4;
						break;
					case 'g':
						kappas[i] = -5;
						break;
					default:
						printf("Sorry, we only have spdfg orbitals\n");
						break;
				}
			}
			else
			{
				switch (orbitalNames[i][idx]) {
					case 's':
						kappas[i] = -1;
						break;
					case 'p':
						if (orbitalNames[i][idx+1] == '1')
							kappas[i] = 1;
						else
							kappas[i] = -2;
						break;
					case 'd':
						if (orbitalNames[i][idx+1] == '3')
							kappas[i] = 2;
						else
							kappas[i] = -3;
						break;
					case 'f':
						if (orbitalNames[i][idx+1] == '5')
							kappas[i] = 3;
						else
							kappas[i] = -4;
						break;
					case 'g':
						if (orbitalNames[i][idx+1] == '7')
							kappas[i] = 4;
						else
							kappas[i] = -5;
						break;
					default:
						printf("Sorry, we only have spdfg orbitals\n");
						break;
				}
			}
		}
	}
	return kappas;
}

uint32_t CreateOperatorAtomicURelativistic(uint32_t const NFermion,const vector<vector<double> > &basisfuns, const vector<double> &grid, const vector<int> &indToK, uint16_t **indices, OperatorType *op, bool Qcontinuum)
{
  if(GlobalTimeOn&0x100) BeginTiming("Create Atomic U Relativistic", &GlobalTime);
  OperatorType optemp;
    SSlaterIntegral *integrals;
    uint32_t NIntegrals, i;

    integrals = GetSlaterIntegrals(basisfuns,grid, indToK, NIntegrals,Qcontinuum);

    int *k;
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

  if(GlobalTimeOn&0x100) EndTiming("Create Atomic U Relativistic", &GlobalTime);
  return 0;
}

uint32_t CreateOperatorAtomicURelativisticOccupationConserving(uint32_t const NFermion,const vector<vector<double> > &basisfuns, const vector<double> &grid, const vector<int> &indToK, uint16_t **indices, OperatorType *op)
{
  if(GlobalTimeOn&0x100) BeginTiming("Create Atomic U Relativistic", &GlobalTime);
  OperatorType optemp;
    SSlaterIntegral *integrals;
    uint32_t NIntegrals, i;

    integrals = GetSlaterIntegralsTwoFilled(basisfuns,grid, indToK, NIntegrals);

    int *k;
    unsigned int *n;
    double *Fk;

    k = integrals[0].tau.kappa;
    n = integrals[0].tau.n;
    Fk = integrals[0].val;

    CreateOperatorAtomicURelativisticOneIntegralOccupationConserving(NFermion,k[0],k[1],k[2],k[3],indices[n[0]],indices[n[1]],indices[n[2]],indices[n[3]],Fk,n[0],n[1],n[2],n[3],op);

    for( i = 1; i < NIntegrals; i++)
    {
        k = integrals[i].tau.kappa;
        n = integrals[i].tau.n;
        Fk = integrals[i].val;


        CreateOperatorAtomicURelativisticOneIntegralOccupationConserving(NFermion,k[0],k[1],k[2],k[3],indices[n[0]],indices[n[1]],indices[n[2]],indices[n[3]],Fk,n[0],n[1],n[2],n[3],&optemp);

        OperatorAdd(op, optemp, 1., 0.);

    }
    for ( i = 0; i < NIntegrals; i++)
    {
    	free(integrals[i].val);
    	free(integrals[i].k);
    }
    free(integrals);

  if(GlobalTimeOn&0x100) EndTiming("Create Atomic U Relativistic", &GlobalTime);
  return 0;
}

uint32_t CreateOperatorAtomicURelativisticOccupationConserving(uint32_t const NFermion,const vector<vector<double> > &basisfuns, const vector<double> &grid, const vector<int> &indToK1,const vector<int> &indToK2, uint16_t **indices,  OperatorType *op)
{
  if(GlobalTimeOn&0x100) BeginTiming("Create Atomic U Relativistic", &GlobalTime);
  OperatorType optemp;
    SSlaterIntegral *integrals;
    uint32_t NIntegrals, i;

    integrals = GetSlaterIntegralsOneFilled(basisfuns,grid, indToK1,indToK2, NIntegrals);

    int *k;
    unsigned int *n;
    double *Fk;

    k = integrals[0].tau.kappa;
    n = integrals[0].tau.n;
    Fk = integrals[0].val;
/*    for (uint32_t j = 0; j < integrals[0].Nval; j++)
        printf("Fk = %.16E\t",Fk[j]);
    printf("\n");*/

    CreateOperatorAtomicURelativisticOneIntegralOccupationConservingOneShell(NFermion,k[0],k[1],k[2],k[3],indices[n[0]],indices[n[1]],indices[n[2]],indices[n[3]],Fk,n[0],n[1],n[2],n[3],op);

    for( i = 1; i < NIntegrals; i++)
    {
        k = integrals[i].tau.kappa;
        n = integrals[i].tau.n;
        Fk = integrals[i].val;
 /*       for (uint32_t j = 0; j < integrals[i].Nval; j++)
        printf("Fk = %.16E\t",Fk[j]);
        printf("\n");*/

        CreateOperatorAtomicURelativisticOneIntegralOccupationConservingOneShell(NFermion,k[0],k[1],k[2],k[3],indices[n[0]],indices[n[1]],indices[n[2]],indices[n[3]],Fk,n[0],n[1],n[2],n[3],&optemp);

        OperatorAdd(op, optemp, 1., 0.);

    }
    for ( i = 0; i < NIntegrals; i++)
    {
    	free(integrals[i].val);
    	free(integrals[i].k);
    }
    free(integrals);

  if(GlobalTimeOn&0x100) EndTiming("Create Atomic U Relativistic", &GlobalTime);
  return 0;
}

uint32_t CreateListOperatorAtomicURelativistic(uint32_t const NFermion,const vector<vector<double> > &basisfuns, const vector<double> &grid, const vector<int> &indToK
                                               ,uint16_t **indices, bool Qcontinuum,lua_State *L, bool IntegralsOnly)
{

    OperatorType *optemp;
    SSlaterIntegral *integrals;
    uint32_t NIntegrals, i;

    integrals = GetSlaterIntegrals(basisfuns,grid, indToK, NIntegrals,Qcontinuum);

    int *k;
    unsigned int *n;
    double *Fk;
    uint32_t countTerms = 0;

    lua_createtable(L,0,0);
    if (IntegralsOnly && lua_istable(L,5))
    {
		vector<uint32_t> idxs;
		lua_check(L,5,idxs);
    	printf("{");
		for( i = 0; i < NIntegrals; i++)
		{
			if (InvolvesOrbital(idxs[0]-1,integrals[i].tau.n) && (InvolvesOrbital(idxs[1]-1,integrals[i].tau.n) || InvolvesOrbital(idxs[2]-1,integrals[i].tau.n)) && (InvolvesOrbital(idxs[3]-1,integrals[i].tau.n) || InvolvesOrbital(idxs[4]-1,integrals[i].tau.n)))
			{
				//printf("%u, ",i+1);
				printf("%u\t%u\t%u\t%u\t\t%u\n, ",integrals[i].tau.n[0],integrals[i].tau.n[1],integrals[i].tau.n[2],integrals[i].tau.n[3],integrals[i].k[0]);
			}
			Lua_Push_SlaterIntegral(L,integrals[i].val,integrals[i].Nval);
			lua_rawseti(L,-2,i+1);
		}
		printf("}\n");
    }
    else
    {
    	for( i = 0; i < NIntegrals; i++)
		{
			Lua_Push_SlaterIntegral(L,integrals[i].val,integrals[i].Nval);
			lua_rawseti(L,-2,i+1);
		}
    }
    if (IntegralsOnly)
    {
    	return 0;
    }
    lua_createtable(L,0,0);
    for( i = 0; i < NIntegrals; i++)
    {
    	k = integrals[i].tau.kappa;
        n = integrals[i].tau.n;
        Fk = integrals[i].val;
        for (uint32_t j = 0; j < integrals[i].Nval; j++)
        {
        	Fk[j] = 0;
        }
        for (uint32_t j = 0; j < integrals[i].Nval; j++)
        {
        	optemp = (OperatorType *) lua_newuserdata(L,sizeof(OperatorType));
        	luaL_setmetatable(L,"Operator_Type");
        	lua_rawseti(L,-2, countTerms +1);
        	countTerms++;

        	Fk[j] = 1;
        	CreateOperatorAtomicURelativisticOneIntegral(NFermion,k[0],k[1],k[2],k[3],indices[n[0]],indices[n[1]],indices[n[2]],indices[n[3]],Fk,n[0],n[1],n[2],n[3],optemp);
        	Fk[j] = 0;
        }

    }
    for ( i = 0; i < NIntegrals; i++)
    {
    	free(integrals[i].val);
    	free(integrals[i].k);
    }
    free(integrals);

    return 0;
}

uint32_t CreateOperatorAtomicUNonRel(uint32_t const NFermion,const vector<Interpolation::InterpolatingFunction *> &ipolfuns, const vector<int> &indToK
                                               ,uint16_t **indices, OperatorType *op, bool Qcontinuum)
{

    OperatorType optemp;
    SSlaterIntegral *integrals;
    uint32_t NIntegrals, i;

    integrals = GetSlaterIntegrals(ipolfuns, indToK, NIntegrals,Qcontinuum);
    //printf("NIntegrals: %d\nn1 = %d\tn2 = %d\tn3 = %d\tn4 = %d\n",NIntegrals,integrals[0].tau.n[0],integrals[0].tau.n[1],integrals[0].tau.n[2],integrals[0].tau.n[3]);
    //printf("%.16f\t%.16f\t",integrals[0].val[0],integrals[0].val[1]);
    int *k;
    unsigned int *n;
    double *Fk;

    k = integrals[0].tau.kappa;
    int l[4];
    l[0] = -k[0]-1;
    l[1] = -k[1]-1;
    l[2] = -k[2]-1;
    l[3] = -k[3]-1;
    n = integrals[0].tau.n;
    Fk = integrals[0].val;

    CreateOperatorAtomicUNonRelOneIntegral(NFermion,l[0],l[1],l[2],l[3],indices[2*n[0]],indices[2*n[0]+1],indices[2*n[1]],indices[2*n[1]+1],indices[2*n[2]],indices[2*n[2]+1],indices[2*n[3]],indices[2*n[3]+1],Fk,n[0],n[1],n[2],n[3],op);

    for( i = 1; i < NIntegrals; i++)
    {
        k = integrals[i].tau.kappa;
        l[0] = -k[0]-1;
        l[1] = -k[1]-1;
        l[2] = -k[2]-1;
        l[3] = -k[3]-1;
        n = integrals[i].tau.n;
        Fk = integrals[i].val;

        CreateOperatorAtomicUNonRelOneIntegral(NFermion,l[0],l[1],l[2],l[3],indices[2*n[0]],indices[2*n[0]+1],indices[2*n[1]],indices[2*n[1]+1],indices[2*n[2]],indices[2*n[2]+1],indices[2*n[3]],indices[2*n[3]+1],Fk,n[0],n[1],n[2],n[3],&optemp);

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
        FreeOperator(&optemp);
    }

    if (n2 != n4)
    {
        CreateOperatorCoulombFourRelativisticShells(NFermion,kappa1,kappa4,kappa3,kappa2, Jindex1,Jindex4,Jindex3,Jindex2,Fk,&optemp);
        OperatorAdd(op, optemp, 1., 0.);
        FreeOperator(&optemp);

        if (n1 != n3)
        {
            CreateOperatorCoulombFourRelativisticShells(NFermion,kappa3,kappa4,kappa1,kappa2, Jindex3,Jindex4,Jindex1,Jindex2,Fk,&optemp);
            OperatorAdd(op, optemp, 1., 0.);
            FreeOperator(&optemp);
        }
    }

    if ( (n1 != n2) || (n3 != n4) )
    {
        CreateOperatorCoulombFourRelativisticShells(NFermion,kappa2,kappa1,kappa4,kappa3, Jindex2,Jindex1,Jindex4,Jindex3,Fk,&optemp);
        OperatorAdd(op, optemp, 1., 0.);
        FreeOperator(&optemp);

        if (n1 != n3)
        {
            CreateOperatorCoulombFourRelativisticShells(NFermion,kappa2,kappa3,kappa4,kappa1, Jindex2,Jindex3,Jindex4,Jindex1,Fk,&optemp);
            OperatorAdd(op, optemp, 1., 0.);
            FreeOperator(&optemp);
        }

        if (n2 != n4)
        {
            CreateOperatorCoulombFourRelativisticShells(NFermion,kappa4,kappa1,kappa2,kappa3, Jindex4,Jindex1,Jindex2,Jindex3,Fk,&optemp);
            OperatorAdd(op, optemp, 1., 0.);
            FreeOperator(&optemp);
            if (n1 != n3)
            {
                CreateOperatorCoulombFourRelativisticShells(NFermion,kappa4,kappa3,kappa2,kappa1, Jindex4,Jindex3,Jindex2,Jindex1,Fk,&optemp);
                OperatorAdd(op, optemp, 1., 0.);
                FreeOperator(&optemp);
            }
        }
    }
    //FreeOperator(&optemp);  //TODO if i uncomment this I get an invalid free - otherwhise heap is lost

}

void CreateOperatorAtomicURelativisticOneIntegralOccupationConserving(uint32_t const NFermion,int const kappa1,int const kappa2, int const kappa3,int const kappa4,uint16_t const *Jindex1,
                                                  uint16_t const *Jindex2,uint16_t const *Jindex3, uint16_t const *Jindex4,double const *Fk, const unsigned int n1,
												  const unsigned int n2, const unsigned int n3, const unsigned int n4, OperatorType *const op)
{
    OperatorType optemp;

    CreateOperatorCoulombFourRelativisticShells(NFermion,kappa1,kappa2,kappa3,kappa4, Jindex1,Jindex2,Jindex3,Jindex4,Fk,op);


    if ( (n1 != n2) && (n3 != n4) )
    {
        CreateOperatorCoulombFourRelativisticShells(NFermion,kappa2,kappa1,kappa4,kappa3, Jindex2,Jindex1,Jindex4,Jindex3,Fk,&optemp);
        OperatorAdd(op, optemp, 1., 0.);
        FreeOperator(&optemp);

    }
    //FreeOperator(&optemp);  //TODO if i uncomment this I get an invalid free - otherwhise heap is lost

}

void CreateOperatorAtomicURelativisticOneIntegralOccupationConservingOneShell(uint32_t const NFermion,int const kappa1,int const kappa2, int const kappa3,int const kappa4,uint16_t const *Jindex1,
                                                  uint16_t const *Jindex2,uint16_t const *Jindex3, uint16_t const *Jindex4,double const *Fk, const unsigned int n1,
												  const unsigned int n2, const unsigned int n3, const unsigned int n4, OperatorType *const op)
{
    OperatorType optemp;

    CreateOperatorCoulombFourRelativisticShells(NFermion,kappa1,kappa2,kappa3,kappa4, Jindex1,Jindex2,Jindex3,Jindex4,Fk,op);
    CreateOperatorCoulombFourRelativisticShells(NFermion,kappa2,kappa1,kappa4,kappa3, Jindex2,Jindex1,Jindex4,Jindex3,Fk,&optemp);
    OperatorAdd(op, optemp, 1., 0.);
    FreeOperator(&optemp);


    if ( (n3 != n2) && (n1 == n4) )
    {
        CreateOperatorCoulombFourRelativisticShells(NFermion,kappa1,kappa3,kappa2,kappa4, Jindex1,Jindex3,Jindex2,Jindex4,Fk,&optemp);
        OperatorAdd(op, optemp, 1., 0.);
        FreeOperator(&optemp);
        CreateOperatorCoulombFourRelativisticShells(NFermion,kappa3,kappa1,kappa4,kappa2, Jindex3,Jindex1,Jindex4,Jindex2,Fk,&optemp);
        OperatorAdd(op, optemp, 1., 0.);
        FreeOperator(&optemp);

    }
    if ( (n4 != n2) && (n1 == n3) )
    {
        CreateOperatorCoulombFourRelativisticShells(NFermion,kappa1,kappa4,kappa3,kappa2, Jindex1,Jindex4,Jindex3,Jindex2,Fk,&optemp);
        OperatorAdd(op, optemp, 1., 0.);
        FreeOperator(&optemp);
        CreateOperatorCoulombFourRelativisticShells(NFermion,kappa4,kappa1,kappa2,kappa3, Jindex4,Jindex1,Jindex2,Jindex3,Fk,&optemp);
        OperatorAdd(op, optemp, 1., 0.);
        FreeOperator(&optemp);
    }
    //FreeOperator(&optemp);  //TODO if i uncomment this I get an invalid free - otherwhise heap is lost

}

void CreateOperatorAtomicUNonRelOneIntegral(uint32_t const NFermion,int const l1,int const l2, int const l3,int const l4,
		uint16_t const *Lindex1Up, uint16_t const *Lindex1Dn, uint16_t const *Lindex2Up, uint16_t const *Lindex2Dn,
		uint16_t const *Lindex3Up, uint16_t const *Lindex3Dn, uint16_t const *Lindex4Up, uint16_t const *Lindex4Dn,
		double const *Fk, const unsigned int n1, const unsigned int n2, const unsigned int n3, const unsigned int n4, OperatorType *const op)
{
    OperatorType optemp;

    CreateOperatorCoulombFourShells(NFermion,l1,l2,l3,l4, Lindex1Up,Lindex1Dn, Lindex2Up,Lindex2Dn, Lindex3Up,Lindex3Dn, Lindex4Up,Lindex4Dn,Fk,op);

    if (n1 != n3)
    {
        CreateOperatorCoulombFourShells(NFermion,l3,l2,l1,l4, Lindex3Up,Lindex3Dn,Lindex2Up,Lindex2Dn,Lindex1Up,Lindex1Dn,Lindex4Up,Lindex4Dn,Fk,&optemp);
        OperatorAdd(op, optemp, 1., 0.);
    }

    if (n2 != n4)
    {
        CreateOperatorCoulombFourShells(NFermion,l1,l4,l3,l2, Lindex1Up,Lindex1Dn,Lindex4Up,Lindex4Dn,Lindex3Up,Lindex3Dn,Lindex2Up,Lindex2Dn,Fk,&optemp);
        OperatorAdd(op, optemp, 1., 0.);

        if (n1 != n3)
        {
            CreateOperatorCoulombFourShells(NFermion,l3,l4,l1,l2, Lindex3Up,Lindex3Dn,Lindex4Up,Lindex4Dn,Lindex1Up,Lindex1Dn,Lindex2Up,Lindex2Dn,Fk,&optemp);
            OperatorAdd(op, optemp, 1., 0.);
        }
    }

    if ( (n1 != n2) || (n3 != n4) )
    {
        CreateOperatorCoulombFourShells(NFermion,l2,l1,l4,l3, Lindex2Up,Lindex2Dn,Lindex1Up,Lindex1Dn,Lindex4Up,Lindex4Dn,Lindex3Up,Lindex3Dn,Fk,&optemp);
        OperatorAdd(op, optemp, 1., 0.);

        if (n1 != n3)
        {
            CreateOperatorCoulombFourShells(NFermion,l2,l3,l4,l1, Lindex2Up,Lindex2Dn,Lindex3Up,Lindex3Dn,Lindex4Up,Lindex4Dn,Lindex1Up,Lindex1Dn,Fk,&optemp);
            OperatorAdd(op, optemp, 1., 0.);
        }

        if (n2 != n4)
        {
            CreateOperatorCoulombFourShells(NFermion,l4,l1,l2,l3, Lindex4Up,Lindex4Dn,Lindex1Up,Lindex1Dn,Lindex2Up,Lindex2Dn,Lindex3Up,Lindex3Dn,Fk,&optemp);
            OperatorAdd(op, optemp, 1., 0.);
            if (n1 != n3)
            {
                CreateOperatorCoulombFourShells(NFermion,l4,l3,l2,l1, Lindex4Up,Lindex4Dn,Lindex3Up,Lindex3Dn,Lindex2Up,Lindex2Dn,Lindex1Up,Lindex1Dn,Fk,&optemp);
                OperatorAdd(op, optemp, 1., 0.);
            }
        }
    }
    //FreeOperator(&optemp);  //TODO if i uncomment this I get an invalid free - otherwhise heap is lost

}

