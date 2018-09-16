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

    if ( len < 3 ||  len > 4)
    {
    	printf("Error in Martin.NewCoulomb: 3 or 4 arguments expected; got %u",len);
    	return -1;
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
			printf("Warning: Martin.Coulomb() number of basis funs (%u) != number of indices (%u)\n",NfunsG , Nindices);
			fflush(stdout);
		}
		else
		{
			if (NfunsG < Nindices)
			{
				return luaL_error(L, "Errorr: Martin.Coulomb() number of basis funs (%u) < number of indices (%u)\n",NfunsG , Nindices);
			}
		}
		if (indToK.size() != Nindices) printf("Warning: Martin.Coulomb() number of kappas (%u) != number of indices (%u)\n",(uint32_t)indToK.size() , Nindices);

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
			printf("Error in LuaOperatorNewCoulombFull: number of functions small part (%u) != large part (%u)\n",NfunsF,NfunsG);
			return -1;
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
		CreateOperatorAtomicURelativistic(NF,wavefuns,x,indToK,indices, op);
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
			printf("Warning: Martin.Coulomb() number of basis funs (%u) != number of indices (%u)\n",NfunsG , Nindices/2);
			fflush(stdout);
		}
		else
		{
			if (NfunsG < Nindices/2)
			{
				return luaL_error(L, "Errorr: Martin.Coulomb() number of basis funs (%u) < number of indices (%u)\n",NfunsG , Nindices/2);
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

    luaL_setmetatable(L,"Operator_Type");


	for(uint32_t i = 0; i < Nindices; i++ )
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
    uint32_t NIntegrals, i;

    integrals = GetSlaterIntegrals(basisfuns,grid, indToK, NIntegrals);

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

    return 0;
}

uint32_t CreateOperatorAtomicUNonRel(uint32_t const NFermion,vector<Interpolation::InterpolatingFunction *> ipolfuns,vector<int> indToK
                                               ,uint16_t **indices, OperatorType *op)
{

    OperatorType optemp;
    SSlaterIntegral *integrals;
    uint32_t NIntegrals, i;

    integrals = GetSlaterIntegrals(ipolfuns, indToK, NIntegrals);

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

