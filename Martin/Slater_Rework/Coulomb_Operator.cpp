/*
 * Coulomb_Operator.cpp
 *
 *  Created on: Jul 21, 2020
 *      Author: brass
 */

#define sgn(x) ((x > 0) - (x < 0))

void Create_Relativistic_Coulomb_Operator(const vector<Interpolation::InterpolatingFunction*> &G,
										  const vector<Interpolation::InterpolatingFunction*> &F,
										  const vector<vector<uint16_t> > &indices,
										  const vector<int> &list_of_kappas,
										  const uint32_t NF,
										  OperatorType &U)
{
	InitNamedOperator(&U, "Coulomb Operator");
	U.NFermion = NF;
	U.NWord = Min(U.NFermion, 1 + (U.NFermion - 1) / 8);
	InitOperatorFixedLength(&U, 2, 8);
	InitOperatorFixedLength(&U, 4, 8);
	vector<double> powx;
	Init_Slater_Kernel(list_of_kappas,G[0]->x,powx);
	uint32_t N = (uint32_t) G.size();
#pragma omp parallel for
	for (uint32_t i1 = 0; i1 < N; i1++)
	{
		int k1 = list_of_kappas[i1];
		uint32_t i2 = i1;
		int k2 = list_of_kappas[i2];
		for (uint32_t i3 = i1; i3 < N; i3++)
		{
			int k3 = list_of_kappas[i3];
			for (uint32_t i4 = i3; i4 < N; i4++)
			{
				int k4 = list_of_kappas[i4];
				int L0 = Max( abs( abs(k1)-abs(k3) ), abs( abs(k2)-abs(k4) ) );
				int L1 = Min( abs(k1)+abs(k3), abs(k2)+abs(k4) ) -1;
				if ( (abs(k1)+abs(k3) + (sgn(k1)+sgn(k3))/2 + L0) % 2 == 0)
				{
					L0++;
				}
				if (L0 <= L1)
				{
					double Fk[(L1-L0)/2 +1];
					for (int L = L0; L <= L1; L +=2)
					{
						Fk[(L - L0)/2] = Slater_Integrate(G[0]->x,powx,L,G[i1]->y,G[i2]->y,G[i3]->y,G[i4]->y,F[i1]->y,F[i2]->y,F[i3]->y,F[i4]->y);
					}
					OperatorType optemp;
					CreateOperatorAtomicURelativisticOneIntegral(NF,k1,k2,k3,k4,&(indices[i1][0]),&(indices[i2][0]),&(indices[i3][0]),&(indices[i4][0]),Fk,i1,i2,i3,i4,&optemp);
#pragma omp critical
					OperatorAdd(&U, optemp, 1., 0.);
					FreeOperator(&optemp);
				}
			}
		}

	}
#pragma omp parallel for
	for (uint32_t i1 = 0; i1 < N; i1++)
	{
		int k1 = list_of_kappas[i1];
		for (uint32_t i2 = i1+1; i2 < N; i2++)
		{
			int k2 = list_of_kappas[i2];
			for (uint32_t i3 = i1; i3 < N; i3++)
			{
				int k3 = list_of_kappas[i3];
				for (uint32_t i4 = i2; i4 < N; i4++)
				{
					int k4 = list_of_kappas[i4];
				    int L0 = Max( abs( abs(k1)-abs(k3) ), abs( abs(k2)-abs(k4) ) );
				    int L1 = Min( abs(k1)+abs(k3), abs(k2)+abs(k4) ) -1;
				    if ( (abs(k1)+abs(k3) + (Sign(k1)+Sign(k3))/2 + L0) % 2 == 0)
				    {
				        L0++;
				    }
				    if (L0 <= L1)
				    {
				    	double Fk[(L1-L0)/2 +1];
				    	for (int L = L0; L <= L1; L +=2)
				    	{
				    		Fk[(L - L0)/2] = Slater_Integrate(G[0]->x,powx,L,G[i1]->y,G[i2]->y,G[i3]->y,G[i4]->y,F[i1]->y,F[i2]->y,F[i3]->y,F[i4]->y);
				    	}
				    	OperatorType optemp;
				    	CreateOperatorAtomicURelativisticOneIntegral(NF,k1,k2,k3,k4,&(indices[i1][0]),&(indices[i2][0]),&(indices[i3][0]),&(indices[i4][0]),Fk,i1,i2,i3,i4,&optemp);
#pragma omp critical
				    	OperatorAdd(&U, optemp, 1., 0.);
						FreeOperator(&optemp);
				    }
				}
			}
		}
	}
}

static int LuaCreate_Relativistic_Coulomb_Operator (lua_State *L)
{
	uint32_t NF;
	vector<vector<uint16_t> > indices;
	vector<int> list_of_kappas;
	vector<Interpolation::InterpolatingFunction*> G,F;
	lua_check(L,1,NF);
	lua_check(L,2,indices);
	lua_getfield(L,2,"kappas");
	if (lua_isnil(L,-1))
	{
		lua_pop(L,1);
		return luaL_error(L,"Error in LuaCreate_Relativistic_Coulomb_Operator: Argument 2: indices.kappas does not exist.\n");
	}
	else
	{
		lua_check(L,-1,list_of_kappas);
		lua_pop(L,1);
		if (indices.size() != list_of_kappas.size()) return luaL_error(L,"Error in LuaCreate_Relativistic_Coulomb_Operator: Argument 2 : indices.kappas and indices do not have same size.\n");
	}
	lua_check(L,3,G);
	lua_check(L,4,F);
	if (G.size() != F.size()) return luaL_error(L,"Error in LuaCreate_Relativistic_Coulomb_Operator: Argument 3 & 4 (1-Particle wave-functions) do not have same size.\n");
	if (G.size() != indices.size()) return luaL_error(L,"Error in LuaCreate_Relativistic_Coulomb_Operator: Argument 2 & 3 do not have same size.\n");

	OperatorType* U = (OperatorType *) lua_newuserdata(L,sizeof(OperatorType));
	luaL_setmetatable(L,"Operator_Type");
	Create_Relativistic_Coulomb_Operator(G,F,indices,list_of_kappas,NF,*U);
	return 1;
}

