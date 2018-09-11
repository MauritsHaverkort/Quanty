/* Some makeshift functions */

#include "MartinHeaders.h"

#include "Quadrature.cpp"
#include "Eigendifferential.cpp"
#include "MartinsLittleHelper.cpp"

#include "CreateOperators.cpp"



static int LuaOperatorNewFullCoulomb (lua_State *L)
{
    OperatorType *op;
    int i;
    uint32_t **tuple, NIntegrals, NGridpoints, Nindices, temp, len, NF;
    uint16_t **indices;
    MatrixType *wavefunctions;
    double *x;

    len = (uint32_t)lua_gettop(L);

    if ( len < 4 ||  len > 6)
    {
    	printf("Error in Martin.NewCoulomb: 4, 5 or 6 arguments expected; got %d",len);
    	return -1;
    }

    NF = (uint32_t)luaL_checkunsigned(L, 1);
    wavefunctions = luaL_checkMatrix(L,2);
    x = luaL_checkarray_double(L,3,&NGridpoints);

    if (wavefunctions->Nj != 2*NGridpoints)
    {
    	printf("Error: relativistic wavefunctions do not match grid. Relativistic wavefunctions need as many values as twice the grid points");
    	fflush(stdout);
    	return -1;
    }

    if (len == 4)
    {
    	vector<int> indToK;
    	lua_getfield(L,4,"kappas");
    	if (lua_isnil(L,-1))
    	{
    		lua_pop(L,1);
    		Nindices = (uint32_t)lua_objlen(L,4);
    		indices = (uint16_t **) malloc(Nindices * sizeof(uint16_t *));
    		if(!indices) printf("malloc fail at indices in LuaOperatorNewFullCoulomb. size = %d\n",Nindices);
    		for( i = 0; i < Nindices; i++ )
    		{
    			lua_rawgeti(L, 4, i+1);
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
    		//PrintArray(indToK);
    		lua_pop(L,1);
    		Nindices = (uint32_t)lua_objlen(L,4);

    		indices = (uint16_t **) malloc(Nindices * sizeof(uint16_t *));
    		if(!indices) printf("malloc fail at indices in LuaOperatorNewFullCoulomb. size = %d\n",Nindices);

    		for( i = 0; i < Nindices; i++ )
    		{
    			lua_rawgeti(L, 4, i+1);
    			indices[i] = luaL_checkarray_unit16_t(L,-1,&temp);
    			lua_pop(L,1);
    		}
    	}
    	if (wavefunctions->Ni > Nindices)
		{
    		printf("Warning: Martin.Coulomb() number of basis funs (%d) != number of indices (%d)\n",wavefunctions->Ni , Nindices);
    		fflush(stdout);
		}
    	else
    	{
    		if (wavefunctions->Ni < Nindices)
    		{
    			return luaL_error(L, "Errorr: Martin.Coulomb() number of basis funs (%d) < number of indices (%d)\n",wavefunctions->Ni , Nindices);
    		}
    	}

    	if (indToK.size() != Nindices) printf("Warning: Martin.Coulomb() number of kappas (%d) != number of indices (%d)\n",(int) indToK.size() , (int) Nindices);
    	op = (OperatorType *) lua_newuserdata(L,sizeof(OperatorType));
    	CreateOperatorFastCoulombRelativistic(NF,wavefunctions->VRe,wavefunctions->Ni,x,NGridpoints,indToK,indices, op);
    }
    else
    {
    	NIntegrals = (uint32_t)lua_objlen(L,4);
    	tuple = (uint32_t **) malloc(NIntegrals*sizeof(uint32_t *));
    	if(!tuple) printf("malloc fail at tuple in LuaOperatorNewFullCoulomb. size = %d\n",NIntegrals);
    	for ( i = 0; i < NIntegrals; i++)
    	{
    		lua_rawgeti(L, 4, i+1);
    		tuple[i] = luaL_checkarray_unit32_t(L,-1,&temp);
    		if ( temp != 8)
    		{
    			printf("Error in Martin.NewCoulomb: Length of tuple[%d] (%d) != 8",i,temp);
    			return -1;
    		}
    		lua_pop(L,1);
    	}

    	luaL_checktype(L, 5, LUA_TTABLE);
    	Nindices = (uint32_t)lua_objlen(L,5);

    	indices = (uint16_t **) malloc(Nindices * sizeof(uint16_t *));
    	if(!indices) printf("malloc fail at indices in LuaOperatorNewFullCoulomb. size = %d\n",Nindices);

    	if (wavefunctions->Ni > Nindices)
		{
    		printf("Warning: Martin.Coulomb() number of basis funs (%d) != number of indices (%d)\n",wavefunctions->Ni , Nindices);
    		fflush(stdout);
		}
    	else
    	{
    		if (wavefunctions->Ni < Nindices)
    		{
    			return luaL_error(L, "Errorr: Martin.Coulomb() number of basis funs (%d) < number of indices (%d)\n",wavefunctions->Ni , Nindices);
    		}
    	}

    	for( i = 0; i < Nindices; i++ )
    	{
    		lua_rawgeti(L, 5, i+1);
    		indices[i] = luaL_checkarray_unit16_t(L,-1,&temp);
    		lua_pop(L,1);
    	}



    	op = (OperatorType *) lua_newuserdata(L,sizeof(OperatorType));


    	if ( len == 5)
    	{
    		CreateOperatorFullCoulombRelativistic(NF,wavefunctions->VRe,wavefunctions->Ni,x,NGridpoints,tuple, NIntegrals,indices, op);
    	}
    	else
    	{
    		uint32_t test = (uint32_t)luaL_checkunsigned(L, 6);
    		if (test == 2)
    		{
    			CreateOperatorSingleCoulombRelativistic(NF,wavefunctions->VRe,wavefunctions->Ni,x,NGridpoints,tuple, NIntegrals,indices, op);
    		}
    		else
    		{
    			CreateOperatorFastCoulombRelativistic(NF,wavefunctions->VRe,wavefunctions->Ni,x,NGridpoints,tuple, NIntegrals,indices, op);
    		}
    	}
    	for( i = 0; i < NIntegrals; i++ )
    	{
    		free(tuple[i]);
    	}
    	free(tuple);
    }

    luaL_setmetatable(L,"Operator_Type");

	FreeMatrix(wavefunctions);
	free(x);

	for( i = 0; i < Nindices; i++ )
	{
		free(indices[i]);
	}
	free(indices);

    return 1;
}

static int LuaOperatorNewDirac(lua_State *L)
{
	OperatorType *op;
	uint32_t **indices;
	uint32_t NF, i, NGroups, NGridpoints, NIndices, temp;
	SKappaGroup *groups;
	MatrixType *basisfuns, *derivfuns;
	int Z;
	double *x;

	NF = (uint32_t)luaL_checkunsigned(L, 1);
	basisfuns = luaL_checkMatrix(L,2);
	derivfuns = luaL_checkMatrix(L,3);
	x = luaL_checkarray_double(L,4,&NGridpoints);
	luaL_checktype(L, 5, LUA_TTABLE);
	NGroups = (uint32_t)lua_objlen(L,5);

	groups = (SKappaGroup *) malloc(NGroups * sizeof(SKappaGroup));
	if (!groups) printf("malloc fail at groups in LuaOperatorNewDirac. size = %d\n",NGroups);

	for( i = 0; i < NGroups; i++ )
	{
		lua_rawgeti(L, 5, i+1);
	   	lua_rawgeti(L, -1, 1);
	   	groups[i].kappa = (int) luaL_checkunsigned(L, -1);
	   	lua_pop(L,1);
	   	lua_rawgeti(L, -1, 2);
		groups[i].indices = luaL_checkarray_unit32_t(L,-1,&(groups[i].n));
	   	lua_pop(L,1);
	   	lua_pop(L,1);
	}
	Z = (int) luaL_checkunsigned(L, 6);
	luaL_checktype(L, 7, LUA_TTABLE);
	NIndices = (uint32_t)lua_objlen(L,7);

	indices = (uint32_t **) malloc(NIndices * sizeof(uint32_t *));
	if(!indices) printf("malloc fail at indices in LuaOperatorNewDirac. size = %d\n",NIndices);
	for( i = 0; i < NIndices; i++ )
	{
	   	lua_rawgeti(L, 7, i+1);
	   	indices[i] = luaL_checkarray_unit32_t(L,-1,&temp);
	   	lua_pop(L,1);
	}

	op = (OperatorType *) lua_newuserdata(L,sizeof(OperatorType));
	CreateOperatorDirac(NF,basisfuns->Ni, basisfuns->VRe,derivfuns->VRe,x,NGridpoints,groups,NGroups,indices,Z,op);

	luaL_setmetatable(L,"Operator_Type");

	FreeMatrix(basisfuns);
	FreeMatrix(derivfuns);
	free(x);
	free(groups);
	for( i = 0; i < NIndices; i++ )
	{
		free(indices[i]);
	}
	free(indices);

	return 1;
}

static int LuaCreateBaseForContinuum(lua_State *L)
{
	uint32_t NGridpoints = 0;
	SKappaGroup group;
	MatrixType *basisfuns, *derivfuns, *contfuns;
	double w;
	//double *x;
	vector<double> x;
	vector<vector<double> > vbasf, vderf,vcontf;

	uint32_t len = (uint32_t)lua_gettop(L);
	if ( len < 5 )
	{
		std::cout << "Error in LuaCreateBaseForContinuum: Wrong number of arguments. Expected 5, 6 or 7, but got " << len << std::endl;
		return 0;
	}

	basisfuns = luaL_checkMatrix(L,1);
	vbasf = ConvertMatrix(basisfuns->VRe,basisfuns->Ni,basisfuns->Nj);

	derivfuns = luaL_checkMatrix(L,2);
	vderf = ConvertMatrix(derivfuns->VRe,derivfuns->Ni,derivfuns->Nj);

	x = luaL_checkarray_double(L,3);
	NGridpoints = (uint32_t) x.size();

	luaL_checktype(L, 4, LUA_TTABLE);
	lua_rawgeti(L, 4, 1);
	group.kappa = (int) luaL_checkunsigned(L, -1);
	lua_pop(L,1);
	lua_rawgeti(L, 4, 2);
	group.indices = luaL_checkarray_unit32_t(L,-1,&(group.n));
	lua_pop(L,1);

	double *eigenval;
	vector<vector<double> > energiesCont;

	if (len == 5)
	{

		contfuns = luaL_checkMatrix(L,5);
		vcontf = ConvertMatrix(contfuns->VRe,contfuns->Ni,contfuns->Nj);

		double dE = 1/2./13.6;
		double Emin = dE;
		energiesCont = CreateEigenDifferentialEnergies((uint32_t) vcontf.size(),Emin,dE);

		eigenval = CreateBaseForContinuum(vbasf,vderf,x,group,vcontf,energiesCont);

		LuaPushArray_double(L,&(eigenval[0]), (int) vcontf.size());
		for (int i = 0; i < vcontf.size(); i++)
		{
			contfuns->VRe[i] = &(vcontf[i][0]);
		}
		contfuns->Ni = (uint32_t) vcontf.size();
		LuaPushTable_MatrixType(L,contfuns);
	}
	if (len == 6)
	{
		w = luaL_checknumber(L, 5);
		luaL_checktype(L, 6, LUA_TTABLE);

		contfuns = luaL_checkMatrix(L,6);
		vcontf = ConvertMatrix(contfuns->VRe,contfuns->Ni,contfuns->Nj);

		eigenval = CreateBaseForContinuum((uint32_t) vbasf.size(), vbasf,vderf,x,NGridpoints,group, (uint32_t) vcontf.size(),vcontf,w);

		LuaPushArray_double(L,&(eigenval[vcontf.size()/2]), (int) vcontf.size()/2);
		for (uint32_t i = 0; i < vcontf.size()/2; i++)
		{
			contfuns->VRe[i] = &(vcontf[vcontf.size()/2+i][0]);
		}
		contfuns->Ni = (uint32_t) vcontf.size()/2;
		LuaPushTable_MatrixType(L,contfuns);
	}
	if (len == 7)
	{
		double Emin = luaL_checknumber(L, 5);
		double dE   = luaL_checknumber(L, 6);
		int NeigenDiffs = (int) luaL_checkunsigned(L, 7);

		contfuns = (MatrixType*)malloc(sizeof(MatrixType));
		contfuns->QComplex=0;
		contfuns->a0=0;
		contfuns->b0=0;
		contfuns->mu=0;
		vcontf = CreateEigenDiffsOnGrid(x,group.kappa,Emin,dE,NeigenDiffs);
		energiesCont = CreateEigenDifferentialEnergies((uint32_t) vcontf.size(),Emin,dE);

		eigenval = CreateBaseForContinuum(vbasf,vderf,x,group,vcontf,energiesCont);

		contfuns->Ni=(uint32_t) vcontf.size();
		contfuns->Nj=(uint32_t) vcontf[0].size();
		contfuns->VRe = (double**) malloc(contfuns->Ni * sizeof(double*));

		LuaPushArray_double(L,&(eigenval[0]), (int) vcontf.size());
		for (uint32_t i = 0; i < vcontf.size(); i++)
		{
			contfuns->VRe[i] = &(vcontf[i][0]);
		}
		LuaPushTable_MatrixType(L,contfuns);


	}


	FreeMatrix(basisfuns);
	FreeMatrix(derivfuns);

	//free(x);


	return 2;
}

static int LuaReadFPLOBasis(lua_State *L)
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
	std::map<string,vector<double> > funs;

	while (!file1.eof())
	{
		name = GetOrbitalName(file1);
		GetOrbitalName(file2);
		y1 = ReadBasisfunction(file1,x);
		y2 = ReadBasisfunction(file2,x);
		VectorJoin(y1,y2);
		funs[name] = y1;
	}
	file1.close();
	file2.close();

	MatrixType *basisfuns = (MatrixType*)malloc(sizeof(MatrixType));
	basisfuns->Ni = Norbitals;
	basisfuns->QComplex = 0;
	basisfuns->VRe = (double **) malloc(sizeof(double*)*Norbitals);
	for (int i = 0; i < Norbitals; i++)
	{
		string charTostring(orbitalNames[i]);
		if (funs.find(charTostring) != funs.end())
		{
			basisfuns->Nj = (uint32_t) funs[charTostring].size();
			basisfuns->VRe[i] = (double *) malloc(sizeof(double)*basisfuns->Nj);
			for (int j = 0; j < basisfuns->Nj; j++)
			{
				basisfuns->VRe[i][j] = funs[charTostring][j];
			}
		}
		else
		{
			std::cout << "LuaReadFPLOBasis: Couldn't find " << orbitalNames[i] << std::endl;
		}

	}

	LuaPushTable_MatrixType(L,basisfuns);
	LuaPushArray_double(L,&x[0],(int) x.size());

	return 2;
}

static int LuaMetricWaveFuns(lua_State *L)
{
	uint32_t NGridpoints = 0;

	MatrixType *basisfuns, *funs2;

	double *x;
	vector<double> vx, metric;
	vector<vector<double> > vbasf, vfuns2;

	basisfuns = luaL_checkMatrix(L,1);
	vbasf = ConvertMatrix(basisfuns->VRe,basisfuns->Ni,basisfuns->Nj);

	funs2 = luaL_checkMatrix(L,2);
	vfuns2 = ConvertMatrix(funs2->VRe,funs2->Ni,funs2->Nj);
	x = luaL_checkarray_double(L,3,&NGridpoints);
	vx = ConvertArray(x,NGridpoints);

	metric = CalcMetric(vbasf,vfuns2,vx);

	LuaPushArray_double(L,&(metric[0]),(int) metric.size());

	FreeMatrix(basisfuns);
	FreeMatrix(funs2);

	free(x);

	return 1;
}

static int LuaPropagator(lua_State *L)
{
	uint32_t n0, n1, nU, nE, nv;
	vector<CompactMatrixType *> U;
	vector<double> w, spec;

	Complex16 *v;
	double  *d0, *d1, *Econt;

	d0 = luaL_checkarray_double(L,1,&n0);
	d1 = luaL_checkarray_double(L,2,&n1);
	nU = (uint32_t)lua_objlen(L,3);
	for (int i = 0; i < nU; i++)
	{
	   	lua_rawgeti(L, 3, i+1);
	   	U.insert(U.end(),(CompactMatrixType *) luaL_checkCompactMatrix(L, -1));
	   	lua_pop(L,1);
	}

	Econt = luaL_checkarray_double(L,4,&nE);
	if (nE != nU)
	{
		printf("errror in LuaPropagator: %d = nE != nU = %d\n",nE,nU);
		return -1;
	}
	v = luaL_checkarray_complex16(L,5,&nv);
	if (nv != n0)
	{
		printf("error in LuaPropagator: %d = nv != n0 = %d\n",nv,n0);
		return -1;
	}
	double Q = luaL_checknumber(L, 6);
	double ww = 0;
	double g0 = luaL_checknumber(L, 7);
	double g1 = luaL_checknumber(L, 8);
	while (ww < Q)
	{
		w.insert(w.end(),ww);
		ww += g1/8.; //TODO adjust
	}
	spec = CalcSpectrumContinuum(w,g0,g1,Q,d0,Econt,d1,U,v);
	LuaPushArray_double(L,&w[0],(int) w.size());
	LuaPushArray_double(L,&spec[0],(int) spec.size());

	free(d0);
	free(d1);
	free(Econt);
	free(v);
	for (int i = 0; i < U.size(); i++)
	{
		FreeCompactMatrix(U[i]);
	}
	return 2;
}

static int LuaPropagator2(lua_State *L)
{
	uint32_t n0, n1, nd1, nU, nE, nv;
	vector<CompactMatrixType *> U;
	vector<double> w, spec;

	Complex16 *v;
	double  *d0, *Econt;
	vector<double *> d1;

	d0 = luaL_checkarray_double(L,1,&n0);
	nd1 = (uint32_t)lua_objlen(L,2);
	d1 = vector<double *>(nd1);
	for (int i = 0; i < nd1; i++)
	{
		lua_rawgeti(L, 2, i+1);
		d1[i] = luaL_checkarray_double(L,-1,&n1);
		lua_pop(L,1);
	}
	nU = (uint32_t)lua_objlen(L,3);
	for (int i = 0; i < nU; i++)
	{
	   	lua_rawgeti(L, 3, i+1);
	   	U.insert(U.end(),(CompactMatrixType *) luaL_checkCompactMatrix(L, -1));
	   	lua_pop(L,1);
	}

	Econt = luaL_checkarray_double(L,4,&nE);
	if (nE != nU)
	{
		printf("errror in LuaPropagator: %d = nE != nU = %d\n",nE,nU);
		return -1;
	}
	v = luaL_checkarray_complex16(L,5,&nv);
	if (nv != n0)
	{
		printf("error in LuaPropagator: %d = nv != n0 = %d\n",nv,n0);
		return -1;
	}
	double Q = luaL_checknumber(L, 6);
	double ww = 0;
	double g0 = luaL_checknumber(L, 7);
	double g1 = luaL_checknumber(L, 8);
	while (ww < Q)
	{
		w.insert(w.end(),ww);
		ww += g1/8.; //TODO adjust
	}
	spec = CalcSpectrumContinuum(w,g0,g1,Q,d0,Econt,d1,U,v);
	LuaPushArray_double(L,&w[0],(int) w.size());
	LuaPushArray_double(L,&spec[0],(int) spec.size());

	free(d0);
	for (int i = 0; i < nd1; i++)
	{
		free(d1[i]);
	}
	free(Econt);
	free(v);
	for (int i = 0; i < U.size(); i++)
	{
		FreeCompactMatrix(U[i]);
	}
	return 2;
}

static const struct luaL_Reg Martin_meta[] = {
  {NULL, NULL}
};

/*
 * Martin.NewCoulomb(NF, basisfuns,x,tuple,indices,...)
 * 	creates all Coulomb terms specified by tupel and calculates
 * 	the needed Slater integrals for the basisfuns on grid x,
 *
 * 		NF: number of fermionic modes
 * 		x: list of points on which the radial wave functions are defined
 *
 * 		basisfuns: list of radial wave functions on x (let N be the number of points)
 * 		basisfuns = { {G1(x0),...,G1(xN),F1(x0),...,F1(xN)},
 * 					  {G2(x0),...,G2(xN),F2(x0),...,F2(xN)}, ...}
 * 				where G (F) stands for the large (small) spinor component
 *
 * 		indices: indices of the mj orbitals corresponding to the basisfuns
 * 		indices = { {m1(m=-j),...,m1(m=+j)}, {m2(m=-j),...,m2(m=+j)}, ... }
 * 			where indices[i] corresponds to basisfuns[i]
 * 		tuple: list of 4-tuple of quantum numbers for slater integrals to be used
 * 		tuple = { {n1,k1,n2,k2,n3,k3,n4,k4}, ... }
 *
 * 		if a 6th argument "fast" is specified, tuple MUST NOT contain duplicates in
 * 		the sense that two 4-tuple are equivalent if you can map one onto the other by
 * 		using any combination of the following permutations: { (1,3), (2,4), (1,2)(3,4) }
 * 		this assures that Slater integrals with the same value are calculated once only,
 * 		which leads to an increase in speed.
 *
 *
 *
 */

static const struct luaL_Reg Martin_func[] = {
  {"NewCoulomb", LuaOperatorNewFullCoulomb},
  {"NewDirac", LuaOperatorNewDirac},
  {"ContinuumBase",LuaCreateBaseForContinuum},
  {"ReadFPLOBasis",LuaReadFPLOBasis},
  {"RB",LuaReadFPLOBasisFunctions},
  {"Propagator",LuaPropagator},
  {"Propagator2",LuaPropagator2},
  {"Metric",LuaMetricWaveFuns},
  {NULL, NULL}
};
