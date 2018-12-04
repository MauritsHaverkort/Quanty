/* Some makeshift functions */

#include "MartinHeaders.h"

#include "Quadrature.cpp"
#include "Eigendifferential.cpp"
#include "MartinsLittleHelper.cpp"
#include "Old_SlaterIntegrals.cpp"
#include "CreateOperators.cpp"
#include "PlotSpectrumWithSelfEnergy.cpp"
#include "CSRMatrix.cpp"
#include "Matrix.cpp"
#include "SparseOperator.cpp"


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
		name = GetOrbitalName(file1,3);
		GetOrbitalName(file2,3);
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

static int LuaMartinTest(lua_State *L)
{
	uint32_t NfunsG = (uint32_t)lua_objlen(L,1);
	vector<Interpolation::InterpolatingFunction *> funsG = vector<Interpolation::InterpolatingFunction *>(NfunsG);
	for(uint32_t i = 0; i < NfunsG; i++)
	{
		lua_rawgeti(L, 1, i+1);
		funsG[i] = (Interpolation::InterpolatingFunction *) luaL_checkudata(L, -1, "InterpolatingFunction_Type");
		lua_pop(L,1);
	}
	vector<Interpolation::InterpolatingFunction *> funsF = vector<Interpolation::InterpolatingFunction *>(NfunsG);
	for(uint32_t i = 0; i < NfunsG; i++)
	{
		lua_rawgeti(L, 2, i+1);
		funsF[i] = (Interpolation::InterpolatingFunction *) luaL_checkudata(L, -1, "InterpolatingFunction_Type");
		lua_pop(L,1);
	}

	uint32_t k = 7;
	vector<double> x = funsG[0]->GetX();
	printf("x = %.8E\n",x[x.size()-1]);
	vector<double> knots = Quad::GetKnots(k,x);
	//PrintArray(&x[0],x.size());
	vector<double> g = (*funsG[0])(knots);
	vector<double> f = (*funsF[0])(knots);
	vector<double> w = Quad::GetWeights(k,x);
	vector<double> dg = w * g;
	vector<double> df = w * f;
	double res = 0;
	for (uint32_t i = 0; i < g.size(); i++)
	{
		res += df[i]*f[i] +dg[i]*g[i];

	}
	printf("%.8E\n",res);
	printf("%.8E\n",Quad::GetMultipoleIntegral(*(Interpolation::CubicSpline*)funsG[0],*(Interpolation::CubicSpline*)funsG[0],0));
	printf("%.8E\n",Quad::GetMultipoleIntegral(*(Interpolation::CubicSpline*)funsG[0],*(Interpolation::CubicSpline*)funsG[0],1));
	printf("%.8E\n\n",Quad::GetMultipoleIntegral(*(Interpolation::CubicSpline*)funsG[0],*(Interpolation::CubicSpline*)funsG[0],2));

	printf("%.8E\n",Quad::GetMultipoleIntegral(*funsG[0],*funsG[0],0,k));
	printf("%.8E\n",Quad::GetMultipoleIntegral(*funsG[0],*funsG[0],1,k));
	printf("%.8E\n",Quad::GetMultipoleIntegral(*funsG[0],*funsG[0],2,k));
	return 0;
}

static int LuaPlotSpectrumWithSelfEnergy(lua_State *L)
{
	double wMin, wMax,dw, gH, gS;
	std::vector<double> psi, H; // TODO make H generic: either Matrix, TriDiagonalMatrix or DiagonalMatrix
	std::vector<ResponsefunctionType *> SigmaRes;
	uint32_t len = lua_gettop(L);
	if (len != 8)
	{
		return luaL_error(L,"PlotSpectrumWithSelfEnergy needs 8 arguments not %d\n",len);
	}
	H   = luaL_checkarray_double(L,1);
	psi = luaL_checkarray_double(L,2);
	if (H.size() != psi.size())
	{
		return luaL_error(L,"PlotSpectrumWithSelfEnergy: Hamiltonian and ground state do not have the same dimensions - dim(H) = %d and dim(psi) = %d\n",H.size(),psi.size());
	}
	luaL_checktype(L, 3, LUA_TTABLE);
	uint32_t lenSigma = lua_objlen(L,3);
	SigmaRes.reserve(lenSigma);
	for(uint32_t i = 1; i <= lenSigma; i++)
	{
		lua_rawgeti(L,3,i);
		SigmaRes.push_back(luaL_checkResponsefunction(L,-1));
		lua_pop(L,1);
	}
	if (((BlockTriDiagonalMatrixType*)SigmaRes[0]->data)->BlockSize != H.size())
	{
		return luaL_error(L,"PlotSpectrumWithSelfEnergy: self energy and Hamiltonian do not have the same dimensions - dim(H) = %d and dim(Sigma) = %d\n",H.size(),((BlockTriDiagonalMatrixType*)SigmaRes[0]->data)->BlockSize);
	}
	gH   = luaL_checknumber(L,4);
	gS   = luaL_checknumber(L,5);
	wMin = luaL_checknumber(L,6);
	wMax = luaL_checknumber(L,7);
	dw   = luaL_checknumber(L,8);
	if (wMin >= wMax || dw <= 0)
	{
		return luaL_error(L,"PlotSpectrumWithSelfEnergy: cannot plot range [%f,%f] in steps of %f\n",wMin,wMax,dw);
	}

	uint32_t lenW = floor((wMax-wMin)/dw) + 1;
	std::vector<std::vector<double> > graph = std::vector<std::vector<double> >(lenW);
	double w = wMin;
	for (uint32_t i = 0; i < lenW; i++)
	{
		//w += dw;
		w = wMin + i*dw;
		graph[i] = std::vector<double>(2);
		graph[i][0] = w;
		graph[i][1] = EvaluateResponseFunction(w,gH,gS,H,SigmaRes,psi);

	}

	lua_createtable(L,graph.size(),0);
	for (uint32_t i = 0; i < graph.size(); i++)
	{
		LuaPushArray_double(L,&(graph[i][0]),2);
		lua_rawseti(L,-2,i+1);
	}

	return 1;
}

static int LuaEvaluateSelfEnergy(lua_State *L)
{
	double wMin, wMax,dw, gS;
	const char *filename;

	std::vector<ResponsefunctionType *> SigmaRes;
	uint32_t len = lua_gettop(L);
/*	if (len != 5)
	{
		return luaL_error(L,"PlotSpectrumWithSelfEnergy needs 8 arguments not %d\n",len);
	}
*/	luaL_checktype(L, 1, LUA_TTABLE);
	uint32_t lenSigma = lua_objlen(L,1);
	SigmaRes.reserve(lenSigma);
	for(uint32_t i = 1; i <= lenSigma; i++)
	{
		lua_rawgeti(L,1,i);
		SigmaRes.push_back(luaL_checkResponsefunction(L,-1));
		lua_pop(L,1);
	}
	gS   = luaL_checknumber(L,2);
	wMin = luaL_checknumber(L,3);
	wMax = luaL_checknumber(L,4);
	dw   = luaL_checknumber(L,5);
	if (wMin >= wMax || dw <= 0)
	{
		return luaL_error(L,"LuaEvaluateSelfEnergy: cannot plot range [%f,%f] in steps of %f\n",wMin,wMax,dw);
	}

	uint32_t lenW = floor((wMax-wMin)/dw) + 1;
	std::vector<CompactMatrixType *> selfenergy = std::vector<CompactMatrixType *>(lenW);
	double w = wMin;
	for (uint32_t i = 0; i < lenW; i++)
	{
		w = wMin + i*dw;
		selfenergy[i] = EvaluateSelfEnergy(w,gS,SigmaRes);
	}
	if (len > 5)
	{
		filename = luaL_checkstring(L,6);
		WriteDataToFile(selfenergy,filename);
	}
/*
	lua_createtable(L,graph.size(),0);
	for (uint32_t i = 0; i < graph.size(); i++)
	{
		LuaPushArray_double(L,&(graph[i][0]),2);
		lua_rawseti(L,-2,i+1);
	}
*/
	return 0;
}

static int LuaReadDataFromFile(lua_State *L)
{
	const char *type = luaL_checkstring(L,1);
	const char *filename = luaL_checkstring(L,2);
	if (strcmp(type, "CompactMatrixType_Table") == 0)
	{
		std::vector<CompactMatrixType *> Ms;
		ReadDataFromFile(Ms,filename);
		lua_push(L,Ms);
		return 1;
	}
	printf("Error: LuaReadDataFromFile: data type not recognized\n");

	return 0;
}

static int LuaWriteDataToFile(lua_State *L)
{
	const char *type = luaL_checkstring(L,1);
	const char *filename = luaL_checkstring(L,2);
	if (strcmp(type, "CompactMatrixType_Table") == 0)
	{
		std::vector<CompactMatrixType *> Ms;
		lua_check(L,3,Ms);
		WriteDataToFile(Ms,filename);
		printf("Table of compact matrices written to file %s\n",filename);

		return 0;
	}
	printf("Error: LuaReadDataFromFile: data type not recognized\n");

	return 0;
}

static int LuaCalcSpectrumWithSelfEnergy(lua_State *L)
{
	double wMin, wMax,dw, gH;
	std::vector<double> psi, H; // TODO make H generic: either Matrix, TriDiagonalMatrix or DiagonalMatrix
	std::vector<CompactMatrixType *> Sigma;
	uint32_t len = lua_gettop(L);

	H   = luaL_checkarray_double(L,1);
	psi = luaL_checkarray_double(L,2);
	if (H.size() != psi.size())
	{
		return luaL_error(L,"CalcSpectrumWithSelfEnergy: Hamiltonian and ground state do not have the same dimensions - dim(H) = %d and dim(psi) = %d\n",H.size(),psi.size());
	}
	lua_check(L,3,Sigma);

	gH   = luaL_checknumber(L,4);
	wMin = luaL_checknumber(L,5);
	wMax = luaL_checknumber(L,6);
	dw   = luaL_checknumber(L,7);
	if (wMin >= wMax || dw <= 0)
	{
		return luaL_error(L,"CalcSpectrumWithSelfEnergy: cannot plot range [%f,%f] in steps of %f\n",wMin,wMax,dw);
	}
	uint32_t lenW = floor((wMax-wMin)/dw) + 1;
	if (lenW != Sigma.size())
	{
		return luaL_error(L,"CalcSpectrumWithSelfEnergy: Grids for Sigma and G do not match\n");
	}
	std::vector<std::vector<double> > graph = std::vector<std::vector<double> >(lenW);
	double w = wMin;
	for (uint32_t i = 0; i < lenW; i++)
	{
		//w += dw;
		w = wMin + i*dw;
		graph[i] = std::vector<double>(2);
		graph[i][0] = w;
		graph[i][1] = EvaluateResponseFunction(w,gH,H,Sigma[i],psi);

	}
	lua_push(L,graph);
	return 1;
}

static int LuaCreateEigenDifferentials(lua_State *L)
{
	double Emin, dE;
	uint32_t NeigenDiffs,Qrelativistic;
	int kappa;
	std::vector<double> grid;
	std::vector<Interpolation::CubicSpline*> eigendifferentials;

	uint32_t len = lua_gettop(L);
	lua_check(L,1,Emin);
	lua_check(L,2,dE);
	lua_check(L,3,NeigenDiffs);
	lua_check(L,4,kappa);
	lua_check(L,5,grid);



	if (len == 6)
	{
		lua_check(L,6,Qrelativistic);
		if (Qrelativistic)
		{
			std::vector<Interpolation::CubicSpline*> eigendifferentialsG = std::vector<Interpolation::CubicSpline*>(NeigenDiffs);
			lua_createtable(L,0,NeigenDiffs);
			for (uint32_t i = 0; i < NeigenDiffs; i++)
			{
				void* memory = lua_newuserdata(L, sizeof(Interpolation::CubicSpline));
				eigendifferentialsG[i] = new(memory) Interpolation::CubicSpline();
				luaL_setmetatable(L, "InterpolatingFunction_Type");
				lua_rawseti(L,-2,i+1);

			}
			std::vector<Interpolation::CubicSpline*> eigendifferentialsF = std::vector<Interpolation::CubicSpline*>(NeigenDiffs);
			lua_createtable(L,0,NeigenDiffs);
			for (uint32_t i = 0; i < NeigenDiffs; i++)
			{
				void* memory = lua_newuserdata(L, sizeof(Interpolation::CubicSpline));
				eigendifferentialsF[i] = new(memory) Interpolation::CubicSpline();
				luaL_setmetatable(L, "InterpolatingFunction_Type");
				lua_rawseti(L,-2,i+1);

			}
			CreateEigenDiffs(grid,kappa,Emin,dE,NeigenDiffs,eigendifferentialsG,eigendifferentialsF);
			return 2;
		}
	}

	eigendifferentials = std::vector<Interpolation::CubicSpline*>(NeigenDiffs);
	lua_createtable(L,0,NeigenDiffs);
	for (uint32_t i = 0; i < NeigenDiffs; i++)
	{
		void* memory = lua_newuserdata(L, sizeof(Interpolation::CubicSpline));
		eigendifferentials[i] = new(memory) Interpolation::CubicSpline();
		luaL_setmetatable(L, "InterpolatingFunction_Type");
		lua_rawseti(L,-2,i+1);

	}
	CreateEigenDiffs(grid,kappa,Emin,dE,NeigenDiffs,eigendifferentials);

	return 1;
}

static int LuaEigenDifferential(lua_State *L)
{
	uint32_t len = lua_gettop(L);
	int kappa;
	double x,q1,q2, result;
	lua_check(L,1,x);
	lua_check(L,2,kappa);
	lua_check(L,3,q1);
	lua_check(L,4,q2);
	if (len == 5 )
	{
		struct
		{
			const char *type;
			uint32_t Nknots;
		} Options;
		LuaOptionsType LuaEDOptions[]={
		    {"Type","char*",&(Options.type)},
		    {"Nknots","uint32_t",&(Options.Nknots)},
		    {NULL, NULL, NULL}};
		Options.Nknots = 4;
		Options.type = "N";
		LuaCheckOptions(L, 5, LuaEDOptions);
		if (strcmp(Options.type,"N")==0)
		{
			result = x * EvaluateEigendifferential(&WaveFunFreeParticleNonRelativistic,x,kappa,q1,q2,Options.Nknots);
		}
		else
		{
			if (strcmp(Options.type,"G")==0)
			{
				result = x * EvaluateEigendifferential(&WaveFunFreeParticleG,x,kappa,q1,q2,Options.Nknots);
			}
			else
			{
				if (strcmp(Options.type,"F")==0)
				{
					result = x * EvaluateEigendifferential(&WaveFunFreeParticleF,x,kappa,q1,q2,Options.Nknots);
				}
				else
				{
					return luaL_error(L,"Type of Eigendifferential not recognized. Possible types are N, G, or F\n");
				}
			}
		}
	}
	else
	{
		result = x * EvaluateEigendifferential(&WaveFunFreeParticleNonRelativistic,x,kappa,q1,q2);
	}
	lua_push(L,result);
	return 1;
}

static int LuaSparseTridiagonalize(lua_State *L)
{
	OperatorType *H = (OperatorType *) luaL_checkudata(L, 1, "Operator_Type");
	std::vector<WaveFunctionType*> psi;
	WaveFunctionType *dets;
	std::vector<double> psiTransformed;
	std::vector<std::vector<double> > TEST;
	LinAlg::CSRMatrix<double> Hsparse;
	LinAlg::Matrix<double> startblock;
	uint32_t N;
	lua_check(L,2,psi);
	lua_check(L,3,dets);
	lua_check(L,4,N);
	startblock.reserve(dets->N,psi.size());
	OperatorToSparseMatrix(H,dets,Hsparse);
	//TEST = std::vector<std::vector<double> >(psi.size());
	printf("Matrix created with %lu nve\n",Hsparse.A.size());
	for(uint32_t i = 0; i < startblock.N; i++)
	{
		WaveFunctionToVector(psi[i],psiTransformed,dets);
		//WaveFunctionToVector(psi[i],TEST[i],dets);
		for (uint32_t j = 0; j < startblock.M; j++)
		{
			startblock.Entries[i*startblock.M + j] = psiTransformed[j];
		}
	}
	//lua_push(L,TEST);
	//return 1; // ----------------------------------------------------------TEST ENDE --------------
	std::vector<LinAlg::Matrix<double> > T,A,B;
	double epsSV = 100*DBL_EPSILON;
	BlockLanczos(Hsparse,startblock,N,epsSV,A,B,T);
	BlockTriDiagonalMatrixType GTri;
	GTri.BlockSize = startblock.N;
    GTri.N = N;
    GTri.QComplex=0;
    InitBlockTridiagonalMatrix(&GTri);
    GTri.mu=0;
    uint32_t bs = GTri.BlockSize;
    for(uint32_t k = 0; k < N; k++)
    {
    	for (uint32_t i = 0; i < A[k].M; i++)
    	{
    		for (uint32_t j = 0; j < A[k].N; j++)
    		{
    			GTri.aRe[(k+1)*bs*bs + i*bs + j] = A[k].Entries[j*A[k].M +i];
    		}
    	}
    	for (uint32_t i = 0; i < B[k].M; i++)
		{
			for (uint32_t j = 0; j < B[k].N; j++)
			{
				GTri.bRe[k*bs*bs + i*bs + j] = B[k].Entries[j*B[k].M +i];
			}
		}
    }

    LuaPushTable_BlockTriDiagonalMatrixType(L, &GTri);

	return 1;

}

static int LuaSparseOperator(lua_State *L)
{
	OperatorType *H = (OperatorType *) luaL_checkudata(L, 1, "Operator_Type");
	WaveFunctionType *psi = (WaveFunctionType *) luaL_checkudata(L, 2, "Wavefunction_Type");
	LinAlg::CSRMatrix<double> M;
	OperatorToSparseMatrix(H,psi,M);
	M.Print();
	//printf("Matrix created with %lu nve\n",M.A.size());
	std::vector<double> vec, a, b, y;
	WaveFunctionToVector(psi,vec,psi);
	uint32_t Mrow = vec.size();
	uint32_t Ncol = 2;
	LinAlg::Matrix<double> startblock = LinAlg::Matrix<double>(Mrow,Ncol);
	/*for(uint32_t i = 0; i < Mrow;i++)
	{
		startblock.Entries[i] = vec[i];
		startblock.Entries[i+Mrow] = vec[i]+i;
	}*/
	startblock.Entries[0] = -1;
	startblock.Entries[1] =  0;
	startblock.Entries[2] =  3;
	startblock.Entries[3] =  1;
	startblock.Entries[4] =  1;
	startblock.Entries[5] =  3;

	startblock.Entries[6] =  2;
	startblock.Entries[7] =  1;
	startblock.Entries[8] =  1;
	startblock.Entries[9] =  0;
	startblock.Entries[10] = 0;
	startblock.Entries[11] = 3;
	std::vector<LinAlg::Matrix<double> > T,A,B;
	double epsSV = 100*DBL_EPSILON;
	uint32_t N = 4;
	BlockLanczos(M,startblock,N,epsSV,A,B,T);
	printf("A:\n");
	for (int i = 0; i < N; i++)
	{
		A[i].print();
		printf("\n");
	}
	printf("B:\n");
	for (int i = 0; i < N; i++)
	{
		B[i].print();
		printf("\n");
	}
	printf("T:\n");
	for (int i = 0; i < N; i++)
	{
		T[i].print();
		printf("\n");
	}

	return 0;
/*	lua_push(L,A);
	lua_push(L,B);
	return 2;*/

/*
	uint32_t dim = 3;
	LinAlg::Matrix<double> A = LinAlg::Matrix<double>(dim,dim,false);
	A.Entries[0] = 1;
	A.Entries[1] = 2;
	A.Entries[2] = 1;
	A.Entries[3] = 2;
	A.Entries[4] = 4;
	A.Entries[5] = 2;
	A.Entries[6] = -1;
	A.Entries[7] = -2;
	A.Entries[8] = 3;
	uint32_t nSV;

	LinAlg::Matrix<double> AtA = LinAlg::Matrix<double>(dim,dim,false);
	BlockSkalarProduct(A,A,AtA);

	LoewdinSVDOrthonormalize(A,0.001,nSV);

	printf("A\n");
	A.print();

	printf("\nSVD\n");
	for (uint32_t i = 0; i < A.M; i++)
	{
		for (uint32_t j = 0; j < A.N; j++)
		{
			printf("%.4f\t",A.Entries[j*3 + i]);
		}
		printf("\n");

	}
*/
/*	printf("\nBlockSkalar\n");
	for (uint32_t i = 0; i < 3; i++)
	{
		for (uint32_t j = 0; j < 3; j++)
		{
			printf("%.4f\t",AtA.Entries[j*AtA.M + i]);
		}
		printf("\n");

	}
	AtA.print();
*/
/*


/*
	Lanczos(M,SkalarProdukt,VectorAdd<double>,VectorMul<double>,vec,6,a,b,T);
	for (uint32_t i = 0; i < a.size(); i++)
	{
		printf("%e\t",a[i]);
	}
	printf("\n");
	for (uint32_t i = 0; i < b.size(); i++)
	{
		printf("%e\t",b[i]);
	}
	printf("\n");
	size_t N = vec.size();
	for (uint32_t i = 0; i < N; i++)
	{
		vec.push_back(2*vec[i]);
	}
	for (uint32_t i = 0; i < N; i++)
	{
		vec.push_back(3*vec[i]);
	}
	M.BlockEvaluate(vec,y,3);
	for (uint32_t i = 0; i < N; i++)
	{
		printf("%.1f\t%.1f\t%.1f\n",y[i],y[i+N],y[i+2*N]);
	}

	std::vector<double> A = std::vector<double>(9);
	A[0] = 1;
	A[1] = 2;
	A[2] = 1;
	A[3] = 2;
	A[4] = 4;
	A[5] = 2;
	A[6] = -1;
	A[7] = -2;
	A[8] = 3;
	uint32_t nSV;

	std::vector<double> AtA = std::vector<double>(6);
	BlockSkalarProduct(A,A,AtA,3,2);

	LoewdinSVDOrthonormalize(A,3,3,0.001,nSV);
	printf("\nSVD\n");
	for (uint32_t i = 0; i < 3; i++)
	{
		for (uint32_t j = 0; j < 3; j++)
		{
			printf("%.4f\t",A[j*3 + i]);
		}
		printf("\n");

	}
	printf("\nBlockSkalar\n");
	for (uint32_t i = 0; i < 2; i++)
	{
		for (uint32_t j = 0; j < 2; j++)
		{
			printf("%.4f\t",AtA[j*2 + i]);
		}
		printf("\n");

	}
*/
    return 0;
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
 *;
 */

static const struct luaL_Reg Martin_func[] = {
  {"NewCoulomb", LuaOperatorNewFullCoulomb},
  {"NewDirac", LuaOperatorNewDirac},
  {"ContinuumBase",LuaCreateBaseForContinuum},
  {"ReadFPLOBasis",LuaReadFPLOBasis},
  {"RB",LuaReadFPLOBasisFunctions},
  {"Test",LuaMartinTest},
  /*{"Propagator",LuaPropagator},
  {"Propagator2",LuaPropagator2},*/
  {"Metric",LuaMetricWaveFuns},
  {"ReadData",LuaReadDataFromFile},
  {"WriteData",LuaWriteDataToFile},
  {"EvaluateSelfEnergy",LuaEvaluateSelfEnergy},
  {"PlotSpectrumWithSelfEnergy",LuaPlotSpectrumWithSelfEnergy},
  {"CalcSpectrumWithSelfEnergy",LuaCalcSpectrumWithSelfEnergy},
  {"CreateEigenDifferentials",LuaCreateEigenDifferentials},
  {"EigenDifferential",LuaEigenDifferential},
  {"ED",LuaEigenDifferential},
  {"SparseOperator",LuaSparseOperator},
  {"SparseTri",LuaSparseTridiagonalize},
  {NULL, NULL}
};

/*
 * static int LuaPropagator(lua_State *L)
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
 */
