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
#include "CouplingToEMField.cpp"
#include "AuxFunctions.cpp"
#include "ResonantSpectraSelfEnergy.cpp"
#include "Continuum-Diagonal-Method.cpp"
#include "Math/Quad.cpp"
#include "Math/BSpline.cpp"
#include "ReduceOperator.cpp"
#include "Auger-Meitner_Selfenergy/Modules.cpp"
#include "Auger-Meitner_Selfenergy/ResonantSpectrum_Selfenergy.cpp"


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
    			double cutoff;
    			lua_check(L,7,cutoff);
    			CreateOperatorSingleCoulombRelativistic(NF,wavefunctions->VRe,wavefunctions->Ni,x,NGridpoints,tuple, NIntegrals,indices, op,cutoff);
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
	uint32_t lenSigma = (uint32_t)lua_objlen(L,3);
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

	lua_createtable(L,(int)graph.size(),0);
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
	uint32_t lenSigma = (uint32_t) lua_objlen(L,1);
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
	if (len == 7)
	{
		std::vector<double> energyshift;
		lua_check(L,7,energyshift);
		if (energyshift.size() != SigmaRes.size())
		{
			return luaL_error(L,"LuaEvaluateSelfEnergy: %d = energyshift.size() != SigmaRes.size() = %d \n",energyshift.size() , SigmaRes.size());
		}
		for (uint32_t i = 0; i < lenW; i++)
		{
			w = wMin + i*dw;
			selfenergy[i] = EvaluateSelfEnergy(w,gS,SigmaRes,energyshift);
		}
	}
	else
	{
		for (uint32_t i = 0; i < lenW; i++)
		{
			w = wMin + i*dw;
			selfenergy[i] = EvaluateSelfEnergy(w,gS,SigmaRes);
		}
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
	if (strcmp(type, "CSRMatrix_Type") == 0)
	{
		LinAlg::CSRMatrix<double> *M = (LinAlg::CSRMatrix<double> *)lua_newuserdata(L, sizeof(LinAlg::CSRMatrix<double>));
		luaL_setmetatable(L,"CSRMatrix_Type");
		ReadDataFromFile(M,filename);
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
	if (strcmp(type, "CSRMatrix_Type") == 0)
	{
		LinAlg::CSRMatrix<double> *M = (LinAlg::CSRMatrix<double> *) luaL_checkudata(L, 3, "CSRMatrix_Type");
		WriteDataToFile(M,filename);
		printf("CSRMatrix written to file %s\n",filename);
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
	//uint32_t len = lua_gettop(L);

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
	std::vector<WaveFunctionType> Opsi;
	WaveFunctionType *dets;
	std::vector<double> psiTransformed;
	LinAlg::CSRMatrix<double> Hsparse;
	LinAlg::Matrix<double> startblock;
	uint32_t N;
	lua_check(L,2,psi);
	lua_check(L,5,dets);
	lua_check(L,4,N);
	OperatorType *OppT = (OperatorType *) luaL_checkudata(L, 3, "Operator_Type");
	Opsi = std::vector<WaveFunctionType>(psi.size());
	startblock.reserve(dets->N,(uint32_t) psi.size());
	OperatorToSparseMatrix(H,dets,Hsparse);
	printf("Matrix created with %lu nve\n",Hsparse.A.size());
	for(uint32_t i = 0; i < startblock.N; i++)
	{
		OperatorPsi(*OppT,*(psi[i]),&(Opsi[i]));
		WaveFunctionToVector(&(Opsi[i]),psiTransformed,dets);
		for (uint32_t j = 0; j < startblock.M; j++)
		{
			startblock.Entries[i*startblock.M + j] = psiTransformed[j];
		}
	}
	for (uint32_t i = 0; i < startblock.N; i++)
	{
		WaveFunctionFree(&(Opsi[i]));
	}
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
    for (uint32_t i = 0; i < bs; i++)
	{
		GTri.bRe[i*bs + i] = 0;
	}
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
    // b0 correspond to B[0]^T since I am stupid and did not use the same convention...
	for (uint32_t i = 0; i < bs; i++)
	{
		for (uint32_t j = i+1; j < bs; j++)
		{
			double help = GTri.bRe[i*bs + j];
			GTri.bRe[i*bs + j] = GTri.bRe[j*bs + i];
			GTri.bRe[j*bs + i] = help;
		}
	}

    LuaPushTable_BlockTriDiagonalMatrixType(L, &GTri);

	return 1;

}

static int LuaCSRTridiagonalize(lua_State *L)
{
	LinAlg::CSRMatrix<double> *Hsparse = (LinAlg::CSRMatrix<double> *) luaL_checkudata(L, 1, "CSRMatrix_Type");
	std::vector<WaveFunctionType*> psi;
	std::vector<WaveFunctionType> Opsi;
	WaveFunctionType *dets;
	std::vector<double> psiTransformed;
	LinAlg::Matrix<double> startblock;
	uint32_t N;
	lua_check(L,2,psi);
	lua_check(L,5,dets);
	lua_check(L,4,N);
	OperatorType *OppT = (OperatorType *) luaL_checkudata(L, 3, "Operator_Type");
	Opsi = std::vector<WaveFunctionType>(psi.size());
	startblock.reserve(dets->N,(uint32_t)psi.size());
	printf("Matrix created with %lu nve\n",Hsparse->A.size());
	for(uint32_t i = 0; i < startblock.N; i++)
	{
		OperatorPsi(*OppT,*(psi[i]),&(Opsi[i]));
		WaveFunctionToVector(&(Opsi[i]),psiTransformed,dets);
		for (uint32_t j = 0; j < startblock.M; j++)
		{
			startblock.Entries[i*startblock.M + j] = psiTransformed[j];
		}
	}
	for (uint32_t i = 0; i < startblock.N; i++)
	{
		WaveFunctionFree(&(Opsi[i]));
	}
	std::vector<LinAlg::Matrix<double> > T,A,B;
	double epsSV = 100*DBL_EPSILON;
	BlockLanczos(*Hsparse,startblock,N,epsSV,A,B,T);
	BlockTriDiagonalMatrixType GTri;
	GTri.BlockSize = startblock.N;
    GTri.N = N;
    GTri.QComplex=0;
    InitBlockTridiagonalMatrix(&GTri);
    GTri.mu=0;
    uint32_t bs = GTri.BlockSize;
    for (uint32_t i = 0; i < bs; i++)
	{
    	GTri.bRe[i*bs + i] = 0;
	}
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
    // b0 correspond to B[0]^T since I am stupid and did not use the same convention...
    for (uint32_t i = 0; i < bs; i++)
	{
    	for (uint32_t j = i+1; j < bs; j++)
    	{
    		double help = GTri.bRe[i*bs + j];
    		GTri.bRe[i*bs + j] = GTri.bRe[j*bs + i];
    		GTri.bRe[j*bs + i] = help;
    	}
	}

    LuaPushTable_BlockTriDiagonalMatrixType(L, &GTri);

	return 1;

}

static int LuaSparseOperator(lua_State *L)
{
	OperatorType *H = (OperatorType *) luaL_checkudata(L, 1, "Operator_Type");
	WaveFunctionType *psi = (WaveFunctionType *) luaL_checkudata(L, 2, "Wavefunction_Type");
	void* memory = lua_newuserdata(L, sizeof(LinAlg::CSRMatrix<double>));
	LinAlg::CSRMatrix<double> *M = new(memory) LinAlg::CSRMatrix<double>;
	luaL_setmetatable(L,"CSRMatrix_Type");
	OperatorToSparseMatrix(H,psi,*M);
	M->Print();
	//printf("Matrix created with %lu nve\n",M.A.size());
	std::vector<double> vec, a, b, y;
	WaveFunctionToVector(psi,vec,psi);
	uint32_t Mrow = (uint32_t) vec.size();
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
	BlockLanczos(*M,startblock,N,epsSV,A,B,T);
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

	return 1;
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

static int LuaOperatorToCSRMatrix(lua_State *L)
{
	OperatorType *H = (OperatorType *) luaL_checkudata(L, 1, "Operator_Type");
	WaveFunctionType *dets;
	lua_check(L,2,dets);

	void* memory = lua_newuserdata(L, sizeof(LinAlg::CSRMatrix<double>));
	LinAlg::CSRMatrix<double> *Hsparse = new(memory) LinAlg::CSRMatrix<double>;
	luaL_setmetatable(L,"CSRMatrix_Type");
	OperatorToSparseMatrix(H,dets,*Hsparse);
	return 1;
}

static int LuaExport(lua_State *L)
{
	std::vector<std::vector<double> > table;
	std::vector<double> list;
	const char *filename = luaL_checkstring(L,2);
	FILE *f = fopen(filename,"w");
	if (luaL_isType(L,1,table))
	{
		lua_check(L,1,table);
		for (size_t i = 0; i < table.size(); i++)
		{
			for (size_t j = 0; j < table[i].size(); j++)
			{
				fprintf(f,"%E\t",table[i][j]);
			}
			fprintf(f,"\n");
		}
	}
	else
	{
		lua_check(L,1,list);
		for (size_t i = 0; i < list.size(); i++)
		{
			fprintf(f,"%E\n",list[i]);
		}
	}
	return 0;
}

static int LuaReadFPLOEnergies(lua_State *L)
{
	string filename; 	//  = "out.scf";
	string orb; 			// = "4s"
	string j;   			// = "1/2"

	vector<string> orbitalLabels;
	vector<double> energies;
	lua_check(L,1,orbitalLabels);
	lua_check(L,2,filename);
	energies.reserve(orbitalLabels.size());
	for (size_t i = 0; i < orbitalLabels.size(); i++)
	{
		orb = orbitalLabels[i].substr(0,2);
		j   = orbitalLabels[i].substr(2,3);
		energies.push_back(ParseFPLOEnergies(filename.c_str(),orb,j));
	}
	lua_push(L,energies);
	return 1;
}

static int LuaRelativisticDipoleOperator(lua_State *L)
{
	vector<Interpolation::InterpolatingFunction* > G, F;
	vector<double> energies;
	vector<int> indToK;
	vector<vector<uint32_t> > indices;
	uint32_t mu, NF;
	lua_check(L,1,G);
	lua_check(L,2,F);
	lua_check(L,3,energies);
	lua_check(L,4,mu);
	lua_getfield(L,5,"kappas");
	if (lua_isnil(L,-1))
	{
		return luaL_error(L,"Error in LuaRelativisticDipoleOperator: indices do not contain field \"kappas\"\n Please use the function CreateAtomicIndicesList to generate indices with option {{\"Kappas\",true}}\n");;
	}
	lua_check(L,-1,indToK);
	lua_pop(L,1);
	lua_check(L,5,indices);
	lua_check(L,6,NF);
	size_t Norbitals = G.size();
	if (Norbitals != F.size() || Norbitals != energies.size() || Norbitals != indToK.size() || Norbitals != indices.size())
	{
		return luaL_error(L,"Error in LuaRelativisticDipoleOperator: Objects of unequal lenght: G : %d\n F : %d\n energies : %d\n indices: %d\n indices[kappas] : %d\n",Norbitals ,F.size(),energies.size(),indices.size(),indToK.size());
	}

	OperatorType *op = (OperatorType *) lua_newuserdata(L,sizeof(OperatorType));
	luaL_setmetatable(L,"Operator_Type");
	CreateRelativisticDipoleOperator(op,G,F,energies,mu,indices,indToK,NF);

	return 1;
}

/*  calculates <psi_i | T_1^+ G_1^+ T_2^+ G_2 T_2 G_1 T_1 | psi_i >
 *  On stack:
 *     H1
 *     H2
 *     T1  - allowed to be a table
 *     T2  - allowed to be a table
 *     psi - allowed to be a table
 *     options
 *        {E1min, E1max, N1, Gamma1}
 *        {E2min, E2max, N2, Gamma2}
 *        NTri1
 *        NTri2
 *        epsilon
 *        restrictions1
 *        restrictions2
 *        E0 E0_1 E0_2
 */
static int LuaCreateResonantSpectraSE(lua_State *L)
{
  struct
  {
    ConvergenceSettingsType Settings;
    double E0, E01, E02;
    uint32_t Npsi;
    uint32_t NT1;
    uint32_t NT2;
    uint32_t NTriSE;
    bool Tensor;
    bool DenseMethod;
  } Options;
  Options.Settings=InitConvergenceSettings();
  Options.E0=DBL_MAX;
  Options.E01=DBL_MAX;
  Options.E02=DBL_MAX;
  Options.Tensor=false;
  Options.DenseMethod=false;
  Options.NTriSE = 5;
  ConvergenceSettingsType SettingsNow;
  SpectraTypeArray *Spec;
  SpectraType *SpecT1psi;
  TriDiagonalMatrixType *SpecTri;
  OperatorType *H1,*H2,*T1, *T2;
  WaveFunctionType *psi, HTpsiDenseBasis,HTHTpsiDenseBasis,Opsi,TpsiSum;
  uint32_t iT,ipsi,len,i;
  LuaOptionsType LuaCreateResonantSpectraOptions[]={
    {"NTri1","uint32_t",&(Options.Settings.NTri1)},
    {"NTri2","uint32_t",&(Options.Settings.NTri2)},
    {"Tensor","bool",&(Options.Tensor)},
    {"restrictions1", "DeterminantRestrictionType", &(Options.Settings.restrictions1)},
    {"restrictions2", "DeterminantRestrictionType", &(Options.Settings.restrictions2)},
    {"Restrictions1", "DeterminantRestrictionType", &(Options.Settings.restrictions1)},
    {"Restrictions2", "DeterminantRestrictionType", &(Options.Settings.restrictions2)},
    {"epsilon","double",&(Options.Settings.epsilon)},
    {"Epsilon","double",&(Options.Settings.epsilon)},
    {"E0","double",&(Options.E0)},
    {"E0_1","double",&(Options.E01)},
    {"E0_2","double",&(Options.E02)},
    {"Emin1","double",&(Options.Settings.Emin1)},
    {"Emax1","double",&(Options.Settings.Emax1)},
    {"NE1","uint32_t",&(Options.Settings.NE1)},
    {"Gamma1","double",&(Options.Settings.Gamma1)},
    {"Emin2","double",&(Options.Settings.Emin2)},
    {"Emax2","double",&(Options.Settings.Emax2)},
    {"NE2","uint32_t",&(Options.Settings.NE2)},
    {"Gamma2","double",&(Options.Settings.Gamma2)},
    {"DenseBorder","uint32_t",&(Options.Settings.DenseBorder)},
    {"ForceDenseMethod","uint32_t",&(Options.DenseMethod)},
	{"NTriSE","uint32_t",&(Options.NTriSE)},
    {NULL, NULL, NULL}};
  /* free garbage full cycle*/
  if(L->l_G->gcrunning) lua_gc(L, LUA_GCCOLLECT, 0);
  /* check number of input parameters */
  len = (uint32_t)lua_gettop(L);
  if(!((len==7)||(len==8)))
  {
    return luaL_error(L, "Error LuaCreateResonantSpectra called with %d arguments,8 or 7 expected (H1, H2, T1, T2^+, psi,U,EnergyContinuum {options})\n",len);
  }
  /* read operator and wavefunctions */
  H1 = (OperatorType *) luaL_checkudata(L, 1, "Operator_Type");
  H2 = (OperatorType *) luaL_checkudata(L, 2, "Operator_Type");
  if(lua_type(L,3)==LUA_TTABLE)
  {
    T1=luaL_checkarray_OperatorType(L, 3, &Options.NT1);
    if(Options.NT1==0) return luaL_error(L, "Errof CreateResonantSpectra called with a table of transition operators (1) of length 0\n");
  }
  else
  {
    Options.NT1 = 1;
    T1 = (OperatorType *) luaL_checkudata(L, 3, "Operator_Type");
  }
  if(lua_type(L,4)==LUA_TTABLE)
  {
    T2=luaL_checkarray_OperatorType(L, 4, &Options.NT2);
    if(Options.NT2==0) return luaL_error(L, "Errof CreateResonantSpectra called with a table of transition operators (2) of length 0\n");
  }
  else
  {
    Options.NT2 = 1;
    T2 = (OperatorType *) luaL_checkudata(L, 4, "Operator_Type");
  }
  if(lua_type(L,5)==LUA_TTABLE)
  {
    psi=luaL_checkarray_WaveFunctionType(L, 5, &Options.Npsi);
    if(Options.Npsi==0) return luaL_error(L, "Errof CreateResonantSpectra called with a table of wave functions of length 0\n");
  }
  else
  {
    Options.Npsi = 1;
    psi = (WaveFunctionType *) luaL_checkudata(L, 5, "Wavefunction_Type");
  }
  std::vector<OperatorType*> U;
  lua_check(L,6,U);
  std::vector<double> EnergyContinuum;
  lua_check(L,7,EnergyContinuum);
  /* read options */
  if (len==8)
  {
    LuaCheckOptions(L, -1, LuaCreateResonantSpectraOptions);
  }
  if(Options.E01==DBL_MAX)Options.E01=Options.E0;
  if(Options.E02==DBL_MAX)Options.E02=Options.E0;

  /* test if operator is Hermitian */
  switch(QOperatorHermitian(*H1, 0.0001*Options.Settings.zero))
  {
    case 3 :
    {
      printf("WARNING: QOperatorHermitian failed in LuaCreateResonantSpectra, this should not happen\n");
      printf("my guess is that you are out of memory and will crash soon\n");
      break;
    }
    case 2 :
    {
      luaL_warning(L, "CreateResonantSpectra called with an operator (H1) that is not fully Hermitian\nStrongly sugested to stop the calculation and repair this\n\n\n\nDo it\n\n\nReally Check this\n\n\nDid you?\n\n\n");
      break;
    }
    case 1 :
    {
      // all fine H is Hermitian upto zero
      break;
    }
    case 0 :
    {
      luaL_error(L, "CreateResonantSpectra called with an operator (H1) that is not Hermitian\n");
      break;
    }
  }
  switch(QOperatorHermitian(*H2, 0.0001*Options.Settings.zero))
  {
    case 3 :
    {
      printf("WARNING: QOperatorHermitian failed in LuaCreateResonantSpectra, this should not happen\n");
      printf("my guess is that you are out of memory and will crash soon\n");
      break;
    }
    case 2 :
    {
      luaL_warning(L, "CreateResonantSpectra called with an operator (H2) that is not fully Hermitian\nStrongly sugested to stop the calculation and repair this\n\n\n\nDo it\n\n\nReally Check this\n\n\nDid you?\n\n\n");
      break;
    }
    case 1 :
    {
      // all fine H is Hermitian upto zero
      break;
    }
    case 0 :
    {
      luaL_error(L, "CreateResonantSpectra called with an operator (H2) that is not Hermitian\n");
      break;
    }
  }
  /* allocate spectra and a pointer to return matrix on lua stack */
  Spec = (SpectraTypeArray *) lua_newuserdata(L, sizeof(SpectraTypeArray));
  luaL_setmetatable(L, "Spectra_Type");
  Spec->NSpectra=0;
  Spec->Spec=NULL;
  /* test if dense methods can be used */
  // at the same time create the wave-functions that define the basis sets of the dense matrices used in the dense methods
  // we need the vectors:
  // HmDenseBasis       : on which the Matrix for the intermediate state can be represented
  // HfDenseBasis       : on which the Matrix fo rthe final states can be represented
		  {
			// psi
			Options.DenseMethod=false;
			goto EndOfUseDenseTest; // TODO Implement DenseMethod
			if (Options.DenseMethod==false) if(psi[0].N>Options.Settings.DenseBorder) goto EndOfUseDenseTest;
			WaveFunctionCopy(psi[0], &HTpsiDenseBasis);
			for(i=1;i<Options.Npsi;i++)
			{
			  WaveFunctionAddUnionBasis(&HTpsiDenseBasis, psi[i]);
			  if (Options.DenseMethod==false) if(HTpsiDenseBasis.N>Options.Settings.DenseBorder)
			  {
				WaveFunctionFree(&HTpsiDenseBasis);
				goto EndOfUseDenseTest;
			  }
			}
			// T psi
			if((Options.Settings.restrictions1.N==0)?(OperatorPsi(T1[0], HTpsiDenseBasis, &TpsiSum)):(OperatorPsiRestricted(T1[0], HTpsiDenseBasis, &TpsiSum, Options.Settings.restrictions1)))
			{
			  WaveFunctionFree(&HTpsiDenseBasis);
			  goto EndOfUseDenseTest;
			}
			for(i=1;i<Options.NT1;i++)
			{
			  if((Options.Settings.restrictions1.N==0)?(OperatorPsi(T1[i], HTpsiDenseBasis, &Opsi)):(OperatorPsiRestricted(T1[i], HTpsiDenseBasis, &Opsi, Options.Settings.restrictions1)))
			  {
				WaveFunctionFree(&TpsiSum);
				WaveFunctionFree(&HTpsiDenseBasis);
				goto EndOfUseDenseTest;
			  }
			  WaveFunctionAddUnionBasis(&TpsiSum, Opsi);
			  WaveFunctionFree(&Opsi);
			  if (Options.DenseMethod==false) if(TpsiSum.N>Options.Settings.DenseBorder)
			  {
				WaveFunctionFree(&TpsiSum);
				WaveFunctionFree(&HTpsiDenseBasis);
				goto EndOfUseDenseTest;
			  }
			}
			WaveFunctionFree(&HTpsiDenseBasis);
			HTpsiDenseBasis=TpsiSum;
			WaveFunctionNormalize(&HTpsiDenseBasis);
			// H T psi
			len=0;
			while(len!=HTpsiDenseBasis.N)
			{
			  len=HTpsiDenseBasis.N;
			  if((Options.Settings.restrictions1.N==0)?(OperatorPsi(*H1, HTpsiDenseBasis, &Opsi)):(OperatorPsiRestricted(*H1, HTpsiDenseBasis, &Opsi, Options.Settings.restrictions1)))
			  {
				WaveFunctionFree(&HTpsiDenseBasis);
				goto EndOfUseDenseTest;
			  }
			  WaveFunctionFree(&HTpsiDenseBasis);
			  if (Options.DenseMethod==false) if(Opsi.N>Options.Settings.DenseBorder)
			  {
				WaveFunctionFree(&Opsi);
				goto EndOfUseDenseTest;
			  }
			  HTpsiDenseBasis=Opsi;
			  WaveFunctionNormalize(&HTpsiDenseBasis);
			}
			// T H T psi
			WaveFunctionCopy(HTpsiDenseBasis, &HTHTpsiDenseBasis);
			if((Options.Settings.restrictions2.N==0)?(OperatorPsi(T2[0], HTHTpsiDenseBasis, &TpsiSum)):(OperatorPsiRestricted(T2[0], HTHTpsiDenseBasis, &TpsiSum, Options.Settings.restrictions2)))
			{
			  WaveFunctionFree(&HTHTpsiDenseBasis);
			  WaveFunctionFree(&HTpsiDenseBasis);
			  goto EndOfUseDenseTest;
			}
			for(i=1;i<Options.NT2;i++)
			{
			  if((Options.Settings.restrictions2.N==0)?(OperatorPsi(T2[i], HTHTpsiDenseBasis, &Opsi)):(OperatorPsiRestricted(T2[i], HTHTpsiDenseBasis, &Opsi, Options.Settings.restrictions2)))
			  {
				WaveFunctionFree(&TpsiSum);
				WaveFunctionFree(&HTHTpsiDenseBasis);
				WaveFunctionFree(&HTpsiDenseBasis);
				goto EndOfUseDenseTest;
			  }
			  WaveFunctionAddUnionBasis(&TpsiSum, Opsi);
			  WaveFunctionFree(&Opsi);
			  if (Options.DenseMethod==false) if(TpsiSum.N>Options.Settings.DenseBorder)
			  {
				WaveFunctionFree(&TpsiSum);
				WaveFunctionFree(&HTHTpsiDenseBasis);
				WaveFunctionFree(&HTpsiDenseBasis);
				goto EndOfUseDenseTest;
			  }
			}
			WaveFunctionFree(&HTHTpsiDenseBasis);
			HTHTpsiDenseBasis=TpsiSum;
			WaveFunctionNormalize(&HTHTpsiDenseBasis);
			// H T H T psi
			len=0;
			while(len!=HTHTpsiDenseBasis.N)
			{
			  len=HTHTpsiDenseBasis.N;
			  if((Options.Settings.restrictions2.N==0)?(OperatorPsi(*H2, HTHTpsiDenseBasis, &Opsi)):(OperatorPsiRestricted(*H2, HTHTpsiDenseBasis, &Opsi, Options.Settings.restrictions2)))
			  {
				WaveFunctionFree(&HTHTpsiDenseBasis);
				WaveFunctionFree(&HTpsiDenseBasis);
				goto EndOfUseDenseTest;
			  }
			  WaveFunctionFree(&HTHTpsiDenseBasis);
			  if (Options.DenseMethod==false) if(Opsi.N>Options.Settings.DenseBorder)
			  {
				WaveFunctionFree(&Opsi);
				WaveFunctionFree(&HTpsiDenseBasis);
				goto EndOfUseDenseTest;
			  }
			  WaveFunctionNormalize(&Opsi);
			  HTHTpsiDenseBasis=Opsi;
			}
			// finsihed all tests and psiDenseBasis still fits
			Options.DenseMethod=true;
		  EndOfUseDenseTest:;
		  }
  if(Options.DenseMethod==true)
  {
    if(Options.Tensor==false)
    {
      ResonantSpectraDense(&HTpsiDenseBasis,&HTHTpsiDenseBasis,*H1, *H2, T1, Options.NT1, T2, Options.NT2, psi, Options.Npsi, &(Spec->Spec), Options.E01, Options.E02, &(Options.Settings));
      Spec->NSpectra=Options.NT1*Options.NT2*Options.Npsi*(Options.Settings.NE1+1);
    }
    else
    {
      return luaL_error(L,"Tensor form for Resonant spectra not yet available");
    }
  }
  else
  {
    // use sparse methods
    if(Options.Tensor==false)
    {
      for(ipsi=0;ipsi<Options.Npsi;ipsi++)
      {
        for(iT=0;iT<Options.NT1;iT++)
        {
          SettingsNow = Options.Settings;
          ResonantSpectraSelfEnergy(*H1, *H2, T1[iT], T2, U, Options.NT2, psi[ipsi], &SpecTri, &SpecT1psi, Options.E01, Options.E02, &(SettingsNow),EnergyContinuum,Options.NTriSE);
          Spec->NSpectra+=(SettingsNow.NE1+1)*Options.NT2;
          Spec->Spec= (SpectraType*)realloc(Spec->Spec, Spec->NSpectra*sizeof(SpectraType));
          for(i=0;i<(SettingsNow.NE1+1)*Options.NT2;i++)
          {
            FreeTriDiagonalMatrix(&(SpecTri[i]));
          }
          free(SpecTri);
          for(i=0;i<(SettingsNow.NE1+1)*Options.NT2;i++)
          {
            sprintf(SpecT1psi[i].Name, "%s T1 %i/%i psi %i/%i",SpecT1psi[i].Name,iT,Options.NT1,ipsi,Options.Npsi);
            Spec->Spec[i+Spec->NSpectra-(SettingsNow.NE1+1)*Options.NT2] = SpecT1psi[i];
          }
          free(SpecT1psi);
        }
      }
    }
    else
    {

      //ResonantSpectraTensor(*H1, *H2, T1, Options.NT1, T2, Options.NT2, psi, Options.Npsi, &SpecTri, &(Spec->Spec[0]), Options.E01, Options.E02, &Options.Settings);
      return luaL_error(L,"Tensor form for Resonant spectra not yet available");
    }
  }
  return 1;
}

static int LuaListU(lua_State *L)
{
	uint32_t Nindices, temp, len, NF;
	uint16_t **indices;
	vector<double> x;
	len = (uint32_t)lua_gettop(L);

	if ( len < 3 ||  len > 5)
	{
		return luaL_error(L,"Error in Martin.ListU: 3 to 5 arguments expected; got %d",len);
	}

	NF = (uint32_t)luaL_checkunsigned(L, 1);
	vector<int> indToK;
	if (len > 3)
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

		if (len == 4)
		{
			CreateListOperatorAtomicURelativistic(NF,wavefuns,x,indToK,indices,true,L,false);
		}
		else
		{
			CreateListOperatorAtomicURelativistic(NF,wavefuns,x,indToK,indices,true,L,true);
			for(uint32_t i = 0; i < Nindices; i++ )
			{
				free(indices[i]);
			}
			free(indices);
			return 1;
		}
	}
	else
	{
		printf("error: nonrelativistic version of ListU not yet implmented\n");
		return 0;
	}


	for(uint32_t i = 0; i < Nindices; i++ )
	{
		free(indices[i]);
	}
	free(indices);

	return 2;
}

static int LuaCDMPlot(lua_State *L)
{
	std::vector<double> poles, SlaterInts;
	std::vector<mypoint> spectrum;
	double wMin, wMax, dw, eps, gamma, shift;

	CompactMatrixType *M;

	uint32_t len = lua_gettop(L);
	lua_check(L,1,wMin);
	lua_check(L,2,wMax);
	lua_check(L,3,dw);
	lua_check(L,4,eps);
	lua_check(L,5,gamma);
	lua_check(L,6,poles);
	lua_check(L,7,SlaterInts);
	lua_check(L,8,M);
	lua_check(L,9,shift);

	uint32_t Npoints = floor((wMax-wMin)/dw) + 1;
	spectrum = std::vector<mypoint>(Npoints);
	for (uint32_t i = 0; i < Npoints; i++)
	{
		spectrum[i].w = wMin + i * dw;
		spectrum[i].y.Re = 0.;
		spectrum[i].y.Im = 0.;
	}
	if (len == 10)
	{
		vector<vector<Complex16> > selfenergy;
		lua_check(L,10,selfenergy);
		PlotSpectrumFromPoles(eps,gamma,shift,poles,SlaterInts,M,spectrum,selfenergy);
	}
	else
	{
		PlotSpectrumFromPoles(eps,gamma,shift,poles,SlaterInts,M,spectrum);
	}

	lua_push(L,spectrum);
	return 1;

}

static int LuaAugerMeitnerSelfEnergy(lua_State *L)
{
	WaveFunctionType * psi;
	vector<OperatorType *> UList;
	OperatorType *H;
	vector<double> poles, SlaterInts, EnergiesUnbound;
	vector<mypoint> spectrum;
	double eps, gamma;
	vector<Interpolation::InterpolatingFunction *> WavefunsBoundG,WavefunsBoundF,WavefunsUnboundG,WavefunsUnboundF;
	uint32_t NF;
	vector<int> indToK;
	uint16_t **indices;
	int arg = 1;

	lua_check(L,arg++,WavefunsBoundG);
	lua_check(L,arg++,WavefunsBoundF);
	lua_check(L,arg++,WavefunsUnboundG);
	lua_check(L,arg++,WavefunsUnboundF);
	 NF = (uint32_t)luaL_checkunsigned(L, arg++);
	 uint32_t Nindices = GetIndices(L,arg++,indToK,indices);
	 lua_check(L,arg++,H);
	 lua_check(L,arg++,psi);
	 lua_check(L,arg++,eps);
	 lua_check(L,arg++,gamma);
	 lua_check(L,arg++,EnergiesUnbound);
	 InitializeSpectrum(L,arg,spectrum);

	 uint32_t NPointsToEvaluateInterpolation = 0;
	 if (lua_gettop(L) == arg)
		 lua_check(L,arg,NPointsToEvaluateInterpolation);

	 AugerMeitnerSelfEnergy(WavefunsBoundG,WavefunsBoundF,WavefunsUnboundG,WavefunsUnboundF,indToK,indices,NF,H,psi,eps,gamma,EnergiesUnbound,spectrum,NPointsToEvaluateInterpolation);
	 lua_push(L,spectrum);

	 for (uint32_t i = 0; i < Nindices; i++)
	 {
		 free(indices[i]);
	 }
	 free(indices);
	 return 1;
}

static int LuaBraketList(lua_State *L)
{
	WaveFunctionType * psi;
	vector<OperatorType *> UList;
	OperatorType *H;
	CompactMatrixType *S, *M;

	lua_check(L,1,psi);
	lua_check(L,2,UList);
	lua_check(L,3,H);

	M = (CompactMatrixType*) lua_newuserdata(L,sizeof(CompactMatrixType));
	luaL_setmetatable(L,"Matrix_Type");
	M->QComplex = psi->QComplex || H->QComplex;
	uint32_t i = 0;
	while (!M->QComplex && i < (uint32_t) UList.size())
	{
		M->QComplex = UList[i]->QComplex;
		i++;
	}
	M->Ni = (uint32_t) UList.size();
	M->Nj = M->Ni;
	InitCompactMatrix(M);
	S = (CompactMatrixType*) lua_newuserdata(L,sizeof(CompactMatrixType));
	luaL_setmetatable(L,"Matrix_Type");
	S->QComplex = 0;
	S->Ni = M->Ni;
	S->Nj = M->Ni;
	InitCompactMatrix(S);
	PsiOppListOperatorOppListPsi(psi,UList,H,M);
	printf("Done\n");
	fflush(stdout);
	//PsiOppListOppListPsi(psi,UList,S);

	return 2;
}

static int LuaReduceOperator(lua_State *L)
{
	vector<uint16_t> ind;
	OperatorType * Opp, *ReducedOpp;

	lua_check(L,1,Opp);
	lua_check(L,2,ind);

	ReducedOpp = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
	luaL_setmetatable(L, "Operator_Type");
	ReduceOperator(*Opp,ReducedOpp,ind);
	return 1;
}

static int LuaMartinTest(lua_State *L)
{

	int    info = 0;
	const char *jobz = "V";
	const char *uplo = "U";
	uint32_t N = 10;
	uint32_t    lwork = 1 + 6 * N + 2 * Sqr(N);
	uint32_t    liwork = 3 + 5 * N;

	double *work = (double*)malloc(lwork*sizeof(double));
	double *m = (double*) malloc(N*N*sizeof(double));
	double *val = (double*) malloc(N*N*sizeof(double));

	for (uint32_t i = 0; i < N*N; i++)
	{
		m[i] = 1;
		val[i] = 0;
	}
	for (uint32_t i = 0; i < lwork; i++)
	{
		work[i] = 0;
	}
	dsyev(jobz, uplo, &N, m, &N, val, work, &lwork, &info);
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
  {"OperatorToCSRMatrix",LuaOperatorToCSRMatrix},
  {"CSRTridiagonalize",LuaCSRTridiagonalize},
  {"Export",LuaExport},
  {"ReadFPLOEnergies",LuaReadFPLOEnergies},
  {"RelativisticDipoleOperator",LuaRelativisticDipoleOperator},
  //{"ResonantSpectrum",LuaCreateResonantSpectraSE},
  {"ResonantSpectrum",LuaCreateResonantSpectrum},
  {"ListU",LuaListU},
  {"PlotSpectrum",LuaCDMPlot},
  {"BraketList",LuaBraketList},
  {"AMSel",LuaAugerMeitnerSelfEnergy},
  {"ReduceOperator",LuaReduceOperator},
  {"CreateContinuumSlaterIntegrals",LuaCreateContinuumSlaterIntegrals},
  {"CreateAugerMeitnerU",LuaCreateAugerMeitnerU},
  {"CNSE",LuaCoolNewSelfEnergy},
  {"CNSEDiag",LuaCoolNewSelfEnergyDiagonal},
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

/* DEBUG QUAD
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
 */
