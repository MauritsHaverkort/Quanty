 
TightBindingDef reduceHilbertSpaceForTightBinding(const TightBindingDef &tb, const vector<AtomType> &consideredAtoms, const vector<AtomType> &unchangedAtoms, const uint filledLigands, const uint emptyLigands)
{
//TODO
	cout << "reduceHilbertSpace for TightBinding not yet fully implemented." << endl;
}

OperatorType reduceHilbertSpaceForOperator(OperatorType &hamil, const uint unchangedStates, const uint filledLigands, const uint emptyLigands)
{
//	TODO
	cout << "reduceHilbertSpace for Operator not yet fully implemented." << endl;
}

CompactMatrixType reducedHilbertSpaceMatrix(const CompactMatrixType &H, const CompactMatrixType &ligandRotMat, const uint unchangedStates, const uint ligandNum)
{
	vector < vector<CompactMatrixType> > fullRotMatMat = CompactnullMatrixMatrix(2,2);
	fullRotMatMat[0][0] = CompactidentityMatrix(unchangedStates);
	fullRotMatMat[1][1] = ligandRotMat;
	CompactMatrixType fullRotMat = MatrixFlatten(fullRotMatMat);
	CompactMatrixType HReduced;
	TransformCompactMatrixcT(&H,&fullRotMat, &HReduced);
	CompactMatrixChop(&HReduced, 100*DBL_EPSILON);
	CompactMatrixType startBlock;
	startBlock.Ni=unchangedStates;
	startBlock.Nj=HReduced.Nj;
	InitCompactIdentityMatrix(&startBlock);
	BlockTriDiagonalMatrixType HTriReduced;
	CompactMatrixBlockBandDiagonalize(&HReduced, &startBlock, &HTriReduced, 1E-6, ligandNum + 2, UINT32_MAX, true, true);
	return BlockTriDiagonalMatrixToCompactMatrix(HTriReduced, 100*DBL_EPSILON);
}

CompactMatrixType reduceHilbertSpaceForMatrix(const CompactMatrixType &H, const uint unchangedStates, const uint filledLigands, const uint emptyLigands)
{
	cout << "reduceHilbertSpace for Matrix not yet fully implemented." << endl;

	uint subMatSize = H.Ni - unchangedStates;
	CompactMatrixType HLigand = Submatrix(H, unchangedStates, subMatSize, unchangedStates, subMatSize);
	CompactMatrixType eigenfuns;
	double eigenvals[HLigand.Ni];
	CompactMatrixEigensystem(&HLigand, &eigenfuns, eigenvals);
	uint numNegVals = HLigand.Ni;
	for (uint i = 0; i < (uint) HLigand.Ni; ++i)
	{
		if (eigenvals[i] > 0)
		{
			numNegVals = i;
			break;
		}
	}
	CompactMatrixType ligandRotMat = Submatrix(eigenfuns, 0, numNegVals, 0, eigenfuns.Ni);
	CompactMatrixType HReducedFilledShells = reducedHilbertSpaceMatrix(H, ligandRotMat, unchangedStates, filledLigands);
	CompactMatrixChop(&HReducedFilledShells, 100 * DBL_EPSILON);
	if (emptyLigands > 0 && numNegVals < HLigand.Ni)
	{
		ligandRotMat = Submatrix(eigenfuns, numNegVals, eigenfuns.Ni - numNegVals, 0, eigenfuns.Ni);
		CompactMatrixType HReducedEmptyShells = reducedHilbertSpaceMatrix(H, ligandRotMat, unchangedStates, emptyLigands);
		if (HReducedEmptyShells.Ni > unchangedStates)
		{
			vector < vector<CompactMatrixType> > HReducedMatMat = CompactnullMatrixMatrix(3, 3);
			HReducedMatMat[0][0] = Submatrix(HReducedFilledShells, 0, unchangedStates, 0, unchangedStates);
			HReducedMatMat[0][1] = Submatrix(HReducedFilledShells, 0, unchangedStates, unchangedStates, HReducedFilledShells.Nj - unchangedStates);
			HReducedMatMat[1][0] = Submatrix(HReducedFilledShells, unchangedStates, HReducedFilledShells.Ni - unchangedStates, 0, unchangedStates);
			HReducedMatMat[1][1] = Submatrix(HReducedFilledShells, unchangedStates, HReducedFilledShells.Ni - unchangedStates, unchangedStates, HReducedFilledShells.Nj - unchangedStates);
			HReducedMatMat[0][2] = Submatrix(HReducedEmptyShells, 0, unchangedStates, unchangedStates, HReducedEmptyShells.Nj - unchangedStates);
			HReducedMatMat[2][0] = Submatrix(HReducedEmptyShells, unchangedStates, HReducedEmptyShells.Ni - unchangedStates, 0, unchangedStates);
			HReducedMatMat[2][2] = Submatrix(HReducedEmptyShells, unchangedStates, HReducedEmptyShells.Ni - unchangedStates, unchangedStates, HReducedEmptyShells.Nj - unchangedStates);
			return MatrixFlatten(HReducedMatMat);
		} else
		{
			return HReducedFilledShells;
		}
	} else
	{
		return HReducedFilledShells;
	}
}





//Input variant 1: TB object, Atoms to be considered, Atoms to be left unchanged, options
//Input variant 2: Hamiltonian, number of states to be left unchagned, options
//Input variant 3: Matrix Representation of Hamiltonian, number of states to be left unchanged, options
//Options: HoppingOrder (int declaring to what order the Ligand orbitals are included, default 200), includeEmpty (bool declaring if I include states with positive energy, default false)
static int LuaReduceHilbertSpace(lua_State * L)
{
	TightBindingDef tb;
	OperatorType op;
	CompactMatrixType mat;
  struct
  {
		uint filledLigands = 200;
		uint emptyLigands = 0;
  } Options;
  LuaOptionsType LuaReduceHilbertSpaceOptions[]={
  	{"filledLigands",  		 "int",  &(Options.filledLigands)},
		{"emptyLigands",   "int", &(Options.emptyLigands)},
    {NULL, NULL, NULL}};

  uint32_t len = (uint32_t) lua_gettop(L);
  if (len < 2 || len > 4)
  	return luaL_error(L,"Error: ReduceHilbertSpace called with %d elements, 2, 3 or 4 expected. \nVariant 1: Tight-Binding Object, list of Atoms to be considered, Atoms to be left unchanged, and optional Options. \nVariant 2: Operator, number of states to be left unchanged, and optional Options. \nVariant 3: Matrix Representation of Operator, number of states to be left unchanged, and optional Options.",len);

  if (LuaTestOptionsList(L, -1))
  			LuaCheckOptions(L, -1, LuaReduceHilbertSpaceOptions);

	if (luaL_isType(L, 1, tb))
	{
		lua_check(L, 1, tb);
//		TODO: QTightBindingHermitian
		vector<AtomType> consideredAtoms;
		lua_check(L, 2, consideredAtoms);
		vector<AtomType> unchangedAtoms;
		lua_check(L, 3, unchangedAtoms);
		lua_push(L, reduceHilbertSpaceForTightBinding(tb, consideredAtoms, unchangedAtoms, Options.filledLigands, Options.emptyLigands));
	} else if (luaL_isType(L, 1, op))
	{
		lua_check(L, 1, op);
		if (!(bool)QOperatorHermitian(op, 10*DBL_EPSILON))
			return luaL_error(L,"Error: ReduceHilbertSpace expects a Hermitian operator as input.");
		uint unchangedStates;
		lua_check(L, 2, unchangedStates);
		lua_push(L, reduceHilbertSpaceForOperator(op, unchangedStates, Options.filledLigands, Options.emptyLigands));
	} else if (luaL_isType(L, 1, mat))
	{
		lua_check(L, 1, mat);
		if (mat.Ni != mat.Nj)
			return luaL_error(L,"Error: ReduceHilbertSpace expects a quadratic matrix as input.");
		if (!(bool)QCompactMatrixHermitian(mat, 10*DBL_EPSILON))
			return luaL_error(L,"Error: ReduceHilbertSpace expects a Hermitian matrix as input.");
		uint unchangedStates;
		lua_check(L, 2, unchangedStates);
		if (unchangedStates >= mat.Ni)
			return luaL_error(L,"Error: Number of unchanged states (%d) must be smaller than dimension of matrix (%d).",unchangedStates, mat.Ni);
		lua_push(L, reduceHilbertSpaceForMatrix(mat, unchangedStates, Options.filledLigands, Options.emptyLigands));
	} else
	{
		return luaL_error(L, "Error: ReduceHilbertSpace only implemented for Tight-Binding Objects, Operators and Matrices.");
	}
	return 1;
}

