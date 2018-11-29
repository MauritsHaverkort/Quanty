 
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

CompactMatrixType reducedHilbertSpaceMatrix(const CompactMatrixType &H, const CompactMatrixType &ligandRotMat, CompactMatrixType &T, const uint unchangedStates, const uint ligandNum)
{
	if (isNullMatrix(ligandRotMat))
	{
		T = CompactnullMatrix();
		return CompactnullMatrix();
	}

	CompactMatrixType out, HReduced, BlockBandTrafo, BlockBandTrafoToDiag;
	CompactMatrixType fullRotMat = SubspaceRotationMatrix(unchangedStates, ligandRotMat, 0);
	TransformCompactMatrixcT(&H, &fullRotMat, &HReduced);
	CompactMatrixChop(&HReduced, 100 * DBL_EPSILON);
	BlockBandTrafo = CompactMatrixType();
	BlockBandTrafo.Ni = unchangedStates;
	BlockBandTrafo.Nj = HReduced.Nj;
	InitCompactIdentityMatrix(&BlockBandTrafo);
	BlockTriDiagonalMatrixType HTriReduced;
	CompactMatrixBlockBandDiagonalize(&HReduced, &BlockBandTrafo, &HTriReduced, 1E-6, ligandNum + 2, UINT32_MAX, true, true);
	CompactMatrixMultiply('T', 'T', &fullRotMat, &BlockBandTrafo, &BlockBandTrafoToDiag);
	CompactMatrixTranspose(&BlockBandTrafoToDiag);
	out = BlockTriDiagonalMatrixToCompactMatrix(HTriReduced, 100 * DBL_EPSILON);
	CompactMatrixCutDimensions(out);
	CompactMatrixCopy(BlockBandTrafoToDiag, &T);
	T = Submatrix(BlockBandTrafoToDiag, out.Ni, BlockBandTrafoToDiag.Nj);
	FreeCompactMatrix(&HReduced);
	FreeCompactMatrix(&BlockBandTrafo);
	FreeCompactMatrix(&BlockBandTrafoToDiag);
	FreeCompactMatrix(&fullRotMat);
	FreeBlockTridiagonalMatrix(&HTriReduced);
	return out;
}

CompactMatrixType reducedHilbertSpaceMatrixFilledAndEmpty(const CompactMatrixType &HRedFilled, const CompactMatrixType &HRedEmpty, const uint unchangedStates)
{
	CompactMatrixType out;
	if (HRedFilled.Ni <= unchangedStates)
	{
		CompactMatrixCopy(HRedEmpty, &out);
	} else if (HRedEmpty.Ni <= unchangedStates)
	{
		CompactMatrixCopy(HRedFilled, &out);
	} else
	{
		vector < vector<CompactMatrixType> > HReducedMatMat = CompactnullMatrixMatrix(3, 3);
		HReducedMatMat[0][0] = Submatrix(HRedFilled, 0, unchangedStates, 0, unchangedStates);
		HReducedMatMat[0][1] = Submatrix(HRedFilled, 0, unchangedStates, unchangedStates, HRedFilled.Nj - unchangedStates);
		HReducedMatMat[1][0] = Submatrix(HRedFilled, unchangedStates, HRedFilled.Ni - unchangedStates, 0, unchangedStates);
		HReducedMatMat[1][1] = Submatrix(HRedFilled, unchangedStates, HRedFilled.Ni - unchangedStates, unchangedStates, HRedFilled.Nj - unchangedStates);
		HReducedMatMat[0][2] = Submatrix(HRedEmpty, 0, unchangedStates, unchangedStates, HRedEmpty.Nj - unchangedStates);
		HReducedMatMat[2][0] = Submatrix(HRedEmpty, unchangedStates, HRedEmpty.Ni - unchangedStates, 0, unchangedStates);
		HReducedMatMat[2][2] = Submatrix(HRedEmpty, unchangedStates, HRedEmpty.Ni - unchangedStates, unchangedStates, HRedEmpty.Nj - unchangedStates);
		out = MatrixFlatten(HReducedMatMat);
		FreeCompactMatrixMatrix (HReducedMatMat);
	}
	return out;
}

CompactMatrixType reducedHilbertSpaceRotMatFilledAndEmpty(const CompactMatrixType &Tfilled, const CompactMatrixType &Tempty, const uint unchangedStates)
{
	CompactMatrixType out;
	if (Tfilled.Ni <= unchangedStates)
	{
		CompactMatrixCopy(Tempty, &out);
	} else if (Tempty.Ni <= unchangedStates)
	{
		CompactMatrixCopy(Tfilled, &out);
	} else
	{
		vector < vector<CompactMatrixType> > TMatMat = CompactnullMatrixMatrix(2, 1);
		TMatMat[0][0] = Submatrix(Tfilled, unchangedStates, Tfilled.Ni - unchangedStates, unchangedStates, Tfilled.Nj - unchangedStates);
		TMatMat[1][0] = Submatrix(Tempty, unchangedStates, Tempty.Ni - unchangedStates, unchangedStates, Tempty.Nj - unchangedStates);
		out = SubspaceRotationMatrix(unchangedStates, MatrixFlatten(TMatMat), 0);
		FreeCompactMatrixMatrix (TMatMat);
	}
	return out;
}

CompactMatrixType reduceHilbertSpaceForMatrix(const CompactMatrixType &H, CompactMatrixType &T, const uint unchangedStates, const uint filledLigands, const uint emptyLigands)
{
	CompactMatrixType out, HLigand, eigenfunsLigand, HRedFilled, Tfilled, HRedEmpty, Tempty;
	uint LigandMatSize = H.Ni - unchangedStates;
	HLigand = Submatrix(H, unchangedStates, LigandMatSize, unchangedStates, LigandMatSize);
	eigenfunsLigand = CompactMatrixType();
	vector<double> eigenvals = vector<double>(LigandMatSize, 0);
	CompactMatrixEigensystem(&HLigand, &eigenfunsLigand, &eigenvals[0]);
	uint numNegVals = LigandMatSize;
	for (uint i = 0; i < (uint) eigenvals.size(); ++i)
	{
		if (eigenvals[i] > 0)
		{
			numNegVals = i;
			break;
		}
	}
	if (numNegVals > 0 && filledLigands > 0)
	{
		CompactMatrixType filledLigRotMatToDiag = Submatrix(eigenfunsLigand, 0, numNegVals, 0, eigenfunsLigand.Ni);
		HRedFilled = reducedHilbertSpaceMatrix(H, filledLigRotMatToDiag, Tfilled, unchangedStates, filledLigands);
		CompactMatrixChop(&HRedFilled, 100 * DBL_EPSILON);
		FreeCompactMatrix(&filledLigRotMatToDiag);
	} else
	{
		HRedFilled = CompactnullMatrix();
		Tfilled = CompactnullMatrix();
	}
	if (numNegVals < eigenfunsLigand.Ni && emptyLigands > 0)
	{
		CompactMatrixType emptyLigRotMatToDiag = Submatrix(eigenfunsLigand, numNegVals, eigenfunsLigand.Ni - numNegVals, 0, eigenfunsLigand.Ni);
		HRedEmpty = reducedHilbertSpaceMatrix(H, emptyLigRotMatToDiag, Tempty, unchangedStates, emptyLigands);
		CompactMatrixChop(&HRedEmpty, 100 * DBL_EPSILON);
		FreeCompactMatrix(&emptyLigRotMatToDiag);
	} else
	{
		HRedEmpty = CompactnullMatrix();
		Tempty = CompactnullMatrix();
	}
	out = reducedHilbertSpaceMatrixFilledAndEmpty(HRedFilled, HRedEmpty, unchangedStates);
	T = reducedHilbertSpaceRotMatFilledAndEmpty(Tfilled, Tempty, unchangedStates);
	FreeCompactMatrix(&HLigand);
	FreeCompactMatrix(&eigenfunsLigand);
	FreeCompactMatrix(&HRedFilled);
	FreeCompactMatrix(&Tfilled);
	FreeCompactMatrix(&HRedEmpty);
	FreeCompactMatrix(&Tempty);
	return out;
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
		CompactMatrixType T;
		lua_push(L, reduceHilbertSpaceForMatrix(mat, T, unchangedStates, Options.filledLigands, Options.emptyLigands));
		lua_push(L, T);
	} else
	{
		return luaL_error(L, "Error: ReduceHilbertSpace only implemented for Tight-Binding Objects, Operators and Matrices.");
	}
	return 2;
}

