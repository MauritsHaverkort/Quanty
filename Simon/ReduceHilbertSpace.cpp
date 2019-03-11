CompactMatrixType reducedHilbertSpaceMatrix(const CompactMatrixType &H, const CompactMatrixType &ligRotMat, CompactMatrixType &T, const uint constOrbs, const uint ligNum)
{
	if (isNullMatrix(ligRotMat))
	{
		T = CompactnullMatrix();
		return CompactnullMatrix();
	}

	CompactMatrixType out, HReduced, BlockBandTrafo, BlockBandTrafoToDiag;
	CompactMatrixType fullRotMat = SubspaceRotationMatrix(constOrbs, ligRotMat, 0);
	TransformCompactMatrixcT(&H, &fullRotMat, &HReduced);
	CompactMatrixChop(&HReduced, 100 * DBL_EPSILON);
	BlockBandTrafo = CompactMatrixType();
	BlockBandTrafo.Ni = constOrbs;
	BlockBandTrafo.Nj = HReduced.Nj;
	InitCompactIdentityMatrix(&BlockBandTrafo);
	BlockTriDiagonalMatrixType HTriReduced;
	CompactMatrixBlockBandDiagonalize(&HReduced, &BlockBandTrafo, &HTriReduced, 1E-6, ligNum + 2, UINT32_MAX, true, true);
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

CompactMatrixType reducedHilbertSpaceMatrixFilledAndEmpty(const CompactMatrixType &HRedFilled, const CompactMatrixType &HRedEmpty, const uint constOrbs)
{
	CompactMatrixType out;
	if (HRedFilled.Ni <= constOrbs)
	{
		CompactMatrixCopy(HRedEmpty, &out);
	} else if (HRedEmpty.Ni <= constOrbs)
	{
		CompactMatrixCopy(HRedFilled, &out);
	} else
	{
		vector < vector<CompactMatrixType> > HReducedMatMat = CompactnullMatrixMatrix(3, 3);
		HReducedMatMat[0][0] = Submatrix(HRedFilled, 0, constOrbs, 0, constOrbs);
		HReducedMatMat[0][1] = Submatrix(HRedFilled, 0, constOrbs, constOrbs, HRedFilled.Nj - constOrbs);
		HReducedMatMat[1][0] = Submatrix(HRedFilled, constOrbs, HRedFilled.Ni - constOrbs, 0, constOrbs);
		HReducedMatMat[1][1] = Submatrix(HRedFilled, constOrbs, HRedFilled.Ni - constOrbs, constOrbs, HRedFilled.Nj - constOrbs);
		HReducedMatMat[0][2] = Submatrix(HRedEmpty, 0, constOrbs, constOrbs, HRedEmpty.Nj - constOrbs);
		HReducedMatMat[2][0] = Submatrix(HRedEmpty, constOrbs, HRedEmpty.Ni - constOrbs, 0, constOrbs);
		HReducedMatMat[2][2] = Submatrix(HRedEmpty, constOrbs, HRedEmpty.Ni - constOrbs, constOrbs, HRedEmpty.Nj - constOrbs);
		out = MatrixFlatten(HReducedMatMat);
		FreeCompactMatrixMatrix (HReducedMatMat);
	}
	return out;
}

CompactMatrixType reducedHilbertSpaceRotMatFilledAndEmpty(const CompactMatrixType &Tfilled, const CompactMatrixType &Tempty, const uint constOrbs, const uint Nj)
{
	CompactMatrixType out;
	if (Tfilled.Ni <= constOrbs && Tempty.Ni <= constOrbs)
	{
		out.Ni = constOrbs;
		out.Nj = Nj;
		out.QComplex = 0;
		InitCompactIdentityMatrix(&out);
	} else if (Tfilled.Ni <= constOrbs)
	{
		CompactMatrixCopy(Tempty, &out);
	} else if (Tempty.Ni <= constOrbs)
	{
		CompactMatrixCopy(Tfilled, &out);
	} else
	{
		vector < vector<CompactMatrixType> > TMatMat = CompactnullMatrixMatrix(2, 1);
		TMatMat[0][0] = Submatrix(Tfilled, constOrbs, Tfilled.Ni - constOrbs, constOrbs, Tfilled.Nj - constOrbs);
		TMatMat[1][0] = Submatrix(Tempty, constOrbs, Tempty.Ni - constOrbs, constOrbs, Tempty.Nj - constOrbs);
		out = SubspaceRotationMatrix(constOrbs, MatrixFlatten(TMatMat), 0);
		FreeCompactMatrixMatrix (TMatMat);
	}
	return out;
}

void addLigandNames(const CompactMatrixType &M, const string &nameBase, vector<string> &names, const uint constOrbs, const bool hasSpin)
{
	uint sameSymLigands = (uint) (M.Ni / constOrbs) - 1;
	uint i;
	for (i = 0; i < sameSymLigands; ++i)
		names.push_back(nameBase + to_string(i + 1) + "_" + multiplicityToNonRelOrb(constOrbs, hasSpin));
	uint restLigandMult = M.Ni % constOrbs;
	if (restLigandMult > 0)
		names.push_back(nameBase + to_string(i + 1) + "_" + multiplicityToNonRelOrb(restLigandMult, hasSpin));
}

CompactMatrixType reduceHilbertSpaceForMatrix(const CompactMatrixType &H, CompactMatrixType &T, const uint constOrbs, const uint filledLig, const uint emptyLig, vector<string> &orbNames,
		const bool hasSpin)
{
	CompactMatrixType out, HLigand, eigenfunsLigand, HRedFilled, Tfilled, HRedEmpty, Tempty;
	uint LigandMatSize = H.Ni - constOrbs;
	orbNames = vector<string>();
	orbNames.push_back("unchanged_" + multiplicityToNonRelOrb(constOrbs, hasSpin));
	HLigand = Submatrix(H, constOrbs, LigandMatSize, constOrbs, LigandMatSize);
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
	if (numNegVals > 0 && filledLig > 0)
	{
		CompactMatrixType filledLigRotMatToDiag = Submatrix(eigenfunsLigand, 0, numNegVals, 0, eigenfunsLigand.Ni);
		HRedFilled = reducedHilbertSpaceMatrix(H, filledLigRotMatToDiag, Tfilled, constOrbs, filledLig);
		CompactMatrixChop(&HRedFilled, 100 * DBL_EPSILON);
		addLigandNames(HRedFilled, "Lval", orbNames, constOrbs, hasSpin);
		FreeCompactMatrix(&filledLigRotMatToDiag);
	} else
	{
		HRedFilled = CompactnullMatrix();
		Tfilled = CompactnullMatrix();
	}
	if (numNegVals < eigenfunsLigand.Ni && emptyLig > 0)
	{
		CompactMatrixType emptyLigRotMatToDiag = Submatrix(eigenfunsLigand, numNegVals, eigenfunsLigand.Ni - numNegVals, 0, eigenfunsLigand.Ni);
		HRedEmpty = reducedHilbertSpaceMatrix(H, emptyLigRotMatToDiag, Tempty, constOrbs, emptyLig);
		CompactMatrixChop(&HRedEmpty, 100 * DBL_EPSILON);
		addLigandNames(HRedEmpty, "Lcon", orbNames, constOrbs, hasSpin);
		FreeCompactMatrix(&emptyLigRotMatToDiag);
	} else
	{
		HRedEmpty = CompactnullMatrix();
		Tempty = CompactnullMatrix();
	}
	out = reducedHilbertSpaceMatrixFilledAndEmpty(HRedFilled, HRedEmpty, constOrbs);
	T = reducedHilbertSpaceRotMatFilledAndEmpty(Tfilled, Tempty, constOrbs, H.Nj);
	FreeCompactMatrix(&HLigand);
	FreeCompactMatrix(&eigenfunsLigand);
	FreeCompactMatrix(&HRedFilled);
	FreeCompactMatrix(&Tfilled);
	FreeCompactMatrix(&HRedEmpty);
	FreeCompactMatrix(&Tempty);
	return out;
}

CompactMatrixType reduceHilbertSpaceForMatrix(const CompactMatrixType &H, CompactMatrixType &T, const uint constOrbs, const uint filledLig, const uint emptyLig)
{
	vector < string > vec = vector<string>();
	return reduceHilbertSpaceForMatrix(H, T, constOrbs, filledLig, emptyLig, vec, true);
}

CompactMatrixType reduceHilbertSpaceForMatrix(const CompactMatrixType &H, const uint constOrbs, const uint filledLig, const uint emptyLig)
{
	CompactMatrixType T;
	vector < string > vec = vector<string>();
	CompactMatrixType out = reduceHilbertSpaceForMatrix(H, T, constOrbs, filledLig, emptyLig, vec, true);
	FreeCompactMatrix(&T);
	return out;
}

OperatorType reduceHilbertSpaceForOperator(const OperatorType &hamil, CompactMatrixType &T, const uint constOrbs, const uint filledLig, const uint emptyLig, vector<string> &orbNames,
		const bool hasSpin)
{
	CompactMatrixType Rho, HMat, HMatRed;
	OperatorType out;
	Rho = CompactMatrixType();
	Rho.Ni = hamil.NFermion;
	Rho.Nj = hamil.NFermion;
	Rho.QComplex = 0;
	InitCompactMatrix(&Rho);
	HMat = CompactMatrixType();
	MeanFieldOperatorOneElectronMatrix(&hamil, &HMat, Rho);
	HMatRed = reduceHilbertSpaceForMatrix(HMat, T, constOrbs, filledLig, emptyLig, orbNames, hasSpin);
	out = CompactMatrixToOperator(HMatRed);
	FreeCompactMatrix(&Rho);
	FreeCompactMatrix(&HMat);
	FreeCompactMatrix(&HMatRed);
	return out;
}

OperatorType reduceHilbertSpaceForOperator(const OperatorType &H, CompactMatrixType &T, const uint constOrbs, const uint filledLig, const uint emptyLig)
{
	vector < string > vec = vector<string>();
	return reduceHilbertSpaceForOperator(H, T, constOrbs, filledLig, emptyLig, vec, true);
}

OperatorType reduceHilbertSpaceForOperator(const OperatorType &H, const uint constOrbs, const uint filledLig, const uint emptyLig)
{
	CompactMatrixType T;
	vector < string > vec = vector<string>();
	OperatorType out = reduceHilbertSpaceForOperator(H, T, constOrbs, filledLig, emptyLig, vec, true);
	FreeCompactMatrix(&T);
	return out;
}

TightBindingDef reduceHilbertSpaceForTightBinding(TightBindingDef &tb, CompactMatrixType &T, const AtomType &constAtom, const uint filledLig, const uint emptyLig, const bool hasSpin)
{
	TightBindingDef tbSorted, out = TightBindingDef();
	MatrixType hMat;
	CompactMatrixType hMatCompact, permMat, rotMat, hMatRed;
	if (!QAtomsEqual(constAtom, tb.Cell.Atom[0]))
	{
		tbSorted = sortTbAtomToFirstPosition(tb, constAtom, permMat);
	} else
	{
		TightBindingCopy(tb, tbSorted);
		permMat = CompactidentityMatrix(tb.NFermion);
	}
	double k[3] = { 0., 0., 0. };
	CreateHamiltonianFromTightBindingDefConstantPhaseInUnitCell(tbSorted, k, &hMat);
	hMatCompact = MatrixTypeToCompactMatrixType(hMat);
	uint constOrbs = constAtom.shelldegeneracy[0];
	vector < string > newNames = vector<string>();
	hMatRed = reduceHilbertSpaceForMatrix(hMatCompact, rotMat, constOrbs, filledLig, emptyLig, newNames, hasSpin);
	CompactMatrixMultiply('T', 'T', &permMat, &rotMat, &T);
	CompactMatrixTranspose(&T);
	TightBindingReadyForRotate(tbSorted, out, OrbsToAtoms(newNames));
	RotateTightBinding(&tbSorted, &rotMat, &out);
	FreeTightBindingDef(&tbSorted);
	FreeMatrix(&hMat);
	FreeCompactMatrix(&hMatCompact);
	FreeCompactMatrix(&permMat);
	FreeCompactMatrix(&rotMat);
	FreeCompactMatrix(&hMatRed);
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
	CompactMatrixType mat, T;
  struct
  {
		uint filledLigands = 1;
		uint emptyLigands = 0;
		bool hasSpin = true;
  } Options;
  LuaOptionsType LuaReduceHilbertSpaceOptions[]={
  	{"filledLigands",  		 "int",  &(Options.filledLigands)},
		{"emptyLigands",   "int", &(Options.emptyLigands)},
		{"hasSpin",   "bool", &(Options.hasSpin)},
    {NULL, NULL, NULL}};

  uint32_t len = (uint32_t) lua_gettop(L);
  if (len < 2 || len > 4)
  	return luaL_error(L,"Error: ReduceHilbertSpace called with %d elements, 2 or 3 expected. \nVariant 1: Tight-Binding Object, Single-Shell Atom to be left unchanged, and optional Options. \nVariant 2: Operator, number of states to be left unchanged, and optional Options. \nVariant 3: Matrix Representation of Operator, number of states to be left unchanged, and optional Options.",len);

  if (LuaTestOptionsList(L, -1))
  			LuaCheckOptions(L, -1, LuaReduceHilbertSpaceOptions);

	if (luaL_isType(L, 1, tb))
	{
		Options.hasSpin = false;
		lua_check(L, 1, tb);
//		TODO: QTightBindingHermitian
		AtomType unchangedAtom;
		lua_check(L, 2, unchangedAtom);
		if (unchangedAtom.Nshell != 1)
			return luaL_error(L, "Error: RedceHilbertSpace only implemented for unchanged atoms with only one shell.");
		lua_push(L, reduceHilbertSpaceForTightBinding(tb, T, unchangedAtom, Options.filledLigands, Options.emptyLigands, Options.hasSpin));
		lua_push(L, T);
		return 2;
	} else if (luaL_isType(L, 1, op))
	{
		lua_check(L, 1, op);
		if (!(bool) QOperatorHermitian(op, 10 * DBL_EPSILON))
			return luaL_error(L, "Error: ReduceHilbertSpace expects a Hermitian operator as input.");
		uint unchangedStates;
		vector < string > orbNames;
		lua_check(L, 2, unchangedStates);
		lua_push(L, reduceHilbertSpaceForOperator(op, T, unchangedStates, Options.filledLigands, Options.emptyLigands, orbNames, Options.hasSpin));
		lua_push(L, T);
		lua_push(L, orbNames);
		return 3;
	} else if (luaL_isType(L, 1, mat))
	{
		lua_check(L, 1, mat);
		if (mat.Ni != mat.Nj)
			return luaL_error(L, "Error: ReduceHilbertSpace expects a quadratic matrix as input.");
		if (!(bool) QCompactMatrixHermitian(mat, 10 * DBL_EPSILON))
			return luaL_error(L, "Error: ReduceHilbertSpace expects a Hermitian matrix as input.");
		uint unchangedStates;
		vector < string > orbNames;
		lua_check(L, 2, unchangedStates);
		if (unchangedStates >= mat.Ni)
			return luaL_error(L, "Error: Number of unchanged states (%d) must be smaller than dimension of matrix (%d).", unchangedStates, mat.Ni);
		lua_push(L, reduceHilbertSpaceForMatrix(mat, T, unchangedStates, Options.filledLigands, Options.emptyLigands, orbNames, Options.hasSpin));
		lua_push(L, T);
		lua_push(L, orbNames);
		return 3;
	} else
	{
		return luaL_error(L, "Error: ReduceHilbertSpace only implemented for Tight-Binding Objects, Operators and Matrices.");
	}
}

