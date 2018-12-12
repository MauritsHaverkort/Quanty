
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

string multiplicityToNonRelOrb(const int mult)
{
	if (mult < 0)
		return "Error";
	switch (mult)
	{
	case 2:
		return "s";
	case 6:
		return "p";
	case 10:
		return "d";
	case 14:
		return "f";
	case 18:
		return "g";
	case 20:
		return "h";
	default:
		return to_string(mult) + "a";
	}
}

void addLigandNames(const CompactMatrixType &M, const string &nameBase, vector<string> &names, const uint unchangedStates)
{
	uint sameSymLigands = (uint) (M.Ni / unchangedStates) - 1;
	uint i;
	for (i = 0; i < sameSymLigands; ++i)
		names.push_back(nameBase + "_" + to_string(i + 1) + "_" + multiplicityToNonRelOrb(2 * unchangedStates));
	uint restLigandMult = M.Ni % unchangedStates;
	if (restLigandMult > 0)
		names.push_back(nameBase + "_" + to_string(i + 1) + "_" + multiplicityToNonRelOrb(2 * restLigandMult));
}

CompactMatrixType reduceHilbertSpaceForMatrix(const CompactMatrixType &H, CompactMatrixType &T, const uint unchangedStates, const uint filledLigands, const uint emptyLigands,
		vector<string> orbNames)
{
	CompactMatrixType out, HLigand, eigenfunsLigand, HRedFilled, Tfilled, HRedEmpty, Tempty;
	uint LigandMatSize = H.Ni - unchangedStates;
	orbNames = vector<string>();
	orbNames.push_back("unchanged_" + multiplicityToNonRelOrb(unchangedStates));
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
		addLigandNames(HRedFilled, "L_val", orbNames, unchangedStates);
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
		addLigandNames(HRedEmpty, "L_con", orbNames, unchangedStates);
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

CompactMatrixType reduceHilbertSpaceForMatrix(const CompactMatrixType &H, CompactMatrixType &T, const uint unchangedStates, const uint filledLigands, const uint emptyLigands)
{
	return reduceHilbertSpaceForMatrix(H, T, unchangedStates, filledLigands, emptyLigands, vector<string>());
}

CompactMatrixType reduceHilbertSpaceForMatrix(const CompactMatrixType &H, const uint unchangedStates, const uint filledLigands, const uint emptyLigands)
{
	CompactMatrixType T;
	CompactMatrixType out = reduceHilbertSpaceForMatrix(H, T, unchangedStates, filledLigands, emptyLigands, vector<string>());
	FreeCompactMatrix(&T);
	return out;
}

OperatorType reduceHilbertSpaceForOperator(const OperatorType &hamil, CompactMatrixType &T, const uint unchangedStates, const uint filledLigands, const uint emptyLigands, vector<string> orbNames)
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
	HMatRed = reduceHilbertSpaceForMatrix(HMat, T, unchangedStates, filledLigands, emptyLigands, orbNames);
	out = CompactMatrixToOperator(HMatRed);
	FreeCompactMatrix(&Rho);
	FreeCompactMatrix(&HMat);
	FreeCompactMatrix(&HMatRed);
	return out;
}

OperatorType reduceHilbertSpaceForOperator(const OperatorType &H, CompactMatrixType &T, const uint unchangedStates, const uint filledLigands, const uint emptyLigands)
{
	return reduceHilbertSpaceForOperator(H, T, unchangedStates, filledLigands, emptyLigands, vector<string>());
}

OperatorType reduceHilbertSpaceForOperator(const OperatorType &H, const uint unchangedStates, const uint filledLigands, const uint emptyLigands)
{
	CompactMatrixType T;
	OperatorType out = reduceHilbertSpaceForOperator(H, T, unchangedStates, filledLigands, emptyLigands, vector<string>());
	FreeCompactMatrix(&T);
	return out;
}


template<typename T>
bool isElementOf(const T element, const vector<T> &set, bool (*eq)(const T&, const T&))
{
	for (uint i = 0; i < (uint) set.size(); ++i)
	{
		if (eq(element, set[i]))
			return true;
	}
	return false;
}

//vector<string> AtomsToOrbs(const vector<AtomType> &atoms)
//{
//	vector < string > out;
//	for (uint i = 0; i < (uint) atoms.size(); ++i)
//		for (uint j = 0; j < atoms[i].Nshell; ++j)
//			out.push_back(to_string(atoms[i].Name) + "_" + to_string(atoms[i].ShellName[j]));
//	return out;
//}


bool QAtomsEqual(const AtomType &a1, const AtomType &a2)
{
	if (strcmp(a1.Name, a2.Name) != 0)
		return false;
	for (uint i = 0; i < 3; ++i)
		if (a1.Position[i] != a2.Position[i])
			return false;
	return true;
}



void sortTbAtomToFirstPosition(TightBindingDef &tb, const AtomType &atom)
{
	MatrixType permMat;
	vector<string> oldOrbs, newOrbs;
	for (uint i = 0; i < (uint) tb.Cell.NAtoms; ++i)
		for (uint j = 0; j < tb.Cell.Atom[i].Nshell; ++j)
		{
			oldOrbs.push_back(to_string(i + 1) + "_" + to_string(tb.Cell.Atom[i].Name) + "_" + to_string(tb.Cell.Atom[i].ShellName[j]));
			if (QAtomsEqual(tb.Cell.Atom[i], atom))
				newOrbs.push_back(oldOrbs[oldOrbs.size() - 1]);
		}
	for (uint i = 0; i < (uint) oldOrbs.size(); ++i)
		if (!isElementOf(oldOrbs[i], newOrbs))
			newOrbs.push_back(oldOrbs[i]);
	permMat = orbitalPermutationMatrix(oldOrbs, newOrbs);

//	TODO:
	debug(oldOrbs);
	debug(newOrbs);

	debug(permMat);
//	debug(atom.Name);
//	debug(atom.AtomNumber);
}

TightBindingDef reduceHilbertSpaceForTightBinding(TightBindingDef &tb, const CompactMatrixType &T, const AtomType &unchangedAtom, const uint filledLigands, const uint emptyLigands)
{
	TightBindingDef out = TightBindingDef();
	MatrixType hMat;
	CompactMatrixType hMatCompact;
//	TODO: Rotate unchanged atoms to beginning
	if (!QAtomsEqual(unchangedAtom, tb.Cell.Atom[0]))
		sortTbAtomToFirstPosition(tb, unchangedAtom);
	double k[3] = {0.,0.,0.};
	CreateHamiltonianFromTightBindingDefConstantPhaseInUnitCell(tb, k, &hMat);
	hMatCompact = MatrixTypeToCompactMatrixType(hMat);



//TODO
	cout << "reduceHilbertSpace for TightBinding not yet fully implemented." << endl;
	FreeMatrix(&hMat);
	FreeCompactMatrix(&hMatCompact);
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
		AtomType unchangedAtom;
		lua_check(L, 2, unchangedAtom);
		lua_push(L, reduceHilbertSpaceForTightBinding(tb, T, unchangedAtom, Options.filledLigands, Options.emptyLigands));
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
		lua_push(L, reduceHilbertSpaceForOperator(op, T, unchangedStates, Options.filledLigands, Options.emptyLigands, orbNames));
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
		lua_push(L, reduceHilbertSpaceForMatrix(mat, T, unchangedStates, Options.filledLigands, Options.emptyLigands, orbNames));
		lua_push(L, T);
		lua_push(L, orbNames);
		return 3;
	} else
	{
		return luaL_error(L, "Error: ReduceHilbertSpace only implemented for Tight-Binding Objects, Operators and Matrices.");
	}
}

