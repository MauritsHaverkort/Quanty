CompactMatrixType CorthogonalShells(const OrcaObject &OrcaMolec, const OrcaObject &OrcaAtom, const double ELower, const double EUpper, const vector<uint32_t> &iAtom, CompactMatrixType &newBasis)
{
	CompactMatrixType out, molecWindowMO, atomicShells, Overlap, atomicShellsOnMO, HMOMat;
	int iLower = 0;
	int iUpper = OrcaMolec.NBasis - 1;
	for (int i = 0; i < OrcaMolec.NBasis; ++i)
	{
		if (OrcaMolec.EigenEnergies[i] < ELower && i >= iLower)
		{
			iLower = i + 1;
		}
		if (OrcaMolec.EigenEnergies[i] > EUpper)
		{
			if (i == 0)
			{
				cout << "Error in CorthogonalShells: Upper energy bound " << EUpper << " eV is below groundstate energy " << OrcaMolec.EigenEnergies[0] << "." << endl;
				throw std::exception();
			}
			iUpper = i - 1;
			break;
		}
	}
	if (iLower > iUpper)
	{
		cout << "Error in CorthogonalShells: Lower index " << iLower << " is above upper index " << iUpper << ". Your energy window is probably too small. Did you use Hartree instead of eV?" << endl;
		throw std::exception();
	}

	vector<uint32_t> iMolec, iMolecFull;
	for (int i = iLower; i <= iUpper; ++i)
	{
		iMolec.push_back(i);
	}
	for (int i = 0; i < OrcaMolec.NBasis; ++i)
	{
		iMolecFull.push_back(i);
	}
	molecWindowMO = Submatrix(OrcaMolec.MolecularOrbitalMatrix, iMolec, iMolecFull);

	vector<uint32_t> iAtomFull;
	for (int i = 0; i < OrcaAtom.NBasis; ++i)
	{
		iAtomFull.push_back(i);
	}
	atomicShells = Submatrix(OrcaAtom.MolecularOrbitalMatrix, iAtom, iAtomFull);

	Overlap = Submatrix(OrcaMolec.OverlapMatrix, iMolecFull, iAtomFull);
	CompactMatrixMultiply3('N','T','T', &atomicShells, &Overlap, &molecWindowMO, &atomicShellsOnMO);
	CompactMatrixSVDOrthogonalizeRow(&atomicShellsOnMO, 1E-12, true);

	newBasis.Ni = atomicShellsOnMO.Ni + atomicShellsOnMO.Nj;
	newBasis.Nj = atomicShellsOnMO.Nj;
	newBasis.QComplex = 0;
	InitCompactMatrix(&newBasis);
	for (int i = 0; i < newBasis.Ni; ++i)
	{
		if (i < atomicShellsOnMO.Ni)
		{
			for (int j = 0; j < atomicShellsOnMO.Nj; ++j)
			{
				newBasis.VRe[i * newBasis.Nj + j] = atomicShellsOnMO.VRe[i * newBasis.Nj + j];
			}
		}
		else
		{
			int j = i - atomicShellsOnMO.Ni;
			newBasis.VRe[i * newBasis.Nj + j] = 1;
		}
	}
	CompactMatrixGramSchmidtOrthogonalizeRow(&newBasis, 1E-12, true);

	HMOMat.Ni = newBasis.Ni;
	HMOMat.Nj = newBasis.Nj;
	HMOMat.QComplex = 0;
	InitCompactMatrix(&HMOMat);
	for (int i = 0; i < HMOMat.Ni; ++i)
	{
		HMOMat.VRe[i * HMOMat.Nj + i] = OrcaMolec.EigenEnergies[iMolec[i]];
	}

	TransformCompactMatrixcT(&HMOMat, &newBasis, &out);

	CompactMatrixChop(&out, DBL_EPSILON);
	CompactMatrixChop(&newBasis, DBL_EPSILON);

	FreeCompactMatrix(&molecWindowMO);
	FreeCompactMatrix(&atomicShells);
	FreeCompactMatrix(&Overlap);
	FreeCompactMatrix(&atomicShellsOnMO);
	FreeCompactMatrix(&HMOMat);
	return out;
}


static int LuaCorthogonalShells(lua_State *L)
{
	OrcaObject OrcaMolec, OrcaAtom;
	double ELower, EUpper;
	vector < uint32_t > iAtom;
	CompactMatrixType newBasis;
	int len = lua_gettop(L);
	if (len < 5 || len > 6)
	{
		return luaL_error(L,
				"CorthogonalShells called with %d elements, 5 or 6 expected: OrcaOutput of a molecule, OrcaOutput of the first atom, lower energy bound [eV], upper energy bound [eV], and either a list of atomic indices or lower and upper atomic index.",
				len);
	}

	lua_check(L, 1, OrcaMolec);
	if (OrcaMolec.Atoms.size() == 0)
	{
		return luaL_error(L, "The molecular OrcaOutput contains no atoms, something has gone sincerely wrong.");
	}

	lua_check(L, 2, OrcaAtom);
	if (OrcaAtom.Atoms.size() == 0)
	{
		return luaL_error(L, "The atomic OrcaOutput contains no atoms, something has gone sincerely wrong.");
	}
	if (OrcaAtom.Atoms.size() > 1)
	{
		return luaL_error(L, "The atomic OrcaOutput is actually that of a molecule with %d atoms. Did you the molecular and the atomic input?", OrcaAtom.Atoms.size());
	}
	if (strcmp(OrcaAtom.Atoms[0].Name, OrcaMolec.Atoms[0].Name) != 0)
	{
		return luaL_error(L, "The atom in the atomic OrcaOutput must be the same as the first atom in the molecular OrcaOutput.");
	}

	lua_check(L,3,ELower);

	lua_check(L,4,EUpper);
	if (ELower > EUpper)
	{
		return luaL_error(L, "Lower energy bound (%f) is above upper energy bound (%f).", ELower, EUpper);
	}

	if (len == 5)
	{
		lua_check(L, 5, iAtom);
		for (int i = 0; i < iAtom.size(); ++i)
		{
			if (iAtom[i] < 1 || iAtom[i] > OrcaAtom.NBasis)
			{
				return luaL_error(L, "Atomic index %d lies outside range [1,%d].", iAtom[i], OrcaAtom.NBasis);
			}
			iAtom[i] = iAtom[i] - 1;
		}
	}
	else if (len == 6)
	{
		uint32_t iLower, iUpper;
		lua_check(L, 5, iLower);
		if (iLower < 1 || iLower > OrcaAtom.NBasis)
		{
			return luaL_error(L, "Lower atomic index %d lies outside range [1,%d].", iLower, OrcaAtom.NBasis);
		}

		lua_check(L, 6, iUpper);
		if (iUpper < 1 || iUpper > OrcaAtom.NBasis)
		{
			return luaL_error(L, "Upper atomic index %d lies outside range [1,%d].", iUpper, OrcaAtom.NBasis);
		}
		if (iLower > iUpper)
		{
			return luaL_error(L, "Lower atomic index (%d) is above upper atomic index (%d).", iLower, iUpper);
		}

		for (int i = iLower - 1; i < iUpper; ++i)
		{
			iAtom.push_back(i);
		}
	}

	lua_push(L, CorthogonalShells(OrcaMolec, OrcaAtom, ELower, EUpper, iAtom, newBasis));
	newBasis.OnLuaStack = true;
	lua_push(L, newBasis);
	return 2;
}
