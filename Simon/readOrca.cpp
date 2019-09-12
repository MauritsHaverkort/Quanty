class FileOpenException: public std::exception
{
public:
	string filename = "";
	FileOpenException(const string &filename)
	{
		this->filename = filename;
	}
	virtual ~FileOpenException() throw () //Promise the compiler that our destructor will not throw any exceptions
	{
	}
};

//void lua_check(lua_State * L, int arg, OrcaObject &in)
//{
//
//}

void lua_push(lua_State * L, const OrcaObject &pushkin)
{
	lua_newtable(L);
	int top = (int) lua_gettop(L);

	lua_push(L, "NBasis");
	lua_push(L, pushkin.NBasis);
	lua_settable(L, top);

	lua_push(L, "OverlapMatrix");
	lua_push(L, pushkin.OverlapMatrix);
	lua_settable(L, top);

	lua_push(L, "MolecularOrbitalMatrix");
	lua_push(L, pushkin.MolecularOrbitalMatrix);
	lua_settable(L, top);

	lua_push(L, "EigenEnergies");
	lua_push(L, pushkin.EigenEnergies);
	lua_settable(L, top);

	lua_push(L, "EigenOccupations");
	lua_push(L, pushkin.EigenOccupations);
	lua_settable(L, top);

	lua_push(L, "BasisOrbitalNames");
	lua_push(L, pushkin.BasisOrbitalNames);
	lua_settable(L, top);

	lua_push(L, "Atoms");
	lua_push(L, pushkin.Atoms);
	lua_settable(L, top);

	lua_push(L, "RadialFunctions");
	lua_push(L, pushkin.RadialFunctions);
	lua_settable(L, top);
}

bool isOrcaLineIgnored(const string &line)
{
	if ((uint)line.size() == 0) return true;
	if (line.find("# Basis set for element") != string::npos) return true;
	if (line.find("NewGTO") != string::npos) return true;
	if (line.find("CARTESIAN COORDINATES") != string::npos) return true;
	if (line.find("MASS") != string::npos) return true;
	if (line.find("OVERLAP MATRIX") != string::npos) return true;
	if (line.find("MOLECULAR ORBITALS") != string::npos) return true;
	if (line.find("---") != string::npos) return true;
	return false;
}

vector<double> stringToDoubleVec(const string &str)
{
	vector<double> out;
	std::stringstream strStr(str);
	while (true)
	{
		double entry;
		strStr >> entry;
		if (!strStr) break;
		out.push_back(entry);
	}
	return out;
}

int stringToLz(const string &str)
{
	if (str.empty())
	{
		return 0;
	}
	else if (str.compare("y") == 0)
	{
		return -1;
	}
	else if (str.compare("z") == 0)
	{
		return 0;
	}
	else if (str.compare("x") == 0)
	{
		return 1;
	}
	else if (str.compare("xy") == 0)
	{
		return -2;
	}
	else if (str.compare("yz") == 0)
	{
		return -1;
	}
	else if (str.compare("z2") == 0)
	{
		return 0;
	}
	else if (str.compare("xz") == 0)
	{
		return 1;
	}
	else if (str.compare("x2y2") == 0)
	{
		return 2;
	}
	else
	{
		return atoi(str.c_str());
	}
}

string constructOrbname(const int siteNum, const string &atomName, const uint n, const string &shell, const int lz, const bool addLz)
{
	//	Expects input like 2_Fe_3p1 or Fe_3p1
	std::stringstream strStr;
	if (siteNum >= 0)
	{
		strStr << siteNum << "_";
	}
	strStr << atomName << "_" << n << shell;
	if (addLz)
	{
		strStr << lz;
	}
	return strStr.str();
}

string constructOrbname(const int siteNum, const string &atomName, const uint n, const string &shell, const int lz)
{
	return constructOrbname(siteNum, atomName, n, shell, lz, true);
}

string constructOrbname(const int siteNum, const string &atomName, const uint n, const string &shell)
{
	return constructOrbname(siteNum, atomName, n, shell, 0, false);
}

string constructOrbname(const string &atomName, const uint n, const string &shell, const int lz)
{
	return constructOrbname(-1, atomName, n, shell, lz, true);
}

string constructOrbname(const string &atomName, const uint n, const string &shell)
{
	return constructOrbname(-1, atomName, n, shell, 0, false);
}

string normalizedOrcaOrbname(const string &orbOrca)
{
	const static string digits = "0123456789";
	//	Expects input like " 0Cr 3dxy "
	uint siteNum;
	string atomName;
	uint n;
	string shell;
	int lz;

	string::size_type lastSiteNumPos = orbOrca.find_first_of(digits);
	lastSiteNumPos =  orbOrca.find_first_not_of(digits, lastSiteNumPos) - 1;
	siteNum = atoi(orbOrca.substr(0, lastSiteNumPos + 1).c_str());

	string::size_type lastAtomNamePos = orbOrca.find_first_of(" ", lastSiteNumPos + 1) - 1;
	atomName = orbOrca.substr(lastSiteNumPos + 1, lastAtomNamePos - lastSiteNumPos);

	string::size_type lastNPos = orbOrca.find_first_of(digits, lastAtomNamePos);
	lastNPos =  orbOrca.find_first_not_of(digits, lastNPos) - 1;
	n = atoi(orbOrca.substr(lastAtomNamePos + 1, lastNPos - lastAtomNamePos).c_str());

	shell = orbOrca.substr(lastNPos + 1, 1);
	uint l = (orbToMultiplicity(shell) / 2 - 1) / 2;

	string lzName = orbOrca.substr(lastNPos + 2);
	if (lzName.find(" ") != string::npos) lzName = lzName.substr(0, lzName.find(" "));
	lz = stringToLz(lzName);

	return constructOrbname(siteNum, atomName, n, shell, lz);
}

void SortOrcaBasis(OrcaObject &orca)
{
	CompactMatrixType tmp;
	vector < string > oldOrbs = orca.BasisOrbitalNames;
	sort(orca.BasisOrbitalNames.begin(), orca.BasisOrbitalNames.end(), atomicStateSorterlnz);
	vector<int> permVec = permutationVector(oldOrbs, orca.BasisOrbitalNames);
	CompactMatrixType permMat;
	permMat.Ni = orca.NBasis;
	permMat.Nj = orca.NBasis;
	permMat.QComplex = 0;
	InitCompactMatrix(&permMat);
	for (uint i = 0; i < (uint) permVec.size(); ++i)
	{
		permMat.VRe[i * permMat.Nj + permVec[i]] = 1;
	}
	CompactMatrixMultiply('N', 'T', &orca.MolecularOrbitalMatrix, &permMat, &tmp);
	bool onStack = orca.MolecularOrbitalMatrix.OnLuaStack;
	orca.MolecularOrbitalMatrix.OnLuaStack = false;
	FreeCompactMatrix(&orca.MolecularOrbitalMatrix);
	orca.MolecularOrbitalMatrix = tmp;
	orca.MolecularOrbitalMatrix.OnLuaStack = onStack;

	TransformCompactMatrixcT(&orca.OverlapMatrix, &permMat, &tmp);
	onStack = orca.OverlapMatrix.OnLuaStack;
	orca.OverlapMatrix.OnLuaStack = false;
	FreeCompactMatrix(&orca.OverlapMatrix);
	orca.OverlapMatrix = tmp;
	orca.OverlapMatrix.OnLuaStack = onStack;
	FreeCompactMatrix(&permMat);
}

double OrcaGaussianNormfactor(const double w, const uint l)
{
	double out = sqrt(4 * M_PI);
	out *= pow(2 * w / M_PI, 3. / 4.);
	out *= pow(4 * w, (double) l / 2.);
	if (l > 0)
	{
		uint den = 2 * l - 1;
		uint tmp = den;
		while (tmp > 1)
		{
			tmp -= 2;
			den *= tmp;
		}
		den *= 2 * l + 1;
		out /= sqrt((double) den);
	}
	return out;
}

OrcaObject FileReadOrca(const string &filename)
{
	const static double HartreeToeV = 27.211391239;
	const static string letters = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
	OrcaObject out;
	ifstream reader;
	string line;
	reader.open(filename.c_str());
	if (!reader.is_open()) throw FileOpenException(filename);

	bool isReadingAtoms = false;
	bool isReadingRadialFunctions = false;
	bool isReadingOverlapMatrix = false;
	bool isReadingMolecularOrbitals = false;
	vector < vector<double> > atomPositions;
	struct
	{
	string currentAtom;
	char currentShell;
	uint shellCounter;
	GaussianRadial currentGauss;
	} basisAux;
	int i, jOffset = 0;
	while (getline(reader, line))
	{
		if (line.find("# of contracted basis functions") != string::npos)
		{
			const string::size_type lastBlank = line.find_last_of(" ");
			out.NBasis = atoi(line.substr(lastBlank + 1).c_str());
		}

		if (line.find("CARTESIAN COORDINATES (A.U.)") != string::npos)
		{
			isReadingAtoms = true;
		}

		if (isReadingAtoms && line.find("INTERNAL COORDINATES") != string::npos)
		{
			isReadingAtoms = false;
		}

		if (isReadingAtoms && !isOrcaLineIgnored(line))
		{
			string::size_type pos = line.find_first_of(letters);
			pos = line.find_first_not_of(letters, pos);
			vector<double> lineVec = stringToDoubleVec(line.substr(pos)); //ZA, FRAG, MASS, X, Y, Z
			vector<double> atomPos(3);
			for (uint i = 0; i < 3; ++i)
			{
				atomPos[i] = lineVec[lineVec.size()-3+i];
			}
			atomPositions.push_back(atomPos);
		}

		if (line.find("# Basis set for element :") != string::npos)
		{
			isReadingRadialFunctions = true;
			string::size_type colonPos = line.find_last_of(':');
			string::size_type atomFirstPos = line.find_first_not_of(' ', colonPos + 1);
			string::size_type atomLastPos = line.find_first_not_of(letters, atomFirstPos);
			basisAux.currentAtom = line.substr(atomFirstPos, atomLastPos - atomFirstPos);
			basisAux.currentShell = 'x';
			basisAux.shellCounter = 0;
		}

		if (isReadingRadialFunctions && !isOrcaLineIgnored(line))
		{
			if (line.find_first_of(letters) != string::npos)
			{
				if (basisAux.shellCounter > 0)
				{
					std::stringstream strstr;
					strstr << basisAux.currentAtom << "_" << basisAux.shellCounter << basisAux.currentShell;
					const string shellNameFull = strstr.str();
					out.RadialFunctions.insert( { shellNameFull, basisAux.currentGauss });
				}
				if (line.find("end") != string::npos)
				{
					isReadingRadialFunctions = false;
				}
				else
				{
					string::size_type shellFirstPos = line.find_first_of(letters);
					if (basisAux.currentShell != tolower(line.at(shellFirstPos)))
					{
						basisAux.shellCounter = 1;
						basisAux.currentShell = tolower(line.at(shellFirstPos));
					}
					else
					{
						basisAux.shellCounter++;
					}
					uint l = 0;
					if (basisAux.currentShell == 's') l = 0;
					if (basisAux.currentShell == 'p') l = 1;
					if (basisAux.currentShell == 'd') l = 2;
					if (basisAux.currentShell == 'f') l = 3;
					if (basisAux.currentShell == 'g') l = 4;
					if (basisAux.currentShell == 'h') l = 5;
					basisAux.currentGauss = GaussianRadial(l);
				}
			}
			else
			{
				vector<double> lineVec = stringToDoubleVec(line);
				basisAux.currentGauss.addGaussian(lineVec[2] * OrcaGaussianNormfactor(lineVec[1], basisAux.currentGauss.getL()), lineVec[1]);
			}
		}

		if (line.find("OVERLAP MATRIX") != string::npos)
		{
			isReadingOverlapMatrix = true;
			out.OverlapMatrix.Ni = out.NBasis;
			out.OverlapMatrix.Nj = out.NBasis;
			out.OverlapMatrix.QComplex = 0;
			InitCompactMatrix(&(out.OverlapMatrix));
			out.OverlapMatrix.OnLuaStack = true;
		}

		if (isReadingOverlapMatrix)
		{
			const string firstCol = line.substr(0, 10);
			if (firstCol.compare(blankFill(10)) == 0)
			{
				const string::size_type firstNumFirstPos = line.find_first_not_of(" ");
				const string::size_type firstNumLastPos = line.find_first_of(" ", firstNumFirstPos);
				jOffset = atoi(line.substr(firstNumFirstPos, firstNumLastPos - firstNumFirstPos).c_str());
			}
			else if (!isOrcaLineIgnored(line))
			{
				vector<double> lineVec = stringToDoubleVec(line);
				i = (int) lineVec[0];
				int j = jOffset - 1;
				for (uint k = 1; k < (uint) lineVec.size(); ++k)
				{
					++j;
					out.OverlapMatrix.VRe[i * out.OverlapMatrix.Nj + j] = lineVec[k];
				}
				if (i == out.OverlapMatrix.Ni - 1 && j == out.OverlapMatrix.Nj - 1) isReadingOverlapMatrix = false;
			}
		}

		if (line.find("MOLECULAR ORBITALS") != string::npos)
		{
			isReadingMolecularOrbitals = true;
			i = -3;
			out.MolecularOrbitalMatrix.Ni = out.NBasis;
			out.MolecularOrbitalMatrix.Nj = out.NBasis;
			out.MolecularOrbitalMatrix.QComplex = 0;
			InitCompactMatrix(&(out.MolecularOrbitalMatrix));
			out.MolecularOrbitalMatrix.OnLuaStack = true;
		}

		if (isReadingMolecularOrbitals && !isOrcaLineIgnored(line))
		{
			if (i == -3)
			{
				std::stringstream strStr(line);
				strStr >> jOffset;
			}
			else if (i == -2)
			{
				vector<double> energies = stringToDoubleVec(line);
				for (uint k = 0; k < (uint) energies.size(); ++k)
				{
					out.EigenEnergies.push_back(energies[k] * HartreeToeV);
				}
			}
			else if (i == -1)
			{
				vector<double> occus = stringToDoubleVec(line);
				for (uint k = 0; k < (uint) occus.size(); ++k)
				{
					out.EigenOccupations.push_back(occus[k]);
				}
			}
			else
			{
				string::size_type orbNameLength = line.find_first_of(" ", 12);
				if ((uint) out.BasisOrbitalNames.size() < out.NBasis)
				{
					const string orbName = normalizedOrcaOrbname(line.substr(0, orbNameLength));
					out.BasisOrbitalNames.push_back(orbName);
				}

				vector<double> entries = stringToDoubleVec(line.substr(orbNameLength));
				double phaseFac = OrcaPhaseConventionFactor(out.BasisOrbitalNames[i]);
				for (uint k = 0; k < (uint) entries.size(); ++k)
				{
					uint j = jOffset + k;
					out.MolecularOrbitalMatrix.VRe[j * out.MolecularOrbitalMatrix.Nj + i] = phaseFac * entries[k]; //we want the matrix to contain row vectors.
				}
			}

			if (i == out.NBasis - 1)
			{
				if ((uint) out.EigenEnergies.size() == out.NBasis)
				{
					isReadingMolecularOrbitals = false;
				}
				else
				{
					i = -3;
				}
			}
			else
			{
				i++;
			}
		}
	}
	reader.close();
	ChangeOrcaOverlapMatrixPhases(out);
	SortOrcaBasis(out);
	out.Atoms = OrbsToAtoms(out.BasisOrbitalNames, false);
	for (uint i = 0; i < (uint) atomPositions.size(); ++i)
	{
		for (uint j = 0; j < 3; ++j)
		{
			out.Atoms[i].Position[j] = atomPositions[i][j];
		}
	}
	return out;
}

void ChangeOrcaOverlapMatrixPhases(OrcaObject &orca)
{
	for (uint i = 0; i < orca.NBasis; ++i)
	{
		for (uint j = 0; j < orca.NBasis; ++j)
		{
			const double faci = OrcaPhaseConventionFactor(orca.BasisOrbitalNames[i]);
			const double facj = OrcaPhaseConventionFactor(orca.BasisOrbitalNames[j]);
			orca.OverlapMatrix.VRe[i * orca.NBasis + j] = orca.OverlapMatrix.VRe[i * orca.NBasis + j] * faci * facj;
		}
	}
}

double OrcaPhaseConventionFactor(const string &orb)
{
	const static string digitsEtc = "+-0123456789";
	const string::size_type pos = orb.find_last_not_of(digitsEtc);
	const string llz = orb.substr(pos);
	if (llz.compare("f-3") == 0)
	{
		return -1.;
	}
	else if (llz.compare("f3") == 0)
	{
		return -1.;
	}
	else if (llz.compare("g-4") == 0)
	{
		return -1. / sqrt(3);
	}
	else if (llz.compare("g-3") == 0)
	{
		return -1. / sqrt(3);
	}
	else if (llz.compare("g-2") == 0)
	{
		return 1. / sqrt(3);
	}
	else if (llz.compare("g-1") == 0)
	{
		return 1. / sqrt(3);
	}
	else if (llz.compare("g0") == 0)
	{
		return 1. / sqrt(3);
	}
	else if (llz.compare("g1") == 0)
	{
		return 1. / sqrt(3);
	}
	else if (llz.compare("g2") == 0)
	{
		return 1. / sqrt(3);
	}
	else if (llz.compare("g3") == 0)
	{
		return -1. / sqrt(3);
	}
	else if (llz.compare("g4") == 0)
	{
		return -1. / sqrt(3);
	}
	else
	{
		return 1.;
	}
}

void OrcaDiagonalizeDegenerateStates(OrcaObject &orca)
{

	CompactMatrixType TrafoFull;
	const uint NFull = orca.NBasis;
	TrafoFull.Ni = NFull;
	TrafoFull.Nj = NFull;
	TrafoFull.QComplex = 0;
	TrafoFull.OnLuaStack = false;
	InitCompactMatrix(&TrafoFull);
	vector<uint> alreadyRotated;
	for (uint k = 0; k < NFull; ++k)
	{
		if (!isElementOf(k, alreadyRotated))
		{
			vector<uint> degOrbIndices;
			degOrbIndices.push_back(k);
			alreadyRotated.push_back(k);

			for (uint l = k + 1; l < NFull; ++l)
			{
				if (fabs(orca.EigenEnergies[k] - orca.EigenEnergies[l]) < 1e-10 && !isElementOf(l, alreadyRotated))
				{

					degOrbIndices.push_back(l);
					alreadyRotated.push_back(l);
				}
			}

			vector<uint> largestEntries;
			for (uint l = 0; l < (uint) degOrbIndices.size(); ++l)
			{
				const uint MOindex = degOrbIndices[l];
				uint iMax = 0;
				double partMax = 0.;
				for (uint m = 0; m < NFull; ++m)
				{
					double part = fabs(orca.MolecularOrbitalMatrix.VRe[MOindex * NFull + m]);
					if (part > partMax && !isElementOf(m, largestEntries))
					{
						iMax = m;
						partMax = part;
					}
				}
				largestEntries.push_back(iMax);
			}

			sort(largestEntries.begin(), largestEntries.end());
			CompactMatrixType subMat;
			subMat.Ni = (uint) degOrbIndices.size();
			subMat.Nj = (uint) degOrbIndices.size();
			subMat.QComplex = 0;
			InitCompactMatrix(&subMat);

			for (uint i = 0; i < (uint) degOrbIndices.size(); ++i)
			{
				const uint iFull = degOrbIndices[i];
				for (uint j = 0; j < (uint) largestEntries.size(); ++j)
				{
					const uint jFull = largestEntries[j];
					subMat.VRe[i * subMat.Nj + j] = orca.MolecularOrbitalMatrix.VRe[iFull * NFull + jFull];
				}
			}

			CompactMatrixType Trafo = GaussEliminationTransformation(subMat);
			for (uint i = 0; i < (uint) degOrbIndices.size(); ++i)
			{
				const uint iFull = degOrbIndices[i];
				for (uint j = 0; j < (uint) degOrbIndices.size(); ++j)
				{
					const uint jFull = degOrbIndices[j];
					TrafoFull.VRe[iFull * TrafoFull.Nj + jFull] = Trafo.VRe[i * Trafo.Nj + j];
				}
			}
			FreeCompactMatrix(&Trafo);
			FreeCompactMatrix(&subMat);
		}
	}

	CompactMatrixType tmpMat;
	CompactMatrixMultiply('N', 'N', &TrafoFull, &orca.MolecularOrbitalMatrix, &tmpMat);
	bool onStack = orca.MolecularOrbitalMatrix.OnLuaStack;
	orca.MolecularOrbitalMatrix.OnLuaStack = false;
	FreeCompactMatrix(&orca.MolecularOrbitalMatrix);
	orca.MolecularOrbitalMatrix = tmpMat;
	orca.MolecularOrbitalMatrix.OnLuaStack = onStack;

	CompactMatrixChop(&orca.MolecularOrbitalMatrix, 1e-14);
	FreeCompactMatrix(&TrafoFull);
}

static int LuaFileReadOrca(lua_State* L)
{
	uint32_t len = (uint32_t) lua_gettop(L);
	if (len < 1 || len > 2) return luaL_error(L, "ReadOrca called with %d arguments, 1 or 2 expected, the path of the Orca output file and optional Options.", len);
	string filename;
  struct
  {
		bool diagonalizeDegenerateStates = false;
  } Options;
  LuaOptionsType LuaFileReadOrcaOptions[]={
		{"diagonalizeDegenerateStates",   "bool", &(Options.diagonalizeDegenerateStates)},
    {NULL, NULL, NULL}};

	lua_check(L, 1, filename);
	if (LuaTestOptionsList(L, -1)) LuaCheckOptions(L, -1, LuaFileReadOrcaOptions);

	try
	{
		OrcaObject Orcaout = FileReadOrca(filename);
		if (Options.diagonalizeDegenerateStates)
		{
			OrcaDiagonalizeDegenerateStates(Orcaout);
		}
		lua_push(L, Orcaout);
		return 1;
	}
	catch (FileOpenException &e)
	{
		return luaL_error(L, "Could not open Orca output file with path %s.", e.filename.c_str());
	}
	catch (ifstream::failure &e)
	{
		return luaL_error(L, "Something in FileReadOrca went wrong while reading through the file %s.", filename.c_str());
	}
	catch (std::exception &e)
	{
		cout << "Exception of type " << e.what() << endl;
		return luaL_error(L, "Something in FileReadOrca went unexpectedly wrong.");
	}
}

CompactMatrixType GaussEliminationTransformation(const CompactMatrixType &mat)
{
	CompactMatrixType out, inCopy, tmp;

	if (mat.Ni != mat.Nj)
	{
		cout << "Error: GaussEliminationTransformation only impelemented for quadratic matrices." << endl;
		throw std::exception();
	}

	out.Ni = mat.Ni;
	out.Nj = mat.Nj;
	out.QComplex = mat.QComplex;
	InitCompactIdentityMatrix(&out);
	out.OnLuaStack = false;
	CompactMatrixCopy(mat, &inCopy);
	inCopy.OnLuaStack = false;

	if (inCopy.VRe[0] == 0)
	{
		uint i;
		for (i = 1; i < inCopy.Ni; ++i)
		{
			if (inCopy.VRe[i * inCopy.Nj] != 0) break;
		}

		if (i == inCopy.Ni)
		{
			cout << "Error: GaussEliminationTransformation called with a zero-column." << endl;
			throw std::exception();
		}

		CompactMatrixType perm;
		perm.Ni = inCopy.Ni;
		perm.Nj = inCopy.Nj;
		perm.QComplex = 0;
		InitCompactIdentityMatrix(&perm);

		perm.VRe[0] = 0;
		perm.VRe[0 * inCopy.Nj + i] = 1;
		perm.VRe[i * inCopy.Nj + 0] = 1;
		perm.VRe[i * inCopy.Nj + i] = 0;

		CompactMatrixMultiply('N', 'N', &perm, &out, &tmp);
		FreeCompactMatrix(&out);
		out = tmp;

		CompactMatrixMultiply('N', 'N', &perm, &inCopy, &tmp);
		FreeCompactMatrix(&inCopy);
		inCopy = tmp;

		FreeCompactMatrix(&perm);
	}

	if (out.QComplex == 0)
	{
		for (uint j = 0; j < mat.Nj - 1; ++j)
		{
			for (uint i = j + 1; i < mat.Ni; ++i)
			{
				if (inCopy.VRe[j * inCopy.Nj + j] == 0)
				{
					cout << "Error in GaussEliminationTransformation: Matrix rows are not linearly independent." << endl;
					throw std::exception();
				}
				double fac = -inCopy.VRe[i * inCopy.Nj + j] / inCopy.VRe[j * inCopy.Nj + j];
				CompactMatrixType TStep;
				TStep.Ni = out.Ni;
				TStep.Nj = out.Nj;
				TStep.QComplex = out.QComplex;
				InitCompactIdentityMatrix(&TStep);
				TStep.VRe[i * TStep.Nj + j] = fac;

				CompactMatrixMultiply('N', 'N', &TStep, &out, &tmp);
				FreeCompactMatrix(&out);
				out = tmp;

				CompactMatrixMultiply('N', 'N', &TStep, &inCopy, &tmp);
				FreeCompactMatrix(&inCopy);
				inCopy = tmp;

				FreeCompactMatrix(&TStep);
			}
		}
		if (inCopy.VRe[(inCopy.Ni - 1) * inCopy.Nj + (inCopy.Nj - 1)] == 0)
		{
			cout << "Error in GaussEliminationTransformation: Matrix rows are not linearly independent." << endl;
			throw std::exception();
		}

		for (int j = mat.Nj - 1; j > 0; --j)
		{
			for (int i = j - 1; i >= 0; --i)
			{
				double fac = -inCopy.VRe[i * inCopy.Nj + j] / inCopy.VRe[j * inCopy.Nj + j];
				CompactMatrixType TStep;
				TStep.Ni = out.Ni;
				TStep.Nj = out.Nj;
				TStep.QComplex = out.QComplex;
				InitCompactIdentityMatrix(&TStep);
				TStep.VRe[i * TStep.Nj + j] = fac;

				CompactMatrixMultiply('N', 'N', &TStep, &out, &tmp);
				FreeCompactMatrix(&out);
				out = tmp;

				CompactMatrixMultiply('N', 'N', &TStep, &inCopy, &tmp);
				FreeCompactMatrix(&inCopy);
				inCopy = tmp;

				FreeCompactMatrix(&TStep);
			}
		}

		for (uint i = 0; i < out.Ni; ++i)
		{
			double newNorm, oldNorm = 0;
			newNorm = inCopy.VRe[i * inCopy.Nj + i]; //no fabs, in order to make all signs positive.
			for (uint j = 0; j < out.Nj; ++j)
			{
				double entry = mat.VRe[i * mat.Nj + j];
				oldNorm += entry * entry;
			}
			oldNorm = sqrt(oldNorm);
			for (uint j = 0; j < out.Nj; ++j)
			{
				out.VRe[i * out.Nj + j] = out.VRe[i * out.Nj + j] * oldNorm / newNorm;
			}
		}
	}
	else
	{
		cout << "Error: GaussEliminationTransformation not yet implemented for complex matrices." << endl;
		throw std::exception();
	}
	FreeCompactMatrix(&inCopy);
	return out;
}

