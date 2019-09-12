 
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

bool QAtomsEqual(const AtomType &a1, const AtomType &a2)
{
	if (strcmp(a1.Name, a2.Name) != 0)
		return false;
	for (uint i = 0; i < 3; ++i)
		if (a1.Position[i] != a2.Position[i])
			return false;
	return true;
}

string multiplicityToNonRelOrb(const int mult, const bool hasSpin)
{
	if (mult < 0)
		return "Error";
	if (hasSpin)
	{
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
		case 22:
			return "h";
		default:
			return to_string(mult) + "a";
		}
	} else
	{
		switch (mult)
		{
		case 1:
			return "s";
		case 3:
			return "p";
		case 5:
			return "d";
		case 7:
			return "f";
		case 9:
			return "g";
		case 11:
			return "h";
		default:
			return to_string(mult) + "a";
		}
	}
}

string multiplicityToNonRelOrb(const int mult)
{
	return multiplicityToNonRelOrb(mult, true);
}

TightBindingDef sortTbAtomToFirstPosition(TightBindingDef &tb, const AtomType &atom, CompactMatrixType &R)
{
	MatrixType permMat;
	TightBindingDef out;
	vector < AtomType > newAtoms;
	newAtoms.push_back(atom);
	vector<string> oldOrbs, newOrbs;
	for (uint i = 0; i < (uint) tb.Cell.NAtoms; ++i)
	{
		if (!QAtomsEqual(tb.Cell.Atom[i], atom))
			newAtoms.push_back(tb.Cell.Atom[i]);
		for (uint j = 0; j < tb.Cell.Atom[i].Nshell; ++j)
		{
			oldOrbs.push_back(to_string(i + 1) + "_" + to_string(tb.Cell.Atom[i].Name) + "_" + to_string(tb.Cell.Atom[i].ShellName[j]));
			if (QAtomsEqual(tb.Cell.Atom[i], atom))
				newOrbs.push_back(oldOrbs[oldOrbs.size() - 1]);
		}
	}
	for (uint i = 0; i < (uint) oldOrbs.size(); ++i)
		if (!isElementOf(oldOrbs[i], newOrbs))
			newOrbs.push_back(oldOrbs[i]);
	permMat = orbitalPermutationMatrix(oldOrbs, newOrbs, false);
	R = MatrixTypeToCompactMatrixType(permMat);
	TightBindingReadyForRotate(tb, out, newAtoms);
	RotateTightBinding(&tb, &R, &out);
	FreeMatrix(&permMat);
	return out;
}

void initAtom(AtomType &a, const vector<string> &shellNames, const vector<vector<string> > &orbitNames, const bool hasSpin)
{
	a.AtomNumber = AtomNameToNumber(a.Name);
	for (uint i = 0; i < 3; ++i)
		a.Position[i] = 0.;
	a.Nshell = (int) shellNames.size();
	if (a.Nshell != 0)
	{
		a.ShellName = (char**) malloc(a.Nshell * sizeof(char*));
		a.shelldegeneracy = (int*) malloc(a.Nshell * sizeof(int));
		a.OrbitalName = (char***) malloc(a.Nshell * sizeof(char**));
		for (int i = 0; i < a.Nshell; ++i)
		{
			a.ShellName[i] = (char*) malloc(MAXStrLength * sizeof(char));
			strcpy(a.ShellName[i], shellNames[i].c_str());
			if ((uint) orbitNames.size() != 0)
			{
				a.shelldegeneracy[i] = (int) orbitNames[i].size();
				if (a.shelldegeneracy[i] != 0)
				{
					a.OrbitalName[i] = (char**) malloc(a.shelldegeneracy[i] * sizeof(char*));
					for (uint j = 0; j < a.shelldegeneracy[i]; ++j)
					{
						a.OrbitalName[i][j] = (char*) malloc(MAXStrLength * sizeof(char));
						strcpy(a.OrbitalName[i][j], orbitNames[i][j].c_str());
					}
				}
			} else
			{
				a.shelldegeneracy[i] = (int) orbToMultiplicity(to_string(a.ShellName[i]), hasSpin);
				if (a.shelldegeneracy[i] > 0)
				{
					a.OrbitalName[i] = (char**) malloc(a.shelldegeneracy[i] * sizeof(char*));
					for (uint j = 0; j < a.shelldegeneracy[i]; ++j)
					{
						a.OrbitalName[i][j] = (char*) malloc(MAXStrLength * sizeof(char));
						if (hasSpin)
						{
							string orbName;
							if (j % 2 == 0)
								orbName = to_string(-(double) (a.shelldegeneracy[i] / 2 - 1) / 2 + (int) j / 2) + "_Dn";
							else
								orbName = to_string(-(double) (a.shelldegeneracy[i] / 2 - 1) / 2 + (int) j / 2) + "_Up";
							strcpy(a.OrbitalName[i][j], orbName.c_str());
						} else
						{
							strcpy(a.OrbitalName[i][j], to_string(-(double) (a.shelldegeneracy[i] - 1) / 2 + (double) j).c_str());
						}
					}
				}
			}
		}
	}
}

vector<string> orbsAddLz(const vector<string> &in)
{
	vector < string > out;
	for (uint i = 0; i < (uint) in.size(); ++i)
	{
		int l = (orbToMultiplicity(in[i], false) - 1) / 2;
		for (int lz = -l; lz <= l; ++lz)
		{
			out.push_back(in[i] + to_string(lz));
		}
	}
	return out;
}

void initAtom(AtomType &a, const vector<string> &shellNames, const vector<vector<string> > &orbitNames)
{
	initAtom(a, shellNames, orbitNames, false);
}

vector<AtomType> OrbsToAtoms(const vector<string> &orbs, const bool hasSpin)
{
	uint lastNum = 0;
	vector < string > atomNames;
	vector < vector<string> > shellNames;
	vector < vector<vector<string> > > orbitNames;
	for (uint i = 0; i < (uint) orbs.size(); ++i)
	{
		uint atomNum, n;
		string atomName, shell;
		int lz;
		deconstructOrbname(orbs[i], atomNum, atomName, n, shell, lz);
		if ((uint) atomNames.size() == 0 || atomNum != lastNum || atomNames[atomNames.size() - 1].compare(atomName) != 0)
		{
			lastNum = atomNum;
			atomNames.push_back(atomName);
			shellNames.push_back(vector<string>());
			orbitNames.push_back(vector<vector<string> >());
		}
		string shellName;
		if (shell.compare("") == 0)
		{
			shellName = "";
		}
		else
		{
			shellName = to_string(n) + shell;
		}
		if ((uint) shellNames[(uint) shellNames.size() - 1].size() == 0 || shellNames[(uint) shellNames.size() - 1][(uint) shellNames[(uint) shellNames.size() - 1].size() - 1].compare(shellName) != 0)
		{
			shellNames[(uint) shellNames.size() - 1].push_back(shellName);
			orbitNames[(uint) orbitNames.size() - 1].push_back(vector<string>());
		}
		orbitNames[(uint) orbitNames.size() - 1][(uint) orbitNames[(uint) orbitNames.size() - 1].size() - 1].push_back(to_string(lz));
	}
	vector < AtomType > out = vector < AtomType > (atomNames.size());
	for (uint i = 0; i < (uint) atomNames.size(); ++i)
	{
		out[i] = AtomType();
		strcpy(out[i].Name, atomNames[i].c_str());
		initAtom(out[i], shellNames[i], orbitNames[i], hasSpin);
	}
	return out;
}

vector<AtomType> OrbsToAtoms(const vector<string> &orbs)
{
	return OrbsToAtoms(orbs, false);
}

