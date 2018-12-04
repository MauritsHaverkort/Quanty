string cutLeadingSpaces(string in)
{
	while (in.find(' ') == 0 || in.find('\t') == 0)
		in = in.substr(1);
	return in;
}

string cutTrailingSpaces(string in)
{
	if (in.length() < 2)
		return in;
	while (in.rfind(' ') == in.length() - 1 || in.rfind('\t') == in.length() - 1)
		in = in.substr(0, in.length() - 1);
	return in;
}

string cutBorderSpaces(string in)
{
	string out = cutLeadingSpaces(in);
	return cutTrailingSpaces(out);
}

bool atomicStateSorter(string s1, string s2)
{
//	Expects input like 2_Fe_3p1 or Fe_3p1
	size_t pos1, pos2;
	pos1 = s1.find('_');
	uint siteNum1 = (uint) strtol(s1.substr(0, pos1).c_str(), NULL, 10);
	pos2 = s2.find('_');
	uint siteNum2 = (uint) strtol(s2.substr(0, pos1).c_str(), NULL, 10);
	if (siteNum1 != siteNum2 && siteNum1 != 0 && siteNum2 != 0)
		return siteNum1 < siteNum2;
	pos1 = s1.find_last_of('_');
	++pos1;
	pos2 = s2.find_last_of('_');
	++pos2;
	uint n1 = (uint) strtol(s1.substr(pos1, 1).c_str(), NULL, 10);
	uint n2 = (uint) strtol(s2.substr(pos2, 1).c_str(), NULL, 10);
	if (n1 != n2)
		return n1 < n2;
	++pos1;
	++pos2;
	uint mult1 = orbToMultiplicity(s1.substr(pos1, 1));
	uint mult2 = orbToMultiplicity(s2.substr(pos2, 1));
	if (mult1 != mult2)
		return mult1 < mult2;
	++pos1;
	++pos2;
	int lz1 = (int) strtol(s1.substr(pos1).c_str(), NULL, 10);
	int lz2 = (int) strtol(s2.substr(pos2).c_str(), NULL, 10);
	return lz1 < lz2;
}

vector<map<string, double> > getWF_coefficients_Filecontent(const string &filename, vector<string> &orbitalNames)
{
	typedef map<string, double> M;
	vector<M> out;
	ifstream file(filename.c_str());
	if (!file.good())
	{
		cout << "Error in getWF_coefficients_Filecontent: Could not open file with path " << filename << "." << endl;
		return out;
	}
	string line, site, orb, state;
	uint oldWanNumber = 1, wanNumber, siteNum;
	int lz;
	double weight;
	size_t pos1, pos2;
	M orbWeightsMap = M();
	while (getline(file, line))
	{
		line = cutBorderSpaces(line);
		if (line.find("WF(") != string::npos)
		{
			pos1 = line.find("WF(");
			pos1 += 3;
			pos2 = line.find(")", pos1);
			wanNumber = (uint) strtol(line.substr(pos1, pos2 - pos1).c_str(), NULL, 10);
			if (wanNumber != oldWanNumber)
			{
				out.push_back(orbWeightsMap);
				orbWeightsMap = M();
				oldWanNumber = wanNumber;
			}

			pos1 = line.find("site ", pos2);
			pos1 += 5;
			pos2 = line.find("(", pos1);
			site = cutBorderSpaces(line.substr(pos1, pos2 - pos1));
			pos1 = pos2 + 1;
			pos2 = line.find(")", pos1);
			siteNum = (uint) strtol(line.substr(pos1, pos2 - pos1).c_str(), NULL, 10);
		} else if (line.find("vec=") != string::npos)
		{
//			Do nothing
		} else if (line.length() <= 3 && line.length() >= 2)
		{
			orb = line;
		} else if (line.find(":") != string::npos)
		{
			pos1 = 0;
			pos2 = line.find(":");
			lz = (int) strtol(line.substr(pos1, pos2 - pos1).c_str(), NULL, 10);
			weight = (double) strtod(line.substr(pos2 + 1).c_str(), NULL);
			state = to_string(siteNum) + "_" + site + "_" + orb + to_string(lz);
			if (!isElementOf(state, orbitalNames))
				orbitalNames.push_back(state);
			orbWeightsMap.insert( { state, weight });
		}
	}
	out.push_back(orbWeightsMap);
	sort(orbitalNames.begin(), orbitalNames.end(), atomicStateSorter);
	return out;
}

CompactMatrixType AtomicToWannierMatrix(const string &filename, vector<string> &orbitalNames)
{
	vector < map<string, double> > filecontent = getWF_coefficients_Filecontent(filename, orbitalNames);
	CompactMatrixType out = CompactMatrixType();
	out.Ni = (uint) filecontent.size();
	out.Nj = (uint) orbitalNames.size();
	out.QComplex = 0;
	InitCompactMatrix(&out);
	for (uint i = 0; i < (uint) out.Ni; ++i)
		for (uint j = 0; j < (uint) out.Nj; ++j)
			if (filecontent[i].find(orbitalNames[j]) != filecontent[i].end())
				out.VRe[i * out.Nj + j] = filecontent[i][orbitalNames[j]];
			else
				out.VRe[i * out.Nj + j] = 0;
	return out;
}

bool isFileAcccessible(const string &name)
{
	ifstream f(name.c_str());
	return f.good();
}

static int LuaAtomicToWannierMatrix(lua_State * L)
{
	string filepath;
	vector < string > orbitalNames;
	uint32_t len = (uint32_t) lua_gettop(L);
	if (len != 1)
		return luaL_error(L, "Error: AtomicToWannierMatrix called with %d elements, only 1, a filepath, expected.", len);
	lua_check(L, 1, filepath);
	if (!isFileAcccessible(filepath))
		return luaL_error(L, "Error: Could not read file with path %s.", filepath.c_str());
	lua_push(L, AtomicToWannierMatrix(filepath, orbitalNames));
	lua_push(L, orbitalNames);
	return 2;
}
