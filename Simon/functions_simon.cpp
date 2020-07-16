
template<typename T>
void debug(const T debugMessage)
{
	cout << endl << endl << "DEBUG: " << debugMessage << endl << endl;
	fflush(stdout);
}

void debug()
{
	debug("");
}

void debug(const Complex16 c)
{
	debug();
	cout << c.Re << " + " << c.Im << endl;
}

template<typename T>
void debug(const vector<T> &vec)
{
	debug();
	for (int i = 0; i < vec.size(); ++i)
		cout << "elem=" << i << ", vec[i]=" << vec[i] << endl;
}

void debug(const vector<Complex16> &vec)
{
	debug();
	for (int i = 0; i < vec.size(); ++i)
		cout << "elem=" << i << ", vec[i]=" << vec[i].Re << " + " << vec[i].Im << endl;
}

template<typename keyType, typename valType>
void debug(const map<keyType, valType> &m)
{
	debug();
	typename map<keyType, valType>::const_iterator iter;
	for (iter = m.begin(); iter != m.end(); ++iter)
		cout << "Key=" << iter->first << ", Val=" << iter->second << endl;
}

void debug(const CompactMatrixType &M, const string &str)
{
	debug(str);
	for (int i = 0; i < M.Ni; ++i)
	{
		for (int j = 0; j < M.Nj; ++j)
		{
			if ((bool) M.QComplex)
				if ((abs(M.V[i * M.Nj + j].Re) < 1E-5) && (abs(M.V[i * M.Nj + j].Im) < 1E-5))
					cout << "0\t";
				else
					cout << numberToStringMaxWidth(M.V[i * M.Nj + j].Re, 4) << "+" << numberToStringMaxWidth(M.V[i * M.Nj + j].Im, 4) << "i" << "\t";
			else
				cout << numberToStringMaxWidth(M.VRe[i * M.Nj + j], 6) << "\t";
		}
		cout << endl;
	}
	cout << endl << endl;
}

void debug(const CompactMatrixType &M, vector<string> vec, const string &str)
{
	debug(str);
	cout << "\t";
	for (int i = 0; i < (int) vec.size(); ++i)
	{
		if ((uint) vec[i].size() > 6) vec[i] = vec[i].substr(vec[i].size() - 6);
		cout << vec[i] << "\t";
	}
	cout << endl;
	for (int i = 0; i < M.Ni; ++i)
	{
		cout << vec[i] << "\t";
		for (int j = 0; j < M.Nj; ++j)
		{
			if ((bool) M.QComplex)
				if ((abs(M.V[i * M.Nj + j].Re) < 1E-5) && (abs(M.V[i * M.Nj + j].Im) < 1E-5))
					cout << "0\t";
				else
					cout << numberToStringMaxWidth(M.V[i * M.Nj + j].Re, 4) << "+" << numberToStringMaxWidth(M.V[i * M.Nj + j].Im, 4) << "i" << "\t";
			else
				cout << numberToStringMaxWidth(M.VRe[i * M.Nj + j], 6) << "\t";
		}
		cout << endl;
	}
	cout << endl << endl;
}

void debug(const MatrixType &M, const string &str)
{
	debug(str);
	cout << "{" << endl;
	for (int i = 0; i < M.Ni; ++i)
	{
		cout << "{\t";
		for (int j = 0; j < M.Nj; ++j)
		{
			cout << numberToStringMaxWidth(M.VRe[i][j], 6) << "\t";
		}
		cout << "}" << endl;
	}
	cout << "}" << endl << endl;
}

void debug(const CompactMatrixType &M)
{
	debug(M, "Matrix-Print");
}

void debug(const CompactMatrixType &M, vector<string> vec)
{
	debug(M, vec, "Matrix- and Vector-Print");
}

void debug(const MatrixType &M)
{
	debug(M, "Matrix-Print");
}

void debug(const doubleVecValPair &vp)
{
	debug(vp.vec);
	debug(vp.val);
}

static int LuaSimontest(lua_State * L)
{
	return 0;
}

static const struct luaL_Reg Simon_meta[] = {
  {NULL, NULL}
};

static const struct luaL_Reg Simon_func[] = {
	{"AtomicNumberOperators", LuaAtomicNumberOperators},
	{"AtomicToWannierMatrix", LuaAtomicToWannierMatrix},
	{"CorthogonalShells", LuaCorthogonalShells},
	{"ExtendedLigandField", LuaExtendedLigandField},
	{"Flatten", LuaFlatten},
	{"IsosurfaceFixedInt", LuaIsosurfaceFixedInt},
	{"IsosurfaceFixedVal", LuaIsosurfaceFixedVal},
	{"MatrixToOperator", LuaMatrixToOperator},
	{"MatrixToSparseMatrix", LuaMatrixToSparseMatrix},
	{"MatrixTrace", LuaMatrixTrace},
	{"OrthogonalizeWithNonOrthoBasis", LuaOrthogonalizeWithNonOrthoBasis},
	{"lsCouplingConstants", LuaLSCouplingConstants},
	{"FileReadOrca", LuaFileReadOrca},
	{"RealspaceDensity", LuaRealspaceDensity},
	{"RemoveSpin", LuaRemoveSpin},
	{"SortByKappa",LuaSortByKappa},
	{"SparseMatrixToMatrix",LuaSparseMatrixToMatrix},
	{"TensorProduct", LuaTensorProduct},
	{"test", LuaSimontest},
	{NULL, NULL}
};
