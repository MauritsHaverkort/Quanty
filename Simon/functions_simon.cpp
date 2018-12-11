
template<typename T>
void debug(T debugMessage)
{
	cout << endl << endl << "DEBUG: " << debugMessage << endl << endl;
	fflush(stdout);
}

void debug()
{
	debug("");
}

template<typename T>
void debug(vector<T> vec)
{
	debug();
	for (int i = 0; i < vec.size(); ++i)
		cout << "elem=" << i << ", vec[i]=" << vec[i] << endl;
}

template<typename keyType, typename valType>
void debug(const map<keyType, valType> &m)
{
	debug();
	typename map<keyType, valType>::const_iterator iter;
	for (iter = m.begin(); iter != m.end(); ++iter)
		cout << "Key=" << iter->first << ", Val=" << iter->second << endl;
}

void debug(CompactMatrixType M, string str)
{
	debug(str);
	cout << "{" << endl;
	for (int i = 0; i < M.Ni; ++i)
	{
		cout << "{\t";
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
		cout << "}" << endl;
	}
	cout << "}" << endl << endl;
}

void debug(MatrixType M, string str)
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

void debug(CompactMatrixType M)
{
	debug(M, "Matrix-Print");
}

void debug(MatrixType M)
{
	debug(M, "Matrix-Print");
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
	{"MatrixToOperator", LuaMatrixToOperator},
	{"ReduceHilbertSpace", LuaReduceHilbertSpace},
  {"test", LuaSimontest},
  {NULL, NULL}
};
