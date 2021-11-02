void lua_push(lua_State *L, const SparseMatrixType *const pushkin)
{
	SparseMatrixType *B = (SparseMatrixType*) lua_newuserdata(L, sizeof(SparseMatrixType));
	luaL_setmetatable(L, "SparseMatrix_Type");
	*B = *pushkin;
}

void lua_push(lua_State *L, const SparseMatrixType &pushkin)
{
	lua_push(L, &pushkin);
}

void lua_check(lua_State *L, int arg, SparseMatrixType &in)
{
	in = *(SparseMatrixType *) luaL_checkudata(L, arg, "SparseMatrix_Type");
}

void lua_check(lua_State *L, int arg, SparseMatrixType* &in)
{
	in = (SparseMatrixType *) luaL_checkudata(L, arg, "SparseMatrix_Type");
}

bool luaL_isType(lua_State *L, int arg, const SparseMatrixType &example)
{
	return luaL_isudata(L, arg, "SparseMatrix_Type");
}

SparseMatrixType CompactMatrixToSparseMatrix(const CompactMatrixType &M)
{
	SparseMatrixType S;
	S.QComplex = M.QComplex;
	S.Ni = M.Ni;
	S.Nj = M.Nj;
	//Assume Number of nonzero Elements ~ Ni*Nj/2 => Advised value 2^NBitsKey ~ Ni*Nj/2
	S.NBitsKey = (uint32_t) log((double) (M.Ni * M.Nj)) / log(2.);
	if (S.NBitsKey == 0) S.NBitsKey = 1;
	S.a0 = M.a0;
	S.b0 = M.b0;
	S.mu = M.mu;
	InitSparseMatrix(&S);

	double epsilon = 10 * DBL_EPSILON;
	if ((bool) S.QComplex)
	{
		for (uint i = 0; i < M.Ni; ++i)
		{
			for (uint j = 0; j < M.Nj; ++j)
			{
				Complex16 elem = M.V[i * M.Nj + j];
				if (fabs(elem.Re) > epsilon || fabs(elem.Im) > epsilon)
				{
					ComplexSparseMatrixAddElement(&S, i, j, elem.Re, elem.Im);
				}
			}
		}
	}
	else
	{
		for (uint i = 0; i < M.Ni; ++i)
		{
			for (uint j = 0; j < M.Nj; ++j)
			{
				double elem = M.VRe[i * M.Nj + j];
				if (fabs(elem) > epsilon)
				{
					RealSparseMatrixAddElement(&S, i, j, elem);
				}
			}
		}
	}
	return S;
}

CompactMatrixType SparseMatrixToCompactMatrix(const SparseMatrixType &S)
{
	CompactMatrixType M;
	M.QComplex = S.QComplex;
	M.Ni = S.Ni;
	M.Nj = S.Nj;
	InitCompactMatrix(&M);
	M.a0 = S.a0;
	M.b0 = S.b0;
	M.mu = S.mu;
	if ((bool) M.QComplex)
	{
		for (uint i = 0; i < M.Ni; ++i)
		{
			for (uint j = 0; j < M.Nj; ++j)
			{
				M.V[i * M.Nj + j].Re = 0;
				M.V[i * M.Nj + j].Im = 0;
			}
		}
		for (uint k = 0; k < S.N; ++k)
		{
			uint ind = S.i[k] * M.Nj + S.j[k];
			M.V[ind].Re = S.VRe[k];
			M.V[ind].Im = S.VIm[k];
		}
	}
	else
	{
		for (uint i = 0; i < M.Ni; ++i)
		{
			for (uint j = 0; j < M.Nj; ++j)
			{
				M.VRe[i * M.Nj + j] = 0;
			}
		}
		for (uint k = 0; k < S.N; ++k)
		{
			uint ind = S.i[k] * M.Nj + S.j[k];
			M.VRe[ind] = S.VRe[k];
		}
	}
	return M;
}


static int LuaMatrixToSparseMatrix(lua_State * L)
{
	CompactMatrixType M;
	uint len = (uint) lua_gettop(L);
	if (len != 1)
		return luaL_error(L, "MatrixToSparseMatrix called with %d arguments, only 1 expected.", len);
	if (!luaL_isType(L, 1, M))
		return luaL_error(L, "MatrixToSparseMatrix called with an argument that is not a matrix.");
	lua_check(L, 1, M);
	lua_push(L, CompactMatrixToSparseMatrix(M));
	FreeCompactMatrix(&M);
	return 1;
}

static int LuaSparseMatrixToMatrix(lua_State * L)
{
	SparseMatrixType M;
	uint len = (uint) lua_gettop(L);
	if (len != 1)
		return luaL_error(L, "SparseMatrixToMatrix called with %d arguments, only 1 expected.", len);
	if (!luaL_isType(L, 1, M))
		return luaL_error(L, "SparseMatrixToMatrix called with an argument that is not a sparse matrix.");
	lua_check(L, 1, M);
	lua_push(L, SparseMatrixToCompactMatrix(M));
	return 1;
}

double RealSparseMatrixTrace(const SparseMatrixType *M)
{
	double trace = 0;
	for (uint k = 0; k < M->N; ++k)
		if (M->i[k] == M->j[k]) trace += M->VRe[k];
	return trace;
}

Complex16 ComplexSparseMatrixTrace(const SparseMatrixType *M)
{
	Complex16 trace;
	trace.Re = 0;
	trace.Im = 0;
	for (uint k = 0; k < M->N; ++k)
	{
		if (M->i[k] == M->j[k])
		{
			trace.Re += M->VRe[k];
			trace.Im += M->VIm[k];
		}
	}
	return trace;
}

static int LuaSparseMatrixUnm(lua_State *L)
{
	SparseMatrixType *Sin, *Sout;
	lua_check(L, 1, Sin);
	Sout = (SparseMatrixType*) lua_newuserdata(L, sizeof(SparseMatrixType));
	luaL_setmetatable(L, "SparseMatrix_Type");
	SparseMatrixCopy(Sout, *Sin);
	if ((bool) Sout->QComplex)
	{
		for (uint k = 0; k < Sout->N; ++k)
		{
			Sout->VRe[k] = -Sout->VRe[k];
			Sout->VIm[k] = -Sout->VIm[k];
		}
	}
	else
	{
		for (uint k = 0; k < Sout->N; ++k)
		{
			Sout->VRe[k] = -Sout->VRe[k];
		}
	}
	return 1;
}

static int LuaSparseMatrixAdd(lua_State *L)
{
	SparseMatrixType *A, *B, *C;
	lua_check(L, 1, A);
	lua_check(L, 2, B);
	if ((A->Ni != B->Ni) || (A->Nj != B->Nj))
	{
		return luaL_error(L, "Error: dimensions of matrix A (%d %d) and B (%d %d) in A+B not equivalent\n", A->Ni, A->Nj, B->Ni, B->Nj);
	}
	C = (SparseMatrixType*) lua_newuserdata(L, sizeof(SparseMatrixType));
	luaL_setmetatable(L, "SparseMatrix_Type");
	SparseMatrixCopy(C, *A);
	SparseMatrixAdd(C, *B, 1, 0);
	//TODO: This gives a segmentation fault:
	//SparseMatrixChop(C, 10 * DBL_EPSILON);
	return 1;
}

static int LuaSparseMatrixSub(lua_State *L)
{
	SparseMatrixType *A, *B, *C;
	lua_check(L, 1, A);
	lua_check(L, 2, B);
	if ((A->Ni != B->Ni) || (A->Nj != B->Nj))
	{
		return luaL_error(L, "Error: dimensions of matrix A (%d %d) and B (%d %d) in A+B not equivalent\n", A->Ni, A->Nj, B->Ni, B->Nj);
	}
	C = (SparseMatrixType*) lua_newuserdata(L, sizeof(SparseMatrixType));
	luaL_setmetatable(L, "SparseMatrix_Type");
	SparseMatrixCopy(C, *A);
	SparseMatrixAdd(C, *B, -1, 0);
	//TODO: This gives a segmentation fault:
	//SparseMatrixChop(C, 10 * DBL_EPSILON);
	return 1;
}

static int LuaSparseMatrixMul(lua_State *L)
{
	SparseMatrixType S1, S2, *Sout;
	double d;
	Complex16 c;
	bool isTimesMat = false, isTimesDouble = false, isTimesComplex = false;
	if (luaL_isType(L, 1, S1) && luaL_isType(L, 2, S2))
	{
		lua_check(L, 1, S1);
		lua_check(L, 2, S2);
		isTimesMat = true;
	}
	else if (luaL_isType(L, 1, d) && luaL_isType(L, 2, S1))
	{
		lua_check(L, 1, d);
		lua_check(L, 2, S1);
		isTimesDouble = true;
	}
	else if (luaL_isType(L, 1, S1) && luaL_isType(L, 2, d))
	{
		lua_check(L, 1, S1);
		lua_check(L, 2, d);
		isTimesDouble = true;
	}
	else if (luaL_isType(L, 1, c) && luaL_isType(L, 2, S1))
	{
		lua_check(L, 1, c);
		lua_check(L, 2, S1);
		isTimesComplex = true;
	}
	else if (luaL_isType(L, 1, S1) && luaL_isType(L, 2, c))
	{
		lua_check(L, 1, S1);
		lua_check(L, 2, c);
		isTimesComplex = true;
	}
	else
	{
		return luaL_error(L, "Error: Sparse Matrix Multiplication called with an unknown type!");
	}


	Sout = (SparseMatrixType*) lua_newuserdata(L, sizeof(SparseMatrixType));
	luaL_setmetatable(L, "SparseMatrix_Type");
	if (isTimesMat)
	{
		SparseMatrixMultiply(Sout, S1, S2);
	}
	else if (isTimesDouble)
	{
		SparseMatrixCopy(Sout, S1);
		SparseMatrixScalarMultiply(Sout, d, 0);
	}
	else if (isTimesComplex)
	{
		SparseMatrixCopy(Sout, S1);
		SparseMatrixScalarMultiply(Sout, c.Re, c.Im);
	}
	return 1;
}

static int LuaSparseMatrixDiv(lua_State *L)
{
	SparseMatrixType *Sin, *Sout;
	double d;
	Complex16 c, cinv;

	lua_check(L, 1, Sin);
	Sout = (SparseMatrixType*) lua_newuserdata(L, sizeof(SparseMatrixType));
	luaL_setmetatable(L, "SparseMatrix_Type");
	SparseMatrixCopy(Sout, *Sin);
	if (luaL_isType(L, 2, d))
	{
		lua_check(L, 2, d);
		SparseMatrixScalarMultiply(Sout, 1 / d, 0);
	}
	else if (luaL_isType(L, 2, c))
	{
		lua_check(L, 2, c);
		cinv.Re = c.Re / (Sqr(c.Re) + Sqr(c.Im));
		cinv.Im = c.Im / (Sqr(c.Re) + Sqr(c.Im));
		SparseMatrixScalarMultiply(Sout, cinv.Re, cinv.Im);
	}
	else
	{
		return luaL_error(L, "Error: Sparse Matrix Division called with an unknown type!");
	}
	return 1;
}

static int LuaSparseMatrixPow(lua_State *L)
{
	cout << "TODO" << endl;
	return 0;
}

static int LuaSparseMatrixEq(lua_State *L)
{
	cout << "TODO" << endl;
	return 0;
}

static int LuaSparseMatrixLen(lua_State *L)
{
	cout << "TODO" << endl;
	return 0;
}

static int LuaSparseMatrixToString(lua_State *L)
{
	SparseMatrixType S;
	lua_check(L, 1, S);
	std::stringstream sstream;
	sstream << "Sparse Matrix print:" << endl;
	for (uint k = 0; k < S.N; ++k)
	{
		sstream << "Element " << S.i[k] + 1 << "," << S.j[k] + 1 << " = ";
		sstream << S.VRe[k];
		if ((bool) S.QComplex)
		{
			sstream << " ";
			if (S.VIm[k] >= 0)
			{
				sstream << "+";
			}
			sstream << S.VIm[k] << " i";
		}
		sstream << endl;
	}
	lua_pushstring(L, sstream.str());
	return 1;
}

static int LuaSparseMatrixIndex(lua_State *L)
{
	cout << "TODO" << endl;
	return 0;
}

static int LuaSparseMatrixDestroy(lua_State *L)
{
	SparseMatrixType *A;
	lua_check(L, 1, A);
	FreeSparseMatrix(A);
	return 0;
}

static const struct luaL_Reg SparseMatrix_meta[] = {
  {"__unm", LuaSparseMatrixUnm},
  {"__add", LuaSparseMatrixAdd},
  {"__sub", LuaSparseMatrixSub},
  {"__mul", LuaSparseMatrixMul},
  {"__div", LuaSparseMatrixDiv},
  {"__pow", LuaSparseMatrixPow},
	{"__eq", LuaSparseMatrixEq},
	{"__len", LuaSparseMatrixLen},
  {"__tostring", LuaSparseMatrixToString},
  {"__index", LuaSparseMatrixIndex},
  {"__gc", LuaSparseMatrixDestroy},
  {NULL, NULL}
};
