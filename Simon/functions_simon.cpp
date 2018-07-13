/* Some makeshift functions */

void debug(char* debugMessage)
{
	printf("\n\nDEBUG: ");
	printf("%s", debugMessage);
	printf("\n\n");
	fflush(stdout);
}



void PartialMeanFieldEigensystem()
{
	//TODO
}

void EigensystemSimon()
{
	//TODO: ruft partialMeanFieldEigensystem auf, wenn die Restrictions dies diktieren.
}

static int LuaPartialMeanFieldEigensystem(lua_State *L)
{
	//TODO
	return 0;
}

//initialises and modifies a list of indices that can, without any approximation, be treated on a mean field level.
void mandatoryMeanFieldList(uint16_t **indices, int * const indSetSize, const DeterminantRestrictionType det, int NFermion, int NBoson) //TODO: rewrite this in C++, using vectors
{
	uint32_t NWord = NBoson + Min(NFermion, 1 + (NFermion - 1) / 8);
	*indices = (uint16_t*) malloc(det.N * NWord * sizeof(uint16_t));
	int i;
	char * detStr;
	for (i = 0; i < det.N; ++i)
	{
		if (det.Min[i] == det.Max[i])
		{
			detStr = DetToString(NFermion, NBoson, det.Determinant[i]);
			debug(detStr);
//		TODO: add every index with 1 to a temporary vector.
//		TODO: compare number of entries to number of max and min electrons. If equal, add those indices to the return vector that are not already in it.
		}
	}
}

static int LuaOrbToMultiplicity(lua_State * L)
{
	char * orb;
	int multiplicity;
	uint32_t len = (uint32_t) lua_gettop(L);
	if (len != 1)
	{
		return luaL_error(L, "OrbToMultiplicity called with %d arguments. Only name of the orbital expected.\n", len);
	}
	orb = (char*) luaL_checkstring(L, 1);
	lua_pop(L, 1);
	multiplicity = nonRelAngularMomentumZNumbers(orb);
	if (multiplicity < 0)
		multiplicity = relAngularMomentumZNumbers(orb);
	if (multiplicity < 0)
	{
		printf("\nCould not interpret orbital \"");
		printf("%s", orb);
		printf("\" as either relativistic or non-relativistic atomic orbital, returning %d.", multiplicity);
	}
	lua_pushnumber(L, multiplicity);
	return 1;
}

static int LuaSimontest(lua_State * L)
{
	DeterminantRestrictionType testDet = *(DeterminantRestrictionType*) luaL_checkarray_DeterminantRestrictionType(L, 1);
	uint16_t * indList;
	int listLength;

//	char * detStr = DetToString(testDet.Determinant[1]);

	mandatoryMeanFieldList(&indList, &listLength, testDet, 3, 0);
	return 0;
}

// accepts a list of Wave Functions as well as the Hamiltonian and a list of other operators. The Hamiltonian is special here because it can be used to calculate the Deviance. It is optional.
static int LuaPrintExpectationValues(lua_State * L)
{
	WaveFunctionType *psiList;
	OperatorType *hamil, *opList;
	uint32_t numberOfWavefunctions, numberOfOperators;
	uint32_t len = (uint32_t) lua_gettop(L);
	if (len < 2 || len > 3)
		return luaL_error(L, "PrintExpectationValues called with %d arguments. 2 or 3 expected, List of Wavefunctions, optional Hamiltonian, and list of other operators.\n", len);
	if (lua_istable(L, 1))
	{
		psiList = luaL_checkarray_WaveFunctionType(L, 1, &numberOfWavefunctions);
	} else
	{
		numberOfWavefunctions = 1;
		psiList = (WaveFunctionType *) luaL_checkudata(L, 1, "Wavefunction_Type");
	}
	if (len == 2)
		hamil = (OperatorType*) luaL_checkudata(L, 2, "Operator_Type");
	if (lua_istable(L, -1))
	{
		opList = luaL_checkarray_OperatorType(L, -1, &numberOfOperators);
	} else
	{
		numberOfOperators = 1;
		opList = (OperatorType *) luaL_checkudata(L, -1, "Operator_Type");
	}

//TODO: weiter, falls BraKet nicht mehr komplett auf Lua Ebene stattfindet

	return 0;
}

//Feed it a table of tables of tables of tables of numbers, and it will return a table of tables.
//TODO: Optional: one Matrix and one or two numbers as inputs, then I embed M in an n times m zero Matrix?
//static int LuaMatrixFlatten(lua_State * L)
//{
//	int len = lua_gettop(L);
//	if (len != 1)
//		return luaL_error(L,
//				"MatrixFlatten called with %d arguments. 1 expected, Table of tables of tables of tables of numbers.\n",
//				len);
//	luaL_checktype(L, 1, LUA_TTABLE);
//	int numOfRowsOuter, numOfColsOuter, numOfAllRows = 0, numOfAllCols = 0;
//	numOfRowsOuter = (int) lua_objlen(L, 1);
//	lua_rawgeti(L, 1, 1);
//	luaL_checktype(L, -1, LUA_TTABLE);
//	numOfColsOuter = (int) lua_objlen(L, -1);
//	lua_pop(L, 1);
//	MatrixType allMs[numOfRowsOuter][numOfColsOuter];
//	uint32_t QComplex = 0;
//	int i, j, k, l;
//	for (i = 0; i < numOfRowsOuter; ++i)
//	{
//		lua_rawgeti(L, 1, i + 1);
//		luaL_checktype(L, -1, LUA_TTABLE);
//		for (j = 0; j < numOfColsOuter; ++j)
//		{
//			lua_rawgeti(L, -1, j + 1);
//			allMs[i][j] = *luaL_checkMatrix(L, -1);
//			lua_pop(L, 1);
//			if (i == 0)
//				numOfAllCols += allMs[i][j].Nj;
//			if (allMs[i][j].Nj != allMs[0][j].Nj)
//				return luaL_error(L,
//						"Matrix at position %d,%d has a different number of columns than matrix at position 1,%d!",
//						i + 1, j + 1, j + 1);
//			if (allMs[i][j].Ni != allMs[i][0].Ni)
//				return luaL_error(L,
//						"Matrix at position %d,%d has a different number of rows than matrix at position %d,%1!",
//						i + 1, j + 1, i + 1);
//			if (!QComplex && allMs[i][j].QComplex)
//				QComplex = allMs[i][j].QComplex;
//		}
//		numOfAllRows += allMs[i][0].Ni;
//		lua_pop(L, 1);
//	}
//	MatrixType fullM;
//	fullM.QComplex = QComplex;
//	fullM.Ni = numOfAllRows;
//	fullM.Nj = numOfAllCols;
//	InitMatrix(&fullM);
//	int rowOffset = 0, colOffset = 0;
//	for (i = 0; i < numOfRowsOuter; ++i)
//	{
//		colOffset = 0;
//		for (j = 0; j < numOfColsOuter; ++j)
//		{
//			for (k = 0; k < allMs[i][j].Ni; ++k)
//			{
//				for (l = 0; l < allMs[i][j].Nj; ++l)
//				{
//					fullM.VRe[rowOffset + k][colOffset + l] = allMs[i][j].VRe[k][l];
//					if (QComplex)
//					{
//						if (allMs[i][j].QComplex)
//							fullM.VIm[rowOffset + k][colOffset + l] = allMs[i][j].VIm[k][l];
//						else
//							fullM.VIm[rowOffset + k][colOffset + l] = 0;
//					}
//				}
//			}
//			colOffset += allMs[i][j].Nj;
//		}
//		rowOffset += allMs[i][0].Ni;
//	}
//	LuaPushTable_MatrixType(L, &fullM);
//	return 1;
//}

MatrixType zeroMatrix(int Ni, int Nj)
{
	MatrixType zeroMat; //TODO: Gibt mir das ein Problem, weil immer die gleiche Matrix initiiert wird? Im Zweifel als extraObjekt auslagern.
	zeroMat.Ni = Ni;
	zeroMat.Nj = Nj;
	if (Ni != 0 && Nj != 0)
		InitMatrix(&zeroMat);
	return zeroMat;
}

MatrixType nullMatrix()
{
	return zeroMatrix(0,0);
}

bool isNullMatrix(MatrixType M)
{
	return (M.Ni == 0 && M.Nj == 0);
}

MatrixType MatrixFlatten(std::vector<std::vector<MatrixType> > allMs)
{
	uint32_t QComplex = 0;
	MatrixType fullM;
	int numOfRowsOuter = allMs.size();
	int numOfColsOuter = allMs[0].size();
	int numOfRowsInner[numOfRowsOuter], numOfColsInner[numOfColsOuter], numOfAllRows = 0, numOfAllCols = 0;
	for (int i = 0; i < numOfRowsOuter; ++i)
	{
		numOfRowsInner[i] = 0;
		for (int j = 0; j < numOfColsOuter; ++j)
		{
			if (!isNullMatrix(allMs[i][j]))
			{
				numOfRowsInner[i] = allMs[i][j].Ni;
				break;
			}
		}
		if (numOfRowsInner[i] == 0)
			printf("\nWarning: Row %d of matrix is empty and will be ignored.\n",i+1);
		numOfAllRows += numOfRowsInner[i];
	}
	for (int j = 0; j < numOfColsOuter; ++j)
	{
		numOfColsInner[j] = 0;
		for (int i = 0; i < numOfRowsOuter; ++i)
		{
			if (!isNullMatrix(allMs[i][j]))
			{
				numOfColsInner[j] = allMs[i][j].Nj;
				break;
			}
		}
		if (numOfColsInner[j] == 0)
			printf("\nWarning: Column %d of matrix is empty and will be ignored.\n",j+1);
		numOfAllCols += numOfColsInner[j];
	}
	for (int i = 0; i < numOfRowsOuter; ++i)
	{
		for (int j = 0; j < numOfColsOuter; ++j)
		{
			if (isNullMatrix(allMs[i][j]))
			{
				allMs[i][j] = zeroMatrix(numOfRowsInner[i],numOfColsInner[j]);
			}
			if (allMs[i][j].Nj != numOfColsInner[j])
			{
				printf("Matrix at position %d,%d has a different number of columns than it should have.", i + 1, j + 1);
				return nullMatrix();
			}
			if (allMs[i][j].Ni != numOfRowsInner[i])
			{
				printf("Matrix at position %d,%d has a different number of rows than it should have.", i + 1, j + 1);
				return nullMatrix();
			}
			if (!QComplex && allMs[i][j].QComplex)
				QComplex = allMs[i][j].QComplex;
		}
	}
	fullM.QComplex = QComplex;
	fullM.Ni = numOfAllRows;
	fullM.Nj = numOfAllCols;
	InitMatrix(&fullM);
	int rowOffset = 0, colOffset = 0;
	for (int i = 0; i < numOfRowsOuter; ++i)
	{
		colOffset = 0;
		for (int j = 0; j < numOfColsOuter; ++j)
		{
			for (int k = 0; k < allMs[i][j].Ni; ++k)
			{
				for (int l = 0; l < allMs[i][j].Nj; ++l)
				{
					fullM.VRe[rowOffset + k][colOffset + l] = allMs[i][j].VRe[k][l];
					if (QComplex)
					{
						if (allMs[i][j].QComplex)
							fullM.VIm[rowOffset + k][colOffset + l] = allMs[i][j].VIm[k][l];
						else
							fullM.VIm[rowOffset + k][colOffset + l] = 0;
					}
				}
			}
			colOffset += allMs[i][j].Nj;
		}
		rowOffset += allMs[i][0].Ni;
	}
	return fullM;
}

//Feed it a table of tables of tables of tables of numbers, and it will return a table of tables.
//TODO: Optional: one Matrix and one or two numbers as inputs, then I embed M in an n times m zero Matrix?
static int LuaMatrixFlatten(lua_State * L)
{
	int len = lua_gettop(L);
	if (len != 1)
		return luaL_error(L, "MatrixFlatten called with %d arguments. 1 expected, Table of tables of tables of tables of numbers.\n", len);
	luaL_checktype(L, 1, LUA_TTABLE);
	int numOfRowsOuter, numOfColsOuter;
	numOfRowsOuter = (int) lua_objlen(L, 1);
	lua_rawgeti(L, 1, 1);
	luaL_checktype(L, -1, LUA_TTABLE);
	numOfColsOuter = (int) lua_objlen(L, -1);
	lua_pop(L, 1);
	std::vector<std::vector<MatrixType> > allMs(numOfRowsOuter);
//	MatrixType allMs[numOfRowsOuter][numOfColsOuter];
	for (int i = 0; i < numOfRowsOuter; ++i)
	{
		lua_rawgeti(L, 1, i + 1);
		luaL_checktype(L, -1, LUA_TTABLE);
		std::vector<MatrixType> col (numOfColsOuter);
		allMs[i] = col;
		for (int j = 0; j < numOfColsOuter; ++j)
		{
			lua_rawgeti(L, -1, j + 1);
			if (lua_isnumber(L, -1))
			{
				if (luaL_checknumber(L, -1) == 0)
					allMs[i][j] = nullMatrix();
				else
					return luaL_error(L, "MatrixFlatten called with entry %d. Matrix or 0 expected.\n", luaL_checknumber(L, -1));
			} else
			{
				allMs[i][j]  = *luaL_checkMatrix(L, -1);
			}
			lua_pop(L, 1);
		}
		lua_pop(L, 1);
	}
	MatrixType fullM = MatrixFlatten(allMs);
	LuaPushTable_MatrixType(L, &fullM);
	return 1;
}


static const struct luaL_Reg Simon_meta[] = {
  {NULL, NULL}
};

static const struct luaL_Reg Simon_func[] = {
  {"test", LuaSimontest},
	{"MatrixFlatten",LuaMatrixFlatten},
	{"OrbToMultiplicity",LuaOrbToMultiplicity},
	{"PartialMeanFieldEigensystem",LuaPartialMeanFieldEigensystem},
	{"PrintExpectationValues",LuaPrintExpectationValues},
  {NULL, NULL}
};
