 
OperatorType CompactMatrixToOperator(const CompactMatrixType &M)
{
	OperatorType op;
	InitOperator(&op);
	InitOperatorFixedLength(&op, 2, 8);
	op.QComplex = M.QComplex;
	OperatorChangeNFandNB(&op, M.Ni, 0);
	uint16_t oadd[2];
	if ((bool) op.QComplex)
	{
		for (uint16_t i = 0; i < (uint16_t) M.Ni; ++i)
		{
			for (uint16_t j = 0; j < (uint16_t) M.Nj; ++j)
			{
				double VRe = M.V[i * M.Nj + j].Re;
				double VIm = M.V[i * M.Nj + j].Im;
				oadd[1] = i + CrAnMask;
				oadd[0] = j;
				OperatorAddLader(&op, oadd, 2, VRe, VIm);
			}
		}
	} else
	{
		for (uint16_t i = 0; i < (uint16_t) M.Ni; ++i)
		{
			for (uint16_t j = 0; j < (uint16_t) M.Nj; ++j)
			{
				double VRe = M.VRe[i * M.Nj + j];
				oadd[1] = i + CrAnMask;
				oadd[0] = j;
				OperatorAddLader(&op, oadd, 2, VRe, 0.);
			}
		}
	}
	return op;
}

static int LuaMatrixToOperator(lua_State * L)
{
	CompactMatrixType M;
	uint len = (uint) lua_gettop(L);
	if (len != 1)
		return luaL_error(L, "MatrixToOperator called with %d arguments, only 1 expected.", len);
	if (!luaL_isType(L, 1, M))
		return luaL_error(L, "MatrixToOperator called with an argument that is not a matrix.");
	lua_check(L, 1, M);
	lua_push(L, CompactMatrixToOperator(M));
	return 1;
}
