/* mwh this function already excisted
  under the same name, using slightly different formatting
   please use the excisting function

CompactMatrixType BlockTriDiagonalMatrixToCompactMatrix(const BlockTriDiagonalMatrixType &M, const double zero)
{
	vector < vector<CompactMatrixType> > outMatMat = CompactnullMatrixMatrix(M.N - 1, M.N - 1);
	for (uint i = 1; i < M.N; ++i)
	{
		CompactMatrixType Atmp, Btmp, BtmpT;
		Atmp.QComplex = M.QComplex;
		Btmp.QComplex = M.QComplex;
		Atmp.Ni = M.BlockSize;
		Atmp.Nj = M.BlockSize;
		Btmp.Ni = M.BlockSize;
		Btmp.Nj = M.BlockSize;
		InitCompactMatrix(&Atmp);
		InitCompactMatrix(&Btmp);
		for (uint k = 0; k < M.BlockSize; ++k)
		{
			for (uint l = 0; l < M.BlockSize; ++l)
			{
				uint indOld = i * M.BlockSize * M.BlockSize + k * M.BlockSize + l;
				uint indNew = k * M.BlockSize + l;
				if ((bool) M.QComplex)
				{
					Atmp.V[indNew].Re = M.a[indOld].Re;
					Atmp.V[indNew].Im = M.a[indOld].Im;
					Btmp.V[indNew].Re = M.b[indOld].Re;
					Btmp.V[indNew].Im = M.b[indOld].Im;
				} else
				{
					Atmp.VRe[indNew] = M.aRe[indOld];
					Btmp.VRe[indNew] = M.bRe[indOld];
				}
			}
		}
		FreeCompactMatrix(&(outMatMat[i - 1][i - 1]));
		outMatMat[i - 1][i - 1] = Atmp;
		if (isZeroMatrix(Btmp, zero))
		{
			FreeCompactMatrix(&Btmp);
			break;
		}
		if (i < M.N - 1)
		{
			FreeCompactMatrix(&(outMatMat[i - 1][i]));
			outMatMat[i - 1][i] = Btmp;
			CompactMatrixCopy(Btmp, &BtmpT);
			if ((bool) Btmp.QComplex)
				CompactMatrixConjugateTranspose(&BtmpT);
			else
				CompactMatrixTranspose(&BtmpT);
			FreeCompactMatrix(&(outMatMat[i][i - 1]));
			outMatMat[i][i - 1] = BtmpT;
		}
	}
	CompactMatrixType out = MatrixFlatten(outMatMat);
	strcpy(out.Name, M.Name);
	out.mu = M.mu;
	FreeCompactMatrixMatrix(outMatMat);
	return out;
}

CompactMatrixType BlockTriDiagonalMatrixToCompactMatrix(const BlockTriDiagonalMatrixType &M)
{
	return BlockTriDiagonalMatrixToCompactMatrix(M, 0);
}
*/

bool isRowZero(const CompactMatrixType &M, const uint row, const double precision)
{
	if ((bool) M.QComplex)
	{
		for (uint j = 0; j < (uint) M.Nj; ++j)
			if ((Abs(M.V[row * M.Nj + j].Re) > precision) || (Abs(M.V[row * M.Nj + j].Im) > precision))
				return false;
	} else
	{
		for (uint j = 0; j < (uint) M.Nj; ++j)
			if (Abs(M.VRe[row * M.Nj + j]) > precision)
				return false;
	}
	return true;
}

bool isRowZero(const CompactMatrixType &M, const uint row)
{
	return isRowZero(M, row, 10 * DBL_EPSILON);
}

bool isColZero(const CompactMatrixType &M, const uint col, const double precision)
{
	if ((bool) M.QComplex)
	{
		for (uint i = 0; i < (uint) M.Ni; ++i)
			if ((Abs(M.V[i * M.Nj + col].Re) > precision) || (Abs(M.V[i * M.Nj + col].Im) > precision))
				return false;
	} else
	{
		for (uint i = 0; i < (uint) M.Ni; ++i)
			if (Abs(M.VRe[i * M.Nj + col]) > precision)
				return false;
	}
	return true;
}

bool isColZero(const CompactMatrixType &M, const uint col)
{
	return isColZero(M, col, 10 * DBL_EPSILON);
}

bool isZeroMatrix(const CompactMatrixType &M, const double zero)
{
	if ((bool) M.QComplex)
	{
		for (uint i = 0; i < (uint) M.Ni; ++i)
			for (uint j = 0; j < (uint) M.Nj; ++j)
				if ((Abs(M.V[i * M.Nj + j].Re) > zero) || (Abs(M.V[i * M.Nj + j].Im) > zero))
					return false;
	} else
	{
		for (uint i = 0; i < (uint) M.Ni; ++i)
			for (uint j = 0; j < (uint) M.Nj; ++j)
				if (Abs(M.VRe[i * M.Nj + j]) > zero)
					return false;
	}
	return true;
}

void FreeCompactMatrixMatrix(vector<vector<CompactMatrixType> > &matMat)
{
	for (uint i = 0; i < (uint) matMat.size(); ++i)
	{
		for (uint j = 0; j < (uint) matMat[i].size(); ++j)
		{
			if (!isNullMatrix(matMat[i][j]))
			{
				FreeCompactMatrix(&matMat[i][j]);
			}
		}
	}
}

void FreeMatrixMatrix(vector<vector<MatrixType> > &matMat)
{
	for (uint i = 0; i < (uint) matMat.size(); ++i)
		for (uint j = 0; j < (uint) matMat.size(); ++j)
			if (!isNullMatrix(matMat[i][j]))
				FreeMatrix(&matMat[i][j]);
}

void CompactMatrixCutDimensions(CompactMatrixType &M)
{
	if (M.Ni != M.Nj)
	{
		cout << "Error in Cut Dimensions: Number of rows (" << M.Ni << ") not equal to number of colums (" << M.Ni << ")." << endl;
		return;
	}
	uint newDim;
	for (newDim = M.Ni; newDim > 0; --newDim)
		if (!isRowZero(M, newDim-1) || !isColZero(M,newDim-1))
			break;
	CompactMatrixType pointerToGarbage = M;
	M = Submatrix(M, newDim);
	FreeCompactMatrix(&pointerToGarbage);
}

CompactMatrixType SubspaceRotationMatrix(const uint dim1, const CompactMatrixType &T, const uint dim2)
{
	CompactMatrixType out, one1, one2;
	one1.QComplex = 0;
	one1.Ni = dim1;
	one1.Nj = dim1;
	one2.QComplex = 0;
	one2.Ni = dim2;
	one2.Nj = dim2;
	InitCompactIdentityMatrix(&one1);
	InitCompactIdentityMatrix(&one2);
	vector < vector<CompactMatrixType> > outMatMat = CompactnullMatrixMatrix(3, 3);
	FreeCompactMatrix(&(outMatMat[0][0]));
	outMatMat[0][0] = one1;
	FreeCompactMatrix(&(outMatMat[1][1]));
	CompactMatrixCopy(T, &(outMatMat[1][1]));
	FreeCompactMatrix(&(outMatMat[2][2]));
	outMatMat[2][2] = one2;
	out = MatrixFlatten(outMatMat);
	FreeCompactMatrixMatrix(outMatMat);
	return out;
}

CompactMatrixType MatrixTypeToCompactMatrixType(const MatrixType &M)
{
	CompactMatrixType out = CompactMatrixType();
	out.Ni = M.Ni;
	out.Nj = M.Nj;
	out.QComplex = M.QComplex;
	InitCompactMatrix(&out);
	strcpy(out.Name, M.Name);
	out.a0 = M.a0;
	out.b0 = M.b0;
	out.mu = M.mu;
	for (uint i = 0; i < M.Ni; ++i)
		for (uint j = 0; j < M.Nj; ++j)
			if (M.QComplex == 0)
			{
				out.VRe[i * out.Nj + j] = M.VRe[i][j];
			} else
			{
				out.V[i * out.Nj + j].Re = M.VRe[i][j];
				out.V[i * out.Nj + j].Im = M.VIm[i][j];
			}
	return out;
}

CompactMatrixType removeSpin(const CompactMatrixType &M)
{
	CompactMatrixType out = CompactMatrixType();
	out.Ni = M.Ni / 2;
	out.Nj = M.Nj / 2;
	out.QComplex = M.QComplex;
	InitCompactMatrix(&out);
	if (out.QComplex == 0)
	{
		for (uint i = 0; i < out.Ni; ++i)
		{
			for (uint j = 0; j < out.Nj; ++j)
			{
				out.VRe[i * out.Nj + j] = M.VRe[2 * i * M.Nj + 2 * j];
			}
		}
	} else
	{
		for (uint i = 0; i < out.Ni; ++i)
		{
			for (uint j = 0; j < out.Nj; ++j)
			{
				out.V[i * out.Nj + j].Re = M.V[2 * i * M.Nj + 2 * j].Re;
				out.V[i * out.Nj + j].Im = M.V[2 * i * M.Nj + 2 * j].Im;
			}
		}
	}
	return out;
}

static int LuaRemoveSpin(lua_State *L)
{
	CompactMatrixType M;
	uint len = (uint) lua_gettop(L);
	if (len != 1)
		return luaL_error(L, "Error: LuaRemoveSpin called with %d argmuents, only 1 expected, a Matrix.", len);
	lua_check(L, 1, M);
	if (M.Ni % 2 != 0 || M.Nj % 2 != 0)
		return luaL_error(L, "Error: LuaRemoveSpin called with a matrix of dimension %d times %d. Number of rows and columns must each be even.", M.Ni, M.Nj);
	lua_push(L, removeSpin(M));
	return 1;
}

CompactMatrixType Submatrix(const CompactMatrixType &M, const vector<uint> &rows, const vector<uint> &cols)
{
	for (uint i = 0; i < (uint) rows.size(); ++i)
	{
		if (rows[i] < 0 || rows[i] >= M.Ni)
		{
			cout << "Error in Submatrix: Row entry " << rows[i] << " lies outside range [0," << (M.Ni - 1) << "]." << endl;
			throw std::exception();
		}
	}
	for (uint i = 0; i < (uint) cols.size(); ++i)
	{
		if (cols[i] < 0 || cols[i] >= M.Nj)
		{
			cout << "Error in Submatrix: Column entry " << cols[i] << " lies outside range [0," << (M.Nj - 1) << "]." << endl;
			throw std::exception();
		}
	}
	CompactMatrixType Msub;
	Msub.Ni = (uint) rows.size();
	Msub.Nj = (uint) cols.size();
	Msub.QComplex = M.QComplex;
	InitCompactMatrix(&Msub);
	Msub.OnLuaStack = M.OnLuaStack;
	for (int k = 0; k < (uint) rows.size(); ++k)
	{
		uint i = rows[k];
		for (int l = 0; l < (uint) cols.size(); ++l)
		{
			uint j = cols[l];
			if (M.QComplex == 0)
			{
				Msub.VRe[k * Msub.Nj + l] = M.VRe[i * M.Nj + j];
			}
			else
			{
				Msub.V[k * Msub.Nj + l] = M.V[i * M.Nj + j];
			}
		}
	}
	return Msub;
}

