bool isRowZero(const CompactMatrixType &M, const uint32_t row, const double precision)
{
	if ((bool) M.QComplex)
	{
		for (int j = 0; j <  M.Nj; ++j)
			if ((Abs(M.V[row * M.Nj + j].Re) > precision) || (Abs(M.V[row * M.Nj + j].Im) > precision))
				return false;
	} else
	{
		for (int j = 0; j < M.Nj; ++j)
			if (Abs(M.VRe[row * M.Nj + j]) > precision)
				return false;
	}
	return true;
}

bool isRowZero(const CompactMatrixType &M, const uint32_t row)
{
	return isRowZero(M, row, 10 * DBL_EPSILON);
}

bool isColZero(const CompactMatrixType &M, const uint32_t col, const double precision)
{
	if ((bool) M.QComplex)
	{
		for (int i = 0; i < M.Ni; ++i)
			if ((Abs(M.V[i * M.Nj + col].Re) > precision) || (Abs(M.V[i * M.Nj + col].Im) > precision))
				return false;
	} else
	{
		for (int i = 0; i < M.Ni; ++i)
			if (Abs(M.VRe[i * M.Nj + col]) > precision)
				return false;
	}
	return true;
}

bool isColZero(const CompactMatrixType &M, const uint32_t col)
{
	return isColZero(M, col, 10 * DBL_EPSILON);
}

bool isZeroMatrix(const CompactMatrixType &M, const double zero)
{
	if ((bool) M.QComplex)
	{
		for (int i = 0; i < M.Ni; ++i)
			for (int j = 0; j < M.Nj; ++j)
				if ((Abs(M.V[i * M.Nj + j].Re) > zero) || (Abs(M.V[i * M.Nj + j].Im) > zero))
					return false;
	} else
	{
		for (int i = 0; i < M.Ni; ++i)
			for (int j = 0; j < M.Nj; ++j)
				if (Abs(M.VRe[i * M.Nj + j]) > zero)
					return false;
	}
	return true;
}

void FreeCompactMatrixMatrix(vector<vector<CompactMatrixType> > &matMat)
{
	for (int i = 0; i < matMat.size(); ++i)
	{
		for (int j = 0; j < matMat[i].size(); ++j)
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
	for (int i = 0; i < matMat.size(); ++i)
		for (int j = 0; j < matMat.size(); ++j)
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
	uint32_t newDim;
	for (newDim = M.Ni; newDim > 0; --newDim)
		if (!isRowZero(M, newDim-1) || !isColZero(M,newDim-1))
			break;
	CompactMatrixType pointerToGarbage = M;
	M = Submatrix(M, newDim);
	FreeCompactMatrix(&pointerToGarbage);
}

CompactMatrixType SubspaceRotationMatrix(const uint32_t dim1, const CompactMatrixType &T, const uint32_t dim2)
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
	for (int i = 0; i < M.Ni; ++i)
		for (int j = 0; j < M.Nj; ++j)
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

CompactMatrixType Submatrix(const CompactMatrixType &M, const vector<uint32_t> &rows, const vector<uint32_t> &cols)
{
	for (int i = 0; i < rows.size(); ++i)
	{
		if (rows[i] < 0 || rows[i] >= M.Ni)
		{
			cout << "Error in Submatrix: Row entry " << rows[i] << " lies outside range [0," << (M.Ni - 1) << "]." << endl;
			throw std::exception();
		}
	}
	for (int i = 0; i < cols.size(); ++i)
	{
		if (cols[i] < 0 || cols[i] >= M.Nj)
		{
			cout << "Error in Submatrix: Column entry " << cols[i] << " lies outside range [0," << (M.Nj - 1) << "]." << endl;
			throw std::exception();
		}
	}
	CompactMatrixType Msub;
	Msub.Ni = (uint32_t) rows.size();
	Msub.Nj = (uint32_t) cols.size();
	Msub.QComplex = M.QComplex;
	InitCompactMatrix(&Msub);
	Msub.OnLuaStack = M.OnLuaStack;
	for (size_t k = 0; k < rows.size(); ++k)
	{
		int i = rows[k];
		for (size_t l = 0; l < cols.size(); ++l)
		{
			int j = cols[l];
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

CompactMatrixType GSOrthogonalizeWithNonOrthoBasis(const CompactMatrixType &in, const CompactMatrixType &Overlap)
{
	if ((bool) in.QComplex or (bool) Overlap.QComplex)
	{
		cout << "Function not implemented for Complex Vectors or Overlap Matrices!" << endl;
		return CompactnullMatrix();
	}
	CompactMatrixType out;
	CompactMatrixCopy(in, &out);
	int Nj = out.Nj;
	for (int i = 0; i < out.Ni; ++i)
	{
		for (int j = 0; j < i; ++j)
		{
			double overl = 0;
			for (int k = 0; k < Nj; ++k)
			{
				for (int l = 0; l < Nj; ++l)
				{
					overl += out.VRe[j * Nj + k] * Overlap.VRe[k * Nj + l] * out.VRe[i * Nj + l];
				}
			}
			if (overl != 0)
			{
				for (int k = 0; k < Nj; ++k)
				{
					out.VRe[i * Nj + k] -= overl * out.VRe[j * Nj + k];
				}
			}
		}
		double norm = 0;
		for (int k = 0; k < Nj; ++k)
		{
			for (int l = 0; l < Nj; ++l)
			{
				norm += out.VRe[i * Nj + k] * Overlap.VRe[k * Nj + l] * out.VRe[i * Nj + l];
			}
		}
		if (norm != 1.0)
		{
			norm = 1.0 / sqrt(norm);
			for (int j = 0; j < Nj; ++j)
				out.VRe[i * Nj + j] = out.VRe[i * Nj + j] * norm;
		}
	}
	return out;
}

/*I did not want to touch Orthonormalize, although I think that on the long run this function and the one it calls should
 * be merged with the existing functions.
 * This is why this Lua Interface Function is written rather sloppily
*/
static int LuaOrthogonalizeWithNonOrthoBasis(lua_State *L)
{
	CompactMatrixType in, Overlap;
	lua_check(L, 1, in);
	lua_check(L, 2, Overlap);
	lua_push(L, GSOrthogonalizeWithNonOrthoBasis(in, Overlap));
	return 1;
}
