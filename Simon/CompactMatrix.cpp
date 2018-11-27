CompactMatrixType BlockTriDiagonalMatrixToCompactMatrix(const BlockTriDiagonalMatrixType &M, const double zero)
{
	vector < vector<CompactMatrixType> > outMatMat = CompactnullMatrixMatrix(M.N - 1, M.N - 1);
	for (uint i = 1; i < M.N; ++i)
	{
		CompactMatrixType Atmp, Btmp, BtmpT;
		Atmp.QComplex = M.QComplex;
		Atmp.Ni = M.BlockSize;
		Atmp.Nj = Atmp.Ni;
		InitCompactMatrix(&Atmp);
		CompactMatrixCopy(Atmp, &Btmp);
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
		outMatMat[i - 1][i - 1] = Atmp;
		if (isZeroMatrix(Btmp, zero))
			break;
		if (i < M.N - 1)
		{
			outMatMat[i - 1][i] = Btmp;
			CompactMatrixCopy(Btmp, &BtmpT);
			if ((bool) Btmp.QComplex)
				CompactMatrixConjugateTranspose(&BtmpT);
			else
				CompactMatrixTranspose(&BtmpT);
			outMatMat[i][i - 1] = BtmpT;
		}
	}
	CompactMatrixType out = MatrixFlatten(outMatMat);
	strcpy(out.Name, M.Name);
	out.mu = M.mu;
	return out;
}

CompactMatrixType BlockTriDiagonalMatrixToCompactMatrix(const BlockTriDiagonalMatrixType &M)
{
	return BlockTriDiagonalMatrixToCompactMatrix(M, 0);
}


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
		for (uint j = 0; j < (uint) matMat.size(); ++j)
			if (!isNullMatrix(matMat[i][j]))
				FreeCompactMatrix(&matMat[i][j]);
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



