
/*
 *     ( A[0][0]=avc0| A[0][1]=bv0 | 0           | 0           | A[0][4]=bc0 | 0           | 0           )
 *   A=( A[1][0]=bv0 | A[1][1]=av1 | A[1][2]=bv1 | 0           | 0           | 0           | 0           )
 *     ( 0           | A[2][1]=bv1 | A[2][2]=av2 | A[2][3]=bv2 | 0           | 0           | 0           )
 *     ( 0           | 0           | A[3][2]=bv2 | A[3][3]=av3 | 0           | 0           | 0           )
 *     ( A[4][0]=bc0 | 0           | 0           | 0           | A[4][4]=ac1 | A[4][5]=bc1 | 0           )
 *     ( 0           | 0           | 0           | 0           | A[5][4]=bc1 | A[5][5]=ac2 | A[5][6]=bc2 )
 *     ( 0           | 0           | 0           | 0           | 0           | A[6][5]=bc2 | A[6][6]=ac3 )
 *
 *   avc0 = a0val + a0con
 *
 *   (A B)[i][j] = \sum_k A[i][k] B[k][j]
 *
 *   Note that the block A[1;;3][1;;3] and the block A[4;;6][4;;6] are normal tri-diagonal matrices
 *
 *   The Green's function representation is
 *   G(w) = a0 + b0 * 1/(w - avc0 - bv0 * 1/(w - av1 - bv1 (...) bv1) * bv0 - bc0 * 1/(w - ac1 - bc1 (...) bc1) * bc0) * b0
 *   G(w) = a0 + b0 * 1/(w - Gval(w) - Gcon(w) ) b0
 *
 */

uint32_t InitBlockDoubleTriDiagonalMatrix(BlockDoubleTriDiagonalMatrixType *const A)
{
  A->val.BlockSize=A->BlockSize;
  A->val.QComplex=A->QComplex;
  InitBlockTridiagonalMatrix(&(A->val));
  A->con.BlockSize=A->BlockSize;
  A->con.QComplex=A->QComplex;
  InitBlockTridiagonalMatrix(&(A->con));
  A->a0.QComplex=A->QComplex;
  A->a0.Ni=A->BlockSize;
  A->a0.Nj=A->BlockSize;
  InitCompactMatrix(&(A->a0));
  A->b0.QComplex=A->QComplex;
  A->b0.Ni=A->BlockSize;
  A->b0.Nj=A->BlockSize;
  InitCompactIdentityMatrix(&(A->b0));
  A->kBT = 0;
  return 0;
}

uint32_t FreeBlockDoubleTriDiagonalMatrix(BlockDoubleTriDiagonalMatrixType *const A)
{
  FreeBlockTridiagonalMatrix(&(A->val));
  FreeBlockTridiagonalMatrix(&(A->con));
  FreeCompactMatrix(&(A->a0));
  FreeCompactMatrix(&(A->b0));
  return 0;
}

uint32_t BlockDoubleTriDiagonalMatrixChop(BlockDoubleTriDiagonalMatrixType *const A, double const Epsilon)
{
  //TODO:
  return 0;
}

uint32_t BlockDoubleTriDiagonalMatrixChopReduceDimension(BlockDoubleTriDiagonalMatrixType *const A, double const Epsilon)
{
  //TODO:
  return 0;
}

/*
 *  B <- A
 *  B not initialized on entry
 */
uint32_t CopyBlockDoubleTriDiagonalMatrix(BlockDoubleTriDiagonalMatrixType A, BlockDoubleTriDiagonalMatrixType * B)
{
  *B=A;
  CopyBlockTridiagonalMatrix(A.val, &(B->val));
  CopyBlockTridiagonalMatrix(A.con, &(B->con));
  CompactMatrixCopy(A.a0, &(B->a0));
  CompactMatrixCopy(A.b0, &(B->b0));
  return 0;
}
