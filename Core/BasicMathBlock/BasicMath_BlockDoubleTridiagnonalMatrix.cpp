
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

/*  QComplex of an object whose parts can each be real or complex: 0 all four real,
 *  1 all four complex, 2 a mixture. Same convention as SetQComplexOperator and as the
 *  lua reader of this type.
 *
 *  A chop can demote one part and not another -- CompactMatrixChop and
 *  BlockTridiagonalMatrixChop each call their own ...ToReal when nothing imaginary
 *  survives -- so after chopping, the container cannot hold a single answer any more.
 *  2 says "ask each part". Every consumer of this type does: Evaluate, Copy, ToSpectra,
 *  ToCompactMatrix and the lua pusher all go through EvaluateBlockTriDiagonalMatrix,
 *  CopyBlockTridiagonalMatrix and the CompactMatrix routines, each of which reads the
 *  QComplex of the operand it is handed. The only place that reads the container's own
 *  flag is ResponseFunctionGetA0, which narrows it with a != 0 test and so treats a
 *  mixture as complex -- conservative and correct.
 */
uint32_t SetBlockDoubleTriDiagonalMatrixQComplex(BlockDoubleTriDiagonalMatrixType *const A)
{
  uint32_t ncomplex = (A->val.QComplex!=0) + (A->con.QComplex!=0) + (A->a0.QComplex!=0) + (A->b0.QComplex!=0);
  if(ncomplex==0)      A->QComplex = 0;
  else if(ncomplex==4) A->QComplex = 1;
  else                 A->QComplex = 2;
  return 0;
}

/*  true if the chop left nothing of b0 at all */
bool BlockDoubleTriDiagonalMatrixB0IsZero(BlockDoubleTriDiagonalMatrixType const *const A)
{
  uint32_t i;
  for(i=0;i<A->b0.Ni*A->b0.Nj;i++)
  {
    if(A->b0.QComplex)
    {
      if((A->b0.V[i].Re!=0)||(A->b0.V[i].Im!=0)) return false;
    }
    else
    {
      if(A->b0.VRe[i]!=0) return false;
    }
  }
  return true;
}

/*  set elements with an absolute value smaller than Epsilon to zero, keeping the
 *  dimensions. The two chains are ordinary block tri-diagonal matrices and a0 and b0 are
 *  ordinary compact matrices, so all four are chopped by their own routine; each of those
 *  demotes its own storage to real if nothing imaginary is left, which is why QComplex
 *  has to be recomputed afterwards and may come out mixed.
 */
uint32_t BlockDoubleTriDiagonalMatrixChop(BlockDoubleTriDiagonalMatrixType *const A, double const Epsilon)
{
  if(BlockTridiagonalMatrixChop(&(A->val), Epsilon))
  {
    printf("BlockTridiagonalMatrixChop failed for val in BlockDoubleTriDiagonalMatrixChop\n");
    fflush(stdout);
    return 1;
  }
  if(BlockTridiagonalMatrixChop(&(A->con), Epsilon))
  {
    printf("BlockTridiagonalMatrixChop failed for con in BlockDoubleTriDiagonalMatrixChop\n");
    fflush(stdout);
    return 1;
  }
  if(CompactMatrixChop(&(A->a0), Epsilon) || CompactMatrixChop(&(A->b0), Epsilon))
  {
    printf("CompactMatrixChop failed in BlockDoubleTriDiagonalMatrixChop\n");
    fflush(stdout);
    return 1;
  }
  SetBlockDoubleTriDiagonalMatrixQComplex(A);
  return 0;
}

/*  as above, and shorten what a chopped hopping has disconnected. The three ways a
 *  dimension collapses are the ones the scalar type has, see
 *  DoubleTriDiagonalMatrixChopReduceDimension: a chain cut at a small b, a branch
 *  disconnected by a small b[0], and a b0 chopped away entirely, which leaves G(w) = a0
 *  with no matrix behind it and so no chains either.
 */
/*  Cut both chains, the valence one to at most Nval and the conduction one to at most
 *  Ncon. a0 and b0 sit in front of the two chains and are left alone.
 */
uint32_t BlockDoubleTriDiagonalMatrixTruncate(BlockDoubleTriDiagonalMatrixType *const A, uint32_t const Nval, uint32_t const Ncon)
{
  if(BlockTridiagonalMatrixTruncate(&(A->val), Nval))
  {
    printf("BlockTridiagonalMatrixTruncate failed for val in BlockDoubleTriDiagonalMatrixTruncate\n");
    fflush(stdout);
    return 1;
  }
  if(BlockTridiagonalMatrixTruncate(&(A->con), Ncon))
  {
    printf("BlockTridiagonalMatrixTruncate failed for con in BlockDoubleTriDiagonalMatrixTruncate\n");
    fflush(stdout);
    return 1;
  }
  return 0;
}

uint32_t BlockDoubleTriDiagonalMatrixChopReduceDimension(BlockDoubleTriDiagonalMatrixType *const A, double const Epsilon)
{
  if(BlockTridiagonalMatrixChopReduceDimension(&(A->val), Epsilon))
  {
    printf("BlockTridiagonalMatrixChopReduceDimension failed for val in BlockDoubleTriDiagonalMatrixChopReduceDimension\n");
    fflush(stdout);
    return 1;
  }
  if(BlockTridiagonalMatrixChopReduceDimension(&(A->con), Epsilon))
  {
    printf("BlockTridiagonalMatrixChopReduceDimension failed for con in BlockDoubleTriDiagonalMatrixChopReduceDimension\n");
    fflush(stdout);
    return 1;
  }
  if(CompactMatrixChop(&(A->a0), Epsilon) || CompactMatrixChop(&(A->b0), Epsilon))
  {
    printf("CompactMatrixChop failed in BlockDoubleTriDiagonalMatrixChopReduceDimension\n");
    fflush(stdout);
    return 1;
  }
  if(BlockDoubleTriDiagonalMatrixB0IsZero(A))
  {
    A->val.N = 0;
    A->con.N = 0;
  }
  SetBlockDoubleTriDiagonalMatrixQComplex(A);
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
