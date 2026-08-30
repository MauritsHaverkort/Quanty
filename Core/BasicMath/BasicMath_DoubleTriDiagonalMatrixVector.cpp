/*
*     ( A[0][0]=a0  | A[0][1]=bv0 | 0           | 0           | A[0][4]=bc0 | 0           | 0           )
*   A=( A[1][0]=bv0 | A[1][1]=av1 | A[1][2]=bv1 | 0           | 0           | 0           | 0           )
*     ( 0           | A[2][1]=bv1 | A[2][2]=av2 | A[2][3]=bv2 | 0           | 0           | 0           )
*     ( 0           | 0           | A[3][2]=bv2 | A[3][3]=av3 | 0           | 0           | 0           )
*     ( A[4][0]=bc0 | 0           | 0           | 0           | A[4][4]=ac1 | A[4][5]=bc1 | 0           )
*     ( 0           | 0           | 0           | 0           | A[5][4]=bc1 | A[5][5]=ac2 | A[5][6]=bc2 )
*     ( 0           | 0           | 0           | 0           | 0           | A[6][5]=bc2 | A[6][6]=ac3 )
*
*   a0 = a0val + a0con
*
*   (A B)[i][j] = \sum_k A[i][k] B[k][j]
*
*   Note that the block A[1;;3][1;;3] and the block A[4;;6][4;;6] are normal tri-diagonal matrices
*/



/* Free a DoubleTriDiagonalMatrix M (the structure is not freed)
 */
void FreeDoubleTriDiagonalMatrix(DoubleTriDiagonalMatrixType * const M)
{
  FreeTriDiagonalMatrix(&M->val);
  FreeTriDiagonalMatrix(&M->con);
  M->mu = 0;
  M->kBT = 0;
  M->a0 = 0;
  M->b0 = 0;
  strcpy(M->Name, "free");
}

/*  Init Tri diagonal matrix
 *  M.val.N M.con.N must be set
 */
uint32_t InitDoubleTriDiagonalMatrix(DoubleTriDiagonalMatrixType * const M)
{
  strcpy(M->Name, "Matrix");
  if ((M->val.N == 0) || (M->con.N == 0))
  {
    printf("Error in InitDoubleTriDiagonalMatrix M->val.N=%i M->con.N=%i\n", M->val.N, M->con.N);
    return 2;
  };
  if (InitTriDiagonalMatrix(&M->val))
  {
    printf("InitTriDiagonalMatrix failed in InitDoubleTriDiagonalMatrix\n");
    fflush(stdout);
    return 1;
  };
  if (InitTriDiagonalMatrix(&M->con))
  {
    FreeTriDiagonalMatrix(&M->val);
    printf("InitTriDiagonalMatrix failed in InitDoubleTriDiagonalMatrix\n");
    fflush(stdout);
    return 1;
  }
  M->a0 = 0;
  M->b0 = 1;
  M->mu = 0.;
  return 0;
}

/*  Init double tri-diagonal matrix and set diagonal (a) equal to one
 *  M.val.N M.con.N must be set
 */
uint32_t InitIdentityDoubleTriDiagonalMatrix(DoubleTriDiagonalMatrixType * const M)
{
  strcpy(M->Name, "I_Matrix");
  if ((M->val.N == 0) || (M->con.N == 0))
  {
    printf("Error in InitIdentityDoubleTriDiagonalMatrix M->val.N=%i M->con.N=%i\n", M->val.N, M->con.N);
    fflush(stdout);
    return 2;
  };
  if (InitIdentityTriDiagonalMatrix(&M->val))
  {
    printf("InitIdentityTriDiagonalMatrix failed in InitIdentityDoubleTriDiagonalMatrix\n");
    fflush(stdout);
    return 1;
  };
  if (InitIdentityTriDiagonalMatrix(&M->con))
  {
    FreeTriDiagonalMatrix(&M->val);
    printf("InitIdentityTriDiagonalMatrix failed in InitIdentityDoubleTriDiagonalMatrix\n");
    fflush(stdout);
    return 1;
  }
  M->a0 = 0;
  M->b0 = 1;
  M->mu = 0.;
  return 0;
}

/* set elements with an absolute value smaller then epsilon to zero
 */
void DoubleTriDiagonalMatrixChop(DoubleTriDiagonalMatrixType * const M, double const epsilon)
{
  uint32_t i;
  for (i = 0; i < M->val.N; i++)
  {
    if ((fabs(M->val.a[i]) < epsilon))
      M->val.a[i] = 0;
    if ((fabs(M->val.b[i]) < epsilon))
      M->val.b[i] = 0;
  }
  if ((fabs(M->val.a[M->val.N]) < epsilon))
    M->val.a[M->val.N] = 0;
  for (i = 0; i < M->con.N; i++)
  {
    if ((fabs(M->con.a[i]) < epsilon))
      M->con.a[i] = 0;
    if ((fabs(M->con.b[i]) < epsilon))
      M->con.b[i] = 0;
  }
  if ((fabs(M->con.a[M->con.N]) < epsilon))
    M->con.a[M->con.N] = 0;
  if ((fabs(M->a0) < epsilon))
    M->a0 = 0;
  if ((fabs(M->b0) < epsilon))
    M->b0 = 0;
}

/*  MI=M^-1
 *  MI not initiallized on call
 */
double DoubleTriDiagonalMatrixInverseOld(DoubleTriDiagonalMatrixType const M, MatrixType * const MI)
{
  uint32_t i;
  MI->Ni = M.val.N + M.con.N + 1;
  MI->Nj = M.val.N + M.con.N + 1;
  MI->QComplex = 0;
  if (InitMatrix(MI))
  {
    printf("InitMatrix failed in DoubleTriDiagonalMatrixInverseOld\n");
    fflush(stdout);
    return -1;
  }
  MI->VRe[0][0] = M.con.a[0] + M.val.a[0];
  for (i = 1; i <= M.val.N; i++)
  {
    MI->VRe[i - 1][i] = M.val.b[i - 1];
    MI->VRe[i][i - 1] = M.val.b[i - 1];
    MI->VRe[i][i] = M.val.a[i];
  }
  i = M.val.N + 1;
  MI->VRe[0][i] = M.val.b[0];
  MI->VRe[i][0] = M.val.b[0];
  MI->VRe[i][i] = M.val.a[1];
  for (i = M.val.N + 2; i <= M.val.N + M.con.N; i++)
  {
    MI->VRe[i - 1][i] = M.val.b[i - M.val.N - 1];
    MI->VRe[i][i - 1] = M.val.b[i - M.val.N - 1];
    MI->VRe[i][i] = M.val.a[i - M.val.N];
  }
  return MatrixInverse(MI);
}

/*  M<- M + alpha B
 */
void DoubleTriDiagonalMatrixAdd(DoubleTriDiagonalMatrixType * const M, DoubleTriDiagonalMatrixType B, double alphaRe)
{
  uint32_t i;
  for (i = 0; i < M->val.N - 1; i++)
  {
    M->val.a[i] += alphaRe * B.val.a[i];
    M->val.b[i] += alphaRe * B.val.b[i];
  }
  M->val.a[M->val.N] += alphaRe * B.val.a[M->val.N];
  for (i = 0; i < M->con.N - 1; i++)
  {
    M->con.a[i] += alphaRe * B.con.a[i];
    M->con.b[i] += alphaRe * B.con.b[i];
  }
  M->con.a[M->con.N] += alphaRe * B.con.a[M->con.N];
}

/*  M<- alpha M
 */
void DoubleTriDiagonalMatrixScalarMultiply(DoubleTriDiagonalMatrixType * const M, double alphaRe)
{
  uint32_t i;
  for (i = 0; i < M->val.N; i++)
  {
    M->val.a[i] *= alphaRe;
    M->val.b[i] *= alphaRe;
  }
  M->val.a[M->val.N] *= alphaRe;
  for (i = 0; i < M->con.N; i++)
  {
    M->con.a[i] *= alphaRe;
    M->con.b[i] *= alphaRe;
  }
  M->con.a[M->con.N] *= alphaRe;
  M->a0 = 0;
  M->b0 = 1;
}

/* Create a new matrix that has dimensions
 * (1+iEndval,1+iEndval)+(1+iEndcon,1+iEndcon)
 * For the valence and conduction part of the matrix
 * and is a copy of a part of A
 * M not allocated on call of the subroutine
 */
uint32_t DoubleTriDiagonalMatrixPart(DoubleTriDiagonalMatrixType * const M, DoubleTriDiagonalMatrixType A, const uint32_t iEndval, const uint32_t iEndcon)
{
  uint32_t i;
  if ((iEndval > A.val.N) || (iEndcon > A.con.N))
    return 1;
  M->val.N = iEndval;
  M->con.N = iEndcon;
  InitDoubleTriDiagonalMatrix(M);
  M->mu = A.mu;
  M->kBT = A.kBT;
  M->a0 = A.a0;
  M->b0 = A.b0;
  for (i = 0; i < M->val.N; i++)
  {
    M->val.a[i] = A.val.a[i];
    M->val.b[i] = A.val.b[i];
  }
  M->val.a[M->val.N] = A.val.a[M->val.N];
  for (i = 0; i < M->con.N; i++)
  {
    M->con.a[i] = A.con.a[i];
    M->con.b[i] = A.con.b[i];
  }
  M->con.a[M->con.N] = A.con.a[M->con.N];
  return 0;
}

/*  M<-A
 */
uint32_t DoubleTriDiagonalMatrixCopy(DoubleTriDiagonalMatrixType * const M, DoubleTriDiagonalMatrixType A)
{
  uint32_t i;
  M->val.N = A.val.N;
  M->con.N = A.con.N;
  InitDoubleTriDiagonalMatrix(M);
  M->mu = A.mu;
  M->kBT = A.kBT;
  M->a0 = A.a0;
  M->b0 = A.b0;
  for (i = 0; i < M->val.N; i++)
  {
    M->val.a[i] = A.val.a[i];
    M->val.b[i] = A.val.b[i];
  }
  M->val.a[M->val.N] = A.val.a[M->val.N];
  for (i = 0; i < M->con.N; i++)
  {
    M->con.a[i] = A.con.a[i];
    M->con.b[i] = A.con.b[i];
  }
  M->con.a[M->con.N] = A.con.a[M->con.N];
  return 0;
}



/* y <- M x
 * y already allocated
 * y[i]=sum_j M[i][j] x[j]
 */
uint32_t RealDoubleTriDiagonalMatrixVectorOld(DoubleTriDiagonalMatrixType const M, const double *const x, double *const y)
{
  uint32_t i;
  y[0] = (M.val.a[0] + M.con.a[0]) * x[0];
  for (i = 1; i <= M.val.N; i++)
  {
    y[i] = M.val.a[i] * x[i];
  }
  for (i = M.val.N + 1; i <= M.val.N + M.con.N; i++)
  {
    y[i] = M.con.a[i - M.val.N] * x[i];
  }
  for (i = 0; i < M.val.N; i++)
  {
    y[i + 1] += M.val.b[i] * x[i];
    y[i] += M.val.b[i] * x[i + 1];
  }
  y[0] += M.con.b[0] * x[M.val.N + 1];
  y[M.val.N + 1] += M.con.b[0] * x[0];
  for (i = M.val.N + 1; i < M.val.N + M.con.N; i++)
  {
    y[i + 1] += M.con.b[i - M.val.N] * x[i];
    y[i] += M.con.b[i - M.val.N] * x[i + 1];
  }
  return 0;
}


/* y <- M x
 * y already allocated
 * y[i]=sum_j M[i][j] x[j]
 */
uint32_t ComplexDoubleTriDiagonalMatrixVectorOld(DoubleTriDiagonalMatrixType const M, const double *const xRe, const double *const xIm, double *const yRe, double *const yIm)
{
  uint32_t i;
  yRe[0] = (M.val.a[0] + M.con.a[0]) * xRe[0];
  yIm[0] = (M.val.a[0] + M.con.a[0]) * xIm[0];
  for (i = 1; i <= M.val.N; i++)
  {
    yRe[i] = M.val.a[i] * xRe[i];
    yIm[i] = M.val.a[i] * xIm[i];
  }
  for (i = M.val.N + 1; i <= M.val.N + M.con.N; i++)
  {
    yRe[i] = M.con.a[i - M.val.N] * xRe[i];
    yIm[i] = M.con.a[i - M.val.N] * xIm[i];
  }
  for (i = 0; i < M.val.N; i++)
  {
    yRe[i + 1] += M.val.b[i] * xRe[i];
    yIm[i + 1] += M.val.b[i] * xIm[i];
    yRe[i] += M.val.b[i] * xRe[i + 1];
    yIm[i] += M.val.b[i] * xIm[i + 1];
  }
  yRe[0] += M.val.b[0] * xRe[M.val.N + 1];
  yIm[0] += M.val.b[0] * xIm[M.val.N + 1];
  yRe[M.val.N + 1] += M.val.b[0] * xRe[0];
  yIm[M.val.N + 1] += M.val.b[0] * xIm[0];
  for (i = M.val.N + 1; i < M.val.N + M.con.N; i++)
  {
    yRe[i + 1] += M.val.b[i - M.val.N] * xRe[i];
    yIm[i + 1] += M.val.b[i - M.val.N] * xIm[i];
    yRe[i] += M.val.b[i - M.val.N] * xRe[i + 1];
    yIm[i] += M.val.b[i - M.val.N] * xIm[i + 1];
  }
  return 0;
}


/* y <- x M x
 * y already allocated
 * y=sum_ij x[i] M[i][j] x[j]
 */
uint32_t RealVectorDoubleTriDiagonalMatrixVectorOld(DoubleTriDiagonalMatrixType const M, const double *const x, double *const y)
{
  uint32_t i;
  *y = x[0] * (M.val.a[0] + M.con.a[0]) * x[0];
  for (i = 1; i <= M.val.N; i++)
  {
    *y += x[i] * M.val.a[i] * x[i];
  }
  for (i = M.val.N + 1; i <= M.val.N + M.con.N; i++)
  {
    *y += x[i] * M.con.a[i - M.val.N] * x[i];
  }
  for (i = 0; i < M.val.N; i++)
  {
    *y += 2 * x[i] * M.val.b[i] * x[i + 1];
  }
  *y += 2 * x[M.val.N + 1] * M.con.b[0] * x[0];
  for (i = M.val.N + 1; i < M.con.N + M.val.N; i++)
  {
    *y += 2 * x[i] * M.con.b[i - M.val.N] * x[i + 1];
  }
  return 0;
}

/* y <- x^* M x
 * y already allocated
 * y=sum_ij x[i]^* M[i][j] x[j]
 */
uint32_t ComplexVectorDoubleTriDiagonalMatrixVectorOld(DoubleTriDiagonalMatrixType const M, const double *const xRe, const double *const xIm, double *const yRe, double *const yIm)
{
  uint32_t i;
  *yRe = xRe[0] * (M.val.a[0] + M.con.a[0]) * xRe[0] + xIm[0] * (M.val.a[0] + M.con.a[0]) * xIm[0];
  *yIm = 0;
  for (i = 1; i <= M.val.N; i++)
  {
    *yRe += xRe[i] * M.val.a[i] * xRe[i] + xIm[i] * M.val.a[i] * xIm[i];
  }
  for (i = M.val.N + 1; i < M.val.N + M.con.N; i++)
  {
    *yRe += xRe[i] * M.con.a[i - M.val.N] * xRe[i] + xIm[i] * M.con.a[i - M.val.N] * xIm[i];
  }
  for (i = 0; i < M.val.N; i++)
  {
    *yRe += 2 * (xRe[i] * M.val.b[i] * xRe[i + 1] + xIm[i] * M.val.b[i] * xIm[i + 1]);
  }
  *yRe += 2 * (xRe[M.val.N + 1] * M.val.b[0] * xRe[0] + xIm[M.val.N + 1] * M.val.b[0] * xIm[0]);
  for (i = 1 + M.val.N; i < M.con.N + M.val.N; i++)
  {
    *yRe += 2 * (xRe[i] * M.val.b[i - M.val.N] * xRe[i + 1] + xIm[i] * M.val.b[i - M.val.N] * xIm[i + 1]);
  }
  return 0;
}
