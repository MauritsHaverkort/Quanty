/*
*     ( A[0][0]=acv0  | A[0][1]=bv0   | 0             | 0           | A[0][4]=bc0   | 0             | 0           )
*   A=( A[1][0]=bv0^+ | A[1][1]=av1   | A[1][2]=bv1   | 0           | 0             | 0             | 0           )
*     ( 0             | A[2][1]=bv1^+ | A[2][2]=av2   | A[2][3]=bv2 | 0             | 0             | 0           )
*     ( 0             | 0             | A[3][2]=bv2^+ | A[3][3]=av3 | 0             | 0             | 0           )
*     ( A[4][0]=bc0^+ | 0             | 0             | 0           | A[4][4]=ac1   | A[4][5]=bc1   | 0           )
*     ( 0             | 0             | 0             | 0           | A[5][4]=bc1^+ | A[5][5]=ac2   | A[5][6]=bc2 )
*     ( 0             | 0             | 0             | 0           | 0             | A[6][5]=bc2^+ | A[6][6]=ac3 )
*
*   acv0 = (a0val + a0con)/2 
*
*   (A B)[i][j] = \sum_k A[i][k] . B[k][j]
*
*   Note that the block A[1;;3][1;;3] and the block A[4;;6][4;;6] are normal blok-tri-diagonal matrices
*/



/*  Remove M from memory, free each matrix M.a[i] and M.b[i] as well
 *  as the list of pointers M.a and M.b
 */
void FreeBlokDoubleTriDiagonalMatrix(BlokDoubleTriDiagonalMatrixType * const M)
{
  FreeBlokTriDiagonalMatrix(&M->val);
  FreeBlokTriDiagonalMatrix(&M->con);
  M->mu = 0;
  M->Epsilon = 0;
  strcpy(M->Name, "free");
}

/*  Free the list of pointers M.val.a M.val.b M.con.a and M.con.b these pointers
 *  should not point to unfreed memory, as then there is a memory leak
 */
void FreeFrameBlokDoubleTriDiagonalMatrix(BlokDoubleTriDiagonalMatrixType * const M)
{
  FreeFrameBlokTriDiagonalMatrix(&M->val);
  FreeFrameBlokTriDiagonalMatrix(&M->con);
  M->mu = 0;
  M->Epsilon = 0;
  strcpy(M->Name, "free");
}

/* Initiallize a matrix M must already be allocated 
*  M.val.N and M.con.N must be set
*  bloksize gives the dimension of the matrices that build the bloks
*/
uint32_t InitBlokDoubleTriDiagonalMatrix(BlokDoubleTriDiagonalMatrixType * const M, uint32_t bloksize,
                                         uint32_t QComplex)
{
  strcpy(M->Name, "Matrix");
  if ((M->val.N == 0) || (M->con.N == 0))
    return 2;
  if (InitBlokTriDiagonalMatrix(&M->val, bloksize, QComplex))
    return 1;
  if (InitBlokTriDiagonalMatrix(&M->con, bloksize, QComplex))
  {
    FreeBlokTriDiagonalMatrix(&M->val);
    return 1;
  }
  return 0;
}

/* Initiallize a matrix M must already be allocated 
*  M.N must be set
*/
uint32_t InitFrameBlokDoubleTriDiagonalMatrix(BlokDoubleTriDiagonalMatrixType * const M)
{
  strcpy(M->Name, "Matrix");
  if ((M->val.N == 0) || (M->con.N == 0))
    return 2;
  if (InitFrameBlokTriDiagonalMatrix(&M->val))
    return 1;
  if (InitFrameBlokTriDiagonalMatrix(&M->con))
  {
    FreeFrameBlokTriDiagonalMatrix(&M->val);
    return 1;
  }
  return 0;
}

/* Initiallize a matrix M must already be allocated 
*  M.N must be set
*  bloksize gives the dimension of the matrices that build the bloks
*/
uint32_t InitIdentityBlokDoubleTriDiagonalMatrix(BlokDoubleTriDiagonalMatrixType * const M, uint32_t bloksize,
                                                 uint32_t QComplex)
{
  strcpy(M->Name, "I_Matrix");
  if ((M->val.N == 0) || (M->con.N == 0))
    return 2;
  if (InitIdentityBlokTriDiagonalMatrix(&M->val, bloksize, QComplex))
    return 1;
  if (InitIdentityBlokTriDiagonalMatrix(&M->con, bloksize, QComplex))
  {
    FreeBlokTriDiagonalMatrix(&M->val);
    return 1;
  }
  return 0;
}

/*  set all ellements with absolute values smaller then epsilon to zero
 */
uint32_t BlokDoubleTriDiagonalMatrixChop(BlokDoubleTriDiagonalMatrixType * const M, double const epsilon)
{
  BlokTriDiagonalMatrixChop(&M->val, epsilon);
  BlokTriDiagonalMatrixChop(&M->con, epsilon);
  return 0;
}

/*
*  M<- M + alpha B
*/
uint32_t BlokDoubleTriDiagonalMatrixAdd(BlokDoubleTriDiagonalMatrixType * const M, BlokDoubleTriDiagonalMatrixType B,
                                        double const alphaRe)
{
  BlokTriDiagonalMatrixAdd(&M->val, B.val, alphaRe);
  BlokTriDiagonalMatrixAdd(&M->con, B.con, alphaRe);
  return 0;
}

/*
*  M<- alpha M
*/
uint32_t BlokDoubleTriDiagonalMatrixScalarMultiply(BlokDoubleTriDiagonalMatrixType * const M, double alphaRe)
{
  BlokTriDiagonalMatrixScalarMultiply(&M->val, alphaRe);
  BlokTriDiagonalMatrixScalarMultiply(&M->con, alphaRe);
  return 0;
}

/*  returns a new blok tridiagonal matrix consisting of bloks 0 to iEndcon and iEndval
 */
uint32_t BlokDoubleTriDiagonalMatrixPart(BlokDoubleTriDiagonalMatrixType * const M, BlokDoubleTriDiagonalMatrixType A,
                                         uint32_t iEndval, uint32_t iEndcon)
{
  M->mu = A.mu;
  M->Epsilon = A.Epsilon;
  strcpy(M->Name, A.Name);
  if ((iEndval >= A.val.N) || (iEndcon >= A.con.N))
    return 1;
  if (BlokTriDiagonalMatrixPart(&M->val, A.val, 0, iEndval))
    return 1;
  if (BlokTriDiagonalMatrixPart(&M->con, A.con, 0, iEndcon))
  {
    FreeBlokTriDiagonalMatrix(&M->val);
    return 1;
  }
  return 0;
}

/*
*  M <- A
*/
uint32_t BlokDoubleTriDiagonalMatrixCopy(BlokDoubleTriDiagonalMatrixType * const M, BlokDoubleTriDiagonalMatrixType A)
{
  M->mu = A.mu;
  M->Epsilon = A.Epsilon;
  strcpy(M->Name, A.Name);
  if (BlokTriDiagonalMatrixCopy(&M->val, A.val))
    return 1;
  if (BlokTriDiagonalMatrixCopy(&M->con, A.con))
  {
    FreeBlokTriDiagonalMatrix(&M->val);
    return 1;
  }
  return 0;
}

/*
*  M <- M^T
*  M.a[i] <- M.a[i]^T
*  M.b[i] <- M.b[i]^*
*/
uint32_t BlokDoubleTriDiagonalMatrixTranspose(BlokDoubleTriDiagonalMatrixType * const M)
{
  BlokTriDiagonalMatrixTranspose(&M->val);
  BlokTriDiagonalMatrixTranspose(&M->con);
  return 0;
}


/*
*  M <- M^*
*  M.a[i] <- M.a[i]^*
*  M.b[i] <- M.b[i]^*
*/
uint32_t BlokDoubleTriDiagonalMatrixComplexConjugate(BlokDoubleTriDiagonalMatrixType * const M)
{
  BlokTriDiagonalMatrixComplexConjugate(&M->val);
  BlokTriDiagonalMatrixComplexConjugate(&M->con);
  return 0;
}


/*
*  M <- M^+
*  M.a[i] <- M.a[i]^+
*  M.b[i] <- M.b[i]
*/
uint32_t BlokDoubleTriDiagonalMatrixConjugateTranspose(BlokDoubleTriDiagonalMatrixType * const M)
{
  BlokTriDiagonalMatrixConjugateTranspose(&M->val);
  BlokTriDiagonalMatrixConjugateTranspose(&M->con);
  return 0;
}

/* y <- M x
 * y already allocated
 * y[i]=sum_j M[i][j] x[j]
 */
uint32_t RealBlokDoubleTriDiagonalMatrixVector(BlokDoubleTriDiagonalMatrixType const M, double *const x,
                                               double *const y)
{
  int const Nblok = M.val.a[0].Ni;
  double temp[Nblok];
  int    i, Nval = 0, Ncon = 0;
  for (i = 1; i < M.val.N; i++)
  {
    Nval += M.val.a[i].Ni;
  }
  for (i = 1; i < M.con.N; i++)
  {
    Ncon += M.con.a[i].Ni;
  }
  for (i = 0; i < Nblok; i++)
  {
    temp[i] = x[i + Nval];
    x[i + Nval] = 0;
  }
  RealBlokTriDiagonalMatrixVector(M.con, &(x[Nval]), &(y[Nval]));
  for (i = 0; i < Nblok; i++)
  {
    x[i + Nval] = temp[i];
    y[i + Nval] = 0;
  }
  RealBlokTriDiagonalMatrixVector(M.val, x, y);
  RealMatrixVector(M.val.b[0], temp, &x[Nval + Nblok]);
  for (i = 0; i < M.val.b[0].Nj; i++)
  {
    y[i] += temp[i];
  }
  return 0;
}


/* y <- M x
 * y already allocated
 * y[i]=sum_j M[i][j] x[j]
 */
uint32_t ComplexBlokDoubleTriDiagonalMatrixVector(BlokDoubleTriDiagonalMatrixType const M, double *const xRe,
                                                  double *const xIm, double *const yRe, double *const yIm)
{
  int const Nblok = M.val.a[0].Ni;
  double tempRe[Nblok], tempIm[Nblok];
  int    i, Nval = 0, Ncon = 0;
  for (i = 1; i < M.val.N; i++)
  {
    Nval += M.val.a[i].Ni;
  }
  for (i = 1; i < M.con.N; i++)
  {
    Ncon += M.con.a[i].Ni;
  }
  for (i = 0; i < Nblok; i++)
  {
    tempRe[i] = xRe[i + Nval];
    tempIm[i] = xIm[i + Nval];
    xRe[i + Nval] = 0;
    xIm[i + Nval] = 0;
  }
  ComplexBlokTriDiagonalMatrixVector(M.con, &(xRe[Nval]), &(xIm[Nval]), &(yRe[Nval]), &(yIm[Nval]));
  for (i = 0; i < Nblok; i++)
  {
    xRe[i + Nval] = tempRe[i];
    xIm[i + Nval] = tempIm[i];
    yRe[i + Nval] = 0;
    yIm[i + Nval] = 0;
  }
  ComplexBlokTriDiagonalMatrixVector(M.val, xRe, xIm, yRe, yIm);
  ComplexMatrixVector(M.val.b[0], tempRe, tempIm, &xRe[Nval + Nblok], &xIm[Nval + Nblok]);
  for (i = 0; i < M.val.b[0].Nj; i++)
  {
    yRe[i] += tempRe[i];
    yIm[i] += tempIm[i];
  }
  return 0;
}

/* y <- x M x
 * y already allocated
 * y=sum_ij x[i] M[i][j] x[j]
 */
uint32_t RealVectorBlokDoubleTriDiagonalMatrixVector(BlokDoubleTriDiagonalMatrixType const M, const double *const x,
                                                     double *const y)
{

  return 0;
}

/* y <- x^* M x
 * y already allocated
 * y=sum_ij x[i]^* M[i][j] x[j]
 */
uint32_t ComplexVectorBlokDoubleTriDiagonalMatrixVector(BlokDoubleTriDiagonalMatrixType const M,
                                                        const double *const xRe, const double *const xIm,
                                                        double *const yRe, double *const yIm)
{
  return 0;
}
