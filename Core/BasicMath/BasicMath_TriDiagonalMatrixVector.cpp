/*
*     ( A[0][0]=a1 | A[0][1]=b1 | 0          | 0          )
*   A=( A[1][0]=b1 | A[1][1]=a2 | A[1][2]=b2 | 0          )
*     ( 0          | A[2][1]=b2 | A[2][2]=a3 | A[2][3]=b3 )
*     ( 0          | 0          | A[3][2]=b3 | A[3][3]=a4 )
*
*   (A B)[i][j] = \sum_k A[i][k] B[k][j]
*
*/

/* Free a TriDiagonalMatrix M (the structure is not freed)
 */
void FreeTriDiagonalMatrix(TriDiagonalMatrixType * const M)
{
  free(M->a);
  free(M->b);
}

/*  Init Tri diagonal matrix
 *  M.N must be set
 */
uint32_t InitTriDiagonalMatrix(TriDiagonalMatrixType * const M)
{
  strcpy(M->Name, "Matrix");
  M->a = (double*)malloc((M->N + 1) * sizeof(double));
  if (M->a == NULL)
  {
    printf("malloc for M->a failed in InitTriDiagonalMatrix\n");
    fflush(stdout);
    return 1;
  };
  M->b = (double*)malloc(M->N * sizeof(double));
  if (M->b == NULL)
  {
    free(M->a);
    printf("malloc for M->b failed in InitTriDiagonalMatrix\n");
    fflush(stdout);
    return 1;
  }
  M->a[0] = 0;
  if(M->N!=0) M->b[0] = 1;
  M->mu = 0;
  return 0;
}

/*  Init tri-diagonal matrix and set diagonal (a) equal to one
 *  M.N must be set
 *
 *  N = 0 is allowed and means there is no chain: a[0], the additive offset, is
 *  the only element, and there is no site energy to set to one and no hopping.
 *  Same as InitTriDiagonalMatrix, which also allocates a[0] and no b for N = 0.
 */
uint32_t InitIdentityTriDiagonalMatrix(TriDiagonalMatrixType * const M)
{
  uint32_t i;
  strcpy(M->Name, "I_Matrix");
  M->a = (double*)malloc((M->N + 1) * sizeof(double));
  if (M->a == NULL)
  {
    printf("malloc for M->a failed in InitTriDiagonalMatrix\n");
    fflush(stdout);
    return 1;
  };
  M->b = (double*)calloc((M->N), sizeof(double));
  if (M->b == NULL)
  {
    free(M->a);
    printf("malloc for M->b failed in InitTriDiagonalMatrix\n");
    fflush(stdout);
    return 1;
  }
  for (i = 1; i <= M->N; i++)
  {
    M->a[i] = 1;
  }
  M->a[0] = 0;
  if(M->N!=0) M->b[0] = 1;
  M->mu = 0;
  return 0;
}

/* set elements with an absolute value smaller then epsilon to zero
 */
void TriDiagonalMatrixChop(TriDiagonalMatrixType * const M, double const epsilon)
{
  uint32_t i;
  for (i = 0; i < M->N; i++)
  {
    if ((fabs(M->a[i]) < epsilon)) M->a[i] = 0;
    if ((fabs(M->b[i]) < epsilon)) M->b[i] = 0;
  }
  if ((fabs(M->a[M->N]) < epsilon)) M->a[M->N] = 0;
}
void TriDiagonalMatrixChopReduceDimension(TriDiagonalMatrixType * const M, double const epsilon)
{
  uint32_t i;
  for (i = 0; i < M->N; i++)
  {
    if ((fabs(M->a[i]) < epsilon)) M->a[i] = 0;
    if ((fabs(M->b[i]) < epsilon)) M->N=i;
  }
  if ((fabs(M->a[M->N]) < epsilon)) M->a[M->N] = 0;
}

/*
 *
 */
void ShiftTriDiagonalMatrix(TriDiagonalMatrixType * const M, double const delta)
{
  uint32_t i;
  for (i = 1; i <= M->N; i++)
  {
    M->a[i] += delta;
  }
}

/*  MI=M^-1
 *  MI not initiallized on call
 */
double TriDiagonalMatrixInverse(TriDiagonalMatrixType const M, MatrixType * const MI)
{
  uint32_t i;
  MI->Ni = M.N;
  MI->Nj = M.N;
  MI->QComplex = 0;
  if (InitMatrix(MI))
  {
    printf("InitMatrix(MI) failed in TriDiagonalMatrixInverse\n");
    fflush(stdout);
    return -1;
  };
  for (i = 1; i < M.N; i++)
  {
    MI->VRe[i - 1][i] = M.b[i];
    MI->VRe[i][i - 1] = M.b[i];
    MI->VRe[i][i] = M.a[i + 1];
  }
  MI->VRe[0][0] = M.a[1];
  return MatrixInverse(MI);
}

/* Given G = a0 + b0 1/(w - a1 - b1 (1/(w - ...)) b1 ) b0
 * return
 *  a1 + b1 (1/(w - ...)) b1
 */
int TriDiagonalMatrixEffectiveHamiltonian(TriDiagonalMatrixType const * const G, TriDiagonalMatrixType * const Heff)
{
  if (G->N==0)
  {
    printf("TriDiagonalMatrixEffectiveHamiltonian failed with on response function with length 0");
    Heff->N=0;
    InitTriDiagonalMatrix(Heff);
    return 1;
  }
  Heff->N = G->N-1;
  InitTriDiagonalMatrix(Heff);
  Heff->mu=G->mu;
  memcpy(Heff->a, &(G->a[1]), (Heff->N+1)*sizeof(double));
  memcpy(Heff->b, &(G->b[1]), (Heff->N  )*sizeof(double));
  return 0;
}
int TriDiagonalMatrixFromEffectiveHamiltonian(TriDiagonalMatrixType const * const Heff, TriDiagonalMatrixType * const G)
{
  G->N = Heff->N+1;
  InitTriDiagonalMatrix(G);
  G->mu=Heff->mu;
  memcpy(&(G->a[1]), Heff->a, (Heff->N+1)*sizeof(double));
  memcpy(&(G->b[1]), Heff->b, (Heff->N  )*sizeof(double));
  return 0;
}


/*  M<- M + alpha B
 */
void TriDiagonalMatrixAdd(TriDiagonalMatrixType * const M, TriDiagonalMatrixType B, double alphaRe)
{
  uint32_t i;
  for (i = 1; i < M->N; i++)
  {
    M->a[i] += alphaRe * B.a[i];
    M->b[i] += alphaRe * B.b[i];
  }
  M->a[M->N] += alphaRe * B.a[M->N];
  M->a[0] = 0;
  M->b[0] = 1;
}

/*  M<- alpha M
 */
void TriDiagonalMatrixScalarMultiply(TriDiagonalMatrixType * const M, double alphaRe)
{
  uint32_t i;
  for (i = 1; i < M->N; i++)
  {
    M->a[i] *= alphaRe;
    M->b[i] *= alphaRe;
  }
  M->a[M->N] *= alphaRe;
  M->a[0] = 0;
  M->b[0] = 1;
}

/* Create a new matrix that has dimensions (1+iEnd-iStart,1+iEnd-iStart)
 * and is a copy of a part of A
 * namely A[{iStart;;iEnd}][{iStart;;iEnd}]
 * a[0] and b[0] are set to A.a[iStart-1] and A.b[iStart-1], the on site energy of the
 * last site before the part and the hopping connecting it to the part, such that the
 * result evaluates as the Green's function of the part as seen from site iStart-1.
 * For iStart=1 this reduces to the offset a[0] and the scaling b[0] of A.
 * M not allocated on call of the subroutine
 */
uint32_t TriDiagonalMatrixPart(TriDiagonalMatrixType * const M, TriDiagonalMatrixType A, const uint32_t iStart, const uint32_t iEnd)
{
  uint32_t i;
  if ((iEnd < iStart) || (iEnd > A.N) || (iStart == 0))
  {
    printf("TriDiagonalMatrixPart called with iStart=%i, iEnd=%i, and A.N=%i.\nThat's not an option\n", iStart, iEnd,A.N);
    return 1;
  }
  M->N = 1 + iEnd - iStart;
  InitTriDiagonalMatrix(M);
  M->a[0] = A.a[iStart-1];
  M->b[0] = A.b[iStart-1];
  for (i = 1; i < M->N; i++)
  {
    M->a[i] = A.a[i + iStart - 1];
    M->b[i] = A.b[i + iStart - 1];
  }
  M->mu = A.mu;
  M->a[M->N] = A.a[iEnd];
  return 0;
}
/*  as TriDiagonalMatrixPart, but each site is doubled in order to represent spin up and spin down
 */
uint32_t TriDiagonalMatrixPartAddSpin(BlockTriDiagonalMatrixType * const M, TriDiagonalMatrixType A, const uint32_t iStart, const uint32_t iEnd)
{
  uint32_t i;
  if ((iEnd < iStart) || (iEnd > A.N) || (iStart == 0))
  {
    printf("TriDiagonalMatrixPartAddSpin called with iStart=%i, iEnd=%i, and A.N=%i.\nThat's not an option\n", iStart, iEnd,A.N);
    return 1;
  }
  M->N = 1 + iEnd - iStart;
  M->QComplex=0;
  M->BlockSize=2;
  InitBlockTridiagonalMatrix(M);
  M->aRe[0] = A.a[iStart-1];
  M->bRe[0] = A.b[iStart-1];
  M->aRe[3] = A.a[iStart-1];
  M->bRe[3] = A.b[iStart-1];
  for (i = 1; i < M->N; i++)
  {
    M->aRe[4*i]   = A.a[i + iStart - 1];
    M->bRe[4*i]   = A.b[i + iStart - 1];
    M->aRe[4*i+3] = A.a[i + iStart - 1];
    M->bRe[4*i+3] = A.b[i + iStart - 1];
  }
  M->mu = A.mu;
  M->aRe[4*M->N]   = A.a[iEnd];
  M->aRe[4*M->N+3] = A.a[iEnd];
  return 0;
}

/*  M<-A
 */
uint32_t TriDiagonalMatrixCopy(TriDiagonalMatrixType * const M, TriDiagonalMatrixType A)
{
  uint32_t i;
  M->N = A.N;
  if (InitTriDiagonalMatrix(M))
  {
    printf("InitTriDiagonalMatrix failed in TriDiagonalMatrixCopy\n");
    fflush(stdout);
    return 1;
  };
  for (i = 0; i < M->N; i++)
  {
    M->a[i] = A.a[i];
    M->b[i] = A.b[i];
  }
  M->mu = A.mu;
  M->a[M->N] = A.a[M->N];
  return 0;
}

/* y <- M x
 * y already allocated
 * y[i]=sum_j M[i][j] x[j]
 */
void RealTriDiagonalMatrixVector(TriDiagonalMatrixType const M, const double *const x, double *const y)
{
  uint32_t i;
  for (i = 0; i < M.N; i++)
  {
    y[i] = M.a[i + 1] * x[i];
  }
  for (i = 1; i < M.N; i++)
  {
    y[i - 1] += M.b[i] * x[i];
    y[i] += M.b[i] * x[i - 1];
  }
}

/* y <- M x
 * y already allocated
 * y[i]=sum_j M[i][j] x[j]
 */
void ComplexTriDiagonalMatrixVector(TriDiagonalMatrixType const M, const double *const xRe, const double *const xIm, double *const yRe, double *const yIm)
{
  uint32_t i;
  for (i = 0; i < M.N; i++)
  {
    yRe[i] = M.a[i + 1] * xRe[i];
    yIm[i] = M.a[i + 1] * xIm[i];
  }
  for (i = 1; i < M.N; i++)
  {
    yRe[i - 1] += M.b[i] * xRe[i];
    yRe[i] += M.b[i] * xRe[i - 1];
    yIm[i - 1] += M.b[i] * xIm[i];
    yIm[i] += M.b[i] * xIm[i - 1];
  }
}

/* y <- x M x
 * y already allocated
 * y=sum_ij x[i] M[i][j] x[j]
 */
void RealVectorTriDiagonalMatrixVector(TriDiagonalMatrixType const M, const double *const x, double *const y)
{
  uint32_t i;
  *y = 0;
  for (i = 0; i < M.N; i++)
  {
    *y += x[i] * M.a[i + 1] * x[i];
  }
  for (i = 1; i < M.N; i++)
  {
    *y += 2 * x[i] * M.b[i] * x[i - 1];
  }
}

/* y <- x^* M x
 * y already allocated
 * y=sum_ij x[i]^* M[i][j] x[j]
 */
void ComplexVectorTriDiagonalMatrixVector(TriDiagonalMatrixType const M, const double *const xRe, const double *const xIm, double *const yRe, double *const yIm)
{
  uint32_t i;
  *yRe = 0;
  *yIm = 0;
  for (i = 0; i < M.N; i++)
  {
    *yRe += xRe[i] * M.a[i + 1] * xRe[i] + xIm[i] * M.a[i + 1] * xIm[i];
  }
  for (i = 1; i < M.N; i++)
  {
    *yRe += 2 * (xRe[i] * M.b[i] * xRe[i - 1] + xIm[i] * M.b[i] * xIm[i - 1]);
  }
}

/*  by Robert
 *  multiply b[blow] to b[bhigh] of a tridiagonal matrix
 */
uint32_t TriDiagonalMatrixBProduct(TriDiagonalMatrixType *SpecTri, uint32_t blow, uint32_t bhigh, double *bProduct)
{
  int i;
  
  if(blow > 0 && bhigh < SpecTri->N)
  {
    if(bhigh < blow)
      *bProduct = 1.0;//SpecTri->b[blow];
    else
    {
      *bProduct = 1.0;
      for(i = blow; i <= bhigh; i++)
        *bProduct *= SpecTri->b[i];
    }
    return 0;
  }
  else
    return 1;
}

/*  by Robert
 *  return the inverse of a tridiagonal matrix
 */
uint32_t TriDiagonalMatrixToGreensFunction(TriDiagonalMatrixType * const SpecTri, double wRes, double energy, double gamma, MatrixType * const invMat)
{
  double thetaRe[SpecTri->N+1], thetaIm[SpecTri->N+1], phiRe[SpecTri->N + 1], phiIm[SpecTri->N + 1], bProduct, denom, shift;
  int i, j;
  
  /* flip sign of tridiagonal matrix (because we want inv(E+w-H+igamma))*/
  //  TriDiagonalMatrixScalarMultiply(SpecTri, -1.0);
  //  /* shift by E+w if necessary*/
  //  if((energy + wRes) != 0)
  //    ShiftTriDiagonalMatrix(SpecTri, energy + wRes);
  shift = energy + wRes;
  /* Initialize theta and phi at start of recurrence. Note that at
   
   http://en.wikipedia.org/wiki/Tridiagonal_matrix#Inversion
   
   theta indices are defined from 0...N (N+1 values), and that is how they
   are here as well. But the matrix indices in the algorithm start at 1,1 so
   the theta indices are shifted w.r.t. the matrix indices during the main
   loop later.  The used range of Phi is defined from 2 ... N+1 (N values),
   so its indices are shifted as well below.
   */
  thetaRe[0] = 1; thetaIm[0] = 0;
  thetaRe[1] = (shift - SpecTri->a[1]); thetaIm[1] = gamma;
  
  phiRe[SpecTri->N-1] = (shift - SpecTri->a[SpecTri->N]); phiIm[SpecTri->N-1] = gamma;
  phiRe[SpecTri->N] = 1.0; phiIm[SpecTri->N] = 0.0;
  /* Calculate the theta and phi recurrences */
  for(i = 2; i < SpecTri->N+1; i++)
  {
    /*Spectri->b should be multiplied by -1 as well, but it cancels out because always squared*/
    thetaRe[i] = (shift - SpecTri->a[i])*thetaRe[i-1] - gamma*thetaIm[i-1] - SpecTri->b[i-1]*SpecTri->b[i-1]*thetaRe[i-2];
    thetaIm[i] = (shift - SpecTri->a[i])*thetaIm[i-1] + gamma * thetaRe[i-1] - SpecTri->b[i-1]*SpecTri->b[i-1]*thetaIm[i-2];
  }
  
  for(i = SpecTri->N - 2; i >= 0; i--)
  {
    phiRe[i] = (shift - SpecTri->a[i+1])*phiRe[i+1] - gamma*phiIm[i+1] - SpecTri->b[i+1]*SpecTri->b[i+1]*phiRe[i+2];
    phiIm[i] = (shift - SpecTri->a[i+1])*phiIm[i+1] + gamma*phiRe[i+1] - SpecTri->b[i+1]*SpecTri->b[i+1]*phiIm[i+2];
  }
  
  /* Compute the matrix inverse */
  invMat->QComplex = 1;
  invMat->Ni = SpecTri->N;
  invMat->Nj = SpecTri->N;
  InitMatrix(invMat);
  for(i = 0; i < SpecTri->N; i++)
  {
    for(j = i; j < SpecTri->N; j++)
    {
      if(TriDiagonalMatrixBProduct(SpecTri, i+1, j - 1+1, &bProduct))
      {
        printf("TriDiagonalMatrixBProduct failed in TriDiagonalMatrixToGreensFunction\n");
        return 1;
      }
      invMat->VRe[i][j] = 0;
      invMat->VIm[i][j] = 0;
      /*Watch for divide by zeros. Shouldn't be a problem for nonsingular matrices.*/
      denom = thetaRe[SpecTri->N]*thetaRe[SpecTri->N] + thetaIm[SpecTri->N]*thetaIm[SpecTri->N];
      if(denom != 0)
      {
        /*Doing real and imaginary parts separately for:
         
         bProduct*theta[i]*phi[j+1]/theta[N]
         */
        invMat->VRe[i][j] = bProduct*(thetaRe[i]*phiRe[j+1]*thetaRe[SpecTri->N] +
                                      thetaRe[i]*phiIm[j+1]*thetaIm[SpecTri->N] -
                                      thetaIm[i]*phiIm[j+1]*thetaRe[SpecTri->N] +
                                      thetaIm[i]*phiRe[j+1]*thetaIm[SpecTri->N])/denom;
        invMat->VIm[i][j] = bProduct*(-thetaRe[i]*phiRe[j+1]*thetaIm[SpecTri->N] +
                                      thetaRe[i]*phiIm[j+1]*thetaRe[SpecTri->N] +
                                      thetaIm[i]*phiRe[j+1]*thetaRe[SpecTri->N] +
                                      thetaIm[i]*phiIm[j+1]*thetaIm[SpecTri->N])/denom;
      }
      /* Handle the (-1)^(i+j) part. I'm assuming this way would
       be faster than pow(i+j,-1), but I'm not positive.  */
      if((i+j)%2 != 0)
      {
        invMat->VRe[i][j] *= -1;
        invMat->VIm[i][j] *= -1;
      }
      /* Symmetric part */
      invMat->VRe[j][i] = invMat->VRe[i][j];
      invMat->VIm[j][i] = invMat->VIm[i][j];
    }
  }
  return 0;
}

