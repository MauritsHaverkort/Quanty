/*
 *     ( A[0][0]=acvB  | A[0][0]=acvBA | A[0][1]=bvB0  | 0             | 0           | A[0][4]=bcB0  | 0             | 0           )
 *     ( A[0][0]=acvBA | A[0][0]=acvA  | A[0][1]=bvA0  | 0             | 0           | A[0][4]=bcA0  | 0             | 0           )
 *     ( A[1][0]=bvB0  | A[1][0]=bvA0  | A[1][1]=av2   | A[1][2]=bv2   | 0           | 0             | 0             | 0           )
 *   A=( 0             | 0             | A[2][1]=bv2   | A[2][2]=av3   | A[2][3]=bv3 | 0             | 0             | 0           )
 *     ( 0             | 0             | 0             | A[3][2]=bv3   | A[3][3]=av4 | 0             | 0             | 0           )
 *     ( A[4][0]=bcB0  | A[4][0]=bcA0  | 0             | 0             | 0           | A[4][4]=ac2   | A[4][5]=bc2   | 0           )
 *     ( 0             | 0             | 0             | 0             | 0           | A[5][4]=bc2   | A[5][5]=ac3   | A[5][6]=bc3 )
 *     ( 0             | 0             | 0             | 0             | 0           | 0             | A[6][5]=bc3   | A[6][6]=ac4 )
 *
 *  A = T Ap T^+
 *
 *     ( A[0][0]=av1   | 0             | A[0][1]=bv1   | 0             | 0           | 0             | 0             | 0           )
 *     ( 0             | A[0][0]=ac1   | 0             | 0             | 0           | A[0][4]=bc1   | 0             | 0           )
 *     ( A[1][0]=bv1   | 0             | A[1][1]=av2   | A[1][2]=bv2   | 0           | 0             | 0             | 0           )
 *  Ap=( 0             | 0             | A[2][1]=bv2   | A[2][2]=av3   | A[2][3]=bv3 | 0             | 0             | 0           )
 *     ( 0             | 0             | 0             | A[3][2]=bv3   | A[3][3]=av4 | 0             | 0             | 0           )
 *     ( 0             | A[4][0]=bc1   | 0             | 0             | 0           | A[4][4]=ac2   | A[4][5]=bc2   | 0           )
 *     ( 0             | 0             | 0             | 0             | 0           | A[5][4]=bc2   | A[5][5]=ac3   | A[5][6]=bc3 )
 *     ( 0             | 0             | 0             | 0             | 0           | 0             | A[6][5]=bc3   | A[6][6]=ac4 )
 *
 *
 *     ( A[0][0]=bv0   | A[0][0]=bc0   | 0             | 0             | 0           | 0             | 0             | 0           )
 *     ( A[0][0]=bc0   | A[0][0]=-bv0  | 0             | 0             | 0           | 0             | 0             | 0           )
 *     ( 0             | 0             | 1             | 0             | 0           | 0             | 0             | 0           )
 *   T=( 0             | 0             | 0             | 1             | 0           | 0             | 0             | 0           )
 *     ( 0             | 0             | 0             | 0             | 1           | 0             | 0             | 0           )
 *     ( 0             | 0             | 0             | 0             | 0           | 1             | 0             | 0           )
 *     ( 0             | 0             | 0             | 0             | 0           | 0             | 1             | 0           )
 *     ( 0             | 0             | 0             | 0             | 0           | 0             | 0             | 1           )
 *
 *   (A B)[i][j] = \sum_k A[i][k] . B[k][j]
 *
 *   Note that the block A[2;;4][2;;4] and the block A[5;;7][5;;7] are normal tri-diagonal matrices
 */

/* Free a Bonding-Antibonding DoubleTriDiagonalMatrix M (the structure is not freed)
 */
void FreeBADoubleTriDiagonalMatrix(BADoubleTriDiagonalMatrixType *const M)
{
  FreeTriDiagonalMatrix(&M->val);
  FreeTriDiagonalMatrix(&M->con);
  M->mu=0;
  M->Epsilon=0;
  M->a0=0;
  M->b0=0;
  M->acvB=0;
  M->acvA=0;
  M->acvBA=0;
  M->bvB0=0;
  M->bvA0=0;
  M->bcB0=0;
  M->bcA0=0;
  strcpy(M->Name,"free");
}

/*  Init BA Double Tri diagonal matrix
 *  M.val.N M.con.N must be set
 */
uint32_t InitBADoubleTriDiagonalMatrix(BADoubleTriDiagonalMatrixType *const M)
{
  strcpy(M->Name,"Matrix");
  if((M->val.N==0)||(M->con.N==0)) {printf("Error in InitBADoubleTriDiagonalMatrix M->val.N=%i M->con.N=%i\n",M->val.N,M->con.N);return 2;};
  if(InitTriDiagonalMatrix(&M->val)){printf("InitTriDiagonalMatrix failed in InitBADoubleTriDiagonalMatrix\n");fflush(stdout);return 1;};
  if(InitTriDiagonalMatrix(&M->con))
  {
    FreeTriDiagonalMatrix(&M->val);
    printf("InitTriDiagonalMatrix failed in InitDoubleTriDiagonalMatrix\n");
    fflush(stdout);
    return 1;
  }
  M->a0=0;
  M->b0=1;
  M->acvB=0;
  M->acvA=0;
  M->acvBA=0;
  M->bvB0=0;
  M->bvA0=0;
  M->bcB0=0;
  M->bcA0=0;
  return 0;
}


void SetBADoubleTriDiagonalMatrixConVal(BADoubleTriDiagonalMatrixType *const M)
{
  M->b0=sqrt(Sqr(M->con.b[0])+Sqr(M->val.b[0]));
  M->con.b[0]/=M->b0;
  M->val.b[0]/=M->b0;
  M->acvB = Sqr(M->val.b[0]) * M->val.a[1] + Sqr(M->con.b[0]) * M->con.a[1];
  M->acvA = Sqr(M->val.b[0]) * M->con.a[1] + Sqr(M->con.b[0]) * M->val.a[1];
  M->acvBA= M->val.b[0] * M->con.b[0] * (M->val.a[1] - M->con.a[1]);
  M->bvB0 = M->val.b[0] * M->val.b[1];
  M->bcB0 = M->con.b[0] * M->con.b[1];
  M->bvA0 = M->con.b[0] * M->val.b[1];
  M->bcA0 =-M->val.b[0] * M->con.b[1];
}


/*  Init Identity BA Doulbe Tri diagonal matrix
 *  M.val.N M.con.N must be set
 */
uint32_t InitIdentityBADoubleTriDiagonalMatrix(BADoubleTriDiagonalMatrixType *const M)
{
  strcpy(M->Name,"Matrix");
  if((M->val.N==0)||(M->con.N==0)) {printf("Error in InitBADoubleTriDiagonalMatrix M->val.N=%i M->con.N=%i\n",M->val.N,M->con.N);return 2;};
  if(InitIdentityTriDiagonalMatrix(&M->val)){printf("InitTriDiagonalMatrix failed in InitBADoubleTriDiagonalMatrix\n");fflush(stdout);return 1;};
  if(InitIdentityTriDiagonalMatrix(&M->con))
  {
    FreeTriDiagonalMatrix(&M->val);
    printf("InitTriDiagonalMatrix failed in InitDoubleTriDiagonalMatrix\n");
    fflush(stdout);
    return 1;
  }
  M->a0=0;
  M->b0=1;
  M->acvB=1;
  M->acvA=1;
  M->acvBA=0;
  M->bvB0=0;
  M->bvA0=0;
  M->bcB0=0;
  M->bcA0=0;
  return 0;
}

/* set elements with an absolute value smaller then epsilon to zero
 */
void BADoubleTriDiagonalMatrixChop(BADoubleTriDiagonalMatrixType *const M, double const epsilon)
{
  uint32_t i;
  for(i=0;i<M->val.N;i++)
  {
    if((fabs(M->val.a[i])<epsilon)) M->val.a[i]=0;
    if((fabs(M->val.b[i])<epsilon)) M->val.b[i]=0;
  }
  if((fabs(M->val.a[M->val.N])<epsilon)) M->val.a[M->val.N]=0;
  for(i=0;i<M->con.N;i++)
  {
    if((fabs(M->con.a[i])<epsilon)) M->con.a[i]=0;
    if((fabs(M->con.b[i])<epsilon)) M->con.b[i]=0;
  }
  if((fabs(M->con.a[M->con.N])<epsilon)) M->con.a[M->con.N]=0;
  if((fabs(M->a0   )<epsilon)) M->a0=0;
  if((fabs(M->b0   )<epsilon)) M->b0=0;
  if((fabs(M->acvB )<epsilon)) M->acvB=0;
  if((fabs(M->acvA )<epsilon)) M->acvA=0;
  if((fabs(M->acvBA)<epsilon)) M->acvBA=0;
  if((fabs(M->bvB0 )<epsilon)) M->bvB0=0;
  if((fabs(M->bvA0 )<epsilon)) M->bvA0=0;
  if((fabs(M->bcB0 )<epsilon)) M->bcB0=0;
  if((fabs(M->bcA0 )<epsilon)) M->bcA0=0;
}

/*  MI=M^-1
 *  MI not initiallized on call
 */
double BADoubleTriDiagonalMatrixInverse(BADoubleTriDiagonalMatrixType const M, MatrixType *const MI)
{
  uint32_t i;
  MI->Ni=M.val.N+M.con.N;
  MI->Nj=M.val.N+M.con.N;
  MI->QComplex=0;
  if(InitMatrix(MI)) {printf("InitMatrix failed in BADoubleTriDiagonalMatrixInverse\n");fflush(stdout);return -1;}
  MI->a0=0;
  if(fabs(M.b0)<=DBL_MIN){printf("Warning inverting a matrix with b0=0\n");}
  MI->b0=1/M.b0;
  MI->VRe[0][0]=M.acvB;
  MI->VRe[1][1]=M.acvA;
  MI->VRe[1][0]=M.acvBA;
  MI->VRe[0][1]=M.acvBA;
  MI->VRe[2][0]=M.bvB0;
  MI->VRe[0][2]=M.bvB0;
  MI->VRe[2][1]=M.bvA0;
  MI->VRe[1][2]=M.bvA0;
  MI->VRe[M.val.N+1][0]=M.bcB0;
  MI->VRe[0][M.val.N+1]=M.bcB0;
  MI->VRe[M.val.N+1][1]=M.bcA0;
  MI->VRe[1][M.val.N+1]=M.bcA0;
  
  
  for(i=2;i<M.val.N;i++)
  {
    MI->VRe[i+1][i]=M.val.b[i];
    MI->VRe[i][i+1]=M.val.b[i];
    MI->VRe[i][i]=M.val.a[i];
  }
  MI->VRe[M.val.N][M.val.N]=M.val.a[M.val.N];
  for(i=2;i<M.con.N;i++)
  {
    MI->VRe[M.val.N+i-1][M.val.N+i]=M.con.b[i];
    MI->VRe[M.val.N+i][M.val.N+i-1]=M.con.b[i];
    MI->VRe[M.val.N+i-1][M.val.N+i-1]=M.con.a[i];
  }
  MI->VRe[M.val.N+M.con.N-1][M.val.N+M.con.N-1]=M.con.a[M.con.N];
  return MatrixInverse(MI);
}

/*  M<- M + alpha B
 */
void BADoubleTriDiagonalMatrixAdd(BADoubleTriDiagonalMatrixType * const M, BADoubleTriDiagonalMatrixType B, double alphaRe)
{
  TriDiagonalMatrixAdd(&(M->val), B.val, alphaRe);
  TriDiagonalMatrixAdd(&(M->con), B.con, alphaRe);
  M->acvB += alphaRe*B.acvB;
  M->acvA += alphaRe*B.acvA;
  M->acvBA+= alphaRe*B.acvBA;
  M->bvB0 += alphaRe*B.bvB0;
  M->bvA0 += alphaRe*B.bvA0;
  M->bcB0 += alphaRe*B.bcB0;
  M->bcA0 += alphaRe*B.bcA0;
}

/*  M<- alpha M
 */
void BADoubleTriDiagonalMatrixScalarMultiply(BADoubleTriDiagonalMatrixType * const M, double alphaRe)
{
  TriDiagonalMatrixScalarMultiply(&(M->val), alphaRe);
  TriDiagonalMatrixScalarMultiply(&(M->con), alphaRe);
  M->acvB *= alphaRe;
  M->acvA *= alphaRe;
  M->acvBA*= alphaRe;
  M->bvB0 *= alphaRe;
  M->bvA0 *= alphaRe;
  M->bcB0 *= alphaRe;
  M->bcA0 *= alphaRe;
}

/* Create a new matrix that has dimensions
 * (iEndval,iEndval)+(iEndcon,iEndcon)
 * For the valence and conduction part of the matrix
 * and is a copy of a part of A
 * M not allocated on call of the subroutine
 */
uint32_t BADoubleTriDiagonalMatrixPart(BADoubleTriDiagonalMatrixType * const M, BADoubleTriDiagonalMatrixType A, const uint32_t iEndval, const uint32_t iEndcon)
{
  TriDiagonalMatrixPart(&(M->val), A.val, 0, iEndval);
  TriDiagonalMatrixPart(&(M->con), A.val, 0, iEndcon);
  M->mu=A.mu;
  M->Epsilon=A.Epsilon;
  M->a0=A.a0;
  M->b0=A.b0;
  M->acvB=A.acvB;
  M->acvA=A.acvA;
  M->acvBA=A.acvBA;
  M->bvB0=A.bvB0;
  M->bvA0=A.bvA0;
  M->bcB0=A.bcB0;
  M->bcA0=A.bcA0;
  return 0;
}

/*  M<-A
 */
uint32_t BADoubleTriDiagonalMatrixCopy(BADoubleTriDiagonalMatrixType * const M, BADoubleTriDiagonalMatrixType A)
{
  BADoubleTriDiagonalMatrixPart(M, A, A.val.N, A.val.N);
  return 0;
}







/* ======================= */









/* y <- M x
 * y already allocated
 * y[i]=sum_j M[i][j] x[j]
 */
uint32_t RealBADoubleTriDiagonalMatrixVector(DoubleTriDiagonalMatrixType const M, const double * const x, double * const y)
{
  uint32_t i;
  y[0]=(M.val.a[0]+M.con.a[0])*x[0]/2.0;
  for(i=1;i<=M.val.N;i++)
  {
    y[i]=M.val.a[i]*x[i];
  }
  for(i=M.val.N+1;i<=M.val.N+M.con.N;i++)
  {
    y[i]=M.con.a[i-M.val.N]*x[i];
  }
  for(i=0;i<M.val.N;i++)
  {
    y[i+1]+=M.val.b[i]*x[i];
    y[i]  +=M.val.b[i]*x[i+1];
  }
  y[0]        +=M.con.b[0]*x[M.val.N+1];
  y[M.val.N+1]+=M.con.b[0]*x[0];
  for(i=M.val.N+1;i<M.val.N+M.con.N;i++)
  {
    y[i+1]+=M.con.b[i-M.val.N]*x[i];
    y[i]  +=M.con.b[i-M.val.N]*x[i+1];
  }
  return 0;
}


/* y <- M x
 * y already allocated
 * y[i]=sum_j M[i][j] x[j]
 */
uint32_t ComplexBADoubleTriDiagonalMatrixVector(DoubleTriDiagonalMatrixType const M, const double * const xRe, const double * const xIm, double * const yRe, double * const yIm)
{
  uint32_t i;
  yRe[0]=(M.val.a[0]+M.con.a[0])*xRe[0]/2.0;
  yIm[0]=(M.val.a[0]+M.con.a[0])*xIm[0]/2.0;
  for(i=1;i<=M.val.N;i++)
  {
    yRe[i]=M.val.a[i]*xRe[i];
    yIm[i]=M.val.a[i]*xIm[i];
  }
  for(i=M.val.N+1;i<=M.val.N+M.con.N;i++)
  {
    yRe[i]=M.con.a[i-M.val.N]*xRe[i];
    yIm[i]=M.con.a[i-M.val.N]*xIm[i];
  }
  for(i=0;i<M.val.N;i++)
  {
    yRe[i+1]+=M.val.b[i]*xRe[i];
    yIm[i+1]+=M.val.b[i]*xIm[i];
    yRe[i]+=M.val.b[i]*xRe[i+1];
    yIm[i]+=M.val.b[i]*xIm[i+1];
  }
  yRe[0]+=M.val.b[0]*xRe[M.val.N+1];
  yIm[0]+=M.val.b[0]*xIm[M.val.N+1];
  yRe[M.val.N+1]+=M.val.b[0]*xRe[0];
  yIm[M.val.N+1]+=M.val.b[0]*xIm[0];
  for(i=M.val.N+1;i<M.val.N+M.con.N;i++)
  {
    yRe[i+1]+=M.val.b[i-M.val.N]*xRe[i];
    yIm[i+1]+=M.val.b[i-M.val.N]*xIm[i];
    yRe[i]+=M.val.b[i-M.val.N]*xRe[i+1];
    yIm[i]+=M.val.b[i-M.val.N]*xIm[i+1];
  }
  return 0;
}


/* y <- x M x
 * y already allocated
 * y=sum_ij x[i] M[i][j] x[j]
 */
uint32_t RealVectorBADoubleTriDiagonalMatrixVector(DoubleTriDiagonalMatrixType const M, const double * const x, double *const y)
{
  uint32_t i;
  *y=x[0]*(M.val.a[0]+M.con.a[0])*x[0]/2.0;
  for(i=1;i<=M.val.N;i++)
  {
    *y+=x[i]*M.val.a[i]*x[i];
  }
  for(i=M.val.N+1;i<=M.val.N+M.con.N;i++)
  {
    *y+=x[i]*M.con.a[i-M.val.N]*x[i];
  }
  for(i=0;i<M.val.N;i++)
  {
    *y+=2*x[i]*M.val.b[i]*x[i+1];
  }
  *y+=2*x[M.val.N+1]*M.con.b[0]*x[0];
  for(i=M.val.N+1;i<M.con.N+M.val.N;i++)
  {
    *y+=2*x[i]*M.con.b[i-M.val.N]*x[i+1];
  }
  return 0;
}

/* y <- x^* M x
 * y already allocated
 * y=sum_ij x[i]^* M[i][j] x[j]
 */
uint32_t ComplexVectorBADoubleTriDiagonalMatrixVector(DoubleTriDiagonalMatrixType const M, const double * const xRe, const double * const xIm, double * const yRe, double * const yIm)
{
  uint32_t i;
  *yRe=xRe[0]*(M.val.a[0]+M.con.a[0])*xRe[0]/2.0+xIm[0]*(M.val.a[0]+M.con.a[0])*xIm[0]/2.0;
  *yIm=0;
  for(i=1;i<=M.val.N;i++)
  {
    *yRe+=xRe[i]*M.val.a[i]*xRe[i]+xIm[i]*M.val.a[i]*xIm[i];
  }
  for(i=M.val.N+1;i<M.val.N+M.con.N;i++)
  {
    *yRe+=xRe[i]*M.con.a[i-M.val.N]*xRe[i]+xIm[i]*M.con.a[i-M.val.N]*xIm[i];
  }
  for(i=0;i<M.val.N;i++)
  {
    *yRe+=2*(xRe[i]*M.val.b[i]*xRe[i+1]+xIm[i]*M.val.b[i]*xIm[i+1]);
  }
  *yRe+=2*(xRe[M.val.N+1]*M.val.b[0]*xRe[0]+xIm[M.val.N+1]*M.val.b[0]*xIm[0]);
  for(i=1+M.val.N;i<M.con.N+M.val.N;i++)
  {
    *yRe+=2*(xRe[i]*M.val.b[i-M.val.N]*xRe[i+1]+xIm[i]*M.val.b[i-M.val.N]*xIm[i+1]);
  }
  return 0;
}
