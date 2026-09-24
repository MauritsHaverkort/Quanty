/*
 *   Initialize a struct of BlockTridiagonalMatrix type
 *   malloc the space for storing N+1 matrices of dimension BlockSize * BlockSize storing A0 to AN
 *   malloc the space for storing N   matrices of dimension BlockSize * BlockSize storing B0 to B(N-1)
 *   malloc the space for storing N+1 integers in BlockSizeList
 *   set B0 to the identity matrix
 */
uint32_t InitBlockTridiagonalMatrix(BlockTriDiagonalMatrixType *const A)
{
  int i;
  strcpy(A->Name, "Block Tridiagonal Matrix");
  A->BlockSizeList= (uint32_t*)malloc((A->N+1)*sizeof(uint32_t));
  if(A->BlockSizeList==NULL)
  {
    printf("malloc failed in InitBlockTridiagonalMatrix\n");
    return 1;
  }
  for(i=0;i<=A->N;i++)
  {
    A->BlockSizeList[i]=A->BlockSize;
  }
  if(A->QComplex)
  {
    A->a= (Complex16*)calloc((A->N+1)*(A->BlockSize*A->BlockSize),sizeof(Complex16));
    if(A->a==NULL)
    {
      free(A->BlockSizeList);
      printf("malloc failed in InitBlockTridiagonalMatrix\n");
      return 1;
    }
    A->b= (Complex16*)calloc((A->N)*(A->BlockSize*A->BlockSize),sizeof(Complex16));
    if(A->b==NULL)
    {
      free(A->BlockSizeList);
      free(A->a);
      printf("malloc failed in InitBlockTridiagonalMatrix\n");
      return 1;
    }
    if(A->N!=0) for(i=0;i<A->BlockSize;i++)
    {
      A->b[i+i*A->BlockSize].Re=1;
    }
  }
  else
  {
    A->aRe= (double*)calloc((A->N+1)*(A->BlockSize*A->BlockSize),sizeof(double));
    if(A->aRe==NULL)
    {
      free(A->BlockSizeList);
      printf("malloc failed in InitBlockTridiagonalMatrix\n");
      return 1;
    }
    A->bRe= (double*)calloc((A->N)*(A->BlockSize*A->BlockSize),sizeof(double));
    if(A->bRe==NULL)
    {
      free(A->BlockSizeList);
      free(A->aRe);
      printf("malloc failed in InitBlockTridiagonalMatrix\n");
      return 1;
    }
    if(A->N!=0) for(i=0;i<A->BlockSize;i++)
    {
      A->bRe[i+i*A->BlockSize]=1;
    }
  }
  A->mu=0;
  return 0;
}

uint32_t FreeBlockTridiagonalMatrix(BlockTriDiagonalMatrixType *A)
{
  if(A->QComplex)
  {
    free(A->a);
    free(A->b);
  }
  else
  {
    free(A->aRe);
    free(A->bRe);
  }
  free(A->BlockSizeList);
  return 0;
}

uint32_t ComplexBlockTridiagonalMatrixToReal(BlockTriDiagonalMatrixType *const A)
{
  uint32_t i;
  if(A->QComplex==0) return 0;
  A->QComplex=0;
  A->aRe=(double*)(A->a);
  A->bRe=(double*)(A->b);
  for(i=0;i<A->BlockSize*A->BlockSize*(A->N+1);i++)
  {
    A->aRe[i]=A->a[i].Re;
  }
  for(i=0;i<A->BlockSize*A->BlockSize*(A->N);i++)
  {
    A->bRe[i]=A->b[i].Re;
  }
  A->aRe= (double*)realloc(A->aRe, A->BlockSize*A->BlockSize*(1+A->N)*sizeof(double));
  A->bRe= (double*)realloc(A->bRe, A->BlockSize*A->BlockSize*(A->N)*sizeof(double));
  A->QComplex=0;
  return 0;
}

uint32_t RealBlockTridiagonalMatrixToComplex(BlockTriDiagonalMatrixType *const A)
{
  uint32_t i;
  if(A->QComplex!=0) return 0;
  A->a= (Complex16*)calloc(A->BlockSize*A->BlockSize*(1+A->N),sizeof(Complex16));
  if(A->a==NULL)
  {
    printf("calloc failed in RealBlockTridiagonalMatrixToComplex\n");
    return 1;
  }
  A->b= (Complex16*)calloc(A->BlockSize*A->BlockSize*(A->N),sizeof(Complex16));
  if(A->b==NULL)
  {
    printf("calloc failed in RealBlockTridiagonalMatrixToComplex\n");
    return 1;
  }
  for(i=0;i<A->BlockSize*A->BlockSize*(A->N+1);i++)
  {
    A->a[i].Re=A->aRe[i];
  }
  for(i=0;i<A->BlockSize*A->BlockSize*(A->N);i++)
  {
    A->b[i].Re=A->bRe[i];
  }
  free(A->aRe);
  free(A->bRe);
  A->QComplex=1;
  return 0;
}

uint32_t BlockTridiagonalMatrixChop(BlockTriDiagonalMatrixType *const A, double const Epsilon)
{
  uint32_t i,KeepComplex=0;
  if(A->QComplex)
  {
    for(i=0;i<A->BlockSize*A->BlockSize*(A->N+1);i++)
    {
      if(Abs(A->a[i].Re)<Epsilon)A->a[i].Re=0;
      if(Abs(A->a[i].Im)<Epsilon)
      {
        A->a[i].Im=0;
      }
      else
      {
        KeepComplex=1;
      }
    }
    for(i=0;i<A->BlockSize*A->BlockSize*(A->N);i++)
    {
      if(Abs(A->b[i].Re)<Epsilon)A->b[i].Re=0;
      if(Abs(A->b[i].Im)<Epsilon)
      {
        A->b[i].Im=0;
      }
      else
      {
        KeepComplex=1;
      }
    }
    if(!KeepComplex) ComplexBlockTridiagonalMatrixToReal(A);
  }
  else
  {
    for(i=0;i<A->BlockSize*A->BlockSize*(A->N+1);i++)
    {
      if(Abs(A->aRe[i])<Epsilon)A->aRe[i]=0;
    }
    for(i=0;i<A->BlockSize*A->BlockSize*(A->N);i++)
    {
      if(Abs(A->bRe[i])<Epsilon)A->bRe[i]=0;
    }
  }
  return 0;
}

/*  Cut the chain to at most N, keeping the blocks a[0..N] and b[0..N-1] and dropping
 *  the rest without looking at the values further down the chain. A chain that is
 *  already short enough is left as it is.
 */
uint32_t BlockTridiagonalMatrixTruncate(BlockTriDiagonalMatrixType *const A, uint32_t const N)
{
  uint32_t const NSqr = A->BlockSize*A->BlockSize;
  if(A->N <= N) return 0;
  A->BlockSizeList = (uint32_t*)realloc(A->BlockSizeList, (N+1)*sizeof(uint32_t));
  if(A->BlockSizeList==NULL)
  {
    printf("realloc for A->BlockSizeList failed in BlockTridiagonalMatrixTruncate\n");
    fflush(stdout);
    return 1;
  }
  if(A->QComplex)
  {
    A->a = (Complex16*)realloc(A->a, (N+1)*NSqr*sizeof(Complex16));
    if(A->a==NULL)
    {
      printf("realloc for A->a failed in BlockTridiagonalMatrixTruncate\n");
      fflush(stdout);
      return 1;
    }
    A->b = (Complex16*)realloc(A->b, (N)*NSqr*sizeof(Complex16));
    if((A->b==NULL)&&(N!=0))
    {
      printf("realloc for A->b failed in BlockTridiagonalMatrixTruncate\n");
      fflush(stdout);
      return 1;
    }
  }
  else
  {
    A->aRe = (double*)realloc(A->aRe, (N+1)*NSqr*sizeof(double));
    if(A->aRe==NULL)
    {
      printf("realloc for A->aRe failed in BlockTridiagonalMatrixTruncate\n");
      fflush(stdout);
      return 1;
    }
    A->bRe = (double*)realloc(A->bRe, (N)*NSqr*sizeof(double));
    if((A->bRe==NULL)&&(N!=0))
    {
      printf("realloc for A->bRe failed in BlockTridiagonalMatrixTruncate\n");
      fflush(stdout);
      return 1;
    }
  }
  A->N = N;
  return 0;
}

void BlockTriDiagonalMatrixSetA0to0(BlockTriDiagonalMatrixType * A)
{
  if(A->QComplex==0)
  {
    memset(A->aRe, 0, A->BlockSize*A->BlockSize*sizeof(double));
  }
  else
  {
    memset(A->a, 0, A->BlockSize*A->BlockSize*sizeof(Complex16));
  }
}



uint32_t BlockTridiagonalMatrixChopReduceDimension(BlockTriDiagonalMatrixType *const A, double const Epsilon)
{
  uint32_t i,j;
  bool KeepComplex=false, KeepBlock;
  if(A->QComplex)
  {
    for(i=0;i<A->N;i++)
    {
      KeepBlock=false;
      for(j=0;j<Sqr(A->BlockSize);j++)
      {
        if(Abs(A->b[i*Sqr(A->BlockSize)+j].Re)<Epsilon)
        {
          A->b[i*Sqr(A->BlockSize)+j].Re=0;
        }
        else
        {
          KeepBlock=true;
        }
        if(Abs(A->b[i*Sqr(A->BlockSize)+j].Im)<Epsilon)
        {
          A->b[i*Sqr(A->BlockSize)+j].Im=0;
        }
        else
        {
          KeepComplex=true;
          KeepBlock=true;
        }
        if(Abs(A->a[i*Sqr(A->BlockSize)+j].Re)<Epsilon)A->a[i*Sqr(A->BlockSize)+j].Re=0;
        if(Abs(A->a[i*Sqr(A->BlockSize)+j].Im)<Epsilon)
        {
          A->a[i*Sqr(A->BlockSize)+j].Im=0;
        }
        else
        {
          KeepComplex=true;
        }
      }
      if(KeepBlock==false)
      {
        A->N=i;
      }
    }
    for(j=0;j<Sqr(A->BlockSize);j++)
    {
      if(Abs(A->a[A->N*Sqr(A->BlockSize)+j].Im)<Epsilon)
      {
        A->a[A->N*Sqr(A->BlockSize)+j].Im=0;
      }
      else
      {
        KeepComplex=true;
      }
    }
    if(KeepComplex==false) ComplexBlockTridiagonalMatrixToReal(A);
  }
  else
  {
    for(i=0;i<A->N;i++)
    {
      KeepBlock=false;
      for(j=0;j<Sqr(A->BlockSize);j++)
      {
        if(Abs(A->bRe[i*Sqr(A->BlockSize)+j])<Epsilon)
        {
          A->bRe[i*Sqr(A->BlockSize)+j]=0;
        }
        else
        {
          KeepBlock=true;
        }
        if(Abs(A->aRe[i*Sqr(A->BlockSize)+j])<Epsilon)A->aRe[i*Sqr(A->BlockSize)+j]=0;
      }
      if(KeepBlock==false)
      {
        A->N=i;
      }
    }
    for(j=0;j<Sqr(A->BlockSize);j++)
    {
      if(Abs(A->aRe[A->N*Sqr(A->BlockSize)+j])<Epsilon)A->aRe[A->N*Sqr(A->BlockSize)+j]=0;
    }
  }
  return 0;
}

void ShiftBlockTridiagonalMatrix(BlockTriDiagonalMatrixType *const A, double const delta)
{
  uint32_t i,k;
  if(A->QComplex)
  {
    for(k=1;k<=A->N;k++)
    {
      for(i=k*A->BlockSize*A->BlockSize;i<(k+1)*A->BlockSize*A->BlockSize;i+=A->BlockSize+1)
      {
        A->a[i].Re+=delta;
      }
    }
  }
  else
  {
    for(k=1;k<=A->N;k++)
    {
      for(i=k*A->BlockSize*A->BlockSize;i<(k+1)*A->BlockSize*A->BlockSize;i+=A->BlockSize+1)
      {
        A->aRe[i]+=delta;
      }
    }
  }
  //A->mu+=delta;
}

/*
 *  B <- A
 *  B not initialized on entry
 */
uint32_t CopyBlockTridiagonalMatrix(BlockTriDiagonalMatrixType A, BlockTriDiagonalMatrixType *B)
{
  *B=A;
  if( (B->BlockSizeList = (uint32_t*)malloc((B->N+1)*sizeof(uint32_t))) == NULL )
  {
    printf("malloc failed in CopyBlockTridiagonalMatrix");
    return 1;
  }
  memcpy(B->BlockSizeList, A.BlockSizeList, (B->N+1)*sizeof(uint32_t));
  
  if(B->QComplex==0)
  {
    B->aRe = (double*)malloc((B->N+1)*B->BlockSize*B->BlockSize*sizeof(double));
    if(B->aRe==NULL)
    {
      printf("malloc failed in CopyBlockTridiagonalMatrix");
      free(B->BlockSizeList);
      return 1;
    }
    memcpy(B->aRe, A.aRe, (B->N+1)*B->BlockSize*B->BlockSize*sizeof(double));
    B->bRe = (double*)malloc((B->N)*B->BlockSize*B->BlockSize*sizeof(double));
    if(B->bRe==NULL)
    {
      free(B->aRe);
      free(B->BlockSizeList);
      printf("malloc failed in CopyBlockTridiagonalMatrix");
      return 1;
    }
    memcpy(B->bRe, A.bRe, (B->N)*B->BlockSize*B->BlockSize*sizeof(double));
  }
  else
  {
    B->a = (Complex16*)malloc((B->N+1)*B->BlockSize*B->BlockSize*sizeof(Complex16));
    if(B->a==NULL)
    {
      free(B->BlockSizeList);
      printf("malloc failed in CopyBlockTridiagonalMatrix");
      return 1;
    }
    memcpy(B->a, A.a, (B->N+1)*B->BlockSize*B->BlockSize*sizeof(Complex16));
    B->b = (Complex16*)malloc((B->N)*B->BlockSize*B->BlockSize*sizeof(Complex16));
    if(B->b==NULL)
    {
      free(B->a);
      free(B->BlockSizeList);
      printf("malloc failed in CopyBlockTridiagonalMatrix");
      return 1;
    }
    memcpy(B->b, A.b, (B->N)*B->BlockSize*B->BlockSize*sizeof(Complex16));
  }
  return 0;
}


/* Given G = a0 + b0^* 1/(w - a1 - b1 (1/(w - ...)) b1^+ ) b0^T
 * return
 *  a1 + b1^* (1/(w - ...)) b1^T
 */
int BlockTriDiagonalMatrixEffectiveHamiltonian(BlockTriDiagonalMatrixType const * const G, BlockTriDiagonalMatrixType * const Heff)
{
  if (G->N==0)
  {
    printf("TriDiagonalMatrixEffectiveHamiltonian failed with on response function with length 0");
    Heff->N=0;
    Heff->QComplex=0;
    Heff->BlockSize = G->BlockSize;
    InitBlockTridiagonalMatrix(Heff);
    return 1;
  }
  Heff->N = G->N-1;
  Heff->QComplex = G->QComplex;
  Heff->BlockSize = G->BlockSize;
  InitBlockTridiagonalMatrix(Heff);
  Heff->mu=G->mu;
  memcpy(Heff->BlockSizeList, &(G->BlockSizeList[1]), (Heff->N+1)*sizeof(uint32_t));
  if(Heff->QComplex==0)
  { // real
    memcpy(Heff->aRe, &(G->aRe[Sqr(G->BlockSize)]), (Heff->N+1)*Sqr(G->BlockSize)*sizeof(double));
    memcpy(Heff->bRe, &(G->bRe[Sqr(G->BlockSize)]), (Heff->N  )*Sqr(G->BlockSize)*sizeof(double));
  }
  else
  { // complex
    memcpy(Heff->a, &(G->a[Sqr(G->BlockSize)]), (Heff->N+1)*Sqr(G->BlockSize)*sizeof(Complex16));
    memcpy(Heff->b, &(G->b[Sqr(G->BlockSize)]), (Heff->N  )*Sqr(G->BlockSize)*sizeof(Complex16));
    for(int i=0; i<Sqr(G->BlockSize); i++)
    {
      Heff->b[i].Im*=-1;
    }
  }
  return 0;
}

int BlockTriDiagonalMatrixFromEffectiveHamiltonian(BlockTriDiagonalMatrixType const * const Heff, BlockTriDiagonalMatrixType * const G)
{
 G->N = Heff->N+1;
 G->QComplex = Heff->QComplex;
 G->BlockSize = Heff->BlockSize;
 InitBlockTridiagonalMatrix(G);
 G->mu=Heff->mu;
 memcpy(&(G->BlockSizeList[1]), Heff->BlockSizeList, (Heff->N+1)*sizeof(uint32_t));
 if(G->QComplex==0)
 { // real
   memcpy(&(G->aRe[Sqr(G->BlockSize)]), Heff->aRe, (Heff->N+1)*Sqr(G->BlockSize)*sizeof(double));
   memcpy(&(G->bRe[Sqr(G->BlockSize)]), Heff->bRe, (Heff->N  )*Sqr(G->BlockSize)*sizeof(double));
 }
 else
 { // complex
   memcpy(&(G->a[Sqr(G->BlockSize)]), Heff->a, (Heff->N+1)*Sqr(G->BlockSize)*sizeof(Complex16));
   memcpy(&(G->b[Sqr(G->BlockSize)]), Heff->b, (Heff->N  )*Sqr(G->BlockSize)*sizeof(Complex16));
   for(int i=0; i<Sqr(G->BlockSize); i++)
   {
     G->b[i+Sqr(G->BlockSize)].Im*=-1;
   }
 }
 return 0;
}


int BlockTriDiagonalMatrixGetTail(BlockTriDiagonalMatrixType const * const A, BlockTriDiagonalMatrixType  * const B, uint32_t N)
{
  if(N>A->N)
  {
    printf("Error in BlockTriDiagonalMatrixGetTail, N (%d) larger then A->N (%d)",N,A->N);
    return 1;
  }
  const uint32_t Nshift = A->N - N;
  B->N = N;
  B->BlockSize = A->BlockSize;
  B->mu = A->mu;
  B->QComplex = A->QComplex;
  InitBlockTridiagonalMatrix(B);
  memcpy(B->BlockSizeList, &(A->BlockSizeList[Nshift]), N*sizeof(uint32_t));
  if(A->QComplex==0)
  {
    memcpy(B->aRe, &(A->aRe[Nshift*A->BlockSize*A->BlockSize]), (N+1)*A->BlockSize*A->BlockSize*sizeof(double));
    memcpy(B->bRe, &(A->bRe[Nshift*A->BlockSize*A->BlockSize]), (N)*A->BlockSize*A->BlockSize*sizeof(double));
  }
  else
  {
    memcpy(B->a, &(A->a[Nshift*A->BlockSize*A->BlockSize]), (N+1)*A->BlockSize*A->BlockSize*sizeof(Complex16));
    memcpy(B->b, &(A->b[Nshift*A->BlockSize*A->BlockSize]), (N)*A->BlockSize*A->BlockSize*sizeof(Complex16));
  }
  return 0;
}
int BlockTriDiagonalMatrixGetTailAndAddSpin(BlockTriDiagonalMatrixType const * const A, BlockTriDiagonalMatrixType  * const B, uint32_t N)
{
  int i,j,k,pA,pB;
  if(N>A->N)
  {
    printf("Error in BlockTriDiagonalMatrixGetTail, N (%d) larger then A->N (%d)",N,A->N);
    return 1;
  }
  const uint32_t Nshift = A->N - N;
  B->N = N;
  B->BlockSize = 2*A->BlockSize;
  B->mu = A->mu;
  B->QComplex = A->QComplex;
  InitBlockTridiagonalMatrix(B);
  for(i=0;i<=N;i++)
  {
    B->BlockSizeList[i] = 2*A->BlockSizeList[Nshift+i];
  }
  pA=Nshift*A->BlockSize*A->BlockSize;
  pB=0;
  if(A->QComplex==0)
  {
    for(k=0;k<N;k++)
    {
      for(i=0;i<A->BlockSize;i++)
      {
        for(j=0;j<A->BlockSize;j++)
        {
          B->aRe[pB]=A->aRe[pA];
          B->aRe[pB+B->BlockSize+1]=A->aRe[pA];
          B->bRe[pB]=A->bRe[pA];
          B->bRe[pB+B->BlockSize+1]=A->bRe[pA];
          pB+=2;
          pA++;
        }
        pB+=B->BlockSize;
      }
    }
    k=N;
    {
      for(i=0;i<A->BlockSize;i++)
      {
        for(j=0;j<A->BlockSize;j++)
        {
          B->aRe[pB]=A->aRe[pA];
          B->aRe[pB+B->BlockSize+1]=A->aRe[pA];
          pB+=2;
          pA++;
        }
        pB+=B->BlockSize;
      }
    }
  }
  else
  {
    for(k=0;k<N;k++)
    {
      for(i=0;i<A->BlockSize;i++)
      {
        for(j=0;j<A->BlockSize;j++)
        {
          B->a[pB]=A->a[pA];
          B->a[pB+B->BlockSize+1]=A->a[pA];
          B->b[pB]=A->b[pA];
          B->b[pB+B->BlockSize+1]=A->b[pA];
          pB+=2;
          pA++;
        }
        pB+=B->BlockSize;
      }
    }
    k=N;
    {
      for(i=0;i<A->BlockSize;i++)
      {
        for(j=0;j<A->BlockSize;j++)
        {
          B->a[pB]=A->a[pA];
          B->a[pB+B->BlockSize+1]=A->a[pA];
          pB+=2;
          pA++;
        }
        pB+=B->BlockSize;
      }
    }
    
  }
  return 0;
}


/* =======================================================================================
 *  The sum of two block tri-diagonal response functions, as a block tri-diagonal one
 * =======================================================================================
 *
 *  With c_X = b0_X^T the two functions are G_X = a0_X + c_X^+ [ (w - H_X)^{-1} ]_{11} c_X,
 *  and their sum is
 *
 *    G_A + G_B = a0_A + a0_B + V^+ [ (w - H_A (+) H_B)^{-1} ] V,    V = ( c_A ; c_B ),
 *
 *  the resolvent of the direct sum of the two chains, seen through the two head couplings.
 *  Its chain is the block Lanczos of H_A (+) H_B started on V. It is found here without a
 *  Lanczos, by unitary transformations that leave the starting block alone:
 *
 *  1) The direct sum is written as a banded Hermitian matrix, the two chains interleaved
 *     block by block, with only the directions their BlockSizeList marks as active. A
 *     unitary on the two head blocks puts an orthonormal basis of the range of V in front.
 *  2) HermitianBlockBandReduce takes the band to block tri-diagonal form, with blocks the
 *     size of that range, by Householder reflections and bulge chasing that never touch the
 *     first block. By the implicit Q theorem the result is the block Lanczos chain, up to a
 *     unitary within every block. It costs O(n^2 r B) for n directions, block size r and
 *     block bandwidth B, where a Lanczos that reorthogonalises costs O(n^3), and it is
 *     backward stable: there is no orthogonality to lose.
 *  3) The unitary within every block is then fixed, block after block, by a polar
 *     decomposition, which makes every coupling Hermitian positive semi-definite, b0
 *     included (b0^T = S^{1/2} with S = V^+ V). That is the gauge the block Lanczos uses
 *     elsewhere (CompactMatrixQRDecompositionEigenHermitian). It makes the chain unique, so
 *     a symmetry of the input survives into every block and degenerate states stay
 *     degenerate.
 *
 *  Deflation. A coupling of rank r' < r means the Krylov space grows by r' directions only.
 *  The chain then continues with blocks of size r': those r' directions are put in front of
 *  the rest of the matrix, which is reduced once more with blocks of size r' (step 2), and
 *  BlockSizeList records r'. The memory stays BlockSize wide and is filled with zeros. A
 *  coupling that is not square can not be Hermitian. It is taken as X Sigma, the left
 *  singular vectors scaled by the singular values, each column with its largest element
 *  real and positive (the same for a head V of rank below BlockSize, row by row). A coupling
 *  of rank zero ends the chain. Every deflation lowers the block size, so at most BlockSize
 *  of these extra reductions happen.
 *
 *  A singular value of a coupling counts as zero below Max(SingularValue, 1000 DBL_EPSILON
 *  |H|), with |H| a bound on the norm of the direct sum, and one of V below Max(SingularValue,
 *  1000 DBL_EPSILON sigma_max(V)). Where the Krylov space is exhausted rounding leaves a
 *  coupling of order DBL_EPSILON |H|, and that has to end the chain.
 *
 *  mu is taken from A, as the list of poles route does.
 */

/*  A Hermitian matrix of NBlock blocks of size r (the last one padded with zeros), of which
 *  the lower block band is stored. Block (I,J) with I>=J and I-J<NDiag is found at
 *  L[(J*NDiag+I-J)*r*r], column major, element (i,j) at [i+j*r]. NDiag is the block
 *  bandwidth plus two, one for the diagonal and one for the bulge of the reduction. Of a
 *  diagonal block only the lower triangle is read.
 */
typedef struct HermitianBlockBandType
{
  uint32_t r;
  uint32_t NBlock;
  uint32_t NDiag;
  Complex16 *L;
}
HermitianBlockBandType;

static int InitHermitianBlockBand(HermitianBlockBandType * const H, uint32_t const n, uint32_t const r, uint32_t const B)
{
  H->r = r;
  H->NBlock = (n+r-1)/r;
  H->NDiag = B+2;
  H->L = (Complex16*)calloc((size_t)(H->NBlock)*(H->NDiag)*r*r, sizeof(Complex16));
  if(H->L==NULL)
  {
    printf("calloc failed in InitHermitianBlockBand\n");
    return 1;
  }
  return 0;
}

static Complex16 * HermitianBlockBandBlock(HermitianBlockBandType const * const H, uint32_t const I, uint32_t const J)
{
  return &(H->L[((size_t)J*(H->NDiag)+(I-J))*(H->r)*(H->r)]);
}

/*  sets element (i,j) of the Hermitian matrix to z, which is element (j,i) set to z^* */
static void HermitianBlockBandSet(HermitianBlockBandType * const H, uint32_t i, uint32_t j, Complex16 z)
{
  uint32_t t;
  if(i<j)
  {
    t=i; i=j; j=t;
    z.Im=-z.Im;
  }
  HermitianBlockBandBlock(H, i/(H->r), j/(H->r))[(i%(H->r))+(j%(H->r))*(H->r)] = z;
}

/*  D = diagonal block I, both triangles, column major */
static void HermitianBlockBandGetDiagonal(HermitianBlockBandType const * const H, uint32_t const I, Complex16 * const D)
{
  uint32_t p,q;
  uint32_t const r=H->r;
  Complex16 const * const blk = HermitianBlockBandBlock(H,I,I);
  for(q=0;q<r;q++)
  {
    for(p=0;p<r;p++)
    {
      if(p>=q)
      {
        D[p+q*r] = blk[p+q*r];
      }
      else
      {
        D[p+q*r].Re =  blk[q+p*r].Re;
        D[p+q*r].Im = -blk[q+p*r].Im;
      }
    }
  }
}


/*  Reduces a Hermitian matrix of block bandwidth NDiag-2 to block tri-diagonal form by a
 *  unitary similarity that leaves block 0 alone.
 *
 *  Sweep Bcur takes the block bandwidth from Bcur to Bcur-1. In block column p the block
 *  (p+Bcur,p) is annihilated by a QR of the two blocks (p+Bcur-1,p) and (p+Bcur,p). The
 *  reflection, applied from both sides to block rows and columns p+Bcur-1 and p+Bcur, fills
 *  in the single block (p+2Bcur,p+Bcur-1). That bulge is annihilated in the same way against
 *  the block above it, which moves it Bcur blocks down, until it falls off the end. The
 *  first rows rotated are p+Bcur-1 >= 1, so block 0 is never touched.
 *
 *  Every step works on a dense copy of the window of blocks j ... i+1+Bcur, which holds all
 *  non-zeros of the two block rows it rotates, and writes the lower band back. A chase has
 *  j = i-Bcur, so the window is at most 2Bcur+2 blocks wide.
 */
static int HermitianBlockBandReduce(HermitianBlockBandType * const H)
{
  uint32_t const r = H->r;
  uint32_t const NBlock = H->NBlock;
  uint32_t const B = H->NDiag-2;
  uint32_t const ldmax = (2*B+2)*r;
  uint32_t r1 = r;
  uint32_t r2 = 2*r;
  uint32_t lwork = 64*r2;
  uint32_t Bcur,p,i,j,lo,hi,ld,I,J,a,b,I0,J0;
  int info=0;
  Complex16 one={1,0}, zero={0,0};
  Complex16 *D, *Q, *T, *tau, *work, *blk;
  if(B<2) return 0;
  D    = (Complex16*)malloc((size_t)ldmax*ldmax*sizeof(Complex16));
  Q    = (Complex16*)malloc((size_t)r2*r2*sizeof(Complex16));
  T    = (Complex16*)malloc((size_t)r2*ldmax*sizeof(Complex16));
  tau  = (Complex16*)malloc((size_t)r2*sizeof(Complex16));
  work = (Complex16*)malloc((size_t)lwork*sizeof(Complex16));
  if((D==NULL)||(Q==NULL)||(T==NULL)||(tau==NULL)||(work==NULL))
  {
    printf("malloc failed in HermitianBlockBandReduce\n");
    free(D); free(Q); free(T); free(tau); free(work);
    return 1;
  }
  for(Bcur=B;Bcur>=2;Bcur--)
  {
    for(p=0;p+Bcur<NBlock;p++)
    {
      for(j=p,i=p+Bcur-1;i+1<NBlock;j=i,i+=Bcur)
      {
        lo = j;
        hi = Min(NBlock-1,i+1+Bcur);
        ld = (hi-lo+1)*r;
        /*  the window, lower band first and then the upper triangle from it */
        for(J=lo;J<=hi;J++)
        {
          for(I=J;I<=hi;I++)
          {
            blk = (I-J<H->NDiag) ? HermitianBlockBandBlock(H,I,J) : NULL;
            for(b=0;b<r;b++)
            {
              for(a=0;a<r;a++)
              {
                D[(I-lo)*r+a+((J-lo)*r+b)*ld] = (blk!=NULL) ? blk[a+b*r] : zero;
              }
            }
          }
        }
        for(b=0;b<ld;b++)
        {
          for(a=0;a<b;a++)
          {
            D[a+b*ld].Re =  D[b+a*ld].Re;
            D[a+b*ld].Im = -D[b+a*ld].Im;
          }
        }
        /*  Q from the QR of the panel, block rows i and i+1 of block column j */
        I0 = (i-lo)*r;
        J0 = (j-lo)*r;
        for(b=0;b<r;b++)
        {
          for(a=0;a<r2;a++)
          {
            Q[a+b*r2] = D[I0+a+(J0+b)*ld];
          }
        }
        zgeqrf(&r2,&r1,Q,&r2,tau,work,&lwork,&info);
        if(info!=0)
        {
          printf("zgeqrf failed in HermitianBlockBandReduce with INFO %i\n",info);
          free(D); free(Q); free(T); free(tau); free(work);
          return 1;
        }
        zungqr(&r2,&r2,&r1,Q,&r2,tau,work,&lwork,&info);
        if(info!=0)
        {
          printf("zungqr failed in HermitianBlockBandReduce with INFO %i\n",info);
          free(D); free(Q); free(T); free(tau); free(work);
          return 1;
        }
        /*  rows I0 ... I0+2r-1 <- Q^+ rows, then columns I0 ... I0+2r-1 <- columns Q */
        zgemm("C","N",&r2,&ld,&r2,&one,Q,&r2,&(D[I0]),&ld,&zero,T,&r2);
        for(b=0;b<ld;b++)
        {
          for(a=0;a<r2;a++)
          {
            D[I0+a+b*ld] = T[a+b*r2];
          }
        }
        zgemm("N","N",&ld,&r2,&r2,&one,&(D[(size_t)I0*ld]),&ld,Q,&r2,&zero,T,&ld);
        memcpy(&(D[(size_t)I0*ld]),T,(size_t)ld*r2*sizeof(Complex16));
        /*  the lower band back, the annihilated block as an exact zero */
        for(J=lo;J<=hi;J++)
        {
          for(I=J;(I<=hi)&&(I-J<H->NDiag);I++)
          {
            blk = HermitianBlockBandBlock(H,I,J);
            for(b=0;b<r;b++)
            {
              for(a=0;a<r;a++)
              {
                blk[a+b*r] = ((I==i+1)&&(J==j)) ? zero : D[(I-lo)*r+a+((J-lo)*r+b)*ld];
              }
            }
          }
        }
      }
    }
  }
  free(D); free(Q); free(T); free(tau); free(work);
  return 0;
}

/*  element (p,q) of block k of a (ab='a') or b (ab='b') of a block tri-diagonal matrix */
static Complex16 BlockTriDiagonalMatrixElement(BlockTriDiagonalMatrixType const * const X, char const ab, uint32_t const k, uint32_t const p, uint32_t const q)
{
  Complex16 z;
  size_t const idx = ((size_t)k*(X->BlockSize)+p)*(X->BlockSize)+q;
  if(X->QComplex)
  {
    return (ab=='a') ? X->a[idx] : X->b[idx];
  }
  z.Re = (ab=='a') ? X->aRe[idx] : X->bRe[idx];
  z.Im = 0;
  return z;
}

/*  the singular value decomposition M = X diag(sv) VT of an m x n column major matrix, with
 *  X m x m and VT n x n. M is overwritten. */
static int ColumnMajorSVD(uint32_t m, uint32_t n, Complex16 * const M, double * const sv, Complex16 * const X, Complex16 * const VT)
{
  uint32_t lwork = 8*(m+n)+64;
  int info=0;
  Complex16 *work;
  double *rwork;
  work  = (Complex16*)malloc((size_t)lwork*sizeof(Complex16));
  rwork = (double*)malloc((size_t)(5*Min(m,n)+5)*sizeof(double));
  if((work==NULL)||(rwork==NULL))
  {
    printf("malloc failed in ColumnMajorSVD\n");
    free(work); free(rwork);
    return 1;
  }
  zgesvd("A","A",&m,&n,M,&m,sv,X,&m,VT,&n,work,&lwork,rwork,&info);
  free(work); free(rwork);
  if(info!=0)
  {
    printf("zgesvd failed in ColumnMajorSVD with INFO %i\n",info);
    return 1;
  }
  return 0;
}

/*  out = (M + M^+)/2 of an r x r column major M, row major into an m x m block */
static void StoreHermitianBlock(Complex16 const * const M, uint32_t const r, Complex16 * const out, uint32_t const m)
{
  uint32_t p,q;
  for(p=0;p<r;p++)
  {
    for(q=0;q<r;q++)
    {
      out[p*m+q].Re = 0.5*(M[p+q*r].Re + M[q+p*r].Re);
      out[p*m+q].Im = 0.5*(M[p+q*r].Im - M[q+p*r].Im);
    }
  }
}

static int GrowChainBlocks(Complex16 ** const a, Complex16 ** const b, uint32_t ** const s, uint32_t * const cap, uint32_t const m)
{
  void *ta, *tb, *ts;
  uint32_t const newcap = 2*(*cap);
  ta = realloc(*a, (size_t)newcap*m*m*sizeof(Complex16));
  if(ta==NULL) {printf("realloc failed in GrowChainBlocks\n"); return 1;}
  *a = (Complex16*)ta;
  tb = realloc(*b, (size_t)newcap*m*m*sizeof(Complex16));
  if(tb==NULL) {printf("realloc failed in GrowChainBlocks\n"); return 1;}
  *b = (Complex16*)tb;
  ts = realloc(*s, (size_t)newcap*sizeof(uint32_t));
  if(ts==NULL) {printf("realloc failed in GrowChainBlocks\n"); return 1;}
  *s = (uint32_t*)ts;
  memset(&((*a)[(size_t)(*cap)*m*m]), 0, (size_t)(newcap-(*cap))*m*m*sizeof(Complex16));
  memset(&((*b)[(size_t)(*cap)*m*m]), 0, (size_t)(newcap-(*cap))*m*m*sizeof(Complex16));
  *cap = newcap;
  return 0;
}

/*  C = A + B, see the comment above HermitianBlockBandType. C is not allocated on entry. */
int BlockTriDiagonalMatrixAdd(BlockTriDiagonalMatrixType const * const A, BlockTriDiagonalMatrixType const * const B, BlockTriDiagonalMatrixType * const C, double const SingularValue)
{
  BlockTriDiagonalMatrixType const * const X[2] = {A, B};
  uint32_t const m = A->BlockSize;
  uint32_t const Kmax = Max(A->N, B->N);
  uint32_t x,k,p,q,g,nh,n,r,rp,w,t,u,K=0,cap=16,ntail,Lb;
  uint32_t s1[2], off1[2];
  uint32_t *s=NULL, *off=NULL, *sOut=NULL;
  double Hnorm=0, fa, fb, fbprev, ZeroH, ZeroV, svmax;
  double *sv=NULL;
  Complex16 one={1,0}, zero={0,0}, z, phase;
  Complex16 *V=NULL, *Xs=NULL, *VT=NULL, *Qh=NULL, *Hh=NULL, *T1=NULL, *T2=NULL, *Cp=NULL;
  Complex16 *Dm=NULL, *U=NULL, *Un=NULL, *M=NULL, *S=NULL, *b0=NULL;
  Complex16 *aOut=NULL, *bOut=NULL;
  HermitianBlockBandType H, Hn;
  int retval=1;
  H.L=NULL;
  if(A->BlockSize!=B->BlockSize)
  {
    printf("Error in BlockTriDiagonalMatrixAdd: block sizes %u and %u differ\n",A->BlockSize,B->BlockSize);
    return 1;
  }

  /*  the active size of block k of each chain, 0 behind its end, and where it sits in the
   *  direct sum: the two first blocks together form the head, the others follow as A_2,
   *  B_2, A_3, B_3, ... */
  s   = (uint32_t*)calloc(2*((size_t)Kmax+2),sizeof(uint32_t));
  off = (uint32_t*)calloc(2*((size_t)Kmax+2),sizeof(uint32_t));
  if((s==NULL)||(off==NULL)) {printf("calloc failed in BlockTriDiagonalMatrixAdd\n"); goto Cleanup;}
  for(x=0;x<2;x++)
  {
    for(k=1;k<=X[x]->N;k++) s[2*k+x] = X[x]->BlockSizeList[k];
  }
  s1[0] = s[2]; s1[1] = s[3];
  off1[0] = 0;  off1[1] = s1[0];
  nh = s1[0]+s1[1];
  n = nh;
  for(k=2;k<=Kmax;k++)
  {
    for(x=0;x<2;x++)
    {
      off[2*k+x] = n;
      n += s[2*k+x];
    }
  }
  /*  |H| <= max over the blocks of |a_k| + |b_k| + |b_{k-1}|, Frobenius norms */
  for(x=0;x<2;x++)
  {
    fbprev=0;
    for(k=1;k<=X[x]->N;k++)
    {
      fa=0;
      for(g=0;g<m*m;g++) {z=BlockTriDiagonalMatrixElement(X[x],'a',k,g/m,g%m); fa+=Sqr(z.Re)+Sqr(z.Im);}
      fb=0;
      if(k<X[x]->N) for(g=0;g<m*m;g++) {z=BlockTriDiagonalMatrixElement(X[x],'b',k,g/m,g%m); fb+=Sqr(z.Re)+Sqr(z.Im);}
      Hnorm = Max(Hnorm, sqrt(fa)+sqrt(fb)+fbprev);
      fbprev = sqrt(fb);
    }
  }

  sOut = (uint32_t*)malloc(cap*sizeof(uint32_t));
  aOut = (Complex16*)calloc((size_t)cap*m*m,sizeof(Complex16));
  bOut = (Complex16*)calloc((size_t)cap*m*m,sizeof(Complex16));
  b0   = (Complex16*)calloc((size_t)m*m,sizeof(Complex16));
  sv   = (double*)malloc((size_t)(2*m+1)*sizeof(double));
  Dm = (Complex16*)malloc((size_t)m*m*sizeof(Complex16));
  U  = (Complex16*)malloc((size_t)m*m*sizeof(Complex16));
  Un = (Complex16*)malloc((size_t)m*m*sizeof(Complex16));
  M  = (Complex16*)malloc((size_t)m*m*sizeof(Complex16));
  S  = (Complex16*)malloc((size_t)m*m*sizeof(Complex16));
  if((sOut==NULL)||(aOut==NULL)||(bOut==NULL)||(b0==NULL)||(sv==NULL)||(Dm==NULL)||(U==NULL)||(Un==NULL)||(M==NULL)||(S==NULL))
  {
    printf("malloc failed in BlockTriDiagonalMatrixAdd\n");
    goto Cleanup;
  }
  if(nh==0) goto Assemble;

  /*  V = ( c_A ; c_B ) on the head, c_X = b0_X^T, and its singular value decomposition */
  V  = (Complex16*)calloc((size_t)nh*m,sizeof(Complex16));
  Xs = (Complex16*)malloc((size_t)nh*nh*sizeof(Complex16));
  VT = (Complex16*)malloc((size_t)m*m*sizeof(Complex16));
  Qh = (Complex16*)malloc((size_t)nh*nh*sizeof(Complex16));
  Hh = (Complex16*)calloc((size_t)nh*nh,sizeof(Complex16));
  T1 = (Complex16*)malloc((size_t)nh*Max(nh,m)*sizeof(Complex16));
  Cp = (Complex16*)malloc((size_t)nh*Max(nh,m)*sizeof(Complex16));
  if((V==NULL)||(Xs==NULL)||(VT==NULL)||(Qh==NULL)||(Hh==NULL)||(T1==NULL)||(Cp==NULL))
  {
    printf("malloc failed in BlockTriDiagonalMatrixAdd\n");
    goto Cleanup;
  }
  for(x=0;x<2;x++)
  {
    for(p=0;p<s1[x];p++)
    {
      for(q=0;q<m;q++) V[off1[x]+p+q*nh] = BlockTriDiagonalMatrixElement(X[x],'b',0,q,p);
    }
  }
  memcpy(Cp,V,(size_t)nh*m*sizeof(Complex16));
  if(ColumnMajorSVD(nh,m,Cp,sv,Xs,VT)) goto Cleanup;
  svmax = sv[0];
  ZeroV = Max(SingularValue, 1000*DBL_EPSILON*svmax);
  ZeroH = Max(SingularValue, 1000*DBL_EPSILON*Hnorm);
  r=0;
  while((r<Min(nh,m))&&(sv[r]>ZeroV)) r++;
  if(r==0) goto Assemble;

  /*  Qh = ( E , E_perp ), E an orthonormal basis of the range of V. For V of full rank
   *  E = X_m VT, the polar factor, and E^+ V = S^{1/2}. Otherwise E = X_r with phases that
   *  make the largest element of every row of E^+ V real and positive. */
  memcpy(Qh,Xs,(size_t)nh*nh*sizeof(Complex16));
  if(r==m)
  {
    zgemm("N","N",&nh,&m,&m,&one,Xs,&nh,VT,&m,&zero,Qh,&nh);
  }
  else
  {
    for(q=0;q<r;q++)
    {
      g=0;
      for(p=1;p<m;p++) if(Sqr(VT[q+p*m].Re)+Sqr(VT[q+p*m].Im) > Sqr(VT[q+g*m].Re)+Sqr(VT[q+g*m].Im)) g=p;
      fa = sqrt(Sqr(VT[q+g*m].Re)+Sqr(VT[q+g*m].Im));
      phase.Re = VT[q+g*m].Re/fa; phase.Im = VT[q+g*m].Im/fa;
      for(p=0;p<nh;p++)
      {
        z = Qh[p+q*nh];
        Qh[p+q*nh].Re = z.Re*phase.Re - z.Im*phase.Im;
        Qh[p+q*nh].Im = z.Re*phase.Im + z.Im*phase.Re;
      }
    }
  }
  /*  b0 = (E^+ V)^T, Hermitian for a V of full rank */
  zgemm("C","N",&r,&m,&nh,&one,Qh,&nh,V,&nh,&zero,Cp,&r);
  if(r==m)
  {
    StoreHermitianBlock(Cp,m,T1,m);
    for(p=0;p<m;p++) for(q=0;q<m;q++) b0[q*m+p] = T1[p*m+q];
  }
  else
  {
    for(p=0;p<r;p++) for(q=0;q<m;q++) b0[q*m+p] = Cp[p+q*r];
  }

  /*  the band: head rotated by Qh, then A_2, B_2, A_3, ... Its half bandwidth w is the
   *  largest distance between two coupled directions. */
  w = nh-1;
  if(Kmax>=2) w = Max(w, off[2*2+1]+s[2*2+1]-1);
  for(k=2;k<Kmax;k++)
  {
    for(x=0;x<2;x++)
    {
      if((s[2*k+x]!=0)&&(s[2*(k+1)+x]!=0)) w = Max(w, off[2*(k+1)+x]+s[2*(k+1)+x]-1-off[2*k+x]);
    }
  }
  if(InitHermitianBlockBand(&H, n, r, Max(1,(w+r-1)/r))) goto Cleanup;
  for(x=0;x<2;x++)
  {
    for(p=0;p<s1[x];p++)
    {
      for(q=0;q<s1[x];q++) Hh[off1[x]+p+(off1[x]+q)*nh] = BlockTriDiagonalMatrixElement(X[x],'a',1,p,q);
    }
  }
  zgemm("N","N",&nh,&nh,&nh,&one,Hh,&nh,Qh,&nh,&zero,T1,&nh);
  zgemm("C","N",&nh,&nh,&nh,&one,Qh,&nh,T1,&nh,&zero,Hh,&nh);
  for(q=0;q<nh;q++) for(p=q;p<nh;p++) HermitianBlockBandSet(&H,p,q,Hh[p+q*nh]);
  for(x=0;x<2;x++)
  {
    if(s[2*2+x]==0) continue;
    /*  Cp = Qh^+ ( b1_X on the rows of X_1 ), the coupling of the head to X_2 */
    for(p=0;p<nh*s[2*2+x];p++) T1[p]=zero;
    for(p=0;p<s1[x];p++)
    {
      for(q=0;q<s[2*2+x];q++) T1[off1[x]+p+q*nh] = BlockTriDiagonalMatrixElement(X[x],'b',1,p,q);
    }
    zgemm("C","N",&nh,&(s[2*2+x]),&nh,&one,Qh,&nh,T1,&nh,&zero,Cp,&nh);
    for(p=0;p<nh;p++)
    {
      for(q=0;q<s[2*2+x];q++) HermitianBlockBandSet(&H,p,off[2*2+x]+q,Cp[p+q*nh]);
    }
  }
  for(k=2;k<=Kmax;k++)
  {
    for(x=0;x<2;x++)
    {
      for(p=0;p<s[2*k+x];p++)
      {
        for(q=0;q<=p;q++) HermitianBlockBandSet(&H,off[2*k+x]+p,off[2*k+x]+q,BlockTriDiagonalMatrixElement(X[x],'a',k,p,q));
        if(k<X[x]->N)
        {
          for(q=0;q<s[2*(k+1)+x];q++) HermitianBlockBandSet(&H,off[2*k+x]+p,off[2*(k+1)+x]+q,BlockTriDiagonalMatrixElement(X[x],'b',k,p,q));
        }
      }
    }
  }
  if(HermitianBlockBandReduce(&H)) goto Cleanup;

  /*  the walk along the chain. U is the basis of the current chain block in the coordinates
   *  of block t of the band. */
  for(p=0;p<r*r;p++) U[p]=zero;
  for(p=0;p<r;p++) U[p+p*r]=one;
  t=0;
  for(;;)
  {
    if((K>=cap)&&GrowChainBlocks(&aOut,&bOut,&sOut,&cap,m)) goto Cleanup;
    /*  a = U^+ D_t U */
    HermitianBlockBandGetDiagonal(&H,t,Dm);
    zgemm("N","N",&r,&r,&r,&one,Dm,&r,U,&r,&zero,S,&r);
    zgemm("C","N",&r,&r,&r,&one,U,&r,S,&r,&zero,M,&r);
    StoreHermitianBlock(M,r,&(aOut[(size_t)K*m*m]),m);
    sOut[K]=r;
    K++;
    if(t+1>=H.NBlock) break;
    /*  the coupling to the next block, U^+ H(t,t+1) = U^+ H(t+1,t)^+ = X diag(sv) VT */
    zgemm("C","C",&r,&r,&r,&one,U,&r,HermitianBlockBandBlock(&H,t+1,t),&r,&zero,M,&r);
    if(ColumnMajorSVD(r,r,M,sv,Dm,S)) goto Cleanup;
    rp=0;
    while((rp<r)&&(sv[rp]>ZeroH)) rp++;
    if(rp==0) break;
    if(rp==r)
    {
      /*  the next basis U VT^+ X^+ makes the coupling X diag(sv) X^+ */
      zgemm("C","C",&r,&r,&r,&one,S,&r,Dm,&r,&zero,Un,&r);
      for(q=0;q<r;q++) for(p=0;p<r;p++) {M[p+q*r].Re = Dm[p+q*r].Re*sv[q]; M[p+q*r].Im = Dm[p+q*r].Im*sv[q];}
      zgemm("N","C",&r,&r,&r,&one,M,&r,Dm,&r,&zero,U,&r);
      StoreHermitianBlock(U,r,&(bOut[(size_t)(K-1)*m*m]),m);
      memcpy(U,Un,(size_t)r*r*sizeof(Complex16));
      t++;
      continue;
    }
    /*  deflation: the coupling X_rp diag(sv), each column with its largest element real
     *  and positive, and the rest of the band again from the rp directions VT^+_rp */
    for(q=0;q<rp;q++)
    {
      g=0;
      for(p=1;p<r;p++) if(Sqr(Dm[p+q*r].Re)+Sqr(Dm[p+q*r].Im) > Sqr(Dm[g+q*r].Re)+Sqr(Dm[g+q*r].Im)) g=p;
      fa = sqrt(Sqr(Dm[g+q*r].Re)+Sqr(Dm[g+q*r].Im));
      phase.Re = Dm[g+q*r].Re/fa; phase.Im = Dm[g+q*r].Im/fa;
      for(p=0;p<r;p++)
      {
        z = Dm[p+q*r];
        Dm[p+q*r].Re = z.Re*phase.Re + z.Im*phase.Im;
        Dm[p+q*r].Im = z.Im*phase.Re - z.Re*phase.Im;
        z = S[q+p*r];
        S[q+p*r].Re = z.Re*phase.Re - z.Im*phase.Im;
        S[q+p*r].Im = z.Re*phase.Im + z.Im*phase.Re;
      }
      for(p=0;p<r;p++)
      {
        bOut[(size_t)(K-1)*m*m+p*m+q].Re = Dm[p+q*r].Re*sv[q];
        bOut[(size_t)(K-1)*m*m+p*m+q].Im = Dm[p+q*r].Im*sv[q];
      }
    }
    /*  Un = VT^+, the new basis of block t+1 in its band coordinates, the rp active
     *  directions first. The right singular vectors of U^+ H(t,t+1) live there. */
    for(q=0;q<r;q++) for(p=0;p<r;p++) {Un[p+q*r].Re = S[q+p*r].Re; Un[p+q*r].Im = -S[q+p*r].Im;}
    ntail = (H.NBlock-t-1)*r;
    if(InitHermitianBlockBand(&Hn, ntail, rp, (2*r-1+rp-1)/rp)) goto Cleanup;
    for(u=t+1;u<H.NBlock;u++)
    {
      Lb = (u-t-1)*r;
      HermitianBlockBandGetDiagonal(&H,u,Dm);
      if(u==t+1)
      {
        zgemm("N","N",&r,&r,&r,&one,Dm,&r,Un,&r,&zero,S,&r);
        zgemm("C","N",&r,&r,&r,&one,Un,&r,S,&r,&zero,Dm,&r);
      }
      for(q=0;q<r;q++) for(p=q;p<r;p++) HermitianBlockBandSet(&Hn,Lb+p,Lb+q,Dm[p+q*r]);
      if(u+1<H.NBlock)
      {
        memcpy(S,HermitianBlockBandBlock(&H,u+1,u),(size_t)r*r*sizeof(Complex16));
        if(u==t+1)
        {
          zgemm("N","N",&r,&r,&r,&one,S,&r,Un,&r,&zero,M,&r);
          memcpy(S,M,(size_t)r*r*sizeof(Complex16));
        }
        for(q=0;q<r;q++) for(p=0;p<r;p++) HermitianBlockBandSet(&Hn,Lb+r+p,Lb+q,S[p+q*r]);
      }
    }
    free(H.L);
    H = Hn;
    if(HermitianBlockBandReduce(&H)) goto Cleanup;
    r = rp;
    for(p=0;p<r*r;p++) U[p]=zero;
    for(p=0;p<r;p++) U[p+p*r]=one;
    t=0;
  }

Assemble:
  C->N = K;
  C->BlockSize = m;
  C->QComplex = 1;
  if(InitBlockTridiagonalMatrix(C)) goto Cleanup;
  C->mu = A->mu;
  C->BlockSizeList[0] = m;
  for(p=0;p<m*m;p++)
  {
    C->a[p] = BlockTriDiagonalMatrixElement(A,'a',0,p/m,p%m);
    z = BlockTriDiagonalMatrixElement(B,'a',0,p/m,p%m);
    C->a[p].Re += z.Re;
    C->a[p].Im += z.Im;
  }
  if(K>0)
  {
    memcpy(&(C->a[(size_t)m*m]),aOut,(size_t)K*m*m*sizeof(Complex16));
    memcpy(C->b,b0,(size_t)m*m*sizeof(Complex16));
    if(K>1) memcpy(&(C->b[(size_t)m*m]),bOut,(size_t)(K-1)*m*m*sizeof(Complex16));
    memcpy(&(C->BlockSizeList[1]),sOut,(size_t)K*sizeof(uint32_t));
  }
  /*  a sum of two real functions is real, and every step above keeps real data real */
  if((A->QComplex==0)&&(B->QComplex==0))
  {
    g=0;
    for(p=0;p<(K+1)*m*m;p++) if(C->a[p].Im!=0) g=1;
    for(p=0;p<K*m*m;p++)     if(C->b[p].Im!=0) g=1;
    if(g==0) ComplexBlockTridiagonalMatrixToReal(C);
  }
  retval=0;
Cleanup:
  free(s); free(off); free(sOut); free(aOut); free(bOut); free(b0); free(sv);
  free(V); free(Xs); free(VT); free(Qh); free(Hh); free(T1); free(Cp);
  free(Dm); free(U); free(Un); free(M); free(S);
  free(H.L);
  if(retval) printf("Error in BlockTriDiagonalMatrixAdd\n");
  return retval;
}

/*  C = A + B for single valued chains, as block chains of block size 1. C is not allocated
 *  on entry. */
int TriDiagonalMatrixAdd(TriDiagonalMatrixType const * const A, TriDiagonalMatrixType const * const B, TriDiagonalMatrixType * const C, double const SingularValue)
{
  TriDiagonalMatrixType const * const X[2] = {A, B};
  BlockTriDiagonalMatrixType BX[2], BC;
  uint32_t x,k;
  for(x=0;x<2;x++)
  {
    BX[x].N = X[x]->N;
    BX[x].BlockSize = 1;
    BX[x].QComplex = 0;
    if(InitBlockTridiagonalMatrix(&(BX[x])))
    {
      if(x==1) FreeBlockTridiagonalMatrix(&(BX[0]));
      printf("InitBlockTridiagonalMatrix failed in TriDiagonalMatrixAdd\n");
      return 1;
    }
    for(k=0;k<=X[x]->N;k++) BX[x].aRe[k] = X[x]->a[k];
    for(k=0;k<X[x]->N;k++)  BX[x].bRe[k] = X[x]->b[k];
    BX[x].mu = X[x]->mu;
  }
  if(BlockTriDiagonalMatrixAdd(&(BX[0]),&(BX[1]),&BC,SingularValue))
  {
    FreeBlockTridiagonalMatrix(&(BX[0]));
    FreeBlockTridiagonalMatrix(&(BX[1]));
    printf("BlockTriDiagonalMatrixAdd failed in TriDiagonalMatrixAdd\n");
    return 1;
  }
  FreeBlockTridiagonalMatrix(&(BX[0]));
  FreeBlockTridiagonalMatrix(&(BX[1]));
  if(BC.QComplex) ComplexBlockTridiagonalMatrixToReal(&BC);
  C->N = BC.N;
  if(InitTriDiagonalMatrix(C))
  {
    FreeBlockTridiagonalMatrix(&BC);
    printf("InitTriDiagonalMatrix failed in TriDiagonalMatrixAdd\n");
    return 1;
  }
  for(k=0;k<=BC.N;k++) C->a[k] = BC.aRe[k];
  for(k=0;k<BC.N;k++)  C->b[k] = BC.bRe[k];
  C->mu = BC.mu;
  FreeBlockTridiagonalMatrix(&BC);
  return 0;
}
