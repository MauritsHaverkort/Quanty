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
