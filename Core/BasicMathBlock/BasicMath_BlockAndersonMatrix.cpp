/*
 *   Initialize a struct of BlockAndersonMatrix type
 *   malloc the space for storing N+1 matrices of dimension BlockSize * BlockSize storing a0 to aN
 *   malloc the space for storing N   matrices of dimension BlockSize * BlockSize storing b0 to bN-1
 */
uint32_t InitBlockAndersonMatrix(BlockAndersonMatrixType *const A)
{
  int i;
  strcpy(A->Name, "Block Anderson Matrix");
  if(A->QComplex)
  {
    A->a= (Complex16*)calloc((A->N+1)*(A->BlockSize*A->BlockSize),sizeof(Complex16));
    if(A->a==NULL)
    {
      printf("malloc failed in InitBlockTridiagonalMatrix\n");
      return 1;
    }
    A->b= (Complex16*)calloc((A->N)*(A->BlockSize*A->BlockSize),sizeof(Complex16));
    if(A->b==NULL)
    {
      free(A->a);
      printf("malloc failed in InitBlockTridiagonalMatrix\n");
      return 1;
    }
    for(i=0;i<A->BlockSize;i++)
    {
      A->b[i+i*A->BlockSize].Re=1;
    }
    
  }
  else
  {
    A->aRe= (double*)calloc((A->N+1)*(A->BlockSize*A->BlockSize),sizeof(double));
    if(A->aRe==NULL)
    {
      printf("malloc failed in InitBlockTridiagonalMatrix\n");
      return 1;
    }
    A->bRe= (double*)calloc((A->N)*(A->BlockSize*A->BlockSize),sizeof(double));
    if(A->bRe==NULL)
    {
      free(A->aRe);
      printf("malloc failed in InitBlockTridiagonalMatrix\n");
      return 1;
    }
    for(i=0;i<A->BlockSize;i++)
    {
      A->bRe[i+i*A->BlockSize]=1;
    }
  }
  A->mu=0;
  return 0;
}

uint32_t FreeBlockAndersonMatrix(BlockAndersonMatrixType *A)
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
  return 0;
}

uint32_t ComplexBlockAndersonMatrixToReal(BlockAndersonMatrixType *const A)
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

uint32_t BlockAndersonMatrixChop(BlockAndersonMatrixType *const A, double const Epsilon)
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
    if(!KeepComplex) ComplexBlockAndersonMatrixToReal(A);
  }
  else
  {
    for(i=0;i<A->BlockSize*A->BlockSize*(A->N+1);i++)
    {
      if(Abs(A->aRe[i])<Epsilon)A->aRe[i]=0;
    }
  }
  return 0;
}

/*
 *  B <- A
 *  B not initialized on entry
 */
uint32_t CopyBlockAndersonMatrix(BlockAndersonMatrixType A, BlockAndersonMatrixType *B)
{
  *B=A;
  if(B->QComplex)
  {
    B->a= (Complex16*)malloc((B->N+1)*(B->BlockSize*B->BlockSize)*sizeof(Complex16));
    if(B->a==NULL)
    {
      printf("malloc failed in CopyBlockAndersonMatrix\n");
      return 1;
    }
    B->b= (Complex16*)malloc((B->N)*(B->BlockSize*B->BlockSize)*sizeof(Complex16));
    if(B->b==NULL)
    {
      free(B->a);
      printf("malloc failed in CopyBlockAndersonMatrix\n");
      return 1;
    }
    memcpy(B->a, A.a, (B->N+1)*(B->BlockSize*B->BlockSize)*sizeof(Complex16));
    memcpy(B->b, A.b, (B->N  )*(B->BlockSize*B->BlockSize)*sizeof(Complex16));
  }
  else
  {
    B->aRe= (double*)malloc((B->N+1)*(B->BlockSize*B->BlockSize)*sizeof(double));
    if(B->aRe==NULL)
    {
      printf("malloc failed in CopyBlockAndersonMatrix\n");
      return 1;
    }
    B->bRe= (double*)malloc((B->N)*(B->BlockSize*B->BlockSize)*sizeof(double));
    if(B->bRe==NULL)
    {
      free(B->aRe);
      printf("malloc failed in CopyBlockAndersonMatrix\n");
      return 1;
    }
    memcpy(B->a, A.a, (B->N+1)*(B->BlockSize*B->BlockSize)*sizeof(double));
    memcpy(B->b, A.b, (B->N  )*(B->BlockSize*B->BlockSize)*sizeof(double));
  }
  return 0;
}
