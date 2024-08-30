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
    if(A->N!=0) for(i=0;i<A->BlockSize;i++)
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
    for(i=0;i<A->BlockSize*A->BlockSize*(A->N);i++)
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
    for(i=A->BlockSize*A->BlockSize*(A->N);i<A->BlockSize*A->BlockSize*(A->N+1);i++)
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
    for(i=0;i<A->BlockSize*A->BlockSize*(A->N);i++)
    {
      if(Abs(A->aRe[i])<Epsilon)A->aRe[i]=0;
      if(Abs(A->bRe[i])<Epsilon)A->bRe[i]=0;
    }
    for(i=A->BlockSize*A->BlockSize*(A->N);i<A->BlockSize*A->BlockSize*(A->N+1);i++)
    {
      if(Abs(A->aRe[i])<Epsilon)A->aRe[i]=0;
    }
  }
  return 0;
}

uint32_t BlockAndersonMatrixChopReduceDimension(BlockAndersonMatrixType *const A, double const Epsilon)
{
  uint32_t i,j;
  bool KeepComplex=false, KeepBlock;
  if(A->QComplex)
  {
    
  }
  else
  {
    for(j=0;j<Sqr(A->BlockSize);j++)
    {
      if(Abs(A->aRe[j])<Epsilon)A->aRe[j]=0;
    }
    if(A->N==0) return 0;
    KeepBlock=false;
    for(j=0;j<Sqr(A->BlockSize);j++)
    {
      if(Abs(A->aRe[Sqr(A->BlockSize)+j])<Epsilon)A->aRe[Sqr(A->BlockSize)+j]=0;
      if(Abs(A->bRe[j])<Epsilon)
      {
        A->bRe[j]=0;
      }
      else
      {
        KeepBlock=true;
      }
    }
    if(KeepBlock==false)
    {
      A->N=0;
    }
    for(i=1;i<A->N;i++)
    {
      KeepBlock=false;
      for(j=0;j<Sqr(A->BlockSize);j++)
      {
        if(Abs(A->aRe[(i+1)*Sqr(A->BlockSize) + j])<Epsilon)A->aRe[(i+1)*Sqr(A->BlockSize) + j]=0;
        if(Abs(A->bRe[i*Sqr(A->BlockSize) + j])<Epsilon)
        {
          A->bRe[i*Sqr(A->BlockSize) + j]=0;
        }
        else
        {
          KeepBlock=true;
        }
      }
      if(KeepBlock==false)
      {
        A->N--;
        if(A->N==i) return 0;
        memcpy(&(A->aRe[(i+1)*Sqr(A->BlockSize)]), &( A->aRe[(i+2)*Sqr(A->BlockSize)]), (A->N - i) * Sqr(A->BlockSize) * sizeof(double));
        memcpy(&(A->bRe[(i  )*Sqr(A->BlockSize)]), &( A->bRe[(i+1)*Sqr(A->BlockSize)]), (A->N - i) * Sqr(A->BlockSize) * sizeof(double));
        i--;
      }
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
    if((B->b==NULL)&&(B->N!=0))
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
    if((B->bRe==NULL)&&(B->N!=0))
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

/*   
 *  The Green's function representation on an Anderson matrix type is defined as
 *  G(w) = A0 + B0^* ( 1 / ( w - A1 - \sum_i ( B_i 1 / ( w - A_{i+1} ) . B_i^+ ) ) ) . B0^T
 *
 *  if A0=0 the inverse is
 *  (B0^T)^-1 . ( w - A1 - \sum_i ( B_i 1 / ( w - A_{i+1} ) . B_i^+ ) ) . (B0^*)^-1
 *
 *  For Block list of poles we now need to transform the matrix A_{i+1} to its eigenbases and split to NBlock different matrices
 *
 *  Take the matrix T_{i+1} as the eigenvectors of A_{i+1} such that
 *  DiagonalMatrix({ a^{i+1}_0, a^{i+1}_1, ..., a^{i+1}_{BlockSize-1} }) = T_{i+1}^* . A_{i+1} . T_{i+1}^T
 *
 */
int BlockAndersonMatrixInverseToListOfPoles(BlockAndersonMatrixType const * const A, BlockListOfPoles * const B, double epsilon)
{
  CompactMatrixType Ai,Aifun,AifunSqr,Bi,Ci;
  double val[A->BlockSize];
  int i,j,ToDiv;
  // we assume A0=0 and B0=1
  B->N=Max(0,(int)(A->N)-1)*A->BlockSize;
  B->QComplex=A->QComplex;
  B->BlockSize=A->BlockSize;
  if(InitBlockListOfPoles(B))
  {
    printf("InitBlockListOfPoles failed in BlockAndersonMatrixInverseToListOfPoles\n");
    return 1;
  }
  strcpy(B->Name, A->Name);
  B->mu=A->mu;
  if(A->N==0)
  {
    printf("BlockAndersonMatrixInverseToListOfPoles failed on inverting a matrix with zero poles\n");
    return 1;
  }
  // add Bi . 1/(w-A_{i+1} . Bi^+
  Ai.Ni=A->BlockSize;
  Ai.Nj=A->BlockSize;
  Ai.QComplex=A->QComplex;
  Ai.Name=(char*)malloc((MAXStrLength+1)*sizeof(char));
  Ai.Name[0]='A'; Ai.Name[1]='\0';
  Bi.Ni=A->BlockSize;
  Bi.Nj=A->BlockSize;
  Bi.QComplex=A->QComplex;
  Bi.Name=(char*)malloc((MAXStrLength+1)*sizeof(char));
  Ai.Name[0]='B'; Ai.Name[1]='\0';
  Ci.Ni=A->BlockSize;
  Ci.Nj=A->BlockSize;
  Ci.QComplex=A->QComplex;
  Ci.Name=(char*)malloc((MAXStrLength+1)*sizeof(char));
  Ci.Name[0]='C'; Ai.Name[1]='\0';
  // we use B->N to count the actual number of poles added. Note that we need to set a0 to the correct pointer in the end
  B->N=0;
  AifunSqr.Ni=A->BlockSize;
  AifunSqr.Nj=A->BlockSize;
  AifunSqr.QComplex=A->QComplex;
  // we reuse the memory for fun_i^T . fun_i^*
  InitCompactMatrix(&AifunSqr);
  // and now the actual loop
  if(A->QComplex==0)
  {
    for(i=1;i<A->N;i++)
    {
      Ai.VRe=&(A->aRe[(i+1)*Sqr(A->BlockSize)]);
      Bi.VRe=&(A->bRe[(i  )*Sqr(A->BlockSize)]);
      CompactMatrixEigensystem(&Ai, &Aifun, val);
      for(j=0;j<A->BlockSize;j++)
      {
        ToDiv=1;
        B->VRe[B->N * (1+Sqr(B->BlockSize))] = val[j];
        SetCompactMatrixToZero(&AifunSqr);
        FunTrDotFunCCPreAllocated(&Aifun, j, &AifunSqr);
        while((j<A->BlockSize)&&(fabs(val[j]-val[j+1])<epsilon))
        {
          j++;
          ToDiv++;
          B->VRe[B->N * (1+Sqr(B->BlockSize))] += val[j];
          FunTrDotFunCCPreAllocated(&Aifun, j, &AifunSqr); // adds spectral weight to the pole
        }
        B->VRe[B->N * (1+Sqr(B->BlockSize))] /= ToDiv;
        Ci.VRe=&(B->VRe[B->N * (1+Sqr(B->BlockSize)) + 1]);
        CompactMatrixMultiply3PreAllocated('N', 'N', 'C', &Bi, &AifunSqr, &Bi, &Ci);
        B->N++;
      }
      FreeCompactMatrix(&Aifun);
    }
    // set a0 to the correct pointer
    B->a0Re = &(B->VRe[B->N * (1+Sqr(B->BlockSize))]);
    //B.a0 becomes A.a1
    memcpy(B->a0Re, &(A->aRe[Sqr(A->BlockSize)]), Sqr(A->BlockSize)*sizeof(double));
  }
  else
  {
    for(i=1;i<A->N;i++)
    {
      Ai.V=&(A->a[(i+1)*Sqr(A->BlockSize)]);
      Bi.V=&(A->b[(i  )*Sqr(A->BlockSize)]);
      CompactMatrixEigensystem(&Ai, &Aifun, val);
      for(j=0;j<A->BlockSize;j++)
      {
        ToDiv=1;
        B->V[B->N * (1+Sqr(B->BlockSize))].Re = val[j];
        B->V[B->N * (1+Sqr(B->BlockSize))].Im = 0;
        SetCompactMatrixToZero(&AifunSqr);
        FunTrDotFunCCPreAllocated(&Aifun, j, &AifunSqr);
        while((j<A->BlockSize)&&(fabs(val[j]-val[j+1])<epsilon))
        {
          j++;
          ToDiv++;
          B->V[B->N * (1+Sqr(B->BlockSize))].Re += val[j];
          FunTrDotFunCCPreAllocated(&Aifun, j, &AifunSqr); // adds spectral weight to the pole
        }
        B->V[B->N * (1+Sqr(B->BlockSize))].Re /= ToDiv;
        Ci.V=&(B->V[B->N * (1+Sqr(B->BlockSize)) + 1]);
        CompactMatrixMultiply3PreAllocated('N', 'N', 'C', &Bi, &AifunSqr, &Bi, &Ci);
        B->N++;
      }
      FreeCompactMatrix(&Aifun);
    }
    // set a0 to the correct pointer
    B->a0 = &(B->V[B->N * (1+Sqr(B->BlockSize))]);
    //B.a0 becomes A.a1
    memcpy(B->a0, &(A->a[Sqr(A->BlockSize)]), Sqr(A->BlockSize)*sizeof(Complex16));
  }
  FreeCompactMatrix(&AifunSqr);
  free(Ai.Name);
  free(Bi.Name);
  free(Ci.Name);
  return 0;
}
