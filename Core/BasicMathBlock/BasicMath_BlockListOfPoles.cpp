/*
 *   Initialize a struct of BlockListOfPoles type
 *   malloc the space for storing N matrices of dimension BlockSize * BlockSize
 *   malloc the space for storing N doubles
 *   QComplex, N and Blocksize need to be set
 */
uint32_t InitBlockListOfPoles(BlockListOfPoles *const A)
{
  strcpy(A->Name, "Block List Of Poles");
  if(A->QComplex)
  {
    A->V= (Complex16*)calloc((A->N+1)*(1+A->BlockSize*A->BlockSize)-1,sizeof(Complex16));
    if(A->V==NULL)
    {
      printf("malloc failed in InitBlockListOfPoles\n");
      return 1;
    }
    A->a0=&(A->V[A->N*(1+A->BlockSize*A->BlockSize)]);
  }
  else
  {
    A->VRe= (double*)calloc((A->N+1)*(1+A->BlockSize*A->BlockSize)-1,sizeof(double));
    if(A->VRe==NULL)
    {
      printf("malloc failed in InitBlockListOfPoles\n");
      return 1;
    }
    A->a0Re=&(A->VRe[A->N*(1+A->BlockSize*A->BlockSize)]);
  }
  A->mu=0;
  return 0;
}

uint32_t FreeBlockListOfPoles(BlockListOfPoles *A)
{
  if(A->QComplex)
  {
    free(A->V);
  }
  else
  {
    free(A->VRe);
  }
  return 0;
}

/*
 *  B <- A
 *  B not initialized on entry
 */
uint32_t CopyBlockListOfPoles(BlockListOfPoles A, BlockListOfPoles *B)
{
  B->N=A.N;
  B->QComplex=A.QComplex;
  B->BlockSize=A.BlockSize;
  if(InitBlockListOfPoles(B)) {printf("InitBlockListOfPoles failed in CopyBlockListOfPoles\n");return 1;}
  strcpy(B->Name, A.Name);
  if(A.QComplex)
  {
    memcpy(B->V, A.V, ((A.N+1)*(1+A.BlockSize*A.BlockSize)-1)*sizeof(Complex16));
  }
  else
  {
    memcpy(B->VRe, A.VRe, ((A.N+1)*(1+A.BlockSize*A.BlockSize)-1)*sizeof(double));
  }
  B->mu=A.mu;
  return 0;
}

uint32_t CopyBlockListOfPoles(BlockListOfPoles A, void *B)
{
	BlockListOfPoles *BTypecast = (BlockListOfPoles*) B;
	return CopyBlockListOfPoles(A, BTypecast);
}

uint32_t CopyBlockListOfPolesCon(BlockListOfPoles A, BlockListOfPoles *B, double const mu, double const T)
{
  uint32_t k,i,kindex=0;
  double occ;
  CopyBlockListOfPoles(A, B);
  if(B->QComplex)
  {
    for(k=0;k<B->N;k++)
    {
      occ=1-FermiFunction(B->V[kindex].Re, mu, T);
      for(i=kindex+1;i<kindex+1+B->BlockSize*B->BlockSize;i++)
      {
        B->V[i].Re*=occ;
        B->V[i].Im*=occ;
      }
      kindex+=1+B->BlockSize*B->BlockSize;
    }
  }
  else
  {
    for(k=0;k<B->N;k++)
    {
      occ=1-FermiFunction(B->VRe[kindex], mu, T);
      for(i=kindex+1;i<kindex+1+B->BlockSize*B->BlockSize;i++)
      {
        B->VRe[i]*=occ;
      }
      kindex+=1+B->BlockSize*B->BlockSize;
    }
  }
  return 0;
}

uint32_t CopyBlockListOfPolesVal(BlockListOfPoles A, BlockListOfPoles *B, double mu, double T)
{
  uint32_t k,i,kindex=0;
  double occ;
  CopyBlockListOfPoles(A, B);
  if(B->QComplex)
  {
    for(k=0;k<B->N;k++)
    {
      occ=FermiFunction(B->V[kindex].Re, mu, T);
      for(i=kindex+1;i<kindex+1+B->BlockSize*B->BlockSize;i++)
      {
        B->V[i].Re*=occ;
        B->V[i].Im*=occ;
      }
      kindex+=1+B->BlockSize*B->BlockSize;
    }
  }
  else
  {
    for(k=0;k<B->N;k++)
    {
      occ=FermiFunction(B->VRe[kindex], mu, T);
      for(i=kindex+1;i<kindex+1+B->BlockSize*B->BlockSize;i++)
      {
        B->VRe[i]*=occ;
      }
      kindex+=1+B->BlockSize*B->BlockSize;
    }
  }
  return 0;
}

/*
 *   A <- scale * A
 */
uint32_t ScaleBlockListOfPoles(BlockListOfPoles *const A, double const scale)
{
  uint32_t i;
  if(A->QComplex)
  {
    for(i=0;i<A->N*(1+A->BlockSize*A->BlockSize);i++)
    {
      if(i%((1+A->BlockSize*A->BlockSize))==0)i++; /* skip the multiplication of a */
      (A->V[i].Re)*=scale;
      (A->V[i].Im)*=scale;
    }
    for(i=0;i<A->BlockSize*A->BlockSize;i++)
    {
      A->a0[i].Re*=scale;
      A->a0[i].Im*=scale;
    }
  }
  else
  {
    for(i=0;i<A->N*(1+A->BlockSize*A->BlockSize);i++)
    {
      if(i%((1+A->BlockSize*A->BlockSize))==0)i++; /* skip the multiplication of a */
      (A->VRe[i])*=scale;
    }
    for(i=0;i<A->BlockSize*A->BlockSize;i++)
    {
      A->a0Re[i]*=scale;
    }
  }
  return 0;
}

uint32_t ScaleBlockListOfPoles(void *const A, double const scale)
{
	BlockListOfPoles *ATypecast = (BlockListOfPoles*) A;
	return ScaleBlockListOfPoles(ATypecast, scale);
}

/*  QComplex changes from 0 to 1
 *  The memory allcoated for doubles is changed to Complex16
 *  The real part is set, the complex part is set to zero
 */
uint32_t MakeRealBlockListOfPolesComplex(BlockListOfPoles *A)
{
  Complex16 *temp;
  uint32_t i;
  if(A->QComplex) return 0;
  temp= (Complex16*)calloc((A->N+1)*(1+A->BlockSize*A->BlockSize)-1,sizeof(Complex16));
  if(temp==NULL)
  {
    printf("malloc failed in MakeRealBlockListOfPolesComplex\n");
    return 1;
  }
  for(i=0;i<(A->N+1)*(1+A->BlockSize*A->BlockSize)-1;i++)
  {
    temp[i].Re=A->VRe[i];
  }
  free(A->VRe);
  A->V=temp;
  A->QComplex=1;
  A->a0=&(A->V[A->N*(1+A->BlockSize*A->BlockSize)]);
  return 0;
}

/*  QComplex changes from 1 to 0
 *  The memory allcoated for Complex is changed to double
 *  The real part is set, the complex part is neglected
 */
uint32_t MakeComplexBlockListOfPolesReal(BlockListOfPoles *A)
{
  double *temp;
  uint32_t i;
  if(A->QComplex==0) return 0;
  temp= (double*)calloc((A->N+1)*(1+A->BlockSize*A->BlockSize)-1,sizeof(double));
  if(temp==NULL)
  {
    printf("malloc failed in MakeComplexBlockListOfPolesReal\n");
    return 1;
  }
  for(i=0;i<(A->N+1)*(1+A->BlockSize*A->BlockSize)-1;i++)
  {
    temp[i]=A->V[i].Re;
  }
  free(A->V);
  A->VRe=temp;
  A->QComplex=0;
  A->a0Re=&(A->VRe[(A->N)*(1+A->BlockSize*A->BlockSize)]);
  return 0;
}

/*  removes (set to zero) all values smaller then epsilon
 *  If there are no complex entries left the memory allocated is changed and a real BlockListOfPoles is returned
 */
uint32_t ChopBlockListOfPoles(BlockListOfPoles *A, double const epsilon)
{
  uint32_t i,StayComplex=0;
  if(A->QComplex)
  {
    for(i=0;i<((A->N+1)*(1+A->BlockSize*A->BlockSize)-1);i++)
    {
      if(fabs(A->V[i].Re)<epsilon)A->V[i].Re=0;
      if(fabs(A->V[i].Im)<epsilon)
      {
        A->V[i].Im=0;
      }
      else
      {
        StayComplex=1;
      }
    }
    if(StayComplex==0)
    {
      return MakeComplexBlockListOfPolesReal(A);
    }
  }
  else
  {
    for(i=0;i<((A->N+1)*(1+A->BlockSize*A->BlockSize)-1);i++)
    {
      if(fabs(A->VRe[i])<epsilon)A->VRe[i]=0;
    }
  }
  return 0;
}

uint32_t RealBlockListOfPolesMergeDegeneratePoles(BlockListOfPoles *A, double const epsilon)
{
  int i,k;
  double WL,WR;
  OrderBlockListOfPoles(A);
  for(i=A->N-2;i>=0;i--)
  {
    if(A->VRe[(i+1) * (1+Sqr(A->BlockSize))] - A->VRe[(i) * (1+Sqr(A->BlockSize))] < epsilon)
    { // merge pole i+1 with i and remove i+1
      WL=0;WR=0;
      for(k=0;k<A->BlockSize;k++)
      {
        WL+=A->VRe[(i  ) * (1+Sqr(A->BlockSize)) + 1 + k + k * A->BlockSize];
        WR+=A->VRe[(i+1) * (1+Sqr(A->BlockSize)) + 1 + k + k * A->BlockSize];
      }
      if(WR<=DBL_EPSILON*(fabs(WL)+fabs(WR)))
      {
        WL=1;
        WR=0;
      }
      else if(WL<=DBL_EPSILON*(fabs(WL)+fabs(WR)))
      {
        WL=0;
        WR=1;
      }
      else
      {
        WL = WL / (WL + WR);
        WR = 1 - WL;
      }
      A->VRe[(i) * (1+Sqr(A->BlockSize))] = WL * A->VRe[(i) * (1+Sqr(A->BlockSize))] + WR * A->VRe[(i+1) * (1+Sqr(A->BlockSize))];
      for(k=0;k<Sqr(A->BlockSize);k++)
      {
        A->VRe[(i) * (1+Sqr(A->BlockSize)) + 1 + k] += A->VRe[(i+1) * (1+Sqr(A->BlockSize)) + 1 + k];
      }
      memmove(&( A->VRe[ (i+1)*(1+Sqr(A->BlockSize)) ] ), &( A->VRe[ (i+2)*(1+Sqr(A->BlockSize)) ] ), ((  (A->N-i-1) * (1+Sqr(A->BlockSize)) ) - 1)*sizeof(double)  );
      A->N--;
      A->a0Re = A->VRe + A->N * (1+Sqr(A->BlockSize));
    }
  }
  return 0;
}
uint32_t ComplexBlockListOfPolesMergeDegeneratePoles(BlockListOfPoles *A, double const epsilon)
{
  int i,k;
  double WL,WR;
  OrderBlockListOfPoles(A);
  for(i=A->N-2;i>=0;i--)
  {
    if(A->V[(i+1) * (1+Sqr(A->BlockSize))].Re - A->V[(i) * (1+Sqr(A->BlockSize))].Re < epsilon)
    { // merge pole i+1 with i and remove i+1
      WL=0;WR=0;
      for(k=0;k<A->BlockSize;k++)
      {
        WL+=A->V[(i  ) * (1+Sqr(A->BlockSize)) + 1 + k + k * A->BlockSize].Re;
        WR+=A->V[(i+1) * (1+Sqr(A->BlockSize)) + 1 + k + k * A->BlockSize].Re;
      }
      if(WR<=DBL_EPSILON*(fabs(WL)+fabs(WR)))
      {
        WL=1;
        WR=0;
      }
      else if(WL<=DBL_EPSILON*(fabs(WL)+fabs(WR)))
      {
        WL=0;
        WR=1;
      }
      else
      {
        WL = WL / (WL + WR);
        WR = 1 - WL;
      }
      A->V[(i) * (1+Sqr(A->BlockSize))].Re = WL * A->V[(i) * (1+Sqr(A->BlockSize))].Re + WR * A->V[(i+1) * (1+Sqr(A->BlockSize))].Re;
      for(k=0;k<Sqr(A->BlockSize);k++)
      {
        A->V[(i) * (1+Sqr(A->BlockSize)) + 1 + k].Re += A->V[(i+1) * (1+Sqr(A->BlockSize)) + 1 + k].Re;
        A->V[(i) * (1+Sqr(A->BlockSize)) + 1 + k].Im += A->V[(i+1) * (1+Sqr(A->BlockSize)) + 1 + k].Im;
      }
      memmove(&( A->V[ (i+1)*(1+Sqr(A->BlockSize)) ] ), &( A->V[ (i+2)*(1+Sqr(A->BlockSize)) ] ), ((  (A->N-i-1) * (1+Sqr(A->BlockSize)) ) - 1)*sizeof(Complex16)  );
      A->N--;
      A->a0 = A->V + A->N * (1+Sqr(A->BlockSize));
    }
  }
  return 0;
}
uint32_t BlockListOfPolesMergeDegeneratePoles(BlockListOfPoles *A, double const epsilon)
{
  if(A->QComplex==0) return RealBlockListOfPolesMergeDegeneratePoles(A, epsilon);
  return ComplexBlockListOfPolesMergeDegeneratePoles(A, epsilon);
}


/*
 *  A <- A + B
 */
uint32_t AddBlockListOfPoles(BlockListOfPoles * const A, BlockListOfPoles const B)
{
  void * temp;
  uint32_t i;
  if(A->BlockSize!=B.BlockSize)
  {
    printf("Tried to add two block list of poles with different matrix sizes %i and %i\n",A->BlockSize,B.BlockSize);
    printf("Error in AddBlockListOfPoles\n");
    return 1;
  }
  if(A->QComplex)
  {
    temp= realloc(A->V, ((A->N+B.N+1)*(1+A->BlockSize*A->BlockSize)-1)*sizeof(Complex16));
  }
  else if(B.QComplex)
  {
    if(MakeRealBlockListOfPolesComplex(A)){printf("MakeRealBlockListOfPolesComplex failed in AddBlockListOfPoles\n");return 1;}
    temp= realloc(A->V, (A->N+B.N+1)*A->BlockSize*A->BlockSize*sizeof(Complex16));
  }
  else
  {
    temp= (double*)realloc(A->VRe, ((A->N+B.N+1)*(1+A->BlockSize*A->BlockSize)-1)*sizeof(double));
  }
  if(temp==NULL)
  {
    printf("realloc failed in AddBlockListOfPoles\n");
    return 1;
  }
  if(A->QComplex)
  {
    A->V=(Complex16*)temp;
    A->a0=&(A->V[A->N*(1+A->BlockSize*A->BlockSize)]);
    memcpy(&(A->V[(A->N+B.N)*(1+A->BlockSize*A->BlockSize)]), A->a0, A->BlockSize*A->BlockSize*sizeof(Complex16));
    A->a0=&(A->V[(A->N+B.N)*(1+A->BlockSize*A->BlockSize)]);
    if(B.QComplex)
    {
      for(i=0;i<A->BlockSize*A->BlockSize;i++)
      {
        A->a0[i].Re+=B.a0[i].Re;
        A->a0[i].Im+=B.a0[i].Im;
      }
      memcpy(&(A->V[A->N*(1+A->BlockSize*A->BlockSize)]), B.V, B.N*(1+A->BlockSize*A->BlockSize)*sizeof(Complex16));
    }
    else
    {
      for(i=0;i<A->BlockSize*A->BlockSize;i++)
      {
        A->a0[i].Re+=B.a0Re[i];
      }
      for(i=0;i<B.N*(1+A->BlockSize*A->BlockSize);i++)
      {
        A->V[i+A->N*(1+A->BlockSize*A->BlockSize)].Re=B.VRe[i];
        A->V[i+A->N*(1+A->BlockSize*A->BlockSize)].Im=0;
      }
    }
  }
  else
  {
    A->VRe=(double*)temp;
    A->a0Re=&(A->VRe[A->N*(1+A->BlockSize*A->BlockSize)]);
    memcpy(&(A->VRe[(A->N+B.N)*(1+A->BlockSize*A->BlockSize)]), A->a0Re, A->BlockSize*A->BlockSize*sizeof(double));
    A->a0Re=&(A->VRe[(A->N+B.N)*(1+A->BlockSize*A->BlockSize)]);
    for(i=0;i<A->BlockSize*A->BlockSize;i++)
    {
      A->a0Re[i]+=B.a0Re[i];
    }
    memcpy(&(A->VRe[A->N*(1+A->BlockSize*A->BlockSize)]), B.VRe, B.N*(1+A->BlockSize*A->BlockSize)*sizeof(double));
  }
  A->N+=B.N;
  return 0;
}

/* A <- alpha A + (1-alpha) B
 *
 */
uint32_t MixBlockListOfPoles(BlockListOfPoles * const A, BlockListOfPoles const B, double const alpha)
{
  if(alpha < 0)
  {
    printf("Can not mix two BlockListOfPole structures with alpha < 0\n");
    return 1;
  }
  else if(alpha > 1)
  {
    printf("Can not mix two BlockListOfPole structures with alpha > 1\n");
    return 1;
  }
  else if(alpha < DBL_EPSILON)
  {
    FreeBlockListOfPoles(A);
    CopyBlockListOfPoles(B, A);
    return 0;
  }
  else if(1-alpha < DBL_EPSILON)
  {
    return 0;
  }
  ScaleBlockListOfPoles(A, alpha / (1-alpha));
  AddBlockListOfPoles(A, B);
  ScaleBlockListOfPoles(A, 1-alpha);
  return 0;
}

/*  A <- 1/2 ( A(w) + A(mu-w) )
 *
 */
uint32_t SymmetrizeBlockListOfPoles(BlockListOfPoles * const A)
{
  void * temp;
  uint32_t i;
  if(A->QComplex)
  {
    temp= realloc(A->V, ((2*A->N+1)*(1+A->BlockSize*A->BlockSize)-1)*sizeof(Complex16));
  }
  else
  {
    temp= realloc(A->VRe, ((2*A->N+1)*(1+A->BlockSize*A->BlockSize)-1)*sizeof(double));
  }
  if(temp==NULL)
  {
    printf("realloc failed in AddBlockListOfPoles\n");
    return 1;
  }
  if(A->QComplex)
  {
    A->V=(Complex16*)temp;
    memcpy(&(A->V[(2*A->N)*(1+A->BlockSize*A->BlockSize)]), A->a0, A->BlockSize*A->BlockSize*sizeof(Complex16));
    A->a0=&(A->V[(2*A->N)*(1+A->BlockSize*A->BlockSize)]);
    memcpy(&(A->V[A->N*(1+A->BlockSize*A->BlockSize)]), A->V, A->N*(1+A->BlockSize*A->BlockSize)*sizeof(Complex16));
    for(i=0;i<A->N*(1+A->BlockSize*A->BlockSize);i++)
    {
      A->V[i].Re*=-1;
    }
    for(i=0;i<A->BlockSize*A->BlockSize;i++)
    {
      A->a0[i].Re*=2;
      A->a0[i].Im*=2;
    }
  }
  else
  {
    A->VRe=(double*)temp;
    memcpy(&(A->VRe[(2*A->N)*(1+A->BlockSize*A->BlockSize)]), A->a0, A->BlockSize*A->BlockSize*sizeof(double));
    A->a0Re=&(A->VRe[(2*A->N)*(1+A->BlockSize*A->BlockSize)]);
    memcpy(&(A->VRe[A->N*(1+A->BlockSize*A->BlockSize)]), A->VRe, A->N*(1+A->BlockSize*A->BlockSize)*sizeof(double));
    for(i=0;i<A->N*(1+A->BlockSize*A->BlockSize);i++)
    {
      A->VRe[i]*=-1;
    }
    for(i=0;i<A->BlockSize*A->BlockSize;i++)
    {
      A->a0Re[i]*=2;
    }
  }
  A->N*=2;
  ScaleBlockListOfPoles(A, 0.5);
  return 0;
}

int OrderBlockListOfPolesCompareReal(double *elem1, double *elem2)
{
  if (elem1[0] < elem2[0])
    return -1;
  else if (elem1[0] > elem2[0])
    return 1;
  return 0;
}

int OrderBlockListOfPolesCompareComplex(Complex16 *elem1, Complex16 *elem2)
{
  if (elem1[0].Re < elem2[0].Re)
    return -1;
  else if (elem1[0].Re > elem2[0].Re)
    return 1;
  return 0;
}

/*  order the list such that a_{i} <= a_{i+1}
 *
 *
 */
uint32_t OrderBlockListOfPoles(BlockListOfPoles * const A)
{
  if(A->QComplex)
  {
    qsort((void *)(A->V), A->N, (A->BlockSize*A->BlockSize+1)*sizeof(Complex16), (compfn) OrderBlockListOfPolesCompareComplex);
  }
  else
  {
    qsort((void *)(A->VRe), A->N, (A->BlockSize*A->BlockSize+1)*sizeof(double), (compfn) OrderBlockListOfPolesCompareReal);
  }
  return 0;
}

/*
 *  Return as a norm the trace of a matrix (doubles)
 */
double ListOfPolesNormFlatWeightReal(double * b, uint32_t Ntau, double a, double amin, double aplus)
{
  double retval;
  uint32_t i;
  retval=0;
  for(i=0;i<Ntau;i++)
  {
    retval+=b[i*Ntau+i];
  }
  return retval;
}
double ListOfPolesNormFlatWeightComplex(Complex16 * b, uint32_t Ntau, double a, double amin, double aplus)
{
  double retval;
  uint32_t i;
  retval=0;
  for(i=0;i<Ntau;i++)
  {
    retval+=b[i*Ntau+i].Re;
  }
  return retval;
}

/*
 *  Return as a norm the trace of a matrix (doubles) times the distance between the pole on the left and right
 */
double ListOfPolesNormFlatWeightDistanceReal(double * b, uint32_t Ntau, double a, double amin, double aplus)
{
  double retval;
  uint32_t i;
  retval=0;
  for(i=0;i<Ntau;i++)
  {
    retval+=b[i*Ntau+i];
  }
  retval*=(aplus-amin);
  return retval;
}
double ListOfPolesNormFlatWeightDistanceComplex(Complex16 * b, uint32_t Ntau, double a, double amin, double aplus)
{
  double retval;
  uint32_t i;
  retval=0;
  for(i=0;i<Ntau;i++)
  {
    retval+=b[i*Ntau+i].Re;
  }
  retval*=(aplus-amin);
  return retval;
}

/*  reduce the number of poles to N
 *
 */
uint32_t ReduceBlockListOfPoles(BlockListOfPoles * const A, uint32_t N)
{
  double const Ef = A->mu;
  double SmallestNorm=0,Normi,Normiplus,Normimin,a,amin,aplus;
  uint32_t iSmallestNorm=0,i;
  if(N==0) {printf("Reducing a list of poles to zero poles does not make sence\n");return 1;};
  if(A->N<=N) {return 0;}
  OrderBlockListOfPoles(A);
  if(A->QComplex)
  {
    Complex16 *iBlock;
    while(A->N>N)
    {
      /* find smallest Norm */
      iBlock=A->V+1;
      a=A->V[0].Re;
      aplus=iBlock[A->BlockSize*A->BlockSize].Re;
      amin=2*a-aplus;
      SmallestNorm=ListOfPolesNormFlatWeightComplex(iBlock,A->BlockSize,a,amin,aplus);
      iSmallestNorm=0;
      for(i=1;i<A->N-1;i++)
      {
        iBlock+=A->BlockSize*A->BlockSize+1;
        amin=a;
        a=aplus;
        aplus=iBlock[A->BlockSize*A->BlockSize].Re;
        Normi=ListOfPolesNormFlatWeightComplex(iBlock,A->BlockSize,a,amin,aplus);
        if(Normi<SmallestNorm)iSmallestNorm=i;
      }
      iBlock+=A->BlockSize*A->BlockSize+1;
      amin=a;
      a=aplus;
      aplus=2*a-amin;
      Normi=ListOfPolesNormFlatWeightComplex(iBlock,A->BlockSize,a,amin,aplus);
      if(Normi<SmallestNorm)iSmallestNorm=A->N-1;
      /* remove smallest Norm */
      if(iSmallestNorm==0)
      {
        Normi=ListOfPolesNormFlatWeightComplex(A->V+1+(1+A->BlockSize*A->BlockSize)*(iSmallestNorm), A->BlockSize, 0, 0, 0);
        Normiplus=ListOfPolesNormFlatWeightComplex(A->V+1+(1+A->BlockSize*A->BlockSize)*(iSmallestNorm+1), A->BlockSize, 0, 0, 0);
        if(Normi+Normiplus<DBL_MIN*1000)
        {
          a=0.5;
        }
        else
        {
          a=Normi/(Normi+Normiplus);
        }
        A->V[1+A->BlockSize*A->BlockSize].Re=(a*A->V[0].Re+(1-a)*A->V[1+A->BlockSize*A->BlockSize].Re);
        for(i=1+(1+A->BlockSize*A->BlockSize)*(iSmallestNorm);i<(1+A->BlockSize*A->BlockSize)*(iSmallestNorm+1);i++)
        {
          A->V[1+A->BlockSize*A->BlockSize+i].Re+=A->V[i].Re;
          A->V[1+A->BlockSize*A->BlockSize+i].Im+=A->V[i].Im;
        }
      }
      else if(iSmallestNorm==A->N-1)
      {
        Normi=ListOfPolesNormFlatWeightComplex(A->V+1+(1+A->BlockSize*A->BlockSize)*(iSmallestNorm), A->BlockSize, 0, 0, 0);
        Normimin=ListOfPolesNormFlatWeightComplex(A->V+1+(1+A->BlockSize*A->BlockSize)*(iSmallestNorm-1), A->BlockSize, 0, 0, 0);
        if(Normi+Normimin<DBL_MIN*1000)
        {
          a=0.5;
        }
        else
        {
          a=Normi/(Normi+Normimin);
        }
        A->V[(1+A->BlockSize*A->BlockSize)*(iSmallestNorm-1)].Re=((1-a)*A->V[(1+A->BlockSize*A->BlockSize)*(iSmallestNorm-1)].Re+a*A->V[(1+A->BlockSize*A->BlockSize)*(iSmallestNorm)].Re);
        for(i=1+(1+A->BlockSize*A->BlockSize)*(iSmallestNorm);i<(1+A->BlockSize*A->BlockSize)*(iSmallestNorm+1);i++)
        {
          A->V[i-1-A->BlockSize*A->BlockSize].Re+=A->V[i].Re;
          A->V[i-1-A->BlockSize*A->BlockSize].Im+=A->V[i].Im;
        }
      }
      else
      {
        a=A->V[(1+A->BlockSize*A->BlockSize)*(iSmallestNorm)].Re;
        aplus=A->V[(1+A->BlockSize*A->BlockSize)*(iSmallestNorm+1)].Re;
        amin=A->V[(1+A->BlockSize*A->BlockSize)*(iSmallestNorm-1)].Re;
        if((aplus>Ef)&&(a<Ef))
        {
          Normi=ListOfPolesNormFlatWeightComplex(A->V+1+(1+A->BlockSize*A->BlockSize)*(iSmallestNorm), A->BlockSize, 0, 0, 0);
          Normimin=ListOfPolesNormFlatWeightComplex(A->V+1+(1+A->BlockSize*A->BlockSize)*(iSmallestNorm-1), A->BlockSize, 0, 0, 0);
          if(Normi+Normimin<DBL_MIN*1000)
          {
            a=0.5;
          }
          else
          {
            a=Normi/(Normi+Normimin);
          }
          A->V[(1+A->BlockSize*A->BlockSize)*(iSmallestNorm-1)].Re=((1-a)*A->V[(1+A->BlockSize*A->BlockSize)*(iSmallestNorm-1)].Re+a*A->V[(1+A->BlockSize*A->BlockSize)*(iSmallestNorm)].Re);
          for(i=1+(1+A->BlockSize*A->BlockSize)*(iSmallestNorm);i<(1+A->BlockSize*A->BlockSize)*(iSmallestNorm+1);i++)
          {
            A->V[i-1-A->BlockSize*A->BlockSize].Re+=A->V[i].Re;
            A->V[i-1-A->BlockSize*A->BlockSize].Im+=A->V[i].Im;
          }
        }
        else if((amin<Ef)&&(a>Ef))
        {
          Normi=ListOfPolesNormFlatWeightComplex(A->V+1+(1+A->BlockSize*A->BlockSize)*(iSmallestNorm), A->BlockSize, 0, 0, 0);
          Normiplus=ListOfPolesNormFlatWeightComplex(A->V+1+(1+A->BlockSize*A->BlockSize)*(iSmallestNorm+1), A->BlockSize, 0, 0, 0);
          if(Normi+Normiplus<DBL_MIN*1000)
          {
            a=0.5;
          }
          else
          {
            a=Normi/(Normi+Normiplus);
          }
          A->V[(1+A->BlockSize*A->BlockSize)*(iSmallestNorm+1)].Re=(a*A->V[(1+A->BlockSize*A->BlockSize)*(iSmallestNorm)].Re+(1-a)*A->V[(1+A->BlockSize*A->BlockSize)*(iSmallestNorm+1)].Re);
          for(i=1+(1+A->BlockSize*A->BlockSize)*(iSmallestNorm);i<(1+A->BlockSize*A->BlockSize)*(iSmallestNorm+1);i++)
          {
            A->V[1+A->BlockSize*A->BlockSize+i].Re+=A->V[i].Re;
            A->V[1+A->BlockSize*A->BlockSize+i].Im+=A->V[i].Im;
          }
        }
        else
        {
          if((aplus-amin)<DBL_MIN*1000)
          {
            Normiplus=0.5;
            Normimin=0.5;
          }
          else
          {
            Normiplus = (a-amin)/(aplus-amin);
            Normimin = 1-Normiplus;
          }
          for(i=1+(1+A->BlockSize*A->BlockSize)*(iSmallestNorm);i<(1+A->BlockSize*A->BlockSize)*(iSmallestNorm+1);i++)
          {
            A->V[i-1-A->BlockSize*A->BlockSize].Re+=Normimin*A->V[i].Re;
            A->V[i-1-A->BlockSize*A->BlockSize].Im+=Normimin*A->V[i].Im;
            A->V[i+1+A->BlockSize*A->BlockSize].Re+=Normiplus*A->V[i].Re;
            A->V[i+1+A->BlockSize*A->BlockSize].Im+=Normiplus*A->V[i].Im;
          }
        }
      }
      memmove(&(A->V[iSmallestNorm*(A->BlockSize*A->BlockSize+1)]), &(A->V[(iSmallestNorm+1)*(A->BlockSize*A->BlockSize+1)]),((A->N-iSmallestNorm)*(A->BlockSize*A->BlockSize+1)-1)*sizeof(Complex16));
      A->N--;
    }
    A->V= (Complex16*)realloc(A->V, ((A->N+1)*(1+A->BlockSize*A->BlockSize)-1)*sizeof(Complex16));
    A->a0=&(A->V[A->N*(1+A->BlockSize*A->BlockSize)]);
  }
  else /* Real */
  {
    double *iBlock;
    while(A->N>N)
    {
      /* find smallest Norm */
      iBlock=A->VRe+1;
      a=A->VRe[0];
      aplus=iBlock[A->BlockSize*A->BlockSize];
      amin=2*a-aplus;
      SmallestNorm=ListOfPolesNormFlatWeightReal(iBlock,A->BlockSize,a,amin,aplus);
      iSmallestNorm=0;
      for(i=1;i<A->N-1;i++)
      {
        iBlock+=A->BlockSize*A->BlockSize+1;
        amin=a;
        a=aplus;
        aplus=iBlock[A->BlockSize*A->BlockSize];
        Normi=ListOfPolesNormFlatWeightReal(iBlock,A->BlockSize,a,amin,aplus);
        if(Normi<SmallestNorm)iSmallestNorm=i;
      }
      iBlock+=A->BlockSize*A->BlockSize+1;
      amin=a;
      a=aplus;
      aplus=2*a-amin;
      Normi=ListOfPolesNormFlatWeightReal(iBlock,A->BlockSize,a,amin,aplus);
      if(Normi<SmallestNorm)iSmallestNorm=A->N-1;
      /* remove smallest Norm */
      if(iSmallestNorm==0)
      {
        Normi     =ListOfPolesNormFlatWeightReal(A->VRe+1+(1+A->BlockSize*A->BlockSize)*(iSmallestNorm)  , A->BlockSize, 0, 0, 0);
        Normiplus =ListOfPolesNormFlatWeightReal(A->VRe+1+(1+A->BlockSize*A->BlockSize)*(iSmallestNorm+1), A->BlockSize, 0, 0, 0);
        if(Normi+Normiplus<DBL_MIN*1000)
        {
          a=0.5;
        }
        else
        {
          a=Normi/(Normi+Normiplus);
        }
        A->VRe[1+A->BlockSize*A->BlockSize]=(a*A->VRe[0]+(1-a)*A->VRe[1+A->BlockSize*A->BlockSize]);
        for(i=1+(1+A->BlockSize*A->BlockSize)*(iSmallestNorm);i<(1+A->BlockSize*A->BlockSize)*(iSmallestNorm+1);i++)
        {
          A->VRe[1+A->BlockSize*A->BlockSize+i]+=A->VRe[i];
        }
      }
      else if(iSmallestNorm==A->N-1)
      {
        Normi=ListOfPolesNormFlatWeightReal(A->VRe+1+(1+A->BlockSize*A->BlockSize)*(iSmallestNorm), A->BlockSize, 0, 0, 0);
        Normimin=ListOfPolesNormFlatWeightReal(A->VRe+1+(1+A->BlockSize*A->BlockSize)*(iSmallestNorm-1), A->BlockSize, 0, 0, 0);
        if(Normi+Normimin<DBL_MIN*1000)
        {
          a=0.5;
        }
        else
        {
          a=Normi/(Normi+Normimin);
        }
        A->VRe[(1+A->BlockSize*A->BlockSize)*(iSmallestNorm-1)]=((1-a)*A->VRe[(1+A->BlockSize*A->BlockSize)*(iSmallestNorm-1)]+a*A->VRe[(1+A->BlockSize*A->BlockSize)*(iSmallestNorm)]);
        for(i=1+(1+A->BlockSize*A->BlockSize)*(iSmallestNorm);i<(1+A->BlockSize*A->BlockSize)*(iSmallestNorm+1);i++)
        {
          A->VRe[i-1-A->BlockSize*A->BlockSize]+=A->VRe[i];
        }
      }
      else
      {
        a=A->VRe[(1+A->BlockSize*A->BlockSize)*(iSmallestNorm)];
        aplus=A->VRe[(1+A->BlockSize*A->BlockSize)*(iSmallestNorm+1)];
        amin=A->VRe[(1+A->BlockSize*A->BlockSize)*(iSmallestNorm-1)];
        if((aplus>Ef)&&(a<Ef))
        {
          Normi=ListOfPolesNormFlatWeightReal(A->VRe+1+(1+A->BlockSize*A->BlockSize)*(iSmallestNorm), A->BlockSize, 0, 0, 0);
          Normimin=ListOfPolesNormFlatWeightReal(A->VRe+1+(1+A->BlockSize*A->BlockSize)*(iSmallestNorm-1), A->BlockSize, 0, 0, 0);
          Normiplus=1-Normimin;
          if(Normi+Normimin<DBL_MIN*1000)
          {
            a=0.5;
          }
          else
          {
            a=Normi/(Normi+Normimin);
          }
          A->VRe[(1+A->BlockSize*A->BlockSize)*(iSmallestNorm-1)]=((1-a)*A->VRe[(1+A->BlockSize*A->BlockSize)*(iSmallestNorm-1)]+a*A->VRe[(1+A->BlockSize*A->BlockSize)*(iSmallestNorm)]);
          for(i=1+(1+A->BlockSize*A->BlockSize)*(iSmallestNorm);i<(1+A->BlockSize*A->BlockSize)*(iSmallestNorm+1);i++)
          {
            A->VRe[i-1-A->BlockSize*A->BlockSize]+=A->VRe[i];
          }
        }
        else if((amin<Ef)&&(a>Ef))
        {
          Normi=ListOfPolesNormFlatWeightReal(A->VRe+1+(1+A->BlockSize*A->BlockSize)*(iSmallestNorm), A->BlockSize, 0, 0, 0);
          Normiplus=ListOfPolesNormFlatWeightReal(A->VRe+1+(1+A->BlockSize*A->BlockSize)*(iSmallestNorm+1), A->BlockSize, 0, 0, 0);
          Normimin=1-Normiplus;
          if(Normi+Normiplus<DBL_MIN*1000)
          {
            a=0.5;
          }
          else
          {
            a=Normi/(Normi+Normiplus);
          }
          A->VRe[(1+A->BlockSize*A->BlockSize)*(iSmallestNorm+1)]=(a*A->VRe[(1+A->BlockSize*A->BlockSize)*(iSmallestNorm)]+(1-a)*A->VRe[(1+A->BlockSize*A->BlockSize)*(iSmallestNorm+1)]);
          for(i=1+(1+A->BlockSize*A->BlockSize)*(iSmallestNorm);i<(1+A->BlockSize*A->BlockSize)*(iSmallestNorm+1);i++)
          {
            A->VRe[1+A->BlockSize*A->BlockSize+i]+=A->VRe[i];
          }
        }
        else
        {
          if((aplus-amin)<DBL_MIN*1000)
          {
            Normiplus=0.5;
            Normimin=0.5;
          }
          else
          {
            Normiplus = (a-amin)/(aplus-amin);
            Normimin = 1-Normiplus;
          }
          for(i=1+(1+A->BlockSize*A->BlockSize)*(iSmallestNorm);i<(1+A->BlockSize*A->BlockSize)*(iSmallestNorm+1);i++)
          {
            A->VRe[i-1-A->BlockSize*A->BlockSize]+=Normimin*A->VRe[i];
            A->VRe[i+1+A->BlockSize*A->BlockSize]+=Normiplus*A->VRe[i];
          }
        }
      }
      memmove(&(A->VRe[iSmallestNorm*(A->BlockSize*A->BlockSize+1)]), &(A->VRe[(iSmallestNorm+1)*(A->BlockSize*A->BlockSize+1)]),((A->N-iSmallestNorm)*(A->BlockSize*A->BlockSize+1)-1)*sizeof(double));
      A->N--;
    }
    A->VRe= (double*)realloc(A->VRe, ((A->N+1)*(1+A->BlockSize*A->BlockSize)-1)*sizeof(double));
    A->a0Re=&(A->VRe[A->N*(1+A->BlockSize*A->BlockSize)]);
  }
  return 0;
}


uint32_t ReduceBlockListOfPoles(BlockListOfPoles * const A, uint32_t Nval, uint32_t Ncon)
{
  BlockListOfPoles Aval, Acon;
  CopyBlockListOfPolesCon(*A, &Acon, A->mu, 0);
  CopyBlockListOfPolesVal(*A, &Aval, A->mu, 0);
  ReduceBlockListOfPoles(&Aval, Nval);
  ReduceBlockListOfPoles(&Acon, Ncon);
  FreeBlockListOfPoles(A);
  CopyBlockListOfPoles(Aval, A);
  AddBlockListOfPoles(A, Acon);
  return 0;
}


uint32_t BlockListOfPolesPlaceOnMesh(BlockListOfPoles * const A, double const *mesh, int const NE)
{
  Complex16 *temp=NULL;
  double *tempRe=NULL;
  uint32_t i, j, k;
  double de, rl, rr;
  int bs = 1+A->BlockSize*A->BlockSize;

  if(NE<=1)
  {
    printf("Error: The mesh should contain at least 2 points.\nSetting N=2.\n");
    return 1;
  }
  for(i=1;i<NE;i++)
  {
    if(mesh[i]-mesh[i-1]<DBL_EPSILON)
    {
      printf("The mesh are two dense or not in correct order:\nNo. %d - %15.10E \nNo. %d - %15.10E\n", i-1, mesh[i-1], i, mesh[i]);
      return 1;
    }
  }

  if(A->QComplex)
  {
    temp = (Complex16*)calloc((NE+1)*bs-1,sizeof(Complex16));
    if(temp==NULL)
    {
      printf("malloc failed in BlockListOfPolesPlaceOnMesh\n");
      return 1;
    }
    for(i=0;i<NE;i++)
    {
      temp[i*bs].Re = mesh[i];
    }
  }
  else
  {
    tempRe = (double*)calloc((NE+1)*bs-1,sizeof(double));
    if(tempRe==NULL)
    {
      printf("malloc failed in BlockListOfPolesPlaceOnMesh\n");
      return 1;
    }
    for(i=0;i<NE;i++)
    {
      tempRe[i*bs] = mesh[i];
    }
  }


  if(A->QComplex)
  {
    for(i=0;i<A->N;i++)
    {
      if(A->V[i*bs].Re<mesh[0])
      {
        for(k=1;k<bs;k++)
        {
          temp[k].Re += A->V[i*bs+k].Re;
          temp[k].Im += A->V[i*bs+k].Im;
        }
      }
      else if(A->V[i*bs].Re>=mesh[NE-1])
      {
        for(k=1;k<bs;k++)
        {
          temp[(NE-1)*bs+k].Re += A->V[i*bs+k].Re;
          temp[(NE-1)*bs+k].Im += A->V[i*bs+k].Im;
        }
      }
      else
      {
        for(j=0;j<NE-1;j++)
        {
          if( (A->V[i*bs].Re>=mesh[j]) && (A->V[i*bs].Re<mesh[j+1]) )
          {
            de = mesh[j+1] - mesh[j];
            rr = (A->V[i*bs].Re-mesh[j])/de;
            rl = 1.0 - rr;
            for(k=1;k<bs;k++)
            {
              temp[j*bs    +k].Re += rl * A->V[i*bs+k].Re;
              temp[(j+1)*bs+k].Re += rr * A->V[i*bs+k].Re;
              temp[j*bs    +k].Im += rl * A->V[i*bs+k].Im;
              temp[(j+1)*bs+k].Im += rr * A->V[i*bs+k].Im;
            }
          }
        }
      }
    }
    memcpy(&(temp[NE*bs]), A->a0, (bs-1)*sizeof(Complex16));
    free(A->V);
    A->N = NE;
    A->V = temp;
    A->a0=&(A->V[A->N*bs]);
  }
  else
  {
    for(i=0;i<A->N;i++)
    {
      if(A->VRe[i*bs]<mesh[0])
      {
        for(k=1;k<bs;k++)
        {
          tempRe[k] += A->VRe[i*bs+k];
        }
      }
      else if(A->VRe[i*bs]>=mesh[NE-1])
      {
        for(k=1;k<bs;k++)
        {
          tempRe[(NE-1)*bs+k] += A->VRe[i*bs+k];
        }
      }
      else
      {
        for(j=0;j<NE-1;j++)
        {
          if( (A->VRe[i*bs]>=mesh[j]) && (A->VRe[i*bs]<mesh[j+1]) )
          {
            if(mesh[j]<A->mu&&mesh[j+1]>=A->mu)
            {
              if(A->VRe[i*bs]>A->mu)
              {
                rr = 1.0;
                rl = 0.0;
              }
              else
              {
                rr = 0.0;
                rl = 1.0;
              }
            }
            else
            {
            de = mesh[j+1] - mesh[j];
            rr = (A->VRe[i*bs]-mesh[j])/de;
            rl = 1.0 - rr;
            }
            for(k=1;k<bs;k++)
            {
              tempRe[j*bs    +k] += rl * A->VRe[i*bs+k];
              tempRe[(j+1)*bs+k] += rr * A->VRe[i*bs+k];
            }
          }
        }
      }
    }
    memcpy(&(tempRe[NE*bs]), A->a0Re, (bs-1)*sizeof(double));
    free(A->VRe);
    A->N = NE;
    A->VRe = tempRe;
    A->a0Re=&(A->VRe[A->N*bs]);
  }

  return 0;
}

/*  removes all negative spectral weigth
 *  For each pole i we calcualte the eigensystem of Bi
 *  We split Bi to negative and positive eigenvalues.
 *  If Bi is negative definite we remove Bi otherwise we keep Bi with the positive spectral weigth
 *  The negative spectral weight of Bi is given to the poles at B(i-1) and B(i+1)
 *    We do not transfer weight through mu.
 *    We start at i=0 and move the negative weight to i=1
 *    We then move the negative weight of i=1 to i=0 and i=2 keeping the norm fixed
 *    and move the negative weigth of i=0 to i=2
 */
uint32_t MakeBlockListOfPolesPhysical(BlockListOfPoles * const A)
{
  char TempName[] = "Mat";
  int i,j,k;
  bool PosDef, NegDef;
  double TrAiNeg,WL,WR;
  CompactMatrixType Ai, AiPos, AiNeg, Bi;
  int iMuL,iMuR; // include 0 to iMuL-1 on the left and iMuR to A->N-1 on the right
  // when there is no pole do nothing
  if(A->N==0) return 0;
  // set some aditional pointers to make life easier
  Ai.Ni=A->BlockSize;
  Ai.Nj=A->BlockSize;
  Ai.QComplex=A->QComplex;
  Ai.Name = TempName;
  Bi.Ni=A->BlockSize;
  Bi.Nj=A->BlockSize;
  Bi.QComplex=A->QComplex;
  Bi.Name = TempName;
  // make sure A is ordered
  OrderBlockListOfPoles(A);
  // do the hard word depending if A is real or complex
  if(A->QComplex==0)
  {
    // when there is one pole remove the negative part
    if(A->N==1)
    {
      Ai.VRe = &( A->VRe[ 1] );
      HermitianCompactMatrixToPosAndNeg(&Ai, &AiPos, &AiNeg, &PosDef, &NegDef, &TrAiNeg);
      if(PosDef==true)
      {
        FreeCompactMatrix(&AiPos);
        FreeCompactMatrix(&AiNeg);
        return 0;
      }
      printf("\n\n WARNING: Found block list of poles %s with only one pole and negative spectral weight (%e). This is unphysical, if the weight is small that might be numerics, if the weight is large there is probably an error in your physical model. I can not move the weight to neighbouring poles so I removed the negative weight, this does not conserve the norm of G.\n\n",A->Name,TrAiNeg);
      if(NegDef==true)
      {
        A->N=0;
        // move a0
        for(k=0;k<Sqr(A->BlockSize);k++)
        {
          A->VRe[k] = A->a0Re[k];
        }
        A->a0Re = A->VRe;
        FreeCompactMatrix(&AiPos);
        FreeCompactMatrix(&AiNeg);
        return 0;
      }
      // set Ai to AiPos
      for(k=0;k<Sqr(A->BlockSize);k++)
      {
        Ai.VRe[k] = AiPos.VRe[k];
      }
    }
    // there is more than one polein A
    // find iMu
    iMuR=0;
    while(iMuR < A->N)
    {
      if(A->VRe[iMuR * (1+Sqr(A->BlockSize))] > A->mu) break;
      iMuR++;
    }
    iMuL=Min(iMuR,A->N-1);
    while(iMuL >= 0)
    {
      if(A->VRe[iMuL * (1+Sqr(A->BlockSize))] < A->mu) break;
      iMuL--;
    }
    // states at mu
    if(iMuR-iMuL!=1)
    {
      // merge all states together
      Ai.VRe = &( A->VRe[ (iMuL+1)*(1+Sqr(A->BlockSize)) + 1] );
      for(i=iMuL+2;i<iMuR;i++)
      {
        Bi.VRe = &( A->VRe[ (i)*(1+Sqr(A->BlockSize)) + 1] );
        for(k=0;k<Sqr(A->BlockSize);k++)
        {
          Ai.VRe[k]+=Bi.VRe[k];
        }
        // state at iMuL+1 contains all spectral weigth at mu, states iMuL+2 to iMuR-1 if present should be removed lateron
      }
      // remove negative spectral weight from state at mu
         // TODO: I have no good option on how to deal with the negative spectral weigth at mu.
      HermitianCompactMatrixToPosAndNeg(&Ai, &AiPos, &AiNeg, &PosDef, &NegDef, &TrAiNeg);
      if(PosDef==false)
      {
        if(fabs(TrAiNeg)>DBL_EPSILON) printf("\n\n WARNING: Found block list of poles %s with negative spectral weight (%e) at the chemical potential. Negative spectral weight is unphysical, if the weight is small that might be numerics, if the weight is large there is probably an error in your physical model. I do not know if the state at the chemical potential belongs to the valence or the conduction states so I can not move the weight to neighbouring poles. Instead I removed the negative weight, this does not conserve the norm of G.\n\n",A->Name,TrAiNeg);
        if(NegDef==false)
        {
          Ai.VRe = &( A->VRe[ (iMuL+1)*(1+Sqr(A->BlockSize)) + 1] );
          for(k=0;k<Sqr(A->BlockSize);k++)
          {
            Ai.VRe[k]=AiPos.VRe[k];
          }
        }
      }
      // remove states from basis
      if(NegDef==true)
      {
        memmove(&( A->VRe[ (iMuL+1)*(1+Sqr(A->BlockSize)) ] ), &( A->VRe[ (iMuR)*(1+Sqr(A->BlockSize)) ] ), (( (1 + A->N - iMuR) * (1+Sqr(A->BlockSize)) ) - 1)*sizeof(double) );
        A->N -= ( iMuR - iMuL -1);
        A->a0Re = A->VRe + A->N * (1+Sqr(A->BlockSize));
        iMuR = iMuL +1;
      }
      else if(iMuR-iMuL>2)
      {
        memmove(&( A->VRe[ (iMuL+2)*(1+Sqr(A->BlockSize)) ] ), &( A->VRe[ (iMuR)*(1+Sqr(A->BlockSize)) ] ), (( (1 + A->N - iMuR) * (1+Sqr(A->BlockSize)) ) - 1)*sizeof(double)  );
        A->N -= ( iMuR - iMuL -2);
        A->a0Re = A->VRe + A->N * (1+Sqr(A->BlockSize));
        iMuR = iMuL +2;
      }
      FreeCompactMatrix(&AiPos);
      FreeCompactMatrix(&AiNeg);
    }
    // states smaller than mu
    labelStartSmallerThanMu:
    if(iMuL==0)
    {
      Ai.VRe = &( A->VRe[1] );
      HermitianCompactMatrixToPosAndNeg(&Ai, &AiPos, &AiNeg, &PosDef, &NegDef, &TrAiNeg);
      if(PosDef==false)
      {
        if(fabs(TrAiNeg)>DBL_EPSILON) printf("\n\n WARNING: Found block list of poles %s with negative spectral weight (%e) below the chemical potential. Negative spectral weight is unphysical, if the weight is small that might be numerics, if the weight is large there is probably an error in your physical model. There is only one pole below the chemical potential so I can not move the weight to neighbouring poles. Instead I removed the negative weight, this does not conserve the norm of %s.\n\n",A->Name,TrAiNeg, A->Name);
        for(k=0;k<Sqr(A->BlockSize);k++)
        {
          Ai.VRe[k]=AiPos.VRe[k];
        }
      }
      if(NegDef==true)
      {
        // remove the state
        memmove(A->VRe, &( A->VRe[ (1+Sqr(A->BlockSize)) ] ), ((  A->N * (1+Sqr(A->BlockSize)) ) - 1)*sizeof(double)  );
        A->N--;
        A->a0Re = A->VRe + A->N * (1+Sqr(A->BlockSize));
        iMuL--;
        iMuR--;
      }
    }
    else if(iMuL>0)
    {
      // move the negative weigth at iMuL to iMuL-1
      Ai.VRe = &( A->VRe[ (iMuL)*(1+Sqr(A->BlockSize)) + 1] );
      HermitianCompactMatrixToPosAndNeg(&Ai, &AiPos, &AiNeg, &PosDef, &NegDef, &TrAiNeg);
      if(NegDef==true)
      {
        Bi.VRe = &( A->VRe[ (iMuL-1)*(1+Sqr(A->BlockSize)) + 1] );
        for(k=0;k<Sqr(A->BlockSize);k++)
        {
          Bi.VRe[k]+=Ai.VRe[k];
        }
        // remove the state
        memmove(&( A->VRe[ (iMuL)*(1+Sqr(A->BlockSize)) ] ), &( A->VRe[ (iMuL+1)*(1+Sqr(A->BlockSize)) ] ), ((  (A->N-iMuL) * (1+Sqr(A->BlockSize))) - 1)*sizeof(double)  );
        A->N--;
        A->a0Re = A->VRe + A->N * (1+Sqr(A->BlockSize));
        iMuL--;
        iMuR--;
        FreeCompactMatrix(&AiPos);
        FreeCompactMatrix(&AiNeg);
        goto labelStartSmallerThanMu;
      }
      if(PosDef==false)
      {
        Bi.VRe = &( A->VRe[ (iMuL-1)*(1+Sqr(A->BlockSize)) + 1] );
        for(k=0;k<Sqr(A->BlockSize);k++)
        {
          Ai.VRe[k]=AiPos.VRe[k];
          Bi.VRe[k]+=AiNeg.VRe[k];
        }
      }
      FreeCompactMatrix(&AiPos);
      FreeCompactMatrix(&AiNeg);
      // loop over states i from iMuL-1 to 1 and move the negative weigth equal to the left and right
      // in an internal sub loop move from j=i+1 to j=iMuL-1 and move the weight to i-1 or j+1
      // as a last step move the negative weight from iMuL to i-1
      for(i=iMuL-1;i>0;i--)
      {
        // move the negative weigth to i-1 and i+1
        Ai.VRe = &( A->VRe[ (i)*(1+Sqr(A->BlockSize)) + 1] );
        HermitianCompactMatrixToPosAndNeg(&Ai, &AiPos, &AiNeg, &PosDef, &NegDef, &TrAiNeg);
        if(PosDef==false)
        {
          for(k=0;k<Sqr(A->BlockSize);k++)
          {
            Ai.VRe[k]=AiPos.VRe[k];
          }
          WL = A->VRe[ (i+1)*(1+Sqr(A->BlockSize))] - A->VRe[ (i)*(1+Sqr(A->BlockSize))];
          WR = A->VRe[ (i)*(1+Sqr(A->BlockSize))] - A->VRe[ (i-1)*(1+Sqr(A->BlockSize))];
          if(WR<=DBL_EPSILON*(fabs(WL)+fabs(WR)))
          {
            WL=1;
            WR=0;
          }
          else if(WL<=DBL_EPSILON*(fabs(WL)+fabs(WR)))
          {
            WL=0;
            WR=1;
          }
          else
          {
            WL = WL / (WL + WR);
            WR = 1 - WL;
          }
          Ai.VRe = &( A->VRe[ (i-1)*(1+Sqr(A->BlockSize)) + 1] );
          Bi.VRe = &( A->VRe[ (i+1)*(1+Sqr(A->BlockSize)) + 1] );
          for(k=0;k<Sqr(A->BlockSize);k++)
          {
            Ai.VRe[k]+=WL*AiNeg.VRe[k];
            Bi.VRe[k]+=WR*AiNeg.VRe[k];
          }
        }
        FreeCompactMatrix(&AiPos);
        FreeCompactMatrix(&AiNeg);
        // loop from i+1 to iMuL-1
        for(j=i+1;j<iMuL;j++)
        { // move negative weith of site j to site i-1 and site j+1
          Ai.VRe = &( A->VRe[ (j)*(1+Sqr(A->BlockSize)) + 1] );
          HermitianCompactMatrixToPosAndNeg(&Ai, &AiPos, &AiNeg, &PosDef, &NegDef, &TrAiNeg);
          if(PosDef==false)
          {
            for(k=0;k<Sqr(A->BlockSize);k++)
            {
              Ai.VRe[k]=AiPos.VRe[k];
            }
            WL = A->VRe[ (j+1)*(1+Sqr(A->BlockSize))] - A->VRe[ (j)*(1+Sqr(A->BlockSize))];
            WR = A->VRe[ (j)*(1+Sqr(A->BlockSize))] - A->VRe[ (i-1)*(1+Sqr(A->BlockSize))];
            if(WR<=DBL_EPSILON*(fabs(WL)+fabs(WR)))
            {
              WL=1;
              WR=0;
            }
            else if(WL<=DBL_EPSILON*(fabs(WL)+fabs(WR)))
            {
              WL=0;
              WR=1;
            }
            else
            {
              WL = WL / (WL + WR);
              WR = 1 - WL;
            }
            Ai.VRe = &( A->VRe[ (i-1)*(1+Sqr(A->BlockSize)) + 1] );
            Bi.VRe = &( A->VRe[ (j+1)*(1+Sqr(A->BlockSize)) + 1] );
            for(k=0;k<Sqr(A->BlockSize);k++)
            {
              Ai.VRe[k]+=WL*AiNeg.VRe[k];
              Bi.VRe[k]+=WR*AiNeg.VRe[k];
            }
          }
          FreeCompactMatrix(&AiPos);
          FreeCompactMatrix(&AiNeg);
        }
        // move negative weith of site iMuL to site i-1
        Ai.VRe = &( A->VRe[ (iMuL)*(1+Sqr(A->BlockSize)) + 1] );
        HermitianCompactMatrixToPosAndNeg(&Ai, &AiPos, &AiNeg, &PosDef, &NegDef, &TrAiNeg);
        if(PosDef==false)
        {
          Bi.VRe = &( A->VRe[ (i-1)*(1+Sqr(A->BlockSize)) + 1] );
          for(k=0;k<Sqr(A->BlockSize);k++)
          {
            Ai.VRe[k]=AiPos.VRe[k];
            Bi.VRe[k]+=AiNeg.VRe[k];
          }
        }
        FreeCompactMatrix(&AiPos);
        FreeCompactMatrix(&AiNeg);
      }
      // move from 0 to iMuL-1 and only move weigth up
      for(i=0;i<iMuL;i++)
      {
        Ai.VRe = &( A->VRe[ (i)*(1+Sqr(A->BlockSize)) + 1] );
        HermitianCompactMatrixToPosAndNeg(&Ai, &AiPos, &AiNeg, &PosDef, &NegDef, &TrAiNeg);
        if(PosDef==false)
        {
          Bi.VRe = &( A->VRe[ (i+1)*(1+Sqr(A->BlockSize)) + 1] );
          if(NegDef==true)
          {
            for(k=0;k<Sqr(A->BlockSize);k++)
            {
              Bi.VRe[k]+=Ai.VRe[k];
            }
          }
          else
          {
            for(k=0;k<Sqr(A->BlockSize);k++)
            {
              Ai.VRe[k]=AiPos.VRe[k];
              Bi.VRe[k]+=AiNeg.VRe[k];
            }
          }
        }
        FreeCompactMatrix(&AiPos);
        FreeCompactMatrix(&AiNeg);
        if(NegDef==true)
        {
          // remove the state
          memmove(&( A->VRe[ (i)*(1+Sqr(A->BlockSize)) ] ), &( A->VRe[ (i+1)*(1+Sqr(A->BlockSize)) ] ), ((  (A->N-i) * (1+Sqr(A->BlockSize)) ) - 1)*sizeof(double)  );
          A->N--;
          A->a0Re = A->VRe + A->N * (1+Sqr(A->BlockSize));
          iMuL--;
          iMuR--;
          i--;
        }
      }
      // remove negative weith from state at iMuL and neglect it
      Ai.VRe = &( A->VRe[ (iMuL)*(1+Sqr(A->BlockSize)) + 1] );
      HermitianCompactMatrixToPosAndNeg(&Ai, &AiPos, &AiNeg, &PosDef, &NegDef, &TrAiNeg);
      if(PosDef==false)
      {
        if(fabs(TrAiNeg)>DBL_EPSILON) printf("\n\n WARNING: Found block list of poles %s with negative spectral weight (%e) below the chemical potential. Negative spectral weight is unphysical, if the weight is small that might be numerics, if the weight is large there is probably an error in your physical model. I tried to move the negative spectral weight to neighbouring poles, but failed as the negative weight is in different subspaces of the matrix. This algorithm can be optimized to find better solutions.\n\n",A->Name,TrAiNeg);
        for(k=0;k<Sqr(A->BlockSize);k++)
        {
          Ai.VRe[k]=AiPos.VRe[k];
        }
      }
      FreeCompactMatrix(&AiPos);
      FreeCompactMatrix(&AiNeg);
      if(NegDef==true)
      {
        // remove the state
        memmove(&( A->VRe[ (iMuL)*(1+Sqr(A->BlockSize)) ] ), &( A->VRe[ (iMuL+1)*(1+Sqr(A->BlockSize)) ] ), ((  (A->N-iMuL) * (1+Sqr(A->BlockSize)) ) - 1)*sizeof(double)  );
        A->N--;
        A->a0Re = A->VRe + A->N * (1+Sqr(A->BlockSize));
        iMuL--;
        iMuR--;
        i--;
      }
    }
    // states larger than mu
    labelStartLargerThanMu:
    if(iMuR==A->N-1)
    {
      Ai.VRe = &( A->VRe[ (A->N-1) * (1+Sqr(A->BlockSize)) + 1] );
      HermitianCompactMatrixToPosAndNeg(&Ai, &AiPos, &AiNeg, &PosDef, &NegDef, &TrAiNeg);
      if(PosDef==false)
      {
        if(fabs(TrAiNeg)>DBL_EPSILON) printf("\n\n WARNING: Found block list of poles %s with negative spectral weight (%e) above the chemical potential. Negative spectral weight is unphysical, if the weight is small that might be numerics, if the weight is large there is probably an error in your physical model. There is only one pole above the chemical potential so I can not move the weight to neighbouring poles. Instead I removed the negative weight, this does not conserve the norm of %s.\n\n",A->Name,TrAiNeg, A->Name);
        for(k=0;k<Sqr(A->BlockSize);k++)
        {
          Ai.VRe[k]=AiPos.VRe[k];
        }
      }
      if(NegDef==true)
      {
        // remove the state
        memmove(&( A->VRe[ (A->N-1) * (1+Sqr(A->BlockSize)) ] ), &( A->VRe[ (A->N) * (1+Sqr(A->BlockSize)) ] ), Sqr(A->BlockSize)*sizeof(double)  );
        A->N--;
        A->a0Re = A->VRe + A->N * (1+Sqr(A->BlockSize));
      }
    }
    else if(iMuR<A->N-1)
    {
      // move the negative weigth at iMuR to iMuR+1
      Ai.VRe = &( A->VRe[ (iMuR)*(1+Sqr(A->BlockSize)) + 1] );
      HermitianCompactMatrixToPosAndNeg(&Ai, &AiPos, &AiNeg, &PosDef, &NegDef, &TrAiNeg);
      if(NegDef==true)
      {
        Bi.VRe = &( A->VRe[ (iMuR+1)*(1+Sqr(A->BlockSize)) + 1] );
        for(k=0;k<Sqr(A->BlockSize);k++)
        {
          Bi.VRe[k]+=Ai.VRe[k];
        }
        // remove the state
        memmove(&( A->VRe[ (iMuR)*(1+Sqr(A->BlockSize)) ] ), &( A->VRe[ (iMuR+1)*(1+Sqr(A->BlockSize)) ] ), ((  (A->N-iMuR) * (1+Sqr(A->BlockSize))) - 1)*sizeof(double)  );
        A->N--;
        A->a0Re = A->VRe + A->N * (1+Sqr(A->BlockSize));
        FreeCompactMatrix(&AiPos);
        FreeCompactMatrix(&AiNeg);
        goto labelStartLargerThanMu;
      }
      if(PosDef==false)
      {
        Bi.VRe = &( A->VRe[ (iMuR+1)*(1+Sqr(A->BlockSize)) + 1] );
        for(k=0;k<Sqr(A->BlockSize);k++)
        {
          Ai.VRe[k]=AiPos.VRe[k];
          Bi.VRe[k]+=AiNeg.VRe[k];
        }
      }
      FreeCompactMatrix(&AiPos);
      FreeCompactMatrix(&AiNeg);
      // loop over states i from iMuR+1 to N-2 and move the negative weigth equal to the left and right
      // in an internal sub loop move from j=i-1 to j=iMuR+1 (down) and move the weight to i+1 or j-1
      // as a last step move the negative weight from iMuR to i+1
      for(i=iMuR+1;i<A->N-1;i++)
      {
        // move the negative weigth to i-1 and i+1
        Ai.VRe = &( A->VRe[ (i)*(1+Sqr(A->BlockSize)) + 1] );
        HermitianCompactMatrixToPosAndNeg(&Ai, &AiPos, &AiNeg, &PosDef, &NegDef, &TrAiNeg);
        if(PosDef==false)
        {
          for(k=0;k<Sqr(A->BlockSize);k++)
          {
            Ai.VRe[k]=AiPos.VRe[k];
          }
          WL = A->VRe[ (i+1)*(1+Sqr(A->BlockSize))] - A->VRe[ (i)*(1+Sqr(A->BlockSize))];
          WR = A->VRe[ (i)*(1+Sqr(A->BlockSize))] - A->VRe[ (i-1)*(1+Sqr(A->BlockSize))];
          if(WR<=DBL_EPSILON*(fabs(WL)+fabs(WR)))
          {
            WL=1;
            WR=0;
          }
          else if(WL<=DBL_EPSILON*(fabs(WL)+fabs(WR)))
          {
            WL=0;
            WR=1;
          }
          else
          {
            WL = WL / (WL + WR);
            WR = 1 - WL;
          }
          Ai.VRe = &( A->VRe[ (i-1)*(1+Sqr(A->BlockSize)) + 1] );
          Bi.VRe = &( A->VRe[ (i+1)*(1+Sqr(A->BlockSize)) + 1] );
          for(k=0;k<Sqr(A->BlockSize);k++)
          {
            Ai.VRe[k]+=WL*AiNeg.VRe[k];
            Bi.VRe[k]+=WR*AiNeg.VRe[k];
          }
        }
        FreeCompactMatrix(&AiPos);
        FreeCompactMatrix(&AiNeg);
        // loop from i-1 to iMuR+1
        for(j=i-1;j>iMuR;j--)
        { // move negative weith of site j to site i+1 and site j-1
          Ai.VRe = &( A->VRe[ (j)*(1+Sqr(A->BlockSize)) + 1] );
          HermitianCompactMatrixToPosAndNeg(&Ai, &AiPos, &AiNeg, &PosDef, &NegDef, &TrAiNeg);
          if(PosDef==false)
          {
            for(k=0;k<Sqr(A->BlockSize);k++)
            {
              Ai.VRe[k]=AiPos.VRe[k];
            }
            WL = A->VRe[ (i+1)*(1+Sqr(A->BlockSize))] - A->VRe[ (j)*(1+Sqr(A->BlockSize))];
            WR = A->VRe[ (j)*(1+Sqr(A->BlockSize))] - A->VRe[ (j-1)*(1+Sqr(A->BlockSize))];
            if(WR<=DBL_EPSILON*(fabs(WL)+fabs(WR)))
            {
              WL=1;
              WR=0;
            }
            else if(WL<=DBL_EPSILON*(fabs(WL)+fabs(WR)))
            {
              WL=0;
              WR=1;
            }
            else
            {
              WL = WL / (WL + WR);
              WR = 1 - WL;
            }
            Ai.VRe = &( A->VRe[ (j-1)*(1+Sqr(A->BlockSize)) + 1] );
            Bi.VRe = &( A->VRe[ (i+1)*(1+Sqr(A->BlockSize)) + 1] );
            for(k=0;k<Sqr(A->BlockSize);k++)
            {
              Ai.VRe[k]+=WL*AiNeg.VRe[k];
              Bi.VRe[k]+=WR*AiNeg.VRe[k];
            }
          }
          FreeCompactMatrix(&AiPos);
          FreeCompactMatrix(&AiNeg);
        }
        // move negative weith of site iMuR to site i+1
        Ai.VRe = &( A->VRe[ (iMuR)*(1+Sqr(A->BlockSize)) + 1] );
        HermitianCompactMatrixToPosAndNeg(&Ai, &AiPos, &AiNeg, &PosDef, &NegDef, &TrAiNeg);
        if(PosDef==false)
        {
          Bi.VRe = &( A->VRe[ (i+1)*(1+Sqr(A->BlockSize)) + 1] );
          for(k=0;k<Sqr(A->BlockSize);k++)
          {
            Ai.VRe[k]=AiPos.VRe[k];
            Bi.VRe[k]+=AiNeg.VRe[k];
          }
        }
        FreeCompactMatrix(&AiPos);
        FreeCompactMatrix(&AiNeg);
      }
      // move from N-1 to iMuR+1 and only move weigth down
      for(i=A->N-1;i>iMuR;i--)
      {
        Ai.VRe = &( A->VRe[ (i)*(1+Sqr(A->BlockSize)) + 1] );
        HermitianCompactMatrixToPosAndNeg(&Ai, &AiPos, &AiNeg, &PosDef, &NegDef, &TrAiNeg);
        if(PosDef==false)
        {
          Bi.VRe = &( A->VRe[ (i-1)*(1+Sqr(A->BlockSize)) + 1] );
          if(NegDef==true)
          {
            for(k=0;k<Sqr(A->BlockSize);k++)
            {
              Bi.VRe[k]+=Ai.VRe[k];
            }
          }
          else
          {
            for(k=0;k<Sqr(A->BlockSize);k++)
            {
              Ai.VRe[k]=AiPos.VRe[k];
              Bi.VRe[k]+=AiNeg.VRe[k];
            }
          }
        }
        FreeCompactMatrix(&AiPos);
        FreeCompactMatrix(&AiNeg);
        if(NegDef==true)
        {
          // remove the state
          memmove(&( A->VRe[ (i)*(1+Sqr(A->BlockSize)) ] ), &( A->VRe[ (i+1)*(1+Sqr(A->BlockSize)) ] ), ((  (A->N-i) * (1+Sqr(A->BlockSize)) ) - 1)*sizeof(double)  );
          A->N--;
          A->a0Re = A->VRe + A->N * (1+Sqr(A->BlockSize));
        }
      }
      // remove negative weith from state at iMuR and neglect it
      Ai.VRe = &( A->VRe[ (iMuR)*(1+Sqr(A->BlockSize)) + 1] );
      HermitianCompactMatrixToPosAndNeg(&Ai, &AiPos, &AiNeg, &PosDef, &NegDef, &TrAiNeg);
      if(PosDef==false)
      {
        if(fabs(TrAiNeg)>DBL_EPSILON) printf("\n\n WARNING: Found block list of poles %s with negative spectral weight (%e) above the chemical potential. Negative spectral weight is unphysical, if the weight is small that might be numerics, if the weight is large there is probably an error in your physical model. I tried to move the negative spectral weight to neighbouring poles, but failed as the negative weight is in different subspaces of the matrix. This algorithm can be optimized to find better solutions.\n\n",A->Name,TrAiNeg);
        for(k=0;k<Sqr(A->BlockSize);k++)
        {
          Ai.VRe[k]=AiPos.VRe[k];
        }
      }
      FreeCompactMatrix(&AiPos);
      FreeCompactMatrix(&AiNeg);
      if(NegDef==true)
      {
        // remove the state
        memmove(&( A->VRe[ (iMuR)*(1+Sqr(A->BlockSize)) ] ), &( A->VRe[ (iMuR+1)*(1+Sqr(A->BlockSize)) ] ), ((  (A->N-iMuR) * (1+Sqr(A->BlockSize)) ) - 1)*sizeof(double)  );
        A->N--;
        A->a0Re = A->VRe + A->N * (1+Sqr(A->BlockSize));
        i--;
      }
    }
  }
  else
  { // Complex
    
    
    
    
    
    
    
    
    
  }
  return 0;
}


/*   A <- A-B
 *   The energy of the poles in B is adjusted to be on the same "grid" as the poles in A
 *   If A pole in A becomes negative weight it is substracted from the remainder of A
 *   If A has only one pole left the returned weight can be negative
 */
uint32_t SubtractBlockListOfPoles(BlockListOfPoles * const A, BlockListOfPoles const * const B)
{
  int i,iMin,j,k;
  double TrA,TrB,Wl,Wr;
  CompactMatrixType Bi,Ai;
  // tests to prevent crashing
  if(A->BlockSize!=B->BlockSize)
  {
    printf("Tried to subtract two BlockListOfPoles with different matrix sizes %i and %i\n",A->BlockSize,B->BlockSize);
    printf("Error in SubtractBlockListOfPoles\n");
    return 1;
  }
  // make A complex if B is complex
  if((A->QComplex==0)&&(B->QComplex!=0))
  {
    MakeRealBlockListOfPolesComplex(A);
  }
  // set some pointers to make life easier
  Bi.Ni=B->BlockSize;
  Bi.Nj=B->BlockSize;
  Bi.QComplex = B->QComplex;
  Ai.Ni=A->BlockSize;
  Ai.Nj=A->BlockSize;
  Ai.QComplex = A->QComplex;
  // A0
  if(A->QComplex==0)
  { // both real
    for(i=0;i<Sqr(A->BlockSize);i++)
    {
      A->a0Re[i]-=B->a0Re[i];
    }
  }
  else
  { // A complex
    if(B->QComplex==0)
    { // B real
      for(i=0;i<Sqr(A->BlockSize);i++)
      {
        A->a0[i].Re-=B->a0Re[i];
      }
    }
    else
    {
      for(i=0;i<Sqr(A->BlockSize);i++)
      {
        A->a0[i].Re-=B->a0[i].Re;
        A->a0[i].Im-=B->a0[i].Im;
      }
    }
  }
  // special case N=0
  if (A->N == 0)
  {
    if (B->N==0)
    {
      return 0;
    }
    // Add one pole as we are going to return a negative weigth (A->N goes from 0 to 1)
    A->N=1;
    if(A->QComplex==0)
    { // A real
      A->VRe = (double*)realloc(A->VRe, (2*(1+Sqr(A->BlockSize))-1)*sizeof(double));
      if(A->VRe==NULL)
      {
        printf("Fatal Error: allocation failed in SubtractBlockListOfPoles\n");
        return 1;
      }
      A->a0Re = &(A->VRe[2*(1+Sqr(A->BlockSize))]);
      for(i=0;i<Sqr(A->BlockSize);i++)
      {
        A->a0Re[i]=A->VRe[i];
        A->VRe[i]=0;
      }
      A->VRe[Sqr(A->BlockSize)]=0;
    }
    else
    { // A complex
      A->V = (Complex16*)realloc(A->V, (2*(1+Sqr(A->BlockSize))-1)*sizeof(Complex16));
      if(A->V==NULL)
      {
        printf("Fatal Error: allocation failed in SubtractBlockListOfPoles\n");
        return 1;
      }
      A->a0 = &(A->V[2*(1+Sqr(A->BlockSize))]);
      for(i=0;i<Sqr(A->BlockSize);i++)
      {
        A->a0[i]=A->V[i];
        A->V[i].Re=0;
        A->V[i].Im=0;
      }
      A->V[Sqr(A->BlockSize)].Re=0;
      A->V[Sqr(A->BlockSize)].Im=0;
    }
  }
  // special case N=1
  if (A->N == 1)
  {
    if(A->QComplex==0)
    { // A and B real
      Ai.VRe = &(A->VRe[1]);
      TrA=0;
      for(k=0;k<A->BlockSize;k++)
      {
        TrA+=Ai.VRe[k*(A->BlockSize+1)];
      }
      for(j=0;j<B->N;j++)
      {
        Bi.VRe = &(B->VRe[(j)*(1+Sqr(B->BlockSize))+1]);
        TrB=0;
        for(k=0;k<B->BlockSize;k++)
        {
          TrB+=Bi.VRe[k*(B->BlockSize+1)];
        }
        for(k=0;k<Sqr(A->BlockSize);k++)
        {
          Ai.VRe[k]-=Bi.VRe[k];
        }
        if(fabs(TrA - TrB) < DBL_EPSILON)
        {
          A->VRe[0] = 0;
        }
        else
        {
          A->VRe[0] = (TrA * (A->VRe[0]) - TrB * (B->VRe[(j)*(1+Sqr(B->BlockSize))])) / (TrA - TrB);
        }
      }
    }
    else
    {  //A complex
      if(B->QComplex==0)
      { // A complex B real
        Ai.V = &(A->V[1]);
        TrA=0;
        for(k=0;k<A->BlockSize;k++)
        {
          TrA+=Ai.V[k*(A->BlockSize+1)].Re;
        }
        for(j=0;j<B->N;j++)
        {
          Bi.VRe = &(B->VRe[(j)*(1+Sqr(B->BlockSize))+1]);
          TrB=0;
          for(k=0;k<B->BlockSize;k++)
          {
            TrB+=Bi.VRe[k*(B->BlockSize+1)];
          }
          for(k=0;k<Sqr(A->BlockSize);k++)
          {
            Ai.V[k].Re-=Bi.VRe[k];
          }
          if(fabs(TrA - TrB) < DBL_EPSILON)
          {
            A->V[0].Re = 0;
            A->V[0].Im = 0;
          }
          else
          {
            A->V[0].Re = (TrA * (A->V[0].Re) - TrB * (B->VRe[(j)*(1+Sqr(B->BlockSize))])) / (TrA - TrB);
          }
        }
      }
      else
      { // A and B Complex
        Ai.V = &(A->V[1]);
        TrA=0;
        for(k=0;k<A->BlockSize;k++)
        {
          TrA+=Ai.V[k*(A->BlockSize+1)].Re;
        }
        for(j=0;j<B->N;j++)
        {
          Bi.V = &(B->V[(j)*(1+Sqr(B->BlockSize))+1]);
          TrB=0;
          for(k=0;k<B->BlockSize;k++)
          {
            TrB+=Bi.V[k*(B->BlockSize+1)].Re;
          }
          for(k=0;k<Sqr(A->BlockSize);k++)
          {
            Ai.V[k]-=Bi.V[k];
          }
          if(fabs(TrA - TrB) < DBL_EPSILON)
          {
            A->V[0].Re = 0;
            A->V[0].Im = 0;
          }
          else
          {
            A->V[0].Re = (TrA * (A->V[0].Re) - TrB * (B->V[(j)*(1+Sqr(B->BlockSize))].Re)) / (TrA - TrB);
          }
        }
      }
    }
    return 0;
  }
  // N>1
  OrderBlockListOfPoles(A);
  if(A->QComplex==0)
  { // A and B real
    for(i=0;i<B->N;i++)
    {
      Bi.VRe = &(B->VRe[i*(Sqr(B->BlockSize)+1)+1]);
      // find position of B[i]
      iMin=0;
      while(iMin<A->N)
      {
        if(B->VRe[i*(Sqr(B->BlockSize)+1)] < A->VRe[iMin*(Sqr(A->BlockSize)+1)] ) break;
        iMin++;
      }
      if(iMin==0)
      {  // below lowest pole
        Ai.VRe = &(A->VRe[1]);
        for(k=0;k<Sqr(A->BlockSize);k++)
        {
          Ai.VRe[k]-=Bi.VRe[k];
        }
      }
      else if(iMin==A->N)
      { // above highest pole
        Ai.VRe = &(A->VRe[(A->N-1)*(Sqr(B->BlockSize)+1)+1]);
        for(k=0;k<Sqr(A->BlockSize);k++)
        {
          Ai.VRe[k]-=Bi.VRe[k];
        }
      }
      else
      { // inbetween pole iMin-1 and iMin
        Wr = B->VRe[i*(Sqr(B->BlockSize)+1)] - A->VRe[(iMin-1)*(Sqr(B->BlockSize)+1)];
        Wl = A->VRe[(iMin)*(Sqr(B->BlockSize)+1)] - B->VRe[i*(Sqr(B->BlockSize)+1)];
        if(Wr <= DBL_MIN + DBL_EPSILON * (fabs(A->VRe[(iMin)*(Sqr(B->BlockSize)+1)]) + fabs(A->VRe[(iMin-1)*(Sqr(B->BlockSize)+1)])))
        {
          if(Wr > Wl)
          { // both are neglegible, equal seperation
            Wr=0.5;
            Wl=0.5;
          }
          else
          { // all weigth to the left
            Wl=1;
            Wr=0;
          }
        }
        else if(Wl <= DBL_MIN + DBL_EPSILON * (fabs(A->VRe[(iMin)*(Sqr(B->BlockSize)+1)]) + fabs(A->VRe[(iMin-1)*(Sqr(B->BlockSize)+1)])))
        { // all weight to the right
          Wl=0;
          Wr=1;
        }
        else
        {
          Wl = Wl / (Wl + Wr);
          Wr = 1 - Wl;
        }
        Ai.VRe = &(A->VRe[(iMin-1)*(Sqr(B->BlockSize)+1)+1]);
        for(k=0;k<Sqr(A->BlockSize);k++)
        {
          Ai.VRe[k]-=Wl * Bi.VRe[k];
        }
        Ai.VRe = &(A->VRe[(iMin)*(Sqr(B->BlockSize)+1)+1]);
        for(k=0;k<Sqr(A->BlockSize);k++)
        {
          Ai.VRe[k]-=Wr * Bi.VRe[k];
        }
      }
    }
  }
  else
  {
    if(B->QComplex==0)
    { // A Complex and B real
      for(i=0;i<B->N;i++)
      {
        Bi.VRe = &(B->VRe[i*(Sqr(B->BlockSize)+1)+1]);
        // find position of B[i]
        iMin=0;
        while(iMin<A->N)
        {
          if(B->VRe[i*(Sqr(B->BlockSize)+1)] < A->V[iMin*(Sqr(A->BlockSize)+1)].Re ) break;
          iMin++;
        }
        if(iMin==0)
        {  // below lowest pole
          Ai.V = &(A->V[1]);
          for(k=0;k<Sqr(A->BlockSize);k++)
          {
            Ai.V[k].Re-=Bi.VRe[k];
          }
        }
        else if(iMin==A->N)
        { // above highest pole
          Ai.V = &(A->V[(A->N-1)*(Sqr(B->BlockSize)+1)+1]);
          for(k=0;k<Sqr(A->BlockSize);k++)
          {
            Ai.V[k].Re-=Bi.VRe[k];
          }
        }
        else
        { // inbetween pole iMin-1 and iMin
          Wr = B->VRe[i*(Sqr(B->BlockSize)+1)] - A->V[(iMin-1)*(Sqr(B->BlockSize)+1)].Re;
          Wl = A->V[(iMin)*(Sqr(B->BlockSize)+1)].Re - B->VRe[i*(Sqr(B->BlockSize)+1)];
          if(Wr <= DBL_MIN + DBL_EPSILON * (fabs(A->V[(iMin)*(Sqr(B->BlockSize)+1)].Re) + fabs(A->V[(iMin-1)*(Sqr(B->BlockSize)+1)].Re)))
          {
            if(Wr > Wl)
            { // both are neglegible, equal seperation
              Wr=0.5;
              Wl=0.5;
            }
            else
            { // all weigth to the left
              Wl=1;
              Wr=0;
            }
          }
          else if(Wl <= DBL_MIN + DBL_EPSILON * (fabs(A->V[(iMin)*(Sqr(B->BlockSize)+1)].Re) + fabs(A->V[(iMin-1)*(Sqr(B->BlockSize)+1)].Re)))
          { // all weight to the right
            Wl=0;
            Wr=1;
          }
          else
          {
            Wl = Wl / (Wl + Wr);
            Wr = 1 - Wl;
          }
          Ai.V = &(A->V[(iMin-1)*(Sqr(B->BlockSize)+1)+1]);
          for(k=0;k<Sqr(A->BlockSize);k++)
          {
            Ai.V[k].Re-=Wl * Bi.VRe[k];
          }
          Ai.V = &(A->V[(iMin)*(Sqr(B->BlockSize)+1)+1]);
          for(k=0;k<Sqr(A->BlockSize);k++)
          {
            Ai.V[k].Re-=Wr * Bi.VRe[k];
          }
        }
      }
    }
    else
    { // A and B Complex
      for(i=0;i<B->N;i++)
      {
        Bi.V = &(B->V[i*(Sqr(B->BlockSize)+1)+1]);
        // find position of B[i]
        iMin=0;
        while(iMin<A->N)
        {
          if(B->V[i*(Sqr(B->BlockSize)+1)].Re < A->V[iMin*(Sqr(A->BlockSize)+1)].Re ) break;
          iMin++;
        }
        if(iMin==0)
        {  // below lowest pole
          Ai.V = &(A->V[1]);
          for(k=0;k<Sqr(A->BlockSize);k++)
          {
            Ai.V[k].Re-=Bi.V[k].Re;
            Ai.V[k].Im-=Bi.V[k].Im;
          }
        }
        else if(iMin==A->N)
        { // above highest pole
          Ai.V = &(A->V[(A->N-1)*(Sqr(B->BlockSize)+1)+1]);
          for(k=0;k<Sqr(A->BlockSize);k++)
          {
            Ai.V[k].Re-=Bi.V[k].Re;
            Ai.V[k].Im-=Bi.V[k].Im;
          }
        }
        else
        { // inbetween pole iMin-1 and iMin
          Wr = B->V[i*(Sqr(B->BlockSize)+1)].Re - A->V[(iMin-1)*(Sqr(B->BlockSize)+1)].Re;
          Wl = A->V[(iMin)*(Sqr(B->BlockSize)+1)].Re - B->V[i*(Sqr(B->BlockSize)+1)].Re;
          if(Wr <= DBL_MIN + DBL_EPSILON * (fabs(A->V[(iMin)*(Sqr(B->BlockSize)+1)].Re) + fabs(A->V[(iMin-1)*(Sqr(B->BlockSize)+1)].Re)))
          {
            if(Wr > Wl)
            { // both are neglegible, equal seperation
              Wr=0.5;
              Wl=0.5;
            }
            else
            { // all weigth to the left
              Wl=1;
              Wr=0;
            }
          }
          else if(Wl <= DBL_MIN + DBL_EPSILON * (fabs(A->V[(iMin)*(Sqr(B->BlockSize)+1)].Re) + fabs(A->V[(iMin-1)*(Sqr(B->BlockSize)+1)].Re)))
          { // all weight to the right
            Wl=0;
            Wr=1;
          }
          else
          {
            Wl = Wl / (Wl + Wr);
            Wr = 1 - Wl;
          }
          Ai.V = &(A->V[(iMin-1)*(Sqr(B->BlockSize)+1)+1]);
          for(k=0;k<Sqr(A->BlockSize);k++)
          {
            Ai.V[k].Re-=Wl * Bi.V[k].Re;
            Ai.V[k].Im-=Wl * Bi.V[k].Im;
          }
          Ai.V = &(A->V[(iMin)*(Sqr(B->BlockSize)+1)+1]);
          for(k=0;k<Sqr(A->BlockSize);k++)
          {
            Ai.V[k].Re-=Wr * Bi.V[k].Re;
            Ai.V[k].Im-=Wr * Bi.V[k].Im;
          }
        }
      }
    }
  }
  // now make A physical
  MakeBlockListOfPolesPhysical(A);
  return 0;
}

uint32_t ShiftBlockListOfPoles(BlockListOfPoles * const A, double delta)
{
  uint32_t i;
  if(A->QComplex==0)
  {
    for(i=0;i<(1+A->BlockSize*A->BlockSize)*A->N;i+=(1+A->BlockSize*A->BlockSize))
    {
      A->VRe[i]+=delta;
    }
  }
  else
  {
    for(i=0;i<(1+A->BlockSize*A->BlockSize)*A->N;i+=(1+A->BlockSize*A->BlockSize))
    {
      A->V[i].Re+=delta;
    }
  }
  //A->mu+=delta;
  return 0;
}

/*  Reduce the dimension of the Green's function represented by a block list of poles
 *  The basis transformation is given by
 *  BasisB = T BasisA
 *  B = T^* A T^T
 */
uint32_t BlockListOfPolesReduceDimcT(BlockListOfPoles const * A, BlockListOfPoles * const B, CompactMatrixType const *T)
{
  uint32_t i,iA,iB;
  CompactMatrixType Bnow, Anow;
  if(T->Nj != A->BlockSize)
  {
    printf("BlockListOfPolesReduceDim failed as T.Nj (%i) != A.BlockSize (%i)\n",T->Nj,A->BlockSize);
    return 1;
  }
  B->QComplex = A->QComplex;
  B->BlockSize = T->Ni;
  B->N=A->N;
  InitBlockListOfPoles(B);
  Anow.QComplex = A->QComplex;
  Anow.Ni=A->BlockSize;
  Anow.Nj=A->BlockSize;
  Bnow.QComplex = B->QComplex;
  Bnow.Ni=B->BlockSize;
  Bnow.Nj=B->BlockSize;
  iA=0;
  iB=0;
  if(B->QComplex)
  {
    for(i=0;i<B->N;i++)
    {
      B->V[iB].Re = A->V[iA].Re;
      iA++;
      iB++;
      Anow.V=&(A->V[iA]);
      Bnow.V=&(B->V[iB]);
      CompactMatrixMultiply3PreAllocated('*', 'N', 'T', T, &Anow, T, &Bnow);
      iA+=A->BlockSize*A->BlockSize;
      iB+=B->BlockSize*B->BlockSize;
    }
  }
  else
  {
    for(i=0;i<B->N;i++)
    {
      B->VRe[iB] = A->VRe[iA];
      iA++;
      iB++;
      Anow.VRe=&(A->VRe[iA]);
      Bnow.VRe=&(B->VRe[iB]);
      CompactMatrixMultiply3PreAllocated('*', 'N', 'T', T, &Anow, T, &Bnow);
      iA+=A->BlockSize*A->BlockSize;
      iB+=B->BlockSize*B->BlockSize;
    }
  }
  return 0;
}
uint32_t BlockListOfPolesTransformcT(BlockListOfPoles const * A, BlockListOfPoles * const B, CompactMatrixType const *T)
{
  return(BlockListOfPolesReduceDimcT(A,B,T));
}
uint32_t BlockListOfPolesTransformInlinecT(BlockListOfPoles * const A, CompactMatrixType const *T)
{
  BlockListOfPoles B;
  if(BlockListOfPolesReduceDimcT(A,&B,T)){return 1;}
  FreeBlockListOfPoles(A);
  *A=B;
  return 0;
}
/*  B = T A T^+
*/
uint32_t BlockListOfPolesReduceDimNC(BlockListOfPoles const * A, BlockListOfPoles * const B, CompactMatrixType const *T)
{
  uint32_t i,iA,iB;
  CompactMatrixType Bnow, Anow;
  if(T->Nj != A->BlockSize)
  {
    printf("BlockListOfPolesReduceDim failed as T.Nj (%i) != A.BlockSize (%i)\n",T->Nj,A->BlockSize);
    return 1;
  }
  B->QComplex = A->QComplex;
  B->BlockSize = T->Ni;
  B->N=A->N;
  InitBlockListOfPoles(B);
  Anow.QComplex = A->QComplex;
  Anow.Ni=A->BlockSize;
  Anow.Nj=A->BlockSize;
  Bnow.QComplex = B->QComplex;
  Bnow.Ni=B->BlockSize;
  Bnow.Nj=B->BlockSize;
  iA=0;
  iB=0;
  if(B->QComplex)
  {
    for(i=0;i<B->N;i++)
    {
      B->V[iB].Re = A->V[iA].Re;
      iA++;
      iB++;
      Anow.V=&(A->V[iA]);
      Bnow.V=&(B->V[iB]);
      CompactMatrixMultiply3PreAllocated('N', 'N', 'C', T, &Anow, T, &Bnow);
      iA+=A->BlockSize*A->BlockSize;
      iB+=B->BlockSize*B->BlockSize;
    }
  }
  else
  {
    for(i=0;i<B->N;i++)
    {
      B->VRe[iB] = A->VRe[iA];
      iA++;
      iB++;
      Anow.VRe=&(A->VRe[iA]);
      Bnow.VRe=&(B->VRe[iB]);
      CompactMatrixMultiply3PreAllocated('*', 'N', 'T', T, &Anow, T, &Bnow);
      iA+=A->BlockSize*A->BlockSize;
      iB+=B->BlockSize*B->BlockSize;
    }
  }
  return 0;
}

uint32_t BlockListOfPolesTransformNC(BlockListOfPoles const * A, BlockListOfPoles * const B, CompactMatrixType const *T)
{
  return(BlockListOfPolesReduceDimNC(A,B,T));
}

uint32_t BlockListOfPolesTransformInlineNC(BlockListOfPoles * const A, CompactMatrixType const *T)
{
  BlockListOfPoles B;
  if(BlockListOfPolesReduceDimNC(A,&B,T)){return 1;}
  FreeBlockListOfPoles(A);
  *A=B;
  return 0;
}


/*  Given G0 and Sigma as Block List of Poles calculate the full Green's function also as a Block List of poles
 *    The energy points of G are given on a fixed mesh
 *    G not initiallized on input
 */
uint32_t BlockListOfPolesG0plusSigmaToGFixedEGrid(BlockListOfPoles const *G0, BlockListOfPoles const *Sigma, BlockListOfPoles * const G, double const Emin, double const Emax, uint32_t const NE)
{
  //TODO: this function is wrong??
  printf("This function \"BlockListOfPolesG0plusSigmaToGFixedEGrid\"is wrong ???why???\n");


  CompactMatrixType S,fun, beta, Char, betasqrt;
  double *val, rest, totnow, totmax, totmin, energy;
  double const dE = (Emax - Emin)/((double)NE);
  int iE;
  uint32_t i,j,k, iS, iLS, iRS, iSigma, iG0, ival, ic,jc;
  if(G0->BlockSize!=Sigma->BlockSize)
  {
    printf("Trying to add a self energy of dimensions %i to a bare Green's function of dimensions %i\n",Sigma->BlockSize,G0->BlockSize);
    return 1;
  }
  val= (double*)malloc((Sigma->N+1)*Sigma->BlockSize*sizeof(double));
  Char.Ni=Sigma->BlockSize;
  Char.Nj=Sigma->BlockSize;
  Char.QComplex=0;
  if((Sigma->QComplex)||(G0->QComplex)) Char.QComplex=1;
  InitCompactMatrix(&Char);
  beta.Ni=Sigma->BlockSize;
  beta.Nj=Sigma->BlockSize;
  beta.QComplex=G0->QComplex;
  beta.Name=Char.Name;
  G->BlockSize=G0->BlockSize;
  G->N=NE+3;
  G->QComplex=0;
  if((Sigma->QComplex)||(G0->QComplex)) G->QComplex=1;
  InitBlockListOfPoles(G);
  energy=Emin-dE;
  if(G->QComplex)
  {
    for(i=0;i<G->N*(1+G->BlockSize*G->BlockSize);i+=(1+G->BlockSize*G->BlockSize))
    {
      G->V[i].Re=energy;
      energy+=dE;
    }
  }
  else
  {
    for(i=0;i<G->N*(1+G->BlockSize*G->BlockSize);i+=(1+G->BlockSize*G->BlockSize))
    {
      G->VRe[i]=energy;
      energy+=dE;
    }
  }
  /* create a compact matrix S from the list of poles Sigma of dimension Sigma->N+1 times BlockSize */
  S.QComplex=Sigma->QComplex;
  S.Ni=(Sigma->N+1)*Sigma->BlockSize;
  S.Nj=(Sigma->N+1)*Sigma->BlockSize;
  InitCompactMatrix(&S);
  iS=(S.Nj+1)*Sigma->BlockSize;
  iLS=Sigma->BlockSize;
  iRS=S.Nj*Sigma->BlockSize;
  iSigma=0;
  if(S.QComplex)
  {
    for(i=0;i<Sigma->N;i++)
    {
      for(j=0;j<Sigma->BlockSize;j++)
      {
        S.V[iS].Re=Sigma->V[iSigma].Re;
        S.V[iS].Im=Sigma->V[iSigma].Im;
        iS+=S.Nj+1;
      }
      iSigma+=1;
      for(j=0;j<Sigma->BlockSize;j++)
      {
        for(k=0;k<Sigma->BlockSize;k++)
        {
          S.V[iRS].Re = Sigma->V[iSigma].Re;
          S.V[iLS].Re = Sigma->V[iSigma].Re;
          S.V[iRS].Im = Sigma->V[iSigma].Im;
          S.V[iLS].Im =-Sigma->V[iSigma].Im;
          iRS++;
          iLS+=S.Nj;
          iSigma++;
        }
        iRS+=S.Nj-Sigma->BlockSize;
        iLS-=(Sigma->BlockSize*S.Nj)-1;
      }
    }
    /* add left top block to S */
    iS=0;
    iSigma=0;
    for(j=0;j<Sigma->BlockSize;j++)
    {
      for(k=0;k<Sigma->BlockSize;k++)
      {
        S.V[iS].Re=Sigma->a0[iSigma].Re;
        S.V[iS].Im=Sigma->a0[iSigma].Im;
        iSigma++;
        iS++;
      }
      iS+=S.Nj-Sigma->BlockSize;
    }
  }
  else
  {
    for(i=0;i<Sigma->N;i++)
    {
      for(j=0;j<Sigma->BlockSize;j++)
      {
        S.VRe[iS]=Sigma->VRe[iSigma];
        iS+=S.Nj+1;
      }
      iSigma+=1;
      for(j=0;j<Sigma->BlockSize;j++)
      {
        for(k=0;k<Sigma->BlockSize;k++)
        {
          S.VRe[iRS] = Sigma->VRe[iSigma];
          S.VRe[iLS] = Sigma->VRe[iSigma];
          iRS++;
          iLS+=S.Nj;
          iSigma++;
        }
        iRS+=S.Nj-Sigma->BlockSize;
        iLS-=(Sigma->BlockSize*S.Nj)-1;
      }
    }
    /* add left top block to S */
    iS=0;
    iSigma=0;
    for(j=0;j<Sigma->BlockSize;j++)
    {
      for(k=0;k<Sigma->BlockSize;k++)
      {
        S.VRe[iS]=Sigma->a0Re[iSigma];
        iSigma++;
        iS++;
      }
      iS+=S.Nj-Sigma->BlockSize;
    }
  }
  /* S is set except for the left top BlockSize by BlockSize part, which depends on G0 */
  /* loop over all poles of G0, set left top of S, diagonalize and add new poles to G  */
  iG0=0;
  totmin=0;
  totmax=0;
  for(i=0;i<G0->N;i++)
  {
    iS=0;
    iSigma=0;
    /* add left top block to S */
    if(Sigma->QComplex)
    {
      if(G0->QComplex)
      {
        for(j=0;j<Sigma->BlockSize;j++)
        {
          S.V[iS].Re=Sigma->a0[iSigma].Re+G0->V[iG0].Re;
          S.V[iS].Im=0; /* diagonal of Hermitian matrix */
          iSigma+=Sigma->BlockSize+1;
          iS+=S.Nj+1;
        }
      }
      else
      {
        for(j=0;j<Sigma->BlockSize;j++)
        {
          S.V[iS].Re=Sigma->a0[iSigma].Re+G0->VRe[iG0];
          S.V[iS].Im=0; /* diagonal of Hermitian matrix */
          iSigma+=Sigma->BlockSize+1;
          iS+=S.Nj+1;
        }
      }
    }
    else
    {
      if(G0->QComplex)
      {
        for(j=0;j<Sigma->BlockSize;j++)
        {
          S.VRe[iS]=Sigma->a0Re[iSigma]+G0->V[iG0].Re;
          iSigma+=Sigma->BlockSize+1;
          iS+=S.Nj+1;
        }
      }
      else
      {
        for(j=0;j<Sigma->BlockSize;j++)
        {
          S.VRe[iS]=Sigma->a0Re[iSigma]+G0->VRe[iG0];
          iSigma+=Sigma->BlockSize+1;
          iS+=S.Nj+1;
        }
      }
    }
    if(CompactMatrixEigensystem(&S, &fun, val))
    {
      printf("CompactMatrixEigensystem failed in BlockListOfPolesG0plusSigmaToGFixedEGrid\n");
      return 1;
    }
    iG0++;
    if(G0->QComplex)
    {
      beta.V=&(G0->V[iG0]);
    }
    else
    {
      beta.VRe=&(G0->VRe[iG0]);
    }
    CompactMatrixSqrt(&beta, &betasqrt);
    for(ival=0;ival<Sigma->BlockSize;ival++)
    {
      /* calculate character matrix for eigen-value ival */
      if(fun.QComplex)
      {
        for(ic=0;ic<G0->BlockSize;ic++)
        {
          for(jc=0;jc<G0->BlockSize;jc++)
          {
            Char.V[ic*Char.Nj+jc].Re = fun.V[ic+ival*fun.Nj].Re*fun.V[jc+ival*fun.Nj].Re - fun.V[ic+ival*fun.Nj].Im*fun.V[jc+ival*fun.Nj].Im;
            Char.V[ic*Char.Nj+jc].Im = fun.V[ic+ival*fun.Nj].Im*fun.V[jc+ival*fun.Nj].Re - fun.V[ic+ival*fun.Nj].Re*fun.V[jc+ival*fun.Nj].Im;
          }
        }
      }
      else
      {
        if(Char.QComplex)
        {
          for(ic=0;ic<G0->BlockSize;ic++)
          {
            for(jc=0;jc<G0->BlockSize;jc++)
            {
              Char.V[ic*Char.Nj+jc].Re = fun.VRe[ic+ival*fun.Nj]*fun.VRe[jc+ival*fun.Nj];
              Char.V[ic*Char.Nj+jc].Im = 0;
            }
          }
        }
        else
        {
          for(ic=0;ic<G0->BlockSize;ic++)
          {
            for(jc=0;jc<G0->BlockSize;jc++)
            {
              Char.VRe[ic*Char.Nj+jc]=fun.VRe[ic+ival*fun.Nj]*fun.VRe[jc+ival*fun.Nj];
            }
          }
        }
      }
      TransformCompactMatrixInlinecT(&Char, &betasqrt); //Need to test
      /* determine the position of the pole in the list of poles */
      rest = (val[ival] - Emin) / dE;
      iE = (int)(floor(rest));
      if(G->QComplex)
      {
        if (iE >= NE)
        {
          totnow=0;
          for(ic=0;ic<G0->BlockSize;ic++)
          {
            totnow+=Char.V[ic+ic*Char.Nj].Re;
          }
          if(totmax+totnow<10*DBL_MIN)
          {
            totmax+=10*DBL_MIN;
            totnow+=10*DBL_MIN;
          }
          G->V[(NE+2)*(1+G->BlockSize*G->BlockSize)].Re = (totmax * G->V[(NE+2)*(1+G->BlockSize*G->BlockSize)].Re + totnow * val[ival]) / (totmax + totnow);
          for(ic=0;ic<G0->BlockSize;ic++)
          {
            for(jc=0;jc<G0->BlockSize;jc++)
            {
              G->V[(NE+2)*(1+G->BlockSize*G->BlockSize)+1+jc+ic*G->BlockSize].Re+=Char.V[ic*Char.Nj+jc].Re;
              G->V[(NE+2)*(1+G->BlockSize*G->BlockSize)+1+jc+ic*G->BlockSize].Im+=Char.V[ic*Char.Nj+jc].Im;
            }
          }
          totmax+=totnow;
        }
        else if (iE < 0)
        {
          totnow=0;
          for(ic=0;ic<G0->BlockSize;ic++)
          {
            totnow+=Char.V[ic+ic*Char.Nj].Re;
          }
          if(totmin+totnow<10*DBL_MIN)
          {
            totmin+=10*DBL_MIN;
            totnow+=10*DBL_MIN;
          }
          G->V[0].Re = (totmin * G->V[0].Re + totnow * val[ival]) / (totmin + totnow);
          for(ic=0;ic<G0->BlockSize;ic++)
          {
            for(jc=0;jc<G0->BlockSize;jc++)
            {
              G->V[1+jc+ic*G0->BlockSize].Re+=Char.V[ic*Char.Nj+jc].Re;
              G->V[1+jc+ic*G0->BlockSize].Im+=Char.V[ic*Char.Nj+jc].Im;
            }
          }
          totmin+=totnow;
        }
        else
        {
          rest -= iE;  /* rest is a number between 0 and 1 pole at iE+1 gets 1-rest pole at iE+2 gest rest */
          for(ic=0;ic<G0->BlockSize;ic++)
          {
            for(jc=0;jc<G0->BlockSize;jc++)
            {
              G->V[(iE+1)*(1+G->BlockSize*G->BlockSize)+1+jc+ic*G->BlockSize].Re+=(1-rest)*Char.V[ic*Char.Nj+jc].Re;
              G->V[(iE+2)*(1+G->BlockSize*G->BlockSize)+1+jc+ic*G->BlockSize].Re+=(  rest)*Char.V[ic*Char.Nj+jc].Re;
              G->V[(iE+1)*(1+G->BlockSize*G->BlockSize)+1+jc+ic*G->BlockSize].Im+=(1-rest)*Char.V[ic*Char.Nj+jc].Im;
              G->V[(iE+2)*(1+G->BlockSize*G->BlockSize)+1+jc+ic*G->BlockSize].Im+=(  rest)*Char.V[ic*Char.Nj+jc].Im;
            }
          }
        }
      }
      else /* all are real */
      {
        if (iE >= NE)
        {
          totnow=0;
          for(ic=0;ic<G0->BlockSize;ic++)
          {
            totnow+=Char.VRe[ic+ic*Char.Nj];
          }
          if(totmax+totnow<10*DBL_MIN)
          {
            totmax+=10*DBL_MIN;
            totnow+=10*DBL_MIN;
          }
          G->VRe[(NE+2)*(1+G->BlockSize*G->BlockSize)] = (totmax * G->VRe[(NE+2)*(1+G->BlockSize*G->BlockSize)] + totnow * val[ival]) / (totmax + totnow);
          for(ic=0;ic<G0->BlockSize;ic++)
          {
            for(jc=0;jc<G0->BlockSize;jc++)
            {
              G->VRe[(NE+2)*(1+G->BlockSize*G->BlockSize)+1+jc+ic*G->BlockSize]+=Char.VRe[ic*Char.Nj+jc];
            }
          }
          totmax+=totnow;
        }
        else if (iE < 0)
        {
          totnow=0;
          for(ic=0;ic<G0->BlockSize;ic++)
          {
            totnow+=Char.VRe[ic+ic*Char.Nj];
          }
          if(totmin+totnow<10*DBL_MIN)
          {
            totmin+=10*DBL_MIN;
            totnow+=10*DBL_MIN;
          }
          G->VRe[0] = (totmin * G->VRe[0] + totnow * val[ival]) / (totmin + totnow);
          for(ic=0;ic<G0->BlockSize;ic++)
          {
            for(jc=0;jc<G0->BlockSize;jc++)
            {
              G->VRe[1+jc+ic*G0->BlockSize]+=Char.VRe[ic*Char.Nj+jc];
            }
          }
          totmin+=totnow;
        }
        else
        {
          rest -= iE;  /* rest is a number between 0 and 1 pole at iE+1 gets 1-rest pole at iE+2 gest rest */
          for(ic=0;ic<G0->BlockSize;ic++)
          {
            for(jc=0;jc<G0->BlockSize;jc++)
            {
              G->VRe[(iE+1)*(1+G->BlockSize*G->BlockSize)+1+jc+ic*G->BlockSize]+=(1.0-rest)*Char.VRe[ic*Char.Nj+jc];
              G->VRe[(iE+2)*(1+G->BlockSize*G->BlockSize)+1+jc+ic*G->BlockSize]+=(    rest)*Char.VRe[ic*Char.Nj+jc];
            }
          }
        }
      }
    }
    FreeCompactMatrix(&betasqrt);
    iG0+=G0->BlockSize*G0->BlockSize;
  }
  ChopBlockListOfPoles(G, 1000*DBL_EPSILON);
  return 0;
}

