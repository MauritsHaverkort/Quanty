/*
 *   Initialize a struct of BlockListOfPoles type
 *   malloc the space for storing N matrices of dimension BlockSize * BlockSize
 *   malloc the space for storing N doubles
 */
uint32_t InitBlockListOfPoles(BlockListOfPoles *const A)
{
  if(A->QComplex)
  {
    A->V=calloc((A->N+1)*(1+A->BlockSize*A->BlockSize)-1,sizeof(Complex16));
    if(A->V==NULL)
    {
      printf("malloc failed in InitBlockListOfPoles\n");
      return 1;
    }
    A->a0=&(A->V[A->N*(1+A->BlockSize*A->BlockSize)]);
  }
  else
  {
    A->VRe=calloc((A->N+1)*(1+A->BlockSize*A->BlockSize)-1,sizeof(double));
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

/*  QComplex changes from 0 to 1
 *  The memory allcoated for doubles is changed to Complex16
 *  The real part is set, the complex part is set to zero
 */
uint32_t MakeRealBlockListOfPolesComplex(BlockListOfPoles *A)
{
  Complex16 *temp;
  uint32_t i;
  temp=calloc((A->N+1)*(1+A->BlockSize*A->BlockSize)-1,sizeof(Complex16));
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
  return 0;
}

/*  QComplex changes from 1 to 0
 *  The memory allcoated for doubles is changed to double
 *  The real part is set, the complex part is neglected
 */
uint32_t MakeComplexBlockListOfPolesReal(BlockListOfPoles *A)
{
  double *temp;
  uint32_t i;
  temp=calloc((A->N+1)*(1+A->BlockSize*A->BlockSize)-1,sizeof(double));
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
    temp=realloc(A->V, ((A->N+B.N+1)*(1+A->BlockSize*A->BlockSize)-1)*sizeof(Complex16));
  }
  else if(B.QComplex)
  {
    if(MakeRealBlockListOfPolesComplex(A)){printf("MakeRealBlockListOfPolesComplex failed in AddBlockListOfPoles\n");return 1;}
    temp=realloc(A->V, (A->N+B.N+1)*A->BlockSize*A->BlockSize*sizeof(Complex16));
  }
  else
  {
    temp=realloc(A->VRe, ((A->N+B.N+1)*(1+A->BlockSize*A->BlockSize)-1)*sizeof(double));
  }
  if(temp==NULL)
  {
    printf("realloc failed in AddBlockListOfPoles\n");
    return 1;
  }
  if(A->QComplex)
  {
    A->V=temp;
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
    A->VRe=temp;
    memcpy(&(A->VRe[(A->N+B.N)*(1+A->BlockSize*A->BlockSize)]), A->a0, A->BlockSize*A->BlockSize*sizeof(double));
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
    temp=realloc(A->V, ((2*A->N+1)*(1+A->BlockSize*A->BlockSize)-1)*sizeof(Complex16));
  }
  else
  {
    temp=realloc(A->VRe, ((2*A->N+1)*(1+A->BlockSize*A->BlockSize)-1)*sizeof(double));
  }
  if(temp==NULL)
  {
    printf("realloc failed in AddBlockListOfPoles\n");
    return 1;
  }
  if(A->QComplex)
  {
    A->V=temp;
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
    A->VRe=temp;
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
  if (Sqr(elem1[0].Re) + Sqr(elem1[0].Im) < Sqr(elem2[0].Re) + Sqr(elem2[0].Im))
    return -1;
  else if (Sqr(elem1[0].Re) + Sqr(elem1[0].Im) > Sqr(elem2[0].Re) + Sqr(elem2[0].Im))
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
    qsort((void *)&(A->V), A->N, (A->BlockSize*A->BlockSize+1)*sizeof(Complex16), (compfn) OrderBlockListOfPolesCompareComplex);
  }
  else
  {
    qsort((void *)&(A->VRe), A->N, (A->BlockSize*A->BlockSize+1)*sizeof(double), (compfn) OrderBlockListOfPolesCompareReal);
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
uint32_t ReduceBlockListOfPoles(BlockListOfPoles * const A, uint32_t N, double const Ef, double (*Norm) (void *, uint32_t, double, double, double) )
{
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
      SmallestNorm=Norm(iBlock,A->BlockSize,a,amin,aplus);
      for(i=1;i<A->N-1;i++)
      {
        iBlock+=A->BlockSize*A->BlockSize+1;
        amin=a;
        a=aplus;
        aplus=iBlock[A->BlockSize*A->BlockSize].Re;
        Normi=Norm(iBlock,A->BlockSize,a,amin,aplus);
        if(Normi<SmallestNorm)iSmallestNorm=i;
      }
      iBlock+=A->BlockSize*A->BlockSize+1;
      amin=a;
      a=aplus;
      aplus=2*a-amin;
      Normi=Norm(iBlock,A->BlockSize,a,amin,aplus);
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
        for(i=1+(1+A->BlockSize*A->BlockSize)*(iSmallestNorm);i<1+(1+A->BlockSize*A->BlockSize)*(iSmallestNorm+1);i++)
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
        for(i=1+(1+A->BlockSize*A->BlockSize)*(iSmallestNorm);i<1+(1+A->BlockSize*A->BlockSize)*(iSmallestNorm+1);i++)
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
          for(i=1+(1+A->BlockSize*A->BlockSize)*(iSmallestNorm);i<1+(1+A->BlockSize*A->BlockSize)*(iSmallestNorm+1);i++)
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
          for(i=1+(1+A->BlockSize*A->BlockSize)*(iSmallestNorm);i<1+(1+A->BlockSize*A->BlockSize)*(iSmallestNorm+1);i++)
          {
            A->V[1+A->BlockSize*A->BlockSize+i].Re+=A->V[i].Re;
            A->V[1+A->BlockSize*A->BlockSize+i].Im+=A->V[i].Im;
          }
        }
        else
        if((aplus-amin)<DBL_MIN*1000)
        {
          Normiplus=0.5;
          Normimin=0.5;
        }
        else
        {
          Normiplus = (a-amin)/(aplus-amin);
          Normimin = 1-Normimin;
        }
        for(i=1+(1+A->BlockSize*A->BlockSize)*(iSmallestNorm);i<1+(1+A->BlockSize*A->BlockSize)*(iSmallestNorm+1);i++)
        {
          A->V[i-1-A->BlockSize*A->BlockSize].Re+=Normimin*A->V[i].Re;
          A->V[i-1-A->BlockSize*A->BlockSize].Im+=Normimin*A->V[i].Im;
          A->V[i+1+A->BlockSize*A->BlockSize].Re+=Normiplus*A->V[i].Re;
          A->V[i+1+A->BlockSize*A->BlockSize].Im+=Normiplus*A->V[i].Im;
        }
      }
      memmove(&(A->V[i*(A->BlockSize*A->BlockSize+1)]), &(A->V[(i+1)*(A->BlockSize*A->BlockSize+1)]),((A->N-i)*(A->BlockSize*A->BlockSize+1)-1)*sizeof(Complex16));
      A->N--;
    }
    A->V=realloc(A->V, ((A->N+1)*(1+A->BlockSize*A->BlockSize)-1)*sizeof(Complex16));
  }
  else
  {
    double *iBlock;
    while(A->N>N)
    {
      /* find smallest Norm */
      iBlock=A->VRe+1;
      a=A->VRe[0];
      aplus=iBlock[A->BlockSize*A->BlockSize];
      amin=2*a-aplus;
      SmallestNorm=Norm(iBlock,A->BlockSize,a,amin,aplus);
      for(i=1;i<A->N-1;i++)
      {
        iBlock+=A->BlockSize*A->BlockSize+1;
        amin=a;
        a=aplus;
        aplus=iBlock[A->BlockSize*A->BlockSize];
        Normi=Norm(iBlock,A->BlockSize,a,amin,aplus);
        if(Normi<SmallestNorm)iSmallestNorm=i;
      }
      iBlock+=A->BlockSize*A->BlockSize+1;
      amin=a;
      a=aplus;
      aplus=2*a-amin;
      Normi=Norm(iBlock,A->BlockSize,a,amin,aplus);
      if(Normi<SmallestNorm)iSmallestNorm=A->N-1;
      /* remove smallest Norm */
      if(iSmallestNorm==0)
      {
        Normi=ListOfPolesNormFlatWeightReal(A->VRe+1+(1+A->BlockSize*A->BlockSize)*(iSmallestNorm), A->BlockSize, 0, 0, 0);
        Normiplus=ListOfPolesNormFlatWeightReal(A->VRe+1+(1+A->BlockSize*A->BlockSize)*(iSmallestNorm+1), A->BlockSize, 0, 0, 0);
        if(Normi+Normiplus<DBL_MIN*1000)
        {
          a=0.5;
        }
        else
        {
          a=Normi/(Normi+Normiplus);
        }
        A->VRe[1+A->BlockSize*A->BlockSize]=(a*A->VRe[0]+(1-a)*A->VRe[1+A->BlockSize*A->BlockSize]);
        for(i=1+(1+A->BlockSize*A->BlockSize)*(iSmallestNorm);i<1+(1+A->BlockSize*A->BlockSize)*(iSmallestNorm+1);i++)
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
        for(i=1+(1+A->BlockSize*A->BlockSize)*(iSmallestNorm);i<1+(1+A->BlockSize*A->BlockSize)*(iSmallestNorm+1);i++)
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
          if(Normi+Normimin<DBL_MIN*1000)
          {
            a=0.5;
          }
          else
          {
            a=Normi/(Normi+Normimin);
          }
          A->VRe[(1+A->BlockSize*A->BlockSize)*(iSmallestNorm-1)]=((1-a)*A->VRe[(1+A->BlockSize*A->BlockSize)*(iSmallestNorm-1)]+a*A->VRe[(1+A->BlockSize*A->BlockSize)*(iSmallestNorm)]);
          for(i=1+(1+A->BlockSize*A->BlockSize)*(iSmallestNorm);i<1+(1+A->BlockSize*A->BlockSize)*(iSmallestNorm+1);i++)
          {
            A->VRe[i-1-A->BlockSize*A->BlockSize]+=A->VRe[i];
          }
        }
        else if((amin<Ef)&&(a>Ef))
        {
          Normi=ListOfPolesNormFlatWeightReal(A->VRe+1+(1+A->BlockSize*A->BlockSize)*(iSmallestNorm), A->BlockSize, 0, 0, 0);
          Normiplus=ListOfPolesNormFlatWeightReal(A->VRe+1+(1+A->BlockSize*A->BlockSize)*(iSmallestNorm+1), A->BlockSize, 0, 0, 0);
          if(Normi+Normiplus<DBL_MIN*1000)
          {
            a=0.5;
          }
          else
          {
            a=Normi/(Normi+Normiplus);
          }
          A->VRe[(1+A->BlockSize*A->BlockSize)*(iSmallestNorm+1)]=(a*A->VRe[(1+A->BlockSize*A->BlockSize)*(iSmallestNorm)]+(1-a)*A->VRe[(1+A->BlockSize*A->BlockSize)*(iSmallestNorm+1)]);
          for(i=1+(1+A->BlockSize*A->BlockSize)*(iSmallestNorm);i<1+(1+A->BlockSize*A->BlockSize)*(iSmallestNorm+1);i++)
          {
            A->VRe[1+A->BlockSize*A->BlockSize+i]+=A->VRe[i];
          }
        }
        else
          if((aplus-amin)<DBL_MIN*1000)
          {
            Normiplus=0.5;
            Normimin=0.5;
          }
          else
          {
            Normiplus = (a-amin)/(aplus-amin);
            Normimin = 1-Normimin;
          }
        for(i=1+(1+A->BlockSize*A->BlockSize)*(iSmallestNorm);i<1+(1+A->BlockSize*A->BlockSize)*(iSmallestNorm+1);i++)
        {
          A->VRe[i-1-A->BlockSize*A->BlockSize]+=Normimin*A->VRe[i];
          A->VRe[i+1+A->BlockSize*A->BlockSize]+=Normiplus*A->VRe[i];
        }
      }
      memmove(&(A->VRe[i*(A->BlockSize*A->BlockSize+1)]), &(A->VRe[(i+1)*(A->BlockSize*A->BlockSize+1)]),((A->N-i)*(A->BlockSize*A->BlockSize+1)-1)*sizeof(double));
      A->N--;
    }
    A->V=realloc(A->V, ((A->N+1)*(1+A->BlockSize*A->BlockSize)-1)*sizeof(double));
  }
  return 0;
}

uint32_t BlockListOfPolesMoveNegativeEnergyToMu(BlockListOfPoles * const A)
{
  /* TO DO */
  return 0;
}

uint32_t BlockListOfPolesMovePositiveEnergyToMu(BlockListOfPoles * const A)
{
  /* TO DO */
  
  return 0;
}


/*   A <- A-B
 *   The energy of the poles in B is adjusted to be on the same "grid" as the poles in A
 *   If A pole in A becomes negative weight it is substracted from the remainder of A
 *   If A has only one pole left the returned weight can be negative
 */
uint32_t SubtractBlockListOfPoles(BlockListOfPoles * const A, BlockListOfPoles const B)
{
  /* TO DO */
  
  return 0;
}