/*    Within this file the Green function is always represented as the resolvent of some matrix
 *    
 *    1) Tri-diagonal matrix representation
 *
 *     ( A[0][0]=a1 | A[0][1]=b1 | 0          | 0          )
 *   A=( A[1][0]=b1 | A[1][1]=a2 | A[1][2]=b2 | 0          )
 *     ( 0          | A[2][1]=b2 | A[2][2]=a3 | A[2][3]=b3 )
 *     ( 0          | 0          | A[3][2]=b3 | A[3][3]=a4 )
 *
 *   G(w) = a0 + b0^2 / (w - a1 - b1^2 / (w - a2 - b2^2 / (w - a3 - b3^2 / (w ... ))))
 *
 *
 *    2) Anderson matrix representation
 *
 *     ( A[0][0]=a1 | A[0][1]=b1 | A[0][2]=b2 | A[0][3]=b3 )
 *   A=( A[1][0]=b1 | A[1][1]=a2 | 0          | 0          )
 *     ( A[2][0]=b2 | 0          | A[2][2]=a3 | 0          )
 *     ( A[3][0]=b3 | 0          | 0          | A[3][3]=a4 )
 *
 *   G(w) = a0 + b0^2 / ( w - a1 - \sum_i ( b_i^2 / ( w - a_{i+1} ) ) )
 *
 *
 *    3) Double Tri-diagonal matrix representation
 *
 *     ( A[0][0]=avc0| A[0][1]=bv0 | 0           | 0           | A[0][4]=bc0 | 0           | 0           )
 *   A=( A[1][0]=bv0 | A[1][1]=av1 | A[1][2]=bv1 | 0           | 0           | 0           | 0           )
 *     ( 0           | A[2][1]=bv1 | A[2][2]=av2 | A[2][3]=bv2 | 0           | 0           | 0           )
 *     ( 0           | 0           | A[3][2]=bv2 | A[3][3]=av3 | 0           | 0           | 0           )
 *     ( A[4][0]=bc0 | 0           | 0           | 0           | A[4][4]=ac1 | A[4][5]=bc1 | 0           )
 *     ( 0           | 0           | 0           | 0           | A[5][4]=bc1 | A[5][5]=ac2 | A[5][6]=bc2 )
 *     ( 0           | 0           | 0           | 0           | 0           | A[6][5]=bc2 | A[6][6]=ac3 )
 *
 *   G(w) = a0 + b0^2 / ( w - avc0 - Gv(w) - Gc(w) )
 *   a0 and b0 stored seperately
 *   Gv(w) = bv0^2 / (w - av1 - bv1^2 / (w - av2 - bv2^2 / (w - av3 - bv3^2 / (w ... ))))
 *   Gc(w) = bc0^2 / (w - ac1 - bc1^2 / (w - ac2 - bc2^2 / (w - ac3 - bc3^2 / (w ... ))))
 */

/*
 *   G_M = G_A + G_B
 *
 *   In order to calculate G_M we define a matrix M' wich has as its left top matrix A and as its right bottom matrix B
 *   This matrix is tri-diagonalized with respect to a start-vector.
 *   The start-vector is given as {A.b[0],0,0,...,B.b[0],0,0,0,...}
 *
 *   input A and B defined, M not initiallized
 */
uint32_t SumGreensfunctionsTriDiagonalRepresentationNonOrtho(TriDiagonalMatrixType const A, TriDiagonalMatrixType const B, TriDiagonalMatrixType *const M)
{
  double vecA[3][A.N],vecB[3][B.N];
  double Atemp,Btemp;
  uint32_t i,j;
  M->N=A.N+B.N;
  if(InitTriDiagonalMatrix(M)){printf("InitTriDiagonalMatrix failed in SumGreensfunctionsTriDiagonalRepresentation\n");fflush(stdout);return 1;};
  M->a[0]=A.a[0]+B.a[0];
  M->b[0]=sqrt(Sqr(A.b[0])+Sqr(B.b[0]));
  for(i=0;i<3;i++)
  {
    for(j=0;j<A.N;j++)
    {
      vecA[i][j]=0;
    }
    for(j=0;j<B.N;j++)
    {
      vecB[i][j]=0;
    }
  }
  if(M->b[0]>10*DBL_MIN)
  {
    vecA[0][0]=A.b[0]/M->b[0];
    vecB[0][0]=B.b[0]/M->b[0];
  }
  for(i=1;i<M->N;i++)
  {
    RealTriDiagonalMatrixVector(A, vecA[(i-1)%3], vecA[i%3]);       /* vec[i] = A vec[i-1] */
    RealTriDiagonalMatrixVector(B, vecB[(i-1)%3], vecB[i%3]);
    RealVectorAdd(vecA[i%3],vecA[(i+1)%3],A.N,-(M->b[i-1]));        /* orthogonalization of vec[i] to vec[i-2], using that vec[i] at this point is H vec[i-1] */
    RealVectorAdd(vecB[i%3],vecB[(i+1)%3],B.N,-(M->b[i-1]));
    RealVectorDotVector(vecA[(i-1)%3],vecA[i%3],A.N,&Atemp);        /* a[i] = vec[i] . vec[i-1] */
    RealVectorDotVector(vecB[(i-1)%3],vecB[i%3],B.N,&Btemp);
    M->a[i]=Atemp+Btemp;
    RealVectorAdd(vecA[i%3],vecA[(i-1)%3],A.N,-M->a[i]);            /* vec[i] -= a[i] vec[i-1] */
    RealVectorAdd(vecB[i%3],vecB[(i-1)%3],B.N,-M->a[i]);
    RealVectorDotVector(vecA[i%3],vecA[i%3],A.N,&Atemp);            /* b[i] = sqrt( vec[i] . vec[i] ) */
    RealVectorDotVector(vecB[i%3],vecB[i%3],B.N,&Btemp);
    M->b[i]=sqrt(Atemp+Btemp);
    if(M->b[i]>100*DBL_MIN)                                         /* vec[i] *= 1/b[i] */
    {
      RealVectorMultiply(vecA[i%3],A.N,1.0/M->b[i]);
      RealVectorMultiply(vecB[i%3],B.N,1.0/M->b[i]);
    }
    else
    {
      M->N=i;
      RealVectorMultiply(vecA[i%3],A.N,0);
      RealVectorMultiply(vecB[i%3],B.N,0);
    }
  }
  RealVectorTriDiagonalMatrixVector(A, vecA[(M->N-1)%3], &Atemp);   /* a[N] = vec[N-1] A vec[N-1] */
  RealVectorTriDiagonalMatrixVector(B, vecB[(M->N-1)%3], &Btemp);
  M->a[M->N]=Atemp+Btemp;
  return 0;

}

uint32_t SumGreensfunctionsTriDiagonalRepresentation(TriDiagonalMatrixType const A, TriDiagonalMatrixType const B, TriDiagonalMatrixType *const M)
{
  double vecA[A.N+B.N][A.N];
  double vecB[A.N+B.N][B.N];
  double Atemp,Btemp;
  uint32_t i,j;
  M->N=A.N+B.N;
  if(InitTriDiagonalMatrix(M)){printf("InitTriDiagonalMatrix failed in SumGreensfunctionsTriDiagonalRepresentation\n");fflush(stdout);return 1;};
  M->a[0]=A.a[0]+B.a[0];
  M->b[0]=sqrt(Sqr(A.b[0])+Sqr(B.b[0]));
  for(i=0;i<3;i++)
  {
    for(j=0;j<A.N;j++)
    {
      vecA[i][j]=0;
    }
    for(j=0;j<B.N;j++)
    {
      vecB[i][j]=0;
    }
  }
  if(M->b[0]>10*DBL_MIN)
  {
    vecA[0][0]=A.b[0]/M->b[0];
    vecB[0][0]=B.b[0]/M->b[0];
  }
  for(i=1;i<M->N;i++)
  {
    RealTriDiagonalMatrixVector(A, vecA[(i-1)], vecA[i]);           /* vec[i] = A vec[i-1] */
    RealTriDiagonalMatrixVector(B, vecB[(i-1)], vecB[i]);
    for(j=i;j>0;)
    {
      j--;
      RealVectorDotVector(vecA[j], vecA[i], A.N, &Atemp);           /* orthogonalization of vec[i] to vec[j<i-1] */
      RealVectorDotVector(vecB[j], vecB[i], B.N, &Btemp);
      RealVectorAdd(vecA[i],vecA[j],A.N,-(Atemp+Btemp));
      RealVectorAdd(vecB[i],vecB[j],B.N,-(Atemp+Btemp));
      if(j==i-1)
      {
        M->a[i]=Atemp+Btemp;
      }
    }
    RealVectorDotVector(vecA[i],vecA[i],A.N,&Atemp);                /* b[i] = sqrt( vec[i] . vec[i] ) */
    RealVectorDotVector(vecB[i],vecB[i],B.N,&Btemp);
    M->b[i]=sqrt(Atemp+Btemp);
    if(M->b[i]>100*(DBL_MIN+M->b[1]*DBL_EPSILON))                   /* vec[i] *= 1/b[i] */
    {
      RealVectorMultiply(vecA[i],A.N,1.0/M->b[i]); 
      RealVectorMultiply(vecB[i],B.N,1.0/M->b[i]);
    }
    else
    {
      M->N=i;
      M->b[i]=0;
      RealVectorMultiply(vecA[i],A.N,0); 
      RealVectorMultiply(vecB[i],B.N,0);
    }
  }
  RealVectorTriDiagonalMatrixVector(A, vecA[(M->N-1)], &Atemp);     /* a[N] = vec[N-1] A vec[N-1] */
  RealVectorTriDiagonalMatrixVector(B, vecB[(M->N-1)], &Btemp);
  M->a[M->N]=Atemp+Btemp;
  return 0;
}

uint32_t SymmetrizeAndersonMatrix(AndersonMatrixType *const M)
{
  uint32_t i;
  double *tmp;
  tmp=realloc(M->a, (2*M->N)*sizeof(double));
  if(tmp==NULL){printf("realloc failed in SymmetrizeAndersonMatrix\n");return 1;}
  M->a=tmp;
  tmp=realloc(M->b, (Max(0,2*M->N-1))*sizeof(double));
  if(tmp==NULL){printf("realloc failed in SymmetrizeAndersonMatrix\n");return 1;}
  M->b=tmp;
  for(i=1;i<M->N;i++)
  {
    M->b[i]*=Sqrtp5;
    M->b[i+M->N-1]=M->b[i];
    M->a[i+M->N]=2*M->mu-M->a[i+1];
  }
  M->N*=2;
  if(M->N!=0)M->N--;
  return 0;
}

uint32_t SymmetrizeBADoubleTriDiagonalMatrix(BADoubleTriDiagonalMatrixType *const M)
{
  printf("SymmetrizeBADoubleTriDiagonalMatrix not jet fully implemented\n");
  double tmp;
  tmp=sqrt(Sqr(M->bvB0)+Sqr(M->bcB0))*Sqrtp5;
  if(M->bvB0>=0)
  {
    M->bvB0=tmp;
  }
  else
  {
    M->bvB0=-tmp;
  }
  if(M->bcB0>=0)
  {
    M->bcB0=tmp;
  }
  else
  {
    M->bcB0=-tmp;
  }
  tmp=sqrt(Sqr(M->bvA0)+Sqr(M->bcA0))*Sqrtp5;
  if(M->bvA0>=0)
  {
    M->bvA0=tmp;
  }
  else
  {
    M->bvA0=-tmp;
  }
  if(M->bcA0>=0)
  {
    M->bcA0=tmp;
  }
  else
  {
    M->bcA0=-tmp;
  }
  return 0;
}

uint32_t SymmetrizeTridiagonalMatrix(TriDiagonalMatrixType *const M)
{
  TriDiagonalMatrixType Mplus,Mmin;
  uint32_t i;
  if(TriDiagonalMatrixCopy(&Mplus,*M)){printf("TriDiagonalMatrixCopy failed in SymmetrizeTridiagonalMatrix\n");fflush(stdout);return 1;};
  if(TriDiagonalMatrixCopy(&Mmin ,*M)){printf("TriDiagonalMatrixCopy failed in SymmetrizeTridiagonalMatrix\n");fflush(stdout);return 1;};
  strcpy(Mplus.Name,M->Name);
  strcpy(Mmin.Name,M->Name);
  FreeTriDiagonalMatrix(M);
  for(i=1;i<=Mmin.N;i++)
  {
    Mmin.a[i]*=-1;
  }
  SumGreensfunctionsTriDiagonalRepresentation(Mplus, Mmin, M);
  strcpy(M->Name,Mmin.Name);
  M->b[0]=Mplus.b[0];
  M->a[0]=Mplus.a[0];
  FreeTriDiagonalMatrix(&Mplus);
  FreeTriDiagonalMatrix(&Mmin);
  return 0;
}



/*  Create a new set of poles, which has lenght N, by removing the M->N-N poles with the smallest weight
 *  This algorithm combines the pole with the smallest weight with its nearest two neighbors
 *  and repeats this untill the number of nonzero poles is N
 */
void ReduceSqrtListOfPolesNumberFermi(ListOfPoles * M, uint32_t N, double const Ef)
{
  uint32_t i;
  for(i=0;i<M->N;i++)
  {
    M->b[i]*=M->b[i];
  }
  ReduceListOfPolesNumberFermi(M, N, Ef);
  for(i=0;i<M->N;i++)
  {
    M->b[i]=sqrt(M->b[i]);
  }
}






/*  Create a new set of poles, which has lenght N, by removing the M->N-N poles with the smallest weight
 *  This algorithm combines the pole with the smallest weight with its nearest two neighbors
 *  and repeats this untill the number of nonzero poles is N
 */
void ReduceSqrtListOfPolesNumber(ListOfPoles * M, uint32_t N)
{
  //TimeType Time;
  if(M->N<N) return;
  //InitializeTiming(&Time);
  //BeginTiming("Total ReduceSqrtListOfPolesNumber",&Time);
  uint32_t i,imin;
  double biminLeft,biminRight,temp;
  //BeginTiming("OrderListOfPoles",&Time);
  OrderListOfPoles(M);
  //EndTiming("OrderListOfPoles",&Time);
  while(M->N>N)
  {
    imin=0;
    for(i=1;i<M->N;i++)
    {
      if(M->b[i]<M->b[imin])
      {
        imin=i;
      }
    }
    if(M->b[imin]>DBL_MIN)
    {
      if(imin==0)
      {
        M->a[imin+1]=(Sqr(M->b[imin]/M->b[imin+1])*M->a[imin]+M->a[imin+1])/(Sqr(M->b[imin]/M->b[imin+1])+1);
        M->b[imin+1]=sqrt(Sqr(M->b[imin])+Sqr(M->b[imin+1]));
      }
      else if(imin==M->N-1)
      {
        M->a[imin-1]=(Sqr(M->b[imin]/M->b[imin-1])*M->a[imin]+M->a[imin-1])/(Sqr(M->b[imin]/M->b[imin-1])+1);
        M->b[imin-1]=sqrt(Sqr(M->b[imin])+Sqr(M->b[imin-1]));
      }
      else
      {
        if(fabs((M->a[imin+1]-M->a[imin-1]))<fabs(10*(DBL_MIN+DBL_EPSILON*M->a[imin])))
        {
          biminLeft=M->b[imin]*sqrt(0.5);
          biminRight=M->b[imin]*sqrt(0.5);
        }
        else if(fabs((M->a[imin+1]-M->a[imin]))<fabs(10*(DBL_MIN+DBL_EPSILON*M->a[imin])))
        {
          biminLeft=M->b[imin];
          biminRight=0;
        }
        else if(fabs((M->a[imin]-M->a[imin-1]))<fabs(10*(DBL_MIN+DBL_EPSILON*M->a[imin])))
        {
          biminLeft=0;
          biminRight=M->b[imin];
        }
        else
        {
          biminLeft =M->b[imin]*sqrt((M->a[imin+1]-M->a[imin  ])/(M->a[imin+1]-M->a[imin-1]));
          biminRight=M->b[imin]*sqrt((M->a[imin  ]-M->a[imin-1])/(M->a[imin+1]-M->a[imin-1]));
        }
        M->a[imin+1]=(Sqr(biminRight/M->b[imin+1])*M->a[imin]+M->a[imin+1])/(Sqr(biminRight/M->b[imin+1])+1);
        M->b[imin+1]=sqrt(Sqr(biminRight)+Sqr(M->b[imin+1]));
        M->a[imin-1]=(Sqr(biminLeft /M->b[imin-1])*M->a[imin]+M->a[imin-1])/(Sqr(biminLeft /M->b[imin-1])+1);
        M->b[imin-1]=sqrt(Sqr(biminLeft )+Sqr(M->b[imin-1]));
      }
    }
    for(i=imin;i<M->N-1;i++)
    {
      M->a[i]=M->a[i+1];
      M->b[i]=M->b[i+1];
    }
    M->N--;
    i=imin;
    while((i>1)&&(M->a[i-1]<M->a[i-2]))
    {
      temp=M->a[i-1];
      M->a[i-1]=M->a[i-2];
      M->a[i-2]=temp;
      temp=M->b[i-1];
      M->b[i-1]=M->b[i-2];
      M->b[i-2]=temp;
      i--;
    }
    i=imin;
    while((i<M->N-1)&&(M->a[i]>M->a[i+1]))
    {
      temp=M->a[i];
      M->a[i]=M->a[i+1];
      M->a[i+1]=temp;
      temp=M->b[i];
      M->b[i]=M->b[i+1];
      M->b[i+1]=temp;
      i++;
    }
  }
  //EndTiming("Total ReduceSqrtListOfPolesNumber",&Time);
  //PrintTiming(Time);
  //FreeTiming(&Time);
}


/*  Create a new set of poles, which has lenght N, by removing the M->N-N poles with the smallest weight
 *  This algorithm combines the pole with the smallest weight with its nearest two neighbors
 *  and repeats this untill the number of nonzero poles is N
 */
void ReduceSqrtListOfPolesMoment(ListOfPoles * M, uint32_t N)
{
  if(M->N<N) return;
  uint32_t i,imin;
  double biminLeft,biminRight,temp,MinMoment,TestMoment;
  OrderListOfPoles(M);
  while(M->N>N)
  {
    imin=0;
    MinMoment=Sqr(M->b[imin])*(M->a[imin+1]-M->a[imin]);
    TestMoment=Sqr(M->b[M->N-1])*(M->a[M->N-1]-M->a[M->N-2]);
    if(TestMoment<MinMoment)
    {
      MinMoment=TestMoment;
      imin=M->N-1;
    }
    for(i=1;i<M->N-1;i++)
    {
      TestMoment=Sqr(M->b[i])*Min(M->a[i+1]-M->a[i],M->a[i]-M->a[i-1]);
      if(TestMoment<MinMoment)
      {
        imin=i;
        MinMoment=TestMoment;
      }
    }
    if(M->b[imin]>DBL_MIN)
    {
      if(imin==0)
      {
        M->a[imin+1]=(Sqr(M->b[imin]/M->b[imin+1])*M->a[imin]+M->a[imin+1])/(Sqr(M->b[imin]/M->b[imin+1])+1);
        M->b[imin+1]=sqrt(Sqr(M->b[imin])+Sqr(M->b[imin+1]));
      }
      else if(imin==M->N-1)
      {
        M->a[imin-1]=(Sqr(M->b[imin]/M->b[imin-1])*M->a[imin]+M->a[imin-1])/(Sqr(M->b[imin]/M->b[imin-1])+1);
        M->b[imin-1]=sqrt(Sqr(M->b[imin])+Sqr(M->b[imin-1]));
      }
      else
      {
        if(fabs((M->a[imin+1]-M->a[imin-1]))<fabs(10*(DBL_MIN+DBL_EPSILON*M->a[imin])))
        {
          biminLeft=M->b[imin]*sqrt(0.5);
          biminRight=M->b[imin]*sqrt(0.5);
        }
        else if(fabs((M->a[imin+1]-M->a[imin]))<fabs(10*(DBL_MIN+DBL_EPSILON*M->a[imin])))
        {
          biminLeft=M->b[imin];
          biminRight=0;
        }
        else if(fabs((M->a[imin]-M->a[imin-1]))<fabs(10*(DBL_MIN+DBL_EPSILON*M->a[imin])))
        {
          biminLeft=0;
          biminRight=M->b[imin];
        }
        else
        {
          biminLeft =M->b[imin]*sqrt((M->a[imin+1]-M->a[imin  ])/(M->a[imin+1]-M->a[imin-1]));
          biminRight=M->b[imin]*sqrt((M->a[imin  ]-M->a[imin-1])/(M->a[imin+1]-M->a[imin-1]));
        }
        M->a[imin+1]=(Sqr(biminRight/M->b[imin+1])*M->a[imin]+M->a[imin+1])/(Sqr(biminRight/M->b[imin+1])+1);
        M->b[imin+1]=sqrt(Sqr(biminRight)+Sqr(M->b[imin+1]));
        M->a[imin-1]=(Sqr(biminLeft /M->b[imin-1])*M->a[imin]+M->a[imin-1])/(Sqr(biminLeft /M->b[imin-1])+1);
        M->b[imin-1]=sqrt(Sqr(biminLeft )+Sqr(M->b[imin-1]));
      }
    }
    for(i=imin;i<M->N-1;i++)
    {
      M->a[i]=M->a[i+1];
      M->b[i]=M->b[i+1];
    }
    M->N--;
    i=imin;
    while((i>1)&&(M->a[i-1]<M->a[i-2]))
    {
      temp=M->a[i-1];
      M->a[i-1]=M->a[i-2];
      M->a[i-2]=temp;
      temp=M->b[i-1];
      M->b[i-1]=M->b[i-2];
      M->b[i-2]=temp;
      i--;
    }
    i=imin;
    while((i<M->N-1)&&(M->a[i]>M->a[i+1]))
    {
      temp=M->a[i];
      M->a[i]=M->a[i+1];
      M->a[i+1]=temp;
      temp=M->b[i];
      M->b[i]=M->b[i+1];
      M->b[i+1]=temp;
      i++;
    }
  }
}

uint32_t ReducePolesAndersonMatrixNumberFermi(AndersonMatrixType *const A, uint32_t N, double const Ef)
{
  ListOfPoles Alp;
  if(A->N<=N) {printf("warning in ReducePolesAndersonMatrixNumber, A->N=%i and N is %i\n",A->N,N);return 0;}
  if(N<=1) {printf("Error in ReducePolesAndersonMatrixNumber N=%i\n",N);return 1;}
  Alp.a=&(A->a[2]);
  Alp.b=&(A->b[1]);
  Alp.N=A->N-1;
  ReduceSqrtListOfPolesNumberFermi(&Alp, N-1, Ef);
  A->N=Alp.N+1;
  return 0;
}

uint32_t ReducePolesAndersonMatrixNumber(AndersonMatrixType *const A, uint32_t N)
{
  ListOfPoles Alp;
  if(A->N<=N) {printf("warning in ReducePolesAndersonMatrixNumber, A->N=%i and N is %i\n",A->N,N);return 0;}
  if(N<=1) {printf("Error in ReducePolesAndersonMatrixNumber N=%i\n",N);return 1;}
  Alp.a=&(A->a[2]);
  Alp.b=&(A->b[1]);
  Alp.N=A->N-1;
  ReduceSqrtListOfPolesNumber(&Alp, N-1);
  A->N=Alp.N+1;
  return 0;
}

uint32_t ReducePolesAndersonMatrixMoment(AndersonMatrixType *const A, uint32_t N)
{
  ListOfPoles Alp;
  if(A->N<=N) {printf("warning in ReducePolesAndersonMatrixNumber, A->N=%i and N is %i\n",A->N,N);return 0;}
  if(N<=1) {printf("Error in ReducePolesAndersonMatrixNumber N=%i\n",N);return 1;}
  Alp.a=&(A->a[2]);
  Alp.b=&(A->b[1]);
  Alp.N=A->N-1;
  ReduceSqrtListOfPolesMoment(&Alp, N-1);
  A->N=Alp.N+1;
  return 0;
}

uint32_t ReducePolesAndersonMatrix(AndersonMatrixType *const A, double tMin, double Ef)
{
  uint32_t i,j;
  for(i=2;i<A->N;i++)
  {
    if(A->b[i-1]<100*DBL_MIN)
    {
      A->b[i-1]=0;
    }
    else if(A->a[i+1]<Ef)
    {
      if(A->b[i-1]<tMin)
      {
        A->a[i+1]=(Sqr(A->b[i-1])*A->a[i]+Sqr(A->b[i])*A->a[i+1])/(Sqr(A->b[i-1])+Sqr(A->b[i]));
        A->a[i]=0;
        A->b[i]=sqrt(Sqr(A->b[i-1])+Sqr(A->b[i]));
        A->b[i-1]=0;
      }
    }
  }
  for(i=A->N;i>2;i--)
  {
    if(A->b[i-1]<100*DBL_MIN)
    {
      A->b[i-1]=0;
    }
    else if(A->a[i-1]>Ef)
    {
      if(A->b[i-1]<tMin)
      {
        A->a[i-1]=(Sqr(A->b[i-1])*A->a[i]+Sqr(A->b[i-2])*A->a[i-1])/(Sqr(A->b[i-1])+Sqr(A->b[i-2]));
        A->a[i]=0;
        A->b[i-2]=sqrt(Sqr(A->b[i-1])+Sqr(A->b[i-2]));
        A->b[i-1]=0;
      }
    }
  }
  i=2;
  while(i<A->N)
  {
    if(A->b[i]==0)
    {
      for(j=i+1;j<A->N;j++)
      {
        A->b[j-1]=A->b[j];
        A->a[j]=A->a[j+1];
      }
      A->N--;
    }
    else
    {
      i++;
    }
  }
  return 0;
}
/*  
 *
 */
uint32_t ReducePolesTriDiagonalMatrix(TriDiagonalMatrixType *const A, uint32_t const NpolesMax, double const Ef)
{
  double tMin;
  AndersonMatrixType An;
  if(TriDiagonalMatrixToAndersonMatrix(*A,&An)){printf("TriDiagonalMatrixToAndersonMatrix failed in ReducePolesTriDiagonalMatrix\n");fflush(stdout);return 1;};
  tMin=A->b[1]/sqrt((double)NpolesMax)+100*DBL_MIN;
  printf("For %s,\nReduce Poles around %lf, from %i ",A->Name,Ef,A->N);
  ReducePolesAndersonMatrix(&An, tMin, Ef);
  FreeTriDiagonalMatrix(A);
  AndersonMatrixToTriDiagonalMatrix(An, A, An.N);
  printf(" to %i by makeing all hybridizations larger than %lf\n",A->N,tMin);
  return 0;
}



/* M=A-B
 *
 */
uint32_t DifferenceGreensFunctionTriDiagonalRepresentation(TriDiagonalMatrixType const A, TriDiagonalMatrixType const B, TriDiagonalMatrixType *const M)
{
  ListOfPoles Ap,Bp;
  if(TriDiagonalMatrixToListOfPoles(A, &Ap)){printf("TriDiagonalMatrixToListOfPoles failed in DifferenceGreensFunctionTriDiagonalRepresentation\n");fflush(stdout);return 1;};
  if(TriDiagonalMatrixToListOfPoles(B, &Bp)){printf("TriDiagonalMatrixToListOfPoles failed in DifferenceGreensFunctionTriDiagonalRepresentation\n");fflush(stdout);return 1;};
  if(SubtractListOfPoles(&Ap, Bp))
  {
    if((Ap.N==1)&&(Ap.b[0]<0))
    {
      printf("SubtractListOfPoles returned a singel pole with weight %24.17E in DifferenceGreensFunctionTriDiagonalRepresentation\n",Ap.b[0]);
      printf("For \"%s\" - \"%s\"\n",A.Name,B.Name);
      printf("Set weight to zero\n");
      Ap.b[0]=0;
    }
    else
    {
      printf("SubtractListOfPoles failed in DifferenceGreensFunctionTriDiagonalRepresentation\n");
      fflush(stdout);
      return 1;
    }
  }
  if(ListOfPolesToTridiagonalMatrix(Ap, M, Ap.N)){printf("ListOfPolesToTridiagonalMatrix failed in DifferenceGreensFunctionTriDiagonalRepresentation\n");fflush(stdout);return 1;};
  M->a[0]=A.a[0]-B.a[0];
  return 0;
}


/* S = G0^{-1} - G^{-1}
 *
 */
uint32_t SelfEnergyTriDiagonalRepresentation(TriDiagonalMatrixType const G0, TriDiagonalMatrixType const G, TriDiagonalMatrixType *const Sigma)
{
  TriDiagonalMatrixType IG0,IG;
  IG0.N=G0.N-1;
  IG.N=G.N-1;
  IG0.a=&(G0.a[1]);
  IG0.b=&(G0.b[1]);
  IG.a=&(G.a[1]);
  IG.b=&(G.b[1]);
  strcpy(IG0.Name,G0.Name);
  strncat(IG0.Name, "^-1 (-1)", MAX(0,(int)((int)MAXStrLength-(int)strlen(IG0.Name))));
  strcpy(IG.Name,G.Name);
  strncat(IG.Name, "^-1 (-1)", MAX(0,(int)((int)MAXStrLength-(int)strlen(IG.Name))));
  return DifferenceGreensFunctionTriDiagonalRepresentation(IG, IG0, Sigma);
}

uint32_t SelfEnergyListOfPolesFromAndersonRepresentation(AndersonMatrixType const G0, AndersonMatrixType const G, ListOfPoles *const Sigma)
{
  uint32_t i;
  ListOfPoles IG0;
  if(AndersonMatrixInverseToListOfPoles(G,Sigma)){printf("AndersonMatrixInverseToListOfPoles failed in SelfEnergyListOfPolesFromAndersonRepresentation\n");return 1;}
  if(AndersonMatrixInverseToListOfPoles(G0, &IG0)){printf("AndersonMatrixInverseToListOfPoles failed in SelfEnergyListOfPolesFromAndersonRepresentation\n");return 1;}
  if(SubtractListOfPoles(Sigma, IG0))
  {
    printf("SubtractListOfPoles failed in SelfEnergyListOfPolesFromAndersonRepresentation\n");
    for(i=0;i<Sigma->N;i++)
    {
      Sigma->a[i]=0;
      Sigma->b[i]=0;
    }
  }
  strcpy(Sigma->Name,"Self energy");
  printf("Self energy moments\n");
  PrintListOfPolesEssentials(*Sigma);
  FreeListOfPoles(&IG0);
  return 0;
}

uint32_t SelfEnergyListOfPolesFromTriDiagonalRepresentation(TriDiagonalMatrixType const G0, TriDiagonalMatrixType const G, ListOfPoles *const Sigma)
{
  ListOfPoles IG0;
  AndersonMatrixType G0And, GAnd;
  if(TriDiagonalMatrixToAndersonMatrix(G0, &G0And)){printf("TriDiagonalMatrixToAndersonMatrix failed in SelfEnergyListOfPolesFromTriDiagonalRepresentation\n");return 1;}
  if(TriDiagonalMatrixToAndersonMatrix(G, &GAnd)){printf("TriDiagonalMatrixToAndersonMatrix failed in SelfEnergyListOfPolesFromTriDiagonalRepresentation\n");return 1;}
  if(AndersonMatrixInverseToListOfPoles(GAnd,Sigma)){printf("AndersonMatrixInverseToListOfPoles failed in SelfEnergyListOfPolesFromAndersonRepresentation\n");return 1;}
  FreeAndersonMatrix(&GAnd);
  if(AndersonMatrixInverseToListOfPoles(G0And, &IG0)){printf("AndersonMatrixInverseToListOfPoles failed in SelfEnergyListOfPolesFromAndersonRepresentation\n");return 1;}
  FreeAndersonMatrix(&G0And);
  if(SubtractListOfPoles(Sigma, IG0)){printf("SubtractListOfPoles failed in SelfEnergyListOfPolesFromAndersonRepresentation\n");return 1;};
  FreeListOfPoles(&IG0);
  return 0;
}


/*  GBath(w) = 1 / ( G_0(w+mu-Sigma(w))^{-1} + Sigma(w) )
 *    G_0(w+mu-Sigma(w))                     = 1 / (w - S - a1 - sum_i b_i^2 / (w - a_{i+1} - S)
 *    G_0(w+mu-Sigma(w))^{-1}                =     (w - S - a1 - sum_i b_i^2 / (w - a_{i+1} - S)
 *    G_0(w+mu-Sigma(w))^{-1} + Sigma        =     (w     - a1 - sum_i b_i^2 / (w - a_{i+1} - S)
 *  ( G_0(w+mu-Sigma(w))^{-1} + Sigma )^{-1} = 1 / (w     - a1 - sum_i b_i^2 / (w - a_{i+1} - S)
 */
uint32_t GBathAndTriAndRepresentation(AndersonMatrixType const G0, double const mu, TriDiagonalMatrixType const Sigma, AndersonMatrixType *const G)
{
  uint32_t i,j;
  ListOfPoles bSigmaPole;
  TriDiagonalMatrixType bSigma;
  G->N=(G0.N-1)*(Sigma.N+1)+1;
  if(InitAndersonMatrix(G)){printf("InitAndersonMatrix failed in GBathAndListAndRepresentation\n");return 1;}
  bSigma.N=Sigma.N+1;
  if(InitTriDiagonalMatrix(&bSigma)){printf("InitTriDiagonalMatrix failed in GBathAndListAndRepresentation\n");return 1;}
  bSigma.a[0]=0;
  for(j=0;j<Sigma.N;j++)
  {
    bSigma.a[j+2]=Sigma.a[j+1];
    bSigma.b[j+1]=Sigma.b[j];
  }
  for(i=1;i<G0.N;i++)
  {
    bSigma.b[0]=G0.b[i];
    bSigma.a[1]=G0.a[i+1]+Sigma.a[0];
    if(TriDiagonalMatrixToListOfPoles(bSigma, &bSigmaPole)){printf("TriDiagonalMatrixToListOfPoles failed in GBathAndListAndRepresentation\n");return 1;}
    for(j=0;j<bSigmaPole.N;j++)
    {
      G->a[j+(i-1)*(Sigma.N+1)+2]=bSigmaPole.a[j];
      G->b[j+(i-1)*(Sigma.N+1)+1]=sqrt(bSigmaPole.b[j]);
    }
    for(j=bSigmaPole.N;j<bSigma.N;j++)
    {
      G->a[j+(i-1)*(Sigma.N+1)+2]=0;
      G->b[j+(i-1)*(Sigma.N+1)+1]=0;
    }
    FreeListOfPoles(&bSigmaPole);
  }
  G->a[0]=G0.a[0];
  G->b[0]=G0.b[0];
  G->a[1]=G0.a[1];
  FreeTriDiagonalMatrix(&bSigma);
  return 0;
}

/*  GBath(w) = 1 / ( G_0(w-Sigma(w))^{-1} + Sigma(w) )
 *    G_0(w-Sigma(w))                     = 1 / (w - S - a1 - b1^2 / (w - S - a2 - b2^2 / (w - S - ...
 *    G_0(w-Sigma(w))^{-1}                =     (w - S - a1 - b1^2 / (w - S - a2 - b2^2 / (w - S - ...
 *    G_0(w-Sigma(w))^{-1} + Sigma        =     (w     - a1 - b1^2 / (w - S - a2 - b2^2 / (w - S - ...
 *  ( G_0(w-Sigma(w))^{-1} + Sigma )^{-1} = 1 / (w     - a1 - b1^2 / (w - S - a2 - b2^2 / (w - S - ...
 */
uint32_t GBathTriDiagonalRepresentation(TriDiagonalMatrixType const G0, TriDiagonalMatrixType const Sigma, TriDiagonalMatrixType *const GBath, uint32_t NTriDiag)
{
  uint32_t i,j;
  MatrixType T;
  SparseMatrixType GS;
  GS.Ni=(G0.N-1)*(Sigma.N+1)+1;
  GS.Nj=(G0.N-1)*(Sigma.N+1)+1;
  GS.QComplex=0;
  GS.NBitsKey=8*((uint32_t)(log((double)GS.Ni)/(8*log(2.0)))+1);
  if(InitSparseMatrix(&GS)){printf("InitSparseMatrix failed in GBathTriDiagonalRepresentation\n");fflush(stdout);return 1;};
  GS.a0=G0.a[0];
  GS.b0=G0.b[0];
  RealSparseMatrixAddElement(&GS,0,0,G0.a[1]);
  RealSparseMatrixAddElement(&GS,0,1,G0.b[1]);
  RealSparseMatrixAddElement(&GS,1,0,G0.b[1]);
  for(i=0;i<G0.N-1;i++)
  {
    RealSparseMatrixAddElement(&GS,1+(i*(Sigma.N+1)),1+(i*(Sigma.N+1)),G0.a[i+2]-Sigma.a[0]);
    for(j=1;j<Sigma.N+1;j++)
    {
      RealSparseMatrixAddElement(&GS,j+1+(i*(Sigma.N+1)) , j+1+(i*(Sigma.N+1)) , Sigma.a[j]  );
      RealSparseMatrixAddElement(&GS,j  +(i*(Sigma.N+1)) , j+1+(i*(Sigma.N+1)) , Sigma.b[j-1]);
      RealSparseMatrixAddElement(&GS,j+1+(i*(Sigma.N+1)) , j  +(i*(Sigma.N+1)) , Sigma.b[j-1]);
    }
  }
  for(i=1;i<G0.N-1;i++)
  {
    RealSparseMatrixAddElement(&GS,1+((i-1)*(Sigma.N+1)) , 1+(i    *(Sigma.N+1)),G0.b[i+1]);
    RealSparseMatrixAddElement(&GS,1+(i    *(Sigma.N+1)) , 1+((i-1)*(Sigma.N+1)),G0.b[i+1]);
  }
  if(SparseMatrixToTriDiagonalMatrix(GS,&T,GBath,NTriDiag)){printf("SparseMatrixToTriDiagonalMatrix failed in GBathTriDiagonalRepresentation\n");fflush(stdout);return 1;};
  FreeMatrix(&T);
  FreeSparseMatrix(&GS);
  return 0;
}

/*  GBath(w) = 1 / ( G_0(w-Sigma(w))^{-1} + Sigma(w) )
 *    G_0(w-Sigma(w))                     = 1 / (w - S - a1 - sum_i b_i^2 / (w - a_{i+1} - S)
 *    G_0(w-Sigma(w))^{-1}                =     (w - S - a1 - sum_i b_i^2 / (w - a_{i+1} - S)
 *    G_0(w-Sigma(w))^{-1} + Sigma        =     (w     - a1 - sum_i b_i^2 / (w - a_{i+1} - S)
 *  ( G_0(w-Sigma(w))^{-1} + Sigma )^{-1} = 1 / (w     - a1 - sum_i b_i^2 / (w - a_{i+1} - S)
 */
uint32_t GBathAndListAndRepresentation(AndersonMatrixType const G0, double const mu, ListOfPoles const Sigma, AndersonMatrixType *const G)
{
  uint32_t i,j;
  ListOfPoles bSigmaPole;
  AndersonMatrixType bSigma;
  G->N=(G0.N-1)*(Sigma.N+1)+1;
  if(InitAndersonMatrix(G)){printf("InitAndersonMatrix failed in GBathAndListAndRepresentation\n");return 1;}
  strcpy(G->Name,"GBath");
  bSigma.N=Sigma.N+1;
  if(InitAndersonMatrix(&bSigma)){printf("InitAndersonMatrix failed in GBathAndListAndRepresentation\n");return 1;}
  bSigma.a[0]=0;
  for(i=1;i<G0.N;i++)
  {
    bSigma.b[0]=G0.b[i];
    bSigma.a[1]=G0.a[i+1]+Sigma.a0;
    for(j=0;j<Sigma.N;j++)
    {
      bSigma.a[j+2]=Sigma.a[j];
      bSigma.b[j+1]=sqrt(Sigma.b[j]);
    }
    if(AndersonMatrixToListOfPoles(bSigma, &bSigmaPole)){printf("AndersonMatrixToListOfPoles failed in GBathAndListAndRepresentation\n");return 1;}
    for(j=0;j<bSigmaPole.N;j++)
    {
      G->a[j+(i-1)*(Sigma.N+1)+2]=bSigmaPole.a[j];
      G->b[j+(i-1)*(Sigma.N+1)+1]=sqrt(bSigmaPole.b[j]);
    }
    for(j=bSigmaPole.N;j<bSigma.N;j++)
    {
      G->a[j+(i-1)*(Sigma.N+1)+2]=0;
      G->b[j+(i-1)*(Sigma.N+1)+1]=0;
    }
    FreeListOfPoles(&bSigmaPole);
  }
  G->a[0]=G0.a[0];
  G->b[0]=G0.b[0];
  G->a[1]=G0.a[1];
  G->mu=mu;
  FreeAndersonMatrix(&bSigma);
  PrintAndersonMatrixEssentials(*G);
  return 0;
}


/*  GBath(w) = 1 / ( G_0(w-Sigma(w))^{-1} + Sigma(w) )
 *    G_0(w-Sigma(w))                     = 1 / (w - S - a1 - sum_i b_i^2 / (w - a_{i+1} - S)
 *    G_0(w-Sigma(w))^{-1}                =     (w - S - a1 - sum_i b_i^2 / (w - a_{i+1} - S)
 *    G_0(w-Sigma(w))^{-1} + Sigma        =     (w     - a1 - sum_i b_i^2 / (w - a_{i+1} - S)
 *  ( G_0(w-Sigma(w))^{-1} + Sigma )^{-1} = 1 / (w     - a1 - sum_i b_i^2 / (w - a_{i+1} - S)
 */
uint32_t GBathAndTriAndRepresentationParallel(AndersonMatrixType const G0, TriDiagonalMatrixType const Sigma, AndersonMatrixType *const G)
{
  int i,j,failed=0;
  G->N=(G0.N-1)*(Sigma.N+1)+1;
  if(InitAndersonMatrix(G)){printf("InitAndersonMatrix (G->N=%i) failed in GBathAndListAndRepresentation\n",G->N);fflush(stdout);return 1;}
  G->a[0]=G0.a[0];
  G->b[0]=G0.b[0];
  G->a[1]=G0.a[1];
  #pragma omp parallel private(i,j)
  {
    ListOfPoles bSigmaPole;
    TriDiagonalMatrixType bSigma;
    bSigma.N=Sigma.N+1;
    if(InitTriDiagonalMatrix(&bSigma))
    {
      printf("InitTriDiagonalMatrix (bSigma.N=%i) failed in GBathAndListAndRepresentation in parallel region\n",bSigma.N);
      fflush(stdout);
      failed=1;
    }
    if(!failed)
    {
      bSigma.a[0]=0;
      for(j=0;j<Sigma.N;j++)
      {
        bSigma.a[j+2]=Sigma.a[j+1];
        bSigma.b[j+1]=Sigma.b[j];
      }
      #pragma omp for
      for(i=1;i<G0.N;i++)
      {
        bSigma.b[0]=G0.b[i];
        bSigma.a[1]=G0.a[i+1]+Sigma.a[0];
        if(TriDiagonalMatrixToListOfPoles(bSigma, &bSigmaPole))
        {
          printf("TriDiagonalMatrixToListOfPoles (bSigma.N=%i) failed in GBathAndListAndRepresentation in parallel region\n",bSigma.N);
          fflush(stdout);
          failed=1;
        }
        if(!failed)
        {
          for(j=0;j<bSigmaPole.N;j++)
          {
            G->a[j+(i-1)*(Sigma.N+1)+2]=bSigmaPole.a[j];
            G->b[j+(i-1)*(Sigma.N+1)+1]=sqrt(bSigmaPole.b[j]);
          }
          for(j=bSigmaPole.N;j<bSigma.N;j++)
          {
            G->a[j+(i-1)*(Sigma.N+1)+2]=0;
            G->b[j+(i-1)*(Sigma.N+1)+1]=0;
          }
          FreeListOfPoles(&bSigmaPole);
        }
      }
      FreeTriDiagonalMatrix(&bSigma);
    }
  }
  if(failed) return 1;
  return 0;
}

/*  GBath(w) = 1 / ( G_0(w-Sigma(w))^{-1} + Sigma(w) )
 *    G_0(w-Sigma(w))                     = 1 / (w - S - a1 - sum_i b_i^2 / (w - a_{i+1} - S)
 *    G_0(w-Sigma(w))^{-1}                =     (w - S - a1 - sum_i b_i^2 / (w - a_{i+1} - S)
 *    G_0(w-Sigma(w))^{-1} + Sigma        =     (w     - a1 - sum_i b_i^2 / (w - a_{i+1} - S)
 *  ( G_0(w-Sigma(w))^{-1} + Sigma )^{-1} = 1 / (w     - a1 - sum_i b_i^2 / (w - a_{i+1} - S)
 */
uint32_t OldGBathAndTriAndRepresentation(AndersonMatrixType const G0, TriDiagonalMatrixType const Sigma, AndersonMatrixType *const G)
{
  uint32_t i,j;
  ListOfPoles bSigmaPole;
  TriDiagonalMatrixType bSigma;
  G->N=(G0.N-1)*(Sigma.N+1)+1;
  if(InitAndersonMatrix(G)){printf("InitAndersonMatrix failed in GBathAndListAndRepresentation\n");return 1;}
  bSigma.N=Sigma.N+1;
  if(InitTriDiagonalMatrix(&bSigma)){printf("InitTriDiagonalMatrix failed in GBathAndListAndRepresentation\n");return 1;}
  bSigma.a[0]=0;
  for(j=0;j<Sigma.N;j++)
  {
    bSigma.a[j+2]=Sigma.a[j+1];
    bSigma.b[j+1]=Sigma.b[j];
  }
  for(i=1;i<G0.N;i++)
  {
    bSigma.b[0]=G0.b[i];
    bSigma.a[1]=G0.a[i+1]+Sigma.a[0]; 
    if(TriDiagonalMatrixToListOfPoles(bSigma, &bSigmaPole)){printf("TriDiagonalMatrixToListOfPoles failed in GBathAndListAndRepresentation\n");return 1;}
    for(j=0;j<bSigmaPole.N;j++)
    {
      G->a[j+(i-1)*(Sigma.N+1)+2]=bSigmaPole.a[j];
      G->b[j+(i-1)*(Sigma.N+1)+1]=sqrt(bSigmaPole.b[j]);
    }
    for(j=bSigmaPole.N;j<bSigma.N;j++)
    {
      G->a[j+(i-1)*(Sigma.N+1)+2]=0;
      G->b[j+(i-1)*(Sigma.N+1)+1]=0;
    }
    FreeListOfPoles(&bSigmaPole);
  }
  G->a[0]=G0.a[0];
  G->b[0]=G0.b[0];
  G->a[1]=G0.a[1];
  FreeTriDiagonalMatrix(&bSigma);
  return 0;
}

/*  GbathNew = (Gcl^{-1}-Sigma^{-1}) ^{-1}
 *
 *  Assume Gcl.a[0]=0;
 *         Gcl.b[0]=1;
 */
uint32_t OldGBathNewTriDiagonalRepresentation(TriDiagonalMatrixType const Gcl, TriDiagonalMatrixType const Sigma, TriDiagonalMatrixType * const GBathNew)
{
  uint32_t i;
  TriDiagonalMatrixType GclInv,SigmaInv,GBathNewInv;
  GclInv.N=Gcl.N-1;
  GclInv.a=&(Gcl.a[1]);
  GclInv.b=&(Gcl.b[1]);
  strcpy(GclInv.Name,Gcl.Name);
  strncat(GclInv.Name, "^-1 (-1)", MAX(0,(int)((int)MAXStrLength-(int)strlen(Gcl.Name))));
  SigmaInv.N=Sigma.N+1;
  InitTriDiagonalMatrix(&SigmaInv);
  SigmaInv.a[0]=0;
  SigmaInv.b[0]=1;
  for(i=0;i<Sigma.N;i++)
  {
    SigmaInv.a[i+1]=Sigma.a[i];
    SigmaInv.b[i+1]=Sigma.b[i];
  }
  SigmaInv.a[Sigma.N+1]=Sigma.a[Sigma.N];
  
  DifferenceGreensFunctionTriDiagonalRepresentation(SigmaInv, GclInv, &GBathNewInv);
  GBathNew->N=GBathNewInv.N+1;
  InitTriDiagonalMatrix(GBathNew);
  GBathNew->a[0]=0;
  GBathNew->b[0]=1;
  for(i=1;i<GBathNew->N;i++)
  {
    GBathNew->a[i]=GBathNewInv.a[i-1];
    GBathNew->b[i]=GBathNewInv.b[i-1];
  }
  GBathNew->a[GBathNew->N]=GBathNewInv.a[GBathNew->N-1];
  FreeTriDiagonalMatrix(&GBathNewInv);
  FreeTriDiagonalMatrix(&SigmaInv);
  return 0;
}


/*  GbathNew = (Sigma + Glocal^{-1}) ^{-1}
 *
 *  Assume Glocal.a[0]=0;
 *         Glocal.b[0]=1;
 */
uint32_t VeryOldGBathNewTriDiagonalRepresentation(TriDiagonalMatrixType const Glocal, TriDiagonalMatrixType const Sigma, TriDiagonalMatrixType * const GBathNew)
{
  uint32_t i;
  TriDiagonalMatrixType GlocInv,GBathNewInv;
  GlocInv.N=Glocal.N-1;
  GlocInv.a=&(Glocal.a[1]);
  GlocInv.b=&(Glocal.b[1]);
  strcpy(GlocInv.Name,Glocal.Name);
  strncat(GlocInv.Name, "^-1 (-1)", MAX(0,(int)((int)MAXStrLength-(int)strlen(Glocal.Name))));
  DifferenceGreensFunctionTriDiagonalRepresentation(GlocInv, Sigma, &GBathNewInv);
  GBathNew->N=GBathNewInv.N+1;
  InitTriDiagonalMatrix(GBathNew);
  GBathNew->a[0]=0;
  GBathNew->b[0]=1;
  for(i=1;i<GBathNew->N;i++)
  {
    GBathNew->a[i]=GBathNewInv.a[i-1];
    GBathNew->b[i]=GBathNewInv.b[i-1];
  }
  GBathNew->a[GBathNew->N]=GBathNewInv.a[GBathNew->N-1];
  FreeTriDiagonalMatrix(&GBathNewInv);
  return 0;
}




/*  This routine returns a tri-diagonal matrix representing a greensfunction as a continued fraction in Gnew
 *  
 */
uint32_t TriDiagonalMatrixDysonScalar(TriDiagonalMatrixType const G, TriDiagonalMatrixType *const Gnew, MatrixType const Rho, uint32_t const NRho, TriDiagonalMatrixType *const Gp)
{
  SparseMatrixType M;
  MatrixType T;
  uint32_t i,j,k,l,NTriDiag,jstuck;
  M.Ni=G.N;
  for(i=0;i<NRho;i++)
  {
    M.Ni+=Gp[i].N*G.N;
  }
  M.Nj=M.Ni;
  M.QComplex=0;
  M.NBitsKey=16;
  InitSparseMatrix(&M);
  M.a0=G.a[0];
  M.b0=G.b[0];
  j=0;
  for(i=0;i<G.N-1;i++)
  {
    RealSparseMatrixAddElement(&M, j, j, G.a[i+1]);
    RealSparseMatrixAddElement(&M, j, j+1, G.b[i+1]);
    RealSparseMatrixAddElement(&M, j+1, j, G.b[i+1]);
    j++;
  }
  RealSparseMatrixAddElement(&M, j, j, G.a[G.N]);
  j++;
  jstuck=j;
  for(k=0;k<NRho;k++)
  {
    for(l=0;l<G.N;l++)
    {
      RealSparseMatrixAddElement(&M, l, j, Gp[k].b[0]*sqrt(fabs(Rho.VRe[l][k])));
      RealSparseMatrixAddElement(&M, j, l, Gp[k].b[0]*sqrt(fabs(Rho.VRe[l][k])));
      for(i=0;i<Gp[k].N-1;i++)
      {
        RealSparseMatrixAddElement(&M, j, j, Gp[k].a[i+1]);
        RealSparseMatrixAddElement(&M, j, j+1, Gp[k].b[i+1]);
        RealSparseMatrixAddElement(&M, j+1, j, Gp[k].b[i+1]);
        j++;
      }
      RealSparseMatrixAddElement(&M, j, j, Gp[k].a[Gp[k].N]);
      j++;
    }
  }
  /*
  printf("\n\n\n\n\n\n\n\n ******* \n\n\n\n\n");
  PrintSparseMatrix(M);
  */
  
  NTriDiag=G.N;
  for(i=0;i<NRho;i++)
  {
    NTriDiag=Max(NTriDiag,Gp[i].N);
  }
  if(SparseMatrixToTriDiagonalMatrix(M, &T, Gnew, NTriDiag)){printf("SparseMatrixToTriDiagonalMatrix failed in TriDiagonalMatrixDysonScalar\n");fflush(stdout);return 1;}
  FreeMatrix(&T);
  FreeSparseMatrix(&M);
  return 0;
}


/*  We start from the following Hamiltonian
 *
 *     ( O       V1      0       0       0         V2      0       0       0         ...       V3      0       0       0   )
 *     ( V1^+    H1+g1a1 g1b1    0       0         0       0       0       0         ...       0       0       0       0   )
 *     ( 0       g1b1    H1+g1a2 g1b2    0         0       0       0       0         ...       0       0       0       0   )
 *     ( 0       0       g1b2    H1+g1a3 ...       0       0       0       0         ...       0       0       0       0   )
 *     ( 0       0       0       ...     ...       0       0       0       0         ...       0       0       0       0   )
 *
 *     ( V2^+    0       0       0       0         H2+g2a1 g2b1    0       0         ...       0       0       0       0   )
 *     ( 0       0       0       0       0         g2b1    H2+g2a2 g2b2    0         ...       0       0       0       0   )
 *     ( 0       0       0       0       0         0       g2b2    H2+g2a3 ...       ...       0       0       0       0   )
 *     ( 0       0       0       0       0         0       0       ...     ...       ...       0       0       0       0   )
 *
 *     ( ...     ...     ...     ...     ...       ...     ...     ...     ...       ...       ...     ...     ...     ... )
 *
 *     ( V3^+    0       0       0       0         0       0       0       0         ...       H3+g3a1 g3b1    0       0   )
 *     ( 0       0       0       0       0         0       0       0       0         ...       g3b1    H3+g3a2 g3b2    0   )
 *     ( 0       0       0       0       0         0       0       0       0         ...       0       g3b2    H3+g3a3 ... )
 *     ( 0       0       0       0       0         0       0       0       0         ...       0       0       ...     ... )
 *
 *  and want to create a tridiagonal matrix on a Krylov basis starting from the first basis function
 *
 *  The matrices included are of special form:
 *  O     - tridiagonal
 *  Vi    - a dense matrix
 *  Hi    - diagonal matrix
 *  gi    - tridiagonal matrix, representing the coupled bath.
 *
 *  gnew  - resulting tridiagonal matrix of length gnew.N (set on entry)
 *
 *  N counts the number of perturbing Operators V
 *
 *  H.Ni Equal to N
 *  H.Nj is Max[V[i].Nj] for all i or larger
 *
 */
/*  This routine needs to be replaced by a standard lanczos from a library!!!!
 *
 */
uint32_t TridiagonalizeImpurityCoupledToBathDyson(TriDiagonalMatrixType const O, MatrixType const *V, MatrixType const H, TriDiagonalMatrixType const *g, uint32_t N, TriDiagonalMatrixType * const gnew)
{
  SparseMatrixType Htot;
  MatrixType T;
  uint32_t i,j,k,l;
  Htot.Ni=O.N;
  for(i=0;i<N;i++)
  {
    Htot.Ni+=g[i].N * V[i].Nj;
  }
  Htot.Nj=Htot.Ni;
  if(Htot.Ni<16000)
  {
    Htot.NBitsKey=16;
  }
  else
  {
    Htot.NBitsKey=24;
  }
  Htot.QComplex = 0;
  Htot.a0=O.a[0];
  Htot.b0=O.b[0];
  for(i=0;i<N;i++)
  {
    if(V[i].QComplex)Htot.QComplex=1;
  }
  if(InitSparseMatrix(&Htot)){printf("Could not initiallize sparse matrix in TridiagonalizeImpurityCoupledToBathDyson\n");return 1;}
  SparseMatrixAddTridiagonalMatrix(&Htot, 0, 0, O);
  j=O.N;
  for(i=0;i<N;i++)
  {
    SparseMatrixAddMatrixTransposed(&Htot, 0, j, V[i], DBL_EPSILON);
    SparseMatrixAddMatrix(&Htot, j, 0, V[i], DBL_EPSILON);
    j+=g[i].N * V[i].Nj;    
  }
  l=O.N;
  if(Htot.QComplex)
  {
    for(i=0;i<N;i++)
    {
      for(j=0;j<g[i].N;j++)
      {
        for(k=0;k<V[i].Nj;k++)
        {
          ComplexSparseMatrixAddElement(&Htot, l, l, H.VRe[i][k]+g[i].a[j+1],0);
          if(j!=g[i].N-1)
          {
            ComplexSparseMatrixAddElement(&Htot, l, l+V[i].Nj, g[i].b[j+1],0);
            ComplexSparseMatrixAddElement(&Htot, l+V[i].Nj, l, g[i].b[j+1],0);
          }
          l++;
        }
      }
    }
  }
  else
  {
    for(i=0;i<N;i++)
    {
      for(j=0;j<g[i].N;j++)
      {
        for(k=0;k<V[i].Nj;k++)
        {
          RealSparseMatrixAddElement(&Htot, l, l, H.VRe[i][k]+g[i].a[j+1]);
          if(j!=g[i].N-1)
          {
            RealSparseMatrixAddElement(&Htot, l, l+V[i].Nj, g[i].b[j+1]);
            RealSparseMatrixAddElement(&Htot, l+V[i].Nj, l, g[i].b[j+1]); 
          }
          l++;
        }
      }
    }
  }
  SparseMatrixToTriDiagonalMatrix(Htot, &T, gnew, gnew->N);
  FreeMatrix(&T);
  FreeSparseMatrix(&Htot);
  return 0;
}


/*  We start from the following Hamiltonian
 *
 *     ( O       V1      0       0       0         V2      0       0       0         ...       V3      0       0       0         0             0             0             0       0             0             0             0       0             0             0             0     0   )
 *
 *     ( V1^+    H1+g1a1 g1b1    0       0         0       0       0       0         ...       0       0       0       0         V12           0             0             0       0             0             0             0       0             0             0             0     0   )
 *     ( 0       g1b1    H1+g1a2 g1b2    0         0       0       0       0         ...       0       0       0       0         0             V12           0             0       0             0             0             0       0             0             0             0     0   )
 *     ( 0       0       g1b2    H1+g1a3 ...       0       0       0       0         ...       0       0       0       0         0             0             V12           0       0             0             0             0       0             0             0             0     0   )
 *     ( 0       0       0       ...     ...       0       0       0       0         ...       0       0       0       0         0             0             0             ...     0             0             0             0       0             0             0             0     0   )
 *
 *     ( V2^+    0       0       0       0         H2+g2a1 g2b1    0       0         ...       0       0       0       0         V21           0             0             0       0             0             0             0       0             0             0             0     0   )
 *     ( 0       0       0       0       0         g2b1    H2+g2a2 g2b2    0         ...       0       0       0       0         0             0             0             0       V21           0             0             0       0             0             0             0     0   )
 *     ( 0       0       0       0       0         0       g2b2    H2+g2a3 ...       ...       0       0       0       0         0             0             0             0       0             0             0             0       V21           0             0             0     0   )
 *     ( 0       0       0       0       0         0       0       ...     ...       ...       0       0       0       0         0             0             0             0       0             0             0             0       0             0             0             0     ... )
 *
 *     ( ...     ...     ...     ...     ...       ...     ...     ...     ...       ...       ...     ...     ...     ...       ...           ...           ...           ...     ...           ...           ...           ...     ...           ...           ...           ...   ... )
 *
 *     ( V3^+    0       0       0       0         0       0       0       0         ...       H3+g3a1 g3b1    0       0         0             0             0             0       0             0             0             0       0             0             0             0     0   )
 *     ( 0       0       0       0       0         0       0       0       0         ...       g3b1    H3+g3a2 g3b2    0         0             0             0             0       0             0             0             0       0             0             0             0     0   )
 *     ( 0       0       0       0       0         0       0       0       0         ...       0       g3b2    H3+g3a3 ...       0             0             0             0       0             0             0             0       0             0             0             0     0   )
 *     ( 0       0       0       0       0         0       0       0       0         ...       0       0       ...     ...       0             0             0             0       0             0             0             0       0             0             0             0     0   )
 *
 *     ( 0       V12^+   0       0       0         V21^+   0       0       0         ...       0       0       0       0         H12+g1a1+g2a1 g1b1          0             0       g2b1          0             0             0       0             0             0             0     0   )
 *     ( 0       0       V12^+   0       0         0       0       0       0         ...       0       0       0       0         g1b1          H12+g1a2+g2a1 g1b2          0       0             g2b1          0             0       0             0             0             0     0   )
 *     ( 0       0       0       V12^+   0         0       0       0       0         ...       0       0       0       0         0             g1b2          H12+g1a3+g2a1 ...     0             0             g2b1          0       0             0             0             0     0   )
 *     ( 0       0       0       0       ...       0       0       0       0         ...       0       0       0       0         0             0             ...           ...     0             0             0             ...     0             0             0             0     0   )
 *
 *     ( 0       0       0       0       0         0       0       0       0         ...       0       0       0       0         g2b1          0             0             0       H12+g1a1+g2a2 g1b1          0             0       g2b2          0             0             0     0   )
 *     ( 0       0       0       0       0         0       V21^+   0       0         ...       0       0       0       0         0             g2b1          0             0       g1b1          H12+g1a2+g2a2 g1b2          0       0             g2b2          0             0     0   )
 *     ( 0       0       0       0       0         0       0       0       0         ...       0       0       0       0         0             0             g2b1          0       0             g1b2          H12+g1a3+g2a2 ...     0             0             g2b2          0     0   )
 *     ( 0       0       0       0       0         0       0       0       0         ...       0       0       0       0         0             0             0             ...     0             0             ...           ...     0             0             0             ...   0   )
 *
 *     ( 0       0       0       0       0         0       0       0       0         ...       0       0       0       0         0             0             0             0       g2b2          0             0             0       H12+g1a1+g2a3 g1b1          0             0     ... )
 *     ( 0       0       0       0       0         0       0       0       0         ...       0       0       0       0         0             0             0             0       0             g2b2          0             0       g1b1          H12+g1a2+g2a3 g1b2          0     ... )
 *     ( 0       0       0       0       0         0       0       V21^+   0         ...       0       0       0       0         0             0             0             0       0             0             g2b2          0       0             g1b2          H12+g1a3+g2a3 ...   ... )
 *     ( 0       0       0       0       0         0       0       0       0         ...       0       0       0       0         0             0             0             0       0             0             0             ...     0             0             ...           ...   ... )
 *
 *     ( 0       0       0       0       0         0       0       0       ...       ...       0       0       0       0         0             0             0             0       0             0             0             0       ...           ...           ...           ...   ... )
 *
 *  and want to create a tridiagonal matrix on a Krylov basis starting from the first basis function
 *
 *  The matrices included are of special form:
 *  O     - tridiagonal
 *  Vi    - a dense matrix
 *  Vii   - a dense matrix
 *  Hi    - diagonal matrix
 *  Hii   - diagonal matrix
 *  gi    - tridiagonal matrix, representing the coupled bath.
 *
 *  gnew  - resulting tridiagonal matrix of length gnew.N (set on entry)
 *
 *  N counts the number of perturbing Operators V
 *
 *  Hi.Ni Equal to N
 *  Hi.Nj is Max[V[i].Nj] for all i or larger
 *
 *  Hii.Ni Equal to (N^2-N)/2
 *  Hii.Nj is Max[V[ij].Nj] for all i,j or larger
 *
 *  Vi,Vii will be labeled V1 and V2 below.
 */
uint32_t TridiagonalizeImpurityCoupledToBathDysonSecondOrderA(TriDiagonalMatrixType const O, MatrixType const *V1, MatrixType const *V2, MatrixType const H1, MatrixType const H2, TriDiagonalMatrixType const *g, uint32_t N, TriDiagonalMatrixType * const gnew)
{
  SparseMatrixType Htot;
  MatrixType T;
  uint32_t i,ii,j,k,l,m,size1,size2,indexi,indexii,indexj;
  
  
  printf("Start of TridiagonalizeImpurityCoupledToBathDysonSecondOrderA\n");
//  for(i=0;i<N;i++)
//  {
//    PrintMatrix(V1[i]);
//  }
//  for(i=0;i<(N*N-N);i++)
//  {
//    PrintMatrix(V2[i]);
//  }
//  PrintMatrix(H1);
//  PrintMatrix(H2);
  
  Htot.Ni=O.N;
  /* Determine dimension of Htot */
  j=0;
  for(i=0;i<N;i++)
  {
    Htot.Ni+=g[i].N * V1[i].Nj;
    for(ii=0;ii<i;ii++)
    {
      Htot.Ni+=g[i].N*g[ii].N*V2[2*j].Nj;
      j++;
    }
  }
  Htot.Nj=Htot.Ni;
  if(Htot.Ni<16000)
  {
    Htot.NBitsKey=16;
  }
  else
  {
    Htot.NBitsKey=24;
  }
  Htot.QComplex = 0;
  Htot.a0=O.a[0];
  Htot.b0=O.b[0];
  j=0;
  for(i=0;i<N;i++)
  {
    if(V1[i].QComplex)Htot.QComplex=1;
    for(ii=0;ii<i;ii++)
    {
      if(V2[2*j  ].QComplex)Htot.QComplex=1;
      if(V2[2*j+1].QComplex)Htot.QComplex=1;
      j++;
    }
  }
  if(InitSparseMatrix(&Htot)){printf("Could not initiallize sparse matrix in TridiagonalizeImpurityCoupledToBathDyson\n");return 1;}
  /* Add matrix O to Htot */
  SparseMatrixAddTridiagonalMatrix(&Htot, 0, 0, O);
  /* Add matrices V1 to Htot */
  j=O.N;
  for(i=0;i<N;i++)
  {
    SparseMatrixAddMatrixTransposed(&Htot, 0, j, V1[i], DBL_EPSILON);
    SparseMatrixAddMatrix(&Htot, j, 0, V1[i], DBL_EPSILON);
    j+=g[i].N * V1[i].Nj;
  }
  size1=j;
  /* Add matrices H1 and g to Htot */
  l=O.N;
  if(Htot.QComplex)
  {
    for(i=0;i<N;i++)
    {
      for(j=0;j<g[i].N;j++)
      {
        for(k=0;k<V1[i].Nj;k++)
        {
          ComplexSparseMatrixAddElement(&Htot, l, l, H1.VRe[i][k]+g[i].a[j+1],0);
          if(j!=g[i].N-1)
          {
            ComplexSparseMatrixAddElement(&Htot, l, l+V1[i].Nj, g[i].b[j+1],0);
            ComplexSparseMatrixAddElement(&Htot, l+V1[i].Nj, l, g[i].b[j+1],0);
          }
          l++;
        }
      }
    }
  }
  else
  {
    for(i=0;i<N;i++)
    {
      for(j=0;j<g[i].N;j++)
      {
        for(k=0;k<V1[i].Nj;k++)
        {
          RealSparseMatrixAddElement(&Htot, l, l, H1.VRe[i][k]+g[i].a[j+1]);
          if(j!=g[i].N-1)
          {
            RealSparseMatrixAddElement(&Htot, l, l+V1[i].Nj, g[i].b[j+1]);
            RealSparseMatrixAddElement(&Htot, l+V1[i].Nj, l, g[i].b[j+1]);
          }
          l++;
        }
      }
    }
  }
  /* Add matrices V2 to Htot */
  j=0;
  indexi=O.N;
  indexj=0;
  for(i=0;i<N;i++)
  {
    indexii=O.N;
    for(ii=0;ii<i;ii++)
    {
      size2=0;
      for(k=0;k<g[i].N;k++)
      {
        SparseMatrixAddMatrixTransposed(&Htot, indexi+k*V2[2*j].Ni, size1+indexj+k*V2[2*j].Nj, V2[2*j], DBL_EPSILON);
        SparseMatrixAddMatrix(&Htot, size1+indexj+k*V2[2*j].Nj, indexi+k*V2[2*j].Ni, V2[2*j], DBL_EPSILON);
        size2+=V2[2*j].Ni;
      }
      for(k=0;k<g[i].N;k++)
      {
        SparseMatrixAddMatrixTransposed(&Htot, indexii+k*V2[2*j+1].Ni, size1+indexj+k*size2, V2[2*j+1], DBL_EPSILON);
        SparseMatrixAddMatrix (&Htot, size1+indexj+k*size2, indexii+k*V2[2*j+1].Ni, V2[2*j+1], DBL_EPSILON);
      }
      indexii+=V1[ii].Nj*g[ii].N;
      indexj+=V2[2*j].Nj*g[i].N*g[ii].N;
      j++;
    }
    indexi+=V1[i].Nj*g[i].N;
  }
  /* Add matrices H2 and g to Htot */
  j=0;
  indexj=size1;
  if(Htot.QComplex)
  {
    
  }
  else
  {
    for(i=0;i<N;i++)
    {
      for(ii=0;ii<i;ii++)
      {
        for(l=0;l<g[ii].N;l++)
        {
          for(m=0;m<g[i].N;m++)
          {
            for(k=0;k<V2[j].Nj;k++)
            {
              RealSparseMatrixAddElement(&Htot, indexj, indexj, H2.VRe[j][k]+g[i].a[m+1]+g[ii].a[l+1]);
              if(m!=g[i].N-1)
              {
                RealSparseMatrixAddElement(&Htot, indexj, indexj+V2[j].Nj, g[i].b[m+1]);
                RealSparseMatrixAddElement(&Htot, indexj+V2[j].Nj, indexj, g[i].b[m+1]);
              }
              if(l!=g[ii].N-1)
              {
                RealSparseMatrixAddElement(&Htot, indexj, indexj+V2[j].Nj*g[i].N, g[ii].b[l+1]);
                RealSparseMatrixAddElement(&Htot, indexj+V2[j].Nj*g[i].N, indexj, g[ii].b[l+1]);
              }
              
              indexj++;
            }
          }
        }
        j++;
      }
    }
  }
  
//  PrintSparseMatrix(Htot);
//  FilePrintSparseMatrixMathematica(Htot, "Htot", "w");
  
  /* tridiagonalize Htot */
  SparseMatrixToTriDiagonalMatrix(Htot, &T, gnew, gnew->N);
  FreeMatrix(&T);
  FreeSparseMatrix(&Htot);
  return 0;
}




/*  We start from the following Hamiltonian
 *
 *     ( H_1     H_12^1         H_12^2          ...   N_12^n        )
 *     ( H_12^1  H_1^2 + H_2^1  0               0     0             )
 *     ( H_12^2  0              H_1^2 + H_2^2   0     0             )
 *     ( ...     0              0               ...   0             )
 *     ( H_12^n  0              0               0     H_1^n + H_2^n )
 *
 *     H_1 is provided on a finite basis in tridiagonal form (G1)
 *     on the same basis the matrix H_12^i H_12^i is given   (M12)
 *     The one particle green function for H2 with the coupling H_12^i is given in tridiagonal form (g2^i)
 *     The expectation value of H_12^i H_1 H_12^i devided by the expectation value of H_12^i H_12^i is given (Emu^i)
 *
 *     Calculated is a tridiaonal matrix representing { 1/(w-H) }_[1,1]
 */
uint32_t TridDiagonalMatrixDysonToMuchSimplified(TriDiagonalMatrixType const G1, MatrixType const *M12, TriDiagonalMatrixType const *g2, MatrixType const *Emu, int32_t const Ni, TriDiagonalMatrixType * const Gnew)
{
  return 0;
}




/*
 *  assume f.a[0]=0
 *         g.a[0]=0;
 *         f.b[0]=1;
 *         g.b[0]=1;
 */
uint32_t GreenFunctionConvoluteTriDiagonalRepresentation(TriDiagonalMatrixType *const fg, TriDiagonalMatrixType const f, TriDiagonalMatrixType const g, uint32_t const NTriDiag)
{
  uint32_t i,j,k;
  SparseMatrixType fgSparse;
  MatrixType T;
  fgSparse.NBitsKey=16;
  fgSparse.QComplex=0;
  fgSparse.Ni=f.N*g.N;
  fgSparse.Nj=f.N*g.N;
  if(InitSparseMatrix(&fgSparse)){printf("InitSparseMatrix failed in GreenFunctionConvoluteTriDiagonalRepresentation\n");return 1;};
  fgSparse.a0=0;
  fgSparse.b0=1;
  k=0;
  for(i=1;i<f.N;i++)
  {
    for(j=1;j<g.N;j++)
    {
      RealSparseMatrixAddElement(&fgSparse, k    , k    , f.a[i]+g.a[j]);
      RealSparseMatrixAddElement(&fgSparse, k    , k+1  , g.b[j]);
      RealSparseMatrixAddElement(&fgSparse, k+1  , k    , g.b[j]);
      RealSparseMatrixAddElement(&fgSparse, k    , k+g.N, f.b[i]);
      RealSparseMatrixAddElement(&fgSparse, k+g.N, k    , f.b[i]);
      k++;
    }
    RealSparseMatrixAddElement(&fgSparse, k    , k    , f.a[i]+g.a[g.N]);
    RealSparseMatrixAddElement(&fgSparse, k    , k+g.N, f.b[i]);
    RealSparseMatrixAddElement(&fgSparse, k+g.N, k    , f.b[i]);
    k++;
  }
  for(j=1;j<g.N;j++)
  {
    RealSparseMatrixAddElement(&fgSparse, k    , k    , f.a[f.N]+g.a[j]);
    RealSparseMatrixAddElement(&fgSparse, k    , k+1  , g.b[j]);
    RealSparseMatrixAddElement(&fgSparse, k+1  , k    , g.b[j]);
    k++;
  }
  RealSparseMatrixAddElement(&fgSparse, k    , k    , f.a[f.N]+g.a[g.N]);
  if(SparseMatrixToTriDiagonalMatrix(fgSparse, &T, fg, NTriDiag)){printf("SparseMatrixToTriDiagonalMatrix failed in GreenFunctionConvoluteTriDiagonalRepresentation\n");return 1;};
  FreeMatrix(&T);
  FreeSparseMatrix(&fgSparse);
  return 0;
}



