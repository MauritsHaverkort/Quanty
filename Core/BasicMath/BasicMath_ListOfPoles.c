uint32_t InitListOfPoles(ListOfPoles *M)
{
  strcpy(M->Name,"List of Poles");
  if(M->N==0) {printf("Initialized a ListOfPoles with size zero.\n");fflush(stdout); return 2;};
  M->a=malloc((M->N)*sizeof(double));
  if(M->a==NULL) {printf("malloc for M->a failed in ListOfPoles\n");fflush(stdout);return 1;};
  M->b=malloc((M->N)*sizeof(double));
  if(M->b==NULL)
  {
    free(M->a);
    printf("malloc for M->b failed in ListOfPoles\n");
    fflush(stdout);
    return 1;
  }
  M->a0=0;
  return 0;
}

void FreeListOfPoles(ListOfPoles *M)
{
  if(M->N!=0)
  {
    free(M->a);
    free(M->b); 
  }
  M->a0=0;
}

void ScaleListOfPoles(ListOfPoles *const A, double const scale)
{
  uint32_t i;
  for(i=0;i<A->N;i++)
  {
    A->b[i]*=scale;
  }
}

uint32_t AddListOfPoles(ListOfPoles *const A, ListOfPoles B)
{
  double *temp;
  uint32_t i;
  temp=realloc(A->a,(A->N+B.N)*sizeof(double));
  if(temp==NULL){printf("realloc of A->a failed in AddListOfPoles\n");fflush(stdout);return 1;}
  A->a=temp;
  temp=realloc(A->b,(A->N+B.N)*sizeof(double));
  if(temp==NULL){printf("realloc of A->b failed in AddListOfPoles\n");fflush(stdout);return 1;}
  A->b=temp;
  for(i=0;i<B.N;i++)
  {
    A->a[i+A->N]=B.a[i];
    A->b[i+A->N]=B.b[i];
  }
  A->a0+=B.a0;
  A->N+=B.N;
  return 0;
}

uint32_t SymmetrizeListOfPoles(ListOfPoles *const A)
{
  double *temp;
  uint32_t i;
  temp=realloc(A->a,(2*A->N)*sizeof(double));
  if(temp==NULL){printf("realloc of A->a failed in AddListOfPoles\n");fflush(stdout);return 1;}
  A->a=temp;
  temp=realloc(A->b,(2*A->N)*sizeof(double));
  if(temp==NULL){printf("realloc of A->b failed in AddListOfPoles\n");fflush(stdout);return 1;}
  A->b=temp;
  for(i=0;i<A->N;i++)
  {
    A->a[i+A->N]=-A->a[i];
    A->b[i+A->N]=A->b[i];
  }
  A->N*=2;
  ScaleListOfPoles(A,0.5);
  return 0;
}

void OrderListOfPoles(ListOfPoles *const M)
{
  uint32_t i,j;
  double temp;
  for(i=1;i<M->N;i++)
  {
    for(j=i;j>0;j--)
    {
      if(M->a[j]<M->a[j-1])
      {
        temp=M->a[j];
        M->a[j]=M->a[j-1];
        M->a[j-1]=temp;
        temp=M->b[j];
        M->b[j]=M->b[j-1];
        M->b[j-1]=temp;
      }
      else
      {
        break;
      }
    }
  }
}

/*  Create a new set of poles, wich all have a minimal weight of bmin
 *
 */
void ReduceListOfPolesWeigth(ListOfPoles *M, double bmin)
{
  printf("ReduceListOfPolesWeigth Still to do/n");
}

/*  Create a new set of poles, which has lenght N, by removing the M->N-N poles with the smallest weight
 *  This algorithm combines the pole with the smallest weight with its nearest two neighbors
 *  and repeats this untill the number of nonzero poles is N
 */
void ReduceListOfPolesNumber(ListOfPoles * M, uint32_t N)
{
  if(M->N<N) return;
  uint32_t i,imin;
  double biminLeft,biminRight,temp;
  OrderListOfPoles(M);
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
        M->a[imin+1]=((M->b[imin]/M->b[imin+1])*M->a[imin]+M->a[imin+1])/((M->b[imin]/M->b[imin+1])+1);
        M->b[imin+1]=M->b[imin]+M->b[imin+1];
      }
      else if(imin==M->N-1)
      {
        M->a[imin-1]=((M->b[imin]/M->b[imin-1])*M->a[imin]+M->a[imin-1])/((M->b[imin]/M->b[imin-1])+1);
        M->b[imin-1]=M->b[imin]+M->b[imin-1];
      }

      else
      {
        if(fabs((M->a[imin+1]-M->a[imin-1]))<fabs(10*(DBL_MIN+DBL_EPSILON*M->a[imin])))
        {
          biminLeft=M->b[imin]*(0.5);
          biminRight=M->b[imin]*(0.5);
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
          biminLeft =M->b[imin]*((M->a[imin+1]-M->a[imin  ])/(M->a[imin+1]-M->a[imin-1]));
          biminRight=M->b[imin]*((M->a[imin  ]-M->a[imin-1])/(M->a[imin+1]-M->a[imin-1]));
        }
        M->a[imin+1]=((biminRight/M->b[imin+1])*M->a[imin]+M->a[imin+1])/((biminRight/M->b[imin+1])+1);
        M->b[imin+1]=biminRight+M->b[imin+1];
        M->a[imin-1]=((biminLeft /M->b[imin-1])*M->a[imin]+M->a[imin-1])/((biminLeft /M->b[imin-1])+1);
        M->b[imin-1]=biminLeft + M->b[imin-1];
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


/*  Create a new set of poles, which has lenght N, by removing the M->N-N poles such that for each reduction weight * (delta E neighbor) is minimum
 *  This algorithm combines the pole with the smallest moment to one neighbor with its nearest two neighbors
 *  and repeats this untill the number of nonzero poles is N
 */
void ReduceListOfPolesMoment(ListOfPoles * M, uint32_t N)
{
  uint32_t i,imin;
  double biminLeft,biminRight,temp,MinMoment,TestMoment;
  if(M->N<N) return;
  if(N==0)
  {
    M->N=0;
    return;
  }
  if(N==1)
  {
    M->a[0]*=M->b[0];
    for(i=1;i<M->N;i++)
    {
      M->a[0]+=M->b[i]*M->a[i];
    }
    for(i=1;i<M->N;i++)
    {
      M->b[0]+=M->b[i];
    }
    M->a[0]/=(M->b[0]);
    M->N=1;
    return;
  }
  OrderListOfPoles(M);
  while(M->N>N)
  {
    imin=0;
    MinMoment=M->b[imin]*(M->a[imin+1]-M->a[imin]);
    TestMoment=M->b[M->N-1]*(M->a[M->N-1]-M->a[M->N-2]);
    if(TestMoment<MinMoment)
    {
      MinMoment=TestMoment;
      imin=M->N-1;
    }
    for(i=1;i<M->N-1;i++)
    {
      TestMoment=M->b[i]*Min(M->a[i+1]-M->a[i],M->a[i]-M->a[i-1]);
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
        M->a[imin+1]=((M->b[imin]/M->b[imin+1])*M->a[imin]+M->a[imin+1])/((M->b[imin]/M->b[imin+1])+1);
        M->b[imin+1]=M->b[imin]+M->b[imin+1];
      }
      else if(imin==M->N-1)
      {
        M->a[imin-1]=((M->b[imin]/M->b[imin-1])*M->a[imin]+M->a[imin-1])/((M->b[imin]/M->b[imin-1])+1);
        M->b[imin-1]=M->b[imin]+M->b[imin-1];
      }
      
      else
      {
        if(fabs((M->a[imin+1]-M->a[imin-1]))<fabs(10*(DBL_MIN+DBL_EPSILON*M->a[imin])))
        {
          biminLeft=M->b[imin]*(0.5);
          biminRight=M->b[imin]*(0.5);
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
          biminLeft =M->b[imin]*((M->a[imin+1]-M->a[imin  ])/(M->a[imin+1]-M->a[imin-1]));
          biminRight=M->b[imin]*((M->a[imin  ]-M->a[imin-1])/(M->a[imin+1]-M->a[imin-1]));
        }
        M->a[imin+1]=((biminRight/M->b[imin+1])*M->a[imin]+M->a[imin+1])/((biminRight/M->b[imin+1])+1);
        M->b[imin+1]=biminRight+M->b[imin+1];
        M->a[imin-1]=((biminLeft /M->b[imin-1])*M->a[imin]+M->a[imin-1])/((biminLeft /M->b[imin-1])+1);
        M->b[imin-1]=biminLeft + M->b[imin-1];
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

/*  Create a new set of poles, wich are minimally spaced by da, by combining poles in a fixed region.
 *
 */
void ReduceListOfPolesSpacing(ListOfPoles * M, double da)
{
  uint32_t iRed,i;
  double PoleBorder;
  OrderListOfPoles(M);
  iRed=0;
  PoleBorder=M->a[iRed]+da/2.0;
  for(i=1;i<M->N;i++)
  {
    if(PoleBorder>M->a[i])
    {
      if(fabs(M->b[i]+M->b[iRed])<10*DBL_MIN)
      {
        M->b[iRed]=0;
        M->a[iRed]=M->a[Min(i+1,M->N-1)];
        M->b[i]=0;
      }
      else
      {
        M->a[iRed]=(M->a[iRed]*M->b[iRed]+M->a[i]*M->b[i])/(M->b[i]+M->b[iRed]);
        M->b[iRed]+=M->b[i];
        M->b[i]=0;
      }
    }
    else
    {
      iRed++;
      if(iRed!=i)
      {
        M->a[iRed]=M->a[i];
        M->b[iRed]=M->b[i];
        M->b[i]=0;
      }
      PoleBorder+=da;
    }
  }
  M->N=iRed+1;
}

/*   A <- A-B
 *   The energy of the poles in B is adjusted to be on the same "grid" as the poles in A
 *   If A pole in A becomes negative weight it is substracted from the remainder of A
 *   If A has only one pole left the returned weight can be negative
 */
uint32_t SubtractListOfPoles(ListOfPoles *const A, const ListOfPoles B)
{
  uint32_t i,NRes,iMin;
  double VRes[A->N+B.N], WRes[A->N+B.N];
  if(A->N==0){printf("SubtractListOfPoles called with NA=0\n");return 1;};
  OrderListOfPoles(A);
  A->a0-=B.a0;
  /* copy B to Res */
  for(i=0;i<B.N;i++)
  {
    VRes[i]=B.a[i];
    WRes[i]=B.b[i];
  }
  NRes=B.N;
  /* subtract B from A */
  while (NRes!=0)
  {
    /* find neighbors in A for VRes[NRes-1] */
    /* special case *NA=1 */
    if(A->N==1)
    {
      if(fabs(A->b[0]-WRes[NRes-1])>(DBL_MIN+DBL_EPSILON*fabs(WRes[NRes-1]+A->b[0])))
      {
        A->a[0]+=WRes[NRes-1]*(A->a[0]-VRes[NRes-1])/(A->b[0]-WRes[NRes-1]);
      }
      A->b[0]-=WRes[NRes-1];
      NRes--;
    }
    else
    {
      if(VRes[NRes-1]>A->a[A->N-2])
      {
        iMin=A->N-2;
      }
      else
      {
        for(i=1;i<A->N;i++)
        {
          if(VRes[NRes-1]<A->a[i])
          {
            iMin=i-1;
            break;
          }
        }
      }
      if(fabs(A->a[iMin+1]-A->a[iMin])<(DBL_MIN+DBL_EPSILON*fabs(A->a[iMin+1]-VRes[NRes-1])))
      {
        A->b[iMin  ]-=0.5*WRes[NRes-1];
        A->b[iMin+1]-=0.5*WRes[NRes-1];
      }
      else
      {
        A->b[iMin+1]-=((VRes[NRes-1]-A->a[iMin  ])/(A->a[iMin+1]-A->a[iMin]))*WRes[NRes-1];
        A->b[iMin  ]-=((A->a[iMin+1]-VRes[NRes-1])/(A->a[iMin+1]-A->a[iMin]))*WRes[NRes-1];
      }
      NRes--;
      if(A->b[iMin+1]<0)
      {
        WRes[NRes]=-A->b[iMin+1];
        VRes[NRes]=A->a[iMin+1];
        NRes++;
        for(i=iMin+2;i<A->N;i++)
        {
          A->b[i-1]=A->b[i];
          A->a[i-1]=A->a[i];
        }
        A->N--;
      }
      if(A->b[iMin]<0)
      {
        WRes[NRes]=-A->b[iMin];
        VRes[NRes]=A->a[iMin];
        NRes++;
        for(i=iMin+1;i<A->N;i++)
        {
          A->b[i-1]=A->b[i];
          A->a[i-1]=A->a[i];
        }
        A->N--;
      }
    }
  }
  if(A->b[0]<0) {printf("SubtractListOfPoles returned with A->b[0]=%24.17E and A->N=%i\n",A->b[0],A->N);return 2;};
  return 0;
}
