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
  return 0;
}

void FreeListOfPoles(ListOfPoles *M)
{
  if(M->N!=0)
  {
    free(M->a);
    free(M->b); 
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
  A->N+=B.N;
  for(i=0;i<B.N;i++)
  {
    A->a[i+A->N]=B.a[i];
    A->b[i+A->N]=B.b[i];
  }
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

/*  Create a new set of poles, wich are minimally spaced by da, by combining poles in a fixed region.
 *
 */
void ReduceListOfPoles(ListOfPoles * M, double da)
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
  M->a=realloc(M->a, M->N*sizeof(double));
  M->b=realloc(M->b, M->N*sizeof(double));
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
