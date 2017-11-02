/*
*     ( A[0][0]=a1 | A[0][1]=b1 | A[0][2]=b2 | A[0][3]=b3 )
*   A=( A[1][0]=b1 | A[1][1]=a2 | 0          | 0          )
*     ( A[2][0]=b2 | 0          | A[2][2]=a3 | 0          )
*     ( A[3][0]=b3 | 0          | 0          | A[3][3]=a4 )
*
*   (A B)[i][j] = \sum_k A[i][k] B[k][j]
*
*/



/* Free an AndersonMatrix M (the structure is not freed)
 */
void FreeAndersonMatrix(AndersonMatrixType * const M)
{
  if (M->N != 0)
  {
    free(M->a);
    free(M->b);
  }
}

/*  Init AndersonMatrix
*   M.N must be set
*/
uint32_t InitAndersonMatrix(AndersonMatrixType * const M)
{
  strcpy(M->Name, "Matrix");
  if (M->N == 0)
  {
    printf("can't initialize an AndersonMatrix of size 0\nError in InitAndersonMatrix\n");
    fflush(stdout);
    return 2;
  };
  M->a = (double*)malloc((M->N + 1) * sizeof(double));
  if (M->a == NULL)
  {
    printf("malloc failed in InitAndersonMatrix\n");
    fflush(stdout);
    return 1;
  };
  M->b = (double*)malloc((M->N) * sizeof(double));
  if (M->b == NULL)
  {
    free(M->a);
    printf("malloc failed in InitAndersonMatrix\n");
    fflush(stdout);
    return 1;
  }
  M->a[0] = 0;
  M->b[0] = 1;
  M->mu = 0;
  return 0;
}

/*  Init Anderson matrix and set diagonal (a) equal to one
 *  M.N must be set
 */
uint32_t InitIdentityAndersonMatrix(AndersonMatrixType * const M)
{
  uint32_t i;
  strcpy(M->Name, "I_Matrix");
  if (M->N == 0)
  {
    printf("can't initialize an AndersonMatrix of size 0\nError in InitIdentityAndersonMatrix\n");
    fflush(stdout);
    return 2;
  };
  M->a = (double*)malloc((M->N + 1) * sizeof(double));
  if (M->a == NULL)
  {
    printf("malloc failed in InitIdentityAndersonMatrix\n");
    fflush(stdout);
    return 1;
  };
  M->b = calloc((M->N), sizeof(double));
  if (M->b == NULL)
  {
    free(M->a);
    printf("malloc failed in InitIdentityAndersonMatrix\n");
    fflush(stdout);
    return 1;
  }
  for (i = 1; i <= M->N; i++)
  {
    M->a[i] = 1;
  }
  M->a[0] = 0;
  M->b[0] = 1;
  M->mu = 0;
  return 0;
}

/* set elements with an absolute value smaller then epsilon to zero
 */
void AndersonMatrixChop(AndersonMatrixType * const M, double const epsilon)
{
  uint32_t i;
  for (i = 0; i < M->N; i++)
  {
    if ((fabs(M->a[i]) < epsilon))
      M->a[i] = 0;
    if ((fabs(M->b[i]) < epsilon))
      M->b[i] = 0;
  }
  if ((fabs(M->a[M->N]) < epsilon))
    M->a[M->N] = 0;
}

/*  MI=M^-1
 *  MI not initiallized on call
 */
double AndersonMatrixInverse(AndersonMatrixType const M, MatrixType * const MI)
{
  uint32_t i;
  MI->Ni = M.N;
  MI->Nj = M.N;
  MI->QComplex = 0;
  if (InitMatrix(MI))
  {
    printf("InitMatrix failed in AndrsonMatrixInverse\n");
    fflush(stdout);
    return -1;
  };
  for (i = 1; i < M.N; i++)
  {
    MI->VRe[0][i] = M.b[i];
    MI->VRe[i][0] = M.b[i];
    MI->VRe[i][i] = M.a[i + 1];
  }
  MI->VRe[0][0] = M.a[1];
  return MatrixInverse(MI);
}


/*
*  M<- M + alpha B
*/
void AndersonMatrixAdd(AndersonMatrixType * const M, AndersonMatrixType B, double alphaRe)
{
  uint32_t i;
  for (i = 1; i < M->N; i++)
  {
    M->a[i] += alphaRe * B.a[i];
    M->b[i] += alphaRe * B.b[i];
  }
  M->a[M->N] += alphaRe * B.a[M->N];
  M->a[0] = 0;
  M->b[0] = 1;
}

/*
*  M<- alpha M
*/
void AndersonMatrixScalarMultiply(AndersonMatrixType * const M, double alphaRe)
{
  uint32_t i;
  for (i = 1; i < M->N; i++)
  {
    M->a[i] *= alphaRe;
    M->b[i] *= alphaRe;
  }
  M->a[M->N] *= alphaRe;
  M->a[0] = 0;
  M->b[0] = 1;
}


/* Create a new matrix that has dimensions (2+iEnd-iStart,2+iEnd-iStart)
 * and is a copy of a part of A
 * namely A[{0,iStart;;iEnd}][{0,iStart;;iEnd}]
 * M not allocated on call of the subroutine
 */
uint32_t AndersonMatrixPart(AndersonMatrixType * const M, AndersonMatrixType A, const uint32_t iStart, const uint32_t iEnd)
{
  uint32_t i;
  if ((iEnd < iStart) || (iEnd >= A.N) || (iStart == 0))
    return 1;
  M->N = 2 + iEnd - iStart;
  InitAndersonMatrix(M);
  M->a[0] = A.a[0];
  M->b[0] = A.b[0];
  M->a[1] = A.a[1];
  for (i = 2; i <= M->N; i++)
  {
    M->a[i] = A.a[i + iStart];
    M->b[i - 1] = A.b[i - 1 + iStart];
  }
  return 0;
}

/*  M<-A
 */
uint32_t AndersonMatrixCopy(AndersonMatrixType * const M, AndersonMatrixType A)
{
  uint32_t i;
  M->N = A.N;
  if (InitAndersonMatrix(M))
    return 1;
  for (i = 0; i < M->N; i++)
  {
    M->a[i] = A.a[i];
    M->b[i] = A.b[i];
  }
  M->a[M->N] = A.a[M->N];
  return 0;
}

/*  M<-A
 *  and M.a[i]<=M.a[i+1] (for i>1)
 */
uint32_t AndersonMatrixOrder(AndersonMatrixType * const M, AndersonMatrixType A)
{
  uint32_t i, j;
  double temp;
  M->N = A.N;
  if (InitAndersonMatrix(M))
    return 1;
  M->a[0] = A.a[0];
  M->b[0] = A.b[0];
  if (A.N == 0)
    return 0;
  M->a[1] = A.a[1];
  if (A.N == 1)
    return 0;
  M->a[2] = A.a[2];
  M->b[1] = A.b[1];
  for (i = 3; i <= M->N; i++)
  {
    M->a[i] = A.a[i];
    M->b[i - 1] = A.b[i - 1];
    for (j = i; j > 2; j--)
    {
      if (M->a[j] < M->a[j - 1])
      {
        temp = M->a[j - 1];
        M->a[j - 1] = M->a[j];
        M->a[j] = temp;
        temp = M->b[j - 2];
        M->b[j - 2] = M->b[j - 1];
        M->b[j - 1] = temp;
      }
      else
        break;
    }
  }
  return 0;
}


/*  M<-A
 *  and M.a[i]<=M.a[i+1] (for i>1)
 *  furthermore M.b[i] = fabs ( M.b[i] )
 */
uint32_t AndersonMatrixOrderAndAbs(AndersonMatrixType * const M, AndersonMatrixType A)
{
  uint32_t i, j;
  double temp;
  M->N = A.N;
  if (InitAndersonMatrix(M))
    return 1;
  M->a[0] = A.a[0];
  M->b[0] = fabs(A.b[0]);
  if (A.N == 0)
    return 0;
  M->a[1] = A.a[1];
  if (A.N == 1)
    return 0;
  M->a[2] = A.a[2];
  M->b[1] = fabs(A.b[1]);
  for (i = 3; i <= M->N; i++)
  {
    M->a[i] = A.a[i];
    M->b[i - 1] = fabs(A.b[i - 1]);
    for (j = i; j > 2; j--)
    {
      if (M->a[j] < M->a[j - 1])
      {
        temp = M->a[j - 1];
        M->a[j - 1] = M->a[j];
        M->a[j] = temp;
        temp = M->b[j - 2];
        M->b[j - 2] = M->b[j - 1];
        M->b[j - 1] = temp;
      }
      else
        break;
    }
  }
  return 0;
}



/* y <- M x
 * y already allocated
 * y[i]=sum_j M[i][j] x[j]
 */
void RealAndersonMatrixVector(AndersonMatrixType const M, const double *const x, double *const y)
{
  uint32_t i;
  for (i = 0; i < M.N; i++)
  {
    y[i] = M.a[i + 1] * x[i];
  }
  for (i = 1; i < M.N; i++)
  {
    y[0] += M.b[i] * x[i];
    y[i] += M.b[i] * x[0];
  }
}


/* y <- M x
 * y already allocated
 * y[i]=sum_j M[i][j] x[j]
 */
void ComplexAndersonMatrixVector(AndersonMatrixType const M, const double *const xRe, const double *const xIm, double *const yRe, double *const yIm)
{
  uint32_t i;
  for (i = 0; i < M.N; i++)
  {
    yRe[i] = M.a[i + 1] * xRe[i];
    yIm[i] = M.a[i + 1] * xIm[i];
  }
  for (i = 1; i < M.N; i++)
  {
    yRe[0] += M.b[i] * xRe[i];
    yRe[i] += M.b[i] * xRe[0];
    yIm[0] += M.b[i] * xIm[i];
    yIm[i] += M.b[i] * xIm[0];
  }
}


/* y <- x M x
 * y already allocated
 * y=sum_ij x[i] M[i][j] x[j]
 */
void RealVectorAndersonMatrixVector(AndersonMatrixType const M, const double *const x, double *const y)
{
  uint32_t i;
  *y = 0;
  for (i = 0; i < M.N; i++)
  {
    *y += x[i] * M.a[i + 1] * x[i];
  }
  for (i = 1; i < M.N; i++)
  {
    *y += 2 * x[i] * M.b[i] * x[0];
  }
}

/* y <- x^* M x
 * y already allocated
 * y=sum_ij x[i]^* M[i][j] x[j]
 */
void ComplexVectorAndersonMatrixVector(AndersonMatrixType const M, const double *const xRe, const double *const xIm, double *const yRe, double *const yIm)
{
  uint32_t i;
  *yRe = 0;
  *yIm = 0;
  for (i = 0; i < M.N; i++)
  {
    *yRe += xRe[i] * M.a[i + 1] * xRe[i] + xIm[i] * M.a[i + 1] * xIm[i];
    *yIm += xRe[i] * M.a[i + 1] * xIm[i] - xIm[i] * M.a[i + 1] * xRe[i];
  }
  for (i = 1; i < M.N; i++)
  {
    *yRe += 2 * (xRe[i] * M.b[i] * xRe[0] + xIm[i] * M.b[i] * xIm[0]);
  }
}




/*
 *
 */
uint32_t AndersonMatrixCalculateChemicalPotential(AndersonMatrixType const M, double const occupation, double *const mu)
{
  ListOfPoles MLP;
  AndersonMatrixToListOfPoles(M, &MLP);
  if (ListOfPolesCalculateChemicalPotential(MLP, occupation, mu))
  {
    printf("ListOfPolesCalculateChemicalPotential failed in AndersonMatrixCalculateChemicalPotential\n");
    FreeListOfPoles(&MLP);
    return 1;
  }
  FreeListOfPoles(&MLP);
  return 0;
}

/* Given an anderson impurity Hammiltonian represented as an one particle matrix and a chemical potential
 * the impurity occupation is calculated
 */
uint32_t CalculateImpurityOccupationAndersonMatrix(AndersonMatrixType const M, double * const occupation)
{
  uint32_t i;
  ListOfPoles MLP;
  if(AndersonMatrixToListOfPoles(M, &MLP))
  {
    printf("AndersonMatrixToListOfPoles failed in CalculateImpurityOccupationAndersonMatrix\n");
  }
  *occupation=0;
  for(i=0;i<MLP.N;i++)
  {
    if(MLP.a[i]==MLP.mu)
    {
      *occupation+=MLP.b[i]/2.0;
    }
    else if(MLP.a[i]<MLP.mu)
    {
      *occupation+=MLP.b[i];
    }
  }
  FreeListOfPoles(&MLP);
  return 0;
}

/*
 *
 */
uint32_t DetermineImpurityShiftFromOccupationAndersonMatrix(AndersonMatrixType M, double const occupation, double *const shift)
{
  uint32_t MaxIt=30, History=10;
  uint32_t i;
  double Shifti[MaxIt+1], occupationi[MaxIt+1], totalb=0;
  Shifti[0]=*shift;
  M.a[1]+=Shifti[0];
  CalculateImpurityOccupationAndersonMatrix(M, &occupationi[0]);
  M.a[1]-=Shifti[0];
  for(i=0;i<M.N;i++)
  {
    totalb+=Sqr(M.b[i]);
  }
  totalb=sqrt(totalb);
  i=0;
  while((i<MaxIt)&&(Abs(occupationi[i]-occupation)>(100*(DBL_MIN+DBL_EPSILON*Abs(occupation)))))
  {
    i++;
    if(i==1)
    {
      if(occupationi[0]>occupation)
      {
        Shifti[i]=Shifti[i-1]+totalb;
      }
      else
      {
        Shifti[i]=Shifti[i-1]-totalb;
      }
    }
    else
    {
      Shifti[i]=InterpolateRootOfset(&Shifti[Max(History,i)-History], &occupationi[Max(History,i)-History], Min(i-1,History), occupation);
    }
    M.a[1]+=Shifti[i];
    CalculateImpurityOccupationAndersonMatrix(M, &occupationi[i]);
    M.a[1]-=Shifti[i];
  }
  *shift=Shifti[i];
  
  MaxIt=i;
  
  
  for(i=0;i<MaxIt+1;i++)
  {
    printf("shift %22.15E occupation %22.15E deviation %22.15E\n",Shifti[i], occupationi[i], occupationi[i]-occupation);
  }
  
  return 0;
}
