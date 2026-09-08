/*   Multiplication and other algorithms accept the following transformations as string
 *   'T' Transpose M^T(i,j) = M(j,i)
 *   'N' Do nothing
 *   '*' Conjugate M^*(i,j).Re = M(i,j).Re  M^*(i,j).Im = -M(i,j).Im
 *   'C' Conjugate Transpose
 */



/* Free a matrix M itselve (i.e. the structure is not freed)
 */
void FreeCompactMatrix(CompactMatrixType * const M)
{
  if(M->OnLuaStack) return;
  free(M->Name);
  if(M->QComplex)
  {
    free(M->V);
  }
  else
  {
    free(M->VRe);
  }
}

void FreeCompactMatrixIncludingPointer(CompactMatrixType * const M)
{
  if(M->OnLuaStack) return;
  FreeCompactMatrix(M);
  free(M);
}

/* Initiallize a matrix M
 *  M.QComplex, M.Ni, M.Nj must be set
 */
uint32_t InitCompactMatrix(CompactMatrixType * const M)
{
  M->Name=(char*)malloc(MAXStrLength*sizeof(char));
  if(M->Name==NULL){printf("malloc for name failed in InitCompactMatrix\n");return 1;};
  strcpy(M->Name, "Matrix");
  M->a0=0;
  M->b0=1;
  M->mu=0;
  if(M->QComplex)
  {
    M->V=(Complex16*)calloc(M->Ni*M->Nj,sizeof(Complex16));
    if(M->V==NULL)
    {
      printf("Malloc failed in InitCompactMatrix for complex allocation\n");
      fflush(stdout);
      return 1;
    }
  }
  else
  {
    M->VRe=(double*)calloc(M->Ni*M->Nj,sizeof(double));
    if(M->VRe==NULL)
    {
      printf("Malloc failed in InitCompactMatrix for real allocation\n");
      fflush(stdout);
      return 1;
    }
  }
  M->OnLuaStack=false;
  return 0;
}

uint32_t ResizeCompactMatrix(CompactMatrixType * const M, uint32_t Ni, uint32_t Nj)
{
  if(M->QComplex)
  {
    M->V=(Complex16*)realloc(M->V,Ni*Nj*sizeof(Complex16));
    if(M->V==NULL)
    {
      printf("Realloc failed in ResizeCompactMatrix for complex allocation\n");
      fflush(stdout);
      return 1;
    }
  }
  else
  {
    M->VRe=(double*)realloc(M->VRe,Ni*Nj*sizeof(double));
    if(M->VRe==NULL)
    {
      printf("Realloc failed in ResizeCompactMatrix for real allocation\n");
      fflush(stdout);
      return 1;
    }
  }
  M->Ni=Ni;
  M->Nj=Nj;
  return 0;
}

/*  Set M to an empty matrix that owns no data. FreeCompactMatrix on such a matrix only calls
 *  free(NULL), so it is safe to leave it to the lua garbage collector. Use it on a freshly
 *  allocated CompactMatrixType before filling it in, when the filling in can fail.
 */
uint32_t SetCompactMatrixEmpty(CompactMatrixType * const M)
{
  M->Name=NULL;
  M->Ni=0;
  M->Nj=0;
  M->QComplex=0;
  M->VRe=NULL;
  M->a0=0;
  M->b0=0;
  M->mu=0;
  M->OnLuaStack=false;
  return 0;
}

uint32_t SetCompactMatrixToZero(CompactMatrixType * const M)
{
  uint32_t i;
  if(M->QComplex)
  {
    for(i=0;i<M->Ni*M->Nj;i++)
    {
      M->V[i].Re=0;
      M->V[i].Im=0;
    }
  }
  else
  {
    for(i=0;i<M->Ni*M->Nj;i++)
    {
      M->VRe[i]=0;
    }
  }
  return 0;
}

/* B <- A
 * B should not be initiallized
 */
uint32_t CompactMatrixCopy(CompactMatrixType const A, CompactMatrixType * const B)
{
  B->Ni=A.Ni;
  B->Nj=A.Nj;
  B->QComplex=A.QComplex;
  if(InitCompactMatrix(B))
  {
    printf("InitCompactMatrix failed in CompactMatrixCopy\n");
    return 1;
  }
  strncpy(B->Name, A.Name, MAXStrLength);
  /*  after InitCompactMatrix, which sets a0=0, b0=1 and mu=0 */
  B->a0=A.a0;
  B->b0=A.b0;
  B->mu=A.mu;
  if(A.QComplex)
  {
    memcpy(B->V,A.V,A.Ni*A.Nj*sizeof(Complex16));
  }
  else
  {
    memcpy(B->VRe,A.VRe,A.Ni*A.Nj*sizeof(double));
  }
  return 0;
}

/* A <-> B
 * Both matrices should be initiallized
 */
uint32_t CompactMatrixSwap(CompactMatrixType * const A, CompactMatrixType * const B)
{
  CompactMatrixType Temp;
  Temp = *A;
  *A = *B;
  *B = *A;
  return 0;
}


/*  M <- A
 *  M should not be allocated
 */
uint32_t CompactMatrixCopyToMatrix(MatrixType * const M, CompactMatrixType const A)
{
  uint32_t i;
  M->Ni = A.Ni;
  M->Nj = A.Nj;
  M->QComplex = A.QComplex;
  if (InitMatrix(M)) return 1;
  strcpy(M->Name, A.Name);
  if (M->QComplex)
  {
    for (i = 0; i < M->Ni*M->Nj; i++)
    {
      M->VRe[0][i] = A.V[i].Re;
      M->VIm[0][i] = A.V[i].Im;
    }
  }
  else
  {
    for (i = 0; i < M->Ni*M->Nj; i++)
    {
      M->VRe[0][i] = A.VRe[i];
    }
  }
  return 0;
}

/* Initiallize a matrix M
 *  M.QComplex, M.Ni, M,Nj must be set
 *  M(i,i) is set to 1
 */
uint32_t InitCompactIdentityMatrix(CompactMatrixType * const M)
{
  uint32_t i;
  M->Name=(char*)malloc(MAXStrLength*sizeof(char));
  if(M->Name==NULL){printf("malloc for name failed in InitCompactMatrix\n");return 1;};
  strcpy(M->Name, "I_Matrix");
  if(M->QComplex)
  {
    M->V=(Complex16*)calloc(M->Ni*M->Nj,sizeof(Complex16));
    if(M->V==NULL)
    {
      printf("Calloc failed in InitCompactMatrix for complex allocation\n");
      fflush(stdout);
      return 1;
    }
    for(i=0;i<Min(M->Ni,M->Nj);i++)
    {
      M->V[i*M->Nj+i].Re=1;
    }
  }
  else
  {
    M->VRe=(double*)calloc(M->Ni*M->Nj,sizeof(double));
    if(M->V==NULL)
    {
      printf("Calloc failed in InitCompactMatrix for real allocation\n");
      fflush(stdout);
      return 1;
    }
    for(i=0;i<Min(M->Ni,M->Nj);i++)
    {
      M->VRe[i*M->Nj+i]=1;
    }
  }
  M->OnLuaStack=false;
  return 0;
}

/* turn a real matrix into a complex matrix i.e. allocate more memory and zero this
 */
uint32_t MakeCompactRealMatrixComplex(CompactMatrixType * const M)
{
  Complex16 *temp;
  uint32_t i;
  if (M->QComplex) return 0;
  temp=(Complex16*)calloc(M->Ni*M->Nj,sizeof(Complex16));
  if(temp==NULL)
  {
    printf("calloc failed in MakeCompactRealMatrixComplex\n");
    return 1;
  }
  for(i=0;i<M->Ni*M->Nj;i++)
  {
    temp[i].Re=M->VRe[i];
  }
  free(M->VRe);
  M->V=temp;
  M->QComplex=1;
  return 0;
}


/* turn a complex matrix into a real matrix i.e. ignore complex part
 */
uint32_t MakeCompactComplexMatrixReal(CompactMatrixType * const M)
{
  double *temp;
  uint32_t i;
  if (!(M->QComplex)) return 0;
  temp=(double*)calloc(M->Ni*M->Nj,sizeof(double));
  if(temp==NULL)
  {
    printf("calloc failed in MakeCompactComplexMatrixReal\n");
    return 1;
  }
  for(i=0;i<M->Ni*M->Nj;i++)
  {
    temp[i]=M->V[i].Re;
  }
  free(M->V);
  M->VRe=temp;
  M->QComplex=0;
  return 0;
}

/* set elements with an absolute value smaller then epsilon to zero
 */
uint32_t CompactMatrixChop(CompactMatrixType * const M, double const epsilon)
{
  int    i, QStillComplex;
  QStillComplex = 0;
  if (M->QComplex)
  {
    for (i = 0; i < M->Ni*M->Nj; i++)
    {
      if (fabs(M->V[i].Re) < epsilon)
      {
        M->V[i].Re = 0;
      }
      if (fabs(M->V[i].Im) < epsilon)
      {
        M->V[i].Im = 0;
      }
      else
      {
          QStillComplex = 1;
      }
    }
    if (QStillComplex == 0)
    {
      if(MakeCompactComplexMatrixReal(M))
      {
        return 1;
      }
    }
  }
  else
  {
    for (i = 0; i < M->Ni*M->Nj; i++)
    {
      if (fabs(M->VRe[i]) < epsilon)
      {
        M->VRe[i] = 0;
      }
    }
  }
  return 0;
}

uint32_t CompactMatrixRemoveRow(CompactMatrixType * const A, uint32_t row)
{
  if(row>=A->Ni)
  {
    printf("Error in CompactMatrixRemoveRow, tried to remove row %u, but there are only %u rows (counting starts at zero)\n",row,A->Ni);
    return 1;
  }
  if(A->QComplex)
  {
    if(!(row+1==A->Ni)) memmove(&(A->V[row*A->Nj]), &(A->V[(row+1)*A->Nj]), (A->Ni-row-1)*A->Nj*sizeof(Complex16));
    A->Ni--;
    A->V= (Complex16*)realloc(A->V, A->Ni*A->Nj*sizeof(Complex16));
  }
  else
  {
    if(!(row+1==A->Ni)) memmove(&(A->VRe[row*A->Nj]), &(A->VRe[(row+1)*A->Nj]), (A->Ni-row-1)*A->Nj*sizeof(double));
    A->Ni--;
    A->VRe= (double*)realloc(A->VRe, A->Ni*A->Nj*sizeof(double));
  }
  return 0;
}


/*   y = A x
 *   y is uninitialized on entry
 */
uint32_t CompactMatrixVector(CompactMatrixType const *A, CompactVectorType *x, CompactVectorType *y)
{
  uint32_t INC=1;
  char TransA='T';
  if(A->Nj!=x->Ni){printf("CompactMatrixVector incompatible matrix (%i by %i) vector (%i) dimensions",A->Ni,A->Nj,x->Ni);return 1;}
  if(A->QComplex&&x->QComplex)
  {
    Complex16 ALPHA,BETA;
    ALPHA.Re=1; ALPHA.Im=0;
    BETA.Re=0;  BETA.Im=0;
    y->V= (Complex16*)malloc(A->Ni*sizeof(Complex16));
    if(y==NULL){printf("malloc failed in CompactMatrixVector\n");return 1;}
    zgemv(&TransA,&(A->Nj),&(A->Ni),&ALPHA,A->V,&(A->Nj),x->V,&INC,&BETA,y->V,&INC);
  }
  else if(A->QComplex)
  {
    int i;
    Complex16 ALPHA,BETA;
    Complex16 xtemp[x->Ni];
    for(i=0;i<x->Ni;i++)
    {
      xtemp[i].Re=x->VRe[i];
      xtemp[i].Im=0;
    }
    ALPHA.Re=1; ALPHA.Im=0;
    BETA.Re=0;  BETA.Im=0;
    y->V= (Complex16*)malloc(A->Ni*sizeof(Complex16));
    if(y==NULL){printf("malloc failed in CompactMatrixVector\n");return 1;}
    zgemv(&TransA,&(A->Nj),&(A->Ni),&ALPHA,A->V,&(A->Nj),xtemp,&INC,&BETA,y->V,&INC);
  }
  else if(x->QComplex)
  {
    int i;
    Complex16 ALPHA,BETA;
    Complex16 Atemp[A->Ni*A->Nj];
    for(i=0;i<A->Ni*A->Nj;i++)
    {
      Atemp[i].Re=A->VRe[i];
      Atemp[i].Im=0;
    }
    ALPHA.Re=1; ALPHA.Im=0;
    BETA.Re=0;  BETA.Im=0;
    y->V= (Complex16*)malloc(A->Ni*sizeof(Complex16));
    if(y==NULL){printf("malloc failed in CompactMatrixVector\n");return 1;}
    zgemv(&TransA,&(A->Nj),&(A->Ni),&ALPHA,Atemp,&(A->Nj),x->V,&INC,&BETA,y->V,&INC);
  }
  else
  {
    double ALPHA=1,BETA=0;
    y->VRe= (double*)malloc(A->Ni*sizeof(double));
    y->Ni=A->Ni;
    if(y==NULL){printf("malloc failed in CompactMatrixVector\n");return 1;}
    dgemv(&TransA,&(A->Nj),&(A->Ni),&ALPHA,A->VRe,&(A->Nj),x->VRe,&INC,&BETA,y->VRe,&INC);
  }
  return 0;
}

/*   B <= A * B
 *   A is diagonal (only A[i,i] is nonzero and stored as a list of doubles
 */
uint32_t MultiplyDiagonalMatrixCompactMatrix(double const *A, CompactMatrixType * const B)
{
  int i,j,index=0;
  if(B->QComplex)
  {
    for(i=0;i<B->Ni;i++)
    {
      for(j=0;j<B->Nj;j++)
      {
        B->V[index].Re*=A[i];
        B->V[index].Im*=A[i];
        index++;
      }
    }
  }
  else
  {
    for(i=0;i<B->Ni;i++)
    {
      for(j=0;j<B->Nj;j++)
      {
        B->VRe[index]*=A[i];
        index++;
      }
    }
  }
  return 0;
}

/*   B <= B * A
 *   A is diagonal (only A[i,i] is nonzero and stored as a list of doubles
 */
uint32_t MultiplyCompactMatrixDiagonalMatrix(CompactMatrixType * const B, double const *A)
{
  int i,j,index=0;
  if(B->QComplex)
  {
    for(i=0;i<B->Ni;i++)
    {
      for(j=0;j<B->Nj;j++)
      {
        B->V[index].Re*=A[j];
        B->V[index].Im*=A[j];
        index++;
      }
    }
  }
  else
  {
    for(i=0;i<B->Ni;i++)
    {
      for(j=0;j<B->Nj;j++)
      {
        B->VRe[index]*=A[j];
        index++;
      }
    }
  }
  return 0;
}

/*  Calculates T A T^dag
 *  dimensions of input T : T.Ni by T.Nj
 *                      A : T.Nj by T.Nj
 *  dimensions of output: T.Ni by T.Ni
 */
uint32_t TransformDiagonalMatrixByCompactMatrixTATdag(double const *A, CompactMatrixType const *T, CompactMatrixType * const TATdag)
{
  uint32_t i,j,k,TATindex,icol,jcol;
  TATdag->Ni=T->Ni;
  TATdag->Nj=T->Ni;
  TATdag->QComplex = T->QComplex;
  InitCompactMatrix(TATdag);
  TATindex=0;
  if(T->QComplex)
  {
    for(i=0;i<TATdag->Ni;i++)
    {
      for(j=0;j<TATdag->Nj;j++)
      {
        icol=i*T->Nj;
        jcol=j*T->Nj;
        TATdag->V[TATindex].Re=0;
        TATdag->V[TATindex].Im=0;
        for(k=0;k<T->Nj;k++)
        {
          TATdag->V[TATindex].Re+=T->V[jcol].Re*A[k]*T->V[icol].Re + T->V[jcol].Im*A[k]*T->V[icol].Im;
          TATdag->V[TATindex].Im+=T->V[jcol].Re*A[k]*T->V[icol].Im - T->V[jcol].Im*A[k]*T->V[icol].Re;
          icol++;
          jcol++;
        }
        TATindex++;
      }
    }
  }
  else
  {
    for(i=0;i<TATdag->Ni;i++)
    {
      for(j=0;j<TATdag->Nj;j++)
      {
        icol=i*T->Nj;
        jcol=j*T->Nj;
        TATdag->VRe[TATindex]=0;
        for(k=0;k<T->Nj;k++)
        {
          TATdag->VRe[TATindex]+=T->VRe[jcol]*A[k]*T->VRe[icol];
          icol++;
          jcol++;
        }
        TATindex++;
      }
    }
  }
  return 0;
}

/*  Calculates T^con A T^tran
 *  dimensions of input T : T.Ni by T.Nj
 *                      A : T.Nj by T.Nj
 *  dimensions of output: T.Ni by T.Ni
 */
uint32_t TransformDiagonalMatrixByCompactMatrixTconATtran(double const * const A, CompactMatrixType const * const T, CompactMatrixType * const TconATtran)
{
  uint32_t i,j,k,TATindex,icol,jcol;
  TconATtran->Ni=T->Ni;
  TconATtran->Nj=T->Ni;
  TconATtran->QComplex = T->QComplex;
  InitCompactMatrix(TconATtran);
  TATindex=0;
  if(T->QComplex)
  {
    for(i=0;i<TconATtran->Ni;i++)
    {
      for(j=0;j<TconATtran->Nj;j++)
      {
        icol=i*T->Nj;
        jcol=j*T->Nj;
        TconATtran->V[TATindex].Re=0;
        TconATtran->V[TATindex].Im=0;
        for(k=0;k<T->Nj;k++)
        {
          TconATtran->V[TATindex].Re+=T->V[jcol].Re*A[k]*T->V[icol].Re + T->V[jcol].Im*A[k]*T->V[icol].Im;
          TconATtran->V[TATindex].Im-=T->V[jcol].Re*A[k]*T->V[icol].Im - T->V[jcol].Im*A[k]*T->V[icol].Re;
          icol++;
          jcol++;
        }
        TATindex++;
      }
    }
  }
  else
  {
    for(i=0;i<TconATtran->Ni;i++)
    {
      for(j=0;j<TconATtran->Nj;j++)
      {
        icol=i*T->Nj;
        jcol=j*T->Nj;
        TconATtran->VRe[TATindex]=0;
        for(k=0;k<T->Nj;k++)
        {
          TconATtran->VRe[TATindex]+=T->VRe[jcol]*A[k]*T->VRe[icol];
          icol++;
          jcol++;
        }
        TATindex++;
      }
    }
  }
  return 0;
}

uint32_t TransformComplexDiagonalMatrixByCompactMatrixTconATtran(Complex16 const *A, CompactMatrixType const * const T, CompactMatrixType * const TconATtran)
{
  uint32_t i,j,k,TATindex,icol,jcol;
  TconATtran->Ni=T->Ni;
  TconATtran->Nj=T->Ni;
  TconATtran->QComplex = 1;
  InitCompactMatrix(TconATtran);
  TATindex=0;
  if(T->QComplex)
  {
    for(i=0;i<TconATtran->Ni;i++)
    {
      for(j=0;j<TconATtran->Nj;j++)
      {
        icol=i*T->Nj;
        jcol=j*T->Nj;
        TconATtran->V[TATindex].Re=0;
        TconATtran->V[TATindex].Im=0;
        for(k=0;k<T->Nj;k++)
        {
          TconATtran->V[TATindex].Re+= T->V[jcol].Re*A[k].Re*T->V[icol].Re + T->V[jcol].Im*A[k].Re*T->V[icol].Im
                                      +T->V[jcol].Re*A[k].Im*T->V[icol].Im - T->V[jcol].Im*A[k].Im*T->V[icol].Re;
          TconATtran->V[TATindex].Im+=-T->V[jcol].Re*A[k].Re*T->V[icol].Im + T->V[jcol].Im*A[k].Re*T->V[icol].Re
                                      +T->V[jcol].Re*A[k].Im*T->V[icol].Re + T->V[jcol].Im*A[k].Im*T->V[icol].Im;
          icol++;
          jcol++;
        }
        TATindex++;
      }
    }
  }
  else
  {
    for(i=0;i<TconATtran->Ni;i++)
    {
      for(j=0;j<TconATtran->Nj;j++)
      {
        icol=i*T->Nj;
        jcol=j*T->Nj;
        TconATtran->V[TATindex].Re=0;
        TconATtran->V[TATindex].Im=0;
        for(k=0;k<T->Nj;k++)
        {
          TconATtran->V[TATindex].Re+=T->VRe[jcol]*A[k].Re*T->VRe[icol];
          TconATtran->V[TATindex].Im+=T->VRe[jcol]*A[k].Im*T->VRe[icol];
          icol++;
          jcol++;
        }
        TATindex++;
      }
    }
  }
  return 0;
}

/*  Calculates TL^con A TR^tran
 *  dimensions of input TL : TL.Ni by TL.Nj
 *                      A  : A.N   by A.N
 *                      TR : TR.Ni by TR.Nj
 *                 A.N == TR.Nj
 *                 TL.Nj = A.N
 *  dimensions of output: TL.Ni by TR.Ni
 */
uint32_t TransformComplexDiagonalMatrixByCompactMatrixTLconATRtran(Complex16 const *A, CompactMatrixType const *TL, CompactMatrixType const *TR, CompactMatrixType * const TLconATRtran)
{
  uint32_t i,j,k,TATindex,icol,jcol;
  if(TL->Nj!=TR->Nj)
  {
    printf("Error in input of TransformComplexDiagonalMatrixByCompactMatrixTLconATRtran, dimension of left and right matrix not compatible\n");
    return 1;
  }
  TLconATRtran->Ni=TL->Ni;
  TLconATRtran->Nj=TR->Ni;
  TLconATRtran->QComplex = 1;
  InitCompactMatrix(TLconATRtran);
  TATindex=0;
  if(TL->QComplex)
  {
    if(TR->QComplex)
    { // L and R complex
      for(i=0;i<TLconATRtran->Ni;i++)
      {
        for(j=0;j<TLconATRtran->Nj;j++)
        {
          icol=i*TL->Nj;
          jcol=j*TR->Nj;
          TLconATRtran->V[TATindex].Re=0;
          TLconATRtran->V[TATindex].Im=0;
          for(k=0;k<TL->Nj;k++)
          {
            TLconATRtran->V[TATindex].Re+=TR->V[jcol].Re*A[k].Re*TL->V[icol].Re+TR->V[jcol].Re*A[k].Im*TL->V[icol].Im
                                         +TR->V[jcol].Im*A[k].Re*TL->V[icol].Im-TR->V[jcol].Im*A[k].Im*TL->V[icol].Re;
            TLconATRtran->V[TATindex].Im+=TR->V[jcol].Re*A[k].Im*TL->V[icol].Re-TR->V[jcol].Re*A[k].Re*TL->V[icol].Im
                                         +TR->V[jcol].Im*A[k].Re*TL->V[icol].Re+TR->V[jcol].Im*A[k].Im*TL->V[icol].Im;
            icol++;
            jcol++;
          }
          TATindex++;
        }
      }
    }
    else
    { // L complex, R real
      for(i=0;i<TLconATRtran->Ni;i++)
      {
        for(j=0;j<TLconATRtran->Nj;j++)
        {
          icol=i*TL->Nj;
          jcol=j*TR->Nj;
          TLconATRtran->V[TATindex].Re=0;
          TLconATRtran->V[TATindex].Im=0;
          for(k=0;k<TL->Nj;k++)
          {
            TLconATRtran->V[TATindex].Re+=TR->VRe[jcol]*A[k].Re*TL->V[icol].Re+TR->VRe[jcol]*A[k].Im*TL->V[icol].Im;
            TLconATRtran->V[TATindex].Im+=TR->VRe[jcol]*A[k].Im*TL->V[icol].Re-TR->VRe[jcol]*A[k].Re*TL->V[icol].Im;
            icol++;
            jcol++;
          }
          TATindex++;
        }
      }
    }
  }
  else
  {
    if(TR->QComplex)
    { // L real, R complex
      for(i=0;i<TLconATRtran->Ni;i++)
      {
        for(j=0;j<TLconATRtran->Nj;j++)
        {
          icol=i*TL->Nj;
          jcol=j*TR->Nj;
          TLconATRtran->V[TATindex].Re=0;
          TLconATRtran->V[TATindex].Im=0;
          for(k=0;k<TL->Nj;k++)
          {
            TLconATRtran->V[TATindex].Re+=TR->V[jcol].Re*A[k].Re*TL->VRe[icol]-TR->V[jcol].Im*A[k].Im*TL->VRe[icol];
            TLconATRtran->V[TATindex].Im+=TR->V[jcol].Re*A[k].Im*TL->VRe[icol]+TR->V[jcol].Im*A[k].Re*TL->VRe[icol];
            icol++;
            jcol++;
          }
          TATindex++;
        }
      }
    }
    else
    { //L and R real
      for(i=0;i<TLconATRtran->Ni;i++)
      {
        for(j=0;j<TLconATRtran->Nj;j++)
        {
          icol=i*TL->Nj;
          jcol=j*TR->Nj;
          TLconATRtran->V[TATindex].Re=0;
          TLconATRtran->V[TATindex].Im=0;
          for(k=0;k<TL->Nj;k++)
          {
            TLconATRtran->V[TATindex].Re+=TR->VRe[jcol]*A[k].Re*TL->VRe[icol];
            TLconATRtran->V[TATindex].Im+=TR->VRe[jcol]*A[k].Im*TL->VRe[icol];
            icol++;
            jcol++;
          }
          TATindex++;
        }
      }
    }
  }
  return 0;
}


/*  Calculates T^dag A T
 *  dimensions of input T : T.Ni by T.Nj
 *                      A : T.Ni by T.Ni
 *  dimensions of output: T.Nj by T.Nj
 */
uint32_t TransformDiagonalMatrixByCompactMatrixTdagAT(double const *A, CompactMatrixType const *T, CompactMatrixType * const TdagAT)
{
  uint32_t i,j,k,TATindex,irow,jrow;
  TdagAT->Ni=T->Nj;
  TdagAT->Nj=T->Nj;
  TdagAT->QComplex = T->QComplex;
  InitCompactMatrix(TdagAT);
  TATindex=0;
  if(T->QComplex)
  {
    for(i=0;i<TdagAT->Ni;i++)
    {
      for(j=0;j<TdagAT->Nj;j++)
      {
        irow=i;
        jrow=j;
        TdagAT->V[TATindex].Re=0;
        TdagAT->V[TATindex].Im=0;
        for(k=0;k<T->Ni;k++)
        {
          TdagAT->V[TATindex].Re+=T->V[jrow].Re*A[k]*T->V[irow].Re + T->V[jrow].Im*A[k]*T->V[irow].Im;
          TdagAT->V[TATindex].Im+=T->V[jrow].Im*A[k]*T->V[irow].Re - T->V[jrow].Re*A[k]*T->V[irow].Im;
          irow+=T->Nj;
          jrow+=T->Nj;
        }
        TATindex++;
      }
    }
  }
  else
  {
    for(i=0;i<TdagAT->Ni;i++)
    {
      for(j=0;j<TdagAT->Nj;j++)
      {
        irow=i;
        jrow=j;
        TdagAT->VRe[TATindex]=0;
        for(k=0;k<T->Ni;k++)
        {
          TdagAT->VRe[TATindex]+=T->VRe[jrow]*A[k]*T->VRe[irow];
          irow+=T->Nj;
          jrow+=T->Nj;
        }
        TATindex++;
      }
    }
  }
  return 0;
}


/*  Calculates T^tran A T^con
 *  dimensions of input T : T.Ni by T.Nj
 *                      A : T.Ni by T.Ni
 *  dimensions of output: T.Nj by T.Nj
 */
uint32_t TransformDiagonalMatrixByCompactMatrixTtranATcon(double const *A, CompactMatrixType const *T, CompactMatrixType * const TtranATcon)
{
  uint32_t i,j,k,TATindex,irow,jrow;
  TtranATcon->Ni=T->Nj;
  TtranATcon->Nj=T->Nj;
  TtranATcon->QComplex = T->QComplex;
  InitCompactMatrix(TtranATcon);
  TATindex=0;
  if(T->QComplex)
  {
    for(i=0;i<TtranATcon->Ni;i++)
    {
      for(j=0;j<TtranATcon->Nj;j++)
      {
        irow=i;
        jrow=j;
        TtranATcon->V[TATindex].Re=0;
        TtranATcon->V[TATindex].Im=0;
        for(k=0;k<T->Ni;k++)
        {
          TtranATcon->V[TATindex].Re+=T->V[jrow].Re*A[k]*T->V[irow].Re + T->V[jrow].Im*A[k]*T->V[irow].Im;
          TtranATcon->V[TATindex].Im-=T->V[jrow].Im*A[k]*T->V[irow].Re - T->V[jrow].Re*A[k]*T->V[irow].Im;
          irow+=T->Nj;
          jrow+=T->Nj;
        }
        TATindex++;
      }
    }
  }
  else
  {
    for(i=0;i<TtranATcon->Ni;i++)
    {
      for(j=0;j<TtranATcon->Nj;j++)
      {
        irow=i;
        jrow=j;
        TtranATcon->VRe[TATindex]=0;
        for(k=0;k<T->Nj;k++)
        {
          TtranATcon->VRe[TATindex]+=T->VRe[jrow]*A[k]*T->VRe[irow];
          irow+=T->Nj;
          jrow+=T->Nj;
        }
        TATindex++;
      }
    }
  }
  return 0;
}


/*  Calculates TR^tran A TL^con
 *  dimensions of input TL : TL.Ni by TL.Nj
 *                      A  : A.N   by A.N
 *                      TR : TR.Ni by TR.Nj
 *                 A.N == TR.Ni
 *                 TL.Ni = A.N
 *  dimensions of output: TR.Nj by TL.Nj
 */
uint32_t TransformComplexDiagonalMatrixByCompactMatrixTRtranATLcon(Complex16 const *A, CompactMatrixType const *TL, CompactMatrixType const *TR, CompactMatrixType * const TRtranATLcon)
{
  uint32_t i,j,k,TATindex,irow,jrow;
  if(TL->Ni!=TR->Ni)
  {
    printf("Error in input of TransformComplexDiagonalMatrixByCompactMatrixTRtranATLcon, dimension of left and right matrix not compatible\n");
    return 1;
  }
  TRtranATLcon->Ni=TR->Nj;
  TRtranATLcon->Nj=TL->Nj;
  TRtranATLcon->QComplex = 1;
  InitCompactMatrix(TRtranATLcon);
  TATindex=0;
  if(TL->QComplex)
  {
    if(TR->QComplex)
    { // L and R complex
      for(i=0;i<TRtranATLcon->Ni;i++)
      {
        for(j=0;j<TRtranATLcon->Nj;j++)
        {
          irow=i;
          jrow=j;
          TRtranATLcon->V[TATindex].Re=0;
          TRtranATLcon->V[TATindex].Im=0;
          for(k=0;k<TRtranATLcon->Nj;k++)
          {
            TRtranATLcon->V[TATindex].Re+=TL->V[jrow].Re*A[k].Re*TR->V[irow].Re+TL->V[jrow].Im*A[k].Im*TR->V[irow].Re
                                         +TL->V[jrow].Im*A[k].Re*TR->V[irow].Im-TL->V[jrow].Re*A[k].Im*TR->V[irow].Im;
            TRtranATLcon->V[TATindex].Im+=TL->V[jrow].Re*A[k].Im*TR->V[irow].Re+TL->V[jrow].Im*A[k].Im*TR->V[irow].Im
                                         +TL->V[jrow].Re*A[k].Re*TR->V[irow].Im-TL->V[jrow].Im*A[k].Re*TR->V[irow].Re;
            irow+=TRtranATLcon->Nj;
            jrow+=TRtranATLcon->Nj;
          }
          TATindex++;
        }
      }
    }
    else
    { // L complex, R real
      for(i=0;i<TRtranATLcon->Ni;i++)
      {
        for(j=0;j<TRtranATLcon->Nj;j++)
        {
          irow=i;
          jrow=j;
          TRtranATLcon->V[TATindex].Re=0;
          TRtranATLcon->V[TATindex].Im=0;
          for(k=0;k<TRtranATLcon->Nj;k++)
          {
            TRtranATLcon->V[TATindex].Re+=TL->V[jrow].Re*A[k].Re*TR->VRe[irow]+TL->V[jrow].Im*A[k].Im*TR->VRe[irow];
            TRtranATLcon->V[TATindex].Im+=TL->V[jrow].Re*A[k].Im*TR->VRe[irow]-TL->V[jrow].Im*A[k].Re*TR->VRe[irow];
            irow+=TRtranATLcon->Nj;
            jrow+=TRtranATLcon->Nj;
          }
          TATindex++;
        }
      }
    }
  }
  else
  {
    if(TR->QComplex)
    { // L real, R complex
      for(i=0;i<TRtranATLcon->Ni;i++)
      {
        for(j=0;j<TRtranATLcon->Nj;j++)
        {
          irow=i;
          jrow=j;
          TRtranATLcon->V[TATindex].Re=0;
          TRtranATLcon->V[TATindex].Im=0;
          for(k=0;k<TRtranATLcon->Nj;k++)
          {
            TRtranATLcon->V[TATindex].Re+=TL->VRe[jrow]*A[k].Re*TR->V[irow].Im-TL->VRe[jrow]*A[k].Im*TR->V[irow].Im;
            TRtranATLcon->V[TATindex].Im+=TL->VRe[jrow]*A[k].Im*TR->V[irow].Re+TL->VRe[jrow]*A[k].Re*TR->V[irow].Im;
            irow+=TRtranATLcon->Nj;
            jrow+=TRtranATLcon->Nj;
          }
          TATindex++;
        }
      }
    }
    else
    { //L and R real
      for(i=0;i<TRtranATLcon->Ni;i++)
      {
        for(j=0;j<TRtranATLcon->Nj;j++)
        {
          irow=i;
          jrow=j;
          TRtranATLcon->V[TATindex].Re=0;
          TRtranATLcon->V[TATindex].Im=0;
          for(k=0;k<TRtranATLcon->Nj;k++)
          {
            TRtranATLcon->V[TATindex].Re+=TL->VRe[jrow]*A[k].Re*TR->VRe[irow];
            TRtranATLcon->V[TATindex].Im+=TL->VRe[jrow]*A[k].Im*TR->VRe[irow];
            irow+=TRtranATLcon->Nj;
            jrow+=TRtranATLcon->Nj;
          }
          TATindex++;
        }
      }
    }
  }
  return 0;
}

/*   C = beta*C + alpha*(A ^{TransA} * B ^{TransB})
 *   TransA can be either:
 *          N for nothing
 *          T for transposed
 *          C for conjugate transposed
 *          * for conjugate
 */
uint32_t CompactMatrixAddMultiply(char const TransA, char const TransB, CompactMatrixType const *A, CompactMatrixType const *B, CompactMatrixType *C, double alphaRe, double alphaIm, double betaRe, double betaIm)
{
  if((C->Ni==0)||(C->Nj==0))return 0;
  CompactMatrixType *Ap,*Bp,*Cp;
  char TransAp=TransA, TransBp=TransB;
  uint32_t Nk,LDA,LDB,LDC,i;
  Complex16 ALPHA,BETA;
  ALPHA.Re=alphaRe; ALPHA.Im=alphaIm;
  BETA.Re =betaRe;  BETA.Im =betaIm;
  LDA=A->Nj;
  LDB=B->Nj;
  LDC=C->Nj;
  // check if C should be real or complex and make complex if needed
  if((C->QComplex)&&(fabs(betaIm)<100*DBL_MIN)&&(fabs(betaRe)<100*DBL_MIN)&&(A->QComplex==0)&&(B->QComplex==0))
  {
    MakeCompactComplexMatrixReal(C);
  }
  if(((A->QComplex)||(B->QComplex)||(fabs(alphaIm)>100*DBL_MIN)||(fabs(betaIm)>100*DBL_MIN))&&(C->QComplex==0))
  {
    if(MakeCompactRealMatrixComplex(C)){printf("Failed to allocate memory in CompactMatrixAddMultiply\n");return 1;};
  }
  // make Ap, Bp, and Cp point to matrices of the same type as C (real or complex)
  Ap=(CompactMatrixType*)A;
  Bp=(CompactMatrixType*)B;
  Cp=C;
  if((Cp->QComplex)&&(A->QComplex==0))
  {
    Ap= (CompactMatrixType*)malloc(sizeof(CompactMatrixType));
    Ap->Ni=A->Ni;
    Ap->Nj=A->Nj;
    Ap->QComplex=1;
    InitCompactMatrix(Ap);
    for(i=0;i<A->Ni*A->Nj;i++)
    {
      Ap->V[i].Re=A->VRe[i];
    }
  }
  if((Cp->QComplex)&&(B->QComplex==0))
  {
    Bp= (CompactMatrixType*)malloc(sizeof(CompactMatrixType));
    Bp->Ni=B->Ni;
    Bp->Nj=B->Nj;
    Bp->QComplex=1;
    InitCompactMatrix(Bp);
    for(i=0;i<B->Ni*B->Nj;i++)
    {
      Bp->V[i].Re=B->VRe[i];
    }
  }
  // replace * by N and conjugate Ap or Bp
  if(TransAp=='*')
  {
    TransAp='N';
    CompactMatrixConjugate(Ap);
  }
  if(TransBp=='*')
  {
    TransBp='N';
    CompactMatrixConjugate(Bp);
  }
  // test if dimensions are ok and set Nk
  if((TransA=='N')||(TransA=='*'))
  {
    if((TransB=='N')||(TransB=='*'))
    {
      Nk=A->Nj;
      if((A->Nj!=B->Ni)||(C->Ni!=A->Ni)||(C->Nj!=B->Nj))
      {
        printf("Trying to multiply matrix (%s)^%c with dimension %i %i with matrix (%s)^%c with dimension %i %i and add it to matrix (%s) with dimension %i %i\n",A->Name,TransA,A->Ni,A->Nj,B->Name,TransB,B->Ni,B->Nj,C->Name,C->Ni,C->Nj);
        printf("Error in dimensions\n");
        return 1;
      }
    }
    else
    {
      Nk=A->Nj;
      if((A->Nj!=B->Nj)||(C->Ni!=A->Ni)||(C->Nj!=B->Ni))
      {
        printf("Trying to multiply matrix (%s)^%c with dimension %i %i with matrix (%s)^%c with dimension %i %i and add it to matrix (%s) with dimension %i %i\n",A->Name,TransA,A->Ni,A->Nj,B->Name,TransB,B->Ni,B->Nj,C->Name,C->Ni,C->Nj);
        printf("Error in dimensions\n");
        return 1;
      }
    }
  }
  else
  {
    if((TransB=='N')||(TransB=='*'))
    {
      Nk=A->Ni;
      if((A->Ni!=B->Ni)||(C->Ni!=A->Nj)||(C->Nj!=B->Nj))
      {
        printf("Trying to multiply matrix (%s)^%c with dimension %i %i with matrix (%s)^%c with dimension %i %i and add it to matrix (%s) with dimension %i %i\n",A->Name,TransA,A->Ni,A->Nj,B->Name,TransB,B->Ni,B->Nj,C->Name,C->Ni,C->Nj);
        printf("Error in dimensions\n");
        return 1;
      }
    }
    else
    {
      Nk=A->Ni;
      if((A->Ni!=B->Nj)||(C->Ni!=A->Nj)||(C->Nj!=B->Ni))
      {
        printf("Trying to multiply matrix (%s)^%c with dimension %i %i with matrix (%s)^%c with dimension %i %i and add it to matrix (%s) with dimension %i %i\n",A->Name,TransA,A->Ni,A->Nj,B->Name,TransB,B->Ni,B->Nj,C->Name,C->Ni,C->Nj);
        printf("Error in dimensions\n");
        return 1;
      }
    }
  }
  // do the multiplication
  if(Cp->QComplex)
  {
    zgemm(&TransBp,&TransAp,&(Cp->Nj),&(Cp->Ni),&Nk,&(ALPHA)   ,Bp->V  ,&LDB,Ap->V  ,&LDA,&(BETA)   ,Cp->V  ,&LDC);
  }
  else
  {
    dgemm(&TransBp,&TransAp,&(Cp->Nj),&(Cp->Ni),&Nk,&(ALPHA.Re),Bp->VRe,&LDB,Ap->VRe,&LDA,&(BETA.Re),Cp->VRe,&LDC);
  }
  // replace * by N and conjugate Ap or Bp
  if((TransA=='*')&&(A->QComplex))
  {
    CompactMatrixConjugate(Ap);
  }
  if((TransB=='*')&&(B->QComplex))
  {
    CompactMatrixConjugate(Bp);
  }
  // free Ap and Bp if needed
  if((Cp->QComplex)&&(A->QComplex==0))
  {
    FreeCompactMatrix(Ap);
    free(Ap);
  }
  if((Cp->QComplex)&&(B->QComplex==0))
  {
    FreeCompactMatrix(Bp);
    free(Bp);
  }
  return 0;
}

/*   C = A ^{TransA} * B ^{TransB}
 *   C should be allocated before the routine is called
 *   TransA can be either:
 *          N for nothing
 *          T for transposed
 *          C for conjugate transposed
 *          * for conjugate
 */
uint32_t CompactMatrixMultiplyPreAllocated(char const TransA, char const TransB, CompactMatrixType const *A, CompactMatrixType const *B, CompactMatrixType *C)
{
  return CompactMatrixAddMultiply(TransA, TransB, A, B, C, 1, 0, 0, 0);
}

/*   C = A ^{TransA} * B ^{TransB}
 *   C should not be allocated before the routine is called
 *   TransA can be either:
 *          N for nothing
 *          T for transposed
 *          C for conjugate transposed
 *          * for conjugate
 */
uint32_t CompactMatrixMultiply(char const TransA, char const TransB, CompactMatrixType const *A, CompactMatrixType const *B, CompactMatrixType *C)
{
  C->QComplex=0;
  if((A->QComplex)||(B->QComplex))
  {
    C->QComplex=1;
  }
  if((TransA=='N')||(TransA=='*'))
  {
    C->Ni=A->Ni;
  }
  else
  {
    C->Ni=A->Nj;
  }
  if((TransB=='N')||(TransB=='*'))
  {
    C->Nj=B->Nj;
  }
  else
  {
    C->Nj=B->Ni;
  }
  if(InitCompactMatrix(C)){printf("Could not initialize matrix in CompactMatrixMultiply\n");return 1;};
  return CompactMatrixAddMultiply(TransA, TransB, A, B, C, 1, 0, 0, 0);
}

/*   D = A ^{TransA} * B ^{TransB} * C ^{TransC}
 *   D should not be allocated before the routine is called
 *   TransA can be either:
 *          N for nothing
 *          T for transposed
 *          C for conjugate transposed
 *          * for conjugate
 */
uint32_t CompactMatrixMultiply3(char const TransA, char const TransB, char const TransC, CompactMatrixType const *A, CompactMatrixType const *B, CompactMatrixType const *C, CompactMatrixType *D)
{
  CompactMatrixType temp;
  if(CompactMatrixMultiply(TransB, TransC, B, C, &temp)){printf("Error in CompactMatrixMultiply3\n"); return 1;}
  if(CompactMatrixMultiply(TransA, 'N', A, &temp, D)){printf("Error in CompactMatrixMultiply3\n"); return 1;}
  FreeCompactMatrix(&temp);
  return 0;
}
uint32_t CompactMatrixMultiply4(char const TransA, char const TransB, char const TransC, char const TransD, CompactMatrixType const *A, CompactMatrixType const *B, CompactMatrixType const *C, CompactMatrixType const *D, CompactMatrixType *E)
{
  CompactMatrixType temp1,temp2;
  if(CompactMatrixMultiply(TransA, TransB, A, B, &temp1)){printf("Error in CompactMatrixMultiply4\n"); return 1;}
  if(CompactMatrixMultiply(TransC, TransD, C, D, &temp2)){printf("Error in CompactMatrixMultiply4\n"); return 1;}
  if(CompactMatrixMultiply('N', 'N', &temp1, &temp2, E)) {printf("Error in CompactMatrixMultiply4\n"); return 1;}
  FreeCompactMatrix(&temp1);
  FreeCompactMatrix(&temp2);
  return 0;
}
uint32_t CompactMatrixMultiply3PreAllocated(char const TransA, char const TransB, char const TransC, CompactMatrixType const *A, CompactMatrixType const *B, CompactMatrixType const *C, CompactMatrixType *D)
{
  CompactMatrixType temp;
  if(CompactMatrixMultiply(TransB, TransC, B, C, &temp)){printf("Error in CompactMatrixMultiply3\n"); return 1;}
  if(CompactMatrixMultiplyPreAllocated(TransA, 'N', A, &temp, D)){printf("Error in CompactMatrixMultiply3\n"); return 1;}
  FreeCompactMatrix(&temp);
  return 0;
}

/*  Calculates TdagAT = T^dag A T
 *
 */
uint32_t TransformCompactMatrixCN(CompactMatrixType const *A, CompactMatrixType const *T, CompactMatrixType * const TdagAT)
{
  return CompactMatrixMultiply3('C', 'N', 'N', T, A, T, TdagAT);
}
/*  Calculates TdagAT = T^dag A T
 *
 */
uint32_t TransformCompactMatrixPreAllocatedCN(CompactMatrixType const *A, CompactMatrixType const *T, CompactMatrixType * const TdagAT)
{
  return CompactMatrixMultiply3PreAllocated('C', 'N', 'N', T, A, T, TdagAT);
}
/*  Calculates A <= T^dag A T
 *
 */
uint32_t TransformCompactMatrixInlineCN(CompactMatrixType * const A, CompactMatrixType const *T)
{
  CompactMatrixType temp;
  if(CompactMatrixMultiply('N', 'N', A, T, &temp)){printf("Error in TransformCompactMatrixInline\n"); return 1;}
  if(CompactMatrixMultiplyPreAllocated('C', 'N', T, &temp, A)){printf("Error in TransformCompactMatrixInline\n"); return 1;}
  FreeCompactMatrix(&temp);
  return 0;
}

/*  Calculates TATdag = T A T^dag
 *
 */
uint32_t TransformCompactMatrixNC(CompactMatrixType const *A, CompactMatrixType const *T, CompactMatrixType * const TATdag)
{
  return CompactMatrixMultiply3('N', 'N', 'C', T, A, T, TATdag);
}
/*  Calculates TdagAT = T^dag A T
 *
 */
uint32_t TransformCompactMatrixPreAllocatedNC(CompactMatrixType const *A, CompactMatrixType const *T, CompactMatrixType * const TATdag)
{
  return CompactMatrixMultiply3PreAllocated('N', 'N', 'C', T, A, T, TATdag);
}
/*  Calculates A <= T A T^dag
 *
 */
uint32_t TransformCompactMatrixInlineNC(CompactMatrixType * const A, CompactMatrixType const *T)
{
  CompactMatrixType temp;
  if(CompactMatrixMultiply('N', 'C', A, T, &temp)){printf("Error in TransformCompactMatrixInlineNNC\n"); return 1;}
  if(CompactMatrixMultiplyPreAllocated('N', 'N', T, &temp, A)){printf("Error in TransformCompactMatrixInlineNNC\n"); return 1;}
  FreeCompactMatrix(&temp);
  return 0;
}

/*  Calculates TTATc = T^T A T^c
 *
 */
uint32_t TransformCompactMatrixTc(CompactMatrixType const *A, CompactMatrixType const *T, CompactMatrixType * const TTATc)
{
  return CompactMatrixMultiply3('T', 'N', '*', T, A, T, TTATc);
}
/*  Calculates TTATc = T^T A T^c
 *
 */
uint32_t TransformCompactMatrixPreAllocatedTc(CompactMatrixType const *A, CompactMatrixType const *T, CompactMatrixType * const TTATc)
{
  return CompactMatrixMultiply3PreAllocated('T', 'N', '*', T, A, T, TTATc);
}
/*  Calculates A <= T^T A T^c
 *
 */
uint32_t TransformCompactMatrixInlineTc(CompactMatrixType * const A, CompactMatrixType const *T)
{
  CompactMatrixType temp;
  if(CompactMatrixMultiply('N', '*', A, T, &temp)){printf("Error in TransformCompactMatrixInline\n"); return 1;}
  if(CompactMatrixMultiplyPreAllocated('T', 'N', T, &temp, A)){printf("Error in TransformCompactMatrixInline\n"); return 1;}
  FreeCompactMatrix(&temp);
  return 0;
}

/*  Calculates TcATT = T^c A T^T
 *
 */
uint32_t TransformCompactMatrixcT(CompactMatrixType const *A, CompactMatrixType const *T, CompactMatrixType * const TcATT)
{
  return CompactMatrixMultiply3('*', 'N', 'T', T, A, T, TcATT);
}
/*  Calculates TdagAT = T^dag A T
 *
 */
uint32_t TransformCompactMatrixPreAllocatedcT(CompactMatrixType const *A, CompactMatrixType const *T, CompactMatrixType * const TcATT)
{
  return CompactMatrixMultiply3PreAllocated('*', 'N', 'T', T, A, T, TcATT);
}
/*  Calculates A <= T^c A T^dag
 *
 */
uint32_t TransformCompactMatrixInlinecT(CompactMatrixType * const A, CompactMatrixType const *T)
{
  CompactMatrixType temp;
  if(CompactMatrixMultiply('N', 'T', A, T, &temp)){printf("Error in TransformCompactMatrixInlineNNC\n"); return 1;}
  if(CompactMatrixMultiplyPreAllocated('*', 'N', T, &temp, A)){printf("Error in TransformCompactMatrixInlineNNC\n"); return 1;}
  FreeCompactMatrix(&temp);
  return 0;
}

/*  For matrix fun
 *    sets all rows except row i to zero and calculates
 *    Out = Out + fun^T . fun^*
 */
int FunTrDotFunCCPreAllocated(CompactMatrixType const * const fun, int i, CompactMatrixType * Out)
{
  int j1,j2;
  if(fun->QComplex==0)
  {
    for(j1=0;j1<fun->Nj;j1++)
    {
      for(j2=0;j2<fun->Nj;j2++)
      {
        Out->VRe[j2 + j1*(fun->Nj)] += fun->VRe[j1 + i * (fun->Nj)] * fun->VRe[j2 + i * (fun->Nj)];
      }
    }
  }
  else
  {
    for(j1=0;j1<fun->Nj;j1++)
    {
      for(j2=0;j2<fun->Nj;j2++)
      {
        Out->V[j2 + j1*(fun->Nj)].Re += fun->V[j1 + i * (fun->Nj)].Re * fun->V[j2 + i * (fun->Nj)].Re + fun->V[j1 + i * (fun->Nj)].Im * fun->V[j2 + i * (fun->Nj)].Im;
        Out->V[j2 + j1*(fun->Nj)].Im += fun->V[j1 + i * (fun->Nj)].Im * fun->V[j2 + i * (fun->Nj)].Re - fun->V[j1 + i * (fun->Nj)].Re * fun->V[j2 + i * (fun->Nj)].Im;
      }
    }
  }
  return 0;
}
/*  For matrix fun
 *    sets all rows except row i to zero and calculates
 *    Out =  fun^T . fun^*
 */
int FunTrDotFunCC(CompactMatrixType const * const fun, int i, CompactMatrixType * Out)
{
  Out->Ni=fun->Nj;
  Out->Nj=fun->Nj;
  Out->QComplex = fun->QComplex;
  InitCompactMatrix(Out);
  return FunTrDotFunCCPreAllocated(fun, i, Out);
}




/* M=M^{T}
 */
uint32_t CompactMatrixTranspose(CompactMatrixType * const M)
{
  uint32_t i,j,iIndex=0,jIndex=0;
  if(M->Ni==M->Nj)
  {
    if(M->QComplex) /* transpose a complex square matrix */
    {
      Complex16 swap;
      for(i=0;i<M->Ni;i++)
      {
        for(j=0;j<M->Nj;j++)
        {
          if(j>i)
          {
            swap=M->V[jIndex];
            M->V[jIndex]=M->V[iIndex];
            M->V[iIndex]=swap;
          }
          jIndex++;
          iIndex+=M->Ni;
        }
        iIndex+=1;
        iIndex-=M->Ni*M->Nj;
      }
    }
    else /* transpose a real square matrix */
    {
      double swap;
      for(i=0;i<M->Ni;i++)
      {
        for(j=0;j<M->Nj;j++)
        {
          if(j>i)
          {
            swap=M->VRe[jIndex];
            M->VRe[jIndex]=M->VRe[iIndex];
            M->VRe[iIndex]=swap;
          }
          jIndex++;
          iIndex+=M->Ni;
        }
        iIndex+=1;
        iIndex-=M->Ni*M->Nj;
      }
    }
  }
  else
  {
    if(M->QComplex) /* transpose a complex rectangular matrix */
    {
      Complex16 *swap;
      swap= (Complex16*)malloc(M->Ni*M->Nj*sizeof(Complex16));
      if(swap==NULL)
      {
        printf("malloc failed in CompactMatrixTranspose\n");
        return 1;
      }
      for(i=0;i<M->Ni;i++)
      {
        for(j=0;j<M->Nj;j++)
        {
          swap[iIndex]=M->V[jIndex];
          jIndex++;
          iIndex+=M->Ni;
        }
        iIndex+=1;
        iIndex-=M->Ni*M->Nj;
      }
      free(M->V);
      M->V=swap;
      i=M->Ni;
      M->Ni=M->Nj;
      M->Nj=i;
    }
    else /* transpose a real rectangular matrix */
    {
      double *swap;
      swap= (double*)malloc(M->Ni*M->Nj*sizeof(double));
      if(swap==NULL)
      {
        printf("malloc failed in CompactMatrixTranspose\n");
        return 1;
      }
      for(i=0;i<M->Ni;i++)
      {
        for(j=0;j<M->Nj;j++)
        {
          swap[iIndex]=M->VRe[jIndex];
          jIndex++;
          iIndex+=M->Ni;
        }
        iIndex+=1;
        iIndex-=M->Ni*M->Nj;
      }
      free(M->VRe);
      M->VRe=swap;
      i=M->Ni;
      M->Ni=M->Nj;
      M->Nj=i;
    }
  }
  return 0;
}

/* M=M^{*}
 */
uint32_t CompactMatrixConjugate(CompactMatrixType * const M)
{
  uint32_t i;
  if(M->QComplex)
  {
    for(i=0;i<M->Ni*M->Nj;i++)
    {
      M->V[i].Im*=-1;
    }
  }
  return 0;
}

/* M=M^{+}
 */
uint32_t CompactMatrixConjugateTranspose(CompactMatrixType * const M)
{
  CompactMatrixConjugate(M);
  CompactMatrixTranspose(M);
  return 0;
}

int SquareCompactMatrixMakeHermitian(CompactMatrixType * const M)
{
  if(M->Ni!=M->Nj) return 1;
  if(M->QComplex==0)
  {
    for(int i=0;i<M->Ni;i++)
    {
      for(int j=i+1;j<M->Ni;j++)
      {
        M->VRe[j + i * M->Nj] += M->VRe[i + j * M->Nj];
        M->VRe[j + i * M->Nj] /=  2;
        M->VRe[i + j * M->Nj] = M->VRe[j + i * M->Nj];
      }
    }
  }
  else
  {
    for(int i=0;i<M->Ni;i++)
    {
      M->V[i + i * M->Nj].Im = 0;
      for(int j=i+1;j<M->Ni;j++)
      {
        M->V[j + i * M->Nj].Re += M->V[i + j * M->Nj].Re;
        M->V[j + i * M->Nj].Im -= M->V[i + j * M->Nj].Im;
        M->V[j + i * M->Nj] /=  2;
        M->V[i + j * M->Nj].Re = M->V[j + i * M->Nj].Re;
        M->V[i + j * M->Nj].Im =-M->V[j + i * M->Nj].Im;
      }
    }
  }
  return 0;
}

int SquareCompactMatrixMakeAntiHermitian(CompactMatrixType * const M)
{
  if(M->Ni!=M->Nj) return 1;
  if(M->QComplex==0)
  {
    for(int i=0;i<M->Ni;i++)
    {
      M->VRe[i + i * M->Nj] = 0;
      for(int j=i+1;j<M->Ni;j++)
      {
        M->VRe[j + i * M->Nj] -= M->VRe[i + j * M->Nj];
        M->VRe[j + i * M->Nj] /=  2;
        M->VRe[i + j * M->Nj] =-M->VRe[j + i * M->Nj];
      }
    }
  }
  else
  {
    for(int i=0;i<M->Ni;i++)
    {
      M->V[i + i * M->Nj].Re = 0;
      for(int j=i+1;j<M->Ni;j++)
      {
        M->V[j + i * M->Nj].Re -= M->V[i + j * M->Nj].Re;
        M->V[j + i * M->Nj].Im += M->V[i + j * M->Nj].Im;
        M->V[j + i * M->Nj] /=  2;
        M->V[i + j * M->Nj].Re =-M->V[j + i * M->Nj].Re;
        M->V[i + j * M->Nj].Im = M->V[j + i * M->Nj].Im;
      }
    }
  }
  return 0;
}

/* M=M^{-1}
 */
double CompactMatrixInverse(CompactMatrixType * const M)
{
  int *IPIV;
  uint32_t LWORK = M->Ni*M->Nj;
  int INFO;
  if(M->Ni!=M->Nj)
  {
    printf("In CompactMatrixInverse tried to inverse a rectangular matrix Ni=%i, Nj=%i\n",M->Ni,M->Nj);
    return -1;
  }
  IPIV = (int*)malloc((M->Ni+1)*sizeof(int));
  if(M->QComplex)
  {
    Complex16 *WORK;
    zgetrf(&M->Ni,&M->Ni,M->V,&M->Ni,IPIV,&INFO);
    if(INFO!=0)
    {
      printf("zgetrf failed in CompactMatrixInverse with INFO %i\n",INFO);
      return 1;
    }
    WORK = (Complex16*)malloc(LWORK*sizeof(Complex16));
    zgetri(&M->Ni,M->V,&M->Ni,IPIV,WORK,&LWORK,&INFO);
    if(INFO!=0)
    {
      printf("zgetri failed in CompactMatrixInverse with INFO %i\n",INFO);
      free(WORK);
      return 1;
    }
    free(WORK);
  }
  else
  {
    double *WORK;
    dgetrf(&M->Ni,&M->Ni,M->VRe,&M->Ni,IPIV,&INFO);
    if(INFO!=0)
    {
      printf("dgetrf failed in CompactMatrixInverse with INFO %i\n",INFO);
      return 1;
    }
    WORK = (double*)malloc(LWORK*sizeof(double));
    dgetri(&M->Ni,M->VRe,&M->Ni,IPIV,WORK,&LWORK,&INFO);
    if(INFO!=0)
    {
      printf("dgetri failed in CompactMatrixInverse with INFO %i\n",INFO);
      free(WORK);
      return 1;
    }
    free(WORK);
  }
  free(IPIV);
  return INFO;
}

/*
 *
 */
uint32_t CompactMatrixScale(CompactMatrixType * const A, double alpha)
{
  uint32_t i;
  if(A->QComplex)
  {
    for(i=0;i<A->Ni*A->Nj;i++)
    {
      A->V[i].Re*=alpha;
      A->V[i].Im*=alpha;
    }
  }
  else
  {
    for(i=0;i<A->Ni*A->Nj;i++)
    {
      A->VRe[i]*=alpha;
    }
  }
  return 0;
}

/*
 *
 */
uint32_t CompactMatrixScaleComplex(CompactMatrixType * const A, Complex16 alpha)
{
  uint32_t i;
  Complex16 v;
  if(alpha.Im==0) return CompactMatrixScale(A, alpha.Re);
  MakeCompactRealMatrixComplex(A);
  for(i=0;i<A->Ni*A->Nj;i++)
  {
    v.Re = A->V[i].Re*alpha.Re - A->V[i].Im*alpha.Im;
    v.Im = A->V[i].Re*alpha.Im + A->V[i].Im*alpha.Re;
    A->V[i].Re=v.Re;
    A->V[i].Im=v.Im;
  }
  return 0;
}

/*  A <- A + B
 *
 */
uint32_t CompactMatrixAdd(CompactMatrixType * const A, CompactMatrixType const * B)
{
  int i;
  if(((A->Ni)!=(B->Ni))||((A->Nj)!=(B->Nj)))
  {
    printf("Error in CompactMatrixAdd can not add matrices with different dimensions\n");
    return 1;
  }
  if(A->QComplex)
  {
    if(B->QComplex)
    {
      for(i=0;i<A->Ni*A->Nj;i++)
      {
        A->V[i].Re+=B->V[i].Re;
        A->V[i].Im+=B->V[i].Im;
      }
    }
    else
    {
      for(i=0;i<A->Ni*A->Nj;i++)
      {
        A->V[i].Re+=B->VRe[i];
      }
    }
  }
  else
  {
    if(B->QComplex)
    {
      MakeCompactRealMatrixComplex(A);
      for(i=0;i<A->Ni*A->Nj;i++)
      {
        A->V[i].Re+=B->V[i].Re;
        A->V[i].Im+=B->V[i].Im;
      }
    }
    else
    {
      for(i=0;i<A->Ni*A->Nj;i++)
      {
        A->VRe[i]+=B->VRe[i];
      }
    }
  }
  return 0;
}

/*  A <- A - B
 *
 */
uint32_t CompactMatrixSub(CompactMatrixType * const A, CompactMatrixType const * B)
{
  int i;
  if(((A->Ni)!=(B->Ni))||((A->Nj)!=(B->Nj)))
  {
    printf("Error in CompactMatrixSub can not subtract matrices with different dimensions\n");
    return 1;
  }
  if(A->QComplex)
  {
    if(B->QComplex)
    {
      for(i=0;i<A->Ni*A->Nj;i++)
      {
        A->V[i].Re-=B->V[i].Re;
        A->V[i].Im-=B->V[i].Im;
      }
    }
    else
    {
      for(i=0;i<A->Ni*A->Nj;i++)
      {
        A->V[i].Re-=B->VRe[i];
      }
    }
  }
  else
  {
    if(B->QComplex)
    {
      MakeCompactRealMatrixComplex(A);
      for(i=0;i<A->Ni*A->Nj;i++)
      {
        A->V[i].Re-=B->V[i].Re;
        A->V[i].Im-=B->V[i].Im;
      }
    }
    else
    {
      for(i=0;i<A->Ni*A->Nj;i++)
      {
        A->VRe[i]-=B->VRe[i];
      }
    }
  }
  return 0;
}


double RealMatrixTrace(const CompactMatrixType *M)
{
  double trace = 0;
  const int N = Min(M->Ni,M->Nj);
  for (int i = 0; i < N; ++i)
    trace += M->VRe[i * M->Nj + i];
  return trace;
}

Complex16 ComplexMatrixTrace(const CompactMatrixType *M)
{
  Complex16 trace;
  trace.Re = 0;
  trace.Im = 0;
  const int N = Min(M->Ni,M->Nj);
  for (int i = 0; i < N; ++i)
  {
    trace.Re += M->V[i * M->Nj + i].Re;
    trace.Im += M->V[i * M->Nj + i].Im;
  }
  return trace;
}

Complex16 CompactMatrixTrace(const CompactMatrixType *M)
{
  if(M->QComplex==0)
  {
    Complex16 trace;
    trace.Re = RealMatrixTrace(M);
    trace.Im = 0;
    return trace;
  }
  else
  {
    return ComplexMatrixTrace(M);
  }
}


/* M[i ... i+B.Ni][j ... j+B.Nj] = M[i ... i+B.Ni][j ... j+B.Nj] + alpha B
 */
uint32_t CompactMatrixAddShifted(CompactMatrixType * const M, const CompactMatrixType B, double const alphaRe, double const alphaIm, int ishift, int jshift)
{
  if ((M->Ni < B.Ni + ishift) || (M->Nj < B.Nj + jshift))
  {
    printf("Matrix B shifted does not fit in matrix M\n");
    printf("M.Ni= %i M.Nj= %i B.Ni= %i B.Nj= %i ishift= %i jshift= %i\n", M->Ni, M->Nj, B.Ni, B.Nj, ishift, jshift);
    return 1;
  }
  int    i, j;
  if ((M->QComplex == 0) && (B.QComplex == 0) && (alphaIm == 0))
  {                             /* all real */
    for (i = 0; i < B.Ni; i++)
    {
      for (j = 0; j < B.Nj; j++)
      {
        M->VRe[(i + ishift)*M->Nj+ j + jshift] += alphaRe * B.VRe[i*B.Nj+j];
      }
    }
    return 0;
  }
  if (M->QComplex == 0)
    if (MakeCompactRealMatrixComplex(M))
      return 1;
  if (B.QComplex == 0)
  {                             /* B real, rest complex */
    for (i = 0; i < B.Ni; i++)
    {
      for (j = 0; j < B.Nj; j++)
      {
        M->V[(i + ishift)*M->Nj+j + jshift].Re += alphaRe * B.VRe[i*B.Nj+j];
        M->V[(i + ishift)*M->Nj+j + jshift].Im += alphaIm * B.VRe[i*B.Nj+j];
      }
    }
    return 0;
  }
  /* all complex */
  for (i = 0; i < B.Ni; i++)
  {
    for (j = 0; j < B.Nj; j++)
    {
      M->V[(i + ishift)*M->Nj+j + jshift].Re += alphaRe * B.V[i*B.Nj+j].Re - alphaIm * B.V[i*B.Nj+j].Im;
      M->V[(i + ishift)*M->Nj+j + jshift].Im += alphaIm * B.V[i*B.Nj+j].Re + alphaRe * B.V[i*B.Nj+j].Im;
    }
  }
  return 0;
}

/*   Matrix A has Ni rows and Nj coulomns with Ni <= Nj
 *   We assume that A^* . A^T = I with I an identity matrix of dimension Ni by Ni
 *   We return matrix NA such that {A,NA} is a square matrix of dimenstion Nj by Nj and {A,NA}^* . {A,NA}^T = I
 */
uint32_t CompactMatrixNullSpace(CompactMatrixType * Ain, CompactMatrixType * NA)
{
  char const *JOBU ="A";  /* in SVD A = U S V^T return full array U */
  char const *JOBVT="N";  /* in SVD A = U S V^T return nothing for VT */
  CompactMatrixType A;
  CompactMatrixCopy(*Ain, &A);
  Complex16 ALPHA,BETA;
  ALPHA.Re=1; ALPHA.Im=0;
  BETA.Re =0; BETA.Im =0;
  uint32_t LWORK;
  int INFO;
  LWORK = Max(5,A.Ni)*Max(5,A.Nj);
  double *WORK;
  Complex16 *WORKIM;
  double *V=NULL;
  Complex16 *Vc=NULL;
  double *S=(double*)calloc(A.Ni,sizeof(double));
  NA->Ni=A.Nj;
  NA->Nj=A.Nj;
  NA->QComplex=A.QComplex;
  InitCompactMatrix(NA);
  if(A.QComplex)
  {
    WORK = (double*)malloc(5*Min(A.Ni,A.Nj)*sizeof(double));
    if(WORK==NULL)
    {
      printf("malloc failed in CompactMatrixSVDOrthogonalize\n");
      return 1;
    }
    WORKIM = (Complex16*)malloc(LWORK*sizeof(Complex16));
    if(WORKIM==NULL)
    {
      free(WORK);
      printf("malloc failed in CompactMatrixSVDOrthogonalize\n");
      return 1;
    }
    zgesvd(JOBU, JOBVT, &A.Nj, &A.Ni, A.V, &A.Nj,S,NA->V, &A.Nj,Vc, &A.Nj, WORKIM, &LWORK, WORK, &INFO);
    if(INFO)
    {
      printf("zgesvd failed in CompactMatrixSVDOrthogonalize with INFO %i\n",INFO);
      return 1;
    }
    memmove(NA->VRe,NA->V+(A.Nj*A.Ni),A.Nj*(A.Nj-A.Ni)*sizeof(Complex16));
    NA->Ni-=A.Ni;
  }
  else
  {
    WORK = (double*)malloc(LWORK*sizeof(double));
    if(WORK==NULL)
    {
      printf("malloc failed in CompactMatrixSVDOrthogonalize\n");
      return 1;
    }
    dgesvd(JOBU, JOBVT, &A.Nj, &A.Ni, A.VRe, &A.Nj,S,NA->VRe, &A.Nj,V, &A.Nj, WORK, &LWORK, &INFO);
    if(INFO)
    {
      printf("dgesvd failed in CompactMatrixSVDOrthogonalize with INFO %i (A->Ni %u, A->Nj %u)\n",INFO,A.Ni,A.Nj);
      return 1;
    }
    memmove(NA->VRe,NA->VRe+(A.Nj*A.Ni),A.Nj*(A.Nj-A.Ni)*sizeof(double));
    NA->Ni-=A.Ni;
  }
  free(S);
  return 0;
}


/*  A is replaced by an orthonormal form of A using Singular value decomposition (eigen value decomposition)
 *  If a SingularValue is smaller than SingularValue this vector is removed from the basis
 *  If ReduceDimension is true the number of Rows in A is allowed to change and A . A^dag = 1
 *  If ReduceDimension is false the number of Rows in A stays constant. The matrix returned will not be fully orthonormal
 *  The trace of A.A^dag is equal to the numer of non-singular elements in A
 *  Diagonalizing A.A^dag wil yield NSingular states with eigen-value zero and all others with eigen-value 1.
 */
uint32_t CompactMatrixSVDOrthogonalizeRow(CompactMatrixType * A, const double SingularValue, const bool ReduceDimension)
{
  uint32_t i,j,NSingular;
  char const *JOBU ="O";  /* in SVD A = U S V^T return the first min(Ni,Nj) collumns of U in matrix A, which is over written */
  char const *JOBVT="A";  /* in SVD A = U S V^T return full array VT */
  char const *TRANSA="N";
  char const *TRANSB="N";
  Complex16 ALPHA,BETA;
  ALPHA.Re=1; ALPHA.Im=0;
  BETA.Re =0; BETA.Im =0;
  uint32_t M,N,K,LDA,LDB,LDC,LWORK;
  int INFO;
  LWORK = Max(5,A->Ni)*Max(5,A->Nj);
  double *WORK;
  Complex16 *WORKIM;
  double *U=NULL;
  Complex16 *Uc=NULL;
  double *S;
  CompactMatrixType V;
  if(A->Ni==0) return 0;
  if(A->Nj==0)
  {
    A->Ni=0;
    return 0;
  }
  V.Ni=A->Ni;
  V.Nj=A->Ni;
  S= (double*)calloc(A->Ni,sizeof(double));
  V.QComplex=A->QComplex;
  InitCompactMatrix(&V);
  if(A->QComplex)
  {
    WORK = (double*)malloc(5*Min(A->Ni,A->Nj)*sizeof(double));
    if(WORK==NULL)
    {
      printf("malloc failed in CompactMatrixSVDOrthogonalize\n");
      return 1;
    }
    WORKIM = (Complex16*)malloc(LWORK*sizeof(Complex16));
    if(WORKIM==NULL)
    {
      free(WORK);
      printf("malloc failed in CompactMatrixSVDOrthogonalize\n");
      return 1;
    }
    zgesvd(JOBU, JOBVT, &A->Nj, &A->Ni, A->V, &A->Nj,S,Uc, &A->Ni,V.V, &A->Ni, WORKIM, &LWORK, WORK, &INFO);
    if(INFO)
    {
      printf("zgesvd failed in CompactMatrixSVDOrthogonalize with INFO %i\n",INFO);
      return 1;
    }
    NSingular=0;
    for(i=0;i<A->Ni;i++)
    {
      if(S[i]<=SingularValue)
      {
        S[i]=0;
        NSingular++;
      }
    }
    
    M=A->Nj;//ok
    N=A->Ni;//?
    K=A->Ni-NSingular;//?
    
    LDA=A->Nj;//ok
    LDB=V.Nj; //ok
    LDC=A->Nj;//ok
    
    zgemm(TRANSA,TRANSB,&M,&N,&K,&ALPHA,A->V,&LDA,V.V,&LDB,&BETA,WORKIM,&LDC);
    
    free(A->V);
    A->V=WORKIM;
    FreeCompactMatrix(&V);
    free(WORK);
    
    if((NSingular!=0)&&(ReduceDimension==true))
    {
      uint32_t Spos =0; // position of vector with smallest norm
      double Smin = 0.0,Snow; // norm of the vector at Spos and norm of the next vector we test
      Smin=DBL_MAX;
      for(i=0;i<A->Ni;i++)
      {
        Snow=0;
        for(j=0;j<A->Nj;j++)
        {
          Snow+=A->V[i*A->Nj+j].Re*A->V[i*A->Nj+j].Re+A->V[i*A->Nj+j].Im*A->V[i*A->Nj+j].Im;
        }
        VerbosityPrintf(0x30000,"Determining which vector to remove, norm found at %d is %24.16e \n",i,Snow);
        if (Snow <= Smin)
        {
          Smin = Snow;
          Spos = i;
        }
      }
      VerbosityPrintf(0x10000,"Singular value removed at position %d \n",Spos);
      CompactMatrixRemoveRow(A, Spos);
    }
  }
  else
  {
    WORK = (double*)malloc(LWORK*sizeof(double));
    if(WORK==NULL)
    {
      printf("malloc failed in CompactMatrixSVDOrthogonalize\n");
      return 1;
    }
    dgesvd(JOBU, JOBVT, &A->Nj, &A->Ni, A->VRe, &A->Nj,S,U, &A->Ni,V.VRe, &A->Ni, WORK, &LWORK, &INFO);
    if(INFO)
    {
      printf("dgesvd failed in CompactMatrixSVDOrthogonalize with INFO %i (A->Ni %u, A->Nj %u)\n",INFO,A->Ni,A->Nj);
      return 1;
    }
    NSingular=0;
    for(i=0;i<A->Ni;i++)
    {
      if(S[i]<=SingularValue)
      {
        VerbosityPrintf(0x10000,"Singular value found at position %d with value %e \n",i,S[i]);
        S[i]=0;
        NSingular++;
      }
    }
    
    M=A->Nj;//ok
    N=A->Ni;//?
    K=A->Ni-NSingular;//?
    
    LDA=A->Nj;//ok
    LDB=V.Nj; //ok
    LDC=A->Nj;//ok
    
    dgemm(TRANSA,TRANSB,&M,&N,&K,&ALPHA.Re,A->VRe,&LDA,V.VRe,&LDB,&BETA.Re,WORK,&LDC);
    
    free(A->VRe);
    A->VRe=WORK;
    FreeCompactMatrix(&V);
    
    if((NSingular!=0)&&(ReduceDimension==true))
    {
      uint32_t Spos = 0; // position of vector with smallest norm
      double Smin = 0.0,Snow; // norm of the vector at Spos and norm of the next vector we test
      for(i=0;i<NSingular;i++)
      Smin=DBL_MAX;
      for(i=0;i<A->Ni;i++)
      {
        Snow=0;
        for(j=0;j<A->Nj;j++)
        {
          Snow+=A->VRe[i*A->Nj+j]*A->VRe[i*A->Nj+j];
        }
        VerbosityPrintf(0x30000,"Determining which vector to remove, norm found at %d is %24.16e \n",i,Snow);
        if (Snow <= Smin)
        {
          Smin = Snow;
          Spos = i;
        }
      }
      VerbosityPrintf(0x10000,"Singular value removed at position %d \n",Spos);
      CompactMatrixRemoveRow(A, Spos);
    }
  }
  free(S);
  if((NSingular>0)&&(ReduceDimension==true)) return CompactMatrixSVDOrthogonalizeRow(A, SingularValue, ReduceDimension);
  return 0;
}
uint32_t CompactMatrixSVDOrthogonalizeRowplusT(CompactMatrixType * A, CompactMatrixType *const Anew, CompactMatrixType *const T, CompactMatrixType *const Tinv, double const SingularValue, bool const ReduceDimension)
{
  //TODO: implement (look at git for older implementation)
  return 1;
}



/*  A is replaced by an orthonormal form of A using Gram Schmidt
 *  If the norm of a vector A_{i} is smaller than SingularValue this vector is set to zero or removed from the basis depending on ReduceDimension
 */
uint32_t CompactMatrixGramSchmidtOrthogonalizeRow(CompactMatrixType * A, double const SingularValue, bool const ReduceDimension)
{
  int i,j,ii,ipos,iipos, NSingular;
  double dot,norm,dotRe,dotIm;
  NSingular=0;
  if(A->QComplex)
  {
    for(i=0;i<A->Ni;i++)
    {
      ipos=i*A->Nj;
      /* set vector i-NSingular equal to the normalized vector i */
      norm=0;
      for(j=ipos;j<ipos+A->Nj;j++)
      {
        norm+=Sqr(A->V[j].Re) + Sqr(A->V[j].Im);
      }
      if(sqrt(norm)<=SingularValue)
      {
        norm=0;
        NSingular++;
      }
      else
      {
        norm=1.0/sqrt(norm);
        for(j=ipos;j<ipos+A->Nj;j++)
        {
          A->V[j].Re*=norm;
          A->V[j].Im*=norm;
          /* shift vector to account for singular vectors */
          A->V[j-NSingular*A->Nj].Re=A->V[j].Re;
          A->V[j-NSingular*A->Nj].Im=A->V[j].Im;
        }
        /* for all vectors after vector i subtract projection on vector i */
        for(ii=i+1+NSingular;ii<A->Ni;ii++)
        {
          iipos = ii*A->Nj;
          dotRe=0;
          dotIm=0;
          for(j=0;j<A->Nj;j++)
          {
            dotRe+=A->V[ipos+j].Re*A->V[iipos+j].Re + A->V[ipos+j].Im*A->V[iipos+j].Im;
            dotIm+=A->V[ipos+j].Re*A->V[iipos+j].Im - A->V[ipos+j].Im*A->V[iipos+j].Re;
          }
          for(j=0;j<A->Nj;j++)
          {
            A->V[iipos+j].Re-=dotRe*A->V[ipos+j].Re - dotIm*A->V[ipos+j].Im;
            A->V[iipos+j].Im-=dotRe*A->V[ipos+j].Im + dotIm*A->V[ipos+j].Re;
          }
        }
      }
    }
    if(NSingular!=0)
    {
      if(ReduceDimension==false)
      {
        for(i=(A->Ni-NSingular)*A->Nj;i<A->Ni*A->Nj;i++)
        {
          A->V[i].Re=0;
          A->V[i].Im=0;
        }
      }
      else
      {
        A->Ni-=NSingular;
        A->V=(Complex16*)realloc(A->V,A->Ni*A->Nj*sizeof(Complex16));
      }
    }
  }
  else
  {
    for(i=0;i<A->Ni;i++)
    {
      ipos=i*A->Nj;
      /* set vector i-NSingular equal to the normalized vector i */
      norm=0;
      for(j=ipos;j<ipos+A->Nj;j++)
      {
        norm+=Sqr(A->VRe[j]);
      }
      if(norm<SingularValue)
      {
        norm=0;
        NSingular++;
      }
      else
      {
        norm=1.0/sqrt(norm);
        for(j=ipos;j<ipos+A->Nj;j++)
        {
          A->VRe[j]*=norm;
          /* shift vector to account for singular vectors */
          A->VRe[j-NSingular*A->Nj]=A->VRe[j];
        }
        /* for all vectors subtract projection on vector i */
        for(ii=i+1;ii<A->Ni;ii++)
        {
          iipos = ii*A->Nj;
          dot=0;
          for(j=0;j<A->Nj;j++)
          {
            dot+=A->VRe[ipos+j]*A->VRe[iipos+j];
          }
          for(j=0;j<A->Nj;j++)
          {
            A->VRe[iipos+j]-=dot*A->VRe[ipos+j];
          }
        }
      }
    }
    if(NSingular!=0)
    {
      if(ReduceDimension==false)
      {
        for(i=(A->Ni-NSingular)*A->Nj;i<A->Ni*A->Nj;i++)
        {
          A->VRe[i]=0;
        }
      }
      else
      {
        A->Ni-=NSingular;
        A->VRe=(double*)realloc(A->VRe,A->Ni*A->Nj*sizeof(double));
      }
    }
  }
  return 0;
}

/*  Anew = T A
 *  A = Tinv Anew
 *  with Anew^+ Anew an identity matrix
 */
uint32_t CompactMatrixGramSchmidtOrthogonalizeRowplusT(CompactMatrixType const * A, CompactMatrixType *const Anew, CompactMatrixType *const T, CompactMatrixType *const Tinv, double const SingularValue, bool const ReduceDimension)
{
  //TODO: add ReduceDimension + debug
  uint32_t i,ii,ipos,iipos, one=1;
  CompactMatrixCopy(*A, Anew);
  T->Ni=A->Ni;
  T->Nj=A->Ni;
  T->QComplex=A->QComplex;
  InitCompactIdentityMatrix(T);
  if(A->QComplex)
  {
    Complex16 dot,norm;
    /* loop over all vectors */
    for(i=0;i<A->Ni;i++)
    {
      ipos=i*A->Nj;
      /* subtract all vectors Anew[ii]<i */
      for(ii=0;ii<i;ii++)
      {
        iipos=ii*A->Nj;
        zdotc(&dot,&(A->Nj), &(Anew->V[iipos]), &one,  &(Anew->V[ipos]), &one);
        dot.Re*=-1;
        dot.Im*=-1;
        zaxpy(&(A->Ni),&dot,&(T->V[ii*A->Ni]),&one,&(T->V[i*A->Ni]) ,&one);
        zaxpy(&(A->Nj),&dot,&(Anew->V[iipos]),&one,&(Anew->V[ipos]),&one);
      }
      /* normalize vector i */
      zdotc(&norm,&(A->Nj), &(Anew->V[ipos]), &one,  &(Anew->V[ipos]), &one);
      if(norm.Re<=SingularValue)
      {
        norm.Re=0;
        norm.Im=0;
      }
      else
      {
        norm.Im=0;
        norm.Re=sqrt(1.0/norm.Re);
      }
      zscal(&(A->Nj),&norm,&(Anew->V[ipos]),&one);
      zscal(&(A->Ni),&norm,&(T->V[i*A->Ni]),&one);
    }
  }
  else /* real */
  {
    double dot,norm;
    /* loop over all vectors */
    for(i=0;i<A->Ni;i++)
    {
      ipos=i*A->Nj;
      /* subtract all vectors Anew[ii]<i */
      for(ii=0;ii<i;ii++)
      {
        iipos=ii*A->Nj;
        dot = -ddot(&(A->Nj), &(Anew->VRe[iipos]), &one,  &(Anew->VRe[ipos]), &one);
        daxpy(&(A->Ni),&dot,&(T->VRe[ii*A->Ni]),&one,&(T->VRe[i*A->Ni]) ,&one);
        daxpy(&(A->Nj),&dot,&(Anew->VRe[iipos]),&one,&(Anew->VRe[ipos]),&one);
      }
      /* normalize vector i */
      norm=ddot(&(A->Nj), &(Anew->VRe[ipos]), &one,  &(Anew->VRe[ipos]), &one);
      if(norm<=SingularValue)
      {
        norm=0;
      }
      else
      {
        norm=sqrt(1.0/norm);
      }
      dscal(&(A->Nj),&norm,&(Anew->VRe[ipos]),&one);
      dscal(&(A->Ni),&norm,&(T->VRe[i*A->Ni]),&one);
    }
  }
  CompactMatrixMultiply('N', 'C', A, Anew, Tinv);
  return 0;
}

/*  See the row version for the implementation, here orthogonalization per column
 *  Slow versions using two aditional transposes
 */
uint32_t CompactMatrixSVDOrthogonalizeCol(CompactMatrixType * A, double const SingularValue, bool const ReduceDimension)
{
  if(CompactMatrixTranspose(A)) {printf("CompactMatrixTranspose failed in CompactMatrixSVDOrthogonalizeCol\n");return 1;};
  if(CompactMatrixSVDOrthogonalizeRow(A,SingularValue,ReduceDimension)) {printf("CompactMatrixSVDOrthogonalizeRow failed in CompactMatrixSVDOrthogonalizeCol\n");return 1;};
  if(CompactMatrixTranspose(A)) {printf("CompactMatrixTranspose failed in CompactMatrixSVDOrthogonalizeCol\n");return 1;};
  return 0;
}
uint32_t CompactMatrixSVDOrthogonalizeColplusT(CompactMatrixType const * A, CompactMatrixType *const Anew, CompactMatrixType *const T, CompactMatrixType *const Tinv, double const SingularValue, bool const ReduceDimension)
{
  CompactMatrixType AT;
  CompactMatrixCopy(*A, &AT);
  if(CompactMatrixTranspose(&AT)) {printf("CompactMatrixTranspose failed in CompactMatrixSVDOrthogonalizeCol\n");return 1;};
  if(CompactMatrixGramSchmidtOrthogonalizeRowplusT(&AT,Anew,T,Tinv,SingularValue,ReduceDimension)) {printf("CompactMatrixSVDOrthogonalizeRowplusT failed in CompactMatrixSVDOrthogonalizeColplusT\n");return 1;};
  if(CompactMatrixTranspose(Anew)) {printf("CompactMatrixTranspose failed in CompactMatrixSVDOrthogonalizeColplusT\n");return 1;};
  FreeCompactMatrix(&AT);
  return 0;
}
uint32_t CompactMatrixGramSchmidtOrthogonalizeCol(CompactMatrixType * A, double const SingularValue, bool const ReduceDimension)
{
  if(CompactMatrixTranspose(A)) {printf("CompactMatrixTranspose failed in CompactMatrixGramSchmidtOrthogonalizeCol\n");return 1;};
  if(CompactMatrixGramSchmidtOrthogonalizeRow(A,SingularValue,ReduceDimension)) {printf("CompactMatrixGramSchmidtOrthogonalizeRow failed in CompactMatrixGramSchmidtOrthogonalizeCol\n");return 1;};
  if(CompactMatrixTranspose(A)) {printf("CompactMatrixTranspose failed in CompactMatrixGramSchmidtOrthogonalizeCol\n");return 1;};
  return 0;
}
uint32_t CompactMatrixGramSchmidtOrthogonalizeColplusT(CompactMatrixType const * A, CompactMatrixType *const Anew, CompactMatrixType *const T, CompactMatrixType *const Tinv, double const SingularValue, bool const ReduceDimension)
{
  CompactMatrixType AT;
  CompactMatrixCopy(*A, &AT);
  if(CompactMatrixTranspose(&AT)) {printf("CompactMatrixTranspose failed in CompactMatrixGramSchmidtOrthogonalizeCol\n");return 1;};
  if(CompactMatrixGramSchmidtOrthogonalizeRowplusT(&AT,Anew,T,Tinv,SingularValue,ReduceDimension)) {printf("CompactMatrixGramSchmidtOrthogonalizeRowplusT failed in CompactMatrixGramSchmidtOrthogonalizeColplusT\n");return 1;};
  if(CompactMatrixTranspose(Anew)) {printf("CompactMatrixTranspose failed in CompactMatrixGramSchmidtOrthogonalizeColplusT\n");return 1;};
  FreeCompactMatrix(&AT);
  return 0;
}

/*  Returns the Generalized Eigensystem of a non-orthnormal matrix.
 *  Given a basis psi_i
 *  Matrix M contains the values psi_i^* . H . psi_j
 *  Matrix S contains the values psi_i^* . psi_j
 *
 *  on exit val defines a diagonal matrix H and
 *  H = T^* . M . T^{\dag}
 *  I = T^* . S . T^{\dag}
 *  S = Sroot . Sroot^{\dag}
 *
 */
uint32_t CompactMatrixGeneralEigensystem(CompactMatrixType const *M, CompactMatrixType * const S, CompactMatrixType * const T, CompactMatrixType * const Sroot, double *val)
{
  int    info;
  uint32_t itype = 1;
  const char *jobz = "V";
  const char *uplo = "U";
  uint32_t   lwork = 1 + 6 * M->Ni + 2 * Sqr(M->Ni);
  uint32_t   liwork = 3 + 5 * M->Ni;
  int    iwork[liwork];
  double work[lwork];
  uint32_t i, j;
  if (M->Ni != M->Nj)
  {
    printf("Generalized Eigensystem called with non-square matrix\n");
    printf("ERORR!!!!!!\n");
    return 1;
  }
  if (CompactMatrixCopy(*S,Sroot))
  {
    printf("CompactMatrixCopy failed to copy S to Sroot in GeneralizedEigensystem\n");
    return 1;
  }
  strcpy(Sroot->Name, "Sqrt of ");
  strncat(Sroot->Name, S->Name, Max(0, MAXStrLength - 8));
  for (i = 0; i < Sroot->Ni; i++)
  {
    for (j = i + 1; j < Sroot->Nj; j++)
    {
      Sroot->VRe[i*Sroot->Nj+j] = 0;
    }
  }
  if (CompactMatrixCopy(*M,T))
  {
    printf("CompactMatrixCopy failed to copy M to T in GeneralizedEigensystem\n");
    FreeCompactMatrix(Sroot);
    return 1;
  }
  strcpy(T->Name, "Eigenvectors of ");
  strncat(T->Name, M->Name, Max(0, MAXStrLength - 16));
  if (M->QComplex == 0)
  {                             /* Real */
    dsygv(&itype, jobz, uplo, &M->Ni, T->VRe, &M->Ni, Sroot->VRe, &M->Ni, val, work, &lwork, &info);
    if (info)
    {
      printf("dsygv failed in Generalized Eigensystem with info value %i\n", info);
      printf("While working on matrix %s\n", M->Name);
      printf("continue diagonalization with dsygvd\n");
      fflush(stdout);
      FreeCompactMatrix(T);
      FreeCompactMatrix(Sroot);
      if (CompactMatrixCopy(*S,Sroot))
      {
        printf("CompactMatrixCopy failed to copy S to Sroot in GeneralizedEigensystem\n");
        return 1;
      }
      strcpy(Sroot->Name, "Sqrt of ");
      strncat(Sroot->Name, S->Name, Max(0, MAXStrLength - 8));
      for (i = 0; i < Sroot->Ni; i++)
      {
        for (j = i + 1; j < Sroot->Nj; j++)
        {
          Sroot->VRe[i*Sroot->Nj+j] = 0;
        }
      }
      if (CompactMatrixCopy(*M,T))
      {
        printf("MatrixCopy failed to copy M to T in GeneralizedEigensystem\n");
        FreeCompactMatrix(Sroot);
        return 1;
      }
      strcpy(T->Name, "Eigenvectors of ");
      strncat(T->Name, M->Name, Max(0, MAXStrLength - 16));
      dsygvd(&itype, jobz, uplo, &M->Ni, T->VRe, &M->Ni, Sroot->VRe, &M->Ni, val, work, &lwork, iwork, &liwork, &info);
      if (info)
      {
        printf("ERORR also dsygvd failed in Eigensystem with info value %i\n", info);
        fflush(stdout);
        return 1;
      }
    }
    return 0;
  }
  else
  {                             /* complex */
    uint32_t  lcwork = 2 * M->Ni + Sqr(M->Ni);
    Complex16 cwork[lcwork];
    if(Sroot->QComplex==0)
    {
      MakeCompactRealMatrixComplex(Sroot);
    }
    zhegv(&itype, jobz, uplo, &M->Ni, M->V, &M->Ni, Sroot->V, &M->Ni, val, cwork, &lcwork, work, &info);
    if (info != 0)
    {
      printf("zhegv failed in Generalized Eigensystem with info value %i\n", info);
      printf("While working on matrix %s\n", M->Name);
      printf("continue diagonalization with zhegvd\n");
      fflush(stdout);
      FreeCompactMatrix(T);
      FreeCompactMatrix(Sroot);
      if (CompactMatrixCopy(*S,Sroot))
      {
        printf("CompactMatrixCopy failed to copy S to Sroot in GeneralizedEigensystem\n");
        return 1;
      }
      strcpy(Sroot->Name, "Sqrt of ");
      strncat(Sroot->Name, S->Name, Max(0, MAXStrLength - 8));
      for (i = 0; i < Sroot->Ni; i++)
      {
        for (j = i + 1; j < Sroot->Nj; j++)
        {
          Sroot->VRe[i*Sroot->Nj+j] = 0;
        }
      }
      if (CompactMatrixCopy(*M,T))
      {
        printf("MatrixCopy failed to copy M to T in GeneralizedEigensystem\n");
        FreeCompactMatrix(Sroot);
        return 1;
      }
      strcpy(T->Name, "Eigenvectors of ");
      strncat(T->Name, M->Name, Max(0, MAXStrLength - 16));
      dsygvd(&itype, jobz, uplo, &M->Ni, T->VRe, &M->Ni, Sroot->VRe, &M->Ni, val, work, &lwork, iwork, &liwork, &info);
      if (info)
      {
        printf("ERORR also dsygvd failed in Eigensystem with info value %i\n", info);
        fflush(stdout);
        return 1;
      }
      if(Sroot->QComplex==0)
      {
        MakeCompactRealMatrixComplex(Sroot);
      }
      zhegvd(&itype, jobz, uplo, &M->Ni, M->V, &M->Ni, Sroot->V, &M->Ni, val, cwork, &lcwork, work, &lwork, iwork, &liwork, &info);
      if (info != 0)
      {
        printf("ERORR, also zhegvd failed in Generalized Eigensystem with info value %i\n", info);
        fflush(stdout);
        return 1;
      }
    }
  }
  return 0;
}














uint32_t CompactMatrixLeftRightNormalize(CompactMatrixType * const funL, CompactMatrixType * const funR)
{
  uint32_t i,j,k;
  Complex16 norm;
  double scale,phase;
  if((funL->Ni!=funR->Ni)||(funL->Nj!=funR->Nj))
  {
    printf("Dimensions of matrices in CompactMatrixLeftRightNormalize not ok to calculate funL^* . funR^T (%i %i)(%i %i)\n",funL->Ni,funL->Nj,funR->Ni,funR->Nj);
    return 1;
  }
  if(funL->QComplex==0) MakeCompactRealMatrixComplex(funL);
  
  for(i=0;i<funL->Ni;i++)
  {
    norm.Re=0;
    norm.Im=0;
    k=i*funL->Nj;
    for(j=0;j<funL->Nj;j++)
    {
      norm.Re+=funL->V[k].Re*funR->V[k].Re + funL->V[k].Im*funR->V[k].Im;
      norm.Im+=funL->V[k].Re*funR->V[k].Im - funL->V[k].Im*funR->V[k].Re;
      k++;
    }
    scale=sqrt(sqrt(Sqr(norm.Re)+Sqr(norm.Im)));
    phase=0.5*atan2(norm.Im,norm.Re);
    if(scale<100*DBL_EPSILON)
    {
      scale=0;
    }
    else
    {
      scale=1/scale;
    }
    for(k=i*funL->Nj;k<(i+1)*funL->Nj;k++)
    {
      norm.Re = funL->V[k].Re;
      norm.Im = funL->V[k].Im;
      funL->V[k].Re = scale * ( cos(phase)*norm.Re - sin(phase)*norm.Im);
      funL->V[k].Im = scale * ( cos(phase)*norm.Im + sin(phase)*norm.Re);
      norm.Re = funR->V[k].Re;
      norm.Im = funR->V[k].Im;
      funR->V[k].Re = scale * ( cos(phase)*norm.Re + sin(phase)*norm.Im);
      funR->V[k].Im = scale * ( cos(phase)*norm.Im - sin(phase)*norm.Re);
    }
  }
  return 0;
}

int SortNonHermitianEigensystemOrder(const void * elem1, const void * elem2)
{
  double f = *((double*)elem1);
  double s = *((double*)elem2);
  if (f > s) return  1;
  if (f < s) return -1;
  return 0;
}

/*  The Lapack standard does not order the eigen values of Non Hermitian matrices
 *  We here order them according to the real part of the eigen values from low to high
 */
uint32_t NonHermitianEigensystemOrder(Complex16 * val, CompactMatrixType * const funL, CompactMatrixType * const funR)
{
  // create a matrix with rows containing val[i] .. funL[i] .. funR[i]
  uint32_t i,j,k;
  double *data;
  if(!(funL->Ni==funR->Ni))
  {
    printf("Number of rows for funL and funR need to be equal in NonHermitianEigensystemOrder\n");
    return 1;
  }
  if(funL->QComplex==0)
  {
    if(funR->QComplex==0)
    { //  all functions real
      // fill data
      data = (double*)malloc((2+funL->Nj+funR->Nj)*funL->Ni*sizeof(double));
      if(data==NULL)
      {
        printf("malloc failed in NonHermitianEigensystemOrder\n");
        return 1;
      }
      k=0;
      for(i=0;i<funL->Ni;i++)
      {
        data[k]   = val[i].Re;
        k++;
        data[k] = val[i].Im;
        k++;
        for(j=0;j<funL->Nj;j++)
        {
          data[k] = funL->VRe[i*funL->Nj+j];
          k++;
        }
        for(j=0;j<funR->Nj;j++)
        {
          data[k] = funR->VRe[i*funR->Nj+j];
          k++;
        }
      }
      // sort data
      qsort(data,funL->Ni, (2+funL->Nj+funR->Nj)*sizeof(double), SortNonHermitianEigensystemOrder);
      // copy back to val, funL and funR
      k=0;
      for(i=0;i<funL->Ni;i++)
      {
        val[i].Re = data[k];
        k++;
        val[i].Im = data[k];
        k++;
        for(j=0;j<funL->Nj;j++)
        {
          funL->VRe[i*funL->Nj+j] = data[k];
          k++;
        }
        for(j=0;j<funR->Nj;j++)
        {
          funR->VRe[i*funR->Nj+j] = data[k];
          k++;
        }
      }
      //
      free(data);
    }
    else
    {  // L real R complex
      // fill data
      data = (double*)malloc((2+funL->Nj+2*funR->Nj)*funL->Ni*sizeof(double));
      if(data==NULL)
      {
        printf("malloc failed in NonHermitianEigensystemOrder\n");
        return 1;
      }
      k=0;
      for(i=0;i<funL->Ni;i++)
      {
        data[k]   = val[i].Re;
        k++;
        data[k] = val[i].Im;
        k++;
        for(j=0;j<funL->Nj;j++)
        {
          data[k] = funL->VRe[i*funL->Nj+j];
          k++;
        }
        for(j=0;j<funR->Nj;j++)
        {
          data[k] = funR->V[i*funR->Nj+j].Re;
          k++;
          data[k] = funR->V[i*funR->Nj+j].Im;
          k++;
        }
      }
      // sort data
      qsort(data,funL->Ni, (2+funL->Nj+2*funR->Nj)*sizeof(double), SortNonHermitianEigensystemOrder);
      // copy back to val, funL and funR
      k=0;
      for(i=0;i<funL->Ni;i++)
      {
        val[i].Re = data[k];
        k++;
        val[i].Im = data[k];
        k++;
        for(j=0;j<funL->Nj;j++)
        {
          funL->VRe[i*funL->Nj+j] = data[k];
          k++;
        }
        for(j=0;j<funR->Nj;j++)
        {
          funR->V[i*funR->Nj+j].Re = data[k];
          k++;
          funR->V[i*funR->Nj+j].Im = data[k];
          k++;
        }
      }
      //
      free(data);
    }
  }
  else
  {
    if(funR->QComplex==0)
    { // L complex R real
      // fill data
      data = (double*)malloc((2+2*funL->Nj+funR->Nj)*funL->Ni*sizeof(double));
      if(data==NULL)
      {
        printf("malloc failed in NonHermitianEigensystemOrder\n");
        return 1;
      }
      k=0;
      for(i=0;i<funL->Ni;i++)
      {
        data[k]   = val[i].Re;
        k++;
        data[k] = val[i].Im;
        k++;
        for(j=0;j<funL->Nj;j++)
        {
          data[k] = funL->V[i*funL->Nj+j].Re;
          k++;
          data[k] = funL->V[i*funL->Nj+j].Im;
          k++;
        }
        for(j=0;j<funR->Nj;j++)
        {
          data[k] = funR->VRe[i*funR->Nj+j];
          k++;
        }
      }
      // sort data
      qsort(data,funL->Ni, (2+2*funL->Nj+funR->Nj)*sizeof(double), SortNonHermitianEigensystemOrder);
      // copy back to val, funL and funR
      k=0;
      for(i=0;i<funL->Ni;i++)
      {
        val[i].Re = data[k];
        k++;
        val[i].Im = data[k];
        k++;
        for(j=0;j<funL->Nj;j++)
        {
          funL->V[i*funL->Nj+j].Re = data[k];
          k++;
          funL->V[i*funL->Nj+j].Im = data[k];
          k++;
        }
        for(j=0;j<funR->Nj;j++)
        {
          funR->VRe[i*funR->Nj+j] = data[k];
          k++;
        }
      }
      //
      free(data);
    }
    else
    { // all functions complex
      // fill data
      data = (double*)malloc((2+2*funL->Nj+2*funR->Nj)*funL->Ni*sizeof(double));
      if(data==NULL)
      {
        printf("malloc failed in NonHermitianEigensystemOrder\n");
        return 1;
      }
      k=0;
      for(i=0;i<funL->Ni;i++)
      {
        data[k]   = val[i].Re;
        k++;
        data[k] = val[i].Im;
        k++;
        for(j=0;j<funL->Nj;j++)
        {
          data[k] = funL->V[i*funL->Nj+j].Re;
          k++;
          data[k] = funL->V[i*funL->Nj+j].Im;
          k++;
        }
        for(j=0;j<funR->Nj;j++)
        {
          data[k] = funR->V[i*funR->Nj+j].Re;
          k++;
          data[k] = funR->V[i*funR->Nj+j].Im;
          k++;
        }
      }
      // sort data
      qsort(data,funL->Ni, (2+2*funL->Nj+2*funR->Nj)*sizeof(double), SortNonHermitianEigensystemOrder);
      // copy back to val, funL and funR
      k=0;
      for(i=0;i<funL->Ni;i++)
      {
        val[i].Re = data[k];
        k++;
        val[i].Im = data[k];
        k++;
        for(j=0;j<funL->Nj;j++)
        {
          funL->V[i*funL->Nj+j].Re = data[k];
          k++;
          funL->V[i*funL->Nj+j].Im = data[k];
          k++;
        }
        for(j=0;j<funR->Nj;j++)
        {
          funR->V[i*funR->Nj+j].Re = data[k];
          k++;
          funR->V[i*funR->Nj+j].Im = data[k];
          k++;
        }
      }
      //
      free(data);
    }
  }
  return 0;
}


/*
 *  D=funL^*.A.funR^{T} with D a diagonal matrix with the values val on the diagonal and
 *  funL^* . funR^T = I
 *  input A defined, fun not initiallized, val allocated
 */
uint32_t NonHermitianCompactMatrixEigensystem(CompactMatrixType const *A, CompactMatrixType * const funL, CompactMatrixType * const funR, Complex16 *val)
{
  int    info;
  const char *jobvl = "V";
  const char *jobvr = "V";
  uint32_t    lwork = 1 + 6 * A->Ni + 2 * Sqr(A->Ni);
  uint32_t i,k;
  
  if (A->Ni != A->Nj)
  {
    printf("NonHermitianCompactMatrixEigensystem called with non-square matrix\n");
    printf("ERORR!!!!!!\n");
    return 1;
  }
  
  funL->Ni=A->Ni;
  funL->Nj=A->Nj;
  funL->QComplex=1;
  if(InitCompactMatrix(funL))
  {
    printf("InitCompactMatrix failed in NonHermitianCompactMatrixEigensystem\n");
    return 1;
  }
  funL->a0=0;
  funL->b0=1;
  funL->mu=0;
  strcpy(funL->Name, "Left Eigen vectors");
  
  funR->Ni=A->Ni;
  funR->Nj=A->Nj;
  funR->QComplex=1;
  if(InitCompactMatrix(funR))
  {
    printf("InitCompactMatrix failed in NonHermitianCompactMatrixEigensystem\n");
    FreeCompactMatrix(funL);
    return 1;
  }
  funR->a0=0;
  funR->b0=1;
  funR->mu=0;
  strcpy(funL->Name, "Right Eigen vectors");
  
  if (A->QComplex == 0)
  {                             /* Real */
    double *work,*valRe,*valIm;
    double *Vl, *Vr, *Acopy;
    work = (double*)malloc(lwork*sizeof(double));
    if(work==NULL)
    {
      printf("malloc for work failed in CompactMatrixEigensystem\n");
      FreeCompactMatrix(funL);
      FreeCompactMatrix(funR);
      return 1;
    }
    valRe= (double*)malloc(A->Ni*sizeof(double));
    if(valRe==NULL)
    {
      printf("malloc failed in NonHermitianCompactMatrixEigensystem\n");
      free(work);
      FreeCompactMatrix(funL);
      FreeCompactMatrix(funR);
      return 1;
    }
    valIm= (double*)malloc(A->Ni*sizeof(double));
    if(valIm==NULL)
    {
      printf("malloc failed in NonHermitianCompactMatrixEigensystem\n");
      free(work);
      free(valRe);
      FreeCompactMatrix(funL);
      FreeCompactMatrix(funR);
      return 1;
    }
    Acopy= (double*)malloc(A->Ni*A->Ni*sizeof(double));
    if(Acopy==NULL)
    {
      printf("malloc failed in NonHermitianCompactMatrixEigensystem\n");
      free(work);
      free(valRe);
      free(valIm);
      FreeCompactMatrix(funL);
      FreeCompactMatrix(funR);
      return 1;
    }
    Vl= (double*)malloc(A->Ni*A->Ni*sizeof(double));
    if(Vl==NULL)
    {
      printf("malloc failed in NonHermitianCompactMatrixEigensystem\n");
      free(work);
      free(valRe);
      free(valIm);
      free(Acopy);
      FreeCompactMatrix(funL);
      FreeCompactMatrix(funR);
      return 1;
    }
    Vr= (double*)malloc(A->Ni*A->Ni*sizeof(double));
    if(Vr==NULL)
    {
      printf("malloc failed in NonHermitianCompactMatrixEigensystem\n");
      free(work);
      free(valRe);
      free(valIm);
      free(Acopy);
      free(Vl);
      FreeCompactMatrix(funL);
      FreeCompactMatrix(funR);
      return 1;
    }
    memcpy(Acopy, A->VRe, A->Ni*A->Ni*sizeof(double));
    dgeev(jobvl, jobvr, &(A->Ni), Acopy, &(A->Ni), valRe, valIm, Vl, &(A->Ni), Vr, &(A->Ni), work, &lwork, &info);
    free(Acopy);
    free(work);
    if (info)
    {
      printf("dgeev failed in NonHermitianCompactMatrixEigensystem with info value %i\n", info);
      printf("While working on matrix %s\n", A->Name);
      fflush(stdout);
      free(valRe);
      free(valIm);
      free(Vl);
      free(Vr);
      return 1;
    }
    for(i=0;i<A->Ni;i++)
    {
      val[i].Re=valRe[i];
      val[i].Im=valIm[i];
      if(val[i].Im==0)
      { // R and L switch due to Fortran to C row to column major
        for(k=i*A->Nj;k<(i+1)*A->Nj;k++)
        {
          funR->V[k].Re=Vl[k];
          funR->V[k].Im=0;
          funL->V[k].Re=Vr[k];
          funL->V[k].Im=0;
        }
      }
      else
      {
        for(k=i*A->Nj;k<(i+1)*A->Nj;k++)
        { // R and L switch and become complex conjugated due to Fortran to C row to column major
          funR->V[k].Re      = Vl[k];
          funR->V[k].Im      =-Vl[k+A->Nj];
          funL->V[k].Re      = Vr[k];
          funL->V[k].Im      =-Vr[k+A->Nj];
          funR->V[k+A->Nj].Re= Vl[k];
          funR->V[k+A->Nj].Im= Vl[k+A->Nj];
          funL->V[k+A->Nj].Re= Vr[k];
          funL->V[k+A->Nj].Im= Vr[k+A->Nj];
        }
        i++;
        val[i].Re=valRe[i];
        val[i].Im=valIm[i];
      }
    }
    free(valRe);
    free(valIm);
    free(Vl);
    free(Vr);
  }
  else
  {                             /* Complex */
    Complex16 *Acopy,*work;
    double *rwork;
    Acopy = (Complex16*)malloc(A->Ni*A->Nj*sizeof(Complex16));
    if(Acopy==NULL)
    {
      printf("malloc failed in NonHermitianCompactMatrixEigensystem\n");
      return 1;
    }
    work  = (Complex16*)malloc(lwork*sizeof(Complex16));
    if(work==NULL)
    {
      free(Acopy);
      printf("malloc failed in NonHermitianCompactMatrixEigensystem\n");
      return 1;
    }
    rwork = (double*)malloc(2*A->Ni*sizeof(double));
    if(rwork==NULL)
    {
      free(Acopy);
      free(work);
      printf("malloc failed in NonHermitianCompactMatrixEigensystem\n");
      return 1;
    }
    memcpy(Acopy, A->VRe, A->Ni*A->Ni*sizeof(Complex16));
    zgeev(jobvl,jobvr,&(A->Ni),Acopy,&(A->Ni),val,funR->V,&(A->Ni),funL->V,&(A->Ni),work,&lwork,rwork,&info);
    CompactMatrixConjugate(funL);
    CompactMatrixConjugate(funR);
    if (info)
    {
      printf("zgeev failed in NonHermitianCompactMatrixEigensystem with info value %i\n", info);
      printf("While working on matrix %s\n", A->Name);
      fflush(stdout);
      return 1;
    }
    free(Acopy);
    free(work);
    free(rwork);
  }
  // normalize vector funR and funL such that funL^* . fulR^T = 1
  CompactMatrixLeftRightNormalize(funL,funR);
  NonHermitianEigensystemOrder(val, funL, funR);
  return 0;
}

/*
 *  D=fun^*.A.fun^{T} with D a diagonal matrix with the values val on the diagonal
 *  input A defined, fun not initiallized, val allocated
 */
uint32_t CompactMatrixEigensystem(CompactMatrixType const *A, CompactMatrixType * const fun, double *val)
{
  if(A->Ni==0)
  {
    CompactMatrixCopy(*A, fun);
    return 0;
  }
  int    info;
  const char *jobz = "V";
  const char *uplo = "U";
  uint32_t    lwork = 1 + 6 * A->Ni + 2 * Sqr(A->Ni);
  uint32_t    liwork = 3 + 5 * A->Ni;
  int    *iwork;
  double *work;
  
  iwork = (int*)malloc(liwork*sizeof(int));
  if(iwork==NULL)
  {
    printf("malloc for iwork failed in CompactMatrixEigensystem\n");
    return 1;
  }
  work = (double*)malloc(lwork*sizeof(double));
  if(work==NULL)
  {
    printf("malloc for work failed in CompactMatrixEigensystem\n");
    return 1;
  }
  if (A->Ni != A->Nj)
  {
    printf("CompactMatrixEigensystem called with non-square matrix\n");
    printf("ERORR!!!!!!\n");
    return 1;
  }
  if (A->QComplex == 0)
  {                             /* Real */
    if (CompactMatrixCopy(*A, fun))
    {
      printf("Error in CompactMatrixCopy in CompactMatrixEigensystem\n");
      return 1;
    }
    strcpy(fun->Name, "Eigenvectors of ");
    strncat(fun->Name, A->Name, Max(0, MAXStrLength - 16));
    dsyev(jobz, uplo, &A->Ni, fun->VRe, &A->Ni, val, work, &lwork, &info);
    if (info)
    {
      printf("dsyev failed in CompactMatrixEigensystem with info value %i\n", info);
      printf("While working on matrix %s\n", A->Name);
      printf("continue diagonalization with dsyvd\n");
      fflush(stdout);
      FreeCompactMatrix(fun);
      if (CompactMatrixCopy(*A,fun))
        return 1;
      strcpy(fun->Name, "Eigenvectors of ");
      strncat(fun->Name, A->Name, Max(0, MAXStrLength - 16));
      dsyevd(jobz, uplo, &A->Ni, fun->VRe, &A->Ni, val, work, &lwork, iwork, &liwork, &info);
      if (info)
      {
        printf("also dsyevd failed in CompactMatrixEigensystem with info value %i\n", info);
        printf("ERORR!!!!!!!!!!\n");
        fflush(stdout);
        return 1;
      }
      free(iwork);
      free(work);
      return 0;
    }
    free(iwork);
    free(work);
    return 0;
  }
  else
  {                             /* complex */
    uint32_t    lcwork = 2 * A->Ni + Sqr(A->Ni);
    Complex16 *cwork;
    cwork = (Complex16*)malloc(lcwork*sizeof(Complex16));
    if(cwork==NULL)
    {
      printf("malloc for cwork failed in CompactMatrixEigensystem\n");
      return 1;
    }
    if (CompactMatrixCopy(*A, fun)) return 1;
    strcpy(fun->Name, "Eigenvectors of ");
    strncat(fun->Name, A->Name, Max(0, MAXStrLength - 16));
    zheev(jobz, uplo, &A->Ni, fun->V, &A->Ni, val, cwork, &lcwork, work, &info);
    if (info != 0)
    {
      printf("zheev failed in CompactMatrixEigensystem with info value %i\n", info);
      printf("While working on matrix %s\n", A->Name);
      printf("continue diagonalization with zheevd\n");
      fflush(stdout);
      FreeCompactMatrix(fun);
      if (CompactMatrixCopy(*A,fun))
        return 1;
      strcpy(fun->Name, "Eigenvectors of ");
      strncat(fun->Name, A->Name, Max(0, MAXStrLength - 16));
      zheevd(jobz, uplo, &A->Ni, fun->V, &A->Ni, val, cwork, &lcwork, work, &lwork, iwork, &liwork, &info);
      if (info != 0)
      {
        printf("also zheevd failed in CompactMatrixEigensystem with info value %i\n", info);
        printf("ERORR!!!!!!!!!!\n");
        fflush(stdout);
        free(iwork);
        free(cwork);
        free(work);
        return 1;
      }
    }
    CompactMatrixConjugate(fun);
    free(iwork);
    free(cwork);
    free(work);
    return 0;
  }
  return 0;
}


/*  Generates a tridiagonal form of matrix H.
 *  We assume Ni==Nj
 *  H needs to be Hermitian, only the upper rigth part of H is read, Hermiticy is not tested
 *  If on entry T==NULL only the resulting tri-diagonal matrix is returned
 *  If T != NULL we assume it points to a compact matrix structure and this structure is filled with the matrix T such that
 *  HTri = T^* H T^T
 */
int CompactMatrixTriDiagonalize(CompactMatrixType *H, CompactMatrixType *T, TriDiagonalMatrixType *HTri)
{
  if(H->Ni==0)
  {
    if(T!=NULL)
    {
      T->Ni=0;
      T->Nj=0;
      T->QComplex=0;
      InitCompactMatrix(T);
    }
    HTri->N=0;
    InitTriDiagonalMatrix(HTri);
    HTri->a[0]=H->a0;
    HTri->mu=H->mu;
    return 0;
  }
  CompactMatrixType *HCopy;
  double *work=NULL, *tau=NULL;
  Complex16 *Cwork=NULL, *Ctau=NULL;
  int info;
  uint32_t lwork= 1 + 6 * H->Ni + 2 * Sqr(H->Ni);
  uint32_t lda=H->Nj;
  char UPLO='L';
  if(T==NULL)
  {
    HCopy = (CompactMatrixType*)malloc(sizeof(CompactMatrixType));
  }
  else
  {
    HCopy = T;
  }
  CompactMatrixCopy(*H, HCopy);
  HTri->N=H->Nj;
  InitTriDiagonalMatrix(HTri);
  if(H->QComplex)
  {
    CompactMatrixConjugate(HCopy);
    Cwork = (Complex16*)malloc(lwork*sizeof(Complex16));
    Ctau = (Complex16*)malloc(H->Nj*sizeof(Complex16));
    zhetrd(&UPLO, &lda, HCopy->V, &lda, &(HTri->a[1]) ,&(HTri->b[1]), Ctau, Cwork, &lwork, &info);
  }
  else
  {
    work = (double*)malloc(lwork*sizeof(double));
    tau  = (double*)malloc(H->Nj*sizeof(double));
    dsytrd(&UPLO, &lda, HCopy->VRe, &lda, &(HTri->a[1]), &(HTri->b[1]), tau, work, &lwork, &info);
  }
  if(T==NULL)
  { // no need to calculate the matrix T. Cleanup and finish.
    if(H->QComplex)
    {
      free(Cwork);
      free(Ctau);
    }
    else
    {
      free(work);
      free(tau);
    }
    FreeCompactMatrix(HCopy);
    return 0;
  }
  // next step, calculte T == HCopy
  if(H->QComplex)
  {
    zungtr(&UPLO, &lda, T->V, &lda, Ctau, Cwork, &lwork, &info);
  }
  else
  {
    dorgtr(&UPLO, &lda, T->VRe, &lda, tau, work, &lwork, &info);
  }
  //Cleanup and finish.
  if(H->QComplex)
  {
    free(Cwork);
    free(Ctau);
  }
  else
  {
    free(work);
    free(tau);
  }
  return 0;
}





/*  On entrance
 *  A is a matrix
 *  SingularValue is the smallest positive eigenvalue of the matrix A^dag A that is considered non-singular
 *
 *  returned:
 *  A = q . r
 *  with q orthnoromal (q^C.q=1) and r hermitian (r=r^C)
 */
uint32_t CompactMatrixQRDecompositionEigenHermitian(CompactMatrixType const *A, double const SingularValue, CompactMatrixType *q, CompactMatrixType *r)
{
  uint32_t i,NSingular=0;
  CompactMatrixType Asqr,fun;
  double val[A->Ni];
  CompactMatrixMultiply('C', 'N', A, A, &Asqr);
  CompactMatrixEigensystem(&Asqr, &fun, val);
  for(i=0;i<A->Ni;i++)
  {
    if(val[i]<SingularValue)
    {
      NSingular++;
      printf("Singular value found (i=%i val=%22.15E) in CompactMatrixQRDecompositionEigenHermitian\n",i,val[i]);
      val[i]=0;
    }
  }
  MakeCompactComplexMatrixReal(&Asqr);
  SetCompactMatrixToZero(&Asqr);
  for(i=0;i<A->Ni;i++)
  {
    Asqr.VRe[i+i*Asqr.Ni]=sqrt(val[i]);
  }
  CompactMatrixMultiply3('T', 'N', '*', &fun, &Asqr, &fun, r);
  for(i=0;i<A->Ni;i++)
  {
    if(val[i]<SingularValue)
    {
      Asqr.VRe[i+i*Asqr.Ni]=0;
    }
    else
    {
      Asqr.VRe[i+i*Asqr.Ni]=(double)1.0/(sqrt(Abs(val[i]))+100*DBL_MIN);
    }
  }
  CompactMatrixMultiply4('T','N','*','C',&fun, &Asqr, &fun, A,q);
  FreeCompactMatrix(&Asqr);
  FreeCompactMatrix(&fun);
  return NSingular;
}

uint32_t CompactMatrixSqrPreAllocated(CompactMatrixType const *A, CompactMatrixType * const SqrA)
{
  return CompactMatrixMultiplyPreAllocated('N', 'N', A, A, SqrA);
}
uint32_t CompactMatrixSqr(CompactMatrixType const *A, CompactMatrixType * const SqrA)
{
  return CompactMatrixMultiply('N', 'N', A, A, SqrA);
}


/*  Given a square matrix calculates PowA = A^n
 *  It is assumed that A->Ni == A->Nj
 *  PowA is allocated but not initiallized
 */
uint32_t CompactMatrixPow(CompactMatrixType const *A, CompactMatrixType * const PowA, const int n)
{
  CompactMatrixType temp;
  if(n==0)
  {
    PowA->Ni = A->Ni;
    PowA->Nj = A->Nj;
    PowA->QComplex=0;
    if(InitCompactIdentityMatrix(PowA))
    {
      printf("CompactMatrixPow faild at InitCompactIdentityMatrix\n");
      return 1;
    }
    return 0;
  }
  if(n==1)
  {
    if(CompactMatrixCopy(*A, PowA))
    {
      printf("CompactMatrixPow faild at CompactMatrixCopy\n");
      return 1;
    }
    return 0;
  }
  if(n<0)
  {
    if(CompactMatrixCopy(*A, &temp))
    {
      printf("CompactMatrixPow faild at CompactMatrixCopy\n");
      return 1;
    }
    if(CompactMatrixInverse(&temp))
    {
      printf("CompactMatrixPow faild at CompactMatrixInverse\n");
      return 1;
    }
    if(CompactMatrixPow(&temp,PowA,-n))
    {
      printf("CompactMatrixPow faild at CompactMatrixPow(%d)\n",-n);
      return 1;
    }
    FreeCompactMatrix(&temp);
    return 0;
  }
  if((n%2)==0) //even
  {
    if(CompactMatrixSqr(A,&temp))
    {
      printf("CompactMatrixPow faild at CompactMatrixSqr\n");
      return 1;
    }
    if(CompactMatrixPow(&temp,PowA,n/2))
    {
      printf("CompactMatrixPow faild at CompactMatrixPow(%d)\n",n/2);
      return 1;
    }
    FreeCompactMatrix(&temp);
    return 0;
  }
  else //odd
  {
    if(CompactMatrixPow(A,&temp,n-1))
    {
      printf("CompactMatrixPow faild at CompactMatrixPow(%d)\n",n-1);
      return 1;
    }
    if(CompactMatrixMultiply('N', 'N', A, &temp, PowA))
    {
      printf("CompactMatrixPow faild at CompactMatrixMultiply\n");
      return 1;
    }
    FreeCompactMatrix(&temp);
    return 0;
  }
}

/*  Given an  hermitian matrix, HermitianCompactMatrixExp calculates its exponent
 *   ExpA = e^A
 *
 */
uint32_t HermitianCompactMatrixExp(CompactMatrixType const *A, CompactMatrixType * const ExpA)
{
  int i;
  double val[A->Ni];
  CompactMatrixType fun;
  if(CompactMatrixEigensystem(A, &fun, val))
  {
    printf("CompactMatrixEigensystem failed in HermitianCompactMatrixExp\n");
    return 1;
  }
  for(i=0;i<fun.Ni;i++)
  {
    val[i]=exp(val[i]);
  }
  TransformDiagonalMatrixByCompactMatrixTtranATcon(val, &fun, ExpA);
  FreeCompactMatrix(&fun);
  return 0;
}

uint32_t NonHermitianCompactMatrixExp(CompactMatrixType const *A, CompactMatrixType * const ExpA)
{
  int i;
  double temp;
  Complex16 val[A->Ni];
  CompactMatrixType funL,funR;
  if(NonHermitianCompactMatrixEigensystem(A, &funL, &funR, val))
  {
    printf("NonHermitianCompactMatrixEigensystem failed in NonHermitianCompactMatrixExp\n");
    return 1;
  }
  for(i=0;i<A->Ni;i++)
  {
    temp = exp(val[i].Re);
    val[i].Re = temp * cos(val[i].Im);
    val[i].Im = temp * sin(val[i].Im);
  }
  TransformComplexDiagonalMatrixByCompactMatrixTRtranATLcon(val, &funL, &funR, ExpA);
  FreeCompactMatrix(&funL);
  FreeCompactMatrix(&funR);
  return 0;
}

/*  Given an semi definite hermitian matrix, CompactMatrixSqrt calculates its sqrt, such that
 *   B * B^{dag} = A
 *
 */
uint32_t CompactMatrixSqrt(CompactMatrixType const *A, CompactMatrixType * const SqrtA)
{
  int i,j,retval=0;
  double val[A->Ni],sumval;
  CompactMatrixType fun;
  CompactMatrixEigensystem(A, &fun, val);
  for(i=0;i<fun.Ni;i++)
  {
    if(val[i]<=0)
    {
      sumval=0;
      for(j=0;j<A->Ni;j++)
      {
        sumval+=Abs(val[j]);
      }
      if(val[i]<-1000*DBL_EPSILON*sumval)
      {
        printf("Negative value in CompactMatrixSqrt, eigen value is %22.15E at pos %i (trace |val| = %22.15E)\n",val[i],i,sumval);
        retval++;
      }
      val[i]=0;
      if(A->QComplex==0)
      {
        for(j=0;j<A->Nj;j++)
        {
          fun.VRe[j+i*A->Nj]=0;
        }
      }
      else
      {
        for(j=0;j<A->Nj;j++)
        {
          fun.V[j+i*A->Nj].Re=0;
          fun.V[j+i*A->Nj].Im=0;
        }
      }
    }
    val[i]=sqrt(val[i]);
  }
  TransformDiagonalMatrixByCompactMatrixTtranATcon(val, &fun, SqrtA);
  FreeCompactMatrix(&fun);
  return retval;
}

/*  A^{-1/2} for a hermitian positive semi definite A.
 *
 *  A direction with an eigenvalue of zero has no inverse. It is left at zero, the same
 *  treatment CompactMatrixSqrt above gives it, so what comes back is the inverse on the
 *  space A actually spans and zero on its null space. A physical caller reaches that case
 *  when A is a spectral weight and some direction carries none of it.
 *
 *  The threshold is relative to the trace of |A|, as a zero eigenvalue of a matrix built
 *  from a product does not come out as exactly zero.
 */
uint32_t CompactMatrixInverseSqrt(CompactMatrixType const *A, CompactMatrixType * const InvSqrtA)
{
  int i,j,retval=0;
  double val[A->Ni],sumval=0;
  CompactMatrixType fun;
  CompactMatrixEigensystem(A, &fun, val);
  for(j=0;j<A->Ni;j++)
  {
    sumval+=Abs(val[j]);
  }
  for(i=0;i<fun.Ni;i++)
  {
    if(val[i]<=1000*DBL_EPSILON*sumval)
    {
      if(val[i]<-1000*DBL_EPSILON*sumval)
      {
        printf("Negative value in CompactMatrixInverseSqrt, eigen value is %22.15E at pos %i (trace |val| = %22.15E)\n",val[i],i,sumval);
        retval++;
      }
      val[i]=0;
      if(A->QComplex==0)
      {
        for(j=0;j<A->Nj;j++)
        {
          fun.VRe[j+i*A->Nj]=0;
        }
      }
      else
      {
        for(j=0;j<A->Nj;j++)
        {
          fun.V[j+i*A->Nj].Re=0;
          fun.V[j+i*A->Nj].Im=0;
        }
      }
    }
    else
    {
      val[i]=1.0/sqrt(val[i]);
    }
  }
  TransformDiagonalMatrixByCompactMatrixTtranATcon(val, &fun, InvSqrtA);
  FreeCompactMatrix(&fun);
  return retval;
}

uint32_t NonHermitianCompactMatrixSqrt(CompactMatrixType const *A, CompactMatrixType * const SqrtA)
{
  //TODO: all
  return 1;
}

/* Given a compact matrix A calculates (w-A)^{-1}_[[1;;BlockSize,1;;BlockSize]] on an energy mesh given by Emin and Emax with N points
 *
 */
uint32_t CompactMatrixToSpectra(CompactMatrixType const *A, SpectraType * const G, uint32_t const BlockSize, double Emin, double Emax, uint32_t N)
{
  //TODO: test
  
  SpectraType *Gtemp;
  CompactMatrixType fun;
  int i,j,k,iEmin,iEplus;
  double *val,wmin,wplus;
  if(A->Ni!=A->Nj)
  {
    printf("CompactMatrixToSpectra expects square matrix as imput Ni=%i Nj=%i\n",A->Ni,A->Nj);
    return 1;
  }
  if(BlockSize>A->Ni)
  {
    printf("CompactMatrixToSpectra want more spectra (%i) then elements in matrix Ni=Nj=%i\n",BlockSize,A->Nj);
    return 1;
  }
  val= (double*)malloc(A->Ni*sizeof(double));
  if(val==NULL)
  {
    printf("malloc failed in CompactMatrixToSpectra\n");
    return 1;
  }
  G[0].N=N;
  G[0].Emin=Emin;
  G[0].Emax=Emax;
  G[0].EminPole=Emin;
  G[0].EmaxPole=Emax;
  G[0].dE=(Emax-Emin)/((double)N);
  G[0].Gamma=0;
  for(i=0;i<BlockSize*BlockSize;i++)
  {
    SpectrumAllocate(G[0], &(G[i]));
  }
  if(CompactMatrixEigensystem(A, &fun, val))
  {
    printf("CompactMatrixEigensystem failed in CompactMatrixToSpectra\n");
    return 1;
  }
  if(A->QComplex)
  {
    Gtemp= (SpectraType*)malloc(BlockSize*BlockSize*sizeof(SpectraType));
    for(i=0;i<BlockSize*BlockSize;i++)
    {
      SpectrumAllocate(G[0], &(Gtemp[i]));
    }
    for(k=0;k<fun.Ni;k++)
    {
      if(val[k]<=Emin)
      {
        iEmin=0;
        iEplus=0;
        wmin=1;
        wplus=0;
      }
      else if(val[k]>=Emax)
      {
        iEmin=N;
        iEplus=N;
        wmin=1;
        wplus=0;
      }
      else
      {
        iEmin=(uint32_t)((val[k]-Emin)/G[0].dE);
        iEplus=iEmin+1;
        wplus=((val[k]-Emin)/G[0].dE)-iEmin;
        wmin=1-wplus;
      }
      /* add the im part of the lorenzian to energy iEmin and iEplus */
      for(i=0;i<BlockSize;i++)
      {
        for(j=0;j<BlockSize;j++)
        {
          Gtemp[i*BlockSize+j].Im[iEmin] -= wmin *(fun.V[k*fun.Nj+j].Re*fun.V[k*fun.Nj+i].Re + fun.V[k*fun.Nj+j].Im*fun.V[k*fun.Nj+i].Im);
          Gtemp[i*BlockSize+j].Im[iEplus]-= wplus*(fun.V[k*fun.Nj+j].Re*fun.V[k*fun.Nj+i].Re + fun.V[k*fun.Nj+j].Im*fun.V[k*fun.Nj+i].Im);
          Gtemp[i*BlockSize+j].Re[iEmin] += wmin *(fun.V[k*fun.Nj+j].Re*fun.V[k*fun.Nj+i].Im - fun.V[k*fun.Nj+j].Im*fun.V[k*fun.Nj+i].Re);
          Gtemp[i*BlockSize+j].Re[iEplus]+= wplus*(fun.V[k*fun.Nj+j].Re*fun.V[k*fun.Nj+i].Im - fun.V[k*fun.Nj+j].Im*fun.V[k*fun.Nj+i].Re);
        }
      }
    }
  }
  else
  {
    Gtemp=G;
    for(k=0;k<fun.Ni;k++)
    {
      if(val[k]<=Emin)
      {
        iEmin=0;
        iEplus=0;
        wmin=1;
        wplus=0;
      }
      else if(val[k]>=Emax)
      {
        iEmin=N;
        iEplus=N;
        wmin=1;
        wplus=0;
      }
      else
      {
        iEmin=(uint32_t)((val[k]-Emin)/G[0].dE);
        iEplus=iEmin+1;
        wplus=((val[k]-Emin)/G[0].dE)-iEmin;
        wmin=1-wplus;
      }
      /* add the im part of the lorenzian to energy iEmin and iEplus */
      for(i=0;i<BlockSize;i++)
      {
        for(j=0;j<BlockSize;j++)
        {
          G[i*BlockSize+j].Im[iEmin] -= wmin *fun.VRe[k*fun.Nj+j]*fun.VRe[k*fun.Nj+i];
          G[i*BlockSize+j].Im[iEplus]-= wplus*fun.VRe[k*fun.Nj+j]*fun.VRe[k*fun.Nj+i];
        }
      }
    }
  }
  /* the im part of the lorenzians is there, we now need to add the real part */
  /* note that in the complex case some lorenzians have a complex pre-factor */
  /* weight w_i at energy E_i contributes w_i / (w - E_i) to all energies w, except at w=E_i */
  for(i=0;i<BlockSize*BlockSize;i++)
  {
    for(j=0;j<=N;j++)
    {
      for(k=0;k<j;k++)
      {
        G[i].Re[k]-=Gtemp[i].Im[j]/(Pi*(k-j));
      }
      for(k=j+1;k<=N;k++)
      {
        G[i].Re[k]-=Gtemp[i].Im[j]/(Pi*(k-j));
      }
    }
  }
  if(A->QComplex)
  {
    for(i=0;i<BlockSize*BlockSize;i++)
    {
      for(j=0;j<=N;j++)
      {
        for(k=0;k<j;k++)
        {
          G[i].Im[k]+=Gtemp[i].Re[j]/(Pi*(k-j));
        }
        for(k=j+1;k<=N;k++)
        {
          G[i].Im[k]+=Gtemp[i].Re[j]/(Pi*(k-j));
        }
      }
      for(j=0;j<=N;j++)
      {
        G[i].Re[j]+=Gtemp[i].Re[j];
        G[i].Im[j]+=Gtemp[i].Im[j];
      }
      SpectrumFree(&(Gtemp[i]));
    }
    free(Gtemp);
  }
  for(i=0;i<BlockSize*BlockSize;i++)
  {
    for(j=0;j<=N;j++)
    {
      G[i].Im[j]*=Pi/G[i].dE;
      G[i].Re[j]*=Pi/G[i].dE;
    }
  }
  return 0;
}

/*  Given matrix A with eigenvalues Aval
 *  split A=Apos+Aneg such that
 *  the eigenvalues of Apos are positive or zero and the eigenvalues of Aneg are positive or negative
 */
uint32_t HermitianCompactMatrixToPosAndNeg(CompactMatrixType const * const A, CompactMatrixType * const Apos, CompactMatrixType * const Aneg, bool *PosDefinite, bool *NegDefinite, double *NegWeight)
{
  int i,ipos,ineg,izero;
  double *val,*valpos,*valneg;
  CompactMatrixType fun;
  if(A->Ni!=A->Nj)
  {
    printf("Error in HermitianCompactMatrixToPosAndNeg. A not a square matrix\n");
    return 1;
  }
  val = (double*)malloc(A->Ni * sizeof(double));
  valpos = (double*)malloc(A->Ni * sizeof(double));
  valneg = (double*)malloc(A->Ni * sizeof(double));
  CompactMatrixEigensystem(A, &fun, val);
  ipos=0;
  ineg=0;
  izero=0;
  *NegWeight=0;
  for(i=0;i<A->Ni;i++)
  {
    if(val[i]<0)
    {
      valpos[i]=0;
      valneg[i]=val[i];
      *NegWeight+=val[i];
      ineg++;
    }
    else if(val[i]>0)
    {
      valpos[i]=val[i];
      valneg[i]=0;
      ipos++;
    }
    else
    {
      valpos[i]=0;
      valneg[i]=0;
      izero++;
    }
  }
  *PosDefinite=false;
  *NegDefinite=false;
  if(izero+ineg==A->Ni)
  {
    *NegDefinite=true;
  }
  else if(izero+ipos==A->Ni)
  {
    *PosDefinite=true;
  }
  TransformDiagonalMatrixByCompactMatrixTtranATcon(valpos, &fun, Apos);
  TransformDiagonalMatrixByCompactMatrixTtranATcon(valneg, &fun, Aneg);
  FreeCompactMatrix(&fun);
  free(val);
  free(valpos);
  free(valneg);
  return 0;
}




/*  Given a compact matrix A calculates  (w-A)^-1 _[[0,0]]
 *  algorithm assumes A to be Hermitian
 */
uint32_t HermitianCompactMatrixToListOfPoles(CompactMatrixType const *A, ListOfPoles *L)
{
  uint32_t i;
  CompactMatrixType fun;
  double *val;
  val= (double*)malloc(A->Ni*sizeof(double));
  CompactMatrixEigensystem(A, &fun, val);
  L->N=A->Ni;
  InitListOfPoles(L);
  L->mu=A->mu;
  L->a0=A->a0;
  memcpy(L->a, val, A->Ni*sizeof(double));
  if(fun.QComplex==0)
  {
    for(i=0;i<A->Nj;i++)
    {
      L->b[i]=Sqr(A->b0)*Sqr(fun.VRe[i*fun.Nj]);
    }
  }
  else
  {
    for(i=0;i<A->Nj;i++)
    {
      L->b[i]=Sqr(A->b0)*(Sqr(fun.V[i*fun.Nj].Re)+Sqr(fun.V[i*fun.Nj].Im));
    }
  }
  free(val);
  FreeCompactMatrix(&fun);
  return 0;
}



/*  Given a compact matrix A calculates vec[i]^* . (w-A)^-1 . vec[i]
 *  returns a list of vec.Ni ListOfPoles
 *  algorithm assumes A to be Hermitian
 */
uint32_t HermitianCompactMatrixToListOfPoles(CompactMatrixType const *A, ListOfPoles *L, CompactMatrixType *vec)
{
  uint32_t i,j;
  CompactMatrixType fun,funvec;
  double *val;
  val= (double*)malloc(A->Ni*sizeof(double));
  CompactMatrixEigensystem(A, &fun, val);
  CompactMatrixMultiply('N', 'C', vec, &fun, &funvec);
  FreeCompactMatrix(&fun);
  for(i=0;i<vec->Ni;i++)
  {
    L[i].mu=A->mu;
    L[i].a0=A->a0;
    L[i].N=A->Ni;
    InitListOfPoles(&(L[i]));
    memcpy(L[i].a, val, A->Ni*sizeof(double));
    if(funvec.QComplex==0)
    {
      for(j=0;j<A->Nj;j++)
      {
        L[i].b[j]=A->b0*Sqr(funvec.VRe[i*funvec.Nj+j]);
      }
    }
    else
    {
      for(j=0;j<A->Nj;j++)
      {
        L[i].b[j]=A->b0*(Sqr(funvec.V[i*funvec.Nj+j].Re)+Sqr(funvec.V[i*funvec.Nj+j].Im));
      }
    }
  }
  free(val);
  FreeCompactMatrix(&funvec);
  return 0;
}
/* fun allocated but not initialized */
uint32_t HermitianCompactMatrixToListOfPolesWithEigenstates(CompactMatrixType const *A, ListOfPoles *L, CompactMatrixType *vec, CompactMatrixType *fun)
{
  uint32_t i,j;
  CompactMatrixType funvec;
  double *val;
  val= (double*)malloc(A->Ni*sizeof(double));
  CompactMatrixEigensystem(A, fun, val);
  CompactMatrixMultiply('N', 'C', vec, fun, &funvec);
  for(i=0;i<vec->Ni;i++)
  {
    L[i].mu=A->mu;
    L[i].a0=A->a0;
    L[i].N=A->Ni;
    InitListOfPoles(&(L[i]));
    memcpy(L[i].a, val, A->Ni*sizeof(double));
    if(funvec.QComplex==0)
    {
      for(j=0;j<A->Nj;j++)
      {
        L[i].b[j]=A->b0*Sqr(funvec.VRe[i*funvec.Nj+j]);
      }
    }
    else
    {
      for(j=0;j<A->Nj;j++)
      {
        L[i].b[j]=A->b0*(Sqr(funvec.V[i*funvec.Nj+j].Re)+Sqr(funvec.V[i*funvec.Nj+j].Im));
      }
    }
  }
  free(val);
  FreeCompactMatrix(&funvec);
  return 0;
}

/*  Given a compact matrix A calculates vec[j]^* . (w-A)^-1 . vec[i]^T
 *  returns a BlockListOfPoles of Blockdimension vec.Ni
 *  algorithm assumes A to be Hermitian
 */
uint32_t HermitianCompactMatrixToBlockListOfPoles(CompactMatrixType const *A, BlockListOfPoles *L, CompactMatrixType *vec)
{
  uint32_t k,l,i,j,ifunvec;
  CompactMatrixType fun,funvec;
  double *val;
  val= (double*)malloc(A->Ni*sizeof(double));
  CompactMatrixEigensystem(A, &fun, val);
  CompactMatrixMultiply('*', 'T', &fun, vec, &funvec);
  FreeCompactMatrix(&fun);
  L->BlockSize=vec->Ni;
  L->N=vec->Nj;
  L->QComplex=funvec.QComplex;
  L->mu=A->mu;
  InitBlockListOfPoles(L);
  l=0;
  ifunvec=0;
  if(funvec.QComplex==0)
  {
    for(k=0;k<funvec.Ni;k++)
    {
      L->VRe[l]=val[k];
      l++;
      for(i=0;i<L->BlockSize;i++)
      {
        for(j=0;j<L->BlockSize;j++)
        {
          L->VRe[l]=funvec.VRe[ifunvec+j]*funvec.VRe[ifunvec+i];
          l++;
        }
      }
      ifunvec+=funvec.Nj;
    }
  }
  else
  {
    for(k=0;k<funvec.Ni;k++)
    {
      L->V[l].Re=val[k];
      l++;
      for(i=0;i<L->BlockSize;i++)
      {
        for(j=0;j<L->BlockSize;j++)
        {
          L->V[l].Re= funvec.V[ifunvec+j].Re*funvec.V[ifunvec+i].Re+funvec.V[ifunvec+j].Im*funvec.V[ifunvec+i].Im;
          L->V[l].Im=-funvec.V[ifunvec+j].Re*funvec.V[ifunvec+i].Im+funvec.V[ifunvec+j].Im*funvec.V[ifunvec+i].Re;
          l++;
        }
      }
      ifunvec+=funvec.Nj;
    }
  }
  free(val);
  FreeCompactMatrix(&funvec);
  return 0;
}
uint32_t HermitianCompactMatrixToBlockListOfPolesWithEigenstates(CompactMatrixType const *A, BlockListOfPoles *L, CompactMatrixType *vec, CompactMatrixType *fun)
{
  uint32_t k,l,i,j,ifunvec;
  CompactMatrixType funvec;
  double *val;
  val= (double*)malloc(A->Ni*sizeof(double));
  CompactMatrixEigensystem(A, fun, val);
  CompactMatrixMultiply('*', 'T', fun, vec, &funvec);
  L->BlockSize=vec->Ni;
  L->N=vec->Nj;
  L->QComplex=funvec.QComplex;
  L->mu=A->mu;
  InitBlockListOfPoles(L);
  l=0;
  ifunvec=0;
  if(funvec.QComplex==0)
  {
    for(k=0;k<funvec.Ni;k++)
    {
      L->VRe[l]=val[k];
      l++;
      for(i=0;i<L->BlockSize;i++)
      {
        for(j=0;j<L->BlockSize;j++)
        {
          L->VRe[l]=funvec.VRe[ifunvec+j]*funvec.VRe[ifunvec+i];
          l++;
        }
      }
      ifunvec+=funvec.Nj;
    }
  }
  else
  {
    for(k=0;k<funvec.Ni;k++)
    {
      L->V[l].Re=val[k];
      l++;
      for(i=0;i<L->BlockSize;i++)
      {
        for(j=0;j<L->BlockSize;j++)
        {
          L->V[l].Re= funvec.V[ifunvec+j].Re*funvec.V[ifunvec+i].Re+funvec.V[ifunvec+j].Im*funvec.V[ifunvec+i].Im;
          L->V[l].Im=-funvec.V[ifunvec+j].Re*funvec.V[ifunvec+i].Im+funvec.V[ifunvec+j].Im*funvec.V[ifunvec+i].Re;
          l++;
        }
      }
      ifunvec+=funvec.Nj;
    }
  }
  free(val);
  FreeCompactMatrix(&funvec);
  return 0;
}
/*  assume that vec is an identity matrix with Nrow rows
 *
 */
uint32_t HermitianCompactMatrixToBlockListOfPoles(CompactMatrixType const *A, BlockListOfPoles *L, uint32_t NRow)
{
  CompactMatrixType vec;
  vec.Ni=NRow;
  vec.Nj=A->Nj;
  vec.QComplex=0;
  InitCompactMatrix(&vec);
  for(int i=0;i<NRow;i++) vec.VRe[i + i *vec.Nj] = 1;
  HermitianCompactMatrixToBlockListOfPoles(A, L, &vec);
  FreeCompactMatrix(&vec);
  return 0;
}







/*
 *
 */
uint32_t ResponseFunctionGetA0(ResponsefunctionType const *A, double *a0)
{
  if(A->Block==true) return 1;
  if((32|A->type)=='t')
  {
    *a0 = ((TriDiagonalMatrixType*)(A->data))->a[0];
  }
  else if((32|A->type)=='a')
  {
    *a0 = ((AndersonMatrixType*)(A->data))->a[0];
  }
  else if((32|A->type)=='l')
  {
    *a0 = ((ListOfPoles*)(A->data))->a0;
  }
  else if((32|A->type)=='n')
  {
    *a0 = ((BADoubleTriDiagonalMatrixType*)(A->data))->a0;
  }
  else if((32|A->type)=='d')
  {
    *a0 = ((DoubleTriDiagonalMatrixType*)(A->data))->a0;
  }
  else
  {
    printf("ERROR: Unrecognized type in ResponseFunctionGetA0\n");
    return 1;
  }
  return 0;
}


uint32_t BlockResponseFunctionGetBlockSize(ResponsefunctionType const *A)
{
  if(A->Block==false) return 1;
  switch ((32|A->type)) 
  {
    case 't':
      return ((BlockTriDiagonalMatrixType*)(A->data))->BlockSize;
    break;
    case 'a':
      return ((BlockAndersonMatrixType*)(A->data))->BlockSize;
    break;
    case 'l':
      return ((BlockListOfPoles*)(A->data))->BlockSize;
    break;
    case 'n':
      return ((BlockNaturalImpurityOrbitalMatrixType*)(A->data))->BlockSize;
    break;
    case 'd':
      return ((BlockDoubleTriDiagonalMatrixType*)(A->data))->BlockSize;
    break;
    default:
      return 0;
    break;
  }
}

/*
 *
 */
uint32_t BlockResponseFunctionGetA0(ResponsefunctionType const *A, CompactMatrixType * const A0)
{
  if(A->Block==false) return 1;
  if((32|A->type)=='t')
  {
    A0->QComplex=((BlockTriDiagonalMatrixType*)A->data)->QComplex;
    A0->Ni=((BlockTriDiagonalMatrixType*)A->data)->BlockSize;
    A0->Nj=A0->Ni;
    InitCompactMatrix(A0);
    A0->mu=((BlockTriDiagonalMatrixType*)A->data)->mu;
    if(A0->QComplex)
    {
      memcpy(A0->V, ((BlockTriDiagonalMatrixType*)A->data)->a, A0->Ni*A0->Nj*sizeof(Complex16));
    }
    else
    {
      memcpy(A0->VRe, ((BlockTriDiagonalMatrixType*)A->data)->aRe, A0->Ni*A0->Nj*sizeof(double));
    }
  }
  else if((32|A->type)=='a')
  {
    A0->QComplex=((BlockAndersonMatrixType*)A->data)->QComplex;
    A0->Ni=((BlockAndersonMatrixType*)A->data)->BlockSize;
    A0->Nj=A0->Ni;
    InitCompactMatrix(A0);
    A0->mu=((BlockAndersonMatrixType*)A->data)->mu;
    if(A0->QComplex)
    {
      memcpy(A0->V, ((BlockAndersonMatrixType*)A->data)->a, A0->Ni*A0->Nj*sizeof(Complex16));
    }
    else
    {
      memcpy(A0->VRe, ((BlockAndersonMatrixType*)A->data)->aRe, A0->Ni*A0->Nj*sizeof(double));
    }
  }
  else if((32|A->type)=='l')
  {
    A0->QComplex=((BlockListOfPoles*)A->data)->QComplex;
    A0->Ni=((BlockListOfPoles*)A->data)->BlockSize;
    A0->Nj=A0->Ni;
    InitCompactMatrix(A0);
    A0->mu=((BlockListOfPoles*)A->data)->mu;
    if(A0->QComplex)
    {
      memcpy(A0->V, ((BlockListOfPoles*)A->data)->a0, A0->Ni*A0->Nj*sizeof(Complex16));
    }
    else
    {
      memcpy(A0->VRe, ((BlockListOfPoles*)A->data)->a0Re, A0->Ni*A0->Nj*sizeof(double));
    }
  }
  else if((32|A->type)=='n')
  {
    uint32_t i;
    CompactMatrixCopy(((BlockNaturalImpurityOrbitalMatrixType*)A->data)->a0, A0);
    //TODO: test if this is correct
    if(A0->QComplex)
    {
      for(i=0;i<A0->Ni*A0->Nj;i++)
      {
        A0->V[i].Re += ((BlockNaturalImpurityOrbitalMatrixType*)A->data)->val.a[i].Re;
        A0->V[i].Re += ((BlockNaturalImpurityOrbitalMatrixType*)A->data)->con.a[i].Re;
        A0->V[i].Im += ((BlockNaturalImpurityOrbitalMatrixType*)A->data)->val.a[i].Im;
        A0->V[i].Im += ((BlockNaturalImpurityOrbitalMatrixType*)A->data)->con.a[i].Im;
      }
    }
    else
    {
      for(i=0;i<A0->Ni*A0->Nj;i++)
      {
        A0->VRe[i] += ((BlockNaturalImpurityOrbitalMatrixType*)A->data)->val.aRe[i];
        A0->VRe[i] += ((BlockNaturalImpurityOrbitalMatrixType*)A->data)->con.aRe[i];
      }
    }
  }
  else if((32|A->type)=='d')
  {
    CompactMatrixCopy(((BlockDoubleTriDiagonalMatrixType*)A->data)->a0, A0);
    A0->mu=((BlockDoubleTriDiagonalMatrixType*)A->data)->mu;
  }
  else
  {
    printf("ERROR: Unrecognized type in BlockResponseFunctionGetA0\n");
    return 1;
  }
  return 0;
}


/*
 *
 */
uint32_t ResponseFunctionGetB0(ResponsefunctionType const *A, double *b0)
{
  if(A->Block==true) return 1;
  /*  if N==0 return b0 = 0 else return the value */
  if((32|A->type)=='t')
  {
    *b0 = (((TriDiagonalMatrixType*)(A->data))->N == 0) ? 0 : ((TriDiagonalMatrixType*)(A->data))->b[0];
  }
  else if((32|A->type)=='a')
  {
    *b0 = (((AndersonMatrixType*)(A->data))->N == 0) ? 0 : ((AndersonMatrixType*)(A->data))->b[0];
  }
  else if((32|A->type)=='l')
  {
    *b0 = 0;
    for(int i=0;i<((ListOfPoles*)(A->data))->N;i++)
    {
      *b0+=((ListOfPoles*)(A->data))->b[i];
    }
    *b0 = sqrt(*b0);
  }
  else if((32|A->type)=='n')
  {
    *b0 = ((BADoubleTriDiagonalMatrixType*)(A->data))->b0;
  }
  else if((32|A->type)=='d')
  {
    *b0 = ((DoubleTriDiagonalMatrixType*)(A->data))->b0;
  }
  else
  {
    printf("ERORR: Unrecognized type in ResponseFunctionGetB0\n");
    return 1;
  }
  return 0;
}

/*
 *
 */
uint32_t BlockResponseFunctionGetB0(ResponsefunctionType const *A, CompactMatrixType * const B0)
{
  if(A->Block==false) return 1;
  if((32|A->type)=='t')
  {
    B0->QComplex=((BlockTriDiagonalMatrixType*)A->data)->QComplex;
    B0->Ni=((BlockTriDiagonalMatrixType*)A->data)->BlockSize;
    B0->Nj=B0->Ni;
    InitCompactMatrix(B0);
    B0->mu=((BlockTriDiagonalMatrixType*)A->data)->mu;
    if(B0->QComplex)
    {
      memcpy(B0->V, ((BlockTriDiagonalMatrixType*)A->data)->b, B0->Ni*B0->Nj*sizeof(Complex16));
    }
    else
    {
      memcpy(B0->VRe, ((BlockTriDiagonalMatrixType*)A->data)->bRe, B0->Ni*B0->Nj*sizeof(double));
    }
  }
  else if((32|A->type)=='a')
  {
    B0->QComplex=((BlockAndersonMatrixType*)A->data)->QComplex;
    B0->Ni=((BlockAndersonMatrixType*)A->data)->BlockSize;
    B0->Nj=B0->Ni;
    InitCompactMatrix(B0);
    B0->mu=((BlockAndersonMatrixType*)A->data)->mu;
    if(B0->QComplex)
    {
      memcpy(B0->V, ((BlockAndersonMatrixType*)A->data)->b, B0->Ni*B0->Nj*sizeof(Complex16));
    }
    else
    {
      memcpy(B0->VRe, ((BlockAndersonMatrixType*)A->data)->bRe, B0->Ni*B0->Nj*sizeof(double));
    }
  }
  else if((32|A->type)=='l')
  {
    BlockListOfPoles const * const L=(BlockListOfPoles*)A->data;
    CompactMatrixType SqrtB0B0;
    uint32_t i,k,iN,NBlock;
    B0->QComplex=L->QComplex;
    B0->Ni=L->BlockSize;
    B0->Nj=B0->Ni;
    InitCompactMatrix(B0);
    B0->mu=L->mu;
    /*  B0^* . B0^T is the total spectral weight, i.e. the sum of the weight of all poles.
     *  This is the same B0 as the one created when transforming to a block tridiagonal form:
     *  BlockListOfPolesToBlockTridiagonalMatrix creates the starting vectors psi as the square root of the weight of
     *  each pole and CompactMatrixBlockBandDiagonalize sets b0 = psi . psi_ortho^C. As psi_ortho is created by an SVD
     *  returning the polar factor (U V^T) this is the Hermitian positive semi definite square root of
     *  psi . psi^C = sum_i B_i^T. (the transpose enters as the starting vectors are created with
     *  CompactMatrixTranspose, which does not conjugate, and is needed for the convention G = A0 + B0^* G1 B0^T)
     *  (compare the non block version in ResponseFunctionGetB0: b0 = sqrt( sum_i b_i ) )
     */
    NBlock=B0->Ni*B0->Nj;
    iN=0;
    for(k=0;k<L->N;k++)
    {
      if(B0->QComplex)
      {
        for(i=0;i<NBlock;i++) B0->V[i]  +=L->V[iN+1+i];
      }
      else
      {
        for(i=0;i<NBlock;i++) B0->VRe[i]+=L->VRe[iN+1+i];
      }
      iN+=1+NBlock;
    }
    CompactMatrixTranspose(B0);
    /*  the weight of a pole of a physical response function is Hermitian. Non Hermitian weights, as created by the
     *  non Hermitian eigenstates in CalculateGBathFromGandSigmaEigensystem, have no Hermitian square root
     */
    SquareCompactMatrixMakeHermitian(B0);
    CompactMatrixSqrt(B0, &SqrtB0B0);
    if(B0->QComplex)
    {
      memcpy(B0->V, SqrtB0B0.V, NBlock*sizeof(Complex16));
    }
    else
    {
      memcpy(B0->VRe, SqrtB0B0.VRe, NBlock*sizeof(double));
    }
    FreeCompactMatrix(&SqrtB0B0);
  }
  else if((32|A->type)=='n')
  {
    BlockNaturalImpurityOrbitalMatrixType * const M=(BlockNaturalImpurityOrbitalMatrixType*)A->data;
    CompactMatrixType T,Bcoup,B0B0,SqrtB0B0;
    B0->QComplex=M->QComplex;
    B0->Ni=M->BlockSize;
    B0->Nj=B0->Ni;
    InitCompactMatrix(B0);
    B0->mu=M->mu;
    GetTFromValandConInBlockNaturalImpurityOrbitalMatrixType(M, &T);
    Bcoup=T;              
    Bcoup.Ni=M->BlockSize;
    CompactMatrixMultiply('N', 'C', &Bcoup, &Bcoup, &B0B0);
    SquareCompactMatrixMakeHermitian(&B0B0);
    CompactMatrixSqrt(&B0B0, &SqrtB0B0);
    if(B0->QComplex)
    { /* val and con can be real while the natural impurity orbital matrix is complex */
      MakeCompactRealMatrixComplex(&SqrtB0B0);
      memcpy(B0->V, SqrtB0B0.V, B0->Ni*B0->Nj*sizeof(Complex16));
    }
    else
    {
      MakeCompactComplexMatrixReal(&SqrtB0B0);
      memcpy(B0->VRe, SqrtB0B0.VRe, B0->Ni*B0->Nj*sizeof(double));
    }
    FreeCompactMatrix(&SqrtB0B0);
    FreeCompactMatrix(&B0B0);
    FreeCompactMatrix(&T);
  }
  else if((32|A->type)=='d')
  {
    CompactMatrixCopy(((BlockDoubleTriDiagonalMatrixType*)A->data)->b0, B0);
    B0->mu=((BlockDoubleTriDiagonalMatrixType*)A->data)->mu;
  }
  else
  {
    printf("ERROR in BlockResponseFunctionGetB0 unrecognized type");
    return 1;
  }
  return 0;
}
/*  convert a response function to a compact matrix
 *  A0 and B0 are neglected
 *
 */
uint32_t ResponseFunctionToCompactMatrix(ResponsefunctionType const * A, CompactMatrixType * const B, const bool deflate, const double SingularValue)
{
  if(A->Block==true)
  { // block version
    if((32|A->type)=='t')
    { // convert block tridiagonal matrix to compact matrix
      if(deflate==true)
      {
        BlockTriDiagonalMatrixToCompactMatrixWithDeflation((BlockTriDiagonalMatrixType*)A->data, B);
      }
      else
      {
        BlockTriDiagonalMatrixToCompactMatrix((BlockTriDiagonalMatrixType*)A->data, B);
      }
      return 0;
    }
    else if((32|A->type)=='a')
    { // convert block anderson matrix to compact matrix
      BlockAndersonMatrixToCompactMatrix((BlockAndersonMatrixType*)A->data, B);      
      return 0;
    }
    else if((32|A->type)=='l')
    { // convert block list of poles response function to compact matrix
      if(deflate==true)
      {
        BlockListOfPolesToCompactMatrixWithDeflation((BlockListOfPoles*)A->data, B, SingularValue);
      }
      else
      {
        BlockListOfPolesToCompactMatrix((BlockListOfPoles*)A->data, B, SingularValue);
      }
      return 0;
    }
    else if((32|A->type)=='n')
    { // convert block natural impurity matrix to compact matrix
      if(deflate==true)
      {
        BlockNaturalImpurityOrbitalMatrixToCompactMatrixWithDeflation((BlockNaturalImpurityOrbitalMatrixType*)A->data, B);
      }
      else
      {
        BlockNaturalImpurityOrbitalMatrixToCompactMatrix((BlockNaturalImpurityOrbitalMatrixType*)A->data, B);
      }
      return 0;
    }
    else if((32|A->type)=='d')
    { // convert block double tridiagonal matrix to compact matrix
      if(deflate==true)
      {
        return BlockDoubleTriDiagonalMatrixToCompactMatrixWithDeflation((BlockDoubleTriDiagonalMatrixType*)(A->data), B);
      }
      else
      {
        return BlockDoubleTriDiagonalMatrixToCompactMatrix((BlockDoubleTriDiagonalMatrixType*)(A->data), B);
      }
    }
  }
  else
  { // non-block version
    if((32|A->type)=='t')
    { // convert tridiagonal matrix to compact matrix
      TriDiagonalMatrixToCompactMatrix( (TriDiagonalMatrixType*)(A->data), B);
      return 0;
    }
    if((32|A->type)=='a')
    { // convert anderson matrix to compact matrix
      AndersonMatrixToCompactMatrix((AndersonMatrixType*)(A->data), B);
      return 0;
    }
    if((32|A->type)=='l')
    { // convert list of poles response function to compact matrix
      ListOfPolesToCompactMatrix((ListOfPoles*)(A->data), B);
      return 0;
    }
    if((32|A->type)=='n')
    { // convert natural impurity matrix to compact matrix
      BADoubleTriDiagonalMatrixTypeToCompactMatrix((BADoubleTriDiagonalMatrixType*)(A->data), B);
      return 0;
    }
    if((32|A->type)=='d')
    {
      return DoubleTriDiagonalMatrixToCompactMatrix( *((DoubleTriDiagonalMatrixType*)(A->data)), B);
    }
  }
  printf("Error: Unidentified type of response function %c in ResponseFunctionToCompactMatrix\n",A->type);
  return 1;
}

/*
 *
 */
uint32_t BlockListOfPolesToCompactMatrix(BlockListOfPoles const *A, CompactMatrixType * const B, const double SingularValue)
{
  BlockTriDiagonalMatrixType temp;
  BlockListOfPolesToBlockTridiagonalMatrix(*A, &temp, A->N, SingularValue);
  BlockTriDiagonalMatrixToCompactMatrix(&temp, B);
  FreeBlockTridiagonalMatrix(&temp);
  return 0;
}
uint32_t BlockListOfPolesToCompactMatrixWithDeflation(BlockListOfPoles const *A, CompactMatrixType * const B, const double SingularValue)
{
  BlockTriDiagonalMatrixType temp;
  BlockListOfPolesToBlockTridiagonalMatrix(*A, &temp, A->N, SingularValue);
  BlockTriDiagonalMatrixToCompactMatrixWithDeflation(&temp, B);
  FreeBlockTridiagonalMatrix(&temp);
  return 0;
}


/*
 *
 */
uint32_t BlockNaturalImpurityOrbitalMatrixToCompactMatrix(BlockNaturalImpurityOrbitalMatrixType const * A, CompactMatrixType * const B)
{
  if((A->val.N==0)&&(A->con.N==0))
  {
    B->Ni=0;
    B->Nj=0;
    B->QComplex=0;
    InitCompactMatrix(B);
    B->mu=A->mu;
    return 0;
  }
  if(A->val.N==0)
  {
    BlockTriDiagonalMatrixToCompactMatrix(&(A->con), B);
    return 0;
  }
  if(A->con.N==0)
  {
    BlockTriDiagonalMatrixToCompactMatrix(&(A->val), B);
    return 0;
  }
  int i,j,k,apos,bposL,bposR;
  B->QComplex=A->QComplex;
  B->Ni=(Max(1,A->val.N)+Max(1,A->con.N))*A->BlockSize;
  B->Nj=(Max(1,A->val.N)+Max(1,A->con.N))*A->BlockSize;
  if(InitCompactMatrix(B))
  {
    printf("InitCompactMatrix failed in BlockNaturalImpurityOrbitalMatrixToCompactMatrix\n");
    return 1;
  }
  B->mu=A->mu;
  if(A->QComplex)
  {
    if((A->val.N==1)&&(A->con.N==1))
    {
      for(i=0;i<2*A->BlockSize;i++)
      {
        memcpy(&(B->V[i*B->Nj]), &(A->TMT.V[i*A->TMT.Nj]), 2*A->BlockSize*sizeof(Complex16));
      }
    }
    else if (A->val.N==1)
    {
      for(i=0;i<2*A->BlockSize;i++)
      {
        memcpy(&(B->V[i*B->Nj]), &(A->TMT.V[i*A->TMT.Nj]), 2*A->BlockSize*sizeof(Complex16));
      }
      for(i=0;i<A->BlockSize;i++)
      {
        memcpy(&(B->V[(i+A->BlockSize*Max(2,A->val.N+1))*B->Nj]), &(A->TMT.V[(i+3*A->BlockSize)*A->TMT.Nj]), 2*A->BlockSize*sizeof(Complex16));
      }
      for(i=0;i<2*A->BlockSize;i++)
      {
        memcpy(&(B->V[(i*B->Nj+A->BlockSize*Max(2,A->val.N+1))]), &(A->TMT.V[(i)*A->TMT.Nj+3*A->BlockSize]), A->BlockSize*sizeof(Complex16));
      }
      for(k=2;k<=A->val.N;k++)
      {
        for(i=0;i<A->BlockSize;i++)
        {
          memcpy(&(B->V[k*A->BlockSize*B->Nj+k*A->BlockSize+i*B->Nj]), &(A->val.a[k*A->BlockSize*A->BlockSize+i*A->BlockSize]), A->BlockSize*sizeof(Complex16));
        }
      }
      for(k=2;k<=A->con.N;k++)
      {
        for(i=0;i<A->BlockSize;i++)
        {
          memcpy(&(B->V[(k+Max(1,A->val.N)-1)*A->BlockSize*B->Nj+(k+Max(1,A->val.N)-1)*A->BlockSize+i*B->Nj]), &(A->con.a[k*A->BlockSize*A->BlockSize+i*A->BlockSize]), A->BlockSize*sizeof(Complex16));
        }
      }
    }
    else if (A->con.N==1)
    {
      for(i=0;i<3*A->BlockSize;i++)
      {
        memcpy(&(B->V[i*B->Nj]), &(A->TMT.V[i*A->TMT.Nj]), 3*A->BlockSize*sizeof(Complex16));
      }
      /*  TMT reaches as far as the second site of each chain, so the onsite energies of the
       *  val chain from there on have to follow. Their couplings are filled by the two
       *  loops after this switch, which run whichever branch was taken. */
      for(k=2;k<=A->val.N;k++)
      {
        for(i=0;i<A->BlockSize;i++)
        {
          memcpy(&(B->V[k*A->BlockSize*B->Nj+k*A->BlockSize+i*B->Nj]), &(A->val.a[k*A->BlockSize*A->BlockSize+i*A->BlockSize]), A->BlockSize*sizeof(Complex16));
        }
      }
    }
    else
    {
      for(i=0;i<3*A->BlockSize;i++)
      {
        memcpy(&(B->V[i*B->Nj]), &(A->TMT.V[i*A->TMT.Nj]), 3*A->BlockSize*sizeof(Complex16));
      }
      for(i=0;i<A->BlockSize;i++)
      {
        memcpy(&(B->V[(i+A->BlockSize*Max(2,A->val.N+1))*B->Nj]), &(A->TMT.V[(i+3*A->BlockSize)*A->TMT.Nj]), 2*A->BlockSize*sizeof(Complex16));
      }
      for(i=0;i<2*A->BlockSize;i++)
      {
        memcpy(&(B->V[(i*B->Nj+A->BlockSize*Max(2,A->val.N+1))]), &(A->TMT.V[(i)*A->TMT.Nj+3*A->BlockSize]), A->BlockSize*sizeof(Complex16));
      }
      for(k=2;k<=A->val.N;k++)
      {
        for(i=0;i<A->BlockSize;i++)
        {
          memcpy(&(B->V[k*A->BlockSize*B->Nj+k*A->BlockSize+i*B->Nj]), &(A->val.a[k*A->BlockSize*A->BlockSize+i*A->BlockSize]), A->BlockSize*sizeof(Complex16));
        }
      }
      for(k=2;k<=A->con.N;k++)
      {
        for(i=0;i<A->BlockSize;i++)
        {
          memcpy(&(B->V[(k+Max(1,A->val.N)-1)*A->BlockSize*B->Nj+(k+Max(1,A->val.N)-1)*A->BlockSize+i*B->Nj]), &(A->con.a[k*A->BlockSize*A->BlockSize+i*A->BlockSize]), A->BlockSize*sizeof(Complex16));
        }
      }
    }
    apos=2*A->BlockSize*A->BlockSize;
    bposL=2*A->BlockSize*B->Nj+3*A->BlockSize;
    bposR=3*A->BlockSize*B->Nj+2*A->BlockSize;
    for(k=2;k<A->val.N;k++)
    {
      for(i=0;i<A->BlockSize;i++)
      {
        for(j=0;j<A->BlockSize;j++)
        {
          B->V[bposL].Re= A->val.b[apos].Re;
          B->V[bposR].Re= A->val.b[apos].Re;
          B->V[bposL].Im=-A->val.b[apos].Im;
          B->V[bposR].Im= A->val.b[apos].Im;
          bposR++;
          bposL+=B->Nj;
          apos++;
        }
        bposR+=B->Nj-A->BlockSize;
        bposL+=(1-A->BlockSize*B->Nj);
      }
      bposR+=A->BlockSize;
      bposL+=(A->BlockSize*B->Nj);
    }
    apos=2*A->BlockSize*A->BlockSize;
    k=1+Max(1,A->val.N);
    bposL=k*A->BlockSize*B->Nj+(k+1)*A->BlockSize;
    bposR=(k+1)*A->BlockSize*B->Nj+k*A->BlockSize;
    for(k=2;k<A->con.N;k++)
    {
      for(i=0;i<A->BlockSize;i++)
      {
        for(j=0;j<A->BlockSize;j++)
        {
          B->V[bposL].Re= A->con.b[apos].Re;
          B->V[bposR].Re= A->con.b[apos].Re;
          B->V[bposL].Im=-A->con.b[apos].Im;
          B->V[bposR].Im= A->con.b[apos].Im;
          bposR++;
          bposL+=B->Nj;
          apos++;
        }
        bposR+=B->Nj-A->BlockSize;
        bposL+=(1-A->BlockSize*B->Nj);
      }
      bposR+=A->BlockSize;
      bposL+=(A->BlockSize*B->Nj);
    }
  }
  else //real
  {
    if((A->val.N==1)&&(A->con.N==1))
    {
      for(i=0;i<2*A->BlockSize;i++)
      {
        memcpy(&(B->VRe[i*B->Nj]), &(A->TMT.VRe[i*A->TMT.Nj]), 2*A->BlockSize*sizeof(double));
      }
    }
    else if (A->val.N==1)
    {
      for(i=0;i<2*A->BlockSize;i++)
      {
        memcpy(&(B->VRe[i*B->Nj]), &(A->TMT.VRe[i*A->TMT.Nj]), 2*A->BlockSize*sizeof(double));
      }
      for(i=0;i<A->BlockSize;i++)
      {
        memcpy(&(B->VRe[(i+A->BlockSize*Max(2,A->val.N+1))*B->Nj]), &(A->TMT.VRe[(i+3*A->BlockSize)*A->TMT.Nj]), 2*A->BlockSize*sizeof(double));
      }
      for(i=0;i<2*A->BlockSize;i++)
      {
        memcpy(&(B->VRe[(i*B->Nj+A->BlockSize*Max(2,A->val.N+1))]), &(A->TMT.VRe[(i)*A->TMT.Nj+3*A->BlockSize]), A->BlockSize*sizeof(double));
      }
      for(k=2;k<=A->val.N;k++)
      {
        for(i=0;i<A->BlockSize;i++)
        {
          memcpy(&(B->VRe[k*A->BlockSize*B->Nj+k*A->BlockSize+i*B->Nj]), &(A->val.aRe[k*A->BlockSize*A->BlockSize+i*A->BlockSize]), A->BlockSize*sizeof(double));
        }
      }
      for(k=2;k<=A->con.N;k++)
      {
        for(i=0;i<A->BlockSize;i++)
        {
          memcpy(&(B->VRe[(k+Max(1,A->val.N)-1)*A->BlockSize*B->Nj+(k+Max(1,A->val.N)-1)*A->BlockSize+i*B->Nj]), &(A->con.aRe[k*A->BlockSize*A->BlockSize+i*A->BlockSize]), A->BlockSize*sizeof(double));
        }
      }
    }
    else if (A->con.N==1)
    {
      for(i=0;i<3*A->BlockSize;i++)
      {
        memcpy(&(B->VRe[i*B->Nj]), &(A->TMT.VRe[i*A->TMT.Nj]), 3*A->BlockSize*sizeof(double));
      }
      /*  TMT reaches as far as the second site of each chain, so the onsite energies of the
       *  val chain from there on have to follow. Their couplings are filled by the two
       *  loops after this switch, which run whichever branch was taken. */
      for(k=2;k<=A->val.N;k++)
      {
        for(i=0;i<A->BlockSize;i++)
        {
          memcpy(&(B->VRe[k*A->BlockSize*B->Nj+k*A->BlockSize+i*B->Nj]), &(A->val.aRe[k*A->BlockSize*A->BlockSize+i*A->BlockSize]), A->BlockSize*sizeof(double));
        }
      }
    }
    else
    {
      for(i=0;i<3*A->BlockSize;i++)
      {
        memcpy(&(B->VRe[i*B->Nj]), &(A->TMT.VRe[i*A->TMT.Nj]), 3*A->BlockSize*sizeof(double));
      }
      for(i=0;i<A->BlockSize;i++)
      {
        memcpy(&(B->VRe[(i+A->BlockSize*Max(2,A->val.N+1))*B->Nj]), &(A->TMT.VRe[(i+3*A->BlockSize)*A->TMT.Nj]), 2*A->BlockSize*sizeof(double));
      }
      for(i=0;i<2*A->BlockSize;i++)
      {
        memcpy(&(B->VRe[(i*B->Nj+A->BlockSize*Max(2,A->val.N+1))]), &(A->TMT.VRe[(i)*A->TMT.Nj+3*A->BlockSize]), A->BlockSize*sizeof(double));
      }
      for(k=2;k<=A->val.N;k++)
      {
        for(i=0;i<A->BlockSize;i++)
        {
          memcpy(&(B->VRe[k*A->BlockSize*B->Nj+k*A->BlockSize+i*B->Nj]), &(A->val.aRe[k*A->BlockSize*A->BlockSize+i*A->BlockSize]), A->BlockSize*sizeof(double));
        }
      }
      for(k=2;k<=A->con.N;k++)
      {
        for(i=0;i<A->BlockSize;i++)
        {
          memcpy(&(B->VRe[(k+Max(1,A->val.N)-1)*A->BlockSize*B->Nj+(k+Max(1,A->val.N)-1)*A->BlockSize+i*B->Nj]), &(A->con.aRe[k*A->BlockSize*A->BlockSize+i*A->BlockSize]), A->BlockSize*sizeof(double));
        }
      }
    }
    apos=2*A->BlockSize*A->BlockSize;
    bposL=2*A->BlockSize*B->Nj+3*A->BlockSize;
    bposR=3*A->BlockSize*B->Nj+2*A->BlockSize;
    for(k=2;k<A->val.N;k++)
    {
      for(i=0;i<A->BlockSize;i++)
      {
        for(j=0;j<A->BlockSize;j++)
        {
          B->VRe[bposL]=A->val.bRe[apos];
          B->VRe[bposR]=A->val.bRe[apos];
          bposR++;
          bposL+=B->Nj;
          apos++;
        }
        bposR+=B->Nj-A->BlockSize;
        bposL+=(1-A->BlockSize*B->Nj);
      }
      bposR+=A->BlockSize;
      bposL+=(A->BlockSize*B->Nj);
    }
    apos=2*A->BlockSize*A->BlockSize;
    k=1+Max(1,A->val.N);
    bposL=k*A->BlockSize*B->Nj+(k+1)*A->BlockSize;
    bposR=(k+1)*A->BlockSize*B->Nj+k*A->BlockSize;
    for(k=2;k<A->con.N;k++)
    {
      for(i=0;i<A->BlockSize;i++)
      {
        for(j=0;j<A->BlockSize;j++)
        {
          B->VRe[bposL]=A->con.bRe[apos];
          B->VRe[bposR]=A->con.bRe[apos];
          bposR++;
          bposL+=B->Nj;
          apos++;
        }
        bposR+=B->Nj-A->BlockSize;
        bposL+=(1-A->BlockSize*B->Nj);
      }
      bposR+=A->BlockSize;
      bposL+=(A->BlockSize*B->Nj);
    }
  }
  return 0;
}


/*
 *
 */
uint32_t BlockNaturalImpurityOrbitalMatrixToCompactMatrixWithDeflation(BlockNaturalImpurityOrbitalMatrixType const * A, CompactMatrixType * const B)
{
  //TODO: BlockNaturalImpurityOrbitalMatrixToCompactMatrixWithDeflation not yet implemented use function without deflation
  BlockNaturalImpurityOrbitalMatrixToCompactMatrix(A,B);
  return 1;
}

/*
 *
 */
uint32_t BlockAndersonMatrixToCompactMatrix(BlockAndersonMatrixType const * A, CompactMatrixType * const B)
{
  uint32_t i,j,Ai,BBi,BAi,BBdagi,BBdagiold,BBiold;
  int k;
  B->QComplex=A->QComplex;
  B->Ni = (A->N) * A->BlockSize;
  B->Nj = (A->N) * A->BlockSize;
  InitCompactMatrix(B);
  B->mu=A->mu;
  strcpy(B->Name, A->Name);
  if(A->N==0) return 0;
  Ai     = A->BlockSize*A->BlockSize;
  BBi    = A->BlockSize;
  BBdagi = A->BlockSize * B->Nj;
  BAi    = 0;
  if(A->QComplex)
  {
    for(k=0;k<((int)A->N)-1;k++)
    {
      BBiold = BBi;
      for(i=0;i<A->BlockSize;i++)
      {
        BBdagiold = BBdagi;
        for(j=0;j<A->BlockSize;j++)
        {
          B->V[BAi]     = A->a[Ai];
          B->V[BBi]     = A->b[Ai];
          B->V[BBdagi].Re  = A->b[Ai].Re;
          B->V[BBdagi].Im  =-A->b[Ai].Im;
          BBi++;
          BAi++;
          BBdagi+=B->Nj;
          Ai++;
        }
        BAi    += B->Nj-A->BlockSize;
        BBi    += B->Nj-A->BlockSize;
        BBdagi  = BBdagiold+1;
      }
      BAi    += A->BlockSize;
      BBi     = BBiold + A->BlockSize;
      BBdagi += A->BlockSize * B->Nj - A->BlockSize;
    }
    for(k=Max(0,((int)A->N)-1);k<A->N;k++)
    {
      BBiold = BBi;
      for(i=0;i<A->BlockSize;i++)
      {
        BBdagiold = BBdagi;
        for(j=0;j<A->BlockSize;j++)
        {
          B->V[BAi]     = A->a[Ai];
          BAi++;
          Ai++;
        }
        BAi    += B->Nj-A->BlockSize;
      }
      BAi    += A->BlockSize;
    }
  }
  else
  {
    for(k=0;k<(((int)A->N)-1);k++)
    {
      BBiold = BBi;
      for(i=0;i<A->BlockSize;i++)
      {
        BBdagiold = BBdagi;
        for(j=0;j<A->BlockSize;j++)
        {
          B->VRe[BAi]     = A->aRe[Ai];
          B->VRe[BBi]     = A->bRe[Ai];
          B->VRe[BBdagi]  = A->bRe[Ai];
          BBi++;
          BAi++;
          BBdagi+=B->Nj;
          Ai++;
        }
        BAi    += B->Nj-A->BlockSize;
        BBi    += B->Nj-A->BlockSize;
        BBdagi  = BBdagiold+1;
      }
      BAi    += A->BlockSize;
      BBi     = BBiold + A->BlockSize;
      BBdagi += A->BlockSize * B->Nj - A->BlockSize;
    }
    for(k=Max(0,((int)A->N)-1);k<A->N;k++)
    {
      BBiold = BBi;
      for(i=0;i<A->BlockSize;i++)
      {
        BBdagiold = BBdagi;
        for(j=0;j<A->BlockSize;j++)
        {
          B->VRe[BAi]     = A->aRe[Ai];
          BAi++;
          Ai++;
        }
        BAi    += B->Nj-A->BlockSize;
      }
      BAi    += A->BlockSize;
    }
  }
  return 0;
}


/*
 *
 */
uint32_t BlockTriDiagonalMatrixToCompactMatrix(BlockTriDiagonalMatrixType const * A, CompactMatrixType * const B)
{
  int i,j,k,apos,bpos,bposL,bposR;
  B->QComplex=A->QComplex;
  B->Ni=A->N * A->BlockSize;
  B->Nj=A->N * A->BlockSize;
  if(InitCompactMatrix(B))
  {
    printf("InitCompactMatrix failed in BlockTriDiagonalMatrixToCompactMatrix\n");
    return 1;
  }
  strcpy(B->Name, A->Name);
  B->mu=A->mu;
  bpos=0;
  apos=A->BlockSize*A->BlockSize;
  if(A->QComplex)
  {
    for(k=1;k<=A->N;k++)
    {
      for(i=0;i<A->BlockSize;i++)
      {
        for(j=0;j<A->BlockSize;j++)
        {
          B->V[bpos].Re=A->a[apos].Re;
          B->V[bpos].Im=A->a[apos].Im;
          bpos++;
          apos++;
        }
        bpos+=B->Nj-A->BlockSize;
      }
      bpos+=A->BlockSize;
    }
    apos=A->BlockSize*A->BlockSize;
    bposL=A->BlockSize;
    bposR=A->BlockSize*B->Nj;
    for(k=1;k<A->N;k++)
    {
      for(i=0;i<A->BlockSize;i++)
      {
        for(j=0;j<A->BlockSize;j++)
        {
          B->V[bposL].Re= A->b[apos].Re;
          B->V[bposR].Re= A->b[apos].Re;
          B->V[bposL].Im= A->b[apos].Im;
          B->V[bposR].Im=-A->b[apos].Im;
          bposL++;
          bposR+=B->Nj;
          apos++;
        }
        bposL+=B->Nj-A->BlockSize;
        bposR+=(1-A->BlockSize*B->Nj);
      }
      bposL+=A->BlockSize;
      bposR+=(A->BlockSize*B->Nj);
    }
  }
  else
  {
    for(k=1;k<=A->N;k++)
    {
      for(i=0;i<A->BlockSize;i++)
      {
        for(j=0;j<A->BlockSize;j++)
        {
          B->VRe[bpos]=A->aRe[apos];
          bpos++;
          apos++;
        }
        bpos+=B->Nj-A->BlockSize;
      }
      bpos+=A->BlockSize;
    }
    apos=A->BlockSize*A->BlockSize;
    bposL=A->BlockSize;
    bposR=A->BlockSize*B->Nj;
    for(k=1;k<A->N;k++)
    {
      for(i=0;i<A->BlockSize;i++)
      {
        for(j=0;j<A->BlockSize;j++)
        {
          B->VRe[bposL]=A->bRe[apos];
          B->VRe[bposR]=A->bRe[apos];
          bposL++;
          bposR+=B->Nj;
          apos++;
        }
        bposL+=B->Nj-A->BlockSize;
        bposR+=(1-A->BlockSize*B->Nj);
      }
      bposL+=A->BlockSize;
      bposR+=(A->BlockSize*B->Nj);
    }
  }
  return 0;
}

/*  Turn a BlockTriDiagonalMatrix into a compact dense matrix
 *  The array BlockSizeList in A contains a list with the size of all sub matrices on the diagonal A[0] to A[N+1]
 *  A0 and B0 are neglected
 *
 */
uint32_t BlockTriDiagonalMatrixToCompactMatrixWithDeflation(BlockTriDiagonalMatrixType const * A, CompactMatrixType * const B)
{
  uint32_t i,j,k,apos,bpos,bposL,bposR;
  B->QComplex=A->QComplex;
  B->Ni=0;
  for(i=1;i<=A->N;i++)
  {
    B->Ni+=A->BlockSizeList[i];
  }
  B->Nj=B->Ni;
  if(InitCompactMatrix(B))
  {
    printf("InitCompactMatrix failed in BlockTriDiagonalMatrixToCompactMatrix\n");
    return 1;
  }
  B->mu=A->mu;
  // set matrix A
  bpos=0;
  apos=A->BlockSize*A->BlockSize; // start A at position A.a[1]
  if(A->QComplex)
  {
    for(k=1;k<=A->N;k++)
    {
      for(i=0;i<A->BlockSizeList[k];i++)
      {
        for(j=0;j<A->BlockSizeList[k];j++) // loop i and j over the block with the correct size for the dense matrix
        {
          B->V[bpos]=A->a[apos];
          bpos++;
          apos++;
        }
        apos+=A->BlockSize-A->BlockSizeList[k]; // at the end of a row for the deflated block move A to the end of the row of the undeflated block
        bpos+=B->Nj-A->BlockSizeList[k];        // move B to the begining of the next row
      }
      apos+=(A->BlockSize-A->BlockSizeList[k])*A->BlockSize; // at the end of the block (rows and columns) move A to the next block (skipe Ndeflated rows)
      bpos+=A->BlockSizeList[k];  // move B to the begining of the next Block
    }
    // set matrix B
    apos=A->BlockSize*A->BlockSize; // Start A at position A.b[1]
    bposL=A->BlockSizeList[1]; // you must love this naming. Set bposL to the first right top block
    bposR=A->BlockSizeList[1]*B->Nj; // Set bposR to the first left botom block
    for(k=1;k<A->N;k++)
    {
      for(i=0;i<A->BlockSizeList[k];i++)  // the right top block has BlockSize[k-1] rows
      {
        for(j=0;j<A->BlockSizeList[k+1];j++)  // the right top block has BlockSize[k] columns
        {
          B->V[bposL].Re= A->b[apos].Re;
          B->V[bposR].Re= A->b[apos].Re;
          B->V[bposL].Im= A->b[apos].Im;
          B->V[bposR].Im=-A->b[apos].Im;
          bposL++;
          bposR+=B->Nj;
          apos++;
        }
        apos+=A->BlockSize-A->BlockSizeList[k+1];
        bposL+=B->Nj-A->BlockSizeList[k+1];
        bposR+=(1-A->BlockSizeList[k+1]*B->Nj);
      }
      apos+=(A->BlockSize-A->BlockSizeList[k])*A->BlockSize;
      bposL+=A->BlockSizeList[k+1];
      bposR+=(A->BlockSizeList[k+1]*B->Nj);
    }
  }
  else
  {
    for(k=1;k<=A->N;k++)
    {
      for(i=0;i<A->BlockSizeList[k];i++)
      {
        for(j=0;j<A->BlockSizeList[k];j++) // loop i and j over the block with the correct size for the dense matrix
        {
          B->VRe[bpos]=A->aRe[apos];
          bpos++;
          apos++;
        }
        apos+=A->BlockSize-A->BlockSizeList[k]; // at the end of a row for the deflated block move A to the end of the row of the undeflated block
        bpos+=B->Nj-A->BlockSizeList[k];        // move B to the begining of the next row
      }
      apos+=(A->BlockSize-A->BlockSizeList[k])*A->BlockSize; // at the end of the block (rows and columns) move A to the next block (skipe Ndeflated rows)
      bpos+=A->BlockSizeList[k];  // move B to the begining of the next Block
    }
    // set matrix B
    apos=A->BlockSize*A->BlockSize; // Start A at position A.b[1]
    bposL=A->BlockSizeList[1]; // you must love this naming. Set bposL to the first right top block
    bposR=A->BlockSizeList[1]*B->Nj; // Set bposR to the first left botom block
    for(k=1;k<A->N;k++)
    {
      for(i=0;i<A->BlockSizeList[k];i++)  // the right top block has BlockSize[k-1] rows
      {
        for(j=0;j<A->BlockSizeList[k+1];j++)  // the right top block has BlockSize[k] columns
        {
          B->VRe[bposL]= A->bRe[apos];
          B->VRe[bposR]= A->bRe[apos];
          bposL++;
          bposR+=B->Nj;
          apos++;
        }
        apos+=A->BlockSize-A->BlockSizeList[k+1];
        bposL+=B->Nj-A->BlockSizeList[k+1];
        bposR+=(1-A->BlockSizeList[k+1]*B->Nj);
      }
      apos+=(A->BlockSize-A->BlockSizeList[k])*A->BlockSize;
      bposL+=A->BlockSizeList[k+1];
      bposR+=(A->BlockSizeList[k+1]*B->Nj);
    }
  }
  return 0;
}

/* Should go in BasicMath_CompactMatrixVector.c
 *
 * !!Not all cases thoroughly tested yet!!
 *
 * Converts the square section enclosed by iStart and iEnd to a CompactMatrix
 *
 * Useful for returning the bath part of a block tridiagonal matrix
 *
 */
uint32_t BlockTriDiagonalMatrixPartToCompactMatrix(BlockTriDiagonalMatrixType const * A, CompactMatrixType * const B, const uint32_t iStart, const uint32_t iEnd)
{
  int i,j,k,apos,bpos,bposL,bposR;
  B->QComplex=A->QComplex;
  B->Ni=(iEnd - iStart +1)* A->BlockSize;
  B->Nj=(iEnd - iStart +1) * A->BlockSize;
  if(InitCompactMatrix(B))
  {
    printf("InitCompactMatrix failed in BlockTriDiagonalMatrixToCompactMatrix\n");
    return 1;
  }
  bpos=0;
  apos=(1+iStart)*A->BlockSize*A->BlockSize;
  if(A->QComplex)
  {
    for(k=iStart;k<=iEnd;k++)
    {
      for(i=0;i<A->BlockSize;i++)
      {
        for(j=0;j<A->BlockSize;j++)
        {
          B->V[bpos].Re=A->a[apos].Re;
          B->V[bpos].Im=A->a[apos].Im;
          bpos++;
          apos++;
        }
        bpos+=B->Nj-A->BlockSize;
      }
      bpos+=A->BlockSize;
    }
    apos=(1+iStart)*A->BlockSize*A->BlockSize;
    bposR=A->BlockSize;
    bposL=A->BlockSize*B->Nj;
    for(k=iStart;k<iEnd;k++)
    {
      for(i=0;i<A->BlockSize;i++)
      {
        for(j=0;j<A->BlockSize;j++)
        {
          B->V[bposL].Re= A->b[apos].Re;
          B->V[bposR].Re= A->b[apos].Re;
          B->V[bposL].Im=-A->b[apos].Im;
          B->V[bposR].Im= A->b[apos].Im;
          bposR++;
          bposL+=B->Nj;
          apos++;
        }
        bposR+=B->Nj-A->BlockSize;
        bposL+=(1-A->BlockSize*B->Nj);
      }
      bposR+=A->BlockSize;
      bposL+=(A->BlockSize*B->Nj);
    }
  }
  else
  {
    for(k=iStart;k<=iEnd;k++)
    {
      for(i=0;i<A->BlockSize;i++)
      {
        for(j=0;j<A->BlockSize;j++)
        {
          B->VRe[bpos]=A->aRe[apos];
          bpos++;
          apos++;
        }
        bpos+=B->Nj-A->BlockSize;
      }
      bpos+=A->BlockSize;
    }
    apos=(1+iStart)*A->BlockSize*A->BlockSize;
    bposR=A->BlockSize;
    bposL=A->BlockSize*B->Nj;
    for(k=iStart;k<iEnd;k++)
    {
      for(i=0;i<A->BlockSize;i++)
      {
        for(j=0;j<A->BlockSize;j++)
        {
          B->VRe[bposL]=A->bRe[apos];
          B->VRe[bposR]=A->bRe[apos];
          bposR++;
          bposL+=B->Nj;
          apos++;
        }
        bposR+=B->Nj-A->BlockSize;
        bposL+=(1-A->BlockSize*B->Nj);
      }
      bposR+=A->BlockSize;
      bposL+=(A->BlockSize*B->Nj);
    }
  }
  return 0;
}

/*  Create the matrix whose resolvent is the Green's function represented by A: one impurity
 *  block with a valence and a conduction chain attached to it.
 *
 *  The dimension is BlockSize + dim(Mval) + dim(Mcon), with Mval and Mcon the matrices that
 *  BlockTriDiagonalMatrixToCompactMatrix[WithDeflation] returns for val and con. The impurity
 *  block is val.a[0] + con.a[0] (both chains share this block, see StructDef.cpp).
 *
 *          ( Aval0+Acon0 | Bval0^* | 0     | Bcon0^* | 0     )
 *      B = ( Bval0^T     |        Mval     | 0       | 0     )
 *          ( 0           |                 | 0       | 0     )
 *          ( Bcon0^T     | 0       | 0     |      Mcon       )
 *          ( 0           | 0       | 0     |                 )
 *
 *  The head coupling plays the role of b0 of its chain, which enters the Green's function as
 *  b0^* . G . b0^T and not as b0 . G . b0^+ the way the couplings inside a chain do (see
 *  BlockTriDiagonalMatrixType in StructDef.cpp). A hermitian matrix with X in the upper right
 *  block gives X . G . X^+, so X has to be b0^*.
 *
 *  A->a0 and A->b0 are the offset and the scaling of the Green's function and are not part of
 *  this matrix. Read them with BlockResponseFunctionGetA0 and ...GetB0, so that
 *    G(w) = A0 + B0^* [ 1/(w - B) ]_{0,0} B0^T
 *
 *  A dense matrix has no mixed representation, so B is real when val and con are both real
 *  and complex otherwise. A mixed A (QComplex 2) therefore lands on QComplex 1 here, unless
 *  it is only a0 or b0 that is complex, as those do not enter this matrix.
 *
 *  With deflation the first block of a chain can be smaller than BlockSize, in which case
 *  only the first BlockSizeList[1] columns of b[0] are non-zero, so only those are copied.
 *
 *  B not allocated on entry.
 */
uint32_t BlockDoubleTriDiagonalMatrixToCompactMatrixDeflate(BlockDoubleTriDiagonalMatrixType const * A, CompactMatrixType * const B, bool const deflate)
{
  uint32_t i,j,bs,nval1,ncon1,offset;
  double re,im;
  CompactMatrixType Mval, Mcon;
  bs = A->BlockSize;
  if (deflate)
  {
    if (BlockTriDiagonalMatrixToCompactMatrixWithDeflation(&(A->val), &Mval)) return 1;
    if (BlockTriDiagonalMatrixToCompactMatrixWithDeflation(&(A->con), &Mcon))
    {
      FreeCompactMatrix(&Mval);
      return 1;
    }
  }
  else
  {
    if (BlockTriDiagonalMatrixToCompactMatrix(&(A->val), &Mval)) return 1;
    if (BlockTriDiagonalMatrixToCompactMatrix(&(A->con), &Mcon))
    {
      FreeCompactMatrix(&Mval);
      return 1;
    }
  }
  /* how many orbitals the first block of each chain contributes */
  nval1 = 0;
  if (A->val.N != 0) nval1 = deflate ? A->val.BlockSizeList[1] : bs;
  ncon1 = 0;
  if (A->con.N != 0) ncon1 = deflate ? A->con.BlockSizeList[1] : bs;
  B->QComplex = ((A->val.QComplex!=0) || (A->con.QComplex!=0)) ? 1 : 0;
  B->Ni = bs + Mval.Ni + Mcon.Ni;
  B->Nj = B->Ni;
  if (InitCompactMatrix(B))
  {
    printf("InitCompactMatrix failed in BlockDoubleTriDiagonalMatrixToCompactMatrix\n");
    FreeCompactMatrix(&Mval);
    FreeCompactMatrix(&Mcon);
    return 1;
  }
  strncpy(B->Name, A->Name, MAXStrLength);
  B->Name[MAXStrLength-1]='\0';
  B->mu = A->mu;
  /* one chain can be real while the other is complex, promote so that both can be copied */
  if ((B->QComplex!=0) && (Mval.QComplex==0)) MakeCompactRealMatrixComplex(&Mval);
  if ((B->QComplex!=0) && (Mcon.QComplex==0)) MakeCompactRealMatrixComplex(&Mcon);
  /*  the impurity block, shared by the two chains. Each chain is read according to its own
   *  QComplex, which for a mixed A need not be the QComplex of A itself.
   */
  for (i = 0; i < bs; i++)
  {
    for (j = 0; j < bs; j++)
    {
      re = 0;
      im = 0;
      if (A->val.QComplex) { re += A->val.a[i*bs + j].Re; im += A->val.a[i*bs + j].Im; }
      else                 { re += A->val.aRe[i*bs + j]; }
      if (A->con.QComplex) { re += A->con.a[i*bs + j].Re; im += A->con.a[i*bs + j].Im; }
      else                 { re += A->con.aRe[i*bs + j]; }
      if (B->QComplex) { B->V[i*B->Nj + j].Re = re; B->V[i*B->Nj + j].Im = im; }
      else             { B->VRe[i*B->Nj + j] = re; }
    }
  }
  /* hopping impurity <-> first valence block, as b0^* in the upper right and b0^T below */
  for (i = 0; i < bs; i++)
  {
    for (j = 0; j < nval1; j++)
    {
      if (A->val.QComplex) { re = A->val.b[i*bs + j].Re; im = A->val.b[i*bs + j].Im; }
      else                 { re = A->val.bRe[i*bs + j];  im = 0; }
      if (B->QComplex)
      {
        B->V[i*B->Nj + bs+j].Re   =  re;
        B->V[i*B->Nj + bs+j].Im   = -im;
        B->V[(bs+j)*B->Nj + i].Re =  re;
        B->V[(bs+j)*B->Nj + i].Im =  im;
      }
      else
      {
        B->VRe[i*B->Nj + bs+j]   = re;
        B->VRe[(bs+j)*B->Nj + i] = re;
      }
    }
  }
  /* hopping impurity <-> first conduction block */
  offset = bs + Mval.Ni;
  for (i = 0; i < bs; i++)
  {
    for (j = 0; j < ncon1; j++)
    {
      if (A->con.QComplex) { re = A->con.b[i*bs + j].Re; im = A->con.b[i*bs + j].Im; }
      else                 { re = A->con.bRe[i*bs + j];  im = 0; }
      if (B->QComplex)
      {
        B->V[i*B->Nj + offset+j].Re   =  re;
        B->V[i*B->Nj + offset+j].Im   = -im;
        B->V[(offset+j)*B->Nj + i].Re =  re;
        B->V[(offset+j)*B->Nj + i].Im =  im;
      }
      else
      {
        B->VRe[i*B->Nj + offset+j]   = re;
        B->VRe[(offset+j)*B->Nj + i] = re;
      }
    }
  }
  /* the two chains */
  for (i = 0; i < Mval.Ni; i++)
  {
    for (j = 0; j < Mval.Nj; j++)
    {
      if (B->QComplex) B->V[(bs+i)*B->Nj + bs+j] = Mval.V[i*Mval.Nj + j];
      else             B->VRe[(bs+i)*B->Nj + bs+j] = Mval.VRe[i*Mval.Nj + j];
    }
  }
  for (i = 0; i < Mcon.Ni; i++)
  {
    for (j = 0; j < Mcon.Nj; j++)
    {
      if (B->QComplex) B->V[(offset+i)*B->Nj + offset+j] = Mcon.V[i*Mcon.Nj + j];
      else             B->VRe[(offset+i)*B->Nj + offset+j] = Mcon.VRe[i*Mcon.Nj + j];
    }
  }
  FreeCompactMatrix(&Mval);
  FreeCompactMatrix(&Mcon);
  return 0;
}

uint32_t BlockDoubleTriDiagonalMatrixToCompactMatrixWithDeflation(BlockDoubleTriDiagonalMatrixType const * A, CompactMatrixType * const B)
{
  return BlockDoubleTriDiagonalMatrixToCompactMatrixDeflate(A, B, true);
}

uint32_t BlockDoubleTriDiagonalMatrixToCompactMatrix(BlockDoubleTriDiagonalMatrixType const * A, CompactMatrixType * const B)
{
  return BlockDoubleTriDiagonalMatrixToCompactMatrixDeflate(A, B, false);
}

/* on input A is a Hermitian matrix, Singular Value is the smallest number larger than 0 that has numerical relavance, BlockSize is the size of the blocks and NBlockMax the maximal number of Blocks
 *  onoutput Matrix B has NBlockMax blocks (a_{0} to a_{NBlockMax+1}. b_{0} to b_{NBlockMax}
 */
uint32_t CompactMatrixBlockBandTridiagonalizeHermitianNoReorthogonalization(CompactMatrixType const *A, BlockTriDiagonalMatrixType *const B, double const SingularValue, uint32_t NBlockMax, uint32_t BlockSize)
{
  //TODO: test after changing compact matrix order
  int i,NSingular;
  CompactMatrixType V,Vm1,Vm2,alpha,beta,R;
  B->N=NBlockMax;
  B->BlockSize=BlockSize;
  B->QComplex=A->QComplex;
  InitBlockTridiagonalMatrix(B);
  
  Vm1.Ni=A->Ni;
  Vm1.Nj=BlockSize;
  Vm1.QComplex=A->QComplex;
  if(InitCompactIdentityMatrix(&Vm1)){printf("Error in InitCompactIdentityMatrix in CompactMatrixBlockBandTridiagonalizeHermitian\n");return 1;};
  CompactMatrixMultiply3('N', 'N', 'C', &Vm1, A, &Vm1, &alpha);
  if(A->QComplex)
  {
    memcpy(&(B->a[BlockSize*BlockSize]), alpha.V, BlockSize*BlockSize*sizeof(Complex16));
  }
  else
  {
    memcpy(&(B->aRe[BlockSize*BlockSize]), alpha.VRe, BlockSize*BlockSize*sizeof(double));
  }
  CompactMatrixMultiply('C', 'C', A, &Vm1, &R);
  CompactMatrixAddMultiply('C', 'C', &Vm1, &alpha, &R, -1.0, 0.0, 1.0, 0.0); //need to check
  
  FreeCompactMatrix(&alpha);
  
  if((NSingular = CompactMatrixQRDecompositionEigenHermitian(&R, SingularValue, &V, &beta)))
  {
    printf("%i Singular values found at block level %i, deflation needed but not implemented\n",NSingular,1);
  }
  
  FreeCompactMatrix(&R);
  if(A->QComplex)
  {
    memcpy(&(B->b[BlockSize*BlockSize]), beta.V, BlockSize*BlockSize*sizeof(Complex16));
  }
  else
  {
    memcpy(&(B->bRe[BlockSize*BlockSize]), beta.VRe, BlockSize*BlockSize*sizeof(double));
  }
  CompactMatrixMultiply3('N', 'N', 'C', &V, A, &V, &alpha);
  if(A->QComplex)
  {
    memcpy(&(B->a[2*BlockSize*BlockSize]), alpha.V, BlockSize*BlockSize*sizeof(Complex16));
  }
  else
  {
    memcpy(&(B->aRe[2*BlockSize*BlockSize]), alpha.VRe, BlockSize*BlockSize*sizeof(double));
  }
  for(i=2;i<B->N;i++)
  {
    Vm2=Vm1;
    Vm1=V;
    CompactMatrixMultiply('C', 'C', A, &Vm1, &R);
    CompactMatrixAddMultiply('C', 'C', &Vm2, &beta, &R, -1.0, 0.0, 1.0, 0.0); //need to check
    CompactMatrixAddMultiply('C', 'C', &V,   &alpha, &R, -1.0, 0.0, 1.0, 0.0); //need to check
    FreeCompactMatrix(&alpha);
    FreeCompactMatrix(&beta);
    if((NSingular=CompactMatrixQRDecompositionEigenHermitian(&R, SingularValue, &V, &beta)))
    {
      printf("%i Singular values found at block level %i, deflation needed but not implemented\n",NSingular,i);
    }
    FreeCompactMatrix(&R);
    if(A->QComplex)
    {
      memcpy(&(B->b[i*BlockSize*BlockSize]), beta.V, BlockSize*BlockSize*sizeof(Complex16));
    }
    else
    {
      memcpy(&(B->bRe[i*BlockSize*BlockSize]), beta.VRe, BlockSize*BlockSize*sizeof(double));
    }
    CompactMatrixMultiply3('N', 'N', 'C', &V, A, &V, &alpha);
    if(A->QComplex)
    {
      memcpy(&(B->a[(i+1)*BlockSize*BlockSize]), alpha.V, BlockSize*BlockSize*sizeof(Complex16));
    }
    else
    {
      memcpy(&(B->aRe[(i+1)*BlockSize*BlockSize]), alpha.VRe, BlockSize*BlockSize*sizeof(double));
    }
    FreeCompactMatrix(&Vm2);
  }
  FreeCompactMatrix(&beta);
  FreeCompactMatrix(&alpha);
  FreeCompactMatrix(&Vm1);
  FreeCompactMatrix(&V);
  return 0;
}

/*  Returns a block tridiagonal form of matrix A in matrix B using the starting vector vec.
 *  B is not initialized on input
 *  vec is destroyed during the calculation. (values and size changed but not freed)
 *  if ReturnT is true vec contains T
 */
uint32_t CompactMatrixBlockBandDiagonalizeReal(CompactMatrixType const * const A, CompactMatrixType * vec, BlockTriDiagonalMatrixType *const B, double const SingularValue, uint32_t const NTri, uint32_t const NOrtho, bool const ReOrthogonalize, bool const ReturnT)
{
  uint32_t i,j,k,Ni;
  char const *charN="N";
  char const *charT="T";
  double one=1;
  double minusone=-1;
  double zero=0;
  uint32_t NSingular;
  double *b,*a,work[vec->Ni*vec->Ni];
  CompactMatrixType *veclist;
  B->BlockSize=vec->Ni;
  B->N=Min(A->Ni,NTri);
  B->QComplex=0;
  InitBlockTridiagonalMatrix(B);
  for(i=0;i<B->BlockSize;i++)
  {
    B->bRe[i+i*B->BlockSize]=0;
  }
  B->mu=A->mu;
  strncpy(B->Name, A->Name, MAXStrLength);
  b=B->bRe;
  a=B->aRe;
  veclist= (CompactMatrixType*)malloc((B->N+2)*sizeof(CompactMatrixType));
  veclist[0]=*vec;
  CompactMatrixCopy(veclist[0], &veclist[1]);
  CompactMatrixSVDOrthogonalizeRow(&veclist[1], SingularValue, true);
  NSingular=veclist[0].Ni-veclist[1].Ni;
  if(NSingular!=0)
  {
    VerbosityPrintf(0x0001,"Reduce dimension from %u to %u in block %u\n",veclist[0].Ni,veclist[1].Ni,0);
    for(k=1;k<=B->N;k++) B->BlockSizeList[k] =veclist[1].Ni;
  }
  if(ReOrthogonalize==true)
  {
    CompactMatrixGramSchmidtOrthogonalizeRow(&veclist[1], SingularValue, false);
  }
  dgemm(charT,charN,&veclist[1].Ni,&veclist[0].Ni,&veclist[0].Nj,&one,veclist[1].VRe,&veclist[1].Nj,veclist[0].VRe,&veclist[0].Nj,&zero,b,&B->BlockSize); // b = v[0] . v[1]^T
  for(i=1;i<=NTri;i++)
  {
    veclist[i+1].Ni=veclist[i].Ni;
    veclist[i+1].Nj=veclist[i].Nj;
    veclist[i+1].QComplex=veclist[i].QComplex;
    if(InitCompactMatrix(&veclist[i+1]))
    {
      // out off memory - could implement freeing some vectors and retry
      printf("InitCompactMatrix failed in CompactMatrixBlockBandDiagonalizeReal\n");
      return 1;
    }
    dgemm(charN,charN,&A->Ni,&veclist[i].Ni,&A->Nj,&one,A->VRe,&A->Nj,veclist[i].VRe,&veclist[i].Nj,&zero,veclist[i+1].VRe,&veclist[i+1].Nj); // v[i+1] = v[i] . A
    a+=B->BlockSize*B->BlockSize;
    dgemm(charT,charN,&veclist[i+1].Ni,&veclist[i].Ni,&veclist[i].Nj,&one,veclist[i+1].VRe,&veclist[i+1].Nj,veclist[i].VRe,&veclist[i].Nj,&zero,a,&B->BlockSize); // a = v[i] . v[i+1]^T
    if(i!=1)
    {
      b+=B->BlockSize*B->BlockSize;
      dgemm(charT,charN,&veclist[i+1].Ni,&veclist[i-1].Ni,&veclist[i-1].Nj,&one,veclist[i+1].VRe,&veclist[i+1].Nj,veclist[i-1].VRe,&veclist[i-1].Nj,&zero,b,&B->BlockSize); // b = v[i-1] . v[i+1]^T
    }
    if(i==NTri) break;
    NSingular=1;
    while(NSingular!=0)
    {
      // orthogonalize veclist[i+1] to veclist[i] ... veclist[Max(1,i+1-NOrtho)]
      // we use matrix multiplication of veclist[i+1] and veclist[k], which is a numerical unstable form of Gramm Schmidt
      // as a soltion we repeat the orthonormalization once more when finished.
      for(j=i;j>((NOrtho>i)?0:i-NOrtho);j--)
      {
        dgemm(charT,charN,&veclist[j].Ni,&veclist[i+1].Ni,&veclist[i+1].Nj,&one,veclist[j].VRe,&veclist[j].Nj,veclist[i+1].VRe,&veclist[i+1].Nj,&zero,work,&B->BlockSize); // w = v[i+1] . v[j]^T
        dgemm(charN,charN,&veclist[j].Nj,&veclist[i+1].Ni,&veclist[j].Ni,&minusone,veclist[j].VRe,&veclist[j].Nj,work,&B->BlockSize,&one,veclist[i+1].VRe,&veclist[i+1].Nj); // v[i+1] <= -1 * w . v[j] + 1 * v[i+1]
        dgemm(charT,charN,&veclist[j].Ni,&veclist[i+1].Ni,&veclist[i+1].Nj,&one,veclist[j].VRe,&veclist[j].Nj,veclist[i+1].VRe,&veclist[i+1].Nj,&zero,work,&B->BlockSize); // w = v[i+1] . v[j]^T
        dgemm(charN,charN,&veclist[j].Nj,&veclist[i+1].Ni,&veclist[j].Ni,&minusone,veclist[j].VRe,&veclist[j].Nj,work,&B->BlockSize,&one,veclist[i+1].VRe,&veclist[i+1].Nj); // v[i+1] <= -1 * w . v[j] + 1 * v[i+1]
      }
      if(ReOrthogonalize==true)
      {
        for(j=1;j<=i;j++)
        {
          dgemm(charT,charN,&veclist[j].Ni,&veclist[i+1].Ni,&veclist[i+1].Nj,&one,veclist[j].VRe,&veclist[j].Nj,veclist[i+1].VRe,&veclist[i+1].Nj,&zero,work,&B->BlockSize); // w = v[i+1] . v[j]^T
          dgemm(charN,charN,&veclist[j].Nj,&veclist[i+1].Ni,&veclist[j].Ni,&minusone,veclist[j].VRe,&veclist[j].Nj,work,&B->BlockSize,&one,veclist[i+1].VRe,&veclist[i+1].Nj); // v[i+1] <= -1 * w . v[j] + 1 * v[i+1]
          dgemm(charT,charN,&veclist[j].Ni,&veclist[i+1].Ni,&veclist[i+1].Nj,&one,veclist[j].VRe,&veclist[j].Nj,veclist[i+1].VRe,&veclist[i+1].Nj,&zero,work,&B->BlockSize); // w = v[i+1] . v[j]^T
          dgemm(charN,charN,&veclist[j].Nj,&veclist[i+1].Ni,&veclist[j].Ni,&minusone,veclist[j].VRe,&veclist[j].Nj,work,&B->BlockSize,&one,veclist[i+1].VRe,&veclist[i+1].Nj); // v[i+1] <= -1 * w . v[j] + 1 * v[i+1]
        }
      }
      NSingular=veclist[i+1].Ni;
      CompactMatrixSVDOrthogonalizeRow(&veclist[i+1], SingularValue, true);
      NSingular-=veclist[i+1].Ni;
      if(NSingular!=0)
      {
        VerbosityPrintf(0x0001,"Reduce dimension from %u to %u in block %u\n",veclist[i].Ni,veclist[i+1].Ni,i);
        for(k=i+1;k<=B->N;k++) B->BlockSizeList[k] =veclist[i+1].Ni;
        if(veclist[i+1].Ni==0) break;
      }
    }
    if(veclist[i+1].Ni==0)
    {
      B->N=i;
      break;
    }
    if(ReOrthogonalize==true)
    {
      CompactMatrixGramSchmidtOrthogonalizeRow(&veclist[i+1], SingularValue, false);
    }
  }
  if(ReturnT==true)
  {
    Ni=0;
    for(j=1;j<=Min(i+1,NTri);j++)
    {
      Ni=Ni+veclist[j].Ni;
    }
    ResizeCompactMatrix(vec, Ni, vec->Nj);
    Ni=0;
    for(j=1;j<=Min(i+1,NTri);j++)
    {
      memcpy(&(vec->VRe[Ni*vec->Nj]),veclist[j].VRe,veclist[j].Ni*veclist[j].Nj*sizeof(double));
      Ni=Ni+veclist[j].Ni;
    }
  }
  for(j=1;j<=Min(i+1,NTri+1);j++)
  {
    FreeCompactMatrix(&veclist[j]);
  }
  free(veclist);
  return 0;
}
uint32_t CompactMatrixBlockBandDiagonalizeComplex(CompactMatrixType const * const A, CompactMatrixType * vec, BlockTriDiagonalMatrixType *const B, double const SingularValue, uint32_t const NTri, uint32_t const NOrtho, bool const ReOrthogonalize,bool const ReturnT)
{
  uint32_t i,j,k,Ni;
  char const *charN="N";
  char const *charT="T";
  char const *charC="C";
  Complex16 one={1,0};
  Complex16 minusone={-1,0};
  Complex16 zero={0,0};
  uint32_t NSingular;
  Complex16 *b,*a,work[vec->Ni*vec->Ni];
  CompactMatrixType *veclist;
  B->BlockSize=vec->Ni;
  B->N=Min(A->Ni,NTri);
  B->QComplex=1;
  InitBlockTridiagonalMatrix(B);
  for(i=0;i<B->BlockSize;i++)
  {
    B->b[i+i*B->BlockSize].Re=0;
  }
  B->mu=A->mu;
  strncpy(B->Name, A->Name, MAXStrLength);
  b=B->b;
  a=B->a;
  veclist= (CompactMatrixType*)malloc((B->N+2)*sizeof(CompactMatrixType));
  veclist[0]=*vec;
  CompactMatrixCopy(veclist[0], &veclist[1]);
  CompactMatrixSVDOrthogonalizeRow(&veclist[1], SingularValue, true);
  NSingular=veclist[0].Ni-veclist[1].Ni;
  if(NSingular!=0)
  {
    VerbosityPrintf(0x0001,"Reduce dimension from %u to %u in block %u\n",veclist[0].Ni,veclist[1].Ni,0);
    for(k=1;k<=B->N;k++) B->BlockSizeList[k] = veclist[1].Ni;
  }
  if(ReOrthogonalize==true)
  {
    CompactMatrixGramSchmidtOrthogonalizeRow(&veclist[1], SingularValue, false);
  }
  zgemm(charC,charN,&veclist[1].Ni,&veclist[0].Ni,&veclist[0].Nj,&one,veclist[1].V,&veclist[1].Nj,veclist[0].V,&veclist[0].Nj,&zero,b,&B->BlockSize); // b = v[0] . v[1]^+
  for(i=1;i<=NTri;i++)
  {
    veclist[i+1].Ni=veclist[i].Ni;
    veclist[i+1].Nj=veclist[i].Nj;
    veclist[i+1].QComplex=veclist[i].QComplex;
    if(InitCompactMatrix(&veclist[i+1]))
    {
      // out off memory - could implement freeing some vectors and retry
      printf("InitCompactMatrix failed in CompactMatrixBlockBandDiagonalizeReal\n");
      return 1;
    }
    zgemm(charT,charN,&A->Nj,&veclist[i].Ni,&A->Ni,&one,A->V,&A->Ni,veclist[i].V,&veclist[i].Nj,&zero,veclist[i+1].V,&veclist[i+1].Nj); // v[i+1] = v[i] . A^* or v[i+1] = v[i] . A^T
    a+=B->BlockSize*B->BlockSize;
    zgemm(charC,charN,&veclist[i+1].Ni,&veclist[i].Ni,&veclist[i].Nj,&one,veclist[i+1].V,&veclist[i+1].Nj,veclist[i].V,&veclist[i].Nj,&zero,a,&B->BlockSize); // a = v[i]^* . v[i+1]^T -> a^* = v[i] . v[i+1]^+
    for(j=0;j<B->BlockSize*B->BlockSize;j++)
    {
      a[j].Im=-a[j].Im;
    }
    if(i!=1)
    {
      b+=B->BlockSize*B->BlockSize;
      zgemm(charC,charN,&veclist[i+1].Ni,&veclist[i-1].Ni,&veclist[i-1].Nj,&one,veclist[i+1].V,&veclist[i+1].Nj,veclist[i-1].V,&veclist[i-1].Nj,&zero,b,&B->BlockSize); // b = v[i-1]^* . v[i+1]^T
      for(j=0;j<B->BlockSize*B->BlockSize;j++)
      {
        b[j].Im=-b[j].Im;
      }
    }
    if(i==NTri) break;
    NSingular=1;
    while(NSingular!=0)
    {
      // orthogonalize veclist[i+1] to veclist[i] ... veclist[Max(1,i+1-NOrtho)]
      // we use matrix multiplication of veclist[i+1] and veclist[k], which is a numerical unstable form of Gramm Schmidt
      // as a soltion we repeat the orthonormalization once more when finished.
      for(j=i;j>((NOrtho>i)?0:i-NOrtho);j--)
      {
        zgemm(charC,charN,&veclist[j].Ni,&veclist[i+1].Ni,&veclist[i+1].Nj,&one,veclist[j].V,&veclist[j].Nj,veclist[i+1].V,&veclist[i+1].Nj,&zero,work,&B->BlockSize); // w = v[i+1] . v[j]^+
        zgemm(charN,charN,&veclist[j].Nj,&veclist[i+1].Ni,&veclist[j].Ni,&minusone,veclist[j].V,&veclist[j].Nj,work,&B->BlockSize,&one,veclist[i+1].V,&veclist[i+1].Nj); // v[i+1] <= -1 * w . v[j] + 1 * v[i+1]
        zgemm(charC,charN,&veclist[j].Ni,&veclist[i+1].Ni,&veclist[i+1].Nj,&one,veclist[j].V,&veclist[j].Nj,veclist[i+1].V,&veclist[i+1].Nj,&zero,work,&B->BlockSize); // w = v[i+1] . v[j]^+
        zgemm(charN,charN,&veclist[j].Nj,&veclist[i+1].Ni,&veclist[j].Ni,&minusone,veclist[j].V,&veclist[j].Nj,work,&B->BlockSize,&one,veclist[i+1].V,&veclist[i+1].Nj); // v[i+1] <= -1 * w . v[j] + 1 * v[i+1]
      }
      if(ReOrthogonalize==true)
      {
        for(j=1;j<=i;j++)
        {
          zgemm(charC,charN,&veclist[j].Ni,&veclist[i+1].Ni,&veclist[i+1].Nj,&one,veclist[j].V,&veclist[j].Nj,veclist[i+1].V,&veclist[i+1].Nj,&zero,work,&B->BlockSize); // w = v[i+1] . v[j]^+
          zgemm(charN,charN,&veclist[j].Nj,&veclist[i+1].Ni,&veclist[j].Ni,&minusone,veclist[j].V,&veclist[j].Nj,work,&B->BlockSize,&one,veclist[i+1].V,&veclist[i+1].Nj); // v[i+1] <= -1 * w . v[j] + 1 * v[i+1]
          zgemm(charC,charN,&veclist[j].Ni,&veclist[i+1].Ni,&veclist[i+1].Nj,&one,veclist[j].V,&veclist[j].Nj,veclist[i+1].V,&veclist[i+1].Nj,&zero,work,&B->BlockSize); // w = v[i+1] . v[j]^+
          zgemm(charN,charN,&veclist[j].Nj,&veclist[i+1].Ni,&veclist[j].Ni,&minusone,veclist[j].V,&veclist[j].Nj,work,&B->BlockSize,&one,veclist[i+1].V,&veclist[i+1].Nj); // v[i+1] <= -1 * w . v[j] + 1 * v[i+1]
        }
      }
      NSingular=veclist[i+1].Ni;
      CompactMatrixSVDOrthogonalizeRow(&veclist[i+1], SingularValue, true);
      NSingular-=veclist[i+1].Ni;
      if(NSingular!=0)
      {
        VerbosityPrintf(0x0001,"Reduce dimension from %u to %u in block %u\n",veclist[i].Ni,veclist[i+1].Ni,i);
        for(k=i+1;k<=B->N;k++) B->BlockSizeList[k] =veclist[i+1].Ni;
        if(veclist[i+1].Ni==0) break;
      }
    }
    if(veclist[i+1].Ni==0)
    {
      B->N = i;
      break;
    }
    if(ReOrthogonalize==true)
    {
      CompactMatrixGramSchmidtOrthogonalizeRow(&veclist[i+1], SingularValue, false);
    }
  }
  if(ReturnT==true)
  {
    Ni=0;
    for(j=1;j<=Min(i+1,NTri);j++)
    {
      Ni=Ni+veclist[j].Ni;
    }
    ResizeCompactMatrix(vec, Ni, vec->Nj);
    Ni=0;
    for(j=1;j<=Min(i+1,NTri);j++)
    {
      memcpy(&(vec->V[Ni*vec->Nj]),veclist[j].V,veclist[j].Ni*veclist[j].Nj*sizeof(Complex16));
      Ni=Ni+veclist[j].Ni;
    }
  }
  for(j=1;j<=Min(i+1,NTri+1);j++)
  {
    FreeCompactMatrix(&veclist[j]);
  }
  free(veclist);
  return 0;
}
uint32_t CompactMatrixBlockBandDiagonalize(CompactMatrixType const * const A, CompactMatrixType * vec, BlockTriDiagonalMatrixType *const B, double const SingularValue, uint32_t const NTri, uint32_t const NOrtho, bool const ReOrthogonalize, bool const ReturnT)
{
  uint32_t ret;
  CompactMatrixType Acopy;
  if(A->Ni==0)
  {
    B->BlockSize=vec->Ni;
    B->N=0;
    B->QComplex=0;
    InitBlockTridiagonalMatrix(B);
    B->mu=A->mu;
    strncpy(B->Name, A->Name, MAXStrLength);
    return 0;
  }
  if(A->QComplex==0)
  {
    if(vec->QComplex==0)
    {
      return CompactMatrixBlockBandDiagonalizeReal(A, vec, B, SingularValue, NTri, NOrtho, ReOrthogonalize, ReturnT);
    }
    else
    {
      CompactMatrixCopy(*A, &Acopy);
      MakeCompactRealMatrixComplex(&Acopy);
      ret = CompactMatrixBlockBandDiagonalizeComplex(&Acopy, vec, B, SingularValue, NTri, NOrtho, ReOrthogonalize, ReturnT);
      FreeCompactMatrix(&Acopy);
      return ret;
    }
  }
  else
  {
    if(vec->QComplex==0)
    {
      MakeCompactRealMatrixComplex(vec);
      return CompactMatrixBlockBandDiagonalizeComplex(A, vec, B, SingularValue, NTri, NOrtho, ReOrthogonalize, ReturnT);
    }
    else
    {
      return CompactMatrixBlockBandDiagonalizeComplex(A, vec, B, SingularValue, NTri, NOrtho, ReOrthogonalize, ReturnT);
    }
  }
}

/* returns 1 if O=O^dag, zero otherwise
 * epsilon gives the accuracy
 * test if | O[i][j] - O[j][i] | > epsilon * | O[i][j] + O[j][i] |
 */

uint32_t QCompactMatrixHermitian(const CompactMatrixType O, const double epsilon)
{
  uint32_t i,j;
  if(O.Ni!=O.Nj)return 0;
  if(O.QComplex)
  {
    for(i=0;i<O.Ni;i++)
    {
      for(j=0;j<=i;j++)
      {
        if(fabs(O.V[i*O.Nj+j].Re-O.V[j*O.Nj+i].Re)>(fabs(O.V[i*O.Nj+j].Re+O.V[j*O.Nj+i].Re) + fabs(O.V[i*O.Nj+j].Im-O.V[j*O.Nj+i].Im))*epsilon+DBL_MIN*100) return 0;
        if(fabs(O.V[i*O.Nj+j].Im+O.V[j*O.Nj+i].Im)>(fabs(O.V[i*O.Nj+j].Re+O.V[j*O.Nj+i].Re) + fabs(O.V[i*O.Nj+j].Im-O.V[j*O.Nj+i].Im))*epsilon+DBL_MIN*100) return 0;
      }
    }
  }
  else
  {
    for(i=0;i<O.Ni;i++)
    {
      for(j=0;j<i;j++)
      {
        if(fabs(O.VRe[i*O.Nj+j]-O.VRe[j*O.Nj+i])>((fabs(O.VRe[i*O.Nj+j]+O.VRe[j*O.Nj+i])*epsilon+DBL_MIN*100))) return 0;
      }
    }
  }
  return 1;
}


CompactMatrixType CompactzeroMatrix(const uint32_t Ni, const uint32_t Nj)
{
  CompactMatrixType zeroMat;
  zeroMat.Ni = Ni;
  zeroMat.Nj = Nj;
  zeroMat.QComplex=0;
  InitCompactMatrix(&zeroMat);
  return zeroMat;
}

CompactMatrixType CompactnullMatrix()
{
  return CompactzeroMatrix(1,0);
}

bool isNullMatrix(const CompactMatrixType M)
{
  return (M.Ni == 0 || M.Nj == 0);
}

CompactMatrixType MatrixFlatten(vector<vector<CompactMatrixType> > allMs, bool printWarnings)
{
  uint32_t QComplex = 0;
  CompactMatrixType fullM;
  uint32_t numOfRowsOuter = (uint32_t)allMs.size();
  uint32_t numOfColsOuter = (uint32_t)allMs[0].size();
  uint32_t numOfRowsInner[numOfRowsOuter], numOfColsInner[numOfColsOuter], numOfAllRows = 0, numOfAllCols = 0;
  for (int i = 0; i < numOfRowsOuter; ++i)
  {
    numOfRowsInner[i] = 0;
    for (int j = 0; j < numOfColsOuter; ++j)
    {
      if (!isNullMatrix(allMs[i][j]))
      {
        numOfRowsInner[i] = allMs[i][j].Ni;
        break;
      }
    }
    if (printWarnings && numOfRowsInner[i] == 0) printf("\nWarning: Row %d of matrix is empty and will be ignored.\n",i+1);
    numOfAllRows += numOfRowsInner[i];
  }
  for (int j = 0; j < numOfColsOuter; ++j)
  {
    numOfColsInner[j] = 0;
    for (int i = 0; i < numOfRowsOuter; ++i)
    {
      if (!isNullMatrix(allMs[i][j]))
      {
        numOfColsInner[j] = allMs[i][j].Nj;
        break;
      }
    }
    if (printWarnings && numOfColsInner[j] == 0) printf("\nWarning: Column %d of matrix is empty and will be ignored.\n",j+1);
    numOfAllCols += numOfColsInner[j];
  }
  for (int i = 0; i < numOfRowsOuter; ++i)
  {
    for (int j = 0; j < numOfColsOuter; ++j)
    {
      if (isNullMatrix(allMs[i][j]))
      {
      	FreeCompactMatrix(&(allMs[i][j]));
        allMs[i][j] = CompactzeroMatrix(numOfRowsInner[i],numOfColsInner[j]);
      }
      if (allMs[i][j].Nj != numOfColsInner[j])
      {
        printf("Matrix at position %d,%d has a different number of columns than it should have.", i + 1, j + 1);
        return CompactnullMatrix();
      }
      if (allMs[i][j].Ni != numOfRowsInner[i])
      {
        printf("Matrix at position %d,%d has a different number of rows than it should have.", i + 1, j + 1);
        return CompactnullMatrix();
      }
      if (!QComplex && allMs[i][j].QComplex) QComplex = allMs[i][j].QComplex;
    }
  }
  fullM.QComplex = QComplex;
  fullM.Ni = numOfAllRows;
  fullM.Nj = numOfAllCols;
  InitCompactMatrix(&fullM);
  int rowOffset = 0, colOffset = 0;
  for (int i = 0; i < numOfRowsOuter; ++i)
  {
    colOffset = 0;
    for (int j = 0; j < numOfColsOuter; ++j)
    {
      for (int k = 0; k < allMs[i][j].Ni; ++k)
      {
        for (int l = 0; l < allMs[i][j].Nj; ++l)
        {
          if(fullM.QComplex==0)
          {
            fullM.VRe[(rowOffset + k)*fullM.Nj+(colOffset + l)] = allMs[i][j].VRe[k*allMs[i][j].Nj+l];
          }
          else
          {
            if(allMs[i][j].QComplex==0)
            {
              fullM.V[(rowOffset + k)*fullM.Nj+(colOffset + l)].Re = allMs[i][j].VRe[k*allMs[i][j].Nj+l];
              fullM.V[(rowOffset + k)*fullM.Nj+(colOffset + l)].Im = 0;
            }
            else
            {
              fullM.V[(rowOffset + k)*fullM.Nj+(colOffset + l)] = allMs[i][j].V[k*allMs[i][j].Nj+l];
            }
          }
        }
      }
      colOffset += allMs[i][j].Nj;
    }
    rowOffset += allMs[i][0].Ni;
  }
  return fullM;
}

CompactMatrixType MatrixFlatten(vector<vector<CompactMatrixType> > allMs)
{
  return MatrixFlatten(allMs, false);
}

CompactMatrixType MatrixFlatten(vector<vector<CompactMatrixType *> > allMs, bool printWarnings)
{
  uint32_t QComplex = 0;
  CompactMatrixType fullM;
  uint32_t numOfRowsOuter = (uint32_t)allMs.size();
  uint32_t numOfColsOuter = (uint32_t)allMs[0].size();
  uint32_t numOfRowsInner[numOfRowsOuter], numOfColsInner[numOfColsOuter], numOfAllRows = 0, numOfAllCols = 0;
  for (int i = 0; i < numOfRowsOuter; ++i)
  {
    numOfRowsInner[i] = 0;
    for (int j = 0; j < numOfColsOuter; ++j)
    {
      if (!isNullMatrix(*allMs[i][j]))
      {
        numOfRowsInner[i] = allMs[i][j]->Ni;
        break;
      }
    }
    if (printWarnings && numOfRowsInner[i] == 0) printf("\nWarning: Row %d of matrix is empty and will be ignored.\n",i+1);
    numOfAllRows += numOfRowsInner[i];
  }
  for (int j = 0; j < numOfColsOuter; ++j)
  {
    numOfColsInner[j] = 0;
    for (int i = 0; i < numOfRowsOuter; ++i)
    {
      if (!isNullMatrix(*allMs[i][j]))
      {
        numOfColsInner[j] = allMs[i][j]->Nj;
        break;
      }
    }
    if (printWarnings && numOfColsInner[j] == 0) printf("\nWarning: Column %d of matrix is empty and will be ignored.\n",j+1);
    numOfAllCols += numOfColsInner[j];
  }
  for (int i = 0; i < numOfRowsOuter; ++i)
  {
    for (int j = 0; j < numOfColsOuter; ++j)
    {
      if (isNullMatrix(*allMs[i][j]))
      {
        allMs[i][j] = (CompactMatrixType*)malloc(sizeof(CompactMatrixType));
        *allMs[i][j] = CompactzeroMatrix(numOfRowsInner[i],numOfColsInner[j]);
      }
      if (allMs[i][j]->Nj != numOfColsInner[j])
      {
        printf("Matrix at position %d,%d has a different number of columns than it should have.", i + 1, j + 1);
        return CompactnullMatrix();
      }
      if (allMs[i][j]->Ni != numOfRowsInner[i])
      {
        printf("Matrix at position %d,%d has a different number of rows than it should have.", i + 1, j + 1);
        return CompactnullMatrix();
      }
      if (!QComplex && allMs[i][j]->QComplex) QComplex = allMs[i][j]->QComplex;
    }
  }
  fullM.QComplex = QComplex;
  fullM.Ni = numOfAllRows;
  fullM.Nj = numOfAllCols;
  InitCompactMatrix(&fullM);
  int rowOffset = 0, colOffset = 0;
  for (int i = 0; i < numOfRowsOuter; ++i)
  {
    colOffset = 0;
    for (int j = 0; j < numOfColsOuter; ++j)
    {
      for (int k = 0; k < allMs[i][j]->Ni; ++k)
      {
        for (int l = 0; l < allMs[i][j]->Nj; ++l)
        {
          if(fullM.QComplex==0)
          {
            fullM.VRe[(rowOffset + k)*fullM.Nj+(colOffset + l)] = allMs[i][j]->VRe[k*allMs[i][j]->Nj+l];
          }
          else
          {
            if(allMs[i][j]->QComplex==0)
            {
              fullM.V[(rowOffset + k)*fullM.Nj+(colOffset + l)].Re = allMs[i][j]->VRe[k*allMs[i][j]->Nj+l];
              fullM.V[(rowOffset + k)*fullM.Nj+(colOffset + l)].Im = 0;
            }
            else
            {
              fullM.V[(rowOffset + k)*fullM.Nj+(colOffset + l)] = allMs[i][j]->V[k*allMs[i][j]->Nj+l];
            }
          }
        }
      }
      colOffset += allMs[i][j]->Nj;
    }
    rowOffset += allMs[i][0]->Ni;
  }
  return fullM;
}

CompactMatrixType MatrixFlatten(vector<vector<CompactMatrixType *> > allMs)
{
  return MatrixFlatten(allMs, false);
}


CompactMatrixType removeSpin(const CompactMatrixType &M)
{
  CompactMatrixType out;
  out.Ni = M.Ni / 2;
  out.Nj = M.Nj / 2;
  out.QComplex = M.QComplex;
  InitCompactMatrix(&out);
  if (out.QComplex == 0)
  {
    for (int i = 0; i < out.Ni; ++i)
    {
      for (int j = 0; j < out.Nj; ++j)
      {
        out.VRe[i * out.Nj + j] = M.VRe[2 * i * M.Nj + 2 * j];
      }
    }
  } else
  {
    for (int i = 0; i < out.Ni; ++i)
    {
      for (int j = 0; j < out.Nj; ++j)
      {
        out.V[i * out.Nj + j].Re = M.V[2 * i * M.Nj + 2 * j].Re;
        out.V[i * out.Nj + j].Im = M.V[2 * i * M.Nj + 2 * j].Im;
      }
    }
  }
  return out;
}

CompactMatrixType ExpandNumberToCompactMatrix(int size, double Re, double Im = 0)
{
  CompactMatrixType M;
  M.QComplex = (Im != 0);
  M.Ni = size;
  M.Nj = size;
  InitCompactMatrix(&M);
  if(M.QComplex==0)
  {
    for (int i = 0; i < size; ++i)
    {
      M.VRe[i*M.Nj+i] = Re;
    }
  }
  else
  {
    for (int i = 0; i < size; ++i)
    {
      M.V[i*M.Nj+i].Re = Re;
      M.V[i*M.Nj+i].Im = Im;
    }
  }
  return M;
}

CompactMatrixType MatrixEnlarge(CompactMatrixType M, int size)
{
	CompactMatrixType out;
	out.Ni = M.Ni * size;
	out.Nj = M.Nj * size;
	out.QComplex = M.QComplex;
	InitCompactMatrix(&out);
	for (uint32_t i = 0; i < M.Ni; ++i)
	{
		for (uint32_t j = 0; j < M.Nj; ++j)
		{
			for (uint32_t k = 0; k < size; ++k)
			{
				uint32_t iNew = i * size + k;
				uint32_t jNew = j * size + k;
				if ((bool) out.QComplex)
				{
					out.V[iNew * out.Nj + jNew].Re = M.V[i * M.Nj + j].Re;
					out.V[iNew * out.Nj + jNew].Im = M.V[i * M.Nj + j].Im;
				}
				else
				{
					out.VRe[iNew * out.Nj + jNew] = M.VRe[i * M.Nj + j];
				}
			}
		}
	}
	return out;
}

CompactMatrixType MatrixExpand(const CompactMatrixType M, const int size, const Complex16 newEntries)
{
  if (size < M.Ni || size < M.Nj)
  {
    cout << "Error in MatrixExpand: New size " << size << "x" << size << " of expanded matrix must be larger than or equal to old size " << M.Ni << "x" << M.Ni << ".\n";
    return CompactnullMatrix();
  }
  if (M.Ni != M.Nj)
  {
    cout << "Error in MatrixExpand: Input of nono-quadratic Matrix (size " << M.Ni << "x" << M.Ni << ".\n";
    return CompactnullMatrix();
  }
  vector < vector<CompactMatrixType> > matMat = CompactnullMatrixMatrix(2, 2);
  matMat[0][0] = M;
  CompactMatrixType identityMat = CompactidentityMatrix(size - M.Ni);
  CompactMatrixScaleComplex(&identityMat, newEntries);
  matMat[1][1] = identityMat;
  return MatrixFlatten(matMat);
}

CompactMatrixType MatrixExpand(const CompactMatrixType M, const int size)
{
  Complex16 z;
  z.Re = 1;
  z.Im = 0;
  return MatrixExpand(M, size, z);
}

void MatrixAbsAverage(CompactMatrixType M, double *AbsAverage)
{
  uint32_t i,j;
  *AbsAverage=0;
  if(M.QComplex==0)
  {
    for(i=0;i<M.Ni;i++)
    {
      for(j=0;j<M.Nj;j++)
      {
        *AbsAverage+=Abs(M.VRe[i*M.Nj+j]);
      }
    }
  }
  else
  {
    for(i=0;i<M.Ni;i++)
    {
      for(j=0;j<M.Nj;j++)
      {
        *AbsAverage+=Abs(M.V[i*M.Nj+j].Re)+Abs(M.V[i*M.Nj+j].Im);
      }
    }
  }
  *AbsAverage/=(M.Ni*M.Nj);
}

vector<vector<CompactMatrixType> > CompactnullMatrixMatrix(const uint32_t rows, const uint32_t cols)
{
  vector < vector<CompactMatrixType> > matMat(rows);
  for (int i = 0; i < rows; ++i)
  {
    vector < CompactMatrixType > col(cols);
    matMat[i] = col;
    for (int j = 0; j < cols; ++j)
      matMat[i][j] = CompactnullMatrix();
  }
  return matMat;
}

CompactMatrixType CompactidentityMatrix(const int rowSize, const int colSize)
{
  CompactMatrixType identityMat;
  identityMat.QComplex = false;
  identityMat.Ni = rowSize;
  identityMat.Nj = colSize;
  InitCompactIdentityMatrix(&identityMat);
  return identityMat;
}

CompactMatrixType CompactidentityMatrix(const int size)
{
  return CompactidentityMatrix(size, size);
}

CompactMatrixType Submatrix(const CompactMatrixType &M, const vector<uint32_t> &inds)
{
	for (uint32_t i = 0; i < (uint32_t) inds.size(); ++i)
	{
		if (inds[i] < 0 || inds[i] >= M.Ni || inds[i] >= M.Nj)
		{
			cout << "Error in Submatrix: index " << inds[i] << " out of range.";
			return CompactnullMatrix();
		}
	}
	CompactMatrixType Msub;
	Msub.Ni = (uint32_t) inds.size();
	Msub.Nj = (uint32_t) inds.size();
	Msub.QComplex = M.QComplex;
	InitCompactMatrix(&Msub);
	Msub.OnLuaStack = M.OnLuaStack;
	for (int i = 0; i < (uint32_t) inds.size(); ++i)
	{
		for (int j = 0; j < (uint32_t) inds.size(); ++j)
		{
			uint32_t ind = inds[i] * M.Nj + inds[j];
			if (M.QComplex == 0)
			{
				Msub.VRe[i * Msub.Nj + j] = M.VRe[ind];
			}
			else
			{
				Msub.V[i * Msub.Nj + j] = M.V[ind];
			}
		}
	}
	return Msub;
}

CompactMatrixType Submatrix(const CompactMatrixType &M, int rowStart, int rows, int colStart, int cols)
{
  if (rowStart < 0 || rows <= 0 || colStart < 0 || cols <= 0 || rowStart + rows > M.Ni || colStart + cols > M.Nj)
  {
    cout << "Error in Submatrix: Invalid parameters:" << endl;
    cout << "rowStart = " << rowStart << ", rows = " << rows << ", colStart = " << colStart << " and cols = " << cols << " is incompatible with a matrix of size " << M.Ni << "x" << M.Nj << "." << endl;
    return CompactnullMatrix();
  }
  CompactMatrixType Msub;
  Msub.Ni = rows;
  Msub.Nj = cols;
  Msub.QComplex = M.QComplex;
  InitCompactMatrix(&Msub);
  Msub.OnLuaStack = M.OnLuaStack;
  for (int i = 0; i < rows; ++i)
  {
    for (int j = 0; j < cols; ++j)
    {
      if(M.QComplex==0)
      {
        Msub.VRe[i*Msub.Nj+j] = M.VRe[(rowStart + i)*M.Nj+(colStart + j)];
      }
      else
      {
        Msub.V[i*Msub.Nj+j] = M.V[(rowStart + i)*M.Nj+(colStart + j)];
      }
    }
  }
  return Msub;
}

CompactMatrixType Submatrix(const CompactMatrixType &M, int rows, int cols)
{
  return Submatrix(M, 0, rows, 0, cols);
}

CompactMatrixType Submatrix(const CompactMatrixType &M, int size)
{
  return Submatrix(M, size, size);
}

OperatorType CompactMatrixToOperator(const CompactMatrixType &M)
{
 OperatorType op;
 InitOperator(&op);
 InitOperatorFixedLength(&op, 2, 8);
 op.QComplex = M.QComplex;
 OperatorChangeNFandNB(&op, M.Ni, 0);
 uint16_t oadd[2];
 if ((bool) op.QComplex)
 {
   for (uint16_t i = 0; i < (uint16_t) M.Ni; ++i)
   {
     for (uint16_t j = 0; j < (uint16_t) M.Nj; ++j)
     {
       double VRe = M.V[i * M.Nj + j].Re;
       double VIm = M.V[i * M.Nj + j].Im;
       if (fabs(VRe) > 0 || fabs(VIm) > 0)
       {
         oadd[1] = i + CrAnMask;
         oadd[0] = j;
         OperatorAddLader(&op, oadd, 2, VRe, VIm);
       }
     }
   }
 }
 else
 {
   for (uint16_t i = 0; i < (uint16_t) M.Ni; ++i)
   {
     for (uint16_t j = 0; j < (uint16_t) M.Nj; ++j)
     {
       double VRe = M.VRe[i * M.Nj + j];
       if (fabs(VRe) > 0)
       {
         oadd[1] = i + CrAnMask;
         oadd[0] = j;
         OperatorAddLader(&op, oadd, 2, VRe, 0.);
       }
     }
   }
 }
 return op;
}
