/* Free a matrix M itselve (i.e. the structure is not freed)
 */
void FreeCompactMatrix(CompactMatrixType * const M)
{
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

/* Initiallize a matrix M must already be allocated
 *  M.QComplex, M.Ni, M,Nj must be set
 */
uint32_t InitCompactMatrix(CompactMatrixType * const M)
{
  M->Name=(char*)malloc(MAXStrLength*sizeof(char));
  if(M->Name==NULL){printf("malloc for name failed in InitCompactMatrix\n");return 1;};
  strcpy(M->Name, "Matrix");
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
    if(M->V==NULL)
    {
      printf("Malloc failed in InitCompactMatrix for real allocation\n");
      fflush(stdout);
      return 1;
    }
  }
  return 0;
}

/*
 *
 */
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

uint32_t CompactMatrixCopy(CompactMatrixType const A, CompactMatrixType * const B)
{
  B->Ni=A.Ni;
  B->Nj=A.Nj;
  B->QComplex=A.QComplex;
  B->a0=A.a0;
  B->b0=A.b0;
  B->mu=A.mu;
  if(InitCompactMatrix(B))
  {
    printf("InitCompactMatrix failed in CompactMatrixCopy\n");
    return 1;
  }
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


/* Initiallize a matrix M must already be allocated
 *  M.QComplex, M.Ni, M,Nj must be set
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

/*   y = A x
 *   TransA can be either:
 *          N for nothing
 *          T for transposed
 *          C for conjugate transposed
 */
uint32_t CompactMatrixVector(char const TransA, CompactMatrixType const *A, void *vx, void *vy)
{
  int INC=1;
  if(A->QComplex)
  {
    Complex16 ALPHA,BETA;
    ALPHA.Re=1; ALPHA.Im=0;
    BETA.Re=0;  BETA.Im=0;
    zgemv(&TransA,&(A->Nj),&(A->Ni),&ALPHA,A->V,&(A->Nj),(Complex16*)vx,&INC,&BETA,(Complex16*)vy,&INC);
  }
  else
  {
    double ALPHA=1,BETA=0;
    dgemv(&TransA,&(A->Nj),&(A->Ni),&ALPHA,A->VRe,&(A->Nj),(double*)vx,&INC,&BETA,(double*)vy,&INC);
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

/*  Calculates T^dag A T
 *  dimensions of input T : T.Ni by T.Nj
 *                      A : T.Nj by T.Nj
 *  dimensions of output: T.Ni by T.Ni
 */
uint32_t TransformDiagonalMatrixByCompactMatrixTdagAT(double const *A, CompactMatrixType const *T, CompactMatrixType * const TdagAT)
{
  uint32_t i,j,k,TATindex,icol,jcol;
  TdagAT->Ni=T->Ni;
  TdagAT->Nj=T->Ni;
  TdagAT->QComplex = T->QComplex;
  InitCompactMatrix(TdagAT);
  TATindex=0;
  if(T->QComplex)
  {
    for(i=0;i<TdagAT->Ni;i++)
    {
      for(j=0;j<TdagAT->Nj;j++)
      {
        icol=i*T->Nj;
        jcol=j*T->Nj;
        TdagAT->V[TATindex].Re=0;
        TdagAT->V[TATindex].Im=0;
        for(k=0;k<T->Nj;k++)
        {
          TdagAT->V[TATindex].Re+=T->V[jcol].Re*A[k]*T->V[icol].Re + T->V[jcol].Im*A[k]*T->V[icol].Im;
          TdagAT->V[TATindex].Im+=T->V[jcol].Re*A[k]*T->V[icol].Im - T->V[jcol].Im*A[k]*T->V[icol].Re;
          icol++;
          jcol++;
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
        icol=i*T->Nj;
        jcol=j*T->Nj;
        TdagAT->VRe[TATindex]=0;
        for(k=0;k<T->Nj;k++)
        {
          TdagAT->VRe[TATindex]+=T->VRe[jcol]*A[k]*T->VRe[icol];
          icol++;
          jcol++;
        }
        TATindex++;
      }
    }
  }
  return 0;
}


/*  Calculates T A T^dag
 *  dimensions of input T : T.Ni by T.Nj
 *                      A : T.Ni by T.Ni
 *  dimensions of output: T.Nj by T.Nj
 */
uint32_t TransformDiagonalMatrixByCompactMatrixTATdag(double const *A, CompactMatrixType const *T, CompactMatrixType * const TATdag)
{
  uint32_t i,j,k,TATindex,irow,jrow;
  TATdag->Ni=T->Nj;
  TATdag->Nj=T->Nj;
  TATdag->QComplex = T->QComplex;
  InitCompactMatrix(TATdag);
  TATindex=0;
  if(T->QComplex)
  {
    for(i=0;i<TATdag->Ni;i++)
    {
      for(j=0;j<TATdag->Nj;j++)
      {
        irow=i;
        jrow=j;
        TATdag->V[TATindex].Re=0;
        TATdag->V[TATindex].Im=0;
        for(k=0;k<T->Ni;k++)
        {
          TATdag->V[TATindex].Re+=T->V[jrow].Re*A[k]*T->V[irow].Re + T->V[jrow].Im*A[k]*T->V[irow].Im;
          TATdag->V[TATindex].Im+=T->V[jrow].Im*A[k]*T->V[irow].Re - T->V[jrow].Re*A[k]*T->V[irow].Im;
          irow+=T->Nj;
          jrow+=T->Nj;
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
        irow=i;
        jrow=j;
        TATdag->VRe[TATindex]=0;
        for(k=0;k<T->Nj;k++)
        {
          TATdag->VRe[TATindex]+=T->VRe[jrow]*A[k]*T->VRe[irow];
          irow+=T->Nj;
          jrow+=T->Nj;
        }
        TATindex++;
      }
    }
  }
  return 0;
}


/*   C = A ^{TransA} * B ^{TransB}
 *   C should not be allocated before the routine is called
 *   TransA can be either:
 *          N for nothing
 *          T for transposed
 *          C for conjugate transposed
 */
uint32_t CompactMatrixMultiply(char const TransA, char const TransB, CompactMatrixType const *A, CompactMatrixType const *B, CompactMatrixType *C)
{
  int Nk,LDA,LDB,LDC,i;
  Complex16 ALPHA,BETA,*temp;
  ALPHA.Re=1; ALPHA.Im=0;
  BETA.Re=0;  BETA.Im=0;
  if((TransA=='N')||(TransA=='n'))
  {
    if((TransB=='N')||(TransB=='n'))
    {
      C->Ni=B->Ni;
      C->Nj=A->Nj;
      Nk=B->Nj;
      if(B->Nj!=A->Ni)
      {
        printf("Trying to multiply matrix (%s)^%c with dimension %i %i with matrix (%s)^%c with dimension %i %i\n",A->Name,TransA,A->Ni,A->Nj,B->Name,TransB,B->Ni,B->Nj);
        printf("Error in dimensions\n");
        return 1;
      }
      LDA=A->Nj;
      LDB=B->Nj;
      LDC=C->Nj;
    }
    else
    {
      C->Ni=B->Nj;
      C->Nj=A->Nj;
      Nk=B->Ni;
      if(B->Ni!=A->Ni)
      {
        printf("Trying to multiply matrix (%s)^%c with dimension %i %i with matrix (%s)^%c with dimension %i %i\n",A->Name,TransA,A->Ni,A->Nj,B->Name,TransB,B->Ni,B->Nj);
        printf("Error in dimensions\n");
        return 1;
      }
      LDA=A->Nj;
      LDB=B->Nj;
      LDC=C->Nj;
    }
  }
  else
  {
    if((TransB=='N')||(TransB=='n'))
    {
      C->Ni=B->Ni;
      C->Nj=A->Ni;
      Nk=B->Nj;
      if(B->Nj!=A->Nj)
      {
        printf("Trying to multiply matrix (%s)^%c with dimension %i %i with matrix (%s)^%c with dimension %i %i\n",A->Name,TransA,A->Ni,A->Nj,B->Name,TransB,B->Ni,B->Nj);
        printf("Error in dimensions\n");
        return 1;
      }
      LDA=A->Nj;
      LDB=B->Nj;
      LDC=C->Nj;
    }
    else
    {
      C->Ni=B->Nj;
      C->Nj=A->Ni;
      Nk=B->Ni;
      if(B->Ni!=A->Nj)
      {
        printf("Trying to multiply matrix (%s)^%c with dimension %i %i with matrix (%s)^%c with dimension %i %i\n",A->Name,TransA,A->Ni,A->Nj,B->Name,TransB,B->Ni,B->Nj);
        printf("Error in dimensions\n");
        return 1;
      }
      LDA=A->Nj;
      LDB=B->Nj;
      LDC=C->Nj;
    }
  }
  C->QComplex=0;
  if((A->QComplex)||(B->QComplex)) C->QComplex=1;
  InitCompactMatrix(C);
  if(A->QComplex)
  {
    if(B->QComplex)
    {
      zgemm(&TransA,&TransB,&C->Nj,&C->Ni,&Nk,&(ALPHA),A->V,&LDA,B->V,&LDB,&(BETA),C->V,&LDC);
    }
    else
    {
      /* A is complex but B is real */
      temp=calloc(B->Ni*B->Nj,sizeof(Complex16));
      for(i=0;i<B->Ni*B->Nj;i++)
      {
        temp[i].Re=B->VRe[i];
      }
      zgemm(&TransA,&TransB,&C->Nj,&C->Ni,&Nk,&(ALPHA),A->V,&LDA,temp,&LDB,&(BETA),C->V,&LDC);
      free(temp);
    }
  }
  else
  {
    if(B->QComplex)
    {
      /* B is complex but A is real */
      temp=calloc(A->Ni*A->Nj,sizeof(Complex16));
      for(i=0;i<A->Ni*A->Nj;i++)
      {
        temp[i].Re=A->VRe[i];
      }
      zgemm(&TransA,&TransB,&C->Nj,&C->Ni,&Nk,&(ALPHA),temp,&LDA,B->V,&LDB,&(BETA),C->V,&LDC);
      free(temp);
    }
    else
    {
      dgemm(&TransA,&TransB,&C->Nj,&C->Ni,&Nk,&(ALPHA.Re),A->VRe,&LDA,B->VRe,&LDB,&(BETA.Re),C->VRe,&LDC);
    }
  }
  return 0;
}

/*   C = A ^{TransA} * B ^{TransB}
 *   C should be allocated before the routine is called
 *   TransA can be either:
 *          N for nothing
 *          T for transposed
 *          C for conjugate transposed
 */
uint32_t CompactMatrixMultiplyPreAllocated(char const TransA, char const TransB, CompactMatrixType const *A, CompactMatrixType const *B, CompactMatrixType *C)
{
  int Nk,LDA,LDB,LDC,i;
  Complex16 ALPHA,BETA,*temp,*tempA,*tempB;
  ALPHA.Re=1.0; ALPHA.Im=0.0;
  BETA.Re=0;  BETA.Im=0;
  if(((A->QComplex)||(B->QComplex))&&(C->QComplex==0))
  {
    printf("Trying to multiply a complex matrix into a real matrix\n");
    return 1;
  }
  if((TransA=='N')||(TransA=='n'))
  {
    if((TransB=='N')||(TransB=='n'))
    {
      Nk=B->Nj;
      if((B->Nj!=A->Ni)||(C->Ni!=B->Ni)||(C->Nj!=A->Nj))
      {
        printf("Trying to multiply matrix (%s)^%c with dimension %i %i with matrix (%s)^%c with dimension %i %i and add it to matrix (%s) with dimension %i %i\n",A->Name,TransA,A->Ni,A->Nj,B->Name,TransB,B->Ni,B->Nj,C->Name,C->Ni,C->Nj);
        printf("Error in dimensions\n");
        return 1;
      }
      LDA=A->Nj;
      LDB=B->Nj;
      LDC=C->Nj;
    }
    else
    {
      Nk=B->Ni;
      if((B->Ni!=A->Ni)||(C->Ni!=B->Nj)||(C->Nj!=A->Nj))
      {
        printf("Trying to multiply matrix (%s)^%c with dimension %i %i with matrix (%s)^%c with dimension %i %i and add it to matrix (%s) with dimension %i %i\n",A->Name,TransA,A->Ni,A->Nj,B->Name,TransB,B->Ni,B->Nj,C->Name,C->Ni,C->Nj);
        printf("Error in dimensions\n");
        return 1;
      }
      LDA=A->Nj;
      LDB=B->Nj;
      LDC=C->Nj;
    }
  }
  else
  {
    if((TransB=='N')||(TransB=='n'))
    {
      Nk=B->Nj;
      if((B->Nj!=A->Nj)||(C->Ni!=B->Ni)||(C->Nj!=A->Ni))
      {
        printf("Trying to multiply matrix (%s)^%c with dimension %i %i with matrix (%s)^%c with dimension %i %i and add it to matrix (%s) with dimension %i %i\n",A->Name,TransA,A->Ni,A->Nj,B->Name,TransB,B->Ni,B->Nj,C->Name,C->Ni,C->Nj);
        printf("Error in dimensions\n");
        return 1;
      }
      LDA=A->Nj;
      LDB=B->Nj;
      LDC=C->Nj;
    }
    else
    {
      Nk=B->Ni;
      if((B->Ni!=A->Nj)||(C->Ni!=B->Nj)||(C->Nj!=A->Ni))
      {
        printf("Trying to multiply matrix (%s)^%c with dimension %i %i with matrix (%s)^%c with dimension %i %i and add it to matrix (%s) with dimension %i %i\n",A->Name,TransA,A->Ni,A->Nj,B->Name,TransB,B->Ni,B->Nj,C->Name,C->Ni,C->Nj);
        printf("Error in dimensions\n");
        return 1;
      }
      LDA=A->Nj;
      LDB=B->Nj;
      LDC=C->Nj;
    }
  }
  if(A->QComplex)
  {
    if(B->QComplex)
    {
      zgemm(&TransA,&TransB,&C->Nj,&C->Ni,&Nk,&(ALPHA),A->V,&LDA,B->V,&LDB,&(BETA),C->V,&LDC);
    }
    else
    {
      /* A is complex but B is real */
      temp=calloc(B->Ni*B->Nj,sizeof(Complex16));
      for(i=0;i<B->Ni*B->Nj;i++)
      {
        temp[i].Re=B->VRe[i];
      }
      zgemm(&TransA,&TransB,&C->Nj,&C->Ni,&Nk,&(ALPHA),A->V,&LDA,temp,&LDB,&(BETA),C->V,&LDC);
      free(temp);
    }
  }
  else
  {
    if(B->QComplex)
    {
      /* B is complex but A is real */
      temp=calloc(A->Ni*A->Nj,sizeof(Complex16));
      for(i=0;i<A->Ni*A->Nj;i++)
      {
        temp[i].Re=A->VRe[i];
      }
      zgemm(&TransA,&TransB,&C->Nj,&C->Ni,&Nk,&(ALPHA),temp,&LDA,B->V,&LDB,&(BETA),C->V,&LDC);
      free(temp);
    }
    else
    {
      /* A and B real */
      if(C->QComplex)
      {
        tempA=calloc(A->Ni*A->Nj,sizeof(Complex16));
        for(i=0;i<A->Ni*A->Nj;i++)
        {
          tempA[i].Re=A->VRe[i];
        }
        tempB=calloc(B->Ni*B->Nj,sizeof(Complex16));
        for(i=0;i<B->Ni*B->Nj;i++)
        {
          tempB[i].Re=B->VRe[i];
        }
        zgemm(&TransA,&TransB,&C->Nj,&C->Ni,&Nk,&(ALPHA),tempA,&LDA,tempB,&LDB,&(BETA),C->V,&LDC);
        free(tempA);
        free(tempB);
      }
      else
      {
        dgemm(&TransA,&TransB,&C->Nj,&C->Ni,&Nk,&(ALPHA.Re),A->VRe,&LDA,B->VRe,&LDB,&(BETA.Re),C->VRe,&LDC);
      }
    }
  }
  return 0;
}

/*   C = C + alpha (A ^{TransA} * B ^{TransB})
 *   C should not be allocated before the routine is called
 *   TransA can be either:
 *          N for nothing
 *          T for transposed
 *          C for conjugate transposed
 */
uint32_t CompactMatrixAddMultiply(char const TransA, char const TransB, CompactMatrixType const *A, CompactMatrixType const *B, CompactMatrixType *C, double alphaRe, double alphaIm)
{
  int Nk,LDA,LDB,LDC,i;
  Complex16 ALPHA,BETA,*temp,*tempA,*tempB;
  ALPHA.Re=alphaRe; ALPHA.Im=alphaIm;
  BETA.Re=1;  BETA.Im=0;
  
  if(((A->QComplex)||(B->QComplex)||(alphaIm!=0))&&(C->QComplex==0))
  {
    MakeCompactRealMatrixComplex(C);
  }
  if((TransA=='N')||(TransA=='n'))
  {
    if((TransB=='N')||(TransB=='n'))
    {
      Nk=B->Nj;
      if((B->Nj!=A->Ni)||(C->Ni!=B->Ni)||(C->Nj!=A->Nj))
      {
        printf("Trying to multiply matrix (%s)^%c with dimension %i %i with matrix (%s)^%c with dimension %i %i and add it to matrix (%s) with dimension %i %i\n",A->Name,TransA,A->Ni,A->Nj,B->Name,TransB,B->Ni,B->Nj,C->Name,C->Ni,C->Nj);
        printf("Error in dimensions\n");
        return 1;
      }
      LDA=A->Nj;
      LDB=B->Nj;
      LDC=C->Nj;
    }
    else
    {
      Nk=B->Ni;
      if((B->Ni!=A->Ni)||(C->Ni!=B->Nj)||(C->Nj!=A->Nj))
      {
        printf("Trying to multiply matrix (%s)^%c with dimension %i %i with matrix (%s)^%c with dimension %i %i and add it to matrix (%s) with dimension %i %i\n",A->Name,TransA,A->Ni,A->Nj,B->Name,TransB,B->Ni,B->Nj,C->Name,C->Ni,C->Nj);
        printf("Error in dimensions\n");
        return 1;
      }
      LDA=A->Nj;
      LDB=B->Nj;
      LDC=C->Nj;
    }
  }
  else
  {
    if((TransB=='N')||(TransB=='n'))
    {
      Nk=B->Nj;
      if((B->Nj!=A->Nj)||(C->Ni!=B->Ni)||(C->Nj!=A->Ni))
      {
        printf("Trying to multiply matrix (%s)^%c with dimension %i %i with matrix (%s)^%c with dimension %i %i and add it to matrix (%s) with dimension %i %i\n",A->Name,TransA,A->Ni,A->Nj,B->Name,TransB,B->Ni,B->Nj,C->Name,C->Ni,C->Nj);
        printf("Error in dimensions\n");
        return 1;
      }
      LDA=A->Nj;
      LDB=B->Nj;
      LDC=C->Nj;
    }
    else
    {
      Nk=B->Ni;
      if((B->Ni!=A->Nj)||(C->Ni!=B->Nj)||(C->Nj!=A->Ni))
      {
        printf("Trying to multiply matrix (%s)^%c with dimension %i %i with matrix (%s)^%c with dimension %i %i and add it to matrix (%s) with dimension %i %i\n",A->Name,TransA,A->Ni,A->Nj,B->Name,TransB,B->Ni,B->Nj,C->Name,C->Ni,C->Nj);
        printf("Error in dimensions\n");
        return 1;
      }
      LDA=A->Nj;
      LDB=B->Nj;
      LDC=C->Nj;
    }
  }
  if(A->QComplex)
  {
    if(B->QComplex)
    {
      zgemm(&TransA,&TransB,&C->Nj,&C->Ni,&Nk,&(ALPHA),A->V,&LDA,B->V,&LDB,&(BETA),C->V,&LDC);
    }
    else
    {
      /* A is complex but B is real */
      temp=calloc(B->Ni*B->Nj,sizeof(Complex16));
      for(i=0;i<B->Ni*B->Nj;i++)
      {
        temp[i].Re=B->VRe[i];
      }
      zgemm(&TransA,&TransB,&C->Nj,&C->Ni,&Nk,&(ALPHA),A->V,&LDA,temp,&LDB,&(BETA),C->V,&LDC);
      free(temp);
    }
  }
  else
  {
    if(B->QComplex)
    {
      /* B is complex but A is real */
      temp=calloc(A->Ni*A->Nj,sizeof(Complex16));
      for(i=0;i<A->Ni*A->Nj;i++)
      {
        temp[i].Re=A->VRe[i];
      }
      zgemm(&TransA,&TransB,&C->Nj,&C->Ni,&Nk,&(ALPHA),temp,&LDA,B->V,&LDB,&(BETA),C->V,&LDC);
      free(temp);
    }
    else
    {
      /* A and B real */
      if(C->QComplex)
      {
        tempA=calloc(A->Ni*A->Nj,sizeof(Complex16));
        for(i=0;i<A->Ni*A->Nj;i++)
        {
          tempA[i].Re=A->VRe[i];
        }
        tempB=calloc(B->Ni*B->Nj,sizeof(Complex16));
        for(i=0;i<B->Ni*B->Nj;i++)
        {
          tempB[i].Re=B->VRe[i];
        }
        zgemm(&TransA,&TransB,&C->Nj,&C->Ni,&Nk,&(ALPHA),tempA,&LDA,tempB,&LDB,&(BETA),C->V,&LDC);
        free(tempA);
        free(tempB);
      }
      else
      {
        dgemm(&TransA,&TransB,&C->Nj,&C->Ni,&Nk,&(ALPHA.Re),A->VRe,&LDA,B->VRe,&LDB,&(BETA.Re),C->VRe,&LDC);
      }
    }
  }
  return 0;
}


/*   D = A ^{TransA} * B ^{TransB} * C ^{TransC}
 *   D should not be allocated before the routine is called
 *   TransA can be either:
 *          N for nothing
 *          T for transposed
 *          C for conjugate transposed
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
uint32_t TransformCompactMatrix(CompactMatrixType const *A, CompactMatrixType const *T, CompactMatrixType * const TdagAT)
{
  return CompactMatrixMultiply3('C', 'N', 'N', T, A, T, TdagAT);
}

/*  Calculates TdagAT = T^dag A T
 *
 */
uint32_t TransformCompactMatrixPreAllocated(CompactMatrixType const *A, CompactMatrixType const *T, CompactMatrixType * const TdagAT)
{
  return CompactMatrixMultiply3PreAllocated('C', 'N', 'N', T, A, T, TdagAT);
}


/*  Calculates A <= T^dag A T
 *
 */
uint32_t TransformCompactMatrixInline(CompactMatrixType * const A, CompactMatrixType const *T)
{
  CompactMatrixType temp;
  if(CompactMatrixMultiply('N', 'N', A, T, &temp)){printf("Error in TransformCompactMatrixInline\n"); return 1;}
  if(CompactMatrixMultiplyPreAllocated('C', 'N', T, &temp, A)){printf("Error in TransformCompactMatrixInline\n"); return 1;}
  FreeCompactMatrix(&temp);
  return 0;
}

/*  Calculates A <= T A T^dag
 *
 */
uint32_t TransformCompactMatrixInlineNNC(CompactMatrixType * const A, CompactMatrixType const *T)
{
  CompactMatrixType temp;
  if(CompactMatrixMultiply('N', 'C', A, T, &temp)){printf("Error in TransformCompactMatrixInlineNNC\n"); return 1;}
  if(CompactMatrixMultiplyPreAllocated('N', 'N', T, &temp, A)){printf("Error in TransformCompactMatrixInlineNNC\n"); return 1;}
  FreeCompactMatrix(&temp);
  return 0;
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
      swap=malloc(M->Ni*M->Nj*sizeof(Complex16));
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
      swap=malloc(M->Ni*M->Nj*sizeof(double));
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


/* M=M^{-1}
 */
double CompactMatrixInverse(CompactMatrixType * const M)
{
  int *IPIV;
  int LWORK = M->Ni*M->Nj;
  int INFO;
  if(M->Ni!=M->Nj)
  {
    printf("In CompactMatrixInverse tried to inverse a rectangular matrix Ni=%i, Nj=%i\n",M->Ni,M->Nj);
    return -1;
  }
  IPIV = malloc((M->Ni+1)*sizeof(int));
  if(M->QComplex)
  {
    Complex16 *WORK;
    zgetrf(&M->Ni,&M->Ni,M->V,&M->Ni,IPIV,&INFO);
    if(INFO!=0)
    {
      printf("zgetrf failed in CompactMatrixInverse with INFO %i\n",INFO);
      return 1;
    }
    WORK = malloc(LWORK*sizeof(Complex16));
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
    WORK = malloc(LWORK*sizeof(double));
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


/*  on input A contains the matrix to be factorized
 *  on output A is destroyed (memory freed)
 *  and A = Q * R
 *  with Q orthogonal and R upper triangular
 */
uint32_t CompactMatrixQRDecomposition(CompactMatrixType *A, CompactMatrixType *Q, CompactMatrixType *R)
{
  int LWORK=A->Ni*A->Nj,INFO;
  
  
  if(A->QComplex)
  {
    
  }
  else
  {
    double *TAU,*WORK;
    TAU=malloc(Min(A->Ni,A->Nj*sizeof(double)));
    WORK=malloc(A->Ni*A->Nj*sizeof(double));
    dgeqrf(&(A->Ni),&(A->Nj),A->VRe,&(A->Ni),TAU,WORK,&LWORK,&INFO);
    
    
  }
  
  printf("To Do\n");
  return 0;
}

/*  A is replaced by an orthonormal form of A using Singular value decomposition (eigen value decomposition)
 *  If a SingularValue is smaller than SingularValue this vector is removed from the basis
 *  Aortho = T A
 *  A = Tinv Aortho
 *  if A has no singular values Tinv . T is the identity matrix
 *  the trance of Tinv . T is equal to the dimension (number of rows) of A minus the number of singular vectors in A.
 */
uint32_t CompactMatrixSVDOrthogonalizeRow(CompactMatrixType * A, double const SingularValue, CompactMatrixType *T, CompactMatrixType *Tinv)
{
  uint32_t i;
  char *JOBU ="A";  /* in SVD A = U S V^T return full array U */
  char *JOBVT="O";  /* in SVD A = U S V^T return the first min(Ni,Nj) rows of V^T in matrix A, which is over written */
  int LWORK, INFO;
  LWORK = Max(5,A->Ni)*Max(5,A->Nj);
  double *WORK;
  Complex16 *WORKIM;
  double *VT=NULL;
  double *S;
  Tinv->Ni=A->Nj;
  Tinv->Nj=A->Nj;
  S=calloc(A->Nj,sizeof(double));
  Tinv->QComplex=A->QComplex;
  InitCompactMatrix(Tinv);
  if(A->QComplex)
  {
    WORK = malloc(5*Min(A->Ni,A->Nj)*sizeof(double));
    if(WORK==NULL)
    {
      printf("malloc failed in CompactMatrixSVDOrthogonalize\n");
      return 1;
    }
    WORKIM = malloc(LWORK*sizeof(Complex16));
    if(WORKIM==NULL)
    {
      free(WORK);
      printf("malloc failed in CompactMatrixSVDOrthogonalize\n");
      return 1;
    }
    zgesvd(JOBU, JOBVT, &A->Nj, &A->Ni, A->V, &A->Nj,S,Tinv->V, &A->Nj,VT, &A->Nj, WORKIM, &LWORK,WORK, &INFO);
    if(INFO)
    {
      printf("zgesvd failed in CompactMatrixSVDOrthogonalize with INFO %i\n",INFO);
      return 1;
    }
    free(WORKIM);
    free(WORK);
  }
  else
  {
    WORK = malloc(LWORK*sizeof(double));
    if(WORK==NULL)
    {
      printf("malloc failed in CompactMatrixSVDOrthogonalize\n");
      return 1;
    }
    dgesvd(JOBU, JOBVT, &A->Nj, &A->Ni, A->VRe, &A->Nj,S,Tinv->VRe, &A->Nj,VT, &A->Nj, WORK, &LWORK, &INFO);
    if(INFO)
    {
      printf("dgesvd failed in CompactMatrixSVDOrthogonalize with INFO %i\n",INFO);
      return 1;
    }
    free(WORK);
  }
  CompactMatrixCopy(*Tinv, T);
  CompactMatrixConjugateTranspose(T);
  for(i=0;i<A->Nj;i++)
  {
    if(S[i]<SingularValue)
    {
      S[i]=0;
    }
  }
  MultiplyCompactMatrixDiagonalMatrix(Tinv, S);
  for(i=0;i<A->Nj;i++)
  {
    if(S[i]<SingularValue)
    {
      S[i]=0;
    }
    else
    {
      S[i]=1.0/S[i];
    }
  }
  MultiplyDiagonalMatrixCompactMatrix(S,T);
  for(i=0;i<A->Nj;i++)
  {
    if(S[i]<100*DBL_MIN)
    {
      S[i]=0;
    }
    else
    {
      S[i]=1;
    }
  }
  MultiplyDiagonalMatrixCompactMatrix(S,A);
  free(S);
  return 0;
}


/*  A is replaced by an orthonormal form of A using Singular value decomposition (eigen value decomposition)
 *  If a SingularValue is smaller than SingularValue this vector is removed from the basis
 *  Aortho = A T
 *  A = Aortho Tinv
 *  if A has no singular values Tinv . T is the identity matrix
 *  the trace of Tinv . T is equal to the dimension (number of columns) of A minus the number of singular vectors in A.
 */
uint32_t CompactMatrixSVDOrthogonalizeCol(CompactMatrixType * A, double const SingularValue, CompactMatrixType *T, CompactMatrixType *Tinv)
{
  uint32_t i;
  char *JOBU ="O";  /* in SVD A = U S V^T return the first min(Ni,Nj) collumns of U in matrix A, which is over written */
  char *JOBVT="A";  /* in SVD A = U S V^T return full array VT */
  int LWORK, INFO;
  LWORK = Max(5,A->Ni)*Max(5,A->Nj);
  double *WORK;
  Complex16 *WORKIM;
  double *U=NULL;
  double *S;
  Tinv->Ni=A->Ni;
  Tinv->Nj=A->Ni;
  S=calloc(A->Ni,sizeof(double));
  Tinv->QComplex=A->QComplex;
  InitCompactMatrix(Tinv);
  if(A->QComplex)
  {
    WORK = malloc(5*Min(A->Ni,A->Nj)*sizeof(double));
    if(WORK==NULL)
    {
      printf("malloc failed in CompactMatrixSVDOrthogonalize\n");
      return 1;
    }
    WORKIM = malloc(LWORK*sizeof(Complex16));
    if(WORKIM==NULL)
    {
      free(WORK);
      printf("malloc failed in CompactMatrixSVDOrthogonalize\n");
      return 1;
    }
    zgesvd(JOBU, JOBVT, &A->Nj, &A->Ni, A->V, &A->Nj,S,U, &A->Ni,Tinv->V, &A->Ni, WORKIM, &LWORK, WORK, &INFO);
    if(INFO)
    {
      printf("zgesvd failed in CompactMatrixSVDOrthogonalize with INFO %i\n",INFO);
      return 1;
    }
    free(WORKIM);
    free(WORK);
  }
  else
  {
    WORK = malloc(LWORK*sizeof(double));
    if(WORK==NULL)
    {
      printf("malloc failed in CompactMatrixSVDOrthogonalize\n");
      return 1;
    }
    dgesvd(JOBU, JOBVT, &A->Nj, &A->Ni, A->VRe, &A->Nj,S,U, &A->Ni,Tinv->VRe, &A->Ni, WORK, &LWORK, &INFO);
    if(INFO)
    {
      printf("dgesvd failed in CompactMatrixSVDOrthogonalize with INFO %i\n",INFO);
      return 1;
    }
    free(WORK);
  }
  CompactMatrixCopy(*Tinv, T);
  CompactMatrixConjugateTranspose(T);
  for(i=0;i<A->Ni;i++)
  {
    if(S[i]<SingularValue)
    {
      S[i]=0;
    }
  }
  MultiplyDiagonalMatrixCompactMatrix(S,Tinv);
  for(i=0;i<A->Ni;i++)
  {
    if(S[i]<SingularValue)
    {
      S[i]=0;
    }
    else
    {
      S[i]=1.0/S[i];
    }
  }
  MultiplyCompactMatrixDiagonalMatrix(T,S);
  for(i=0;i<A->Ni;i++)
  {
    if(S[i]<100*DBL_MIN)
    {
      S[i]=0;
    }
    else
    {
      S[i]=1;
    }
  }
  MultiplyCompactMatrixDiagonalMatrix(A,S);
  
  
  
  free(S);
  return 0;
}

/*  A is replaced by an orthonormal form of A using Gram Schmidt
 *  If the norm of a vector A_{i} is smaller than SingularValue this vector is set to zero
 */
uint32_t CompactMatrixGramSchmidtOrthogonalizeCol(CompactMatrixType * A, double const SingularValue)
{
  int i,j,ii,ipos,iipos;
  double dot,norm,dotRe,dotIm;
  if(A->QComplex)
  {
    for(i=0;i<A->Ni;i++)
    {
      ipos=i*A->Nj;
      /* normalize vector i */
      norm=0;
      for(j=0;j<A->Nj;j++)
      {
        norm+=Sqr(A->V[ipos+j].Re) + Sqr(A->V[ipos+j].Im);
      }
      if(norm<SingularValue)
      {
        norm=0;
      }
      else
      {
        norm=1.0/sqrt(norm);
      }
      for(j=0;j<A->Nj;j++)
      {
        A->V[ipos+j].Re*=norm;
        A->V[ipos+j].Im*=norm;
      }
      /* for all vectors subtract projection on vector i */
      for(ii=i+1;ii<A->Ni;ii++)
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
  else
  {
    for(i=0;i<A->Ni;i++)
    {
      ipos=i*A->Nj;
      /* normalize vector i */
      norm=0;
      for(j=0;j<A->Nj;j++)
      {
        norm+=Sqr(A->VRe[ipos+j]);
      }
      if(norm<SingularValue)
      {
        norm=0;
      }
      else
      {
        norm=1.0/sqrt(norm);
      }
      for(j=0;j<A->Nj;j++)
      {
        A->VRe[ipos+j]*=norm;
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
  return 0;
}

/*  Anew = A T
 *  A = Anew Tinv
 *  with Anew^+ Anew an identity matrix
 */
uint32_t CompactMatrixGramSchmidtOrthogonalizeColplusT(CompactMatrixType const * A, CompactMatrixType *const Anew, CompactMatrixType *const T, CompactMatrixType *const Tinv, double const SingularValue)
{
  int i,ii,ipos,iipos, one=1;
  CompactMatrixCopy(*A, Anew);
  T->Ni=A->Ni;
  T->Nj=A->Ni;
  T->QComplex=A->QComplex;
  InitCompactIdentityMatrix(T);
  if(A->QComplex)
  {
    
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
      if(norm<SingularValue)
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
  CompactMatrixMultiply('C', 'N', Anew, A, Tinv);return 0;
}


uint32_t CompactMatrixLowdinOrthogonalize(CompactMatrixType * A)
{
  printf("To Do\n");
  return 1;
}

/*
 *  A=fun.H.fun^{\dag} with H a diagonal matrix with the values val on the diagonal
 *  input M defined, T not initiallized, val allocated
 */
uint32_t CompactMatrixEigensystem(CompactMatrixType const *A, CompactMatrixType * const fun, double *val)
{
  int    info;
  const char *jobz = "V";
  const char *uplo = "U";
  int    lwork = 1 + 6 * A->Ni + 2 * Sqr(A->Ni);
  int    liwork = 3 + 5 * A->Ni;
  int    *iwork;
  double *work;
  
  iwork = malloc(liwork*sizeof(int));
  if(iwork==NULL)
  {
    printf("malloc for iwork failed in CompactMatrixEigensystem\n");
    return 1;
  }
  work = malloc(lwork*sizeof(double));
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
      printf("dsyv failed in CompactMatrixEigensystem with info value %i\n", info);
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
        printf("also dsyvd failed in CompactMatrixEigensystem with info value %i\n", info);
        printf("ERORR!!!!!!!!!!\n");
        fflush(stdout);
        return 1;
      }
      return 0;
    }
    return 0;
  }
  else
  {                             /* complex */
    int    lcwork = 2 * A->Ni + Sqr(A->Ni);
    Complex16 *cwork;
    cwork =malloc(lcwork*sizeof(Complex16));
    if(cwork==NULL)
    {
      printf("malloc for cwork failed in CompactMatrixEigensystem\n");
      return 1;
    }
    if (CompactMatrixCopy(*A, fun))
      return 1;
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
        return 1;
      }
    }
    return 0;
  }
  return 0;
}

/*  On entrance A is ...
 *  SingularValue is the smallest positive eigenvalue of the matrix A^dag A that is considered non-singular
 *
 *  returned:
 *  A = q . r
 *  with q orthnoromal and r hermitian
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
  CompactMatrixMultiply3('N', 'N', 'C', &fun, &Asqr, &fun, r);
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
  CompactMatrixMultiply4('N','N','C','C',&fun, &Asqr, &fun, A,q);
  FreeCompactMatrix(&Asqr);
  FreeCompactMatrix(&fun);
  return NSingular;
}


/*  Given an semi definite hermitian matrix, CompactMatrixSqrt calculates its sqrt, such that
 *   B * B^{dag} = A
 *
 */
uint32_t CompactMatrixSqrt(CompactMatrixType const *A, CompactMatrixType * const SqrtA)
{
  int i;
  double val[A->Ni];
  CompactMatrixType fun;
  
  CompactMatrixEigensystem(A, &fun, val);
  for(i=0;i<fun.Ni;i++)
  {
    if(val[i]<0)
    {
      printf("Negative value in CompactMatrixSqrt, eigen value is %22.15E at pos %i\n",val[i],i);
      val[i]=0;
    }
    val[i]=sqrt(val[i]);
  }
  TransformDiagonalMatrixByCompactMatrixTATdag(val, &fun, SqrtA);
  FreeCompactMatrix(&fun);
  return 0;
}


/*  Given an semi definite hermitian matrix, CompactMatrixSqrt calculates its sqrt, such that
 *   B * B^{dag} = A
 *
 */
uint32_t CompactMatrixSqrtold(CompactMatrixType const *A, CompactMatrixType *SqrtA)
{
  int i,j;
  double val[A->Ni],sqrtval;
  CompactMatrixType fun,B;
  
  CompactMatrixEigensystem(A, &fun, val);
  CompactMatrixCopy(fun, &B);
  
  if(B.QComplex)
  {
    for(i=0;i<fun.Ni;i++)
    {
      if(val[i]<0)
      {
        printf("Negative value in CompactMatrixSqrt, eigen value is %22.15E at pos %i\n",val[i],i);
        val[i]=0;
      }
      sqrtval=sqrt(val[i]);
      {
        for(j=i*B.Nj;j<fun.Nj+i*B.Nj;j++)
        {
          B.V[j].Re*=sqrtval;
          B.V[j].Im*=sqrtval;
        }
      }
    }
  }
  else
  {
    for(i=0;i<fun.Ni;i++)
    {
      if(val[i]<0)
      {
        printf("Negative value in CompactMatrixSqrt, eigen value is %22.15E at pos %i\n",val[i],i);
        val[i]=0;
      }
      sqrtval=sqrt(val[i]);
      {
        for(j=i*B.Nj;j<fun.Nj+i*B.Nj;j++)
        {
          B.VRe[j]*=sqrtval;
        }
      }
    }
  }
  
  CompactMatrixMultiply('N', 'C', &B, &fun, SqrtA);
  
  FreeCompactMatrix(&fun);
  FreeCompactMatrix(&B);
  
  return 0;
}

/*
 *
 */
uint32_t BlockNaturalImpurityOrbitalMatrixToCompactMatrix(BlockNaturalImpurityOrbitalMatrixType const * A, CompactMatrixType * const B)
{
  int i,j,k,apos,bposL,bposR;
  B->QComplex=A->QComplex;
  B->Ni=(Max(1,A->val.N)+Max(1,A->con.N))*A->BlockSize;
  B->Nj=(A->val.N+A->con.N)*A->BlockSize;
  if(InitCompactMatrix(B))
  {
    printf("InitCompactMatrix failed in BlockNaturalImpurityOrbitalMatrixToCompactMatrix\n");
    return 1;
  }
  if(A->QComplex)
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
uint32_t CompactMatrixToSpectra(CompactMatrixType const *A, SpectraType * const G, uint32_t const BlockSize, double Emin, double Emax, uint32_t N)
{
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
  val=malloc(A->Ni*sizeof(double));
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
    Gtemp=malloc(BlockSize*BlockSize*sizeof(SpectraType));
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
          Gtemp[i*BlockSize+j].Im[iEmin] -= wmin *(fun.V[k*fun.Nj+i].Re*fun.V[k*fun.Nj+j].Re + fun.V[k*fun.Nj+i].Im*fun.V[k*fun.Nj+j].Im);
          Gtemp[i*BlockSize+j].Im[iEplus]-= wplus*(fun.V[k*fun.Nj+i].Re*fun.V[k*fun.Nj+j].Re + fun.V[k*fun.Nj+i].Im*fun.V[k*fun.Nj+j].Im);
          Gtemp[i*BlockSize+j].Re[iEmin] += wmin *(fun.V[k*fun.Nj+i].Re*fun.V[k*fun.Nj+j].Im - fun.V[k*fun.Nj+i].Im*fun.V[k*fun.Nj+j].Re);
          Gtemp[i*BlockSize+j].Re[iEplus]+= wplus*(fun.V[k*fun.Nj+i].Re*fun.V[k*fun.Nj+j].Im - fun.V[k*fun.Nj+i].Im*fun.V[k*fun.Nj+j].Re);
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
          G[i*BlockSize+j].Im[iEmin] -= wmin *fun.VRe[k*fun.Nj+i]*fun.VRe[k*fun.Nj+j];
          G[i*BlockSize+j].Im[iEplus]-= wplus*fun.VRe[k*fun.Nj+i]*fun.VRe[k*fun.Nj+j];
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

/* on input A is a Hermitian matrix, Singular Value is the smallest number larger than 0 that has numerical relavance, BlockSize is the size of the blocks and NBlockMax the maximal number of Blocks
 *  onoutput Matrix B has NBlockMax blocks (a_{0} to a_{NBlockMax+1}. b_{0} to b_{NBlockMax}
 */
uint32_t CompactMatrixBlockBandTridiagonalizeHermitianNoReorthogonalization(CompactMatrixType const *A, BlockTriDiagonalMatrixType *const B, double const SingularValue, uint32_t NBlockMax, uint32_t BlockSize)
{
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
  CompactMatrixAddMultiply('C', 'C', &Vm1, &alpha, &R, -1.0, 0.0);
  
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
    CompactMatrixAddMultiply('C', 'C', &Vm2, &beta, &R, -1.0, 0.0);
    CompactMatrixAddMultiply('C', 'C', &V,   &alpha, &R, -1.0, 0.0);
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


