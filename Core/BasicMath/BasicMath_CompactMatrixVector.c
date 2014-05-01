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


/* M=M^{T}
 */
uint32_t CompactMatrixTranspose(CompactMatrixType * const M)
{
  printf("to do\n");
  return 0;
}

/* M=M^{*}
 */
uint32_t CompactMatrixConjugate(CompactMatrixType * const M)
{
  printf("to do\n");
  
  return 0;
}

/* M=M^{*}
 */
uint32_t CompactMatrixConjugateTranspose(CompactMatrixType * const M)
{
  printf("to do\n");
  
  return 0;
}


/* M=M^{-1}
 */
double CompactMatrixInverse(CompactMatrixType * const M)
{
  printf("to do\n");
  
  return 0;
}


/* M[i ... i+B.Ni][j ... j+B.Nj] = M[i ... i+B.Ni][j ... j+B.Nj] + alpha B
 */
uint32_t CompactMatrixAddShifted(CompactMatrixType * const M, const CompactMatrixType B, double const alphaRe, double const alphaIm,
                          int ishift, int jshift)
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
















