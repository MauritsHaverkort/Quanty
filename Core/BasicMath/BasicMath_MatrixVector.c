/*
 *     ( A[0][0] | A[0][1] | A[0][2] )
 *   A=( A[1][0] | A[1][1] | A[1][2] )
 *     ( A[2][0] | A[2][1] | A[2][2] )
 *
 *   (A B)[i][j] = \sum_k A[i][k] B[k][j]
 *
 *   y[j] = \sum_i A[i][j] x[i]
 *
 *   A[i][j]=A[0][i*A.Nj+j]
 */

/* allocate the memory for a matrix in pointer format 
 * M[i] is a pointer to a vector representing row i 
 * Ni rows of lenght Nj 
 */
double **MatrixMalloc(const int Ni, const int Nj)
{
  double **M;
  int    i;
  M = (double**)malloc(Ni * sizeof(double *));
  if (M == NULL)
  {
    printf("MatrixMalloc Failed for pointers of matrix of size %i by %i\n", Ni, Nj);
    fflush(stdout);
    return NULL;
  }
  else
  {
    M[0] = (double*)malloc(Ni * Nj * sizeof(double));
    if (M[0] == NULL)
    {
      printf("MatrixMalloc Failed for matrix of size %i by %i\n", Ni, Nj);
      fflush(stdout);
      free(M);
      return NULL;
    }
    else
    {
      for (i = 1; i < Ni; i++)
      {
        M[i] = &M[0][i * Nj];
      }
      return M;
    }
  }
}

/* like MatrixMalloc, but added zero fill
 */
double **MatrixCalloc(const int Ni, const int Nj)
{
  double **M;
  int    i;
  M = (double**)malloc(Ni * sizeof(double *));
  if (M == NULL)
  {
    printf("MatrixCalloc Failed for pointers of matrix of size %i by %i\n", Ni, Nj);
    fflush(stdout);
    return NULL;
  }
  else
  {
    M[0] = (double*)calloc(Ni * Nj, sizeof(double));
    if (M[0] == NULL)
    {
      printf("MatrixCalloc Failed for matrix of size %i by %i\n", Ni, Nj);
      fflush(stdout);
      free(M);
      return NULL;
    }
    else
    {
      for (i = 1; i < Ni; i++)
      {
        M[i] = &M[0][i * Nj];
      }
      return M;
    }
  }
}

/* Free a matrix M itselve (i.e. the structure is not freed)
 */
void FreeMatrix(MatrixType * const M)
{
  free(M->VRe[0]);
  free(M->VRe);
  if (M->QComplex)
  {
    free(M->VIm[0]);
    free(M->VIm);
  }
  M->Ni = 0;
  M->Nj = 0;
}

/* Initiallize a matrix M must already be allocated 
*  M.QComplex, M.Ni, M,Nj must be set
*/
uint32_t InitMatrix(MatrixType * const M)
{
  strcpy(M->Name, "Matrix");
  M->mu = 0.;
  M->VRe = MatrixCalloc(M->Ni, M->Nj);
  if (M->VRe == NULL)
  {
    printf("MatrixCalloc failed in InitMatrix\n");
    fflush(stdout);
    return 1;
  }
  if (M->QComplex)
  {
    M->VIm = MatrixCalloc(M->Ni, M->Nj);
    if (M->VIm == NULL)
    {
      printf("MatrixCalloc failed in InitMatrix for complex allocation\n");
      fflush(stdout);
      M->QComplex = 0;
      FreeMatrix(M);
      M->QComplex = 1;
      return 1;
    }
  }
  return 0;
}

/* Initiallize a matrix M must already be allocated 
*  M.QComplex, M.Ni, M,Nj must be set
*/
uint32_t InitIdentityMatrix(MatrixType * const M)
{
  uint32_t i;
  strcpy(M->Name, "I_Matrix");
  M->mu = 0.;
  M->VRe = MatrixCalloc(M->Ni, M->Nj);
  if (M->VRe == NULL)
  {
    return 1;
  }
  if (M->QComplex)
  {
    M->VIm = MatrixCalloc(M->Ni, M->Nj);
    if (M->VIm == NULL)
    {
      M->QComplex = 0;
      FreeMatrix(M);
      M->QComplex = 1;
      return 1;
    }
  }
  for (i = 0; i < Min(M->Ni, M->Nj); i++)
  {
    M->VRe[i][i] = 1;
  }
  return 0;
}

/* turn a real matrix into a complex matrix i.e. allocate more memory and zero this
 */
uint32_t MakeRealMatrixComplex(MatrixType * const M)
{
  if (M->QComplex)
    return 0;
  M->VIm = MatrixCalloc(M->Ni, M->Nj);
  if (M->VIm == NULL)
    return 1;
  M->QComplex = 1;
  return 0;
}

/* turn a complex matrix into a real matrix i.e. free memory
 */
uint32_t MakeComplexMatrixReal(MatrixType * const M)
{
  if (M->QComplex==0)
    return 0;
  free(M->VIm);
  M->QComplex = 0;
  return 0;
}

/* set elements with an absolute value smaller then epsilon to zero
 */
void MatrixChop(MatrixType * const M, double const epsilon)
{
  int    i, j, QStillComplex;
  QStillComplex = 0;
  if (M->QComplex)
  {
    for (i = 0; i < M->Ni; i++)
    {
      for (j = 0; j < M->Nj; j++)
      {
        if (fabs(M->VRe[i][j]) < epsilon)
        {
          M->VRe[i][j] = 0;
        };
        if (fabs(M->VIm[i][j]) < epsilon)
        {
          M->VIm[i][j] = 0;
        }
        else
          QStillComplex = 1;
      }
    }
    if (QStillComplex == 0)
    {
      free(M->VIm[0]);
      free(M->VIm);
      M->QComplex = 0;
    }
  }
  else
  {
    for (i = 0; i < M->Ni; i++)
    {
      for (j = 0; j < M->Nj; j++)
      {
        if (fabs(M->VRe[i][j]) < epsilon)
          M->VRe[i][j] = 0;
      }
    }
  }
}

/* M=M^{T}
 */
uint32_t MatrixTranspose(MatrixType * const M)
{
  int    i, j;
  double temp;
  strncat(M->Name, "^T", MAX(0, (int)((int)MAXStrLength - (int)strlen(M->Name))));
  if (M->Ni == M->Nj)
  {                             /* Ni==Nj */
    if (M->QComplex)
    {                           /* Ni==Nj complex */
      for (i = 1; i < M->Ni; i++)
      {
        for (j = 0; j < i; j++)
        {
          temp = M->VRe[i][j];
          M->VRe[i][j] = M->VRe[j][i];
          M->VRe[j][i] = temp;
          temp = M->VIm[i][j];
          M->VIm[i][j] = M->VIm[j][i];
          M->VIm[j][i] = temp;
        }
      }
    }
    else
    {                           /* Ni==Nj Real */
      for (i = 1; i < M->Ni; i++)
      {
        for (j = 0; j < i; j++)
        {
          temp = M->VRe[i][j];
          M->VRe[i][j] = M->VRe[j][i];
          M->VRe[j][i] = temp;
        }
      }
    }
  }
  else
  {                             /* Ni != Nj */
    MatrixType MT;
    MT.Ni = M->Nj;
    MT.Nj = M->Ni;
    strcpy(MT.Name, M->Name);
    MT.QComplex = M->QComplex;
    if (InitMatrix(&MT))
      return 1;
    if (MT.QComplex)
    {
      for (i = 0; i < M->Ni; i++)
      {
        for (j = 0; j < M->Nj; j++)
        {
          MT.VRe[j][i] = M->VRe[i][j];
          MT.VIm[j][i] = M->VIm[i][j];
        }
      }
    }
    else
    {
      for (i = 0; i < M->Ni; i++)
      {
        for (j = 0; j < M->Nj; j++)
        {
          MT.VRe[j][i] = M->VRe[i][j];
        }
      }
    }
    FreeMatrix(M);
    M->Ni = MT.Ni;
    M->Nj = MT.Nj;
    M->QComplex = MT.QComplex;
    M->VRe = MT.VRe;
    M->VIm = MT.VIm;
  }
  return 0;
}

/* M=M^{*}
 */
void MatrixComplexConjugate(MatrixType * const M)
{
  strncat(M->Name, "^*", MAX(0, (int)((int)MAXStrLength - (int)strlen(M->Name))));
  if (M->QComplex)
  {
    int    i, j;
    for (i = 0; i < M->Ni; i++)
    {
      for (j = 0; j < M->Nj; j++)
      {
        M->VIm[i][j] = -M->VIm[i][j];
      }
    }
  }
}

/* M=M^{\dag}
 */
uint32_t MatrixConjugateTranspose(MatrixType * M)
{
  int    i, j;
  double temp;
  strncat(M->Name, "^+", MAX(0, (int)((int)MAXStrLength - (int)strlen(M->Name))));
  if (M->Ni == M->Nj)
  {                             /* Ni==Nj */
    if (M->QComplex)
    {                           /* Ni==Nj complex */
      for (i = 1; i < M->Ni; i++)
      {
        for (j = 0; j < i; j++)
        {
          temp = M->VRe[i][j];
          M->VRe[i][j] = M->VRe[j][i];
          M->VRe[j][i] = temp;
          temp = M->VIm[i][j];
          M->VIm[i][j] = -M->VIm[j][i];
          M->VIm[j][i] = -temp;
        }
        M->VIm[i][i] *= -1;
      }
      M->VIm[0][0] *= -1;
    }
    else
    {                           /* Ni==Nj Real */
      for (i = 1; i < M->Ni; i++)
      {
        for (j = 0; j < i; j++)
        {
          temp = M->VRe[i][j];
          M->VRe[i][j] = M->VRe[j][i];
          M->VRe[j][i] = temp;
        }
      }
    }
  }
  else
  {                             /* Ni != Nj */
    MatrixType MT;
    MT.Ni = M->Nj;
    MT.Nj = M->Ni;
    MT.QComplex = M->QComplex;
    if (InitMatrix(&MT))
      return 1;
    if (MT.QComplex)
    {
      for (i = 0; i < M->Ni; i++)
      {
        for (j = 0; j < M->Nj; j++)
        {
          MT.VRe[j][i] = M->VRe[i][j];
          MT.VIm[j][i] = -M->VIm[i][j];
        }
      }
    }
    else
    {
      for (i = 0; i < M->Ni; i++)
      {
        for (j = 0; j < M->Nj; j++)
        {
          MT.VRe[j][i] = M->VRe[i][j];
        }
      }
    }
    FreeMatrix(M);
    M->Ni = MT.Ni;
    M->Nj = MT.Nj;
    M->QComplex = MT.QComplex;
    M->VRe = MT.VRe;
    M->VIm = MT.VIm;
  }
  return 0;
}

/* M=M^{-1}
 */
double MatrixInverse(MatrixType * M)
{
  if (M->Ni != M->Nj)
  {
    printf("Tried to invert a non square matrix Ni=%i Nj=%i\n", M->Ni, M->Nj);
    return -1;
  }
  strncat(M->Name, "^-1", MAX(0, (int)((int)MAXStrLength - (int)strlen(M->Name))));
  if (M->QComplex)
  {
    int    i, j, k, jmaxcol;
    double singular;
    double MInvRe[M->Ni][M->Ni];
    double MInvIm[M->Ni][M->Ni];
    double temp, TIm, TRe;
    for (i = 0; i < M->Ni; i++) /* Initiallize MInv to an identity Matrix */
    {
      for (j = 0; j < M->Ni; j++)
      {
        MInvRe[i][j] = 0.0;
        MInvIm[i][j] = 0.0;
      }
      MInvRe[i][i] = 1.0;
    }
    singular = 1.0;
    for (i = 0; i < M->Ni; i++) /* loop over all rows */
    {
      jmaxcol = i;
      for (j = i + 1; j < M->Ni; j++) /* find jmaxcol */
      {
        if ((Sqr(M->VRe[j][i]) + Sqr(M->VIm[j][i])) > (Sqr(M->VRe[jmaxcol][i]) + Sqr(M->VIm[jmaxcol][i])))
        {
          jmaxcol = j;
        }
      }
      if (jmaxcol != i)         /* swapp row i and row jmaxcol */
      {
        for (j = 0; j < M->Ni; j++)
        {
          temp = M->VRe[i][j];
          M->VRe[i][j] = M->VRe[jmaxcol][j];
          M->VRe[jmaxcol][j] = temp;
          temp = M->VIm[i][j];
          M->VIm[i][j] = M->VIm[jmaxcol][j];
          M->VIm[jmaxcol][j] = temp;
          temp = MInvRe[i][j];
          MInvRe[i][j] = MInvRe[jmaxcol][j];
          MInvRe[jmaxcol][j] = temp;
          temp = MInvIm[i][j];
          MInvIm[i][j] = MInvIm[jmaxcol][j];
          MInvIm[jmaxcol][j] = temp;
        }
      }
      temp = Sqr(M->VIm[i][i]) + Sqr(M->VRe[i][i]);
      singular = Min(temp, singular);
      if (temp < 100 * DBL_MIN)
      {
        temp = 100 * DBL_MIN;
        printf("Singular matrix inversion of matrix %s\n", M->Name);
      }
      TIm = -M->VIm[i][i] / temp;
      TRe = M->VRe[i][i] / temp;
      for (j = 0; j < M->Ni; j++) /* multipli row i such that element [i][i] becomes one */
      {
        temp = TRe * M->VRe[i][j] - TIm * M->VIm[i][j];
        M->VIm[i][j] = TRe * M->VIm[i][j] + TIm * M->VRe[i][j];
        M->VRe[i][j] = temp;
        temp = TRe * MInvRe[i][j] - TIm * MInvIm[i][j];
        MInvIm[i][j] = TRe * MInvIm[i][j] + TIm * MInvRe[i][j];
        MInvRe[i][j] = temp;
      }
      for (j = 0; j < M->Ni; j++) /* subtract row i from all other rows such that collum i is only non-zero for row i */
      {
        if (j != i)
        {
          TRe = M->VRe[j][i];
          TIm = M->VIm[j][i];
          for (k = 0; k < M->Ni; k++)
          {
            M->VRe[j][k] -= TRe * M->VRe[i][k] - TIm * M->VIm[i][k];
            M->VIm[j][k] -= TRe * M->VIm[i][k] + TIm * M->VRe[i][k];
            MInvRe[j][k] -= TRe * MInvRe[i][k] - TIm * MInvIm[i][k];
            MInvIm[j][k] -= TRe * MInvIm[i][k] + TIm * MInvRe[i][k];
          }
        }
      }
    }
    for (i = 0; i < M->Ni; i++) /* The matrix M now is the identity matrix and MInv contains the inverse matrix */
    {
      for (j = 0; j < M->Ni; j++)
      {
        M->VRe[i][j] = MInvRe[i][j];
        M->VIm[i][j] = MInvIm[i][j];
      }
    }
    return singular;
  }
  else
  {
    int    i, j, k, jmaxcol;
    double singular;
    double MInvRe[M->Ni][M->Ni];
    double temp, TRe;
    for (i = 0; i < M->Ni; i++) /* Initiallize MInv to an identity Matrix */
    {
      for (j = 0; j < M->Ni; j++)
      {
        MInvRe[i][j] = 0.0;
      }
      MInvRe[i][i] = 1.0;
    }
    singular = 1.0;
    for (i = 0; i < M->Ni; i++) /* loop over all rows */
    {
      jmaxcol = i;
      for (j = i + 1; j < M->Ni; j++) /* find jmaxcol */
      {
        if (Sqr(M->VRe[j][i]) > Sqr(M->VRe[jmaxcol][i]))
        {
          jmaxcol = j;
        }
      }
      if (jmaxcol != i)         /* swapp row i and row jmaxcol */
      {
        for (j = 0; j < M->Ni; j++)
        {
          temp = M->VRe[i][j];
          M->VRe[i][j] = M->VRe[jmaxcol][j];
          M->VRe[jmaxcol][j] = temp;
          temp = MInvRe[i][j];
          MInvRe[i][j] = MInvRe[jmaxcol][j];
          MInvRe[jmaxcol][j] = temp;
        }
      }
      temp = Sqr(M->VRe[i][i]);
      singular = Min(temp, singular);
      if (temp < 100 * DBL_MIN)
      {
        temp = 100 * DBL_MIN;
        printf("Singular matrix inversion of matrix %s\n", M->Name);
      }
      TRe = M->VRe[i][i] / temp;
      for (j = 0; j < M->Ni; j++) /* multipli row i such that element [i][i] becomes one */
      {
        M->VRe[i][j] = TRe * M->VRe[i][j];
        MInvRe[i][j] = TRe * MInvRe[i][j];
      }
      for (j = 0; j < M->Ni; j++) /* subtract row i from all other rows such that collum i is only non-zero for row i */
      {
        if (j != i)
        {
          TRe = M->VRe[j][i];
          for (k = 0; k < M->Ni; k++)
          {
            M->VRe[j][k] -= TRe * M->VRe[i][k];
            MInvRe[j][k] -= TRe * MInvRe[i][k];
          }
        }
      }
    }
    for (i = 0; i < M->Ni; i++) /* The matrix M now is the identity matrix and MInv contains the inverse matrix */
    {
      for (j = 0; j < M->Ni; j++)
      {
        M->VRe[i][j] = MInvRe[i][j];
      }
    }
    return singular;
  }
  return 0;
}

/* M[i ... i+B.Ni][j ... j+B.Nj] = M[i ... i+B.Ni][j ... j+B.Nj] + alpha B
 */
uint32_t MatrixAddShifted(MatrixType * const M, const MatrixType B, double const alphaRe, double const alphaIm,
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
        M->VRe[i + ishift][j + jshift] += alphaRe * B.VRe[i][j];
      }
    }
    return 0;
  }
  if (M->QComplex == 0)
    if (MakeRealMatrixComplex(M))
      return 1;
  if (B.QComplex == 0)
  {                             /* B real, rest complex */
    for (i = 0; i < B.Ni; i++)
    {
      for (j = 0; j < B.Nj; j++)
      {
        M->VRe[i + ishift][j + jshift] += alphaRe * B.VRe[i][j];
        M->VIm[i + ishift][j + jshift] += alphaIm * B.VRe[i][j];
      }
    }
    return 0;
  }
  /* all complex */
  for (i = 0; i < B.Ni; i++)
  {
    for (j = 0; j < B.Nj; j++)
    {
      M->VRe[i + ishift][j + jshift] += alphaRe * B.VRe[i][j] - alphaIm * B.VIm[i][j];
      M->VIm[i + ishift][j + jshift] += alphaIm * B.VRe[i][j] + alphaRe * B.VIm[i][j];
    }
  }
  return 0;
}


/* M[i ... i+B.Nj][j ... j+B.Ni] = M[i ... i+B.Nj][j ... j+B.Ni] + alpha B^{\dag}
 */
uint32_t MatrixAddShiftedConjugateTranspose(MatrixType * const M, const MatrixType B, double const alphaRe,
                                            double const alphaIm, int ishift, int jshift)
{
  if ((M->Ni < B.Nj + ishift) || (M->Nj < B.Ni + jshift))
  {
    printf("Matrix B^{\\dag} shifted does not fit in matrix M\n");
    printf("M.Ni= %i M.Nj= %i B.Ni= %i B.Nj= %i ishift= %i jshift= %i\n", M->Ni, M->Nj, B.Ni, B.Nj, ishift, jshift);
    return 1;
  }
  int    i, j;
  if ((M->QComplex == 0) && (B.QComplex == 0) && (alphaIm == 0))
  {                             /* all real */
    for (i = 0; i < B.Nj; i++)
    {
      for (j = 0; j < B.Ni; j++)
      {
        M->VRe[i + ishift][j + jshift] += alphaRe * B.VRe[j][i];
      }
    }
    return 0;
  }
  if (M->QComplex == 0)
    if (MakeRealMatrixComplex(M))
      return 1;
  if (B.QComplex == 0)
  {                             /* B real, rest complex */
    for (i = 0; i < B.Nj; i++)
    {
      for (j = 0; j < B.Ni; j++)
      {
        M->VRe[i + ishift][j + jshift] += alphaRe * B.VRe[j][i];
        M->VIm[i + ishift][j + jshift] += alphaIm * B.VRe[j][i];
      }
    }
    return 0;
  }
  /* all complex */
  for (i = 0; i < B.Nj; i++)
  {
    for (j = 0; j < B.Ni; j++)
    {
      M->VRe[i + ishift][j + jshift] += alphaRe * B.VRe[j][i] + alphaIm * B.VIm[j][i];
      M->VIm[i + ishift][j + jshift] += alphaIm * B.VRe[j][i] - alphaRe * B.VIm[j][i];
    }
  }
  return 0;
}

/* M=M + alpha B
 */
uint32_t MatrixAdd(MatrixType * const M, const MatrixType B, double const alphaRe, double const alphaIm)
{
  if ((M->Ni != B.Ni) || (M->Nj != B.Nj))
  {
    printf("Can only add matrices of the same dimension M+B failed\n");
    printf("M.Ni= %i M.Nj= %i B.Ni= %i B.Nj= %i\n", M->Ni, M->Nj, B.Ni, B.Nj);
    fflush(stdout);
    return 1;
  }
  int    i, j;
  if ((M->QComplex == 0) && (B.QComplex == 0) && (alphaIm == 0))
  {                             /* all real */
    for (i = 0; i < M->Ni; i++)
    {
      for (j = 0; j < M->Nj; j++)
      {
        M->VRe[i][j] += alphaRe * B.VRe[i][j];
      }
    }
    return 0;
  }
  if (M->QComplex == 0)
    if (MakeRealMatrixComplex(M))
      return 1;
  if (B.QComplex == 0)
  {                             /* B real, rest complex */
    for (i = 0; i < M->Ni; i++)
    {
      for (j = 0; j < M->Nj; j++)
      {
        M->VRe[i][j] += alphaRe * B.VRe[i][j];
        M->VIm[i][j] += alphaIm * B.VRe[i][j];
      }
    }
    return 0;
  }
  /* all complex */
  for (i = 0; i < M->Ni; i++)
  {
    for (j = 0; j < M->Nj; j++)
    {
      M->VRe[i][j] += alphaRe * B.VRe[i][j] - alphaIm * B.VIm[i][j];
      M->VIm[i][j] += alphaIm * B.VRe[i][j] + alphaRe * B.VIm[i][j];
    }
  }
  return 0;
}

/* M = alpha M
 */
uint32_t MatrixScalarMultiply(MatrixType * const M, double const alphaRe, double const alphaIm)
{
  int    i, j;
  if ((M->QComplex == 0) && (alphaIm == 0))
  {                             /* real */
    for (i = 0; i < M->Ni; i++)
    {
      for (j = 0; j < M->Nj; j++)
      {
        M->VRe[i][j] *= alphaRe;
      }
    }
    return 0;
  }
  if (M->QComplex == 0)
    if (MakeRealMatrixComplex(M))
      return 1;
  double temp;
  /* all complex */
  for (i = 0; i < M->Ni; i++)
  {
    for (j = 0; j < M->Nj; j++)
    {
      temp = M->VRe[i][j];
      M->VRe[i][j] = alphaRe * M->VRe[i][j] - alphaIm * M->VIm[i][j];
      M->VIm[i][j] = alphaRe * M->VIm[i][j] + alphaIm * temp;
    }
  }
  return 0;
}

/* M = T^* A T^{transpose}
 * M[i][j]=\sum_k,l T[i][l]^* A[l][k] T[j][k]
 * (A.Ni==A.Nj==T.Nj)
 * M.Ni=M.Nj=T.Ni;
 * M on entrance not initiallized
 */
uint32_t MatrixTransform(MatrixType * const M, MatrixType const A, MatrixType const T)
{
  if ((A.Ni != A.Nj) || (A.Ni != T.Nj))
  {
    printf("Invallid dimensions for matrix transformation A.Ni should be T.Nj and A should be square\n");
    printf("A.Ni=%i A.Nj=%i T.Ni=%i T.Nj=%i \n", A.Ni, A.Nj, T.Ni, T.Nj);
    return 1;
  }
  uint32_t i, j, k, l;
  M->Ni = T.Ni;
  M->Nj = T.Ni;
  M->QComplex = 1;
  if ((A.QComplex == 0) && (T.QComplex == 0))
    M->QComplex = 0;
  if (InitMatrix(M))
    return 1;
  strcpy(M->Name, "U^* A U^T U=(");
  strncat(M->Name, T.Name, MAX(0, (int)((int)MAXStrLength - (int)24)) / 2);
  strcat(M->Name, ") A=(");
  strncat(M->Name, A.Name, MAX(0, (int)((int)MAXStrLength - (int)24)) / 2);
  strcat(M->Name, ")");
  if ((A.QComplex) && (T.QComplex))
  {                             /* all complex */
    for (i = 0; i < M->Ni; i++)
    {
      for (j = 0; j < M->Nj; j++)
      {
        M->VRe[i][j] = 0;
        M->VIm[i][j] = 0;
        for (k = 0; k < A.Nj; k++)
        {
          for (l = 0; l < A.Ni; l++)
          {
            M->VRe[i][j] +=
              T.VRe[i][l] * A.VRe[l][k] * T.VRe[j][k] + T.VIm[i][l] * A.VIm[l][k] * T.VRe[j][k] +
              T.VIm[i][l] * A.VRe[l][k] * T.VIm[j][k] - T.VRe[i][l] * A.VIm[l][k] * T.VIm[j][k];
            M->VIm[i][j] +=
              T.VIm[i][l] * A.VIm[l][k] * T.VIm[j][k] + T.VRe[i][l] * A.VRe[l][k] * T.VIm[j][k] +
              T.VRe[i][l] * A.VIm[l][k] * T.VRe[j][k] - T.VIm[i][l] * A.VRe[l][k] * T.VRe[j][k];
          }
        }
      }
    }
    return 0;
  }
  if (A.QComplex)
  {                             /* A complex, T real */
    for (i = 0; i < M->Ni; i++)
    {
      for (j = 0; j < M->Nj; j++)
      {
        M->VRe[i][j] = 0;
        M->VIm[i][j] = 0;
        for (k = 0; k < A.Nj; k++)
        {
          for (l = 0; l < A.Ni; l++)
          {
            M->VRe[i][j] += T.VRe[i][l] * A.VRe[l][k] * T.VRe[j][k];
            M->VIm[i][j] += T.VRe[i][l] * A.VIm[l][k] * T.VRe[j][k];
          }
        }
      }
    }
    return 0;
  }
  if (T.QComplex)
  {                             /* T complex, A real */
    for (i = 0; i < M->Ni; i++)
    {
      for (j = 0; j < M->Nj; j++)
      {
        M->VRe[i][j] = 0;
        M->VIm[i][j] = 0;
        for (k = 0; k < A.Nj; k++)
        {
          for (l = 0; l < A.Ni; l++)
          {
            M->VRe[i][j] += T.VRe[i][l] * A.VRe[l][k] * T.VRe[j][k] + T.VIm[i][l] * A.VRe[l][k] * T.VIm[j][k];
            M->VIm[i][j] += T.VRe[i][l] * A.VRe[l][k] * T.VIm[j][k] - T.VIm[i][l] * A.VRe[l][k] * T.VRe[j][k];
          }
        }
      }
    }
    return 0;
  }
  /* all real */
  for (i = 0; i < M->Ni; i++)
  {
    for (j = 0; j < M->Nj; j++)
    {
      M->VRe[i][j] = 0;
      for (k = 0; k < A.Nj; k++)
      {
        for (l = 0; l < A.Ni; l++)
        {
          M->VRe[i][j] += T.VRe[i][l] * A.VRe[l][k] * T.VRe[j][k];
        }
      }
    }
  }
  return 0;
}


/* M = T^* A T^{transpose}
 * M[i][j]=\sum_k,l T[i][l]^* A[l][k] T[j][k]
 * (A.Ni==A.Nj==T.Nj)
 * M.Ni=M.Nj=T.Ni;
 * M on entrance not initiallized
 */
uint32_t MatrixTransformPart(MatrixType * const M, MatrixType const A, MatrixType const T, const uint32_t MNiMin,
                             const uint32_t MNiMax, const uint32_t MNjMin, const uint32_t MNjMax)
{
  if ((A.Ni != A.Nj) || (A.Ni != T.Nj))
  {
    printf("Invallid dimensions for matrix transformation A.Ni should be T.Nj and A should be square\n");
    printf("A.Ni=%i A.Nj=%i T.Ni=%i T.Nj=%i \n", A.Ni, A.Nj, T.Ni, T.Nj);
    return 1;
  }
  if (MNiMin > MNiMax)
  {
    printf("MNiMin should be smaller then MNiMax\n");
    return 1;
  }
  if (MNjMin > MNjMax)
  {
    printf("MNjMin should be smaller then MNjMax\n");
    return 1;
  }
  if (MNjMax >= T.Ni)
  {
    printf("MNjMax=%i should be smaller then T.Ni=%i\n", MNjMax, T.Ni);
    return 1;
  }
  if (MNiMax >= T.Ni)
  {
    printf("MNiMax=%i should be smaller then T.Ni=%i", MNiMax, T.Ni);
    return 1;
  }
  uint32_t i, j, k, l;
  M->Ni = 1 + MNiMax - MNiMin;
  M->Nj = 1 + MNjMax - MNjMin;
  M->QComplex = 1;
  if ((A.QComplex == 0) && (T.QComplex == 0))
    M->QComplex = 0;
  if (InitMatrix(M))
    return 1;
  strcpy(M->Name, "U^* A U^T U=(");
  strncat(M->Name, T.Name, MAX(0, (int)((int)MAXStrLength - (int)24)) / 2);
  strcat(M->Name, ") A=(");
  strncat(M->Name, A.Name, MAX(0, (int)((int)MAXStrLength - (int)24)) / 2);
  strcat(M->Name, ")");
  if ((A.QComplex) && (T.QComplex))
  {                             /* all complex */
    for (i = MNiMin; i <= MNiMax; i++)
    {
      for (j = MNjMin; j <= MNjMax; j++)
      {
        M->VRe[i - MNiMin][j - MNjMin] = 0;
        M->VIm[i - MNiMin][j - MNjMin] = 0;
        for (k = 0; k < A.Nj; k++)
        {
          for (l = 0; l < A.Ni; l++)
          {
            M->VRe[i - MNiMin][j - MNjMin] +=
              T.VRe[i][l] * A.VRe[l][k] * T.VRe[j][k] + T.VIm[i][l] * A.VIm[l][k] * T.VRe[j][k] +
              T.VIm[i][l] * A.VRe[l][k] * T.VIm[j][k] - T.VRe[i][l] * A.VIm[l][k] * T.VIm[j][k];
            M->VIm[i - MNiMin][j - MNjMin] +=
              T.VIm[i][l] * A.VIm[l][k] * T.VIm[j][k] + T.VRe[i][l] * A.VRe[l][k] * T.VIm[j][k] +
              T.VRe[i][l] * A.VIm[l][k] * T.VRe[j][k] - T.VIm[i][l] * A.VRe[l][k] * T.VRe[j][k];
          }
        }
      }
    }
    return 0;
  }
  if (A.QComplex)
  {                             /* A complex, T real */
    for (i = MNiMin; i <= MNiMax; i++)
    {
      for (j = MNjMin; j <= MNjMax; j++)
      {
        M->VRe[i - MNiMin][j - MNjMin] = 0;
        M->VIm[i - MNiMin][j - MNjMin] = 0;
        for (k = 0; k < A.Nj; k++)
        {
          for (l = 0; l < A.Ni; l++)
          {
            M->VRe[i - MNiMin][j - MNjMin] += T.VRe[i][l] * A.VRe[l][k] * T.VRe[j][k];
            M->VIm[i - MNiMin][j - MNjMin] += T.VRe[i][l] * A.VIm[l][k] * T.VRe[j][k];
          }
        }
      }
    }
    return 0;
  }
  if (T.QComplex)
  {                             /* T complex, A real */
    for (i = MNiMin; i <= MNiMax; i++)
    {
      for (j = MNjMin; j <= MNjMax; j++)
      {
        M->VRe[i - MNiMin][j - MNjMin] = 0;
        M->VIm[i - MNiMin][j - MNjMin] = 0;
        for (k = 0; k < A.Nj; k++)
        {
          for (l = 0; l < A.Ni; l++)
          {
            M->VRe[i - MNiMin][j - MNjMin] +=
              T.VRe[i][l] * A.VRe[l][k] * T.VRe[j][k] + T.VIm[i][l] * A.VRe[l][k] * T.VIm[j][k];
            M->VIm[i - MNiMin][j - MNjMin] +=
              T.VRe[i][l] * A.VRe[l][k] * T.VIm[j][k] - T.VIm[i][l] * A.VRe[l][k] * T.VRe[j][k];
          }
        }
      }
    }
    return 0;
  }
  /* all real */
  for (i = MNiMin; i <= MNiMax; i++)
  {
    for (j = MNjMin; j <= MNjMax; j++)
    {
      M->VRe[i - MNiMin][j - MNjMin] = 0;
      for (k = 0; k < A.Nj; k++)
      {
        for (l = 0; l < A.Ni; l++)
        {
          M->VRe[i - MNiMin][j - MNjMin] += T.VRe[i][l] * A.VRe[l][k] * T.VRe[j][k];
        }
      }
    }
  }
  return 0;
}

/* M = A B
 * M[i][j]=\sum_k A[i][k] B[k][j]
 * if ((A.Nj!=B.Ni)) then error
 * (M.Ni=A.Ni)&&(M.Nj=B.Nj)
 * M on entrance not initiallized
 */
uint32_t MatrixMultiply(MatrixType * const M, MatrixType const A, MatrixType const B)
{
  if (A.Nj != B.Ni)
  {
    printf("Invallid dimensions for matrix multiplication A.Nj should be B.Ni\n");
    printf("A.Ni=%i A.Nj=%i B.Ni=%i B.Nj=%i \n", A.Ni, A.Nj, B.Ni, B.Nj);
    return 1;
  }
  uint32_t i, j, k;
  strcpy(M->Name, A.Name);
  strncat(M->Name, " * ", MAX(0, (int)((int)MAXStrLength - (int)strlen(A.Name))));
  strncat(M->Name, B.Name, MAX(0, (int)((int)MAXStrLength - (int)strlen(A.Name) + (int)4)));
  M->Ni = A.Ni;
  M->Nj = B.Nj;
  if ((A.QComplex) && (B.QComplex))
  {                             /* all complex */
    M->QComplex = 1;
    if (InitMatrix(M))
      return 1;
    for (i = 0; i < M->Ni; i++)
    {
      for (j = 0; j < M->Nj; j++)
      {
        M->VRe[i][j] = 0;
        M->VIm[i][j] = 0;
        for (k = 0; k < A.Nj; k++)
        {
          M->VRe[i][j] += A.VRe[i][k] * B.VRe[k][j] - A.VIm[i][k] * B.VIm[k][j];
          M->VIm[i][j] += A.VRe[i][k] * B.VIm[k][j] + A.VIm[i][k] * B.VRe[k][j];
        }
      }
    }
    return 0;
  }
  if (A.QComplex)
  {                             /* A complex, B real */
    M->QComplex = 1;
    if (InitMatrix(M))
      return 1;
    for (i = 0; i < M->Ni; i++)
    {
      for (j = 0; j < M->Nj; j++)
      {
        M->VRe[i][j] = 0;
        M->VIm[i][j] = 0;
        for (k = 0; k < A.Nj; k++)
        {
          M->VRe[i][j] += A.VRe[i][k] * B.VRe[k][j];
          M->VIm[i][j] += A.VIm[i][k] * B.VRe[k][j];
        }
      }
    }
    return 0;
  }
  if (B.QComplex)
  {                             /* B complex, A real */
    M->QComplex = 1;
    if (InitMatrix(M))
      return 1;
    for (i = 0; i < M->Ni; i++)
    {
      for (j = 0; j < M->Nj; j++)
      {
        M->VRe[i][j] = 0;
        M->VIm[i][j] = 0;
        for (k = 0; k < A.Nj; k++)
        {
          M->VRe[i][j] += A.VRe[i][k] * B.VRe[k][j];
          M->VIm[i][j] += A.VRe[i][k] * B.VIm[k][j];
        }
      }
    }
    return 0;
  }
  /* all real */
  M->QComplex = 0;
  if (InitMatrix(M))
    return 1;
  for (i = 0; i < M->Ni; i++)
  {
    for (j = 0; j < M->Nj; j++)
    {
      M->VRe[i][j] = 0;
      for (k = 0; k < A.Nj; k++)
      {
        M->VRe[i][j] += A.VRe[i][k] * B.VRe[k][j];
      }
    }
  }
  return 0;
}

/* Create a new matrix that has dimensions (1+iEnd-iStart,1+jEnd-jStart) 
 * and is a copy of a part of A
 * M not allocated on call of the subroutine
 */
uint32_t MatrixPart(MatrixType * const M, MatrixType const A, const uint32_t iStart, const uint32_t iEnd,
                    const uint32_t jStart, const uint32_t jEnd)
{
  uint32_t i, j;
  if ((iEnd < iStart) || (jEnd < jStart) || (iEnd >= A.Ni) || (jEnd >= A.Nj))
    return 1;
  M->Ni = 1 + iEnd - iStart;
  M->Nj = 1 + jEnd - jStart;
  M->QComplex = A.QComplex;
  if (InitMatrix(M))
    return 1;
  strcpy(M->Name, "Part of ");
  strncat(M->Name, A.Name, MAX(0, (int)((int)MAXStrLength - (int)10)));
  if (M->QComplex)
  {
    for (i = iStart; i <= iEnd; i++)
    {
      for (j = jStart; j <= jEnd; j++)
      {
        M->VRe[i - iStart][j - jStart] = A.VRe[i][j];
        M->VIm[i - iStart][j - jStart] = A.VIm[i][j];
      }
    }
  }
  else
  {
    for (i = iStart; i <= iEnd; i++)
    {
      for (j = jStart; j <= jEnd; j++)
      {
        M->VRe[i - iStart][j - jStart] = A.VRe[i][j];
      }
    }
  }
  return 0;
}

/*
 * M<-A
 */
uint32_t MatrixCopy(MatrixType * const M, MatrixType const A)
{
  uint32_t i, j;
  M->Ni = A.Ni;
  M->Nj = A.Nj;
  M->QComplex = A.QComplex;
  if (InitMatrix(M))
    return 1;
  strcpy(M->Name, A.Name);
  if (M->QComplex)
  {
    for (i = 0; i < M->Ni; i++)
    {
      for (j = 0; j < M->Nj; j++)
      {
        M->VRe[i][j] = A.VRe[i][j];
        M->VIm[i][j] = A.VIm[i][j];
      }
    }
  }
  else
  {
    for (i = 0; i < M->Ni; i++)
    {
      for (j = 0; j < M->Nj; j++)
      {
        M->VRe[i][j] = A.VRe[i][j];
      }
    }
  }
  return 0;
}

/*
 * M<-A
 */
uint32_t MatrixCopyToCompact(CompactMatrixType * const M, MatrixType const A)
{
  uint32_t i;
  M->Ni = A.Ni;
  M->Nj = A.Nj;
  M->QComplex = A.QComplex;
  if (InitCompactMatrix(M)) return 1;
  strcpy(M->Name, A.Name);
  if (M->QComplex)
  {
    for (i = 0; i < M->Ni*M->Nj; i++)
    {
      M->V[i].Re = A.VRe[0][i];
      M->V[i].Im = A.VIm[0][i];
    }
  }
  else
  {
    for (i = 0; i < M->Ni*M->Nj; i++)
    {
      M->VRe[i] = A.VRe[0][i];
    }
  }
  return 0;
}

/* b <- M a
 * b already allocated
 * b[i]=sum_j M[i][j] a[j]
 */
uint32_t RealMatrixVector(MatrixType const M, const double *const a, double *const b)
{
  uint32_t i, j;
  if (M.QComplex)
  {
    printf("tried to use the function RealMatrixVector on a complex Matrix\n");
    printf("Name of the matrix %s\n", M.Name);
    return 1;
  }
  for (i = 0; i < M.Ni; i++)
  {
    b[i] = 0;
    for (j = 0; j < M.Nj; j++)
    {
      b[i] += M.VRe[i][j] * a[j];
    }
  }
  return 0;
}

/* b <- M^T a
 * b already allocated
 * b[i]=sum_j M[j][i] a[j]
 */
uint32_t RealTransposeMatrixVector(MatrixType const M, const double *const a, double *const b)
{
  uint32_t i, j;
  if (M.QComplex)
  {
    printf("tried to use the function RealMatrixVector on a complex Matrix\n");
    printf("Name of the matrix %s\n", M.Name);
    return 1;
  }
  for (i = 0; i < M.Ni; i++)
  {
    b[i] = 0;
    for (j = 0; j < M.Nj; j++)
    {
      b[i] += M.VRe[j][i] * a[j];
    }
  }
  return 0;
}

/* b <- M a
 * b already allocated
 * b[i]=sum_j M[i][j] a[j]
 */
uint32_t ComplexMatrixVector(MatrixType const M, const double *const aRe, const double *const aIm, double *const bRe,
                             double *const bIm)
{
  uint32_t i, j;
  if (M.QComplex)
  {                             /* M Complex */
    for (i = 0; i < M.Ni; i++)
    {
      bRe[i] = 0;
      bIm[i] = 0;
      for (j = 0; j < M.Nj; j++)
      {
        bRe[i] += M.VRe[i][j] * aRe[j] - M.VIm[i][j] * aIm[j];
        bIm[i] += M.VRe[i][j] * aIm[j] + M.VIm[i][j] * aRe[j];
      }
    }
  }
  else
  {                             /* M Real */
    for (i = 0; i < M.Ni; i++)
    {
      bRe[i] = 0;
      bIm[i] = 0;
      for (j = 0; j < M.Nj; j++)
      {
        bRe[i] += M.VRe[i][j] * aRe[j];
        bIm[i] += M.VRe[i][j] * aIm[j];
      }
    }
  }
  return 0;
}


/* b <- M^+ a
 * b already allocated
 * b[i]=sum_j M[j][i]^* a[j]
 */
uint32_t ComplexConjugateTransposeMatrixVector(MatrixType const M, const double *const aRe, const double *const aIm,
                                               double *const bRe, double *const bIm)
{
  uint32_t i, j;
  if (M.QComplex)
  {                             /* M Complex */
    for (i = 0; i < M.Ni; i++)
    {
      bRe[i] = 0;
      bIm[i] = 0;
      for (j = 0; j < M.Nj; j++)
      {
        bRe[i] += M.VRe[j][i] * aRe[j] + M.VIm[j][i] * aIm[j];
        bIm[i] += M.VRe[j][i] * aIm[j] - M.VIm[j][i] * aRe[j];
      }
    }
  }
  else
  {                             /* M Real */
    for (i = 0; i < M.Ni; i++)
    {
      bRe[i] = 0;
      bIm[i] = 0;
      for (j = 0; j < M.Nj; j++)
      {
        bRe[i] += M.VRe[j][i] * aRe[j];
        bIm[i] += M.VRe[j][i] * aIm[j];
      }
    }
  }
  return 0;
}

/* b <- a M a
 * b already allocated
 * b=sum_ij a[i] M[i][j] a[j]
 */
uint32_t RealVectorMatrixVector(MatrixType const M, const double *const a, double *const b)
{
  uint32_t i, j;
  if (M.QComplex)
  {
    printf("tried to use the function RealVectorMatrixVector on a complex Matrix\n");
    printf("Name of the matrix %s\n", M.Name);
    return 1;
  }
  *b = 0;
  for (i = 0; i < M.Ni; i++)
  {
    for (j = 0; j < M.Nj; j++)
    {
      *b += a[i] * M.VRe[i][j] * a[j];
    }
  }
  return 0;
}

/* b <- a M a
 * b already allocated
 * b=sum_ij a[i] M[i][j] a[j]
 */
uint32_t ComplexVectorMatrixVector(MatrixType const M, const double *const aRe, const double *const aIm,
                                   double *const bRe, double *const bIm)
{
  uint32_t i, j;
  *bRe = 0;
  *bIm = 0;
  if (M.QComplex)
  {                             /* M Complex */
    for (i = 0; i < M.Ni; i++)
    {
      for (j = 0; j < M.Nj; j++)
      {
        *bRe +=
          aRe[i] * M.VRe[i][j] * aRe[j] + aIm[i] * M.VIm[i][j] * aRe[j] + aIm[i] * M.VRe[i][j] * aIm[j] -
          aRe[i] * M.VIm[i][j] * aIm[j];
        *bIm +=
          aIm[i] * M.VIm[i][j] * aIm[j] + aRe[i] * M.VRe[i][j] * aIm[j] + aRe[i] * M.VIm[i][j] * aRe[j] -
          aIm[i] * M.VRe[i][j] * aRe[j];
      }
    }
  }
  else
  {                             /* M Real */
    for (i = 0; i < M.Ni; i++)
    {
      for (j = 0; j < M.Nj; j++)
      {
        *bRe += aRe[i] * M.VRe[i][j] * aRe[j] + aIm[i] * M.VRe[i][j] * aIm[j];
        *bIm += aRe[i] * M.VRe[i][j] * aIm[j] - aIm[i] * M.VRe[i][j] * aRe[j];
      }
    }
  }
  return 0;
}

/* M transformed such that
 * M[i]^* . M[j] = delta_{i,j}
 * fully dependent vectors are set to zero
 * At the same time the overlap between the new and old M is kept "as large as possible"
 * Whereby M[0] is first optimized, then M[1], etc.
 * if M[1]=1/2 (M[0]+M[2]) then M[1] becomes M[2]
 * and M[2] is set to zero
 */
void MatrixOrthonormalize(MatrixType * const M)
{
  uint32_t i, j, Repeat;
  double tempRe, tempIm;
  if (M->QComplex)
  {                             /* complex */
    for (i = 0; i < M->Ni; i++)
    {
      Repeat = 1;
      while (Repeat)
      {
        Repeat = 0;
        for (j = 0; j < i; j++)
        {
          ComplexVectorDotVector(M->VRe[j], M->VIm[j], M->VRe[i], M->VIm[i], M->Nj, &tempRe, &tempIm);
          ComplexVectorAdd(M->VRe[i], M->VIm[i], M->VRe[j], M->VIm[j], M->Nj, -tempRe, -tempIm);
        }
        ComplexVectorDotVector(M->VRe[i], M->VIm[i], M->VRe[i], M->VIm[i], M->Nj, &tempRe, &tempIm);
        if (i >= M->Nj)
        {
          ComplexVectorMultiply(M->VRe[i], M->VIm[i], M->Nj, 0, 0);
        }
        else
        {
          if (tempRe < 0.5)
            Repeat = 1;
          tempRe = sqrt(tempRe);
          ComplexVectorMultiply(M->VRe[i], M->VIm[i], M->Nj, 1.0 / tempRe, 0);
        }
      }
    }
    return;
  }
  /* Real */
  for (i = 0; i < M->Ni; i++)
  {
    Repeat = 1;
    while (Repeat)
    {
      Repeat = 0;
      for (j = 0; j < i; j++)
      {
        RealVectorDotVector(M->VRe[j], M->VRe[i], M->Nj, &tempRe);
        RealVectorAdd(M->VRe[i], M->VRe[j], M->Nj, -tempRe);
      }
      RealVectorDotVector(M->VRe[i], M->VRe[i], M->Nj, &tempRe);
      if (i >= M->Nj)
      {
        RealVectorMultiply(M->VRe[i], M->Nj, 0);
      }
      else
      {
        if (tempRe < 0.5)
          Repeat = 1;
        tempRe = sqrt(tempRe);
        RealVectorMultiply(M->VRe[i], M->Nj, 1.0 / tempRe);
      }
    }
  }
  return;
}

/*  S contains the overlap matrix <psi_i psi_j>
 *  on return T contans a transformation matrix returning an orthonormal basis
 */
uint32_t MatrixBasisOrthonormalize(MatrixType const S, MatrixType * const T)
{

  return 0;
}
