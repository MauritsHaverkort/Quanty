#include "BasicMath_BlokAndersonMatrixVector.h"

/*
*     ( A[0][0]=a1   | A[0][1]=b1   | A[1][2]=b2   | A[2][3]=b3 )
*   A=( A[1][0]=b1^+ | A[1][1]=a2   | 0            | 0          )
*     ( A[2][1]=b2^+ | 0            | A[2][2]=a3   | 0          )
*     ( A[3][2]=b3^+ | 0            | 0            | A[3][3]=a4 )
*
*   with a and b matrices
*
*   (A B)[i][j] = \sum_k A[i][k] . B[k][j]
*
*/



void FreeBlokAndersonMatrix(BlokAndersonMatrixType * const M)
{
  uint32_t i;
  if (M->N != 0)
  {
    for (i = 0; i < M->N - 1; i++)
    {
      FreeMatrix(&(M->a[i]));
      FreeMatrix(&(M->b[i]));
    }
    FreeMatrix(&(M->a[M->N - 1]));
  }
  free(M->a);
  free(M->b);
}

void FreeFrameBlokAndersonMatrix(BlokAndersonMatrixType * const M)
{
  free(M->a);
  free(M->b);
}

/* Initiallize a matrix M must already be allocated 
*  M.N must be set
*  bloksize gives the dimension of the matrices that build the bloks
*/
uint32_t InitBlokAndersonMatrix(BlokAndersonMatrixType * const M, uint32_t bloksize, uint32_t QComplex)
{
  uint32_t i, fail = 0;
  strcpy(M->Name, "Matrix");
  if (M->N == 0)
    return 2;
  M->a = malloc(M->N * sizeof(MatrixType));
  if (M->a == NULL)
    return 1;
  M->b = malloc((M->N - 1) * sizeof(MatrixType));
  if (M->b == NULL)
  {
    free(M->a);
    return 1;
  }
  for (i = 0; i < M->N - 1; i++)
  {
    M->a[i].QComplex = QComplex;
    M->a[i].Ni = bloksize;
    M->a[i].Nj = bloksize;
    M->b[i].QComplex = QComplex;
    M->b[i].Ni = bloksize;
    M->b[i].Nj = bloksize;
    if (InitMatrix(&(M->a[i])))
    {
      fail = i;
      break;
    }
    if (InitMatrix(&(M->b[i])))
    {
      FreeMatrix(&(M->a[i]));
      fail = i;
      break;
    };
  }
  if (!fail)
  {
    M->a[M->N].QComplex = QComplex;
    M->a[M->N].Ni = bloksize;
    M->a[M->N].Nj = bloksize;
    if (InitMatrix(&(M->a[M->N])))
    {
      fail = M->N;
    }
  }
  if (fail)
  {
    for (i = 0; i < fail; i++)
    {
      FreeMatrix(&(M->a[i]));
      FreeMatrix(&(M->b[i]));
    }
    free(M->a);
    free(M->b);
  }
  return 0;
}

/* Initiallize a matrix M must already be allocated 
*  M.N must be set
*  bloksize gives the dimension of the matrices that build the bloks
*/
uint32_t InitFrameBlokAndersonMatrix(BlokAndersonMatrixType * const M)
{
  strcpy(M->Name, "Matrix");
  if (M->N == 0)
    return 2;
  M->a = malloc(M->N * sizeof(MatrixType));
  if (M->a == NULL)
    return 1;
  M->b = malloc((M->N - 1) * sizeof(MatrixType));
  if (M->b == NULL)
  {
    free(M->a);
    return 1;
  }
  return 0;
}

/* Initiallize a matrix M must already be allocated 
*  M.N must be set
*  bloksize gives the dimension of the matrices that build the bloks
*/
uint32_t InitIdentityBlokAndersonMatrix(BlokAndersonMatrixType * const M, uint32_t bloksize, uint32_t QComplex)
{
  uint32_t i, fail = 0;
  strcpy(M->Name, "I_Matrix");
  if (M->N == 0)
    return 2;
  M->a = malloc(M->N * sizeof(MatrixType));
  if (M->a == NULL)
    return 1;
  M->b = malloc((M->N - 1) * sizeof(MatrixType));
  if (M->b == NULL)
  {
    free(M->a);
    return 1;
  }
  for (i = 0; i < M->N - 1; i++)
  {
    M->a[i].QComplex = QComplex;
    M->a[i].Ni = bloksize;
    M->a[i].Nj = bloksize;
    M->b[i].QComplex = QComplex;
    M->b[i].Ni = bloksize;
    M->b[i].Nj = bloksize;
    if (InitIdentityMatrix(&(M->a[i])))
    {
      fail = i;
      break;
    }
    if (InitMatrix(&(M->b[i])))
    {
      FreeMatrix(&(M->a[i]));
      fail = i;
      break;
    };
  }
  if (!fail)
  {
    M->a[M->N - 1].QComplex = QComplex;
    M->a[M->N - 1].Ni = bloksize;
    M->a[M->N - 1].Nj = bloksize;
    if (InitIdentityMatrix(&(M->a[M->N - 1])))
    {
      fail = M->N - 1;
    }
  }
  if (fail)
  {
    for (i = 0; i < fail; i++)
    {
      FreeMatrix(&(M->a[i]));
      FreeMatrix(&(M->b[i]));
    }
    free(M->a);
    free(M->b);
  }
  return 0;
}

/*  set all ellements with absolute values smaller then epsilon to zero
 */
void BlokAndersonMatrixChop(BlokAndersonMatrixType * const M, double const epsilon)
{
  uint32_t i;
  for (i = 0; i < M->N - 1; i++)
  {
    MatrixChop(&(M->a[i]), epsilon);
    MatrixChop(&(M->b[i]), epsilon);
  }
  MatrixChop(&(M->a[M->N]), epsilon);
}


/*
*  M<- M + alpha B
*/
void BlokAndersonMatrixAdd(BlokAndersonMatrixType * const M, BlokAndersonMatrixType B, double alphaRe)
{
  uint32_t i;
  for (i = 0; i < M->N - 1; i++)
  {
    MatrixAdd(&(M->a[i]), B.a[i], alphaRe, 0);
    MatrixAdd(&(M->b[i]), B.b[i], alphaRe, 0);
  }
  MatrixAdd(&(M->a[M->N]), B.a[M->N], alphaRe, 0);
}

/*
*  M<- alpha M
*/
void BlokAndersonMatrixScalarMultiply(BlokAndersonMatrixType * const M, double alphaRe)
{
  uint32_t i;
  for (i = 0; i < M->N - 1; i++)
  {
    MatrixScalarMultiply(&(M->a[i]), alphaRe, 0);
    MatrixScalarMultiply(&(M->b[i]), alphaRe, 0);
  }
  MatrixScalarMultiply(&(M->a[M->N]), alphaRe, 0);
}

/* Create a new BlokAndersonMatrix that has dimensions (2+iEnd-iStart,2+iEnd-iStart)
 * and is a copy of a part of A
 * namely A[{0,iStart;;iEnd}][{0,iStart;;iEnd}]
 * M not allocated on call of the subroutine
 */
uint32_t BlokAndersonMatrixPart(BlokAndersonMatrixType * const M, BlokAndersonMatrixType A, const uint32_t iStart, const uint32_t iEnd)
{
  uint32_t i;
  if ((iEnd < iStart) || (iEnd >= A.N) || (iStart == 0))
    return 1;
  M->N = 2 + iEnd - iStart;
  InitFrameBlokAndersonMatrix(M);
  strcpy(M->Name, A.Name);
  if (MatrixCopy(&(M->a[0]), A.a[0]))
  {
    FreeFrameBlokAndersonMatrix(M);
    return 1;
  }
  for (i = 1; i < M->N; i++)
  {
    MatrixCopy(&(M->a[i]), A.a[i + iStart]);
    MatrixCopy(&(M->b[i - 1]), A.b[i - 1 + iStart]);
  }
  return 0;
}

/*
*  M <- A
*/
uint32_t BlokAndersonMatrixCopy(BlokAndersonMatrixType * const M, BlokAndersonMatrixType A)
{
  uint32_t i;
  M->N = A.N;
  if (InitBlokAndersonMatrix(M, A.a[0].Ni, A.a[0].QComplex))
    return 1;
  for (i = 0; i < M->N - 1; i++)
  {
    MatrixCopy(&(M->a[i]), A.a[i]);
    MatrixCopy(&(M->b[i]), A.b[i]);
  }
  MatrixCopy(&(M->a[M->N]), A.a[M->N]);
  return 0;
}

/*
*  M <- M^T
*  M.a[i] <- M.a[i]^T
*  M.b[i] <- M.b[i]^+
*/
void BlokAndersonMatrixTranspose(BlokAndersonMatrixType * const M)
{
  uint32_t i;
  for (i = 0; i < M->N - 1; i++)
  {
    MatrixTranspose(&(M->a[i]));
    MatrixComplexConjugate(&(M->b[i]));
  }
  MatrixTranspose(&(M->a[M->N]));
}


/*
*  M <- M^*
*  M.a[i] <- M.a[i]^*
*  M.b[i] <- M.b[i]^*
*/
void BlokAndersonMatrixComplexConjugate(BlokAndersonMatrixType * const M)
{
  uint32_t i;
  for (i = 0; i < M->N - 1; i++)
  {
    MatrixComplexConjugate(&(M->a[i]));
    MatrixComplexConjugate(&(M->b[i]));
  }
  MatrixComplexConjugate(&(M->a[M->N]));
}


/*
*  M <- M^+
*  M.a[i] <- M.a[i]^+
*  M.b[i] <- M.b[i]^T
*/
void BlokAndersonMatrixConjugateTranspose(BlokAndersonMatrixType * const M)
{
  uint32_t i;
  for (i = 0; i < M->N - 1; i++)
  {
    MatrixConjugateTranspose(&(M->a[i]));
    MatrixTranspose(&(M->b[i]));
  }
  MatrixConjugateTranspose(&(M->a[M->N]));
}



/* y <- M x
 * y already allocated
 * y[i]=sum_j M[i][j] x[j]
 */
uint32_t RealBlokAndersonMatrixVector(BlokAndersonMatrixType const M, const double *const x, double *const y)
{
  uint32_t i, StartOfBloki[M.N], StartOfBlokj[M.N], Ni, Nj, MaxNi, MaxNj;
  StartOfBloki[0] = 0;
  StartOfBlokj[0] = 0;
  MaxNi = M.a[0].Ni;
  MaxNj = M.a[0].Nj;
  for (i = 1; i < M.N; i++)
  {
    StartOfBloki[i] = StartOfBloki[i - 1] + M.a[i - 1].Ni;
    StartOfBlokj[i] = StartOfBlokj[i - 1] + M.a[i - 1].Nj;
    if (M.a[i].Ni > MaxNi)
      MaxNi = M.a[i].Ni;
    if (M.a[i].Nj > MaxNj)
      MaxNj = M.a[i].Nj;
  }
  Ni = StartOfBloki[M.N - 1] + M.a[M.N - 1].Ni;
  Nj = StartOfBlokj[M.N - 1] + M.a[M.N - 1].Nj;
  for (i = 0; i < M.N; i++)
  {
    if (RealMatrixVector(M.a[i], &x[StartOfBlokj[i]], &y[StartOfBloki[i]]))
    {
      printf("Function RealMatrixVector nags in RealBlokAndersonMatrixVector\n");
      printf("Called with matrix name %s\n", M.Name);
      printf("for blok M.a[%i]\n", i);
      return 1;
    }
  }
  double temp[Max(MaxNi, MaxNj)];
  for (i = 1; i < M.N; i++)
  {
    if (RealMatrixVector(M.b[i - 1], &x[StartOfBlokj[i]], temp))
    {
      printf("Function RealMatrixVector nags in RealBlokAndersonMatrixVector\n");
      printf("Called with matrix name %s\n", M.Name);
      printf("for blok M.b[%i-1]\n", i);
      return 1;
    }
    RealVectorAdd(&y[StartOfBloki[0]], temp, M.b[i - 1].Ni, 1);
    if (RealTransposeMatrixVector(M.b[i - 1], &x[StartOfBloki[0]], temp))
    {
      printf("Function RealTransposeMatrixVector nags in RealBlokAndersonMatrixVector\n");
      printf("Called with matrix name %s\n", M.Name);
      printf("for blok M.b[%i-1]\n", i);
      return 1;
    }
    RealVectorAdd(&y[StartOfBlokj[i - 1]], temp, M.b[i - 1].Nj, 1);
  }
  return 0;
}

/* y <- M x
 * y already allocated
 * y[i]=sum_j M[i][j] x[j]
 */
uint32_t ComplexBlokAndersonMatrixVector(BlokAndersonMatrixType const M, const double *const xRe, const double *const xIm, double *const yRe, double *const yIm)
{
  uint32_t i, StartOfBloki[M.N], StartOfBlokj[M.N], Ni, Nj, MaxNi, MaxNj;
  StartOfBloki[0] = 0;
  StartOfBlokj[0] = 0;
  MaxNi = M.a[0].Ni;
  MaxNj = M.a[0].Nj;
  for (i = 1; i < M.N; i++)
  {
    StartOfBloki[i] = StartOfBloki[i - 1] + M.a[i - 1].Ni;
    StartOfBlokj[i] = StartOfBlokj[i - 1] + M.a[i - 1].Nj;
    if (M.a[i].Ni > MaxNi)
      MaxNi = M.a[i].Ni;
    if (M.a[i].Nj > MaxNj)
      MaxNj = M.a[i].Nj;
  }
  Ni = StartOfBloki[M.N - 1] + M.a[M.N - 1].Ni;
  Nj = StartOfBlokj[M.N - 1] + M.a[M.N - 1].Nj;
  for (i = 0; i < M.N; i++)
  {
    if (ComplexMatrixVector
        (M.a[i], &xRe[StartOfBlokj[i]], &xIm[StartOfBlokj[i]], &yRe[StartOfBloki[i]], &yIm[StartOfBloki[i]]))
    {
      printf("Function ComplexMatrixVector nags in ComplexBlokAndersonMatrixVector\n");
      printf("Called with matrix name %s\n", M.Name);
      printf("for blok M.a[%i]\n", i);
      return 1;
    }
  }
  double tempRe[Max(MaxNi, MaxNj)];
  double tempIm[Max(MaxNi, MaxNj)];
  for (i = 1; i < M.N; i++)
  {
    if (ComplexMatrixVector(M.b[i - 1], &xRe[StartOfBlokj[i]], &xIm[StartOfBlokj[i]], tempRe, tempIm))
    {
      printf("Function ComplexMatrixVector nags in ComplexBlokAndersonMatrixVector\n");
      printf("Called with matrix name %s\n", M.Name);
      printf("for blok M.b[%i-1]\n", i);
      return 1;
    }
    ComplexVectorAdd(&yRe[StartOfBloki[0]], &yIm[StartOfBloki[0]], tempRe, tempIm, M.b[i - 1].Ni, 1, 0);
    if (ComplexConjugateTransposeMatrixVector(M.b[i - 1], &xRe[StartOfBloki[0]], &xIm[StartOfBloki[0]], tempRe, tempIm))
    {
      printf("Function ComplexConjugateTransposeMatrixVector nags in ComplexBlokAndersonMatrixVector\n");
      printf("Called with matrix name %s\n", M.Name);
      printf("for blok M.b[%i-1]\n", i);
      return 1;
    }
    ComplexVectorAdd(&yRe[StartOfBlokj[i - 1]], &yIm[StartOfBlokj[i - 1]], tempRe, tempIm, M.b[i - 1].Nj, 1, 0);
  }
  return 0;
}

/* y <- x M x
 * y already allocated
 * y=sum_ij x[i] M[i][j] x[j]
 */
uint32_t RealVectorBlokAndersonMatrixVector(BlokAndersonMatrixType const M, const double *const x, double *const y)
{
  uint32_t i, StartOfBloki[M.N], StartOfBlokj[M.N], Ni, Nj, MaxNi, MaxNj;
  double ytemp;
  *y = 0;
  StartOfBloki[0] = 0;
  StartOfBlokj[0] = 0;
  MaxNi = M.a[0].Ni;
  MaxNj = M.a[0].Nj;
  for (i = 1; i < M.N; i++)
  {
    StartOfBloki[i] = StartOfBloki[i - 1] + M.a[i - 1].Ni;
    StartOfBlokj[i] = StartOfBlokj[i - 1] + M.a[i - 1].Nj;
    if (M.a[i].Ni > MaxNi)
      MaxNi = M.a[i].Ni;
    if (M.a[i].Nj > MaxNj)
      MaxNj = M.a[i].Nj;
  }
  double temp[Max(MaxNi, MaxNj)];
  Ni = StartOfBloki[M.N - 1] + M.a[M.N - 1].Ni;
  Nj = StartOfBlokj[M.N - 1] + M.a[M.N - 1].Nj;
  for (i = 0; i < M.N; i++)
  {
    if (RealVectorMatrixVector(M.a[i], &x[StartOfBlokj[i]], &ytemp))
    {
      printf("Function RealVectorMatrixVector nags in RealVectorBlokAndersonMatrixVector\n");
      printf("Called with matrix name %s\n", M.Name);
      printf("for blok M.a[%i]\n", i);
      return 1;
    }
    *y += ytemp;
  }
  for (i = 1; i < M.N; i++)
  {
    if (RealMatrixVector(M.b[i - 1], &x[StartOfBlokj[i]], temp))
    {
      printf("Function RealMatrixVector nags in RealVectorBlokAndersonMatrixVector\n");
      printf("Called with matrix name %s\n", M.Name);
      printf("for blok M.b[%i-1]\n", i);
      return 1;
    }
    RealVectorDotVector(&x[StartOfBloki[0]], temp, M.b[i - 1].Ni, &ytemp);
    *y += ytemp;
    if (RealTransposeMatrixVector(M.b[i - 1], &x[StartOfBloki[0]], temp))
    {
      printf("Function RealTransposeMatrixVector nags in RealVectorBlokAndersonMatrixVector\n");
      printf("Called with matrix name %s\n", M.Name);
      printf("for blok M.b[%i-1]\n", i);
      return 1;
    }
    RealVectorDotVector(&x[StartOfBlokj[i - 1]], temp, M.b[i - 1].Nj, &ytemp);
    *y += ytemp;
  }
  return 0;
}

/* y <- x^* M x
 * y already allocated
 * y=sum_ij x[i]^* M[i][j] x[j]
 */
uint32_t ComplexVectorBlokAndersonMatrixVector(BlokAndersonMatrixType const M, const double *const xRe,const double *const xIm, double *const yRe, double *const yIm)
{
  uint32_t i, StartOfBloki[M.N], StartOfBlokj[M.N], Ni, Nj, MaxNi, MaxNj;
  double ytempRe, ytempIm;
  *yRe = 0;
  *yIm = 0;
  StartOfBloki[0] = 0;
  StartOfBlokj[0] = 0;
  MaxNi = M.a[0].Ni;
  MaxNj = M.a[0].Nj;
  for (i = 1; i < M.N; i++)
  {
    StartOfBloki[i] = StartOfBloki[i - 1] + M.a[i - 1].Ni;
    StartOfBlokj[i] = StartOfBlokj[i - 1] + M.a[i - 1].Nj;
    if (M.a[i].Ni > MaxNi)
      MaxNi = M.a[i].Ni;
    if (M.a[i].Nj > MaxNj)
      MaxNj = M.a[i].Nj;
  }
  double tempRe[Max(MaxNi, MaxNj)], tempIm[Max(MaxNi, MaxNj)];
  Ni = StartOfBloki[M.N - 1] + M.a[M.N - 1].Ni;
  Nj = StartOfBlokj[M.N - 1] + M.a[M.N - 1].Nj;
  for (i = 0; i < M.N; i++)
  {
    if (ComplexVectorMatrixVector(M.a[i], &xRe[StartOfBlokj[i]], &xIm[StartOfBlokj[i]], &ytempRe, &ytempIm))
    {
      printf("Function ComplexVectorMatrixVector nags in ComplexVectorBlokAndersonMatrixVector\n");
      printf("Called with matrix name %s\n", M.Name);
      printf("for blok M.a[%i]\n", i);
      return 1;
    }
    *yRe += ytempRe;
    *yIm += ytempIm;
  }
  for (i = 1; i < M.N; i++)
  {
    if (ComplexMatrixVector(M.b[i - 1], &xRe[StartOfBlokj[i]], &xIm[StartOfBlokj[i]], tempRe, tempIm))
    {
      printf("Function ComplexMatrixVector nags in ComplexVectorBlokAndersonMatrixVector\n");
      printf("Called with matrix name %s\n", M.Name);
      printf("for blok M.b[%i-1]\n", i);
      return 1;
    }
    ComplexVectorDotVector(&xRe[StartOfBloki[0]], &xIm[StartOfBloki[0]], tempRe, tempIm, M.b[i - 1].Ni, &ytempRe,
                           &ytempIm);
    *yRe += ytempRe;
    *yIm += ytempIm;
    if (ComplexConjugateTransposeMatrixVector(M.b[i - 1], &xRe[StartOfBloki[0]], &xIm[StartOfBloki[0]], tempRe, tempIm))
    {
      printf("Function ComplexConjugateTransposeMatrixVector nags in ComplexVectorBlokAndersonMatrixVector\n");
      printf("Called with matrix name %s\n", M.Name);
      printf("for blok M.b[%i-1]\n", i);
      return 1;
    }
    ComplexVectorDotVector(&xRe[StartOfBlokj[i - 1]], &xIm[StartOfBlokj[i - 1]], tempRe, tempIm, M.b[i - 1].Nj,
                           &ytempRe, &ytempIm);
    *yRe += ytempRe;
    *yIm += ytempIm;
  }
  return 0;
}
