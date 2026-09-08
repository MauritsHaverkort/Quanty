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
uint32_t SumGreensfunctionsTriDiagonalRepresentationNonOrtho(TriDiagonalMatrixType const A, TriDiagonalMatrixType const B, TriDiagonalMatrixType * const M)
{
  double vecA[3][A.N], vecB[3][B.N];
  double Atemp, Btemp;
  uint32_t i, j;
  M->N = A.N + B.N;
  if (InitTriDiagonalMatrix(M))
  {
    printf("InitTriDiagonalMatrix failed in SumGreensfunctionsTriDiagonalRepresentation\n");
    fflush(stdout);
    return 1;
  };
  M->a[0] = A.a[0] + B.a[0];
  M->b[0] = sqrt(Sqr(A.b[0]) + Sqr(B.b[0]));
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < A.N; j++)
    {
      vecA[i][j] = 0;
    }
    for (j = 0; j < B.N; j++)
    {
      vecB[i][j] = 0;
    }
  }
  if (M->b[0] > 10 * DBL_MIN)
  {
    vecA[0][0] = A.b[0] / M->b[0];
    vecB[0][0] = B.b[0] / M->b[0];
  }
  for (i = 1; i < M->N; i++)
  {
    RealTriDiagonalMatrixVector(A, vecA[(i - 1) % 3], vecA[i % 3]); /* vec[i] = A vec[i-1] */
    RealTriDiagonalMatrixVector(B, vecB[(i - 1) % 3], vecB[i % 3]);
    RealVectorAdd(vecA[i % 3], vecA[(i + 1) % 3], A.N, -(M->b[i - 1])); /* orthogonalization of vec[i] to vec[i-2], using that vec[i] at this point is H vec[i-1] */
    RealVectorAdd(vecB[i % 3], vecB[(i + 1) % 3], B.N, -(M->b[i - 1]));
    RealVectorDotVector(vecA[(i - 1) % 3], vecA[i % 3], A.N, &Atemp); /* a[i] = vec[i] . vec[i-1] */
    RealVectorDotVector(vecB[(i - 1) % 3], vecB[i % 3], B.N, &Btemp);
    M->a[i] = Atemp + Btemp;
    RealVectorAdd(vecA[i % 3], vecA[(i - 1) % 3], A.N, -M->a[i]); /* vec[i] -= a[i] vec[i-1] */
    RealVectorAdd(vecB[i % 3], vecB[(i - 1) % 3], B.N, -M->a[i]);
    RealVectorDotVector(vecA[i % 3], vecA[i % 3], A.N, &Atemp); /* b[i] = sqrt( vec[i] . vec[i] ) */
    RealVectorDotVector(vecB[i % 3], vecB[i % 3], B.N, &Btemp);
    M->b[i] = sqrt(Atemp + Btemp);
    if (M->b[i] > 100 * DBL_MIN)  /* vec[i] *= 1/b[i] */
    {
      RealVectorMultiply(vecA[i % 3], A.N, 1.0 / M->b[i]);
      RealVectorMultiply(vecB[i % 3], B.N, 1.0 / M->b[i]);
    }
    else
    {
      M->N = i;
      RealVectorMultiply(vecA[i % 3], A.N, 0);
      RealVectorMultiply(vecB[i % 3], B.N, 0);
    }
  }
  RealVectorTriDiagonalMatrixVector(A, vecA[(M->N - 1) % 3], &Atemp); /* a[N] = vec[N-1] A vec[N-1] */
  RealVectorTriDiagonalMatrixVector(B, vecB[(M->N - 1) % 3], &Btemp);
  M->a[M->N] = Atemp + Btemp;
  return 0;

}

uint32_t SumGreensfunctionsTriDiagonalRepresentation(TriDiagonalMatrixType const A, TriDiagonalMatrixType const B, TriDiagonalMatrixType * const M)
{
  double vecA[A.N + B.N][A.N];
  double vecB[A.N + B.N][B.N];
  double Atemp, Btemp;
  uint32_t i, j;
  M->N = A.N + B.N;
  if (InitTriDiagonalMatrix(M))
  {
    printf("InitTriDiagonalMatrix failed in SumGreensfunctionsTriDiagonalRepresentation\n");
    fflush(stdout);
    return 1;
  };
  M->a[0] = A.a[0] + B.a[0];
  M->b[0] = sqrt(Sqr(A.b[0]) + Sqr(B.b[0]));
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < A.N; j++)
    {
      vecA[i][j] = 0;
    }
    for (j = 0; j < B.N; j++)
    {
      vecB[i][j] = 0;
    }
  }
  if (M->b[0] > 10 * DBL_MIN)
  {
    vecA[0][0] = A.b[0] / M->b[0];
    vecB[0][0] = B.b[0] / M->b[0];
  }
  for (i = 1; i < M->N; i++)
  {
    RealTriDiagonalMatrixVector(A, vecA[(i - 1)], vecA[i]); /* vec[i] = A vec[i-1] */
    RealTriDiagonalMatrixVector(B, vecB[(i - 1)], vecB[i]);
    for (j = i; j > 0;)
    {
      j--;
      RealVectorDotVector(vecA[j], vecA[i], A.N, &Atemp); /* orthogonalization of vec[i] to vec[j<i-1] */
      RealVectorDotVector(vecB[j], vecB[i], B.N, &Btemp);
      RealVectorAdd(vecA[i], vecA[j], A.N, -(Atemp + Btemp));
      RealVectorAdd(vecB[i], vecB[j], B.N, -(Atemp + Btemp));
      if (j == i - 1)
      {
        M->a[i] = Atemp + Btemp;
      }
    }
    RealVectorDotVector(vecA[i], vecA[i], A.N, &Atemp); /* b[i] = sqrt( vec[i] . vec[i] ) */
    RealVectorDotVector(vecB[i], vecB[i], B.N, &Btemp);
    M->b[i] = sqrt(Atemp + Btemp);
    if (M->b[i] > 100 * (DBL_MIN + M->b[1] * DBL_EPSILON))  /* vec[i] *= 1/b[i] */
    {
      RealVectorMultiply(vecA[i], A.N, 1.0 / M->b[i]);
      RealVectorMultiply(vecB[i], B.N, 1.0 / M->b[i]);
    }
    else
    {
      M->N = i;
      M->b[i] = 0;
      RealVectorMultiply(vecA[i], A.N, 0);
      RealVectorMultiply(vecB[i], B.N, 0);
    }
  }
  RealVectorTriDiagonalMatrixVector(A, vecA[(M->N - 1)], &Atemp); /* a[N] = vec[N-1] A vec[N-1] */
  RealVectorTriDiagonalMatrixVector(B, vecB[(M->N - 1)], &Btemp);
  M->a[M->N] = Atemp + Btemp;
  return 0;
}

uint32_t SymmetrizeAndersonMatrix(AndersonMatrixType * const M)
{
  uint32_t i;
  double *tmp;
  tmp = (double*)realloc(M->a, (2 * M->N) * sizeof(double));
  if (tmp == NULL)
  {
    printf("realloc failed in SymmetrizeAndersonMatrix\n");
    return 1;
  }
  M->a = tmp;
  tmp = (double*)realloc(M->b, (Max(0, 2 * M->N - 1)) * sizeof(double));
  if (tmp == NULL)
  {
    printf("realloc failed in SymmetrizeAndersonMatrix\n");
    return 1;
  }
  M->b = tmp;
  for (i = 1; i < M->N; i++)
  {
    M->b[i] *= Sqrtp5;
    M->b[i + M->N - 1] = M->b[i];
    M->a[i + M->N] = 2 * M->mu - M->a[i + 1];
  }
  M->N *= 2;
  if (M->N != 0)
    M->N--;
  return 0;
}

uint32_t SymmetrizeBADoubleTriDiagonalMatrix(BADoubleTriDiagonalMatrixType * const M)
{
  printf("SymmetrizeBADoubleTriDiagonalMatrix not jet fully implemented\n");
  double tmp;
  tmp = sqrt(Sqr(M->bvB0) + Sqr(M->bcB0)) * Sqrtp5;
  if (M->bvB0 >= 0)
  {
    M->bvB0 = tmp;
  }
  else
  {
    M->bvB0 = -tmp;
  }
  if (M->bcB0 >= 0)
  {
    M->bcB0 = tmp;
  }
  else
  {
    M->bcB0 = -tmp;
  }
  tmp = sqrt(Sqr(M->bvA0) + Sqr(M->bcA0)) * Sqrtp5;
  if (M->bvA0 >= 0)
  {
    M->bvA0 = tmp;
  }
  else
  {
    M->bvA0 = -tmp;
  }
  if (M->bcA0 >= 0)
  {
    M->bcA0 = tmp;
  }
  else
  {
    M->bcA0 = -tmp;
  }
  return 0;
}

uint32_t SymmetrizeTridiagonalMatrix(TriDiagonalMatrixType * const M)
{
  TriDiagonalMatrixType Mplus, Mmin;
  uint32_t i;
  if (TriDiagonalMatrixCopy(&Mplus, *M))
  {
    printf("TriDiagonalMatrixCopy failed in SymmetrizeTridiagonalMatrix\n");
    fflush(stdout);
    return 1;
  };
  if (TriDiagonalMatrixCopy(&Mmin, *M))
  {
    printf("TriDiagonalMatrixCopy failed in SymmetrizeTridiagonalMatrix\n");
    fflush(stdout);
    return 1;
  };
  strcpy(Mplus.Name, M->Name);
  strcpy(Mmin.Name, M->Name);
  FreeTriDiagonalMatrix(M);
  for (i = 1; i <= Mmin.N; i++)
  {
    Mmin.a[i] *= -1;
  }
  SumGreensfunctionsTriDiagonalRepresentation(Mplus, Mmin, M);
  strcpy(M->Name, Mmin.Name);
  M->b[0] = Mplus.b[0];
  M->a[0] = Mplus.a[0];
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
  for (i = 0; i < M->N; i++)
  {
    M->b[i] *= M->b[i];
  }
  ReduceListOfPolesNumberFermi(M, N, Ef);
  for (i = 0; i < M->N; i++)
  {
    M->b[i] = sqrt(M->b[i]);
  }
}






/*  Create a new set of poles, which has lenght N, by removing the M->N-N poles with the smallest weight
 *  This algorithm combines the pole with the smallest weight with its nearest two neighbors
 *  and repeats this untill the number of nonzero poles is N
 */
void ReduceSqrtListOfPolesNumber(ListOfPoles * M, uint32_t N)
{
  //TimeType Time;
  if (M->N < N)
    return;
  //InitializeTiming(&Time);
  //BeginTiming("Total ReduceSqrtListOfPolesNumber",&Time);
  uint32_t i, imin;
  double biminLeft, biminRight, temp;
  //BeginTiming("OrderListOfPoles",&Time);
  OrderListOfPoles(M);
  //EndTiming("OrderListOfPoles",&Time);
  while (M->N > N)
  {
    imin = 0;
    for (i = 1; i < M->N; i++)
    {
      if (M->b[i] < M->b[imin])
      {
        imin = i;
      }
    }
    if (M->b[imin] > DBL_MIN)
    {
      if (imin == 0)
      {
        M->a[imin + 1] =
          (Sqr(M->b[imin] / M->b[imin + 1]) * M->a[imin] + M->a[imin + 1]) / (Sqr(M->b[imin] / M->b[imin + 1]) + 1);
        M->b[imin + 1] = sqrt(Sqr(M->b[imin]) + Sqr(M->b[imin + 1]));
      }
      else if (imin == M->N - 1)
      {
        M->a[imin - 1] =
          (Sqr(M->b[imin] / M->b[imin - 1]) * M->a[imin] + M->a[imin - 1]) / (Sqr(M->b[imin] / M->b[imin - 1]) + 1);
        M->b[imin - 1] = sqrt(Sqr(M->b[imin]) + Sqr(M->b[imin - 1]));
      }
      else
      {
        if (fabs((M->a[imin + 1] - M->a[imin - 1])) < fabs(10 * (DBL_MIN + DBL_EPSILON * M->a[imin])))
        {
          biminLeft = M->b[imin] * sqrt(0.5);
          biminRight = M->b[imin] * sqrt(0.5);
        }
        else if (fabs((M->a[imin + 1] - M->a[imin])) < fabs(10 * (DBL_MIN + DBL_EPSILON * M->a[imin])))
        {
          biminLeft = M->b[imin];
          biminRight = 0;
        }
        else if (fabs((M->a[imin] - M->a[imin - 1])) < fabs(10 * (DBL_MIN + DBL_EPSILON * M->a[imin])))
        {
          biminLeft = 0;
          biminRight = M->b[imin];
        }
        else
        {
          biminLeft = M->b[imin] * sqrt((M->a[imin + 1] - M->a[imin]) / (M->a[imin + 1] - M->a[imin - 1]));
          biminRight = M->b[imin] * sqrt((M->a[imin] - M->a[imin - 1]) / (M->a[imin + 1] - M->a[imin - 1]));
        }
        M->a[imin + 1] =
          (Sqr(biminRight / M->b[imin + 1]) * M->a[imin] + M->a[imin + 1]) / (Sqr(biminRight / M->b[imin + 1]) + 1);
        M->b[imin + 1] = sqrt(Sqr(biminRight) + Sqr(M->b[imin + 1]));
        M->a[imin - 1] =
          (Sqr(biminLeft / M->b[imin - 1]) * M->a[imin] + M->a[imin - 1]) / (Sqr(biminLeft / M->b[imin - 1]) + 1);
        M->b[imin - 1] = sqrt(Sqr(biminLeft) + Sqr(M->b[imin - 1]));
      }
    }
    for (i = imin; i < M->N - 1; i++)
    {
      M->a[i] = M->a[i + 1];
      M->b[i] = M->b[i + 1];
    }
    M->N--;
    i = imin;
    while ((i > 1) && (M->a[i - 1] < M->a[i - 2]))
    {
      temp = M->a[i - 1];
      M->a[i - 1] = M->a[i - 2];
      M->a[i - 2] = temp;
      temp = M->b[i - 1];
      M->b[i - 1] = M->b[i - 2];
      M->b[i - 2] = temp;
      i--;
    }
    i = imin;
    while ((i < M->N - 1) && (M->a[i] > M->a[i + 1]))
    {
      temp = M->a[i];
      M->a[i] = M->a[i + 1];
      M->a[i + 1] = temp;
      temp = M->b[i];
      M->b[i] = M->b[i + 1];
      M->b[i + 1] = temp;
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
  if (M->N < N)
    return;
  uint32_t i, imin;
  double biminLeft, biminRight, temp, MinMoment, TestMoment;
  OrderListOfPoles(M);
  while (M->N > N)
  {
    imin = 0;
    MinMoment = Sqr(M->b[imin]) * (M->a[imin + 1] - M->a[imin]);
    TestMoment = Sqr(M->b[M->N - 1]) * (M->a[M->N - 1] - M->a[M->N - 2]);
    if (TestMoment < MinMoment)
    {
      MinMoment = TestMoment;
      imin = M->N - 1;
    }
    for (i = 1; i < M->N - 1; i++)
    {
      TestMoment = Sqr(M->b[i]) * Min(M->a[i + 1] - M->a[i], M->a[i] - M->a[i - 1]);
      if (TestMoment < MinMoment)
      {
        imin = i;
        MinMoment = TestMoment;
      }
    }
    if (M->b[imin] > DBL_MIN)
    {
      if (imin == 0)
      {
        M->a[imin + 1] =
          (Sqr(M->b[imin] / M->b[imin + 1]) * M->a[imin] + M->a[imin + 1]) / (Sqr(M->b[imin] / M->b[imin + 1]) + 1);
        M->b[imin + 1] = sqrt(Sqr(M->b[imin]) + Sqr(M->b[imin + 1]));
      }
      else if (imin == M->N - 1)
      {
        M->a[imin - 1] =
          (Sqr(M->b[imin] / M->b[imin - 1]) * M->a[imin] + M->a[imin - 1]) / (Sqr(M->b[imin] / M->b[imin - 1]) + 1);
        M->b[imin - 1] = sqrt(Sqr(M->b[imin]) + Sqr(M->b[imin - 1]));
      }
      else
      {
        if (fabs((M->a[imin + 1] - M->a[imin - 1])) < fabs(10 * (DBL_MIN + DBL_EPSILON * M->a[imin])))
        {
          biminLeft = M->b[imin] * sqrt(0.5);
          biminRight = M->b[imin] * sqrt(0.5);
        }
        else if (fabs((M->a[imin + 1] - M->a[imin])) < fabs(10 * (DBL_MIN + DBL_EPSILON * M->a[imin])))
        {
          biminLeft = M->b[imin];
          biminRight = 0;
        }
        else if (fabs((M->a[imin] - M->a[imin - 1])) < fabs(10 * (DBL_MIN + DBL_EPSILON * M->a[imin])))
        {
          biminLeft = 0;
          biminRight = M->b[imin];
        }
        else
        {
          biminLeft = M->b[imin] * sqrt((M->a[imin + 1] - M->a[imin]) / (M->a[imin + 1] - M->a[imin - 1]));
          biminRight = M->b[imin] * sqrt((M->a[imin] - M->a[imin - 1]) / (M->a[imin + 1] - M->a[imin - 1]));
        }
        M->a[imin + 1] =
          (Sqr(biminRight / M->b[imin + 1]) * M->a[imin] + M->a[imin + 1]) / (Sqr(biminRight / M->b[imin + 1]) + 1);
        M->b[imin + 1] = sqrt(Sqr(biminRight) + Sqr(M->b[imin + 1]));
        M->a[imin - 1] =
          (Sqr(biminLeft / M->b[imin - 1]) * M->a[imin] + M->a[imin - 1]) / (Sqr(biminLeft / M->b[imin - 1]) + 1);
        M->b[imin - 1] = sqrt(Sqr(biminLeft) + Sqr(M->b[imin - 1]));
      }
    }
    for (i = imin; i < M->N - 1; i++)
    {
      M->a[i] = M->a[i + 1];
      M->b[i] = M->b[i + 1];
    }
    M->N--;
    i = imin;
    while ((i > 1) && (M->a[i - 1] < M->a[i - 2]))
    {
      temp = M->a[i - 1];
      M->a[i - 1] = M->a[i - 2];
      M->a[i - 2] = temp;
      temp = M->b[i - 1];
      M->b[i - 1] = M->b[i - 2];
      M->b[i - 2] = temp;
      i--;
    }
    i = imin;
    while ((i < M->N - 1) && (M->a[i] > M->a[i + 1]))
    {
      temp = M->a[i];
      M->a[i] = M->a[i + 1];
      M->a[i + 1] = temp;
      temp = M->b[i];
      M->b[i] = M->b[i + 1];
      M->b[i + 1] = temp;
      i++;
    }
  }
}

uint32_t ReducePolesAndersonMatrixNumberFermi(AndersonMatrixType * const A, uint32_t N, double const Ef)
{
  ListOfPoles Alp;
  if (A->N <= N)
  {
    printf("warning in ReducePolesAndersonMatrixNumber, A->N=%i and N is %i\n", A->N, N);
    return 0;
  }
  if (N <= 1)
  {
    printf("Error in ReducePolesAndersonMatrixNumber N=%i\n", N);
    return 1;
  }
  Alp.a = &(A->a[2]);
  Alp.b = &(A->b[1]);
  Alp.N = A->N - 1;
  ReduceSqrtListOfPolesNumberFermi(&Alp, N - 1, Ef);
  A->N = Alp.N + 1;
  return 0;
}

uint32_t ReducePolesAndersonMatrixNumber(AndersonMatrixType * const A, uint32_t N)
{
  ListOfPoles Alp;
  if (A->N <= N)
  {
    printf("warning in ReducePolesAndersonMatrixNumber, A->N=%i and N is %i\n", A->N, N);
    return 0;
  }
  if (N <= 1)
  {
    printf("Error in ReducePolesAndersonMatrixNumber N=%i\n", N);
    return 1;
  }
  Alp.a = &(A->a[2]);
  Alp.b = &(A->b[1]);
  Alp.N = A->N - 1;
  ReduceSqrtListOfPolesNumber(&Alp, N - 1);
  A->N = Alp.N + 1;
  return 0;
}

uint32_t ReducePolesAndersonMatrixMoment(AndersonMatrixType * const A, uint32_t N)
{
  ListOfPoles Alp;
  if (A->N <= N)
  {
    printf("warning in ReducePolesAndersonMatrixNumber, A->N=%i and N is %i\n", A->N, N);
    return 0;
  }
  if (N <= 1)
  {
    printf("Error in ReducePolesAndersonMatrixNumber N=%i\n", N);
    return 1;
  }
  Alp.a = &(A->a[2]);
  Alp.b = &(A->b[1]);
  Alp.N = A->N - 1;
  ReduceSqrtListOfPolesMoment(&Alp, N - 1);
  A->N = Alp.N + 1;
  return 0;
}

uint32_t ReducePolesAndersonMatrix(AndersonMatrixType * const A, double tMin, double Ef)
{
  uint32_t i, j;
  for (i = 2; i < A->N; i++)
  {
    if (A->b[i - 1] < 100 * DBL_MIN)
    {
      A->b[i - 1] = 0;
    }
    else if (A->a[i + 1] < Ef)
    {
      if (A->b[i - 1] < tMin)
      {
        A->a[i + 1] = (Sqr(A->b[i - 1]) * A->a[i] + Sqr(A->b[i]) * A->a[i + 1]) / (Sqr(A->b[i - 1]) + Sqr(A->b[i]));
        A->a[i] = 0;
        A->b[i] = sqrt(Sqr(A->b[i - 1]) + Sqr(A->b[i]));
        A->b[i - 1] = 0;
      }
    }
  }
  for (i = A->N; i > 2; i--)
  {
    if (A->b[i - 1] < 100 * DBL_MIN)
    {
      A->b[i - 1] = 0;
    }
    else if (A->a[i - 1] > Ef)
    {
      if (A->b[i - 1] < tMin)
      {
        A->a[i - 1] =
          (Sqr(A->b[i - 1]) * A->a[i] + Sqr(A->b[i - 2]) * A->a[i - 1]) / (Sqr(A->b[i - 1]) + Sqr(A->b[i - 2]));
        A->a[i] = 0;
        A->b[i - 2] = sqrt(Sqr(A->b[i - 1]) + Sqr(A->b[i - 2]));
        A->b[i - 1] = 0;
      }
    }
  }
  i = 2;
  while (i < A->N)
  {
    if (A->b[i] == 0)
    {
      for (j = i + 1; j < A->N; j++)
      {
        A->b[j - 1] = A->b[j];
        A->a[j] = A->a[j + 1];
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
uint32_t ReducePolesTriDiagonalMatrix(TriDiagonalMatrixType * const A, uint32_t const NpolesMax, double const Ef)
{
  double tMin;
  AndersonMatrixType An;
  if (TriDiagonalMatrixToAndersonMatrix(*A, &An))
  {
    printf("TriDiagonalMatrixToAndersonMatrix failed in ReducePolesTriDiagonalMatrix\n");
    fflush(stdout);
    return 1;
  };
  tMin = A->b[1] / sqrt((double)NpolesMax) + 100 * DBL_MIN;
  printf("For %s,\nReduce Poles around %lf, from %i ", A->Name, Ef, A->N);
  ReducePolesAndersonMatrix(&An, tMin, Ef);
  FreeTriDiagonalMatrix(A);
  AndersonMatrixToTriDiagonalMatrix(An, A, An.N);
  printf(" to %i by makeing all hybridizations larger than %lf\n", A->N, tMin);
  return 0;
}



/* M=A-B
 *
 */
uint32_t DifferenceGreensFunctionTriDiagonalRepresentation(TriDiagonalMatrixType const A, TriDiagonalMatrixType const B, TriDiagonalMatrixType * const M)
{
  ListOfPoles Ap, Bp;
  if (TriDiagonalMatrixToListOfPoles(A, &Ap))
  {
    printf("TriDiagonalMatrixToListOfPoles failed in DifferenceGreensFunctionTriDiagonalRepresentation\n");
    fflush(stdout);
    return 1;
  };
  if (TriDiagonalMatrixToListOfPoles(B, &Bp))
  {
    printf("TriDiagonalMatrixToListOfPoles failed in DifferenceGreensFunctionTriDiagonalRepresentation\n");
    fflush(stdout);
    return 1;
  };
  if (SubtractListOfPoles(&Ap, &Bp))
  {
    if ((Ap.N == 1) && (Ap.b[0] < 0))
    {
      printf
        ("SubtractListOfPoles returned a singel pole with weight %24.17E in DifferenceGreensFunctionTriDiagonalRepresentation\n",
         Ap.b[0]);
      printf("For \"%s\" - \"%s\"\n", A.Name, B.Name);
      printf("Set weight to zero\n");
      Ap.b[0] = 0;
    }
    else
    {
      printf("SubtractListOfPoles failed in DifferenceGreensFunctionTriDiagonalRepresentation\n");
      fflush(stdout);
      return 1;
    }
  }
  if (ListOfPolesToTridiagonalMatrix(Ap, M, Ap.N))
  {
    printf("ListOfPolesToTridiagonalMatrix failed in DifferenceGreensFunctionTriDiagonalRepresentation\n");
    fflush(stdout);
    return 1;
  };
  M->a[0] = A.a[0] - B.a[0];
  return 0;
}


/* S = Gbath^{-1} - Gcluster^{-1}
 *
 */
uint32_t SelfEnergyTriDiagonalRepresentation(TriDiagonalMatrixType const G0, TriDiagonalMatrixType const G, TriDiagonalMatrixType * const Sigma)
{
  TriDiagonalMatrixType IG0, IG;
  IG0.N = G0.N - 1;
  IG.N = G.N - 1;
  IG0.a = &(G0.a[1]);
  IG0.b = &(G0.b[1]);
  IG.a = &(G.a[1]);
  IG.b = &(G.b[1]);
  strcpy(IG0.Name, G0.Name);
  strncat(IG0.Name, "^-1 (-1)", Max(0, (int)((int)MAXStrLength - (int)strlen(IG0.Name))));
  strcpy(IG.Name, G.Name);
  strncat(IG.Name, "^-1 (-1)", Max(0, (int)((int)MAXStrLength - (int)strlen(IG.Name))));
  return DifferenceGreensFunctionTriDiagonalRepresentation(IG, IG0, Sigma);
}

uint32_t SelfEnergyListOfPolesFromAndersonRepresentation(AndersonMatrixType const G0, AndersonMatrixType const G, ListOfPoles * const Sigma)
{
  uint32_t i;
  ListOfPoles IG0;
  if (AndersonMatrixInverseToListOfPoles(G, Sigma))
  {
    printf("AndersonMatrixInverseToListOfPoles failed in SelfEnergyListOfPolesFromAndersonRepresentation\n");
    return 1;
  }
  if (AndersonMatrixInverseToListOfPoles(G0, &IG0))
  {
    printf("AndersonMatrixInverseToListOfPoles failed in SelfEnergyListOfPolesFromAndersonRepresentation\n");
    return 1;
  }
  if (SubtractListOfPoles(Sigma, &IG0))
  {
    printf("SubtractListOfPoles failed in SelfEnergyListOfPolesFromAndersonRepresentation\n");
    for (i = 0; i < Sigma->N; i++)
    {
      Sigma->a[i] = 0;
      Sigma->b[i] = 0;
    }
  }
  Sigma->mu = G.mu;
  strcpy(Sigma->Name, "Self energy");
  //printf("Self energy moments\n");
  //PrintListOfPolesEssentials(*Sigma);
  FreeListOfPoles(&IG0);
  return 0;
}

uint32_t SelfEnergyListOfPolesFromTriDiagonalRepresentation(TriDiagonalMatrixType const G0, TriDiagonalMatrixType const G, ListOfPoles * const Sigma)
{
  ListOfPoles IG0;
  AndersonMatrixType G0And, GAnd;
  if (TriDiagonalMatrixToAndersonMatrix(G0, &G0And))
  {
    printf("TriDiagonalMatrixToAndersonMatrix failed in SelfEnergyListOfPolesFromTriDiagonalRepresentation\n");
    return 1;
  }
  if (TriDiagonalMatrixToAndersonMatrix(G, &GAnd))
  {
    printf("TriDiagonalMatrixToAndersonMatrix failed in SelfEnergyListOfPolesFromTriDiagonalRepresentation\n");
    return 1;
  }
  if (AndersonMatrixInverseToListOfPoles(GAnd, Sigma))
  {
    printf("AndersonMatrixInverseToListOfPoles failed in SelfEnergyListOfPolesFromAndersonRepresentation\n");
    return 1;
  }
  FreeAndersonMatrix(&GAnd);
  if (AndersonMatrixInverseToListOfPoles(G0And, &IG0))
  {
    printf("AndersonMatrixInverseToListOfPoles failed in SelfEnergyListOfPolesFromAndersonRepresentation\n");
    return 1;
  }
  FreeAndersonMatrix(&G0And);
  if (SubtractListOfPoles(Sigma, &IG0))
  {
    printf("SubtractListOfPoles failed in SelfEnergyListOfPolesFromAndersonRepresentation\n");
    return 1;
  };
  FreeListOfPoles(&IG0);
  return 0;
}


uint32_t BlockSelfEnergyListOfPolesFromAndersonRepresentation(BlockAndersonMatrixType const G0, BlockAndersonMatrixType const G, BlockListOfPoles * const Sigma, double epsilon)
{
  BlockListOfPoles IG0;
  if (BlockAndersonMatrixInverseToListOfPoles(&G, Sigma, epsilon))
  {
    printf("BlockAndersonMatrixInverseToBlockListOfPoles failed in BlockSelfEnergyListOfPolesFromAndersonRepresentation\n");
    return 1;
  }
  if (BlockAndersonMatrixInverseToListOfPoles(&G0, &IG0, epsilon))
  {
    printf("BlockAndersonMatrixInverseToBlockListOfPoles failed in BlockSelfEnergyListOfPolesFromAndersonRepresentation\n");
    return 1;
  }
  if (SubtractBlockListOfPoles(Sigma, &IG0))
  {
    printf("SubtractBlockListOfPoles failed in BlockSelfEnergyListOfPolesFromAndersonRepresentation\n");
    return 1;
  }
  Sigma->mu = G.mu;
  strcpy(Sigma->Name, "Self energy");
  //printf("Self energy moments\n");
  //PrintListOfPolesEssentials(*Sigma);
  FreeBlockListOfPoles(&IG0);
  return 0;
}
/*
 *  Sigma must point to a Responsefunction object that is not initiallized.
 *  either both G and G0 are Block response functions or they are both single valued response functions
 *  epsilon the smallest energy difference for which two poles in the same block of the anderson matrix with different eigenenergies will be considered different
 */
int CalculateSelfenergy(ResponsefunctionType  * const G0, ResponsefunctionType  * const G, ResponsefunctionType * const Sigma, double epsilon)
{
  if(G0->Block==false)
  {
    ResponsefunctionType G0And, GAnd;
    AndersonMatrixType *G0AndData,*GAndData;
    if(!((G0->type=='a')||(G0->type=='A')))
    {
      ConvertResponsefunction(G0, &G0And, 'A', "Default", 0, 0, 0, 0, false, false, false);
      G0AndData = (AndersonMatrixType*)G0And.data;
    }
    else
    {
      G0AndData = (AndersonMatrixType*)G0->data;
    }
    if(!((G->type=='a')||(G->type=='A')))
    {
      ConvertResponsefunction(G, &GAnd, 'A', "Default", 0, 0, 0, 0, false, false, false);
      GAndData = (AndersonMatrixType*)GAnd.data;
    }
    else
    {
      GAndData = (AndersonMatrixType*)G->data;
    }
    Sigma->type='l';
    Sigma->Block=false;
    Sigma->NSuggestedGridPoints=0;
    Sigma->SuggestedGrid=NULL;
    Sigma->data = (void*)malloc(sizeof(ListOfPoles));
    SelfEnergyListOfPolesFromAndersonRepresentation(*G0AndData, *GAndData, (ListOfPoles*)Sigma->data);
    if(!((G0->type=='a')||(G0->type=='A')))
    {
      FreeResponsefunction(&G0And);
    }
    if(!((G->type=='a')||(G->type=='A')))
    {
      FreeResponsefunction(&GAnd);
    }
  }
  else
  {
    ResponsefunctionType G0And, GAnd;
    BlockAndersonMatrixType *G0AndData,*GAndData;
    if(!((G0->type=='a')||(G0->type=='A')))
    {
      ConvertResponsefunction(G0, &G0And, 'A', "Default", 0, 0, 0, 0, false, false, false);
      G0AndData = (BlockAndersonMatrixType*)G0And.data;
    }
    else
    {
      G0AndData = (BlockAndersonMatrixType*)G0->data;
    }
    if(!((G->type=='a')||(G->type=='A')))
    {
      ConvertResponsefunction(G, &GAnd, 'A', "Default", 0, 0, 0, 0, false, false, false);
      GAndData = (BlockAndersonMatrixType*)GAnd.data;
    }
    else
    {
      GAndData = (BlockAndersonMatrixType*)G->data;
    }
    Sigma->type='l';
    Sigma->Block=true;
    Sigma->NSuggestedGridPoints=0;
    Sigma->SuggestedGrid=NULL;
    Sigma->data = (void*)malloc(sizeof(BlockListOfPoles));
    BlockSelfEnergyListOfPolesFromAndersonRepresentation(*G0AndData, *GAndData, (BlockListOfPoles*)Sigma->data,epsilon);
    if(!((G0->type=='a')||(G0->type=='A')))
    {
      FreeResponsefunction(&G0And);
    }
    if(!((G->type=='a')||(G->type=='A')))
    {
      FreeResponsefunction(&GAnd);
    }
  }
  return 0;
}

/*  One chain of a double tri-diagonal matrix under the energy inversion
 *  G(mu+w) -> -G(mu-w)^*.
 *
 *  The chain stands in the denominator of G = a0 + b0^* ( w - Gval - Gcon )^-1 b0^T, so
 *  its head has to follow, where the head of a chain of a natural impurity orbital matrix
 *  does not: those two chains are added and their heads are a constant of G, which this
 *  map leaves alone. a[0] of the val and of the con chain together hold the onsite energy
 *  of the observed site, half of it in each, and that energy maps to 2 mu - a[0]val
 *  - a[0]con, so mu - a[0] for each of the two keeps the half and half split.
 *
 *  Below the head a[i] -> 2 mu - a[i]. The b[i] do not follow: flipping the sign of all of
 *  them is a similarity transformation with diag(1,-1,1,-1,...) which leaves the head
 *  element of the resolvent, and with it G, alone.
 */
static void DoubleTriDiagonalMatrixChainInvertEnergy(TriDiagonalMatrixType * const T, double const Mu)
{
  T->a[0] = Mu - T->a[0];
  for(int i=1;i<=T->N;i++)
  {
    T->a[i] = 2 * Mu - T->a[i];
  }
  return;
}
static void BlockDoubleTriDiagonalMatrixChainInvertEnergy(BlockTriDiagonalMatrixType * const T, double const Mu)
{
  if(T->QComplex==true)
  {
    for(int i=0;i<(T->N+1)*Sqr(T->BlockSize);i++)
    {
      T->a[i].Re*=-1;
      T->a[i].Im*=-1;
    }
    for(int j=0;j<T->BlockSizeList[0];j++)
    {
      T->a[j*(1+T->BlockSize)].Re+=Mu;
    }
    for(int i=1;i<=T->N;i++)
    {
      for(int j=0;j<T->BlockSizeList[i];j++)
      {
        T->a[i*Sqr(T->BlockSize)+j*(1+T->BlockSize)].Re+=2*Mu;
      }
    }
  }
  else
  {
    for(int i=0;i<(T->N+1)*Sqr(T->BlockSize);i++)
    {
      T->aRe[i]*=-1;
    }
    for(int j=0;j<T->BlockSizeList[0];j++)
    {
      T->aRe[j*(1+T->BlockSize)]+=Mu;
    }
    for(int i=1;i<=T->N;i++)
    {
      for(int j=0;j<T->BlockSizeList[i];j++)
      {
        T->aRe[i*Sqr(T->BlockSize)+j*(1+T->BlockSize)]+=2*Mu;
      }
    }
  }
  return;
}
/*  G(mu+w) -> -G(mu-w)^*, applied to the pole part of G only: with G = A0 + P(w), where
 *  A0 is the constant term and P vanishes at infinity, the result is A0 - P(2 mu - w)^+.
 *  The constant stays where it is, so the line above holds as it stands for A0 = 0.
 *
 *  On the data that is a -> 2 mu - a for every pole energy or chain diagonal, with the
 *  weights and the couplings left alone; for a chain, flipping the sign of all the b is a
 *  similarity transformation with diag(1,-1,1,-1,...) which leaves G alone. The formats
 *  that keep a valence and a conduction chain exchange the two, as what was below mu is
 *  above it now.
 */
int ResponseFunctionInvertEnergy(ResponsefunctionType * const Gin, ResponsefunctionType * const Gout)
{
  ResponsefunctionCopy(Gin,Gout);
  
  int i;
  double TwoMu;
  if(Gout->Block==true)
  {
    switch (Gout->type|32)
    {
      case 'a':
      {
        TwoMu = 2 * ((BlockAndersonMatrixType*)(Gout->data))->mu;
        if((((BlockAndersonMatrixType*)(Gout->data))->QComplex)==0)
        {
          for(int i = Sqr( ((BlockAndersonMatrixType*)(Gout->data))->BlockSize  ); i < Sqr( ((BlockTriDiagonalMatrixType*)(Gout->data))->BlockSize  ) + ( ((BlockAndersonMatrixType*)(Gout->data))->N) * Sqr( ((BlockAndersonMatrixType*)(Gout->data))->BlockSize  ) ; i++)
          {
            (((BlockAndersonMatrixType*)(Gout->data))->aRe[i])*=-1;
          }
          for(int i=0;i<( ((BlockAndersonMatrixType*)(Gout->data))->N);i++)
          {
            for(int j=0;j<( ((BlockAndersonMatrixType*)(Gout->data))->BlockSize);j++)
            {
              (((BlockAndersonMatrixType*)(Gout->data))->aRe[Sqr( ((BlockAndersonMatrixType*)(Gout->data))->BlockSize  ) + i * Sqr( ((BlockAndersonMatrixType*)(Gout->data))->BlockSize ) + j * (1+((BlockAndersonMatrixType*)(Gout->data))->BlockSize) ])+=TwoMu;
            }
          }
        }
        else
        {
          for(int i = Sqr( ((BlockAndersonMatrixType*)(Gout->data))->BlockSize  ); i < Sqr( ((BlockTriDiagonalMatrixType*)(Gout->data))->BlockSize  ) + ( ((BlockAndersonMatrixType*)(Gout->data))->N) * Sqr( ((BlockAndersonMatrixType*)(Gout->data))->BlockSize  ) ; i++)
          {
            (((BlockAndersonMatrixType*)(Gout->data))->a[i].Re)*=-1;
            (((BlockAndersonMatrixType*)(Gout->data))->a[i].Im)*=-1;
          }
          for(int i=0;i<( ((BlockAndersonMatrixType*)(Gout->data))->N);i++)
          {
            for(int j=0;j<( ((BlockAndersonMatrixType*)(Gout->data))->BlockSize);j++)
            {
              (((BlockAndersonMatrixType*)(Gout->data))->a[Sqr( ((BlockAndersonMatrixType*)(Gout->data))->BlockSize  ) + i * Sqr( ((BlockAndersonMatrixType*)(Gout->data))->BlockSize ) + j * (1+((BlockAndersonMatrixType*)(Gout->data))->BlockSize) ].Re)+=TwoMu;
            }
          }
        }
      }
      break;
      case 'l':
      {
        TwoMu = 2 * ((BlockListOfPoles*)(Gout->data))->mu;
        //   a_i can be found at position (i-1)*(1+BlockSize*BlockSize)
        if((((BlockListOfPoles*)(Gout->data))->QComplex)==0)
        {
          for(int i = 0; i < ( ((BlockListOfPoles*)(Gout->data))->N) * (1+Sqr( ((BlockListOfPoles*)(Gout->data))->BlockSize  )) ; i+=(1+Sqr( ((BlockListOfPoles*)(Gout->data))->BlockSize  )))
          {
            (((BlockListOfPoles*)(Gout->data))->VRe[i])*=-1;
            (((BlockListOfPoles*)(Gout->data))->VRe[i])+=TwoMu;
          }
        }
        else
        {
          for(int i = 0; i < ( ((BlockListOfPoles*)(Gout->data))->N) * (1+Sqr( ((BlockListOfPoles*)(Gout->data))->BlockSize  )) ; i+=(1+Sqr( ((BlockListOfPoles*)(Gout->data))->BlockSize  )))
          {
            (((BlockListOfPoles*)(Gout->data))->V[i].Re)*=-1;
            (((BlockListOfPoles*)(Gout->data))->V[i].Im)*=-1;
            (((BlockListOfPoles*)(Gout->data))->V[i].Re)+=TwoMu;
          }
        }
      }
      break;
      case 't':
      {
        TwoMu = 2 * ((BlockTriDiagonalMatrixType*)(Gout->data))->mu;
        if((((BlockTriDiagonalMatrixType*)(Gout->data))->QComplex)==0)
        {
          for(int i = Sqr( ((BlockTriDiagonalMatrixType*)(Gout->data))->BlockSize  ); i < Sqr( ((BlockTriDiagonalMatrixType*)(Gout->data))->BlockSize  ) + ( ((BlockTriDiagonalMatrixType*)(Gout->data))->N) * Sqr( ((BlockTriDiagonalMatrixType*)(Gout->data))->BlockSize  ) ; i++)
          {
            (((BlockTriDiagonalMatrixType*)(Gout->data))->aRe[i])*=-1;
          }
          for(int i=0;i<( ((BlockTriDiagonalMatrixType*)(Gout->data))->N);i++)
          {
            for(int j=0;j<( ((BlockTriDiagonalMatrixType*)(Gout->data))->BlockSizeList[i]);j++)
            {
              (((BlockTriDiagonalMatrixType*)(Gout->data))->aRe[Sqr( ((BlockTriDiagonalMatrixType*)(Gout->data))->BlockSize  ) + i * Sqr( ((BlockTriDiagonalMatrixType*)(Gout->data))->BlockSize ) + j * (1+((BlockTriDiagonalMatrixType*)(Gout->data))->BlockSize) ])+=TwoMu;
            }
          }
        }
        else
        {
          for(int i = Sqr( ((BlockTriDiagonalMatrixType*)(Gout->data))->BlockSize  ); i < Sqr( ((BlockTriDiagonalMatrixType*)(Gout->data))->BlockSize  ) + ( ((BlockTriDiagonalMatrixType*)(Gout->data))->N) * Sqr( ((BlockTriDiagonalMatrixType*)(Gout->data))->BlockSize  ) ; i++)
          {
            (((BlockTriDiagonalMatrixType*)(Gout->data))->a[i].Re)*=-1;
            (((BlockTriDiagonalMatrixType*)(Gout->data))->a[i].Im)*=-1;
          }
          for(int i=0;i<( ((BlockTriDiagonalMatrixType*)(Gout->data))->N);i++)
          {
            for(int j=0;j<( ((BlockTriDiagonalMatrixType*)(Gout->data))->BlockSizeList[i]);j++)
            {
              (((BlockTriDiagonalMatrixType*)(Gout->data))->a[Sqr( ((BlockTriDiagonalMatrixType*)(Gout->data))->BlockSize  ) + i * Sqr( ((BlockTriDiagonalMatrixType*)(Gout->data))->BlockSize ) + j * (1+((BlockTriDiagonalMatrixType*)(Gout->data))->BlockSize) ].Re)+=TwoMu;
            }
          }
        }
      }
      break;
      case 'n':
      {
        TwoMu = 2 * ((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data))->mu;
        
        if((((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data))->QComplex)==0)
        {
          for(int i = Sqr( ((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data))->BlockSize  ); i < Sqr( ((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data))->BlockSize  ) + ( ((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data))->val.N) * Sqr( ((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data))->BlockSize  ) ; i++)
          {
            (((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data))->val.aRe[i])*=-1;
          }
          for(int i=0;i< ( ((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data))->val.N);i++)
          {
            for(int j=0;j<( ((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data))->val.BlockSizeList[i]);j++)
            {
              (((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data))->val.aRe[Sqr( ((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data))->BlockSize  ) + i * Sqr( ((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data))->BlockSize ) + j * (1+((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data))->BlockSize) ])+=TwoMu;
            }
          }
        }
        else
        {
          for(int i = Sqr( ((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data))->BlockSize  ); i < Sqr( ((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data))->BlockSize  ) + ( ((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data))->val.N) * Sqr( ((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data))->BlockSize  ) ; i++)
          {
            (((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data))->val.a[i].Re)*=-1;
            (((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data))->val.a[i].Im)*=-1;
          }
          for(int i=0;i< ( ((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data))->val.N);i++)
          {
            for(int j=0;j<( ((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data))->val.BlockSizeList[i]);j++)
            {
              (((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data))->val.a[Sqr( ((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data))->BlockSize  ) + i * Sqr( ((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data))->BlockSize ) + j * (1+((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data))->BlockSize) ].Re)+=TwoMu;
            }
          }
        }
        if((((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data))->QComplex)==0)
        {
          for(int i = Sqr( ((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data))->BlockSize  ); i < Sqr( ((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data))->BlockSize  ) + ( ((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data))->con.N) * Sqr( ((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data))->BlockSize  ) ; i++)
          {
            (((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data))->con.aRe[i])*=-1;
          }
          for(int i=0;i< ( ((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data))->con.N);i++)
          {
            for(int j=0;j<( ((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data))->con.BlockSizeList[i]);j++)
            {
              (((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data))->con.aRe[Sqr( ((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data))->BlockSize  ) + i * Sqr( ((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data))->BlockSize ) + j * (1+((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data))->BlockSize) ])+=TwoMu;
            }
          }
        }
        else
        {
          for(int i = Sqr( ((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data))->BlockSize  ); i < Sqr( ((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data))->BlockSize  ) + ( ((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data))->con.N) * Sqr( ((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data))->BlockSize  ) ; i++)
          {
            (((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data))->con.a[i].Re)*=-1;
            (((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data))->con.a[i].Im)*=-1;
          }
          for(int i=0;i< ( ((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data))->con.N);i++)
          {
            for(int j=0;j<( ((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data))->con.BlockSizeList[i]);j++)
            {
              (((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data))->con.a[Sqr( ((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data))->BlockSize  ) + i * Sqr( ((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data))->BlockSize ) + j * (1+((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data))->BlockSize) ].Re)+=TwoMu;
            }
          }
        }
        BlockTriDiagonalMatrixType temp = ((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data))->con;
        ((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data))->con = ((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data))->val;
        ((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data))->val = temp;
        //TODO: one could implement this by swapping and multiplying by -1 sub matrices of Gin.TMT.
        SetTMTFromValandConInBlockNaturalImpurityOrbitalMatrixType( ((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data)));
      }
      break;
      case 'd':
      {
        BlockDoubleTriDiagonalMatrixType *D = (BlockDoubleTriDiagonalMatrixType*)(Gout->data);
        BlockTriDiagonalMatrixType temp;
        BlockDoubleTriDiagonalMatrixChainInvertEnergy(&(D->val), D->mu);
        BlockDoubleTriDiagonalMatrixChainInvertEnergy(&(D->con), D->mu);
        /*  a0 is the constant of G, which this map leaves where it is, and b0 enters
         *  twice with a conjugate between the two, so neither of them follows. The
         *  states that were below mu are above it now, so the two chains change place. */
        temp   = D->con;
        D->con = D->val;
        D->val = temp;
      }
      break;
      default:
      {
        printf("Error in ResponseFunctionInvertEnergy unrecognized type: %c\n",Gout->type);
        fflush(stdout);
        return 1;
      }
      break;
    }
  }
  else
  {
    switch (Gout->type|32)
    {
      case 'a':
      {
        TwoMu = 2 * ((AndersonMatrixType*)(Gout->data))->mu;
        for(i=1;i<=((AndersonMatrixType*)(Gout->data))->N;i++)
        {
          ((AndersonMatrixType*)(Gout->data))->a[i]*=- 1;
        }
        for(i=1;i<=((AndersonMatrixType*)(Gout->data))->N;i++)
        {
          ((AndersonMatrixType*)(Gout->data))->a[i]+= TwoMu;
        }
      }
      break;
      case 'l':
      {
        TwoMu = 2 * ((ListOfPoles*)(Gout->data))->mu;
        for(i=0;i<((ListOfPoles*)(Gout->data))->N;i++)
        {
          ((ListOfPoles*)(Gout->data))->a[i]*=-1;
          ((ListOfPoles*)(Gout->data))->a[i]+=TwoMu;
        }
      }
      break;
      case 't':
      {
        TwoMu = 2 * ((TriDiagonalMatrixType*)(Gout->data))->mu;
        for(i=1;i<=((TriDiagonalMatrixType*)(Gout->data))->N;i++)
        {
          ((TriDiagonalMatrixType*)(Gout->data))->a[i]*=-1;
        }
        for(i=1;i<=((TriDiagonalMatrixType*)(Gout->data))->N;i++)
        {
          ((TriDiagonalMatrixType*)(Gout->data))->a[i]+=TwoMu;
        }
      }
      break;
      case 'n':
      {
        BADoubleTriDiagonalMatrixType *B = (BADoubleTriDiagonalMatrixType*)(Gout->data);
        TriDiagonalMatrixType temp;
        TwoMu = 2 * B->mu;
        /*  a[0] of the two chains is a constant of G, which this map leaves where it is,
         *  the two chains being added and not inverted. */
        for(i=1;i<=B->val.N;i++)
        {
          B->val.a[i]*=-1;
          B->val.a[i]+=TwoMu;
        }
        for(i=1;i<=B->con.N;i++)
        {
          B->con.a[i]*=-1;
          B->con.a[i]+=TwoMu;
        }
        /*  the states that were below mu are above it now, so the two chains change
         *  place, whole: their lengths need not be equal. The seven cached numbers that
         *  describe the observed combination of the two head sites follow from the
         *  chains and are derived again, which is what the block case does as well. */
        temp   = B->con;
        B->con = B->val;
        B->val = temp;
        SetBADoubleTriDiagonalMatrixConValFromChains(B);
      }
      break;
      case 'd':
      {
        DoubleTriDiagonalMatrixType *D = (DoubleTriDiagonalMatrixType*)(Gout->data);
        TriDiagonalMatrixType temp;
        DoubleTriDiagonalMatrixChainInvertEnergy(&(D->val), D->mu);
        DoubleTriDiagonalMatrixChainInvertEnergy(&(D->con), D->mu);
        /*  a0 is the constant of G and b0 enters squared, so neither of them follows */
        temp   = D->con;
        D->con = D->val;
        D->val = temp;
      }
      break;
      default:
      {
        printf("Error in ResponseFunctionInvertEnergy unrecognized type: %c\n",Gout->type);
        fflush(stdout);
        return 1;
      }
      break;
    }
  }
  return 0;
}

/*  G(w) -> G(w)^T
 *  G(w)^T = (A0 + B0^* ( 1 / ( w - A1 - \sum_i ( B_i 1 / ( w - A_{i+1} ) . B_i^+ ) ) ) . B0^T)^T
 *         = A0^T + B0 ( 1 / ( w - A1^T - \sum_i ( B_i^* 1 / ( w - A_{i+1}^T ) . B_i^T ) ) ) . B0^+
 * B_i -> B_i^*
 * A_i -> A_i^T
 */
void BlockAndersonMatrixTypeTranspose(BlockAndersonMatrixType * const G)
{
  CompactMatrixType M;
  M.Ni = G->BlockSize;
  M.Nj = G->BlockSize;
  M.QComplex = G->QComplex;
  if(G->QComplex==true)
  {
    for(int i=0;i<=G->N;i++)
    {
      M.V = &( G->a[i*(G->BlockSize*G->BlockSize)] );
      CompactMatrixTranspose(&M);
    }
    for(int i=0;i<G->N;i++)
    {
      M.V = &( G->b[i*(G->BlockSize*G->BlockSize)] );
      CompactMatrixConjugate(&M);
    }
  }
  else
  {
    for(int i=0;i<=G->N;i++)
    {
      M.VRe = &( G->aRe[i*(G->BlockSize*G->BlockSize)] );
      CompactMatrixTranspose(&M);
    }
  }
  return;
}
void BlockListOfPolesTranspose(BlockListOfPoles * const G)
{
  CompactMatrixType M;
  M.Ni = G->BlockSize;
  M.Nj = G->BlockSize;
  M.QComplex = G->QComplex;
  if(G->QComplex==true)
  {
    M.V = G->a0;
    CompactMatrixTranspose(&M);
    for(int i=0;i<G->N;i++)
    {
      M.V = &( G->V[i*(1+G->BlockSize*G->BlockSize)+1] );
      CompactMatrixTranspose(&M);
    }
  }
  else
  {
    M.VRe = G->a0Re;
    CompactMatrixTranspose(&M);
    for(int i=0;i<G->N;i++)
    {
      M.VRe = &( G->VRe[i*(1+G->BlockSize*G->BlockSize)+1] );
      CompactMatrixTranspose(&M);
    }
  }
  return;
}
/*
 *   G0 =        A0 + B0^* . ( G1 ) . B0^T
 *   G1 = 1/(w - A1 - B1 . ( G2 ) . B1^+)
 *   G2 = 1/(w - A2 - B2 . ( G3 ) . B2^+)
 *   ...
 *   Gn = 1/(w - An)
 *   G0^T =  A0^T + B0 . ( G1^T ) . B0^*
 *   G1^T = 1/(w - A1^T - B1^* . ( G2^T ) . B1^T)
 * B_i -> B_i^*
 * A_i -> A_i^T
 */
void BlockTriDiagonalMatrixTypeTranspose(BlockTriDiagonalMatrixType * const G)
{
  CompactMatrixType M;
  M.Ni = G->BlockSize;
  M.Nj = G->BlockSize;
  M.QComplex = G->QComplex;
  if(G->QComplex==true)
  {
    for(int i=0;i<=G->N;i++)
    {
      M.V = &( G->a[i*(G->BlockSize*G->BlockSize)] );
      CompactMatrixTranspose(&M);
    }
    for(int i=0;i<G->N;i++)
    {
      M.V = &( G->b[i*(G->BlockSize*G->BlockSize)] );
      CompactMatrixConjugate(&M);
    }
  }
  else
  {
    for(int i=0;i<=G->N;i++)
    {
      M.VRe = &( G->aRe[i*(G->BlockSize*G->BlockSize)] );
      CompactMatrixTranspose(&M);
    }
  }
  return;
}
/*
 *
 */
void BlockNaturalImpurityOrbitalMatrixTypeTranspose(BlockNaturalImpurityOrbitalMatrixType * const G)
{
  BlockTriDiagonalMatrixTypeTranspose(&(G->val));
  BlockTriDiagonalMatrixTypeTranspose(&(G->con));
  CompactMatrixTranspose(&(G->a0));
  CompactMatrixTranspose(&(G->TMT));
  return;
}
/*  G(w) -> G(w)^T
 *
 *   G   = A0   + B0^* ( w - Gval(w)   - Gcon(w)   )^-1 B0^T
 *   G^T = A0^T + B0   ( w - Gval(w)^T - Gcon(w)^T )^-1 B0^+
 *
 *  the same A_i -> A_i^T and B_i -> B_i^* as everywhere else. B0 stands in the head
 *  convention, so with B0 -> B0^* the head reads B0 ( ) B0^+, and the two chains are added
 *  in the denominator, where each of them transposes as the tri-diagonal response function
 *  it is. a0 and b0 carry a QComplex of their own, which for a mixed matrix need not be the
 *  QComplex of the matrix itself; CompactMatrixConjugate reads it and does nothing to a
 *  real matrix.
 */
void BlockDoubleTriDiagonalMatrixTypeTranspose(BlockDoubleTriDiagonalMatrixType * const G)
{
  BlockTriDiagonalMatrixTypeTranspose(&(G->val));
  BlockTriDiagonalMatrixTypeTranspose(&(G->con));
  CompactMatrixTranspose(&(G->a0));
  CompactMatrixConjugate(&(G->b0));
  return;
}


int ResponseFunctionTranspose(ResponsefunctionType * const Gin, ResponsefunctionType * const Gout)
{
  ResponsefunctionCopy(Gin,Gout);
  if(Gout->Block==false) return 0;
  switch (Gout->type|32)
  {
    case 'a':
    {
      BlockAndersonMatrixTypeTranspose(((BlockAndersonMatrixType*)(Gout->data)));
    }
    break;
    case 'l':
    {
      BlockListOfPolesTranspose(((BlockListOfPoles*)(Gout->data)));
    }
    break;
    case 't':
    {
      BlockTriDiagonalMatrixTypeTranspose(((BlockTriDiagonalMatrixType*)(Gout->data)));
    }
    break;
    case 'n':
    {
      BlockNaturalImpurityOrbitalMatrixTypeTranspose(((BlockNaturalImpurityOrbitalMatrixType*)(Gout->data)));
    }
    break;
    case 'd':
    {
      BlockDoubleTriDiagonalMatrixTypeTranspose(((BlockDoubleTriDiagonalMatrixType*)(Gout->data)));
    }
    break;
    default:
    {
      printf("Error in ResponseFunctionTranspose unrecognized type: %c\n",Gout->type);
      fflush(stdout);
      return 1;
    }
    break;
  }
  return 0;
}

/*  Split a response function into the effective Hamiltonian standing in its denominator and
 *  the offset and scaling in front of it,
 *
 *    G(w) = A0 + B0^* ( w - Heff(w) )^-1 B0^T
 *
 *  Gout holds Heff, which is itself a response function; A0 and B0 come back as matrices,
 *  1 by 1 for a single valued G.
 *
 *  Which type Heff comes back as is chosen per input type by what is numerically cheapest,
 *  or most numerically stable. (Future options can influence the method used and type returned):
 *
 *    t -> t   no numerical work at all, only a shift in the data structure
 *    a -> l   one diagonalization per bath block and nothing more
 *    l -> l   through a
 *    n -> l   through a
 *    d -> l   through a
 *
 *  The three that need a transformation take the standard route, ConvertResponsefunction or
 *  ConvertBlockMatrix, so that they behave exactly as ChangeType would and take the same
 *  options.
 */
int ResponseFunctionEffectiveHamiltonian(ResponsefunctionType * const Gin, ResponsefunctionType * const Gout, CompactMatrixType * const A0, CompactMatrixType * const B0, double const T, double const SingularValue, double const MinimalPoleWeigth, double const MinimalPoleDistance)
{
  ResponsefunctionType Ga;
  bool converted = false;
  char type = Gin->type|32;
  if((type!='t')&&(type!='a')&&(type!='l')&&(type!='n')&&(type!='d'))
  {
    printf("Error in ResponseFunctionEffectiveHamiltonian: unrecognized type %c\n",Gin->type);
    fflush(stdout);
    return 1;
  }
  /*  A0 and B0 first, from Gin, whatever the route below does */
  if(Gin->Block==true)
  {
    if(BlockResponseFunctionGetA0(Gin, A0)||BlockResponseFunctionGetB0(Gin, B0))
    {
      printf("Error in ResponseFunctionEffectiveHamiltonian: could not read A0 or B0\n");
      fflush(stdout);
      return 1;
    }
  }
  else
  {
    A0->Ni=1; A0->Nj=1; A0->QComplex=0;
    B0->Ni=1; B0->Nj=1; B0->QComplex=0;
    if(InitCompactMatrix(A0)||InitCompactMatrix(B0))
    {
      printf("InitCompactMatrix failed in ResponseFunctionEffectiveHamiltonian\n");
      fflush(stdout);
      return 1;
    }
    if(ResponseFunctionGetA0(Gin, A0->VRe)||ResponseFunctionGetB0(Gin, B0->VRe))
    {
      printf("Error in ResponseFunctionEffectiveHamiltonian: could not read a0 or b0\n");
      fflush(stdout);
      FreeCompactMatrix(A0);
      FreeCompactMatrix(B0);
      return 1;
    }
  }
  /*  SingularValue and T are magnitudes and a temperature, so a negative one cannot be
   *  meant: Set standard values
   */
  double sv = (SingularValue < 0) ? 1000*DBL_EPSILON : SingularValue;
  double kt = (T < 0) ? sqrt(DBL_MIN) : T;
  Gout->Block = Gin->Block;
  Gout->NSuggestedGridPoints = 0;
  Gout->SuggestedGrid = NULL;
  /*  t keeps its own type: the effective Hamiltonian of a chain is that chain read one
   *  level down, so there is nothing to compute */
  if(type=='t')
  {
    Gout->type = 'T';
    if(Gin->Block==true)
    {
      Gout->data = malloc(sizeof(BlockTriDiagonalMatrixType));
      if(BlockTriDiagonalMatrixEffectiveHamiltonian((BlockTriDiagonalMatrixType*)(Gin->data), (BlockTriDiagonalMatrixType*)(Gout->data)))
      {
        printf("BlockTriDiagonalMatrixEffectiveHamiltonian failed in ResponseFunctionEffectiveHamiltonian\n");
        fflush(stdout);
        free(Gout->data);
        Gout->data = NULL;
        FreeCompactMatrix(A0);
        FreeCompactMatrix(B0);
        return 1;
      }
    }
    else
    {
      Gout->data = malloc(sizeof(TriDiagonalMatrixType));
      if(TriDiagonalMatrixEffectiveHamiltonian((TriDiagonalMatrixType*)(Gin->data), (TriDiagonalMatrixType*)(Gout->data)))
      {
        printf("TriDiagonalMatrixEffectiveHamiltonian failed in ResponseFunctionEffectiveHamiltonian\n");
        fflush(stdout);
        free(Gout->data);
        Gout->data = NULL;
        FreeCompactMatrix(A0);
        FreeCompactMatrix(B0);
        return 1;
      }
    }
    return 0;
  }
  /*  everything else ends at an Anderson matrix, whose effective Hamiltonian is its bath
   *  as a list of poles. l, n and d get there the standard way */
  if(type!='a')
  {
    if(ConvertResponsefunction(Gin, &Ga, 'A', "Default", kt, sv, MinimalPoleWeigth, MinimalPoleDistance, false, false, false))
    {
      printf("ConvertResponsefunction to an Anderson matrix failed in ResponseFunctionEffectiveHamiltonian\n");
      fflush(stdout);
      FreeCompactMatrix(A0);
      FreeCompactMatrix(B0);
      return 1;
    }
    converted = true;
  }
  Gout->type = 'L';
  if(Gin->Block==true)
  {
    BlockAndersonMatrixType *Aptr = (BlockAndersonMatrixType*)(converted ? Ga.data : Gin->data);
    Gout->data = malloc(sizeof(BlockListOfPoles));
    /*  the rank one weights and then the merge, rather than the grouping that
     *  BlockAndersonMatrixInverseToListOfPoles does with a tolerance of its own: this way
     *  MinimalPoleDistance means here what it means everywhere else, a negative value
     *  asking BlockListOfPolesMergeDegeneratePoles for its standard value, and degenerate
     *  poles are merged across bath blocks and not only inside one */
    if(BlockAndersonMatrixInverseToBlockListOfPoles(Aptr, (BlockListOfPoles*)(Gout->data)))
    {
      printf("BlockAndersonMatrixInverseToBlockListOfPoles failed in ResponseFunctionEffectiveHamiltonian\n");
      fflush(stdout);
      free(Gout->data);
      Gout->data = NULL;
      if(converted) FreeResponsefunction(&Ga);
      FreeCompactMatrix(A0);
      FreeCompactMatrix(B0);
      return 1;
    }
    if(BlockListOfPolesMergeDegeneratePoles((BlockListOfPoles*)(Gout->data), MinimalPoleDistance, false))
    {
      printf("BlockListOfPolesMergeDegeneratePoles failed in ResponseFunctionEffectiveHamiltonian\n");
      fflush(stdout);
      FreeBlockListOfPoles((BlockListOfPoles*)(Gout->data));
      free(Gout->data);
      Gout->data = NULL;
      if(converted) FreeResponsefunction(&Ga);
      FreeCompactMatrix(A0);
      FreeCompactMatrix(B0);
      return 1;
    }
  }
  else
  {
    AndersonMatrixType *Aptr = (AndersonMatrixType*)(converted ? Ga.data : Gin->data);
    Gout->data = malloc(sizeof(ListOfPoles));
    if(AndersonMatrixEffectiveHamiltonian(Aptr, (ListOfPoles*)(Gout->data)))
    {
      printf("AndersonMatrixEffectiveHamiltonian failed in ResponseFunctionEffectiveHamiltonian\n");
      fflush(stdout);
      free(Gout->data);
      Gout->data = NULL;
      if(converted) FreeResponsefunction(&Ga);
      FreeCompactMatrix(A0);
      FreeCompactMatrix(B0);
      return 1;
    }
  }
  if(converted) FreeResponsefunction(&Ga);
  return 0;
}

/*  The inverse of ResponseFunctionEffectiveHamiltonian: given the effective Hamiltonian
 *  standing in a denominator, return the response function it describes,
 *
 *    G(w) = ( w - Heff(w) )^-1
 *
 *  Which type G comes back as is chosen per input type by what is numerically cheapest,
 *  or most numerically stable. (Future options can influence the method used and type
 *  returned):
 *
 *    t -> t   no numerical work at all, only a shift in the data structure
 *    l -> a   one square root per pole and nothing more
 *    a -> a   through l
 *    n -> a   through l
 *    d -> a   through l
 *
 *  The three that need a transformation take the standard route,
 *  ConvertResponsefunction or ConvertBlockMatrix, so that they behave exactly as ChangeType
 *  would and take the same options.
 *
 *  t -> t is a shift by one block: a[0] and b[0] of Heff become a[1] and b[1] of G, and G
 *  keeps a[0] zero and b[0] the identity. b[0] of Heff sits at the top level, where a
 *  coupling enters as b0^* ( ) b0^T, and lands one level down in G, where it enters as
 *  b ( ) b^+, so b[1] of G is the conjugate of b[0] of Heff.
 *
 *  Gout is not allocated on entry.
 */
int ResponseFunctionFromEffectiveHamiltonian(ResponsefunctionType * const Gin, ResponsefunctionType * const Gout, double const T, double const SingularValue, double const MinimalPoleWeigth, double const MinimalPoleDistance)
{
  ResponsefunctionType Gl;
  bool converted = false;
  char type = Gin->type|32;
  if((type!='t')&&(type!='a')&&(type!='l')&&(type!='n')&&(type!='d'))
  {
    printf("Error in ResponseFunctionFromEffectiveHamiltonian: unrecognized type %c\n",Gin->type);
    fflush(stdout);
    return 1;
  }
  /*  SingularValue and T are a magnitude and a temperature, so a negative one cannot be
   *  meant: Set standard values
   */
  double sv = (SingularValue < 0) ? 1000*DBL_EPSILON : SingularValue;
  double kt = (T < 0) ? sqrt(DBL_MIN) : T;
  Gout->Block = Gin->Block;
  Gout->NSuggestedGridPoints = 0;
  Gout->SuggestedGrid = NULL;
  if(type=='t')
  {
    Gout->type = 'T';
    if(Gin->Block==true)
    {
      Gout->data = malloc(sizeof(BlockTriDiagonalMatrixType));
      if(BlockTriDiagonalMatrixFromEffectiveHamiltonian((BlockTriDiagonalMatrixType*)(Gin->data), (BlockTriDiagonalMatrixType*)(Gout->data)))
      {
        printf("BlockTriDiagonalMatrixFromEffectiveHamiltonian failed in ResponseFunctionFromEffectiveHamiltonian\n");
        fflush(stdout);
        free(Gout->data);
        Gout->data = NULL;
        return 1;
      }
    }
    else
    {
      Gout->data = malloc(sizeof(TriDiagonalMatrixType));
      if(TriDiagonalMatrixFromEffectiveHamiltonian((TriDiagonalMatrixType*)(Gin->data), (TriDiagonalMatrixType*)(Gout->data)))
      {
        printf("TriDiagonalMatrixFromEffectiveHamiltonian failed in ResponseFunctionFromEffectiveHamiltonian\n");
        fflush(stdout);
        free(Gout->data);
        Gout->data = NULL;
        return 1;
      }
    }
    return 0;
  }
  /*  everything else starts from a list of poles, whose a0 becomes the impurity block of an
   *  Anderson matrix and whose poles become one bath block each */
  if(type!='l')
  {
    if(ConvertResponsefunction(Gin, &Gl, 'L', "Default", kt, sv, MinimalPoleWeigth, MinimalPoleDistance, false, false, false))
    {
      printf("ConvertResponsefunction to a list of poles failed in ResponseFunctionFromEffectiveHamiltonian\n");
      fflush(stdout);
      return 1;
    }
    converted = true;
  }
  Gout->type = 'A';
  if(Gin->Block==true)
  {
    BlockListOfPoles *Lptr = (BlockListOfPoles*)(converted ? Gl.data : Gin->data);
    Gout->data = malloc(sizeof(BlockAndersonMatrixType));
    if(BlockAndersonMatrixFromListOfPolesAsSelfEnergy(Lptr, (BlockAndersonMatrixType*)(Gout->data)))
    {
      printf("BlockAndersonMatrixFromListOfPolesAsSelfEnergy failed in ResponseFunctionFromEffectiveHamiltonian\n");
      fflush(stdout);
      free(Gout->data);
      Gout->data = NULL;
      if(converted) FreeResponsefunction(&Gl);
      return 1;
    }
  }
  else
  {
    ListOfPoles *Lptr = (ListOfPoles*)(converted ? Gl.data : Gin->data);
    Gout->data = malloc(sizeof(AndersonMatrixType));
    if(AndersonMatrixFromEffectiveHamiltonian(Lptr, (AndersonMatrixType*)(Gout->data)))
    {
      printf("AndersonMatrixFromEffectiveHamiltonian failed in ResponseFunctionFromEffectiveHamiltonian\n");
      fflush(stdout);
      free(Gout->data);
      Gout->data = NULL;
      if(converted) FreeResponsefunction(&Gl);
      return 1;
    }
  }
  if(converted) FreeResponsefunction(&Gl);
  return 0;
}


int ResponseFunctionAdd(ResponsefunctionType * const GA, ResponsefunctionType * const GB, ResponsefunctionType * const Gsum)
{
  ResponsefunctionType GBLP;
  ConvertResponsefunction(GA, Gsum, 'l', "Default", 0, 0, 0, 0, false, false, false);
  ConvertResponsefunction(GB, &GBLP, 'l', "Default", 0, 0, 0, 0, false, false, false);
  
  if(GA->Block)
  {
    AddBlockListOfPoles((BlockListOfPoles*)(Gsum->data), *(BlockListOfPoles*)(GBLP.data));
  }
  else
  {
    AddListOfPoles((ListOfPoles*)(Gsum->data), *(ListOfPoles*)(GBLP.data));
  }
  FreeResponsefunction(&GBLP);
  return 0;
}


int ResponseFunctionSubtract(ResponsefunctionType * const GA, ResponsefunctionType * const GB, ResponsefunctionType * const Gdif)
{
  int retval=0;
  ResponsefunctionType GBLP;
  ConvertResponsefunction(GA, Gdif,  'l', "Default", 0, 0, 0, 1000*DBL_EPSILON, false, false, false);
  ConvertResponsefunction(GB, &GBLP, 'l', "Default", 0, 0, 0, 1000*DBL_EPSILON, false, false, false);
  
  if(GA->Block)
  {
    retval=SubtractBlockListOfPoles((BlockListOfPoles*)(Gdif->data), (BlockListOfPoles*)(GBLP.data));
  }
  else
  {
    retval=SubtractListOfPoles((ListOfPoles*)(Gdif->data), (ListOfPoles*)(GBLP.data) );
  }
  FreeResponsefunction(&GBLP);
  return retval;
}


/*  Adds a self energy to all sites except for the impurity site
 *  G0 is converted to an Anderson Impurity Matrix
 *  Sigma to a lsit of poles type
 *
 */
int CalculateGBath(ResponsefunctionType  * const G0, ResponsefunctionType  * const Sigma, ResponsefunctionType * const G, double Epsilon)
{
 if(G0->Block==false)
 {
   ResponsefunctionType G0And, SigmaLP;
   AndersonMatrixType *G0AndData;
   ListOfPoles        *SigmaLPData;
   if(!((G0->type=='a')||(G0->type=='A')))
   {
     ConvertResponsefunction(G0, &G0And, 'A', "Default", 0, 0, 0, 0, false, false, false);
     G0AndData = (AndersonMatrixType*)G0And.data;
   }
   else
   {
     G0AndData = (AndersonMatrixType*)G0->data;
   }
   if(!((Sigma->type=='l')||(Sigma->type=='L')))
   {
     ConvertResponsefunction(Sigma, &SigmaLP, 'L', "Default", 0, 0, 0, 0, false, false, false);
     SigmaLPData = (ListOfPoles*)SigmaLP.data;
   }
   else
   {
     SigmaLPData = (ListOfPoles*)Sigma->data;
   }
   G->type='a';
   G->Block=false;
   G->NSuggestedGridPoints=0;
   G->SuggestedGrid=NULL;
   G->data = (void*)malloc(sizeof(AndersonMatrixType));
   GBathAndListAndRepresentation(*G0AndData, G0AndData->mu, *SigmaLPData, (AndersonMatrixType *)(G->data) );
   if(!((G0->type=='a')||(G0->type=='A')))
   {
     FreeResponsefunction(&G0And);
   }
   if(!((Sigma->type=='l')||(Sigma->type=='L')))
   {
     FreeResponsefunction(&SigmaLP);
   }
 }
 else
 {
   ResponsefunctionType G0And, SigmaLP;
   BlockAndersonMatrixType *G0AndData;
   BlockListOfPoles        *SigmaLPData;
   if(!((G0->type=='a')||(G0->type=='A')))
   {
     ConvertResponsefunction(G0, &G0And, 'A', "Default", 0, 0, 0, 0, false, false, false);
     G0AndData = (BlockAndersonMatrixType*)G0And.data;
   }
   else
   {
     G0AndData = (BlockAndersonMatrixType*)G0->data;
   }
   if(!((Sigma->type=='l')||(Sigma->type=='L')))
   {
     ConvertResponsefunction(Sigma, &SigmaLP, 'L', "Default", 0, 0, 0, 0, false, false, false);
     SigmaLPData = (BlockListOfPoles*)SigmaLP.data;
   }
   else
   {
     SigmaLPData = (BlockListOfPoles*)Sigma->data;
   }
   G->type='a';
   G->Block=true;
   G->NSuggestedGridPoints=0;
   G->SuggestedGrid=NULL;
   G->data = (void*)malloc(sizeof(BlockAndersonMatrixType));
   BlockGBathAndListAndRepresentation(G0AndData, G0AndData->mu, SigmaLPData, (BlockAndersonMatrixType *)(G->data), Epsilon);
   if(!((G0->type=='a')||(G0->type=='A')))
   {
     FreeResponsefunction(&G0And);
   }
   if(!((Sigma->type=='l')||(Sigma->type=='L')))
   {
     FreeResponsefunction(&SigmaLP);
   }
 }
 return 0;
}


int CalculateGBathOnGrid(ResponsefunctionType  * const G0, ResponsefunctionType  * const Sigma, CompactVectorType *const Grid, ResponsefunctionType * const G, double Epsilon)
{
 if(G0->Block==false)
 {
   ResponsefunctionType G0And, SigmaLP;
   AndersonMatrixType *G0AndData;
   ListOfPoles        *SigmaLPData;
   if(!((G0->type=='a')||(G0->type=='A')))
   {
     ConvertResponsefunction(G0, &G0And, 'A', "Default", 0, 0, 0, 0, false, false, false);
     G0AndData = (AndersonMatrixType*)G0And.data;
   }
   else
   {
     G0AndData = (AndersonMatrixType*)G0->data;
   }
   if(!((Sigma->type=='l')||(Sigma->type=='L')))
   {
     ConvertResponsefunction(Sigma, &SigmaLP, 'L', "Default", 0, 0, 0, 0, false, false, false);
     SigmaLPData = (ListOfPoles*)SigmaLP.data;
   }
   else
   {
     SigmaLPData = (ListOfPoles*)Sigma->data;
   }
   G->type='a';
   G->Block=false;
   G->NSuggestedGridPoints=0;
   G->SuggestedGrid=NULL;
   G->data = (void*)malloc(sizeof(AndersonMatrixType));
   GBathAndListAndRepresentationOnGrid(*G0AndData, G0AndData->mu, *SigmaLPData, Grid->VRe, Grid->Ni, (AndersonMatrixType *)(G->data) );
   if(!((G0->type=='a')||(G0->type=='A')))
   {
     FreeResponsefunction(&G0And);
   }
   if(!((Sigma->type=='l')||(Sigma->type=='L')))
   {
     FreeResponsefunction(&SigmaLP);
   }
 }
 else
 {
   ResponsefunctionType G0And, SigmaLP;
   BlockAndersonMatrixType *G0AndData;
   BlockListOfPoles        *SigmaLPData;
   if(!((G0->type=='a')||(G0->type=='A')))
   {
     ConvertResponsefunction(G0, &G0And, 'A', "Default", 0, 0, 0, 0, false, false, false);
     G0AndData = (BlockAndersonMatrixType*)G0And.data;
   }
   else
   {
     G0AndData = (BlockAndersonMatrixType*)G0->data;
   }
   if(!((Sigma->type=='l')||(Sigma->type=='L')))
   {
     ConvertResponsefunction(Sigma, &SigmaLP, 'L', "Default", 0, 0, 0, 0, false, false, false);
     SigmaLPData = (BlockListOfPoles*)SigmaLP.data;
   }
   else
   {
     SigmaLPData = (BlockListOfPoles*)Sigma->data;
   }
   G->type='a';
   G->Block=true;
   G->NSuggestedGridPoints=0;
   G->SuggestedGrid=NULL;
   G->data = (void*)malloc(sizeof(BlockAndersonMatrixType));
   BlockGBathAndListAndRepresentationOnGrid(G0AndData, G0AndData->mu, SigmaLPData, Grid->VRe, Grid->Ni, (BlockAndersonMatrixType *)(G->data) );
   if(!((G0->type=='a')||(G0->type=='A')))
   {
     FreeResponsefunction(&G0And);
   }
   if(!((Sigma->type=='l')||(Sigma->type=='L')))
   {
     FreeResponsefunction(&SigmaLP);
   }
 }
 return 0;
}

/*  Calculate GBath(w) = ( G^-1 - Simga(w) )^-1
 *
 *  Calcualte G^-1
 *  Convert this to a list of poles
 *  Subtract list of poles to list of poles, mergin negative weight ( https://doi.org/10.1103/PhysRevB.90.085102 )
 *  Invert list of poles to GBath in Anderson Matrix form.
 */

int   CalculateGBathFromGandSigmaDirectSubtract(ResponsefunctionType  * const G, ResponsefunctionType  * const Sigma, ResponsefunctionType  * const GBath)
{
  ResponsefunctionType Heff, HeffLP, SigmaLP;
  CompactMatrixType A0, B0;
  /*  no options here: -DBL_MAX is the unset marker, which for MinimalPoleDistance means
   *  the standard merge and for the others is read by nothing */
  if(ResponseFunctionEffectiveHamiltonian(G, &Heff, &A0, &B0, -DBL_MAX, -DBL_MAX, -DBL_MAX, -DBL_MAX))
  {
    printf("ResponseFunctionEffectiveHamiltonian failed in CalculateGBathFromGandSigmaDirectSubtract\n");
    fflush(stdout);
    return 1;
  }
  FreeCompactMatrix(&A0);
  FreeCompactMatrix(&B0);
  ConvertResponsefunction(&Heff, &HeffLP,  'L', "Default", 0, 0, 0, 0, false, false, false); // capital L as C has controle of pointer
  FreeResponsefunction(&Heff);
  ConvertResponsefunction(Sigma, &SigmaLP, 'L', "Default", 0, 0, 0, 0, false, false, false); // capital L as C has controle of pointer
  if(HeffLP.Block==true)
  {
    if(SigmaLP.Block==true)
    { // Block Block
      SubtractBlockListOfPoles((BlockListOfPoles*)(HeffLP.data), (BlockListOfPoles*)(SigmaLP.data));
    }
    else
    { // Block Single
      BlockListOfPoles SigmaLPBlock;
      CopyListOfPolesToBlockListOfPoles((ListOfPoles *)(SigmaLP.data), &SigmaLPBlock);
      SubtractBlockListOfPoles((BlockListOfPoles*)(HeffLP.data), &SigmaLPBlock);
      FreeBlockListOfPoles(&SigmaLPBlock);
    }
  }
  else
  {
    if(SigmaLP.Block==true)
    { // Single Block
      BlockListOfPoles *HeffLPBlock;
      HeffLPBlock = (BlockListOfPoles *)malloc(sizeof(BlockListOfPoles));
      CopyListOfPolesToBlockListOfPoles((ListOfPoles *)(HeffLP.data), HeffLPBlock);
      SubtractBlockListOfPoles(HeffLPBlock, (BlockListOfPoles*)(SigmaLP.data));
      FreeListOfPoles((ListOfPoles *)(HeffLP.data));
      free(HeffLP.data);
      HeffLP.data = HeffLPBlock;
      HeffLP.Block=true;
    }
    else
    { // Single Single
      SubtractListOfPoles((ListOfPoles *)(HeffLP.data), (ListOfPoles *)(SigmaLP.data));
    }
  }
  FreeResponsefunction(&SigmaLP);
  if(ResponseFunctionFromEffectiveHamiltonian(&HeffLP, GBath, -DBL_MAX, -DBL_MAX, -DBL_MAX, -DBL_MAX))
  {
    printf("ResponseFunctionFromEffectiveHamiltonian failed in CalculateGBathFromGandSigmaDirectSubtract\n");
    fflush(stdout);
    FreeResponsefunction(&HeffLP);
    FreeResponsefunction(&SigmaLP);
    return 1;
  }
  FreeResponsefunction(&HeffLP);
  return 0;
}


int   CalculateGBathFromGandSigmaEigensystem(ResponsefunctionType  * const G, ResponsefunctionType  * const Sigma, ResponsefunctionType  * const GBath, double Epsilon, CompactVectorType *const Grid, uint32_t Nv, uint32_t Nc, ResponsefunctionType  * const Gv, ResponsefunctionType  * const Gc)
{
  CompactMatrixType MGBath, MG, MSigma, funL, funR;
  Complex16 *val;
  ResponseFunctionToCompactMatrix(G, &MG, true, Epsilon);
  ResponseFunctionToCompactMatrix(Sigma, &MSigma, true, Epsilon);
  MGBath.Ni = MG.Ni + MSigma.Ni;
  MGBath.Nj = MG.Nj + MSigma.Nj;
  MGBath.QComplex = 1;
  InitCompactMatrix(&MGBath);
  CompactMatrixAddShifted(&MGBath, MG, 1, 0, 0, 0);
  CompactMatrixAddShifted(&MGBath, MSigma, 1, 0, MG.Ni, MG.Nj);
  SquareCompactMatrixMakeHermitian(&MGBath);
  if(Sigma->Block==false)
  {
    double a0, b0;
    ResponseFunctionGetA0(Sigma, &a0);
    ResponseFunctionGetB0(Sigma, &b0);
    MGBath.V[0].Re-=a0;
    // Note we add -I*b0 on the off diagonal making the matrix non hermitian. The coupling enters as I^2 = -1 which we need to get the subtraction
    MGBath.V[MG.Ni].Im =-b0;
    MGBath.V[MG.Ni * MGBath.Nj].Im =-b0;
  }
  else
  {
    CompactMatrixType A0, B0;
    BlockResponseFunctionGetA0(Sigma, &A0);
    BlockResponseFunctionGetB0(Sigma, &B0);
    SquareCompactMatrixMakeHermitian(&A0);
    SquareCompactMatrixMakeHermitian(&B0);
    CompactMatrixAddShifted(&MGBath, A0, -1, 0, 0, 0);
    // Note we add I*B0 on the off diagonal making the matrix non hermitian. The coupling enters as I^2 = -1 which we need to get the subtraction
    CompactMatrixAddShifted(&MGBath, B0,  0,-1, MG.Ni, 0);
    CompactMatrixAddShifted(&MGBath, B0,  0,-1, 0, MG.Nj);
    FreeCompactMatrix(&A0);
    FreeCompactMatrix(&B0);
  }
  val = (Complex16*)malloc(MGBath.Ni * sizeof(Complex16));
  NonHermitianCompactMatrixEigensystem(&MGBath, &funL, &funR, val);
  GBath->type='L';
  GBath->Block=G->Block;
  GBath->SuggestedGrid=NULL;
  GBath->NSuggestedGridPoints=0;
  if(GBath->Block==true)
  {
    uint32_t pos=0;
    int i,j,k;
    BlockListOfPoles *GBathLP;
    GBath->data=malloc(sizeof(BlockListOfPoles));
    GBathLP = (BlockListOfPoles*)GBath->data;
    GBathLP->BlockSize = BlockResponseFunctionGetBlockSize(G);
    GBathLP->N = MGBath.Ni;
    GBathLP->QComplex=1;
    InitBlockListOfPoles(GBathLP);
    GBathLP->mu=ResponseFunctionGetMu(G);
    for(k=0;k<GBathLP->N;k++)
    {
      GBathLP->V[pos].Re = val[k].Re;
      pos++;
      for(i=0;i<GBathLP->BlockSize;i++)
      {
        for(j=0;j<GBathLP->BlockSize;j++)
        {
          GBathLP->V[pos].Re = funL.V[k*funL.Nj+i].Re * funR.V[k*funR.Nj+j].Re + funL.V[k*funL.Nj+i].Im * funR.V[k*funR.Nj+j].Im;
          GBathLP->V[pos].Im = funL.V[k*funL.Nj+i].Re * funR.V[k*funR.Nj+j].Im - funL.V[k*funL.Nj+i].Im * funR.V[k*funR.Nj+j].Re;
          pos++;
        }
      }
    }
    BlockListOfPolesMergeDegeneratePoles(GBathLP, Epsilon, false);
    MakeBlockListOfPolesPhysical(GBathLP);
  }
  else
  {
    ListOfPoles *GBathLP;
    GBath->data=malloc(sizeof(ListOfPoles));
    GBathLP = (ListOfPoles*)GBath->data;
    GBathLP->N=MGBath.Ni;
    InitListOfPoles(GBathLP);
    GBathLP->mu=ResponseFunctionGetMu(G);
    for(int i=0;i<GBathLP->N;i++)
    {
      GBathLP->a[i] = val[i].Re;
      GBathLP->b[i] = funL.V[i*GBathLP->N].Re * funR.V[i*GBathLP->N].Re + funL.V[i*GBathLP->N].Im * funR.V[i*GBathLP->N].Im;
    }
    ReduceListOfPolesSpacing(GBathLP, Epsilon);
    MakeListOfPolesPhysical(GBathLP);
  }
  
/*
  PrintCompactMatrixMinimal(MGBath);
  PrintCompactMatrixMinimal(funL);
  PrintCompactMatrixMinimal(funR);
  Complex16 Weight;
  for(int i=0;i<MGBath.Ni;i++)
  {
    Weight.Re = funL.V[i * MGBath.Ni].Re * funR.V[i * MGBath.Ni].Re + funL.V[i * MGBath.Ni].Im * funR.V[i * MGBath.Ni].Im;
    Weight.Im = funL.V[i * MGBath.Ni].Re * funR.V[i * MGBath.Ni].Im - funL.V[i * MGBath.Ni].Im * funR.V[i * MGBath.Ni].Re;
    printf("val ( %17.10E %17.10E ) weight ( %17.10E %17.10E ) \n",val[i].Re,val[i].Im,Weight.Re,Weight.Im);
  }
*/
  free(val);
  FreeCompactMatrix(&MGBath);
  FreeCompactMatrix(&funL);
  FreeCompactMatrix(&funR);
  return 0;
}


int   CalculateGBathFromGandSigma(ResponsefunctionType  * const G, ResponsefunctionType  * const Sigma, ResponsefunctionType  * const GBath, double Epsilon, CompactVectorType *const Grid, uint32_t Nv, uint32_t Nc, ResponsefunctionType  * const Gv, ResponsefunctionType  * const Gc, char const * const Method, char const * const Type)
{
  char GBathType;
  int retval;
  if(Method==NULL)                               retval = CalculateGBathFromGandSigmaDirectSubtract(G, Sigma, GBath);
  else if (strncmp("Subtract",    Method, 8)==0) retval = CalculateGBathFromGandSigmaDirectSubtract(G, Sigma, GBath);
  else if (strncmp("Eigensystem", Method,11)==0) retval = CalculateGBathFromGandSigmaEigensystem(G, Sigma, GBath, Epsilon, Grid, Nv, Nc, Gv, Gc);
  else return 1;
  if(retval!=0) return retval;
  
  // check if GBath->type equals Type, if not convert to the requested type. if Type is Null accept all types
  if(Type==NULL)                     GBathType=GBath->type;
  else if (strcmp("and",Type) == 0 ) GBathType='a';
  else if (strcmp("And",Type) == 0 ) GBathType='a';
  else if (strcmp("LP", Type) == 0 ) GBathType='l';
  else if (strcmp("lp", Type) == 0 ) GBathType='l';
  else if (strcmp("Tri",Type) == 0 ) GBathType='t';
  else if (strcmp("tri",Type) == 0 ) GBathType='t';
  else if (strcmp("Nat",Type) == 0 ) GBathType='n';
  else if (strcmp("nat",Type) == 0 ) GBathType='n';
  else
  {
    printf("WARNING: Unrecognzied type of required bath green's function in CalculateGBathFromGandSigma\n");
    printf("Returned arbetrary type %c \n",GBath->type);
    GBathType=GBath->type;
  }
  // TODO: change type to requested type
  
  // TODO: check if we need to set Gv and or Gc
  
  return 0;
}



/*  calculate G = G0( omega - Sigma(omega) )
 *
 *  represent G0 and S as a list of poles then
 *   G = sum_i B_i 1/(w-a_i - S) is for each i a response function on Anderson form
 *     we can transform this for each i to a list of poles and sum list of poles
 */
int CalculateGfromG0andSigma(ResponsefunctionType  * const G0, ResponsefunctionType  * const Sigma, ResponsefunctionType * const G, double Epsilon)
{
  if(G0->Block==false)
  {
    ResponsefunctionType G0LP, SigmaLP;
    ListOfPoles          *SigmaLPData, *G0LPData;
    if(!((G0->type=='l')||(G0->type=='L')))
    {
      ConvertResponsefunction(G0, &G0LP, 'L', "Default", 0, 0, 0, 0, false, false, false);
      G0LPData = (ListOfPoles*)G0LP.data;
    }
    else
    {
      G0LPData = (ListOfPoles*)G0->data;
    }
    if(!((Sigma->type=='l')||(Sigma->type=='L')))
    {
      ConvertResponsefunction(Sigma, &SigmaLP, 'L', "Default", 0, 0, 0, 0, false, false, false);
      SigmaLPData = (ListOfPoles*)SigmaLP.data;
    }
    else
    {
      SigmaLPData = (ListOfPoles*)Sigma->data;
    }
    G->type='l';
    G->Block=false;
    G->NSuggestedGridPoints=0;
    G->SuggestedGrid=NULL;
    G->data = (void*)malloc(sizeof(ListOfPoles));
    // We find G0 as a list of poles in G0LPData
    // We find Sigma as a list of poles in SigmaLPData
    // We need to calculate G as a list of poles in GLPData
    GListListListRepresentation(G0LPData, SigmaLPData, (ListOfPoles*)G->data);
    if(!((G0->type=='l')||(G0->type=='L')))
    {
      FreeResponsefunction(&G0LP);
    }
    if(!((Sigma->type=='l')||(Sigma->type=='L')))
    {
      FreeResponsefunction(&SigmaLP);
    }
  }
  else
  {
    ResponsefunctionType G0LP, SigmaLP;
    BlockListOfPoles     *SigmaLPData, *G0LPData;
    if(!((G0->type=='l')||(G0->type=='L')))
    {
      ConvertResponsefunction(G0, &G0LP, 'L', "Default", 0, 0, 0, 0, false, false, false);
      G0LPData = (BlockListOfPoles*)G0LP.data;
    }
    else
    {
      G0LPData = (BlockListOfPoles*)G0->data;
    }
    if(!((Sigma->type=='l')||(Sigma->type=='L')))
    {
      ConvertResponsefunction(Sigma, &SigmaLP, 'L', "Default", 0, 0, 0, 0, false, false, false);
      SigmaLPData = (BlockListOfPoles*)SigmaLP.data;
    }
    else
    {
      SigmaLPData = (BlockListOfPoles*)Sigma->data;
    }
    G->type='l';
    G->Block=true;
    G->NSuggestedGridPoints=0;
    G->SuggestedGrid=NULL;
    G->data = (void*)malloc(sizeof(BlockListOfPoles));
    // We find G0 as a list of poles in G0LPData
    // We find Sigma as a list of poles in SigmaLPData
    // We need to calculate G as a list of poles in GLPData
    BlockGListListListRepresentation(G0LPData, SigmaLPData, (BlockListOfPoles*)G->data);
    if(!((G0->type=='l')||(G0->type=='L')))
    {
      FreeResponsefunction(&G0LP);
    }
    if(!((Sigma->type=='l')||(Sigma->type=='L')))
    {
      FreeResponsefunction(&SigmaLP);
    }
  }
  return 1;
}
int CalculateGfromG0andSigmaOnGrid(ResponsefunctionType  * const G0, ResponsefunctionType  * const Sigma, CompactVectorType *const Grid, ResponsefunctionType * const G, double Epsilon)
{
  if(G0->Block==false)
  {
    ResponsefunctionType G0LP, SigmaLP;
    ListOfPoles        *SigmaLPData, *G0LPData;
    if(!((G0->type=='l')||(G0->type=='L')))
    {
      ConvertResponsefunction(G0, &G0LP, 'L', "Default", 0, 0, 0, 0, false, false, false);
      G0LPData = (ListOfPoles*)G0LP.data;
    }
    else
    {
      G0LPData = (ListOfPoles*)G0->data;
    }
    if(!((Sigma->type=='l')||(Sigma->type=='L')))
    {
      ConvertResponsefunction(Sigma, &SigmaLP, 'L', "Default", 0, 0, 0, 0, false, false, false);
      SigmaLPData = (ListOfPoles*)SigmaLP.data;
    }
    else
    {
      SigmaLPData = (ListOfPoles*)Sigma->data;
    }
    G->type='l';
    G->Block=false;
    G->NSuggestedGridPoints=0;
    G->SuggestedGrid=NULL;
    G->data = (void*)malloc(sizeof(ListOfPoles));
    // We find G0 as a list of poles in G0LPData
    // We find Sigma as a list of poles in SigmaLPData
    // We need to calculate G as a list of poles in GLPData
    GListListListRepresentationOnGrid(G0LPData, SigmaLPData, Grid->VRe, Grid->Ni, (ListOfPoles*)G->data);
    if(!((G0->type=='l')||(G0->type=='L')))
    {
      FreeResponsefunction(&G0LP);
    }
    if(!((Sigma->type=='l')||(Sigma->type=='L')))
    {
      FreeResponsefunction(&SigmaLP);
    }
  }
  else
  {
    ResponsefunctionType G0LP, SigmaLP;
    BlockListOfPoles     *SigmaLPData, *G0LPData;
    if(!((G0->type=='l')||(G0->type=='L')))
    {
      ConvertResponsefunction(G0, &G0LP, 'L', "Default", 0, 0, 0, 0, false, false, false);
      G0LPData = (BlockListOfPoles*)G0LP.data;
    }
    else
    {
      G0LPData = (BlockListOfPoles*)G0->data;
    }
    if(!((Sigma->type=='l')||(Sigma->type=='L')))
    {
      ConvertResponsefunction(Sigma, &SigmaLP, 'L', "Default", 0, 0, 0, 0, false, false, false);
      SigmaLPData = (BlockListOfPoles*)SigmaLP.data;
    }
    else
    {
      SigmaLPData = (BlockListOfPoles*)Sigma->data;
    }
    G->type='l';
    G->Block=true;
    G->NSuggestedGridPoints=0;
    G->SuggestedGrid=NULL;
    G->data = (void*)malloc(sizeof(BlockListOfPoles));
    // We find G0 as a list of poles in G0LPData
    // We find Sigma as a list of poles in SigmaLPData
    // We need to calculate G as a list of poles in GLPData
    BlockGListListListRepresentationOnGrid(G0LPData, SigmaLPData, Grid->VRe, Grid->Ni, (BlockListOfPoles*)G->data);
    if(!((G0->type=='l')||(G0->type=='L')))
    {
      FreeResponsefunction(&G0LP);
    }
    if(!((Sigma->type=='l')||(Sigma->type=='L')))
    {
      FreeResponsefunction(&SigmaLP);
    }
  }
  return 1;
}



/*  GBath(w) = 1 / ( G_0(w+mu-Sigma(w))^{-1} + Sigma(w) )
 *    G_0(w+mu-Sigma(w))                     = 1 / (w - S - a1 - sum_i b_i^2 / (w - a_{i+1} - S)
 *    G_0(w+mu-Sigma(w))^{-1}                =     (w - S - a1 - sum_i b_i^2 / (w - a_{i+1} - S)
 *    G_0(w+mu-Sigma(w))^{-1} + Sigma        =     (w     - a1 - sum_i b_i^2 / (w - a_{i+1} - S)
 *  ( G_0(w+mu-Sigma(w))^{-1} + Sigma )^{-1} = 1 / (w     - a1 - sum_i b_i^2 / (w - a_{i+1} - S)
 */
uint32_t GBathAndTriAndRepresentation(AndersonMatrixType const G0, double const mu, TriDiagonalMatrixType const Sigma, AndersonMatrixType * const G)
{
  uint32_t i, j;
  ListOfPoles bSigmaPole;
  TriDiagonalMatrixType bSigma;
  G->N = (G0.N - 1) * (Sigma.N + 1) + 1;
  if (InitAndersonMatrix(G))
  {
    printf("InitAndersonMatrix failed in GBathAndListAndRepresentation\n");
    return 1;
  }
  bSigma.N = Sigma.N + 1;
  if (InitTriDiagonalMatrix(&bSigma))
  {
    printf("InitTriDiagonalMatrix failed in GBathAndListAndRepresentation\n");
    return 1;
  }
  bSigma.a[0] = 0;
  for (j = 0; j < Sigma.N; j++)
  {
    bSigma.a[j + 2] = Sigma.a[j + 1];
    bSigma.b[j + 1] = Sigma.b[j];
  }
  for (i = 1; i < G0.N; i++)
  {
    bSigma.b[0] = G0.b[i];
    bSigma.a[1] = G0.a[i + 1] + Sigma.a[0];
    if (TriDiagonalMatrixToListOfPoles(bSigma, &bSigmaPole))
    {
      printf("TriDiagonalMatrixToListOfPoles failed in GBathAndListAndRepresentation\n");
      return 1;
    }
    for (j = 0; j < bSigmaPole.N; j++)
    {
      G->a[j + (i - 1) * (Sigma.N + 1) + 2] = bSigmaPole.a[j];
      G->b[j + (i - 1) * (Sigma.N + 1) + 1] = sqrt(bSigmaPole.b[j]);
    }
    for (j = bSigmaPole.N; j < bSigma.N; j++)
    {
      G->a[j + (i - 1) * (Sigma.N + 1) + 2] = 0;
      G->b[j + (i - 1) * (Sigma.N + 1) + 1] = 0;
    }
    FreeListOfPoles(&bSigmaPole);
  }
  G->a[0] = G0.a[0];
  G->b[0] = G0.b[0];
  G->a[1] = G0.a[1];
  FreeTriDiagonalMatrix(&bSigma);
  return 0;
}

/*  GBath(w) = 1 / ( G_0(w-Sigma(w))^{-1} + Sigma(w) )
 *    G_0(w-Sigma(w))                     = 1 / (w - S - a1 - b1^2 / (w - S - a2 - b2^2 / (w - S - ...
 *    G_0(w-Sigma(w))^{-1}                =     (w - S - a1 - b1^2 / (w - S - a2 - b2^2 / (w - S - ...
 *    G_0(w-Sigma(w))^{-1} + Sigma        =     (w     - a1 - b1^2 / (w - S - a2 - b2^2 / (w - S - ...
 *  ( G_0(w-Sigma(w))^{-1} + Sigma )^{-1} = 1 / (w     - a1 - b1^2 / (w - S - a2 - b2^2 / (w - S - ...
 */
uint32_t GBathTriDiagonalRepresentation(TriDiagonalMatrixType const G0, TriDiagonalMatrixType const Sigma, TriDiagonalMatrixType * const GBath, uint32_t NTriDiag)
{
  uint32_t i, j;
  MatrixType T;
  SparseMatrixType GS;
  GS.Ni = (G0.N - 1) * (Sigma.N + 1) + 1;
  GS.Nj = (G0.N - 1) * (Sigma.N + 1) + 1;
  GS.QComplex = 0;
  GS.NBitsKey = 8 * ((uint32_t) (log((double)GS.Ni) / (8 * log(2.0))) + 1);
  if (InitSparseMatrix(&GS))
  {
    printf("InitSparseMatrix failed in GBathTriDiagonalRepresentation\n");
    fflush(stdout);
    return 1;
  };
  GS.a0 = G0.a[0];
  GS.b0 = G0.b[0];
  RealSparseMatrixAddElement(&GS, 0, 0, G0.a[1]);
  RealSparseMatrixAddElement(&GS, 0, 1, G0.b[1]);
  RealSparseMatrixAddElement(&GS, 1, 0, G0.b[1]);
  for (i = 0; i < G0.N - 1; i++)
  {
    RealSparseMatrixAddElement(&GS, 1 + (i * (Sigma.N + 1)), 1 + (i * (Sigma.N + 1)), G0.a[i + 2] - Sigma.a[0]);
    for (j = 1; j < Sigma.N + 1; j++)
    {
      RealSparseMatrixAddElement(&GS, j + 1 + (i * (Sigma.N + 1)), j + 1 + (i * (Sigma.N + 1)), Sigma.a[j]);
      RealSparseMatrixAddElement(&GS, j + (i * (Sigma.N + 1)), j + 1 + (i * (Sigma.N + 1)), Sigma.b[j - 1]);
      RealSparseMatrixAddElement(&GS, j + 1 + (i * (Sigma.N + 1)), j + (i * (Sigma.N + 1)), Sigma.b[j - 1]);
    }
  }
  for (i = 1; i < G0.N - 1; i++)
  {
    RealSparseMatrixAddElement(&GS, 1 + ((i - 1) * (Sigma.N + 1)), 1 + (i * (Sigma.N + 1)), G0.b[i + 1]);
    RealSparseMatrixAddElement(&GS, 1 + (i * (Sigma.N + 1)), 1 + ((i - 1) * (Sigma.N + 1)), G0.b[i + 1]);
  }
  if (SparseMatrixToTriDiagonalMatrix(GS, &T, GBath, NTriDiag))
  {
    printf("SparseMatrixToTriDiagonalMatrix failed in GBathTriDiagonalRepresentation\n");
    fflush(stdout);
    return 1;
  };
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
uint32_t GBathAndListAndRepresentation(AndersonMatrixType const G0, double const mu, ListOfPoles const Sigma, AndersonMatrixType * const G)
{
  if(G0.N==0)
  {
    AndersonMatrixCopy(G, G0);
    return 0;
  }
  uint32_t i, j;
  ListOfPoles bSigmaPole;
  AndersonMatrixType bSigma;
  G->N = (G0.N - 1) * (Sigma.N + 1) + 1;
  if (InitAndersonMatrix(G))
  {
    printf("InitAndersonMatrix failed in GBathAndListAndRepresentation\n");
    return 1;
  }
  strcpy(G->Name, "GBath");
  bSigma.N = Sigma.N + 1;
  if (InitAndersonMatrix(&bSigma))
  {
    printf("InitAndersonMatrix failed in GBathAndListAndRepresentation\n");
    return 1;
  }
  bSigma.a[0] = 0;
  for (i = 1; i < G0.N; i++)
  {
    bSigma.b[0] = G0.b[i];
    bSigma.a[1] = G0.a[i + 1] + Sigma.a0;
    for (j = 0; j < Sigma.N; j++)
    {
      bSigma.a[j + 2] = Sigma.a[j];
      bSigma.b[j + 1] = sqrt(Sigma.b[j]);
    }
    if (AndersonMatrixToListOfPoles(bSigma, &bSigmaPole))
    {
      printf("AndersonMatrixToListOfPoles failed in GBathAndListAndRepresentation\n");
      return 1;
    }
    for (j = 0; j < bSigmaPole.N; j++)
    {
      G->a[j + (i - 1) * (Sigma.N + 1) + 2] = bSigmaPole.a[j];
      G->b[j + (i - 1) * (Sigma.N + 1) + 1] = sqrt(bSigmaPole.b[j]);
    }
    for (j = bSigmaPole.N; j < bSigma.N; j++)
    {
      G->a[j + (i - 1) * (Sigma.N + 1) + 2] = 0;
      G->b[j + (i - 1) * (Sigma.N + 1) + 1] = 0;
    }
    FreeListOfPoles(&bSigmaPole);
  }
  G->a[0] = G0.a[0];
  G->b[0] = G0.b[0];
  G->a[1] = G0.a[1];
  G->mu = mu;
  FreeAndersonMatrix(&bSigma);
  OrderAndersonMatrix(G);
  return 0;
}
/*
 *
 */
uint32_t GBathAndListAndRepresentationOnGrid(AndersonMatrixType const G0, double const mu, ListOfPoles const Sigma, double *Grid, uint32_t NGridPoints, AndersonMatrixType * const G)
{
  if(G0.N==0)
  {
    AndersonMatrixCopy(G, G0);
    return 0;
  }
  if(NGridPoints==0)
  {
    printf("Zero gridpoints in GBathAndListAndRepresentationOnGrid. That does not make sence.\n");
    return 1;
  }
  uint32_t i, j;
  ListOfPoles bSigmaPole,Ghyb;
  AndersonMatrixType bSigma;
  // we take the hybridization part as a list of poles on a grid to which we add
  // once all added we bring it to an anderson form
  Ghyb.N=NGridPoints;
  InitListOfPoles(&Ghyb);
  for(i=0;i<Ghyb.N;i++)
  {
    Ghyb.a[i]=Grid[i];
    Ghyb.b[i]=0;
  }
  bSigma.N = Sigma.N + 1;
  if (InitAndersonMatrix(&bSigma))
  {
    printf("InitAndersonMatrix failed in GBathAndListAndRepresentation\n");
    return 1;
  }
  bSigma.a[0] = 0;
  for (i = 1; i < G0.N; i++)
  {
    bSigma.b[0] = G0.b[i];
    bSigma.a[1] = G0.a[i + 1] + Sigma.a0;
    for (j = 0; j < Sigma.N; j++)
    {
      bSigma.a[j + 2] =  Sigma.a[j];
      bSigma.b[j + 1] = sqrt(Sigma.b[j]);
    }
    if (AndersonMatrixToListOfPoles(bSigma, &bSigmaPole))
    {
      printf("AndersonMatrixToListOfPoles failed in GBathAndListAndRepresentation\n");
      return 1;
    }
    // add bSigmaPole to the Ghyb
    // both are sorted
    
    AddToListOfPolesOnGridOrdered(&Ghyb, &bSigmaPole, Grid[0], Grid[NGridPoints-1]);
    FreeListOfPoles(&bSigmaPole);
  }
  G->N = NGridPoints+1;
  if (InitAndersonMatrix(G))
  {
    printf("InitAndersonMatrix failed in GBathAndListAndRepresentation\n");
    return 1;
  }
  for(i=0;i<NGridPoints;i++)
  {
    G->a[i+2] = Ghyb.a[i];
    G->b[i+1] = sqrt(Ghyb.b[i]);
  }
  G->a[0] = G0.a[0];
  G->b[0] = G0.b[0];
  G->a[1] = G0.a[1];
  G->mu = mu;
  strcpy(G->Name, "GBath");
  FreeListOfPoles(&Ghyb);
  FreeAndersonMatrix(&bSigma);
  OrderAndersonMatrix(G);
  return 0;
}



/*  GBath(w) = 1 / ( G_0(w-Sigma(w))^{-1} + Sigma(w) )
 *    G_0(w-Sigma(w))                     = 1 / (w - S - A1 - sum_i B_i . ( 1 / (w - A_{i+1} - S) ) . B_i^+ )
 *    G_0(w-Sigma(w))^{-1}                =     (w - S - A1 - sum_i B_i . ( 1 / (w - A_{i+1} - S) ) . B_i^+ )
 *    G_0(w-Sigma(w))^{-1} + Sigma        =     (w     - A1 - sum_i B_i . ( 1 / (w - A_{i+1} - S) ) . B_i^+ )
 *  ( G_0(w-Sigma(w))^{-1} + Sigma )^{-1} = 1 / (w     - A1 - sum_i B_i . ( 1 / (w - A_{i+1} - S) ) . B_i^+ )
 *  ( G_0(w-Sigma(w))^{-1} + Sigma )^{-1} = 1 / (w     - A1 - sum_i B_i . ( 1 / (w - A_{i+1} - AS0 - sum_j BS_j^2 /(w-aS_j) ) . B_i^+ )
 *
 *  Epsilon determines if two poles in the same block have the same energy or not.
 */
uint32_t BlockGBathAndListAndRepresentation(BlockAndersonMatrixType const * const G0, double const mu, BlockListOfPoles const * const Sigma, BlockAndersonMatrixType * const G, double Epsilon)
{
  if(G0->N==0)
  {
    CopyBlockAndersonMatrix(*G0, G);
    return 0;
  }
  int i,j,k,l;
  CompactMatrixType SigmaInv,G0Bi,GBiSqr,GBi;
  BlockListOfPoles G0iWithSigma;
  // write the inverse of Sigma as a matrix
  BlockListOfPolesInversedToCompactMatrix(Sigma, &SigmaInv);
  if((G0->QComplex!=0)&&(Sigma->QComplex==0)) MakeCompactRealMatrixComplex(&SigmaInv);
  // we use G0Bi and GBiSqr as a dummie pointer Compact matrix to elements in G0
  G0Bi.Ni=G0->BlockSize;
  G0Bi.Nj=G0->BlockSize;
  G0Bi.QComplex=G0->QComplex;
  G0Bi.Name=SigmaInv.Name;
  GBiSqr.Ni=G0->BlockSize;
  GBiSqr.Nj=G0->BlockSize;
  GBiSqr.QComplex=SigmaInv.QComplex;
  GBiSqr.Name=SigmaInv.Name;
  // initialize G -- this can be very large
  G->BlockSize = G0->BlockSize;
  G->QComplex = SigmaInv.QComplex;
  G->N = 1 + (G0->N-1) * (Sigma->N+1) * G0->BlockSize;
  InitBlockAndersonMatrix(G);
  // set G.A0 and G.B0
  if(G->QComplex==0)
  {
    memcpy(G->aRe, G0->aRe, 2*Sqr(G0->BlockSize)*sizeof(double));
    memcpy(G->bRe, G0->bRe, Sqr(G0->BlockSize)*sizeof(double));
  }
  else if(G0->QComplex==0)
  {
    for(i=0;i<Sqr(G0->BlockSize);i++)
    {
      G->b[i].Re = G0->bRe[i];
      G->b[i].Im = 0;
    }
    for(i=0;i<2*Sqr(G0->BlockSize);i++)
    {
      G->a[i].Re = G0->aRe[i];
      G->a[i].Im = 0;
    }
  }
  else
  {
    memcpy(G->a, G0->a, 2*Sqr(G0->BlockSize)*sizeof(Complex16));
    memcpy(G->b, G0->b, Sqr(G0->BlockSize)*sizeof(Complex16));
  }
  G->mu=G0->mu;
  G->N=1;
  // loop over all blocks in G0 (i=1;i<G0.N;i++ )
  for(i=1;i<G0->N;i++)
  {
    // Set the left top block of SigmaInv equal to G0.A[i+1] + Sigma.A[0]
    if(SigmaInv.QComplex==0)
    {
      for(k=0;k<G0->BlockSize;k++)
      {
        for(l=0;l<G0->BlockSize;l++)
        {
          SigmaInv.VRe[l + k*SigmaInv.Nj] = G0->aRe[(i+1) * Sqr(G0->BlockSize) + l + k*(G0->BlockSize)] + Sigma->a0Re[l + k*(G0->BlockSize)];
        }
      }
    }
    else
    {
      if(Sigma->QComplex==0)
      {
        for(k=0;k<G0->BlockSize;k++)
        {
          for(l=0;l<G0->BlockSize;l++)
          {
            SigmaInv.V[l + k*SigmaInv.Nj].Re = G0->a[(i+1) * Sqr(G0->BlockSize) + l + k*(G0->BlockSize)].Re + Sigma->a0Re[l + k*(G0->BlockSize)];
            SigmaInv.V[l + k*SigmaInv.Nj].Im = G0->a[(i+1) * Sqr(G0->BlockSize) + l + k*(G0->BlockSize)].Im;
          }
        }
      }
      else if(G0->QComplex==0)
      {
        for(k=0;k<G0->BlockSize;k++)
        {
          for(l=0;l<G0->BlockSize;l++)
          {
            SigmaInv.V[l + k*SigmaInv.Nj].Re = G0->aRe[(i+1) * Sqr(G0->BlockSize) + l + k*(G0->BlockSize)] + Sigma->a0[l + k*(G0->BlockSize)].Re;
            SigmaInv.V[l + k*SigmaInv.Nj].Im =                                                               Sigma->a0[l + k*(G0->BlockSize)].Im;
          }
        }
      }
      else
      { // all complex
        for(k=0;k<G0->BlockSize;k++)
        {
          for(l=0;l<G0->BlockSize;l++)
          {
            SigmaInv.V[l + k*SigmaInv.Nj].Re = G0->a[(i+1) * Sqr(G0->BlockSize) + l + k*(G0->BlockSize)].Re + Sigma->a0[l + k*(G0->BlockSize)].Re;
            SigmaInv.V[l + k*SigmaInv.Nj].Im = G0->a[(i+1) * Sqr(G0->BlockSize) + l + k*(G0->BlockSize)].Im + Sigma->a0[l + k*(G0->BlockSize)].Im;
          }
        }
      }
    }
    // SigmaInv to BlockListOfPoles
    HermitianCompactMatrixToBlockListOfPoles(&SigmaInv, &G0iWithSigma, G0->BlockSize);
    // Merge poles with energy differences smaller than Epsilon
    BlockListOfPolesMergeDegeneratePoles(&G0iWithSigma, Epsilon, false);
    // rotate BlockListOpPoles by G0.B[i]
    if(G0->QComplex==0)
    {
      G0Bi.VRe = &(G0->bRe[i*Sqr(G0->BlockSize)]);
    }
    else
    {
      G0Bi.V = &(G0->b[i*Sqr(G0->BlockSize)]);
    }
    BlockListOfPolesTransformInlineNC(&G0iWithSigma, &G0Bi);
    // List of Poles to AndersonImpurity matrix part
    if(G->QComplex==0)
    {
      for(j=0;j<G0iWithSigma.N;j++)
      {
        GBiSqr.VRe = &(G0iWithSigma.VRe[1 + j * (Sqr(G0iWithSigma.BlockSize)+1)]);
        CompactMatrixSqrt(&GBiSqr, &GBi);
        memcpy(&(G->bRe[G->N * Sqr(G->BlockSize)]), GBi.VRe, Sqr(G->BlockSize)*sizeof(double));
        FreeCompactMatrix(&GBi);
        G->N++;
        for(k=0;k<G->BlockSize;k++)
        {
          G->aRe[G->N*Sqr(G->BlockSize) + k * (1+G->BlockSize)] = G0iWithSigma.VRe[j * (Sqr(G0iWithSigma.BlockSize)+1)];
        }
      }
    }
    else
    {
      for(j=0;j<G0iWithSigma.N;j++)
      {
        GBiSqr.V = &(G0iWithSigma.V[1 + j * (Sqr(G0iWithSigma.BlockSize)+1)]);
        CompactMatrixSqrt(&GBiSqr, &GBi);
        memcpy(&(G->b[G->N * Sqr(G->BlockSize)]), GBi.V, Sqr(G->BlockSize)*sizeof(double));
        FreeCompactMatrix(&GBi);
        G->N++;
        for(k=0;k<G->BlockSize;k++)
        {
          G->a[G->N*Sqr(G->BlockSize) + k * (1+G->BlockSize)] = G0iWithSigma.V[j * (Sqr(G0iWithSigma.BlockSize)+1)];
        }
      }
    }
    FreeBlockListOfPoles(&G0iWithSigma);
  }
  FreeCompactMatrix(&SigmaInv);
  return 0;
}
/*
 *
 */
uint32_t BlockGBathAndListAndRepresentationOnGrid(BlockAndersonMatrixType const * const G0, double const mu, BlockListOfPoles const * const Sigma, double const * const Grid, uint32_t NGridPoints, BlockAndersonMatrixType * const G)
{
  int i,k,l;
  BlockListOfPoles G0iWithSigma,GInv;
  CompactMatrixType G0Bi,GBiSqr,GBi;
  if(G0->N==0)
  {
    CopyBlockAndersonMatrix(*G0, G);
    return 0;
  }
  CompactMatrixType SigmaInv;
  // write the inverse of Sigma as a matrix
  BlockListOfPolesInversedToCompactMatrix(Sigma, &SigmaInv);
  if((G0->QComplex!=0)&&(Sigma->QComplex==0)) MakeCompactRealMatrixComplex(&SigmaInv);
  // we use G0Bi and GBiSqr as a dummie pointer Compact matrix to elements in G0 and G
  G0Bi.Ni=G0->BlockSize;
  G0Bi.Nj=G0->BlockSize;
  G0Bi.QComplex=G0->QComplex;
  G0Bi.Name=SigmaInv.Name;
  // initialize G
  G->BlockSize = G0->BlockSize;
  G->QComplex = SigmaInv.QComplex;
  G->N = 1 + NGridPoints;
  InitBlockAndersonMatrix(G);
  GBiSqr.Ni=G->BlockSize;
  GBiSqr.Nj=G->BlockSize;
  GBiSqr.QComplex=G->QComplex;
  GBiSqr.Name=G->Name;
  // initialize Ginv
  GInv.N = NGridPoints;
  GInv.QComplex = G->QComplex;
  GInv.BlockSize = G->BlockSize;
  InitBlockListOfPoles(&GInv);
  if(GInv.QComplex==0)
  {
    for(i=0;i<NGridPoints;i++)
    {
      GInv.VRe[i*(1+Sqr(GInv.BlockSize))] = Grid[i];
    }
  }
  else
  {
    for(i=0;i<NGridPoints;i++)
    {
      GInv.V[i*(1+Sqr(GInv.BlockSize))].Re = Grid[i];
    }
  }
  // set G.A0 and G.B0
  if(G->QComplex==0)
  {
    memcpy(G->aRe, G0->aRe, 2*Sqr(G0->BlockSize)*sizeof(double));
    memcpy(G->bRe, G0->bRe, Sqr(G0->BlockSize)*sizeof(double));
  }
  else if(G0->QComplex==0)
  {
    for(i=0;i<Sqr(G0->BlockSize);i++)
    {
      G->b[i].Re = G0->bRe[i];
      G->b[i].Im = 0;
    }
    for(i=0;i<2*Sqr(G0->BlockSize);i++)
    {
      G->a[i].Re = G0->aRe[i];
      G->a[i].Im = 0;
    }
  }
  else
  {
    memcpy(G->a, G0->a, 2*Sqr(G0->BlockSize)*sizeof(Complex16));
    memcpy(G->b, G0->b, Sqr(G0->BlockSize)*sizeof(Complex16));
  }
  G->mu=G0->mu;
  // loop over all blocks in G0 (i=1;i<G0.N;i++ )
  for(i=1;i<G0->N;i++)
  {
    // Set the left top block of SigmaInv equal to G0.A[i+1] + Sigma.A[0]
    if(SigmaInv.QComplex==0)
    {
      for(k=0;k<G0->BlockSize;k++)
      {
        for(l=0;l<G0->BlockSize;l++)
        {
          SigmaInv.VRe[l + k*SigmaInv.Nj] = G0->aRe[(i+1) * Sqr(G0->BlockSize) + l + k*(G0->BlockSize)] + Sigma->a0Re[l + k*(G0->BlockSize)];
        }
      }
    }
    else
    {
      if(Sigma->QComplex==0)
      {
        for(k=0;k<G0->BlockSize;k++)
        {
          for(l=0;l<G0->BlockSize;l++)
          {
            SigmaInv.V[l + k*SigmaInv.Nj].Re = G0->a[(i+1) * Sqr(G0->BlockSize) + l + k*(G0->BlockSize)].Re + Sigma->a0Re[l + k*(G0->BlockSize)];
            SigmaInv.V[l + k*SigmaInv.Nj].Im = G0->a[(i+1) * Sqr(G0->BlockSize) + l + k*(G0->BlockSize)].Im;
          }
        }
      }
      else if(G0->QComplex==0)
      {
        for(k=0;k<G0->BlockSize;k++)
        {
          for(l=0;l<G0->BlockSize;l++)
          {
            SigmaInv.V[l + k*SigmaInv.Nj].Re = G0->aRe[(i+1) * Sqr(G0->BlockSize) + l + k*(G0->BlockSize)] + Sigma->a0[l + k*(G0->BlockSize)].Re;
            SigmaInv.V[l + k*SigmaInv.Nj].Im =                                                               Sigma->a0[l + k*(G0->BlockSize)].Im;
          }
        }
      }
      else
      { // all complex
        for(k=0;k<G0->BlockSize;k++)
        {
          for(l=0;l<G0->BlockSize;l++)
          {
            SigmaInv.V[l + k*SigmaInv.Nj].Re = G0->a[(i+1) * Sqr(G0->BlockSize) + l + k*(G0->BlockSize)].Re + Sigma->a0[l + k*(G0->BlockSize)].Re;
            SigmaInv.V[l + k*SigmaInv.Nj].Im = G0->a[(i+1) * Sqr(G0->BlockSize) + l + k*(G0->BlockSize)].Im + Sigma->a0[l + k*(G0->BlockSize)].Im;
          }
        }
      }
    }
    // SigmaInv to BlockListOfPoles
    HermitianCompactMatrixToBlockListOfPoles(&SigmaInv, &G0iWithSigma, G0->BlockSize);
    // rotate BlockListOpPoles by G0.B[i]
    if(G0->QComplex==0)
    {
      G0Bi.VRe = &(G0->bRe[i*Sqr(G0->BlockSize)]);
    }
    else
    {
      G0Bi.V = &(G0->b[i*Sqr(G0->BlockSize)]);
    }
    BlockListOfPolesTransformInlineNC(&G0iWithSigma, &G0Bi);
    // Add ListofPoles on Grid
    AddBlockListOfPolesOnMeshOfA(&GInv, &G0iWithSigma);
    FreeBlockListOfPoles(&G0iWithSigma);
  }
  FreeCompactMatrix(&SigmaInv);
  // add GInv into G (GInv.A0 = 0 by construction)
  if(G->QComplex==0)
  {
    for(i=0;i<NGridPoints;i++)
    {
      GBiSqr.VRe = &(GInv.VRe[i * (1+Sqr(G->BlockSize)) + 1]);
      CompactMatrixSqrt(&GBiSqr, &GBi);
      for(k=0;k<G->BlockSize;k++) G->aRe[k*(1+G->BlockSize) + (i+2)*Sqr(G->BlockSize)] = GInv.VRe[i * (1+Sqr(G->BlockSize))];
      memcpy(&(G->bRe[(i+1)*Sqr(G->BlockSize)]), GBi.VRe, Sqr(G->BlockSize)*sizeof(double));
      FreeCompactMatrix(&GBi);
    }
  }
  else
  {
    for(i=0;i<NGridPoints;i++)
    {
      GBiSqr.V = &(GInv.V[i * (1+Sqr(G->BlockSize)) + 1]);
      CompactMatrixSqrt(&GBiSqr, &GBi);
      for(k=0;k<G->BlockSize;k++) G->a[k*(1+G->BlockSize) + (i+2)*Sqr(G->BlockSize)].Re = GInv.V[i * (1+Sqr(G->BlockSize))].Re;
      memcpy(&(G->b[(i+1)*Sqr(G->BlockSize)]), GBi.V, Sqr(G->BlockSize)*sizeof(double));
      FreeCompactMatrix(&GBi);
    }
  }
  return 0;
}




/*  GBath(w) = 1 / ( G_0(w-Sigma(w))^{-1} + Sigma(w) )
 *    G_0(w-Sigma(w))                     = 1 / (w - S - a1 - sum_i b_i^2 / (w - a_{i+1} - S)
 *    G_0(w-Sigma(w))^{-1}                =     (w - S - a1 - sum_i b_i^2 / (w - a_{i+1} - S)
 *    G_0(w-Sigma(w))^{-1} + Sigma        =     (w     - a1 - sum_i b_i^2 / (w - a_{i+1} - S)
 *  ( G_0(w-Sigma(w))^{-1} + Sigma )^{-1} = 1 / (w     - a1 - sum_i b_i^2 / (w - a_{i+1} - S)
 */
uint32_t GBathAndTriAndRepresentationParallel(AndersonMatrixType const G0, TriDiagonalMatrixType const Sigma, AndersonMatrixType * const G)
{
  int    i, j, failed = 0;
  G->N = (G0.N - 1) * (Sigma.N + 1) + 1;
  if (InitAndersonMatrix(G))
  {
    printf("InitAndersonMatrix (G->N=%i) failed in GBathAndListAndRepresentation\n", G->N);
    fflush(stdout);
    return 1;
  }
  G->a[0] = G0.a[0];
  G->b[0] = G0.b[0];
  G->a[1] = G0.a[1];
#pragma omp parallel private(i,j)
  {
    ListOfPoles bSigmaPole;
    TriDiagonalMatrixType bSigma;
    bSigma.N = Sigma.N + 1;
    if (InitTriDiagonalMatrix(&bSigma))
    {
      printf("InitTriDiagonalMatrix (bSigma.N=%i) failed in GBathAndListAndRepresentation in parallel region\n",
             bSigma.N);
      fflush(stdout);
      failed = 1;
    }
    if (!failed)
    {
      bSigma.a[0] = 0;
      for (j = 0; j < Sigma.N; j++)
      {
        bSigma.a[j + 2] = Sigma.a[j + 1];
        bSigma.b[j + 1] = Sigma.b[j];
      }
#pragma omp for
      for (i = 1; i < G0.N; i++)
      {
        bSigma.b[0] = G0.b[i];
        bSigma.a[1] = G0.a[i + 1] + Sigma.a[0];
        if (TriDiagonalMatrixToListOfPoles(bSigma, &bSigmaPole))
        {
          printf
            ("TriDiagonalMatrixToListOfPoles (bSigma.N=%i) failed in GBathAndListAndRepresentation in parallel region\n",
             bSigma.N);
          fflush(stdout);
          failed = 1;
        }
        if (!failed)
        {
          for (j = 0; j < bSigmaPole.N; j++)
          {
            G->a[j + (i - 1) * (Sigma.N + 1) + 2] = bSigmaPole.a[j];
            G->b[j + (i - 1) * (Sigma.N + 1) + 1] = sqrt(bSigmaPole.b[j]);
          }
          for (j = bSigmaPole.N; j < bSigma.N; j++)
          {
            G->a[j + (i - 1) * (Sigma.N + 1) + 2] = 0;
            G->b[j + (i - 1) * (Sigma.N + 1) + 1] = 0;
          }
          FreeListOfPoles(&bSigmaPole);
        }
      }
      FreeTriDiagonalMatrix(&bSigma);
    }
  }
  if (failed)
    return 1;
  return 0;
}

/*  GBath(w) = 1 / ( G_0(w-Sigma(w))^{-1} + Sigma(w) )
 *    G_0(w-Sigma(w))                     = 1 / (w - S - a1 - sum_i b_i^2 / (w - a_{i+1} - S)
 *    G_0(w-Sigma(w))^{-1}                =     (w - S - a1 - sum_i b_i^2 / (w - a_{i+1} - S)
 *    G_0(w-Sigma(w))^{-1} + Sigma        =     (w     - a1 - sum_i b_i^2 / (w - a_{i+1} - S)
 *  ( G_0(w-Sigma(w))^{-1} + Sigma )^{-1} = 1 / (w     - a1 - sum_i b_i^2 / (w - a_{i+1} - S)
 */
uint32_t OldGBathAndTriAndRepresentation(AndersonMatrixType const G0, TriDiagonalMatrixType const Sigma, AndersonMatrixType * const G)
{
  uint32_t i, j;
  ListOfPoles bSigmaPole;
  TriDiagonalMatrixType bSigma;
  G->N = (G0.N - 1) * (Sigma.N + 1) + 1;
  if (InitAndersonMatrix(G))
  {
    printf("InitAndersonMatrix failed in GBathAndListAndRepresentation\n");
    return 1;
  }
  bSigma.N = Sigma.N + 1;
  if (InitTriDiagonalMatrix(&bSigma))
  {
    printf("InitTriDiagonalMatrix failed in GBathAndListAndRepresentation\n");
    return 1;
  }
  bSigma.a[0] = 0;
  for (j = 0; j < Sigma.N; j++)
  {
    bSigma.a[j + 2] = Sigma.a[j + 1];
    bSigma.b[j + 1] = Sigma.b[j];
  }
  for (i = 1; i < G0.N; i++)
  {
    bSigma.b[0] = G0.b[i];
    bSigma.a[1] = G0.a[i + 1] + Sigma.a[0];
    if (TriDiagonalMatrixToListOfPoles(bSigma, &bSigmaPole))
    {
      printf("TriDiagonalMatrixToListOfPoles failed in GBathAndListAndRepresentation\n");
      return 1;
    }
    for (j = 0; j < bSigmaPole.N; j++)
    {
      G->a[j + (i - 1) * (Sigma.N + 1) + 2] = bSigmaPole.a[j];
      G->b[j + (i - 1) * (Sigma.N + 1) + 1] = sqrt(bSigmaPole.b[j]);
    }
    for (j = bSigmaPole.N; j < bSigma.N; j++)
    {
      G->a[j + (i - 1) * (Sigma.N + 1) + 2] = 0;
      G->b[j + (i - 1) * (Sigma.N + 1) + 1] = 0;
    }
    FreeListOfPoles(&bSigmaPole);
  }
  G->a[0] = G0.a[0];
  G->b[0] = G0.b[0];
  G->a[1] = G0.a[1];
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
  TriDiagonalMatrixType GclInv, SigmaInv, GBathNewInv;
  GclInv.N = Gcl.N - 1;
  GclInv.a = &(Gcl.a[1]);
  GclInv.b = &(Gcl.b[1]);
  strcpy(GclInv.Name, Gcl.Name);
  strncat(GclInv.Name, "^-1 (-1)", Max(0, (int)((int)MAXStrLength - (int)strlen(Gcl.Name))));
  SigmaInv.N = Sigma.N + 1;
  InitTriDiagonalMatrix(&SigmaInv);
  SigmaInv.a[0] = 0;
  SigmaInv.b[0] = 1;
  for (i = 0; i < Sigma.N; i++)
  {
    SigmaInv.a[i + 1] = Sigma.a[i];
    SigmaInv.b[i + 1] = Sigma.b[i];
  }
  SigmaInv.a[Sigma.N + 1] = Sigma.a[Sigma.N];

  DifferenceGreensFunctionTriDiagonalRepresentation(SigmaInv, GclInv, &GBathNewInv);
  GBathNew->N = GBathNewInv.N + 1;
  InitTriDiagonalMatrix(GBathNew);
  GBathNew->a[0] = 0;
  GBathNew->b[0] = 1;
  for (i = 1; i < GBathNew->N; i++)
  {
    GBathNew->a[i] = GBathNewInv.a[i - 1];
    GBathNew->b[i] = GBathNewInv.b[i - 1];
  }
  GBathNew->a[GBathNew->N] = GBathNewInv.a[GBathNew->N - 1];
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
  TriDiagonalMatrixType GlocInv, GBathNewInv;
  GlocInv.N = Glocal.N - 1;
  GlocInv.a = &(Glocal.a[1]);
  GlocInv.b = &(Glocal.b[1]);
  strcpy(GlocInv.Name, Glocal.Name);
  strncat(GlocInv.Name, "^-1 (-1)", Max(0, (int)((int)MAXStrLength - (int)strlen(Glocal.Name))));
  DifferenceGreensFunctionTriDiagonalRepresentation(GlocInv, Sigma, &GBathNewInv);
  GBathNew->N = GBathNewInv.N + 1;
  InitTriDiagonalMatrix(GBathNew);
  GBathNew->a[0] = 0;
  GBathNew->b[0] = 1;
  for (i = 1; i < GBathNew->N; i++)
  {
    GBathNew->a[i] = GBathNewInv.a[i - 1];
    GBathNew->b[i] = GBathNewInv.b[i - 1];
  }
  GBathNew->a[GBathNew->N] = GBathNewInv.a[GBathNew->N - 1];
  FreeTriDiagonalMatrix(&GBathNewInv);
  return 0;
}

/*  G(w) = G_0(w-Sigma(w))
 *       = Ga0 + sum_i Gb_i 1/(w - Ga_i - Sigma(w))
 *       = Ga0 + sum_i Gb_i 1/(w - Ga_i - Sa0 - sum_j Sb_j (1/(w-Sa_j) ) )
 */
int GListListListRepresentation(ListOfPoles * G0, ListOfPoles * Sigma, ListOfPoles *G)
{
  if(G0->N==0) // special case when G0 has no poles
  {
    ListOfPolesCopy(G, *G0);
    return 0;
  }
  if(Sigma->N==0) // specia case when Sigma has no poles, still need to shift the Hartree term
  {
    ListOfPolesCopy(G, *G0);
    for(int i=0; i<G->N;i++)
    {
      G->a[i]-=Sigma->a0;
    }
    return 0;
  }
  // general case. The number of poles in G is the product of the number of poles in G0 and (Sigma+1)
  G->N = (G0->N) * (Sigma->N + 1);
  if(InitListOfPoles(G))
  {
    printf("InitListOfPoles failed in GListListListRepresentation\n");
    return 1;
  }
  strcpy(G->Name, "G");
  // for each pole in G0 we create an anderson representation of G0 with sigma and rotate this to a list of pole representation
  ListOfPoles bSigmaPole;
  AndersonMatrixType bSigma;
  bSigma.N = Sigma->N + 1;
  if (InitAndersonMatrix(&bSigma))
  {
    printf("InitAndersonMatrix failed in GBathAndListAndRepresentation\n");
    return 1;
  }
  bSigma.a[0] = 0;
  bSigma.a[1] = 0; // Sigma->a0 + G0 pole energy
  bSigma.b[0] = 0; // sqrt of G0 pole weigth
  for (int i=0; i<Sigma->N; i++)
  {
    bSigma.a[i+2] = Sigma->a[i];
    bSigma.b[i+1] = sqrt(Sigma->b[i]);
  }
  for (int i = 0; i < G0->N; i++)
  {
    bSigma.b[0] = sqrt(G0->b[i]);
    bSigma.a[1] = G0->a[i] + Sigma->a0;
    if (AndersonMatrixToListOfPoles(bSigma, &bSigmaPole))
    {
      printf("AndersonMatrixToListOfPoles failed in GListListListRepresentation\n");
      return 1;
    }
    for (int j = 0; j < bSigmaPole.N; j++)
    {
      G->a[j + i * (Sigma->N + 1) ] = bSigmaPole.a[j];
      G->b[j + i * (Sigma->N + 1) ] = bSigmaPole.b[j];
    }
    for (int j = bSigmaPole.N; j < bSigma.N; j++)
    {
      G->a[j + i * (Sigma->N + 1) ] = 0;
      G->b[j + i * (Sigma->N + 1) ] = 0;
    }
    FreeListOfPoles(&bSigmaPole);
  }
  OrderListOfPoles(G);
  G->a0 = G0->a0;
  G->mu = G0->mu;
  FreeAndersonMatrix(&bSigma);
  return 0;
}


int GListListListRepresentationOnGrid(ListOfPoles * G0, ListOfPoles * Sigma, double *Grid, uint32_t NGridPoints,  ListOfPoles *G)
{
  if(G0->N==0) // special case when G0 has no poles
  {
    ListOfPolesCopy(G, *G0);
    return 0;
  }
  if(NGridPoints==0)
  {
    printf("Zero gridpoints in GListListListRepresentationOnGrid. That does not make sence.\n");
    return 1;
  }
  // general case. The number of poles in G is determined by the grid
  G->N = NGridPoints;
  if(InitListOfPoles(G))
  {
    printf("InitListOfPoles failed in GListListListRepresentationOnGrid\n");
    return 1;
  }
  strcpy(G->Name, "G");
  G->a0 = G0->a0;
  G->mu = G0->mu;
  for(int i=0;i<G->N;i++)
  {
    G->a[i]=Grid[i];
    G->b[i]=0;
  }
  // for each pole in G0 we create an anderson representation of G0 with sigma and rotate this to a list of pole representation
  ListOfPoles bSigmaPole;
  AndersonMatrixType bSigma;
  bSigma.N = Sigma->N + 1;
  if (InitAndersonMatrix(&bSigma))
  {
    printf("InitAndersonMatrix failed in GBathAndListAndRepresentation\n");
    return 1;
  }
  bSigma.a[0] = 0;
  bSigma.a[1] = 0; // Sigma->a0 + G0 pole energy
  bSigma.b[0] = 0; // sqrt of G0 pole weigth
  for (int i=0; i<Sigma->N; i++)
  {
    bSigma.a[i+2] = Sigma->a[i];
    bSigma.b[i+1] = sqrt(Sigma->b[i]);
  }
  for (int i = 0; i < G0->N; i++)
  {
    bSigma.b[0] = sqrt(G0->b[i]);
    bSigma.a[1] = G0->a[i] + Sigma->a0;
    if (AndersonMatrixToListOfPoles(bSigma, &bSigmaPole))
    {
      printf("AndersonMatrixToListOfPoles failed in GListListListRepresentation\n");
      return 1;
    }
    AddToListOfPolesOnGridOrdered(G, &bSigmaPole, Grid[0], Grid[NGridPoints-1]);
    FreeListOfPoles(&bSigmaPole);
  }
  FreeAndersonMatrix(&bSigma);
  return 0;
}

/*
 *  Warning this is not the general form of G.
 *  with G0 and S represented as block list of poles
 *  Blocksize of G0 and Sigma need to be the same, not tested
 */
int BlockGListListListRepresentation(BlockListOfPoles * G0, BlockListOfPoles * Sigma, BlockListOfPoles *G)
{
  if(G0->N==0) // special case when G0 has no poles
  {
    CopyBlockListOfPoles(*G0, G);
    return 0;
  }
  // general case. The maximal number of poles in G is the product of the number of poles in G0 and (Sigma+1) times the size of the block.
  G->N = (G0->N) * (Sigma->N + 1) * G0->BlockSize;
  G->QComplex = G0->QComplex | Sigma->QComplex;
  G->BlockSize = G0->BlockSize;
  if(InitBlockListOfPoles(G))
  {
    printf("InitBlockListOfPoles failed in GListListListRepresentation\n");
    return 1;
  }
  strcpy(G->Name, "G");
  G->N=0;
  // for each pole in G0 we create an anderson representation of G0 with sigma and rotate this to a list of pole representation
  BlockListOfPoles bSigmaPole;
  BlockAndersonMatrixType bSigma;
  BlockAndersonMatrixFromListOfPolesAsSelfEnergy(Sigma, &bSigma);
  if(G->QComplex==0)
  {
    CompactMatrixType B,BSqrt;
    B.Ni=G0->BlockSize;
    B.Nj=G0->BlockSize;
    B.QComplex=G0->QComplex;
    B.Name=(char*)G0->Name;
    for (int i = 0; i < G0->N; i++)
    {
      // set B0 of bSigma to sqrt of G0->B[i]
      B.VRe = &(G0->VRe[i*(Sqr(G0->BlockSize)+1)+1]);
      CompactMatrixSqrt(&B, &BSqrt);
      memcpy(bSigma.bRe, BSqrt.VRe, Sqr(G0->BlockSize)*sizeof(double));
      FreeCompactMatrix(&BSqrt);
      // set A1 of bSigma to A0 of Sigma + G0->a[i] I
      memcpy(&(bSigma.aRe[Sqr(G0->BlockSize)]), Sigma->a0Re, Sqr(G0->BlockSize)*sizeof(double));
      for(int j=0; j<G0->BlockSize; j++)
      {
        bSigma.aRe[Sqr(G0->BlockSize) + j + j * G0->BlockSize] += G0->VRe[i*(Sqr(G0->BlockSize)+1)];
      }
      // transform to list of poles
      if (BlockAndersonMatrixToBlockListOfPoles(&bSigma, &bSigmaPole, 1000*DBL_EPSILON, false))
      {
        printf("BlockAndersonMatrixToBlockListOfPoles failed in GListListListRepresentation\n");
        return 1;
      }
      // add poles from bSigmaPole to G
      memcpy(&(G->VRe[(G->N)*(Sqr(G->BlockSize)+1)]), bSigmaPole.VRe, (Sqr(G->BlockSize)+1)*bSigmaPole.N*sizeof(double));
      G->N+=bSigmaPole.N;
      FreeBlockListOfPoles(&bSigmaPole);
    }
    memcpy(G->a0, G0->a0, Sqr(G0->BlockSize)*sizeof(double));
  }
  else
  {
    if(bSigma.QComplex==0)
    {
      RealBlockAndersonMatrixToComplex(&bSigma);
    }
    if(G0->QComplex==0)
    { // G0 real Sigma complex
      CompactMatrixType B,BSqrt;
      B.Ni=G0->BlockSize;
      B.Nj=G0->BlockSize;
      B.QComplex=G0->QComplex;
      B.Name=(char*)G0->Name;
      for (int i = 0; i < G0->N; i++)
      {
        // set B0 of bSigma to sqrt of G0->B[i]
        B.VRe = &(G0->VRe[i*(Sqr(G0->BlockSize)+1)+1]);
        CompactMatrixSqrt(&B, &BSqrt);
        for(int j=0;j<Sqr(G0->BlockSize);j++)
        {
          bSigma.b[j].Re = BSqrt.VRe[j];
          bSigma.b[j].Im = 0;
        }
        FreeCompactMatrix(&BSqrt);
        // set A1 of bSigma to A0 of Sigma + G0->a[i] I
        memcpy(&(bSigma.a[Sqr(G0->BlockSize)]), Sigma->a0, Sqr(G0->BlockSize)*sizeof(Complex16));
        for(int j=0; j<G0->BlockSize; j++)
        {
          bSigma.a[Sqr(G0->BlockSize) + j + j * G0->BlockSize].Re += G0->VRe[i*(Sqr(G0->BlockSize)+1)];
        }
        // transform to list of poles
        if (BlockAndersonMatrixToBlockListOfPoles(&bSigma, &bSigmaPole, 1000*DBL_EPSILON, false))
        {
          printf("BlockAndersonMatrixToBlockListOfPoles failed in GListListListRepresentation\n");
          return 1;
        }
        // add poles from bSigmaPole to G
        memcpy(&(G->V[(G->N)*(Sqr(G->BlockSize)+1)]), bSigmaPole.V, (Sqr(G->BlockSize)+1)*bSigmaPole.N*sizeof(Complex16));
        G->N+=bSigmaPole.N;
      }
      for(int j=0; j<Sqr(G0->BlockSize); j++)
      {
        G->a0[j].Re = G0->a0Re[j];
        G->a0[j].Im = 0;
      }
    }
    else
    {
      CompactMatrixType B,BSqrt;
      B.Ni=G0->BlockSize;
      B.Nj=G0->BlockSize;
      B.QComplex=G0->QComplex;
      B.Name=(char*)G0->Name;
      for (int i = 0; i < G0->N; i++)
      {
        // set B0 of bSigma to sqrt of G0->B[i]
        B.V = &(G0->V[i*(Sqr(G0->BlockSize)+1)+1]);
        CompactMatrixSqrt(&B, &BSqrt);
        CompactMatrixConjugate(&BSqrt);
        memcpy(bSigma.b, BSqrt.V, Sqr(G0->BlockSize)*sizeof(Complex16));
        FreeCompactMatrix(&BSqrt);
        // set A1 of bSigma to A0 of Sigma + G0->a[i] I
        memcpy(&(bSigma.a[Sqr(G0->BlockSize)]), Sigma->a0, Sqr(G0->BlockSize)*sizeof(Complex16));
        for(int j=0; j<G0->BlockSize; j++)
        {
          bSigma.a[Sqr(G0->BlockSize) + j + j * G0->BlockSize].Re += G0->V[i*(Sqr(G0->BlockSize)+1)].Re;
          bSigma.a[Sqr(G0->BlockSize) + j + j * G0->BlockSize].Im += G0->V[i*(Sqr(G0->BlockSize)+1)].Im;
        }
        // transform to list of poles
        if (BlockAndersonMatrixToBlockListOfPoles(&bSigma, &bSigmaPole, 1000*DBL_EPSILON, false))
        {
          printf("BlockAndersonMatrixToBlockListOfPoles failed in GListListListRepresentation\n");
          return 1;
        }
        // add poles from bSigmaPole to G
        memcpy(&(G->V[(G->N)*(Sqr(G->BlockSize)+1)]), bSigmaPole.V, (Sqr(G->BlockSize)+1)*bSigmaPole.N*sizeof(Complex16));
        G->N+=bSigmaPole.N;
        FreeBlockListOfPoles(&bSigmaPole);
      }
      memcpy(G->a0, G0->a0, Sqr(G0->BlockSize)*sizeof(Complex16));
    }
  }
  OrderBlockListOfPoles(G);
  G->mu = G0->mu;
  FreeBlockAndersonMatrix(&bSigma);
  return 0;
}


int BlockGListListListRepresentationOnGrid(BlockListOfPoles * G0, BlockListOfPoles * Sigma, double *Grid, uint32_t NGridPoints,  BlockListOfPoles *G)
{
  if(G0->N==0) // special case when G0 has no poles
  {
    CopyBlockListOfPoles(*G0, G);
    return 0;
  }
  // general case. The number of poles in G is determined by the grid
  if(NGridPoints==0)
  {
    printf("Zero gridpoints in BlockGListListListRepresentationOnGrid. That does not make sence.\n");
    return 1;
  }
  // general case. The number of poles in G is determined by the grid
  G->N = NGridPoints;
  G->QComplex = G0->QComplex | Sigma->QComplex;
  G->BlockSize = G0->BlockSize;
  if(InitBlockListOfPoles(G))
  {
    printf("InitBlockListOfPoles failed in GListListListRepresentation\n");
    return 1;
  }
  strcpy(G->Name, "G");
  G->a0 = G0->a0;
  G->mu = G0->mu;
  // initialize G
  if(G->QComplex==0)
  {
    for(int i=0;i<G->N;i++)
    {
      G->VRe[i*(Sqr(G->BlockSize)+1)]=Grid[i];
    }
  }
  else
  {
    for(int i=0;i<G->N;i++)
    {
      G->V[i*(Sqr(G->BlockSize)+1)].Re=Grid[i];
    }
  }
  // for each pole in G0 we create an anderson representation of G0 with sigma and rotate this to a list of pole representation
  BlockListOfPoles bSigmaPole;
  BlockAndersonMatrixType bSigma;
  BlockAndersonMatrixFromListOfPolesAsSelfEnergy(Sigma, &bSigma);
  if(G->QComplex==0)
  {
    CompactMatrixType B,BSqrt;
    B.Ni=G0->BlockSize;
    B.Nj=G0->BlockSize;
    B.QComplex=G0->QComplex;
    B.Name=(char*)G0->Name;
    for (int i = 0; i < G0->N; i++)
    {
      // set B0 of bSigma to sqrt of G0->B[i]
      B.VRe = &(G0->VRe[i*(Sqr(G0->BlockSize)+1)+1]);
      CompactMatrixSqrt(&B, &BSqrt);
      memcpy(bSigma.bRe, BSqrt.VRe, Sqr(G0->BlockSize)*sizeof(double));
      FreeCompactMatrix(&BSqrt);
      // set A1 of bSigma to A0 of Sigma + G0->a[i] I
      memcpy(&(bSigma.aRe[Sqr(G0->BlockSize)]), Sigma->a0Re, Sqr(G0->BlockSize)*sizeof(double));
      for(int j=0; j<G0->BlockSize; j++)
      {
        bSigma.aRe[Sqr(G0->BlockSize) + j + j * G0->BlockSize] += G0->VRe[i*(Sqr(G0->BlockSize)+1)];
      }
      // transform to list of poles
      if (BlockAndersonMatrixToBlockListOfPoles(&bSigma, &bSigmaPole, 1000*DBL_EPSILON, false))
      {
        printf("BlockAndersonMatrixToBlockListOfPoles failed in GListListListRepresentation\n");
        return 1;
      }
      // add poles from bSigmaPole to G
      AddBlockListOfPolesOnMeshOfA(G, &bSigmaPole);
      FreeBlockListOfPoles(&bSigmaPole);
    }
    memcpy(G->a0, G0->a0, Sqr(G0->BlockSize)*sizeof(double));
  }
  else
  {
    if(bSigma.QComplex==0)
    {
      RealBlockAndersonMatrixToComplex(&bSigma);
    }
    if(G0->QComplex==0)
    { // G0 real Sigma complex
      CompactMatrixType B,BSqrt;
      B.Ni=G0->BlockSize;
      B.Nj=G0->BlockSize;
      B.QComplex=G0->QComplex;
      B.Name=(char*)G0->Name;
      for (int i = 0; i < G0->N; i++)
      {
        // set B0 of bSigma to sqrt of G0->B[i]
        B.VRe = &(G0->VRe[i*(Sqr(G0->BlockSize)+1)+1]);
        CompactMatrixSqrt(&B, &BSqrt);
        for(int j=0;j<Sqr(G0->BlockSize);j++)
        {
          bSigma.b[j].Re = BSqrt.VRe[j];
          bSigma.b[j].Im = 0;
        }
        FreeCompactMatrix(&BSqrt);
        // set A1 of bSigma to A0 of Sigma + G0->a[i] I
        memcpy(&(bSigma.a[Sqr(G0->BlockSize)]), Sigma->a0, Sqr(G0->BlockSize)*sizeof(Complex16));
        for(int j=0; j<G0->BlockSize; j++)
        {
          bSigma.a[Sqr(G0->BlockSize) + j + j * G0->BlockSize].Re += G0->VRe[i*(Sqr(G0->BlockSize)+1)];
        }
        // transform to list of poles
        if (BlockAndersonMatrixToBlockListOfPoles(&bSigma, &bSigmaPole,DBL_MIN,false))
        {
          printf("BlockAndersonMatrixToBlockListOfPoles failed in GListListListRepresentation\n");
          return 1;
        }
        // add poles from bSigmaPole to G
        AddBlockListOfPolesOnMeshOfA(G, &bSigmaPole);
        FreeBlockListOfPoles(&bSigmaPole);
      }
      for(int j=0; j<Sqr(G0->BlockSize); j++)
      {
        G->a0[j].Re = G0->a0Re[j];
        G->a0[j].Im = 0;
      }
    }
    else
    {
      CompactMatrixType B,BSqrt;
      B.Ni=G0->BlockSize;
      B.Nj=G0->BlockSize;
      B.QComplex=G0->QComplex;
      B.Name=(char*)G0->Name;
      for (int i = 0; i < G0->N; i++)
      {
        // set B0 of bSigma to sqrt of G0->B[i]
        B.V = &(G0->V[i*(Sqr(G0->BlockSize)+1)+1]);
        CompactMatrixSqrt(&B, &BSqrt);
        CompactMatrixConjugate(&BSqrt);
        memcpy(bSigma.b, BSqrt.V, Sqr(G0->BlockSize)*sizeof(Complex16));
        FreeCompactMatrix(&BSqrt);
        // set A1 of bSigma to A0 of Sigma + G0->a[i] I
        memcpy(&(bSigma.a[Sqr(G0->BlockSize)]), Sigma->a0, Sqr(G0->BlockSize)*sizeof(Complex16));
        for(int j=0; j<G0->BlockSize; j++)
        {
          bSigma.a[Sqr(G0->BlockSize) + j + j * G0->BlockSize].Re += G0->V[i*(Sqr(G0->BlockSize)+1)].Re;
          bSigma.a[Sqr(G0->BlockSize) + j + j * G0->BlockSize].Im += G0->V[i*(Sqr(G0->BlockSize)+1)].Im;
        }
        // transform to list of poles
        if (BlockAndersonMatrixToBlockListOfPoles(&bSigma, &bSigmaPole, 1000*DBL_EPSILON, false))
        {
          printf("BlockAndersonMatrixToBlockListOfPoles failed in GListListListRepresentation\n");
          return 1;
        }
        // add poles from bSigmaPole to G
        AddBlockListOfPolesOnMeshOfA(G, &bSigmaPole);
        FreeBlockListOfPoles(&bSigmaPole);
      }
      memcpy(G->a0, G0->a0, Sqr(G0->BlockSize)*sizeof(Complex16));
    }
  }
  G->mu = G0->mu;
  FreeBlockAndersonMatrix(&bSigma);
  return 0;
}








/*  This routine returns a tri-diagonal matrix representing a greensfunction as a continued fraction in Gnew
 *  
 */
uint32_t TriDiagonalMatrixDysonScalar(TriDiagonalMatrixType const G, TriDiagonalMatrixType * const Gnew,  MatrixType const Rho, uint32_t const NRho, TriDiagonalMatrixType * const Gp)
{
  SparseMatrixType M;
  MatrixType T;
  uint32_t i, j, k, l, NTriDiag, jstuck;
  M.Ni = G.N;
  for (i = 0; i < NRho; i++)
  {
    M.Ni += Gp[i].N * G.N;
  }
  M.Nj = M.Ni;
  M.QComplex = 0;
  M.NBitsKey = 16;
  InitSparseMatrix(&M);
  M.a0 = G.a[0];
  M.b0 = G.b[0];
  j = 0;
  for (i = 0; i < G.N - 1; i++)
  {
    RealSparseMatrixAddElement(&M, j, j, G.a[i + 1]);
    RealSparseMatrixAddElement(&M, j, j + 1, G.b[i + 1]);
    RealSparseMatrixAddElement(&M, j + 1, j, G.b[i + 1]);
    j++;
  }
  RealSparseMatrixAddElement(&M, j, j, G.a[G.N]);
  j++;
  jstuck = j;
  for (k = 0; k < NRho; k++)
  {
    for (l = 0; l < G.N; l++)
    {
      RealSparseMatrixAddElement(&M, l, j, Gp[k].b[0] * sqrt(fabs(Rho.VRe[l][k])));
      RealSparseMatrixAddElement(&M, j, l, Gp[k].b[0] * sqrt(fabs(Rho.VRe[l][k])));
      for (i = 0; i < Gp[k].N - 1; i++)
      {
        RealSparseMatrixAddElement(&M, j, j, Gp[k].a[i + 1]);
        RealSparseMatrixAddElement(&M, j, j + 1, Gp[k].b[i + 1]);
        RealSparseMatrixAddElement(&M, j + 1, j, Gp[k].b[i + 1]);
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

  NTriDiag = G.N;
  for (i = 0; i < NRho; i++)
  {
    NTriDiag = Max(NTriDiag, Gp[i].N);
  }
  if (SparseMatrixToTriDiagonalMatrix(M, &T, Gnew, NTriDiag))
  {
    printf("SparseMatrixToTriDiagonalMatrix failed in TriDiagonalMatrixDysonScalar\n");
    fflush(stdout);
    return 1;
  }
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
  uint32_t i, j, k, l;
  Htot.Ni = O.N;
  for (i = 0; i < N; i++)
  {
    Htot.Ni += g[i].N * V[i].Nj;
  }
  Htot.Nj = Htot.Ni;
  if (Htot.Ni < 16000)
  {
    Htot.NBitsKey = 16;
  }
  else
  {
    Htot.NBitsKey = 24;
  }
  Htot.QComplex = 0;
  Htot.a0 = O.a[0];
  Htot.b0 = O.b[0];
  for (i = 0; i < N; i++)
  {
    if (V[i].QComplex)
      Htot.QComplex = 1;
  }
  if (InitSparseMatrix(&Htot))
  {
    printf("Could not initiallize sparse matrix in TridiagonalizeImpurityCoupledToBathDyson\n");
    return 1;
  }
  SparseMatrixAddTridiagonalMatrix(&Htot, 0, 0, O);
  j = O.N;
  for (i = 0; i < N; i++)
  {
    SparseMatrixAddMatrixTransposed(&Htot, 0, j, V[i], DBL_EPSILON);
    SparseMatrixAddMatrix(&Htot, j, 0, V[i], DBL_EPSILON);
    j += g[i].N * V[i].Nj;
  }
  l = O.N;
  if (Htot.QComplex)
  {
    for (i = 0; i < N; i++)
    {
      for (j = 0; j < g[i].N; j++)
      {
        for (k = 0; k < V[i].Nj; k++)
        {
          ComplexSparseMatrixAddElement(&Htot, l, l, H.VRe[i][k] + g[i].a[j + 1], 0);
          if (j != g[i].N - 1)
          {
            ComplexSparseMatrixAddElement(&Htot, l, l + V[i].Nj, g[i].b[j + 1], 0);
            ComplexSparseMatrixAddElement(&Htot, l + V[i].Nj, l, g[i].b[j + 1], 0);
          }
          l++;
        }
      }
    }
  }
  else
  {
    for (i = 0; i < N; i++)
    {
      for (j = 0; j < g[i].N; j++)
      {
        for (k = 0; k < V[i].Nj; k++)
        {
          RealSparseMatrixAddElement(&Htot, l, l, H.VRe[i][k] + g[i].a[j + 1]);
          if (j != g[i].N - 1)
          {
            RealSparseMatrixAddElement(&Htot, l, l + V[i].Nj, g[i].b[j + 1]);
            RealSparseMatrixAddElement(&Htot, l + V[i].Nj, l, g[i].b[j + 1]);
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
  uint32_t i, ii, j, k, l, m, size1, size2, indexi, indexii, indexj;


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

  Htot.Ni = O.N;
  /* Determine dimension of Htot */
  j = 0;
  for (i = 0; i < N; i++)
  {
    Htot.Ni += g[i].N * V1[i].Nj;
    for (ii = 0; ii < i; ii++)
    {
      Htot.Ni += g[i].N * g[ii].N * V2[2 * j].Nj;
      j++;
    }
  }
  Htot.Nj = Htot.Ni;
  if (Htot.Ni < 16000)
  {
    Htot.NBitsKey = 16;
  }
  else
  {
    Htot.NBitsKey = 24;
  }
  Htot.QComplex = 0;
  Htot.a0 = O.a[0];
  Htot.b0 = O.b[0];
  j = 0;
  for (i = 0; i < N; i++)
  {
    if (V1[i].QComplex)
      Htot.QComplex = 1;
    for (ii = 0; ii < i; ii++)
    {
      if (V2[2 * j].QComplex)
        Htot.QComplex = 1;
      if (V2[2 * j + 1].QComplex)
        Htot.QComplex = 1;
      j++;
    }
  }
  if (InitSparseMatrix(&Htot))
  {
    printf("Could not initiallize sparse matrix in TridiagonalizeImpurityCoupledToBathDyson\n");
    return 1;
  }
  /* Add matrix O to Htot */
  SparseMatrixAddTridiagonalMatrix(&Htot, 0, 0, O);
  /* Add matrices V1 to Htot */
  j = O.N;
  for (i = 0; i < N; i++)
  {
    SparseMatrixAddMatrixTransposed(&Htot, 0, j, V1[i], DBL_EPSILON);
    SparseMatrixAddMatrix(&Htot, j, 0, V1[i], DBL_EPSILON);
    j += g[i].N * V1[i].Nj;
  }
  size1 = j;
  /* Add matrices H1 and g to Htot */
  l = O.N;
  if (Htot.QComplex)
  {
    for (i = 0; i < N; i++)
    {
      for (j = 0; j < g[i].N; j++)
      {
        for (k = 0; k < V1[i].Nj; k++)
        {
          ComplexSparseMatrixAddElement(&Htot, l, l, H1.VRe[i][k] + g[i].a[j + 1], 0);
          if (j != g[i].N - 1)
          {
            ComplexSparseMatrixAddElement(&Htot, l, l + V1[i].Nj, g[i].b[j + 1], 0);
            ComplexSparseMatrixAddElement(&Htot, l + V1[i].Nj, l, g[i].b[j + 1], 0);
          }
          l++;
        }
      }
    }
  }
  else
  {
    for (i = 0; i < N; i++)
    {
      for (j = 0; j < g[i].N; j++)
      {
        for (k = 0; k < V1[i].Nj; k++)
        {
          RealSparseMatrixAddElement(&Htot, l, l, H1.VRe[i][k] + g[i].a[j + 1]);
          if (j != g[i].N - 1)
          {
            RealSparseMatrixAddElement(&Htot, l, l + V1[i].Nj, g[i].b[j + 1]);
            RealSparseMatrixAddElement(&Htot, l + V1[i].Nj, l, g[i].b[j + 1]);
          }
          l++;
        }
      }
    }
  }
  /* Add matrices V2 to Htot */
  j = 0;
  indexi = O.N;
  indexj = 0;
  for (i = 0; i < N; i++)
  {
    indexii = O.N;
    for (ii = 0; ii < i; ii++)
    {
      size2 = 0;
      for (k = 0; k < g[i].N; k++)
      {
        SparseMatrixAddMatrixTransposed(&Htot, indexi + k * V2[2 * j].Ni, size1 + indexj + k * V2[2 * j].Nj, V2[2 * j],
                                        DBL_EPSILON);
        SparseMatrixAddMatrix(&Htot, size1 + indexj + k * V2[2 * j].Nj, indexi + k * V2[2 * j].Ni, V2[2 * j],
                              DBL_EPSILON);
        size2 += V2[2 * j].Ni;
      }
      for (k = 0; k < g[i].N; k++)
      {
        SparseMatrixAddMatrixTransposed(&Htot, indexii + k * V2[2 * j + 1].Ni, size1 + indexj + k * size2,
                                        V2[2 * j + 1], DBL_EPSILON);
        SparseMatrixAddMatrix(&Htot, size1 + indexj + k * size2, indexii + k * V2[2 * j + 1].Ni, V2[2 * j + 1],
                              DBL_EPSILON);
      }
      indexii += V1[ii].Nj * g[ii].N;
      indexj += V2[2 * j].Nj * g[i].N * g[ii].N;
      j++;
    }
    indexi += V1[i].Nj * g[i].N;
  }
  /* Add matrices H2 and g to Htot */
  j = 0;
  indexj = size1;
  if (Htot.QComplex)
  {

  }
  else
  {
    for (i = 0; i < N; i++)
    {
      for (ii = 0; ii < i; ii++)
      {
        for (l = 0; l < g[ii].N; l++)
        {
          for (m = 0; m < g[i].N; m++)
          {
            for (k = 0; k < V2[j].Nj; k++)
            {
              RealSparseMatrixAddElement(&Htot, indexj, indexj, H2.VRe[j][k] + g[i].a[m + 1] + g[ii].a[l + 1]);
              if (m != g[i].N - 1)
              {
                RealSparseMatrixAddElement(&Htot, indexj, indexj + V2[j].Nj, g[i].b[m + 1]);
                RealSparseMatrixAddElement(&Htot, indexj + V2[j].Nj, indexj, g[i].b[m + 1]);
              }
              if (l != g[ii].N - 1)
              {
                RealSparseMatrixAddElement(&Htot, indexj, indexj + V2[j].Nj * g[i].N, g[ii].b[l + 1]);
                RealSparseMatrixAddElement(&Htot, indexj + V2[j].Nj * g[i].N, indexj, g[ii].b[l + 1]);
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
uint32_t TridDiagonalMatrixDysonToMuchSimplified(TriDiagonalMatrixType const G1, MatrixType const *M12,
                                                 TriDiagonalMatrixType const *g2, MatrixType const *Emu,
                                                 int32_t const Ni, TriDiagonalMatrixType * const Gnew)
{
  return 0;
}




/*
 *  assume f.a[0]=0
 *         g.a[0]=0;
 *         f.b[0]=1;
 *         g.b[0]=1;
 */
uint32_t GreenFunctionConvoluteTriDiagonalRepresentation(TriDiagonalMatrixType * const fg,
                                                         TriDiagonalMatrixType const f, TriDiagonalMatrixType const g,
                                                         uint32_t const NTriDiag)
{
  uint32_t i, j, k;
  SparseMatrixType fgSparse;
  MatrixType T;
  fgSparse.NBitsKey = 16;
  fgSparse.QComplex = 0;
  fgSparse.Ni = f.N * g.N;
  fgSparse.Nj = f.N * g.N;
  if (InitSparseMatrix(&fgSparse))
  {
    printf("InitSparseMatrix failed in GreenFunctionConvoluteTriDiagonalRepresentation\n");
    return 1;
  };
  fgSparse.a0 = 0;
  fgSparse.b0 = 1;
  k = 0;
  for (i = 1; i < f.N; i++)
  {
    for (j = 1; j < g.N; j++)
    {
      RealSparseMatrixAddElement(&fgSparse, k, k, f.a[i] + g.a[j]);
      RealSparseMatrixAddElement(&fgSparse, k, k + 1, g.b[j]);
      RealSparseMatrixAddElement(&fgSparse, k + 1, k, g.b[j]);
      RealSparseMatrixAddElement(&fgSparse, k, k + g.N, f.b[i]);
      RealSparseMatrixAddElement(&fgSparse, k + g.N, k, f.b[i]);
      k++;
    }
    RealSparseMatrixAddElement(&fgSparse, k, k, f.a[i] + g.a[g.N]);
    RealSparseMatrixAddElement(&fgSparse, k, k + g.N, f.b[i]);
    RealSparseMatrixAddElement(&fgSparse, k + g.N, k, f.b[i]);
    k++;
  }
  for (j = 1; j < g.N; j++)
  {
    RealSparseMatrixAddElement(&fgSparse, k, k, f.a[f.N] + g.a[j]);
    RealSparseMatrixAddElement(&fgSparse, k, k + 1, g.b[j]);
    RealSparseMatrixAddElement(&fgSparse, k + 1, k, g.b[j]);
    k++;
  }
  RealSparseMatrixAddElement(&fgSparse, k, k, f.a[f.N] + g.a[g.N]);
  if (SparseMatrixToTriDiagonalMatrix(fgSparse, &T, fg, NTriDiag))
  {
    printf("SparseMatrixToTriDiagonalMatrix failed in GreenFunctionConvoluteTriDiagonalRepresentation\n");
    return 1;
  };
  FreeMatrix(&T);
  FreeSparseMatrix(&fgSparse);
  return 0;
}



/****************************************************************************/
/********************  Evaluate Matrix to Complex16      ********************/
/****************************************************************************/


Complex16 EvaluateTriDiagonalMatrix(TriDiagonalMatrixType const * A, const double omega, const double Gamma)
{
  Complex16 res;
  double denom;
  uint32_t j;
  if(A->N==0)
  {
    res.Re = A->a[0];
    res.Im=0;
    return res;
  }
  denom = (Sqr(Gamma * 0.5) + Sqr(omega - A->a[A->N]));
  res.Re = (omega - A->a[A->N]) / denom;
  res.Im = (-Gamma * 0.5) / denom;
  for (j = A->N - 1; j >= 1; j--)
  {
    denom = Sqr(Gamma * 0.5 - res.Im * Sqr(A->b[j])) + Sqr(omega - A->a[j] - res.Re * Sqr(A->b[j]));
    res.Re = (omega - A->a[j] - res.Re * Sqr(A->b[j])) / denom;
    res.Im = (res.Im * Sqr(A->b[j]) - Gamma * 0.5) / denom;
  }
  res.Re = A->a[0] + res.Re * Sqr(A->b[0]);
  res.Im = res.Im * Sqr(A->b[0]);
  return res;
}

Complex16 EvaluateTriDiagonalMatrix(void const * A, const double omega, const double Gamma)
{
	TriDiagonalMatrixType *ATypecast = (TriDiagonalMatrixType*) A;
	return EvaluateTriDiagonalMatrix(ATypecast, omega, Gamma);
}

Complex16 EvaluateListOfPoles(ListOfPoles const * A, const double omega, const double Gamma)
{
  Complex16 res;
  int i;
  double bdivdenom,resIm,resRe,omina;
  double const Gh = Gamma * 0.5;
  resIm=0;
  resRe=A->a0;
#pragma omp parallel for private(bdivdenom,omina) reduction(+ : resIm,resRe)
  for(i=0;i<A->N;i++)
  {
    omina = omega - A->a[i];
    bdivdenom = A->b[i] / ( Sqr(Gh) + Sqr(omina) );
    resIm -= Gh * bdivdenom;
    resRe += omina * bdivdenom;
  }
  res.Im=resIm;
  res.Re=resRe;
  return res;
}

Complex16 EvaluateListOfPoles(void const * A, const double omega, const double Gamma)
{
	ListOfPoles *ATypecast = (ListOfPoles*) A;
	return EvaluateListOfPoles(ATypecast, omega, Gamma);
}

Complex16 EvaluateAndersonMatrix(AndersonMatrixType const * A, const double omega, const double Gamma)
{
  Complex16 res;
  uint32_t j;
  double denom;
  res.Im=0;
  res.Re=0;
  if(A->N==0)
  { /* no observed site and no bath, so G(w) = a[0]. There is no a[1] and no b[0] to read */
    res.Re = A->a[0];
    res.Im = 0;
    return res;
  }
  for (j = A->N; j >= 2; j--)
  {
    denom = Sqr(Gamma * 0.5) + Sqr(omega - A->a[j]);
    res.Re += Sqr(A->b[j - 1]) * (omega - A->a[j]) / denom;
    res.Im += Sqr(A->b[j - 1]) * (-Gamma * 0.5) / denom;
  }
  denom = (Sqr(Gamma * 0.5 - res.Im) + Sqr(omega - A->a[1] - res.Re));
  res.Re = A->a[0] + (omega - A->a[1] - res.Re) * Sqr(A->b[0]) / denom;
  res.Im = (res.Im - Gamma * 0.5) * Sqr(A->b[0]) / denom;
  return res;
}

Complex16 EvaluateAndersonMatrix(void const * A, const double omega, const double Gamma)
{
	AndersonMatrixType *ATypecast = (AndersonMatrixType*) A;
	return EvaluateAndersonMatrix(ATypecast, omega, Gamma);
}

Complex16 EvaluateBADoubleTriDiagonalMatrix(BADoubleTriDiagonalMatrixType const * A, const double omega, const double Gamma)
{
  Complex16 res,rescon;
  rescon = EvaluateTriDiagonalMatrix(&(A->con), omega, Gamma);
  res    = EvaluateTriDiagonalMatrix(&(A->val), omega, Gamma);
  res.Re += rescon.Re;
  res.Im += rescon.Im;
  res.Re *= Sqr(A->b0);
  res.Im *= Sqr(A->b0);
  res.Re += A->a0;
  return res;
}

Complex16 EvaluateBADoubleTriDiagonalMatrix(void const * A, const double omega, const double Gamma)
{
	BADoubleTriDiagonalMatrixType *Atypecast = (BADoubleTriDiagonalMatrixType*) A;
	return EvaluateBADoubleTriDiagonalMatrix(Atypecast, omega, Gamma);
}

Complex16 EvaluateDoubleTriDiagonalMatrix(DoubleTriDiagonalMatrixType const * A, const double omega, const double Gamma)
{
  Complex16 res;
  res.Re = omega;
  res.Im = Gamma/2;
  res -= EvaluateTriDiagonalMatrix(&(A->val),omega,Gamma);
  res -= EvaluateTriDiagonalMatrix(&(A->con),omega,Gamma);
  res = A->a0 + A->b0 * A->b0 / res;
  return res;
}

Complex16 EvaluateDoubleTriDiagonalMatrix(void const * A, const double omega, const double Gamma)
{
  DoubleTriDiagonalMatrixType *ATypecast = (DoubleTriDiagonalMatrixType*) A;
  return EvaluateDoubleTriDiagonalMatrix(ATypecast, omega, Gamma);
}

Complex16 EvaluateResponsefunctionType(ResponsefunctionType const *A, const double omega, const double Gamma)
{
  const Complex16 res={0.0,0.0};
  switch (A->type)
  {
    case 't':
    case 'T':
    {
      return EvaluateTriDiagonalMatrix(((TriDiagonalMatrixType*)(A->data)), omega, Gamma);
    }
    case 'a':
    case 'A':
    {
      return EvaluateAndersonMatrix(A->data, omega, Gamma);
    }
    case 'l':
    case 'L':
    {
      return EvaluateListOfPoles(A->data, omega, Gamma);
    }
    case 'n':
    case 'N':
    {
      return EvaluateBADoubleTriDiagonalMatrix(A->data, omega, Gamma);
    }
    case 'd':
    case 'D':
    {
      return EvaluateDoubleTriDiagonalMatrix(A->data, omega, Gamma);
    }
    default:
    {
      printf("Error in EvaluateResponseFunction: Type of Response function should be either T/t,A/a,L/l, N/n, or D/d and not %c\n",A->type);
      return res;
    }
  }
  return res;
}


/*
 *
 */
int ConvertResponsefunction(const ResponsefunctionType *const Min, ResponsefunctionType *Mout, const char type, const char * const Method, const double T, const double SingularValue, const double MinimalPoleWeigth, const double MinimalPoleDistance, bool deflate, bool relativeDistance, bool relativeWeigth)
{
  if(Min->Block==true)
  {
    return ConvertBlockMatrix(Min,Mout,type,Method, T,SingularValue, MinimalPoleWeigth, MinimalPoleDistance, deflate, relativeDistance, relativeWeigth);
  }
  Mout->NSuggestedGridPoints = 0;
  Mout->type = type;
  Mout->Block = false;
  /*  every implemented cell below allocates Mout->data. A cell that is not
   *  implemented only prints, so data staying NULL is what marks the failure */
  Mout->data = NULL;
  Mout->NSuggestedGridPoints=Min->NSuggestedGridPoints;
  if(Mout->NSuggestedGridPoints==0)
  {
    Mout->SuggestedGrid=NULL;
  }
  else
  {
    Mout->SuggestedGrid= (double*)malloc(Mout->NSuggestedGridPoints*sizeof(double));
    memcpy(Mout->SuggestedGrid, Min->SuggestedGrid, Mout->NSuggestedGridPoints*sizeof(double));
  }
  switch (Min->type|32)      // original type
  {
    case 'l': // original list of poles
    {
      switch (type|32)     // target type
      {
        case 'l':
        {
          Mout->data = (void*) malloc(sizeof(ListOfPoles));
          ListOfPolesCopy((ListOfPoles*)Mout->data, *((ListOfPoles*)Min->data ));
        }
        break;
        case 't':
        {
          Mout->data = (void*) malloc(sizeof(TriDiagonalMatrixType));
          ListOfPolesToTridiagonalMatrix(*(ListOfPoles*)Min->data, (TriDiagonalMatrixType*)Mout->data, ((BlockListOfPoles*)(Min->data))->N);
        }
        break;
        case 'a':
        {
          Mout->data = (void*) malloc(sizeof(AndersonMatrixType));
          ListOfPolesToAndersonMatrix(*(ListOfPoles*)Min->data, (AndersonMatrixType*)Mout->data);
        }
        break;
        case 'n':
        {
          Mout->data = (void*) malloc(sizeof(BADoubleTriDiagonalMatrixType));
          ListOfPolesToBADoubleTridiagonalMatrix(*(ListOfPoles*)Min->data, (BADoubleTriDiagonalMatrixType*)Mout->data, ((BlockListOfPoles*)(Min->data))->N, ((BlockListOfPoles*)(Min->data))->N, ((BlockListOfPoles*)(Min->data))->mu, T);
        }
        break;
        case 'd':
        {
          Mout->data = (void*) malloc(sizeof(DoubleTriDiagonalMatrixType));
          ListOfPolesToDoubleTridiagonalMatrix(*(ListOfPoles*)Min->data, (DoubleTriDiagonalMatrixType*)Mout->data, ((ListOfPoles*)(Min->data))->mu, T);
        }
        break;
        default:
        {
          printf("Error in ConvertResponsefunction unrecognized type to convert to: %c\n", type);
        }
        break;
      }
    }
    break;
    case 't': // original tri-diagonal
    {
      switch (type|32)
      {
        case 'l':
        {
          Mout->data = (void*) malloc(sizeof(ListOfPoles));
          TriDiagonalMatrixToListOfPoles(*((TriDiagonalMatrixType*)Min->data), (ListOfPoles*)Mout->data);
        }
        break;
        case 't':
        {
          Mout->data = (void*) malloc(sizeof(TriDiagonalMatrixType));
          TriDiagonalMatrixCopy((TriDiagonalMatrixType*)Mout->data, *(TriDiagonalMatrixType*)Min->data);
        }
        break;
        case 'a':
        {
          Mout->data = (void*) malloc(sizeof(AndersonMatrixType));
          TriDiagonalMatrixToAndersonMatrix(*(TriDiagonalMatrixType*)Min->data, (AndersonMatrixType*)Mout->data);
        }
        break;
        case 'n':
        {
          Mout->data = (void*) malloc(sizeof(BADoubleTriDiagonalMatrixType));
          TriDiagonalMatrixToBADoubleTriDiagonalMatrix(*(TriDiagonalMatrixType*)Min->data, (BADoubleTriDiagonalMatrixType*)Mout->data, ((TriDiagonalMatrixType*)Min->data)->N, ((TriDiagonalMatrixType*)Min->data)->N, ((TriDiagonalMatrixType*)Min->data)->mu, T);
        }
        break;
        case 'd':
        {
          Mout->data = (void*) malloc(sizeof(DoubleTriDiagonalMatrixType));
          TriDiagonalMatrixToDoubleTriDiagonalMatrix(*(TriDiagonalMatrixType*)Min->data, (DoubleTriDiagonalMatrixType*)Mout->data, ((TriDiagonalMatrixType*)Min->data)->mu, T);
        }
        break;
        default:
        {
          printf("Error in ConvertResponsefunction unrecognized type to convert to: %c\n",type);
        }
        break;
      }
    }
    break;
    case 'a': // original anderson matrix
    {
      switch (type|32)
      {
        case 'l':
        {
          Mout->data = (void*) malloc(sizeof(ListOfPoles));
          AndersonMatrixToListOfPoles(*((AndersonMatrixType*)(Min->data)), (ListOfPoles*)Mout->data);
        }
        break;
        case 't':
        {
          Mout->data = (void*) malloc(sizeof(TriDiagonalMatrixType));
          AndersonMatrixToTriDiagonalMatrix(*(AndersonMatrixType*)Min->data, (TriDiagonalMatrixType*)Mout->data, ((AndersonMatrixType*)Min->data)->N);
        }
        break;
        case 'a':
        {
          Mout->data = (void*) malloc(sizeof(AndersonMatrixType));
          AndersonMatrixCopy((AndersonMatrixType*)Mout->data,*(AndersonMatrixType*)Min->data);
        }
        break;
        case 'n':
        {
          Mout->data = (void*) malloc(sizeof(BADoubleTriDiagonalMatrixType));
          AndersonMatrixToBADoubleTriDiagonalMatrix(*(AndersonMatrixType*)Min->data, (BADoubleTriDiagonalMatrixType*)Mout->data, ((AndersonMatrixType*)Min->data)->N, ((AndersonMatrixType*)Min->data)->N, ((AndersonMatrixType*)Min->data)->mu, T);
        }
        break;
        case 'd':
        {
          Mout->data = (void*) malloc(sizeof(DoubleTriDiagonalMatrixType));
          AndersonMatrixToDoubleTriDiagonalMatrix(*(AndersonMatrixType*)Min->data, (DoubleTriDiagonalMatrixType*)Mout->data, ((AndersonMatrixType*)Min->data)->mu, T);
        }
        break;
        default:
        {
          printf("Error in ConvertResponsefunction unrecognized type of M: %c",Min->type);
        }
        break;
      }
    }
    break;
    case 'n': // original natural impurity orbitals (BADoubleTriDiagonalMatrixType)
    {
      switch (type|32)
      {
        case 'l':
        {
          Mout->data = (void*) malloc(sizeof(ListOfPoles));
          BADoubleTriDiagonalMatrixToListOfPoles((BADoubleTriDiagonalMatrixType*)Min->data, (ListOfPoles*)Mout->data);
        }
        break;
        case 't':
        {
          Mout->data = (void*) malloc(sizeof(TriDiagonalMatrixType));
          BADoubleTriDiagonalMatrixToTriDiagonalMatrix((BADoubleTriDiagonalMatrixType*)Min->data, (TriDiagonalMatrixType*)Mout->data);
        }
        break;
        case 'a':
        {
          Mout->data = (void*) malloc(sizeof(AndersonMatrixType));
          BADoubleTriDiagonalMatrixToAndersonMatrix((BADoubleTriDiagonalMatrixType*)Min->data, (AndersonMatrixType*)Mout->data);
        }
        break;
        case 'n':
        {
          Mout->data = (void*) malloc(sizeof(BADoubleTriDiagonalMatrixType));
          BADoubleTriDiagonalMatrixCopy((BADoubleTriDiagonalMatrixType*)Mout->data, *(BADoubleTriDiagonalMatrixType*)Min->data);
        }
        break;
        case 'd':
        {
          Mout->data = (void*) malloc(sizeof(DoubleTriDiagonalMatrixType));
          BADoubleTriDiagonalMatrixToDoubleTriDiagonalMatrix((BADoubleTriDiagonalMatrixType*)Min->data, (DoubleTriDiagonalMatrixType*)Mout->data, ((BADoubleTriDiagonalMatrixType*)(Min->data))->mu, T);
        }
        break;
        default:
        {
          printf("Error in ConvertResponsefunction unrecognized type to convert to: %c\n",type);
        }
        break;
      }
    }
    break;
    case 'd': // original double tri-diagonal (DoubleTriDiagonalMatrixType)
    {
      switch (type|32)
      {
        case 'l':
        {
          Mout->data = (void*) malloc(sizeof(ListOfPoles));
          DoubleTriDiagonalMatrixToListOfPoles(*(DoubleTriDiagonalMatrixType*)Min->data, (ListOfPoles*)Mout->data);
        }
        break;
        case 't':
        {
          Mout->data = (void*) malloc(sizeof(TriDiagonalMatrixType));
          DoubleTriDiagonalMatrixToTriDiagonalMatrix(*(DoubleTriDiagonalMatrixType*)Min->data, (TriDiagonalMatrixType*)Mout->data);
        }
        break;
        case 'n':
        {
          Mout->data = (void*) malloc(sizeof(BADoubleTriDiagonalMatrixType));
          DoubleTriDiagonalMatrixToBADoubleTriDiagonalMatrix(*(DoubleTriDiagonalMatrixType*)Min->data, (BADoubleTriDiagonalMatrixType*)Mout->data, ((DoubleTriDiagonalMatrixType*)(Min->data))->mu, T);
        }
        break;
        case 'a':
        {
          Mout->data = (void*) malloc(sizeof(AndersonMatrixType));
          DoubleTriDiagonalMatrixToAndersonMatrix(*(DoubleTriDiagonalMatrixType*)Min->data, (AndersonMatrixType*)Mout->data);
        }
        break;
        case 'd':
        {
          Mout->data = (void*) malloc(sizeof(DoubleTriDiagonalMatrixType));
          DoubleTriDiagonalMatrixCopy((DoubleTriDiagonalMatrixType*)Mout->data, *(DoubleTriDiagonalMatrixType*)Min->data);
        }
        break;
        default:
        {
          printf("Error in ConvertResponsefunction unrecognized type to convert to: %c\n",type);
        }
        break;
      }
    }
    break;
    default:
    {
      printf("Error in ConvertResponsefunction unrecognized type of M: %c\n",Min->type);
    }
    break;
  }
  if(Mout->data==NULL)
  {
    if(Mout->NSuggestedGridPoints!=0)
    {
      free(Mout->SuggestedGrid);
      Mout->SuggestedGrid=NULL;
      Mout->NSuggestedGridPoints=0;
    }
    return 1;
  }
  return 0;
}


/*
 *
 */
int ResponsefunctionCopy(ResponsefunctionType *Min, ResponsefunctionType *Mout)
{
  return ConvertResponsefunction(Min,Mout,Min->type, "Default", 0, 0, 0, 0, false, false, false);
}

/*  set elements smaller than epsilon in absolute value to zero
 *  If Deflate is true also remove a block when all elements are zero
 */
int ResponsefunctionChop(ResponsefunctionType *M, double epsilon, bool Deflate)
{
  if(M->Block==true)
  {
    switch (M->type|32)
    {
      case 'l':
      {
        if(Deflate==true)
        {
          ChopBlockListOfPolesReduceDimension((BlockListOfPoles*)(M->data), epsilon);
        }
        else
        {
          ChopBlockListOfPoles((BlockListOfPoles*)(M->data), epsilon);
        }
      }
        break;
      case 't':
      {
        if(Deflate==true)
        {
          BlockTridiagonalMatrixChopReduceDimension((BlockTriDiagonalMatrixType*)(M->data), epsilon);
        }
        else
        {
          BlockTridiagonalMatrixChop((BlockTriDiagonalMatrixType*)(M->data), epsilon);
        }
      }
        break;
      case 'a':
      {
        if(Deflate==true)
        {
          BlockAndersonMatrixChopReduceDimension((BlockAndersonMatrixType*)(M->data), epsilon);
        }
        else
        {
          BlockAndersonMatrixChop((BlockAndersonMatrixType*)(M->data), epsilon);
        }
      }
        break;
      case 'n':
      {
        if(Deflate==true)
        {
          BlockNaturalImpurityOrbitalMatrixChopReduceDimension((BlockNaturalImpurityOrbitalMatrixType*)(M->data), epsilon);
        }
        else
        {
          BlockNaturalImpurityOrbitalMatrixChop((BlockNaturalImpurityOrbitalMatrixType*)(M->data), epsilon);
        }
      }
        break;
      case 'd':
      {
        if(Deflate==true)
        {
          BlockDoubleTriDiagonalMatrixChopReduceDimension((BlockDoubleTriDiagonalMatrixType*)(M->data), epsilon);
        }
        else
        {
          BlockDoubleTriDiagonalMatrixChop((BlockDoubleTriDiagonalMatrixType*)(M->data), epsilon);
        }
      }
        break;
      default:
      {
        printf("Error in ResponsefunctionChop unrecognized type to chop: %c\n",M->type);
      }
        break;
    }
  }
  else
  {
    switch (M->type|32)
    {
      case 'l':
      {
        if(Deflate==true)
        {
          ListOfPoleChopReduceDimension((ListOfPoles*)(M->data), epsilon);
        }
        else
        {
          ListOfPoleChop((ListOfPoles*)(M->data), epsilon);
        }
      }
        break;
      case 't':
      {
        if(Deflate==true)
        {
          TriDiagonalMatrixChopReduceDimension((TriDiagonalMatrixType *)(M->data), epsilon);
        }
        else
        {
          TriDiagonalMatrixChop((TriDiagonalMatrixType *)(M->data), epsilon);
        }
      }
        break;
      case 'a':
      {
        if(Deflate==true)
        {
          AndersonMatrixChopReduceDimension((AndersonMatrixType *)(M->data), epsilon);
        }
        else
        {
          AndersonMatrixChop((AndersonMatrixType *)(M->data), epsilon);
        }
      }
        break;
      case 'n':
      {
        if(Deflate==true)
        {
          BADoubleTriDiagonalMatrixChopReduceDimension((BADoubleTriDiagonalMatrixType *)(M->data), epsilon);
        }
        else
        {
          BADoubleTriDiagonalMatrixChop((BADoubleTriDiagonalMatrixType *)(M->data), epsilon);
        }
      }
        break;
      case 'd':
      {
        if(Deflate==true)
        {
          DoubleTriDiagonalMatrixChopReduceDimension((DoubleTriDiagonalMatrixType *)(M->data), epsilon);
        }
        else
        {
          DoubleTriDiagonalMatrixChop((DoubleTriDiagonalMatrixType *)(M->data), epsilon);
        }
      }
        break;
      default:
      {
        printf("Error in ResponsefunctionChop unrecognized type to chop: %c\n",M->type);
      }
        break;
    }
  }
  return 0;
}


/*
 *
 */
int ResponsefunctionMul(ResponsefunctionType *M, double x)
{
  if(M->Block==true)
  {
    switch (M->type|32)
    {
      case 'l':
      {
        if( (((BlockListOfPoles*)(M->data))->QComplex) ==0 )
        {
          int k=0;
          for(int i=0;i<(((BlockListOfPoles*)(M->data))->N); i++)
          {
            k++;
            for(int j=0;j<Sqr( (((BlockListOfPoles*)(M->data))->BlockSize) ); j++)
            {
              ((BlockListOfPoles*)(M->data))->VRe[k]*=x;
              k++;
            }
          }
          for(int j=0;j<Sqr( (((BlockListOfPoles*)(M->data))->BlockSize) ); j++)
          {
            ((BlockListOfPoles*)(M->data))->a0Re[j]*=x;
          }
        }
        else
        {
          int k=0;
          for(int i=0;i<(((BlockListOfPoles*)(M->data))->N); i++)
          {
            k++;
            for(int j=0;j<Sqr( (((BlockListOfPoles*)(M->data))->BlockSize) ); j++)
            {
              ((BlockListOfPoles*)(M->data))->V[k].Re*=x;
              ((BlockListOfPoles*)(M->data))->V[k].Im*=x;
              k++;
            }
          }
          for(int j=0;j<Sqr( (((BlockListOfPoles*)(M->data))->BlockSize) ); j++)
          {
            ((BlockListOfPoles*)(M->data))->a0[j].Re*=x;
            ((BlockListOfPoles*)(M->data))->a0[j].Im*=x;
          }
        }
      }
        break;
      case 't':
      {
        double xSqr = x;
        x = sqrt(x);
        if( (((BlockTriDiagonalMatrixType*)(M->data))->QComplex) ==0 )
        {
          for(int i=0;i<Sqr(((BlockTriDiagonalMatrixType*)(M->data))->BlockSize);i++)
          {
            (((BlockTriDiagonalMatrixType*)(M->data))->bRe[i]*=x);
            (((BlockTriDiagonalMatrixType*)(M->data))->aRe[i]*=xSqr);
          }
        }
        else
        {
          for(int i=0;i<Sqr(((BlockTriDiagonalMatrixType*)(M->data))->BlockSize);i++)
          {
            (((BlockTriDiagonalMatrixType*)(M->data))->b[i].Re*=x);
            (((BlockTriDiagonalMatrixType*)(M->data))->b[i].Im*=x);
            (((BlockTriDiagonalMatrixType*)(M->data))->a[i].Re*=xSqr);
            (((BlockTriDiagonalMatrixType*)(M->data))->a[i].Im*=xSqr);
          }
        }
      }
        break;
      case 'a':
      {
        double xSqr = x;
        x = sqrt(x);
        if( (((BlockAndersonMatrixType*)(M->data))->QComplex) ==0 )
        {
          for(int i=0;i<Sqr(((BlockAndersonMatrixType*)(M->data))->BlockSize);i++)
          {
            (((BlockAndersonMatrixType*)(M->data))->bRe[i]*=x);
            (((BlockAndersonMatrixType*)(M->data))->aRe[i]*=xSqr);
          }
        }
        else
        {
          for(int i=0;i<Sqr(((BlockAndersonMatrixType*)(M->data))->BlockSize);i++)
          {
            (((BlockAndersonMatrixType*)(M->data))->b[i].Re*=x);
            (((BlockAndersonMatrixType*)(M->data))->b[i].Im*=x);
            (((BlockAndersonMatrixType*)(M->data))->a[i].Re*=xSqr);
            (((BlockAndersonMatrixType*)(M->data))->a[i].Im*=xSqr);
          }
        }
      }
        break;
      case 'n':
      {
        double xSqr = x;
        x = sqrt(x);
        if( (((BlockNaturalImpurityOrbitalMatrixType*)(M->data))->QComplex) ==0 )
        {
          for(int i=0;i<Sqr(((BlockNaturalImpurityOrbitalMatrixType*)(M->data))->BlockSize);i++)
          {
            (((BlockNaturalImpurityOrbitalMatrixType*)(M->data))->val.bRe[i]*=x);
            (((BlockNaturalImpurityOrbitalMatrixType*)(M->data))->con.bRe[i]*=x);
            (((BlockNaturalImpurityOrbitalMatrixType*)(M->data))->a0.VRe[i]*=xSqr);
          }
        }
        else
        {
          for(int i=0;i<Sqr(((BlockNaturalImpurityOrbitalMatrixType*)(M->data))->BlockSize);i++)
          {
            (((BlockNaturalImpurityOrbitalMatrixType*)(M->data))->val.b[i].Re*=x);
            (((BlockNaturalImpurityOrbitalMatrixType*)(M->data))->val.b[i].Im*=x);
            (((BlockNaturalImpurityOrbitalMatrixType*)(M->data))->con.b[i].Re*=x);
            (((BlockNaturalImpurityOrbitalMatrixType*)(M->data))->con.b[i].Im*=x);
            (((BlockNaturalImpurityOrbitalMatrixType*)(M->data))->a0.V[i].Re*=xSqr);
            (((BlockNaturalImpurityOrbitalMatrixType*)(M->data))->a0.V[i].Im*=xSqr);
          }
        }      }
        break;
      case 'd':
      {
        /*  G = a0 + b0^* ( w - Gval - Gcon )^-1 b0^T, so a0 takes x and b0 the square root
         *  of it, b0 appearing twice. The two chains stand inside the denominator and are
         *  left alone. a0 and b0 carry a QComplex each, which for a mixed matrix need not
         *  be the QComplex of the matrix itself. */
        double xSqrt = sqrt(x);
        BlockDoubleTriDiagonalMatrixType *D = (BlockDoubleTriDiagonalMatrixType*)(M->data);
        for(int i=0;i<Sqr(D->BlockSize);i++)
        {
          if(D->a0.QComplex) { D->a0.V[i].Re*=x;     D->a0.V[i].Im*=x;     }
          else                 D->a0.VRe[i]*=x;
          if(D->b0.QComplex) { D->b0.V[i].Re*=xSqrt; D->b0.V[i].Im*=xSqrt; }
          else                 D->b0.VRe[i]*=xSqrt;
        }
      }
        break;
      default:
      {
        printf("Error in ResponsefunctionMul unrecognized type: %c\n",M->type);
        fflush(stdout);
        return 1;
      }
        break;
    }
  }
  else
  {
    switch (M->type|32)
    {
      case 'l':
      {
        for(int i=0;i<((ListOfPoles*)(M->data))->N;i++)
        {
          ((ListOfPoles*)(M->data))->b[i]*=x;
        }
        ((ListOfPoles*)(M->data))->a0*= x;
      }
        break;
      case 't':
      {
        ((TriDiagonalMatrixType*)(M->data))->b[0] *=sqrt(x);
        ((TriDiagonalMatrixType*)(M->data))->a[0] *= x;
      }
        break;
      case 'a':
      {
        ((AndersonMatrixType*)(M->data))->b[0] *=sqrt(x);
        ((AndersonMatrixType*)(M->data))->a[0] *= x;
      }
        break;
      case 'n':
      {
        ((BADoubleTriDiagonalMatrixType*)(M->data))->b0 *=sqrt(x);
        ((BADoubleTriDiagonalMatrixType*)(M->data))->a0 *= x;
      }
        break;
      case 'd':
      {
        /*  as the block case: a0 takes x, b0 the square root of it, the two chains nothing */
        ((DoubleTriDiagonalMatrixType*)(M->data))->b0 *=sqrt(x);
        ((DoubleTriDiagonalMatrixType*)(M->data))->a0 *= x;
      }
        break;
      default:
      {
        printf("Error in ResponsefunctionMul unrecognized type: %c\n",M->type);
        fflush(stdout);
        return 1;
      }
        break;
    }
  }
  return 0;
}




/*   FreeResponsefunction frees the content of a response function
 *   Response functions of types with a capital letter are freed
 *   Response functions of types with small letters are controled by a different instance and not freed
 *   The original variable of ResponsefunctionType is not freed
 */
uint32_t FreeResponsefunction(ResponsefunctionType *M)
{
  if(M->Block==true) return FreeBlockMatrixFromLua(M);
  if((M->type>='A') && (M->type<='Z') && (M->NSuggestedGridPoints!=0))
  {
    free(M->SuggestedGrid);
    M->SuggestedGrid=NULL;
    M->NSuggestedGridPoints=0;
  }
  switch (M->type) {
    case 'l':
    {
      return 0;
    }
    case 'L':
    {
      FreeListOfPoles((ListOfPoles *)M->data);
      free((ListOfPoles *)(M->data));
      return 0;
    }
    case 'a':
    {
      return 0;
    }
    case 'A':
    {
      FreeAndersonMatrix((AndersonMatrixType *)M->data);
      free((AndersonMatrixType *)(M->data));
      return 0;
    }
    case 't':
    {
      return 0;
    }
    case 'T':
    {
      FreeTriDiagonalMatrix((TriDiagonalMatrixType *)M->data);
      free((TriDiagonalMatrixType *)(M->data));
      return 0;
    }
    case 'n':
    {
      return 0;
    }
    case 'N':
    {
      FreeBADoubleTriDiagonalMatrix((BADoubleTriDiagonalMatrixType *)M->data);
      free((BADoubleTriDiagonalMatrixType *)(M->data));
      return 0;
    }
    case 'd':
    {
      return 0;
    }
    case 'D':
    {
      FreeDoubleTriDiagonalMatrix((DoubleTriDiagonalMatrixType *)M->data);
      free((DoubleTriDiagonalMatrixType *)(M->data));
      return 0;
    }
    default:
    {
      printf("Error in FreeResponsefunction unrecognized type %c\n",M->type);
      return 1;
    }
  }
}


/*   FreeResponsefunction frees the content of a response function
 *   Response functions of types with a capital letter are freed
 *   Response functions of types with small letters are controled by a different instance and not freed
 *   The original variable of ResponsefunctionType is freed if the content was freed. Otherwise it is not freed
 */
uint32_t FreeResponsefunctionIncludingOriginalPointer(ResponsefunctionType *M)
{
  if((M->type>='A') && (M->type<='Z'))
  {
    uint32_t retval;
    retval=FreeResponsefunction(M);
    free(M);
    return retval;
  }
  else
  {
    return FreeResponsefunction(M);
  }
}
