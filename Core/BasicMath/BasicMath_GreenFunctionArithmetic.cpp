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

/*  Superseded by TriDiagonalMatrixChainInvertEnergy and by
 *  BlockTriDiagonalMatrixChainInvertEnergy, which do the same and transpose as they go:
 *  A_i -> 2 mu - A_i^T rather than 2 mu - A_i, and B_i -> B_i^*. For a single valued chain
 *  and for a real symmetric one the two agree; for a complex matrix valued one they do not,
 *  and it is the transposing pair that is the energy inversion. These two also leave the
 *  constant term of G alone, where the new ones map it to -A0^T.
 *
 *  One chain of a double tri-diagonal matrix under the energy inversion as it was defined
 *  then, G(mu+w) -> -G(mu-w)^+, with ^+ the conjugate transpose.
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
static void DoubleTriDiagonalMatrixChainInvertEnergyOld(TriDiagonalMatrixType * const T, double const Mu)
{
  T->a[0] = Mu - T->a[0];
  for(int i=1;i<=T->N;i++)
  {
    T->a[i] = 2 * Mu - T->a[i];
  }
  return;
}
static void BlockDoubleTriDiagonalMatrixChainInvertEnergyOld(BlockTriDiagonalMatrixType * const T, double const Mu)
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
/*  One square block of a chain, seen as a CompactMatrixType so that the matrix routines can
 *  be used on it. V points into the flat array of blocks and has to be taken from aRe or
 *  from a according to QComplex, which is the callers job: the two are a union and index
 *  with a different stride.
 */
/*  A -> Shift 1 - A^T, the energy inversion of one block. The shift sits on the diagonal,
 *  so it does not matter that it is added after the transpose, and NShift is the number of
 *  rows that carry it: a block that was deflated keeps its storage but only part of it is
 *  in use. With NShift = 0 this is A -> -A^T, which is what the constant term of a response
 *  function follows.
 */
static void CompactMatrixInvertEnergy(CompactMatrixType * const M, uint32_t const NShift, double const Shift)
{
  CompactMatrixTranspose(M);
  if(M->QComplex!=0)
  {
    for(uint32_t k=0;k<M->Ni*M->Nj;k++)
    {
      M->V[k].Re *= -1;
      M->V[k].Im *= -1;
    }
    for(uint32_t k=0;k<NShift;k++) M->V[k*(1+M->Nj)].Re += Shift;
  }
  else
  {
    for(uint32_t k=0;k<M->Ni*M->Nj;k++) M->VRe[k] *= -1;
    for(uint32_t k=0;k<NShift;k++) M->VRe[k*(1+M->Nj)] += Shift;
  }
  return;
}
/*  One chain of blocks:  a[i] -> 2 mu 1 - a[i]^T  and  b[i] -> b[i]^*.
 *
 *  a[0] is not part of the chain and what it is depends on the format, so the caller says
 *  what happens to it. For a tri-diagonal and an Anderson response function it is the
 *  constant term A0 and maps to -A0^T, which is ReflectHead with HeadShift = 0. The two
 *  chains of a double tri-diagonal one hold half of the onsite energy of the observed block
 *  each, so each half maps to mu - a[0] and HeadShift is mu. The two head sites of a
 *  natural impurity orbital matrix are not reflected at all and only follow the transpose.
 *
 *  M is a view on one block of the chain, which is stored in a flat array of them and not
 *  as a CompactMatrixType. The blocks are square, so CompactMatrixTranspose swaps in place
 *  and does not reallocate; it would free the pointer of a rectangular one, and that
 *  pointer stands in the middle of the array of the chain.
 */
static void BlockTriDiagonalMatrixChainInvertEnergy(BlockTriDiagonalMatrixType * const T, double const Mu, bool const ReflectHead, double const HeadShift)
{
  CompactMatrixType M;
  uint32_t const Bs2 = Sqr(T->BlockSize);
  M.Ni       = T->BlockSize;
  M.Nj       = T->BlockSize;
  M.QComplex = T->QComplex;
  for(uint32_t i=0;i<=T->N;i++)
  {
    if(T->QComplex!=0) M.V   = &(T->a[i*Bs2]);
    else               M.VRe = &(T->aRe[i*Bs2]);
    if(i==0)
    {
      if(ReflectHead==true) CompactMatrixInvertEnergy(&M, T->BlockSizeList[0], HeadShift);
      else                  CompactMatrixTranspose(&M);
    }
    else
    {
      CompactMatrixInvertEnergy(&M, T->BlockSizeList[i], 2*Mu);
    }
  }
  for(uint32_t i=0;i<T->N;i++)
  {
    if(T->QComplex!=0) M.V   = &(T->b[i*Bs2]);
    else               M.VRe = &(T->bRe[i*Bs2]);
    CompactMatrixConjugate(&M);
  }
  return;
}
/*  the same for a single valued chain, where the transpose is nothing and only a[0], the
 *  constant term, has changed compared to what this used to do */
static void TriDiagonalMatrixChainInvertEnergy(TriDiagonalMatrixType * const T, double const Mu, bool const ReflectHead, double const HeadShift)
{
  if(ReflectHead==true) T->a[0] = HeadShift - T->a[0];
  for(uint32_t i=1;i<=T->N;i++) T->a[i] = 2 * Mu - T->a[i];
  return;
}
/*  The energy inversion. Throughout, ^* is the element by element complex conjugate, ^T the
 *  transpose and ^+ the two combined. On the stored numbers the map reads
 *
 *      A0  -> -A0^T            the constant term, the value of G at infinity
 *      e_n -> 2 mu - e_n       the poles are reflected about mu
 *      w_n -> w_n^T            the residues
 *
 *  and in the continued fraction formats the same statement is A_i -> 2 mu - A_i^T for the
 *  blocks of a chain and B_i -> B_i^* for the blocks that couple them. As a function of a
 *  complex w that is
 *
 *      Ginv(mu + w) = -G(mu - w)^T
 *
 *  the transpose alone, and mu - w the plain reflection of the argument: nothing is
 *  conjugated, neither the matrix nor w. Evaluated the way Quanty evaluates, at a real
 *  frequency with +i Gamma added on both sides, the same statement reads
 *
 *      Ginv(mu + w, Gamma) = -[ G(mu - w, Gamma) ]^*
 *
 *  where the conjugate takes care of the sign of Gamma as well as of the matrix. That
 *  second form is a corollary and not the definition: it needs w_n^* = w_n^T, which holds
 *  because the residues are Hermitian. It is the form the test checks, since Quanty has no
 *  negative Gamma to offer.
 *
 *  It is ^T and not ^+ that makes this the map one wants. A response function is
 *  G_ij(w) = <0| A_i ( w - H + E0 )^-1 A_j^+ |0>, so its residues are w_n = v v^+ with
 *  v_i = <0|A_i|n> and are Hermitian. The other side of the Lehmann sum, which is what
 *  reflecting the energy axis asks for, is <0| A_j^+ ... A_i |0>: the index order swaps
 *  together with the reflection, w_n -> w_n^T. Written with ^+ the map would leave the
 *  Hermitian residues alone and give something that is not the response function of
 *  anything unless G is symmetric -- invisible for every real G and silently wrong for a
 *  complex one.
 *
 *  A0 follows the same rule rather than staying where it is, so the line above holds as it
 *  stands and not only for A0 = 0. G and a self energy without its Hartree part both have
 *  A0 = 0 and do not notice. Applying the map twice is still the identity: -(-A0^T)^T = A0.
 *
 *  The transpose is folded into the loops so that each of the five formats is touched once.
 */
int ResponseFunctionInvertEnergy(ResponsefunctionType * const Gin, ResponsefunctionType * const Gout)
{
  ResponsefunctionCopy(Gin,Gout);

  double Mu;
  if(Gout->Block==true)
  {
    switch (Gout->type|32)
    {
      case 'a':
      {
        BlockAndersonMatrixType *A = (BlockAndersonMatrixType*)(Gout->data);
        CompactMatrixType M;
        uint32_t const Bs2 = Sqr(A->BlockSize);
        Mu = A->mu;
        M.Ni       = A->BlockSize;
        M.Nj       = A->BlockSize;
        M.QComplex = A->QComplex;
        /*  a[0] is the constant term and maps to -A0^T, a[1] and up are the poles */
        for(uint32_t i=0;i<=A->N;i++)
        {
          if(A->QComplex!=0) M.V   = &(A->a[i*Bs2]);
          else               M.VRe = &(A->aRe[i*Bs2]);
          if(i==0) CompactMatrixInvertEnergy(&M, 0,             0.0);
          else     CompactMatrixInvertEnergy(&M, A->BlockSize, 2*Mu);
        }
        for(uint32_t i=0;i<A->N;i++)
        {
          if(A->QComplex!=0) M.V   = &(A->b[i*Bs2]);
          else               M.VRe = &(A->bRe[i*Bs2]);
          CompactMatrixConjugate(&M);
        }
      }
      break;
      case 'l':
      {
        BlockListOfPoles *L = (BlockListOfPoles*)(Gout->data);
        uint32_t const Bs2 = Sqr(L->BlockSize);
        CompactMatrixType M;
        Mu = L->mu;
        M.Ni       = L->BlockSize;
        M.Nj       = L->BlockSize;
        M.QComplex = L->QComplex;
        M.V        = L->a0;
        CompactMatrixInvertEnergy(&M, 0, 0.0);
        /*  pole i stands at (1+BlockSize*BlockSize)*i, its residue follows it */
        for(uint32_t i=0;i<L->N;i++)
        {
          if(L->QComplex!=0)
          {
            L->V[i*(1+Bs2)].Re  = 2*Mu - L->V[i*(1+Bs2)].Re;
            L->V[i*(1+Bs2)].Im *= -1;
            M.V = &(L->V[i*(1+Bs2)+1]);
          }
          else
          {
            L->VRe[i*(1+Bs2)] = 2*Mu - L->VRe[i*(1+Bs2)];
            M.VRe = &(L->VRe[i*(1+Bs2)+1]);
          }
          CompactMatrixTranspose(&M);
        }
      }
      break;
      case 't':
      {
        BlockTriDiagonalMatrixType *T = (BlockTriDiagonalMatrixType*)(Gout->data);
        /*  a[0] is the constant term: reflected about 0, that is negated, and transposed */
        BlockTriDiagonalMatrixChainInvertEnergy(T, T->mu, true, 0.0);
      }
      break;
      case 'n':
      {
        BlockNaturalImpurityOrbitalMatrixType *B = (BlockNaturalImpurityOrbitalMatrixType*)(Gout->data);
        BlockTriDiagonalMatrixType temp;
        /*  the head sites of the two chains are the natural orbitals themselves and this
         *  map leaves them where they are */
        BlockTriDiagonalMatrixChainInvertEnergy(&(B->val), B->mu, false, 0.0);
        BlockTriDiagonalMatrixChainInvertEnergy(&(B->con), B->mu, false, 0.0);
        CompactMatrixInvertEnergy(&(B->a0), 0, 0.0);
        temp   = B->con;
        B->con = B->val;
        B->val = temp;
        /*  TMT describes the two head sites and follows from the chains, so it is derived
         *  again from the chains as they now stand */
        SetTMTFromValandConInBlockNaturalImpurityOrbitalMatrixType(B);
      }
      break;
      case 'd':
      {
        BlockDoubleTriDiagonalMatrixType *D = (BlockDoubleTriDiagonalMatrixType*)(Gout->data);
        BlockTriDiagonalMatrixType temp;
        /*  a[0] of the two chains holds half of the onsite energy of the observed block
         *  each, and that energy maps to 2 mu - a[0]val - a[0]con, so mu - a[0] for each of
         *  the two keeps the half and half split. The b[i] would flip sign as well, which
         *  is a similarity transformation with diag(1,-1,1,-1,...) and leaves the head
         *  element of the resolvent, and with it G, alone; only the conjugate of the
         *  transpose is left. */
        BlockTriDiagonalMatrixChainInvertEnergy(&(D->val), D->mu, true, D->mu);
        BlockTriDiagonalMatrixChainInvertEnergy(&(D->con), D->mu, true, D->mu);
        CompactMatrixInvertEnergy(&(D->a0), 0, 0.0);
        CompactMatrixConjugate(&(D->b0));
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
    /*  a transpose is nothing here, so only the constant term has changed compared to what
     *  this used to do */
    switch (Gout->type|32)
    {
      case 'a':
      {
        AndersonMatrixType *A = (AndersonMatrixType*)(Gout->data);
        A->a[0] *= -1;
        for(uint32_t i=1;i<=A->N;i++) A->a[i] = 2 * A->mu - A->a[i];
      }
      break;
      case 'l':
      {
        ListOfPoles *L = (ListOfPoles*)(Gout->data);
        L->a0 *= -1;
        for(uint32_t i=0;i<L->N;i++) L->a[i] = 2 * L->mu - L->a[i];
      }
      break;
      case 't':
      {
        TriDiagonalMatrixType *T = (TriDiagonalMatrixType*)(Gout->data);
        TriDiagonalMatrixChainInvertEnergy(T, T->mu, true, 0.0);
      }
      break;
      case 'n':
      {
        BADoubleTriDiagonalMatrixType *B = (BADoubleTriDiagonalMatrixType*)(Gout->data);
        TriDiagonalMatrixType temp;
        TriDiagonalMatrixChainInvertEnergy(&(B->val), B->mu, false, 0.0);
        TriDiagonalMatrixChainInvertEnergy(&(B->con), B->mu, false, 0.0);
        B->a0 *= -1;
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
        TriDiagonalMatrixChainInvertEnergy(&(D->val), D->mu, true, D->mu);
        TriDiagonalMatrixChainInvertEnergy(&(D->con), D->mu, true, D->mu);
        D->a0 *= -1;
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

/*  The inverse of a Hermitian positive semi definite matrix, with the directions that carry
 *  no weight left at zero instead of blown up. It is CompactMatrixInverseSqrt one power
 *  further along: there the caller holds the square of the matrix, here the matrix itself,
 *  and inverting the matrix itself rather than its square saves squaring the condition
 *  number.
 */
static uint32_t HermitianCompactMatrixPseudoInverse(CompactMatrixType const * const A, CompactMatrixType * const Ainv, double const SingularValue)
{
  uint32_t i,j;
  double val[A->Ni],sumval=0;
  CompactMatrixType fun;
  if(CompactMatrixEigensystem(A, &fun, val)) return 1;
  for(i=0;i<A->Ni;i++)
  {
    sumval+=Abs(val[i]);
  }
  for(i=0;i<fun.Ni;i++)
  {
    if(Abs(val[i]) <= SingularValue*sumval)
    {
      val[i]=0;
      for(j=0;j<A->Nj;j++)
      {
        if(A->QComplex==0)
        {
          fun.VRe[j+i*A->Nj]=0;
        }
        else
        {
          fun.V[j+i*A->Nj].Re=0;
          fun.V[j+i*A->Nj].Im=0;
        }
      }
    }
    else
    {
      val[i]=1.0/val[i];
    }
  }
  TransformDiagonalMatrixByCompactMatrixTtranATcon(val, &fun, Ainv);
  FreeCompactMatrix(&fun);
  return 0;
}

/*  The largest element of a matrix, the crude magnitude used below to decide whether a
 *  matrix is zero and to set the scale of that decision.
 */
static double CompactMatrixMaxElement(CompactMatrixType const * const M)
{
  uint32_t i;
  double max=0,v;
  for(i=0;i<M->Ni*M->Nj;i++)
  {
    if(M->QComplex==0) v = Abs(M->VRe[i]);
    else               v = Abs(M->V[i].Re) + Abs(M->V[i].Im);
    if(v>max) max=v;
  }
  return max;
}

/*  Step (ii) of ResponseFunctionDownfold for an effective Hamiltonian that is a chain,
 *  Appendix C1 of Zacinskis et al. rather than its Appendix C3.
 *
 *  Heff is the chain Atilde1, Btilde1, Atilde2, ... and the congruence with the isometry P
 *  is wanted at a block size of Nsub rather than the Ntilde the chain is written in:
 *
 *    P Heff(w) P^+ = P Atilde1 P^+ + ( P Btilde1^* ) 1/( w - Ht ) ( P Btilde1^* )^+
 *
 *  with Ht the chain from its second block on. The head is one triple product. The rest is
 *  the resolvent of Ht seen through the Nsub row vectors D = P Btilde1^*, which is a block
 *  Lanczos of Ht started on those rows -- Eq. (C7), with the projection onto the subspace as
 *  the starting block vector. CompactMatrixBlockBandDiagonalize does exactly that and hands
 *  back a chain of block size Nsub whose head coupling b[0] is the S^{1/2} of that Lanczos
 *  and whose a[0] it leaves at zero for the caller to fill with the head.
 *
 *  That routine wants the conjugate of the left factor as its starting block, D^* and not
 *  D: it walks the Krylov space with v -> v A^T, so its row vectors live in the conjugated
 *  space, and it closes the chain with b[0] = v[0] v[1]^+. The list of poles route obeys
 *  the same rule, BlockListOfPolesToBlockTridiagonalMatrix hands it beta^T and not beta.
 *  A real response function cannot tell the two apart and a complex one can, which is what
 *  the complex Tri case of the test is there for.
 *
 *  Why this is worth a route of its own: the general step (ii) rotates a list of poles, and
 *  a chain of M blocks of size Ntilde only becomes a list of poles by diagonalising it, so
 *  the general route pays (M Ntilde)^3 and answers with a star of M Ntilde poles. Here
 *  nothing is diagonalised and the answer is a chain again, of at most M Ntilde / Nsub
 *  blocks of size Nsub.
 *
 *  The Krylov space of Ht started on D cannot leave the first i blocks of the chain after i
 *  steps, so ceil( dim Ht / Nsub ) steps exhaust it and the projection is then exact.
 *  NTriMax cuts the chain shorter than that, which is an approximation that keeps the low
 *  moments of the projected function and drops the high ones, as cutting a chain always is.
 *
 *  HeffP is not allocated on entry.
 */
static int BlockTriDiagonalMatrixDownfoldChain(BlockTriDiagonalMatrixType const * const Heff, CompactMatrixType const * const P, BlockTriDiagonalMatrixType * const HeffP, double const SingularValue, uint32_t const NTriMax)
{
  char name[] = "block of a chain";
  CompactMatrixType Ht, A1, D, Dfull, view;
  uint32_t const Ntilde = Heff->BlockSize;
  uint32_t const Nsub   = P->Ni;
  uint32_t i,j,NTri;
  int retval = 1;
  SetCompactMatrixEmpty(&Ht); SetCompactMatrixEmpty(&A1);
  SetCompactMatrixEmpty(&D);  SetCompactMatrixEmpty(&Dfull);
  /*  view is a window on a block of the chain and owns nothing, so it is never freed */
  SetCompactMatrixEmpty(&view);
  view.Ni = Ntilde; view.Nj = Ntilde; view.QComplex = Heff->QComplex; view.Name = name;

  /*  the head, A1 = P Atilde1 P^+ */
  if(Heff->QComplex) view.V = Heff->a; else view.VRe = Heff->aRe;
  if(CompactMatrixMultiply3('N','N','C',P,&view,P,&A1))
  {
    printf("CompactMatrixMultiply3 failed in BlockTriDiagonalMatrixDownfoldChain\n");
    goto Cleanup;
  }
  SquareCompactMatrixMakeHermitian(&A1);

  HeffP->BlockSize = Nsub;
  HeffP->QComplex  = A1.QComplex;
  if(Heff->N!=0)
  {
    /*  D^* = ( P Btilde1^* )^*, the starting block: the Nsub row vectors the chain behind
     *  the head is seen through, conjugated as the block Lanczos wants them */
    if(Heff->QComplex) view.V = Heff->b; else view.VRe = Heff->bRe;
    if(CompactMatrixMultiply('*','N',P,&view,&D))
    {
      printf("CompactMatrixMultiply failed in BlockTriDiagonalMatrixDownfoldChain\n");
      goto Cleanup;
    }
    /*  Ht, the chain from its second block on, with D^* in its first block */
    if(BlockTriDiagonalMatrixToCompactMatrix(Heff,&Ht))
    {
      printf("BlockTriDiagonalMatrixToCompactMatrix failed in BlockTriDiagonalMatrixDownfoldChain\n");
      goto Cleanup;
    }
    Dfull.Ni = Nsub; Dfull.Nj = Ht.Ni; Dfull.QComplex = D.QComplex;
    if(InitCompactMatrix(&Dfull)||SetCompactMatrixToZero(&Dfull))
    {
      printf("InitCompactMatrix failed in BlockTriDiagonalMatrixDownfoldChain\n");
      goto Cleanup;
    }
    for(i=0;i<Nsub;i++)
    {
      for(j=0;j<Ntilde;j++)
      {
        if(Dfull.QComplex) Dfull.V[i*Dfull.Nj+j]   = D.V[i*D.Nj+j];
        else               Dfull.VRe[i*Dfull.Nj+j] = D.VRe[i*D.Nj+j];
      }
    }
    NTri = (Ht.Ni + Nsub - 1)/Nsub;
    if(NTriMax!=NTriMaxNoLimit) NTri = Min(NTri,NTriMax);
    /*  CompactMatrixBlockBandDiagonalize eats its starting block, which is why Dfull is a
     *  copy of D and not D itself. NOrtho and ReOrthogonalize are what its other caller,
     *  DiagonalMatrixBlockTridiagonalize, passes: orthogonalise against every Krylov block
     *  and not only the last two. A chain that runs to the length asked for here is long
     *  enough that a Lanczos without it loses orthogonality and with it the poles at the
     *  far end of the chain. */
    if(CompactMatrixBlockBandDiagonalize(&Ht,&Dfull,HeffP,SingularValue,NTri,Ht.Ni,true,false))
    {
      printf("CompactMatrixBlockBandDiagonalize failed in BlockTriDiagonalMatrixDownfoldChain\n");
      goto Cleanup;
    }
    if((A1.QComplex!=0)&&(HeffP->QComplex==0)) RealBlockTridiagonalMatrixToComplex(HeffP);
  }
  else
  {
    /*  a head and nothing behind it */
    HeffP->N = 0;
    if(InitBlockTridiagonalMatrix(HeffP))
    {
      printf("InitBlockTridiagonalMatrix failed in BlockTriDiagonalMatrixDownfoldChain\n");
      goto Cleanup;
    }
  }
  /*  the head energy, which the Lanczos leaves at zero */
  for(i=0;i<Sqr(Nsub);i++)
  {
    if(HeffP->QComplex==0)
    {
      HeffP->aRe[i] = A1.VRe[i];
    }
    else
    {
      HeffP->a[i].Re = (A1.QComplex==0) ? A1.VRe[i] : A1.V[i].Re;
      HeffP->a[i].Im = (A1.QComplex==0) ? 0         : A1.V[i].Im;
    }
  }
  HeffP->mu = Heff->mu;
  retval = 0;
Cleanup:
  FreeCompactMatrix(&Ht);
  FreeCompactMatrix(&A1);
  FreeCompactMatrix(&D);
  FreeCompactMatrix(&Dfull);
  return retval;
}

/*  Downfold a response function onto a subspace: the sub block of its inverse, inverted
 *  again,
 *
 *    Gsub(w) = ( R G^-1(w) R^+ )^-1 ,
 *
 *  with R a matrix of Nsub rows and BlockSize columns and of full row rank. R = ( 0 | 1 )
 *  picks the lower right block of a two by two block structure; any other R picks the
 *  subspace its rows span. The rows need not be orthonormal, the identity above holds for
 *  any R of full row rank.
 *
 *  This is the symmetric improved estimator for the self energy of a discrete many body
 *  system. With G the augmented propagator of the operators a_m and q_m = [a_m,H1] and R
 *  selecting the q half, Gsub is Sigma(w) - Sigma^HF. Sigma^HF is not part of it: it is the
 *  off diagonal block of the spectral weight of G, which the caller reads from B0 as
 *  B0^* B0^T.
 *
 *  S comes back as the spectral weight of the result, sum_i W_i, which for the self energy
 *  is the sum rule S = <{q,q^+}> - Sigma^HF <{a,a^+}>^-1 Sigma^HF.
 *
 *  The route is the one of Appendix C of Zacinskis et al., written for a response function
 *  held as G = A0 + B0^* ( w - Heff(w) )^-1 B0^T, in which B0^* is the Hermitian root
 *  S~^{1/2} of the spectral weight of G:
 *
 *    (i)   inversion    G^-1(w) = S~^{-1/2} ( w - Heff(w) ) S~^{-1/2}, which is Heff and
 *                       costs nothing at all, ResponseFunctionEffectiveHamiltonian has it
 *    (ii)  projection   R G^-1(w) R^+ = S^{-1/2} ( w - P Heff(w) P^+ ) S^{-1/2}
 *    (iii) reinversion  Gsub(w)       = S^{1/2} ( w - P Heff(w) P^+ )^-1 S^{1/2}
 *
 *  with
 *
 *    S^-1 = R S~^-1 R^+ ,    P = S^{1/2} R S~^{-1/2} .
 *
 *  Appendix C carries the congruences of (ii) and (iii) separately, one with S~^{-1/2} and
 *  one with S^{1/2}. Composing them into the single P above is what leaves one rotation of
 *  Heff and nothing else to do. P is an isometry, P P^+ = S^{1/2} S^-1 S^{1/2} = 1,
 *  whatever S~ does, so the projection is perfectly conditioned even where S~^{-1/2} on its
 *  own is not.
 *
 *  Causality survives by construction: every weight of Heff is Hermitian and positive semi
 *  definite, P W P^+ is a congruence and stays so, and so does S^{1/2} ( ) S^{1/2}.
 *
 *  S~ is the Gram matrix of the states behind G and is near singular whenever two of them
 *  are near dependent -- for the self energy at small interaction the q states vanish
 *  outright. What that costs is accuracy in the directions of P, not a large answer: a small
 *  eigenvalue of S~ is a large one of S^-1 and so a small weight in S. Directions of S~
 *  below SingularValue are deflated away and carry no weight at all.
 *
 *  A0 has to vanish. A constant added to G is not a constant added to G^-1 and the three
 *  steps above do not hold for it; physical Green's functions have A0 = 0.
 *
 *  B0 has to be Hermitian, as it is for every representation the library builds, where it
 *  comes out of a symmetric orthonormalisation. A b0 written down by hand is only defined up
 *  to a gauge, B0 -> B0 V with V unitary leaves G alone but hides a unitary in Heff, and
 *  the caller has to take that gauge out before asking for a subspace of the inverse.
 *
 *  There are two routes through step (ii) and the type of Gin picks one. An effective
 *  Hamiltonian that is a chain, which is what a tri-diagonal Gin has, takes Appendix C1:
 *  BlockTriDiagonalMatrixDownfoldChain projects it with a block Lanczos, nothing is
 *  diagonalised, and Gout comes back a chain of at most (M Ntilde / Nsub) blocks. Every
 *  other type reaches step (ii) as a list of poles and takes Appendix C3, the congruence
 *  one weight at a time, and Gout comes back a star of as many poles as Heff has. NTriMax
 *  is a ceiling on the length of the chain the first route produces and the second route
 *  refuses it, there being no chain to cut; NTriMaxNoLimit asks for no ceiling, which
 *  leaves the projection exact.
 *
 *  Gout is therefore a chain for a tri-diagonal Gin and an Anderson matrix otherwise, in
 *  both cases where step (iii) leaves it and at no cost beyond it. ChangeType turns either
 *  into the discrete spectral representation Gsub(w) = sum_i W_i / ( w - a_i ), at the
 *  price of one dense diagonalisation. Gout and S are not allocated on entry.
 */
int ResponseFunctionDownfold(ResponsefunctionType * const Gin, CompactMatrixType const * const R, ResponsefunctionType * const Gout, CompactMatrixType * const S, double const T, double const SingularValue, double const MinimalPoleWeigth, double const MinimalPoleDistance, uint32_t const NTriMax)
{
  ResponsefunctionType Heff, HeffP;
  CompactMatrixType A0, B0, InvSqrtStilde, Y, Sinv, SqrtS, P, Pc;
  uint32_t *QComplexOut;
  double *bReOut;
  Complex16 *bOut;
  uint32_t BlockSizeOut;
  double sv = (SingularValue < 0) ? 1000*DBL_EPSILON : SingularValue;
  double a0max, b0max;
  uint32_t i;
  int retval = 1;
  /*  everything the cleanup at the end touches, empty, so that a failure anywhere can jump
   *  straight to it. A lower case type is one this routine does not own and does not free */
  SetCompactMatrixEmpty(&A0);   SetCompactMatrixEmpty(&B0);  SetCompactMatrixEmpty(&InvSqrtStilde);
  SetCompactMatrixEmpty(&Y);    SetCompactMatrixEmpty(&Sinv); SetCompactMatrixEmpty(&SqrtS);
  SetCompactMatrixEmpty(&P);    SetCompactMatrixEmpty(&Pc);
  Heff.type='l';  Heff.Block=true;  Heff.data=NULL;  Heff.SuggestedGrid=NULL;  Heff.NSuggestedGridPoints=0;
  HeffP.type='l'; HeffP.Block=true; HeffP.data=NULL; HeffP.SuggestedGrid=NULL; HeffP.NSuggestedGridPoints=0;

  if(Gin->Block==false)
  {
    printf("Error in ResponseFunctionDownfold: a single valued response function has no subspace to downfold onto\n");
    fflush(stdout);
    return 1;
  }
  if(ResponseFunctionEffectiveHamiltonian(Gin, &Heff, &A0, &B0, T, SingularValue, MinimalPoleWeigth, MinimalPoleDistance))
  {
    printf("ResponseFunctionEffectiveHamiltonian failed in ResponseFunctionDownfold\n");
    fflush(stdout);
    return 1;
  }
  if(R->Nj != B0.Ni)
  {
    printf("Error in ResponseFunctionDownfold: the matrix selecting the subspace has %i columns and the response function has a block size of %i\n",R->Nj,B0.Ni);
    goto Cleanup;
  }
  a0max = CompactMatrixMaxElement(&A0);
  b0max = CompactMatrixMaxElement(&B0);
  if(a0max > sv*(1.0+Sqr(b0max)))
  {
    printf("Error in ResponseFunctionDownfold: the response function has a constant term A0, its largest element is %E. A constant added to G is not a constant added to G^-1 and a subspace of the inverse is not defined by the three steps this routine takes\n",a0max);
    goto Cleanup;
  }
  /*  the Hermitian root of the spectral weight of G is B0^*: the block convention is
   *  G = A0 + B0^* ( w - Heff )^-1 B0^T and ( B0^* )^+ = B0^T. Everything below is written
   *  in terms of that root, so conjugate once here and be done with it */
  CompactMatrixConjugate(&B0);
  if(QCompactMatrixHermitian(B0,1E-8)==0)
  {
    printf("Error in ResponseFunctionDownfold: B0 of the response function is not Hermitian. It is then defined only up to a unitary gauge, which sits in Heff and not in B0, and the subspace of the inverse would be taken in the wrong basis\n");
    goto Cleanup;
  }
  SquareCompactMatrixMakeHermitian(&B0);
  if(HermitianCompactMatrixPseudoInverse(&B0,&InvSqrtStilde,sv))
  {
    printf("HermitianCompactMatrixPseudoInverse failed in ResponseFunctionDownfold\n");
    goto Cleanup;
  }
  /*  Y = R S~^{-1/2}, from which both S and P follow: Y Y^+ = R S~^-1 R^+ = S^-1 */
  if(CompactMatrixMultiply('N','N',R,&InvSqrtStilde,&Y))
  {
    printf("CompactMatrixMultiply failed in ResponseFunctionDownfold\n");
    goto Cleanup;
  }
  if(CompactMatrixMultiply('N','C',&Y,&Y,&Sinv))
  {
    printf("CompactMatrixMultiply failed in ResponseFunctionDownfold\n");
    goto Cleanup;
  }
  SquareCompactMatrixMakeHermitian(&Sinv);
  CompactMatrixInverseSqrt(&Sinv,&SqrtS);
  if(CompactMatrixMultiply('N','N',&SqrtS,&Y,&P))
  {
    printf("CompactMatrixMultiply failed in ResponseFunctionDownfold\n");
    goto Cleanup;
  }
  /*  Step (ii). An effective Hamiltonian that is a chain takes the route of Appendix C1,
   *  which keeps it a chain and diagonalises nothing; everything else is a list of poles by
   *  the time it gets here and takes the congruence of Appendix C3, one weight at a time. */
  if((Heff.type|32)=='t')
  {
    HeffP.type  = 'T';
    HeffP.Block = true;
    HeffP.data  = malloc(sizeof(BlockTriDiagonalMatrixType));
    if(HeffP.data==NULL)
    {
      printf("malloc failed in ResponseFunctionDownfold\n");
      goto Cleanup;
    }
    if(BlockTriDiagonalMatrixDownfoldChain((BlockTriDiagonalMatrixType*)(Heff.data),&P,(BlockTriDiagonalMatrixType*)(HeffP.data),sv,NTriMax))
    {
      printf("BlockTriDiagonalMatrixDownfoldChain failed in ResponseFunctionDownfold\n");
      free(HeffP.data);
      HeffP.data = NULL;
      HeffP.type = 'l';
      goto Cleanup;
    }
  }
  else
  {
    if(NTriMax!=NTriMaxNoLimit)
    {
      printf("Error in ResponseFunctionDownfold: NTriMax is a ceiling on the length of the chain the projection produces, and a response function of type %c has no chain to produce -- it is downfolded one pole at a time and the result is a star\n",Gin->type);
      goto Cleanup;
    }
    /*  Rotate computes U^* X U^T, so U = P^* is the congruence P X P^+ */
    if(CompactMatrixCopy(P,&Pc))
    {
      printf("CompactMatrixCopy failed in ResponseFunctionDownfold\n");
      goto Cleanup;
    }
    CompactMatrixConjugate(&Pc);
    if(BlockResponseFunctionRotate(&Heff,&Pc,&HeffP,T,SingularValue,MinimalPoleWeigth,MinimalPoleDistance))
    {
      printf("BlockResponseFunctionRotate failed in ResponseFunctionDownfold\n");
      goto Cleanup;
    }
  }
  if(ResponseFunctionFromEffectiveHamiltonian(&HeffP,Gout,T,SingularValue,MinimalPoleWeigth,MinimalPoleDistance))
  {
    printf("ResponseFunctionFromEffectiveHamiltonian failed in ResponseFunctionDownfold\n");
    goto Cleanup;
  }
  /*  What FromEffectiveHamiltonian returns is ( w - P Heff(w) P^+ )^-1, with a head energy
   *  a[0] of zero and a head coupling b[0] of one. The S^{1/2} that step (iii) puts in front
   *  of it is exactly that head coupling, and writing it there is free. Turning the result
   *  into a list of poles instead, the pole expansion of Eq. (20), costs a dense
   *  diagonalisation of the whole matrix: the caller who wants that asks ChangeType for it
   *  and pays for it there.
   *
   *  A chain stays a chain and everything else has become a star, and the two hold their
   *  head energy and head coupling the same way, as the first block of a and of b.
   *
   *  b0 enters as b0^* ( ) b0^T, so the head coupling to write down is the conjugate of the
   *  Hermitian root S^{1/2}. P carries the QComplex of SqrtS into Heff and from there into
   *  the result, but a conversion on the way is free to notice that everything it holds is
   *  real, so promote when it did. */
  if((Gout->type|32)=='a')
  {
    BlockAndersonMatrixType * const Aout = (BlockAndersonMatrixType*)(Gout->data);
    if((Aout->QComplex==0)&&(SqrtS.QComplex!=0)) RealBlockAndersonMatrixToComplex(Aout);
    QComplexOut = &(Aout->QComplex); bReOut = Aout->bRe; bOut = Aout->b; BlockSizeOut = Aout->BlockSize;
  }
  else if((Gout->type|32)=='t')
  {
    BlockTriDiagonalMatrixType * const Tout = (BlockTriDiagonalMatrixType*)(Gout->data);
    if((Tout->QComplex==0)&&(SqrtS.QComplex!=0)) RealBlockTridiagonalMatrixToComplex(Tout);
    QComplexOut = &(Tout->QComplex); bReOut = Tout->bRe; bOut = Tout->b; BlockSizeOut = Tout->BlockSize;
  }
  else
  {
    printf("Error in ResponseFunctionDownfold: ResponseFunctionFromEffectiveHamiltonian returned a %c, which is neither a chain nor a star, and there is then no head coupling to put S^{1/2} in\n",Gout->type);
    FreeResponsefunction(Gout);
    goto Cleanup;
  }
  for(i=0;i<Sqr(BlockSizeOut);i++)
  {
    if(*QComplexOut==0)
    {
      bReOut[i] = SqrtS.VRe[i];
    }
    else
    {
      bOut[i].Re = (SqrtS.QComplex==0) ?  SqrtS.VRe[i] :  SqrtS.V[i].Re;
      bOut[i].Im = (SqrtS.QComplex==0) ?  0            : -SqrtS.V[i].Im;
    }
  }
  if(CompactMatrixMultiply('N','N',&SqrtS,&SqrtS,S))
  {
    printf("CompactMatrixMultiply failed in ResponseFunctionDownfold\n");
    FreeResponsefunction(Gout);
    goto Cleanup;
  }
  retval = 0;
Cleanup:
  FreeResponsefunction(&Heff);
  FreeResponsefunction(&HeffP);
  FreeCompactMatrix(&A0);
  FreeCompactMatrix(&B0);
  FreeCompactMatrix(&InvSqrtStilde);
  FreeCompactMatrix(&Y);
  FreeCompactMatrix(&Sinv);
  FreeCompactMatrix(&SqrtS);
  FreeCompactMatrix(&P);
  FreeCompactMatrix(&Pc);
  if(retval) printf("Error in ResponseFunctionDownfold\n");
  fflush(stdout);
  return retval;
}


int ResponseFunctionAdd(ResponsefunctionType * const GA, ResponsefunctionType * const GB, ResponsefunctionType * const Gsum)
{
  ResponsefunctionType GBLP;
  /*  two chains add up to a chain, without going through a list of poles. Lower case, as
   *  the list of poles below: the memory goes to Lua, which pushes it as a userdata. */
  if(((GA->type|32)=='t')&&((GB->type|32)=='t'))
  {
    int retval;
    Gsum->type = 't';
    Gsum->Block = GA->Block;
    Gsum->NSuggestedGridPoints = 0;
    Gsum->SuggestedGrid = NULL;
    Gsum->data = GA->Block ? malloc(sizeof(BlockTriDiagonalMatrixType)) : malloc(sizeof(TriDiagonalMatrixType));
    if(Gsum->data==NULL)
    {
      printf("malloc failed in ResponseFunctionAdd\n");
      return 1;
    }
    if(GA->Block)
    {
      retval = BlockTriDiagonalMatrixAdd((BlockTriDiagonalMatrixType*)(GA->data), (BlockTriDiagonalMatrixType*)(GB->data), (BlockTriDiagonalMatrixType*)(Gsum->data), 0);
    }
    else
    {
      retval = TriDiagonalMatrixAdd((TriDiagonalMatrixType*)(GA->data), (TriDiagonalMatrixType*)(GB->data), (TriDiagonalMatrixType*)(Gsum->data), 0);
    }
    if(retval)
    {
      free(Gsum->data);
      printf("Error in ResponseFunctionAdd\n");
    }
    return retval;
  }
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
  ConvertResponsefunction(&Heff, &HeffLP,  'L', "Default", -DBL_MAX, -DBL_MAX, -DBL_MAX, -DBL_MAX, false, false, false); // capital L as C has controle of pointer
  FreeResponsefunction(&Heff);
  ConvertResponsefunction(Sigma, &SigmaLP, 'L', "Default", -DBL_MAX, -DBL_MAX, -DBL_MAX, -DBL_MAX, false, false, false); // capital L as C has controle of pointer
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


/*  The smaller of a ceiling on the chain length and the length the conversion would
 *  produce on its own. Spelled out rather than written with the Min macro because the
 *  conversions below call their input Min.
 */
uint32_t NTriMaxOr(uint32_t const NTriMax, uint32_t const N)
{
  return (NTriMax<N)?NTriMax:N;
}

/*
 *
 */
int ConvertResponsefunction(const ResponsefunctionType *const Min, ResponsefunctionType *Mout, const char type, const char * const Method, const double T, const double SingularValue, const double MinimalPoleWeigth, const double MinimalPoleDistance, bool deflate, bool relativeDistance, bool relativeWeigth)
{
  return ConvertResponsefunction(Min, Mout, type, Method, T, SingularValue, MinimalPoleWeigth, MinimalPoleDistance, deflate, relativeDistance, relativeWeigth, NTriMaxNoLimit, NTriMaxNoLimit);
}

/*  NTriMaxVal and NTriMaxCon put a ceiling on the length of the chain(s) the conversion
 *  produces: NTriMaxVal on the one chain of a tri-diagonal matrix, and the two of them on
 *  the valence and the conduction chain of a NaturalImpurityOrbital or a DoubleTri matrix.
 *  A list of poles and an Anderson matrix have no chain and ignore both.
 *
 *  Where the routine that does the work takes a length of its own the ceiling is handed to
 *  it, so that it stops there rather than building a chain that is then thrown away. The
 *  cut at the end of this function is what makes the ceiling hold for the routes that have
 *  no such argument, and is a no-op for the ones that do.
 */
int ConvertResponsefunction(const ResponsefunctionType *const Min, ResponsefunctionType *Mout, const char type, const char * const Method, const double T, const double SingularValue, const double MinimalPoleWeigth, const double MinimalPoleDistance, bool deflate, bool relativeDistance, bool relativeWeigth, uint32_t const NTriMaxVal, uint32_t const NTriMaxCon)
{
  if(Min->Block==true)
  {
    return ConvertBlockMatrix(Min,Mout,type,Method, T,SingularValue, MinimalPoleWeigth, MinimalPoleDistance, deflate, relativeDistance, relativeWeigth, NTriMaxVal, NTriMaxCon);
  }
  Mout->NSuggestedGridPoints = 0;
  Mout->type = type;
  Mout->Block = false;
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
          ListOfPolesToTridiagonalMatrix(*(ListOfPoles*)Min->data, (TriDiagonalMatrixType*)Mout->data, NTriMaxOr(NTriMaxVal, ((ListOfPoles*)(Min->data))->N));
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
          /*  the two length arguments of ListOfPolesToBADoubleTridiagonalMatrix reduce the
           *  number of poles before the chains are built, which is a merge and not the cut
           *  NTriMax asks for, so they keep their old value and the cut happens at the end */
          ListOfPolesToBADoubleTridiagonalMatrix(*(ListOfPoles*)Min->data, (BADoubleTriDiagonalMatrixType*)Mout->data, ((ListOfPoles*)(Min->data))->N, ((ListOfPoles*)(Min->data))->N, ((ListOfPoles*)(Min->data))->mu, T);
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
          /*  a pole merge rather than a cut, as for the list of poles above */
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
          AndersonMatrixToTriDiagonalMatrix(*(AndersonMatrixType*)Min->data, (TriDiagonalMatrixType*)Mout->data, NTriMaxOr(NTriMaxVal, ((AndersonMatrixType*)Min->data)->N));
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
          /*  a pole merge rather than a cut, as for the list of poles above */
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
  return ResponsefunctionTruncate(Mout, NTriMaxVal, NTriMaxCon);
}


/*  Copy a response function in its single valued (non-block) form into the block form of the same type, with
 *  BlockSize=1. The numbers are copied, nothing is transformed, and the result evaluates to the same function:
 *    l -> L  ListOfPoles                 to BlockListOfPoles
 *    t -> T  TriDiagonalMatrixType       to BlockTriDiagonalMatrixType
 *    a -> A  AndersonMatrixType          to BlockAndersonMatrixType
 *    d -> D  DoubleTriDiagonalMatrixType to BlockDoubleTriDiagonalMatrixType
 *    n -> N  BADoubleTriDiagonalMatrixType to BlockNaturalImpurityOrbitalMatrixType
 *  The last one is not a pure copy, the two forms differ:
 *    single valued  G = a0 + b0^2 (Gval + Gcon)
 *    block          G = a0 +       Gval + Gcon
 *  so b0 is multiplied into the head of both chains, a[0] -> b0^2 a[0] and b[0] -> b0 b[0]. The derived numbers
 *  acvB, acvA, acvBA, bvB0, bcB0, bvA0 and bcA0 of the single valued form are not needed, the block form keeps
 *  their equivalent in TMT, which is calculated from val and con.
 *
 *  On exit Mout is a block response function of upper case type, i.e. C controls the allocation. Mout->data is
 *  allocated here.
 */
int ResponseFunctionToBlockResponseFunction(ResponsefunctionType const * const Min, ResponsefunctionType * const Mout)
{
  uint32_t i;
  if(Min->Block==true)
  {
    printf("Error in ResponseFunctionToBlockResponseFunction: the response function is already a block response function\n");
    return 1;
  }
  Mout->type = Min->type & ~32;  // upper case, C controls the allocation
  Mout->Block = true;
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
  switch (Min->type|32)
  {
    case 'l':
    {
      Mout->data = malloc(sizeof(BlockListOfPoles));
      if(CopyListOfPolesToBlockListOfPoles((ListOfPoles*)Min->data, (BlockListOfPoles*)Mout->data)) return 1;
    }
    break;
    case 't':
    {
      TriDiagonalMatrixType const * const A = (TriDiagonalMatrixType*)Min->data;
      BlockTriDiagonalMatrixType * const B = (BlockTriDiagonalMatrixType*)malloc(sizeof(BlockTriDiagonalMatrixType));
      Mout->data = B;
      B->N = A->N;
      B->BlockSize = 1;
      B->QComplex = 0;
      if(InitBlockTridiagonalMatrix(B)) return 1;
      strcpy(B->Name, A->Name);
      B->mu = A->mu;
      memcpy(B->aRe, A->a, (A->N+1)*sizeof(double));
      memcpy(B->bRe, A->b, (A->N)*sizeof(double));
    }
    break;
    case 'a':
    {
      AndersonMatrixType const * const A = (AndersonMatrixType*)Min->data;
      BlockAndersonMatrixType * const B = (BlockAndersonMatrixType*)malloc(sizeof(BlockAndersonMatrixType));
      Mout->data = B;
      B->N = A->N;
      B->BlockSize = 1;
      B->QComplex = 0;
      if(InitBlockAndersonMatrix(B)) return 1;
      strcpy(B->Name, A->Name);
      B->mu = A->mu;
      memcpy(B->aRe, A->a, (A->N+1)*sizeof(double));
      memcpy(B->bRe, A->b, (A->N)*sizeof(double));
    }
    break;
    case 'd':
    {
      DoubleTriDiagonalMatrixType const * const A = (DoubleTriDiagonalMatrixType*)Min->data;
      BlockDoubleTriDiagonalMatrixType * const B = (BlockDoubleTriDiagonalMatrixType*)malloc(sizeof(BlockDoubleTriDiagonalMatrixType));
      Mout->data = B;
      B->BlockSize = 1;
      B->QComplex = 0;
      B->val.N = A->val.N;
      B->con.N = A->con.N;
      if(InitBlockDoubleTriDiagonalMatrix(B)) return 1;
      strcpy(B->Name, A->Name);
      B->mu = A->mu;
      B->kBT = A->kBT;
      B->a0.VRe[0] = A->a0;
      B->b0.VRe[0] = A->b0;
      strcpy(B->val.Name, A->val.Name);
      B->val.mu = A->val.mu;
      memcpy(B->val.aRe, A->val.a, (A->val.N+1)*sizeof(double));
      memcpy(B->val.bRe, A->val.b, (A->val.N)*sizeof(double));
      strcpy(B->con.Name, A->con.Name);
      B->con.mu = A->con.mu;
      memcpy(B->con.aRe, A->con.a, (A->con.N+1)*sizeof(double));
      memcpy(B->con.bRe, A->con.b, (A->con.N)*sizeof(double));
    }
    break;
    case 'n':
    {
      BADoubleTriDiagonalMatrixType const * const A = (BADoubleTriDiagonalMatrixType*)Min->data;
      BlockNaturalImpurityOrbitalMatrixType * const B = (BlockNaturalImpurityOrbitalMatrixType*)malloc(sizeof(BlockNaturalImpurityOrbitalMatrixType));
      Mout->data = B;
      B->BlockSize = 1;
      B->QComplex = 0;
      B->val.N = A->val.N;
      B->con.N = A->con.N;
      if(InitBlockNaturalImpurityOrbitalMatrix(B)) return 1;
      strcpy(B->Name, A->Name);
      B->mu = A->mu;
      B->Epsilon = A->Epsilon;
      B->a0.VRe[0] = A->a0;
      strcpy(B->val.Name, A->val.Name);
      B->val.mu = A->val.mu;
      memcpy(B->val.aRe, A->val.a, (A->val.N+1)*sizeof(double));
      memcpy(B->val.bRe, A->val.b, (A->val.N)*sizeof(double));
      strcpy(B->con.Name, A->con.Name);
      B->con.mu = A->con.mu;
      memcpy(B->con.aRe, A->con.a, (A->con.N+1)*sizeof(double));
      memcpy(B->con.bRe, A->con.b, (A->con.N)*sizeof(double));
      // the overall scaling b0 of the single valued form goes into the head of both chains
      B->val.aRe[0] *= Sqr(A->b0);
      B->con.aRe[0] *= Sqr(A->b0);
      if(A->val.N>0) B->val.bRe[0] *= A->b0;
      if(A->con.N>0) B->con.bRe[0] *= A->b0;
      // TMT follows from val and con
      FreeCompactMatrix(&(B->TMT));
      if(SetTMTFromValandConInBlockNaturalImpurityOrbitalMatrixType(B)) return 1;
    }
    break;
    default:
    {
      printf("Error in ResponseFunctionToBlockResponseFunction: unrecognized type of response function: %c\n", Min->type);
      return 1;
    }
  }
  return 0;
}


/*  The direct sum of m block response functions A[i] of one type (NULL is a zero response function), each
 *  placed at the rows and columns Offset[i] ... Offset[i]+n_i-1 of a response function of block size
 *  BlockSize. See BlockListOfPolesDirectSum and its siblings. The A[i] must all be block response functions
 *  of the same type, which is checked, and have the same mu, which is left to the caller.
 *
 *  On exit Mout is a block response function of upper case type, i.e. C controls the allocation.
 */
int ResponseFunctionDirectSum(ResponsefunctionType const * const * const A, uint32_t const * const Offset, uint32_t const m, uint32_t const BlockSize, ResponsefunctionType * const Mout)
{
  uint32_t i, ret=1;
  char type=0;
  for(i=0;i<m;i++)
  {
    if(A[i]==NULL) continue;
    if(A[i]->Block!=true)
    {
      printf("Error in ResponseFunctionDirectSum: response function %i is not a block response function\n", i+1);
      return 1;
    }
    if(type==0) type = A[i]->type|32;
    if((A[i]->type|32)!=type)
    {
      printf("Error in ResponseFunctionDirectSum: response function %i is of type %c and not of type %c as the ones before\n", i+1, A[i]->type, type);
      return 1;
    }
  }
  if(type==0)
  {
    printf("Error in ResponseFunctionDirectSum: all response functions are zero, the type of the sum is not defined\n");
    return 1;
  }
  Mout->type = type & ~32;  // upper case, C controls the allocation
  Mout->Block = true;
  Mout->data = NULL;
  Mout->NSuggestedGridPoints = 0;
  Mout->SuggestedGrid = NULL;
  /*  Each case collects the data pointers in an array of its own type. An array of void* read back as an
   *  array of BlockListOfPoles* (or any other type) breaks the aliasing rules, and the Intel compiler
   *  then drops the stores and reads garbage. */
  switch(type)
  {
    case 'l':
    {
      BlockListOfPoles const **data = (BlockListOfPoles const **)malloc(m*sizeof(BlockListOfPoles*));
      for(i=0;i<m;i++) data[i] = (A[i]==NULL) ? (NULL) : ((BlockListOfPoles const *)A[i]->data);
      Mout->data = malloc(sizeof(BlockListOfPoles));
      ret = BlockListOfPolesDirectSum(data, Offset, m, BlockSize, (BlockListOfPoles*)Mout->data);
      free(data);
    }
    break;
    case 't':
    {
      BlockTriDiagonalMatrixType const **data = (BlockTriDiagonalMatrixType const **)malloc(m*sizeof(BlockTriDiagonalMatrixType*));
      for(i=0;i<m;i++) data[i] = (A[i]==NULL) ? (NULL) : ((BlockTriDiagonalMatrixType const *)A[i]->data);
      Mout->data = malloc(sizeof(BlockTriDiagonalMatrixType));
      ret = BlockTriDiagonalMatrixDirectSum(data, Offset, m, BlockSize, (BlockTriDiagonalMatrixType*)Mout->data);
      free(data);
    }
    break;
    case 'a':
    {
      BlockAndersonMatrixType const **data = (BlockAndersonMatrixType const **)malloc(m*sizeof(BlockAndersonMatrixType*));
      for(i=0;i<m;i++) data[i] = (A[i]==NULL) ? (NULL) : ((BlockAndersonMatrixType const *)A[i]->data);
      Mout->data = malloc(sizeof(BlockAndersonMatrixType));
      ret = BlockAndersonMatrixDirectSum(data, Offset, m, BlockSize, (BlockAndersonMatrixType*)Mout->data);
      free(data);
    }
    break;
    case 'n':
    {
      BlockNaturalImpurityOrbitalMatrixType const **data = (BlockNaturalImpurityOrbitalMatrixType const **)malloc(m*sizeof(BlockNaturalImpurityOrbitalMatrixType*));
      for(i=0;i<m;i++) data[i] = (A[i]==NULL) ? (NULL) : ((BlockNaturalImpurityOrbitalMatrixType const *)A[i]->data);
      Mout->data = malloc(sizeof(BlockNaturalImpurityOrbitalMatrixType));
      ret = BlockNaturalImpurityOrbitalMatrixDirectSum(data, Offset, m, BlockSize, (BlockNaturalImpurityOrbitalMatrixType*)Mout->data);
      free(data);
    }
    break;
    case 'd':
    {
      BlockDoubleTriDiagonalMatrixType const **data = (BlockDoubleTriDiagonalMatrixType const **)malloc(m*sizeof(BlockDoubleTriDiagonalMatrixType*));
      for(i=0;i<m;i++) data[i] = (A[i]==NULL) ? (NULL) : ((BlockDoubleTriDiagonalMatrixType const *)A[i]->data);
      Mout->data = malloc(sizeof(BlockDoubleTriDiagonalMatrixType));
      ret = BlockDoubleTriDiagonalMatrixDirectSum(data, Offset, m, BlockSize, (BlockDoubleTriDiagonalMatrixType*)Mout->data);
      free(data);
    }
    break;
    default:
      printf("Error in ResponseFunctionDirectSum: unrecognized type of response function: %c\n", type);
  }
  return ret;
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

/*  Cut the chain of a response function that carries one, so that it is at most NTriMaxVal
 *  (Tri), or its valence chain at most NTriMaxVal and its conduction chain at most
 *  NTriMaxCon (NaturalImpurityOrbital and DoubleTri). A chain that is already short enough
 *  is left alone, and NTriMaxNoLimit asks for no cut at all. What is dropped is dropped
 *  without looking at the values further down the chain.
 *
 *  A list of poles and an Anderson matrix have no chain: they are left as they are.
 */
int ResponsefunctionTruncate(ResponsefunctionType *M, uint32_t const NTriMaxVal, uint32_t const NTriMaxCon)
{
  if((NTriMaxVal==NTriMaxNoLimit)&&(NTriMaxCon==NTriMaxNoLimit)) return 0;
  if(M->Block==true)
  {
    switch (M->type|32)
    {
      case 't':
        return BlockTridiagonalMatrixTruncate((BlockTriDiagonalMatrixType*)(M->data), NTriMaxVal);
      case 'n':
        return BlockNaturalImpurityOrbitalMatrixTruncate((BlockNaturalImpurityOrbitalMatrixType*)(M->data), NTriMaxVal, NTriMaxCon);
      case 'd':
        return BlockDoubleTriDiagonalMatrixTruncate((BlockDoubleTriDiagonalMatrixType*)(M->data), NTriMaxVal, NTriMaxCon);
      case 'l':
      case 'a':
        return 0;
      default:
        printf("Error in ResponsefunctionTruncate unrecognized type to truncate: %c\n",M->type);
        return 1;
    }
  }
  switch (M->type|32)
  {
    case 't':
      return TriDiagonalMatrixTruncate((TriDiagonalMatrixType*)(M->data), NTriMaxVal);
    case 'n':
      return BADoubleTriDiagonalMatrixTruncate((BADoubleTriDiagonalMatrixType*)(M->data), NTriMaxVal, NTriMaxCon);
    case 'd':
      return DoubleTriDiagonalMatrixTruncate((DoubleTriDiagonalMatrixType*)(M->data), NTriMaxVal, NTriMaxCon);
    case 'l':
    case 'a':
      return 0;
    default:
      printf("Error in ResponsefunctionTruncate unrecognized type to truncate: %c\n",M->type);
      return 1;
  }
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
