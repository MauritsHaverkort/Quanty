
/*
*   MatrixToSparseMatrix         (MatrixType const M, SparseMatrixType * const U, double const epsilon)          M=U for M[i,j]>epsilon
*   MatrixToTriDiagonalMatrix    (MatrixType const M, MatrixType *const T, TriDiagonalMatrixType *const H)       M=T H T^{\dag}
*   MatrixToAndersonMatrix       (MatrixType const M, MatrixType *const T, AndersonMatrixType *const H)          M=T H T^{\dag}
*   MatrixToBlokTriDiagonalMatrix(MatrixType const M, MatrixType *const T, BlokTriDiagonalMatrixType *const H)   M=T H T^{\dag}
*   MatrixToBlokAndersonMatrix   (MatrixType const M, MatrixType *const T, BlokAndersonMatrixType *const H)      M=T H T^{\dag}
*   MatrixToEigensystem          (MatrixType const M, MatrixType *const T, double *val)                          M=T H T^{\dag} with H diagonal and elements val on the diagonal
*/


/****************************************************************************/
/********************         Matrix conversions         ********************/
/********************         AndersonMatrix To ...      ********************/
/****************************************************************************/

uint32_t AndersonMatrixToMatrix(AndersonMatrixType const M, MatrixType * const Mfull)
{
  int    i;
  Mfull->QComplex = 0;
  Mfull->Ni = M.N;
  Mfull->Nj = M.N;
  if (InitMatrix(Mfull))
  {
    printf("Could not initiallize Matrix in AndersonMatrixToMatrix\n");
    return 1;
  };
  strcpy(Mfull->Name, M.Name);
  for (i = 0; i < M.N; i++)
  {
    Mfull->VRe[i][i] = M.a[i + 1];
  }
  for (i = 1; i < M.N; i++)
  {
    Mfull->VRe[i][0] = M.b[i];
    Mfull->VRe[0][i] = M.b[i];
  }
  return 0;
}

uint32_t AndersonMatrixToTriDiagonalMatrixNonOrthogonal(AndersonMatrixType const M, TriDiagonalMatrixType * const H,
                                                        uint32_t NTriDiag)
{
  uint32_t i;
  double *vec[3];
  if ((vec[0] = malloc(M.N * sizeof(double))) == NULL)
  {
    printf("malloc 1 failed in AndersonMatrixToTriDiagonalMatrixNonOrthogonal\n");
    fflush(stdout);
    return 1;
  };
  if ((vec[1] = malloc(M.N * sizeof(double))) == NULL)
  {
    printf("malloc 1 failed in AndersonMatrixToTriDiagonalMatrixNonOrthogonal\n");
    free(vec[0]);
    fflush(stdout);
    return 1;
  };
  if ((vec[2] = malloc(M.N * sizeof(double))) == NULL)
  {
    printf("malloc 1 failed in AndersonMatrixToTriDiagonalMatrixNonOrthogonal\n");
    free(vec[0]);
    free(vec[1]);
    fflush(stdout);
    return 1;
  };
  H->N = NTriDiag;
  if (InitTriDiagonalMatrix(H))
  {
    printf("InitTriDiagonalMatrix failed in AndersonMatrixToTriDiagonalMatrixNonOrthogonal\n");
    free(vec[0]);
    free(vec[1]);
    free(vec[2]);
    fflush(stdout);
    return 1;
  };
  for (i = 0; i < M.N; i++)
  {
    vec[0][i] = 0;
    vec[1][i] = 0;
    vec[2][i] = 0;
  }
  vec[0][0] = 1;
  H->a[0] = M.a[0];
  H->b[0] = M.b[0];
  for (i = 1; i < NTriDiag; i++)
  {
    RealAndersonMatrixVector(M, vec[(i - 1) % 3], vec[i % 3]);  /* vec[i] = M vec[i-1]            */
    RealVectorAdd(vec[i % 3], vec[(i + 1) % 3], M.N, -(H->b[i - 1])); /* orthogonalization of vec[i] to vec[i-2], using that vec[i] at this point is H vec[i-1] */
    RealVectorDotVector(vec[(i - 1) % 3], vec[i % 3], M.N, &(H->a[i])); /* a[i] = vec[i] . vec[i-1]       */
    RealVectorAdd(vec[i % 3], vec[(i - 1) % 3], M.N, -H->a[i]); /* vec[i] -= a[i] vec[i-1]        */
    RealVectorDotVector(vec[i % 3], vec[i % 3], M.N, &(H->b[i])); /* b[i] = sqrt( vec[i] . vec[i] ) */
    H->b[i] = sqrt(H->b[i]);
    if (H->b[i] > H->b[1] * 100 * DBL_EPSILON)  /* vec[i] *= 1/b[i]               */
    {
      RealVectorMultiply(vec[i % 3], M.N, 1.0 / H->b[i]);
    }
    else
    {
      RealVectorMultiply(vec[i % 3], M.N, 0);
      H->b[i] = 0;
      H->a[i + 1] = 0;
      H->N = i;
      i = NTriDiag;
    }
  }
  RealVectorAndersonMatrixVector(M, vec[(H->N - 1) % 3], &(H->a[H->N]));  /* a[N] = vec[N-1] M vec[N-1]    */
  free(vec[0]);
  free(vec[1]);
  free(vec[2]);
  return 0;
}

/*
 *  M=T.H.T^{\dag}
 *  input M defined, H and T not initiallized
 */
uint32_t AndersonMatrixToTriDiagonalMatrixPlusT(AndersonMatrixType const M, MatrixType * const T,
                                                TriDiagonalMatrixType * const H, uint32_t NTriDiag)
{
  uint32_t i, j;
  double tempRe;
  T->QComplex = 0;
  T->Ni = NTriDiag;
  T->Nj = M.N;
  if (InitMatrix(T))
  {
    printf("InitMatrix failed in AndersonMatrixToTriDiagonalMatrix\n");
    fflush(stdout);
    return 1;
  };
  H->N = NTriDiag;
  if (InitTriDiagonalMatrix(H))
  {
    FreeMatrix(T);
    printf("InitTriDiagonalMatrix failed in AndersonMatrixToTriDiagonalMatrix\n");
    return 1;
  }
  strcpy(T->Name, "TriDiagonalizer of ");
  strncat(T->Name, M.Name, Max(0, MAXStrLength - 20));
  strcpy(H->Name, "TriDiagonal form of ");
  strncat(H->Name, M.Name, Max(0, MAXStrLength - 22));
  H->a[0] = M.a[0];
  H->b[0] = M.b[0];
  H->mu = M.mu;
  for (i = 0; i < M.N; i++)
  {
    T->VRe[0][i] = 0;
  }
  T->VRe[0][0] = 1;
  for (i = 1; i < NTriDiag; i++)
  {
    RealAndersonMatrixVector(M, T->VRe[i - 1], T->VRe[i]);
    for (j = 0; j < i - 1; j++)
    {
      RealVectorDotVector(T->VRe[j], T->VRe[i], T->Nj, &tempRe);
      RealVectorAdd(T->VRe[i], T->VRe[j], T->Nj, -tempRe);
    }
    RealVectorDotVector(T->VRe[i - 1], T->VRe[i], T->Nj, &(H->a[i]));
    RealVectorAdd(T->VRe[i], T->VRe[i - 1], T->Nj, -H->a[i]);
    RealVectorDotVector(T->VRe[i], T->VRe[i], T->Nj, &(H->b[i]));
    H->b[i] = sqrt(H->b[i]);
    if (H->b[i] > H->b[1] * 100 * DBL_EPSILON)
    {
      RealVectorMultiply(T->VRe[i], T->Nj, 1.0 / H->b[i]);
      for (j = 0; j < i; j++)
      {
        RealVectorDotVector(T->VRe[j], T->VRe[i], T->Nj, &tempRe);
        RealVectorAdd(T->VRe[i], T->VRe[j], T->Nj, -tempRe);
      }
      RealVectorDotVector(T->VRe[i], T->VRe[i], T->Nj, &tempRe);
      tempRe = sqrt(tempRe);
      H->b[i] *= tempRe;
      if (tempRe > 100 * DBL_EPSILON)
      {
        RealVectorMultiply(T->VRe[i], T->Nj, 1.0 / tempRe);
      }
      else
      {
        RealVectorMultiply(T->VRe[i], T->Nj, 0);
        H->b[i] = 0;
        H->a[i + 1] = 0;
        H->N = i + 1;
        T->Ni = i + 1;
        i = NTriDiag;
      }
    }
    else
    {
      RealVectorMultiply(T->VRe[i], T->Nj, 0);
      H->b[i] = 0;
      H->a[i + 1] = 0;
      H->N = i + 1;
      T->Ni = i + 1;
      i = NTriDiag;
    }
  }
  RealVectorAndersonMatrixVector(M, T->VRe[NTriDiag - 1], &(H->a[NTriDiag]));
  return 0;
}

uint32_t AndersonMatrixToTriDiagonalMatrix(AndersonMatrixType const M, TriDiagonalMatrixType * const H,
                                           uint32_t NTriDiag)
{
  MatrixType T;
  if (AndersonMatrixToTriDiagonalMatrixPlusT(M, &T, H, NTriDiag))
  {
    printf("AndersonMatrixToTriDiagonalMatrixPlusT failed in AndersonMatrixToTriDiagonalMatrix\n");
  };
  FreeMatrix(&T);
  return 0;
}

uint32_t AndersonMatrixToDoubleTriDiagonalMatrixNonOrthogonal(AndersonMatrixType const M,
                                                              DoubleTriDiagonalMatrixType * const H, double const Ef,
                                                              double const epsilon)
{
  /* To DO */
  return 1;
}

/*
*  M=T.H.T^{\dag}
*  input M defined, H and T not initiallized
*  note that when epsilon != 0 it can happen that M.N > H.N
*  in that case M is only roughly T.H.T^{\dag} as the transformation made is not unitary
*  The transformation made requires an aditional multiplication of Sqrt(0.5) to the doubled elements in H (only for the coupling)
*  and M.V and M.C are not orthonormal after the transformation. The overlap should be set to zero.
*  The eigenvalues of M contain the eigenvalues of H.
*  The aditional eigenvalues of M have zero values for the first element of their corresponding eigenfunctions.
*
*  diagonal energies in the AndersonMatrix larger  then Ef+epsilon belong to the conduction part of the double tri diagonal matrix
*  diagonal energies in the AndersonMatrix smaller then Ef-epsilon belong to the valence    part of the double tri diagonal matrix
*  diagonal energies in the AndersonMatrix in the interval Ef-epsilon to Ef+epsilon contribute partly to the valence, partly to the conductin part, weighted according to distance with Ef
*/
uint32_t AndersonMatrixToDoubleTriDiagonalMatrixPlusT(AndersonMatrixType const M, MatrixType * const T,
                                                      DoubleTriDiagonalMatrixType * const H, double const mu,
                                                      double const epsilon)
{
  AndersonMatrixType MV, MC;
  TriDiagonalMatrixType MVTri, MCTri;
  MatrixType AtoTV, AtoTC;
  uint32_t mapV[M.N + 1], mapC[M.N + 1];
  int    i, j;
  if (M.N <= 1)
  {
    printf("Can't create a double tri diagonal matrix of a constant\n");
    return 2;
  };
  MV.N = M.N;
  MC.N = M.N;
  if (InitAndersonMatrix(&MV))
  {
    printf("InitAndersonMatrix failed in AndersonMatrixToDoubleTriDiagonalMatrix\n");
    fflush(stdout);
    return 1;
  };
  if (InitAndersonMatrix(&MC))
  {
    FreeAndersonMatrix(&MV);
    printf("InitAndersonMatrix failed in AndersonMatrixToDoubleTriDiagonalMatrix\n");
    fflush(stdout);
    return 1;
  }
  mapV[0] = 0;
  mapC[0] = 0;
  MV.N = 1;
  MC.N = 1;
  MV.a[0] = 0;
  MV.b[0] = 1;
  MC.a[0] = 0;
  MC.b[0] = 1;
  MV.a[1] = M.a[1];
  MC.a[1] = M.a[1];
  for (i = 2; i <= M.N; i++)
  {
    if (fabs(M.b[i - 1]) > 100 * DBL_MIN)
    {
      if (M.a[i] < mu - epsilon)
      {
        MV.N++;
        MV.a[MV.N] = M.a[i];
        MV.b[MV.N - 1] = M.b[i - 1];
        mapV[MV.N - 1] = i - 1;
      }
      else if (M.a[i] > mu + epsilon)
      {
        MC.N++;
        MC.a[MC.N] = M.a[i];
        MC.b[MC.N - 1] = M.b[i - 1];
        mapC[MC.N - 1] = i - 1;
      }
      else
      {
        MV.N++;
        MV.a[MV.N] = M.a[i];
        MV.b[MV.N - 1] = M.b[i - 1] * sqrt(0.5 - (M.a[i] - mu) / (DBL_MIN + 2 * epsilon));
        mapV[MV.N - 1] = i - 1;
        MC.N++;
        MC.a[MC.N] = M.a[i];
        MC.b[MC.N - 1] = M.b[i - 1] * sqrt(0.5 + (M.a[i] - mu) / (DBL_MIN + 2 * epsilon));
        mapC[MC.N - 1] = i - 1;
      }
    }
  }
  if (MV.N == 1)
  {
    MV.N++;
    MV.a[MV.N] = mu - epsilon - (10 * fabs(mu) * DBL_EPSILON);
    MV.b[MV.N - 1] = 0;
    mapV[MV.N - 1] = 0;
  }
  if (MC.N == 1)
  {
    MC.N++;
    MC.a[MC.N] = mu + epsilon + (10 * fabs(mu) * DBL_EPSILON);
    MC.b[MC.N - 1] = 0;
    mapC[MC.N - 1] = 0;
  }
  AndersonMatrixToTriDiagonalMatrixPlusT(MV, &AtoTV, &MVTri, MV.N);
  AndersonMatrixToTriDiagonalMatrixPlusT(MC, &AtoTC, &MCTri, MC.N);
  T->Ni = MV.N + MC.N - 1;
  T->Nj = M.N;
  T->QComplex = 0;
  InitMatrix(T);
  for (i = 1; i < T->Ni; i++)
  {
    T->VRe[i][0] = 0;
  }
  T->VRe[0][0] = 1;
  for (i = 1; i < AtoTV.Ni; i++)
  {
    for (j = 1; j < AtoTV.Nj; j++)
    {
      T->VRe[i][mapV[j]] = AtoTV.VRe[i][j];
    }
  }
  for (i = 1; i < AtoTC.Ni; i++)
  {
    for (j = 1; j < AtoTC.Nj; j++)
    {
      T->VRe[i + AtoTV.Ni - 1][mapC[j]] = AtoTC.VRe[i][j];
    }
  }
  H->val.N = MVTri.N - 1;
  H->con.N = MCTri.N - 1;
  InitDoubleTriDiagonalMatrix(H);
  H->mu = mu;
  H->Epsilon = epsilon;
  H->a0 = M.a[0];
  H->b0 = M.b[0];
  for (i = 0; i < H->val.N; i++)
  {
    H->val.a[i] = MVTri.a[i + 1];
    H->val.b[i] = MVTri.b[i + 1];
  }
  H->val.a[H->val.N] = MVTri.a[H->val.N + 1];
  for (i = 0; i < H->con.N; i++)
  {
    H->con.a[i] = MCTri.a[i + 1];
    H->con.b[i] = MCTri.b[i + 1];
  }
  H->con.a[H->con.N] = MCTri.a[H->con.N + 1];
  FreeMatrix(&AtoTV);
  FreeMatrix(&AtoTC);
  FreeAndersonMatrix(&MV);
  FreeAndersonMatrix(&MC);
  FreeTriDiagonalMatrix(&MVTri);
  FreeTriDiagonalMatrix(&MCTri);
  return 0;
}

uint32_t AndersonMatrixToDoubleTriDiagonalMatrix(AndersonMatrixType const M, DoubleTriDiagonalMatrixType * const H,
                                                 double const Ef, double const epsilon)
{
  MatrixType T;
  if (AndersonMatrixToDoubleTriDiagonalMatrixPlusT(M, &T, H, Ef, epsilon))
  {
    printf("AndersonMatrixToDoubleTriDiagonalMatrixPlusT failed in AndersonMatrixToDoubleTriDiagonalMatrix\n");
    return 1;
  };
  FreeMatrix(&T);
  return 0;
}


/*  This algorithm assumes M.a[0]=0 and M.b[0]=1
 *  Gand = a0 + b0^2 1 / (w - a1 - \sum_i b_i^2 / ( w - a_i+1 ))
 *  a0=0 and b0=1 assumed
 *  -Gand^-1 = a1 + \sum_i b_i^2 / ( w - a_i+1 )
 *  Notice the minus sign
 */
uint32_t AndersonMatrixInverseToListOfPoles(AndersonMatrixType const M, ListOfPoles * const H)
{
  uint32_t i;
  H->N = M.N - 1;
  if (InitListOfPoles(H))
  {
    printf("InitListOfPoles failed in AndersonMatrixInverseToListOfPoles\n");
    return 1;
  }
  H->a0 = M.a[1];
  for (i = 0; i < H->N; i++)
  {
    H->b[i] = Sqr(M.b[i + 1]);
    H->a[i] = M.a[i + 2];
  }
  return 0;
}

/*
 *
 */
uint32_t AndersonMatrixToEigensystemOnsiteShift(AndersonMatrixType const M, double const shift, MatrixType * const T,
                                                double *val)
{
  MatrixType Tone, Ttwo;
  TriDiagonalMatrixType Mtri;
  if (AndersonMatrixToTriDiagonalMatrixPlusT(M, &Tone, &Mtri, M.N))
  {
    printf("AndersonMatrixToTriDiagonalMatrix failed in AndersonMatrixToEigensystem\n");
    return 1;
  };
  if (M.N > 0)
  {
    Mtri.a[1] += shift;
  }
  if (TriDiagonalMatrixToEigensystem(Mtri, &Ttwo, val))
  {
    printf("TriDiagonalMatrixToEigensystem failed in AndersonMatrixToEigensystem\n");
    return 1;
  };
  if (MatrixMultiply(T, Ttwo, Tone))
  {
    printf("MatrixMultiply failed in AndersonMatrixToEigensystem\n");
    return 1;
  };
  FreeMatrix(&Tone);
  FreeMatrix(&Ttwo);
  FreeTriDiagonalMatrix(&Mtri);
  return 0;
}

/*
 *
 */
uint32_t AndersonMatrixToEigensystem(AndersonMatrixType const M, MatrixType * const T, double *val)
{
  MatrixType Tone, Ttwo;
  TriDiagonalMatrixType Mtri;
  if (AndersonMatrixToTriDiagonalMatrixPlusT(M, &Tone, &Mtri, M.N))
  {
    printf("AndersonMatrixToTriDiagonalMatrix failed in AndersonMatrixToEigensystem\n");
    return 1;
  }
  if (TriDiagonalMatrixToEigensystem(Mtri, &Ttwo, val))
  {
    printf("TriDiagonalMatrixToEigensystem failed in AndersonMatrixToEigensystem\n");
    return 1;
  };
  if (MatrixMultiply(T, Ttwo, Tone))
  {
    printf("MatrixMultiply failed in AndersonMatrixToEigensystem\n");
    return 1;
  };
  FreeMatrix(&Tone);
  FreeMatrix(&Ttwo);
  FreeTriDiagonalMatrix(&Mtri);
  return 0;
}

/*
 *
 */
uint32_t AndersonMatrixToListOfPoles(AndersonMatrixType const M, ListOfPoles * const A)
{
  MatrixType T;
  uint32_t i;
  A->N = M.N + 1;
  if (InitListOfPoles(A))
  {
    printf("InitListOfPoles failed in TriDiagonalMatrixToListOfPoles\n");
    fflush(stdout);
    return 1;
  };
  if (AndersonMatrixToEigensystem(M, &T, A->a))
  {
    printf("AndersonMatrixToEigensystem failed in TriDiagonalMatrixToListOfPoles\n");
    fflush(stdout);
    return 1;
  };
  A->N = T.Ni;
  for (i = 0; i < A->N; i++)
  {
    A->b[i] = Sqr(M.b[0] * T.VRe[i][0]);
  }
  FreeMatrix(&T);
  A->mu=M.mu;
  return 0;
}

uint32_t TriDiagonalMatrixFindOnsiteOccupation(TriDiagonalMatrixType const M, double *occupation)
{
  uint32_t i;
  MatrixType Eigenfun;
  double Eigenval[M.N];
  TriDiagonalMatrixToEigensystem(M, &Eigenfun, Eigenval);
  *occupation = 0;
  for (i = 0; i < M.N; i++)
  {
    if (Eigenval[i] <= M.mu)
    {
      *occupation += Sqr(Eigenfun.VRe[i][0]);
    }
  }
  FreeMatrix(&Eigenfun);
  return 0;
}


/*   input:
 *     AndersonMatrix M           : the matrix for which an onsite shift is found such that the eigensystem has an occupation of "occupation"
 *     double         occupation  : the wanted occupation
 *     double *       onsiteshift : the starting value for the shift
 *
 *   output:
 *     double *       onsiteshift : the shift that brings the eigensystem of M to have an occupation of "occupation" on the impurity
 */
uint32_t AndersonMatrixFindOnsiteEnergyShift(AndersonMatrixType const M, double const occupation,
                                             double *const onsiteshift)
{
  uint32_t const maxit = 30;
  uint32_t i;
  double shifti[maxit], occi[maxit];
  shifti[0] = *onsiteshift;
  TriDiagonalMatrixType MTri;
  AndersonMatrixToTriDiagonalMatrix(M, &MTri, M.N);
  MTri.a[1] = M.a[1] + shifti[0];
  TriDiagonalMatrixFindOnsiteOccupation(MTri, &(occi[0]));
  shifti[1] = shifti[0] - (occupation - occi[0]) * Abs(MTri.b[1]);
  i = 0;
  while ((Abs(occupation - occi[i]) > 100 * DBL_EPSILON) && (i < maxit - 1))
  {
    i++;
    if (i != 1)
    {
      shifti[i] =
        shifti[i - 1] + (occupation - occi[i - 1]) * (shifti[i - 2] - shifti[i - 1]) / (occi[i - 2] - occi[i - 1]);
      shifti[i] = InterpolateRootOfset(shifti, occi, i, occupation);
    }
    MTri.a[1] = M.a[1] + shifti[i];
    TriDiagonalMatrixFindOnsiteOccupation(MTri, &(occi[i]));
  }
  *onsiteshift = shifti[i];
  if (i >= maxit - 1)
  {
    printf("In AndersonMatrixFindOnsiteEnergyShift onsite shift did not converge in %i itterations\n", maxit);
    printf("loop    shift                  occupation            difference\n");
    i++;
    while (i > 0)
    {
      i--;
      printf("%4i %22.15E %22.15E %22.15E\n", i, shifti[i], occi[i], occupation - occi[i]);
    }
  }
  return 0;
}

uint32_t EigenstateToTridiagonalMatrix(double const *psi0, double const *val, uint32_t const N,
                                       TriDiagonalMatrixType * const H)
{
  uint32_t i, j;
  MatrixType T;
  double tempRe;
  T.QComplex = 0;
  T.Ni = N;
  T.Nj = N;
  if (InitMatrix(&T))
    return 1;
  for (i = 0; i < N; i++)
  {
    T.VRe[0][i] = psi0[i];
  }
  H->N = N;
  if (InitTriDiagonalMatrix(H))
  {
    FreeMatrix(&T);
    return 1;
  }
  for (i = 1; i < H->N; i++)
  {
    for (j = 0; j < N; j++)
    {
      T.VRe[i][j] = val[j] * T.VRe[i - 1][j];
    }
    for (j = 0; j < i - 1; j++)
    {
      RealVectorDotVector(T.VRe[j], T.VRe[i], T.Nj, &tempRe);
      RealVectorAdd(T.VRe[i], T.VRe[j], T.Nj, -tempRe);
    }
    RealVectorDotVector(T.VRe[i - 1], T.VRe[i], T.Nj, &(H->a[i]));
    RealVectorAdd(T.VRe[i], T.VRe[i - 1], T.Nj, -H->a[i]);
    RealVectorDotVector(T.VRe[i], T.VRe[i], T.Ni, &(H->b[i]));
    H->b[i] = sqrt(H->b[i]);
    if (H->b[i] < 1000 * DBL_EPSILON)
    {
      RealVectorMultiply(T.VRe[i], T.Ni, 0);
      H->N = i;
    }
    else
    {
      RealVectorMultiply(T.VRe[i], T.Ni, 1.0 / H->b[i]);
    }
    for (j = 0; j < i; j++)
    {
      RealVectorDotVector(T.VRe[j], T.VRe[i], T.Nj, &tempRe);
      RealVectorAdd(T.VRe[i], T.VRe[j], T.Nj, -tempRe);
    }
    RealVectorDotVector(T.VRe[i], T.VRe[i], T.Ni, &tempRe);
    tempRe = sqrt(tempRe);
    H->b[i] *= tempRe;
    if (tempRe < 1000 * DBL_EPSILON)
    {
      RealVectorMultiply(T.VRe[i], T.Ni, 0);
      H->N = i;
    }
    else
    {
      RealVectorMultiply(T.VRe[i], T.Ni, 1.0 / tempRe);
    }
  }
  H->a[N] = 0;
  for (i = 0; i < N; i++)
  {
    H->a[N] += T.VRe[N - 1][i] * val[i] * T.VRe[N - 1][i];
  }
  FreeMatrix(&T);
  return 0;
}

/*
 *
 */
uint32_t AndersonMatrixToBADoubleTriDiagonalMatrixSimple(AndersonMatrixType const M,
                                                         BADoubleTriDiagonalMatrixType * const A)
{
  TriDiagonalMatrixType MTri;
  SecondDoubleTriDiagonalMatrixType MSecTri;
  AndersonMatrixToTriDiagonalMatrix(M, &MTri, M.N);
  TriDiagonalMatrixToSecondDoubleTriDiagonalMatrix(MTri, &MSecTri, M.mu, 0.0000001);
  MapSecondoubleTridiagonalMatrixToBA(MSecTri, A);
  FreeTriDiagonalMatrix(&MTri);
  FreeDoubleTriDiagonalMatrix(&(MSecTri.Dt));
  return 0;
}



/*
 *
 */
uint32_t AndersonMatrixToBADoubleTriDiagonalMatrix(AndersonMatrixType const M, double const onsiteshift,
                                                   BADoubleTriDiagonalMatrixType * const A, double *const occupation)
{
  uint32_t i, j, nVal;
  MatrixType Eigenfun;
  double Eigenval[M.N], psival[M.N], psicon[M.N];
  if (M.N == 0)
  {
    printf("AndersonMatrixToBADoubleTriDiagonalMatrix called with a matrix of size 0\n");
    *occupation = 0;
    A->val.N = 0;
    A->con.N = 0;
    InitBADoubleTriDiagonalMatrix(A);
    strcpy(A->Name, M.Name);
    A->a0 = M.a[0];
    A->b0 = 0;
    return 0;
  }
  /* find eigenbasis of the Anderson matrix */
  AndersonMatrixToEigensystemOnsiteShift(M, onsiteshift, &Eigenfun, Eigenval);
  for (i = 0; i < Eigenfun.Ni; i++)
  {
    if (Eigenfun.VRe[i][0] < 0)
    {
      for (j = 0; j < Eigenfun.Nj; j++)
      {
        Eigenfun.VRe[i][j] *= -1;
      }
    }
  }
  nVal = 0;
  *occupation = 0;
  while ((nVal < Eigenfun.Ni) && (Eigenval[nVal] <= M.mu))
  {
    *occupation += Sqr(Eigenfun.VRe[nVal][0]);
    nVal++;
  }
  if (*occupation <= DBL_EPSILON)
  {
    AndersonMatrixToTriDiagonalMatrix(M, &(A->con), M.N);
    A->val.N = 2;
    InitTriDiagonalMatrix(&(A->val));
    A->val.a[0] = 0;
    A->val.b[0] = 0;
    A->val.a[1] = -Abs(A->con.a[1]) - 2 * Abs(A->con.b[1]) - 100 * DBL_EPSILON;
    A->val.b[1] = 0;
    A->val.a[2] = -Abs(A->con.a[1]) - 2 * Abs(A->con.b[1]) - 100 * DBL_EPSILON;
    A->a0 = M.a[0];
    A->b0 = M.b[0];
    A->mu = M.mu;
    strcpy(A->Name, M.Name);
    SetBADoubleTriDiagonalMatrixConVal(A);
  }
  else if (*occupation >= (1 - DBL_EPSILON))
  {
    AndersonMatrixToTriDiagonalMatrix(M, &(A->val), M.N);
    A->con.N = 2;
    InitTriDiagonalMatrix(&(A->con));
    A->con.a[0] = 0;
    A->con.b[0] = 0;
    A->con.a[1] = -Abs(A->val.a[1]) - 2 * Abs(A->val.b[1]) - 100 * DBL_EPSILON;
    A->con.b[1] = 0;
    A->con.a[2] = -Abs(A->val.a[1]) - 2 * Abs(A->val.b[1]) - 100 * DBL_EPSILON;
    A->a0 = M.a[0];
    A->b0 = M.b[0];
    A->mu = M.mu;
    strcpy(A->Name, M.Name);
    SetBADoubleTriDiagonalMatrixConVal(A);
  }
  else
  {
    printf
      ("In AndersonMatrixToBADoubleTriDiagonalMatrix:\n Impurity occupation is %22.15E for shift %22.15E at chemical potential %22.15E\n",
       *occupation, onsiteshift, M.mu);
    if (nVal == 0)
    {
      printf("no states found in valence band; lowest eigen-state has energy %22.15E\n", Eigenval[0]);
    }
    if (nVal == Eigenfun.Ni)
    {
      printf("no states found in conduction band; highest eigen-state has energy %22.15E\n", Eigenval[Eigenfun.Ni - 1]);
    }
    /* create psival and psicon */
    for (i = 0; i < nVal; i++)
    {
      psival[i] = sqrt(1.0 / *occupation) * Eigenfun.VRe[i][0];
      psicon[i] = 0;
    }
    for (i = nVal; i < Eigenfun.Ni; i++)
    {
      psival[i] = 0;
      psicon[i] = sqrt(1.0 / (1.0 - *occupation)) * Eigenfun.VRe[i][0];
    }
    /* create A->val and A->con */
    EigenstateToTridiagonalMatrix(psival, &(Eigenval[0]), Eigenfun.Ni, &(A->val));
    A->val.b[0] = sqrt(*occupation);
    EigenstateToTridiagonalMatrix(psicon, &(Eigenval[0]), Eigenfun.Ni, &(A->con));
    A->con.b[0] = sqrt(1.0 - *occupation);
    A->a0 = M.a[0];
    A->b0 = M.b[0];
    A->mu = M.mu;
    strcpy(A->Name, M.Name);
    SetBADoubleTriDiagonalMatrixConVal(A);
    A->acvB -= onsiteshift;
  }
  FreeMatrix(&Eigenfun);
  return 0;
}

/*
 *
 */
uint32_t AndersonMatrixToBADoubleTriDiagonalMatrixOld(AndersonMatrixType const M, double const onsiteshift,
                                                      BADoubleTriDiagonalMatrixType * const A, double *const occupation)
{
  uint32_t i, j, nVal;
  MatrixType Eigenfun;
  double Eigenval[M.N], psival[M.N], psicon[M.N], psib[M.N];
  if (M.N == 0)
  {
    printf("AndersonMatrixToBADoubleTriDiagonalMatrix called with a matrix of size 0\n");
    *occupation = 0;
    A->val.N = 0;
    A->con.N = 0;
    InitBADoubleTriDiagonalMatrix(A);
    strcpy(A->Name, M.Name);
    A->a0 = M.a[0];
    A->b0 = 0;
    return 0;
  }
  /* find eigenbasis of the Anderson matrix */
  AndersonMatrixToEigensystemOnsiteShift(M, onsiteshift, &Eigenfun, Eigenval);
  for (i = 0; i < Eigenfun.Ni; i++)
  {
    if (Eigenfun.VRe[i][0] < 0)
    {
      for (j = 0; j < Eigenfun.Nj; j++)
      {
        Eigenfun.VRe[i][j] *= -1;
      }
    }
  }
  nVal = 0;
  *occupation = 0;
  while ((nVal < Eigenfun.Ni) && (Eigenval[nVal] <= M.mu))
  {
    *occupation += Sqr(Eigenfun.VRe[nVal][0]);
    nVal++;
  }
  if (*occupation <= DBL_EPSILON)
  {
    AndersonMatrixToTriDiagonalMatrix(M, &(A->con), M.N);
    A->val.N = 2;
    InitTriDiagonalMatrix(&(A->val));
    A->val.a[0] = 0;
    A->val.b[0] = 0;
    A->val.a[1] = -Abs(A->con.a[1]) - 2 * Abs(A->con.b[1]) - 100 * DBL_EPSILON;
    A->val.b[1] = 0;
    A->val.a[2] = -Abs(A->con.a[1]) - 2 * Abs(A->con.b[1]) - 100 * DBL_EPSILON;
    A->a0 = M.a[0];
    A->b0 = M.b[0];
    A->mu = M.mu;
    strcpy(A->Name, M.Name);
    SetBADoubleTriDiagonalMatrixConVal(A);
  }
  else if (*occupation >= (1 - DBL_EPSILON))
  {
    AndersonMatrixToTriDiagonalMatrix(M, &(A->val), M.N);
    A->con.N = 2;
    InitTriDiagonalMatrix(&(A->con));
    A->con.a[0] = 0;
    A->con.b[0] = 0;
    A->con.a[1] = -Abs(A->val.a[1]) - 2 * Abs(A->val.b[1]) - 100 * DBL_EPSILON;
    A->con.b[1] = 0;
    A->con.a[2] = -Abs(A->val.a[1]) - 2 * Abs(A->val.b[1]) - 100 * DBL_EPSILON;
    A->a0 = M.a[0];
    A->b0 = M.b[0];
    A->mu = M.mu;
    strcpy(A->Name, M.Name);
    SetBADoubleTriDiagonalMatrixConVal(A);
  }
  else
  {
    printf
      ("In AndersonMatrixToBADoubleTriDiagonalMatrix:\n Impurity occupation is %22.15E for shift %22.15E at chemical potential %22.15E\n",
       *occupation, onsiteshift, M.mu);
    if (nVal == 0)
    {
      printf("no states found in valence band; lowest eigen-state has energy %22.15E\n", Eigenval[0]);
    }
    if (nVal == Eigenfun.Ni)
    {
      printf("no states found in conduction band; highest eigen-state has energy %22.15E\n", Eigenval[Eigenfun.Ni - 1]);
    }
    /* create psib */
    for (i = 0; i < nVal; i++)
    {
      psib[i] = sqrt((1.0 - *occupation) / *occupation) * Eigenfun.VRe[i][0];
    }
    for (i = nVal; i < Eigenfun.Ni; i++)
    {
      psib[i] = -sqrt(*occupation / (1.0 - *occupation)) * Eigenfun.VRe[i][0];
    }
    /* create psival and psicon */
    for (i = 0; i < Eigenfun.Ni; i++)
    {
      psival[i] = sqrt(*occupation) * Eigenfun.VRe[i][0] + sqrt(1.0 - *occupation) * psib[i];
      psicon[i] = sqrt(1.0 - *occupation) * Eigenfun.VRe[i][0] - sqrt(*occupation) * psib[i];
    }
    /* create A->val and A->con */
    EigenstateToTridiagonalMatrix(psival, &(Eigenval[0]), Eigenfun.Ni, &(A->val));
    A->val.b[0] = sqrt(*occupation);
    EigenstateToTridiagonalMatrix(psicon, &(Eigenval[0]), Eigenfun.Ni, &(A->con));
    A->con.b[0] = sqrt(1.0 - *occupation);
    A->a0 = M.a[0];
    A->b0 = M.b[0];
    A->mu = M.mu;
    strcpy(A->Name, M.Name);
    SetBADoubleTriDiagonalMatrixConVal(A);
    A->acvB -= onsiteshift;
  }
  FreeMatrix(&Eigenfun);
  return 0;
}



/*
 *  Given GBath in anderson form, return GBath in BA form.
 *
 *  an onsite energy (onsiteshift) is added to the impurity site
 *  the one particle Hamiltonian is diagonalized
 *  using the chemical potential M.mu and temperature kbT the occupation is determined
 *  a valence and conduction function is defined such that the impurity function is a linear combinaiton of the two.
 *  a tridiagonal form of the valence and conduction states is created starting from the valence and conduction functions.
 *  The valence and conduction functions are linear combined to the impurity and extra site in the BA form.
 */
uint32_t AndersonMatrixToBADoubleTriDiagonalMatrixTemperature(AndersonMatrixType const M, double const onsiteshift,
                                                              double const kbT, BADoubleTriDiagonalMatrixType * const A,
                                                              double *const occupation)
{
  uint32_t i, j;
  MatrixType Eigenfun, Eigenfundbl;
  double Eigenval[2 * M.N], psival[2 * M.N], psicon[2 * M.N];
  if (M.N == 0)
  {
    printf("AndersonMatrixToBADoubleTriDiagonalMatrix called with a matrix of size 0\n");
    *occupation = 0;
    A->val.N = 0;
    A->con.N = 0;
    InitBADoubleTriDiagonalMatrix(A);
    strcpy(A->Name, M.Name);
    A->a0 = M.a[0];
    A->b0 = 0;
    return 0;
  }
  /* find eigenbasis of the Anderson matrix */
  AndersonMatrixToEigensystemOnsiteShift(M, onsiteshift, &Eigenfun, Eigenval);
  for (i = 0; i < Eigenfun.Ni; i++)
  {
    if (Eigenfun.VRe[i][0] < 0)
    {
      for (j = 0; j < Eigenfun.Nj; j++)
      {
        Eigenfun.VRe[i][j] *= -1;
      }
    }
  }
  *occupation = 0;
  for (i = 0; i < Eigenfun.Ni; i++)
  {
    *occupation += Sqr(Eigenfun.VRe[i][0]) * FermiFunction(Eigenval[i], M.mu, kbT);
  }
  if (*occupation <= 10 * DBL_EPSILON)
  {
    AndersonMatrixToTriDiagonalMatrix(M, &(A->con), M.N);
    A->val.N = 2;
    InitTriDiagonalMatrix(&(A->val));
    A->val.a[0] = 0;
    A->val.b[0] = 0;
    A->val.a[1] = -Abs(A->con.a[1]) - 2 * Abs(A->con.b[1]) - 100 * DBL_EPSILON;
    A->val.b[1] = 0;
    A->val.a[2] = -Abs(A->con.a[1]) - 2 * Abs(A->con.b[1]) - 100 * DBL_EPSILON;
    A->a0 = M.a[0];
    A->b0 = M.b[0];
    A->mu = M.mu;
    strcpy(A->Name, M.Name);
    SetBADoubleTriDiagonalMatrixConVal(A);
  }
  else if (*occupation >= (1 - 10 * DBL_EPSILON))
  {
    AndersonMatrixToTriDiagonalMatrix(M, &(A->val), M.N);
    A->con.N = 2;
    InitTriDiagonalMatrix(&(A->con));
    A->con.a[0] = 0;
    A->con.b[0] = 0;
    A->con.a[1] = -Abs(A->val.a[1]) - 2 * Abs(A->val.b[1]) - 100 * DBL_EPSILON;
    A->con.b[1] = 0;
    A->con.a[2] = -Abs(A->val.a[1]) - 2 * Abs(A->val.b[1]) - 100 * DBL_EPSILON;
    A->a0 = M.a[0];
    A->b0 = M.b[0];
    A->mu = M.mu;
    strcpy(A->Name, M.Name);
    SetBADoubleTriDiagonalMatrixConVal(A);
  }
  else
  {
    printf
      ("In AndersonMatrixToBADoubleTriDiagonalMatrix:\n Impurity occupation is %22.15E for shift %22.15E at chemical potential %22.15E\n",
       *occupation, onsiteshift, M.mu);
    /* double Eigenfun and Eigenval */
    Eigenfundbl.Ni = 2 * Eigenfun.Ni;
    Eigenfundbl.Nj = 2 * Eigenfun.Nj;
    Eigenfundbl.QComplex = 0; // real
    InitMatrix(&Eigenfundbl);
    for (i = 0; i < Eigenfun.Ni; i++)
    {
      Eigenval[i + Eigenfun.Ni] = Eigenval[i];
      for (j = 0; j < Eigenfun.Nj; j++)
      {
        Eigenfundbl.VRe[i][j] = Eigenfun.VRe[i][j] * sqrt(FermiFunction(Eigenval[i], M.mu, kbT));
        Eigenfundbl.VRe[i][j + Eigenfun.Nj] = 0;
        Eigenfundbl.VRe[i + Eigenfun.Ni][j] = Eigenfun.VRe[i][j] * sqrt(1 - FermiFunction(Eigenval[i], M.mu, kbT));
        Eigenfundbl.VRe[i + Eigenfun.Ni][j + Eigenfun.Nj] = 0;
      }
      Eigenfundbl.VRe[i][i + Eigenfun.Nj] = sqrt(1 - FermiFunction(Eigenval[i], M.mu, kbT));
      Eigenfundbl.VRe[i + Eigenfun.Ni][i + Eigenfun.Nj] = -sqrt(FermiFunction(Eigenval[i], M.mu, kbT));
    }

    /* create psival and psicon */
    for (i = 0; i < Eigenfun.Ni; i++)
    {
      psival[i] = sqrt(1.0 / *occupation) * Eigenfundbl.VRe[i][0];
      psival[i + Eigenfun.Ni] = 0;
      psicon[i + Eigenfun.Ni] = sqrt(1.0 / (1 - *occupation)) * Eigenfundbl.VRe[i + Eigenfun.Ni][0];
      psicon[i] = 0;
    }
    /* create A->val and A->con */
    EigenstateToTridiagonalMatrix(psival, &(Eigenval[0]), Eigenfundbl.Ni, &(A->val));
    A->val.b[0] = sqrt(*occupation);
    EigenstateToTridiagonalMatrix(psicon, &(Eigenval[0]), Eigenfundbl.Ni, &(A->con));
    A->con.b[0] = sqrt(1.0 - *occupation);
    A->a0 = M.a[0];
    A->b0 = M.b[0];
    A->mu = M.mu;
    strcpy(A->Name, M.Name);
    SetBADoubleTriDiagonalMatrixConVal(A);
    A->acvB -= onsiteshift;
  }
  FreeMatrix(&Eigenfun);
  return 0;
}

/*
 *  Determine mu such that the occupation of the first state in Eigenfun is equal to Nelec
 */
uint32_t DetermineMuFromEigensystem(MatrixType const Eigenfun, double const *Eigenval, double const Nelec,
                                    double const kbT, double *const mu)
{
  uint32_t const maxitter = 30;
  uint32_t i, converged = 0, itter = 0;
  double occupation[maxitter + 1], mulist[maxitter + 1];
  occupation[0] = 0;
  if ((Nelec <= 0) || (Nelec >= 1))
  {
    return 0;
  }
  mulist[0] = *mu;
  while (!converged)
  {
    for (i = 0; i < Eigenfun.Ni; i++)
    {
      occupation[itter] += Sqr(Eigenfun.VRe[i][0]) * FermiFunction(Eigenval[i], mulist[itter], kbT);
    }
    if (Abs(occupation[itter] - Nelec) < 100 * DBL_EPSILON)
    {
      converged = 1;
    }
    if (itter == 0)
    {
      mulist[1] = *mu + 0.1;
    }
    else
    {
      mulist[itter + 1] = InterpolateRootOfset(mulist, occupation, itter + 1, Nelec);
    }
    if (itter == maxitter)
    {
      printf("in DetermineMuFromEigensystem reached maximum itterations\n");

      converged = 1;
    }
    else if (!converged)
    {
      itter++;
    }
  }
  return 0;
}

/*
 *
 */
uint32_t AndersonMatrixToBADoubleTriDiagonalMatrixWrong(AndersonMatrixType const M,
                                                        BADoubleTriDiagonalMatrixType * const A, uint32_t NTriDiagVal,
                                                        uint32_t NTriDiagCon, double Ef, double Epsilon)
{
  ListOfPoles LP;
  if (AndersonMatrixToListOfPoles(M, &LP))
  {
    printf("AndersonMatrixToListOfPoles failed in AndersonMatrixToBADoubleTriDiagonalMatrix\n");
    return 1;
  };
  if (ListOfPolesToBADoubleTridiagonalMatrix(LP, A, NTriDiagVal, NTriDiagCon, Ef, Epsilon))
  {
    printf("ListOfPolesToBADoubleTridiagonalMatrix failed in AndersonMatrixToBADoubleTriDiagonalMatrix\n");
    return 1;
  };
  FreeListOfPoles(&LP);
  strcpy(A->Name, "GBath");
  return 0;
}

/****************************************************************************/
/********************         Matrix conversions         ********************/
/********************         BlokAndersonMatrix To ...  ********************/
/****************************************************************************/






/****************************************************************************/
/********************         Matrix conversions         ********************/
/********************         TriDiagonalMatrix To ...   ********************/
/****************************************************************************/

uint32_t TriDiagonalMatrixToEigenvalues(TriDiagonalMatrixType const M, double *val)
{
  int    i, info;
  const char *jobz = "N";
  double e[M.N];
  int    lwork = (1 + 4 * M.N + Sqr(M.N));
  int    liwork = 3 + 5 * M.N;
  double *work,*Z=NULL;
  int   *iwork;
  work = malloc(lwork * sizeof(double));
  if (work == NULL)
  {
    printf("malloc of work failed in TriDiagonalMatrixToEigensystem\n");
    return 1;
  };
  iwork = malloc(liwork * sizeof(int));
  if (iwork == NULL)
  {
    printf("malloc of iwork failed in TriDiagonalMatrixToEigensystem\n");
    return 1;
  };
  
  for (i = 0; i < M.N - 1; i++) e[i] = M.b[i + 1];
  e[M.N - 1] = 0.0;
  for (i = 0; i < M.N; i++) val[i] = M.a[i + 1];
  
  dstevd(jobz, &M.N, val, e, Z, &M.N, work, &lwork, iwork, &liwork, &info);
  if (info != 0)
  {
    printf("dstevd failed in TriDiagonalMatrixToEigensystem with info value %i\n", info);
    printf("While working on matrix %s\n", M.Name);
    printf("continue diagonalization with dstev\n");
    fflush(stdout);
    dstev(jobz, &M.N, val, e, Z, &M.N, work, &info);
    if (info != 0)
    {
      printf("Also dstev failed with info value %i\n", info);
      printf("I don't know what to do\n");
      printf("ERROR!!!!!!!!!!!\n");
      fflush(stdout);
      return 1;
    }
  }
  free(work);
  free(iwork);
  return 0;
  
}


uint32_t TriDiagonalMatrixToEigensystem(TriDiagonalMatrixType const M, MatrixType * const T, double *val)
{
  int    i, info;
  const char *jobz = "V";
  double e[M.N];
  int    lwork = (1 + 4 * M.N + Sqr(M.N));
  int    liwork = 3 + 5 * M.N;
  double *work;
  int   *iwork;
  work = malloc(lwork * sizeof(double));
  if (work == NULL)
  {
    printf("malloc of work failed in TriDiagonalMatrixToEigensystem\n");
    return 1;
  };
  iwork = malloc(liwork * sizeof(int));
  if (iwork == NULL)
  {
    printf("malloc of iwork failed in TriDiagonalMatrixToEigensystem\n");
    return 1;
  };

  for (i = 0; i < M.N - 1; i++)
    e[i] = M.b[i + 1];
  e[M.N - 1] = 0.0;
  for (i = 0; i < M.N; i++)
    val[i] = M.a[i + 1];

  T->Ni = M.N;
  T->Nj = M.N;
  T->QComplex = 0;
  strcpy(T->Name, "Eigenvectors of ");
  strncat(T->Name, M.Name, Max(0, MAXStrLength - 18));
  if (InitMatrix(T))
  {
    printf("InitMatrix(T) failed in TriDiagonalMatrixToEigensystem\n");
    fflush(stdout);
    free(work);
    free(iwork);
    return 1;
  }

  dstevd(jobz, &M.N, val, e, T->VRe[0], &M.N, work, &lwork, iwork, &liwork, &info);
  if (info != 0)
  {
    printf("dstevd failed in TriDiagonalMatrixToEigensystem with info value %i\n", info);
    printf("While working on matrix %s\n", M.Name);
    printf("continue diagonalization with dstev\n");
    fflush(stdout);
    dstev(jobz, &M.N, val, e, T->VRe[0], &M.N, work, &info);
    if (info != 0)
    {
      printf("Also dstev failed with info value %i\n", info);
      printf("I don't know what to do\n");
      printf("ERROR!!!!!!!!!!!\n");
      fflush(stdout);
      return 1;
    }
  }
  free(work);
  free(iwork);
  return 0;
}

uint32_t TriDiagonalMatrixToListOfPoles(TriDiagonalMatrixType const M, ListOfPoles * const A)
{
  MatrixType T;
  uint32_t i;
  A->N = M.N;
  if (InitListOfPoles(A))
  {
    printf("InitListOfPoles failed in TriDiagonalMatrixToListOfPoles\n");
    fflush(stdout);
    return 1;
  };
  if (TriDiagonalMatrixToEigensystem(M, &T, A->a))
  {
    printf("TriDiagonalMatrixToEigensystem failed in TriDiagonalMatrixToListOfPoles\n");
    fflush(stdout);
    return 1;
  };
  for (i = 0; i < M.N; i++)
  {
    A->b[i] = Sqr(M.b[0] * T.VRe[i][0]);
  }
  FreeMatrix(&T);
  A->a0 = M.a[0];
  A->mu = M.mu;
  strcpy(A->Name, M.Name);
  return 0;
}

uint32_t TriDiagonalMatrixToAndersonMatrix(TriDiagonalMatrixType const M, AndersonMatrixType * const A)
{
  TriDiagonalMatrixType Mminor;
  MatrixType Tminor;
  uint32_t i;
  A->N = M.N;
  InitAndersonMatrix(A);
  if (M.N <= 2)
  {
    for (i = 0; i < M.N; i++)
    {
      A->a[i] = M.a[i];
      A->b[i] = M.b[i];
    }
    A->a[M.N] = M.a[M.N];
    return 0;
  }
  A->a[0] = M.a[0];
  A->b[0] = M.b[0];
  A->a[1] = M.a[1];
  TriDiagonalMatrixPart(&Mminor, M, 2, M.N);
  TriDiagonalMatrixToEigensystem(Mminor, &Tminor, &(A->a[2]));
  for (i = 1; i < M.N; i++)
  {
    A->b[i] = fabs(M.b[1] * Tminor.VRe[i - 1][0]);
  }
  FreeMatrix(&Tminor);
  FreeTriDiagonalMatrix(&Mminor);
  return 0;
}


uint32_t TriDiagonalMatrixToAndersonMatrixPlusT(TriDiagonalMatrixType const M, MatrixType * const T,
                                                AndersonMatrixType * const A)
{
  TriDiagonalMatrixType Mminor;
  MatrixType Tminor;
  uint32_t i, j;
  int    sign;
  A->N = M.N;
  InitAndersonMatrix(A);
  T->QComplex = 0;
  T->Ni = M.N;
  T->Nj = M.N;
  InitIdentityMatrix(T);
  strcpy(T->Name, "T_Matrix");
  if (M.N <= 2)
  {
    for (i = 0; i < M.N; i++)
    {
      A->a[i] = M.a[i];
      A->b[i] = M.b[i];
    }
    A->a[M.N] = M.a[M.N];
    return 0;
  }
  A->a[0] = M.a[0];
  A->b[0] = M.b[0];
  A->a[1] = M.a[1];
  TriDiagonalMatrixPart(&Mminor, M, 2, M.N);
  TriDiagonalMatrixToEigensystem(Mminor, &Tminor, &(A->a[2]));
  for (i = 1; i < M.N; i++)
  {
    A->b[i] = fabs(M.b[1] * Tminor.VRe[i - 1][0]);
  }
  for (i = 1; i < M.N; i++)
  {
    sign = ((M.b[1] * Tminor.VRe[i - 1][0]) > 0) - ((M.b[1] * Tminor.VRe[i - 1][0]) < 0);
    for (j = 1; j < M.N; j++)
    {
      T->VRe[i][j] = sign * Tminor.VRe[i - 1][j - 1];
    }
  }
  FreeMatrix(&Tminor);
  FreeTriDiagonalMatrix(&Mminor);
  return 0;
}

uint32_t TriDiagonalMatrixToDoubleTriDiagonalMatrix(TriDiagonalMatrixType const M,
                                                    DoubleTriDiagonalMatrixType * const A, double Ef, double epsilon)
{
  AndersonMatrixType An;
  if (TriDiagonalMatrixToAndersonMatrix(M, &An))
  {
    printf("TriDiagonalMatrixToAndersonMatrix failed in TriDiagonalMatrixToDoubleTriDiagonalMatrix\n");
    fflush(stdout);
    return 1;
  };
  if (AndersonMatrixToDoubleTriDiagonalMatrix(An, A, Ef, epsilon))
  {
    printf("AndersonMatrixToDoubleTriDiagonalMatrix failed in TriDiagonalMatrixToDoubleTriDiagonalMatrix\n");
    fflush(stdout);
    return 1;
  };
  FreeAndersonMatrix(&An);
  return 0;
}


uint32_t TriDiagonalMatrixToSecondDoubleTriDiagonalMatrix(TriDiagonalMatrixType const M,
                                                          SecondDoubleTriDiagonalMatrixType * const A, double Ef,
                                                          double epsilon)
{
  TriDiagonalMatrixType SecTri;
  if (M.N < 2)
  {
    printf("Help error in TriDiagonalMatrixToSecondDoubleTriDiagonalMatrix catch not programmed yet\n");
  };
  A->a0 = M.a[0];
  A->b0 = M.b[0];
  A->mu = M.mu;
  SecTri.N = M.N - 1;
  SecTri.a = &(M.a[1]);
  SecTri.b = &(M.b[1]);
  return TriDiagonalMatrixToDoubleTriDiagonalMatrix(SecTri, &(A->Dt), Ef, epsilon);
}

uint32_t TriDiagonalMatrixToDoubleTriDiagonalMatrixPlusT(TriDiagonalMatrixType const M, MatrixType * const T,
                                                         DoubleTriDiagonalMatrixType * const A, double Ef,
                                                         double epsilon)
{
  AndersonMatrixType An;
  MatrixType Tone, Ttwo;
  if (TriDiagonalMatrixToAndersonMatrixPlusT(M, &Tone, &An))
  {
    printf("TriDiagonalMatrixToAndersonMatrix failed in TriDiagonalMatrixToDoubleTriDiagonalMatrix\n");
    fflush(stdout);
    return 1;
  };
  if (AndersonMatrixToDoubleTriDiagonalMatrixPlusT(An, &Ttwo, A, Ef, epsilon))
  {
    printf("AndersonMatrixToDoubleTriDiagonalMatrix failed in TriDiagonalMatrixToDoubleTriDiagonalMatrix\n");
    fflush(stdout);
    return 1;
  };
  FreeAndersonMatrix(&An);
  if (MatrixMultiply(T, Ttwo, Tone))
  {
    printf("MatrixMultiply failed in TriDiagonalMatrixToDoubleTriDiagonalMatrix\n");
    fflush(stdout);
    return 1;
  };
  FreeMatrix(&Tone);
  FreeMatrix(&Ttwo);
  return 0;
}


uint32_t TriDiagonalMatrixToSecondDoubleTriDiagonalMatrixPlusT(TriDiagonalMatrixType const M, MatrixType * const T,
                                                               SecondDoubleTriDiagonalMatrixType * const A, double Ef,
                                                               double epsilon)
{
  TriDiagonalMatrixType SecTri;
  if (M.N < 2)
  {
    printf("Help error in TriDiagonalMatrixToSecondDoubleTriDiagonalMatrix catch not programmed yet\n");
  };
  A->a0 = M.a[0];
  A->b0 = M.b[0];
  SecTri.N = M.N - 1;
  SecTri.a = &(M.a[1]);
  SecTri.b = &(M.b[1]);
  return TriDiagonalMatrixToDoubleTriDiagonalMatrixPlusT(SecTri, T, &(A->Dt), Ef, epsilon);
}


uint32_t TriDiagonalMatrixToMatrix(TriDiagonalMatrixType const M, MatrixType * const A)
{
  uint32_t i;
  A->QComplex = 0;
  A->Ni = M.N;
  A->Nj = M.N;
  if (InitMatrix(A))
  {
    printf("InitMatrix failed in TriDiagonalMatrixToMatrix\n");
    return 1;
  }
  A->a0 = M.a[0];
  A->b0 = M.b[0];
  for (i = 1; i < M.N; i++)
  {
    A->VRe[i - 1][i - 1] = M.a[i];
    A->VRe[i][i - 1] = M.b[i];
    A->VRe[i - 1][i] = M.b[i];
  }
  A->VRe[M.N - 1][M.N - 1] = M.a[M.N];
  return 0;
}

/*
 *
 */
uint32_t TriDiagonalMatrixToBADoubleTriDiagonalMatrix(TriDiagonalMatrixType const M,
                                                      BADoubleTriDiagonalMatrixType * const A, uint32_t NTriDiagVal,
                                                      uint32_t NTriDiagCon, double Ef, double Epsilon)
{
  ListOfPoles LP;
  if (TriDiagonalMatrixToListOfPoles(M, &LP))
  {
    printf("TriDiagonalMatrixToListOfPoles failed in TriDiagonalMatrixToBADoubleTriDiagonalMatrix\n");
    return 1;
  };
  if (ListOfPolesToBADoubleTridiagonalMatrix(LP, A, NTriDiagVal, NTriDiagCon, Ef, Epsilon))
  {
    printf("ListOfPolesToBADoubleTridiagonalMatrix failed in TriDiagonalMatrixToBADoubleTriDiagonalMatrix\n");
    return 1;
  };
  FreeListOfPoles(&LP);
  strcpy(A->Name, M.Name);
  return 0;
}


/*
 *
 */
uint32_t TriDiagonalMatrixToBADoubleTriDiagonalMatrixOccupation(TriDiagonalMatrixType const M,
                                                                BADoubleTriDiagonalMatrixType * const A,
                                                                uint32_t NTriDiagVal, uint32_t NTriDiagCon, double Ef,
                                                                double Occupation)
{
  ListOfPoles LP;
  if (TriDiagonalMatrixToListOfPoles(M, &LP))
  {
    printf("TriDiagonalMatrixToListOfPoles failed in TriDiagonalMatrixToBADoubleTriDiagonalMatrix\n");
    return 1;
  };
  if (ListOfPolesToBADoubleTridiagonalMatrixOccupation(LP, A, NTriDiagVal, NTriDiagCon, Ef, Occupation))
  {
    printf("ListOfPolesToBADoubleTridiagonalMatrix failed in TriDiagonalMatrixToBADoubleTriDiagonalMatrix\n");
    return 1;
  };
  FreeListOfPoles(&LP);
  strcpy(A->Name, M.Name);
  return 0;
}

/******************************************************************************/
/********************         Matrix conversions           ********************/
/********************         BlokTriDiagonalMatrix To ... ********************/
/******************************************************************************/






/********************************************************************************/
/********************         Matrix conversions             ********************/
/********************         DoubleTriDiagonalMatrix To ... ********************/
/********************************************************************************/

uint32_t DoubleTriDiagonalMatrixToTriDiagonalMatrix(DoubleTriDiagonalMatrixType const M,
                                                    TriDiagonalMatrixType * const H, uint32_t NTriDiag)
{
  uint32_t i, j;
  double tempRe;
  MatrixType T;
  T.QComplex = 0;
  T.Ni = NTriDiag;
  T.Nj = M.val.N + M.con.N + 1;
  InitMatrix(&T);
  H->N = NTriDiag;
  if (InitTriDiagonalMatrix(H))
  {
    FreeMatrix(&T);
    printf("InitTriDiagonalMatrix failed in DoubleTriDiagonalMatrixToTriDiagonalMatrix\n");
    return 1;
  }
  strcpy(T.Name, "TriDiagonalizer of ");
  strncat(T.Name, M.Name, Max(0, MAXStrLength - 20));
  strcpy(H->Name, "TriDiagonal form of ");
  strncat(H->Name, M.Name, Max(0, MAXStrLength - 22));
  H->a[0] = M.a0;
  H->b[0] = M.b0;
  for (i = 0; i < M.val.N + M.con.N + 1; i++)
  {
    T.VRe[0][i] = 0;
  }
  T.VRe[0][0] = 1;
  for (i = 1; i < NTriDiag; i++)
  {
    RealDoubleTriDiagonalMatrixVector(M, T.VRe[i - 1], T.VRe[i]);
    for (j = 0; j < i - 1; j++)
    {
      RealVectorDotVector(T.VRe[j], T.VRe[i], T.Nj, &tempRe);
      RealVectorAdd(T.VRe[i], T.VRe[j], T.Nj, -tempRe);
    }
    RealVectorDotVector(T.VRe[i - 1], T.VRe[i], T.Nj, &(H->a[i]));
    RealVectorAdd(T.VRe[i], T.VRe[i - 1], T.Nj, -H->a[i]);
    RealVectorDotVector(T.VRe[i], T.VRe[i], T.Nj, &(H->b[i]));
    H->b[i] = sqrt(H->b[i]);
    if (H->b[i] > 100 * DBL_MIN)
    {
      RealVectorMultiply(T.VRe[i], T.Nj, 1.0 / H->b[i]);
      for (j = 0; j < i; j++)
      {
        RealVectorDotVector(T.VRe[j], T.VRe[i], T.Nj, &tempRe);
        RealVectorAdd(T.VRe[i], T.VRe[j], T.Nj, -tempRe);
      }
      RealVectorDotVector(T.VRe[i], T.VRe[i], T.Nj, &tempRe);
      tempRe = sqrt(tempRe);
      H->b[i] *= tempRe;
      if (tempRe > 100 * DBL_MIN)
      {
        RealVectorMultiply(T.VRe[i], T.Nj, 1.0 / tempRe);
      }
      else
      {
        RealVectorMultiply(T.VRe[i], T.Nj, 0);
      }
    }
    else
    {
      RealVectorMultiply(T.VRe[i], T.Nj, 0);
    }
  }
  RealVectorDoubleTriDiagonalMatrixVector(M, T.VRe[NTriDiag - 1], &(H->a[NTriDiag]));
  FreeMatrix(&T);
  return 0;
}



/**************************************************************************************/
/********************         Matrix conversions                   ********************/
/********************         SecondDoubleTriDiagonalMatrix To ... ********************/
/**************************************************************************************/

uint32_t SecondDoubleTriDiagonalMatrixToTriDiagonalMatrix(SecondDoubleTriDiagonalMatrixType const M,
                                                          TriDiagonalMatrixType * const H, uint32_t const NTriDiag)
{
  uint32_t i;
  TriDiagonalMatrixType Hlow;
  if (NTriDiag > 1)
    if (DoubleTriDiagonalMatrixToTriDiagonalMatrix(M.Dt, &Hlow, NTriDiag - 1))
    {
      printf("DoubleTriDiagonalMatrixToTriDiagonalMatrix failed in SecondDoubleTriDiagonalMatrixToTriDiagonalMatrix\n");
      return 1;
    };
  H->N = NTriDiag;
  InitTriDiagonalMatrix(H);
  H->a[0] = M.a0;
  H->b[0] = M.b0;
  H->a[1] = Hlow.a[0];
  for (i = 0; i < NTriDiag - 1; i++)
  {
    H->b[i + 1] = Hlow.b[i];
    H->a[i + 2] = Hlow.a[i + 1];
  }
  FreeTriDiagonalMatrix(&Hlow);
  return 0;
}



/************************************************************************************/
/********************         Matrix conversions                 ********************/
/********************         BlokDoubleTriDiagonalMatrix To ... ********************/
/************************************************************************************/





/****************************************************************************/
/********************         Matrix conversions         ********************/
/********************         SparseMatrix To ...        ********************/
/****************************************************************************/

uint32_t SparseMatrixToTriDiagonalMatrixNonOrthogonal(SparseMatrixType M, TriDiagonalMatrixType * const H,
                                                      uint32_t NTriDiag)
{
  /* to do */
  return 1;
}

/*
 *
 */
uint32_t SparseMatrixToTriDiagonalMatrix(SparseMatrixType M, MatrixType * const T, TriDiagonalMatrixType * const H,
                                         uint32_t NTriDiag)
{
  uint32_t i, j;
  double tempRe;
  if (M.QComplex)
  {
    printf("SparseMatrixToTriDiagonalMatrix not implemented for complex sparse matricess");
    fflush(stdout);
    return 1;
  };
  T->QComplex = M.QComplex;
  T->Ni = NTriDiag;
  T->Nj = M.Nj;
  if (InitMatrix(T))
  {
    printf("InitMatrix failed in SparseMatrixToTriDiagonalMatrix\n");
    fflush(stdout);
    return 1;
  };
  H->N = NTriDiag;
  if (InitTriDiagonalMatrix(H))
  {
    FreeMatrix(T);
    printf("InitTriDiagonalMatrix failed in SparseMatrixToTriDiagonalMatrix\n");
    fflush(stdout);
    return 1;
  }
  strcpy(T->Name, "TriDiagonalizer of ");
  strncat(T->Name, M.Name, Max(0, MAXStrLength - 20));
  strcpy(H->Name, "TriDiagonal form of ");
  strncat(H->Name, M.Name, Max(0, MAXStrLength - 22));
  H->a[0] = M.a0;
  H->b[0] = M.b0;
  for (i = 0; i < T->Nj; i++)
  {
    T->VRe[0][i] = 0;
  }
  T->VRe[0][0] = 1;
  for (i = 1; i < H->N; i++)
  {
    RealSparseMatrixVector(M, T->VRe[i - 1], T->VRe[i]);
    for (j = 0; j < i - 1; j++)
    {
      RealVectorDotVector(T->VRe[j], T->VRe[i], T->Nj, &tempRe);
      RealVectorAdd(T->VRe[i], T->VRe[j], T->Nj, -tempRe);
    }
    RealVectorDotVector(T->VRe[i - 1], T->VRe[i], T->Nj, &(H->a[i]));
    RealVectorAdd(T->VRe[i], T->VRe[i - 1], T->Nj, -H->a[i]);
    RealVectorDotVector(T->VRe[i], T->VRe[i], T->Nj, &(H->b[i]));
    H->b[i] = sqrt(H->b[i]);
    if (H->b[i] < (100 * DBL_MIN + DBL_EPSILON * H->b[1]))
    {
      RealVectorMultiply(T->VRe[i], T->Nj, 0);
      H->N = i;
    }
    else
    {
      RealVectorMultiply(T->VRe[i], T->Nj, 1.0 / H->b[i]);
    }
    for (j = 0; j < i; j++)
    {
      RealVectorDotVector(T->VRe[j], T->VRe[i], T->Nj, &tempRe);
      RealVectorAdd(T->VRe[i], T->VRe[j], T->Nj, -tempRe);
    }
    RealVectorDotVector(T->VRe[i], T->VRe[i], T->Nj, &tempRe);
    tempRe = sqrt(tempRe);
    H->b[i] *= tempRe;
    if (tempRe < (100 * DBL_MIN))
    {
      RealVectorMultiply(T->VRe[i], T->Nj, 0);
    }
    else
    {
      RealVectorMultiply(T->VRe[i], T->Nj, 1.0 / tempRe);
    }
  }
  RealVectorSparseMatrixVector(M, T->VRe[NTriDiag - 1], &(H->a[NTriDiag]));
  return 0;
}


/****************************************************************************/
/********************         Matrix conversions         ********************/
/********************         Matrix To ...              ********************/
/****************************************************************************/

/*
*   U->NBitsKey must be set
*   addvised value 2^NBitsKey ~ number of nonzero elements
*/
uint32_t MatrixToSparseMatrix(MatrixType const M, SparseMatrixType * const U, double const epsilon)
{
  uint32_t i, j;
  U->QComplex = M.QComplex;
  U->Ni = M.Ni;
  U->Nj = M.Nj;
  if (InitSparseMatrix(U))
    return 1;
  if (M.QComplex)
  {                             /* complex */
    for (i = 0; i < M.Ni; i++)
    {
      for (j = 0; j < M.Ni; j++)
      {
        if ((Sqr(M.VRe[i][j]) + Sqr(M.VIm[i][j])) > Sqr(epsilon))
        {
          if (ComplexSparseMatrixAddElement(U, i, j, M.VRe[i][j], M.VIm[i][j]))
          {
            FreeSparseMatrix(U);
            return 1;
          }
        }
      }
    }
  }
  /* Real */
  for (i = 0; i < M.Ni; i++)
  {
    for (j = 0; j < M.Ni; j++)
    {
      if ((Sqr(M.VRe[i][j])) > Sqr(epsilon))
      {
        if (RealSparseMatrixAddElement(U, i, j, M.VRe[i][j]))
        {
          FreeSparseMatrix(U);
          return 1;
        }
      }
    }
  }
  return 0;
}

/*
*  M=T.H.T^{\dag}
*  input M defined, H and T not initiallized
*/
uint32_t TriDiagonalize(MatrixType const M, MatrixType * const T, TriDiagonalMatrixType * const H)
{
  uint32_t i, j;
  double tempRe, tempIm;
  if (M.Ni != M.Nj)
    return 2;
  T->QComplex = M.QComplex;
  T->Ni = M.Ni;
  T->Nj = M.Nj;
  if (InitMatrix(T))
    return 1;
  H->N = M.Ni;
  if (InitTriDiagonalMatrix(H))
  {
    FreeMatrix(T);
    return 1;
  }
  strcpy(T->Name, "TriDiagonalizer of ");
  strncat(T->Name, M.Name, Max(0, MAXStrLength - 20));
  strcpy(H->Name, "TriDiagonal form of ");
  strncat(H->Name, M.Name, Max(0, MAXStrLength - 22));
  H->a[0] = 0;
  H->b[0] = 1;
  if (M.QComplex)
  {                             /* complex tri-diagonalization */
    for (i = 0; i < M.Ni; i++)
    {
      T->VRe[0][i] = 0;
      T->VIm[0][i] = 0;
    }
    T->VRe[0][0] = 1;
    for (i = 1; i < M.Ni; i++)
    {
      ComplexMatrixVector(M, T->VRe[i - 1], T->VIm[i - 1], T->VRe[i], T->VIm[i]); /* T[i]=M T[i-1]      */
      for (j = 0; j < i - 1; j++)
      {
        ComplexVectorDotVector(T->VRe[j], T->VIm[j], T->VRe[i], T->VIm[i], T->Nj, &tempRe, &tempIm);  /* temp=T[j]^* . T[i] */
        ComplexVectorAdd(T->VRe[i], T->VIm[i], T->VRe[j], T->VIm[j], T->Nj, -tempRe, -tempIm);  /* T[i]+=-temp T[j]   */
      }
      ComplexVectorDotVector(T->VRe[i - 1], T->VIm[i - 1], T->VRe[i], T->VIm[i], T->Nj, &(H->a[i]), &tempIm);
      ComplexVectorAdd(T->VRe[i], T->VIm[i], T->VRe[j], T->VIm[j], T->Nj, -H->a[i], -tempIm);
      ComplexVectorDotVector(T->VRe[i], T->VIm[i], T->VRe[i], T->VIm[i], T->Nj, &(H->b[i]), &tempIm);
      H->b[i] = sqrt(H->b[i]);
      ComplexVectorMultiply(T->VRe[i], T->VIm[i], T->Nj, 1.0 / H->b[i], 0);
      for (j = 0; j < i; j++)
      {
        ComplexVectorDotVector(T->VRe[j], T->VIm[j], T->VRe[i], T->VIm[i], T->Nj, &tempRe, &tempIm);
        ComplexVectorAdd(T->VRe[i], T->VIm[i], T->VRe[j], T->VIm[j], T->Nj, -tempRe, -tempIm);
      }
      ComplexVectorDotVector(T->VRe[i], T->VIm[i], T->VRe[i], T->VIm[i], T->Nj, &tempRe, &tempIm);
      tempRe = sqrt(tempRe);
      H->b[i] *= tempRe;
      ComplexVectorMultiply(T->VRe[i], T->VIm[i], T->Nj, 1.0 / tempRe, 0);
    }
    ComplexVectorMatrixVector(M, T->VRe[M.Ni - 1], T->VIm[M.Ni - 1], &(H->a[M.Ni]), &tempIm);
    return 0;
  }
  /* real tri-diagonalization */
  for (i = 0; i < M.Ni; i++)
  {
    T->VRe[0][i] = 0;
  }
  T->VRe[0][0] = 1;
  for (i = 1; i < M.Ni; i++)
  {
    RealMatrixVector(M, T->VRe[i - 1], T->VRe[i]);
    for (j = 0; j < i - 1; j++)
    {
      RealVectorDotVector(T->VRe[j], T->VRe[i], T->Nj, &tempRe);
      RealVectorAdd(T->VRe[i], T->VRe[j], T->Nj, -tempRe);
    }
    RealVectorDotVector(T->VRe[i - 1], T->VRe[i], T->Nj, &(H->a[i]));
    RealVectorAdd(T->VRe[i], T->VRe[i - 1], T->Nj, -H->a[i]);
    RealVectorDotVector(T->VRe[i], T->VRe[i], T->Ni, &(H->b[i]));
    H->b[i] = sqrt(H->b[i]);
    RealVectorMultiply(T->VRe[i], T->Ni, 1.0 / H->b[i]);
    for (j = 0; j < i; j++)
    {
      RealVectorDotVector(T->VRe[j], T->VRe[i], T->Nj, &tempRe);
      RealVectorAdd(T->VRe[i], T->VRe[j], T->Nj, -tempRe);
    }
    RealVectorDotVector(T->VRe[i], T->VRe[i], T->Ni, &tempRe);
    tempRe = sqrt(tempRe);
    H->b[i] *= tempRe;
    RealVectorMultiply(T->VRe[i], T->Ni, 1.0 / tempRe);
  }
  RealVectorMatrixVector(M, T->VRe[M.Ni - 1], &(H->a[M.Ni]));
  return 0;
}

uint32_t MatrixToTriDiagonalMatrix(MatrixType const M, MatrixType * const T, TriDiagonalMatrixType * const H)
{
  return TriDiagonalize(M, T, H);
}

uint32_t MatrixToTriDiagonalMatrixNonOrthogonal(MatrixType const M, TriDiagonalMatrixType * const H)
{
  /* to do */
  return 1;
}

/*
*  M=T^{*}.H.T^{T} with H a diagonal matrix with the values val on the diagonal
*  input M defined, T not initiallized
*/
uint32_t Eigensystem(MatrixType const M, MatrixType * const T, double *val)
{
  int    info;
  const char *jobz = "V";
  const char *uplo = "U";
  int    lwork = 1 + 6 * M.Ni + 2 * Sqr(M.Ni);
  int    liwork = 3 + 5 * M.Ni;
  int    iwork[liwork];
  double work[lwork];
  if (M.Ni != M.Nj)
  {
    printf("Eigensystem called with non-square matrix\n");
    printf("ERORR!!!!!!\n");
    return 1;
  }
  if (M.QComplex == 0)
  {                             /* Real */
    if (MatrixCopy(T, M))
      return 1;
    strcpy(T->Name, "Eigenvectors of ");
    strncat(T->Name, M.Name, Max(0, MAXStrLength - 16));
    dsyev(jobz, uplo, &M.Ni, T->VRe[0], &M.Ni, val, work, &lwork, &info);
    if (info)
    {
      printf("dsyv failed in Eigensystem with info value %i\n", info);
      printf("While working on matrix %s\n", M.Name);
      printf("continue diagonalization with dsyvd\n");
      fflush(stdout);
      FreeMatrix(T);
      if (MatrixCopy(T, M))
        return 1;
      strcpy(T->Name, "Eigenvectors of ");
      strncat(T->Name, M.Name, Max(0, MAXStrLength - 16));
      dsyevd(jobz, uplo, &M.Ni, T->VRe[0], &M.Ni, val, work, &lwork, iwork, &liwork, &info);
      if (info)
      {
        printf("also dsyvd failed in Eigensystem with info value %i\n", info);
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
    Complex16 MReIm[M.Ni * M.Ni];
    int    lcwork = 2 * M.Ni + Sqr(M.Ni);
    Complex16 cwork[lcwork];
    int    i;
    T->Ni = M.Ni;
    T->Nj = M.Nj;
    T->QComplex = M.QComplex;
    if (InitMatrix(T)) return 1;
    strcpy(T->Name, "Eigenvectors of ");
    strncat(T->Name, M.Name, Max(0, MAXStrLength - 16));
    for (i = 0; i < M.Ni * M.Ni; i++)
    {
      MReIm[i].Re = M.VRe[0][i];
      MReIm[i].Im = M.VIm[0][i];
    }
    zheev(jobz, uplo, &M.Ni, MReIm, &M.Ni, val, cwork, &lcwork, work, &info);
    if (info != 0)
    {
      printf("zheev failed in Eigensystem with info value %i\n", info);
      printf("While working on matrix %s\n", M.Name);
      printf("continue diagonalization with zheevd\n");
      fflush(stdout);
      for (i = 0; i < M.Ni * M.Ni; i++)
      {
        MReIm[i].Re = M.VRe[0][i];
        MReIm[i].Im = M.VIm[0][i];
      }
      zheevd(jobz, uplo, &M.Ni, MReIm, &M.Ni, val, cwork, &lcwork, work, &lwork, iwork, &liwork, &info);
      if (info != 0)
      {
        printf("also zheevd failed in Eigensystem with info value %i\n", info);
        printf("ERORR!!!!!!!!!!\n");
        fflush(stdout);
        return 1;
      }
    }
    for (i = 0; i < M.Ni * M.Ni; i++)
    {
      T->VRe[0][i] = MReIm[i].Re;
      T->VIm[0][i] = -MReIm[i].Im;
    }
    return 0;
  }
  return 0;
}

uint32_t MatrixToEigensystem(MatrixType const M, MatrixType * const T, double *val)
{
  return Eigensystem(M, T, val);
}


/*  Calculates the eigen values (val) and eigen functions (T) of matrix M
 *  H=T^*.M.T^T with H a diagonal matrix with the values val on the diagonal
 *  input M defined, H and T not initiallized
 */
uint32_t EigensystemToCompact(MatrixType const M, CompactMatrixType * const T, double *val)
{
  int    info;
  const char *jobz = "V";
  const char *uplo = "U";
  int    lwork = 1 + 6 * M.Ni + 2 * Sqr(M.Ni);
  int    liwork = 3 + 5 * M.Ni;
  int    iwork[liwork];
  double work[lwork];
  if (M.Ni != M.Nj)
  {
    printf("Eigensystem called with non-square matrix\n");
    printf("ERORR!!!!!!\n");
    return 1;
  }
  if (M.QComplex == 0)
  {                             /* Real */
    if (MatrixCopyToCompact(T, M)) return 1;
    strcpy(T->Name, "Eigenvectors of ");
    strncat(T->Name, M.Name, Max(0, MAXStrLength - 16));
    dsyev(jobz, uplo, &M.Ni, T->VRe, &M.Ni, val, work, &lwork, &info);
    if (info)
    {
      printf("dsyv failed in Eigensystem with info value %i\n", info);
      printf("While working on matrix %s\n", M.Name);
      printf("continue diagonalization with dsyvd\n");
      fflush(stdout);
      FreeCompactMatrix(T);
      if (MatrixCopyToCompact(T, M)) return 1;
      strcpy(T->Name, "Eigenvectors of ");
      strncat(T->Name, M.Name, Max(0, MAXStrLength - 16));
      dsyevd(jobz, uplo, &M.Ni, T->VRe, &M.Ni, val, work, &lwork, iwork, &liwork, &info);
      if (info)
      {
        printf("also dsyvd failed in Eigensystem with info value %i\n", info);
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
    int    lcwork = 2 * M.Ni + Sqr(M.Ni);
    Complex16 cwork[lcwork];
    int    i;
    if (MatrixCopyToCompact(T, M)) return 1;
    zheev(jobz, uplo, &M.Ni, T->V, &M.Ni, val, cwork, &lcwork, work, &info);
    if (info != 0)
    {
      printf("zheev failed in Eigensystem with info value %i\n", info);
      printf("While working on matrix %s\n", M.Name);
      printf("continue diagonalization with zheevd\n");
      fflush(stdout);
      FreeCompactMatrix(T);
      if (MatrixCopyToCompact(T, M)) return 1;
      zheevd(jobz, uplo, &M.Ni, T->V, &M.Ni, val, cwork, &lcwork, work, &lwork, iwork, &liwork, &info);
      if (info != 0)
      {
        printf("also zheevd failed in Eigensystem with info value %i\n", info);
        printf("ERORR!!!!!!!!!!\n");
        fflush(stdout);
        return 1;
      }
    }
    for (i = 0; i < M.Ni * M.Ni; i++)
    {
      T->V[i].Im *= -1;
    }
    return 0;
  }
  return 0;
}

/*  Returns the Generalized Eigensystem of a non-orthnormal matrix.
 *  Given a basis psi_i
 *  Matrix M contains the values psi_i^* . H . psi_j
 *  Matrix S contains the values psi_i^* . psi_j
 *
 *  on exit val defines a diagonal matrix H and
 *  H = T . M . T^{\dag}
 *  I = T . S . T^{\dag}
 *  S = Sroot . Sroot^{\dag}
 *
 */
uint32_t GeneralizedEigensystem(MatrixType const M, MatrixType const S, MatrixType * const T, MatrixType * const Sroot,  double *val)
{
  //TODO: H = T . M . T^dag seems wrong should be H = T^* . M . T^T
  int    info, itype = 1;
  const char *jobz = "V";
  const char *uplo = "U";
  int    lwork = 1 + 6 * M.Ni + 2 * Sqr(M.Ni);
  int    liwork = 3 + 5 * M.Ni;
  int    iwork[liwork];
  double work[lwork];
  uint32_t i, j;
  if (M.Ni != M.Nj)
  {
    printf("Generalized Eigensystem called with non-square matrix\n");
    printf("ERORR!!!!!!\n");
    return 1;
  }
  if (M.QComplex == 0)
  {                             /* Real */
    if (MatrixCopy(Sroot, S))
    {
      printf("MatrixCopy failed to copy S to Sroot in GeneralizedEigensystem\n");
      return 1;
    }
    strcpy(Sroot->Name, "Sqrt of ");
    strncat(Sroot->Name, S.Name, Max(0, MAXStrLength - 8));
    for (i = 0; i < Sroot->Ni; i++)
    {
      for (j = i + 1; j < Sroot->Nj; j++)
      {
        Sroot->VRe[i][j] = 0;
      }
    }
    if (MatrixCopy(T, M))
    {
      printf("MatrixCopy failed to copy M to T in GeneralizedEigensystem\n");
      FreeMatrix(Sroot);
      return 1;
    }
    strcpy(T->Name, "Eigenvectors of ");
    strncat(T->Name, M.Name, Max(0, MAXStrLength - 16));
    dsygv(&itype, jobz, uplo, &M.Ni, T->VRe[0], &M.Ni, Sroot->VRe[0], &M.Ni, val, work, &lwork, &info);
    if (info)
    {
      printf("dsygv failed in Generalized Eigensystem with info value %i\n", info);
      printf("While working on matrix %s\n", M.Name);
      printf("continue diagonalization with dsygvd\n");
      fflush(stdout);
      FreeMatrix(T);
      FreeMatrix(Sroot);
      if (MatrixCopy(Sroot, S))
      {
        printf("MatrixCopy failed to copy S to Sroot in GeneralizedEigensystem\n");
        return 1;
      }
      strcpy(Sroot->Name, "Sqrt of ");
      strncat(Sroot->Name, S.Name, Max(0, MAXStrLength - 8));
      for (i = 0; i < Sroot->Ni; i++)
      {
        for (j = i + 1; j < Sroot->Nj; j++)
        {
          Sroot->VRe[i][j] = 0;
        }
      }
      if (MatrixCopy(T, M))
      {
        printf("MatrixCopy failed to copy M to T in GeneralizedEigensystem\n");
        FreeMatrix(Sroot);
        return 1;
      }
      strcpy(T->Name, "Eigenvectors of ");
      strncat(T->Name, M.Name, Max(0, MAXStrLength - 16));
      dsygvd(&itype, jobz, uplo, &M.Ni, T->VRe[0], &M.Ni, Sroot->VRe[0], &M.Ni, val, work, &lwork, iwork, &liwork,
             &info);
      if (info)
      {
        printf("also dsygvd failed in Eigensystem with info value %i\n", info);
        printf("ERORR!!!!!!!!!!\n");
        printf("ERORR!!!!!!!!!!\n");
        printf("ERORR!!!!!!!!!!\n");
        printf("ERORR!!!!!!!!!!\n");
        printf("ERORR!!!!!!!!!!\n");
        fflush(stdout);
        return 1;
      }
    }
    return 0;
  }
  else
  {                             /* complex */
    Complex16 MReIm[M.Ni * M.Ni];
    Complex16 SReIm[M.Ni * M.Ni];
    int    lcwork = 2 * M.Ni + Sqr(M.Ni);
    Complex16 cwork[lcwork];
    T->Ni = M.Ni;
    T->Nj = M.Nj;
    T->QComplex = M.QComplex;
    Sroot->Ni = M.Ni;
    Sroot->Nj = M.Nj;
    Sroot->QComplex = M.QComplex;
    if (InitMatrix(T))
      return 1;
    if (InitMatrix(Sroot))
      return 1;
    strcpy(T->Name, "Eigenvectors of ");
    strncat(T->Name, M.Name, Max(0, MAXStrLength - 16));
    strcpy(Sroot->Name, "Sqrt of ");
    strncat(Sroot->Name, S.Name, Max(0, MAXStrLength - 8));
    for (i = 0; i < M.Ni * M.Ni; i++)
    {
      MReIm[i].Re = M.VRe[0][i];
      MReIm[i].Im = M.VIm[0][i];
      SReIm[i].Re = S.VRe[0][i];
      SReIm[i].Im = S.VIm[0][i];
    }
    zhegv(&itype, jobz, uplo, &M.Ni, MReIm, &M.Ni, SReIm, &M.Ni, val, cwork, &lcwork, work, &info);
    if (info != 0)
    {
      printf("zhegv failed in Generalized Eigensystem with info value %i\n", info);
      printf("While working on matrix %s\n", M.Name);
      printf("continue diagonalization with zhegvd\n");
      fflush(stdout);
      for (i = 0; i < M.Ni * M.Ni; i++)
      {
        MReIm[i].Re = M.VRe[0][i];
        MReIm[i].Im = M.VIm[0][i];
        SReIm[i].Re = S.VRe[0][i];
        SReIm[i].Im = S.VIm[0][i];
      }
      zhegvd(&itype, jobz, uplo, &M.Ni, MReIm, &M.Ni, SReIm, &M.Ni, val, cwork, &lcwork, work, &lwork, iwork, &liwork,
             &info);
      if (info != 0)
      {
        printf("also zhegvd failed in Generalized Eigensystem with info value %i\n", info);
        printf("ERORR!!!!!!!!!!\n");
        printf("ERORR!!!!!!!!!!\n");
        printf("ERORR!!!!!!!!!!\n");
        printf("ERORR!!!!!!!!!!\n");
        printf("ERORR!!!!!!!!!!\n");
        printf("ERORR!!!!!!!!!!\n");
        fflush(stdout);
        return 1;
      }
    }
    for (i = 0; i < M.Ni * M.Ni; i++)
    {
      T->VRe[0][i] = MReIm[i].Re;
      T->VIm[0][i] = -MReIm[i].Im;
      Sroot->VRe[0][i] = SReIm[i].Re;
      Sroot->VIm[0][i] = -SReIm[i].Im;
    }
    for (i = 0; i < Sroot->Ni; i++)
    {
      for (j = i + 1; j < Sroot->Nj; j++)
      {
        Sroot->VRe[i][j] = 0;
        Sroot->VIm[i][j] = 0;
      }
    }
  }
  return 0;
}

uint32_t MatrixToGeneralizedEigensystem(MatrixType const M, MatrixType const S, MatrixType * const T, MatrixType * const Sroot, double *val)
{
  return GeneralizedEigensystem(M, S, T, Sroot, val);
}


/*  Returns the Generalized Eigenvalues of a non-orthnormal matrix.
 *  Given a basis psi_i
 *  Matrix M contains the values psi_i^* . H . psi_j
 *  Matrix S contains the values psi_i^* . psi_j
 *
 *  on exit val defines a diagonal matrix H and
 *  H = T . M . T^{\dag}
 *  I = T . S . T^{\dag}
 *  S = Sroot . Sroot^{\dag}
 *
 */
uint32_t GeneralizedEigenvalues(MatrixType const M, MatrixType const S, double *val)
{
  //TODO: H = T . M . T^dag seems wrong should be H = T^* . M . T^T
  int    info, itype = 1;
  const char *jobz = "N";
  const char *uplo = "U";
  int    lwork = 1 + 6 * M.Ni + 2 * Sqr(M.Ni);
  int    liwork = 3 + 5 * M.Ni;
  int    iwork[liwork];
  double work[lwork];
  uint32_t i, j;
  MatrixType Sroot, T;
  if (M.Ni != M.Nj)
  {
    printf("Generalized Eigenvalues called with non-square matrix\n");
    printf("ERORR!!!!!!\n");
    return 1;
  }
  if (M.QComplex == 0)
  {                             /* Real */
    if (MatrixCopy(&Sroot, S))
    {
      printf("MatrixCopy failed to copy S to Sroot in GeneralizedEigenvalues\n");
      return 1;
    }
    strcpy(Sroot.Name, "Sqrt of ");
    strncat(Sroot.Name, S.Name, Max(0, MAXStrLength - 8));
    for (i = 0; i < Sroot.Ni; i++)
    {
      for (j = i + 1; j < Sroot.Nj; j++)
      {
        Sroot.VRe[i][j] = 0;
      }
    }
    if (MatrixCopy(&T, M))
    {
      printf("MatrixCopy failed to copy M to T in GeneralizedEigenvalues\n");
      FreeMatrix(&Sroot);
      return 1;
    }
    strcpy(T.Name, "Eigenvectors of ");
    strncat(T.Name, M.Name, Max(0, MAXStrLength - 16));
    dsygv(&itype, jobz, uplo, &M.Ni, T.VRe[0], &M.Ni, Sroot.VRe[0], &M.Ni, val, work, &lwork, &info);
    if (info)
    {
      printf("dsygv failed in Generalized Eigenvalues with info value %i\n", info);
      printf("While working on matrix %s\n", M.Name);
      printf("continue diagonalization with dsygvd\n");
      fflush(stdout);
      FreeMatrix(&T);
      FreeMatrix(&Sroot);
      if (MatrixCopy(&Sroot, S))
      {
        printf("MatrixCopy failed to copy S to Sroot in GeneralizedEigenvalues\n");
        return 1;
      }
      strcpy(Sroot.Name, "Sqrt of ");
      strncat(Sroot.Name, S.Name, Max(0, MAXStrLength - 8));
      for (i = 0; i < Sroot.Ni; i++)
      {
        for (j = i + 1; j < Sroot.Nj; j++)
        {
          Sroot.VRe[i][j] = 0;
        }
      }
      if (MatrixCopy(&T, M))
      {
        printf("MatrixCopy failed to copy M to T in GeneralizedEigenvalues\n");
        FreeMatrix(&Sroot);
        return 1;
      }
      strcpy(T.Name, "Eigenvectors of ");
      strncat(T.Name, M.Name, Max(0, MAXStrLength - 16));
      dsygvd(&itype, jobz, uplo, &M.Ni, T.VRe[0], &M.Ni, Sroot.VRe[0], &M.Ni, val, work, &lwork, iwork, &liwork, &info);
      if (info)
      {
        printf("also dsygvd failed in Eigenvalues with info value %i\n", info);
        printf("ERORR!!!!!!!!!!\n");
        printf("ERORR!!!!!!!!!!\n");
        printf("ERORR!!!!!!!!!!\n");
        printf("ERORR!!!!!!!!!!\n");
        printf("ERORR!!!!!!!!!!\n");
        fflush(stdout);
        return 1;
      }
    }
    FreeMatrix(&Sroot);
    FreeMatrix(&T);
    return 0;
  }
  else
  {                             /* complex */
    Complex16 MReIm[M.Ni * M.Ni];
    Complex16 SReIm[M.Ni * M.Ni];
    int    lcwork = 2 * M.Ni + Sqr(M.Ni);
    Complex16 cwork[lcwork];
    T.Ni = M.Ni;
    T.Nj = M.Nj;
    T.QComplex = M.QComplex;
    Sroot.Ni = M.Ni;
    Sroot.Nj = M.Nj;
    Sroot.QComplex = M.QComplex;
    if (InitMatrix(&T))
      return 1;
    if (InitMatrix(&Sroot))
      return 1;
    strcpy(T.Name, "Eigenvectors of ");
    strncat(T.Name, M.Name, Max(0, MAXStrLength - 16));
    strcpy(Sroot.Name, "Sqrt of ");
    strncat(Sroot.Name, S.Name, Max(0, MAXStrLength - 8));
    for (i = 0; i < M.Ni * M.Ni; i++)
    {
      MReIm[i].Re = M.VRe[0][i];
      MReIm[i].Im = M.VIm[0][i];
      SReIm[i].Re = S.VRe[0][i];
      SReIm[i].Im = S.VIm[0][i];
    }
    zhegv(&itype, jobz, uplo, &M.Ni, MReIm, &M.Ni, SReIm, &M.Ni, val, cwork, &lcwork, work, &info);
    if (info != 0)
    {
      printf("zhegv failed in Generalized Eigenvalues with info value %i\n", info);
      printf("While working on matrix %s\n", M.Name);
      printf("continue diagonalization with zhegvd\n");
      fflush(stdout);
      for (i = 0; i < M.Ni * M.Ni; i++)
      {
        MReIm[i].Re = M.VRe[0][i];
        MReIm[i].Im = M.VIm[0][i];
        SReIm[i].Re = S.VRe[0][i];
        SReIm[i].Im = S.VIm[0][i];
      }
      zhegvd(&itype, jobz, uplo, &M.Ni, MReIm, &M.Ni, SReIm, &M.Ni, val, cwork, &lcwork, work, &lwork, iwork, &liwork,
             &info);
      if (info != 0)
      {
        printf("also zhegvd failed in Generalized Eigenvalues with info value %i\n", info);
        printf("ERORR!!!!!!!!!!\n");
        printf("ERORR!!!!!!!!!!\n");
        printf("ERORR!!!!!!!!!!\n");
        printf("ERORR!!!!!!!!!!\n");
        printf("ERORR!!!!!!!!!!\n");
        printf("ERORR!!!!!!!!!!\n");
        fflush(stdout);
        return 1;
      }
    }
    for (i = 0; i < M.Ni * M.Ni; i++)
    {
      T.VRe[0][i] = MReIm[i].Re;
      T.VIm[0][i] = -MReIm[i].Im;
      Sroot.VRe[0][i] = SReIm[i].Re;
      Sroot.VIm[0][i] = -SReIm[i].Im;
    }
    for (i = 0; i < Sroot.Ni; i++)
    {
      for (j = i + 1; j < Sroot.Nj; j++)
      {
        Sroot.VRe[i][j] = 0;
        Sroot.VIm[i][j] = 0;
      }
    }
  }
  FreeMatrix(&Sroot);
  FreeMatrix(&T);
  return 0;
}


/*  Returns the Generalized Eigensystem of a non-orthnormal matrix.
 *  Given a basis psi_i
 *  Matrix M contains the values psi_i^* . H . psi_j
 *  Matrix S contains the values psi_i^* . psi_j
 *
 *  on exit val defines a diagonal matrix H and
 *  H = T . M . T^{\dag}
 *  I = T . S . T^{\dag}
 *  S = Sroot . Sroot^{\dag}
 *
 */
uint32_t GeneralizedEigensystemSingular(MatrixType const M, MatrixType const S, MatrixType * const T, MatrixType * const Sroot, double *val)
{
  //TODO: test H = T . M . T^{\dag} does not seem right, should be H = T^* . M . T^T
  uint32_t i, j, NewDim, retval = 0;
  double Sval[S.Ni], bound;
  MatrixType EigenSysS, TMT, TST, Ttemp;
  MatrixCopy(&TST, S);
  MatrixScalarMultiply(&TST, -1, 0);
  Eigensystem(TST, &EigenSysS, Sval);
  MatrixScalarMultiply(&TST, 0, 0);
  MatrixTransform(&TMT, M, EigenSysS);
  NewDim = 1;
  bound = 0;
  for (i = 0; i < M.Ni; i++)
  {
    bound -= Sval[i];
  }
  bound *= DBL_EPSILON * 100;
  bound += DBL_MIN * 100;
  for (i = 0; i < M.Ni; i++)
  {
    if (-Sval[i] > bound)
    {
      TST.VRe[i][i] = -Sval[i];
      NewDim = i + 1;
    }
    else
    {
      for (j = 0; j < M.Nj; j++)
      {
        TMT.VRe[i][j] = 0;
        TMT.VRe[j][i] = 0;
      }
      TST.VRe[i][i] = 1;
    }
  }
  printf("Singularity border =%10.3E : Number of states kept = %i\n", bound, NewDim);
  printf("Norm of basis ( 0...7): ");
  for (i = 0; i < Min(8, S.Ni); i++)
  {
    printf("%10.3E ", -Sval[i]);
  }
  printf("..........\n");
  if (GeneralizedEigenvalues(TMT, TST, val))
  {
    printf("GeneralizedEigenvalues failed in GeneralizedEigensystemSingular\n");
    retval = 1;
  }
  for (i = NewDim; i < M.Ni; i++)
  {
    TMT.VRe[i][i] = val[M.Ni - 1] + 100 * DBL_EPSILON * Abs(val[M.Ni - 1]);
  }
  if (GeneralizedEigensystem(TMT, TST, &Ttemp, Sroot, val))
  {
    printf("GeneralizedEigensystem failed in GeneralizedEigensystemSingular\n");
    retval = 1;
  }
  for (i = NewDim; i < M.Ni; i++)
  {
    for (j = 0; j < M.Ni; j++)
    {
      Ttemp.VRe[j][i] = 0;
    }
  }
  MatrixMultiply(T, Ttemp, EigenSysS);
  FreeMatrix(&EigenSysS);
  FreeMatrix(&TMT);
  FreeMatrix(&TST);
  FreeMatrix(&Ttemp);
  return retval;
}

/*
 *
 */
uint32_t RemoveStatesWithNegligibleNorm(MatrixType * const H, MatrixType * const S, MatrixType * const T)
{
  MatrixType temp;
  double val[S->Ni];
  MatrixScalarMultiply(S, -1, 0);
  Eigensystem(*S, T, val);
  MatrixScalarMultiply(S, -1, 0);
  MatrixCopy(&temp, *H);
  MatrixTransform(H, temp, *T);
  FreeMatrix(&temp);
  MatrixCopy(&temp, *S);
  MatrixTransform(S, temp, *T);
  FreeMatrix(&temp);
  return 0;
}


/*
*  M=T.H.T^{\dag}
*  input M defined, H and T not initiallized
*/
uint32_t MatrixToAndersonMatrix(MatrixType const M, MatrixType * const T, AndersonMatrixType * const H)
{
  //TODO: M=T.H.T^dag does not seem right, should be H=T^*.M.T^T
  int    i, j, info;
  const char *jobz = "V";
  const char *uplo = "U";
  int    lwork = 1 + 6 * M.Ni + 2 * Sqr(M.Ni);
  int    liwork = 3 + 5 * M.Ni;
  int    iwork[liwork];
  double work[lwork];
  MatrixType U;
  if (M.Ni != M.Nj)
  {
    printf("MatrixToAndersonMatrix called with non-square matrix\n");
    printf("ERORR!!!!!!\n");
    return 1;
  }
  H->N = M.Ni;
  if (InitAndersonMatrix(H))
    return 1;
  H->a[0] = 0;
  H->b[0] = 1;
  if (M.QComplex == 0)
  {                             /* Real */
    if (MatrixCopy(T, M))
      return 1;
    for (i = 1; i < T->Ni; i++)
    {
      T->VRe[0][i] = 0;
      T->VRe[i][0] = 0;
    }
    H->a[1] = M.VRe[0][0];
    T->VRe[0][0] = 1;
    strcpy(T->Name, "Eigenvectors of Minor of ");
    strncat(T->Name, M.Name, Max(0, MAXStrLength - 30));
    i = T->Ni - 1;
    dsyev(jobz, uplo, &i, &(T->VRe[1][1]), &M.Ni, &(H->a[2]), work, &lwork, &info);
    if (info)
    {
      printf("dseyv failed in MatrixToAndersonMatrix with info value %i\n", info);
      printf("While working on matrix %s\n", M.Name);
      printf("continue diagonalization with dsyvd\n");
      fflush(stdout);
      dsyevd(jobz, uplo, &i, &(T->VRe[1][1]), &M.Ni, &(H->a[2]), work, &lwork, iwork, &liwork, &info);
      if (info)
      {
        printf("also dsyevd failed in MatrixToAndersonMatrix with info value %i\n", info);
        printf("ERORR!!!!!!!!!!\n");
        fflush(stdout);
        return 1;
      }
    }
    if (MatrixTransform(&U, M, *T))
      return 1;
    for (i = 1; i < T->Ni; i++)
    {
      if (U.VRe[0][i] < 0)
      {
        H->b[i] = -U.VRe[0][i];
        for (j = 0; j < T->Ni; j++)
        {
          T->VRe[i][j] *= -1;
        }
      }
      else
      {
        H->b[i] = U.VRe[0][i];
      }
    }
    return 0;
  }
  else
  {                             /* Complex */
    Complex16 MReIm[M.Ni * M.Ni];
    Complex16 phase;
    int    lcwork = 2 * M.Ni + Sqr(M.Ni);
    Complex16 cwork[lcwork];
    T->Ni = M.Ni;
    T->Nj = M.Nj;
    T->QComplex = M.QComplex;
    strcpy(T->Name, "Eigenvectors of Minor of ");
    strncat(T->Name, M.Name, Max(0, MAXStrLength - 30));
    if (InitMatrix(T))
      return 1;
    for (i = 0; i < M.Ni * M.Ni; i++)
    {
      MReIm[i].Re = M.VRe[0][i];
      MReIm[i].Im = M.VIm[0][i];
    }
    for (i = 1; i < T->Ni; i++)
    {
      MReIm[i].Re = 0;
      MReIm[i].Im = 0;
      MReIm[i * M.Ni].Re = 0;
      MReIm[i * M.Ni].Im = 0;
    }
    MReIm[0].Re = 1;
    MReIm[0].Im = 0;
    i = T->Ni - 1;
    zheev(jobz, uplo, &i, &MReIm[1 + M.Ni], &M.Ni, &(H->a[2]), cwork, &lcwork, work, &info);
    if (info != 0)
    {
      printf("zheev failed in Eigensystem with info value %i\n", info);
      printf("While working on matrix %s\n", M.Name);
      printf("continue diagonalization with zheevd\n");
      fflush(stdout);
      zheevd(jobz, uplo, &i, &MReIm[1 + M.Ni], &M.Ni, &(H->a[2]), cwork, &lcwork, work, &lwork, iwork, &liwork, &info);
      if (info != 0)
      {
        printf("also zheevd failed in Eigensystem with info value %i\n", info);
        printf("ERORR!!!!!!!!!!\n");
        fflush(stdout);
        return 1;
      }
    }
    H->a[1] = M.VRe[0][0];
    for (i = 0; i < M.Ni * M.Ni; i++)
    {
      T->VRe[0][i] = MReIm[i].Re;
      T->VIm[0][i] = -MReIm[i].Im;
    }
    if (MatrixTransform(&U, M, *T))
      return 1;
    for (i = 1; i < T->Ni; i++)
    {
      if (fabs(U.VRe[0][i]) > fabs(U.VIm[0][i]))
      {
        phase.Im = 1;
        phase.Re = -U.VRe[0][i] / U.VIm[0][i];
      }
      else
      {
        phase.Re = 1;
        phase.Im = -U.VIm[0][i] / U.VRe[0][i];
      }
      work[0] = sqrt(Sqr(phase.Re) + Sqr(phase.Im));
      phase.Im /= work[0];
      phase.Re /= work[0];
      if ((U.VRe[0][i] * phase.Re - U.VIm[0][i] * phase.Im) < 0)
      {
        phase.Re *= -1;
        phase.Im *= -1;
      }
      for (j = 0; j < T->Ni; j++)
      {
        work[0] = T->VRe[i][j] * phase.Re - T->VIm[i][j] * phase.Im;
        T->VIm[i][j] = T->VRe[i][j] * phase.Im + T->VIm[i][j] * phase.Re;
        T->VRe[i][j] = work[0];
      }
      H->b[i] = U.VRe[0][i] * phase.Re - U.VIm[0][i] * phase.Im;
    }
    return 0;
  }
  return 0;
}

/*
*  M=T.H.T^{\dag}
*  input M defined, H and T not initiallized
******************* still to do ***************************
*  a[i>1] is diagonal
*  b[i] is random
*  one could make b[i] such that sum_j b[i][k][j] is real and positive
*/
uint32_t BlokTriDiagonalize(MatrixType const M, MatrixType * const T, BlokTriDiagonalMatrixType * const H, const uint32_t BlokDim)
{
  //TODO: M=T.H.T^Dag does not seem right, should be H=T^*.M.T^T
  int    i, j, k, l, CurrentDim, NBlok;
  const char *jobu = "N";
  const char *jobvt = "A";
  const char *jobz = "V";
  const char *uplo = "U";
  double SingularVal[M.Ni], val[M.Ni];
  double *U;
  const int lwork = 1 + 6 * M.Ni + 2 * Sqr(M.Ni);
  double work[lwork];
  const int liwork = 3 + 5 * M.Ni;
  int    iwork[liwork];
  double diagval[BlokDim];
  int    info;
  MatrixType THT;
  Complex16 *UC, *THTC;
  if (M.Ni != M.Nj)
  {
    printf("BlokTriDiagonalize called with non-square matrix\n");
    printf("ERORR!!!!!!\n");
    return 1;
  }
  T->Ni = M.Ni;
  T->Nj = M.Nj;
  T->QComplex = M.QComplex;
  if (InitIdentityMatrix(T))
    return 1;
  strcpy(T->Name, "BlokTridiagonalizer of ");
  strncat(T->Name, M.Name, Max(0, MAXStrLength - 28));
  NBlok = ceil(((double)M.Ni) / ((double)BlokDim));
  H->N = NBlok;
  if (InitFrameBlokTriDiagonalMatrix(H))
  {
    FreeMatrix(T);
    return 1;
  }

  if (M.QComplex == 0)
  {                             /* Real */
    CurrentDim = BlokDim;
    for (i = BlokDim; i < M.Ni - CurrentDim; i += CurrentDim)
    {
      /* calculate the values on the diagonal of the previous created blok */
      if (MatrixTransformPart(&THT, M, *T, i - CurrentDim, i - 1, i - CurrentDim, i - 1))
      {
        FreeMatrix(T);
        FreeBlokTriDiagonalMatrix(H);
        return 1;
      };
      for (j = 0; j < CurrentDim; j++)
      {
        diagval[j] = THT.VRe[j][j];
      }
      FreeMatrix(&THT);
      /* make another CurrentDim row's zero */
      if (MatrixTransformPart(&THT, M, *T, i, M.Ni - 1, i - CurrentDim, i - 1))
      {
        FreeMatrix(T);
        FreeBlokTriDiagonalMatrix(H);
        return 1;
      };
      U = malloc(THT.Ni * THT.Ni * sizeof(double));
      if (U == NULL)
      {
        FreeMatrix(T);
        FreeMatrix(&THT);
        FreeBlokTriDiagonalMatrix(H);
        return 1;
      }
      dgesvd(jobu, jobvt, &(THT.Nj), &(THT.Ni), THT.VRe[0], &(THT.Nj), SingularVal, U, &(THT.Ni), U, &(THT.Ni), work,&lwork, &info);
      FreeMatrix(&THT);
      if (info)
      {
        printf("Possible error in dgesvd in BlokTriDiagonalize\n");
        printf("info is %i", info);
      }
      for (j = 0; j < M.Nj; j++)
      {
        for (k = i; k < M.Ni; k++)
        {
          work[k] = 0;
          for (l = i; l < M.Ni; l++)
          {
            work[k] += U[(l - i) * (M.Ni - i) + k - i] * T->VRe[l][j];
          }
        }
        for (l = i; l < M.Ni; l++)
        {
          T->VRe[l][j] = work[l];
        }
      }
      free(U);
      while ((SingularVal[CurrentDim - 1] < fabs(diagval[CurrentDim - 1]) * DBL_EPSILON * 100) && (CurrentDim > 1))
      {
        CurrentDim--;
      }
      if ((SingularVal[CurrentDim - 1] < fabs(diagval[CurrentDim - 1]) * DBL_EPSILON * 100))
      {                         /* Diagonalize the right bottom block of the matrix and leave */

        i = M.Ni + 1;
      }
      else
      {                         /* Diagonalize the right bottom block of the matrix, the new block on the diagonal and continue */
        if (MatrixTransformPart(&THT, M, *T, i + CurrentDim, M.Ni - 1, i + CurrentDim, M.Ni - 1))
        {
          FreeMatrix(T);
          FreeBlokTriDiagonalMatrix(H);
          return 1;
        };
        dsyev(jobz, uplo, &THT.Ni, THT.VRe[0], &THT.Ni, val, work, &lwork, &info);
        if (info)
        {
          printf("dsyv failed in BlokTriDiagonalize with info value %i\n", info);
          printf("While working on matrix %s\n", M.Name);
          printf("continue diagonalization with dsyvd\n");
          fflush(stdout);
          dsyevd(jobz, uplo, &THT.Ni, THT.VRe[0], &THT.Ni, val, work, &lwork, iwork, &liwork, &info);
          if (info)
          {
            printf("also dsyvd failed in BlokTriDiagonalize with info value %i\n", info);
            printf("ERORR!!!!!!!!!!\n");
            fflush(stdout);
            return 1;
          }
        }
        for (j = 0; j < M.Nj; j++)
        {
          for (k = i + CurrentDim; k < M.Ni; k++)
          {
            work[k] = 0;
            for (l = i + CurrentDim; l < M.Ni; l++)
            {
              work[k] += THT.VRe[(k - i - CurrentDim)][l - i - CurrentDim] * T->VRe[l][j];
            }
          }
          for (l = i + CurrentDim; l < M.Ni; l++)
          {
            T->VRe[l][j] = work[l];
          }
        }
        FreeMatrix(&THT);
        /* Diagonalize the newly formed blok */
        if (MatrixTransformPart
            (&THT, M, *T, i, Min(M.Ni - 1, i + CurrentDim - 1), i, Min(M.Ni - 1, i + CurrentDim - 1)))
        {
          FreeMatrix(T);
          FreeBlokTriDiagonalMatrix(H);
          return 1;
        };
        dsyev(jobz, uplo, &THT.Ni, THT.VRe[0], &THT.Ni, val, work, &lwork, &info);
        if (info)
        {
          printf("dsyv failed in BlokTriDiagonalize with info value %i\n", info);
          printf("While working on matrix %s\n", M.Name);
          printf("continue diagonalization with dsyvd\n");
          fflush(stdout);
          dsyevd(jobz, uplo, &THT.Ni, THT.VRe[0], &THT.Ni, val, work, &lwork, iwork, &liwork, &info);
          if (info)
          {
            printf("also dsyvd failed in BlokTriDiagonalize with info value %i\n", info);
            printf("ERORR!!!!!!!!!!\n");
            fflush(stdout);
            return 1;
          }
        }
        for (j = 0; j < M.Nj; j++)
        {
          for (k = i; k < Min(M.Ni, i + CurrentDim); k++)
          {
            work[k] = 0;
            for (l = i; l < Min(M.Ni, i + CurrentDim); l++)
            {
              work[k] += THT.VRe[(k - i)][l - i] * T->VRe[l][j];
            }
          }
          for (l = i; l < Min(M.Ni, i + CurrentDim); l++)
          {
            T->VRe[l][j] = work[l];
          }
        }
        FreeMatrix(&THT);
      }
      /* additional orhtonormalization of T */
      MatrixOrthonormalize(T);
    }                           /* end loop over all bloks (i) */
    /* create the blok tri-diagonal matrix H */
    for (i = 0; i < H->N; i++)
    {
      if (MatrixTransformPart
          (&(H->a[i]), M, *T, i * BlokDim, Min(M.Ni, (i + 1) * (BlokDim) - 1), i * BlokDim,
           Min(M.Nj, (i + 1) * (BlokDim) - 1)))
      {
        printf("Unexpected out of memory in BlokTriDiagonalize\n");
        return 1;
      }
    }
    for (i = 0; i < H->N - 1; i++)
    {
      if (MatrixTransformPart
          (&(H->b[i]), M, *T, i * BlokDim, Min(M.Ni, (i + 1) * (BlokDim) - 1), (i + 1) * BlokDim,
           Min(M.Nj, (i + 2) * (BlokDim) - 1)))
      {
        printf("Unexpected out of memory in BlokTriDiagonalize\n");
        return 1;
      }
    }
  }                             /* end Real */
  else
  {                             /* Complex */
    const int Clwork = 1 + 2 * M.Ni + Sqr(M.Ni);
    Complex16 Cwork[Clwork];
    CurrentDim = BlokDim;
    for (i = BlokDim; i < M.Ni - CurrentDim; i += CurrentDim)
    {
      /* calculate the values on the diagonal of the previous created blok */
      if (MatrixTransformPart(&THT, M, *T, i - CurrentDim, i - 1, i - CurrentDim, i - 1))
      {
        printf("MatrixTransformPart failed in BlokTriDiagonalize C 10\n");
        fflush(stdout);
        FreeMatrix(T);
        FreeBlokTriDiagonalMatrix(H);
        return 1;
      };
      for (j = 0; j < CurrentDim; j++)
      {
        diagval[j] = THT.VRe[j][j];
      }
      FreeMatrix(&THT);
      /* make another CurrentDim row's zero */
      if (MatrixTransformPart(&THT, M, *T, i, M.Ni - 1, i - CurrentDim, i - 1))
      {
        printf("MatrixTransformPart failed in BlokTriDiagonalize C 11\n");
        fflush(stdout);
        FreeMatrix(T);
        FreeBlokTriDiagonalMatrix(H);
        return 1;
      };
      UC = malloc(THT.Ni * THT.Ni * sizeof(Complex16));
      if (UC == NULL)
      {
        printf("Allocation of UC failed in BlokTriDiagonalize C 20\n");
        fflush(stdout);
        FreeMatrix(T);
        FreeMatrix(&THT);
        FreeBlokTriDiagonalMatrix(H);
        return 1;
      }
      THTC = malloc(THT.Ni * THT.Nj * sizeof(Complex16));
      if (THTC == NULL)
      {
        printf("Allocation of THTC failed in BlokTriDiagonalize C 21\n");
        fflush(stdout);
        FreeMatrix(T);
        FreeMatrix(&THT);
        FreeBlokTriDiagonalMatrix(H);
        free(UC);
        return 1;
      }
      for (j = 0; j < THT.Ni * THT.Nj; j++)
      {
        THTC[j].Re = THT.VRe[0][j];
        THTC[j].Im = THT.VIm[0][j];
      }
      zgesvd(jobu, jobvt, &(THT.Nj), &(THT.Ni), THTC, &(THT.Nj), SingularVal, UC, &(THT.Ni), UC, &(THT.Ni), Cwork,
             &Clwork, work, &info);
      FreeMatrix(&THT);
      free(THTC);
      if (info)
      {
        printf("Possible error in dgesvd in BlokTriDiagonalize\n");
        printf("info is %i", info);
        fflush(stdout);
      }
      for (j = 0; j < M.Nj; j++)
      {
        for (k = i; k < M.Ni; k++)
        {
          Cwork[k].Re = 0;
          Cwork[k].Im = 0;
          for (l = i; l < M.Ni; l++)
          {
            Cwork[k].Re +=
              UC[(l - i) * (M.Ni - i) + k - i].Re * T->VRe[l][j] - UC[(l - i) * (M.Ni - i) + k - i].Im * T->VIm[l][j];
            Cwork[k].Im +=
              UC[(l - i) * (M.Ni - i) + k - i].Re * T->VIm[l][j] + UC[(l - i) * (M.Ni - i) + k - i].Im * T->VRe[l][j];
          }
        }
        for (l = i; l < M.Ni; l++)
        {
          T->VRe[l][j] = Cwork[l].Re;
          T->VIm[l][j] = Cwork[l].Im;
        }
      }
      free(UC);
      while ((SingularVal[CurrentDim - 1] < fabs(diagval[CurrentDim - 1]) * DBL_EPSILON * 100) && (CurrentDim > 1))
      {
        CurrentDim--;
      }
      if ((SingularVal[CurrentDim - 1] < fabs(diagval[CurrentDim - 1]) * DBL_EPSILON * 100))
      {                         /* Diagonalize the right bottom block of the matrix and leave */

        i = M.Ni + 1;
      }
      else
      {                         /* Diagonalize the right bottom block of the matrix, the new block on the diagonal and continue */
        if (MatrixTransformPart(&THT, M, *T, i + CurrentDim, M.Ni - 1, i + CurrentDim, M.Ni - 1))
        {
          FreeMatrix(T);
          FreeBlokTriDiagonalMatrix(H);
          return 1;
        };
        THTC = malloc(THT.Ni * THT.Nj * sizeof(Complex16));
        if (THTC == NULL)
        {
          printf("Allocation of THTC failed in BlokTriDiagonalize C 40\n");
          fflush(stdout);
          FreeMatrix(T);
          FreeMatrix(&THT);
          FreeBlokTriDiagonalMatrix(H);
          return 1;
        }
        for (j = 0; j < THT.Ni * THT.Nj; j++)
        {
          THTC[j].Re = THT.VRe[0][j];
          THTC[j].Im = THT.VIm[0][j];
        }
        zheev(jobz, uplo, &THT.Ni, THTC, &THT.Ni, val, Cwork, &Clwork, work, &info);
        if (info)
        {
          printf("zheev failed in BlokTriDiagonalize with info value %i\n", info);
          printf("While working on matrix %s\n", M.Name);
          printf("continue diagonalization with zheevd\n");
          fflush(stdout);
          zheevd(jobz, uplo, &THT.Ni, THTC, &THT.Ni, val, Cwork, &Clwork, work, &lwork, iwork, &liwork, &info);
          if (info)
          {
            printf("also zheev failed in BlokTriDiagonalize with info value %i\n", info);
            printf("ERORR!!!!!!!!!!\n");
            fflush(stdout);
            return 1;
          }
        }
        for (j = 0; j < M.Nj; j++)
        {
          for (k = i + CurrentDim; k < M.Ni; k++)
          {
            Cwork[k].Re = 0;
            Cwork[k].Im = 0;
            for (l = i + CurrentDim; l < M.Ni; l++)
            {
              Cwork[k].Re +=
                THTC[(k - i - CurrentDim) * THT.Ni + l - i - CurrentDim].Re * T->VRe[l][j] +
                THTC[(k - i - CurrentDim) * THT.Ni + l - i - CurrentDim].Im * T->VIm[l][j];
              Cwork[k].Im +=
                THTC[(k - i - CurrentDim) * THT.Ni + l - i - CurrentDim].Re * T->VIm[l][j] -
                THTC[(k - i - CurrentDim) * THT.Ni + l - i - CurrentDim].Im * T->VRe[l][j];
            }
          }
          for (l = i + CurrentDim; l < M.Ni; l++)
          {
            T->VRe[l][j] = Cwork[l].Re;
            T->VIm[l][j] = Cwork[l].Im;
          }
        }
        FreeMatrix(&THT);
        free(THTC);
        /* Diagonalize the newly formed blok */
        if (MatrixTransformPart
            (&THT, M, *T, i, Min(M.Ni - 1, i + CurrentDim - 1), i, Min(M.Ni - 1, i + CurrentDim - 1)))
        {
          FreeMatrix(T);
          FreeBlokTriDiagonalMatrix(H);
          return 1;
        };
        THTC = malloc(THT.Ni * THT.Nj * sizeof(Complex16));
        if (THTC == NULL)
        {
          printf("Allocation of THTC failed in BlokTriDiagonalize C 40\n");
          fflush(stdout);
          FreeMatrix(T);
          FreeMatrix(&THT);
          FreeBlokTriDiagonalMatrix(H);
          return 1;
        }
        for (j = 0; j < THT.Ni * THT.Nj; j++)
        {
          THTC[j].Re = THT.VRe[0][j];
          THTC[j].Im = THT.VIm[0][j];
        }
        zheev(jobz, uplo, &THT.Ni, THTC, &THT.Ni, val, Cwork, &Clwork, work, &info);
        if (info)
        {
          printf("zheev failed in BlokTriDiagonalize with info value %i\n", info);
          printf("While working on matrix %s\n", M.Name);
          printf("continue diagonalization with zheevd\n");
          fflush(stdout);
          zheevd(jobz, uplo, &THT.Ni, THTC, &THT.Ni, val, Cwork, &Clwork, work, &lwork, iwork, &liwork, &info);
          if (info)
          {
            printf("also zheevd failed in BlokTriDiagonalize with info value %i\n", info);
            printf("ERORR!!!!!!!!!!\n");
            fflush(stdout);
            return 1;
          }
        }
        for (j = 0; j < M.Nj; j++)
        {
          for (k = i; k < Min(M.Ni, i + CurrentDim); k++)
          {
            Cwork[k].Re = 0;
            Cwork[k].Im = 0;
            for (l = i; l < Min(M.Ni, i + CurrentDim); l++)
            {
              Cwork[k].Re +=
                THTC[(k - i) * THT.Ni + l - i].Re * T->VRe[l][j] + THTC[(k - i) * THT.Ni + l - i].Im * T->VIm[l][j];
              Cwork[k].Im +=
                THTC[(k - i) * THT.Ni + l - i].Re * T->VIm[l][j] - THTC[(k - i) * THT.Ni + l - i].Im * T->VRe[l][j];
            }
          }
          for (l = i; l < Min(M.Ni, i + CurrentDim); l++)
          {
            T->VRe[l][j] = Cwork[l].Re;
            T->VIm[l][j] = Cwork[l].Im;
          }
        }
        FreeMatrix(&THT);
        free(THTC);
      }
      /* additional orhtonormalization of T */
      MatrixOrthonormalize(T);
    }                           /* end loop over all bloks (i) */
    /* create the blok tri-diagonal matrix H */
    for (i = 0; i < H->N; i++)
    {
      if (MatrixTransformPart
          (&(H->a[i]), M, *T, i * BlokDim, Min(M.Ni, (i + 1) * (BlokDim) - 1), i * BlokDim,
           Min(M.Nj, (i + 1) * (BlokDim) - 1)))
      {
        printf("Unexpected out of memory in BlokTriDiagonalize 30\n");
        return 1;
      }
    }
    for (i = 0; i < H->N - 1; i++)
    {
      if (MatrixTransformPart
          (&(H->b[i]), M, *T, i * BlokDim, Min(M.Ni, (i + 1) * (BlokDim) - 1), (i + 1) * BlokDim,
           Min(M.Nj, (i + 2) * (BlokDim) - 1)))
      {
        printf("Unexpected out of memory in BlokTriDiagonalize 31\n");
        return 1;
      }
    }
  }                             /* end complex */
  return 0;
}

uint32_t MatrixToBlokTriDiagonalMatrix(MatrixType const M, MatrixType * const T, BlokTriDiagonalMatrixType * const H, const uint32_t NBlok)
{
  return BlokTriDiagonalize(M, T, H, NBlok);
}

/*
*  M=T.H.T^{\dag}
*  input M defined, H and T not initiallized
*/
uint32_t MatrixToBlokAndersonMatrix(MatrixType const M, MatrixType * const T, BlokAndersonMatrixType * const H, const uint32_t BlokDim)
{
  //TODO: M=T.H.T^Dag does not seem right, should be H=T^*.M.T^T
  int    i, j, info, NBlok;
  const char *jobz = "V";
  const char *uplo = "U";
  int    lwork = 1 + 6 * M.Ni + 2 * Sqr(M.Ni);
  int    liwork = 3 + 5 * M.Ni;
  int    iwork[liwork];
  double work[lwork];
  double val[M.Ni];
  MatrixType U;
  Complex16 phase;
  if (M.Ni != M.Nj)
  {
    printf("MatrixToAndersonMatrix called with non-square matrix\n");
    printf("ERORR!!!!!!\n");
    return 1;
  }
  NBlok = ceil(((double)M.Ni) / ((double)BlokDim));
  H->N = NBlok;
  if (InitFrameBlokAndersonMatrix(H))
    return 1;
  if (M.QComplex == 0)
  {                             /* Real */
    if (MatrixCopy(T, M))
    {
      FreeFrameBlokAndersonMatrix(H);
      return 1;
    }
    for (i = 0; i < T->Ni; i++)
    {
      for (j = 0; j < BlokDim; j++)
      {
        T->VRe[j][i] = 0;
        T->VRe[i][j] = 0;
      }
    }
    for (i = 0; i < BlokDim; i++)
    {
      T->VRe[i][i] = 1;
    }
    strcpy(T->Name, "Eigenvectors of Minor of ");
    strncat(T->Name, M.Name, Max(0, MAXStrLength - 30));
    i = T->Ni - BlokDim;
    dsyev(jobz, uplo, &i, &(T->VRe[BlokDim][BlokDim]), &M.Ni, val, work, &lwork, &info);
    if (info)
    {
      printf("dseyv failed in MatrixToAndersonMatrix with info value %i\n", info);
      printf("While working on matrix %s\n", M.Name);
      printf("continue diagonalization with dsyvd\n");
      fflush(stdout);
      dsyevd(jobz, uplo, &i, &(T->VRe[BlokDim][BlokDim]), &M.Ni, val, work, &lwork, iwork, &liwork, &info);
      if (info)
      {
        printf("also dsyevd failed in MatrixToAndersonMatrix with info value %i\n", info);
        printf("ERORR!!!!!!!!!!\n");
        fflush(stdout);
        return 1;
      }
    }
    if (MatrixTransform(&U, M, *T))
      return 1;
    for (i = BlokDim; i < T->Ni; i++)
    {
      phase.Re = 0;
      for (j = 0; j < BlokDim; j++)
      {
        phase.Re += U.VRe[j][i];
      }
      if (phase.Re < 0)
      {
        for (j = 0; j < T->Ni; j++)
        {
          T->VRe[i][j] *= -1;
        }
      }
    }
  }
  else
  {                             /* Complex */
    Complex16 MReIm[M.Ni * M.Ni];
    int    lcwork = 2 * M.Ni + Sqr(M.Ni);
    Complex16 cwork[lcwork];
    T->Ni = M.Ni;
    T->Nj = M.Nj;
    T->QComplex = M.QComplex;
    strcpy(T->Name, "Eigenvectors of Minor of ");
    strncat(T->Name, M.Name, Max(0, MAXStrLength - 30));
    if (InitMatrix(T))
    {
      FreeFrameBlokAndersonMatrix(H);
      return 1;
    }
    for (i = 0; i < M.Ni * M.Ni; i++)
    {
      MReIm[i].Re = M.VRe[0][i];
      MReIm[i].Im = M.VIm[0][i];
    }
    for (i = 0; i < T->Ni; i++)
    {
      for (j = 0; j < BlokDim; j++)
      {
        MReIm[i + M.Ni * j].Re = 0;
        MReIm[i + M.Ni * j].Im = 0;
        MReIm[i * M.Ni + j].Re = 0;
        MReIm[i * M.Ni + j].Im = 0;
      }
    }
    for (i = 0; i < BlokDim; i++)
    {
      MReIm[i * (M.Ni + 1)].Re = 1;
      MReIm[i * (M.Ni + 1)].Im = 1;
    }
    i = T->Ni - BlokDim;
    zheev(jobz, uplo, &i, &(MReIm[BlokDim + M.Ni * BlokDim]), &M.Ni, val, cwork, &lcwork, work, &info);
    if (info != 0)
    {
      printf("zheev failed in Eigensystem with info value %i\n", info);
      printf("While working on matrix %s\n", M.Name);
      printf("continue diagonalization with zheevd\n");
      fflush(stdout);
      zheevd(jobz, uplo, &i, &(MReIm[BlokDim + M.Ni * BlokDim]), &M.Ni, val, cwork, &lcwork, work, &lwork, iwork,
             &liwork, &info);
      if (info != 0)
      {
        printf("also zheevd failed in Eigensystem with info value %i\n", info);
        printf("ERORR!!!!!!!!!!\n");
        fflush(stdout);
        return 1;
      }
    }
    for (i = 0; i < M.Ni * M.Ni; i++)
    {
      T->VRe[0][i] = MReIm[i].Re;
      T->VIm[0][i] = -MReIm[i].Im;
    }
    if (MatrixTransform(&U, M, *T))
      return 1;
    for (i = BlokDim; i < T->Ni; i++)
    {
      cwork[0].Re = 0;
      cwork[0].Im = 0;
      for (j = 0; j < BlokDim; j++)
      {
        cwork[0].Re += U.VRe[j][i];
        cwork[0].Im += U.VIm[j][i];
      }
      if (fabs(cwork[0].Re) > fabs(cwork[0].Im))
      {
        phase.Im = 1;
        phase.Re = -cwork[0].Re / cwork[0].Im;
      }
      else
      {
        phase.Re = 1;
        phase.Im = -cwork[0].Im / cwork[0].Re;
      }
      work[0] = sqrt(Sqr(phase.Re) + Sqr(phase.Im));
      phase.Im /= work[0];
      phase.Re /= work[0];
      if ((cwork[0].Re * phase.Re - cwork[0].Im * phase.Im) < 0)
      {
        phase.Re *= -1;
        phase.Im *= -1;
      }
      for (j = 0; j < T->Ni; j++)
      {
        work[0] = T->VRe[i][j] * phase.Re - T->VIm[i][j] * phase.Im;
        T->VIm[i][j] = T->VRe[i][j] * phase.Im + T->VIm[i][j] * phase.Re;
        T->VRe[i][j] = work[0];
      }
    }
  }
  /* create the blok anderson matrix H by transforming M with the use of T */
  FreeMatrix(&U);
  for (i = 0; i < H->N; i++)
  {
    if (MatrixTransformPart
        (&(H->a[i]), M, *T, i * BlokDim, Min(M.Ni, (i + 1) * (BlokDim) - 1), i * BlokDim,
         Min(M.Nj, (i + 1) * (BlokDim) - 1)))
    {
      printf("Unexpected out of memory in BlokTriDiagonalize 30\n");
      return 1;
    }
  }
  for (i = 0; i < H->N - 1; i++)
  {
    if (MatrixTransformPart(&(H->b[i]), M, *T, 0, (BlokDim) - 1, (i + 1) * BlokDim, Min(M.Nj, (i + 2) * (BlokDim) - 1)))
    {
      printf("Unexpected out of memory in BlokTriDiagonalize 31\n");
      return 1;
    }
  }
  return 0;
}


/* 
*
*/
uint32_t MatrixToDoubleTriDiagonalMatrixNoOverlap(MatrixType const M, MatrixType * const Tval, MatrixType * const Tcon, DoubleTriDiagonalMatrixType * const H, double const Ef)
{
  int    i, iEf;
  MatrixType TMtoAnd;
  AndersonMatrixType MAnd;
  if (MatrixToAndersonMatrix(M, &TMtoAnd, &MAnd))
    return 1;
  iEf = M.Ni;
  /* val goes from 1 to iEf-1 */
  /* con goes from iEf to M.Ni-1 */
  for (i = 1; i < M.Ni; i++)
  {
    if (MAnd.a[i] > Ef)
    {
      iEf = i;
      i = M.Ni;
    }
  }


  //TODO: general coding

  return 0;
}


uint32_t MatrixToBlokDoubleTriDiagonalMatrixNoOverlap(MatrixType const M, MatrixType * const Tval,
                                                      MatrixType * const Tcon,
                                                      BlokDoubleTriDiagonalMatrixType * const H, const uint32_t BlokDim,
                                                      double const Ef)
{
  //TODO: General coding
  return 0;
}


/****************************************************************************/
/********************         Matrix conversions         ********************/
/********************         ListOfPoles To ...         ********************/
/****************************************************************************/

/*
 *
 */
uint32_t ListOfPolesToTridiagonalMatrixNonOrthogonal(const ListOfPoles A, TriDiagonalMatrixType * const M,
                                                     uint32_t NTriDiag)
{
  uint32_t i;
  AndersonMatrixType AM;
  AM.N = A.N + 1;
  if (InitAndersonMatrix(&AM))
  {
    printf("InitAndersonMatrix failed in ListOfPolesToTridiagonalMatrix\n;");
    fflush(stdout);
    return 1;
  };
  AM.a[0] = 0;
  AM.a[1] = 0;
  AM.b[0] = 1;
  for (i = 0; i < A.N; i++)
  {
    AM.b[i + 1] = sqrt(fabs(A.b[i]));
    AM.a[i + 2] = A.a[i];
  }
  AndersonMatrixToTriDiagonalMatrixNonOrthogonal(AM, M, NTriDiag + 1);
  for (i = 0; i < M->N - 1; i++)
  {
    M->a[i] = M->a[i + 1];
    M->b[i] = M->b[i + 1];
  }
  M->a[M->N - 1] = M->a[M->N];
  M->N--;
  FreeAndersonMatrix(&AM);
  return 0;
}

/*
 *
 */
uint32_t ListOfPolesToTridiagonalMatrix(const ListOfPoles A, TriDiagonalMatrixType * const M, uint32_t NTriDiag)
{
  uint32_t i;
  AndersonMatrixType AM;
  AM.N = A.N + 1;
  if (InitAndersonMatrix(&AM))
  {
    printf("InitAndersonMatrix failed in ListOfPolesToTridiagonalMatrix\n;");
    fflush(stdout);
    return 1;
  };
  AM.a[0] = 0;
  AM.a[1] = A.a0;
  AM.b[0] = 1;
  for (i = 0; i < A.N; i++)
  {
    AM.b[i + 1] = sqrt(fabs(A.b[i]));
    AM.a[i + 2] = A.a[i];
  }
  AndersonMatrixToTriDiagonalMatrix(AM, M, NTriDiag + 1);
  for (i = 0; i < M->N - 1; i++)
  {
    M->a[i] = M->a[i + 1];
    M->b[i] = M->b[i + 1];
  }
  M->a[M->N - 1] = M->a[M->N];
  M->N--;
  FreeAndersonMatrix(&AM);
  return 0;
}

/*
 *
 */
uint32_t ListOfPolesToAndersonMatrix(const ListOfPoles A, AndersonMatrixType * const M)
{
  TriDiagonalMatrixType Atri;
  ListOfPolesToTridiagonalMatrix(A, &Atri, A.N);
  TriDiagonalMatrixToAndersonMatrix(Atri, M);
  FreeTriDiagonalMatrix(&Atri);
  M->mu = A.mu;
  return 0;
}

/*
 *
 */
uint32_t ListOfPolesToBADoubleTridiagonalMatrix(const ListOfPoles A, BADoubleTriDiagonalMatrixType * const M, uint32_t NTriDiagval, uint32_t NTriDiagcon, double mu, double Epsilon)
{
  ListOfPoles LPcon, LPval;
  ListOfPolesCopyCon(&LPcon, A, mu, Epsilon);
  ListOfPolesCopyVal(&LPval, A, mu, Epsilon);
  ReduceListOfPolesNumber(&LPcon, NTriDiagcon);
  ReduceListOfPolesNumber(&LPval, NTriDiagval);
  ListOfPolesToTridiagonalMatrix(LPcon, &(M->con), NTriDiagcon);
  ListOfPolesToTridiagonalMatrix(LPval, &(M->val), NTriDiagval);
  M->a0 = A.a0;
  M->mu = mu;
  M->Epsilon = Epsilon;
  SetBADoubleTriDiagonalMatrixConVal(M);
  FreeListOfPoles(&LPval);
  FreeListOfPoles(&LPcon);
  return 0;
}

/*
 *
 */
uint32_t ListOfPolesToBADoubleTridiagonalMatrixOccupation(const ListOfPoles A, BADoubleTriDiagonalMatrixType * const M,
                                                          uint32_t NTriDiagval, uint32_t NTriDiagcon, double mu,
                                                          double occupation)
{
  double muOccupation;
  ListOfPoles Achop;
  ListOfPoles LPcon, LPval;
  ListOfPolesCopyConOccupation(&LPcon, A, &muOccupation, occupation);
  ListOfPolesCopyValOccupation(&LPval, A, &muOccupation, occupation);
  printf("\nIn ListOfPolesToBADoubleTridiagonalMatrixOccupation mu input is: %22.15E,\n", mu);
  printf("for an occupation of %22.15E      mu needs to be: %22.15E\n", occupation, muOccupation);
  ListOfPolesDetermineOccupation(A, mu, &muOccupation);
  printf("An occupation of     %22.15E is found at mu=      %22.15E\n\n", muOccupation, mu);

  ListOfPoleCopyChoped(A, &Achop, 0.01);
  PrintListOfPoles(Achop);
  FreeListOfPoles(&Achop);

  ReduceListOfPolesNumber(&LPcon, NTriDiagcon);
  ReduceListOfPolesNumber(&LPval, NTriDiagval);
  ListOfPolesToTridiagonalMatrix(LPcon, &(M->con), NTriDiagcon);
  ListOfPolesToTridiagonalMatrix(LPval, &(M->val), NTriDiagval);
  M->a0 = A.a0;
  M->mu = mu;
  M->Epsilon = Abs(mu - muOccupation);
  SetBADoubleTriDiagonalMatrixConVal(M);
  FreeListOfPoles(&LPval);
  FreeListOfPoles(&LPcon);
  return 0;
}

/****************************************************************************/
/********************  Matrix conversions                ********************/
/********************  BADoubleTridiagonalMatrix To ...  ********************/
/****************************************************************************/
