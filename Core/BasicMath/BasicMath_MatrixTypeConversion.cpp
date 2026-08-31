/****************************************************************************/
/********************         Matrix conversions         ********************/
/********************         Poles To ...               ********************/
/****************************************************************************/

/*  Recurrence coefficients of the discrete measure  sum_i Wpole[i] delta(x-Epole[i])
 *
 *  A list of poles, the bath of an Anderson matrix and one branch of a double
 *  tri-diagonal matrix are all this same measure, and turning any of them into
 *  a chain is the same calculation:
 *
 *    sum_i Wpole[i] / (z - Epole[i]) = W / ( z - a[0] - b[0]^2/( z - a[1] - ... ) )
 *
 *  with W = sum_i Wpole[i] the total weight. A caller holding a
 *  TriDiagonalMatrixType H writes sqrt(W) as the hopping into the chain and
 *  passes &H->a[1] and &H->b[1], or &H->a[2] and &H->b[2] when the head of H is
 *  a site of its own.
 *
 *  Epole and Wpole are NPoles energies and NPoles weights, in any order; poles
 *  without weight are ignored, weights are not allowed to be negative. a and b
 *  are supplied by the caller and must hold NMax and NMax-1 elements. On return
 *  *N is the number of chain sites written, which is at most NMax and at most
 *  the number of distinct poles that carry weight, and *W the total weight.
 *
 *  The method is that of Rutishauser, Kahan, Pal and Walker; see Gragg and
 *  Harrod, Numer. Math. 44 (1984) 317. One pole is absorbed at a time by a
 *  sequence of plane rotations, which costs O(NPoles*NMax) operations and
 *  O(NPoles) memory and needs no re-orthogonalization.
 *
 *  The obvious alternative, Lanczos on the diagonal matrix of the poles with
 *  sqrt(Wpole) as the starting vector, is not usable for a long chain. It loses
 *  orthogonality completely, a relative error of 1e-2 on G for 1000 poles,
 *  unless every Lanczos vector is kept and re-orthogonalized against all
 *  previous ones. That costs O(NPoles*NMax^2) operations and O(NPoles*NMax)
 *  memory: 1.5 seconds and 8 Mb for 1000 poles and 100 seconds and 128 Mb for
 *  4000, where this routine needs 7 and 120 milliseconds at the same accuracy
 *  of 1e-14.
 *
 *  Poles at the same energy are merged first. That is not cosmetic: the method
 *  is not robust against a measure whose support is smaller than the number of
 *  poles handed to it, and degenerate bath levels are common.
 */
uint32_t PolesToTriDiagonalCoefficients(uint32_t const NPoles, double const *const Epole, double const *const Wpole,
                                        uint32_t const NMax, double *const a, double *const b,
                                        uint32_t *const N, double *const W)
{
  uint32_t i, j, NKept, NChain;
  double *EW;                   /* the poles as (energy,weight) pairs, so that they can be sorted */
  double *alpha, *betasqr;      /* the running Jacobi matrix; betasqr[0] ends up as the total weight */
  double p, gamma, sigma, sigmaold, t, tnew, rho, tmp, Elam, bmin;
  *N = 0;
  *W = 0;
  if ((NPoles == 0) || (NMax == 0))
  {
    return 0;
  }
  if ((EW = (double*)malloc(2 * NPoles * sizeof(double))) == NULL)
  {
    printf("malloc for EW failed in PolesToTriDiagonalCoefficients\n");
    fflush(stdout);
    return 1;
  }
  NKept = 0;
  for (i = 0; i < NPoles; i++)
  {
    if (Wpole[i] > 0)
    {
      EW[2 * NKept    ] = Epole[i];
      EW[2 * NKept + 1] = Wpole[i];
      NKept++;
    }
  }
  if (NKept == 0)
  {
    free(EW);
    return 0;
  }
  qsort((void *)EW, NKept, 2 * sizeof(double), (compfn) OrderListOfPolesCompare);
  /* merge poles that sit at the same energy, comparing against the first of a
   * run rather than the previous pole so that a long run can not drift */
  j = 0;
  for (i = 1; i < NKept; i++)
  {
    if ((EW[2 * i] - EW[2 * j]) <= 4 * DBL_EPSILON * (Abs(EW[2 * i]) + Abs(EW[2 * j])))
    {
      EW[2 * j + 1] += EW[2 * i + 1];
    }
    else
    {
      j++;
      EW[2 * j    ] = EW[2 * i    ];
      EW[2 * j + 1] = EW[2 * i + 1];
    }
  }
  NKept = j + 1;
  NChain = Min(NMax, NKept);
  if ((alpha = (double*)malloc(NChain * sizeof(double))) == NULL)
  {
    printf("malloc for alpha failed in PolesToTriDiagonalCoefficients\n");
    fflush(stdout);
    free(EW);
    return 1;
  }
  if ((betasqr = (double*)malloc(NChain * sizeof(double))) == NULL)
  {
    printf("malloc for betasqr failed in PolesToTriDiagonalCoefficients\n");
    fflush(stdout);
    free(EW);
    free(alpha);
    return 1;
  }
  for (j = 0; j < NChain; j++)
  {
    alpha[j] = EW[2 * j];
    betasqr[j] = 0;
  }
  betasqr[0] = EW[1];
  for (i = 1; i < NKept; i++)
  {
    Elam = EW[2 * i];
    p = EW[2 * i + 1];
    gamma = 1;
    sigma = 0;
    t = 0;
    for (j = 0; (j < NChain) && (j <= i); j++)
    {
      rho = betasqr[j] + p;
      tmp = gamma * rho;
      sigmaold = sigma;
      if (rho <= 0)
      {
        gamma = 1;
        sigma = 0;
      }
      else
      {
        gamma = betasqr[j] / rho;
        sigma = p / rho;
      }
      tnew = sigma * (alpha[j] - Elam) - gamma * t;
      alpha[j] -= tnew - t;
      t = tnew;
      if (sigma <= 0)
      {
        p = sigmaold * betasqr[j];
      }
      else
      {
        p = t * t / sigma;
      }
      betasqr[j] = tmp;
    }
  }
  *W = betasqr[0];
  /*  the chain ends where the hopping vanishes: what comes after it is
   *  decoupled from the head and does not contribute to the Green's function */
  bmin = sqrt(*W) * 100 * DBL_EPSILON;
  a[0] = alpha[0];
  *N = 1;
  for (j = 1; j < NChain; j++)
  {
    b[j - 1] = sqrt(Abs(betasqr[j]));
    if (b[j - 1] <= bmin)
    {
      b[j - 1] = 0;
      break;
    }
    a[j] = alpha[j];
    (*N)++;
  }
  free(EW);
  free(alpha);
  free(betasqr);
  return 0;
}



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
  Mfull->mu = M.mu;
  Mfull->a0 = M.a[0];
  Mfull->b0 = M.b[0];
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

uint32_t AndersonMatrixToCompactMatrix(AndersonMatrixType * const M, CompactMatrixType * const Mfull)
{
  int    i;
  Mfull->QComplex = 0;
  Mfull->Ni = M->N;
  Mfull->Nj = M->N;
  if (InitCompactMatrix(Mfull))
  {
    printf("Could not initiallize CompactMatrix in AndersonMatrixToCompactMatrix\n");
    return 1;
  };
  Mfull->mu = M->mu;
  Mfull->a0 = M->a[0];
  Mfull->b0 = M->b[0];
  strcpy(Mfull->Name, M->Name);
  for (i = 0; i < M->N; i++)
  {
    Mfull->VRe[i * M->N + i] = M->a[i + 1];
  }
  for (i = 1; i < M->N; i++)
  {
    Mfull->VRe[i * M->N] = M->b[i];
    Mfull->VRe[       i] = M->b[i];
  }
  return 0;
}

uint32_t AndersonMatrixToTriDiagonalMatrixNonOrthogonal(AndersonMatrixType const M, TriDiagonalMatrixType * const H, uint32_t NTriDiag)
{
  uint32_t i;
  double *vec[3];
  if ((vec[0] = (double*)malloc(M.N * sizeof(double))) == NULL)
  {
    printf("malloc 1 failed in AndersonMatrixToTriDiagonalMatrixNonOrthogonal\n");
    fflush(stdout);
    return 1;
  };
  if ((vec[1] = (double*)malloc(M.N * sizeof(double))) == NULL)
  {
    printf("malloc 1 failed in AndersonMatrixToTriDiagonalMatrixNonOrthogonal\n");
    free(vec[0]);
    fflush(stdout);
    return 1;
  };
  if ((vec[2] = (double*)malloc(M.N * sizeof(double))) == NULL)
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
  strncpy(H->Name, M.Name, MAXStrLength);
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
 *  M=T.H.T^T
 *  input M defined, H and T not initiallized
 */
uint32_t AndersonMatrixToTriDiagonalMatrixPlusT(AndersonMatrixType const M, MatrixType * const T, TriDiagonalMatrixType * const H, uint32_t NTriDiag)
{
  if(M.N==0)
  {
    T->QComplex = 0;
    T->Ni = 0;
    T->Nj = 0;
    if (InitMatrix(T))
    {
      printf("InitMatrix failed in AndersonMatrixToTriDiagonalMatrix\n");
      fflush(stdout);
      return 1;
    };
    H->N = 0;
    if (InitTriDiagonalMatrix(H))
    {
      FreeMatrix(T);
      printf("InitTriDiagonalMatrix failed in AndersonMatrixToTriDiagonalMatrix\n");
      return 1;
    }
    strncpy(H->Name, M.Name, MAXStrLength);
    strcpy(T->Name, "TriDiagonalizer of ");
    strncat(T->Name, M.Name, Max(0, MAXStrLength - 20));
    H->a[0] = M.a[0];
    H->mu = M.mu;
    return 0;
  }
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
  strncpy(H->Name, M.Name, MAXStrLength);
  strcpy(T->Name, "TriDiagonalizer of ");
  strncat(T->Name, M.Name, Max(0, MAXStrLength - 20));
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

/*  M = T.H.T^T with T orthogonal and H tri-diagonal; T itself is not returned.
 *
 *  The bath of an Anderson matrix is diagonal by construction: site 0 is the
 *  impurity with energy a[1], coupled to the bath levels a[2...N] by the
 *  hoppings b[1...N-1]. The chain is therefore the Jacobi matrix of the bath as
 *  seen from the impurity and no eigenvalue decomposition is needed. H->a[1] is
 *  the impurity, H->b[1] the hopping into the chain, and H->a[2...] / H->b[2...]
 *  the chain itself.
 *
 *  NTriDiag is an upper bound on the dimension of H. The chain is shorter when
 *  the bath has fewer distinct levels than that, or when a hopping vanishes.
 */
uint32_t AndersonMatrixToTriDiagonalMatrix(AndersonMatrixType const M, TriDiagonalMatrixType * const H, uint32_t NTriDiag)
{
  uint32_t i, NAlloc, NBath, NChain;
  double *Ebath, *Wbath, Whyb;
  NAlloc = Min(NTriDiag, M.N);
  H->N = NAlloc;
  if (InitTriDiagonalMatrix(H))
  {
    printf("InitTriDiagonalMatrix failed in AndersonMatrixToTriDiagonalMatrix\n");
    fflush(stdout);
    return 1;
  }
  strncpy(H->Name, M.Name, MAXStrLength);
  H->mu = M.mu;
  H->a[0] = M.a[0];
  if (NAlloc == 0)
  {
    return 0;
  }
  H->b[0] = M.b[0];
  H->a[1] = M.a[1];
  H->N = 1;
  NBath = M.N - 1;
  if ((NBath == 0) || (NAlloc == 1))
  {
    return 0;
  }
  if ((Ebath = (double*)malloc(NBath * sizeof(double))) == NULL)
  {
    printf("malloc for Ebath failed in AndersonMatrixToTriDiagonalMatrix\n");
    fflush(stdout);
    return 1;
  }
  if ((Wbath = (double*)malloc(NBath * sizeof(double))) == NULL)
  {
    printf("malloc for Wbath failed in AndersonMatrixToTriDiagonalMatrix\n");
    fflush(stdout);
    free(Ebath);
    return 1;
  }
  for (i = 0; i < NBath; i++)
  {
    Ebath[i] = M.a[i + 2];
    Wbath[i] = Sqr(M.b[i + 1]);
  }
  if (PolesToTriDiagonalCoefficients(NBath, Ebath, Wbath, NAlloc - 1, &(H->a[2]), &(H->b[2]), &NChain, &Whyb))
  {
    printf("PolesToTriDiagonalCoefficients failed in AndersonMatrixToTriDiagonalMatrix\n");
    fflush(stdout);
    free(Ebath);
    free(Wbath);
    return 1;
  }
  free(Ebath);
  free(Wbath);
  if (NChain != 0)
  {
    H->b[1] = sqrt(Whyb);
    H->N = 1 + NChain;
  }
  return 0;
}

uint32_t AndersonMatrixToDoubleTriDiagonalMatrixNonOrthogonal(AndersonMatrixType const M, DoubleTriDiagonalMatrixType * const H, double const Ef, double const epsilon)
{
  /* To DO */
  return 1;
}

/*  Kept for its history and to be reimplemented; not called by anything.
 *
 *  This is the version that returned the transformation matrix T. It follows the
 *  same convention for the shared impurity site as everything else, half of the
 *  energy in each head slot, but it still differs from
 *  AndersonMatrixToDoubleTriDiagonalMatrix below: it splits with a linear ramp
 *  over [mu-epsilon,mu+epsilon] rather than with a Fermi function, it invents a
 *  bath site with zero coupling when a branch would be empty, and it
 *  tri-diagonalizes with re-orthogonalized Lanczos, which costs the O(N^2)
 *  memory that PolesToTriDiagonalCoefficients exists to avoid. Use the routine
 *  below until this one is rewritten on top of it.
 */
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
*  TODO
*  epsilon is replaced in name by kBT, not yet changed the implementation
*/
uint32_t AndersonMatrixToDoubleTriDiagonalMatrixPlusT(AndersonMatrixType const M, MatrixType * const T, DoubleTriDiagonalMatrixType * const H, double const mu, double const epsilon)
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
  /*  the impurity site is shared by the two chains and its energy is the sum of
   *  the two head slots, so each chain is built with half of it */
  MV.a[1] = 0.5 * M.a[1];
  MC.a[1] = 0.5 * M.a[1];
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
  // TODO implementation needs to be changed
  H->kBT = epsilon;
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


/*  Split the bath of an Anderson matrix at the chemical potential and put the
 *  two halves in the two chains of a double tri-diagonal matrix,
 *
 *    G(w) = a0 + b0^2 / ( w - Gval(w) - Gcon(w) )
 *
 *  The bath level a[i+1] couples to the impurity with b[i] and therefore
 *  carries a weight b[i]^2 in the hybridization. The Fermi function divides
 *  that weight over the two chains, f to the valence and 1-f to the conduction
 *  chain, so the two chains together reproduce the original function exactly.
 *  FermiFunction returns exactly 1 below mu-40 kBT and exactly 0 above
 *  mu+40 kBT, so a level outside that window ends up in one chain only and the
 *  levels inside it are the only ones that appear twice. The chains are as long
 *  as they need to be to be exact; nothing is truncated here.
 *
 *  The impurity energy a[1] is shared by the two chains: val.a[0] + con.a[0] is
 *  what enters the denominator, and each chain is given half of it. How it is
 *  split between the two is a gauge with no canonical value, so equal halves is
 *  a convention and not an invariant.
 *
 *  A chain can come out empty, with mu below or above the whole bath. That is
 *  val.N = 0 with val.b[0] = 0, and the impurity keeps its half of a[1].
 */
uint32_t AndersonMatrixToDoubleTriDiagonalMatrix(AndersonMatrixType const M, DoubleTriDiagonalMatrixType * const H, double const mu, double const kBT)
{
  uint32_t i, NBath, NVal, NCon;
  double *Ebath, *Wval, *Wcon, Wvaltot, Wcontot, x, f, g;
  if (M.N == 0)
  {
    printf("AndersonMatrixToDoubleTriDiagonalMatrix called with a matrix of dimension 0\n");
    fflush(stdout);
    return 1;
  }
  NBath = M.N - 1;
  /*  each chain can hold at most one site per bath level, so both are allocated
   *  for NBath sites here and their length set to what the Lanczos returns */
  H->val.N = NBath;
  H->con.N = NBath;
  if (InitDoubleTriDiagonalMatrix(H))
  {
    printf("InitDoubleTriDiagonalMatrix failed in AndersonMatrixToDoubleTriDiagonalMatrix\n");
    fflush(stdout);
    return 1;
  }
  strncpy(H->Name, M.Name, MAXStrLength);
  H->mu = mu;
  H->kBT = kBT;
  H->a0 = M.a[0];
  H->b0 = M.b[0];
  H->val.a[0] = 0.5 * M.a[1];
  H->con.a[0] = 0.5 * M.a[1];
  H->val.N = 0;
  H->con.N = 0;
  if (NBath == 0)
  {
    return 0;
  }
  if ((Ebath = (double*)malloc(NBath * sizeof(double))) == NULL)
  {
    printf("malloc for Ebath failed in AndersonMatrixToDoubleTriDiagonalMatrix\n");
    fflush(stdout);
    return 1;
  }
  if ((Wval = (double*)malloc(NBath * sizeof(double))) == NULL)
  {
    printf("malloc for Wval failed in AndersonMatrixToDoubleTriDiagonalMatrix\n");
    fflush(stdout);
    free(Ebath);
    return 1;
  }
  if ((Wcon = (double*)malloc(NBath * sizeof(double))) == NULL)
  {
    printf("malloc for Wcon failed in AndersonMatrixToDoubleTriDiagonalMatrix\n");
    fflush(stdout);
    free(Ebath);
    free(Wval);
    return 1;
  }
  for (i = 0; i < NBath; i++)
  {
    Ebath[i] = M.a[i + 2];
    x = M.a[i + 2] - mu;
    f = FermiFunction(M.a[i + 2], mu, kBT);
    /*  1-f loses its relative accuracy where f is close to one, so the
     *  conduction weight is formed from the Fermi function of -x instead,
     *  except where FermiFunction itself returns an exact value */
    if ((x == 0) || (kBT <= 0) || (Abs(x) >= 40 * kBT))
    {
      g = 1.0 - f;
    }
    else
    {
      g = 1.0 / (1.0 + exp(-x / kBT));
    }
    Wval[i] = Sqr(M.b[i + 1]) * f;
    Wcon[i] = Sqr(M.b[i + 1]) * g;
  }
  if (PolesToTriDiagonalCoefficients(NBath, Ebath, Wval, NBath, &(H->val.a[1]), &(H->val.b[1]), &NVal, &Wvaltot))
  {
    printf("PolesToTriDiagonalCoefficients failed for val in AndersonMatrixToDoubleTriDiagonalMatrix\n");
    fflush(stdout);
    free(Ebath);
    free(Wval);
    free(Wcon);
    return 1;
  }
  if (PolesToTriDiagonalCoefficients(NBath, Ebath, Wcon, NBath, &(H->con.a[1]), &(H->con.b[1]), &NCon, &Wcontot))
  {
    printf("PolesToTriDiagonalCoefficients failed for con in AndersonMatrixToDoubleTriDiagonalMatrix\n");
    fflush(stdout);
    free(Ebath);
    free(Wval);
    free(Wcon);
    return 1;
  }
  free(Ebath);
  free(Wval);
  free(Wcon);
  H->val.b[0] = sqrt(Wvaltot);
  H->val.N = NVal;
  H->con.b[0] = sqrt(Wcontot);
  H->con.N = NCon;
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
uint32_t AndersonMatrixToEigensystemOnsiteShift(AndersonMatrixType const M, double const shift, MatrixType * const T, double *val)
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

uint32_t AndersonMatrixToEigenSystemCompact(const AndersonMatrixType * const M, double *val, CompactMatrixType * const T)
{
  //TODO all
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
  A->a0 = M.a[0];
  strncpy(A->Name, M.Name, MAXStrLength);
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
uint32_t AndersonMatrixFindOnsiteEnergyShift(AndersonMatrixType const M, double const occupation, double *const onsiteshift)
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

uint32_t EigenstateToTridiagonalMatrix(double const *psi0, double const *val, uint32_t const N, TriDiagonalMatrixType * const H)
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
uint32_t AndersonMatrixToBADoubleTriDiagonalMatrixSimple(AndersonMatrixType const M, BADoubleTriDiagonalMatrixType * const A)
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
 *  Transform a tridiagonal matrix to Natural orbital format using chemical potential Ef / mu with the soft edge given by
 *  the interval Ef-Epsilon to Ef+Epsilon, scaled by ( Epole - (Ef-Epsilon) ) / (2 Epsilon))
 *  Reduce the number of poles in for the valence and conduction part to NTriDiagVal and NTriDiagCon
 *
 *  Go via ListOfPoles
 */
uint32_t AndersonMatrixToBADoubleTriDiagonalMatrix(AndersonMatrixType const M, BADoubleTriDiagonalMatrixType * const A, uint32_t NTriDiagVal, uint32_t NTriDiagCon, double Ef, double Epsilon)
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
  strcpy(A->Name, M.Name);
  return 0;
}

/*
 *
 */
uint32_t AndersonMatrixToBADoubleTriDiagonalMatrixOccupation(AndersonMatrixType const M, double const onsiteshift, BADoubleTriDiagonalMatrixType * const A, double *const occupation)
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
uint32_t AndersonMatrixToBADoubleTriDiagonalMatrixTemperature(AndersonMatrixType const M, double const onsiteshift, double const kbT, BADoubleTriDiagonalMatrixType * const A, double *const occupation)
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
uint32_t DetermineMuFromEigensystem(MatrixType const Eigenfun, double const *Eigenval, double const Nelec, double const kbT, double *const mu)
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
uint32_t AndersonMatrixToBADoubleTriDiagonalMatrixWrong(AndersonMatrixType const M,BADoubleTriDiagonalMatrixType * const A, uint32_t NTriDiagVal,uint32_t NTriDiagCon, double Ef, double Epsilon)
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
/********************         TriDiagonalMatrix To ...   ********************/
/****************************************************************************/

uint32_t TriDiagonalMatrixToEigenvalues(TriDiagonalMatrixType const M, double *val)
{
  int    i, info;
  const char *jobz = "N";
  double e[M.N];
  uint32_t    lwork = (1 + 4 * M.N + Sqr(M.N));
  uint32_t    liwork = 3 + 5 * M.N;
  double *work,*Z=NULL;
  int   *iwork;
  work = (double*)malloc(lwork * sizeof(double));
  if (work == NULL)
  {
    printf("malloc of work failed in TriDiagonalMatrixToEigensystem\n");
    return 1;
  };
  iwork = (int*)malloc(liwork * sizeof(int));
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
  if(M.N==0)
  {
    T->Ni = 0;
    T->Nj = 0;
    T->QComplex = 0;
    strcpy(T->Name, "Eigenvectors of ");
    strncat(T->Name, M.Name, Max(0, MAXStrLength - 18));
    if (InitMatrix(T))
    {
      printf("InitMatrix(T) failed in TriDiagonalMatrixToEigensystem\n");
      fflush(stdout);
      return 1;
    }
    return 0;
  }
  int    i, info;
  const char *jobz = "V";
  double e[M.N];
  uint32_t    lwork = (1 + 4 * M.N + Sqr(M.N));
  uint32_t    liwork = 3 + 5 * M.N;
  double *work;
  int   *iwork;
  work = (double*)malloc(lwork * sizeof(double));
  if (work == NULL)
  {
    printf("malloc of work failed in TriDiagonalMatrixToEigensystem\n");
    return 1;
  };
  iwork = (int*)malloc(liwork * sizeof(int));
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

uint32_t TriDiagonalMatrixToEigensystemCompact(const TriDiagonalMatrixType * const M, double *val, CompactMatrixType * const T)
{
  int    i, info;
  const char *jobz = "V";
  double e[M->N];
  uint32_t    lwork = (1 + 4 * M->N + Sqr(M->N));
  uint32_t    liwork = 3 + 5 * M->N;
  double *work;
  int   *iwork;
  work = (double*)malloc(lwork * sizeof(double));
  if (work == NULL)
  {
    printf("malloc of work failed in TriDiagonalMatrixToEigensystem\n");
    return 1;
  };
  iwork = (int*)malloc(liwork * sizeof(int));
  if (iwork == NULL)
  {
    printf("malloc of iwork failed in TriDiagonalMatrixToEigensystem\n");
    return 1;
  };
  
  for (i = 0; i < M->N - 1; i++) e[i] = M->b[i + 1];
  e[M->N - 1] = 0.0;
  for (i = 0; i < M->N; i++) val[i] = M->a[i + 1];
  
  T->Ni = M->N;
  T->Nj = M->N;
  T->QComplex = 0;
  if (InitCompactMatrix(T))
  {
    printf("InitCompactMatrix(T) failed in TriDiagonalMatrixToEigensystemCompact\n");
    fflush(stdout);
    free(work);
    free(iwork);
    return 1;
  }
  strcpy(T->Name, "Eigenvectors of ");
  strncat(T->Name, M->Name, Max(0, MAXStrLength - 18));
  dstevd(jobz, &M->N, val, e, T->VRe, &M->N, work, &lwork, iwork, &liwork, &info);
  if (info != 0)
  {
    printf("dstevd failed in TriDiagonalMatrixToEigensystemCompact with info value %i\n", info);
    printf("While working on matrix %s\n", M->Name);
    printf("continue diagonalization with dstev\n");
    fflush(stdout);
    dstev(jobz, &M->N, val, e, T->VRe, &M->N, work, &info);
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
  
  A->mu = M.mu;
  strncpy(A->Name, M.Name, MAXStrLength);
  return 0;
}


uint32_t TriDiagonalMatrixToAndersonMatrixPlusT(TriDiagonalMatrixType const M, MatrixType * const T, AndersonMatrixType * const A)
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

uint32_t TriDiagonalMatrixToDoubleTriDiagonalMatrix(TriDiagonalMatrixType const M, DoubleTriDiagonalMatrixType * const A, double Ef, double epsilon)
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


uint32_t TriDiagonalMatrixToSecondDoubleTriDiagonalMatrix(TriDiagonalMatrixType const M, SecondDoubleTriDiagonalMatrixType * const A, double Ef, double epsilon)
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

/*  Kept for its history and to be reimplemented; not called by anything outside
 *  Test/. It routes through AndersonMatrixToDoubleTriDiagonalMatrixPlusT and so
 *  inherits every difference listed there.
 */
uint32_t TriDiagonalMatrixToDoubleTriDiagonalMatrixPlusT(TriDiagonalMatrixType const M, MatrixType * const T, DoubleTriDiagonalMatrixType * const A, double Ef,  double epsilon)
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

uint32_t TriDiagonalMatrixToSecondDoubleTriDiagonalMatrixPlusT(TriDiagonalMatrixType const M, MatrixType * const T,  SecondDoubleTriDiagonalMatrixType * const A, double Ef, double epsilon)
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


uint32_t TriDiagonalMatrixToCompactMatrix(const TriDiagonalMatrixType * const M, CompactMatrixType * const A)
{
  uint32_t i;
  A->QComplex = 0;
  A->Ni = M->N;
  A->Nj = M->N;
  if (InitCompactMatrix(A))
  {
    printf("InitMatrix failed in TriDiagonalMatrixToMatrix\n");
    return 1;
  }
  A->a0 = M->a[0];
  if(M->N!=0)
  {
    A->b0 = M->b[0];
  }
  else
  {
    A->b0 = 0;
  }
  A->mu = M->mu;
  for (i = 1; i < M->N; i++)
  {
    A->VRe[(i - 1)*M->N + i - 1] = M->a[i];
    A->VRe[(i    )*M->N + i - 1] = M->b[i];
    A->VRe[(i - 1)*M->N + i    ] = M->b[i];
  }
  if(M->N!=0) A->VRe[(M->N - 1) * M->N + M->N - 1] = M->a[M->N];
  return 0;
}



/*
 *  Transform a tridiagonal matrix to Natural orbital format using chemical potential Ef / mu with the soft edge given by
 *  the interval Ef-Epsilon to Ef+Epsilon, scaled by ( Epole - (Ef-Epsilon) ) / (2 Epsilon))
 *  Reduce the number of poles in for the valence and conduction part to NTriDiagVal and NTriDiagCon
 *
 *  Go via ListOfPoles
 */
uint32_t TriDiagonalMatrixToBADoubleTriDiagonalMatrix(TriDiagonalMatrixType const M, BADoubleTriDiagonalMatrixType * const A, uint32_t NTriDiagVal, uint32_t NTriDiagCon, double Ef, double Epsilon)
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
uint32_t TriDiagonalMatrixToBADoubleTriDiagonalMatrixOccupation(TriDiagonalMatrixType const M,BADoubleTriDiagonalMatrixType * const A,uint32_t NTriDiagVal, uint32_t NTriDiagCon, double Ef,double Occupation)
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

/*  Convert a tridiagonal matrix to a Blocktridiagonal matrix
 *
 */
uint32_t TriDiagonalMatrixToBlockTriDiagonalMatrix(TriDiagonalMatrixType const * const Min, BlockTriDiagonalMatrixType * Mout)
{
  Mout->BlockSize = 1;
  Mout->N = Min->N;
  Mout->QComplex = 0;
  if(InitBlockTridiagonalMatrix(Mout)) return 1;
  memcpy(Mout->Name, Min->Name, sizeof(char)*MAXStrLength);
  memcpy(Mout->aRe, Min->a, sizeof(double)*(Min->N+1));
  memcpy(Mout->bRe, Min->b, sizeof(double)*(Min->N));
  return 0;
}
uint32_t ConvertTriDiagonalMatrixToBlockTriDiagonalMatrixNoCopy(TriDiagonalMatrixType * Min, BlockTriDiagonalMatrixType * Mout)
{
  Mout->BlockSize = 1;
  Mout->N = Min->N;
  Mout->QComplex = 0;
  if(InitBlockTridiagonalMatrix(Mout)) return 1;
  memcpy(Mout->Name, Min->Name, sizeof(char)*MAXStrLength);
  Mout->aRe = Min->a;
  Min->a=NULL;
  Mout->bRe = Min->b;
  Min->b=NULL;
  return 0;
}


/********************************************************************************/
/********************         Matrix conversions             ********************/
/********************         DoubleTriDiagonalMatrix To ... ********************/
/********************************************************************************/

/*  A double tri-diagonal matrix is a band matrix of half bandwidth 2
 *
 *  Tri-diagonalizing G from the impurity has its first step for free: the format
 *  already gives a[1] = val.a[0]+con.a[0] and b[1] = sqrt(val.b[0]^2+con.b[0]^2).
 *  What is left is the minor, two chains that do not touch, seen from
 *  (val.b[0] val_1 + con.b[0] con_1)/b[1] -- the Jacobi matrix of a sum of two
 *  measures.
 *
 *  Read the matrix as val_N ... val_1, impurity, con_1 ... con_N and it is
 *  tri-diagonal, but that puts the impurity in the middle and Lanczos has to
 *  start at it. Interleave the two chains behind the impurity instead,
 *
 *      impurity, val_1, con_1, val_2, con_2, val_3, con_3, ...
 *          1       2      3      4      5      6      7
 *
 *  and every coupling sits at distance 1 or 2: impurity to val_1 and to con_1,
 *  then val_k to val_k+1 and con_k to con_k+1. So the whole matrix is a symmetric
 *  band matrix of half bandwidth 2, block tri-diagonal with 2x2 blocks, and the
 *  starting vector is e_1 with no rotation needed. Chains of unequal length are
 *  no different: once the shorter one runs out the rest of the longer one is
 *  consecutive, at distance 1.
 *
 *  dsbtrd then reduces it with Givens rotations, annihilating only in planes 2
 *  and beyond, so it leaves e_1 alone and what comes out is the Lanczos chain
 *  from the impurity -- computed by orthogonal similarity rather than by the
 *  three term recurrence. That is the point: the recurrence loses orthogonality
 *  and needs every vector kept and re-orthogonalized, O(N^2) of memory, while a
 *  band reduction is backward stable by construction, T = Q^T (A+E) Q with
 *  |E| = O(eps |A|), and needs 4N doubles. It also costs O(N^2 kd) rather than
 *  the two eigenvalue problems a route through the Anderson form would need.
 *
 *  The signs of the off diagonals are a gauge -- flipping a basis vector flips
 *  the two hoppings that touch it, and taking them in turn makes all of them
 *  positive -- so the result is made positive.
 */
uint32_t DoubleTriDiagonalMatrixToTriDiagonalMatrix(DoubleTriDiagonalMatrixType const M, TriDiagonalMatrixType * const H)
{
  char const *VECT = "N";
  char const *UPLO = "U";
  const uint32_t KD = 2;
  const uint32_t LDAB = 3;
  uint32_t N, k, m, p, q;
  double *AB, *WORK;
  int INFO = 0;
  N = 1 + M.val.N + M.con.N;
  H->N = N;
  if (InitTriDiagonalMatrix(H))
  {
    printf("InitTriDiagonalMatrix failed in DoubleTriDiagonalMatrixToTriDiagonalMatrix\n");
    fflush(stdout);
    return 1;
  }
  strncpy(H->Name, M.Name, MAXStrLength);
  H->mu   = M.mu;
  H->a[0] = M.a0;
  H->b[0] = M.b0;
  H->a[1] = M.val.a[0] + M.con.a[0];   /* the impurity site is shared, its energy is the sum */
  if (N == 1)
  {
    return 0;                          /* both chains empty, only the impurity is left */
  }
  if ((M.val.N == 0) || (M.con.N == 0))
  {
    /*  one chain empty: the minor is a single chain and already tri-diagonal,
     *  so the answer is the impurity in front of it */
    TriDiagonalMatrixType const * const C = (M.val.N == 0) ? &(M.con) : &(M.val);
    H->b[1] = fabs(C->b[0]);
    for (k = 1; k <= C->N; k++)
    {
      H->a[1 + k] = C->a[k];
    }
    for (k = 1; k < C->N; k++)
    {
      H->b[1 + k] = fabs(C->b[k]);
    }
    return 0;
  }
  if ((AB = (double*)calloc(LDAB * N, sizeof(double))) == NULL)
  {
    printf("calloc for AB failed in DoubleTriDiagonalMatrixToTriDiagonalMatrix\n");
    fflush(stdout);
    return 1;
  }
  if ((WORK = (double*)malloc(N * sizeof(double))) == NULL)
  {
    printf("malloc for WORK failed in DoubleTriDiagonalMatrixToTriDiagonalMatrix\n");
    fflush(stdout);
    free(AB);
    return 1;
  }
  /*  band storage, UPLO='U': the element A(i,j) with i<=j of the fortran matrix
   *  lives at AB(KD+1+i-j,j), which is AB[(KD+i-j) + (j-1)*LDAB] here. Site
   *  numbers are 1 based, val_k at 2k and con_k at 2k+1 while both chains last,
   *  and consecutive afterwards */
  m = Min(M.val.N, M.con.N);
  #define DoubleTriPosVal(k) (((k) <= m) ? (2 * (k)) : (m + (k) + 1))
  #define DoubleTriPosCon(k) (((k) <= m) ? (2 * (k) + 1) : (m + (k) + 1))
  #define DoubleTriBandSet(i,j,v) AB[(KD + (i) - (j)) + ((j) - 1) * LDAB] = (v)
  DoubleTriBandSet(1, 1, M.val.a[0] + M.con.a[0]);
  for (k = 1; k <= M.val.N; k++)
  {
    p = DoubleTriPosVal(k);
    DoubleTriBandSet(p, p, M.val.a[k]);
  }
  for (k = 1; k <= M.con.N; k++)
  {
    p = DoubleTriPosCon(k);
    DoubleTriBandSet(p, p, M.con.a[k]);
  }
  DoubleTriBandSet(1, DoubleTriPosVal(1), M.val.b[0]);
  DoubleTriBandSet(1, DoubleTriPosCon(1), M.con.b[0]);
  for (k = 1; k < M.val.N; k++)
  {
    p = DoubleTriPosVal(k);
    q = DoubleTriPosVal(k + 1);
    DoubleTriBandSet(p, q, M.val.b[k]);
  }
  for (k = 1; k < M.con.N; k++)
  {
    p = DoubleTriPosCon(k);
    q = DoubleTriPosCon(k + 1);
    DoubleTriBandSet(p, q, M.con.b[k]);
  }
  #undef DoubleTriPosVal
  #undef DoubleTriPosCon
  #undef DoubleTriBandSet
  dsbtrd(VECT, UPLO, &N, &KD, AB, &LDAB, &(H->a[1]), &(H->b[1]), NULL, &N, WORK, &INFO);
  free(AB);
  free(WORK);
  if (INFO != 0)
  {
    printf("dsbtrd failed in DoubleTriDiagonalMatrixToTriDiagonalMatrix with INFO %i\n", INFO);
    fflush(stdout);
    return 1;
  }
  for (k = 1; k < N; k++)
  {
    H->b[k] = fabs(H->b[k]);
  }
  return 0;
}

/*  The bath of a double tri-diagonal matrix is its two chains, side by side
 *
 *  Read the matrix in the order val_N ... val_1, impurity, con_1 ... con_N and it
 *  is one ordinary tri-diagonal chain with the impurity in the middle. Deleting
 *  the impurity, which is what the Anderson form asks for, therefore leaves two
 *  chains that do not touch: the valence chain and the conduction chain. The
 *  bath is the two spectra side by side, and the impurity reaches eigenvector k
 *  of a chain only through that chain's head, so the hopping is b[0]*c_k with
 *  c_k the first component of that eigenvector.
 *
 *  The same thing falls out of the algebra. EvaluateTriDiagonalMatrix returns
 *  a[0] + b[0]^2 g, with g the continued fraction over a[1...N], so
 *
 *    w - Gval - Gcon = w - (val.a[0]+con.a[0]) - [ val.b[0]^2 g_val + con.b[0]^2 g_con ]
 *
 *  which is the Anderson denominator with the impurity at val.a[0]+con.a[0] and
 *  a hybridization whose poles are the two chain spectra.
 *
 *  So this costs two tri-diagonal eigenvalue problems of size val.N and con.N.
 *  Going through a list of poles instead would need one of size 1+val.N+con.N to
 *  build G's poles, a tri-diagonalization to get back to a chain and a second
 *  eigenvalue problem on its minor, which is about four times the work for
 *  chains of equal length and forms G's poles only to undo them again.
 *
 *  The dimension is preserved, 1+val.N+con.N: this is a rotation of the bath
 *  into its eigenbasis and nothing else. Eigenvectors whose head component
 *  vanishes are kept, with a hopping of zero; Chop removes them.
 */
uint32_t DoubleTriDiagonalMatrixToAndersonMatrix(DoubleTriDiagonalMatrixType const M, AndersonMatrixType * const H)
{
  CompactMatrixType T;
  uint32_t k, NBath;
  double *EB;                   /* the bath as (energy,hopping) pairs, so that it can be sorted */
  NBath = M.val.N + M.con.N;
  H->N = 1 + NBath;
  if (InitAndersonMatrix(H))
  {
    printf("InitAndersonMatrix failed in DoubleTriDiagonalMatrixToAndersonMatrix\n");
    fflush(stdout);
    return 1;
  }
  strncpy(H->Name, M.Name, MAXStrLength);
  H->mu   = M.mu;
  H->a[0] = M.a0;
  H->b[0] = M.b0;
  H->a[1] = M.val.a[0] + M.con.a[0];   /* the impurity site is shared, its energy is the sum */
  if (NBath == 0)
  {
    return 0;
  }
  if (M.val.N != 0)
  {
    if (TriDiagonalMatrixToEigensystemCompact(&(M.val), &(H->a[2]), &T))
    {
      printf("TriDiagonalMatrixToEigensystemCompact failed for val in DoubleTriDiagonalMatrixToAndersonMatrix\n");
      fflush(stdout);
      return 1;
    }
    /*  T.VRe[k*N+j] is component j of eigenvector k. The sign of a hopping is a
     *  gauge on the bath orbital, so only the absolute value is kept */
    for (k = 0; k < M.val.N; k++)
    {
      H->b[1 + k] = fabs(M.val.b[0] * T.VRe[k * M.val.N]);
    }
    FreeCompactMatrix(&T);       /* freed before the second one, so the peak is the larger chain, not the sum */
  }
  if (M.con.N != 0)
  {
    if (TriDiagonalMatrixToEigensystemCompact(&(M.con), &(H->a[2 + M.val.N]), &T))
    {
      printf("TriDiagonalMatrixToEigensystemCompact failed for con in DoubleTriDiagonalMatrixToAndersonMatrix\n");
      fflush(stdout);
      return 1;
    }
    for (k = 0; k < M.con.N; k++)
    {
      H->b[1 + M.val.N + k] = fabs(M.con.b[0] * T.VRe[k * M.con.N]);
    }
    FreeCompactMatrix(&T);
  }
  /*  each half leaves LAPACK in ascending order, but the two halves interleave
   *  wherever the Fermi window put a level in both chains, so order the bath */
  if ((EB = (double*)malloc(2 * NBath * sizeof(double))) == NULL)
  {
    printf("malloc for EB failed in DoubleTriDiagonalMatrixToAndersonMatrix\n");
    fflush(stdout);
    return 1;
  }
  for (k = 0; k < NBath; k++)
  {
    EB[2 * k    ] = H->a[2 + k];
    EB[2 * k + 1] = H->b[1 + k];
  }
  qsort((void *)EB, NBath, 2 * sizeof(double), (compfn) OrderListOfPolesCompare);
  for (k = 0; k < NBath; k++)
  {
    H->a[2 + k] = EB[2 * k    ];
    H->b[1 + k] = EB[2 * k + 1];
  }
  free(EB);
  return 0;
}

/*  The natural impurity orbital form adds its two chains where the double
 *  tri-diagonal form nests them (see the struct comments), so there is no
 *  rotation from one to the other -- the two split the same G differently. What
 *  they agree on are the poles of G, so that is the route, and it is the same one
 *  BADoubleTriDiagonalMatrixToDoubleTriDiagonalMatrix takes the other way.
 *
 *  Going through the poles also means d -> n splits at mu exactly as l -> n does,
 *  with ListOfPolesCopyVal and ListOfPolesCopyCon, so the two agree by
 *  construction.
 */
uint32_t DoubleTriDiagonalMatrixToBADoubleTriDiagonalMatrix(DoubleTriDiagonalMatrixType const M, BADoubleTriDiagonalMatrixType * const H, double const mu, double const kBT)
{
  ListOfPoles LP;
  if (DoubleTriDiagonalMatrixToListOfPoles(M, &LP))
  {
    printf("DoubleTriDiagonalMatrixToListOfPoles failed in DoubleTriDiagonalMatrixToBADoubleTriDiagonalMatrix\n");
    fflush(stdout);
    return 1;
  }
  if (ListOfPolesToBADoubleTridiagonalMatrix(LP, H, LP.N, LP.N, mu, kBT))
  {
    printf("ListOfPolesToBADoubleTridiagonalMatrix failed in DoubleTriDiagonalMatrixToBADoubleTriDiagonalMatrix\n");
    fflush(stdout);
    FreeListOfPoles(&LP);
    return 1;
  }
  FreeListOfPoles(&LP);
  strncpy(H->Name, M.Name, MAXStrLength);
  return 0;
}

/*  The whole matrix is tri-diagonal, with the impurity in the middle
 *
 *  The poles of G are the eigenvalues of the whole matrix, with weights the
 *  square of what the impurity sees of each eigenvector -- not the two chain
 *  spectra, which are the poles of the hybridization and are what
 *  DoubleTriDiagonalMatrixToAndersonMatrix wants.
 *
 *  So this needs the whole matrix diagonalized, and it is cheap: read it as
 *
 *      val_N ... val_2, val_1, impurity, con_1, con_2, ... con_N
 *
 *  and it is an ordinary tri-diagonal matrix, half bandwidth 1. Only the observed
 *  site sits in the middle rather than at an end, which costs nothing here --
 *  the weights are read from the middle row of the eigenvectors instead of the
 *  first. One tri-diagonal eigenvalue problem of size 1+val.N+con.N and no more.
 *
 *  (The interleaved ordering that DoubleTriDiagonalMatrixToTriDiagonalMatrix uses
 *  is only needed when the observed site has to come first, which a Lanczos chain
 *  requires and an eigenvalue decomposition does not.)
 *
 *  Poles of zero weight are kept, as everywhere else in the d row; ReducePoles
 *  and Chop remove them. dstevd returns the eigenvalues ascending, so the list
 *  comes out ordered.
 */
uint32_t DoubleTriDiagonalMatrixToListOfPoles(DoubleTriDiagonalMatrixType const M, ListOfPoles * const H)
{
  TriDiagonalMatrixType T;
  CompactMatrixType U;
  uint32_t N, k, imp;
  N = 1 + M.val.N + M.con.N;
  imp = M.val.N;                 /* 0 based index of the impurity, site M.val.N+1 */
  T.N = N;
  if (InitTriDiagonalMatrix(&T))
  {
    printf("InitTriDiagonalMatrix failed in DoubleTriDiagonalMatrixToListOfPoles\n");
    fflush(stdout);
    return 1;
  }
  strncpy(T.Name, M.Name, MAXStrLength);
  /*  the valence chain backwards, then the impurity, then the conduction chain.
   *  b[i] is the bond between site i and site i+1 */
  for (k = 1; k <= M.val.N; k++)
  {
    T.a[M.val.N + 1 - k] = M.val.a[k];
  }
  T.a[M.val.N + 1] = M.val.a[0] + M.con.a[0];
  for (k = 1; k <= M.con.N; k++)
  {
    T.a[M.val.N + 1 + k] = M.con.a[k];
  }
  for (k = 1; k < M.val.N; k++)
  {
    T.b[M.val.N - k] = M.val.b[k];
  }
  if (M.val.N != 0) T.b[M.val.N    ] = M.val.b[0];
  if (M.con.N != 0) T.b[M.val.N + 1] = M.con.b[0];
  for (k = 1; k < M.con.N; k++)
  {
    T.b[M.val.N + 1 + k] = M.con.b[k];
  }
  H->N = N;
  if (InitListOfPoles(H))
  {
    printf("InitListOfPoles failed in DoubleTriDiagonalMatrixToListOfPoles\n");
    fflush(stdout);
    FreeTriDiagonalMatrix(&T);
    return 1;
  }
  if (TriDiagonalMatrixToEigensystemCompact(&T, H->a, &U))
  {
    printf("TriDiagonalMatrixToEigensystemCompact failed in DoubleTriDiagonalMatrixToListOfPoles\n");
    fflush(stdout);
    FreeTriDiagonalMatrix(&T);
    return 1;
  }
  FreeTriDiagonalMatrix(&T);
  /*  U.VRe[k*N+j] is component j of eigenvector k */
  for (k = 0; k < N; k++)
  {
    H->b[k] = Sqr(M.b0 * U.VRe[k * N + imp]);
  }
  FreeCompactMatrix(&U);
  H->a0 = M.a0;
  H->mu = M.mu;
  strncpy(H->Name, M.Name, MAXStrLength);
  return 0;
}
/*  Create the matrix whose resolvent is the Green's function represented by M: one impurity
 *  site with a valence and a conduction chain attached to it.
 *
 *  The dimension is 1 + val.N + con.N. The two chains are the matrices that
 *  TriDiagonalMatrixToCompactMatrix returns for val and con, the impurity site energy is
 *  val.a[0] + con.a[0] (both chains share this site, see StructDef.cpp) and the hopping from
 *  the impurity into the first site of a chain is val.b[0] resp. con.b[0]. Those four numbers
 *  are the a0 and b0 that TriDiagonalMatrixToCompactMatrix leaves in Mval and Mcon.
 *
 *          ( aval0+acon0 | bval0 | 0     | bcon0 | 0     )
 *      H = ( bval0       |       Mval    | 0     | 0     )
 *          ( bcon0       | 0     | 0     |     Mcon      )
 *
 *  H.a0 and H.b0 are the a0 and b0 of M itself, so that
 *    G(w) = H.a0 + H.b0^2 * [ 1/(w - H) ]_{0,0}
 *
 *  H not allocated on entry.
 */
uint32_t DoubleTriDiagonalMatrixToCompactMatrix(DoubleTriDiagonalMatrixType const M, CompactMatrixType * const H)
{
  uint32_t i,j,offset;
  CompactMatrixType Mval, Mcon;
  if (TriDiagonalMatrixToCompactMatrix(&(M.val), &Mval))
  {
    printf("TriDiagonalMatrixToCompactMatrix failed for val in DoubleTriDiagonalMatrixToCompactMatrix\n");
    return 1;
  }
  if (TriDiagonalMatrixToCompactMatrix(&(M.con), &Mcon))
  {
    printf("TriDiagonalMatrixToCompactMatrix failed for con in DoubleTriDiagonalMatrixToCompactMatrix\n");
    FreeCompactMatrix(&Mval);
    return 1;
  }
  H->QComplex = 0;
  H->Ni = 1 + Mval.Ni + Mcon.Ni;
  H->Nj = H->Ni;
  if (InitCompactMatrix(H))
  {
    printf("InitCompactMatrix failed in DoubleTriDiagonalMatrixToCompactMatrix\n");
    FreeCompactMatrix(&Mval);
    FreeCompactMatrix(&Mcon);
    return 1;
  }
  strncpy(H->Name, M.Name, MAXStrLength);
  H->Name[MAXStrLength-1]='\0';
  H->a0 = M.a0;
  H->b0 = M.b0;
  H->mu = M.mu;
  /* the impurity site, shared by the two chains */
  H->VRe[0] = Mval.a0 + Mcon.a0;
  /* the valence chain and its hopping to the impurity */
  if (Mval.Ni != 0)
  {
    H->VRe[0*H->Nj + 1] = Mval.b0;
    H->VRe[1*H->Nj + 0] = Mval.b0;
    for (i = 0; i < Mval.Ni; i++)
    {
      for (j = 0; j < Mval.Nj; j++)
      {
        H->VRe[(1+i)*H->Nj + 1+j] = Mval.VRe[i*Mval.Nj + j];
      }
    }
  }
  /* the conduction chain and its hopping to the impurity */
  if (Mcon.Ni != 0)
  {
    offset = 1 + Mval.Ni;
    H->VRe[0*H->Nj + offset] = Mcon.b0;
    H->VRe[offset*H->Nj + 0] = Mcon.b0;
    for (i = 0; i < Mcon.Ni; i++)
    {
      for (j = 0; j < Mcon.Nj; j++)
      {
        H->VRe[(offset+i)*H->Nj + offset+j] = Mcon.VRe[i*Mcon.Nj + j];
      }
    }
  }
  FreeCompactMatrix(&Mval);
  FreeCompactMatrix(&Mcon);
  return 0;
}

/*  A non block matrix has no BlockSizeList and thus nothing to deflate. Present so that the
 *  block and the non block branch of ResponseFunctionToCompactMatrix can be written the same
 *  way.
 */
uint32_t DoubleTriDiagonalMatrixToCompactMatrixWithDeflation(DoubleTriDiagonalMatrixType const M, CompactMatrixType * const H)
{
  return DoubleTriDiagonalMatrixToCompactMatrix(M, H);
}

/********************************************************************************/
/********************         Matrix conversions             ********************/
/********************         BADoubleTriDiagonalMatrix To ... ******************/
/********************************************************************************/

uint32_t BADoubleTriDiagonalMatrixTypeToCompactMatrix(const BADoubleTriDiagonalMatrixType * const M, CompactMatrixType * const A)
{
  if((M->val.N==0)&&(M->con.N==0))
  {
    A->Ni=0;
    A->Nj=0;
    A->QComplex=0;
    InitCompactMatrix(A);
    A->a0 = M->a0;
    A->b0 = 0;
    A->mu = M->mu;
    return 0;
  }
  if(M->val.N==0)
  {
    TriDiagonalMatrixToCompactMatrix(&(M->con), A);
    A->a0 = M->a0;
    A->b0*= M->b0;
    A->mu = M->mu;
    return 0;
  }
  if(M->con.N==0)
  {
    TriDiagonalMatrixToCompactMatrix(&(M->val), A);
    A->a0 = M->a0;
    A->b0*= M->b0;
    A->mu = M->mu;
    return 0;
  }
  uint32_t i,j;
  A->QComplex = 0;
  A->Ni = M->val.N+M->con.N;
  A->Nj = M->val.N+M->con.N;
  if (InitCompactMatrix(A))
  {
    printf("InitMatrix failed in TriDiagonalMatrixToMatrix\n");
    return 1;
  }
  A->a0 = M->a0;
  A->b0 = M->b0;
  A->mu = M->mu;
  A->VRe[0*A->Nj + 0] = M->acvB;
  A->VRe[0*A->Nj + 1] = M->acvBA;
  A->VRe[1*A->Nj + 0] = M->acvBA;
  A->VRe[1*A->Nj + 1] = M->acvA;
  
  if(M->val.N>1)
  {
    A->VRe[0*A->Nj + 2] = M->bvB0;
    A->VRe[1*A->Nj + 2] = M->bvA0;
    A->VRe[2*A->Nj + 0] = M->bvB0;
    A->VRe[2*A->Nj + 1] = M->bvA0;
  }
  if(M->con.N>1)
  {
    A->VRe[0*A->Nj + M->val.N+1] = M->bcB0;
    A->VRe[1*A->Nj + M->val.N+1] = M->bcA0;
    A->VRe[(M->val.N+1)*A->Nj + 0] = M->bcB0;
    A->VRe[(M->val.N+1)*A->Nj + 1] = M->bcA0;
  }
  for (i = 2; i < M->val.N; i++)
  {
    j=i+1;
    A->VRe[(j - 1)*A->Nj + j - 1] = M->val.a[i];
    A->VRe[(j    )*A->Nj + j - 1] = M->val.b[i];
    A->VRe[(j - 1)*A->Nj + j    ] = M->val.b[i];
  }
  if(M->val.N>1) A->VRe[(M->val.N  ) * A->Nj + M->val.N ] = M->val.a[M->val.N];
  
  for (i = 2; i < M->con.N; i++)
  {
    j=i+M->val.N;
    A->VRe[(j - 1)*A->Nj + j - 1] = M->con.a[i];
    A->VRe[(j    )*A->Nj + j - 1] = M->con.b[i];
    A->VRe[(j - 1)*A->Nj + j    ] = M->con.b[i];
  }
  if(M->con.N>1) A->VRe[(M->con.N +M->val.N - 1) * A->Nj + M->con.N + M->val.N -1] = M->con.a[M->con.N];
  return 0;
}

/*  On input M defined
 *  T not initiallized
 *  val allocated
 */
uint32_t BADoubleTriDiagonalMatrixTypeToEigenSystemCompact(const BADoubleTriDiagonalMatrixType * const M, double * val, CompactMatrixType * const T)
{
  CompactMatrixType Mc;
  BADoubleTriDiagonalMatrixTypeToCompactMatrix(M, &Mc);
  CompactMatrixEigensystem(&Mc, T, val);
  FreeCompactMatrix(&Mc);
  return 0;
}

uint32_t BADoubleTriDiagonalMatrixToTriDiagonalMatrix(BADoubleTriDiagonalMatrixType * const M, TriDiagonalMatrixType * const H)
{
  CompactMatrixType Mc, *T=NULL;
  BADoubleTriDiagonalMatrixTypeToCompactMatrix(M, &Mc);
  CompactMatrixTriDiagonalize(&Mc, T, H);
  FreeCompactMatrix(&Mc);
  return 0;
}

/*  The natural impurity orbital form adds its two chains, the double tri-diagonal
 *  form nests them, so there is no rotation from one to the other. The route is
 *  through the Anderson form, whose bath is one tri-diagonal eigenvalue problem
 *  away (see BADoubleTriDiagonalMatrixToAndersonMatrix) and is already the list
 *  of levels and hoppings that the split at mu needs.
 */
uint32_t BADoubleTriDiagonalMatrixToDoubleTriDiagonalMatrix(BADoubleTriDiagonalMatrixType * const M, DoubleTriDiagonalMatrixType * const H, double const mu, double const kBT)
{
  AndersonMatrixType An;
  if (BADoubleTriDiagonalMatrixToAndersonMatrix(M, &An))
  {
    printf("BADoubleTriDiagonalMatrixToAndersonMatrix failed in BADoubleTriDiagonalMatrixToDoubleTriDiagonalMatrix\n");
    fflush(stdout);
    return 1;
  }
  if (AndersonMatrixToDoubleTriDiagonalMatrix(An, H, mu, kBT))
  {
    printf("AndersonMatrixToDoubleTriDiagonalMatrix failed in BADoubleTriDiagonalMatrixToDoubleTriDiagonalMatrix\n");
    fflush(stdout);
    FreeAndersonMatrix(&An);
    return 1;
  }
  FreeAndersonMatrix(&An);
  return 0;
}

uint32_t BADoubleTriDiagonalMatrixToListOfPoles(BADoubleTriDiagonalMatrixType * const M, ListOfPoles * const H)
{
  ListOfPoles Hcon;
  TriDiagonalMatrixToListOfPoles(M->val, H);
  TriDiagonalMatrixToListOfPoles(M->con, &Hcon);
  AddListOfPoles(H, Hcon);
  FreeListOfPoles(&Hcon);
  ScaleListOfPoles(H, M->b0*M->b0);
  H->a0 = M->a0;
  return 0;
}

/*  The bath of an n object, seen from the observed orbital, is tri-diagonal
 *
 *  In the B/A basis the matrix holds the observed orbital B, the extra natural
 *  orbital A, and the two chains from their second site onwards. B couples to A
 *  and to both chain heads, and so does A:
 *
 *      av_N ... av2 - B - ac2 ... ac_N            av2 - A - ac2
 *                     |                                 |
 *                     A                                 B
 *
 *  Drop B, which is what the Anderson form asks for, and what is left is a
 *  single chain read from one end to the other,
 *
 *      av_N - ... - av3 - av2 - A - ac2 - ac3 - ... - ac_N
 *
 *  because A is left with exactly two neighbours, the head of each chain. So the
 *  bath of the Anderson form is one tri-diagonal eigenvalue problem, and the
 *  hoppings are what B sees of each eigenvector: bvB0 on av2, acvBA on A and
 *  bcB0 on ac2.
 *
 *  This used to go through the dense matrix, which costs a full symmetric
 *  eigenvalue problem and 2*N^2 of workspace; the route through the poles of G
 *  costs three tri-diagonal eigenvalue problems and a tri-diagonalization.
 *
 *  val.a[0] and con.a[0] are dropped, as they are by the dense matrix and by
 *  BADoubleTriDiagonalMatrixToListOfPoles: an additive constant on G has no home
 *  in an auxiliary matrix. Only Evaluate honours them.
 */
uint32_t BADoubleTriDiagonalMatrixToAndersonMatrix(BADoubleTriDiagonalMatrixType * const M, AndersonMatrixType * const H)
{
  TriDiagonalMatrixType Mminor;
  CompactMatrixType Tminor;
  uint32_t k, p, nv, nc, n;
  double *Eminor, V;
  /*  a chain of length zero on one side leaves an ordinary tri-diagonal matrix,
   *  the same three cases BADoubleTriDiagonalMatrixTypeToCompactMatrix has */
  if ((M->val.N == 0) && (M->con.N == 0))
  {
    H->N = 0;
    if (InitAndersonMatrix(H)) return 1;
    H->a[0] = M->a0;
    H->mu = M->mu;
    strncpy(H->Name, M->Name, MAXStrLength);
    return 0;
  }
  if ((M->val.N == 0) || (M->con.N == 0))
  {
    if (TriDiagonalMatrixToAndersonMatrix((M->val.N == 0) ? M->con : M->val, H))
    {
      printf("TriDiagonalMatrixToAndersonMatrix failed in BADoubleTriDiagonalMatrixToAndersonMatrix\n");
      fflush(stdout);
      return 1;
    }
    H->a[0] = M->a0;
    H->b[0] *= M->b0;
    H->mu = M->mu;
    strncpy(H->Name, M->Name, MAXStrLength);
    return 0;
  }
  /*  a chain of length one contributes nothing beyond the rotation that made B
   *  and A, so it has no site of its own here and its hopping is not read */
  nv = (M->val.N > 1) ? M->val.N - 1 : 0;
  nc = (M->con.N > 1) ? M->con.N - 1 : 0;
  n = 1 + nv + nc;
  Mminor.N = n;
  if (InitTriDiagonalMatrix(&Mminor))
  {
    printf("InitTriDiagonalMatrix failed in BADoubleTriDiagonalMatrixToAndersonMatrix\n");
    fflush(stdout);
    return 1;
  }
  Mminor.a[0] = 0;
  Mminor.b[0] = 1;
  /*  the valence chain backwards, then A, then the conduction chain forwards */
  for (p = 0; p < nv; p++)
  {
    Mminor.a[1 + p] = M->val.a[M->val.N - p];
  }
  Mminor.a[1 + nv] = M->acvA;
  for (p = 0; p < nc; p++)
  {
    Mminor.a[2 + nv + p] = M->con.a[2 + p];
  }
  for (p = 0; p + 1 < nv; p++)
  {
    Mminor.b[1 + p] = M->val.b[M->val.N - 1 - p];
  }
  if (nv > 0) Mminor.b[nv] = M->bvA0;
  if (nc > 0) Mminor.b[1 + nv] = M->bcA0;
  for (p = 0; p + 1 < nc; p++)
  {
    Mminor.b[2 + nv + p] = M->con.b[2 + p];
  }
  if ((Eminor = (double*)malloc(n * sizeof(double))) == NULL)
  {
    printf("malloc for Eminor failed in BADoubleTriDiagonalMatrixToAndersonMatrix\n");
    fflush(stdout);
    FreeTriDiagonalMatrix(&Mminor);
    return 1;
  }
  if (TriDiagonalMatrixToEigensystemCompact(&Mminor, Eminor, &Tminor))
  {
    printf("TriDiagonalMatrixToEigensystemCompact failed in BADoubleTriDiagonalMatrixToAndersonMatrix\n");
    fflush(stdout);
    FreeTriDiagonalMatrix(&Mminor);
    free(Eminor);
    return 1;
  }
  FreeTriDiagonalMatrix(&Mminor);
  H->N = 1 + n;
  if (InitAndersonMatrix(H))
  {
    printf("InitAndersonMatrix failed in BADoubleTriDiagonalMatrixToAndersonMatrix\n");
    fflush(stdout);
    FreeCompactMatrix(&Tminor);
    free(Eminor);
    return 1;
  }
  H->a[0] = M->a0;
  H->b[0] = M->b0;
  H->a[1] = M->acvB;
  H->mu = M->mu;
  strncpy(H->Name, M->Name, MAXStrLength);
  /*  Tminor.VRe[k*n+j] is component j of eigenvector k. The sign of a hopping is
   *  a gauge on the bath orbital, so only the absolute value is kept */
  for (k = 0; k < n; k++)
  {
    H->a[2 + k] = Eminor[k];
    V = M->acvBA * Tminor.VRe[k * n + nv];
    if (nv > 0) V += M->bvB0 * Tminor.VRe[k * n + nv - 1];
    if (nc > 0) V += M->bcB0 * Tminor.VRe[k * n + nv + 1];
    H->b[1 + k] = fabs(V);
  }
  FreeCompactMatrix(&Tminor);
  free(Eminor);
  return 0;
}


/**************************************************************************************/
/********************         Matrix conversions                   ********************/
/********************         SecondDoubleTriDiagonalMatrix To ... ********************/
/**************************************************************************************/

uint32_t SecondDoubleTriDiagonalMatrixToTriDiagonalMatrix(SecondDoubleTriDiagonalMatrixType const M, TriDiagonalMatrixType * const H, uint32_t const NTriDiag)
{
  uint32_t i, Ncopy;
  TriDiagonalMatrixType Hlow;
  Hlow.N = 0;
  if (NTriDiag > 1)
    if (DoubleTriDiagonalMatrixToTriDiagonalMatrix(M.Dt, &Hlow))
    {
      printf("DoubleTriDiagonalMatrixToTriDiagonalMatrix failed in SecondDoubleTriDiagonalMatrixToTriDiagonalMatrix\n");
      return 1;
    };
  /*  the inner conversion returns the whole chain now, so take what fits */
  Ncopy = Min(NTriDiag - 1, Hlow.N);
  H->N = NTriDiag;
  InitTriDiagonalMatrix(H);
  H->a[0] = M.a0;
  H->b[0] = M.b0;
  H->a[1] = Hlow.a[0];
  for (i = 0; i < Ncopy; i++)
  {
    H->b[i + 1] = Hlow.b[i];
    H->a[i + 2] = Hlow.a[i + 1];
  }
  FreeTriDiagonalMatrix(&Hlow);
  return 0;
}





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
uint32_t SparseMatrixToTriDiagonalMatrix(SparseMatrixType M, MatrixType * const T, TriDiagonalMatrixType * const H, uint32_t NTriDiag)
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

/*  DO NOT USE THIS FUNCTION
 *  use CompactMatrixTriDiagonalize(CompactMatrixType *H, CompactMatrixType *T, TriDiagonalMatrixType *HTri)
 *
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
/*  DO NOT USE THIS FUNCTION
 *  use CompactMatrixTriDiagonalize(CompactMatrixType *H, CompactMatrixType *T, TriDiagonalMatrixType *HTri)
 */
uint32_t MatrixToTriDiagonalMatrix(MatrixType const M, MatrixType * const T, TriDiagonalMatrixType * const H)
{
  return TriDiagonalize(M, T, H);
}
uint32_t MatrixToTriDiagonalMatrixNonOrthogonal(MatrixType const M, TriDiagonalMatrixType * const H)
{
  /* to do */
  return 1;
}




/*  This function is outdated use the functions acting on compact matrix instead
 *  you find them in the file BasicMatch_CompactMatrixVector
 *
*  M=T^{*}.H.T^{T} with H a diagonal matrix with the values val on the diagonal
*  input M defined, T not initiallized
*
*/
uint32_t Eigensystem(MatrixType const M, MatrixType * const T, double *val)
{
  int    info;
  const char *jobz = "V";
  const char *uplo = "U";
  uint32_t    lwork = 1 + 6 * M.Ni + 2 * Sqr(M.Ni);
  uint32_t    liwork = 3 + 5 * M.Ni;
  double *work;
  int   *iwork;
  work = (double*)malloc(lwork * sizeof(double));
  if (work == NULL)
  {
    printf("malloc of work failed in Eigensystem\n");
    return 1;
  };
  iwork = (int*)malloc(liwork * sizeof(int));
  if (iwork == NULL)
  {
    printf("malloc of iwork failed in Eigensystem\n");
    free(work);
    return 1;
  };
  if (M.Ni != M.Nj)
  {
    printf("Eigensystem called with non-square matrix\n");
    printf("ERORR!!!!!!\n");
    free(work);
    free(iwork);
    return 1;
  }
  if (M.QComplex == 0)
  {                             /* Real */
    if (MatrixCopy(T, M)) return 1;
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
      {
        free(work);
        free(iwork);
        return 1;
      }
      strcpy(T->Name, "Eigenvectors of ");
      strncat(T->Name, M.Name, Max(0, MAXStrLength - 16));
      dsyevd(jobz, uplo, &M.Ni, T->VRe[0], &M.Ni, val, work, &lwork, iwork, &liwork, &info);
      if (info)
      {
        printf("also dsyvd failed in Eigensystem with info value %i\n", info);
        printf("ERORR!!!!!!!!!!\n");
        fflush(stdout);
        free(work);
        free(iwork);
        return 1;
      }
    }
  }
  else
  {                             /* complex */
    Complex16 MReIm[M.Ni * M.Ni];
    uint32_t    lcwork = 2 * M.Ni + Sqr(M.Ni);
    Complex16 cwork[lcwork];
    int    i;
    T->Ni = M.Ni;
    T->Nj = M.Nj;
    T->QComplex = M.QComplex;
    if (InitMatrix(T))
    {
      free(work);
      free(iwork);
      return 1;
    }
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
        free(work);
        free(iwork);
        return 1;
      }
    }
    for (i = 0; i < M.Ni * M.Ni; i++)
    {
      T->VRe[0][i] = MReIm[i].Re;
      T->VIm[0][i] = -MReIm[i].Im;
    }
  }
  free(work);
  free(iwork);
  return 0;
}
/*  This function is outdated use the functions acting on compact matrix instead
 *  you find them in the file BasicMatch_CompactMatrixVector
 */
uint32_t MatrixToEigensystem(MatrixType const M, MatrixType * const T, double *val)
{
  return Eigensystem(M, T, val);
}

/*  This function is outdated use the functions acting on compact matrix instead
 *  you find them in the file BasicMatch_CompactMatrixVector
 *
 *  Returns only the eigenvalues of M, not the eigenvectors.
 *  input M defined; val must point to an array of at least M.Ni doubles.
 *
 */
uint32_t Eigenvalues(MatrixType const M, double *val)
{
  int    info;
  const char *jobz = "N";
  const char *uplo = "U";
  uint32_t    lwork = 1 + 6 * M.Ni + 2 * Sqr(M.Ni);
  uint32_t    liwork = 3 + 5 * M.Ni;
  double *work;
  int   *iwork;
  if (M.Ni != M.Nj)
  {
    printf("Eigenvalues called with non-square matrix\n");
    printf("ERORR!!!!!!\n");
    return 1;
  }
  work = (double*)malloc(lwork * sizeof(double));
  if (work == NULL)
  {
    printf("malloc of work failed in Eigenvalues\n");
    return 1;
  };
  iwork = (int*)malloc(liwork * sizeof(int));
  if (iwork == NULL)
  {
    printf("malloc of iwork failed in Eigenvalues\n");
    free(work);
    return 1;
  };
  if (M.QComplex == 0)
  {                             /* Real */
    double *MRe = (double*)malloc(M.Ni * M.Ni * sizeof(double));
    if (MRe == NULL)
    {
      printf("malloc of MRe failed in Eigenvalues\n");
      free(work);
      free(iwork);
      return 1;
    }
    memcpy(MRe, M.VRe[0], M.Ni * M.Ni * sizeof(double));
    dsyev(jobz, uplo, &M.Ni, MRe, &M.Ni, val, work, &lwork, &info);
    if (info)
    {
      printf("dsyv failed in Eigenvalues with info value %i\n", info);
      printf("While working on matrix %s\n", M.Name);
      printf("continue diagonalization with dsyvd\n");
      fflush(stdout);
      memcpy(MRe, M.VRe[0], M.Ni * M.Ni * sizeof(double));
      dsyevd(jobz, uplo, &M.Ni, MRe, &M.Ni, val, work, &lwork, iwork, &liwork, &info);
      if (info)
      {
        printf("also dsyvd failed in Eigenvalues with info value %i\n", info);
        printf("ERORR!!!!!!!!!!\n");
        fflush(stdout);
        free(MRe);
        free(work);
        free(iwork);
        return 1;
      }
    }
    free(MRe);
  }
  else
  {                             /* complex */
    uint32_t    lcwork = 2 * M.Ni + Sqr(M.Ni);
    Complex16 *MReIm = (Complex16*)malloc(M.Ni * M.Ni * sizeof(Complex16));
    Complex16 *cwork = (Complex16*)malloc(lcwork * sizeof(Complex16));
    int    i;
    if (MReIm == NULL || cwork == NULL)
    {
      printf("malloc of MReIm/cwork failed in Eigenvalues\n");
      free(MReIm);
      free(cwork);
      free(work);
      free(iwork);
      return 1;
    }
    for (i = 0; i < M.Ni * M.Ni; i++)
    {
      MReIm[i].Re = M.VRe[0][i];
      MReIm[i].Im = M.VIm[0][i];
    }
    zheev(jobz, uplo, &M.Ni, MReIm, &M.Ni, val, cwork, &lcwork, work, &info);
    if (info != 0)
    {
      printf("zheev failed in Eigenvalues with info value %i\n", info);
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
        printf("also zheevd failed in Eigenvalues with info value %i\n", info);
        printf("ERORR!!!!!!!!!!\n");
        fflush(stdout);
        free(MReIm);
        free(cwork);
        free(work);
        free(iwork);
        return 1;
      }
    }
    free(MReIm);
    free(cwork);
  }
  free(work);
  free(iwork);
  return 0;
}




/*  This function is outdated use the functions acting on compact matrix instead
 *  you find them in the file BasicMatch_CompactMatrixVector
 *  Calculates the eigen values (val) and eigen functions (T) of matrix M
 *  H=T^*.M.T^T with H a diagonal matrix with the values val on the diagonal
 *  input M defined, H and T not initiallized
 */
uint32_t EigensystemToCompact(MatrixType const M, CompactMatrixType * const T, double *val)
{
  int    info;
  const char *jobz = "V";
  const char *uplo = "U";
  uint32_t    lwork = 1 + 6 * M.Ni + 2 * Sqr(M.Ni);
  uint32_t    liwork = 3 + 5 * M.Ni;
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
    uint32_t  lcwork = 2 * M.Ni + Sqr(M.Ni);
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


/*  This function is outdated use the functions acting on compact matrix instead
 *  you find them in the file BasicMatch_CompactMatrixVector
 */
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
  int    info;
  uint32_t itype = 1;
  const char *jobz = "V";
  const char *uplo = "U";
  uint32_t   lwork = 1 + 6 * M.Ni + 2 * Sqr(M.Ni);
  uint32_t   liwork = 3 + 5 * M.Ni;
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
      dsygvd(&itype, jobz, uplo, &M.Ni, T->VRe[0], &M.Ni, Sroot->VRe[0], &M.Ni, val, work, &lwork, iwork, &liwork, &info);
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
    uint32_t  lcwork = 2 * M.Ni + Sqr(M.Ni);
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
      zhegvd(&itype, jobz, uplo, &M.Ni, MReIm, &M.Ni, SReIm, &M.Ni, val, cwork, &lcwork, work, &lwork, iwork, &liwork, &info);
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
  int    info;
  uint32_t itype = 1;
  const char *jobz = "N";
  const char *uplo = "U";
  uint32_t lwork = 1 + 6 * M.Ni + 2 * Sqr(M.Ni);
  uint32_t liwork = 3 + 5 * M.Ni;
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
    uint32_t  lcwork = 2 * M.Ni + Sqr(M.Ni);
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

/*  Rotate degenerate states of an eigensystem
 *  States are degenerate if |val[i+1]-val[i]| < epsilon
 *  Method determines how the states are rotated
 */
uint32_t RotateStates(CompactMatrixType *fun, uint32_t start, uint32_t end, char *method)
{
  if(strcmp(method, "QRDecomposition")==0)
  {
    // for N degenerate states, find the N columns with largest elements, and use a QR decomposite on this matrix
    CompactMatrixType subfun;
    subfun.Ni=end-start+1;
    subfun.Nj=end-start+1;
    if(fun->QComplex==0)
    { // real
      //TODO: implement real QR rotation
    }
    else
    { // complex
      //TODO: implement complex QR rotation
    }
  }
  else
  {
    printf("ERROR: Unknown method in Rotate: %s",method);
    return 1;
  }
  return 0;
}
uint32_t RotateDegenerateStates(CompactMatrixType *fun, double *val, double epsilon, char *method)
{
  uint32_t degenerateStart=0, degenerateEnd=0;
  while(degenerateEnd<(fun->Ni-1))
  {
    if(Abs(val[degenerateStart]-val[degenerateEnd+1])<=epsilon)
    {
      degenerateEnd++;
    }
    else
    {
      if(degenerateStart==degenerateEnd)
      {
        degenerateStart++;
        degenerateEnd++;
        continue;
      }
      else
      {
        RotateStates(fun, degenerateStart, degenerateEnd, method);
        degenerateEnd++;
        degenerateStart=degenerateEnd;
      }
    }
    if((degenerateEnd==(fun->Ni-1))&&(degenerateEnd>degenerateStart))
    {
      RotateStates(fun, degenerateStart, degenerateEnd, method);
      degenerateEnd++;
      degenerateStart=degenerateEnd;
    }
  }
  return 0;
}

/*  turns a matrix into a compact matrix. The matrix A is destroyed (freed) on access
 *
 */
uint32_t MatrixToCompactMatrix(MatrixType *A, CompactMatrixType *B)
{
  uint32_t i;
  B->Ni=A->Ni;
  B->Nj=A->Nj;
  B->QComplex=A->QComplex;
  if(A->QComplex==0)
  {
    B->VRe=A->VRe[0];
    if((A->Ni!=0)&&(A->Nj!=0)) free(A->VRe);
  }
  else
  {
    InitCompactMatrix(B);
    for(i=0;i<B->Ni*B->Nj;i++)
    {
      B->V[i].Re=A->VRe[0][i];
      B->V[i].Im=A->VIm[0][i];
    }
    if((A->Ni!=0)&&(A->Nj!=0))
    {
      free(A->VRe[0]);
      free(A->VRe);
      free(A->VIm[0]);
      free(A->VIm);
    }
  }
  B->a0=A->a0;
  B->b0=A->b0;
  B->mu=A->mu;
  B->OnLuaStack=0;
  B->Name=(char*)malloc(MAXStrLength*sizeof(char));
  strcpy(B->Name,A->Name);
  return 0;
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


/*  DO NOT USE THIS FUNCTION, USE CompactMatrixToAndersonMatrix instead
 *
*  M=T.H.T^{\dag}
*  input M defined, H and T not initiallized
*/
uint32_t MatrixToAndersonMatrix(MatrixType const M, MatrixType * const T, AndersonMatrixType * const H)
{
  //TODO: M=T.H.T^dag does not seem right, should be H=T^*.M.T^T
  int    info;
  uint32_t i,j;
  const char *jobz = "V";
  const char *uplo = "U";
  const uint32_t lwork = 1 + 6 * M.Ni + 2 * Sqr(M.Ni);
  const uint32_t liwork = 3 + 5 * M.Ni;
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
    uint32_t  lcwork = 2 * M.Ni + Sqr(M.Ni);
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
*  M=T^*.H.T^T
*  input M defined, H and T not initiallized
*/
uint32_t CompactMatrixToAndersonMatrix(CompactMatrixType * const M, CompactMatrixType * const T, AndersonMatrixType * const H)
{
  if(M->Ni==0)
  {
    T->Ni=0;
    T->Nj=0;
    T->QComplex=0;
    InitCompactMatrix(T);
    H->N=0;
    InitAndersonMatrix(H);
    H->a[0]=M->a0;
    H->mu=M->mu;
    return 0;
  }
  int    info;
  uint32_t i,j;
  const char jobz  = 'V';
  const char trans = 'T'; // Difference between C and Fortran Matrix_C . Vector_C = Matrix_{Fortran}^T . Vector_{Fortran}
  const char uplo  = 'U';
  /*  dsyev and zheev are QR iterations and want a workspace of O(N): (NB+2)*N
   *  doubles and (NB+1)*N Complex16 for the block size NB the library picks,
   *  plus 3N-2 doubles of rwork for zheev. What used to stand here was
   *  1+6N+2N^2 doubles and 2N+N^2 Complex16, which is what dsyevd and zheevd,
   *  the divide and conquer solvers, want -- and they were stack arrays. That is
   *  30 Mb at N=1000 and 3 Gb at N=10000, so this function died with SIGSEGV
   *  above about N=500.
   *
   *  NB is 32 in Accelerate; 128 is taken here as an upper bound over libraries.
   *  Being generous costs nothing, LAPACK uses no more than the blocked size, so
   *  this runs at exactly the speed the oversized workspace did.
   */
  const uint32_t NBmax = 128;
  uint32_t lwork, lrwork;
  double *work=NULL, *rwork=NULL;
  double pRe, pIm, vRe, vIm;
  Complex16 *cwork=NULL;
  Complex16  alpha;
  Complex16  beta;
  alpha.Re=1; alpha.Im=0;
  beta.Re=0;  beta.Im=0;
  if (M->Ni != M->Nj)
  {
    printf("MatrixToAndersonMatrix called with non-square matrix\n");
    printf("ERORR!!!!!!\n");
    return 1;
  }
  H->N = M->Ni;
  if (InitAndersonMatrix(H)) return 1;
  CompactMatrixCopy(*M, T);
  H->a[0] = T->a0;
  H->b[0] = T->b0;
  H->mu   = T->mu;
  i=(M->Ni)-1;
  j=1;
  if (T->QComplex == 0)
  {
    H->a[1] = M->VRe[0];
    lwork = (NBmax + 2) * i + 1;                  /* the documented minimum is 3N-1 */
    work = (double*)malloc(lwork * sizeof(double));
    if (work == NULL)
    {
      printf("malloc of work failed in CompactMatrixToAndersonMatrix\n");
      fflush(stdout);
      return 1;
    }
    dsyev(&jobz, &uplo, &i, &(T->VRe[1+M->Ni]), &M->Ni, &(H->a[2]), work, &lwork, &info);
    if (info != 0) printf("dsyev failed in CompactMatrixToAndersonMatrix with info value %i\n", info);
    free(work);
    T->VRe[0]=1;
    dgemv(&trans, &i, &i, &(alpha.Re), &(T->VRe[1+M->Ni]), &M->Ni, &(M->VRe[1]), &j, &(beta.Re), &(H->b[1]), &j);
    for(i=1;i<M->Ni;i++)
    {
      T->VRe[i]=0;
      T->VRe[i*M->Ni]=0;
      if(H->b[i]<0)
      {
        H->b[i]*=-1;
        for(j=1;j<M->Ni;j++)
        {
          T->VRe[i*T->Nj+j]*=-1;
        }
      }
    }
  }
  else
  {
    H->a[1] = M->V[0].Re;
    /*  zheev keeps its real workspace apart from its complex one, and wants
     *  only 3N-2 of it */
    lrwork = 3 * i + 1;
    rwork = (double*)malloc(lrwork * sizeof(double));
    if (rwork == NULL)
    {
      printf("malloc of rwork failed in CompactMatrixToAndersonMatrix\n");
      fflush(stdout);
      return 1;
    }
    /*  the documented minimum is 2N-1, and zgemv below writes i elements into
     *  cwork, so (NBmax+1)*i+1 covers both */
    lwork = (NBmax + 1) * i + 1;
    cwork = (Complex16*)malloc(lwork * sizeof(Complex16));
    if (cwork == NULL)
    {
      printf("malloc of cwork failed in CompactMatrixToAndersonMatrix\n");
      fflush(stdout);
      free(rwork);
      return 1;
    }
    zheev(&jobz, &uplo, &i, &(T->V[1+M->Ni]),   &M->Ni, &(H->a[2]), cwork, &lwork, rwork, &info);
    if (info != 0) printf("zheev failed in CompactMatrixToAndersonMatrix with info value %i\n", info);
    free(rwork);
    T->V[0].Re=1;
    T->V[0].Im=0;
    /*  the impurity couples to the bath through the column below it, M[r][0].
     *  For a real symmetric matrix that is the same as the row M[0][r] the real
     *  branch reads, but for a Hermitian one it is its complex conjugate, and
     *  the conjugate is the wrong one here: a row major C array read column
     *  major by fortran is conj(M), so zheev returns the eigenvectors of conj(M)
     *  and the coupling U^dagger w comes out as a plain, unconjugated sum over
     *  them. Hence the column, with M->Ni as the stride, and not the row.
     */
    zgemv(&trans, &i, &i, &alpha, &(T->V[1+M->Ni]), &M->Ni, &(M->V[M->Ni]), &M->Ni, &beta, cwork, &j);
    /*  cwork[i-1] is the coupling of the impurity to eigenvector i-1, which sits
     *  in row i of T. An AndersonMatrixType keeps its hoppings real, and it may:
     *  the phase of a coupling is a gauge on the bath orbital it points at, so
     *  rotating that orbital by the conjugate phase makes the hopping real and
     *  positive and leaves the matrix Hermitian. This is the complex version of
     *  the sign flip the real branch does just above.
     */
    for(i=1;i<M->Ni;i++)
    {
      T->V[i].Re=0;
      T->V[i*M->Ni].Re=0;
      T->V[i].Im=0;
      T->V[i*M->Ni].Im=0;
      H->b[i] = sqrt(Sqr(cwork[i-1].Re) + Sqr(cwork[i-1].Im));
      if(H->b[i]>0)
      {
        pRe =  cwork[i-1].Re / H->b[i];
        pIm = -cwork[i-1].Im / H->b[i];
        for(j=1;j<M->Ni;j++)
        {
          vRe = T->V[i*T->Nj+j].Re;
          vIm = T->V[i*T->Nj+j].Im;
          T->V[i*T->Nj+j].Re = vRe*pRe - vIm*pIm;
          T->V[i*T->Nj+j].Im = vRe*pIm + vIm*pRe;
        }
      }
    }
    free(cwork);
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
  //M=T.H.T^Dag does not seem right, should be H=T^*.M.T^T
  int    i, j, k, l, CurrentDim, NBlok;
  const char *jobu = "N";
  const char *jobvt = "A";
  const char *jobz = "V";
  const char *uplo = "U";
  double SingularVal[M.Ni], val[M.Ni];
  double *U;
  const uint32_t lwork = 1 + 6 * M.Ni + 2 * Sqr(M.Ni);
  double work[lwork];
  const uint32_t liwork = 3 + 5 * M.Ni;
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
      U = (double*)malloc(THT.Ni * THT.Ni * sizeof(double));
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
    const uint32_t Clwork = 1 + 2 * M.Ni + Sqr(M.Ni);
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
      UC = (Complex16*)malloc(THT.Ni * THT.Ni * sizeof(Complex16));
      if (UC == NULL)
      {
        printf("Allocation of UC failed in BlokTriDiagonalize C 20\n");
        fflush(stdout);
        FreeMatrix(T);
        FreeMatrix(&THT);
        FreeBlokTriDiagonalMatrix(H);
        return 1;
      }
      THTC = (Complex16*)malloc(THT.Ni * THT.Nj * sizeof(Complex16));
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
      zgesvd(jobu, jobvt, &(THT.Nj), &(THT.Ni), THTC, &(THT.Nj), SingularVal, UC, &(THT.Ni), UC, &(THT.Ni), Cwork, &Clwork, work, &info);
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
        THTC = (Complex16*)malloc(THT.Ni * THT.Nj * sizeof(Complex16));
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
        THTC = (Complex16*)malloc(THT.Ni * THT.Nj * sizeof(Complex16));
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
  //M=T.H.T^Dag does not seem right, should be H=T^*.M.T^T
  uint32_t i, j, NBlok;
  int info;
  const char *jobz = "V";
  const char *uplo = "U";
  uint32_t    lwork = 1 + 6 * M.Ni + 2 * Sqr(M.Ni);
  uint32_t    liwork = 3 + 5 * M.Ni;
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
    uint32_t  lcwork = 2 * M.Ni + Sqr(M.Ni);
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


  //general coding

  return 1;
}


uint32_t MatrixToBlokDoubleTriDiagonalMatrixNoOverlap(MatrixType const M, MatrixType * const Tval, MatrixType * const Tcon, BlokDoubleTriDiagonalMatrixType * const H, const uint32_t BlokDim, double const Ef)
{
  //general coding
  return 1;
}


/****************************************************************************/
/********************         Matrix conversions         ********************/
/********************         ListOfPoles To ...         ********************/
/****************************************************************************/


/*  Given a Green's function in a list of poles representation
 *    g = a0 + sum_i b_i / (w - a_i)
 *  the occupation is
 *    rho = sum_i f(a_i, mu, kbT) b_i
 *  with f the Fermi function. A pole exactly at mu is counted half.
 *  The constant offset a0 carries no spectral weight and does not contribute to rho.
 *  rho is returned as a 1 by 1 matrix such that the scalar and the block case return the same type
 */
uint32_t ListOfPolesToDensityMatrix(ListOfPoles const M, CompactMatrixType * const rho, double const mu, double const kbT)
{
  uint32_t i;
  rho->Ni=1;
  rho->Nj=1;
  rho->QComplex=0;
  InitCompactMatrix(rho);
  rho->mu=mu;
  strncpy(rho->Name, M.Name, MAXStrLength);
  for(i=0;i<M.N;i++)
  {
    rho->VRe[0] += FermiFunction(M.a[i], mu, kbT) * M.b[i];
  }
  return 0;
}


/*
 *
 */
uint32_t ListOfPolesToCompactMatrix(const ListOfPoles * const A, CompactMatrixType * const AM)
{
  TriDiagonalMatrixType ATri;
  ListOfPolesToTridiagonalMatrix(*A, &ATri, A->N);
  TriDiagonalMatrixToCompactMatrix(&ATri, AM);
  FreeTriDiagonalMatrix(&ATri);
  return 0;
}


/*
 *
 */
uint32_t ListOfPolesToEigensystemCompact(const ListOfPoles * const A, double * val, CompactMatrixType * const T)
{
  memcpy(val, &(A->a), A->N*sizeof(double));
  T->Ni=A->N;
  T->Nj=A->N;
  T->QComplex=0;
  InitCompactIdentityMatrix(T);
  return 0;
}
/*
 *
 */
uint32_t ListOfPolesToTridiagonalMatrixNonOrthogonal(const ListOfPoles A, TriDiagonalMatrixType * const M, uint32_t NTriDiag)
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
  strncpy(AM.Name, A.Name, MAXStrLength);
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

/*  g(w) = a0 + sum_i b[i]/(w - a[i]) becomes
 *  a0 + b[0]^2/(w - a[1] - b[1]^2/(w - a[2] - ...)), with b[0]^2 the total
 *  weight of the poles.
 *
 *  NTriDiag is an upper bound on the dimension. The chain is shorter when the
 *  poles are fewer or degenerate, or when a hopping vanishes.
 */
uint32_t ListOfPolesToTridiagonalMatrix(const ListOfPoles A, TriDiagonalMatrixType * const M, uint32_t NTriDiag)
{
  uint32_t i, NChain;
  double *Wpole, Wtot;
  M->N = Min(NTriDiag, A.N);
  if (InitTriDiagonalMatrix(M))
  {
    printf("InitTriDiagonalMatrix failed in ListOfPolesToTridiagonalMatrix\n");
    fflush(stdout);
    return 1;
  }
  snprintf(M->Name, MAXStrLength, "%s", A.Name);
  M->mu = A.mu;
  M->a[0] = A.a0;
  if (M->N == 0)
  {
    return 0;
  }
  if ((Wpole = (double*)malloc(A.N * sizeof(double))) == NULL)
  {
    printf("malloc for Wpole failed in ListOfPolesToTridiagonalMatrix\n");
    fflush(stdout);
    return 1;
  }
  for (i = 0; i < A.N; i++)
  {
    Wpole[i] = Abs(A.b[i]);
  }
  if (PolesToTriDiagonalCoefficients(A.N, A.a, Wpole, M->N, &(M->a[1]), &(M->b[1]), &NChain, &Wtot))
  {
    printf("PolesToTriDiagonalCoefficients failed in ListOfPolesToTridiagonalMatrix\n");
    fflush(stdout);
    free(Wpole);
    return 1;
  }
  free(Wpole);
  M->b[0] = sqrt(Wtot);
  M->N = NChain;
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
  ListOfPolesToTridiagonalMatrix(LPcon, &(M->con), Min(NTriDiagcon,LPcon.N));
  ListOfPolesToTridiagonalMatrix(LPval, &(M->val), Min(NTriDiagval,LPval.N));
  if((M->val.N==0)&&(M->con.N!=0))
  {
    M->val.a = (double*)realloc(M->val.a,2*sizeof(double));
    M->val.b = (double*)realloc(M->val.b,sizeof(double));
    M->val.a[1]=0;
    M->val.b[0]=0;
    M->val.N=1;
  }
  if((M->con.N==0)&&(M->val.N!=0))
  {
    M->con.a = (double*)realloc(M->con.a,2*sizeof(double));
    M->con.b = (double*)realloc(M->con.b,sizeof(double));
    M->con.a[1]=0;
    M->con.b[0]=0;
    M->con.N=1;
  }
  M->con.a[0]=0;
  M->val.a[0]=0;
  M->a0 = A.a0;
  M->mu = mu;
  M->Epsilon = Epsilon;
  if(M->con.N > LPcon.N)
  {
    if(M->con.N>1) M->con.N--;
  }
  if(M->val.N > LPval.N)
  {
    if(M->val.N>1) M->val.N--;
  }
  snprintf(M->Name, MAXStrLength, "%s",A.Name);
  SetBADoubleTriDiagonalMatrixConVal(M);
  FreeListOfPoles(&LPval);
  FreeListOfPoles(&LPcon);
  return 0;
}

/*
 *
 */
uint32_t ListOfPolesToBADoubleTridiagonalMatrixOccupation(const ListOfPoles A, BADoubleTriDiagonalMatrixType * const M, uint32_t NTriDiagval, uint32_t NTriDiagcon, double mu, double occupation)
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

  ListOfPoleCopyChopReduceDimension(A, &Achop, 0.01);
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

/*  A list of poles is the function itself, a0 + sum_i b[i]/(w-a[i]), while the
 *  double tri-diagonal form puts its two chains in the denominator. The poles
 *  can therefore not simply be split at mu and tri-diagonalized -- that is what
 *  the natural impurity orbital form does, and it works there only because that
 *  form adds its two chains instead of nesting them. The poles first have to
 *  become the chain whose spectrum they are, which is the Anderson form.
 */
uint32_t ListOfPolesToDoubleTridiagonalMatrix(const ListOfPoles A, DoubleTriDiagonalMatrixType * const M, double const mu, double const kBT)
{
  AndersonMatrixType An;
  if (ListOfPolesToAndersonMatrix(A, &An))
  {
    printf("ListOfPolesToAndersonMatrix failed in ListOfPolesToDoubleTridiagonalMatrix\n");
    fflush(stdout);
    return 1;
  }
  if (AndersonMatrixToDoubleTriDiagonalMatrix(An, M, mu, kBT))
  {
    printf("AndersonMatrixToDoubleTriDiagonalMatrix failed in ListOfPolesToDoubleTridiagonalMatrix\n");
    fflush(stdout);
    FreeAndersonMatrix(&An);
    return 1;
  }
  FreeAndersonMatrix(&An);
  return 0;
}

/****************************************************************************/
/********************  Matrix conversions                ********************/
/********************  BADoubleTridiagonalMatrix To ...  ********************/
/****************************************************************************/






