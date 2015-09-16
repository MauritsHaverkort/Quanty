/*
 *
 */
uint32_t FindGroundStateFromG0andSigma2(DMFTSettings const O, AndersonMatrixType const G0, ListOfPoles const Sigma,
                                       double *const mu, double *const Nelectron, WaveFunctionType * const Psi,
                                       BADoubleTriDiagonalMatrixType * const GBathBARet,
                                       AndersonMatrixType * const GBathAndRet, TimeType * const Time,
                                       double *const onsiteshiftstart)
{
  uint32_t i, j, converged = 0, iiter;
  uint32_t const nitermax = 50, nitermem = 60;
  uint16_t FermionIncluded[2];
  AndersonMatrixType GBathAnd;
  BADoubleTriDiagonalMatrixType GBathBA, GBathBATrunc;
  double WaveFunctionAccuracy, PsiEnergy, PsiQuality, onsiteshift[nitermax + 1], occupation, Rho[4 + 4 * O.NBath],
  occlist[nitermax + 1], ABocclist[nitermax + 1], mulist[nitermax + 1];
  OperatorType HAnd, Sx, Sy, Sz, S2;
  uint16_t *SpinUp, *SpinDn;
  WaveFunctionAccuracy = O.zeroLancWavefunction;
  onsiteshift[0] = *onsiteshiftstart;
  mulist[0] = *mu;
  iiter = 0;
  GBathAndListAndRepresentation(G0, *mu, Sigma, &GBathAnd);
  /*
   The bath Greens function is not symmetric around mu
   if((O.SymOpt)&2)
   {
   printf("Symmetrize the bath greens function around mu=%22.15E\n",GBathAnd.mu);
   SymmetrizeAndersonMatrix(&GBathAnd);
   }
   */
  ReducePolesAndersonMatrixNumberFermi(&GBathAnd, O.NGBath, *mu);
  while (!converged)
  {
    GBathAnd.mu = mulist[iiter];
    AndersonMatrixToBADoubleTriDiagonalMatrixTemperature(GBathAnd, onsiteshift[iiter], O.Teff, &GBathBA, &occupation);
    GBathBATrunc = GBathBA;
    GBathBATrunc.con.N = Min(GBathBA.con.N, O.NBath + 1);
    GBathBATrunc.val.N = Min(GBathBA.val.N, O.NBath + 1);
    TestIfPsiHasCorrectDimensions(GBathBATrunc, Psi);
    BADoubleTriDiagonalMatrixToAndersonImpurityOperator(&HAnd, GBathBATrunc, GBathBATrunc, O.U, 0);
    if ((O.SymOpt) & 4)
    {
      printf("Restricted HAnd for the ground-state to the first bath site\n");
      HAnd.l[2].N = 8;
    }
    /*
     The bath Greens function is not symmetric around mu
     if((O.SymOpt)&2)
     {
     printf("Symmetrize the bath greens function around mu=%22.15E\n",GBathBA.mu);
     SymmetrizeBADoubleTriDiagonalMatrix(&GBathBA);
     }
     */
    /* if singlet set is true add the (U+1) * S^2 to the operator HAndGrd */
    if ((O.SymOpt) & 1)
    {
      printf("Add %22.15E * S^2 to HAnd for the calculation of the ground-state wave-function\n", 1 + O.U);
      printf("This stabilizes the ground-state by enlarging the singlet triplet energy\n");
      SpinUp = calloc(HAnd.NFermion / 2, sizeof(uint16_t));
      SpinDn = calloc(HAnd.NFermion / 2, sizeof(uint16_t));
      SpinUp[0] = 0;
      SpinUp[1] = 2;
      SpinDn[0] = 1;
      SpinDn[1] = 3;
      for (i = 4; i < 2 + HAnd.NFermion / 2; i++)
      {
        SpinUp[i - 2] = i;
        SpinDn[i - 2] = i + HAnd.NFermion / 2 - 2;
      }
      CreateSpinOperators(HAnd.NFermion, HAnd.NFermion / 2, SpinUp, SpinDn, &Sx, &Sy, &Sz, &S2);
      FreeOperator(&Sx);
      FreeOperator(&Sy);
      FreeOperator(&Sz);
      OperatorAddNormalOrder(&HAnd, S2, fabs(O.U) + 1, 0, &QDetNotOrderedFcFaBcBa);
      FreeOperator(&S2);
    }
    /* create ground-state wave-function */
    BeginTiming("Lanczos GroundState", Time);
    LanczosGroundState(HAnd, Psi, &WaveFunctionAccuracy, O.zeroLancEnergy, O.U + *mu, Psi, 0);
    EndTiming("Lanczos GroundState", Time);
    /* calculate NelectronNew */
    CalculateDensity(Rho, *Psi);
    ABocclist[iiter] = Rho[0] + Rho[1] + Rho[2] + Rho[3];
    occlist[iiter] = Rho[0] + Rho[1];
    /* determine if ground-state calculation has to be restarted */
    if (O.NmuSym & 2)
    {
      /* calculations at constant mu, but still need to make sure the bath Greens function is calculated at the same filling as the impurity greens function */
      /* this is done by optimizing NelectronNew such that the occupation of the impurity and the selected bath site is integer */
      if (Abs(ABocclist[iiter] - O.NSite) > 100 * DBL_EPSILON)
      {
        mulist[iiter + 1] = *mu;
        if (iiter == 0)
        {
          onsiteshift[iiter + 1] = onsiteshift[iiter] - (ABocclist[iiter] - O.NSite) * (Abs(G0.b[1]) + Abs(O.U / 2));
        }
        else if (iiter == nitermax)
        {
          printf("occupation of bath and impurity green function did not converge in %i itterations\n", iiter);
          printf("loop   onsiteshift            total occupation       mu                     impurity occupation\n");
          for (j = 0; j < iiter; j++)
          {
            printf("%3i %22.15E %22.15E %22.15E %22.15E\n", j, onsiteshift[j], ABocclist[j], mulist[j], occlist[j]);
          }
          converged = 1;
          /* pass the converged onsite shift back to onsiteshift start as the next starting point */
          *onsiteshiftstart = onsiteshift[iiter];
        }
        else
        {
          onsiteshift[iiter + 1] =
          InterpolateRootOfset(&(onsiteshift[Max(iiter, nitermem) - nitermem]),
                               &(ABocclist[Max(iiter, nitermem) - nitermem]), Min(iiter + 1, nitermem), O.NSite);
        }
        iiter++;
      }
      else
      {
        converged = 1;
        /* pass the converged onsite shift back to onsiteshift start as the next starting point */
        *onsiteshiftstart = onsiteshift[iiter];
      }
    }
    else
    {
      /* calculations at constant occupation, find the correct mu and onsiteshift */
      if ((Abs(ABocclist[iiter] - O.NSite) > 100 * DBL_EPSILON)
          || (Abs(occlist[iiter] - *Nelectron) > 100 * DBL_EPSILON))
      {
        if (iiter == 0)
        {
          onsiteshift[iiter + 1] = onsiteshift[iiter] - (ABocclist[iiter] - O.NSite) * (Abs(G0.b[1]) + Abs(O.U / 2));
          mulist[iiter + 1] = mulist[iiter];
        }
        else if (iiter == nitermax)
        {
          printf("occupation of bath and impurity green function did not converge in %i itterations\n", iiter);
          printf("loop   onsiteshift            total occupation       mu                     impurity occupation\n");
          for (j = 0; j < iiter; j++)
          {
            printf("%3i %22.15E %22.15E %22.15E %22.15E\n", j, onsiteshift[j], ABocclist[j], mulist[j], occlist[j]);
          }
          converged = 1;
          *mu = mulist[nitermax];
        }
        else
        {
          InterpolateRoot2DOfset(&(onsiteshift[Max(iiter, nitermem) - nitermem]),
                                 &(mulist[Max(iiter, nitermem) - nitermem]),
                                 &(ABocclist[Max(iiter, nitermem) - nitermem]),
                                 &(occlist[Max(iiter, nitermem) - nitermem]), Min(iiter + 1, nitermem), O.NSite,
                                 *Nelectron, &(onsiteshift[iiter + 1]), &(mulist[iiter + 1]));
        }
        iiter++;
      }
      else
      {
        converged = 1;
        *mu = mulist[iiter];
      }
    }
    if (!converged)
      FreeOperator(&HAnd);
  }
  /* print info on groundstate */
  PrintInfoOnGBathAndG0InDMFTFindGroundstate(G0, GBathBA);
  PsiHermitianOperatorPsi(HAnd, *Psi, &PsiEnergy);
  PsiHermitianOperatorSqrPsi(HAnd, *Psi, &PsiQuality);
  PsiQuality /= Sqr(PsiEnergy);
  printf("Energy                  <psi|H^2|psi>/(<psi|H|psi>^2)\n");
  printf("%22.15E  %22.15E\n", PsiEnergy, PsiQuality);
  PrintWaveFunctionDensity(*Psi);
  PrintWaveFunctionPrefactorHistogram(*Psi, 2.0, 20);
  FermionIncluded[0]=0;
  FermionIncluded[1]=1;
  PrintWaveFunctionPartialCorrelationTensor(*Psi, 2, FermionIncluded);
  *GBathBARet = GBathBA;
  *GBathAndRet = GBathAnd;
  return 0;
}