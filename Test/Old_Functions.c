

uint32_t LanczosGroundStateRestrictedRR(const OperatorType H, WaveFunctionType * const psi, DeterminantRestrictionType const Restrictions, double *epsilon, const double zero, double shiftin, WaveFunctionType * OtherPsi, uint32_t NOtherPsi)
{
  printf("\nStart of LanczosGroundStateRestrictedRR\n\n");
  fflush(stdout);
  const uint32_t NKrylovMax = 100, NKrylovStart = 25;
  uint32_t NKrylov;
  WaveFunctionType v[3];
  double energy = 4 * zero, energy1 = 2 * zero, energy2;
  double shift;
  uint32_t cycles, i, vi, vip, vim, minloop = 3;
  MatrixType Eigen;
  TriDiagonalMatrixType M;
  double val[NKrylovMax];
  M.N = NKrylovMax;
  InitTriDiagonalMatrix(&M);
  // Shift H down to twice the energy of the starting wave-function
  // for the case that this is zero, add an aditional shift of 1
  // Make sure at least 2 loops are done
  PsiOperatorPsiRR(H, *psi, &energy1);
  if (energy1 < 0)
  {
    shift = -Abs(energy) - 1;
  }
  else
  {
    shift = -2 * energy - 1;
  }
  shift -= shiftin;
  energy1 = 0;
  energy = 1;
  // loop until the relative energy gain is smaller than zero.
  while (((energy - energy1) > (zero * (1.0 + fabs(energy)))) || minloop)
  {
    if (minloop)
      minloop--;
    cycles = 0;
    energy = energy1;
    // Calculate psi_new = H^3 psi_old whereby the first time H acts on psi_old the basis size is enlarged
    //  This leaves three functions allocated which can be used in the next part
  Repeat:
    if (OperatorPsiRestrictedRR(H, *psi, &v[0], Restrictions))
      goto OutOffMemv0;
    if (RealWaveFunctionAddMultipliedUnionBasis(&v[0], *psi, shift))
      goto OutOffMemv0add;
    if (RealWaveFunctionOrthogonalize(&v[0], OtherPsi, NOtherPsi))
      goto OutOffMemv0Ortho;
    RealWaveFunctionNormalize(&v[0]);
    printf("BasisSetSize After Expansion: %12i Truncation: ", v[0].N);
    RealWaveFunctionTruncateBasis(&v[0], *epsilon);
    printf("%12i\n", v[0].N);
    fflush(stdout);
    if (OperatorPsiConserveBasisRR(H, v[0], &v[1]))
      goto OutOffMemv1;
    RealWaveFunctionAddMultipliedSameBasis(v[1], v[0], shift);
    RealWaveFunctionNormalize(&v[1]);
    if (OperatorPsiConserveBasisRR(H, v[1], &v[2]))
      goto OutOffMemv2;
    RealWaveFunctionAddMultipliedSameBasis(v[2], v[1], shift);
    RealWaveFunctionNormalize(&v[2]);
    RealWaveFunctionFree(psi);
    if (OperatorPsiConserveBasisRR(H, v[2], psi))
      goto OutOffMempsi;
    RealWaveFunctionAddMultipliedSameBasis(*psi, v[2], shift);
    RealWaveFunctionOrthogonalizeConserveBasis(psi, OtherPsi, NOtherPsi);
    RealWaveFunctionNormalize(psi);
    goto EnoughMem;
  OutOffMempsi:
    RealWaveFunctionFreeValues(&v[2]);
  OutOffMemv2:
    RealWaveFunctionFreeValues(&v[1]);
  OutOffMemv1:
  OutOffMemv0Ortho:
  OutOffMemv0add:
    RealWaveFunctionFreeValues(&v[0]);
  OutOffMemv0:
    printf("Out off memory increase epsilon from %22.15E to %22.15E\n", *epsilon, *epsilon * 2);
    *epsilon *= 2.0;
    if (psi->NBlok == 1)
    {
      printf("Outof memory in LanczosGroundStateRR,\npsi->NBLok=1 so no way to recover.\n");
      AllocFailed();
    }
    printf("BasisSetSize Before Truncation: %12i ", psi->N);
    RealWaveFunctionTruncateBasis(psi, *epsilon);
    printf("After Truncation: %12i\n", psi->N);
    fflush(stdout);
    goto Repeat;
  EnoughMem:
    // For a given basis set determine the lowest state
    //  In the previous part the memory for 3 functions v[0], v[1] and v[2] is allocated on this basis set.
    energy1 = 2 * zero;
    energy2 = 4 * zero;
    // loop until the relative energy gain is smaller than zero, loop a minimum of 3 times
    while ((cycles < 3) || ((energy2 - energy1) > (zero * (1.0 + fabs(energy1)))))
    {
      // Calculate the Krylov basis of size NKrylov, starting from the sofar best approximation to the ground-state
      NKrylov = NKrylovStart + cycles % (NKrylovMax - NKrylovStart);
      M.N = NKrylov;
      cycles++;
      if ((cycles & 15) == 15)
      {
        printf("Reached %i-th cycle with energy %22.15E\n", cycles, energy1);
        fflush(stdout);
      }
      energy2 = energy1;
      RealWaveFunctionMultiply(v[0], 0);
      RealWaveFunctionAddSameBasis(v[0], *psi);
      OperatorPsiConserveBasisPreAllocatedRR(&H, &(v[0]), &(v[1]));
      RealWaveFunctionAddMultipliedSameBasis(v[1], v[0], shift);
      M.a[1] = RealWaveFunctionDotProduct(&v[1], &v[0]);
      for (i = 1; i < M.N; i++)
      {
        vip = (i + 1) % 3;
        vi = (i) % 3;
        vim = (i - 1) % 3;
        RealWaveFunctionAddMultipliedConserveBasis(&(v[vi]), v[vim], -M.a[i]);
        M.b[i] = RealWaveFunctionCalculateNorm(v[vi]);
        if (M.b[i] > DBL_MIN * 10)
          RealWaveFunctionMultiply(v[vi], (double)1.0 / M.b[i]);
        else
          RealWaveFunctionMultiply(v[vi], 0.0);
        RealWaveFunctionOrthogonalizeConserveBasis(&v[vi], OtherPsi, NOtherPsi);
        OperatorPsiConserveBasisPreAllocatedRR(&H, &(v[vi]), &(v[vip]));
        RealWaveFunctionAddMultipliedSameBasis(v[vip], v[vi], shift);
        M.a[i + 1] = RealWaveFunctionDotProduct(&v[vip], &v[vi]);
        RealWaveFunctionAddMultipliedConserveBasis(&(v[vip]), v[vim], -M.b[i]);
      }
      // Calculate the eigen-states in the Krylov basis
      TriDiagonalMatrixToEigensystem(M, &Eigen, val);
      // Given the eigen-states in the Krylov basis reconstruct the ground-state in the original basis
      RealWaveFunctionMultiply(v[0], 0);
      RealWaveFunctionAddSameBasis(v[0], *psi);
      RealWaveFunctionMultiply(*psi, Eigen.VRe[0][0]);
      OperatorPsiConserveBasisPreAllocatedRR(&H, &(v[0]), &(v[1]));
      RealWaveFunctionAddMultipliedSameBasis(v[1], v[0], shift);
      for (i = 1; i < M.N; i++)
      {
        vip = (i + 1) % 3;
        vi = (i) % 3;
        vim = (i - 1) % 3;
        RealWaveFunctionAddMultipliedSameBasis(v[vi], v[vim], -M.a[i]);
        if (M.b[i] > DBL_MIN * 10)
          RealWaveFunctionMultiply(v[vi], (double)1.0 / M.b[i]);
        else
          RealWaveFunctionMultiply(v[vi], 0.0);
        RealWaveFunctionOrthogonalizeConserveBasis(&v[vi], OtherPsi, NOtherPsi);
        OperatorPsiConserveBasisPreAllocatedRR(&H, &(v[vi]), &(v[vip]));
        RealWaveFunctionAddMultipliedSameBasis(v[vip], v[vi], shift);
        RealWaveFunctionAddMultipliedSameBasis(v[vip], v[vim], -M.b[i]);
        RealWaveFunctionAddMultipliedSameBasis(*psi, v[vi], Eigen.VRe[0][i]);
      }
      // Set the shift such that the largest eigenvalue found becomes zero
      shift -= val[M.N - 1];
      shift -= shiftin;
      // Orthogonalize to previously found functions
      RealWaveFunctionOrthogonalizeConserveBasis(psi, OtherPsi, NOtherPsi);
      RealWaveFunctionNormalize(psi);
      PsiOperatorPsiRR(H, *psi, &energy1);
    }
    RealWaveFunctionFreeValues(&v[0]);
    RealWaveFunctionFreeValues(&v[1]);
    RealWaveFunctionFreeValues(&v[2]);
    RealWaveFunctionTruncateBasis(psi, DBL_EPSILON);
    printf("                                           Eigenstate: %12i\n", psi->N);
    RealWaveFunctionNormalize(psi);
    FreeMatrix(&Eigen);
    printf("Energy: %22.15E  ", energy1);
    printf("NumberOfLaczosCycles: %i\n", cycles);
    fflush(stdout);
  }
  return 0;
}

uint32_t LanczosGroundStateRestrictedRC(const OperatorType H, WaveFunctionType * const psi, DeterminantRestrictionType const Restrictions, double *epsilon, const double zero, double shiftin, WaveFunctionType * OtherPsi, uint32_t NOtherPsi)
{
  printf("\nStart of LanczosGroundStateRestrictedRC\n\n");
  fflush(stdout);
  const uint32_t NKrylovMax = 100;
  uint32_t NKrylov;
  WaveFunctionType v[3];
  double energy = 4 * zero, energy1 = 2 * zero, energy2;
  double Imzero;
  double shift;
  uint32_t cycles, i, vi, vip, vim;
  MatrixType Eigen;
  TriDiagonalMatrixType M;
  double val[NKrylovMax];
  M.N = NKrylovMax;
  InitTriDiagonalMatrix(&M);
  while ((energy - energy1) > (zero * (1.0 + fabs(energy))))
  {
    cycles = 0;
    energy = energy1;
    shift = -shiftin - energy * 0.49;
  Repeat:
    if (OperatorPsiRestrictedRC(H, *psi, &v[0], Restrictions))
      goto OutOffMemv0;
    ComplexWaveFunctionMultiply(*psi, shift);
    if (ComplexWaveFunctionAddUnionBasis(&v[0], *psi))
      goto OutOffMemv0add;
    if (ComplexWaveFunctionOrthogonalize(&v[0], OtherPsi, NOtherPsi))
      goto OutOffMemv0Ortho;
    ComplexWaveFunctionNormalize(&v[0]);
    printf("BasisSetSize After Expansion: %12i Truncation: ", v[0].N);
    ComplexWaveFunctionTruncateBasis(&v[0], *epsilon);
    printf("%12i\n", v[0].N);
    fflush(stdout);
    if (OperatorPsiConserveBasisRC(H, v[0], &v[1]))
      goto OutOffMemv1;
    ComplexWaveFunctionMultiply(v[0], shift);
    ComplexWaveFunctionAddSameBasis(v[1], v[0]);
    ComplexWaveFunctionNormalize(&v[1]);
    if (OperatorPsiConserveBasisRC(H, v[1], &v[2]))
      goto OutOffMemv2;
    ComplexWaveFunctionMultiply(v[1], shift);
    ComplexWaveFunctionAddSameBasis(v[2], v[1]);
    ComplexWaveFunctionNormalize(&v[2]);
    ComplexWaveFunctionFree(psi);
    if (OperatorPsiConserveBasisRC(H, v[2], psi))
      goto OutOffMempsi;
    ComplexWaveFunctionMultiply(v[2], shift);
    ComplexWaveFunctionAddSameBasis(*psi, v[2]);
    ComplexWaveFunctionOrthogonalizeConserveBasis(psi, OtherPsi, NOtherPsi);
    ComplexWaveFunctionNormalize(psi);
    goto EnoughMem;
  OutOffMempsi:
    ComplexWaveFunctionFreeValues(&v[2]);
  OutOffMemv2:
    ComplexWaveFunctionFreeValues(&v[1]);
  OutOffMemv1:
  OutOffMemv0Ortho:
  OutOffMemv0add:
    ComplexWaveFunctionFreeValues(&v[0]);
    ComplexWaveFunctionMultiply(*psi, 1.0 / shift);
  OutOffMemv0:
    printf("Out off memory increase epsilon from %22.15E to %22.15E\n", *epsilon, *epsilon * 2);
    *epsilon *= 2.0;
    if (psi->NBlok == 1)
    {
      printf("Outof memory in LanczosGroundStateRC,\npsi->NBLok=1 so no way to recover.\n");
      AllocFailed();
    }
    printf("BasisSetSize Before Truncation: %12i ", psi->N);
    ComplexWaveFunctionTruncateBasis(psi, *epsilon);
    printf("After Truncation: %12i\n", psi->N);
    fflush(stdout);
    goto Repeat;
  EnoughMem:
    
    energy1 = 2 * zero;
    energy2 = 4 * zero;
    while ((cycles < 3) || ((energy2 - energy1) > (zero * (1.0 + fabs(energy1)))))
    {
      NKrylov = 5 + cycles % (NKrylovMax - 5);
      M.N = NKrylov;
      
      cycles++;
      if ((cycles & 15) == 15)
      {
        printf("Reached %i-th cycle with energy %22.15E\n", cycles, energy1);
        fflush(stdout);
      }
      energy2 = energy1;
      ComplexWaveFunctionMultiply(v[0], 0);
      ComplexWaveFunctionAddSameBasis(v[0], *psi);
      OperatorPsiConserveBasisPreAllocatedRC(&H, &(v[0]), &(v[1]));
      ComplexWaveFunctionDotProduct(&v[1], &v[0], &(M.a[1]), &Imzero);
      for (i = 1; i < M.N; i++)
      {
        vip = (i + 1) % 3;
        vi = (i) % 3;
        vim = (i - 1) % 3;
        ComplexWaveFunctionAddMultipliedConserveBasis(&(v[vi]), v[vim], -M.a[i]);
        M.b[i] = ComplexWaveFunctionCalculateNorm(v[vi]);
        if (M.b[i] > DBL_MIN * 10)
          ComplexWaveFunctionMultiply(v[vi], (double)1.0 / M.b[i]);
        else
          ComplexWaveFunctionMultiply(v[vi], 0.0);
        ComplexWaveFunctionOrthogonalizeConserveBasis(&v[vi], OtherPsi, NOtherPsi);
        OperatorPsiConserveBasisPreAllocatedRC(&H, &(v[vi]), &(v[vip]));
        ComplexWaveFunctionDotProduct(&v[vip], &v[vi], &(M.a[i + 1]), &Imzero);
        ComplexWaveFunctionAddMultipliedConserveBasis(&(v[vip]), v[vim], -M.b[i]);
      }
      TriDiagonalMatrixToEigensystem(M, &Eigen, val);
      ComplexWaveFunctionMultiply(v[0], 0);
      ComplexWaveFunctionAddSameBasis(v[0], *psi);
      ComplexWaveFunctionMultiply(*psi, Eigen.VRe[0][0]);
      OperatorPsiConserveBasisPreAllocatedRC(&H, &(v[0]), &(v[1]));
      for (i = 1; i < M.N; i++)
      {
        vip = (i + 1) % 3;
        vi = (i) % 3;
        vim = (i - 1) % 3;
        ComplexWaveFunctionAddMultipliedSameBasis(v[vi], v[vim], -M.a[i]);
        if (M.b[i] > DBL_MIN * 10)
          RealWaveFunctionMultiply(v[vi], (double)1.0 / M.b[i]);
        else
          ComplexWaveFunctionMultiply(v[vi], 0.0);
        ComplexWaveFunctionOrthogonalizeConserveBasis(&v[vi], OtherPsi, NOtherPsi);
        OperatorPsiConserveBasisPreAllocated(&H, &(v[vi]), &(v[vip]));
        ComplexWaveFunctionAddMultipliedSameBasis(v[vip], v[vim], -M.b[i]);
        ComplexWaveFunctionAddMultipliedSameBasis(*psi, v[vi], Eigen.VRe[0][i]);
      }
      ComplexWaveFunctionOrthogonalizeConserveBasis(psi, OtherPsi, NOtherPsi);
      ComplexWaveFunctionNormalize(psi);
      PsiOperatorPsiRHC(H, *psi, &energy1);
    }
    ComplexWaveFunctionFreeValues(&v[0]);
    ComplexWaveFunctionFreeValues(&v[1]);
    ComplexWaveFunctionFreeValues(&v[2]);
    FreeMatrix(&Eigen);
    printf("Energy: %22.15E  ", energy1);
    printf("NumberOfLaczosCycles: %i\n", cycles);
    fflush(stdout);
  }
  return 0;
}

uint32_t LanczosGroundStateRestrictedMC(const OperatorType H, WaveFunctionType * const psi, DeterminantRestrictionType const Restrictions, double *epsilon, const double zero, double shiftin, WaveFunctionType * OtherPsi, uint32_t NOtherPsi)
{
  printf("\nStart of LanczosGroundStateRestrictedMC\n\n");
  fflush(stdout);
  const uint32_t NKrylovMax = 100;
  uint32_t NKrylov;
  WaveFunctionType v[3];
  double energy = 4 * zero, energy1 = 2 * zero, energy2;
  double Imzero;
  double shift;
  uint32_t cycles, i, vi, vip, vim;
  MatrixType Eigen;
  TriDiagonalMatrixType M;
  double val[NKrylovMax];
  M.N = NKrylovMax;
  InitTriDiagonalMatrix(&M);
  while ((energy - energy1) > (zero * (1.0 + fabs(energy))))
  {
    cycles = 0;
    energy = energy1;
    shift = -shiftin - energy * 0.49;
  Repeat:
    if (OperatorPsiRestrictedMC(H, *psi, &v[0], Restrictions))
      goto OutOffMemv0;
    ComplexWaveFunctionMultiply(*psi, shift);
    if (ComplexWaveFunctionAddUnionBasis(&v[0], *psi))
      goto OutOffMemv0add;
    if (ComplexWaveFunctionOrthogonalize(&v[0], OtherPsi, NOtherPsi))
      goto OutOffMemv0Ortho;
    ComplexWaveFunctionNormalize(&v[0]);
    printf("BasisSetSize After Expansion: %12i Truncation: ", v[0].N);
    ComplexWaveFunctionTruncateBasis(&v[0], *epsilon);
    printf("%12i\n", v[0].N);
    fflush(stdout);
    if (OperatorPsiConserveBasisMC(H, v[0], &v[1]))
      goto OutOffMemv1;
    ComplexWaveFunctionMultiply(v[0], shift);
    ComplexWaveFunctionAddSameBasis(v[1], v[0]);
    ComplexWaveFunctionNormalize(&v[1]);
    if (OperatorPsiConserveBasisMC(H, v[1], &v[2]))
      goto OutOffMemv2;
    ComplexWaveFunctionMultiply(v[1], shift);
    ComplexWaveFunctionAddSameBasis(v[2], v[1]);
    ComplexWaveFunctionNormalize(&v[2]);
    ComplexWaveFunctionFree(psi);
    if (OperatorPsiConserveBasisMC(H, v[2], psi))
      goto OutOffMempsi;
    ComplexWaveFunctionMultiply(v[2], shift);
    ComplexWaveFunctionAddSameBasis(*psi, v[2]);
    ComplexWaveFunctionOrthogonalizeConserveBasis(psi, OtherPsi, NOtherPsi);
    ComplexWaveFunctionNormalize(psi);
    goto EnoughMem;
  OutOffMempsi:
    ComplexWaveFunctionFreeValues(&v[2]);
  OutOffMemv2:
    ComplexWaveFunctionFreeValues(&v[1]);
  OutOffMemv1:
  OutOffMemv0Ortho:
  OutOffMemv0add:
    ComplexWaveFunctionFreeValues(&v[0]);
    ComplexWaveFunctionMultiply(*psi, 1.0 / shift);
  OutOffMemv0:
    printf("Out off memory increase epsilon from %22.15E to %22.15E\n", *epsilon, *epsilon * 2);
    *epsilon *= 2.0;
    if (psi->NBlok == 1)
    {
      printf("Outof memory in LanczosGroundStateMC,\npsi->NBLok=1 so no way to recover.\n");
      AllocFailed();
    }
    printf("BasisSetSize Before Truncation: %12i ", psi->N);
    ComplexWaveFunctionTruncateBasis(psi, *epsilon);
    printf("After Truncation: %12i\n", psi->N);
    fflush(stdout);
    goto Repeat;
  EnoughMem:
    
    energy1 = 2 * zero;
    energy2 = 4 * zero;
    while ((cycles < 3) || ((energy2 - energy1) > (zero * (1.0 + fabs(energy1)))))
    {
      NKrylov = 5 + cycles % (NKrylovMax - 5);
      M.N = NKrylov;
      
      cycles++;
      if ((cycles & 15) == 15)
      {
        printf("Reached %i-th cycle with energy %22.15E\n", cycles, energy1);
        fflush(stdout);
      }
      energy2 = energy1;
      ComplexWaveFunctionMultiply(v[0], 0);
      ComplexWaveFunctionAddSameBasis(v[0], *psi);
      OperatorPsiConserveBasisPreAllocatedMC(&H, &(v[0]), &(v[1]));
      ComplexWaveFunctionDotProduct(&v[1], &v[0], &(M.a[1]), &Imzero);
      for (i = 1; i < M.N; i++)
      {
        vip = (i + 1) % 3;
        vi = (i) % 3;
        vim = (i - 1) % 3;
        ComplexWaveFunctionAddMultipliedConserveBasis(&(v[vi]), v[vim], -M.a[i]);
        M.b[i] = ComplexWaveFunctionCalculateNorm(v[vi]);
        if (M.b[i] > DBL_MIN * 10)
          ComplexWaveFunctionMultiply(v[vi], (double)1.0 / M.b[i]);
        else
          ComplexWaveFunctionMultiply(v[vi], 0.0);
        ComplexWaveFunctionOrthogonalizeConserveBasis(&v[vi], OtherPsi, NOtherPsi);
        OperatorPsiConserveBasisPreAllocatedMC(&H, &(v[vi]), &(v[vip]));
        ComplexWaveFunctionDotProduct(&v[vip], &v[vi], &(M.a[i + 1]), &Imzero);
        ComplexWaveFunctionAddMultipliedConserveBasis(&(v[vip]), v[vim], -M.b[i]);
      }
      TriDiagonalMatrixToEigensystem(M, &Eigen, val);
      ComplexWaveFunctionMultiply(v[0], 0);
      ComplexWaveFunctionAddSameBasis(v[0], *psi);
      ComplexWaveFunctionMultiply(*psi, Eigen.VRe[0][0]);
      OperatorPsiConserveBasisPreAllocatedMC(&H, &(v[0]), &(v[1]));
      for (i = 1; i < M.N; i++)
      {
        vip = (i + 1) % 3;
        vi = (i) % 3;
        vim = (i - 1) % 3;
        ComplexWaveFunctionAddMultipliedSameBasis(v[vi], v[vim], -M.a[i]);
        if (M.b[i] > DBL_MIN * 10)
          ComplexWaveFunctionMultiply(v[vi], (double)1.0 / M.b[i]);
        else
          ComplexWaveFunctionMultiply(v[vi], 0.0);
        ComplexWaveFunctionOrthogonalizeConserveBasis(&v[vi], OtherPsi, NOtherPsi);
        OperatorPsiConserveBasisPreAllocatedMC(&H, &(v[vi]), &(v[vip]));
        ComplexWaveFunctionAddMultipliedSameBasis(v[vip], v[vim], -M.b[i]);
        ComplexWaveFunctionAddMultipliedSameBasis(*psi, v[vi], Eigen.VRe[0][i]);
      }
      ComplexWaveFunctionOrthogonalizeConserveBasis(psi, OtherPsi, NOtherPsi);
      ComplexWaveFunctionNormalize(psi);
      PsiOperatorPsiMHC(H, *psi, &energy1);
    }
    ComplexWaveFunctionFreeValues(&v[0]);
    ComplexWaveFunctionFreeValues(&v[1]);
    ComplexWaveFunctionFreeValues(&v[2]);
    FreeMatrix(&Eigen);
    printf("Energy: %22.15E  ", energy1);
    printf("NumberOfLaczosCycles: %i\n", cycles);
    fflush(stdout);
  }
  return 0;
}

uint32_t LanczosGroundStateRestrictedCC(const OperatorType H, WaveFunctionType * const psi, DeterminantRestrictionType const Restrictions, double *epsilon, const double zero, double shiftin, WaveFunctionType * OtherPsi, uint32_t NOtherPsi)
{
  printf("\nStart of LanczosGroundStateRestrictedCC\n\n");
  fflush(stdout);
  const uint32_t NKrylovMax = 100;
  uint32_t NKrylov;
  WaveFunctionType v[3];
  double energy = 4 * zero, energy1 = 2 * zero, energy2;
  double Imzero;
  double shift;
  uint32_t cycles, i, vi, vip, vim;
  MatrixType Eigen;
  TriDiagonalMatrixType M;
  double val[NKrylovMax];
  M.N = NKrylovMax;
  InitTriDiagonalMatrix(&M);
  while ((energy - energy1) > (zero * (1.0 + fabs(energy))))
  {
    cycles = 0;
    energy = energy1;
    shift = -shiftin - energy * 0.49;
  Repeat:
    if (OperatorPsiRestrictedCC(H, *psi, &v[0], Restrictions))
      goto OutOffMemv0;
    ComplexWaveFunctionMultiply(*psi, shift);
    if (ComplexWaveFunctionAddUnionBasis(&v[0], *psi))
      goto OutOffMemv0add;
    if (ComplexWaveFunctionOrthogonalize(&v[0], OtherPsi, NOtherPsi))
      goto OutOffMemv0Ortho;
    ComplexWaveFunctionNormalize(&v[0]);
    printf("BasisSetSize After Expansion: %12i Truncation: ", v[0].N);
    ComplexWaveFunctionTruncateBasis(&v[0], *epsilon);
    printf("%12i\n", v[0].N);
    fflush(stdout);
    if (OperatorPsiConserveBasisCC(H, v[0], &v[1]))
      goto OutOffMemv1;
    ComplexWaveFunctionMultiply(v[0], shift);
    ComplexWaveFunctionAddSameBasis(v[1], v[0]);
    ComplexWaveFunctionNormalize(&v[1]);
    if (OperatorPsiConserveBasisCC(H, v[1], &v[2]))
      goto OutOffMemv2;
    ComplexWaveFunctionMultiply(v[1], shift);
    ComplexWaveFunctionAddSameBasis(v[2], v[1]);
    ComplexWaveFunctionNormalize(&v[2]);
    ComplexWaveFunctionFree(psi);
    if (OperatorPsiConserveBasisCC(H, v[2], psi))
      goto OutOffMempsi;
    ComplexWaveFunctionMultiply(v[2], shift);
    ComplexWaveFunctionAddSameBasis(*psi, v[2]);
    ComplexWaveFunctionOrthogonalizeConserveBasis(psi, OtherPsi, NOtherPsi);
    ComplexWaveFunctionNormalize(psi);
    goto EnoughMem;
  OutOffMempsi:
    ComplexWaveFunctionFreeValues(&v[2]);
  OutOffMemv2:
    ComplexWaveFunctionFreeValues(&v[1]);
  OutOffMemv1:
  OutOffMemv0Ortho:
  OutOffMemv0add:
    ComplexWaveFunctionFreeValues(&v[0]);
    ComplexWaveFunctionMultiply(*psi, 1.0 / shift);
  OutOffMemv0:
    printf("Out off memory increase epsilon from %22.15E to %22.15E\n", *epsilon, *epsilon * 2);
    *epsilon *= 2.0;
    if (psi->NBlok == 1)
    {
      printf("Outof memory in LanczosGroundStateCC,\npsi->NBLok=1 so no way to recover.\n");
      AllocFailed();
    }
    printf("BasisSetSize Before Truncation: %12i ", psi->N);
    ComplexWaveFunctionTruncateBasis(psi, *epsilon);
    printf("After Truncation: %12i\n", psi->N);
    fflush(stdout);
    goto Repeat;
  EnoughMem:
    
    energy1 = 2 * zero;
    energy2 = 4 * zero;
    while ((cycles < 3) || ((energy2 - energy1) > (zero * (1.0 + fabs(energy1)))))
    {
      NKrylov = 5 + cycles % (NKrylovMax - 5);
      M.N = NKrylov;
      
      cycles++;
      if ((cycles & 15) == 15)
      {
        printf("Reached %i-th cycle with energy %22.15E\n", cycles, energy1);
        fflush(stdout);
      }
      energy2 = energy1;
      ComplexWaveFunctionMultiply(v[0], 0);
      ComplexWaveFunctionAddSameBasis(v[0], *psi);
      OperatorPsiConserveBasisPreAllocatedCC(&H, &(v[0]), &(v[1]));
      ComplexWaveFunctionDotProduct(&v[1], &v[0], &(M.a[0]), &Imzero);
      for (i = 1; i < M.N; i++)
      {
        vip = (i + 1) % 3;
        vi = (i) % 3;
        vim = (i - 1) % 3;
        ComplexWaveFunctionAddMultipliedConserveBasis(&(v[vi]), v[vim], -M.a[i]);
        M.b[i] = ComplexWaveFunctionCalculateNorm(v[vi]);
        if (M.b[i] > DBL_MIN * 10)
          ComplexWaveFunctionMultiply(v[vi], (double)1.0 / M.b[i]);
        else
          ComplexWaveFunctionMultiply(v[vi], 0.0);
        ComplexWaveFunctionOrthogonalizeConserveBasis(&v[vi], OtherPsi, NOtherPsi);
        OperatorPsiConserveBasisPreAllocatedCC(&H, &(v[vi]), &(v[vip]));
        ComplexWaveFunctionDotProduct(&v[vip], &v[vi], &(M.a[i + 1]), &Imzero);
        ComplexWaveFunctionAddMultipliedConserveBasis(&(v[vip]), v[vim], -M.b[i]);
      }
      TriDiagonalMatrixToEigensystem(M, &Eigen, val);
      ComplexWaveFunctionMultiply(v[0], 0);
      ComplexWaveFunctionAddSameBasis(v[0], *psi);
      ComplexWaveFunctionMultiply(*psi, Eigen.VRe[0][0]);
      OperatorPsiConserveBasisPreAllocatedCC(&H, &(v[0]), &(v[1]));
      for (i = 1; i < M.N; i++)
      {
        vip = (i + 1) % 3;
        vi = (i) % 3;
        vim = (i - 1) % 3;
        ComplexWaveFunctionAddMultipliedSameBasis(v[vi], v[vim], -M.a[i]);
        if (M.b[i] > DBL_MIN * 10)
          ComplexWaveFunctionMultiply(v[vi], (double)1.0 / M.b[i]);
        else
          ComplexWaveFunctionMultiply(v[vi], 0.0);
        ComplexWaveFunctionOrthogonalizeConserveBasis(&v[vi], OtherPsi, NOtherPsi);
        OperatorPsiConserveBasisPreAllocatedCC(&H, &(v[vi]), &(v[vip]));
        ComplexWaveFunctionAddMultipliedSameBasis(v[vip], v[vim], -M.b[i]);
        ComplexWaveFunctionAddMultipliedSameBasis(*psi, v[vi], Eigen.VRe[0][i]);
      }
      ComplexWaveFunctionOrthogonalizeConserveBasis(psi, OtherPsi, NOtherPsi);
      ComplexWaveFunctionNormalize(psi);
      PsiOperatorPsiCHC(H, *psi, &energy1);
    }
    ComplexWaveFunctionFreeValues(&v[0]);
    ComplexWaveFunctionFreeValues(&v[1]);
    ComplexWaveFunctionFreeValues(&v[2]);
    FreeMatrix(&Eigen);
    printf("Energy: %22.15E  ", energy1);
    printf("NumberOfLaczosCycles: %i\n", cycles);
    fflush(stdout);
  }
  return 0;
}


uint32_t LanczosGroundStateRRold2(const OperatorType H, WaveFunctionType * const psi, double *epsilon, const double zero, double shiftin, WaveFunctionType * OtherPsi, uint32_t NOtherPsi)
{
  printf("\nStart of LanczosGroundStateRR\n\n");
  fflush(stdout);
  const uint32_t NKrylovMax = 100;
  uint32_t NKrylov;
  WaveFunctionType v[3];
  double energy = 4 * zero, energy1 = 2 * zero, energy2;
  double shift;
  uint32_t cycles, i, vi, vip, vim;
  MatrixType Eigen;
  TriDiagonalMatrixType M;
  double val[NKrylovMax];
  M.N = NKrylovMax;
  InitTriDiagonalMatrix(&M);
  while ((energy - energy1) > (zero * (1.0 + fabs(energy))))
  {
    cycles = 0;
    energy = energy1;
    shift = -shiftin - energy * 0.49;
  Repeat:
    if (OperatorPsiRR(H, *psi, &v[0]))
      goto OutOffMemv0;
    RealWaveFunctionMultiply(*psi, shift);
    if (RealWaveFunctionAddUnionBasis(&v[0], *psi))
      goto OutOffMemv0add;
    if (RealWaveFunctionOrthogonalize(&v[0], OtherPsi, NOtherPsi))
      goto OutOffMemv0Ortho;
    RealWaveFunctionNormalize(&v[0]);
    printf("BasisSetSize After Expansion: %12i Truncation: ", v[0].N);
    RealWaveFunctionTruncateBasis(&v[0], *epsilon);
    printf("%12i\n", v[0].N);
    fflush(stdout);
    if (OperatorPsiConserveBasisRR(H, v[0], &v[1]))
      goto OutOffMemv1;
    RealWaveFunctionMultiply(v[0], shift);
    RealWaveFunctionAddSameBasis(v[1], v[0]);
    RealWaveFunctionNormalize(&v[1]);
    if (OperatorPsiConserveBasisRR(H, v[1], &v[2]))
      goto OutOffMemv2;
    RealWaveFunctionMultiply(v[1], shift);
    RealWaveFunctionAddSameBasis(v[2], v[1]);
    RealWaveFunctionNormalize(&v[2]);
    RealWaveFunctionFree(psi);
    if (OperatorPsiConserveBasisRR(H, v[2], psi))
      goto OutOffMempsi;
    RealWaveFunctionMultiply(v[2], shift);
    RealWaveFunctionAddSameBasis(*psi, v[2]);
    RealWaveFunctionOrthogonalizeConserveBasis(psi, OtherPsi, NOtherPsi);
    RealWaveFunctionNormalize(psi);
    goto EnoughMem;
  OutOffMempsi:
    RealWaveFunctionFreeValues(&v[2]);
  OutOffMemv2:
    RealWaveFunctionFreeValues(&v[1]);
  OutOffMemv1:
  OutOffMemv0Ortho:
  OutOffMemv0add:
    RealWaveFunctionFreeValues(&v[0]);
    RealWaveFunctionMultiply(*psi, 1.0 / shift);
  OutOffMemv0:
    printf("Out off memory increase epsilon from %22.15E to %22.15E\n", *epsilon, *epsilon * 2);
    *epsilon *= 2.0;
    if (psi->NBlok == 1)
    {
      printf("Outof memory in LanczosGroundStateRR,\npsi->NBLok=1 so no way to recover.\n");
      AllocFailed();
    }
    printf("BasisSetSize Before Truncation: %12i ", psi->N);
    RealWaveFunctionTruncateBasis(psi, *epsilon);
    printf("After Truncation: %12i\n", psi->N);
    fflush(stdout);
    goto Repeat;
  EnoughMem:
    
    energy1 = 2 * zero;
    energy2 = 4 * zero;
    while ((cycles < 3) || ((energy2 - energy1) > (zero * (1.0 + fabs(energy1)))))
    {
      NKrylov = 5 + cycles % (NKrylovMax - 5);
      M.N = NKrylov;
      
      cycles++;
      if ((cycles & 15) == 15)
      {
        printf("Reached %i-th cycle with energy %22.15E\n", cycles, energy1);
        fflush(stdout);
      }
      energy2 = energy1;
      RealWaveFunctionMultiply(v[0], 0);
      RealWaveFunctionAddSameBasis(v[0], *psi);
      OperatorPsiConserveBasisPreAllocatedRR(&H, &(v[0]), &(v[1]));
      M.a[1] = RealWaveFunctionDotProduct(&v[1], &v[0]);
      for (i = 1; i < M.N; i++)
      {
        vip = (i + 1) % 3;
        vi = (i) % 3;
        vim = (i - 1) % 3;
        RealWaveFunctionAddMultipliedConserveBasis(&(v[vi]), v[vim], -M.a[i]);
        M.b[i] = RealWaveFunctionCalculateNorm(v[vi]);
        if (M.b[i] > DBL_MIN * 10)
          RealWaveFunctionMultiply(v[vi], (double)1.0 / M.b[i]);
        else
          RealWaveFunctionMultiply(v[vi], 0.0);
        RealWaveFunctionOrthogonalizeConserveBasis(&v[vi], OtherPsi, NOtherPsi);
        OperatorPsiConserveBasisPreAllocatedRR(&H, &(v[vi]), &(v[vip]));
        M.a[i + 1] = RealWaveFunctionDotProduct(&v[vip], &v[vi]);
        RealWaveFunctionAddMultipliedConserveBasis(&(v[vip]), v[vim], -M.b[i]);
      }
      TriDiagonalMatrixToEigensystem(M, &Eigen, val);
      RealWaveFunctionMultiply(v[0], 0);
      RealWaveFunctionAddSameBasis(v[0], *psi);
      RealWaveFunctionMultiply(*psi, Eigen.VRe[0][0]);
      OperatorPsiConserveBasisPreAllocatedRR(&H, &(v[0]), &(v[1]));
      for (i = 1; i < M.N; i++)
      {
        vip = (i + 1) % 3;
        vi = (i) % 3;
        vim = (i - 1) % 3;
        RealWaveFunctionAddMultipliedSameBasis(v[vi], v[vim], -M.a[i]);
        if (M.b[i] > DBL_MIN * 10)
          RealWaveFunctionMultiply(v[vi], (double)1.0 / M.b[i]);
        else
          RealWaveFunctionMultiply(v[vi], 0.0);
        RealWaveFunctionOrthogonalizeConserveBasis(&v[vi], OtherPsi, NOtherPsi);
        OperatorPsiConserveBasisPreAllocatedRR(&H, &(v[vi]), &(v[vip]));
        RealWaveFunctionAddMultipliedSameBasis(v[vip], v[vim], -M.b[i]);
        RealWaveFunctionAddMultipliedSameBasis(*psi, v[vi], Eigen.VRe[0][i]);
      }
      RealWaveFunctionOrthogonalizeConserveBasis(psi, OtherPsi, NOtherPsi);
      RealWaveFunctionNormalize(psi);
      PsiOperatorPsiRR(H, *psi, &energy1);
    }
    RealWaveFunctionFreeValues(&v[0]);
    RealWaveFunctionFreeValues(&v[1]);
    RealWaveFunctionFreeValues(&v[2]);
    FreeMatrix(&Eigen);
    printf("Energy: %22.15E  ", energy1);
    printf("NumberOfLaczosCycles: %i\n", cycles);
    fflush(stdout);
  }
  return 0;
}


/*  add the scaled normal ordered product of two operators to AB
 *  AB <- AB + alpha A B
 *  AB must be initiallized, i.e AB.Maxlength must be A.Maxlength+B.Maxlength
 *  it is advised to set AB.l[il].NBitsKey to sensible values
 */
uint32_t OperatorMultiplyScaledNormalOrder(OperatorType * const AB, OperatorType const A, OperatorType const B,
                                           const double alphaRe, const double alphaIm,
                                           uint32_t(*QDetNotOrdered) (uint16_t, uint16_t, uint32_t))
{
  fflush(stdout);
  uint32_t ilA, ilB, iNA, iNB, k;
  uint16_t det[A.MaxLength + B.MaxLength];
  
  if (AB->MaxLength != A.MaxLength + B.MaxLength)
  {
    printf("Length of output operator does not match in OperatorMultiplyScaledNormalOrder\n");
    fflush(stdout);
    return 1;
  }
  
  for (ilA = 0; ilA <= A.MaxLength; ilA++)
  {
    for (ilB = 0; ilB <= B.MaxLength; ilB++)
    {
      if ((AB->l[ilA + ilB].QComplex) || (A.l[ilA].QComplex) || (B.l[ilB].QComplex) || (alphaIm != 0))
      {                         /* at least one complex */
        if (!(AB->l[ilA + ilB].QComplex))
          if (MakeRealOperatorFixedLengthComplex(AB, ilA + ilB))
          {
            printf("MakeRealOperatorFixedLengthComplex failed in OperatorMultiplyNormalOrder\n");
            fflush(stdout);
            return 1;
          }
        if (A.l[ilA].QComplex)
        {                       /* A complex */
          if (B.l[ilB].QComplex)
          {                     /* A and B complex */
            for (iNA = 0; iNA < A.l[ilA].N; iNA++)
            {
              for (k = 0; k < ilA; k++)
              {
                det[k] = A.l[ilA].o[iNA * ilA + k];
              }
              for (iNB = 0; iNB < B.l[ilB].N; iNB++)
              {
                for (k = 0; k < ilB; k++)
                {
                  det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                }
                if (ComplexOperatorAddLaderNormalOrder
                    (AB, det, ilA + ilB,
                     alphaRe * (A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Re - A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Im) -
                     alphaIm * (A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Im + A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Re),
                     alphaRe * (A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Im + A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Re) +
                     alphaIm * (A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Re - A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Im),
                     QDetNotOrdered))
                {
                  printf("ComplexOperatorAddLaderNormalOrder failed in OperatorMultiplyNormalOrder\n");
                  fflush(stdout);
                  return 1;
                }
              }
            }
          }
          else
          {                     /* A complex B real */
            for (iNA = 0; iNA < A.l[ilA].N; iNA++)
            {
              for (k = 0; k < ilA; k++)
              {
                det[k] = A.l[ilA].o[iNA * ilA + k];
              }
              for (iNB = 0; iNB < B.l[ilB].N; iNB++)
              {
                for (k = 0; k < ilB; k++)
                {
                  det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                }
                if (ComplexOperatorAddLaderNormalOrder
                    (AB, det, ilA + ilB,
                     alphaRe * (A.l[ilA].V[iNA].Re * B.l[ilB].VRe[iNB]) -
                     alphaIm * (A.l[ilA].V[iNA].Im * B.l[ilB].VRe[iNB]),
                     alphaRe * (A.l[ilA].V[iNA].Im * B.l[ilB].VRe[iNB]) +
                     alphaIm * (A.l[ilA].V[iNA].Re * B.l[ilB].VRe[iNB]), QDetNotOrdered))
                {
                  printf("ComplexOperatorAddLaderNormalOrder failed in OperatorMultiplyNormalOrder\n");
                  fflush(stdout);
                  return 1;
                }
              }
            }
          }
        }
        else
        {                       /* A real */
          if (B.l[ilB].QComplex)
          {                     /* A real and B complex */
            for (iNA = 0; iNA < A.l[ilA].N; iNA++)
            {
              for (k = 0; k < ilA; k++)
              {
                det[k] = A.l[ilA].o[iNA * ilA + k];
              }
              for (iNB = 0; iNB < B.l[ilB].N; iNB++)
              {
                for (k = 0; k < ilB; k++)
                {
                  det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                }
                if (ComplexOperatorAddLaderNormalOrder
                    (AB, det, ilA + ilB,
                     alphaRe * (A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Re) -
                     alphaIm * (A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Im),
                     alphaRe * (A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Im) +
                     alphaIm * (A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Re), QDetNotOrdered))
                {
                  printf("ComplexOperatorAddLaderNormalOrder failed in OperatorMultiplyNormalOrder\n");
                  fflush(stdout);
                  return 1;
                }
              }
            }
          }
          else
          {                     /* A and B real */
            for (iNA = 0; iNA < A.l[ilA].N; iNA++)
            {
              for (k = 0; k < ilA; k++)
              {
                det[k] = A.l[ilA].o[iNA * ilA + k];
              }
              for (iNB = 0; iNB < B.l[ilB].N; iNB++)
              {
                for (k = 0; k < ilB; k++)
                {
                  det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                }
                if (ComplexOperatorAddLaderNormalOrder
                    (AB, det, ilA + ilB, alphaRe * (A.l[ilA].VRe[iNA] * B.l[ilB].VRe[iNB]),
                     alphaIm * (A.l[ilA].VRe[iNA] * B.l[ilB].VRe[iNB]), QDetNotOrdered))
                {
                  printf("ComplexOperatorAddLaderNormalOrder failed in OperatorMultiplyNormalOrder\n");
                  fflush(stdout);
                  return 1;
                }
              }
            }
          }
        }
      }
      else
      {                         /* All real */
        for (iNA = 0; iNA < A.l[ilA].N; iNA++)
        {
          for (k = 0; k < ilA; k++)
          {
            det[k] = A.l[ilA].o[iNA * ilA + k];
          }
          for (iNB = 0; iNB < B.l[ilB].N; iNB++)
          {
            for (k = 0; k < ilB; k++)
            {
              det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
            }
            if (RealOperatorAddLaderNormalOrder
                (AB, det, ilA + ilB, alphaRe * A.l[ilA].VRe[iNA] * B.l[ilB].VRe[iNB], QDetNotOrdered))
            {
              printf("RealOperatorAddLaderNormalOrder failed in OperatorMultiplyNormalOrder\n");
              fflush(stdout);
              return 1;
            }
          }
        }
      }
    }
  }
  SetQComplexOperator(AB);
  return 0;
}





/*  add the product of two operators to AB
 *  AB <- AB + A B
 *  AB must be initiallized, i.e AB.Maxlength must be A.Maxlength+B.Maxlength
 *  it is advised to set AB.l[il].NBitsKey to sensible values
 */
uint32_t OperatorMultiply(OperatorType * const AB, OperatorType const A, OperatorType const B)
{
  uint32_t ilA, ilB, iNA, iNB, k;
  uint16_t det[A.MaxLength + B.MaxLength];
  
  if (AB->MaxLength != A.MaxLength + B.MaxLength)
  {
    printf("Length of output operator does not match in OperatorMultiply\n");
    fflush(stdout);
    return 1;
  }
  
  for (ilA = 0; ilA <= A.MaxLength; ilA++)
  {
    for (ilB = 0; ilB <= B.MaxLength; ilB++)
    {
      if ((AB->l[ilA + ilB].QComplex) || (A.l[ilA].QComplex) || (B.l[ilB].QComplex))
      {                         /* at least one complex */
        if (!(AB->l[ilA + ilB].QComplex))
          if (MakeRealOperatorFixedLengthComplex(AB, ilA + ilB))
          {
            printf("MakeRealOperatorFixedLengthComplex at length=%i failed in OperatorMultiply\n", ilA + ilB);
            fflush(stdout);
            return 1;
          }
        if (A.l[ilA].QComplex)
        {                       /* A complex */
          if (B.l[ilB].QComplex)
          {                     /* A and B complex */
            for (iNA = 0; iNA < A.l[ilA].N; iNA++)
            {
              for (k = 0; k < ilA; k++)
              {
                det[k] = A.l[ilA].o[iNA * ilA + k];
              }
              for (iNB = 0; iNB < B.l[ilB].N; iNB++)
              {
                for (k = 0; k < ilB; k++)
                {
                  det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                }
                if (ComplexOperatorAddLader
                    (AB, det, ilA + ilB,
                     A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Re - A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Im,
                     A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Im + A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Re))
                {
                  printf("ComplexOperatorAddLader failed in OperatorMultiply\n");
                  fflush(stdout);
                  return 1;
                }
              }
            }
          }
          else
          {                     /* A complex B real */
            for (iNA = 0; iNA < A.l[ilA].N; iNA++)
            {
              for (k = 0; k < ilA; k++)
              {
                det[k] = A.l[ilA].o[iNA * ilA + k];
              }
              for (iNB = 0; iNB < B.l[ilB].N; iNB++)
              {
                for (k = 0; k < ilB; k++)
                {
                  det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                }
                if (ComplexOperatorAddLader
                    (AB, det, ilA + ilB, A.l[ilA].V[iNA].Re * B.l[ilB].VRe[iNB],
                     A.l[ilA].V[iNA].Im * B.l[ilB].VRe[iNB]))
                {
                  printf("ComplexOperatorAddLader failed in OperatorMultiply\n");
                  fflush(stdout);
                  return 1;
                }
              }
            }
          }
        }
        else
        {                       /* A real */
          if (B.l[ilB].QComplex)
          {                     /* A real and B complex */
            for (iNA = 0; iNA < A.l[ilA].N; iNA++)
            {
              for (k = 0; k < ilA; k++)
              {
                det[k] = A.l[ilA].o[iNA * ilA + k];
              }
              for (iNB = 0; iNB < B.l[ilB].N; iNB++)
              {
                for (k = 0; k < ilB; k++)
                {
                  det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                }
                if (ComplexOperatorAddLader
                    (AB, det, ilA + ilB, A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Re,
                     A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Im))
                {
                  printf("ComplexOperatorAddLader failed in OperatorMultiply\n");
                  fflush(stdout);
                  return 1;
                }
              }
            }
          }
          else
          {                     /* A and B real */
            for (iNA = 0; iNA < A.l[ilA].N; iNA++)
            {
              for (k = 0; k < ilA; k++)
              {
                det[k] = A.l[ilA].o[iNA * ilA + k];
              }
              for (iNB = 0; iNB < B.l[ilB].N; iNB++)
              {
                for (k = 0; k < ilB; k++)
                {
                  det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                }
                if (ComplexOperatorAddLader(AB, det, ilA + ilB, A.l[ilA].VRe[iNA] * B.l[ilB].VRe[iNB], 0))
                {
                  printf("ComplexOperatorAddLader failed in OperatorMultiply\n");
                  fflush(stdout);
                  return 1;
                }
              }
            }
          }
        }
      }
      else
      {                         /* All real */
        for (iNA = 0; iNA < A.l[ilA].N; iNA++)
        {
          for (k = 0; k < ilA; k++)
          {
            det[k] = A.l[ilA].o[iNA * ilA + k];
          }
          for (iNB = 0; iNB < B.l[ilB].N; iNB++)
          {
            for (k = 0; k < ilB; k++)
            {
              det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
            }
            if (RealOperatorAddLader(AB, det, ilA + ilB, A.l[ilA].VRe[iNA] * B.l[ilB].VRe[iNB]))
            {
              printf("RealOperatorAddLader failed in OperatorMultiply\n");
              fflush(stdout);
              return 1;
            }
          }
        }
      }
    }
  }
  SetQComplexOperator(AB);
  return 0;
}


/*  add the normal ordered product of two operators to AB
 *  AB <- AB + A B
 *  AB must be initiallized, i.e AB.Maxlength must be A.Maxlength+B.Maxlength
 *  it is advised to set AB.l[il].NBitsKey to sensible values
 */
uint32_t OperatorMultiplyNormalOrder(OperatorType * const AB, OperatorType const A, OperatorType const B,
                                     uint32_t(*QDetNotOrdered) (uint16_t, uint16_t, uint32_t))
{
  fflush(stdout);
  uint32_t ilA, ilB, iNA, iNB, k;
  uint16_t det[A.MaxLength + B.MaxLength];
  if (AB->MaxLength != A.MaxLength + B.MaxLength)
  {
    printf("Length of output operator does not match in OperatorMultiplyNormalOrder\n");
    fflush(stdout);
    return 1;
  }
  
  for (ilA = 0; ilA <= A.MaxLength; ilA++)
  {
    for (ilB = 0; ilB <= B.MaxLength; ilB++)
    {
      if ((AB->l[ilA + ilB].QComplex) || (A.l[ilA].QComplex) || (B.l[ilB].QComplex))
      {                         /* at least one complex */
        if (!(AB->l[ilA + ilB].QComplex))
          if (MakeRealOperatorFixedLengthComplex(AB, ilA + ilB))
          {
            printf("MakeRealOperatorFixedLengthComplex failed in OperatorMultiplyNormalOrder\n");
            fflush(stdout);
            return 1;
          }
        if (A.l[ilA].QComplex)
        {                       /* A complex */
          if (B.l[ilB].QComplex)
          {                     /* A and B complex */
            for (iNA = 0; iNA < A.l[ilA].N; iNA++)
            {
              for (k = 0; k < ilA; k++)
              {
                det[k] = A.l[ilA].o[iNA * ilA + k];
              }
              for (iNB = 0; iNB < B.l[ilB].N; iNB++)
              {
                for (k = 0; k < ilB; k++)
                {
                  det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                }
                if (ComplexOperatorAddLaderNormalOrder
                    (AB, det, ilA + ilB,
                     A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Re - A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Im,
                     A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Im + A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Re, QDetNotOrdered))
                {
                  printf("ComplexOperatorAddLaderNormalOrder failed in OperatorMultiplyNormalOrder\n");
                  fflush(stdout);
                  return 1;
                }
              }
            }
          }
          else
          {                     /* A complex B real */
            for (iNA = 0; iNA < A.l[ilA].N; iNA++)
            {
              for (k = 0; k < ilA; k++)
              {
                det[k] = A.l[ilA].o[iNA * ilA + k];
              }
              for (iNB = 0; iNB < B.l[ilB].N; iNB++)
              {
                for (k = 0; k < ilB; k++)
                {
                  det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                }
                if (ComplexOperatorAddLaderNormalOrder
                    (AB, det, ilA + ilB, A.l[ilA].V[iNA].Re * B.l[ilB].VRe[iNB], A.l[ilA].V[iNA].Im * B.l[ilB].VRe[iNB],
                     QDetNotOrdered))
                {
                  printf("ComplexOperatorAddLaderNormalOrder failed in OperatorMultiplyNormalOrder\n");
                  fflush(stdout);
                  return 1;
                }
              }
            }
          }
        }
        else
        {                       /* A real */
          if (B.l[ilB].QComplex)
          {                     /* A real and B complex */
            for (iNA = 0; iNA < A.l[ilA].N; iNA++)
            {
              for (k = 0; k < ilA; k++)
              {
                det[k] = A.l[ilA].o[iNA * ilA + k];
              }
              for (iNB = 0; iNB < B.l[ilB].N; iNB++)
              {
                for (k = 0; k < ilB; k++)
                {
                  det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                }
                if (ComplexOperatorAddLaderNormalOrder
                    (AB, det, ilA + ilB, A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Re, A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Im,
                     QDetNotOrdered))
                {
                  printf("ComplexOperatorAddLaderNormalOrder failed in OperatorMultiplyNormalOrder\n");
                  fflush(stdout);
                  return 1;
                }
              }
            }
          }
          else
          {                     /* A and B real */
            for (iNA = 0; iNA < A.l[ilA].N; iNA++)
            {
              for (k = 0; k < ilA; k++)
              {
                det[k] = A.l[ilA].o[iNA * ilA + k];
              }
              for (iNB = 0; iNB < B.l[ilB].N; iNB++)
              {
                for (k = 0; k < ilB; k++)
                {
                  det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                }
                if (ComplexOperatorAddLaderNormalOrder
                    (AB, det, ilA + ilB, A.l[ilA].VRe[iNA] * B.l[ilB].VRe[iNB], 0, QDetNotOrdered))
                {
                  printf("ComplexOperatorAddLaderNormalOrder failed in OperatorMultiplyNormalOrder\n");
                  fflush(stdout);
                  return 1;
                }
              }
            }
          }
        }
      }
      else
      {                         /* All real */
        for (iNA = 0; iNA < A.l[ilA].N; iNA++)
        {
          for (k = 0; k < ilA; k++)
          {
            det[k] = A.l[ilA].o[iNA * ilA + k];
          }
          for (iNB = 0; iNB < B.l[ilB].N; iNB++)
          {
            for (k = 0; k < ilB; k++)
            {
              det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
            }
            if (RealOperatorAddLaderNormalOrder
                (AB, det, ilA + ilB, A.l[ilA].VRe[iNA] * B.l[ilB].VRe[iNB], QDetNotOrdered))
            {
              printf("RealOperatorAddLaderNormalOrder failed in OperatorMultiplyNormalOrder\n");
              fflush(stdout);
              return 1;
            }
          }
        }
      }
    }
  }
  SetQComplexOperator(AB);
  return 0;
}



int mainOld2(int argc, const char *argv[])
{
  const uint32_t imuMax = 20, iHistMax = 20, StringLengthMax = 255;
  TimeType Time;
  DMFTSettings HSet;
  WaveFunctionType Psi, OPsi, PsiTemp;
  TriDiagonalMatrixType G0, GBath, GclPup, GclPdn, GclMup, GclMdn, Gcl, GBathNew, GBathOld, gPerturbation[4], Sigma;
  BADoubleTriDiagonalMatrixType GAnd, GAndTrunk, BAG0;
  AndersonMatrixType GBathAnd, GclAnd, G0And, GBathNewAnd;
  ListOfPoles SigmaPoles, GclList, GListTemp;
  OperatorType HAnd, HAndGrd, OPESdn, OPESup, OIPESdn, OIPESup, VPerturbation[4], Sx, Sy, Sz, S2;
  MatrixType T;
  SpectraType Spectra[10];
  char   GBathFileName[StringLengthMax], GClFileName[StringLengthMax], GSigmaFileName[StringLengthMax],
    PsiFileName[StringLengthMax], SpecFileName[StringLengthMax], Buffer[StringLengthMax];
  uint32_t iOuter, i, PsiNBathV, PsiNBathC, imu;
  uint16_t oadd[1], *SpinUp, *SpinDn;
  double Ef = 0.0, epsilonAnd = 0.00001, mu[imuMax + 1], G0Occupation, GClOccupation[imuMax], G0OccupationReal;
  double PsiEnergy, PsiQuality, PsiHist[iHistMax], *PsiRho;

  InitializeTiming(&Time);
  BeginTiming("Total", &Time);
  TimeStampAll();

  /* read the DMFT settings */
  {
    argc -= 1;
    if (FileReadDMFTSettings(&HSet, argv[argc]))
    {
      printf("Fatal Error in main\n Failed to read DMFT Settings from file %s\n", argv[argc]);
      return 1;
    };
    PrintDMFTSettings(HSet);
  }
  /* read G0 */
  {
    if (FileReadPlainTriDiagonalMatrixType(&G0, HSet.G0FileName))
    {
      printf("Fatal Error in main\n Failed to read G0 from file %s\n", HSet.G0FileName);
      return 1;
    };
    strcpy(G0.Name, "G0");
    if ((HSet.SymOpt) & 2)
    {
      printf("Symmetrize G0(w) -= G0(-w)^*\n");
      SymmetrizeTridiagonalMatrix(&G0);
    }
    TriDiagonalMatrixToBADoubleTriDiagonalMatrix(G0, &BAG0, G0.N, G0.N, Ef, epsilonAnd);
    if ((HSet.SymOpt) & 2)
    {
      printf("Symmetrize G0(w) -= G0(-w)^*\n");
      SymmetrizeBADoubleTriDiagonalMatrix(&BAG0);
    }
    TriDiagonalMatrixToAndersonMatrix(G0, &T, &G0And);
    if ((HSet.SymOpt) & 2)
    {
      printf("Symmetrize G0(w) -= G0(-w)^*\n");
      SymmetrizeAndersonMatrix(&G0And);
    }
    FreeMatrix(&T);
    //G0Occupation=2*Sqr(BAG0.val.b[0]);
    G0Occupation = HSet.Nelectron;
    printf("\nRead G0 from file \"%s\"\n", HSet.G0FileName);
    printf("The occupation of the cluster is %22.15E electrons\n", G0Occupation);
  }
  /* read G_Bath_iStart and psi_iStart */
  {
    if (HSet.iStart == 0)
    {
      sprintf(GBathFileName, "%s000", HSet.GRootFileName);
      sprintf(PsiFileName, "%s000", HSet.WaveFunctionRootFileName);
    }
    else if (HSet.iStart < 10)
    {
      sprintf(GBathFileName, "%s00%1i", HSet.GRootFileName, HSet.iStart);
      sprintf(PsiFileName, "%s00%1i", HSet.WaveFunctionRootFileName, HSet.iStart);
    }
    else if (HSet.iStart < 100)
    {
      sprintf(GBathFileName, "%s0%2i", HSet.GRootFileName, HSet.iStart);
      sprintf(PsiFileName, "%s0%2i", HSet.WaveFunctionRootFileName, HSet.iStart);
    }
    else
    {
      sprintf(GBathFileName, "%s%i", HSet.GRootFileName, HSet.iStart);
      sprintf(PsiFileName, "%s%i", HSet.WaveFunctionRootFileName, HSet.iStart);
    }
    if (FileReadPlainTriDiagonalMatrixType(&GBath, GBathFileName))
    {
      printf("Fatal Error in main\n Failed to read GBath from file %s !\n", GBathFileName);
      return 1;
    };
    if (file_exists(PsiFileName))
    {
      Psi.NBitsKey = HSet.NBitsKey;
      if (FileReadWaveFunctionPlainBinary(&Psi, PsiFileName))
      {
        printf("Fatal Error in main\n Failed to read Psi from file %s !\n", PsiFileName);
        return 1;
      };
      printf("\nRead the starting wave function from the file %s\n", PsiFileName);
      sscanf(Psi.Name, "%s %i %i", Buffer, &PsiNBathC, &PsiNBathV);
      printf
        ("Wave-function \"%s\" found has %i conduction bath sites, %i valence bath sites and %i Fermions in basis\n",
         Psi.Name, PsiNBathC, PsiNBathV, Psi.NFermion);
      if ((2 * (PsiNBathC + PsiNBathV + 2)) != Psi.NFermion)
      {
        printf("The total number of Fermions is incompatible witht the bath sizes\n");
        printf("This wave-function can not be used to start from\n");
        return 1;
      }
    }
    else
    {
      printf("\nWave-function file: '%s' does not excist\n Started from new function\n", PsiFileName);
      CreateWaveFunctionDMFTSecond(&Psi, HSet.NBath, HSet.NBath, HSet.NBitsKey);
    }
    strcpy(GBath.Name, "Gbath");
  }


  /* do the DMFT loops */
  for (iOuter = HSet.iStart; iOuter < HSet.iStart + HSet.Ni; iOuter++)
  {
    BeginTiming("DMFTLoop", &Time);
    printf("\n\n***********************************************************************\n");
    printf("Start of loop %i from %i to %i\n", iOuter + 1, HSet.iStart, HSet.iStart + HSet.Ni);
    printf("***********************************************************************\n\n\n");
    /* Create splitted anderson matrix */
    {
      //TriDiagonalMatrixToBADoubleTriDiagonalMatrixOccupation(GBath, &GAnd, GBath.N, GBath.N, Ef, G0Occupation/2);
      TriDiagonalMatrixToBADoubleTriDiagonalMatrix(GBath, &GAnd, GBath.N, GBath.N, 0, epsilonAnd);
      if ((HSet.SymOpt) & 2)
      {
        printf("Symmetrize GAnd(w) -= GAnd(-w)^*\n");
        printf("Size of GAnd went from %i + %i ", GAnd.val.N, GAnd.con.N);
        SymmetrizeBADoubleTriDiagonalMatrix(&GAnd);
        printf("to %i + %i ", GAnd.val.N, GAnd.con.N);
        GAnd.val.N = Min(GAnd.val.N, HSet.NGBath);
        GAnd.con.N = Min(GAnd.con.N, HSet.NGBath);
        printf("to %i + %i\n", GAnd.val.N, GAnd.con.N);
      }
      if ((HSet.SymOpt) & 16)
      {
        printf("Enforce GAnd to fullfill the knonw sumrules\n");
        printf("Not yet implemented\n");
      }
    }
    /* Create GAnd and GAndTrunk as a copy pointing to the same memory */
    {
      GAndTrunk.Ef = GAnd.Ef;
      GAndTrunk.Epsilon = GAnd.Epsilon;
      GAndTrunk.a0 = GAnd.a0;
      GAndTrunk.b0 = GAnd.b0;
      GAndTrunk.acvB = GAnd.acvB;
      GAndTrunk.acvA = GAnd.acvA;
      GAndTrunk.acvBA = GAnd.acvBA;
      GAndTrunk.bvB0 = GAnd.bvB0;
      GAndTrunk.bcB0 = GAnd.bcB0;
      GAndTrunk.bvA0 = GAnd.bvA0;
      GAndTrunk.bcA0 = GAnd.bcA0;
      GAndTrunk.val.N = GAnd.val.N;
      GAndTrunk.val.a = GAnd.val.a;
      GAndTrunk.val.b = GAnd.val.b;
      GAndTrunk.con.N = GAnd.con.N;
      GAndTrunk.con.a = GAnd.con.a;
      GAndTrunk.con.b = GAnd.con.b;
    }
    /* Create VPerturbation */
    {
      if (HSet.Dyson)
      {
        for (i = 0; i < 4; i++)
        {
          InitOperator(&(VPerturbation[i]));
          VPerturbation[i].NFermion = 4 + 2 * (Min(-1 + GAnd.val.N, HSet.NBath) + Min(-1 + GAnd.con.N, HSet.NBath));
          VPerturbation[i].NBoson = 0;
          VPerturbation[i].NWord = Min(VPerturbation[i].NFermion, 1 + (VPerturbation[i].NFermion - 1) / 8);
          VPerturbation[i].QComplex = 0;
          InitOperatorFixedLength(&(VPerturbation[i]), 1, 1);
        }
        if (GAnd.val.N > (HSet.NBath + 1))
        {
          oadd[0] = CrAnMask + 3 + HSet.NBath;
          RealOperatorAddLader(&(VPerturbation[0]), oadd, 1, GAnd.val.b[HSet.NBath + 1]);
          oadd[0] = CrAnMask + 3 + Min(-1 + GAnd.con.N, HSet.NBath) + 2 * HSet.NBath;
          RealOperatorAddLader(&(VPerturbation[1]), oadd, 1, GAnd.val.b[HSet.NBath + 1]);
        }
        if (GAnd.con.N > HSet.NBath + 1)
        {
          oadd[0] = 3 + Min(-1 + GAnd.val.N, HSet.NBath) + HSet.NBath;
          RealOperatorAddLader(&(VPerturbation[2]), oadd, 1, GAnd.con.b[HSet.NBath + 1]);
          oadd[0] = 3 + 2 * Min(-1 + GAnd.val.N, HSet.NBath) + 2 * HSet.NBath;
          RealOperatorAddLader(&(VPerturbation[3]), oadd, 1, GAnd.con.b[HSet.NBath + 1]);
        }
        /* Create gPerturbation */
        if (GAnd.val.N > HSet.NBath + 1)
        {
          gPerturbation[0].N = GAnd.val.N - HSet.NBath - 1;
          gPerturbation[1].N = GAnd.val.N - HSet.NBath - 1;
          InitTriDiagonalMatrix(&(gPerturbation[0]));
          InitTriDiagonalMatrix(&(gPerturbation[1]));
          for (i = 1; i < GAnd.val.N - HSet.NBath - 1; i++)
          {
            gPerturbation[0].a[i] = -GAnd.val.a[i + HSet.NBath + 1];
            gPerturbation[0].b[i] = GAnd.val.b[i + HSet.NBath + 1];
            gPerturbation[1].a[i] = -GAnd.val.a[i + HSet.NBath + 1];
            gPerturbation[1].b[i] = GAnd.val.b[i + HSet.NBath + 1];
          }
          gPerturbation[0].a[GAnd.val.N - HSet.NBath - 1] = -GAnd.val.a[GAnd.val.N];
          gPerturbation[1].a[GAnd.val.N - HSet.NBath - 1] = -GAnd.val.a[GAnd.val.N];
        }
        else
        {
          gPerturbation[0].N = 1;
          gPerturbation[1].N = 1;
          InitTriDiagonalMatrix(&(gPerturbation[0]));
          InitTriDiagonalMatrix(&(gPerturbation[1]));
          gPerturbation[0].a[1] = 0;
          gPerturbation[1].a[1] = 0;
        }
        if (GAnd.con.N > HSet.NBath + 1)
        {
          gPerturbation[2].N = GAnd.con.N - HSet.NBath - 1;
          gPerturbation[3].N = GAnd.con.N - HSet.NBath - 1;
          InitTriDiagonalMatrix(&(gPerturbation[2]));
          InitTriDiagonalMatrix(&(gPerturbation[3]));
          for (i = 1; i < GAnd.con.N - HSet.NBath - 1; i++)
          {
            gPerturbation[2].a[i] = GAnd.con.a[i + HSet.NBath + 1];
            gPerturbation[2].b[i] = GAnd.con.b[i + HSet.NBath + 1];
            gPerturbation[3].a[i] = GAnd.con.a[i + HSet.NBath + 1];
            gPerturbation[3].b[i] = GAnd.con.b[i + HSet.NBath + 1];
          }
          gPerturbation[2].a[GAnd.con.N - HSet.NBath - 1] = GAnd.con.a[GAnd.con.N];
          gPerturbation[3].a[GAnd.con.N - HSet.NBath - 1] = GAnd.con.a[GAnd.con.N];
        }
        else
        {
          gPerturbation[2].N = 1;
          gPerturbation[3].N = 1;
          InitTriDiagonalMatrix(&(gPerturbation[2]));
          InitTriDiagonalMatrix(&(gPerturbation[3]));
          gPerturbation[2].a[1] = 0;
          gPerturbation[3].a[1] = 0;
        }
      }
    }
    /* Create GAndTrunk */
    {
      if (GAnd.con.N < HSet.NBath + 1)
      {
        printf("Warning: GAnd.con.N(%i)<HSet.NBath+1(%i+1)\n", GAnd.con.N, HSet.NBath);
      }
      else
      {
        GAndTrunk.con.N = HSet.NBath + 1;
      }
      if (GAnd.val.N < HSet.NBath + 1)
      {
        printf("Warning: GAnd.val.N(%i)<HSet.NBath+1(%i+1)\n", GAnd.val.N, HSet.NBath);
      }
      else
      {
        GAndTrunk.val.N = HSet.NBath + 1;
      }
    }
    /* test if psi has correct dimensions and correct if needed */
    {
      sscanf(Psi.Name, "%s %i %i", Buffer, &PsiNBathC, &PsiNBathV);
      if (PsiNBathC != GAndTrunk.con.N - 1)
      {
        printf("\nBefore diagonalizing:\nThe wavefunction has %i conduction bath sites\nThe Anderson Hamiltonian %i\n",
               PsiNBathC, GAndTrunk.con.N - 1);
        if (PsiNBathC > GAndTrunk.con.N - 1)
        {
          printf("Removing %i unoccupied fermions from the wave-function\n", PsiNBathC - GAndTrunk.con.N + 1);
          for (i = 0; i < PsiNBathC - GAndTrunk.con.N + 1; i++)
          {
            WaveFunctionRemoveFermion(Psi, &PsiTemp, Psi.NFermion - 1);
            WaveFunctionFree(&Psi);
            WaveFunctionRemoveFermion(PsiTemp, &Psi, 4 + PsiNBathV + PsiNBathC - i - 1);
            WaveFunctionFree(&PsiTemp);
          }
        }
        else
        {
          printf("Adding %i unoccupied fermions to the wave-function\n", GAndTrunk.con.N - 1 - PsiNBathC);
          for (i = 0; i < GAndTrunk.con.N - 1 - PsiNBathC; i++)
          {
            WaveFunctionAddFermion(Psi, &PsiTemp, Psi.NFermion, 0);
            WaveFunctionFree(&Psi);
            WaveFunctionAddFermion(PsiTemp, &Psi, 4 + PsiNBathV + PsiNBathC + i, 0);
            WaveFunctionFree(&PsiTemp);
          }
        }
      }
      if (PsiNBathV != GAndTrunk.val.N - 1)
      {
        printf("\nBefore diagonalizing:\nThe wavefunction has %i valence bath sites\nThe Anderson Hamiltonian %i\n",
               PsiNBathV, GAndTrunk.val.N - 1);
        if (PsiNBathV > GAndTrunk.val.N - 1)
        {
          printf("Removing %i occupied fermions from the wave-function\n", PsiNBathV - GAndTrunk.val.N + 1);
          for (i = 0; i < PsiNBathV - GAndTrunk.val.N + 1; i++)
          {
            WaveFunctionRemoveFermion(Psi, &PsiTemp, Psi.NFermion - GAndTrunk.con.N);
            WaveFunctionFree(&Psi);
            WaveFunctionRemoveFermion(PsiTemp, &Psi, 4 + PsiNBathV - i - 1);
            WaveFunctionFree(&PsiTemp);
          }
        }
        else
        {
          printf("Adding %i occupied fermions to the wave-function\n", GAndTrunk.val.N - 1 - PsiNBathV);
          for (i = 0; i < GAndTrunk.val.N - 1 - PsiNBathV; i++)
          {
            WaveFunctionAddFermion(Psi, &PsiTemp, Psi.NFermion - GAndTrunk.val.N + 1, 1);
            WaveFunctionFree(&Psi);
            WaveFunctionAddFermion(PsiTemp, &Psi, 4 + PsiNBathV + i, 1);
            WaveFunctionFree(&PsiTemp);
          }
        }
        printf("Corrected basis\n");
        PrintWaveFunctionReadableTruncated(Psi, 0.1);
      }
      sprintf(Psi.Name, "DMFT %i %i", GAndTrunk.con.N - 1, GAndTrunk.val.N - 1);
    }
    /* Print some info on G0 and GAnd */
    {
      printf
        ("                               G0v                    G0c                    GAndv                  Gandc\n");
      printf("b0                          %22.15E %22.15E %22.15E %22.15E\n", BAG0.val.b[0], BAG0.con.b[0],
             GAnd.val.b[0], GAnd.con.b[0]);
      printf("b1                          %22.15E %22.15E %22.15E %22.15E\n", BAG0.val.b[1], BAG0.con.b[1],
             GAnd.val.b[1], GAnd.con.b[1]);
      printf("b2                          %22.15E %22.15E %22.15E %22.15E\n", BAG0.val.b[2], BAG0.con.b[2],
             GAnd.val.b[2], GAnd.con.b[2]);
      printf("a1                          %22.15E %22.15E %22.15E %22.15E\n", BAG0.val.a[1], BAG0.con.a[1],
             GAnd.val.a[1], GAnd.con.a[1]);
      printf("a2                          %22.15E %22.15E %22.15E %22.15E\n", BAG0.val.a[2], BAG0.con.a[2],
             GAnd.val.a[2], GAnd.con.a[2]);
      printf("e_i                         %22.15E                        %22.15E\n", BAG0.acvB, GAnd.acvB);
      printf("e_i^p                       %22.15E                        %22.15E\n", BAG0.acvA, GAnd.acvA);
      printf("t_b                         %22.15E                        %22.15E\n", BAG0.acvBA, GAnd.acvBA);
      printf("t_c                         %22.15E                        %22.15E\n", BAG0.bcB0, GAnd.bcB0);
      printf("t_c^p                       %22.15E                        %22.15E\n", BAG0.bcA0, GAnd.bcA0);
      printf("t_v                         %22.15E                        %22.15E\n", BAG0.bvB0, GAnd.bvB0);
      printf("t_v^p                       %22.15E                        %22.15E\n", BAG0.bvA0, GAnd.bvA0);
      printf("(t_c)^2+(t_c^p)^2           %22.15E                        %22.15E\n", Sqr(BAG0.bcB0) + Sqr(BAG0.bcA0),
             Sqr(GAnd.bcB0) + Sqr(GAnd.bcA0));
      printf("(t_v)^2+(t_v^p)^2           %22.15E                        %22.15E\n", Sqr(BAG0.bvB0) + Sqr(BAG0.bvA0),
             Sqr(GAnd.bvB0) + Sqr(GAnd.bvA0));
      printf("(t_c)^2+(t_v)^2+(t_b)^2     %22.15E                        %22.15E\n",
             Sqr(BAG0.bcB0) + Sqr(BAG0.bvB0) + Sqr(BAG0.acvBA), Sqr(GAnd.bcB0) + Sqr(GAnd.bvB0) + Sqr(GAnd.acvBA));
      printf("(t_c^p)^2+(t_v^p)^2+(t_b)^2 %22.15E                        %22.15E\n",
             Sqr(BAG0.bcA0) + Sqr(BAG0.bvA0) + Sqr(BAG0.acvBA), Sqr(GAnd.bcA0) + Sqr(GAnd.bvA0) + Sqr(GAnd.acvBA));
    }
    /* find ground-state with correct occupation by adjusting mu */
    {
      BeginTiming("Adjust mu", &Time);
      imu = 0;
      if ((iOuter == 0) || ((HSet.SymOpt) & 8))
      {
        mu[0] = Sqr(G0Occupation) * HSet.U / 2;
      }
      while (imu < imuMax)
      {
        /* Create HAnd */
        BADoubleTriDiagonalMatrixToAndersonImpurityOperator(&HAnd, GAndTrunk, GAndTrunk, HSet.U, mu[imu]);
        /* Set HAndGrd to HAnd */
        OperatorCopy(&HAndGrd, HAnd);
        /* if SymOpt (4) is set keep only H on the impurity and first bath site */
        if ((HSet.SymOpt) & 4)
        {
          printf("Restricted HAnd for the ground-state to the first bath site\n");
          HAndGrd.l[2].N = 8;
        }
        /* if singlet set is true add the (U+1) * S^2 to the operator HAndGrd */
        if ((HSet.SymOpt) & 1)
        {
          printf("Add %22.15E * S^2 to HAnd for the calculation of the ground-state wave-function\n", 1 + HSet.U);
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
          OperatorAddNormalOrder(&HAndGrd, S2, fabs(HSet.U) + 1, 0, &QDetNotOrderedFcFaBcBa);
          FreeOperator(&S2);
        }
        /* create ground-state wave-function */
        BeginTiming("Lanczos GroundState", &Time);
        LanczosGroundState(HAndGrd, &Psi, &HSet.zeroLancWavefunction, HSet.zeroLancEnergy, HSet.U, &Psi, 0);
        FreeOperator(&HAndGrd);
        EndTiming("Lanczos GroundState", &Time);
        PsiHermitianOperatorPsi(HAnd, Psi, &PsiEnergy);
        PsiHermitianOperatorSqrPsi(HAnd, Psi, &PsiQuality);
        PsiQuality /= Sqr(PsiEnergy);
        printf("Energy                  <psi|H^2|psi>/(<psi|H|psi>^2)\n");
        printf("%22.15E  %22.15E\n", PsiEnergy, PsiQuality);
        /* Start adjusting mu */
        if ((HSet.SymOpt) & 8)
        {
          imu = imuMax;
        }
        else
        {
          /* calculate cluster occupation */
          PsiRho = calloc(Psi.NFermion, sizeof(double));
          if (PsiRho == NULL)
          {
            printf("calloc failed in main for PsiRho\n Fatal Error\n");
            fflush(stdout);
            return 1;
          }
          CalculateDensity(PsiRho, Psi);
          GClOccupation[imu] = PsiRho[0] + PsiRho[1];
          free(PsiRho);
          /* Stop if occupation is close enough to occupation of G0 or if imu==imuMax-1 */
          if ((imu == imuMax - 1) || (Abs(GClOccupation[imu] - G0Occupation) < (1000 * DBL_EPSILON))
              || (imu == 0 ? 0 : (Abs(GClOccupation[imu] - GClOccupation[imu - 1]) < (1000 * DBL_EPSILON))))
          {                     /* Stop */
            /* Warning if mu not converged when imu reached imuMax */
            if (Abs(GClOccupation[imu] - G0Occupation) > (1000 * DBL_EPSILON))
            {
              printf("WARNING! stop changing imu without converging\n");
            }
            /* print occupation v.s. mu */
            printf("occupation v.s. mu\n");
            printf("  i    mu[i]                  Occupation             Difference\n");
            for (i = 0; i <= imu; i++)
            {
              printf("%3i %22.15E %22.15E %22.15E\n", i, mu[i], GClOccupation[i], GClOccupation[i] - G0Occupation);
            }
            mu[0] = mu[imu];
            imu = imuMax;
          }
          else
          {                     /* make another loop */
            FreeOperator(&HAnd);
            /* interpolate mu[imu+1] */
            if (imu == 0)
            {
              if (GClOccupation[0] > G0Occupation)
              {
                mu[1] = mu[0] - 0.1 * Abs(HSet.U);
              }
              else
              {
                mu[1] = mu[0] + 0.1 * Abs(HSet.U);
              }
            }
            else
            {
              if (Abs(GClOccupation[imu] - GClOccupation[imu - 1]) < (100 * DBL_EPSILON))
              {
                printf("Convergence creterium for occupation not reached, two mu's give the same occupation\n");
                mu[imu + 1] = mu[imu];
                imu = imuMax;
              }
              else
              {
                mu[imu + 1] =
                  mu[imu] + (G0Occupation - GClOccupation[imu]) * (mu[imu] - mu[imu - 1]) / (GClOccupation[imu] -
                                                                                             GClOccupation[imu - 1]);
              }
            }
            printf("\nFilling should be        : %22.15E.\nThe current occupation is: %22.15E.\n", G0Occupation,
                   GClOccupation[imu]);
            printf("mu adjusted from %22.15E to %22.15E\n", mu[imu], mu[imu + 1]);
            printf("Restarting Lanczos Loop %i for the %i-th time\n", iOuter, imu + 1);
            fflush(stdout);
          }
        }
        imu++;
      }
      EndTiming("Adjust mu", &Time);
      PsiHermitianOperatorPsi(HAnd, Psi, &PsiEnergy);
      PsiHermitianOperatorSqrPsi(HAnd, Psi, &PsiQuality);
      PsiQuality /= Sqr(PsiEnergy);
      PrintWaveFunctionReadableTruncated(Psi, 0.01);
      CalculatePrefactorHistogram(PsiHist, Psi, 10.0, 20);
      printf("\nHistogram of prefactors in Wave-function\n");
      for (i = 0; i < 19; i++)
      {
        printf("10^(-%2i) to 10^(-%2i) | %22.15E\n", i, i + 1, PsiHist[i]);
      }
      printf("10^(-%2i) to 0        | %22.15E\n", i, PsiHist[19]);
      printf("\n");
      printf("Energy                  <psi|H^2|psi>/(<psi|H|psi>^2)\n");
      printf("%22.15E  %22.15E\n", PsiEnergy, PsiQuality);
      PsiRho = calloc(Psi.NFermion, sizeof(double));
      if (PsiRho == NULL)
      {
        printf("calloc failed in main for PsiRho\n Fatal Error\n");
        fflush(stdout);
        return 1;
      }
      CalculateDensity(PsiRho, Psi);
      printf("Occupation \n");
      for (i = 0; i < Psi.NFermion; i++)
      {
        printf("%5i %22.15E\n", i, PsiRho[i]);
      }
      free(PsiRho);
    }
    /*                                           */
    /* create spin up and spin down PES and IPES */
    /*                                           */
    {
      CreatePESdnOperatorDMFT(&OPESdn, HAnd.NFermion);
      CreatePESupOperatorDMFT(&OPESup, HAnd.NFermion);
      CreateIPESdnOperatorDMFT(&OIPESdn, HAnd.NFermion);
      CreateIPESupOperatorDMFT(&OIPESup, HAnd.NFermion);
      /* in the last extra loop call DMFTSpectra */
      if ((iOuter == HSet.iStart + HSet.Ni - 1) && (HSet.Spectra))
      {
        if (DMFTSpec(HSet, PsiEnergy, Psi, HAnd, VPerturbation, gPerturbation, OPESdn, OPESup, OIPESdn, OIPESup))
        {
          printf("DMFTSpec failed\n");
          return 1;
        }
      }
      /* if not the last loop continue */
      BeginTiming("Lanczos Tri Diagonalization", &Time);
      OperatorPsi(OPESup, Psi, &OPsi);
      switch (HSet.Dyson)
      {
      case 1:
        GclMup.N = Max(HSet.NTriDiag, HSet.NGCluster);
        InitTriDiagonalMatrix(&GclMup);
        LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation,
                                         gPerturbation, 4, &GclMup);
        break;
      default:
        GclMup.N = HSet.NTriDiag;
        InitTriDiagonalMatrix(&GclMup);
        LanczosTriDiagonalize(HAnd, OPsi, HSet.zeroLancWavefunction, &GclMup, HSet.Gamma);
        break;
      }
      EndTiming("Lanczos Tri Diagonalization", &Time);
      BeginTiming("Lanczos Tri Diagonalization", &Time);
      OperatorPsi(OPESdn, Psi, &OPsi);
      switch (HSet.Dyson)
      {
      case 1:
        GclMdn.N = Max(HSet.NTriDiag, HSet.NGCluster);
        InitTriDiagonalMatrix(&GclMdn);
        LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation,
                                         gPerturbation, 4, &GclMdn);
        break;
      default:
        GclMdn.N = HSet.NTriDiag;
        InitTriDiagonalMatrix(&GclMdn);
        LanczosTriDiagonalize(HAnd, OPsi, HSet.zeroLancWavefunction, &GclMdn, HSet.Gamma);
        break;
      }
      EndTiming("Lanczos Tri Diagonalization", &Time);
      BeginTiming("Lanczos Tri Diagonalization", &Time);
      OperatorPsi(OIPESup, Psi, &OPsi);
      switch (HSet.Dyson)
      {
      case 1:
        GclPup.N = Max(HSet.NTriDiag, HSet.NGCluster);
        InitTriDiagonalMatrix(&GclPup);
        LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation,
                                         gPerturbation, 4, &GclPup);
        break;
      default:
        GclPup.N = HSet.NTriDiag;
        InitTriDiagonalMatrix(&GclPup);
        LanczosTriDiagonalize(HAnd, OPsi, HSet.zeroLancWavefunction, &GclPup, HSet.Gamma);
        break;
      }
      EndTiming("Lanczos Tri Diagonalization", &Time);
      BeginTiming("Lanczos Tri Diagonalization", &Time);
      OperatorPsi(OIPESdn, Psi, &OPsi);
      switch (HSet.Dyson)
      {
      case 1:
        GclPdn.N = Max(HSet.NTriDiag, HSet.NGCluster);
        InitTriDiagonalMatrix(&GclPdn);
        LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation,
                                         gPerturbation, 4, &GclPdn);
        break;
      default:
        GclPdn.N = HSet.NTriDiag;
        InitTriDiagonalMatrix(&GclPdn);
        LanczosTriDiagonalize(HAnd, OPsi, HSet.zeroLancWavefunction, &GclPdn, HSet.Gamma);
        break;
      }
      EndTiming("Lanczos Tri Diagonalization", &Time);
      /* print occupations from spectra */
      printf("\nWeight of PES and IPES Greens functions\n");
      printf("occupation:   dn                      up\n");
      printf("occupied  : %22.15E %22.15E\n", Sqr(GclMdn.b[0]), Sqr(GclMup.b[0]));
      printf("unoccupied: %22.15E %22.15E\n", Sqr(GclPdn.b[0]), Sqr(GclPup.b[0]));
    }
    printf("=== Start G math ===\n");
    fflush(stdout);
    BeginTiming("G math", &Time);
    /* G math - create Gcl as the sum of PES - IPES */
    {
      ShiftTriDiagonalMatrix(&GclMdn, -PsiEnergy);
      ShiftTriDiagonalMatrix(&GclMup, -PsiEnergy);
      ShiftTriDiagonalMatrix(&GclPdn, -PsiEnergy);
      ShiftTriDiagonalMatrix(&GclPup, -PsiEnergy);
      for (i = 1; i <= GclMdn.N; i++)
      {
        GclMdn.a[i] *= -1;
      }
      for (i = 1; i <= GclMup.N; i++)
      {
        GclMup.a[i] *= -1;
      }
      strcpy(GclPup.Name, "Gcluster-Plus-Up");
      strcpy(GclPdn.Name, "Gcluster-Plus-Dn");
      strcpy(GclMup.Name, "Gcluster-Min-Up");
      strcpy(GclMdn.Name, "Gcluster-Min-Dn");
      printf("Start SumGreensfunctionsTriDiagonalRepresentation\n");
      fflush(stdout);
      TriDiagonalMatrixToListOfPoles(GclPup, &GclList);
      TriDiagonalMatrixToListOfPoles(GclPdn, &GListTemp);
      AddListOfPoles(&GclList, GListTemp);
      FreeListOfPoles(&GListTemp);
      TriDiagonalMatrixToListOfPoles(GclMup, &GListTemp);
      AddListOfPoles(&GclList, GListTemp);
      FreeListOfPoles(&GListTemp);
      TriDiagonalMatrixToListOfPoles(GclMdn, &GListTemp);
      AddListOfPoles(&GclList, GListTemp);
      FreeListOfPoles(&GListTemp);
      ScaleListOfPoles(&GclList, 0.5);
      ListOfPolesDetermineOccupation(GclList, Ef, &G0OccupationReal);
      printf("\nOccupation of Cluster Greens function for Ef=%22.15E: %22.15E\n", Ef, G0OccupationReal);
      ListOfPolesDetermineOccupation(GclList, mu[0], &G0OccupationReal);
      printf("\nOccupation of Cluster Greens function for Ef=%22.15E: %22.15E\n\n", mu[0], G0OccupationReal);
      /* symmetrize */
      if ((HSet.SymOpt) & 2)
      {
        SymmetrizeListOfPoles(&GclList);
      }
      printf("Cluster Green function created of size %i ", GclList.N);
      ReduceListOfPolesNumberFermi(&GclList, HSet.NGCluster, Ef);
      printf("reduced to %i\n", GclList.N);
      ListOfPolesToTridiagonalMatrix(GclList, &Gcl, GclList.N);
      strcpy(Gcl.Name, "Gcluster");
      Gcl.b[0] = 1.0;
      /* make Gcl a physical solution of the full lattice model */
      /* Create self energy (SDMFT) of G0 and Gcl SDMFT = - Gcl^-1 + G0^-1 */
      /* Create Gclprime from SDMFT Gcl = (G0^-1 - SDMFT)^-1 */
      /* G math - create Cluster Self energy From Gcl and GBath */
    }
    /* G math - create Self energy */
    {
      BeginTiming("G math Selfenergy", &Time);
      printf("Start SelfEnergyTriDiagonalRepresentation\n");
      fflush(stdout);
      TriDiagonalMatrixToAndersonMatrix(Gcl, &T, &GclAnd);
      FreeMatrix(&T);
      TriDiagonalMatrixToAndersonMatrix(GBath, &T, &GBathAnd);
      FreeMatrix(&T);
      SelfEnergyListOfPolesFromAndersonRepresentation(GBathAnd, GclAnd, &SigmaPoles);
      /*
         GaussianBroadenListOfPoles(&SigmaPoles, 0.1, 10);
         ReduceListOfPolesNumber(&SigmaPoles, SigmaPoles.N/5);
       */
      printf("Cluster self energy created of size %i with shift %22.15E\n", SigmaPoles.N, SigmaPoles.a0);
      strcpy(SigmaPoles.Name, "Sigma");
      /* copy sigma to tridiagonal format */
      ListOfPolesToTridiagonalMatrix(SigmaPoles, &Sigma, SigmaPoles.N);
      printf("Cluster self energy created in tridiagonal form of size %i with weight %22.15E and shift %22.15E\n",
             Sigma.N, Sigma.b[0], Sigma.a[0]);
      EndTiming("G math Selfenergy", &Time);
    }
    /* G math - create GBath */
    {
      BeginTiming("G math GBath", &Time);
      BeginTiming("G math GBath AndListAndRepresentation", &Time);
      printf("Start GBathAndListAndRepresentation\n");
      fflush(stdout);
      printf("mu=%22.15E\n", mu[0]);
      GBathAndTriAndRepresentation(G0And, mu[0], Sigma, &GBathNewAnd);
      EndTiming("G math GBath AndListAndRepresentation", &Time);
      printf("New bath Green function created of size %i ", GBathNewAnd.N);
      strcpy(GBathNewAnd.Name, "Gbath-New");
      BeginTiming("G math GBath ReducePolesAndersonMatrixNumber", &Time);
      ReducePolesAndersonMatrixNumberFermi(&GBathNewAnd, HSet.NGBath, Ef);
      EndTiming("G math GBath ReducePolesAndersonMatrixNumber", &Time);
      printf("reduced to %i\n", GBathNewAnd.N);
      BeginTiming("G math GBath AndersonMatrixToTriDiagonalMatrix", &Time);
      AndersonMatrixToTriDiagonalMatrix(GBathNewAnd, &T, &GBathNew, GBathNewAnd.N);
      FreeMatrix(&T);
      EndTiming("G math GBath AndersonMatrixToTriDiagonalMatrix", &Time);
      /* symmetrize */
      if ((HSet.SymOpt) & 2)
      {
        SymmetrizeTridiagonalMatrix(&GBathNew);
      }
      EndTiming("G math GBath", &Time);
    }
    /* calculate GBath for the next loop by mixing GBath with GBathNew */
    {
      printf("Start mixing GBath with GBathNew\n");
      TriDiagonalMatrixCopy(&GBathOld, GBath);
      FreeTriDiagonalMatrix(&GBath);
      GBathOld.N = Min(GBathOld.N, HSet.NGBath);
      strcpy(GBathOld.Name, "Gbath-Old");
      GBathOld.b[0] = sqrt(HSet.Mix);
      GBathNew.b[0] = sqrt(1 - HSet.Mix);
      SumGreensfunctionsTriDiagonalRepresentation(GBathNew, GBathOld, &GBath);
      GBath.a[0] = HSet.Mix * GBathOld.a[0] + (1 - HSet.Mix) * GBathNew.a[0];
      strcpy(GBath.Name, "Gbath");
      GBathOld.b[0] = 1;
      GBathNew.b[0] = 1;
      if ((HSet.SymOpt) & 16)
      {
        printf("Enforce GBath to fullfill the knonw sumrules\n");
        printf("Not yet implemented\n");
      }
    }
    EndTiming("G math", &Time);
    /* create spectra */
    {
      printf("Start create spectra\n");
      fflush(stdout);
      Spectra[0].N = HSet.NE;
      Spectra[0].Emin = HSet.Emin;
      Spectra[0].Emax = HSet.Emax;
      Spectra[0].dE = (HSet.Emax - HSet.Emin) / HSet.NE;
      Spectra[0].EminPole = HSet.Emin;
      Spectra[0].EmaxPole = HSet.Emax;
      Spectra[0].Gamma = HSet.Gamma;
      for (i = 0; i < 10; i++)
      {
        SpectrumAllocate(Spectra[0], &Spectra[i]);
      }
      TriDiagonalMatrixToSpectrum(0, GBath, &(Spectra[0])); /*  3 */
      TriDiagonalMatrixToSpectrum(0, GBathOld, &(Spectra[1]));  /*  5 */
      TriDiagonalMatrixToSpectrum(0, GBathNew, &(Spectra[2]));  /*  7 */
      BADoubleTriDiagonalMatrixToSpectrum(0, GAnd, &(Spectra[3]));  /*  9 */
      TriDiagonalMatrixToSpectrum(0, G0, &(Spectra[4]));  /* 11 */
      TriDiagonalMatrixToSpectrum(0, Sigma, &(Spectra[5])); /* 13 */
      AndersonImpurityMatrixToSpectrum(0, GBathNewAnd, &(Spectra[6]));  /* 15 */
      TriDiagonalMatrixToSpectrum(0, GclPup, &(Spectra[7]));  /* 17 */
      TriDiagonalMatrixToSpectrum(0, GclMup, &(Spectra[8]));  /* 19 */
      TriDiagonalMatrixToSpectrum(0, Gcl, &(Spectra[9])); /* 21 */
    }
    /* write to disk */
    {
      printf("Start write to dist\n");
      fflush(stdout);
      if (iOuter == 0)
      {
        sprintf(GBathFileName, "%s001", HSet.GRootFileName);
        sprintf(GClFileName, "Cl%s001", HSet.GRootFileName);
        sprintf(GSigmaFileName, "Sigma%s001", HSet.GRootFileName);
        sprintf(PsiFileName, "%s001", HSet.WaveFunctionRootFileName);
        sprintf(SpecFileName, "%s001", HSet.SpecRootFileName);
      }
      else if (iOuter < 9)
      {
        sprintf(GBathFileName, "%s00%1i", HSet.GRootFileName, iOuter + 1);
        sprintf(GClFileName, "Cl%s00%1i", HSet.GRootFileName, iOuter + 1);
        sprintf(GSigmaFileName, "Sigma%s00%1i", HSet.GRootFileName, iOuter + 1);
        sprintf(PsiFileName, "%s00%1i", HSet.WaveFunctionRootFileName, iOuter + 1);
        sprintf(SpecFileName, "%s00%1i", HSet.SpecRootFileName, iOuter + 1);
      }
      else if (iOuter < 99)
      {
        sprintf(GBathFileName, "%s0%2i", HSet.GRootFileName, iOuter + 1);
        sprintf(GClFileName, "Cl%s0%2i", HSet.GRootFileName, iOuter + 1);
        sprintf(GSigmaFileName, "Sigma%s0%2i", HSet.GRootFileName, iOuter + 1);
        sprintf(PsiFileName, "%s0%2i", HSet.WaveFunctionRootFileName, iOuter + 1);
        sprintf(SpecFileName, "%s0%2i", HSet.SpecRootFileName, iOuter + 1);
      }
      else if (iOuter < 999)
      {
        sprintf(GBathFileName, "%s%3i", HSet.GRootFileName, iOuter + 1);
        sprintf(GClFileName, "Cl%s%3i", HSet.GRootFileName, iOuter + 1);
        sprintf(GSigmaFileName, "Sigma%s%3i", HSet.GRootFileName, iOuter + 1);
        sprintf(PsiFileName, "%s%3i", HSet.WaveFunctionRootFileName, iOuter + 1);
        sprintf(SpecFileName, "%s%3i", HSet.SpecRootFileName, iOuter + 1);
      }
      FilePrintSpectra(Spectra, 10, SpecFileName, "w");
      FilePrintPlainTriDiagonalMatrix(GBath, GBathFileName, "w");
      FilePrintPlainTriDiagonalMatrix(Gcl, GClFileName, "w");
      FilePrintPlainTriDiagonalMatrix(Sigma, GSigmaFileName, "w");
      FilePrintWaveFunctionPlainBinary(Psi, PsiFileName, "w");
    }
    /* free variables */
    {
      printf("Start free variables\n");
      fflush(stdout);
      FreeTriDiagonalMatrix(&GclMup);
      FreeTriDiagonalMatrix(&GclMdn);
      FreeTriDiagonalMatrix(&GclPup);
      FreeTriDiagonalMatrix(&GclPdn);
      FreeTriDiagonalMatrix(&Gcl);
      FreeListOfPoles(&SigmaPoles);
      FreeTriDiagonalMatrix(&Sigma);
      FreeTriDiagonalMatrix(&GBathNew);
      FreeTriDiagonalMatrix(&GBathOld);
      FreeBADoubleTriDiagonalMatrix(&GAnd);
      FreeAndersonMatrix(&GBathNewAnd);
      FreeAndersonMatrix(&GBathAnd);
      EndTiming("DMFTLoop", &Time);
      PrintTiming(Time);
      if (HSet.Dyson)
        for (i = 0; i < 4; i++)
        {
          FreeOperator(&(VPerturbation[i]));
          FreeTriDiagonalMatrix(&(gPerturbation[i]));
        }
      FreeOperator(&HAnd);
      FreeOperator(&OPESdn);
      FreeOperator(&OPESup);
      FreeOperator(&OIPESdn);
      FreeOperator(&OIPESup);
    }
  }
  EndTiming("Total", &Time);
  PrintTiming(Time);
  return 0;
}

int mainOld(int argc, const char *argv[])
{
  const uint32_t imuMax = 20, iHistMax = 20, StringLengthMax = 255;
  TimeType Time;
  DMFTSettings HSet;
  InitializeTiming(&Time);
  BeginTiming("Total", &Time);
  TimeStampAll();
  argc -= 1;
  /* read the DMFT settings */
  if (FileReadDMFTSettings(&HSet, argv[argc]))
  {
    printf("Fatal Error in main\n Failed to read DMFT Settings from file %s\n", argv[argc]);
    return 1;
  };
  PrintDMFTSettings(HSet);
  WaveFunctionType Psi, OPsi, PsiTemp;
  TriDiagonalMatrixType G0, GBath, GclPup, GclPdn, GclMup, GclMdn, Gcl, GBathNew, GBathOld, gPerturbation[4], Sigma;
  BADoubleTriDiagonalMatrixType GAnd, GAndTrunk, BAG0;
  AndersonMatrixType GBathAnd, GclAnd, G0And, GBathNewAnd;
  ListOfPoles SigmaPoles, GclList, GListTemp;
  OperatorType HAnd, HAndGrd, OPESdn, OPESup, OIPESdn, OIPESup, VPerturbation[4], Sx, Sy, Sz, S2;
  MatrixType T;
  SpectraType Spectra[10];
  char   GBathFileName[StringLengthMax], GClFileName[StringLengthMax], GSigmaFileName[StringLengthMax],
    PsiFileName[StringLengthMax], SpecFileName[StringLengthMax], Buffer[StringLengthMax];
  uint32_t iOuter, i, PsiNBathV, PsiNBathC, imu;
  uint16_t oadd[1], *SpinUp, *SpinDn;
  double Ef = 0.0, epsilonAnd = 0.00001, mu[imuMax + 1], G0Occupation, GClOccupation[imuMax], G0OccupationReal;
  double PsiEnergy, PsiQuality, PsiHist[iHistMax], *PsiRho;
  /* read G0 */
  if (FileReadPlainTriDiagonalMatrixType(&G0, HSet.G0FileName))
  {
    printf("Fatal Error in main\n Failed to read G0 from file %s\n", HSet.G0FileName);
    return 1;
  };
  strcpy(G0.Name, "G0");
  if ((HSet.SymOpt) & 2)
  {
    printf("Symmetrize G0(w) -= G0(-w)^*\n");
    SymmetrizeTridiagonalMatrix(&G0);
  }
  TriDiagonalMatrixToBADoubleTriDiagonalMatrix(G0, &BAG0, G0.N, G0.N, Ef, epsilonAnd);
  if ((HSet.SymOpt) & 2)
  {
    printf("Symmetrize G0(w) -= G0(-w)^*\n");
    SymmetrizeBADoubleTriDiagonalMatrix(&BAG0);
  }
  TriDiagonalMatrixToAndersonMatrix(G0, &T, &G0And);
  if ((HSet.SymOpt) & 2)
  {
    printf("Symmetrize G0(w) -= G0(-w)^*\n");
    SymmetrizeAndersonMatrix(&G0And);
  }
  FreeMatrix(&T);
  G0Occupation = 2 * Sqr(BAG0.val.b[0]);
  printf("\nRead G0 from file \"%s\"\n", HSet.G0FileName);
  printf("The occupation of the cluster is %22.15E electrons\n", G0Occupation);
  printf("If you don't like this, change G0\n\n");
  /* read G_Bath_iStart and psi_iStart */
  if (HSet.iStart == 0)
  {
    sprintf(GBathFileName, "%s000", HSet.GRootFileName);
    sprintf(PsiFileName, "%s000", HSet.WaveFunctionRootFileName);
  }
  else if (HSet.iStart < 10)
  {
    sprintf(GBathFileName, "%s00%1i", HSet.GRootFileName, HSet.iStart);
    sprintf(PsiFileName, "%s00%1i", HSet.WaveFunctionRootFileName, HSet.iStart);
  }
  else if (HSet.iStart < 100)
  {
    sprintf(GBathFileName, "%s0%2i", HSet.GRootFileName, HSet.iStart);
    sprintf(PsiFileName, "%s0%2i", HSet.WaveFunctionRootFileName, HSet.iStart);
  }
  else
  {
    sprintf(GBathFileName, "%s%i", HSet.GRootFileName, HSet.iStart);
    sprintf(PsiFileName, "%s%i", HSet.WaveFunctionRootFileName, HSet.iStart);
  }
  if (FileReadPlainTriDiagonalMatrixType(&GBath, GBathFileName))
  {
    printf("Fatal Error in main\n Failed to read GBath from file %s !\n", GBathFileName);
    return 1;
  };
  if (file_exists(PsiFileName))
  {
    Psi.NBitsKey = HSet.NBitsKey;
    if (FileReadWaveFunctionPlainBinary(&Psi, PsiFileName))
    {
      printf("Fatal Error in main\n Failed to read Psi from file %s !\n", PsiFileName);
      return 1;
    };
    printf("\nRead the starting wave function from the file %s\n", PsiFileName);
    sscanf(Psi.Name, "%s %i %i", Buffer, &PsiNBathC, &PsiNBathV);
    printf("Wave-function \"%s\" found has %i conduction bath sites, %i valence bath sites and %i Fermions in basis\n",
           Psi.Name, PsiNBathC, PsiNBathV, Psi.NFermion);
    if ((2 * (PsiNBathC + PsiNBathV + 2)) != Psi.NFermion)
    {
      printf("The total number of Fermions is incompatible witht the bath sizes\n");
      printf("This wave-function can not be used to start from\n");
      return 1;
    }
  }
  else
  {
    printf("\nWave-function file: '%s' does not excist\n Started from new function\n", PsiFileName);
    CreateWaveFunctionDMFTSecond(&Psi, HSet.NBath, HSet.NBath, HSet.NBitsKey);
  }
  strcpy(GBath.Name, "Gbath");
  /* do the DMFT loops */
  for (iOuter = HSet.iStart; iOuter < HSet.iStart + HSet.Ni; iOuter++)
  {
    BeginTiming("DMFTLoop", &Time);
    printf("\n\n***********************************************************************\n");
    printf("Start of loop %i from %i to %i\n", iOuter + 1, HSet.iStart, HSet.iStart + HSet.Ni);
    printf("***********************************************************************\n\n\n");
    /* Create splitted anderson matrix */
    //TriDiagonalMatrixToBADoubleTriDiagonalMatrixOccupation(GBath, &GAnd, GBath.N, GBath.N, Ef, G0Occupation/2);
    TriDiagonalMatrixToBADoubleTriDiagonalMatrix(GBath, &GAnd, GBath.N, GBath.N, Ef, epsilonAnd);
    if ((HSet.SymOpt) & 2)
    {
      printf("Symmetrize GAnd(w) -= GAnd(-w)^*\n");
      printf("Size of GAnd went from %i + %i ", GAnd.val.N, GAnd.con.N);
      SymmetrizeBADoubleTriDiagonalMatrix(&GAnd);
      printf("to %i + %i ", GAnd.val.N, GAnd.con.N);
      GAnd.val.N = Min(GAnd.val.N, HSet.NGBath);
      GAnd.con.N = Min(GAnd.con.N, HSet.NGBath);
      printf("to %i + %i\n", GAnd.val.N, GAnd.con.N);
    }
    if ((HSet.SymOpt) & 16)
    {
      printf("Enforce GAnd to fullfill the knonw sumrules\n");
      printf("Not yet implemented\n");
    }
    /* Create GAnd and GAndTrunk as a copy pointing to the same memory */
    GAndTrunk.Ef = GAnd.Ef;
    GAndTrunk.Epsilon = GAnd.Epsilon;
    GAndTrunk.a0 = GAnd.a0;
    GAndTrunk.b0 = GAnd.b0;
    GAndTrunk.acvB = GAnd.acvB;
    GAndTrunk.acvA = GAnd.acvA;
    GAndTrunk.acvBA = GAnd.acvBA;
    GAndTrunk.bvB0 = GAnd.bvB0;
    GAndTrunk.bcB0 = GAnd.bcB0;
    GAndTrunk.bvA0 = GAnd.bvA0;
    GAndTrunk.bcA0 = GAnd.bcA0;
    GAndTrunk.val.N = GAnd.val.N;
    GAndTrunk.val.a = GAnd.val.a;
    GAndTrunk.val.b = GAnd.val.b;
    GAndTrunk.con.N = GAnd.con.N;
    GAndTrunk.con.a = GAnd.con.a;
    GAndTrunk.con.b = GAnd.con.b;
    if (HSet.Dyson)
    {
      /* Create VPerturbation */
      for (i = 0; i < 4; i++)
      {
        InitOperator(&(VPerturbation[i]));
        VPerturbation[i].NFermion = 4 + 2 * (Min(-1 + GAnd.val.N, HSet.NBath) + Min(-1 + GAnd.con.N, HSet.NBath));
        VPerturbation[i].NBoson = 0;
        VPerturbation[i].NWord = Min(VPerturbation[i].NFermion, 1 + (VPerturbation[i].NFermion - 1) / 8);
        VPerturbation[i].QComplex = 0;
        InitOperatorFixedLength(&(VPerturbation[i]), 1, 1);
      }
      if (GAnd.val.N > (HSet.NBath + 1))
      {
        oadd[0] = CrAnMask + 3 + HSet.NBath;
        RealOperatorAddLader(&(VPerturbation[0]), oadd, 1, GAnd.val.b[HSet.NBath + 1]);
        oadd[0] = CrAnMask + 3 + Min(-1 + GAnd.con.N, HSet.NBath) + 2 * HSet.NBath;
        RealOperatorAddLader(&(VPerturbation[1]), oadd, 1, GAnd.val.b[HSet.NBath + 1]);
      }
      if (GAnd.con.N > HSet.NBath + 1)
      {
        oadd[0] = 3 + Min(-1 + GAnd.val.N, HSet.NBath) + HSet.NBath;
        RealOperatorAddLader(&(VPerturbation[2]), oadd, 1, GAnd.con.b[HSet.NBath + 1]);
        oadd[0] = 3 + 2 * Min(-1 + GAnd.val.N, HSet.NBath) + 2 * HSet.NBath;
        RealOperatorAddLader(&(VPerturbation[3]), oadd, 1, GAnd.con.b[HSet.NBath + 1]);
      }
      /* Create gPerturbation */
      if (GAnd.val.N > HSet.NBath + 1)
      {
        gPerturbation[0].N = GAnd.val.N - HSet.NBath - 1;
        gPerturbation[1].N = GAnd.val.N - HSet.NBath - 1;
        InitTriDiagonalMatrix(&(gPerturbation[0]));
        InitTriDiagonalMatrix(&(gPerturbation[1]));
        for (i = 1; i < GAnd.val.N - HSet.NBath - 1; i++)
        {
          gPerturbation[0].a[i] = -GAnd.val.a[i + HSet.NBath + 1];
          gPerturbation[0].b[i] = GAnd.val.b[i + HSet.NBath + 1];
          gPerturbation[1].a[i] = -GAnd.val.a[i + HSet.NBath + 1];
          gPerturbation[1].b[i] = GAnd.val.b[i + HSet.NBath + 1];
        }
        gPerturbation[0].a[GAnd.val.N - HSet.NBath - 1] = -GAnd.val.a[GAnd.val.N];
        gPerturbation[1].a[GAnd.val.N - HSet.NBath - 1] = -GAnd.val.a[GAnd.val.N];
      }
      else
      {
        gPerturbation[0].N = 1;
        gPerturbation[1].N = 1;
        InitTriDiagonalMatrix(&(gPerturbation[0]));
        InitTriDiagonalMatrix(&(gPerturbation[1]));
        gPerturbation[0].a[1] = 0;
        gPerturbation[1].a[1] = 0;
      }
      if (GAnd.con.N > HSet.NBath + 1)
      {
        gPerturbation[2].N = GAnd.con.N - HSet.NBath - 1;
        gPerturbation[3].N = GAnd.con.N - HSet.NBath - 1;
        InitTriDiagonalMatrix(&(gPerturbation[2]));
        InitTriDiagonalMatrix(&(gPerturbation[3]));
        for (i = 1; i < GAnd.con.N - HSet.NBath - 1; i++)
        {
          gPerturbation[2].a[i] = GAnd.con.a[i + HSet.NBath + 1];
          gPerturbation[2].b[i] = GAnd.con.b[i + HSet.NBath + 1];
          gPerturbation[3].a[i] = GAnd.con.a[i + HSet.NBath + 1];
          gPerturbation[3].b[i] = GAnd.con.b[i + HSet.NBath + 1];
        }
        gPerturbation[2].a[GAnd.con.N - HSet.NBath - 1] = GAnd.con.a[GAnd.con.N];
        gPerturbation[3].a[GAnd.con.N - HSet.NBath - 1] = GAnd.con.a[GAnd.con.N];
      }
      else
      {
        gPerturbation[2].N = 1;
        gPerturbation[3].N = 1;
        InitTriDiagonalMatrix(&(gPerturbation[2]));
        InitTriDiagonalMatrix(&(gPerturbation[3]));
        gPerturbation[2].a[1] = 0;
        gPerturbation[3].a[1] = 0;
      }
    }
    /* Create GAndTrunk */
    if (GAnd.con.N < HSet.NBath + 1)
    {
      printf("Warning: GAnd.con.N(%i)<HSet.NBath+1(%i+1)\n", GAnd.con.N, HSet.NBath);
    }
    else
    {
      GAndTrunk.con.N = HSet.NBath + 1;
    }
    if (GAnd.val.N < HSet.NBath + 1)
    {
      printf("Warning: GAnd.val.N(%i)<HSet.NBath+1(%i+1)\n", GAnd.val.N, HSet.NBath);
    }
    else
    {
      GAndTrunk.val.N = HSet.NBath + 1;
    }
    /* test if psi has correct dimensions and correct if needed */
    sscanf(Psi.Name, "%s %i %i", Buffer, &PsiNBathC, &PsiNBathV);
    if (PsiNBathC != GAndTrunk.con.N - 1)
    {
      printf("\nBefore diagonalizing:\nThe wavefunction has %i conduction bath sites\nThe Anderson Hamiltonian %i\n",
             PsiNBathC, GAndTrunk.con.N - 1);
      if (PsiNBathC > GAndTrunk.con.N - 1)
      {
        printf("Removing %i unoccupied fermions from the wave-function\n", PsiNBathC - GAndTrunk.con.N + 1);
        for (i = 0; i < PsiNBathC - GAndTrunk.con.N + 1; i++)
        {
          WaveFunctionRemoveFermion(Psi, &PsiTemp, Psi.NFermion - 1);
          WaveFunctionFree(&Psi);
          WaveFunctionRemoveFermion(PsiTemp, &Psi, 4 + PsiNBathV + PsiNBathC - i - 1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      else
      {
        printf("Adding %i unoccupied fermions to the wave-function\n", GAndTrunk.con.N - 1 - PsiNBathC);
        for (i = 0; i < GAndTrunk.con.N - 1 - PsiNBathC; i++)
        {
          WaveFunctionAddFermion(Psi, &PsiTemp, Psi.NFermion, 0);
          WaveFunctionFree(&Psi);
          WaveFunctionAddFermion(PsiTemp, &Psi, 4 + PsiNBathV + PsiNBathC + i, 0);
          WaveFunctionFree(&PsiTemp);
        }
      }
    }
    if (PsiNBathV != GAndTrunk.val.N - 1)
    {
      printf("\nBefore diagonalizing:\nThe wavefunction has %i valence bath sites\nThe Anderson Hamiltonian %i\n",
             PsiNBathV, GAndTrunk.val.N - 1);
      if (PsiNBathV > GAndTrunk.val.N - 1)
      {
        printf("Removing %i occupied fermions from the wave-function\n", PsiNBathV - GAndTrunk.val.N + 1);
        for (i = 0; i < PsiNBathV - GAndTrunk.val.N + 1; i++)
        {
          WaveFunctionRemoveFermion(Psi, &PsiTemp, Psi.NFermion - GAndTrunk.con.N);
          WaveFunctionFree(&Psi);
          WaveFunctionRemoveFermion(PsiTemp, &Psi, 4 + PsiNBathV - i - 1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      else
      {
        printf("Adding %i occupied fermions to the wave-function\n", GAndTrunk.val.N - 1 - PsiNBathV);
        for (i = 0; i < GAndTrunk.val.N - 1 - PsiNBathV; i++)
        {
          WaveFunctionAddFermion(Psi, &PsiTemp, Psi.NFermion - GAndTrunk.val.N + 1, 1);
          WaveFunctionFree(&Psi);
          WaveFunctionAddFermion(PsiTemp, &Psi, 4 + PsiNBathV + i, 1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      printf("Corrected basis\n");
      PrintWaveFunctionReadableTruncated(Psi, 0.1);
    }
    sprintf(Psi.Name, "DMFT %i %i", GAndTrunk.con.N - 1, GAndTrunk.val.N - 1);
    /* Print some info on G0 and GAnd */
    printf
      ("                               G0v                    G0c                    GAndv                  Gandc\n");
    printf("b0                          %22.15E %22.15E %22.15E %22.15E\n", BAG0.val.b[0], BAG0.con.b[0], GAnd.val.b[0],
           GAnd.con.b[0]);
    printf("b1                          %22.15E %22.15E %22.15E %22.15E\n", BAG0.val.b[1], BAG0.con.b[1], GAnd.val.b[1],
           GAnd.con.b[1]);
    printf("b2                          %22.15E %22.15E %22.15E %22.15E\n", BAG0.val.b[2], BAG0.con.b[2], GAnd.val.b[2],
           GAnd.con.b[2]);
    printf("a1                          %22.15E %22.15E %22.15E %22.15E\n", BAG0.val.a[1], BAG0.con.a[1], GAnd.val.a[1],
           GAnd.con.a[1]);
    printf("a2                          %22.15E %22.15E %22.15E %22.15E\n", BAG0.val.a[2], BAG0.con.a[2], GAnd.val.a[2],
           GAnd.con.a[2]);
    printf("e_i                         %22.15E                        %22.15E\n", BAG0.acvB, GAnd.acvB);
    printf("e_i^p                       %22.15E                        %22.15E\n", BAG0.acvA, GAnd.acvA);
    printf("t_b                         %22.15E                        %22.15E\n", BAG0.acvBA, GAnd.acvBA);
    printf("t_c                         %22.15E                        %22.15E\n", BAG0.bcB0, GAnd.bcB0);
    printf("t_c^p                       %22.15E                        %22.15E\n", BAG0.bcA0, GAnd.bcA0);
    printf("t_v                         %22.15E                        %22.15E\n", BAG0.bvB0, GAnd.bvA0);
    printf("t_v^p                       %22.15E                        %22.15E\n", BAG0.bvA0, GAnd.bvA0);
    printf("(t_c)^2+(t_c^p)^2           %22.15E                        %22.15E\n", Sqr(BAG0.bcB0) + Sqr(BAG0.bcA0),
           Sqr(GAnd.bcB0) + Sqr(GAnd.bcA0));
    printf("(t_v)^2+(t_v^p)^2           %22.15E                        %22.15E\n", Sqr(BAG0.bvB0) + Sqr(BAG0.bvA0),
           Sqr(GAnd.bvB0) + Sqr(GAnd.bvA0));
    printf("(t_c)^2+(t_v)^2+(t_b)^2     %22.15E                        %22.15E\n",
           Sqr(BAG0.bcB0) + Sqr(BAG0.bvB0) + Sqr(BAG0.acvBA), Sqr(GAnd.bcB0) + Sqr(GAnd.bvB0) + Sqr(GAnd.acvBA));
    printf("(t_c^p)^2+(t_v^p)^2+(t_b)^2 %22.15E                        %22.15E\n",
           Sqr(BAG0.bcA0) + Sqr(BAG0.bvA0) + Sqr(BAG0.acvBA), Sqr(GAnd.bcA0) + Sqr(GAnd.bvA0) + Sqr(GAnd.acvBA));
    /* find ground-state with correct occupation by adjusting mu */
    BeginTiming("Adjust mu", &Time);
    imu = 0;
    if ((iOuter == 0) || ((HSet.SymOpt) & 8))
    {
      mu[0] = Sqr(G0Occupation) * HSet.U / 2;
    }
    while (imu < imuMax)
    {
      /* Create HAnd */
      BADoubleTriDiagonalMatrixToAndersonImpurityOperator(&HAnd, GAndTrunk, GAndTrunk, HSet.U, mu[imu]);
      /* Set HAndGrd to HAnd */
      OperatorCopy(&HAndGrd, HAnd);
      /* if SymOpt (4) is set keep only H on the impurity and first bath site */
      if ((HSet.SymOpt) & 4)
      {
        printf("Restricted HAnd for the ground-state to the first bath site\n");
        HAndGrd.l[2].N = 8;
      }
      /* if singlet set is true add the (U+1) * S^2 to the operator HAndGrd */
      if ((HSet.SymOpt) & 1)
      {
        printf("Add %22.15E * S^2 to HAnd for the calculation of the ground-state wave-function\n", 1 + HSet.U);
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
        OperatorAddNormalOrder(&HAndGrd, S2, fabs(HSet.U) + 1, 0, &QDetNotOrderedFcFaBcBa);
        FreeOperator(&S2);
      }
      /* create ground-state wave-function */
      BeginTiming("Lanczos GroundState", &Time);
      LanczosGroundState(HAndGrd, &Psi, &HSet.zeroLancWavefunction, HSet.zeroLancEnergy, HSet.U, &Psi, 0);
      FreeOperator(&HAndGrd);
      EndTiming("Lanczos GroundState", &Time);
      PsiHermitianOperatorPsi(HAnd, Psi, &PsiEnergy);
      PsiHermitianOperatorSqrPsi(HAnd, Psi, &PsiQuality);
      PsiQuality /= Sqr(PsiEnergy);
      printf("Energy                  <psi|H^2|psi>/(<psi|H|psi>^2)\n");
      printf("%22.15E  %22.15E\n", PsiEnergy, PsiQuality);
      /* Start adjusting mu */
      if ((HSet.SymOpt) & 8)
      {
        imu = imuMax;
      }
      else
      {
        /* calculate cluster occupation */
        PsiRho = calloc(Psi.NFermion, sizeof(double));
        if (PsiRho == NULL)
        {
          printf("calloc failed in main for PsiRho\n Fatal Error\n");
          fflush(stdout);
          return 1;
        }
        CalculateDensity(PsiRho, Psi);
        GClOccupation[imu] = PsiRho[0] + PsiRho[1];
        free(PsiRho);
        /* Stop if occupation is close enough to occupation of G0 or if imu==imuMax-1 */
        if ((imu == imuMax - 1) || (Abs(GClOccupation[imu] - G0Occupation) < (1000 * DBL_EPSILON))
            || (imu == 0 ? 0 : (Abs(GClOccupation[imu] - GClOccupation[imu - 1]) < (1000 * DBL_EPSILON))))
        {                       /* Stop */
          /* Warning if mu not converged when imu reached imuMax */
          if (Abs(GClOccupation[imu] - G0Occupation) > (1000 * DBL_EPSILON))
          {
            printf("WARNING! stop changing imu without converging\n");
          }
          /* print occupation v.s. mu */
          printf("occupation v.s. mu\n");
          printf("  i    mu[i]                  Occupation             Difference\n");
          for (i = 0; i <= imu; i++)
          {
            printf("%3i %22.15E %22.15E %22.15E\n", i, mu[i], GClOccupation[i], GClOccupation[i] - G0Occupation);
          }
          mu[0] = mu[imu];
          imu = imuMax;
        }
        else
        {                       /* make another loop */
          FreeOperator(&HAnd);
          /* interpolate mu[imu+1] */
          if (imu == 0)
          {
            if (GClOccupation[0] > G0Occupation)
            {
              mu[1] = mu[0] - 0.1 * Abs(HSet.U);
            }
            else
            {
              mu[1] = mu[0] + 0.1 * Abs(HSet.U);
            }
          }
          else
          {
            if (Abs(GClOccupation[imu] - GClOccupation[imu - 1]) < (100 * DBL_EPSILON))
            {
              printf("Convergence creterium for occupation not reached, two mu's give the same occupation\n");
              mu[imu + 1] = mu[imu];
              imu = imuMax;
            }
            else
            {
              mu[imu + 1] =
                mu[imu] + (G0Occupation - GClOccupation[imu]) * (mu[imu] - mu[imu - 1]) / (GClOccupation[imu] -
                                                                                           GClOccupation[imu - 1]);
            }
          }
          printf("\nFilling should be        : %22.15E.\nThe current occupation is: %22.15E.\n", G0Occupation,
                 GClOccupation[imu]);
          printf("mu adjusted from %22.15E to %22.15E\n", mu[imu], mu[imu + 1]);
          printf("Restarting Lanczos Loop %i for the %i-th time\n", iOuter, imu + 1);
          fflush(stdout);
        }
      }
      imu++;
    }
    EndTiming("Adjust mu", &Time);
    PsiHermitianOperatorPsi(HAnd, Psi, &PsiEnergy);
    PsiHermitianOperatorSqrPsi(HAnd, Psi, &PsiQuality);
    PsiQuality /= Sqr(PsiEnergy);
    PrintWaveFunctionReadableTruncated(Psi, 0.01);
    CalculatePrefactorHistogram(PsiHist, Psi, 10.0, 20);
    printf("\nHistogram of prefactors in Wave-function\n");
    for (i = 0; i < 19; i++)
    {
      printf("10^(-%2i) to 10^(-%2i) | %22.15E\n", i, i + 1, PsiHist[i]);
    }
    printf("10^(-%2i) to 0        | %22.15E\n", i, PsiHist[19]);
    printf("\n");
    printf("Energy                  <psi|H^2|psi>/(<psi|H|psi>^2)\n");
    printf("%22.15E  %22.15E\n", PsiEnergy, PsiQuality);
    PsiRho = calloc(Psi.NFermion, sizeof(double));
    if (PsiRho == NULL)
    {
      printf("calloc failed in main for PsiRho\n Fatal Error\n");
      fflush(stdout);
      return 1;
    }
    CalculateDensity(PsiRho, Psi);
    printf("Occupation \n");
    for (i = 0; i < Psi.NFermion; i++)
    {
      printf("%5i %22.15E\n", i, PsiRho[i]);
    }
    free(PsiRho);
    /*                                           */
    /* create spin up and spin down PES and IPES */
    /*                                           */
    CreatePESdnOperatorDMFT(&OPESdn, HAnd.NFermion);
    CreatePESupOperatorDMFT(&OPESup, HAnd.NFermion);
    CreateIPESdnOperatorDMFT(&OIPESdn, HAnd.NFermion);
    CreateIPESupOperatorDMFT(&OIPESup, HAnd.NFermion);
    /* in the last extra loop call DMFTSpectra */
    if ((iOuter == HSet.iStart + HSet.Ni - 1) && (HSet.Spectra))
    {
      if (DMFTSpec(HSet, PsiEnergy, Psi, HAnd, VPerturbation, gPerturbation, OPESdn, OPESup, OIPESdn, OIPESup))
      {
        printf("DMFTSpec failed\n");
        return 1;
      }
    }
    /* if not the last loop continue */
    BeginTiming("Lanczos Tri Diagonalization", &Time);
    OperatorPsi(OPESup, Psi, &OPsi);
    switch (HSet.Dyson)
    {
    case 1:
      GclMup.N = Max(HSet.NTriDiag, HSet.NGCluster);
      InitTriDiagonalMatrix(&GclMup);
      LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation,
                                       gPerturbation, 4, &GclMup);
      break;
    default:
      GclMup.N = HSet.NTriDiag;
      InitTriDiagonalMatrix(&GclMup);
      LanczosTriDiagonalize(HAnd, OPsi, HSet.zeroLancWavefunction, &GclMup, HSet.Gamma);
      break;
    }
    EndTiming("Lanczos Tri Diagonalization", &Time);
    BeginTiming("Lanczos Tri Diagonalization", &Time);
    OperatorPsi(OPESdn, Psi, &OPsi);
    switch (HSet.Dyson)
    {
    case 1:
      GclMdn.N = Max(HSet.NTriDiag, HSet.NGCluster);
      InitTriDiagonalMatrix(&GclMdn);
      LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation,
                                       gPerturbation, 4, &GclMdn);
      break;
    default:
      GclMdn.N = HSet.NTriDiag;
      InitTriDiagonalMatrix(&GclMdn);
      LanczosTriDiagonalize(HAnd, OPsi, HSet.zeroLancWavefunction, &GclMdn, HSet.Gamma);
      break;
    }
    EndTiming("Lanczos Tri Diagonalization", &Time);
    BeginTiming("Lanczos Tri Diagonalization", &Time);
    OperatorPsi(OIPESup, Psi, &OPsi);
    switch (HSet.Dyson)
    {
    case 1:
      GclPup.N = Max(HSet.NTriDiag, HSet.NGCluster);
      InitTriDiagonalMatrix(&GclPup);
      LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation,
                                       gPerturbation, 4, &GclPup);
      break;
    default:
      GclPup.N = HSet.NTriDiag;
      InitTriDiagonalMatrix(&GclPup);
      LanczosTriDiagonalize(HAnd, OPsi, HSet.zeroLancWavefunction, &GclPup, HSet.Gamma);
      break;
    }
    EndTiming("Lanczos Tri Diagonalization", &Time);
    BeginTiming("Lanczos Tri Diagonalization", &Time);
    OperatorPsi(OIPESdn, Psi, &OPsi);
    switch (HSet.Dyson)
    {
    case 1:
      GclPdn.N = Max(HSet.NTriDiag, HSet.NGCluster);
      InitTriDiagonalMatrix(&GclPdn);
      LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation,
                                       gPerturbation, 4, &GclPdn);
      break;
    default:
      GclPdn.N = HSet.NTriDiag;
      InitTriDiagonalMatrix(&GclPdn);
      LanczosTriDiagonalize(HAnd, OPsi, HSet.zeroLancWavefunction, &GclPdn, HSet.Gamma);
      break;
    }
    EndTiming("Lanczos Tri Diagonalization", &Time);
    /* print occupations from spectra */
    printf("/nWeight of PES and IPES Greens functions\n");
    printf("occupation:   dn                      up\n");
    printf("occupied  : %22.15E %22.15E\n", Sqr(GclMdn.b[0]), Sqr(GclMup.b[0]));
    printf("unoccupied: %22.15E %22.15E\n", Sqr(GclPdn.b[0]), Sqr(GclPup.b[0]));
    /* G math - create Gcl as the sum of PES - IPES */
    printf("=== Start G math ===\n");
    fflush(stdout);
    BeginTiming("G math", &Time);
    ShiftTriDiagonalMatrix(&GclMdn, -PsiEnergy);
    ShiftTriDiagonalMatrix(&GclMup, -PsiEnergy);
    ShiftTriDiagonalMatrix(&GclPdn, -PsiEnergy);
    ShiftTriDiagonalMatrix(&GclPup, -PsiEnergy);
    for (i = 1; i <= GclMdn.N; i++)
    {
      GclMdn.a[i] *= -1;
    }
    for (i = 1; i <= GclMup.N; i++)
    {
      GclMup.a[i] *= -1;
    }
    strcpy(GclPup.Name, "Gcluster-Plus-Up");
    strcpy(GclPdn.Name, "Gcluster-Plus-Dn");
    strcpy(GclMup.Name, "Gcluster-Min-Up");
    strcpy(GclMdn.Name, "Gcluster-Min-Dn");
    printf("Start SumGreensfunctionsTriDiagonalRepresentation\n");
    fflush(stdout);
    TriDiagonalMatrixToListOfPoles(GclPup, &GclList);
    TriDiagonalMatrixToListOfPoles(GclPdn, &GListTemp);
    AddListOfPoles(&GclList, GListTemp);
    FreeListOfPoles(&GListTemp);
    TriDiagonalMatrixToListOfPoles(GclMup, &GListTemp);
    AddListOfPoles(&GclList, GListTemp);
    FreeListOfPoles(&GListTemp);
    TriDiagonalMatrixToListOfPoles(GclMdn, &GListTemp);
    AddListOfPoles(&GclList, GListTemp);
    FreeListOfPoles(&GListTemp);
    ScaleListOfPoles(&GclList, 0.5);

    ListOfPolesDetermineOccupation(GclList, Ef, &G0OccupationReal);
    printf("\nOccupation of Cluster Greens function: %22.15E\n\n", G0OccupationReal);

    /* symmetrize */
    if ((HSet.SymOpt) & 2)
    {
      SymmetrizeListOfPoles(&GclList);
    }
    printf("Cluster Green function created of size %i ", GclList.N);
    ReduceListOfPolesNumberFermi(&GclList, HSet.NGCluster, Ef);
    printf("reduced to %i\n", GclList.N);
    ListOfPolesToTridiagonalMatrix(GclList, &Gcl, GclList.N);
    strcpy(Gcl.Name, "Gcluster");
    Gcl.b[0] = 1.0;
    /* make Gcl a physical solution of the full lattice model */
    /* Create self energy (SDMFT) of G0 and Gcl SDMFT = - Gcl^-1 + G0^-1 */
    /* Create Gclprime from SDMFT Gcl = (G0^-1 - SDMFT)^-1 */
    /* G math - create Cluster Self energy From Gcl and GBath */
    BeginTiming("G math Selfenergy", &Time);
    printf("Start SelfEnergyTriDiagonalRepresentation\n");
    fflush(stdout);
    TriDiagonalMatrixToAndersonMatrix(Gcl, &T, &GclAnd);
    FreeMatrix(&T);
    TriDiagonalMatrixToAndersonMatrix(GBath, &T, &GBathAnd);
    FreeMatrix(&T);
    SelfEnergyListOfPolesFromAndersonRepresentation(GBathAnd, GclAnd, &SigmaPoles);
    /*
       GaussianBroadenListOfPoles(&SigmaPoles, 0.1, 10);
       ReduceListOfPolesNumber(&SigmaPoles, SigmaPoles.N/5);
     */
    printf("Cluster self energy created of size %i with shift %22.15E\n", SigmaPoles.N, SigmaPoles.a0);
    strcpy(SigmaPoles.Name, "Sigma");
    /* copy sigma to tridiagonal format */
    ListOfPolesToTridiagonalMatrix(SigmaPoles, &Sigma, SigmaPoles.N);
    printf("Cluster self energy created in tridiagonal form of size %i with weight %22.15E and shift %22.15E\n",
           Sigma.N, Sigma.b[0], Sigma.a[0]);
    EndTiming("G math Selfenergy", &Time);
    /* G math - create GBath */
    BeginTiming("G math GBath", &Time);
    BeginTiming("G math GBath AndListAndRepresentation", &Time);
    printf("Start GBathAndListAndRepresentation\n");
    fflush(stdout);
    GBathAndTriAndRepresentation(G0And, mu[imu], Sigma, &GBathNewAnd);
    EndTiming("G math GBath AndListAndRepresentation", &Time);
    printf("New bath Green function created of size %i ", GBathNewAnd.N);
    strcpy(GBathNewAnd.Name, "Gbath-New");
    BeginTiming("G math GBath ReducePolesAndersonMatrixNumber", &Time);
    ReducePolesAndersonMatrixNumberFermi(&GBathNewAnd, HSet.NGBath, Ef);
    EndTiming("G math GBath ReducePolesAndersonMatrixNumber", &Time);
    printf("reduced to %i\n", GBathNewAnd.N);
    BeginTiming("G math GBath AndersonMatrixToTriDiagonalMatrix", &Time);
    AndersonMatrixToTriDiagonalMatrix(GBathNewAnd, &T, &GBathNew, GBathNewAnd.N);
    FreeMatrix(&T);
    EndTiming("G math GBath AndersonMatrixToTriDiagonalMatrix", &Time);
    /* symmetrize */
    if ((HSet.SymOpt) & 2)
    {
      SymmetrizeTridiagonalMatrix(&GBathNew);
    }
    EndTiming("G math GBath", &Time);
    /* calculate GBath for the next loop by mixing GBath with GBathNew */
    printf("Start mixing GBath with GBathNew\n");
    TriDiagonalMatrixCopy(&GBathOld, GBath);
    FreeTriDiagonalMatrix(&GBath);
    GBathOld.N = Min(GBathOld.N, HSet.NGBath);
    strcpy(GBathOld.Name, "Gbath-Old");
    GBathOld.b[0] = sqrt(HSet.Mix);
    GBathNew.b[0] = sqrt(1 - HSet.Mix);
    SumGreensfunctionsTriDiagonalRepresentation(GBathNew, GBathOld, &GBath);
    GBath.a[0] = HSet.Mix * GBathOld.a[0] + (1 - HSet.Mix) * GBathNew.a[0];
    strcpy(GBath.Name, "Gbath");
    GBathOld.b[0] = 1;
    GBathNew.b[0] = 1;
    if ((HSet.SymOpt) & 16)
    {
      printf("Enforce GBath to fullfill the knonw sumrules\n");
      printf("Not yet implemented\n");
    }
    EndTiming("G math", &Time);
    /* create spectra */
    printf("Start create spectra\n");
    fflush(stdout);
    Spectra[0].N = HSet.NE;
    Spectra[0].Emin = HSet.Emin;
    Spectra[0].Emax = HSet.Emax;
    Spectra[0].dE = (HSet.Emax - HSet.Emin) / HSet.NE;
    Spectra[0].EminPole = HSet.Emin;
    Spectra[0].EmaxPole = HSet.Emax;
    Spectra[0].Gamma = HSet.Gamma;
    for (i = 0; i < 10; i++)
    {
      SpectrumAllocate(Spectra[0], &Spectra[i]);
    }
    TriDiagonalMatrixToSpectrum(0, GBath, &(Spectra[0])); /*  3 */
    TriDiagonalMatrixToSpectrum(0, GBathOld, &(Spectra[1]));  /*  5 */
    TriDiagonalMatrixToSpectrum(0, GBathNew, &(Spectra[2]));  /*  7 */
    BADoubleTriDiagonalMatrixToSpectrum(0, GAnd, &(Spectra[3]));  /*  9 */
    TriDiagonalMatrixToSpectrum(0, G0, &(Spectra[4]));  /* 11 */
    TriDiagonalMatrixToSpectrum(0, Sigma, &(Spectra[5])); /* 13 */
    AndersonImpurityMatrixToSpectrum(0, GBathNewAnd, &(Spectra[6]));  /* 15 */
    TriDiagonalMatrixToSpectrum(0, GclPup, &(Spectra[7]));  /* 17 */
    TriDiagonalMatrixToSpectrum(0, GclMup, &(Spectra[8]));  /* 19 */
    TriDiagonalMatrixToSpectrum(0, Gcl, &(Spectra[9])); /* 21 */
    /* write to disk */
    printf("Start write to dist\n");
    fflush(stdout);
    if (iOuter == 0)
    {
      sprintf(GBathFileName, "%s001", HSet.GRootFileName);
      sprintf(GClFileName, "Cl%s001", HSet.GRootFileName);
      sprintf(GSigmaFileName, "Sigma%s001", HSet.GRootFileName);
      sprintf(PsiFileName, "%s001", HSet.WaveFunctionRootFileName);
      sprintf(SpecFileName, "%s001", HSet.SpecRootFileName);
    }
    else if (iOuter < 9)
    {
      sprintf(GBathFileName, "%s00%1i", HSet.GRootFileName, iOuter + 1);
      sprintf(GClFileName, "Cl%s00%1i", HSet.GRootFileName, iOuter + 1);
      sprintf(GSigmaFileName, "Sigma%s00%1i", HSet.GRootFileName, iOuter + 1);
      sprintf(PsiFileName, "%s00%1i", HSet.WaveFunctionRootFileName, iOuter + 1);
      sprintf(SpecFileName, "%s00%1i", HSet.SpecRootFileName, iOuter + 1);
    }
    else if (iOuter < 99)
    {
      sprintf(GBathFileName, "%s0%2i", HSet.GRootFileName, iOuter + 1);
      sprintf(GClFileName, "Cl%s0%2i", HSet.GRootFileName, iOuter + 1);
      sprintf(GSigmaFileName, "Sigma%s0%2i", HSet.GRootFileName, iOuter + 1);
      sprintf(PsiFileName, "%s0%2i", HSet.WaveFunctionRootFileName, iOuter + 1);
      sprintf(SpecFileName, "%s0%2i", HSet.SpecRootFileName, iOuter + 1);
    }
    else if (iOuter < 999)
    {
      sprintf(GBathFileName, "%s%3i", HSet.GRootFileName, iOuter + 1);
      sprintf(GClFileName, "Cl%s%3i", HSet.GRootFileName, iOuter + 1);
      sprintf(GSigmaFileName, "Sigma%s%3i", HSet.GRootFileName, iOuter + 1);
      sprintf(PsiFileName, "%s%3i", HSet.WaveFunctionRootFileName, iOuter + 1);
      sprintf(SpecFileName, "%s%3i", HSet.SpecRootFileName, iOuter + 1);
    }
    FilePrintSpectra(Spectra, 10, SpecFileName, "w");
    FilePrintPlainTriDiagonalMatrix(GBath, GBathFileName, "w");
    FilePrintPlainTriDiagonalMatrix(Gcl, GClFileName, "w");
    FilePrintPlainTriDiagonalMatrix(Sigma, GSigmaFileName, "w");
    FilePrintWaveFunctionPlainBinary(Psi, PsiFileName, "w");
    /* free variables */
    printf("Start free variables\n");
    fflush(stdout);
    FreeTriDiagonalMatrix(&GclMup);
    FreeTriDiagonalMatrix(&GclMdn);
    FreeTriDiagonalMatrix(&GclPup);
    FreeTriDiagonalMatrix(&GclPdn);
    FreeTriDiagonalMatrix(&Gcl);
    FreeListOfPoles(&SigmaPoles);
    FreeTriDiagonalMatrix(&Sigma);
    FreeTriDiagonalMatrix(&GBathNew);
    FreeTriDiagonalMatrix(&GBathOld);
    FreeBADoubleTriDiagonalMatrix(&GAnd);
    FreeAndersonMatrix(&GBathNewAnd);
    FreeAndersonMatrix(&GBathAnd);
    EndTiming("DMFTLoop", &Time);
    PrintTiming(Time);
    if (HSet.Dyson)
      for (i = 0; i < 4; i++)
      {
        FreeOperator(&(VPerturbation[i]));
        FreeTriDiagonalMatrix(&(gPerturbation[i]));
      }
    FreeOperator(&HAnd);
    FreeOperator(&OPESdn);
    FreeOperator(&OPESup);
    FreeOperator(&OIPESdn);
    FreeOperator(&OIPESup);
  }
  EndTiming("Total", &Time);
  PrintTiming(Time);
  return 0;
}





/* =========================================== old Lanczos functions ============================================= */




/*
*   Eigenstates of a 2 by 2 matrix
*   The lowest eigen vector and energy can be found in Energy[0] and {Functions[0][0],Functions[0][1]}
*   Exept when the matrix on entry is already diagonal or close to diagonal |b1|<(DBL_EPSILON*|a0-a1|)
*     in that case the returned eigen-system are the entered basis functions
*/
void Eigensystem2(const double ain0, const double ain1, const double b1, double Energy[], double Functions[][2])
{
  double shift, a0, a1;
  shift = (ain0 + ain1) / 2.0;
  a0 = ain0 - shift;
  a1 = ain1 - shift;
  if (fabs(b1) < (10 * (DBL_EPSILON * fabs(a0 - a1) + DBL_MIN)))
  {
    if (a0 < a1)
    {
      Energy[0] = a0;
      Energy[1] = a1;
      Functions[0][0] = 1;
      Functions[0][1] = 0;
      Functions[1][0] = 0;
      Functions[1][1] = 1;
    }
    else
    {
      Energy[0] = a1;
      Energy[1] = a0;
      Functions[0][0] = 0;
      Functions[0][1] = 1;
      Functions[1][0] = 1;
      Functions[1][1] = 0;
    }
  }
  else
  {
    Energy[0] = 0.5 * (a0 + a1 - sqrt(Sqr(a0 - a1) + 4 * Sqr(b1)));
    Energy[1] = 0.5 * (a0 + a1 + sqrt(Sqr(a0 - a1) + 4 * Sqr(b1)));
    /* use only the vector for wich |E-a| is largest */
    if (fabs((Energy[0] - a1)) > fabs(Energy[1] - a1))
    {                           /* fun 0 stable */
      Functions[0][0] = (Energy[0] - a1) / sqrt(Sqr(Energy[0] - a1) + Sqr(b1));
      Functions[0][1] = (b1) / sqrt(Sqr(Energy[0] - a1) + Sqr(b1));
      Functions[1][1] = Functions[0][0];
      Functions[1][0] = -Functions[0][1];
    }
    else
    {                           /* fun 1 stable */
      Functions[1][0] = (Energy[1] - a1) / sqrt(Sqr(Energy[1] - a1) + Sqr(b1));
      Functions[1][1] = (b1) / sqrt(Sqr(Energy[1] - a1) + Sqr(b1));
      Functions[0][1] = Functions[1][0];
      Functions[0][0] = -Functions[1][1];
    }
  }
  Energy[0] += shift;
  Energy[1] += shift;
  return;
}

uint32_t Eigensystem3(const double a0, const double a1, const double a2, const double b1, const double b2,
                      double Energy[3], double Functions[3][3])
{
  double e[3] = { b1, b2, 0 };
  double g, r, p, f, b, s, c, t, shift;
  int    nIter;
  int    m, l, i, k;

  shift = (a0 + a1 + a2) / 3.0;
  Energy[0] = a0 - shift;
  Energy[1] = a1 - shift;
  Energy[2] = a2 - shift;
  Functions[0][0] = 1;
  Functions[0][1] = 0;
  Functions[0][2] = 0;
  Functions[1][0] = 0;
  Functions[1][1] = 1;
  Functions[1][2] = 0;
  Functions[2][0] = 0;
  Functions[2][1] = 0;
  Functions[2][2] = 1;

  for (l = 0; l < 2; l++)
  {
    nIter = 0;
    while (1)
    {
      // Check for convergence and exit iteration loop if off-diagonal
      // element e(l) is zero
      for (m = l; m <= 1; m++)
      {
        g = fabs(Energy[m]) + fabs(Energy[m + 1]);
        if (fabs(e[m]) + g == g)
          break;
      }
      if (m == l)
        break;

      if (nIter++ >= 30)
        return 1;

      // Calculate g = d_m - k
      g = (Energy[l + 1] - Energy[l]) / (e[l] + e[l]);
      r = sqrt(Sqr(g) + 1.0);
      if (g > 0)
        g = Energy[m] - Energy[l] + e[l] / (g + r);
      else
        g = Energy[m] - Energy[l] + e[l] / (g - r);

      s = c = 1.0;
      p = 0.0;
      for (i = m - 1; i >= l; i--)
      {
        f = s * e[i];
        b = c * e[i];
        if (fabs(f) > fabs(g))
        {
          c = g / f;
          r = sqrt(Sqr(c) + 1.0);
          e[i + 1] = f * r;
          c *= (s = 1.0 / r);
        }
        else
        {
          s = f / g;
          r = sqrt(Sqr(s) + 1.0);
          e[i + 1] = g * r;
          s *= (c = 1.0 / r);
        }

        g = Energy[i + 1] - p;
        r = (Energy[i] - g) * s + 2.0 * c * b;
        p = s * r;
        Energy[i + 1] = g + p;
        g = c * r - b;

        for (k = 0; k < 3; k++)
        {
          t = Functions[i + 1][k];
          Functions[i + 1][k] = s * Functions[i][k] + c * t;
          Functions[i][k] = c * Functions[i][k] - s * t;
        }

      }
      Energy[l] -= p;
      e[l] = g;
      e[m] = 0.0;
    }
  }
  Energy[0] += shift;
  Energy[1] += shift;
  Energy[2] += shift;
  if ((Energy[0] <= Energy[1]) && (Energy[1] <= Energy[2]))
  {
    return 0;
  }
  if ((Energy[0] <= Energy[2]) && (Energy[2] <= Energy[1]))
  {
    r = Energy[1];
    Energy[1] = Energy[2];
    Energy[2] = r;
    e[0] = Functions[1][0];
    e[1] = Functions[1][1];
    e[2] = Functions[1][2];
    Functions[1][0] = Functions[2][0];
    Functions[1][1] = Functions[2][1];
    Functions[1][2] = Functions[2][2];
    Functions[2][0] = e[0];
    Functions[2][1] = e[1];
    Functions[2][2] = e[2];
    return 0;
  }
  if ((Energy[1] <= Energy[0]) && (Energy[0] <= Energy[2]))
  {
    r = Energy[1];
    Energy[1] = Energy[0];
    Energy[0] = r;
    e[0] = Functions[1][0];
    e[1] = Functions[1][1];
    e[2] = Functions[1][2];
    Functions[1][0] = Functions[0][0];
    Functions[1][1] = Functions[0][1];
    Functions[1][2] = Functions[0][2];
    Functions[0][0] = e[0];
    Functions[0][1] = e[1];
    Functions[0][2] = e[2];
    return 0;
  }
  if ((Energy[2] <= Energy[1]) && (Energy[1] <= Energy[0]))
  {
    r = Energy[2];
    Energy[2] = Energy[0];
    Energy[0] = r;
    e[0] = Functions[2][0];
    e[1] = Functions[2][1];
    e[2] = Functions[2][2];
    Functions[2][0] = Functions[0][0];
    Functions[2][1] = Functions[0][1];
    Functions[2][2] = Functions[0][2];
    Functions[0][0] = e[0];
    Functions[0][1] = e[1];
    Functions[0][2] = e[2];
    return 0;
  }
  if ((Energy[1] <= Energy[2]) && (Energy[2] <= Energy[0]))
  {
    r = Energy[1];
    Energy[1] = Energy[2];
    Energy[2] = Energy[0];
    Energy[0] = r;
    e[0] = Functions[1][0];
    e[1] = Functions[1][1];
    e[2] = Functions[1][2];
    Functions[1][0] = Functions[2][0];
    Functions[1][1] = Functions[2][1];
    Functions[1][2] = Functions[2][2];
    Functions[2][0] = Functions[0][0];
    Functions[2][1] = Functions[0][1];
    Functions[2][2] = Functions[0][2];
    Functions[0][0] = e[0];
    Functions[0][1] = e[1];
    Functions[0][2] = e[2];
    return 0;
  }
  if ((Energy[1] <= Energy[2]) && (Energy[2] <= Energy[0]))
  {
    r = Energy[2];
    Energy[2] = Energy[1];
    Energy[1] = Energy[0];
    Energy[0] = r;
    e[0] = Functions[2][0];
    e[1] = Functions[2][1];
    e[2] = Functions[2][2];
    Functions[2][0] = Functions[1][0];
    Functions[2][1] = Functions[1][1];
    Functions[2][2] = Functions[1][2];
    Functions[1][0] = Functions[0][0];
    Functions[1][1] = Functions[0][1];
    Functions[1][2] = Functions[0][2];
    Functions[0][0] = e[0];
    Functions[0][1] = e[1];
    Functions[0][2] = e[2];
    return 0;
  }
  return 2;
}



/*
*   Returns a single wave function orthogonal to the functions given in OhterPsi with the lowest energy, given
*     H         Hamiltonian
*     psi       Starting wave function
*     epsilon   Value to wich psi is truncated before a single excitation is added
*     zero      Smallest (positive) value that should be distinguised from zero
*     shiftin   Positive value larger then the largest eigenvalue
*     OtherPsi  Wave-functions the solution has to be orthogonal to
*     NOtherPsi Number of wave-functions to orthogonalize to
*
*
*   The values in psi are distroyed and the returned vector has the "same" storage location as the input wave function
*
*   algorithm:
*   while(notconverged)
*     psi = H psi
*     for(i=0;i<NOtherPsi;i++)
*     {
*       psi-=(OtherPsi_i . psi) OtherPsi_i
*     }
*     psi = psi/|psi|
*   end while
*   return psi
*/
void LanczosGroundStateFullRestarted(const OperatorType H, WaveFunctionType * const psi, double *epsilon,
                                     const double zero, double shiftin, WaveFunctionType * OtherPsi, uint32_t NOtherPsi)
{
  printf("Start of LanczosGroundState\n\n");
  fflush(stdout);
  WaveFunctionType Opsi;
  double energy = 4 * zero, energy1 = 2 * zero, energy2;
  double shift;
  uint32_t cycles;
  while ((energy - energy1) > (zero * (1.0 + fabs(energy))))
  {
    cycles = 0;
    energy = energy1;
    shift = -shiftin - energy * 0.49;
    while (OperatorPsi(H, *psi, &Opsi)) /* Newly allocate memory for Opsi both in Value, Determinant and Next */
    {                           /* Allocation Failed */
      printf("Allocation of new memory failed in LanczosGroundState for the function OperatorPsi\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", *epsilon,
             *epsilon * 2.0);
      fflush(stdout);
      *epsilon *= 2.0;
      if (psi->NBlok == 1)
        AllocFailed();
      RealWaveFunctionTruncateBasis(psi, *epsilon);
    }
    printf("BasisSetSizeAfter Expansion/Truncation : %12i ", Opsi.N);
    RealWaveFunctionMultiply(*psi, shift);
    while (RealWaveFunctionAddUnionBasis(&Opsi, *psi))  /* Possibly expand allocated memory of Opsi */
    {                           /* Allocation Failed */
      printf("Allocation of new memory failed in LanczosGroundState for the function RealWaveFunctionAddUnionBasis\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", *epsilon,
             *epsilon * 2.0);
      fflush(stdout);
      *epsilon *= 2.0;
      RealWaveFunctionFreeValues(&Opsi);
      RealWaveFunctionMultiply(*psi, 1.0 / shift);
      if (psi->NBlok == 1)
        AllocFailed();
      RealWaveFunctionTruncateBasis(psi, *epsilon);
      OperatorPsi(H, *psi, &Opsi);
      RealWaveFunctionMultiply(*psi, shift);
    }
    RealWaveFunctionFree(psi);  /* Remove allocated memory for psi from memory */
    while (RealWaveFunctionOrthogonalize(&Opsi, OtherPsi, NOtherPsi))
    {                           /* Allocation Failed */
      printf("Allocation of new memory failed in LanczosGroundState for the function RealWaveFunctionOrthogonalize\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", *epsilon,
             *epsilon * 2.0);
      fflush(stdout);
      *epsilon *= 2.0;
      if (Opsi.NBlok == 1)
        AllocFailed();
      RealWaveFunctionTruncateBasis(&Opsi, *epsilon);
    }
    RealWaveFunctionNormalize(&Opsi);
    RealWaveFunctionTruncateBasis(&Opsi, *epsilon); /* Possibly remove allocated memory for Opsi */
    printf("%12i\n", Opsi.N);
    while (OperatorPsiConserveBasis(H, Opsi, psi))  /* Newly allocate memory for psi in Value. Determinant and Next point to same positions as Opsi */
    {                           /* Allocation Failed */
      printf("Allocation of new memory failed in LanczosGroundState for the function OperatorPsiConserveBasis\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", *epsilon,
             *epsilon * 2.0);
      fflush(stdout);
      *epsilon *= 2.0;
      if (Opsi.NBlok == 1)
        AllocFailed();
      RealWaveFunctionTruncateBasis(&Opsi, *epsilon);
    }
    RealWaveFunctionMultiply(Opsi, shift);
    RealWaveFunctionAddSameBasis(*psi, Opsi);
    RealWaveFunctionNormalize(psi);
    energy1 = 2 * zero;
    energy2 = 4 * zero;
    while ((energy2 - energy1) > (zero * (1.0 + fabs(energy1))))
    {
      cycles++;
      if ((cycles & 255) == 255)
      {
        printf("Reached %i-th cycle with energy %22.15E\n", cycles, energy1);
        fflush(stdout);
      }
      OperatorPsiConserveBasisPreAllocated(H, *psi, &Opsi);
      RealWaveFunctionOrthogonalizeConserveBasis(&Opsi, OtherPsi, NOtherPsi);
      energy2 = RealWaveFunctionDotProductSameBasis(*psi, Opsi);
      RealWaveFunctionMultiply(*psi, shift);
      RealWaveFunctionAddSameBasis(Opsi, *psi);
      RealWaveFunctionNormalize(&Opsi);
      OperatorPsiConserveBasisPreAllocated(H, Opsi, psi);
      RealWaveFunctionOrthogonalizeConserveBasis(psi, OtherPsi, NOtherPsi);
      energy1 = RealWaveFunctionDotProductSameBasis(*psi, Opsi);
      RealWaveFunctionMultiply(Opsi, shift);
      RealWaveFunctionAddSameBasis(*psi, Opsi);
      RealWaveFunctionNormalize(psi);
    }
    RealWaveFunctionFreeValues(&Opsi);  /* Free Opsi.ValueRe, but don't free Opsi.Determinant as this is the same space as where psi.Determinant points to */
    printf("Energy: %22.15E  ", energy1);
    printf("NumberOfLaczosCycles: %i\n", cycles);
  }
  printf("\n");
  fflush(stdout);
}

void LanczosGroundStateDouble(const OperatorType H, WaveFunctionType * const psi, double *epsilon, const double zero,
                              double shiftin, WaveFunctionType * OtherPsi, uint32_t NOtherPsi)
{
  printf("Start of LanczosGroundState\n\n");
  fflush(stdout);
  WaveFunctionType Opsi;
  double energy = 4 * zero, energy1 = 2 * zero, energy2;
  double shift, a0, a1, b1;
  uint32_t cycles;
  while ((energy - energy1) > (0.5 * zero * (1.0 + fabs(energy))))
  {
    cycles = 0;
    energy = energy1;
    shift = -shiftin - energy * 0.49;
    while (OperatorPsi(H, *psi, &Opsi)) /* Newly allocate memory for Opsi both in Value, Determinant and Next */
    {                           /* Allocation Failed */
      printf("Allocation of new memory failed in LanczosGroundState for the function OperatorPsi\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", *epsilon,
             *epsilon * 2.0);
      fflush(stdout);
      *epsilon *= 2.0;
      if (psi->NBlok == 1)
        AllocFailed();
      RealWaveFunctionTruncateBasis(psi, *epsilon);
    }
    printf("BasisSetSizeAfter Expansion/Truncation : %12i ", Opsi.N);
    RealWaveFunctionMultiply(*psi, shift);
    while (RealWaveFunctionAddUnionBasis(&Opsi, *psi))  /* Possibly expand allocated memory of Opsi */
    {                           /* Allocation Failed */
      printf("Allocation of new memory failed in LanczosGroundState for the function RealWaveFunctionAddUnionBasis\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", *epsilon,
             *epsilon * 2.0);
      fflush(stdout);
      *epsilon *= 2.0;
      RealWaveFunctionFreeValues(&Opsi);
      RealWaveFunctionMultiply(*psi, 1.0 / shift);
      if (psi->NBlok == 1)
        AllocFailed();
      RealWaveFunctionTruncateBasis(psi, *epsilon);
      OperatorPsi(H, *psi, &Opsi);
      RealWaveFunctionMultiply(*psi, shift);
    }
    RealWaveFunctionFree(psi);  /* Remove allocated memory for psi from memory */
    while (RealWaveFunctionOrthogonalize(&Opsi, OtherPsi, NOtherPsi))
    {                           /* Allocation Failed */
      printf("Allocation of new memory failed in LanczosGroundState for the function RealWaveFunctionOrthogonalize\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", *epsilon,
             *epsilon * 2.0);
      fflush(stdout);
      *epsilon *= 2.0;
      if (Opsi.NBlok == 1)
        AllocFailed();
      RealWaveFunctionTruncateBasis(&Opsi, *epsilon);
    }
    RealWaveFunctionNormalize(&Opsi);
    RealWaveFunctionTruncateBasis(&Opsi, *epsilon); /* Possibly remove allocated memory for Opsi */
    printf("%12i\n", Opsi.N);
    while (OperatorPsiConserveBasis(H, Opsi, psi))  /* Newly allocate memory for psi in Value. Determinant and Next point to same positions as Opsi */
    {                           /* Allocation Failed */
      printf("Allocation of new memory failed in LanczosGroundState for the function OperatorPsiConserveBasis\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", *epsilon,
             *epsilon * 2.0);
      fflush(stdout);
      *epsilon *= 2.0;
      if (Opsi.NBlok == 1)
        AllocFailed();
      RealWaveFunctionTruncateBasis(&Opsi, *epsilon);
    }
    RealWaveFunctionMultiply(Opsi, shift);
    RealWaveFunctionAddSameBasis(*psi, Opsi);
    RealWaveFunctionNormalize(psi);
    energy1 = 2 * zero;
    energy2 = 4 * zero;
    while ((energy2 - energy1) > (zero * (1.0 + fabs(energy1))))
    {
      cycles++;
      if ((cycles & 255) == 255)
      {
        printf("Reached %i-th cycle with energy %22.15E\n", cycles, energy1);
        fflush(stdout);
      }
      energy2 = energy1;
      OperatorPsiConserveBasisPreAllocated(H, *psi, &Opsi);
      a0 = RealWaveFunctionDotProductSameBasis(*psi, Opsi);
      RealWaveFunctionAddMultipliedSameBasis(Opsi, *psi, -a0);
      b1 = RealWaveFunctionCalculateNorm(Opsi);
      if (b1 == 0)
      {
        energy1 = a0;
      }
      else
      {
        RealWaveFunctionMultiply(Opsi, (double)1.0 / b1);
        PsiOperatorPsiRR(H, Opsi, &a1);
        energy1 = 0.5 * (a0 + a1 - sqrt(a0 * a0 - 2 * a0 * a1 + a1 * a1 + 4 * b1 * b1));
        RealWaveFunctionMultiply(*psi, (energy1 - a1) * 2);
        RealWaveFunctionAddMultipliedSameBasis(*psi, Opsi, 2 * b1);
        RealWaveFunctionOrthogonalizeConserveBasis(psi, OtherPsi, NOtherPsi);
        RealWaveFunctionNormalize(psi);
      }
    }
    RealWaveFunctionFreeValues(&Opsi);  /* Free  Opsi.ValueRe, but don't free  Opsi.Determinant as this is the same space as where psi.Determinant points to */
    printf("Energy: %22.15E  ", energy1);
    printf("NumberOfLaczosCycles: %i\n", cycles);
  }
  printf("\n");
  fflush(stdout);
}

uint32_t LanczosGroundStateRRTripple(const OperatorType H, WaveFunctionType * const psi, double *epsilon,
                                     const double zero, double shiftin, WaveFunctionType * OtherPsi, uint32_t NOtherPsi)
{
  printf("Start of LanczosGroundStateRR\n\n");
  fflush(stdout);
  WaveFunctionType Opsi, OOpsi;
  double energy = 4 * zero, energy1 = 2 * zero, energy2;
  double shift, a0, a1, a2, b1, b2;
  double EnT[3], FunT[3][3];
  uint32_t cycles;

  while (RealWaveFunctionOrthogonalize(psi, OtherPsi, NOtherPsi))
  {                             /* Allocation Failed */
    printf("Allocation of new memory failed in LanczosGroundState for the function RealWaveFunctionOrthogonalize\n");
    printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", *epsilon,
           *epsilon * 2.0);
    fflush(stdout);
    *epsilon *= 2.0;
    if (Opsi.NBlok == 1)
      AllocFailed();
    RealWaveFunctionTruncateBasis(psi, *epsilon);
  }
  RealWaveFunctionNormalize(psi);
  RealWaveFunctionTruncateBasis(psi, *epsilon);

  while ((energy - energy1) > (zero * (1.0 + fabs(energy))))
  {
    cycles = 0;
    energy = energy1;
    shift = -shiftin - energy * 0.49;
    while (OperatorPsiRR(H, *psi, &Opsi)) /* Newly allocate memory for Opsi both in Value, Determinant and Next */
    {                           /* Allocation Failed */
      printf("Allocation of new memory failed in LanczosGroundState for the function OperatorPsi\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", *epsilon,
             *epsilon * 2.0);
      fflush(stdout);
      *epsilon *= 2.0;
      if (psi->NBlok == 1)
        AllocFailed();
      RealWaveFunctionTruncateBasis(psi, *epsilon);
    }
    printf("BasisSetSizeAfter Expansion/Truncation : %12i ", Opsi.N);
    RealWaveFunctionMultiply(*psi, shift);
    while (RealWaveFunctionAddUnionBasis(&Opsi, *psi))  /* Possibly expand allocated memory of Opsi */
    {                           /* Allocation Failed */
      printf("Allocation of new memory failed in LanczosGroundState for the function RealWaveFunctionAddUnionBasis\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", *epsilon,
             *epsilon * 2.0);
      fflush(stdout);
      *epsilon *= 2.0;
      RealWaveFunctionFreeValues(&Opsi);
      RealWaveFunctionMultiply(*psi, 1.0 / shift);
      if (psi->NBlok == 1)
        AllocFailed();
      RealWaveFunctionTruncateBasis(psi, *epsilon);
      OperatorPsiRR(H, *psi, &Opsi);
      RealWaveFunctionMultiply(*psi, shift);
    }
    RealWaveFunctionFree(psi);  /* Remove allocated memory for psi from memory */
    while (RealWaveFunctionOrthogonalize(&Opsi, OtherPsi, NOtherPsi))
    {                           /* Allocation Failed */
      printf("Allocation of new memory failed in LanczosGroundState for the function RealWaveFunctionOrthogonalize\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", *epsilon,
             *epsilon * 2.0);
      fflush(stdout);
      *epsilon *= 2.0;
      if (Opsi.NBlok == 1)
        AllocFailed();
      RealWaveFunctionTruncateBasis(&Opsi, *epsilon);
    }
    RealWaveFunctionNormalize(&Opsi);
    RealWaveFunctionTruncateBasis(&Opsi, *epsilon); /* Possibly remove allocated memory for Opsi */
    printf("%12i\n", Opsi.N);
    while (OperatorPsiConserveBasisRR(H, Opsi, &OOpsi)) /* Newly allocate memory for OOpsi in Value. Determinant and Next point to same positions as Opsi */
    {                           /* Allocation Failed */
      printf("Allocation of new memory failed in LanczosGroundState for the function OperatorPsiConserveBasis\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", *epsilon,
             *epsilon * 2.0);
      fflush(stdout);
      *epsilon *= 2.0;
      if (Opsi.NBlok == 1)
        AllocFailed();
      RealWaveFunctionTruncateBasis(&Opsi, *epsilon);
    }
    RealWaveFunctionMultiply(Opsi, shift);
    RealWaveFunctionAddSameBasis(OOpsi, Opsi);
    RealWaveFunctionNormalize(&OOpsi);
    while (OperatorPsiConserveBasisRR(H, OOpsi, psi)) /* Newly allocate memory for psi in Value. Determinant and Next point to same positions as Opsi */
    {                           /* Allocation Failed */
      printf("Allocation of new memory failed in LanczosGroundState for the function OperatorPsiConserveBasis\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", *epsilon,
             *epsilon * 2.0);
      fflush(stdout);
      *epsilon *= 2.0;
      if (Opsi.NBlok == 1)
        AllocFailed();
      RealWaveFunctionTruncateBasis(&Opsi, *epsilon);
      RealWaveFunctionFreeValues(&OOpsi);
      OperatorPsiConserveBasisRR(H, Opsi, &OOpsi);
    }
    RealWaveFunctionMultiply(OOpsi, shift);
    RealWaveFunctionAddSameBasis(*psi, OOpsi);
    RealWaveFunctionNormalize(psi);
    energy1 = 2 * zero;
    energy2 = 4 * zero;
    while ((cycles < 4) || ((energy2 - energy1) > (zero * (1.0 + fabs(energy1)))))
    {
      cycles++;
      if ((cycles & 255) == 255)
      {
        printf("Reached %i-th cycle with energy %22.15E\n", cycles, energy1);
        fflush(stdout);
      }
      energy2 = energy1;
      OperatorPsiConserveBasisPreAllocatedRR(H, *psi, &Opsi);
      a0 = RealWaveFunctionDotProductSameBasis(*psi, Opsi);
      RealWaveFunctionAddMultipliedSameBasis(Opsi, *psi, -a0);
      b1 = RealWaveFunctionCalculateNorm(Opsi);
      if (b1 == 0)
      {
        a1 = fabs(a0) + 2.0;
        a2 = fabs(a0) + fabs(a1) + 2.0;
        b2 = 0;
      }
      else
      {
        RealWaveFunctionMultiply(Opsi, (double)1.0 / b1);
        OperatorPsiConserveBasisPreAllocatedRR(H, Opsi, &OOpsi);
        a1 = RealWaveFunctionDotProductSameBasis(Opsi, OOpsi);
        RealWaveFunctionAddMultipliedSameBasis(OOpsi, Opsi, -a1);
        b2 = RealWaveFunctionCalculateNorm(OOpsi);
        if (b2 == 0)
        {
          a2 = fabs(a0) + fabs(a1) + 2.0;
        }
        else
        {
          RealWaveFunctionMultiply(OOpsi, (double)1.0 / b2);
          PsiOperatorPsiRR(H, OOpsi, &a2);
        }
      }
      Eigensystem3(a0, a1, a2, b1, b2, EnT, FunT);
      energy1 = EnT[0];
      RealWaveFunctionMultiply(*psi, FunT[0][0]);
      RealWaveFunctionAddMultipliedSameBasis(*psi, Opsi, FunT[0][1]);
      RealWaveFunctionAddMultipliedSameBasis(*psi, OOpsi, FunT[0][2]);
      RealWaveFunctionOrthogonalizeConserveBasis(psi, OtherPsi, NOtherPsi);
      RealWaveFunctionNormalize(psi);
    }
    RealWaveFunctionFreeValues(&Opsi);  /* Free  Opsi.ValueRe, but don't free  Opsi.Determinant as this is the same space as where psi.Determinant points to */
    RealWaveFunctionFreeValues(&OOpsi); /* Free OOpsi.ValueRe, but don't free OOpsi.Determinant as this is the same space as where psi.Determinant points to */
    printf("Energy: %22.15E  ", energy1);
    printf("NumberOfLaczosCycles: %i\n", cycles);
  }
  printf("\n");
  fflush(stdout);
  return 0;
}


uint32_t LanczosGroundStateRCTripple(const OperatorType H, WaveFunctionType * const psi, double *epsilon,
                                     const double zero, double shiftin, WaveFunctionType * OtherPsi, uint32_t NOtherPsi)
{
  printf("Start of LanczosGroundStateRC\n\n");
  fflush(stdout);
  WaveFunctionType Opsi, OOpsi;
  double energy = 4 * zero, energy1 = 2 * zero, energy2;
  double shift, a0, a1, a2, b1, b2, a0Im, a1Im;
  double EnT[3], FunT[3][3];
  uint32_t cycles;
  while ((energy - energy1) > (zero * (1.0 + fabs(energy))))
  {
    cycles = 0;
    energy = energy1;
    shift = -shiftin - energy * 0.49;
    while (OperatorPsiRC(H, *psi, &Opsi)) /* Newly allocate memory for Opsi both in Value, Determinant and Next */
    {                           /* Allocation Failed */
      printf("Allocation of new memory failed in LanczosGroundState for the function OperatorPsi\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", *epsilon,
             *epsilon * 2.0);
      fflush(stdout);
      *epsilon *= 2.0;
      if (psi->NBlok == 1)
        AllocFailed();
      ComplexWaveFunctionTruncateBasis(psi, *epsilon);
    }
    printf("BasisSetSizeAfter Expansion/Truncation : %12i ", Opsi.N);
    ComplexWaveFunctionMultiply(*psi, shift);
    while (ComplexWaveFunctionAddUnionBasis(&Opsi, *psi)) /* Possibly expand allocated memory of Opsi */
    {                           /* Allocation Failed */
      printf("Allocation of new memory failed in LanczosGroundState for the function RealWaveFunctionAddUnionBasis\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", *epsilon,
             *epsilon * 2.0);
      fflush(stdout);
      *epsilon *= 2.0;
      ComplexWaveFunctionFreeValues(&Opsi);
      ComplexWaveFunctionMultiply(*psi, 1.0 / shift);
      if (psi->NBlok == 1)
        AllocFailed();
      ComplexWaveFunctionTruncateBasis(psi, *epsilon);
      OperatorPsiRC(H, *psi, &Opsi);
      ComplexWaveFunctionMultiply(*psi, shift);
    }
    ComplexWaveFunctionFree(psi); /* Remove allocated memory for psi from memory */
    while (ComplexWaveFunctionOrthogonalize(&Opsi, OtherPsi, NOtherPsi))
    {                           /* Allocation Failed */
      printf("Allocation of new memory failed in LanczosGroundState for the function RealWaveFunctionOrthogonalize\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", *epsilon,
             *epsilon * 2.0);
      fflush(stdout);
      *epsilon *= 2.0;
      if (Opsi.NBlok == 1)
        AllocFailed();
      ComplexWaveFunctionTruncateBasis(&Opsi, *epsilon);
    }
    ComplexWaveFunctionNormalize(&Opsi);
    ComplexWaveFunctionTruncateBasis(&Opsi, *epsilon);  /* Possibly remove allocated memory for Opsi */
    printf("%12i\n", Opsi.N);
    while (OperatorPsiConserveBasisRC(H, Opsi, &OOpsi)) /* Newly allocate memory for OOpsi in Value. Determinant and Next point to same positions as Opsi */
    {                           /* Allocation Failed */
      printf("Allocation of new memory failed in LanczosGroundState for the function OperatorPsiConserveBasis\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", *epsilon,
             *epsilon * 2.0);
      fflush(stdout);
      *epsilon *= 2.0;
      if (Opsi.NBlok == 1)
        AllocFailed();
      ComplexWaveFunctionTruncateBasis(&Opsi, *epsilon);
    }
    ComplexWaveFunctionMultiply(Opsi, shift);
    ComplexWaveFunctionAddSameBasis(OOpsi, Opsi);
    ComplexWaveFunctionNormalize(&OOpsi);
    while (OperatorPsiConserveBasisRC(H, OOpsi, psi)) /* Newly allocate memory for psi in Value. Determinant and Next point to same positions as Opsi */
    {                           /* Allocation Failed */
      printf("Allocation of new memory failed in LanczosGroundState for the function OperatorPsiConserveBasis\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", *epsilon,
             *epsilon * 2.0);
      fflush(stdout);
      *epsilon *= 2.0;
      if (Opsi.NBlok == 1)
        AllocFailed();
      ComplexWaveFunctionTruncateBasis(&Opsi, *epsilon);
      ComplexWaveFunctionFreeValues(&OOpsi);
      OperatorPsiConserveBasisRC(H, Opsi, &OOpsi);
    }
    ComplexWaveFunctionMultiply(OOpsi, shift);
    ComplexWaveFunctionAddSameBasis(*psi, OOpsi);
    ComplexWaveFunctionNormalize(psi);
    energy1 = 2 * zero;
    energy2 = 4 * zero;
    while ((cycles < 4) || ((energy2 - energy1) > (zero * (1.0 + fabs(energy1)))))
    {
      cycles++;
      if ((cycles & 255) == 255)
      {
        printf("Reached %i-th cycle with energy %22.15E\n", cycles, energy1);
        fflush(stdout);
      }
      energy2 = energy1;
      OperatorPsiConserveBasisPreAllocatedRC(H, *psi, &Opsi);
      ComplexWaveFunctionDotProductSameBasis(*psi, Opsi, &a0, &a0Im);
      ComplexWaveFunctionAddMultipliedSameBasis(Opsi, *psi, -a0);
      b1 = ComplexWaveFunctionCalculateNorm(Opsi);
      if (b1 == 0)
      {
        a1 = fabs(a0) + 2.0;
        a2 = fabs(a0) + fabs(a1) + 2.0;
        b2 = 0;
      }
      else
      {
        ComplexWaveFunctionMultiply(Opsi, (double)1.0 / b1);
        OperatorPsiConserveBasisPreAllocatedRC(H, Opsi, &OOpsi);
        ComplexWaveFunctionDotProductSameBasis(Opsi, OOpsi, &a1, &a1Im);
        ComplexWaveFunctionAddMultipliedSameBasis(OOpsi, Opsi, -a1);
        b2 = ComplexWaveFunctionCalculateNorm(OOpsi);
        if (b2 == 0)
        {
          a2 = fabs(a0) + fabs(a1) + 2.0;
        }
        else
        {
          ComplexWaveFunctionMultiply(OOpsi, (double)1.0 / b2);
          PsiOperatorPsiRHC(H, OOpsi, &a2);
        }
      }
      Eigensystem3(a0, a1, a2, b1, b2, EnT, FunT);
      energy1 = EnT[0];
      ComplexWaveFunctionMultiply(*psi, FunT[0][0]);
      ComplexWaveFunctionAddMultipliedSameBasis(*psi, Opsi, FunT[0][1]);
      ComplexWaveFunctionAddMultipliedSameBasis(*psi, OOpsi, FunT[0][2]);
      ComplexWaveFunctionOrthogonalizeConserveBasis(psi, OtherPsi, NOtherPsi);
      ComplexWaveFunctionNormalize(psi);
    }
    ComplexWaveFunctionFreeValues(&Opsi); /* Free  Opsi.ValueRe, but don't free  Opsi.Determinant as this is the same space as where psi.Determinant points to */
    ComplexWaveFunctionFreeValues(&OOpsi);  /* Free OOpsi.ValueRe, but don't free OOpsi.Determinant as this is the same space as where psi.Determinant points to */
    printf("Energy: %22.15E  ", energy1);
    printf("NumberOfLaczosCycles: %i\n", cycles);
  }
  printf("\n");
  fflush(stdout);
  return 0;
}

uint32_t LanczosGroundStateCCTripple(const OperatorType H, WaveFunctionType * const psi, double *epsilon,
                                     const double zero, double shiftin, WaveFunctionType * OtherPsi, uint32_t NOtherPsi)
{
  printf("Start of LanczosGroundStateCC\n\n");
  fflush(stdout);
  WaveFunctionType Opsi, OOpsi;
  double energy = 4 * zero, energy1 = 2 * zero, energy2;
  double shift, a0, a1, a2, b1, b2, a0Im, a1Im;
  double EnT[3], FunT[3][3];
  uint32_t cycles;
  while ((energy - energy1) > (zero * (1.0 + fabs(energy))))
  {
    cycles = 0;
    energy = energy1;
    shift = -shiftin - energy * 0.49;
    while (OperatorPsiCC(H, *psi, &Opsi)) /* Newly allocate memory for Opsi both in Value, Determinant and Next */
    {                           /* Allocation Failed */
      printf("Allocation of new memory failed in LanczosGroundState for the function OperatorPsi\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", *epsilon,
             *epsilon * 2.0);
      fflush(stdout);
      *epsilon *= 2.0;
      if (psi->NBlok == 1)
        AllocFailed();
      ComplexWaveFunctionTruncateBasis(psi, *epsilon);
    }
    printf("BasisSetSizeAfter Expansion/Truncation : %12i ", Opsi.N);
    ComplexWaveFunctionMultiply(*psi, shift);
    while (ComplexWaveFunctionAddUnionBasis(&Opsi, *psi)) /* Possibly expand allocated memory of Opsi */
    {                           /* Allocation Failed */
      printf("Allocation of new memory failed in LanczosGroundState for the function RealWaveFunctionAddUnionBasis\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", *epsilon,
             *epsilon * 2.0);
      fflush(stdout);
      *epsilon *= 2.0;
      ComplexWaveFunctionFreeValues(&Opsi);
      ComplexWaveFunctionMultiply(*psi, 1.0 / shift);
      if (psi->NBlok == 1)
        AllocFailed();
      ComplexWaveFunctionTruncateBasis(psi, *epsilon);
      OperatorPsiCC(H, *psi, &Opsi);
      ComplexWaveFunctionMultiply(*psi, shift);
    }
    ComplexWaveFunctionFree(psi); /* Remove allocated memory for psi from memory */
    while (ComplexWaveFunctionOrthogonalize(&Opsi, OtherPsi, NOtherPsi))
    {                           /* Allocation Failed */
      printf("Allocation of new memory failed in LanczosGroundState for the function RealWaveFunctionOrthogonalize\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", *epsilon,
             *epsilon * 2.0);
      fflush(stdout);
      *epsilon *= 2.0;
      if (Opsi.NBlok == 1)
        AllocFailed();
      ComplexWaveFunctionTruncateBasis(&Opsi, *epsilon);
    }
    ComplexWaveFunctionNormalize(&Opsi);
    ComplexWaveFunctionTruncateBasis(&Opsi, *epsilon);  /* Possibly remove allocated memory for Opsi */
    printf("%12i\n", Opsi.N);
    while (OperatorPsiConserveBasisCC(H, Opsi, &OOpsi)) /* Newly allocate memory for OOpsi in Value. Determinant and Next point to same positions as Opsi */
    {                           /* Allocation Failed */
      printf("Allocation of new memory failed in LanczosGroundState for the function OperatorPsiConserveBasis\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", *epsilon,
             *epsilon * 2.0);
      fflush(stdout);
      *epsilon *= 2.0;
      if (Opsi.NBlok == 1)
        AllocFailed();
      ComplexWaveFunctionTruncateBasis(&Opsi, *epsilon);
    }
    ComplexWaveFunctionMultiply(Opsi, shift);
    ComplexWaveFunctionAddSameBasis(OOpsi, Opsi);
    ComplexWaveFunctionNormalize(&OOpsi);
    while (OperatorPsiConserveBasisCC(H, OOpsi, psi)) /* Newly allocate memory for psi in Value. Determinant and Next point to same positions as Opsi */
    {                           /* Allocation Failed */
      printf("Allocation of new memory failed in LanczosGroundState for the function OperatorPsiConserveBasis\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", *epsilon,
             *epsilon * 2.0);
      fflush(stdout);
      *epsilon *= 2.0;
      if (Opsi.NBlok == 1)
        AllocFailed();
      ComplexWaveFunctionTruncateBasis(&Opsi, *epsilon);
      ComplexWaveFunctionFreeValues(&OOpsi);
      OperatorPsiConserveBasisCC(H, Opsi, &OOpsi);
    }
    ComplexWaveFunctionMultiply(OOpsi, shift);
    ComplexWaveFunctionAddSameBasis(*psi, OOpsi);
    ComplexWaveFunctionNormalize(psi);
    energy1 = 2 * zero;
    energy2 = 4 * zero;
    while ((cycles < 4) || ((energy2 - energy1) > (zero * (1.0 + fabs(energy1)))))
    {
      cycles++;
      if ((cycles & 255) == 255)
      {
        printf("Reached %i-th cycle with energy %22.15E\n", cycles, energy1);
        fflush(stdout);
      }
      energy2 = energy1;
      OperatorPsiConserveBasisPreAllocatedCC(H, *psi, &Opsi);
      ComplexWaveFunctionDotProductSameBasis(*psi, Opsi, &a0, &a0Im);
      ComplexWaveFunctionAddMultipliedSameBasis(Opsi, *psi, -a0);
      b1 = ComplexWaveFunctionCalculateNorm(Opsi);
      if (b1 == 0)
      {
        a1 = fabs(a0) + 2.0;
        a2 = fabs(a0) + fabs(a1) + 2.0;
        b2 = 0;
      }
      else
      {
        ComplexWaveFunctionMultiply(Opsi, (double)1.0 / b1);
        OperatorPsiConserveBasisPreAllocatedCC(H, Opsi, &OOpsi);
        ComplexWaveFunctionDotProductSameBasis(Opsi, OOpsi, &a1, &a1Im);
        ComplexWaveFunctionAddMultipliedSameBasis(OOpsi, Opsi, -a1);
        b2 = ComplexWaveFunctionCalculateNorm(OOpsi);
        if (b2 == 0)
        {
          a2 = fabs(a0) + fabs(a1) + 2.0;
        }
        else
        {
          ComplexWaveFunctionMultiply(OOpsi, (double)1.0 / b2);
          PsiOperatorPsiCHC(H, OOpsi, &a2);
        }
      }
      Eigensystem3(a0, a1, a2, b1, b2, EnT, FunT);
      energy1 = EnT[0];
      ComplexWaveFunctionMultiply(*psi, FunT[0][0]);
      ComplexWaveFunctionAddMultipliedSameBasis(*psi, Opsi, FunT[0][1]);
      ComplexWaveFunctionAddMultipliedSameBasis(*psi, OOpsi, FunT[0][2]);
      ComplexWaveFunctionOrthogonalizeConserveBasis(psi, OtherPsi, NOtherPsi);
      ComplexWaveFunctionNormalize(psi);
    }
    ComplexWaveFunctionFreeValues(&Opsi); /* Free  Opsi.ValueRe, but don't free  Opsi.Determinant as this is the same space as where psi.Determinant points to */
    ComplexWaveFunctionFreeValues(&OOpsi);  /* Free OOpsi.ValueRe, but don't free OOpsi.Determinant as this is the same space as where psi.Determinant points to */
    printf("Energy: %22.15E  ", energy1);
    printf("NumberOfLaczosCycles: %i\n", cycles);
  }
  printf("\n");
  fflush(stdout);
  return 0;
}


void LanczosTriDiagonalizeConserveBasisPreAllocated(const OperatorType H, double shift, WaveFunctionType v[],
                                                    TriDiagonalMatrixType * const M)
{
  uint32_t i, vi, vip, vim;
  M->b[0] = 1.0;
  OperatorPsiConserveBasisPreAllocated(H, v[0], &(v[1]));
  RealWaveFunctionAddMultipliedConserveBasis(&(v[1]), v[0], shift);
  M->a[0] = RealWaveFunctionDotProduct(v[1], v[0]);
  for (i = 1; i < M->N; i++)
  {
    vip = (i + 1) % 3;
    vi = (i) % 3;
    vim = (i - 1) % 3;
    RealWaveFunctionAddMultipliedConserveBasis(&(v[vi]), v[vim], -M->a[i - 1]);
    M->b[i] = RealWaveFunctionCalculateNorm(v[vi]);
    RealWaveFunctionMultiply(v[vi], (double)1.0 / M->b[i]);
    OperatorPsiConserveBasisPreAllocated(H, v[vi], &(v[vip]));
    RealWaveFunctionAddMultipliedConserveBasis(&(v[vip]), v[vim], shift);
    M->a[i] = RealWaveFunctionDotProduct(v[vip], v[vi]);
    RealWaveFunctionAddMultipliedConserveBasis(&(v[vip]), v[vim], -M->b[i]);
  }
}


/* =========================================== old DMFT functions ============================================= */



int OldmainWorkingDMFT(int argc, const char *argv[])
{
  TimeType Time;
  DMFTSettings HSet;
  InitializeTiming(&Time);
  BeginTiming("Total", &Time);
  TimeStampAll();
  argc -= 1;
  /* read the DMFT settings */
  if (FileReadDMFTSettings(&HSet, argv[argc]))
  {
    printf("Fatal Error in main\n Failed to read DMFT Settings from file %s\n", argv[argc]);
    return 1;
  };
  PrintDMFTSettings(HSet);

  WaveFunctionType Psi, OPsi, PsiTemp;
  TriDiagonalMatrixType G0, GBath, GclPup, GclPdn, GclMup, GclMdn, Gclup, Gcldn, Gcl, Sigma, GBathNew, GBathOld,
    gPerturbation[4];
  SecondDoubleTriDiagonalMatrixType GAnd, GAndTrunk;
  OperatorType HAnd, OPESdn, OPESup, OIPESdn, OIPESup, VPerturbation[4];
  MatrixType T;
  SpectraType Spectra[10];
  char   GBathFileName[255], PsiFileName[255], SpecFileName[25], Buffer[255];
  uint32_t iOuter, i, PsiNBathV, PsiNBathC;
  uint16_t oadd[1];
  double Ef = 0.0, epsilonAnd = 0.00001, mu;
  double PsiEnergy, PsiQuality, PsiHist[20], *PsiRho;
  /* read G0 */
  if (FileReadPlainTriDiagonalMatrixType(&G0, HSet.G0FileName))
  {
    printf("Fatal Error in main\n Failed to read G0 from file %s\n", HSet.G0FileName);
    return 1;
  };
  strcpy(G0.Name, "G0");
  /* read G_Bath_iStart and psi_iStart */
  if (HSet.iStart == 0)
  {
    sprintf(GBathFileName, "%s000", HSet.GRootFileName);
    sprintf(PsiFileName, "%s000", HSet.WaveFunctionRootFileName);
  }
  else if (HSet.iStart < 10)
  {
    sprintf(GBathFileName, "%s00%1i", HSet.GRootFileName, HSet.iStart);
    sprintf(PsiFileName, "%s00%1i", HSet.WaveFunctionRootFileName, HSet.iStart);
  }
  else if (HSet.iStart < 100)
  {
    sprintf(GBathFileName, "%s0%2i", HSet.GRootFileName, HSet.iStart);
    sprintf(PsiFileName, "%s0%2i", HSet.WaveFunctionRootFileName, HSet.iStart);
  }
  else
  {
    sprintf(GBathFileName, "%s%i", HSet.GRootFileName, HSet.iStart);
    sprintf(PsiFileName, "%s%i", HSet.WaveFunctionRootFileName, HSet.iStart);
  }
  if (FileReadPlainTriDiagonalMatrixType(&GBath, GBathFileName))
  {
    printf("Fatal Error in main\n Failed to read GBath from file %s !\n", GBathFileName);
    return 1;
  };
  if (file_exists(PsiFileName))
  {
    Psi.NBitsKey = HSet.NBitsKey;
    if (FileReadWaveFunctionPlainBinary(&Psi, PsiFileName))
    {
      printf("Fatal Error in main\n Failed to read Psi from file %s !\n", PsiFileName);
      return 1;
    };
    printf("\nRead the starting wave function from the file %s\n", PsiFileName);
    sscanf(Psi.Name, "%s %i %i", Buffer, &PsiNBathC, &PsiNBathV);
    printf("Wave-function \"%s\" found has %i conduction bath sites, %i valence bath sites and %i Fermions in basis\n",
           Psi.Name, PsiNBathC, PsiNBathV, Psi.NFermion);
    if ((2 * (PsiNBathC + PsiNBathV + 2)) != Psi.NFermion)
    {
      printf("The total number of Fermions is incompatible witht the bath sizes\n");
      printf("This wave-function can not be used to start from\n");
      return 1;
    }
  }
  else
  {
    printf("\nWave-function file: '%s' does not excist\n Started from new function\n", PsiFileName);
    CreateWaveFunctionDMFTSecond(&Psi, HSet.NBath, HSet.NBath, HSet.NBitsKey);
  }
  strcpy(GBath.Name, "Gbath");
  /* do the DMFT loops */
  for (iOuter = HSet.iStart; iOuter < HSet.iStart + HSet.Ni; iOuter++)
  {
    BeginTiming("DMFTLoop", &Time);
    printf("\n\n***********************************************************************\n");
    printf("Start of loop %i from %i to %i\n", iOuter + 1, HSet.iStart, HSet.iStart + HSet.Ni);
    printf("***********************************************************************\n\n\n");
    /* Create splitted anderson matrix */
    TriDiagonalMatrixToSecondDoubleTriDiagonalMatrix(GBath, &T, &GAnd, Ef, epsilonAnd);
    /* Create GAnd and GAndTrunk as a copy pointing to the same memory */
    GAndTrunk.Ef = GAnd.Ef;
    GAndTrunk.Epsilon = GAnd.Epsilon;
    GAndTrunk.a0 = GAnd.a0;
    GAndTrunk.b0 = GAnd.b0;
    GAndTrunk.Dt.Ef = GAnd.Dt.Ef;
    GAndTrunk.Dt.Epsilon = GAnd.Dt.Epsilon;
    GAndTrunk.Dt.a0 = GAnd.Dt.a0;
    GAndTrunk.Dt.b0 = GAnd.Dt.b0;
    GAndTrunk.Dt.val.N = GAnd.Dt.val.N;
    GAndTrunk.Dt.val.a = GAnd.Dt.val.a;
    GAndTrunk.Dt.val.b = GAnd.Dt.val.b;
    GAndTrunk.Dt.con.N = GAnd.Dt.con.N;
    GAndTrunk.Dt.con.a = GAnd.Dt.con.a;
    GAndTrunk.Dt.con.b = GAnd.Dt.con.b;
    /* Create VPerturbation */
    for (i = 0; i < 4; i++)
    {
      InitOperator(&(VPerturbation[i]));
      VPerturbation[i].NFermion = 4 + 2 * (Min(GAnd.Dt.val.N, HSet.NBath) + Min(GAnd.Dt.con.N, HSet.NBath));
      VPerturbation[i].NBoson = 0;
      VPerturbation[i].NWord = Min(VPerturbation[i].NFermion, 1 + (VPerturbation[i].NFermion - 1) / 8);
      VPerturbation[i].QComplex = 0;
      InitOperatorFixedLength(&(VPerturbation[i]), 1, 1);
    }
    if (GAnd.Dt.val.N > HSet.NBath)
    {
      oadd[0] = CrAnMask + 3 + HSet.NBath;
      RealOperatorAddLader(&(VPerturbation[0]), oadd, 1, GAnd.Dt.val.b[HSet.NBath]);
      oadd[0] = CrAnMask + 3 + Min(GAnd.Dt.con.N, HSet.NBath) + 2 * HSet.NBath;
      RealOperatorAddLader(&(VPerturbation[1]), oadd, 1, GAnd.Dt.val.b[HSet.NBath]);
    }
    if (GAnd.Dt.con.N > HSet.NBath)
    {
      oadd[0] = 3 + Min(GAnd.Dt.val.N, HSet.NBath) + HSet.NBath;
      RealOperatorAddLader(&(VPerturbation[2]), oadd, 1, GAnd.Dt.con.b[HSet.NBath]);
      oadd[0] = 3 + 2 * Min(GAnd.Dt.val.N, HSet.NBath) + 2 * HSet.NBath;
      RealOperatorAddLader(&(VPerturbation[3]), oadd, 1, GAnd.Dt.con.b[HSet.NBath]);
    }
    /* Create gPerturbation */
    if (GAnd.Dt.val.N > HSet.NBath)
    {
      gPerturbation[0].N = GAnd.Dt.val.N - HSet.NBath;
      gPerturbation[1].N = GAnd.Dt.val.N - HSet.NBath;
      InitTriDiagonalMatrix(&(gPerturbation[0]));
      InitTriDiagonalMatrix(&(gPerturbation[1]));
      for (i = 1; i < GAnd.Dt.val.N - HSet.NBath; i++)
      {
        gPerturbation[0].a[i] = -GAnd.Dt.val.a[i + HSet.NBath];
        gPerturbation[0].b[i] = GAnd.Dt.val.b[i + HSet.NBath];
        gPerturbation[1].a[i] = -GAnd.Dt.val.a[i + HSet.NBath];
        gPerturbation[1].b[i] = GAnd.Dt.val.b[i + HSet.NBath];
      }
      gPerturbation[0].a[GAnd.Dt.val.N - HSet.NBath] = -GAnd.Dt.val.a[GAnd.Dt.val.N];
      gPerturbation[1].a[GAnd.Dt.val.N - HSet.NBath] = -GAnd.Dt.val.a[GAnd.Dt.val.N];
    }
    else
    {
      gPerturbation[0].N = 1;
      gPerturbation[1].N = 1;
      InitTriDiagonalMatrix(&(gPerturbation[0]));
      InitTriDiagonalMatrix(&(gPerturbation[1]));
      gPerturbation[0].a[1] = 0;
      gPerturbation[1].a[1] = 0;
    }
    if (GAnd.Dt.con.N > HSet.NBath)
    {
      gPerturbation[2].N = GAnd.Dt.con.N - HSet.NBath;
      gPerturbation[3].N = GAnd.Dt.con.N - HSet.NBath;
      InitTriDiagonalMatrix(&(gPerturbation[2]));
      InitTriDiagonalMatrix(&(gPerturbation[3]));
      for (i = 1; i < GAnd.Dt.con.N - HSet.NBath; i++)
      {
        gPerturbation[2].a[i] = GAnd.Dt.con.a[i + HSet.NBath];
        gPerturbation[2].b[i] = GAnd.Dt.con.b[i + HSet.NBath];
        gPerturbation[3].a[i] = GAnd.Dt.con.a[i + HSet.NBath];
        gPerturbation[3].b[i] = GAnd.Dt.con.b[i + HSet.NBath];
      }
      gPerturbation[2].a[GAnd.Dt.con.N - HSet.NBath] = GAnd.Dt.con.a[GAnd.Dt.con.N];
      gPerturbation[3].a[GAnd.Dt.con.N - HSet.NBath] = GAnd.Dt.con.a[GAnd.Dt.con.N];
    }
    else
    {
      gPerturbation[2].N = 1;
      gPerturbation[3].N = 1;
      InitTriDiagonalMatrix(&(gPerturbation[2]));
      InitTriDiagonalMatrix(&(gPerturbation[3]));
      gPerturbation[2].a[1] = 0;
      gPerturbation[3].a[1] = 0;
    };
    /* Create HAnd */
    FreeMatrix(&T);
    if (GAnd.Dt.con.N < HSet.NBath)
    {
      printf("Warning: GAnd.Dt.con.N(%i)<HSet.NBath(%i)\n", GAnd.Dt.con.N, HSet.NBath);
    }
    else
    {
      GAndTrunk.Dt.con.N = HSet.NBath;
    }
    if (GAnd.Dt.val.N < HSet.NBath)
    {
      printf("Warning: GAnd.Dt.val.N(%i)<HSet.NBath(%i)\n", GAnd.Dt.val.N, HSet.NBath);
    }
    else
    {
      GAndTrunk.Dt.val.N = HSet.NBath;
    }
    mu = -HSet.U / 2;
    SecondDoubleTriDiagonalMatrixToAndersonImpurityOperator(&HAnd, GAndTrunk, GAndTrunk, HSet.U, mu);
    /* test if psi has correct dimensions and correct if needed */
    sscanf(Psi.Name, "%s %i %i", Buffer, &PsiNBathC, &PsiNBathV);
    if (PsiNBathC != GAndTrunk.Dt.con.N)
    {
      printf("\nBefore diagonalizing:\nThe wavefunction has %i conduction bath sites\nThe Anderson Hamiltonian %i\n",
             PsiNBathC, GAndTrunk.Dt.con.N);
      if (PsiNBathC > GAndTrunk.Dt.con.N)
      {
        printf("Removing %i unoccupied fermions from the wave-function\n", PsiNBathC - GAndTrunk.Dt.con.N);
        for (i = 0; i < PsiNBathC - GAndTrunk.Dt.con.N; i++)
        {
          WaveFunctionRemoveFermion(Psi, &PsiTemp, Psi.NFermion - 1);
          WaveFunctionFree(&Psi);
          WaveFunctionRemoveFermion(PsiTemp, &Psi, 4 + PsiNBathV + PsiNBathC - i - 1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      else
      {
        printf("Adding %i unoccupied fermions to the wave-function\n", GAndTrunk.Dt.con.N - PsiNBathC);
        for (i = 0; i < GAndTrunk.Dt.con.N - PsiNBathC; i++)
        {
          WaveFunctionAddFermion(Psi, &PsiTemp, Psi.NFermion, 0);
          WaveFunctionFree(&Psi);
          WaveFunctionAddFermion(PsiTemp, &Psi, 4 + PsiNBathV + PsiNBathC + i, 0);
          WaveFunctionFree(&PsiTemp);
        }
      }
    }
    if (PsiNBathV != GAndTrunk.Dt.val.N)
    {
      printf("\nBefore diagonalizing:\nThe wavefunction has %i valence bath sites\nThe Anderson Hamiltonian %i\n",
             PsiNBathV, GAndTrunk.Dt.val.N);
      if (PsiNBathV > GAndTrunk.Dt.val.N)
      {
        printf("Removing %i occupied fermions from the wave-function\n", PsiNBathV - GAndTrunk.Dt.val.N);
        for (i = 0; i < PsiNBathV - GAndTrunk.Dt.val.N; i++)
        {
          WaveFunctionRemoveFermion(Psi, &PsiTemp, Psi.NFermion - GAndTrunk.Dt.con.N - 1);
          WaveFunctionFree(&Psi);
          WaveFunctionRemoveFermion(PsiTemp, &Psi, 4 + PsiNBathV - i - 1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      else
      {
        printf("Adding %i occupied fermions to the wave-function\n", GAndTrunk.Dt.val.N - PsiNBathV);
        for (i = 0; i < GAndTrunk.Dt.val.N - PsiNBathV; i++)
        {
          WaveFunctionAddFermion(Psi, &PsiTemp, Psi.NFermion - GAndTrunk.Dt.val.N, 1);
          WaveFunctionFree(&Psi);
          WaveFunctionAddFermion(PsiTemp, &Psi, 4 + PsiNBathV + i, 1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      printf("Corrected basis\n");
      PrintWaveFunctionReadableTruncated(Psi, 0.1);
    }
    sprintf(Psi.Name, "DMFT %i %i", GAndTrunk.Dt.con.N, GAndTrunk.Dt.val.N);
    /* create ground-state wave-function */
    BeginTiming("Lanczos GroundState", &Time);
    LanczosGroundState(HAnd, &Psi, &HSet.zeroLancWavefunction, HSet.zeroLancEnergy, HSet.U, &Psi, 0);
    LanczosGroundState(HAnd, &Psi, &HSet.zeroLancWavefunction, HSet.zeroLancEnergy, HSet.U, &Psi, 0);
    EndTiming("Lanczos GroundState", &Time);
    PsiHermitianOperatorPsi(HAnd, Psi, &PsiEnergy);
    PsiHermitianOperatorSqrPsi(HAnd, Psi, &PsiQuality);
    PsiQuality /= Sqr(PsiEnergy);
    PrintWaveFunctionReadableTruncated(Psi, 0.1);
    CalculatePrefactorHistogram(PsiHist, Psi, 10.0, 20);
    printf("\nHistogram of prefactors in Wave-function\n");
    for (i = 0; i < 19; i++)
    {
      printf("10^(-%2i) to 10^(-%2i) | %22.15E\n", i, i + 1, PsiHist[i]);
    }
    printf("10^(-%2i) to 0        | %22.15E\n", i, PsiHist[19]);
    printf("\n");
    printf("Energy                  <psi|H^2|psi>/(<psi|H|psi>^2)\n");
    printf("%22.15E  %22.15E\n", PsiEnergy, PsiQuality);
    PsiRho = calloc(Psi.NFermion, sizeof(double));
    if (PsiRho == NULL)
    {
      printf("calloc failed in main for PsiRho\n Fatal Error\n");
      fflush(stdout);
      return 1;
    }
    CalculateDensity(PsiRho, Psi);
    printf("Occupation \n");
    for (i = 0; i < Psi.NFermion; i++)
    {
      printf("%5i %22.15E\n", i, PsiRho[i]);
    }
    free(PsiRho);
    /*                                           */
    /* create spin up and spin down PES and IPES */
    /*                                           */
    CreatePESdnOperatorDMFT(&OPESdn, HAnd.NFermion);
    CreatePESupOperatorDMFT(&OPESup, HAnd.NFermion);
    CreateIPESdnOperatorDMFT(&OIPESdn, HAnd.NFermion);
    CreateIPESupOperatorDMFT(&OIPESup, HAnd.NFermion);
    GclMup.N = Max(HSet.NTriDiag, G0.N);
    InitTriDiagonalMatrix(&GclMup);
    GclMdn.N = Max(HSet.NTriDiag, G0.N);
    InitTriDiagonalMatrix(&GclMdn);
    GclPup.N = Max(HSet.NTriDiag, G0.N);
    InitTriDiagonalMatrix(&GclPup);
    GclPdn.N = Max(HSet.NTriDiag, G0.N);
    InitTriDiagonalMatrix(&GclPdn);
    BeginTiming("Lanczos Tri Diagonalization Dyson", &Time);
    OperatorPsi(OPESup, Psi, &OPsi);
    LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation, gPerturbation,
                                     4, &GclMup);
    EndTiming("Lanczos Tri Diagonalization Dyson", &Time);
    BeginTiming("Lanczos Tri Diagonalization Dyson", &Time);
    OperatorPsi(OPESdn, Psi, &OPsi);
    LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation, gPerturbation,
                                     4, &GclMdn);
    EndTiming("Lanczos Tri Diagonalization Dyson", &Time);
    BeginTiming("Lanczos Tri Diagonalization Dyson", &Time);
    OperatorPsi(OIPESup, Psi, &OPsi);
    LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation, gPerturbation,
                                     4, &GclPup);
    EndTiming("Lanczos Tri Diagonalization Dyson", &Time);
    BeginTiming("Lanczos Tri Diagonalization Dyson", &Time);
    OperatorPsi(OIPESdn, Psi, &OPsi);
    LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation, gPerturbation,
                                     4, &GclPdn);
    EndTiming("Lanczos Tri Diagonalization Dyson", &Time);
    /* G math - create Gcl as the sum of PES - IPES */
    printf("=== Start G math ===\n");
    BeginTiming("G math", &Time);
    ShiftTriDiagonalMatrix(&GclMdn, -PsiEnergy);
    ShiftTriDiagonalMatrix(&GclMup, -PsiEnergy);
    ShiftTriDiagonalMatrix(&GclPdn, -PsiEnergy);
    ShiftTriDiagonalMatrix(&GclPup, -PsiEnergy);
    for (i = 1; i <= GclMdn.N; i++)
    {
      GclMdn.a[i] *= -1;
    }
    for (i = 1; i <= GclMup.N; i++)
    {
      GclMup.a[i] *= -1;
    }
    strcpy(GclPup.Name, "Gcluster-Plus-Up");
    strcpy(GclPdn.Name, "Gcluster-Plus-Dn");
    strcpy(GclMup.Name, "Gcluster-Min-Up");
    strcpy(GclMdn.Name, "Gcluster-Min-Dn");
    printf("Start SumGreensfunctionsTriDiagonalRepresentation\n");
    fflush(stdout);
    SumGreensfunctionsTriDiagonalRepresentation(GclMup, GclPup, &Gclup);
    strcpy(Gclup.Name, "Gcluster-Up");
    SumGreensfunctionsTriDiagonalRepresentation(GclMdn, GclPdn, &Gcldn);
    strcpy(Gcldn.Name, "Gcluster-Dn");
    SumGreensfunctionsTriDiagonalRepresentation(Gcldn, Gclup, &Gcl);
    strcpy(Gcl.Name, "Gcluster");
    Gcl.b[0] = 1.0;
    ReducePolesTriDiagonalMatrix(&Gcl, Max(4 * HSet.NTriDiag, 6 * G0.N), 0.0);
    /* symmetrize */
    SymmetrizeTridiagonalMatrix(&Gcl);
    /* G math - create Self energy From Gcl and GBath */
    BeginTiming("G math Selfenergy", &Time);
    printf("Start SelfEnergyTriDiagonalRepresentation\n");
    SelfEnergyTriDiagonalRepresentation(GBath, Gcl, &Sigma);
    strcpy(Sigma.Name, "Sigma");
    EndTiming("G math Selfenergy", &Time);


    printf("Start GBathTriDiagonalRepresentation\n");
    GBathTriDiagonalRepresentation(G0, Sigma, &GBathNew, Max(4 * HSet.NTriDiag, 8 * G0.N));
    strcpy(GBathNew.Name, "Gbath-New");

    /* calculate GBath for the next loop by mixing GBath with GBathNew */
    printf("  Start mixing GBath with GBathNew\n");
    TriDiagonalMatrixCopy(&GBathOld, GBath);
    strcpy(GBathOld.Name, "Gbath-Old");
    FreeTriDiagonalMatrix(&GBath);
    GBathOld.b[0] = sqrt(HSet.Mix);
    GBathNew.b[0] = sqrt(1 - HSet.Mix);
    SumGreensfunctionsTriDiagonalRepresentation(GBathNew, GBathOld, &GBath);
    strcpy(GBath.Name, "Gbath");
    GBathOld.b[0] = 1;
    GBathNew.b[0] = 1;
    ReducePolesTriDiagonalMatrix(&GBath, Max(4 * HSet.NTriDiag, 2 * G0.N), 0.0);
    /* symmetrize */
    SymmetrizeTridiagonalMatrix(&GBath);
    EndTiming("G math", &Time);
    /* create spectra */
    printf("Start create spectra\n");
    Spectra[0].N = HSet.NE;
    Spectra[0].Emin = HSet.Emin;
    Spectra[0].Emax = HSet.Emax;
    Spectra[0].dE = (HSet.Emax - HSet.Emin) / HSet.NE;
    Spectra[0].EminPole = HSet.Emin;
    Spectra[0].EmaxPole = HSet.Emax;
    Spectra[0].Gamma = HSet.Gamma;
    for (i = 0; i < 10; i++)
    {
      SpectrumAllocate(Spectra[0], &Spectra[i]);
    }
    TriDiagonalMatrixToSpectrum(0, GBath, &(Spectra[0])); /*  3 */
    TriDiagonalMatrixToSpectrum(0, GBathOld, &(Spectra[1]));  /*  5 */
    TriDiagonalMatrixToSpectrum(0, GBathNew, &(Spectra[2]));  /*  7 */
    SecondDoubleTriDiagonalMatrixToSpectrum(0, GAnd, &(Spectra[3]));  /*  9 */
    TriDiagonalMatrixToSpectrum(0, G0, &(Spectra[4]));  /* 11 */
    TriDiagonalMatrixToSpectrum(0, Sigma, &(Spectra[5])); /* 13 */
    //    DoubleTriDiagonalMatrixToSpectrum(0,SigmaDouble,&(Spectra[6]));  /* 15 */
    TriDiagonalMatrixToSpectrum(0, GclPup, &(Spectra[7]));  /* 17 */
    TriDiagonalMatrixToSpectrum(0, GclMup, &(Spectra[8]));  /* 19 */
    TriDiagonalMatrixToSpectrum(0, Gcl, &(Spectra[9])); /* 21 */
    /* write to disk */
    printf("Start write to dist\n");
    if (iOuter == 0)
    {
      sprintf(GBathFileName, "%s001", HSet.GRootFileName);
      sprintf(PsiFileName, "%s001", HSet.WaveFunctionRootFileName);
      sprintf(SpecFileName, "%s001", HSet.SpecRootFileName);
    }
    else if (iOuter < 9)
    {
      sprintf(GBathFileName, "%s00%1i", HSet.GRootFileName, iOuter + 1);
      sprintf(PsiFileName, "%s00%1i", HSet.WaveFunctionRootFileName, iOuter + 1);
      sprintf(SpecFileName, "%s00%1i", HSet.SpecRootFileName, iOuter + 1);
    }
    else if (iOuter < 99)
    {
      sprintf(GBathFileName, "%s0%2i", HSet.GRootFileName, iOuter + 1);
      sprintf(PsiFileName, "%s0%2i", HSet.WaveFunctionRootFileName, iOuter + 1);
      sprintf(SpecFileName, "%s0%2i", HSet.SpecRootFileName, iOuter + 1);
    }
    else if (iOuter < 999)
    {
      sprintf(GBathFileName, "%s%3i", HSet.GRootFileName, iOuter + 1);
      sprintf(PsiFileName, "%s%3i", HSet.WaveFunctionRootFileName, iOuter + 1);
      sprintf(SpecFileName, "%s%3i", HSet.SpecRootFileName, iOuter + 1);
    }
    FilePrintSpectra(Spectra, 10, SpecFileName, "w");
    FilePrintPlainTriDiagonalMatrix(GBath, GBathFileName, "w");
    FilePrintWaveFunctionPlainBinary(Psi, PsiFileName, "w");
    /* free variables */
    printf("Start free variables\n");
    for (i = 0; i < 4; i++)
    {
      FreeOperator(&(VPerturbation[i]));
      FreeTriDiagonalMatrix(&(gPerturbation[i]));
    }
    FreeOperator(&HAnd);
    FreeOperator(&OPESdn);
    FreeOperator(&OPESup);
    FreeOperator(&OIPESdn);
    FreeOperator(&OIPESup);
    FreeTriDiagonalMatrix(&GclMup);
    FreeTriDiagonalMatrix(&GclMdn);
    FreeTriDiagonalMatrix(&GclPup);
    FreeTriDiagonalMatrix(&GclPdn);
    FreeTriDiagonalMatrix(&Gcl);
    FreeTriDiagonalMatrix(&Sigma);
    FreeTriDiagonalMatrix(&GBathNew);
    FreeTriDiagonalMatrix(&GBathOld);
    FreeDoubleTriDiagonalMatrix(&(GAnd.Dt));
    EndTiming("DMFTLoop", &Time);

    PrintTiming(Time);
  }



  EndTiming("Total", &Time);
  PrintTiming(Time);
  return 0;
}



/*   GLocal(w) = G0( w - Sigma(w) )
 *             = 1/( w - H0 - Sigma(w) )
 *
 *   GLocal(w)^{-1} = w - H0 - Sigma(w)
 *                  = - ( -G0^{-1} + Sigma(w) )
 *
 */
uint32_t GLocalTriDiagonalRepresentationOld2(TriDiagonalMatrixType const G0, TriDiagonalMatrixType const Sigma,
                                             TriDiagonalMatrixType * const G, uint32_t NTriDiag)
{
  TriDiagonalMatrixType G0Inv;
  uint32_t i;

  G0Inv.N = G0.N - 1;
  G0Inv.a = &(G0.a[1]);
  G0Inv.b = &(G0.b[1]);
  if (SumGreensfunctionsTriDiagonalRepresentation(G0Inv, Sigma, G))
  {
    printf("SumGreensfunctionsTriDiagonalRepresentation failed in GLocalTriDiagonalRepresentation\n");
    return 1;
  };
  G->a[G->N] = G->a[G->N - 1];
  for (i = G->N - 1; i > 0; i--)
  {
    G->a[i] = G->a[i - 1];
    G->b[i] = G->b[i - 1];
  }
  G->a[0] = G0.a[0];
  G->b[0] = G0.b[0];

  if (NTriDiag > (G->N))
  {
    printf("Warning GLocalTriDiagonalRepresentation called with NTriDiag=%i, but G0.N=%i and Sigma.N=%i\n", NTriDiag,
           G0.N, Sigma.N);
    printf("Returned Glocal with Glocal.N=%i\n", (G->N));
  }
  else
  {
    G->N = NTriDiag;
  }
  return 0;
}


/*   GLocal(w) = G0( w - Sigma(w) )
 *
 *             = 1/( w - H0 - Sigma(w) )
 *
 *   Sigma(w) = sa0 + (sb0^2)/(w - sa1 - sb1^2 (w-sa2 ...
 *   G0(w)    = 0   +    1   /(w - ga1 - gb1^2 (w-ga2 ...
 *
 */
uint32_t GLocalTriDiagonalRepresentationOld(TriDiagonalMatrixType const G0, TriDiagonalMatrixType const Sigma,
                                            TriDiagonalMatrixType * const G, uint32_t NTriDiag)
{
  uint32_t i, j;
  MatrixType T;
  SparseMatrixType GS;
  GS.Ni = G0.N * (Sigma.N + 1);
  GS.Nj = G0.N * (Sigma.N + 1);
  GS.QComplex = 0;
  GS.NBitsKey = 8 * ((uint32_t) (log((double)GS.Ni) / (8 * log(2.0))) + 1);
  if (InitSparseMatrix(&GS))
  {
    printf("InitSparseMatrix failed in GLocalTriDiagonalRepresentation\n");
    fflush(stdout);
    return 1;
  };
  GS.a0 = G0.a[0];
  GS.b0 = G0.b[0];
  for (i = 0; i < G0.N; i++)
  {
    RealSparseMatrixAddElement(&GS, i * (Sigma.N + 1), i * (Sigma.N + 1), G0.a[i + 1] - Sigma.a[0]);
    for (j = 1; j < Sigma.N + 1; j++)
    {
      RealSparseMatrixAddElement(&GS, i * (Sigma.N + 1) + j, i * (Sigma.N + 1) + j, Sigma.a[j]);
      RealSparseMatrixAddElement(&GS, i * (Sigma.N + 1) + j - 1, i * (Sigma.N + 1) + j, Sigma.b[j - 1]);
      RealSparseMatrixAddElement(&GS, i * (Sigma.N + 1) + j, i * (Sigma.N + 1) + j - 1, Sigma.b[j - 1]);
    }
  }
  for (i = 1; i < G0.N; i++)
  {
    RealSparseMatrixAddElement(&GS, (i - 1) * (Sigma.N + 1), (i) * (Sigma.N + 1), G0.b[i]);
    RealSparseMatrixAddElement(&GS, (i) * (Sigma.N + 1), (i - 1) * (Sigma.N + 1), G0.b[i]);
  }
  if (SparseMatrixToTriDiagonalMatrix(GS, &T, G, NTriDiag))
  {
    printf("SparseMatrixToTriDiagonalMatrix failed in GLocalTriDiagonalRepresentation\n");
    fflush(stdout);
    return 1;
  };
  FreeMatrix(&T);
  FreeSparseMatrix(&GS);
  return 0;
}

uint32_t OldGBathTriDiagonalRepresentation(TriDiagonalMatrixType const G0, TriDiagonalMatrixType const Sigma,
                                           TriDiagonalMatrixType * const GBath)
{
  uint32_t i;
  TriDiagonalMatrixType G0Inv, G;
  G0Inv.N = G0.N - 2;
  G0Inv.a = &(G0.a[2]);
  G0Inv.b = &(G0.b[2]);
  if (SumGreensfunctionsTriDiagonalRepresentation(G0Inv, Sigma, &G))
  {
    printf("SumGreensfunctionsTriDiagonalRepresentation failed in GLocalTriDiagonalRepresentation\n");
    return 1;
  };
  GBath->N = G.N + 2;
  InitTriDiagonalMatrix(GBath);
  GBath->a[0] = 0;
  GBath->b[0] = 1;
  GBath->a[1] = G0.a[1];
  GBath->b[1] = G0.b[1];
  for (i = 2; i < GBath->N; i++)
  {
    GBath->a[i] = G.a[i - 2];
    GBath->b[i] = G.b[i - 2];
  }
  GBath->a[GBath->N] = G.a[GBath->N - 2];
  FreeTriDiagonalMatrix(&G);
  return 0;
}




int oldmainDyson(DMFTSettings HSet, TimeType * Time)
{
  WaveFunctionType Psi, OPsi, PsiTemp;
  TriDiagonalMatrixType G0, GBath, GclPup, GclPdn, GclMup, GclMdn, Gclup, Gcldn, Gcl, Sigma, GBathNew, GBathOld,
    gPerturbation[4];
  DoubleTriDiagonalMatrixType SigmaDouble;
  SecondDoubleTriDiagonalMatrixType GAnd, GAndTrunk;
  OperatorType HAnd, OPESdn, OPESup, OIPESdn, OIPESup, VPerturbation[4];
  MatrixType T;
  SpectraType Spectra[10];
  char   GBathFileName[255], PsiFileName[255], SpecFileName[25], Buffer[255];
  uint32_t iOuter, i, PsiNBathV, PsiNBathC;
  uint16_t oadd[1];
  double Ef = 0.0, epsilonAnd = 0.00001, mu;
  double PsiEnergy, PsiQuality, PsiHist[20], *PsiRho;
  /* read G0 */
  if (FileReadPlainTriDiagonalMatrixType(&G0, HSet.G0FileName))
  {
    printf("Fatal Error in main\n Failed to read G0 from file %s\n", HSet.G0FileName);
    return 1;
  };
  strcpy(G0.Name, "G0");
  /* read G_Bath_iStart and psi_iStart */
  if (HSet.iStart == 0)
  {
    sprintf(GBathFileName, "%s000", HSet.GRootFileName);
    sprintf(PsiFileName, "%s000", HSet.WaveFunctionRootFileName);
  }
  else if (HSet.iStart < 10)
  {
    sprintf(GBathFileName, "%s00%1i", HSet.GRootFileName, HSet.iStart);
    sprintf(PsiFileName, "%s00%1i", HSet.WaveFunctionRootFileName, HSet.iStart);
  }
  else if (HSet.iStart < 100)
  {
    sprintf(GBathFileName, "%s0%2i", HSet.GRootFileName, HSet.iStart);
    sprintf(PsiFileName, "%s0%2i", HSet.WaveFunctionRootFileName, HSet.iStart);
  }
  else
  {
    sprintf(GBathFileName, "%s%i", HSet.GRootFileName, HSet.iStart);
    sprintf(PsiFileName, "%s%i", HSet.WaveFunctionRootFileName, HSet.iStart);
  }
  if (FileReadPlainTriDiagonalMatrixType(&GBath, GBathFileName))
  {
    printf("Fatal Error in main\n Failed to read GBath from file %s !\n", GBathFileName);
    return 1;
  };
  if (file_exists(PsiFileName))
  {
    Psi.NBitsKey = HSet.NBitsKey;
    if (FileReadWaveFunctionPlainBinary(&Psi, PsiFileName))
    {
      printf("Fatal Error in main\n Failed to read Psi from file %s !\n", PsiFileName);
      return 1;
    };
    printf("\nRead the starting wave function from the file %s\n", PsiFileName);
    sscanf(Psi.Name, "%s %i %i", Buffer, &PsiNBathC, &PsiNBathV);
    printf("Wave-function \"%s\" found has %i conduction bath sites, %i valence bath sites and %i Fermions in basis\n",
           Psi.Name, PsiNBathC, PsiNBathV, Psi.NFermion);
    if ((2 * (PsiNBathC + PsiNBathV + 2)) != Psi.NFermion)
    {
      printf("The total number of Fermions is incompatible witht the bath sizes\n");
      printf("This wave-function can not be used to start from\n");
      return 1;
    }
  }
  else
  {
    printf("\nWave-function file: '%s' does not excist\n Started from new function\n", PsiFileName);
    CreateWaveFunctionDMFTSecond(&Psi, HSet.NBath, HSet.NBath, HSet.NBitsKey);
  }
  strcpy(GBath.Name, "Gbath");
  /* do the DMFT loops */
  for (iOuter = HSet.iStart; iOuter < HSet.iStart + HSet.Ni; iOuter++)
  {
    BeginTiming("DMFTLoop", Time);
    printf("\n\n***********************************************************************\n");
    printf("Start of loop %i from %i to %i\n", iOuter + 1, HSet.iStart, HSet.iStart + HSet.Ni);
    printf("***********************************************************************\n\n\n");
    /* Create splitted anderson matrix */
    TriDiagonalMatrixToSecondDoubleTriDiagonalMatrix(GBath, &T, &GAnd, Ef, epsilonAnd);
    /* Create GAnd and GAndTrunk as a copy pointing to the same memory */
    GAndTrunk.Ef = GAnd.Ef;
    GAndTrunk.Epsilon = GAnd.Epsilon;
    GAndTrunk.a0 = GAnd.a0;
    GAndTrunk.b0 = GAnd.b0;
    GAndTrunk.Dt.Ef = GAnd.Dt.Ef;
    GAndTrunk.Dt.Epsilon = GAnd.Dt.Epsilon;
    GAndTrunk.Dt.a0 = GAnd.Dt.a0;
    GAndTrunk.Dt.b0 = GAnd.Dt.b0;
    GAndTrunk.Dt.val.N = GAnd.Dt.val.N;
    GAndTrunk.Dt.val.a = GAnd.Dt.val.a;
    GAndTrunk.Dt.val.b = GAnd.Dt.val.b;
    GAndTrunk.Dt.con.N = GAnd.Dt.con.N;
    GAndTrunk.Dt.con.a = GAnd.Dt.con.a;
    GAndTrunk.Dt.con.b = GAnd.Dt.con.b;
    /* Create VPerturbation */
    for (i = 0; i < 4; i++)
    {
      InitOperator(&(VPerturbation[i]));
      VPerturbation[i].NFermion = 4 + 2 * (Min(GAnd.Dt.val.N, HSet.NBath) + Min(GAnd.Dt.con.N, HSet.NBath));
      VPerturbation[i].NBoson = 0;
      VPerturbation[i].NWord = Min(VPerturbation[i].NFermion, 1 + (VPerturbation[i].NFermion - 1) / 8);
      VPerturbation[i].QComplex = 0;
      InitOperatorFixedLength(&(VPerturbation[i]), 1, 1);
    }
    if (GAnd.Dt.val.N > HSet.NBath)
    {
      oadd[0] = CrAnMask + 3 + HSet.NBath;
      RealOperatorAddLader(&(VPerturbation[0]), oadd, 1, GAnd.Dt.val.b[HSet.NBath]);
      oadd[0] = CrAnMask + 3 + Min(GAnd.Dt.con.N, HSet.NBath) + 2 * HSet.NBath;
      RealOperatorAddLader(&(VPerturbation[1]), oadd, 1, GAnd.Dt.val.b[HSet.NBath]);
    }
    if (GAnd.Dt.con.N > HSet.NBath)
    {
      oadd[0] = 3 + Min(GAnd.Dt.val.N, HSet.NBath) + HSet.NBath;
      RealOperatorAddLader(&(VPerturbation[2]), oadd, 1, GAnd.Dt.con.b[HSet.NBath]);
      oadd[0] = 3 + 2 * Min(GAnd.Dt.val.N, HSet.NBath) + 2 * HSet.NBath;
      RealOperatorAddLader(&(VPerturbation[3]), oadd, 1, GAnd.Dt.con.b[HSet.NBath]);
    }
    /* Create gPerturbation */
    if (GAnd.Dt.val.N > HSet.NBath)
    {
      gPerturbation[0].N = GAnd.Dt.val.N - HSet.NBath;
      gPerturbation[1].N = GAnd.Dt.val.N - HSet.NBath;
      InitTriDiagonalMatrix(&(gPerturbation[0]));
      InitTriDiagonalMatrix(&(gPerturbation[1]));
      for (i = 1; i < GAnd.Dt.val.N - HSet.NBath; i++)
      {
        gPerturbation[0].a[i] = -GAnd.Dt.val.a[i + HSet.NBath];
        gPerturbation[0].b[i] = GAnd.Dt.val.b[i + HSet.NBath];
        gPerturbation[1].a[i] = -GAnd.Dt.val.a[i + HSet.NBath];
        gPerturbation[1].b[i] = GAnd.Dt.val.b[i + HSet.NBath];
      }
      gPerturbation[0].a[GAnd.Dt.val.N - HSet.NBath] = -GAnd.Dt.val.a[GAnd.Dt.val.N];
      gPerturbation[1].a[GAnd.Dt.val.N - HSet.NBath] = -GAnd.Dt.val.a[GAnd.Dt.val.N];
    }
    else
    {
      gPerturbation[0].N = 1;
      gPerturbation[1].N = 1;
      InitTriDiagonalMatrix(&(gPerturbation[0]));
      InitTriDiagonalMatrix(&(gPerturbation[1]));
      gPerturbation[0].a[1] = 0;
      gPerturbation[1].a[1] = 0;
    }
    if (GAnd.Dt.con.N > HSet.NBath)
    {
      gPerturbation[2].N = GAnd.Dt.con.N - HSet.NBath;
      gPerturbation[3].N = GAnd.Dt.con.N - HSet.NBath;
      InitTriDiagonalMatrix(&(gPerturbation[2]));
      InitTriDiagonalMatrix(&(gPerturbation[3]));
      for (i = 1; i < GAnd.Dt.con.N - HSet.NBath; i++)
      {
        gPerturbation[2].a[i] = GAnd.Dt.con.a[i + HSet.NBath];
        gPerturbation[2].b[i] = GAnd.Dt.con.b[i + HSet.NBath];
        gPerturbation[3].a[i] = GAnd.Dt.con.a[i + HSet.NBath];
        gPerturbation[3].b[i] = GAnd.Dt.con.b[i + HSet.NBath];
      }
      gPerturbation[2].a[GAnd.Dt.con.N - HSet.NBath] = GAnd.Dt.con.a[GAnd.Dt.con.N];
      gPerturbation[3].a[GAnd.Dt.con.N - HSet.NBath] = GAnd.Dt.con.a[GAnd.Dt.con.N];
    }
    else
    {
      gPerturbation[2].N = 1;
      gPerturbation[3].N = 1;
      InitTriDiagonalMatrix(&(gPerturbation[2]));
      InitTriDiagonalMatrix(&(gPerturbation[3]));
      gPerturbation[2].a[1] = 0;
      gPerturbation[3].a[1] = 0;
    };
    /* Create HAnd */
    FreeMatrix(&T);
    if (GAnd.Dt.con.N < HSet.NBath)
    {
      printf("Warning: GAnd.Dt.con.N(%i)<HSet.NBath(%i)\n", GAnd.Dt.con.N, HSet.NBath);
    }
    else
    {
      GAndTrunk.Dt.con.N = HSet.NBath;
    }
    if (GAnd.Dt.val.N < HSet.NBath)
    {
      printf("Warning: GAnd.Dt.val.N(%i)<HSet.NBath(%i)\n", GAnd.Dt.val.N, HSet.NBath);
    }
    else
    {
      GAndTrunk.Dt.val.N = HSet.NBath;
    }
    mu = -HSet.U / 2;
    SecondDoubleTriDiagonalMatrixToAndersonImpurityOperator(&HAnd, GAndTrunk, GAndTrunk, HSet.U, mu);
    /* test if psi has correct dimensions and correct if needed */
    sscanf(Psi.Name, "%s %i %i", Buffer, &PsiNBathC, &PsiNBathV);
    if (PsiNBathC != GAndTrunk.Dt.con.N)
    {
      printf("\nBefore diagonalizing:\nThe wavefunction has %i conduction bath sites\nThe Anderson Hamiltonian %i\n",
             PsiNBathC, GAndTrunk.Dt.con.N);
      if (PsiNBathC > GAndTrunk.Dt.con.N)
      {
        printf("Removing %i unoccupied fermions from the wave-function\n", PsiNBathC - GAndTrunk.Dt.con.N);
        for (i = 0; i < PsiNBathC - GAndTrunk.Dt.con.N; i++)
        {
          WaveFunctionRemoveFermion(Psi, &PsiTemp, Psi.NFermion - 1);
          WaveFunctionFree(&Psi);
          WaveFunctionRemoveFermion(PsiTemp, &Psi, 4 + PsiNBathV + PsiNBathC - i - 1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      else
      {
        printf("Adding %i unoccupied fermions to the wave-function\n", GAndTrunk.Dt.con.N - PsiNBathC);
        for (i = 0; i < GAndTrunk.Dt.con.N - PsiNBathC; i++)
        {
          WaveFunctionAddFermion(Psi, &PsiTemp, Psi.NFermion, 0);
          WaveFunctionFree(&Psi);
          WaveFunctionAddFermion(PsiTemp, &Psi, 4 + PsiNBathV + PsiNBathC + i, 0);
          WaveFunctionFree(&PsiTemp);
        }
      }
    }
    if (PsiNBathV != GAndTrunk.Dt.val.N)
    {
      printf("\nBefore diagonalizing:\nThe wavefunction has %i valence bath sites\nThe Anderson Hamiltonian %i\n",
             PsiNBathV, GAndTrunk.Dt.val.N);
      if (PsiNBathV > GAndTrunk.Dt.val.N)
      {
        printf("Removing %i occupied fermions from the wave-function\n", PsiNBathV - GAndTrunk.Dt.val.N);
        for (i = 0; i < PsiNBathV - GAndTrunk.Dt.val.N; i++)
        {
          WaveFunctionRemoveFermion(Psi, &PsiTemp, Psi.NFermion - GAndTrunk.Dt.con.N - 1);
          WaveFunctionFree(&Psi);
          WaveFunctionRemoveFermion(PsiTemp, &Psi, 4 + PsiNBathV - i - 1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      else
      {
        printf("Adding %i occupied fermions to the wave-function\n", GAndTrunk.Dt.val.N - PsiNBathV);
        for (i = 0; i < GAndTrunk.Dt.val.N - PsiNBathV; i++)
        {
          WaveFunctionAddFermion(Psi, &PsiTemp, Psi.NFermion - GAndTrunk.Dt.val.N, 1);
          WaveFunctionFree(&Psi);
          WaveFunctionAddFermion(PsiTemp, &Psi, 4 + PsiNBathV + i, 1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      printf("Corrected basis\n");
      PrintWaveFunctionReadableTruncated(Psi, 0.1);
    }
    sprintf(Psi.Name, "DMFT %i %i", GAndTrunk.Dt.con.N, GAndTrunk.Dt.val.N);
    /* create ground-state wave-function */
    BeginTiming("Lanczos GroundState", Time);
    LanczosGroundState(HAnd, &Psi, &HSet.zeroLancWavefunction, HSet.zeroLancEnergy, HSet.U, &Psi, 0);
    LanczosGroundState(HAnd, &Psi, &HSet.zeroLancWavefunction, HSet.zeroLancEnergy, HSet.U, &Psi, 0);
    EndTiming("Lanczos GroundState", Time);
    PsiHermitianOperatorPsi(HAnd, Psi, &PsiEnergy);
    PsiHermitianOperatorSqrPsi(HAnd, Psi, &PsiQuality);
    PsiQuality /= Sqr(PsiEnergy);
    PrintWaveFunctionReadableTruncated(Psi, 0.1);
    CalculatePrefactorHistogram(PsiHist, Psi, 10.0, 20);
    printf("\nHistogram of prefactors in Wave-function\n");
    for (i = 0; i < 19; i++)
    {
      printf("10^(-%2i) to 10^(-%2i) | %22.15E\n", i, i + 1, PsiHist[i]);
    }
    printf("10^(-%2i) to 0        | %22.15E\n", i, PsiHist[19]);
    printf("\n");
    printf("Energy                  <psi|H^2|psi>/(<psi|H|psi>^2)\n");
    printf("%22.15E  %22.15E\n", PsiEnergy, PsiQuality);
    PsiRho = calloc(Psi.NFermion, sizeof(double));
    if (PsiRho == NULL)
    {
      printf("calloc failed in main for PsiRho\n Fatal Error\n");
      fflush(stdout);
      return 1;
    }
    CalculateDensity(PsiRho, Psi);
    printf("Occupation \n");
    for (i = 0; i < Psi.NFermion; i++)
    {
      printf("%5i %22.15E\n", i, PsiRho[i]);
    }
    free(PsiRho);
    /*                                           */
    /* create spin up and spin down PES and IPES */
    /*                                           */
    CreatePESdnOperatorDMFT(&OPESdn, HAnd.NFermion);
    CreatePESupOperatorDMFT(&OPESup, HAnd.NFermion);
    CreateIPESdnOperatorDMFT(&OIPESdn, HAnd.NFermion);
    CreateIPESupOperatorDMFT(&OIPESup, HAnd.NFermion);
    GclMup.N = Max(HSet.NTriDiag, G0.N);
    InitTriDiagonalMatrix(&GclMup);
    GclMdn.N = Max(HSet.NTriDiag, G0.N);
    InitTriDiagonalMatrix(&GclMdn);
    GclPup.N = Max(HSet.NTriDiag, G0.N);
    InitTriDiagonalMatrix(&GclPup);
    GclPdn.N = Max(HSet.NTriDiag, G0.N);
    InitTriDiagonalMatrix(&GclPdn);
    BeginTiming("Lanczos Tri Diagonalization Dyson", Time);
    OperatorPsi(OPESup, Psi, &OPsi);
    LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation, gPerturbation,
                                     4, &GclMup);
    EndTiming("Lanczos Tri Diagonalization Dyson", Time);
    BeginTiming("Lanczos Tri Diagonalization Dyson", Time);
    OperatorPsi(OPESdn, Psi, &OPsi);
    LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation, gPerturbation,
                                     4, &GclMdn);
    EndTiming("Lanczos Tri Diagonalization Dyson", Time);
    BeginTiming("Lanczos Tri Diagonalization Dyson", Time);
    OperatorPsi(OIPESup, Psi, &OPsi);
    LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation, gPerturbation,
                                     4, &GclPup);
    EndTiming("Lanczos Tri Diagonalization Dyson", Time);
    BeginTiming("Lanczos Tri Diagonalization Dyson", Time);
    OperatorPsi(OIPESdn, Psi, &OPsi);
    LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation, gPerturbation,
                                     4, &GclPdn);
    EndTiming("Lanczos Tri Diagonalization Dyson", Time);
    /* G math - create Gcl as the sum of PES - IPES */
    printf("Start G math\n");
    BeginTiming("G math", Time);
    ShiftTriDiagonalMatrix(&GclMdn, -PsiEnergy);
    ShiftTriDiagonalMatrix(&GclMup, -PsiEnergy);
    ShiftTriDiagonalMatrix(&GclPdn, -PsiEnergy);
    ShiftTriDiagonalMatrix(&GclPup, -PsiEnergy);
    for (i = 1; i <= GclMdn.N; i++)
    {
      GclMdn.a[i] *= -1;
    }
    for (i = 1; i <= GclMup.N; i++)
    {
      GclMup.a[i] *= -1;
    }
    strcpy(GclPup.Name, "Gcluster-Plus-Up");
    strcpy(GclPdn.Name, "Gcluster-Plus-Dn");
    strcpy(GclMup.Name, "Gcluster-Min-Up");
    strcpy(GclMdn.Name, "Gcluster-Min-Dn");
    printf("  Start SumGreensfunctionsTriDiagonalRepresentation\n");
    fflush(stdout);
    SumGreensfunctionsTriDiagonalRepresentation(GclMup, GclPup, &Gclup);
    strcpy(Gclup.Name, "Gcluster-Up");
    SumGreensfunctionsTriDiagonalRepresentation(GclMdn, GclPdn, &Gcldn);
    strcpy(Gcldn.Name, "Gcluster-Dn");
    SumGreensfunctionsTriDiagonalRepresentation(Gcldn, Gclup, &Gcl);
    strcpy(Gcl.Name, "Gcluster");
    Gcl.b[0] = 1.0;
    printf("Gcl.N = %i ", Gcl.N);
    SymmetrizeTridiagonalMatrix(&Gcl);
    printf("=> %i\n", Gcl.N);
    /* G math - create Self energy From Gcl and GBath */
    BeginTiming("G math Selfenergy", Time);
    printf("  Start SelfEnergyTriDiagonalRepresentation\n");
    printf("G0.N = %i ", G0.N);
    SelfEnergyTriDiagonalRepresentation(GBath, Gcl, &Sigma);
    printf("Sigma.N = %i ", Sigma.N);
    TriDiagonalMatrixToDoubleTriDiagonalMatrix(Sigma, &T, &SigmaDouble, Ef, epsilonAnd);
    FreeMatrix(&T);
    Sigma.N = Min(Sigma.N, Max(4 * HSet.NTriDiag, 2 * G0.N));
    printf("=> %i ", Sigma.N);
    strcpy(Sigma.Name, "Sigma");
    SymmetrizeTridiagonalMatrix(&Sigma);
    printf("=> %i ", Sigma.N);
    Sigma.N = Min(Sigma.N, Max(2 * HSet.NTriDiag, 2 * G0.N));
    printf("=> %i\n", Sigma.N);
    EndTiming("G math Selfenergy", Time);


    printf("  Start GBathNewTriDiagonalRepresentation\n");
    GBathTriDiagonalRepresentation(G0, Sigma, &GBathNew, Max(2 * HSet.NTriDiag, 2 * G0.N));
    GBathNew.N = Min(GBathNew.N, Max(4 * HSet.NTriDiag, 2 * G0.N));
    strcpy(GBathNew.Name, "Gbath-New");
    SymmetrizeTridiagonalMatrix(&GBathNew);
    GBathNew.N = Min(GBathNew.N, Max(4 * HSet.NTriDiag, 2 * G0.N));


    /* calculate GBath for the next loop by mixing GBath with GBathNew */
    printf("Start mixing GBath with GBathNew\n");
    TriDiagonalMatrixCopy(&GBathOld, GBath);
    strcpy(GBathOld.Name, "Gbath-Old");
    FreeTriDiagonalMatrix(&GBath);
    GBathOld.b[0] = sqrt(HSet.Mix);
    GBathNew.b[0] = sqrt(1 - HSet.Mix);
    SumGreensfunctionsTriDiagonalRepresentation(GBathNew, GBathOld, &GBath);
    strcpy(GBath.Name, "Gbath");
    GBathOld.b[0] = 1;
    GBathNew.b[0] = 1;
    SymmetrizeTridiagonalMatrix(&GBath);
    GBath.N = Min(GBath.N, Max(4 * HSet.NTriDiag, 2 * G0.N));
    EndTiming("G math", Time);
    /* create spectra */
    printf("Start create spectra\n");
    Spectra[0].N = HSet.NE;
    Spectra[0].Emin = HSet.Emin;
    Spectra[0].Emax = HSet.Emax;
    Spectra[0].dE = (HSet.Emax - HSet.Emin) / HSet.NE;
    Spectra[0].EminPole = HSet.Emin;
    Spectra[0].EmaxPole = HSet.Emax;
    Spectra[0].Gamma = HSet.Gamma;
    for (i = 0; i < 10; i++)
    {
      SpectrumAllocate(Spectra[0], &Spectra[i]);
    }
    TriDiagonalMatrixToSpectrum(0, GBath, &(Spectra[0])); /*  3 */
    TriDiagonalMatrixToSpectrum(0, GBathOld, &(Spectra[1]));  /*  5 */
    TriDiagonalMatrixToSpectrum(0, GBathNew, &(Spectra[2]));  /*  7 */
    SecondDoubleTriDiagonalMatrixToSpectrum(0, GAnd, &(Spectra[3]));  /*  9 */
    TriDiagonalMatrixToSpectrum(0, G0, &(Spectra[4]));  /* 11 */
    TriDiagonalMatrixToSpectrum(0, Sigma, &(Spectra[5])); /* 13 */
    DoubleTriDiagonalMatrixToSpectrum(0, SigmaDouble, &(Spectra[6])); /* 15 */
    TriDiagonalMatrixToSpectrum(0, GclPup, &(Spectra[7]));  /* 17 */
    TriDiagonalMatrixToSpectrum(0, GclMup, &(Spectra[8]));  /* 19 */
    TriDiagonalMatrixToSpectrum(0, Gcl, &(Spectra[9])); /* 21 */
    /* write to disk */
    printf("Start write to dist\n");
    if (iOuter == 0)
    {
      sprintf(GBathFileName, "%s001", HSet.GRootFileName);
      sprintf(PsiFileName, "%s001", HSet.WaveFunctionRootFileName);
      sprintf(SpecFileName, "%s001", HSet.SpecRootFileName);
    }
    else if (iOuter < 9)
    {
      sprintf(GBathFileName, "%s00%1i", HSet.GRootFileName, iOuter + 1);
      sprintf(PsiFileName, "%s00%1i", HSet.WaveFunctionRootFileName, iOuter + 1);
      sprintf(SpecFileName, "%s00%1i", HSet.SpecRootFileName, iOuter + 1);
    }
    else if (iOuter < 99)
    {
      sprintf(GBathFileName, "%s0%2i", HSet.GRootFileName, iOuter + 1);
      sprintf(PsiFileName, "%s0%2i", HSet.WaveFunctionRootFileName, iOuter + 1);
      sprintf(SpecFileName, "%s0%2i", HSet.SpecRootFileName, iOuter + 1);
    }
    else if (iOuter < 999)
    {
      sprintf(GBathFileName, "%s%3i", HSet.GRootFileName, iOuter + 1);
      sprintf(PsiFileName, "%s%3i", HSet.WaveFunctionRootFileName, iOuter + 1);
      sprintf(SpecFileName, "%s%3i", HSet.SpecRootFileName, iOuter + 1);
    }
    FilePrintSpectra(Spectra, 10, SpecFileName, "w");
    FilePrintPlainTriDiagonalMatrix(GBath, GBathFileName, "w");
    FilePrintWaveFunctionPlainBinary(Psi, PsiFileName, "w");
    /* free variables */
    printf("Start free variables\n");
    for (i = 0; i < 4; i++)
    {
      FreeOperator(&(VPerturbation[i]));
      FreeTriDiagonalMatrix(&(gPerturbation[i]));
    }
    FreeOperator(&HAnd);
    FreeOperator(&OPESdn);
    FreeOperator(&OPESup);
    FreeOperator(&OIPESdn);
    FreeOperator(&OIPESup);
    FreeTriDiagonalMatrix(&GclMup);
    FreeTriDiagonalMatrix(&GclMdn);
    FreeTriDiagonalMatrix(&GclPup);
    FreeTriDiagonalMatrix(&GclPdn);
    FreeTriDiagonalMatrix(&Gcl);
    FreeTriDiagonalMatrix(&Sigma);
    FreeTriDiagonalMatrix(&GBathNew);
    FreeTriDiagonalMatrix(&GBathOld);
    FreeDoubleTriDiagonalMatrix(&(GAnd.Dt));
    EndTiming("DMFTLoop", Time);

    PrintTiming(*Time);
  }
  return 0;
};


int oldmainNoDyson(DMFTSettings HSet, TimeType * Time)
{
  WaveFunctionType Psi, OPsi, PsiTemp;
  TriDiagonalMatrixType G0, GBath, GclPup, GclPdn, GclMup, GclMdn, Gclup, Gcldn, Gcl, Sigma, GBathNew, GBathOld;
  SecondDoubleTriDiagonalMatrixType GAnd;
  OperatorType HAnd, OPESdn, OPESup, OIPESdn, OIPESup;
  MatrixType T;
  SpectraType Spectra[10];
  char   GBathFileName[255], PsiFileName[255], SpecFileName[25], Buffer[255];
  uint32_t iOuter, i, PsiNBathV, PsiNBathC;
  double Ef = 0.0, epsilonAnd = 0.0000001, mu;
  double PsiEnergy, PsiQuality, PsiHist[20], *PsiRho;
  /* read G0 */
  if (FileReadPlainTriDiagonalMatrixType(&G0, HSet.G0FileName))
  {
    printf("Fatal Error in main\n Failed to read G0 from file %s\n", HSet.G0FileName);
    return 1;
  };
  strcpy(G0.Name, "G0");
  /* read G_Bath_iStart and psi_iStart */
  if (HSet.iStart == 0)
  {
    sprintf(GBathFileName, "%s000", HSet.GRootFileName);
    sprintf(PsiFileName, "%s000", HSet.WaveFunctionRootFileName);
  }
  else if (HSet.iStart < 10)
  {
    sprintf(GBathFileName, "%s00%1i", HSet.GRootFileName, HSet.iStart);
    sprintf(PsiFileName, "%s00%1i", HSet.WaveFunctionRootFileName, HSet.iStart);
  }
  else if (HSet.iStart < 100)
  {
    sprintf(GBathFileName, "%s0%2i", HSet.GRootFileName, HSet.iStart);
    sprintf(PsiFileName, "%s0%2i", HSet.WaveFunctionRootFileName, HSet.iStart);
  }
  else
  {
    sprintf(GBathFileName, "%s%i", HSet.GRootFileName, HSet.iStart);
    sprintf(PsiFileName, "%s%i", HSet.WaveFunctionRootFileName, HSet.iStart);
  }
  if (FileReadPlainTriDiagonalMatrixType(&GBath, GBathFileName))
  {
    printf("Fatal Error in main\n Failed to read GBath from file %s !\n", GBathFileName);
    return 1;
  };
  if (file_exists(PsiFileName))
  {
    Psi.NBitsKey = HSet.NBitsKey;
    if (FileReadWaveFunctionPlainBinary(&Psi, PsiFileName))
    {
      printf("Fatal Error in main\n Failed to read Psi from file %s !\n", PsiFileName);
      return 1;
    };
    printf("\nRead the starting wave function from the file %s\n", PsiFileName);
    sscanf(Psi.Name, "%s %i %i", Buffer, &PsiNBathC, &PsiNBathV);
    printf("Wave-function \"%s\" found has %i conduction bath sites, %i valence bath sites and %i Fermions in basis\n",
           Psi.Name, PsiNBathC, PsiNBathV, Psi.NFermion);
    if ((2 * (PsiNBathC + PsiNBathV + 2)) != Psi.NFermion)
    {
      printf("The total number of Fermions is incompatible witht the bath sizes\n");
      printf("This wave-function can not be used to start from\n");
      return 1;
    }
  }
  else
  {
    printf("\nWave-function file: '%s' does not excist\n Started from new function\n", PsiFileName);
    CreateWaveFunctionDMFTSecond(&Psi, HSet.NBath, HSet.NBath, HSet.NBitsKey);
  }
  strcpy(GBath.Name, "Gbath");
  /* do the DMFT loops */
  for (iOuter = HSet.iStart; iOuter < HSet.iStart + HSet.Ni; iOuter++)
  {
    BeginTiming("DMFTLoop", Time);
    printf("\n\n***********************************************************************\n");
    printf("Start of loop %i from %i to %i\n", iOuter + 1, HSet.iStart, HSet.iStart + HSet.Ni);
    printf("***********************************************************************\n\n\n");
    TriDiagonalMatrixToSecondDoubleTriDiagonalMatrix(GBath, &T, &GAnd, Ef, epsilonAnd);

    PrintTriDiagonalMatrix(GBath);
    PrintDoubleTriDiagonalMatrix(GAnd.Dt);

    FreeMatrix(&T);
    if (GAnd.Dt.con.N < HSet.NBath)
    {
      printf("Warning: GAnd.Dt.con.N(%i)<HSet.NBath(%i)\n", GAnd.Dt.con.N, HSet.NBath);
    }
    else
    {
      GAnd.Dt.con.N = HSet.NBath;
    }
    if (GAnd.Dt.val.N < HSet.NBath)
    {
      printf("Warning: GAnd.Dt.val.N(%i)<HSet.NBath(%i)\n", GAnd.Dt.val.N, HSet.NBath);
    }
    else
    {
      GAnd.Dt.val.N = HSet.NBath;
    }
    mu = -HSet.U / 2;
    SecondDoubleTriDiagonalMatrixToAndersonImpurityOperator(&HAnd, GAnd, GAnd, HSet.U, mu);
    PrintOperatorReadable(HAnd);
    /* test if psi has correct dimensions and correct if needed */
    sscanf(Psi.Name, "%s %i %i", Buffer, &PsiNBathC, &PsiNBathV);
    if (PsiNBathC != GAnd.Dt.con.N)
    {
      printf("\nBefore diagonalizing:\nThe wavefunction has %i conduction bath sites\nThe Anderson Hamiltonian %i\n",
             PsiNBathC, GAnd.Dt.con.N);
      if (PsiNBathC > GAnd.Dt.con.N)
      {
        printf("Removing %i unoccupied fermions from the wave-function\n", PsiNBathC - GAnd.Dt.con.N);
        for (i = 0; i < PsiNBathC - GAnd.Dt.con.N; i++)
        {
          WaveFunctionRemoveFermion(Psi, &PsiTemp, Psi.NFermion - 1);
          WaveFunctionFree(&Psi);
          WaveFunctionRemoveFermion(PsiTemp, &Psi, 4 + PsiNBathV + PsiNBathC - i - 1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      else
      {
        printf("Adding %i unoccupied fermions to the wave-function\n", GAnd.Dt.con.N - PsiNBathC);
        for (i = 0; i < GAnd.Dt.con.N - PsiNBathC; i++)
        {
          WaveFunctionAddFermion(Psi, &PsiTemp, Psi.NFermion, 0);
          WaveFunctionFree(&Psi);
          WaveFunctionAddFermion(PsiTemp, &Psi, 4 + PsiNBathV + PsiNBathC + i, 0);
          WaveFunctionFree(&PsiTemp);
        }
      }
    }
    if (PsiNBathV != GAnd.Dt.val.N)
    {
      printf("\nBefore diagonalizing:\nThe wavefunction has %i valence bath sites\nThe Anderson Hamiltonian %i\n",
             PsiNBathV, GAnd.Dt.val.N);
      if (PsiNBathV > GAnd.Dt.val.N)
      {
        printf("Removing %i occupied fermions from the wave-function\n", PsiNBathV - GAnd.Dt.val.N);
        for (i = 0; i < PsiNBathV - GAnd.Dt.val.N; i++)
        {
          WaveFunctionRemoveFermion(Psi, &PsiTemp, Psi.NFermion - GAnd.Dt.con.N - 1);
          WaveFunctionFree(&Psi);
          WaveFunctionRemoveFermion(PsiTemp, &Psi, 4 + PsiNBathV - i - 1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      else
      {
        printf("Adding %i occupied fermions to the wave-function\n", GAnd.Dt.val.N - PsiNBathV);
        for (i = 0; i < GAnd.Dt.val.N - PsiNBathV; i++)
        {
          WaveFunctionAddFermion(Psi, &PsiTemp, Psi.NFermion - GAnd.Dt.val.N, 1);
          WaveFunctionFree(&Psi);
          WaveFunctionAddFermion(PsiTemp, &Psi, 4 + PsiNBathV + i, 1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      printf("Corrected basis\n");
      PrintWaveFunctionReadableTruncated(Psi, 0.1);
    }
    sprintf(Psi.Name, "DMFT %i %i", GAnd.Dt.con.N, GAnd.Dt.val.N);
    /* create ground-state wave-funciton */
    BeginTiming("Lanczos GroundState", Time);
    LanczosGroundState(HAnd, &Psi, &HSet.zeroLancWavefunction, HSet.zeroLancEnergy, HSet.U, &Psi, 0);
    LanczosGroundState(HAnd, &Psi, &HSet.zeroLancWavefunction, HSet.zeroLancEnergy, HSet.U, &Psi, 0);
    EndTiming("Lanczos GroundState", Time);
    PsiHermitianOperatorPsi(HAnd, Psi, &PsiEnergy);
    PsiHermitianOperatorSqrPsi(HAnd, Psi, &PsiQuality);
    PsiQuality /= Sqr(PsiEnergy);
    PrintWaveFunctionReadableTruncated(Psi, 0.1);
    CalculatePrefactorHistogram(PsiHist, Psi, 10.0, 20);
    printf("\nHistogram of prefactors in Wave-function\n");
    for (i = 0; i < 19; i++)
    {
      printf("10^(-%2i) to 10^(-%2i) | %22.15E\n", i, i + 1, PsiHist[i]);
    }
    printf("10^(-%2i) to 0        | %22.15E\n", i, PsiHist[19]);
    printf("\n");
    printf("Energy                  <psi|H^2|psi>/(<psi|H|psi>^2)\n");
    printf("%22.15E  %22.15E\n", PsiEnergy, PsiQuality);
    PsiRho = calloc(Psi.NFermion, sizeof(double));
    if (PsiRho == NULL)
    {
      printf("calloc failed in main for PsiRho\n Fatal Error\n");
      fflush(stdout);
      return 1;
    }
    CalculateDensity(PsiRho, Psi);
    printf("Occupation \n");
    for (i = 0; i < Psi.NFermion; i++)
    {
      printf("%5i %22.15E\n", i, PsiRho[i]);
    }
    free(PsiRho);
    /* create spin up and spin down PES and IPES */
    CreatePESdnOperatorDMFT(&OPESdn, HAnd.NFermion);
    CreatePESupOperatorDMFT(&OPESup, HAnd.NFermion);
    CreateIPESdnOperatorDMFT(&OIPESdn, HAnd.NFermion);
    CreateIPESupOperatorDMFT(&OIPESup, HAnd.NFermion);
    GclMup.N = HSet.NTriDiag;
    InitTriDiagonalMatrix(&GclMup);
    GclMdn.N = HSet.NTriDiag;
    InitTriDiagonalMatrix(&GclMdn);
    GclPup.N = HSet.NTriDiag;
    InitTriDiagonalMatrix(&GclPup);
    GclPdn.N = HSet.NTriDiag;
    InitTriDiagonalMatrix(&GclPdn);
    BeginTiming("Lanczos Tri Diagonalization PES up", Time);
    OperatorPsi(OPESup, Psi, &OPsi);
    LanczosTriDiagonalize(HAnd, OPsi, HSet.zeroLancWavefunction, &GclMup, HSet.Gamma);
    EndTiming("Lanczos Tri Diagonalization PES up", Time);
    BeginTiming("Lanczos Tri Diagonalization PES dn", Time);
    OperatorPsi(OPESdn, Psi, &OPsi);
    LanczosTriDiagonalize(HAnd, OPsi, HSet.zeroLancWavefunction, &GclMdn, HSet.Gamma);
    EndTiming("Lanczos Tri Diagonalization PES dn", Time);
    BeginTiming("Lanczos Tri Diagonalization IPES up", Time);
    OperatorPsi(OIPESup, Psi, &OPsi);
    LanczosTriDiagonalize(HAnd, OPsi, HSet.zeroLancWavefunction, &GclPup, HSet.Gamma);
    EndTiming("Lanczos Tri Diagonalization IPES up", Time);
    BeginTiming("Lanczos Tri Diagonalization IPES dn", Time);
    OperatorPsi(OIPESdn, Psi, &OPsi);
    LanczosTriDiagonalize(HAnd, OPsi, HSet.zeroLancWavefunction, &GclPdn, HSet.Gamma);
    EndTiming("Lanczos Tri Diagonalization IPES dn", Time);
    /* G math */
    printf("Start G math\n");
    ShiftTriDiagonalMatrix(&GclMdn, -PsiEnergy);
    ShiftTriDiagonalMatrix(&GclMup, -PsiEnergy);
    ShiftTriDiagonalMatrix(&GclPdn, -PsiEnergy);
    ShiftTriDiagonalMatrix(&GclPup, -PsiEnergy);
    for (i = 1; i <= GclMdn.N; i++)
    {
      GclMdn.a[i] *= -1;
    }
    for (i = 1; i <= GclMup.N; i++)
    {
      GclMup.a[i] *= -1;
    }
    strcpy(GclPup.Name, "Gcluster-Plus-Up");
    strcpy(GclPdn.Name, "Gcluster-Plus-Dn");
    strcpy(GclMup.Name, "Gcluster-Min-Up");
    strcpy(GclMdn.Name, "Gcluster-Min-Dn");
    printf("  Start SumGreensfunctionsTriDiagonalRepresentation\n");
    fflush(stdout);
    SumGreensfunctionsTriDiagonalRepresentation(GclPup, GclMup, &Gclup);
    strcpy(Gclup.Name, "Gcluster-Up");
    SumGreensfunctionsTriDiagonalRepresentation(GclPdn, GclMdn, &Gcldn);
    strcpy(Gcldn.Name, "Gcluster-Dn");
    SumGreensfunctionsTriDiagonalRepresentation(Gclup, Gcldn, &Gcl);
    strcpy(Gcl.Name, "Gcluster");
    Gcl.b[0] = 1.0;
    // SymmetrizeTridiagonalMatrix(&Gcl);
    Gcl.N = Min(Gcl.N, Max(4 * HSet.NTriDiag, G0.N));
    printf("  Start SelfEnergyTriDiagonalRepresentation\n");
    SelfEnergyTriDiagonalRepresentation(G0, Gcl, &Sigma);
    Sigma.N = Min(Sigma.N, Max(4 * HSet.NTriDiag, G0.N));
    strcpy(Sigma.Name, "Sigma");
    // SymmetrizeTridiagonalMatrix(&Sigma);
    Sigma.N = Min(Sigma.N, Max(4 * HSet.NTriDiag, G0.N));



    printf("  Start GBathNewTriDiagonalRepresentation\n");
    // SumGreensfunctionsTriDiagonalRepresentation(Gcl, Sigma, &GBathNew);
    GBathTriDiagonalRepresentation(G0, Sigma, &GBathNew, Max(4 * HSet.NTriDiag, 4 * G0.N));



    GBathNew.N = Min(GBathNew.N, Max(4 * HSet.NTriDiag, G0.N));
    strcpy(GBathNew.Name, "Gbath-New");
    // SymmetrizeTridiagonalMatrix(&GBathNew);
    GBathNew.N = Min(GBathNew.N, Max(4 * HSet.NTriDiag, G0.N));
    /* calculate GBath for the next loop by mixing GBath with GBathNew */
    printf("Start mixing GBath with GBathNew\n");
    TriDiagonalMatrixCopy(&GBathOld, GBath);
    strcpy(GBathOld.Name, "Gbath-Old");
    FreeTriDiagonalMatrix(&GBath);
    GBathOld.b[0] = sqrt(HSet.Mix);
    GBathNew.b[0] = sqrt(1 - HSet.Mix);
    SumGreensfunctionsTriDiagonalRepresentation(GBathNew, GBathOld, &GBath);
    strcpy(GBath.Name, "Gbath");
    GBathOld.b[0] = 1;
    GBathNew.b[0] = 1;
    //  SymmetrizeTridiagonalMatrix(&GBath);
    GBath.N = Min(GBath.N, Max(4 * HSet.NTriDiag, G0.N));
    /* create spectra */
    printf("Start create spectra\n");
    Spectra[0].N = HSet.NE;
    Spectra[0].Emin = HSet.Emin;
    Spectra[0].Emax = HSet.Emax;
    Spectra[0].dE = (HSet.Emax - HSet.Emin) / HSet.NE;
    Spectra[0].EminPole = HSet.Emin;
    Spectra[0].EmaxPole = HSet.Emax;
    Spectra[0].Gamma = HSet.Gamma;
    for (i = 0; i < 10; i++)
    {
      SpectrumAllocate(Spectra[0], &Spectra[i]);
    }
    TriDiagonalMatrixToSpectrum(0, GBath, &(Spectra[0]));
    TriDiagonalMatrixToSpectrum(0, GBathOld, &(Spectra[1]));
    TriDiagonalMatrixToSpectrum(0, GBathNew, &(Spectra[2]));
    SecondDoubleTriDiagonalMatrixToSpectrum(0, GAnd, &(Spectra[3]));
    TriDiagonalMatrixToSpectrum(0, G0, &(Spectra[4]));
    TriDiagonalMatrixToSpectrum(0, Sigma, &(Spectra[5]));
    //    TriDiagonalMatrixToSpectrum(0,GLocal,&(Spectra[6]));
    TriDiagonalMatrixToSpectrum(0, GclPup, &(Spectra[7]));
    TriDiagonalMatrixToSpectrum(0, GclMup, &(Spectra[8]));
    TriDiagonalMatrixToSpectrum(0, Gcl, &(Spectra[9]));
    /* write to disk */
    printf("Start write to dist\n");
    if (iOuter == 0)
    {
      sprintf(GBathFileName, "%s001", HSet.GRootFileName);
      sprintf(PsiFileName, "%s001", HSet.WaveFunctionRootFileName);
      sprintf(SpecFileName, "%s001", HSet.SpecRootFileName);
    }
    else if (iOuter < 9)
    {
      sprintf(GBathFileName, "%s00%1i", HSet.GRootFileName, iOuter + 1);
      sprintf(PsiFileName, "%s00%1i", HSet.WaveFunctionRootFileName, iOuter + 1);
      sprintf(SpecFileName, "%s00%1i", HSet.SpecRootFileName, iOuter + 1);
    }
    else if (iOuter < 99)
    {
      sprintf(GBathFileName, "%s0%2i", HSet.GRootFileName, iOuter + 1);
      sprintf(PsiFileName, "%s0%2i", HSet.WaveFunctionRootFileName, iOuter + 1);
      sprintf(SpecFileName, "%s0%2i", HSet.SpecRootFileName, iOuter + 1);
    }
    else if (iOuter < 999)
    {
      sprintf(GBathFileName, "%s%3i", HSet.GRootFileName, iOuter + 1);
      sprintf(PsiFileName, "%s%3i", HSet.WaveFunctionRootFileName, iOuter + 1);
      sprintf(SpecFileName, "%s%3i", HSet.SpecRootFileName, iOuter + 1);
    }
    FilePrintSpectra(Spectra, 10, SpecFileName, "w");
    FilePrintPlainTriDiagonalMatrix(GBath, GBathFileName, "w");
    FilePrintWaveFunctionPlainBinary(Psi, PsiFileName, "w");
    /* free variables */
    printf("Start free variables\n");
    FreeOperator(&HAnd);
    FreeOperator(&OPESdn);
    FreeOperator(&OPESup);
    FreeOperator(&OIPESdn);
    FreeOperator(&OIPESup);
    FreeTriDiagonalMatrix(&GclMup);
    FreeTriDiagonalMatrix(&GclMdn);
    FreeTriDiagonalMatrix(&GclPup);
    FreeTriDiagonalMatrix(&GclPdn);
    FreeTriDiagonalMatrix(&Gcl);
    //   FreeTriDiagonalMatrix(&GLocal);
    FreeTriDiagonalMatrix(&Sigma);
    FreeTriDiagonalMatrix(&GBathNew);
    FreeTriDiagonalMatrix(&GBathOld);
    FreeDoubleTriDiagonalMatrix(&(GAnd.Dt));
    EndTiming("DMFTLoop", Time);
    PrintTiming(*Time);
  }
  return 0;
};


int old2mainDyson(DMFTSettings HSet, TimeType * Time)
{
  WaveFunctionType Psi, OPsi, PsiTemp;
  TriDiagonalMatrixType G0, GBath, GclPup, GclPdn, GclMup, GclMdn, Gclup, Gcldn, Gcl, Sigma, GLocal, GBathNew, GBathOld,
    gPerturbation[4];
  DoubleTriDiagonalMatrixType GAnd, GAndTrunk;
  OperatorType HAnd, OPESdn, OPESup, OIPESdn, OIPESup, VPerturbation[4];
  MatrixType T;
  SpectraType Spectra[10];
  char   GBathFileName[255], PsiFileName[255], SpecFileName[25], Buffer[255];
  uint32_t iOuter, i, PsiNBathV, PsiNBathC;
  uint16_t oadd[1];
  double Ef = 0.0, epsilonAnd = 0.1, mu;
  double PsiEnergy, PsiQuality, PsiHist[20], *PsiRho;
  /* read G0 */
  if (FileReadPlainTriDiagonalMatrixType(&G0, HSet.G0FileName))
  {
    printf("Fatal Error in main\n Failed to read G0 from file %s\n", HSet.G0FileName);
    return 1;
  };
  strcpy(G0.Name, "G0");
  /* read G_Bath_iStart and psi_iStart */
  if (HSet.iStart == 0)
  {
    sprintf(GBathFileName, "%s000", HSet.GRootFileName);
    sprintf(PsiFileName, "%s000", HSet.WaveFunctionRootFileName);
  }
  else if (HSet.iStart < 10)
  {
    sprintf(GBathFileName, "%s00%1i", HSet.GRootFileName, HSet.iStart);
    sprintf(PsiFileName, "%s00%1i", HSet.WaveFunctionRootFileName, HSet.iStart);
  }
  else if (HSet.iStart < 100)
  {
    sprintf(GBathFileName, "%s0%2i", HSet.GRootFileName, HSet.iStart);
    sprintf(PsiFileName, "%s0%2i", HSet.WaveFunctionRootFileName, HSet.iStart);
  }
  else
  {
    sprintf(GBathFileName, "%s%i", HSet.GRootFileName, HSet.iStart);
    sprintf(PsiFileName, "%s%i", HSet.WaveFunctionRootFileName, HSet.iStart);
  }
  if (FileReadPlainTriDiagonalMatrixType(&GBath, GBathFileName))
  {
    printf("Fatal Error in main\n Failed to read GBath from file %s !\n", GBathFileName);
    return 1;
  };
  if (file_exists(PsiFileName))
  {
    Psi.NBitsKey = HSet.NBitsKey;
    if (FileReadWaveFunctionPlainBinary(&Psi, PsiFileName))
    {
      printf("Fatal Error in main\n Failed to read Psi from file %s !\n", PsiFileName);
      return 1;
    };
    printf("\nRead the starting wave function from the file %s\n", PsiFileName);
    sscanf(Psi.Name, "%s %i %i", Buffer, &PsiNBathC, &PsiNBathV);
    printf("Wave-function \"%s\" found has %i conduction bath sites, %i valence bath sites and %i Fermions in basis\n",
           Psi.Name, PsiNBathC, PsiNBathV, Psi.NFermion);
    if ((2 * (PsiNBathC + PsiNBathV + 1)) != Psi.NFermion)
    {
      printf("The total number of Fermions is incompatible witht the bath sizes\n");
      printf("This wave-function can not be used to start from\n");
      return 1;
    }
  }
  else
  {
    printf("\nWave-function file: '%s' does not excist\n Started from new function\n", PsiFileName);
    CreateWaveFunctionDMFT(&Psi, HSet.NBath, HSet.NBath, HSet.NBitsKey);
  }
  strcpy(GBath.Name, "Gbath");
  /* do the DMFT loops */
  for (iOuter = HSet.iStart; iOuter < HSet.iStart + HSet.Ni; iOuter++)
  {
    BeginTiming("DMFTLoop", Time);
    printf("\n\n***********************************************************************\n");
    printf("Start of loop %i from %i to %i\n", iOuter + 1, HSet.iStart, HSet.iStart + HSet.Ni);
    printf("***********************************************************************\n\n\n");
    /* Create splitted anderson matrix */
    TriDiagonalMatrixToDoubleTriDiagonalMatrix(GBath, &T, &GAnd, Ef, epsilonAnd);
    /* Create GAnd and GAndTrunk as a copy pointing to the same memory */
    GAndTrunk.Ef = GAnd.Ef;
    GAndTrunk.Epsilon = GAnd.Epsilon;
    GAndTrunk.a0 = GAnd.a0;
    GAndTrunk.b0 = GAnd.b0;
    GAndTrunk.val.N = GAnd.val.N;
    GAndTrunk.val.a = GAnd.val.a;
    GAndTrunk.val.b = GAnd.val.b;
    GAndTrunk.con.N = GAnd.con.N;
    GAndTrunk.con.a = GAnd.con.a;
    GAndTrunk.con.b = GAnd.con.b;
    /* Create VPerturbation */
    for (i = 0; i < 4; i++)
    {
      InitOperator(&(VPerturbation[i]));
      VPerturbation[i].NFermion = 2 + 2 * (Min(GAnd.val.N, HSet.NBath) + Min(GAnd.con.N, HSet.NBath));
      VPerturbation[i].NBoson = 0;
      VPerturbation[i].NWord = Min(VPerturbation[i].NFermion, 1 + (VPerturbation[i].NFermion - 1) / 8);
      VPerturbation[i].QComplex = 0;
      InitOperatorFixedLength(&(VPerturbation[i]), 1, 1);
    }
    if (GAnd.val.N > HSet.NBath)
    {
      oadd[0] = CrAnMask + 1 + HSet.NBath;
      RealOperatorAddLader(&(VPerturbation[0]), oadd, 1, GAnd.val.b[HSet.NBath]);
      oadd[0] = CrAnMask + 1 + Min(GAnd.con.N, HSet.NBath) + 2 * HSet.NBath;
      RealOperatorAddLader(&(VPerturbation[1]), oadd, 1, GAnd.val.b[HSet.NBath]);
    }
    if (GAnd.con.N > HSet.NBath)
    {
      oadd[0] = 1 + Min(GAnd.val.N, HSet.NBath) + HSet.NBath;
      RealOperatorAddLader(&(VPerturbation[2]), oadd, 1, GAnd.con.b[HSet.NBath]);
      oadd[0] = 1 + 2 * Min(GAnd.val.N, HSet.NBath) + 2 * HSet.NBath;
      RealOperatorAddLader(&(VPerturbation[3]), oadd, 1, GAnd.con.b[HSet.NBath]);
    }
    /* Create gPerturbation */
    if (GAnd.val.N > HSet.NBath)
    {
      gPerturbation[0].N = GAnd.val.N - HSet.NBath;
      gPerturbation[1].N = GAnd.val.N - HSet.NBath;
      InitTriDiagonalMatrix(&(gPerturbation[0]));
      InitTriDiagonalMatrix(&(gPerturbation[1]));
      for (i = 1; i < GAnd.val.N - HSet.NBath; i++)
      {
        gPerturbation[0].a[i] = -GAnd.val.a[i + HSet.NBath];
        gPerturbation[0].b[i] = GAnd.val.b[i + HSet.NBath];
        gPerturbation[1].a[i] = -GAnd.val.a[i + HSet.NBath];
        gPerturbation[1].b[i] = GAnd.val.b[i + HSet.NBath];
      }
      gPerturbation[0].a[GAnd.val.N - HSet.NBath] = -GAnd.val.a[GAnd.val.N];
      gPerturbation[1].a[GAnd.val.N - HSet.NBath] = -GAnd.val.a[GAnd.val.N];
    }
    else
    {
      gPerturbation[0].N = 1;
      gPerturbation[1].N = 1;
      InitTriDiagonalMatrix(&(gPerturbation[0]));
      InitTriDiagonalMatrix(&(gPerturbation[1]));
      gPerturbation[0].a[1] = 0;
      gPerturbation[1].a[1] = 0;
    }
    if (GAnd.con.N > HSet.NBath)
    {
      gPerturbation[2].N = GAnd.con.N - HSet.NBath;
      gPerturbation[3].N = GAnd.con.N - HSet.NBath;
      InitTriDiagonalMatrix(&(gPerturbation[2]));
      InitTriDiagonalMatrix(&(gPerturbation[3]));
      for (i = 1; i < GAnd.con.N - HSet.NBath; i++)
      {
        gPerturbation[2].a[i] = GAnd.con.a[i + HSet.NBath];
        gPerturbation[2].b[i] = GAnd.con.b[i + HSet.NBath];
        gPerturbation[3].a[i] = GAnd.con.a[i + HSet.NBath];
        gPerturbation[3].b[i] = GAnd.con.b[i + HSet.NBath];
      }
      gPerturbation[2].a[GAnd.con.N - HSet.NBath] = GAnd.con.a[GAnd.con.N];
      gPerturbation[3].a[GAnd.con.N - HSet.NBath] = GAnd.con.a[GAnd.con.N];
    }
    else
    {
      gPerturbation[2].N = 1;
      gPerturbation[3].N = 1;
      InitTriDiagonalMatrix(&(gPerturbation[2]));
      InitTriDiagonalMatrix(&(gPerturbation[3]));
      gPerturbation[2].a[1] = 0;
      gPerturbation[3].a[1] = 0;
    };
    /* Create HAnd */
    FreeMatrix(&T);
    if (GAnd.con.N < HSet.NBath)
    {
      printf("Warning: GAnd.con.N(%i)<HSet.NBath(%i)\n", GAnd.con.N, HSet.NBath);
    }
    else
    {
      GAndTrunk.con.N = HSet.NBath;
    }
    if (GAnd.val.N < HSet.NBath)
    {
      printf("Warning: GAnd.val.N(%i)<HSet.NBath(%i)\n", GAnd.val.N, HSet.NBath);
    }
    else
    {
      GAndTrunk.val.N = HSet.NBath;
    }
    mu = -HSet.U / 2;
    DoubleTriDiagonalMatrixToAndersonImpurityOperator(&HAnd, GAndTrunk, GAndTrunk, HSet.U, mu);
    /* test if psi has correct dimensions and correct if needed */
    sscanf(Psi.Name, "%s %i %i", Buffer, &PsiNBathC, &PsiNBathV);
    if (PsiNBathC != GAndTrunk.con.N)
    {
      printf("\nBefore diagonalizing:\nThe wavefunction has %i conduction bath sites\nThe Anderson Hamiltonian %i\n",
             PsiNBathC, GAndTrunk.con.N);
      if (PsiNBathC > GAndTrunk.con.N)
      {
        printf("Removing %i unoccupied fermions from the wave-function\n", PsiNBathC - GAndTrunk.con.N);
        for (i = 0; i < PsiNBathC - GAndTrunk.con.N; i++)
        {
          WaveFunctionRemoveFermion(Psi, &PsiTemp, Psi.NFermion - 1);
          WaveFunctionFree(&Psi);
          WaveFunctionRemoveFermion(PsiTemp, &Psi, 2 + PsiNBathV + PsiNBathC - i - 1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      else
      {
        printf("Adding %i unoccupied fermions to the wave-function\n", GAndTrunk.con.N - PsiNBathC);
        for (i = 0; i < GAndTrunk.con.N - PsiNBathC; i++)
        {
          WaveFunctionAddFermion(Psi, &PsiTemp, Psi.NFermion, 0);
          WaveFunctionFree(&Psi);
          WaveFunctionAddFermion(PsiTemp, &Psi, 2 + PsiNBathV + PsiNBathC + i, 0);
          WaveFunctionFree(&PsiTemp);
        }
      }
    }
    if (PsiNBathV != GAndTrunk.val.N)
    {
      printf("\nBefore diagonalizing:\nThe wavefunction has %i valence bath sites\nThe Anderson Hamiltonian %i\n",
             PsiNBathV, GAndTrunk.val.N);
      if (PsiNBathV > GAndTrunk.val.N)
      {
        printf("Removing %i occupied fermions from the wave-function\n", PsiNBathV - GAndTrunk.val.N);
        for (i = 0; i < PsiNBathV - GAndTrunk.val.N; i++)
        {
          WaveFunctionRemoveFermion(Psi, &PsiTemp, Psi.NFermion - GAndTrunk.con.N - 1);
          WaveFunctionFree(&Psi);
          WaveFunctionRemoveFermion(PsiTemp, &Psi, 2 + PsiNBathV - i - 1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      else
      {
        printf("Adding %i occupied fermions to the wave-function\n", GAndTrunk.val.N - PsiNBathV);
        for (i = 0; i < GAndTrunk.val.N - PsiNBathV; i++)
        {
          WaveFunctionAddFermion(Psi, &PsiTemp, Psi.NFermion - GAndTrunk.val.N, 1);
          WaveFunctionFree(&Psi);
          WaveFunctionAddFermion(PsiTemp, &Psi, 2 + PsiNBathV + i, 1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      printf("Corrected basis\n");
      PrintWaveFunctionReadableTruncated(Psi, 0.1);
    }
    sprintf(Psi.Name, "DMFT %i %i", GAndTrunk.con.N, GAndTrunk.val.N);
    /* create ground-state wave-function */
    BeginTiming("Lanczos GroundState", Time);
    LanczosGroundState(HAnd, &Psi, &HSet.zeroLancWavefunction, HSet.zeroLancEnergy, HSet.U, &Psi, 0);
    LanczosGroundState(HAnd, &Psi, &HSet.zeroLancWavefunction, HSet.zeroLancEnergy, HSet.U, &Psi, 0);
    EndTiming("Lanczos GroundState", Time);
    PsiHermitianOperatorPsi(HAnd, Psi, &PsiEnergy);
    PsiHermitianOperatorSqrPsi(HAnd, Psi, &PsiQuality);
    PsiQuality /= Sqr(PsiEnergy);
    PrintWaveFunctionReadableTruncated(Psi, 0.1);
    CalculatePrefactorHistogram(PsiHist, Psi, 10.0, 20);
    printf("\nHistogram of prefactors in Wave-function\n");
    for (i = 0; i < 19; i++)
    {
      printf("10^(-%2i) to 10^(-%2i) | %22.15E\n", i, i + 1, PsiHist[i]);
    }
    printf("10^(-%2i) to 0        | %22.15E\n", i, PsiHist[19]);
    printf("\n");
    printf("Energy                  <psi|H^2|psi>/(<psi|H|psi>^2)\n");
    printf("%22.15E  %22.15E\n", PsiEnergy, PsiQuality);
    PsiRho = calloc(Psi.NFermion, sizeof(double));
    if (PsiRho == NULL)
    {
      printf("calloc failed in main for PsiRho\n Fatal Error\n");
      fflush(stdout);
      return 1;
    }
    CalculateDensity(PsiRho, Psi);
    printf("Occupation \n");
    for (i = 0; i < Psi.NFermion; i++)
    {
      printf("%5i %22.15E\n", i, PsiRho[i]);
    }
    free(PsiRho);
    /*                                           */
    /* create spin up and spin down PES and IPES */
    /*                                           */
    CreatePESdnOperatorDMFT(&OPESdn, HAnd.NFermion);
    CreatePESupOperatorDMFT(&OPESup, HAnd.NFermion);
    CreateIPESdnOperatorDMFT(&OIPESdn, HAnd.NFermion);
    CreateIPESupOperatorDMFT(&OIPESup, HAnd.NFermion);
    GclMup.N = Max(HSet.NTriDiag, G0.N);
    InitTriDiagonalMatrix(&GclMup);
    GclMdn.N = Max(HSet.NTriDiag, G0.N);
    InitTriDiagonalMatrix(&GclMdn);
    GclPup.N = Max(HSet.NTriDiag, G0.N);
    InitTriDiagonalMatrix(&GclPup);
    GclPdn.N = Max(HSet.NTriDiag, G0.N);
    InitTriDiagonalMatrix(&GclPdn);
    BeginTiming("Lanczos Tri Diagonalization Dyson", Time);
    OperatorPsi(OPESup, Psi, &OPsi);
    LanczosTriDiagonalizeDyson(HAnd, OPsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation, gPerturbation, 4,
                               &GclMup);
    EndTiming("Lanczos Tri Diagonalization Dyson", Time);
    BeginTiming("Lanczos Tri Diagonalization Dyson", Time);
    OperatorPsi(OPESdn, Psi, &OPsi);
    LanczosTriDiagonalizeDyson(HAnd, OPsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation, gPerturbation, 4,
                               &GclMdn);
    EndTiming("Lanczos Tri Diagonalization Dyson", Time);
    BeginTiming("Lanczos Tri Diagonalization Dyson", Time);
    OperatorPsi(OIPESup, Psi, &OPsi);
    LanczosTriDiagonalizeDyson(HAnd, OPsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation, gPerturbation, 4,
                               &GclPup);
    EndTiming("Lanczos Tri Diagonalization Dyson", Time);
    BeginTiming("Lanczos Tri Diagonalization Dyson", Time);
    OperatorPsi(OIPESdn, Psi, &OPsi);
    LanczosTriDiagonalizeDyson(HAnd, OPsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation, gPerturbation, 4,
                               &GclPdn);
    EndTiming("Lanczos Tri Diagonalization Dyson", Time);
    /* G math */
    printf("Start G math\n");
    ShiftTriDiagonalMatrix(&GclMdn, -PsiEnergy);
    ShiftTriDiagonalMatrix(&GclMup, -PsiEnergy);
    ShiftTriDiagonalMatrix(&GclPdn, -PsiEnergy);
    ShiftTriDiagonalMatrix(&GclPup, -PsiEnergy);
    for (i = 1; i <= GclMdn.N; i++)
    {
      GclMdn.a[i] *= -1;
    }
    for (i = 1; i <= GclMup.N; i++)
    {
      GclMup.a[i] *= -1;
    }
    strcpy(GclPup.Name, "Gcluster-Plus-Up");
    strcpy(GclPdn.Name, "Gcluster-Plus-Dn");
    strcpy(GclMup.Name, "Gcluster-Min-Up");
    strcpy(GclMdn.Name, "Gcluster-Min-Dn");
    printf("  Start SumGreensfunctionsTriDiagonalRepresentation\n");
    fflush(stdout);
    SumGreensfunctionsTriDiagonalRepresentation(GclPup, GclMup, &Gclup);
    strcpy(Gclup.Name, "Gcluster-Up");
    SumGreensfunctionsTriDiagonalRepresentation(GclPdn, GclMdn, &Gcldn);
    strcpy(Gcldn.Name, "Gcluster-Dn");
    SumGreensfunctionsTriDiagonalRepresentation(Gclup, Gcldn, &Gcl);
    strcpy(Gcl.Name, "Gcluster");
    Gcl.b[0] = 1.0;
    printf("Gcl.N = %i ", Gcl.N);
    SymmetrizeTridiagonalMatrix(&Gcl);
    Gcl.N = Min(Gcl.N, Max(4 * HSet.NTriDiag, 2 * G0.N));
    printf("=> %i\n", Gcl.N);
    printf("  Start SelfEnergyTriDiagonalRepresentation\n");
    printf("GBath.N = %i ", GBath.N);
    SelfEnergyTriDiagonalRepresentation(GBath, Gcl, &Sigma);
    printf("Sigma.N = %i ", Sigma.N);
    Sigma.N = Min(Sigma.N, Max(4 * HSet.NTriDiag, 2 * G0.N));
    printf("=> %i ", Sigma.N);
    strcpy(Sigma.Name, "Sigma");
    SymmetrizeTridiagonalMatrix(&Sigma);
    printf("=> %i ", Sigma.N);
    Sigma.N = Min(Sigma.N, Max(4 * HSet.NTriDiag, 2 * G0.N));
    printf("=> %i\n", Sigma.N);
    printf("  Start GLocalTriDiagonalRepresentation\n");
    GLocalTriDiagonalRepresentation(G0, Sigma, &GLocal, Sigma.N + G0.N - 1);
    printf("Glocal.N = %i ", GLocal.N);
    GLocal.N = Min(GLocal.N, Max(4 * HSet.NTriDiag, 2 * G0.N));
    printf("=> %i ", GLocal.N);
    strcpy(GLocal.Name, "Glocal");
    SymmetrizeTridiagonalMatrix(&GLocal);
    printf("=> %i ", GLocal.N);
    GLocal.N = Min(GLocal.N, Max(4 * HSet.NTriDiag, 2 * G0.N));
    printf("=> %i\n", GLocal.N);
    printf("  Start GBathNewTriDiagonalRepresentation\n");
    GBathNewTriDiagonalRepresentation(GLocal, Sigma, &GBathNew);
    GBathNew.N = Min(GBathNew.N, Max(4 * HSet.NTriDiag, 2 * G0.N));
    strcpy(GBathNew.Name, "Gbath-New");
    SymmetrizeTridiagonalMatrix(&GBathNew);
    GBathNew.N = Min(GBathNew.N, Max(4 * HSet.NTriDiag, 2 * G0.N));
    /* calculate GBath for the next loop by mixing GBath with GBathNew */
    printf("Start mixing GBath with GBathNew\n");
    TriDiagonalMatrixCopy(&GBathOld, GBath);
    strcpy(GBathOld.Name, "Gbath-Old");
    FreeTriDiagonalMatrix(&GBath);
    GBathOld.b[0] = sqrt(HSet.Mix);
    GBathNew.b[0] = sqrt(1 - HSet.Mix);
    SumGreensfunctionsTriDiagonalRepresentation(GBathNew, GBathOld, &GBath);
    strcpy(GBath.Name, "Gbath");
    GBathOld.b[0] = 1;
    GBathNew.b[0] = 1;
    SymmetrizeTridiagonalMatrix(&GBath);
    GBath.N = Min(GBath.N, Max(4 * HSet.NTriDiag, 2 * G0.N));
    /* create spectra */
    printf("Start create spectra\n");
    Spectra[0].N = HSet.NE;
    Spectra[0].Emin = HSet.Emin;
    Spectra[0].Emax = HSet.Emax;
    Spectra[0].dE = (HSet.Emax - HSet.Emin) / HSet.NE;
    Spectra[0].EminPole = HSet.Emin;
    Spectra[0].EmaxPole = HSet.Emax;
    Spectra[0].Gamma = HSet.Gamma;
    for (i = 0; i < 10; i++)
    {
      SpectrumAllocate(Spectra[0], &Spectra[i]);
    }
    TriDiagonalMatrixToSpectrum(0, GBath, &(Spectra[0])); /*  3 */
    TriDiagonalMatrixToSpectrum(0, GBathOld, &(Spectra[1]));  /*  5 */
    TriDiagonalMatrixToSpectrum(0, GBathNew, &(Spectra[2]));  /*  7 */
    DoubleTriDiagonalMatrixToSpectrum(0, GAnd, &(Spectra[3]));  /*  9 */
    TriDiagonalMatrixToSpectrum(0, G0, &(Spectra[4]));  /* 11 */
    TriDiagonalMatrixToSpectrum(0, Sigma, &(Spectra[5])); /* 13 */
    TriDiagonalMatrixToSpectrum(0, GLocal, &(Spectra[6]));  /* 15 */
    TriDiagonalMatrixToSpectrum(0, Gclup, &(Spectra[7])); /* 17 */
    TriDiagonalMatrixToSpectrum(0, Gcldn, &(Spectra[8])); /* 19 */
    TriDiagonalMatrixToSpectrum(0, Gcl, &(Spectra[9])); /* 21 */
    /* write to disk */
    printf("Start write to dist\n");
    if (iOuter == 0)
    {
      sprintf(GBathFileName, "%s001", HSet.GRootFileName);
      sprintf(PsiFileName, "%s001", HSet.WaveFunctionRootFileName);
      sprintf(SpecFileName, "%s001", HSet.SpecRootFileName);
    }
    else if (iOuter < 9)
    {
      sprintf(GBathFileName, "%s00%1i", HSet.GRootFileName, iOuter + 1);
      sprintf(PsiFileName, "%s00%1i", HSet.WaveFunctionRootFileName, iOuter + 1);
      sprintf(SpecFileName, "%s00%1i", HSet.SpecRootFileName, iOuter + 1);
    }
    else if (iOuter < 99)
    {
      sprintf(GBathFileName, "%s0%2i", HSet.GRootFileName, iOuter + 1);
      sprintf(PsiFileName, "%s0%2i", HSet.WaveFunctionRootFileName, iOuter + 1);
      sprintf(SpecFileName, "%s0%2i", HSet.SpecRootFileName, iOuter + 1);
    }
    else if (iOuter < 999)
    {
      sprintf(GBathFileName, "%s%3i", HSet.GRootFileName, iOuter + 1);
      sprintf(PsiFileName, "%s%3i", HSet.WaveFunctionRootFileName, iOuter + 1);
      sprintf(SpecFileName, "%s%3i", HSet.SpecRootFileName, iOuter + 1);
    }
    FilePrintSpectra(Spectra, 10, SpecFileName, "w");
    FilePrintPlainTriDiagonalMatrix(GBath, GBathFileName, "w");
    FilePrintWaveFunctionPlainBinary(Psi, PsiFileName, "w");
    /* free variables */
    printf("Start free variables\n");
    for (i = 0; i < 4; i++)
    {
      FreeOperator(&(VPerturbation[i]));
      FreeTriDiagonalMatrix(&(gPerturbation[i]));
    }
    FreeOperator(&HAnd);
    FreeOperator(&OPESdn);
    FreeOperator(&OPESup);
    FreeOperator(&OIPESdn);
    FreeOperator(&OIPESup);
    FreeTriDiagonalMatrix(&GclMup);
    FreeTriDiagonalMatrix(&GclMdn);
    FreeTriDiagonalMatrix(&GclPup);
    FreeTriDiagonalMatrix(&GclPdn);
    FreeTriDiagonalMatrix(&Gcl);
    FreeTriDiagonalMatrix(&GLocal);
    FreeTriDiagonalMatrix(&Sigma);
    FreeTriDiagonalMatrix(&GBathNew);
    FreeTriDiagonalMatrix(&GBathOld);
    FreeDoubleTriDiagonalMatrix(&GAnd);
    EndTiming("DMFTLoop", Time);

    PrintTiming(*Time);
  }
  return 0;
};

int oldmainDysonRho(int argc, const char *argv[])
{
  DMFTSettings HSet;
  WaveFunctionType Psi, OPsi, PsiTemp;
  TriDiagonalMatrixType G0, GBath, GclPup, GclPdn, GclMup, GclMdn, Gclup, Gcldn, Gcl, Sigma, GLocal, GBathNew, GBathOld,
    Gtemp, Gperturb[4];
  DoubleTriDiagonalMatrixType GAnd, GAndOrig;
  OperatorType HAnd, OPESdn, OPESup, OIPESdn, OIPESup;
  TimeType Time;
  MatrixType T, Rho, RhoBound;
  SpectraType Spectra[10];
  char   GBathFileName[255], PsiFileName[255], SpecFileName[25], Buffer[255];
  uint32_t iOuter, i, PsiNBathV, PsiNBathC;
  double Ef = 0.0, epsilonAnd = 0.0000001, mu;
  double PsiEnergy, PsiQuality, PsiHist[20], *PsiRho;
  InitializeTiming(&Time);
  BeginTiming("Total", &Time);
  TimeStampAll();
  argc -= 1;
  /* read the DMFT settings */
  if (FileReadDMFTSettings(&HSet, argv[argc]))
  {
    printf("Fatal Error in main\n Failed to read DMFT Settings from file %s\n", argv[argc]);
    return 1;
  };
  PrintDMFTSettings(HSet);
  /* read G0 */
  if (FileReadPlainTriDiagonalMatrixType(&G0, HSet.G0FileName))
  {
    printf("Fatal Error in main\n Failed to read G0 from file %s\n", HSet.G0FileName);
    return 1;
  };
  strcpy(G0.Name, "G0");
  /* read G_Bath_iStart and psi_iStart */
  if (HSet.iStart == 0)
  {
    sprintf(GBathFileName, "%s000", HSet.GRootFileName);
    sprintf(PsiFileName, "%s000", HSet.WaveFunctionRootFileName);
  }
  else if (HSet.iStart < 10)
  {
    sprintf(GBathFileName, "%s00%1i", HSet.GRootFileName, HSet.iStart);
    sprintf(PsiFileName, "%s00%1i", HSet.WaveFunctionRootFileName, HSet.iStart);
  }
  else if (HSet.iStart < 100)
  {
    sprintf(GBathFileName, "%s0%2i", HSet.GRootFileName, HSet.iStart);
    sprintf(PsiFileName, "%s0%2i", HSet.WaveFunctionRootFileName, HSet.iStart);
  }
  else
  {
    sprintf(GBathFileName, "%s%i", HSet.GRootFileName, HSet.iStart);
    sprintf(PsiFileName, "%s%i", HSet.WaveFunctionRootFileName, HSet.iStart);
  }
  if (FileReadPlainTriDiagonalMatrixType(&GBath, GBathFileName))
  {
    printf("Fatal Error in main\n Failed to read GBath from file %s !\n", GBathFileName);
    return 1;
  };
  if (file_exists(PsiFileName))
  {
    Psi.NBitsKey = HSet.NBitsKey;
    if (FileReadWaveFunctionPlainBinary(&Psi, PsiFileName))
    {
      printf("Fatal Error in main\n Failed to read Psi from file %s !\n", PsiFileName);
      return 1;
    };
    printf("\nRead the starting wave function from the file %s\n", PsiFileName);
    sscanf(Psi.Name, "%s %i %i", Buffer, &PsiNBathC, &PsiNBathV);
    printf("Wave-function \"%s\" found has %i conduction bath sites, %i valence bath sites and %i Fermions in basis\n",
           Psi.Name, PsiNBathC, PsiNBathV, Psi.NFermion);
    if ((2 * (PsiNBathC + PsiNBathV + 1)) != Psi.NFermion)
    {
      printf("The total number of Fermions is incompatible witht the bath sizes\n");
      printf("This wave-function can not be used to start from\n");
      return 1;
    }
  }
  else
  {
    printf("\nWave-function file: '%s' does not excist\n Started from new function\n", PsiFileName);
    CreateWaveFunctionDMFT(&Psi, HSet.NBath, HSet.NBath, HSet.NBitsKey);
  }
  strcpy(GBath.Name, "Gbath");
  /* allocate Rho matrices */
  Rho.QComplex = 0;
  Rho.Nj = 2 * (HSet.NBath + 1);
  Rho.Ni = HSet.NTriDiag;
  if (InitMatrix(&Rho))
  {
    printf("InitMatrix(&Rho) failed in main\n");
    fflush(stdout);
    return 1;
  };
  RhoBound.QComplex = 0;
  RhoBound.Nj = 4;
  RhoBound.Ni = HSet.NTriDiag;
  if (InitMatrix(&RhoBound))
  {
    printf("InitMatrix(&RhoBound) failed in main\n");
    fflush(stdout);
    return 1;
  };
  /* do the DMFT loops */
  for (iOuter = HSet.iStart; iOuter < HSet.iStart + HSet.Ni; iOuter++)
  {
    BeginTiming("DMFTLoop", &Time);
    printf("\n\n***********************************************************************\n");
    printf("Start of loop %i from %i to %i\n", iOuter + 1, HSet.iStart, HSet.iStart + HSet.Ni);
    printf("***********************************************************************\n\n\n");
    /* create GBath as a double tri-diagonal matrix */
    TriDiagonalMatrixToDoubleTriDiagonalMatrix(GBath, &T, &GAnd, Ef, epsilonAnd);
    FreeMatrix(&T);
    /* make GAndOrig point to GAnd with original size */
    GAndOrig = GAnd;
    /* reduce size of GBath for exact diagonalization */
    if (GAnd.con.N < HSet.NBath)
    {
      printf("Warning: GAnd.con.N(%i)<HSet.NBath(%i)\n", GAnd.con.N, HSet.NBath);
    }
    else
    {
      GAnd.con.N = HSet.NBath;
    }
    if (GAnd.val.N < HSet.NBath)
    {
      printf("Warning: GAnd.val.N(%i)<HSet.NBath(%i)\n", GAnd.val.N, HSet.NBath);
    }
    else
    {
      GAnd.val.N = HSet.NBath;
    }
    mu = -HSet.U / 2;
    DoubleTriDiagonalMatrixToAndersonImpurityOperator(&HAnd, GAnd, GAnd, HSet.U, mu);
    /* test if psi has correct dimensions and correct if needed */
    sscanf(Psi.Name, "%s %i %i", Buffer, &PsiNBathC, &PsiNBathV);
    if (PsiNBathC != GAnd.con.N)
    {
      printf("\nBefore diagonalizing:\nThe wavefunction has %i conduction bath sites\nThe Anderson Hamiltonian %i\n",
             PsiNBathC, GAnd.con.N);
      if (PsiNBathC > GAnd.con.N)
      {
        printf("Removing %i unoccupied fermions from the wave-function\n", PsiNBathC - GAnd.con.N);
        for (i = 0; i < PsiNBathC - GAnd.con.N; i++)
        {
          WaveFunctionRemoveFermion(Psi, &PsiTemp, Psi.NFermion - 1);
          WaveFunctionFree(&Psi);
          WaveFunctionRemoveFermion(PsiTemp, &Psi, 2 + PsiNBathV + PsiNBathC - i - 1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      else
      {
        printf("Adding %i unoccupied fermions to the wave-function\n", GAnd.con.N - PsiNBathC);
        for (i = 0; i < GAnd.con.N - PsiNBathC; i++)
        {
          WaveFunctionAddFermion(Psi, &PsiTemp, Psi.NFermion, 0);
          WaveFunctionFree(&Psi);
          WaveFunctionAddFermion(PsiTemp, &Psi, 2 + PsiNBathV + PsiNBathC + i, 0);
          WaveFunctionFree(&PsiTemp);
        }
      }
    }
    if (PsiNBathV != GAnd.val.N)
    {
      printf("\nBefore diagonalizing:\nThe wavefunction has %i valence bath sites\nThe Anderson Hamiltonian %i\n",
             PsiNBathV, GAnd.val.N);
      if (PsiNBathV > GAnd.val.N)
      {
        printf("Removing %i occupied fermions from the wave-function\n", PsiNBathV - GAnd.val.N);
        for (i = 0; i < PsiNBathV - GAnd.val.N; i++)
        {
          WaveFunctionRemoveFermion(Psi, &PsiTemp, Psi.NFermion - GAnd.con.N - 1);
          WaveFunctionFree(&Psi);
          WaveFunctionRemoveFermion(PsiTemp, &Psi, 2 + PsiNBathV - i - 1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      else
      {
        printf("Adding %i occupied fermions to the wave-function\n", GAnd.val.N - PsiNBathV);
        for (i = 0; i < GAnd.val.N - PsiNBathV; i++)
        {
          WaveFunctionAddFermion(Psi, &PsiTemp, Psi.NFermion - GAnd.val.N, 1);
          WaveFunctionFree(&Psi);
          WaveFunctionAddFermion(PsiTemp, &Psi, 2 + PsiNBathV + i, 1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      printf("Corrected basis\n");
      PrintWaveFunctionReadableTruncated(Psi, 0.1);
    }
    sprintf(Psi.Name, "DMFT %i %i", GAnd.con.N, GAnd.val.N);
    /* create ground-state wave-function */
    BeginTiming("Lanczos GroundState", &Time);
    LanczosGroundState(HAnd, &Psi, &HSet.zeroLancWavefunction, HSet.zeroLancEnergy, HSet.U, &Psi, 0);
    LanczosGroundState(HAnd, &Psi, &HSet.zeroLancWavefunction, HSet.zeroLancEnergy, HSet.U, &Psi, 0);
    EndTiming("Lanczos GroundState", &Time);
    PsiHermitianOperatorPsi(HAnd, Psi, &PsiEnergy);
    PsiHermitianOperatorSqrPsi(HAnd, Psi, &PsiQuality);
    PsiQuality /= Sqr(PsiEnergy);
    PrintWaveFunctionReadableTruncated(Psi, 0.1);
    CalculatePrefactorHistogram(PsiHist, Psi, 10.0, 20);
    printf("\nHistogram of prefactors in Wave-function\n");
    for (i = 0; i < 19; i++)
    {
      printf("10^(-%2i) to 10^(-%2i) | %22.15E\n", i, i + 1, PsiHist[i]);
    }
    printf("10^(-%2i) to 0        | %22.15E\n", i, PsiHist[19]);
    printf("\n");
    printf("Energy                  <psi|H^2|psi>/(<psi|H|psi>^2)\n");
    printf("%22.15E  %22.15E\n", PsiEnergy, PsiQuality);
    PsiRho = calloc(Psi.NFermion, sizeof(double));
    if (PsiRho == NULL)
    {
      printf("calloc failed in main for PsiRho\n Fatal Error\n");
      fflush(stdout);
      return 1;
    }
    CalculateDensity(PsiRho, Psi);
    printf("Occupation \n");
    for (i = 0; i < Psi.NFermion; i++)
    {
      printf("%5i %22.15E\n", i, PsiRho[i]);
    }
    free(PsiRho);
    /* create Gperturb */
    if (GAndOrig.val.N <= HSet.NBath)
    {
      Gperturb[0].N = 0;
      Gperturb[2].N = 0;
    }
    else
    {
      Gperturb[0].N = GAndOrig.val.N - HSet.NBath;
      Gperturb[0].a = &GAndOrig.val.a[HSet.NBath];
      Gperturb[0].b = &GAndOrig.val.b[HSet.NBath];
      Gperturb[2].N = GAndOrig.val.N - HSet.NBath;
      Gperturb[2].a = &GAndOrig.val.a[HSet.NBath];
      Gperturb[2].b = &GAndOrig.val.b[HSet.NBath];
    }
    if (GAndOrig.con.N <= HSet.NBath)
    {
      Gperturb[1].N = 0;
      Gperturb[3].N = 0;
    }
    else
    {
      Gperturb[1].N = GAndOrig.con.N - HSet.NBath;
      Gperturb[1].a = &GAndOrig.con.a[HSet.NBath];
      Gperturb[1].b = &GAndOrig.con.b[HSet.NBath];
      Gperturb[3].N = GAndOrig.con.N - HSet.NBath;
      Gperturb[3].a = &GAndOrig.con.a[HSet.NBath];
      Gperturb[3].b = &GAndOrig.con.b[HSet.NBath];
    }



    /* create spin up and spin down PES and IPES */
    CreatePESdnOperatorDMFT(&OPESdn, HAnd.NFermion);
    CreatePESupOperatorDMFT(&OPESup, HAnd.NFermion);
    CreateIPESdnOperatorDMFT(&OIPESdn, HAnd.NFermion);
    CreateIPESupOperatorDMFT(&OIPESup, HAnd.NFermion);
    GclMup.N = HSet.NTriDiag;
    InitTriDiagonalMatrix(&GclMup);
    GclMdn.N = HSet.NTriDiag;
    InitTriDiagonalMatrix(&GclMdn);
    GclPup.N = HSet.NTriDiag;
    InitTriDiagonalMatrix(&GclPup);
    GclPdn.N = HSet.NTriDiag;
    InitTriDiagonalMatrix(&GclPdn);
    Gtemp.N = HSet.NTriDiag;
    InitTriDiagonalMatrix(&Gtemp);
    BeginTiming("Lanczos Tri Diagonalization PES up", &Time);
    OperatorPsi(OPESup, Psi, &OPsi);
    LanczosTriDiagonalizePlusDensity(HAnd, OPsi, HSet.zeroLancWavefunction, &Gtemp, &Rho, HSet.Gamma);
    /* copy the relavant elements from rho to rhobound */
    for (i = 0; i < Rho.Ni; i++)
    {
      RhoBound.VRe[i][0] = 1 - Rho.VRe[i][GAnd.val.N + 1];
      RhoBound.VRe[i][1] = Rho.VRe[i][GAnd.val.N + GAnd.con.N + 1];
      RhoBound.VRe[i][2] = 1 - Rho.VRe[i][GAnd.val.N + GAnd.con.N + GAnd.val.N + 1];
      RhoBound.VRe[i][3] = Rho.VRe[i][GAnd.val.N + GAnd.con.N + GAnd.val.N + GAnd.con.N + 1];
    }
    PrintMatrix(RhoBound);
    /* use Dyson to add more bath sites to spectrum than ground-state */
    TriDiagonalMatrixDysonScalar(Gtemp, &GclMup, RhoBound, 4, Gperturb);
    EndTiming("Lanczos Tri Diagonalization PES up", &Time);
    BeginTiming("Lanczos Tri Diagonalization PES dn", &Time);
    OperatorPsi(OPESdn, Psi, &OPsi);
    LanczosTriDiagonalizePlusDensity(HAnd, OPsi, HSet.zeroLancWavefunction, &Gtemp, &Rho, HSet.Gamma);
    /* copy the relavant elements from rho to rhobound */
    for (i = 0; i < Rho.Ni; i++)
    {
      RhoBound.VRe[i][0] = 1 - Rho.VRe[i][GAnd.val.N + 1];
      RhoBound.VRe[i][1] = Rho.VRe[i][GAnd.val.N + GAnd.con.N + 1];
      RhoBound.VRe[i][2] = 1 - Rho.VRe[i][GAnd.val.N + GAnd.con.N + GAnd.val.N + 1];
      RhoBound.VRe[i][3] = Rho.VRe[i][GAnd.val.N + GAnd.con.N + GAnd.val.N + GAnd.con.N + 1];
    }
    PrintMatrix(RhoBound);
    /* use Dyson to add more bath sites to spectrum than ground-state */
    TriDiagonalMatrixDysonScalar(Gtemp, &GclMdn, RhoBound, 4, Gperturb);
    PrintTriDiagonalMatrix(Gtemp);
    PrintTriDiagonalMatrix(GclMdn);
    EndTiming("Lanczos Tri Diagonalization PES dn", &Time);
    BeginTiming("Lanczos Tri Diagonalization IPES up", &Time);
    OperatorPsi(OIPESup, Psi, &OPsi);
    LanczosTriDiagonalizePlusDensity(HAnd, OPsi, HSet.zeroLancWavefunction, &Gtemp, &Rho, HSet.Gamma);
    /* copy the relavant elements from rho to rhobound */
    for (i = 0; i < Rho.Ni; i++)
    {
      RhoBound.VRe[i][0] = 1 - Rho.VRe[i][GAnd.val.N + 1];
      RhoBound.VRe[i][1] = Rho.VRe[i][GAnd.val.N + GAnd.con.N + 1];
      RhoBound.VRe[i][2] = 1 - Rho.VRe[i][GAnd.val.N + GAnd.con.N + GAnd.val.N + 1];
      RhoBound.VRe[i][3] = Rho.VRe[i][GAnd.val.N + GAnd.con.N + GAnd.val.N + GAnd.con.N + 1];
    }
    PrintMatrix(RhoBound);
    /* use Dyson to add more bath sites to spectrum than ground-state */
    TriDiagonalMatrixDysonScalar(Gtemp, &GclPup, RhoBound, 4, Gperturb);
    EndTiming("Lanczos Tri Diagonalization IPES up", &Time);
    BeginTiming("Lanczos Tri Diagonalization IPES dn", &Time);
    OperatorPsi(OIPESdn, Psi, &OPsi);
    LanczosTriDiagonalizePlusDensity(HAnd, OPsi, HSet.zeroLancWavefunction, &Gtemp, &Rho, HSet.Gamma);
    /* copy the relavant elements from rho to rhobound */
    for (i = 0; i < Rho.Ni; i++)
    {
      RhoBound.VRe[i][0] = 1 - Rho.VRe[i][GAnd.val.N + 1];
      RhoBound.VRe[i][1] = Rho.VRe[i][GAnd.val.N + GAnd.con.N + 1];
      RhoBound.VRe[i][2] = 1 - Rho.VRe[i][GAnd.val.N + GAnd.con.N + GAnd.val.N + 1];
      RhoBound.VRe[i][3] = Rho.VRe[i][GAnd.val.N + GAnd.con.N + GAnd.val.N + GAnd.con.N + 1];
    }
    PrintMatrix(RhoBound);
    /* use Dyson to add more bath sites to spectrum than ground-state */
    TriDiagonalMatrixDysonScalar(Gtemp, &GclPdn, RhoBound, 4, Gperturb);
    PrintTriDiagonalMatrix(Gtemp);
    PrintTriDiagonalMatrix(GclPdn);
    EndTiming("Lanczos Tri Diagonalization IPES dn", &Time);
    /* G math */
    printf("Start G math\n");
    ShiftTriDiagonalMatrix(&GclMdn, -PsiEnergy);
    ShiftTriDiagonalMatrix(&GclMup, -PsiEnergy);
    ShiftTriDiagonalMatrix(&GclPdn, -PsiEnergy);
    ShiftTriDiagonalMatrix(&GclPup, -PsiEnergy);
    for (i = 1; i <= GclMdn.N; i++)
    {
      GclMdn.a[i] *= -1;
    }
    for (i = 1; i <= GclMup.N; i++)
    {
      GclMup.a[i] *= -1;
    }
    strcpy(GclPup.Name, "Gcluster-Plus-Up");
    strcpy(GclPdn.Name, "Gcluster-Plus-Dn");
    strcpy(GclMup.Name, "Gcluster-Min-Up");
    strcpy(GclMdn.Name, "Gcluster-Min-Dn");
    printf("  Start SumGreensfunctionsTriDiagonalRepresentation\n");
    fflush(stdout);
    SumGreensfunctionsTriDiagonalRepresentation(GclPup, GclMup, &Gclup);
    strcpy(Gclup.Name, "Gcluster-Up");
    SumGreensfunctionsTriDiagonalRepresentation(GclPdn, GclMdn, &Gcldn);
    strcpy(Gcldn.Name, "Gcluster-Dn");
    SumGreensfunctionsTriDiagonalRepresentation(Gclup, Gcldn, &Gcl);
    strcpy(Gcl.Name, "Gcluster");
    Gcl.b[0] = 1.0;
    SymmetrizeTridiagonalMatrix(&Gcl);
    Gcl.N = Min(Gcl.N, 4 * HSet.NTriDiag);
    printf("  Start SelfEnergyTriDiagonalRepresentation\n");
    SelfEnergyTriDiagonalRepresentation(GBath, Gcl, &Sigma);
    Sigma.N = Min(Sigma.N, 4 * HSet.NTriDiag);
    strcpy(Sigma.Name, "Sigma");
    SymmetrizeTridiagonalMatrix(&Sigma);
    Sigma.N = Min(Sigma.N, 4 * HSet.NTriDiag);
    printf("  Start GLocalTriDiagonalRepresentation\n");
    GLocalTriDiagonalRepresentation(G0, Sigma, &GLocal, HSet.NTriDiag);
    GLocal.N = Min(GLocal.N, 4 * HSet.NTriDiag);
    strcpy(GLocal.Name, "Glocal");
    SymmetrizeTridiagonalMatrix(&GLocal);
    GLocal.N = Min(GLocal.N, 4 * HSet.NTriDiag);
    printf("  Start GBathNewTriDiagonalRepresentation\n");
    GBathNewTriDiagonalRepresentation(GLocal, Sigma, &GBathNew);
    GBathNew.N = Min(GBathNew.N, 4 * HSet.NTriDiag);
    strcpy(GBathNew.Name, "Gbath-New");
    SymmetrizeTridiagonalMatrix(&GBathNew);
    GBathNew.N = Min(GBathNew.N, 4 * HSet.NTriDiag);
    /* calculate GBath for the next loop by mixing GBath with GBathNew */
    printf("Start mixing GBath with GBathNew\n");
    TriDiagonalMatrixCopy(&GBathOld, GBath);
    strcpy(GBathOld.Name, "Gbath-Old");
    FreeTriDiagonalMatrix(&GBath);
    GBathOld.b[0] = sqrt(HSet.Mix);
    GBathNew.b[0] = sqrt(1 - HSet.Mix);
    SumGreensfunctionsTriDiagonalRepresentation(GBathNew, GBathOld, &GBath);
    strcpy(GBath.Name, "Gbath");
    GBathOld.b[0] = 1;
    GBathNew.b[0] = 1;
    SymmetrizeTridiagonalMatrix(&GBath);
    GBath.N = Min(GBath.N, 4 * HSet.NTriDiag);
    /* create spectra */
    printf("Start create spectra\n");
    Spectra[0].N = HSet.NE;
    Spectra[0].Emin = HSet.Emin;
    Spectra[0].Emax = HSet.Emax;
    Spectra[0].dE = (HSet.Emax - HSet.Emin) / HSet.NE;
    Spectra[0].EminPole = HSet.Emin;
    Spectra[0].EmaxPole = HSet.Emax;
    Spectra[0].Gamma = HSet.Gamma;
    for (i = 0; i < 10; i++)
    {
      SpectrumAllocate(Spectra[0], &Spectra[i]);
    }
    TriDiagonalMatrixToSpectrum(0, GBath, &(Spectra[0]));
    TriDiagonalMatrixToSpectrum(0, GBathOld, &(Spectra[1]));
    TriDiagonalMatrixToSpectrum(0, GBathNew, &(Spectra[2]));
    DoubleTriDiagonalMatrixToSpectrum(0, GAnd, &(Spectra[3]));
    TriDiagonalMatrixToSpectrum(0, G0, &(Spectra[4]));
    TriDiagonalMatrixToSpectrum(0, Sigma, &(Spectra[5]));
    TriDiagonalMatrixToSpectrum(0, GLocal, &(Spectra[6]));
    TriDiagonalMatrixToSpectrum(0, Gclup, &(Spectra[7]));
    TriDiagonalMatrixToSpectrum(0, Gcldn, &(Spectra[8]));
    TriDiagonalMatrixToSpectrum(0, Gcl, &(Spectra[9]));
    /* write to disk */
    printf("Start write to dist\n");
    if (iOuter == 0)
    {
      sprintf(GBathFileName, "%s001", HSet.GRootFileName);
      sprintf(PsiFileName, "%s001", HSet.WaveFunctionRootFileName);
      sprintf(SpecFileName, "%s001", HSet.SpecRootFileName);
    }
    else if (iOuter < 9)
    {
      sprintf(GBathFileName, "%s00%1i", HSet.GRootFileName, iOuter + 1);
      sprintf(PsiFileName, "%s00%1i", HSet.WaveFunctionRootFileName, iOuter + 1);
      sprintf(SpecFileName, "%s00%1i", HSet.SpecRootFileName, iOuter + 1);
    }
    else if (iOuter < 99)
    {
      sprintf(GBathFileName, "%s0%2i", HSet.GRootFileName, iOuter + 1);
      sprintf(PsiFileName, "%s0%2i", HSet.WaveFunctionRootFileName, iOuter + 1);
      sprintf(SpecFileName, "%s0%2i", HSet.SpecRootFileName, iOuter + 1);
    }
    else if (iOuter < 999)
    {
      sprintf(GBathFileName, "%s%3i", HSet.GRootFileName, iOuter + 1);
      sprintf(PsiFileName, "%s%3i", HSet.WaveFunctionRootFileName, iOuter + 1);
      sprintf(SpecFileName, "%s%3i", HSet.SpecRootFileName, iOuter + 1);
    }
    FilePrintSpectra(Spectra, 10, SpecFileName, "w");
    FilePrintPlainTriDiagonalMatrix(GBath, GBathFileName, "w");
    FilePrintWaveFunctionPlainBinary(Psi, PsiFileName, "w");
    /* free variables */
    printf("Start free variables\n");
    FreeOperator(&HAnd);
    FreeOperator(&OPESdn);
    FreeOperator(&OPESup);
    FreeOperator(&OIPESdn);
    FreeOperator(&OIPESup);
    FreeTriDiagonalMatrix(&GclMup);
    FreeTriDiagonalMatrix(&GclMdn);
    FreeTriDiagonalMatrix(&GclPup);
    FreeTriDiagonalMatrix(&GclPdn);
    FreeTriDiagonalMatrix(&Gtemp);
    FreeTriDiagonalMatrix(&Gcl);
    FreeTriDiagonalMatrix(&GLocal);
    FreeTriDiagonalMatrix(&Sigma);
    FreeTriDiagonalMatrix(&GBathNew);
    FreeTriDiagonalMatrix(&GBathOld);
    FreeDoubleTriDiagonalMatrix(&GAnd);
    EndTiming("DMFTLoop", &Time);
  }
  EndTiming("Total", &Time);
  PrintTiming(Time);
  return 0;
};

int old2mainNoDyson(DMFTSettings HSet, TimeType * Time)
{
  WaveFunctionType Psi, OPsi, PsiTemp;
  TriDiagonalMatrixType G0, GBath, GclPup, GclPdn, GclMup, GclMdn, Gclup, Gcldn, Gcl, Sigma, GLocal, GBathNew, GBathOld;
  DoubleTriDiagonalMatrixType GAnd;
  OperatorType HAnd, OPESdn, OPESup, OIPESdn, OIPESup;
  MatrixType T;
  SpectraType Spectra[10];
  char   GBathFileName[255], PsiFileName[255], SpecFileName[25], Buffer[255];
  uint32_t iOuter, i, PsiNBathV, PsiNBathC;
  double Ef = 0.0, epsilonAnd = 0.0000001, mu;
  double PsiEnergy, PsiQuality, PsiHist[20], *PsiRho;
  /* read G0 */
  if (FileReadPlainTriDiagonalMatrixType(&G0, HSet.G0FileName))
  {
    printf("Fatal Error in main\n Failed to read G0 from file %s\n", HSet.G0FileName);
    return 1;
  };
  strcpy(G0.Name, "G0");
  /* read G_Bath_iStart and psi_iStart */
  if (HSet.iStart == 0)
  {
    sprintf(GBathFileName, "%s000", HSet.GRootFileName);
    sprintf(PsiFileName, "%s000", HSet.WaveFunctionRootFileName);
  }
  else if (HSet.iStart < 10)
  {
    sprintf(GBathFileName, "%s00%1i", HSet.GRootFileName, HSet.iStart);
    sprintf(PsiFileName, "%s00%1i", HSet.WaveFunctionRootFileName, HSet.iStart);
  }
  else if (HSet.iStart < 100)
  {
    sprintf(GBathFileName, "%s0%2i", HSet.GRootFileName, HSet.iStart);
    sprintf(PsiFileName, "%s0%2i", HSet.WaveFunctionRootFileName, HSet.iStart);
  }
  else
  {
    sprintf(GBathFileName, "%s%i", HSet.GRootFileName, HSet.iStart);
    sprintf(PsiFileName, "%s%i", HSet.WaveFunctionRootFileName, HSet.iStart);
  }
  if (FileReadPlainTriDiagonalMatrixType(&GBath, GBathFileName))
  {
    printf("Fatal Error in main\n Failed to read GBath from file %s !\n", GBathFileName);
    return 1;
  };
  if (file_exists(PsiFileName))
  {
    Psi.NBitsKey = HSet.NBitsKey;
    if (FileReadWaveFunctionPlainBinary(&Psi, PsiFileName))
    {
      printf("Fatal Error in main\n Failed to read Psi from file %s !\n", PsiFileName);
      return 1;
    };
    printf("\nRead the starting wave function from the file %s\n", PsiFileName);
    sscanf(Psi.Name, "%s %i %i", Buffer, &PsiNBathC, &PsiNBathV);
    printf("Wave-function \"%s\" found has %i conduction bath sites, %i valence bath sites and %i Fermions in basis\n",
           Psi.Name, PsiNBathC, PsiNBathV, Psi.NFermion);
    if ((2 * (PsiNBathC + PsiNBathV + 1)) != Psi.NFermion)
    {
      printf("The total number of Fermions is incompatible witht the bath sizes\n");
      printf("This wave-function can not be used to start from\n");
      return 1;
    }
  }
  else
  {
    printf("\nWave-function file: '%s' does not excist\n Started from new function\n", PsiFileName);
    CreateWaveFunctionDMFT(&Psi, HSet.NBath, HSet.NBath, HSet.NBitsKey);
  }
  strcpy(GBath.Name, "Gbath");
  /* do the DMFT loops */
  for (iOuter = HSet.iStart; iOuter < HSet.iStart + HSet.Ni; iOuter++)
  {
    BeginTiming("DMFTLoop", Time);
    printf("\n\n***********************************************************************\n");
    printf("Start of loop %i from %i to %i\n", iOuter + 1, HSet.iStart, HSet.iStart + HSet.Ni);
    printf("***********************************************************************\n\n\n");
    TriDiagonalMatrixToDoubleTriDiagonalMatrix(GBath, &T, &GAnd, Ef, epsilonAnd);
    FreeMatrix(&T);
    if (GAnd.con.N < HSet.NBath)
    {
      printf("Warning: GAnd.con.N(%i)<HSet.NBath(%i)\n", GAnd.con.N, HSet.NBath);
    }
    else
    {
      GAnd.con.N = HSet.NBath;
    }
    if (GAnd.val.N < HSet.NBath)
    {
      printf("Warning: GAnd.val.N(%i)<HSet.NBath(%i)\n", GAnd.val.N, HSet.NBath);
    }
    else
    {
      GAnd.val.N = HSet.NBath;
    }
    mu = -HSet.U / 2;
    DoubleTriDiagonalMatrixToAndersonImpurityOperator(&HAnd, GAnd, GAnd, HSet.U, mu);
    /* test if psi has correct dimensions and correct if needed */
    sscanf(Psi.Name, "%s %i %i", Buffer, &PsiNBathC, &PsiNBathV);
    if (PsiNBathC != GAnd.con.N)
    {
      printf("\nBefore diagonalizing:\nThe wavefunction has %i conduction bath sites\nThe Anderson Hamiltonian %i\n",
             PsiNBathC, GAnd.con.N);
      if (PsiNBathC > GAnd.con.N)
      {
        printf("Removing %i unoccupied fermions from the wave-function\n", PsiNBathC - GAnd.con.N);
        for (i = 0; i < PsiNBathC - GAnd.con.N; i++)
        {
          WaveFunctionRemoveFermion(Psi, &PsiTemp, Psi.NFermion - 1);
          WaveFunctionFree(&Psi);
          WaveFunctionRemoveFermion(PsiTemp, &Psi, 2 + PsiNBathV + PsiNBathC - i - 1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      else
      {
        printf("Adding %i unoccupied fermions to the wave-function\n", GAnd.con.N - PsiNBathC);
        for (i = 0; i < GAnd.con.N - PsiNBathC; i++)
        {
          WaveFunctionAddFermion(Psi, &PsiTemp, Psi.NFermion, 0);
          WaveFunctionFree(&Psi);
          WaveFunctionAddFermion(PsiTemp, &Psi, 2 + PsiNBathV + PsiNBathC + i, 0);
          WaveFunctionFree(&PsiTemp);
        }
      }
    }
    if (PsiNBathV != GAnd.val.N)
    {
      printf("\nBefore diagonalizing:\nThe wavefunction has %i valence bath sites\nThe Anderson Hamiltonian %i\n",
             PsiNBathV, GAnd.val.N);
      if (PsiNBathV > GAnd.val.N)
      {
        printf("Removing %i occupied fermions from the wave-function\n", PsiNBathV - GAnd.val.N);
        for (i = 0; i < PsiNBathV - GAnd.val.N; i++)
        {
          WaveFunctionRemoveFermion(Psi, &PsiTemp, Psi.NFermion - GAnd.con.N - 1);
          WaveFunctionFree(&Psi);
          WaveFunctionRemoveFermion(PsiTemp, &Psi, 2 + PsiNBathV - i - 1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      else
      {
        printf("Adding %i occupied fermions to the wave-function\n", GAnd.val.N - PsiNBathV);
        for (i = 0; i < GAnd.val.N - PsiNBathV; i++)
        {
          WaveFunctionAddFermion(Psi, &PsiTemp, Psi.NFermion - GAnd.val.N, 1);
          WaveFunctionFree(&Psi);
          WaveFunctionAddFermion(PsiTemp, &Psi, 2 + PsiNBathV + i, 1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      printf("Corrected basis\n");
      PrintWaveFunctionReadableTruncated(Psi, 0.1);
    }
    sprintf(Psi.Name, "DMFT %i %i", GAnd.con.N, GAnd.val.N);
    /* create ground-state wave-funciton */
    BeginTiming("Lanczos GroundState", Time);
    LanczosGroundState(HAnd, &Psi, &HSet.zeroLancWavefunction, HSet.zeroLancEnergy, HSet.U, &Psi, 0);
    LanczosGroundState(HAnd, &Psi, &HSet.zeroLancWavefunction, HSet.zeroLancEnergy, HSet.U, &Psi, 0);
    EndTiming("Lanczos GroundState", Time);
    PsiHermitianOperatorPsi(HAnd, Psi, &PsiEnergy);
    PsiHermitianOperatorSqrPsi(HAnd, Psi, &PsiQuality);
    PsiQuality /= Sqr(PsiEnergy);
    PrintWaveFunctionReadableTruncated(Psi, 0.1);
    CalculatePrefactorHistogram(PsiHist, Psi, 10.0, 20);
    printf("\nHistogram of prefactors in Wave-function\n");
    for (i = 0; i < 19; i++)
    {
      printf("10^(-%2i) to 10^(-%2i) | %22.15E\n", i, i + 1, PsiHist[i]);
    }
    printf("10^(-%2i) to 0        | %22.15E\n", i, PsiHist[19]);
    printf("\n");
    printf("Energy                  <psi|H^2|psi>/(<psi|H|psi>^2)\n");
    printf("%22.15E  %22.15E\n", PsiEnergy, PsiQuality);
    PsiRho = calloc(Psi.NFermion, sizeof(double));
    if (PsiRho == NULL)
    {
      printf("calloc failed in main for PsiRho\n Fatal Error\n");
      fflush(stdout);
      return 1;
    }
    CalculateDensity(PsiRho, Psi);
    printf("Occupation \n");
    for (i = 0; i < Psi.NFermion; i++)
    {
      printf("%5i %22.15E\n", i, PsiRho[i]);
    }
    free(PsiRho);
    /* create spin up and spin down PES and IPES */
    CreatePESdnOperatorDMFT(&OPESdn, HAnd.NFermion);
    CreatePESupOperatorDMFT(&OPESup, HAnd.NFermion);
    CreateIPESdnOperatorDMFT(&OIPESdn, HAnd.NFermion);
    CreateIPESupOperatorDMFT(&OIPESup, HAnd.NFermion);
    GclMup.N = HSet.NTriDiag;
    InitTriDiagonalMatrix(&GclMup);
    GclMdn.N = HSet.NTriDiag;
    InitTriDiagonalMatrix(&GclMdn);
    GclPup.N = HSet.NTriDiag;
    InitTriDiagonalMatrix(&GclPup);
    GclPdn.N = HSet.NTriDiag;
    InitTriDiagonalMatrix(&GclPdn);
    BeginTiming("Lanczos Tri Diagonalization PES up", Time);
    OperatorPsi(OPESup, Psi, &OPsi);
    LanczosTriDiagonalize(HAnd, OPsi, HSet.zeroLancWavefunction, &GclMup, HSet.Gamma);
    EndTiming("Lanczos Tri Diagonalization PES up", Time);
    BeginTiming("Lanczos Tri Diagonalization PES dn", Time);
    OperatorPsi(OPESdn, Psi, &OPsi);
    LanczosTriDiagonalize(HAnd, OPsi, HSet.zeroLancWavefunction, &GclMdn, HSet.Gamma);
    EndTiming("Lanczos Tri Diagonalization PES dn", Time);
    BeginTiming("Lanczos Tri Diagonalization IPES up", Time);
    OperatorPsi(OIPESup, Psi, &OPsi);
    LanczosTriDiagonalize(HAnd, OPsi, HSet.zeroLancWavefunction, &GclPup, HSet.Gamma);
    EndTiming("Lanczos Tri Diagonalization IPES up", Time);
    BeginTiming("Lanczos Tri Diagonalization IPES dn", Time);
    OperatorPsi(OIPESdn, Psi, &OPsi);
    LanczosTriDiagonalize(HAnd, OPsi, HSet.zeroLancWavefunction, &GclPdn, HSet.Gamma);
    EndTiming("Lanczos Tri Diagonalization IPES dn", Time);
    /* G math */
    printf("Start G math\n");
    ShiftTriDiagonalMatrix(&GclMdn, -PsiEnergy);
    ShiftTriDiagonalMatrix(&GclMup, -PsiEnergy);
    ShiftTriDiagonalMatrix(&GclPdn, -PsiEnergy);
    ShiftTriDiagonalMatrix(&GclPup, -PsiEnergy);
    for (i = 1; i <= GclMdn.N; i++)
    {
      GclMdn.a[i] *= -1;
    }
    for (i = 1; i <= GclMup.N; i++)
    {
      GclMup.a[i] *= -1;
    }
    strcpy(GclPup.Name, "Gcluster-Plus-Up");
    strcpy(GclPdn.Name, "Gcluster-Plus-Dn");
    strcpy(GclMup.Name, "Gcluster-Min-Up");
    strcpy(GclMdn.Name, "Gcluster-Min-Dn");
    printf("  Start SumGreensfunctionsTriDiagonalRepresentation\n");
    fflush(stdout);
    SumGreensfunctionsTriDiagonalRepresentation(GclPup, GclMup, &Gclup);
    strcpy(Gclup.Name, "Gcluster-Up");
    SumGreensfunctionsTriDiagonalRepresentation(GclPdn, GclMdn, &Gcldn);
    strcpy(Gcldn.Name, "Gcluster-Dn");
    SumGreensfunctionsTriDiagonalRepresentation(Gclup, Gcldn, &Gcl);
    strcpy(Gcl.Name, "Gcluster");
    Gcl.b[0] = 1.0;
    SymmetrizeTridiagonalMatrix(&Gcl);
    Gcl.N = Min(Gcl.N, Max(4 * HSet.NTriDiag, G0.N));
    printf("  Start SelfEnergyTriDiagonalRepresentation\n");
    SelfEnergyTriDiagonalRepresentation(GBath, Gcl, &Sigma);
    Sigma.N = Min(Sigma.N, Max(4 * HSet.NTriDiag, G0.N));
    strcpy(Sigma.Name, "Sigma");
    SymmetrizeTridiagonalMatrix(&Sigma);
    Sigma.N = Min(Sigma.N, Max(4 * HSet.NTriDiag, G0.N));
    printf("  Start GLocalTriDiagonalRepresentation\n");
    GLocalTriDiagonalRepresentation(G0, Sigma, &GLocal, HSet.NTriDiag);
    GLocal.N = Min(GLocal.N, Max(4 * HSet.NTriDiag, G0.N));
    strcpy(GLocal.Name, "Glocal");
    SymmetrizeTridiagonalMatrix(&GLocal);
    GLocal.N = Min(GLocal.N, Max(4 * HSet.NTriDiag, G0.N));
    printf("  Start GBathNewTriDiagonalRepresentation\n");
    GBathNewTriDiagonalRepresentation(GLocal, Sigma, &GBathNew);
    GBathNew.N = Min(GBathNew.N, Max(4 * HSet.NTriDiag, G0.N));
    strcpy(GBathNew.Name, "Gbath-New");
    SymmetrizeTridiagonalMatrix(&GBathNew);
    GBathNew.N = Min(GBathNew.N, Max(4 * HSet.NTriDiag, G0.N));
    /* calculate GBath for the next loop by mixing GBath with GBathNew */
    printf("Start mixing GBath with GBathNew\n");
    TriDiagonalMatrixCopy(&GBathOld, GBath);
    strcpy(GBathOld.Name, "Gbath-Old");
    FreeTriDiagonalMatrix(&GBath);
    GBathOld.b[0] = sqrt(HSet.Mix);
    GBathNew.b[0] = sqrt(1 - HSet.Mix);
    SumGreensfunctionsTriDiagonalRepresentation(GBathNew, GBathOld, &GBath);
    strcpy(GBath.Name, "Gbath");
    GBathOld.b[0] = 1;
    GBathNew.b[0] = 1;
    SymmetrizeTridiagonalMatrix(&GBath);
    GBath.N = Min(GBath.N, Max(4 * HSet.NTriDiag, G0.N));
    /* create spectra */
    printf("Start create spectra\n");
    Spectra[0].N = HSet.NE;
    Spectra[0].Emin = HSet.Emin;
    Spectra[0].Emax = HSet.Emax;
    Spectra[0].dE = (HSet.Emax - HSet.Emin) / HSet.NE;
    Spectra[0].EminPole = HSet.Emin;
    Spectra[0].EmaxPole = HSet.Emax;
    Spectra[0].Gamma = HSet.Gamma;
    for (i = 0; i < 10; i++)
    {
      SpectrumAllocate(Spectra[0], &Spectra[i]);
    }
    TriDiagonalMatrixToSpectrum(0, GBath, &(Spectra[0]));
    TriDiagonalMatrixToSpectrum(0, GBathOld, &(Spectra[1]));
    TriDiagonalMatrixToSpectrum(0, GBathNew, &(Spectra[2]));
    DoubleTriDiagonalMatrixToSpectrum(0, GAnd, &(Spectra[3]));
    TriDiagonalMatrixToSpectrum(0, G0, &(Spectra[4]));
    TriDiagonalMatrixToSpectrum(0, Sigma, &(Spectra[5]));
    TriDiagonalMatrixToSpectrum(0, GLocal, &(Spectra[6]));
    TriDiagonalMatrixToSpectrum(0, Gclup, &(Spectra[7]));
    TriDiagonalMatrixToSpectrum(0, Gcldn, &(Spectra[8]));
    TriDiagonalMatrixToSpectrum(0, Gcl, &(Spectra[9]));
    /* write to disk */
    printf("Start write to dist\n");
    if (iOuter == 0)
    {
      sprintf(GBathFileName, "%s001", HSet.GRootFileName);
      sprintf(PsiFileName, "%s001", HSet.WaveFunctionRootFileName);
      sprintf(SpecFileName, "%s001", HSet.SpecRootFileName);
    }
    else if (iOuter < 9)
    {
      sprintf(GBathFileName, "%s00%1i", HSet.GRootFileName, iOuter + 1);
      sprintf(PsiFileName, "%s00%1i", HSet.WaveFunctionRootFileName, iOuter + 1);
      sprintf(SpecFileName, "%s00%1i", HSet.SpecRootFileName, iOuter + 1);
    }
    else if (iOuter < 99)
    {
      sprintf(GBathFileName, "%s0%2i", HSet.GRootFileName, iOuter + 1);
      sprintf(PsiFileName, "%s0%2i", HSet.WaveFunctionRootFileName, iOuter + 1);
      sprintf(SpecFileName, "%s0%2i", HSet.SpecRootFileName, iOuter + 1);
    }
    else if (iOuter < 999)
    {
      sprintf(GBathFileName, "%s%3i", HSet.GRootFileName, iOuter + 1);
      sprintf(PsiFileName, "%s%3i", HSet.WaveFunctionRootFileName, iOuter + 1);
      sprintf(SpecFileName, "%s%3i", HSet.SpecRootFileName, iOuter + 1);
    }
    FilePrintSpectra(Spectra, 10, SpecFileName, "w");
    FilePrintPlainTriDiagonalMatrix(GBath, GBathFileName, "w");
    FilePrintWaveFunctionPlainBinary(Psi, PsiFileName, "w");
    /* free variables */
    printf("Start free variables\n");
    FreeOperator(&HAnd);
    FreeOperator(&OPESdn);
    FreeOperator(&OPESup);
    FreeOperator(&OIPESdn);
    FreeOperator(&OIPESup);
    FreeTriDiagonalMatrix(&GclMup);
    FreeTriDiagonalMatrix(&GclMdn);
    FreeTriDiagonalMatrix(&GclPup);
    FreeTriDiagonalMatrix(&GclPdn);
    FreeTriDiagonalMatrix(&Gcl);
    FreeTriDiagonalMatrix(&GLocal);
    FreeTriDiagonalMatrix(&Sigma);
    FreeTriDiagonalMatrix(&GBathNew);
    FreeTriDiagonalMatrix(&GBathOld);
    FreeDoubleTriDiagonalMatrix(&GAnd);
    EndTiming("DMFTLoop", Time);
    PrintTiming(*Time);
  }
  return 0;
};

int oldmainNoDysonSecond(DMFTSettings HSet, TimeType * Time)
{
  WaveFunctionType Psi, OPsi, PsiTemp;
  TriDiagonalMatrixType G0, GBath, GclPup, GclPdn, GclMup, GclMdn, Gclup, Gcldn, Gcl, Sigma, GLocal, GBathNew, GBathOld;
  SecondDoubleTriDiagonalMatrixType GAnd;
  OperatorType HAnd, OPESdn, OPESup, OIPESdn, OIPESup;
  MatrixType T;
  SpectraType Spectra[10];
  char   GBathFileName[255], PsiFileName[255], SpecFileName[25], Buffer[255];
  uint32_t iOuter, i, PsiNBathV, PsiNBathC;
  double Ef = 0.0, epsilonAnd = 0.0000001, mu;
  double PsiEnergy, PsiQuality, PsiHist[20], *PsiRho;
  /* read G0 */
  if (FileReadPlainTriDiagonalMatrixType(&G0, HSet.G0FileName))
  {
    printf("Fatal Error in main\n Failed to read G0 from file %s\n", HSet.G0FileName);
    return 1;
  };
  strcpy(G0.Name, "G0");
  /* read G_Bath_iStart and psi_iStart */
  if (HSet.iStart == 0)
  {
    sprintf(GBathFileName, "%s000", HSet.GRootFileName);
    sprintf(PsiFileName, "%s000", HSet.WaveFunctionRootFileName);
  }
  else if (HSet.iStart < 10)
  {
    sprintf(GBathFileName, "%s00%1i", HSet.GRootFileName, HSet.iStart);
    sprintf(PsiFileName, "%s00%1i", HSet.WaveFunctionRootFileName, HSet.iStart);
  }
  else if (HSet.iStart < 100)
  {
    sprintf(GBathFileName, "%s0%2i", HSet.GRootFileName, HSet.iStart);
    sprintf(PsiFileName, "%s0%2i", HSet.WaveFunctionRootFileName, HSet.iStart);
  }
  else
  {
    sprintf(GBathFileName, "%s%i", HSet.GRootFileName, HSet.iStart);
    sprintf(PsiFileName, "%s%i", HSet.WaveFunctionRootFileName, HSet.iStart);
  }
  if (FileReadPlainTriDiagonalMatrixType(&GBath, GBathFileName))
  {
    printf("Fatal Error in main\n Failed to read GBath from file %s !\n", GBathFileName);
    return 1;
  };
  if (file_exists(PsiFileName))
  {
    Psi.NBitsKey = HSet.NBitsKey;
    if (FileReadWaveFunctionPlainBinary(&Psi, PsiFileName))
    {
      printf("Fatal Error in main\n Failed to read Psi from file %s !\n", PsiFileName);
      return 1;
    };
    printf("\nRead the starting wave function from the file %s\n", PsiFileName);
    sscanf(Psi.Name, "%s %i %i", Buffer, &PsiNBathC, &PsiNBathV);
    printf("Wave-function \"%s\" found has %i conduction bath sites, %i valence bath sites and %i Fermions in basis\n",
           Psi.Name, PsiNBathC, PsiNBathV, Psi.NFermion);
    if ((2 * (PsiNBathC + PsiNBathV + 2)) != Psi.NFermion)
    {
      printf("The total number of Fermions is incompatible witht the bath sizes\n");
      printf("This wave-function can not be used to start from\n");
      return 1;
    }
  }
  else
  {
    printf("\nWave-function file: '%s' does not excist\n Started from new function\n", PsiFileName);
    CreateWaveFunctionDMFTSecond(&Psi, HSet.NBath, HSet.NBath, HSet.NBitsKey);
  }
  strcpy(GBath.Name, "Gbath");
  /* do the DMFT loops */
  for (iOuter = HSet.iStart; iOuter < HSet.iStart + HSet.Ni; iOuter++)
  {
    BeginTiming("DMFTLoop", Time);
    printf("\n\n***********************************************************************\n");
    printf("Start of loop %i from %i to %i\n", iOuter + 1, HSet.iStart, HSet.iStart + HSet.Ni);
    printf("***********************************************************************\n\n\n");
    TriDiagonalMatrixToSecondDoubleTriDiagonalMatrix(GBath, &T, &GAnd, Ef, epsilonAnd);

    PrintTriDiagonalMatrix(GBath);
    PrintDoubleTriDiagonalMatrix(GAnd.Dt);

    FreeMatrix(&T);
    if (GAnd.Dt.con.N < HSet.NBath)
    {
      printf("Warning: GAnd.Dt.con.N(%i)<HSet.NBath(%i)\n", GAnd.Dt.con.N, HSet.NBath);
    }
    else
    {
      GAnd.Dt.con.N = HSet.NBath;
    }
    if (GAnd.Dt.val.N < HSet.NBath)
    {
      printf("Warning: GAnd.Dt.val.N(%i)<HSet.NBath(%i)\n", GAnd.Dt.val.N, HSet.NBath);
    }
    else
    {
      GAnd.Dt.val.N = HSet.NBath;
    }
    mu = -HSet.U / 2;
    SecondDoubleTriDiagonalMatrixToAndersonImpurityOperator(&HAnd, GAnd, GAnd, HSet.U, mu);
    PrintOperatorReadable(HAnd);
    /* test if psi has correct dimensions and correct if needed */
    sscanf(Psi.Name, "%s %i %i", Buffer, &PsiNBathC, &PsiNBathV);
    if (PsiNBathC != GAnd.Dt.con.N)
    {
      printf("\nBefore diagonalizing:\nThe wavefunction has %i conduction bath sites\nThe Anderson Hamiltonian %i\n",
             PsiNBathC, GAnd.Dt.con.N);
      if (PsiNBathC > GAnd.Dt.con.N)
      {
        printf("Removing %i unoccupied fermions from the wave-function\n", PsiNBathC - GAnd.Dt.con.N);
        for (i = 0; i < PsiNBathC - GAnd.Dt.con.N; i++)
        {
          WaveFunctionRemoveFermion(Psi, &PsiTemp, Psi.NFermion - 1);
          WaveFunctionFree(&Psi);
          WaveFunctionRemoveFermion(PsiTemp, &Psi, 4 + PsiNBathV + PsiNBathC - i - 1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      else
      {
        printf("Adding %i unoccupied fermions to the wave-function\n", GAnd.Dt.con.N - PsiNBathC);
        for (i = 0; i < GAnd.Dt.con.N - PsiNBathC; i++)
        {
          WaveFunctionAddFermion(Psi, &PsiTemp, Psi.NFermion, 0);
          WaveFunctionFree(&Psi);
          WaveFunctionAddFermion(PsiTemp, &Psi, 4 + PsiNBathV + PsiNBathC + i, 0);
          WaveFunctionFree(&PsiTemp);
        }
      }
    }
    if (PsiNBathV != GAnd.Dt.val.N)
    {
      printf("\nBefore diagonalizing:\nThe wavefunction has %i valence bath sites\nThe Anderson Hamiltonian %i\n",
             PsiNBathV, GAnd.Dt.val.N);
      if (PsiNBathV > GAnd.Dt.val.N)
      {
        printf("Removing %i occupied fermions from the wave-function\n", PsiNBathV - GAnd.Dt.val.N);
        for (i = 0; i < PsiNBathV - GAnd.Dt.val.N; i++)
        {
          WaveFunctionRemoveFermion(Psi, &PsiTemp, Psi.NFermion - GAnd.Dt.con.N - 1);
          WaveFunctionFree(&Psi);
          WaveFunctionRemoveFermion(PsiTemp, &Psi, 4 + PsiNBathV - i - 1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      else
      {
        printf("Adding %i occupied fermions to the wave-function\n", GAnd.Dt.val.N - PsiNBathV);
        for (i = 0; i < GAnd.Dt.val.N - PsiNBathV; i++)
        {
          WaveFunctionAddFermion(Psi, &PsiTemp, Psi.NFermion - GAnd.Dt.val.N, 1);
          WaveFunctionFree(&Psi);
          WaveFunctionAddFermion(PsiTemp, &Psi, 4 + PsiNBathV + i, 1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      printf("Corrected basis\n");
      PrintWaveFunctionReadableTruncated(Psi, 0.1);
    }
    sprintf(Psi.Name, "DMFT %i %i", GAnd.Dt.con.N, GAnd.Dt.val.N);
    /* create ground-state wave-funciton */
    BeginTiming("Lanczos GroundState", Time);
    LanczosGroundState(HAnd, &Psi, &HSet.zeroLancWavefunction, HSet.zeroLancEnergy, HSet.U, &Psi, 0);
    LanczosGroundState(HAnd, &Psi, &HSet.zeroLancWavefunction, HSet.zeroLancEnergy, HSet.U, &Psi, 0);
    EndTiming("Lanczos GroundState", Time);
    PsiHermitianOperatorPsi(HAnd, Psi, &PsiEnergy);
    PsiHermitianOperatorSqrPsi(HAnd, Psi, &PsiQuality);
    PsiQuality /= Sqr(PsiEnergy);
    PrintWaveFunctionReadableTruncated(Psi, 0.1);
    CalculatePrefactorHistogram(PsiHist, Psi, 10.0, 20);
    printf("\nHistogram of prefactors in Wave-function\n");
    for (i = 0; i < 19; i++)
    {
      printf("10^(-%2i) to 10^(-%2i) | %22.15E\n", i, i + 1, PsiHist[i]);
    }
    printf("10^(-%2i) to 0        | %22.15E\n", i, PsiHist[19]);
    printf("\n");
    printf("Energy                  <psi|H^2|psi>/(<psi|H|psi>^2)\n");
    printf("%22.15E  %22.15E\n", PsiEnergy, PsiQuality);
    PsiRho = calloc(Psi.NFermion, sizeof(double));
    if (PsiRho == NULL)
    {
      printf("calloc failed in main for PsiRho\n Fatal Error\n");
      fflush(stdout);
      return 1;
    }
    CalculateDensity(PsiRho, Psi);
    printf("Occupation \n");
    for (i = 0; i < Psi.NFermion; i++)
    {
      printf("%5i %22.15E\n", i, PsiRho[i]);
    }
    free(PsiRho);
    /* create spin up and spin down PES and IPES */
    CreatePESdnOperatorDMFT(&OPESdn, HAnd.NFermion);
    CreatePESupOperatorDMFT(&OPESup, HAnd.NFermion);
    CreateIPESdnOperatorDMFT(&OIPESdn, HAnd.NFermion);
    CreateIPESupOperatorDMFT(&OIPESup, HAnd.NFermion);
    GclMup.N = HSet.NTriDiag;
    InitTriDiagonalMatrix(&GclMup);
    GclMdn.N = HSet.NTriDiag;
    InitTriDiagonalMatrix(&GclMdn);
    GclPup.N = HSet.NTriDiag;
    InitTriDiagonalMatrix(&GclPup);
    GclPdn.N = HSet.NTriDiag;
    InitTriDiagonalMatrix(&GclPdn);
    BeginTiming("Lanczos Tri Diagonalization PES up", Time);
    OperatorPsi(OPESup, Psi, &OPsi);
    LanczosTriDiagonalize(HAnd, OPsi, HSet.zeroLancWavefunction, &GclMup, HSet.Gamma);
    EndTiming("Lanczos Tri Diagonalization PES up", Time);
    BeginTiming("Lanczos Tri Diagonalization PES dn", Time);
    OperatorPsi(OPESdn, Psi, &OPsi);
    LanczosTriDiagonalize(HAnd, OPsi, HSet.zeroLancWavefunction, &GclMdn, HSet.Gamma);
    EndTiming("Lanczos Tri Diagonalization PES dn", Time);
    BeginTiming("Lanczos Tri Diagonalization IPES up", Time);
    OperatorPsi(OIPESup, Psi, &OPsi);
    LanczosTriDiagonalize(HAnd, OPsi, HSet.zeroLancWavefunction, &GclPup, HSet.Gamma);
    EndTiming("Lanczos Tri Diagonalization IPES up", Time);
    BeginTiming("Lanczos Tri Diagonalization IPES dn", Time);
    OperatorPsi(OIPESdn, Psi, &OPsi);
    LanczosTriDiagonalize(HAnd, OPsi, HSet.zeroLancWavefunction, &GclPdn, HSet.Gamma);
    EndTiming("Lanczos Tri Diagonalization IPES dn", Time);
    /* G math */
    printf("Start G math\n");
    ShiftTriDiagonalMatrix(&GclMdn, -PsiEnergy);
    ShiftTriDiagonalMatrix(&GclMup, -PsiEnergy);
    ShiftTriDiagonalMatrix(&GclPdn, -PsiEnergy);
    ShiftTriDiagonalMatrix(&GclPup, -PsiEnergy);
    for (i = 1; i <= GclMdn.N; i++)
    {
      GclMdn.a[i] *= -1;
    }
    for (i = 1; i <= GclMup.N; i++)
    {
      GclMup.a[i] *= -1;
    }
    strcpy(GclPup.Name, "Gcluster-Plus-Up");
    strcpy(GclPdn.Name, "Gcluster-Plus-Dn");
    strcpy(GclMup.Name, "Gcluster-Min-Up");
    strcpy(GclMdn.Name, "Gcluster-Min-Dn");
    printf("  Start SumGreensfunctionsTriDiagonalRepresentation\n");
    fflush(stdout);
    SumGreensfunctionsTriDiagonalRepresentation(GclPup, GclMup, &Gclup);
    strcpy(Gclup.Name, "Gcluster-Up");
    SumGreensfunctionsTriDiagonalRepresentation(GclPdn, GclMdn, &Gcldn);
    strcpy(Gcldn.Name, "Gcluster-Dn");
    SumGreensfunctionsTriDiagonalRepresentation(Gclup, Gcldn, &Gcl);
    strcpy(Gcl.Name, "Gcluster");
    Gcl.b[0] = 1.0;
    // SymmetrizeTridiagonalMatrix(&Gcl);
    Gcl.N = Min(Gcl.N, Max(4 * HSet.NTriDiag, G0.N));
    printf("  Start SelfEnergyTriDiagonalRepresentation\n");
    SelfEnergyTriDiagonalRepresentation(GBath, Gcl, &Sigma);
    Sigma.N = Min(Sigma.N, Max(4 * HSet.NTriDiag, G0.N));
    strcpy(Sigma.Name, "Sigma");
    // SymmetrizeTridiagonalMatrix(&Sigma);
    Sigma.N = Min(Sigma.N, Max(4 * HSet.NTriDiag, G0.N));
    printf("  Start GLocalTriDiagonalRepresentation\n");
    GLocalTriDiagonalRepresentation(G0, Sigma, &GLocal, HSet.NTriDiag);
    GLocal.N = Min(GLocal.N, Max(4 * HSet.NTriDiag, G0.N));
    strcpy(GLocal.Name, "Glocal");
    // SymmetrizeTridiagonalMatrix(&GLocal);
    GLocal.N = Min(GLocal.N, Max(4 * HSet.NTriDiag, G0.N));
    printf("  Start GBathNewTriDiagonalRepresentation\n");
    GBathNewTriDiagonalRepresentation(GLocal, Sigma, &GBathNew);
    GBathNew.N = Min(GBathNew.N, Max(4 * HSet.NTriDiag, G0.N));
    strcpy(GBathNew.Name, "Gbath-New");
    // SymmetrizeTridiagonalMatrix(&GBathNew);
    GBathNew.N = Min(GBathNew.N, Max(4 * HSet.NTriDiag, G0.N));
    /* calculate GBath for the next loop by mixing GBath with GBathNew */
    printf("Start mixing GBath with GBathNew\n");
    TriDiagonalMatrixCopy(&GBathOld, GBath);
    strcpy(GBathOld.Name, "Gbath-Old");
    FreeTriDiagonalMatrix(&GBath);
    GBathOld.b[0] = sqrt(HSet.Mix);
    GBathNew.b[0] = sqrt(1 - HSet.Mix);
    SumGreensfunctionsTriDiagonalRepresentation(GBathNew, GBathOld, &GBath);
    strcpy(GBath.Name, "Gbath");
    GBathOld.b[0] = 1;
    GBathNew.b[0] = 1;
    //  SymmetrizeTridiagonalMatrix(&GBath);
    GBath.N = Min(GBath.N, Max(4 * HSet.NTriDiag, G0.N));
    /* create spectra */
    printf("Start create spectra\n");
    Spectra[0].N = HSet.NE;
    Spectra[0].Emin = HSet.Emin;
    Spectra[0].Emax = HSet.Emax;
    Spectra[0].dE = (HSet.Emax - HSet.Emin) / HSet.NE;
    Spectra[0].EminPole = HSet.Emin;
    Spectra[0].EmaxPole = HSet.Emax;
    Spectra[0].Gamma = HSet.Gamma;
    for (i = 0; i < 10; i++)
    {
      SpectrumAllocate(Spectra[0], &Spectra[i]);
    }
    TriDiagonalMatrixToSpectrum(0, GBath, &(Spectra[0]));
    TriDiagonalMatrixToSpectrum(0, GBathOld, &(Spectra[1]));
    TriDiagonalMatrixToSpectrum(0, GBathNew, &(Spectra[2]));
    SecondDoubleTriDiagonalMatrixToSpectrum(0, GAnd, &(Spectra[3]));
    TriDiagonalMatrixToSpectrum(0, G0, &(Spectra[4]));
    TriDiagonalMatrixToSpectrum(0, Sigma, &(Spectra[5]));
    TriDiagonalMatrixToSpectrum(0, GLocal, &(Spectra[6]));
    TriDiagonalMatrixToSpectrum(0, GclPup, &(Spectra[7]));
    TriDiagonalMatrixToSpectrum(0, GclMup, &(Spectra[8]));
    TriDiagonalMatrixToSpectrum(0, Gcl, &(Spectra[9]));
    /* write to disk */
    printf("Start write to dist\n");
    if (iOuter == 0)
    {
      sprintf(GBathFileName, "%s001", HSet.GRootFileName);
      sprintf(PsiFileName, "%s001", HSet.WaveFunctionRootFileName);
      sprintf(SpecFileName, "%s001", HSet.SpecRootFileName);
    }
    else if (iOuter < 9)
    {
      sprintf(GBathFileName, "%s00%1i", HSet.GRootFileName, iOuter + 1);
      sprintf(PsiFileName, "%s00%1i", HSet.WaveFunctionRootFileName, iOuter + 1);
      sprintf(SpecFileName, "%s00%1i", HSet.SpecRootFileName, iOuter + 1);
    }
    else if (iOuter < 99)
    {
      sprintf(GBathFileName, "%s0%2i", HSet.GRootFileName, iOuter + 1);
      sprintf(PsiFileName, "%s0%2i", HSet.WaveFunctionRootFileName, iOuter + 1);
      sprintf(SpecFileName, "%s0%2i", HSet.SpecRootFileName, iOuter + 1);
    }
    else if (iOuter < 999)
    {
      sprintf(GBathFileName, "%s%3i", HSet.GRootFileName, iOuter + 1);
      sprintf(PsiFileName, "%s%3i", HSet.WaveFunctionRootFileName, iOuter + 1);
      sprintf(SpecFileName, "%s%3i", HSet.SpecRootFileName, iOuter + 1);
    }
    FilePrintSpectra(Spectra, 10, SpecFileName, "w");
    FilePrintPlainTriDiagonalMatrix(GBath, GBathFileName, "w");
    FilePrintWaveFunctionPlainBinary(Psi, PsiFileName, "w");
    /* free variables */
    printf("Start free variables\n");
    FreeOperator(&HAnd);
    FreeOperator(&OPESdn);
    FreeOperator(&OPESup);
    FreeOperator(&OIPESdn);
    FreeOperator(&OIPESup);
    FreeTriDiagonalMatrix(&GclMup);
    FreeTriDiagonalMatrix(&GclMdn);
    FreeTriDiagonalMatrix(&GclPup);
    FreeTriDiagonalMatrix(&GclPdn);
    FreeTriDiagonalMatrix(&Gcl);
    FreeTriDiagonalMatrix(&GLocal);
    FreeTriDiagonalMatrix(&Sigma);
    FreeTriDiagonalMatrix(&GBathNew);
    FreeTriDiagonalMatrix(&GBathOld);
    FreeDoubleTriDiagonalMatrix(&(GAnd.Dt));
    EndTiming("DMFTLoop", Time);
    PrintTiming(*Time);
  }
  return 0;
};

int oldmainDysonSecond(DMFTSettings HSet, TimeType * Time)
{
  WaveFunctionType Psi, OPsi, PsiTemp;
  TriDiagonalMatrixType G0, GBath, GclPup, GclPdn, GclMup, GclMdn, Gclup, Gcldn, Gcl, Sigma, GLocal, GBathNew, GBathOld,
    gPerturbation[4];
  SecondDoubleTriDiagonalMatrixType GAnd, GAndTrunk;
  OperatorType HAnd, OPESdn, OPESup, OIPESdn, OIPESup, VPerturbation[4];
  MatrixType T;
  SpectraType Spectra[10];
  char   GBathFileName[255], PsiFileName[255], SpecFileName[25], Buffer[255];
  uint32_t iOuter, i, PsiNBathV, PsiNBathC;
  uint16_t oadd[1];
  double Ef = 0.0, epsilonAnd = 0.00001, mu;
  double PsiEnergy, PsiQuality, PsiHist[20], *PsiRho;
  /* read G0 */
  if (FileReadPlainTriDiagonalMatrixType(&G0, HSet.G0FileName))
  {
    printf("Fatal Error in main\n Failed to read G0 from file %s\n", HSet.G0FileName);
    return 1;
  };
  strcpy(G0.Name, "G0");
  /* read G_Bath_iStart and psi_iStart */
  if (HSet.iStart == 0)
  {
    sprintf(GBathFileName, "%s000", HSet.GRootFileName);
    sprintf(PsiFileName, "%s000", HSet.WaveFunctionRootFileName);
  }
  else if (HSet.iStart < 10)
  {
    sprintf(GBathFileName, "%s00%1i", HSet.GRootFileName, HSet.iStart);
    sprintf(PsiFileName, "%s00%1i", HSet.WaveFunctionRootFileName, HSet.iStart);
  }
  else if (HSet.iStart < 100)
  {
    sprintf(GBathFileName, "%s0%2i", HSet.GRootFileName, HSet.iStart);
    sprintf(PsiFileName, "%s0%2i", HSet.WaveFunctionRootFileName, HSet.iStart);
  }
  else
  {
    sprintf(GBathFileName, "%s%i", HSet.GRootFileName, HSet.iStart);
    sprintf(PsiFileName, "%s%i", HSet.WaveFunctionRootFileName, HSet.iStart);
  }
  if (FileReadPlainTriDiagonalMatrixType(&GBath, GBathFileName))
  {
    printf("Fatal Error in main\n Failed to read GBath from file %s !\n", GBathFileName);
    return 1;
  };
  if (file_exists(PsiFileName))
  {
    Psi.NBitsKey = HSet.NBitsKey;
    if (FileReadWaveFunctionPlainBinary(&Psi, PsiFileName))
    {
      printf("Fatal Error in main\n Failed to read Psi from file %s !\n", PsiFileName);
      return 1;
    };
    printf("\nRead the starting wave function from the file %s\n", PsiFileName);
    sscanf(Psi.Name, "%s %i %i", Buffer, &PsiNBathC, &PsiNBathV);
    printf("Wave-function \"%s\" found has %i conduction bath sites, %i valence bath sites and %i Fermions in basis\n",
           Psi.Name, PsiNBathC, PsiNBathV, Psi.NFermion);
    if ((2 * (PsiNBathC + PsiNBathV + 2)) != Psi.NFermion)
    {
      printf("The total number of Fermions is incompatible witht the bath sizes\n");
      printf("This wave-function can not be used to start from\n");
      return 1;
    }
  }
  else
  {
    printf("\nWave-function file: '%s' does not excist\n Started from new function\n", PsiFileName);
    CreateWaveFunctionDMFTSecond(&Psi, HSet.NBath, HSet.NBath, HSet.NBitsKey);
  }
  strcpy(GBath.Name, "Gbath");
  /* do the DMFT loops */
  for (iOuter = HSet.iStart; iOuter < HSet.iStart + HSet.Ni; iOuter++)
  {
    BeginTiming("DMFTLoop", Time);
    printf("\n\n***********************************************************************\n");
    printf("Start of loop %i from %i to %i\n", iOuter + 1, HSet.iStart, HSet.iStart + HSet.Ni);
    printf("***********************************************************************\n\n\n");
    /* Create splitted anderson matrix */
    TriDiagonalMatrixToSecondDoubleTriDiagonalMatrix(GBath, &T, &GAnd, Ef, epsilonAnd);
    /* Create GAnd and GAndTrunk as a copy pointing to the same memory */
    GAndTrunk.Ef = GAnd.Ef;
    GAndTrunk.Epsilon = GAnd.Epsilon;
    GAndTrunk.a0 = GAnd.a0;
    GAndTrunk.b0 = GAnd.b0;
    GAndTrunk.Dt.Ef = GAnd.Dt.Ef;
    GAndTrunk.Dt.Epsilon = GAnd.Dt.Epsilon;
    GAndTrunk.Dt.a0 = GAnd.Dt.a0;
    GAndTrunk.Dt.b0 = GAnd.Dt.b0;
    GAndTrunk.Dt.val.N = GAnd.Dt.val.N;
    GAndTrunk.Dt.val.a = GAnd.Dt.val.a;
    GAndTrunk.Dt.val.b = GAnd.Dt.val.b;
    GAndTrunk.Dt.con.N = GAnd.Dt.con.N;
    GAndTrunk.Dt.con.a = GAnd.Dt.con.a;
    GAndTrunk.Dt.con.b = GAnd.Dt.con.b;
    /* Create VPerturbation */
    for (i = 0; i < 4; i++)
    {
      InitOperator(&(VPerturbation[i]));
      VPerturbation[i].NFermion = 4 + 2 * (Min(GAnd.Dt.val.N, HSet.NBath) + Min(GAnd.Dt.con.N, HSet.NBath));
      VPerturbation[i].NBoson = 0;
      VPerturbation[i].NWord = Min(VPerturbation[i].NFermion, 1 + (VPerturbation[i].NFermion - 1) / 8);
      VPerturbation[i].QComplex = 0;
      InitOperatorFixedLength(&(VPerturbation[i]), 1, 1);
    }
    if (GAnd.Dt.val.N > HSet.NBath)
    {
      oadd[0] = CrAnMask + 3 + HSet.NBath;
      RealOperatorAddLader(&(VPerturbation[0]), oadd, 1, GAnd.Dt.val.b[HSet.NBath]);
      oadd[0] = CrAnMask + 3 + Min(GAnd.Dt.con.N, HSet.NBath) + 2 * HSet.NBath;
      RealOperatorAddLader(&(VPerturbation[1]), oadd, 1, GAnd.Dt.val.b[HSet.NBath]);
    }
    if (GAnd.Dt.con.N > HSet.NBath)
    {
      oadd[0] = 3 + Min(GAnd.Dt.val.N, HSet.NBath) + HSet.NBath;
      RealOperatorAddLader(&(VPerturbation[2]), oadd, 1, GAnd.Dt.con.b[HSet.NBath]);
      oadd[0] = 3 + 2 * Min(GAnd.Dt.val.N, HSet.NBath) + 2 * HSet.NBath;
      RealOperatorAddLader(&(VPerturbation[3]), oadd, 1, GAnd.Dt.con.b[HSet.NBath]);
    }
    /* Create gPerturbation */
    if (GAnd.Dt.val.N > HSet.NBath)
    {
      gPerturbation[0].N = GAnd.Dt.val.N - HSet.NBath;
      gPerturbation[1].N = GAnd.Dt.val.N - HSet.NBath;
      InitTriDiagonalMatrix(&(gPerturbation[0]));
      InitTriDiagonalMatrix(&(gPerturbation[1]));
      for (i = 1; i < GAnd.Dt.val.N - HSet.NBath; i++)
      {
        gPerturbation[0].a[i] = -GAnd.Dt.val.a[i + HSet.NBath];
        gPerturbation[0].b[i] = GAnd.Dt.val.b[i + HSet.NBath];
        gPerturbation[1].a[i] = -GAnd.Dt.val.a[i + HSet.NBath];
        gPerturbation[1].b[i] = GAnd.Dt.val.b[i + HSet.NBath];
      }
      gPerturbation[0].a[GAnd.Dt.val.N - HSet.NBath] = -GAnd.Dt.val.a[GAnd.Dt.val.N];
      gPerturbation[1].a[GAnd.Dt.val.N - HSet.NBath] = -GAnd.Dt.val.a[GAnd.Dt.val.N];
    }
    else
    {
      gPerturbation[0].N = 1;
      gPerturbation[1].N = 1;
      InitTriDiagonalMatrix(&(gPerturbation[0]));
      InitTriDiagonalMatrix(&(gPerturbation[1]));
      gPerturbation[0].a[1] = 0;
      gPerturbation[1].a[1] = 0;
    }
    if (GAnd.Dt.con.N > HSet.NBath)
    {
      gPerturbation[2].N = GAnd.Dt.con.N - HSet.NBath;
      gPerturbation[3].N = GAnd.Dt.con.N - HSet.NBath;
      InitTriDiagonalMatrix(&(gPerturbation[2]));
      InitTriDiagonalMatrix(&(gPerturbation[3]));
      for (i = 1; i < GAnd.Dt.con.N - HSet.NBath; i++)
      {
        gPerturbation[2].a[i] = GAnd.Dt.con.a[i + HSet.NBath];
        gPerturbation[2].b[i] = GAnd.Dt.con.b[i + HSet.NBath];
        gPerturbation[3].a[i] = GAnd.Dt.con.a[i + HSet.NBath];
        gPerturbation[3].b[i] = GAnd.Dt.con.b[i + HSet.NBath];
      }
      gPerturbation[2].a[GAnd.Dt.con.N - HSet.NBath] = GAnd.Dt.con.a[GAnd.Dt.con.N];
      gPerturbation[3].a[GAnd.Dt.con.N - HSet.NBath] = GAnd.Dt.con.a[GAnd.Dt.con.N];
    }
    else
    {
      gPerturbation[2].N = 1;
      gPerturbation[3].N = 1;
      InitTriDiagonalMatrix(&(gPerturbation[2]));
      InitTriDiagonalMatrix(&(gPerturbation[3]));
      gPerturbation[2].a[1] = 0;
      gPerturbation[3].a[1] = 0;
    };
    /* Create HAnd */
    FreeMatrix(&T);
    if (GAnd.Dt.con.N < HSet.NBath)
    {
      printf("Warning: GAnd.Dt.con.N(%i)<HSet.NBath(%i)\n", GAnd.Dt.con.N, HSet.NBath);
    }
    else
    {
      GAndTrunk.Dt.con.N = HSet.NBath;
    }
    if (GAnd.Dt.val.N < HSet.NBath)
    {
      printf("Warning: GAnd.Dt.val.N(%i)<HSet.NBath(%i)\n", GAnd.Dt.val.N, HSet.NBath);
    }
    else
    {
      GAndTrunk.Dt.val.N = HSet.NBath;
    }
    mu = -HSet.U / 2;
    SecondDoubleTriDiagonalMatrixToAndersonImpurityOperator(&HAnd, GAndTrunk, GAndTrunk, HSet.U, mu);
    /* test if psi has correct dimensions and correct if needed */
    sscanf(Psi.Name, "%s %i %i", Buffer, &PsiNBathC, &PsiNBathV);
    if (PsiNBathC != GAndTrunk.Dt.con.N)
    {
      printf("\nBefore diagonalizing:\nThe wavefunction has %i conduction bath sites\nThe Anderson Hamiltonian %i\n",
             PsiNBathC, GAndTrunk.Dt.con.N);
      if (PsiNBathC > GAndTrunk.Dt.con.N)
      {
        printf("Removing %i unoccupied fermions from the wave-function\n", PsiNBathC - GAndTrunk.Dt.con.N);
        for (i = 0; i < PsiNBathC - GAndTrunk.Dt.con.N; i++)
        {
          WaveFunctionRemoveFermion(Psi, &PsiTemp, Psi.NFermion - 1);
          WaveFunctionFree(&Psi);
          WaveFunctionRemoveFermion(PsiTemp, &Psi, 4 + PsiNBathV + PsiNBathC - i - 1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      else
      {
        printf("Adding %i unoccupied fermions to the wave-function\n", GAndTrunk.Dt.con.N - PsiNBathC);
        for (i = 0; i < GAndTrunk.Dt.con.N - PsiNBathC; i++)
        {
          WaveFunctionAddFermion(Psi, &PsiTemp, Psi.NFermion, 0);
          WaveFunctionFree(&Psi);
          WaveFunctionAddFermion(PsiTemp, &Psi, 4 + PsiNBathV + PsiNBathC + i, 0);
          WaveFunctionFree(&PsiTemp);
        }
      }
    }
    if (PsiNBathV != GAndTrunk.Dt.val.N)
    {
      printf("\nBefore diagonalizing:\nThe wavefunction has %i valence bath sites\nThe Anderson Hamiltonian %i\n",
             PsiNBathV, GAndTrunk.Dt.val.N);
      if (PsiNBathV > GAndTrunk.Dt.val.N)
      {
        printf("Removing %i occupied fermions from the wave-function\n", PsiNBathV - GAndTrunk.Dt.val.N);
        for (i = 0; i < PsiNBathV - GAndTrunk.Dt.val.N; i++)
        {
          WaveFunctionRemoveFermion(Psi, &PsiTemp, Psi.NFermion - GAndTrunk.Dt.con.N - 1);
          WaveFunctionFree(&Psi);
          WaveFunctionRemoveFermion(PsiTemp, &Psi, 4 + PsiNBathV - i - 1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      else
      {
        printf("Adding %i occupied fermions to the wave-function\n", GAndTrunk.Dt.val.N - PsiNBathV);
        for (i = 0; i < GAndTrunk.Dt.val.N - PsiNBathV; i++)
        {
          WaveFunctionAddFermion(Psi, &PsiTemp, Psi.NFermion - GAndTrunk.Dt.val.N, 1);
          WaveFunctionFree(&Psi);
          WaveFunctionAddFermion(PsiTemp, &Psi, 4 + PsiNBathV + i, 1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      printf("Corrected basis\n");
      PrintWaveFunctionReadableTruncated(Psi, 0.1);
    }
    sprintf(Psi.Name, "DMFT %i %i", GAndTrunk.Dt.con.N, GAndTrunk.Dt.val.N);
    /* create ground-state wave-function */
    BeginTiming("Lanczos GroundState", Time);
    LanczosGroundState(HAnd, &Psi, &HSet.zeroLancWavefunction, HSet.zeroLancEnergy, HSet.U, &Psi, 0);
    LanczosGroundState(HAnd, &Psi, &HSet.zeroLancWavefunction, HSet.zeroLancEnergy, HSet.U, &Psi, 0);
    EndTiming("Lanczos GroundState", Time);
    PsiHermitianOperatorPsi(HAnd, Psi, &PsiEnergy);
    PsiHermitianOperatorSqrPsi(HAnd, Psi, &PsiQuality);
    PsiQuality /= Sqr(PsiEnergy);
    PrintWaveFunctionReadableTruncated(Psi, 0.1);
    CalculatePrefactorHistogram(PsiHist, Psi, 10.0, 20);
    printf("\nHistogram of prefactors in Wave-function\n");
    for (i = 0; i < 19; i++)
    {
      printf("10^(-%2i) to 10^(-%2i) | %22.15E\n", i, i + 1, PsiHist[i]);
    }
    printf("10^(-%2i) to 0        | %22.15E\n", i, PsiHist[19]);
    printf("\n");
    printf("Energy                  <psi|H^2|psi>/(<psi|H|psi>^2)\n");
    printf("%22.15E  %22.15E\n", PsiEnergy, PsiQuality);
    PsiRho = calloc(Psi.NFermion, sizeof(double));
    if (PsiRho == NULL)
    {
      printf("calloc failed in main for PsiRho\n Fatal Error\n");
      fflush(stdout);
      return 1;
    }
    CalculateDensity(PsiRho, Psi);
    printf("Occupation \n");
    for (i = 0; i < Psi.NFermion; i++)
    {
      printf("%5i %22.15E\n", i, PsiRho[i]);
    }
    free(PsiRho);
    /*                                           */
    /* create spin up and spin down PES and IPES */
    /*                                           */
    CreatePESdnOperatorDMFT(&OPESdn, HAnd.NFermion);
    CreatePESupOperatorDMFT(&OPESup, HAnd.NFermion);
    CreateIPESdnOperatorDMFT(&OIPESdn, HAnd.NFermion);
    CreateIPESupOperatorDMFT(&OIPESup, HAnd.NFermion);
    GclMup.N = Max(HSet.NTriDiag, G0.N);
    InitTriDiagonalMatrix(&GclMup);
    GclMdn.N = Max(HSet.NTriDiag, G0.N);
    InitTriDiagonalMatrix(&GclMdn);
    GclPup.N = Max(HSet.NTriDiag, G0.N);
    InitTriDiagonalMatrix(&GclPup);
    GclPdn.N = Max(HSet.NTriDiag, G0.N);
    InitTriDiagonalMatrix(&GclPdn);
    BeginTiming("Lanczos Tri Diagonalization Dyson", Time);
    OperatorPsi(OPESup, Psi, &OPsi);
    LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation, gPerturbation,
                                     4, &GclMup);
    EndTiming("Lanczos Tri Diagonalization Dyson", Time);
    BeginTiming("Lanczos Tri Diagonalization Dyson", Time);
    OperatorPsi(OPESdn, Psi, &OPsi);
    LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation, gPerturbation,
                                     4, &GclMdn);
    EndTiming("Lanczos Tri Diagonalization Dyson", Time);
    BeginTiming("Lanczos Tri Diagonalization Dyson", Time);
    OperatorPsi(OIPESup, Psi, &OPsi);
    LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation, gPerturbation,
                                     4, &GclPup);
    EndTiming("Lanczos Tri Diagonalization Dyson", Time);
    BeginTiming("Lanczos Tri Diagonalization Dyson", Time);
    OperatorPsi(OIPESdn, Psi, &OPsi);
    LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation, gPerturbation,
                                     4, &GclPdn);
    EndTiming("Lanczos Tri Diagonalization Dyson", Time);
    /* G math */
    printf("Start G math\n");
    BeginTiming("G math", Time);
    ShiftTriDiagonalMatrix(&GclMdn, -PsiEnergy);
    ShiftTriDiagonalMatrix(&GclMup, -PsiEnergy);
    ShiftTriDiagonalMatrix(&GclPdn, -PsiEnergy);
    ShiftTriDiagonalMatrix(&GclPup, -PsiEnergy);
    for (i = 1; i <= GclMdn.N; i++)
    {
      GclMdn.a[i] *= -1;
    }
    for (i = 1; i <= GclMup.N; i++)
    {
      GclMup.a[i] *= -1;
    }
    strcpy(GclPup.Name, "Gcluster-Plus-Up");
    strcpy(GclPdn.Name, "Gcluster-Plus-Dn");
    strcpy(GclMup.Name, "Gcluster-Min-Up");
    strcpy(GclMdn.Name, "Gcluster-Min-Dn");
    printf("  Start SumGreensfunctionsTriDiagonalRepresentation\n");
    fflush(stdout);
    SumGreensfunctionsTriDiagonalRepresentation(GclPup, GclMup, &Gclup);
    strcpy(Gclup.Name, "Gcluster-Up");
    SumGreensfunctionsTriDiagonalRepresentation(GclPdn, GclMdn, &Gcldn);
    strcpy(Gcldn.Name, "Gcluster-Dn");
    SumGreensfunctionsTriDiagonalRepresentation(Gclup, Gcldn, &Gcl);
    strcpy(Gcl.Name, "Gcluster");
    Gcl.b[0] = 1.0;
    printf("Gcl.N = %i ", Gcl.N);
    // SymmetrizeTridiagonalMatrix(&Gcl);
    Gcl.N = Min(Gcl.N, Max(4 * HSet.NTriDiag, 2 * G0.N));
    printf("=> %i\n", Gcl.N);
    BeginTiming("G math Selfenergy", Time);
    printf("  Start SelfEnergyTriDiagonalRepresentation\n");
    printf("GBath.N = %i ", GBath.N);
    SelfEnergyTriDiagonalRepresentation(GBath, Gcl, &Sigma);
    printf("Sigma.N = %i ", Sigma.N);
    Sigma.N = Min(Sigma.N, Max(4 * HSet.NTriDiag, 2 * G0.N));
    printf("=> %i ", Sigma.N);
    strcpy(Sigma.Name, "Sigma");
    // SymmetrizeTridiagonalMatrix(&Sigma);
    printf("=> %i ", Sigma.N);
    Sigma.N = Min(Sigma.N, Max(4 * HSet.NTriDiag, 2 * G0.N));
    printf("=> %i\n", Sigma.N);
    EndTiming("G math Selfenergy", Time);
    printf("  Start GLocalTriDiagonalRepresentation\n");
    GLocalTriDiagonalRepresentation(G0, Sigma, &GLocal, Sigma.N + G0.N - 1);
    printf("Glocal.N = %i ", GLocal.N);
    GLocal.N = Min(GLocal.N, Max(4 * HSet.NTriDiag, 2 * G0.N));
    printf("=> %i ", GLocal.N);
    strcpy(GLocal.Name, "Glocal");
    // SymmetrizeTridiagonalMatrix(&GLocal);
    printf("=> %i ", GLocal.N);
    GLocal.N = Min(GLocal.N, Max(4 * HSet.NTriDiag, 2 * G0.N));
    printf("=> %i\n", GLocal.N);
    printf("  Start GBathNewTriDiagonalRepresentation\n");
    GBathNewTriDiagonalRepresentation(GLocal, Sigma, &GBathNew);
    GBathNew.N = Min(GBathNew.N, Max(4 * HSet.NTriDiag, 2 * G0.N));
    strcpy(GBathNew.Name, "Gbath-New");
    // SymmetrizeTridiagonalMatrix(&GBathNew);
    GBathNew.N = Min(GBathNew.N, Max(4 * HSet.NTriDiag, 2 * G0.N));
    /* calculate GBath for the next loop by mixing GBath with GBathNew */
    printf("Start mixing GBath with GBathNew\n");
    TriDiagonalMatrixCopy(&GBathOld, GBath);
    strcpy(GBathOld.Name, "Gbath-Old");
    FreeTriDiagonalMatrix(&GBath);
    GBathOld.b[0] = sqrt(HSet.Mix);
    GBathNew.b[0] = sqrt(1 - HSet.Mix);
    SumGreensfunctionsTriDiagonalRepresentation(GBathNew, GBathOld, &GBath);
    strcpy(GBath.Name, "Gbath");
    GBathOld.b[0] = 1;
    GBathNew.b[0] = 1;
    // SymmetrizeTridiagonalMatrix(&GBath);
    GBath.N = Min(GBath.N, Max(4 * HSet.NTriDiag, 2 * G0.N));
    EndTiming("G math", Time);
    /* create spectra */
    printf("Start create spectra\n");
    Spectra[0].N = HSet.NE;
    Spectra[0].Emin = HSet.Emin;
    Spectra[0].Emax = HSet.Emax;
    Spectra[0].dE = (HSet.Emax - HSet.Emin) / HSet.NE;
    Spectra[0].EminPole = HSet.Emin;
    Spectra[0].EmaxPole = HSet.Emax;
    Spectra[0].Gamma = HSet.Gamma;
    for (i = 0; i < 10; i++)
    {
      SpectrumAllocate(Spectra[0], &Spectra[i]);
    }
    TriDiagonalMatrixToSpectrum(0, GBath, &(Spectra[0])); /*  3 */
    TriDiagonalMatrixToSpectrum(0, GBathOld, &(Spectra[1]));  /*  5 */
    TriDiagonalMatrixToSpectrum(0, GBathNew, &(Spectra[2]));  /*  7 */
    SecondDoubleTriDiagonalMatrixToSpectrum(0, GAnd, &(Spectra[3]));  /*  9 */
    TriDiagonalMatrixToSpectrum(0, G0, &(Spectra[4]));  /* 11 */
    TriDiagonalMatrixToSpectrum(0, Sigma, &(Spectra[5])); /* 13 */
    TriDiagonalMatrixToSpectrum(0, GLocal, &(Spectra[6]));  /* 15 */
    TriDiagonalMatrixToSpectrum(0, GclPup, &(Spectra[7]));  /* 17 */
    TriDiagonalMatrixToSpectrum(0, GclMup, &(Spectra[8]));  /* 19 */
    TriDiagonalMatrixToSpectrum(0, Gcl, &(Spectra[9])); /* 21 */
    /* write to disk */
    printf("Start write to dist\n");
    if (iOuter == 0)
    {
      sprintf(GBathFileName, "%s001", HSet.GRootFileName);
      sprintf(PsiFileName, "%s001", HSet.WaveFunctionRootFileName);
      sprintf(SpecFileName, "%s001", HSet.SpecRootFileName);
    }
    else if (iOuter < 9)
    {
      sprintf(GBathFileName, "%s00%1i", HSet.GRootFileName, iOuter + 1);
      sprintf(PsiFileName, "%s00%1i", HSet.WaveFunctionRootFileName, iOuter + 1);
      sprintf(SpecFileName, "%s00%1i", HSet.SpecRootFileName, iOuter + 1);
    }
    else if (iOuter < 99)
    {
      sprintf(GBathFileName, "%s0%2i", HSet.GRootFileName, iOuter + 1);
      sprintf(PsiFileName, "%s0%2i", HSet.WaveFunctionRootFileName, iOuter + 1);
      sprintf(SpecFileName, "%s0%2i", HSet.SpecRootFileName, iOuter + 1);
    }
    else if (iOuter < 999)
    {
      sprintf(GBathFileName, "%s%3i", HSet.GRootFileName, iOuter + 1);
      sprintf(PsiFileName, "%s%3i", HSet.WaveFunctionRootFileName, iOuter + 1);
      sprintf(SpecFileName, "%s%3i", HSet.SpecRootFileName, iOuter + 1);
    }
    FilePrintSpectra(Spectra, 10, SpecFileName, "w");
    FilePrintPlainTriDiagonalMatrix(GBath, GBathFileName, "w");
    FilePrintWaveFunctionPlainBinary(Psi, PsiFileName, "w");
    /* free variables */
    printf("Start free variables\n");
    for (i = 0; i < 4; i++)
    {
      FreeOperator(&(VPerturbation[i]));
      FreeTriDiagonalMatrix(&(gPerturbation[i]));
    }
    FreeOperator(&HAnd);
    FreeOperator(&OPESdn);
    FreeOperator(&OPESup);
    FreeOperator(&OIPESdn);
    FreeOperator(&OIPESup);
    FreeTriDiagonalMatrix(&GclMup);
    FreeTriDiagonalMatrix(&GclMdn);
    FreeTriDiagonalMatrix(&GclPup);
    FreeTriDiagonalMatrix(&GclPdn);
    FreeTriDiagonalMatrix(&Gcl);
    FreeTriDiagonalMatrix(&GLocal);
    FreeTriDiagonalMatrix(&Sigma);
    FreeTriDiagonalMatrix(&GBathNew);
    FreeTriDiagonalMatrix(&GBathOld);
    FreeDoubleTriDiagonalMatrix(&(GAnd.Dt));
    EndTiming("DMFTLoop", Time);

    PrintTiming(*Time);
  }
  return 0;
};


int oldmainWorkingOnDoubleTriDiagonalBathFormat(int argc, const char *argv[])
{
  TimeType Time;
  DMFTSettings HSet;
  InitializeTiming(&Time);
  BeginTiming("Total", &Time);
  TimeStampAll();
  argc -= 1;
  /* read the DMFT settings */
  if (FileReadDMFTSettings(&HSet, argv[argc]))
  {
    printf("Fatal Error in main\n Failed to read DMFT Settings from file %s\n", argv[argc]);
    return 1;
  };
  PrintDMFTSettings(HSet);
  WaveFunctionType Psi, OPsi, PsiTemp;
  TriDiagonalMatrixType G0, GBath, GclPup, GclPdn, GclMup, GclMdn, Gcl, GBathNew, GBathOld, gPerturbation[4], Sigma;
  SecondDoubleTriDiagonalMatrixType GAnd, GAndTrunk;
  AndersonMatrixType GBathAnd, GclAnd, G0And, GBathNewAnd;
  ListOfPoles SigmaPoles, GclList, GListTemp;
  OperatorType HAnd, OPESdn, OPESup, OIPESdn, OIPESup, VPerturbation[4];
  MatrixType T;
  SpectraType Spectra[10];
  char   GBathFileName[255], GClFileName[255], GSigmaFileName[255], PsiFileName[255], SpecFileName[255], Buffer[255];
  uint32_t iOuter, i, PsiNBathV, PsiNBathC;
  uint16_t oadd[1];
  double Ef = 0.0, epsilonAnd = 0.00001, mu;
  double PsiEnergy, PsiQuality, PsiHist[20], *PsiRho;
  /* read G0 */
  if (FileReadPlainTriDiagonalMatrixType(&G0, HSet.G0FileName))
  {
    printf("Fatal Error in main\n Failed to read G0 from file %s\n", HSet.G0FileName);
    return 1;
  };
  strcpy(G0.Name, "G0");
  /* read G_Bath_iStart and psi_iStart */
  if (HSet.iStart == 0)
  {
    sprintf(GBathFileName, "%s000", HSet.GRootFileName);
    sprintf(PsiFileName, "%s000", HSet.WaveFunctionRootFileName);
  }
  else if (HSet.iStart < 10)
  {
    sprintf(GBathFileName, "%s00%1i", HSet.GRootFileName, HSet.iStart);
    sprintf(PsiFileName, "%s00%1i", HSet.WaveFunctionRootFileName, HSet.iStart);
  }
  else if (HSet.iStart < 100)
  {
    sprintf(GBathFileName, "%s0%2i", HSet.GRootFileName, HSet.iStart);
    sprintf(PsiFileName, "%s0%2i", HSet.WaveFunctionRootFileName, HSet.iStart);
  }
  else
  {
    sprintf(GBathFileName, "%s%i", HSet.GRootFileName, HSet.iStart);
    sprintf(PsiFileName, "%s%i", HSet.WaveFunctionRootFileName, HSet.iStart);
  }
  if (FileReadPlainTriDiagonalMatrixType(&GBath, GBathFileName))
  {
    printf("Fatal Error in main\n Failed to read GBath from file %s !\n", GBathFileName);
    return 1;
  };
  if (file_exists(PsiFileName))
  {
    Psi.NBitsKey = HSet.NBitsKey;
    if (FileReadWaveFunctionPlainBinary(&Psi, PsiFileName))
    {
      printf("Fatal Error in main\n Failed to read Psi from file %s !\n", PsiFileName);
      return 1;
    };
    printf("\nRead the starting wave function from the file %s\n", PsiFileName);
    sscanf(Psi.Name, "%s %i %i", Buffer, &PsiNBathC, &PsiNBathV);
    printf("Wave-function \"%s\" found has %i conduction bath sites, %i valence bath sites and %i Fermions in basis\n",
           Psi.Name, PsiNBathC, PsiNBathV, Psi.NFermion);
    if ((2 * (PsiNBathC + PsiNBathV + 2)) != Psi.NFermion)
    {
      printf("The total number of Fermions is incompatible witht the bath sizes\n");
      printf("This wave-function can not be used to start from\n");
      return 1;
    }
  }
  else
  {
    printf("\nWave-function file: '%s' does not excist\n Started from new function\n", PsiFileName);
    CreateWaveFunctionDMFTSecond(&Psi, HSet.NBath, HSet.NBath, HSet.NBitsKey);
  }
  strcpy(GBath.Name, "Gbath");
  TriDiagonalMatrixToAndersonMatrix(G0, &T, &G0And);
  FreeMatrix(&T);
  /* do the DMFT loops */
  for (iOuter = HSet.iStart; iOuter < HSet.iStart + HSet.Ni; iOuter++)
  {
    BeginTiming("DMFTLoop", &Time);
    printf("\n\n***********************************************************************\n");
    printf("Start of loop %i from %i to %i\n", iOuter + 1, HSet.iStart, HSet.iStart + HSet.Ni);
    printf("***********************************************************************\n\n\n");
    /* Create splitted anderson matrix */
    TriDiagonalMatrixToSecondDoubleTriDiagonalMatrix(GBath, &T, &GAnd, Ef, epsilonAnd);
    /* Create GAnd and GAndTrunk as a copy pointing to the same memory */
    GAndTrunk.Ef = GAnd.Ef;
    GAndTrunk.Epsilon = GAnd.Epsilon;
    GAndTrunk.a0 = GAnd.a0;
    GAndTrunk.b0 = GAnd.b0;
    GAndTrunk.Dt.Ef = GAnd.Dt.Ef;
    GAndTrunk.Dt.Epsilon = GAnd.Dt.Epsilon;
    GAndTrunk.Dt.a0 = GAnd.Dt.a0;
    GAndTrunk.Dt.b0 = GAnd.Dt.b0;
    GAndTrunk.Dt.val.N = GAnd.Dt.val.N;
    GAndTrunk.Dt.val.a = GAnd.Dt.val.a;
    GAndTrunk.Dt.val.b = GAnd.Dt.val.b;
    GAndTrunk.Dt.con.N = GAnd.Dt.con.N;
    GAndTrunk.Dt.con.a = GAnd.Dt.con.a;
    GAndTrunk.Dt.con.b = GAnd.Dt.con.b;
    if (HSet.Dyson)
    {
      /* Create VPerturbation */
      for (i = 0; i < 4; i++)
      {
        InitOperator(&(VPerturbation[i]));
        VPerturbation[i].NFermion = 4 + 2 * (Min(GAnd.Dt.val.N, HSet.NBath) + Min(GAnd.Dt.con.N, HSet.NBath));
        VPerturbation[i].NBoson = 0;
        VPerturbation[i].NWord = Min(VPerturbation[i].NFermion, 1 + (VPerturbation[i].NFermion - 1) / 8);
        VPerturbation[i].QComplex = 0;
        InitOperatorFixedLength(&(VPerturbation[i]), 1, 1);
      }
      if (GAnd.Dt.val.N > HSet.NBath)
      {
        oadd[0] = CrAnMask + 3 + HSet.NBath;
        RealOperatorAddLader(&(VPerturbation[0]), oadd, 1, GAnd.Dt.val.b[HSet.NBath]);
        oadd[0] = CrAnMask + 3 + Min(GAnd.Dt.con.N, HSet.NBath) + 2 * HSet.NBath;
        RealOperatorAddLader(&(VPerturbation[1]), oadd, 1, GAnd.Dt.val.b[HSet.NBath]);
      }
      if (GAnd.Dt.con.N > HSet.NBath)
      {
        oadd[0] = 3 + Min(GAnd.Dt.val.N, HSet.NBath) + HSet.NBath;
        RealOperatorAddLader(&(VPerturbation[2]), oadd, 1, GAnd.Dt.con.b[HSet.NBath]);
        oadd[0] = 3 + 2 * Min(GAnd.Dt.val.N, HSet.NBath) + 2 * HSet.NBath;
        RealOperatorAddLader(&(VPerturbation[3]), oadd, 1, GAnd.Dt.con.b[HSet.NBath]);
      }
      /* Create gPerturbation */
      if (GAnd.Dt.val.N > HSet.NBath)
      {
        gPerturbation[0].N = GAnd.Dt.val.N - HSet.NBath;
        gPerturbation[1].N = GAnd.Dt.val.N - HSet.NBath;
        InitTriDiagonalMatrix(&(gPerturbation[0]));
        InitTriDiagonalMatrix(&(gPerturbation[1]));
        for (i = 1; i < GAnd.Dt.val.N - HSet.NBath; i++)
        {
          gPerturbation[0].a[i] = -GAnd.Dt.val.a[i + HSet.NBath];
          gPerturbation[0].b[i] = GAnd.Dt.val.b[i + HSet.NBath];
          gPerturbation[1].a[i] = -GAnd.Dt.val.a[i + HSet.NBath];
          gPerturbation[1].b[i] = GAnd.Dt.val.b[i + HSet.NBath];
        }
        gPerturbation[0].a[GAnd.Dt.val.N - HSet.NBath] = -GAnd.Dt.val.a[GAnd.Dt.val.N];
        gPerturbation[1].a[GAnd.Dt.val.N - HSet.NBath] = -GAnd.Dt.val.a[GAnd.Dt.val.N];
      }
      else
      {
        gPerturbation[0].N = 1;
        gPerturbation[1].N = 1;
        InitTriDiagonalMatrix(&(gPerturbation[0]));
        InitTriDiagonalMatrix(&(gPerturbation[1]));
        gPerturbation[0].a[1] = 0;
        gPerturbation[1].a[1] = 0;
      }
      if (GAnd.Dt.con.N > HSet.NBath)
      {
        gPerturbation[2].N = GAnd.Dt.con.N - HSet.NBath;
        gPerturbation[3].N = GAnd.Dt.con.N - HSet.NBath;
        InitTriDiagonalMatrix(&(gPerturbation[2]));
        InitTriDiagonalMatrix(&(gPerturbation[3]));
        for (i = 1; i < GAnd.Dt.con.N - HSet.NBath; i++)
        {
          gPerturbation[2].a[i] = GAnd.Dt.con.a[i + HSet.NBath];
          gPerturbation[2].b[i] = GAnd.Dt.con.b[i + HSet.NBath];
          gPerturbation[3].a[i] = GAnd.Dt.con.a[i + HSet.NBath];
          gPerturbation[3].b[i] = GAnd.Dt.con.b[i + HSet.NBath];
        }
        gPerturbation[2].a[GAnd.Dt.con.N - HSet.NBath] = GAnd.Dt.con.a[GAnd.Dt.con.N];
        gPerturbation[3].a[GAnd.Dt.con.N - HSet.NBath] = GAnd.Dt.con.a[GAnd.Dt.con.N];
      }
      else
      {
        gPerturbation[2].N = 1;
        gPerturbation[3].N = 1;
        InitTriDiagonalMatrix(&(gPerturbation[2]));
        InitTriDiagonalMatrix(&(gPerturbation[3]));
        gPerturbation[2].a[1] = 0;
        gPerturbation[3].a[1] = 0;
      }
    }
    /* Create HAnd */
    FreeMatrix(&T);
    if (GAnd.Dt.con.N < HSet.NBath)
    {
      printf("Warning: GAnd.Dt.con.N(%i)<HSet.NBath(%i)\n", GAnd.Dt.con.N, HSet.NBath);
    }
    else
    {
      GAndTrunk.Dt.con.N = HSet.NBath;
    }
    if (GAnd.Dt.val.N < HSet.NBath)
    {
      printf("Warning: GAnd.Dt.val.N(%i)<HSet.NBath(%i)\n", GAnd.Dt.val.N, HSet.NBath);
    }
    else
    {
      GAndTrunk.Dt.val.N = HSet.NBath;
    }
    mu = -HSet.U / 2;
    SecondDoubleTriDiagonalMatrixToAndersonImpurityOperator(&HAnd, GAndTrunk, GAndTrunk, HSet.U, mu);
    /* test if psi has correct dimensions and correct if needed */
    sscanf(Psi.Name, "%s %i %i", Buffer, &PsiNBathC, &PsiNBathV);
    if (PsiNBathC != GAndTrunk.Dt.con.N)
    {
      printf("\nBefore diagonalizing:\nThe wavefunction has %i conduction bath sites\nThe Anderson Hamiltonian %i\n",
             PsiNBathC, GAndTrunk.Dt.con.N);
      if (PsiNBathC > GAndTrunk.Dt.con.N)
      {
        printf("Removing %i unoccupied fermions from the wave-function\n", PsiNBathC - GAndTrunk.Dt.con.N);
        for (i = 0; i < PsiNBathC - GAndTrunk.Dt.con.N; i++)
        {
          WaveFunctionRemoveFermion(Psi, &PsiTemp, Psi.NFermion - 1);
          WaveFunctionFree(&Psi);
          WaveFunctionRemoveFermion(PsiTemp, &Psi, 4 + PsiNBathV + PsiNBathC - i - 1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      else
      {
        printf("Adding %i unoccupied fermions to the wave-function\n", GAndTrunk.Dt.con.N - PsiNBathC);
        for (i = 0; i < GAndTrunk.Dt.con.N - PsiNBathC; i++)
        {
          WaveFunctionAddFermion(Psi, &PsiTemp, Psi.NFermion, 0);
          WaveFunctionFree(&Psi);
          WaveFunctionAddFermion(PsiTemp, &Psi, 4 + PsiNBathV + PsiNBathC + i, 0);
          WaveFunctionFree(&PsiTemp);
        }
      }
    }
    if (PsiNBathV != GAndTrunk.Dt.val.N)
    {
      printf("\nBefore diagonalizing:\nThe wavefunction has %i valence bath sites\nThe Anderson Hamiltonian %i\n",
             PsiNBathV, GAndTrunk.Dt.val.N);
      if (PsiNBathV > GAndTrunk.Dt.val.N)
      {
        printf("Removing %i occupied fermions from the wave-function\n", PsiNBathV - GAndTrunk.Dt.val.N);
        for (i = 0; i < PsiNBathV - GAndTrunk.Dt.val.N; i++)
        {
          WaveFunctionRemoveFermion(Psi, &PsiTemp, Psi.NFermion - GAndTrunk.Dt.con.N - 1);
          WaveFunctionFree(&Psi);
          WaveFunctionRemoveFermion(PsiTemp, &Psi, 4 + PsiNBathV - i - 1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      else
      {
        printf("Adding %i occupied fermions to the wave-function\n", GAndTrunk.Dt.val.N - PsiNBathV);
        for (i = 0; i < GAndTrunk.Dt.val.N - PsiNBathV; i++)
        {
          WaveFunctionAddFermion(Psi, &PsiTemp, Psi.NFermion - GAndTrunk.Dt.val.N, 1);
          WaveFunctionFree(&Psi);
          WaveFunctionAddFermion(PsiTemp, &Psi, 4 + PsiNBathV + i, 1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      printf("Corrected basis\n");
      PrintWaveFunctionReadableTruncated(Psi, 0.1);
    }
    sprintf(Psi.Name, "DMFT %i %i", GAndTrunk.Dt.con.N, GAndTrunk.Dt.val.N);
    /* create ground-state wave-function */
    BeginTiming("Lanczos GroundState", &Time);
    LanczosGroundState(HAnd, &Psi, &HSet.zeroLancWavefunction, HSet.zeroLancEnergy, HSet.U, &Psi, 0);
    LanczosGroundState(HAnd, &Psi, &HSet.zeroLancWavefunction, HSet.zeroLancEnergy, HSet.U, &Psi, 0);
    EndTiming("Lanczos GroundState", &Time);
    PsiHermitianOperatorPsi(HAnd, Psi, &PsiEnergy);
    PsiHermitianOperatorSqrPsi(HAnd, Psi, &PsiQuality);
    PsiQuality /= Sqr(PsiEnergy);
    PrintWaveFunctionReadableTruncated(Psi, 0.1);
    CalculatePrefactorHistogram(PsiHist, Psi, 10.0, 20);
    printf("\nHistogram of prefactors in Wave-function\n");
    for (i = 0; i < 19; i++)
    {
      printf("10^(-%2i) to 10^(-%2i) | %22.15E\n", i, i + 1, PsiHist[i]);
    }
    printf("10^(-%2i) to 0        | %22.15E\n", i, PsiHist[19]);
    printf("\n");
    printf("Energy                  <psi|H^2|psi>/(<psi|H|psi>^2)\n");
    printf("%22.15E  %22.15E\n", PsiEnergy, PsiQuality);
    PsiRho = calloc(Psi.NFermion, sizeof(double));
    if (PsiRho == NULL)
    {
      printf("calloc failed in main for PsiRho\n Fatal Error\n");
      fflush(stdout);
      return 1;
    }
    CalculateDensity(PsiRho, Psi);
    printf("Occupation \n");
    for (i = 0; i < Psi.NFermion; i++)
    {
      printf("%5i %22.15E\n", i, PsiRho[i]);
    }
    free(PsiRho);
    /*                                           */
    /* create spin up and spin down PES and IPES */
    /*                                           */
    CreatePESdnOperatorDMFT(&OPESdn, HAnd.NFermion);
    CreatePESupOperatorDMFT(&OPESup, HAnd.NFermion);
    CreateIPESdnOperatorDMFT(&OIPESdn, HAnd.NFermion);
    CreateIPESupOperatorDMFT(&OIPESup, HAnd.NFermion);
    /* in the last extra loop call DMFTSpectra */
    if ((iOuter == HSet.iStart + HSet.Ni - 1) && (HSet.Spectra))
    {
      if (DMFTSpec(HSet, PsiEnergy, Psi, HAnd, VPerturbation, gPerturbation, OPESdn, OPESup, OIPESdn, OIPESup))
      {
        printf("DMFTSpec failed\n");
        return 1;
      }
    }
    /* if not the last loop continue */
    BeginTiming("Lanczos Tri Diagonalization", &Time);
    OperatorPsi(OPESup, Psi, &OPsi);
    switch (HSet.Dyson)
    {
    case 1:
      GclMup.N = Max(HSet.NTriDiag, HSet.NGCluster);
      InitTriDiagonalMatrix(&GclMup);
      LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation,
                                       gPerturbation, 4, &GclMup);
      break;
    default:
      GclMup.N = HSet.NTriDiag;
      InitTriDiagonalMatrix(&GclMup);
      LanczosTriDiagonalize(HAnd, OPsi, HSet.zeroLancWavefunction, &GclMup, HSet.Gamma);
      break;
    }
    EndTiming("Lanczos Tri Diagonalization", &Time);
    BeginTiming("Lanczos Tri Diagonalization", &Time);
    OperatorPsi(OPESdn, Psi, &OPsi);
    switch (HSet.Dyson)
    {
    case 1:
      GclMdn.N = Max(HSet.NTriDiag, HSet.NGCluster);
      InitTriDiagonalMatrix(&GclMdn);
      LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation,
                                       gPerturbation, 4, &GclMdn);
      break;
    default:
      GclMdn.N = HSet.NTriDiag;
      InitTriDiagonalMatrix(&GclMdn);
      LanczosTriDiagonalize(HAnd, OPsi, HSet.zeroLancWavefunction, &GclMdn, HSet.Gamma);
      break;
    }
    EndTiming("Lanczos Tri Diagonalization", &Time);
    BeginTiming("Lanczos Tri Diagonalization", &Time);
    OperatorPsi(OIPESup, Psi, &OPsi);
    switch (HSet.Dyson)
    {
    case 1:
      GclPup.N = Max(HSet.NTriDiag, HSet.NGCluster);
      InitTriDiagonalMatrix(&GclPup);
      LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation,
                                       gPerturbation, 4, &GclPup);
      break;
    default:
      GclPup.N = HSet.NTriDiag;
      InitTriDiagonalMatrix(&GclPup);
      LanczosTriDiagonalize(HAnd, OPsi, HSet.zeroLancWavefunction, &GclPup, HSet.Gamma);
      break;
    }
    EndTiming("Lanczos Tri Diagonalization", &Time);
    BeginTiming("Lanczos Tri Diagonalization", &Time);
    OperatorPsi(OIPESdn, Psi, &OPsi);
    switch (HSet.Dyson)
    {
    case 1:
      GclPdn.N = Max(HSet.NTriDiag, HSet.NGCluster);
      InitTriDiagonalMatrix(&GclPdn);
      LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation,
                                       gPerturbation, 4, &GclPdn);
      break;
    default:
      GclPdn.N = HSet.NTriDiag;
      InitTriDiagonalMatrix(&GclPdn);
      LanczosTriDiagonalize(HAnd, OPsi, HSet.zeroLancWavefunction, &GclPdn, HSet.Gamma);
      break;
    }
    EndTiming("Lanczos Tri Diagonalization", &Time);
    /* G math - create Gcl as the sum of PES - IPES */
    printf("=== Start G math ===\n");
    fflush(stdout);
    BeginTiming("G math", &Time);
    ShiftTriDiagonalMatrix(&GclMdn, -PsiEnergy);
    ShiftTriDiagonalMatrix(&GclMup, -PsiEnergy);
    ShiftTriDiagonalMatrix(&GclPdn, -PsiEnergy);
    ShiftTriDiagonalMatrix(&GclPup, -PsiEnergy);
    for (i = 1; i <= GclMdn.N; i++)
    {
      GclMdn.a[i] *= -1;
    }
    for (i = 1; i <= GclMup.N; i++)
    {
      GclMup.a[i] *= -1;
    }
    strcpy(GclPup.Name, "Gcluster-Plus-Up");
    strcpy(GclPdn.Name, "Gcluster-Plus-Dn");
    strcpy(GclMup.Name, "Gcluster-Min-Up");
    strcpy(GclMdn.Name, "Gcluster-Min-Dn");
    printf("Start SumGreensfunctionsTriDiagonalRepresentation\n");
    fflush(stdout);
    TriDiagonalMatrixToListOfPoles(GclPup, &GclList);
    TriDiagonalMatrixToListOfPoles(GclPdn, &GListTemp);
    AddListOfPoles(&GclList, GListTemp);
    FreeListOfPoles(&GListTemp);
    TriDiagonalMatrixToListOfPoles(GclMup, &GListTemp);
    AddListOfPoles(&GclList, GListTemp);
    FreeListOfPoles(&GListTemp);
    TriDiagonalMatrixToListOfPoles(GclMdn, &GListTemp);
    AddListOfPoles(&GclList, GListTemp);
    FreeListOfPoles(&GListTemp);
    ScaleListOfPoles(&GclList, 0.5);
    /* symmetrize */
    SymmetrizeListOfPoles(&GclList);
    printf("Cluster Green function created of size %i ", GclList.N);
    ReduceListOfPolesNumber(&GclList, HSet.NGCluster);
    printf("reduced to %i\n", GclList.N);
    ListOfPolesToTridiagonalMatrix(GclList, &Gcl, GclList.N);
    strcpy(Gcl.Name, "Gcluster");
    Gcl.b[0] = 1.0;
    /* make Gcl a physical solution of the full lattice model */
    /* Create self energy (SDMFT) of G0 and Gcl SDMFT = - Gcl^-1 + G0^-1 */

    /* Create Gclprime from SDMFT Gcl = (G0^-1 - SDMFT)^-1 */

    /* G math - create Cluster Self energy From Gcl and GBath */
    BeginTiming("G math Selfenergy", &Time);
    printf("Start SelfEnergyTriDiagonalRepresentation\n");
    fflush(stdout);
    TriDiagonalMatrixToAndersonMatrix(Gcl, &T, &GclAnd);
    FreeMatrix(&T);
    TriDiagonalMatrixToAndersonMatrix(GBath, &T, &GBathAnd);
    FreeMatrix(&T);
    SelfEnergyListOfPolesFromAndersonRepresentation(GBathAnd, GclAnd, &SigmaPoles);

    /*
       GaussianBroadenListOfPoles(&SigmaPoles, 0.1, 10);
       ReduceListOfPolesNumber(&SigmaPoles, SigmaPoles.N/5);
     */

    printf("Cluster self energy created of size %i\n", SigmaPoles.N);
    strcpy(SigmaPoles.Name, "Sigma");
    /* copy sigma to tridiagonal format */
    ListOfPolesToTridiagonalMatrix(SigmaPoles, &Sigma, SigmaPoles.N);
    EndTiming("G math Selfenergy", &Time);
    /* G math - create GBath */
    BeginTiming("G math GBath", &Time);
    printf("Start GBathAndListAndRepresentation\n");
    fflush(stdout);
    GBathAndTriAndRepresentation(G0And, Sigma, &GBathNewAnd);
    printf("New bath Green function created of size %i ", GBathNewAnd.N);
    strcpy(GBathNewAnd.Name, "Gbath-New");
    ReducePolesAndersonMatrixNumber(&GBathNewAnd, HSet.NGBath);
    printf("reduced to %i\n", GBathNewAnd.N);
    AndersonMatrixToTriDiagonalMatrix(GBathNewAnd, &T, &GBathNew, GBathNewAnd.N);
    FreeMatrix(&T);
    EndTiming("G math GBath", &Time);
    /* calculate GBath for the next loop by mixing GBath with GBathNew */
    printf("Start mixing GBath with GBathNew\n");
    TriDiagonalMatrixCopy(&GBathOld, GBath);
    FreeTriDiagonalMatrix(&GBath);
    GBathOld.N = Min(GBathOld.N, HSet.NGBath);
    strcpy(GBathOld.Name, "Gbath-Old");
    GBathOld.b[0] = sqrt(HSet.Mix);
    GBathNew.b[0] = sqrt(1 - HSet.Mix);
    SumGreensfunctionsTriDiagonalRepresentation(GBathNew, GBathOld, &GBath);
    strcpy(GBath.Name, "Gbath");
    GBathOld.b[0] = 1;
    GBathNew.b[0] = 1;
    /* symmetrize */
    SymmetrizeTridiagonalMatrix(&GBath);
    EndTiming("G math", &Time);
    /* create spectra */
    printf("Start create spectra\n");
    fflush(stdout);
    Spectra[0].N = HSet.NE;
    Spectra[0].Emin = HSet.Emin;
    Spectra[0].Emax = HSet.Emax;
    Spectra[0].dE = (HSet.Emax - HSet.Emin) / HSet.NE;
    Spectra[0].EminPole = HSet.Emin;
    Spectra[0].EmaxPole = HSet.Emax;
    Spectra[0].Gamma = HSet.Gamma;
    for (i = 0; i < 10; i++)
    {
      SpectrumAllocate(Spectra[0], &Spectra[i]);
    }
    TriDiagonalMatrixToSpectrum(0, GBath, &(Spectra[0])); /*  3 */
    TriDiagonalMatrixToSpectrum(0, GBathOld, &(Spectra[1]));  /*  5 */
    TriDiagonalMatrixToSpectrum(0, GBathNew, &(Spectra[2]));  /*  7 */
    SecondDoubleTriDiagonalMatrixToSpectrum(0, GAnd, &(Spectra[3]));  /*  9 */
    TriDiagonalMatrixToSpectrum(0, G0, &(Spectra[4]));  /* 11 */
    TriDiagonalMatrixToSpectrum(0, Sigma, &(Spectra[5])); /* 13 */
    AndersonImpurityMatrixToSpectrum(0, GBathNewAnd, &(Spectra[6]));  /* 15 */
    TriDiagonalMatrixToSpectrum(0, GclPup, &(Spectra[7]));  /* 17 */
    TriDiagonalMatrixToSpectrum(0, GclMup, &(Spectra[8]));  /* 19 */
    TriDiagonalMatrixToSpectrum(0, Gcl, &(Spectra[9])); /* 21 */
    /* write to disk */
    printf("Start write to dist\n");
    fflush(stdout);
    if (iOuter == 0)
    {
      sprintf(GBathFileName, "%s001", HSet.GRootFileName);
      sprintf(GClFileName, "Cl%s001", HSet.GRootFileName);
      sprintf(GSigmaFileName, "Sigma%s001", HSet.GRootFileName);
      sprintf(PsiFileName, "%s001", HSet.WaveFunctionRootFileName);
      sprintf(SpecFileName, "%s001", HSet.SpecRootFileName);
    }
    else if (iOuter < 9)
    {
      sprintf(GBathFileName, "%s00%1i", HSet.GRootFileName, iOuter + 1);
      sprintf(GClFileName, "Cl%s00%1i", HSet.GRootFileName, iOuter + 1);
      sprintf(GSigmaFileName, "Sigma%s00%1i", HSet.GRootFileName, iOuter + 1);
      sprintf(PsiFileName, "%s00%1i", HSet.WaveFunctionRootFileName, iOuter + 1);
      sprintf(SpecFileName, "%s00%1i", HSet.SpecRootFileName, iOuter + 1);
    }
    else if (iOuter < 99)
    {
      sprintf(GBathFileName, "%s0%2i", HSet.GRootFileName, iOuter + 1);
      sprintf(GClFileName, "Cl%s0%2i", HSet.GRootFileName, iOuter + 1);
      sprintf(GSigmaFileName, "Sigma%s0%2i", HSet.GRootFileName, iOuter + 1);
      sprintf(PsiFileName, "%s0%2i", HSet.WaveFunctionRootFileName, iOuter + 1);
      sprintf(SpecFileName, "%s0%2i", HSet.SpecRootFileName, iOuter + 1);
    }
    else if (iOuter < 999)
    {
      sprintf(GBathFileName, "%s%3i", HSet.GRootFileName, iOuter + 1);
      sprintf(GClFileName, "Cl%s%3i", HSet.GRootFileName, iOuter + 1);
      sprintf(GSigmaFileName, "Sigma%s%3i", HSet.GRootFileName, iOuter + 1);
      sprintf(PsiFileName, "%s%3i", HSet.WaveFunctionRootFileName, iOuter + 1);
      sprintf(SpecFileName, "%s%3i", HSet.SpecRootFileName, iOuter + 1);
    }
    FilePrintSpectra(Spectra, 10, SpecFileName, "w");
    FilePrintPlainTriDiagonalMatrix(GBath, GBathFileName, "w");
    FilePrintPlainTriDiagonalMatrix(Gcl, GClFileName, "w");
    FilePrintPlainTriDiagonalMatrix(Sigma, GSigmaFileName, "w");
    FilePrintWaveFunctionPlainBinary(Psi, PsiFileName, "w");
    /* free variables */
    printf("Start free variables\n");
    fflush(stdout);
    FreeTriDiagonalMatrix(&GclMup);
    FreeTriDiagonalMatrix(&GclMdn);
    FreeTriDiagonalMatrix(&GclPup);
    FreeTriDiagonalMatrix(&GclPdn);
    FreeTriDiagonalMatrix(&Gcl);
    FreeListOfPoles(&SigmaPoles);
    FreeTriDiagonalMatrix(&Sigma);
    FreeTriDiagonalMatrix(&GBathNew);
    FreeTriDiagonalMatrix(&GBathOld);
    FreeDoubleTriDiagonalMatrix(&(GAnd.Dt));
    FreeAndersonMatrix(&GBathNewAnd);
    FreeAndersonMatrix(&GBathAnd);
    EndTiming("DMFTLoop", &Time);
    PrintTiming(Time);
    if (HSet.Dyson)
      for (i = 0; i < 4; i++)
      {
        FreeOperator(&(VPerturbation[i]));
        FreeTriDiagonalMatrix(&(gPerturbation[i]));
      }
    FreeOperator(&HAnd);
    FreeOperator(&OPESdn);
    FreeOperator(&OPESup);
    FreeOperator(&OIPESdn);
    FreeOperator(&OIPESup);
  }
  EndTiming("Total", &Time);
  PrintTiming(Time);
  return 0;
}
