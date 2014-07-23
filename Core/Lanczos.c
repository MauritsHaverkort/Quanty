/*
*   Lanczos Subroutines
*   Functions Specified:
*     LanczosGroundState returns the lowest eigenvector of a given hamiltonian
*     LanczosTriDiagonalize returns a list containing a and b of a tridiagonal matrix given an hamiltonian and initial wave function
*
*/


/*  Different versions of LanczosGroundState depending on QComplex of H and psi
 *
 */

uint32_t LanczosGroundStateRR(const OperatorType H, WaveFunctionType * const psi, double *epsilon, const double zero, double shiftin, WaveFunctionType * OtherPsi, uint32_t NOtherPsi)
{
  printf("\nStart of LanczosGroundStateRR\n\n");
  fflush(stdout);
  const uint32_t NKrylovMax = 100, NKrylovStart = 25;
  uint32_t NKrylov;
  WaveFunctionType v[3];
  double energy = 4 * zero, energy1 = 2 * zero, energy2, energysqr;
  double shift;
  uint32_t cycles, i, vi, vip, vim, minloop = 3, statetopick;
  MatrixType Eigen;
  TriDiagonalMatrixType M;
  double val[NKrylovMax];
  M.N = NKrylovMax;
  InitTriDiagonalMatrix(&M);
  // Orthogonalize to previously found functions
  RealWaveFunctionOrthogonalizeConserveBasis(psi, OtherPsi, NOtherPsi);
  RealWaveFunctionNormalize(psi);
  // Shift H down to twice the energy of the starting wave-function
  // for the case that this is zero, add an aditional shift of 1
  // Make sure at least 2 loops are done
  PsiOperatorPsiRR(H, *psi, &energy1);
  PsiOperatorSqrPsiRR(H, *psi, &energysqr);
  if (energy1 < 0)
  {
    shift = -Abs(energy) - 1;
  }
  else
  {
    shift = -2 * energy - 1;
  }

  // loop until the relative energy gain is smaller than zero.
  printf("Energy: %22.15E variance %22.15E Number of determinants %d\n", energy1,energysqr-Sqr(energy1),psi->N);
  while ((((energy - energy1) > (zero * (1.0 + fabs(energy)))) || minloop)&&((energysqr-Sqr(energy1))>zero*(fabs(energy1)+DBL_EPSILON)))
  {
    if (minloop)
      minloop--;
    cycles = 0;
    energy = energy1;
    // Calculate psi_new = H^3 psi_old whereby the first time H acts on psi_old the basis size is enlarged
    //  This leaves three functions allocated which can be used in the next part
  Repeat:
    if (OperatorPsiRR(H, *psi, &v[0]))
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
    
    // Orthogonalize to previously found functions
    RealWaveFunctionOrthogonalizeConserveBasis(psi, OtherPsi, NOtherPsi);
    RealWaveFunctionNormalize(psi);
    
    // For a given basis set determine the lowest state
    //  In the previous part the memory for 3 functions v[0], v[1] and v[2] is allocated on this basis set.
    energy1 = 2 * zero;
    energy2 = 4 * zero;
    // loop until the relative energy gain is smaller than zero, loop a minimum of 3 times
    
    while ((energy2 - energy1) > (zero * (1.0 + fabs(energy1))))
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
      
      // test
      //PrintMatrix(Eigen);
      
      
      // Given the eigen-states in the Krylov basis reconstruct the ground-state in the original basis
      // pick the lowest state which has an overlap larger than sqrt(DBL_EPSILON) with the previous state
      statetopick=0;
      for(i=0;i<M.N;i++)
      {
        if(Sqr(Eigen.VRe[i][0])>DBL_EPSILON)
        {
          statetopick=i;
          break;
        }
      }
      
      
      
      RealWaveFunctionMultiply(v[0], 0);
      RealWaveFunctionAddSameBasis(v[0], *psi);
      RealWaveFunctionMultiply(*psi, Eigen.VRe[statetopick][0]);
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
        RealWaveFunctionAddMultipliedSameBasis(*psi, v[vi], Eigen.VRe[statetopick][i]);
      }
      // Set the shift such that the largest eigenvalue found becomes zero
      shift -= val[M.N - 1];
      shift -= shiftin;
      // Orthogonalize to previously found functions
      
      RealWaveFunctionOrthogonalizeConserveBasis(psi, OtherPsi, NOtherPsi);
      RealWaveFunctionNormalize(psi);
      PsiOperatorPsiRR(H, *psi, &energy1);
      PsiOperatorSqrPsiRR(H, *psi, &energysqr);
    }
    RealWaveFunctionFreeValues(&v[0]);
    RealWaveFunctionFreeValues(&v[1]);
    RealWaveFunctionFreeValues(&v[2]);
    RealWaveFunctionTruncateBasis(psi, DBL_EPSILON);
    RealWaveFunctionNormalize(psi);
    FreeMatrix(&Eigen);
    printf("Energy: %22.15E variance %22.15E ", energy1,energysqr-Sqr(energy1));
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

uint32_t LanczosGroundStateRC(const OperatorType H, WaveFunctionType * const psi, double *epsilon, const double zero, double shiftin, WaveFunctionType * OtherPsi, uint32_t NOtherPsi)
{
  printf("\nStart of LanczosGroundStateRC\n\n");
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
    if (OperatorPsiRC(H, *psi, &v[0]))
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

uint32_t LanczosGroundStateMC(const OperatorType H, WaveFunctionType * const psi, double *epsilon, const double zero, double shiftin, WaveFunctionType * OtherPsi, uint32_t NOtherPsi)
{
  printf("\nStart of LanczosGroundStateMC\n\n");
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
    if (OperatorPsiMC(H, *psi, &v[0]))
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

uint32_t LanczosGroundStateCC(const OperatorType H, WaveFunctionType * const psi, double *epsilon, const double zero, double shiftin, WaveFunctionType * OtherPsi, uint32_t NOtherPsi)
{
  printf("\nStart of LanczosGroundStateCC\n\n");
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
    if (OperatorPsiCC(H, *psi, &v[0]))
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

uint32_t LanczosGroundStateCR(const OperatorType H, WaveFunctionType * const psi, double *epsilon, const double zero, double shiftin, WaveFunctionType * OtherPsi, uint32_t NOtherPsi)
{
  if (MakeComplexWaveFunctionFromReal(psi, *psi))
  {
    printf
      ("Lanczos called with a real wave-function and a complex Hamiltonian.\nTurned the wave-function complex in step 1.\nBut failed.\n");
    fflush(stdout);
    return 1;
  }
  return LanczosGroundStateCC(H, psi, epsilon, zero, shiftin, OtherPsi, NOtherPsi);
}

uint32_t LanczosGroundStateMR(const OperatorType H, WaveFunctionType * const psi, double *epsilon, const double zero, double shiftin, WaveFunctionType * OtherPsi, uint32_t NOtherPsi)
{
  if (MakeComplexWaveFunctionFromReal(psi, *psi))
  {
    printf
      ("Lanczos called with a real wave-function and a mixed complex Hamiltonian.\nTurned the wave-function complex in step 1.\nBut failed.\n");
    fflush(stdout);
    return 1;
  }
  return LanczosGroundStateMC(H, psi, epsilon, zero, shiftin, OtherPsi, NOtherPsi);
}






/*  Different versions of LanczosGroundStateRestricted depending on QComplex of H and psi
 *
 */
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

uint32_t LanczosGroundStateRestrictedCR(const OperatorType H, WaveFunctionType * const psi, DeterminantRestrictionType const Restrictions, double *epsilon, const double zero, double shiftin, WaveFunctionType * OtherPsi, uint32_t NOtherPsi)
{
  if (MakeComplexWaveFunctionFromReal(psi, *psi))
  {
    printf
    ("LanczosGroundStateRestrictedCR called with a real wave-function and a complex Hamiltonian.\nTurned the wave-function complex in step 1.\nBut failed.\n");
    fflush(stdout);
    return 1;
  }
  return LanczosGroundStateRestrictedCC(H, psi, Restrictions, epsilon, zero, shiftin, OtherPsi, NOtherPsi);
}

uint32_t LanczosGroundStateRestrictedMR(const OperatorType H, WaveFunctionType * const psi, DeterminantRestrictionType const Restrictions, double *epsilon, const double zero, double shiftin, WaveFunctionType * OtherPsi, uint32_t NOtherPsi)
{
  if (MakeComplexWaveFunctionFromReal(psi, *psi))
  {
    printf
    ("LanczosGroundStateRestrictedMR called with a real wave-function and a mixed complex Hamiltonian.\nTurned the wave-function complex in step 1.\nBut failed.\n");
    fflush(stdout);
    return 1;
  }
  return LanczosGroundStateRestrictedMC(H, psi, Restrictions, epsilon, zero, shiftin, OtherPsi, NOtherPsi);
}





/*  Different versions of LanczosTriDiagonalizedepending on QComplex of H and psi
 *
 */
uint32_t LanczosTriDiagonalizeRR(const OperatorType H, WaveFunctionType psi, double epsilon, TriDiagonalMatrixType * const M, double Gamma)
{
  printf("\nStart of LanczosTriDiagonalizeRR\n\n");
  fflush(stdout);
  WaveFunctionType v[3];
  uint32_t i, vi=0, vip=0, vim=0;
  double temp, errorestimate = 1.0;
  double QRe[3], QIm[3], PRe[3], PIm[3], ValEst[3];
  QRe[0] = 1;
  QRe[1] = 0;
  QRe[2] = 0;
  QIm[0] = 0;
  QIm[1] = 0;
  QIm[2] = 0;
  PRe[0] = 0;
  PRe[1] = 0;
  PRe[2] = 1;
  PIm[0] = 0;
  PIm[1] = 0;
  PIm[2] = 0;
  ValEst[0] = 0.0;
  ValEst[1] = 0.0;
  ValEst[2] = 0.0;
  /* copy the name of psi to v[0] in order to do the loop more cleanely  */
  v[0].N = psi.N;
  v[0].NBlok = psi.NBlok;
  v[0].NFermion = psi.NFermion;
  v[0].NBoson = psi.NBoson;
  v[0].NWord = psi.NWord;
  v[0].NBitsKey = psi.NBitsKey;
  v[0].ValueRe = psi.ValueRe;
  v[0].Determinant = psi.Determinant;
  v[0].Next = psi.Next;
  v[0].KeyList = psi.KeyList;
  v[0].QComplex = psi.QComplex;
  /* finished copying */
  M->a[0] = 0;
  M->b[0] = RealWaveFunctionCalculateNorm(v[0]);
  if (M->b[0] == 0)
  {
    RealWaveFunctionMultiply(v[0], (double)0.0);
  }
  else
  {
    RealWaveFunctionMultiply(v[0], (double)1.0 / M->b[0]);
  }
  printf("BasisSetSize: %12i ", v[0].N);
  RealWaveFunctionTruncateBasis(&(v[0]), epsilon);
  printf("%12i\n", v[0].N);
  while (OperatorPsiRR(H, v[0], &(v[1])))
  {                             /* Allocation of Memory failed */
    printf("Allocation of new memory failed in LanczosTriDiagonalize for the function OperatorPsi\n");
    printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
           epsilon * 2.0);
    fflush(stdout);
    epsilon *= 2.0;
    if (v[0].NBlok == 1)
      AllocFailed();
    RealWaveFunctionTruncateBasis(&(v[0]), epsilon);
  }
  M->a[1] = RealWaveFunctionDotProduct(&v[1], &v[0]);
  printf("BasisSetSize: %12i ", v[1].N);
  for (i = 1; i < M->N; i++)
  {
    vip = (i + 1) % 3;
    vi = (i) % 3;
    vim = (i - 1) % 3;
    while (RealWaveFunctionAddMultipliedUnionBasis(&(v[vi]), v[vim], -M->a[i]))
    {                           /*Allocation of Memory failed */
      printf
        ("Allocation of new memory failed in LanczosTriDiagonalize for the function RealWaveFunctionAddMultipliedUnionBasis(&(v[vi]),v[vim],-M->a[i])\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
             epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      temp = RealWaveFunctionCalculateNorm(v[vi]);
      if (v[vi].NBlok == 1)
        AllocFailed();
      RealWaveFunctionTruncateBasis(&(v[vi]), epsilon * temp / errorestimate);
    }
    M->b[i] = RealWaveFunctionCalculateNorm(v[vi]);
    if (M->b[i] == 0)
    {
      RealWaveFunctionMultiply(v[vi], (double)0.0);
    }
    else
    {
      RealWaveFunctionMultiply(v[vi], (double)1.0 / M->b[i]);
    }
    RealWaveFunctionTruncateBasis(&(v[vi]), epsilon / errorestimate);
    printf("%12i ", v[vi].N);
    printf("| M.a[%3i] M.b[%3i]: %12.5E %12.5E | Err: %12.5E\n", i, i - 1, M->a[i], M->b[i - 1], errorestimate);
    /*
       QRe[vi]=M->b[i-1]*QRe[vip]-0.5*Gamma*QIm[vim];
       QIm[vi]=M->b[i-1]*QIm[vip]+0.5*Gamma*QRe[vim];
       PRe[vi]=M->b[i-1]*PRe[vip]-0.5*Gamma*PIm[vim];
       PIm[vi]=M->b[i-1]*PIm[vip]+0.5*Gamma*PRe[vim];
       ValEst[vi]=fabs((PIm[vi]*QRe[vi]-PRe[vi]*QIm[vi])/(QIm[vi]*QIm[vi]+QRe[vi]*QRe[vi]));
       if(i>5)errorestimate=fabs(2*(ValEst[vi]-ValEst[vim])/(Gamma));
     */
    errorestimate *= M->b[i] * M->b[i] / (M->b[i] * M->b[i] + 4 * Gamma * Gamma);
    while (OperatorPsiRR(H, v[vi], &(v[vip])))
    {                           /*Allocation of Memory failed */
      printf("Allocation of new memory failed in LanczosTriDiagonalize for the function OperatorPsi\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
             epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      if (v[vi].NBlok == 1)
        AllocFailed();
      RealWaveFunctionTruncateBasis(&(v[vi]), epsilon / errorestimate);
    }
    printf("BasisSetSize: %12i ", v[vip].N);
    M->a[i + 1] = RealWaveFunctionDotProduct(&v[vip], &v[vi]);
    while (RealWaveFunctionAddMultipliedUnionBasis(&(v[vip]), v[vim], -M->b[i]))
    {                           /*Allocation of Memory failed */
      printf
        ("Allocation of new memory failed in LanczosTriDiagonalize for the function RealWaveFunctionAddMultipliedUnionBasis(&(v[vip]),v[vim],-M->b[i])\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
             epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      temp = RealWaveFunctionCalculateNorm(v[vip]);
      if (v[vip].NBlok == 1)
        AllocFailed();
      RealWaveFunctionTruncateBasis(&(v[vip]), epsilon * temp / errorestimate);
    }
    RealWaveFunctionFree(&v[vim]);
  }
  printf("             | M.a[%3i] M.b[%3i]: %12.5E %12.5E | Err: %12.5E\n\n", i, i - 1, M->a[i], M->b[i - 1],
         errorestimate);
  if (M->N > 1)
  {
    RealWaveFunctionFree(&v[vi]);
    RealWaveFunctionFree(&v[vip]);
  }
  fflush(stdout);
  return 0;
}

uint32_t LanczosTriDiagonalizeCR(const OperatorType H, WaveFunctionType psi, double epsilon,
                                 TriDiagonalMatrixType * const M, double Gamma)
{
  uint32_t retval;
  WaveFunctionType Cpsi;
  if (MakeComplexWaveFunctionFromReal(&Cpsi, psi))
  {
    printf
      ("LanczosTriDiagonalizeCR called with a real wave-function and a complex Hamiltonian.\nTurned the wave-function complex in step 1.\nBut failed.\n");
    fflush(stdout);
    return 1;
  }
  retval = LanczosTriDiagonalizeCC(H, Cpsi, epsilon, M, Gamma);
  MakeRealWaveFunctionFromComplex(Cpsi, &Cpsi);
  return retval;
}

uint32_t LanczosTriDiagonalizeMR(const OperatorType H, WaveFunctionType psi, double epsilon,
                                 TriDiagonalMatrixType * const M, double Gamma)
{
  uint32_t retval;
  WaveFunctionType Cpsi;
  if (MakeComplexWaveFunctionFromReal(&Cpsi, psi))
  {
    printf
      ("LanczosTriDiagonalizeMR called with a real wave-function and a mixed complex Hamiltonian.\nTurned the wave-function complex in step 1.\nBut failed.\n");
    fflush(stdout);
    return 1;
  }
  retval = LanczosTriDiagonalizeMC(H, Cpsi, epsilon, M, Gamma);
  MakeRealWaveFunctionFromComplex(Cpsi, &Cpsi);
  return retval;
}

uint32_t LanczosTriDiagonalizeRC(const OperatorType H, WaveFunctionType psi, double epsilon,
                                 TriDiagonalMatrixType * const M, double Gamma)
{
  printf("\nStart of LanczosTriDiagonalizeRC\n\n");
  fflush(stdout);
  WaveFunctionType v[3];
  uint32_t i, vi=0, vip=0, vim=0;
  double temp, errorestimate = 1.0, junk;
  /* copy the name of psi to v[0] in order to do the loop more cleanely  */
  v[0].N = psi.N;
  v[0].NBlok = psi.NBlok;
  v[0].NFermion = psi.NFermion;
  v[0].NBoson = psi.NBoson;
  v[0].NWord = psi.NWord;
  v[0].NBitsKey = psi.NBitsKey;
  v[0].ValueRe = psi.ValueRe;
  v[0].ValueIm = psi.ValueIm;
  v[0].Determinant = psi.Determinant;
  v[0].Next = psi.Next;
  v[0].KeyList = psi.KeyList;
  v[0].QComplex = psi.QComplex;
  /* finished copying */
  M->a[0] = 0;
  M->b[0] = ComplexWaveFunctionCalculateNorm(v[0]);
  ComplexWaveFunctionMultiply(v[0], (double)1.0 / M->b[0]);
  printf("BasisSetSize: %12i ", v[0].N);
  ComplexWaveFunctionTruncateBasis(&(v[0]), epsilon);
  printf("%12i\n", v[0].N);
  while (OperatorPsiRC(H, v[0], &(v[1])))
  {                             /* Allocation of Memory failed */
    printf("Allocation of new memory failed in LanczosTriDiagonalize for the function OperatorPsi\n");
    printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
           epsilon * 2.0);
    fflush(stdout);
    epsilon *= 2.0;
    if (v[0].NBlok == 1)
      AllocFailed();
    ComplexWaveFunctionTruncateBasis(&(v[0]), epsilon);
  }
  ComplexWaveFunctionDotProduct(&v[1], &v[0], &M->a[1], &junk);
  printf("BasisSetSize: %12i ", v[1].N);
  for (i = 1; i < M->N; i++)
  {
    vip = (i + 1) % 3;
    vi = (i) % 3;
    vim = (i - 1) % 3;
    while (ComplexWaveFunctionAddMultipliedUnionBasis(&(v[vi]), v[vim], -M->a[i]))
    {                           /*Allocation of Memory failed */
      printf
        ("Allocation of new memory failed in LanczosTriDiagonalize for the function RealWaveFunctionAddMultipliedUnionBasis(&(v[vi]),v[vim],-M->a[i])\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
             epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      temp = ComplexWaveFunctionCalculateNorm(v[vi]);
      if (v[vi].NBlok == 1)
        AllocFailed();
      ComplexWaveFunctionTruncateBasis(&(v[vi]), epsilon * temp / errorestimate);
    }
    M->b[i] = ComplexWaveFunctionCalculateNorm(v[vi]);
    ComplexWaveFunctionMultiply(v[vi], (double)1.0 / M->b[i]);
    ComplexWaveFunctionTruncateBasis(&(v[vi]), epsilon / errorestimate);
    printf("%12i ", v[vi].N);
    printf("| M.a[%3i] M.b[%3i]: %12.5E %12.5E | Err: %12.5E\n", i, i - 1, M->a[i], M->b[i - 1], errorestimate);
    errorestimate *= M->b[i] * M->b[i] / (M->b[i] * M->b[i] + 4 * Gamma * Gamma);
    while (OperatorPsiRC(H, v[vi], &(v[vip])))
    {                           /*Allocation of Memory failed */
      printf("Allocation of new memory failed in LanczosTriDiagonalize for the function OperatorPsi\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
             epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      if (v[vi].NBlok == 1)
        AllocFailed();
      ComplexWaveFunctionTruncateBasis(&(v[vi]), epsilon / errorestimate);
    }
    printf("BasisSetSize: %12i ", v[vip].N);
    ComplexWaveFunctionDotProduct(&v[vip], &v[vi], &M->a[i + 1], &junk);
    while (ComplexWaveFunctionAddMultipliedUnionBasis(&(v[vip]), v[vim], -M->b[i]))
    {                           /*Allocation of Memory failed */
      printf
        ("Allocation of new memory failed in LanczosTriDiagonalize for the function RealWaveFunctionAddMultipliedUnionBasis(&(v[vip]),v[vim],-M->b[i])\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
             epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      temp = ComplexWaveFunctionCalculateNorm(v[vip]);
      if (v[vip].NBlok == 1)
        AllocFailed();
      ComplexWaveFunctionTruncateBasis(&(v[vip]), epsilon * temp / errorestimate);
    }
    ComplexWaveFunctionFree(&v[vim]);
  }
  printf("             | M.a[%3i] M.b[%3i]: %12.5E %12.5E | Err: %12.5E\n\n", i, i - 1, M->a[i], M->b[i - 1],
         errorestimate);
  if (M->N > 1)
  {
    ComplexWaveFunctionFree(&v[vi]);
    ComplexWaveFunctionFree(&v[vip]);
  }
  fflush(stdout);
  return 0;
}

uint32_t LanczosTriDiagonalizeCC(const OperatorType H, WaveFunctionType psi, double epsilon,
                                 TriDiagonalMatrixType * const M, double Gamma)
{
  printf("\nStart of LanczosTriDiagonalizeCC\n\n");
  fflush(stdout);
  WaveFunctionType v[3];
  uint32_t i, vi=0, vip=0, vim=0;
  double temp, errorestimate = 1.0, junk;
  /* copy the name of psi to v[0] in order to do the loop more cleanely  */
  v[0].N = psi.N;
  v[0].NBlok = psi.NBlok;
  v[0].NFermion = psi.NFermion;
  v[0].NBoson = psi.NBoson;
  v[0].NWord = psi.NWord;
  v[0].NBitsKey = psi.NBitsKey;
  v[0].ValueRe = psi.ValueRe;
  v[0].ValueIm = psi.ValueIm;
  v[0].Determinant = psi.Determinant;
  v[0].Next = psi.Next;
  v[0].KeyList = psi.KeyList;
  v[0].QComplex = psi.QComplex;
  /* finished copying */
  M->a[0] = 0;
  M->b[0] = ComplexWaveFunctionCalculateNorm(v[0]);
  ComplexWaveFunctionMultiply(v[0], (double)1.0 / M->b[0]);
  printf("BasisSetSize: %12i ", v[0].N);
  ComplexWaveFunctionTruncateBasis(&(v[0]), epsilon);
  printf("%12i\n", v[0].N);
  while (OperatorPsiCC(H, v[0], &(v[1])))
  {                             /* Allocation of Memory failed */
    printf("Allocation of new memory failed in LanczosTriDiagonalize for the function OperatorPsi\n");
    printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
           epsilon * 2.0);
    fflush(stdout);
    epsilon *= 2.0;
    if (v[0].NBlok == 1)
      AllocFailed();
    ComplexWaveFunctionTruncateBasis(&(v[0]), epsilon);
  }
  ComplexWaveFunctionDotProduct(&v[1], &v[0], &M->a[1], &junk);
  printf("BasisSetSize: %12i ", v[1].N);
  for (i = 1; i < M->N; i++)
  {
    vip = (i + 1) % 3;
    vi = (i) % 3;
    vim = (i - 1) % 3;
    while (ComplexWaveFunctionAddMultipliedUnionBasis(&(v[vi]), v[vim], -M->a[i]))
    {                           /*Allocation of Memory failed */
      printf
        ("Allocation of new memory failed in LanczosTriDiagonalize for the function RealWaveFunctionAddMultipliedUnionBasis(&(v[vi]),v[vim],-M->a[i])\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
             epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      temp = ComplexWaveFunctionCalculateNorm(v[vi]);
      if (v[vi].NBlok == 1)
        AllocFailed();
      ComplexWaveFunctionTruncateBasis(&(v[vi]), epsilon * temp / errorestimate);
    }
    M->b[i] = ComplexWaveFunctionCalculateNorm(v[vi]);
    ComplexWaveFunctionMultiply(v[vi], (double)1.0 / M->b[i]);
    ComplexWaveFunctionTruncateBasis(&(v[vi]), epsilon / errorestimate);
    printf("%12i ", v[vi].N);
    printf("| M.a[%3i] M.b[%3i]: %12.5E %12.5E | Err: %12.5E\n", i, i - 1, M->a[i], M->b[i - 1], errorestimate);
    errorestimate *= M->b[i] * M->b[i] / (M->b[i] * M->b[i] + 4 * Gamma * Gamma);
    while (OperatorPsiCC(H, v[vi], &(v[vip])))
    {                           /*Allocation of Memory failed */
      printf("Allocation of new memory failed in LanczosTriDiagonalize for the function OperatorPsi\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
             epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      if (v[vi].NBlok == 1)
        AllocFailed();
      ComplexWaveFunctionTruncateBasis(&(v[vi]), epsilon / errorestimate);
    }
    printf("BasisSetSize: %12i ", v[vip].N);
    ComplexWaveFunctionDotProduct(&v[vip], &v[vi], &M->a[i + 1], &junk);
    while (ComplexWaveFunctionAddMultipliedUnionBasis(&(v[vip]), v[vim], -M->b[i]))
    {                           /*Allocation of Memory failed */
      printf
        ("Allocation of new memory failed in LanczosTriDiagonalize for the function RealWaveFunctionAddMultipliedUnionBasis(&(v[vip]),v[vim],-M->b[i])\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
             epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      temp = ComplexWaveFunctionCalculateNorm(v[vip]);
      if (v[vip].NBlok == 1)
        AllocFailed();
      ComplexWaveFunctionTruncateBasis(&(v[vip]), epsilon * temp / errorestimate);
    }
    ComplexWaveFunctionFree(&v[vim]);
  }
  printf("             | M.a[%3i] M.b[%3i]: %12.5E %12.5E | Err: %12.5E\n\n", i, i - 1, M->a[i], M->b[i - 1],
         errorestimate);
  if (M->N > 1)
  {
    ComplexWaveFunctionFree(&v[vi]);
    ComplexWaveFunctionFree(&v[vip]);
  }
  fflush(stdout);
  return 0;
}

uint32_t LanczosTriDiagonalizeMC(const OperatorType H, WaveFunctionType psi, double epsilon,
                                 TriDiagonalMatrixType * const M, double Gamma)
{
  printf("\nStart of LanczosTriDiagonalizeMC\n\n");
  fflush(stdout);
  WaveFunctionType v[3];
  uint32_t i, vi=0, vip=0, vim=0;
  double temp, errorestimate = 1.0, junk;
  /* copy the name of psi to v[0] in order to do the loop more cleanely  */
  v[0].N = psi.N;
  v[0].NBlok = psi.NBlok;
  v[0].NFermion = psi.NFermion;
  v[0].NBoson = psi.NBoson;
  v[0].NWord = psi.NWord;
  v[0].NBitsKey = psi.NBitsKey;
  v[0].ValueRe = psi.ValueRe;
  v[0].ValueIm = psi.ValueIm;
  v[0].Determinant = psi.Determinant;
  v[0].Next = psi.Next;
  v[0].KeyList = psi.KeyList;
  v[0].QComplex = psi.QComplex;
  /* finished copying */
  M->a[0] = 0;
  M->b[0] = ComplexWaveFunctionCalculateNorm(v[0]);
  ComplexWaveFunctionMultiply(v[0], (double)1.0 / M->b[0]);
  printf("BasisSetSize: %12i ", v[0].N);
  ComplexWaveFunctionTruncateBasis(&(v[0]), epsilon);
  printf("%12i\n", v[0].N);
  while (OperatorPsiMC(H, v[0], &(v[1])))
  {                             /* Allocation of Memory failed */
    printf("Allocation of new memory failed in LanczosTriDiagonalize for the function OperatorPsi\n");
    printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
           epsilon * 2.0);
    fflush(stdout);
    epsilon *= 2.0;
    if (v[0].NBlok == 1)
      AllocFailed();
    ComplexWaveFunctionTruncateBasis(&(v[0]), epsilon);
  }
  ComplexWaveFunctionDotProduct(&v[1], &v[0], &M->a[1], &junk);
  printf("BasisSetSize: %12i ", v[1].N);
  for (i = 1; i < M->N; i++)
  {
    vip = (i + 1) % 3;
    vi = (i) % 3;
    vim = (i - 1) % 3;
    while (ComplexWaveFunctionAddMultipliedUnionBasis(&(v[vi]), v[vim], -M->a[i]))
    {                           /*Allocation of Memory failed */
      printf
        ("Allocation of new memory failed in LanczosTriDiagonalize for the function RealWaveFunctionAddMultipliedUnionBasis(&(v[vi]),v[vim],-M->a[i])\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
             epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      temp = ComplexWaveFunctionCalculateNorm(v[vi]);
      if (v[vi].NBlok == 1)
        AllocFailed();
      ComplexWaveFunctionTruncateBasis(&(v[vi]), epsilon * temp / errorestimate);
    }
    M->b[i] = ComplexWaveFunctionCalculateNorm(v[vi]);
    ComplexWaveFunctionMultiply(v[vi], (double)1.0 / M->b[i]);
    ComplexWaveFunctionTruncateBasis(&(v[vi]), epsilon / errorestimate);
    printf("%12i ", v[vi].N);
    printf("| M.a[%3i] M.b[%3i]: %12.5E %12.5E | Err: %12.5E\n", i, i - 1, M->a[i], M->b[i - 1], errorestimate);
    errorestimate *= M->b[i] * M->b[i] / (M->b[i] * M->b[i] + 4 * Gamma * Gamma);
    while (OperatorPsiMC(H, v[vi], &(v[vip])))
    {                           /*Allocation of Memory failed */
      printf("Allocation of new memory failed in LanczosTriDiagonalize for the function OperatorPsi\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
             epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      if (v[vi].NBlok == 1)
        AllocFailed();
      ComplexWaveFunctionTruncateBasis(&(v[vi]), epsilon / errorestimate);
    }
    printf("BasisSetSize: %12i ", v[vip].N);
    ComplexWaveFunctionDotProduct(&v[vip], &v[vi], &M->a[i + 1], &junk);
    while (ComplexWaveFunctionAddMultipliedUnionBasis(&(v[vip]), v[vim], -M->b[i]))
    {                           /*Allocation of Memory failed */
      printf
        ("Allocation of new memory failed in LanczosTriDiagonalize for the function RealWaveFunctionAddMultipliedUnionBasis(&(v[vip]),v[vim],-M->b[i])\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
             epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      temp = ComplexWaveFunctionCalculateNorm(v[vip]);
      if (v[vip].NBlok == 1)
        AllocFailed();
      ComplexWaveFunctionTruncateBasis(&(v[vip]), epsilon * temp / errorestimate);
    }
    ComplexWaveFunctionFree(&v[vim]);
  }
  printf("             | M.a[%3i] M.b[%3i]: %12.5E %12.5E | Err: %12.5E\n\n", i, i - 1, M->a[i], M->b[i - 1],
         errorestimate);
  if (M->N > 1)
  {
    ComplexWaveFunctionFree(&v[vi]);
    ComplexWaveFunctionFree(&v[vip]);
  }
  fflush(stdout);
  return 0;
}




/*  Different versions of LanczosTriDiagonalizeRestricted depending on QComplex of H and psi
 *
 */
uint32_t LanczosTriDiagonalizeRestrictedRR(const OperatorType H, WaveFunctionType psi, DeterminantRestrictionType const Restrictions, double epsilon, TriDiagonalMatrixType * const M, double Gamma)
{
  printf("\nStart of LanczosTriDiagonalizeRestrictedRR\n\n");
  fflush(stdout);
  WaveFunctionType v[3];
  uint32_t i, vi=0, vip=0, vim=0;
  double temp, errorestimate = 1.0;
  double QRe[3], QIm[3], PRe[3], PIm[3], ValEst[3];
  QRe[0] = 1;
  QRe[1] = 0;
  QRe[2] = 0;
  QIm[0] = 0;
  QIm[1] = 0;
  QIm[2] = 0;
  PRe[0] = 0;
  PRe[1] = 0;
  PRe[2] = 1;
  PIm[0] = 0;
  PIm[1] = 0;
  PIm[2] = 0;
  ValEst[0] = 0.0;
  ValEst[1] = 0.0;
  ValEst[2] = 0.0;
  /* copy the name of psi to v[0] in order to do the loop more cleanely  */
  v[0].N = psi.N;
  v[0].NBlok = psi.NBlok;
  v[0].NFermion = psi.NFermion;
  v[0].NBoson = psi.NBoson;
  v[0].NWord = psi.NWord;
  v[0].NBitsKey = psi.NBitsKey;
  v[0].ValueRe = psi.ValueRe;
  v[0].Determinant = psi.Determinant;
  v[0].Next = psi.Next;
  v[0].KeyList = psi.KeyList;
  v[0].QComplex = psi.QComplex;
  /* finished copying */
  M->a[0] = 0;
  M->b[0] = RealWaveFunctionCalculateNorm(v[0]);
  if (M->b[0] == 0)
  {
    RealWaveFunctionMultiply(v[0], (double)0.0);
  }
  else
  {
    RealWaveFunctionMultiply(v[0], (double)1.0 / M->b[0]);
  }
  printf("BasisSetSize: %12i ", v[0].N);
  RealWaveFunctionTruncateBasis(&(v[0]), epsilon);
  printf("%12i\n", v[0].N);
  while (OperatorPsiRestrictedRR(H, v[0], &(v[1]), Restrictions))
  {                             /* Allocation of Memory failed */
    printf("Allocation of new memory failed in LanczosTriDiagonalize for the function OperatorPsi\n");
    printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
           epsilon * 2.0);
    fflush(stdout);
    epsilon *= 2.0;
    if (v[0].NBlok == 1)
      AllocFailed();
    RealWaveFunctionTruncateBasis(&(v[0]), epsilon);
  }
  M->a[1] = RealWaveFunctionDotProduct(&v[1], &v[0]);
  printf("BasisSetSize: %12i ", v[1].N);
  for (i = 1; i < M->N; i++)
  {
    vip = (i + 1) % 3;
    vi = (i) % 3;
    vim = (i - 1) % 3;
    while (RealWaveFunctionAddMultipliedUnionBasis(&(v[vi]), v[vim], -M->a[i]))
    {                           /*Allocation of Memory failed */
      printf
      ("Allocation of new memory failed in LanczosTriDiagonalize for the function RealWaveFunctionAddMultipliedUnionBasis(&(v[vi]),v[vim],-M->a[i])\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
             epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      temp = RealWaveFunctionCalculateNorm(v[vi]);
      if (v[vi].NBlok == 1)
        AllocFailed();
      RealWaveFunctionTruncateBasis(&(v[vi]), epsilon * temp / errorestimate);
    }
    M->b[i] = RealWaveFunctionCalculateNorm(v[vi]);
    if (M->b[i] == 0)
    {
      RealWaveFunctionMultiply(v[vi], (double)0.0);
    }
    else
    {
      RealWaveFunctionMultiply(v[vi], (double)1.0 / M->b[i]);
    }
    RealWaveFunctionTruncateBasis(&(v[vi]), epsilon / errorestimate);
    printf("%12i ", v[vi].N);
    printf("| M.a[%3i] M.b[%3i]: %12.5E %12.5E | Err: %12.5E\n", i, i - 1, M->a[i], M->b[i - 1], errorestimate);
    /*
     QRe[vi]=M->b[i-1]*QRe[vip]-0.5*Gamma*QIm[vim];
     QIm[vi]=M->b[i-1]*QIm[vip]+0.5*Gamma*QRe[vim];
     PRe[vi]=M->b[i-1]*PRe[vip]-0.5*Gamma*PIm[vim];
     PIm[vi]=M->b[i-1]*PIm[vip]+0.5*Gamma*PRe[vim];
     ValEst[vi]=fabs((PIm[vi]*QRe[vi]-PRe[vi]*QIm[vi])/(QIm[vi]*QIm[vi]+QRe[vi]*QRe[vi]));
     if(i>5)errorestimate=fabs(2*(ValEst[vi]-ValEst[vim])/(Gamma));
     */
    errorestimate *= M->b[i] * M->b[i] / (M->b[i] * M->b[i] + 4 * Gamma * Gamma);
    while (OperatorPsiRestrictedRR(H, v[vi], &(v[vip]),Restrictions))
    {                           /*Allocation of Memory failed */
      printf("Allocation of new memory failed in LanczosTriDiagonalize for the function OperatorPsi\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
             epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      if (v[vi].NBlok == 1)
        AllocFailed();
      RealWaveFunctionTruncateBasis(&(v[vi]), epsilon / errorestimate);
    }
    printf("BasisSetSize: %12i ", v[vip].N);
    M->a[i + 1] = RealWaveFunctionDotProduct(&v[vip], &v[vi]);
    while (RealWaveFunctionAddMultipliedUnionBasis(&(v[vip]), v[vim], -M->b[i]))
    {                           /*Allocation of Memory failed */
      printf
      ("Allocation of new memory failed in LanczosTriDiagonalize for the function RealWaveFunctionAddMultipliedUnionBasis(&(v[vip]),v[vim],-M->b[i])\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
             epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      temp = RealWaveFunctionCalculateNorm(v[vip]);
      if (v[vip].NBlok == 1)
        AllocFailed();
      RealWaveFunctionTruncateBasis(&(v[vip]), epsilon * temp / errorestimate);
    }
    RealWaveFunctionFree(&v[vim]);
  }
  printf("             | M.a[%3i] M.b[%3i]: %12.5E %12.5E | Err: %12.5E\n\n", i, i - 1, M->a[i], M->b[i - 1],
         errorestimate);
  if (M->N > 1)
  {
    RealWaveFunctionFree(&v[vi]);
    RealWaveFunctionFree(&v[vip]);
  }
  fflush(stdout);
  return 0;
}

uint32_t LanczosTriDiagonalizeRestrictedCC(const OperatorType H, WaveFunctionType psi, DeterminantRestrictionType const Restrictions, double epsilon, TriDiagonalMatrixType * const M, double Gamma)
{
  printf("\nStart of LanczosTriDiagonalizeRestrictedCC\n\n");
  fflush(stdout);
  WaveFunctionType v[3];
  uint32_t i, vi=0, vip=0, vim=0;
  double temp, errorestimate = 1.0, junk;
  /* copy the name of psi to v[0] in order to do the loop more cleanely  */
  v[0].N = psi.N;
  v[0].NBlok = psi.NBlok;
  v[0].NFermion = psi.NFermion;
  v[0].NBoson = psi.NBoson;
  v[0].NWord = psi.NWord;
  v[0].NBitsKey = psi.NBitsKey;
  v[0].ValueRe = psi.ValueRe;
  v[0].ValueIm = psi.ValueIm;
  v[0].Determinant = psi.Determinant;
  v[0].Next = psi.Next;
  v[0].KeyList = psi.KeyList;
  v[0].QComplex = psi.QComplex;
  /* finished copying */
  M->a[0] = 0;
  M->b[0] = ComplexWaveFunctionCalculateNorm(v[0]);
  ComplexWaveFunctionMultiply(v[0], (double)1.0 / M->b[0]);
  printf("BasisSetSize: %12i ", v[0].N);
  ComplexWaveFunctionTruncateBasis(&(v[0]), epsilon);
  printf("%12i\n", v[0].N);
  while (OperatorPsiRestrictedCC(H, v[0], &(v[1]), Restrictions))
  {                             /* Allocation of Memory failed */
    printf("Allocation of new memory failed in LanczosTriDiagonalize for the function OperatorPsi\n");
    printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
           epsilon * 2.0);
    fflush(stdout);
    epsilon *= 2.0;
    if (v[0].NBlok == 1)
      AllocFailed();
    ComplexWaveFunctionTruncateBasis(&(v[0]), epsilon);
  }
  ComplexWaveFunctionDotProduct(&v[1], &v[0], &M->a[1], &junk);
  printf("BasisSetSize: %12i ", v[1].N);
  for (i = 1; i < M->N; i++)
  {
    vip = (i + 1) % 3;
    vi = (i) % 3;
    vim = (i - 1) % 3;
    while (ComplexWaveFunctionAddMultipliedUnionBasis(&(v[vi]), v[vim], -M->a[i]))
    {                           /*Allocation of Memory failed */
      printf
      ("Allocation of new memory failed in LanczosTriDiagonalize for the function RealWaveFunctionAddMultipliedUnionBasis(&(v[vi]),v[vim],-M->a[i])\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
             epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      temp = ComplexWaveFunctionCalculateNorm(v[vi]);
      if (v[vi].NBlok == 1)
        AllocFailed();
      ComplexWaveFunctionTruncateBasis(&(v[vi]), epsilon * temp / errorestimate);
    }
    M->b[i] = ComplexWaveFunctionCalculateNorm(v[vi]);
    ComplexWaveFunctionMultiply(v[vi], (double)1.0 / M->b[i]);
    ComplexWaveFunctionTruncateBasis(&(v[vi]), epsilon / errorestimate);
    printf("%12i ", v[vi].N);
    printf("| M.a[%3i] M.b[%3i]: %12.5E %12.5E | Err: %12.5E\n", i, i - 1, M->a[i], M->b[i - 1], errorestimate);
    errorestimate *= M->b[i] * M->b[i] / (M->b[i] * M->b[i] + 4 * Gamma * Gamma);
    while (OperatorPsiRestrictedCC(H, v[vi], &(v[vip]), Restrictions))
    {                           /*Allocation of Memory failed */
      printf("Allocation of new memory failed in LanczosTriDiagonalize for the function OperatorPsi\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
             epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      if (v[vi].NBlok == 1)
        AllocFailed();
      ComplexWaveFunctionTruncateBasis(&(v[vi]), epsilon / errorestimate);
    }
    printf("BasisSetSize: %12i ", v[vip].N);
    ComplexWaveFunctionDotProduct(&v[vip], &v[vi], &M->a[i + 1], &junk);
    while (ComplexWaveFunctionAddMultipliedUnionBasis(&(v[vip]), v[vim], -M->b[i]))
    {                           /*Allocation of Memory failed */
      printf
      ("Allocation of new memory failed in LanczosTriDiagonalize for the function RealWaveFunctionAddMultipliedUnionBasis(&(v[vip]),v[vim],-M->b[i])\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
             epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      temp = ComplexWaveFunctionCalculateNorm(v[vip]);
      if (v[vip].NBlok == 1)
        AllocFailed();
      ComplexWaveFunctionTruncateBasis(&(v[vip]), epsilon * temp / errorestimate);
    }
    ComplexWaveFunctionFree(&v[vim]);
  }
  printf("             | M.a[%3i] M.b[%3i]: %12.5E %12.5E | Err: %12.5E\n\n", i, i - 1, M->a[i], M->b[i - 1],
         errorestimate);
  if (M->N > 1)
  {
    ComplexWaveFunctionFree(&v[vi]);
    ComplexWaveFunctionFree(&v[vip]);
  }
  fflush(stdout);
  return 0;
}

uint32_t LanczosTriDiagonalizeRestrictedCR(const OperatorType H, WaveFunctionType psi, DeterminantRestrictionType const Restrictions, double epsilon, TriDiagonalMatrixType * const M, double Gamma)
{
  uint32_t retval;
  WaveFunctionType Cpsi;
  if (MakeComplexWaveFunctionFromReal(&Cpsi, psi))
  {
    printf
    ("LanczosTriDiagonalizeRestrictedCR called with a real wave-function and a complex Hamiltonian.\nTurned the wave-function complex in step 1.\nBut failed.\n");
    fflush(stdout);
    return 1;
  }
  retval = LanczosTriDiagonalizeRestrictedCC(H, Cpsi, Restrictions, epsilon, M, Gamma);
  MakeRealWaveFunctionFromComplex(Cpsi, &Cpsi);
  return retval;
}

uint32_t LanczosTriDiagonalizeRestrictedRC(const OperatorType H, WaveFunctionType psi, DeterminantRestrictionType const Restrictions, double epsilon, TriDiagonalMatrixType * const M, double Gamma)
{
  printf("\nStart of LanczosTriDiagonalizeRestrictedRC\n\n");
  fflush(stdout);
  WaveFunctionType v[3];
  uint32_t i, vi=0, vip=0, vim=0;
  double temp, errorestimate = 1.0, junk;
  /* copy the name of psi to v[0] in order to do the loop more cleanely  */
  v[0].N = psi.N;
  v[0].NBlok = psi.NBlok;
  v[0].NFermion = psi.NFermion;
  v[0].NBoson = psi.NBoson;
  v[0].NWord = psi.NWord;
  v[0].NBitsKey = psi.NBitsKey;
  v[0].ValueRe = psi.ValueRe;
  v[0].ValueIm = psi.ValueIm;
  v[0].Determinant = psi.Determinant;
  v[0].Next = psi.Next;
  v[0].KeyList = psi.KeyList;
  v[0].QComplex = psi.QComplex;
  /* finished copying */
  M->a[0] = 0;
  M->b[0] = ComplexWaveFunctionCalculateNorm(v[0]);
  ComplexWaveFunctionMultiply(v[0], (double)1.0 / M->b[0]);
  printf("BasisSetSize: %12i ", v[0].N);
  ComplexWaveFunctionTruncateBasis(&(v[0]), epsilon);
  printf("%12i\n", v[0].N);
  while (OperatorPsiRestrictedRC(H, v[0], &(v[1]),Restrictions))
  {                             /* Allocation of Memory failed */
    printf("Allocation of new memory failed in LanczosTriDiagonalize for the function OperatorPsi\n");
    printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
           epsilon * 2.0);
    fflush(stdout);
    epsilon *= 2.0;
    if (v[0].NBlok == 1)
      AllocFailed();
    ComplexWaveFunctionTruncateBasis(&(v[0]), epsilon);
  }
  ComplexWaveFunctionDotProduct(&v[1], &v[0], &M->a[1], &junk);
  printf("BasisSetSize: %12i ", v[1].N);
  for (i = 1; i < M->N; i++)
  {
    vip = (i + 1) % 3;
    vi = (i) % 3;
    vim = (i - 1) % 3;
    while (ComplexWaveFunctionAddMultipliedUnionBasis(&(v[vi]), v[vim], -M->a[i]))
    {                           /*Allocation of Memory failed */
      printf
      ("Allocation of new memory failed in LanczosTriDiagonalize for the function RealWaveFunctionAddMultipliedUnionBasis(&(v[vi]),v[vim],-M->a[i])\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
             epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      temp = ComplexWaveFunctionCalculateNorm(v[vi]);
      if (v[vi].NBlok == 1)
        AllocFailed();
      ComplexWaveFunctionTruncateBasis(&(v[vi]), epsilon * temp / errorestimate);
    }
    M->b[i] = ComplexWaveFunctionCalculateNorm(v[vi]);
    ComplexWaveFunctionMultiply(v[vi], (double)1.0 / M->b[i]);
    ComplexWaveFunctionTruncateBasis(&(v[vi]), epsilon / errorestimate);
    printf("%12i ", v[vi].N);
    printf("| M.a[%3i] M.b[%3i]: %12.5E %12.5E | Err: %12.5E\n", i, i - 1, M->a[i], M->b[i - 1], errorestimate);
    errorestimate *= M->b[i] * M->b[i] / (M->b[i] * M->b[i] + 4 * Gamma * Gamma);
    while (OperatorPsiRestrictedRC(H, v[vi], &(v[vip]), Restrictions))
    {                           /*Allocation of Memory failed */
      printf("Allocation of new memory failed in LanczosTriDiagonalize for the function OperatorPsi\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
             epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      if (v[vi].NBlok == 1)
        AllocFailed();
      ComplexWaveFunctionTruncateBasis(&(v[vi]), epsilon / errorestimate);
    }
    printf("BasisSetSize: %12i ", v[vip].N);
    ComplexWaveFunctionDotProduct(&v[vip], &v[vi], &M->a[i + 1], &junk);
    while (ComplexWaveFunctionAddMultipliedUnionBasis(&(v[vip]), v[vim], -M->b[i]))
    {                           /*Allocation of Memory failed */
      printf
      ("Allocation of new memory failed in LanczosTriDiagonalize for the function RealWaveFunctionAddMultipliedUnionBasis(&(v[vip]),v[vim],-M->b[i])\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
             epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      temp = ComplexWaveFunctionCalculateNorm(v[vip]);
      if (v[vip].NBlok == 1)
        AllocFailed();
      ComplexWaveFunctionTruncateBasis(&(v[vip]), epsilon * temp / errorestimate);
    }
    ComplexWaveFunctionFree(&v[vim]);
  }
  printf("             | M.a[%3i] M.b[%3i]: %12.5E %12.5E | Err: %12.5E\n\n", i, i - 1, M->a[i], M->b[i - 1],
         errorestimate);
  if (M->N > 1)
  {
    ComplexWaveFunctionFree(&v[vi]);
    ComplexWaveFunctionFree(&v[vip]);
  }
  fflush(stdout);
  return 0;
}

uint32_t LanczosTriDiagonalizeRestrictedMC(const OperatorType H, WaveFunctionType psi, DeterminantRestrictionType const Restrictions, double epsilon, TriDiagonalMatrixType * const M, double Gamma)
{
  printf("\nStart of LanczosTriDiagonalizeRestrictedMC\n\n");
  fflush(stdout);
  WaveFunctionType v[3];
  uint32_t i, vi=0, vip=0, vim=0;
  double temp, errorestimate = 1.0, junk;
  /* copy the name of psi to v[0] in order to do the loop more cleanely  */
  v[0].N = psi.N;
  v[0].NBlok = psi.NBlok;
  v[0].NFermion = psi.NFermion;
  v[0].NBoson = psi.NBoson;
  v[0].NWord = psi.NWord;
  v[0].NBitsKey = psi.NBitsKey;
  v[0].ValueRe = psi.ValueRe;
  v[0].ValueIm = psi.ValueIm;
  v[0].Determinant = psi.Determinant;
  v[0].Next = psi.Next;
  v[0].KeyList = psi.KeyList;
  v[0].QComplex = psi.QComplex;
  /* finished copying */
  M->a[0] = 0;
  M->b[0] = ComplexWaveFunctionCalculateNorm(v[0]);
  ComplexWaveFunctionMultiply(v[0], (double)1.0 / M->b[0]);
  printf("BasisSetSize: %12i ", v[0].N);
  ComplexWaveFunctionTruncateBasis(&(v[0]), epsilon);
  printf("%12i\n", v[0].N);
  while (OperatorPsiRestrictedMC(H, v[0], &(v[1]), Restrictions))
  {                             /* Allocation of Memory failed */
    printf("Allocation of new memory failed in LanczosTriDiagonalize for the function OperatorPsi\n");
    printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
           epsilon * 2.0);
    fflush(stdout);
    epsilon *= 2.0;
    if (v[0].NBlok == 1)
      AllocFailed();
    ComplexWaveFunctionTruncateBasis(&(v[0]), epsilon);
  }
  ComplexWaveFunctionDotProduct(&v[1], &v[0], &M->a[1], &junk);
  printf("BasisSetSize: %12i ", v[1].N);
  for (i = 1; i < M->N; i++)
  {
    vip = (i + 1) % 3;
    vi = (i) % 3;
    vim = (i - 1) % 3;
    while (ComplexWaveFunctionAddMultipliedUnionBasis(&(v[vi]), v[vim], -M->a[i]))
    {                           /*Allocation of Memory failed */
      printf
      ("Allocation of new memory failed in LanczosTriDiagonalize for the function RealWaveFunctionAddMultipliedUnionBasis(&(v[vi]),v[vim],-M->a[i])\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
             epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      temp = ComplexWaveFunctionCalculateNorm(v[vi]);
      if (v[vi].NBlok == 1)
        AllocFailed();
      ComplexWaveFunctionTruncateBasis(&(v[vi]), epsilon * temp / errorestimate);
    }
    M->b[i] = ComplexWaveFunctionCalculateNorm(v[vi]);
    ComplexWaveFunctionMultiply(v[vi], (double)1.0 / M->b[i]);
    ComplexWaveFunctionTruncateBasis(&(v[vi]), epsilon / errorestimate);
    printf("%12i ", v[vi].N);
    printf("| M.a[%3i] M.b[%3i]: %12.5E %12.5E | Err: %12.5E\n", i, i - 1, M->a[i], M->b[i - 1], errorestimate);
    errorestimate *= M->b[i] * M->b[i] / (M->b[i] * M->b[i] + 4 * Gamma * Gamma);
    while (OperatorPsiRestrictedMC(H, v[vi], &(v[vip]), Restrictions))
    {                           /*Allocation of Memory failed */
      printf("Allocation of new memory failed in LanczosTriDiagonalize for the function OperatorPsi\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
             epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      if (v[vi].NBlok == 1)
        AllocFailed();
      ComplexWaveFunctionTruncateBasis(&(v[vi]), epsilon / errorestimate);
    }
    printf("BasisSetSize: %12i ", v[vip].N);
    ComplexWaveFunctionDotProduct(&v[vip], &v[vi], &M->a[i + 1], &junk);
    while (ComplexWaveFunctionAddMultipliedUnionBasis(&(v[vip]), v[vim], -M->b[i]))
    {                           /*Allocation of Memory failed */
      printf
      ("Allocation of new memory failed in LanczosTriDiagonalize for the function RealWaveFunctionAddMultipliedUnionBasis(&(v[vip]),v[vim],-M->b[i])\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
             epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      temp = ComplexWaveFunctionCalculateNorm(v[vip]);
      if (v[vip].NBlok == 1)
        AllocFailed();
      ComplexWaveFunctionTruncateBasis(&(v[vip]), epsilon * temp / errorestimate);
    }
    ComplexWaveFunctionFree(&v[vim]);
  }
  printf("             | M.a[%3i] M.b[%3i]: %12.5E %12.5E | Err: %12.5E\n\n", i, i - 1, M->a[i], M->b[i - 1],
         errorestimate);
  if (M->N > 1)
  {
    ComplexWaveFunctionFree(&v[vi]);
    ComplexWaveFunctionFree(&v[vip]);
  }
  fflush(stdout);
  return 0;
}

uint32_t LanczosTriDiagonalizeRestrictedMR(const OperatorType H, WaveFunctionType psi, DeterminantRestrictionType const Restrictions, double epsilon, TriDiagonalMatrixType * const M, double Gamma)
{
  uint32_t retval;
  WaveFunctionType Cpsi;
  if (MakeComplexWaveFunctionFromReal(&Cpsi, psi))
  {
    printf
    ("LanczosTriDiagonalizeRestrictedMR called with a real wave-function and a mixed complex Hamiltonian.\nTurned the wave-function complex in step 1.\nBut failed.\n");
    fflush(stdout);
    return 1;
  }
  retval = LanczosTriDiagonalizeRestrictedMC(H, Cpsi, Restrictions, epsilon, M, Gamma);
  MakeRealWaveFunctionFromComplex(Cpsi, &Cpsi);
  return retval;
}










/*  Different versions of LanczosTriDiagonalizePlusDensity depending on QComplex of H and psi
 *
 */
uint32_t LanczosTriDiagonalizePlusDensityRR(const OperatorType H, WaveFunctionType psi, double epsilon,
                                            TriDiagonalMatrixType * const M, MatrixType * const Rho, double Gamma)
{
  printf("\nStart of LanczosTriDiagonalizePlusDensityRR\n\n");
  fflush(stdout);
  WaveFunctionType v[3];
  uint32_t i, vi=0, vip=0, vim=0;
  double temp, errorestimate = 1.0;
  double QRe[3], QIm[3], PRe[3], PIm[3], ValEst[3];
  QRe[0] = 1;
  QRe[1] = 0;
  QRe[2] = 0;
  QIm[0] = 0;
  QIm[1] = 0;
  QIm[2] = 0;
  PRe[0] = 0;
  PRe[1] = 0;
  PRe[2] = 1;
  PIm[0] = 0;
  PIm[1] = 0;
  PIm[2] = 0;
  ValEst[0] = 0.0;
  ValEst[1] = 0.0;
  ValEst[2] = 0.0;
  /* copy the name of psi to v[0] in order to do the loop more cleanely  */
  v[0].N = psi.N;
  v[0].NBlok = psi.NBlok;
  v[0].NFermion = psi.NFermion;
  v[0].NBoson = psi.NBoson;
  v[0].NWord = psi.NWord;
  v[0].NBitsKey = psi.NBitsKey;
  v[0].ValueRe = psi.ValueRe;
  v[0].Determinant = psi.Determinant;
  v[0].Next = psi.Next;
  v[0].KeyList = psi.KeyList;
  v[0].QComplex = psi.QComplex;
  /* finished copying */
  M->a[0] = 0;
  M->b[0] = RealWaveFunctionCalculateNorm(v[0]);
  if (M->b[0] == 0)
  {
    RealWaveFunctionMultiply(v[0], (double)0.0);
  }
  else
  {
    RealWaveFunctionMultiply(v[0], (double)1.0 / M->b[0]);
  }
  CalculateDensity(Rho->VRe[0], v[0]);
  printf("BasisSetSize: %12i ", v[0].N);
  RealWaveFunctionTruncateBasis(&(v[0]), epsilon);
  printf("%12i\n", v[0].N);
  while (OperatorPsiRR(H, v[0], &(v[1])))
  {                             /* Allocation of Memory failed */
    printf("Allocation of new memory failed in LanczosTriDiagonalize for the function OperatorPsi\n");
    printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
           epsilon * 2.0);
    fflush(stdout);
    epsilon *= 2.0;
    if (v[0].NBlok == 1)
      AllocFailed();
    RealWaveFunctionTruncateBasis(&(v[0]), epsilon);
  }
  M->a[1] = RealWaveFunctionDotProduct(&v[1], &v[0]);
  printf("BasisSetSize: %12i ", v[1].N);
  for (i = 1; i < M->N; i++)
  {
    vip = (i + 1) % 3;
    vi = (i) % 3;
    vim = (i - 1) % 3;
    while (RealWaveFunctionAddMultipliedUnionBasis(&(v[vi]), v[vim], -M->a[i]))
    {                           /*Allocation of Memory failed */
      printf
        ("Allocation of new memory failed in LanczosTriDiagonalize for the function RealWaveFunctionAddMultipliedUnionBasis(&(v[vi]),v[vim],-M->a[i])\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
             epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      temp = RealWaveFunctionCalculateNorm(v[vi]);
      if (v[vi].NBlok == 1)
        AllocFailed();
      RealWaveFunctionTruncateBasis(&(v[vi]), epsilon * temp / errorestimate);
    }
    M->b[i] = RealWaveFunctionCalculateNorm(v[vi]);
    if (M->b[i] == 0)
    {
      RealWaveFunctionMultiply(v[vi], (double)0.0);
    }
    else
    {
      RealWaveFunctionMultiply(v[vi], (double)1.0 / M->b[i]);
    }
    CalculateDensity(Rho->VRe[i], v[vi]);
    RealWaveFunctionTruncateBasis(&(v[vi]), epsilon / errorestimate);
    printf("%12i ", v[vi].N);
    printf("| M.a[%3i] M.b[%3i]: %12.5E %12.5E | Err: %12.5E\n", i, i - 1, M->a[i], M->b[i - 1], errorestimate);
    /*
       QRe[vi]=M->b[i-1]*QRe[vip]-0.5*Gamma*QIm[vim];
       QIm[vi]=M->b[i-1]*QIm[vip]+0.5*Gamma*QRe[vim];
       PRe[vi]=M->b[i-1]*PRe[vip]-0.5*Gamma*PIm[vim];
       PIm[vi]=M->b[i-1]*PIm[vip]+0.5*Gamma*PRe[vim];
       ValEst[vi]=fabs((PIm[vi]*QRe[vi]-PRe[vi]*QIm[vi])/(QIm[vi]*QIm[vi]+QRe[vi]*QRe[vi]));
       if(i>5)errorestimate=fabs(2*(ValEst[vi]-ValEst[vim])/(Gamma));
     */
    errorestimate *= M->b[i] * M->b[i] / (M->b[i] * M->b[i] + 4 * Gamma * Gamma);
    while (OperatorPsiRR(H, v[vi], &(v[vip])))
    {                           /*Allocation of Memory failed */
      printf("Allocation of new memory failed in LanczosTriDiagonalize for the function OperatorPsi\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
             epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      if (v[vi].NBlok == 1)
        AllocFailed();
      RealWaveFunctionTruncateBasis(&(v[vi]), epsilon / errorestimate);
    }
    printf("BasisSetSize: %12i ", v[vip].N);
    M->a[i + 1] = RealWaveFunctionDotProduct(&v[vip], &v[vi]);
    while (RealWaveFunctionAddMultipliedUnionBasis(&(v[vip]), v[vim], -M->b[i]))
    {                           /*Allocation of Memory failed */
      printf
        ("Allocation of new memory failed in LanczosTriDiagonalize for the function RealWaveFunctionAddMultipliedUnionBasis(&(v[vip]),v[vim],-M->b[i])\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
             epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      temp = RealWaveFunctionCalculateNorm(v[vip]);
      if (v[vip].NBlok == 1)
        AllocFailed();
      RealWaveFunctionTruncateBasis(&(v[vip]), epsilon * temp / errorestimate);
    }
    RealWaveFunctionFree(&v[vim]);
  }
  printf("             | M.a[%3i] M.b[%3i]: %12.5E %12.5E | Err: %12.5E\n\n", i, i - 1, M->a[i], M->b[i - 1],
         errorestimate);
  if (M->N > 1)
  {
    RealWaveFunctionFree(&v[vi]);
    RealWaveFunctionFree(&v[vip]);
  }
  fflush(stdout);
  return 0;
}

uint32_t LanczosTriDiagonalizePlusDensityCR(const OperatorType H, WaveFunctionType psi, double epsilon,
                                            TriDiagonalMatrixType * const M, MatrixType * const Rho, double Gamma)
{
  uint32_t retval;
  WaveFunctionType Cpsi;
  if (MakeComplexWaveFunctionFromReal(&Cpsi, psi))
  {
    printf
      ("LanczosTriDiagonalizePlusDensityCR called with a real wave-function and a complex Hamiltonian.\nTurned the wave-function complex in step 1.\nBut failed.\n");
    return 1;
  }
  retval = LanczosTriDiagonalizePlusDensityCC(H, Cpsi, epsilon, M, Rho, Gamma);
  MakeRealWaveFunctionFromComplex(Cpsi, &Cpsi);
  return retval;
}

uint32_t LanczosTriDiagonalizePlusDensityMR(const OperatorType H, WaveFunctionType psi, double epsilon,
                                            TriDiagonalMatrixType * const M, MatrixType * const Rho, double Gamma)
{
  uint32_t retval;
  WaveFunctionType Cpsi;
  if (MakeComplexWaveFunctionFromReal(&Cpsi, psi))
  {
    printf
      ("LanczosTriDiagonalizePlusDensityMR called with a real wave-function and a mixed complex Hamiltonian.\nTurned the wave-function complex in step 1.\nBut failed.\n");
    return 1;
  }
  retval = LanczosTriDiagonalizePlusDensityMC(H, Cpsi, epsilon, M, Rho, Gamma);
  MakeRealWaveFunctionFromComplex(Cpsi, &Cpsi);
  return retval;
}

uint32_t LanczosTriDiagonalizePlusDensityRC(const OperatorType H, WaveFunctionType psi, double epsilon,
                                            TriDiagonalMatrixType * const M, MatrixType * const Rho, double Gamma)
{
  printf("\nStart of LanczosTriDiagonalizePlusDensityRC\n\n");
  fflush(stdout);
  WaveFunctionType v[3];
  uint32_t i, vi=0, vip=0, vim=0;
  double temp, errorestimate = 1.0, junk;
  /* copy the name of psi to v[0] in order to do the loop more cleanely  */
  v[0].N = psi.N;
  v[0].NBlok = psi.NBlok;
  v[0].NFermion = psi.NFermion;
  v[0].NBoson = psi.NBoson;
  v[0].NWord = psi.NWord;
  v[0].NBitsKey = psi.NBitsKey;
  v[0].ValueRe = psi.ValueRe;
  v[0].ValueIm = psi.ValueIm;
  v[0].Determinant = psi.Determinant;
  v[0].Next = psi.Next;
  v[0].KeyList = psi.KeyList;
  v[0].QComplex = psi.QComplex;
  /* finished copying */
  M->a[0] = 0;
  M->b[0] = ComplexWaveFunctionCalculateNorm(v[0]);
  ComplexWaveFunctionMultiply(v[0], (double)1.0 / M->b[0]);
  printf("BasisSetSize: %12i ", v[0].N);
  ComplexWaveFunctionTruncateBasis(&(v[0]), epsilon);
  printf("%12i\n", v[0].N);
  while (OperatorPsiRC(H, v[0], &(v[1])))
  {                             /* Allocation of Memory failed */
    printf("Allocation of new memory failed in LanczosTriDiagonalize for the function OperatorPsi\n");
    printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
           epsilon * 2.0);
    fflush(stdout);
    epsilon *= 2.0;
    if (v[0].NBlok == 1)
      AllocFailed();
    ComplexWaveFunctionTruncateBasis(&(v[0]), epsilon);
  }
  ComplexWaveFunctionDotProduct(&v[1], &v[0], &M->a[1], &junk);
  printf("BasisSetSize: %12i ", v[1].N);
  for (i = 1; i < M->N; i++)
  {
    vip = (i + 1) % 3;
    vi = (i) % 3;
    vim = (i - 1) % 3;
    while (ComplexWaveFunctionAddMultipliedUnionBasis(&(v[vi]), v[vim], -M->a[i]))
    {                           /*Allocation of Memory failed */
      printf
        ("Allocation of new memory failed in LanczosTriDiagonalize for the function RealWaveFunctionAddMultipliedUnionBasis(&(v[vi]),v[vim],-M->a[i])\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
             epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      temp = ComplexWaveFunctionCalculateNorm(v[vi]);
      if (v[vi].NBlok == 1)
        AllocFailed();
      ComplexWaveFunctionTruncateBasis(&(v[vi]), epsilon * temp / errorestimate);
    }
    M->b[i] = ComplexWaveFunctionCalculateNorm(v[vi]);
    ComplexWaveFunctionMultiply(v[vi], (double)1.0 / M->b[i]);
    ComplexWaveFunctionTruncateBasis(&(v[vi]), epsilon / errorestimate);
    printf("%12i ", v[vi].N);
    printf("| M.a[%3i] M.b[%3i]: %12.5E %12.5E | Err: %12.5E\n", i, i - 1, M->a[i], M->b[i - 1], errorestimate);
    errorestimate *= M->b[i] * M->b[i] / (M->b[i] * M->b[i] + 4 * Gamma * Gamma);
    while (OperatorPsiRC(H, v[vi], &(v[vip])))
    {                           /*Allocation of Memory failed */
      printf("Allocation of new memory failed in LanczosTriDiagonalize for the function OperatorPsi\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
             epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      if (v[vi].NBlok == 1)
        AllocFailed();
      ComplexWaveFunctionTruncateBasis(&(v[vi]), epsilon / errorestimate);
    }
    printf("BasisSetSize: %12i ", v[vip].N);
    ComplexWaveFunctionDotProduct(&v[vip], &v[vi], &M->a[i + 1], &junk);
    while (ComplexWaveFunctionAddMultipliedUnionBasis(&(v[vip]), v[vim], -M->b[i]))
    {                           /*Allocation of Memory failed */
      printf
        ("Allocation of new memory failed in LanczosTriDiagonalize for the function RealWaveFunctionAddMultipliedUnionBasis(&(v[vip]),v[vim],-M->b[i])\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
             epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      temp = ComplexWaveFunctionCalculateNorm(v[vip]);
      if (v[vip].NBlok == 1)
        AllocFailed();
      ComplexWaveFunctionTruncateBasis(&(v[vip]), epsilon * temp / errorestimate);
    }
    ComplexWaveFunctionFree(&v[vim]);
  }
  printf("             | M.a[%3i] M.b[%3i]: %12.5E %12.5E | Err: %12.5E\n\n", i, i - 1, M->a[i], M->b[i - 1],
         errorestimate);
  if (M->N > 1)
  {
    ComplexWaveFunctionFree(&v[vi]);
    ComplexWaveFunctionFree(&v[vip]);
  }
  fflush(stdout);
  return 0;
}

uint32_t LanczosTriDiagonalizePlusDensityCC(const OperatorType H, WaveFunctionType psi, double epsilon,
                                            TriDiagonalMatrixType * const M, MatrixType * const Rho, double Gamma)
{
  printf("\nStart of LanczosTriDiagonalizePlusDensityCC\n\n");
  fflush(stdout);
  WaveFunctionType v[3];
  uint32_t i, vi=0, vip=0, vim=0;
  double temp, errorestimate = 1.0, junk;
  /* copy the name of psi to v[0] in order to do the loop more cleanely  */
  v[0].N = psi.N;
  v[0].NBlok = psi.NBlok;
  v[0].NFermion = psi.NFermion;
  v[0].NBoson = psi.NBoson;
  v[0].NWord = psi.NWord;
  v[0].NBitsKey = psi.NBitsKey;
  v[0].ValueRe = psi.ValueRe;
  v[0].ValueIm = psi.ValueIm;
  v[0].Determinant = psi.Determinant;
  v[0].Next = psi.Next;
  v[0].KeyList = psi.KeyList;
  v[0].QComplex = psi.QComplex;
  /* finished copying */
  M->a[0] = 0;
  M->b[0] = ComplexWaveFunctionCalculateNorm(v[0]);
  ComplexWaveFunctionMultiply(v[0], (double)1.0 / M->b[0]);
  printf("BasisSetSize: %12i ", v[0].N);
  ComplexWaveFunctionTruncateBasis(&(v[0]), epsilon);
  printf("%12i\n", v[0].N);
  while (OperatorPsiCC(H, v[0], &(v[1])))
  {                             /* Allocation of Memory failed */
    printf("Allocation of new memory failed in LanczosTriDiagonalize for the function OperatorPsi\n");
    printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
           epsilon * 2.0);
    fflush(stdout);
    epsilon *= 2.0;
    if (v[0].NBlok == 1)
      AllocFailed();
    ComplexWaveFunctionTruncateBasis(&(v[0]), epsilon);
  }
  ComplexWaveFunctionDotProduct(&v[1], &v[0], &M->a[1], &junk);
  printf("BasisSetSize: %12i ", v[1].N);
  for (i = 1; i < M->N; i++)
  {
    vip = (i + 1) % 3;
    vi = (i) % 3;
    vim = (i - 1) % 3;
    while (ComplexWaveFunctionAddMultipliedUnionBasis(&(v[vi]), v[vim], -M->a[i]))
    {                           /*Allocation of Memory failed */
      printf
        ("Allocation of new memory failed in LanczosTriDiagonalize for the function RealWaveFunctionAddMultipliedUnionBasis(&(v[vi]),v[vim],-M->a[i])\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
             epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      temp = ComplexWaveFunctionCalculateNorm(v[vi]);
      if (v[vi].NBlok == 1)
        AllocFailed();
      ComplexWaveFunctionTruncateBasis(&(v[vi]), epsilon * temp / errorestimate);
    }
    M->b[i] = ComplexWaveFunctionCalculateNorm(v[vi]);
    ComplexWaveFunctionMultiply(v[vi], (double)1.0 / M->b[i]);
    ComplexWaveFunctionTruncateBasis(&(v[vi]), epsilon / errorestimate);
    printf("%12i ", v[vi].N);
    printf("| M.a[%3i] M.b[%3i]: %12.5E %12.5E | Err: %12.5E\n", i, i - 1, M->a[i], M->b[i - 1], errorestimate);
    errorestimate *= M->b[i] * M->b[i] / (M->b[i] * M->b[i] + 4 * Gamma * Gamma);
    while (OperatorPsiCC(H, v[vi], &(v[vip])))
    {                           /*Allocation of Memory failed */
      printf("Allocation of new memory failed in LanczosTriDiagonalize for the function OperatorPsi\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
             epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      if (v[vi].NBlok == 1)
        AllocFailed();
      ComplexWaveFunctionTruncateBasis(&(v[vi]), epsilon / errorestimate);
    }
    printf("BasisSetSize: %12i ", v[vip].N);
    ComplexWaveFunctionDotProduct(&v[vip], &v[vi], &M->a[i + 1], &junk);
    while (ComplexWaveFunctionAddMultipliedUnionBasis(&(v[vip]), v[vim], -M->b[i]))
    {                           /*Allocation of Memory failed */
      printf
        ("Allocation of new memory failed in LanczosTriDiagonalize for the function RealWaveFunctionAddMultipliedUnionBasis(&(v[vip]),v[vim],-M->b[i])\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
             epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      temp = ComplexWaveFunctionCalculateNorm(v[vip]);
      if (v[vip].NBlok == 1)
        AllocFailed();
      ComplexWaveFunctionTruncateBasis(&(v[vip]), epsilon * temp / errorestimate);
    }
    ComplexWaveFunctionFree(&v[vim]);
  }
  printf("             | M.a[%3i] M.b[%3i]: %12.5E %12.5E | Err: %12.5E\n\n", i, i - 1, M->a[i], M->b[i - 1],
         errorestimate);
  if (M->N > 1)
  {
    ComplexWaveFunctionFree(&v[vi]);
    ComplexWaveFunctionFree(&v[vip]);
  }
  fflush(stdout);
  return 0;
}

uint32_t LanczosTriDiagonalizePlusDensityMC(const OperatorType H, WaveFunctionType psi, double epsilon,
                                            TriDiagonalMatrixType * const M, MatrixType * const Rho, double Gamma)
{
  printf("\nStart of LanczosTriDiagonalizePlusDensityMC\n\n");
  fflush(stdout);
  WaveFunctionType v[3];
  uint32_t i, vi=0, vip=0, vim=0;
  double temp, errorestimate = 1.0, junk;
  /* copy the name of psi to v[0] in order to do the loop more cleanely  */
  v[0].N = psi.N;
  v[0].NBlok = psi.NBlok;
  v[0].NFermion = psi.NFermion;
  v[0].NBoson = psi.NBoson;
  v[0].NWord = psi.NWord;
  v[0].NBitsKey = psi.NBitsKey;
  v[0].ValueRe = psi.ValueRe;
  v[0].ValueIm = psi.ValueIm;
  v[0].Determinant = psi.Determinant;
  v[0].Next = psi.Next;
  v[0].KeyList = psi.KeyList;
  v[0].QComplex = psi.QComplex;
  /* finished copying */
  M->a[0] = 0;
  M->b[0] = ComplexWaveFunctionCalculateNorm(v[0]);
  ComplexWaveFunctionMultiply(v[0], (double)1.0 / M->b[0]);
  printf("BasisSetSize: %12i ", v[0].N);
  ComplexWaveFunctionTruncateBasis(&(v[0]), epsilon);
  printf("%12i\n", v[0].N);
  while (OperatorPsiMC(H, v[0], &(v[1])))
  {                             /* Allocation of Memory failed */
    printf("Allocation of new memory failed in LanczosTriDiagonalize for the function OperatorPsi\n");
    printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
           epsilon * 2.0);
    fflush(stdout);
    epsilon *= 2.0;
    if (v[0].NBlok == 1)
      AllocFailed();
    ComplexWaveFunctionTruncateBasis(&(v[0]), epsilon);
  }
  ComplexWaveFunctionDotProduct(&v[1], &v[0], &M->a[1], &junk);
  printf("BasisSetSize: %12i ", v[1].N);
  for (i = 1; i < M->N; i++)
  {
    vip = (i + 1) % 3;
    vi = (i) % 3;
    vim = (i - 1) % 3;
    while (ComplexWaveFunctionAddMultipliedUnionBasis(&(v[vi]), v[vim], -M->a[i]))
    {                           /*Allocation of Memory failed */
      printf
        ("Allocation of new memory failed in LanczosTriDiagonalize for the function RealWaveFunctionAddMultipliedUnionBasis(&(v[vi]),v[vim],-M->a[i])\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
             epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      temp = ComplexWaveFunctionCalculateNorm(v[vi]);
      if (v[vi].NBlok == 1)
        AllocFailed();
      ComplexWaveFunctionTruncateBasis(&(v[vi]), epsilon * temp / errorestimate);
    }
    M->b[i] = ComplexWaveFunctionCalculateNorm(v[vi]);
    ComplexWaveFunctionMultiply(v[vi], (double)1.0 / M->b[i]);
    ComplexWaveFunctionTruncateBasis(&(v[vi]), epsilon / errorestimate);
    printf("%12i ", v[vi].N);
    printf("| M.a[%3i] M.b[%3i]: %12.5E %12.5E | Err: %12.5E\n", i, i - 1, M->a[i], M->b[i - 1], errorestimate);
    errorestimate *= M->b[i] * M->b[i] / (M->b[i] * M->b[i] + 4 * Gamma * Gamma);
    while (OperatorPsiMC(H, v[vi], &(v[vip])))
    {                           /*Allocation of Memory failed */
      printf("Allocation of new memory failed in LanczosTriDiagonalize for the function OperatorPsi\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
             epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      if (v[vi].NBlok == 1)
        AllocFailed();
      ComplexWaveFunctionTruncateBasis(&(v[vi]), epsilon / errorestimate);
    }
    printf("BasisSetSize: %12i ", v[vip].N);
    ComplexWaveFunctionDotProduct(&v[vip], &v[vi], &M->a[i + 1], &junk);
    while (ComplexWaveFunctionAddMultipliedUnionBasis(&(v[vip]), v[vim], -M->b[i]))
    {                           /*Allocation of Memory failed */
      printf
        ("Allocation of new memory failed in LanczosTriDiagonalize for the function RealWaveFunctionAddMultipliedUnionBasis(&(v[vip]),v[vim],-M->b[i])\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
             epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      temp = ComplexWaveFunctionCalculateNorm(v[vip]);
      if (v[vip].NBlok == 1)
        AllocFailed();
      ComplexWaveFunctionTruncateBasis(&(v[vip]), epsilon * temp / errorestimate);
    }
    ComplexWaveFunctionFree(&v[vim]);
  }
  printf("             | M.a[%3i] M.b[%3i]: %12.5E %12.5E | Err: %12.5E\n\n", i, i - 1, M->a[i], M->b[i - 1],
         errorestimate);
  if (M->N > 1)
  {
    ComplexWaveFunctionFree(&v[vi]);
    ComplexWaveFunctionFree(&v[vip]);
  }
  fflush(stdout);
  return 0;
}


uint32_t LanczosTriDiagonalizeKrylovRRNonOrtho(const OperatorType H, WaveFunctionType psi, double epsilon,
                                               TriDiagonalMatrixType * const M, WaveFunctionType * PsiBasis,
                                               double Gamma)
{
  printf("\nStart of LanczosTriDiagonalizeKrylovNonOrthoRR\n\n");
  fflush(stdout);
  uint32_t i;
  double temp, errorestimate = 1.0;
  /* Allocate memory for *PsiBasis */
  if ((PsiBasis = calloc(M->N + 1, sizeof(WaveFunctionType))) == NULL)
  {
    printf("calloc failed in LanczosTriDiagonalizeKrylovRR 00\n");
    fflush(stdout);
    return 1;
  };
  /* copy the name of psi to PsiBasis[0] in order to do the loop more cleanely  */
  PsiBasis[0].N = psi.N;
  PsiBasis[0].NBlok = psi.NBlok;
  PsiBasis[0].NFermion = psi.NFermion;
  PsiBasis[0].NBoson = psi.NBoson;
  PsiBasis[0].NWord = psi.NWord;
  PsiBasis[0].NBitsKey = psi.NBitsKey;
  PsiBasis[0].ValueRe = psi.ValueRe;
  PsiBasis[0].Determinant = psi.Determinant;
  PsiBasis[0].Next = psi.Next;
  PsiBasis[0].KeyList = psi.KeyList;
  PsiBasis[0].QComplex = psi.QComplex;
  /* calculate a[0] and b[0]  */
  M->a[0] = 0;
  M->b[0] = RealWaveFunctionCalculateNorm(PsiBasis[0]);
  /* Normalize PsiBasis[0] */
  if (M->b[0] == 0)
  {
    RealWaveFunctionMultiply(PsiBasis[0], (double)0.0);
  }
  else
  {
    RealWaveFunctionMultiply(PsiBasis[0], (double)1.0 / M->b[0]);
  }
  printf("BasisSetSize: %12i ", PsiBasis[0].N);
  RealWaveFunctionTruncateBasis(&(PsiBasis[0]), epsilon);
  printf("%12i\n", PsiBasis[0].N);
  /* Create PsiBasis[1] */
  while (OperatorPsiRR(H, PsiBasis[0], &(PsiBasis[1])))
  {                             /* Allocation of Memory failed */
    printf("Allocation of new memory failed in LanczosTriDiagonalizeKrylovRR for the function OperatorPsi\n");
    printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
           epsilon * 2.0);
    fflush(stdout);
    epsilon *= 2.0;
    if (PsiBasis[0].NBlok == 1)
      AllocFailed();
    RealWaveFunctionTruncateBasis(&(PsiBasis[0]), epsilon);
  }
  M->a[1] = RealWaveFunctionDotProduct(&PsiBasis[1], &PsiBasis[0]);
  printf("BasisSetSize: %12i ", PsiBasis[1].N);
  /* Create PsiBasis[2] upto PsiBasis[M->N-1] */
  for (i = 1; i < M->N; i++)
  {
    while (RealWaveFunctionAddMultipliedUnionBasis(&(PsiBasis[i]), PsiBasis[i - 1], -M->a[i]))
    {                           /*Allocation of Memory failed */
      printf
        ("Allocation of new memory failed in LanczosTriDiagonalizeKrylovRR for the function RealWaveFunctionAddMultipliedUnionBasis(&(v[vi]),v[vim],-M->a[i])\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
             epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      temp = RealWaveFunctionCalculateNorm(PsiBasis[i]);
      if (PsiBasis[i].NBlok == 1)
        AllocFailed();
      RealWaveFunctionTruncateBasis(&(PsiBasis[i]), epsilon * temp / errorestimate);
    }
    M->b[i] = RealWaveFunctionCalculateNorm(PsiBasis[i]);
    if (M->b[i] == 0)
    {
      RealWaveFunctionMultiply(PsiBasis[i], (double)0.0);
    }
    else
    {
      RealWaveFunctionMultiply(PsiBasis[i], (double)1.0 / M->b[i]);
    }
    RealWaveFunctionTruncateBasis(&(PsiBasis[i]), epsilon);
    printf("%12i ", PsiBasis[i].N);
    printf("| M.a[%3i] M.b[%3i]: %12.5E %12.5E | Err: %12.5E\n", i, i - 1, M->a[i], M->b[i - 1], errorestimate);
    errorestimate *= M->b[i] * M->b[i] / (M->b[i] * M->b[i] + 4 * Gamma * Gamma);
    while (OperatorPsiRR(H, PsiBasis[i], &(PsiBasis[i + 1])))
    {                           /*Allocation of Memory failed */
      printf("Allocation of new memory failed in LanczosTriDiagonalize for the function OperatorPsi\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
             epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      if (PsiBasis[i].NBlok == 1)
        AllocFailed();
      RealWaveFunctionTruncateBasis(&(PsiBasis[i]), epsilon / errorestimate);
    }
    printf("BasisSetSize: %12i ", PsiBasis[i + 1].N);
    M->a[i + 1] = RealWaveFunctionDotProduct(&PsiBasis[i + 1], &PsiBasis[i]);
    while (RealWaveFunctionAddMultipliedUnionBasis(&(PsiBasis[i + 1]), PsiBasis[i - 1], -M->b[i]))
    {                           /*Allocation of Memory failed */
      printf
        ("Allocation of new memory failed in LanczosTriDiagonalizeKrylovRR for the function RealWaveFunctionAddMultipliedUnionBasis(&(v[vip]),v[vim],-M->b[i])\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
             epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      temp = RealWaveFunctionCalculateNorm(PsiBasis[i + 1]);
      if (PsiBasis[i + 1].NBlok == 1)
        AllocFailed();
      RealWaveFunctionTruncateBasis(&(PsiBasis[i + 1]), epsilon * temp / errorestimate);
    }
  }
  printf("             | M.a[%3i] M.b[%3i]: %12.5E %12.5E | Err: %12.5E\n\n", i, i - 1, M->a[i], M->b[i - 1],
         errorestimate);
  RealWaveFunctionFree(&PsiBasis[M->N]);
  fflush(stdout);
  return 0;
}


uint32_t LanczosTriDiagonalizeKrylovRR(const OperatorType H, WaveFunctionType psi, double epsilon,
                                       TriDiagonalMatrixType * const M, WaveFunctionType * const PsiBasis)
{
  printf("\nStart of LanczosTriDiagonalizeKrylovRR\n\n");
  fflush(stdout);
  uint32_t i, j, k, stop;
  double temp, dot;
  double errorestimate = 1;
  /* copy the name of psi to PsiBasis[0] in order to do the loop more cleanely  */
  PsiBasis[0].N = psi.N;
  PsiBasis[0].NBlok = psi.NBlok;
  PsiBasis[0].NFermion = psi.NFermion;
  PsiBasis[0].NBoson = psi.NBoson;
  PsiBasis[0].NWord = psi.NWord;
  PsiBasis[0].NBitsKey = psi.NBitsKey;
  PsiBasis[0].ValueRe = psi.ValueRe;
  PsiBasis[0].Determinant = psi.Determinant;
  PsiBasis[0].Next = psi.Next;
  PsiBasis[0].KeyList = psi.KeyList;
  PsiBasis[0].QComplex = psi.QComplex;
  /* calculate a[0] and b[0]  */
  M->a[0] = 0;
  M->b[0] = RealWaveFunctionCalculateNorm(PsiBasis[0]);
  /* Normalize PsiBasis[0] */
  if (M->b[0] == 0)
  {
    RealWaveFunctionMultiply(PsiBasis[0], (double)0.0);
  }
  else
  {
    RealWaveFunctionMultiply(PsiBasis[0], (double)1.0 / M->b[0]);
  }
  printf("BasisSetSize: %12i ", PsiBasis[0].N);
  RealWaveFunctionTruncateBasis(&(PsiBasis[0]), epsilon / errorestimate);
  printf("%12i ", PsiBasis[0].N);
  printf("| M.a[  0]         : %12.5E\n", M->a[0]);
  /* Create PsiBasis[1] upto PsiBasis[M-N-1] */
  for (i = 1; i < M->N; i++)
  {
    while (OperatorPsiRR(H, PsiBasis[i - 1], &(PsiBasis[i])))
    {                           /* Allocation of Memory failed */
      printf("Allocation of new memory failed in LanczosTriDiagonalizeKrylovRR for the function OperatorPsi\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
             epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      stop = 1;
      for (j = 0; j < i; j++)
      {
        if (PsiBasis[j].NBlok != 1)
        {
          RealWaveFunctionTruncateBasis(&(PsiBasis[j]), epsilon / errorestimate);
          stop = 0;
        }
      }
      if (stop == 1)
        AllocFailed();
    }
    for (j = i; j > 0;)
    {
      j--;
      dot = RealWaveFunctionDotProduct(&PsiBasis[j], &PsiBasis[i]);
      if (j == i - 1)
      {
        M->a[i] = dot;
      }
      while (RealWaveFunctionAddMultipliedUnionBasis(&(PsiBasis[i]), PsiBasis[j], -dot))
      {                         /* Allocation of Memory failed */
        printf
          ("Allocation of new memory failed in LanczosTriDiagonalizeKrylovRR for the function RealWaveFunctionAddMultipliedUnionBasis(&(PsiBasis[i]),PsiBasis[j],-dot)\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
               epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        temp = RealWaveFunctionCalculateNorm(PsiBasis[i]);
        stop = 1;
        if (PsiBasis[i].NBlok != 1)
        {
          RealWaveFunctionTruncateBasis(&(PsiBasis[i]), epsilon * temp / errorestimate);
          stop = 0;
        }
        for (k = 0; k < i; k++)
        {
          if (PsiBasis[k].NBlok != 1)
          {
            RealWaveFunctionTruncateBasis(&(PsiBasis[k]), epsilon / errorestimate);
            stop = 0;
          }
        }
        if (stop == 1)
          AllocFailed();
      }
    }
    dot = RealWaveFunctionCalculateNorm(PsiBasis[i]);
    RealWaveFunctionMultiply(PsiBasis[i], 1 / dot);
    M->b[i] = dot;
    errorestimate *= Sqr(M->b[i]) / (Sqr(M->b[i]) + DBL_MIN + DBL_EPSILON * Sqr(M->a[i]));
    printf("BasisSetSize: %12i ", PsiBasis[i].N);
    RealWaveFunctionTruncateBasis(&(PsiBasis[i]), epsilon / errorestimate);
    printf("%12i ", PsiBasis[i].N);
    printf("| M.a[%3i] M.b[%3i]: %12.5E %12.5E | Err: %12.5E\n", i, i - 1, M->a[i], M->b[i - 1], errorestimate);
  }
  PsiOperatorPsiRR(H, PsiBasis[M->N - 1], &(M->a[M->N]));
  printf("BasisSetSize: %12i ", PsiBasis[M->N - 1].N);
  printf("             | M.a[%3i] M.b[%3i]: %12.5E %12.5E | Err: %12.5E\n\n", M->N, M->N - 1, M->a[M->N],
         M->b[M->N - 1], errorestimate);
  fflush(stdout);
  return 0;
}

uint32_t LanczosTriDiagonalizeKrylovRC(const OperatorType H, WaveFunctionType psi, double epsilon,
                                       TriDiagonalMatrixType * const M, WaveFunctionType * const PsiBasis)
{
  printf("\nStart of LanczosTriDiagonalizeKrylovRC\n\n");
  printf("Still to do\n");
  fflush(stdout);
  return 0;
}

uint32_t LanczosTriDiagonalizeKrylovCC(const OperatorType H, WaveFunctionType psi, double epsilon,
                                       TriDiagonalMatrixType * const M, WaveFunctionType * const PsiBasis)
{
  printf("\nStart of LanczosTriDiagonalizeKrylovCC\n\n");
  printf("Still to do\n");
  fflush(stdout);
  return 0;
}

uint32_t LanczosTriDiagonalizeKrylovMC(const OperatorType H, WaveFunctionType psi, double epsilon,
                                       TriDiagonalMatrixType * const M, WaveFunctionType * const PsiBasis)
{
  printf("\nStart of LanczosTriDiagonalizeKrylovMC\n\n");
  printf("Still to do\n");
  fflush(stdout);
  return 0;
}


uint32_t ArnoldiTriDiagonalizeKrylovRR(const OperatorType H, WaveFunctionType psi, double epsilon, MatrixType * const M,
                                       WaveFunctionType * const PsiBasis)
{
  printf("\nStart of ArnoldiTriDiagonalizeKrylovRR\n\n");
  fflush(stdout);
  uint32_t i, j, stop;
  double temp;
  double errorestimate = 1;
  MatrixType lastcolumn;
  /* copy the name of psi to PsiBasis[0] in order to do the loop more cleanely  */
  PsiBasis[0].N = psi.N;
  PsiBasis[0].NBlok = psi.NBlok;
  PsiBasis[0].NFermion = psi.NFermion;
  PsiBasis[0].NBoson = psi.NBoson;
  PsiBasis[0].NWord = psi.NWord;
  PsiBasis[0].NBitsKey = psi.NBitsKey;
  PsiBasis[0].ValueRe = psi.ValueRe;
  PsiBasis[0].Determinant = psi.Determinant;
  PsiBasis[0].Next = psi.Next;
  PsiBasis[0].KeyList = psi.KeyList;
  PsiBasis[0].QComplex = psi.QComplex;
  /* calculate a[0] and b[0]  */
  M->a0 = 0;
  M->b0 = RealWaveFunctionCalculateNorm(PsiBasis[0]);
  /* Normalize PsiBasis[0] */
  if (M->b0 == 0)
  {
    RealWaveFunctionMultiply(PsiBasis[0], (double)0.0);
  }
  else
  {
    RealWaveFunctionMultiply(PsiBasis[0], (double)1.0 / M->b0);
  }
  printf("BasisSetSize: %12i ", PsiBasis[0].N);
  RealWaveFunctionTruncateBasis(&(PsiBasis[0]), epsilon);
  printf("%12i ", PsiBasis[0].N);
  printf("| M.a[  0] M.b[  0]: %12.5E %12.5E | Err: %12.5E\n", M->a0, M->b0, errorestimate);
  /* Create PsiBasis[1] upto PsiBasis[M-N-1] */
  for (i = 1; i < M->Ni; i++)
  {
    while (OperatorPsiRR(H, PsiBasis[i - 1], &(PsiBasis[i])))
    {                           /* Allocation of Memory failed */
      printf("Allocation of new memory failed in LanczosTriDiagonalizeKrylovRR for the function OperatorPsi\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
             epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      stop = 1;
      for (j = 0; j < i; j++)
      {
        if (PsiBasis[j].NBlok != 1)
        {
          RealWaveFunctionTruncateBasis(&(PsiBasis[j]), epsilon);
          stop = 0;
        }
      }
      if (stop == 1)
        AllocFailed();
    }
    for (j = i; j > 0;)
    {
      j--;
      M->VRe[j][i - 1] = RealWaveFunctionDotProduct(&PsiBasis[j], &PsiBasis[i]);
      while (RealWaveFunctionAddMultipliedUnionBasis(&(PsiBasis[i]), PsiBasis[j], -M->VRe[j][i - 1]))
      {                         /* Allocation of Memory failed */
        printf
          ("Allocation of new memory failed in LanczosTriDiagonalizeKrylovRR for the function RealWaveFunctionAddMultipliedUnionBasis(&(PsiBasis[i]),PsiBasis[j],-dot)\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
               epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        temp = RealWaveFunctionCalculateNorm(PsiBasis[i]);
        stop = 1;
        if (PsiBasis[i].NBlok != 1)
        {
          RealWaveFunctionTruncateBasis(&(PsiBasis[i]), epsilon * temp);
          stop = 0;
        }
        for (j = 0; j < i; j++)
        {
          if (PsiBasis[j].NBlok != 1)
          {
            RealWaveFunctionTruncateBasis(&(PsiBasis[j]), epsilon);
            stop = 0;
          }
        }
        if (stop == 1)
          AllocFailed();
      }
    }
    M->VRe[i][i - 1] = RealWaveFunctionCalculateNorm(PsiBasis[i]);
    RealWaveFunctionMultiply(PsiBasis[i], 1 / M->VRe[i][i - 1]);
    printf("BasisSetSize: %12i ", PsiBasis[i].N);
    RealWaveFunctionTruncateBasis(&(PsiBasis[i]), epsilon);
    printf("%12i ", PsiBasis[i].N);
    printf("| M.a[%3i] M.b[%3i]: %12.5E %12.5E | Err: %12.5E\n", i, i, M->VRe[i - 1][i - 1], M->VRe[i][i - 1],
           errorestimate);
  }
  /* Calculate column M->Ni */
  lastcolumn.Ni = M->Ni;
  lastcolumn.Nj = 1;
  lastcolumn.QComplex = 0;
  InitMatrix(&lastcolumn);
  PsiPrimeOperatorPsiMatrix(H, &(PsiBasis[0]), &(PsiBasis[M->Ni - 1]), &lastcolumn);
  for (i = 0; i < M->Ni; i++)
  {
    M->VRe[i][M->Nj - 1] = lastcolumn.VRe[i][0];
  }
  FreeMatrix(&lastcolumn);
  fflush(stdout);
  return 0;
}

uint32_t ArnoldiTriDiagonalizeKrylovRC(const OperatorType H, WaveFunctionType psi, double epsilon, MatrixType * const M,
                                       WaveFunctionType * const PsiBasis)
{
  printf("\nStart of ArnoldiTriDiagonalizeKrylovRC\n\n");
  fflush(stdout);
  uint32_t i, j, stop;
  double temp;
  double errorestimate = 1;
  MatrixType lastcolumn;
  /* test input */
  if (M->QComplex == 0)
  {
    printf
      ("ArnoldiTriDiagonalizeKrylovRC entered with real matrix for results/n Although this might be intentional the alghoritms assumes a complex matrix to store the results\n Matrix turned into complex matrix");
    MakeRealMatrixComplex(M);
  }
  /* copy the name of psi to PsiBasis[0] in order to do the loop more cleanely  */
  PsiBasis[0].N = psi.N;
  PsiBasis[0].NBlok = psi.NBlok;
  PsiBasis[0].NFermion = psi.NFermion;
  PsiBasis[0].NBoson = psi.NBoson;
  PsiBasis[0].NWord = psi.NWord;
  PsiBasis[0].NBitsKey = psi.NBitsKey;
  PsiBasis[0].ValueRe = psi.ValueRe;
  PsiBasis[0].ValueIm = psi.ValueIm;
  PsiBasis[0].Determinant = psi.Determinant;
  PsiBasis[0].Next = psi.Next;
  PsiBasis[0].KeyList = psi.KeyList;
  PsiBasis[0].QComplex = psi.QComplex;
  /* calculate a[0] and b[0]  */
  M->a0 = 0;
  M->b0 = ComplexWaveFunctionCalculateNorm(PsiBasis[0]);
  /* Normalize PsiBasis[0] */
  if (M->b0 == 0)
  {
    ComplexWaveFunctionMultiply(PsiBasis[0], (double)0.0);
  }
  else
  {
    ComplexWaveFunctionMultiply(PsiBasis[0], (double)1.0 / M->b0);
  }
  printf("BasisSetSize: %12i ", PsiBasis[0].N);
  ComplexWaveFunctionTruncateBasis(&(PsiBasis[0]), epsilon);
  printf("%12i ", PsiBasis[0].N);
  printf("| M.a[  0] M.b[  0]: %12.5E %12.5E | Err: %12.5E\n", M->a0, M->b0, errorestimate);
  /* Create PsiBasis[1] upto PsiBasis[M-N-1] */
  for (i = 1; i < M->Ni; i++)
  {
    while (OperatorPsiRC(H, PsiBasis[i - 1], &(PsiBasis[i])))
    {                           /* Allocation of Memory failed */
      printf("Allocation of new memory failed in LanczosTriDiagonalizeKrylovRC for the function OperatorPsi\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
             epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      stop = 1;
      for (j = 0; j < i; j++)
      {
        if (PsiBasis[j].NBlok != 1)
        {
          ComplexWaveFunctionTruncateBasis(&(PsiBasis[j]), epsilon);
          stop = 0;
        }
      }
      if (stop == 1)
        AllocFailed();
    }
    for (j = i; j > 0;)
    {
      j--;
      ComplexWaveFunctionDotProduct(&PsiBasis[j], &PsiBasis[i], &(M->VRe[j][i - 1]), &(M->VIm[j][i - 1]));
      while (ComplexWaveFunctionAddComplexMultipliedUnionBasis
             (&(PsiBasis[i]), PsiBasis[j], -M->VRe[j][i - 1], -M->VIm[j][i - 1]))
      {                         /* Allocation of Memory failed */
        printf
          ("Allocation of new memory failed in LanczosTriDiagonalizeKrylovRC for the function RealWaveFunctionAddMultipliedUnionBasis(&(PsiBasis[i]),PsiBasis[j],-dot)\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
               epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        temp = ComplexWaveFunctionCalculateNorm(PsiBasis[i]);
        stop = 1;
        if (PsiBasis[i].NBlok != 1)
        {
          ComplexWaveFunctionTruncateBasis(&(PsiBasis[i]), epsilon * temp);
          stop = 0;
        }
        for (j = 0; j < i; j++)
        {
          if (PsiBasis[j].NBlok != 1)
          {
            ComplexWaveFunctionTruncateBasis(&(PsiBasis[j]), epsilon);
            stop = 0;
          }
        }
        if (stop == 1)
          AllocFailed();
      }
    }
    M->VIm[i][i - 1] = 0;
    M->VRe[i][i - 1] = ComplexWaveFunctionCalculateNorm(PsiBasis[i]);
    ComplexWaveFunctionMultiply(PsiBasis[i], 1 / M->VRe[i][i - 1]);
    printf("BasisSetSize: %12i ", PsiBasis[i].N);
    ComplexWaveFunctionTruncateBasis(&(PsiBasis[i]), epsilon);
    printf("%12i ", PsiBasis[i].N);
    printf("| M.a[%3i] M.b[%3i]: %12.5E %12.5E | Err: %12.5E\n", i, i, M->VRe[i - 1][i - 1], M->VRe[i][i - 1],
           errorestimate);
  }
  /* Calculate column M->Ni */
  lastcolumn.Ni = M->Ni;
  lastcolumn.Nj = 1;
  lastcolumn.QComplex = 1;
  InitMatrix(&lastcolumn);
  PsiPrimeOperatorPsiMatrix(H, &(PsiBasis[0]), &(PsiBasis[M->Ni - 1]), &lastcolumn);
  for (i = 0; i < M->Ni; i++)
  {
    M->VRe[i][M->Nj - 1] = lastcolumn.VRe[i][0];
    M->VIm[i][M->Nj - 1] = lastcolumn.VIm[i][0];
  }
  FreeMatrix(&lastcolumn);
  fflush(stdout);
  return 0;
}

uint32_t ArnoldiTriDiagonalizeKrylovCC(const OperatorType H, WaveFunctionType psi, double epsilon, MatrixType * const M,
                                       WaveFunctionType * const PsiBasis)
{
  printf("\nStart of ArnoldiTriDiagonalizeKrylovCC\n\n");
  fflush(stdout);
  uint32_t i, j, stop;
  double temp;
  double errorestimate = 1;
  MatrixType lastcolumn;
  /* test input */
  if (M->QComplex == 0)
  {
    printf
      ("ArnoldiTriDiagonalizeKrylovCC entered with real matrix for results/n Although this might be intentional the alghoritms assumes a complex matrix to store the results\n Matrix turned into complex matrix");
    MakeRealMatrixComplex(M);
  }
  /* copy the name of psi to PsiBasis[0] in order to do the loop more cleanely  */
  PsiBasis[0].N = psi.N;
  PsiBasis[0].NBlok = psi.NBlok;
  PsiBasis[0].NFermion = psi.NFermion;
  PsiBasis[0].NBoson = psi.NBoson;
  PsiBasis[0].NWord = psi.NWord;
  PsiBasis[0].NBitsKey = psi.NBitsKey;
  PsiBasis[0].ValueRe = psi.ValueRe;
  PsiBasis[0].ValueIm = psi.ValueIm;
  PsiBasis[0].Determinant = psi.Determinant;
  PsiBasis[0].Next = psi.Next;
  PsiBasis[0].KeyList = psi.KeyList;
  PsiBasis[0].QComplex = psi.QComplex;
  /* calculate a[0] and b[0]  */
  M->a0 = 0;
  M->b0 = ComplexWaveFunctionCalculateNorm(PsiBasis[0]);
  /* Normalize PsiBasis[0] */
  if (M->b0 == 0)
  {
    ComplexWaveFunctionMultiply(PsiBasis[0], (double)0.0);
  }
  else
  {
    ComplexWaveFunctionMultiply(PsiBasis[0], (double)1.0 / M->b0);
  }
  printf("BasisSetSize: %12i ", PsiBasis[0].N);
  ComplexWaveFunctionTruncateBasis(&(PsiBasis[0]), epsilon);
  printf("%12i ", PsiBasis[0].N);
  printf("| M.a[  0] M.b[  0]: %12.5E %12.5E | Err: %12.5E\n", M->a0, M->b0, errorestimate);
  /* Create PsiBasis[1] upto PsiBasis[M-N-1] */
  for (i = 1; i < M->Ni; i++)
  {
    while (OperatorPsiCC(H, PsiBasis[i - 1], &(PsiBasis[i])))
    {                           /* Allocation of Memory failed */
      printf("Allocation of new memory failed in LanczosTriDiagonalizeKrylovCC for the function OperatorPsi\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
             epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      stop = 1;
      for (j = 0; j < i; j++)
      {
        if (PsiBasis[j].NBlok != 1)
        {
          ComplexWaveFunctionTruncateBasis(&(PsiBasis[j]), epsilon);
          stop = 0;
        }
      }
      if (stop == 1)
        AllocFailed();
    }
    for (j = i; j > 0;)
    {
      j--;
      ComplexWaveFunctionDotProduct(&PsiBasis[j], &PsiBasis[i], &(M->VRe[j][i - 1]), &(M->VIm[j][i - 1]));
      while (ComplexWaveFunctionAddComplexMultipliedUnionBasis
             (&(PsiBasis[i]), PsiBasis[j], -M->VRe[j][i - 1], -M->VIm[j][i - 1]))
      {                         /* Allocation of Memory failed */
        printf
          ("Allocation of new memory failed in LanczosTriDiagonalizeKrylovCC for the function RealWaveFunctionAddMultipliedUnionBasis(&(PsiBasis[i]),PsiBasis[j],-dot)\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
               epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        temp = ComplexWaveFunctionCalculateNorm(PsiBasis[i]);
        stop = 1;
        if (PsiBasis[i].NBlok != 1)
        {
          ComplexWaveFunctionTruncateBasis(&(PsiBasis[i]), epsilon * temp);
          stop = 0;
        }
        for (j = 0; j < i; j++)
        {
          if (PsiBasis[j].NBlok != 1)
          {
            ComplexWaveFunctionTruncateBasis(&(PsiBasis[j]), epsilon);
            stop = 0;
          }
        }
        if (stop == 1)
          AllocFailed();
      }
    }
    M->VIm[i][i - 1] = 0;
    M->VRe[i][i - 1] = ComplexWaveFunctionCalculateNorm(PsiBasis[i]);
    ComplexWaveFunctionMultiply(PsiBasis[i], 1 / M->VRe[i][i - 1]);
    printf("BasisSetSize: %12i ", PsiBasis[i].N);
    ComplexWaveFunctionTruncateBasis(&(PsiBasis[i]), epsilon);
    printf("%12i ", PsiBasis[i].N);
    printf("| M.a[%3i] M.b[%3i]: %12.5E %12.5E | Err: %12.5E\n", i, i, M->VRe[i - 1][i - 1], M->VRe[i][i - 1],
           errorestimate);
  }
  /* Calculate column M->Ni */
  lastcolumn.Ni = M->Ni;
  lastcolumn.Nj = 1;
  lastcolumn.QComplex = 1;
  InitMatrix(&lastcolumn);
  PsiPrimeOperatorPsiMatrix(H, &(PsiBasis[0]), &(PsiBasis[M->Ni - 1]), &lastcolumn);
  for (i = 0; i < M->Ni; i++)
  {
    M->VRe[i][M->Nj - 1] = lastcolumn.VRe[i][0];
    M->VIm[i][M->Nj - 1] = lastcolumn.VIm[i][0];
  }
  FreeMatrix(&lastcolumn);
  fflush(stdout);
  return 0;
}

uint32_t ArnoldiTriDiagonalizeKrylovMC(const OperatorType H, WaveFunctionType psi, double epsilon, MatrixType * const M,
                                       WaveFunctionType * const PsiBasis)
{
  printf("\nStart of ArnoldiTriDiagonalizeKrylovMC\n\n");
  fflush(stdout);
  uint32_t i, j, stop;
  double temp;
  double errorestimate = 1;
  MatrixType lastcolumn;
  /* test input */
  if (M->QComplex == 0)
  {
    printf
      ("ArnoldiTriDiagonalizeKrylovMC entered with real matrix for results/n Although this might be intentional the alghoritms assumes a complex matrix to store the results\n Matrix turned into complex matrix");
    MakeRealMatrixComplex(M);
  }
  /* copy the name of psi to PsiBasis[0] in order to do the loop more cleanely  */
  PsiBasis[0].N = psi.N;
  PsiBasis[0].NBlok = psi.NBlok;
  PsiBasis[0].NFermion = psi.NFermion;
  PsiBasis[0].NBoson = psi.NBoson;
  PsiBasis[0].NWord = psi.NWord;
  PsiBasis[0].NBitsKey = psi.NBitsKey;
  PsiBasis[0].ValueRe = psi.ValueRe;
  PsiBasis[0].ValueIm = psi.ValueIm;
  PsiBasis[0].Determinant = psi.Determinant;
  PsiBasis[0].Next = psi.Next;
  PsiBasis[0].KeyList = psi.KeyList;
  PsiBasis[0].QComplex = psi.QComplex;
  /* calculate a[0] and b[0]  */
  M->a0 = 0;
  M->b0 = ComplexWaveFunctionCalculateNorm(PsiBasis[0]);
  /* Normalize PsiBasis[0] */
  if (M->b0 == 0)
  {
    ComplexWaveFunctionMultiply(PsiBasis[0], (double)0.0);
  }
  else
  {
    ComplexWaveFunctionMultiply(PsiBasis[0], (double)1.0 / M->b0);
  }
  printf("BasisSetSize: %12i ", PsiBasis[0].N);
  ComplexWaveFunctionTruncateBasis(&(PsiBasis[0]), epsilon);
  printf("%12i ", PsiBasis[0].N);
  printf("| M.a[  0] M.b[  0]: %12.5E %12.5E | Err: %12.5E\n", M->a0, M->b0, errorestimate);
  /* Create PsiBasis[1] upto PsiBasis[M-N-1] */
  for (i = 1; i < M->Ni; i++)
  {
    while (OperatorPsiMC(H, PsiBasis[i - 1], &(PsiBasis[i])))
    {                           /* Allocation of Memory failed */
      printf("Allocation of new memory failed in LanczosTriDiagonalizeKrylovMC for the function OperatorPsi\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
             epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      stop = 1;
      for (j = 0; j < i; j++)
      {
        if (PsiBasis[j].NBlok != 1)
        {
          ComplexWaveFunctionTruncateBasis(&(PsiBasis[j]), epsilon);
          stop = 0;
        }
      }
      if (stop == 1)
        AllocFailed();
    }
    for (j = i; j > 0;)
    {
      j--;
      ComplexWaveFunctionDotProduct(&PsiBasis[j], &PsiBasis[i], &(M->VRe[j][i - 1]), &(M->VIm[j][i - 1]));
      while (ComplexWaveFunctionAddComplexMultipliedUnionBasis
             (&(PsiBasis[i]), PsiBasis[j], -M->VRe[j][i - 1], -M->VIm[j][i - 1]))
      {                         /* Allocation of Memory failed */
        printf
          ("Allocation of new memory failed in LanczosTriDiagonalizeKrylovMC for the function RealWaveFunctionAddMultipliedUnionBasis(&(PsiBasis[i]),PsiBasis[j],-dot)\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
               epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        temp = ComplexWaveFunctionCalculateNorm(PsiBasis[i]);
        stop = 1;
        if (PsiBasis[i].NBlok != 1)
        {
          ComplexWaveFunctionTruncateBasis(&(PsiBasis[i]), epsilon * temp);
          stop = 0;
        }
        for (j = 0; j < i; j++)
        {
          if (PsiBasis[j].NBlok != 1)
          {
            ComplexWaveFunctionTruncateBasis(&(PsiBasis[j]), epsilon);
            stop = 0;
          }
        }
        if (stop == 1)
          AllocFailed();
      }
    }
    M->VIm[i][i - 1] = 0;
    M->VRe[i][i - 1] = ComplexWaveFunctionCalculateNorm(PsiBasis[i]);
    ComplexWaveFunctionMultiply(PsiBasis[i], 1 / M->VRe[i][i - 1]);
    printf("BasisSetSize: %12i ", PsiBasis[i].N);
    ComplexWaveFunctionTruncateBasis(&(PsiBasis[i]), epsilon);
    printf("%12i ", PsiBasis[i].N);
    printf("| M.a[%3i] M.b[%3i]: %12.5E %12.5E | Err: %12.5E\n", i, i, M->VRe[i - 1][i - 1], M->VRe[i][i - 1],
           errorestimate);
  }
  /* Calculate column M->Ni */
  lastcolumn.Ni = M->Ni;
  lastcolumn.Nj = 1;
  lastcolumn.QComplex = 1;
  InitMatrix(&lastcolumn);
  PsiPrimeOperatorPsiMatrix(H, &(PsiBasis[0]), &(PsiBasis[M->Ni - 1]), &lastcolumn);
  for (i = 0; i < M->Ni; i++)
  {
    M->VRe[i][M->Nj - 1] = lastcolumn.VRe[i][0];
    M->VIm[i][M->Nj - 1] = lastcolumn.VIm[i][0];
  }
  FreeMatrix(&lastcolumn);
  fflush(stdout);
  return 0;
}





uint32_t LanczosGroundState(const OperatorType H, WaveFunctionType * const psi, double *epsilon, const double zero,
                            double shiftin, WaveFunctionType * OtherPsi, uint32_t NOtherPsi)
{
  if (psi->QComplex)
  {
    if (H.QComplex == 0)
    {
      return LanczosGroundStateRC(H, psi, epsilon, zero, shiftin, OtherPsi, NOtherPsi);
    }
    else if (H.QComplex == 1)
    {
      return LanczosGroundStateCC(H, psi, epsilon, zero, shiftin, OtherPsi, NOtherPsi);
    }
    else
    {
      return LanczosGroundStateMC(H, psi, epsilon, zero, shiftin, OtherPsi, NOtherPsi);
    }
  }
  else
  {
    if (H.QComplex == 0)
    {
      return LanczosGroundStateRR(H, psi, epsilon, zero, shiftin, OtherPsi, NOtherPsi);
    }
    else if (H.QComplex == 1)
    {
      if (MakeComplexWaveFunctionFromReal(psi, *psi))
      {
        printf("MakeComplexWaveFunctionFromReal in LanczosGroundState failed\n");
        return 1;
      }
      return LanczosGroundStateCC(H, psi, epsilon, zero, shiftin, OtherPsi, NOtherPsi);
    }
    else
    {
      if (MakeComplexWaveFunctionFromReal(psi, *psi))
      {
        printf("MakeComplexWaveFunctionFromReal in LanczosGroundState failed\n");
        return 1;
      }
      return LanczosGroundStateMC(H, psi, epsilon, zero, shiftin, OtherPsi, NOtherPsi);
    }
  }
}


uint32_t LanczosGroundStateRestricted(const OperatorType H, WaveFunctionType * const psi, DeterminantRestrictionType const Restrictions, double *epsilon, const double zero, double shiftin, WaveFunctionType * OtherPsi, uint32_t NOtherPsi)
{
  if (psi->QComplex)
  {
    if (H.QComplex == 0)
    {
      return LanczosGroundStateRestrictedRC(H, psi, Restrictions, epsilon, zero, shiftin, OtherPsi, NOtherPsi);
    }
    else if (H.QComplex == 1)
    {
      return LanczosGroundStateRestrictedCC(H, psi, Restrictions, epsilon, zero, shiftin, OtherPsi, NOtherPsi);
    }
    else
    {
      return LanczosGroundStateRestrictedMC(H, psi, Restrictions, epsilon, zero, shiftin, OtherPsi, NOtherPsi);
    }
  }
  else
  {
    if (H.QComplex == 0)
    {
      return LanczosGroundStateRestrictedRR(H, psi, Restrictions, epsilon, zero, shiftin, OtherPsi, NOtherPsi);
    }
    else if (H.QComplex == 1)
    {
      if (MakeComplexWaveFunctionFromReal(psi, *psi))
      {
        printf("MakeComplexWaveFunctionFromReal in LanczosGroundState failed\n");
        return 1;
      }
      return LanczosGroundStateRestrictedCC(H, psi, Restrictions, epsilon, zero, shiftin, OtherPsi, NOtherPsi);
    }
    else
    {
      if (MakeComplexWaveFunctionFromReal(psi, *psi))
      {
        printf("MakeComplexWaveFunctionFromReal in LanczosGroundState failed\n");
        return 1;
      }
      return LanczosGroundStateRestrictedMC(H, psi, Restrictions, epsilon, zero, shiftin, OtherPsi, NOtherPsi);
    }
  }
}

uint32_t LanczosEigenSystem(const OperatorType H, WaveFunctionType const psiStart, WaveFunctionType * const psiList, double *epsilon, const double zero, uint32_t Npsi)
{
  uint32_t i;
  for (i = 0; i < Npsi; i++)
  {
    WaveFunctionCopy(psiStart, &(psiList[i]));
    WaveFunctionOrthogonalize(&(psiList[i]), psiList, i);
    WaveFunctionNormalize(&(psiList[i]));
    LanczosGroundState(H, &(psiList[i]), epsilon, zero, 0, psiList, i);
  }
  return 0;
}

uint32_t LanczosEigenSystemRestricted(const OperatorType H, WaveFunctionType const psiStart, WaveFunctionType * const psiList, DeterminantRestrictionType const Restrictions, double *epsilon, const double zero, uint32_t Npsi)
{
  uint32_t i;
  for (i = 0; i < Npsi; i++)
  {
    WaveFunctionCopy(psiStart, &(psiList[i]));
    WaveFunctionOrthogonalize(&(psiList[i]), psiList, i);
    WaveFunctionNormalize(&(psiList[i]));
    LanczosGroundStateRestricted(H, &(psiList[i]), Restrictions, epsilon, zero, 0, psiList, i);
  }
  return 0;
}




uint32_t LanczosTriDiagonalize(const OperatorType H, WaveFunctionType psi, double epsilon, TriDiagonalMatrixType * const M, double Gamma)
{
  WaveFunctionType Cpsi;
  if (psi.QComplex)
  {
    if (H.QComplex == 0)
    {
      return LanczosTriDiagonalizeRC(H, psi, epsilon, M, Gamma);
    }
    else if (H.QComplex == 1)
    {
      return LanczosTriDiagonalizeCC(H, psi, epsilon, M, Gamma);
    }
    else
    {
      return LanczosTriDiagonalizeMC(H, psi, epsilon, M, Gamma);
    }
  }
  else
  {
    if (H.QComplex == 0)
    {
      return LanczosTriDiagonalizeRR(H, psi, epsilon, M, Gamma);
    }
    else if (H.QComplex == 1)
    {
      if (MakeComplexWaveFunctionFromReal(&Cpsi, psi))
      {
        printf("MakeComplexWaveFunctionFromReal failed in LanczosTriDiagonalize\n");
        return 1;
      }
      return LanczosTriDiagonalizeCC(H, Cpsi, epsilon, M, Gamma);
    }
    else
    {
      if (MakeComplexWaveFunctionFromReal(&Cpsi, psi))
      {
        printf("MakeComplexWaveFunctionFromReal failed in LanczosTriDiagonalize\n");
        return 1;
      }
      return LanczosTriDiagonalizeMC(H, Cpsi, epsilon, M, Gamma);
    }
  }
}


uint32_t LanczosTriDiagonalizeRestricted(const OperatorType H, WaveFunctionType psi, DeterminantRestrictionType Restrictions, double epsilon, TriDiagonalMatrixType * const M, double Gamma)
{
  WaveFunctionType Cpsi;
  if (psi.QComplex)
  {
    if (H.QComplex == 0)
    {
      return LanczosTriDiagonalizeRestrictedRC(H, psi, Restrictions, epsilon, M, Gamma);
    }
    else if (H.QComplex == 1)
    {
      return LanczosTriDiagonalizeRestrictedCC(H, psi, Restrictions, epsilon, M, Gamma);
    }
    else
    {
      return LanczosTriDiagonalizeRestrictedMC(H, psi, Restrictions, epsilon, M, Gamma);
    }
  }
  else
  {
    if (H.QComplex == 0)
    {
      return LanczosTriDiagonalizeRestrictedRR(H, psi, Restrictions, epsilon, M, Gamma);
    }
    else if (H.QComplex == 1)
    {
      if (MakeComplexWaveFunctionFromReal(&Cpsi, psi))
      {
        printf("MakeComplexWaveFunctionFromReal failed in LanczosTriDiagonalize\n");
        return 1;
      }
      return LanczosTriDiagonalizeRestrictedCC(H, Cpsi, Restrictions, epsilon, M, Gamma);
    }
    else
    {
      if (MakeComplexWaveFunctionFromReal(&Cpsi, psi))
      {
        printf("MakeComplexWaveFunctionFromReal failed in LanczosTriDiagonalize\n");
        return 1;
      }
      return LanczosTriDiagonalizeRestrictedMC(H, Cpsi, Restrictions, epsilon, M, Gamma);
    }
  }
}

uint32_t LanczosTriDiagonalizePlusDensity(const OperatorType H, WaveFunctionType psi, double epsilon,
                                          TriDiagonalMatrixType * const M, MatrixType * const Rho, double Gamma)
{
  WaveFunctionType Cpsi;
  if (psi.QComplex)
  {
    if (H.QComplex == 0)
    {
      return LanczosTriDiagonalizePlusDensityRC(H, psi, epsilon, M, Rho, Gamma);
    }
    else if (H.QComplex == 1)
    {
      return LanczosTriDiagonalizePlusDensityCC(H, psi, epsilon, M, Rho, Gamma);
    }
    else
    {
      return LanczosTriDiagonalizePlusDensityMC(H, psi, epsilon, M, Rho, Gamma);
    }
  }
  else
  {
    if (H.QComplex == 0)
    {
      return LanczosTriDiagonalizePlusDensityRR(H, psi, epsilon, M, Rho, Gamma);
    }
    else if (H.QComplex == 1)
    {
      if (MakeComplexWaveFunctionFromReal(&Cpsi, psi))
      {
        printf("MakeComplexWaveFunctionFromReal failed in LanczosTriDiagonalize\n");
        return 1;
      }
      return LanczosTriDiagonalizePlusDensityCC(H, Cpsi, epsilon, M, Rho, Gamma);
    }
    else
    {
      if (MakeComplexWaveFunctionFromReal(&Cpsi, psi))
      {
        printf("MakeComplexWaveFunctionFromReal failed in LanczosTriDiagonalize\n");
        return 1;
      }
      return LanczosTriDiagonalizePlusDensityMC(H, Cpsi, epsilon, M, Rho, Gamma);
    }
  }
}


/*   A subroutine to create the krylov space for a given Matrix H and starting function psi. 
 *   Used is a Lanczos routine with reorthogonalization
 *
 *   On input are given:
 *   OperatorType                H : The matrix
 *   WaveFunctionType          psi : The starting vector
 *   double                epsilon : size of smallest element in psi considered different from zero
 *   double                  Gamma : accuracy needed (assumed broadening when looked at the resulting continued fraction - distance from the real axis where G should be accurate
 *   TriDiagonalMatrixType       M : M->N sets the size of the Krylov space.
 *
 *   On output created are:
 *   TriDiagonalMatrixType       M : M->a and M->b contain the tridiagonal matrix
 *   WaveFunctionType    *PsiBasis : pointer to a vector of length M->N wavefunctions representing the Krylov basis
 *                                 : the wavefunctions are not allocated before the calculation starts
 *
 *   The function psi will be changed and become part of *PsiBasis
 *
 *
 *   The resulting Matrix has the form:
 *
 *  The element a[0] stores a possilbe real ofset when used to represent a continued fraction
 *  The element b[0] stores a possile weighting when used to represent a continued fraction
 *  The element a[i] stores the elements M[i-1][i-1] of the full matrix for i>=1
 *  The element b[i] stores the elements M[i][i-1]==M[i-1][i] of the full matrix for i>=1
 *  a is a vector of lenght N+1
 *  b is a vector of length N
 *  The matrix represented has a dimension of N*N
 *
 *     ( A[0][0]=a1 | A[0][1]=b1 | 0          | 0          )
 *   A=( A[1][0]=b1 | A[1][1]=a2 | A[1][2]=b2 | 0          )
 *     ( 0          | A[2][1]=b2 | A[2][2]=a3 | A[2][3]=b3 )
 *     ( 0          | 0          | A[3][2]=b3 | A[3][3]=a4 )
 *
 * 
 */
uint32_t LanczosTriDiagonalizeKrylov(const OperatorType H, WaveFunctionType psi, double epsilon,
                                     TriDiagonalMatrixType * const M, WaveFunctionType * const PsiBasis)
{
  WaveFunctionType Cpsi;
  if (psi.QComplex)
  {
    if (H.QComplex == 0)
    {
      return LanczosTriDiagonalizeKrylovRC(H, psi, epsilon, M, PsiBasis);
    }
    else if (H.QComplex == 1)
    {
      return LanczosTriDiagonalizeKrylovCC(H, psi, epsilon, M, PsiBasis);
    }
    else
    {
      return LanczosTriDiagonalizeKrylovMC(H, psi, epsilon, M, PsiBasis);
    }
  }
  else
  {
    if (H.QComplex == 0)
    {
      return LanczosTriDiagonalizeKrylovRR(H, psi, epsilon, M, PsiBasis);
    }
    else if (H.QComplex == 1)
    {
      if (MakeComplexWaveFunctionFromReal(&Cpsi, psi))
      {
        printf("MakeComplexWaveFunctionFromReal failed in LanczosTriDiagonalizeKrylov\n");
        return 1;
      }
      return LanczosTriDiagonalizeKrylovCC(H, Cpsi, epsilon, M, PsiBasis);
    }
    else
    {
      if (MakeComplexWaveFunctionFromReal(&Cpsi, psi))
      {
        printf("MakeComplexWaveFunctionFromReal failed in LanczosTriDiagonalizeKrylov\n");
        return 1;
      }
      return LanczosTriDiagonalizeKrylovMC(H, Cpsi, epsilon, M, PsiBasis);
    }
  }
}

/*
 *   
 */
uint32_t LanczosTriDiagonalizeKrylovPerturbation(const OperatorPerturbationType H, WaveFunctionPerturbationType psi,
                                                 double epsilon, TriDiagonalMatrixType * const M,
                                                 WaveFunctionPerturbationType * const PsiBasis)
{
  printf("\nStart of LanczosTriDiagonalizeKrylovRR\n\n");
  printf("Still to correct for complex systems\n");
  fflush(stdout);
  uint32_t i, j, k, stop;
  double temp, dotRe, dotIm;
  double errorestimate = 1;
  /* copy the name of psi to PsiBasis[0] in order to do the loop more cleanely  */
  PsiBasis[0].NV = psi.NV;
  PsiBasis[0].order = psi.order;
  PsiBasis[0].psiOrder = psi.psiOrder;
  /* calculate a[0] and b[0]  */
  M->a[0] = 0;
  WaveFunctionPerturbationCalculateNorm(PsiBasis[0], &(M->b[0]));
  /* Normalize PsiBasis[0] */
  if (M->b[0] == 0)
  {
    WaveFunctionPerturbationMultiply(PsiBasis[0], (double)0.0, (double)0.0);
  }
  else
  {
    WaveFunctionPerturbationMultiply(PsiBasis[0], (double)1.0 / M->b[0], (double)0.0);
  }
  printf("BasisSetSize: %12i ", PsiBasis[0].psiOrder[0][0].N);
  WaveFunctionPerturbationTruncateBasis(&(PsiBasis[0]), epsilon / errorestimate);
  printf("%12i ", PsiBasis[0].psiOrder[0][0].N);
  printf("| M.a[  0]         : %12.5E\n", M->a[0]);
  /* Create PsiBasis[1] upto PsiBasis[M-N-1] */
  for (i = 1; i < M->N; i++)
  {
    while (OperatorPsiPerturbation(H, PsiBasis[i - 1], &(PsiBasis[i])))
    {                           /* Allocation of Memory failed */
      printf("Allocation of new memory failed in LanczosTriDiagonalizeKrylovRR for the function OperatorPsi\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
             epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      stop = 1;
      for (j = 0; j < i; j++)
      {
        /* Only truncate if there are parts of psi that have more than one block allocated */
        if (!QWaveFunctionPerturbationMinimumBloksAlocated(PsiBasis[j]))
        {
          WaveFunctionPerturbationTruncateBasis(&(PsiBasis[j]), epsilon / errorestimate);
          stop = 0;
        }
      }
      if (stop == 1)
        AllocFailed();
    }
    for (j = i; j > 0;)
    {
      j--;
      WaveFunctionPerturbationDotProduct(PsiBasis[j], PsiBasis[i], &dotRe, &dotIm);
      if (j == i - 1)
      {
        M->a[i] = dotRe;
      }
      while (WaveFunctionPerturbationAddMultipliedUnionBasis(&(PsiBasis[i]), PsiBasis[j], -dotRe))
      {                         /* Allocation of Memory failed */
        printf
          ("Allocation of new memory failed in LanczosTriDiagonalizeKrylovRR for the function RealWaveFunctionAddMultipliedUnionBasis(&(PsiBasis[i]),PsiBasis[j],-dot)\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon,
               epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        WaveFunctionPerturbationCalculateNorm(PsiBasis[i], &temp);
        stop = 1;
        if (!QWaveFunctionPerturbationMinimumBloksAlocated(PsiBasis[j]))
        {
          WaveFunctionPerturbationTruncateBasis(&(PsiBasis[i]), epsilon * temp / errorestimate);
          stop = 0;
        }
        for (k = 0; k < i; k++)
        {
          if (!QWaveFunctionPerturbationMinimumBloksAlocated(PsiBasis[j]))
          {
            WaveFunctionPerturbationTruncateBasis(&(PsiBasis[k]), epsilon / errorestimate);
            stop = 0;
          }
        }
        if (stop == 1)
          AllocFailed();
      }
    }
    WaveFunctionPerturbationCalculateNorm(PsiBasis[i], &dotRe);
    WaveFunctionPerturbationMultiply(PsiBasis[i], 1 / dotRe, 0);
    M->b[i] = dotRe;
    errorestimate *= Sqr(M->b[i]) / (Sqr(M->b[i]) + DBL_MIN + DBL_EPSILON * Sqr(M->a[i]));
    printf("BasisSetSize: %12i ", PsiBasis[i].psiOrder[0][0].N);
    WaveFunctionPerturbationTruncateBasis(&(PsiBasis[i]), epsilon / errorestimate);
    printf("%12i ", PsiBasis[i].psiOrder[0][0].N);
    printf("| M.a[%3i] M.b[%3i]: %12.5E %12.5E | Err: %12.5E\n", i, i - 1, M->a[i], M->b[i - 1], errorestimate);
  }
  PsiHermitianOperatorPsiPerturbation(H, PsiBasis[M->N - 1], &(M->a[M->N]));
  printf("BasisSetSize: %12i ", PsiBasis[M->N - 1].psiOrder[0][0].N);
  printf("             | M.a[%3i] M.b[%3i]: %12.5E %12.5E | Err: %12.5E\n\n", M->N, M->N - 1, M->a[M->N],
         M->b[M->N - 1], errorestimate);
  fflush(stdout);
  return 0;
}



/*   A subroutine to kreate the krylov space for a given Matrix H and starting function psi.
 *   Used is an Arnoldi routine
 *   For Hermitian matrices H one should use the Lanczos equivalent
 *   Arnoldi is only usefull for non Hermitian matrices H
 *
 *   On input are given:
 *   OperatorType                H : The matrix
 *   WaveFunctionType          psi : The starting vector
 *   double                epsilon : size of smallest element in psi considered different from zero
 *   double                  Gamma : accuracy needed (assumed broadening when looked at the resulting continued fraction - distance from the real axis where G should be accurate
 *   MatrixType                  M : M->Ni=M->Nj sets the size of the Krylov space.
 *
 *   On output created are:
 *   MatrixType                  M : Uper triangle j>=i-1 contains the upper Hessenberg matrix, lower triangle is zero.
 *   WaveFunctionType    *PsiBasis : pointer to a vector of length M->N wavefunctions representing the Krylov basis
 *                                 : the wavefunctions are not allocated before the calculation starts
 *
 *   The function psi will be changed and become part of *PsiBasis
 *
 *
 *   The resulting Matrix has the form:
 *
 *  The element a0 stores a possilbe real ofset when used to represent a continued fraction
 *  The element b0 stores a possile weighting when used to represent a continued fraction
 *  The matrix represented has a dimension of N*N
 *
 *     ( A[0][0]=h11 | A[0][1]=h12 | A[1][2]=h13 | A[2][3]=h14 )
 *   A=( A[1][0]=h21 | A[1][1]=h22 | A[1][2]=h23 | A[2][3]=h24 )
 *     ( 0           | A[2][1]=h23 | A[2][2]=h33 | A[2][3]=h34 )
 *     ( 0           | 0           | A[3][2]=h43 | A[3][3]=h44 )
 *
 *
 */
uint32_t ArnoldiTriDiagonalizeKrylov(const OperatorType H, WaveFunctionType psi, double epsilon, MatrixType * const M,
                                     WaveFunctionType * const PsiBasis)
{
  WaveFunctionType Cpsi;
  if (psi.QComplex)
  {
    if (H.QComplex == 0)
    {
      return ArnoldiTriDiagonalizeKrylovRC(H, psi, epsilon, M, PsiBasis);
    }
    else if (H.QComplex == 1)
    {
      return ArnoldiTriDiagonalizeKrylovCC(H, psi, epsilon, M, PsiBasis);
    }
    else
    {
      return ArnoldiTriDiagonalizeKrylovMC(H, psi, epsilon, M, PsiBasis);
    }
  }
  else
  {
    if (H.QComplex == 0)
    {
      return ArnoldiTriDiagonalizeKrylovRR(H, psi, epsilon, M, PsiBasis);
    }
    else if (H.QComplex == 1)
    {
      if (MakeComplexWaveFunctionFromReal(&Cpsi, psi))
      {
        printf("MakeComplexWaveFunctionFromReal failed in LanczosTriDiagonalizeKrylov\n");
        return 1;
      }
      return ArnoldiTriDiagonalizeKrylovCC(H, Cpsi, epsilon, M, PsiBasis);
    }
    else
    {
      if (MakeComplexWaveFunctionFromReal(&Cpsi, psi))
      {
        printf("MakeComplexWaveFunctionFromReal failed in LanczosTriDiagonalizeKrylov\n");
        return 1;
      }
      return ArnoldiTriDiagonalizeKrylovMC(H, Cpsi, epsilon, M, PsiBasis);
    }
  }
}



/*    Given a basis set Psi_i for i=0 to NPsi in the vector PsiBasis
 *    Given an operator H12=h12 e12 such that h12 acts within PsiBasis and e12 creates or annialates something somerwhereeles and therefore
 *    <Psi_i | H12 | Psi_j >=<Psi_i | h12 | Psi_j > is zero for all i and j
 *    Create an orthonormal basis that spanns the set of vectors h12 Psi_i
 */
uint32_t ExtendKrylovBasis(const OperatorType h12, WaveFunctionType * const PsiBasis,
                           WaveFunctionType * const PsiBasisExtended, uint32_t NPsi, double epsilon)
{
  uint32_t i, j, stop;
  double norm;
  for (i = 0; i < NPsi; i++)
  {
    while (OperatorPsi(h12, PsiBasis[i], &(PsiBasisExtended[i])))
    {
      printf("OperatorPsi failed in ExtendKrylovBasis\n");
      printf("Reduce accuracy from %g to %g", epsilon, epsilon * 2);
      fflush(stdout);
      epsilon *= 2;
      stop = 1;
      for (j = 0; j < NPsi; j++)
      {
        if (PsiBasis[j].NBlok != 1)
        {
          WaveFunctionTruncateBasis(&(PsiBasis[j]), epsilon);
          stop = 0;
        }
      }
      for (j = 0; j < i; j++)
      {
        if (PsiBasisExtended[j].NBlok != 1)
        {
          WaveFunctionTruncateBasis(&(PsiBasisExtended[j]), epsilon);
          stop = 0;
        }
      }
      if (stop)
      {
        printf("OperatorPsi failed in ExtendKrylovBasis\n");
        printf("Reducing accuracy does not help as all functions have NBlok=1\n");
        printf("Stop at i=%i in ExtendKrylovBasis\n", i);
        AllocFailed();
      }
    }


    WaveFunctionCalculateNorm(PsiBasisExtended[i], &norm);
    printf("In ExtendKrylov Norm%i=%lf\n", i, norm);

    while (WaveFunctionOrthogonalize(&(PsiBasisExtended[i]), PsiBasisExtended, i))
    {
      printf("WaveFunctionOrthogonalize failed in ExtendKrylovBasis\n");
      printf("Reduce accuracy from %g to %g", epsilon, epsilon * 2);
      fflush(stdout);
      epsilon *= 2;
      stop = 1;
      for (j = 0; j < NPsi; j++)
      {
        if (PsiBasis[j].NBlok != 1)
        {
          WaveFunctionTruncateBasis(&(PsiBasis[j]), epsilon);
          stop = 0;
        }
      }
      for (j = 0; j < i; j++)
      {
        if (PsiBasisExtended[j].NBlok != 1)
        {
          WaveFunctionTruncateBasis(&(PsiBasisExtended[j]), epsilon);
          stop = 0;
        }
      }
      if (stop)
      {
        printf("OperatorPsi failed in ExtendKrylovBasis\n");
        printf("Reducing accuracy does not help as all functions have NBlok=1\n");
        printf("Stop at i=%i in ExtendKrylovBasis\n", i);
        AllocFailed();
      }
    }

    WaveFunctionNormalize(&(PsiBasisExtended[i]));
    while (WaveFunctionOrthogonalize(&(PsiBasisExtended[i]), PsiBasisExtended, i))
    {
      printf("WaveFunctionOrthogonalize failed in ExtendKrylovBasis\n");
      printf("Reduce accuracy from %g to %g", epsilon, epsilon * 2);
      fflush(stdout);
      epsilon *= 2;
      stop = 1;
      for (j = 0; j < NPsi; j++)
      {
        if (PsiBasis[j].NBlok != 1)
        {
          WaveFunctionTruncateBasis(&(PsiBasis[j]), epsilon);
          stop = 0;
        }
      }
      for (j = 0; j < i; j++)
      {
        if (PsiBasisExtended[j].NBlok != 1)
        {
          WaveFunctionTruncateBasis(&(PsiBasisExtended[j]), epsilon);
          stop = 0;
        }
      }
      if (stop)
      {
        printf("OperatorPsi failed in ExtendKrylovBasis\n");
        printf("Reducing accuracy does not help as all functions have NBlok=1\n");
        printf("Stop at i=%i in ExtendKrylovBasis\n", i);
        AllocFailed();
      }
    }
    WaveFunctionCalculateNorm(PsiBasisExtended[i], &norm);
    WaveFunctionNormalize(&(PsiBasisExtended[i]));
    if (norm < 0.99)
    {
      while (WaveFunctionOrthogonalize(&(PsiBasisExtended[i]), PsiBasisExtended, i))
      {
        printf("WaveFunctionOrthogonalize failed in ExtendKrylovBasis\n");
        printf("Reduce accuracy from %g to %g", epsilon, epsilon * 2);
        fflush(stdout);
        epsilon *= 2;
        stop = 1;
        for (j = 0; j < NPsi; j++)
        {
          if (PsiBasis[j].NBlok != 1)
          {
            WaveFunctionTruncateBasis(&(PsiBasis[j]), epsilon);
            stop = 0;
          }
        }
        for (j = 0; j < i; j++)
        {
          if (PsiBasisExtended[j].NBlok != 1)
          {
            WaveFunctionTruncateBasis(&(PsiBasisExtended[j]), epsilon);
            stop = 0;
          }
        }
        if (stop)
        {
          printf("OperatorPsi failed in ExtendKrylovBasis\n");
          printf("Reducing accuracy does not help as all functions have NBlok=1\n");
          printf("Stop at i=%i in ExtendKrylovBasis\n", i);
          AllocFailed();
        }
      }
      WaveFunctionCalculateNorm(PsiBasisExtended[i], &norm);
      if (norm < 0.99)
      {
        WaveFunctionMultiply(&PsiBasisExtended[i], 0, 0);
      }
      else
      {
        WaveFunctionMultiply(&PsiBasisExtended[i], 1.0 / norm, 0);
      }
    }
    WaveFunctionTruncateBasis(&(PsiBasisExtended[i]), epsilon);
  }
  return 0;
}


/* Given wavefunction psi and operator O a tridiagonal matrix and its Krylov basis is created of length NKrylov
 * For each operator V[i], with i from 0 to NV-1 the following steps are taken
 *   The extended basis is defined as O^n V[i] PsiKrilov[k], but will never be calculated. (for n element 0 to nmax)
 *   Stemp contains the overlap of the basis functions, representing the non-orthogonality of the basis
 *   Htemp contains the Hamiltonian on this non-orthogonal basis
 *   Calculate the matrix Stemp with elements stemp[j,k] = <PsiKrylov[j] | (O^n V[i])^dag   (O^n' V[i]) PsiKrilov[k]>
 *   Calculate the matrix Htemp with elements htemp[j,k] = <PsiKrylov[j] | (O^n V[i])^dag O (O^n' V[i]) PsiKrilov[k]>
 *   Diagonalize Htemp with respect to the overlap matrix Stemp and store eigenvalues in Hval[i], Eigenvectors in HtempEigenVec
 *   Create Vmat[i] as Vmat[i] = Stemp . ConjuagetTranspose[HtempEigenVec]
 *
 *  Once Hval and VMat are know one can call TridiagonalizeImpurityCoupledToBathDyson
 *
 *  Gout->N contains the length of the newly created tri-diagonal matrix and has to be set on input
 *  Gout will be initiallized within the routine
 *
 *  In order to calculate Stemp and Htemp one first calculates the matrices VOV(m) = V^dag O^m V with m element 0 to 2 nmax+1
 */
uint32_t LanczosTriDiagonalizeDysonFullLowMem(OperatorType const O, uint32_t nmax, WaveFunctionType psi, double epsilon,
                                              uint32_t const NKrylov, OperatorType const *V,
                                              TriDiagonalMatrixType const *g, uint32_t const NV,
                                              TriDiagonalMatrixType * Gout)
{
  uint32_t i, j, n, ni, nj;
  TriDiagonalMatrixType OKrylov;
  MatrixType Stemp, Htemp, HtempEigenFun, StempRoot, VOV[2 * nmax + 2];
  MatrixType Hval, *VMat;
  WaveFunctionType *PsiKrylov;
  OperatorType Otemp, Otempold, VOVtemp;
  printf("\nStart of LanczosTriDiagonalizeDysonLowMem\n\n");
  /* Initiallize matrices */
  Htemp.Ni = (2 * nmax + 1) * NKrylov;
  Htemp.Nj = (2 * nmax + 1) * NKrylov;
  Htemp.QComplex = (O.QComplex || psi.QComplex);
  for (i = 0; i < NV; i++)
  {
    Htemp.QComplex = (Htemp.QComplex || V[i].QComplex);
  }
  if (InitMatrix(&Htemp))
  {
    printf("InitMatrix failed in LanczosTriDiagonalizeDyson\n");
    fflush(stdout);
    return 1;
  };
  Stemp.Ni = (2 * nmax + 1) * NKrylov;
  Stemp.Nj = (2 * nmax + 1) * NKrylov;
  Stemp.QComplex = Htemp.QComplex;
  if (InitMatrix(&Stemp))
  {
    printf("InitMatrix failed in LanczosTriDiagonalizeDyson\n");
    fflush(stdout);
    return 1;
  };
  for (n = 0; n < 2 * nmax + 2; n++)
  {
    VOV[n].Ni = NKrylov;
    VOV[n].Nj = NKrylov;
    VOV[n].QComplex = (O.QComplex || psi.QComplex);
    for (i = 0; i < NV; i++)
    {
      VOV[n].QComplex = (VOV[n].QComplex || V[i].QComplex);
    }
    if (InitMatrix(&VOV[n]))
    {
      printf("InitMatrix failed in LanczosTriDiagonalizeDyson\n");
      fflush(stdout);
      return 1;
    };
  }
  Hval.Ni = NV;
  Hval.Nj = (2 * nmax + 1) * NKrylov;
  Hval.QComplex = 0;
  if (InitMatrix(&Hval))
  {
    printf("InitMatrix failed in LanczosTriDiagonalizeDyson\n");
    fflush(stdout);
    return 1;
  };
  if ((VMat = calloc(NV, sizeof(MatrixType))) == NULL)
  {
    printf("calloc failed in LanczosTriDiagonalizeDyson\n");
    fflush(stdout);
    return 1;
  };
  /* Create Krylov basis of O */
  OKrylov.N = NKrylov;
  InitTriDiagonalMatrix(&OKrylov);
  if ((PsiKrylov = calloc(NKrylov, sizeof(WaveFunctionType))) == NULL)
  {
    printf("calloc failed in LanczosTriDiagonalizeDyson\n");
    fflush(stdout);
    return 1;
  };
  LanczosTriDiagonalizeKrylov(O, psi, epsilon, &OKrylov, PsiKrylov);
  /* loop over all V[i] */
  for (i = 0; i < NV; i++)
  {
    printf("Create perturbation %i from %i\n", i + 1, NV);
    /* create matrices VOV[n] */
    /* initiallize VOVtemp for V^{\dag} V */
    InitOperatorForMult(&VOVtemp, 2, &V[i], &V[i]);
    /* create VOV[0] */
    OperatorMultiply2NormalOrder(&VOVtemp, V[i], V[i], 'D', 'N', &QDetNotOrderedFcFaBcBa);
    PsiOperatorPsiMatrix(VOVtemp, PsiKrylov, &VOV[0]);
    FreeOperator(&VOVtemp);
    for (n = 1; n <= 2 * nmax + 1; n++)
    {
      if (n == 1)
      {
        OperatorCopy(&Otemp, O);
      }
      else
      {
        OperatorCopy(&Otempold, Otemp);
        FreeOperator(&Otemp);
        InitOperatorForMult(&Otemp, 2, &Otempold, &O);
        OperatorMultiplyNormalOrder(&Otemp, O, Otempold, &QDetNotOrderedFcFaBcBa);
        FreeOperator(&Otempold);
      }
      InitOperatorForMult(&VOVtemp, 3, &V[i], &Otemp, &V[i]);
      OperatorMultiply3NormalOrder(&VOVtemp, V[i], Otemp, V[i], 'D', 'N', 'N', &QDetNotOrderedFcFaBcBa);
      PsiOperatorPsiMatrix(VOVtemp, PsiKrylov, &VOV[i]);
      FreeOperator(&VOVtemp);
    }
    FreeOperator(&Otemp);
    /* create Htemp and Stemp from VOVtemp */
    for (ni = 0; ni < 2 * nmax + 1; ni++)
    {
      MatrixAddShifted(&Htemp, VOV[2 * ni + 1], 1, 0, ni * NKrylov, ni * NKrylov);
      MatrixAddShifted(&Stemp, VOV[2 * ni], 1, 0, ni * NKrylov, ni * NKrylov);
      for (nj = ni + 1; nj < 2 * nmax + 1; nj++)
      {
        MatrixAddShifted(&Htemp, VOV[ni + nj + 1], 1, 0, ni * NKrylov, nj * NKrylov);
        MatrixAddShifted(&Stemp, VOV[ni + nj], 1, 0, ni * NKrylov, nj * NKrylov);
        MatrixAddShiftedConjugateTranspose(&Htemp, VOV[ni + nj + 1], 1, 0, nj * NKrylov, ni * NKrylov);
        MatrixAddShiftedConjugateTranspose(&Stemp, VOV[ni + nj], 1, 0, nj * NKrylov, ni * NKrylov);
      }
    }
    /* diagonalize Htemp with repsect to Stemp */
    GeneralizedEigensystemSingular(Htemp, Stemp, &HtempEigenFun, &StempRoot, Hval.VRe[i]);
    MatrixMultiply(&(VMat[i]), HtempEigenFun, Stemp);
    /* Free HtempEigenFun */
    FreeMatrix(&HtempEigenFun);
  }
  /* Free PsiKrylov */
  for (j = 0; j < NKrylov; j++)
  {
    WaveFunctionFree(&(PsiKrylov[j]));
  }
  free(PsiKrylov);
  FreeMatrix(&Htemp);
  FreeMatrix(&Stemp);
  /* create Gout */
  if (TridiagonalizeImpurityCoupledToBathDyson(OKrylov, VMat, Hval, g, NV, Gout))
  {
    printf("TridiagonalizeImpurityCoupledToBathDyson failed in LanczosTriDiagonalizeDyson\n");
    FreeTriDiagonalMatrix(&OKrylov);
    FreeMatrix(&Hval);
    for (i = 0; i < NV; i++)
    {
      FreeMatrix(&(VMat[i]));
    }
    free(VMat);
    return 1;
  }
  /* Free things allocated */
  FreeTriDiagonalMatrix(&OKrylov);
  FreeMatrix(&Hval);
  for (i = 0; i < NV; i++)
  {
    FreeMatrix(&(VMat[i]));
  }
  free(VMat);
  return 0;
}





/* Given wavefunction psi and operator O a tridiagonal matrix and its Krylov basis is created of length NKrylov
 * For each operator V[i], with i from 0 to NV-1 the following steps are taken
 *   The extended basis is defined as V[i] PsiKrilov[k], but will never be calculated.
 *   Stemp contains the overlap of the basis functions, representing the non-orthogonality of the basis
 *   Htemp contains the Hamiltonian on this non-orthogonal basis
 *   Calculate the matrix Stemp with elements stemp[j,k] = <PsiKrylov[j] | V[i]^dag   V[i] PsiKrilov[k]>
 *   Calculate the matrix Htemp with elements htemp[j,k] = <PsiKrylov[j] | V[i]^dag O V[i] PsiKrilov[k]>
 *   Diagonalize Htemp with respect to the overlap matrix Stemp and store eigenvalues in Hval[i], Eigenvectors in HtempEigenVec
 *   Create Vmat[i] as Vmat[i] = Stemp . ConjuagetTranspose[HtempEigenVec]
 *
 *  Once Hval and VMat are know one can call TridiagonalizeImpurityCoupledToBathDyson
 *
 *  Gout->N contains the length of the newly created tri-diagonal matrix and has to be set on input
 *  Gout will be initiallized within the routine
 */
uint32_t LanczosTriDiagonalizeDysonLowMem(OperatorType const O, WaveFunctionType psi, double epsilon,
                                          uint32_t const NKrylov, OperatorType const *V, TriDiagonalMatrixType const *g,
                                          uint32_t const NV, TriDiagonalMatrixType * Gout)
{
  uint32_t i, j;
  TriDiagonalMatrixType OKrylov;
  MatrixType Stemp, Htemp, HtempEigenFun, StempRoot;
  MatrixType Hval, *VMat;
  WaveFunctionType *PsiKrylov;
  OperatorType Otemp;
  TimeType Time;
  InitializeTiming(&Time);
  BeginTiming("Total", &Time);
  printf("\nStart of LanczosTriDiagonalizeDysonLowMem\n\n");
  /* Initiallize matrices */
  Htemp.Ni = NKrylov;
  Htemp.Nj = NKrylov;
  Htemp.QComplex = (O.QComplex || psi.QComplex);
  for (i = 0; i < NV; i++)
  {
    Htemp.QComplex = (Htemp.QComplex || V[i].QComplex);
  }
  if (InitMatrix(&Htemp))
  {
    printf("InitMatrix failed in LanczosTriDiagonalizeDyson\n");
    fflush(stdout);
    return 1;
  };
  Stemp.Ni = NKrylov;
  Stemp.Nj = NKrylov;
  Stemp.QComplex = Htemp.QComplex;
  if (InitMatrix(&Stemp))
  {
    printf("InitMatrix failed in LanczosTriDiagonalizeDyson\n");
    fflush(stdout);
    return 1;
  };
  Hval.Ni = NV;
  Hval.Nj = NKrylov;
  Hval.QComplex = 0;
  if (InitMatrix(&Hval))
  {
    printf("InitMatrix failed in LanczosTriDiagonalizeDyson\n");
    fflush(stdout);
    return 1;
  };
  if ((VMat = calloc(NV, sizeof(MatrixType))) == NULL)
  {
    printf("calloc failed in LanczosTriDiagonalizeDyson\n");
    fflush(stdout);
    return 1;
  };
  /* Create Krylov basis of O */
  OKrylov.N = NKrylov;
  InitTriDiagonalMatrix(&OKrylov);
  if ((PsiKrylov = calloc(NKrylov, sizeof(WaveFunctionType))) == NULL)
  {
    printf("calloc failed in LanczosTriDiagonalizeDyson\n");
    fflush(stdout);
    return 1;
  };
  BeginTiming("Krylov", &Time);
  LanczosTriDiagonalizeKrylov(O, psi, epsilon, &OKrylov, PsiKrylov);
  EndTiming("Krylov", &Time);
  /* loop over all V[i] */
  BeginTiming("Creation of matrix elements", &Time);
  for (i = 0; i < NV; i++)
  {
    printf("Create perturbation %i from %i\n", i + 1, NV);
    /* create extended non orthonormal matrices on original basis */
    /* initiallize Otemp for V^{\dag} V */
    InitOperatorForMult(&Otemp, 2, &V[i], &V[i]);
    /* create Stemp */
    OperatorMultiply2NormalOrder(&Otemp, V[i], V[i], 'D', 'N', &QDetNotOrderedFcFaBcBa);
    PsiOperatorPsiMatrix(Otemp, PsiKrylov, &Stemp);
    FreeOperator(&Otemp);
    /* initiallize Otemp for V^{\dag} O V */
    InitOperatorForMult(&Otemp, 3, &V[i], &O, &V[i]);
    /* create Htemp */
    OperatorMultiply3NormalOrder(&Otemp, V[i], O, V[i], 'D', 'N', 'N', &QDetNotOrderedFcFaBcBa);
    PsiOperatorPsiMatrix(Otemp, PsiKrylov, &Htemp);
    FreeOperator(&Otemp);
    /* diagonalize Htemp with repsect to Stemp */
    GeneralizedEigensystemSingular(Htemp, Stemp, &HtempEigenFun, &StempRoot, Hval.VRe[i]);
    MatrixMultiply(&(VMat[i]), HtempEigenFun, Stemp);
    /* Free HtempEigenFun */
    FreeMatrix(&HtempEigenFun);
  }
  /* Free PsiKrylov */
  for (j = 0; j < NKrylov; j++)
  {
    WaveFunctionFree(&(PsiKrylov[j]));
  }
  free(PsiKrylov);
  FreeMatrix(&Htemp);
  FreeMatrix(&Stemp);
  EndTiming("Creation of matrix elements", &Time);
  /* create Gout */
  BeginTiming("Final Tridiagonalize", &Time);
  if (TridiagonalizeImpurityCoupledToBathDyson(OKrylov, VMat, Hval, g, NV, Gout))
  {
    printf("TridiagonalizeImpurityCoupledToBathDyson failed in LanczosTriDiagonalizeDyson\n");
    FreeTriDiagonalMatrix(&OKrylov);
    FreeMatrix(&Hval);
    for (i = 0; i < NV; i++)
    {
      FreeMatrix(&(VMat[i]));
    }
    free(VMat);
    return 1;
  }
  EndTiming("Final Tridiagonalize", &Time);
  /* Free things allocated */
  FreeTriDiagonalMatrix(&OKrylov);
  FreeMatrix(&Hval);
  for (i = 0; i < NV; i++)
  {
    FreeMatrix(&(VMat[i]));
  }
  free(VMat);
  EndTiming("Total", &Time);
  PrintTiming(Time);
  FreeTiming(&Time);
  return 0;
}


/* Given wavefunction psi and operator O a tridiagonal matrix and its Krylov basis is created of length NKrylov
 * For each operator V[i], with i from 0 to NV-1 the following steps are taken
 *   The first order extended basis is defined as V[i] PsiKrilov[k], but will never be calculated.
 *   Stemp contains the overlap of the basis functions, representing the non-orthogonality of the basis
 *   Htemp contains the Hamiltonian on this non-orthogonal basis
 *   Calculate the matrix Stemp with elements stemp[j,k] = < PsiKrylov[j] | V[i]^dag   V[i] | PsiKrilov[k] >
 *   Calculate the matrix Htemp with elements htemp[j,k] = < PsiKrylov[j] | V[i]^dag O V[i] | PsiKrilov[k] >
 *   Diagonalize Htemp with respect to the overlap matrix Stemp and store eigenvalues in Hval[i], Eigenvectors in Tmat[i]
 *   Create Vmat[i] as Vmat[i] = Stemp . ConjuagetTranspose[Tmat[i]]
 *
 * For each operator V[ii] with ii<i the following steps are taken
 *   The second order extended basis is defined as V[ii] V[i] PsiKrilov[k], but will never be calculated.
 *   Stemp contains the overlap of the basis functions, representing the non-orthogonality of the basis
 *   Htemp contains the Hamiltonian on this non-orthogonal basis
 *   Calculate the matrix Stemp with elements stemp[j,k] = < PsiKrylov[j] | V[i]^dag V[ii]^dag V[ii] V[i] | PsiKrilov[k] >
 *   Calculate the matrix Htemp with elements htemp[j,k] = < PsiKrylov[j] | V[i]^dag V[ii]^dag O V[ii] V[i] | PsiKrilov[k] >
 *   Diagonalize Htemp with respect to the overlap matrix Stemp and store eigenvalues in Hval[i], Eigenvectors in HtempEigenVec
 *   Create Vmat[i] as Vmat[i] = Stemp . ConjuagetTranspose[HtempEigenVec]
 *
 *  Once Hval and VMat are know one can call TridiagonalizeImpurityCoupledToBathDysonSecondOrder
 *
 *  Gout->N contains the length of the newly created tri-diagonal matrix and has to be set on input
 *  Gout will be initiallized within the routine
 */
uint32_t LanczosTriDiagonalizeDysonSecondOrderALowMem(OperatorType const O, WaveFunctionType psi, double epsilon,
                                                      uint32_t const NKrylov, OperatorType const *V,
                                                      TriDiagonalMatrixType const *g, uint32_t const NV,
                                                      TriDiagonalMatrixType * Gout)
{
  uint32_t i, ii, j;
  TriDiagonalMatrixType OKrylov;
  MatrixType Stemp, Htemp, StempRoot, HtempEigenFun, VMat2Temp;
  MatrixType Hval, *VMat, *TMat, Hval2, *VMat2;
  WaveFunctionType *PsiKrylov;
  OperatorType Otemp;
  printf("\nStart of LanczosTriDiagonalizeDysonSecondOrderALowMem\n\n");
  /* Initiallize matrices */
  Htemp.Ni = NKrylov;
  Htemp.Nj = NKrylov;
  Htemp.QComplex = (O.QComplex || psi.QComplex);
  for (i = 0; i < NV; i++)
  {
    Htemp.QComplex = (Htemp.QComplex || V[i].QComplex);
  }
  if (InitMatrix(&Htemp))
  {
    printf("InitMatrix failed in LanczosTriDiagonalizeDyson\n");
    fflush(stdout);
    return 1;
  };
  Stemp.Ni = NKrylov;
  Stemp.Nj = NKrylov;
  Stemp.QComplex = Htemp.QComplex;
  if (InitMatrix(&Stemp))
  {
    printf("InitMatrix failed in LanczosTriDiagonalizeDyson\n");
    fflush(stdout);
    return 1;
  };
  Hval.Ni = NV;
  Hval.Nj = NKrylov;
  Hval.QComplex = 0;
  if (InitMatrix(&Hval))
  {
    printf("InitMatrix failed in LanczosTriDiagonalizeDyson\n");
    fflush(stdout);
    return 1;
  };
  if ((VMat = calloc(NV, sizeof(MatrixType))) == NULL)
  {
    printf("calloc failed in LanczosTriDiagonalizeDyson\n");
    fflush(stdout);
    return 1;
  };
  if ((TMat = calloc(NV, sizeof(MatrixType))) == NULL)
  {
    printf("calloc failed in LanczosTriDiagonalizeDyson\n");
    fflush(stdout);
    return 1;
  };
  Hval2.Ni = (NV * NV - NV) / 2;
  Hval2.Nj = NKrylov;
  Hval2.QComplex = 0;
  if (InitMatrix(&Hval2))
  {
    printf("InitMatrix failed in LanczosTriDiagonalizeDyson\n");
    fflush(stdout);
    return 1;
  };
  if ((VMat2 = calloc((NV * NV - NV), sizeof(MatrixType))) == NULL)
  {
    printf("calloc failed in LanczosTriDiagonalizeDyson\n");
    fflush(stdout);
    return 1;
  };
  /* Create Krylov basis of O */
  OKrylov.N = NKrylov;
  InitTriDiagonalMatrix(&OKrylov);
  if ((PsiKrylov = calloc(NKrylov, sizeof(WaveFunctionType))) == NULL)
  {
    printf("calloc failed in LanczosTriDiagonalizeDyson\n");
    fflush(stdout);
    return 1;
  };
  LanczosTriDiagonalizeKrylov(O, psi, epsilon, &OKrylov, PsiKrylov);
  /* loop over all V[i] */
  for (i = 0; i < NV; i++)
  {
    printf("Create perturbation %i from %i\n", i + 1, NV);
    /* create extended non orthonormal matrices on original basis */
    /* initiallize Otemp for V^{\dag} V */
    InitOperatorForMult(&Otemp, 2, &V[i], &V[i]);
    /* create Stemp */
    OperatorMultiply2NormalOrder(&Otemp, V[i], V[i], 'D', 'N', &QDetNotOrderedFcFaBcBa);
    PsiOperatorPsiMatrix(Otemp, PsiKrylov, &Stemp);
    FreeOperator(&Otemp);
    /* initiallize Otemp for V^{\dag} O V */
    InitOperatorForMult(&Otemp, 3, &V[i], &O, &V[i]);
    /* create Htemp */
    OperatorMultiply3NormalOrder(&Otemp, V[i], O, V[i], 'D', 'N', 'N', &QDetNotOrderedFcFaBcBa);
    PsiOperatorPsiMatrix(Otemp, PsiKrylov, &Htemp);
    FreeOperator(&Otemp);
    GeneralizedEigensystemSingular(Htemp, Stemp, &TMat[i], &StempRoot, Hval.VRe[i]);
    /* create V on orthonormal eigenbasis of Htemp */
    MatrixMultiply(&(VMat[i]), TMat[i], Stemp);
    MatrixTranspose(&(TMat[i]));
  }
  /* loop over all V[ii,i] */
  j = 0;                        /* j counts from 0 to (NV*NV-NV)/2 */
  for (i = 0; i < NV; i++)
  {
    for (ii = 0; ii < i; ii++)
    {
      printf("Create 2nd order perturbation %i %i from %i\n", i + 1, ii + 1, (NV * NV - NV) / 2);
      /* create extended non orthonormal matrices on original basis */
      /* initiallize Otemp for V^{\dag} V^{\dag} V V */
      InitOperatorForMult(&Otemp, 4, &V[i], &V[ii], &V[ii], &V[i]);
      OperatorMultiply4NormalOrder(&Otemp, V[i], V[ii], V[ii], V[i], 'D', 'D', 'N', 'N', &QDetNotOrderedFcFaBcBa);
      PsiOperatorPsiMatrix(Otemp, PsiKrylov, &Stemp);
      FreeOperator(&Otemp);
      /* initiallize Otemp for V^{\dag} V^{\dag} O V V */
      InitOperatorForMult(&Otemp, 5, &V[i], &V[ii], &O, &V[ii], &V[i]);
      /* create Htemp */
      OperatorMultiply5NormalOrder(&Otemp, V[i], V[ii], O, V[ii], V[i], 'D', 'D', 'N', 'N', 'N', &QDetNotOrderedFcFaBcBa);
      PsiOperatorPsiMatrix(Otemp, PsiKrylov, &Htemp);
      FreeOperator(&Otemp);
      /* diagonalize Htemp with repsect to Stemp */
      GeneralizedEigensystemSingular(Htemp, Stemp, &HtempEigenFun, &StempRoot, Hval2.VRe[j]);
      /* create V on orthonormal eigenbasis of Htemp */
      MatrixMultiply(&VMat2Temp, HtempEigenFun, Stemp);
      MatrixMultiply(&(VMat2[2 * j]), VMat2Temp, TMat[i]);
      MatrixMultiply(&(VMat2[2 * j + 1]), VMat2Temp, TMat[ii]);
      /* Free HtempEigenFun */
      FreeMatrix(&HtempEigenFun);
      FreeMatrix(&VMat2Temp);
      j++;                      /* j counts from 0 to (NV*NV-NV)/2 */
    }
  }
  /* Free PsiKrylov */
  for (j = 0; j < NKrylov; j++)
  {
    WaveFunctionFree(&(PsiKrylov[j]));
  }
  free(PsiKrylov);
  FreeMatrix(&Htemp);
  FreeMatrix(&Stemp);
  /* create Gout */
  if (TridiagonalizeImpurityCoupledToBathDysonSecondOrderA(OKrylov, VMat, VMat2, Hval, Hval2, g, NV, Gout))
  {
    printf("TridiagonalizeImpurityCoupledToBathDyson failed in LanczosTriDiagonalizeDyson\n");
    FreeTriDiagonalMatrix(&OKrylov);
    FreeMatrix(&Hval);
    FreeMatrix(&Hval2);
    for (i = 0; i < NV; i++)
    {
      FreeMatrix(&(VMat[i]));
      FreeMatrix(&(TMat[i]));
    }
    for (i = 0; i < (NV * NV - NV); i++)
    {
      FreeMatrix(&(VMat2[i]));
    }
    free(VMat);
    free(VMat2);
    free(TMat);
    return 1;
  }
  /* Free things allocated */
  FreeTriDiagonalMatrix(&OKrylov);
  FreeMatrix(&Hval);
  FreeMatrix(&Hval2);
  for (i = 0; i < NV; i++)
  {
    FreeMatrix(&(VMat[i]));
    FreeMatrix(&(TMat[i]));
  }
  for (i = 0; i < (NV * NV - NV); i++)
  {
    FreeMatrix(&(VMat2[i]));
  }
  free(VMat);
  free(VMat2);
  free(TMat);
  return 0;
}






/*  Given wavefunction psi and operator O a tridiagonal matrix and its Krylov basis is created of length NKrylov
 *  For each operator V[i], with i from 0 to NV-1 the following steps are taken
 *    Create an extended Krylov basis as V[i] PsiKrylov[j]
 *    On this basis calculate the operator O store in matrix Htemp
 *    Between the original and the extended Krylov basis calculate the operator V store in matrix Vtemp
 *    Diagonalize Htemp store eigenvalues in Hval[i]
 *    Transform Vtemp to eigenbasis of Htemp and store new matrix in VMat[i]
 *    Remove Vtemp and Htemp from memory
 *
 *  Once Hval and VMat are known one can call TridiagonalizeImpurityCoupledToBathDyson
 *
 *  Gout->N contains the length of the newly created tri-diagonal matrix and has to be set on input
 *  Gout will be initiallized within the routine
 */
uint32_t LanczosTriDiagonalizeDyson(OperatorType const O, WaveFunctionType psi, double epsilon, uint32_t const NKrylov,
                                    OperatorType const *V, TriDiagonalMatrixType const *g, uint32_t const NV,
                                    TriDiagonalMatrixType * Gout)
{
  uint32_t i, j;
  TriDiagonalMatrixType OKrylov;
  WaveFunctionType *PsiKrylov, *VPsiKrylov;
  MatrixType Htemp, Vtemp, Hval, HtempEigenFun, *VMat;
  /* Create Krylov basis of O */
  OKrylov.N = NKrylov;
  InitTriDiagonalMatrix(&OKrylov);
  if ((PsiKrylov = calloc(NKrylov, sizeof(WaveFunctionType))) == NULL)
  {
    printf("calloc failed in LanczosTriDiagonalizeDyson\n");
    fflush(stdout);
    return 1;
  };
  LanczosTriDiagonalizeKrylov(O, psi, epsilon, &OKrylov, PsiKrylov);
  /* Create extended matrices initiallize */
  if ((VPsiKrylov = calloc(NKrylov, sizeof(WaveFunctionType))) == NULL)
  {
    printf("calloc failed in LanczosTriDiagonalizeDyson\n");
    fflush(stdout);
    return 1;
  };
  Htemp.Ni = NKrylov;
  Htemp.Nj = NKrylov;
  Htemp.QComplex = (O.QComplex || psi.QComplex);
  for (i = 0; i < NV; i++)
  {
    Htemp.QComplex = (Htemp.QComplex || V[i].QComplex);
  }
  if (InitMatrix(&Htemp))
  {
    printf("InitMatrix failed in LanczosTriDiagonalizeDyson\n");
    fflush(stdout);
    return 1;
  };
  Vtemp.Ni = NKrylov;
  Vtemp.Nj = NKrylov;
  Vtemp.QComplex = Htemp.QComplex;
  if (InitMatrix(&Vtemp))
  {
    printf("InitMatrix failed in LanczosTriDiagonalizeDyson\n");
    fflush(stdout);
    return 1;
  };
  Hval.Ni = NV;
  Hval.Nj = NKrylov;
  Hval.QComplex = 0;
  if (InitMatrix(&Hval))
  {
    printf("InitMatrix failed in LanczosTriDiagonalizeDyson\n");
    fflush(stdout);
    return 1;
  };
  if ((VMat = calloc(NV, sizeof(MatrixType))) == NULL)
  {
    printf("calloc failed in LanczosTriDiagonalizeDyson\n");
    fflush(stdout);
    return 1;
  };
  /* Create extended matrices loop over all operators */
  for (i = 0; i < NV; i++)
  {
    /* Create extended matrices on original basis */
    ExtendKrylovBasis(V[i], PsiKrylov, VPsiKrylov, NKrylov, epsilon);
    PsiPrimeOperatorPsiMatrix(V[i], VPsiKrylov, PsiKrylov, &Vtemp);
    PsiOperatorPsiMatrix(O, VPsiKrylov, &Htemp);
    /* Free VPsiKrylov */
    for (j = 0; j < NKrylov; j++)
    {
      WaveFunctionFree(&(VPsiKrylov[j]));
    }
    /* Basis transform to eigenbasis of Htemp */
    Eigensystem(Htemp, &HtempEigenFun, Hval.VRe[i]);
    MatrixMultiply(&(VMat[i]), HtempEigenFun, Vtemp);
    /* Free HtempEigenFun */
    FreeMatrix(&HtempEigenFun);
  }
  /* Free PsiKrylov */
  for (j = 0; j < NKrylov; j++)
  {
    WaveFunctionFree(&(PsiKrylov[j]));
  }
  free(PsiKrylov);
  free(VPsiKrylov);
  FreeMatrix(&Htemp);
  FreeMatrix(&Vtemp);
  /* create Gout */
  if (TridiagonalizeImpurityCoupledToBathDyson(OKrylov, VMat, Hval, g, NV, Gout))
  {
    printf("TridiagonalizeImpurityCoupledToBathDyson failed in LanczosTriDiagonalizeDyson\n");
    FreeTriDiagonalMatrix(&OKrylov);
    FreeMatrix(&Hval);
    for (i = 0; i < NV; i++)
    {
      FreeMatrix(&(VMat[i]));
    }
    free(VMat);
    return 1;
  }
  /* Free things allocated */
  FreeTriDiagonalMatrix(&OKrylov);
  FreeMatrix(&Hval);
  for (i = 0; i < NV; i++)
  {
    FreeMatrix(&(VMat[i]));
  }
  free(VMat);
  return 0;
}
