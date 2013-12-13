void TestGFBA()
{
  uint32_t i, NBath = 4;
  double epsilon, PsiEnergy, PsiQuality, PsiHist[20];
  double *PsiRho;
  TriDiagonalMatrixType G0Tri;
  BADoubleTriDiagonalMatrixType G0BADoubleTri;
  WaveFunctionType Psi;
  OperatorType H;
  printf(" ************** Test G0 to Hamiltonian to G0 *************\n");
  if (FileReadPlainTriDiagonalMatrixType(&G0Tri, "G0_2D"))
  {
    printf("Fatal Error in main\n Failed to read G0 from file\n");
    return;
  };
  strcpy(G0Tri.Name, "G0");
  for (i = 1; i <= G0Tri.N; i++)
  {
    G0Tri.a[i] += 0.1;
  }
  TriDiagonalMatrixToBADoubleTriDiagonalMatrix(G0Tri, &G0BADoubleTri, NBath + 1, NBath + 1, 0.0, 0.0000001);
  PrintTriDiagonalMatrix(G0Tri);
  PrintBADoubleTriDiagonalMatrix(G0BADoubleTri);
  BADoubleTriDiagonalMatrixToAndersonImpurityOperator(&H, G0BADoubleTri, G0BADoubleTri, 0.0, 0.0);
  CreateWaveFunctionDMFTSecond(&Psi, NBath, NBath, 16);

  PrintOperatorReadable(H);
  PrintWaveFunctionReadable(Psi);

  epsilon = 1E-12;
  LanczosGroundState(H, &Psi, &epsilon, 1E-12, 0, &Psi, 0);

  PrintWaveFunctionReadableTruncated(Psi, 0.01);


  PsiHermitianOperatorPsi(H, Psi, &PsiEnergy);
  PsiHermitianOperatorSqrPsi(H, Psi, &PsiQuality);
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
    return;
  }
  CalculateDensity(PsiRho, Psi);
  printf("Occupation \n");
  for (i = 0; i < Psi.NFermion; i++)
  {
    printf("%5i %22.15E\n", i, PsiRho[i]);
  }
  free(PsiRho);
}
