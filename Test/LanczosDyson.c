/*
 *  LanczosDyson.c
 *  Model
 *
 *  Created by Maurits Haverkort on 3/23/12.
 *
 */

void TestLanczosDyson3()
{
  TightBindingDef TB;
  SpectraType DOS,Gspec[1];
  BlockListOfPoles G;
  int Nk[]={200,200,1};
  // create G0
  FileReadTightBindingDef(&TB, "two_band_TB");
  PrintTightBindingDef(TB);
  DOS.N=100;
  DOS.Emin=-4;
  DOS.Emax=6;
  DOS.dE=(DOS.Emax-DOS.Emin)/DOS.N;
  SpectrumAllocate(DOS, &DOS);
  TightBindingDOS(TB, Nk, &DOS);
  
  FilePrintSpectra(&DOS, 1, "TBDOS", "w");
  
  TightBindingHamiltonianToBlockListOfPoles(TB, Nk, &G, -3, 5.5, 1000);
  
  BlockListOfPolesToSpectra(G, Gspec, -4, 6, 100);
  
  PrintBlockListOfPoles(G);
  
  FilePrintSpectra(Gspec, TB.NFermion*TB.NFermion, "TBDOS2", "w");
  
  
  // turn G0 into impurity natural orbitals
  
  // no spin PES and IPES
  
  // with spin PES and IPES
  
}


void TestLanczosDyson2()
{
  printf("Start Test Lanczos Dyson2\n");

  MatrixType T, RhoMat;
  TriDiagonalMatrixType G0, Gp[2], GPES, GIPES, GPESD1, GIPESD1, GPESD2, GIPESD2;
  SecondDoubleTriDiagonalMatrixType G0Double;
  OperatorType V[2], O, TPES, TIPES;
  WaveFunctionType psi, Opsi;
  uint32_t i, NTriDiag = 10, NCellV = 4, NCellC = 4, NKrylov = 150;
  uint16_t oadd[2];
  uint8_t *Determinant;
  double PsiEnergy, PsiQuality, PsiHist[20], *PsiRho;
  double epsilon, zero;
  SpectraType Spectra[10];
  G0.N = 150;
  InitTriDiagonalMatrix(&G0);
  for (i = 1; i < G0.N; i++)
  {
    G0.a[i] = 0.0;
    G0.b[i] = 1.0;
  }
  G0.a[G0.N] = 0.0;
  TriDiagonalMatrixToSecondDoubleTriDiagonalMatrixPlusT(G0, &T, &G0Double, 0.0, 0.00000001);
  FreeMatrix(&T);
  PrintTriDiagonalMatrix(G0);
  PrintSecondDoubleTriDiagonalMatrix(G0Double);

  InitOperator(&O);
  O.NFermion = 2 + NCellV + NCellC;
  O.NBoson = 0;
  O.NWord = Min(O.NFermion, 1 + (O.NFermion - 1) / 8);
  O.QComplex = 0;
  InitOperatorFixedLength(&O, 2, 8);

  oadd[1] = CrAnMask + 0;
  oadd[0] = 0;
  OperatorAddLaderNormalOrder(&O, oadd, 2, G0Double.Dt.a0, 0.0, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 1;
  oadd[0] = 1;
  OperatorAddLaderNormalOrder(&O, oadd, 2, (G0Double.Dt.val.a[0] + G0Double.Dt.con.a[0]) / 2, 0.0,
                              &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 1;
  oadd[0] = 0;
  OperatorAddLaderNormalOrder(&O, oadd, 2, G0Double.Dt.b0, 0.0, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 0;
  oadd[0] = 1;
  OperatorAddLaderNormalOrder(&O, oadd, 2, G0Double.Dt.b0, 0.0, &QDetNotOrderedFcFaBcBa);

  oadd[1] = CrAnMask + 2;
  oadd[0] = 2;
  OperatorAddLaderNormalOrder(&O, oadd, 2, G0Double.Dt.val.a[1], 0.0, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 2;
  oadd[0] = 1;
  OperatorAddLaderNormalOrder(&O, oadd, 2, G0Double.Dt.val.b[0], 0.0, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 1;
  oadd[0] = 2;
  OperatorAddLaderNormalOrder(&O, oadd, 2, G0Double.Dt.val.b[0], 0.0, &QDetNotOrderedFcFaBcBa);

  for (i = 1; i < NCellV; i++)
  {
    oadd[1] = CrAnMask + i + 2;
    oadd[0] = i + 2;
    OperatorAddLaderNormalOrder(&O, oadd, 2, G0Double.Dt.val.a[i + 1], 0.0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i + 2;
    oadd[0] = i + 1;
    OperatorAddLaderNormalOrder(&O, oadd, 2, G0Double.Dt.val.b[i], 0.0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i + 1;
    oadd[0] = i + 2;
    OperatorAddLaderNormalOrder(&O, oadd, 2, G0Double.Dt.val.b[i], 0.0, &QDetNotOrderedFcFaBcBa);
  }

  oadd[1] = CrAnMask + 2 + NCellV;
  oadd[0] = 2 + NCellV;
  OperatorAddLaderNormalOrder(&O, oadd, 2, G0Double.Dt.con.a[1], 0.0, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 2 + NCellV;
  oadd[0] = 1;
  OperatorAddLaderNormalOrder(&O, oadd, 2, G0Double.Dt.con.b[0], 0.0, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 1;
  oadd[0] = 2 + NCellV;
  OperatorAddLaderNormalOrder(&O, oadd, 2, G0Double.Dt.con.b[0], 0.0, &QDetNotOrderedFcFaBcBa);

  for (i = 1; i < NCellC; i++)
  {
    oadd[1] = CrAnMask + i + 2 + NCellV;
    oadd[0] = i + 2 + NCellV;
    OperatorAddLaderNormalOrder(&O, oadd, 2, G0Double.Dt.con.a[i + 1], 0.0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i + 2 + NCellV;
    oadd[0] = i + 1 + NCellV;
    OperatorAddLaderNormalOrder(&O, oadd, 2, G0Double.Dt.con.b[i], 0.0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i + 1 + NCellV;
    oadd[0] = i + 2 + NCellV;
    OperatorAddLaderNormalOrder(&O, oadd, 2, G0Double.Dt.con.b[i], 0.0, &QDetNotOrderedFcFaBcBa);
  }

  PrintOperatorReadable(O);

  for (i = 0; i < 2; i++)
  {
    InitOperator(&V[i]);
    V[i].NFermion = 2 + NCellV + NCellC;
    V[i].NBoson = 0;
    V[i].NWord = Min(V[i].NFermion, 1 + (V[i].NFermion - 1) / 8);
    V[i].QComplex = 0;
    InitOperatorFixedLength(&V[i], 1, 1);
  }
  oadd[0] = CrAnMask + 1 + NCellV;
  OperatorAddLaderNormalOrder(&V[0], oadd, 1, G0Double.Dt.val.b[NCellV], 0.0, &QDetNotOrderedFcFaBcBa);
  oadd[0] = 1 + NCellV + NCellC;
  OperatorAddLaderNormalOrder(&V[1], oadd, 1, G0Double.Dt.con.b[NCellC], 0.0, &QDetNotOrderedFcFaBcBa);

  PrintOperatorReadable(V[0]);
  PrintOperatorReadable(V[1]);

  Gp[0].N = G0Double.Dt.val.N - NCellV;
  InitTriDiagonalMatrix(&Gp[0]);
  for (i = 1; i < Gp[0].N; i++)
  {
    Gp[0].a[i] = -G0Double.Dt.val.a[i + NCellV];
    Gp[0].b[i] = G0Double.Dt.val.b[i + NCellV];
  }
  Gp[0].a[Gp[0].N] = -G0Double.Dt.val.a[Gp[0].N + NCellV];

  Gp[1].N = G0Double.Dt.con.N - NCellC;
  InitTriDiagonalMatrix(&Gp[1]);
  for (i = 1; i < Gp[1].N; i++)
  {
    Gp[1].a[i] = G0Double.Dt.con.a[i + NCellC];
    Gp[1].b[i] = G0Double.Dt.con.b[i + NCellC];
  }
  Gp[1].a[Gp[1].N] = G0Double.Dt.con.a[Gp[1].N + NCellC];

  PrintTriDiagonalMatrix(Gp[0]);
  PrintTriDiagonalMatrix(Gp[1]);

  InitOperator(&TPES);
  TPES.NFermion = 2 + NCellV + NCellC;
  TPES.NBoson = 0;
  TPES.NWord = Min(TPES.NFermion, 1 + (TPES.NFermion - 1) / 8);
  TPES.QComplex = 0;
  InitOperatorFixedLength(&TPES, 1, 1);
  oadd[0] = 0;
  RealOperatorAddLader(&TPES, oadd, 1, 1.0);

  InitOperator(&TIPES);
  TIPES.NFermion = 2 + NCellV + NCellC;
  TIPES.NBoson = 0;
  TIPES.NWord = Min(TIPES.NFermion, 1 + (TIPES.NFermion - 1) / 8);
  TIPES.QComplex = 0;
  InitOperatorFixedLength(&TIPES, 1, 1);
  oadd[0] = CrAnMask + 0;
  RealOperatorAddLader(&TIPES, oadd, 1, 1.0);

  PrintOperatorReadable(TPES);
  PrintOperatorReadable(TIPES);

  psi.QComplex = 0;
  psi.NFermion = 2 + NCellV + NCellC;
  psi.NBoson = 0;
  psi.NBitsKey = 16;
  WaveFunctionInit(&psi);
  Determinant = calloc(psi.NWord, sizeof(uint8_t));
  for (i = 0; i < psi.NWord; i++)
  {
    Determinant[i] = 0;
  }
  for (i = 1; i < 2 + NCellV; i++)
  {
    OperatorCr(Determinant, i, psi.NFermion);
  }
  RealWaveFunctionAddElement(&psi, Determinant, 1.0);
  PrintWaveFunctionFullReadable(psi);

  epsilon = 1E-15;
  zero = 1E-15;
  LanczosGroundState(O, &psi, &epsilon, zero, 0.0, &psi, 0);
  LanczosGroundState(O, &psi, &epsilon, zero, 0.0, &psi, 0);






  PrintWaveFunctionReadableTruncated(psi, 0.1);
  PsiHermitianOperatorPsi(O, psi, &PsiEnergy);
  PsiHermitianOperatorSqrPsi(O, psi, &PsiQuality);
  PsiQuality /= Sqr(PsiEnergy);
  CalculatePrefactorHistogram(PsiHist, psi, 10.0, 20);
  printf("\nHistogram of prefactors in Wave-function\n");
  for (i = 0; i < 19; i++)
  {
    printf("10^(-%2i) to 10^(-%2i) | %22.15E\n", i, i + 1, PsiHist[i]);
  }
  printf("10^(-%2i) to 0        | %22.15E\n", i, PsiHist[19]);
  printf("\n");
  printf("Energy                  <psi|H^2|psi>/(<psi|H|psi>^2)\n");
  printf("%22.15E  %22.15E\n", PsiEnergy, PsiQuality);
  PsiRho = calloc(psi.NFermion, sizeof(double));
  if (PsiRho == NULL)
  {
    printf("calloc failed in main for PsiRho\n Fatal Error\n");
    fflush(stdout);
    return;
  }
  CalculateDensity(PsiRho, psi);
  printf("Occupation \n");
  for (i = 0; i < psi.NFermion; i++)
  {
    printf("%5i %22.15E\n", i, PsiRho[i]);
  }
  free(PsiRho);

  CalculateDensityCorrelationMatrix(&RhoMat, psi);
  PrintMatrix(RhoMat);



  fflush(stdout);

  GPES.N = NTriDiag;
  GIPES.N = NTriDiag;
  GPESD1.N = NKrylov;
  GIPESD1.N = NKrylov;
  GPESD2.N = NKrylov;
  GIPESD2.N = NKrylov;
  InitTriDiagonalMatrix(&GPES);
  InitTriDiagonalMatrix(&GIPES);
  InitTriDiagonalMatrix(&GPESD1);
  InitTriDiagonalMatrix(&GIPESD1);

  OperatorPsi(TPES, psi, &Opsi);
  LanczosTriDiagonalize(O, Opsi, epsilon, &GPES, 0);

  OperatorPsi(TIPES, psi, &Opsi);
  LanczosTriDiagonalize(O, Opsi, epsilon, &GIPES, 0);

  OperatorPsi(TPES, psi, &Opsi);
  LanczosTriDiagonalizeDysonLowMem(O, Opsi, epsilon, NTriDiag, V, Gp, 2, &GPESD1);

  OperatorPsi(TIPES, psi, &Opsi);
  LanczosTriDiagonalizeDysonLowMem(O, Opsi, epsilon, NTriDiag, V, Gp, 2, &GIPESD1);

  OperatorPsi(TPES, psi, &Opsi);
  LanczosTriDiagonalizeDysonFullLowMem(O, 0, Opsi, epsilon, NTriDiag, V, Gp, 2, &GPESD2);

  OperatorPsi(TIPES, psi, &Opsi);
  LanczosTriDiagonalizeDysonFullLowMem(O, 0, Opsi, epsilon, NTriDiag, V, Gp, 2, &GIPESD2);


  printf("End Tridiagonalization\n");

  fflush(stdout);

  Spectra[0].N = 6000;
  Spectra[0].Emin = -3;
  Spectra[0].Emax = 3;
  Spectra[0].dE = (Spectra[0].Emax - Spectra[0].Emin) / Spectra[0].N;
  Spectra[0].EminPole = -1;
  Spectra[0].EmaxPole = 1;
  Spectra[0].Gamma = 0.08;
  for (i = 0; i < 11; i++)
  {
    SpectrumAllocate(Spectra[0], &Spectra[i]);
  }

  TriDiagonalMatrixToSpectrum(PsiEnergy, GPES, &Spectra[0]);
  TriDiagonalMatrixToSpectrum(PsiEnergy, GIPES, &Spectra[1]);
  TriDiagonalMatrixToSpectrum(PsiEnergy, GPESD1, &Spectra[2]);
  TriDiagonalMatrixToSpectrum(PsiEnergy, GIPESD1, &Spectra[3]);
  TriDiagonalMatrixToSpectrum(PsiEnergy, GPESD2, &Spectra[4]);
  TriDiagonalMatrixToSpectrum(PsiEnergy, GIPESD2, &Spectra[5]);
  ShiftTriDiagonalMatrix(&GPESD1, -PsiEnergy);
  ShiftTriDiagonalMatrix(&GIPESD1, -PsiEnergy);
  ShiftTriDiagonalMatrix(&GPESD2, -PsiEnergy);
  ShiftTriDiagonalMatrix(&GIPESD2, -PsiEnergy);
  for (i = 1; i <= GPESD1.N; i++)
  {
    GPESD1.a[i] *= -1;
  }
  for (i = 1; i <= GPESD2.N; i++)
  {
    GPESD2.a[i] *= -1;
  }
  TriDiagonalMatrixToSpectrum(0, GPESD1, &Spectra[6]);
  TriDiagonalMatrixToSpectrum(0, GIPESD1, &Spectra[7]);
  TriDiagonalMatrixToSpectrum(0, GPESD2, &Spectra[8]);
  TriDiagonalMatrixToSpectrum(0, GIPESD2, &Spectra[9]);
  TriDiagonalMatrixToSpectrum(0, G0, &Spectra[10]);
  FilePrintSpectra(Spectra, 11, "/Users/haverkort/ProgramOutput/LancDyson2", "w");
}


void TestLanczosDyson()
{
  printf("Start Test Lanczos Dyson\n");

  uint32_t i;
  uint16_t oadd[6];
  uint32_t NCell = 5, NCelltot = 200, NTriDiag = 5;
  OperatorType O, TPES, h12, h12Dag;
  double epsilon, zero, ZeroEneryShift = -30;
  double *OMatrixExtendedEigenVal;
  WaveFunctionType psi, Opsi, *PsiKrylov, *PsiKrylovExtended;
  uint8_t *Determinant;
  double PsiEnergy, PsiQuality, PsiHist[20], *PsiRho;
  MatrixType Rho, h12matrix, h12matrixEigenOExtended, OMatrixExdended, OMatrixExtendedEigenFun;
  MatrixType test1, test2, test3, OExtendedEigenVal, V[1];
  TriDiagonalMatrixType G, Gp[1], GKrylov, GwithBath, GwithBathSec, GwithBathSecLowMem;
  SpectraType Spectra[6];

  psi.QComplex = 0;
  psi.NFermion = NCell;
  psi.NBoson = 0;
  psi.NBitsKey = 16;
  WaveFunctionInit(&psi);
  Determinant = calloc(psi.NWord, sizeof(uint8_t));
  for (i = 0; i < psi.NWord; i++)
  {
    Determinant[i] = 0;
  }
  OperatorCr(Determinant, 0, psi.NFermion);
  OperatorCr(Determinant, 1, psi.NFermion);
  OperatorCr(Determinant, 2, psi.NFermion);
  RealWaveFunctionAddElement(&psi, Determinant, 1.0);
  PrintWaveFunctionFullReadable(psi);

  Rho.QComplex = 0;
  Rho.Nj = NCell;
  Rho.Ni = NTriDiag;
  InitMatrix(&Rho);

  InitOperator(&TPES);
  TPES.NFermion = NCell;
  TPES.NBoson = 0;
  TPES.NWord = Min(TPES.NFermion, 1 + (TPES.NFermion - 1) / 8);
  TPES.QComplex = 0;
  InitOperatorFixedLength(&TPES, 1, 1);
  oadd[0] = 1;
  RealOperatorAddLader(&TPES, oadd, 1, 1.0);

  InitOperator(&h12);
  h12.NFermion = NCell;
  h12.NBoson = 0;
  h12.NWord = Min(h12.NFermion, 1 + (h12.NFermion - 1) / 8);
  h12.QComplex = 0;
  InitOperatorFixedLength(&h12, 1, 1);
  oadd[0] = NCell - 1;
  RealOperatorAddLader(&h12, oadd, 1, 0.25);

  InitOperator(&O);
  O.NFermion = NCell;
  O.NBoson = 0;
  O.NWord = Min(O.NFermion, 1 + (O.NFermion - 1) / 8);
  O.QComplex = 0;
  InitOperatorFixedLength(&O, 2, 8);
  InitOperatorFixedLength(&O, 4, 8);
  for (i = 1; i < NCell; i++)
  {
    oadd[1] = CrAnMask + i - 1;
    oadd[0] = i;
    OperatorAddLaderNormalOrder(&O, oadd, 2, 0.25, 0.0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i;
    oadd[0] = i - 1;
    OperatorAddLaderNormalOrder(&O, oadd, 2, 0.25, 0.0, &QDetNotOrderedFcFaBcBa);
  }
  oadd[3] = CrAnMask + 0;
  oadd[2] = 0;
  oadd[1] = CrAnMask + 1;
  oadd[0] = 1;
  OperatorAddLaderNormalOrder(&O, oadd, 4, -1.0, 0, &QDetNotOrderedFcFaBcBa);
//  oadd[3]=CrAnMask+0  ;  oadd[2]=0  ;     oadd[1]=CrAnMask+2;   oadd[0]=2;     OperatorAddLaderNormalOrder(&O, oadd, 4,-1.0, 0, &QDetNotOrderedFcFaBcBa);
  oadd[3] = CrAnMask + 1;
  oadd[2] = 1;
  oadd[1] = CrAnMask + 2;
  oadd[0] = 2;
  OperatorAddLaderNormalOrder(&O, oadd, 4, -0.5, 0, &QDetNotOrderedFcFaBcBa);

//  oadd[1]=CrAnMask+0  ;  oadd[0]=0  ;     OperatorAddLaderNormalOrder(&O, oadd, 2,-0.5, 0.0, &QDetNotOrderedFcFaBcBa);
//  oadd[1]=CrAnMask+1  ;  oadd[0]=1  ;     OperatorAddLaderNormalOrder(&O, oadd, 2,-0.5, 0.0, &QDetNotOrderedFcFaBcBa);
//  oadd[1]=CrAnMask+2  ;  oadd[0]=2  ;     OperatorAddLaderNormalOrder(&O, oadd, 2,-0.5, 0.0, &QDetNotOrderedFcFaBcBa);
  OperatorAddLader(&O, oadd, 0, ZeroEneryShift, 0);
  SetQComplexOperator(&O);

  epsilon = 1E-15;
  zero = 1E-15;
  LanczosGroundState(O, &psi, &epsilon, zero, 0.0, &psi, 0);
  LanczosGroundState(O, &psi, &epsilon, zero, 0.0, &psi, 0);

  PrintWaveFunctionReadableTruncated(psi, 0.1);
  PsiHermitianOperatorPsi(O, psi, &PsiEnergy);
  PsiHermitianOperatorSqrPsi(O, psi, &PsiQuality);
  PsiQuality /= Sqr(PsiEnergy);
  CalculatePrefactorHistogram(PsiHist, psi, 10.0, 20);
  printf("\nHistogram of prefactors in Wave-function\n");
  for (i = 0; i < 19; i++)
  {
    printf("10^(-%2i) to 10^(-%2i) | %22.15E\n", i, i + 1, PsiHist[i]);
  }
  printf("10^(-%2i) to 0        | %22.15E\n", i, PsiHist[19]);
  printf("\n");
  printf("Energy                  <psi|H^2|psi>/(<psi|H|psi>^2)\n");
  printf("%22.15E  %22.15E\n", PsiEnergy, PsiQuality);
  PsiRho = calloc(psi.NFermion, sizeof(double));
  if (PsiRho == NULL)
  {
    printf("calloc failed in main for PsiRho\n Fatal Error\n");
    fflush(stdout);
    return;
  }
  CalculateDensity(PsiRho, psi);
  printf("Occupation \n");
  for (i = 0; i < psi.NFermion; i++)
  {
    printf("%5i %22.15E\n", i, PsiRho[i]);
  }
  free(PsiRho);

  fflush(stdout);

  printf("Start normal Lanczos\n");

  G.N = NTriDiag;
  InitTriDiagonalMatrix(&G);
  OperatorPsi(TPES, psi, &Opsi);
  PrintWaveFunctionReadableTruncated(Opsi, 0.1);
  LanczosTriDiagonalizePlusDensity(O, Opsi, epsilon, &G, &Rho, epsilon);

  fflush(stdout);

  printf("Start Krylov Tridiagonalization\n");

  GKrylov.N = NTriDiag;

  /* Allocate memory for *PsiBasis */
  if ((PsiKrylov = calloc(NTriDiag, sizeof(WaveFunctionType))) == NULL)
  {
    printf("calloc failed in main\n");
    fflush(stdout);
    return;
  };

  InitTriDiagonalMatrix(&GKrylov);
  OperatorPsi(TPES, psi, &Opsi);
  PrintWaveFunctionReadableTruncated(Opsi, 0.0);
  LanczosTriDiagonalizeKrylov(O, Opsi, epsilon, &GKrylov, PsiKrylov);

  printf("End Krylov Tridiagonalization\n");

  fflush(stdout);

  Spectra[0].N = 5000;
  Spectra[0].Emin = -1;
  Spectra[0].Emax = 4;
  Spectra[0].dE = (Spectra[0].Emax - Spectra[0].Emin) / Spectra[0].N;
  Spectra[0].EminPole = -1;
  Spectra[0].EmaxPole = 1;
  Spectra[0].Gamma = 0.01;
  for (i = 0; i < 6; i++)
  {
    SpectrumAllocate(Spectra[0], &Spectra[i]);
  }


  printf("Start creation of aditional basis given by h12 PsiKrylov\n");


  if ((PsiKrylovExtended = calloc(NTriDiag, sizeof(WaveFunctionType))) == NULL)
  {
    printf("calloc failed in main\n");
    fflush(stdout);
    return;
  };

  ExtendKrylovBasis(h12, PsiKrylov, PsiKrylovExtended, NTriDiag, epsilon);

  printf("Start calculation of H on extended basis\n");

  h12matrix.Ni = NTriDiag;
  h12matrix.Nj = NTriDiag;
  h12matrix.QComplex = 0;
  InitMatrix(&h12matrix);
  OMatrixExdended.Ni = NTriDiag;
  OMatrixExdended.Nj = NTriDiag;
  OMatrixExdended.QComplex = 0;
  InitMatrix(&OMatrixExdended);

  PsiPrimeOperatorPsiMatrix(h12, PsiKrylovExtended, PsiKrylov, &h12matrix);

  // for(i=0;i<5;i++) PrintWaveFunctionPlain(PsiKrylovExtended[i]);
  // PrintOperatorFullReadable(O);


  PsiOperatorPsiMatrix(O, PsiKrylovExtended, &OMatrixExdended);

  OMatrixExtendedEigenVal = malloc(sizeof(double) * NTriDiag);

  Eigensystem(OMatrixExdended, &OMatrixExtendedEigenFun, OMatrixExtendedEigenVal);

  printf("Eigenfun and eigenval OExtended\n");
  for (i = 0; i < NTriDiag; i++)
  {
    printf("%lf ", OMatrixExtendedEigenVal[i]);
  }
  printf("\n");
  PrintMatrix(OMatrixExtendedEigenFun);

  printf("\n");

  MatrixMultiply(&h12matrixEigenOExtended, OMatrixExtendedEigenFun, h12matrix);

  MatrixCopy(&test1, OMatrixExtendedEigenFun);
  MatrixConjugateTranspose(&test1);
  MatrixMultiply(&test2, OMatrixExdended, test1);
  MatrixMultiply(&test3, OMatrixExtendedEigenFun, test2);
  PrintMatrix(test3);
  PrintMatrix(OMatrixExdended);


  PrintMatrix(h12matrix);


  PrintMatrix(h12matrixEigenOExtended);



  printf("Start creation of Gp\n");

  Gp[0].N = NCelltot - NCell;
  InitTriDiagonalMatrix(&Gp[0]);
  for (i = 0; i < Gp[0].N; i++)
  {
    Gp[0].a[i + 1] = 0;
    Gp[0].b[i] = 0.25;
  }
  Gp[0].a[0] = 0;


  printf("Start creation of GwithBath\n");

  MatrixCopy(&(V[0]), h12matrixEigenOExtended);


  GwithBath.N = 200;
  OExtendedEigenVal.Ni = 1;
  OExtendedEigenVal.Nj = OMatrixExtendedEigenFun.Ni;
  OExtendedEigenVal.QComplex = 0;
  InitMatrix(&OExtendedEigenVal);
  for (i = 0; i < OMatrixExtendedEigenFun.Ni; i++)
  {
    OExtendedEigenVal.VRe[0][i] = OMatrixExtendedEigenVal[i];
  }
  TridiagonalizeImpurityCoupledToBathDyson(GKrylov, V, OExtendedEigenVal, Gp, 1, &GwithBath);


  OperatorPsi(TPES, psi, &Opsi);
  GwithBathSec.N = 200;
  LanczosTriDiagonalizeDyson(O, Opsi, epsilon, NTriDiag, &h12, Gp, 1, &GwithBathSec);

  OperatorTranspose(h12, &h12Dag);
  OperatorPsi(TPES, psi, &Opsi);
  GwithBathSecLowMem.N = 200;
  LanczosTriDiagonalizeDysonLowMem(O, Opsi, epsilon, NTriDiag, &h12, Gp, 1, &GwithBathSecLowMem);


  printf("Start creation of Spectra\n");

  TriDiagonalMatrixToSpectrum(0, Gp[0], &Spectra[0]);
  TriDiagonalMatrixToSpectrum(PsiEnergy, G, &Spectra[1]);
  TriDiagonalMatrixToSpectrum(PsiEnergy, GKrylov, &Spectra[2]);
  TriDiagonalMatrixToSpectrum(PsiEnergy, GwithBath, &Spectra[3]);
  TriDiagonalMatrixToSpectrum(PsiEnergy, GwithBathSec, &Spectra[4]);
  TriDiagonalMatrixToSpectrum(PsiEnergy, GwithBathSecLowMem, &Spectra[5]);

  FilePrintSpectra(Spectra, 6, "/Users/haverkort/ProgramOutput/LancDyson", "w");

  printf("End of Test Lanczos Dyson\n");
}


void TestLanczosDysonOld()
{
  printf("Start Test Lanczos Dyson\n");

  uint32_t i;
  uint16_t oadd[6];
  uint32_t NCell = 10, NCelltot = 400;
  OperatorType O, TPES, Odens, OplusG;
  double epsilon, zero, ZeroEneryShift = 0;
  WaveFunctionType psi, Opsi, *PsiKrylov, *PsiKrylovAr;
  uint8_t *Determinant;
  double PsiEnergy, PsiQuality, PsiHist[20], *PsiRho;
  MatrixType Rho, RhoBound, RhoMatrix, GKrylovNonHer, HonKrylovBasis;
  TriDiagonalMatrixType G, Gp[1], Gnew, GKrylov;
  SpectraType Spectra[6];
  psi.QComplex = 0;
  psi.NFermion = NCell;
  psi.NBoson = 0;
  psi.NBitsKey = 16;
  if (WaveFunctionInit(&psi))
  {
    printf("WaveFunctionInit failed in TestLanczosDyson\n");
    fflush(stdout);
    return;
  };
  Determinant = calloc(psi.NWord, sizeof(uint8_t));
  for (i = 0; i < psi.NWord; i++)
  {
    Determinant[i] = 0;
  }
  OperatorCr(Determinant, 0, psi.NFermion);
  OperatorCr(Determinant, 1, psi.NFermion);
//  OperatorCr(Determinant,2,psi.NFermion);
  RealWaveFunctionAddElement(&psi, Determinant, 1.0);
  PrintWaveFunctionFullReadable(psi);

  Rho.QComplex = 0;
  Rho.Nj = NCell;
  Rho.Ni = 50;
  InitMatrix(&Rho);

  InitOperator(&TPES);
  TPES.NFermion = NCell;
  TPES.NBoson = 0;
  TPES.NWord = Min(TPES.NFermion, 1 + (TPES.NFermion - 1) / 8);
  TPES.QComplex = 0;
  InitOperatorFixedLength(&TPES, 1, 1);
  oadd[0] = CrAnMask + 0;
  RealOperatorAddLader(&TPES, oadd, 1, 1.0);

  InitOperator(&Odens);
  Odens.NFermion = NCell;
  Odens.NBoson = 0;
  Odens.NWord = Min(Odens.NFermion, 1 + (Odens.NFermion - 1) / 8);
  Odens.QComplex = 0;
  InitOperatorFixedLength(&Odens, 2, 1);
  oadd[1] = CrAnMask + NCell - 1;
  oadd[0] = NCell - 1;
  RealOperatorAddLader(&Odens, oadd, 2, 1.0);

  InitOperator(&O);
  O.NFermion = NCell;
  O.NBoson = 0;
  O.NWord = Min(O.NFermion, 1 + (O.NFermion - 1) / 8);
  O.QComplex = 0;
  InitOperatorFixedLength(&O, 2, 8);
  InitOperatorFixedLength(&O, 4, 8);
  for (i = 1; i < NCell; i++)
  {
    oadd[1] = CrAnMask + i - 1;
    oadd[0] = i;
    OperatorAddLaderNormalOrder(&O, oadd, 2, 0.25, 0.0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i;
    oadd[0] = i - 1;
    OperatorAddLaderNormalOrder(&O, oadd, 2, 0.25, 0.0, &QDetNotOrderedFcFaBcBa);
  }
//  oadd[3]=CrAnMask+0  ;  oadd[2]=0  ;     oadd[1]=CrAnMask+1;   oadd[0]=1;     OperatorAddLaderNormalOrder(&O, oadd, 4,-0.5, 0, &QDetNotOrderedFcFaBcBa);
//  oadd[3]=CrAnMask+0  ;  oadd[2]=0  ;     oadd[1]=CrAnMask+2;   oadd[0]=2;     OperatorAddLaderNormalOrder(&O, oadd, 4,-0.5, 0, &QDetNotOrderedFcFaBcBa);
//  oadd[3]=CrAnMask+1  ;  oadd[2]=1  ;     oadd[1]=CrAnMask+2;   oadd[0]=2;     OperatorAddLaderNormalOrder(&O, oadd, 4,-0.5, 0, &QDetNotOrderedFcFaBcBa);

  oadd[1] = CrAnMask + 0;
  oadd[0] = 0;
  OperatorAddLaderNormalOrder(&O, oadd, 2, -0.25, 0.0, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 1;
  oadd[0] = 1;
  OperatorAddLaderNormalOrder(&O, oadd, 2, -0.25, 0.0, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 2;
  oadd[0] = 2;
  OperatorAddLaderNormalOrder(&O, oadd, 2, -0.25, 0.0, &QDetNotOrderedFcFaBcBa);
  OperatorAddLader(&O, oadd, 0, ZeroEneryShift, 0);
  SetQComplexOperator(&O);
  PrintOperatorFullReadable(O);
  PrintWaveFunctionFullReadable(psi);

  PrintWaveFunctionReadableTruncated(psi, 0.1);
  PsiHermitianOperatorPsi(O, psi, &PsiEnergy);
  PsiHermitianOperatorSqrPsi(O, psi, &PsiQuality);
  PsiQuality /= Sqr(PsiEnergy);
  CalculatePrefactorHistogram(PsiHist, psi, 10.0, 20);
  printf("\nHistogram of prefactors in Wave-function\n");
  for (i = 0; i < 19; i++)
  {
    printf("10^(-%2i) to 10^(-%2i) | %22.15E\n", i, i + 1, PsiHist[i]);
  }
  printf("10^(-%2i) to 0        | %22.15E\n", i, PsiHist[19]);
  printf("\n");
  printf("Energy                  <psi|H^2|psi>/(<psi|H|psi>^2)\n");
  printf("%22.15E  %22.15E\n", PsiEnergy, PsiQuality);
  PsiRho = calloc(psi.NFermion, sizeof(double));
  if (PsiRho == NULL)
  {
    printf("calloc failed in main for PsiRho\n Fatal Error\n");
    fflush(stdout);
    return;
  }
  CalculateDensity(PsiRho, psi);
  printf("Occupation \n");
  for (i = 0; i < psi.NFermion; i++)
  {
    printf("%5i %22.15E\n", i, PsiRho[i]);
  }
  free(PsiRho);

  printf("Hello world A\n");


  epsilon = 0.000000001;
  zero = 0.000000001;
  LanczosGroundState(O, &psi, &epsilon, zero, 0.0, &psi, 0);
  LanczosGroundState(O, &psi, &epsilon, zero, 0.0, &psi, 0);

  printf("Hello world B\n");

  PrintWaveFunctionReadableTruncated(psi, 0.1);
  PsiHermitianOperatorPsi(O, psi, &PsiEnergy);
  PsiHermitianOperatorSqrPsi(O, psi, &PsiQuality);
  PsiQuality /= Sqr(PsiEnergy);
  CalculatePrefactorHistogram(PsiHist, psi, 10.0, 20);
  printf("\nHistogram of prefactors in Wave-function\n");
  for (i = 0; i < 19; i++)
  {
    printf("10^(-%2i) to 10^(-%2i) | %22.15E\n", i, i + 1, PsiHist[i]);
  }
  printf("10^(-%2i) to 0        | %22.15E\n", i, PsiHist[19]);
  printf("\n");
  printf("Energy                  <psi|H^2|psi>/(<psi|H|psi>^2)\n");
  printf("%22.15E  %22.15E\n", PsiEnergy, PsiQuality);
  PsiRho = calloc(psi.NFermion, sizeof(double));
  if (PsiRho == NULL)
  {
    printf("calloc failed in main for PsiRho\n Fatal Error\n");
    fflush(stdout);
    return;
  }
  CalculateDensity(PsiRho, psi);
  printf("Occupation \n");
  for (i = 0; i < psi.NFermion; i++)
  {
    printf("%5i %22.15E\n", i, PsiRho[i]);
  }
  free(PsiRho);

  fflush(stdout);

  G.N = 15;
  InitTriDiagonalMatrix(&G);
  OperatorPsi(TPES, psi, &Opsi);
  PrintWaveFunctionReadableTruncated(Opsi, 0.1);
  LanczosTriDiagonalizePlusDensity(O, Opsi, epsilon, &G, &Rho, epsilon);


  printf("Start Krylov Tridiagonalization\n");
  fflush(stdout);

  GKrylov.N = G.N;

  /* Allocate memory for *PsiBasis */
  if ((PsiKrylov = calloc(GKrylov.N, sizeof(WaveFunctionType))) == NULL)
  {
    printf("calloc failed in main\n");
    fflush(stdout);
    return;
  };

  InitTriDiagonalMatrix(&GKrylov);
  OperatorPsi(TPES, psi, &Opsi);
  PrintWaveFunctionReadableTruncated(Opsi, 0.0);
  LanczosTriDiagonalizeKrylov(O, Opsi, epsilon, &GKrylov, PsiKrylov);



  /* for(i=0;i<GKrylov.N;i++)
     {
     PrintWaveFunctionReadableTruncated(PsiKrylov[i], 0.0);
     } */


  fflush(stdout);

  RhoMatrix.Ni = GKrylov.N;
  RhoMatrix.Nj = GKrylov.N;
  RhoMatrix.QComplex = 0;
  InitMatrix(&RhoMatrix);

  PsiOperatorPsiMatrix(Odens, PsiKrylov, &RhoMatrix);

  PrintMatrix(RhoMatrix);

  printf("End Krylov Tridiagonalization\n");


  printf("Start Krylov Tridiagonalization Arnoldi\n");
  fflush(stdout);

  GKrylovNonHer.Ni = G.N;
  GKrylovNonHer.Nj = G.N;
  GKrylovNonHer.QComplex = 1;
  InitMatrix(&GKrylovNonHer);


  InitOperator(&OplusG);
  OplusG.NFermion = NCell;
  OplusG.NBoson = 0;
  OplusG.NWord = Min(OplusG.NFermion, 1 + (OplusG.NFermion - 1) / 8);
  OplusG.QComplex = 0;
  OperatorAdd(&OplusG, O, 1, 0);
  OperatorAdd(&OplusG, Odens, 0, 1);

  PrintOperatorReadable(OplusG);

  /* Allocate memory for *PsiBasis */
  if ((PsiKrylovAr = calloc(GKrylovNonHer.Ni, sizeof(WaveFunctionType))) == NULL)
  {
    printf("calloc failed in main\n");
    fflush(stdout);
    return;
  };


  OperatorPsi(TPES, psi, &Opsi);

  /*MakeComplexWaveFunctionFromReal(&COpsi, Opsi);
     ComplexWaveFunctionComplexMultiply(COpsi, sqrt(0.5), sqrt(0.5)); */

  PrintWaveFunctionReadableTruncated(Opsi, 0.0);
  ArnoldiTriDiagonalizeKrylov(OplusG, Opsi, epsilon, &GKrylovNonHer, PsiKrylovAr);

  PrintMatrix(GKrylovNonHer);

  printf("End Krylov Tridiagonalization Arnoldi\n");


  PrintMatrix(Rho);

  Spectra[0].N = 5000;
  Spectra[0].Emin = -1;
  Spectra[0].Emax = 4;
  Spectra[0].dE = (Spectra[0].Emax - Spectra[0].Emin) / Spectra[0].N;
  Spectra[0].EminPole = -1;
  Spectra[0].EmaxPole = 1;
  Spectra[0].Gamma = 0.01;
  for (i = 0; i < 6; i++)
  {
    SpectrumAllocate(Spectra[0], &Spectra[i]);
  }

  RhoBound.QComplex = 0;
  RhoBound.Nj = 1;
  RhoBound.Ni = 50;
  InitMatrix(&RhoBound);

  for (i = 0; i < Rho.Ni; i++)
  {
    RhoBound.VRe[0][i] = Rho.VRe[i][NCell - 1];
  }

  PrintMatrix(RhoBound);


  Gp[0].N = NCelltot - NCell;
  InitTriDiagonalMatrix(&Gp[0]);
  for (i = 0; i < Gp[0].N; i++)
  {
    Gp[0].a[i + 1] = ZeroEneryShift;
    Gp[0].b[i] = 0.25;
  }
  Gp[0].a[0] = 0;

  printf("Hello world C\n");

  TriDiagonalMatrixDysonScalar(G, &Gnew, RhoBound, 1, Gp);



  TriDiagonalMatrixToSpectrum(0 + 0 * PsiEnergy, Gp[0], &Spectra[0]);
  TriDiagonalMatrixToSpectrum(PsiEnergy, G, &Spectra[1]);
  TriDiagonalMatrixToSpectrum(PsiEnergy, Gnew, &Spectra[2]);


  TriDiagonalMatrixToMatrix(GKrylov, &HonKrylovBasis);

  PrintMatrix(HonKrylovBasis);
  PrintMatrix(RhoMatrix);

  MatrixToSpectrum(HonKrylovBasis, PsiEnergy, Spectra[0].Gamma, &RhoMatrix, &(Spectra[0]), 1, &(Spectra[3]));



  FilePrintSpectra(Spectra, 4, "/Users/haverkort/ProgramOutput/LancDyson", "w");


}
