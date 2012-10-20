/*
 *  LanczosDyson.c
 *  Model
 *
 *  Created by Maurits Haverkort on 3/23/12.
 *  Copyright 2012 Max Planck Institut. All rights reserved.
 *
 */


void TestLanczosDyson()
{
  printf("Start Test Lanczos Dyson\n");
  
  uint32_t i;
  uint16_t oadd[6];
  uint32_t NCell=5,NCelltot=300;
  OperatorType O,TPES;
  double epsilon,zero;
  WaveFunctionType psi,Opsi;
  uint8_t *Determinant;
  double PsiEnergy,PsiQuality,PsiHist[20],*PsiRho;  
  MatrixType Rho,RhoBound;
  TriDiagonalMatrixType G,Gp[1],Gnew;
  SpectraType Spectra[2];
  psi.QComplex=0;
  psi.NFermion=NCell;
  psi.NBoson=0;
  psi.NBitsKey=16;
  if(WaveFunctionInit(&psi)){printf("WaveFunctionInit failed in TestLanczosDyson\n");fflush(stdout);return;};
  Determinant=calloc(psi.NWord,sizeof(uint8_t));
  for(i=0;i<psi.NWord;i++)
  {
    Determinant[i]=0;
  }
  OperatorCr(Determinant,0,psi.NFermion);
  OperatorCr(Determinant,1,psi.NFermion);
  RealWaveFunctionAddElement(&psi,Determinant,1.0);
  PrintWaveFunctionFullReadable(psi);
  Rho.QComplex=0;
  Rho.Nj=NCell;
  Rho.Ni=50;
  InitMatrix(&Rho);
  
  InitOperator(&TPES);
  TPES.NFermion=NCell;
  TPES.NBoson=0;
  TPES.NWord=Min(TPES.NFermion,1+(TPES.NFermion-1)/8);
  TPES.QComplex=0;
  InitOperatorFixedLength(&TPES,1,1);
  oadd[0]=0;
  RealOperatorAddLader(&TPES, oadd, 1, 1.0);
  
  InitOperator(&O);
  O.NFermion=NCell;
  O.NBoson=0;
  O.NWord=Min(O.NFermion,1+(O.NFermion-1)/8);
  O.QComplex=0;
  InitOperatorFixedLength(&O,2,8);
  InitOperatorFixedLength(&O,4,8);
  for(i=1;i<NCell;i++)
  {
    oadd[1]=CrAnMask+i-1;     oadd[0]=i  ;     OperatorAddLaderNormalOrder(&O, oadd, 2, 0.25, 0.0, &QDetNotOrderedFcFaBcBa);
    oadd[1]=CrAnMask+i  ;     oadd[0]=i-1;     OperatorAddLaderNormalOrder(&O, oadd, 2, 0.25, 0.0, &QDetNotOrderedFcFaBcBa);
  }
  oadd[3]=CrAnMask+0  ;  oadd[2]=0  ;     oadd[1]=CrAnMask+1;   oadd[0]=1;     OperatorAddLaderNormalOrder(&O, oadd, 4,-1.0, 0, &QDetNotOrderedFcFaBcBa);
  oadd[1]=CrAnMask+0  ;  oadd[0]=0  ;     OperatorAddLaderNormalOrder(&O, oadd, 2,-0.0, 0.0, &QDetNotOrderedFcFaBcBa);
  //OperatorAddLader(&O, oadd, 0, -20.0, 0);
  SetQComplexOperator(&O);
  PrintOperatorFullReadable(O);
  PrintWaveFunctionFullReadable(psi);
  
  PrintWaveFunctionReadableTruncated(psi, 0.1);
  PsiHermitianOperatorPsi(O, psi, &PsiEnergy);
  PsiHermitianOperatorSqrPsi(O,psi,&PsiQuality);
  PsiQuality/=Sqr(PsiEnergy);
  CalculatePrefactorHistogram(PsiHist, psi, 10.0, 20);
  printf("\nHistogram of prefactors in Wave-function\n");
  for(i=0;i<19;i++)
  {
    printf("10^(-%2i) to 10^(-%2i) | %22.15E\n",i,i+1,PsiHist[i]);
  }
  printf("10^(-%2i) to 0        | %22.15E\n",i,PsiHist[19]);
  printf("\n");
  printf("Energy                  <psi|H^2|psi>/(<psi|H|psi>^2)\n");
  printf("%22.15E  %22.15E\n",PsiEnergy,PsiQuality);
  PsiRho=calloc(psi.NFermion,sizeof(double));
  if(PsiRho==NULL){printf("calloc failed in main for PsiRho\n Fatal Error\n");fflush(stdout);return;}
  CalculateDensity(PsiRho,psi);
  printf("Occupation \n");
  for(i=0;i<psi.NFermion;i++)
  {
    printf("%5i %22.15E\n",i,PsiRho[i]);
  }
  free(PsiRho);
  
  printf("Hello world A\n");
  
  
  epsilon=0.000001;
  zero=   0.000001;
  LanczosGroundState(O, &psi, &epsilon, zero, 0.0, &psi, 0);  
  
  printf("Hello world B\n");
  
  PrintWaveFunctionReadableTruncated(psi, 0.1);
    PsiHermitianOperatorPsi(O, psi, &PsiEnergy);
    PsiHermitianOperatorSqrPsi(O,psi,&PsiQuality);
    PsiQuality/=Sqr(PsiEnergy);
  CalculatePrefactorHistogram(PsiHist, psi, 10.0, 20);
  printf("\nHistogram of prefactors in Wave-function\n");
  for(i=0;i<19;i++)
  {
    printf("10^(-%2i) to 10^(-%2i) | %22.15E\n",i,i+1,PsiHist[i]);
  }
  printf("10^(-%2i) to 0        | %22.15E\n",i,PsiHist[19]);
  printf("\n");
  printf("Energy                  <psi|H^2|psi>/(<psi|H|psi>^2)\n");
  printf("%22.15E  %22.15E\n",PsiEnergy,PsiQuality);
  PsiRho=calloc(psi.NFermion,sizeof(double));
  if(PsiRho==NULL){printf("calloc failed in main for PsiRho\n Fatal Error\n");fflush(stdout);return;}
  CalculateDensity(PsiRho,psi);
  printf("Occupation \n");
  for(i=0;i<psi.NFermion;i++)
  {
    printf("%5i %22.15E\n",i,PsiRho[i]);
  }
  free(PsiRho);
  
  fflush(stdout);
  
  G.N=20;
  InitTriDiagonalMatrix(&G);
  OperatorPsi(TPES,psi,&Opsi);
  LanczosTriDiagonalizePlusDensity(O, Opsi, epsilon,&G,&Rho, epsilon);
   
  PrintMatrix(Rho);
  
    Spectra[0].N       =2000;
    Spectra[0].Emin    =-1;
    Spectra[0].Emax    =1;
    Spectra[0].dE      =(Spectra[0].Emax-Spectra[0].Emin)/Spectra[0].N;
    Spectra[0].EminPole=-1;
    Spectra[0].EmaxPole=1;
    Spectra[0].Gamma   =0.01;
    for(i=0;i<2;i++)
    {
      SpectrumAllocate(Spectra[0], &Spectra[i]);
    }
  
  RhoBound.QComplex=0;
  RhoBound.Nj=1;
  RhoBound.Ni=50;
  InitMatrix(&RhoBound);
  
  for(i=0;i<Rho.Ni;i++)
  {
    RhoBound.VRe[0][i]=Rho.VRe[NCell][i];
  }
  
  Gp[0].N=NCelltot-NCell;
  InitTriDiagonalMatrix(&Gp[0]);
  for(i=0;i<Gp[0].N;i++)
  {
    Gp[0].a[i]=0;
    Gp[0].b[i]=0.25;
  }
  Gp[0].a[Gp[0].N]=0;
  
  TriDiagonalMatrixDysonScalar(G, &Gnew, RhoBound, 1, Gp);
  
  
  
  TriDiagonalMatrixToSpectrum(0, G, &Spectra[0]);
  TriDiagonalMatrixToSpectrum(0, Gnew, &Spectra[1]);
  FilePrintSpectra(Spectra, 2, "LancDyson", "w");
  
  
}