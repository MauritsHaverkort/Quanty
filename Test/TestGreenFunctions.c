void GreenTransforms();
void TestSelfenergy();
void TestGBath();

void TestGreenFunctionArithmetic()
{
  printf("***************************************************************************\n");
  printf("***********************       Green Transforms       **********************\n");
  printf("***************************************************************************\n");
  GreenTransforms();
  printf("***************************************************************************\n");
  printf("***********************       Test Selfenergy        **********************\n");
  printf("***************************************************************************\n");
  TestSelfenergy();
  printf("***************************************************************************\n");
  printf("***********************       Test Bath Green        **********************\n");
  printf("***************************************************************************\n");
  TestGBath();
}

void TestGBath()
{
  uint32_t i;
  MatrixType T;
  AndersonMatrixType G0And, GAnd, GBathAnd;
  TriDiagonalMatrixType G0Tri, GTri, SigmaTri, GBathTri, GBathTriFromAnd, GBathTriFromALA;
  ListOfPoles SigmaList;
  GAnd.N=6;
  G0And.N=5;
  InitAndersonMatrix(&G0And);
  InitAndersonMatrix(&GAnd);
  GAnd.a[0]=0;
  GAnd.b[0]=1;
  for(i=1;i<GAnd.N;i++)
  {
    GAnd.a[i]=i*i-i-1.5;
    GAnd.b[i]=i*i/2.0+i/4.0;
  }
  GAnd.a[GAnd.N]=2.0;
  G0And.a[0]=0;
  G0And.b[0]=1;
  for(i=1;i<G0And.N;i++)
  {
    G0And.a[i]=GAnd.a[i];
    G0And.b[i]=GAnd.b[i]/(1+GAnd.b[i+1]);
  }
  G0And.a[G0And.N]=GAnd.a[G0And.N];
  AndersonMatrixToTriDiagonalMatrixPlusT(G0And, &T, &G0Tri, G0And.N);
  FreeMatrix(&T);
  AndersonMatrixToTriDiagonalMatrixPlusT(GAnd, &T, &GTri, GAnd.N);
  FreeMatrix(&T);
  SelfEnergyTriDiagonalRepresentation(G0Tri, GTri, &SigmaTri);
  TriDiagonalMatrixToListOfPoles(SigmaTri, &SigmaList);
  
  GBathTriDiagonalRepresentation(G0Tri, SigmaTri, &GBathTri, 30);
  
  GBathAndTriAndRepresentation(G0And, 0, SigmaTri, &GBathAnd);
  AndersonMatrixToTriDiagonalMatrixPlusT(GBathAnd, &T, &GBathTriFromAnd, GBathAnd.N);
  FreeMatrix(&T);
  FreeAndersonMatrix(&GBathAnd);
  
  GBathAndListAndRepresentation(G0And, 0, SigmaList, &GBathAnd);
  AndersonMatrixToTriDiagonalMatrixPlusT(GBathAnd, &T, &GBathTriFromALA, GBathAnd.N);
  FreeMatrix(&T);
  
  
  PrintTriDiagonalMatrix(G0Tri);
  PrintTriDiagonalMatrix(SigmaTri);
  PrintTriDiagonalMatrix(GBathTri);
  PrintTriDiagonalMatrix(GBathTriFromAnd);
  PrintTriDiagonalMatrix(GBathTriFromALA);
  
}

void TestSelfenergy()
{
  uint32_t i;
  MatrixType T;
  AndersonMatrixType G0And, GAnd;
  TriDiagonalMatrixType G0Tri, GTri, SigmaTri;
  ListOfPoles SigmaList;
  GAnd.N=6;
  G0And.N=5;
  InitAndersonMatrix(&G0And);
  InitAndersonMatrix(&GAnd);
  GAnd.a[0]=0;
  GAnd.b[0]=1;
  for(i=1;i<GAnd.N;i++)
  {
    GAnd.a[i]=i*i-i-1.5;
    GAnd.b[i]=i*i/2.0+i/4.0;
  }
  GAnd.a[GAnd.N]=2.0;
  G0And.a[0]=0;
  G0And.b[0]=1;
  for(i=1;i<G0And.N;i++)
  {
    G0And.a[i]=GAnd.a[i];
    G0And.b[i]=GAnd.b[i]/(1+GAnd.b[i+1]);
  }
  G0And.a[G0And.N]=GAnd.a[G0And.N];
  AndersonMatrixToTriDiagonalMatrixPlusT(G0And, &T, &G0Tri, G0And.N);
  FreeMatrix(&T);
  AndersonMatrixToTriDiagonalMatrixPlusT(GAnd, &T, &GTri, GAnd.N);
  FreeMatrix(&T);
  
  
  PrintAndersonMatrix(G0And);
  PrintAndersonMatrix(GAnd);
  PrintTriDiagonalMatrix(G0Tri);
  PrintTriDiagonalMatrix(GTri);
  SelfEnergyTriDiagonalRepresentation(G0Tri, GTri, &SigmaTri);
  TriDiagonalMatrixToListOfPoles(SigmaTri, &SigmaList);
  PrintTriDiagonalMatrix(SigmaTri);
  PrintListOfPoles(SigmaList);
  FreeListOfPoles(&SigmaList);
  SelfEnergyListOfPolesFromAndersonRepresentation(G0And, GAnd, &SigmaList);
  PrintListOfPoles(SigmaList);
  FreeListOfPoles(&SigmaList);
  SelfEnergyListOfPolesFromTriDiagonalRepresentation(G0Tri, GTri, &SigmaList);
  PrintListOfPoles(SigmaList);
  printf("\n");
  for(i=0;i<SigmaList.N;i++)
  {
    printf("%i a=%22.15E b=%22.15E sqrt(b)=%22.15E\n",i,SigmaList.a[i],SigmaList.b[i],sqrt(SigmaList.b[i]));
  }
  FreeListOfPoles(&SigmaList);
}

void GreenTransforms()
{
  uint32_t i;
  TriDiagonalMatrixType MTri,MTriBack;
  AndersonMatrixType MAnd;
  ListOfPoles MList;
  DoubleTriDiagonalMatrixType MDoubleTri;
  SecondDoubleTriDiagonalMatrixType MSecDoubleTri;
  MatrixType T;
  MTri.N=7;
  InitTriDiagonalMatrix(&MTri);
  for(i=0;i<MTri.N;i++)
  {
    MTri.a[i]=(i+1.)/2-i*i/4.;
    MTri.b[i]=(i+1.)/(2i+2);
  }
  MTri.a[MTri.N]=MTri.N-1;
  TriDiagonalMatrixToAndersonMatrixPlusT(MTri, &T, &MAnd);
  FreeMatrix(&T);
  TriDiagonalMatrixToListOfPoles(MTri, &MList);

  TriDiagonalMatrixToDoubleTriDiagonalMatrixPlusT(MTri, &T, &MDoubleTri, 0.0, 0.0);
  FreeMatrix(&T);
  TriDiagonalMatrixToSecondDoubleTriDiagonalMatrix(MTri, &T, &MSecDoubleTri, 0.0, 0.0);
  FreeMatrix(&T);
  
  printf("===============\n");
  PrintTriDiagonalMatrix(MTri);
  printf("===============\n");
  PrintAndersonMatrix(MAnd);
  printf("===============\n");
  PrintListOfPoles(MList);
  printf("===============\n");
  PrintDoubleTriDiagonalMatrix(MDoubleTri);
  printf("===============\n");
  PrintSecondDoubleTriDiagonalMatrix(MSecDoubleTri);
  printf("===============\n");
  
  AndersonMatrixToTriDiagonalMatrixPlusT(MAnd, &T, &MTriBack, MAnd.N);
  printf("=============== Back from And ============\n");
  PrintTriDiagonalMatrix(MTriBack);
  FreeMatrix(&T);
  FreeTriDiagonalMatrix(&MTriBack);
  
  ListOfPolesToTridiagonalMatrix(MList, &MTriBack, MList.N);
  printf("=============== Back from List ============\n");
  PrintTriDiagonalMatrix(MTriBack);
  FreeTriDiagonalMatrix(&MTriBack);
  
  DoubleTriDiagonalMatrixToTriDiagonalMatrix(MDoubleTri, &MTriBack, 7);
  printf("=============== Back from DoubleTri ============\n");
  PrintTriDiagonalMatrix(MTriBack);
  FreeTriDiagonalMatrix(&MTriBack);
  
  SecondDoubleTriDiagonalMatrixToTriDiagonalMatrix(MSecDoubleTri, &MTriBack, 7);
  printf("=============== Back from Second DoubleTri ============\n");
  PrintTriDiagonalMatrix(MTriBack);
  FreeTriDiagonalMatrix(&MTriBack);
  
  
}



void GreenSpec()
{
  
  
  uint32_t i,Dim=10;
  TriDiagonalMatrixType At,Bt,Ct,ApBt,BpCt,ApCt,Apt,Bpt,Cpt;
  SpectraType Spec[15];
  Spec[0].N=4000;
  Spec[0].Emin=-10.0;
  Spec[0].Emax= 10.0;
  Spec[0].EminPole=-10.0;
  Spec[0].EmaxPole= 10.0;
  Spec[0].dE=20.0/4000.0;
  Spec[0].Gamma=0.1;
  Spec[0].phase=0.0;
  for(i=0;i<15;i++)
  {
    SpectrumAllocate(Spec[0], &(Spec[i]));
  }
  At.N=Dim;
  Bt.N=Dim;
  Ct.N=Dim;
  InitTriDiagonalMatrix(&At);
  InitTriDiagonalMatrix(&Bt);
  InitTriDiagonalMatrix(&Ct);
  for(i=0;i<Dim;i++)
  {
    At.a[i+1]=(double)i;
    At.b[i  ]=sqrt(1.0/(double)(i+1));
    Bt.a[i+1]=-(double)i;
    Bt.b[i  ]=sqrt(1.0/(double)(i+1));
    Ct.a[i+1]=0;
    Ct.b[i  ]=sqrt(1.0/(double)(i+1));
  }
  At.a[0]=0;
  Bt.a[0]=0;
  Ct.a[0]=0;
  SumGreensfunctionsTriDiagonalRepresentation(At, Bt, &ApBt);
  SumGreensfunctionsTriDiagonalRepresentation(Bt, Ct, &BpCt);
  SumGreensfunctionsTriDiagonalRepresentation(At, Ct, &ApCt);
  DifferenceGreensFunctionTriDiagonalRepresentation(ApBt, Bt, &Apt);
  DifferenceGreensFunctionTriDiagonalRepresentation(BpCt, Ct, &Bpt);
  DifferenceGreensFunctionTriDiagonalRepresentation(ApCt, At, &Cpt);
  TriDiagonalMatrixToSpectrum(0.0, At, &Spec[0]);
  TriDiagonalMatrixToSpectrum(0.0, Bt, &Spec[1]);
  TriDiagonalMatrixToSpectrum(0.0, Ct, &Spec[2]);
  TriDiagonalMatrixToSpectrum(0.0, ApBt, &Spec[3]);
  TriDiagonalMatrixToSpectrum(0.0, BpCt, &Spec[4]);
  TriDiagonalMatrixToSpectrum(0.0, ApCt, &Spec[5]);
  TriDiagonalMatrixToSpectrum(0.0, Apt, &Spec[9]);
  TriDiagonalMatrixToSpectrum(0.0, Bpt, &Spec[10]);
  TriDiagonalMatrixToSpectrum(0.0, Cpt, &Spec[11]);
  CopySpectra(Spec[0], &Spec[6]);
  CopySpectra(Spec[1], &Spec[7]);
  CopySpectra(Spec[2], &Spec[8]);
  AddMultipliedSpectra(&Spec[6], Spec[1], 1.0, 0.0);
  AddMultipliedSpectra(&Spec[7], Spec[2], 1.0, 0.0);
  AddMultipliedSpectra(&Spec[8], Spec[0], 1.0, 0.0);
  FreeTriDiagonalMatrix(&Apt);
  FreeTriDiagonalMatrix(&Bpt);
  FreeTriDiagonalMatrix(&Cpt);
  ApCt.N=15;
  DifferenceGreensFunctionTriDiagonalRepresentation(ApCt, Ct, &Apt);
  ApCt.N=10;
  DifferenceGreensFunctionTriDiagonalRepresentation(ApCt, Ct, &Bpt);
  ApCt.N=5;
  DifferenceGreensFunctionTriDiagonalRepresentation(ApCt, Ct, &Cpt);
  TriDiagonalMatrixToSpectrum(0.0, Apt, &Spec[12]);
  TriDiagonalMatrixToSpectrum(0.0, Bpt, &Spec[13]);
  TriDiagonalMatrixToSpectrum(0.0, Cpt, &Spec[14]);
  
  FilePrintSpectra(Spec,15,"TestSpec_01","w");
  
}

