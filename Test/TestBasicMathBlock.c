void TestBasicMathBlock()
{
  int Nk[3];
  uint32_t i,j,k;
  BlockListOfPoles                      BLP,B2LP,B3LP,Sigma;
  BlockTriDiagonalMatrixType            BTD;
  BlockAndersonMatrixType               BAM;
  BlockNaturalImpurityOrbitalMatrixType BNA;
  CompactMatrixType T,psi,M1compact,TMT;
  SpectraType G1[64];
  TightBindingDef TB;
  double val[32];

  
  printf("\n============= Structure types ============\n");
  
  printf("\n============= BlockListOfPoles ============\n");
  BLP.N=4;
  BLP.BlockSize=3;
  BLP.QComplex=0;
  InitBlockListOfPoles(&BLP);
  strcpy(BLP.Name, "BLP");
  BLP.VRe[0*(1+Sqr(BLP.BlockSize))]= 1;
  BLP.VRe[1*(1+Sqr(BLP.BlockSize))]= 2;
  BLP.VRe[2*(1+Sqr(BLP.BlockSize))]=-2;
  BLP.VRe[3*(1+Sqr(BLP.BlockSize))]= 3;
  for(i=1;i<(1+Sqr(BLP.BlockSize));i++)
  {
    BLP.VRe[i]=i;
  }
  BLP.VRe[1*(1+Sqr(BLP.BlockSize))+2]=1;
  BLP.VRe[1*(1+Sqr(BLP.BlockSize))+4]=1;
  BLP.VRe[1*(1+Sqr(BLP.BlockSize))+9]=2;
  BLP.VRe[3*(1+Sqr(BLP.BlockSize))+2]=1;
  BLP.VRe[3*(1+Sqr(BLP.BlockSize))+4]=2;
  BLP.VRe[3*(1+Sqr(BLP.BlockSize))+9]=2;
  
  PrintBlockListOfPoles(BLP);
  
  printf("\n============= BlockTriDiagonalMatrixType ============\n");
  BTD.N=4;
  BTD.BlockSize=3;
  BTD.QComplex=0;
  InitBlockTridiagonalMatrix(&BTD);
  PrintBlockTriDiagonalMatrix(BTD);
  
  printf("\n============= BlockAndersonMatrixType ============\n");
  BAM.N=4;
  BAM.BlockSize=3;
  BAM.QComplex=0;
  InitBlockAndersonMatrix(&BAM);
  PrintBlockAndersonMatrix(BAM);
  
  printf("\n============= BlockNaturalImpurityOrbitalMatrixType ============\n");
  BNA.BlockSize=3;
  BNA.QComplex=0;
  BNA.val.N=2;
  BNA.con.N=3;
  InitBlockNaturalImpurityOrbitalMatrix(&BNA);
  PrintBlockNaturalImpurityOrbitalMatrix(BNA);
  
  
  printf("\n============= BlockListOfPoles Real ============\n");
  
  CopyBlockListOfPoles(BLP, &B2LP);
  OrderBlockListOfPoles(&B2LP);
  PrintBlockListOfPoles(B2LP);
  ReduceBlockListOfPoles(&B2LP, 2, 0);
  PrintBlockListOfPoles(B2LP);
  
  T.Ni=2;
  T.Nj=3;
  T.QComplex=0;
  InitCompactMatrix(&T);
  T.VRe[2]=1;
  T.VRe[3]= sqrt(0.5);
  T.VRe[4]=-sqrt(0.5);
  PrintCompactMatrixMinimal(T);
  BlockListOfPolesReduceDim(&BLP, &B3LP, &T);
  PrintBlockListOfPoles(B3LP);
  
  FreeBlockListOfPoles(&B2LP);
  FreeBlockListOfPoles(&B3LP);
  
  printf("\n============= BlockListOfPoles Complex ============\n");
  
  
  CopyBlockListOfPoles(BLP, &B2LP);
  MakeRealBlockListOfPolesComplex(&B2LP);
  
  B2LP.V[0*(1+Sqr(BLP.BlockSize))+2].Im=1;
  B2LP.V[0*(1+Sqr(BLP.BlockSize))+4].Im=1;
  B2LP.V[0*(1+Sqr(BLP.BlockSize))+9].Im=2;
  B2LP.V[1*(1+Sqr(BLP.BlockSize))+2].Im=1;
  B2LP.V[1*(1+Sqr(BLP.BlockSize))+4].Im=2;
  B2LP.V[1*(1+Sqr(BLP.BlockSize))+9].Im=2;
  
  PrintBlockListOfPoles(B2LP);
  
  OrderBlockListOfPoles(&B2LP);
  PrintBlockListOfPoles(B2LP);
  ReduceBlockListOfPoles(&B2LP, 2, 0);
  PrintBlockListOfPoles(B2LP);
  
  T.Ni=2;
  T.Nj=3;
  T.QComplex=0;
  InitCompactMatrix(&T);
  T.VRe[2]=1;
  T.VRe[3]= sqrt(0.5);
  T.VRe[4]=-sqrt(0.5);
  PrintCompactMatrixMinimal(T);
  BlockListOfPolesReduceDim(&B2LP, &B3LP, &T);
  PrintBlockListOfPoles(B3LP);
  
  FreeBlockListOfPoles(&B2LP);
  FreeBlockListOfPoles(&B3LP);
  FreeCompactMatrix(&T);
  
  printf("\n============= BlockAndersonMatrix ============\n");
  // nothing to do
  printf("\n============= BlockTridiagonalMatrix ============\n");
  // nothing to do
  printf("\n============= BlocknaturalImpurityOrbital ============\n");
  // nothing to do
  
  
  printf("\n============= BlockTypeConversions ============\n");
  
  printf("\n============= BlockListOfPolesTo... ============\n");
  
  ScaleBlockListOfPoles(&BLP, 0);
  BLP.VRe[0*(1+Sqr(BLP.BlockSize))+1]=1;
  BLP.VRe[0*(1+Sqr(BLP.BlockSize))+5]=1;
  BLP.VRe[0*(1+Sqr(BLP.BlockSize))+9]=2;
  BLP.VRe[0*(1+Sqr(BLP.BlockSize))+2]=0.5;
  BLP.VRe[0*(1+Sqr(BLP.BlockSize))+4]=0.5;
  
  BLP.VRe[1*(1+Sqr(BLP.BlockSize))+1]=1;
  BLP.VRe[1*(1+Sqr(BLP.BlockSize))+5]=1;
  BLP.VRe[1*(1+Sqr(BLP.BlockSize))+9]=2;
  BLP.VRe[1*(1+Sqr(BLP.BlockSize))+2]=-0.5;
  BLP.VRe[1*(1+Sqr(BLP.BlockSize))+4]=-0.5;
  
  BLP.VRe[2*(1+Sqr(BLP.BlockSize))+1]=1;
  BLP.VRe[2*(1+Sqr(BLP.BlockSize))+5]=1;
  BLP.VRe[2*(1+Sqr(BLP.BlockSize))+9]=2;
  BLP.VRe[2*(1+Sqr(BLP.BlockSize))+2]=0.5;
  BLP.VRe[2*(1+Sqr(BLP.BlockSize))+4]=0.5;
  BLP.VRe[2*(1+Sqr(BLP.BlockSize))+6]=1;
  BLP.VRe[2*(1+Sqr(BLP.BlockSize))+8]=1;
  
  BLP.VRe[3*(1+Sqr(BLP.BlockSize))+1]=0;
  BLP.VRe[3*(1+Sqr(BLP.BlockSize))+5]=1;
  BLP.VRe[3*(1+Sqr(BLP.BlockSize))+9]=0;
  
  printf("\n============= BlockListOfPolesToSpectra ============\n");
  
  PrintBlockListOfPoles(BLP);
  BlockListOfPolesToSpectra(BLP, G1, -5, 5, 100);
  FilePrintSpectra(G1, 9,"G1","w");
  for(i=0;i<9;i++)
  {
    SpectrumFree(&G1[i]);
  }
  
  printf("\n============= BlockListOfPolesInversedToCompactMatrix ============\n");
  
  BlockListOfPolesInversedToCompactMatrix(&BLP, &T);
  PrintCompactMatrixMinimal(T);
  FreeCompactMatrix(&T);
  FreeBlockListOfPoles(&BLP);
  
  printf("\n============= Eigensystem to compact ============\n");
  
  MatrixType M1;
  M1.QComplex=1;
  M1.Ni=5;
  M1.Nj=5;
  InitMatrix(&M1);
  for(i=0;i<5;i++)
  {
    for(j=0;j<5;j++)
    {
      M1.VRe[i][j]=sqrt((double)i+(double)j);
      M1.VIm[i][j]=(double)i-(double)j;
    }
  }
  PrintMatrix(M1);
  EigensystemToCompact(M1, &T, val);
  PrintCompactMatrixMinimal(T);
  for(i=0;i<5;i++)
  {
    printf("val %i = %22.15F\n",i,val[i]);
  }
  MatrixCopyToCompact(&M1compact, M1);
  PrintCompactMatrixMinimal(M1compact);
  CompactMatrixMultiply3('*', 'N', 'T', &T, &M1compact, &T, &TMT);
  CompactMatrixChop(&TMT, 100*DBL_EPSILON);
  PrintCompactMatrixMinimal(TMT);
  
  printf("\n============= TightBindingHamiltonianToBlockListOfPoles ============\n");
  
  
  FileReadTightBindingDef(&TB, "two_band_TB_1D_dimer_real");
  PrintTightBindingDef(TB);
  Nk[0]=5000;
  Nk[1]=1;
  Nk[2]=1;
  TightBindingHamiltonianToBlockListOfPoles(TB, Nk, &BLP, -2, 2, 4000);
  BlockListOfPolesToSpectra(BLP, G1, -2, 2, 4000);
  FilePrintSpectra(G1, 4,"G2","w");
  for(i=0;i<4;i++)
  {
    SpectrumFree(&G1[i]);
  }
  FreeBlockListOfPoles(&BLP);
  
  printf("\n============= TightBindingHamiltonianMultiSiteToBlockListOfPoles ============\n");
  int vec[9]={0,0,0,1,0,0,0,0,0};
  TightBindingHamiltonianMultiSiteToBlockListOfPoles(TB, Nk, 2, vec, &BLP, -2, 2, 4000);
  BlockListOfPolesToSpectra(BLP, G1, -2, 2, 4000);
  FilePrintSpectra(G1, 16,"G3","w");
  for(i=0;i<16;i++)
  {
    SpectrumFree(&G1[i]);
  }
  FreeBlockListOfPoles(&BLP);

  printf("\n============= TightBindingHamiltonianPlusSelfEnergyToBlockListOfPoles ============\n");
  
  Sigma.N=0;
  Sigma.BlockSize=2;
  Sigma.QComplex=0;
  InitBlockListOfPoles(&Sigma);
  Sigma.a0Re[0]=-0.25;
  Sigma.a0Re[3]= 0.25;
  TightBindingHamiltonianPlusSelfEnergyToBlockListOfPoles(TB, Nk, Sigma, &BLP, -2, 2, 4000);
  BlockListOfPolesToSpectra(BLP, G1, -2, 2, 4000);
  FilePrintSpectra(G1, 4,"G4","w");
  for(i=0;i<4;i++)
  {
    SpectrumFree(&G1[i]);
  }
  FreeBlockListOfPoles(&BLP);
  FreeBlockListOfPoles(&Sigma);
  
  Sigma.N=1;
  Sigma.BlockSize=2;
  Sigma.QComplex=0;
  InitBlockListOfPoles(&Sigma);
  Sigma.VRe[0]=0;
  Sigma.VRe[1]= 0.25;
  Sigma.VRe[4]= 0.25;
  strcpy(Sigma.Name, "Sigma");
  printf("The self energy is:\n");
  PrintBlockListOfPoles(Sigma);
  TightBindingHamiltonianPlusSelfEnergyToBlockListOfPoles(TB, Nk, Sigma, &BLP, -2, 2, 4000);
  BlockListOfPolesToSpectra(BLP, G1, -2, 2, 4000);
  FilePrintSpectra(G1, 4,"G5","w");
  for(i=0;i<4;i++)
  {
    SpectrumFree(&G1[i]);
  }
  FreeBlockListOfPoles(&BLP);
  FreeBlockListOfPoles(&Sigma);
  
  printf("\n============= TightBindingHamiltonianPlusSelfEnergyMultiSiteToBlockListOfPoles ============\n");
  
  Sigma.N=0;
  Sigma.BlockSize=2;
  Sigma.QComplex=0;
  InitBlockListOfPoles(&Sigma);
  Sigma.a0Re[ 0]=-0.25;
  Sigma.a0Re[ 3]= 0.25;
  TightBindingHamiltonianPlusSelfEnergyMultiSiteToBlockListOfPoles(TB, Nk, 2, vec, Sigma, &BLP, -2, 2, 4000);
  
  BlockListOfPolesToSpectra(BLP, G1, -2, 2, 4000);
  FilePrintSpectra(G1, 16,"G6","w");
  for(i=0;i<16;i++)
  {
    SpectrumFree(&G1[i]);
  }
  FreeBlockListOfPoles(&BLP);
  FreeBlockListOfPoles(&Sigma);

  
  printf("\n============= TightBindingHamiltonianAddSpinPlusSelfEnergyToBlockListOfPoles ============\n");
  
  Sigma.N=0;
  Sigma.BlockSize=4;
  Sigma.QComplex=0;
  InitBlockListOfPoles(&Sigma);
  Sigma.a0Re[ 2]= 0.25;
  Sigma.a0Re[ 8]= 0.25;
  Sigma.a0Re[ 7]= 0.25;
  Sigma.a0Re[13]= 0.25;
  TightBindingHamiltonianAddSpinPlusSelfEnergyToBlockListOfPoles(TB, Nk, Sigma, &BLP, -2, 2, 4000);
  
  BlockListOfPolesToSpectra(BLP, G1, -2, 2, 4000);
  FilePrintSpectra(G1, 16,"G7","w");
  for(i=0;i<16;i++)
  {
    SpectrumFree(&G1[i]);
  }
  FreeBlockListOfPoles(&BLP);
  FreeBlockListOfPoles(&Sigma);
  
  printf("\n============= TightBindingHamiltonianAddSpinPlusSelfEnergyMultiSiteToBlockListOfPoles ============\n");
  
  Sigma.N=0;
  Sigma.BlockSize=4;
  Sigma.QComplex=0;
  InitBlockListOfPoles(&Sigma);
  Sigma.a0Re[ 0]= 0.25;
  Sigma.a0Re[ 5]=-0.25;
  Sigma.a0Re[10]=-0.25;
  Sigma.a0Re[15]= 0.25;
  TightBindingHamiltonianAddSpinPlusSelfEnergyMultiSiteToBlockListOfPoles(TB, Nk, 2, vec, Sigma, &BLP, -2, 2, 4000);
  
  BlockListOfPolesToSpectra(BLP, G1, -2, 2, 4000);
  FilePrintSpectra(G1, 64,"G8","w");
  for(i=0;i<64;i++)
  {
    SpectrumFree(&G1[i]);
  }
  FreeBlockListOfPoles(&BLP);
  FreeBlockListOfPoles(&Sigma);
  
  printf("\n============= DiagonalMatrixBlockTridiagonalize ============\n");
  
  FreeBlockTridiagonalMatrix(&BTD);
  
  val[0]=1;
  val[1]=2;
  val[2]=3;
  val[3]=4;
  val[4]=5;
  val[5]=6;
  val[6]=7;
  val[7]=8;
  
  psi.Ni=3;
  psi.Nj=8;
  psi.QComplex=1;
  InitCompactMatrix(&psi);
  for(j=0;j<psi.Nj;j++)
  {
    psi.V[0*psi.Nj+j].Re=1;
    psi.V[0*psi.Nj+j].Im=1;
  }
  for(j=0;j<psi.Nj;j++)
  {
    psi.V[1*psi.Nj+j].Re=sqrt((double)j);
    psi.V[1*psi.Nj+j].Im=1;
  }
  for(j=0;j<psi.Nj;j++)
  {
    psi.V[2*psi.Nj+j].Re=psi.Nj-j;
    psi.V[2*psi.Nj+j].Im=psi.Nj-j;
  }
  strcpy(psi.Name,"beta-list");
  printf("beta-list\n");
  PrintCompactMatrixMinimal(psi);
  DiagonalMatrixBlockTridiagonalize(val, &psi, &T, &BTD, 2, DBL_EPSILON);
  
  printf("Transformation matrix\n");
  CompactMatrixChop(&T, 1000*DBL_EPSILON);
  BlockTridiagonalMatrixChop(&BTD, 1000*DBL_EPSILON);
  PrintCompactMatrixMinimal(T);
  PrintBlockTriDiagonalMatrix(BTD);
  FreeBlockTridiagonalMatrix(&BTD);
  
  
  
  
  exit(1);
  
  
  
  
  
  
  
  
  
  
  printf("\n============= BlockListOfPolesToBlockTridiagonalMatrix ============\n");
  printf("Real\n");
  BLP.N=3;
  BLP.BlockSize=4;
  BLP.QComplex=0;
  InitBlockListOfPoles(&BLP);
  BLP.VRe[0 * (1+BLP.BlockSize*BLP.BlockSize)] = -0.2;
  BLP.VRe[1 * (1+BLP.BlockSize*BLP.BlockSize)] = -0.1;
  BLP.VRe[2 * (1+BLP.BlockSize*BLP.BlockSize)] =  0.3;
  
  for(k=0;k<BLP.N;k++)
  {
    for(i=0;i<BLP.BlockSize;i++)
    {
      for(j=0;j<BLP.BlockSize;j++)
      {
        BLP.VRe[k * (1+BLP.BlockSize*BLP.BlockSize) + 1 + i + j * BLP.BlockSize] = sqrt((double)(k+i+j));
      }
      BLP.VRe[k * (1+BLP.BlockSize*BLP.BlockSize) + 1 + i * (BLP.BlockSize+1) ]+=1;
    }
  }
  PrintBlockListOfPoles(BLP);
  
  BlockListOfPolesToBlockTridiagonalMatrix(BLP, &BTD, 4, DBL_EPSILON);
  BlockTridiagonalMatrixChop(&BTD, 100*DBL_EPSILON);
  
  PrintBlockTriDiagonalMatrix(BTD);
  
  BlockTriDiagonalMatrixToCompactMatrix(&BTD, &M1compact);
  
  PrintCompactMatrixMinimal(M1compact);
  
  CompactMatrixEigensystem(&M1compact, &T, val);
  
  PrintCompactMatrixMinimal(T);
  for(i=0;i<12;i++)
  {
    printf("%22.15E\n",val[i]);
  }
  
  FreeBlockTridiagonalMatrix(&BTD);
  FreeBlockListOfPoles(&BLP);
  FreeCompactMatrix(&M1compact);
  FreeCompactMatrix(&T);
  
  printf("Complex\n");
  
  BLP.N=3;
  BLP.BlockSize=4;
  BLP.QComplex=1;
  InitBlockListOfPoles(&BLP);
  BLP.V[0 * (1+BLP.BlockSize*BLP.BlockSize)].Re = -0.2;
  BLP.V[1 * (1+BLP.BlockSize*BLP.BlockSize)].Re = -0.1;
  BLP.V[2 * (1+BLP.BlockSize*BLP.BlockSize)].Re =  0.3;
  
  for(k=0;k<BLP.N;k++)
  {
    for(i=0;i<BLP.BlockSize;i++)
    {
      for(j=0;j<BLP.BlockSize;j++)
      {
        BLP.V[k * (1+BLP.BlockSize*BLP.BlockSize) + 1 + i + j * BLP.BlockSize].Re = sqrt((double)(k+i+j));
        BLP.V[k * (1+BLP.BlockSize*BLP.BlockSize) + 1 + i + j * BLP.BlockSize].Im = ((double)i-(double)j)*0.1*sqrt((double)(k+i+j));
      }
      BLP.V[k * (1+BLP.BlockSize*BLP.BlockSize) + 1 + i * (BLP.BlockSize+1) ].Re+=2;
    }
  }
  PrintBlockListOfPoles(BLP);
  
  BlockListOfPolesToBlockTridiagonalMatrix(BLP, &BTD, 4, DBL_EPSILON);
  BlockTridiagonalMatrixChop(&BTD, 100*DBL_EPSILON);
  
  PrintBlockTriDiagonalMatrix(BTD);
  
  BlockTriDiagonalMatrixToCompactMatrix(&BTD, &M1compact);
  
  PrintCompactMatrixMinimal(M1compact);
  
  CompactMatrixEigensystem(&M1compact, &T, val);
  
  PrintCompactMatrixMinimal(T);
  for(i=0;i<12;i++)
  {
    printf("%22.15E\n",val[i]);
  }

  
  printf("\n============= BlockListOfPolesToBlockTridiagonalMatrixReduceDimension ============\n");
  
  printf("\n============= BlockListOfPolesToBlockNaturalImpurityOrbitalMatrix ============\n");
  
  printf("\n============= BlockListOfPolesToBlockNaturalImpurityOrbitalMatrixReduceDimension ============\n");
  
  printf("\n============= BlockTriDiagonalMatrixToBlockListOfPoles ============\n");
  
  printf("\n============= BlockTriDiagonalMatrixToSpectra ============\n");
  
  printf("\n============= BlockNaturalImpurityOrbitalToSpectra ============\n");
  
  printf("\n============= BlockTriDiagonalMatrixToBlockAndersonMatrix ============\n");
  
  printf("\n============= BlockListOfPolesToBlockAndersonMatrix ============\n");
  
  printf("\n============= BlockListOfPolesToBlockAndersonMatrixReduceDimension ============\n");
  
  
  
  printf("\n============= BlockMatrixToOperator ============\n");
  
  
}