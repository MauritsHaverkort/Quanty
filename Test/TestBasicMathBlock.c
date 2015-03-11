void TestBasicMathBlock()
{
  uint32_t i;
  BlockListOfPoles                      BLP,B2LP,B3LP;
  BlockTriDiagonalMatrixType            BTD;
  BlockAndersonMatrixType               BAM;
  BlockNaturalImpurityOrbitalMatrixType BNA;
  CompactMatrixType T;
  
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
  
  
  printf("\n============= BlockListOfPoles ============\n");
  
  OrderBlockListOfPoles(&BLP);
  PrintBlockListOfPoles(BLP);
  CopyBlockListOfPoles(BLP, &B2LP);
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
  
  
  printf("\n============= BlockAndersonMatrix ============\n");
  printf("\n============= BlockTridiagonalMatrix ============\n");
  printf("\n============= BlocknaturalImpurityOrbital ============\n");
  printf("\n============= BlockTypeConversions ============\n");
  printf("\n============= BlockMatrixToOperator ============\n");
  
}