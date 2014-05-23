void TestOppMult()
{
  printf("Test multiplication of operators\n");
  OperatorType A, B, AB, BA;
  WaveFunctionType psi, Apsi, Bpsi, ABpsi, BApsi;
  uint8_t *Determinant;
  uint16_t oadd[10];
  psi.NBitsKey = 8;
  psi.QComplex = 0;
  psi.NBoson = 0;
  psi.NFermion = 10;
  if (WaveFunctionInit(&psi))
  {
    printf("WaveFunctionInit failed in CreateWaveFunctionDMFT\n");
    fflush(stdout);
    return;
  };
  Determinant = calloc(psi.NWord, sizeof(uint8_t *));
  sprintf(psi.Name, "psi test");
  OperatorCr(Determinant, 0, psi.NFermion);
  OperatorCr(Determinant, 1, psi.NFermion);
  RealWaveFunctionAddElement(&psi, Determinant, 1.0);
  free(Determinant);
  PrintWaveFunctionReadable(psi);

  InitNamedOperator(&A, "OppA");
  InitNamedOperator(&B, "OppB");
  InitNamedOperator(&BA, "OppBA");
  InitNamedOperator(&AB, "OppAB");
  A.NFermion = psi.NFermion;
  A.NBoson = psi.NBoson;
  A.NWord = Min(A.NFermion, 1 + (A.NFermion - 1) / 8);
  B.NFermion = psi.NFermion;
  B.NBoson = psi.NBoson;
  B.NWord = Min(A.NFermion, 1 + (B.NFermion - 1) / 8);
  AB.NFermion = psi.NFermion;
  AB.NBoson = psi.NBoson;
  AB.NWord = Min(AB.NFermion, 1 + (AB.NFermion - 1) / 8);
  BA.NFermion = psi.NFermion;
  BA.NBoson = psi.NBoson;
  BA.NWord = Min(A.NFermion, 1 + (BA.NFermion - 1) / 8);
  InitOperatorFixedLength(&A, 1, 8);
  InitOperatorFixedLength(&A, 2, 8);
  InitOperatorFixedLength(&B, 1, 8);
  InitOperatorFixedLength(&B, 2, 8);
  InitOperatorFixedLength(&AB, 1, 8);
  InitOperatorFixedLength(&AB, 2, 8);
  InitOperatorFixedLength(&AB, 3, 8);
  InitOperatorFixedLength(&AB, 4, 8);
  InitOperatorFixedLength(&BA, 1, 8);
  InitOperatorFixedLength(&BA, 2, 8);
  InitOperatorFixedLength(&BA, 3, 8);
  InitOperatorFixedLength(&BA, 4, 8);


  oadd[0] = CrAnMask + 0;
  oadd[1] = CrAnMask + 2;
  OperatorAddLader(&A, oadd, 2, 1.0, 0.0);


  oadd[0] = 1;
  oadd[1] = 0;
  OperatorAddLader(&B, oadd, 2, 1.0, 0.0);



  OperatorMultiply2NormalOrder(&AB, A, B, 'N', 'N', &QDetNotOrderedFcFaBcBa);
  OperatorMultiply2NormalOrder(&BA, B, A, 'N', 'N', &QDetNotOrderedFcFaBcBa);

  OperatorPsi(A, psi, &Apsi);
  OperatorPsi(B, psi, &Bpsi);
  OperatorPsi(AB, psi, &ABpsi);
  OperatorPsi(BA, psi, &BApsi);

  PrintOperatorReadable(A);
  PrintOperatorReadable(B);
  PrintOperatorReadable(AB);
  PrintOperatorReadable(BA);

  PrintWaveFunctionReadable(psi);
  PrintWaveFunctionReadable(Apsi);
  PrintWaveFunctionReadable(Bpsi);
  PrintWaveFunctionReadable(ABpsi);
  PrintWaveFunctionReadable(BApsi);

  printf("*********************************\n\n\n\n\n\n\n\n");





  OperatorType Sx, Sy, Sz, S2, Stest;
  uint16_t SpinUp[6], SpinDn[6];
  SpinUp[0] = 0;
  SpinUp[1] = 2;
  SpinUp[2] = 4;
  SpinDn[0] = 1;
  SpinDn[1] = 3;
  SpinDn[2] = 5;
  CreateSpinOperators(6, 3, SpinUp, SpinDn, &Sx, &Sy, &Sz, &S2);
  PrintOperatorReadable(Sx);
  PrintOperatorReadable(Sy);
  PrintOperatorReadable(Sz);
  PrintOperatorReadable(S2);

  InitOperatorCopySettings(&Stest, S2);
  //OperatorMultiplyNormalOrder(&Stest, Sy, Sx,&QDetNotOrderedFcFaBcBa);
  OperatorMultiply2NormalOrder(&Stest, Sy, Sx, 'N', 'N', &QDetNotOrderedFcFaBcBa);
  OperatorScale(&Stest, -1.0, 0.0);
  //OperatorMultiplyNormalOrder(&Stest, Sx, Sy,&QDetNotOrderedFcFaBcBa);
  OperatorMultiply2NormalOrder(&Stest, Sx, Sy, 'N', 'N', &QDetNotOrderedFcFaBcBa);
  strcpy(Stest.Name, "i Sz = Sx.Sy - Sy.Sx");
  PrintOperatorReadable(Stest);
  FreeOperator(&Stest);

  InitOperatorCopySettings(&Stest, S2);
  OperatorMultiplyNormalOrder(&Stest, Sz, Sy, &QDetNotOrderedFcFaBcBa);
  OperatorScale(&Stest, -1.0, 0.0);
  OperatorMultiplyNormalOrder(&Stest, Sy, Sz, &QDetNotOrderedFcFaBcBa);
  strcpy(Stest.Name, "i Sx = Sy.Sz - Sz.Sy");
  PrintOperatorReadable(Stest);
  FreeOperator(&Stest);

  InitOperatorCopySettings(&Stest, S2);
  OperatorMultiplyNormalOrder(&Stest, Sx, Sz, &QDetNotOrderedFcFaBcBa);
  OperatorScale(&Stest, -1.0, 0.0);
  OperatorMultiplyNormalOrder(&Stest, Sz, Sx, &QDetNotOrderedFcFaBcBa);
  strcpy(Stest.Name, "i Sy = Sz.Sx - Sx.Sz");
  PrintOperatorReadable(Stest);
  FreeOperator(&Stest);

  printf("Test mult\n");

  InitOperatorCopySettings(&Stest, S2);
  OperatorMultiplyNormalOrder(&Stest, Sx, Sz, &QDetNotOrderedFcFaBcBa);
  PrintOperatorReadable(Stest);
  OperatorScale(&Stest, 0.0, 1.0);
  PrintOperatorReadable(Stest);


}


void TestRotateAtomic()
{
  uint32_t i;
  OperatorType H,U,Sx,Sy,Sz,Ssqr,Lx,Ly,Lz,Lsqr,Jx,Jy,Jz,Jsqr,ldots;
  WaveFunctionType psi, Rotpsi, Rpsi;
  double epsilon=0.00000001, zero=100*DBL_EPSILON, RhoDiag[14];
  MatrixType Rho, EigenRho;
  OperatorType RH;
  uint32_t NFermion=14,L=3,Nelectron=3;
  uint16_t LindexUp[]={1,3,5,7,9,11,13};
  uint16_t LindexDn[]={0,2,4,6,8,10,12};
  double F[]={0,10.0,6.0,4.0};
  CreateOperatorCoulombSingleShell(NFermion, L, LindexUp, LindexDn, F, &U);
  CreateOperatorSx(NFermion, 2*L+1, LindexUp, LindexDn, &Sx);
  CreateOperatorSy(NFermion, 2*L+1, LindexUp, LindexDn, &Sy);
  CreateOperatorSz(NFermion, 2*L+1, LindexUp, LindexDn, &Sz);
  CreateOperatorSsqr(NFermion, 2*L+1, LindexUp, LindexDn, &Ssqr);
  CreateOperatorLx(NFermion, L, LindexUp, LindexDn, &Lx);
  CreateOperatorLy(NFermion, L, LindexUp, LindexDn, &Ly);
  CreateOperatorLz(NFermion, L, LindexUp, LindexDn, &Lz);
  CreateOperatorLsqr(NFermion, L, LindexUp, LindexDn, &Lsqr);
  CreateOperatorJx(NFermion, L, LindexUp, LindexDn, &Jx);
  CreateOperatorJy(NFermion, L, LindexUp, LindexDn, &Jy);
  CreateOperatorJz(NFermion, L, LindexUp, LindexDn, &Jz);
  CreateOperatorJsqr(NFermion, L, LindexUp, LindexDn, &Jsqr);
  CreateOperatorldots(NFermion, L, LindexUp, LindexDn, &ldots);
  OperatorCopy(&H, U);
  OperatorAdd(&H, ldots, 0.01, 0.0);
  OperatorAdd(&H, Sx, 0.00002, 0.0);
  OperatorAdd(&H, Lx, 0.00001, 0.0);
  CreateStartWaveFunctionOnsiteEnergy(H, &psi, Nelectron, 0, 16);
  LanczosGroundState(H, &psi, &epsilon, zero, 0, &psi, 0);
  printf("\n================\n Operator \n================\n ");
  PrintOperatorReadable(H);
  printf("\n================\n Groundstate \n================\n ");
  PrintWaveFunctionReadable(psi);
  printf("\n================\n Density Matrix \n================\n ");
  CalculateDensityMatrix(&Rho, psi);
  PrintMatrix(Rho);
  Eigensystem(Rho, &EigenRho, RhoDiag);
  PrintMatrix(EigenRho);
  if(EigenRho.QComplex==0)
  {
    MatrixTranspose(&EigenRho);
  }
  else
  {
    MatrixConjugateTranspose(&EigenRho);
  }

  for(i=0;i<NFermion;i++)
  {
    printf("Density of orbital %i: %lf\n",i,RhoDiag[i]);
  }
  
  printf("\n================\n Rotate psi \n================\n ");
  RotateWaveFunction(&Rotpsi, psi, EigenRho);
  PrintWaveFunctionReadable(Rotpsi);
  
  
  printf("\n================\n Rotate Operator \n================\n ");
  RotateOperatorNormalOrder(&RH, H, EigenRho, &QDetNotOrderedFcFaBcBa);
  OperatorChop(&RH, 0.000000001);
  PrintOperatorReadable(RH);
  
  
  printf("\n================\n Eigensystem of RH \n================\n ");
  WaveFunctionCopy(Rotpsi, &Rpsi);
  LanczosGroundState(RH, &Rpsi, &epsilon, zero, 0, &Rotpsi, 0);
  
  PrintWaveFunctionReadable(Rpsi);
  
}



void TestRotate()
{
  
  
  TestRotateAtomic();
  return;
  
  MatrixType Rho, RotMat, RhoRot, Rho2;
  OperatorType Opp, OppRot;
  WaveFunctionType Psi, PsiRot, PsiRot2;
  uint32_t NFermion=10,QComplex=1,NElectron=5,NBitsKey=16;
  uint32_t i;
  double epsilon=0.0000001,zero=DBL_EPSILON;
  double RhoVal[NFermion];
  CreateRandomTwoParticleOperator(NFermion, QComplex, &Opp);
  
  PrintOperatorReadable(Opp);
  CreateStartWaveFunctionOnsiteEnergy(Opp, &Psi, NElectron, 0, NBitsKey);
  
  printf("==========\nGroundstate on original basis\n==========\n\n");
  
  LanczosGroundState(Opp, &Psi, &epsilon, zero, 0, &Psi, 0);
  
  PrintWaveFunctionReadableTruncated(Psi, 0.01);
  
  CalculateDensityMatrix(&Rho, Psi);
  
  
  
  PrintMatrix(Rho);
  
  Eigensystem(Rho, &RotMat, RhoVal);
  
  
  PrintMatrix(RotMat);
  for(i=0;i<NFermion;i++)
  {
    printf("Density is %lf\n",RhoVal[i]);
  }
  
  printf("==========\nRotate operator\n==========\n\n");
  
  
  if(RotMat.QComplex==0)
  {
    MatrixTranspose(&RotMat);
  }
  else
  {
    MatrixConjugateTranspose(&RotMat);
  }
  
  
  RotateOperator(&OppRot, Opp, RotMat);
  
  PrintOperatorReadable(OppRot);
  
  printf("==========\nGroundstte rotated operator\n==========\n\n");
  
  CreateStartWaveFunctionOnsiteEnergy(OppRot, &PsiRot, NElectron, 0, NBitsKey);
  
  
  LanczosGroundState(OppRot, &PsiRot, &epsilon, zero, 0, &PsiRot, 0);
  
  PrintWaveFunctionReadableTruncated(PsiRot, 0.01);
  
  CalculateDensityMatrix(&RhoRot, PsiRot);
  
  PrintMatrix(RhoRot);
  
  
  printf("==========\nRotate original eigenstate\n==========\n\n");
  
  RotateWaveFunction(&PsiRot2, Psi, RotMat);
  
  PrintWaveFunctionReadable(PsiRot2);
  
  CalculateDensityMatrix(&Rho2, PsiRot2);
  
  PrintMatrix(Rho2);
  
}
