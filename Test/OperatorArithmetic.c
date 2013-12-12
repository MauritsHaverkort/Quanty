void TestOppMult()
{
  printf("Test multiplication of operators\n");
  OperatorType A,B,AB,BA;
  WaveFunctionType psi,Apsi,Bpsi,ABpsi,BApsi;
  uint8_t *Determinant;
  uint16_t oadd[10];
  psi.NBitsKey=8;
  psi.QComplex=0;
  psi.NBoson=0;
  psi.NFermion=10;
  if(WaveFunctionInit(&psi)){printf("WaveFunctionInit failed in CreateWaveFunctionDMFT\n");fflush(stdout);return;};
  Determinant=calloc(psi.NWord,sizeof(uint8_t *));
  sprintf(psi.Name,"psi test");
  OperatorCr(Determinant,0,psi.NFermion);
  OperatorCr(Determinant,1,psi.NFermion);
  RealWaveFunctionAddElement(&psi, Determinant, 1.0);
  free(Determinant);
  PrintWaveFunctionReadable(psi);
  
  InitNamedOperator(&A,"OppA");
  InitNamedOperator(&B,"OppB");
  InitNamedOperator(&BA,"OppBA");
  InitNamedOperator(&AB,"OppAB");
  A.NFermion=psi.NFermion;
  A.NBoson=psi.NBoson;
  A.NWord=Min(A.NFermion,1+(A.NFermion-1)/8);
  B.NFermion=psi.NFermion;
  B.NBoson=psi.NBoson;
  B.NWord=Min(A.NFermion,1+(B.NFermion-1)/8);
  AB.NFermion=psi.NFermion;
  AB.NBoson=psi.NBoson;
  AB.NWord=Min(AB.NFermion,1+(AB.NFermion-1)/8);
  BA.NFermion=psi.NFermion;
  BA.NBoson=psi.NBoson;
  BA.NWord=Min(A.NFermion,1+(BA.NFermion-1)/8);
  InitOperatorFixedLength(&A,1,8);
  InitOperatorFixedLength(&A,2,8);
  InitOperatorFixedLength(&B,1,8);
  InitOperatorFixedLength(&B,2,8);
  InitOperatorFixedLength(&AB,1,8);
  InitOperatorFixedLength(&AB,2,8);
  InitOperatorFixedLength(&AB,3,8);
  InitOperatorFixedLength(&AB,4,8);
  InitOperatorFixedLength(&BA,1,8);
  InitOperatorFixedLength(&BA,2,8);
  InitOperatorFixedLength(&BA,3,8);
  InitOperatorFixedLength(&BA,4,8);
  
  
  oadd[0]=CrAnMask+0;
  oadd[1]=CrAnMask+2;
  OperatorAddLader(&A, oadd, 2, 1.0, 0.0);
  
  
  oadd[0]=1;
  oadd[1]=0;
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
  
  
  
  
  
  OperatorType Sx,Sy,Sz,S2,Stest;
  uint16_t SpinUp[6],SpinDn[6];
  SpinUp[0]=0;
  SpinUp[1]=2;
  SpinUp[2]=4;
  SpinDn[0]=1;
  SpinDn[1]=3;
  SpinDn[2]=5;
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
  OperatorMultiplyNormalOrder(&Stest, Sz, Sy,&QDetNotOrderedFcFaBcBa);
  OperatorScale(&Stest, -1.0, 0.0);
  OperatorMultiplyNormalOrder(&Stest, Sy, Sz,&QDetNotOrderedFcFaBcBa);
  strcpy(Stest.Name, "i Sx = Sy.Sz - Sz.Sy");
  PrintOperatorReadable(Stest);
  FreeOperator(&Stest);
  
  InitOperatorCopySettings(&Stest, S2);
  OperatorMultiplyNormalOrder(&Stest, Sx, Sz,&QDetNotOrderedFcFaBcBa);
  OperatorScale(&Stest, -1.0, 0.0);
  OperatorMultiplyNormalOrder(&Stest, Sz, Sx,&QDetNotOrderedFcFaBcBa);
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