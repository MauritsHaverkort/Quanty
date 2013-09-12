void TestS()
{
  OperatorType Sx,Sy,Sz,S2,Stest;
  uint16_t SpinUp[2],SpinDn[2];
  SpinUp[0]=0;
  SpinUp[1]=2;
  SpinDn[0]=1;
  SpinDn[1]=3;
  CreateSpinOperators(2, 1, SpinUp, SpinDn, &Sx, &Sy, &Sz, &S2);
  PrintOperatorReadable(Sx);
  PrintOperatorReadable(Sy);
  PrintOperatorReadable(Sz);
  PrintOperatorReadable(S2);
  
  InitOperatorCopySettings(&Stest, S2);
  OperatorMultiplyNormalOrder(&Stest, Sy, Sx,&QDetNotOrderedFcFaBcBa);
  OperatorScale(&Stest, -1.0, 0.0);
  OperatorMultiplyNormalOrder(&Stest, Sx, Sy,&QDetNotOrderedFcFaBcBa);
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