
uint32_t ImpurityModelFromTB()
{
  StuttgartCTRL CTRLin;
  StuttgartHAMR HAMRin;
  TightBindingDef TB;

  CompactMatrixType AT;
  TimeType Time;
  InitializeTiming(&Time);
  BeginTiming("Total", &Time);
  SpectraType DOS,Gspec[64], GTriSpec[64], Gspec2[64], GTriSpecreduced[25];
  BlockListOfPoles G,GLP2;
  BlockTriDiagonalMatrixType Gtri, Gtrireduced;
  int Nk[]={200,200,1};



  /* Read StuttgartNMTO files*/
  if(FileReadStuttgartCTRL(&CTRLin, "CTRL"))
  {
    printf("Error reading StuttgartNMTO CTRL file.\n");
    return 1;
  }
  if(FileReadStuttgartHAMR(&HAMRin, "HAMR"))
  {
    printf("Error reading StuttgartNMTO HAMR file.\n");
    return 1;
  }
  if(TightBindingDefFromStuttgartCTRLHAMR(CTRLin, HAMRin, &TB))
  {
    printf("Error creating TB definition from StuttgartNMTO CTRL and HAMR.");
    return 1;
  }

  /* Create G0 from tight binding definition */
  DOS.N=100;
  DOS.Emin=-3.0;
  DOS.Emax=6.5;
  DOS.dE=(DOS.Emax-DOS.Emin)/DOS.N;
  SpectrumAllocate(DOS, &DOS);



printf("Point 1\n"); fflush(stdout);
  BeginTiming("TightBindingDOS", &Time);
  TightBindingDOS(TB, Nk, &DOS);
  EndTiming("TightBindingDOS", &Time);
printf("Point 2\n"); fflush(stdout);

  FilePrintSpectra(&DOS, 1, "TBDOS1", "w");

printf("Point 3\n"); fflush(stdout);
  BeginTiming("TightBindingHamiltonianToBlockListOfPoles", &Time);
  TightBindingHamiltonianToBlockListOfPoles(TB, Nk, &G, -2.1, 5.5, 100);
  EndTiming("TightBindingHamiltonianToBlockListOfPoles", &Time);
printf("Point 4\n"); fflush(stdout);

  PrintBlockListOfPoles(G);

printf("Point 5\n"); fflush(stdout);
  BeginTiming("BlockListOfPolesToSpectra", &Time);
  BlockListOfPolesToSpectra(G, Gspec, -3.0, 6.5, 100);
  EndTiming("BlockListOfPolesToSpectra", &Time);
printf("Point 6\n"); fflush(stdout);


  PrintBlockListOfPoles(G);
  FilePrintSpectra(Gspec, TB.NFermion*TB.NFermion, "TBDOS2", "w");
  

printf("Point 7\n"); fflush(stdout);  
  BeginTiming("BlockListOfPolesToBlockTridiagonalMatrix", &Time);
  BlockListOfPolesToBlockTridiagonalMatrix(G, &Gtri, 4, DBL_EPSILON);
  EndTiming("BlockListOfPolesToBlockTridiagonalMatrix", &Time);
printf("Point 8\n"); fflush(stdout);


  PrintBlockTriDiagonalMatrix(Gtri);
  

printf("Point 9\n"); fflush(stdout);  
  BeginTiming("BlockTriDiagonalMatrixToSpectra", &Time);
  BlockTriDiagonalMatrixToSpectra(&Gtri, GTriSpec, -3.0, 6.5, 100, 0.05);
  EndTiming("BlockTriDiagonalMatrixToSpectra", &Time);
printf("Point 10\n"); fflush(stdout);  

  
  FilePrintSpectra(GTriSpec, TB.NFermion*TB.NFermion, "TBDOS3", "w");
  

printf("Point 11\n"); fflush(stdout);  
  BeginTiming("BlockTriDiagonalMatrixToBlockListOfPoles", &Time);
  BlockTriDiagonalMatrixToBlockListOfPoles(&Gtri, &GLP2, -2.1, 5.5, 100);
  EndTiming("BlockTriDiagonalMatrixToBlockListOfPoles", &Time);
printf("Point 12\n"); fflush(stdout);  


printf("Point 13\n"); fflush(stdout);
  BeginTiming("BlockListOfPolesToSpectra", &Time);
  BlockListOfPolesToSpectra(GLP2, Gspec2, -3, 6.5, 100);
  EndTiming("BlockListOfPolesToSpectra", &Time);
printf("Point 14\n"); fflush(stdout);  


  FilePrintSpectra(Gspec2, TB.NFermion*TB.NFermion, "TBDOS4", "w");
  

  AT.Ni=5;
  AT.Nj=8;
  AT.QComplex=0;
  InitCompactIdentityMatrix(&AT);
//  AT.VRe[0]=sqrt(0.5);
//  AT.VRe[1]=sqrt(0.5);
  
printf("Point 15\n"); fflush(stdout);  
  BeginTiming("BlockListOfPolesToBlockTridiagonalMatrixReduced", &Time);
  BlockListOfPolesToBlockTridiagonalMatrixReduceDimension(G, &Gtrireduced, &AT, 500, DBL_EPSILON);
  EndTiming("BlockListOfPolesToBlockTridiagonalMatrixReduced", &Time);
printf("Point 16\n"); fflush(stdout);  
  BeginTiming("BlockTriDiagonalMatrixToSpectra", &Time);
  BlockTriDiagonalMatrixToSpectra(&Gtrireduced, GTriSpecreduced, -3.0, 6.5, 100, 0.05);
  EndTiming("BlockTriDiagonalMatrixToSpectra", &Time);
  
printf("Point 17\n"); fflush(stdout);  
  FilePrintSpectra(GTriSpecreduced, TB.NFermion*TB.NFermion, "TBDOS3reduced", "w");
  
  // turn G0 into impurity natural orbitals
  
  
  
  // no spin PES and IPES
  
  // with spin PES and IPES
  
  EndTiming("Total", &Time);
  printf("\n\n***********************************************************************\n\n\n");
  PrintTiming(Time);
  return 0;
}


