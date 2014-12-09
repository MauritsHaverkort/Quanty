
uint32_t ImpurityHFromTB(OperatorType *const Himp,uint32_t NBath,uint32_t NFermionShift, uint32_t NFermion)
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
  int Nk[]={20,20,20};



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

  BeginTiming("TightBindingDOS", &Time);
  TightBindingDOS(TB, Nk, &DOS);
  EndTiming("TightBindingDOS", &Time);

  FilePrintSpectra(&DOS, 1, "TBDOS1", "w");

  BeginTiming("TightBindingHamiltonianToBlockListOfPoles", &Time);
  TightBindingHamiltonianToBlockListOfPoles(TB, Nk, &G, -2.1, 5.5, 100);
  EndTiming("TightBindingHamiltonianToBlockListOfPoles", &Time);

  BeginTiming("BlockListOfPolesToSpectra", &Time);
  BlockListOfPolesToSpectra(G, Gspec, -3.0, 6.5, 100);
  EndTiming("BlockListOfPolesToSpectra", &Time);

  FilePrintSpectra(Gspec, TB.NFermion*TB.NFermion, "TBDOS2", "w");
  

  BeginTiming("BlockListOfPolesToBlockTridiagonalMatrix", &Time);
  BlockListOfPolesToBlockTridiagonalMatrix(G, &Gtri, 100, DBL_EPSILON);
  EndTiming("BlockListOfPolesToBlockTridiagonalMatrix", &Time);
  
  BeginTiming("BlockTriDiagonalMatrixToSpectra", &Time);
  BlockTriDiagonalMatrixToSpectra(&Gtri, GTriSpec, -3.0, 6.5, 100, 0.05);
  EndTiming("BlockTriDiagonalMatrixToSpectra", &Time);

  AT.Ni=5; AT.Nj=8; AT.QComplex=1;
  InitCompactIdentityMatrix(&AT);

  BeginTiming("BlockListOfPolesToBlockTridiagonalMatrixReduced", &Time);
  BlockListOfPolesToBlockTridiagonalMatrixReduceDimension(G, &Gtrireduced, &AT, 100, DBL_EPSILON);
  EndTiming("BlockListOfPolesToBlockTridiagonalMatrixReduced", &Time);
  
  BeginTiming("BlockTriDiagonalMatrixToSpectra", &Time);
  BlockTriDiagonalMatrixToSpectra(&Gtrireduced, GTriSpecreduced, -3.0, 6.5, 100, 0.05);
  EndTiming("BlockTriDiagonalMatrixToSpectra", &Time);
  
//  PrintBlockTriDiagonalMatrix(Gtrireduced);
  
  FilePrintSpectra(GTriSpecreduced, 25, "TBDOS3reduced", "w");
  
  BlockTridiagonalMatrixToImpurityHamiltonian(Himp, &Gtrireduced, NBath+1, 1, NFermionShift, NFermion);
  
  EndTiming("Total", &Time);
  printf("\n\n***********************************************************************\n\n\n");
  PrintTiming(Time); 
  return 0;
}

uint32_t ImpurityModelFromTB(uint32_t Z, uint32_t dCount)
{
  TimeType Time;
  InitializeTiming(&Time);
  BeginTiming("Total", &Time);
  printf("\n========================================================");
  printf("\n====       start impurity model from TB      ===========");
  printf("\n========================================================\n");
  fflush(stdout);

  uint32_t i, j;
  double Bmag[3], Hex[3], zeta2p, zeta3d;
  double red2P3D, red3D3D, Fdd[3], Fpd[2], Gpd[2];
  double Udd, Upd, Delta, ep, ed, el;
  uint32_t NFermion;
  uint16_t Index2Pdn[3] = {0, 1, 2};
  uint16_t Index2Pup[3] = {3, 4, 5};
  uint16_t Index3Ddn[5] = {6, 7, 8, 9, 10};
  uint16_t Index3Dup[5] = {11, 12, 13, 14, 15};
  uint16_t Index3D[10] = {6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
  HarmonicsExpansion D4hAkma1g, D4hAkmb1g, D4hAkmb2g, D4hAkmeg;
  double HoppingEnergiesa1g[4]={1., 0., 0., 0.};
  double HoppingEnergiesb1g[4]={0., 1., 0., 0.};
  double HoppingEnergiesb2g[4]={0., 0., 1., 0.};
  double HoppingEnergieseg[4] ={0., 0., 0., 1.};
  HarmonicsExpansion Vx, Vy, Vz, Vleft, Vright;
  int Vl[]={1, 1, 1};
  int Vm[]={-1, 0, 1};
  double VxRe[]={sqrt(0.5), 0.0, sqrt(0.5)};
  double VyRe[]={sqrt(0.5), 0.0, -sqrt(0.5)};
  double VzRe[]={0.0, 1.0, 0.0};
  double VlRe[]={1.0, 0.0, 0.0};
  double VrRe[]={0.0, 0.0, 1.0};
  OperatorType OppH, OppSx, OppSy, OppSz, OppSsqr;
  OperatorType OppLx, OppLy, OppLz, OppLsqr;
  OperatorType OppJx, OppJy, OppJz, OppJsqr;
  OperatorType Oppldots2p, Oppldots3d;
  OperatorType OppNthreehalf, OppNfivehalf;
  OperatorType OppUdd, OppUpd, OppNd, OppNp;
  OperatorType OppNda1g, OppNdb1g, OppNdb2g, OppNdeg;
  OperatorType Opptemp, Opptemp2;
  OperatorType TXASx, TXASy, TXASz, TXASl, TXASr, TRIXS[3];
  double epsilon, zero;
  double energy, energysqr, Esoc, ECoulomb;
  double Sx, Sy, Sz, Ssqr, Lx, Ly, Lz, Lsqr, Jx, Jy, Jz, Jsqr;
  double Nd, Ndsqr, Nthreehalf, Nfivehalf, Nthreehalfsqr, Nfivehalfsqr;
  double Nda1g, Ndb1g, Ndb2g, Ndeg;
  uint32_t NBitsKey;
  WaveFunctionType psi, Opsi;
  uint16_t *CreateString;
  TriDiagonalMatrixType SpecTri;
  SpectraType Spectra[35];
  MatrixType Density;
  double Hist[30], OccHist[10];
  uint32_t doXAS, doRIXS, doFY;
  uint16_t Nelectrons;
  MatrixType Rotated;
  uint32_t NBathv;

  uint32_t numwRes = 20, numwResFY = 300;
  double gamma, gammaRes, NKrylovXAS, NKrylovRIXS, wRes[numwRes],wResFY[numwResFY],Gfwhm;
  SpectraType RIXSSpectra[3][numwRes];
  

  
  printf("\n========================================================");
  printf("\n==========          set parameters           ===========");
  printf("\n========================================================\n");
  
  doXAS=1; doRIXS=1; doFY=1;
  
  Bmag[0]= 0.000000; Bmag[1]= 0.000000;  Bmag[2]= 0.000000; /* x, y, z (units of eV) */
  Hex[0] = 0.000000; Hex[1] = 0.000000; Hex[2] = 0.000010; /* x, y, z (units of eV) */
  
  if(get2p3dAtomicInfo(Z,dCount,Fdd,Fpd,Gpd,&zeta2p,&zeta3d))
  {
    return 1;
  }
  red2P3D = 0.8; red3D3D = 0.8;
  Fdd[1] *= red3D3D; Fdd[2] *= red3D3D; Fpd[1] *= red2P3D; 
  Gpd[0] *= red2P3D; Gpd[1] *= red2P3D;

  Udd = 7.0;
  Upd = 8.5;
  Delta = 4.5;      /* energy difference between dn and dn+1 ligand hole */
  
  Fdd[0] = Udd + 2.0/63.0 * (Fdd[1] + Fdd[2]);
  Fpd[0] = Upd + 1.0/15.0 * Gpd[0] + 3.0/70.0 * Gpd[1];
  
  if(dCount == 1)
    ed = Delta - 1. * Udd - 6. * Upd;
  else if(dCount == 2)
    ed = Delta - 2. * Udd - 6. * Upd;
  else if(dCount == 3)
    ed = Delta - 3. * Udd - 6. * Upd;
  else if(dCount == 4)
    ed = Delta - 4. * Udd - 6. * Upd;
  else if(dCount == 5)
    ed = Delta - 5. * Udd - 6. * Upd;
  else if(dCount == 7)
    ed = Delta - 6. * Udd - 6. * Upd;
  else if(dCount == 7)
    ed = Delta - 7. * Udd - 6. * Upd;
  else if(dCount == 8)
    ed = Delta - 8. * Udd - 6. * Upd;
  else if(dCount == 9)
    ed = Delta - 9. * Udd - 6. * Upd;
  else if(dCount == 10)
    ed = Delta - 10. * Udd - 6. * Upd;
  else
    ed = 0;
  
  ep = 0;  el = 0.0;


  ed += 4.0; //accounting for shift of HAMR
  
  
  ed -= 3.0; //temporary for testing...Ni
//  ed -= 10.0; //temporary for testing...Co
  
  
  
  /* Choose the number of bath orbitals, and set
   * the number of electrons and fermions accordingly */
  NBathv = 15;
  Nelectrons = 6 + dCount + 10*NBathv;
  NFermion = 6 + 10 + 10*(NBathv);
  
  epsilon = 1.0e-3; zero = 1.0e-3;

  NBitsKey=20;
  
  /* Create a starting wavefunction. High spin state. */
  CreateString=calloc(6+10+10*NBathv, sizeof(uint16_t));
  
  /*2p electrons*/
  for (i=0; i<6; i++)
    CreateString[i] = i;
  
  /*make a highspin 3d config and fill the ligand orbitals*/
  j = 0;
  for (i=6; i<Nelectrons; i++)
  {
    if((i+j)%2 == 0 && i < (6+dCount) && j < (10-dCount) )
    {
      j++; CreateString[i] = i+j;
    }
    else if(i < (6+dCount))
    {
      CreateString[i] = i+j;
    }
    else
    {
      CreateString[i] = i + (10-dCount);
    }
  }
 
  /* Used for number operators */
  PotentialExpandedOnC("D4h", 2, HoppingEnergiesa1g, &D4hAkma1g);
  PotentialExpandedOnC("D4h", 2, HoppingEnergiesb1g, &D4hAkmb1g);
  PotentialExpandedOnC("D4h", 2, HoppingEnergiesb2g, &D4hAkmb2g);
  PotentialExpandedOnC("D4h", 2, HoppingEnergieseg,  &D4hAkmeg );

  
  
  /* Polarization dependent transition operators for XAS and RIXS.
     These are created from potentials expanded on spherical harmonics.
     Here the expanded potentials are given.  */  
  Vx.N=3; Vx.QComplex=0; Vx.l=Vl; Vx.m=Vm; Vx.Re=VxRe;
  Vy.N=3; Vy.QComplex=0; Vy.l=Vl; Vy.m=Vm; Vy.Re=VyRe;
  Vz.N=3; Vz.QComplex=0; Vz.l=Vl; Vz.m=Vm; Vz.Re=VzRe;
  Vleft.N=3; Vleft.QComplex=0; Vleft.l=Vl; Vleft.m=Vm; Vleft.Re=VlRe;
  Vright.N=3; Vright.QComplex=0; Vright.l=Vl; Vright.m=Vm; Vright.Re=VrRe;
  
  
  printf("\n========================================================");
  printf("\n==========         define operators          ===========");
  printf("\n========================================================\n");
  BeginTiming("define operators", &Time);
  
  
  ImpurityHFromTB(&OppH,NBathv,6,NFermion);
  OperatorChop(&OppH, epsilon);
  //PrintOperatorReadable(OppH);
  
  /*
  OperatorConjugateTranspose(OppH, &Opptemp);
  OperatorAdd(&Opptemp, OppH, -1.0, 0.0);
  OperatorChop(&Opptemp, DBL_EPSILON*1000);
  printf("\nPrinting Opptemp = H - H* \n");
  PrintOperatorReadable(Opptemp);
  FreeOperator(&Opptemp);
  */
  
  
  CreateOperatorSx(NFermion, 5, Index3Dup, Index3Ddn, &OppSx);
  CreateOperatorSy(NFermion, 5, Index3Dup, Index3Ddn, &OppSy);
  CreateOperatorSz(NFermion, 5, Index3Dup, Index3Ddn, &OppSz);
  CreateOperatorSsqr(NFermion, 5, Index3Dup, Index3Ddn, &OppSsqr);
  
  CreateOperatorLx(NFermion, 2, Index3Dup, Index3Ddn, &OppLx);
  CreateOperatorLy(NFermion, 2, Index3Dup, Index3Ddn, &OppLy);
  CreateOperatorLz(NFermion, 2, Index3Dup, Index3Ddn, &OppLz);
  CreateOperatorLsqr(NFermion, 2, Index3Dup, Index3Ddn, &OppLsqr);
  
  CreateOperatorJx(NFermion, 2, Index3Dup, Index3Ddn, &OppJx);
  CreateOperatorJy(NFermion, 2, Index3Dup, Index3Ddn, &OppJy);
  CreateOperatorJz(NFermion, 2, Index3Dup, Index3Ddn, &OppJz);
  CreateOperatorJsqr(NFermion, 2, Index3Dup, Index3Ddn, &OppJsqr);

  CreateOperatorldots(NFermion, 1, Index2Pup, Index2Pdn, &Oppldots2p);
  CreateOperatorldots(NFermion, 2, Index3Dup, Index3Ddn, &Oppldots3d);
  
  CreateOperatorCoulombSingleShell(NFermion, 2, Index3Dup, Index3Ddn, Fdd, &OppUdd);
  CreateOperatorCoulombTwoShellsConstantOccupation(NFermion, 1, 2, Index2Pup, Index2Pdn, Index3Dup, Index3Ddn, Fpd, Gpd, &OppUpd);
  
  /*
  OperatorConjugateTranspose(OppUdd, &Opptemp);
  OperatorAddNormalOrder(&Opptemp, OppUdd, -1.0, 0.0, &QDetNotOrderedFcFaBcBa);
  OperatorChop(&Opptemp, DBL_EPSILON*1000);
  printf("\nPrinting Opptemp = OppUdd - OppUdd* \n");
  PrintOperatorReadable(Opptemp);
  FreeOperator(&Opptemp);
  
  OperatorConjugateTranspose(OppUpd, &Opptemp);
  OperatorAddNormalOrder(&Opptemp, OppUpd, -1.0, 0.0, &QDetNotOrderedFcFaBcBa);
  OperatorChop(&Opptemp, DBL_EPSILON*1000);
  printf("\nPrinting Opptemp = OppUpd - OppUpd* \n");
  PrintOperatorReadable(Opptemp);
  FreeOperator(&Opptemp);
  */
  //CreateOperatorCF(NFermion, 2, 2, Index3Dup, Index3Ddn, Index3Dup, Index3Ddn, D4hAkmCF, &OppCF);
  
  CreateOperatorNumberShell(NFermion, 10, Index3D, &OppNd);
  ////CreateOperatorNumberShell(NFermion, 10, IndexLD, &OppNL);
  CreateOperatorCF(NFermion, 2, 2, Index3Dup, Index3Ddn, Index3Dup, Index3Ddn, D4hAkma1g, &OppNda1g);
  CreateOperatorCF(NFermion, 2, 2, Index3Dup, Index3Ddn, Index3Dup, Index3Ddn, D4hAkmb1g, &OppNdb1g);
  CreateOperatorCF(NFermion, 2, 2, Index3Dup, Index3Ddn, Index3Dup, Index3Ddn, D4hAkmb2g, &OppNdb2g);
  CreateOperatorCF(NFermion, 2, 2, Index3Dup, Index3Ddn, Index3Dup, Index3Ddn, D4hAkmeg,  &OppNdeg );
  
  
  
  
  OperatorCopy(&OppNthreehalf, Oppldots3d); /* <ls> = -3/2n(j=3/2)+n(j=5/2) */
  OperatorAdd(&OppNthreehalf, OppNd, -1.0, 0.0);
  OperatorScale(&OppNthreehalf, -0.4, 0.0);
  OperatorCopy(&OppNfivehalf, OppNd);
  OperatorAdd(&OppNfivehalf, OppNthreehalf, -1.0, 0.0);
  
  
    /* create the transition operator for XAS */
  CreateOperatorCF(NFermion, 2, 1, Index3Dup, Index3Ddn, Index2Pup, Index2Pdn, Vx, &TXASx);
  CreateOperatorCF(NFermion, 2, 1, Index3Dup, Index3Ddn, Index2Pup, Index2Pdn, Vy, &TXASy);
  CreateOperatorCF(NFermion, 2, 1, Index3Dup, Index3Ddn, Index2Pup, Index2Pdn, Vz, &TXASz);
  CreateOperatorCF(NFermion, 2, 1, Index3Dup, Index3Ddn, Index2Pup, Index2Pdn, Vleft, &TXASl);
  CreateOperatorCF(NFermion, 2, 1, Index3Dup, Index3Ddn, Index2Pup, Index2Pdn, Vright, &TXASr);

  /*Create transition operators for RIXS */
  CreateOperatorCF(NFermion, 1, 2, Index2Pup, Index2Pdn, Index3Dup, Index3Ddn, Vx, &TRIXS[0]);
  CreateOperatorCF(NFermion, 1, 2, Index2Pup, Index2Pdn, Index3Dup, Index3Ddn, Vy, &TRIXS[1]);
  CreateOperatorCF(NFermion, 1, 2, Index2Pup, Index2Pdn, Index3Dup, Index3Ddn, Vz, &TRIXS[2]);
  
  
  /*Create decay operators for RIXS*/
 /* for(i = 0; i < 3; i++)
  {
    OperatorConjugateTranspose(TRIXS[i],&Opptemp);
    FreeOperator(TRIXS[i]);
    OperatorCopy(&TRIXS[i], Opptemp);
    FreeOperator(&Opptemp);
  }*/
//  OperatorConjugateTranspose(TXASx,&TRIXS[0]);
//  OperatorTranspose(TXASy,&TRIXS[1]);
//  OperatorTranspose(TXASz,&TRIXS[2]);

    /* Define rotation matrix */
  
  printf("\n========================================================");
  printf("\n====  Define rotation matrix to 'real' orbitals   ======");
  printf("\n========================================================\n");
  
  Rotated.QComplex = 0;
  Rotated.Ni = NFermion;
  Rotated.Nj = NFermion;
  InitMatrix(&Rotated);
  for (i=0; i<Rotated.Ni; i++)
  {
    for (j=0; j<Rotated.Nj; j++)
    {
      Rotated.VRe[i][j] = 0.0;
    }
  }
//    for (i=0; i<6; i++)
//    {
//      Rotated.VRe[i][i] = 1.0;
//    }
  Rotated.VRe[0][0] =-sqrt(0.5);    Rotated.VRe[0][2] = sqrt(0.5);
  Rotated.VRe[3][3] =-sqrt(0.5);    Rotated.VRe[3][5] = sqrt(0.5);

  Rotated.VRe[1][1] = 1.0;
  Rotated.VRe[4][4] = 1.0;

  Rotated.VRe[2][0] =-sqrt(0.5);    Rotated.VRe[2][2] =-sqrt(0.5);
  Rotated.VRe[5][3] =-sqrt(0.5);    Rotated.VRe[5][5] =-sqrt(0.5);
  
  for (i=0; i<NBathv+1; i++)
  {
    Rotated.VRe[ 6+i*10][ 6+i*10] =-sqrt(0.5);   Rotated.VRe[ 6+i*10][10+i*10] = sqrt(0.5);
    Rotated.VRe[ 11+i*10][ 11+i*10] =-sqrt(0.5);   Rotated.VRe[ 11+i*10][15+i*10] = sqrt(0.5);
    
    Rotated.VRe[ 7+i*10][ 7+i*10] =-sqrt(0.5);   Rotated.VRe[ 7+i*10][9+i*10] = sqrt(0.5);
    Rotated.VRe[ 12+i*10][ 12+i*10] =-sqrt(0.5);   Rotated.VRe[ 12+i*10][14+i*10] = sqrt(0.5);
    
    Rotated.VRe[8+i*10][8+i*10] = 1.0;
    Rotated.VRe[13+i*10][13+i*10] = 1.0;
    
    Rotated.VRe[9+i*10][ 7+i*10] =-sqrt(0.5);   Rotated.VRe[9+i*10][9+i*10] =-sqrt(0.5);
    Rotated.VRe[14+i*10][ 12+i*10] =-sqrt(0.5);   Rotated.VRe[14+i*10][14+i*10] =-sqrt(0.5);
    
    Rotated.VRe[10+i*10][ 6+i*10] = sqrt(0.5);   Rotated.VRe[10+i*10][10+i*10] = sqrt(0.5);
    Rotated.VRe[15+i*10][ 11+i*10] = sqrt(0.5);   Rotated.VRe[15+i*10][15+i*10] = sqrt(0.5);
  }
  
  printf("\n========================================================");
  printf("\n==== Rotate the operators from lm to real orbital ======");
  printf("\n========================================================\n");
  
  /* Rotate Operators */
  RotateOperatorNormalOrder(&Opptemp, OppSx, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&OppSx); OperatorCopy(&OppSx, Opptemp); FreeOperator(&Opptemp);
  RotateOperatorNormalOrder(&Opptemp, OppSy, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&OppSy); OperatorCopy(&OppSy, Opptemp); FreeOperator(&Opptemp);
  RotateOperatorNormalOrder(&Opptemp, OppSz, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&OppSz); OperatorCopy(&OppSz, Opptemp); FreeOperator(&Opptemp);
  RotateOperatorNormalOrder(&Opptemp, OppSsqr, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&OppSsqr); OperatorCopy(&OppSsqr, Opptemp); FreeOperator(&Opptemp);
  
  RotateOperatorNormalOrder(&Opptemp, OppLx, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&OppLx); OperatorCopy(&OppLx, Opptemp); FreeOperator(&Opptemp);
  RotateOperatorNormalOrder(&Opptemp, OppLy, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&OppLy); OperatorCopy(&OppLy, Opptemp); FreeOperator(&Opptemp);
  RotateOperatorNormalOrder(&Opptemp, OppLz, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&OppLz); OperatorCopy(&OppLz, Opptemp); FreeOperator(&Opptemp);
  RotateOperatorNormalOrder(&Opptemp, OppLsqr, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&OppLsqr); OperatorCopy(&OppLsqr, Opptemp); FreeOperator(&Opptemp);
  
  RotateOperatorNormalOrder(&Opptemp, OppJx, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&OppJx); OperatorCopy(&OppJx, Opptemp); FreeOperator(&Opptemp);
  RotateOperatorNormalOrder(&Opptemp, OppJy, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&OppJy); OperatorCopy(&OppJy, Opptemp); FreeOperator(&Opptemp);
  RotateOperatorNormalOrder(&Opptemp, OppJz, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&OppJz); OperatorCopy(&OppJz, Opptemp); FreeOperator(&Opptemp);
  RotateOperatorNormalOrder(&Opptemp, OppJsqr, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&OppJsqr); OperatorCopy(&OppJsqr, Opptemp); FreeOperator(&Opptemp);
  
  RotateOperatorNormalOrder(&Opptemp, Oppldots2p, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&Oppldots2p); OperatorCopy(&Oppldots2p, Opptemp); FreeOperator(&Opptemp);
  RotateOperatorNormalOrder(&Opptemp, Oppldots3d, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&Oppldots3d); OperatorCopy(&Oppldots3d, Opptemp); FreeOperator(&Opptemp);
  
  RotateOperatorNormalOrder(&Opptemp, OppUdd, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&OppUdd); OperatorCopy(&OppUdd, Opptemp); FreeOperator(&Opptemp);
  RotateOperatorNormalOrder(&Opptemp, OppUpd, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&OppUpd); OperatorCopy(&OppUpd, Opptemp); FreeOperator(&Opptemp);

  RotateOperatorNormalOrder(&Opptemp, OppNd, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&OppNd); OperatorCopy(&OppNd, Opptemp); FreeOperator(&Opptemp);
  RotateOperatorNormalOrder(&Opptemp, OppNda1g, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&OppNda1g); OperatorCopy(&OppNda1g, Opptemp); FreeOperator(&Opptemp);
  RotateOperatorNormalOrder(&Opptemp, OppNdb1g, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&OppNdb1g); OperatorCopy(&OppNdb1g, Opptemp); FreeOperator(&Opptemp);
  RotateOperatorNormalOrder(&Opptemp, OppNdb2g, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&OppNdb2g); OperatorCopy(&OppNdb2g, Opptemp); FreeOperator(&Opptemp);
  RotateOperatorNormalOrder(&Opptemp, OppNdeg, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&OppNdeg); OperatorCopy(&OppNdeg, Opptemp); FreeOperator(&Opptemp);
  RotateOperatorNormalOrder(&Opptemp, OppNthreehalf, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&OppNthreehalf); OperatorCopy(&OppNthreehalf, Opptemp); FreeOperator(&Opptemp);
  RotateOperatorNormalOrder(&Opptemp, OppNfivehalf, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&OppNfivehalf); OperatorCopy(&OppNfivehalf, Opptemp); FreeOperator(&Opptemp);
  
  /* Rotate spectral operators */
  
  RotateOperatorNormalOrder(&Opptemp, TXASx, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&TXASx); OperatorCopy(&TXASx, Opptemp); FreeOperator(&Opptemp);
  RotateOperatorNormalOrder(&Opptemp, TXASy, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&TXASy); OperatorCopy(&TXASy, Opptemp); FreeOperator(&Opptemp);
  RotateOperatorNormalOrder(&Opptemp, TXASz, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&TXASz); OperatorCopy(&TXASz, Opptemp); FreeOperator(&Opptemp);
  RotateOperatorNormalOrder(&Opptemp, TXASl, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&TXASl); OperatorCopy(&TXASl, Opptemp); FreeOperator(&Opptemp);
  RotateOperatorNormalOrder(&Opptemp, TXASr, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&TXASr); OperatorCopy(&TXASr, Opptemp); FreeOperator(&Opptemp);
  
  for(i = 0; i < 3; i++)
  {
    RotateOperatorNormalOrder(&Opptemp, TRIXS[i], Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&TRIXS[i]); OperatorCopy(&TRIXS[i], Opptemp); FreeOperator(&Opptemp);
  }
  
  
  OperatorConjugateTranspose(OppUdd, &Opptemp);
  OperatorAdd(&Opptemp, OppUdd, -1.0, 0.0);
  InitOperatorCopySettings(&Opptemp2, Opptemp);
  NormalOrderOperator(Opptemp, &Opptemp2, &QDetNotOrderedFcFaBaBc); // different orderings possible, all need to be the same to get annihilation
  OperatorChop(&Opptemp2, 1000*DBL_EPSILON);
  printf("\nPrinting Opptemp = Udd - Udd^+ \n");
  PrintOperatorReadable(Opptemp2);
  FreeOperator(&Opptemp);
  FreeOperator(&Opptemp2);
  
  /*
  OperatorConjugateTranspose(OppUdd, &Opptemp);
  OperatorAddNormalOrder(&Opptemp, OppUdd, -1.0, 0.0, &QDetNotOrderedFcFaBcBa);
  OperatorChop(&Opptemp, DBL_EPSILON*1000);
  printf("\nPrinting Opptemp = Udd - Udd* \n");
  PrintOperatorReadable(Opptemp);
  FreeOperator(&Opptemp);
  */
  /*                                          */
  /* sum operators to create the Hamiltonian */
  /*                                          */
  OperatorAdd(&OppH, OppUdd, 1.0, 0.0);
  OperatorAdd(&OppH, OppUpd, 1.0, 0.0);
  OperatorAdd(&OppH, Oppldots2p, zeta2p, 0.0);
  OperatorAdd(&OppH, Oppldots3d, zeta3d, 0.0);  
  OperatorAdd(&OppH, OppNd, ed, 0.0);
  OperatorAdd(&OppH, OppSx, 2*Bmag[0]+Hex[0], 0.0);
  OperatorAdd(&OppH, OppSy, 2*Bmag[1]+Hex[1], 0.0);
  OperatorAdd(&OppH, OppSz, 2*Bmag[2]+Hex[2], 0.0);
  OperatorAdd(&OppH, OppLx, Bmag[0], 0.0);
  OperatorAdd(&OppH, OppLy, Bmag[1], 0.0);
  OperatorAdd(&OppH, OppLz, Bmag[2], 0.0);
  
  OperatorChop(&OppH, DBL_EPSILON*1000);

  
  
  OperatorConjugateTranspose(OppH, &Opptemp);
  OperatorAdd(&Opptemp, OppH, -1.0, 0.0);
  InitOperatorCopySettings(&Opptemp2, Opptemp);
  NormalOrderOperator(Opptemp, &Opptemp2, &QDetNotOrderedFcFaBaBc); // different orderings possible, all need to be the same to get annihilation
  OperatorChop(&Opptemp2, 1000*DBL_EPSILON);
  printf("\nPrinting Opptemp = H - H^+ \n");
  PrintOperatorReadable(Opptemp2);
  FreeOperator(&Opptemp);
  FreeOperator(&Opptemp2);
  
  
  EndTiming("define operators", &Time);
  
  
  printf("\n========================================================");
  printf("\n==========         find groundstate          ===========");
  printf("\n========================================================\n");
  BeginTiming("find groundstate", &Time);
  
 
  
  /*
   * Shouldn't need this for manganese, because only spin up electrons
  CreateString[6] = 6;
  CreateString[7] = 8;
  CreateString[8] = 10;
  CreateString[9] = 12;
  CreateString[10] = 14;
  CreateString[11] = 15;
  CreateString[12] = 16;
  CreateString[13] = 17;
  CreateString[14] = 18;
  CreateString[15] = 19;
  */
  
  printf("\n\n\n\n");
  printf("Wavefunction string: ");
  for(i = 0; i < Nelectrons; i++)
  {
    printf("%d ", CreateString[i]);
  }
  printf("\n\n\n\n");

  

  CreateWaveFunctionFromString(&psi, CreateString, Nelectrons, NFermion, 0, NBitsKey);
  //  WaveFunctionRandomize(&psi);

   epsilon = 1.0e-4; zero = 1.0e-4;
  BeginTiming("LanczosGroundState", &Time);
  LanczosGroundState(OppH, &psi, &epsilon, zero, 0, &psi, 0);

  WaveFunctionRandomize(&psi);

  LanczosGroundState(OppH, &psi, &epsilon, zero, 0, &psi, 0);
  
  WaveFunctionRandomize(&psi);
  
  LanczosGroundState(OppH, &psi, &epsilon, zero, 0, &psi, 0);
  LanczosGroundState(OppH, &psi, &epsilon, zero, 0, &psi, 0);
  LanczosGroundState(OppH, &psi, &epsilon, zero, 0, &psi, 0);
  LanczosGroundState(OppH, &psi, &epsilon, zero, 0, &psi, 0);
  LanczosGroundState(OppH, &psi, &epsilon, zero, 0, &psi, 0);
  LanczosGroundState(OppH, &psi, &epsilon, zero, 0, &psi, 0);
  EndTiming("LanczosGroundState", &Time);
  
  CalculatePrefactorHistogram(Hist, psi, 2.0, 30);
  printf("\nPrefactor Histogram of groundstate\n");
  printf("Interval                                 | percentage of states   | alpha^2 * number of states \n");
  for(i=0;i<30;i++)
  {
    printf("2^-%2i (%10.3E) to 2^-%2i (%10.3E) | %22.15E | %22.15E \n",
           i, pow(0.5,(double)i), i+1, pow(0.5,(double)(i+1)),
           Hist[i], Hist[i]*pow(0.5,(double)(2*i)));
  }
  
  CalculateOccupationHistogram(OccHist, psi, Index3D, 10);
  
  PsiHermitianOperatorPsi(OppH, psi, &energy);
  PsiHermitianOperatorSqrPsi(OppH, psi, &energysqr);
  PsiHermitianOperatorPsi(Oppldots3d, psi, &Esoc);
  Esoc *= zeta3d;
  PsiHermitianOperatorPsi(OppUdd, psi, &ECoulomb);
//  PsiHermitianOperatorPsi(OppHop, psi, &EKin);
  PsiHermitianOperatorPsi(OppSx, psi, &Sx);
  PsiHermitianOperatorPsi(OppSy, psi, &Sy);
  PsiHermitianOperatorPsi(OppSz, psi, &Sz);
  PsiHermitianOperatorPsi(OppSsqr, psi, &Ssqr);
  PsiHermitianOperatorPsi(OppLx, psi, &Lx);
  PsiHermitianOperatorPsi(OppLy, psi, &Ly);
  PsiHermitianOperatorPsi(OppLz, psi, &Lz);
  PsiHermitianOperatorPsi(OppLsqr, psi, &Lsqr);
  PsiHermitianOperatorPsi(OppJx, psi, &Jx);
  PsiHermitianOperatorPsi(OppJy, psi, &Jy);
  PsiHermitianOperatorPsi(OppJz, psi, &Jz);
  PsiHermitianOperatorPsi(OppJsqr, psi, &Jsqr);
  PsiHermitianOperatorPsi(OppNd, psi, &Nd);
  PsiHermitianOperatorPsi(OppNda1g, psi, &Nda1g);
  PsiHermitianOperatorPsi(OppNdb1g, psi, &Ndb1g);
  PsiHermitianOperatorPsi(OppNdb2g, psi, &Ndb2g);
  PsiHermitianOperatorPsi(OppNdeg,  psi, &Ndeg );
  PsiHermitianOperatorPsi(OppNthreehalf, psi, &Nthreehalf);
  PsiHermitianOperatorPsi(OppNfivehalf, psi, &Nfivehalf);
  PsiHermitianOperatorSqrPsi(OppNd, psi, &Ndsqr);
  PsiHermitianOperatorSqrPsi(OppNthreehalf, psi, &Nthreehalfsqr);
  PsiHermitianOperatorSqrPsi(OppNfivehalf, psi, &Nfivehalfsqr);
  EndTiming("find groundstate", &Time);
  PrintTiming(Time);
  
  printf("\n========================================================");
  printf("\n==========   Groundstate expectation values   ===========");
  printf("\n========================================================\n");
  printf("Energy              %22.15E  (variance %22.15E)\n", energy, energysqr-energy*energy);
//  printf("Kinetic energy      %22.15E\n", EKin);
  printf("Relativistic energy %22.15E\n", Esoc);
  printf("Coulomb energy (3d) %22.15E\n", ECoulomb);
  printf("Sx %22.15E Sy %22.15E Sz %22.15E Ssqr %22.15E\n", Sx, Sy, Sz, Ssqr);
  printf("Lx %22.15E Ly %22.15E Lz %22.15E Lsqr %22.15E\n", Lx, Ly, Lz, Lsqr);
  printf("Jx %22.15E Jy %22.15E Jz %22.15E Jsqr %22.15E\n", Jx, Jy, Jz, Jsqr);
  printf("Nd   %22.15E <Nd^2>   %22.15E\n", Nd, Ndsqr);
  printf("Nda1g %22.14E\n", Nda1g);
  printf("Ndb1g %22.14E\n", Ndb1g);
  printf("Ndb2g %22.14E\n", Ndb2g);
  printf("Ndeg  %22.14E\n", Ndeg);
  printf("N3/2 %22.15E <N3/2^2> %22.15E\n", Nthreehalf, Nthreehalfsqr);
  printf("N5/2 %22.15E <N5/2^2> %22.15E\n", Nfivehalf, Nfivehalfsqr);
  printf("\n");
  
  printf("Occupation of d-shell\n");
  for(i=0;i<=10;i++)
  {
    printf("d^%2i | %22.15E\n", i, OccHist[i]);
  }
  printf("\n");
  
  printf("diagonal part of density matrix\n");
  CalculateDensityMatrix(&Density, psi);
  for(i=0;i<Density.Ni;i++)
  {
    printf("%i %22.15E\n",i,Density.VRe[i][i]);
  }
  printf("\n");
  //PrintMatrix(Density);
  epsilon = 1.0e-3; zero = 1.0e-3;
  if(doXAS)
  {
    printf("\n========================================================");
    printf("\n==========                XAS                ===========");
    printf("\n========================================================\n");
    BeginTiming("XAS", &Time);

    SpecTri.N = 150;
    Gfwhm = 0.2;
    InitTriDiagonalMatrix(&SpecTri);
     
    for(i=26;i<31;i++)
    {
      Spectra[i].N=7000;
      Spectra[i].Emin=-75.0;
      Spectra[i].Emax=-30.0;
      Spectra[i].dE = (Spectra[i].Emax - Spectra[i].Emin) / Spectra[i].N;
      Spectra[i].Gamma = 0.4;
      SpectrumAllocate(Spectra[i], &Spectra[i]);
    }
    
    BeginTiming("spectrum", &Time);
    OperatorPsi(TXASx, psi, &Opsi);
    LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[26].Gamma);
    TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[26]);
    //FilePrintPlainTriDiagonalMatrix(SpecTri, "XASTri_x", "w");
    EndTiming("spectrum", &Time);
    
    BeginTiming("spectrum", &Time);
    OperatorPsi(TXASy, psi, &Opsi);
    LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[27].Gamma);
    TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[27]);
    //FilePrintPlainTriDiagonalMatrix(SpecTri, "XASTri_y", "w");
    EndTiming("spectrum", &Time);
    
    BeginTiming("spectrum", &Time);
    OperatorPsi(TXASz, psi, &Opsi);
    LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[28].Gamma);
    TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[28]);
    //FilePrintPlainTriDiagonalMatrix(SpecTri, "XASTri_z", "w");
    EndTiming("spectrum", &Time);
    /*
    BeginTiming("spectrum", &Time);
    OperatorPsi(TXASl, psi, &Opsi);
    LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[29].Gamma);
    TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[29]);
    EndTiming("spectrum", &Time);
    
    BeginTiming("spectrum", &Time);
    OperatorPsi(TXASr, psi, &Opsi);
    LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[30].Gamma);
    TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[30]);
    EndTiming("spectrum", &Time);*/
    BroadenG(&Spectra[26], Gfwhm);
    BroadenG(&Spectra[27], Gfwhm);
    BroadenG(&Spectra[28], Gfwhm);
    FilePrintSpectra(&Spectra[26], 3, "XAS_Spectra","w");
    EndTiming("XAS", &Time);
    PrintTiming(Time);
  }
  if(doFY)
  {
    printf("\n========================================================");
    printf("\n==========       Fluorescence Yield          ===========");
    printf("\n========================================================\n");
    BeginTiming("FY", &Time);

    gamma = 0.1;
    Gfwhm = 0.2;
    gammaRes = 0.4;
    NKrylovXAS = 150;
    numwResFY = 300;
    for(i = 0; i < numwResFY; i++)
      wResFY[i] = -66.0 + ((double)i)*0.1;
  
    epsilon = 1e-05; zero=1e-05;
    /* Just doing xx for now */
    BeginTiming("spectrum", &Time);
    FluorescenceYield(OppH,TXASx,TRIXS,1,psi,energy,epsilon,numwResFY,&Spectra[0],wResFY,gammaRes,Gfwhm,NKrylovXAS);
    EndTiming("spectrum", &Time);

    FilePrintSpectra(&Spectra[0], 1, "Fluorescence_Yield","w");

    EndTiming("FY", &Time);
    PrintTiming(Time);
  }
  if(doRIXS)
  {
    printf("\n========================================================");
    printf("\n==========               RIXS                ===========");
    printf("\n========================================================\n");
    BeginTiming("RIXS", &Time);

    gamma = 0.1;
    Gfwhm = 0.5;
    gammaRes = 0.4;
    NKrylovXAS = 150;
    NKrylovRIXS = 150;
    
    for(i = 0; i < numwRes; i++)
      wRes[i] = -60.0 + ((double)i)*0.6;
  
    epsilon = 1e-04; zero=1e-04;
    for(i=0;i<3;i++)
    {
      for(j=0; j < numwRes; j++)
      {
        RIXSSpectra[i][j].N=4601;
        RIXSSpectra[i][j].Emin=-5.0;
        RIXSSpectra[i][j].Emax=18.0;
        RIXSSpectra[i][j].dE = (RIXSSpectra[i][j].Emax - RIXSSpectra[i][j].Emin) / RIXSSpectra[i][j].N;
        RIXSSpectra[i][j].Gamma = gamma;
        SpectrumAllocate(RIXSSpectra[i][j], &RIXSSpectra[i][j]);
      }
    }

    /* Just set to do xx for testing.  To do xx, xy, xz, set the parameter to 3 in the call, and
       uncomment the file printing lines*/
    BeginTiming("spectrum", &Time);
    ResonantSpectraOrig(OppH,TXASx,TRIXS,1,psi,energy,epsilon,numwRes,RIXSSpectra,wRes,gammaRes,gamma,Gfwhm,NKrylovXAS,NKrylovRIXS);
    EndTiming("spectrum", &Time);

    FilePrintRIXSSpectrum(numwRes,RIXSSpectra,wRes, 1,0, "RIXS_xx_Real","w");
    FilePrintRIXSSpectrum(numwRes,RIXSSpectra,wRes, 0,0, "RIXS_xx_Imag","w");
//    FilePrintRIXSSpectrum(numwRes,RIXSSpectra,wRes, 1,1, "RIXS_xy_Real","w");
//    FilePrintRIXSSpectrum(numwRes,RIXSSpectra,wRes, 0,1, "RIXS_xy_Imag","w");
//    FilePrintRIXSSpectrum(numwRes,RIXSSpectra,wRes, 1,2, "RIXS_xz_Real","w");
//    FilePrintRIXSSpectrum(numwRes,RIXSSpectra,wRes, 0,2, "RIXS_xz_Imag","w");

    EndTiming("RIXS", &Time);
    PrintTiming(Time);
  }
  
  EndTiming("Total", &Time);
  PrintTiming(Time);
  return 0;
  
  
}
