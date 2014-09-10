
void MnOLigandField()
{
  TimeType Time;
  InitializeTiming(&Time);
  BeginTiming("Total", &Time);
  printf("\n========================================================");
  printf("\n==========     start ligand field theory     ===========");
  printf("\n========================================================\n");
  uint32_t i, j;
  double Bmag[3], Hex[3], zeta2p, zeta3d;
  double Rydberg, red2P3D, red3D3D, Fdd[3], Fpd[2], Gpd[2];
  double Udd, Upd, Delta, ep, ed, el;
  uint32_t NFermion = 26;
  uint16_t Index2Pdn[3] = {0, 2, 4};
  uint16_t Index2Pup[3] = {1, 3, 5};
  uint16_t Index3Ddn[5] = {6, 8, 10, 12, 14};
  uint16_t Index3Dup[5] = {7, 9, 11, 13, 15};
  uint16_t IndexLDdn[5] = {16, 18, 20, 22, 24};
  uint16_t IndexLDup[5] = {17, 19, 21, 23, 25};
  uint16_t Index3D[10] = {6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
  uint16_t IndexLD[10] = {16, 17, 18, 19, 20, 21, 22, 23, 24, 25};
  uint16_t jlist[10];
  double jRe[10], jIm[10];
  HarmonicsExpansion D4hAkmHop, D4hAkmCF;
  HarmonicsExpansion D4hAkma1g, D4hAkmb1g, D4hAkmb2g, D4hAkmeg;
  double HoppingEnergiesa1g[4]={1., 0., 0., 0.};
  double HoppingEnergiesb1g[4]={0., 1., 0., 0.};
  double HoppingEnergiesb2g[4]={0., 0., 1., 0.};
  double HoppingEnergieseg[4] ={0., 0., 0., 1.};
  double HoppingEnergies[4], CFEnergies[4];
  double Ea1g, Eb1g, Eb2g, Eeg, Hopa1g, Hopb1g, Hopb2g, Hopeg, pds, pdp;
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
  OperatorType OppUdd, OppUpd, OppNd, OppNL;
  OperatorType OppNda1g, OppNdb1g, OppNdb2g, OppNdeg;
  OperatorType Opptemp, OppHop, OppCF;
  OperatorType TPES[10], TIPES[10], TCPES[6], TXASx, TXASy, TXASz, TXASl, TXASr, TCharge, TRIXS[3];
  double epsilon, zero;
  double energy, energysqr, Esoc, ECoulomb, EKin;
  double Sx, Sy, Sz, Ssqr, Lx, Ly, Lz, Lsqr, Jx, Jy, Jz, Jsqr;
  double Nd, Ndsqr, Nthreehalf, Nfivehalf, Nthreehalfsqr, Nfivehalfsqr;
  double Nda1g, Ndb1g, Ndb2g, Ndeg;
  uint32_t NBitsKey;
  WaveFunctionType psi, Opsi;
  uint16_t CreateString[21];
  TriDiagonalMatrixType SpecTri;
  SpectraType Spectra[35];
  MatrixType Density;
  double Hist[30], OccHist[10];
  uint16_t oadd = 0;
  uint32_t doPES, doIPES, doCPES, doXAS, doXiSpin, doXiCharge, doRIXS;
  uint32_t rotate=0;
  MatrixType Rotated;


  uint32_t numwRes = 20;
  double gamma, gammaRes, NKrylovXAS, NKrylovRIXS, wRes[numwRes];
  SpectraType RIXSSpectra[3][numwRes];


  
  printf("\n========================================================");
  printf("\n==========          set parameters           ===========");
  printf("\n========================================================\n");
  
  doPES=0;
  doIPES=0;
  doCPES=0;
  doXAS=1;
  doXiSpin=0;
  doXiCharge=0;
  doRIXS=1;
  
  Rydberg = 13.6057;
  
  Bmag[0]= 0.000000; /* magnetic field in the x direction (units of eV) */
  Bmag[1]= 0.000000; /* magnetic field in the y direction (units of eV) */
  Bmag[2]= 0.000000; /* magnetic field in the z direction (units of eV) */
  Hex[0] = 0.000000; /* exchange field in the x direction (units of eV) */
  Hex[1] = 0.000000; /* exchange field in the y direction (units of eV) */
  Hex[2] = 0.000000; /* exchange field in the z direction (units of eV) */
  
  zeta2p = 6.8462;
  zeta3d = 0.040;
  
  red2P3D = 0.8;
  red3D3D = 0.8;
  
  Udd = 7.0;
  Upd = 8.5;
  Delta = 5.5;      /* energy differnece between d8 and d9 ligand hole */
  
  Fdd[1] = 10.316 * red3D3D;    /* F2 */
  Fdd[2] = 6.414 * red3D3D;     /* F4 */
  
  Fpd[1] = 6.321 * red2P3D;     /* F2 */
  
  Gpd[0] = 4.606 * red2P3D;     /* G1 */
  Gpd[1] = 2.618 * red2P3D;      /* G3 */
  
  Fdd[0] = Udd + 2.0/63.0 * (Fdd[1] + Fdd[2]);
  Fpd[0] = Upd + 1.0/15.0 * Gpd[0] + 3.0/70.0 * Gpd[1];
  
/* First, solve

    6ep +  5ed + 10eL + 10Udd + 30Upd = 0
    6ep +  6ed +  9eL + 15Udd + 36Upd = Delta
    5ep +  6ed + 10eL + 15Udd + 30Upd = 0
    5ep +  7ed +  9eL + 21Udd + 35Upd = Delta + Udd - Upd

   for ep, ed, eL, to get:

    ep = (1/21)( 10Delta + 15Udd - 90Upd)
    ed = (1/21)( 10Delta - 90Udd - 90Upd)
    eL = (1/21)(-11Delta + 15Udd + 36Upd)

  Now shift eL to zero, to get:

    ed = (1/21)(  21Delta - 105Udd - 126Upd)
       =  Delta - 5Udd - 6Upd
    eL = 0.0
*/
 
  ed = Delta - 5. * Udd - 6. * Upd;
  el = 0.0;
  
  epsilon = 1.0e-3;
  zero = 1.0e-3;
  SpecTri.N = 35;
  InitTriDiagonalMatrix(&SpecTri);
  
  Ea1g = 0.3;
  Eb1g = 0.3;
  Eb2g = -0.2;
  Eeg  = -0.2;
  
  CFEnergies[0] = Ea1g;
  CFEnergies[1] = Eb1g;
  CFEnergies[2] = Eb2g;
  CFEnergies[3] = Eeg;
  
  PotentialExpandedOnY("D4h", 2, CFEnergies, &D4hAkmCF);
  
  /*
  pds = -1.88;
  pdp = pds/2.17;
  Hopa1g = sqrt(3.0) * pds;
  Hopb1g = sqrt(3.0) * pds;
  Hopb2g = 2.0 * pdp;
  Hopeg  = 2.0 * pdp;
  */
  Hopa1g = -1.92;
  Hopb1g = -1.92;
  Hopb2g = 1.15;
  Hopeg  = 1.15;
  
  
  HoppingEnergies[0] = Hopa1g;
  HoppingEnergies[1] = Hopb1g;
  HoppingEnergies[2] = Hopb2g;
  HoppingEnergies[3] = Hopeg;
  
  PotentialExpandedOnY("D4h", 2, HoppingEnergies, &D4hAkmHop);
  PotentialExpandedOnY("D4h", 2, HoppingEnergiesa1g, &D4hAkma1g);
  PotentialExpandedOnY("D4h", 2, HoppingEnergiesb1g, &D4hAkmb1g);
  PotentialExpandedOnY("D4h", 2, HoppingEnergiesb2g, &D4hAkmb2g);
  PotentialExpandedOnY("D4h", 2, HoppingEnergieseg,  &D4hAkmeg );

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
  
  /* first define some operators */
  
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
  
  CreateOperatorCF(NFermion, 2, 2, Index3Dup, Index3Ddn, Index3Dup, Index3Ddn, D4hAkmCF, &OppCF);
  
  CreateOperatorCF(NFermion, 2, 2, IndexLDup, IndexLDdn, Index3Dup, Index3Ddn, D4hAkmHop, &OppHop);
  CreateOperatorCF(NFermion, 2, 2, Index3Dup, Index3Ddn, IndexLDup, IndexLDdn, D4hAkmHop, &Opptemp);
  OperatorAdd(&OppHop, Opptemp, 1.0, 0.0); FreeOperator(&Opptemp);
  
  CreateOperatorNumberShell(NFermion, 10, Index3D, &OppNd);
  CreateOperatorNumberShell(NFermion, 10, IndexLD, &OppNL);
  CreateOperatorCF(NFermion, 2, 2, Index3Dup, Index3Ddn, Index3Dup, Index3Ddn, D4hAkma1g, &OppNda1g);
  CreateOperatorCF(NFermion, 2, 2, Index3Dup, Index3Ddn, Index3Dup, Index3Ddn, D4hAkmb1g, &OppNdb1g);
  CreateOperatorCF(NFermion, 2, 2, Index3Dup, Index3Ddn, Index3Dup, Index3Ddn, D4hAkmb2g, &OppNdb2g);
  CreateOperatorCF(NFermion, 2, 2, Index3Dup, Index3Ddn, Index3Dup, Index3Ddn, D4hAkmeg,  &OppNdeg );
  OperatorCopy(&OppNthreehalf, Oppldots3d); /* <ls> = -3/2n(j=3/2)+n(j=5/2) */
  OperatorAdd(&OppNthreehalf, OppNd, -1.0, 0.0);
  OperatorScale(&OppNthreehalf, -0.4, 0.0);
  OperatorCopy(&OppNfivehalf, OppNd);
  OperatorAdd(&OppNfivehalf, OppNthreehalf, -1.0, 0.0);
  
  /* sum these operators to create the Hamiltonian */
  
  OperatorCopy(&OppH, OppUdd);
  OperatorAdd(&OppH, OppUpd, 1.0, 0.0);
  OperatorAdd(&OppH, Oppldots2p, zeta2p, 0.0);
  OperatorAdd(&OppH, Oppldots3d, zeta3d, 0.0);
  OperatorAdd(&OppH, OppCF, 1.0, 0.0);
  OperatorAdd(&OppH, OppHop, 1.0, 0.0);
  OperatorAdd(&OppH, OppNd, ed, 0.0);
  OperatorAdd(&OppH, OppNL, el, 0.0);
  OperatorAdd(&OppH, OppSx, 2*Bmag[0]+Hex[0], 0.0);
  OperatorAdd(&OppH, OppSy, 2*Bmag[1]+Hex[1], 0.0);
  OperatorAdd(&OppH, OppSz, 2*Bmag[2]+Hex[2], 0.0);
  OperatorAdd(&OppH, OppLx, Bmag[0], 0.0);
  OperatorAdd(&OppH, OppLy, Bmag[1], 0.0);
  OperatorAdd(&OppH, OppLz, Bmag[2], 0.0);
  
  OperatorChop(&OppH, DBL_EPSILON*1000);
  
  PrintOperatorReadable(OppH);
  
  /* create the transition operator for XAS */
  
  CreateOperatorCF(NFermion, 2, 1, Index3Dup, Index3Ddn, Index2Pup, Index2Pdn, Vx, &TXASx);
  CreateOperatorCF(NFermion, 2, 1, Index3Dup, Index3Ddn, Index2Pup, Index2Pdn, Vy, &TXASy);
  CreateOperatorCF(NFermion, 2, 1, Index3Dup, Index3Ddn, Index2Pup, Index2Pdn, Vz, &TXASz);
  CreateOperatorCF(NFermion, 2, 1, Index3Dup, Index3Ddn, Index2Pup, Index2Pdn, Vleft, &TXASl);
  CreateOperatorCF(NFermion, 2, 1, Index3Dup, Index3Ddn, Index2Pup, Index2Pdn, Vright, &TXASr);

  /*Create decay operators for RIXS*/
  for(i = 0; i < 3; i++)
  {
    InitOperator(&TRIXS[i]);
  }
  OperatorConjugateTranspose(TXASx,&TRIXS[0]);
  OperatorConjugateTranspose(TXASy,&TRIXS[1]);
  OperatorConjugateTranspose(TXASz,&TRIXS[2]);


  EndTiming("define operators", &Time);
  
  if (rotate)
  {
    printf("\n========================================================");
    printf("\n==== Rotate the operators from lm to real orbital ======");
    printf("\n========================================================\n");
    
    Rotated.QComplex = 0; //1;
    Rotated.Ni = 26;
    Rotated.Nj = 26;
    InitMatrix(&Rotated);
    for (i=0; i<Rotated.Ni; i++)
    {
      for (j=0; j<Rotated.Nj; j++)
      {
        Rotated.VRe[i][j] = 0.0;
//        Rotated.VIm[i][j] = 0.0;
      }
    }



    for (i=0; i<6; i++)
    {
      Rotated.VRe[i][i] = 1.0;
    }
/*
    Rotated.VIm[0][0] =-sqrt(0.5);    Rotated.VRe[0][4] = sqrt(0.5);
    Rotated.VIm[1][1] =-sqrt(0.5);    Rotated.VRe[1][5] = sqrt(0.5);
    
    Rotated.VRe[2][2] = 1.0;
    Rotated.VRe[3][3] = 1.0;
    
    Rotated.VIm[4][0] =-sqrt(0.5);    Rotated.VRe[4][4] =-sqrt(0.5);
    Rotated.VIm[5][1] =-sqrt(0.5);    Rotated.VRe[5][5] =-sqrt(0.5);
  */  
    for (i=0; i<2; i++)
    {
	/*
      Rotated.VIm[ 6+i*10][ 6+i*10] =-sqrt(0.5);   Rotated.VRe[ 6+i*10][14+i*10] = sqrt(0.5);
      Rotated.VIm[ 7+i*10][ 7+i*10] =-sqrt(0.5);   Rotated.VRe[ 7+i*10][15+i*10] = sqrt(0.5);
      
      Rotated.VIm[ 8+i*10][ 8+i*10] =-sqrt(0.5);   Rotated.VRe[ 8+i*10][12+i*10] = sqrt(0.5);
      Rotated.VIm[ 9+i*10][ 9+i*10] =-sqrt(0.5);   Rotated.VRe[ 9+i*10][13+i*10] = sqrt(0.5);
      
      Rotated.VRe[10+i*10][10+i*10] = 1.0;
      Rotated.VRe[11+i*10][11+i*10] = 1.0;
      
      Rotated.VIm[12+i*10][ 8+i*10] =-sqrt(0.5);   Rotated.VRe[12+i*10][12+i*10] =-sqrt(0.5);
      Rotated.VIm[13+i*10][ 9+i*10] =-sqrt(0.5);   Rotated.VRe[13+i*10][13+i*10] =-sqrt(0.5);
      
      Rotated.VIm[14+i*10][ 6+i*10] = sqrt(0.5);   Rotated.VRe[14+i*10][14+i*10] = sqrt(0.5);
      Rotated.VIm[15+i*10][ 7+i*10] = sqrt(0.5);   Rotated.VRe[15+i*10][15+i*10] = sqrt(0.5);
	  */
	Rotated.VRe[ 6+i*10][ 6+i*10] =-sqrt(0.5);   Rotated.VRe[ 6+i*10][14+i*10] = sqrt(0.5);
    Rotated.VRe[ 7+i*10][ 7+i*10] =-sqrt(0.5);   Rotated.VRe[ 7+i*10][15+i*10] = sqrt(0.5);
    
    Rotated.VRe[ 8+i*10][ 8+i*10] =-sqrt(0.5);   Rotated.VRe[ 8+i*10][12+i*10] = sqrt(0.5);
    Rotated.VRe[ 9+i*10][ 9+i*10] =-sqrt(0.5);   Rotated.VRe[ 9+i*10][13+i*10] = sqrt(0.5);
    
    Rotated.VRe[10+i*10][10+i*10] = 1.0;
    Rotated.VRe[11+i*10][11+i*10] = 1.0;
    
    Rotated.VRe[12+i*10][ 8+i*10] =-sqrt(0.5);   Rotated.VRe[12+i*10][12+i*10] =-sqrt(0.5);
    Rotated.VRe[13+i*10][ 9+i*10] =-sqrt(0.5);   Rotated.VRe[13+i*10][13+i*10] =-sqrt(0.5);
    
    Rotated.VRe[14+i*10][ 6+i*10] = sqrt(0.5);   Rotated.VRe[14+i*10][14+i*10] = sqrt(0.5);
    Rotated.VRe[15+i*10][ 7+i*10] = sqrt(0.5);   Rotated.VRe[15+i*10][15+i*10] = sqrt(0.5);
    }
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
    RotateOperatorNormalOrder(&Opptemp, OppCF, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&OppCF); OperatorCopy(&OppCF, Opptemp); FreeOperator(&Opptemp);
    RotateOperatorNormalOrder(&Opptemp, OppHop, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&OppHop); OperatorCopy(&OppHop, Opptemp); FreeOperator(&Opptemp);
    
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
    RotateOperatorNormalOrder(&Opptemp, OppNL, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&OppNL); OperatorCopy(&OppNL, Opptemp); FreeOperator(&Opptemp);
    RotateOperatorNormalOrder(&Opptemp, OppH, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&OppH); OperatorCopy(&OppH, Opptemp); FreeOperator(&Opptemp);
    
    
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


    for(j = 0; j < 3; j++)
    {
      RotateOperatorNormalOrder(&Opptemp, TRIXS[j], Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
      FreeOperator(&TRIXS[j]); OperatorCopy(&TRIXS[j], Opptemp); FreeOperator(&Opptemp);
    }
  }
  
  printf("\n========================================================");
  printf("\n==========         find groundstate          ===========");
  printf("\n========================================================\n");
  BeginTiming("find groundstate", &Time);
  
  NBitsKey=16;
  for (i=0; i<7; i++)
  {
    CreateString[i] = i;
  }
  for (i=7; i<14; i++)
  {
    CreateString[i] = i+5;
  }
  for (i=14; i<21; i++)
  {
    CreateString[i] = i+5;
  }
  CreateWaveFunctionFromString(&psi, CreateString, 21, NFermion, 0, NBitsKey);
   
  
  BeginTiming("LanczosGroundState", &Time);
  LanczosGroundState(OppH, &psi, &epsilon, zero, 0, &psi, 0);

  WaveFunctionRandomize(&psi);

  LanczosGroundState(OppH, &psi, &epsilon, zero, 0, &psi, 0);
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
  PsiHermitianOperatorPsi(OppHop, psi, &EKin);
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
  printf("\n==========   Groundstte expectation values   ===========");
  printf("\n========================================================\n");
  printf("Energy              %22.15E  (variance %22.15E)\n", energy, energysqr-energy*energy);
  printf("Kinetic energy      %22.15E\n", EKin);
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
  PrintMatrix(Density);
  
  if(doXAS)
  {
    printf("\n========================================================");
    printf("\n==========                XAS                ===========");
    printf("\n========================================================\n");
    BeginTiming("XAS", &Time);
    
    for(i=26;i<31;i++)
    {
      Spectra[i].N=7000;
      Spectra[i].Emin=-65.0;
      Spectra[i].Emax=-30.0;
      Spectra[i].dE = (Spectra[i].Emax - Spectra[i].Emin) / Spectra[i].N;
      Spectra[i].Gamma = 0.3;
      SpectrumAllocate(Spectra[i], &Spectra[i]);
    }
    
    BeginTiming("spectrum", &Time);
    OperatorPsi(TXASx, psi, &Opsi);
    LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[26].Gamma);
    TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[26]);
//    FilePrintPlainTriDiagonalMatrix(SpecTri, "NiLigandField_XASTri_x", "w");
    EndTiming("spectrum", &Time);
    
    BeginTiming("spectrum", &Time);
    OperatorPsi(TXASy, psi, &Opsi);
    LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[27].Gamma);
    TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[27]);
//    FilePrintPlainTriDiagonalMatrix(SpecTri, "NiLigandField_XASTri_y", "w");
    EndTiming("spectrum", &Time);
    
    BeginTiming("spectrum", &Time);
    OperatorPsi(TXASz, psi, &Opsi);
    LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[28].Gamma);
    TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[28]);
//    FilePrintPlainTriDiagonalMatrix(SpecTri, "NiLigandField_XASTri_z", "w");
    EndTiming("spectrum", &Time);
    
    BeginTiming("spectrum", &Time);
    OperatorPsi(TXASl, psi, &Opsi);
    LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[29].Gamma);
    TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[29]);
    EndTiming("spectrum", &Time);
    
    BeginTiming("spectrum", &Time);
    OperatorPsi(TXASr, psi, &Opsi);
    LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[30].Gamma);
    TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[30]);
    EndTiming("spectrum", &Time);
    
    FilePrintSpectra(&Spectra[26], 5, "MnLigandField_XAS","w");
    EndTiming("XAS", &Time);
    PrintTiming(Time);
  }
  if(doRIXS)
  {
    printf("\n========================================================");
    printf("\n==========               RIXS                ===========");
    printf("\n========================================================\n");
    BeginTiming("RIXS", &Time);
//    for(i = 0; i < 3; i++)
//    {
//      InitOperator(&TRIXS[i]);
//    }
//    OperatorConjugateTranspose(TXASx,&TRIXS[0]);
//    OperatorConjugateTranspose(TXASy,&TRIXS[1]);
//    OperatorConjugateTranspose(TXASz,&TRIXS[2]);

    gamma = 0.1;
    gammaRes = 0.4;
    NKrylovXAS = 75;
    NKrylovRIXS = 75;

    for(i = 0; i < numwRes; i++)
      wRes[i] = -50.0 + ((double)i)*0.1;
  
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

    BeginTiming("spectrum", &Time);
//    ResonantSpectra(OppH,TXASx,TRIXS,3,psi,energy,epsilon,&Spectra[26],-77.4,0.4,0.1,70,50);
    ResonantSpectra(OppH,TXASx,TRIXS,3,psi,energy,epsilon,numwRes,RIXSSpectra,wRes,gammaRes,gamma,NKrylovXAS,NKrylovRIXS);



    EndTiming("spectrum", &Time);
    FilePrintRIXSSpectra(numwRes,RIXSSpectra,wRes, 3, "MnLigandField_RIXS","w");
    EndTiming("RIXS", &Time);
    PrintTiming(Time);
  }
  
  EndTiming("Total", &Time);
  PrintTiming(Time);
}

