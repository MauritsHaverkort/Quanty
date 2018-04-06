#include "LigandField.h"

void CeLigandField()
{
  TimeType Time;
  InitializeTiming(&Time);
  BeginTiming("Total", &Time);
  printf("\n========================================================");
  printf("\n==========     start ligand field theory     ===========");
  printf("\n========================================================\n");
  uint32_t i;
  double Bmag[3], Hex[3], t, zeta3d, zeta4f;
  double Rydberg, red3D4F, red4F4F, Fff[4], Fdf[3], Gdf[3];
  double Uff, Udf, Delta, epsilonf;
  uint16_t Index3Ddn[5]={0,2,4,6,8};
  uint16_t Index3Dup[5]={1,3,5,7,9};
  uint16_t Index4Fdn[7]={10,12,14,16,18,20,22};
  uint16_t Index4Fup[7]={11,13,15,17,19,21,23};
  uint16_t IndexLFdn[7]={24,26,28,30,32,34,36};
  uint16_t IndexLFup[7]={25,27,29,31,33,35,37};
  uint16_t Index4F[14]={10,11,12,13,14,15,16,17,18,19,20,21,22,23};
  uint16_t jlist[14];
  double jRe[14],jIm[14];
  HarmonicsExpansion SphericalAkm,Vx,Vy,Vz,Vleft,Vright;
  int Vl[]={1,1,1};
  int Vm[]={-1,0,1};
  double VxRe[]={sqrt(0.5),0.0, sqrt(0.5)};
  double VyRe[]={sqrt(0.5),0.0,-sqrt(0.5)};
  double VzRe[]={0.0,1.0,0.0};
  double VlRe[]={1.0,0.0,0.0};
  double VrRe[]={0.0,0.0,1.0};
  OperatorType OppH, OppSx, OppSy, OppSz, OppSsqr, OppLx, OppLy, OppLz, OppLsqr, OppJx, OppJy, OppJz, OppJsqr, Oppldots4f, Oppldots3d, OppUff, OppUdf, OppDelta, OppHyb, OppNf, OppNfivehalf, OppNsevenhalf, Opptemp;
  OperatorType TPES[14], TIPES[14], TCPES[10], TXASx, TXASy, TXASz, TXASl, TXASr, TCharge;
  double epsilon, zero;
  double energy, EKin, Esoc, ECoulomb, Sx, Sy, Sz, Ssqr, Lx, Ly, Lz, Lsqr, Jx, Jy, Jz, Jsqr, Nf, Nfivehalf, Nsevenhalf, energysqr, Nfsqr, Nfivehalfsqr, Nsevenhalfsqr;
  uint32_t NBitsKey;
  WaveFunctionType psi, Opsi;
  uint16_t CreateString[24];
  TriDiagonalMatrixType SpecTri;
  SpectraType Spectra[47];
  MatrixType Density;
  double Hist[40],OccHist[15];
  uint16_t oadd=0;
  DeterminantRestrictionType Restrictions;
  uint32_t doPES, doIPES, doCPES, doXAS, doXiSpin, doXiCharge;
  
  /* allocate space for restrictions */
  Restrictions.N=1;
  Restrictions.Min=malloc(Restrictions.N*sizeof(uint32_t));
  Restrictions.Max=malloc(Restrictions.N*sizeof(uint32_t));
  Restrictions.Determinant=malloc(Restrictions.N*sizeof(uint8_t*));
  for(i=0;i<Restrictions.N;i++)
  {
    Restrictions.Determinant[i]=malloc(5*sizeof(uint8_t)); /* we have 38 fermions in the basis, Determinants are coded per byte i.e. blocks of 8 fermions. 8*5=40 and can store the 38 fermions */
  }
  
  printf("\n========================================================");
  printf("\n==========          set parameters           ===========");
  printf("\n========================================================\n");
  
  doPES=0;
  doIPES=0;
  doCPES=1;
  doXAS=0;
  doXiSpin=0;
  doXiCharge=0;
  
  Rydberg=13.6057; /* calculations are done in eV this alows one to imput in mixed eV and Rydberg */
  
  Bmag[0]=0.000000; /* magnetic field in the x direction (units of eV) */
  Bmag[1]=0.000000; /* magnetic field in the y direction (units of eV) */
  Bmag[2]=0.000000; /* magnetic field in the z direction (units of eV) */
  Hex[0]= 0.000000; /* exchange field in the x direction (units of eV) */
  Hex[1]= 0.000000; /* exchange field in the y direction (units of eV) */
  Hex[2]= 0.000000; /* exchange field in the z direction (units of eV) */
  
  t=0.27;            /* hopping between f shel and ligand shell */
  
  zeta3d=0.54721 * Rydberg;       /* spin orbit of the d core shell */
  zeta4f=0.00789 * Rydberg;       /* spin orbit of the f shell */
  
  red3D4F=0.76;     /* reduction of 3d-4f slater integrals from atomic HF values */
  red4F4F=0.57;     /* reduction of 4f-4f slater integrals from atomic HF values */
  
  Uff=8.0;          /* Uff - shperical part of the Coulmb interactin */
  Udf=10.3;         /* Udf - shperical part of the Coulmb interactin */
  Delta=-1.9;       /* energy differnece between f0 and f1 ligand hole */
  
  Fff[1] = 0.9277208 * Rydberg * red4F4F;  /* F2 */
  Fff[2] = 0.5832808 * Rydberg * red4F4F;  /* F4 */
  Fff[3] = 0.4199624 * Rydberg * red4F4F;  /* F6 */
  
  Fdf[1] = 0.5497914 * Rydberg * red3D4F;  /* F2 */
  Fdf[2] = 0.2484883 * Rydberg * red3D4F;  /* F4 */
  
  Gdf[0] = 0.3724648 * Rydberg * red3D4F;  /* G1 */
  Gdf[1] = 0.2178949 * Rydberg * red3D4F;  /* G3 */
  Gdf[2] = 0.1503945 * Rydberg * red3D4F;  /* G5 */
  
  Fff[0] = Uff + (572.0 * Fff[1] + 390.0 * Fff[2] + 500.0 * Fff[3])/27885.0;  /* F0 and U are not defined in the same manner here the correction */
  Fdf[0] = Udf + (99.0 * Gdf[0] + 44.0 * Gdf[1] + 50.0 * Gdf[2])/2310.0;      /* F0 and U are not defined in the same manner here the correction */
  
  epsilonf = Delta - 10 * Udf;
  
  epsilon = 0.00001;           /* smallest value stored as prefactor in wavefunction <== enlarge to speed up calculation !!!                       */
  zero = 10 * DBL_EPSILON;        /* smallest value to be considered different from zero in convergnece tests <== enlarge to speed up calculation !!! */
  SpecTri.N=50;                  /* number of poles in spectra  <== reduce here to speed up calculation of spectra !!!!                              */
  InitTriDiagonalMatrix(&SpecTri);
  
  PotentialExpandedOnC("C1",3,&t,&SphericalAkm); /* the hybridization is taken to be spherical */
  
  Vx.N=3;             /* for XAS I need polarization dependent transition operators */
  Vx.QComplex=0;      /* these are created from potentials expanded on spherical harmonics */
  Vx.l=Vl;            /* here the expanded potentials are given */
  Vx.m=Vm;
  Vx.Re=VxRe;
  Vy.N=3;
  Vy.QComplex=0;
  Vy.l=Vl;
  Vy.m=Vm;
  Vy.Re=VyRe;
  Vz.N=3;
  Vz.QComplex=0;
  Vz.l=Vl;
  Vz.m=Vm;
  Vz.Re=VzRe;
  Vleft.N=3;
  Vleft.QComplex=0;
  Vleft.l=Vl;
  Vleft.m=Vm;
  Vleft.Re=VlRe;
  Vright.N=3;
  Vright.QComplex=0;
  Vright.l=Vl;
  Vright.m=Vm;
  Vright.Re=VrRe;
  
  Restrictions.Min[0]=0;
  Restrictions.Max[0]=3;
  Restrictions.Determinant[0][0] = 0x00; /* bit 1 to 8 are set to zero (d states) */
  Restrictions.Determinant[0][1] = 0x3F; /* bit 9 and 10 are set to zero, (d states) bit 11 to 16 are set to one (f states) */
  Restrictions.Determinant[0][2] = 0xFF; /* bit 17 to 24 are set to one (f states) */
  Restrictions.Determinant[0][3] = 0x00; /* bit 25 to 32 are set to zero (Lf states) */
  Restrictions.Determinant[0][4] = 0x00; /* bit 32 to 38 are set to zero (Lf states), bit 39 and 40 are set to zero, to supress junk in det !- must be zero -! */
  
  printf("\n========================================================");
  printf("\n==========         define operators          ===========");
  printf("\n========================================================\n");
  BeginTiming("define operators", &Time);
  
  /* first define some operators */
  
  CreateOperatorSx(38, 7, Index4Fup, Index4Fdn, &OppSx);
  CreateOperatorSy(38, 7, Index4Fup, Index4Fdn, &OppSy);
  CreateOperatorSz(38, 7, Index4Fup, Index4Fdn, &OppSz);
  CreateOperatorSsqr(38, 7, Index4Fup, Index4Fdn, &OppSsqr);
  CreateOperatorLx(38, 3, Index4Fup, Index4Fdn, &OppLx);
  CreateOperatorLy(38, 3, Index4Fup, Index4Fdn, &OppLy);
  CreateOperatorLz(38, 3, Index4Fup, Index4Fdn, &OppLz);
  CreateOperatorLsqr(38, 3, Index4Fup, Index4Fdn, &OppLsqr);
  CreateOperatorJx(38, 3, Index4Fup, Index4Fdn, &OppJx);
  CreateOperatorJy(38, 3, Index4Fup, Index4Fdn, &OppJy);
  CreateOperatorJz(38, 3, Index4Fup, Index4Fdn, &OppJz);
  CreateOperatorJsqr(38, 3, Index4Fup, Index4Fdn, &OppJsqr);
  CreateOperatorldots(38, 3, Index4Fup, Index4Fdn, &Oppldots4f);
  CreateOperatorldots(38, 2, Index3Dup, Index3Ddn, &Oppldots3d);
  CreateOperatorCoulombSingleShell(38, 3, Index4Fup, Index4Fdn, Fff, &OppUff);
  CreateOperatorCoulombTwoShellsConstantOccupation(38, 2, 3, Index3Dup, Index3Ddn, Index4Fup, Index4Fdn, Fdf, Gdf, &OppUdf);
  CreateOperatorCF(38, 3, 3, IndexLFup, IndexLFdn, Index4Fup, Index4Fdn, SphericalAkm, &OppHyb);
  CreateOperatorCF(38, 3, 3, Index4Fup, Index4Fdn, IndexLFup, IndexLFdn, SphericalAkm, &Opptemp);
  OperatorAdd(&OppHyb, Opptemp, 1.0, 0.0);
  FreeOperator(&Opptemp);
  CreateOperatorNumberShell(38, 14, Index4F, &OppNf);
  OperatorCopy(&OppNfivehalf, Oppldots4f);
  OperatorAdd(&OppNfivehalf, OppNf, -1.5, 0.0);
  OperatorScale(&OppNfivehalf,-0.285714285714285714,0.0);
  OperatorCopy(&OppNsevenhalf, OppNf);
  OperatorAdd(&OppNsevenhalf, OppNfivehalf, -1.0, 0.0);
  OperatorCopy(&OppDelta, OppNf);
  OperatorScale(&OppDelta, epsilonf, 0);
  
  /* summ these operators to create the Hamiltonian */
  
  OperatorCopy(&OppH, OppUff);
  OperatorAdd(&OppH, OppUdf, 1.0, 0.0);
  OperatorAdd(&OppH, Oppldots3d, zeta3d, 0.0);
  OperatorAdd(&OppH, Oppldots4f, zeta4f, 0.0);
  OperatorAdd(&OppH, OppHyb, 1.0, 0.0);
  OperatorAdd(&OppH, OppDelta, 1.0, 0.0);
  OperatorAdd(&OppH, OppSx, 2*Bmag[0]+Hex[0], 0.0);
  OperatorAdd(&OppH, OppSy, 2*Bmag[1]+Hex[1], 0.0);
  OperatorAdd(&OppH, OppSz, 2*Bmag[2]+Hex[2], 0.0);
  OperatorAdd(&OppH, OppLx, Bmag[0], 0.0);
  OperatorAdd(&OppH, OppLy, Bmag[1], 0.0);
  OperatorAdd(&OppH, OppLz, Bmag[2], 0.0);
  
  OperatorChop(&OppH, DBL_EPSILON*1000);
  
  /* create the transition operator for pes */
  
  jIm[0]=0;   jIm[1]=0;
  jlist[0]=11; jlist[1]=12; jRe[0]=-sqrt(6.0/7.0); jRe[1]= sqrt(1.0/7.0);
  CreateOperatorAnnihilateList(38, jlist, jRe, jIm, 2, &TPES[0]);
  jlist[0]=13; jlist[1]=14; jRe[0]=-sqrt(5.0/7.0); jRe[1]= sqrt(2.0/7.0);
  CreateOperatorAnnihilateList(38, jlist, jRe, jIm, 2, &TPES[1]);
  jlist[0]=15; jlist[1]=16; jRe[0]=-sqrt(4.0/7.0); jRe[1]= sqrt(3.0/7.0);
  CreateOperatorAnnihilateList(38, jlist, jRe, jIm, 2, &TPES[2]);
  jlist[0]=17; jlist[1]=18; jRe[0]=-sqrt(3.0/7.0); jRe[1]= sqrt(4.0/7.0);
  CreateOperatorAnnihilateList(38, jlist, jRe, jIm, 2, &TPES[3]);
  jlist[0]=19; jlist[1]=20; jRe[0]=-sqrt(2.0/7.0); jRe[1]= sqrt(5.0/7.0);
  CreateOperatorAnnihilateList(38, jlist, jRe, jIm, 2, &TPES[4]);
  jlist[0]=21; jlist[1]=22; jRe[0]=-sqrt(1.0/7.0); jRe[1]= sqrt(6.0/7.0);
  CreateOperatorAnnihilateList(38, jlist, jRe, jIm, 2, &TPES[5]);
  jlist[0]=10;              jRe[0]=1.0;
  CreateOperatorAnnihilateList(38, jlist, jRe, jIm, 1, &TPES[6]);
  jlist[0]=11; jlist[1]=12; jRe[0]= sqrt(1.0/7.0); jRe[1]= sqrt(6.0/7.0);
  CreateOperatorAnnihilateList(38, jlist, jRe, jIm, 2, &TPES[7]);
  jlist[0]=13; jlist[1]=14; jRe[0]= sqrt(2.0/7.0); jRe[1]= sqrt(5.0/7.0);
  CreateOperatorAnnihilateList(38, jlist, jRe, jIm, 2, &TPES[8]);
  jlist[0]=15; jlist[1]=16; jRe[0]= sqrt(3.0/7.0); jRe[1]= sqrt(4.0/7.0);
  CreateOperatorAnnihilateList(38, jlist, jRe, jIm, 2, &TPES[9]);
  jlist[0]=17; jlist[1]=18; jRe[0]= sqrt(4.0/7.0); jRe[1]= sqrt(3.0/7.0);
  CreateOperatorAnnihilateList(38, jlist, jRe, jIm, 2, &TPES[10]);
  jlist[0]=19; jlist[1]=20; jRe[0]= sqrt(5.0/7.0); jRe[1]= sqrt(2.0/7.0);
  CreateOperatorAnnihilateList(38, jlist, jRe, jIm, 2, &TPES[11]);
  jlist[0]=21; jlist[1]=22; jRe[0]= sqrt(6.0/7.0); jRe[1]= sqrt(1.0/7.0);
  CreateOperatorAnnihilateList(38, jlist, jRe, jIm, 2, &TPES[12]);
  jlist[0]=23;              jRe[0]=1.0;
  CreateOperatorAnnihilateList(38, jlist, jRe, jIm, 1, &TPES[13]);
  
  /* create the transition operator for ipes */
  
  jIm[0]=0;   jIm[1]=0;
  jlist[0]=11; jlist[1]=12; jRe[0]=-sqrt(6.0/7.0); jRe[1]= sqrt(1.0/7.0);
  CreateOperatorCreateList(38, jlist, jRe, jIm, 2, &TIPES[0]);
  jlist[0]=13; jlist[1]=14; jRe[0]=-sqrt(5.0/7.0); jRe[1]= sqrt(2.0/7.0);
  CreateOperatorCreateList(38, jlist, jRe, jIm, 2, &TIPES[1]);
  jlist[0]=15; jlist[1]=16; jRe[0]=-sqrt(4.0/7.0); jRe[1]= sqrt(3.0/7.0);
  CreateOperatorCreateList(38, jlist, jRe, jIm, 2, &TIPES[2]);
  jlist[0]=17; jlist[1]=18; jRe[0]=-sqrt(3.0/7.0); jRe[1]= sqrt(4.0/7.0);
  CreateOperatorCreateList(38, jlist, jRe, jIm, 2, &TIPES[3]);
  jlist[0]=19; jlist[1]=20; jRe[0]=-sqrt(2.0/7.0); jRe[1]= sqrt(5.0/7.0);
  CreateOperatorCreateList(38, jlist, jRe, jIm, 2, &TIPES[4]);
  jlist[0]=21; jlist[1]=22; jRe[0]=-sqrt(1.0/7.0); jRe[1]= sqrt(6.0/7.0);
  CreateOperatorCreateList(38, jlist, jRe, jIm, 2, &TIPES[5]);
  jlist[0]=10;              jRe[0]=1.0;
  CreateOperatorCreateList(38, jlist, jRe, jIm, 1, &TIPES[6]);
  jlist[0]=11; jlist[1]=12; jRe[0]= sqrt(1.0/7.0); jRe[1]= sqrt(6.0/7.0);
  CreateOperatorCreateList(38, jlist, jRe, jIm, 2, &TIPES[7]);
  jlist[0]=13; jlist[1]=14; jRe[0]= sqrt(2.0/7.0); jRe[1]= sqrt(5.0/7.0);
  CreateOperatorCreateList(38, jlist, jRe, jIm, 2, &TIPES[8]);
  jlist[0]=15; jlist[1]=16; jRe[0]= sqrt(3.0/7.0); jRe[1]= sqrt(4.0/7.0);
  CreateOperatorCreateList(38, jlist, jRe, jIm, 2, &TIPES[9]);
  jlist[0]=17; jlist[1]=18; jRe[0]= sqrt(4.0/7.0); jRe[1]= sqrt(3.0/7.0);
  CreateOperatorCreateList(38, jlist, jRe, jIm, 2, &TIPES[10]);
  jlist[0]=19; jlist[1]=20; jRe[0]= sqrt(5.0/7.0); jRe[1]= sqrt(2.0/7.0);
  CreateOperatorCreateList(38, jlist, jRe, jIm, 2, &TIPES[11]);
  jlist[0]=21; jlist[1]=22; jRe[0]= sqrt(6.0/7.0); jRe[1]= sqrt(1.0/7.0);
  CreateOperatorCreateList(38, jlist, jRe, jIm, 2, &TIPES[12]);
  jlist[0]=23;              jRe[0]=1.0;
  CreateOperatorCreateList(38, jlist, jRe, jIm, 1, &TIPES[13]);
  
  /* create the transition operator for core level pes */
  
  jIm[0]=0;   jIm[1]=0;
  jlist[0]=1; jlist[1]=2; jRe[0]=-sqrt(4.0/5.0); jRe[1]= sqrt(1.0/5.0);
  CreateOperatorAnnihilateList(38, jlist, jRe, jIm, 2, &TCPES[0]);
  jlist[0]=3; jlist[1]=4; jRe[0]=-sqrt(3.0/5.0); jRe[1]= sqrt(2.0/5.0);
  CreateOperatorAnnihilateList(38, jlist, jRe, jIm, 2, &TCPES[1]);
  jlist[0]=5; jlist[1]=6; jRe[0]=-sqrt(2.0/5.0); jRe[1]= sqrt(3.0/5.0);
  CreateOperatorAnnihilateList(38, jlist, jRe, jIm, 2, &TCPES[2]);
  jlist[0]=7; jlist[1]=8; jRe[0]=-sqrt(1.0/5.0); jRe[1]= sqrt(4.0/5.0);
  CreateOperatorAnnihilateList(38, jlist, jRe, jIm, 2, &TCPES[3]);
  jlist[0]=0;             jRe[0]=1.0;
  CreateOperatorAnnihilateList(38, jlist, jRe, jIm, 1, &TCPES[4]);
  jlist[0]=1; jlist[1]=2; jRe[0]= sqrt(1.0/5.0); jRe[1]= sqrt(4.0/5.0);
  CreateOperatorAnnihilateList(38, jlist, jRe, jIm, 2, &TCPES[5]);
  jlist[0]=3; jlist[1]=4; jRe[0]= sqrt(2.0/5.0); jRe[1]= sqrt(3.0/5.0);
  CreateOperatorAnnihilateList(38, jlist, jRe, jIm, 2, &TCPES[6]);
  jlist[0]=5; jlist[1]=6; jRe[0]= sqrt(3.0/5.0); jRe[1]= sqrt(2.0/5.0);
  CreateOperatorAnnihilateList(38, jlist, jRe, jIm, 2, &TCPES[7]);
  jlist[0]=7; jlist[1]=8; jRe[0]= sqrt(4.0/5.0); jRe[1]= sqrt(1.0/5.0);
  CreateOperatorAnnihilateList(38, jlist, jRe, jIm, 2, &TCPES[8]);
  jlist[0]=9;             jRe[0]=1.0;
  CreateOperatorAnnihilateList(38, jlist, jRe, jIm, 1, &TCPES[9]);
  
  /* create the transition operator for XAS */
  
  CreateOperatorCF(38, 3, 2, Index4Fup, Index4Fdn, Index3Dup, Index3Ddn, Vx, &TXASx);
  CreateOperatorCF(38, 3, 2, Index4Fup, Index4Fdn, Index3Dup, Index3Ddn, Vy, &TXASy);
  CreateOperatorCF(38, 3, 2, Index4Fup, Index4Fdn, Index3Dup, Index3Ddn, Vz, &TXASz);
  CreateOperatorCF(38, 3, 2, Index4Fup, Index4Fdn, Index3Dup, Index3Ddn, Vleft, &TXASl);
  CreateOperatorCF(38, 3, 2, Index4Fup, Index4Fdn, Index3Dup, Index3Ddn, Vright, &TXASr);
  EndTiming("define operators", &Time);
  
  printf("\n========================================================");
  printf("\n==========         find groundstate          ===========");
  printf("\n========================================================\n");
  BeginTiming("find groundstate", &Time);

  NBitsKey=16;
  CreateString[ 0]= 0;
  CreateString[ 1]= 1;
  CreateString[ 2]= 2;
  CreateString[ 3]= 3;
  CreateString[ 4]= 4;
  CreateString[ 5]= 5;
  CreateString[ 6]= 6;
  CreateString[ 7]= 7;
  CreateString[ 8]= 8;
  CreateString[ 9]= 9;
  CreateString[10]=24;
  CreateString[11]=25;
  CreateString[12]=26;
  CreateString[13]=27;
  CreateString[14]=28;
  CreateString[15]=29;
  CreateString[16]=30;
  CreateString[17]=31;
  CreateString[18]=32;
  CreateString[19]=33;
  CreateString[20]=34;
  CreateString[21]=35;
  CreateString[22]=36;
  CreateString[23]=37;
  CreateWaveFunctionFromString(&psi, CreateString, 24, 38, 0, NBitsKey);
  
  BeginTiming("LanczosGroundState", &Time);
  LanczosGroundStateRestricted(OppH, &psi, Restrictions, &epsilon, zero, 0, &psi, 0);
  EndTiming("LanczosGroundState", &Time);
  
  CalculatePrefactorHistogram(Hist, psi, 2.0, 30);
  printf("\nPrefactor Histogram of groundstate\n");
  printf("Interval                                 | percentage of states   | alpha^2 * number of states \n");
  for(i=0;i<30;i++)
  {
    printf("2^-%2i (%10.3E) to 2^-%2i (%10.3E) | %22.15E | %22.15E \n",i,pow(0.5,(double)i),i+1,pow(0.5,(double)(i+1)),Hist[i],Hist[i]*pow(0.5,(double)(2*i)));
  }
  
  CalculateOccupationHistogram(OccHist, psi, Index4F, 14);
  
  PsiHermitianOperatorPsi(OppH, psi, &energy);
  PsiHermitianOperatorSqrPsi(OppH, psi, &energysqr);
  PsiHermitianOperatorPsi(Oppldots4f, psi, &Esoc);
  Esoc*=zeta4f;
  PsiHermitianOperatorPsi(OppUff, psi, &ECoulomb);
  PsiHermitianOperatorPsi(OppHyb, psi, &EKin);
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
  PsiHermitianOperatorPsi(OppNf, psi, &Nf);
  PsiHermitianOperatorPsi(OppNfivehalf, psi, &Nfivehalf);
  PsiHermitianOperatorPsi(OppNsevenhalf, psi, &Nsevenhalf);
  PsiHermitianOperatorSqrPsi(OppNf, psi, &Nfsqr);
  PsiHermitianOperatorSqrPsi(OppNfivehalf, psi, &Nfivehalfsqr);
  PsiHermitianOperatorSqrPsi(OppNsevenhalf, psi, &Nsevenhalfsqr);
  EndTiming("find groundstate", &Time);
  PrintTiming(Time);
  
  printf("\n========================================================");
  printf("\n==========   Groundstte expectation values   ===========");
  printf("\n========================================================\n");
  printf("Energy              %22.15E  (variance %22.15E)\n",energy,energysqr-energy*energy);
  printf("Kinetic energy      %22.15E\n",EKin);
  printf("Relativistic energy %22.15E\n",Esoc);
  printf("Coulomb energy (4f) %22.15E\n",ECoulomb);
  printf("Sx %22.15E Sy %22.15E Sz %22.15E Ssqr %22.15E\n",Sx,Sy,Sz,Ssqr);
  printf("Lx %22.15E Ly %22.15E Lz %22.15E Lsqr %22.15E\n",Lx,Ly,Lz,Lsqr);
  printf("Jx %22.15E Jy %22.15E Jz %22.15E Jsqr %22.15E\n",Jx,Jy,Jz,Jsqr);
  printf("Nf   %22.15E <Nf^2>   %22.15E\n",Nf,Nfsqr);
  printf("N5/2 %22.15E <N5/2^2> %22.15E\n",Nfivehalf,Nfivehalfsqr);
  printf("N7/2 %22.15E <N7/2^2> %22.15E\n",Nsevenhalf,Nsevenhalfsqr);
  printf("\n");
  
  printf("Occupation of f-shell\n");
  for(i=0;i<=14;i++)
  {
    printf("f^%2i | %22.15E\n",i,OccHist[i]);
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
  
  if(doPES)
  {
    printf("\n========================================================");
    printf("\n==========                PES                ===========");
    printf("\n========================================================\n");
    BeginTiming("PES", &Time);
    
    for(i=0;i<14;i++)
    {
      Spectra[i].N=6000;
      Spectra[i].Emin=-10.0;
      Spectra[i].Emax= 50.0;
      Spectra[i].dE = (Spectra[i].Emax - Spectra[i].Emin) / Spectra[i].N;
      Spectra[i].Gamma = 0.08;
      SpectrumAllocate(Spectra[i], &Spectra[i]);
    }
    
    for(i=0;i<14;i++)
    {
      BeginTiming("spectrum", &Time);
      OperatorPsi(TPES[i], psi, &Opsi);
      LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[i].Gamma);
      TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[i]);
      EndTiming("spectrum", &Time);
    }
    
    FilePrintSpectra(&Spectra[0], 14, "CeLigandField_PES","w");
    EndTiming("PES", &Time);
    PrintTiming(Time);
  }
  
  if(doIPES)
  {
    printf("\n========================================================");
    printf("\n==========               IPES                ===========");
    printf("\n========================================================\n");
    BeginTiming("IPES", &Time);
    
    for(i=14;i<28;i++)
    {
      Spectra[i].N=6000;
      Spectra[i].Emin=-10.0;
      Spectra[i].Emax= 50.0;
      Spectra[i].dE = (Spectra[i].Emax - Spectra[i].Emin) / Spectra[i].N;
      Spectra[i].Gamma = 0.08;
      SpectrumAllocate(Spectra[i], &Spectra[i]);
    }
    
    for(i=0;i<14;i++)
    {
      BeginTiming("spectrum", &Time);
      OperatorPsi(TIPES[i], psi, &Opsi);
      LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[14+i].Gamma);
      TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[14+i]);
      EndTiming("spectrum", &Time);
    }
    
    FilePrintSpectra(&Spectra[14], 14, "CeLigandField_IPES","w");
    EndTiming("IPES", &Time);
    PrintTiming(Time);
  }
  
  if(doCPES)
  {
    printf("\n========================================================");
    printf("\n==========               CPES                ===========");
    printf("\n========================================================\n");
    BeginTiming("CPES", &Time);
    
    for(i=28;i<38;i++)
    {
      Spectra[i].N=7000;
      Spectra[i].Emin=-40.0;
      Spectra[i].Emax= 30.0;
      Spectra[i].dE = (Spectra[i].Emax - Spectra[i].Emin) / Spectra[i].N;
      Spectra[i].Gamma = 0.08;
      SpectrumAllocate(Spectra[i], &Spectra[i]);
    }
    
    for(i=0;i<10;i++)
    {
      BeginTiming("spectrum", &Time);
      OperatorPsi(TCPES[i], psi, &Opsi);
      LanczosTriDiagonalizeRestricted(OppH, Opsi, Restrictions, epsilon, &SpecTri,  Spectra[28+i].Gamma);
      TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[28+i]);
      EndTiming("spectrum", &Time);
    }
    
    FilePrintSpectra(&Spectra[28], 10, "CeLigandField_CPES","w");
    EndTiming("CPES", &Time);
    PrintTiming(Time);
  }
  
  if(doXAS)
  {
    printf("\n========================================================");
    printf("\n==========                XAS                ===========");
    printf("\n========================================================\n");
    BeginTiming("XAS", &Time);
    
    for(i=38;i<43;i++)
    {
      Spectra[i].N=7000;
      Spectra[i].Emin=-40.0;
      Spectra[i].Emax= 30.0;
      Spectra[i].dE = (Spectra[i].Emax - Spectra[i].Emin) / Spectra[i].N;
      Spectra[i].Gamma = 0.08;
      SpectrumAllocate(Spectra[i], &Spectra[i]);
    }
    
    BeginTiming("spectrum", &Time);
    OperatorPsi(TXASx, psi, &Opsi);
    LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[38].Gamma);
    TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[38]);
    EndTiming("spectrum", &Time);
    
    BeginTiming("spectrum", &Time);
    OperatorPsi(TXASy, psi, &Opsi);
    LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[39].Gamma);
    TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[39]);
    EndTiming("spectrum", &Time);
    
    BeginTiming("spectrum", &Time);
    OperatorPsi(TXASz, psi, &Opsi);
    LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[40].Gamma);
    TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[40]);
    EndTiming("spectrum", &Time);
    
    BeginTiming("spectrum", &Time);
    OperatorPsi(TXASl, psi, &Opsi);
    LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[41].Gamma);
    TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[41]);
    EndTiming("spectrum", &Time);
    
    BeginTiming("spectrum", &Time);
    OperatorPsi(TXASr, psi, &Opsi);
    LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[42].Gamma);
    TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[42]);
    EndTiming("spectrum", &Time);
    
    FilePrintSpectra(&Spectra[38], 5, "CeLigandField_XAS","w");
    EndTiming("XAS", &Time);
    PrintTiming(Time);
  }
  
  if(doXiSpin)
  {
    printf("\n========================================================");
    printf("\n==========               Spin                ===========");
    printf("\n========================================================\n");
    BeginTiming("Spin", &Time);
    
    for(i=43;i<46;i++)
    {
      Spectra[i].N=6000;
      Spectra[i].Emin=-1.0;
      Spectra[i].Emax= 11.0;
      Spectra[i].dE = (Spectra[i].Emax - Spectra[i].Emin) / Spectra[i].N;
      Spectra[i].Gamma = 0.01;
      SpectrumAllocate(Spectra[i], &Spectra[i]);
    }
    
    BeginTiming("spectrum", &Time);
    OperatorPsi(OppSx, psi, &Opsi);
    LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[43].Gamma);
    TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[43]);
    EndTiming("spectrum", &Time);
    
    BeginTiming("spectrum", &Time);
    OperatorPsi(OppSy, psi, &Opsi);
    LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[44].Gamma);
    TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[44]);
    EndTiming("spectrum", &Time);
    
    BeginTiming("spectrum", &Time);
    OperatorPsi(OppSz, psi, &Opsi);
    LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[45].Gamma);
    TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[45]);
    EndTiming("spectrum", &Time);
    
    FilePrintSpectra(&Spectra[43], 3, "CeLigandField_Spin","w");
    
    EndTiming("Spin", &Time);
    PrintTiming(Time);
  }
  
  if(doXiCharge)
  {
    printf("\n========================================================");
    printf("\n==========              charge               ===========");
    printf("\n========================================================\n");
    BeginTiming("Charge", &Time);
    
    Spectra[46].N=6000;
    Spectra[46].Emin=-1.0;
    Spectra[46].Emax= 11.0;
    Spectra[46].dE = (Spectra[46].Emax - Spectra[46].Emin) / Spectra[46].N;
    Spectra[46].Gamma = 0.01;
    SpectrumAllocate(Spectra[46], &Spectra[46]);
    
    OperatorCopy(&TCharge, OppNf);
    OperatorAddLader(&TCharge, &oadd, 0, -1, 0);
    
    BeginTiming("spectrum", &Time);
    OperatorPsi(TCharge, psi, &Opsi);
    LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[46].Gamma);
    TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[46]);
    EndTiming("spectrum", &Time);
    
    FilePrintSpectra(&Spectra[46], 1, "CeLigandField_Charge","w");
    
    EndTiming("Charge", &Time);
    PrintTiming(Time);
  }
  
  
  EndTiming("Total", &Time);
  PrintTiming(Time);
}