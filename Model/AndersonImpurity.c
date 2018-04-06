#include "AndersonImpurity.h"

void CeImpurity()
{
  TimeType Time;
  InitializeTiming(&Time);
  BeginTiming("Total", &Time);
  printf("\n========================================================");
  printf("\n==========      start andersen impurity      ===========");
  printf("\n========================================================\n");
  uint32_t i,j;
  double Bmag[3], Hex[3], zeta3d, zeta4f;
  double Rydberg, red3D4F, red4F4F, Fff[4], Fdf[3], Gdf[3];
  double Uff, Udf, Delta, epsilonf;
  uint16_t Index3Ddn[5]={0,2,4,6,8};
  uint16_t Index3Dup[5]={1,3,5,7,9};
  uint16_t Index4Fdn[7]={10,12,14,16,18,20,22};
  uint16_t Index4Fup[7]={11,13,15,17,19,21,23};
  uint16_t IndexLFdn[7]={24,26,28,30,32,34,36};
  uint16_t IndexLFup[7]={25,27,29,31,33,35,37};
  uint16_t IndexLFvcup[7], IndexLFvcdn[7], IndexLFvcnup[7], IndexLFvcndn[7];
  uint16_t Index4F[14]={10,11,12,13,14,15,16,17,18,19,20,21,22,23};
  uint16_t IndexLF[14]={24,25,26,27,28,29,30,31,32,33,34,35,36,37};
  uint16_t IndexF[14];
  uint16_t jlist[14];
  uint16_t IndexCrTemp[16],IndexAnTemp[16];
  double OneList[]= {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
  double ZeroList[]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double jRe[14],jIm[14];
  HarmonicsExpansion SphericalAkm,Vx,Vy,Vz,Vleft,Vright;
  int Vl[]={1,1,1};
  int Vm[]={-1,0,1};
  double VxRe[]={sqrt(0.5),0.0, sqrt(0.5)};
  double VyRe[]={sqrt(0.5),0.0,-sqrt(0.5)};
  double VzRe[]={0.0,1.0,0.0};
  double VlRe[]={1.0,0.0,0.0};
  double VrRe[]={0.0,0.0,1.0};
  OperatorType Opptemp;
  OperatorType OppH, OppHpart, OppSx, OppSy, OppSz, OppSsqr, OppLx, OppLy, OppLz, OppLsqr, OppJx, OppJy, OppJz, OppJsqr, Oppldots4f, Oppldots3d, OppUff, OppUdf, OppDelta, OppHyb, OppHyb52, OppHyb72, OppNf, OppNL, *OppNc, *OppNv, OppNfivehalf, OppNsevenhalf;
  OperatorType TPES[14], TIPES[14], TCPES[10], TXASx, TXASy, TXASz, TXASl, TXASr, TCharge;
  OperatorType OppHybIL,OppHybIv,OppHybIc,OppHybLv,OppHybLc,*OppHybc,*OppHybv;
  OperatorType OppHybIL52,OppHybIv52,OppHybIc52,OppHybLv52,OppHybLc52,*OppHybc52,*OppHybv52;
  OperatorType OppHybIL72,OppHybIv72,OppHybIc72,OppHybLv72,OppHybLc72,*OppHybc72,*OppHybv72;
  OperatorType OppNL52, OppNL72, *OppNc52, *OppNv52, *OppNc72, *OppNv72;
  double epsilon, zero;
  double energy, EKin, Esoc, ECoulomb, Sx, Sy, Sz, Ssqr, Lx, Ly, Lz, Lsqr, Jx, Jy, Jz, Jsqr, Nf, Nfivehalf, Nsevenhalf, energysqr, Nfsqr, Nfivehalfsqr, Nsevenhalfsqr;
  uint32_t NBitsKey;
  WaveFunctionType psi, Opsi;
  uint16_t *CreateString;
  TriDiagonalMatrixType SpecTri;
  SpectraType Spectra[47];
  MatrixType Density;
  double Hist[40],OccHist[15];
  uint16_t oadd=0;
  AndersonMatrixType GBath;
  BADoubleTriDiagonalMatrixType GBathBA, GBathBA52, GBathBA72;
  uint32_t NBathv, NBathc, NFermion;
  double NfStart, Nf52Start, Nf72Start, MaxNf, t=1.0;
  double occupation, occupation52, occupation72, onsiteshift, onsiteshift52, onsiteshift72, kbTeff;
  DeterminantRestrictionType Restrictions;
  uint32_t doPES, doIPES, doCPES, doCPESShort, doXAS, doXiSpin, doXiCharge, NGroundStateLoops, iGroundStateLoops;
  MatrixType jjzRotationMatrix;
  
  /* allocate space for restrictions */
  Restrictions.N=1;
  Restrictions.Min=malloc(Restrictions.N*sizeof(uint32_t));
  Restrictions.Max=malloc(Restrictions.N*sizeof(uint32_t));
  Restrictions.Determinant=malloc(Restrictions.N*sizeof(uint8_t*));
  for(i=0;i<Restrictions.N;i++)
  {
    Restrictions.Determinant[i]=malloc(16*sizeof(uint8_t)); /* we have 38+? fermions in the basis, Determinants are coded per byte i.e. blocks of 8 fermions. 8*16=128 and can store the fermions (note that after transformation the number of sites is less than Nbath) */
  }
  
  
  printf("\n========================================================");
  printf("\n==========           Read Gbath              ===========");
  printf("\n========================================================\n");
  
  FileReadAndersonMatrixTypeUndeterminedPlusMu(&GBath, "CeGBath", &GBath.mu);
  
  printf("\n========================================================");
  printf("\n==========          set parameters           ===========");
  printf("\n========================================================\n");
  
  doPES=0;
  doIPES=0;
  doCPES=0;
  doCPESShort=1;
  doXAS=0;
  doXiSpin=0;
  doXiCharge=0;
  NGroundStateLoops=3;
  
  Rydberg=13.6057; /* calculations are done in eV this alows one to imput in mixed eV and Rydberg */
  
  Bmag[0]=0.000000; /* magnetic field in the x direction (units of eV) */
  Bmag[1]=0.000000; /* magnetic field in the y direction (units of eV) */
  Bmag[2]=0.000000; /* magnetic field in the z direction (units of eV) */
  Hex[0]= 0.000000; /* exchange field in the x direction (units of eV) */
  Hex[1]= 0.000000; /* exchange field in the y direction (units of eV) */
  Hex[2]= 0.000000; /* exchange field in the z direction (units of eV) */
  
  zeta3d=0.54721 * Rydberg;       /* spin orbit of the d core shell */
  zeta4f=0.00789 * Rydberg;       /* spin orbit of the f shell */
  
  red3D4F=0.76;    /* reduction of 3d-4f slater integrals from atomic HF values */
  red4F4F=0.57;     /* reduction of 4f-4f slater integrals from atomic HF values */
  
  Uff=8.0;          /* Uff - shperical part of the Coulmb interactin */
  Udf=10.3;         /* Udf - shperical part of the Coulmb interactin */
  Delta=-1.9;       /* energy differnece between f0 and f1 ligand hole ???? */
  epsilonf = Delta - 10 * Udf; /* onsite energy of f orbital. Onsite energy of ligand orbitals given by GBath */
  
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
  
  NBathv = 2;
  NBathc = 0;
  MaxNf = 1;
  
  NfStart = 0.8344785458087937;     /* numbers obtained with MaxMf = 1 */
  Nf52Start = 0.6029731957814112;
  Nf72Start = 0.2315053500273825;
  
  /*
    Nf52Start = NfStart * (6.0/14.0);
    Nf72Start = NfStart * (8.0/14.0);
  */
  /*
   Nf52Start = NfStart;
   Nf72Start = 0.000001;
  */
  
  NBitsKey=16;
  
  epsilon = 0.0001;           /* smallest value stored as prefactor in wavefunction <== enlarge to speed up calculation !!!                       */
  zero = 10 * DBL_EPSILON;        /* smallest value to be considered different from zero in convergnece tests <== enlarge to speed up calculation !!! */
  SpecTri.N=50;                  /* number of poles in spectra  <== reduce here to speed up calculation of spectra !!!!                              */
  InitTriDiagonalMatrix(&SpecTri);
  
  
  /* some settings one should not change */
  
  NFermion=10+14+14+14*(NBathv+NBathc);
  
  CreateString=calloc(10+14+14*NBathv, sizeof(uint16_t));
  for(i=0;i<10;i++)
  {
    CreateString[i]=(uint16_t)i;
  }
  for(i=10;i<(10+14*(NBathv+1));i++)
  {
    CreateString[i]=(uint16_t)i+14;
  }
  
  if(NBathc!=0)
  {
    OppHybc=malloc(NBathc*sizeof(OperatorType));
    OppHybc52=malloc(NBathc*sizeof(OperatorType));
    OppHybc72=malloc(NBathc*sizeof(OperatorType));
    OppNc=malloc(NBathc*sizeof(OperatorType));
    OppNc52=malloc(NBathc*sizeof(OperatorType));
    OppNc72=malloc(NBathc*sizeof(OperatorType));
  }
  if(NBathv!=0)
  {
    OppHybv=malloc(NBathv*sizeof(OperatorType));
    OppHybv52=malloc(NBathv*sizeof(OperatorType));
    OppHybv72=malloc(NBathv*sizeof(OperatorType));
    OppNv=malloc(NBathv*sizeof(OperatorType));
    OppNv52=malloc(NBathv*sizeof(OperatorType));
    OppNv72=malloc(NBathv*sizeof(OperatorType));
  }
  
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
  Restrictions.Max[0]=MaxNf;
  Restrictions.Determinant[0][0] = 0x00; /* bit 1 to 8 are set to zero (d states) */
  Restrictions.Determinant[0][1] = 0x3F; /* bit 9 and 10 are set to zero, (d states) bit 11 to 16 are set to one (f states) */
  Restrictions.Determinant[0][2] = 0xFF; /* bit 17 to 24 are set to one (f states) */
  Restrictions.Determinant[0][3] = 0x00; /* bit 25 to 32 are set to zero (Lf states) */
  Restrictions.Determinant[0][4] = 0x00;
  Restrictions.Determinant[0][5] = 0x00;
  Restrictions.Determinant[0][6] = 0x00;
  Restrictions.Determinant[0][7] = 0x00;
  Restrictions.Determinant[0][8] = 0x00;
  Restrictions.Determinant[0][9] = 0x00;
  Restrictions.Determinant[0][10] = 0x00;
  Restrictions.Determinant[0][11] = 0x00;
  Restrictions.Determinant[0][12] = 0x00;
  Restrictions.Determinant[0][13] = 0x00;
  Restrictions.Determinant[0][14] = 0x00;
  Restrictions.Determinant[0][15] = 0x00;
  
  
  printf("\n========================================================");
  printf("\n==========         define operators          ===========");
  printf("\n========================================================\n");
  BeginTiming("define operators", &Time);
  
  /* define some operators */
  
  CreateOperatorSx(NFermion, 7, Index4Fup, Index4Fdn, &OppSx);
  CreateOperatorSy(NFermion, 7, Index4Fup, Index4Fdn, &OppSy);
  CreateOperatorSz(NFermion, 7, Index4Fup, Index4Fdn, &OppSz);
  CreateOperatorSsqr(NFermion, 7, Index4Fup, Index4Fdn, &OppSsqr);
  CreateOperatorLx(NFermion, 3, Index4Fup, Index4Fdn, &OppLx);
  CreateOperatorLy(NFermion, 3, Index4Fup, Index4Fdn, &OppLy);
  CreateOperatorLz(NFermion, 3, Index4Fup, Index4Fdn, &OppLz);
  CreateOperatorLsqr(NFermion, 3, Index4Fup, Index4Fdn, &OppLsqr);
  CreateOperatorJx(NFermion, 3, Index4Fup, Index4Fdn, &OppJx);
  CreateOperatorJy(NFermion, 3, Index4Fup, Index4Fdn, &OppJy);
  CreateOperatorJz(NFermion, 3, Index4Fup, Index4Fdn, &OppJz);
  CreateOperatorJsqr(NFermion, 3, Index4Fup, Index4Fdn, &OppJsqr);
  CreateOperatorldots(NFermion, 3, Index4Fup, Index4Fdn, &Oppldots4f);
  CreateOperatorldots(NFermion, 2, Index3Dup, Index3Ddn, &Oppldots3d);
  CreateOperatorCoulombSingleShell(NFermion, 3, Index4Fup, Index4Fdn, Fff, &OppUff);
  CreateOperatorCoulombTwoShellsConstantOccupation(NFermion, 2, 3, Index3Dup, Index3Ddn, Index4Fup, Index4Fdn, Fdf, Gdf, &OppUdf);
  CreateOperatorNumberShell(NFermion, 14, Index4F, &OppNf);
  CreateOperatorNumberShell(NFermion, 14, IndexLF, &OppNL);
  OperatorCopy(&OppNfivehalf, Oppldots4f);
  OperatorAdd(&OppNfivehalf, OppNf, -1.5, 0.0);
  OperatorScale(&OppNfivehalf,-0.285714285714285714,0.0);
  OperatorCopy(&OppNsevenhalf, OppNf);
  OperatorAdd(&OppNsevenhalf, OppNfivehalf, -1.0, 0.0);
  OperatorCopy(&OppDelta, OppNf);
  OperatorScale(&OppDelta, epsilonf, 0);
  
  /* summ these operators to create the Hamiltonian */
  
  OperatorCopy(&OppHpart, OppUff);
  OperatorAdd(&OppHpart, OppUdf, 1.0, 0.0);
  OperatorAdd(&OppHpart, Oppldots3d, zeta3d, 0.0);
  OperatorAdd(&OppHpart, Oppldots4f, zeta4f, 0.0);
  OperatorAdd(&OppHpart, OppDelta, 1.0, 0.0);
  
  /* define hopping operators */
  
  CreateOperatorCF(NFermion, 3, 3, IndexLFup, IndexLFdn, Index4Fup, Index4Fdn, SphericalAkm, &OppHybIL);
  CreateOperatorCF(NFermion, 3, 3, Index4Fup, Index4Fdn, IndexLFup, IndexLFdn, SphericalAkm, &Opptemp);
  OperatorAdd(&OppHybIL, Opptemp, 1.0, 0.0);
  FreeOperator(&Opptemp);
  
  if(NBathv!=0)
  {
    for(i=0;i<14;i++)
    {
      IndexF[i]=10+14+14+i;
    }
    CreateOperatorNumberShell(NFermion, 14, IndexF, &OppNv[0]);
    for(i=0;i<7;i++)
    {
      IndexLFvcdn[i]=10+14+14+2*i;
      IndexLFvcup[i]=10+14+14+2*i+1;
    }
    CreateOperatorCF(NFermion, 3, 3, IndexLFvcup, IndexLFvcdn, Index4Fup, Index4Fdn, SphericalAkm, &OppHybIv);
    CreateOperatorCF(NFermion, 3, 3, Index4Fup, Index4Fdn, IndexLFvcup, IndexLFvcdn, SphericalAkm, &Opptemp);
    OperatorAdd(&OppHybIv, Opptemp, 1.0, 0.0);
    FreeOperator(&Opptemp);
    CreateOperatorCF(NFermion, 3, 3, IndexLFvcup, IndexLFvcdn, IndexLFup, IndexLFdn, SphericalAkm, &OppHybLv);
    CreateOperatorCF(NFermion, 3, 3, IndexLFup, IndexLFdn, IndexLFvcup, IndexLFvcdn, SphericalAkm, &Opptemp);
    OperatorAdd(&OppHybLv, Opptemp, 1.0, 0.0);
    FreeOperator(&Opptemp);
  }
  for(j=1;j<NBathv;j++)
  {
    for(i=0;i<14;i++)
    {
      IndexF[i]=10+14+14+14*j+i;
    }
    CreateOperatorNumberShell(NFermion, 14, IndexF, &OppNv[j]);
    for(i=0;i<7;i++)
    {
      IndexLFvcdn[i]=10+14+14+14*(j-1)+2*i;
      IndexLFvcup[i]=10+14+14+14*(j-1)+2*i+1;
      IndexLFvcndn[i]=10+14+14+14*j+2*i;
      IndexLFvcnup[i]=10+14+14+14*j+2*i+1;
    }
    CreateOperatorCF(NFermion, 3, 3, IndexLFvcup, IndexLFvcdn, IndexLFvcnup, IndexLFvcndn, SphericalAkm, &OppHybv[j-1]);
    CreateOperatorCF(NFermion, 3, 3, IndexLFvcnup, IndexLFvcndn, IndexLFvcup, IndexLFvcdn, SphericalAkm, &Opptemp);
    OperatorAdd(&OppHybv[j-1], Opptemp, 1.0, 0.0);
    FreeOperator(&Opptemp);
    
  }
  if(NBathc!=0)
  {
    for(i=0;i<14;i++)
    {
      IndexF[i]=10+14+14+14*NBathv+i;
    }
    CreateOperatorNumberShell(NFermion, 14, IndexF, &OppNc[0]);
    for(i=0;i<7;i++)
    {
      IndexLFvcdn[i]=10+14+14+14*NBathv+2*i;
      IndexLFvcup[i]=10+14+14+14*NBathv+2*i+1;
    }
    CreateOperatorCF(NFermion, 3, 3, IndexLFvcup, IndexLFvcdn, Index4Fup, Index4Fdn, SphericalAkm, &OppHybIc);
    CreateOperatorCF(NFermion, 3, 3, Index4Fup, Index4Fdn, IndexLFvcup, IndexLFvcdn, SphericalAkm, &Opptemp);
    OperatorAdd(&OppHybIc, Opptemp, 1.0, 0.0);
    FreeOperator(&Opptemp);
    CreateOperatorCF(NFermion, 3, 3, IndexLFvcup, IndexLFvcdn, IndexLFup, IndexLFdn, SphericalAkm, &OppHybLc);
    CreateOperatorCF(NFermion, 3, 3, IndexLFup, IndexLFdn, IndexLFvcup, IndexLFvcdn, SphericalAkm, &Opptemp);
    OperatorAdd(&OppHybLc, Opptemp, 1.0, 0.0);
    FreeOperator(&Opptemp);
  }
  for(j=1;j<NBathc;j++)
  {
    for(i=0;i<14;i++)
    {
      IndexF[i]=10+14+14+1414*NBathv*j+i;
    }
    CreateOperatorNumberShell(NFermion, 14, IndexF, &OppNc[j]);
    for(i=0;i<7;i++)
    {
      IndexLFvcdn[i]=10+14+14+14*NBathv+14*(j-1)+2*i;
      IndexLFvcup[i]=10+14+14+14*NBathv+14*(j-1)+2*i+1;
      IndexLFvcndn[i]=10+14+14+14*NBathv+14*j+2*i;
      IndexLFvcnup[i]=10+14+14+14*NBathv+14*j+2*i+1;
    }
    CreateOperatorCF(NFermion, 3, 3, IndexLFvcup, IndexLFvcdn, IndexLFvcnup, IndexLFvcndn, SphericalAkm, &OppHybc[j-1]);
    CreateOperatorCF(NFermion, 3, 3, IndexLFvcnup, IndexLFvcndn, IndexLFvcup, IndexLFvcdn, SphericalAkm, &Opptemp);
    OperatorAdd(&OppHybc[j-1], Opptemp, 1.0, 0.0);
    FreeOperator(&Opptemp);
  }
  
   /* define hopping operators for 5/2 on jjz basis */
  
  for(i=0;i<6;i++)
  {
    IndexCrTemp[i]=(uint16_t)(10+14+i);
    IndexAnTemp[i]=(uint16_t)(10+14+i);
  }
  CreateOperatorCreateAnnihilateList(NFermion, IndexCrTemp, IndexAnTemp, OneList, ZeroList, 6, &OppNL52);
  for(i=0;i<6;i++)
  {
    IndexCrTemp[i]=(uint16_t)(10+i);
    IndexAnTemp[i]=(uint16_t)(10+14+i);
    IndexCrTemp[i+6]=(uint16_t)(10+14+i);
    IndexAnTemp[i+6]=(uint16_t)(10+i);
  }
  CreateOperatorCreateAnnihilateList(NFermion, IndexCrTemp, IndexAnTemp, OneList, ZeroList, 12, &OppHybIL52);
  if(NBathv!=0)
  {
    for(i=0;i<6;i++)
    {
      IndexCrTemp[i]=(uint16_t)(10+i);
      IndexAnTemp[i]=(uint16_t)(10+14+14+i);
      IndexCrTemp[i+6]=(uint16_t)(10+14+14+i);
      IndexAnTemp[i+6]=(uint16_t)(10+i);
    }
    CreateOperatorCreateAnnihilateList(NFermion, IndexCrTemp, IndexAnTemp, OneList, ZeroList, 12, &OppHybIv52);
    for(i=0;i<6;i++)
    {
      IndexCrTemp[i]=(uint16_t)(10+14+i);
      IndexAnTemp[i]=(uint16_t)(10+14+14+i);
      IndexCrTemp[i+6]=(uint16_t)(10+14+14+i);
      IndexAnTemp[i+6]=(uint16_t)(10+14+i);
    }
    CreateOperatorCreateAnnihilateList(NFermion, IndexCrTemp, IndexAnTemp, OneList, ZeroList, 12, &OppHybLv52);
  }
  for(j=1;j<NBathv;j++)
  {
    for(i=0;i<6;i++)
    {
      IndexCrTemp[i]=(uint16_t)(10+14+14+14*(j-1)+i);
      IndexAnTemp[i]=(uint16_t)(10+14+14+14*(j-0)+i);
      IndexCrTemp[i+6]=(uint16_t)(10+14+14+14*(j-0)+i);
      IndexAnTemp[i+6]=(uint16_t)(10+14+14+14*(j-1)+i);
    }
    CreateOperatorCreateAnnihilateList(NFermion, IndexCrTemp, IndexAnTemp, OneList, ZeroList, 12, &OppHybv52[j-1]);
  }
  for(j=0;j<NBathv;j++)
  {
    for(i=0;i<6;i++)
    {
      IndexCrTemp[i]=(uint16_t)(10+14+14+14*j+i);
      IndexAnTemp[i]=(uint16_t)(10+14+14+14*j+i);
    }
    CreateOperatorCreateAnnihilateList(NFermion, IndexCrTemp, IndexAnTemp, OneList, ZeroList, 6, &OppNv52[j]);
  }
  if(NBathc!=0)
  {
    for(i=0;i<6;i++)
    {
      IndexCrTemp[i]=(uint16_t)(10+i);
      IndexAnTemp[i]=(uint16_t)(10+14+14+14*NBathv+i);
      IndexCrTemp[i+6]=(uint16_t)(10+14+14+14*NBathv+i);
      IndexAnTemp[i+6]=(uint16_t)(10+i);
    }
    CreateOperatorCreateAnnihilateList(NFermion, IndexCrTemp, IndexAnTemp, OneList, ZeroList, 12, &OppHybIc52);
    for(i=0;i<6;i++)
    {
      IndexCrTemp[i]=(uint16_t)(10+14+14*NBathv+i);
      IndexAnTemp[i]=(uint16_t)(10+14+14+14*NBathv+i);
      IndexCrTemp[i+6]=(uint16_t)(10+14+14+14*NBathv+i);
      IndexAnTemp[i+6]=(uint16_t)(10+14+14*NBathv+i);
    }
    CreateOperatorCreateAnnihilateList(NFermion, IndexCrTemp, IndexAnTemp, OneList, ZeroList, 12, &OppHybLc52);
  }
  for(j=1;j<NBathc;j++)
  {
    for(i=0;i<6;i++)
    {
      IndexCrTemp[i]=(uint16_t)(10+14+14+14+14*NBathv*(j-1)+i);
      IndexAnTemp[i]=(uint16_t)(10+14+14+14+14*NBathv*(j-0)+i);
      IndexCrTemp[i+6]=(uint16_t)(10+14+14+14+14*NBathv*(j-0)+i);
      IndexAnTemp[i+6]=(uint16_t)(10+14+14+14+14*NBathv*(j-1)+i);
    }
    CreateOperatorCreateAnnihilateList(NFermion, IndexCrTemp, IndexAnTemp, OneList, ZeroList, 12, &OppHybc52[j-1]);
  }
  for(j=0;j<NBathc;j++)
  {
    for(i=0;i<6;i++)
    {
      IndexCrTemp[i]=(uint16_t)(10+14+14+14*NBathv+14*j+i);
      IndexAnTemp[i]=(uint16_t)(10+14+14+14*NBathv+14*j+i);
    }
    CreateOperatorCreateAnnihilateList(NFermion, IndexCrTemp, IndexAnTemp, OneList, ZeroList, 6, &OppNc52[j]);
  }
  
  
  /* define hopping operators for 7/2 on jjz basis */
  
  for(i=0;i<8;i++)
  {
    IndexCrTemp[i]=(uint16_t)(16+14+i);
    IndexAnTemp[i]=(uint16_t)(16+14+i);
  }
  CreateOperatorCreateAnnihilateList(NFermion, IndexCrTemp, IndexAnTemp, OneList, ZeroList, 8, &OppNL72);
  for(i=0;i<8;i++)
  {
    IndexCrTemp[i]=(uint16_t)(16+i);
    IndexAnTemp[i]=(uint16_t)(16+14+i);
    IndexCrTemp[i+8]=(uint16_t)(16+14+i);
    IndexAnTemp[i+8]=(uint16_t)(16+i);
  }
  CreateOperatorCreateAnnihilateList(NFermion, IndexCrTemp, IndexAnTemp, OneList, ZeroList, 16, &OppHybIL72);
  if(NBathv!=0)
  {
    for(i=0;i<8;i++)
    {
      IndexCrTemp[i]=(uint16_t)(16+i);
      IndexAnTemp[i]=(uint16_t)(16+14+14+i);
      IndexCrTemp[i+8]=(uint16_t)(16+14+14+i);
      IndexAnTemp[i+8]=(uint16_t)(16+i);
    }
    CreateOperatorCreateAnnihilateList(NFermion, IndexCrTemp, IndexAnTemp, OneList, ZeroList, 16, &OppHybIv72);
    for(i=0;i<8;i++)
    {
      IndexCrTemp[i]=(uint16_t)(16+14+i);
      IndexAnTemp[i]=(uint16_t)(16+14+14+i);
      IndexCrTemp[i+8]=(uint16_t)(16+14+14+i);
      IndexAnTemp[i+8]=(uint16_t)(16+14+i);
    }
    CreateOperatorCreateAnnihilateList(NFermion, IndexCrTemp, IndexAnTemp, OneList, ZeroList, 16, &OppHybLv72);
  }
  for(j=1;j<NBathv;j++)
  {
    for(i=0;i<8;i++)
    {
      IndexCrTemp[i]=(uint16_t)(16+14+14+14*(j-1)+i);
      IndexAnTemp[i]=(uint16_t)(16+14+14+14*(j-0)+i);
      IndexCrTemp[i+8]=(uint16_t)(16+14+14+14*(j-0)+i);
      IndexAnTemp[i+8]=(uint16_t)(16+14+14+14*(j-1)+i);
    }
    CreateOperatorCreateAnnihilateList(NFermion, IndexCrTemp, IndexAnTemp, OneList, ZeroList, 16, &OppHybv72[j-1]);
  }
  for(j=0;j<NBathv;j++)
  {
    for(i=0;i<8;i++)
    {
      IndexCrTemp[i]=(uint16_t)(16+14+14+14*j+i);
      IndexAnTemp[i]=(uint16_t)(16+14+14+14*j+i);
    }
    CreateOperatorCreateAnnihilateList(NFermion, IndexCrTemp, IndexAnTemp, OneList, ZeroList, 8, &OppNv72[j]);
  }
  if(NBathc!=0)
  {
    for(i=0;i<8;i++)
    {
      IndexCrTemp[i]=(uint16_t)(16+i);
      IndexAnTemp[i]=(uint16_t)(16+14+14+14*NBathv+i);
      IndexCrTemp[i+8]=(uint16_t)(16+14+14+14*NBathv+i);
      IndexAnTemp[i+8]=(uint16_t)(16+i);
    }
    CreateOperatorCreateAnnihilateList(NFermion, IndexCrTemp, IndexAnTemp, OneList, ZeroList, 16, &OppHybIc72);
    for(i=0;i<8;i++)
    {
      IndexCrTemp[i]=(uint16_t)(16+14+14*NBathv+i);
      IndexAnTemp[i]=(uint16_t)(16+14+14+14*NBathv+i);
      IndexCrTemp[i+8]=(uint16_t)(16+14+14+14*NBathv+i);
      IndexAnTemp[i+8]=(uint16_t)(16+14+14*NBathv+i);
    }
    CreateOperatorCreateAnnihilateList(NFermion, IndexCrTemp, IndexAnTemp, OneList, ZeroList, 16, &OppHybLc72);
  }
  for(j=1;j<NBathc;j++)
  {
    for(i=0;i<8;i++)
    {
      IndexCrTemp[i]=(uint16_t)(16+14+14+14+14*NBathv*(j-1)+i);
      IndexAnTemp[i]=(uint16_t)(16+14+14+14+14*NBathv*(j-0)+i);
      IndexCrTemp[i+8]=(uint16_t)(16+14+14+14+14*NBathv*(j-0)+i);
      IndexAnTemp[i+8]=(uint16_t)(16+14+14+14+14*NBathv*(j-1)+i);
    }
    CreateOperatorCreateAnnihilateList(NFermion, IndexCrTemp, IndexAnTemp, OneList, ZeroList, 16, &OppHybc72[j-1]);
  }
  for(j=0;j<NBathc;j++)
  {
    for(i=0;i<8;i++)
    {
      IndexCrTemp[i]=(uint16_t)(16+14+14+14*NBathv+14*j+i);
      IndexAnTemp[i]=(uint16_t)(16+14+14+14*NBathv+14*j+i);
    }
    CreateOperatorCreateAnnihilateList(NFermion, IndexCrTemp, IndexAnTemp, OneList, ZeroList, 8, &OppNc72[j]);
  }

  
  
  /* create the transition operator for pes */
  
  jIm[0]=0;   jIm[1]=0;
  jlist[0]=11; jlist[1]=12; jRe[0]=-sqrt(6.0/7.0); jRe[1]= sqrt(1.0/7.0);
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 2, &TPES[0]);
  jlist[0]=13; jlist[1]=14; jRe[0]=-sqrt(5.0/7.0); jRe[1]= sqrt(2.0/7.0);
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 2, &TPES[1]);
  jlist[0]=15; jlist[1]=16; jRe[0]=-sqrt(4.0/7.0); jRe[1]= sqrt(3.0/7.0);
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 2, &TPES[2]);
  jlist[0]=17; jlist[1]=18; jRe[0]=-sqrt(3.0/7.0); jRe[1]= sqrt(4.0/7.0);
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 2, &TPES[3]);
  jlist[0]=19; jlist[1]=20; jRe[0]=-sqrt(2.0/7.0); jRe[1]= sqrt(5.0/7.0);
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 2, &TPES[4]);
  jlist[0]=21; jlist[1]=22; jRe[0]=-sqrt(1.0/7.0); jRe[1]= sqrt(6.0/7.0);
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 2, &TPES[5]);
  jlist[0]=10;              jRe[0]=1.0;
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 1, &TPES[6]);
  jlist[0]=11; jlist[1]=12; jRe[0]= sqrt(1.0/7.0); jRe[1]= sqrt(6.0/7.0);
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 2, &TPES[7]);
  jlist[0]=13; jlist[1]=14; jRe[0]= sqrt(2.0/7.0); jRe[1]= sqrt(5.0/7.0);
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 2, &TPES[8]);
  jlist[0]=15; jlist[1]=16; jRe[0]= sqrt(3.0/7.0); jRe[1]= sqrt(4.0/7.0);
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 2, &TPES[9]);
  jlist[0]=17; jlist[1]=18; jRe[0]= sqrt(4.0/7.0); jRe[1]= sqrt(3.0/7.0);
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 2, &TPES[10]);
  jlist[0]=19; jlist[1]=20; jRe[0]= sqrt(5.0/7.0); jRe[1]= sqrt(2.0/7.0);
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 2, &TPES[11]);
  jlist[0]=21; jlist[1]=22; jRe[0]= sqrt(6.0/7.0); jRe[1]= sqrt(1.0/7.0);
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 2, &TPES[12]);
  jlist[0]=23;              jRe[0]=1.0;
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 1, &TPES[13]);
  
  /* create the transition operator for ipes */
  
  jIm[0]=0;   jIm[1]=0;
  jlist[0]=11; jlist[1]=12; jRe[0]=-sqrt(6.0/7.0); jRe[1]= sqrt(1.0/7.0);
  CreateOperatorCreateList(NFermion, jlist, jRe, jIm, 2, &TIPES[0]);
  jlist[0]=13; jlist[1]=14; jRe[0]=-sqrt(5.0/7.0); jRe[1]= sqrt(2.0/7.0);
  CreateOperatorCreateList(NFermion, jlist, jRe, jIm, 2, &TIPES[1]);
  jlist[0]=15; jlist[1]=16; jRe[0]=-sqrt(4.0/7.0); jRe[1]= sqrt(3.0/7.0);
  CreateOperatorCreateList(NFermion, jlist, jRe, jIm, 2, &TIPES[2]);
  jlist[0]=17; jlist[1]=18; jRe[0]=-sqrt(3.0/7.0); jRe[1]= sqrt(4.0/7.0);
  CreateOperatorCreateList(NFermion, jlist, jRe, jIm, 2, &TIPES[3]);
  jlist[0]=19; jlist[1]=20; jRe[0]=-sqrt(2.0/7.0); jRe[1]= sqrt(5.0/7.0);
  CreateOperatorCreateList(NFermion, jlist, jRe, jIm, 2, &TIPES[4]);
  jlist[0]=21; jlist[1]=22; jRe[0]=-sqrt(1.0/7.0); jRe[1]= sqrt(6.0/7.0);
  CreateOperatorCreateList(NFermion, jlist, jRe, jIm, 2, &TIPES[5]);
  jlist[0]=10;              jRe[0]=1.0;
  CreateOperatorCreateList(NFermion, jlist, jRe, jIm, 1, &TIPES[6]);
  jlist[0]=11; jlist[1]=12; jRe[0]= sqrt(1.0/7.0); jRe[1]= sqrt(6.0/7.0);
  CreateOperatorCreateList(NFermion, jlist, jRe, jIm, 2, &TIPES[7]);
  jlist[0]=13; jlist[1]=14; jRe[0]= sqrt(2.0/7.0); jRe[1]= sqrt(5.0/7.0);
  CreateOperatorCreateList(NFermion, jlist, jRe, jIm, 2, &TIPES[8]);
  jlist[0]=15; jlist[1]=16; jRe[0]= sqrt(3.0/7.0); jRe[1]= sqrt(4.0/7.0);
  CreateOperatorCreateList(NFermion, jlist, jRe, jIm, 2, &TIPES[9]);
  jlist[0]=17; jlist[1]=18; jRe[0]= sqrt(4.0/7.0); jRe[1]= sqrt(3.0/7.0);
  CreateOperatorCreateList(NFermion, jlist, jRe, jIm, 2, &TIPES[10]);
  jlist[0]=19; jlist[1]=20; jRe[0]= sqrt(5.0/7.0); jRe[1]= sqrt(2.0/7.0);
  CreateOperatorCreateList(NFermion, jlist, jRe, jIm, 2, &TIPES[11]);
  jlist[0]=21; jlist[1]=22; jRe[0]= sqrt(6.0/7.0); jRe[1]= sqrt(1.0/7.0);
  CreateOperatorCreateList(NFermion, jlist, jRe, jIm, 2, &TIPES[12]);
  jlist[0]=23;              jRe[0]=1.0;
  CreateOperatorCreateList(NFermion, jlist, jRe, jIm, 1, &TIPES[13]);
  
  /* create the transition operator for core level pes */
  
  jIm[0]=0;   jIm[1]=0;
  jlist[0]=1; jlist[1]=2; jRe[0]=-sqrt(4.0/5.0); jRe[1]= sqrt(1.0/5.0);
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 2, &TCPES[0]);
  jlist[0]=3; jlist[1]=4; jRe[0]=-sqrt(3.0/5.0); jRe[1]= sqrt(2.0/5.0);
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 2, &TCPES[1]);
  jlist[0]=5; jlist[1]=6; jRe[0]=-sqrt(2.0/5.0); jRe[1]= sqrt(3.0/5.0);
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 2, &TCPES[2]);
  jlist[0]=7; jlist[1]=8; jRe[0]=-sqrt(1.0/5.0); jRe[1]= sqrt(4.0/5.0);
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 2, &TCPES[3]);
  jlist[0]=0;             jRe[0]=1.0;
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 1, &TCPES[4]);
  jlist[0]=1; jlist[1]=2; jRe[0]= sqrt(1.0/5.0); jRe[1]= sqrt(4.0/5.0);
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 2, &TCPES[5]);
  jlist[0]=3; jlist[1]=4; jRe[0]= sqrt(2.0/5.0); jRe[1]= sqrt(3.0/5.0);
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 2, &TCPES[6]);
  jlist[0]=5; jlist[1]=6; jRe[0]= sqrt(3.0/5.0); jRe[1]= sqrt(2.0/5.0);
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 2, &TCPES[7]);
  jlist[0]=7; jlist[1]=8; jRe[0]= sqrt(4.0/5.0); jRe[1]= sqrt(1.0/5.0);
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 2, &TCPES[8]);
  jlist[0]=9;             jRe[0]=1.0;
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 1, &TCPES[9]);
  
  /* create the transition operator for XAS */
  
  CreateOperatorCF(NFermion, 3, 2, Index4Fup, Index4Fdn, Index3Dup, Index3Ddn, Vx, &TXASx);
  CreateOperatorCF(NFermion, 3, 2, Index4Fup, Index4Fdn, Index3Dup, Index3Ddn, Vy, &TXASy);
  CreateOperatorCF(NFermion, 3, 2, Index4Fup, Index4Fdn, Index3Dup, Index3Ddn, Vz, &TXASz);
  CreateOperatorCF(NFermion, 3, 2, Index4Fup, Index4Fdn, Index3Dup, Index3Ddn, Vleft, &TXASl);
  CreateOperatorCF(NFermion, 3, 2, Index4Fup, Index4Fdn, Index3Dup, Index3Ddn, Vright, &TXASr);
  EndTiming("define operators", &Time);
  
  
  printf("\n========================================================");
  printf("\n==== Rotate the operators from lmssz to jjz basis ======");
  printf("\n========================================================\n");
  
  jjzRotationMatrix.QComplex=0;
  jjzRotationMatrix.Ni=10+14+14+NBathc*14+NBathv*14;
  jjzRotationMatrix.Nj=10+14+14+NBathc*14+NBathv*14;
  InitIdentityMatrix(&jjzRotationMatrix);
  for(i=10;i<24;i++)
  {
    jjzRotationMatrix.VRe[i][i]=0;
  }
  jjzRotationMatrix.VRe[10][11]=-sqrt(6.0/7.0); jjzRotationMatrix.VRe[10][12]= sqrt(1.0/7.0);
  jjzRotationMatrix.VRe[11][13]=-sqrt(5.0/7.0); jjzRotationMatrix.VRe[11][14]= sqrt(2.0/7.0);
  jjzRotationMatrix.VRe[12][15]=-sqrt(4.0/7.0); jjzRotationMatrix.VRe[12][16]= sqrt(3.0/7.0);
  jjzRotationMatrix.VRe[13][17]=-sqrt(3.0/7.0); jjzRotationMatrix.VRe[13][18]= sqrt(4.0/7.0);
  jjzRotationMatrix.VRe[14][19]=-sqrt(2.0/7.0); jjzRotationMatrix.VRe[14][20]= sqrt(5.0/7.0);
  jjzRotationMatrix.VRe[15][21]=-sqrt(1.0/7.0); jjzRotationMatrix.VRe[15][22]= sqrt(6.0/7.0);
  
  jjzRotationMatrix.VRe[16][10]=1.0;
  jjzRotationMatrix.VRe[17][11]= sqrt(1.0/7.0); jjzRotationMatrix.VRe[17][12]= sqrt(6.0/7.0);
  jjzRotationMatrix.VRe[18][13]= sqrt(2.0/7.0); jjzRotationMatrix.VRe[18][14]= sqrt(5.0/7.0);
  jjzRotationMatrix.VRe[19][15]= sqrt(3.0/7.0); jjzRotationMatrix.VRe[19][16]= sqrt(4.0/7.0);
  jjzRotationMatrix.VRe[20][17]= sqrt(4.0/7.0); jjzRotationMatrix.VRe[20][18]= sqrt(3.0/7.0);
  jjzRotationMatrix.VRe[21][19]= sqrt(5.0/7.0); jjzRotationMatrix.VRe[21][20]= sqrt(2.0/7.0);
  jjzRotationMatrix.VRe[22][21]= sqrt(6.0/7.0); jjzRotationMatrix.VRe[22][22]= sqrt(1.0/7.0);
  jjzRotationMatrix.VRe[23][23]=1.0;
  
  MatrixTranspose(&jjzRotationMatrix);
  
  
  /* rotate some operators */
  
  RotateOperatorNormalOrder(&Opptemp, OppSx, jjzRotationMatrix, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&OppSx); OperatorCopy(&OppSx, Opptemp); FreeOperator(&Opptemp);
  RotateOperatorNormalOrder(&Opptemp, OppSy, jjzRotationMatrix, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&OppSy); OperatorCopy(&OppSy, Opptemp); FreeOperator(&Opptemp);
  RotateOperatorNormalOrder(&Opptemp, OppSz, jjzRotationMatrix, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&OppSz); OperatorCopy(&OppSz, Opptemp); FreeOperator(&Opptemp);
  RotateOperatorNormalOrder(&Opptemp, OppSsqr, jjzRotationMatrix, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&OppSsqr); OperatorCopy(&OppSsqr, Opptemp); FreeOperator(&Opptemp);
  RotateOperatorNormalOrder(&Opptemp, OppLx, jjzRotationMatrix, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&OppLx); OperatorCopy(&OppLx, Opptemp); FreeOperator(&Opptemp);
  RotateOperatorNormalOrder(&Opptemp, OppLy, jjzRotationMatrix, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&OppLy); OperatorCopy(&OppLy, Opptemp); FreeOperator(&Opptemp);
  RotateOperatorNormalOrder(&Opptemp, OppLz, jjzRotationMatrix, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&OppLz); OperatorCopy(&OppLz, Opptemp); FreeOperator(&Opptemp);
  RotateOperatorNormalOrder(&Opptemp, OppLsqr, jjzRotationMatrix, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&OppLsqr); OperatorCopy(&OppLsqr, Opptemp); FreeOperator(&Opptemp);
  RotateOperatorNormalOrder(&Opptemp, OppJx, jjzRotationMatrix, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&OppJx); OperatorCopy(&OppJx, Opptemp); FreeOperator(&Opptemp);
  RotateOperatorNormalOrder(&Opptemp, OppJy, jjzRotationMatrix, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&OppJy); OperatorCopy(&OppJy, Opptemp); FreeOperator(&Opptemp);
  RotateOperatorNormalOrder(&Opptemp, OppJz, jjzRotationMatrix, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&OppJz); OperatorCopy(&OppJz, Opptemp); FreeOperator(&Opptemp);
  RotateOperatorNormalOrder(&Opptemp, OppJsqr, jjzRotationMatrix, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&OppJsqr); OperatorCopy(&OppJsqr, Opptemp); FreeOperator(&Opptemp);
  RotateOperatorNormalOrder(&Opptemp, Oppldots3d, jjzRotationMatrix, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&Oppldots3d); OperatorCopy(&Oppldots3d, Opptemp); FreeOperator(&Opptemp);
  RotateOperatorNormalOrder(&Opptemp, Oppldots4f, jjzRotationMatrix, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&Oppldots4f); OperatorCopy(&Oppldots4f, Opptemp); FreeOperator(&Opptemp);
  RotateOperatorNormalOrder(&Opptemp, OppUdf, jjzRotationMatrix, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&OppUdf); OperatorCopy(&OppUdf, Opptemp); FreeOperator(&Opptemp);
  RotateOperatorNormalOrder(&Opptemp, OppUff, jjzRotationMatrix, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&OppUff); OperatorCopy(&OppUff, Opptemp); FreeOperator(&Opptemp);
  RotateOperatorNormalOrder(&Opptemp, OppNf, jjzRotationMatrix, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&OppNf); OperatorCopy(&OppNf, Opptemp); FreeOperator(&Opptemp);
  RotateOperatorNormalOrder(&Opptemp, OppNL, jjzRotationMatrix, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&OppNL); OperatorCopy(&OppNL, Opptemp); FreeOperator(&Opptemp);
  RotateOperatorNormalOrder(&Opptemp, OppNfivehalf, jjzRotationMatrix, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&OppNfivehalf); OperatorCopy(&OppNfivehalf, Opptemp); FreeOperator(&Opptemp);
  RotateOperatorNormalOrder(&Opptemp, OppNsevenhalf, jjzRotationMatrix, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&OppNsevenhalf); OperatorCopy(&OppNsevenhalf, Opptemp); FreeOperator(&Opptemp);
  RotateOperatorNormalOrder(&Opptemp, OppDelta, jjzRotationMatrix, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&OppDelta); OperatorCopy(&OppDelta, Opptemp); FreeOperator(&Opptemp);
  
  RotateOperatorNormalOrder(&Opptemp, OppHpart, jjzRotationMatrix, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&OppHpart); OperatorCopy(&OppHpart, Opptemp); FreeOperator(&Opptemp);
  
  
  
  
  /* rotate hopping operators */
    /* nothing to do */
  
  /* rotate pes */
  
  for(i=0;i<14;i++)
  {
    RotateOperatorNormalOrder(&Opptemp, TPES[i], jjzRotationMatrix, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&TPES[i]); OperatorCopy(&TPES[i], Opptemp); FreeOperator(&Opptemp);
  }
  
  /* rotate ipes */
  
  for(i=0;i<14;i++)
  {
    RotateOperatorNormalOrder(&Opptemp, TIPES[i], jjzRotationMatrix, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&TIPES[i]); OperatorCopy(&TIPES[i], Opptemp); FreeOperator(&Opptemp);
  }
  
  /* rotate core level PES */
  
    /* nothing to do */
  
  /* rotate TXAS */
  
  RotateOperatorNormalOrder(&Opptemp, TXASx, jjzRotationMatrix, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&TXASx); OperatorCopy(&TXASx, Opptemp); FreeOperator(&Opptemp);
  RotateOperatorNormalOrder(&Opptemp, TXASy, jjzRotationMatrix, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&TXASy); OperatorCopy(&TXASy, Opptemp); FreeOperator(&Opptemp);
  RotateOperatorNormalOrder(&Opptemp, TXASz, jjzRotationMatrix, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&TXASz); OperatorCopy(&TXASz, Opptemp); FreeOperator(&Opptemp);
  RotateOperatorNormalOrder(&Opptemp, TXASl, jjzRotationMatrix, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&TXASl); OperatorCopy(&TXASl, Opptemp); FreeOperator(&Opptemp);
  RotateOperatorNormalOrder(&Opptemp, TXASr, jjzRotationMatrix, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&TXASr); OperatorCopy(&TXASr, Opptemp); FreeOperator(&Opptemp);
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  printf("\n========================================================");
  printf("\n========== occupation dependent Hamiltonian  ===========");
  printf("\n========================================================\n");
  
 

  occupation=NfStart/14.0;
  occupation52=Nf52Start/6.0;
  occupation72=Nf72Start/8.0;
  
  
  CreateWaveFunctionFromString(&psi, CreateString, 10+14*(1+NBathv), NFermion, 0, NBitsKey);
  
  
  for(iGroundStateLoops=0;iGroundStateLoops<NGroundStateLoops;iGroundStateLoops++)
  {
    printf("\n========================================================");
    printf("\n==========       Start of loop %3i / %3i     ===========",iGroundStateLoops,NGroundStateLoops);
    printf("\n========================================================\n");
    
    
    onsiteshift=Uff/2.0+Delta;
    onsiteshift52=onsiteshift;
    onsiteshift72=onsiteshift;
    
    
    printf("guess the rotation of the bath by assuming an occupation of (total) %lf/14.\n",occupation*14.0);
    DetermineImpurityShiftFromOccupationAndersonMatrix(GBath, occupation, &onsiteshift);
    printf("guess the rotation of the bath by assuming an occupation of (5/2)   %lf/6.\n",occupation52*6.0);
    DetermineImpurityShiftFromOccupationAndersonMatrix(GBath, occupation52, &onsiteshift52);
    printf("guess the rotation of the bath by assuming an occupation of (7/2)   %lf/8.\n",occupation72*8.0);
    DetermineImpurityShiftFromOccupationAndersonMatrix(GBath, occupation72, &onsiteshift72);
    
    printf("Onsite shift is %22.15E\n",onsiteshift);
    printf("Onsite shift for 5/2 is %22.15E\n",onsiteshift52);
    printf("Onsite shift for 7/2 is %22.15E\n",onsiteshift72);
    
    kbTeff=0.0001;
    AndersonMatrixToBADoubleTriDiagonalMatrixTemperature(GBath, onsiteshift, kbTeff, &GBathBA, &occupation);
    AndersonMatrixToBADoubleTriDiagonalMatrixTemperature(GBath, onsiteshift52, kbTeff, &GBathBA52, &occupation52);
    AndersonMatrixToBADoubleTriDiagonalMatrixTemperature(GBath, onsiteshift72, kbTeff, &GBathBA72, &occupation72);
    
    /* create OppHyb based on GBathBA */
    
    OperatorCopy(&OppHyb, OppHybIL);
    OperatorScale(&OppHyb, GBathBA.acvBA, 0.0);
    OperatorAdd(&OppHyb, OppNL, GBathBA.acvB, 0.0);
    if((NBathv!=0)&&(GBathBA.val.N>1))
    {
      OperatorAdd(&OppHyb, OppHybIv, GBathBA.bvA0, 0.0);
      OperatorAdd(&OppHyb, OppHybLv, GBathBA.bvB0, 0.0);
      OperatorAdd(&OppHyb, OppNv[0], GBathBA.val.a[2], 0.0);
    }
    if((NBathc!=0)&&(GBathBA.con.N>1))
    {
      OperatorAdd(&OppHyb, OppHybIc, GBathBA.bcA0, 0.0);
      OperatorAdd(&OppHyb, OppHybLc, GBathBA.bcB0, 0.0);
      OperatorAdd(&OppHyb, OppNc[0], GBathBA.con.a[2], 0.0);
    }
    for(j=1;j<Min(NBathv,GBathBA.val.N-1);j++)
    {
      OperatorAdd(&OppHyb, OppHybv[j-1], GBathBA.val.b[j+1], 0.0);
      OperatorAdd(&OppHyb, OppNv[j],     GBathBA.val.a[j+2], 0.0);
    }
    for(j=1;j<Min(NBathc,GBathBA.con.N-1);j++)
    {
      OperatorAdd(&OppHyb, OppHybc[j-1], GBathBA.con.b[j+1], 0.0);
      OperatorAdd(&OppHyb, OppNc[j],     GBathBA.con.a[j+2], 0.0);
    }
    
    /* create OppHyb52 based on GBathBA52 */
    
    OperatorCopy(&OppHyb52, OppHybIL52);
    OperatorScale(&OppHyb52, GBathBA52.acvBA, 0.0);
    OperatorAdd(&OppHyb52, OppNL52, GBathBA52.acvB, 0.0);
    if((NBathv!=0)&&(GBathBA52.val.N>1))
    {
      OperatorAdd(&OppHyb52, OppHybIv52, GBathBA52.bvA0, 0.0);
      OperatorAdd(&OppHyb52, OppHybLv52, GBathBA52.bvB0, 0.0);
      OperatorAdd(&OppHyb52, OppNv52[0], GBathBA52.val.a[2], 0.0);
    }
    if((NBathc!=0)&&(GBathBA52.con.N>1))
    {
      OperatorAdd(&OppHyb52, OppHybIc52, GBathBA52.bcA0, 0.0);
      OperatorAdd(&OppHyb52, OppHybLc52, GBathBA52.bcB0, 0.0);
      OperatorAdd(&OppHyb52, OppNc52[0], GBathBA52.con.a[2], 0.0);
    }
    for(j=1;j<Min(NBathv,GBathBA52.val.N-1);j++)
    {
      OperatorAdd(&OppHyb52, OppHybv52[j-1], GBathBA52.val.b[j+1], 0.0);
      OperatorAdd(&OppHyb52, OppNv52[j],     GBathBA52.val.a[j+2], 0.0);
    }
    for(j=1;j<Min(NBathc,GBathBA52.con.N-1);j++)
    {
      OperatorAdd(&OppHyb52, OppHybc52[j-1], GBathBA.con.b[j+1], 0.0);
      OperatorAdd(&OppHyb52, OppNc52[j],     GBathBA.con.a[j+2], 0.0);
    }
    
    /* create OppHyb72 based on GBathBA */
    
    OperatorCopy(&OppHyb72, OppHybIL72);
    OperatorScale(&OppHyb72, GBathBA72.acvBA, 0.0);
    OperatorAdd(&OppHyb72, OppNL72, GBathBA72.acvB, 0.0);
    if((NBathv!=0)&&(GBathBA72.val.N>1))
    {
      OperatorAdd(&OppHyb72, OppHybIv72, GBathBA72.bvA0, 0.0);
      OperatorAdd(&OppHyb72, OppHybLv72, GBathBA72.bvB0, 0.0);
      OperatorAdd(&OppHyb72, OppNv72[0], GBathBA72.val.a[2], 0.0);
    }
    if((NBathc!=0)&&(GBathBA72.con.N>1))
    {
      OperatorAdd(&OppHyb72, OppHybIc72, GBathBA72.bcA0, 0.0);
      OperatorAdd(&OppHyb72, OppHybLc72, GBathBA72.bcB0, 0.0);
      OperatorAdd(&OppHyb72, OppNc72[0], GBathBA72.con.a[2], 0.0);
    }
    for(j=1;j<Min(NBathv,GBathBA72.val.N-1);j++)
    {
      OperatorAdd(&OppHyb72, OppHybv72[j-1], GBathBA72.val.b[j+1], 0.0);
      OperatorAdd(&OppHyb72, OppNv72[j],     GBathBA72.val.a[j+2], 0.0);
    }
    for(j=1;j<Min(NBathc,GBathBA72.con.N-1);j++)
    {
      OperatorAdd(&OppHyb72, OppHybc72[j-1], GBathBA72.con.b[j+1], 0.0);
      OperatorAdd(&OppHyb72, OppNc72[j],     GBathBA72.con.a[j+2], 0.0);
    }
    
    /* add OppHyb to Hamiltonian */
    
    OperatorCopy(&OppH, OppHpart);
    OperatorAdd(&OppH, OppHyb52, 1.0, 0.0);
    OperatorAdd(&OppH, OppHyb72, 1.0, 0.0);
    
    printf("\n========================================================");
    printf("\n==========   find groundstate %3i-st time    ===========",iGroundStateLoops);
    printf("\n========================================================\n");
    BeginTiming("find groundstate", &Time);
    
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
    
    
    if(iGroundStateLoops!=NGroundStateLoops-1)
    {
      FreeOperator(&OppHyb);
      FreeOperator(&OppHyb52);
      FreeOperator(&OppHyb72);
      FreeOperator(&OppH);
      FreeBADoubleTriDiagonalMatrix(&GBathBA);
      FreeBADoubleTriDiagonalMatrix(&GBathBA52);
      FreeBADoubleTriDiagonalMatrix(&GBathBA72);
      
      occupation=Nf/14.0;
      occupation52=Nfivehalf/6.0;
      occupation72=Nsevenhalf/8.0;
      
      
      
    }
    
    
  }
  
  
  
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
  
  if(doCPESShort)
  {
    printf("\n========================================================");
    printf("\n==========            CPES Short             ===========");
    printf("\n========================================================\n");
    BeginTiming("CPESShort", &Time);
    
    for(i=28;i<30;i++)
    {
      Spectra[i].N=7000;
      Spectra[i].Emin=-40.0;
      Spectra[i].Emax= 30.0;
      Spectra[i].dE = (Spectra[i].Emax - Spectra[i].Emin) / Spectra[i].N;
      Spectra[i].Gamma = 0.08;
      SpectrumAllocate(Spectra[i], &Spectra[i]);
    }
    
    for(i=0;i<2;i++)
    {
      BeginTiming("spectrum", &Time);
      OperatorPsi(TCPES[4*i], psi, &Opsi);
      LanczosTriDiagonalizeRestricted(OppH, Opsi, Restrictions, epsilon, &SpecTri,  Spectra[28+i].Gamma);
      TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[28+i]);
      EndTiming("spectrum", &Time);
    }
    
    FilePrintSpectra(&Spectra[28], 2, "CeLigandField_CPESShort","w");
    EndTiming("CPESShort", &Time);
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