#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <unistd.h>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#define omp_get_num_procs() 0
#define omp_get_max_threads() 0
#endif

#include "../Core/StructDef.c"
#include "../header.h"
#include "../Core/TimeStamp.c"
#include "../Core/Hash.c"
#include "../Core/BasicMath.c"
#include "../IO/Print.c"
#include "../IO/Read.c"
#include "../IO/BitMap.c"
#include "../Core/Error.c"
#include "../Core/WaveFunctionArithmetic.c"
#include "../Core/OperatorArithmetic.c"
#include "../Core/OperatorWaveFunction.c"
#include "../Core/Lanczos.c"
#include "../Core/Spectral.c"
#include "../Core/DensityMatrix.c"
#include "../Core/Random.c"
#include "../Core/MatrixOperator.c"
#include "../Model/StandardOperators.c"
#include "../Hubbard_Model/CellFunctions.c"
#include "../Tight_Binding/Tight_Binding_lib.c"
#include "../Tight_Binding/ImpurityUnfolding.c"
#include "../DMFT/SingleSiteFunctions.c"

int mainNoDyson();
int oldmainNoDysonSecond();
int mainDyson();
int oldmainDysonSecond();

int main(int argc, const char * argv[])
{
  TimeType Time;
  DMFTSettings HSet;
  InitializeTiming(&Time);
  BeginTiming("Total",&Time);
  TimeStampAll();
  argc-=1;
  /* read the DMFT settings */
  if(FileReadDMFTSettings(&HSet,argv[argc])){printf("Fatal Error in main\n Failed to read DMFT Settings from file %s\n",argv[argc]);return 1;};
  PrintDMFTSettings(HSet);
  
  WaveFunctionType Psi,OPsi,PsiTemp;
  TriDiagonalMatrixType G0,GBath,GclPup,GclPdn,GclMup,GclMdn,Gclup,Gcldn,Gcl,Sigma,GBathNew,GBathOld,gPerturbation[4];
  SecondDoubleTriDiagonalMatrixType GAnd,GAndTrunk;
  OperatorType HAnd, OPESdn,OPESup,OIPESdn,OIPESup,VPerturbation[4];
  MatrixType T;
  SpectraType Spectra[10];
  char GBathFileName[255],PsiFileName[255],SpecFileName[25],Buffer[255];
  uint32_t iOuter,i,PsiNBathV,PsiNBathC;
  uint16_t oadd[1];
  double Ef=0.0,epsilonAnd=0.00001,mu;
  double PsiEnergy,PsiQuality,PsiHist[20],*PsiRho;
  /* read G0 */
  if(FileReadPlainTriDiagonalMatrixType(&G0,HSet.G0FileName)){printf("Fatal Error in main\n Failed to read G0 from file %s\n",HSet.G0FileName);return 1;};
  strcpy(G0.Name,"G0");
  /* read G_Bath_iStart and psi_iStart*/
  if(HSet.iStart==0)
  {
    sprintf(GBathFileName,"%s000",HSet.GRootFileName);
    sprintf(PsiFileName,"%s000",HSet.WaveFunctionRootFileName);
  }
  else if(HSet.iStart<10)
  {
    sprintf(GBathFileName,"%s00%1i",HSet.GRootFileName,HSet.iStart);
    sprintf(PsiFileName,"%s00%1i",HSet.WaveFunctionRootFileName,HSet.iStart);
  }
  else if(HSet.iStart<100)
  {
    sprintf(GBathFileName,"%s0%2i",HSet.GRootFileName,HSet.iStart);
    sprintf(PsiFileName,"%s0%2i",HSet.WaveFunctionRootFileName,HSet.iStart);
  }
  else
  {
    sprintf(GBathFileName,"%s%i",HSet.GRootFileName,HSet.iStart);
    sprintf(PsiFileName,"%s%i",HSet.WaveFunctionRootFileName,HSet.iStart);
  }
  if(FileReadPlainTriDiagonalMatrixType(&GBath,GBathFileName)){printf("Fatal Error in main\n Failed to read GBath from file %s !\n",GBathFileName);return 1;};
  if(file_exists(PsiFileName))
  {
    Psi.NBitsKey=HSet.NBitsKey;
    if(FileReadWaveFunctionPlainBinary(&Psi, PsiFileName)){printf("Fatal Error in main\n Failed to read Psi from file %s !\n",PsiFileName);return 1;};
    printf("\nRead the starting wave function from the file %s\n",PsiFileName);
    sscanf(Psi.Name,"%s %i %i",Buffer,&PsiNBathC,&PsiNBathV);
    printf("Wave-function \"%s\" found has %i conduction bath sites, %i valence bath sites and %i Fermions in basis\n",Psi.Name,PsiNBathC,PsiNBathV,Psi.NFermion);
    if((2*(PsiNBathC+PsiNBathV+2))!=Psi.NFermion)
    {
      printf("The total number of Fermions is incompatible witht the bath sizes\n");
      printf("This wave-function can not be used to start from\n");
      return 1;
    }
  }
  else
  {
    printf("\nWave-function file: '%s' does not excist\n Started from new function\n",PsiFileName);
    CreateWaveFunctionDMFTSecond(&Psi, HSet.NBath, HSet.NBath, HSet.NBitsKey);
  }
  strcpy(GBath.Name,"Gbath");
  /* do the DMFT loops */
  for(iOuter=HSet.iStart;iOuter<HSet.iStart+HSet.Ni;iOuter++)
  {
    BeginTiming("DMFTLoop",&Time);
    printf("\n\n***********************************************************************\n");
    printf("Start of loop %i from %i to %i\n",iOuter+1,HSet.iStart,HSet.iStart+HSet.Ni);
    printf("***********************************************************************\n\n\n");
    /* Create splitted anderson matrix */
    TriDiagonalMatrixToSecondDoubleTriDiagonalMatrix(GBath, &T, &GAnd, Ef, epsilonAnd);
    /* Create GAnd and GAndTrunk as a copy pointing to the same memory */
    GAndTrunk.Ef=GAnd.Ef;
    GAndTrunk.Epsilon=GAnd.Epsilon;
    GAndTrunk.a0=GAnd.a0;
    GAndTrunk.b0=GAnd.b0;
    GAndTrunk.Dt.Ef=GAnd.Dt.Ef;
    GAndTrunk.Dt.Epsilon=GAnd.Dt.Epsilon;
    GAndTrunk.Dt.a0=GAnd.Dt.a0;
    GAndTrunk.Dt.b0=GAnd.Dt.b0;
    GAndTrunk.Dt.val.N=GAnd.Dt.val.N;
    GAndTrunk.Dt.val.a=GAnd.Dt.val.a;
    GAndTrunk.Dt.val.b=GAnd.Dt.val.b;
    GAndTrunk.Dt.con.N=GAnd.Dt.con.N;
    GAndTrunk.Dt.con.a=GAnd.Dt.con.a;
    GAndTrunk.Dt.con.b=GAnd.Dt.con.b;
    /* Create VPerturbation */
    for(i=0;i<4;i++)
    {
      InitOperator(&(VPerturbation[i]));
      VPerturbation[i].NFermion=4+2*(Min(GAnd.Dt.val.N,HSet.NBath)+Min(GAnd.Dt.con.N,HSet.NBath));
      VPerturbation[i].NBoson=0;
      VPerturbation[i].NWord=Min(VPerturbation[i].NFermion,1+(VPerturbation[i].NFermion-1)/8);
      VPerturbation[i].QComplex=0;
      InitOperatorFixedLength(&(VPerturbation[i]),1,1);
    }
    if(GAnd.Dt.val.N>HSet.NBath)
    {
      oadd[0]=CrAnMask+3+HSet.NBath;
      RealOperatorAddLader(&(VPerturbation[0]), oadd, 1, GAnd.Dt.val.b[HSet.NBath]);
      oadd[0]=CrAnMask+3+Min(GAnd.Dt.con.N,HSet.NBath)+2*HSet.NBath;
      RealOperatorAddLader(&(VPerturbation[1]), oadd, 1, GAnd.Dt.val.b[HSet.NBath]);
    }
    if(GAnd.Dt.con.N>HSet.NBath)
    {
      oadd[0]=3+Min(GAnd.Dt.val.N,HSet.NBath)+HSet.NBath;
      RealOperatorAddLader(&(VPerturbation[2]), oadd, 1, GAnd.Dt.con.b[HSet.NBath]);
      oadd[0]=3+2*Min(GAnd.Dt.val.N,HSet.NBath)+2*HSet.NBath;
      RealOperatorAddLader(&(VPerturbation[3]), oadd, 1, GAnd.Dt.con.b[HSet.NBath]);
    }
    /* Create gPerturbation */
    if(GAnd.Dt.val.N>HSet.NBath)
    {
      gPerturbation[0].N=GAnd.Dt.val.N-HSet.NBath;
      gPerturbation[1].N=GAnd.Dt.val.N-HSet.NBath;
      InitTriDiagonalMatrix(&(gPerturbation[0]));
      InitTriDiagonalMatrix(&(gPerturbation[1]));
      for(i=1;i<GAnd.Dt.val.N-HSet.NBath;i++)
      {
        gPerturbation[0].a[i]=-GAnd.Dt.val.a[i+HSet.NBath];
        gPerturbation[0].b[i]=GAnd.Dt.val.b[i+HSet.NBath];
        gPerturbation[1].a[i]=-GAnd.Dt.val.a[i+HSet.NBath];
        gPerturbation[1].b[i]=GAnd.Dt.val.b[i+HSet.NBath];
      }
      gPerturbation[0].a[GAnd.Dt.val.N-HSet.NBath]=-GAnd.Dt.val.a[GAnd.Dt.val.N];
      gPerturbation[1].a[GAnd.Dt.val.N-HSet.NBath]=-GAnd.Dt.val.a[GAnd.Dt.val.N];
    }
    else
    {
      gPerturbation[0].N=1;
      gPerturbation[1].N=1;
      InitTriDiagonalMatrix(&(gPerturbation[0]));
      InitTriDiagonalMatrix(&(gPerturbation[1]));
      gPerturbation[0].a[1]=0;
      gPerturbation[1].a[1]=0;
    }
    if(GAnd.Dt.con.N>HSet.NBath)
    {
      gPerturbation[2].N=GAnd.Dt.con.N-HSet.NBath;
      gPerturbation[3].N=GAnd.Dt.con.N-HSet.NBath;
      InitTriDiagonalMatrix(&(gPerturbation[2]));
      InitTriDiagonalMatrix(&(gPerturbation[3]));
      for(i=1;i<GAnd.Dt.con.N-HSet.NBath;i++)
      {
        gPerturbation[2].a[i]=GAnd.Dt.con.a[i+HSet.NBath];
        gPerturbation[2].b[i]=GAnd.Dt.con.b[i+HSet.NBath];
        gPerturbation[3].a[i]=GAnd.Dt.con.a[i+HSet.NBath];
        gPerturbation[3].b[i]=GAnd.Dt.con.b[i+HSet.NBath];
      }
      gPerturbation[2].a[GAnd.Dt.con.N-HSet.NBath]=GAnd.Dt.con.a[GAnd.Dt.con.N];
      gPerturbation[3].a[GAnd.Dt.con.N-HSet.NBath]=GAnd.Dt.con.a[GAnd.Dt.con.N];
    }
    else
    {
      gPerturbation[2].N=1;
      gPerturbation[3].N=1;
      InitTriDiagonalMatrix(&(gPerturbation[2]));
      InitTriDiagonalMatrix(&(gPerturbation[3]));
      gPerturbation[2].a[1]=0;
      gPerturbation[3].a[1]=0;
    };
    /* Create HAnd */
    FreeMatrix(&T);
    if(GAnd.Dt.con.N<HSet.NBath)
    {
      printf("Warning: GAnd.Dt.con.N(%i)<HSet.NBath(%i)\n",GAnd.Dt.con.N,HSet.NBath);
    }
    else
    {
      GAndTrunk.Dt.con.N=HSet.NBath;
    }
    if(GAnd.Dt.val.N<HSet.NBath)
    {
      printf("Warning: GAnd.Dt.val.N(%i)<HSet.NBath(%i)\n",GAnd.Dt.val.N,HSet.NBath);
    }
    else
    {
      GAndTrunk.Dt.val.N=HSet.NBath;
    }
    mu=-HSet.U/2;
    SecondDoubleTriDiagonalMatrixToAndersonImpurityOperator(&HAnd, GAndTrunk, GAndTrunk, HSet.U, mu);
    /* test if psi has correct dimensions and correct if needed */
    sscanf(Psi.Name,"%s %i %i",Buffer,&PsiNBathC,&PsiNBathV);
    if(PsiNBathC!=GAndTrunk.Dt.con.N)
    {
      printf("\nBefore diagonalizing:\nThe wavefunction has %i conduction bath sites\nThe Anderson Hamiltonian %i\n",PsiNBathC,GAndTrunk.Dt.con.N);
      if(PsiNBathC>GAndTrunk.Dt.con.N)
      {
        printf("Removing %i unoccupied fermions from the wave-function\n",PsiNBathC-GAndTrunk.Dt.con.N);
        for(i=0;i<PsiNBathC-GAndTrunk.Dt.con.N;i++)
        {
          WaveFunctionRemoveFermion(Psi, &PsiTemp, Psi.NFermion-1);
          WaveFunctionFree(&Psi);
          WaveFunctionRemoveFermion(PsiTemp, &Psi, 4+PsiNBathV+PsiNBathC-i-1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      else
      {
        printf("Adding %i unoccupied fermions to the wave-function\n",GAndTrunk.Dt.con.N-PsiNBathC);
        for(i=0;i<GAndTrunk.Dt.con.N-PsiNBathC;i++)
        {
          WaveFunctionAddFermion(Psi, &PsiTemp, Psi.NFermion,0);
          WaveFunctionFree(&Psi);
          WaveFunctionAddFermion(PsiTemp, &Psi, 4+PsiNBathV+PsiNBathC+i,0);
          WaveFunctionFree(&PsiTemp);
        }
      }
    }
    if(PsiNBathV!=GAndTrunk.Dt.val.N)
    {
      printf("\nBefore diagonalizing:\nThe wavefunction has %i valence bath sites\nThe Anderson Hamiltonian %i\n",PsiNBathV,GAndTrunk.Dt.val.N);
      if(PsiNBathV>GAndTrunk.Dt.val.N)
      {
        printf("Removing %i occupied fermions from the wave-function\n",PsiNBathV-GAndTrunk.Dt.val.N);
        for(i=0;i<PsiNBathV-GAndTrunk.Dt.val.N;i++)
        {
          WaveFunctionRemoveFermion(Psi, &PsiTemp, Psi.NFermion-GAndTrunk.Dt.con.N-1);
          WaveFunctionFree(&Psi);
          WaveFunctionRemoveFermion(PsiTemp, &Psi, 4+PsiNBathV-i-1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      else
      {
        printf("Adding %i occupied fermions to the wave-function\n",GAndTrunk.Dt.val.N-PsiNBathV);
        for(i=0;i<GAndTrunk.Dt.val.N-PsiNBathV;i++)
        {
          WaveFunctionAddFermion(Psi, &PsiTemp, Psi.NFermion-GAndTrunk.Dt.val.N,1);
          WaveFunctionFree(&Psi);
          WaveFunctionAddFermion(PsiTemp, &Psi, 4+PsiNBathV+i,1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      printf("Corrected basis\n");
      PrintWaveFunctionReadableTruncated(Psi, 0.1);
    }
    sprintf(Psi.Name,"DMFT %i %i",GAndTrunk.Dt.con.N,GAndTrunk.Dt.val.N);
    /* create ground-state wave-function */
    BeginTiming("Lanczos GroundState",&Time);
    LanczosGroundState(HAnd,&Psi,&HSet.zeroLancWavefunction,HSet.zeroLancEnergy,HSet.U,&Psi,0);
    LanczosGroundState(HAnd,&Psi,&HSet.zeroLancWavefunction,HSet.zeroLancEnergy,HSet.U,&Psi,0);
    EndTiming("Lanczos GroundState",&Time);
    PsiHermitianOperatorPsi(HAnd, Psi, &PsiEnergy);
    PsiHermitianOperatorSqrPsi(HAnd,Psi,&PsiQuality);
    PsiQuality/=Sqr(PsiEnergy);
    PrintWaveFunctionReadableTruncated(Psi, 0.1);
    CalculatePrefactorHistogram(PsiHist, Psi, 10.0, 20);
    printf("\nHistogram of prefactors in Wave-function\n");
    for(i=0;i<19;i++)
    {
      printf("10^(-%2i) to 10^(-%2i) | %22.15E\n",i,i+1,PsiHist[i]);
    }
    printf("10^(-%2i) to 0        | %22.15E\n",i,PsiHist[19]);
    printf("\n");
    printf("Energy                  <psi|H^2|psi>/(<psi|H|psi>^2)\n");
    printf("%22.15E  %22.15E\n",PsiEnergy,PsiQuality);
    PsiRho=calloc(Psi.NFermion,sizeof(double));
    if(PsiRho==NULL){printf("calloc failed in main for PsiRho\n Fatal Error\n");fflush(stdout);return 1;}
    CalculateDensity(PsiRho,Psi);
    printf("Occupation \n");
    for(i=0;i<Psi.NFermion;i++)
    {
      printf("%5i %22.15E\n",i,PsiRho[i]);
    }
    free(PsiRho);
    /*                                           */
    /* create spin up and spin down PES and IPES */
    /*                                           */
    CreatePESdnOperatorDMFT(&OPESdn, HAnd.NFermion);
    CreatePESupOperatorDMFT(&OPESup, HAnd.NFermion);
    CreateIPESdnOperatorDMFT(&OIPESdn, HAnd.NFermion);
    CreateIPESupOperatorDMFT(&OIPESup, HAnd.NFermion);
    GclMup.N=Max(HSet.NTriDiag,G0.N);
    InitTriDiagonalMatrix(&GclMup);
    GclMdn.N=Max(HSet.NTriDiag,G0.N);
    InitTriDiagonalMatrix(&GclMdn);
    GclPup.N=Max(HSet.NTriDiag,G0.N);
    InitTriDiagonalMatrix(&GclPup);
    GclPdn.N=Max(HSet.NTriDiag,G0.N);
    InitTriDiagonalMatrix(&GclPdn);
    BeginTiming("Lanczos Tri Diagonalization Dyson",&Time);
    OperatorPsi(OPESup,Psi,&OPsi);
    LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation, gPerturbation, 4, &GclMup);
    EndTiming("Lanczos Tri Diagonalization Dyson",&Time);
    BeginTiming("Lanczos Tri Diagonalization Dyson",&Time);
    OperatorPsi(OPESdn,Psi,&OPsi);
    LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation, gPerturbation, 4, &GclMdn);
    EndTiming("Lanczos Tri Diagonalization Dyson",&Time);
    BeginTiming("Lanczos Tri Diagonalization Dyson",&Time);
    OperatorPsi(OIPESup,Psi,&OPsi);
    LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation, gPerturbation, 4, &GclPup);
    EndTiming("Lanczos Tri Diagonalization Dyson",&Time);
    BeginTiming("Lanczos Tri Diagonalization Dyson",&Time);
    OperatorPsi(OIPESdn,Psi,&OPsi);
    LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation, gPerturbation, 4, &GclPdn);
    EndTiming("Lanczos Tri Diagonalization Dyson",&Time);
    /* G math - create Gcl as the sum of PES - IPES */
    printf("=== Start G math ===\n");
    BeginTiming("G math",&Time);
    ShiftTriDiagonalMatrix(&GclMdn, -PsiEnergy);
    ShiftTriDiagonalMatrix(&GclMup, -PsiEnergy);
    ShiftTriDiagonalMatrix(&GclPdn, -PsiEnergy);
    ShiftTriDiagonalMatrix(&GclPup, -PsiEnergy);
    for(i=1;i<=GclMdn.N;i++)
    {
      GclMdn.a[i]*=-1;
    }
    for(i=1;i<=GclMup.N;i++)
    {
      GclMup.a[i]*=-1;
    }
    strcpy(GclPup.Name,"Gcluster-Plus-Up");
    strcpy(GclPdn.Name,"Gcluster-Plus-Dn");
    strcpy(GclMup.Name,"Gcluster-Min-Up");
    strcpy(GclMdn.Name,"Gcluster-Min-Dn");
    printf("Start SumGreensfunctionsTriDiagonalRepresentation\n");fflush(stdout);
    SumGreensfunctionsTriDiagonalRepresentation(GclMup, GclPup, &Gclup);
    strcpy(Gclup.Name,"Gcluster-Up");
    SumGreensfunctionsTriDiagonalRepresentation(GclMdn, GclPdn, &Gcldn);
    strcpy(Gcldn.Name,"Gcluster-Dn");
    SumGreensfunctionsTriDiagonalRepresentation(Gcldn, Gclup, &Gcl);
    strcpy(Gcl.Name,"Gcluster");
    Gcl.b[0]=1.0;
    ReducePolesTriDiagonalMatrix(&Gcl,Max(4*HSet.NTriDiag,6*G0.N),0.0);
    /* symmetrize */
    SymmetrizeTridiagonalMatrix(&Gcl);
    /* G math - create Self energy From Gcl and GBath */
    BeginTiming("G math Selfenergy",&Time);
    printf("Start SelfEnergyTriDiagonalRepresentation\n");
    SelfEnergyTriDiagonalRepresentation(GBath, Gcl, &Sigma);
    strcpy(Sigma.Name,"Sigma");
    EndTiming("G math Selfenergy",&Time);
    
    
    printf("Start GBathTriDiagonalRepresentation\n");
    GBathTriDiagonalRepresentation(G0, Sigma, &GBathNew,Max(4*HSet.NTriDiag,8*G0.N));
    strcpy(GBathNew.Name,"Gbath-New");
    
    /* calculate GBath for the next loop by mixing GBath with GBathNew */
    printf("  Start mixing GBath with GBathNew\n");
    TriDiagonalMatrixCopy(&GBathOld, GBath);
    strcpy(GBathOld.Name,"Gbath-Old");
    FreeTriDiagonalMatrix(&GBath);
    GBathOld.b[0]=sqrt(HSet.Mix);
    GBathNew.b[0]=sqrt(1-HSet.Mix);
    SumGreensfunctionsTriDiagonalRepresentation(GBathNew, GBathOld, &GBath);
    strcpy(GBath.Name,"Gbath");
    GBathOld.b[0]=1;
    GBathNew.b[0]=1;
    ReducePolesTriDiagonalMatrix(&GBath,Max(4*HSet.NTriDiag,2*G0.N),0.0);
    /* symmetrize */
    SymmetrizeTridiagonalMatrix(&GBath);
    EndTiming("G math",&Time);
    /* create spectra */
    printf("Start create spectra\n");
    Spectra[0].N       =HSet.NE;
    Spectra[0].Emin    =HSet.Emin;
    Spectra[0].Emax    =HSet.Emax;
    Spectra[0].dE      =(HSet.Emax-HSet.Emin)/HSet.NE;
    Spectra[0].EminPole=HSet.Emin;
    Spectra[0].EmaxPole=HSet.Emax;
    Spectra[0].Gamma   =HSet.Gamma;
    for(i=0;i<10;i++)
    {
      SpectrumAllocate(Spectra[0], &Spectra[i]);
    }
    TriDiagonalMatrixToSpectrum(0,GBath,&(Spectra[0]));              /*  3 */
    TriDiagonalMatrixToSpectrum(0,GBathOld,&(Spectra[1]));           /*  5 */
    TriDiagonalMatrixToSpectrum(0,GBathNew,&(Spectra[2]));           /*  7 */
    SecondDoubleTriDiagonalMatrixToSpectrum(0,GAnd,&(Spectra[3]));   /*  9 */
    TriDiagonalMatrixToSpectrum(0,G0,&(Spectra[4]));                 /* 11 */
    TriDiagonalMatrixToSpectrum(0,Sigma,&(Spectra[5]));              /* 13 */
//    DoubleTriDiagonalMatrixToSpectrum(0,SigmaDouble,&(Spectra[6]));  /* 15 */
    TriDiagonalMatrixToSpectrum(0,GclPup,&(Spectra[7]));             /* 17 */
    TriDiagonalMatrixToSpectrum(0,GclMup,&(Spectra[8]));             /* 19 */
    TriDiagonalMatrixToSpectrum(0,Gcl,&(Spectra[9]));                /* 21 */
    /* write to disk */
    printf("Start write to dist\n");
    if(iOuter==0)
    {
      sprintf(GBathFileName,"%s001",HSet.GRootFileName);
      sprintf(PsiFileName,"%s001",HSet.WaveFunctionRootFileName);
      sprintf(SpecFileName,"%s001",HSet.SpecRootFileName);
    }
    else if(iOuter<9)
    {
      sprintf(GBathFileName,"%s00%1i",HSet.GRootFileName,iOuter+1);
      sprintf(PsiFileName,"%s00%1i",HSet.WaveFunctionRootFileName,iOuter+1);
      sprintf(SpecFileName,"%s00%1i",HSet.SpecRootFileName,iOuter+1);
    }
    else if(iOuter<99)
    {
      sprintf(GBathFileName,"%s0%2i",HSet.GRootFileName,iOuter+1);
      sprintf(PsiFileName,"%s0%2i",HSet.WaveFunctionRootFileName,iOuter+1);
      sprintf(SpecFileName,"%s0%2i",HSet.SpecRootFileName,iOuter+1);
    }
    else if(iOuter<999)
    {
      sprintf(GBathFileName,"%s%3i",HSet.GRootFileName,iOuter+1);
      sprintf(PsiFileName,"%s%3i",HSet.WaveFunctionRootFileName,iOuter+1);
      sprintf(SpecFileName,"%s%3i",HSet.SpecRootFileName,iOuter+1);
    }
    FilePrintSpectra(Spectra, 10, SpecFileName, "w");
    FilePrintPlainTriDiagonalMatrix(GBath,GBathFileName,"w");
    FilePrintWaveFunctionPlainBinary(Psi, PsiFileName, "w");
    /* free variables */
    printf("Start free variables\n");
    for(i=0;i<4;i++)
    {
      FreeOperator(&(VPerturbation[i]));
      FreeTriDiagonalMatrix(&(gPerturbation[i]));
    }
    FreeOperator(&HAnd);
    FreeOperator(&OPESdn);
    FreeOperator(&OPESup);
    FreeOperator(&OIPESdn);
    FreeOperator(&OIPESup);
    FreeTriDiagonalMatrix(&GclMup);
    FreeTriDiagonalMatrix(&GclMdn);
    FreeTriDiagonalMatrix(&GclPup);
    FreeTriDiagonalMatrix(&GclPdn);
    FreeTriDiagonalMatrix(&Gcl);
    FreeTriDiagonalMatrix(&Sigma);
    FreeTriDiagonalMatrix(&GBathNew);
    FreeTriDiagonalMatrix(&GBathOld);
    FreeDoubleTriDiagonalMatrix(&(GAnd.Dt));
    EndTiming("DMFTLoop",&Time);
    
    PrintTiming(Time);
  }
  
  
  
  EndTiming("Total",&Time);
  PrintTiming(Time);
  return 0;
}







