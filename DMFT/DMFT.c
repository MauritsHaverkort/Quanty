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
#include "../Core/WaveFunctionPerturbationArithmetic.c"
#include "../Core/OperatorArithmetic.c"
#include "../Core/OperatorWaveFunction.c"
#include "../Core/OperatorWaveFunctionPerturbation.c"
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
  TriDiagonalMatrixType G0,GBath,GclPup,GclPdn,GclMup,GclMdn,Gcl,GBathNew,GBathOld,gPerturbation[4],Sigma;
  BADoubleTriDiagonalMatrixType GAnd,GAndTrunk;
  AndersonMatrixType GBathAnd, GclAnd,G0And,GBathNewAnd;
  ListOfPoles SigmaPoles,GclList,GListTemp;
  OperatorType HAnd, OPESdn,OPESup,OIPESdn,OIPESup,VPerturbation[4];
  MatrixType T;
  SpectraType Spectra[10];
  char GBathFileName[255],GClFileName[255],GSigmaFileName[255],PsiFileName[255],SpecFileName[255],Buffer[255];
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
  TriDiagonalMatrixToAndersonMatrix(G0, &T, &G0And);
  FreeMatrix(&T);
  /* do the DMFT loops */
  for(iOuter=HSet.iStart;iOuter<HSet.iStart+HSet.Ni;iOuter++)
  {
    BeginTiming("DMFTLoop",&Time);
    printf("\n\n***********************************************************************\n");
    printf("Start of loop %i from %i to %i\n",iOuter+1,HSet.iStart,HSet.iStart+HSet.Ni);
    printf("***********************************************************************\n\n\n");
    /* Create splitted anderson matrix */
    TriDiagonalMatrixToBADoubleTriDiagonalMatrix(GBath, &GAnd, 2+GBath.N/2, 2+GBath.N/2, Ef, epsilonAnd);
    /* Create GAnd and GAndTrunk as a copy pointing to the same memory */
    GAndTrunk.Ef=GAnd.Ef;
    GAndTrunk.Epsilon=GAnd.Epsilon;
    GAndTrunk.a0=GAnd.a0;
    GAndTrunk.b0=GAnd.b0;
    GAndTrunk.acvB=GAnd.acvB;
    GAndTrunk.acvA=GAnd.acvA;
    GAndTrunk.acvBA=GAnd.acvBA;
    GAndTrunk.bvB0=GAnd.bvB0;
    GAndTrunk.bcB0=GAnd.bcB0;
    GAndTrunk.bvA0=GAnd.bvA0;
    GAndTrunk.bcA0=GAnd.bcA0;
    GAndTrunk.val.N=GAnd.val.N;
    GAndTrunk.val.a=GAnd.val.a;
    GAndTrunk.val.b=GAnd.val.b;
    GAndTrunk.con.N=GAnd.con.N;
    GAndTrunk.con.a=GAnd.con.a;
    GAndTrunk.con.b=GAnd.con.b;
    if(HSet.Dyson)
    {
      /* Create VPerturbation */
      for(i=0;i<4;i++)
      {
        InitOperator(&(VPerturbation[i]));
        VPerturbation[i].NFermion=4+2*(Min(-1+GAnd.val.N,HSet.NBath)+Min(-1+GAnd.con.N,HSet.NBath));
        VPerturbation[i].NBoson=0;
        VPerturbation[i].NWord=Min(VPerturbation[i].NFermion,1+(VPerturbation[i].NFermion-1)/8);
        VPerturbation[i].QComplex=0;
        InitOperatorFixedLength(&(VPerturbation[i]),1,1);
      }
      if(GAnd.val.N>(HSet.NBath+1))
      {
        oadd[0]=CrAnMask+3+HSet.NBath;
        RealOperatorAddLader(&(VPerturbation[0]), oadd, 1, GAnd.val.b[HSet.NBath+1]);
        oadd[0]=CrAnMask+3+Min(-1+GAnd.con.N,HSet.NBath)+2*HSet.NBath;
        RealOperatorAddLader(&(VPerturbation[1]), oadd, 1, GAnd.val.b[HSet.NBath+1]);
      }
      if(GAnd.con.N>HSet.NBath+1)
      {
        oadd[0]=3+Min(-1+GAnd.val.N,HSet.NBath)+HSet.NBath;
        RealOperatorAddLader(&(VPerturbation[2]), oadd, 1, GAnd.con.b[HSet.NBath+1]);
        oadd[0]=3+2*Min(-1+GAnd.val.N,HSet.NBath)+2*HSet.NBath;
        RealOperatorAddLader(&(VPerturbation[3]), oadd, 1, GAnd.con.b[HSet.NBath+1]);
      }
      /* Create gPerturbation */
      if(GAnd.val.N>HSet.NBath+1)
      {
        gPerturbation[0].N=GAnd.val.N-HSet.NBath-1;
        gPerturbation[1].N=GAnd.val.N-HSet.NBath-1;
        InitTriDiagonalMatrix(&(gPerturbation[0]));
        InitTriDiagonalMatrix(&(gPerturbation[1]));
        for(i=1;i<GAnd.val.N-HSet.NBath-1;i++)
        {
          gPerturbation[0].a[i]=-GAnd.val.a[i+HSet.NBath+1];
          gPerturbation[0].b[i]=GAnd.val.b[i+HSet.NBath+1];
          gPerturbation[1].a[i]=-GAnd.val.a[i+HSet.NBath+1];
          gPerturbation[1].b[i]=GAnd.val.b[i+HSet.NBath+1];
        }
        gPerturbation[0].a[GAnd.val.N-HSet.NBath-1]=-GAnd.val.a[GAnd.val.N];
        gPerturbation[1].a[GAnd.val.N-HSet.NBath-1]=-GAnd.val.a[GAnd.val.N];
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
      if(GAnd.con.N>HSet.NBath+1)
      {
        gPerturbation[2].N=GAnd.con.N-HSet.NBath-1;
        gPerturbation[3].N=GAnd.con.N-HSet.NBath-1;
        InitTriDiagonalMatrix(&(gPerturbation[2]));
        InitTriDiagonalMatrix(&(gPerturbation[3]));
        for(i=1;i<GAnd.con.N-HSet.NBath-1;i++)
        {
          gPerturbation[2].a[i]=GAnd.con.a[i+HSet.NBath+1];
          gPerturbation[2].b[i]=GAnd.con.b[i+HSet.NBath+1];
          gPerturbation[3].a[i]=GAnd.con.a[i+HSet.NBath+1];
          gPerturbation[3].b[i]=GAnd.con.b[i+HSet.NBath+1];
        }
        gPerturbation[2].a[GAnd.con.N-HSet.NBath-1]=GAnd.con.a[GAnd.con.N];
        gPerturbation[3].a[GAnd.con.N-HSet.NBath-1]=GAnd.con.a[GAnd.con.N];
      }
      else
      {
        gPerturbation[2].N=1;
        gPerturbation[3].N=1;
        InitTriDiagonalMatrix(&(gPerturbation[2]));
        InitTriDiagonalMatrix(&(gPerturbation[3]));
        gPerturbation[2].a[1]=0;
        gPerturbation[3].a[1]=0;
      }
    }
    /* Create HAnd */
    if(GAnd.con.N<HSet.NBath+1)
    {
      printf("Warning: GAnd.con.N(%i)<HSet.NBath+1(%i+1)\n",GAnd.con.N,HSet.NBath);
    }
    else
    {
      GAndTrunk.con.N=HSet.NBath+1;
    }
    if(GAnd.val.N<HSet.NBath+1)
    {
      printf("Warning: GAnd.val.N(%i)<HSet.NBath+1(%i+1)\n",GAnd.val.N,HSet.NBath);
    }
    else
    {
      GAndTrunk.val.N=HSet.NBath+1;
    }
    mu=HSet.U/2;
    BADoubleTriDiagonalMatrixToAndersonImpurityOperator(&HAnd, GAndTrunk, GAndTrunk, HSet.U, mu);
    /* test if psi has correct dimensions and correct if needed */
    sscanf(Psi.Name,"%s %i %i",Buffer,&PsiNBathC,&PsiNBathV);
    if(PsiNBathC!=GAndTrunk.con.N-1)
    {
      printf("\nBefore diagonalizing:\nThe wavefunction has %i conduction bath sites\nThe Anderson Hamiltonian %i\n",PsiNBathC,GAndTrunk.con.N-1);
      if(PsiNBathC>GAndTrunk.con.N-1)
      {
        printf("Removing %i unoccupied fermions from the wave-function\n",PsiNBathC-GAndTrunk.con.N+1);
        for(i=0;i<PsiNBathC-GAndTrunk.con.N+1;i++)
        {
          WaveFunctionRemoveFermion(Psi, &PsiTemp, Psi.NFermion-1);
          WaveFunctionFree(&Psi);
          WaveFunctionRemoveFermion(PsiTemp, &Psi, 4+PsiNBathV+PsiNBathC-i-1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      else
      {
        printf("Adding %i unoccupied fermions to the wave-function\n",GAndTrunk.con.N-1-PsiNBathC);
        for(i=0;i<GAndTrunk.con.N-1-PsiNBathC;i++)
        {
          WaveFunctionAddFermion(Psi, &PsiTemp, Psi.NFermion,0);
          WaveFunctionFree(&Psi);
          WaveFunctionAddFermion(PsiTemp, &Psi, 4+PsiNBathV+PsiNBathC+i,0);
          WaveFunctionFree(&PsiTemp);
        }
      }
    }
    if(PsiNBathV!=GAndTrunk.val.N-1)
    {
      printf("\nBefore diagonalizing:\nThe wavefunction has %i valence bath sites\nThe Anderson Hamiltonian %i\n",PsiNBathV,GAndTrunk.val.N-1);
      if(PsiNBathV>GAndTrunk.val.N-1)
      {
        printf("Removing %i occupied fermions from the wave-function\n",PsiNBathV-GAndTrunk.val.N+1);
        for(i=0;i<PsiNBathV-GAndTrunk.val.N+1;i++)
        {
          WaveFunctionRemoveFermion(Psi, &PsiTemp, Psi.NFermion-GAndTrunk.con.N);
          WaveFunctionFree(&Psi);
          WaveFunctionRemoveFermion(PsiTemp, &Psi, 4+PsiNBathV-i-1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      else
      {
        printf("Adding %i occupied fermions to the wave-function\n",GAndTrunk.val.N-1-PsiNBathV);
        for(i=0;i<GAndTrunk.val.N-1-PsiNBathV;i++)
        {
          WaveFunctionAddFermion(Psi, &PsiTemp, Psi.NFermion-GAndTrunk.val.N+1,1);
          WaveFunctionFree(&Psi);
          WaveFunctionAddFermion(PsiTemp, &Psi, 4+PsiNBathV+i,1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      printf("Corrected basis\n");
      PrintWaveFunctionReadableTruncated(Psi, 0.1);
    }
    sprintf(Psi.Name,"DMFT %i %i",GAndTrunk.con.N-1,GAndTrunk.val.N-1);
    /* In debug mode print info */
    // PrintBADoubleTriDiagonalMatrix(GAnd);
    // PrintTriDiagonalMatrix(gPerturbation[0]);
    // PrintTriDiagonalMatrix(gPerturbation[1]);
    // PrintTriDiagonalMatrix(gPerturbation[2]);
    // PrintTriDiagonalMatrix(gPerturbation[3]);
    // PrintBADoubleTriDiagonalMatrix(GAndTrunk);
    // PrintOperatorReadable(HAnd);
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
    /* in the last extra loop call DMFTSpectra */
    if((iOuter==HSet.iStart+HSet.Ni-1)&&(HSet.Spectra))
    {
      if(DMFTSpec(HSet, PsiEnergy, Psi, HAnd, VPerturbation, gPerturbation, OPESdn, OPESup, OIPESdn, OIPESup)){printf("DMFTSpec failed\n");return 1;}
    }
    /* if not the last loop continue */
    BeginTiming("Lanczos Tri Diagonalization",&Time);
    OperatorPsi(OPESup,Psi,&OPsi);
    switch (HSet.Dyson) {
      case 1:
        GclMup.N=Max(HSet.NTriDiag,HSet.NGCluster);
        InitTriDiagonalMatrix(&GclMup);
        LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation, gPerturbation, 4, &GclMup);
      break;
      default:
        GclMup.N=HSet.NTriDiag;
        InitTriDiagonalMatrix(&GclMup);
        LanczosTriDiagonalize(HAnd, OPsi, HSet.zeroLancWavefunction, &GclMup, HSet.Gamma);
      break;
    }
    EndTiming("Lanczos Tri Diagonalization",&Time);
    BeginTiming("Lanczos Tri Diagonalization",&Time);
    OperatorPsi(OPESdn,Psi,&OPsi);
    switch (HSet.Dyson) {
      case 1:
        GclMdn.N=Max(HSet.NTriDiag,HSet.NGCluster);
        InitTriDiagonalMatrix(&GclMdn);
        LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation, gPerturbation, 4, &GclMdn);
      break;
      default:
        GclMdn.N=HSet.NTriDiag;
        InitTriDiagonalMatrix(&GclMdn);
        LanczosTriDiagonalize(HAnd, OPsi, HSet.zeroLancWavefunction, &GclMdn, HSet.Gamma);
      break;
    }
    EndTiming("Lanczos Tri Diagonalization",&Time);
    BeginTiming("Lanczos Tri Diagonalization",&Time);
    OperatorPsi(OIPESup,Psi,&OPsi);
    switch (HSet.Dyson) {
      case 1:
        GclPup.N=Max(HSet.NTriDiag,HSet.NGCluster);
        InitTriDiagonalMatrix(&GclPup);
        LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation, gPerturbation, 4, &GclPup);
      break;
      default:
        GclPup.N=HSet.NTriDiag;
        InitTriDiagonalMatrix(&GclPup);
        LanczosTriDiagonalize(HAnd, OPsi, HSet.zeroLancWavefunction, &GclPup, HSet.Gamma);
      break;
    }
    EndTiming("Lanczos Tri Diagonalization",&Time);
    BeginTiming("Lanczos Tri Diagonalization",&Time);
    OperatorPsi(OIPESdn,Psi,&OPsi);
    switch (HSet.Dyson) {
      case 1:
        GclPdn.N=Max(HSet.NTriDiag,HSet.NGCluster);
        InitTriDiagonalMatrix(&GclPdn);
        LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation, gPerturbation, 4, &GclPdn);
      break;
      default:
        GclPdn.N=HSet.NTriDiag;
        InitTriDiagonalMatrix(&GclPdn);
        LanczosTriDiagonalize(HAnd, OPsi, HSet.zeroLancWavefunction, &GclPdn, HSet.Gamma);
      break;
    }
    EndTiming("Lanczos Tri Diagonalization",&Time);
    /* G math - create Gcl as the sum of PES - IPES */
    printf("=== Start G math ===\n");fflush(stdout);
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
    TriDiagonalMatrixToListOfPoles(GclPup, &GclList);
    TriDiagonalMatrixToListOfPoles(GclPdn, &GListTemp);
    AddListOfPoles(&GclList, GListTemp);
    FreeListOfPoles(&GListTemp);
    TriDiagonalMatrixToListOfPoles(GclMup, &GListTemp);
    AddListOfPoles(&GclList, GListTemp);
    FreeListOfPoles(&GListTemp);
    TriDiagonalMatrixToListOfPoles(GclMdn, &GListTemp);
    AddListOfPoles(&GclList, GListTemp);
    FreeListOfPoles(&GListTemp);
    ScaleListOfPoles(&GclList,0.5);
    /* symmetrize */
    SymmetrizeListOfPoles(&GclList);
    printf("Cluster Green function created of size %i ",GclList.N);
    ReduceListOfPolesNumber(&GclList, HSet.NGCluster);
    printf("reduced to %i\n",GclList.N);
    ListOfPolesToTridiagonalMatrix(GclList, &Gcl, GclList.N);
    strcpy(Gcl.Name,"Gcluster");
    Gcl.b[0]=1.0;
    /* make Gcl a physical solution of the full lattice model */
    /* Create self energy (SDMFT) of G0 and Gcl SDMFT = - Gcl^-1 + G0^-1 */
    
    /* Create Gclprime from SDMFT Gcl = (G0^-1 - SDMFT)^-1 */
    
    /* G math - create Cluster Self energy From Gcl and GBath */
    BeginTiming("G math Selfenergy",&Time);
    printf("Start SelfEnergyTriDiagonalRepresentation\n");fflush(stdout);
    TriDiagonalMatrixToAndersonMatrix(Gcl, &T, &GclAnd);
    FreeMatrix(&T);
    TriDiagonalMatrixToAndersonMatrix(GBath, &T, &GBathAnd);
    FreeMatrix(&T);
    SelfEnergyListOfPolesFromAndersonRepresentation(GBathAnd, GclAnd, &SigmaPoles);
    
    /*
    GausianBroadenListOfPoles(&SigmaPoles, 0.1, 10);
    ReduceListOfPolesNumber(&SigmaPoles, SigmaPoles.N/5);
    */
    
    printf("Cluster self energy created of size %i\n",SigmaPoles.N);
    strcpy(SigmaPoles.Name,"Sigma");
    /* copy sigma to tridiagonal format */
    ListOfPolesToTridiagonalMatrix(SigmaPoles, &Sigma, SigmaPoles.N);
    EndTiming("G math Selfenergy",&Time);
    /* G math - create GBath */
    BeginTiming("G math GBath",&Time);
    printf("Start GBathAndListAndRepresentation\n");fflush(stdout);    
    GBathAndTriAndRepresentation(G0And, Sigma, &GBathNewAnd);
    printf("New bath Green function created of size %i ",GBathNewAnd.N);
    strcpy(GBathNewAnd.Name,"Gbath-New");
    ReducePolesAndersonMatrixNumber(&GBathNewAnd, HSet.NGBath);
    printf("reduced to %i\n",GBathNewAnd.N);
    AndersonMatrixToTriDiagonalMatrix(GBathNewAnd, &T, &GBathNew, GBathNewAnd.N);
    FreeMatrix(&T);
    EndTiming("G math GBath",&Time);
    /* calculate GBath for the next loop by mixing GBath with GBathNew */
    printf("Start mixing GBath with GBathNew\n");
    TriDiagonalMatrixCopy(&GBathOld, GBath);
    FreeTriDiagonalMatrix(&GBath);
    GBathOld.N=Min(GBathOld.N,HSet.NGBath);
    strcpy(GBathOld.Name,"Gbath-Old");
    GBathOld.b[0]=sqrt(HSet.Mix);
    GBathNew.b[0]=sqrt(1-HSet.Mix);
    SumGreensfunctionsTriDiagonalRepresentation(GBathNew, GBathOld, &GBath);
    strcpy(GBath.Name,"Gbath");
    GBathOld.b[0]=1;
    GBathNew.b[0]=1;
    /* symmetrize */
    SymmetrizeTridiagonalMatrix(&GBath);
    EndTiming("G math",&Time);
    /* create spectra */
    printf("Start create spectra\n");fflush(stdout);
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
    BADoubleTriDiagonalMatrixToSpectrum(0,GAnd,&(Spectra[3]));       /*  9 */
    TriDiagonalMatrixToSpectrum(0,G0,&(Spectra[4]));                 /* 11 */
    TriDiagonalMatrixToSpectrum(0,Sigma,&(Spectra[5]));              /* 13 */
    AndersonImpurityMatrixToSpectrum(0,GBathNewAnd,&(Spectra[6]));   /* 15 */
    TriDiagonalMatrixToSpectrum(0,GclPup,&(Spectra[7]));             /* 17 */
    TriDiagonalMatrixToSpectrum(0,GclMup,&(Spectra[8]));             /* 19 */
    TriDiagonalMatrixToSpectrum(0,Gcl,&(Spectra[9]));                /* 21 */
    /* write to disk */
    printf("Start write to dist\n");fflush(stdout);
    if(iOuter==0)
    {
      sprintf(GBathFileName,"%s001",HSet.GRootFileName);
      sprintf(GClFileName,"Cl%s001",HSet.GRootFileName);
      sprintf(GSigmaFileName,"Sigma%s001",HSet.GRootFileName);
      sprintf(PsiFileName,"%s001",HSet.WaveFunctionRootFileName);
      sprintf(SpecFileName,"%s001",HSet.SpecRootFileName);
    }
    else if(iOuter<9)
    {
      sprintf(GBathFileName,"%s00%1i",HSet.GRootFileName,iOuter+1);
      sprintf(GClFileName,"Cl%s00%1i",HSet.GRootFileName,iOuter+1);
      sprintf(GSigmaFileName,"Sigma%s00%1i",HSet.GRootFileName,iOuter+1);
      sprintf(PsiFileName,"%s00%1i",HSet.WaveFunctionRootFileName,iOuter+1);
      sprintf(SpecFileName,"%s00%1i",HSet.SpecRootFileName,iOuter+1);
    }
    else if(iOuter<99)
    {
      sprintf(GBathFileName,"%s0%2i",HSet.GRootFileName,iOuter+1);
      sprintf(GClFileName,"Cl%s0%2i",HSet.GRootFileName,iOuter+1);
      sprintf(GSigmaFileName,"Sigma%s0%2i",HSet.GRootFileName,iOuter+1);
      sprintf(PsiFileName,"%s0%2i",HSet.WaveFunctionRootFileName,iOuter+1);
      sprintf(SpecFileName,"%s0%2i",HSet.SpecRootFileName,iOuter+1);
    }
    else if(iOuter<999)
    {
      sprintf(GBathFileName,"%s%3i",HSet.GRootFileName,iOuter+1);
      sprintf(GClFileName,"Cl%s%3i",HSet.GRootFileName,iOuter+1);
      sprintf(GSigmaFileName,"Sigma%s%3i",HSet.GRootFileName,iOuter+1);
      sprintf(PsiFileName,"%s%3i",HSet.WaveFunctionRootFileName,iOuter+1);
      sprintf(SpecFileName,"%s%3i",HSet.SpecRootFileName,iOuter+1);
    }
    FilePrintSpectra(Spectra, 10, SpecFileName, "w");
    FilePrintPlainTriDiagonalMatrix(GBath,GBathFileName,"w");
    FilePrintPlainTriDiagonalMatrix(Gcl,GClFileName,"w");
    FilePrintPlainTriDiagonalMatrix(Sigma,GSigmaFileName,"w");
    FilePrintWaveFunctionPlainBinary(Psi, PsiFileName, "w");
    /* free variables */
    printf("Start free variables\n");fflush(stdout);
    FreeTriDiagonalMatrix(&GclMup);
    FreeTriDiagonalMatrix(&GclMdn);
    FreeTriDiagonalMatrix(&GclPup);
    FreeTriDiagonalMatrix(&GclPdn);
    FreeTriDiagonalMatrix(&Gcl);
    FreeListOfPoles(&SigmaPoles);
    FreeTriDiagonalMatrix(&Sigma);
    FreeTriDiagonalMatrix(&GBathNew);
    FreeTriDiagonalMatrix(&GBathOld);
    FreeBADoubleTriDiagonalMatrix(&GAnd);
    FreeAndersonMatrix(&GBathNewAnd);
    FreeAndersonMatrix(&GBathAnd);
    EndTiming("DMFTLoop",&Time);
    PrintTiming(Time);
    if(HSet.Dyson) for(i=0;i<4;i++)
    {
      FreeOperator(&(VPerturbation[i]));
      FreeTriDiagonalMatrix(&(gPerturbation[i]));
    }
    FreeOperator(&HAnd);
    FreeOperator(&OPESdn);
    FreeOperator(&OPESup);
    FreeOperator(&OIPESdn);
    FreeOperator(&OIPESup);
  }
  EndTiming("Total",&Time);
  PrintTiming(Time);
  return 0;
}

uint32_t DMFTSpec(DMFTSettings const HSet, double const Egrd, WaveFunctionType const psiGrd, OperatorType HAnd, OperatorType const *VPerturbation, TriDiagonalMatrixType const *gPerturbation, OperatorType const OAndn, OperatorType const OAnup, OperatorType const OCrdn, OperatorType const OCrup)
{
  uint16_t oadd[4];
  uint32_t i,iQ,iOuter=HSet.Ni+HSet.iStart-1;
  char SpecFileName[255];
  OperatorType T,Nup,Ndn;
  WaveFunctionType Tpsi;
  SpectraType Spectra[HSet.NQ],SpectraTemp;
  TriDiagonalMatrixType SpectraTri;
  printf("Start of DMFT Spec\n");
  fflush(stdout);
  InitOperator(&Nup);
  Nup.QComplex=0;
  Nup.NFermion=HAnd.NFermion;
  Nup.NBoson=HAnd.NBoson;
  Nup.NWord=HAnd.NWord;
  InitOperatorFixedLength(&Nup, 2, 2);
  InitOperator(&Ndn);
  Ndn.QComplex=0;
  Ndn.NFermion=HAnd.NFermion;
  Ndn.NBoson=HAnd.NBoson;
  Ndn.NWord=HAnd.NWord;
  InitOperatorFixedLength(&Ndn, 2, 2);
  OperatorMultiply2NormalOrder(&Nup, OCrup, OAnup, 'N', 'N', &QDetNotOrderedFcFaBcBa);
  OperatorMultiply2NormalOrder(&Ndn, OCrdn, OAndn, 'N', 'N', &QDetNotOrderedFcFaBcBa);
  /* 1 - spin susceptibility */
  if(HSet.Spectra&1)
  {
    printf("\nStart of DMFT Spec spin susceptibility\n\n");
    fflush(stdout);
    OperatorCopy(&T, Nup);
    OperatorAdd(&T, Ndn, -1.0, 0);
    OperatorScale(&T, 0.5, 0);
    OperatorPsi(T, psiGrd, &Tpsi);
    switch (HSet.Dyson)
    {
      case 1:
        SpectraTri.N=Max(HSet.NTriDiag,HSet.NGCluster);
        InitTriDiagonalMatrix(&SpectraTri);
        LanczosTriDiagonalizeDysonLowMem(HAnd, Tpsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation, gPerturbation, 4, &SpectraTri);
      break;
      default:
        SpectraTri.N=HSet.NTriDiag;
        InitTriDiagonalMatrix(&SpectraTri);
        LanczosTriDiagonalize(HAnd, Tpsi, HSet.zeroLancWavefunction, &SpectraTri, HSet.Gamma);
      break;
    }
    Spectra[0].N       =HSet.NE;
    Spectra[0].Emin    =HSet.Emin;
    Spectra[0].Emax    =HSet.Emax;
    Spectra[0].dE      =(HSet.Emax-HSet.Emin)/HSet.NE;
    Spectra[0].EminPole=HSet.Emin;
    Spectra[0].EmaxPole=HSet.Emax;
    Spectra[0].Gamma   =HSet.Gamma;
    SpectrumAllocate(Spectra[0], &Spectra[0]);
    TriDiagonalMatrixToSpectrum(Egrd, SpectraTri, &Spectra[0]);
    if(iOuter==0)
    {
      sprintf(SpecFileName,"%s001_Spin",HSet.SpecRootFileName);
    }
    else if(iOuter<9)
    {
      sprintf(SpecFileName,"%s00%1i_Spin",HSet.SpecRootFileName,iOuter+1);
    }
    else if(iOuter<99)
    {
      sprintf(SpecFileName,"%s0%2i_Spin",HSet.SpecRootFileName,iOuter+1);
    }
    else if(iOuter<999)
    {
      sprintf(SpecFileName,"%s%3i_Spin",HSet.SpecRootFileName,iOuter+1);
    }
    FilePrintSpectra(Spectra, 1, SpecFileName, "w");
    FreeTriDiagonalMatrix(&SpectraTri);
    SpectrumFree(&Spectra[0]);
    FreeOperator(&T);
  }
  /* 2 - charge susceptibility */
  if(HSet.Spectra&2)
  {
    printf("\nStart of DMFT Spec charge susceptibility\n\n");
    fflush(stdout);
    OperatorCopy(&T, Nup);
    OperatorAdd(&T, Ndn, 1.0, 0);
    OperatorPsi(T, psiGrd, &Tpsi);
    switch (HSet.Dyson)
    {
      case 1:
        SpectraTri.N=Max(HSet.NTriDiag,HSet.NGCluster);
        InitTriDiagonalMatrix(&SpectraTri);
        LanczosTriDiagonalizeDysonLowMem(HAnd, Tpsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation, gPerturbation, 4, &SpectraTri);
        break;
      default:
        SpectraTri.N=HSet.NTriDiag;
        InitTriDiagonalMatrix(&SpectraTri);
        LanczosTriDiagonalize(HAnd, Tpsi, HSet.zeroLancWavefunction, &SpectraTri, HSet.Gamma);
        break;
    }
    Spectra[0].N       =HSet.NE;
    Spectra[0].Emin    =HSet.Emin;
    Spectra[0].Emax    =HSet.Emax;
    Spectra[0].dE      =(HSet.Emax-HSet.Emin)/HSet.NE;
    Spectra[0].EminPole=HSet.Emin;
    Spectra[0].EmaxPole=HSet.Emax;
    Spectra[0].Gamma   =HSet.Gamma;
    SpectrumAllocate(Spectra[0], &Spectra[0]);
    TriDiagonalMatrixToSpectrum(Egrd, SpectraTri, &Spectra[0]);
    if(iOuter==0)
    {
      sprintf(SpecFileName,"%s001_Charge",HSet.SpecRootFileName);
    }
    else if(iOuter<9)
    {
      sprintf(SpecFileName,"%s00%1i_Charge",HSet.SpecRootFileName,iOuter+1);
    }
    else if(iOuter<99)
    {
      sprintf(SpecFileName,"%s0%2i_Charge",HSet.SpecRootFileName,iOuter+1);
    }
    else if(iOuter<999)
    {
      sprintf(SpecFileName,"%s%3i_Charge",HSet.SpecRootFileName,iOuter+1);
    }
    FilePrintSpectra(Spectra, 1, SpecFileName, "w");
    FreeTriDiagonalMatrix(&SpectraTri);
    SpectrumFree(&Spectra[0]);
    FreeOperator(&T);
  }
  /* 4 - XAS */
  if(HSet.Spectra&4)
  {
    printf("\nStart of DMFT Spec XAS\n\n");
    fflush(stdout);
    Spectra[0].N       =HSet.NE;
    Spectra[0].Emin    =HSet.Emin;
    Spectra[0].Emax    =HSet.Emax;
    Spectra[0].dE      =(HSet.Emax-HSet.Emin)/HSet.NE;
    Spectra[0].EminPole=HSet.Emin;
    Spectra[0].EmaxPole=HSet.Emax;
    Spectra[0].Gamma   =HSet.Gamma;
    for(iQ=0;iQ<HSet.NQ;iQ++)
    { /* init spectraTri */
      SpectraTri.N=Max(HSet.NTriDiag,HSet.NGCluster);
      InitTriDiagonalMatrix(&SpectraTri);
      /* shift Onsite energy of Anderson operator */
      oadd[1]=CrAnMask+0; oadd[0]=0; RealOperatorAddLaderNormalOrder(&HAnd,oadd,2,-HSet.Q[iQ],&QDetNotOrderedFcFaBcBa);
      oadd[1]=CrAnMask+1; oadd[0]=1; RealOperatorAddLaderNormalOrder(&HAnd,oadd,2,-HSet.Q[iQ],&QDetNotOrderedFcFaBcBa);
      /* calculate XAS spin up */
      OperatorPsi(OCrup, psiGrd, &Tpsi);
      LanczosTriDiagonalizeDysonLowMem(HAnd, Tpsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation, gPerturbation, 4, &SpectraTri);
      /* transform tridiagonal to spectra */
      SpectrumAllocate(Spectra[0], &Spectra[iQ]);
      TriDiagonalMatrixToSpectrum(Egrd, SpectraTri, &Spectra[iQ]);
      /* calculate XAS spin dn */
      OperatorPsi(OCrdn, psiGrd, &Tpsi);
      LanczosTriDiagonalizeDysonLowMem(HAnd, Tpsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation, gPerturbation, 4, &SpectraTri);
      /* transform tridiagonal to spectra */
      SpectrumAllocate(Spectra[0], &SpectraTemp);
      TriDiagonalMatrixToSpectrum(Egrd, SpectraTri, &SpectraTemp);
      /* sum spectra */
      for(i=0;i<Spectra[0].N;i++)
      {
        Spectra[iQ].Re[i]+=SpectraTemp.Re[i];
        Spectra[iQ].Im[i]+=SpectraTemp.Im[i];
      }
      /* shift Onsite energy of Anderson operator back */
      oadd[1]=CrAnMask+0; oadd[0]=0; RealOperatorAddLaderNormalOrder(&HAnd,oadd,2,HSet.Q[iQ],&QDetNotOrderedFcFaBcBa);
      oadd[1]=CrAnMask+1; oadd[0]=1; RealOperatorAddLaderNormalOrder(&HAnd,oadd,2,HSet.Q[iQ],&QDetNotOrderedFcFaBcBa);
      /* free variables */
      SpectrumFree(&SpectraTemp);
      FreeTriDiagonalMatrix(&SpectraTri);
    }
    /* write spectra to disk */
    if(iOuter==0)
    {
      sprintf(SpecFileName,"%s001_XAS",HSet.SpecRootFileName);
    }
    else if(iOuter<9)
    {
      sprintf(SpecFileName,"%s00%1i_XAS",HSet.SpecRootFileName,iOuter+1);
    }
    else if(iOuter<99)
    {
      sprintf(SpecFileName,"%s0%2i_XAS",HSet.SpecRootFileName,iOuter+1);
    }
    else if(iOuter<999)
    {
      sprintf(SpecFileName,"%s%3i_XAS",HSet.SpecRootFileName,iOuter+1);
    }
    FilePrintSpectra(Spectra, HSet.NQ, SpecFileName, "w");
    /* free variables */
    for(iQ=0;iQ<HSet.NQ;iQ++)
    {
      SpectrumFree(&(Spectra[iQ]));
    }
  }
  /* 8 - core PES */
  if(HSet.Spectra&8)
  {
    printf("\nStart of DMFT Spec core PES\n\n");
    fflush(stdout);
    /* create Identity Operator */
    InitOperator(&T);
    T.QComplex=0;
    T.NFermion=HAnd.NFermion;
    T.NBoson=HAnd.NBoson;
    T.NWord=HAnd.NWord;
    InitOperatorFixedLength(&T, 0, 0);
    RealOperatorAddLaderNormalOrder(&T,oadd,0,1,&QDetNotOrderedFcFaBcBa);
    /* set settings for spectra */
    Spectra[0].N       =HSet.NE;
    Spectra[0].Emin    =HSet.Emin;
    Spectra[0].Emax    =HSet.Emax;
    Spectra[0].dE      =(HSet.Emax-HSet.Emin)/HSet.NE;
    Spectra[0].EminPole=HSet.Emin;
    Spectra[0].EmaxPole=HSet.Emax;
    Spectra[0].Gamma   =HSet.Gamma;
    for(iQ=0;iQ<HSet.NQ;iQ++)
    { /* init spectraTri */
      SpectraTri.N=Max(HSet.NTriDiag,HSet.NGCluster);
      InitTriDiagonalMatrix(&SpectraTri);
      /* shift Onsite energy of Anderson operator */
      oadd[1]=CrAnMask+0; oadd[0]=0; RealOperatorAddLaderNormalOrder(&HAnd,oadd,2,-HSet.Q[iQ],&QDetNotOrderedFcFaBcBa);
      oadd[1]=CrAnMask+1; oadd[0]=1; RealOperatorAddLaderNormalOrder(&HAnd,oadd,2,-HSet.Q[iQ],&QDetNotOrderedFcFaBcBa);
      /* calculate cPES */
      OperatorPsi(T, psiGrd, &Tpsi);
      LanczosTriDiagonalizeDysonLowMem(HAnd, Tpsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation, gPerturbation, 4, &SpectraTri);
      /* transform tridiagonal to spectra */
      SpectrumAllocate(Spectra[0], &Spectra[iQ]);
      TriDiagonalMatrixToSpectrum(Egrd, SpectraTri, &Spectra[iQ]);
      /* shift Onsite energy of Anderson operator back */
      oadd[1]=CrAnMask+0; oadd[0]=0; RealOperatorAddLaderNormalOrder(&HAnd,oadd,2,HSet.Q[iQ],&QDetNotOrderedFcFaBcBa);
      oadd[1]=CrAnMask+1; oadd[0]=1; RealOperatorAddLaderNormalOrder(&HAnd,oadd,2,HSet.Q[iQ],&QDetNotOrderedFcFaBcBa);
      /* free variables */
      FreeTriDiagonalMatrix(&SpectraTri);
    }
    /* write spectra to disk */
    if(iOuter==0)
    {
      sprintf(SpecFileName,"%s001_cPES",HSet.SpecRootFileName);
    }
    else if(iOuter<9)
    {
      sprintf(SpecFileName,"%s00%1i_cPES",HSet.SpecRootFileName,iOuter+1);
    }
    else if(iOuter<99)
    {
      sprintf(SpecFileName,"%s0%2i_cPES",HSet.SpecRootFileName,iOuter+1);
    }
    else if(iOuter<999)
    {
      sprintf(SpecFileName,"%s%3i_cPES",HSet.SpecRootFileName,iOuter+1);
    }
    FilePrintSpectra(Spectra, HSet.NQ, SpecFileName, "w");
    /* free variables */
    for(iQ=0;iQ<HSet.NQ;iQ++)
    {
      SpectrumFree(&(Spectra[iQ]));
    }
  }
  return 0;
}


