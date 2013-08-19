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
  const uint32_t imuMax=20, iHistMax=20, StringLengthMax=255;
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
  BADoubleTriDiagonalMatrixType GAnd,GAndTrunk,BAG0;
  AndersonMatrixType GBathAnd, GclAnd,G0And,GBathNewAnd;
  ListOfPoles SigmaPoles,GclList,GListTemp;
  OperatorType HAnd, HAndGrd, OPESdn,OPESup,OIPESdn,OIPESup,VPerturbation[4],Sx,Sy,Sz,S2;
  MatrixType T;
  SpectraType Spectra[10];
  char GBathFileName[StringLengthMax],GClFileName[StringLengthMax],GSigmaFileName[StringLengthMax],PsiFileName[StringLengthMax],SpecFileName[StringLengthMax],Buffer[StringLengthMax];
  uint32_t iOuter,i,PsiNBathV,PsiNBathC,imu;
  uint16_t oadd[1],*SpinUp,*SpinDn;
  double Ef=0.0,epsilonAnd=0.00001,mu[imuMax+1],G0Occupation,GClOccupation[imuMax];
  double PsiEnergy,PsiQuality,PsiHist[iHistMax],*PsiRho;
  /* read G0 */
  if(FileReadPlainTriDiagonalMatrixType(&G0,HSet.G0FileName)){printf("Fatal Error in main\n Failed to read G0 from file %s\n",HSet.G0FileName);return 1;};
  strcpy(G0.Name,"G0");
  if((HSet.SymOpt)&2)
  {
    printf("Symmetrize G0(w) -= G0(-w)^*\n");
    SymmetrizeTridiagonalMatrix(&G0);
  }
  TriDiagonalMatrixToBADoubleTriDiagonalMatrix(G0, &BAG0, G0.N, G0.N, Ef, epsilonAnd);
  if((HSet.SymOpt)&2)
  {
    printf("Symmetrize G0(w) -= G0(-w)^*\n");
    SymmetrizeBADoubleTriDiagonalMatrix(&BAG0);
  }
  TriDiagonalMatrixToAndersonMatrix(G0, &T, &G0And);
  if((HSet.SymOpt)&2)
  {
    printf("Symmetrize G0(w) -= G0(-w)^*\n");
    SymmetrizeAndersonMatrix(&G0And);
  }
  FreeMatrix(&T);
  G0Occupation=2*Sqr(BAG0.val.b[0]);
  printf("\nRead G0 from file \"%s\"\n",HSet.G0FileName);
  printf("The occupation of the cluster is %22.15E electrons\n",G0Occupation);
  printf("If you don't like this, change G0\n\n");
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
    TriDiagonalMatrixToBADoubleTriDiagonalMatrix(GBath, &GAnd, GBath.N, GBath.N, Ef, epsilonAnd);
    if((HSet.SymOpt)&2)
    {
      printf("Symmetrize GAnd(w) -= GAnd(-w)^*\n");
      printf("Size of GAnd went from %i + %i ",GAnd.val.N,GAnd.con.N);
      SymmetrizeBADoubleTriDiagonalMatrix(&GAnd);
      printf("to %i + %i ",GAnd.val.N,GAnd.con.N);
      GAnd.val.N=Min(GAnd.val.N,HSet.NGBath);
      GAnd.con.N=Min(GAnd.con.N,HSet.NGBath);
      printf("to %i + %i\n",GAnd.val.N,GAnd.con.N);
    }
    if((HSet.SymOpt)&16)
    {
      printf("Enforce GAnd to fullfill the knonw sumrules\n");
      printf("Not yet implemented\n");
    }
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
    /* Create GAndTrunk */
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
    /* Print some info on G0 and GAnd */
    printf("                               G0v                    G0c                    GAndv                  Gandc\n");
    printf("b0                          %22.15E %22.15E %22.15E %22.15E\n",BAG0.val.b[0],BAG0.con.b[0],GAnd.val.b[0],GAnd.con.b[0]);
    printf("b1                          %22.15E %22.15E %22.15E %22.15E\n",BAG0.val.b[1],BAG0.con.b[1],GAnd.val.b[1],GAnd.con.b[1]);
    printf("b2                          %22.15E %22.15E %22.15E %22.15E\n",BAG0.val.b[2],BAG0.con.b[2],GAnd.val.b[2],GAnd.con.b[2]);
    printf("a1                          %22.15E %22.15E %22.15E %22.15E\n",BAG0.val.a[1],BAG0.con.a[1],GAnd.val.a[1],GAnd.con.a[1]);
    printf("a2                          %22.15E %22.15E %22.15E %22.15E\n",BAG0.val.a[2],BAG0.con.a[2],GAnd.val.a[2],GAnd.con.a[2]);
    printf("e_i                         %22.15E                        %22.15E\n",BAG0.acvB,GAnd.acvB);
    printf("e_i^p                       %22.15E                        %22.15E\n",BAG0.acvA,GAnd.acvA);
    printf("t_b                         %22.15E                        %22.15E\n",BAG0.acvBA,GAnd.acvBA);
    printf("t_c                         %22.15E                        %22.15E\n",BAG0.bcB0,GAnd.bcB0);
    printf("t_c^p                       %22.15E                        %22.15E\n",BAG0.bcA0,GAnd.bcA0);
    printf("t_v                         %22.15E                        %22.15E\n",BAG0.bvB0,GAnd.bvA0);
    printf("t_v^p                       %22.15E                        %22.15E\n",BAG0.bvA0,GAnd.bvA0);
    printf("(t_c)^2+(t_c^p)^2           %22.15E                        %22.15E\n",Sqr(BAG0.bcB0)+Sqr(BAG0.bcA0),Sqr(GAnd.bcB0)+Sqr(GAnd.bcA0));
    printf("(t_v)^2+(t_v^p)^2           %22.15E                        %22.15E\n",Sqr(BAG0.bvB0)+Sqr(BAG0.bvA0),Sqr(GAnd.bvB0)+Sqr(GAnd.bvA0));
    printf("(t_c)^2+(t_v)^2+(t_b)^2     %22.15E                        %22.15E\n",Sqr(BAG0.bcB0)+Sqr(BAG0.bvB0)+Sqr(BAG0.acvBA),Sqr(GAnd.bcB0)+Sqr(GAnd.bvB0)+Sqr(GAnd.acvBA));
    printf("(t_c^p)^2+(t_v^p)^2+(t_b)^2 %22.15E                        %22.15E\n",Sqr(BAG0.bcA0)+Sqr(BAG0.bvA0)+Sqr(BAG0.acvBA),Sqr(GAnd.bcA0)+Sqr(GAnd.bvA0)+Sqr(GAnd.acvBA));
    /* find ground-state with correct occupation by adjusting mu */
    BeginTiming("Adjust mu",&Time);
    imu=0;
    if((iOuter==0)||((HSet.SymOpt)&8))
    {
      mu[0]=Sqr(G0Occupation) * HSet.U/2;
    }
    while(imu<imuMax)
    {
      /* Create HAnd */
      BADoubleTriDiagonalMatrixToAndersonImpurityOperator(&HAnd, GAndTrunk, GAndTrunk, HSet.U, mu[imu]);
      /* Set HAndGrd to HAnd */
      OperatorCopy(&HAndGrd,HAnd);
      /* if SymOpt (4) is set keep only H on the impurity and first bath site */
      if((HSet.SymOpt)&4)
      {
        printf("Restricted HAnd for the ground-state to the first bath site\n");
        HAndGrd.l[2].N=8;
      }
      /* if singlet set is true add the (U+1) * S^2 to the operator HAndGrd */
      if((HSet.SymOpt)&1)
      {
        printf("Add %22.15E * S^2 to HAnd for the calculation of the ground-state wave-function\n",1+HSet.U);
        printf("This stabilizes the ground-state by enlarging the singlet triplet energy\n");
        SpinUp=calloc(HAnd.NFermion/2, sizeof(uint16_t));
        SpinDn=calloc(HAnd.NFermion/2, sizeof(uint16_t));
        SpinUp[0]=0;
        SpinUp[1]=2;
        SpinDn[0]=1;
        SpinDn[1]=3;
        for(i=4;i<2+HAnd.NFermion/2;i++)
        {
          SpinUp[i-2]=i;
          SpinDn[i-2]=i+HAnd.NFermion/2-2;
        }
        CreateSpinOperators(HAnd.NFermion, HAnd.NFermion/2, SpinUp, SpinDn, &Sx, &Sy, &Sz, &S2);
        FreeOperator(&Sx);
        FreeOperator(&Sy);
        FreeOperator(&Sz);
        OperatorAddNormalOrder(&HAndGrd, S2, fabs(HSet.U)+1, 0, &QDetNotOrderedFcFaBcBa);
        FreeOperator(&S2);
      }
      /* create ground-state wave-function */
      BeginTiming("Lanczos GroundState",&Time);
      LanczosGroundState(HAndGrd,&Psi,&HSet.zeroLancWavefunction,HSet.zeroLancEnergy,HSet.U,&Psi,0);
      FreeOperator(&HAndGrd);
      EndTiming("Lanczos GroundState",&Time);
      PsiHermitianOperatorPsi(HAnd, Psi, &PsiEnergy);
      PsiHermitianOperatorSqrPsi(HAnd,Psi,&PsiQuality);
      PsiQuality/=Sqr(PsiEnergy);
      printf("Energy                  <psi|H^2|psi>/(<psi|H|psi>^2)\n");
      printf("%22.15E  %22.15E\n",PsiEnergy,PsiQuality);
      /* Start adjusting mu */
      if((HSet.SymOpt)&8)
      {
        imu=imuMax;
      }
      else
      {
        /* calculate cluster occupation */
        PsiRho=calloc(Psi.NFermion,sizeof(double));
        if(PsiRho==NULL){printf("calloc failed in main for PsiRho\n Fatal Error\n");fflush(stdout);return 1;}
        CalculateDensity(PsiRho,Psi);
        GClOccupation[imu]=PsiRho[0]+PsiRho[1];
        free(PsiRho);
        /* Stop if occupation is close enough to occupation of G0 or if imu==imuMax-1 */
        if((imu==imuMax-1)||(Abs(GClOccupation[imu]-G0Occupation)<(1000*DBL_EPSILON))||(imu==0?0:(Abs(GClOccupation[imu]-GClOccupation[imu-1])<(1000*DBL_EPSILON))))
        { /* Stop */
          /* Warning if mu not converged when imu reached imuMax */
          if(Abs(GClOccupation[imu]-G0Occupation)>(1000*DBL_EPSILON))
          {
            printf("WARNING! stop changing imu without converging\n");
          }
          /* print occupation v.s. mu */
          printf("occupation v.s. mu\n");
          printf("  i    mu[i]                  Occupation             Difference\n");
          for(i=0;i<=imu;i++)
          {
            printf("%3i %22.15E %22.15E %22.15E\n",i,mu[i],GClOccupation[i],GClOccupation[i]-G0Occupation);
          }
          mu[0]=mu[imu];
          imu=imuMax;
        }
        else
        { /* make another loop */
          FreeOperator(&HAnd);
          /* interpolate mu[imu+1] */
          if(imu==0)
          {
            if(GClOccupation[0]>G0Occupation)
            {
              mu[1]=mu[0]-0.1*Abs(HSet.U);
            }
            else
            {
              mu[1]=mu[0]+0.1*Abs(HSet.U);
            }
          }
          else
          {
            if(Abs(GClOccupation[imu]-GClOccupation[imu-1])<(100*DBL_EPSILON))
            {
              printf("Convergence creterium for occupation not reached, two mu's give the same occupation\n");
              mu[imu+1]=mu[imu];
              imu=imuMax;
            }
            else
            {
              mu[imu+1]=mu[imu]+(G0Occupation-GClOccupation[imu])*(mu[imu]-mu[imu-1])/(GClOccupation[imu]-GClOccupation[imu-1]);
            }
          }
          printf("\nFilling should be        : %22.15E.\nThe current occupation is: %22.15E.\n", G0Occupation, GClOccupation[imu]);
          printf("mu adjusted from %22.15E to %22.15E\n", mu[imu], mu[imu+1]);
          printf("Restarting Lanczos Loop %i for the %i-th time\n", iOuter,imu+1);
          fflush(stdout);
        }
      }
      imu++;
    }
    EndTiming("Adjust mu", &Time);
    PsiHermitianOperatorPsi(HAnd, Psi, &PsiEnergy);
    PsiHermitianOperatorSqrPsi(HAnd,Psi,&PsiQuality);
    PsiQuality/=Sqr(PsiEnergy);
    PrintWaveFunctionReadableTruncated(Psi, 0.01);
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
    /* print occupations from spectra */
    printf("occupation:   dn                      up\n");
    printf("occupied  : %22.15E %22.15E\n",Sqr(GclMdn.b[0]),Sqr(GclMup.b[0]));
    printf("unoccupied: %22.15E %22.15E\n",Sqr(GclPdn.b[0]),Sqr(GclPup.b[0]));
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
    if((HSet.SymOpt)&2)
    {
      SymmetrizeListOfPoles(&GclList);
    }
    printf("Cluster Green function created of size %i ",GclList.N);
    ReduceListOfPolesNumberFermi(&GclList, HSet.NGCluster,Ef);
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
    printf("Cluster self energy created of size %i with shift %22.15E\n",SigmaPoles.N,SigmaPoles.a0);
    strcpy(SigmaPoles.Name,"Sigma");
    /* copy sigma to tridiagonal format */
    ListOfPolesToTridiagonalMatrix(SigmaPoles, &Sigma, SigmaPoles.N);
    printf("Cluster self energy created in tridiagonal form of size %i with weight %22.15E and shift %22.15E\n",Sigma.N,Sigma.b[0],Sigma.a[0]);
    EndTiming("G math Selfenergy",&Time);
    /* G math - create GBath */
    BeginTiming("G math GBath",&Time);
    BeginTiming("G math GBath AndListAndRepresentation",&Time);
    printf("Start GBathAndListAndRepresentation\n");fflush(stdout);    
    GBathAndTriAndRepresentation(G0And, Sigma, &GBathNewAnd);
    EndTiming("G math GBath AndListAndRepresentation",&Time);
    printf("New bath Green function created of size %i ",GBathNewAnd.N);
    strcpy(GBathNewAnd.Name,"Gbath-New");
    BeginTiming("G math GBath ReducePolesAndersonMatrixNumber",&Time);
    ReducePolesAndersonMatrixNumberFermi(&GBathNewAnd, HSet.NGBath, Ef);
    EndTiming("G math GBath ReducePolesAndersonMatrixNumber",&Time);
    printf("reduced to %i\n",GBathNewAnd.N);
    BeginTiming("G math GBath AndersonMatrixToTriDiagonalMatrix",&Time);
    AndersonMatrixToTriDiagonalMatrix(GBathNewAnd, &T, &GBathNew, GBathNewAnd.N);
    FreeMatrix(&T);
    EndTiming("G math GBath AndersonMatrixToTriDiagonalMatrix",&Time);
    /* symmetrize */
    if((HSet.SymOpt)&2)
    {
      SymmetrizeTridiagonalMatrix(&GBathNew);
    }
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
    GBath.a[0]=HSet.Mix*GBathOld.a[0]+(1-HSet.Mix)*GBathNew.a[0];
    strcpy(GBath.Name,"Gbath");
    GBathOld.b[0]=1;
    GBathNew.b[0]=1;
    if((HSet.SymOpt)&16)
    {
      printf("Enforce GBath to fullfill the knonw sumrules\n");
      printf("Not yet implemented\n");
    }
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
  OperatorType T,Nup,Ndn,Projectd0,Projectd1,Projectd2;
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
  InitOperatorCopySettings(&Projectd0, Nup);
  InitOperatorCopySettings(&Projectd1, Nup);
  InitOperatorCopySettings(&Projectd2, Nup);
  InitOperatorFixedLength(&Projectd0, 4, 2);
  InitOperatorFixedLength(&Projectd1, 4, 2);
  InitOperatorFixedLength(&Projectd2, 4, 2);
  OperatorMultiply2NormalOrder(&Projectd2, Nup, Ndn, 'N', 'N', &QDetNotOrderedFcFaBcBa);
  OperatorAdd(&Projectd1, Nup, 1, 0);
  OperatorAdd(&Projectd1, Ndn, 1, 0);
  OperatorAdd(&Projectd1, Projectd2, -2, 0);
  OperatorAddLader(&Projectd0, oadd, 0, 1, 0);
  OperatorAdd(&Projectd0, Projectd1, -1, 0);
  OperatorAdd(&Projectd0, Projectd2, -1, 0);
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
  /* 4 - charge susceptibility minus Identity */
  if(HSet.Spectra&4)
  {
    printf("\nStart of DMFT Spec charge susceptibility minus Identity Operator \n\n");
    fflush(stdout);
    OperatorCopy(&T, Nup);
    OperatorAdd(&T, Ndn, 1.0, 0);
    OperatorAddLader(&T, oadd, 0, -1, 0);
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
      sprintf(SpecFileName,"%s001_Inelastic_Charge",HSet.SpecRootFileName);
    }
    else if(iOuter<9)
    {
      sprintf(SpecFileName,"%s00%1i_Inelastic_Charge",HSet.SpecRootFileName,iOuter+1);
    }
    else if(iOuter<99)
    {
      sprintf(SpecFileName,"%s0%2i_Inelastic_Charge",HSet.SpecRootFileName,iOuter+1);
    }
    else if(iOuter<999)
    {
      sprintf(SpecFileName,"%s%3i_Inelastic_Charge",HSet.SpecRootFileName,iOuter+1);
    }
    FilePrintSpectra(Spectra, 1, SpecFileName, "w");
    FreeTriDiagonalMatrix(&SpectraTri);
    SpectrumFree(&Spectra[0]);
    FreeOperator(&T);
  }
  /* 8 - charge susceptibility squared */
  if(HSet.Spectra&8)
  {
    printf("\nStart of DMFT Spec NupNdn \n\n");
    fflush(stdout);
    InitOperator(&T);
    T.QComplex=0;
    T.NFermion=HAnd.NFermion;
    T.NBoson=HAnd.NBoson;
    T.NWord=HAnd.NWord;
    InitOperatorFixedLength(&T, 4, 2);
    OperatorMultiply2NormalOrder(&T, Nup, Ndn, 'N', 'N', &QDetNotOrderedFcFaBcBa);
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
      sprintf(SpecFileName,"%s001_NupNdn",HSet.SpecRootFileName);
    }
    else if(iOuter<9)
    {
      sprintf(SpecFileName,"%s00%1i_NupNdn",HSet.SpecRootFileName,iOuter+1);
    }
    else if(iOuter<99)
    {
      sprintf(SpecFileName,"%s0%2i_NupNdn",HSet.SpecRootFileName,iOuter+1);
    }
    else if(iOuter<999)
    {
      sprintf(SpecFileName,"%s%3i_NupNdn",HSet.SpecRootFileName,iOuter+1);
    }
    FilePrintSpectra(Spectra, 1, SpecFileName, "w");
    FreeTriDiagonalMatrix(&SpectraTri);
    SpectrumFree(&Spectra[0]);
    FreeOperator(&T);
  }
  /* 16 - XAS */
  if(HSet.Spectra&16)
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
    {
      printf("Start of XAS itteration for Q=%22.15E (%i/%i)\n",HSet.Q[iQ],iQ,HSet.NQ);fflush(stdout);
      /* init spectraTri */
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
  /* 32 - core PES */
  if(HSet.Spectra&32)
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
    {
      printf("Start of core PES itteration for Q=%22.15E (%i/%i)\n",HSet.Q[iQ],iQ,HSet.NQ);fflush(stdout);
      /* init spectraTri */
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
  
  
  
  /* 64 - XAS - d0 */
  if(HSet.Spectra&64)
  {
    printf("\nStart of DMFT Spec XAS d0\n\n");
    fflush(stdout);
    Spectra[0].N       =HSet.NE;
    Spectra[0].Emin    =HSet.Emin;
    Spectra[0].Emax    =HSet.Emax;
    Spectra[0].dE      =(HSet.Emax-HSet.Emin)/HSet.NE;
    Spectra[0].EminPole=HSet.Emin;
    Spectra[0].EmaxPole=HSet.Emax;
    Spectra[0].Gamma   =HSet.Gamma;
    for(iQ=0;iQ<HSet.NQ;iQ++)
    {
      printf("Start of XAS itteration for Q=%22.15E (%i/%i)\n",HSet.Q[iQ],iQ,HSet.NQ);fflush(stdout);
      /* init spectraTri */
      SpectraTri.N=Max(HSet.NTriDiag,HSet.NGCluster);
      InitTriDiagonalMatrix(&SpectraTri);
      /* shift Onsite energy of Anderson operator */
      oadd[1]=CrAnMask+0; oadd[0]=0; RealOperatorAddLaderNormalOrder(&HAnd,oadd,2,-HSet.Q[iQ],&QDetNotOrderedFcFaBcBa);
      oadd[1]=CrAnMask+1; oadd[0]=1; RealOperatorAddLaderNormalOrder(&HAnd,oadd,2,-HSet.Q[iQ],&QDetNotOrderedFcFaBcBa);
      /* calculate XAS spin up */
      InitOperator(&T);
      T.QComplex=0;
      T.NFermion=HAnd.NFermion;
      T.NBoson=HAnd.NBoson;
      T.NWord=HAnd.NWord;
      InitOperatorFixedLength(&T, 2, 2);
      InitOperatorFixedLength(&T, 4, 2);
      InitOperatorFixedLength(&T, 6, 2);
      OperatorMultiply2NormalOrder(&T, OCrup, Projectd0, 'N', 'N', &QDetNotOrderedFcFaBcBa);
      OperatorPsi(T, psiGrd, &Tpsi);
      FreeOperator(&T);
      LanczosTriDiagonalizeDysonLowMem(HAnd, Tpsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation, gPerturbation, 4, &SpectraTri);
      /* transform tridiagonal to spectra */
      SpectrumAllocate(Spectra[0], &Spectra[iQ]);
      TriDiagonalMatrixToSpectrum(Egrd, SpectraTri, &Spectra[iQ]);
      /* calculate XAS spin dn */
      InitOperator(&T);
      T.QComplex=0;
      T.NFermion=HAnd.NFermion;
      T.NBoson=HAnd.NBoson;
      T.NWord=HAnd.NWord;
      InitOperatorFixedLength(&T, 2, 2);
      InitOperatorFixedLength(&T, 4, 2);
      InitOperatorFixedLength(&T, 6, 2);
      OperatorMultiply2NormalOrder(&T, OCrdn, Projectd0, 'N', 'N', &QDetNotOrderedFcFaBcBa);
      OperatorPsi(T, psiGrd, &Tpsi);
      FreeOperator(&T);
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
      sprintf(SpecFileName,"%s001_XAS_d0",HSet.SpecRootFileName);
    }
    else if(iOuter<9)
    {
      sprintf(SpecFileName,"%s00%1i_XAS_d0",HSet.SpecRootFileName,iOuter+1);
    }
    else if(iOuter<99)
    {
      sprintf(SpecFileName,"%s0%2i_XAS_d0",HSet.SpecRootFileName,iOuter+1);
    }
    else if(iOuter<999)
    {
      sprintf(SpecFileName,"%s%3i_XAS_d0",HSet.SpecRootFileName,iOuter+1);
    }
    FilePrintSpectra(Spectra, HSet.NQ, SpecFileName, "w");
    /* free variables */
    for(iQ=0;iQ<HSet.NQ;iQ++)
    {
      SpectrumFree(&(Spectra[iQ]));
    }
  }
  
  
  
  /* 128 - XAS - d1 */
  if(HSet.Spectra&128)
  {
    printf("\nStart of DMFT Spec XAS d1\n\n");
    fflush(stdout);
    Spectra[0].N       =HSet.NE;
    Spectra[0].Emin    =HSet.Emin;
    Spectra[0].Emax    =HSet.Emax;
    Spectra[0].dE      =(HSet.Emax-HSet.Emin)/HSet.NE;
    Spectra[0].EminPole=HSet.Emin;
    Spectra[0].EmaxPole=HSet.Emax;
    Spectra[0].Gamma   =HSet.Gamma;
    for(iQ=0;iQ<HSet.NQ;iQ++)
    {
      printf("Start of XAS itteration for Q=%22.15E (%i/%i)\n",HSet.Q[iQ],iQ,HSet.NQ);fflush(stdout);
      /* init spectraTri */
      SpectraTri.N=Max(HSet.NTriDiag,HSet.NGCluster);
      InitTriDiagonalMatrix(&SpectraTri);
      /* shift Onsite energy of Anderson operator */
      oadd[1]=CrAnMask+0; oadd[0]=0; RealOperatorAddLaderNormalOrder(&HAnd,oadd,2,-HSet.Q[iQ],&QDetNotOrderedFcFaBcBa);
      oadd[1]=CrAnMask+1; oadd[0]=1; RealOperatorAddLaderNormalOrder(&HAnd,oadd,2,-HSet.Q[iQ],&QDetNotOrderedFcFaBcBa);
      /* calculate XAS spin up */
      InitOperator(&T);
      T.QComplex=0;
      T.NFermion=HAnd.NFermion;
      T.NBoson=HAnd.NBoson;
      T.NWord=HAnd.NWord;
      InitOperatorFixedLength(&T, 2, 2);
      InitOperatorFixedLength(&T, 4, 2);
      InitOperatorFixedLength(&T, 6, 2);
      OperatorMultiply2NormalOrder(&T, OCrup, Projectd1, 'N', 'N', &QDetNotOrderedFcFaBcBa);
      OperatorPsi(T, psiGrd, &Tpsi);
      FreeOperator(&T);
      LanczosTriDiagonalizeDysonLowMem(HAnd, Tpsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation, gPerturbation, 4, &SpectraTri);
      /* transform tridiagonal to spectra */
      SpectrumAllocate(Spectra[0], &Spectra[iQ]);
      TriDiagonalMatrixToSpectrum(Egrd, SpectraTri, &Spectra[iQ]);
      /* calculate XAS spin dn */
      InitOperator(&T);
      T.QComplex=0;
      T.NFermion=HAnd.NFermion;
      T.NBoson=HAnd.NBoson;
      T.NWord=HAnd.NWord;
      InitOperatorFixedLength(&T, 2, 2);
      InitOperatorFixedLength(&T, 4, 2);
      InitOperatorFixedLength(&T, 6, 2);
      OperatorMultiply2NormalOrder(&T, OCrdn, Projectd1, 'N', 'N', &QDetNotOrderedFcFaBcBa);
      OperatorPsi(T, psiGrd, &Tpsi);
      FreeOperator(&T);
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
      sprintf(SpecFileName,"%s001_XAS_d1",HSet.SpecRootFileName);
    }
    else if(iOuter<9)
    {
      sprintf(SpecFileName,"%s00%1i_XAS_d1",HSet.SpecRootFileName,iOuter+1);
    }
    else if(iOuter<99)
    {
      sprintf(SpecFileName,"%s0%2i_XAS_d1",HSet.SpecRootFileName,iOuter+1);
    }
    else if(iOuter<999)
    {
      sprintf(SpecFileName,"%s%3i_XAS_d1",HSet.SpecRootFileName,iOuter+1);
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


