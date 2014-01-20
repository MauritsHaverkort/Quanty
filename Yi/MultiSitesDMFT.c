#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <unistd.h>
#include <stdarg.h>
#include <mkl_cblas.h>
//#include <Accelerate/Accelerate.h>


#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#define omp_get_num_procs() 0
#define omp_get_max_threads() 0
#endif

#include "../Core/StructDef.c"
#include "../Headers/header.h"
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
#include "../Yi/MultiSitesFunctions.c"



int main(int argc, const char *argv[])
{
  TimeType Time;
  MultiSitesDMFTSettings DMFTSet;
  FileReadMultiSitesDMFTSettings(&DMFTSet, argv[argc-1]);
  int NSite = DMFTSet.NSite;
  AndersonMatrixType G0And[NSite], GBathupAnd[NSite], GBathdnAnd[NSite], GClupAnd[NSite], GCldnAnd[NSite];
  ListOfPoles SigmaupLP[NSite], SigmaupLPNew[NSite], SigmaupLPOld[NSite];
  ListOfPoles SigmadnLP[NSite], SigmadnLPNew[NSite], SigmadnLPOld[NSite];
  BADoubleTriDiagonalMatrixType GBathupBA[NSite], GBathupBATrunc[NSite];
  BADoubleTriDiagonalMatrixType GBathdnBA[NSite], GBathdnBATrunc[NSite];
  WaveFunctionType Psi;
  OperatorType HLoc, HAnd, *VPerturbation;
  TriDiagonalMatrixType *gPerturbation;
  uint32_t iDMFT, i;
  double mu, Nelectron, onsiteshiftstart;
  SpectraType Spectra[12*NSite];
  char GBathFileName[MAXStrLength], GClFileName[MAXStrLength], GSigmaFileName[MAXStrLength];
  char PsiFileName[MAXStrLength], SpecFileName[MAXStrLength];
  
  InitializeTiming(&Time);
  BeginTiming("Total", &Time);
  TimeStampAll();
  
  BeginTiming("Init", &Time);
  if (initializeMultiSitesDMFT(argv[argc-1], &DMFTSet, G0And, SigmaupLP, SigmadnLP, &Psi))
  {
    printf("initializeMultiSitesDMFT failed in main\n");
    return 1;
  }
  
  VPerturbation=malloc(4*DMFTSet.NSite*sizeof(OperatorType));
  gPerturbation=malloc(4*DMFTSet.NSite*sizeof(TriDiagonalMatrixType));
  mu = DMFTSet.mu;
  Nelectron = DMFTSet.Nelectron;
  EndTiming("Init", &Time);
  
  if (file_exists("onsiteshift")) {
    onsiteshiftstart = read_oneline("onsiteshift");
    printf("\non site shift read from previous run: %22.15E\n", onsiteshiftstart);
  }
  else
  {
    onsiteshiftstart = mu;
  }
  
  for(iDMFT=DMFTSet.iStart; iDMFT<DMFTSet.iStart+DMFTSet.Ni; iDMFT++)
  {
    BeginTiming("DMFTLoop", &Time);
    {
      printf("\n\n***********************************************************************\n");
      printf("Start of loop %i from %i to %i\n",
             iDMFT+1, DMFTSet.iStart, DMFTSet.iStart+DMFTSet.Ni);
      printf("***********************************************************************\n\n\n");
      
      BeginTiming("Find Groundstate and mu", &Time);
      MultiSitesFindGroundStateFromG0andSigma(&DMFTSet, G0And, SigmaupLP, SigmadnLP,
                                              &mu, &Nelectron, &Psi,
                                              GBathupBA, GBathdnBA,
                                              GBathupAnd, GBathdnAnd,
                                              &Time, &onsiteshiftstart);
      EndTiming("Find Groundstate and mu", &Time);
      
      BeginTiming("Set up the V and g", &Time);
      MultiSitesCreateGBathBATruncandVPerturbation(&DMFTSet, GBathupBA, GBathdnBA,
                                                   GBathupBATrunc, GBathdnBATrunc,
                                                   VPerturbation, gPerturbation);
      EndTiming("Set up the V and g", &Time);
          
      BeginTiming("Calculate G cluster", &Time);
      {
        MultiSitesTestPsiDimensions(&DMFTSet, GBathupBATrunc, GBathdnBATrunc, &Psi);
        //LocalHamiltonian(&HLoc, GBathupBATrunc, GBathdnBATrunc, DMFTSet.NSite, DMFTSet.U, DMFTSet.J, 0);
        DCAHamiltonian(&HLoc, GBathupBATrunc, GBathdnBATrunc, DMFTSet.NSite, DMFTSet.U, 0);
        MultiSitesBADoubleTriDiagonalMatrixToAndersonImpurityOperator(&HAnd, &HLoc, GBathupBATrunc, GBathdnBATrunc, DMFTSet.NSite);
        CreateMultiSitesGclDMFT(&DMFTSet, HAnd, Psi, VPerturbation, gPerturbation, mu, GClupAnd, GCldnAnd, &Time);
      }
      EndTiming("Calculate G cluster", &Time);
      
      BeginTiming("Calculate Selfenergy", &Time);
      for (i=0; i<DMFTSet.NSite; i++)
      {
        SelfEnergyListOfPolesFromAndersonRepresentation(GBathupAnd[i], GClupAnd[i], &SigmaupLPNew[i]);
        ListOfPolesCopy(&SigmaupLPOld[i], SigmaupLP[i]);
        MixListOfPoles(&SigmaupLP[i], SigmaupLPNew[i], DMFTSet.Mix);
        ReduceListOfPolesNumberFermi(&SigmaupLP[i], DMFTSet.NGBath, SigmaupLP[i].mu);
        
        SelfEnergyListOfPolesFromAndersonRepresentation(GBathdnAnd[i], GCldnAnd[i], &SigmadnLPNew[i]);
        ListOfPolesCopy(&SigmadnLPOld[i], SigmadnLP[i]);
        MixListOfPoles(&SigmadnLP[i], SigmadnLPNew[i], DMFTSet.Mix);
        ReduceListOfPolesNumberFermi(&SigmadnLP[i], DMFTSet.NGBath, SigmadnLP[i].mu);
      }
      EndTiming("Calculate Selfenergy", &Time);
      
      printf("Start create spectra\n");fflush(stdout);
      Spectra[0].N        = DMFTSet.NE;
      Spectra[0].Emin     = DMFTSet.Emin;
      Spectra[0].Emax     = DMFTSet.Emax;
      Spectra[0].dE       = (DMFTSet.Emax-DMFTSet.Emin)/DMFTSet.NE;
      Spectra[0].EminPole = DMFTSet.Emin;
      Spectra[0].EmaxPole = DMFTSet.Emax;
      Spectra[0].Gamma    = DMFTSet.Gamma;
      for (i=0;i<12*NSite;i++)
      {
        SpectrumAllocate(Spectra[0], &Spectra[i]);
      }
      for (i=0; i<DMFTSet.NSite; i++)
      {
        AndersonImpurityMatrixToSpectrum(0, GBathupAnd[i],      &(Spectra[12*i+ 0]));                  /*  3 + 24i*/
        AndersonImpurityMatrixToSpectrum(0, GBathdnAnd[i],      &(Spectra[12*i+ 1]));                  /*  5 + 24i*/
        AndersonImpurityMatrixToSpectrum(0, GClupAnd[i],        &(Spectra[12*i+ 2]));                  /*  7 + 24i*/
        AndersonImpurityMatrixToSpectrum(0, GCldnAnd[i],        &(Spectra[12*i+ 3]));                  /*  9 + 24i*/
        AndersonImpurityMatrixToSpectrum(DMFTSet.mu0, G0And[i], &(Spectra[12*i+ 4]));                  /* 11 + 24i*/
        ListOfPolesToSpectrum(0, SigmaupLP[i],                  &(Spectra[12*i+ 5]));                  /* 13 + 24i*/
        ListOfPolesToSpectrum(0, SigmadnLP[i],                  &(Spectra[12*i+ 6]));                  /* 15 + 24i*/
        ListOfPolesToSpectrum(0, SigmaupLPNew[i],               &(Spectra[12*i+ 7]));                  /* 17 + 24i*/
        ListOfPolesToSpectrum(0, SigmadnLPNew[i],               &(Spectra[12*i+ 8]));                  /* 19 + 24i*/
        ListOfPolesToSpectrum(0, SigmaupLPOld[i],               &(Spectra[12*i+ 9]));                  /* 21 + 24i*/
        ListOfPolesToSpectrum(0, SigmadnLPOld[i],               &(Spectra[12*i+10]));                  /* 23 + 24i*/
        AndersonImpurityMatrixToSpectrum(0, G0And[i],           &(Spectra[12*i+11]));                  /* 25 + 24i*/
      }
      
      
      printf("Start write to disk\n"); fflush(stdout);
      
      sprintf(SpecFileName, "%s%03d", DMFTSet.SpecRootFileName, iDMFT+1);
      FilePrintSpectra(Spectra, 12*NSite, SpecFileName, "w");
      
      for (i=0; i<DMFTSet.NSite; i++)
      {
        sprintf(GBathFileName,"%sBathup%03d", DMFTSet.RootFileName[i], iDMFT+1);
        sprintf(GClFileName,"%sClusterup%03d", DMFTSet.RootFileName[i], iDMFT+1);
        sprintf(GSigmaFileName,"%sSigmaup%03d", DMFTSet.RootFileName[i], iDMFT+1);
        FilePrintPlainAndersonMatrix(GBathupAnd[i], GBathFileName, "w");
        FilePrintPlainAndersonMatrix(GClupAnd[i], GClFileName, "w");
        FilePrintPlainListOfPoles(SigmaupLP[i], GSigmaFileName, "w");
        
        sprintf(GBathFileName,"%sBathdn%03d", DMFTSet.RootFileName[i], iDMFT+1);
        sprintf(GClFileName,"%sClusterdn%03d", DMFTSet.RootFileName[i], iDMFT+1);
        sprintf(GSigmaFileName,"%sSigmadn%03d", DMFTSet.RootFileName[i], iDMFT+1);
        FilePrintPlainAndersonMatrix(GBathdnAnd[i], GBathFileName, "w");
        FilePrintPlainAndersonMatrix(GCldnAnd[i], GClFileName, "w");
        FilePrintPlainListOfPoles(SigmadnLP[i], GSigmaFileName, "w");
      }
      
      sprintf(PsiFileName, "%s%03d", DMFTSet.WaveFunctionRootFileName, iDMFT+1);
      FilePrintWaveFunctionPlainBinary(Psi, PsiFileName, "w");
      
      for (i=0; i<DMFTSet.NSite; i++)
      {
        FreeAndersonMatrix(&GBathupAnd[i]);
        FreeAndersonMatrix(&GClupAnd[i]);
        FreeListOfPoles(&SigmaupLPNew[i]);
        FreeListOfPoles(&SigmaupLPOld[i]);
        FreeBADoubleTriDiagonalMatrix(&GBathupBA[i]);
        
        FreeAndersonMatrix(&GBathdnAnd[i]);
        FreeAndersonMatrix(&GCldnAnd[i]);
        FreeListOfPoles(&SigmadnLPNew[i]);
        FreeListOfPoles(&SigmadnLPOld[i]);
        FreeBADoubleTriDiagonalMatrix(&GBathdnBA[i]);
      }
      
      /* subroutine to calculate the spectra
      if ((iDMFT==DMFTSet.iStart+DMFTSet.Ni-1) && (DMFTSet.Spectra)) {
        if (DMFTSpec(DMFTSet, HAnd, Psi, VPerturbation, gPerturbation)) {
          printf("\n========== ERROR: DMFTSpec failed! ==========\n");
          return 1;
        }
      }
      */
      
      for(i=0; i<4*DMFTSet.NSite; i++)
      {
        FreeTriDiagonalMatrix(&gPerturbation[i]);
        FreeOperator(&VPerturbation[i]);
      }
    }
    
    EndTiming("DMFTLoop", &Time);
  }
  
  EndTiming("Total", &Time);
  printf("\n\n***********************************************************************\n\n\n");
  PrintTiming(Time);
  return 0;
}