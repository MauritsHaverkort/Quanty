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
#include "../Yi/DMFTFunctions.c"




int main(int argc, const char *argv[])
{
  TimeType Time;
  DMFTSettings DMFTSet;
  AndersonMatrixType G0And, GBathAnd, GClAnd;
  ListOfPoles SigmaLP, SigmaLPNew, SigmaLPOld;
  BADoubleTriDiagonalMatrixType GBathBA, GBathBATrunc;
  WaveFunctionType Psi;
  OperatorType *VPerturbation, HAnd;
  TriDiagonalMatrixType *gPerturbation;
  uint32_t iOuter, i;
  double mu, Nelectron, onsiteshiftstart;
  SpectraType Spectra[10];
  char   GBathFileName[MAXStrLength], GClFileName[MAXStrLength], GSigmaFileName[MAXStrLength],
    PsiFileName[MAXStrLength], SpecFileName[MAXStrLength];
  double PsiEnergy;
  OperatorType OPESdn, OPESup,OIPESdn,OIPESup;

  InitializeTiming(&Time);
  BeginTiming("Total", &Time);
  TimeStampAll();
  BeginTiming("Init", &Time);
  if (initializeDMFT(argv[argc - 1], &DMFTSet, &G0And, &SigmaLP, &Psi))
  {
    printf("initializeDMFT failed in main\n");
    return 1;
  }
  VPerturbation = malloc(2 * DMFTSet.NSite * sizeof(OperatorType));
  gPerturbation = malloc(2 * DMFTSet.NSite * sizeof(TriDiagonalMatrixType));
  mu = DMFTSet.mu;
  Nelectron = DMFTSet.Nelectron;
  EndTiming("Init", &Time);
  onsiteshiftstart = mu;
  for (iOuter = DMFTSet.iStart; iOuter < DMFTSet.iStart + DMFTSet.Ni; iOuter++)
  {
    BeginTiming("DMFTLoop", &Time);
    {
      printf("\n\n***********************************************************************\n");
      printf("Start of loop %i from %i to %i\n", iOuter + 1, DMFTSet.iStart, DMFTSet.iStart + DMFTSet.Ni);
      printf("***********************************************************************\n\n\n");

      BeginTiming("Find Groundstate and mu", &Time);
      FindGroundStateFromG0andSigma2(DMFTSet, G0And, SigmaLP, &mu, &Nelectron, &Psi, &GBathBA, &GBathAnd, &Time,
                                    &onsiteshiftstart);
      EndTiming("Find Groundstate and mu", &Time);
      BeginTiming("Set up the V and g", &Time);
      createGBathBATruncandVPerturbation(DMFTSet, GBathBA, &GBathBATrunc, VPerturbation, gPerturbation);
      EndTiming("Set up the V and g", &Time);

      BeginTiming("Calculate G cluster", &Time);
      {
        TestIfPsiHasCorrectDimensions(GBathBATrunc, &Psi);
        BADoubleTriDiagonalMatrixToAndersonImpurityOperator(&HAnd, GBathBATrunc, GBathBATrunc, DMFTSet.U, 0);
        CreateGclDMFT(DMFTSet, HAnd, Psi, VPerturbation, gPerturbation, mu, &GClAnd, &Time);
      }
      EndTiming("Calculate G cluster", &Time);
      BeginTiming("Calculate Selfenergy", &Time);
      SelfEnergyListOfPolesFromAndersonRepresentation(GBathAnd, GClAnd, &SigmaLPNew);
      ListOfPolesCopy(&SigmaLPOld, SigmaLP);
      MixListOfPoles(&SigmaLP, SigmaLPNew, DMFTSet.Mix);

      ReduceListOfPolesNumberFermi(&SigmaLP, DMFTSet.NGBath, SigmaLP.mu);

      EndTiming("Calculate Selfenergy", &Time);

      printf("Start create spectra\n");
      fflush(stdout);
      Spectra[0].N = DMFTSet.NE;
      Spectra[0].Emin = DMFTSet.Emin;
      Spectra[0].Emax = DMFTSet.Emax;
      Spectra[0].dE = (DMFTSet.Emax - DMFTSet.Emin) / DMFTSet.NE;
      Spectra[0].EminPole = DMFTSet.Emin;
      Spectra[0].EmaxPole = DMFTSet.Emax;
      Spectra[0].Gamma = DMFTSet.Gamma;
      for (i = 0; i < 7; i++)
      {
        SpectrumAllocate(Spectra[0], &Spectra[i]);
      }
      AndersonImpurityMatrixToSpectrum(0, GBathAnd, &(Spectra[0])); /*  3 */
      AndersonImpurityMatrixToSpectrum(0, GClAnd, &(Spectra[1])); /*  5 */
      AndersonImpurityMatrixToSpectrum(DMFTSet.mu0, G0And, &(Spectra[2]));  /*  7 */
      ListOfPolesToSpectrum(0, SigmaLP, &(Spectra[3])); /*  9 */
      ListOfPolesToSpectrum(0, SigmaLPNew, &(Spectra[4]));  /* 11 */
      ListOfPolesToSpectrum(0, SigmaLPOld, &(Spectra[5]));  /* 13 */
      AndersonImpurityMatrixToSpectrum(0, G0And, &(Spectra[6]));  /* 15 */

      printf("Start write to disk\n");
      fflush(stdout);
      sprintf(GBathFileName, "%sBath%03i", DMFTSet.RootFileName, iOuter + 1);
      sprintf(GClFileName, "%sCluster%03i", DMFTSet.RootFileName, iOuter + 1);
      sprintf(GSigmaFileName, "%sSigma%03i", DMFTSet.RootFileName, iOuter + 1);
      sprintf(PsiFileName, "%s%03i", DMFTSet.WaveFunctionRootFileName, iOuter + 1);
      sprintf(SpecFileName, "%s%03i", DMFTSet.SpecRootFileName, iOuter + 1);
      FilePrintSpectra(Spectra, 7, SpecFileName, "w");
      FilePrintPlainAndersonMatrix(GBathAnd, GBathFileName, "w");
      FilePrintPlainAndersonMatrix(GClAnd, GClFileName, "w");
      FilePrintPlainListOfPoles(SigmaLP, GSigmaFileName, "w");
      FilePrintWaveFunctionPlainBinary(Psi, PsiFileName, "w");

      FreeAndersonMatrix(&GBathAnd);
      FreeAndersonMatrix(&GClAnd);
      FreeListOfPoles(&SigmaLPNew);
      FreeListOfPoles(&SigmaLPOld);
    }
    EndTiming("DMFTLoop", &Time);
    if((DMFTSet.Spectra&1)||(DMFTSet.Spectra&&(iOuter == DMFTSet.iStart + DMFTSet.Ni -1)))
    {
      printf("\n\nCalculate extra spectra after the DMFT loops\n\n");
      BeginTiming("ExtraSpectra", &Time);
      /* Get ground-state energy needed for shift */
      PsiHermitianOperatorPsi(HAnd, Psi, &PsiEnergy);
      /* create operators */
      CreatePESdnOperatorDMFT(&OPESdn, HAnd.NFermion);
      CreatePESupOperatorDMFT(&OPESup, HAnd.NFermion);
      CreateIPESdnOperatorDMFT(&OIPESdn, HAnd.NFermion);
      CreateIPESupOperatorDMFT(&OIPESup, HAnd.NFermion);
      if (DMFTSpec(DMFTSet, iOuter, PsiEnergy, Psi, HAnd, VPerturbation, gPerturbation, OPESdn, OPESup, OIPESdn, OIPESup))
      {
        printf("DMFTSpec failed\n");
        return 1;
      }
      FreeOperator(&OPESdn);
      FreeOperator(&OPESup);
      FreeOperator(&OIPESdn);
      FreeOperator(&OIPESup);
      EndTiming("ExtraSpectra", &Time);
    }
    FreeOperator(&HAnd);
    FreeBADoubleTriDiagonalMatrix(&GBathBA);
    for (i = 0; i < 2 * DMFTSet.NSite; i++)
    {
      FreeTriDiagonalMatrix(&gPerturbation[i]);
      FreeOperator(&VPerturbation[i]);
    }
  }
  EndTiming("Total", &Time);
  printf("\n\n***********************************************************************\n\n\n");
  PrintTiming(Time);
  return 0;
}



uint32_t DMFTSpec(DMFTSettings const HSet, uint32_t const iOuter, double const Egrd, WaveFunctionType const psiGrd, OperatorType HAnd,
                  OperatorType const *VPerturbation, TriDiagonalMatrixType const *gPerturbation,
                  OperatorType const OAndn, OperatorType const OAnup, OperatorType const OCrdn,
                  OperatorType const OCrup)
{
  uint16_t oadd[4];
  uint32_t i, iQ;
  char   SpecFileName[255];
  OperatorType T, Nup, Ndn, Projectd0, Projectd1, Projectd2;
  WaveFunctionType Tpsi;
  SpectraType Spectra[HSet.NQ], SpectraTemp;
  TriDiagonalMatrixType SpectraTri;
  printf("Start of DMFT Spec\n");
  fflush(stdout);
  InitOperator(&Nup);
  Nup.QComplex = 0;
  Nup.NFermion = HAnd.NFermion;
  Nup.NBoson = HAnd.NBoson;
  Nup.NWord = HAnd.NWord;
  InitOperatorFixedLength(&Nup, 2, 2);
  InitOperator(&Ndn);
  Ndn.QComplex = 0;
  Ndn.NFermion = HAnd.NFermion;
  Ndn.NBoson = HAnd.NBoson;
  Ndn.NWord = HAnd.NWord;
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
  /* 2 - spin susceptibility */
  if (HSet.Spectra & 2)
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
      SpectraTri.N = Max(HSet.NTriDiag, HSet.NGCluster);
      InitTriDiagonalMatrix(&SpectraTri);
      LanczosTriDiagonalizeDysonLowMem(HAnd, Tpsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation,
                                       gPerturbation, 4, &SpectraTri);
      break;
    default:
      SpectraTri.N = HSet.NTriDiag;
      InitTriDiagonalMatrix(&SpectraTri);
      LanczosTriDiagonalize(HAnd, Tpsi, HSet.zeroLancWavefunction, &SpectraTri, HSet.Gamma);
      break;
    }
    Spectra[0].N = HSet.NE;
    Spectra[0].Emin = HSet.Emin;
    Spectra[0].Emax = HSet.Emax;
    Spectra[0].dE = (HSet.Emax - HSet.Emin) / HSet.NE;
    Spectra[0].EminPole = HSet.Emin;
    Spectra[0].EmaxPole = HSet.Emax;
    Spectra[0].Gamma = HSet.Gamma;
    SpectrumAllocate(Spectra[0], &Spectra[0]);
    TriDiagonalMatrixToSpectrum(Egrd, SpectraTri, &Spectra[0]);
    sprintf(SpecFileName, "%s%03i_Spin", HSet.SpecRootFileName, iOuter + 1);
    FilePrintSpectra(Spectra, 1, SpecFileName, "w");
    FreeTriDiagonalMatrix(&SpectraTri);
    SpectrumFree(&Spectra[0]);
    FreeOperator(&T);
  }
  /* 4 - charge susceptibility */
  if (HSet.Spectra & 4)
  {
    printf("\nStart of DMFT Spec charge susceptibility\n\n");
    fflush(stdout);
    OperatorCopy(&T, Nup);
    OperatorAdd(&T, Ndn, 1.0, 0);
    OperatorPsi(T, psiGrd, &Tpsi);
    switch (HSet.Dyson)
    {
    case 1:
      SpectraTri.N = Max(HSet.NTriDiag, HSet.NGCluster);
      InitTriDiagonalMatrix(&SpectraTri);
      LanczosTriDiagonalizeDysonLowMem(HAnd, Tpsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation,
                                       gPerturbation, 4, &SpectraTri);
      break;
    default:
      SpectraTri.N = HSet.NTriDiag;
      InitTriDiagonalMatrix(&SpectraTri);
      LanczosTriDiagonalize(HAnd, Tpsi, HSet.zeroLancWavefunction, &SpectraTri, HSet.Gamma);
      break;
    }
    Spectra[0].N = HSet.NE;
    Spectra[0].Emin = HSet.Emin;
    Spectra[0].Emax = HSet.Emax;
    Spectra[0].dE = (HSet.Emax - HSet.Emin) / HSet.NE;
    Spectra[0].EminPole = HSet.Emin;
    Spectra[0].EmaxPole = HSet.Emax;
    Spectra[0].Gamma = HSet.Gamma;
    SpectrumAllocate(Spectra[0], &Spectra[0]);
    TriDiagonalMatrixToSpectrum(Egrd, SpectraTri, &Spectra[0]);
    sprintf(SpecFileName, "%s%03i_Charge", HSet.SpecRootFileName, iOuter + 1);
    FilePrintSpectra(Spectra, 1, SpecFileName, "w");
    FreeTriDiagonalMatrix(&SpectraTri);
    SpectrumFree(&Spectra[0]);
    FreeOperator(&T);
  }
  /* 8 - charge susceptibility minus Identity */
  if (HSet.Spectra & 8)
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
      SpectraTri.N = Max(HSet.NTriDiag, HSet.NGCluster);
      InitTriDiagonalMatrix(&SpectraTri);
      LanczosTriDiagonalizeDysonLowMem(HAnd, Tpsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation,
                                       gPerturbation, 4, &SpectraTri);
      break;
    default:
      SpectraTri.N = HSet.NTriDiag;
      InitTriDiagonalMatrix(&SpectraTri);
      LanczosTriDiagonalize(HAnd, Tpsi, HSet.zeroLancWavefunction, &SpectraTri, HSet.Gamma);
      break;
    }
    Spectra[0].N = HSet.NE;
    Spectra[0].Emin = HSet.Emin;
    Spectra[0].Emax = HSet.Emax;
    Spectra[0].dE = (HSet.Emax - HSet.Emin) / HSet.NE;
    Spectra[0].EminPole = HSet.Emin;
    Spectra[0].EmaxPole = HSet.Emax;
    Spectra[0].Gamma = HSet.Gamma;
    SpectrumAllocate(Spectra[0], &Spectra[0]);
    TriDiagonalMatrixToSpectrum(Egrd, SpectraTri, &Spectra[0]);
    sprintf(SpecFileName, "%s%03i_Inelastic_Charge", HSet.SpecRootFileName, iOuter + 1);
    FilePrintSpectra(Spectra, 1, SpecFileName, "w");
    FreeTriDiagonalMatrix(&SpectraTri);
    SpectrumFree(&Spectra[0]);
    FreeOperator(&T);
  }
  /* 16 - charge susceptibility squared */
  if (HSet.Spectra & 16)
  {
    printf("\nStart of DMFT Spec NupNdn \n\n");
    fflush(stdout);
    InitOperator(&T);
    T.QComplex = 0;
    T.NFermion = HAnd.NFermion;
    T.NBoson = HAnd.NBoson;
    T.NWord = HAnd.NWord;
    InitOperatorFixedLength(&T, 4, 2);
    OperatorMultiply2NormalOrder(&T, Nup, Ndn, 'N', 'N', &QDetNotOrderedFcFaBcBa);
    OperatorPsi(T, psiGrd, &Tpsi);
    switch (HSet.Dyson)
    {
    case 1:
      SpectraTri.N = Max(HSet.NTriDiag, HSet.NGCluster);
      InitTriDiagonalMatrix(&SpectraTri);
      LanczosTriDiagonalizeDysonLowMem(HAnd, Tpsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation,
                                       gPerturbation, 4, &SpectraTri);
      break;
    default:
      SpectraTri.N = HSet.NTriDiag;
      InitTriDiagonalMatrix(&SpectraTri);
      LanczosTriDiagonalize(HAnd, Tpsi, HSet.zeroLancWavefunction, &SpectraTri, HSet.Gamma);
      break;
    }
    Spectra[0].N = HSet.NE;
    Spectra[0].Emin = HSet.Emin;
    Spectra[0].Emax = HSet.Emax;
    Spectra[0].dE = (HSet.Emax - HSet.Emin) / HSet.NE;
    Spectra[0].EminPole = HSet.Emin;
    Spectra[0].EmaxPole = HSet.Emax;
    Spectra[0].Gamma = HSet.Gamma;
    SpectrumAllocate(Spectra[0], &Spectra[0]);
    TriDiagonalMatrixToSpectrum(Egrd, SpectraTri, &Spectra[0]);
    sprintf(SpecFileName, "%s%03i_NupNdn", HSet.SpecRootFileName, iOuter + 1);
    FilePrintSpectra(Spectra, 1, SpecFileName, "w");
    FreeTriDiagonalMatrix(&SpectraTri);
    SpectrumFree(&Spectra[0]);
    FreeOperator(&T);
  }
  /* 32 - XAS */
  if (HSet.Spectra & 32)
  {
    printf("\nStart of DMFT Spec XAS\n\n");
    fflush(stdout);
    Spectra[0].N = HSet.NE;
    Spectra[0].Emin = HSet.Emin;
    Spectra[0].Emax = HSet.Emax;
    Spectra[0].dE = (HSet.Emax - HSet.Emin) / HSet.NE;
    Spectra[0].EminPole = HSet.Emin;
    Spectra[0].EmaxPole = HSet.Emax;
    Spectra[0].Gamma = HSet.Gamma;
    for (iQ = 0; iQ < HSet.NQ; iQ++)
    {
      printf("Start of XAS itteration for Q=%22.15E (%i/%i)\n", HSet.Q[iQ], iQ, HSet.NQ);
      fflush(stdout);
      /* init spectraTri */
      SpectraTri.N = Max(HSet.NTriDiag, HSet.NGCluster);
      InitTriDiagonalMatrix(&SpectraTri);
      /* shift Onsite energy of Anderson operator */
      oadd[1] = CrAnMask + 0;
      oadd[0] = 0;
      RealOperatorAddLaderNormalOrder(&HAnd, oadd, 2, -HSet.Q[iQ], &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + 1;
      oadd[0] = 1;
      RealOperatorAddLaderNormalOrder(&HAnd, oadd, 2, -HSet.Q[iQ], &QDetNotOrderedFcFaBcBa);
      /* calculate XAS spin up */
      OperatorPsi(OCrup, psiGrd, &Tpsi);
      LanczosTriDiagonalizeDysonLowMem(HAnd, Tpsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation,
                                       gPerturbation, 4, &SpectraTri);
      /* transform tridiagonal to spectra */
      SpectrumAllocate(Spectra[0], &Spectra[iQ]);
      TriDiagonalMatrixToSpectrum(Egrd, SpectraTri, &Spectra[iQ]);
      /* calculate XAS spin dn */
      OperatorPsi(OCrdn, psiGrd, &Tpsi);
      LanczosTriDiagonalizeDysonLowMem(HAnd, Tpsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation,
                                       gPerturbation, 4, &SpectraTri);
      /* transform tridiagonal to spectra */
      SpectrumAllocate(Spectra[0], &SpectraTemp);
      TriDiagonalMatrixToSpectrum(Egrd, SpectraTri, &SpectraTemp);
      /* sum spectra */
      for (i = 0; i < Spectra[0].N; i++)
      {
        Spectra[iQ].Re[i] += SpectraTemp.Re[i];
        Spectra[iQ].Im[i] += SpectraTemp.Im[i];
      }
      /* shift Onsite energy of Anderson operator back */
      oadd[1] = CrAnMask + 0;
      oadd[0] = 0;
      RealOperatorAddLaderNormalOrder(&HAnd, oadd, 2, HSet.Q[iQ], &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + 1;
      oadd[0] = 1;
      RealOperatorAddLaderNormalOrder(&HAnd, oadd, 2, HSet.Q[iQ], &QDetNotOrderedFcFaBcBa);
      /* free variables */
      SpectrumFree(&SpectraTemp);
      FreeTriDiagonalMatrix(&SpectraTri);
    }
    /* write spectra to disk */
    sprintf(SpecFileName, "%s%03i_XAS", HSet.SpecRootFileName, iOuter + 1);
    FilePrintSpectra(Spectra, HSet.NQ, SpecFileName, "w");
    /* free variables */
    for (iQ = 0; iQ < HSet.NQ; iQ++)
    {
      SpectrumFree(&(Spectra[iQ]));
    }
  }
  /* 64 - core PES */
  if (HSet.Spectra & 64)
  {
    printf("\nStart of DMFT Spec core PES\n\n");
    fflush(stdout);
    /* create Identity Operator */
    InitOperator(&T);
    T.QComplex = 0;
    T.NFermion = HAnd.NFermion;
    T.NBoson = HAnd.NBoson;
    T.NWord = HAnd.NWord;
    InitOperatorFixedLength(&T, 0, 0);
    RealOperatorAddLaderNormalOrder(&T, oadd, 0, 1, &QDetNotOrderedFcFaBcBa);
    /* set settings for spectra */
    Spectra[0].N = HSet.NE;
    Spectra[0].Emin = HSet.Emin;
    Spectra[0].Emax = HSet.Emax;
    Spectra[0].dE = (HSet.Emax - HSet.Emin) / HSet.NE;
    Spectra[0].EminPole = HSet.Emin;
    Spectra[0].EmaxPole = HSet.Emax;
    Spectra[0].Gamma = HSet.Gamma;
    for (iQ = 0; iQ < HSet.NQ; iQ++)
    {
      printf("Start of core PES itteration for Q=%22.15E (%i/%i)\n", HSet.Q[iQ], iQ, HSet.NQ);
      fflush(stdout);
      /* init spectraTri */
      SpectraTri.N = Max(HSet.NTriDiag, HSet.NGCluster);
      InitTriDiagonalMatrix(&SpectraTri);
      /* shift Onsite energy of Anderson operator */
      oadd[1] = CrAnMask + 0;
      oadd[0] = 0;
      RealOperatorAddLaderNormalOrder(&HAnd, oadd, 2, -HSet.Q[iQ], &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + 1;
      oadd[0] = 1;
      RealOperatorAddLaderNormalOrder(&HAnd, oadd, 2, -HSet.Q[iQ], &QDetNotOrderedFcFaBcBa);
      /* calculate cPES */
      OperatorPsi(T, psiGrd, &Tpsi);
      LanczosTriDiagonalizeDysonLowMem(HAnd, Tpsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation,
                                       gPerturbation, 4, &SpectraTri);
      /* transform tridiagonal to spectra */
      SpectrumAllocate(Spectra[0], &Spectra[iQ]);
      TriDiagonalMatrixToSpectrum(Egrd, SpectraTri, &Spectra[iQ]);
      /* shift Onsite energy of Anderson operator back */
      oadd[1] = CrAnMask + 0;
      oadd[0] = 0;
      RealOperatorAddLaderNormalOrder(&HAnd, oadd, 2, HSet.Q[iQ], &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + 1;
      oadd[0] = 1;
      RealOperatorAddLaderNormalOrder(&HAnd, oadd, 2, HSet.Q[iQ], &QDetNotOrderedFcFaBcBa);
      /* free variables */
      FreeTriDiagonalMatrix(&SpectraTri);
    }
    /* write spectra to disk */
    sprintf(SpecFileName, "%s%03i_cPES", HSet.SpecRootFileName, iOuter + 1);
    FilePrintSpectra(Spectra, HSet.NQ, SpecFileName, "w");
    /* free variables */
    for (iQ = 0; iQ < HSet.NQ; iQ++)
    {
      SpectrumFree(&(Spectra[iQ]));
    }
  }
  /* 128 - XAS - d0 */
  if (HSet.Spectra & 128)
  {
    printf("\nStart of DMFT Spec XAS d0\n\n");
    fflush(stdout);
    Spectra[0].N = HSet.NE;
    Spectra[0].Emin = HSet.Emin;
    Spectra[0].Emax = HSet.Emax;
    Spectra[0].dE = (HSet.Emax - HSet.Emin) / HSet.NE;
    Spectra[0].EminPole = HSet.Emin;
    Spectra[0].EmaxPole = HSet.Emax;
    Spectra[0].Gamma = HSet.Gamma;
    for (iQ = 0; iQ < HSet.NQ; iQ++)
    {
      printf("Start of XAS itteration for Q=%22.15E (%i/%i)\n", HSet.Q[iQ], iQ, HSet.NQ);
      fflush(stdout);
      /* init spectraTri */
      SpectraTri.N = Max(HSet.NTriDiag, HSet.NGCluster);
      InitTriDiagonalMatrix(&SpectraTri);
      /* shift Onsite energy of Anderson operator */
      oadd[1] = CrAnMask + 0;
      oadd[0] = 0;
      RealOperatorAddLaderNormalOrder(&HAnd, oadd, 2, -HSet.Q[iQ], &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + 1;
      oadd[0] = 1;
      RealOperatorAddLaderNormalOrder(&HAnd, oadd, 2, -HSet.Q[iQ], &QDetNotOrderedFcFaBcBa);
      /* calculate XAS spin up */
      InitOperator(&T);
      T.QComplex = 0;
      T.NFermion = HAnd.NFermion;
      T.NBoson = HAnd.NBoson;
      T.NWord = HAnd.NWord;
      InitOperatorFixedLength(&T, 2, 2);
      InitOperatorFixedLength(&T, 4, 2);
      InitOperatorFixedLength(&T, 6, 2);
      OperatorMultiply2NormalOrder(&T, OCrup, Projectd0, 'N', 'N', &QDetNotOrderedFcFaBcBa);
      OperatorPsi(T, psiGrd, &Tpsi);
      FreeOperator(&T);
      LanczosTriDiagonalizeDysonLowMem(HAnd, Tpsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation,
                                       gPerturbation, 4, &SpectraTri);
      /* transform tridiagonal to spectra */
      SpectrumAllocate(Spectra[0], &Spectra[iQ]);
      TriDiagonalMatrixToSpectrum(Egrd, SpectraTri, &Spectra[iQ]);
      /* calculate XAS spin dn */
      InitOperator(&T);
      T.QComplex = 0;
      T.NFermion = HAnd.NFermion;
      T.NBoson = HAnd.NBoson;
      T.NWord = HAnd.NWord;
      InitOperatorFixedLength(&T, 2, 2);
      InitOperatorFixedLength(&T, 4, 2);
      InitOperatorFixedLength(&T, 6, 2);
      OperatorMultiply2NormalOrder(&T, OCrdn, Projectd0, 'N', 'N', &QDetNotOrderedFcFaBcBa);
      OperatorPsi(T, psiGrd, &Tpsi);
      FreeOperator(&T);
      LanczosTriDiagonalizeDysonLowMem(HAnd, Tpsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation,
                                       gPerturbation, 4, &SpectraTri);
      /* transform tridiagonal to spectra */
      SpectrumAllocate(Spectra[0], &SpectraTemp);
      TriDiagonalMatrixToSpectrum(Egrd, SpectraTri, &SpectraTemp);
      /* sum spectra */
      for (i = 0; i < Spectra[0].N; i++)
      {
        Spectra[iQ].Re[i] += SpectraTemp.Re[i];
        Spectra[iQ].Im[i] += SpectraTemp.Im[i];
      }
      /* shift Onsite energy of Anderson operator back */
      oadd[1] = CrAnMask + 0;
      oadd[0] = 0;
      RealOperatorAddLaderNormalOrder(&HAnd, oadd, 2, HSet.Q[iQ], &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + 1;
      oadd[0] = 1;
      RealOperatorAddLaderNormalOrder(&HAnd, oadd, 2, HSet.Q[iQ], &QDetNotOrderedFcFaBcBa);
      /* free variables */
      SpectrumFree(&SpectraTemp);
      FreeTriDiagonalMatrix(&SpectraTri);
    }
    /* write spectra to disk */
    sprintf(SpecFileName, "%s%03i_XAS_d0", HSet.SpecRootFileName, iOuter + 1);
    FilePrintSpectra(Spectra, HSet.NQ, SpecFileName, "w");
    /* free variables */
    for (iQ = 0; iQ < HSet.NQ; iQ++)
    {
      SpectrumFree(&(Spectra[iQ]));
    }
  }
  /* 256 - XAS - d1 */
  if (HSet.Spectra & 256)
  {
    printf("\nStart of DMFT Spec XAS d1\n\n");
    fflush(stdout);
    Spectra[0].N = HSet.NE;
    Spectra[0].Emin = HSet.Emin;
    Spectra[0].Emax = HSet.Emax;
    Spectra[0].dE = (HSet.Emax - HSet.Emin) / HSet.NE;
    Spectra[0].EminPole = HSet.Emin;
    Spectra[0].EmaxPole = HSet.Emax;
    Spectra[0].Gamma = HSet.Gamma;
    for (iQ = 0; iQ < HSet.NQ; iQ++)
    {
      printf("Start of XAS itteration for Q=%22.15E (%i/%i)\n", HSet.Q[iQ], iQ, HSet.NQ);
      fflush(stdout);
      /* init spectraTri */
      SpectraTri.N = Max(HSet.NTriDiag, HSet.NGCluster);
      InitTriDiagonalMatrix(&SpectraTri);
      /* shift Onsite energy of Anderson operator */
      oadd[1] = CrAnMask + 0;
      oadd[0] = 0;
      RealOperatorAddLaderNormalOrder(&HAnd, oadd, 2, -HSet.Q[iQ], &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + 1;
      oadd[0] = 1;
      RealOperatorAddLaderNormalOrder(&HAnd, oadd, 2, -HSet.Q[iQ], &QDetNotOrderedFcFaBcBa);
      /* calculate XAS spin up */
      InitOperator(&T);
      T.QComplex = 0;
      T.NFermion = HAnd.NFermion;
      T.NBoson = HAnd.NBoson;
      T.NWord = HAnd.NWord;
      InitOperatorFixedLength(&T, 2, 2);
      InitOperatorFixedLength(&T, 4, 2);
      InitOperatorFixedLength(&T, 6, 2);
      OperatorMultiply2NormalOrder(&T, OCrup, Projectd1, 'N', 'N', &QDetNotOrderedFcFaBcBa);
      OperatorPsi(T, psiGrd, &Tpsi);
      FreeOperator(&T);
      LanczosTriDiagonalizeDysonLowMem(HAnd, Tpsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation,
                                       gPerturbation, 4, &SpectraTri);
      /* transform tridiagonal to spectra */
      SpectrumAllocate(Spectra[0], &Spectra[iQ]);
      TriDiagonalMatrixToSpectrum(Egrd, SpectraTri, &Spectra[iQ]);
      /* calculate XAS spin dn */
      InitOperator(&T);
      T.QComplex = 0;
      T.NFermion = HAnd.NFermion;
      T.NBoson = HAnd.NBoson;
      T.NWord = HAnd.NWord;
      InitOperatorFixedLength(&T, 2, 2);
      InitOperatorFixedLength(&T, 4, 2);
      InitOperatorFixedLength(&T, 6, 2);
      OperatorMultiply2NormalOrder(&T, OCrdn, Projectd1, 'N', 'N', &QDetNotOrderedFcFaBcBa);
      OperatorPsi(T, psiGrd, &Tpsi);
      FreeOperator(&T);
      LanczosTriDiagonalizeDysonLowMem(HAnd, Tpsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation,
                                       gPerturbation, 4, &SpectraTri);
      /* transform tridiagonal to spectra */
      SpectrumAllocate(Spectra[0], &SpectraTemp);
      TriDiagonalMatrixToSpectrum(Egrd, SpectraTri, &SpectraTemp);
      /* sum spectra */
      for (i = 0; i < Spectra[0].N; i++)
      {
        Spectra[iQ].Re[i] += SpectraTemp.Re[i];
        Spectra[iQ].Im[i] += SpectraTemp.Im[i];
      }
      /* shift Onsite energy of Anderson operator back */
      oadd[1] = CrAnMask + 0;
      oadd[0] = 0;
      RealOperatorAddLaderNormalOrder(&HAnd, oadd, 2, HSet.Q[iQ], &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + 1;
      oadd[0] = 1;
      RealOperatorAddLaderNormalOrder(&HAnd, oadd, 2, HSet.Q[iQ], &QDetNotOrderedFcFaBcBa);
      /* free variables */
      SpectrumFree(&SpectraTemp);
      FreeTriDiagonalMatrix(&SpectraTri);
    }
    /* write spectra to disk */
    sprintf(SpecFileName, "%s%03i_XAS_d1", HSet.SpecRootFileName, iOuter + 1);
    FilePrintSpectra(Spectra, HSet.NQ, SpecFileName, "w");
    /* free variables */
    for (iQ = 0; iQ < HSet.NQ; iQ++)
    {
      SpectrumFree(&(Spectra[iQ]));
    }
  }
  /* 512 - core PES - d0 */
  if (HSet.Spectra & 512)
  {
    printf("\nStart of DMFT Spec core PES\n\n");
    fflush(stdout);
    /* set settings for spectra */
    Spectra[0].N = HSet.NE;
    Spectra[0].Emin = HSet.Emin;
    Spectra[0].Emax = HSet.Emax;
    Spectra[0].dE = (HSet.Emax - HSet.Emin) / HSet.NE;
    Spectra[0].EminPole = HSet.Emin;
    Spectra[0].EmaxPole = HSet.Emax;
    Spectra[0].Gamma = HSet.Gamma;
    for (iQ = 0; iQ < HSet.NQ; iQ++)
    {
      printf("Start of core PES itteration for Q=%22.15E (%i/%i)\n", HSet.Q[iQ], iQ, HSet.NQ);
      fflush(stdout);
      /* init spectraTri */
      SpectraTri.N = Max(HSet.NTriDiag, HSet.NGCluster);
      InitTriDiagonalMatrix(&SpectraTri);
      /* shift Onsite energy of Anderson operator */
      oadd[1] = CrAnMask + 0;
      oadd[0] = 0;
      RealOperatorAddLaderNormalOrder(&HAnd, oadd, 2, -HSet.Q[iQ], &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + 1;
      oadd[0] = 1;
      RealOperatorAddLaderNormalOrder(&HAnd, oadd, 2, -HSet.Q[iQ], &QDetNotOrderedFcFaBcBa);
      /* calculate cPES */
      OperatorPsi(Projectd0, psiGrd, &Tpsi);
      LanczosTriDiagonalizeDysonLowMem(HAnd, Tpsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation,
                                       gPerturbation, 4, &SpectraTri);
      /* transform tridiagonal to spectra */
      SpectrumAllocate(Spectra[0], &Spectra[iQ]);
      TriDiagonalMatrixToSpectrum(Egrd, SpectraTri, &Spectra[iQ]);
      /* shift Onsite energy of Anderson operator back */
      oadd[1] = CrAnMask + 0;
      oadd[0] = 0;
      RealOperatorAddLaderNormalOrder(&HAnd, oadd, 2, HSet.Q[iQ], &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + 1;
      oadd[0] = 1;
      RealOperatorAddLaderNormalOrder(&HAnd, oadd, 2, HSet.Q[iQ], &QDetNotOrderedFcFaBcBa);
      /* free variables */
      FreeTriDiagonalMatrix(&SpectraTri);
    }
    /* write spectra to disk */
    if (iOuter == 0)
    {
      sprintf(SpecFileName, "%s001_cPES_d0", HSet.SpecRootFileName);
    }
    else if (iOuter < 9)
    {
      sprintf(SpecFileName, "%s00%1i_cPES_d0", HSet.SpecRootFileName, iOuter + 1);
    }
    else if (iOuter < 99)
    {
      sprintf(SpecFileName, "%s0%2i_cPES_d0", HSet.SpecRootFileName, iOuter + 1);
    }
    else if (iOuter < 999)
    {
      sprintf(SpecFileName, "%s%3i_cPES_d0", HSet.SpecRootFileName, iOuter + 1);
    }
    FilePrintSpectra(Spectra, HSet.NQ, SpecFileName, "w");
    /* free variables */
    for (iQ = 0; iQ < HSet.NQ; iQ++)
    {
      SpectrumFree(&(Spectra[iQ]));
    }
  }
  /* 1024 - core PES - d1 */
  if (HSet.Spectra & 1024)
  {
    printf("\nStart of DMFT Spec core PES\n\n");
    fflush(stdout);
    /* set settings for spectra */
    Spectra[0].N = HSet.NE;
    Spectra[0].Emin = HSet.Emin;
    Spectra[0].Emax = HSet.Emax;
    Spectra[0].dE = (HSet.Emax - HSet.Emin) / HSet.NE;
    Spectra[0].EminPole = HSet.Emin;
    Spectra[0].EmaxPole = HSet.Emax;
    Spectra[0].Gamma = HSet.Gamma;
    for (iQ = 0; iQ < HSet.NQ; iQ++)
    {
      printf("Start of core PES itteration for Q=%22.15E (%i/%i)\n", HSet.Q[iQ], iQ, HSet.NQ);
      fflush(stdout);
      /* init spectraTri */
      SpectraTri.N = Max(HSet.NTriDiag, HSet.NGCluster);
      InitTriDiagonalMatrix(&SpectraTri);
      /* shift Onsite energy of Anderson operator */
      oadd[1] = CrAnMask + 0;
      oadd[0] = 0;
      RealOperatorAddLaderNormalOrder(&HAnd, oadd, 2, -HSet.Q[iQ], &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + 1;
      oadd[0] = 1;
      RealOperatorAddLaderNormalOrder(&HAnd, oadd, 2, -HSet.Q[iQ], &QDetNotOrderedFcFaBcBa);
      /* calculate cPES */
      OperatorPsi(Projectd1, psiGrd, &Tpsi);
      LanczosTriDiagonalizeDysonLowMem(HAnd, Tpsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation,
                                       gPerturbation, 4, &SpectraTri);
      /* transform tridiagonal to spectra */
      SpectrumAllocate(Spectra[0], &Spectra[iQ]);
      TriDiagonalMatrixToSpectrum(Egrd, SpectraTri, &Spectra[iQ]);
      /* shift Onsite energy of Anderson operator back */
      oadd[1] = CrAnMask + 0;
      oadd[0] = 0;
      RealOperatorAddLaderNormalOrder(&HAnd, oadd, 2, HSet.Q[iQ], &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + 1;
      oadd[0] = 1;
      RealOperatorAddLaderNormalOrder(&HAnd, oadd, 2, HSet.Q[iQ], &QDetNotOrderedFcFaBcBa);
      /* free variables */
      FreeTriDiagonalMatrix(&SpectraTri);
    }
    /* write spectra to disk */
    if (iOuter == 0)
    {
      sprintf(SpecFileName, "%s001_cPES_d1", HSet.SpecRootFileName);
    }
    else if (iOuter < 9)
    {
      sprintf(SpecFileName, "%s00%1i_cPES_d1", HSet.SpecRootFileName, iOuter + 1);
    }
    else if (iOuter < 99)
    {
      sprintf(SpecFileName, "%s0%2i_cPES_d1", HSet.SpecRootFileName, iOuter + 1);
    }
    else if (iOuter < 999)
    {
      sprintf(SpecFileName, "%s%3i_cPES_d1", HSet.SpecRootFileName, iOuter + 1);
    }
    FilePrintSpectra(Spectra, HSet.NQ, SpecFileName, "w");
    /* free variables */
    for (iQ = 0; iQ < HSet.NQ; iQ++)
    {
      SpectrumFree(&(Spectra[iQ]));
    }
  }
  /* 2048 - core PES - d2 */
  if (HSet.Spectra & 2048)
  {
    printf("\nStart of DMFT Spec core PES\n\n");
    fflush(stdout);
    /* set settings for spectra */
    Spectra[0].N = HSet.NE;
    Spectra[0].Emin = HSet.Emin;
    Spectra[0].Emax = HSet.Emax;
    Spectra[0].dE = (HSet.Emax - HSet.Emin) / HSet.NE;
    Spectra[0].EminPole = HSet.Emin;
    Spectra[0].EmaxPole = HSet.Emax;
    Spectra[0].Gamma = HSet.Gamma;
    for (iQ = 0; iQ < HSet.NQ; iQ++)
    {
      printf("Start of core PES itteration for Q=%22.15E (%i/%i)\n", HSet.Q[iQ], iQ, HSet.NQ);
      fflush(stdout);
      /* init spectraTri */
      SpectraTri.N = Max(HSet.NTriDiag, HSet.NGCluster);
      InitTriDiagonalMatrix(&SpectraTri);
      /* shift Onsite energy of Anderson operator */
      oadd[1] = CrAnMask + 0;
      oadd[0] = 0;
      RealOperatorAddLaderNormalOrder(&HAnd, oadd, 2, -HSet.Q[iQ], &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + 1;
      oadd[0] = 1;
      RealOperatorAddLaderNormalOrder(&HAnd, oadd, 2, -HSet.Q[iQ], &QDetNotOrderedFcFaBcBa);
      /* calculate cPES */
      OperatorPsi(Projectd2, psiGrd, &Tpsi);
      LanczosTriDiagonalizeDysonLowMem(HAnd, Tpsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation,
                                       gPerturbation, 4, &SpectraTri);
      /* transform tridiagonal to spectra */
      SpectrumAllocate(Spectra[0], &Spectra[iQ]);
      TriDiagonalMatrixToSpectrum(Egrd, SpectraTri, &Spectra[iQ]);
      /* shift Onsite energy of Anderson operator back */
      oadd[1] = CrAnMask + 0;
      oadd[0] = 0;
      RealOperatorAddLaderNormalOrder(&HAnd, oadd, 2, HSet.Q[iQ], &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + 1;
      oadd[0] = 1;
      RealOperatorAddLaderNormalOrder(&HAnd, oadd, 2, HSet.Q[iQ], &QDetNotOrderedFcFaBcBa);
      /* free variables */
      FreeTriDiagonalMatrix(&SpectraTri);
    }
    /* write spectra to disk */
    if (iOuter == 0)
    {
      sprintf(SpecFileName, "%s001_cPES_d2", HSet.SpecRootFileName);
    }
    else if (iOuter < 9)
    {
      sprintf(SpecFileName, "%s00%1i_cPES_d2", HSet.SpecRootFileName, iOuter + 1);
    }
    else if (iOuter < 99)
    {
      sprintf(SpecFileName, "%s0%2i_cPES_d2", HSet.SpecRootFileName, iOuter + 1);
    }
    else if (iOuter < 999)
    {
      sprintf(SpecFileName, "%s%3i_cPES_d2", HSet.SpecRootFileName, iOuter + 1);
    }
    FilePrintSpectra(Spectra, HSet.NQ, SpecFileName, "w");
    /* free variables */
    for (iQ = 0; iQ < HSet.NQ; iQ++)
    {
      SpectrumFree(&(Spectra[iQ]));
    }
  }





  return 0;
}
