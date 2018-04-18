#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "../Core/StructDef.c"
#include "../header.h"
#include "../Core/TimeStamp.c"
#include "../Core/Hash.c"
#include "../Core/BasicMath.c"
#include "../IO/Print.c"
#include "../IO/Read.c"
#include "../Core/Error.c"
#include "../Core/WaveFunctionArithmetic.c"
#include "../Core/ComplexWaveFunctionArithmetic.c"
#include "../Core/OperatorArithmetic.c"
#include "../Core/OperatorWaveFunction.c"
#include "../Core/Lanczos.c"
#include "../Core/Spectral.c"
#include "../Core/DensityMatrix.c"
#include "../Model/StandardOperators.c"
#include "../Hubbard_Model/CellFunctions.c"


/*  Creates an Dimerized Spin Chain Hamiltonian
 *  Exchange between site 2 i     and 2 i + 1 is JA
 *  Exchange between site 2 i + 1 and 2 i + 2 is JB ( Between site 2 i - 1 and 1 is JB )
 *
 *  Si . Sj = Szi Szj + 1/2 ( Si+ Sj- + Si- Sj+ )
 *  Sz = 1/2 (a+_{up} a_{up} - a+_{dn} a_{dn})
 *  S+ = a+_{up} a_{dn}
 *  S- = a+_{dn} a_{up}
 */
uint32_t CreateDimerizedSpinChainHamiltonian(DimerizedSpinChainSettings Set, OperatorType * H)
{
  uint32_t i, j;
  uint16_t oadd[4];
  InitOperator(H);
  H->NFermion = 4 * Set.Ndimer;
  H->NWord = Min(H->NFermion, 1 + (H->NFermion - 1) / 8);
  InitOperatorFixedLength(H, 2, 8);
  InitOperatorFixedLength(H, 4, 8);
  for (i = 0; i < 2 * Set.Ndimer; i += 2)
  {                             /* add S2i . S2i+1 to the Hamiltonian with JA */
    j = i + 1;
    if (j >= 2 * Set.Ndimer)
      j -= 2 * Set.Ndimer;
    oadd[0] = CrAnMask + 2 * i + 1;
    oadd[1] = 2 * i + 1;
    oadd[2] = CrAnMask + 2 * j + 1;
    oadd[3] = 2 * j + 1;
    RealOperatorAddLaderNormalOrder(H, oadd, 4, 0.25 * Set.JA, &QDetNotOrderedFcFaBcBa);
    oadd[0] = CrAnMask + 2 * i;
    oadd[1] = 2 * i;
    oadd[2] = CrAnMask + 2 * j;
    oadd[3] = 2 * j;
    RealOperatorAddLaderNormalOrder(H, oadd, 4, 0.25 * Set.JA, &QDetNotOrderedFcFaBcBa);
    oadd[0] = CrAnMask + 2 * i;
    oadd[1] = 2 * i;
    oadd[2] = CrAnMask + 2 * j + 1;
    oadd[3] = 2 * j + 1;
    RealOperatorAddLaderNormalOrder(H, oadd, 4, -0.25 * Set.JA, &QDetNotOrderedFcFaBcBa);
    oadd[0] = CrAnMask + 2 * i + 1;
    oadd[1] = 2 * i + 1;
    oadd[2] = CrAnMask + 2 * j;
    oadd[3] = 2 * j;
    RealOperatorAddLaderNormalOrder(H, oadd, 4, -0.25 * Set.JA, &QDetNotOrderedFcFaBcBa);
    oadd[0] = CrAnMask + 2 * i + 1;
    oadd[1] = 2 * i;
    oadd[2] = CrAnMask + 2 * j;
    oadd[3] = 2 * j + 1;
    RealOperatorAddLaderNormalOrder(H, oadd, 4, 0.5 * Set.JA, &QDetNotOrderedFcFaBcBa);
    oadd[0] = CrAnMask + 2 * i;
    oadd[1] = 2 * i + 1;
    oadd[2] = CrAnMask + 2 * j + 1;
    oadd[3] = 2 * j;
    RealOperatorAddLaderNormalOrder(H, oadd, 4, 0.5 * Set.JA, &QDetNotOrderedFcFaBcBa);
  }
  for (i = 1; i < 2 * Set.Ndimer; i += 2)
  {                             /* add S2i+1 . S2i+2 to the Hamiltonian with JB */
    j = i + 1;
    if (j >= 2 * Set.Ndimer)
      j -= 2 * Set.Ndimer;
    oadd[0] = CrAnMask + 2 * i + 1;
    oadd[1] = 2 * i + 1;
    oadd[2] = CrAnMask + 2 * j + 1;
    oadd[3] = 2 * j + 1;
    RealOperatorAddLaderNormalOrder(H, oadd, 4, 0.25 * Set.JB, &QDetNotOrderedFcFaBcBa);
    oadd[0] = CrAnMask + 2 * i;
    oadd[1] = 2 * i;
    oadd[2] = CrAnMask + 2 * j;
    oadd[3] = 2 * j;
    RealOperatorAddLaderNormalOrder(H, oadd, 4, 0.25 * Set.JB, &QDetNotOrderedFcFaBcBa);
    oadd[0] = CrAnMask + 2 * i;
    oadd[1] = 2 * i;
    oadd[2] = CrAnMask + 2 * j + 1;
    oadd[3] = 2 * j + 1;
    RealOperatorAddLaderNormalOrder(H, oadd, 4, -0.25 * Set.JB, &QDetNotOrderedFcFaBcBa);
    oadd[0] = CrAnMask + 2 * i + 1;
    oadd[1] = 2 * i + 1;
    oadd[2] = CrAnMask + 2 * j;
    oadd[3] = 2 * j;
    RealOperatorAddLaderNormalOrder(H, oadd, 4, -0.25 * Set.JB, &QDetNotOrderedFcFaBcBa);
    oadd[0] = CrAnMask + 2 * i + 1;
    oadd[1] = 2 * i;
    oadd[2] = CrAnMask + 2 * j;
    oadd[3] = 2 * j + 1;
    RealOperatorAddLaderNormalOrder(H, oadd, 4, 0.5 * Set.JB, &QDetNotOrderedFcFaBcBa);
    oadd[0] = CrAnMask + 2 * i;
    oadd[1] = 2 * i + 1;
    oadd[2] = CrAnMask + 2 * j + 1;
    oadd[3] = 2 * j;
    RealOperatorAddLaderNormalOrder(H, oadd, 4, 0.5 * Set.JB, &QDetNotOrderedFcFaBcBa);
  }
  OperatorRemoveSmallValues(H, DBL_MIN * 100);
  return 0;
}


uint32_t CreateDimerizedSpinChainStartWaveFunctionOperator(DimerizedSpinChainSettings Set, OperatorType * O)
{
  uint16_t oadd[2 * Set.Ndimer];
  uint32_t i;
  InitOperator(O);
  O->NFermion = 4 * Set.Ndimer;
  O->NWord = Min(O->NFermion, 1 + (O->NFermion - 1) / 8);
  InitOperatorFixedLength(O, 2 * Set.Ndimer, 1);
  for (i = 0; i < 2 * Set.Ndimer; i++)
  {
    oadd[i] = CrAnMask + (uint16_t) (2 * i + (i) % 2);
  }
  RealOperatorAddLaderNormalOrder(O, oadd, 2 * Set.Ndimer, 1.0, &QDetNotOrderedFcFaBcBa);
  return 0;
}

uint32_t CreateDimerizedSpinChainStartWaveFunction(DimerizedSpinChainSettings Set, WaveFunctionType * Psi)
{
  OperatorType O;
  WaveFunctionType f;
  uint8_t *Determinant;
  CreateDimerizedSpinChainStartWaveFunctionOperator(Set, &O);
  f.QComplex = 0;
  f.NFermion = 4 * Set.Ndimer;
  f.NBoson = 0;
  f.NBitsKey = Set.NBitsKey;
  if(WaveFunctionInit(&f))
  {
    printf("WaveFunctionInit failed in CreateDimerizedSpinChainStartWaveFunction\n");
    return 1;
  }
  Determinant = calloc(f.NWord, sizeof(uint8_t));
  WaveFunctionAddElement(&f, Determinant, 1);
  OperatorPsi(O, f, Psi);
  FreeOperator(&O);
  return 0;
}


uint32_t CreateDimerizedSpinChainSzOperatorMomentumSpace(DimerizedSpinChainSettings Set, OperatorType * O)
{
  uint16_t i, k;
  uint16_t oadd[2];
  double unitcellration = 2.0 * Pi / (2 * Set.Ndimer);
  for (i = 0; i < 2 * Set.Ndimer; i++)
  {
    InitOperator(&O[i]);
    O[i].NFermion = 4 * Set.Ndimer;
    O[i].NWord = Min(O[i].NFermion, 1 + (O[i].NFermion - 1) / 8);
    InitOperatorFixedLength(&O[i], 2, 8);
    sprintf(O[i].Name, "Operator Sz at k: %i.", i);
  }
  for (k = 0; k < 2 * Set.Ndimer; k++)
  {
    for (i = 0; i < 2 * Set.Ndimer; i++)
    {
      oadd[0] = CrAnMask + 1 + 2 * i;
      oadd[1] = 1 + 2 * i;
      ComplexOperatorAddLaderNormalOrder(&O[k], oadd, 2,
                                         cos(unitcellration * k * i) * sqrt(0.125 / ((double)Set.Ndimer)),
                                         sin(unitcellration * k * i) * sqrt(0.125 / ((double)Set.Ndimer)),
                                         &QDetNotOrderedFcFaBcBa);
      oadd[0] = CrAnMask + 2 * i;
      oadd[1] = 2 * i;
      ComplexOperatorAddLaderNormalOrder(&O[k], oadd, 2,
                                         -cos(unitcellration * k * i) * sqrt(0.125 / ((double)Set.Ndimer)),
                                         -sin(unitcellration * k * i) * sqrt(0.125 / ((double)Set.Ndimer)),
                                         &QDetNotOrderedFcFaBcBa);
    }
    OperatorRemoveSmallValues(&O[k], DBL_MIN * 100);
  }
  return 0;
}

int main(int argc, const char *argv[])
{
  argc -= 1;
  double TotEmax, *EnergyList, *WaveFunctionQuality;
  char   SpecFileName[255];
  double Hist[20];
  uint32_t jpsi;
  uint32_t i, j;
  MatrixType Rho, RhoCorr;
  MatrixType RhoTot, RhoCorrTot;
  WaveFunctionType *PsiList;
  DimerizedSpinChainSettings Set;
  WaveFunctionType OPsi;
  OperatorType H;
  TimeType Time;
  TriDiagonalMatrixType s;
  OperatorType *T;
  SpectraType *Spec;

  InitializeTiming(&Time);
  BeginTiming("Total", &Time);
  TimeStampAll();
  FileReadDimerizedSpinChainSettings(&Set, argv[argc]);
  PrintDimerizedSpinChainSettings(Set);
  BeginTiming("Create Hamiltonian", &Time);
  CreateDimerizedSpinChainHamiltonian(Set, &H);
  EndTiming("Create Hamiltonian", &Time);
  PrintOperatorReadable(H);
  BeginTiming("Create Transition Matrix", &Time);
  s.N = Set.NTriDiag;
  InitTriDiagonalMatrix(&s);
  T = malloc(Set.Ndimer * 2 * sizeof(OperatorType));
  Spec = malloc(Set.Ndimer * 2 * sizeof(SpectraType));
  Spec[0].N = Set.NE;
  Spec[0].Emin = Set.Emin;
  Spec[0].Emax = Set.Emax;
  Spec[0].dE = (Set.Emax - Set.Emin) / Set.NE;
  Spec[0].EminPole = Set.Emin;
  Spec[0].EmaxPole = Set.Emax;
  Spec[0].Gamma = Set.Gamma;
  for (i = 0; i < Set.Ndimer * 2; i++)
  {
    SpectrumAllocate(Spec[0], &Spec[i]);
  }
  CreateDimerizedSpinChainSzOperatorMomentumSpace(Set, T);
  EndTiming("Create Transition Matrix", &Time);
  for (i = 0; i < Set.Ndimer * 2; i++)
  {
    PrintOperatorReadable(T[i]);
  }
  PsiList = malloc(Set.NPsi * sizeof(WaveFunctionType));
  EnergyList = malloc(Set.NPsi * sizeof(double));
  WaveFunctionQuality = malloc(Set.NPsi * sizeof(double));
  for (jpsi = 0; jpsi < Set.NPsi; jpsi++)
  {
    CreateDimerizedSpinChainStartWaveFunction(Set, &PsiList[jpsi]);
    if (jpsi != 0)
      WaveFunctionOrthogonalize(&PsiList[jpsi], PsiList, jpsi);
    WaveFunctionRandomize(&PsiList[jpsi]);
    PrintWaveFunctionReadableTruncated(PsiList[jpsi], 0.1);
    TotEmax = 0.75 * Set.Ndimer * (fabs(Set.JA) + fabs(Set.JB));
    BeginTiming("LanczosGroundState", &Time);
    LanczosGroundState(H, &PsiList[jpsi], &Set.zeroLancWavefunction, Set.zeroLancEnergy, TotEmax, PsiList, jpsi);
    EndTiming("LanczosGroundState", &Time);
    TotEmax *= 2;
    BeginTiming("LanczosGroundState", &Time);
    LanczosGroundState(H, &PsiList[jpsi], &Set.zeroLancWavefunction, Set.zeroLancEnergy, TotEmax, PsiList, jpsi);
    EndTiming("LanczosGroundState", &Time);
    sprintf(PsiList[jpsi].Name, "Eigenstate #%i", jpsi);
    PsiHermitianOperatorPsi(H, PsiList[jpsi], &EnergyList[jpsi]);
    PsiHermitianOperatorSqrPsi(H, PsiList[jpsi], &WaveFunctionQuality[jpsi]);
    WaveFunctionQuality[jpsi] /= Sqr(EnergyList[jpsi]);
    PrintWaveFunctionReadableTruncated(PsiList[jpsi], 0.1);
    CalculatePrefactorHistogram(Hist, PsiList[jpsi], 10.0, 20);
    printf("\nHistogram of prefactors in Wave-function\n");
    for (i = 0; i < 19; i++)
    {
      printf("10^(-%2i) to 10^(-%2i) | %22.15E\n", i, i + 1, Hist[i]);
    }
    printf("10^(-%2i) to 0        | %22.15E\n", i, Hist[19]);
    printf("\nOrthonormality test of wave-functions\n");
    for (i = 0; i <= jpsi; i++)
    {
      for (j = 0; j <= jpsi; j++)
      {
        printf(" %22.15E", WaveFunctionDotProduct(PsiList[i], PsiList[j]));
      }
      printf("\n");
    }
    printf("\n");
    printf("State#   Energy                  <psi|H^2|psi>/(<psi|H|psi>^2)\n");
    for (i = 0; i <= jpsi; i++)
    {
      printf("%2i      %22.15E  %22.15E\n", i, EnergyList[i], WaveFunctionQuality[i]);
    }
    CalculateDensityMatrix(&Rho, PsiList[jpsi]);
    PrintMatrix(Rho);
    CalculateDensityCorrelationMatrix(&RhoCorr, PsiList[jpsi]);
    PrintMatrix(RhoCorr);
    if ((EnergyList[jpsi] - EnergyList[0]) > 1E-4)
    {
      printf("Average of Density matrices for state 0 to %i\n", jpsi - 1);
      MatrixScalarMultiply(&RhoTot, 1.0 / ((double)(jpsi)), 0);
      MatrixScalarMultiply(&RhoCorrTot, 1.0 / ((double)(jpsi)), 0);
      PrintMatrix(RhoTot);
      PrintMatrix(RhoCorrTot);
      EndTiming("Total", &Time);
      PrintTiming(Time);
      return 0;
    }
    if (jpsi == 0)
    {
      MatrixCopy(&RhoTot, Rho);
      MatrixCopy(&RhoCorrTot, RhoCorr);
    }
    else
    {
      MatrixAdd(&RhoTot, Rho, 1, 0);
      MatrixAdd(&RhoCorrTot, RhoCorr, 1, 0);
    }
    FreeMatrix(&Rho);
    FreeMatrix(&RhoCorr);
    for (i = 0; i < Set.Ndimer * 2; i++)
    {
      BeginTiming("Lanczos Tri Diagonalization", &Time);
      OperatorPsi(T[i], PsiList[jpsi], &OPsi);
      LanczosTriDiagonalize(H, OPsi, Set.zeroLancWavefunction, &s, Set.Gamma);
      TriDiagonalMatrixToSpectrum(EnergyList[jpsi], s, &Spec[i]);
      EndTiming("Lanczos Tri Diagonalization", &Time);
    }
    if (jpsi < 10)
      sprintf(SpecFileName, "%s00%1i", Set.SpecFile, jpsi);
    else if (jpsi < 100)
      sprintf(SpecFileName, "%s0%2i", Set.SpecFile, jpsi);
    else if (jpsi < 1000)
      sprintf(SpecFileName, "%s%3i", Set.SpecFile, jpsi);
    else
      sprintf(SpecFileName, "%s%i", Set.SpecFile, jpsi);
    FilePrintSpectra(Spec, 2 * Set.Ndimer, SpecFileName, "w");
  }
  EndTiming("Total", &Time);
  PrintTiming(Time);
  return 0;
}
