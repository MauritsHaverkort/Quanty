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

/*  Creates an Heisenberg Hamiltonian
 *  the fermions at position 2 i     correspond to the real-space position i with spin down
 *  the fermions at position 2 i + 1 correspont to the real-space position i wiht spin up
 *
 *  Si . Sj = Szi Szj + 1/2 ( Si+ Sj- + Si- Sj+ )
 *  Sz = 1/2 (a+_{up} a_{up} - a+_{dn} a_{dn})
 *  S+ = a+_{up} a_{dn}
 *  S- = a+_{dn} a_{up}
 */
uint32_t CreateHeisenbergModelHamiltonian(HeisenbergModelSettings Set, OperatorType * H)
{
  uint32_t i, j;
  int32_t C[3];
  uint16_t oadd[4];
  InitOperator(H);
  H->NFermion = 2 * Set.Cluster.NCell;
  H->NWord = Min(H->NFermion, 1 + (H->NFermion - 1) / 8);
  InitOperatorFixedLength(H, 2, 8);
  InitOperatorFixedLength(H, 4, 8);
  for (i = 0; i < Set.Cluster.NCell; i++)
  {
    for (j = 0; j < Set.Cluster.NCell; j++)
    {
      ShortesConnectingVectorPeriodicSuperCell(Set.Cluster.CellPos[i], Set.Cluster.CellPos[j], Set.Cluster, C);
      if (Dot(C, C) == 1)
      {                         /* add Si . Sj to the Hamiltonian */
        oadd[0] = CrAnMask + 2 * i + 1;
        oadd[1] = 2 * i + 1;
        oadd[2] = CrAnMask + 2 * j + 1;
        oadd[3] = 2 * j + 1;
        RealOperatorAddLaderNormalOrder(H, oadd, 4, 0.25 * Set.J, &QDetNotOrderedFcFaBcBa);
        oadd[0] = CrAnMask + 2 * i;
        oadd[1] = 2 * i;
        oadd[2] = CrAnMask + 2 * j;
        oadd[3] = 2 * j;
        RealOperatorAddLaderNormalOrder(H, oadd, 4, 0.25 * Set.J, &QDetNotOrderedFcFaBcBa);
        oadd[0] = CrAnMask + 2 * i;
        oadd[1] = 2 * i;
        oadd[2] = CrAnMask + 2 * j + 1;
        oadd[3] = 2 * j + 1;
        RealOperatorAddLaderNormalOrder(H, oadd, 4, -0.25 * Set.J, &QDetNotOrderedFcFaBcBa);
        oadd[0] = CrAnMask + 2 * i + 1;
        oadd[1] = 2 * i + 1;
        oadd[2] = CrAnMask + 2 * j;
        oadd[3] = 2 * j;
        RealOperatorAddLaderNormalOrder(H, oadd, 4, -0.25 * Set.J, &QDetNotOrderedFcFaBcBa);
        oadd[0] = CrAnMask + 2 * i + 1;
        oadd[1] = 2 * i;
        oadd[2] = CrAnMask + 2 * j;
        oadd[3] = 2 * j + 1;
        RealOperatorAddLaderNormalOrder(H, oadd, 4, 0.5 * Set.J, &QDetNotOrderedFcFaBcBa);
        oadd[0] = CrAnMask + 2 * i;
        oadd[1] = 2 * i + 1;
        oadd[2] = CrAnMask + 2 * j + 1;
        oadd[3] = 2 * j;
        RealOperatorAddLaderNormalOrder(H, oadd, 4, 0.5 * Set.J, &QDetNotOrderedFcFaBcBa);
      }
    }
  }
  return 0;
}

/*  Creates Set.Cluster.NCell operators
 *  O[0] flips the spin at site 0
 *  O[i] flips the spin at site 0 and at site i (sum normalized to sqrt(0.5))
 *  The 'flip' operator is the Sz operator
 */
uint32_t CreateHeisenbergModelSpinFlipOperatorRealSpace(HeisenbergModelSettings Set, OperatorType * O)
{
  uint16_t i;
  uint16_t oadd[2];
  for (i = 0; i < Set.Cluster.NCell; i++)
  {
    InitOperator(&O[i]);
    O[i].NFermion = 2 * Set.Cluster.NCell;
    O[i].NWord = Min(O[i].NFermion, 1 + (O[i].NFermion - 1) / 8);
    InitOperatorFixedLength(&O[i], 2, 8);
    sprintf(O[i].Name, "Operator Sz on site 0 and %i.", i);
  }
  oadd[0] = CrAnMask + 1;
  oadd[1] = 1;
  RealOperatorAddLaderNormalOrder(&O[0], oadd, 2, 0.5, &QDetNotOrderedFcFaBcBa);
  oadd[0] = CrAnMask;
  oadd[1] = 0;
  RealOperatorAddLaderNormalOrder(&O[0], oadd, 2, -0.5, &QDetNotOrderedFcFaBcBa);
  for (i = 1; i < Set.Cluster.NCell; i++)
  {
    oadd[0] = CrAnMask + 1;
    oadd[1] = 1;
    RealOperatorAddLaderNormalOrder(&O[i], oadd, 2, sqrt(0.125), &QDetNotOrderedFcFaBcBa);
    oadd[0] = CrAnMask;
    oadd[1] = 0;
    RealOperatorAddLaderNormalOrder(&O[i], oadd, 2, -sqrt(0.125), &QDetNotOrderedFcFaBcBa);
    oadd[0] = CrAnMask + 1 + 2 * i;
    oadd[1] = 1 + 2 * i;
    RealOperatorAddLaderNormalOrder(&O[i], oadd, 2, sqrt(0.125), &QDetNotOrderedFcFaBcBa);
    oadd[0] = CrAnMask + 2 * i;
    oadd[1] = 2 * i;
    RealOperatorAddLaderNormalOrder(&O[i], oadd, 2, -sqrt(0.125), &QDetNotOrderedFcFaBcBa);
  }
  return 0;
}

/*  Creates Set.Cluster.NCell operators
 *  O[0] flips the spin at site 0
 *  O[i] flips the spin at site 0 and at site i (sum normalized to sqrt(0.5))
 *  The 'flip' operator is the Sx operator
 */
uint32_t CreateHeisenbergModelSpinFlipOperatorSxRealSpace(HeisenbergModelSettings Set, OperatorType * O)
{
  uint16_t i;
  uint16_t oadd[2];
  for (i = 0; i < Set.Cluster.NCell; i++)
  {
    InitOperator(&O[i]);
    O[i].NFermion = 2 * Set.Cluster.NCell;
    O[i].NWord = Min(O[i].NFermion, 1 + (O[i].NFermion - 1) / 8);
    InitOperatorFixedLength(&O[i], 2, 8);
    sprintf(O[i].Name, "Operator Sx on site 0 and %i.", i);
  }
  oadd[0] = CrAnMask + 1;
  oadd[1] = 0;
  RealOperatorAddLaderNormalOrder(&O[0], oadd, 2, 0.5, &QDetNotOrderedFcFaBcBa);
  oadd[0] = CrAnMask;
  oadd[1] = 1;
  RealOperatorAddLaderNormalOrder(&O[0], oadd, 2, 0.5, &QDetNotOrderedFcFaBcBa);
  for (i = 1; i < Set.Cluster.NCell; i++)
  {
    oadd[0] = CrAnMask + 1;
    oadd[1] = 0;
    RealOperatorAddLaderNormalOrder(&O[i], oadd, 2, sqrt(0.125), &QDetNotOrderedFcFaBcBa);
    oadd[0] = CrAnMask;
    oadd[1] = 1;
    RealOperatorAddLaderNormalOrder(&O[i], oadd, 2, sqrt(0.125), &QDetNotOrderedFcFaBcBa);
    oadd[0] = CrAnMask + 1 + 2 * i;
    oadd[1] = 2 * i;
    RealOperatorAddLaderNormalOrder(&O[i], oadd, 2, sqrt(0.125), &QDetNotOrderedFcFaBcBa);
    oadd[0] = CrAnMask + 2 * i;
    oadd[1] = 1 + 2 * i;
    RealOperatorAddLaderNormalOrder(&O[i], oadd, 2, sqrt(0.125), &QDetNotOrderedFcFaBcBa);
  }
  return 0;
}

uint32_t CreateHeisenbergModelSpinFlipOperatorMomentumSpace(HeisenbergModelSettings Set, OperatorType * O)
{
  uint16_t i, k;
  uint16_t oadd[2];
  for (i = 0; i < Set.Cluster.NCell; i++)
  {
    InitOperator(&O[i]);
    O[i].NFermion = 2 * Set.Cluster.NCell;
    O[i].NWord = Min(O[i].NFermion, 1 + (O[i].NFermion - 1) / 8);
    InitOperatorFixedLength(&O[i], 2, 8);
    sprintf(O[i].Name, "Operator Sz at k: %i.", i);
  }
  for (k = 0; k < Set.Cluster.NCell; k++)
  {
    for (i = 0; i < Set.Cluster.NCell; i++)
    {
      oadd[0] = CrAnMask + 1 + 2 * i;
      oadd[1] = 1 + 2 * i;
      ComplexOperatorAddLaderNormalOrder(&O[k], oadd, 2,
                                         cos(Dot(Set.Cluster.kVector[k], Set.Cluster.CellPos[i])) * sqrt(0.25 /
                                                                                                         ((double)Set.
                                                                                                          Cluster.
                                                                                                          NCell)),
                                         sin(Dot(Set.Cluster.kVector[k], Set.Cluster.CellPos[i])) * sqrt(0.25 /
                                                                                                         ((double)Set.
                                                                                                          Cluster.
                                                                                                          NCell)),
                                         &QDetNotOrderedFcFaBcBa);
      oadd[0] = CrAnMask + 2 * i;
      oadd[1] = 2 * i;
      ComplexOperatorAddLaderNormalOrder(&O[k], oadd, 2,
                                         -cos(Dot(Set.Cluster.kVector[k], Set.Cluster.CellPos[i])) * sqrt(0.25 /
                                                                                                          ((double)Set.
                                                                                                           Cluster.
                                                                                                           NCell)),
                                         -sin(Dot(Set.Cluster.kVector[k], Set.Cluster.CellPos[i])) * sqrt(0.25 /
                                                                                                          ((double)Set.
                                                                                                           Cluster.
                                                                                                           NCell)),
                                         &QDetNotOrderedFcFaBcBa);
    }
  }
  return 0;
}

uint32_t CreateHeisenbergModelSpinFlipOperatorSxMomentumSpace(HeisenbergModelSettings Set, OperatorType * O)
{
  uint16_t i, k;
  uint16_t oadd[2];
  for (i = 0; i < Set.Cluster.NCell; i++)
  {
    InitOperator(&O[i]);
    O[i].NFermion = 2 * Set.Cluster.NCell;
    O[i].NWord = Min(O[i].NFermion, 1 + (O[i].NFermion - 1) / 8);
    InitOperatorFixedLength(&O[i], 2, 8);
    sprintf(O[i].Name, "Operator Sz at k: %i.", i);
  }
  for (k = 0; k < Set.Cluster.NCell; k++)
  {
    for (i = 0; i < Set.Cluster.NCell; i++)
    {
      oadd[0] = CrAnMask + 1 + 2 * i;
      oadd[1] = 2 * i;
      ComplexOperatorAddLaderNormalOrder(&O[k], oadd, 2,
                                         cos(Dot(Set.Cluster.kVector[k], Set.Cluster.CellPos[i])) * sqrt(0.25 /
                                                                                                         ((double)Set.
                                                                                                          Cluster.
                                                                                                          NCell)),
                                         sin(Dot(Set.Cluster.kVector[k], Set.Cluster.CellPos[i])) * sqrt(0.25 /
                                                                                                         ((double)Set.
                                                                                                          Cluster.
                                                                                                          NCell)),
                                         &QDetNotOrderedFcFaBcBa);
      oadd[0] = CrAnMask + 2 * i;
      oadd[1] = 1 + 2 * i;
      ComplexOperatorAddLaderNormalOrder(&O[k], oadd, 2,
                                         cos(Dot(Set.Cluster.kVector[k], Set.Cluster.CellPos[i])) * sqrt(0.25 /
                                                                                                         ((double)Set.
                                                                                                          Cluster.
                                                                                                          NCell)),
                                         sin(Dot(Set.Cluster.kVector[k], Set.Cluster.CellPos[i])) * sqrt(0.25 /
                                                                                                         ((double)Set.
                                                                                                          Cluster.
                                                                                                          NCell)),
                                         &QDetNotOrderedFcFaBcBa);
    }
  }
  return 0;
}

uint32_t CreateHeisenbergModelStartWaveFunctionOperator(HeisenbergModelSettings Set, OperatorType * O)
{
  uint16_t oadd[Set.Cluster.NCell];
  uint32_t i;
  InitOperator(O);
  O->NFermion = 2 * Set.Cluster.NCell;
  O->NWord = Min(O->NFermion, 1 + (O->NFermion - 1) / 8);
  InitOperatorFixedLength(O, Set.Cluster.NCell, 1);
  for (i = 0; i < Set.Cluster.NCell; i++)
  {
    oadd[i] =
      CrAnMask + (uint16_t) (2 * i +
                             ((Set.Cluster.CellPos[i][0] + Set.Cluster.CellPos[i][1] + Set.Cluster.CellPos[i][2]) % 2));
  }
  RealOperatorAddLaderNormalOrder(O, oadd, Set.Cluster.NCell, 1.0, &QDetNotOrderedFcFaBcBa);
  return 0;
}

uint32_t CreateHeisenbergModelStartWaveFunction(HeisenbergModelSettings Set, WaveFunctionType * Psi)
{
  OperatorType O;
  WaveFunctionType f;
  uint8_t *Determinant;
  CreateHeisenbergModelStartWaveFunctionOperator(Set, &O);
  f.QComplex = 0;
  f.NFermion = 2 * Set.Cluster.NCell;
  f.NBoson = 0;
  f.NBitsKey = Set.NBitsKey;
  WaveFunctionInit(&f);
  Determinant = calloc(f.NWord, sizeof(uint8_t));
  WaveFunctionAddElement(&f, Determinant, 1);
  OperatorPsi(O, f, Psi);
  FreeOperator(&O);
  return 0;
}

uint32_t CreateHeisenbergModelStartWaveFunctionOperatorFM(HeisenbergModelSettings Set, OperatorType * O)
{
  uint16_t oadd[Set.Cluster.NCell];
  uint32_t i;
  InitOperator(O);
  O->NFermion = 2 * Set.Cluster.NCell;
  O->NWord = Min(O->NFermion, 1 + (O->NFermion - 1) / 8);
  InitOperatorFixedLength(O, Set.Cluster.NCell, 1);
  for (i = 0; i < Set.Cluster.NCell; i++)
  {
    oadd[i] = CrAnMask + (uint16_t) (2 * i + 1);
  }
  RealOperatorAddLaderNormalOrder(O, oadd, Set.Cluster.NCell, 1.0, &QDetNotOrderedFcFaBcBa);
  return 0;
}

uint32_t CreateHeisenbergModelStartWaveFunctionFM(HeisenbergModelSettings Set, WaveFunctionType * Psi)
{
  OperatorType O;
  WaveFunctionType f;
  uint8_t *Determinant;
  CreateHeisenbergModelStartWaveFunctionOperatorFM(Set, &O);
  f.QComplex = 0;
  f.NFermion = 2 * Set.Cluster.NCell;
  f.NBoson = 0;
  f.NBitsKey = Set.NBitsKey;
  WaveFunctionInit(&f);
  Determinant = calloc(f.NWord, sizeof(uint8_t));
  WaveFunctionAddElement(&f, Determinant, 1);
  OperatorPsi(O, f, Psi);
  FreeOperator(&O);
  return 0;
}





int main(int argc, const char *argv[])
{
  double TotEmax, *EnergyList;
  uint32_t jpsi;
  uint32_t i;
  MatrixType SpinRho, Rho;
  WaveFunctionType *PsiList;
  HeisenbergModelSettings Set;
  WaveFunctionType OPsi;
  OperatorType H, *T;
  TriDiagonalMatrixType s;
  SpectraType *Spec;
  TimeType Time;
  InitializeTiming(&Time);
  BeginTiming("Total", &Time);
  TimeStampAll();
  FileReadHeisenbergModelSettings(&Set,
                                  "/Users/maurits/Documents/Programs/Model_lib/Model/build/Debug/Heisenberg.Settings");
  PrintHeisenbergModelSettings(Set);
  BeginTiming("Create Hamiltonian", &Time);
  CreateHeisenbergModelHamiltonian(Set, &H);
  EndTiming("Create Hamiltonian", &Time);
  PrintOperatorReadable(H);
  TotEmax = fabs(Set.Cluster.NCell * Set.J * 8);
  PsiList = malloc(Set.NPsi * sizeof(WaveFunctionType));
  EnergyList = malloc(Set.NPsi * sizeof(double));
  for (jpsi = 0; jpsi < Set.NPsi; jpsi++)
  {
    if (Set.J > 0)
      CreateHeisenbergModelStartWaveFunction(Set, &PsiList[jpsi]);
    else
      CreateHeisenbergModelStartWaveFunctionFM(Set, &PsiList[jpsi]);
    PrintWaveFunctionReadable(PsiList[jpsi]);
    BeginTiming("LanczosGroundState", &Time);
    LanczosGroundState(H, &PsiList[jpsi], &Set.zeroLancWavefunction, Set.zeroLancEnergy, TotEmax, PsiList, jpsi);
    EndTiming("LanczosGroundState", &Time);
    BeginTiming("LanczosGroundState", &Time);
    LanczosGroundState(H, &PsiList[jpsi], &Set.zeroLancWavefunction, Set.zeroLancEnergy, TotEmax, PsiList, jpsi);
    EndTiming("LanczosGroundState", &Time);
    if (OperatorPsiConserveBasis(H, PsiList[jpsi], &OPsi))
      AllocFailed();
    EnergyList[jpsi] = WaveFunctionDotProduct(OPsi, PsiList[jpsi]);
    PrintWaveFunctionReadableTruncated(PsiList[jpsi], 0.1);
  }
  printf("Orthonormality test of wave-functions\n");
  for (jpsi = 0; jpsi < Set.NPsi; jpsi++)
  {
    for (i = 0; i < Set.NPsi; i++)
    {
      printf(" %22.15E", WaveFunctionDotProduct(PsiList[i], PsiList[jpsi]));
    }
    printf("\n");
  }
  printf("State#   Energy\n");
  for (jpsi = 0; jpsi < Set.NPsi; jpsi++)
  {
    printf("%2i      %22.15E\n", jpsi, EnergyList[jpsi]);
  }
  for (jpsi = 0; jpsi < Set.NPsi; jpsi++)
  {
    CalculateSpinCorrelationMatrix(&SpinRho, PsiList[jpsi]);
    PrintMatrix(SpinRho);
    FreeMatrix(&SpinRho);
    CalculateDensityMatrix(&Rho, PsiList[jpsi]);
    PrintMatrix(Rho);
    FreeMatrix(&Rho);
  }
  s.N = Set.NTriDiag;
  InitTriDiagonalMatrix(&s);
  T = malloc(Set.Cluster.NCell * sizeof(OperatorType));
  Spec = malloc(Set.Cluster.NCell * sizeof(SpectraType));
  Spec[0].N = Set.NE;
  Spec[0].Emin = Set.Emin;
  Spec[0].Emax = Set.Emax;
  Spec[0].dE = (Set.Emax - Set.Emin) / Set.NE;
  Spec[0].EminPole = Set.Emin;
  Spec[0].EmaxPole = Set.Emax;
  Spec[0].Gamma = Set.Gamma;
  for (i = 0; i < Set.Cluster.NCell; i++)
  {
    SpectrumAllocate(Spec[0], &Spec[i]);
  }
  if (Set.LocalFourierBasis == 1)
  {
    if (Set.J > 0)
      CreateHeisenbergModelSpinFlipOperatorRealSpace(Set, T);
    else
      CreateHeisenbergModelSpinFlipOperatorSxRealSpace(Set, T);
  }
  else
  {
    if (Set.J > 0)
      CreateHeisenbergModelSpinFlipOperatorMomentumSpace(Set, T);
    else
      CreateHeisenbergModelSpinFlipOperatorSxMomentumSpace(Set, T);
  }


  for (i = 0; i < Set.Cluster.NCell; i++)
  {
    PrintOperatorReadable(T[i]);
  }
  for (jpsi = 0; jpsi < Set.NPsi; jpsi++)
  {
    for (i = 0; i < Set.Cluster.NCell; i++)
    {
      BeginTiming("Lanczos Tri Diagonalization", &Time);
      OperatorPsi(T[i], PsiList[jpsi], &OPsi);

      LanczosTriDiagonalize(H, OPsi, Set.zeroLancWavefunction, &s, Set.Gamma);
      TriDiagonalMatrixToSpectrum(EnergyList[jpsi], s, &Spec[i]);
      EndTiming("Lanczos Tri Diagonalization", &Time);
    }
    if (Set.LocalFourierBasis == 1)
    {
      SiteFourierTransform(Set.Cluster.NCell, Set.Cluster.CellPos, Set.Cluster.rSuperCell, Set.Cluster.kVector, Spec);
    }
    FilePrintSpectra(Spec, Set.Cluster.NCell, Set.SpecFile, "w");
  }
  EndTiming("Total", &Time);
  PrintTiming(Time);

  return 0;
}
