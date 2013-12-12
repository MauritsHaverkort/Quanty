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

int main(int argc, const char *argv[])
{
  TimeStampAll();
  argc -= 1;
  DMFTXASSettings Set;
  AndersonMatrixType A;
  DoubleTriDiagonalMatrixType BathFull, Bath;
  TriDiagonalMatrixType s;
  MatrixType Ttot;
  OperatorType OdTri, OdTriQ, Tup, Tdn;
  WaveFunctionType psidTri, Tpsi;
  SpectraType Spec[2];
  uint16_t oadd[4];
  double mu = 0.0, epsilon, Energy, Hist[20], shiftin, RhodTri[100], HistTot;
  uint32_t i, nFermion, nBoson;

  FileReadDMFTXASSettings(&Set, argv[argc]);
  PrintDMFTXASSettings(Set);

  FileReadAndersonMatrixType(&A, Set.BathFilename);
  strcpy(A.Name, "Input Anderson Impurity Matrix");
  PrintAndersonMatrix(A);
  AndersonMatrixToDoubleTriDiagonalMatrix(A, &Ttot, &BathFull, 0.0, 0.01);
  DoubleTriDiagonalMatrixPart(&Bath, BathFull, Min(BathFull.val.N, Set.NBathMax), Min(BathFull.con.N, Set.NBathMax));
  PrintDoubleTriDiagonalMatrix(Bath);

  DoubleTriDiagonalMatrixToAndersonImpurityOperator(&OdTri, Bath, Bath, Set.U, mu);
  PrintOperatorReadable(OdTri);

  Spec[0].N = Set.NE;
  Spec[0].Emin = Set.Emin;
  Spec[0].Emax = Set.Emax;
  Spec[0].dE = (Spec[0].Emax - Spec[0].Emin) / Spec[0].N;
  Spec[0].Gamma = Set.Gamma;
  SpectrumAllocate(Spec[0], &Spec[0]);
  SpectrumAllocate(Spec[0], &Spec[1]);

  InitOperatorCopySettings(&OdTriQ, OdTri);
  InitOperatorCopySettings(&Tup, OdTri);
  InitOperatorCopySettings(&Tdn, OdTri);
  OperatorAdd(&OdTriQ, OdTri, 1.0, 0.0);
  oadd[0] = CrAnMask + 0;
  RealOperatorAddLaderNormalOrder(&Tup, oadd, 1, 1.0, &QDetNotOrderedFcFaBcBa);
  oadd[0] = CrAnMask + 1;
  RealOperatorAddLaderNormalOrder(&Tdn, oadd, 1, 1.0, &QDetNotOrderedFcFaBcBa);

  s.N = Set.NTriDiag;
  InitTriDiagonalMatrix(&s);

  oadd[1] = CrAnMask + 0;
  oadd[0] = 0;
  RealOperatorAddLaderNormalOrder(&OdTriQ, oadd, 2, -Set.Q, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 1;
  oadd[0] = 1;
  RealOperatorAddLaderNormalOrder(&OdTriQ, oadd, 2, -Set.Q, &QDetNotOrderedFcFaBcBa);

  nFermion = OdTri.NFermion / 2;
  nBoson = 0;
  CreateStartWaveFunctionOnsiteEnergy(OdTri, &psidTri, nFermion, nBoson, Set.NBitsKey);

  CalculateDensity(RhodTri, psidTri);
  printf("RhodTri\n");
  for (i = 0; i < OdTri.NFermion; i++)
  {
    printf("%i %lf\n", i, RhodTri[i]);
  }

  shiftin = 0.0;
  epsilon = Set.zeroLancWavefunction;
  LanczosGroundState(OdTri, &psidTri, &epsilon, Set.zeroLancEnergy, shiftin, &psidTri, 0);
  LanczosGroundState(OdTri, &psidTri, &epsilon, Set.zeroLancEnergy, shiftin, &psidTri, 0);
  PsiHermitianOperatorPsi(OdTri, psidTri, &Energy);
  printf("Ground-state energy is: %lf", Energy);
  PrintWaveFunctionReadableTruncated(psidTri, 0.1);

  CalculateDensity(RhodTri, psidTri);

  printf("RhodTri\n");
  for (i = 0; i < OdTri.NFermion; i++)
  {
    printf("%i %lf\n", i, RhodTri[i]);
  }

  CalculatePrefactorHistogram(Hist, psidTri, 10, 20);
  printf("Logarithmic Prefactor Histogram on base 10\n");
  HistTot = 0;
  for (i = 0; i < 20; i++)
  {
    HistTot += Hist[i];
    printf("%3i %lf\n", i, Hist[i]);
  }
  printf("tot %lf\n", HistTot);

  OperatorPsi(Tup, psidTri, &Tpsi);
  LanczosTriDiagonalize(OdTriQ, Tpsi, epsilon, &s, Set.Gamma);
  TriDiagonalMatrixToSpectrum(Energy, s, &Spec[0]);

  OperatorPsi(Tdn, psidTri, &Tpsi);
  LanczosTriDiagonalize(OdTriQ, Tpsi, epsilon, &s, Set.Gamma);
  TriDiagonalMatrixToSpectrum(Energy, s, &Spec[1]);

  FilePrintSpectra(Spec, 2, Set.SpecFile, "w");

  return 0;
}
