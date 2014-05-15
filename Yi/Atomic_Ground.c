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
#include <Accelerate/Accelerate.h>

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
#include "../Tight_Binding/TightBindingToG.c"
#include "../DMFT/SingleSiteFunctions.c"
#include "./AtomicFunctions.c"

int main()
{
# define Naccuracy 3
  WaveFunctionType initPsi, PsiGrd[Naccuracy];
  OperatorType AtomicHam, OSz, OSsqr;
  MatrixType Rho;
  uint32_t NOrb, NBasis, Nelectron;
  uint32_t i;
  uint16_t *IndexUp, *IndexDn;
  double EGrd[Naccuracy], VarEGrd[Naccuracy], Sz[Naccuracy], Ssqr[Naccuracy];
  double accuracy, zero=1.e-15;

  Nelectron = 2;                                          /* atomic number Z */
  NOrb = 6;                                               /* see function AtomicCI for the ordering of orbitals */
  AtomicCI("RadMEListZ2.dat", Nelectron, NOrb, &NBasis, &AtomicHam);

  IndexUp=malloc(NBasis*sizeof(uint16_t));
  IndexDn=malloc(NBasis*sizeof(uint16_t));
  for (i=0; i<NBasis; i++)
  {
    IndexUp[i] = 2*i;
    IndexDn[i] = 2*i+1;
  }
  InitWavefunctionRandom(&initPsi, NBasis*2, Nelectron);

  /* loop over accuracies */
  for(i=0; i<Naccuracy; i++)
  {
    accuracy=1/pow(10, (double)(17-2*(Naccuracy-i)));
    printf("\nCalculating the ground-state with truncation of basis states smaller than : %22.15E\n", accuracy);
    /* create starting point wavefunction */
    if (i==0)
    {
      WaveFunctionCopy(initPsi, &PsiGrd[i]);
    }
    else
    {
      WaveFunctionCopy(PsiGrd[i-1], &PsiGrd[i]);
    }
    /* Find groundstate */
    LanczosGroundState(AtomicHam, &(PsiGrd[i]), &accuracy, zero, 0, &(PsiGrd[i]), 0);

//    PrintWaveFunctionReadable(PsiGrd[i]);

    PsiHermitianOperatorPsi(AtomicHam, PsiGrd[i], &EGrd[i]);
    PsiHermitianOperatorSqrPsi(AtomicHam, PsiGrd[i], &VarEGrd[i]);
    VarEGrd[i]-=Sqr(EGrd[i]);
    CreateOperatorSz(NBasis*2, NBasis, IndexUp, IndexDn, &OSz);
    CreateOperatorSsqr(NBasis*2, NBasis, IndexUp, IndexDn, &OSsqr);
    PsiHermitianOperatorPsi(OSz, PsiGrd[i], &Sz[i]);
    PsiHermitianOperatorPsi(OSsqr, PsiGrd[i], &Ssqr[i]);

    printf("\n E                      VarE                   Sz                   S^2                   Nstates\n");
    printf("%22.15E %22.15E %22.15E %22.15E %10i\n", EGrd[i], VarEGrd[i], Sz[i], Ssqr[i], PsiGrd[i].N);
  }
  /* print an overview of groundstates */
  printf("\n\n==================================\nOriginal Hamiltonian\n==================================\n\n");
  printf("\n E                      VarE                   Sz                   S^2                   Nstates       Accuracy\n");
  for(i=0; i<Naccuracy; i++)
  {
    accuracy=1/pow(10, (double)(17-2*(Naccuracy-i)));
    printf("%22.15E %22.15E %22.15E %22.15E %10i %22.15E\n", EGrd[i], VarEGrd[i], Sz[i], Ssqr[i], PsiGrd[i].N, accuracy);
  }
  /* Density matrix */
  CalculateDensityMatrix(&Rho, PsiGrd[Naccuracy-1]);
  PrintMatrix(Rho);

  return 0;
}
