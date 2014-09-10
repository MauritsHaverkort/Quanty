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
#include <complex.h>

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
#include "../Core/MeanField.c"
#include "../Core/Spectral.c"
#include "../Core/DensityMatrix.c"
#include "../Core/Random.c"
#include "../Core/MatrixOperator.c"
#include "../Model/StandardOperators.c"
#include "../Hubbard_Model/CellFunctions.c"
#include "../Tight_Binding/Tight_Binding_lib.c"
#include "../Tight_Binding/TightBindingToG.c"
#include "../Tight_Binding/ImpurityUnfolding.c"
#include "../Tight_Binding/StuttgartNMTO.c"
#include "../DMFT/SingleSiteFunctions.c"
#include "../IO/StuttgartNMTO.c"

#include "./ResonantSpectra.c"
#include "./NiO_Calcs.c"
#include "./MnO_Calcs.c"

int main(int argc, const char *argv[])
{
/*  TriDiagonalMatrixType SpecTri;
  double wRes, energy, gamma;
  MatrixType invMat;
  int i, j;

  SpecTri.N = 5;
  InitTriDiagonalMatrix(&SpecTri);
  for(i = 1; i <= 5; i++)
  {
    SpecTri.a[i] = i;
    if(i <= 4)
      SpecTri.b[i] = 1;
  }

  wRes = 10.0;
  energy = 20.0;
  gamma = 1.0;
  printf("Original matrix:\n");
  PrintPlainTriDiagonalMatrix(SpecTri);

  if(TriDiagonalMatrixToGreensFunction(&SpecTri, wRes, energy, gamma, &invMat))
  {
    printf("Problem with tridiagonalmatrixtogreensfunction\n");

  }
  else
  {
    printf("Successful run of tridiagonalmatrixtogreensfunction\n");
    printf("Greens function matrix is:\n");
    
    for(i = 0; i < invMat.Ni; i++)
    {    
      for(j = 0; j < invMat.Nj; j++)
      {
        printf("%f   %f    ", invMat.VRe[i][j], invMat.VIm[i][j]);

      }
      printf("\n");
    }

  }  
  */
//  NiOLigandField();
  MnOLigandField();
  return 0;

}
