/*
*   Test.c contains a simple interface to test routines that test simple parts of the program
*   In principle I try to suply each subroutine with a test routine that tests it expected
*   behaviour on correct input and on unexpected input
*   As alwyas, I legg behind in test subroutines....
*
*   Test.c basically links everything togehter, the real subroutines are in other files
*/
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

#include "../Test/TestSpectral.c"
#include "../Test/TestBasicMath.c"
#include "../Test/TestCore.c"
#include "../Test/PointerAritmethic.c"
#include "../Test/TestFinitsizeClusterModels.c"
#include "../Test/TestIO.c"
#include "../Test/LanczosDyson.c"
#include "../Test/TestGreenFunctions.c"
#include "../Test/TestGFBA.c"
#include "../Test/TestStandardOperators.c"
#include "../Test/TestAtomicPhysics.c"
#include "../Test/OperatorArithmetic.c"
#include "../Test/StuttgartNMTO.c"
#include "../Test/HarmonicOscillator.c"
#include "../Test/BasisOptimization.c"

int main(int argc, const char *argv[])
{
  TimeStampAll();
  argc -= 1;
  for (; argc != 0; argc--)
  {
    /* functions in TestCore */
    if (!strncmp("Struct", argv[argc], strlen(argv[argc])))
      TestStruct();
    if (!strncmp("Hash", argv[argc], strlen(argv[argc])))
      TestHash();
    if (!strncmp("Wave", argv[argc], strlen(argv[argc])))
      TestWaveFunctionArithmetic();
    if (!strncmp("WaveAdd", argv[argc], strlen(argv[argc])))
      TestAddFermion();
    if (!strncmp("Opar", argv[argc], strlen(argv[argc])))
      TestOperatorArithmetic();
    if (!strncmp("OOpar", argv[argc], strlen(argv[argc])))
      TestOperatorArithmeticNormalOrder();
    if (!strncmp("Cell", argv[argc], strlen(argv[argc])))
      TestCellFunctions();



    if (!strncmp("Cr", argv[argc], strlen(argv[argc])))
      TestCreate();
    if (!strncmp("VMP", argv[argc], strlen(argv[argc])))
      TestVectorMatrixParallel();
    if (!strncmp("point", argv[argc], strlen(argv[argc])))
      TestPointerArithmetic();
    if (!strncmp("math", argv[argc], strlen(argv[argc])))
      TestBasicMath();
    if (!strncmp("Complex", argv[argc], strlen(argv[argc])))
      TestComplex();

    if (!strncmp("AllMatrix", argv[argc], strlen(argv[argc])))
      TestAllMatrix();
    if (!strncmp("MatrixVector", argv[argc], strlen(argv[argc])))
      TestMatrixVector();
    if (!strncmp("CompactMatrixVector", argv[argc], strlen(argv[argc])))
      TestCompactMatrixVector();
    if (!strncmp("ConvertMatrix", argv[argc], strlen(argv[argc])))
      TestMatrixTypeConversion();
    if (!strncmp("Combi", argv[argc], strlen(argv[argc])))
      TestCombinatorial();
    if (!strncmp("3J", argv[argc], strlen(argv[argc])))
      TestThreeJSymbol();

    if (!strncmp("DGESVD", argv[argc], strlen(argv[argc])))
      TestDGESVD();
    if (!strncmp("Lanczos", argv[argc], strlen(argv[argc])))
      TestLanczos();

    if (!strncmp("GFSpec", argv[argc], strlen(argv[argc])))
      TestGreenFunctionArithmetic();
    if (!strncmp("GFCSpec", argv[argc], strlen(argv[argc])))
      TestGreenFunctionCreation();
    if (!strncmp("Spec", argv[argc], strlen(argv[argc])))
      TestSpec();
    if (!strncmp("SpecSelf", argv[argc], strlen(argv[argc])))
      TestSpecSelf();
    if (!strncmp("KKTrans", argv[argc], strlen(argv[argc])))
      TestKKTrans();
    if (!strncmp("interpolate", argv[argc], strlen(argv[argc])))
      TestInterpolate();

    if (!strncmp("Inverse", argv[argc], strlen(argv[argc])))
      TestMatrixInverse();

    if (!strncmp("bmp", argv[argc], strlen(argv[argc])))
      TestBMP();
    if (!strncmp("waveIO", argv[argc], strlen(argv[argc])))
      TestWaveIO();

    if (!strncmp("TriMatrix", argv[argc], strlen(argv[argc])))
      TestTriDiagonalMatrix();

    if (!strncmp("LancDy", argv[argc], strlen(argv[argc])))
      TestLanczosDyson();
    if (!strncmp("LancDy2", argv[argc], strlen(argv[argc])))
      TestLanczosDyson2();
    if (!strncmp("LancDy3", argv[argc], strlen(argv[argc])))
      TestLanczosDyson3();

    if (!strncmp("GenEig", argv[argc], strlen(argv[argc])))
      TestGeneralizedEigensystem();

    if (!strncmp("Green", argv[argc], strlen(argv[argc])))
      TestGreenFunctionArithmetic();

    if (!strncmp("GFBA", argv[argc], strlen(argv[argc])))
      TestGFBA();

    if (!strncmp("StdOppS", argv[argc], strlen(argv[argc])))
      TestS();

    if (!strncmp("OppMult", argv[argc], strlen(argv[argc])))
      TestOppMult();
    if (!strncmp("Rotate", argv[argc], strlen(argv[argc])))
      TestRotate();
    
    if (!strncmp("AtomicPhysics", argv[argc], strlen(argv[argc])))
      TestAtomicPhysics();

    if (!strncmp("HarmOsc", argv[argc], strlen(argv[argc])))
      TestHarmonicOscillator();
    
    if (!strncmp("NMTO", argv[argc], strlen(argv[argc])))
      TestStuttgartNMTO();
    
    if (!strncmp("BasisOpt", argv[argc], strlen(argv[argc])))
      BasisOptimization();
    
    if (!strncmp("MF", argv[argc], strlen(argv[argc])))
      TestMeanField();
  }
  
  
  return 0;
}
