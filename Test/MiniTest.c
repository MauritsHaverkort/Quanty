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




int main(int argc, const char *argv[])
{
  TimeStampAll();
  uint32_t i, j;
  uint32_t nElec = 7;
  uint32_t nOrb = 10;
  uint32_t nConfig[BinomialCoefficient(nOrb, nElec)];
  nConfig[0] = (1 << nElec) - 1;
  printf("#%4i | %4X\n", 0, nConfig[0]);
  j = 0;
  for (i = 1 << (nElec - 1); i > 0; i = i >> 1)
  {
    j++;
    nConfig[j] = nConfig[j - 1] + i;
    printf("#%4i | %4X\n", j, nConfig[j]);
  }


  return 0;
}
