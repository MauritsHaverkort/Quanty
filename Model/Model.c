#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <stdint.h> 
#include <string.h>
#include <math.h>
#include <unistd.h>

#include "../Core/StructDef.c"
#include "../Core/BasicMath.c"
#include "../Core/Error.c"
#include "../Core/Hash.c"
#include "../IO/Print.c"
#include "../IO/Read.c"
#include "../Core/VectorMatrix.c"
#include "../Core/ComplexVectorMatrix.c"
#include "../Core/Spectral.c"
#include "../Core/Lanczos.c"
#include "../Tools/TemperatureAvverage.c"

int main (int argc, const char * argv[]) 
{
  if (argc<=2)
  {
    printf("Need at least two arguments in order to do something\n");
    return 1;
  };
  if (!strncmp("Temperature",argv[1],strlen(argv[1]))) TemperatureAvverage(argv[2]);
  return 0;

}
