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
#include "../Model/StandardOperators.c"
#include "../Hubbard_Model/CellFunctions.c"
#include "../Tight_Binding/Tight_Binding_lib.c"
#include "../Tight_Binding/ImpurityUnfolding.c"

int main (int argc, const char * argv[]) 
{
  TimeStampAll();
  argc-=1;
  char comments[1000];
  char filenamepre[255];
  FILE *stream;
  int Dimension,ClusterSize,Nrepetitions;
  double Rho,Vpot;
  if((stream=fopen(argv[argc],"r"))==NULL)
  {
    printf("ERROR\nCan't open %s for reading.\n",argv[argc]);
    return 1;
  }
  fscanf(stream,"%i" ,&(Dimension));               fgets(comments,1000,stream);
  fscanf(stream,"%i" ,&(ClusterSize));             fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(Rho));                     fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(Vpot));                    fgets(comments,1000,stream);
  fscanf(stream,"%i" ,&(Nrepetitions));            fgets(comments,1000,stream);
  fclose(stream);
  if(Rho<0)
  {
    printf("Negative impurity denisty in input file\n");
    printf("Use single impurity calculation\n");
    sprintf(filenamepre,"D%i_N%i_V%5.3f_SingleImpurity",Dimension,ClusterSize,Vpot);
    if(Dimension==1) TightBindingImpurityUnfolding1Impurity1D(ClusterSize, Vpot, filenamepre);
    if(Dimension==2) TightBindingImpurityUnfolding1Impurity2D(ClusterSize, Vpot, filenamepre);
    if(Dimension==3) TightBindingImpurityUnfolding1Impurity3D(ClusterSize, Vpot, filenamepre);
  }
  else
  {
    sprintf(filenamepre,"D%i_N%i_Rho%5.3f_V%5.3f_Avv%i",Dimension,ClusterSize,Rho,Vpot,Nrepetitions);
    TightBindingImpurityUnfoldingRandomImpuritySeries(Dimension, ClusterSize, Rho, Vpot, Nrepetitions, filenamepre);
  }
  
//  TightBindingImpurityUnfolding1Impurity1D(3, 0.5, "1D_N02_V3.0");


//  TightBindingImpurityUnfoldingRandomImpuritySeries(1, 100, 0.5, 3.0, 10, "D1_N02_Rho0.5_V6.0_Avv10");
//  TightBindingImpurityUnfoldingRandomImpuritySeries(2,  10, 0.5, 6.0, 10, "D2_N02_Rho0.5_V6.0_Avv10");
//  TightBindingImpurityUnfoldingRandomImpuritySeries(3,   4, 0.5, 9.0, 10, "D3_N02_Rho0.5_V6.0_Avv10");


  return 0;
}
