#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <unistd.h>

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


/*
*   Operators are given as structures as defined in StructDef.c
*   The one particle basis is given by "orbitals", which stands for sites, spin, local-degeneracies, k-points, etc.
*   The many particle basis is given by Slater integrals that can be seen as a string of occupied(1) and unoccupied(0) orbitals
*   Operators are given by two matricess M_{ij} and M_{ijkl} whereby the operator is:
*   \sum_{i,j} M_{i,j} a^{dag}_{i}a_{j} + \sum_{i,j,k,l} a^{dag}_{i}a^{dag}_{j}a_{k}a_{l}
*   with a^{dag} and a the usual creation and annihilation operators respectively
*
*   The folowing function creates the operator for a one band Hubbard model with nearest neighbor hopping t, onsite repulsion U
*   with 'atoms' at all integer index positions and a unit cell given by the vectors a, b, and c.
*      a vector (0,0,1) does not result in onsite hopping.
*      1D can be realized by setting a=(n,0,0) and b=(0,1,0) c=(0,0,1)
*      2D etc.
*/
/*
OperatorType CreateSimpelHubbardModel(const HubbardModelSettings HSet)
{
  OperatorType H;
  int32_t i,j,ri,rj,rk;
  double  r[3];
  H.QComplex=0;
  H.NFermion=HSet.NSite;
  H.NBoson=0;
  InitOperator(&H);
  H.NCrCrAnAn=H.NFermion;
  H.VCrCrAnAnRe=malloc(H.NFermion*sizeof(double));
  H.iCrCrAnAn=malloc(H.NFermion*sizeof(uint32_t));
  H.jCrCrAnAn=malloc(H.NFermion*sizeof(uint32_t));
  H.kCrCrAnAn=malloc(H.NFermion*sizeof(uint32_t));
  H.lCrCrAnAn=malloc(H.NFermion*sizeof(uint32_t));
  for(i=0;i<H.NFermion;i++)
  {
    H.VCrCrAnAnRe[i]=HSet.U;
    H.iCrCrAnAn[i]=i;
    H.jCrCrAnAn[i]=i+H.NFermion;
    H.kCrCrAnAn[i]=i+H.NFermion;
    H.lCrCrAnAn[i]=i;
  }
  H.NCrAn=0;
  H.VCrAnRe=malloc(12*H.NFermion*sizeof(double));
  H.iCrAn=malloc(12*H.NFermion*sizeof(uint32_t));
  H.jCrAn=malloc(12*H.NFermion*sizeof(uint32_t));
  for(i=0;i<H.NFermion;i++)
  {
    H.VCrAnRe[H.NCrAn]=-0.5*HSet.U;
    H.iCrAn[H.NCrAn]=i;
    H.jCrAn[H.NCrAn]=i;
    H.NCrAn++;
    H.VCrAnRe[H.NCrAn]=-0.5*HSet.U;
    H.iCrAn[H.NCrAn]=i+H.NFermion;
    H.jCrAn[H.NCrAn]=i+H.NFermion;
    H.NCrAn++;
    for(j=0;j<i;j++)
    {
      for(ri=-1;ri<=1;ri++)
      {
        for(rj=-1;rj<=1;rj++)
        {
          for(rk=-1;rk<=1;rk++)
          {
            r[0]=HSet.AtomPos[j][0]+ri*HSet.rUnitCell[0][0]+rj*HSet.rUnitCell[1][0]+rk*HSet.rUnitCell[2][0];
            r[1]=HSet.AtomPos[j][1]+ri*HSet.rUnitCell[0][1]+rj*HSet.rUnitCell[1][1]+rk*HSet.rUnitCell[2][1];
            r[2]=HSet.AtomPos[j][2]+ri*HSet.rUnitCell[0][2]+rj*HSet.rUnitCell[1][2]+rk*HSet.rUnitCell[2][2];
            if(DistanceSqr(HSet.AtomPos[i],r)==1)
            {
              H.VCrAnRe[H.NCrAn]=HSet.t;
              H.iCrAn[H.NCrAn]=i;
              H.jCrAn[H.NCrAn]=j;
              H.NCrAn++;
              H.VCrAnRe[H.NCrAn]=HSet.t;
              H.iCrAn[H.NCrAn]=i+H.NFermion;
              H.jCrAn[H.NCrAn]=j+H.NFermion;
              H.NCrAn++;
              H.VCrAnRe[H.NCrAn]=HSet.t;
              H.iCrAn[H.NCrAn]=j;
              H.jCrAn[H.NCrAn]=i;
              H.NCrAn++;
              H.VCrAnRe[H.NCrAn]=HSet.t;
              H.iCrAn[H.NCrAn]=j+H.NFermion;
              H.jCrAn[H.NCrAn]=i+H.NFermion;
              H.NCrAn++;
            }
          }
        }
      }
    }
  }
  H.VCrAnRe=realloc(H.VCrAnRe,H.NCrAn*sizeof(double));
  H.iCrAn=realloc(H.iCrAn,H.NCrAn*sizeof(uint32_t));
  H.jCrAn=realloc(H.jCrAn,H.NCrAn*sizeof(uint32_t));
  H.NFermion*=2;
  H.NWord=Min(H.NFermion,1+(H.NFermion-1)/8);
  printf("\nHubbard Model Hamiltonian created\n");
  printf("%i atoms created\n",H.NFermion/2);
  printf("Within the unit cell spanned by\n");
  printf("%4i %4i %4i\n",HSet.rUnitCell[0][0],HSet.rUnitCell[0][1],HSet.rUnitCell[0][2]);
  printf("%4i %4i %4i\n",HSet.rUnitCell[1][0],HSet.rUnitCell[1][1],HSet.rUnitCell[1][2]);
  printf("%4i %4i %4i\n",HSet.rUnitCell[2][0],HSet.rUnitCell[2][1],HSet.rUnitCell[2][2]);
  printf("At positions\n");
  printf("   x    y    z  spin index\n");
  for(i=0;i<H.NFermion/2;i++)
  {
    printf("%4i %4i %4i    dn  %4i\n",HSet.AtomPos[i][0],HSet.AtomPos[i][1],HSet.AtomPos[i][2],i);
  }
  for(i=0;i<H.NFermion/2;i++)
  {
    printf("%4i %4i %4i    up  %4i\n",HSet.AtomPos[i][0],HSet.AtomPos[i][1],HSet.AtomPos[i][2],i+H.NFermion/2);
  }
  printf("The Hamiltonian has %i 'orbitals' (cites*spin*local_degeneracy*etc.)\n",H.NFermion);
  printf("The two particle part of H:\n");
  for(i=0;i<H.NCrCrAnAn;i++)
  {
    printf("%4i %4i %4i %4i %10.6f\n",H.iCrCrAnAn[i],H.jCrCrAnAn[i],H.kCrCrAnAn[i],H.lCrCrAnAn[i],H.VCrCrAnAnRe[i]);
  }
  printf("The one particle part of H:\n");
  for(i=0;i<H.NCrAn;i++)
  {
    printf("%4i %4i %10.6f\n",H.iCrAn[i],H.jCrAn[i],H.VCrAnRe[i]);
  }
  printf("\n");
  return H;
}
*/
/*
OperatorType CreateSimpelTwistedHubbardModel(const HubbardModelSettings HSet, const uint32_t iq[])
{
  OperatorType H;
  int32_t AtomPos[100][3];
  int32_t i,j,k,ri,rj,rk;
  double   a,b,c,r[3],g[3][3],rUnitCell[3][3],phase;
  H.NFermion=0;
  InitOperator(&H);
  // Create reciprocal lattice vectors 
  rUnitCell[0][0]=(double)HSet.rUnitCell[0][0];
  rUnitCell[0][1]=(double)HSet.rUnitCell[0][1];
  rUnitCell[0][2]=(double)HSet.rUnitCell[0][2];
  rUnitCell[1][0]=(double)HSet.rUnitCell[1][0];
  rUnitCell[1][1]=(double)HSet.rUnitCell[1][1];
  rUnitCell[1][2]=(double)HSet.rUnitCell[1][2];
  rUnitCell[2][0]=(double)HSet.rUnitCell[2][0];
  rUnitCell[2][1]=(double)HSet.rUnitCell[2][1];
  rUnitCell[2][2]=(double)HSet.rUnitCell[2][2];
  ResiprocalLatticeVector(rUnitCell, g);
  // determine the unit cell by a loop over i=min(x) to max(x), j=min(y) to max(y) and k=min(z) to max(z) 
  // and then selecting only the points that lay witnin the unit cell 
  for(i=-(abs(HSet.rUnitCell[0][0])+abs(HSet.rUnitCell[1][0])+abs(HSet.rUnitCell[2][0]));i<(abs(HSet.rUnitCell[0][0])+abs(HSet.rUnitCell[1][0])+abs(HSet.rUnitCell[2][0]));i++)
  {
    for(j=-(abs(HSet.rUnitCell[0][1])+abs(HSet.rUnitCell[1][1])+abs(HSet.rUnitCell[2][1]));j<(abs(HSet.rUnitCell[0][1])+abs(HSet.rUnitCell[1][1])+abs(HSet.rUnitCell[2][1]));j++)
    {
      for(k=-(abs(HSet.rUnitCell[0][2])+abs(HSet.rUnitCell[1][2])+abs(HSet.rUnitCell[2][2]));k<(abs(HSet.rUnitCell[0][2])+abs(HSet.rUnitCell[1][2])+abs(HSet.rUnitCell[2][2]));k++)
      {
        a=-((double)(HSet.rUnitCell[1][2]*HSet.rUnitCell[2][1]*i - HSet.rUnitCell[1][1]*HSet.rUnitCell[2][2]*i - 
          HSet.rUnitCell[0][2]*HSet.rUnitCell[2][1]*j + HSet.rUnitCell[0][1]*HSet.rUnitCell[2][2]*j + HSet.rUnitCell[0][2]*HSet.rUnitCell[1][1]*k - 
          HSet.rUnitCell[0][1]*HSet.rUnitCell[1][2]*k)/(double)
        (-(HSet.rUnitCell[0][2]*HSet.rUnitCell[1][1]*HSet.rUnitCell[2][0]) + HSet.rUnitCell[0][1]*HSet.rUnitCell[1][2]*HSet.rUnitCell[2][0] + 
          HSet.rUnitCell[0][2]*HSet.rUnitCell[1][0]*HSet.rUnitCell[2][1] - HSet.rUnitCell[0][0]*HSet.rUnitCell[1][2]*HSet.rUnitCell[2][1] - 
          HSet.rUnitCell[0][1]*HSet.rUnitCell[1][0]*HSet.rUnitCell[2][2] + HSet.rUnitCell[0][0]*HSet.rUnitCell[1][1]*HSet.rUnitCell[2][2]));
        b=-((double)(HSet.rUnitCell[1][2]*HSet.rUnitCell[2][0]*i - HSet.rUnitCell[1][0]*HSet.rUnitCell[2][2]*i - 
          HSet.rUnitCell[0][2]*HSet.rUnitCell[2][0]*j + HSet.rUnitCell[0][0]*HSet.rUnitCell[2][2]*j + HSet.rUnitCell[0][2]*HSet.rUnitCell[1][0]*k - 
          HSet.rUnitCell[0][0]*HSet.rUnitCell[1][2]*k)/(double)
        (HSet.rUnitCell[0][2]*HSet.rUnitCell[1][1]*HSet.rUnitCell[2][0] - HSet.rUnitCell[0][1]*HSet.rUnitCell[1][2]*HSet.rUnitCell[2][0] - 
          HSet.rUnitCell[0][2]*HSet.rUnitCell[1][0]*HSet.rUnitCell[2][1] + HSet.rUnitCell[0][0]*HSet.rUnitCell[1][2]*HSet.rUnitCell[2][1] + 
          HSet.rUnitCell[0][1]*HSet.rUnitCell[1][0]*HSet.rUnitCell[2][2] - HSet.rUnitCell[0][0]*HSet.rUnitCell[1][1]*HSet.rUnitCell[2][2]));
        c=-((double)(-(HSet.rUnitCell[1][1]*HSet.rUnitCell[2][0]*i) + HSet.rUnitCell[1][0]*HSet.rUnitCell[2][1]*i + 
          HSet.rUnitCell[0][1]*HSet.rUnitCell[2][0]*j - HSet.rUnitCell[0][0]*HSet.rUnitCell[2][1]*j - HSet.rUnitCell[0][1]*HSet.rUnitCell[1][0]*k + 
          HSet.rUnitCell[0][0]*HSet.rUnitCell[1][1]*k)/(double)
        (HSet.rUnitCell[0][2]*HSet.rUnitCell[1][1]*HSet.rUnitCell[2][0] - HSet.rUnitCell[0][1]*HSet.rUnitCell[1][2]*HSet.rUnitCell[2][0] - 
          HSet.rUnitCell[0][2]*HSet.rUnitCell[1][0]*HSet.rUnitCell[2][1] + HSet.rUnitCell[0][0]*HSet.rUnitCell[1][2]*HSet.rUnitCell[2][1] + 
          HSet.rUnitCell[0][1]*HSet.rUnitCell[1][0]*HSet.rUnitCell[2][2] - HSet.rUnitCell[0][0]*HSet.rUnitCell[1][1]*HSet.rUnitCell[2][2]));      
        if((a>=0)&&(a<1.0)&&(b>=0)&&(b<1.0)&&(c>=0)&&(c<1.0))
        {
          AtomPos[H.NFermion][0]=i;
          AtomPos[H.NFermion][1]=j;
          AtomPos[H.NFermion][2]=k;
          H.NFermion++;
        }
      }
    }
  }
  // setup memory for two electron part of the Hamiltonian //
  H.NCrCrAnAn=H.NFermion;
  H.VCrCrAnAnRe=malloc(H.NFermion*sizeof(double));
  H.VCrCrAnAnIm=malloc(H.NFermion*sizeof(double));
  H.iCrCrAnAn=malloc(H.NFermion*sizeof(uint32_t));
  H.jCrCrAnAn=malloc(H.NFermion*sizeof(uint32_t));
  H.kCrCrAnAn=malloc(H.NFermion*sizeof(uint32_t));
  H.lCrCrAnAn=malloc(H.NFermion*sizeof(uint32_t));
  // fill the two electron part of the Hamiltonian //
  for(i=0;i<H.NFermion;i++)
  {
    H.VCrCrAnAnRe[i]=HSet.U;
    H.VCrCrAnAnIm[i]=0.0;
    H.iCrCrAnAn[i]=i;
    H.jCrCrAnAn[i]=i+H.NFermion;
    H.kCrCrAnAn[i]=i+H.NFermion;
    H.lCrCrAnAn[i]=i;
  }
  // setup memory fo the one electron part of the Hamiltonian //
  H.NCrAn=0;
  H.VCrAnRe=malloc(12*H.NFermion*sizeof(double));
  H.VCrAnIm=malloc(12*H.NFermion*sizeof(double));
  H.iCrAn=malloc(12*H.NFermion*sizeof(uint32_t));
  H.jCrAn=malloc(12*H.NFermion*sizeof(uint32_t));
  for(i=0;i<H.NFermion;i++)
  {
    // onsite energy equal to -0.5 U //
    H.VCrAnRe[H.NCrAn]=-0.5*HSet.U;
    H.VCrAnIm[H.NCrAn]=0;
    H.iCrAn[H.NCrAn]=i;
    H.jCrAn[H.NCrAn]=i;
    H.NCrAn++;
    H.VCrAnRe[H.NCrAn]=-0.5*HSet.U;
    H.VCrAnIm[H.NCrAn]=0;
    H.iCrAn[H.NCrAn]=i+H.NFermion;
    H.jCrAn[H.NCrAn]=i+H.NFermion;
    H.NCrAn++;
    // set hopping //
    // loop over all other orbitals different from i (loop over j) //
    // loop over hopping to the unit cell on the left, the same unit cell, or the unit cell on the right (ri,rj,rk) //
    for(j=0;j<i;j++)
    {
      for(ri=-1;ri<=1;ri++)
      {
        for(rj=-1;rj<=1;rj++)
        {
          for(rk=-1;rk<=1;rk++)
          {
            r[0]=AtomPos[j][0]+ri*HSet.rUnitCell[0][0]+rj*HSet.rUnitCell[1][0]+rk*HSet.rUnitCell[2][0];
            r[1]=AtomPos[j][1]+ri*HSet.rUnitCell[0][1]+rj*HSet.rUnitCell[1][1]+rk*HSet.rUnitCell[2][1];
            r[2]=AtomPos[j][2]+ri*HSet.rUnitCell[0][2]+rj*HSet.rUnitCell[1][2]+rk*HSet.rUnitCell[2][2];
            if(DistanceSqr(AtomPos[i],r)==1)
            {
              r[0]-=AtomPos[i][0];
              r[1]-=AtomPos[i][1];
              r[2]-=AtomPos[i][2];
              phase=((double)Dot(g[0],r)*iq[0])/((double)HSet.Nq[0])+((double)Dot(g[1],r)*iq[1])/((double)HSet.Nq[1])+((double)Dot(g[2],r)*iq[2])/((double)HSet.Nq[2]);
              H.VCrAnRe[H.NCrAn]=HSet.t*cos(phase);
              H.VCrAnIm[H.NCrAn]=HSet.t*sin(phase);
              H.iCrAn[H.NCrAn]=i;
              H.jCrAn[H.NCrAn]=j;
              H.NCrAn++;
              H.VCrAnRe[H.NCrAn]=HSet.t*cos(phase);
              H.VCrAnIm[H.NCrAn]=HSet.t*sin(phase);
              H.iCrAn[H.NCrAn]=i+H.NFermion;
              H.jCrAn[H.NCrAn]=j+H.NFermion;
              H.NCrAn++;
              H.VCrAnRe[H.NCrAn]=HSet.t*cos(phase);
              H.VCrAnIm[H.NCrAn]=-HSet.t*sin(phase);
              H.iCrAn[H.NCrAn]=j;
              H.jCrAn[H.NCrAn]=i;
              H.NCrAn++;
              H.VCrAnRe[H.NCrAn]=HSet.t*cos(phase);
              H.VCrAnIm[H.NCrAn]=-HSet.t*sin(phase);
              H.iCrAn[H.NCrAn]=j+H.NFermion;
              H.jCrAn[H.NCrAn]=i+H.NFermion;
              H.NCrAn++;
            }
          }
        }
      }
    }
  }
  H.VCrAnRe=realloc(H.VCrAnRe,H.NCrAn*sizeof(double));
  H.VCrAnIm=realloc(H.VCrAnIm,H.NCrAn*sizeof(double));
  H.iCrAn=realloc(H.iCrAn,H.NCrAn*sizeof(uint32_t));
  H.jCrAn=realloc(H.jCrAn,H.NCrAn*sizeof(uint32_t));
  H.NFermion*=2;
  H.NWord=Min(H.NFermion,1+(H.NFermion-1)/8);
  printf("\nHubbard Model Hamiltonian created\n");
  printf("%i atoms created\n",H.NFermion/2);
  printf("Within the unit cell spanned by\n");
  printf("%4i %4i %4i\n",HSet.rUnitCell[0][0],HSet.rUnitCell[0][1],HSet.rUnitCell[0][2]);
  printf("%4i %4i %4i\n",HSet.rUnitCell[1][0],HSet.rUnitCell[1][1],HSet.rUnitCell[1][2]);
  printf("%4i %4i %4i\n",HSet.rUnitCell[2][0],HSet.rUnitCell[2][1],HSet.rUnitCell[2][2]);
  printf("At positions\n");
  printf("   x    y    z  spin index\n");
  for(i=0;i<H.NFermion/2;i++)
  {
    printf("%4i %4i %4i    dn  %4i\n",AtomPos[i][0],AtomPos[i][1],AtomPos[i][2],i);
  }
  for(i=0;i<H.NFermion/2;i++)
  {
    printf("%4i %4i %4i    up  %4i\n",AtomPos[i][0],AtomPos[i][1],AtomPos[i][2],i+H.NFermion/2);
  }
  printf("For the pseudo k point %10.6f %10.6f %10.6f (units of rec. lat. vec. (2 pi))\n",iq[0]/((double)HSet.Nq[0]),iq[1]/((double)HSet.Nq[1]),iq[2]/((double)HSet.Nq[2]));
  printf("For the pseudo k point %10.6f %10.6f %10.6f (carthesian coordinates  (2 pi))\n",(g[0][0]*iq[0]/((double)HSet.Nq[0])+g[1][0]*iq[1]/((double)HSet.Nq[1])+g[2][0]*iq[2]/((double)HSet.Nq[2]))/(2*Pi),(g[0][1]*iq[0]/((double)HSet.Nq[0])+g[1][1]*iq[1]/((double)HSet.Nq[1])+g[2][1]*iq[2]/((double)HSet.Nq[2]))/(2*Pi),(g[0][2]*iq[0]/((double)HSet.Nq[0])+g[1][2]*iq[1]/((double)HSet.Nq[1])+g[2][2]*iq[2]/((double)HSet.Nq[2]))/(2*Pi));
  printf("With pseudo reciprocal lattice vectors (units of 2 Pi)\n");
  printf("%10.6f %10.6f %10.6f\n",g[0][0]/(2.0*Pi),g[0][1]/(2.0*Pi),g[0][2]/(2.0*Pi));
  printf("%10.6f %10.6f %10.6f\n",g[1][0]/(2.0*Pi),g[1][1]/(2.0*Pi),g[1][2]/(2.0*Pi));
  printf("%10.6f %10.6f %10.6f\n",g[2][0]/(2.0*Pi),g[2][1]/(2.0*Pi),g[2][2]/(2.0*Pi));
  printf("The Hamiltonian has %i 'orbitals' (cites*spin*local_degeneracy*etc.)\n",H.NFermion);
  printf("The two particle part of H:\n");
  for(i=0;i<H.NCrCrAnAn;i++)
  {
    printf("%4i %4i %4i %4i %10.6f\n",H.iCrCrAnAn[i],H.jCrCrAnAn[i],H.kCrCrAnAn[i],H.lCrCrAnAn[i],H.VCrCrAnAnRe[i]);
  }
  printf("The one particle part of H:\n");
  for(i=0;i<H.NCrAn;i++)
  {
    if(H.VCrAnIm[i]<0)
    {
      printf("%4i %4i %10.6f - I%10.6f\n",H.iCrAn[i],H.jCrAn[i],H.VCrAnRe[i],-H.VCrAnIm[i]);
    }
    else
    {
      printf("%4i %4i %10.6f + I%10.6f\n",H.iCrAn[i],H.jCrAn[i],H.VCrAnRe[i],H.VCrAnIm[i]);
    }
  }
  printf("\n");
  return H;
}
*/

/*
void CreateSimpelHubbardModelOpticsT(const HubbardModelSettings HSet, OperatorType *const H)
{
  InitOperator(H);
  int32_t AtomPos[100][3];
  int32_t i,j,k,ri,rj,rk;
  double   a,b,c,r[3];
  H->NFermion=0;
  for(i=-(abs(HSet.rUnitCell[0][0])+abs(HSet.rUnitCell[1][0])+abs(HSet.rUnitCell[2][0]));i<(abs(HSet.rUnitCell[0][0])+abs(HSet.rUnitCell[1][0])+abs(HSet.rUnitCell[2][0]));i++)
  {
    for(j=-(abs(HSet.rUnitCell[0][1])+abs(HSet.rUnitCell[1][1])+abs(HSet.rUnitCell[2][1]));j<(abs(HSet.rUnitCell[0][1])+abs(HSet.rUnitCell[1][1])+abs(HSet.rUnitCell[2][1]));j++)
    {
      for(k=-(abs(HSet.rUnitCell[0][2])+abs(HSet.rUnitCell[1][2])+abs(HSet.rUnitCell[2][2]));k<(abs(HSet.rUnitCell[0][2])+abs(HSet.rUnitCell[1][2])+abs(HSet.rUnitCell[2][2]));k++)
      {
        a=-((double)(HSet.rUnitCell[1][2]*HSet.rUnitCell[2][1]*i - HSet.rUnitCell[1][1]*HSet.rUnitCell[2][2]*i - 
          HSet.rUnitCell[0][2]*HSet.rUnitCell[2][1]*j + HSet.rUnitCell[0][1]*HSet.rUnitCell[2][2]*j + HSet.rUnitCell[0][2]*HSet.rUnitCell[1][1]*k - 
          HSet.rUnitCell[0][1]*HSet.rUnitCell[1][2]*k)/(double)
        (-(HSet.rUnitCell[0][2]*HSet.rUnitCell[1][1]*HSet.rUnitCell[2][0]) + HSet.rUnitCell[0][1]*HSet.rUnitCell[1][2]*HSet.rUnitCell[2][0] + 
          HSet.rUnitCell[0][2]*HSet.rUnitCell[1][0]*HSet.rUnitCell[2][1] - HSet.rUnitCell[0][0]*HSet.rUnitCell[1][2]*HSet.rUnitCell[2][1] - 
          HSet.rUnitCell[0][1]*HSet.rUnitCell[1][0]*HSet.rUnitCell[2][2] + HSet.rUnitCell[0][0]*HSet.rUnitCell[1][1]*HSet.rUnitCell[2][2]));
        b=-((double)(HSet.rUnitCell[1][2]*HSet.rUnitCell[2][0]*i - HSet.rUnitCell[1][0]*HSet.rUnitCell[2][2]*i - 
          HSet.rUnitCell[0][2]*HSet.rUnitCell[2][0]*j + HSet.rUnitCell[0][0]*HSet.rUnitCell[2][2]*j + HSet.rUnitCell[0][2]*HSet.rUnitCell[1][0]*k - 
          HSet.rUnitCell[0][0]*HSet.rUnitCell[1][2]*k)/(double)
        (HSet.rUnitCell[0][2]*HSet.rUnitCell[1][1]*HSet.rUnitCell[2][0] - HSet.rUnitCell[0][1]*HSet.rUnitCell[1][2]*HSet.rUnitCell[2][0] - 
          HSet.rUnitCell[0][2]*HSet.rUnitCell[1][0]*HSet.rUnitCell[2][1] + HSet.rUnitCell[0][0]*HSet.rUnitCell[1][2]*HSet.rUnitCell[2][1] + 
          HSet.rUnitCell[0][1]*HSet.rUnitCell[1][0]*HSet.rUnitCell[2][2] - HSet.rUnitCell[0][0]*HSet.rUnitCell[1][1]*HSet.rUnitCell[2][2]));
        c=-((double)(-(HSet.rUnitCell[1][1]*HSet.rUnitCell[2][0]*i) + HSet.rUnitCell[1][0]*HSet.rUnitCell[2][1]*i + 
          HSet.rUnitCell[0][1]*HSet.rUnitCell[2][0]*j - HSet.rUnitCell[0][0]*HSet.rUnitCell[2][1]*j - HSet.rUnitCell[0][1]*HSet.rUnitCell[1][0]*k + 
          HSet.rUnitCell[0][0]*HSet.rUnitCell[1][1]*k)/(double)
        (HSet.rUnitCell[0][2]*HSet.rUnitCell[1][1]*HSet.rUnitCell[2][0] - HSet.rUnitCell[0][1]*HSet.rUnitCell[1][2]*HSet.rUnitCell[2][0] - 
          HSet.rUnitCell[0][2]*HSet.rUnitCell[1][0]*HSet.rUnitCell[2][1] + HSet.rUnitCell[0][0]*HSet.rUnitCell[1][2]*HSet.rUnitCell[2][1] + 
          HSet.rUnitCell[0][1]*HSet.rUnitCell[1][0]*HSet.rUnitCell[2][2] - HSet.rUnitCell[0][0]*HSet.rUnitCell[1][1]*HSet.rUnitCell[2][2]));      
        if((a>=0)&&(a<1.0)&&(b>=0)&&(b<1.0)&&(c>=0)&&(c<1.0))
        {
          AtomPos[H->NFermion][0]=i;
          AtomPos[H->NFermion][1]=j;
          AtomPos[H->NFermion][2]=k;
          H->NFermion++;
        }
      }
    }
  }
  H->NCrAn=0;
  H->VCrAnRe=malloc(12*H->NFermion*sizeof(double));
  H->iCrAn=malloc(12*H->NFermion*sizeof(uint32_t));
  H->jCrAn=malloc(12*H->NFermion*sizeof(uint32_t));
  for(i=0;i<H->NFermion;i++)
  {
    for(j=0;j<i;j++)
    {
      for(ri=-1;ri<=1;ri++)
      {
        for(rj=-1;rj<=1;rj++)
        {
          for(rk=-1;rk<=1;rk++)
          {
            r[0]=AtomPos[j][0]+ri*HSet.rUnitCell[0][0]+rj*HSet.rUnitCell[1][0]+rk*HSet.rUnitCell[2][0];
            r[1]=AtomPos[j][1]+ri*HSet.rUnitCell[0][1]+rj*HSet.rUnitCell[1][1]+rk*HSet.rUnitCell[2][1];
            r[2]=AtomPos[j][2]+ri*HSet.rUnitCell[0][2]+rj*HSet.rUnitCell[1][2]+rk*HSet.rUnitCell[2][2];
            if((DistanceSqr(AtomPos[i],r)==1)&&((AtomPos[i][0]-r[0])*HSet.polarization[0]+(AtomPos[i][1]-r[1])*HSet.polarization[1]+(AtomPos[i][2]-r[2])*HSet.polarization[2])!=0)
            {
              H->VCrAnRe[H->NCrAn]=HSet.t*((AtomPos[i][0]-r[0])*HSet.polarization[0]+(AtomPos[i][1]-r[1])*HSet.polarization[1]+(AtomPos[i][2]-r[2])*HSet.polarization[2]);
              H->iCrAn[H->NCrAn]=i;
              H->jCrAn[H->NCrAn]=j;
              H->NCrAn++;
              H->VCrAnRe[H->NCrAn]=HSet.t*((AtomPos[i][0]-r[0])*HSet.polarization[0]+(AtomPos[i][1]-r[1])*HSet.polarization[1]+(AtomPos[i][2]-r[2])*HSet.polarization[2]);
              H->iCrAn[H->NCrAn]=i+H->NFermion;
              H->jCrAn[H->NCrAn]=j+H->NFermion;
              H->NCrAn++;
              H->VCrAnRe[H->NCrAn]=HSet.t*((r[0]-AtomPos[i][0])*HSet.polarization[0]+(r[1]-AtomPos[i][1])*HSet.polarization[1]+(r[2]-AtomPos[i][2])*HSet.polarization[2]);
              H->iCrAn[H->NCrAn]=j;
              H->jCrAn[H->NCrAn]=i;
              H->NCrAn++;
              H->VCrAnRe[H->NCrAn]=HSet.t*((r[0]-AtomPos[i][0])*HSet.polarization[0]+(r[1]-AtomPos[i][1])*HSet.polarization[1]+(r[2]-AtomPos[i][2])*HSet.polarization[2]);
              H->iCrAn[H->NCrAn]=j+H->NFermion;
              H->jCrAn[H->NCrAn]=i+H->NFermion;
              H->NCrAn++;
            }
          }
        }
      }
    }
  }
  H->VCrAnRe=realloc(H->VCrAnRe,H->NCrAn*sizeof(double));
  H->iCrAn=realloc(H->iCrAn,H->NCrAn*sizeof(uint32_t));
  H->jCrAn=realloc(H->jCrAn,H->NCrAn*sizeof(uint32_t));
  H->NFermion*=2;
  H->NWord=Min(H->NFermion,1+(H->NFermion-1)/8);
  printf("\nHubbard Model Optical transition matrix created\n");
  printf("The one particle part of T_{Opt}:\n");
  for(i=0;i<H->NCrAn;i++)
  {
    printf("%4i %4i %10.6f\n",H->iCrAn[i],H->jCrAn[i],H->VCrAnRe[i]);
  }
  printf("\n");
}
*/

/*
void CreateSimpelTwistedHubbardModelOpticsT(const HubbardModelSettings HSet, OperatorType *const H, const uint32_t iq[])
{
  InitOperator(H);
  int32_t AtomPos[100][3];
  int32_t i,j,k,ri,rj,rk;
  double   a,b,c,r[3],g[3][3],rUnitCell[3][3],phase;
  H->NFermion=0;
  // Create reciprocal lattice vectors //
  rUnitCell[0][0]=(double)HSet.rUnitCell[0][0];
  rUnitCell[0][1]=(double)HSet.rUnitCell[0][1];
  rUnitCell[0][2]=(double)HSet.rUnitCell[0][2];
  rUnitCell[1][0]=(double)HSet.rUnitCell[1][0];
  rUnitCell[1][1]=(double)HSet.rUnitCell[1][1];
  rUnitCell[1][2]=(double)HSet.rUnitCell[1][2];
  rUnitCell[2][0]=(double)HSet.rUnitCell[2][0];
  rUnitCell[2][1]=(double)HSet.rUnitCell[2][1];
  rUnitCell[2][2]=(double)HSet.rUnitCell[2][2];
  ResiprocalLatticeVector(rUnitCell, g);
  for(i=-(abs(HSet.rUnitCell[0][0])+abs(HSet.rUnitCell[1][0])+abs(HSet.rUnitCell[2][0]));i<(abs(HSet.rUnitCell[0][0])+abs(HSet.rUnitCell[1][0])+abs(HSet.rUnitCell[2][0]));i++)
  {
    for(j=-(abs(HSet.rUnitCell[0][1])+abs(HSet.rUnitCell[1][1])+abs(HSet.rUnitCell[2][1]));j<(abs(HSet.rUnitCell[0][1])+abs(HSet.rUnitCell[1][1])+abs(HSet.rUnitCell[2][1]));j++)
    {
      for(k=-(abs(HSet.rUnitCell[0][2])+abs(HSet.rUnitCell[1][2])+abs(HSet.rUnitCell[2][2]));k<(abs(HSet.rUnitCell[0][2])+abs(HSet.rUnitCell[1][2])+abs(HSet.rUnitCell[2][2]));k++)
      {
        a=-((double)(HSet.rUnitCell[1][2]*HSet.rUnitCell[2][1]*i - HSet.rUnitCell[1][1]*HSet.rUnitCell[2][2]*i - 
          HSet.rUnitCell[0][2]*HSet.rUnitCell[2][1]*j + HSet.rUnitCell[0][1]*HSet.rUnitCell[2][2]*j + HSet.rUnitCell[0][2]*HSet.rUnitCell[1][1]*k - 
          HSet.rUnitCell[0][1]*HSet.rUnitCell[1][2]*k)/(double)
        (-(HSet.rUnitCell[0][2]*HSet.rUnitCell[1][1]*HSet.rUnitCell[2][0]) + HSet.rUnitCell[0][1]*HSet.rUnitCell[1][2]*HSet.rUnitCell[2][0] + 
          HSet.rUnitCell[0][2]*HSet.rUnitCell[1][0]*HSet.rUnitCell[2][1] - HSet.rUnitCell[0][0]*HSet.rUnitCell[1][2]*HSet.rUnitCell[2][1] - 
          HSet.rUnitCell[0][1]*HSet.rUnitCell[1][0]*HSet.rUnitCell[2][2] + HSet.rUnitCell[0][0]*HSet.rUnitCell[1][1]*HSet.rUnitCell[2][2]));
        b=-((double)(HSet.rUnitCell[1][2]*HSet.rUnitCell[2][0]*i - HSet.rUnitCell[1][0]*HSet.rUnitCell[2][2]*i - 
          HSet.rUnitCell[0][2]*HSet.rUnitCell[2][0]*j + HSet.rUnitCell[0][0]*HSet.rUnitCell[2][2]*j + HSet.rUnitCell[0][2]*HSet.rUnitCell[1][0]*k - 
          HSet.rUnitCell[0][0]*HSet.rUnitCell[1][2]*k)/(double)
        (HSet.rUnitCell[0][2]*HSet.rUnitCell[1][1]*HSet.rUnitCell[2][0] - HSet.rUnitCell[0][1]*HSet.rUnitCell[1][2]*HSet.rUnitCell[2][0] - 
          HSet.rUnitCell[0][2]*HSet.rUnitCell[1][0]*HSet.rUnitCell[2][1] + HSet.rUnitCell[0][0]*HSet.rUnitCell[1][2]*HSet.rUnitCell[2][1] + 
          HSet.rUnitCell[0][1]*HSet.rUnitCell[1][0]*HSet.rUnitCell[2][2] - HSet.rUnitCell[0][0]*HSet.rUnitCell[1][1]*HSet.rUnitCell[2][2]));
        c=-((double)(-(HSet.rUnitCell[1][1]*HSet.rUnitCell[2][0]*i) + HSet.rUnitCell[1][0]*HSet.rUnitCell[2][1]*i + 
          HSet.rUnitCell[0][1]*HSet.rUnitCell[2][0]*j - HSet.rUnitCell[0][0]*HSet.rUnitCell[2][1]*j - HSet.rUnitCell[0][1]*HSet.rUnitCell[1][0]*k + 
          HSet.rUnitCell[0][0]*HSet.rUnitCell[1][1]*k)/(double)
        (HSet.rUnitCell[0][2]*HSet.rUnitCell[1][1]*HSet.rUnitCell[2][0] - HSet.rUnitCell[0][1]*HSet.rUnitCell[1][2]*HSet.rUnitCell[2][0] - 
          HSet.rUnitCell[0][2]*HSet.rUnitCell[1][0]*HSet.rUnitCell[2][1] + HSet.rUnitCell[0][0]*HSet.rUnitCell[1][2]*HSet.rUnitCell[2][1] + 
          HSet.rUnitCell[0][1]*HSet.rUnitCell[1][0]*HSet.rUnitCell[2][2] - HSet.rUnitCell[0][0]*HSet.rUnitCell[1][1]*HSet.rUnitCell[2][2]));      
        if((a>=0)&&(a<1.0)&&(b>=0)&&(b<1.0)&&(c>=0)&&(c<1.0))
        {
          AtomPos[H->NFermion][0]=i;
          AtomPos[H->NFermion][1]=j;
          AtomPos[H->NFermion][2]=k;
          H->NFermion++;
        }
      }
    }
  }
  H->NCrAn=0;
  H->VCrAnRe=malloc(12*H->NFermion*sizeof(double));
  H->VCrAnIm=malloc(12*H->NFermion*sizeof(double));
  H->iCrAn=malloc(12*H->NFermion*sizeof(uint32_t));
  H->jCrAn=malloc(12*H->NFermion*sizeof(uint32_t));
  for(i=0;i<H->NFermion;i++)
  {
    for(j=0;j<i;j++)
    {
      for(ri=-1;ri<=1;ri++)
      {
        for(rj=-1;rj<=1;rj++)
        {
          for(rk=-1;rk<=1;rk++)
          {
            r[0]=AtomPos[j][0]+ri*HSet.rUnitCell[0][0]+rj*HSet.rUnitCell[1][0]+rk*HSet.rUnitCell[2][0];
            r[1]=AtomPos[j][1]+ri*HSet.rUnitCell[0][1]+rj*HSet.rUnitCell[1][1]+rk*HSet.rUnitCell[2][1];
            r[2]=AtomPos[j][2]+ri*HSet.rUnitCell[0][2]+rj*HSet.rUnitCell[1][2]+rk*HSet.rUnitCell[2][2];
            if((DistanceSqr(AtomPos[i],r)==1)&&((AtomPos[i][0]-r[0])*HSet.polarization[0]+(AtomPos[i][1]-r[1])*HSet.polarization[1]+(AtomPos[i][2]-r[2])*HSet.polarization[2])!=0)
            {
              r[0]-=AtomPos[i][0];
              r[1]-=AtomPos[i][1];
              r[2]-=AtomPos[i][2];
              phase=Dot(g[0],r)*iq[0]/((double)HSet.Nq[0])+Dot(g[1],r)*iq[1]/((double)HSet.Nq[1])+Dot(g[2],r)*iq[2]/((double)HSet.Nq[2]);
              r[0]+=AtomPos[i][0];
              r[1]+=AtomPos[i][1];
              r[2]+=AtomPos[i][2];
              H->VCrAnRe[H->NCrAn]=-sin(phase)*HSet.t*((AtomPos[i][0]-r[0])*HSet.polarization[0]+(AtomPos[i][1]-r[1])*HSet.polarization[1]+(AtomPos[i][2]-r[2])*HSet.polarization[2]);
              H->VCrAnIm[H->NCrAn]= cos(phase)*HSet.t*((AtomPos[i][0]-r[0])*HSet.polarization[0]+(AtomPos[i][1]-r[1])*HSet.polarization[1]+(AtomPos[i][2]-r[2])*HSet.polarization[2]);
              H->iCrAn[H->NCrAn]=i;
              H->jCrAn[H->NCrAn]=j;
              H->NCrAn++;
              H->VCrAnRe[H->NCrAn]=-sin(phase)*HSet.t*((AtomPos[i][0]-r[0])*HSet.polarization[0]+(AtomPos[i][1]-r[1])*HSet.polarization[1]+(AtomPos[i][2]-r[2])*HSet.polarization[2]);
              H->VCrAnIm[H->NCrAn]= cos(phase)*HSet.t*((AtomPos[i][0]-r[0])*HSet.polarization[0]+(AtomPos[i][1]-r[1])*HSet.polarization[1]+(AtomPos[i][2]-r[2])*HSet.polarization[2]);
              H->iCrAn[H->NCrAn]=i+H->NFermion;
              H->jCrAn[H->NCrAn]=j+H->NFermion;
              H->NCrAn++;
              H->VCrAnRe[H->NCrAn]=-sin(phase)*HSet.t*((r[0]-AtomPos[i][0])*HSet.polarization[0]+(r[1]-AtomPos[i][1])*HSet.polarization[1]+(r[2]-AtomPos[i][2])*HSet.polarization[2]);
              H->VCrAnIm[H->NCrAn]= cos(phase)*HSet.t*((r[0]-AtomPos[i][0])*HSet.polarization[0]+(r[1]-AtomPos[i][1])*HSet.polarization[1]+(r[2]-AtomPos[i][2])*HSet.polarization[2]);
              H->iCrAn[H->NCrAn]=j;
              H->jCrAn[H->NCrAn]=i;
              H->NCrAn++;
              H->VCrAnRe[H->NCrAn]=-sin(phase)*HSet.t*((r[0]-AtomPos[i][0])*HSet.polarization[0]+(r[1]-AtomPos[i][1])*HSet.polarization[1]+(r[2]-AtomPos[i][2])*HSet.polarization[2]);
              H->VCrAnIm[H->NCrAn]= cos(phase)*HSet.t*((r[0]-AtomPos[i][0])*HSet.polarization[0]+(r[1]-AtomPos[i][1])*HSet.polarization[1]+(r[2]-AtomPos[i][2])*HSet.polarization[2]);
              H->iCrAn[H->NCrAn]=j+H->NFermion;
              H->jCrAn[H->NCrAn]=i+H->NFermion;
              H->NCrAn++;
            }
          }
        }
      }
    }
  }
  H->VCrAnRe=realloc(H->VCrAnRe,H->NCrAn*sizeof(double));
  H->VCrAnIm=realloc(H->VCrAnIm,H->NCrAn*sizeof(double));
  H->iCrAn=realloc(H->iCrAn,H->NCrAn*sizeof(uint32_t));
  H->jCrAn=realloc(H->jCrAn,H->NCrAn*sizeof(uint32_t));
  H->NFermion*=2;
  H->NWord=Min(H->NFermion,1+(H->NFermion-1)/8);
  printf("\nHubbard Model Optical transition matrix created\n");
  printf("The one particle part of T_{Opt}:\n");
  for(i=0;i<H->NCrAn;i++)
  {
    if(H->VCrAnIm[i]<0)
    {
      printf("%4i %4i %10.6f - I %10.6f\n",H->iCrAn[i],H->jCrAn[i],H->VCrAnRe[i],-H->VCrAnIm[i]);
    }
    else
    {
      printf("%4i %4i %10.6f + I %10.6f\n",H->iCrAn[i],H->jCrAn[i],H->VCrAnRe[i],H->VCrAnIm[i]);
    }
  }
  printf("\n");
}
*/


/*
WaveFunctionType CreateStartWaveFunitionNeelState(const HubbardModelSettings HSet, const OperatorType H)
{
  int32_t AtomPos[100][3];
  int32_t i,j,k,iorb,Nup,Ndn;
  double   a,b,c;
  WaveFunctionType psi,Opsi;
  uint8_t *Determinant;
  OperatorType T;
  iorb=0;

  for(i=-(abs(HSet.rUnitCell[0][0])+abs(HSet.rUnitCell[1][0])+abs(HSet.rUnitCell[2][0]));i<(abs(HSet.rUnitCell[0][0])+abs(HSet.rUnitCell[1][0])+abs(HSet.rUnitCell[2][0]));i++)
  {
    for(j=-(abs(HSet.rUnitCell[0][1])+abs(HSet.rUnitCell[1][1])+abs(HSet.rUnitCell[2][1]));j<(abs(HSet.rUnitCell[0][1])+abs(HSet.rUnitCell[1][1])+abs(HSet.rUnitCell[2][1]));j++)
    {
      for(k=-(abs(HSet.rUnitCell[0][2])+abs(HSet.rUnitCell[1][2])+abs(HSet.rUnitCell[2][2]));k<(abs(HSet.rUnitCell[0][2])+abs(HSet.rUnitCell[1][2])+abs(HSet.rUnitCell[2][2]));k++)
      {
        a=-((double)(HSet.rUnitCell[1][2]*HSet.rUnitCell[2][1]*i - HSet.rUnitCell[1][1]*HSet.rUnitCell[2][2]*i - 
          HSet.rUnitCell[0][2]*HSet.rUnitCell[2][1]*j + HSet.rUnitCell[0][1]*HSet.rUnitCell[2][2]*j + HSet.rUnitCell[0][2]*HSet.rUnitCell[1][1]*k - 
          HSet.rUnitCell[0][1]*HSet.rUnitCell[1][2]*k)/(double)
        (-(HSet.rUnitCell[0][2]*HSet.rUnitCell[1][1]*HSet.rUnitCell[2][0]) + HSet.rUnitCell[0][1]*HSet.rUnitCell[1][2]*HSet.rUnitCell[2][0] + 
          HSet.rUnitCell[0][2]*HSet.rUnitCell[1][0]*HSet.rUnitCell[2][1] - HSet.rUnitCell[0][0]*HSet.rUnitCell[1][2]*HSet.rUnitCell[2][1] - 
          HSet.rUnitCell[0][1]*HSet.rUnitCell[1][0]*HSet.rUnitCell[2][2] + HSet.rUnitCell[0][0]*HSet.rUnitCell[1][1]*HSet.rUnitCell[2][2]));
        b=-((double)(HSet.rUnitCell[1][2]*HSet.rUnitCell[2][0]*i - HSet.rUnitCell[1][0]*HSet.rUnitCell[2][2]*i - 
          HSet.rUnitCell[0][2]*HSet.rUnitCell[2][0]*j + HSet.rUnitCell[0][0]*HSet.rUnitCell[2][2]*j + HSet.rUnitCell[0][2]*HSet.rUnitCell[1][0]*k - 
          HSet.rUnitCell[0][0]*HSet.rUnitCell[1][2]*k)/(double)
        (HSet.rUnitCell[0][2]*HSet.rUnitCell[1][1]*HSet.rUnitCell[2][0] - HSet.rUnitCell[0][1]*HSet.rUnitCell[1][2]*HSet.rUnitCell[2][0] - 
          HSet.rUnitCell[0][2]*HSet.rUnitCell[1][0]*HSet.rUnitCell[2][1] + HSet.rUnitCell[0][0]*HSet.rUnitCell[1][2]*HSet.rUnitCell[2][1] + 
          HSet.rUnitCell[0][1]*HSet.rUnitCell[1][0]*HSet.rUnitCell[2][2] - HSet.rUnitCell[0][0]*HSet.rUnitCell[1][1]*HSet.rUnitCell[2][2]));
        c=-((double)(-(HSet.rUnitCell[1][1]*HSet.rUnitCell[2][0]*i) + HSet.rUnitCell[1][0]*HSet.rUnitCell[2][1]*i + 
          HSet.rUnitCell[0][1]*HSet.rUnitCell[2][0]*j - HSet.rUnitCell[0][0]*HSet.rUnitCell[2][1]*j - HSet.rUnitCell[0][1]*HSet.rUnitCell[1][0]*k + 
          HSet.rUnitCell[0][0]*HSet.rUnitCell[1][1]*k)/(double)
        (HSet.rUnitCell[0][2]*HSet.rUnitCell[1][1]*HSet.rUnitCell[2][0] - HSet.rUnitCell[0][1]*HSet.rUnitCell[1][2]*HSet.rUnitCell[2][0] - 
          HSet.rUnitCell[0][2]*HSet.rUnitCell[1][0]*HSet.rUnitCell[2][1] + HSet.rUnitCell[0][0]*HSet.rUnitCell[1][2]*HSet.rUnitCell[2][1] + 
          HSet.rUnitCell[0][1]*HSet.rUnitCell[1][0]*HSet.rUnitCell[2][2] - HSet.rUnitCell[0][0]*HSet.rUnitCell[1][1]*HSet.rUnitCell[2][2]));      
        if((a>=0)&&(a<1.0)&&(b>=0)&&(b<1.0)&&(c>=0)&&(c<1.0))
        {
          AtomPos[iorb][0]=i;
          AtomPos[iorb][1]=j;
          AtomPos[iorb][2]=k;
          iorb++;
        }
      }
    }
  }
  psi.NBitsKey=HSet.NBitsKey;
  psi.NFermion=H.NFermion;
  psi.NWord=H.NWord;
  psi.N=0;
  psi.NBlok=1;
  psi.KeyList=calloc((uint32_t)1<<(psi.NBitsKey),sizeof(uint32_t));
  psi.ValueRe=malloc(psi.NBlok*sizeof(double *));
  psi.Next=malloc(psi.NBlok*sizeof(uint32_t *));
  psi.Determinant=malloc(psi.NBlok*sizeof(uint8_t *));
  if((psi.ValueRe==NULL)||(psi.Next==NULL)||(psi.Determinant==NULL)) AllocFailed();
  psi.ValueRe[psi.NBlok-1]=malloc(BlokSize*sizeof(double));
  psi.Next[psi.NBlok-1]=malloc(BlokSize*sizeof(uint32_t));
  psi.Determinant[psi.NBlok-1]=malloc(BlokSize*psi.NWord*sizeof(uint8_t));
  if((psi.ValueRe[psi.NBlok-1]==NULL)||(psi.Next[psi.NBlok-1]==NULL)||(psi.Determinant[psi.NBlok-1]==NULL)) AllocFailed();
  Determinant=calloc(psi.NWord,sizeof(uint8_t));
  // Step one create a half filled Neel state //
  Nup=0;
  Ndn=0;
  for(i=0;i<H.NFermion/2;i++)
  {
    if((AtomPos[i][0]+AtomPos[i][1]+AtomPos[i][2])&1==1)
    {
      OperatorCr(Determinant,i,H.NFermion);
      Nup++;
    }
    else
    {
      OperatorCr(Determinant,i+H.NFermion/2,H.NFermion);
      Ndn++;
    }
  }
  RealWaveFunctionAddElement(&psi, Determinant, 1.0);
  // Remove electrons from spin up states if HSet.NelecUp<Nup //
  T.NCr=0;T.NAn=0;T.NCrAn=0;T.NCrCrAnAn=0;
  if(HSet.NelecUp<Nup)
  {
    T.NAn=H.NFermion/2.0;
    T.iAn=malloc(T.NAn*sizeof(uint32_t));
    T.VAnRe=malloc(T.NAn*sizeof(double));
    for(i=0;i<T.NAn;i++)
    {
      T.iAn[i]=i;
      T.VAnRe[i]=1.0;
    }
    while(HSet.NelecUp<Nup)
    {
      Nup--;
      OperatorPsi(T,psi,&Opsi);
      RealWaveFunctionFree(&psi);
      psi=Opsi;
    }
    FreeOperator(&T);
  }
  // Remove electrons from spin dn states if HSet.NelecDn<Ndn //
  if(HSet.NelecDn<Ndn)
  {
    T.NAn=H.NFermion/2;
    T.iAn=malloc(T.NAn*sizeof(uint32_t));
    T.VAnRe=malloc(T.NAn*sizeof(double));
    for(i=0;i<T.NAn;i++)
    {
      T.iAn[i]=i+H.NFermion/2;
      T.VAnRe[i]=1.0;
    }
    while(HSet.NelecDn<Ndn)
    {
      Ndn--;
      OperatorPsi(T,psi,&Opsi);
      RealWaveFunctionFree(&psi);
      psi=Opsi;
    }
    FreeOperator(&T);
  }
  // Add electrons to spin up states if HSet.NelecUp>Nup //
  if(HSet.NelecUp>Nup)
  {
    T.NCr=H.NFermion/2.0;
    T.iCr=malloc(T.NCr*sizeof(uint32_t));
    T.VCrRe=malloc(T.NCr*sizeof(double));
    for(i=0;i<T.NCr;i++)
    {
      T.iCr[i]=i;
      T.VCrRe[i]=1.0;
    }
    while(HSet.NelecUp>Nup)
    {
      Nup++;
      OperatorPsi(T,psi,&Opsi);
      RealWaveFunctionFree(&psi);
      psi=Opsi;
    }
    FreeOperator(&T);
  }
  // Add electrons to spin dn states if HSet.NelecDn>Ndn //
  if(HSet.NelecDn>Ndn)
  {
    T.NCr=H.NFermion/2;
    T.iCr=malloc(T.NCr*sizeof(uint32_t));
    T.VCrRe=malloc(T.NCr*sizeof(double));
    for(i=0;i<T.NCr;i++)
    {
      T.iCr[i]=i+H.NFermion/2;
      T.VCrRe[i]=1.0;
    }
    while(HSet.NelecDn>Ndn)
    {
      Ndn++;
      OperatorPsi(T,psi,&Opsi);
      RealWaveFunctionFree(&psi);
      psi=Opsi;
    }
    FreeOperator(&T);
  }
  RealWaveFunctionRandomize(&psi);
  printf("Starting Wave Function Created\n\n");
  return psi;
}
*/
/*
WaveFunctionType CreateStartComplexWaveFunitionNeelState(const HubbardModelSettings HSet)
{
  int32_t AtomPos[100][3];
  int32_t i,j,k,iorb,Nup,Ndn;
  double   a,b,c;
  WaveFunctionType psi,Opsi;
  uint8_t *Determinant;
  OperatorType T;
  iorb=0;

  for(i=-(abs(HSet.rUnitCell[0][0])+abs(HSet.rUnitCell[1][0])+abs(HSet.rUnitCell[2][0]));i<(abs(HSet.rUnitCell[0][0])+abs(HSet.rUnitCell[1][0])+abs(HSet.rUnitCell[2][0]));i++)
  {
    for(j=-(abs(HSet.rUnitCell[0][1])+abs(HSet.rUnitCell[1][1])+abs(HSet.rUnitCell[2][1]));j<(abs(HSet.rUnitCell[0][1])+abs(HSet.rUnitCell[1][1])+abs(HSet.rUnitCell[2][1]));j++)
    {
      for(k=-(abs(HSet.rUnitCell[0][2])+abs(HSet.rUnitCell[1][2])+abs(HSet.rUnitCell[2][2]));k<(abs(HSet.rUnitCell[0][2])+abs(HSet.rUnitCell[1][2])+abs(HSet.rUnitCell[2][2]));k++)
      {
        a=-((double)(HSet.rUnitCell[1][2]*HSet.rUnitCell[2][1]*i - HSet.rUnitCell[1][1]*HSet.rUnitCell[2][2]*i - 
          HSet.rUnitCell[0][2]*HSet.rUnitCell[2][1]*j + HSet.rUnitCell[0][1]*HSet.rUnitCell[2][2]*j + HSet.rUnitCell[0][2]*HSet.rUnitCell[1][1]*k - 
          HSet.rUnitCell[0][1]*HSet.rUnitCell[1][2]*k)/(double)
        (-(HSet.rUnitCell[0][2]*HSet.rUnitCell[1][1]*HSet.rUnitCell[2][0]) + HSet.rUnitCell[0][1]*HSet.rUnitCell[1][2]*HSet.rUnitCell[2][0] + 
          HSet.rUnitCell[0][2]*HSet.rUnitCell[1][0]*HSet.rUnitCell[2][1] - HSet.rUnitCell[0][0]*HSet.rUnitCell[1][2]*HSet.rUnitCell[2][1] - 
          HSet.rUnitCell[0][1]*HSet.rUnitCell[1][0]*HSet.rUnitCell[2][2] + HSet.rUnitCell[0][0]*HSet.rUnitCell[1][1]*HSet.rUnitCell[2][2]));
        b=-((double)(HSet.rUnitCell[1][2]*HSet.rUnitCell[2][0]*i - HSet.rUnitCell[1][0]*HSet.rUnitCell[2][2]*i - 
          HSet.rUnitCell[0][2]*HSet.rUnitCell[2][0]*j + HSet.rUnitCell[0][0]*HSet.rUnitCell[2][2]*j + HSet.rUnitCell[0][2]*HSet.rUnitCell[1][0]*k - 
          HSet.rUnitCell[0][0]*HSet.rUnitCell[1][2]*k)/(double)
        (HSet.rUnitCell[0][2]*HSet.rUnitCell[1][1]*HSet.rUnitCell[2][0] - HSet.rUnitCell[0][1]*HSet.rUnitCell[1][2]*HSet.rUnitCell[2][0] - 
          HSet.rUnitCell[0][2]*HSet.rUnitCell[1][0]*HSet.rUnitCell[2][1] + HSet.rUnitCell[0][0]*HSet.rUnitCell[1][2]*HSet.rUnitCell[2][1] + 
          HSet.rUnitCell[0][1]*HSet.rUnitCell[1][0]*HSet.rUnitCell[2][2] - HSet.rUnitCell[0][0]*HSet.rUnitCell[1][1]*HSet.rUnitCell[2][2]));
        c=-((double)(-(HSet.rUnitCell[1][1]*HSet.rUnitCell[2][0]*i) + HSet.rUnitCell[1][0]*HSet.rUnitCell[2][1]*i + 
          HSet.rUnitCell[0][1]*HSet.rUnitCell[2][0]*j - HSet.rUnitCell[0][0]*HSet.rUnitCell[2][1]*j - HSet.rUnitCell[0][1]*HSet.rUnitCell[1][0]*k + 
          HSet.rUnitCell[0][0]*HSet.rUnitCell[1][1]*k)/(double)
        (HSet.rUnitCell[0][2]*HSet.rUnitCell[1][1]*HSet.rUnitCell[2][0] - HSet.rUnitCell[0][1]*HSet.rUnitCell[1][2]*HSet.rUnitCell[2][0] - 
          HSet.rUnitCell[0][2]*HSet.rUnitCell[1][0]*HSet.rUnitCell[2][1] + HSet.rUnitCell[0][0]*HSet.rUnitCell[1][2]*HSet.rUnitCell[2][1] + 
          HSet.rUnitCell[0][1]*HSet.rUnitCell[1][0]*HSet.rUnitCell[2][2] - HSet.rUnitCell[0][0]*HSet.rUnitCell[1][1]*HSet.rUnitCell[2][2]));      
        if((a>=0)&&(a<1.0)&&(b>=0)&&(b<1.0)&&(c>=0)&&(c<1.0))
        {
          AtomPos[iorb][0]=i;
          AtomPos[iorb][1]=j;
          AtomPos[iorb][2]=k;
          iorb++;
        }
      }
    }
  }
  psi.NBitsKey=HSet.NBitsKey;
  psi.NFermion=2*iorb;
  psi.NWord=Min(psi.NFermion,1+(psi.NFermion-1)/8);
  psi.N=0;
  psi.NBlok=1;
  psi.KeyList=calloc((uint32_t)1<<(psi.NBitsKey),sizeof(uint32_t));
  psi.ValueRe=malloc(psi.NBlok*sizeof(double *));
  psi.ValueIm=malloc(psi.NBlok*sizeof(double *));
  psi.Next=malloc(psi.NBlok*sizeof(uint32_t *));
  psi.Determinant=malloc(psi.NBlok*sizeof(uint8_t *));
  if((psi.ValueRe==NULL)||(psi.ValueIm==NULL)||(psi.Next==NULL)||(psi.Determinant==NULL)) AllocFailed();
  psi.ValueRe[psi.NBlok-1]=malloc(BlokSize*sizeof(double));
  psi.ValueIm[psi.NBlok-1]=malloc(BlokSize*sizeof(double));
  psi.Next[psi.NBlok-1]=malloc(BlokSize*sizeof(uint32_t));
  psi.Determinant[psi.NBlok-1]=malloc(BlokSize*psi.NWord*sizeof(uint8_t));
  if((psi.ValueRe[psi.NBlok-1]==NULL)||(psi.ValueIm[psi.NBlok-1]==NULL)||(psi.Next[psi.NBlok-1]==NULL)||(psi.Determinant[psi.NBlok-1]==NULL)) AllocFailed();
  Determinant=calloc(psi.NWord,sizeof(uint8_t));
  // Step one create a half filled Neel state //
  Nup=0;
  Ndn=0;
  for(i=0;i<psi.NFermion/2;i++)
  {
    if((AtomPos[i][0]+AtomPos[i][1]+AtomPos[i][2])&1==1)
    {
      OperatorCr(Determinant,i,psi.NFermion);
      Nup++;
    }
    else
    {
      OperatorCr(Determinant,i+psi.NFermion/2,psi.NFermion);
      Ndn++;
    }
  }
  ComplexWaveFunctionAddElement(&psi, Determinant, 1.0, 0.0);
  // Remove electrons from spin up states if HSet.NelecUp<Nup //
  T.NCr=0;T.NAn=0;T.NCrAn=0;T.NCrCrAnAn=0;
  if(HSet.NelecUp<Nup)
  {
    T.NAn=psi.NFermion/2.0;
    T.iAn=malloc(T.NAn*sizeof(uint32_t));
    T.VAnRe=malloc(T.NAn*sizeof(double));
    T.VAnIm=malloc(T.NAn*sizeof(double));
    for(i=0;i<T.NAn;i++)
    {
      T.iAn[i]=i;
      T.VAnRe[i]=1.0;
      T.VAnIm[i]=0.0;
    }
    while(HSet.NelecUp<Nup)
    {
      Nup--;
      ComplexOperatorPsi(T,psi,&Opsi);
      ComplexWaveFunctionFree(&psi);
      psi=Opsi;
    }
    FreeOperator(&T);
  }
  // Remove electrons from spin dn states if HSet.NelecDn<Ndn //
  if(HSet.NelecDn<Ndn)
  {
    T.NAn=psi.NFermion/2;
    T.iAn=malloc(T.NAn*sizeof(uint32_t));
    T.VAnRe=malloc(T.NAn*sizeof(double));
    T.VAnIm=malloc(T.NAn*sizeof(double));
    for(i=0;i<T.NAn;i++)
    {
      T.iAn[i]=i+psi.NFermion/2;
      T.VAnRe[i]=1.0;
      T.VAnIm[i]=0.0;
    }
    while(HSet.NelecDn<Ndn)
    {
      Ndn--;
      ComplexOperatorPsi(T,psi,&Opsi);
      ComplexWaveFunctionFree(&psi);
      psi=Opsi;
    }
    FreeOperator(&T);
  }
  // Add electrons to spin up states if HSet.NelecUp>Nup //
  if(HSet.NelecUp>Nup)
  {
    T.NCr=psi.NFermion/2.0;
    T.iCr=malloc(T.NCr*sizeof(uint32_t));
    T.VCrRe=malloc(T.NCr*sizeof(double));
    T.VCrIm=malloc(T.NCr*sizeof(double));
    for(i=0;i<T.NCr;i++)
    {
      T.iCr[i]=i;
      T.VCrRe[i]=1.0;
      T.VCrIm[i]=0.0;
    }
    while(HSet.NelecUp>Nup)
    {
      Nup++;
      ComplexOperatorPsi(T,psi,&Opsi);
      ComplexWaveFunctionFree(&psi);
      psi=Opsi;
    }
    FreeOperator(&T);
  }
  // Add electrons to spin dn states if HSet.NelecDn>Ndn //
  if(HSet.NelecDn>Ndn)
  {
    T.NCr=psi.NFermion/2;
    T.iCr=malloc(T.NCr*sizeof(uint32_t));
    T.VCrRe=malloc(T.NCr*sizeof(double));
    T.VCrIm=malloc(T.NCr*sizeof(double));
    for(i=0;i<T.NCr;i++)
    {
      T.iCr[i]=i+psi.NFermion/2;
      T.VCrRe[i]=1.0;
      T.VCrIm[i]=0.0;
    }
    while(HSet.NelecDn>Ndn)
    {
      Ndn++;
      ComplexOperatorPsi(T,psi,&Opsi);
      ComplexWaveFunctionFree(&psi);
      psi=Opsi;
    }
    FreeOperator(&T);
  }
  ComplexWaveFunctionRandomize(&psi);
  printf("Starting Wave Function Created with %i elements\n\n",psi.N);
  return psi;
}
*/
/*
WaveFunctionType CreateStartWaveFunctionSimpelHubbardModel(const HubbardModelSettings HSet, const OperatorType H)
{
  WaveFunctionType psi,Opsi;
  uint8_t *Determinant;
  uint32_t i,phase;
  psi.NBitsKey=HSet.NBitsKey;
  psi.NFermion=H.NFermion;
  psi.NWord=H.NWord;
  psi.N=0;
  psi.NBlok=1;
  psi.KeyList=calloc((uint32_t)1<<(psi.NBitsKey),sizeof(uint32_t));
  psi.ValueRe=malloc(psi.NBlok*sizeof(double *));
  psi.Next=malloc(psi.NBlok*sizeof(uint32_t *));
  psi.Determinant=malloc(psi.NBlok*sizeof(uint8_t *));
  if((psi.ValueRe==NULL)||(psi.Next==NULL)||(psi.Determinant==NULL)) AllocFailed();
  psi.ValueRe[psi.NBlok-1]=malloc(BlokSize*sizeof(double));
  psi.Next[psi.NBlok-1]=malloc(BlokSize*sizeof(uint32_t));
  psi.Determinant[psi.NBlok-1]=malloc(BlokSize*psi.NWord*sizeof(uint8_t));
  if((psi.ValueRe[psi.NBlok-1]==NULL)||(psi.Next[psi.NBlok-1]==NULL)||(psi.Determinant[psi.NBlok-1]==NULL)) AllocFailed();
  Determinant=calloc(psi.NWord,sizeof(uint8_t));
  for(i=0;i<HSet.NelecUp;i++) OperatorCr(Determinant,i,psi.NFermion);
  for(i=psi.NFermion/2;i<(psi.NFermion/2+HSet.NelecDn);i++) OperatorCr(Determinant,i,psi.NFermion);
  RealWaveFunctionAddElement(&psi, Determinant, 5.0);
  phase=OperatorCr(Determinant,psi.NFermion-1,psi.NFermion);
  phase*=OperatorAn(Determinant,psi.NFermion/2,psi.NFermion);
  if(phase) RealWaveFunctionAddElement(&psi, Determinant, 3.0);
  free(Determinant);
  Determinant=calloc(psi.NWord,sizeof(uint8_t));
  for(i=0;i<HSet.NelecUp;i++) OperatorCr(Determinant,i,psi.NFermion);
  for(i=psi.NFermion/2;i<(psi.NFermion/2+HSet.NelecDn);i++) OperatorCr(Determinant,i,psi.NFermion);
  phase=OperatorCr(Determinant,psi.NFermion/2-1,psi.NFermion);
  phase*=OperatorAn(Determinant,HSet.NelecUp-1,psi.NFermion);
  if(phase) RealWaveFunctionAddElement(&psi, Determinant, 1.0);
  free(Determinant);
  OperatorPsi(H,psi,&Opsi);
  RealWaveFunctionRandomize(&Opsi);
  RealWaveFunctionNormalize(&Opsi);
  RealWaveFunctionFree(&psi);
  OperatorPsi(H,Opsi,&psi);
  RealWaveFunctionRandomize(&psi);
  RealWaveFunctionNormalize(&psi);
  RealWaveFunctionFree(&Opsi);
  printf("Starting Wave Function Created\n\n");
  return psi;  
}
*/

/*
void mainHubbardModel (HubbardModelSettings HSet) 
{
  WaveFunctionType psi,Opsi,*PsiList;
  OperatorType H,T;
  TriDiagonalMatrixType s;
  SpectraType *SpecPES,*SpecIPES,SpecOpt,*SpecSPINUp,*SpecSPINDn;
  double TotEmax,Energy,*EnergyList;
  uint32_t i,j,jacc,jpsi;
  char flname[64];
  FILE *stream;
  time_t t0,t1,t2,tTotal=0,tLanczosGroundState=0,tLanczosTriDiagonalize=0,tTriDiagonalMatrixToSpectrum=0;
  
  time(&t0);

  s.N=HSet.NTriDiag;
  s.a=malloc(s.N*sizeof(double));
  s.b=malloc(s.N*sizeof(double));
  PsiList=malloc(HSet.NPsi*sizeof(WaveFunctionType));
  EnergyList=malloc(HSet.NPsi*sizeof(double));



  printf("Hubbard Model with only nearest neighbor hopping:\n");

  H=CreateSimpelHubbardModel(HSet);
  fflush(stdout);

  SpecOpt.N=HSet.NE;
  SpecOpt.Emin=HSet.Emin;
  SpecOpt.Emax=HSet.Emax;
  SpecOpt.dE=(HSet.Emax-HSet.Emin)/((double)HSet.NE);
  SpecOpt.EminPole=HSet.Emin;
  SpecOpt.EmaxPole=HSet.Emax;
  SpecOpt.Gamma=HSet.Gamma;

  SpectrumAllocate(SpecOpt,&SpecOpt);
  SpecPES=malloc(H.NFermion/2*sizeof(SpectraType));
  SpecIPES=malloc(H.NFermion/2*sizeof(SpectraType));
  SpecSPINUp=malloc(H.NFermion/2*sizeof(SpectraType));
  SpecSPINDn=malloc(H.NFermion/2*sizeof(SpectraType));
  for(j=0;j<H.NFermion/2;++j)
  {
    SpectrumAllocate(SpecOpt,&SpecPES[j]);
    SpectrumAllocate(SpecOpt,&SpecIPES[j]);
    SpectrumAllocate(SpecOpt,&SpecSPINUp[j]);
    SpectrumAllocate(SpecOpt,&SpecSPINDn[j]);
  }
  for(jpsi=0;jpsi<HSet.NPsi;jpsi++)
  {
    PsiList[jpsi]=CreateStartWaveFunitionNeelState(HSet, H);
  }
  TotEmax=Min(HSet.NelecUp,HSet.NelecDn)*HSet.U+(HSet.NelecUp+HSet.NelecDn)*2*fabs(HSet.t)-0.5*(HSet.NelecUp+HSet.NelecDn)*HSet.U;
  for(jacc=0;jacc<=HSet.Ni;jacc++)
  {
    if(jacc!=0)
    {
      HSet.zeroLancWavefunction/=10;
      HSet.zeroLancEnergy/=10;
      if((HSet.zeroLancEnergy/(TotEmax-EnergyList[0]))<1.0E-14)HSet.zeroLancEnergy*=10;
      printf("\n=================\nIncreased accuracy to %22.15E %22.15E \n=================\n",HSet.zeroLancWavefunction,HSet.zeroLancEnergy);
      fflush(stdout);
    }
    for(jpsi=0;jpsi<HSet.NPsi;jpsi++)
    {
      sprintf(flname,"Opt_PES_IPES_SPIN_%i_%i.dat",jacc,jpsi);
      psi=PsiList[jpsi];
      printf("Starting density\n");
      PrintWaveFunctionDensity(psi);
      printf("Number of determinants used for starting function %i\n",psi.N);
      printf("Get ground-state\n");
      fflush(stdout);
      time(&t1);
      LanczosGroundState(H,&psi,&HSet.zeroLancWavefunction,HSet.zeroLancEnergy,TotEmax,PsiList,jpsi);
      time(&t2);
      tLanczosGroundState+=t2-t1;
      PsiList[jpsi]=psi;
      PrintWaveFunctionDensity(psi);
      printf("Number of determinants used for eigen function %i\n",psi.N);
      printf("Time used for LanczosGroundState %i\n",(int)(t2-t1));
      fflush(stdout);
      if(OperatorPsiConserveBasis(H,psi,&Opsi))AllocFailed();
      Energy=RealWaveFunctionDotProduct(Opsi,psi);
      EnergyList[jpsi]=Energy;
      printf("Ground state energy is %22.15E\n",Energy);
      RealWaveFunctionFreeValues(&Opsi);
      fflush(stdout);
  
      printf("Optical conductivety \n\n");
      CreateSimpelHubbardModelOpticsT(HSet,&T);
      OperatorPsi(T,psi,&Opsi);
      time(&t1);
      LanczosTriDiagonalize(H, Opsi, HSet.zeroLancWavefunction,&s,HSet.Gamma);
      time(&t2);
      tLanczosTriDiagonalize+=t2-t1;
      printf("Time used for LanczosTriDiagonalize %i\n",(int)(t2-t1));
      time(&t1);
      TriDiagonalMatrixToSpectrum(Energy,s,&SpecOpt);
      time(&t2);
      tTriDiagonalMatrixToSpectrum+=t2-t1;
      printf("Time used for TriDiagonalMatrixToSpectrum %i\n",(int)(t2-t1));
      fflush(stdout);
      FreeOperator(&T);
      
      printf("\n\nPES\n\n");
      T.NAn=1;
      T.iAn=malloc(2*sizeof(uint32_t));
      T.VAnRe=malloc(2*sizeof(double));
      T.iAn[0]=0;
      T.VAnRe[0]=1.0;
      OperatorPsi(T,psi,&Opsi);
      time(&t1);
      LanczosTriDiagonalize(H, Opsi, HSet.zeroLancWavefunction,&s,HSet.Gamma);
      time(&t2);
      tLanczosTriDiagonalize+=t2-t1;
      printf("Time used for LanczosTriDiagonalize %i\n",(int)(t2-t1));
      time(&t1);
      TriDiagonalMatrixToSpectrum(Energy,s,&SpecPES[0]);
      time(&t2);
      tTriDiagonalMatrixToSpectrum+=t2-t1;
      printf("Time used for TriDiagonalMatrixToSpectrum %i\n",(int)(t2-t1));
      fflush(stdout);
      T.NAn=2;
      T.VAnRe[0]=sqrt(0.5);
      T.VAnRe[1]=sqrt(0.5);
      for(j=1;j<psi.NFermion/2;++j)
      {
        T.iAn[1]=j;
        OperatorPsi(T,psi,&Opsi);
        LanczosTriDiagonalize(H, Opsi, HSet.zeroLancWavefunction,&s,HSet.Gamma);
        TriDiagonalMatrixToSpectrum(Energy,s,&SpecPES[j]);
      }
      FreeOperator(&T);
  
  
      printf("\n\nIPES\n\n");
      T.NCr=1;
      T.iCr=malloc(2*sizeof(uint32_t));
      T.VCrRe=malloc(2*sizeof(double));
      T.iCr[0]=0;
      T.VCrRe[0]=1.0;
      OperatorPsi(T,psi,&Opsi);
      time(&t1);
      LanczosTriDiagonalize(H, Opsi, HSet.zeroLancWavefunction,&s,HSet.Gamma);
      time(&t2);
      tLanczosTriDiagonalize+=t2-t1;
      printf("Time used for LanczosTriDiagonalize %i\n",(int)(t2-t1));
      time(&t1);
      TriDiagonalMatrixToSpectrum(Energy,s,&SpecIPES[0]);
      time(&t2);
      tTriDiagonalMatrixToSpectrum+=t2-t1;
      printf("Time used for TriDiagonalMatrixToSpectrum %i\n",(int)(t2-t1));
      fflush(stdout);
      T.NCr=2;
      T.VCrRe[0]=sqrt(0.5);
      T.VCrRe[1]=sqrt(0.5);
      for(j=1;j<psi.NFermion/2;++j)
      {
        T.iCr[1]=j;
        OperatorPsi(T,psi,&Opsi);
        LanczosTriDiagonalize(H, Opsi, HSet.zeroLancWavefunction,&s,HSet.Gamma);
        TriDiagonalMatrixToSpectrum(Energy,s,&SpecIPES[j]);
      }
      FreeOperator(&T);
  
      printf("\n\nSPIN Up\n\n");
      T.NCrAn=1;
      T.iCrAn=malloc(2*sizeof(uint32_t));
      T.jCrAn=malloc(2*sizeof(uint32_t));
      T.VCrAnRe=malloc(2*sizeof(double));
      T.iCrAn[0]=0;
      T.jCrAn[0]=H.NFermion/2;
      T.VCrAnRe[0]=1.0;
      OperatorPsi(T,psi,&Opsi);
      time(&t1);
      LanczosTriDiagonalize(H, Opsi, HSet.zeroLancWavefunction,&s,HSet.Gamma);
      time(&t2);
      tLanczosTriDiagonalize+=t2-t1;
      printf("Time used for LanczosTriDiagonalize %i\n",(int)(t2-t1));
      time(&t1);
      TriDiagonalMatrixToSpectrum(Energy,s,&SpecSPINUp[0]);
      time(&t2);
      tTriDiagonalMatrixToSpectrum+=t2-t1;
      printf("Time used for TriDiagonalMatrixToSpectrum %i\n",(int)(t2-t1));
      fflush(stdout);
      T.NCrAn=2;
      T.VCrAnRe[0]=sqrt(0.5);
      T.VCrAnRe[1]=sqrt(0.5);
      for(j=1;j<psi.NFermion/2;++j)
      {
        T.iCrAn[1]=j;
        T.jCrAn[1]=j+H.NFermion/2;
        OperatorPsi(T,psi,&Opsi);
        LanczosTriDiagonalize(H, Opsi, HSet.zeroLancWavefunction,&s,HSet.Gamma);
        TriDiagonalMatrixToSpectrum(Energy,s,&SpecSPINUp[j]);
      }
      FreeOperator(&T);
      
        
      printf("\n\nSPIN Dn\n\n");
      T.NCrAn=1;
      T.iCrAn=malloc(2*sizeof(uint32_t));
      T.jCrAn=malloc(2*sizeof(uint32_t));
      T.VCrAnRe=malloc(2*sizeof(double));
      T.jCrAn[0]=0;
      T.iCrAn[0]=H.NFermion/2;
      T.VCrAnRe[0]=1.0;
      OperatorPsi(T,psi,&Opsi);
      time(&t1);
      LanczosTriDiagonalize(H, Opsi, HSet.zeroLancWavefunction,&s,HSet.Gamma);
      time(&t2);
      tLanczosTriDiagonalize+=t2-t1;
      printf("Time used for LanczosTriDiagonalize %i\n",(int)(t2-t1));
      time(&t1);
      TriDiagonalMatrixToSpectrum(Energy,s,&SpecSPINDn[0]);
      time(&t2);
      tTriDiagonalMatrixToSpectrum+=t2-t1;
      printf("Time used for TriDiagonalMatrixToSpectrum %i\n",(int)(t2-t1));
      fflush(stdout);
      T.NCrAn=2;
      T.VCrAnRe[0]=sqrt(0.5);
      T.VCrAnRe[1]=sqrt(0.5);
      for(j=1;j<psi.NFermion/2;++j)
      {
        T.jCrAn[1]=j;
        T.iCrAn[1]=j+H.NFermion/2;
        OperatorPsi(T,psi,&Opsi);
        LanczosTriDiagonalize(H, Opsi, HSet.zeroLancWavefunction,&s,HSet.Gamma);
        TriDiagonalMatrixToSpectrum(Energy,s,&SpecSPINDn[j]);
      }
      FreeOperator(&T);

      
      printf("\n Fourier transform spectra within cluster\n");
      SiteFourierTransform(HSet.NSite, HSet.NE, HSet.AtomPos, HSet.rUnitCell, HSet.kVector, SpecSPINUp);
      SiteFourierTransform(HSet.NSite, HSet.NE, HSet.AtomPos, HSet.rUnitCell, HSet.kVector, SpecSPINDn);
      SiteFourierTransform(HSet.NSite, HSet.NE, HSet.AtomPos, HSet.rUnitCell, HSet.kVector, SpecPES);
      SiteFourierTransform(HSet.NSite, HSet.NE, HSet.AtomPos, HSet.rUnitCell, HSet.kVector, SpecIPES);
      
      
      
      printf("\n\nSpectra written to file %s\n",flname);

      stream=fopen(flname,"w");
      if (stream==NULL)
      {
        printf("Can't open file %s\n",flname);
        printf("That's bad program terminates\n");
        exit(1);
      }  
      for(i=0;i<=HSet.NE;++i)
      {
        fprintf(stream,"%22.15E        %14.7E %14.7E ",HSet.Emin+(HSet.Emax-HSet.Emin)*i/((double)HSet.NE),SpecOpt.Re[i],SpecOpt.Im[i]);
        for(j=0;j<psi.NFermion/2;++j)
        {
          fprintf(stream,"%14.7E %14.7E ",SpecPES[j].Re[i],SpecPES[j].Im[i]);
        }
        for(j=0;j<psi.NFermion/2;++j)
        {
          fprintf(stream,"%14.7E %14.7E ",SpecIPES[j].Re[i],SpecIPES[j].Im[i]);
        }
        for(j=0;j<psi.NFermion/2;++j)
        {
          fprintf(stream,"%14.7E %14.7E ",SpecIPES[j].Re[i],SpecSPINUp[j].Im[i]);
        }
        for(j=0;j<psi.NFermion/2;++j)
        {
          fprintf(stream,"%14.7E %14.7E ",SpecIPES[j].Re[i],SpecSPINDn[j].Im[i]);
        }
        fprintf(stream,"\n");
      }
      fprintf(stream,"\n");
      fclose(stream);
    }// End of loop over different intial states //
    printf("\nEnergy overview\n");
    for(i=0;i<HSet.NPsi;i++)
    {
      printf("EigenState %6i  Energy %22.15E  N(BasisSet) %12i\n",i,EnergyList[i],PsiList[i].N);
    }
    printf("\nOrthogonality test\n");
    for(i=0;i<HSet.NPsi;i++)
    {
      for(j=0;j<HSet.NPsi;j++)
      {
        printf("%10.3E ",RealWaveFunctionDotProduct(PsiList[i],PsiList[j]));
      }
      printf("\n");
    }
    
    fflush(stdout);
    
  }// End of loop over different accuracies //
  
  time(&tTotal);
  printf("\nTiming Results\n");
  printf("Total time                                : %20i seconds\n",(int)(tTotal-t0));
  printf("Time spent in LanczosGroundState          : %20i seconds\n",(int)(tLanczosGroundState));
  printf("Time spent in LanczosTriDiagonalize       : %20i seconds\n",(int)(tLanczosTriDiagonalize));
  printf("Time spent in TriDiagonalMatrixToSpectrum : %20i seconds\n",(int)(tTriDiagonalMatrixToSpectrum));
};
*/


/*
void mainTwistedHubbardModel (HubbardModelSettings HSet) 
{
  WaveFunctionType psi,Opsi,*PsiList;
  OperatorType H,T;
  TriDiagonalMatrixType s;
  double TotEmax,Energy,*EnergyList,junk,dotRe,dotIm;
  SpectraType *SpecPES,*SpecIPES,SpecOpt,*SpecSPINUp,*SpecSPINDn;
  uint32_t i,j,jacc,jpsi,iq[3];
  char flname[64];
  FILE *stream;
  time_t t0,t1,t2,tTotal=0,tLanczosGroundState=0,tLanczosTriDiagonalize=0,tTriDiagonalMatrixToSpectrum=0;
  
  time(&t0);

  s.N=HSet.NTriDiag;
  H.NFermion=(-(HSet.rUnitCell[0][2]*HSet.rUnitCell[1][1]) + HSet.rUnitCell[0][1]*HSet.rUnitCell[1][2])*HSet.rUnitCell[2][0] + (HSet.rUnitCell[0][2]*HSet.rUnitCell[1][0] - HSet.rUnitCell[0][0]*HSet.rUnitCell[1][2])*HSet.rUnitCell[2][1] + (-(HSet.rUnitCell[0][1]*HSet.rUnitCell[1][0]) + HSet.rUnitCell[0][0]*HSet.rUnitCell[1][1])*HSet.rUnitCell[2][2];
  printf("Number of sites in Hubbard Model %i\n",H.NFermion);
  H.NFermion*=2;
  // allocate memory //
  
  SpecOpt.N=HSet.NE;
  SpecOpt.Emin=HSet.Emin;
  SpecOpt.Emax=HSet.Emax;
  SpecOpt.dE=(HSet.Emax-HSet.Emin)/((double)HSet.NE);
  SpecOpt.EminPole=HSet.Emin;
  SpecOpt.EmaxPole=HSet.Emax;
  SpecOpt.Gamma=HSet.Gamma;

  SpectrumAllocate(SpecOpt,&SpecOpt);
  
  SpecPES=malloc(H.NFermion/2*sizeof(SpectraType));
  SpecIPES=malloc(H.NFermion/2*sizeof(SpectraType));
  SpecSPINUp=malloc(H.NFermion/2*sizeof(SpectraType));
  SpecSPINDn=malloc(H.NFermion/2*sizeof(SpectraType));
  for(j=0;j<H.NFermion;++j)
  {
    SpectrumAllocate(SpecOpt,&SpecPES[j]);
    SpectrumAllocate(SpecOpt,&SpecIPES[j]);
    SpectrumAllocate(SpecOpt,&SpecSPINUp[j]);
    SpectrumAllocate(SpecOpt,&SpecSPINDn[j]);
  }
  
  s.a=malloc(s.N*sizeof(double));
  s.b=malloc(s.N*sizeof(double));
  PsiList=malloc(HSet.NPsi*sizeof(WaveFunctionType));
  EnergyList=malloc(HSet.NPsi*sizeof(double));

  for(jpsi=0;jpsi<HSet.NPsi;jpsi++)
  {
    PsiList[jpsi]=CreateStartComplexWaveFunitionNeelState(HSet);
  }
  // loop over all pseudo k points the twisted boundary conditions //
  for(iq[0]=0;iq[0]<HSet.Nq[0];iq[0]++)
  { for(iq[1]=0;iq[1]<HSet.Nq[1];iq[1]++)
    { for(iq[2]=0;iq[2]<HSet.Nq[2];iq[2]++)
      {
        printf("TwistedHubbard Model with only nearest neighbor hopping:\n");
        H=CreateSimpelTwistedHubbardModel(HSet,iq);
        fflush(stdout);
        TotEmax=Min(HSet.NelecUp,HSet.NelecDn)*HSet.U+(HSet.NelecUp+HSet.NelecDn)*2*fabs(HSet.t)-0.5*(HSet.NelecUp+HSet.NelecDn)*HSet.U;
        // loop over different accuracies //
        for(jacc=0;jacc<=HSet.Ni;jacc++)
        {
          if(jacc!=0)
          {
            HSet.zeroLancWavefunction/=10;
            HSet.zeroLancEnergy/=10;
            if((HSet.zeroLancEnergy/(TotEmax-EnergyList[0]))<1.0E-14)HSet.zeroLancEnergy*=10;
            printf("\n=================\nIncreased accuracy to %22.15E %22.15E \n=================\n",HSet.zeroLancWavefunction,HSet.zeroLancEnergy);
            fflush(stdout);
          }
          // loop over number of eigen-states //
          for(jpsi=0;jpsi<HSet.NPsi;jpsi++)
          {
            sprintf(flname,"Opt_PES_IPES_SPIN_%i_%i_q%iq%iq%i.dat",jacc,jpsi,iq[0],iq[1],iq[2]);
            psi=PsiList[jpsi];
            printf("Starting density\n");
            ComplexWaveFunctionRandomize(&psi);
            PrintWaveFunctionDensity(psi);
            printf("Number of determinants used for starting function %i\n",psi.N);
            printf("Get ground-state\n");
            fflush(stdout);
            time(&t1);
            ComplexLanczosGroundState(H,&psi,&HSet.zeroLancWavefunction,HSet.zeroLancEnergy,TotEmax,PsiList,jpsi);
            time(&t2);
            tLanczosGroundState+=t2-t1;
            PsiList[jpsi]=psi;
            PrintWaveFunctionDensity(psi);
            printf("Number of determinants used for eigen function %i\n",psi.N);
            printf("Time used for LanczosGroundState %i\n",(int)(t2-t1));
            fflush(stdout);
            if(ComplexOperatorPsiConserveBasis(H,psi,&Opsi))AllocFailed();
            ComplexWaveFunctionDotProduct(Opsi,psi,&Energy,&junk);
            EnergyList[jpsi]=Energy;
            printf("Ground state energy is %22.15E + I (%22.15E)\n",Energy,junk);
            ComplexWaveFunctionFreeValues(&Opsi);
            fflush(stdout);
  
            printf("Optical conductivety \n\n");
            CreateSimpelTwistedHubbardModelOpticsT(HSet,&T,iq);
            ComplexOperatorPsi(T,psi,&Opsi);
            time(&t1);
            ComplexLanczosTriDiagonalize(H, Opsi, HSet.zeroLancWavefunction,&s,HSet.Gamma);
            time(&t2);
            tLanczosTriDiagonalize+=t2-t1;
            printf("Time used for LanczosTriDiagonalize %i\n",(int)(t2-t1));
            time(&t1);
            TriDiagonalMatrixToSpectrum(Energy,s,&SpecOpt);
            time(&t2);
            tTriDiagonalMatrixToSpectrum+=t2-t1;
            printf("Time used for TriDiagonalMatrixToSpectrum %i\n",(int)(t2-t1));
            fflush(stdout);
            FreeOperator(&T);
      
            printf("\n\nPES\n\n");
            T.NAn=HSet.NSite;
            T.iAn=malloc(T.NAn*sizeof(uint32_t));
            T.VAnRe=malloc(T.NAn*sizeof(double));
            T.VAnIm=malloc(T.NAn*sizeof(double));
            for(j=0;j<psi.NFermion/2;j++)
            {
              for(i=0;i<psi.NFermion/2;i++)
              {
                T.iAn[i]=i;
                T.VAnRe[i]=(1.0/sqrt((double)HSet.NSite))*cos((double)Dot(HSet.AtomPos[i],HSet.kVector[j]));
                T.VAnIm[i]=-(1.0/sqrt((double)HSet.NSite))*sin((double)Dot(HSet.AtomPos[i],HSet.kVector[j]));
              }
              ComplexOperatorPsi(T,psi,&Opsi);
              ComplexLanczosTriDiagonalize(H, Opsi, HSet.zeroLancWavefunction,&s,HSet.Gamma);
              TriDiagonalMatrixToSpectrum(Energy,s,&SpecPES[j]);
            }
            FreeOperator(&T);

            printf("\n\nIPES\n\n");
            T.NCr=HSet.NSite;
            T.iCr=malloc(T.NCr*sizeof(uint32_t));
            T.VCrRe=malloc(T.NCr*sizeof(double));
            T.VCrIm=malloc(T.NCr*sizeof(double));
            for(j=0;j<psi.NFermion/2;j++)
            {
              for(i=0;i<psi.NFermion/2;i++)
              {
                T.iCr[i]=i;
                T.VCrRe[i]=(1.0/sqrt((double)HSet.NSite))*cos((double)Dot(HSet.AtomPos[i],HSet.kVector[j]));
                T.VCrIm[i]=(1.0/sqrt((double)HSet.NSite))*sin((double)Dot(HSet.AtomPos[i],HSet.kVector[j]));
              }
              ComplexOperatorPsi(T,psi,&Opsi);
              ComplexLanczosTriDiagonalize(H, Opsi, HSet.zeroLancWavefunction,&s,HSet.Gamma);
              TriDiagonalMatrixToSpectrum(Energy,s,&SpecIPES[j]);
            }
            FreeOperator(&T);


            printf("\n\nSPIN Up\n\n");
            T.NCrAn=HSet.NSite;
            T.iCrAn=malloc(T.NCrAn*sizeof(uint32_t));
            T.jCrAn=malloc(T.NCrAn*sizeof(uint32_t));
            T.VCrAnRe=malloc(T.NCrAn*sizeof(double));
            T.VCrAnIm=malloc(T.NCrAn*sizeof(double));
            for(j=0;j<psi.NFermion/2;j++)
            {
              for(i=0;i<psi.NFermion/2;i++)
              {
                T.iCrAn[i]=0;
                T.jCrAn[i]=i+psi.NFermion/2;
                T.VCrAnRe[i]=(1.0/sqrt((double)HSet.NSite))*cos((double)(Dot(HSet.AtomPos[i],HSet.kVector[j])-Dot(HSet.AtomPos[0],HSet.kVector[j])));
                T.VCrAnIm[i]=-(1.0/sqrt((double)HSet.NSite))*sin((double)(Dot(HSet.AtomPos[i],HSet.kVector[j])-Dot(HSet.AtomPos[0],HSet.kVector[j])));
              }
              ComplexOperatorPsi(T,psi,&Opsi);
              ComplexLanczosTriDiagonalize(H, Opsi, HSet.zeroLancWavefunction,&s,HSet.Gamma);
              TriDiagonalMatrixToSpectrum(Energy,s,&SpecSPINUp[j]);
            }
            FreeOperator(&T); 
            
                         
            printf("\n\nSPIN Dn\n\n");
            T.NCrAn=2*HSet.NSite;
            T.iCrAn=malloc(T.NCrAn*sizeof(uint32_t));
            T.jCrAn=malloc(T.NCrAn*sizeof(uint32_t));
            T.VCrAnRe=malloc(T.NCrAn*sizeof(double));
            T.VCrAnIm=malloc(T.NCrAn*sizeof(double));
            for(j=0;j<psi.NFermion/2;j++)
            {
              for(i=0;i<psi.NFermion/2;i++)
              {
                T.iCrAn[i]=i;
                T.jCrAn[i]=i;
                T.iCrAn[i+psi.NFermion/2]=i+psi.NFermion/2;
                T.jCrAn[i+psi.NFermion/2]=i+psi.NFermion/2;
                T.VCrAnRe[i]=(1.0/sqrt((double)HSet.NSite))*cos((double)Dot(HSet.AtomPos[i],HSet.kVector[j]));
                T.VCrAnIm[i]=(1.0/sqrt((double)HSet.NSite))*sin((double)Dot(HSet.AtomPos[i],HSet.kVector[j]));
                T.VCrAnRe[i+psi.NFermion/2]=-(1.0/sqrt((double)HSet.NSite))*cos((double)Dot(HSet.AtomPos[i],HSet.kVector[j]));
                T.VCrAnIm[i+psi.NFermion/2]=-(1.0/sqrt((double)HSet.NSite))*sin((double)Dot(HSet.AtomPos[i],HSet.kVector[j]));
              }
              ComplexOperatorPsi(T,psi,&Opsi);
              ComplexLanczosTriDiagonalize(H, Opsi, HSet.zeroLancWavefunction,&s,HSet.Gamma);
              TriDiagonalMatrixToSpectrum(Energy,s,&SpecSPINDn[j]);
            }
            FreeOperator(&T);              
  
            printf("\n\nSpectra written to file %s\n",flname);

            stream=fopen(flname,"w");
            if (stream==NULL)
            {
              printf("Can't open file %s\n",flname);
              printf("That's bad program terminates\n");
              exit(1);
            }

            for(i=0;i<=HSet.NE;++i)
            {
              fprintf(stream,"%22.15E        %14.7E %14.7E ",HSet.Emin+(HSet.Emax-HSet.Emin)*i/((double)HSet.NE),SpecOpt.Re[i],SpecOpt.Im[i]);
              for(j=0;j<psi.NFermion/2;++j)
              {
                fprintf(stream,"%14.7E %14.7E ",SpecPES[j].Re[i],SpecPES[j].Im[i]);
              }
              for(j=0;j<psi.NFermion/2;++j)
              {
                fprintf(stream,"%14.7E %14.7E ",SpecIPES[j].Re[i],SpecIPES[j].Im[i]);
              }
              for(j=0;j<psi.NFermion/2;++j)
              {
                fprintf(stream,"%14.7E %14.7E ",SpecSPINUp[j].Re[i],SpecSPINUp[j].Im[i]);
              }
              for(j=0;j<psi.NFermion/2;++j)
              {
                fprintf(stream,"%14.7E %14.7E ",SpecSPINDn[j].Re[i],SpecSPINDn[j].Im[i]);
              }
              fprintf(stream,"\n");
            }
            fprintf(stream,"\n");
            fclose(stream);
          }// End of loop over different intial states //
          printf("\nEnergy overview\n");
          for(i=0;i<HSet.NPsi;i++)
          {
            printf("EigenState %6i  Energy %22.15E  N(BasisSet) %12i\n",i,EnergyList[i],PsiList[i].N);
          }
          printf("\nOrthogonality test\n");
          for(i=0;i<HSet.NPsi;i++)
          {
            for(j=0;j<HSet.NPsi;j++)
            {
              ComplexWaveFunctionDotProduct(PsiList[i],PsiList[j],&dotRe,&dotIm);
              printf("%10.3E %10.3E",dotRe,dotIm);
            }
            printf("\n");
          }
    
          fflush(stdout);
    
        }// End of loop over different accuracies //
        FreeOperator(&H);
      } // End of loop over iq[2] //
    } // End of loop over iq[1] //
  } // End of loop over iq[0] //

  
  time(&tTotal);
  printf("\nTiming Results\n");
  printf("Total time                                : %20i seconds\n",(int)(tTotal-t0));
  printf("Time spent in LanczosGroundState          : %20i seconds\n",(int)(tLanczosGroundState));
  printf("Time spent in LanczosTriDiagonalize       : %20i seconds\n",(int)(tLanczosTriDiagonalize));
  printf("Time spent in TriDiagonalMatrixToSpectrum : %20i seconds\n",(int)(tTriDiagonalMatrixToSpectrum));
};
*/


uint32_t HubbardModelCreateStartWaveFunction(HubbardModelSettings Set, WaveFunctionType * psi)
{
  uint8_t *Determinant;
  uint32_t iatom, ishell, ifermion, icell, ielectron, i;
  psi->QComplex = 0;
  psi->NBitsKey = Set.NBitsKey;
  psi->NFermion = Set.USuperCell.NFermion;
  psi->NBoson = Set.USuperCell.NBoson;
  psi->NWord = Set.USuperCell.NWord;
  psi->N = 0;
  psi->NBlok = 1;
  psi->KeyList = calloc((uint32_t) 1 << (psi->NBitsKey), sizeof(uint32_t));
  psi->ValueRe = malloc(psi->NBlok * sizeof(double *));
  psi->Next = malloc(psi->NBlok * sizeof(uint32_t *));
  psi->Determinant = malloc(psi->NBlok * sizeof(uint8_t *));
  if ((psi->ValueRe == NULL) || (psi->Next == NULL) || (psi->Determinant == NULL))
    AllocFailed();
  psi->ValueRe[psi->NBlok - 1] = malloc(BlokSize * sizeof(double));
  psi->Next[psi->NBlok - 1] = malloc(BlokSize * sizeof(uint32_t));
  psi->Determinant[psi->NBlok - 1] = malloc(BlokSize * psi->NWord * sizeof(uint8_t));
  if ((psi->ValueRe[psi->NBlok - 1] == NULL) || (psi->Next[psi->NBlok - 1] == NULL)
      || (psi->Determinant[psi->NBlok - 1] == NULL))
    AllocFailed();
  Determinant = calloc(psi->NWord, sizeof(uint8_t));
  for (iatom = 0; iatom < Set.SuperCell.Crystal.NAtoms; iatom++)
  {
    for (ishell = 0; ishell < Set.SuperCell.Crystal.Atom[iatom].Nshell; ishell++)
    {
      for (ifermion = 0; ifermion < Set.SuperCell.Crystal.Atom[iatom].shelldegeneracy[ishell]; ifermion++)
      {
        for (ielectron = 0; ielectron < Set.NElectron[iatom][ishell][ifermion]; ielectron++)
        {
          icell = ielectron;    /* needs to be changed */
          i = ifermion + Set.TightBinding.FermionIndex[iatom][ishell] + icell * Set.TightBinding.NFermion;
          OperatorCr(Determinant, i, psi->NFermion);
        }
      }
    }
  }
  RealWaveFunctionAddElement(psi, Determinant, 1.0);
  sprintf(psi->Name, "Start Wavefunction");
  return 0;
}

uint32_t HubbardModelEnlargeStartWaveFunction(OperatorType H, WaveFunctionType * psi)
{
  WaveFunctionType OPsi;
  OperatorPsi(H, *psi, &OPsi);
  RealWaveFunctionFree(psi);
  RealWaveFunctionRandomize(&OPsi);
  RealWaveFunctionNormalize(&OPsi);
  OperatorPsi(H, OPsi, psi);
  RealWaveFunctionRandomize(psi);
  RealWaveFunctionNormalize(psi);
  return 0;
}


int main(int argc, const char *argv[])
{
  uint32_t iq[3], jpsi, i, j, Npsi;
  double q[3], Hist[20], TotEmax, *EnergyList, *WaveFunctionQuality, tempRe, tempIm;
  OperatorType H;
  WaveFunctionType *PsiList;
  MatrixType Rho, RhoTot, RhoCorr, RhoCorrTot;
  TimeType Time;
  InitializeTiming(&Time);
  BeginTiming("Total", &Time);
  TimeStampAll();
  argc -= 1;
  HubbardModelSettings HSet;
  FileReadHubbardModelSettings(&HSet, argv[argc]);
  PrintHubbardModelSettings(HSet);
  PsiList = malloc(HSet.NPsi * sizeof(WaveFunctionType));
  EnergyList = malloc(HSet.NPsi * sizeof(double));
  WaveFunctionQuality = malloc(HSet.NPsi * sizeof(double));
  for (iq[0] = 0; iq[0] < HSet.Nq[0]; iq[0]++)
  {
    for (iq[1] = 0; iq[1] < HSet.Nq[1]; iq[1]++)
    {
      for (iq[2] = 0; iq[2] < HSet.Nq[2]; iq[2]++)
      {
        BeginTiming("Total_new_q", &Time);
        q[0] =
          ((double)iq[0]) / ((double)HSet.Nq[0]) * HSet.SuperCell.gSuperCell[0][0] +
          ((double)iq[1]) / ((double)HSet.Nq[1]) * HSet.SuperCell.gSuperCell[1][0] +
          ((double)iq[2]) / ((double)HSet.Nq[2]) * HSet.SuperCell.gSuperCell[2][0];
        q[1] =
          ((double)iq[0]) / ((double)HSet.Nq[0]) * HSet.SuperCell.gSuperCell[0][1] +
          ((double)iq[1]) / ((double)HSet.Nq[1]) * HSet.SuperCell.gSuperCell[1][1] +
          ((double)iq[2]) / ((double)HSet.Nq[2]) * HSet.SuperCell.gSuperCell[2][1];
        q[2] =
          ((double)iq[0]) / ((double)HSet.Nq[0]) * HSet.SuperCell.gSuperCell[0][2] +
          ((double)iq[1]) / ((double)HSet.Nq[1]) * HSet.SuperCell.gSuperCell[1][2] +
          ((double)iq[2]) / ((double)HSet.Nq[2]) * HSet.SuperCell.gSuperCell[2][2];
        printf("\n\nStarted with q vector %i %i %i is %22.15lf %22.15lf %22.15lf\n", iq[0], iq[1], iq[2], q[0], q[1],
               q[2]);
        CreateHamiltonianOperatorFromTightBindingDef(HSet.TightBindingSuperCell, q, &H);
        OperatorAdd(&H, HSet.USuperCell, 1, 0);
        PrintOperatorReadable(H);
        if (QOperatorHermitian(H, DBL_EPSILON))
        {
          printf("Operator is Hermitian\n");
        }
        else
        {
          printf("Operator is NOT Hermitian\n");
        }
        Npsi = 0;
        for (jpsi = 0; jpsi < HSet.NPsi; jpsi++)
        {
          Npsi++;
          HubbardModelCreateStartWaveFunction(HSet, &PsiList[jpsi]);
          HubbardModelEnlargeStartWaveFunction(H, &PsiList[jpsi]);
          if (PsiList[jpsi].QComplex == 0)
          {
            if (jpsi != 0)
              RealWaveFunctionOrthogonalize(&PsiList[jpsi], PsiList, jpsi);
            RealWaveFunctionRandomize(&PsiList[jpsi]);
          }
          else
          {
            if (jpsi != 0)
              ComplexWaveFunctionOrthogonalize(&PsiList[jpsi], PsiList, jpsi);
            ComplexWaveFunctionRandomize(&PsiList[jpsi]);
          }
          TotEmax = 10.0;       /* Have to change this */
          BeginTiming("LanczosGroundState", &Time);
          LanczosGroundState(H, &PsiList[jpsi], &HSet.zeroLancWavefunction, HSet.zeroLancEnergy, TotEmax, PsiList,
                             jpsi);
          EndTiming("LanczosGroundState", &Time);
          TotEmax *= 2;
          BeginTiming("LanczosGroundState", &Time);
          LanczosGroundState(H, &PsiList[jpsi], &HSet.zeroLancWavefunction, HSet.zeroLancEnergy, TotEmax, PsiList,
                             jpsi);
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

              WaveFunctionDotProduct(PsiList[i], PsiList[j], &tempRe, &tempIm);
              if (fabs(tempRe) < DBL_EPSILON)
                tempRe = 0;
              printf(" %22.15E", tempRe);
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
          if ((EnergyList[jpsi] - EnergyList[0]) > HSet.DeltaEMax)
          {
            printf("Average of Density matrices for state 0 to %i\n", jpsi - 1);
            MatrixScalarMultiply(&RhoTot, 1.0 / ((double)(jpsi)), 0);
            MatrixScalarMultiply(&RhoCorrTot, 1.0 / ((double)(jpsi)), 0);
            PrintMatrix(RhoTot);
            PrintMatrix(RhoCorrTot);
            goto SkipExcitedStateCalculation;
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
//          for(i=0;i<Set.NSite*4;i++)
//          {
//            BeginTiming("Lanczos Tri Diagonalization",&Time);
//            OperatorPsi(T[i],PsiList[jpsi],&OPsi);
//            LanczosTriDiagonalize(H, OPsi, Set.zeroLancWavefunction,&s,Set.Gamma);
//            TriDiagonalMatrixToSpectrum(EnergyList[jpsi],s,&Spec[i]);
//            EndTiming("Lanczos Tri Diagonalization",&Time);
//          }
//          if(jpsi<10)        sprintf(SpecFileName,"%s00%1i",Set.SpecFile,jpsi);
//          else if(jpsi<100)  sprintf(SpecFileName,"%s0%2i" ,Set.SpecFile,jpsi);
//          else if(jpsi<1000) sprintf(SpecFileName,"%s%3i"  ,Set.SpecFile,jpsi);
//          else               sprintf(SpecFileName,"%s%i"   ,Set.SpecFile,jpsi);
//          FilePrintSpectra(Spec, 4*Set.NSite,SpecFileName,"w");

        }

      SkipExcitedStateCalculation:
        FreeOperator(&H);
        for (jpsi = 0; jpsi < Npsi; jpsi++)
        {
          RealWaveFunctionFree(&PsiList[jpsi]);
        }
        FreeMatrix(&RhoTot);
        FreeMatrix(&RhoCorrTot);

        EndTiming("Total_new_q", &Time);
        PrintTiming(Time);
      }
    }
  }


  EndTiming("Total", &Time);
  PrintTiming(Time);
  return 0;
};
