#include "tjModel.h"
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <omp.h>
#include <limits.h>
#include <unistd.h>
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



/*  Creates a tJ-model Hamiltonian
 *  the fermions at position 2 i     correspond to the real-space position i with spin down
 *  the fermions at position 2 i + 1 correspont to the real-space position i wiht spin up
 *
 *  Si . Sj = Szi Szj + 1/2 ( Si+ Sj- + Si- Sj+ )
 *  Sz = 1/2 (a+_{up} a_{up} - a+_{dn} a_{dn})
 *  S+ = a+_{up} a_{dn}
 *  S- = a+_{dn} a_{up}
 */
uint32_t CreatetJModelHamiltonian(tJModelSettings Set, OperatorType * H)
{
  uint32_t i, j;
  int32_t C[3], A[3], B[3];
  uint16_t oadd[8];
  InitOperator(H);
  H->NFermion = 2 * Set.Cluster.NCell;
  H->NWord = Min(H->NFermion, 1 + (H->NFermion - 1) / 8);
  InitOperatorFixedLength(H, 2, 8);
  InitOperatorFixedLength(H, 4, 8);
  InitOperatorFixedLength(H, 6, 8);
  for (i = 0; i < Set.Cluster.NCell; i++)
  {
    A[0] = Set.Cluster.CellPos[i][0];
    A[1] = Set.Cluster.CellPos[i][1];
    A[2] = Set.Cluster.CellPos[i][2];
    for (j = 0; j < Set.Cluster.NCell; j++)
    {
      B[0] = Set.Cluster.CellPos[j][0];
      B[1] = Set.Cluster.CellPos[j][1];
      B[2] = Set.Cluster.CellPos[j][2];
      ShortesConnectingVectorPeriodicSuperCell(A, B, Set.Cluster, C);
      if (Dot(C, C) == 1)
      {                         /* add Si . Sj to the Hamiltonian */
        oadd[3] = CrAnMask + 2 * i + 1;
        oadd[2] = 2 * i + 1;
        oadd[1] = CrAnMask + 2 * j + 1;
        oadd[0] = 2 * j + 1;
        RealOperatorAddLaderNormalOrder(H, oadd, 4, 0.25 * Set.J, &QDetNotOrderedFcFaBcBa);
        oadd[3] = CrAnMask + 2 * i;
        oadd[2] = 2 * i;
        oadd[1] = CrAnMask + 2 * j;
        oadd[0] = 2 * j;
        RealOperatorAddLaderNormalOrder(H, oadd, 4, 0.25 * Set.J, &QDetNotOrderedFcFaBcBa);
        oadd[3] = CrAnMask + 2 * i;
        oadd[2] = 2 * i;
        oadd[1] = CrAnMask + 2 * j + 1;
        oadd[0] = 2 * j + 1;
        RealOperatorAddLaderNormalOrder(H, oadd, 4, -0.25 * Set.J, &QDetNotOrderedFcFaBcBa);
        oadd[3] = CrAnMask + 2 * i + 1;
        oadd[2] = 2 * i + 1;
        oadd[1] = CrAnMask + 2 * j;
        oadd[0] = 2 * j;
        RealOperatorAddLaderNormalOrder(H, oadd, 4, -0.25 * Set.J, &QDetNotOrderedFcFaBcBa);
        oadd[3] = CrAnMask + 2 * i + 1;
        oadd[2] = 2 * i;
        oadd[1] = CrAnMask + 2 * j;
        oadd[0] = 2 * j + 1;
        RealOperatorAddLaderNormalOrder(H, oadd, 4, 0.5 * Set.J, &QDetNotOrderedFcFaBcBa);
        oadd[3] = CrAnMask + 2 * i;
        oadd[2] = 2 * i + 1;
        oadd[1] = CrAnMask + 2 * j + 1;
        oadd[0] = 2 * j;
        RealOperatorAddLaderNormalOrder(H, oadd, 4, 0.5 * Set.J, &QDetNotOrderedFcFaBcBa);
        /* subtract density density term */
        oadd[3] = CrAnMask + 2 * i + 1;
        oadd[2] = 2 * i + 1;
        oadd[1] = CrAnMask + 2 * j + 1;
        oadd[0] = 2 * j + 1;
        RealOperatorAddLaderNormalOrder(H, oadd, 4, -0.25 * Set.J, &QDetNotOrderedFcFaBcBa);
        oadd[3] = CrAnMask + 2 * i + 1;
        oadd[2] = 2 * i + 1;
        oadd[1] = CrAnMask + 2 * j;
        oadd[0] = 2 * j;
        RealOperatorAddLaderNormalOrder(H, oadd, 4, -0.25 * Set.J, &QDetNotOrderedFcFaBcBa);
        oadd[3] = CrAnMask + 2 * i;
        oadd[2] = 2 * i;
        oadd[1] = CrAnMask + 2 * j + 1;
        oadd[0] = 2 * j + 1;
        RealOperatorAddLaderNormalOrder(H, oadd, 4, -0.25 * Set.J, &QDetNotOrderedFcFaBcBa);
        oadd[3] = CrAnMask + 2 * i;
        oadd[2] = 2 * i;
        oadd[1] = CrAnMask + 2 * j;
        oadd[0] = 2 * j;
        RealOperatorAddLaderNormalOrder(H, oadd, 4, -0.25 * Set.J, &QDetNotOrderedFcFaBcBa);
        /* add hopping keeping double occupation constant */
        /* i to j spin dn */
        oadd[5] = CrAnMask + 2 * i;
        oadd[4] = CrAnMask + 2 * i + 1;
        oadd[3] = 2 * i + 1;
        oadd[2] = CrAnMask + 2 * j + 1;
        oadd[1] = 2 * j + 1;
        oadd[0] = 2 * j;
        RealOperatorAddLaderNormalOrder(H, oadd, 6, Set.t, &QDetNotOrderedFcFaBcBa);
        oadd[5] = CrAnMask + 2 * i;
        oadd[4] = 2 * i + 1;
        oadd[3] = CrAnMask + 2 * i + 1;
        oadd[2] = 2 * j + 1;
        oadd[1] = CrAnMask + 2 * j + 1;
        oadd[0] = 2 * j;
        RealOperatorAddLaderNormalOrder(H, oadd, 6, Set.t, &QDetNotOrderedFcFaBcBa);
        /* i to j spin up */
        oadd[5] = CrAnMask + 2 * i + 1;
        oadd[4] = CrAnMask + 2 * i;
        oadd[3] = 2 * i;
        oadd[2] = CrAnMask + 2 * j;
        oadd[1] = 2 * j;
        oadd[0] = 2 * j + 1;
        RealOperatorAddLaderNormalOrder(H, oadd, 6, Set.t, &QDetNotOrderedFcFaBcBa);
        oadd[5] = CrAnMask + 2 * i + 1;
        oadd[4] = 2 * i;
        oadd[3] = CrAnMask + 2 * i;
        oadd[2] = 2 * j;
        oadd[1] = CrAnMask + 2 * j;
        oadd[0] = 2 * j + 1;
        RealOperatorAddLaderNormalOrder(H, oadd, 6, Set.t, &QDetNotOrderedFcFaBcBa);
        /* j to i spin dn */
        oadd[5] = CrAnMask + 2 * j;
        oadd[4] = CrAnMask + 2 * j + 1;
        oadd[3] = 2 * j + 1;
        oadd[2] = CrAnMask + 2 * i + 1;
        oadd[1] = 2 * i + 1;
        oadd[0] = 2 * i;
        RealOperatorAddLaderNormalOrder(H, oadd, 6, Set.t, &QDetNotOrderedFcFaBcBa);
        oadd[5] = CrAnMask + 2 * j;
        oadd[4] = 2 * j + 1;
        oadd[3] = CrAnMask + 2 * j + 1;
        oadd[2] = 2 * i + 1;
        oadd[1] = CrAnMask + 2 * i + 1;
        oadd[0] = 2 * i;
        RealOperatorAddLaderNormalOrder(H, oadd, 6, Set.t, &QDetNotOrderedFcFaBcBa);
        /* j to i spin up */
        oadd[5] = CrAnMask + 2 * j + 1;
        oadd[4] = CrAnMask + 2 * j;
        oadd[3] = 2 * j;
        oadd[2] = CrAnMask + 2 * i;
        oadd[1] = 2 * i;
        oadd[0] = 2 * i + 1;
        RealOperatorAddLaderNormalOrder(H, oadd, 6, Set.t, &QDetNotOrderedFcFaBcBa);
        oadd[5] = CrAnMask + 2 * j + 1;
        oadd[4] = 2 * j;
        oadd[3] = CrAnMask + 2 * j;
        oadd[2] = 2 * i;
        oadd[1] = CrAnMask + 2 * i;
        oadd[0] = 2 * i + 1;
        RealOperatorAddLaderNormalOrder(H, oadd, 6, Set.t, &QDetNotOrderedFcFaBcBa);
      }
    }
  }
  OperatorRemoveSmallValues(H, DBL_MIN * 100);
  return 0;
}

uint32_t CreatetJModelHamiltonianLessThenHalfFilled(tJModelSettings Set, OperatorType * H)
{
  uint32_t i, j;
  int32_t C[3];
  uint16_t oadd[6];
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
        oadd[3] = CrAnMask + 2 * i + 1;
        oadd[2] = 2 * i + 1;
        oadd[1] = CrAnMask + 2 * j + 1;
        oadd[0] = 2 * j + 1;
        RealOperatorAddLaderNormalOrder(H, oadd, 4, 0.25 * Set.J, &QDetNotOrderedFcFaBcBa);
        oadd[3] = CrAnMask + 2 * i;
        oadd[2] = 2 * i;
        oadd[1] = CrAnMask + 2 * j;
        oadd[0] = 2 * j;
        RealOperatorAddLaderNormalOrder(H, oadd, 4, 0.25 * Set.J, &QDetNotOrderedFcFaBcBa);
        oadd[3] = CrAnMask + 2 * i;
        oadd[2] = 2 * i;
        oadd[1] = CrAnMask + 2 * j + 1;
        oadd[0] = 2 * j + 1;
        RealOperatorAddLaderNormalOrder(H, oadd, 4, -0.25 * Set.J, &QDetNotOrderedFcFaBcBa);
        oadd[3] = CrAnMask + 2 * i + 1;
        oadd[2] = 2 * i + 1;
        oadd[1] = CrAnMask + 2 * j;
        oadd[0] = 2 * j;
        RealOperatorAddLaderNormalOrder(H, oadd, 4, -0.25 * Set.J, &QDetNotOrderedFcFaBcBa);
        oadd[3] = CrAnMask + 2 * i + 1;
        oadd[2] = 2 * i;
        oadd[1] = CrAnMask + 2 * j;
        oadd[0] = 2 * j + 1;
        RealOperatorAddLaderNormalOrder(H, oadd, 4, 0.5 * Set.J, &QDetNotOrderedFcFaBcBa);
        oadd[3] = CrAnMask + 2 * i;
        oadd[2] = 2 * i + 1;
        oadd[1] = CrAnMask + 2 * j + 1;
        oadd[0] = 2 * j;
        RealOperatorAddLaderNormalOrder(H, oadd, 4, 0.5 * Set.J, &QDetNotOrderedFcFaBcBa);
        /* subtract density density term */
        oadd[3] = CrAnMask + 2 * i + 1;
        oadd[2] = 2 * i + 1;
        oadd[1] = CrAnMask + 2 * j + 1;
        oadd[0] = 2 * j + 1;
        RealOperatorAddLaderNormalOrder(H, oadd, 4, -0.25 * Set.J, &QDetNotOrderedFcFaBcBa);
        oadd[3] = CrAnMask + 2 * i;
        oadd[2] = 2 * i;
        oadd[1] = CrAnMask + 2 * j + 1;
        oadd[0] = 2 * j + 1;
        RealOperatorAddLaderNormalOrder(H, oadd, 4, -0.25 * Set.J, &QDetNotOrderedFcFaBcBa);
        oadd[3] = CrAnMask + 2 * i + 1;
        oadd[2] = 2 * i + 1;
        oadd[1] = CrAnMask + 2 * j;
        oadd[0] = 2 * j;
        RealOperatorAddLaderNormalOrder(H, oadd, 4, -0.25 * Set.J, &QDetNotOrderedFcFaBcBa);
        oadd[3] = CrAnMask + 2 * i;
        oadd[2] = 2 * i;
        oadd[1] = CrAnMask + 2 * j;
        oadd[0] = 2 * j;
        RealOperatorAddLaderNormalOrder(H, oadd, 4, -0.25 * Set.J, &QDetNotOrderedFcFaBcBa);
        /* add hopping keeping double occupation constant */
        /* hop from j_up to i_up if i_dn is not occupied (1-st line rest permutation in i and j and spin */
        oadd[3] = CrAnMask + 2 * i + 1;
        oadd[2] = 2 * j + 1;
        oadd[1] = 2 * i;
        oadd[0] = CrAnMask + 2 * i;
        RealOperatorAddLader(H, oadd, 4, Set.t);
        oadd[3] = CrAnMask + 2 * i;
        oadd[2] = 2 * j;
        oadd[1] = 2 * i + 1;
        oadd[0] = CrAnMask + 2 * i + 1;
        RealOperatorAddLader(H, oadd, 4, Set.t);
        oadd[3] = CrAnMask + 2 * j + 1;
        oadd[2] = 2 * i + 1;
        oadd[1] = 2 * j;
        oadd[0] = CrAnMask + 2 * j;
        RealOperatorAddLader(H, oadd, 4, Set.t);
        oadd[3] = CrAnMask + 2 * j;
        oadd[2] = 2 * i;
        oadd[1] = 2 * j + 1;
        oadd[0] = CrAnMask + 2 * j + 1;
        RealOperatorAddLader(H, oadd, 4, Set.t);
      }
    }
  }
  OperatorRemoveSmallValues(H, DBL_MIN * 100);
  return 0;
}

uint32_t CreatetJModelHamiltonianMoreThenHalfFilled(tJModelSettings Set, OperatorType * H)
{
  uint32_t i, j;
  int32_t C[3];
  uint16_t oadd[6];
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
        oadd[3] = CrAnMask + 2 * i + 1;
        oadd[2] = 2 * i + 1;
        oadd[1] = CrAnMask + 2 * j + 1;
        oadd[0] = 2 * j + 1;
        RealOperatorAddLaderNormalOrder(H, oadd, 4, 0.25 * Set.J, &QDetNotOrderedFcFaBcBa);
        oadd[3] = CrAnMask + 2 * i;
        oadd[2] = 2 * i;
        oadd[1] = CrAnMask + 2 * j;
        oadd[0] = 2 * j;
        RealOperatorAddLaderNormalOrder(H, oadd, 4, 0.25 * Set.J, &QDetNotOrderedFcFaBcBa);
        oadd[3] = CrAnMask + 2 * i;
        oadd[2] = 2 * i;
        oadd[1] = CrAnMask + 2 * j + 1;
        oadd[0] = 2 * j + 1;
        RealOperatorAddLaderNormalOrder(H, oadd, 4, -0.25 * Set.J, &QDetNotOrderedFcFaBcBa);
        oadd[3] = CrAnMask + 2 * i + 1;
        oadd[2] = 2 * i + 1;
        oadd[1] = CrAnMask + 2 * j;
        oadd[0] = 2 * j;
        RealOperatorAddLaderNormalOrder(H, oadd, 4, -0.25 * Set.J, &QDetNotOrderedFcFaBcBa);
        oadd[3] = CrAnMask + 2 * i + 1;
        oadd[2] = 2 * i;
        oadd[1] = CrAnMask + 2 * j;
        oadd[0] = 2 * j + 1;
        RealOperatorAddLaderNormalOrder(H, oadd, 4, 0.5 * Set.J, &QDetNotOrderedFcFaBcBa);
        oadd[3] = CrAnMask + 2 * i;
        oadd[2] = 2 * i + 1;
        oadd[1] = CrAnMask + 2 * j + 1;
        oadd[0] = 2 * j;
        RealOperatorAddLaderNormalOrder(H, oadd, 4, 0.5 * Set.J, &QDetNotOrderedFcFaBcBa);
        /* subtract density density term */
        oadd[3] = CrAnMask + 2 * i + 1;
        oadd[2] = 2 * i + 1;
        oadd[1] = CrAnMask + 2 * j + 1;
        oadd[0] = 2 * j + 1;
        RealOperatorAddLaderNormalOrder(H, oadd, 4, -0.25 * Set.J, &QDetNotOrderedFcFaBcBa);
        oadd[3] = CrAnMask + 2 * i + 1;
        oadd[2] = 2 * i + 1;
        oadd[1] = CrAnMask + 2 * j;
        oadd[0] = 2 * j;
        RealOperatorAddLaderNormalOrder(H, oadd, 4, -0.25 * Set.J, &QDetNotOrderedFcFaBcBa);
        oadd[3] = CrAnMask + 2 * i;
        oadd[2] = 2 * i;
        oadd[1] = CrAnMask + 2 * j + 1;
        oadd[0] = 2 * j + 1;
        RealOperatorAddLaderNormalOrder(H, oadd, 4, -0.25 * Set.J, &QDetNotOrderedFcFaBcBa);
        oadd[3] = CrAnMask + 2 * i;
        oadd[2] = 2 * i;
        oadd[1] = CrAnMask + 2 * j;
        oadd[0] = 2 * j;
        RealOperatorAddLaderNormalOrder(H, oadd, 4, -0.25 * Set.J, &QDetNotOrderedFcFaBcBa);
        /* add hopping keeping double occupation constant */
        /* hop from j_up to i_up if j_dn is occupied (1-st line rest permutation in i and j and spin */
        oadd[3] = CrAnMask + 2 * i + 1;
        oadd[2] = 2 * j + 1;
        oadd[1] = CrAnMask + 2 * j;
        oadd[0] = 2 * j;
        RealOperatorAddLader(H, oadd, 4, Set.t);
        oadd[3] = CrAnMask + 2 * i;
        oadd[2] = 2 * j;
        oadd[1] = CrAnMask + 2 * j + 1;
        oadd[0] = 2 * j + 1;
        RealOperatorAddLader(H, oadd, 4, Set.t);
        oadd[3] = CrAnMask + 2 * j + 1;
        oadd[2] = 2 * i + 1;
        oadd[1] = CrAnMask + 2 * i;
        oadd[0] = 2 * i;
        RealOperatorAddLader(H, oadd, 4, Set.t);
        oadd[3] = CrAnMask + 2 * j;
        oadd[2] = 2 * i;
        oadd[1] = CrAnMask + 2 * i + 1;
        oadd[0] = 2 * i + 1;
        RealOperatorAddLader(H, oadd, 4, Set.t);
      }
    }
  }
  OperatorRemoveSmallValues(H, DBL_MIN * 100);
  return 0;
}


/*  Creates Set.Cluster.NCell operators
 *  O[0] flips the spin at site 0
 *  O[i] flips the spin at site 0 and at site i (sum normalized to sqrt(0.5))
 *  The 'flip' operator is the Sz operator
 */
uint32_t CreatetJModelSpinFlipOperatorRealSpace(tJModelSettings Set, OperatorType * O)
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
  oadd[1] = CrAnMask + 1;
  oadd[0] = 1;
  RealOperatorAddLaderNormalOrder(&O[0], oadd, 2, 0.5, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask;
  oadd[0] = 0;
  RealOperatorAddLaderNormalOrder(&O[0], oadd, 2, -0.5, &QDetNotOrderedFcFaBcBa);
  for (i = 1; i < Set.Cluster.NCell; i++)
  {
    oadd[1] = CrAnMask + 1;
    oadd[0] = 1;
    RealOperatorAddLaderNormalOrder(&O[i], oadd, 2, sqrt(0.125), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask;
    oadd[0] = 0;
    RealOperatorAddLaderNormalOrder(&O[i], oadd, 2, -sqrt(0.125), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 1 + 2 * i;
    oadd[0] = 1 + 2 * i;
    RealOperatorAddLaderNormalOrder(&O[i], oadd, 2, sqrt(0.125), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 2 * i;
    oadd[0] = 2 * i;
    RealOperatorAddLaderNormalOrder(&O[i], oadd, 2, -sqrt(0.125), &QDetNotOrderedFcFaBcBa);
  }
  return 0;
}

/*  Creates Set.Cluster.NCell operators
 *  O[0] flips the spin at site 0
 *  O[i] flips the spin at site 0 and at site i (sum normalized to sqrt(0.5))
 *  The 'flip' operator is the Sx operator
 */
uint32_t CreatetJModelSpinFlipOperatorSxRealSpace(tJModelSettings Set, OperatorType * O)
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
  oadd[1] = CrAnMask + 1;
  oadd[0] = 0;
  RealOperatorAddLaderNormalOrder(&O[0], oadd, 2, 0.5, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask;
  oadd[0] = 1;
  RealOperatorAddLaderNormalOrder(&O[0], oadd, 2, 0.5, &QDetNotOrderedFcFaBcBa);
  for (i = 1; i < Set.Cluster.NCell; i++)
  {
    oadd[1] = CrAnMask + 1;
    oadd[0] = 0;
    RealOperatorAddLaderNormalOrder(&O[i], oadd, 2, sqrt(0.125), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask;
    oadd[0] = 1;
    RealOperatorAddLaderNormalOrder(&O[i], oadd, 2, sqrt(0.125), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 1 + 2 * i;
    oadd[0] = 2 * i;
    RealOperatorAddLaderNormalOrder(&O[i], oadd, 2, sqrt(0.125), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 2 * i;
    oadd[0] = 1 + 2 * i;
    RealOperatorAddLaderNormalOrder(&O[i], oadd, 2, sqrt(0.125), &QDetNotOrderedFcFaBcBa);
  }
  return 0;
}

uint32_t CreatetJModelSpinFlipOperatorMomentumSpace(tJModelSettings Set, OperatorType * O)
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
      oadd[1] = CrAnMask + 1 + 2 * i;
      oadd[0] = 1 + 2 * i;
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
      oadd[1] = CrAnMask + 2 * i;
      oadd[0] = 2 * i;
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

uint32_t CreatetJModelSpinFlipOperatorSxMomentumSpace(tJModelSettings Set, OperatorType * O)
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
      oadd[1] = CrAnMask + 1 + 2 * i;
      oadd[0] = 2 * i;
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
      oadd[1] = CrAnMask + 2 * i;
      oadd[0] = 1 + 2 * i;
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

uint32_t CreatetJModelStartWaveFunctionOperator(tJModelSettings Set, OperatorType * O)
{
  uint16_t oadd[Set.Cluster.NCell];
  uint32_t i, j;
  InitOperator(O);
  O->NFermion = 2 * Set.Cluster.NCell;
  O->NWord = Min(O->NFermion, 1 + (O->NFermion - 1) / 8);
  if (Set.nElectron <= Set.Cluster.NCell)
  {
    InitOperatorFixedLength(O, Set.nElectron, 1);
    for (j = 0; j <= Set.Cluster.NCell; j++)
    {
      for (i = 0; i < Set.nElectron; i++)
      {
        oadd[i] = CrAnMask + (uint16_t) ((2 * (i + j) + (i % 2)) % (2 * Set.Cluster.NCell));
      }
      RealOperatorAddLaderNormalOrder(O, oadd, Set.nElectron, 1.0, &QDetNotOrderedFcFaBcBa);
    }
  }
  else
  {
    InitOperatorFixedLength(O, Set.nElectron, 1);
    for (i = 0; i < Set.Cluster.NCell; i++)
    {
      oadd[i] = CrAnMask + (uint16_t) (2 * i + ((1 + i) % 2));
    }
    for (j = 0; j < Set.Cluster.NCell; j++)
    {
      for (i = 0; i < Set.nElectron - Set.Cluster.NCell; i++)
      {
        oadd[i + Set.Cluster.NCell] = CrAnMask + (uint16_t) ((2 * (i + j) + ((i) % 2)) % (2 * Set.Cluster.NCell));
      }
      RealOperatorAddLaderNormalOrder(O, oadd, Set.nElectron, 1.0, &QDetNotOrderedFcFaBcBa);
    }
    for (i = 0; i < Set.Cluster.NCell; i++)
    {
      oadd[i] = CrAnMask + (uint16_t) (2 * i + ((i) % 2));
    }
    for (j = 0; j < Set.Cluster.NCell; j++)
    {
      for (i = 0; i < Set.nElectron - Set.Cluster.NCell; i++)
      {
        oadd[i + Set.Cluster.NCell] = CrAnMask + (uint16_t) ((2 * (i + j) + ((i) % 2)) % (2 * Set.Cluster.NCell));
      }
      RealOperatorAddLaderNormalOrder(O, oadd, Set.nElectron, 1.0, &QDetNotOrderedFcFaBcBa);
    }
  }
  return 0;
}

uint32_t CreatetJModelStartWaveFunction(tJModelSettings Set, WaveFunctionType * Psi)
{
  OperatorType O;
  WaveFunctionType f;
  uint8_t *Determinant;
  CreatetJModelStartWaveFunctionOperator(Set, &O);
  f.QComplex = 0;
  f.NFermion = 2 * Set.Cluster.NCell;
  f.NBoson = 0;
  f.NBitsKey = Set.NBitsKey;
  WaveFunctionInit(&f);
  Determinant = calloc(f.NWord, sizeof(uint8_t));
  WaveFunctionAddElement(&f, Determinant, 1.0);
  OperatorPsi(O, f, Psi);
  WaveFunctionRandomize(Psi);
  WaveFunctionNormalize(Psi);
  FreeOperator(&O);
  sprintf(Psi->Name, "Start Wavefunction");
  return 0;
}

uint32_t CreatetJModelStartWaveFunctionOperatorFM(tJModelSettings Set, OperatorType * O)
{
  uint16_t oadd[Set.nElectron];
  uint32_t i, j;
  InitOperator(O);
  O->NFermion = 2 * Set.Cluster.NCell;
  O->NWord = Min(O->NFermion, 1 + (O->NFermion - 1) / 8);
  if (Set.nElectron <= Set.Cluster.NCell)
  {
    InitOperatorFixedLength(O, Set.nElectron, 1);
    for (j = 0; j <= (Set.Cluster.NCell - Set.nElectron); j++)
    {
      for (i = 0; i < Set.nElectron; i++)
      {
        oadd[i] = CrAnMask + (uint16_t) (2 * (i + j) + 1);
      }
      RealOperatorAddLaderNormalOrder(O, oadd, Set.nElectron, 1.0, &QDetNotOrderedFcFaBcBa);
    }
  }
  else
  {
    InitOperatorFixedLength(O, Set.nElectron, 1);
    for (i = 0; i < Set.Cluster.NCell; i++)
    {
      oadd[i] = CrAnMask + (uint16_t) (2 * i + 1);
    }
    for (j = 0; j <= (2 * Set.Cluster.NCell - Set.nElectron); j++)
    {
      for (i = 0; i < Set.nElectron - Set.Cluster.NCell; i++)
      {
        oadd[i + Set.Cluster.NCell] = CrAnMask + (uint16_t) (2 * (i + j));
      }
      RealOperatorAddLaderNormalOrder(O, oadd, Set.nElectron, 1.0, &QDetNotOrderedFcFaBcBa);
    }
  }
  return 0;
}

uint32_t CreatetJModelStartWaveFunctionFM(tJModelSettings Set, WaveFunctionType * Psi)
{
  OperatorType O;
  WaveFunctionType f;
  uint8_t *Determinant;
  CreatetJModelStartWaveFunctionOperatorFM(Set, &O);
  f.QComplex = 0;
  f.NFermion = 2 * Set.Cluster.NCell;
  f.NBoson = 0;
  f.NBitsKey = Set.NBitsKey;
  WaveFunctionInit(&f);
  Determinant = calloc(f.NWord, sizeof(uint8_t));
  WaveFunctionAddElement(&f, Determinant, 1);
  OperatorPsi(O, f, Psi);
  WaveFunctionRandomize(Psi);
  WaveFunctionNormalize(Psi);
  FreeOperator(&O);
  sprintf(Psi->Name, "Start Wavefunction");
  return 0;
}





int main(int argc, const char *argv[])
{
  double TotEmax, *EnergyList, *WaveFunctionQuality;
  char   SpecFileName[255];
  double Hist[20];
  uint32_t jpsi;
  uint32_t i, j;
  argc -= 1;
  MatrixType SpinRho, Rho, RhoCorr;
  MatrixType SpinRhoTot, RhoTot, RhoCorrTot;
  WaveFunctionType *PsiList;
  tJModelSettings Set;
  WaveFunctionType OPsi;
  OperatorType H, *T;
  TriDiagonalMatrixType s;
  SpectraType *Spec;
  TimeType Time;
  InitializeTiming(&Time);
  BeginTiming("Total", &Time);
  TimeStampAll();
  FileReadtJModelSettings(&Set, argv[argc]);
  PrinttJModelSettings(Set);
  BeginTiming("Create Hamiltonian", &Time);
  if (Set.nElectron > Set.Cluster.NCell)
  {
    CreatetJModelHamiltonianMoreThenHalfFilled(Set, &H);
  }
  else if (Set.nElectron < Set.Cluster.NCell)
  {
    CreatetJModelHamiltonianLessThenHalfFilled(Set, &H);
  }
  else
  {
    CreatetJModelHamiltonian(Set, &H);
  }
  EndTiming("Create Hamiltonian", &Time);
  PrintOperatorReadable(H);
  BeginTiming("Create Transition Matrix", &Time);
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
      CreatetJModelSpinFlipOperatorRealSpace(Set, T);
    else
      CreatetJModelSpinFlipOperatorSxRealSpace(Set, T);
  }
  else
  {
    if (Set.J > 0)
      CreatetJModelSpinFlipOperatorMomentumSpace(Set, T);
    else
      CreatetJModelSpinFlipOperatorSxMomentumSpace(Set, T);
  }
  EndTiming("Create Transition Matrix", &Time);
  for (i = 0; i < Set.Cluster.NCell; i++)
  {
    PrintOperatorReadable(T[i]);
  }
  PsiList = malloc(Set.NPsi * sizeof(WaveFunctionType));
  EnergyList = malloc(Set.NPsi * sizeof(double));
  WaveFunctionQuality = malloc(Set.NPsi * sizeof(double));
  for (jpsi = 0; jpsi < Set.NPsi; jpsi++)
  {
    if (Set.J > 0)
      CreatetJModelStartWaveFunction(Set, &PsiList[jpsi]);
    else
      CreatetJModelStartWaveFunctionFM(Set, &PsiList[jpsi]);
    if (jpsi != 0)
      WaveFunctionOrthogonalize(&PsiList[jpsi], PsiList, jpsi);
    WaveFunctionRandomize(&PsiList[jpsi]);
    PrintWaveFunctionReadableTruncated(PsiList[jpsi], 0.1);
    TotEmax = 16 * Set.NE * (fabs(Set.J) + fabs(Set.t));
    BeginTiming("LanczosGroundState", &Time);
    LanczosGroundState(H, &PsiList[jpsi], &Set.zeroLancWavefunction, Set.zeroLancEnergy, TotEmax, PsiList, jpsi);
    EndTiming("LanczosGroundState", &Time);
    TotEmax *= 2;
    BeginTiming("LanczosGroundState", &Time);
    LanczosGroundState(H, &PsiList[jpsi], &Set.zeroLancWavefunction, Set.zeroLancEnergy, TotEmax, PsiList, jpsi);
    EndTiming("LanczosGroundState", &Time);
    sprintf(PsiList[jpsi].Name, "Eigenstate #%i", jpsi);
    if (OperatorPsiConserveBasis(H, PsiList[jpsi], &OPsi))
      AllocFailed();
    EnergyList[jpsi] = WaveFunctionDotProduct(OPsi, PsiList[jpsi]);
    PsiOperatorSqrPsiRR(H, PsiList[jpsi], &WaveFunctionQuality[jpsi]);
    WaveFunctionQuality[jpsi] /= Sqr(EnergyList[jpsi]);
    WaveFunctionFreeValues(&OPsi);
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
    CalculateSpinCorrelationMatrix(&SpinRho, PsiList[jpsi]);
    PrintMatrix(SpinRho);
    CalculateDensityMatrix(&Rho, PsiList[jpsi]);
    PrintMatrix(Rho);
    CalculateDensityCorrelationMatrix(&RhoCorr, PsiList[jpsi]);
    PrintMatrix(RhoCorr);
    if ((EnergyList[jpsi] - EnergyList[0]) > 1E-5)
    {
      printf("Average of Density matrices for state 0 to %i\n", jpsi - 1);
      MatrixScalarMultiply(&SpinRhoTot, 1.0 / ((double)(jpsi)), 0);
      MatrixScalarMultiply(&RhoTot, 1.0 / ((double)(jpsi)), 0);
      MatrixScalarMultiply(&RhoCorrTot, 1.0 / ((double)(jpsi)), 0);
      PrintMatrix(SpinRhoTot);
      PrintMatrix(RhoTot);
      PrintMatrix(RhoCorrTot);
      EndTiming("Total", &Time);
      PrintTiming(Time);
      return 0;
    }
    if (jpsi == 0)
    {
      MatrixCopy(&SpinRhoTot, SpinRho);
      MatrixCopy(&RhoTot, Rho);
      MatrixCopy(&RhoCorrTot, RhoCorr);
    }
    else
    {
      MatrixAdd(&SpinRhoTot, SpinRho, 1, 0);
      MatrixAdd(&RhoTot, Rho, 1, 0);
      MatrixAdd(&RhoCorrTot, RhoCorr, 1, 0);
    }
    FreeMatrix(&SpinRho);
    FreeMatrix(&Rho);
    FreeMatrix(&RhoCorr);
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
    if (jpsi < 10)
      sprintf(SpecFileName, "%s00%1i", Set.SpecFile, jpsi);
    else if (jpsi < 100)
      sprintf(SpecFileName, "%s0%2i", Set.SpecFile, jpsi);
    else if (jpsi < 1000)
      sprintf(SpecFileName, "%s%3i", Set.SpecFile, jpsi);
    else
      sprintf(SpecFileName, "%s%i", Set.SpecFile, jpsi);
    FilePrintSpectra(Spec, Set.Cluster.NCell, SpecFileName, "w");
  }
  EndTiming("Total", &Time);
  PrintTiming(Time);
  return 0;
}
