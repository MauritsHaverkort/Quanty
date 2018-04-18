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

/*  Basis on a rung
*   0 (ud-du) * sqrt(1/2)
*   1 dd
*   2 (ud+du) * sqrt(1/2)
*   3 uu
*
*   Some shorthand
*   d1 == a^{\dag}_1
*   a1 == a_1
*
*   some spin operator relations
*   S^2 |S,M> = s(s+1) |S,M>
*   Sz  |S,M> = M      |S,M>
*   S+  |S,M> = sqrt[s(s+1)-m(m+1)] |S,M+1>
*   S-  |S,M> = sqrt[s(s+1)-m(m-1)] |S,M-1>
*   S+ = Sx + i Sy
*   S- = Sx - i Sy
*   Sx = 1/2   (S- + S+)
*   Sy = 1/2 i (S- - S+)
*   Si . Sj = Szi Szj + 1/2 (Si+ Sj- + Si- Sj+)
*   
*   Some operators
*
*   2 Sz1   =  0  0  1  0  =  d0 a2
*              0 -1  0  0  -  d1 a1
*              1  0  0  0  +  d2 a0
*              0  0  0  1  +  d3 a3
*
*   2 Sz2   =  0  0 -1  0  = -d0 a2
*              0 -1  0  0  -  d1 a1
*             -1  0  0  0  -  d2 a0
*              0  0  0  1  +  d3 a3
*
*   t=sqrt(1/2)
*
*   S-1     =  0  0  0 -t  = -d0 a3
*              t  0  t  0  +  d1 a0 + d1 a2              
*              0  0  0  t  +  d2 a3
*              0  0  0  0
*
*   S+1     =  0  t  0  0  =  d0 a1
*              0  0  0  0
*              0  t  0  0  +  d2 a1
*             -t  0  t  0  -  d3 a0 + d3 a2
*
*   S-2     =  0  0  0  t  =  d0 a3
*             -t  0  t  0  -  d1 a0 + d1 a2
*              0  0  0  t  +  d2 a3
*              0  0  0  0
*
*   S+2     =  0 -t  0  0  = -d0 a1
*              0  0  0  0
*              0  t  0  0  +  d2 a1
*              t  0  t  0  +  d3 a0 + d3 a2
*
*   Sz1+Sz2 =  0  0  0  0  =
*              0 -1  0  0  -  d1 a1
*              0  0  0  0
*              0  0  0  1  +  d3 a3
*
*   Sz1-Sz2 =  0  0  1  0  =  d0 a2
*              0  0  0  0
*              1  0  0  0  +  d2 a0
*              0  0  0  0
*
*   Composite Operators
*
*   Si1 . Sj1 = 1/4 ( d0i a2i - d1i a1i + d2i a0i + d3i a3i )( d0j a2j - d1j a1j + d2j a0j + d3j a3j )
*             + 1/4 (-d0i a3i + d1i a0i + d1i a2i + d2i a3i )( d0j a1j + d2j a1j - d3j a0j + d3j a2j )
*             + 1/4 ( d0i a1i + d2i a1i - d3i a0i + d3i a2i )(-d0j a3j + d1j a0j + d1j a2j + d2j a3j )
*/


uint32_t CreateSpinLaderSpinOperatorsSz(SpinLaderSettings Set, OperatorType * S)
{
  uint32_t i;
  uint16_t oadd[2];
  for (i = 0; i < 2 * Set.NSite; i++)
  {
    InitNamedOperator(&S[i], "Sz");
    S[i].NFermion = 4 * Set.NSite;
    S[i].NWord = Min(S[i].NFermion, (1 + S[i].NFermion - 1) / 8);
    InitOperatorFixedLength(&S[i], 2, 8);
  }
  for (i = 0; i < Set.NSite; i++)
  {                             /* acting on left spin on a rung */
    oadd[1] = CrAnMask + 4 * i + 2;
    oadd[0] = 4 * i + 0;
    RealOperatorAddLaderNormalOrder(&S[2 * i], oadd, 2, 0.5, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4 * i + 0;
    oadd[0] = 4 * i + 2;
    RealOperatorAddLaderNormalOrder(&S[2 * i], oadd, 2, 0.5, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4 * i + 3;
    oadd[0] = 4 * i + 3;
    RealOperatorAddLaderNormalOrder(&S[2 * i], oadd, 2, 0.5, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4 * i + 1;
    oadd[0] = 4 * i + 1;
    RealOperatorAddLaderNormalOrder(&S[2 * i], oadd, 2, -0.5, &QDetNotOrderedFcFaBcBa);
    /* acting on right spin on a rung */
    oadd[1] = CrAnMask + 4 * i + 2;
    oadd[0] = 4 * i + 0;
    RealOperatorAddLaderNormalOrder(&S[2 * i + 1], oadd, 2, -0.5, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4 * i + 0;
    oadd[0] = 4 * i + 2;
    RealOperatorAddLaderNormalOrder(&S[2 * i + 1], oadd, 2, -0.5, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4 * i + 3;
    oadd[0] = 4 * i + 3;
    RealOperatorAddLaderNormalOrder(&S[2 * i + 1], oadd, 2, 0.5, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4 * i + 1;
    oadd[0] = 4 * i + 1;
    RealOperatorAddLaderNormalOrder(&S[2 * i + 1], oadd, 2, -0.5, &QDetNotOrderedFcFaBcBa);
  }
  return 0;
}

uint32_t CreateSpinLaderSpinOperatorsSplus(SpinLaderSettings Set, OperatorType * S)
{
  uint32_t i;
  uint16_t oadd[2];
  for (i = 0; i < 2 * Set.NSite; i++)
  {
    InitNamedOperator(&S[i], "S plus");
    S[i].NFermion = 4 * Set.NSite;
    S[i].NWord = Min(S[i].NFermion, (1 + S[i].NFermion - 1) / 8);
    InitOperatorFixedLength(&S[i], 2, 8);
  }
  for (i = 0; i < Set.NSite; i++)
  {                             /* acting on left spin on a rung */
    oadd[1] = CrAnMask + 4 * i + 3;
    oadd[0] = 4 * i + 0;
    RealOperatorAddLaderNormalOrder(&S[2 * i], oadd, 2, -sqrt(0.5), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4 * i + 3;
    oadd[0] = 4 * i + 2;
    RealOperatorAddLaderNormalOrder(&S[2 * i], oadd, 2, sqrt(0.5), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4 * i + 0;
    oadd[0] = 4 * i + 1;
    RealOperatorAddLaderNormalOrder(&S[2 * i], oadd, 2, sqrt(0.5), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4 * i + 2;
    oadd[0] = 4 * i + 1;
    RealOperatorAddLaderNormalOrder(&S[2 * i], oadd, 2, sqrt(0.5), &QDetNotOrderedFcFaBcBa);
    /* acting on right spin on a rung */
    oadd[1] = CrAnMask + 4 * i + 3;
    oadd[0] = 4 * i + 0;
    RealOperatorAddLaderNormalOrder(&S[2 * i + 1], oadd, 2, sqrt(0.5), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4 * i + 3;
    oadd[0] = 4 * i + 2;
    RealOperatorAddLaderNormalOrder(&S[2 * i + 1], oadd, 2, sqrt(0.5), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4 * i + 0;
    oadd[0] = 4 * i + 1;
    RealOperatorAddLaderNormalOrder(&S[2 * i + 1], oadd, 2, -sqrt(0.5), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4 * i + 2;
    oadd[0] = 4 * i + 1;
    RealOperatorAddLaderNormalOrder(&S[2 * i + 1], oadd, 2, sqrt(0.5), &QDetNotOrderedFcFaBcBa);
  }
  return 0;
}

uint32_t CreateSpinLaderSpinOperatorsSminus(SpinLaderSettings Set, OperatorType * S)
{
  uint32_t i;
  uint16_t oadd[2];
  for (i = 0; i < 2 * Set.NSite; i++)
  {
    InitNamedOperator(&S[i], "Sminus");
    S[i].NFermion = 4 * Set.NSite;
    S[i].NWord = Min(S[i].NFermion, (1 + S[i].NFermion - 1) / 8);
    InitOperatorFixedLength(&S[i], 2, 8);
  }
  for (i = 0; i < Set.NSite; i++)
  {                             /* acting on left spin on a rung */
    oadd[1] = CrAnMask + 4 * i + 1;
    oadd[0] = 4 * i + 0;
    RealOperatorAddLaderNormalOrder(&S[2 * i], oadd, 2, sqrt(0.5), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4 * i + 1;
    oadd[0] = 4 * i + 2;
    RealOperatorAddLaderNormalOrder(&S[2 * i], oadd, 2, sqrt(0.5), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4 * i + 0;
    oadd[0] = 4 * i + 3;
    RealOperatorAddLaderNormalOrder(&S[2 * i], oadd, 2, -sqrt(0.5), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4 * i + 2;
    oadd[0] = 4 * i + 3;
    RealOperatorAddLaderNormalOrder(&S[2 * i], oadd, 2, sqrt(0.5), &QDetNotOrderedFcFaBcBa);
    /* acting on right spin on a rung */
    oadd[1] = CrAnMask + 4 * i + 1;
    oadd[0] = 4 * i + 0;
    RealOperatorAddLaderNormalOrder(&S[2 * i + 1], oadd, 2, -sqrt(0.5), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4 * i + 1;
    oadd[0] = 4 * i + 2;
    RealOperatorAddLaderNormalOrder(&S[2 * i + 1], oadd, 2, sqrt(0.5), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4 * i + 0;
    oadd[0] = 4 * i + 3;
    RealOperatorAddLaderNormalOrder(&S[2 * i + 1], oadd, 2, sqrt(0.5), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4 * i + 2;
    oadd[0] = 4 * i + 3;
    RealOperatorAddLaderNormalOrder(&S[2 * i + 1], oadd, 2, sqrt(0.5), &QDetNotOrderedFcFaBcBa);
  }
  return 0;
}

/*  Within the basis set fermion 0 to 3 is single occupied
 *  Operators of the form A1 A0 are thus needles
 *  remove them by setting their value to zero and runing remove opp. zero values
 */
uint32_t SetValuesOfUnnededOperatorsToZero(OperatorType * H)
{
  uint32_t i, j, k, l;
  for (i = 0; i <= H->MaxLength; i++)
  {
    for (j = 0; j < H->l[i].N; j++)
    {
      for (k = 0; k < i; k++)
      {
        for (l = k + 1; l < i; l++)
        {
          if ((H->l[i].o[i * j + k] >> 2) == (H->l[i].o[i * j + l] >> 2))
          {
            H->l[i].VRe[j] = 0;
          }
        }
      }
    }
  }
  return 0;
}

uint32_t CreateSpinLaderHamiltonian(SpinLaderSettings Set, OperatorType * H)
{
  uint32_t i, j;
  uint16_t oadd[2];
  OperatorType *Sz, *Splus, *Smin;
  OperatorType Hrung, Hleg, Hcyc;
  OperatorType OtempA, OtempB;
  Sz = malloc(Set.NSite * 2 * sizeof(OperatorType));
  Splus = malloc(Set.NSite * 2 * sizeof(OperatorType));
  Smin = malloc(Set.NSite * 2 * sizeof(OperatorType));
  CreateSpinLaderSpinOperatorsSz(Set, Sz);
  CreateSpinLaderSpinOperatorsSplus(Set, Splus);
  CreateSpinLaderSpinOperatorsSminus(Set, Smin);

  InitNamedOperator(H, "Hamiltonian");
  H->NFermion = 4 * Set.NSite;
  H->NWord = Min(H->NFermion, 1 + (H->NFermion - 1) / 8);
  InitOperatorFixedLength(H, 2, 8);
  InitOperatorFixedLength(H, 4, 8);
  InitOperatorFixedLength(H, 6, 8);
  InitOperatorFixedLength(H, 8, 8);
  InitNamedOperator(&Hleg, "H leg");
  Hleg.NFermion = 4 * Set.NSite;
  Hleg.NWord = Min(Hleg.NFermion, 1 + (Hleg.NFermion - 1) / 8);
  InitOperatorFixedLength(&Hleg, 2, 8);
  InitOperatorFixedLength(&Hleg, 4, 8);
  InitNamedOperator(&Hcyc, "H cyc");
  Hcyc.NFermion = 4 * Set.NSite;
  Hcyc.NWord = Min(Hcyc.NFermion, 1 + (Hcyc.NFermion - 1) / 8);
  InitOperatorFixedLength(&Hcyc, 2, 8);
  InitOperatorFixedLength(&Hcyc, 4, 8);
  InitOperatorFixedLength(&Hcyc, 6, 8);
  InitOperatorFixedLength(&Hcyc, 8, 8);
  InitNamedOperator(&Hrung, "H rung");
  Hrung.NFermion = 4 * Set.NSite;
  Hrung.NWord = Min(Hrung.NFermion, 1 + (Hrung.NFermion - 1) / 8);
  InitOperatorFixedLength(&Hrung, 2, 8);
  InitOperatorFixedLength(&Hrung, 4, 8);
  for (i = 0; i < Set.NSite; i++)
  {
    j = (i + 1) % Set.NSite;
    /* add Si_1 . Si_2 to the Hamiltonian */
    /* The following definition is correct, but leaves to much junk, i.e. operators that do not do anything
       OperatorMultiplyNormalOrder(&Hrung, Sz[2*i], Sz[2*i+1], &QDetNotOrderedFcFaBcBa);
       OperatorMultiplyNormalOrder(&Hrung, Sz[2*i], Sz[2*i+1], &QDetNotOrderedFcFaBcBa);
       OperatorMultiplyNormalOrder(&Hrung, Splus[2*i], Smin[2*i+1], &QDetNotOrderedFcFaBcBa);
       OperatorMultiplyNormalOrder(&Hrung, Smin[2*i], Splus[2*i+1], &QDetNotOrderedFcFaBcBa); */
    oadd[1] = CrAnMask + 4 * i;
    oadd[0] = 4 * i;
    RealOperatorAddLaderNormalOrder(&Hrung, oadd, 2, -1.00 * Set.Jrung, &QDetNotOrderedFcFaBcBa);
    RealOperatorAddLaderNormalOrder(&Hrung, oadd, 0, 0.25 * Set.Jrung, &QDetNotOrderedFcFaBcBa);
    /* add 2 Si_1 . Sj_1 to the Hamiltonian */
    OperatorMultiplyNormalOrder(&Hleg, Sz[2 * i], Sz[2 * j], &QDetNotOrderedFcFaBcBa);
    OperatorMultiplyNormalOrder(&Hleg, Sz[2 * i], Sz[2 * j], &QDetNotOrderedFcFaBcBa);
    OperatorMultiplyNormalOrder(&Hleg, Splus[2 * i], Smin[2 * j], &QDetNotOrderedFcFaBcBa);
    OperatorMultiplyNormalOrder(&Hleg, Smin[2 * i], Splus[2 * j], &QDetNotOrderedFcFaBcBa);
    /* add 2 Si_2 . Sj_2 to the Hamiltonian */
    OperatorMultiplyNormalOrder(&Hleg, Sz[2 * i + 1], Sz[2 * j + 1], &QDetNotOrderedFcFaBcBa);
    OperatorMultiplyNormalOrder(&Hleg, Sz[2 * i + 1], Sz[2 * j + 1], &QDetNotOrderedFcFaBcBa);
    OperatorMultiplyNormalOrder(&Hleg, Splus[2 * i + 1], Smin[2 * j + 1], &QDetNotOrderedFcFaBcBa);
    OperatorMultiplyNormalOrder(&Hleg, Smin[2 * i + 1], Splus[2 * j + 1], &QDetNotOrderedFcFaBcBa);
    /* add (Si_1 . Sj_1)(Si_2 . Sj_2) to the Hamiltonian */
    InitNamedOperator(&OtempA, "OtempA");
    OtempA.NFermion = 4 * Set.NSite;
    OtempA.NWord = Min(OtempA.NFermion, 1 + (OtempA.NFermion - 1) / 8);
    InitOperatorFixedLength(&OtempA, 2, 8);
    InitOperatorFixedLength(&OtempA, 4, 8);
    OperatorMultiplyNormalOrder(&OtempA, Sz[2 * i], Sz[2 * j], &QDetNotOrderedFcFaBcBa);
    OperatorMultiplyScaledNormalOrder(&OtempA, Splus[2 * i], Smin[2 * j], 0.5, 0, &QDetNotOrderedFcFaBcBa);
    OperatorMultiplyScaledNormalOrder(&OtempA, Smin[2 * i], Splus[2 * j], 0.5, 0, &QDetNotOrderedFcFaBcBa);
    InitNamedOperator(&OtempB, "OtempB");
    OtempB.NFermion = 4 * Set.NSite;
    OtempB.NWord = Min(OtempB.NFermion, 1 + (OtempB.NFermion - 1) / 8);
    InitOperatorFixedLength(&OtempB, 2, 8);
    InitOperatorFixedLength(&OtempB, 4, 8);
    oadd[1] = CrAnMask + 4 * i;
    oadd[0] = 4 * i;
    OperatorMultiplyNormalOrder(&OtempB, Sz[2 * i + 1], Sz[2 * j + 1], &QDetNotOrderedFcFaBcBa);
    OperatorMultiplyScaledNormalOrder(&OtempB, Splus[2 * i + 1], Smin[2 * j + 1], 0.5, 0, &QDetNotOrderedFcFaBcBa);
    OperatorMultiplyScaledNormalOrder(&OtempB, Smin[2 * i + 1], Splus[2 * j + 1], 0.5, 0, &QDetNotOrderedFcFaBcBa);
    OperatorMultiplyNormalOrder(&Hcyc, OtempA, OtempB, &QDetNotOrderedFcFaBcBa);
    FreeOperator(&OtempA);
    FreeOperator(&OtempB);
    /* add (Si_1 . Si_2)(Sj_1 . Sj_2) to the Hamiltonian */
    InitNamedOperator(&OtempA, "OtempA");
    OtempA.NFermion = 4 * Set.NSite;
    OtempA.NWord = Min(OtempA.NFermion, 1 + (OtempA.NFermion - 1) / 8);
    InitOperatorFixedLength(&OtempA, 2, 8);
    InitOperatorFixedLength(&OtempA, 4, 8);
    oadd[1] = CrAnMask + 4 * i;
    oadd[0] = 4 * i;
    RealOperatorAddLaderNormalOrder(&OtempA, oadd, 2, -1.00, &QDetNotOrderedFcFaBcBa);
    RealOperatorAddLaderNormalOrder(&OtempA, oadd, 0, 0.25, &QDetNotOrderedFcFaBcBa);
    InitNamedOperator(&OtempB, "OtempB");
    OtempB.NFermion = 4 * Set.NSite;
    OtempB.NWord = Min(OtempB.NFermion, 1 + (OtempB.NFermion - 1) / 8);
    InitOperatorFixedLength(&OtempB, 2, 8);
    InitOperatorFixedLength(&OtempB, 4, 8);
    oadd[1] = CrAnMask + 4 * i;
    oadd[0] = 4 * i;
    RealOperatorAddLaderNormalOrder(&OtempB, oadd, 2, -1.00, &QDetNotOrderedFcFaBcBa);
    RealOperatorAddLaderNormalOrder(&OtempB, oadd, 0, 0.25, &QDetNotOrderedFcFaBcBa);
    OperatorMultiplyNormalOrder(&Hcyc, OtempA, OtempB, &QDetNotOrderedFcFaBcBa);
    FreeOperator(&OtempA);
    FreeOperator(&OtempB);
    /* add -(Si_1 . Sj_2)(Si_2 . Sj_1) to the Hamiltonian */
    InitNamedOperator(&OtempA, "OtempA");
    OtempA.NFermion = 4 * Set.NSite;
    OtempA.NWord = Min(OtempA.NFermion, 1 + (OtempA.NFermion - 1) / 8);
    InitOperatorFixedLength(&OtempA, 2, 8);
    InitOperatorFixedLength(&OtempA, 4, 8);
    OperatorMultiplyScaledNormalOrder(&OtempA, Sz[2 * i], Sz[2 * j + 1], -1.0, 0, &QDetNotOrderedFcFaBcBa);
    OperatorMultiplyScaledNormalOrder(&OtempA, Splus[2 * i], Smin[2 * j + 1], -0.5, 0, &QDetNotOrderedFcFaBcBa);
    OperatorMultiplyScaledNormalOrder(&OtempA, Smin[2 * i], Splus[2 * j + 1], -0.5, 0, &QDetNotOrderedFcFaBcBa);
    InitNamedOperator(&OtempB, "OtempB");
    OtempB.NFermion = 4 * Set.NSite;
    OtempB.NWord = Min(OtempB.NFermion, 1 + (OtempB.NFermion - 1) / 8);
    InitOperatorFixedLength(&OtempB, 2, 8);
    InitOperatorFixedLength(&OtempB, 4, 8);
    oadd[1] = CrAnMask + 4 * i;
    oadd[0] = 4 * i;
    OperatorMultiplyScaledNormalOrder(&OtempB, Sz[2 * i + 1], Sz[2 * j], -1.0, 0, &QDetNotOrderedFcFaBcBa);
    OperatorMultiplyScaledNormalOrder(&OtempB, Splus[2 * i + 1], Smin[2 * j], -0.5, 0, &QDetNotOrderedFcFaBcBa);
    OperatorMultiplyScaledNormalOrder(&OtempB, Smin[2 * i + 1], Splus[2 * j], -0.5, 0, &QDetNotOrderedFcFaBcBa);
    OperatorMultiplyNormalOrder(&Hcyc, OtempA, OtempB, &QDetNotOrderedFcFaBcBa);
    FreeOperator(&OtempA);
    FreeOperator(&OtempB);
  }
  OperatorScale(&Hcyc, Set.Jcyc, 0);
  OperatorScale(&Hleg, 0.5 * Set.Jleg, 0);
  OperatorAdd(H, Hrung, 1, 0);
  OperatorAdd(H, Hleg, 1, 0);
  OperatorAdd(H, Hcyc, 1, 0);
  SetValuesOfUnnededOperatorsToZero(H);
  OperatorRemoveSmallValues(H, DBL_MIN * 100);
  FreeOperator(&Hrung);
  FreeOperator(&Hleg);
  FreeOperator(&Hcyc);
  return 0;
}

uint32_t CreateSpinLaderStartWaveFunction(const SpinLaderSettings Set, WaveFunctionType * psi)
{
  uint8_t *Determinant;
  int    i;
  psi->QComplex = 0;
  psi->NBitsKey = Set.NBitsKey;
  psi->NFermion = Set.NSite * 4;
  psi->NBoson = 0;
  psi->NWord = Min(psi->NFermion, 1 + (psi->NFermion - 1) / 8);
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

  for (i = 0; i < psi->NFermion; i += 4)
  {
    OperatorCr(Determinant, i, psi->NFermion);
  }
  WaveFunctionAddElement(psi, Determinant, 1.0);

  sprintf(psi->Name, "Start Wavefunction");
  return 0;
}

/*  Create for each qz 4 operators:
 *  Sz_1 + Sz_2
 *  Sz_1 - Sz_2
 *  Sz_1 + i Sz_2
 *  Sz_1 - i Sz_2
 */
uint32_t CreateSpinLaderSpinFlipOperator(SpinLaderSettings Set, OperatorType * O)
{
  uint32_t i, j, k;
  uint16_t oadd[2];
  double phase;
  for (j = 0; j < 4; j++)
  {
    for (i = 0; i < Set.NSite; i++)
    {
      InitOperator(&O[i + j * Set.NSite]);
      O[i + j * Set.NSite].NFermion = 4 * Set.NSite;
      O[i + j * Set.NSite].NWord = Min(O[i + j * Set.NSite].NFermion, 1 + (O[i + j * Set.NSite].NFermion - 1) / 8);
      InitOperatorFixedLength(&O[i + j * Set.NSite], 2, 8);
      switch (j)
      {
      case 0:
        sprintf(O[i + j * Set.NSite].Name, "Operator Sz_1 +   Sz_2 at k: %i ", i);
        break;
      case 1:
        sprintf(O[i + j * Set.NSite].Name, "Operator Sz_1 + i Sz_2 at k: %i ", i);
        break;
      case 2:
        sprintf(O[i + j * Set.NSite].Name, "Operator Sz_1 -   Sz_2 at k: %i ", i);
        break;
      case 3:
        sprintf(O[i + j * Set.NSite].Name, "Operator Sz_1 - i Sz_2 at k: %i ", i);
        break;
      }
    }
  }
  for (j = 0; j < 4; j++)
  {
    for (k = 0; k < Set.NSite; k++)
    {
      for (i = 0; i < Set.NSite; i++)
      {
        phase = 2 * Pi * (((double)(i * k)) / ((double)Set.NSite));
        oadd[1] = CrAnMask + 4 * i + 2;
        oadd[0] = 4 * i + 0;
        OperatorAddLaderNormalOrder(&O[k + j * Set.NSite], oadd, 2, 0.5 * cos(phase), 0.5 * sin(phase),
                                    &QDetNotOrderedFcFaBcBa);
        oadd[1] = CrAnMask + 4 * i + 0;
        oadd[0] = 4 * i + 2;
        OperatorAddLaderNormalOrder(&O[k + j * Set.NSite], oadd, 2, 0.5 * cos(phase), 0.5 * sin(phase),
                                    &QDetNotOrderedFcFaBcBa);
        oadd[1] = CrAnMask + 4 * i + 3;
        oadd[0] = 4 * i + 3;
        OperatorAddLaderNormalOrder(&O[k + j * Set.NSite], oadd, 2, 0.5 * cos(phase), 0.5 * sin(phase),
                                    &QDetNotOrderedFcFaBcBa);
        oadd[1] = CrAnMask + 4 * i + 1;
        oadd[0] = 4 * i + 1;
        OperatorAddLaderNormalOrder(&O[k + j * Set.NSite], oadd, 2, -0.5 * cos(phase), -0.5 * sin(phase),
                                    &QDetNotOrderedFcFaBcBa);
        /* acting on right spin on a rung */
        phase = 2 * Pi * (((double)(i * k)) / ((double)Set.NSite) + ((double)j) / 4.0);
        oadd[1] = CrAnMask + 4 * i + 2;
        oadd[0] = 4 * i + 0;
        OperatorAddLaderNormalOrder(&O[k + j * Set.NSite], oadd, 2, -0.5 * cos(phase), -0.5 * sin(phase),
                                    &QDetNotOrderedFcFaBcBa);
        oadd[1] = CrAnMask + 4 * i + 0;
        oadd[0] = 4 * i + 2;
        OperatorAddLaderNormalOrder(&O[k + j * Set.NSite], oadd, 2, -0.5 * cos(phase), -0.5 * sin(phase),
                                    &QDetNotOrderedFcFaBcBa);
        oadd[1] = CrAnMask + 4 * i + 3;
        oadd[0] = 4 * i + 3;
        OperatorAddLaderNormalOrder(&O[k + j * Set.NSite], oadd, 2, 0.5 * cos(phase), 0.5 * sin(phase),
                                    &QDetNotOrderedFcFaBcBa);
        oadd[1] = CrAnMask + 4 * i + 1;
        oadd[0] = 4 * i + 1;
        OperatorAddLaderNormalOrder(&O[k + j * Set.NSite], oadd, 2, -0.5 * cos(phase), -0.5 * sin(phase),
                                    &QDetNotOrderedFcFaBcBa);
      }
      OperatorRemoveSmallValues(&O[k + j * Set.NSite], DBL_EPSILON * 100);
    }
  }
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
  MatrixType Rho, RhoCorr;
  MatrixType RhoTot, RhoCorrTot;
  WaveFunctionType *PsiList;
  SpinLaderSettings Set;
  WaveFunctionType OPsi;
  OperatorType H;
  TimeType Time;
  TriDiagonalMatrixType s;
  OperatorType *T;
  SpectraType *Spec;
  InitializeTiming(&Time);
  BeginTiming("Total", &Time);
  TimeStampAll();
  FileReadSpinLaderSettings(&Set, argv[argc]);
  PrintSpinLaderSettings(Set);
  BeginTiming("Create Hamiltonian", &Time);
  CreateSpinLaderHamiltonian(Set, &H);
  EndTiming("Create Hamiltonian", &Time);
  PrintOperatorReadable(H);
  BeginTiming("Create Transition Matrix", &Time);
  s.N = Set.NTriDiag;
  InitTriDiagonalMatrix(&s);
  T = malloc(Set.NSite * 4 * sizeof(OperatorType));
  Spec = malloc(Set.NSite * 4 * sizeof(SpectraType));
  Spec[0].N = Set.NE;
  Spec[0].Emin = Set.Emin;
  Spec[0].Emax = Set.Emax;
  Spec[0].dE = (Set.Emax - Set.Emin) / Set.NE;
  Spec[0].EminPole = Set.Emin;
  Spec[0].EmaxPole = Set.Emax;
  Spec[0].Gamma = Set.Gamma;
  for (i = 0; i < Set.NSite * 4; i++)
  {
    SpectrumAllocate(Spec[0], &Spec[i]);
  }
  CreateSpinLaderSpinFlipOperator(Set, T);
  EndTiming("Create Transition Matrix", &Time);
  for (i = 0; i < Set.NSite * 4; i++)
  {
    PrintOperatorReadable(T[i]);
  }
  PsiList = malloc(Set.NPsi * sizeof(WaveFunctionType));
  EnergyList = malloc(Set.NPsi * sizeof(double));
  WaveFunctionQuality = malloc(Set.NPsi * sizeof(double));
  for (jpsi = 0; jpsi < Set.NPsi; jpsi++)
  {
    CreateSpinLaderStartWaveFunction(Set, &PsiList[jpsi]);
    if (jpsi != 0)
      WaveFunctionOrthogonalize(&PsiList[jpsi], PsiList, jpsi);
    WaveFunctionRandomize(&PsiList[jpsi]);
    PrintWaveFunctionReadableTruncated(PsiList[jpsi], 0.1);
    TotEmax = 4 * Set.NSite * (fabs(Set.Jcyc) + fabs(Set.Jrung) + fabs(Set.Jleg));
    BeginTiming("LanczosGroundState", &Time);
    LanczosGroundState(H, &PsiList[jpsi], &Set.zeroLancWavefunction, Set.zeroLancEnergy, TotEmax, PsiList, jpsi);
    EndTiming("LanczosGroundState", &Time);
    TotEmax *= 2;
    BeginTiming("LanczosGroundState", &Time);
    LanczosGroundState(H, &PsiList[jpsi], &Set.zeroLancWavefunction, Set.zeroLancEnergy, TotEmax, PsiList, jpsi);
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
    CalculateDensityMatrix(&Rho, PsiList[jpsi]);
    PrintMatrix(Rho);
    CalculateDensityCorrelationMatrix(&RhoCorr, PsiList[jpsi]);
    PrintMatrix(RhoCorr);
    if ((EnergyList[jpsi] - EnergyList[0]) > 1E-5)
    {
      printf("Average of Density matrices for state 0 to %i\n", jpsi - 1);
      MatrixScalarMultiply(&RhoTot, 1.0 / ((double)(jpsi)), 0);
      MatrixScalarMultiply(&RhoCorrTot, 1.0 / ((double)(jpsi)), 0);
      PrintMatrix(RhoTot);
      PrintMatrix(RhoCorrTot);
      EndTiming("Total", &Time);
      PrintTiming(Time);
      return 0;
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
    for (i = 0; i < Set.NSite * 4; i++)
    {
      BeginTiming("Lanczos Tri Diagonalization", &Time);
      OperatorPsi(T[i], PsiList[jpsi], &OPsi);
      LanczosTriDiagonalize(H, OPsi, Set.zeroLancWavefunction, &s, Set.Gamma);
      TriDiagonalMatrixToSpectrum(EnergyList[jpsi], s, &Spec[i]);
      EndTiming("Lanczos Tri Diagonalization", &Time);
    }
    if (jpsi < 10)
      sprintf(SpecFileName, "%s00%1i", Set.SpecFile, jpsi);
    else if (jpsi < 100)
      sprintf(SpecFileName, "%s0%2i", Set.SpecFile, jpsi);
    else if (jpsi < 1000)
      sprintf(SpecFileName, "%s%3i", Set.SpecFile, jpsi);
    else
      sprintf(SpecFileName, "%s%i", Set.SpecFile, jpsi);
    FilePrintSpectra(Spec, 4 * Set.NSite, SpecFileName, "w");

  }
  EndTiming("Total", &Time);
  PrintTiming(Time);
  return 0;
}
