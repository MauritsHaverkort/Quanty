/* ===========================================================================
 *
 *   G(w) for a Hamiltonian H in a space P plus hybridisation functions, with
 *   up to NExcitations excitations at once in the baths they stand for
 *
 *   CalculateGforHwithSigmaAsOperator allows one excitation.  This file lifts
 *   that restriction.  The routine below reduces to the old one for
 *   NExcitations = 1 (same space, same answer), but CreateSpectra keeps calling
 *   the old routine there, so that nothing changes for existing scripts.
 *
 *   Reference: Y. Lu, X. Cao, P. Hansmann and M. W. Haverkort, Phys. Rev. B
 *   100, 115134 (2019), section II C.  Their p is NExcitations here.  Their
 *   chain sites beyond L are the poles of the hybridisation functions here, and
 *   their Krylov space K^M of the subsystem is PsiBasis.
 *
 * =========================================================================== */

/*  the bath modes: every pole of every hybridisation function factorized into
 *  modes of rank one.  Mode mu belongs to entry Entry[mu], sits at energy E[mu]
 *  and couples to channel j of its entry with L[mu * BMax + j].  The modes of
 *  entry i are First[i] .. First[i+1]-1, so the global index is entry major    */
typedef struct BathExcModesType
{
  uint32_t N;
  uint32_t BMax;
  uint32_t *Entry;
  uint32_t *First;
  uint32_t *Parity;       /* per entry, 1 fermionic, 0 bosonic */
  double *E;
  Complex16 *L;
}
BathExcModesType;

/*  one system space X_Level^S and the block of the auxiliary Hamiltonian built
 *  on it.  Level 0 is P itself.  For Level > 0 the space is spanned by
 *
 *    V[s_Level] ... V[s_2] V[s_1] | e_p >
 *
 *  for every string s of (entry,channel) pairs whose entries form the multiset
 *  S = Entry[0..Level-1] (sorted), and every orthonormal basis function e_p of
 *  P.  It is deflated and diagonalized: E are its energies.
 *
 *  The block is X_Level^S tensor the bath states Tuple.  A tuple is Level mode
 *  indices whose entries form S.  Element (t,lambda) of the block sits at
 *  Offset + t*Dim + lambda of the full vector and has energy
 *  E[lambda] + sum of the mode energies of tuple t.
 *
 *  Couplings: for every distinct entry i of S the space is reached from the
 *  parent space S\i by V[i][j].  Z[c][j] = < lambda | V[i][j] | eta >, Dim x
 *  Dim of the parent.  A link (t, tp, mu, f) says: tuple t of this block is
 *  tuple tp of the parent block with mode mu added, with amplitude
 *  f * sum_j L[mu][j] Z[c][j].                                                */
typedef struct BathExcSpaceType
{
  uint32_t Level;
  uint32_t *Entry;
  uint32_t NString;
  uint32_t *String;      /* NString * Level pair indices, position r is applied r-th */
  uint32_t Dim;
  double *E;
  CompactMatrixType U;   /* Dim x NString*North, < n | w_(s,p) >                     */
  CompactMatrixType C;   /* Dim x NString*North, | n > = sum_nu C_(n,nu) | w_nu >    */
  uint32_t NCoupling;
  uint32_t *CouplingEntry;
  uint32_t *CouplingParent;
  CompactMatrixType **Z; /* NCoupling arrays of Sigma[i].BlockSize matrices          */
  uint32_t NTuple;
  uint32_t *Tuple;       /* NTuple * Level                                           */
  double *TupleEnergy;   /* NTuple                                                   */
  uint32_t NLink;
  uint32_t *LinkTuple, *LinkParentTuple, *LinkMode, *LinkCoupling;
  double *LinkFactor;
  size_t Offset;
}
BathExcSpaceType;

/*  relative cut on the eigenvalues of the overlap matrix of a space X_m^S        */
static double const BathExcEpsDefl = 1e-10;
/*  relative cut when a residue of a hybridisation function is factorized       */
static double const BathExcEpsChol = 1e-10;
/*  relative cut in the re-orthogonalization of the block Lanczos               */
static double const BathExcEpsLanczos = 1e-12;
/*  memory (in Complex16) the stored Lanczos vectors may take before the full
 *  re-orthogonalization is given up for a local one                            */
static double const BathExcMaxLanczosMemory = 1.25e8;
/*  hard cap on the number of raw vectors of one space X_m^S                    */
static uint32_t const BathExcMaxRaw = 20000;
/*  hard cap on the number of strings enumerated for one space                  */
static double const BathExcMaxStringEnumeration = 1e7;

/*  The parity of the number of fermion ladder operators in every term of O.
 *  Returns 0 (even) or 1 (odd), -1 for an operator without terms and -2 when
 *  the terms do not agree.  Boson ladder operators carry an index of NFermion
 *  or above and do not count.
 */
int OperatorFermionParity(OperatorType const *const O)
{
  uint32_t il, n, k, nf;
  int parity = -1;
  for (il = 0; il <= O->MaxLength; il++)
  {
    for (n = 0; n < O->l[il].N; n++)
    {
      nf = 0;
      for (k = 0; k < il; k++)
      {
        if ((uint32_t)(O->l[il].o[n * il + k] & NotCrAnMask) < O->NFermion) nf++;
      }
      if (parity == -1) parity = (int)(nf % 2);
      else if (parity != (int)(nf % 2)) return -2;
    }
  }
  return parity;
}

/*  A complex Ni x Nj view on memory owned by someone else                       */
static CompactMatrixType BathExcView(Complex16 *const V, uint32_t const Ni, uint32_t const Nj)
{
  CompactMatrixType v;
  CompactMatrixSetUnallocated(&v);
  v.Ni = Ni;
  v.Nj = Nj;
  v.QComplex = 1;
  v.V = V;
  v.OnLuaStack = true;
  return v;
}

/*  lexicographic comparison of two tuples of length n                          */
static int BathExcTupleCompare(uint32_t const *const a, uint32_t const *const b, uint32_t const n)
{
  uint32_t r;
  for (r = 0; r < n; r++)
  {
    if (a[r] < b[r]) return -1;
    if (a[r] > b[r]) return 1;
  }
  return 0;
}

/*  index of tuple t (length n) in the lexicographically sorted list, or
 *  UINT32_MAX when it is not there                                             */
static uint32_t BathExcTupleFind(uint32_t const *const List, uint32_t const NList, uint32_t const *const t, uint32_t const n)
{
  uint32_t lo = 0, hi = NList, mid;
  int c;
  if (n == 0) return (NList > 0) ? 0 : UINT32_MAX;
  while (lo < hi)
  {
    mid = lo + (hi - lo) / 2;
    c = BathExcTupleCompare(&(List[(size_t)mid * n]), t, n);
    if (c == 0) return mid;
    if (c < 0) lo = mid + 1; else hi = mid;
  }
  return UINT32_MAX;
}

static void BathExcSpaceInit(BathExcSpaceType *const S)
{
  memset(S, 0, sizeof(BathExcSpaceType));
  CompactMatrixSetUnallocated(&(S->U));
  CompactMatrixSetUnallocated(&(S->C));
}

static void BathExcSpaceFree(BathExcSpaceType *const S, BlockListOfPoles const *const Sigma)
{
  uint32_t c, j;
  free(S->Entry);
  free(S->String);
  free(S->E);
  FreeCompactMatrix(&(S->U));
  FreeCompactMatrix(&(S->C));
  if (S->Z != NULL)
  {
    for (c = 0; c < S->NCoupling; c++)
    {
      if (S->Z[c] == NULL) continue;
      for (j = 0; j < Sigma[S->CouplingEntry[c]].BlockSize; j++) FreeCompactMatrix(&(S->Z[c][j]));
      free(S->Z[c]);
    }
    free(S->Z);
  }
  free(S->CouplingEntry);
  free(S->CouplingParent);
  free(S->Tuple);
  free(S->TupleEnergy);
  free(S->LinkTuple);
  free(S->LinkParentTuple);
  free(S->LinkMode);
  free(S->LinkCoupling);
  free(S->LinkFactor);
  BathExcSpaceInit(S);
}

/*  The bath modes.  Every residue W of every pole is factorized as L^dag L and
 *  every row of L is one mode, exactly as the old routine builds its bath.  With
 *  a grid the poles are first binned on it, which conserves the zeroth and the
 *  first moment of every bin and caps the number of modes of an entry at
 *  BlockSize times the number of grid points.  That matters here: the number of
 *  bath states grows as the number of modes to the power NExcitations.
 */
static int BathExcModesBuild(BlockListOfPoles const *const Sigma, uint32_t const NSigma, uint32_t const *const Parity, double const *const SuggestedGrid, uint32_t const NSuggestedGridPoints, BathExcModesType *const M, double *const DiscardedWeight, double *const NegativeWeight)
{
  uint32_t i, j, k, n, c, r, Ng, NModeMax = 0, BMax = 0, clow, chigh;
  double en, uleft, uright, wmax, norm;
  double *Grid = NULL, *wval = NULL;
  CompactMatrixType Wm, Wfun, *Wgrid = NULL;
  bool UseGrid = (SuggestedGrid != NULL) && (NSuggestedGridPoints >= 2);
  memset(M, 0, sizeof(BathExcModesType));
  for (i = 0; i < NSigma; i++) BMax = Max(BMax, Sigma[i].BlockSize);
  if (UseGrid)
  {
    double Emin = DBL_MAX, Emax = -DBL_MAX;
    for (k = 1; k < NSuggestedGridPoints; k++)
    {
      if (SuggestedGrid[k] <= SuggestedGrid[k - 1])
      {
        printf("CalculateGforHwithSigmaAsOperatorBathExcitations: the grid is not strictly increasing at position %u (%g <= %g)\n", k, SuggestedGrid[k], SuggestedGrid[k - 1]);
        return 1;
      }
    }
    for (i = 0; i < NSigma; i++)
    {
      for (n = 0; n < Sigma[i].N; n++)
      {
        en = Sigma[i].QComplex ? Sigma[i].V[n * (1 + Sqr(Sigma[i].BlockSize))].Re : Sigma[i].VRe[n * (1 + Sqr(Sigma[i].BlockSize))];
        Emin = Min(Emin, en);
        Emax = Max(Emax, en);
      }
    }
    Ng = NSuggestedGridPoints + ((Emin < SuggestedGrid[0]) ? 1 : 0) + ((Emax > SuggestedGrid[NSuggestedGridPoints - 1]) ? 1 : 0);
    Grid = (double*)malloc(Ng * sizeof(double));
    if (Grid == NULL) { printf("malloc failed in BathExcModesBuild\n"); return 1; }
    k = 0;
    if (Emin < SuggestedGrid[0]) Grid[k++] = Emin;
    for (c = 0; c < NSuggestedGridPoints; c++) Grid[k++] = SuggestedGrid[c];
    if (Emax > SuggestedGrid[NSuggestedGridPoints - 1]) Grid[k++] = Emax;
    for (i = 0; i < NSigma; i++) NModeMax += Ng * Sigma[i].BlockSize;
  }
  else
  {
    Ng = 0;
    for (i = 0; i < NSigma; i++) NModeMax += Sigma[i].N * Sigma[i].BlockSize;
  }
  M->BMax = BMax;
  M->First = (uint32_t*)calloc(NSigma + 1, sizeof(uint32_t));
  M->Parity = (uint32_t*)calloc(NSigma + 1, sizeof(uint32_t));
  M->Entry = (uint32_t*)calloc(Max(NModeMax, 1), sizeof(uint32_t));
  M->E = (double*)calloc(Max(NModeMax, 1), sizeof(double));
  M->L = (Complex16*)calloc((size_t)Max(NModeMax, 1) * BMax, sizeof(Complex16));
  wval = (double*)malloc(Max(BMax, 1) * sizeof(double));
  if ((M->First == NULL) || (M->Parity == NULL) || (M->Entry == NULL) || (M->E == NULL) || (M->L == NULL) || (wval == NULL))
  {
    printf("malloc failed in BathExcModesBuild\n");
    free(Grid); free(wval);
    return 1;
  }
  for (i = 0; i < NSigma; i++) M->Parity[i] = Parity[i];
  M->N = 0;
  for (i = 0; i < NSigma; i++)
  {
    uint32_t bi = Sigma[i].BlockSize, NPole;
    M->First[i] = M->N;
    CompactMatrixSetUnallocated(&Wm);
    Wm.Ni = bi; Wm.Nj = bi; Wm.QComplex = 1;
    if (InitCompactMatrix(&Wm)) { printf("InitCompactMatrix failed in BathExcModesBuild\n"); free(Grid); free(wval); return 1; }
    if (UseGrid)
    {
      Wgrid = (CompactMatrixType*)calloc(Ng, sizeof(CompactMatrixType));
      if (Wgrid == NULL) { printf("malloc failed in BathExcModesBuild\n"); FreeCompactMatrix(&Wm); free(Grid); free(wval); return 1; }
      for (c = 0; c < Ng; c++)
      {
        CompactMatrixSetUnallocated(&(Wgrid[c]));
        Wgrid[c].Ni = bi; Wgrid[c].Nj = bi; Wgrid[c].QComplex = 1;
        if (InitCompactMatrix(&(Wgrid[c]))) { printf("InitCompactMatrix failed in BathExcModesBuild\n"); FreeCompactMatrix(&Wm); free(Grid); free(wval); return 1; }
      }
      for (n = 0; n < Sigma[i].N; n++)
      {
        BlockListOfPolesGetPole(&(Sigma[i]), n, &Wm, &en);
        clow = 0;
        chigh = Ng - 1;
        while (chigh - clow > 1)
        {
          c = (clow + chigh) / 2;
          if (Grid[c] <= en) clow = c; else chigh = c;
        }
        FindWeigth(&(Grid[clow]), &en, &(Grid[chigh]), &uleft, &uright);
        for (k = 0; k < bi * bi; k++)
        {
          Wgrid[clow].V[k].Re  += uleft  * Wm.V[k].Re;  Wgrid[clow].V[k].Im  += uleft  * Wm.V[k].Im;
          Wgrid[chigh].V[k].Re += uright * Wm.V[k].Re;  Wgrid[chigh].V[k].Im += uright * Wm.V[k].Im;
        }
      }
      NPole = Ng;
    }
    else
    {
      NPole = Sigma[i].N;
    }
    for (n = 0; n < NPole; n++)
    {
      if (UseGrid)
      {
        en = Grid[n];
        memcpy(Wm.V, Wgrid[n].V, bi * bi * sizeof(Complex16));
        SquareCompactMatrixMakeHermitian(&Wm);
      }
      else
      {
        BlockListOfPolesGetPole(&(Sigma[i]), n, &Wm, &en);
      }
      if (HermitianCompactMatrixEigensystemRowVectors(&Wm, &Wfun, wval)) { FreeCompactMatrix(&Wm); free(Grid); free(wval); return 1; }
      wmax = Max(wval[bi - 1], 0.0);
      for (k = 0; k < bi; k++)
      {
        if (wval[k] > BathExcEpsChol * wmax)
        {
          norm = sqrt(wval[k]);
          r = M->N;
          M->Entry[r] = i;
          M->E[r] = en;
          for (j = 0; j < bi; j++)
          {
            M->L[(size_t)r * BMax + j].Re =  norm * Wfun.V[k * bi + j].Re;
            M->L[(size_t)r * BMax + j].Im = -norm * Wfun.V[k * bi + j].Im;
          }
          M->N++;
        }
        else if (wval[k] < -BathExcEpsChol * Max(wmax, DBL_MIN)) *NegativeWeight -= wval[k];
        else if (wval[k] > 0) *DiscardedWeight += wval[k];
      }
      FreeCompactMatrix(&Wfun);
    }
    if (UseGrid)
    {
      for (c = 0; c < Ng; c++) FreeCompactMatrix(&(Wgrid[c]));
      free(Wgrid);
      Wgrid = NULL;
    }
    FreeCompactMatrix(&Wm);
  }
  M->First[NSigma] = M->N;
  free(Grid);
  free(wval);
  return 0;
}

static void BathExcModesFree(BathExcModesType *const M)
{
  free(M->Entry);
  free(M->First);
  free(M->Parity);
  free(M->E);
  free(M->L);
  memset(M, 0, sizeof(BathExcModesType));
}

/*  Everything the matrix vector product needs                                   */
typedef struct BathExcSystemType
{
  uint32_t North;
  CompactMatrixType *HO;          /* North x North, the Hamiltonian in P        */
  uint32_t NSpace;
  BathExcSpaceType *Space;        /* Space[0] is P                              */
  BathExcModesType *Modes;
  BlockListOfPoles const *Sigma;
  size_t Dtot;
  Complex16 *WorkU;               /* scratch, WorkSize each                     */
  Complex16 *WorkW;
  size_t WorkSize;
}
BathExcSystemType;

/*  y = M x for the auxiliary Hamiltonian
 *
 *      M = (  HO    A_1^dag                )
 *          (  A_1   D_1      A_2^dag       )
 *          (        A_2      D_2      ...  )
 *
 *  with D_m diagonal and A_m the couplings of the links.  M is never formed.
 */
static int BathExcMatVec(BathExcSystemType const *const Sys, Complex16 const *const x, Complex16 *const y)
{
  uint32_t s, t, l, c, j, lambda, bi;
  size_t k;
  BathExcSpaceType const *S, *Sp;
  CompactMatrixType xv, yv, Uj, Wj;
  Complex16 coef, cc, *ys, *xsrow, *Ujrow, *Wjrow;
  double e;
  memset(y, 0, Sys->Dtot * sizeof(Complex16));
  /*  P  */
  xv = BathExcView((Complex16*)x, 1, Sys->North);
  yv = BathExcView(y, 1, Sys->North);
  CompactMatrixAddMultiply('N', 'T', &xv, Sys->HO, &yv, 1, 0, 0, 0);
  /*  the diagonal of the bath blocks  */
  for (s = 1; s < Sys->NSpace; s++)
  {
    S = &(Sys->Space[s]);
    for (t = 0; t < S->NTuple; t++)
    {
      for (lambda = 0; lambda < S->Dim; lambda++)
      {
        k = S->Offset + (size_t)t * S->Dim + lambda;
        e = S->E[lambda] + S->TupleEnergy[t];
        y[k].Re = e * x[k].Re;
        y[k].Im = e * x[k].Im;
      }
    }
  }
  /*  the couplings, space by space and parent by parent  */
  for (s = 1; s < Sys->NSpace; s++)
  {
    S = &(Sys->Space[s]);
    for (c = 0; c < S->NCoupling; c++)
    {
      Sp = &(Sys->Space[S->CouplingParent[c]]);
      bi = Sys->Sigma[S->CouplingEntry[c]].BlockSize;
      /*  U_j = X_parent Z_j^T  and  W_j = 0  */
      xv = BathExcView((Complex16*)&(x[Sp->Offset]), Sp->NTuple, Sp->Dim);
      for (j = 0; j < bi; j++)
      {
        Uj = BathExcView(&(Sys->WorkU[(size_t)j * Sp->NTuple * S->Dim]), Sp->NTuple, S->Dim);
        CompactMatrixAddMultiply('N', 'T', &xv, &(S->Z[c][j]), &Uj, 1, 0, 0, 0);
        memset(&(Sys->WorkW[(size_t)j * Sp->NTuple * S->Dim]), 0, (size_t)Sp->NTuple * S->Dim * sizeof(Complex16));
      }
      ys = &(y[S->Offset]);
      for (l = 0; l < S->NLink; l++)
      {
        if (S->LinkCoupling[l] != c) continue;
        t = S->LinkTuple[l];
        xsrow = (Complex16*)&(x[S->Offset + (size_t)t * S->Dim]);
        for (j = 0; j < bi; j++)
        {
          cc = Sys->Modes->L[(size_t)S->LinkMode[l] * Sys->Modes->BMax + j];
          coef.Re = S->LinkFactor[l] * cc.Re;
          coef.Im = S->LinkFactor[l] * cc.Im;
          if ((coef.Re == 0) && (coef.Im == 0)) continue;
          Ujrow = &(Sys->WorkU[((size_t)j * Sp->NTuple + S->LinkParentTuple[l]) * S->Dim]);
          Wjrow = &(Sys->WorkW[((size_t)j * Sp->NTuple + S->LinkParentTuple[l]) * S->Dim]);
          for (lambda = 0; lambda < S->Dim; lambda++)
          {
            /*  y_S[t] += coef U_j[tp]  */
            ys[(size_t)t * S->Dim + lambda].Re += coef.Re * Ujrow[lambda].Re - coef.Im * Ujrow[lambda].Im;
            ys[(size_t)t * S->Dim + lambda].Im += coef.Re * Ujrow[lambda].Im + coef.Im * Ujrow[lambda].Re;
            /*  W_j[tp] += conj(coef) x_S[t]  */
            Wjrow[lambda].Re += coef.Re * xsrow[lambda].Re + coef.Im * xsrow[lambda].Im;
            Wjrow[lambda].Im += coef.Re * xsrow[lambda].Im - coef.Im * xsrow[lambda].Re;
          }
        }
      }
      /*  y_parent += W_j conj(Z_j)  */
      yv = BathExcView(&(y[Sp->Offset]), Sp->NTuple, Sp->Dim);
      for (j = 0; j < bi; j++)
      {
        Wj = BathExcView(&(Sys->WorkW[(size_t)j * Sp->NTuple * S->Dim]), Sp->NTuple, S->Dim);
        CompactMatrixAddMultiply('N', '*', &Wj, &(S->Z[c][j]), &yv, 1, 0, 1, 0);
      }
    }
  }
  return 0;
}

/*  What the front end hands to the back ends: P orthonormalized (HO, Cseed),
 *  the system spaces X_m^S with their energies and couplings, and the parity
 *  of every hybridisation.  Everything in it is owned by it and released by
 *  BathExcFrontEndFree.                                                        */
typedef struct BathExcFrontType
{
  uint32_t North, NFull, QComplexOut, NSpace, NSpaceAlloc;
  uint32_t *Parity;
  CompactMatrixType HO, Cseed;
  BathExcSpaceType *Space;
  OperatorType Heff;
  bool QHeffAllocated;
}
BathExcFrontType;

static void BathExcFrontEndInit(BathExcFrontType *const F)
{
  memset(F, 0, sizeof(BathExcFrontType));
  CompactMatrixSetUnallocated(&(F->HO));
  CompactMatrixSetUnallocated(&(F->Cseed));
}

static void BathExcFrontEndFree(BathExcFrontType *const F, BlockListOfPoles const *const Sigma)
{
  uint32_t s;
  FreeCompactMatrix(&(F->HO));
  FreeCompactMatrix(&(F->Cseed));
  if (F->Space != NULL)
  {
    for (s = 0; s < F->NSpaceAlloc; s++) BathExcSpaceFree(&(F->Space[s]), Sigma);
    free(F->Space);
  }
  free(F->Parity);
  if (F->QHeffAllocated) FreeOperator(&(F->Heff));
  BathExcFrontEndInit(F);
}

/*  Steps A, B and C of CalculateGforHwithSigmaAsOperatorBathExcitations, shared
 *  with the energy mesh back end CalculateGforHwithSigmaAsOperatorBathExcitationsMesh.
 *  On failure everything is released and 1 is returned.                        */
static int BathExcFrontEndBuild(const OperatorType *const H, const WaveFunctionType *const PsiBasis, const uint32_t NPsiBasis, const uint32_t BlockSize, OperatorType ** V, const BlockListOfPoles *const Sigma, const uint32_t NSigma, const uint32_t NExcitations, const bool NonCrossing, const uint32_t NTriFull, const CompactMatrixType *const B0, BathExcFrontType *const F)
{
  uint32_t i, j, jp, k, p, q, a, s, sp, r, lev, c;
  uint32_t Npsi = NPsiBasis, m = BlockSize, NFull, North = 0;
  uint32_t QComplexOut, NPair = 0, NSpaceAlloc = 0, NSpace = 0;
  uint32_t *PairEntry = NULL, *PairChannel = NULL, *Parity = NULL;
  int ierr = 0, par;
  double *Sval = NULL, norm;
  bool QHeffAllocated = false;
  CompactMatrixType SMat, HMat, Sfun, T, TS, HO, Cseed;
  OperatorType Heff, *Hptr;
  BathExcSpaceType *Space = NULL;

  BathExcFrontEndInit(F);
  CompactMatrixSetUnallocated(&SMat);  CompactMatrixSetUnallocated(&HMat);
  CompactMatrixSetUnallocated(&Sfun);  CompactMatrixSetUnallocated(&T);
  CompactMatrixSetUnallocated(&TS);    CompactMatrixSetUnallocated(&HO);
  CompactMatrixSetUnallocated(&Cseed);

  /* ---------------------------------------------------------------------- */
  /*  A   input                                                             */
  /* ---------------------------------------------------------------------- */
  if ((H == NULL) || (PsiBasis == NULL))
  {
    printf("CalculateGforHwithSigmaAsOperatorBathExcitations called with a NULL pointer for H or PsiBasis\n");
    return 1;
  }
  if ((Npsi == 0) || (m == 0) || (m > Npsi))
  {
    printf("CalculateGforHwithSigmaAsOperatorBathExcitations called with NPsiBasis=%u and BlockSize=%u, needs 0 < BlockSize <= NPsiBasis\n", Npsi, m);
    return 1;
  }
  if ((NTriFull == 0) || (NExcitations == 0))
  {
    printf("CalculateGforHwithSigmaAsOperatorBathExcitations called with NTriFull=%u and NExcitations=%u, both need to be at least 1\n", NTriFull, NExcitations);
    return 1;
  }
  if (B0 == NULL)
  {
    NFull = m;
  }
  else
  {
    NFull = B0->Ni;
    if ((B0->Nj != m) || (NFull < m))
    {
      printf("CalculateGforHwithSigmaAsOperatorBathExcitations called with B0 of size %u x %u, needs BlockSize=%u columns and at least as many rows\n", B0->Ni, B0->Nj, m);
      return 1;
    }
  }
  Parity = (uint32_t*)calloc(Max(NSigma, 1), sizeof(uint32_t));
  if (Parity == NULL) { printf("malloc failed in CalculateGforHwithSigmaAsOperatorBathExcitations\n"); return 1; }
  for (i = 0; i < NSigma; i++)
  {
    int pi = -1;
    if (Sigma[i].BlockSize == 0)
    {
      printf("CalculateGforHwithSigmaAsOperatorBathExcitations called with Sigma[%u].BlockSize=0\n", i);
      free(Parity); return 1;
    }
    if ((V == NULL) || (V[i] == NULL))
    {
      printf("CalculateGforHwithSigmaAsOperatorBathExcitations called with no operators for Sigma[%u]\n", i);
      free(Parity); return 1;
    }
    /*  the statistics of the bath particle is the parity of its operators     */
    for (j = 0; j < Sigma[i].BlockSize; j++)
    {
      par = OperatorFermionParity(&(V[i][j]));
      if (par == -2)
      {
        printf("CalculateGforHwithSigmaAsOperatorBathExcitations: operator %u of hybridisation %u mixes terms with an even and an odd number of fermion operators.\n", j + 1, i + 1);
        printf("  Such an operator can not couple to a bath of free particles, the coupling would not conserve fermion parity.\n");
        free(Parity); return 1;
      }
      if (par == -1) continue;
      if ((pi != -1) && (pi != par))
      {
        printf("CalculateGforHwithSigmaAsOperatorBathExcitations: the operators of hybridisation %u do not all have the same parity.\n", i + 1);
        printf("  One hybridisation function stands for one kind of particle, fermion (odd operators) or boson (even operators).\n");
        free(Parity); return 1;
      }
      pi = par;
    }
    Parity[i] = (pi == 1) ? 1 : 0;
  }
  QComplexOut = (H->QComplex != 0);
  for (a = 0; a < Npsi; a++) if (PsiBasis[a].QComplex) QComplexOut = 1;
  for (i = 0; i < NSigma; i++)
  {
    if (Sigma[i].QComplex) QComplexOut = 1;
    for (j = 0; j < Sigma[i].BlockSize; j++) if (V[i][j].QComplex) QComplexOut = 1;
  }
  if ((B0 != NULL) && (B0->QComplex)) QComplexOut = 1;
  VerbosityPrintf(0x0001, "CalculateGforHwithSigmaAsOperatorBathExcitations: NPsiBasis=%u BlockSize=%u NSigma=%u NExcitations=%u NonCrossing=%s NTriFull=%u\n", Npsi, m, NSigma, NExcitations, NonCrossing ? "true" : "false", NTriFull);

  /* ---------------------------------------------------------------------- */
  /*  A2  the constant term of every hybridisation function goes into H      */
  /* ---------------------------------------------------------------------- */
  Hptr = (OperatorType*)H;
  for (i = 0; i < NSigma; i++)
  {
    uint32_t bi = Sigma[i].BlockSize;
    CompactMatrixType A0;
    OperatorType Oprod;
    CompactMatrixSetUnallocated(&A0);
    A0.Ni = bi; A0.Nj = bi; A0.QComplex = 1;
    if (InitCompactMatrix(&A0)) { printf("InitCompactMatrix failed in CalculateGforHwithSigmaAsOperatorBathExcitations\n"); ierr = 1; goto cleanup; }
    BlockListOfPolesGetA0(&(Sigma[i]), &A0);
    for (j = 0; j < bi; j++)
    {
      for (jp = 0; jp < bi; jp++)
      {
        if ((fabs(A0.V[j * bi + jp].Re) < 100 * DBL_MIN) && (fabs(A0.V[j * bi + jp].Im) < 100 * DBL_MIN)) continue;
        if (QHeffAllocated == false)
        {
          if (OperatorCopy(&Heff, *H)) { printf("OperatorCopy failed in CalculateGforHwithSigmaAsOperatorBathExcitations\n"); FreeCompactMatrix(&A0); ierr = 1; goto cleanup; }
          QHeffAllocated = true;
          Hptr = &Heff;
        }
        InitOperatorForMult(&Oprod, 2, &(V[i][j]), &(V[i][jp]));
        OperatorMultiply2NormalOrder(&Oprod, V[i][j], V[i][jp], 'D', 'N', &QDetNotOrderedFcFaBcBa);
        OperatorAddNormalOrder(&Heff, Oprod, A0.V[j * bi + jp].Re, A0.V[j * bi + jp].Im, &QDetNotOrderedFcFaBcBa);
        FreeOperator(&Oprod);
      }
    }
    FreeCompactMatrix(&A0);
  }

  /* ---------------------------------------------------------------------- */
  /*  B   overlap and Hamiltonian in the given basis, then orthonormalize    */
  /*      (the same as step B of CalculateGforHwithSigmaAsOperator)          */
  /* ---------------------------------------------------------------------- */
  SMat.Ni = Npsi; SMat.Nj = Npsi; SMat.QComplex = 1;
  HMat.Ni = Npsi; HMat.Nj = Npsi; HMat.QComplex = 1;
  if (InitCompactMatrix(&SMat) || InitCompactMatrix(&HMat)) { printf("InitCompactMatrix failed in CalculateGforHwithSigmaAsOperatorBathExcitations\n"); ierr = 1; goto cleanup; }
  if (PsiListPsiListCompactMatrix(PsiBasis, Npsi, &SMat)) { ierr = 1; goto cleanup; }
  if (PsiListOperatorPsiListCompactMatrix(Hptr, PsiBasis, Npsi, &HMat)) { ierr = 1; goto cleanup; }
  SquareCompactMatrixMakeHermitian(&HMat);
  Sval = (double*)malloc(Npsi * sizeof(double));
  if (Sval == NULL) { printf("malloc failed in CalculateGforHwithSigmaAsOperatorBathExcitations\n"); ierr = 1; goto cleanup; }
  if (HermitianCompactMatrixEigensystemRowVectors(&SMat, &Sfun, Sval)) { ierr = 1; goto cleanup; }
  North = 0;
  for (p = 0; p < Npsi; p++) if (Sval[p] > SigmaAsOperatorEpsOrtho * Sval[Npsi - 1]) North++;
  if (North < m)
  {
    printf("CalculateGforHwithSigmaAsOperatorBathExcitations: the %u basis functions span only %u directions, fewer than the %u seeds\n", Npsi, North, m);
    ierr = 1; goto cleanup;
  }
  T.Ni = Npsi; T.Nj = North; T.QComplex = 1;
  if (InitCompactMatrix(&T)) { printf("InitCompactMatrix failed in CalculateGforHwithSigmaAsOperatorBathExcitations\n"); ierr = 1; goto cleanup; }
  q = 0;
  for (p = 0; p < Npsi; p++)
  {
    if (Sval[p] <= SigmaAsOperatorEpsOrtho * Sval[Npsi - 1]) continue;
    norm = 1.0 / sqrt(Sval[p]);
    for (a = 0; a < Npsi; a++)
    {
      T.V[a * North + q].Re = norm * Sfun.V[p * Npsi + a].Re;
      T.V[a * North + q].Im = norm * Sfun.V[p * Npsi + a].Im;
    }
    q++;
  }
  FreeCompactMatrix(&Sfun); CompactMatrixSetUnallocated(&Sfun);
  TS.Ni = North; TS.Nj = Npsi; TS.QComplex = 1;
  HO.Ni = North; HO.Nj = North; HO.QComplex = 1;
  if (InitCompactMatrix(&TS) || InitCompactMatrix(&HO)) { printf("InitCompactMatrix failed in CalculateGforHwithSigmaAsOperatorBathExcitations\n"); ierr = 1; goto cleanup; }
  CompactMatrixAddMultiply('C', 'N', &T, &HMat, &TS, 1, 0, 0, 0);
  CompactMatrixAddMultiply('N', 'N', &TS, &T, &HO, 1, 0, 0, 0);
  SquareCompactMatrixMakeHermitian(&HO);
  CompactMatrixAddMultiply('C', 'N', &T, &SMat, &TS, 1, 0, 0, 0);
  Cseed.Ni = North; Cseed.Nj = m; Cseed.QComplex = 1;
  if (InitCompactMatrix(&Cseed)) { printf("InitCompactMatrix failed in CalculateGforHwithSigmaAsOperatorBathExcitations\n"); ierr = 1; goto cleanup; }
  for (p = 0; p < North; p++)
  {
    for (a = 0; a < m; a++) Cseed.V[p * m + a] = TS.V[p * Npsi + a];
  }
  FreeCompactMatrix(&HMat); CompactMatrixSetUnallocated(&HMat);
  FreeCompactMatrix(&SMat); CompactMatrixSetUnallocated(&SMat);
  FreeCompactMatrix(&TS); CompactMatrixSetUnallocated(&TS);

  /* ---------------------------------------------------------------------- */
  /*  C   the system spaces X_m^S                                            */
  /* ---------------------------------------------------------------------- */
  for (i = 0; i < NSigma; i++) NPair += Sigma[i].BlockSize;
  PairEntry = (uint32_t*)malloc(Max(NPair, 1) * sizeof(uint32_t));
  PairChannel = (uint32_t*)malloc(Max(NPair, 1) * sizeof(uint32_t));
  if ((PairEntry == NULL) || (PairChannel == NULL)) { printf("malloc failed in CalculateGforHwithSigmaAsOperatorBathExcitations\n"); ierr = 1; goto cleanup; }
  k = 0;
  for (i = 0; i < NSigma; i++)
  {
    for (j = 0; j < Sigma[i].BlockSize; j++) { PairEntry[k] = i; PairChannel[k] = j; k++; }
  }
  /*  the number of multisets of size <= NExcitations out of NSigma entries    */
  {
    double count = 1, total = 1;
    for (lev = 1; lev <= NExcitations; lev++)
    {
      count = count * (double)(NSigma + lev - 1) / (double)lev;
      total += count;
    }
    if (total > 1e6)
    {
      printf("CalculateGforHwithSigmaAsOperatorBathExcitations: %g system spaces for %u hybridisations and %u excitations, that is not a sensible calculation\n", total, NSigma, NExcitations);
      ierr = 1; goto cleanup;
    }
    NSpaceAlloc = (uint32_t)total;
  }
  Space = (BathExcSpaceType*)malloc(NSpaceAlloc * sizeof(BathExcSpaceType));
  if (Space == NULL) { printf("malloc failed in CalculateGforHwithSigmaAsOperatorBathExcitations\n"); ierr = 1; goto cleanup; }
  for (s = 0; s < NSpaceAlloc; s++) BathExcSpaceInit(&(Space[s]));
  /*  Space[0] is P: one empty tuple, the raw vectors are e_p themselves       */
  Space[0].Level = 0;
  Space[0].Dim = North;
  Space[0].NString = 1;
  Space[0].NTuple = 1;
  Space[0].TupleEnergy = (double*)calloc(1, sizeof(double));
  if (Space[0].TupleEnergy == NULL) { printf("malloc failed in CalculateGforHwithSigmaAsOperatorBathExcitations\n"); ierr = 1; goto cleanup; }
  NSpace = 1;

  for (lev = 1; lev <= NExcitations; lev++)
  {
    uint32_t *ms = (uint32_t*)calloc(lev, sizeof(uint32_t));
    uint32_t *seq = (uint32_t*)calloc(lev, sizeof(uint32_t));
    uint32_t *sorted = (uint32_t*)calloc(lev, sizeof(uint32_t));
    uint32_t NSpaceBefore = NSpace;
    bool more;
    if ((ms == NULL) || (seq == NULL) || (sorted == NULL)) { free(ms); free(seq); free(sorted); printf("malloc failed in CalculateGforHwithSigmaAsOperatorBathExcitations\n"); ierr = 1; goto cleanup; }
    if (pow((double)NPair, (double)lev) > BathExcMaxStringEnumeration)
    {
      printf("CalculateGforHwithSigmaAsOperatorBathExcitations: %u channels to the power %u is too many strings to enumerate\n", NPair, lev);
      free(ms); free(seq); free(sorted); ierr = 1; goto cleanup;
    }
    /*  loop over the non decreasing sequences ms of entries                   */
    more = (NSigma > 0);
    while (more)
    {
      BathExcSpaceType *S = &(Space[NSpace]);
      uint32_t NRaw, NStringMax, d, di, pc;
      OperatorType *VsOp = NULL;
      CompactMatrixType SextO, HextO, Mtmp, MtmpO, TM, Xdefl, Ws, Hred, Sxfun, ufun;
      double *sxval = NULL;
      OperatorType Oprod;
      bool ok = true;

      BathExcSpaceInit(S);
      S->Level = lev;
      S->Entry = (uint32_t*)malloc(lev * sizeof(uint32_t));
      NStringMax = (uint32_t)pow((double)NPair, (double)lev);
      S->String = (uint32_t*)malloc((size_t)NStringMax * lev * sizeof(uint32_t));
      if ((S->Entry == NULL) || (S->String == NULL)) { free(ms); free(seq); free(sorted); printf("malloc failed in CalculateGforHwithSigmaAsOperatorBathExcitations\n"); ierr = 1; goto cleanup; }
      memcpy(S->Entry, ms, lev * sizeof(uint32_t));
      /*  all strings of channel pairs whose entries are the multiset ms        */
      S->NString = 0;
      for (r = 0; r < lev; r++) seq[r] = 0;
      while (true)
      {
        /*  sort the entries of seq and compare with ms                        */
        for (r = 0; r < lev; r++) sorted[r] = PairEntry[seq[r]];
        for (r = 1; r < lev; r++)
        {
          uint32_t v = sorted[r], w = r;
          while ((w > 0) && (sorted[w - 1] > v)) { sorted[w] = sorted[w - 1]; w--; }
          sorted[w] = v;
        }
        if (memcmp(sorted, ms, lev * sizeof(uint32_t)) == 0)
        {
          memcpy(&(S->String[(size_t)S->NString * lev]), seq, lev * sizeof(uint32_t));
          S->NString++;
        }
        r = lev;
        while (r > 0)
        {
          r--;
          seq[r]++;
          if (seq[r] < NPair) break;
          seq[r] = 0;
          if (r == 0) { r = UINT32_MAX; break; }
        }
        if (r == UINT32_MAX) break;
      }
      NRaw = S->NString * North;
      if (NRaw > BathExcMaxRaw)
      {
        printf("CalculateGforHwithSigmaAsOperatorBathExcitations: space %u of level %u has %u strings times %u basis functions = %u raw vectors,\n", NSpace, lev, S->NString, North, NRaw);
        printf("  which is above the compiled in limit BathExcMaxRaw = %u. Reduce NTri or the number of channels.\n", BathExcMaxRaw);
        free(ms); free(seq); free(sorted); ierr = 1; goto cleanup;
      }
      /*  a parent has to exist for every distinct entry, otherwise the space
       *  could still be reached, but its couplings to that parent are zero      */
      /*  the string operators V[s_lev] ... V[s_1]                              */
      VsOp = (OperatorType*)calloc(S->NString, sizeof(OperatorType));
      if (VsOp == NULL) { free(ms); free(seq); free(sorted); printf("malloc failed in CalculateGforHwithSigmaAsOperatorBathExcitations\n"); ierr = 1; goto cleanup; }
      for (sp = 0; sp < S->NString; sp++)
      {
        uint32_t *str = &(S->String[(size_t)sp * lev]);
        OperatorCopy(&(VsOp[sp]), V[PairEntry[str[0]]][PairChannel[str[0]]]);
        for (r = 1; r < lev; r++)
        {
          OperatorType *Vr = &(V[PairEntry[str[r]]][PairChannel[str[r]]]);
          InitOperatorForMult(&Oprod, 2, Vr, &(VsOp[sp]));
          OperatorMultiply2NormalOrder(&Oprod, *Vr, VsOp[sp], 'N', 'N', &QDetNotOrderedFcFaBcBa);
          FreeOperator(&(VsOp[sp]));
          VsOp[sp] = Oprod;
        }
      }
      /*  Sext and Hext of the raw vectors, from operators that act in P       */
      CompactMatrixSetUnallocated(&SextO); CompactMatrixSetUnallocated(&HextO);
      CompactMatrixSetUnallocated(&Mtmp);  CompactMatrixSetUnallocated(&MtmpO);
      CompactMatrixSetUnallocated(&TM);
      SextO.Ni = NRaw;  SextO.Nj = NRaw;  SextO.QComplex = 1;
      HextO.Ni = NRaw;  HextO.Nj = NRaw;  HextO.QComplex = 1;
      Mtmp.Ni  = Npsi;  Mtmp.Nj  = Npsi;  Mtmp.QComplex  = 1;
      MtmpO.Ni = North; MtmpO.Nj = North; MtmpO.QComplex = 1;
      TM.Ni    = North; TM.Nj    = Npsi;  TM.QComplex    = 1;
      if (InitCompactMatrix(&SextO) || InitCompactMatrix(&HextO) || InitCompactMatrix(&Mtmp) || InitCompactMatrix(&MtmpO) || InitCompactMatrix(&TM))
      {
        printf("InitCompactMatrix failed in CalculateGforHwithSigmaAsOperatorBathExcitations for space %u\n", NSpace);
        free(ms); free(seq); free(sorted); ierr = 1; goto cleanup;
      }
      for (sp = 0; sp < S->NString; sp++)
      {
        for (pc = sp; pc < S->NString; pc++)
        {
          uint32_t pass;
          for (pass = 0; pass < 2; pass++)
          {
            CompactMatrixType *Dest = (pass == 0) ? &SextO : &HextO;
            if (pass == 0)
            {
              InitOperatorForMult(&Oprod, 2, &(VsOp[sp]), &(VsOp[pc]));
              OperatorMultiply2NormalOrder(&Oprod, VsOp[sp], VsOp[pc], 'D', 'N', &QDetNotOrderedFcFaBcBa);
            }
            else
            {
              InitOperatorForMult(&Oprod, 3, &(VsOp[sp]), Hptr, &(VsOp[pc]));
              OperatorMultiply3NormalOrder(&Oprod, VsOp[sp], *Hptr, VsOp[pc], 'D', 'N', 'N', &QDetNotOrderedFcFaBcBa);
            }
            if (PsiListOperatorPsiListCompactMatrix(&Oprod, PsiBasis, Npsi, &Mtmp)) { FreeOperator(&Oprod); ok = false; break; }
            FreeOperator(&Oprod);
            CompactMatrixAddMultiply('C', 'N', &T, &Mtmp, &TM, 1, 0, 0, 0);
            CompactMatrixAddMultiply('N', 'N', &TM, &T, &MtmpO, 1, 0, 0, 0);
            for (p = 0; p < North; p++)
            {
              for (q = 0; q < North; q++)
              {
                Dest->V[((size_t)sp * North + p) * NRaw + pc * North + q] = MtmpO.V[p * North + q];
                Dest->V[((size_t)pc * North + q) * NRaw + sp * North + p].Re =  MtmpO.V[p * North + q].Re;
                Dest->V[((size_t)pc * North + q) * NRaw + sp * North + p].Im = -MtmpO.V[p * North + q].Im;
              }
            }
          }
          if (!ok) break;
        }
        if (!ok) break;
      }
      for (sp = 0; sp < S->NString; sp++) FreeOperator(&(VsOp[sp]));
      free(VsOp);
      FreeCompactMatrix(&Mtmp);
      FreeCompactMatrix(&MtmpO);
      FreeCompactMatrix(&TM);
      if (!ok) { FreeCompactMatrix(&SextO); FreeCompactMatrix(&HextO); free(ms); free(seq); free(sorted); ierr = 1; goto cleanup; }
      /*  deflate, diagonalize (the same as step C of the old routine)          */
      sxval = (double*)malloc(NRaw * sizeof(double));
      if (sxval == NULL) { printf("malloc failed in CalculateGforHwithSigmaAsOperatorBathExcitations\n"); free(ms); free(seq); free(sorted); ierr = 1; goto cleanup; }
      if (HermitianCompactMatrixEigensystemRowVectors(&SextO, &Sxfun, sxval)) { free(sxval); free(ms); free(seq); free(sorted); ierr = 1; goto cleanup; }
      FreeCompactMatrix(&SextO);
      di = 0;
      for (k = 0; k < NRaw; k++) if (sxval[k] > BathExcEpsDefl * sxval[NRaw - 1]) di++;
      S->Dim = di;
      if (di == 0)
      {
        VerbosityPrintf(0x0001, "CalculateGforHwithSigmaAsOperatorBathExcitations: space %u of level %u is empty, skipped\n", NSpace, lev);
        free(sxval);
        FreeCompactMatrix(&Sxfun);
        FreeCompactMatrix(&HextO);
        BathExcSpaceFree(S, Sigma);
      }
      else
      {
        CompactMatrixSetUnallocated(&Xdefl); CompactMatrixSetUnallocated(&Ws);
        CompactMatrixSetUnallocated(&Hred);  CompactMatrixSetUnallocated(&TM);
        Xdefl.Ni = NRaw; Xdefl.Nj = di; Xdefl.QComplex = 1;
        Ws.Ni = di; Ws.Nj = NRaw; Ws.QComplex = 1;
        Hred.Ni = di; Hred.Nj = di; Hred.QComplex = 1;
        TM.Ni = di; TM.Nj = NRaw; TM.QComplex = 1;
        if (InitCompactMatrix(&Xdefl) || InitCompactMatrix(&Ws) || InitCompactMatrix(&Hred) || InitCompactMatrix(&TM))
        {
          printf("InitCompactMatrix failed in CalculateGforHwithSigmaAsOperatorBathExcitations for space %u\n", NSpace);
          free(sxval); free(ms); free(seq); free(sorted); ierr = 1; goto cleanup;
        }
        q = 0;
        for (k = 0; k < NRaw; k++)
        {
          if (sxval[k] <= BathExcEpsDefl * sxval[NRaw - 1]) continue;
          norm = 1.0 / sqrt(sxval[k]);
          for (p = 0; p < NRaw; p++)
          {
            Xdefl.V[(size_t)p * di + q].Re = norm * Sxfun.V[(size_t)k * NRaw + p].Re;
            Xdefl.V[(size_t)p * di + q].Im = norm * Sxfun.V[(size_t)k * NRaw + p].Im;
            Ws.V[(size_t)q * NRaw + p].Re =  sxval[k] * norm * Sxfun.V[(size_t)k * NRaw + p].Re;
            Ws.V[(size_t)q * NRaw + p].Im = -sxval[k] * norm * Sxfun.V[(size_t)k * NRaw + p].Im;
          }
          q++;
        }
        FreeCompactMatrix(&Sxfun);
        free(sxval);
        CompactMatrixAddMultiply('C', 'N', &Xdefl, &HextO, &TM, 1, 0, 0, 0);
        CompactMatrixAddMultiply('N', 'N', &TM, &Xdefl, &Hred, 1, 0, 0, 0);
        SquareCompactMatrixMakeHermitian(&Hred);
        FreeCompactMatrix(&TM);
        FreeCompactMatrix(&HextO);
        S->E = (double*)malloc(di * sizeof(double));
        if (S->E == NULL) { printf("malloc failed in CalculateGforHwithSigmaAsOperatorBathExcitations\n"); free(ms); free(seq); free(sorted); ierr = 1; goto cleanup; }
        if (HermitianCompactMatrixEigensystemRowVectors(&Hred, &ufun, S->E)) { free(ms); free(seq); free(sorted); ierr = 1; goto cleanup; }
        FreeCompactMatrix(&Hred);
        /*  U_{n,nu} = < n | w_nu > and C_{n,nu}, | n > = sum_nu C_{n,nu} | w_nu > */
        S->U.Ni = di; S->U.Nj = NRaw; S->U.QComplex = 1;
        S->C.Ni = di; S->C.Nj = NRaw; S->C.QComplex = 1;
        if (InitCompactMatrix(&(S->U)) || InitCompactMatrix(&(S->C))) { printf("InitCompactMatrix failed in CalculateGforHwithSigmaAsOperatorBathExcitations\n"); FreeCompactMatrix(&ufun); free(ms); free(seq); free(sorted); ierr = 1; goto cleanup; }
        CompactMatrixAddMultiply('*', 'N', &ufun, &Ws, &(S->U), 1, 0, 0, 0);
        CompactMatrixAddMultiply('N', 'T', &ufun, &Xdefl, &(S->C), 1, 0, 0, 0);
        FreeCompactMatrix(&ufun);
        FreeCompactMatrix(&Ws);
        FreeCompactMatrix(&Xdefl);
        /*  the couplings.  For every distinct entry i of S the parent is S\i and
         *  V[i][j] maps its raw vector (s', q) to the raw vector (s' + (i,j), q)
         *  of S, so Z_j = U[:, map] C_parent^T.                               */
        S->CouplingEntry = (uint32_t*)malloc(lev * sizeof(uint32_t));
        S->CouplingParent = (uint32_t*)malloc(lev * sizeof(uint32_t));
        S->Z = (CompactMatrixType**)calloc(lev, sizeof(CompactMatrixType*));
        if ((S->CouplingEntry == NULL) || (S->CouplingParent == NULL) || (S->Z == NULL)) { printf("malloc failed in CalculateGforHwithSigmaAsOperatorBathExcitations\n"); free(ms); free(seq); free(sorted); ierr = 1; goto cleanup; }
        S->NCoupling = 0;
        for (r = 0; r < lev; r++)
        {
          uint32_t *pms, ps = UINT32_MAX, w;
          BathExcSpaceType *Pp;
          if ((r > 0) && (ms[r] == ms[r - 1])) continue;
          i = ms[r];
          /*  the multiset S \ i                                                */
          pms = sorted;
          w = 0;
          for (k = 0; k < lev; k++) if (k != r) pms[w++] = ms[k];
          if (lev == 1)
          {
            ps = 0;
          }
          else
          {
            for (sp = 1; sp < NSpaceBefore; sp++)
            {
              if ((Space[sp].Level == lev - 1) && (memcmp(Space[sp].Entry, pms, (lev - 1) * sizeof(uint32_t)) == 0)) { ps = sp; break; }
            }
          }
          if (ps == UINT32_MAX) continue;
          Pp = &(Space[ps]);
          c = S->NCoupling;
          S->CouplingEntry[c] = i;
          S->CouplingParent[c] = ps;
          S->Z[c] = (CompactMatrixType*)calloc(Sigma[i].BlockSize, sizeof(CompactMatrixType));
          if (S->Z[c] == NULL) { printf("malloc failed in CalculateGforHwithSigmaAsOperatorBathExcitations\n"); free(ms); free(seq); free(sorted); ierr = 1; goto cleanup; }
          S->NCoupling++;
          for (j = 0; j < Sigma[i].BlockSize; j++)
          {
            uint32_t pair = 0, NRawP = Pp->NString * North, sidx;
            CompactMatrixType Usub;
            for (k = 0; k < NPair; k++) if ((PairEntry[k] == i) && (PairChannel[k] == j)) pair = k;
            CompactMatrixSetUnallocated(&(S->Z[c][j]));
            CompactMatrixSetUnallocated(&Usub);
            S->Z[c][j].Ni = di; S->Z[c][j].Nj = Pp->Dim; S->Z[c][j].QComplex = 1;
            Usub.Ni = di; Usub.Nj = NRawP; Usub.QComplex = 1;
            if (InitCompactMatrix(&(S->Z[c][j])) || InitCompactMatrix(&Usub)) { printf("InitCompactMatrix failed in CalculateGforHwithSigmaAsOperatorBathExcitations\n"); free(ms); free(seq); free(sorted); ierr = 1; goto cleanup; }
            for (sp = 0; sp < Pp->NString; sp++)
            {
              /*  the string of the parent with (i,j) applied last              */
              for (w = 0; w + 1 < lev; w++) seq[w] = Pp->String[(size_t)sp * (lev - 1) + w];
              seq[lev - 1] = pair;
              sidx = UINT32_MAX;
              for (pc = 0; pc < S->NString; pc++)
              {
                if (memcmp(&(S->String[(size_t)pc * lev]), seq, lev * sizeof(uint32_t)) == 0) { sidx = pc; break; }
              }
              if (sidx == UINT32_MAX) continue;
              for (d = 0; d < di; d++)
              {
                for (q = 0; q < North; q++)
                {
                  Usub.V[(size_t)d * NRawP + sp * North + q] = S->U.V[(size_t)d * NRaw + sidx * North + q];
                }
              }
            }
            if (lev == 1)
            {
              memcpy(S->Z[c][j].V, Usub.V, (size_t)di * North * sizeof(Complex16));
            }
            else
            {
              CompactMatrixAddMultiply('N', 'T', &Usub, &(Pp->C), &(S->Z[c][j]), 1, 0, 0, 0);
            }
            FreeCompactMatrix(&Usub);
          }
        }
        FreeCompactMatrix(&(S->U)); CompactMatrixSetUnallocated(&(S->U));
        if (S->NCoupling == 0)
        {
          VerbosityPrintf(0x0001, "CalculateGforHwithSigmaAsOperatorBathExcitations: space %u of level %u has no parent, skipped\n", NSpace, lev);
          BathExcSpaceFree(S, Sigma);
        }
        else
        {
          VerbosityPrintf(0x0001, "CalculateGforHwithSigmaAsOperatorBathExcitations: level %u space %u, %u strings, %u raw vectors deflated to %u\n", lev, NSpace, S->NString, NRaw, di);
          NSpace++;
        }
      }
      /*  next non decreasing sequence of entries                                */
      r = lev;
      more = false;
      while (r > 0)
      {
        r--;
        if (ms[r] + 1 < NSigma)
        {
          ms[r]++;
          for (k = r + 1; k < lev; k++) ms[k] = ms[r];
          more = true;
          break;
        }
      }
    }
    /*  the C matrices of level lev-1 are no longer needed                      */
    for (sp = 1; sp < NSpaceBefore; sp++)
    {
      if (Space[sp].Level + 1 == lev) { FreeCompactMatrix(&(Space[sp].C)); CompactMatrixSetUnallocated(&(Space[sp].C)); }
    }
    free(ms); free(seq); free(sorted);
  }
  for (sp = 1; sp < NSpace; sp++) { FreeCompactMatrix(&(Space[sp].C)); CompactMatrixSetUnallocated(&(Space[sp].C)); }

  /*  hand over                                                                */
  F->North = North;
  F->NFull = NFull;
  F->QComplexOut = QComplexOut;
  F->NSpace = NSpace;
  F->NSpaceAlloc = NSpaceAlloc;
  F->Parity = Parity;             Parity = NULL;
  F->HO = HO;                     CompactMatrixSetUnallocated(&HO);
  F->Cseed = Cseed;               CompactMatrixSetUnallocated(&Cseed);
  F->Space = Space;               Space = NULL;
  F->QHeffAllocated = QHeffAllocated;
  if (QHeffAllocated) F->Heff = Heff;
  QHeffAllocated = false;

cleanup:
  FreeCompactMatrix(&SMat);
  FreeCompactMatrix(&HMat);
  FreeCompactMatrix(&Sfun);
  FreeCompactMatrix(&T);
  FreeCompactMatrix(&TS);
  FreeCompactMatrix(&HO);
  FreeCompactMatrix(&Cseed);
  if (Space != NULL)
  {
    for (s = 0; s < NSpaceAlloc; s++) BathExcSpaceFree(&(Space[s]), Sigma);
    free(Space);
  }
  free(Sval);
  free(PairEntry);
  free(PairChannel);
  free(Parity);
  if (QHeffAllocated) FreeOperator(&Heff);
  return ierr;
}

/*
 *  Given an operator H and basis PsiBasis[i] with i in [0,NPsiBasis)
 *  Given hybridisation functions Sigma[i] with channel operators V[i][j]
 *  calculate
 *    G_{lm}(w) = < PsiBasis[l] | 1 / ( w - H_aux ) | PsiBasis[m] >
 *    with l,m in [0,BlockSize)
 *  where H_aux is H plus the baths the hybridisation functions stand for, in
 *  the space of at most NExcitations bath excitations, and return it in block
 *  tridiagonal form in G.
 *
 *  --------------------------------------------------------------------------
 *  What is calculated
 *  --------------------------------------------------------------------------
 *
 *  Sigma[i](w) = a0[i] + sum_n W[i][n] / ( w - a[i][n] ) is read as a bath of
 *  free particles.  W = L^dag L, and every row mu of L is a mode b_mu at energy
 *  a[i][n] coupled by
 *
 *    T = sum_mu sum_j ( L_{mu j} b_mu^dag V[i][j] + h.c. )
 *
 *  The coupling is a term of a Hamiltonian, so it is even in fermion parity and
 *  the bath particle has the parity of V: an odd V (c, c^dag c^dag c) couples
 *  to a fermion, an even V (c^dag c, a boson operator) to a boson.  All
 *  operators of one entry must have the same parity; that is checked.  At one
 *  excitation the statistics never shows, at two it sets the exchange sign of
 *  the crossed processes and whether a mode may hold two excitations.
 *
 *  The states are B(tau) | lambda > with lambda in X_m^S (see
 *  BathExcSpaceType) and B(tau) the bath creation operators of the tuple tau,
 *  the one of the largest mode index leftmost, normalized with 1/sqrt(n!) for a
 *  mode that holds n bosons.  Adding mode mu to tau' gives
 *
 *    b_mu^dag V_j B(tau') | eta > = (-1)^{pi_mu N_odd(tau',< mu)} sqrt(n_mu(tau))
 *                                   B(tau) V_j | eta >
 *
 *  with pi_mu the parity of the mode and N_odd(tau',< mu) the number of
 *  fermionic modes of tau' with an index below mu.  V_j | eta > is projected
 *  on X_m^S, which is the one approximation of the method (Rayleigh-Ritz, as
 *  in the old routine and in the paper), and the order of the system operators
 *  is the business of the operator algebra.
 *
 *  NonCrossing == false: tuples are multisets of modes and every mode of a
 *  tuple can be the one that is removed.  A state with two excitations is then
 *  reached along both orders, which gives the crossed processes, and the space
 *  is exact within at most NExcitations excitations.
 *
 *  NonCrossing == true: tuples are ordered histories and only the excitation
 *  added last can leave.  The excitations are distinguishable, there is no
 *  exchange sign, no Pauli exclusion and no Bose factor between them.  This is
 *  the nested (rainbow) approximation.  Its only use here is to measure how
 *  much the crossed processes contribute; the energy mesh method computes the
 *  same approximation far more cheaply.
 *
 *  --------------------------------------------------------------------------
 *  How it is calculated
 *  --------------------------------------------------------------------------
 *
 *  A,B  as in CalculateGforHwithSigmaAsOperator: a0 into H, P orthonormalized.
 *
 *  C    For m = 1 .. NExcitations and every multiset S of m entries: the
 *       overlap and Hamiltonian matrices of the raw vectors
 *       V[s_m]...V[s_1] | e_p > from operator products that act in P, so no
 *       wave function outside P is ever created.  Deflate, diagonalize.
 *       V[i][j] maps a raw vector of S\i into a raw vector of S, so the
 *       couplings Z follow from the eigenvectors without any further operator.
 *
 *  D    The modes, optionally binned on the grid.
 *
 *  E    The tuples and links of every block, dropping tuples that no link
 *       reaches.
 *
 *  F    Block Lanczos on the auxiliary Hamiltonian from the BlockSize seeds.
 *       Full re-orthogonalization while the vectors fit in
 *       BathExcMaxLanczosMemory, a local one against the two previous blocks
 *       beyond it.  The local one keeps the moments right and can split a
 *       converged pole into copies that share its weight.
 *
 *  G    The chain padded up to the block size of B0, as in the old routine.
 *
 *  --------------------------------------------------------------------------
 *  Arguments
 *  --------------------------------------------------------------------------
 *
 *  As CalculateGforHwithSigmaAsOperator, and
 *
 *  NExcitations         at most this many bath excitations at once, >= 1
 *  NonCrossing          see above
 *  SuggestedGrid        grid on which the poles of the hybridisation functions
 *                       are binned before they become modes, NULL for none.
 *                       (The old routine bins the poles of the self energy
 *                       after they are split by the intermediate states; that
 *                       can not be done here, the crossed processes need to
 *                       know which mode is which.)
 *
 *  returns 0 on success, 1 on failure
 */
int CalculateGforHwithSigmaAsOperatorBathExcitations(const OperatorType *const H, const WaveFunctionType *const PsiBasis, const uint32_t NPsiBasis, const uint32_t BlockSize, OperatorType ** V, const BlockListOfPoles *const Sigma, const uint32_t NSigma, const uint32_t NExcitations, const bool NonCrossing, const uint32_t NTriFull, const double *const SuggestedGrid, const uint32_t NSuggestedGridPoints, const CompactMatrixType *const B0, BlockTriDiagonalMatrixType * G)
{
  uint32_t i, j, k, p, q, s, t, r, c;
  uint32_t m = BlockSize, NFull, North, NTriUse, NKeep;
  uint32_t QComplexOut, NSpace;
  uint32_t *Parity;
  int ierr = 0;
  size_t Dtot, WorkSize = 0;
  double norm, Re;
  double DiscardedWeight = 0, NegativeWeight = 0;
  bool FullReorth;
  CompactMatrixType HO, Cseed, Ablk, Bblk, Lblk;
  BlockTriDiagonalMatrixType chain;
  BathExcSpaceType *Space;
  BathExcModesType Modes;
  BathExcSystemType Sys;
  BathExcFrontType Front;
  Complex16 *Qall = NULL, *Y = NULL, *Reorth = NULL;

  memset(&Modes, 0, sizeof(BathExcModesType));
  memset(&Sys, 0, sizeof(BathExcSystemType));
  CompactMatrixSetUnallocated(&Ablk);  CompactMatrixSetUnallocated(&Bblk);
  CompactMatrixSetUnallocated(&Lblk);
  chain.a = NULL;

  /* ---------------------------------------------------------------------- */
  /*  A, B, C  input, P, the system spaces                                  */
  /* ---------------------------------------------------------------------- */
  BathExcFrontEndInit(&Front);
  if (G == NULL) { printf("CalculateGforHwithSigmaAsOperatorBathExcitations called with a NULL pointer for G\n"); return 1; }
  if (BathExcFrontEndBuild(H, PsiBasis, NPsiBasis, BlockSize, V, Sigma, NSigma, NExcitations, NonCrossing, NTriFull, B0, &Front)) return 1;
  North = Front.North;
  NFull = Front.NFull;
  QComplexOut = Front.QComplexOut;
  NSpace = Front.NSpace;
  Parity = Front.Parity;
  Space = Front.Space;
  HO = Front.HO;          /* views, owned by Front */
  HO.OnLuaStack = true;
  Cseed = Front.Cseed;
  Cseed.OnLuaStack = true;

  /* ---------------------------------------------------------------------- */
  /*  D   the modes                                                          */
  /* ---------------------------------------------------------------------- */
  if (BathExcModesBuild(Sigma, NSigma, Parity, SuggestedGrid, NSuggestedGridPoints, &Modes, &DiscardedWeight, &NegativeWeight)) { ierr = 1; goto cleanup; }
  if (NegativeWeight > 0)
  {
    printf("CalculateGforHwithSigmaAsOperatorBathExcitations: %g negative spectral weight found in the residues of the hybridisation functions\n", NegativeWeight);
    printf("  A residue that is not positive semi definite has no Hamiltonian representation. The negative part is dropped.\n");
  }
  for (i = 0; i < NSigma; i++)
  {
    VerbosityPrintf(0x0001, "CalculateGforHwithSigmaAsOperatorBathExcitations: hybridisation %u, %s, %u modes\n", i + 1, Parity[i] ? "fermionic" : "bosonic", Modes.First[i + 1] - Modes.First[i]);
  }

  /* ---------------------------------------------------------------------- */
  /*  E   tuples and links                                                   */
  /* ---------------------------------------------------------------------- */
  Dtot = North;
  for (s = 1; s < NSpace; s++)
  {
    BathExcSpaceType *S = &(Space[s]);
    uint32_t L = S->Level, *tup, *lo = NULL, *hi = NULL, NTupleMax, NLinkMax, nkeep;
    size_t cnt;
    bool more;
    tup = (uint32_t*)calloc(L, sizeof(uint32_t));
    lo = (uint32_t*)calloc(L, sizeof(uint32_t));
    hi = (uint32_t*)calloc(L, sizeof(uint32_t));
    if ((tup == NULL) || (lo == NULL) || (hi == NULL)) { free(tup); free(lo); free(hi); printf("malloc failed in CalculateGforHwithSigmaAsOperatorBathExcitations\n"); ierr = 1; goto cleanup; }
    /*  count first.  Exact: position r runs over the modes of entry S[r],
     *  non decreasing, strictly increasing within a fermionic entry.
     *  NonCrossing: every position runs over all modes of all entries of S and
     *  the multiset of entries has to be S.                                   */
    cnt = 0;
    {
      uint32_t npass;
      for (npass = 0; npass < 2; npass++)
      {
        cnt = 0;
        /*  initialise the odometer                                           */
        more = true;
        for (r = 0; r < L; r++)
        {
          if (NonCrossing) { lo[r] = 0; hi[r] = Modes.N; }
          else { lo[r] = Modes.First[S->Entry[r]]; hi[r] = Modes.First[S->Entry[r] + 1]; }
        }
        r = 0;
        /*  depth first enumeration in lexicographic order                    */
        {
          int32_t pos = 0;
          uint32_t *cur = tup;
          cur[0] = lo[0];
          if (lo[0] >= hi[0]) more = false;
          while (more && (pos >= 0))
          {
            if (cur[pos] >= hi[pos])
            {
              pos--;
              if (pos >= 0) cur[pos]++;
              continue;
            }
            if (NonCrossing)
            {
              uint32_t en = Modes.Entry[cur[pos]];
              bool okentry = false;
              for (k = 0; k < L; k++) if (S->Entry[k] == en) okentry = true;
              if (!okentry) { cur[pos]++; continue; }
            }
            if (pos + 1 < L)
            {
              uint32_t start = lo[pos + 1];
              if (!NonCrossing && (S->Entry[pos + 1] == S->Entry[pos]))
              {
                start = cur[pos] + (Modes.Parity[S->Entry[pos]] ? 1 : 0);
              }
              pos++;
              cur[pos] = start;
              continue;
            }
            /*  a complete tuple                                                */
            {
              bool accept = true;
              if (NonCrossing)
              {
                /*  the multiset of entries has to be S                           */
                uint32_t e1[64];
                if (L > 64) accept = false;
                else
                {
                  for (k = 0; k < L; k++) e1[k] = Modes.Entry[cur[k]];
                  for (k = 1; k < L; k++)
                  {
                    uint32_t v = e1[k], w = k;
                    while ((w > 0) && (e1[w - 1] > v)) { e1[w] = e1[w - 1]; w--; }
                    e1[w] = v;
                  }
                  if (memcmp(e1, S->Entry, L * sizeof(uint32_t)) != 0) accept = false;
                }
              }
              if (accept)
              {
                if (npass == 1) memcpy(&(S->Tuple[cnt * L]), cur, L * sizeof(uint32_t));
                cnt++;
              }
            }
            cur[pos]++;
          }
        }
        if (npass == 0)
        {
          NTupleMax = (uint32_t)cnt;
          if (cnt > (size_t)UINT32_MAX / 2) { printf("CalculateGforHwithSigmaAsOperatorBathExcitations: too many bath states in space %u\n", s); free(tup); free(lo); free(hi); ierr = 1; goto cleanup; }
          S->Tuple = (uint32_t*)malloc(Max(cnt, (size_t)1) * L * sizeof(uint32_t));
          if (S->Tuple == NULL) { printf("malloc failed in CalculateGforHwithSigmaAsOperatorBathExcitations for %zu tuples\n", cnt); free(tup); free(lo); free(hi); ierr = 1; goto cleanup; }
        }
      }
    }
    S->NTuple = (uint32_t)cnt;
    free(tup); free(lo); free(hi);
    /*  links.  Exact: one per distinct mode of the tuple.  NonCrossing: the
     *  last one only.                                                          */
    NLinkMax = S->NTuple * (NonCrossing ? 1 : L);
    S->LinkTuple = (uint32_t*)malloc(Max(NLinkMax, 1) * sizeof(uint32_t));
    S->LinkParentTuple = (uint32_t*)malloc(Max(NLinkMax, 1) * sizeof(uint32_t));
    S->LinkMode = (uint32_t*)malloc(Max(NLinkMax, 1) * sizeof(uint32_t));
    S->LinkCoupling = (uint32_t*)malloc(Max(NLinkMax, 1) * sizeof(uint32_t));
    S->LinkFactor = (double*)malloc(Max(NLinkMax, 1) * sizeof(double));
    S->TupleEnergy = (double*)malloc(Max(S->NTuple, 1) * sizeof(double));
    if ((S->LinkTuple == NULL) || (S->LinkParentTuple == NULL) || (S->LinkMode == NULL) || (S->LinkCoupling == NULL) || (S->LinkFactor == NULL) || (S->TupleEnergy == NULL))
    {
      printf("malloc failed in CalculateGforHwithSigmaAsOperatorBathExcitations for the links of space %u\n", s);
      ierr = 1; goto cleanup;
    }
    S->NLink = 0;
    nkeep = 0;
    for (t = 0; t < S->NTuple; t++)
    {
      uint32_t *tt = &(S->Tuple[(size_t)t * L]), parent[64], rfirst, rlast, nlinkt = 0;
      double esum = 0;
      if (L > 64) { printf("CalculateGforHwithSigmaAsOperatorBathExcitations: more than 64 excitations\n"); ierr = 1; goto cleanup; }
      for (r = 0; r < L; r++) esum += Modes.E[tt[r]];
      rfirst = NonCrossing ? L - 1 : 0;
      rlast = L;
      for (r = rfirst; r < rlast; r++)
      {
        uint32_t mu = tt[r], w, nodd = 0, nmu = 0, ptup, ps;
        if (!NonCrossing && (r > 0) && (tt[r - 1] == mu)) continue;  /* one link per distinct mode */
        w = 0;
        for (k = 0; k < L; k++)
        {
          if (k == r) continue;
          parent[w++] = tt[k];
        }
        /*  the coupling of the entry of mu                                     */
        c = UINT32_MAX;
        for (k = 0; k < S->NCoupling; k++) if (S->CouplingEntry[k] == Modes.Entry[mu]) c = k;
        if (c == UINT32_MAX) continue;
        ps = S->CouplingParent[c];
        if (L == 1) ptup = 0;
        else ptup = BathExcTupleFind(Space[ps].Tuple, Space[ps].NTuple, parent, L - 1);
        if (ptup == UINT32_MAX) continue;
        if (NonCrossing)
        {
          S->LinkFactor[S->NLink] = 1;
        }
        else
        {
          for (k = 0; k < L; k++) if (tt[k] == mu) nmu++;
          for (k = 0; k < L - 1; k++) if ((parent[k] < mu) && Modes.Parity[Modes.Entry[parent[k]]]) nodd++;
          S->LinkFactor[S->NLink] = ((Modes.Parity[Modes.Entry[mu]] && (nodd % 2)) ? -1.0 : 1.0) * sqrt((double)nmu);
        }
        S->LinkTuple[S->NLink] = nkeep;
        S->LinkParentTuple[S->NLink] = ptup;
        S->LinkMode[S->NLink] = mu;
        S->LinkCoupling[S->NLink] = c;
        S->NLink++;
        nlinkt++;
      }
      /*  a tuple no link reaches never takes part, it is dropped               */
      if (nlinkt == 0) continue;
      if (nkeep != t) memmove(&(S->Tuple[(size_t)nkeep * L]), tt, L * sizeof(uint32_t));
      S->TupleEnergy[nkeep] = esum;
      nkeep++;
    }
    if (nkeep != S->NTuple)
    {
      VerbosityPrintf(0x0001, "CalculateGforHwithSigmaAsOperatorBathExcitations: space %u, %u of %u bath states are not reached and dropped\n", s, S->NTuple - nkeep, S->NTuple);
    }
    S->NTuple = nkeep;
    S->Offset = Dtot;
    Dtot += (size_t)S->NTuple * S->Dim;
    VerbosityPrintf(0x0001, "CalculateGforHwithSigmaAsOperatorBathExcitations: level %u space %u, %u x %u = %zu states, %u links\n", S->Level, s, S->Dim, S->NTuple, (size_t)S->NTuple * S->Dim, S->NLink);
    for (c = 0; c < S->NCoupling; c++)
    {
      BathExcSpaceType *Pp = &(Space[S->CouplingParent[c]]);
      WorkSize = Max(WorkSize, (size_t)Sigma[S->CouplingEntry[c]].BlockSize * Pp->NTuple * S->Dim);
    }
  }
  Space[0].Offset = 0;
  VerbosityPrintf(0x0001, "CalculateGforHwithSigmaAsOperatorBathExcitations: the auxiliary Hamiltonian has dimension %zu\n", Dtot);
  if (Dtot > (size_t)INT32_MAX)
  {
    printf("CalculateGforHwithSigmaAsOperatorBathExcitations: the auxiliary space has %zu states, more than the %d the linear algebra can index.\n", Dtot, INT32_MAX);
    printf("  Bin the hybridisation functions on a coarser grid (GridPoints) or reduce NTri.\n");
    ierr = 1; goto cleanup;
  }
  Sys.North = North;
  Sys.HO = &HO;
  Sys.NSpace = NSpace;
  Sys.Space = Space;
  Sys.Modes = &Modes;
  Sys.Sigma = Sigma;
  Sys.Dtot = Dtot;
  Sys.WorkSize = Max(WorkSize, (size_t)1);
  Sys.WorkU = (Complex16*)malloc(Sys.WorkSize * sizeof(Complex16));
  Sys.WorkW = (Complex16*)malloc(Sys.WorkSize * sizeof(Complex16));
  if ((Sys.WorkU == NULL) || (Sys.WorkW == NULL)) { printf("malloc failed in CalculateGforHwithSigmaAsOperatorBathExcitations for the work space\n"); ierr = 1; goto cleanup; }

  /* ---------------------------------------------------------------------- */
  /*  F   block Lanczos on the auxiliary Hamiltonian                         */
  /* ---------------------------------------------------------------------- */
  NTriUse = NTriFull;
  if ((double)NTriUse * (double)m > (double)Dtot) NTriUse = (uint32_t)((Dtot + m - 1) / m);
  if (NTriUse == 0) NTriUse = 1;
  if (NTriUse != NTriFull)
  {
    VerbosityPrintf(0x0001, "CalculateGforHwithSigmaAsOperatorBathExcitations: NTriFull reduced from %u to %u, the space has only %zu dimensions\n", NTriFull, NTriUse, Dtot);
  }
  FullReorth = ((double)NTriUse * (double)m * (double)Dtot <= BathExcMaxLanczosMemory);
  NKeep = FullReorth ? NTriUse : 2;
  if (!FullReorth)
  {
    VerbosityPrintf(0x0001, "CalculateGforHwithSigmaAsOperatorBathExcitations: %u Lanczos vectors of %zu states do not fit, re-orthogonalization only against the two previous blocks\n", NTriUse * m, Dtot);
  }
  Qall = (Complex16*)calloc((size_t)NKeep * m * Dtot, sizeof(Complex16));
  Y = (Complex16*)calloc((size_t)m * Dtot, sizeof(Complex16));
  Reorth = (Complex16*)calloc((size_t)m * NKeep * m, sizeof(Complex16));
  Ablk.Ni = m; Ablk.Nj = m; Ablk.QComplex = 1;
  Bblk.Ni = m; Bblk.Nj = m; Bblk.QComplex = 1;
  Lblk.Ni = m; Lblk.Nj = m; Lblk.QComplex = 1;
  if ((Qall == NULL) || (Y == NULL) || (Reorth == NULL) || InitCompactMatrix(&Ablk) || InitCompactMatrix(&Bblk) || InitCompactMatrix(&Lblk))
  {
    printf("malloc failed in CalculateGforHwithSigmaAsOperatorBathExcitations for %u Lanczos vectors of length %zu\n", NKeep * m, Dtot);
    ierr = 1; goto cleanup;
  }
  chain.N = NTriUse;
  chain.BlockSize = m;
  chain.QComplex = 1;
  if (InitBlockTridiagonalMatrix(&chain)) { printf("InitBlockTridiagonalMatrix failed in CalculateGforHwithSigmaAsOperatorBathExcitations\n"); chain.a = NULL; ierr = 1; goto cleanup; }
  chain.mu = 0;
  {
    CompactMatrixType Qk, Qj, Yv, Empty;
    uint32_t l, slot, slotprev, pass2;
    size_t rowlen = Dtot;
    Empty = BathExcView(NULL, m, 0);
    for (l = 0; l < m; l++)
    {
      for (p = 0; p < North; p++) Qall[(size_t)l * rowlen + p] = Cseed.V[p * m + l];
    }
    Qk = BathExcView(Qall, m, (uint32_t)Dtot);
    RowBlockPairOrthonormalizeWithTransform(&Qk, &Empty, &Lblk, BathExcEpsLanczos);
    memcpy(chain.b, Lblk.V, Sqr(m) * sizeof(Complex16));
    Yv = BathExcView(Y, m, (uint32_t)Dtot);
    for (k = 0; k < NTriUse; k++)
    {
      slot = k % NKeep;
      Qk = BathExcView(&(Qall[(size_t)slot * m * rowlen]), m, (uint32_t)Dtot);
      /*  Y = M Q_k                                                             */
      for (l = 0; l < m; l++) BathExcMatVec(&Sys, &(Qk.V[(size_t)l * rowlen]), &(Y[(size_t)l * rowlen]));
      /*  A_k = Q_k^dag M Q_k.  Y Q_k^dag is its transpose; written this way
       *  zgemm conjugates on the fly, where '*' would conjugate Q_k in place,
       *  twice, for every product                                            */
      CompactMatrixAddMultiply('N', 'C', &Yv, &Qk, &Bblk, 1, 0, 0, 0);
      for (p = 0; p < m; p++)
      {
        for (q = 0; q < m; q++) Ablk.V[p * m + q] = Bblk.V[q * m + p];
      }
      SquareCompactMatrixMakeHermitian(&Ablk);
      memcpy(&(chain.a[(size_t)(k + 1) * Sqr(m)]), Ablk.V, Sqr(m) * sizeof(Complex16));
      if (k + 1 == NTriUse) break;
      /*  the residual                                                          */
      CompactMatrixAddMultiply('T', 'N', &Ablk, &Qk, &Yv, -1, 0, 1, 0);
      if (k != 0)
      {
        slotprev = (k - 1) % NKeep;
        Qj = BathExcView(&(Qall[(size_t)slotprev * m * rowlen]), m, (uint32_t)Dtot);
        memcpy(Bblk.V, &(chain.b[(size_t)k * Sqr(m)]), Sqr(m) * sizeof(Complex16));
        CompactMatrixAddMultiply('T', 'N', &Bblk, &Qj, &Yv, -1, 0, 1, 0);
      }
      /*  re-orthogonalization, twice.  With full re-orthogonalization the
       *  blocks 0..k are contiguous in Qall and one product against all of
       *  them replaces k+1 small ones:  C = Y Q^dag,  Y -= C Q               */
      for (pass2 = 0; pass2 < 2; pass2++)
      {
        uint32_t jfirst = FullReorth ? 0 : ((k > 0) ? k - 1 : 0), nb;
        if (FullReorth)
        {
          CompactMatrixType Qprev, Coef;
          nb = (k + 1) * m;
          Qprev = BathExcView(Qall, nb, (uint32_t)Dtot);
          Coef = BathExcView(Reorth, m, nb);
          CompactMatrixAddMultiply('N', 'C', &Yv, &Qprev, &Coef, 1, 0, 0, 0);
          CompactMatrixAddMultiply('N', 'N', &Coef, &Qprev, &Yv, -1, 0, 1, 0);
        }
        else
        {
          for (j = jfirst; j <= k; j++)
          {
            CompactMatrixType Coef;
            Qj = BathExcView(&(Qall[(size_t)(j % NKeep) * m * rowlen]), m, (uint32_t)Dtot);
            Coef = BathExcView(Reorth, m, m);
            CompactMatrixAddMultiply('N', 'C', &Yv, &Qj, &Coef, 1, 0, 0, 0);
            CompactMatrixAddMultiply('N', 'N', &Coef, &Qj, &Yv, -1, 0, 1, 0);
          }
        }
      }
      norm = 0;
      for (l = 0; l < m; l++)
      {
        size_t z;
        Re = 0;
        for (z = 0; z < Dtot; z++) Re += Sqr(Y[(size_t)l * rowlen + z].Re) + Sqr(Y[(size_t)l * rowlen + z].Im);
        norm = Max(norm, sqrt(Re));
      }
      RowBlockPairOrthonormalizeWithTransform(&Yv, &Empty, &Lblk, BathExcEpsLanczos * Max(norm, DBL_MIN));
      slot = (k + 1) % NKeep;
      memcpy(&(Qall[(size_t)slot * m * rowlen]), Y, (size_t)m * rowlen * sizeof(Complex16));
      for (p = 0; p < m; p++)
      {
        for (q = 0; q < m; q++)
        {
          chain.b[(size_t)(k + 1) * Sqr(m) + p * m + q].Re =  Lblk.V[p * m + q].Re;
          chain.b[(size_t)(k + 1) * Sqr(m) + p * m + q].Im = -Lblk.V[p * m + q].Im;
        }
      }
    }
  }
  free(Qall); Qall = NULL;
  free(Y); Y = NULL;
  free(Reorth); Reorth = NULL;

  /* ---------------------------------------------------------------------- */
  /*  G   the chain padded back up to the full block size                    */
  /*      (the same as step G of CalculateGforHwithSigmaAsOperator)          */
  /* ---------------------------------------------------------------------- */
  G->N = chain.N;
  G->BlockSize = NFull;
  G->QComplex = 1;
  if (InitBlockTridiagonalMatrix(G)) { printf("InitBlockTridiagonalMatrix failed in CalculateGforHwithSigmaAsOperatorBathExcitations\n"); ierr = 1; goto cleanup; }
  G->mu = 0;
  strcpy(G->Name, "G with self energy");
  G->BlockSizeList[0] = NFull;
  for (k = 1; k <= G->N; k++) G->BlockSizeList[k] = m;
  for (k = 1; k <= G->N; k++)
  {
    for (p = 0; p < m; p++)
    {
      for (q = 0; q < m; q++) G->a[(size_t)k * Sqr(NFull) + p * NFull + q] = chain.a[(size_t)k * Sqr(m) + p * m + q];
    }
  }
  for (k = 1; k < G->N; k++)
  {
    for (p = 0; p < m; p++)
    {
      for (q = 0; q < m; q++) G->b[(size_t)k * Sqr(NFull) + p * NFull + q] = chain.b[(size_t)k * Sqr(m) + p * m + q];
    }
  }
  for (k = 0; k < Sqr(NFull); k++) { G->b[k].Re = 0; G->b[k].Im = 0; }
  for (p = 0; p < NFull; p++)
  {
    for (q = 0; q < m; q++)
    {
      for (k = 0; k < m; k++)
      {
        Complex16 b0pk;
        if (B0 == NULL) { b0pk.Re = (p == k) ? 1 : 0; b0pk.Im = 0; }
        else if (B0->QComplex) { b0pk = B0->V[p * m + k]; }
        else { b0pk.Re = B0->VRe[p * m + k]; b0pk.Im = 0; }
        G->b[p * NFull + q].Re += b0pk.Re * chain.b[k * m + q].Re - b0pk.Im * chain.b[k * m + q].Im;
        G->b[p * NFull + q].Im += b0pk.Re * chain.b[k * m + q].Im + b0pk.Im * chain.b[k * m + q].Re;
      }
    }
  }
  if (QComplexOut == 0)
  {
    norm = 0;
    Re = 0;
    for (k = 0; k < (G->N + 1) * Sqr(NFull); k++) { norm = Max(norm, fabs(G->a[k].Im)); Re = Max(Re, fabs(G->a[k].Re)); }
    for (k = 0; k < G->N * Sqr(NFull); k++)       { norm = Max(norm, fabs(G->b[k].Im)); Re = Max(Re, fabs(G->b[k].Re)); }
    if (norm > 1e-10 * Max(Re, 1.0))
    {
      printf("CalculateGforHwithSigmaAsOperatorBathExcitations: all input is real but the result has an imaginary part of %g, kept complex\n", norm);
    }
    else
    {
      ComplexBlockTridiagonalMatrixToReal(G);
    }
  }

cleanup:
  FreeCompactMatrix(&Ablk);
  FreeCompactMatrix(&Bblk);
  FreeCompactMatrix(&Lblk);
  if (chain.a != NULL) FreeBlockTridiagonalMatrix(&chain);
  BathExcModesFree(&Modes);
  free(Sys.WorkU);
  free(Sys.WorkW);
  free(Qall);
  free(Y);
  free(Reorth);
  BathExcFrontEndFree(&Front, Sigma);
  return ierr;
}

/* ===========================================================================
 *
 *   The same Green function on a fixed energy mesh, in the nested (non
 *   crossing) approximation
 *
 * =========================================================================== */

/*  the same for B series at once, stored as a[t*B + b]: every butterfly runs
 *  over the B contiguous series, and the twiddles tw[k] = exp(sign 2 pi i k/n),
 *  k < n/2, are computed once by the caller                                  */
static void BathExcFFTBatch(Complex16 *const a, uint32_t const n, size_t const B, Complex16 const *const tw)
{
  uint32_t i, j, k, len, half, step;
  size_t b;
  Complex16 *row = (Complex16*)malloc(B * sizeof(Complex16));
  if (row == NULL) return;
  for (i = 1, j = 0; i < n; i++)
  {
    uint32_t bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j)
    {
      memcpy(row, &(a[(size_t)i * B]), B * sizeof(Complex16));
      memcpy(&(a[(size_t)i * B]), &(a[(size_t)j * B]), B * sizeof(Complex16));
      memcpy(&(a[(size_t)j * B]), row, B * sizeof(Complex16));
    }
  }
  free(row);
  for (len = 2; len <= n; len <<= 1)
  {
    half = len >> 1;
    step = n / len;
    for (i = 0; i < n; i += len)
    {
      for (k = 0; k < half; k++)
      {
        Complex16 w = tw[(size_t)k * step];
        Complex16 *u = &(a[(size_t)(i + k) * B]), *v = &(a[(size_t)(i + k + half) * B]);
        for (b = 0; b < B; b++)
        {
          double vr = v[b].Re * w.Re - v[b].Im * w.Im;
          double vi = v[b].Re * w.Im + v[b].Im * w.Re;
          v[b].Re = u[b].Re - vr;  v[b].Im = u[b].Im - vi;
          u[b].Re += vr;           u[b].Im += vi;
        }
      }
    }
  }
}

/*  the smallest power of two >= n                                              */
static uint32_t BathExcPow2(size_t const n)
{
  uint32_t p = 1;
  while (p < n) p <<= 1;
  return p;
}

/*  M = z - Hdiag - Sigma, inverted in place.  Hdiag is either a full Dim x Dim
 *  matrix (Hfull != NULL) or a list of Dim energies                          */
static int BathExcResolvent(Complex16 const z, uint32_t const Dim, CompactMatrixType const *const Hfull, double const *const Ediag, Complex16 const *const Sig, CompactMatrixType *const Mat)
{
  uint32_t a, b;
  for (a = 0; a < Dim; a++)
  {
    for (b = 0; b < Dim; b++)
    {
      Complex16 v;
      v.Re = -Sig[a * Dim + b].Re;
      v.Im = -Sig[a * Dim + b].Im;
      if (Hfull != NULL) { v.Re -= Hfull->V[a * Dim + b].Re; v.Im -= Hfull->V[a * Dim + b].Im; }
      if (a == b)
      {
        v.Re += z.Re;
        v.Im += z.Im;
        if (Ediag != NULL) v.Re -= Ediag[a];
      }
      Mat->V[a * Dim + b] = v;
    }
  }
  if (CompactMatrixInverse(Mat) < 0) return 1;
  return 0;
}

/*  what the evaluation at one y needs, all of it read only                     */
typedef struct BathExcMeshType
{
  uint32_t p, m, North, NSpace;
  BathExcSpaceType const *Space;
  BlockListOfPoles const *Sigma;
  CompactMatrixType const *HOrel, *Cseed;
  Complex16 ***Wbin;
  int64_t const *kL, *kH;
  int64_t jL, jH, kOutL;
  size_t NOut;
  double dw;
}
BathExcMeshType;

/*  The nested evaluation at z = w_k + i y for every mesh point k, level by
 *  level from the top down, and its contribution wq y [A(k-1) - 2A(k) + A(k+1)]
 *  to Racc.  Allocates its own work space, so several y can run at once.      */
static int BathExcMeshOneY(BathExcMeshType const *const X, double const y, double const wq, Complex16 *const Racc)
{
  uint32_t i, l, lp, k, b, s, c, n, lev, a;
  uint32_t p = X->p, m = X->m, North = X->North, NSpace = X->NSpace, DimMax = North;
  BathExcSpaceType const *Space = X->Space;
  BlockListOfPoles const *Sigma = X->Sigma;
  int64_t const *kL = X->kL, *kH = X->kH;
  int64_t jL = X->jL, jH = X->jH;
  double dw = X->dw;
  size_t GmeshMax = 0;
  int ierr = 0;
  Complex16 **Gmesh = NULL, *Gseed = NULL, *Sig = NULL, *Kser = NULL, *Fhat = NULL, *What = NULL, *Xbuf = NULL, *Kt = NULL;
  CompactMatrixType Mat, Tmp;
  CompactMatrixSetUnallocated(&Mat);
  CompactMatrixSetUnallocated(&Tmp);
  for (s = 1; s < NSpace; s++) DimMax = Max(DimMax, Space[s].Dim);
  Gmesh = (Complex16**)calloc(NSpace, sizeof(Complex16*));
  if (Gmesh == NULL) return 1;
  for (s = 1; s < NSpace; s++)
  {
    lev = Space[s].Level;
    if (lev >= p) continue;
    Gmesh[s] = (Complex16*)malloc((size_t)(kH[lev] - kL[lev] + 1) * Sqr(Space[s].Dim) * sizeof(Complex16));
    if (Gmesh[s] == NULL) { printf("malloc failed in CalculateGforHwithSigmaAsOperatorBathExcitationsMesh for the mesh of space %u\n", s); ierr = 1; goto done; }
    GmeshMax = Max(GmeshMax, (size_t)(kH[lev] - kL[lev] + 1));
  }
  Gseed = (Complex16*)malloc((size_t)(kH[0] - kL[0] + 1) * m * m * sizeof(Complex16));
  Sig = (Complex16*)malloc((size_t)Max((size_t)(kH[0] - kL[0] + 1), GmeshMax) * Sqr(DimMax) * sizeof(Complex16));
  Mat.Ni = DimMax; Mat.Nj = DimMax; Mat.QComplex = 1;
  Tmp.Ni = North; Tmp.Nj = m; Tmp.QComplex = 1;
  if ((Gseed == NULL) || (Sig == NULL) || InitCompactMatrix(&Mat) || InitCompactMatrix(&Tmp)) { printf("malloc failed in CalculateGforHwithSigmaAsOperatorBathExcitationsMesh\n"); ierr = 1; goto done; }
  for (int32_t levi = (int32_t)p - 1; levi >= 0; levi--)
  {
    lev = (uint32_t)levi;
    for (s = 0; s < NSpace; s++)
    {
      BathExcSpaceType const *S = &(Space[s]);
      uint32_t DimS = (s == 0) ? North : S->Dim;
      size_t Nk = (size_t)(kH[lev] - kL[lev] + 1), kk;
      if (S->Level != lev) continue;
      memset(Sig, 0, Nk * DimS * DimS * sizeof(Complex16));
      /*  the self energy from every space one level up                      */
      for (uint32_t sc = 1; sc < NSpace; sc++)
      {
        BathExcSpaceType const *C = &(Space[sc]);
        if (C->Level != lev + 1) continue;
        for (c = 0; c < C->NCoupling; c++)
        {
          uint32_t bi, DimC = C->Dim;
          if (C->CouplingParent[c] != s) continue;
          i = C->CouplingEntry[c];
          bi = Sigma[i].BlockSize;
          if (C->Level == p)
          {
            /*  top level: exact.  Sigma += sum_l Z_l^dag T_l with
             *  T_l[lambda,:] = sum_l' Delta_ll'(z - E_C[lambda]) Z_l'[lambda,:]  */
            uint32_t stride = 1 + bi * bi;
            Complex16 *Tl = (Complex16*)malloc((size_t)bi * DimC * DimS * sizeof(Complex16));
            Complex16 *dl = (Complex16*)malloc(bi * bi * sizeof(Complex16));
            if ((Tl == NULL) || (dl == NULL)) { free(Tl); free(dl); ierr = 1; goto done; }
            for (kk = 0; kk < Nk; kk++)
            {
              double x = (double)(kL[lev] + (int64_t)kk) * dw;
              CompactMatrixType Sv = BathExcView(&(Sig[kk * DimS * DimS]), DimS, DimS);
              memset(Tl, 0, (size_t)bi * DimC * DimS * sizeof(Complex16));
              for (uint32_t lam = 0; lam < DimC; lam++)
              {
                memset(dl, 0, bi * bi * sizeof(Complex16));
                for (n = 0; n < Sigma[i].N; n++)
                {
                  double en = Sigma[i].QComplex ? Sigma[i].V[n * stride].Re : Sigma[i].VRe[n * stride];
                  double dr = x - C->E[lam] - en, den = dr * dr + y * y;
                  Complex16 inv;
                  inv.Re = dr / den;
                  inv.Im = -y / den;
                  for (k = 0; k < bi * bi; k++)
                  {
                    Complex16 w;
                    if (Sigma[i].QComplex) w = Sigma[i].V[n * stride + 1 + k];
                    else { w.Re = Sigma[i].VRe[n * stride + 1 + k]; w.Im = 0; }
                    dl[k].Re += w.Re * inv.Re - w.Im * inv.Im;
                    dl[k].Im += w.Re * inv.Im + w.Im * inv.Re;
                  }
                }
                for (l = 0; l < bi; l++)
                {
                  for (lp = 0; lp < bi; lp++)
                  {
                    Complex16 d = dl[l * bi + lp];
                    if ((d.Re == 0) && (d.Im == 0)) continue;
                    Complex16 const *zr = &(C->Z[c][lp].V[(size_t)lam * DimS]);
                    Complex16 *tr = &(Tl[((size_t)l * DimC + lam) * DimS]);
                    for (b = 0; b < DimS; b++)
                    {
                      tr[b].Re += d.Re * zr[b].Re - d.Im * zr[b].Im;
                      tr[b].Im += d.Re * zr[b].Im + d.Im * zr[b].Re;
                    }
                  }
                }
              }
              for (l = 0; l < bi; l++)
              {
                CompactMatrixType Tv = BathExcView(&(Tl[(size_t)l * DimC * DimS]), DimC, DimS);
                CompactMatrixAddMultiply('C', 'N', &(C->Z[c][l]), &Tv, &Sv, 1, 0, 1, 0);
              }
            }
            free(Tl);
            free(dl);
          }
          else
          {
            /*  below the top level: the binned hybridisation convolved with
             *  K_ll'(x) = Z_l^dag G_C(x) Z_l' along the real axis, by FFT      */
            size_t Nc = (size_t)(kH[lev + 1] - kL[lev + 1] + 1), Nw = (size_t)(jH - jL + 1), t, DD = (size_t)DimS * DimS;
            uint32_t Lf = BathExcPow2(Nc + Nw - 1), ab;
            Complex16 *twf, *twb;
            Kser = (Complex16*)malloc((size_t)Lf * DD * sizeof(Complex16));
            Fhat = (Complex16*)calloc((size_t)Lf * DD, sizeof(Complex16));
            What = (Complex16*)malloc((size_t)Lf * sizeof(Complex16));
            Xbuf = (Complex16*)malloc(Nc * DimC * DimS * sizeof(Complex16));
            Kt = (Complex16*)malloc((size_t)Lf * sizeof(Complex16));   /* the twiddles, both signs */
            if ((Kser == NULL) || (Fhat == NULL) || (What == NULL) || (Xbuf == NULL) || (Kt == NULL)) { printf("malloc failed in CalculateGforHwithSigmaAsOperatorBathExcitationsMesh for the convolution\n"); ierr = 1; goto done; }
            twf = Kt;
            twb = &(Kt[Lf / 2]);
            for (t = 0; t < Lf / 2; t++)
            {
              twf[t].Re = cos(2 * M_PI * t / Lf);  twf[t].Im = -sin(2 * M_PI * t / Lf);
              twb[t].Re = twf[t].Re;               twb[t].Im = -twf[t].Im;
            }
            for (lp = 0; lp < bi; lp++)
            {
              /*  X(x) = G_C(x) Z_l'                                              */
              for (t = 0; t < Nc; t++)
              {
                CompactMatrixType Gv = BathExcView(&(Gmesh[sc][t * DimC * DimC]), DimC, DimC);
                CompactMatrixType Xv = BathExcView(&(Xbuf[t * DimC * DimS]), DimC, DimS);
                CompactMatrixAddMultiply('N', 'N', &Gv, &(C->Z[c][lp]), &Xv, 1, 0, 0, 0);
              }
              for (l = 0; l < bi; l++)
              {
                bool any = false;
                for (t = 0; t < Nw; t++) if ((X->Wbin[i][t][l * bi + lp].Re != 0) || (X->Wbin[i][t][l * bi + lp].Im != 0)) any = true;
                if (!any) continue;   /* e.g. the spin off diagonal blocks */
                /*  K_ll'(x) = Z_l^dag X(x), block t of the series, zero padded  */
                for (t = 0; t < Nc; t++)
                {
                  CompactMatrixType Kv = BathExcView(&(Kser[t * DD]), DimS, DimS);
                  CompactMatrixType Xv = BathExcView(&(Xbuf[t * DimC * DimS]), DimC, DimS);
                  CompactMatrixAddMultiply('C', 'N', &(C->Z[c][l]), &Xv, &Kv, 1, 0, 0, 0);
                }
                memset(&(Kser[Nc * DD]), 0, (Lf - Nc) * DD * sizeof(Complex16));
                memset(What, 0, (size_t)Lf * sizeof(Complex16));
                for (t = 0; t < Nw; t++) What[t] = X->Wbin[i][t][l * bi + lp];
                BathExcFFTBatch(What, Lf, 1, twf);
                BathExcFFTBatch(Kser, Lf, DD, twf);
                for (t = 0; t < Lf; t++)
                {
                  Complex16 wt = What[t], *ks = &(Kser[t * DD]), *fh = &(Fhat[t * DD]);
                  for (ab = 0; ab < DD; ab++)
                  {
                    fh[ab].Re += ks[ab].Re * wt.Re - ks[ab].Im * wt.Im;
                    fh[ab].Im += ks[ab].Re * wt.Im + ks[ab].Im * wt.Re;
                  }
                }
              }
            }
            BathExcFFTBatch(Fhat, Lf, DD, twb);
            for (kk = 0; kk < Nk; kk++)
            {
              /*  F(k) sits at t = k - kL[lev+1] - jL                          */
              t = (size_t)(kL[lev] + (int64_t)kk - kL[lev + 1] - jL);
              for (ab = 0; ab < DD; ab++)
              {
                Sig[kk * DD + ab].Re += Fhat[t * DD + ab].Re / Lf;
                Sig[kk * DD + ab].Im += Fhat[t * DD + ab].Im / Lf;
              }
            }
            free(Kser); Kser = NULL;
            free(Fhat); Fhat = NULL;
            free(What); What = NULL;
            free(Xbuf); Xbuf = NULL;
            free(Kt); Kt = NULL;
          }
        }
      }
      /*  the resolvent at every mesh point                                  */
      Mat.Ni = DimS; Mat.Nj = DimS;
      for (kk = 0; kk < Nk; kk++)
      {
        Complex16 z;
        z.Re = (double)(kL[lev] + (int64_t)kk) * dw;
        z.Im = y;
        if (BathExcResolvent(z, DimS, (s == 0) ? X->HOrel : NULL, (s == 0) ? NULL : S->E, &(Sig[kk * DimS * DimS]), &Mat))
        {
          printf("CalculateGforHwithSigmaAsOperatorBathExcitationsMesh: the resolvent is singular at w = %g + i %g\n", z.Re, z.Im);
          ierr = 1; goto done;
        }
        if (s == 0)
        {
          /*  the seed block Cseed^dag G_P Cseed                             */
          CompactMatrixType Gs = BathExcView(&(Gseed[kk * m * m]), m, m);
          CompactMatrixAddMultiply('N', 'N', &Mat, X->Cseed, &Tmp, 1, 0, 0, 0);
          CompactMatrixAddMultiply('C', 'N', X->Cseed, &Tmp, &Gs, 1, 0, 0, 0);
        }
        else
        {
          memcpy(&(Gmesh[s][kk * DimS * DimS]), Mat.V, (size_t)DimS * DimS * sizeof(Complex16));
        }
      }
    }
  }
  /*  R_k += w_q y [A(k-1) - 2 A(k) + A(k+1)],  A = -(1/pi)(G - G^dag)/(2i)  */
  for (size_t ko = 0; ko < X->NOut; ko++)
  {
    size_t kk = (size_t)(X->kOutL + (int64_t)ko - kL[0]);
    for (a = 0; a < m; a++)
    {
      for (b = 0; b < m; b++)
      {
        double sr = 0, si = 0, cf;
        int dk;
        for (dk = -1; dk <= 1; dk++)
        {
          Complex16 gab = Gseed[(kk + dk) * m * m + a * m + b];
          Complex16 gba = Gseed[(kk + dk) * m * m + b * m + a];
          double dr = gab.Re - gba.Re, di = gab.Im + gba.Im;   /* G_ab - conj(G_ba) */
          cf = (dk == 0) ? -2.0 : 1.0;
          sr += cf * (-di / (2 * M_PI));
          si += cf * ( dr / (2 * M_PI));
        }
        Racc[ko * m * m + a * m + b].Re += wq * y * sr;
        Racc[ko * m * m + a * m + b].Im += wq * y * si;
      }
    }
  }

done:
  if (Gmesh != NULL)
  {
    for (s = 0; s < NSpace; s++) free(Gmesh[s]);
    free(Gmesh);
  }
  free(Gseed);
  free(Sig);
  free(Kser);
  free(Fhat);
  free(What);
  free(Xbuf);
  free(Kt);
  FreeCompactMatrix(&Mat);
  FreeCompactMatrix(&Tmp);
  return ierr;
}

/*
 *  Given the same input as CalculateGforHwithSigmaAsOperatorBathExcitations
 *  calculate the seed block of G(w) on the uniform mesh w_k = k dw,
 *  dw = (Emax - Emin)/NE, energies measured from E0, in the nested (non
 *  crossing) approximation with up to NExcitations bath excitations, and
 *  return it as a block list of poles that sit on the mesh points.
 *
 *  --------------------------------------------------------------------------
 *  What is calculated
 *  --------------------------------------------------------------------------
 *
 *  With distinguishable bath excitations of which the last one added leaves
 *  first, the resolvent below a state of X_m^S with bath energy E_tau depends
 *  on tau only through z - E_tau.  So one matrix function per system space
 *  describes everything:
 *
 *    G_S(z)     = ( z - E_S - Sigma_S(z) )^-1                (G_P with H_P)
 *    Sigma_S(z) = sum_C sum_{l,l'} sum_n W_i[n]_{ll'} Z_l^dag G_C(z - e_i[n]) Z_l'
 *
 *  the sum over the spaces C = S + i one level up, Z_l = < C | V[i][l] | S >.
 *  At the top level G_C is diagonal and the sum over the poles of the
 *  hybridisation is done exactly, with the hybridisation evaluated at
 *  z - E_C[lambda].  At one excitation this is the whole calculation and the
 *  result is exact.  Below the top level the poles are binned on the lattice
 *  j dw (conserving the zeroth and the first moment of every bin), so
 *  z - e_i[n] is again a lattice point and the sum over n is a convolution
 *  along the real axis, done by FFT.  The binning keeps the moments M0 .. M3
 *  of G_P.  This is the same approximation as
 *  CalculateGforHwithSigmaAsOperatorBathExcitations with NonCrossing = true.
 *
 *  The result is the list of poles at w_k with residues
 *
 *    R_k = integral tent_k(w) A(w) dw ,   A = -(1/pi) Im G(w + i0)
 *
 *  with tent_k the triangle of width 2 dw at w_k: every pole of the exact G is
 *  split over its two neighbouring mesh points, conserving its weight and its
 *  first moment, exactly as the binning on a grid elsewhere in Quanty.  M0 and
 *  M1 are exact, M2 is off by at most dw^2/4 per unit of weight.
 *
 *  --------------------------------------------------------------------------
 *  How the residues are obtained without a broadening
 *  --------------------------------------------------------------------------
 *
 *  A_y(x) = -(1/pi) Im G(x + i y) is the spectrum broadened by a Lorentzian of
 *  half width y.  It is harmonic in (x,y), so d_x^2 A_y = - d_y^2 A_y, and two
 *  integrations by parts give int_0^inf y d_x^2 A_y dy = - A(x).  With the
 *  second difference on the mesh in place of d_x^2 this is exact:
 *
 *    R_k = -(1/dw) int_0^inf y [ A_y(w_{k-1}) - 2 A_y(w_k) + A_y(w_{k+1}) ] dy
 *
 *  (for one pole: y A_y = (W/pi)(1 - d^2/(d^2+y^2)), the 1 drops out of the
 *  second difference and the rest integrates to the tent.)  Only Im G enters.
 *  The integrand is bounded, W/pi at most, and a pole at distance d from a
 *  mesh point gives structure at y ~ d, so the integral is done on a log grid
 *  from y_min = BathExcMeshYmin dw, which errs by about 0.2 y_min / dw
 *  pointwise.  Beyond Y = BathExcMeshYmaxFactor times the spectral radius the
 *  moment expansion G = M0/z + M1/z^2 + M2/z^3 is integrated analytically.
 *  M0, M1 and M2 are exact:  S,  < H_P >  and  < H_P^2 > + < V^dag V >.
 *  Measured on 300 random poles with some on and next to mesh points: 2e-5
 *  pointwise, 1e-6 in M0 and M1, with 64 nodes.
 *
 *  For a matrix G, Im G is the Hermitian part (G - G^dag)/(2i); the residues
 *  come out Hermitian, and what the quadrature leaves negative is removed and
 *  reported.
 *
 *  --------------------------------------------------------------------------
 *  The mesh
 *  --------------------------------------------------------------------------
 *
 *  Every pole of the result has to sit inside the mesh, so the mesh covers a
 *  bound on the spectrum as well as [Emin, Emax]:  the extremes of the
 *  diagonal blocks widened by the norm of the coupling, which for the nested
 *  tree is at most twice max_S || sum_C Z^dag W_tot Z ||^(1/2).  The spaces one
 *  level up are needed at z - e and live on the mesh shifted by the range of
 *  the bins.
 *
 *  --------------------------------------------------------------------------
 *  Cost
 *  --------------------------------------------------------------------------
 *
 *  Per node of the y quadrature: one inversion per mesh point and system space
 *  below the top level, the exact top level sums (mesh points x poles of the
 *  hybridisation x dimension of the top spaces), and the FFT convolutions.
 *  Memory: the matrix functions of one y at a time.
 *
 *  returns 0 on success, 1 on failure
 */
static double const BathExcMeshYmin = 1e-4;
static double const BathExcMeshYmaxFactor = 100;
static double const BathExcMeshPanelWidth = 2.5;
static double const BathExcMeshKeepWeight = 1e-14;

int CalculateGforHwithSigmaAsOperatorBathExcitationsMesh(const OperatorType *const H, const WaveFunctionType *const PsiBasis, const uint32_t NPsiBasis, const uint32_t BlockSize, OperatorType ** V, const BlockListOfPoles *const Sigma, const uint32_t NSigma, const uint32_t NExcitations, const double E0, const double Emin, const double Emax, const uint32_t NE, const CompactMatrixType *const B0, BlockListOfPoles * G)
{
  /*  8 point Gauss-Legendre on [-1,1]                                        */
  static double const GLx[8] = {-0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498,
                                 0.1834346424956498,  0.5255324099163290,  0.7966664774136267,  0.9602898564975363};
  static double const GLw[8] = { 0.1012285362903763,  0.2223810344533745,  0.3137066458778873,  0.3626837833783620,
                                 0.3626837833783620,  0.3137066458778873,  0.2223810344533745,  0.1012285362903763};
  uint32_t i, j, l, lp, k, a, b, s, c, lev, q, n;
  uint32_t m = BlockSize, NFull, North, NSpace, p = NExcitations;
  int ierr = 0;
  int64_t kL[64], kH[64], jL = 0, jH = 0, kOutL, kOutH;
  double dw, lo, hi, Rspec, ymin, ymax, *emin = NULL, *emax = NULL;
  double dlo[64], dhi[64], bnorm[65];
  double ClippedWeight = 0;
  BathExcFrontType Front;
  BathExcSpaceType *Space;
  CompactMatrixType HOrel, Cseed, Tmp, M0, M1, M2;
  Complex16 *Rsum = NULL;
  Complex16 ***Wbin = NULL;         /* Wbin[i][j - jL] : bi x bi                 */
  double *ynode = NULL, *yweight = NULL;
  uint32_t Ny = 0;
  size_t NOut;

  BathExcFrontEndInit(&Front);
  CompactMatrixSetUnallocated(&HOrel);
  CompactMatrixSetUnallocated(&Tmp);   CompactMatrixSetUnallocated(&M0);
  CompactMatrixSetUnallocated(&M1);    CompactMatrixSetUnallocated(&M2);
  if (G == NULL) { printf("CalculateGforHwithSigmaAsOperatorBathExcitationsMesh called with a NULL pointer for G\n"); return 1; }
  if ((NE == 0) || !(Emax > Emin) || (p == 0) || (p > 60))
  {
    printf("CalculateGforHwithSigmaAsOperatorBathExcitationsMesh called with Emin=%g Emax=%g NE=%u NExcitations=%u\n", Emin, Emax, NE, p);
    return 1;
  }
  dw = (Emax - Emin) / (double)NE;

  /* ---------------------------------------------------------------------- */
  /*  A, B, C  input, P, the system spaces                                  */
  /* ---------------------------------------------------------------------- */
  if (BathExcFrontEndBuild(H, PsiBasis, NPsiBasis, BlockSize, V, Sigma, NSigma, NExcitations, true, 1, B0, &Front)) return 1;
  North = Front.North;
  NFull = Front.NFull;
  NSpace = Front.NSpace;
  Space = Front.Space;
  Cseed = Front.Cseed;
  Cseed.OnLuaStack = true;
  /*  energies from E0 on                                                     */
  HOrel.Ni = North; HOrel.Nj = North; HOrel.QComplex = 1;
  if (InitCompactMatrix(&HOrel)) { ierr = 1; goto cleanup; }
  memcpy(HOrel.V, Front.HO.V, (size_t)North * North * sizeof(Complex16));
  for (a = 0; a < North; a++) HOrel.V[a * North + a].Re -= E0;
  for (s = 1; s < NSpace; s++)
  {
    for (a = 0; a < Space[s].Dim; a++) Space[s].E[a] -= E0;
  }

  /* ---------------------------------------------------------------------- */
  /*  D   the poles of the hybridisations binned on the lattice j dw         */
  /* ---------------------------------------------------------------------- */
  emin = (double*)malloc(Max(NSigma, 1) * sizeof(double));
  emax = (double*)malloc(Max(NSigma, 1) * sizeof(double));
  if ((emin == NULL) || (emax == NULL)) { printf("malloc failed in CalculateGforHwithSigmaAsOperatorBathExcitationsMesh\n"); ierr = 1; goto cleanup; }
  jL = INT64_MAX; jH = INT64_MIN;
  for (i = 0; i < NSigma; i++)
  {
    uint32_t bi = Sigma[i].BlockSize, stride = 1 + bi * bi;
    emin[i] = DBL_MAX; emax[i] = -DBL_MAX;
    for (n = 0; n < Sigma[i].N; n++)
    {
      double en = Sigma[i].QComplex ? Sigma[i].V[n * stride].Re : Sigma[i].VRe[n * stride];
      emin[i] = Min(emin[i], en);
      emax[i] = Max(emax[i], en);
    }
    if (Sigma[i].N == 0) { emin[i] = 0; emax[i] = 0; }
    jL = Min(jL, (int64_t)floor(emin[i] / dw));
    jH = Max(jH, (int64_t)floor(emax[i] / dw) + 1);
  }
  if (NSigma == 0) { jL = 0; jH = 0; }
  Wbin = (Complex16***)calloc(Max(NSigma, 1), sizeof(Complex16**));
  if (Wbin == NULL) { printf("malloc failed in CalculateGforHwithSigmaAsOperatorBathExcitationsMesh\n"); ierr = 1; goto cleanup; }
  for (i = 0; i < NSigma; i++)
  {
    uint32_t bi = Sigma[i].BlockSize, stride = 1 + bi * bi;
    Wbin[i] = (Complex16**)calloc((size_t)(jH - jL + 1), sizeof(Complex16*));
    if (Wbin[i] == NULL) { printf("malloc failed in CalculateGforHwithSigmaAsOperatorBathExcitationsMesh\n"); ierr = 1; goto cleanup; }
    for (j = 0; j < (uint32_t)(jH - jL + 1); j++)
    {
      Wbin[i][j] = (Complex16*)calloc(bi * bi, sizeof(Complex16));
      if (Wbin[i][j] == NULL) { printf("malloc failed in CalculateGforHwithSigmaAsOperatorBathExcitationsMesh\n"); ierr = 1; goto cleanup; }
    }
    for (n = 0; n < Sigma[i].N; n++)
    {
      double en = Sigma[i].QComplex ? Sigma[i].V[n * stride].Re : Sigma[i].VRe[n * stride];
      double u = en / dw, f;
      int64_t j0 = (int64_t)floor(u);
      f = u - (double)j0;
      for (k = 0; k < bi * bi; k++)
      {
        Complex16 w;
        if (Sigma[i].QComplex) w = Sigma[i].V[n * stride + 1 + k];
        else { w.Re = Sigma[i].VRe[n * stride + 1 + k]; w.Im = 0; }
        Wbin[i][j0 - jL][k].Re += (1 - f) * w.Re;  Wbin[i][j0 - jL][k].Im += (1 - f) * w.Im;
        Wbin[i][j0 + 1 - jL][k].Re += f * w.Re;    Wbin[i][j0 + 1 - jL][k].Im += f * w.Im;
      }
    }
  }

  /* ---------------------------------------------------------------------- */
  /*  E   the moments of the seed block and a bound on the spectrum          */
  /* ---------------------------------------------------------------------- */
  M0.Ni = m; M0.Nj = m; M0.QComplex = 1;
  M1.Ni = m; M1.Nj = m; M1.QComplex = 1;
  M2.Ni = m; M2.Nj = m; M2.QComplex = 1;
  Tmp.Ni = North; Tmp.Nj = m; Tmp.QComplex = 1;
  if (InitCompactMatrix(&M0) || InitCompactMatrix(&M1) || InitCompactMatrix(&M2) || InitCompactMatrix(&Tmp)) { ierr = 1; goto cleanup; }
  CompactMatrixAddMultiply('C', 'N', &Cseed, &Cseed, &M0, 1, 0, 0, 0);
  CompactMatrixAddMultiply('N', 'N', &HOrel, &Cseed, &Tmp, 1, 0, 0, 0);
  CompactMatrixAddMultiply('C', 'N', &Cseed, &Tmp, &M1, 1, 0, 0, 0);
  CompactMatrixAddMultiply('C', 'N', &Tmp, &Tmp, &M2, 1, 0, 0, 0);
  /*  the extremes of the diagonal blocks and the norm of the couplings.  The
   *  coupling out of space S is sum_C sum_ll' Wtot_ll' Z_l^dag Z_l'; it is
   *  added to M2 for S = P                                                     */
  {
    double *val = (double*)malloc(Max(North, 1) * sizeof(double));
    CompactMatrixType fun;
    if (val == NULL) { ierr = 1; goto cleanup; }
    if (HermitianCompactMatrixEigensystemRowVectors(&HOrel, &fun, val)) { free(val); ierr = 1; goto cleanup; }
    FreeCompactMatrix(&fun);
    for (lev = 0; lev <= p; lev++) { dlo[lev] = DBL_MAX; dhi[lev] = -DBL_MAX; bnorm[lev] = 0; }
    dlo[0] = val[0];
    dhi[0] = val[North - 1];
    free(val);
  }
  for (s = 0; s < NSpace; s++)
  {
    uint32_t DimS = Space[s].Dim;
    CompactMatrixType Cpl, ZW;
    double *val;
    CompactMatrixSetUnallocated(&Cpl);
    Cpl.Ni = DimS; Cpl.Nj = DimS; Cpl.QComplex = 1;
    if (InitCompactMatrix(&Cpl)) { ierr = 1; goto cleanup; }
    if (s > 0)
    {
      double elo = 0, ehi = 0;
      for (k = 0; k < Space[s].Level; k++) { elo += emin[Space[s].Entry[k]]; ehi += emax[Space[s].Entry[k]]; }
      for (a = 0; a < DimS; a++)
      {
        dlo[Space[s].Level] = Min(dlo[Space[s].Level], Space[s].E[a] + elo - Space[s].Level * dw);
        dhi[Space[s].Level] = Max(dhi[Space[s].Level], Space[s].E[a] + ehi + Space[s].Level * dw);
      }
    }
    for (uint32_t sc = 1; sc < NSpace; sc++)
    {
      BathExcSpaceType *C = &(Space[sc]);
      if (C->Level != Space[s].Level + 1) continue;
      for (c = 0; c < C->NCoupling; c++)
      {
        uint32_t bi;
        if (C->CouplingParent[c] != s) continue;
        i = C->CouplingEntry[c];
        bi = Sigma[i].BlockSize;
        CompactMatrixSetUnallocated(&ZW);
        ZW.Ni = C->Dim; ZW.Nj = DimS; ZW.QComplex = 1;
        if (InitCompactMatrix(&ZW)) { FreeCompactMatrix(&Cpl); ierr = 1; goto cleanup; }
        for (l = 0; l < bi; l++)
        {
          /*  ZW = sum_l' Wtot_ll' Z_l'                                          */
          SetCompactMatrixToZero(&ZW);
          for (lp = 0; lp < bi; lp++)
          {
            Complex16 wt = {0, 0};
            for (j = 0; j < (uint32_t)(jH - jL + 1); j++) { wt.Re += Wbin[i][j][l * bi + lp].Re; wt.Im += Wbin[i][j][l * bi + lp].Im; }
            for (a = 0; a < C->Dim * DimS; a++)
            {
              ZW.V[a].Re += wt.Re * C->Z[c][lp].V[a].Re - wt.Im * C->Z[c][lp].V[a].Im;
              ZW.V[a].Im += wt.Re * C->Z[c][lp].V[a].Im + wt.Im * C->Z[c][lp].V[a].Re;
            }
          }
          CompactMatrixAddMultiply('C', 'N', &(C->Z[c][l]), &ZW, &Cpl, 1, 0, 1, 0);
        }
        FreeCompactMatrix(&ZW);
      }
    }
    SquareCompactMatrixMakeHermitian(&Cpl);
    val = (double*)malloc(Max(DimS, 1) * sizeof(double));
    {
      CompactMatrixType fun;
      if ((val == NULL) || HermitianCompactMatrixEigensystemRowVectors(&Cpl, &fun, val)) { free(val); FreeCompactMatrix(&Cpl); ierr = 1; goto cleanup; }
      FreeCompactMatrix(&fun);
    }
    bnorm[Space[s].Level] = Max(bnorm[Space[s].Level], sqrt(Max(val[DimS - 1], 0.0)));
    free(val);
    if (s == 0)
    {
      CompactMatrixType CT;
      CompactMatrixSetUnallocated(&CT);
      CT.Ni = North; CT.Nj = m; CT.QComplex = 1;
      if (InitCompactMatrix(&CT)) { FreeCompactMatrix(&Cpl); ierr = 1; goto cleanup; }
      CompactMatrixAddMultiply('N', 'N', &Cpl, &Cseed, &CT, 1, 0, 0, 0);
      CompactMatrixAddMultiply('C', 'N', &Cseed, &CT, &M2, 1, 0, 1, 0);
      FreeCompactMatrix(&CT);
    }
    FreeCompactMatrix(&Cpl);
  }
  /*  block Gershgorin: level m couples to m-1 and m+1 only, so every
   *  eigenvalue lies within the diagonal range of some level widened by the
   *  norms of its two couplings                                              */
  lo = DBL_MAX;
  hi = -DBL_MAX;
  for (lev = 0; lev <= p; lev++)
  {
    double bb = bnorm[lev] + ((lev > 0) ? bnorm[lev - 1] : 0);
    if (dlo[lev] == DBL_MAX) continue;
    lo = Min(lo, dlo[lev] - bb);
    hi = Max(hi, dhi[lev] + bb);
  }
  kOutL = (int64_t)floor(Min(lo, Emin) / dw) - 1;
  kOutH = (int64_t)ceil(Max(hi, Emax) / dw) + 1;
  kL[0] = kOutL - 1;
  kH[0] = kOutH + 1;
  for (lev = 1; lev < p; lev++)
  {
    kL[lev] = kL[lev - 1] - jH;
    kH[lev] = kH[lev - 1] - jL;
  }
  NOut = (size_t)(kOutH - kOutL + 1);
  Rspec = Max(Max(fabs(lo), fabs(hi)), Max(fabs(Emin), fabs(Emax)));
  VerbosityPrintf(0x0001, "CalculateGforHwithSigmaAsOperatorBathExcitationsMesh: spectrum within [%g, %g], %zu mesh points of %g, bins %lld .. %lld\n", lo, hi, NOut, dw, (long long)jL, (long long)jH);

  /*  the y quadrature                                                         */
  ymin = BathExcMeshYmin * dw;
  ymax = BathExcMeshYmaxFactor * Max(Rspec, dw);
  {
    uint32_t npanel = (uint32_t)ceil((log(ymax) - log(ymin)) / BathExcMeshPanelWidth), pn;
    double h = (log(ymax) - log(ymin)) / npanel;
    Ny = 8 * npanel;
    ynode = (double*)malloc(Ny * sizeof(double));
    yweight = (double*)malloc(Ny * sizeof(double));
    if ((ynode == NULL) || (yweight == NULL)) { ierr = 1; goto cleanup; }
    for (pn = 0; pn < npanel; pn++)
    {
      for (q = 0; q < 8; q++)
      {
        double u = log(ymin) + h * pn + h * (GLx[q] + 1) / 2;
        ynode[pn * 8 + q] = exp(u);
        yweight[pn * 8 + q] = GLw[q] * h / 2 * exp(u);    /* dy = y du */
      }
    }
    VerbosityPrintf(0x0001, "CalculateGforHwithSigmaAsOperatorBathExcitationsMesh: %u nodes in y from %g to %g\n", Ny, ymin, ymax);
  }

  /* ---------------------------------------------------------------------- */
  /*  F   the nested evaluation, one y at a time, the y in parallel          */
  /* ---------------------------------------------------------------------- */
  Rsum = (Complex16*)calloc(NOut * m * m, sizeof(Complex16));
  if (Rsum == NULL) { printf("malloc failed in CalculateGforHwithSigmaAsOperatorBathExcitationsMesh\n"); ierr = 1; goto cleanup; }
  {
    BathExcMeshType X;
    int failed = 0;
    X.p = p; X.m = m; X.North = North; X.NSpace = NSpace; X.Space = Space; X.Sigma = Sigma;
    X.HOrel = &HOrel; X.Cseed = &Cseed; X.Wbin = Wbin;
    X.kL = kL; X.kH = kH; X.jL = jL; X.jH = jH; X.kOutL = kOutL; X.NOut = NOut; X.dw = dw;
#pragma omp parallel
    {
      Complex16 *Rloc = (Complex16*)calloc(NOut * m * m, sizeof(Complex16));
      int qq;
      if (Rloc == NULL) failed = 1;
#pragma omp for schedule(dynamic)
      for (qq = 0; qq < (int)Ny; qq++)
      {
        if (failed) continue;
        if (BathExcMeshOneY(&X, ynode[qq], yweight[qq], Rloc)) failed = 1;
      }
#pragma omp critical
      {
        if (Rloc != NULL)
        {
          for (size_t z = 0; z < NOut * m * m; z++) { Rsum[z].Re += Rloc[z].Re; Rsum[z].Im += Rloc[z].Im; }
        }
      }
      free(Rloc);
    }
    if (failed) { ierr = 1; goto cleanup; }
  }

  /* ---------------------------------------------------------------------- */
  /*  G   the analytic tail beyond ymax, the residues, the list of poles     */
  /* ---------------------------------------------------------------------- */
  {
    CompactMatrixType *Mom[3] = {&M0, &M1, &M2};
    CompactMatrixType Rk, fun, Rfull, BR;
    double *val = (double*)malloc(m * sizeof(double)), wmaxall = 0;
    uint32_t Nkeep = 0, np;
    int order, jj;
    if (val == NULL) { ierr = 1; goto cleanup; }
    SquareCompactMatrixMakeHermitian(&M0);
    SquareCompactMatrixMakeHermitian(&M1);
    SquareCompactMatrixMakeHermitian(&M2);
    CompactMatrixSetUnallocated(&Rk);
    CompactMatrixSetUnallocated(&Rfull);
    CompactMatrixSetUnallocated(&BR);
    Rk.Ni = m; Rk.Nj = m; Rk.QComplex = 1;
    if (InitCompactMatrix(&Rk)) { free(val); ierr = 1; goto cleanup; }
    for (size_t ko = 0; ko < NOut; ko++)
    {
      double x = (double)(kOutL + (int64_t)ko) * dw, tail;
      for (n = 0; n < 3; n++)
      {
        tail = 0;
        for (order = 2; order <= 4; order += 2)
        {
          double coef = (order == 2) ? dw * dw : pow(dw, 4) / 12, cc = 1;
          int kk2 = (int)n + 1 + order;
          Complex16 w0, p1, p2, I;
          for (jj = 0; jj < order; jj++) cc *= -(double)(n + 1 + jj);
          /*  int_Y^inf y (x+iy)^-kk2 dy = w0^(2-kk2)/(2-kk2) - x w0^(1-kk2)/(1-kk2) */
          w0.Re = x; w0.Im = ymax;
          p1.Re = 1; p1.Im = 0;
          for (jj = 0; jj < kk2 - 1; jj++)
          {
            /*  p1 = w0^-(jj+1)                                                */
            double den = w0.Re * w0.Re + w0.Im * w0.Im;
            Complex16 t2;
            t2.Re = (p1.Re * w0.Re + p1.Im * w0.Im) / den;
            t2.Im = (p1.Im * w0.Re - p1.Re * w0.Im) / den;
            p1 = t2;
          }
          /*  p1 = w0^(1-kk2),  p2 = w0^(2-kk2) = p1 w0                          */
          p2.Re = p1.Re * w0.Re - p1.Im * w0.Im;
          p2.Im = p1.Re * w0.Im + p1.Im * w0.Re;
          I.Re = p2.Re / (2 - kk2) - x * p1.Re / (1 - kk2);
          I.Im = p2.Im / (2 - kk2) - x * p1.Im / (1 - kk2);
          tail += coef * cc * (-1.0 / M_PI) * I.Im;
        }
        for (a = 0; a < m * m; a++)
        {
          Rsum[ko * m * m + a].Re += tail * Mom[n]->V[a].Re;
          Rsum[ko * m * m + a].Im += tail * Mom[n]->V[a].Im;
        }
      }
      for (a = 0; a < m * m; a++)
      {
        Rk.V[a].Re = -Rsum[ko * m * m + a].Re / dw;
        Rk.V[a].Im = -Rsum[ko * m * m + a].Im / dw;
      }
      SquareCompactMatrixMakeHermitian(&Rk);
      /*  what the quadrature leaves negative is removed                     */
      if (HermitianCompactMatrixEigensystemRowVectors(&Rk, &fun, val)) { free(val); FreeCompactMatrix(&Rk); ierr = 1; goto cleanup; }
      SetCompactMatrixToZero(&Rk);
      for (n = 0; n < m; n++)
      {
        if (val[n] <= 0) { ClippedWeight -= val[n]; continue; }
        for (a = 0; a < m; a++)
        {
          for (b = 0; b < m; b++)
          {
            Complex16 fa = fun.V[n * m + a], fb = fun.V[n * m + b];
            Rk.V[a * m + b].Re += val[n] * (fa.Re * fb.Re + fa.Im * fb.Im);
            Rk.V[a * m + b].Im += val[n] * (fa.Im * fb.Re - fa.Re * fb.Im);
          }
        }
      }
      FreeCompactMatrix(&fun);
      memcpy(&(Rsum[ko * m * m]), Rk.V, m * m * sizeof(Complex16));
      {
        double tr = 0;
        for (a = 0; a < m; a++) tr += Rk.V[a * m + a].Re;
        wmaxall = Max(wmaxall, tr);
      }
    }
    free(val);
    FreeCompactMatrix(&Rk);
    for (size_t ko = 0; ko < NOut; ko++)
    {
      double tr = 0;
      for (a = 0; a < m; a++) tr += Rsum[ko * m * m + a * m + a].Re;
      if (tr > BathExcMeshKeepWeight * wmaxall) Nkeep++;
    }
    {
      double tr0 = 0;
      for (a = 0; a < m; a++) tr0 += M0.V[a * m + a].Re;
      if (ClippedWeight > 1e-6 * tr0)
      {
        printf("CalculateGforHwithSigmaAsOperatorBathExcitationsMesh: %g negative weight (of %g) left by the quadrature was removed\n", ClippedWeight, tr0);
      }
    }
    /*  the list of poles, residues B0 R B0^dag                                */
    G->N = Max(Nkeep, 1);
    G->BlockSize = NFull;
    G->QComplex = 1;
    if (InitBlockListOfPoles(G)) { ierr = 1; goto cleanup; }
    G->N = Nkeep;
    G->mu = 0;
    Rk = BathExcView(NULL, m, m);
    BR.Ni = NFull; BR.Nj = m; BR.QComplex = 1;
    Rfull.Ni = NFull; Rfull.Nj = NFull; Rfull.QComplex = 1;
    if (InitCompactMatrix(&BR) || InitCompactMatrix(&Rfull)) { ierr = 1; goto cleanup; }
    np = 0;
    for (size_t ko = 0; ko < NOut; ko++)
    {
      double tr = 0;
      Complex16 *dst;
      for (a = 0; a < m; a++) tr += Rsum[ko * m * m + a * m + a].Re;
      if (!(tr > BathExcMeshKeepWeight * wmaxall)) continue;
      dst = &(G->V[(size_t)np * (1 + NFull * NFull)]);
      dst[0].Re = (double)(kOutL + (int64_t)ko) * dw;
      dst[0].Im = 0;
      Rk.V = &(Rsum[ko * m * m]);
      if (B0 == NULL)
      {
        memcpy(&(dst[1]), Rk.V, m * m * sizeof(Complex16));
      }
      else
      {
        CompactMatrixType B0c;
        CompactMatrixSetUnallocated(&B0c);
        B0c.Ni = NFull; B0c.Nj = m; B0c.QComplex = 1;
        InitCompactMatrix(&B0c);
        for (a = 0; a < NFull * m; a++)
        {
          if (B0->QComplex) B0c.V[a] = B0->V[a];
          else { B0c.V[a].Re = B0->VRe[a]; B0c.V[a].Im = 0; }
        }
        CompactMatrixAddMultiply('N', 'N', &B0c, &Rk, &BR, 1, 0, 0, 0);
        CompactMatrixAddMultiply('N', 'C', &BR, &B0c, &Rfull, 1, 0, 0, 0);
        memcpy(&(dst[1]), Rfull.V, NFull * NFull * sizeof(Complex16));
        FreeCompactMatrix(&B0c);
      }
      np++;
    }
    FreeCompactMatrix(&BR);
    FreeCompactMatrix(&Rfull);
    strcpy(G->Name, "G on a mesh");
    VerbosityPrintf(0x0001, "CalculateGforHwithSigmaAsOperatorBathExcitationsMesh: %u poles on the mesh\n", Nkeep);
  }

cleanup:
  FreeCompactMatrix(&HOrel);
  FreeCompactMatrix(&Tmp);
  FreeCompactMatrix(&M0);
  FreeCompactMatrix(&M1);
  FreeCompactMatrix(&M2);
  if (Wbin != NULL)
  {
    for (i = 0; i < NSigma; i++)
    {
      if (Wbin[i] == NULL) continue;
      for (j = 0; j < (uint32_t)(jH - jL + 1); j++) free(Wbin[i][j]);
      free(Wbin[i]);
    }
    free(Wbin);
  }
  free(Rsum);
  free(ynode);
  free(yweight);
  free(emin);
  free(emax);
  BathExcFrontEndFree(&Front, Sigma);
  return ierr;
}
