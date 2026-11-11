/*
 * AtomHF.cpp
 *
 * The configuration-average Dirac-Hartree-Fock atom (DFT stage A4, design section 8), and
 * its nonrelativistic limit, Cowan's average-of-configuration Hartree-Fock.
 *
 * The energy is the average over all determinants of the (nonrelativistic) configuration,
 *   E_av = sum_a qbar_a I_a + sum_a N_aa [F^0(aa) - g_a/(g_a-1) sum_{k>0} Gamma^k_aa F^k(aa)]
 *          + sum_{a<b} N_ab [F^0(ab) - sum_k Gamma^k_ab G^k(ab)],
 * over the subshells a (Dirac: n kappa; Schroedinger: n l), with qbar_a = q_s g_a/G_s, the
 * pair numbers N of design section 8.1 and
 *   Gamma^k_ab = (j_a k j_b; 1/2 0 -1/2)^2   (Dirac),   1/2 (l_a k l_b; 0 0 0)^2   (Cowan),
 * both 0 unless l_a + l_b + k is even. The radial integrals
 *   R^k(ab;cd) = int rho_ac Y^k[rho_bd]/r dr,   rho_ab = P_a P_b + Q_a Q_b,
 * use the Galerkin solution of the Poisson equation (BSplinePoissonSolve). Galerkin makes
 * R^k(ab;cd) = (2k+1) b_ac^T A^-1 b_bd, b_i = int B_i rho/r, exactly symmetric, so the
 * discrete energy is a quadratic form whose exact derivative is the Fock operator:
 *   dE_av/dc_a = 2 qbar_a (F_a c_a - eps_a S c_a),
 *   F_a = H0 + M[V_N] + sum_b C_ab (M[Y^0_bb/r] - sum_k Gamma^k_ab K^k[b])   (HFCoupling),
 *   K^k[b]_alpha,beta = int rho_alpha,b Y^k[rho_beta,b]/r dr = (2k+1) B^T A^-1 B,
 * B_i,alpha = int B_i rho_alpha,b / r over the basis functions alpha of the group of a; B is
 * banded (the DKB functions are local), so one banded solve with all alpha as right-hand
 * sides gives the exchange matrix.
 *
 * All routines return 0 on success and 1 on error and do not print.
 */

#include "AtomHF.h"

/* Gamma^k_ab of the energy functional: kappaA = kappaB = 0 for Cowan's nonrelativistic
 * shells (1/2 (l_a k l_b; 0 0 0)^2), otherwise (j_a k j_b; 1/2 0 -1/2)^2 with j = |kappa| - 1/2;
 * 0 unless l_a + l_b + k is even (l of the large component; the small components have the
 * same parity). */
double HFGamma(int la, int kappaA, int lb, int kappaB, int k)
{
  double t;
  if ((la + lb + k) % 2 != 0) return 0.;
  if (kappaA == 0 || kappaB == 0)
  {
    t = ThreeJSymbol(la, 0, k, 0, lb, 0);
    return 0.5*t*t;
  }
  t = ThreeJSymbolHalf(2*abs(kappaA) - 1, 1, 2*k, 0, 2*abs(kappaB) - 1, -1);
  return t*t;
}

/* The subshells, pair numbers and angular factors of the atom's configuration, and the
 * Poisson factorizations for k = 0 .. 2 lmax + 1. Dirac: one subshell per orbital of the
 * atom. Schroedinger: one per orbital (Cowan) or, with jjAngular, two (j = l -+ 1/2) with
 * the same radial function and the relativistic angular factors: then E_av must equal
 * Cowan's (design section 8.6, check 5); no Fock matrices in that mode.
 * The ratios N_ab/qbar_a are taken in closed form (no division by qbar):
 *   2N_aa/qbar_a = (q_s-1)(g_a-1)/(G_s-1),
 *   N_ab/qbar_a = (q_s-1) g_b/(G_s-1) (same shell),   q_t g_b/G_t = qbar_b (other shells). */
int HFSetup(DiracHFType *hf, DiracAtomType *atom, int jjAngular, uint32_t nShell, const int *shellN, const int *shellL, const double *shellQ)
{
  RadialBasisType *b = atom->basis;
  uint32_t o, s, a, c, k, nMax, lMax = 0, ld = b->degree + 1;
  int t, nt;
  memset(hf, 0, sizeof(DiracHFType));
  hf->atom = atom;
  hf->nShell = nShell;
  hf->shellQ = (double *) malloc(nShell*sizeof(double));
  hf->shellG = (double *) malloc(nShell*sizeof(double));
  for (s = 0; s < nShell; s++)
  {
    hf->shellQ[s] = shellQ[s];
    hf->shellG[s] = 4.*shellL[s] + 2.;
    if ((uint32_t) shellL[s] > lMax) lMax = (uint32_t) shellL[s];
  }
  nMax = 2*atom->nOrb;
  hf->sub = (HFSubshellType *) calloc(nMax, sizeof(HFSubshellType));
  for (o = 0; o < atom->nOrb; o++)
  {
    DiracOrbitalType *orb = atom->orb + o;
    for (s = 0; s < nShell; s++) if (shellN[s] == orb->n && shellL[s] == orb->l) break;
    if (s == nShell) { HFFree(hf); return 1; }
    nt = (!atom->relativistic && jjAngular && orb->l > 0) ? 2 : 1;
    for (t = 0; t < nt; t++)
    {
      HFSubshellType *sub = hf->sub + hf->nSub;
      sub->n = orb->n;
      sub->l = orb->l;
      sub->shell = s;
      sub->coef = orb->coef;
      sub->orbital = o;
      if (atom->relativistic) sub->kappa = orb->kappa;
      else if (jjAngular) sub->kappa = (t == 0) ? -(orb->l + 1) : orb->l;
      else sub->kappa = 0;
      sub->g = (sub->kappa == 0) ? 4.*orb->l + 2. : 2.*abs(sub->kappa);
      sub->qbar = shellQ[s]*sub->g/hf->shellG[s];
      sub->closed = (shellQ[s] == hf->shellG[s]);
      hf->nSub++;
    }
  }
  hf->kMax = 2*lMax + 1;
  hf->Gam = (double *) calloc(hf->nSub*hf->nSub*(hf->kMax + 1), sizeof(double));
  hf->Rself = (double *) calloc(hf->nSub, sizeof(double));
  hf->Rpair = (double *) calloc(hf->nSub*hf->nSub, sizeof(double));
  hf->Npair = (double *) calloc(hf->nSub*hf->nSub, sizeof(double));
  for (a = 0; a < hf->nSub; a++)
    for (c = 0; c < hf->nSub; c++)
    {
      HFSubshellType *A = hf->sub + a, *C = hf->sub + c;
      double qs = hf->shellQ[A->shell], Gs = hf->shellG[A->shell], qt = hf->shellQ[C->shell], Gt = hf->shellG[C->shell];
      for (k = 0; k <= hf->kMax; k++) hf->Gam[(a*hf->nSub + c)*(hf->kMax + 1) + k] = HFGamma(A->l, A->kappa, C->l, C->kappa, (int) k);
      if (a == c)
      {
        hf->Npair[a*hf->nSub + a] = qs*(qs - 1.)/2.*A->g*(A->g - 1.)/(Gs*(Gs - 1.));
        hf->Rself[a] = (qs - 1.)*(A->g - 1.)/(Gs - 1.);
      }
      else if (A->shell == C->shell)
      {
        hf->Npair[a*hf->nSub + c] = qs*(qs - 1.)*A->g*C->g/(Gs*(Gs - 1.));
        hf->Rpair[a*hf->nSub + c] = (qs - 1.)*C->g/(Gs - 1.);
      }
      else
      {
        hf->Npair[a*hf->nSub + c] = qs*qt*A->g*C->g/(Gs*Gt);
        hf->Rpair[a*hf->nSub + c] = qt*C->g/Gt;
      }
    }
  // Poisson factorizations, k = 0 .. kMax
  hf->nPois = b->nSpline - 1;
  hf->band = (double **) calloc(hf->kMax + 1, sizeof(double *));
  for (k = 0; k <= hf->kMax; k++)
  {
    hf->band[k] = (double *) malloc(ld*hf->nPois*sizeof(double));
    if (BSplinePoissonSetup(b->knot, b->nKnot, b->degree, b->rq, b->wq, b->nPoint, k, hf->band[k], ld)) { HFFree(hf); return 1; }
  }
  hf->P = (double *) calloc(hf->nSub*b->nPoint, sizeof(double));
  hf->Q = (double *) calloc(hf->nSub*b->nPoint, sizeof(double));
  return 0;
}

void HFFree(DiracHFType *hf)
{
  uint32_t k;
  if (hf->band != NULL) for (k = 0; k <= hf->kMax; k++) free(hf->band[k]);
  free(hf->band);
  free(hf->sub);
  free(hf->shellQ);
  free(hf->shellG);
  free(hf->Gam);
  free(hf->Rself);
  free(hf->Rpair);
  free(hf->Npair);
  free(hf->P);
  free(hf->Q);
  memset(hf, 0, sizeof(DiracHFType));
}

// P and Q of every subshell on the quadrature points
void HFOrbitalsOnPoints(DiracHFType *hf)
{
  RadialBasisType *b = hf->atom->basis;
  uint32_t a, nP = b->nPoint;
  double *dP = (double *) malloc(nP*sizeof(double)), *dQ = (double *) malloc(nP*sizeof(double));
  for (a = 0; a < hf->nSub; a++)
  {
    DiracOrbitalType *orb = hf->atom->orb + hf->sub[a].orbital;
    AtomOrbitalOnPoints(b, hf->atom->relativistic, orb->kappa, hf->atom->c, hf->sub[a].coef, hf->P + a*nP, dP, hf->Q + a*nP, dQ);
  }
  free(dP);
  free(dQ);
}

/* Y^k[rho](r)/r on the quadrature points of basis b, rho the radial pair density on them,
 * band the factorized Poisson matrix of BSplinePoissonSetup for this k */
int RadialYkOnPoints(const RadialBasisType *b, uint32_t k, const double *band, const double *rho, double *Yr)
{
  uint32_t q;
  double *Y = (double *) malloc(b->nSpline*sizeof(double));
  if (BSplinePoissonSolve(b->knot, b->nKnot, b->degree, b->rq, b->wq, b->nPoint, k, band, b->degree + 1, rho, Y)) { free(Y); return 1; }
  for (q = 0; q < b->nPoint; q++) Yr[q] = BSplineDerivative(b->rq[q], b->knot, b->nKnot, b->degree, Y, 0)/b->rq[q];
  free(Y);
  return 0;
}

// Y^k[rho](r)/r on the quadrature points, rho the radial pair density on them
int HFYk(DiracHFType *hf, uint32_t k, const double *rho, double *Yr)
{
  return RadialYkOnPoints(hf->atom->basis, k, hf->band[k], rho, Yr);
}

// int rho Yr dr on the quadrature points, compensated summation
double HFIntegral(const RadialBasisType *b, const double *rho, const double *Yr)
{
  DFTReal sum = 0., comp = 0.;
  uint32_t q;
  for (q = 0; q < b->nPoint; q++) NeumaierAdd(&sum, &comp, b->wq[q]*rho[q]*Yr[q]);
  return sum + comp;
}

/* E_av of the subshells' radial functions (HFOrbitalsOnPoints first), which must be
 * normalized: Etot, and its parts Ekin = sum qbar <a|H0|a>, Eenuc = sum qbar <a|V_N|a>,
 * Edirect (the F^0 terms) and Eexchange (the rest). */
int HFEnergy(DiracHFType *hf, double *Etot, double *Ekin, double *Eenuc, double *Edirect, double *Eexchange)
{
  DiracAtomType *atom = hf->atom;
  RadialBasisType *b = atom->basis;
  uint32_t a, c, k, q, nP = b->nPoint, K1 = hf->kMax + 1;
  double *rho = (double *) malloc(nP*sizeof(double)), *Yr = (double *) malloc(nP*sizeof(double));
  double *Y0 = (double *) malloc(hf->nSub*nP*sizeof(double)), I, kin, Fk, gam, sumD = 0., sumX = 0.;
  int error = 0;
  *Ekin = 0.;
  *Eenuc = 0.;
  // one-electron part; Y^0 of every subshell
  for (a = 0; a < hf->nSub && !error; a++)
  {
    DiracGroupType *gr = atom->group + atom->orb[hf->sub[a].orbital].group;
    I = AtomBilinear(hf->sub[a].coef, gr->Hfix, hf->sub[a].coef, atom->Nf);
    kin = AtomBilinear(hf->sub[a].coef, gr->H0, hf->sub[a].coef, atom->Nf);
    *Ekin += hf->sub[a].qbar*kin;
    *Eenuc += hf->sub[a].qbar*(I - kin);
    for (q = 0; q < nP; q++) rho[q] = hf->P[a*nP + q]*hf->P[a*nP + q] + hf->Q[a*nP + q]*hf->Q[a*nP + q];
    error |= HFYk(hf, 0, rho, Y0 + a*nP);
  }
  for (a = 0; a < hf->nSub && !error; a++)
  {
    HFSubshellType *A = hf->sub + a;
    for (q = 0; q < nP; q++) rho[q] = hf->P[a*nP + q]*hf->P[a*nP + q] + hf->Q[a*nP + q]*hf->Q[a*nP + q];
    // aa: N_aa [F^0 - g/(g-1) sum_k>0 Gamma^k F^k]
    sumD += hf->Npair[a*hf->nSub + a]*HFIntegral(b, rho, Y0 + a*nP);
    for (k = 1; k <= hf->kMax; k++)
    {
      gam = hf->Gam[(a*hf->nSub + a)*K1 + k];
      if (gam == 0.) continue;
      error |= HFYk(hf, k, rho, Yr);
      sumX -= hf->Npair[a*hf->nSub + a]*A->g/(A->g - 1.)*gam*HFIntegral(b, rho, Yr);
    }
    for (c = a + 1; c < hf->nSub; c++)
    {
      double Nac = hf->Npair[a*hf->nSub + c], *rac;
      sumD += Nac*HFIntegral(b, rho, Y0 + c*nP);
      rac = (double *) malloc(nP*sizeof(double));
      for (q = 0; q < nP; q++) rac[q] = hf->P[a*nP + q]*hf->P[c*nP + q] + hf->Q[a*nP + q]*hf->Q[c*nP + q];
      for (k = 0; k <= hf->kMax; k++)
      {
        gam = hf->Gam[(a*hf->nSub + c)*K1 + k];
        if (gam == 0.) continue;
        error |= HFYk(hf, k, rac, Yr);
        Fk = HFIntegral(b, rac, Yr);
        sumX -= Nac*gam*Fk;
      }
      free(rac);
    }
  }
  *Edirect = sumD;
  *Eexchange = sumX;
  *Etot = *Ekin + *Eenuc + sumD + sumX;
  free(rho); free(Yr); free(Y0);
  return error;
}

/* K^k[b] in the basis of symmetry group `group`: K_alpha,beta = int rho_alpha,b Y^k[rho_beta,b]/r
 * = (2k+1) B^T A^-1 B with B_i,alpha = int B_i (P_alpha P_b + Q_alpha Q_b)/r dr, i over the
 * Poisson functions B_1 .. B_{N-1} (Nf x Nf, dense, column major). */
int HFExchangeMatrix(DiracHFType *hf, uint32_t group, uint32_t bsub, uint32_t k, double *K)
{
  DiracAtomType *atom = hf->atom;
  RadialBasisType *bs = atom->basis;
  uint32_t Nf = atom->Nf, nP = bs->nPoint, nPois = hf->nPois, d = bs->degree, q, j, m, s, mu, i, idx[2*(BSplineMaxDegree+1)], kd = d, ld = d + 1;
  double P[2*(BSplineMaxDegree+1)], dP[2*(BSplineMaxDegree+1)], Q[2*(BSplineMaxDegree+1)], dQ[2*(BSplineMaxDegree+1)];
  double val[BSplineMaxDegree+1], *Bm, *X, w, Pb, Qb, scale = 2.*k + 1., zero = 0.;
  int kappa = atom->group[group].kappa, info = 0;
  char uplo = 'L', tr = 'T', nt = 'N';
  Bm = (double *) calloc(nPois*Nf, sizeof(double));
  for (q = 0; q < nP; q++)
  {
    m = AtomLocalFunctions(bs, atom->relativistic, kappa, atom->c, q, idx, P, dP, Q, dQ);
    if (m == 0) continue;
    mu = BSplineFindInterval(bs->rq[q], bs->knot, bs->nKnot);
    if (mu >= bs->nKnot) continue;
    BSplineValues(bs->rq[q], mu, bs->knot, bs->nKnot, d, 0, val);
    Pb = hf->P[bsub*nP + q];
    Qb = hf->Q[bsub*nP + q];
    w = bs->wq[q]/bs->rq[q];
    for (s = 0; s <= d; s++)
    {
      int64_t sp = (int64_t) mu - d + s;            // B-spline index; Poisson index sp-1
      if (sp < 1 || sp >= (int64_t) bs->nSpline) continue;
      i = (uint32_t) (sp - 1);
      for (j = 0; j < m; j++) Bm[i + idx[j]*nPois] += w*val[s]*(P[j]*Pb + Q[j]*Qb);
    }
  }
  X = (double *) malloc(nPois*Nf*sizeof(double));
  for (q = 0; q < nPois*Nf; q++) X[q] = Bm[q];
  dpbtrs(&uplo, &nPois, &kd, &Nf, hf->band[k], &ld, X, &nPois, &info);
  if (info == 0) dgemm(&tr, &nt, &Nf, &Nf, &nPois, &scale, Bm, &nPois, X, &nPois, &zero, K, &Nf);
  // symmetric to rounding; make it exactly so
  for (j = 0; j < Nf && info == 0; j++)
    for (i = j + 1; i < Nf; i++) K[i + j*Nf] = K[j + i*Nf] = 0.5*(K[i + j*Nf] + K[j + i*Nf]);
  free(Bm);
  free(X);
  return (info != 0);
}

/* The coefficient of J_b - sum_k Gamma^k_ab K^k[b] in the Fock operator of subshell a:
 * N_ab/qbar_a for b != a, and 2N_aa/qbar_a g_a/(g_a-1) for b = a. The self term
 * N_aa [F^0(aa) - g/(g-1) sum_{k>0} Gamma^k F^k(aa)] equals, with Gamma^0_aa = 1/g,
 * N_aa g/(g-1) [F^0(aa) - sum_{k>=0} Gamma^k_aa F^k(aa)]; written so, its derivative is
 * the same operator form as the pair terms (J - K, with K^k[a] c_a = M[Y^k_aa/r] c_a).
 * Then every closed subshell of a symmetry group has the same Fock operator, the
 * closed-shell (Roothaan) one: the coefficient is g_b for every b. */
double HFCoupling(const DiracHFType *hf, uint32_t a, uint32_t b)
{
  if (a != b) return hf->Rpair[a*hf->nSub + b];
  return hf->Rself[a]*hf->sub[a].g/(hf->sub[a].g - 1.);
}

/* The Fock matrices F[a] (Nf x Nf each, allocated by the caller) of every subshell, from
 * the radial functions on the quadrature points (HFOrbitalsOnPoints first):
 *   F_a = H0 + M[V_N] + sum_b HFCoupling(a, b) (M[Y^0_bb/r] - sum_k Gamma^k_ab K^k[b]),
 * so that dE_av/dc_a = 2 qbar_a F_a c_a for unnormalized coefficients (design check 6).
 * Not for the jjAngular mode (its subshells share radial functions). */
int HFFockMatrices(DiracHFType *hf, double **F)
{
  DiracAtomType *atom = hf->atom;
  RadialBasisType *b = atom->basis;
  uint32_t a, c, k, q, i, gr, nP = b->nPoint, Nf = atom->Nf, K1 = hf->kMax + 1;
  double *rho = (double *) malloc(nP*sizeof(double)), *v = (double *) malloc(nP*sizeof(double));
  double *Y0 = (double *) malloc(hf->nSub*nP*sizeof(double)), *M = (double *) malloc(Nf*Nf*sizeof(double)), coefK;
  int error = 0, need;
  for (a = 0; a < hf->nSub && !error; a++)
  {
    for (q = 0; q < nP; q++) rho[q] = hf->P[a*nP + q]*hf->P[a*nP + q] + hf->Q[a*nP + q]*hf->Q[a*nP + q];
    error |= HFYk(hf, 0, rho, Y0 + a*nP);
  }
  // local part
  for (a = 0; a < hf->nSub && !error; a++)
  {
    HFSubshellType *A = hf->sub + a;
    DiracGroupType *g = atom->group + atom->orb[A->orbital].group;
    for (q = 0; q < nP; q++) v[q] = 0.;
    for (c = 0; c < hf->nSub; c++)
      for (q = 0; q < nP; q++) v[q] += HFCoupling(hf, a, c)*Y0[c*nP + q];
    AtomPotentialMatrix(b, atom->relativistic, atom->orb[A->orbital].kappa, atom->c, v, NULL, M);
    for (i = 0; i < Nf*Nf; i++) F[a][i] = g->Hfix[i] + M[i];
  }
  // exchange: K^k[c] once per group, partner c and k, added to every subshell of the group
  for (gr = 0; gr < atom->nGroup && !error; gr++)
    for (c = 0; c < hf->nSub && !error; c++)
      for (k = 0; k <= hf->kMax && !error; k++)
      {
        need = 0;
        for (a = 0; a < hf->nSub; a++)
          if (atom->orb[hf->sub[a].orbital].group == gr && HFCoupling(hf, a, c)*hf->Gam[(a*hf->nSub + c)*K1 + k] != 0.) need = 1;
        if (!need) continue;
        error |= HFExchangeMatrix(hf, gr, c, k, M);
        for (a = 0; a < hf->nSub; a++)
        {
          if (atom->orb[hf->sub[a].orbital].group != gr) continue;
          coefK = -HFCoupling(hf, a, c)*hf->Gam[(a*hf->nSub + c)*K1 + k];
          if (coefK == 0.) continue;
          for (i = 0; i < Nf*Nf; i++) F[a][i] += coefK*M[i];
        }
      }
  free(rho); free(v); free(Y0); free(M);
  return error;
}

/* The gradient of E_av on the orbitals, orthonormal per symmetry group (F from
 * HFFockMatrices of the same orbitals); sets eps_a = c_a^T F_a c_a.
 *  - Towards the unoccupied space: G_a = 2 qbar_a r_a with
 *      r_a = F_a c_a - sum_b S c_b (c_b^T F_a c_a)   (b over the subshells of a's group),
 *    *gradient = max_a |G_a| in the dual norm sqrt(G^T S^-1 G) (independent of the basis).
 *  - Rotations between two occupied subshells a, b of one group (c_a cos t + c_b sin t,
 *    c_b cos t - c_a sin t): dE/dt = g_ab = 2 (qbar_a c_b^T F_a c_a - qbar_b c_a^T F_b c_b),
 *    *pairGradient = max |g_ab| (zero by symmetry for two closed subshells).
 * Both vanish at the stationary point. Returns 1 if LAPACK fails. */
int HFGradient(DiracHFType *hf, double **F, double *gradient, double *pairGradient)
{
  DiracAtomType *atom = hf->atom;
  uint32_t Nf = atom->Nf, a, b, i, j, N, nrhs = 1;
  double *Fc = (double *) malloc(Nf*sizeof(double)), *Sc = (double *) malloc(Nf*sizeof(double)), *r = (double *) malloc(Nf*sizeof(double));
  double *Sa = (double *) malloc(Nf*Nf*sizeof(double)), *ra = (double *) malloc(Nf*sizeof(double)), lam, norm, gab;
  int info = 0;
  char uplo = 'L';
  *gradient = 0.;
  *pairGradient = 0.;
  for (a = 0; a < hf->nSub && info == 0; a++)
  {
    HFSubshellType *A = hf->sub + a;
    uint32_t gr = atom->orb[A->orbital].group;
    DiracGroupType *g = atom->group + gr;
    for (i = 0; i < Nf; i++) { Fc[i] = 0.; for (j = 0; j < Nf; j++) Fc[i] += F[a][i + j*Nf]*A->coef[j]; }
    A->energy = 0.;
    for (i = 0; i < Nf; i++) { A->energy += A->coef[i]*Fc[i]; r[i] = Fc[i]; }
    for (b = 0; b < hf->nSub; b++)
    {
      HFSubshellType *B = hf->sub + b;
      if (atom->orb[B->orbital].group != gr) continue;
      for (i = 0; i < Nf; i++) { Sc[i] = 0.; for (j = 0; j < Nf; j++) Sc[i] += g->S[i + j*Nf]*B->coef[j]; }
      lam = 0.;
      for (i = 0; i < Nf; i++) lam += B->coef[i]*Fc[i];
      for (i = 0; i < Nf; i++) r[i] -= lam*Sc[i];
      if (b > a && !(A->closed && B->closed))
      {
        gab = 2.*(A->qbar*lam - B->qbar*AtomBilinear(A->coef, F[b], B->coef, Nf));
        if (fabs(gab) > *pairGradient) *pairGradient = fabs(gab);
      }
    }
    // dual norm on the active functions: r^T S^-1 r
    N = g->N;
    for (i = 0; i < N; i++)
    {
      ra[i] = r[g->active[i]];
      for (j = 0; j < N; j++) Sa[i + j*N] = g->S[g->active[i] + g->active[j]*Nf];
    }
    dpotrf(&uplo, &N, Sa, &N, &info);
    if (info != 0) break;
    for (i = 0; i < N; i++) Sc[i] = ra[i];
    dpotrs(&uplo, &N, &nrhs, Sa, &N, Sc, &N, &info);
    norm = 0.;
    for (i = 0; i < N; i++) norm += ra[i]*Sc[i];
    norm = 2.*A->qbar*sqrt(fabs(norm));
    if (norm > *gradient) *gradient = norm;
  }
  free(Fc); free(Sc); free(r); free(Sa); free(ra);
  return (info != 0);
}

/* Orthonormalize the orbitals of every symmetry group (Lowdin, two passes: the smallest
 * change of the set). Returns 1 if they are linearly dependent or LAPACK fails. */
int HFOrthonormalize(DiracHFType *hf)
{
  DiracAtomType *atom = hf->atom;
  uint32_t Nf = atom->Nf, gr, a, m, i;
  double *C = (double *) malloc(hf->nSub*Nf*sizeof(double)), minEig;
  int error = 0;
  for (gr = 0; gr < atom->nGroup && !error; gr++)
  {
    m = 0;
    for (a = 0; a < hf->nSub; a++)
      if (atom->orb[hf->sub[a].orbital].group == gr) { for (i = 0; i < Nf; i++) C[i + m*Nf] = hf->sub[a].coef[i]; m++; }
    if (m == 0) continue;
    error = DiracOrthonormalize(atom->group[gr].S, Nf, C, m, DiracLowdin, 0., &minEig);
    m = 0;
    for (a = 0; a < hf->nSub && !error; a++)
      if (atom->orb[hf->sub[a].orbital].group == gr) { for (i = 0; i < Nf; i++) hf->sub[a].coef[i] = C[i + m*Nf]; m++; }
  }
  free(C);
  return error;
}

/* New orbitals of one class of subshells of symmetry group gr with the same Fock matrix
 * Fa: the closed subshells of the group together (their Fock matrices are equal, see
 * HFCoupling), or one open subshell. The other occupied orbitals of the group (Co, the
 * input orbitals, S-orthonormal) are projected out exactly,
 *   F'' = Q^T Fa Q + sigma S Co Co^T S,   Q = 1 - Co Co^T S:
 * F'' is block diagonal in span(Co) + its S-orthogonal complement, with eigenvalue sigma
 * on span(Co) and Q^T Fa Q on the complement, so its eigenvectors below sigma are
 * orthogonal to Co and satisfy Fa c - eps S c in span(S Co): the gradient towards the
 * unoccupied space vanishes at the fixed point. sigma only has to lie above the occupied
 * eigenvalues. Each subshell of the class takes the root of largest overlap with its input
 * (distinct roots), with the sign of the input. out[a] receives the new coefficients. */
int HFClassOrbitals(DiracHFType *hf, uint32_t gr, const double *Fa, const int *inClass, double sigma, int rayleigh, double **out)
{
  DiracAtomType *atom = hf->atom;
  DiracGroupType *g = atom->group + gr;
  uint32_t Nf = atom->Nf, a, b, i, j, mo = 0, first, nRoot, best;
  double *Co = (double *) malloc(hf->nSub*Nf*sizeof(double)), *A, *B, *Mo, *H, *E, *C, ov, bestOv, sgn;
  int *taken, error;
  for (b = 0; b < hf->nSub; b++)
  {
    if (atom->orb[hf->sub[b].orbital].group != gr || inClass[b]) continue;
    for (i = 0; i < Nf; i++) Co[i + mo*Nf] = hf->sub[b].coef[i];
    mo++;
  }
  A = (double *) calloc(Nf*(mo + 1), sizeof(double));     // S Co
  B = (double *) calloc(Nf*(mo + 1), sizeof(double));     // Fa Co
  Mo = (double *) calloc((mo + 1)*(mo + 1), sizeof(double));
  H = (double *) malloc(Nf*Nf*sizeof(double));
  E = (double *) malloc(Nf*sizeof(double));
  C = (double *) malloc(Nf*Nf*sizeof(double));
  for (b = 0; b < mo; b++)
    for (i = 0; i < Nf; i++)
      for (j = 0; j < Nf; j++)
      {
        A[i + b*Nf] += g->S[i + j*Nf]*Co[j + b*Nf];
        B[i + b*Nf] += Fa[i + j*Nf]*Co[j + b*Nf];
      }
  for (a = 0; a < mo; a++)
    for (b = 0; b < mo; b++)
      for (i = 0; i < Nf; i++) Mo[a + b*mo] += Co[i + a*Nf]*B[i + b*Nf];
  for (j = 0; j < Nf; j++)
    for (i = 0; i < Nf; i++)
    {
      double h = Fa[i + j*Nf];
      for (a = 0; a < mo; a++)
      {
        h += sigma*A[i + a*Nf]*A[j + a*Nf] - A[i + a*Nf]*B[j + a*Nf] - B[i + a*Nf]*A[j + a*Nf];
        for (b = 0; b < mo; b++) h += A[i + a*Nf]*Mo[a + b*mo]*A[j + b*Nf];
      }
      H[i + j*Nf] = h;
    }
  for (j = 0; j < Nf; j++) for (i = j + 1; i < Nf; i++) H[i + j*Nf] = H[j + i*Nf] = 0.5*(H[i + j*Nf] + H[j + i*Nf]);
  error = DiracGeneralizedEigen(H, g->S, Nf, g->active, g->N, rayleigh, E, C, NULL, NULL);
  first = atom->relativistic ? DiracElectronicFirst(E, g->N, atom->c) : 0;
  nRoot = 0;
  while (first + nRoot < g->N && E[first + nRoot] < sigma/2.) nRoot++;
  taken = (int *) calloc(nRoot + 1, sizeof(int));
  for (a = 0; a < hf->nSub && !error; a++)
  {
    if (!inClass[a]) continue;
    best = nRoot;
    bestOv = -1.;
    for (j = 0; j < nRoot; j++)
    {
      if (taken[j]) continue;
      ov = fabs(AtomBilinear(C + (first + j)*Nf, g->S, hf->sub[a].coef, Nf));
      if (ov > bestOv) { bestOv = ov; best = j; }
    }
    if (best == nRoot) { error = 1; break; }
    taken[best] = 1;
    sgn = (AtomBilinear(C + (first + best)*Nf, g->S, hf->sub[a].coef, Nf) < 0.) ? -1. : 1.;
    for (i = 0; i < Nf; i++) out[a][i] = sgn*C[i + (first + best)*Nf];
  }
  free(Co); free(A); free(B); free(Mo); free(H); free(E); free(C); free(taken);
  return error;
}

/* The configuration-average Hartree-Fock SCF (design section 8.4), from the orbitals in
 * hf (e.g. a converged LDA atom), orthonormal per symmetry group. One iteration:
 *  1. Fock matrices and E_av of the input orbitals x (HFFockMatrices, HFEnergy) and the
 *     gradient (HFGradient, reported only);
 *  2. new orbitals per class (HFClassOrbitals: the closed subshells of a group in one
 *     eigenproblem, every open subshell in its own, the others projected out),
 *     orthonormalized (Lowdin);
 *  3. one sweep of Brillouin rotations over the pairs of a group that are not both closed
 *     (their Fock matrices differ, so step 2 does not fix the rotation between them):
 *     Newton t = -g_ab/E'' with the fixed-Fock second derivative
 *       E'' = 2 [qbar_a (c_b^T F_a c_b - c_a^T F_a c_a) + qbar_b (c_a^T F_b c_a - c_b^T F_b c_b)],
 *     skipped if E'' <= 0;
 *  4. converged when |dE| < energyTolerance and the step max_a |R_a|_S < residualTolerance,
 *     R = x_out - x; else Pulay's DIIS on the stacked coefficients (S metric, the last
 *     `history` iterations, x_new = sum_i c_i (x_i + mixing R_i)), orthonormalized.
 * The fixed point is stationary: step 2 makes the gradient towards the unoccupied space
 * vanish, step 3 the occupied-occupied one; so the step R measures both. The gradient
 * itself is no criterion: in the dual norm it has a rounding floor of a few 1e-8 (the
 * steep functions near the nucleus, F ~ 1/r1^2), while the step goes to 1e-13 (Be, Ne,
 * nonrelativistic) with E stable to 1e-14. On return the orbitals are those of the last
 * energy; hf holds the energies, eps_a, and E from 1/2 sum qbar (eps + I), which equals
 * Etot for any orbitals (E is homogeneous of degree 2 and 4 in the one- and two-electron
 * parts): a check of the Fock matrices against the energy, not of convergence.
 * Returns 1 if a step fails; hf->converged tells whether it converged. */
int HFSCF(DiracHFType *hf, uint32_t maxIterations, uint32_t history, double mixing, double energyTolerance,
          double residualTolerance, int rayleigh)
{
  DiracAtomType *atom = hf->atom;
  uint32_t Nf = atom->Nf, nS = hf->nSub, nx = nS*Nf, it, a, b, i, j, k, gr, m, slot, nB, nrhs = 1;
  double **F, **out, *x, *R, *X, *RR, *SR, *Bm, *rhs, Ebefore = 0., sigma, sum, cs, I, Eone;
  double cf, sf, gab, d2, t, ca, cb;
  int *inClass, *ipiv, info, error = 0, any;
  if (history == 0) history = 1;
  F = (double **) malloc(nS*sizeof(double *));
  out = (double **) malloc(nS*sizeof(double *));
  for (a = 0; a < nS; a++) { F[a] = (double *) malloc(Nf*Nf*sizeof(double)); out[a] = (double *) malloc(Nf*sizeof(double)); }
  x = (double *) malloc(nx*sizeof(double));
  R = (double *) malloc(nx*sizeof(double));
  X = (double *) malloc(history*nx*sizeof(double));
  RR = (double *) malloc(history*nx*sizeof(double));
  SR = (double *) malloc(nx*sizeof(double));
  Bm = (double *) malloc((history + 1)*(history + 1)*sizeof(double));
  rhs = (double *) malloc((history + 1)*sizeof(double));
  ipiv = (int *) malloc((history + 1)*sizeof(int));
  inClass = (int *) malloc(nS*sizeof(int));
  sigma = 2.*atom->Z*atom->Z + 10.;
  hf->converged = 0;
  error = HFOrthonormalize(hf);
  for (it = 1; it <= maxIterations && !error; it++)
  {
    // 1. Fock matrices, energy, gradient of the input orbitals
    HFOrbitalsOnPoints(hf);
    if (HFFockMatrices(hf, F) || HFEnergy(hf, &hf->Etot, &hf->Ekin, &hf->Eenuc, &hf->Edirect, &hf->Eexchange)) { error = 1; break; }
    if (HFGradient(hf, F, &hf->gradient, &hf->pairGradient)) { error = 1; break; }
    Eone = 0.;
    for (a = 0; a < nS; a++)
    {
      DiracGroupType *g = atom->group + atom->orb[hf->sub[a].orbital].group;
      I = AtomBilinear(hf->sub[a].coef, g->Hfix, hf->sub[a].coef, Nf);
      Eone += 0.5*hf->sub[a].qbar*(hf->sub[a].energy + I);
    }
    hf->Eidentity = Eone;
    hf->dE = hf->Etot - Ebefore;
    hf->iterations = it;
    Ebefore = hf->Etot;
    for (a = 0; a < nS; a++) for (i = 0; i < Nf; i++) x[a*Nf + i] = hf->sub[a].coef[i];
    // 2. new orbitals per class
    for (gr = 0; gr < atom->nGroup && !error; gr++)
    {
      any = -1;
      for (a = 0; a < nS; a++)
      {
        inClass[a] = (atom->orb[hf->sub[a].orbital].group == gr && hf->sub[a].closed);
        if (inClass[a] && any < 0) any = (int) a;
      }
      if (any >= 0) error |= HFClassOrbitals(hf, gr, F[any], inClass, sigma, rayleigh, out);
      for (a = 0; a < nS && !error; a++)
      {
        if (atom->orb[hf->sub[a].orbital].group != gr || hf->sub[a].closed) continue;
        for (b = 0; b < nS; b++) inClass[b] = (b == a);
        error |= HFClassOrbitals(hf, gr, F[a], inClass, sigma, rayleigh, out);
      }
    }
    if (error) break;
    for (a = 0; a < nS; a++) for (i = 0; i < Nf; i++) hf->sub[a].coef[i] = out[a][i];
    if (HFOrthonormalize(hf)) { error = 1; break; }
    // 3. Brillouin rotations
    for (a = 0; a < nS; a++)
      for (b = a + 1; b < nS; b++)
      {
        HFSubshellType *A = hf->sub + a, *B = hf->sub + b;
        if (atom->orb[A->orbital].group != atom->orb[B->orbital].group || (A->closed && B->closed)) continue;
        gab = 2.*(A->qbar*AtomBilinear(B->coef, F[a], A->coef, Nf) - B->qbar*AtomBilinear(A->coef, F[b], B->coef, Nf));
        d2 = 2.*(A->qbar*(AtomBilinear(B->coef, F[a], B->coef, Nf) - AtomBilinear(A->coef, F[a], A->coef, Nf))
               + B->qbar*(AtomBilinear(A->coef, F[b], A->coef, Nf) - AtomBilinear(B->coef, F[b], B->coef, Nf)));
        if (!(d2 > 0.)) continue;
        t = -gab/d2;
        cf = cos(t);
        sf = sin(t);
        for (i = 0; i < Nf; i++)
        {
          ca = A->coef[i];
          cb = B->coef[i];
          A->coef[i] = cf*ca + sf*cb;
          B->coef[i] = cf*cb - sf*ca;
        }
      }
    // 4. DIIS on the coefficients
    for (a = 0; a < nS; a++) for (i = 0; i < Nf; i++) R[a*Nf + i] = hf->sub[a].coef[i] - x[a*Nf + i];
    hf->residual = 0.;
    for (a = 0; a < nS; a++)
    {
      double ra = sqrt(fabs(AtomBilinear(R + a*Nf, atom->group[atom->orb[hf->sub[a].orbital].group].S, R + a*Nf, Nf)));
      if (ra > hf->residual) hf->residual = ra;
    }
    if ((it > 1 && fabs(hf->dE) < energyTolerance && hf->residual < residualTolerance) || it == maxIterations)
    {
      // keep the orbitals of the energy
      hf->converged = (it < maxIterations || (fabs(hf->dE) < energyTolerance && hf->residual < residualTolerance));
      for (a = 0; a < nS; a++) for (i = 0; i < Nf; i++) hf->sub[a].coef[i] = x[a*Nf + i];
      break;
    }
    slot = (it - 1) % history;
    for (i = 0; i < nx; i++) { X[slot*nx + i] = x[i]; RR[slot*nx + i] = R[i]; }
    m = (it < history) ? it : history;
    nB = m + 1;
    for (j = 0; j < m; j++)
    {
      for (a = 0; a < nS; a++)
      {
        DiracGroupType *g = atom->group + atom->orb[hf->sub[a].orbital].group;
        for (i = 0; i < Nf; i++) { SR[a*Nf + i] = 0.; for (k = 0; k < Nf; k++) SR[a*Nf + i] += g->S[i + k*Nf]*RR[j*nx + a*Nf + k]; }
      }
      for (k = 0; k <= j; k++)
      {
        sum = 0.;
        for (i = 0; i < nx; i++) sum += SR[i]*RR[k*nx + i];
        Bm[j + k*nB] = Bm[k + j*nB] = sum;
      }
    }
    for (j = 0; j < m; j++) { Bm[m + j*nB] = Bm[j + m*nB] = -1.; rhs[j] = 0.; }
    Bm[m + m*nB] = 0.;
    rhs[m] = -1.;
    dgesv(&nB, &nrhs, Bm, &nB, ipiv, rhs, &nB, &info);
    for (i = 0; i < nx; i++)
    {
      if (info != 0) cs = x[i] + mixing*R[i];
      else { cs = 0.; for (j = 0; j < m; j++) cs += rhs[j]*(X[j*nx + i] + mixing*RR[j*nx + i]); }
      hf->sub[i/Nf].coef[i % Nf] = cs;
    }
    if (HFOrthonormalize(hf)) { error = 1; break; }
  }
  for (a = 0; a < nS; a++) { free(F[a]); free(out[a]); }
  free(F); free(out); free(x); free(R); free(X); free(RR); free(SR); free(Bm); free(rhs); free(ipiv); free(inClass);
  return error;
}
