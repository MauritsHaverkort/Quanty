/*
 * AtomSCF.cpp
 *
 * The self-consistent Dirac-Kohn-Sham atom (DFT stage A3, design section 9), and its
 * nonrelativistic (Schroedinger) limit. Spherical, spin unpolarized; the density comes from
 * the mean occupations of the subshells,
 *     n(r) = 1/(4 pi r^2) sum_a qbar_a (P_a^2 + Q_a^2),   qbar_a = q_s g_a/G_s,
 * which is also the j-averaging of the NIST atomic reference data.
 *
 * Every quantity lives on the quadrature points of one radial basis. The Hamiltonian of a
 * symmetry group (kappa, or l without relativity) is
 *     H = H0 + M[V_N] + M[v, g],
 * H0 the kinetic part (V = 0), M[v, g]_ab = int v (P_a P_b + Q_a Q_b) + g (P_a P_b + Q_a Q_b)' dr,
 * with the local potential v = V_H + de/dn - 4 (de/dsigma) n'/r and, for a GGA, the
 * coupling g = 2 (de/dsigma) n' of the weak form of design section 9. H0 and M[V_N] are
 * built once; M[v, g] each iteration, by direct quadrature on the DKB functions.
 * Mixing: Pulay (DIIS) on (v, g) on the quadrature points.
 *
 * All routines return 0 on success and 1 on error and do not print.
 */

#include "AtomSCF.h"

/* The shells of a configuration such as "1s2 2s2 2p6", "[Ne] 3s2 3p6", "[Ar] 3d7.5 4s0.5" or
 * "[Kr] 4d3/2(4) 4d5/2(5)": n, l, kappa and the (possibly fractional) occupation q of each,
 * at most maxShell, appended to shell[*nShell ..]. An nl entry ("3d8") is a nonrelativistic
 * shell (kappa = 0; Cowan's average over its 4l+2 spinors), an n l 2j/2 entry with the
 * occupation in parentheses ("4d5/2(5)") one relativistic subshell (kappa = -(l+1) for
 * j = l+1/2, l for j = l-1/2; the average over its 2j+1 spinors). "[He]", "[Ne]", "[Ar]",
 * "[Kr]", "[Xe]" and "[Rn]" stand for the closed shells of the noble gases. Returns 1 for a
 * syntax error, a shell or subshell given twice (also as nl and as n kappa), n <= l, l > 6,
 * j not l -+ 1/2, q < 0 or q above the degeneracy. */
int AtomParseConfiguration(const char *configuration, uint32_t maxShell, uint32_t *nShell, AtomShellType *shell)
{
  const char *cores[6][2] = {{"[He]", "1s2"}, {"[Ne]", "[He] 2s2 2p6"}, {"[Ar]", "[Ne] 3s2 3p6"}, {"[Kr]", "[Ar] 3d10 4s2 4p6"},
                             {"[Xe]", "[Kr] 4d10 5s2 5p6"}, {"[Rn]", "[Xe] 4f14 5d10 6s2 6p6"}};
  const char *letters = "spdfghi";
  const char *p = configuration, *where;
  char *end;
  uint32_t i, j;
  int n, l, kappa, twoJ, found;
  double q, G;
  while (*p != '\0')
  {
    while (*p == ' ' || *p == '\t' || *p == ',') p++;
    if (*p == '\0') break;
    if (*p == '[')
    {
      found = 0;
      for (i = 0; i < 6 && !found; i++)
        if (strncmp(p, cores[i][0], 4) == 0)
        {
          if (AtomParseConfiguration(cores[i][1], maxShell, nShell, shell)) return 1;
          p += 4;
          found = 1;
        }
      if (!found) return 1;
      continue;
    }
    n = (int) strtol(p, &end, 10);
    if (end == p || n < 1) return 1;
    p = end;
    where = strchr(letters, *p);
    if (where == NULL || *p == '\0') return 1;
    l = (int) (where - letters);
    p++;
    kappa = 0;
    twoJ = (int) strtol(p, &end, 10);
    if (end != p && end[0] == '/' && end[1] == '2')
    {
      // n l 2j/2 (q)
      if (twoJ == 2*l + 1) kappa = -(l + 1);
      else if (twoJ == 2*l - 1 && l > 0) kappa = l;
      else return 1;
      p = end + 2;
      if (*p != '(') return 1;
      q = strtod(p + 1, &end);
      if (end == p + 1 || *end != ')') return 1;
      p = end + 1;
      G = twoJ + 1.;
    }
    else
    {
      q = strtod(p, &end);
      if (end == p) return 1;
      p = end;
      G = 2.*(2*l + 1);
    }
    if (n <= l || q < 0. || q > G) return 1;
    for (j = 0; j < *nShell; j++)
      if (shell[j].n == n && shell[j].l == l && (kappa == 0 || shell[j].kappa == 0 || shell[j].kappa == kappa)) return 1;
    if (*nShell >= maxShell) return 1;
    shell[*nShell].n = n;
    shell[*nShell].l = l;
    shell[*nShell].kappa = kappa;
    shell[*nShell].q = q;
    (*nShell)++;
  }
  return 0;
}

/* The ground-state configuration of the neutral atom Z = 1 .. 103, from the NIST Atomic Spectra
 * Database (ionization energies, "Ground Shells" of the neutral atoms, fetched 11 November 2026;
 * the cores [Cd] and [Hg] written out), in the syntax of AtomParseConfiguration; NULL for Z
 * outside 1 .. 103 or not an integer. */
const char *AtomGroundStateConfiguration(double Z)
{
  static const char *table[103] = {
    "1s1",
    "1s2",
    "1s2 2s1",
    "1s2 2s2",
    "1s2 2s2 2p1",
    "1s2 2s2 2p2",
    "1s2 2s2 2p3",
    "1s2 2s2 2p4",
    "1s2 2s2 2p5",
    "1s2 2s2 2p6",
    "[Ne] 3s1",
    "[Ne] 3s2",
    "[Ne] 3s2 3p1",
    "[Ne] 3s2 3p2",
    "[Ne] 3s2 3p3",
    "[Ne] 3s2 3p4",
    "[Ne] 3s2 3p5",
    "[Ne] 3s2 3p6",
    "[Ar] 4s1",
    "[Ar] 4s2",
    "[Ar] 3d1 4s2",
    "[Ar] 3d2 4s2",
    "[Ar] 3d3 4s2",
    "[Ar] 3d5 4s1",
    "[Ar] 3d5 4s2",
    "[Ar] 3d6 4s2",
    "[Ar] 3d7 4s2",
    "[Ar] 3d8 4s2",
    "[Ar] 3d10 4s1",
    "[Ar] 3d10 4s2",
    "[Ar] 3d10 4s2 4p1",
    "[Ar] 3d10 4s2 4p2",
    "[Ar] 3d10 4s2 4p3",
    "[Ar] 3d10 4s2 4p4",
    "[Ar] 3d10 4s2 4p5",
    "[Ar] 3d10 4s2 4p6",
    "[Kr] 5s1",
    "[Kr] 5s2",
    "[Kr] 4d1 5s2",
    "[Kr] 4d2 5s2",
    "[Kr] 4d4 5s1",
    "[Kr] 4d5 5s1",
    "[Kr] 4d5 5s2",
    "[Kr] 4d7 5s1",
    "[Kr] 4d8 5s1",
    "[Kr] 4d10",
    "[Kr] 4d10 5s1",
    "[Kr] 4d10 5s2",
    "[Kr] 4d10 5s2 5p1",
    "[Kr] 4d10 5s2 5p2",
    "[Kr] 4d10 5s2 5p3",
    "[Kr] 4d10 5s2 5p4",
    "[Kr] 4d10 5s2 5p5",
    "[Kr] 4d10 5s2 5p6",
    "[Xe] 6s1",
    "[Xe] 6s2",
    "[Xe] 5d1 6s2",
    "[Xe] 4f1 5d1 6s2",
    "[Xe] 4f3 6s2",
    "[Xe] 4f4 6s2",
    "[Xe] 4f5 6s2",
    "[Xe] 4f6 6s2",
    "[Xe] 4f7 6s2",
    "[Xe] 4f7 5d1 6s2",
    "[Xe] 4f9 6s2",
    "[Xe] 4f10 6s2",
    "[Xe] 4f11 6s2",
    "[Xe] 4f12 6s2",
    "[Xe] 4f13 6s2",
    "[Xe] 4f14 6s2",
    "[Xe] 4f14 5d1 6s2",
    "[Xe] 4f14 5d2 6s2",
    "[Xe] 4f14 5d3 6s2",
    "[Xe] 4f14 5d4 6s2",
    "[Xe] 4f14 5d5 6s2",
    "[Xe] 4f14 5d6 6s2",
    "[Xe] 4f14 5d7 6s2",
    "[Xe] 4f14 5d9 6s1",
    "[Xe] 4f14 5d10 6s1",
    "[Xe] 4f14 5d10 6s2",
    "[Xe] 4f14 5d10 6s2 6p1",
    "[Xe] 4f14 5d10 6s2 6p2",
    "[Xe] 4f14 5d10 6s2 6p3",
    "[Xe] 4f14 5d10 6s2 6p4",
    "[Xe] 4f14 5d10 6s2 6p5",
    "[Xe] 4f14 5d10 6s2 6p6",
    "[Rn] 7s1",
    "[Rn] 7s2",
    "[Rn] 6d1 7s2",
    "[Rn] 6d2 7s2",
    "[Rn] 5f2 6d1 7s2",
    "[Rn] 5f3 6d1 7s2",
    "[Rn] 5f4 6d1 7s2",
    "[Rn] 5f6 7s2",
    "[Rn] 5f7 7s2",
    "[Rn] 5f7 6d1 7s2",
    "[Rn] 5f9 7s2",
    "[Rn] 5f10 7s2",
    "[Rn] 5f11 7s2",
    "[Rn] 5f12 7s2",
    "[Rn] 5f13 7s2",
    "[Rn] 5f14 7s2",
    "[Rn] 5f14 7s2 7p1"};
  int z = (int) Z;
  if ((double) z != Z || z < 1 || z > 103) return NULL;
  return table[z - 1];
}

/* The mass number used when a finite nucleus is asked for without A or Rrms: that of the most
 * abundant isotope (NIST Atomic Weights and Isotopic Compositions, fetched 11 November 2026),
 * for elements without a stable isotope the mass number NIST gives in brackets (Tc, Pm, Po ..
 * Ac, Np, Pu), and for Am .. Lr, where NIST gives none, that of the longest-lived isotope.
 * 0 for Z outside 1 .. 103 or not an integer. */
int AtomDefaultMassNumber(double Z)
{
  static const int table[103] = {
    1, 4, 7, 9, 11, 12, 14, 16, 19, 20,
    23, 24, 27, 28, 31, 32, 35, 40, 39, 40,
    45, 48, 51, 52, 55, 56, 59, 58, 63, 64,
    69, 74, 75, 80, 79, 84, 85, 88, 89, 90,
    93, 98, 98, 102, 103, 106, 107, 114, 115, 120,
    121, 130, 127, 132, 133, 138, 139, 140, 141, 142,
    145, 152, 153, 158, 159, 164, 165, 166, 169, 174,
    175, 180, 181, 184, 187, 192, 193, 195, 197, 202,
    205, 208, 209, 209, 210, 222, 223, 226, 227, 232,
    231, 238, 237, 244, 243, 247, 247, 251, 252, 257,
    258, 259, 262};
  int z = (int) Z;
  if ((double) z != Z || z < 1 || z > 103) return 0;
  return table[z - 1];
}

// the degeneracy of a shell: 4l+2 for an nl shell, 2j+1 = 2|kappa| for a relativistic subshell
double AtomShellDegeneracy(const AtomShellType *shell)
{
  return (shell->kappa == 0) ? 4.*shell->l + 2. : 2.*abs(shell->kappa);
}

/* P, P', Q, Q' at the quadrature points of the orbital with coefficients coef (Dirac: 2n DKB
 * coefficients; Schroedinger: n B-spline coefficients, Q = Q' = 0). */
void AtomOrbitalOnPoints(const RadialBasisType *basis, int relativistic, int kappa, double c, const double *coef,
                         double *P, double *dP, double *Q, double *dQ)
{
  uint32_t n = basis->nSpline - 2, d = basis->degree, q, j, mu;
  int64_t i;
  double val[3*(BSplineMaxDegree+1)], h = 1./(2.*c), r, B, dB, d2B, a, b;
  for (q = 0; q < basis->nPoint; q++)
  {
    P[q] = dP[q] = Q[q] = dQ[q] = 0.;
    r = basis->rq[q];
    mu = BSplineFindInterval(r, basis->knot, basis->nKnot);
    if (mu >= basis->nKnot) continue;
    BSplineValues(r, mu, basis->knot, basis->nKnot, d, 2, val);
    for (j = 0; j <= d; j++)
    {
      i = (int64_t) mu - d + j - 1;
      if (i < 0 || i >= n) continue;
      B = val[j];
      dB = val[d + 1 + j];
      d2B = val[2*(d + 1) + j];
      a = coef[i];
      if (!relativistic)
      {
        P[q] += a*B;
        dP[q] += a*dB;
        continue;
      }
      b = coef[n + i];
      P[q]  += a*B + h*b*(dB - kappa*B/r);
      dP[q] += a*dB + h*b*(d2B - kappa*dB/r + kappa*B/(r*r));
      Q[q]  += h*a*(dB + kappa*B/r) + b*B;
      dQ[q] += h*a*(d2B + kappa*dB/r - kappa*B/(r*r)) + b*dB;
    }
  }
}

/* The basis functions that are nonzero at quadrature point q, with their values: idx[m] the
 * index (Dirac: 0 .. 2n-1, u^(1) then u^(2); Schroedinger: 0 .. n-1), P, P', Q, Q' (Q = 0
 * without relativity), for m = 0 .. (return value)-1, at most 2(degree+1). The same
 * pointwise expressions as DiracDKBMatrices (design section 5). */
uint32_t AtomLocalFunctions(const RadialBasisType *basis, int relativistic, int kappa, double c, uint32_t q,
                            uint32_t *idx, double *P, double *dP, double *Q, double *dQ)
{
  uint32_t n = basis->nSpline - 2, d = basis->degree, j, mu, m = 0;
  int64_t i;
  double val[3*(BSplineMaxDegree+1)], h = 1./(2.*c), r = basis->rq[q], B, dB, d2B;
  mu = BSplineFindInterval(r, basis->knot, basis->nKnot);
  if (mu >= basis->nKnot) return 0;
  BSplineValues(r, mu, basis->knot, basis->nKnot, d, 2, val);
  for (j = 0; j <= d; j++)
  {
    i = (int64_t) mu - d + j - 1;
    if (i < 0 || i >= n) continue;
    B = val[j];
    dB = val[d + 1 + j];
    d2B = val[2*(d + 1) + j];
    if (!relativistic)
    {
      idx[m] = (uint32_t) i; P[m] = B; dP[m] = dB; Q[m] = 0.; dQ[m] = 0.; m++;
      continue;
    }
    idx[m] = (uint32_t) i;                              // u_i^(1)
    P[m] = B;
    dP[m] = dB;
    Q[m] = h*(dB + kappa*B/r);
    dQ[m] = h*(d2B + kappa*dB/r - kappa*B/(r*r));
    m++;
    idx[m] = n + (uint32_t) i;                          // u_i^(2)
    P[m] = h*(dB - kappa*B/r);
    dP[m] = h*(d2B - kappa*dB/r + kappa*B/(r*r));
    Q[m] = B;
    dQ[m] = dB;
    m++;
  }
  return m;
}

/* M_ab = int [ v (P_a P_b + Q_a Q_b) + g (P_a P_b + Q_a Q_b)' ] dr over the basis functions
 * (Dirac: the 2n DKB functions of design section 5, M is 2n x 2n; Schroedinger: the n
 * B-splines B_1 .. B_n, Q = 0), dense, column major, both triangles. v and g on the
 * quadrature points; g may be NULL. Direct quadrature: on every point the 2(degree+1)
 * nonzero functions (AtomLocalFunctions), so that M[v, NULL] equals the V part of
 * DiracDKBMatrices to rounding. */
void AtomPotentialMatrix(const RadialBasisType *basis, int relativistic, int kappa, double c, const double *v, const double *g, double *M)
{
  uint32_t n = basis->nSpline - 2, Nf = relativistic ? 2*n : n, q, j, k, m, idx[2*(BSplineMaxDegree+1)];
  double P[2*(BSplineMaxDegree+1)], dP[2*(BSplineMaxDegree+1)], Q[2*(BSplineMaxDegree+1)], dQ[2*(BSplineMaxDegree+1)], wv, wg;
  for (j = 0; j < Nf*Nf; j++) M[j] = 0.;
  for (q = 0; q < basis->nPoint; q++)
  {
    m = AtomLocalFunctions(basis, relativistic, kappa, c, q, idx, P, dP, Q, dQ);
    wv = basis->wq[q]*v[q];
    wg = (g != NULL) ? basis->wq[q]*g[q] : 0.;
    for (k = 0; k < m; k++)
      for (j = 0; j < m; j++)
        M[idx[j] + idx[k]*Nf] += wv*(P[j]*P[k] + Q[j]*Q[k]) + wg*(dP[j]*P[k] + P[j]*dP[k] + dQ[j]*Q[k] + Q[j]*dQ[k]);
  }
}

/* The kinetic part H0 and the overlap S of one symmetry group (dense, Nf x Nf): Dirac,
 * DiracDKBMatrices with V = 0; Schroedinger, H0 = 1/2 int B_i' B_j' + l(l+1)/2 int B_i B_j/r^2,
 * S = int B_i B_j, over B_1 .. B_n. */
int AtomKineticMatrices(const RadialBasisType *basis, int relativistic, int kappa, int l, double c, double *H0, double *S)
{
  uint32_t n = basis->nSpline - 2, q, i;
  double *zero, *f, *M;
  int error = 0;
  zero = (double *) calloc(basis->nPoint, sizeof(double));
  if (relativistic) error = DiracDKBMatrices(basis, kappa, c, zero, NULL, H0, S);
  else
  {
    f = (double *) malloc(basis->nPoint*sizeof(double));
    M = (double *) malloc(n*n*sizeof(double));
    for (q = 0; q < basis->nPoint; q++) f[q] = 0.5*l*(l + 1.)/(basis->rq[q]*basis->rq[q]);
    error |= DiracDenseMatrix(basis, NULL, 1, 1, H0);
    for (i = 0; i < n*n; i++) H0[i] *= 0.5;
    error |= DiracDenseMatrix(basis, f, 0, 0, M);
    for (i = 0; i < n*n; i++) H0[i] += M[i];
    error |= DiracDenseMatrix(basis, NULL, 0, 0, S);
    free(f);
    free(M);
  }
  free(zero);
  return error;
}

// sum_ij x_i A_ij y_j for the Nf x Nf matrix A
double AtomBilinear(const double *x, const double *A, const double *y, uint32_t Nf)
{
  double s = 0., t;
  uint32_t i, j;
  for (j = 0; j < Nf; j++)
  {
    t = 0.;
    for (i = 0; i < Nf; i++) t += x[i]*A[i + j*Nf];
    s += t*y[j];
  }
  return s;
}

/* Thomas-Fermi screening, the start of the iteration: the potential of the electrons
 * Z (1 - phi(x))/r, x = r/(0.8853 Z^(-1/3)), with Latter's fit of the Thomas-Fermi function
 * phi. Only a starting point; its details do not enter the result. */
DFTReal AtomThomasFermiScreening(DFTReal Z, DFTReal r)
{
  DFTReal x = r*pow(Z, 1./3.)/0.8853, s = sqrt(x);
  DFTReal phi = 1./(1. + 0.02747*s + 1.243*x - 0.1486*x*s + 0.2302*x*x + 0.007298*x*x*s + 0.006944*x*x*x);
  return Z*(1. - phi)/r;
}

/* The density of the orbitals, as the radial density f = 4 pi r^2 n = sum qbar (P^2 + Q^2)
 * and its derivative df, on the quadrature points (atom->f, atom->df). */
void AtomDensity(DiracAtomType *atom)
{
  RadialBasisType *b = atom->basis;
  uint32_t q, o, nP = b->nPoint;
  double *P, *dP, *Q, *dQ;
  P = (double *) malloc(nP*sizeof(double));
  dP = (double *) malloc(nP*sizeof(double));
  Q = (double *) malloc(nP*sizeof(double));
  dQ = (double *) malloc(nP*sizeof(double));
  for (q = 0; q < nP; q++) atom->f[q] = atom->df[q] = 0.;
  for (o = 0; o < atom->nOrb; o++)
  {
    if (atom->orb[o].q == 0.) continue;
    AtomOrbitalOnPoints(b, atom->relativistic, atom->orb[o].kappa, atom->c, atom->orb[o].coef, P, dP, Q, dQ);
    for (q = 0; q < nP; q++)
    {
      atom->f[q] += atom->orb[o].q*(P[q]*P[q] + Q[q]*Q[q]);
      atom->df[q] += 2.*atom->orb[o].q*(P[q]*dP[q] + Q[q]*dQ[q]);
    }
  }
  free(P); free(dP); free(Q); free(dQ);
}

/* The potentials of the density (f, df) on the quadrature points: V_H = Y^0[f]/r, the
 * exchange-correlation energy density e, v = V_H + de/dn - 4 (de/dsigma) n'/r and
 * g = 2 (de/dsigma) n' (0 for an LDA), with n = f/(4 pi r^2), n' = (df - 2f/r)/(4 pi r^2).
 * Also the energies Ecoul = 1/2 int f V_H and Exc = int 4 pi r^2 e. */
int AtomPotential(DiracAtomType *atom, const double *f, const double *df, double *v, double *g, double *Ecoul, double *Exc)
{
  RadialBasisType *b = atom->basis;
  uint32_t q, nP = b->nPoint, ld = b->degree + 1;
  double *Y, n, dn, e, en, es, r, sumC = 0., cC = 0., sumX = 0., cX = 0.;
  Y = (double *) malloc(b->nSpline*sizeof(double));
  if (BSplinePoissonSolve(b->knot, b->nKnot, b->degree, b->rq, b->wq, nP, 0, atom->poisson, ld, f, Y)) { free(Y); return 1; }
  for (q = 0; q < nP; q++)
  {
    r = b->rq[q];
    v[q] = BSplineDerivative(r, b->knot, b->nKnot, b->degree, Y, 0)/r;      // V_H
    NeumaierAdd(&sumC, &cC, 0.5*b->wq[q]*f[q]*v[q]);
    n = f[q]/(4.*M_PI*r*r);
    dn = (df[q] - 2.*f[q]/r)/(4.*M_PI*r*r);
    XCEvaluate(atom->method, atom->c, n, dn*dn, atom->cutoff, &e, &en, &es);
    NeumaierAdd(&sumX, &cX, b->wq[q]*4.*M_PI*r*r*e);
    v[q] += en - 4.*es*dn/r;
    g[q] = 2.*es*dn;
  }
  *Ecoul = sumC + cC;
  *Exc = sumX + cX;
  free(Y);
  return 0;
}

// AtomDensity, then AtomPotential of that density
int AtomDensityPotential(DiracAtomType *atom, double *v, double *g, double *Ecoul, double *Exc)
{
  AtomDensity(atom);
  return AtomPotential(atom, atom->f, atom->df, v, g, Ecoul, Exc);
}

/* Set up the atom for the shells (AtomParseConfiguration; copied) on the basis with the
 * nuclear potential VN (on the quadrature points; copied): the subshells (Dirac: an nl
 * shell gives kappa = -(l+1) and l with qbar = q (2j+1)/(4l+2), an n kappa subshell one;
 * Schroedinger: one per nl shell), their symmetry groups with H0, S and H0 + M[V_N], and
 * the Poisson factorization (k = 0). Shells with q = 0 give empty orbitals (no density;
 * the SCF still solves for them); a configuration of empty shells only is the bare nucleus.
 * Returns 1 for an n kappa subshell without relativity or a failing setup. */
int AtomSetup(DiracAtomType *atom, RadialBasisType *basis, const DFTReal *VN, double Z, int relativistic, double c, int method,
              int originBoundary, DFTReal cutoff, uint32_t nShell, const AtomShellType *shell)
{
  uint32_t s, o, k, nf, nP = basis->nPoint, n = basis->nSpline - 2, nLarge;
  int found, kap, l, t;
  double *M;
  memset(atom, 0, sizeof(DiracAtomType));
  atom->basis = basis;
  atom->Z = Z;
  atom->relativistic = relativistic;
  atom->c = c;
  atom->method = method;
  atom->cutoff = cutoff;
  atom->Nf = relativistic ? 2*n : n;
  nf = atom->Nf;
  atom->orb = (DiracOrbitalType *) calloc(2*nShell + 1, sizeof(DiracOrbitalType));
  atom->group = (DiracGroupType *) calloc(2*nShell + 1, sizeof(DiracGroupType));
  atom->nShell = nShell;
  atom->shell = (AtomShellType *) malloc((nShell + 1)*sizeof(AtomShellType));
  for (s = 0; s < nShell; s++) atom->shell[s] = shell[s];
  for (s = 0; s < nShell; s++)
  {
    l = shell[s].l;
    if (shell[s].kappa != 0 && !relativistic) return 1;
    for (t = 0; t < (relativistic && l > 0 && shell[s].kappa == 0 ? 2 : 1); t++)
    {
      DiracOrbitalType *orb = atom->orb + atom->nOrb;
      orb->n = shell[s].n;
      orb->l = l;
      orb->shell = s;
      if (relativistic)
      {
        if (shell[s].kappa != 0) kap = shell[s].kappa;
        else kap = (t == 0) ? -(l + 1) : l;               // j = l+1/2 first, then l-1/2
        orb->kappa = kap;
        orb->q = shell[s].q*2.*abs(kap)/AtomShellDegeneracy(shell + s);
      }
      else
      {
        orb->kappa = 0;
        orb->q = shell[s].q;
      }
      orb->coef = (double *) calloc(nf, sizeof(double));
      // its symmetry group
      found = 0;
      for (k = 0; k < atom->nGroup && !found; k++)
        if (atom->group[k].kappa == orb->kappa && atom->group[k].l == l) { orb->group = k; found = 1; }
      if (!found)
      {
        DiracGroupType *gr = atom->group + atom->nGroup;
        gr->kappa = orb->kappa;
        gr->l = l;
        gr->H0 = (double *) malloc(nf*nf*sizeof(double));
        gr->S = (double *) malloc(nf*nf*sizeof(double));
        gr->Hfix = (double *) malloc(nf*nf*sizeof(double));
        gr->active = (uint32_t *) malloc(nf*sizeof(uint32_t));
        if (AtomKineticMatrices(basis, relativistic, orb->kappa, l, c, gr->H0, gr->S)) return 1;
        M = (double *) malloc(nf*nf*sizeof(double));
        AtomPotentialMatrix(basis, relativistic, orb->kappa, c, VN, NULL, M);
        for (k = 0; k < nf*nf; k++) gr->Hfix[k] = gr->H0[k] + M[k];
        free(M);
        if (relativistic) DiracActiveFunctions(basis, orb->kappa, originBoundary, gr->active, &(gr->N), &nLarge);
        else
        {
          gr->N = n;
          for (k = 0; k < n; k++) gr->active[k] = k;
        }
        orb->group = atom->nGroup;
        atom->nGroup++;
      }
      atom->nElectron += orb->q;
      atom->nOrb++;
    }
  }
  // the roots: n - l - 1 within the electronic branch of the group
  for (o = 0; o < atom->nOrb; o++) atom->orb[o].root = (uint32_t) (atom->orb[o].n - atom->orb[o].l - 1);
  atom->VN = (DFTReal *) malloc(nP*sizeof(DFTReal));
  for (k = 0; k < nP; k++) atom->VN[k] = VN[k];
  atom->f = (DFTReal *) calloc(nP, sizeof(DFTReal));
  atom->df = (DFTReal *) calloc(nP, sizeof(DFTReal));
  atom->v = (DFTReal *) calloc(nP, sizeof(DFTReal));
  atom->g = (DFTReal *) calloc(nP, sizeof(DFTReal));
  atom->poisson = (double *) malloc((basis->degree + 1)*(basis->nSpline - 1)*sizeof(double));
  return BSplinePoissonSetup(basis->knot, basis->nKnot, basis->degree, basis->rq, basis->wq, nP, 0, atom->poisson, basis->degree + 1);
}

void AtomFree(DiracAtomType *atom)
{
  uint32_t k;
  for (k = 0; k < atom->nOrb; k++) { free(atom->orb[k].coef); free(atom->orb[k].extra); }
  free(atom->shell);
  for (k = 0; k < atom->nGroup; k++)
  {
    free(atom->group[k].H0);
    free(atom->group[k].S);
    free(atom->group[k].Hfix);
    free(atom->group[k].active);
  }
  free(atom->orb);
  free(atom->group);
  free(atom->VN);
  free(atom->f);
  free(atom->df);
  free(atom->v);
  free(atom->g);
  free(atom->poisson);
}

/* One diagonalization of every group in the potential (v, g): the orbitals (coefficients,
 * energies). Returns 1 if an eigenproblem fails or a group has too few electronic states. */
int AtomSolveOrbitals(DiracAtomType *atom, const double *v, const double *g, int rayleigh)
{
  uint32_t k, o, i, nf = atom->Nf, first;
  double *H, *M, *E, *C;
  int error = 0;
  H = (double *) malloc(nf*nf*sizeof(double));
  M = (double *) malloc(nf*nf*sizeof(double));
  E = (double *) malloc(nf*sizeof(double));
  C = (double *) malloc(nf*nf*sizeof(double));
  for (k = 0; k < atom->nGroup && !error; k++)
  {
    DiracGroupType *gr = atom->group + k;
    AtomPotentialMatrix(atom->basis, atom->relativistic, gr->kappa, atom->c, v, XCIsGGA(atom->method) ? g : NULL, M);
    for (i = 0; i < nf*nf; i++) H[i] = gr->Hfix[i] + M[i];
    if (DiracGeneralizedEigen(H, gr->S, nf, gr->active, gr->N, rayleigh, E, C, NULL, NULL)) { error = 1; break; }
    first = atom->relativistic ? DiracElectronicFirst(E, gr->N, atom->c) : 0;
    for (o = 0; o < atom->nOrb; o++)
    {
      if (atom->orb[o].group != k) continue;
      if (first + atom->orb[o].root >= gr->N) { error = 1; break; }
      atom->orb[o].energy = E[first + atom->orb[o].root];
      for (i = 0; i < nf; i++) atom->orb[o].coef[i] = C[i + (first + atom->orb[o].root)*nf];
    }
  }
  free(H); free(M); free(E); free(C);
  return error;
}

/* The self-consistent field (design section 9). Starts from Thomas-Fermi screening and
 * mixes (v, g) on the quadrature points by Pulay's DIIS with the last `history` iterations
 * and the step mixing (x_new = sum_i c_i (x_i + mixing R_i)). The residual R = out - in is
 * weighted with the density, as it enters the energy (first order int f dv):
 *     |R|^2 = sum_q w_q f_q (R_v^2 + R_g^2/r_q^2),   f = 4 pi r^2 n,
 * (g enters as int g (P_a P_b)' dr). An r^2 weight instead stalled the Dirac atoms at 1e-6:
 * beyond r ~ 10 the density falls to its rounding level (1e-24 for Ne) and v_xc ~ n^(1/3)
 * there is noise without effect on anything. Converged when the change of the total
 * energy is below energyTolerance and the residual below residualTolerance, at most
 * maxIterations. On return the orbitals belong to the last input potential; the energies
 * (DiracAtomType) to the density of these orbitals:
 *   Etot = Ekin + Eenuc + Ecoul + Exc with Ekin = sum qbar <a|H0|a>, Eenuc = int f V_N,
 *   Eeig = sum qbar eps_a - int (f v_out + df g_out) + Ecoul + Exc, which equals Etot at
 *   self-consistency (their difference is first order in the residual).
 * Returns 1 if a step fails; atom->converged tells whether it converged. */
int AtomSCF(DiracAtomType *atom, uint32_t maxIterations, uint32_t history, double mixing, double energyTolerance,
            double residualTolerance, int rayleigh, int mixDensity)
{
  RadialBasisType *b = atom->basis;
  uint32_t nP = b->nPoint, nx = 2*nP, it, i, j, k, o, m = 0, slot, nrhs = 1, nB;
  double *x, *R, *X, *RR, *vout, *gout, *Bm, *rhs, *W, *W2, *vin, *gin, Ecoul, Exc, Ekin, Eenuc, Eeig, Ebefore = 0., norm, sum, cs;
  int *ipiv, info, error = 0;
  if (history == 0) history = 1;
  x = (double *) malloc(nx*sizeof(double));                 // input (v, g)
  R = (double *) malloc(nx*sizeof(double));
  X = (double *) malloc(history*nx*sizeof(double));          // stored inputs
  RR = (double *) malloc(history*nx*sizeof(double));         // stored residuals
  vout = (double *) malloc(nP*sizeof(double));
  W = (double *) malloc(nP*sizeof(double));                   // weights of the residual norm
  W2 = (double *) malloc(nP*sizeof(double));
  gout = (double *) malloc(nP*sizeof(double));
  Bm = (double *) malloc((history + 1)*(history + 1)*sizeof(double));
  rhs = (double *) malloc((history + 1)*sizeof(double));
  ipiv = (int *) malloc((history + 1)*sizeof(int));
  vin = (double *) malloc(nP*sizeof(double));
  gin = (double *) calloc(nP, sizeof(double));
  for (i = 0; i < nP; i++) vin[i] = AtomThomasFermiScreening(atom->Z, b->rq[i]);
  if (mixDensity)
  {
    // the density of the Thomas-Fermi orbitals is the first input
    if (AtomSolveOrbitals(atom, vin, gin, rayleigh)) error = 1;
    AtomDensity(atom);
    for (i = 0; i < nP; i++) { x[i] = atom->f[i]; x[nP + i] = atom->df[i]; }
  }
  else
    for (i = 0; i < nP; i++) { x[i] = vin[i]; x[nP + i] = 0.; }
  atom->converged = 0;
  for (it = 1; it <= maxIterations && !error; it++)
  {
    if (mixDensity && AtomPotential(atom, x, x + nP, vin, gin, &Ecoul, &Exc)) { error = 1; break; }
    if (AtomSolveOrbitals(atom, mixDensity ? vin : x, mixDensity ? gin : x + nP, rayleigh)) { error = 1; break; }
    if (AtomDensityPotential(atom, vout, gout, &Ecoul, &Exc)) { error = 1; break; }
    // energies
    Ekin = 0.;
    Eeig = 0.;
    for (o = 0; o < atom->nOrb; o++)
    {
      Ekin += atom->orb[o].q*AtomBilinear(atom->orb[o].coef, atom->group[atom->orb[o].group].H0, atom->orb[o].coef, atom->Nf);
      Eeig += atom->orb[o].q*atom->orb[o].energy;
    }
    Eenuc = 0.;
    for (i = 0; i < nP; i++)
    {
      Eenuc += b->wq[i]*atom->f[i]*atom->VN[i];
      Eeig -= b->wq[i]*(atom->f[i]*vout[i] + atom->df[i]*gout[i]);
    }
    if (!mixDensity) for (i = 0; i < nP; i++) { vin[i] = x[i]; gin[i] = x[nP + i]; }
    atom->Ekin = Ekin;
    atom->Eenuc = Eenuc;
    atom->Ecoul = Ecoul;
    atom->Exc = Exc;
    atom->Etot = Ekin + Eenuc + Ecoul + Exc;
    atom->Eeig = Eeig + Ecoul + Exc;
    // residual
    // residual of the mixed quantity: potential (v, g) weighted with the density, or the
    // density (f, df) with weights w and w r^2
    norm = 0.;
    for (i = 0; i < nP; i++)
    {
      if (mixDensity)
      {
        R[i] = atom->f[i] - x[i];
        R[nP + i] = atom->df[i] - x[nP + i];
        W[i] = b->wq[i];
        W2[i] = b->wq[i]*b->rq[i]*b->rq[i];
      }
      else
      {
        R[i] = vout[i] - x[i];
        R[nP + i] = gout[i] - x[nP + i];
        W[i] = b->wq[i]*atom->f[i];
        W2[i] = W[i]/(b->rq[i]*b->rq[i]);
      }
      norm += W[i]*R[i]*R[i] + W2[i]*R[nP + i]*R[nP + i];
    }
    atom->residual = sqrt(norm);
    // the potential residual in both modes, weighted with the output density: it bounds the
    // difference of the two energy expressions, |Eeig - Etot| <= sqrt(N) |R_v| + ... (to
    // first order)
    norm = 0.;
    for (i = 0; i < nP; i++)
    {
      double rv = vout[i] - (mixDensity ? vin[i] : x[i]), rg = gout[i] - (mixDensity ? gin[i] : x[nP + i]);
      norm += b->wq[i]*atom->f[i]*(rv*rv + rg*rg/(b->rq[i]*b->rq[i]));
    }
    atom->potentialResidual = sqrt(norm);
    atom->dE = atom->Etot - Ebefore;
    atom->iterations = it;
    Ebefore = atom->Etot;
    if (it > 1 && fabs(atom->dE) < energyTolerance && atom->residual < residualTolerance)
    {
      atom->converged = 1;
      break;
    }
    if (it == maxIterations) break;
    // Pulay mixing
    slot = (it - 1) % history;
    for (i = 0; i < nx; i++)
    {
      X[slot*nx + i] = x[i];
      RR[slot*nx + i] = R[i];
    }
    m = (it < history) ? it : history;
    nB = m + 1;
    for (j = 0; j < m; j++)
      for (k = 0; k <= j; k++)
      {
        sum = 0.;
        for (i = 0; i < nP; i++) sum += W[i]*RR[j*nx + i]*RR[k*nx + i] + W2[i]*RR[j*nx + nP + i]*RR[k*nx + nP + i];
        Bm[j + k*nB] = Bm[k + j*nB] = sum;
      }
    for (j = 0; j < m; j++) { Bm[m + j*nB] = Bm[j + m*nB] = -1.; rhs[j] = 0.; }
    Bm[m + m*nB] = 0.;
    rhs[m] = -1.;
    dgesv(&nB, &nrhs, Bm, &nB, ipiv, rhs, &nB, &info);
    if (info != 0)
    {
      // singular history: plain mixing
      for (i = 0; i < nx; i++) x[i] += mixing*R[i];
      continue;
    }
    for (i = 0; i < nx; i++)
    {
      cs = 0.;
      for (j = 0; j < m; j++) cs += rhs[j]*(X[j*nx + i] + mixing*RR[j*nx + i]);
      x[i] = cs;
    }
  }
  for (i = 0; i < nP; i++)
  {
    atom->v[i] = vin[i];
    atom->g[i] = gin[i];
  }
  free(x); free(R); free(X); free(RR); free(vout); free(gout); free(Bm); free(rhs); free(ipiv); free(W); free(W2); free(vin); free(gin);
  return error;
}

/* Sign convention of orbital o, P > 0 on its innermost lobe (InnermostLobeSign), for Dirac and
 * Schroedinger orbitals alike. */
void AtomSignConvention(DiracAtomType *atom, uint32_t o)
{
  uint32_t q, nP = atom->basis->nPoint;
  double *P = (double *) malloc(nP*sizeof(double)), *dP = (double *) malloc(nP*sizeof(double));
  double *Q = (double *) malloc(nP*sizeof(double)), *dQ = (double *) malloc(nP*sizeof(double));
  AtomOrbitalOnPoints(atom->basis, atom->relativistic, atom->orb[o].kappa, atom->c, atom->orb[o].coef, P, dP, Q, dQ);
  if (InnermostLobeSign(P, nP) < 0.) for (q = 0; q < atom->Nf; q++) atom->orb[o].coef[q] = -atom->orb[o].coef[q];
  free(P); free(dP); free(Q); free(dQ);
}
