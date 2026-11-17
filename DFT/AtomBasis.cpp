/*
 * AtomBasis.cpp
 *
 * Atomic basis functions (DFT stage A5, design sections 6 and 6.1; decisions of 2026-11-11
 * in section 0.2): every orbital (n, kappa) of the atom may have an extra nuclear charge Q
 * and a confinement radius rc. Its function comes from its own generating calculation, in
 * which only its own equation has the nuclear charge Z+Q and the Lorentz-scalar
 * confinement V_beta = (r/rc)^N; the other electrons see Z and no confinement:
 *  - an occupied orbital: the self-consistent field (HFSCF, Hartree-Fock or Kohn-Sham) with
 *    the extra one-electron term Delta H = M[(Q/Z) V_N] + M_beta[V_beta] in its equation,
 *    the stationary point of E + qbar <a|Delta H|a> with the orbitals orthonormal (the
 *    target its own class, the other occupied orbitals of its kappa projected out);
 *  - an empty orbital (q = 0): the root n-l-1 of its ground-state operator plus Delta H
 *    (Hartree-Fock: the Fock operator of V^N; Kohn-Sham: the ground-state potential).
 * The other orbitals are those of the ground state. The functions of different generating
 * calculations are not orthogonal; AtomOrthonormalizeBasis makes them orthonormal per kappa,
 * Gram-Schmidt from the lowest n upwards (decided by Maurits).
 *
 * All routines return 0 on success and 1 on error and do not print.
 */

#include "AtomBasis.h"

/* The extra one-electron term of orbital o with nuclear charge Z+Q and the confinement
 * Vbeta (on the quadrature points, NULL for none): extra = M[(Q/Z) V_N] + M_beta[V_beta]
 * (Dirac: the DKB matrices of design section 5, the scalar enters as beta V_beta;
 * Schroedinger: both as potentials). beta (NULL allowed) receives M_beta[V_beta] alone, for
 * <beta V_beta> = int V_beta (P^2 - Q^2) dr. Nf x Nf each. */
int AtomTargetMatrices(DiracAtomType *atom, uint32_t o, double Q, const DFTReal *Vbeta, double *extra, double *beta)
{
  RadialBasisType *b = atom->basis;
  DiracGroupType *g = atom->group + atom->orb[o].group;
  uint32_t i, Nf = atom->Nf, nP = b->nPoint;
  DFTReal *V = (DFTReal *) malloc(nP*sizeof(DFTReal)), *zero = (DFTReal *) calloc(nP, sizeof(DFTReal));
  double *H = (double *) malloc(Nf*Nf*sizeof(double)), *S = (double *) malloc(Nf*Nf*sizeof(double));
  int error = 0;
  for (i = 0; i < nP; i++) V[i] = (Q/atom->Z)*atom->VN[i];
  if (atom->relativistic)
  {
    error |= DiracDKBMatrices(b, atom->orb[o].kappa, atom->c, V, Vbeta, H, S);
    for (i = 0; i < Nf*Nf; i++) extra[i] = H[i] - g->H0[i];
    if (beta != NULL)
    {
      if (Vbeta == NULL) for (i = 0; i < Nf*Nf; i++) beta[i] = 0.;
      else
      {
        error |= DiracDKBMatrices(b, atom->orb[o].kappa, atom->c, zero, Vbeta, H, S);
        for (i = 0; i < Nf*Nf; i++) beta[i] = H[i] - g->H0[i];
      }
    }
  }
  else
  {
    if (Vbeta != NULL) for (i = 0; i < nP; i++) V[i] += Vbeta[i];
    AtomPotentialMatrix(b, 0, 0, atom->c, V, NULL, extra);
    if (beta != NULL)
    {
      if (Vbeta == NULL) for (i = 0; i < Nf*Nf; i++) beta[i] = 0.;
      else AtomPotentialMatrix(b, 0, 0, atom->c, Vbeta, NULL, beta);
    }
  }
  free(V); free(zero); free(H); free(S);
  return error;
}

/* The root `root` (n-l-1) of the electronic branch of F (Nf x Nf) in symmetry group `group`:
 * its coefficients and eigenvalue. Returns 1 if the eigenproblem fails or has too few
 * electronic states. */
int AtomRootSolve(DiracAtomType *atom, uint32_t group, const double *F, uint32_t root, int rayleigh, double *coef, double *energy)
{
  DiracGroupType *g = atom->group + group;
  uint32_t Nf = atom->Nf, i, first;
  double *E = (double *) malloc(Nf*sizeof(double)), *C = (double *) malloc(Nf*Nf*sizeof(double));
  int error = DiracGeneralizedEigen(F, g->S, Nf, g->active, g->N, rayleigh, E, C, NULL, NULL);
  first = atom->relativistic ? DiracElectronicFirst(E, g->N, atom->c) : 0;
  if (!error && first + root < g->N)
  {
    for (i = 0; i < Nf; i++) coef[i] = C[i + (first + root)*Nf];
    *energy = E[first + root];
  }
  else error = 1;
  free(E); free(C);
  return error;
}

void AtomBasisFree(AtomBasisType *res, uint32_t nOrb)
{
  uint32_t o;
  if (res->Fground != NULL) for (o = 0; o < nOrb; o++) free(res->Fground[o]);
  free(res->Fground);
  free(res->coef); free(res->energy); free(res->energyGround); free(res->beta); free(res->Egenerating);
  free(res->generated); free(res->converged);
  memset(res, 0, sizeof(AtomBasisType));
}

/* The basis functions of the atom (set up by AtomSetup with every orbital wanted; the
 * method of the atom is the Kohn-Sham functional, or the LDA start of Hartree-Fock):
 *  1. the ground state: AtomSCF, then (hartreeFock) HFSCF; the empty orbitals from the
 *     ground-state operators (V^N Fock operator, or the Kohn-Sham potential);
 *  2. per orbital o with Q[o] != 0 or rc[o] > 0 its generating calculation (file header),
 *     from the ground-state orbitals;
 *  3. <a|F_ground|a> of every function is left to the caller after orthonormalization
 *     (res->Fground holds the ground-state matrices).
 * On return the atom holds the ground state; res the functions (not yet orthonormal).
 * Returns 1 if a step fails; res->groundConverged and res->converged[o] tell convergence. */
int AtomBasisFunctions(DiracAtomType *atom, int hartreeFock, const double *Q, const double *rc, double Nconf,
                       const AtomSCFOptionsType *opt, AtomBasisType *res)
{
  RadialBasisType *b = atom->basis;
  DiracHFType hf;
  uint32_t o, i, Nf = atom->Nf, nOrb = atom->nOrb, nP = b->nPoint;
  double *C0, *Mbeta, Eground;
  DFTReal *Vbeta = NULL;
  int error = 0;
  memset(res, 0, sizeof(AtomBasisType));
  res->coef = (double *) calloc(nOrb*Nf, sizeof(double));
  res->energy = (double *) calloc(nOrb, sizeof(double));
  res->energyGround = (double *) calloc(nOrb, sizeof(double));
  res->beta = (double *) calloc(nOrb, sizeof(double));
  res->Egenerating = (double *) calloc(nOrb, sizeof(double));
  res->generated = (int *) calloc(nOrb, sizeof(int));
  res->converged = (int *) calloc(nOrb, sizeof(int));
  res->Fground = (double **) calloc(nOrb, sizeof(double *));
  for (o = 0; o < nOrb; o++) res->Fground[o] = (double *) malloc(Nf*Nf*sizeof(double));
  C0 = (double *) malloc(nOrb*Nf*sizeof(double));
  Mbeta = (double *) malloc(Nf*Nf*sizeof(double));

  // 1. the ground state; without electrons the bare nucleus: every orbital a root of H0 + V_N
  if (!(atom->nElectron > 0.))
    for (o = 0; o < nOrb && !error; o++)
      error = AtomRootSolve(atom, atom->orb[o].group, atom->group[atom->orb[o].group].Hfix, atom->orb[o].root, opt->rayleigh,
                            atom->orb[o].coef, &(atom->orb[o].energy));
  else if (hartreeFock) error = AtomSCF(atom, opt->maxIterations, 8, 0.3, 1e-8, 1e-6, opt->rayleigh, 0);   // the LDA start
  else error = AtomSCF(atom, opt->maxIterations, opt->history, opt->mixing, opt->energyTolerance, opt->residualTolerance,
                       opt->rayleigh, opt->mixDensity);
  if (!error) error = HFSetup(&hf, atom, 0);
  if (error) { free(C0); free(Mbeta); return 1; }
  hf.ks = !hartreeFock;
  if (!(atom->nElectron > 0.)) res->groundConverged = 1;   // nothing to iterate; the energies are 0
  else if (hartreeFock)
  {
    error = HFSCF(&hf, opt->maxIterations, opt->history, opt->mixing, opt->energyTolerance, opt->orbitalTolerance, opt->rayleigh);
    res->groundConverged = hf.converged;
    res->iterations = hf.iterations;
    res->dE = hf.dE;
    res->residual = hf.residual;
    res->gradient = hf.gradient;
    res->pairGradient = hf.pairGradient;
    res->Eeig = hf.Eidentity;
    if (!error && !atom->relativistic)
    {
      // Cowan's E_av again with the relativistic angular factors of j = l -+ 1/2 (design check 5)
      DiracHFType jj;
      double e1, e2, e3, e4;
      error = HFSetup(&jj, atom, 1);
      if (!error)
      {
        HFOrbitalsOnPoints(&jj);
        error = HFEnergy(&jj, &res->EavJJ, &e1, &e2, &e3, &e4);
        HFFree(&jj);
      }
    }
  }
  else
  {
    res->groundConverged = atom->converged;
    res->iterations = atom->iterations;
    res->dE = atom->dE;
    res->residual = atom->residual;
    res->Eeig = atom->Eeig;
    res->potentialResidual = atom->potentialResidual;
  }
  // the ground-state operators of the ground-state orbitals; the empty orbitals of Hartree-Fock
  if (!error)
  {
    HFOrbitalsOnPoints(&hf);
    error = HFFockMatrices(&hf, res->Fground);
  }
  for (o = 0; o < nOrb && !error; o++)
    if (hartreeFock && hf.sub[o].empty)
      error = AtomRootSolve(atom, atom->orb[o].group, res->Fground[o], atom->orb[o].root, opt->rayleigh, atom->orb[o].coef, &(atom->orb[o].energy));
    else if (hartreeFock) atom->orb[o].energy = AtomBilinear(atom->orb[o].coef, res->Fground[o], atom->orb[o].coef, Nf);
  if (!error)
  {
    HFOrbitalsOnPoints(&hf);
    error = HFEnergy(&hf, &res->Etot, &res->Ekin, &res->Eenuc, &res->Edirect, &res->Eexchange);
  }
  Eground = res->Etot;
  for (o = 0; o < nOrb; o++)
  {
    for (i = 0; i < Nf; i++) res->coef[o*Nf + i] = C0[o*Nf + i] = atom->orb[o].coef[i];
    res->energy[o] = atom->orb[o].energy;
    res->Egenerating[o] = Eground;
    res->converged[o] = 1;
  }

  // 2. the generating calculations
  if (Nconf > 0.) Vbeta = (DFTReal *) malloc(nP*sizeof(DFTReal));
  for (o = 0; o < nOrb && !error; o++)
  {
    double Etot, e1, e2, e3, e4;
    if (Q[o] == 0. && !(rc[o] > 0.)) continue;
    res->generated[o] = 1;
    if (rc[o] > 0. && (Vbeta == NULL || ConfinementPotential(b, rc[o], Nconf, Vbeta))) { error = 1; break; }
    atom->orb[o].extra = (double *) malloc(Nf*Nf*sizeof(double));
    error = AtomTargetMatrices(atom, o, Q[o], rc[o] > 0. ? Vbeta : NULL, atom->orb[o].extra, Mbeta);
    for (uint32_t a = 0; a < nOrb; a++) for (i = 0; i < Nf; i++) atom->orb[a].coef[i] = C0[a*Nf + i];
    if (!error && !hf.sub[o].empty)
    {
      error = HFSCF(&hf, opt->maxIterations, opt->history, opt->mixing, opt->energyTolerance, opt->orbitalTolerance, opt->rayleigh);
      res->converged[o] = hf.converged;
      res->energy[o] = hf.sub[o].energy;
      res->Egenerating[o] = hf.Etot;
    }
    else if (!error)
    {
      double **F = (double **) malloc(nOrb*sizeof(double *));
      for (uint32_t a = 0; a < nOrb; a++) F[a] = (double *) malloc(Nf*Nf*sizeof(double));
      HFOrbitalsOnPoints(&hf);
      error = HFFockMatrices(&hf, F);
      if (!error) error = AtomRootSolve(atom, atom->orb[o].group, F[o], atom->orb[o].root, opt->rayleigh, atom->orb[o].coef, &(res->energy[o]));
      if (!error) error = HFEnergy(&hf, &Etot, &e1, &e2, &e3, &e4);
      res->Egenerating[o] = Eground;    // an empty orbital does not change the energy
      for (uint32_t a = 0; a < nOrb; a++) free(F[a]);
      free(F);
    }
    for (i = 0; i < Nf; i++) res->coef[o*Nf + i] = atom->orb[o].coef[i];
    res->beta[o] = AtomBilinear(atom->orb[o].coef, Mbeta, atom->orb[o].coef, Nf);
    free(atom->orb[o].extra);
    atom->orb[o].extra = NULL;
  }
  // the atom keeps its ground state
  for (o = 0; o < nOrb; o++) for (i = 0; i < Nf; i++) atom->orb[o].coef[i] = C0[o*Nf + i];
  HFFree(&hf);
  free(C0); free(Mbeta); free(Vbeta);
  return error;
}

/* Orthonormalize the functions coef (nOrb x Nf) per symmetry group, in the order of
 * increasing n (method DiracGramSchmidt: the lowest n of each kappa stays as it is; or
 * DiracLowdin), two passes (DiracOrthonormalize). *minEigenvalue receives the smallest
 * eigenvalue of the overlap matrices; below redundancy the functions are (nearly) linearly
 * dependent and 1 is returned. *orthoError receives max |C^T S C - 1| afterwards. */
int AtomOrthonormalizeBasis(DiracAtomType *atom, double *coef, int method, double redundancy, double *minEigenvalue, double *orthoError)
{
  uint32_t gr, o, a, m, i, j, Nf = atom->Nf, *list = (uint32_t *) malloc(atom->nOrb*sizeof(uint32_t));
  double *C = (double *) malloc(atom->nOrb*Nf*sizeof(double)), minEig, err;
  int error = 0;
  *minEigenvalue = 1.;
  *orthoError = 0.;
  for (gr = 0; gr < atom->nGroup && !error; gr++)
  {
    m = 0;
    for (o = 0; o < atom->nOrb; o++)
      if (atom->orb[o].group == gr)
      {
        // insert in the order of n
        for (a = m; a > 0 && atom->orb[list[a-1]].n > atom->orb[o].n; a--) list[a] = list[a-1];
        list[a] = o;
        m++;
      }
    for (j = 0; j < m; j++) for (i = 0; i < Nf; i++) C[i + j*Nf] = coef[list[j]*Nf + i];
    error = DiracOrthonormalize(atom->group[gr].S, Nf, C, m, method, redundancy, &minEig);
    if (minEig < *minEigenvalue) *minEigenvalue = minEig;
    if (error) break;
    for (j = 0; j < m; j++) for (i = 0; i < Nf; i++) coef[list[j]*Nf + i] = C[i + j*Nf];
    err = DiracOrthonormalityError(atom->group[gr].S, Nf, C, m);
    if (err > *orthoError) *orthoError = err;
  }
  free(list); free(C);
  return error;
}
