/*
 * LuaAtomicBasis.cpp
 *
 *     G, F, x, Atom = CalculateAtomicBasisFunctions(Z [, configuration] [, {options}])
 *
 * Atomic radial basis functions (DFT stage A5, design sections 6, 6.1, 10 and 11.1), the
 * drop-in replacement of ReadFPLOBasisFunctions: G (large components P) and F (small
 * components Q) as DKBSplines, orthonormal, ready for NewOperator("Dirac", NF, Z, indices,
 * G, F, Atom.Nucleus) and NewOperator("AU", NF, indices, G, F); x the quadrature points.
 * Schroedinger ({"Relativistic", false}): G as BSplines, F empty.
 *
 * Z: the nuclear charge (a real number). configuration: e.g. "[Ar] 3d8 4s2" or
 * "[Kr] 4d3/2(4) 4d5/2(5)"; an orbital with occupation 0 ("4p0", "5p1/2(0)") is calculated as
 * an empty orbital. Without it, the ground-state configuration of the neutral atom
 * (AtomGroundStateConfiguration, NIST, Z = 1 .. 103). Every orbital of the configuration and
 * all n <= NMax are calculated and returned, ordered by n, l, j (Atom.Labels lists them; a
 * selection is made in Lua). Each orbital (n, kappa) may have its own extra nuclear charge Q
 * and confinement scale P; it then comes from its own generating calculation, in which only
 * its own equation has Z+Q and V_beta = (r/(P rho))^N (AtomBasis.cpp). All are made
 * orthonormal per kappa, Gram-Schmidt from the lowest n upwards.
 * Options (with the nucleus and basis options of DiracHartreeFock); none is required:
 *     Method            "DHF" (default; configuration-average Hartree-Fock), "LDA",
 *                       "LDA-VWN", "RLDA-MV", "PBE", "PBE-VWN"
 *     Relativistic      true (default): Dirac; false: Schroedinger
 *     NMax              {s = 7, p = 6, d = 6, f = 5}, per l, with per-kappa entries such
 *                       as ["p1/2"] = 7 (Schroedinger: per l only)
 *     Q                 {["4p"] = 1, ["5s1/2"] = 2}: the extra nuclear charge of the equation
 *                       of an orbital (default 0); an nl key sets both j (n kappa wins); only
 *                       orbitals that are calculated
 *     P                 {["4d5/2"] = 0.8}: the confinement scale of an orbital (default 1),
 *                       likewise; needs Confinement
 *     Confinement       {Rho = 6, N = 10}: rc = P(n kappa) Rho for every orbital; without it
 *                       there is no confinement
 *     Nucleus           as DiracHartreeFock; a finite nucleus without A or Rrms takes the mass
 *                       number of AtomDefaultMassNumber (the most abundant isotope)
 *     ConfinementVmax   with a confinement and no rMax: rMax = max rc Vmax^(1/N), 1e5 Ha
 *     FourierBessel     {KMax = 30, NK = 1001}: the transforms of the returned orbitals on
 *                       k = 0 .. KMax (NK points), Atom.FourierBessel = {k, P, Q}
 *                       (FourierBessel.cpp; Q empty without relativity)
 *     Orthonormalize    "GramSchmidt" (default, lowest n first) or "Lowdin"
 *     RedundancyTolerance  smallest eigenvalue of an overlap matrix, 1e-8 (error below)
 *     Precision         "Normal" (default) or "High" (design 0.1): EnergyTolerance 1e-13 |E|
 *                       (|E| of the Thomas-Fermi atom), OrbitalTolerance 1e-9 (the eigensolver
 *                       floor of the orbital step: 1.6e-10 for Hg with a finite nucleus, 2e-9
 *                       with a point nucleus), NQuadrature Degree+9, and ConvergenceCheck;
 *                       ResidualTolerance (Kohn-Sham) stays: the DIIS residual stalls at
 *                       1.7e-9 (Re), 7.5e-9 (U, finite nucleus), 8.7e-8 (U, point nucleus), and
 *                       Etot is second order in it. Options given explicitly win
 *     ConvergenceCheck  true: repeat with twice the knot intervals (same rMin, rMax, Degree,
 *                       NQuadrature) and report Atom.ConvergenceCheck = {NIntervals, dEtot,
 *                       dEnergy = max |change of an orbital energy|}
 *     MaxIterations 200, History 8, Mixing (Kohn-Sham ground state 0.3), OrbitalMixing
 *     (every orbital SCF: Hartree-Fock and the generating ones, 0.5), EnergyTolerance 1e-10,
 *     ResidualTolerance (Kohn-Sham potential residual) 1e-8, OrbitalTolerance (orbital step)
 *     1e-8, DensityCutoff 1e-30, MixDensity (default for the GGAs), RayleighQuotient,
 *     OriginBoundary
 * Atom: the ground state (Etot, Ekin, Eenuc, and Edirect, Eexchange for DHF or Ecoul, Exc
 * for Kohn-Sham; Converged, Iterations, Residual, dE; DHF: Gradient, PairGradient;
 * Schroedinger: VirialRatio -(Etot - Ekin)/Ekin), Potential = {rVN, rVH, rVxc (Kohn-Sham)},
 * r times the nuclear, Hartree and exchange-correlation potentials of the ground state, as
 * BSplines through the points x (degree of the orbitals; rVxc the local potential, also for a
 * GGA; f(r, 1) = V + r V', so V' = (f'(r) - f(r)/r)/r), Method, Relativistic, c, Z,
 * NElectron, Configuration, A, rMin, rMax, Nucleus (an options list for NewOperator("Dirac")), Labels,
 * MinOverlapEigenvalue and OrthonormalityError (max |C^T S C - 1|) of the
 * orthonormalization, for a "Function3D" nucleus NuclearMultipoleMoments[i] = {L, M, Value =
 * int r^L S_LM rho d^3r} (L >= 2; S_LM as in NucleusMultipole.cpp), and Orbitals[i] = {Label, n, l, kappa,
 * j, Occupation, Q, P, rc, Energy (eigenvalue of its generating equation), EnergyGround
 * (<a|F|a> of the returned function with the ground-state operator), Vbeta (<beta V_beta>),
 * RExpectation {[-3] (l > 0), [-1], [1], [2]} = <r^k> of the returned function,
 * EtotGenerating, Generated, Converged}; SlaterIntegrals in the format of GetSlaterIntegrals,
 * t["a b a b"] = {[k] = F^k(ab)} and t["a b b a"] = {[k] = G^k(ab)} for the returned orbitals
 * (jj-resolved in Dirac); Dirac also SlaterIntegralsNL, the same for the nl shells, the
 * radial integrals averaged over their j subshells with the weights (2j+1)/(4l+2), and
 * Zeta[nl] = 2 (<F>_{l+1/2} - <F>_{l-1/2})/(2l+1) of the EnergyGround values, for shells
 * returned with both j (design section 10, open point 1).
 */

/* n, l and kappa of an orbital label: Dirac "4d5/2" (kappa = -(l+1) for j = l+1/2, l for
 * j = l-1/2), Schroedinger "4d" (kappa = 0); with allowNL also "4d" in Dirac (kappa = 0, both
 * j). Returns 1 for anything else. */
int DFTParseOrbitalLabel(const char *label, int relativistic, int allowNL, int *n, int *l, int *kappa)
{
  const char *letters = "spdfghi", *where;
  char *end;
  int twoJ;
  *n = (int) strtol(label, &end, 10);
  if (end == label || *n < 1 || *end == '\0') return 1;
  where = strchr(letters, *end);
  if (where == NULL) return 1;
  *l = (int) (where - letters);
  if (*n <= *l) return 1;
  label = end + 1;
  if (*label == '\0')
  {
    *kappa = 0;
    return (relativistic && !allowNL);
  }
  if (!relativistic) return 1;
  twoJ = (int) strtol(label, &end, 10);
  if (end == label || strcmp(end, "/2") != 0) return 1;
  if (twoJ == 2*(*l) + 1) *kappa = -(*l + 1);
  else if (twoJ == 2*(*l) - 1 && *l > 0) *kappa = *l;
  else return 1;
  return 0;
}

// the label of orbital (n, l, kappa)
static void DFTOrbitalLabel(char *label, size_t size, int n, int l, int kappa)
{
  const char *letters = "spdfghi";
  if (kappa == 0) snprintf(label, size, "%d%c", n, letters[l]);
  else snprintf(label, size, "%d%c%d/2", n, letters[l], 2*abs(kappa) - 1);
}

typedef struct { int n, l, kappa; } DFTOrbitalKeyType;

static uint32_t DFTAddOrbital(DFTOrbitalKeyType *list, uint32_t *n, uint32_t max, int nn, int l, int kappa)
{
  uint32_t i;
  for (i = 0; i < *n; i++) if (list[i].n == nn && list[i].l == l && list[i].kappa == kappa) return i;
  if (*n >= max) return max;
  list[*n].n = nn;
  list[*n].l = l;
  list[*n].kappa = kappa;
  return (*n)++;
}

// a number field of the table at the top of the stack, or def if it is absent
// push a BSpline of degree n through the points (x_i, y_i), x increasing
static void DFTLuaPushBSplineOnPoints(lua_State *L, const std::vector<double> &x, const std::vector<double> &y, uint32_t n)
{
  void *memory = lua_newuserdata(L, sizeof(BSpline::TBSpline));
  BSpline::TBSpline *f = new(memory) BSpline::TBSpline;
  luaL_setmetatable(L, "InterpolatingFunction_Type");
  f->Create(x, y, n);
}

static double DFTLuaFieldNumber(lua_State *L, const char *who, const char *table, const char *name, double def)
{
  double v = def;
  lua_getfield(L, -1, name);
  if (!lua_isnil(L, -1))
  {
    if (lua_type(L, -1) != LUA_TNUMBER) { luaL_error(L, "Error in %s: %s.%s must be a number\n", who, table, name); return def; }
    v = lua_tonumber(L, -1);
  }
  lua_pop(L, 1);
  return v;
}

/* R^k(xy;xy) (type 0) or R^k(xy;yx) (type 1) of the orbitals x, y with P, Q on the quadrature
 * points ([x*nPoint + q]); band the Poisson factorizations per k; rho, rho2, Y workspace. */
static double DFTRadialRk(const RadialBasisType *b, double **band, const double *P, const double *Q, uint32_t x, uint32_t y,
                          int type, int k, double *rho, double *rho2, double *Y, int *error)
{
  uint32_t q, nP = b->nPoint;
  for (q = 0; q < nP; q++)
  {
    if (type == 0) { rho[q] = P[x*nP + q]*P[x*nP + q] + Q[x*nP + q]*Q[x*nP + q]; rho2[q] = P[y*nP + q]*P[y*nP + q] + Q[y*nP + q]*Q[y*nP + q]; }
    else rho[q] = rho2[q] = P[x*nP + q]*P[y*nP + q] + Q[x*nP + q]*Q[y*nP + q];
  }
  if (RadialYkOnPoints(b, (uint32_t) k, band[k], rho2, Y)) { *error = 1; return 0.; }
  return HFIntegral(b, rho, Y);
}

/* Slater integrals of the orbitals list[0 .. n-1] (coefficients coef, labels labels), in the
 * format of GetSlaterIntegrals, into the table at index tJJ:
 *   t["a b a b"] = {[k] = F^k(ab)},   t["a b b a"] = {[k] = G^k(ab)}   (a before b in list),
 * F^k(ab) = R^k(ab;ab) = int rho_aa Y^k[rho_bb]/r, G^k(ab) = R^k(ab;ba) = int rho_ab Y^k[rho_ab]/r,
 * for the k the angular factors allow (l_a + l_b + k even; triangle with j_a, j_b, or with
 * l_a, l_b without relativity). Dirac, with tNL > 0: the same for the nl shells into the
 * table at tNL, the radial integrals averaged over the j subshells returned with the
 * weights (2j+1)/(4l+2) (renormalized if a shell is returned with one j only), for the
 * nonrelativistic k: F^k with k = 0, 2, .. 2 min(l, l'), G^k with |l-l'| <= k <= l+l'.
 * Returns 1 if the Poisson equation fails. */
static int DFTSlaterTables(lua_State *L, DiracAtomType *atom, const double *coef, const uint32_t *list, uint32_t n, char labels[][32], int tJJ, int tNL)
{
  RadialBasisType *b = atom->basis;
  uint32_t nP = b->nPoint, Nf = atom->Nf, i, j, k, kMax = 0, nPois = b->nSpline - 1, ld = b->degree + 1;
  double *P = (double *) malloc(n*nP*sizeof(double)), *Q = (double *) malloc(n*nP*sizeof(double)), *dP = (double *) malloc(nP*sizeof(double));
  double *dQ = (double *) malloc(nP*sizeof(double)), *rho = (double *) malloc(nP*sizeof(double)), *rho2 = (double *) malloc(nP*sizeof(double));
  double *Y = (double *) malloc(nP*sizeof(double)), **band;
  int error = 0;
  char key[160];
  for (i = 0; i < n; i++)
  {
    DiracOrbitalType *o = atom->orb + list[i];
    AtomOrbitalOnPoints(b, atom->relativistic, o->kappa, atom->c, coef + list[i]*Nf, P + i*nP, dP, Q + i*nP, dQ);
    if ((uint32_t) (2*o->l + 1) > kMax) kMax = (uint32_t) (2*o->l + 1);
  }
  band = (double **) calloc(kMax + 1, sizeof(double *));
  for (k = 0; k <= kMax && !error; k++)
  {
    band[k] = (double *) malloc(ld*nPois*sizeof(double));
    error = BSplinePoissonSetup(b->knot, b->nKnot, b->degree, b->rq, b->wq, nP, k, band[k], ld);
  }
  for (i = 0; i < n && !error; i++)
    for (j = i; j < n && !error; j++)
    {
      DiracOrbitalType *A = atom->orb + list[i], *B = atom->orb + list[j];
      int ta = atom->relativistic ? 2*abs(A->kappa) - 1 : 2*A->l, tb = atom->relativistic ? 2*abs(B->kappa) - 1 : 2*B->l;
      for (int type = 0; type < 2 && !error; type++)
      {
        if (type == 1 && i == j) break;        // G^k(aa) = F^k(aa)
        int k0 = (type == 0) ? 0 : abs(ta - tb)/2, k1 = (type == 0) ? Min(ta, tb) : (ta + tb)/2;
        snprintf(key, sizeof(key), "%s %s %s %s", labels[i], labels[j], type == 0 ? labels[i] : labels[j], type == 0 ? labels[j] : labels[i]);
        lua_pushstring(L, key);
        lua_newtable(L);
        for (int kk = k0; kk <= k1 && !error; kk++)
        {
          if ((type == 0 && kk % 2 != 0) || (type == 1 && (A->l + B->l + kk) % 2 != 0)) continue;
          lua_pushnumber(L, DFTRadialRk(b, band, P, Q, i, j, type, kk, rho, rho2, Y, &error));
          lua_rawseti(L, -2, kk);
        }
        lua_settable(L, tJJ);
      }
    }
  // the j averages of the nl shells
  for (i = 0; i < n && !error && tNL > 0; i++)
    for (j = i; j < n && !error; j++)
    {
      DiracOrbitalType *A = atom->orb + list[i], *B = atom->orb + list[j];
      uint32_t t;
      int firstA = 1, firstB = 1, same = (A->n == B->n && A->l == B->l);
      for (t = 0; t < i; t++) if (atom->orb[list[t]].n == A->n && atom->orb[list[t]].l == A->l) firstA = 0;
      for (t = 0; t < j; t++) if (atom->orb[list[t]].n == B->n && atom->orb[list[t]].l == B->l) firstB = 0;
      if (!firstA || !firstB || (same && i != j)) continue;
      char la[16], lb[16];
      DFTOrbitalLabel(la, 16, A->n, A->l, 0);
      DFTOrbitalLabel(lb, 16, B->n, B->l, 0);
      for (int type = 0; type < 2 && !error; type++)
      {
        if (type == 1 && same) break;
        int k0 = (type == 0) ? 0 : abs(A->l - B->l), k1 = (type == 0) ? 2*Min(A->l, B->l) : A->l + B->l;
        snprintf(key, sizeof(key), "%s %s %s %s", la, lb, type == 0 ? la : lb, type == 0 ? lb : la);
        lua_pushstring(L, key);
        lua_newtable(L);
        for (int kk = k0; kk <= k1 && !error; kk += 2)
        {
          double sum = 0., wsum = 0.;
          for (uint32_t x = 0; x < n; x++)
            for (uint32_t y = 0; y < n; y++)
            {
              DiracOrbitalType *X = atom->orb + list[x], *W = atom->orb + list[y];
              if (X->n != A->n || X->l != A->l || W->n != B->n || W->l != B->l) continue;
              double w = (2.*abs(X->kappa)/(4.*A->l + 2.))*(2.*abs(W->kappa)/(4.*B->l + 2.));
              sum += w*DFTRadialRk(b, band, P, Q, x, y, type, kk, rho, rho2, Y, &error);
              wsum += w;
            }
          lua_pushnumber(L, sum/wsum);
          lua_rawseti(L, -2, kk);
        }
        lua_settable(L, tNL);
      }
    }
  for (k = 0; k <= kMax; k++) free(band[k]);
  free(band); free(P); free(Q); free(dP); free(dQ); free(rho); free(rho2); free(Y);
  return error;
}

static int LuaCalculateAtomicBasisFunctions(lua_State *L)
{
  struct
  {
    DFTLuaNucleusOptions nuc;
    double Z;
    char *Configuration, *Method, *Orthonormalize;
    bool Relativistic, RayleighQuotient, OriginBoundary, MixDensity, ConvergenceCheck;
    char *Precision;
    int NMax, Qtab, Ptab, Confinement, FourierBessel;
    uint32_t MaxIterations, History;
    double Mixing, OrbitalMixing, EnergyTolerance, ResidualTolerance, OrbitalTolerance, DensityCutoff, RedundancyTolerance, ConfinementVmax;
  } Options;
  LuaOptionsType LuaBasisOptions[] = {
    {"Method", "char*", &(Options.Method)},
    {"Relativistic", "bool", &(Options.Relativistic)},
    {"NMax", "table", &(Options.NMax)},
    {"Q", "table", &(Options.Qtab)},
    {"P", "table", &(Options.Ptab)},
    {"Confinement", "table", &(Options.Confinement)},
    {"FourierBessel", "table", &(Options.FourierBessel)},
    {"ConfinementVmax", "double", &(Options.ConfinementVmax)},
    {"Orthonormalize", "char*", &(Options.Orthonormalize)},
    {"RedundancyTolerance", "double", &(Options.RedundancyTolerance)},
    {"Nucleus", "char*", &(Options.nuc.Nucleus)},
    {"A", "double", &(Options.nuc.A)},
    {"Rrms", "double", &(Options.nuc.Rrms)},
    {"FermiThickness", "double", &(Options.nuc.FermiThickness)},
    {"Density", "function", &(Options.nuc.Density)},
    {"DensityTolerance", "double", &(Options.nuc.DensityTolerance)},
    {"NucleusLmax", "uint32_t", &(Options.nuc.NucleusLmax)},
    {"NucleusNTheta", "uint32_t", &(Options.nuc.NucleusNTheta)},
    {"NucleusNPhi", "uint32_t", &(Options.nuc.NucleusNPhi)},
    {"Degree", "uint32_t", &(Options.nuc.Degree)},
    {"NIntervals", "uint32_t", &(Options.nuc.NIntervals)},
    {"rMin", "double", &(Options.nuc.rMin)},
    {"rMax", "double", &(Options.nuc.rMax)},
    {"NQuadrature", "uint32_t", &(Options.nuc.NQuadrature)},
    {"RayleighQuotient", "bool", &(Options.RayleighQuotient)},
    {"OriginBoundary", "bool", &(Options.OriginBoundary)},
    {"MaxIterations", "uint32_t", &(Options.MaxIterations)},
    {"History", "uint32_t", &(Options.History)},
    {"Mixing", "double", &(Options.Mixing)},
    {"OrbitalMixing", "double", &(Options.OrbitalMixing)},
    {"EnergyTolerance", "double", &(Options.EnergyTolerance)},
    {"ResidualTolerance", "double", &(Options.ResidualTolerance)},
    {"OrbitalTolerance", "double", &(Options.OrbitalTolerance)},
    {"DensityCutoff", "double", &(Options.DensityCutoff)},
    {"MixDensity", "bool", &(Options.MixDensity)},
    {"Precision", "char*", &(Options.Precision)},
    {"ConvergenceCheck", "bool", &(Options.ConvergenceCheck)},
    {NULL, NULL, NULL}};
  const char *who = "CalculateAtomicBasisFunctions";
  const uint32_t maxShell = 256;
  const char *letters = "spdfghi";
  AtomShellType shell[256];
  DFTOrbitalKeyType want[256];
  uint32_t nShell = 0, nWant = 0, len = (uint32_t) lua_gettop(L), i, k, o, s, nRet = 0, ret[256], map[256];
  double c, a0fm, *V, charge, Q[256], Pnk[256], rc[256], Rho = 0., Nconf = 0., rcMax = 0., minEig = 1., orthoError = 0.;
  int method, hartreeFock, error, n, l, kappa, nmaxL[7], nmaxK[7][2], densityIndex = 0, optArg = 0;
  RadialBasisType basis;
  NucleusType nucleus;
  DiracAtomType atom;
  AtomBasisType res;
  AtomSCFOptionsType scf;
  char label[32];

  // CalculateAtomicBasisFunctions(Z), (Z, configuration), (Z, configuration, {options}), (Z, {options})
  if (len < 1 || len > 3) return luaL_error(L, "Error in %s: 1 to 3 arguments expected (Z, configuration, {options}), got %d\n", who, (int) len);
  Options.Z = luaL_checknumber(L, 1);
  Options.Configuration = NULL;
  if (len >= 2 && lua_type(L, 2) == LUA_TSTRING) Options.Configuration = (char *) lua_tostring(L, 2);
  else if (len == 2 && lua_istable(L, 2)) optArg = 2;
  else if (len >= 2 && !lua_isnil(L, 2)) return luaL_error(L, "Error in %s: argument 2 must be the configuration (a string), or the options\n", who);
  if (len == 3)
  {
    if (!lua_istable(L, 3)) return luaL_error(L, "Error in %s: argument 3 must be the options\n", who);
    optArg = 3;
  }
  c = LuaGetPhysicalConstant(L, "cAU");
  a0fm = LuaGetPhysicalConstant(L, "a0")*1e15;
  DFTLuaNucleusDefaults(&Options.nuc, a0fm);
  Options.Method = (char *) "DHF";
  Options.Orthonormalize = (char *) "GramSchmidt";
  Options.Relativistic = true;
  Options.RayleighQuotient = true;
  Options.OriginBoundary = true;
  Options.MixDensity = false;
  Options.NMax = Options.Qtab = Options.Ptab = Options.Confinement = Options.FourierBessel = LUA_NOREF;
  Options.MaxIterations = 200;
  Options.History = 8;
  Options.Mixing = 0.3;
  Options.OrbitalMixing = 0.5;
  Options.EnergyTolerance = 1e-10;
  Options.ResidualTolerance = 1e-8;
  Options.OrbitalTolerance = 1e-8;
  Options.DensityCutoff = 1e-30;
  Options.RedundancyTolerance = 1e-8;
  Options.ConfinementVmax = 1e5;
  Options.Precision = (char *) "Normal";
  Options.ConvergenceCheck = false;
  if (optArg > 0) LuaCheckOptions(L, optArg, LuaBasisOptions);
  // the preset {"Precision", "High"} (design 0.1): tighter defaults, explicit options win
  if (strcicmp(Options.Precision, "High") == 0)
  {
    double scale = 0.7687*pow(fabs(Options.Z), 7./3.);     // |E| of the Thomas-Fermi atom
    if (!(optArg > 0 && DFTLuaHasOption(L, optArg, "EnergyTolerance"))) Options.EnergyTolerance = 1e-13*scale;
    if (!(optArg > 0 && DFTLuaHasOption(L, optArg, "OrbitalTolerance"))) Options.OrbitalTolerance = 1e-9;
    if (!(optArg > 0 && DFTLuaHasOption(L, optArg, "NQuadrature"))) Options.nuc.NQuadrature = Options.nuc.Degree + 9;
    if (!(optArg > 0 && DFTLuaHasOption(L, optArg, "ConvergenceCheck"))) Options.ConvergenceCheck = true;
  }
  else if (strcicmp(Options.Precision, "Normal") != 0)
    return luaL_error(L, "Error in %s: Precision \"%s\" is not Normal or High\n", who, Options.Precision);
  if (!(Options.Z > 0.)) return luaL_error(L, "Error in %s: the nuclear charge Z = %f must be positive\n", who, Options.Z);
  if (Options.Configuration == NULL)
  {
    Options.Configuration = (char *) AtomGroundStateConfiguration(Options.Z);
    if (Options.Configuration == NULL)
      return luaL_error(L, "Error in %s: no ground-state configuration known for Z = %f (Z = 1 .. 103); give the configuration\n", who, Options.Z);
  }
  if (strcicmp(Options.Method, "DHF") == 0) { hartreeFock = 1; method = XCLDAVWN; }
  else
  {
    hartreeFock = 0;
    method = XCMethodFromName(Options.Method);
    if (method < 0) return luaL_error(L, "Error in %s: unknown Method \"%s\" (DHF, LDA, LDA-VWN, RLDA-MV, PBE, PBE-VWN)\n", who, Options.Method);
    if (!(optArg > 0 && DFTLuaHasOption(L, optArg, "MixDensity"))) Options.MixDensity = XCIsGGA(method) ? true : false;
  }
  if (AtomParseConfiguration(Options.Configuration, maxShell, &nShell, shell))
    return luaL_error(L, "Error in %s: can not read the configuration \"%s\"\n", who, Options.Configuration);
  if (!Options.Relativistic)
    for (s = 0; s < nShell; s++)
      if (shell[s].kappa != 0) return luaL_error(L, "Error in %s: n kappa occupations need {\"Relativistic\", true}\n", who);

  // the orbitals: the configuration
  for (s = 0; s < nShell; s++)
  {
    if (shell[s].kappa != 0) DFTAddOrbital(want, &nWant, maxShell, shell[s].n, shell[s].l, shell[s].kappa);
    else if (!Options.Relativistic) DFTAddOrbital(want, &nWant, maxShell, shell[s].n, shell[s].l, 0);
    else
    {
      if (shell[s].l > 0) DFTAddOrbital(want, &nWant, maxShell, shell[s].n, shell[s].l, shell[s].l);
      DFTAddOrbital(want, &nWant, maxShell, shell[s].n, shell[s].l, -(shell[s].l + 1));
    }
  }
  // NMax
  for (l = 0; l < 7; l++) { nmaxL[l] = 0; nmaxK[l][0] = nmaxK[l][1] = -1; }
  if (Options.NMax != LUA_NOREF)
  {
    lua_rawgeti(L, LUA_REGISTRYINDEX, Options.NMax);
    lua_pushnil(L);
    while (lua_next(L, -2) != 0)
    {
      const char *key = (lua_type(L, -2) == LUA_TSTRING) ? lua_tostring(L, -2) : NULL;
      const char *where = (key != NULL && key[0] != '\0') ? strchr(letters, key[0]) : NULL;
      if (where == NULL || lua_type(L, -1) != LUA_TNUMBER)
        return luaL_error(L, "Error in %s: NMax entries are l = n (s = 5) or kappa = n ([\"p1/2\"] = 5)\n", who);
      l = (int) (where - letters);
      n = (int) lua_tointeger(L, -1);
      if (key[1] == '\0') nmaxL[l] = n;
      else
      {
        char full[32];
        snprintf(full, sizeof(full), "%d%s", l + 1, key);
        if (!Options.Relativistic || DFTParseOrbitalLabel(full, 1, 0, &n, &l, &kappa))
          return luaL_error(L, "Error in %s: NMax entry \"%s\" is not a kappa such as p1/2 (Dirac only)\n", who, key);
        nmaxK[l][kappa > 0 ? 0 : 1] = (int) lua_tointeger(L, -1);
      }
      lua_pop(L, 1);
    }
    lua_pop(L, 1);
    luaL_unref(L, LUA_REGISTRYINDEX, Options.NMax);
    for (l = 0; l < 7; l++)
      for (int t = 0; t < (Options.Relativistic ? 2 : 1); t++)
      {
        int nm = (nmaxK[l][t] >= 0) ? nmaxK[l][t] : nmaxL[l];
        if (Options.Relativistic && t == 0 && l == 0) continue;     // no kappa = 0
        kappa = !Options.Relativistic ? 0 : (t == 0 ? l : -(l + 1));
        for (n = l + 1; n <= nm; n++)
          if (DFTAddOrbital(want, &nWant, maxShell, n, l, kappa) == maxShell) return luaL_error(L, "Error in %s: more than %d orbitals\n", who, (int) maxShell);
      }
  }
  // Q and P per orbital (options "Q" and "P"): nl keys first, then n kappa keys (they win);
  // only orbitals that are calculated (the configuration and NMax)
  for (i = 0; i < maxShell; i++) { Q[i] = 0.; Pnk[i] = 1.; rc[i] = 0.; }
  for (int which = 0; which < 2; which++)
  {
    int ref = (which == 0) ? Options.Qtab : Options.Ptab;
    const char *name = (which == 0) ? "Q" : "P";
    if (ref == LUA_NOREF) continue;
    lua_rawgeti(L, LUA_REGISTRYINDEX, ref);
    for (int pass = 0; pass < 2; pass++)
    {
      lua_pushnil(L);
      while (lua_next(L, -2) != 0)
      {
        if (lua_type(L, -2) != LUA_TSTRING || lua_type(L, -1) != LUA_TNUMBER
            || DFTParseOrbitalLabel(lua_tostring(L, -2), Options.Relativistic, 1, &n, &l, &kappa))
          return luaL_error(L, "Error in %s: the entries of %s are [\"4d5/2\"] = number or [\"4d\"] = number\n", who, name);
        if ((pass == 0) == (kappa == 0 || !Options.Relativistic))
        {
          double v = lua_tonumber(L, -1);
          if (which == 1 && !(v > 0.)) return luaL_error(L, "Error in %s: P[\"%s\"] must be positive\n", who, lua_tostring(L, -2));
          for (int t = 0; t < 2; t++)
          {
            int kap = kappa;
            if (kappa == 0 && Options.Relativistic) { if (t == 0 && l == 0) continue; kap = (t == 0) ? l : -(l + 1); }
            else if (t == 1) break;
            for (k = 0; k < nWant; k++) if (want[k].n == n && want[k].l == l && want[k].kappa == kap) break;
            if (k == nWant)
              return luaL_error(L, "Error in %s: %s[\"%s\"]: the orbital is not calculated (add it to the configuration, with occupation 0 if empty, or to NMax)\n",
                                who, name, lua_tostring(L, -2));
            if (which == 0) Q[k] = v; else Pnk[k] = v;
          }
        }
        lua_pop(L, 1);
      }
    }
    lua_pop(L, 1);
    luaL_unref(L, LUA_REGISTRYINDEX, ref);
  }
  if (Options.Confinement != LUA_NOREF)
  {
    lua_rawgeti(L, LUA_REGISTRYINDEX, Options.Confinement);
    Rho = DFTLuaFieldNumber(L, who, "Confinement", "Rho", 0.);
    Nconf = DFTLuaFieldNumber(L, who, "Confinement", "N", 0.);
    lua_getfield(L, -1, "P");
    int hasAtomP = !lua_isnil(L, -1);
    lua_pop(L, 2);
    luaL_unref(L, LUA_REGISTRYINDEX, Options.Confinement);
    if (hasAtomP) return luaL_error(L, "Error in %s: Confinement takes Rho and N; the radius is P Rho with P per orbital (option P, default 1)\n", who);
    if (!(Rho > 0.) || !(Nconf > 0.)) return luaL_error(L, "Error in %s: Confinement needs Rho > 0 and N > 0\n", who);
    for (k = 0; k < nWant; k++) { rc[k] = Pnk[k]*Rho; if (rc[k] > rcMax) rcMax = rc[k]; }
  }
  else if (Options.Ptab != LUA_NOREF) return luaL_error(L, "Error in %s: the confinement scale P needs the option Confinement {Rho = .., N = ..}\n", who);

  // the shells of the atom: the configuration and an empty one for every other orbital
  for (k = 0; k < nWant; k++)
  {
    int found = 0;
    for (s = 0; s < nShell && !found; s++)
      if (shell[s].n == want[k].n && shell[s].l == want[k].l && (shell[s].kappa == 0 || shell[s].kappa == want[k].kappa)) found = 1;
    if (found) continue;
    if (nShell >= maxShell) return luaL_error(L, "Error in %s: more than %d shells\n", who, (int) maxShell);
    shell[nShell].n = want[k].n;
    shell[nShell].l = want[k].l;
    shell[nShell].kappa = want[k].kappa;
    shell[nShell].q = 0.;
    nShell++;
  }

  // basis, nucleus, atom
  if (Options.nuc.rMax <= 0.) Options.nuc.rMax = (rcMax > 0.) ? ConfinementRMax(rcMax, Nconf, Options.ConfinementVmax) : 40.;
  // a finite nucleus without A or Rrms: the mass number of AtomDefaultMassNumber
  if ((strcicmp(Options.nuc.Nucleus, "Gaussian") == 0 || strcicmp(Options.nuc.Nucleus, "Uniform") == 0 || strcicmp(Options.nuc.Nucleus, "Fermi") == 0)
      && !(Options.nuc.A > 0.) && !(Options.nuc.Rrms > 0.))
    Options.nuc.A = AtomDefaultMassNumber(Options.Z);
  if (Options.nuc.Density != LUA_NOREF)
  {
    lua_rawgeti(L, LUA_REGISTRYINDEX, Options.nuc.Density);     // kept for Atom.Nucleus
    densityIndex = lua_gettop(L);
  }
  double moments[1024];
  if (strcicmp(Options.nuc.Nucleus, "Function3D") == 0)
  {
    if (NucleusMultipoleCount(Options.nuc.NucleusLmax) > 1024) return luaL_error(L, "Error in %s: NucleusLmax %d is too large\n", who, (int) Options.nuc.NucleusLmax);
    Options.nuc.multipoleMoment = moments;
  }
  DFTLuaNucleusBasis(L, who, Options.Z, &Options.nuc, a0fm, &nucleus, &basis, &V, &charge);
  error = AtomSetup(&atom, &basis, V, Options.Z, Options.Relativistic ? 1 : 0, c, method, Options.OriginBoundary ? 1 : 0,
                    Options.DensityCutoff, nShell, shell);
  free(V);
  if (error) { AtomFree(&atom); RadialBasisFree(&basis); return luaL_error(L, "Error in %s: the atom could not be set up\n", who); }
  // the wanted orbitals among those of the atom
  for (k = 0; k < nWant; k++)
  {
    map[k] = atom.nOrb;
    for (o = 0; o < atom.nOrb; o++)
      if (atom.orb[o].n == want[k].n && atom.orb[o].l == want[k].l && atom.orb[o].kappa == want[k].kappa) map[k] = o;
    if (map[k] == atom.nOrb) { AtomFree(&atom); RadialBasisFree(&basis); return luaL_error(L, "Error in %s: an orbital is missing (internal)\n", who); }
  }
  {
    double Qo[256], rco[256];
    for (o = 0; o < atom.nOrb; o++) { Qo[o] = 0.; rco[o] = 0.; }
    for (k = 0; k < nWant; k++) { Qo[map[k]] = Q[k]; rco[map[k]] = rc[k]; }
    scf.maxIterations = Options.MaxIterations;
    scf.history = Options.History;
    scf.mixing = hartreeFock ? Options.OrbitalMixing : Options.Mixing;
    scf.energyTolerance = Options.EnergyTolerance;
    scf.residualTolerance = Options.ResidualTolerance;
    scf.orbitalTolerance = Options.OrbitalTolerance;
    scf.rayleigh = Options.RayleighQuotient ? 1 : 0;
    scf.mixDensity = Options.MixDensity ? 1 : 0;
    error = AtomBasisFunctions(&atom, hartreeFock, Qo, rco, Nconf, &scf, &res);
  }
  if (!error)
  {
    int m = (strcicmp(Options.Orthonormalize, "Lowdin") == 0) ? DiracLowdin : DiracGramSchmidt;
    if (strcicmp(Options.Orthonormalize, "Lowdin") != 0 && strcicmp(Options.Orthonormalize, "GramSchmidt") != 0) error = 2;
    else if (AtomOrthonormalizeBasis(&atom, res.coef, m, Options.RedundancyTolerance, &minEig, &orthoError)) error = 3;
  }
  if (error)
  {
    AtomBasisFree(&res, atom.nOrb);
    AtomFree(&atom);
    RadialBasisFree(&basis);
    if (error == 2) return luaL_error(L, "Error in %s: Orthonormalize \"%s\" is not GramSchmidt or Lowdin\n", who, Options.Orthonormalize);
    if (error == 3) return luaL_error(L, "Error in %s: the basis functions are (nearly) linearly dependent: smallest overlap eigenvalue %f < RedundancyTolerance %f\n",
                                      who, minEig, Options.RedundancyTolerance);
    return luaL_error(L, "Error in %s: a self-consistent field or eigenproblem failed\n", who);
  }
  if (!res.groundConverged)
    luaL_warning(L, "Warning in %s: the ground state is not converged in %d iterations (dE = %f, residual = %f)\n", who, (int) res.iterations, res.dE, res.residual);
  for (o = 0; o < atom.nOrb; o++)
  {
    res.energyGround[o] = AtomBilinear(res.coef + o*atom.Nf, res.Fground[o], res.coef + o*atom.Nf, atom.Nf);
    if (!res.converged[o])
      luaL_warning(L, "Warning in %s: the generating calculation of orbital %d%c is not converged\n", who, atom.orb[o].n, letters[atom.orb[o].l]);
  }
  // the returned orbitals: all, ordered by n, l, j
  {
    nRet = 0;
    for (k = 0; k < nWant; k++)
    {
      for (i = nRet; i > 0; i--)
      {
        DFTOrbitalKeyType *a = want + ret[i-1], *b = want + k;
        int before = (b->n < a->n) || (b->n == a->n && (b->l < a->l || (b->l == a->l && abs(b->kappa) < abs(a->kappa))));
        if (!before) break;
        ret[i] = ret[i-1];
      }
      ret[i] = k;
      nRet++;
    }
  }

  // G, F, x
  {
    uint32_t list[256];
    for (o = 0; o < atom.nOrb; o++) for (i = 0; i < atom.Nf; i++) atom.orb[o].coef[i] = res.coef[o*atom.Nf + i];
    for (k = 0; k < nRet; k++) list[k] = map[ret[k]];
    DFTLuaPushRadialFunctions(L, &atom, list, nRet);
    LuaPushArray_double(L, basis.rq, (int) basis.nPoint);
  }
  // Atom
  lua_newtable(L);
  lua_pushnumber(L, res.Etot);   lua_setfield(L, -2, "Etot");
  lua_pushnumber(L, res.Ekin);   lua_setfield(L, -2, "Ekin");
  lua_pushnumber(L, res.Eenuc);  lua_setfield(L, -2, "Eenuc");
  lua_pushnumber(L, res.Edirect);   lua_setfield(L, -2, hartreeFock ? "Edirect" : "Ecoul");
  lua_pushnumber(L, res.Eexchange); lua_setfield(L, -2, hartreeFock ? "Eexchange" : "Exc");
  lua_pushboolean(L, res.groundConverged); lua_setfield(L, -2, "Converged");
  lua_pushinteger(L, res.iterations); lua_setfield(L, -2, "Iterations");
  lua_pushnumber(L, res.residual); lua_setfield(L, -2, "Residual");
  lua_pushnumber(L, res.dE);     lua_setfield(L, -2, "dE");
  lua_pushstring(L, Options.Method); lua_setfield(L, -2, "Method");
  lua_pushboolean(L, Options.Relativistic); lua_setfield(L, -2, "Relativistic");
  lua_pushnumber(L, c);          lua_setfield(L, -2, "c");
  lua_pushnumber(L, Options.Z);  lua_setfield(L, -2, "Z");
  lua_pushstring(L, Options.Configuration); lua_setfield(L, -2, "Configuration");
  if (Options.nuc.A > 0.) { lua_pushnumber(L, Options.nuc.A); lua_setfield(L, -2, "A"); }
  lua_pushnumber(L, atom.nElectron); lua_setfield(L, -2, "NElectron");
  lua_pushnumber(L, Options.nuc.rMin); lua_setfield(L, -2, "rMin");
  lua_pushnumber(L, Options.nuc.rMax); lua_setfield(L, -2, "rMax");
  lua_pushnumber(L, minEig);     lua_setfield(L, -2, "MinOverlapEigenvalue");
  lua_pushnumber(L, orthoError); lua_setfield(L, -2, "OrthonormalityError");
  if (hartreeFock)
  {
    lua_pushnumber(L, res.gradient);     lua_setfield(L, -2, "Gradient");
    lua_pushnumber(L, res.pairGradient); lua_setfield(L, -2, "PairGradient");
  }
  if (!Options.Relativistic) { lua_pushnumber(L, -(res.Etot - res.Ekin)/res.Ekin); lua_setfield(L, -2, "VirialRatio"); }
  /* potentials of the ground state: r V_N, r V_H (of the ground-state density) and r V_xc
   * (Kohn-Sham), each a BSpline through its values on the points x, of the degree of the
   * orbitals (the same points and degree, so the same knots). The GGA potential of the solver
   * is in weak form, v + g d/dr (AtomPotential); the local potential is v - g', g' the
   * derivative of the BSpline of g through the points. */
  {
    uint32_t nP = basis.nPoint;
    double *VH = (double *) malloc(nP*sizeof(double)), *v = (double *) malloc(nP*sizeof(double)), *g = (double *) malloc(nP*sizeof(double)), Ec, Ex;
    std::vector<double> r(basis.rq, basis.rq + nP), rV(nP), dg(nP, 0.);
    AtomDensity(&atom);
    RadialYkOnPoints(&basis, 0, atom.poisson, atom.f, VH);
    lua_newtable(L);
    for (i = 0; i < nP; i++) rV[i] = r[i]*atom.VN[i];
    DFTLuaPushBSplineOnPoints(L, r, rV, basis.degree); lua_setfield(L, -2, "rVN");
    for (i = 0; i < nP; i++) rV[i] = r[i]*VH[i];
    DFTLuaPushBSplineOnPoints(L, r, rV, basis.degree); lua_setfield(L, -2, "rVH");
    if (!hartreeFock && !AtomPotential(&atom, atom.f, atom.df, v, g, &Ec, &Ex))
    {
      if (XCIsGGA(atom.method))
      {
        BSpline::TBSpline gs;
        gs.Create(r, std::vector<double>(g, g + nP), basis.degree);
        dg = gs(r, 1);
      }
      for (i = 0; i < nP; i++) rV[i] = r[i]*(v[i] - VH[i] - dg[i]);
      DFTLuaPushBSplineOnPoints(L, r, rV, basis.degree); lua_setfield(L, -2, "rVxc");
    }
    lua_setfield(L, -2, "Potential");
    free(VH); free(v); free(g);
  }
  // Nucleus, as options for NewOperator("Dirac")
  lua_newtable(L);
  lua_newtable(L); lua_pushstring(L, "Nucleus"); lua_rawseti(L, -2, 1);
  lua_pushstring(L, nucleus.model == NucleusPoint ? "Point" : nucleus.model == NucleusGaussian ? "Gaussian" : nucleus.model == NucleusFermi ? "Fermi"
                  : nucleus.model == NucleusUniform ? "Uniform" : nucleus.model == NucleusFunction ? "Function" : "Function3D");
  lua_rawseti(L, -2, 2); lua_rawseti(L, -2, 1);
  k = 2;
  if (nucleus.model == NucleusGaussian || nucleus.model == NucleusFermi || nucleus.model == NucleusUniform)
  {
    lua_newtable(L); lua_pushstring(L, "Rrms"); lua_rawseti(L, -2, 1); lua_pushnumber(L, nucleus.Rrms); lua_rawseti(L, -2, 2); lua_rawseti(L, -2, k++);
  }
  if (nucleus.model == NucleusFermi)
  {
    lua_newtable(L); lua_pushstring(L, "FermiThickness"); lua_rawseti(L, -2, 1); lua_pushnumber(L, Options.nuc.FermiThickness); lua_rawseti(L, -2, 2); lua_rawseti(L, -2, k++);
  }
  if ((nucleus.model == NucleusFunction || nucleus.model == NucleusFunction3D) && densityIndex > 0)
  {
    lua_newtable(L); lua_pushstring(L, "Density"); lua_rawseti(L, -2, 1); lua_pushvalue(L, densityIndex); lua_rawseti(L, -2, 2); lua_rawseti(L, -2, k++);
  }
  if (nucleus.model == NucleusFunction3D)
  {
    const char *nm[3] = {"NucleusLmax", "NucleusNTheta", "NucleusNPhi"};
    uint32_t nv[3] = {Options.nuc.NucleusLmax, Options.nuc.NucleusNTheta, Options.nuc.NucleusNPhi};
    for (int t = 0; t < 3; t++)
    {
      lua_newtable(L); lua_pushstring(L, nm[t]); lua_rawseti(L, -2, 1); lua_pushinteger(L, nv[t]); lua_rawseti(L, -2, 2); lua_rawseti(L, -2, k++);
    }
  }
  lua_setfield(L, -2, "Nucleus");
  // Function3D: the multipole moments q_LM = int r^L S_LM rho d^3r of L >= 2
  if (nucleus.model == NucleusFunction3D)
  {
    uint32_t nLM = NucleusMultipoleCount(Options.nuc.NucleusLmax), t = 0;
    lua_newtable(L);
    for (i = 0; i < nLM; i++)
    {
      int Lq, Mq;
      NucleusMultipoleLM(i, &Lq, &Mq);
      if (Lq < 2) continue;
      lua_newtable(L);
      lua_pushinteger(L, Lq); lua_setfield(L, -2, "L");
      lua_pushinteger(L, Mq); lua_setfield(L, -2, "M");
      lua_pushnumber(L, moments[i]); lua_setfield(L, -2, "Value");
      lua_rawseti(L, -2, ++t);
    }
    lua_setfield(L, -2, "NuclearMultipoleMoments");
  }
  // Labels and Orbitals
  lua_createtable(L, nRet, 0);
  for (k = 0; k < nRet; k++)
  {
    DFTOrbitalLabel(label, sizeof(label), want[ret[k]].n, want[ret[k]].l, want[ret[k]].kappa);
    lua_pushstring(L, label);
    lua_rawseti(L, -2, k+1);
  }
  lua_setfield(L, -2, "Labels");
  lua_createtable(L, nRet, 0);
  for (k = 0; k < nRet; k++)
  {
    uint32_t w = ret[k];
    o = map[w];
    DiracOrbitalType *orb = atom.orb + o;
    lua_newtable(L);
    DFTOrbitalLabel(label, sizeof(label), orb->n, orb->l, orb->kappa);
    lua_pushstring(L, label);         lua_setfield(L, -2, "Label");
    lua_pushinteger(L, orb->n);       lua_setfield(L, -2, "n");
    lua_pushinteger(L, orb->l);       lua_setfield(L, -2, "l");
    if (Options.Relativistic)
    {
      lua_pushinteger(L, orb->kappa); lua_setfield(L, -2, "kappa");
      lua_pushnumber(L, abs(orb->kappa) - 0.5); lua_setfield(L, -2, "j");
    }
    lua_pushnumber(L, orb->q);        lua_setfield(L, -2, "Occupation");
    lua_pushnumber(L, Q[w]);          lua_setfield(L, -2, "Q");
    lua_pushnumber(L, Pnk[w]);        lua_setfield(L, -2, "P");
    lua_pushnumber(L, rc[w]);         lua_setfield(L, -2, "rc");
    lua_pushnumber(L, res.energy[o]); lua_setfield(L, -2, "Energy");
    lua_pushnumber(L, res.energyGround[o]); lua_setfield(L, -2, "EnergyGround");
    lua_pushnumber(L, res.beta[o]);   lua_setfield(L, -2, "Vbeta");
    lua_pushnumber(L, res.Egenerating[o]); lua_setfield(L, -2, "EtotGenerating");
    lua_pushboolean(L, res.generated[o]); lua_setfield(L, -2, "Generated");
    lua_pushboolean(L, res.converged[o]); lua_setfield(L, -2, "Converged");
    {
      // <r^k> = int (P^2 + Q^2) r^k dr of the returned function, k = -3 (l > 0), -1, 1, 2
      uint32_t nP = basis.nPoint, q;
      double *P = (double *) malloc(nP*sizeof(double)), *dP = (double *) malloc(nP*sizeof(double));
      double *Qs = (double *) malloc(nP*sizeof(double)), *dQ = (double *) malloc(nP*sizeof(double));
      int ks[4] = {-3, -1, 1, 2};
      AtomOrbitalOnPoints(&basis, atom.relativistic, orb->kappa, c, res.coef + o*atom.Nf, P, dP, Qs, dQ);
      lua_newtable(L);
      for (int t = 0; t < 4; t++)
      {
        DFTReal sum = 0., comp = 0.;
        if (ks[t] == -3 && orb->l == 0) continue;
        for (q = 0; q < nP; q++) NeumaierAdd(&sum, &comp, basis.wq[q]*(P[q]*P[q] + Qs[q]*Qs[q])*pow(basis.rq[q], ks[t]));
        lua_pushnumber(L, sum + comp);
        lua_rawseti(L, -2, ks[t]);
      }
      lua_setfield(L, -2, "RExpectation");
      free(P); free(dP); free(Qs); free(dQ);
    }
    lua_rawseti(L, -2, k+1);
  }
  lua_setfield(L, -2, "Orbitals");
  // Fourier-Bessel transforms of the returned orbitals on k = 0 .. KMax (NK points)
  if (Options.FourierBessel != LUA_NOREF)
  {
    double kMax, *kk, *Pk, *Qk;
    uint32_t nk;
    lua_rawgeti(L, LUA_REGISTRYINDEX, Options.FourierBessel);
    kMax = DFTLuaFieldNumber(L, who, "FourierBessel", "KMax", 0.);
    nk = (uint32_t) DFTLuaFieldNumber(L, who, "FourierBessel", "NK", 1001.);
    lua_pop(L, 1);
    luaL_unref(L, LUA_REGISTRYINDEX, Options.FourierBessel);
    if (!(kMax > 0.) || nk < 2)
    {
      AtomBasisFree(&res, atom.nOrb); AtomFree(&atom); RadialBasisFree(&basis);
      return luaL_error(L, "Error in %s: FourierBessel needs KMax > 0 and NK >= 2\n", who);
    }
    kk = (double *) malloc(nk*sizeof(double));
    Pk = (double *) malloc(nk*sizeof(double));
    Qk = (double *) malloc(nk*sizeof(double));
    for (i = 0; i < nk; i++) kk[i] = kMax*i/(nk - 1.);
    lua_newtable(L);
    LuaPushArray_double(L, kk, (int) nk);
    lua_setfield(L, -2, "k");
    lua_createtable(L, nRet, 0);
    lua_createtable(L, nRet, 0);
    for (k = 0; k < nRet && !error; k++)
    {
      o = map[ret[k]];
      error = FourierBesselTransform(&basis, atom.relativistic, atom.orb[o].kappa, atom.orb[o].l, c, atom.orb[o].coef, nk, kk, Pk, Qk);
      LuaPushArray_double(L, Pk, (int) nk);
      lua_rawseti(L, -3, k + 1);
      if (atom.relativistic) { LuaPushArray_double(L, Qk, (int) nk); lua_rawseti(L, -2, k + 1); }
    }
    lua_setfield(L, -3, "Q");
    lua_setfield(L, -2, "P");
    lua_setfield(L, -2, "FourierBessel");
    free(kk); free(Pk); free(Qk);
    if (error)
    {
      AtomBasisFree(&res, atom.nOrb); AtomFree(&atom); RadialBasisFree(&basis);
      return luaL_error(L, "Error in %s: the Fourier-Bessel transform failed\n", who);
    }
  }
  // Slater integrals of the returned orbitals, jj-resolved (Dirac) or nl (Schroedinger), the
  // j averages, and zeta_nl = 2 (<F>_{l+1/2} - <F>_{l-1/2})/(2l+1) of EnergyGround
  {
    uint32_t list[256];
    static char labs[256][32];
    int tAtom = lua_gettop(L);
    for (k = 0; k < nRet; k++)
    {
      list[k] = map[ret[k]];
      DFTOrbitalLabel(labs[k], 32, atom.orb[list[k]].n, atom.orb[list[k]].l, atom.orb[list[k]].kappa);
    }
    lua_newtable(L);
    lua_newtable(L);
    if (DFTSlaterTables(L, &atom, res.coef, list, nRet, labs, tAtom + 1, Options.Relativistic ? tAtom + 2 : 0))
    {
      AtomBasisFree(&res, atom.nOrb); AtomFree(&atom); RadialBasisFree(&basis);
      return luaL_error(L, "Error in %s: the Slater integrals failed (Poisson equation)\n", who);
    }
    if (Options.Relativistic) lua_setfield(L, tAtom, "SlaterIntegralsNL"); else lua_pop(L, 1);
    lua_setfield(L, tAtom, "SlaterIntegrals");
    if (Options.Relativistic)
    {
      lua_newtable(L);
      for (k = 0; k < nRet; k++)
      {
        DiracOrbitalType *A = atom.orb + list[k];
        if (A->l == 0 || A->kappa < 0) continue;                 // from the j = l-1/2 member
        for (uint32_t t = 0; t < nRet; t++)
        {
          DiracOrbitalType *X = atom.orb + list[t];
          if (X->n == A->n && X->l == A->l && X->kappa < 0)
          {
            char la[16];
            DFTOrbitalLabel(la, 16, A->n, A->l, 0);
            lua_pushnumber(L, 2.*(res.energyGround[list[t]] - res.energyGround[list[k]])/(2.*A->l + 1.));
            lua_setfield(L, -2, la);
          }
        }
      }
      lua_setfield(L, tAtom, "Zeta");
    }
  }
  AtomBasisFree(&res, atom.nOrb);
  AtomFree(&atom);
  RadialBasisFree(&basis);
  if (densityIndex > 0) lua_remove(L, densityIndex);   // the reference was released by DFTLuaNuclearPotential
  // ConvergenceCheck: the same calculation with twice the knot intervals (the same rMin, rMax,
  // Degree, NQuadrature); Atom.ConvergenceCheck = {NIntervals, dEtot, dEnergy = max |d eps|}
  if (Options.ConvergenceCheck)
  {
    int tAtom = lua_gettop(L), nOpt = (optArg > 0) ? (int) lua_objlen(L, optArg) : 0, m = 0;
    lua_pushcfunction(L, LuaCalculateAtomicBasisFunctions);
    lua_pushnumber(L, Options.Z);
    lua_pushstring(L, Options.Configuration);
    lua_newtable(L);
    for (int t = 1; t <= nOpt; t++)
    {
      lua_rawgeti(L, optArg, t);
      lua_rawgeti(L, -1, 1);
      const char *name = lua_tostring(L, -1);
      int skip = (name != NULL) && (strcicmp(name, "NIntervals") == 0 || strcicmp(name, "ConvergenceCheck") == 0 || strcicmp(name, "FourierBessel") == 0
                                    || strcicmp(name, "rMin") == 0 || strcicmp(name, "rMax") == 0 || strcicmp(name, "NQuadrature") == 0);
      lua_pop(L, 1);
      if (skip) lua_pop(L, 1); else lua_rawseti(L, -2, ++m);
    }
    const char *nm[5] = {"NIntervals", "ConvergenceCheck", "rMin", "rMax", "NQuadrature"};
    double nv[5] = {2.*Options.nuc.NIntervals, 0., Options.nuc.rMin, Options.nuc.rMax, (double) Options.nuc.NQuadrature};
    for (int t = 0; t < 5; t++)
    {
      lua_newtable(L);
      lua_pushstring(L, nm[t]); lua_rawseti(L, -2, 1);
      if (t == 1) lua_pushboolean(L, 0); else lua_pushnumber(L, nv[t]);
      lua_rawseti(L, -2, 2);
      lua_rawseti(L, -2, ++m);
    }
    lua_call(L, 3, 4);
    // the refined Atom at the top, the original at tAtom
    double dE, dEps = 0.;
    lua_getfield(L, -1, "Etot"); lua_getfield(L, tAtom, "Etot");
    dE = lua_tonumber(L, -2) - lua_tonumber(L, -1);
    lua_pop(L, 2);
    lua_getfield(L, -1, "Orbitals"); lua_getfield(L, tAtom, "Orbitals");
    for (k = 1; k <= (uint32_t) lua_objlen(L, -1); k++)
    {
      lua_rawgeti(L, -2, k); lua_getfield(L, -1, "Energy");
      lua_rawgeti(L, -3, k); lua_getfield(L, -1, "Energy");
      dEps = fmax(dEps, fabs(lua_tonumber(L, -3) - lua_tonumber(L, -1)));
      lua_pop(L, 4);
    }
    lua_pop(L, 2 + 4);
    lua_newtable(L);
    lua_pushinteger(L, 2*Options.nuc.NIntervals); lua_setfield(L, -2, "NIntervals");
    lua_pushnumber(L, dE);   lua_setfield(L, -2, "dEtot");
    lua_pushnumber(L, dEps); lua_setfield(L, -2, "dEnergy");
    lua_setfield(L, tAtom, "ConvergenceCheck");
  }
  return 4;
}
