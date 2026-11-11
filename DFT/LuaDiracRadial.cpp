/*
 * LuaDiracRadial.cpp
 *
 * Temporary Lua entry point for testing the radial Dirac solver of DFT stage A1:
 *
 *     E, G, F, Info = DiracRadialSolve(Z, kappa, {options})
 *
 * The final interface is CalculateAtomicBasisFunctions (task A5, design section 11.1).
 *
 * One electron with relativistic quantum number kappa in the field of a nucleus of charge
 * Z (a real number). Returns the NStates lowest electronic eigenvalues E (E[i] belongs to
 * n = l(kappa) + i), the large components G = P and the small components F = Q of these
 * states as interpolating functions of type "DKBSpline" (exact evaluation, 0 outside
 * [0, rMax]), and a table Info. Options (lengths in bohr):
 *     Nucleus          "Point", "Gaussian" (default), "Fermi" or "Function"
 *     A, Rrms          mass number, or the root-mean-square radius; Rrms = (0.836 A^(1/3)
 *                      + 0.570) fm if only A is given (Visscher and Dyall)
 *     FermiThickness   skin thickness t of the Fermi nucleus, default 2.30 fm
 *     Density          Lua function rho(r) (charge per volume) for Nucleus "Function";
 *                      sampled here on the quadrature points (no call back from C)
 *     DensityTolerance |4 pi int r^2 rho - Z| / Z above which a warning is printed,
 *                      default 1e-10; the density is always rescaled to Z
 *     Degree           B-spline degree, default 8 (order 9)
 *     NIntervals       knot intervals, default 120
 *     rMin, rMax       first interior knot and the end of the knots; default rMin =
 *                      Rrms/50 for a finite nucleus and 1e-6/Z for a point nucleus,
 *                      rMax = 40, or from the confinement (ConfinementVmax)
 *     NQuadrature      Gauss points per knot interval, default Degree+5, at most 127
 *     NStates          number of states returned, default 3
 *     RayleighQuotient true (default): the eigenvalues are the Rayleigh quotients of the
 *                      eigenvectors (DiracRadialEigen); false: those of the eigensolver
 *                      (dsygv), always in Info.EEigensolver
 *     OriginBoundary   true (default): the DKB functions on B_1 that are nonzero at r = 0
 *                      are left out for kappa != -1 (DiracActiveFunctions), which removes
 *                      the spurious state of a point nucleus for kappa = +1; false: all 2n
 *     ConfinementRadius rc > 0 switches on the confinement V_conf = (r/rc)^N (design
 *                      section 6), default none
 *     ConfinementPower N, default 8
 *     ConfinementType  "Scalar" (default): Lorentz scalar, beta V_conf; "Vector": added
 *                      to V, for comparison only (negative-continuum intruders)
 *     ConfinementVmax  with a confinement and no rMax given, rMax = rc Vmax^(1/N), where
 *                      V_conf reaches Vmax; default 1e5 Ha (rMax = 4.2 rc for N = 8)
 * Info.NBasis is the number of DKB functions in use, Info.NLarge and Info.NSmall those of
 * type u^(1) and u^(2); the negative continuum has NSmall states.
 * Info.BetaVbeta[i] is <beta V_conf> = int V_conf (P^2 - Q^2) dr for "Scalar" and
 * int V_conf (P^2 + Q^2) dr for "Vector"; Info.NormQ[i] = int Q^2 dr.
 * c is PhysicalConstants.cAU.value, 1 bohr = PhysicalConstants.a0.value m.
 */

/* The options of the DFT Lua entries that fix the nucleus and the radial basis. */
typedef struct
{
  char *Nucleus;
  double A, Rrms, FermiThickness, DensityTolerance, rMin, rMax;
  uint32_t Degree, NIntervals, NQuadrature;
  int Density;
} DFTLuaNucleusOptions;

int DFTLuaNucleus(lua_State *L, const char *who, double Z, DFTLuaNucleusOptions *o, double a0fm, NucleusType *nucleus);
int DFTLuaNuclearPotential(lua_State *L, const char *who, double Z, DFTLuaNucleusOptions *o, const NucleusType *nucleus,
                           RadialBasisType *basis, double **V, double *charge);

void DFTLuaNucleusDefaults(DFTLuaNucleusOptions *o, double a0fm)
{
  o->Nucleus = (char *) "Gaussian";
  o->A = 0.;
  o->Rrms = 0.;
  o->FermiThickness = 2.30/a0fm;
  o->DensityTolerance = 1e-10;
  o->rMin = 0.;
  o->rMax = 0.;
  o->Degree = 8;
  o->NIntervals = 120;
  o->NQuadrature = 0;
  o->Density = LUA_NOREF;
}

/* The nucleus, the radial basis and the nuclear potential *V (on the quadrature points,
 * allocated here) from the options of a DFT Lua entry (who names it in the messages);
 * o->rMin and o->rMax (rMax must be set) receive the values used, *charge the norm of a
 * "Function" density before rescaling. Raises a Lua error on failure. */
int DFTLuaNucleusBasis(lua_State *L, const char *who, double Z, DFTLuaNucleusOptions *o, double a0fm,
                       NucleusType *nucleus, RadialBasisType *basis, double **V, double *charge)
{
  if (o->NQuadrature == 0) o->NQuadrature = o->Degree + 5;
  if (o->NQuadrature > DiracMaxQuadraturePoints)
    return luaL_error(L, "Error in %s: NQuadrature = %d, at most %d points per knot interval (the 34-digit Gauss-Legendre table)\n",
                      who, (int) o->NQuadrature, DiracMaxQuadraturePoints);
  DFTLuaNucleus(L, who, Z, o, a0fm, nucleus);
  if (o->rMin <= 0.) o->rMin = (nucleus->Rrms > 0.) ? nucleus->Rrms/50. : 1e-6/fabs(Z);

  if (RadialBasisCreate(o->Degree, o->NIntervals, o->rMin, o->rMax, o->NQuadrature, basis))
  {
    if (o->Density != LUA_NOREF) luaL_unref(L, LUA_REGISTRYINDEX, o->Density);
    return luaL_error(L, "Error in %s: no B-spline basis for Degree %d, NIntervals %d, rMin %f, rMax %f\n",
                      who, (int) o->Degree, (int) o->NIntervals, o->rMin, o->rMax);
  }
  return DFTLuaNuclearPotential(L, who, Z, o, nucleus, basis, V, charge);
}

/* The nucleus model of a DFT Lua entry from its options (who names it in the messages).
 * Raises a Lua error for an unknown or incomplete model. */
int DFTLuaNucleus(lua_State *L, const char *who, double Z, DFTLuaNucleusOptions *o, double a0fm, NucleusType *nucleus)
{
  nucleus->Z = Z;
  nucleus->A = o->A;
  nucleus->Rrms = 0.;
  nucleus->cFermi = 0.;
  nucleus->aFermi = 0.;
  if (strcicmp(o->Nucleus, "Point") == 0) nucleus->model = NucleusPoint;
  else if (strcicmp(o->Nucleus, "Gaussian") == 0) nucleus->model = NucleusGaussian;
  else if (strcicmp(o->Nucleus, "Fermi") == 0) nucleus->model = NucleusFermi;
  else if (strcicmp(o->Nucleus, "Function") == 0) nucleus->model = NucleusFunction;
  else return luaL_error(L, "Error in %s: unknown Nucleus \"%s\" (Point, Gaussian, Fermi or Function)\n", who, o->Nucleus);
  if ((nucleus->model == NucleusFunction) != (o->Density != LUA_NOREF))
  {
    if (o->Density != LUA_NOREF) luaL_unref(L, LUA_REGISTRYINDEX, o->Density);
    return luaL_error(L, "Error in %s: the option Density goes with Nucleus \"Function\", and only with it\n", who);
  }
  if (nucleus->model == NucleusGaussian || nucleus->model == NucleusFermi)
  {
    if (o->Rrms > 0.) nucleus->Rrms = o->Rrms;
    else if (o->A > 0.) nucleus->Rrms = (0.836*pow(o->A, 1./3.) + 0.570)/a0fm;
    else return luaL_error(L, "Error in %s: a %s nucleus needs the option A or Rrms\n", who, o->Nucleus);
  }
  if (nucleus->model == NucleusFermi && NucleusFermiParameters(nucleus->Rrms, o->FermiThickness, &nucleus->cFermi, &nucleus->aFermi))
    return luaL_error(L, "Error in %s: Fermi nucleus with Rrms = %f bohr is too small for the skin thickness %f bohr (Rrms^2 <= 7/5 pi^2 a^2)\n",
                      who, nucleus->Rrms, o->FermiThickness);
  return 0;
}

/* The nuclear potential *V (allocated here) on the quadrature points of basis, for the
 * nucleus of DFTLuaNucleus; a "Function" density (o->Density, released here) is sampled
 * on them and rescaled to Z, *charge receives its norm before that. On a Lua error the
 * basis is freed. */
int DFTLuaNuclearPotential(lua_State *L, const char *who, double Z, DFTLuaNucleusOptions *o, const NucleusType *nucleus,
                           RadialBasisType *basis, double **V, double *charge)
{
  double *rho = NULL;
  uint32_t q;
  *charge = 0.;
  // sample the nuclear density of a Lua function on the quadrature points, normalize to Z
  if (nucleus->model == NucleusFunction)
  {
    rho = (double *) malloc(basis->nPoint*sizeof(double));
    for (q = 0; q < basis->nPoint; q++)
    {
      lua_rawgeti(L, LUA_REGISTRYINDEX, o->Density);
      lua_pushnumber(L, basis->rq[q]);
      lua_call(L, 1, 1);
      if (lua_type(L, -1) != LUA_TNUMBER)
      {
        free(rho);
        RadialBasisFree(basis);
        return luaL_error(L, "Error in %s: Density(%f) is not a number\n", who, basis->rq[q]);
      }
      rho[q] = lua_tonumber(L, -1);
      lua_pop(L, 1);
    }
    luaL_unref(L, LUA_REGISTRYINDEX, o->Density);
    *charge = NucleusCharge(basis, rho);
    if (!(*charge > 0.))
    {
      free(rho);
      RadialBasisFree(basis);
      return luaL_error(L, "Error in %s: the nuclear charge 4 pi int r^2 rho dr = %f is not positive\n", who, *charge);
    }
    if (fabs(*charge - Z) > o->DensityTolerance*fabs(Z))
      luaL_warning(L, "Warning in %s: 4 pi int r^2 rho dr = %f differs from Z = %f by %f relative; the density is rescaled to Z\n",
                   who, *charge, Z, fabs(*charge - Z)/fabs(Z));
    for (q = 0; q < basis->nPoint; q++) rho[q] *= Z/(*charge);
  }
  *V = (double *) malloc(basis->nPoint*sizeof(double));
  if (NucleusPotential(nucleus, basis, rho, *V))
  {
    free(*V);
    free(rho);
    RadialBasisFree(basis);
    return luaL_error(L, "Error in %s: the nuclear potential failed\n", who);
  }
  free(rho);
  return 0;
}

static int LuaDiracRadialSolve(lua_State *L)
{
  struct
  {
    DFTLuaNucleusOptions nuc;
    uint32_t NStates;
    bool RayleighQuotient, OriginBoundary;
    double ConfinementRadius, ConfinementPower, ConfinementVmax;
    char *ConfinementType;
  } Options;
  LuaOptionsType LuaDiracRadialOptions[] = {
    {"Nucleus", "char*", &(Options.nuc.Nucleus)},
    {"A", "double", &(Options.nuc.A)},
    {"Rrms", "double", &(Options.nuc.Rrms)},
    {"FermiThickness", "double", &(Options.nuc.FermiThickness)},
    {"Density", "function", &(Options.nuc.Density)},
    {"DensityTolerance", "double", &(Options.nuc.DensityTolerance)},
    {"Degree", "uint32_t", &(Options.nuc.Degree)},
    {"NIntervals", "uint32_t", &(Options.nuc.NIntervals)},
    {"rMin", "double", &(Options.nuc.rMin)},
    {"rMax", "double", &(Options.nuc.rMax)},
    {"NQuadrature", "uint32_t", &(Options.nuc.NQuadrature)},
    {"NStates", "uint32_t", &(Options.NStates)},
    {"RayleighQuotient", "bool", &(Options.RayleighQuotient)},
    {"OriginBoundary", "bool", &(Options.OriginBoundary)},
    {"ConfinementRadius", "double", &(Options.ConfinementRadius)},
    {"ConfinementPower", "double", &(Options.ConfinementPower)},
    {"ConfinementVmax", "double", &(Options.ConfinementVmax)},
    {"ConfinementType", "char*", &(Options.ConfinementType)},
    {NULL, NULL, NULL}};
  RadialBasisType basis;
  NucleusType nucleus;
  double c, a0fm, Z, *V, *E, *Eeig2, *C, *res, *coef, charge = 0., *Vconf = NULL, *ones;
  int vectorConfinement = 0;
  uint32_t len = (uint32_t) lua_gettop(L), N, first, i, q, nStates, nE, nLarge, *active;
  int kappa, l;

  if (len < 2 || len > 3) return luaL_error(L, "Error in DiracRadialSolve: 2 or 3 arguments expected (Z, kappa, {options}), got %d\n", (int) len);
  Z = luaL_checknumber(L, 1);
  kappa = (int) luaL_checkinteger(L, 2);
  if (kappa == 0) return luaL_error(L, "Error in DiracRadialSolve: kappa = 0\n");
  l = DiracKappaL(kappa);
  c = LuaGetPhysicalConstant(L, "cAU");
  a0fm = LuaGetPhysicalConstant(L, "a0")*1e15;   // 1 bohr in fm

  DFTLuaNucleusDefaults(&Options.nuc, a0fm);
  Options.NStates = 3;
  Options.RayleighQuotient = true;
  Options.OriginBoundary = true;
  Options.ConfinementRadius = 0.;
  Options.ConfinementPower = 8.;
  Options.ConfinementVmax = 1e5;
  Options.ConfinementType = (char *) "Scalar";
  if (len == 3) LuaCheckOptions(L, 3, LuaDiracRadialOptions);
  if (strcicmp(Options.ConfinementType, "Vector") == 0) vectorConfinement = 1;
  else if (strcicmp(Options.ConfinementType, "Scalar") != 0)
    return luaL_error(L, "Error in DiracRadialSolve: unknown ConfinementType \"%s\" (Scalar or Vector)\n", Options.ConfinementType);
  if (Options.ConfinementRadius < 0. || !(Options.ConfinementPower > 0.) || !(Options.ConfinementVmax > 0.))
    return luaL_error(L, "Error in DiracRadialSolve: ConfinementRadius, ConfinementPower and ConfinementVmax must be positive\n");
  if (Options.nuc.rMax <= 0.)
    Options.nuc.rMax = (Options.ConfinementRadius > 0.) ? ConfinementRMax(Options.ConfinementRadius, Options.ConfinementPower, Options.ConfinementVmax) : 40.;

  DFTLuaNucleusBasis(L, "DiracRadialSolve", Z, &Options.nuc, a0fm, &nucleus, &basis, &V, &charge);
  N = 2*(basis.nSpline - 2);

  if (Options.ConfinementRadius > 0.)
  {
    Vconf = (double *) malloc(basis.nPoint*sizeof(double));
    ConfinementPotential(&basis, Options.ConfinementRadius, Options.ConfinementPower, Vconf);
    if (vectorConfinement) for (q = 0; q < basis.nPoint; q++) V[q] += Vconf[q];
  }

  E = (double *) malloc(N*sizeof(double));
  Eeig2 = (double *) malloc(N*sizeof(double));
  res = (double *) malloc(N*sizeof(double));
  C = (double *) malloc(N*N*sizeof(double));
  if (DiracRadialEigen(&basis, kappa, c, V, vectorConfinement ? NULL : Vconf, Options.RayleighQuotient ? 1 : 0,
                       Options.OriginBoundary ? 1 : 0, E, C, &nE, res, Eeig2))
  {
    free(V); free(Vconf); free(E); free(Eeig2); free(res); free(C);
    RadialBasisFree(&basis);
    return luaL_error(L, "Error in DiracRadialSolve: the generalized eigenvalue problem (dsygv) failed\n");
  }
  first = DiracElectronicFirst(E, nE, c);
  nStates = Options.NStates;
  if (first + nStates > nE) nStates = nE - first;
  active = (uint32_t *) malloc(N*sizeof(uint32_t));
  DiracActiveFunctions(&basis, kappa, Options.OriginBoundary ? 1 : 0, active, &nE, &nLarge);
  free(active);

  // E, G, F
  lua_createtable(L, nStates, 0);
  for (i = 0; i < nStates; i++)
  {
    lua_pushnumber(L, E[first + i]);
    lua_rawseti(L, -2, i+1);
  }
  coef = (double *) malloc(N*sizeof(double));
  for (int comp = 0; comp < 2; comp++)
  {
    lua_createtable(L, nStates, 0);
    for (i = 0; i < nStates; i++)
    {
      for (q = 0; q < N; q++) coef[q] = C[q + (first + i)*N];
      DiracSignConvention(&basis, kappa, c, coef);
      void *memory = lua_newuserdata(L, sizeof(Interpolation::DKBSpline));
      Interpolation::DKBSpline *f = new(memory) Interpolation::DKBSpline;
      luaL_setmetatable(L, "InterpolatingFunction_Type");
      f->Create(basis.knot, basis.nKnot, basis.degree, kappa, (uint32_t) comp, c, coef, basis.nQ);
      lua_rawseti(L, -2, i+1);
    }
  }
  free(coef);

  // Info
  lua_newtable(L);
  lua_pushnumber(L, c);                      lua_setfield(L, -2, "c");
  lua_pushinteger(L, kappa);                 lua_setfield(L, -2, "kappa");
  lua_pushinteger(L, l);                     lua_setfield(L, -2, "l");
  lua_pushinteger(L, first);                 lua_setfield(L, -2, "NNegativeContinuum");
  lua_pushinteger(L, nE);                    lua_setfield(L, -2, "NBasis");
  lua_pushinteger(L, nLarge);                lua_setfield(L, -2, "NLarge");
  lua_pushinteger(L, nE - nLarge);           lua_setfield(L, -2, "NSmall");
  lua_pushinteger(L, basis.nKnot);           lua_setfield(L, -2, "NKnots");
  lua_pushnumber(L, Options.nuc.rMin);       lua_setfield(L, -2, "rMin");
  lua_pushnumber(L, Options.nuc.rMax);       lua_setfield(L, -2, "rMax");
  lua_createtable(L, nStates, 0);
  for (i = 0; i < nStates; i++) { lua_pushnumber(L, Eeig2[first + i]); lua_rawseti(L, -2, i+1); }
  lua_setfield(L, -2, "EEigensolver");
  lua_createtable(L, nStates, 0);
  for (i = 0; i < nStates; i++) { lua_pushnumber(L, res[first + i]); lua_rawseti(L, -2, i+1); }
  lua_setfield(L, -2, "Residual");
  ones = (double *) malloc(basis.nPoint*sizeof(double));
  for (q = 0; q < basis.nPoint; q++) ones[q] = 1.;
  lua_createtable(L, nStates, 0);
  for (i = 0; i < nStates; i++)
  {
    // int Q^2 = (int (P^2+Q^2) - int (P^2-Q^2))/2
    lua_pushnumber(L, 0.5*(DiracExpectation(&basis, kappa, c, C + (first + i)*N, ones, 1) - DiracExpectation(&basis, kappa, c, C + (first + i)*N, ones, -1)));
    lua_rawseti(L, -2, i+1);
  }
  lua_setfield(L, -2, "NormQ");
  free(ones);
  if (Vconf != NULL)
  {
    lua_createtable(L, nStates, 0);
    for (i = 0; i < nStates; i++)
    {
      lua_pushnumber(L, DiracExpectation(&basis, kappa, c, C + (first + i)*N, Vconf, vectorConfinement ? 1 : -1));
      lua_rawseti(L, -2, i+1);
    }
    lua_setfield(L, -2, "BetaVbeta");
    lua_newtable(L);
    lua_pushnumber(L, Options.ConfinementRadius);   lua_setfield(L, -2, "Radius");
    lua_pushnumber(L, Options.ConfinementPower);    lua_setfield(L, -2, "Power");
    lua_pushnumber(L, Options.ConfinementVmax);     lua_setfield(L, -2, "Vmax");
    lua_pushstring(L, vectorConfinement ? "Vector" : "Scalar"); lua_setfield(L, -2, "Type");
    lua_setfield(L, -2, "Confinement");
  }
  LuaPushArray_double(L, basis.rq, (int) basis.nPoint);
  lua_setfield(L, -2, "r");
  LuaPushArray_double(L, basis.wq, (int) basis.nPoint);
  lua_setfield(L, -2, "w");
  LuaPushArray_double(L, V, (int) basis.nPoint);
  lua_setfield(L, -2, "V");
  lua_newtable(L);
  lua_pushstring(L, Options.nuc.Nucleus);    lua_setfield(L, -2, "Model");
  lua_pushnumber(L, nucleus.Z);              lua_setfield(L, -2, "Z");
  lua_pushnumber(L, nucleus.A);              lua_setfield(L, -2, "A");
  lua_pushnumber(L, nucleus.Rrms);           lua_setfield(L, -2, "Rrms");
  lua_pushnumber(L, nucleus.cFermi);         lua_setfield(L, -2, "cFermi");
  lua_pushnumber(L, nucleus.aFermi);         lua_setfield(L, -2, "aFermi");
  if (nucleus.model == NucleusFunction) { lua_pushnumber(L, charge); lua_setfield(L, -2, "DensityCharge"); }
  lua_setfield(L, -2, "Nucleus");

  free(V); free(Vconf); free(E); free(Eeig2); free(res); free(C);
  RadialBasisFree(&basis);
  return 4;
}

/*
 * Temporary Lua entry point for testing the radial Poisson solver (DFT stage A2):
 *
 *     Y = RadialPoisson(rho, k, {options})
 *
 * Y^k(r) = r int_0^rMax rho(r') r<^k / r>^(k+1) dr' for the radial density rho(r) (a Lua
 * function, per dr, e.g. P_a P_b + Q_a Q_b; it should vanish beyond rMax), from the
 * Galerkin solution of Y'' - k(k+1) Y/r^2 = -(2k+1) rho/r (design section 7). rho is
 * sampled on the quadrature points. Y is returned exactly, as an interpolating function of
 * type "BSpline" on the exponential knots. Options: Degree (default 8), NIntervals (120),
 * rMin (1e-6), rMax (40), NQuadrature (Degree+5, at most 127).
 */
static int LuaRadialPoisson(lua_State *L)
{
  struct
  {
    double rMin, rMax;
    uint32_t Degree, NIntervals, NQuadrature;
  } Options;
  LuaOptionsType LuaRadialPoissonOptions[] = {
    {"Degree", "uint32_t", &(Options.Degree)},
    {"NIntervals", "uint32_t", &(Options.NIntervals)},
    {"rMin", "double", &(Options.rMin)},
    {"rMax", "double", &(Options.rMax)},
    {"NQuadrature", "uint32_t", &(Options.NQuadrature)},
    {NULL, NULL, NULL}};
  RadialBasisType basis;
  uint32_t len = (uint32_t) lua_gettop(L), q, i, k, ld;
  double *rho, *band, *Y;
  int error;

  if (len < 2 || len > 3) return luaL_error(L, "Error in RadialPoisson: 2 or 3 arguments expected (rho, k, {options}), got %d\n", (int) len);
  luaL_checktype(L, 1, LUA_TFUNCTION);
  if (luaL_checkinteger(L, 2) < 0) return luaL_error(L, "Error in RadialPoisson: k < 0\n");
  k = (uint32_t) luaL_checkinteger(L, 2);
  Options.Degree = 8;
  Options.NIntervals = 120;
  Options.rMin = 1e-6;
  Options.rMax = 40.;
  Options.NQuadrature = 0;
  if (len == 3) LuaCheckOptions(L, 3, LuaRadialPoissonOptions);
  if (Options.NQuadrature == 0) Options.NQuadrature = Options.Degree + 5;
  if (Options.NQuadrature > DiracMaxQuadraturePoints)
    return luaL_error(L, "Error in RadialPoisson: NQuadrature = %d, at most %d points per knot interval\n", (int) Options.NQuadrature, DiracMaxQuadraturePoints);
  if (RadialBasisCreate(Options.Degree, Options.NIntervals, Options.rMin, Options.rMax, Options.NQuadrature, &basis))
    return luaL_error(L, "Error in RadialPoisson: no B-spline basis for Degree %d, NIntervals %d, rMin %f, rMax %f\n",
                      (int) Options.Degree, (int) Options.NIntervals, Options.rMin, Options.rMax);
  rho = (double *) malloc(basis.nPoint*sizeof(double));
  for (q = 0; q < basis.nPoint; q++)
  {
    lua_pushvalue(L, 1);
    lua_pushnumber(L, basis.rq[q]);
    lua_call(L, 1, 1);
    if (lua_type(L, -1) != LUA_TNUMBER)
    {
      free(rho);
      RadialBasisFree(&basis);
      return luaL_error(L, "Error in RadialPoisson: rho(%f) is not a number\n", basis.rq[q]);
    }
    rho[q] = lua_tonumber(L, -1);
    lua_pop(L, 1);
  }
  ld = basis.degree + 1;
  band = (double *) malloc(ld*(basis.nSpline - 1)*sizeof(double));
  Y = (double *) malloc(basis.nSpline*sizeof(double));
  error = BSplinePoissonSetup(basis.knot, basis.nKnot, basis.degree, basis.rq, basis.wq, basis.nPoint, k, band, ld);
  if (!error) error = BSplinePoissonSolve(basis.knot, basis.nKnot, basis.degree, basis.rq, basis.wq, basis.nPoint, k, band, ld, rho, Y);
  if (error)
  {
    free(rho); free(band); free(Y);
    RadialBasisFree(&basis);
    return luaL_error(L, "Error in RadialPoisson: the Poisson solver failed\n");
  }
  // the breakpoints, for x (Integrate is exact on them)
  std::vector<double> x, knots(basis.knot, basis.knot + basis.nKnot), coef(Y, Y + basis.nSpline);
  for (i = 0; i < basis.nKnot; i++) if (i == 0 || basis.knot[i] > basis.knot[i-1]) x.push_back(basis.knot[i]);
  void *memory = lua_newuserdata(L, sizeof(BSpline::TBSpline));
  new(memory) BSpline::TBSpline(x, knots, coef, basis.degree);
  luaL_setmetatable(L, "InterpolatingFunction_Type");
  free(rho); free(band); free(Y);
  RadialBasisFree(&basis);
  return 1;
}

/*
 * Temporary Lua entry point for testing the on-site orthonormalization (DFT stage A2,
 * design section 6.1):
 *
 *     G2, F2, Info = DiracOrthonormalize(G, F, {options})
 *
 * G and F are tables of the large and small components (type "DKBSpline") of m spinors
 * with the same kappa, knots and c, e.g. from DiracRadialSolve with different potentials
 * on the same knots. Returns the orthonormalized spinors (in the overlap of the DKB basis,
 * exact in coefficient space) and Info with the smallest eigenvalue of the overlap matrix
 * O = C^T S C (MinOverlapEigenvalue) and max |C^T S C - 1| before and after
 * (OrthonormalityErrorBefore, OrthonormalityError). Options:
 *     Method                   "GramSchmidt" (default; ordered, in the order of G) or
 *                              "Lowdin" (symmetric)
 *     RedundancyThreshold      error if an eigenvalue of O is below it (the orbitals are
 *                              nearly linearly dependent); default 1e-8
 *     OrthonormalityTolerance  error if max |C^T S C - 1| exceeds it afterwards;
 *                              default 1e-13
 */
static int LuaDiracOrthonormalize(lua_State *L)
{
  struct
  {
    char *Method;
    double RedundancyThreshold, OrthonormalityTolerance;
  } Options;
  LuaOptionsType LuaDiracOrthonormalizeOptions[] = {
    {"Method", "char*", &(Options.Method)},
    {"RedundancyThreshold", "double", &(Options.RedundancyThreshold)},
    {"OrthonormalityTolerance", "double", &(Options.OrthonormalityTolerance)},
    {NULL, NULL, NULL}};
  RadialBasisType basis;
  uint32_t len = (uint32_t) lua_gettop(L), m, i, j, n, N;
  Interpolation::DKBSpline *g0, *f;
  std::vector<Interpolation::DKBSpline *> Gs, Fs;
  double *S, *C, minEig = 0., before, after;
  int method;

  if (len < 2 || len > 3) return luaL_error(L, "Error in DiracOrthonormalize: 2 or 3 arguments expected (G, F, {options}), got %d\n", (int) len);
  luaL_checktype(L, 1, LUA_TTABLE);
  luaL_checktype(L, 2, LUA_TTABLE);
  Options.Method = (char *) "GramSchmidt";
  Options.RedundancyThreshold = 1e-8;
  Options.OrthonormalityTolerance = 1e-13;
  if (len == 3) LuaCheckOptions(L, 3, LuaDiracOrthonormalizeOptions);
  if (strcicmp(Options.Method, "GramSchmidt") == 0) method = DiracGramSchmidt;
  else if (strcicmp(Options.Method, "Lowdin") == 0) method = DiracLowdin;
  else return luaL_error(L, "Error in DiracOrthonormalize: unknown Method \"%s\" (GramSchmidt or Lowdin)\n", Options.Method);
  m = (uint32_t) lua_objlen(L, 1);
  if (m == 0 || m != (uint32_t) lua_objlen(L, 2)) return luaL_error(L, "Error in DiracOrthonormalize: G and F must hold the same number (> 0) of functions\n");
  for (i = 0; i < m; i++)
  {
    lua_rawgeti(L, 1, i+1);
    lua_check(L, -1, f);
    Gs.push_back(f);
    lua_pop(L, 1);
    lua_rawgeti(L, 2, i+1);
    lua_check(L, -1, f);
    Fs.push_back(f);
    lua_pop(L, 1);
  }
  // all on one knot set, one kappa and one c; G[i] and F[i] two components of one spinor
  g0 = Gs[0];
  for (i = 0; i < m; i++)
  {
    for (j = 0; j < 2; j++)
    {
      f = (j == 0) ? Gs[i] : Fs[i];
      if (f->component != j)
        return luaL_error(L, "Error in DiracOrthonormalize: %s[%d] is not a %s component\n", j ? "F" : "G", (int) i+1, j ? "small" : "large");
      if (f->kappa != g0->kappa || f->degree != g0->degree || f->knots != g0->knots || f->speedOfLight != g0->speedOfLight || f->nQ != g0->nQ)
        return luaL_error(L, "Error in DiracOrthonormalize: all functions must have the same kappa, knots, c and quadrature (%s[%d] differs from G[1])\n",
                          j ? "F" : "G", (int) i+1);
    }
    if (Gs[i]->a != Fs[i]->a || Gs[i]->b != Fs[i]->b)
      return luaL_error(L, "Error in DiracOrthonormalize: G[%d] and F[%d] are not the two components of one spinor\n", (int) i+1, (int) i+1);
  }
  if (RadialBasisFromKnots(g0->knots.data(), (uint32_t) g0->knots.size(), g0->degree, g0->nQ, &basis))
    return luaL_error(L, "Error in DiracOrthonormalize: invalid knots\n");
  n = basis.nSpline - 2;
  N = 2*n;
  S = (double *) malloc(N*N*sizeof(double));
  C = (double *) malloc(N*m*sizeof(double));
  if (DiracDKBOverlap(&basis, g0->kappa, g0->speedOfLight, S))
  {
    free(S); free(C);
    RadialBasisFree(&basis);
    return luaL_error(L, "Error in DiracOrthonormalize: the DKB overlap matrix failed\n");
  }
  for (j = 0; j < m; j++)
    for (i = 0; i < n; i++)
    {
      C[i + j*N] = Gs[j]->a[i+1];
      C[n + i + j*N] = Gs[j]->b[i+1];
    }
  before = DiracOrthonormalityError(S, N, C, m);
  if (DiracOrthonormalize(S, N, C, m, method, Options.RedundancyThreshold, &minEig))
  {
    free(S); free(C);
    RadialBasisFree(&basis);
    return luaL_error(L, "Error in DiracOrthonormalize: the overlap matrix has the eigenvalue %f < RedundancyThreshold %f; the orbitals are (nearly) linearly dependent\n",
                      minEig, Options.RedundancyThreshold);
  }
  after = DiracOrthonormalityError(S, N, C, m);
  if (!(after <= Options.OrthonormalityTolerance))
  {
    free(S); free(C);
    RadialBasisFree(&basis);
    return luaL_error(L, "Error in DiracOrthonormalize: max |C^T S C - 1| = %f after orthonormalization exceeds OrthonormalityTolerance %f\n",
                      after, Options.OrthonormalityTolerance);
  }
  for (j = 0; j < 2; j++)
  {
    lua_createtable(L, m, 0);
    for (i = 0; i < m; i++)
    {
      void *memory = lua_newuserdata(L, sizeof(Interpolation::DKBSpline));
      Interpolation::DKBSpline *h = new(memory) Interpolation::DKBSpline;
      luaL_setmetatable(L, "InterpolatingFunction_Type");
      h->Create(basis.knot, basis.nKnot, basis.degree, g0->kappa, j, g0->speedOfLight, C + i*N, basis.nQ);
      lua_rawseti(L, -2, i+1);
    }
  }
  lua_newtable(L);
  lua_pushnumber(L, minEig);   lua_setfield(L, -2, "MinOverlapEigenvalue");
  lua_pushnumber(L, before);   lua_setfield(L, -2, "OrthonormalityErrorBefore");
  lua_pushnumber(L, after);    lua_setfield(L, -2, "OrthonormalityError");
  free(S); free(C);
  RadialBasisFree(&basis);
  return 3;
}
