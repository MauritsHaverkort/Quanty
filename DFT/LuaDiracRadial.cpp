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
 *                      rMax = 40
 *     NQuadrature      Gauss points per knot interval, default Degree+5, at most 127
 *     NStates          number of states returned, default 3
 *     RayleighQuotient true (default): the eigenvalues are the Rayleigh quotients of the
 *                      eigenvectors (DiracRadialEigen); false: those of dsygvd
 * c is PhysicalConstants.cAU.value, 1 bohr = PhysicalConstants.a0.value m.
 */

static int LuaDiracRadialSolve(lua_State *L)
{
  struct
  {
    char *Nucleus;
    double A, Rrms, FermiThickness, DensityTolerance, rMin, rMax;
    uint32_t Degree, NIntervals, NQuadrature, NStates;
    bool RayleighQuotient;
    int Density;
  } Options;
  LuaOptionsType LuaDiracRadialOptions[] = {
    {"Nucleus", "char*", &(Options.Nucleus)},
    {"A", "double", &(Options.A)},
    {"Rrms", "double", &(Options.Rrms)},
    {"FermiThickness", "double", &(Options.FermiThickness)},
    {"Density", "function", &(Options.Density)},
    {"DensityTolerance", "double", &(Options.DensityTolerance)},
    {"Degree", "uint32_t", &(Options.Degree)},
    {"NIntervals", "uint32_t", &(Options.NIntervals)},
    {"rMin", "double", &(Options.rMin)},
    {"rMax", "double", &(Options.rMax)},
    {"NQuadrature", "uint32_t", &(Options.NQuadrature)},
    {"NStates", "uint32_t", &(Options.NStates)},
    {"RayleighQuotient", "bool", &(Options.RayleighQuotient)},
    {NULL, NULL, NULL}};
  RadialBasisType basis;
  NucleusType nucleus;
  double c, a0fm, Z, *V, *rho = NULL, *E, *Edsygvd, *C, *res, *coef, charge = 0.;
  uint32_t len = (uint32_t) lua_gettop(L), N, first, i, q, nStates;
  int kappa, l;

  if (len < 2 || len > 3) return luaL_error(L, "Error in DiracRadialSolve: 2 or 3 arguments expected (Z, kappa, {options}), got %d\n", len);
  Z = luaL_checknumber(L, 1);
  kappa = (int) luaL_checkinteger(L, 2);
  if (kappa == 0) return luaL_error(L, "Error in DiracRadialSolve: kappa = 0\n");
  l = DiracKappaL(kappa);
  c = LuaGetPhysicalConstant(L, "cAU");
  a0fm = LuaGetPhysicalConstant(L, "a0")*1e15;   // 1 bohr in fm

  Options.Nucleus = (char *) "Gaussian";
  Options.A = 0.;
  Options.Rrms = 0.;
  Options.FermiThickness = 2.30/a0fm;
  Options.DensityTolerance = 1e-10;
  Options.Degree = 8;
  Options.NIntervals = 120;
  Options.rMin = 0.;
  Options.rMax = 40.;
  Options.NQuadrature = 0;
  Options.NStates = 3;
  Options.RayleighQuotient = true;
  Options.Density = LUA_NOREF;
  if (len == 3) LuaCheckOptions(L, 3, LuaDiracRadialOptions);
  if (Options.NQuadrature == 0) Options.NQuadrature = Options.Degree + 5;
  if (Options.NQuadrature > DiracMaxQuadraturePoints)
    return luaL_error(L, "Error in DiracRadialSolve: NQuadrature = %d, at most %d points per knot interval (the 34-digit Gauss-Legendre table)\n",
                      (int) Options.NQuadrature, DiracMaxQuadraturePoints);

  // the nucleus
  nucleus.Z = Z;
  nucleus.A = Options.A;
  nucleus.Rrms = 0.;
  nucleus.cFermi = 0.;
  nucleus.aFermi = 0.;
  if (strcicmp(Options.Nucleus, "Point") == 0) nucleus.model = NucleusPoint;
  else if (strcicmp(Options.Nucleus, "Gaussian") == 0) nucleus.model = NucleusGaussian;
  else if (strcicmp(Options.Nucleus, "Fermi") == 0) nucleus.model = NucleusFermi;
  else if (strcicmp(Options.Nucleus, "Function") == 0) nucleus.model = NucleusFunction;
  else return luaL_error(L, "Error in DiracRadialSolve: unknown Nucleus \"%s\" (Point, Gaussian, Fermi or Function)\n", Options.Nucleus);
  if ((nucleus.model == NucleusFunction) != (Options.Density != LUA_NOREF))
  {
    if (Options.Density != LUA_NOREF) luaL_unref(L, LUA_REGISTRYINDEX, Options.Density);
    return luaL_error(L, "Error in DiracRadialSolve: the option Density goes with Nucleus \"Function\", and only with it\n");
  }
  if (nucleus.model == NucleusGaussian || nucleus.model == NucleusFermi)
  {
    if (Options.Rrms > 0.) nucleus.Rrms = Options.Rrms;
    else if (Options.A > 0.) nucleus.Rrms = (0.836*pow(Options.A, 1./3.) + 0.570)/a0fm;
    else return luaL_error(L, "Error in DiracRadialSolve: a %s nucleus needs the option A or Rrms\n", Options.Nucleus);
  }
  if (nucleus.model == NucleusFermi && NucleusFermiParameters(nucleus.Rrms, Options.FermiThickness, &nucleus.cFermi, &nucleus.aFermi))
    return luaL_error(L, "Error in DiracRadialSolve: Fermi nucleus with Rrms = %g bohr is too small for the skin thickness %g bohr (Rrms^2 <= 7/5 pi^2 a^2)\n",
                      nucleus.Rrms, Options.FermiThickness);
  if (Options.rMin <= 0.) Options.rMin = (nucleus.Rrms > 0.) ? nucleus.Rrms/50. : 1e-6/fabs(Z);

  if (RadialBasisCreate(Options.Degree, Options.NIntervals, Options.rMin, Options.rMax, Options.NQuadrature, &basis))
    return luaL_error(L, "Error in DiracRadialSolve: no B-spline basis for Degree %d, NIntervals %d, rMin %g, rMax %g\n",
                      (int) Options.Degree, (int) Options.NIntervals, Options.rMin, Options.rMax);
  N = 2*(basis.nSpline - 2);

  // sample the nuclear density of a Lua function on the quadrature points, normalize to Z
  if (nucleus.model == NucleusFunction)
  {
    rho = (double *) malloc(basis.nPoint*sizeof(double));
    for (q = 0; q < basis.nPoint; q++)
    {
      lua_rawgeti(L, LUA_REGISTRYINDEX, Options.Density);
      lua_pushnumber(L, basis.rq[q]);
      lua_call(L, 1, 1);
      if (lua_type(L, -1) != LUA_TNUMBER)
      {
        free(rho);
        RadialBasisFree(&basis);
        return luaL_error(L, "Error in DiracRadialSolve: Density(%g) is not a number\n", basis.rq[q]);
      }
      rho[q] = lua_tonumber(L, -1);
      lua_pop(L, 1);
    }
    luaL_unref(L, LUA_REGISTRYINDEX, Options.Density);
    charge = NucleusCharge(&basis, rho);
    if (!(charge > 0.))
    {
      free(rho);
      RadialBasisFree(&basis);
      return luaL_error(L, "Error in DiracRadialSolve: the nuclear charge 4 pi int r^2 rho dr = %g is not positive\n", charge);
    }
    if (fabs(charge - Z) > Options.DensityTolerance*fabs(Z))
    {
      char message[256];
      snprintf(message, sizeof(message), "Warning in DiracRadialSolve: 4 pi int r^2 rho dr = %.15g differs from Z = %.15g by %.3e relative; the density is rescaled to Z\n",
               charge, Z, fabs(charge - Z)/fabs(Z));
      luaL_warning(L, "%s", message);
    }
    for (q = 0; q < basis.nPoint; q++) rho[q] *= Z/charge;
  }
  V = (double *) malloc(basis.nPoint*sizeof(double));
  if (NucleusPotential(&nucleus, &basis, rho, V))
  {
    free(V);
    free(rho);
    RadialBasisFree(&basis);
    return luaL_error(L, "Error in DiracRadialSolve: the nuclear potential failed\n");
  }

  E = (double *) malloc(N*sizeof(double));
  Edsygvd = (double *) malloc(N*sizeof(double));
  res = (double *) malloc(N*sizeof(double));
  C = (double *) malloc(N*N*sizeof(double));
  if (DiracRadialEigen(&basis, kappa, c, V, NULL, Options.RayleighQuotient ? 1 : 0, E, C, res, Edsygvd))
  {
    free(V); free(rho); free(E); free(Edsygvd); free(res); free(C);
    RadialBasisFree(&basis);
    return luaL_error(L, "Error in DiracRadialSolve: the generalized eigenvalue problem (dsygvd) failed\n");
  }
  first = DiracElectronicFirst(E, N, c);
  nStates = Options.NStates;
  if (first + nStates > N) nStates = N - first;

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
  lua_pushinteger(L, N);                     lua_setfield(L, -2, "NBasis");
  lua_pushinteger(L, basis.nKnot);           lua_setfield(L, -2, "NKnots");
  lua_pushnumber(L, Options.rMin);           lua_setfield(L, -2, "rMin");
  lua_pushnumber(L, Options.rMax);           lua_setfield(L, -2, "rMax");
  lua_createtable(L, nStates, 0);
  for (i = 0; i < nStates; i++) { lua_pushnumber(L, Edsygvd[first + i]); lua_rawseti(L, -2, i+1); }
  lua_setfield(L, -2, "Edsygvd");
  lua_createtable(L, nStates, 0);
  for (i = 0; i < nStates; i++) { lua_pushnumber(L, res[first + i]); lua_rawseti(L, -2, i+1); }
  lua_setfield(L, -2, "Residual");
  LuaPushArray_double(L, basis.rq, (int) basis.nPoint);
  lua_setfield(L, -2, "r");
  LuaPushArray_double(L, basis.wq, (int) basis.nPoint);
  lua_setfield(L, -2, "w");
  LuaPushArray_double(L, V, (int) basis.nPoint);
  lua_setfield(L, -2, "V");
  lua_newtable(L);
  lua_pushstring(L, Options.Nucleus);        lua_setfield(L, -2, "Model");
  lua_pushnumber(L, nucleus.Z);              lua_setfield(L, -2, "Z");
  lua_pushnumber(L, nucleus.A);              lua_setfield(L, -2, "A");
  lua_pushnumber(L, nucleus.Rrms);           lua_setfield(L, -2, "Rrms");
  lua_pushnumber(L, nucleus.cFermi);         lua_setfield(L, -2, "cFermi");
  lua_pushnumber(L, nucleus.aFermi);         lua_setfield(L, -2, "aFermi");
  if (nucleus.model == NucleusFunction) { lua_pushnumber(L, charge); lua_setfield(L, -2, "DensityCharge"); }
  lua_setfield(L, -2, "Nucleus");

  free(V); free(rho); free(E); free(Edsygvd); free(res); free(C);
  RadialBasisFree(&basis);
  return 4;
}
