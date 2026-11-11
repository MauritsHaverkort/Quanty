/*
 * LuaAtom.cpp
 *
 * Temporary Lua entry points of DFT stages A3 and A4, until CalculateAtomicBasisFunctions (A5):
 *
 *     Atom, G, F = DiracKohnSham(Z, configuration, {options})
 *     Atom, G, F = DiracHartreeFock(Z, configuration, {options})
 *     e, dedn, dedsigma = XCFunctional(method, n, sigma)
 */

/*
 * The self-consistent Kohn-Sham atom (AtomSCF.cpp): Dirac (default) or Schroedinger,
 * spherical and spin unpolarized, with the mean occupations of the configuration, e.g.
 * "[Ne] 3s2 3p6" or "1s2 2s2 2p6" (fractional occupations allowed; Dirac: each nl shell is
 * shared between j = l-1/2 and l+1/2 in proportion to 2j+1).
 * Returns the table Atom with the energies (Etot = Ekin + Ecoul + Eenuc + Exc, the
 * components as in the NIST atomic reference data, and EtotEigenvalueSum, the expression
 * through the eigenvalues, equal to Etot at self-consistency), Converged, Iterations,
 * Residual (of the mixed quantity), PotentialResidual (sqrt(int f |v_out - v_in|^2 + ...),
 * which bounds |EtotEigenvalueSum - Etot| by about sqrt(NElectron) times it), dE, and Orbitals[i] = {Label, n, l, kappa, j, Occupation, Energy}; and the
 * radial functions: Dirac, G = P and F = Q as DKBSplines; Schroedinger, G = P as a BSpline
 * and F empty. Options: those of the nucleus and the basis as in DiracRadialSolve
 * (Nucleus, A, Rrms, FermiThickness, Density, DensityTolerance, Degree, NIntervals, rMin,
 * rMax with default 40, NQuadrature), and
 *     Method            "LDA" (Slater + PW92, default), "LDA-VWN" (Slater + VWN5, NIST's
 *                       LDA), "RLDA-MV" (VWN5 and the MacDonald-Vosko factor on exchange,
 *                       NIST's RLDA), "PBE", "PBE-VWN" (PBE on VWN5)
 *     Relativistic      true (default): Dirac; false: Schroedinger (the c -> infinity
 *                       limit, which c = PhysicalConstants.cAU.value large can not give in
 *                       double precision)
 *     MaxIterations     200
 *     History           Pulay history, 8
 *     Mixing            step of the Pulay mixing, 0.3
 *     EnergyTolerance   change of Etot between iterations, 1e-10 Ha
 *     ResidualTolerance norm of the potential residual, 1e-8
 *     DensityCutoff     density below which e_xc = 0, 1e-30
 *     MixDensity        Pulay mixing of the density (f and f') instead of the potential
 *                       (v and g); default true for the GGAs (PBE, PBE-VWN), whose Dirac
 *                       atoms do not converge with potential mixing, false for the LDAs
 *                       (decided by Maurits, 2026-11-11)
 *     RayleighQuotient, OriginBoundary   as in DiracRadialSolve, true
 */

/* 1 if the options table at index arg has an entry {name, value} (case insensitive). */
int DFTLuaHasOption(lua_State *L, int arg, const char *name)
{
  uint32_t i, len;
  int found = 0;
  if (!lua_istable(L, arg)) return 0;
  len = (uint32_t) lua_objlen(L, arg);
  for (i = 1; i <= len && !found; i++)
  {
    lua_rawgeti(L, arg, i);
    if (lua_istable(L, -1))
    {
      lua_rawgeti(L, -1, 1);
      if (lua_type(L, -1) == LUA_TSTRING && strcicmp(lua_tostring(L, -1), name) == 0) found = 1;
      lua_pop(L, 1);
    }
    lua_pop(L, 1);
  }
  return found;
}
/* Push the table Orbitals[i] = {Label, n, l, kappa, j, Occupation, Energy} of the atom
 * (kappa and j only for Dirac). */
void DFTLuaPushOrbitals(lua_State *L, DiracAtomType *atom)
{
  const char *letters = "spdfghi";
  char label[32];
  uint32_t o;
  lua_createtable(L, atom->nOrb, 0);
  for (o = 0; o < atom->nOrb; o++)
  {
    DiracOrbitalType *orb = atom->orb + o;
    lua_newtable(L);
    if (atom->relativistic) snprintf(label, sizeof(label), "%d%c%d/2", orb->n, letters[orb->l], 2*abs(orb->kappa) - 1);
    else snprintf(label, sizeof(label), "%d%c", orb->n, letters[orb->l]);
    lua_pushstring(L, label);       lua_setfield(L, -2, "Label");
    lua_pushinteger(L, orb->n);     lua_setfield(L, -2, "n");
    lua_pushinteger(L, orb->l);     lua_setfield(L, -2, "l");
    if (atom->relativistic)
    {
      lua_pushinteger(L, orb->kappa); lua_setfield(L, -2, "kappa");
      lua_pushnumber(L, abs(orb->kappa) - 0.5); lua_setfield(L, -2, "j");
    }
    lua_pushnumber(L, orb->q);      lua_setfield(L, -2, "Occupation");
    lua_pushnumber(L, orb->energy); lua_setfield(L, -2, "Energy");
    lua_rawseti(L, -2, o+1);
  }
}

/* Push G and F, the radial functions of the atom's orbitals in the sign convention of
 * AtomSignConvention (applied here): Dirac, G = P and F = Q as DKBSplines; Schroedinger,
 * G = P as a BSpline and F empty. */
void DFTLuaPushRadialFunctions(lua_State *L, DiracAtomType *atom)
{
  RadialBasisType *basis = atom->basis;
  uint32_t o, i;
  lua_createtable(L, atom->nOrb, 0);
  lua_createtable(L, atom->relativistic ? atom->nOrb : 0, 0);
  for (o = 0; o < atom->nOrb; o++)
  {
    AtomSignConvention(atom, o);
    if (atom->relativistic)
    {
      for (int comp = 0; comp < 2; comp++)
      {
        void *memory = lua_newuserdata(L, sizeof(Interpolation::DKBSpline));
        Interpolation::DKBSpline *f = new(memory) Interpolation::DKBSpline;
        luaL_setmetatable(L, "InterpolatingFunction_Type");
        f->Create(basis->knot, basis->nKnot, basis->degree, atom->orb[o].kappa, (uint32_t) comp, atom->c, atom->orb[o].coef, basis->nQ);
        lua_rawseti(L, comp == 0 ? -3 : -2, o+1);
      }
    }
    else
    {
      // P = sum_i coef_i B_{i+1}, a B-spline on the knots (the first and the last coefficient 0)
      std::vector<double> x, knots(basis->knot, basis->knot + basis->nKnot), y(basis->nSpline, 0.);
      for (i = 0; i < basis->nKnot; i++) if (i == 0 || basis->knot[i] > basis->knot[i-1]) x.push_back(basis->knot[i]);
      for (i = 0; i + 2 < basis->nSpline; i++) y[i+1] = atom->orb[o].coef[i];
      void *memory = lua_newuserdata(L, sizeof(BSpline::TBSpline));
      new(memory) BSpline::TBSpline(x, knots, y, basis->degree);
      luaL_setmetatable(L, "InterpolatingFunction_Type");
      lua_rawseti(L, -3, o+1);
    }
  }
}

static int LuaDiracKohnSham(lua_State *L)
{
  struct
  {
    DFTLuaNucleusOptions nuc;
    char *Method;
    bool Relativistic, RayleighQuotient, OriginBoundary, MixDensity;
    uint32_t MaxIterations, History;
    double Mixing, EnergyTolerance, ResidualTolerance, DensityCutoff;
  } Options;
  LuaOptionsType LuaDiracKohnShamOptions[] = {
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
    {"Method", "char*", &(Options.Method)},
    {"Relativistic", "bool", &(Options.Relativistic)},
    {"RayleighQuotient", "bool", &(Options.RayleighQuotient)},
    {"OriginBoundary", "bool", &(Options.OriginBoundary)},
    {"MaxIterations", "uint32_t", &(Options.MaxIterations)},
    {"History", "uint32_t", &(Options.History)},
    {"Mixing", "double", &(Options.Mixing)},
    {"EnergyTolerance", "double", &(Options.EnergyTolerance)},
    {"ResidualTolerance", "double", &(Options.ResidualTolerance)},
    {"DensityCutoff", "double", &(Options.DensityCutoff)},
    {"MixDensity", "bool", &(Options.MixDensity)},
    {NULL, NULL, NULL}};
  const uint32_t maxShell = 64;
  RadialBasisType basis;
  NucleusType nucleus;
  DiracAtomType atom;
  double c, a0fm, Z, *V, charge, shellQ[64];
  int shellN[64], shellL[64], method, error;
  uint32_t len = (uint32_t) lua_gettop(L), nShell = 0;
  const char *configuration;

  if (len < 2 || len > 3) return luaL_error(L, "Error in DiracKohnSham: 2 or 3 arguments expected (Z, configuration, {options}), got %d\n", (int) len);
  Z = luaL_checknumber(L, 1);
  configuration = luaL_checkstring(L, 2);
  c = LuaGetPhysicalConstant(L, "cAU");
  a0fm = LuaGetPhysicalConstant(L, "a0")*1e15;
  DFTLuaNucleusDefaults(&Options.nuc, a0fm);
  Options.Method = (char *) "LDA";
  Options.Relativistic = true;
  Options.RayleighQuotient = true;
  Options.OriginBoundary = true;
  Options.MaxIterations = 200;
  Options.History = 8;
  Options.Mixing = 0.3;
  Options.EnergyTolerance = 1e-10;
  Options.ResidualTolerance = 1e-8;
  Options.DensityCutoff = 1e-30;
  Options.MixDensity = false;
  if (len == 3) LuaCheckOptions(L, 3, LuaDiracKohnShamOptions);
  method = XCMethodFromName(Options.Method);
  if (method < 0) return luaL_error(L, "Error in DiracKohnSham: unknown Method \"%s\" (LDA, LDA-VWN, RLDA-MV, PBE, PBE-VWN)\n", Options.Method);
  if (len < 3 || !DFTLuaHasOption(L, 3, "MixDensity")) Options.MixDensity = XCIsGGA(method) ? true : false;
  if (AtomParseConfiguration(configuration, maxShell, &nShell, shellN, shellL, shellQ))
    return luaL_error(L, "Error in DiracKohnSham: can not read the configuration \"%s\"\n", configuration);
  if (Options.nuc.rMax <= 0.) Options.nuc.rMax = 40.;
  DFTLuaNucleusBasis(L, "DiracKohnSham", Z, &Options.nuc, a0fm, &nucleus, &basis, &V, &charge);

  error = AtomSetup(&atom, &basis, V, Z, Options.Relativistic ? 1 : 0, c, method, Options.OriginBoundary ? 1 : 0, Options.DensityCutoff,
                    nShell, shellN, shellL, shellQ);
  if (!error) error = AtomSCF(&atom, Options.MaxIterations, Options.History, Options.Mixing, Options.EnergyTolerance,
                              Options.ResidualTolerance, Options.RayleighQuotient ? 1 : 0, Options.MixDensity ? 1 : 0);
  free(V);
  if (error)
  {
    AtomFree(&atom);
    RadialBasisFree(&basis);
    return luaL_error(L, "Error in DiracKohnSham: the self-consistent field failed (an eigenproblem or the Poisson equation)\n");
  }
  if (!atom.converged)
    luaL_warning(L, "Warning in DiracKohnSham: not converged in %d iterations (dE = %f, residual = %f)\n",
                 (int) atom.iterations, atom.dE, atom.residual);

  // Atom
  lua_newtable(L);
  lua_pushnumber(L, atom.Etot);   lua_setfield(L, -2, "Etot");
  lua_pushnumber(L, atom.Ekin);   lua_setfield(L, -2, "Ekin");
  lua_pushnumber(L, atom.Ecoul);  lua_setfield(L, -2, "Ecoul");
  lua_pushnumber(L, atom.Eenuc);  lua_setfield(L, -2, "Eenuc");
  lua_pushnumber(L, atom.Exc);    lua_setfield(L, -2, "Exc");
  lua_pushnumber(L, atom.Eeig);   lua_setfield(L, -2, "EtotEigenvalueSum");
  lua_pushboolean(L, atom.converged); lua_setfield(L, -2, "Converged");
  lua_pushinteger(L, atom.iterations); lua_setfield(L, -2, "Iterations");
  lua_pushnumber(L, atom.residual); lua_setfield(L, -2, "Residual");
  lua_pushnumber(L, atom.potentialResidual); lua_setfield(L, -2, "PotentialResidual");
  lua_pushnumber(L, atom.dE);     lua_setfield(L, -2, "dE");
  lua_pushnumber(L, atom.nElectron); lua_setfield(L, -2, "NElectron");
  lua_pushnumber(L, c);           lua_setfield(L, -2, "c");
  lua_pushstring(L, Options.Method); lua_setfield(L, -2, "Method");
  lua_pushboolean(L, Options.Relativistic); lua_setfield(L, -2, "Relativistic");
  lua_pushboolean(L, Options.MixDensity); lua_setfield(L, -2, "MixDensity");
  lua_pushnumber(L, Options.nuc.rMin); lua_setfield(L, -2, "rMin");
  lua_pushnumber(L, Options.nuc.rMax); lua_setfield(L, -2, "rMax");
  DFTLuaPushOrbitals(L, &atom);
  lua_setfield(L, -2, "Orbitals");
  DFTLuaPushRadialFunctions(L, &atom);
  AtomFree(&atom);
  RadialBasisFree(&basis);
  return 3;
}

/*
 * The configuration-average Hartree-Fock atom (AtomHF.cpp, design section 8): Dirac
 * (default; subshells n kappa) or Schroedinger (Cowan's average of configuration), for
 * the configuration as in DiracKohnSham, every shell with q > 0. Starts from the
 * self-consistent LDA-VWN atom of the configuration.
 * Returns the table Atom with Etot = Ekin + Eenuc + Edirect + Eexchange (Ekin = sum qbar
 * <a|H0|a>, H0 without the nucleus; Edirect the F^0 terms; Eexchange the rest),
 * EtotEigenvalueSum = 1/2 sum qbar (eps + I) (equal to Etot for any orbitals, a check of
 * the Fock operators against the energy), Converged, Iterations, Residual (the orbital
 * step max_a |c_out - c_in|_S of the last iteration), Gradient (max_a |dE/dc_a| towards
 * the unoccupied orbitals, dual norm; it has a rounding floor of 1e-8 .. 1e-6, so it is
 * reported, not used), PairGradient (max |dE/dt| over rotations of two occupied subshells
 * of one kappa, Brillouin), dE, NElectron, c, Relativistic, rMin, rMax, Orbitals as in
 * DiracKohnSham (Energy = eps_a = <a|F_a|a>); Schroedinger only: EavJJ, E_av of the
 * converged radial functions with the relativistic angular factors and pair counts of the
 * j = l -+ 1/2 subshells, equal to Etot (design check 5). G and F as in DiracKohnSham.
 * Options: those of the nucleus and the basis as in DiracKohnSham, and
 *     Relativistic      true (default): Dirac; false: Schroedinger
 *     MaxIterations     200
 *     History           Pulay history, 8
 *     Mixing            step of the Pulay mixing, 0.5
 *     EnergyTolerance   change of Etot between iterations, 1e-10 Ha
 *     ResidualTolerance the orbital step, 1e-8 (the eigensolver noise of a point-nucleus
 *                       Dirac basis puts its floor at 2e-9 for Hg; E then varies by 4e-11)
 *     RayleighQuotient, OriginBoundary   as in DiracRadialSolve, true
 */
static int LuaDiracHartreeFock(lua_State *L)
{
  struct
  {
    DFTLuaNucleusOptions nuc;
    bool Relativistic, RayleighQuotient, OriginBoundary;
    uint32_t MaxIterations, History;
    double Mixing, EnergyTolerance, ResidualTolerance;
  } Options;
  LuaOptionsType LuaDiracHartreeFockOptions[] = {
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
    {"Relativistic", "bool", &(Options.Relativistic)},
    {"RayleighQuotient", "bool", &(Options.RayleighQuotient)},
    {"OriginBoundary", "bool", &(Options.OriginBoundary)},
    {"MaxIterations", "uint32_t", &(Options.MaxIterations)},
    {"History", "uint32_t", &(Options.History)},
    {"Mixing", "double", &(Options.Mixing)},
    {"EnergyTolerance", "double", &(Options.EnergyTolerance)},
    {"ResidualTolerance", "double", &(Options.ResidualTolerance)},
    {NULL, NULL, NULL}};
  const uint32_t maxShell = 64;
  RadialBasisType basis;
  NucleusType nucleus;
  DiracAtomType atom;
  DiracHFType hf, jj;
  double c, a0fm, Z, *V, charge, shellQ[64], EavJJ = 0., e1, e2, e3, e4;
  int shellN[64], shellL[64], error;
  uint32_t len = (uint32_t) lua_gettop(L), nShell = 0, s, a;
  const char *configuration;

  if (len < 2 || len > 3) return luaL_error(L, "Error in DiracHartreeFock: 2 or 3 arguments expected (Z, configuration, {options}), got %d\n", (int) len);
  Z = luaL_checknumber(L, 1);
  configuration = luaL_checkstring(L, 2);
  c = LuaGetPhysicalConstant(L, "cAU");
  a0fm = LuaGetPhysicalConstant(L, "a0")*1e15;
  DFTLuaNucleusDefaults(&Options.nuc, a0fm);
  Options.Relativistic = true;
  Options.RayleighQuotient = true;
  Options.OriginBoundary = true;
  Options.MaxIterations = 200;
  Options.History = 8;
  Options.Mixing = 0.5;
  Options.EnergyTolerance = 1e-10;
  Options.ResidualTolerance = 1e-8;
  if (len == 3) LuaCheckOptions(L, 3, LuaDiracHartreeFockOptions);
  if (AtomParseConfiguration(configuration, maxShell, &nShell, shellN, shellL, shellQ))
    return luaL_error(L, "Error in DiracHartreeFock: can not read the configuration \"%s\"\n", configuration);
  for (s = 0; s < nShell; s++)
    if (!(shellQ[s] > 0.))
      return luaL_error(L, "Error in DiracHartreeFock: shell %d%c of \"%s\" is empty (every shell needs q > 0)\n", shellN[s], "spdfghi"[shellL[s]], configuration);
  if (Options.nuc.rMax <= 0.) Options.nuc.rMax = 40.;
  DFTLuaNucleusBasis(L, "DiracHartreeFock", Z, &Options.nuc, a0fm, &nucleus, &basis, &V, &charge);

  // the LDA atom as the start, then Hartree-Fock
  error = AtomSetup(&atom, &basis, V, Z, Options.Relativistic ? 1 : 0, c, XCLDAVWN, Options.OriginBoundary ? 1 : 0, 1e-30,
                    nShell, shellN, shellL, shellQ);
  if (!error) error = AtomSCF(&atom, Options.MaxIterations, 8, 0.3, 1e-8, 1e-6, Options.RayleighQuotient ? 1 : 0, 0);
  free(V);
  if (!error) error = HFSetup(&hf, &atom, 0, nShell, shellN, shellL, shellQ);
  if (!error)
  {
    error = HFSCF(&hf, Options.MaxIterations, Options.History, Options.Mixing, Options.EnergyTolerance, Options.ResidualTolerance,
                  Options.RayleighQuotient ? 1 : 0);
    if (!error && !Options.Relativistic)
    {
      error = HFSetup(&jj, &atom, 1, nShell, shellN, shellL, shellQ);
      if (!error)
      {
        HFOrbitalsOnPoints(&jj);
        error = HFEnergy(&jj, &EavJJ, &e1, &e2, &e3, &e4);
        HFFree(&jj);
      }
    }
    if (error) HFFree(&hf);
  }
  if (error)
  {
    AtomFree(&atom);
    RadialBasisFree(&basis);
    return luaL_error(L, "Error in DiracHartreeFock: the self-consistent field failed (an eigenproblem, the Poisson equation or linearly dependent orbitals)\n");
  }
  if (!hf.converged)
    luaL_warning(L, "Warning in DiracHartreeFock: not converged in %d iterations (dE = %f, residual = %f)\n",
                 (int) hf.iterations, hf.dE, hf.residual);
  for (a = 0; a < hf.nSub; a++) atom.orb[hf.sub[a].orbital].energy = hf.sub[a].energy;

  // Atom
  lua_newtable(L);
  lua_pushnumber(L, hf.Etot);      lua_setfield(L, -2, "Etot");
  lua_pushnumber(L, hf.Ekin);      lua_setfield(L, -2, "Ekin");
  lua_pushnumber(L, hf.Eenuc);     lua_setfield(L, -2, "Eenuc");
  lua_pushnumber(L, hf.Edirect);   lua_setfield(L, -2, "Edirect");
  lua_pushnumber(L, hf.Eexchange); lua_setfield(L, -2, "Eexchange");
  lua_pushnumber(L, hf.Eidentity); lua_setfield(L, -2, "EtotEigenvalueSum");
  if (!Options.Relativistic) { lua_pushnumber(L, EavJJ); lua_setfield(L, -2, "EavJJ"); }
  lua_pushboolean(L, hf.converged); lua_setfield(L, -2, "Converged");
  lua_pushinteger(L, hf.iterations); lua_setfield(L, -2, "Iterations");
  lua_pushnumber(L, hf.residual);  lua_setfield(L, -2, "Residual");
  lua_pushnumber(L, hf.gradient);  lua_setfield(L, -2, "Gradient");
  lua_pushnumber(L, hf.pairGradient); lua_setfield(L, -2, "PairGradient");
  lua_pushnumber(L, hf.dE);        lua_setfield(L, -2, "dE");
  lua_pushnumber(L, atom.nElectron); lua_setfield(L, -2, "NElectron");
  lua_pushnumber(L, c);            lua_setfield(L, -2, "c");
  lua_pushboolean(L, Options.Relativistic); lua_setfield(L, -2, "Relativistic");
  lua_pushnumber(L, Options.nuc.rMin); lua_setfield(L, -2, "rMin");
  lua_pushnumber(L, Options.nuc.rMax); lua_setfield(L, -2, "rMax");
  DFTLuaPushOrbitals(L, &atom);
  lua_setfield(L, -2, "Orbitals");
  DFTLuaPushRadialFunctions(L, &atom);
  HFFree(&hf);
  AtomFree(&atom);
  RadialBasisFree(&basis);
  return 3;
}

/*
 * e, dedn, dedsigma = XCFunctional(method, n, sigma)
 * The exchange-correlation energy per volume e(n, sigma) of a spin-unpolarized density n
 * with sigma = |grad n|^2 (default 0), and its derivatives, for method "LDA", "LDA-VWN",
 * "RLDA-MV" (with c = PhysicalConstants.cAU.value), "PBE" or "PBE-VWN" (XCFunctional.cpp).
 */
static int LuaXCFunctional(lua_State *L)
{
  uint32_t len = (uint32_t) lua_gettop(L);
  const char *name;
  int method;
  double n, sigma = 0., e, en, es;
  if (len < 2 || len > 3) return luaL_error(L, "Error in XCFunctional: 2 or 3 arguments expected (method, n, sigma), got %d\n", (int) len);
  name = luaL_checkstring(L, 1);
  method = XCMethodFromName(name);
  if (method < 0) return luaL_error(L, "Error in XCFunctional: unknown method \"%s\"\n", name);
  n = luaL_checknumber(L, 2);
  if (len == 3) sigma = luaL_checknumber(L, 3);
  XCEvaluate(method, LuaGetPhysicalConstant(L, "cAU"), n, sigma, 0., &e, &en, &es);
  lua_pushnumber(L, e);
  lua_pushnumber(L, en);
  lua_pushnumber(L, es);
  return 3;
}
