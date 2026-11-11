/*
 * LuaAtomicOperators.cpp
 *
 * The exact paths of NewOperator("Dirac", ...) and NewOperator("AU", ...) for radial Dirac
 * spinors of type "DKBSpline" (DFT stage A5, design section 11.2), reached from
 * AtomicOperators/AtomicDirac.cpp and AtomicU.cpp:
 *
 *     H0 = NewOperator("Dirac", NF, Z, indices, G, F, {options})
 *     U  = NewOperator("AU", NF, indices, G, F)
 *
 * All functions must share one knot set, degree and number of Gauss points, and must
 * have been made with the present c = PhysicalConstants.cAU.value (the DKB functions
 * depend on c through h = 1/(2c)). G[i] and F[i] are the two components of one spinor.
 *  - Dirac: <a|h_D|b> = c_a^T H_kappa c_b with H_kappa the DKB matrix of design section 5
 *    for the nucleus of the options (Nucleus, A, Rrms, FermiThickness, Density,
 *    DensityTolerance as in DiracRadialSolve; default "Point" with charge Z), without
 *    confinement. For the solver's orbitals, nucleus and c this is the solver's matrix.
 *  - AU: R^k(13;24) = int rho_13 Y^k[rho_24]/r dr with the Galerkin solution Y^k of design
 *    section 7 on the shared knots (rho_ab = P_a P_b + Q_a Q_b), assembled exactly as
 *    Create_Relativistic_Coulomb_Operator (Martin/Slater_Rework/Coulomb_Operator.cpp) does
 *    with its grid integrals.
 */

/* 1 if the table at index arg is a list with an interpolating function of type DKBSpline
 * (then all must be, which DFTLuaDKBSpinors checks: they never reach the grid paths) */
int DFTLuaFunctionsAreDKB(lua_State *L, int arg)
{
  uint32_t i;
  int is = 0;
  if (!lua_istable(L, arg)) return 0;
  for (i = 1; i <= (uint32_t) lua_objlen(L, arg) && !is; i++)
  {
    lua_rawgeti(L, arg, i);
    Interpolation::InterpolatingFunction *f = (Interpolation::InterpolatingFunction *) luaL_testudata(L, -1, "InterpolatingFunction_Type");
    if (f != NULL && f->type == "DKBSpline") is = 1;
    lua_pop(L, 1);
  }
  return is;
}

/* The n spinors G[i], F[i] (tables at argG, argF) as coefficient vectors coef[i*Nf ..]
 * (allocated here, Nf = 2(nSpline-2)) with their kappa, and the radial basis of their knots
 * (created here). Raises a Lua error if they are no DKB spinors of one knot set, or were
 * made with another c. */
int DFTLuaDKBSpinors(lua_State *L, const char *who, int argG, int argF, uint32_t n, double c, RadialBasisType *basis,
                     int *kappa, double **coef)
{
  std::vector<Interpolation::DKBSpline *> large(n), small(n);   // not G, F: G is a Lua macro
  uint32_t i, j, ns, Nf;
  if (!lua_istable(L, argG) || !lua_istable(L, argF) || lua_objlen(L, argG) < n || lua_objlen(L, argF) < n)
    return luaL_error(L, "Error in %s: %d large and small components (G, F) expected\n", who, (int) n);
  for (i = 0; i < n; i++)
  {
    for (int comp = 0; comp < 2; comp++)
    {
      lua_rawgeti(L, comp == 0 ? argG : argF, i + 1);
      Interpolation::InterpolatingFunction *f = (Interpolation::InterpolatingFunction *) luaL_testudata(L, -1, "InterpolatingFunction_Type");
      lua_pop(L, 1);
      if (f == NULL || f->type != "DKBSpline")
        return luaL_error(L, "Error in %s: %s[%d] is not an interpolating function of type DKBSpline (all or none must be)\n", who, comp == 0 ? "G" : "F", (int) i + 1);
      if (comp == 0) large[i] = (Interpolation::DKBSpline *) f; else small[i] = (Interpolation::DKBSpline *) f;
    }
    if (large[i]->component != 0 || small[i]->component != 1)
      return luaL_error(L, "Error in %s: G[%d] must be a large and F[%d] a small component\n", who, (int) i + 1, (int) i + 1);
    if (large[i]->kappa != small[i]->kappa || large[i]->a != small[i]->a || large[i]->b != small[i]->b)
      return luaL_error(L, "Error in %s: G[%d] and F[%d] are components of different spinors\n", who, (int) i + 1, (int) i + 1);
    if (large[i]->knots != large[0]->knots || large[i]->degree != large[0]->degree || large[i]->nQ != large[0]->nQ)
      return luaL_error(L, "Error in %s: orbital %d has another knot set, degree or quadrature than orbital 1 (one knot set is required)\n", who, (int) i + 1);
    if (large[i]->speedOfLight != c)
      return luaL_error(L, "Error in %s: orbital %d was made with c = %f, PhysicalConstants.cAU.value is %f\n", who, (int) i + 1, large[i]->speedOfLight, c);
  }
  if (RadialBasisFromKnots(&(large[0]->knots[0]), (uint32_t) large[0]->knots.size(), large[0]->degree, large[0]->nQ, basis))
    return luaL_error(L, "Error in %s: no radial basis for the knots of the functions\n", who);
  ns = basis->nSpline - 2;
  Nf = 2*ns;
  *coef = (double *) malloc(n*Nf*sizeof(double));
  for (i = 0; i < n; i++)
  {
    kappa[i] = large[i]->kappa;
    for (j = 0; j < ns; j++)
    {
      (*coef)[i*Nf + j] = large[i]->a[j + 1];
      (*coef)[i*Nf + ns + j] = large[i]->b[j + 1];
    }
  }
  return 0;
}

// the indices argument of NewOperator: index lists per orbital and the field kappas
static void DFTLuaIndices(lua_State *L, const char *who, int arg, uint32_t ***indices, uint32_t *n, std::vector<int> &indToK)
{
  uint32_t i, temp;
  luaL_checktype(L, arg, LUA_TTABLE);
  *n = (uint32_t) lua_objlen(L, arg);
  lua_getfield(L, arg, "kappas");
  if (lua_isnil(L, -1)) { luaL_error(L, "Error in %s: the indices have no field kappas\n", who); return; }
  lua_check(L, -1, indToK);
  lua_pop(L, 1);
  if (indToK.size() != *n) { luaL_error(L, "Error in %s: indices.kappas and indices differ in length\n", who); return; }
  *indices = (uint32_t **) malloc((*n)*sizeof(uint32_t *));
  for (i = 0; i < *n; i++)
  {
    lua_rawgeti(L, arg, i + 1);
    (*indices)[i] = luaL_checkarray_unit32_t(L, -1, &temp);
    lua_pop(L, 1);
    if (temp != (uint32_t) (2*abs(indToK[i]))) { luaL_error(L, "Error in %s: orbital %d has %d indices, kappa = %d needs %d\n", who, (int) i + 1, (int) temp, indToK[i], 2*abs(indToK[i])); return; }
  }
}

int DFTLuaNewOperatorDiracDKB(lua_State *L)
{
  struct { DFTLuaNucleusOptions nuc; } Options;
  LuaOptionsType LuaDiracDKBOptions[] = {
    {"Nucleus", "char*", &(Options.nuc.Nucleus)},
    {"A", "double", &(Options.nuc.A)},
    {"Rrms", "double", &(Options.nuc.Rrms)},
    {"FermiThickness", "double", &(Options.nuc.FermiThickness)},
    {"Density", "function", &(Options.nuc.Density)},
    {"DensityTolerance", "double", &(Options.nuc.DensityTolerance)},
    {NULL, NULL, NULL}};
  const char *who = "NewOperator(\"Dirac\",...)";
  RadialBasisType basis;
  NucleusType nucleus;
  OperatorType *op;
  uint32_t NF, n, i, j, gr, Nf, **indices = NULL, len = (uint32_t) lua_gettop(L);
  std::vector<int> indToK;
  std::vector<SKappaGroup> groups;
  double Z, c, a0fm, *V, *coef, *H, *S, charge;
  int *kappa;
  if (len < 5 || len > 6) return luaL_error(L, "Error in %s: 5 or 6 arguments expected (NF, Z, indices, G, F, {options}), got %d\n", who, (int) len);
  NF = (uint32_t) luaL_checkunsigned(L, 1);
  Z = luaL_checknumber(L, 2);
  c = LuaGetPhysicalConstant(L, "cAU");
  a0fm = LuaGetPhysicalConstant(L, "a0")*1e15;
  DFTLuaNucleusDefaults(&Options.nuc, a0fm);
  Options.nuc.Nucleus = (char *) "Point";
  if (len == 6) LuaCheckOptions(L, 6, LuaDiracDKBOptions);
  DFTLuaIndices(L, who, 3, &indices, &n, indToK);
  kappa = (int *) malloc(n*sizeof(int));
  DFTLuaDKBSpinors(L, who, 4, 5, n, c, &basis, kappa, &coef);
  for (i = 0; i < n; i++)
    if (kappa[i] != indToK[i])
      return luaL_error(L, "Error in %s: orbital %d has kappa = %d, indices.kappas says %d\n", who, (int) i + 1, kappa[i], indToK[i]);
  DFTLuaNucleus(L, who, Z, &Options.nuc, a0fm, &nucleus);
  DFTLuaNuclearPotential(L, who, Z, &Options.nuc, &nucleus, &basis, &V, &charge);
  GroupIndicesAccordingToKappa(indToK, groups);
  Nf = 2*(basis.nSpline - 2);
  H = (double *) malloc(Nf*Nf*sizeof(double));
  S = (double *) malloc(Nf*Nf*sizeof(double));
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  InitOperator(op);
  InitOperatorFixedLength(op, 2, 8);
  op->NFermion = NF;
  op->NWord = Min(op->NFermion, 1 + (op->NFermion - 1)/8);
  for (gr = 0; gr < (uint32_t) groups.size(); gr++)
  {
    std::vector<std::vector<double> > elements(groups[gr].n, std::vector<double>(groups[gr].n, 0.));
    if (DiracDKBMatrices(&basis, groups[gr].kappa, c, V, NULL, H, S))
      return luaL_error(L, "Error in %s: the DKB matrices for kappa = %d failed\n", who, groups[gr].kappa);
    for (i = 0; i < groups[gr].n; i++)
      for (j = i; j < groups[gr].n; j++)
        elements[i][j] = AtomBilinear(coef + groups[gr].indices[i]*Nf, H, coef + groups[gr].indices[j]*Nf, Nf);
    OperatorAddKappaBlock(op, &groups[gr], indices, elements);
  }
  luaL_setmetatable(L, "Operator_Type");
  for (gr = 0; gr < (uint32_t) groups.size(); gr++) free(groups[gr].indices);
  for (i = 0; i < n; i++) free(indices[i]);
  free(indices); free(kappa); free(coef); free(V); free(H); free(S);
  RadialBasisFree(&basis);
  return 1;
}

int DFTLuaNewOperatorUDKB(lua_State *L)
{
  const char *who = "NewOperator(\"AU\",...)";
  RadialBasisType basis;
  uint32_t NF, n, i, i1, i2, i3, i4, q, k, nP, nPois, kMax = 0, len = (uint32_t) lua_gettop(L);
  std::vector<std::vector<uint16_t> > indices;
  std::vector<int> indToK;
  double c, *coef, *P, *Q, *dP, *dQ, **band, *rho, **Ycache;
  int *kappa;
  if (len != 4) return luaL_error(L, "Error in %s: for radial functions of type DKBSpline 4 arguments (NF, indices, G, F), got %d (no options)\n", who, (int) len);
  NF = (uint32_t) luaL_checkunsigned(L, 1);
  c = LuaGetPhysicalConstant(L, "cAU");
  lua_check(L, 2, indices);
  lua_getfield(L, 2, "kappas");
  if (lua_isnil(L, -1)) return luaL_error(L, "Error in %s: the indices have no field kappas\n", who);
  lua_check(L, -1, indToK);
  lua_pop(L, 1);
  n = (uint32_t) indices.size();
  if (indToK.size() != n) return luaL_error(L, "Error in %s: indices.kappas and indices differ in length\n", who);
  if (lua_objlen(L, 3) != n || lua_objlen(L, 4) != n) return luaL_error(L, "Error in %s: G and F need one function per orbital of the indices (%d)\n", who, (int) n);
  kappa = (int *) malloc(n*sizeof(int));
  DFTLuaDKBSpinors(L, who, 3, 4, n, c, &basis, kappa, &coef);
  for (i = 0; i < n; i++)
  {
    if (kappa[i] != indToK[i])
      return luaL_error(L, "Error in %s: orbital %d has kappa = %d, indices.kappas says %d\n", who, (int) i + 1, kappa[i], indToK[i]);
    if ((uint32_t) (2*abs(kappa[i])) > kMax) kMax = (uint32_t) (2*abs(kappa[i]));
  }
  // P, Q of every orbital on the quadrature points; Poisson factorizations k = 0 .. kMax-1
  nP = basis.nPoint;
  nPois = basis.nSpline - 1;
  P = (double *) malloc(n*nP*sizeof(double));
  Q = (double *) malloc(n*nP*sizeof(double));
  dP = (double *) malloc(nP*sizeof(double));
  dQ = (double *) malloc(nP*sizeof(double));
  rho = (double *) malloc(nP*sizeof(double));
  for (i = 0; i < n; i++) AtomOrbitalOnPoints(&basis, 1, kappa[i], c, coef + i*2*(basis.nSpline - 2), P + i*nP, dP, Q + i*nP, dQ);
  band = (double **) calloc(kMax, sizeof(double *));
  for (k = 0; k < kMax; k++)
  {
    band[k] = (double *) malloc((basis.degree + 1)*nPois*sizeof(double));
    if (BSplinePoissonSetup(basis.knot, basis.nKnot, basis.degree, basis.rq, basis.wq, nP, k, band[k], basis.degree + 1))
      return luaL_error(L, "Error in %s: the Poisson matrix for k = %d failed\n", who, (int) k);
  }
  // Y^k[rho_24]/r on the points, computed once per pair (i2 <= i4) and k
  Ycache = (double **) calloc(n*n*kMax, sizeof(double *));
  OperatorType *U = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  InitNamedOperator(U, "Coulomb Operator");
  U->NFermion = NF;
  U->NWord = Min(U->NFermion, 1 + (U->NFermion - 1)/8);
  InitOperatorFixedLength(U, 2, 8);
  InitOperatorFixedLength(U, 4, 8);
  // the tuples of Create_Relativistic_Coulomb_Operator: i2 = i1, and i2 > i1
  for (i1 = 0; i1 < n; i1++)
    for (i2 = i1; i2 < n; i2++)
      for (i3 = i1; i3 < n; i3++)
        for (i4 = (i2 == i1) ? i3 : i2; i4 < n; i4++)
        {
          int k1 = kappa[i1], k2 = kappa[i2], k3 = kappa[i3], k4 = kappa[i4];
          int L0 = Max(abs(abs(k1) - abs(k3)), abs(abs(k2) - abs(k4)));
          int L1 = Min(abs(k1) + abs(k3), abs(k2) + abs(k4)) - 1;
          if ((abs(k1) + abs(k3) + (Sign(k1) + Sign(k3))/2 + L0) % 2 == 0) L0++;
          if (L0 > L1) continue;
          double Fk[(L1 - L0)/2 + 1];
          uint32_t a = Min(i2, i4), b = Max(i2, i4);
          for (int Lk = L0; Lk <= L1; Lk += 2)
          {
            double **Y = Ycache + (a*n + b)*kMax + Lk;
            if (*Y == NULL)
            {
              *Y = (double *) malloc(nP*sizeof(double));
              for (q = 0; q < nP; q++) rho[q] = P[a*nP + q]*P[b*nP + q] + Q[a*nP + q]*Q[b*nP + q];
              if (RadialYkOnPoints(&basis, (uint32_t) Lk, band[Lk], rho, *Y))
                return luaL_error(L, "Error in %s: the Poisson equation for k = %d failed\n", who, Lk);
            }
            for (q = 0; q < nP; q++) rho[q] = P[i1*nP + q]*P[i3*nP + q] + Q[i1*nP + q]*Q[i3*nP + q];
            Fk[(Lk - L0)/2] = HFIntegral(&basis, rho, *Y);
          }
          OperatorType optemp;
          CreateOperatorAtomicURelativisticOneIntegral(NF, k1, k2, k3, k4, &(indices[i1][0]), &(indices[i2][0]), &(indices[i3][0]), &(indices[i4][0]),
                                                       Fk, i1, i2, i3, i4, &optemp);
          OperatorAdd(U, optemp, 1., 0.);
          FreeOperator(&optemp);
        }
  luaL_setmetatable(L, "Operator_Type");
  for (i = 0; i < n*n*kMax; i++) free(Ycache[i]);
  for (k = 0; k < kMax; k++) free(band[k]);
  free(Ycache); free(band); free(P); free(Q); free(dP); free(dQ); free(rho); free(coef); free(kappa);
  RadialBasisFree(&basis);
  return 1;
}
