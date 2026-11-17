/*
 * LuaAtom.cpp
 *
 * Lua helpers of the DFT atom (used by CalculateAtomicBasisFunctions, LuaAtomicBasis.cpp)
 * and the entry point of the exchange-correlation functionals:
 *
 *     e, dedn, dedsigma = XCFunctional(method, n, sigma)
 *
 * The temporary entry points DiracKohnSham and DiracHartreeFock of stages A3 and A4 are
 * merged into CalculateAtomicBasisFunctions ({"Method", ...} and the fields of Atom).
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

/* Push G and F, the radial functions of the atom's orbitals in the sign convention of
 * AtomSignConvention (applied here): Dirac, G = P and F = Q as DKBSplines; Schroedinger,
 * G = P as a BSpline and F empty; both give the quadrature points as f.x. list[0 .. nList-1]
 * picks the orbitals and their order (NULL: all, in the atom's order). */
void DFTLuaPushRadialFunctions(lua_State *L, DiracAtomType *atom, const uint32_t *list, uint32_t nList)
{
  RadialBasisType *basis = atom->basis;
  uint32_t k, o, i;
  if (list == NULL) nList = atom->nOrb;
  lua_createtable(L, nList, 0);
  lua_createtable(L, atom->relativistic ? nList : 0, 0);
  for (k = 0; k < nList; k++)
  {
    o = (list == NULL) ? k : list[k];
    AtomSignConvention(atom, o);
    if (atom->relativistic)
    {
      for (int comp = 0; comp < 2; comp++)
      {
        void *memory = lua_newuserdata(L, sizeof(Interpolation::DKBSpline));
        Interpolation::DKBSpline *f = new(memory) Interpolation::DKBSpline;
        luaL_setmetatable(L, "InterpolatingFunction_Type");
        f->Create(basis->knot, basis->nKnot, basis->degree, atom->orb[o].kappa, (uint32_t) comp, atom->c, atom->orb[o].coef, basis->nQ);
        lua_rawseti(L, comp == 0 ? -3 : -2, k+1);
      }
    }
    else
    {
      // P = sum_i coef_i B_{i+1}, a B-spline on the knots (the first and the last coefficient 0)
      std::vector<double> x, knots(basis->knot, basis->knot + basis->nKnot), y(basis->nSpline, 0.);
      for (i = 0; i < basis->nKnot; i++) if (i == 0 || basis->knot[i] > basis->knot[i-1]) x.push_back(basis->knot[i]);
      for (i = 0; i + 2 < basis->nSpline; i++) y[i+1] = atom->orb[o].coef[i];
      void *memory = lua_newuserdata(L, sizeof(BSpline::TBSpline));
      BSpline::TBSpline *f = new(memory) BSpline::TBSpline(x, knots, y, basis->degree);
      f->points.assign(basis->rq, basis->rq + basis->nPoint);   // f.x: the quadrature points, as for a DKBSpline
      luaL_setmetatable(L, "InterpolatingFunction_Type");
      lua_rawseti(L, -3, k+1);
    }
  }
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
