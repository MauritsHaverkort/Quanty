/*
 * LuaDiracRadial.cpp
 *
 * The options of the DFT Lua entries (CalculateAtomicBasisFunctions, NewOperator("Dirac"),
 * NewOperator("NuclearMultipole")) that fix the nucleus and the radial basis, and their
 * setup (DFTLuaNucleusDefaults, DFTLuaNucleus, DFTLuaNuclearPotential, DFTLuaNucleusBasis).
 * Options (lengths in bohr):
 *     Nucleus          "Point", "Gaussian" (default), "Uniform" (sphere, R = sqrt(5/3) Rrms),
 *                      "Fermi", "Function" or "Function3D"
 *     A, Rrms          mass number, or the root-mean-square radius; Rrms = (0.836 A^(1/3)
 *                      + 0.570) fm if only A is given (Visscher and Dyall)
 *     FermiThickness   skin thickness t of the Fermi nucleus, default 2.30 fm
 *     Density          Lua function rho(r) (charge per volume) for Nucleus "Function", or
 *                      rho(x, y, z) for "Function3D"; sampled on the quadrature points (no
 *                      call back from C)
 *     DensityTolerance |4 pi int r^2 rho - Z| / Z above which a warning is printed,
 *                      default 1e-10; the density is always rescaled to Z
 *     NucleusLmax, NucleusNTheta, NucleusNPhi   Function3D: the even L <= Lmax of the
 *                      multipole expansion and its angular grid
 *     Degree           B-spline degree, default 8 (order 9)
 *     NIntervals       knot intervals, default 120
 *     rMin, rMax       first interior knot and the end of the knots; default rMin =
 *                      Rrms/50 for a finite nucleus and 1e-6/Z for a point nucleus,
 *                      rMax = 40 (or from the confinement of CalculateAtomicBasisFunctions)
 *     NQuadrature      Gauss points per knot interval, default Degree+5, at most 127
 * c is PhysicalConstants.cAU.value, 1 bohr = PhysicalConstants.a0.value m.
 * The temporary test entries of stages A1 and A2 (DiracRadialSolve, RadialPoisson,
 * DiracOrthonormalize) are merged into CalculateAtomicBasisFunctions: a configuration
 * without electrons is the bare nucleus, Atom.Potential.rVH the Poisson solution, and the
 * option Orthonormalize the on-site orthonormalization.
 */

/* The options of the DFT Lua entries that fix the nucleus and the radial basis. */
typedef struct
{
  char *Nucleus;
  double A, Rrms, FermiThickness, DensityTolerance, rMin, rMax;
  uint32_t Degree, NIntervals, NQuadrature;
  int Density;
  uint32_t NucleusLmax, NucleusNTheta, NucleusNPhi;   // Function3D: even L <= Lmax; the angular grid (0: Lmax+8, 2 Lmax+16)
  double *multipoleMoment;   // Function3D, set by the caller to receive q_LM (NucleusMultipoleCount(Lmax)), or NULL
  double *multipoleV;        // Function3D, set by the caller to receive V_LM on the basis points, or NULL
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
  o->NucleusLmax = 4;
  o->NucleusNTheta = 0;
  o->NucleusNPhi = 0;
  o->multipoleMoment = NULL;
  o->multipoleV = NULL;
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
  else if (strcicmp(o->Nucleus, "Function3D") == 0) nucleus->model = NucleusFunction3D;
  else if (strcicmp(o->Nucleus, "Uniform") == 0) nucleus->model = NucleusUniform;
  else return luaL_error(L, "Error in %s: unknown Nucleus \"%s\" (Point, Gaussian, Uniform, Fermi, Function or Function3D)\n", who, o->Nucleus);
  if ((nucleus->model == NucleusFunction || nucleus->model == NucleusFunction3D) != (o->Density != LUA_NOREF))
  {
    if (o->Density != LUA_NOREF) luaL_unref(L, LUA_REGISTRYINDEX, o->Density);
    return luaL_error(L, "Error in %s: the option Density goes with Nucleus \"Function\" or \"Function3D\", and only with them\n", who);
  }
  if (nucleus->model == NucleusFunction3D)
  {
    if (o->NucleusNTheta == 0) o->NucleusNTheta = o->NucleusLmax + 8;
    if (o->NucleusNPhi == 0) o->NucleusNPhi = 2*o->NucleusLmax + 16;
    if (o->NucleusNTheta > DiracMaxQuadraturePoints)
      return luaL_error(L, "Error in %s: NucleusNTheta = %d, at most %d (Gauss-Legendre in cos theta)\n", who, (int) o->NucleusNTheta, DiracMaxQuadraturePoints);
  }
  if (nucleus->model == NucleusGaussian || nucleus->model == NucleusFermi || nucleus->model == NucleusUniform)
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
  }
  // Function3D: rho(x, y, z) on the product grid (Gauss-Legendre in cos theta, uniform in
  // phi) at every quadrature point, projected on S_LM (NucleusMultipole.cpp); the spherical
  // part rho_00/sqrt(4 pi) is the nucleus of the atom
  if (nucleus->model == NucleusFunction3D)
  {
    uint32_t nT = o->NucleusNTheta, nF = o->NucleusNPhi, nP = basis->nPoint, i, j, nLM = NucleusMultipoleCount(o->NucleusLmax);
    double *xT = (double *) malloc(nT*sizeof(double)), *wT = (double *) malloc(nT*sizeof(double));
    double *rho3 = (double *) malloc(nP*nT*nF*sizeof(double)), *rhoLM = (double *) malloc(nLM*nP*sizeof(double));
    GaussLegendreNodesWeights(nT, -1., 1., xT, wT);
    for (q = 0; q < nP; q++)
      for (i = 0; i < nT; i++)
        for (j = 0; j < nF; j++)
        {
          double r = basis->rq[q], st = sqrt(fmax(0., 1. - xT[i]*xT[i])), ph = 2.*M_PI*j/nF;
          lua_rawgeti(L, LUA_REGISTRYINDEX, o->Density);
          lua_pushnumber(L, r*st*cos(ph));
          lua_pushnumber(L, r*st*sin(ph));
          lua_pushnumber(L, r*xT[i]);
          lua_call(L, 3, 1);
          if (lua_type(L, -1) != LUA_TNUMBER)
          {
            free(xT); free(wT); free(rho3); free(rhoLM);
            RadialBasisFree(basis);
            return luaL_error(L, "Error in %s: Density(%f, %f, %f) is not a number\n", who, r*st*cos(ph), r*st*sin(ph), r*xT[i]);
          }
          rho3[(q*nT + i)*nF + j] = lua_tonumber(L, -1);
          lua_pop(L, 1);
        }
    luaL_unref(L, LUA_REGISTRYINDEX, o->Density);
    NucleusMultipoleProject(rho3, nP, nT, nF, xT, wT, o->NucleusLmax, rhoLM);
    rho = (double *) malloc(nP*sizeof(double));
    for (q = 0; q < nP; q++) rho[q] = rhoLM[q]/sqrt(4.*M_PI);
    *charge = NucleusCharge(basis, rho);
    if (*charge > 0.)
    {
      for (i = 0; i < nLM*nP; i++) rhoLM[i] *= Z/(*charge);
      if (o->multipoleMoment != NULL)
        for (i = 0; i < nLM; i++)
        {
          int Lq, Mq;
          NucleusMultipoleLM(i, &Lq, &Mq);
          o->multipoleMoment[i] = NucleusMultipoleMoment(basis, Lq, rhoLM + i*nP);
        }
      if (o->multipoleV != NULL && NucleusMultipolePotentials(basis, o->NucleusLmax, rhoLM, o->multipoleV))
      {
        free(xT); free(wT); free(rho3); free(rhoLM); free(rho);
        RadialBasisFree(basis);
        return luaL_error(L, "Error in %s: the multipole potentials failed\n", who);
      }
    }
    free(xT); free(wT); free(rho3); free(rhoLM);
  }
  if (nucleus->model == NucleusFunction || nucleus->model == NucleusFunction3D)
  {
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
