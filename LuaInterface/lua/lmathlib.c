/*
** $Id: lmathlib.c,v 1.83.1.1 2013/04/12 18:48:47 roberto Exp $
** Standard mathematical library
** See Copyright Notice in lua.h
*/


#include <stdlib.h>
#include <math.h>

#define lmathlib_c
#define LUA_LIB

#include "lua.h"

#include "lauxlib.h"
#include "lualib.h"


#undef PI
#define PI	((lua_Number)(3.1415926535897932384626433832795))
#define RADIANS_PER_DEGREE	((lua_Number)(PI/180.0))



static int math_abs (lua_State *L)
{
  if (luaL_testudata(L, 1, "Complex_Type") != NULL)
  {
    Complex16 *a = (Complex16 *) luaL_checkudata(L, 1, "Complex_Type");
    lua_pushnumber(L, sqrt(Sqr(a->Im)+Sqr(a->Re)));
  }
  else
  {
    lua_pushnumber(L, l_mathop(fabs)(luaL_checknumber(L, 1)));
  }
  return 1;
}

static int math_sin (lua_State *L)
{
  if (luaL_testudata(L, 1, "Complex_Type") != NULL)
  {
    Complex16 *res;
    Complex16 *a = (Complex16 *) luaL_checkudata(L, 1, "Complex_Type");
    res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
    res->Re = sin(a->Re) * cosh(a->Im);
    res->Im = cos(a->Re) * sinh(a->Im);
    luaL_setmetatable(L, "Complex_Type");
  }
  else
  {
    lua_pushnumber(L, l_mathop(sin)(luaL_checknumber(L, 1)));
  }
  return 1;
}

static int math_sinh (lua_State *L) {
  if (luaL_testudata(L, 1, "Complex_Type") != NULL)
  {
    Complex16 *res;
    Complex16 *a = (Complex16 *) luaL_checkudata(L, 1, "Complex_Type");
    res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
    res->Re = sinh(a->Re) * cos(a->Im);
    res->Im = cosh(a->Re) * sin(a->Im);
    luaL_setmetatable(L, "Complex_Type");
  }
  else
  {
    lua_pushnumber(L, l_mathop(sinh)(luaL_checknumber(L, 1)));
  }
  return 1;
}

static int math_cos (lua_State *L) {
  if (luaL_testudata(L, 1, "Complex_Type") != NULL)
  {
    Complex16 *res;
    Complex16 *a = (Complex16 *) luaL_checkudata(L, 1, "Complex_Type");
    res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
    res->Re = cos(a->Re) * cosh(a->Im);
    res->Im =-sin(a->Re) * sinh(a->Im);
    luaL_setmetatable(L, "Complex_Type");
  }
  else
  {
    lua_pushnumber(L, l_mathop(cos)(luaL_checknumber(L, 1)));
  }
  return 1;
}

static int math_cosh (lua_State *L) {
  if (luaL_testudata(L, 1, "Complex_Type") != NULL)
  {
    Complex16 *res;
    Complex16 *a = (Complex16 *) luaL_checkudata(L, 1, "Complex_Type");
    res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
    res->Re = cosh(a->Re) * cos(a->Im);
    res->Im = sinh(a->Re) * sin(a->Im);
    luaL_setmetatable(L, "Complex_Type");
  }
  else
  {
    lua_pushnumber(L, l_mathop(cosh)(luaL_checknumber(L, 1)));
  }
  return 1;
}

static int math_tan (lua_State *L) {
  if (luaL_testudata(L, 1, "Complex_Type") != NULL)
  {
    Complex16 *res;
    Complex16 *a = (Complex16 *) luaL_checkudata(L, 1, "Complex_Type");
    res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
    res->Re = sin(2*a->Re)  / (cos(2*a->Re)+cosh(2*a->Im));
    res->Im = sinh(2*a->Im) / (cos(2*a->Re)+cosh(2*a->Im));
    luaL_setmetatable(L, "Complex_Type");
  }
  else
  {
    lua_pushnumber(L, l_mathop(tan)(luaL_checknumber(L, 1)));
  }
  return 1;
}

static int math_tanh (lua_State *L) {
  if (luaL_testudata(L, 1, "Complex_Type") != NULL)
  {
    Complex16 *res;
    Complex16 *a = (Complex16 *) luaL_checkudata(L, 1, "Complex_Type");
    res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
    res->Re = sinh(2*a->Re) / (cosh(2*a->Re)+cos(2*a->Im));
    res->Im = sin(2*a->Im)  / (cosh(2*a->Re)+cos(2*a->Im));
    luaL_setmetatable(L, "Complex_Type");
  }
  else
  {
    lua_pushnumber(L, l_mathop(tanh)(luaL_checknumber(L, 1)));
  }
  return 1;
}

static int math_asin (lua_State *L) {
  lua_pushnumber(L, l_mathop(asin)(luaL_checknumber(L, 1)));
  return 1;
}

static int math_acos (lua_State *L) {
  lua_pushnumber(L, l_mathop(acos)(luaL_checknumber(L, 1)));
  return 1;
}

static int math_atan (lua_State *L) {
  lua_pushnumber(L, l_mathop(atan)(luaL_checknumber(L, 1)));
  return 1;
}

static int math_atan2 (lua_State *L) {
  lua_pushnumber(L, l_mathop(atan2)(luaL_checknumber(L, 1),
                                luaL_checknumber(L, 2)));
  return 1;
}

static int math_ceil (lua_State *L) {
  lua_pushnumber(L, l_mathop(ceil)(luaL_checknumber(L, 1)));
  return 1;
}

static int math_floor (lua_State *L) {
  lua_pushnumber(L, l_mathop(floor)(luaL_checknumber(L, 1)));
  return 1;
}

static int math_fmod (lua_State *L) {
  lua_pushnumber(L, l_mathop(fmod)(luaL_checknumber(L, 1),
                               luaL_checknumber(L, 2)));
  return 1;
}

static int math_modf (lua_State *L) {
  lua_Number ip;
  lua_Number fp = l_mathop(modf)(luaL_checknumber(L, 1), &ip);
  lua_pushnumber(L, ip);
  lua_pushnumber(L, fp);
  return 2;
}

static int math_sqrt (lua_State *L) {
  double s=luaL_checknumber(L, 1);
  if(s<0) luaL_error(L, "Calculating sqrt(%d). Can not take the root of a negative number\n",s);
  lua_pushnumber(L, l_mathop(sqrt)(s));
  return 1;
}

static int math_pow (lua_State *L) {
  lua_Number x = luaL_checknumber(L, 1);
  lua_Number y = luaL_checknumber(L, 2);
  lua_pushnumber(L, l_mathop(pow)(x, y));
  return 1;
}

static int math_log (lua_State *L)
{
  if (luaL_testudata(L, 1, "Complex_Type") != NULL)
  {
    Complex16 *res;
    Complex16 *a = (Complex16 *) luaL_checkudata(L, 1, "Complex_Type");
    res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
    luaL_setmetatable(L, "Complex_Type");
    res->Re = 0.5*log(Sqr(a->Re)+Sqr(a->Im));
    res->Im = atan2(a->Im,a->Re);
    if (!lua_isnoneornil(L, 2))
    {
      lua_Number base = luaL_checknumber(L, 2);
      res->Re /= l_mathop(log)(base);
      res->Im /= l_mathop(log)(base);
    }
  }
  else
  {
    lua_Number x = luaL_checknumber(L, 1);
    lua_Number res;
    if (lua_isnoneornil(L, 2))
      res = l_mathop(log)(x);
    else
    {
      lua_Number base = luaL_checknumber(L, 2);
      if (base == (lua_Number)10.0) res = l_mathop(log10)(x);
      else res = l_mathop(log)(x)/l_mathop(log)(base);
    }
    lua_pushnumber(L, res);
  }
  return 1;
}

#if defined(LUA_COMPAT_LOG10)
static int math_log10 (lua_State *L)
{
  if (luaL_testudata(L, 1, "Complex_Type") != NULL)
  {
    Complex16 *res;
    Complex16 *a = (Complex16 *) luaL_checkudata(L, 1, "Complex_Type");
    res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
    luaL_setmetatable(L, "Complex_Type");
    res->Re = 0.5*log10(Sqr(a->Re)+Sqr(a->Im));
    res->Im = atan2(a->Im,a->Re)/log(10.0);
  }
  else
  {
    lua_pushnumber(L, l_mathop(log10)(luaL_checknumber(L, 1)));
  }
  return 1;
}
#endif

static int math_exp (lua_State *L)
{
  if (luaL_testudata(L, 1, "Complex_Type") != NULL)
  {
    Complex16 *res;
    Complex16 *a = (Complex16 *) luaL_checkudata(L, 1, "Complex_Type");
    res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
    luaL_setmetatable(L, "Complex_Type");
    res->Re = exp(a->Re)*cos(a->Im);
    res->Im = exp(a->Re)*sin(a->Im);
  }
  else if (luaL_testudata(L, 1, "Operator_Type") != NULL)
  {
    return LuaOperatorExp(L);
  }
  else
  {
    lua_pushnumber(L, l_mathop(exp)(luaL_checknumber(L, 1)));
  }
  return 1;
}

static int math_erf (lua_State *L)
{
  if (luaL_testudata(L, 1, "Complex_Type") != NULL)
  {
    Complex16 *res;
    Complex16 *a = (Complex16 *) luaL_checkudata(L, 1, "Complex_Type");
    res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
    luaL_setmetatable(L, "Complex_Type");
    *res = Cerf(*a);
  }
  else
  {
    lua_pushnumber(L, l_mathop(erf)(luaL_checknumber(L, 1)));
  }
  return 1;
}

static int math_erfc (lua_State *L)
{
  if (luaL_testudata(L, 1, "Complex_Type") != NULL)
  {
    Complex16 *res;
    Complex16 *a = (Complex16 *) luaL_checkudata(L, 1, "Complex_Type");
    res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
    luaL_setmetatable(L, "Complex_Type");
    *res = Cerfc(*a);
  }
  else
  {
    lua_pushnumber(L, l_mathop(erfc)(luaL_checknumber(L, 1)));
  }
  return 1;
}

static int math_deg (lua_State *L) {
  lua_pushnumber(L, luaL_checknumber(L, 1)/RADIANS_PER_DEGREE);
  return 1;
}

static int math_rad (lua_State *L) {
  lua_pushnumber(L, luaL_checknumber(L, 1)*RADIANS_PER_DEGREE);
  return 1;
}

static int math_frexp (lua_State *L) {
  int e;
  lua_pushnumber(L, l_mathop(frexp)(luaL_checknumber(L, 1), &e));
  lua_pushinteger(L, e);
  return 2;
}

static int math_ldexp (lua_State *L) {
  lua_Number x = luaL_checknumber(L, 1);
  int ep = luaL_checkint(L, 2);
  lua_pushnumber(L, l_mathop(ldexp)(x, ep));
  return 1;
}



static int math_min (lua_State *L) {
  int n = lua_gettop(L);  /* number of arguments */
  lua_Number dmin = luaL_checknumber(L, 1);
  int i;
  for (i=2; i<=n; i++) {
    lua_Number d = luaL_checknumber(L, i);
    if (d < dmin)
      dmin = d;
  }
  lua_pushnumber(L, dmin);
  return 1;
}


static int math_max (lua_State *L) {
  int n = lua_gettop(L);  /* number of arguments */
  lua_Number dmax = luaL_checknumber(L, 1);
  int i;
  for (i=2; i<=n; i++) {
    lua_Number d = luaL_checknumber(L, i);
    if (d > dmax)
      dmax = d;
  }
  lua_pushnumber(L, dmax);
  return 1;
}


static int math_random (lua_State *L) {
  /* the `%' avoids the (rare) case of r==1, and is needed also because on
     some systems (SunOS!) `rand()' may return a value larger than RAND_MAX */
  lua_Number r = (lua_Number)(rand()%RAND_MAX) / (lua_Number)RAND_MAX;
  switch (lua_gettop(L)) {  /* check number of arguments */
    case 0: {  /* no arguments */
      lua_pushnumber(L, r);  /* Number between 0 and 1 */
      break;
    }
    case 1: {  /* only upper limit */
      lua_Number u = luaL_checknumber(L, 1);
      luaL_argcheck(L, (lua_Number)1.0 <= u, 1, "interval is empty");
      lua_pushnumber(L, l_mathop(floor)(r*u) + (lua_Number)(1.0));  /* [1, u] */
      break;
    }
    case 2: {  /* lower and upper limits */
      lua_Number l = luaL_checknumber(L, 1);
      lua_Number u = luaL_checknumber(L, 2);
      luaL_argcheck(L, l <= u, 2, "interval is empty");
      lua_pushnumber(L, l_mathop(floor)(r*(u-l+1)) + l);  /* [l, u] */
      break;
    }
    default: return luaL_error(L, "wrong number of arguments");
  }
  return 1;
}


static int math_randomseed (lua_State *L) {
  srand(luaL_checkunsigned(L, 1));
  (void)rand(); /* discard first value to avoid undesirable correlations */
  return 0;
}



/* on stack l,m,theta,phi
 *  returns a complex number
 */
static int math_SphericalHarmonicY(lua_State * L)
{
  int l,m;
  double theta,phi;
  Complex16 *res;
  if(lua_gettop(L)!=4) return luaL_error(L, "math.SphericalHarmonicsY expects 4 parameters l,m,theta,phi");
  l = (int) luaL_checknumber(L, 1);
  if((l<0)||(l>12)) return luaL_error(L, "Spherical Harmonics only implemented for 0<=l<=12 l=%d",l);
  m = (int) luaL_checknumber(L, 2);
  if(Abs(m)>l) return luaL_error(L, "Spherical Harmonics only defined for -l<m<l l=%d m=%d",l,m);
  theta = (double) luaL_checknumber(L,3);
  phi   = (double) luaL_checknumber(L,4);
  res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
  res->Re = ReSphericalHarmonicY(l, m, theta, phi);
  res->Im = ImSphericalHarmonicY(l, m, theta, phi);
  luaL_setmetatable(L, "Complex_Type");
  return 1;
}
static int math_SphericalHarmonicC(lua_State * L)
{
  int l,m;
  double theta,phi;
  Complex16 *res;
  if(lua_gettop(L)!=4) return luaL_error(L, "math.SphericalHarmonicsC expects 4 parameters l,m,theta,phi");
  l = (int) luaL_checknumber(L, 1);
  if((l<0)||(l>12)) return luaL_error(L, "Spherical Harmonics only implemented for 0<=l<=12 l=%d",l);
  m = (int) luaL_checknumber(L, 2);
  if(Abs(m)>l) return luaL_error(L, "Spherical Harmonics only defined for -l<m<l l=%d m=%d",l,m);
  theta = (double) luaL_checknumber(L,3);
  phi   = (double) luaL_checknumber(L,4);
  res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
  res->Re = sqrt(4*Pi/((double)(2*l+1))) * ReSphericalHarmonicY(l, m, theta, phi);
  res->Im = sqrt(4*Pi/((double)(2*l+1))) * ImSphericalHarmonicY(l, m, theta, phi);
  luaL_setmetatable(L, "Complex_Type");
  return 1;
}

/* on stack l,m,theta,phi
 *  returns a real number
 */
static int math_KubicHarmonicK(lua_State * L)
{
  int l,m;
  double theta,phi;
  double res;
  if(lua_gettop(L)!=4) return luaL_error(L, "math.SphericalHarmonicsY expects 4 parameters l,m,theta,phi");
  l = (int) luaL_checknumber(L, 1);
  if((l<0)||(l>12)) return luaL_error(L, "Spherical Harmonics only implemented for 0<=l<=6 l=%d",l);
  m = (int) luaL_checknumber(L, 2);
  if(Abs(m)>l) return luaL_error(L, "Spherical Harmonics only defined for -l<m<l l=%d m=%d",l,m);
  theta = (double) luaL_checknumber(L,3);
  phi   = (double) luaL_checknumber(L,4);
  res=CubicHarmonicK(l, m, theta, phi);
  lua_pushnumber(L, res);
  return 1;
}
static int math_KubicHarmonicE(lua_State * L)
{
  int l,m;
  double theta,phi;
  double res;
  if(lua_gettop(L)!=4) return luaL_error(L, "math.SphericalHarmonicsY expects 4 parameters l,m,theta,phi");
  l = (int) luaL_checknumber(L, 1);
  if((l<0)||(l>12)) return luaL_error(L, "Spherical Harmonics only implemented for 0<=l<=6 l=%d",l);
  m = (int) luaL_checknumber(L, 2);
  if(Abs(m)>l) return luaL_error(L, "Spherical Harmonics only defined for -l<m<l l=%d m=%d",l,m);
  theta = (double) luaL_checknumber(L,3);
  phi   = (double) luaL_checknumber(L,4);
  res=sqrt(4*Pi/((double)(2*l+1))) * CubicHarmonicK(l, m, theta, phi);
  lua_pushnumber(L, res);
  return 1;
}

/* on stack l,m,theta,phi
 *  returns a real number
 */
static int math_TesseralHarmonicZ(lua_State * L)
{
  int l,m;
  double theta,phi;
  double res;
  if(lua_gettop(L)!=4) return luaL_error(L, "math.SphericalHarmonicsY expects 4 parameters l,m,theta,phi");
  l = (int) luaL_checknumber(L, 1);
  if((l<0)||(l>12)) return luaL_error(L, "Spherical Harmonics only implemented for 0<l<12 l=%d",l);
  m = (int) luaL_checknumber(L, 2);
  if(Abs(m)>l) return luaL_error(L, "Spherical Harmonics only defined for -l<m<l l=%d m=%d",l,m);
  theta = (double) luaL_checknumber(L,3);
  phi   = (double) luaL_checknumber(L,4);
  res=TesseralHarmonicZ(l, m, theta, phi);
  lua_pushnumber(L, res);
  return 1;
}
static int math_TesseralHarmonicD(lua_State * L)
{
  int l,m;
  double theta,phi;
  double res;
  if(lua_gettop(L)!=4) return luaL_error(L, "math.SphericalHarmonicsY expects 4 parameters l,m,theta,phi");
  l = (int) luaL_checknumber(L, 1);
  if((l<0)||(l>12)) return luaL_error(L, "Spherical Harmonics only implemented for 0<l<12 l=%d",l);
  m = (int) luaL_checknumber(L, 2);
  if(Abs(m)>l) return luaL_error(L, "Spherical Harmonics only defined for -l<m<l l=%d m=%d",l,m);
  theta = (double) luaL_checknumber(L,3);
  phi   = (double) luaL_checknumber(L,4);
  res=sqrt(4*Pi/((double)(2*l+1))) * TesseralHarmonicZ(l, m, theta, phi);
  lua_pushnumber(L, res);
  return 1;
}

/* on stack l,x
 *  returns a real number
 */
static int math_SphericalBesselJ(lua_State * L)
{
  int l;
  double x;
  double res;
  if(lua_gettop(L)!=2) return luaL_error(L, "math.SphericalBesselJ expects 2 parameters l,x");
  l = (int) luaL_checknumber(L, 1);
  x = (double) luaL_checknumber(L, 2);
  res = SphericalBesselJ(x, l);
  lua_pushnumber(L, res);
  return 1;
}


/* on stack n,l,x
 *  returns a real number
 */
static int math_LaguerreL(lua_State * L)
{
  int n,l;
  double x;
  double res;
  if(lua_gettop(L)!=3) return luaL_error(L, "math.LaguerreL expects 3 parameters n,l,x");
  n = (int) luaL_checknumber(L, 1);
  l = (int) luaL_checknumber(L, 2);
  x = (double) luaL_checknumber(L, 3);
  res = LaguerrePolynomial(n,l,x);
  lua_pushnumber(L, res);
  return 1;
}

/* on stack w,GFWHM,LFWHM
 *  returns a real number
 */
static int math_Voight(lua_State * L)
{
  double w,GFWHM,LFWHM;
  double res;
  if(lua_gettop(L)!=3) return luaL_error(L, "math.Voight expects 3 parameters w,GFWMH,LFWHM");
  w     = (double) luaL_checknumber(L, 1);
  GFWHM = (double) luaL_checknumber(L, 2);
  LFWHM = (double) luaL_checknumber(L, 3);
  res = Voight(w, GFWHM, LFWHM);
  lua_pushnumber(L, res);
  return 1;
}

/* on stack n
 *  returns a real number
 */
static int math_fact(lua_State * L)
{
  int n;
  double res;
  if(lua_gettop(L)!=1) return luaL_error(L, "math.fact expects 1 parameter n");
  n = (int) luaL_checknumber(L, 1);
  res = Factorial(n);
  lua_pushnumber(L, res);
  return 1;
}

/* on stack n,k
 *  returns a real number
 */
static int math_binomial(lua_State * L)
{
  uint32_t n,k;
  double res;
  if(lua_gettop(L)!=2) return luaL_error(L, "math.Binomial expects 2 parameters n and k");
  n = (uint32_t) luaL_checknumber(L, 1);
  k = (uint32_t) luaL_checknumber(L, 2);
  res = BinomialCoefficientDouble(n, k);
  lua_pushnumber(L, res);
  return 1;
}

static int math_ThreeJSymbol(lua_State *L)
{
  double l1,l2,l3,m1,m2,m3;
  double threeJ;
  int il1,il2,il3,im1,im2,im3;
  if(lua_gettop(L)!=3) return luaL_error(L, "math.ThreeJSymbol expects 6 parameters grouped in 3 lists of length 2:\n{l1,m1}, {l2,m2}, {l3,m3}");
  luaL_checktype(L, 1, LUA_TTABLE);
  if(lua_objlen(L,1)!=2) return luaL_error(L, "math.ThreeJSymbol expects 6 parameters grouped in 3 lists of length 2:\n{l1,m1}, {l2,m2}, {l3,m3}");
  luaL_checktype(L, 2, LUA_TTABLE);
  if(lua_objlen(L,2)!=2) return luaL_error(L, "math.ThreeJSymbol expects 6 parameters grouped in 3 lists of length 2:\n{l1,m1}, {l2,m2}, {l3,m3}");
  luaL_checktype(L, 3, LUA_TTABLE);
  if(lua_objlen(L,3)!=2) return luaL_error(L, "math.ThreeJSymbol expects 6 parameters grouped in 3 lists of length 2:\n{l1,m1}, {l2,m2}, {l3,m3}");
  lua_rawgeti(L,1,1);
  l1=(double)luaL_checknumber(L, -1);
  lua_pop(L,1);
  lua_rawgeti(L,1,2);
  m1=(double)luaL_checknumber(L, -1);
  lua_pop(L,1);
  lua_rawgeti(L,2,1);
  l2=(double)luaL_checknumber(L, -1);
  lua_pop(L,1);
  lua_rawgeti(L,2,2);
  m2=(double)luaL_checknumber(L, -1);
  lua_pop(L,1);
  lua_rawgeti(L,3,1);
  l3=(double)luaL_checknumber(L, -1);
  lua_pop(L,1);
  lua_rawgeti(L,3,2);
  m3=(double)luaL_checknumber(L, -1);
  lua_pop(L,1);
  if((Abs(round(l1)-l1)<10*DBL_EPSILON*Abs(l1)+DBL_MIN)&&
     (Abs(round(m1)-m1)<10*DBL_EPSILON*Abs(m1)+DBL_MIN)&&
     (Abs(round(l2)-l2)<10*DBL_EPSILON*Abs(l2)+DBL_MIN)&&
     (Abs(round(m2)-m2)<10*DBL_EPSILON*Abs(m2)+DBL_MIN)&&
     (Abs(round(l3)-l3)<10*DBL_EPSILON*Abs(l3)+DBL_MIN)&&
     (Abs(round(m3)-m3)<10*DBL_EPSILON*Abs(m3)+DBL_MIN))
  {
    il1 = (int)round(l1);
    im1 = (int)round(m1);
    il2 = (int)round(l2);
    im2 = (int)round(m2);
    il3 = (int)round(l3);
    im3 = (int)round(m3);
    threeJ=ThreeJSymbol(il1, im1, il2, im2, il3, im3);
    lua_pushnumber(L, threeJ);
  }
  else if((Abs(round(2*l1)-2*l1)<10*DBL_EPSILON*Abs(2*l1)+DBL_MIN)&&
          (Abs(round(2*m1)-2*m1)<10*DBL_EPSILON*Abs(2*m1)+DBL_MIN)&&
          (Abs(round(2*l2)-2*l2)<10*DBL_EPSILON*Abs(2*l2)+DBL_MIN)&&
          (Abs(round(2*m2)-2*m2)<10*DBL_EPSILON*Abs(2*m2)+DBL_MIN)&&
          (Abs(round(2*l3)-2*l3)<10*DBL_EPSILON*Abs(2*l3)+DBL_MIN)&&
          (Abs(round(2*m3)-2*m3)<10*DBL_EPSILON*Abs(2*m3)+DBL_MIN))
  {
    il1 = (int)round(2*l1);
    im1 = (int)round(2*m1);
    il2 = (int)round(2*l2);
    im2 = (int)round(2*m2);
    il3 = (int)round(2*l3);
    im3 = (int)round(2*m3);
    threeJ=ThreeJSymbolHalf(il1, im1, il2, im2, il3, im3);
    lua_pushnumber(L, threeJ);
  }
  else
  {
    return luaL_error(L, "ThreeJSymbol expects 6 angular momentum as input that need to be integer or half integer values");
  }
  return 1;
}


static int math_SlaterCoefficientC(lua_State *L)
{
  double l1,l2,l3,m1,m2,m3;
  double Ck;
  int il1,il2,il3,im1,im2,im3;
  if(lua_gettop(L)!=3) return luaL_error(L, "math.SlaterCoefficientC expects 6 parameters grouped in 3 lists of length 2:\n{l1,m1}, {l2,m2}, {l3,m3}");
  if(lua_istable(L, 1))
  {
    if(lua_objlen(L,1)!=2) return luaL_error(L, "math.SlaterCoefficientC expects 6 parameters grouped in 3 lists of length 2:\n{l1,m1}, {l2,m2}, {l3,m3}");
    luaL_checktype(L, 2, LUA_TTABLE);
    if(lua_objlen(L,2)!=2) return luaL_error(L, "math.SlaterCoefficientC expects 6 parameters grouped in 3 lists of length 2:\n{l1,m1}, {l2,m2}, {l3,m3}");
    luaL_checktype(L, 3, LUA_TTABLE);
    if(lua_objlen(L,3)!=2) return luaL_error(L, "math.SlaterCoefficientC expects 6 parameters grouped in 3 lists of length 2:\n{l1,m1}, {l2,m2}, {l3,m3}");
    lua_rawgeti(L,1,1);
    l1=(double)luaL_checknumber(L, -1);
    lua_pop(L,1);
    lua_rawgeti(L,1,2);
    m1=(double)luaL_checknumber(L, -1);
    lua_pop(L,1);
    lua_rawgeti(L,2,1);
    l2=(double)luaL_checknumber(L, -1);
    lua_pop(L,1);
    lua_rawgeti(L,2,2);
    m2=(double)luaL_checknumber(L, -1);
    lua_pop(L,1);
    lua_rawgeti(L,3,1);
    l3=(double)luaL_checknumber(L, -1);
    lua_pop(L,1);
    lua_rawgeti(L,3,2);
    m3=(double)luaL_checknumber(L, -1);
    lua_pop(L,1);
  }
  else
  {
    luaL_checktype(L, 2, LUA_TTABLE);
    if(lua_objlen(L,2)!=2) return luaL_error(L, "math.SlaterCoefficientC expects 5 parameters grouped in one integer and 2 lists of length 2:\nk, {l1,m1}, {l2,m2}");
    luaL_checktype(L, 3, LUA_TTABLE);
    if(lua_objlen(L,3)!=2) return luaL_error(L, "math.SlaterCoefficientC expects 5 parameters grouped in one integer and 2 lists of length 2:\nk, {l1,m1}, {l2,m2}");
    lua_rawgeti(L,2,1);
    l1=(double)luaL_checknumber(L, -1);
    lua_pop(L,1);
    lua_rawgeti(L,2,2);
    m1=(double)luaL_checknumber(L, -1);
    lua_pop(L,1);
    lua_rawgeti(L,3,1);
    l3=(double)luaL_checknumber(L, -1);
    lua_pop(L,1);
    lua_rawgeti(L,3,2);
    m3=(double)luaL_checknumber(L, -1);
    lua_pop(L,1);
    l2=(double)luaL_checknumber(L, 1);
    m2=m1-m3;
  }
  if((Abs(round(l1)-l1)<10*DBL_EPSILON*Abs(l1)+DBL_MIN)&&
     (Abs(round(m1)-m1)<10*DBL_EPSILON*Abs(m1)+DBL_MIN)&&
     (Abs(round(l2)-l2)<10*DBL_EPSILON*Abs(l2)+DBL_MIN)&&
     (Abs(round(m2)-m2)<10*DBL_EPSILON*Abs(m2)+DBL_MIN)&&
     (Abs(round(l3)-l3)<10*DBL_EPSILON*Abs(l3)+DBL_MIN)&&
     (Abs(round(m3)-m3)<10*DBL_EPSILON*Abs(m3)+DBL_MIN))
  {
    il1 = (int)round(l1);
    im1 = (int)round(m1);
    il2 = (int)round(l2);
    im2 = (int)round(m2);
    il3 = (int)round(l3);
    im3 = (int)round(m3);
    if(((-im1+im2+im3)!=0)||(il2>il1+il3)||(il2<Abs(il1-il3))||(((il1+il2+il3)&1)==1))
    {
      Ck=0;
    }
    else
    {
      Ck=SlaterCoefficientC(il2, il1, im1, il3, im3);
    }
    lua_pushnumber(L, Ck);
  }
  else
  {
    return luaL_error(L, "math.SlaterCoefficientC expects 6 angular momentum as input that need to be integer values");
  }
  return 1;
}

static int math_IntegerQ(lua_State *L)
{
  double nd;
  if(lua_type(L,-1)==LUA_TNUMBER)
  {
    nd=luaL_checknumber(L, -1);
    if(Abs(round(nd)-nd)<10*DBL_EPSILON*Abs(nd)+DBL_MIN)
    {
      lua_pushboolean(L, 1);
      return 1;
    }
  }
  lua_pushboolean(L, 0);
  return 1;
}

static int math_EvenQ(lua_State *L)
{
  double nd;
  if(lua_type(L,-1)==LUA_TNUMBER)
  {
    nd=luaL_checknumber(L, -1)/2.;
    if(Abs(round(nd)-nd)<(10*DBL_EPSILON*Abs(nd)+10*DBL_MIN))
    {
      lua_pushboolean(L, 1);
      return 1;
    }
  }
  else
  {
    printf("A\n");
  }
  lua_pushboolean(L, 0);
  return 1;
}

static int math_OddQ(lua_State *L)
{
  double nd;
  if(lua_type(L,-1)==LUA_TNUMBER)
  {
    nd=(1.+luaL_checknumber(L, -1))/2.;
    if(Abs(round(nd)-nd)<10*DBL_EPSILON*Abs(nd))
    {
      lua_pushboolean(L, 1);
    }
    return 1;
  }
  lua_pushboolean(L, 0);
  return 1;
}

static int math_Dot(lua_State *L)
{
  CompactMatrixType *A,*B,C;
  luaL_checktype(L, 1, LUA_TTABLE);
  luaL_checktype(L, 2, LUA_TTABLE);
  A = (CompactMatrixType *) luaL_checkCompactMatrix(L, 1);
  B = (CompactMatrixType *) luaL_checkCompactMatrix(L, 2);
  CompactMatrixMultiply('N', 'N', A, B, &C);
  LuaPushTable_CompactMatrixType(L, &C);
  FreeCompactMatrix(A);
  free(A);
  FreeCompactMatrix(B);
  free(B);
  FreeCompactMatrix(&C);
  return 1;
}

static int math_MatrixTranspose(lua_State *L)
{
  CompactMatrixType *A;
  luaL_checktype(L, 1, LUA_TTABLE);
  A = (CompactMatrixType *) luaL_checkCompactMatrix(L, 1);
  CompactMatrixTranspose(A);
  LuaPushTable_CompactMatrixType(L, A);
  FreeCompactMatrix(A);
  free(A);
  return 1;
}

static int math_MatrixConjugateTranspose(lua_State *L)
{
  CompactMatrixType *A;
  luaL_checktype(L, 1, LUA_TTABLE);
  A = (CompactMatrixType *) luaL_checkCompactMatrix(L, 1);
  CompactMatrixConjugateTranspose(A);
  LuaPushTable_CompactMatrixType(L, A);
  FreeCompactMatrix(A);
  free(A);
  return 1;
}

static int math_MatrixConjugate(lua_State *L)
{
  CompactMatrixType *A;
  luaL_checktype(L, 1, LUA_TTABLE);
  A = (CompactMatrixType *) luaL_checkCompactMatrix(L, 1);
  CompactMatrixConjugate(A);
  LuaPushTable_CompactMatrixType(L, A);
  FreeCompactMatrix(A);
  free(A);
  return 1;
}

static const luaL_Reg mathlib[] = {
  {"abs",   math_abs},{"Abs",   math_abs},
  {"acos",  math_acos},{"Acos",  math_acos},
  {"asin",  math_asin},{"Asin",  math_asin},
  {"atan2", math_atan2},{"Atan2", math_atan2},
  {"atan",  math_atan},{"Atan",  math_atan},
  {"ceil",  math_ceil},{"Ceil",  math_ceil},
  {"cosh",   math_cosh},{"Cosh",   math_cosh},
  {"cos",   math_cos},{"Cos",   math_cos},
  {"deg",   math_deg},{"Deg",   math_deg},
  {"exp",   math_exp},{"Exp",   math_exp},
  {"erf",   math_erf},{"Erf",   math_erf},
  {"erfc",   math_erfc},{"Erfc",   math_erfc},
  {"floor", math_floor},{"Floor", math_floor},
  {"fmod",   math_fmod},{"Fmod",   math_fmod},
  {"frexp", math_frexp},{"Frexp", math_frexp},
  {"ldexp", math_ldexp},{"Ldexp", math_ldexp},
#if defined(LUA_COMPAT_LOG10)
  {"log10", math_log10},{"Log10", math_log10},
#endif
  {"log",   math_log},{"Log",   math_log},
  {"max",   math_max},{"Max",   math_max},
  {"min",   math_min},{"Min",   math_min},
  {"modf",   math_modf},{"Modf",   math_modf},
  {"pow",   math_pow},{"Pow",   math_pow},
  {"rad",   math_rad},{"Rad",   math_rad},
  {"random",     math_random},{"Random",     math_random},
  {"randomseed", math_randomseed},{"RandomSeed", math_randomseed},
  {"sinh",   math_sinh},{"Sinh",   math_sinh},
  {"sin",   math_sin},{"Sin",   math_sin},
  {"sqrt",  math_sqrt},{"Sqrt",  math_sqrt},
  {"tanh",   math_tanh},{"Tanh",   math_tanh},
  {"tan",   math_tan},{"Tan",   math_tan},
  {"fact",   math_fact},{"Fact",   math_fact},{"factorial",   math_fact},{"Factorial",   math_fact},
  {"Binomial",   math_binomial},{"BinomialCoefficient",   math_binomial},
  {"ThreeJ",   math_ThreeJSymbol},{"ThreeJSymbol",   math_ThreeJSymbol},
  {"SlaterCoefficientC", math_SlaterCoefficientC},
  {"SphericalHarmonicY", math_SphericalHarmonicY},
  {"SphericalHarmonicC", math_SphericalHarmonicC},
  {"TesseralHarmonicZ",  math_TesseralHarmonicZ},
  {"TesseralHarmonicD",  math_TesseralHarmonicD},
  {"KubicHarmonicK",     math_KubicHarmonicK},
  {"KubicHarmonicE",     math_KubicHarmonicE},
  {"SphericalBesselJ",   math_SphericalBesselJ},
  {"LaguerreL",          math_LaguerreL},
  {"Voight",             math_Voight},
  {"IntegerQ",           math_IntegerQ},
  {"EvenQ",              math_EvenQ},
  {"OddQ",               math_OddQ},
  {"Dot",                math_Dot},
  {"MatrixTranspose",    math_MatrixTranspose},
  {"MatrixConjugateTranspose",    math_MatrixConjugateTranspose},
  {"MatrixConjugate",    math_MatrixConjugate},
  {NULL, NULL}
};

/*
** Open math library
*/
LUAMOD_API int luaopen_math (lua_State *L) {
  luaL_newlib(L, mathlib);
  lua_pushnumber(L, PI);
  lua_setfield(L, -2, "pi");
  lua_pushnumber(L, PI);
  lua_setfield(L, -2, "Pi");
  lua_pushnumber(L, HUGE_VAL);
  lua_setfield(L, -2, "huge");
  lua_pushnumber(L, HUGE_VAL);
  lua_setfield(L, -2, "Huge");
  return 1;
}

