/*   calculates < psi | psi > or < psi | O | psi >
 
 */
static int LuaBraKet(lua_State * L)
{
  uint32_t len = (uint32_t)lua_gettop(L);
  if (len==2)
  {
    WaveFunctionType *psiA, *psiB;
    Complex16 res;
    psiA = (WaveFunctionType*) luaL_checkudata(L, 1, "Wavefunction_Type");
    psiB = (WaveFunctionType*) luaL_checkudata(L, 2, "Wavefunction_Type");
    if(QWaveFunctionSameBasis(psiA, psiB))
    {
      WaveFunctionDotProductSameBasis(*psiA, *psiB, &(res.Re), &(res.Im));
    }
    else
    {
      WaveFunctionDotProduct(*psiA, *psiB, &(res.Re), &(res.Im));
    }
    if(Abs(res.Im)<10*(DBL_EPSILON*Abs(res.Re)+DBL_MIN))
    {
      lua_pushnumber(L, res.Re);
    }
    else
    {
      Complex16 *temp;
      temp = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
      temp->Re = res.Re;
      temp->Im = res.Im;
      luaL_setmetatable(L, "Complex_Type");
    }
    return 1;
  }
  else if (len==3)
  {
    WaveFunctionType *psiA, *psiB;
    OperatorType *O;
    Complex16 res;
    psiA = (WaveFunctionType*) luaL_checkudata(L, 1, "Wavefunction_Type");
    O    = (OperatorType*)     luaL_checkudata(L, 2, "Operator_Type");
    psiB = (WaveFunctionType*) luaL_checkudata(L, 3, "Wavefunction_Type");
    if(QWaveFunctionSame(psiA, psiB))
    {
      PsiOperatorPsi(*O, *psiA, &(res.Re), &(res.Im));
    }
    else
    {
      PsiPrimeOperatorPsi(*O, *psiA, *psiB, &(res.Re), &(res.Im));
    }
    if(Abs(res.Im)<10*(DBL_EPSILON*Abs(res.Re)+DBL_MIN))
    {
      lua_pushnumber(L, res.Re);
    }
    else
    {
      Complex16 *temp;
      temp = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
      temp->Re = res.Re;
      temp->Im = res.Im;
      luaL_setmetatable(L, "Complex_Type");
    }
    return 1;
  }
  return 0;
}