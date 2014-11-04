/*  define the ssp.NewOperator functions
 *  depending on the first argument of ssp.NewOperator different types of standard operators are created. Here the implementation
 */

static int LuaOperatorNewCoulomb(lua_State * L)
{
  OperatorType *op;
  uint32_t Ntemp,NShell,l,lA,lB;
  uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t len = (uint32_t)lua_gettop(L);
  if(len==4) /* creates the coulomb operator in a single shell of angular momentum L */
  {
    uint16_t *LindexUp, *LindexDn;
    double *F;
    LindexUp = luaL_checkarray_unit16_t(L, 2, &NShell);
    LindexDn = luaL_checkarray_unit16_t(L, 3, &Ntemp);
    if(NShell!=Ntemp)
    {
      return luaL_error(L, "Error in LuaOperatorNewCoulomb, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
    }
    l=(NShell-1)/2;
    if(2*l+1!=NShell)
    {
      return luaL_error(L, "Error in LuaOperatorNewCoulomb, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
    }
    F = luaL_checkarray_double(L, 4, &Ntemp);
    if(Ntemp!=l+1)
    {
      return luaL_error(L, "Error in LuaOperatorNewCoulomb, Lenght of F (%d) is not equal to l+1\n",Ntemp);
    }
    op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    CreateOperatorCoulombSingleShell(NFermion, l, LindexUp, LindexDn, F, op);
    luaL_setmetatable(L, "ssp.operator");
    free(LindexUp);
    free(LindexDn);
    free(F);
    return 1;
  }
  else if (len==7)
  {
    uint16_t *LindexUpA, *LindexDnA, *LindexUpB, *LindexDnB;
    double *F,*G;
    LindexUpA = luaL_checkarray_unit16_t(L, 2, &NShell);
    LindexDnA = luaL_checkarray_unit16_t(L, 3, &Ntemp);
    if(NShell!=Ntemp)
    {
      return luaL_error(L, "Error in LuaOperatorNewCoulomb, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
    }
    lA=(NShell-1)/2;
    if(2*lA+1!=NShell)
    {
      return luaL_error(L, "Error in LuaOperatorNewCoulomb, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
    }
    LindexUpB = luaL_checkarray_unit16_t(L, 4, &NShell);
    LindexDnB = luaL_checkarray_unit16_t(L, 5, &Ntemp);
    if(NShell!=Ntemp)
    {
      return luaL_error(L, "Error in LuaOperatorNewCoulomb, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
    }
    lB=(NShell-1)/2;
    if(2*lB+1!=NShell)
    {
      return luaL_error(L, "Error in LuaOperatorNewCoulomb, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
    }
    F = luaL_checkarray_double(L, 6, &Ntemp);
    if(Ntemp!=(Min(lA,lB)+1))
    {
      return luaL_error(L, "Error in LuaOperatorNewCoulomb, Lenght of F (%d) is not equal to (Min(lA,lB)+1)1\n",Ntemp);
    }
    G = luaL_checkarray_double(L, 7, &Ntemp);
    if(Ntemp!=(Min(lA,lB)+1))
    {
      return luaL_error(L, "Error in LuaOperatorNewCoulomb, Lenght of G (%d) is not equal to (Min(lA,lB)+1)1\n",Ntemp);
    }
    op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    CreateOperatorCoulombTwoShellsConstantOccupation(NFermion, lA, lB, LindexUpA, LindexDnA, LindexUpB, LindexDnB, F, G, op);
    luaL_setmetatable(L, "ssp.operator");
    free(LindexUpA);
    free(LindexDnA);
    free(LindexUpB);
    free(LindexDnB);
    free(F);
    free(G);
    return 1;
  }
  else
  {
    return luaL_error(L, "Error in LuaOperatorNewCoulomb, expected 4 or 7 arguments not %d",len);
  }
  return 1;
}


static int LuaOperatorNewSx(lua_State * L)
{
  OperatorType *op;
  uint32_t Ntemp;
  uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t NShell;
  uint16_t *LindexUp, *LindexDn;
  LindexUp = luaL_checkarray_unit16_t(L, 2, &NShell);
  LindexDn = luaL_checkarray_unit16_t(L, 3, &Ntemp);
  if(NShell!=Ntemp)
  {
    return luaL_error(L, "Error in LuaOperatorNewSx, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorSx(NFermion, NShell, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "ssp.operator");
  free(LindexUp);
  free(LindexDn);
  return 1;
}

static int LuaOperatorNewSy(lua_State * L)
{
  OperatorType *op;
  uint32_t Ntemp;
  uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t NShell;
  uint16_t *LindexUp, *LindexDn;
  LindexUp = luaL_checkarray_unit16_t(L, 2, &NShell);
  LindexDn = luaL_checkarray_unit16_t(L, 3, &Ntemp);
  if(NShell!=Ntemp)
  {
    return luaL_error(L, "Error in LuaOperatorNewSy, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorSy(NFermion, NShell, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "ssp.operator");
  free(LindexUp);
  free(LindexDn);
  return 1;
}

static int LuaOperatorNewSz(lua_State * L)
{
  OperatorType *op;
  uint32_t Ntemp;
  uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t NShell;
  uint16_t *LindexUp, *LindexDn;
  LindexUp = luaL_checkarray_unit16_t(L, 2, &NShell);
  LindexDn = luaL_checkarray_unit16_t(L, 3, &Ntemp);
  if(NShell!=Ntemp)
  {
    return luaL_error(L, "Error in LuaOperatorNewSz, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorSz(NFermion, NShell, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "ssp.operator");
  free(LindexUp);
  free(LindexDn);
  return 1;
}

static int LuaOperatorNewSsqr(lua_State * L)
{
  OperatorType *op;
  uint32_t Ntemp;
  uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t NShell;
  uint16_t *LindexUp, *LindexDn;
  LindexUp = luaL_checkarray_unit16_t(L, 2, &NShell);
  LindexDn = luaL_checkarray_unit16_t(L, 3, &Ntemp);
  if(NShell!=Ntemp)
  {
    return luaL_error(L, "Error in LuaOperatorNewSsqr, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorSsqr(NFermion, NShell, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "ssp.operator");
  free(LindexUp);
  free(LindexDn);
  return 1;
}

static int LuaOperatorNewSplus(lua_State * L)
{
  OperatorType *op;
  uint32_t Ntemp;
  uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t NShell;
  uint16_t *LindexUp, *LindexDn;
  LindexUp = luaL_checkarray_unit16_t(L, 2, &NShell);
  LindexDn = luaL_checkarray_unit16_t(L, 3, &Ntemp);
  if(NShell!=Ntemp)
  {
    return luaL_error(L, "Error in LuaOperatorNewSplus, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorSplus(NFermion, NShell, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "ssp.operator");
  free(LindexUp);
  free(LindexDn);
  return 1;
}

static int LuaOperatorNewSmin(lua_State * L)
{
  OperatorType *op;
  uint32_t Ntemp;
  uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t NShell;
  uint16_t *LindexUp, *LindexDn;
  LindexUp = luaL_checkarray_unit16_t(L, 2, &NShell);
  LindexDn = luaL_checkarray_unit16_t(L, 3, &Ntemp);
  if(NShell!=Ntemp)
  {
    return luaL_error(L, "Error in LuaOperatorNewSmin, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorSmin(NFermion, NShell, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "ssp.operator");
  free(LindexUp);
  free(LindexDn);
  return 1;
}



static int LuaOperatorNewLx(lua_State * L)
{
  OperatorType *op;
  uint32_t Ntemp;
  uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t NShell,l;
  uint16_t *LindexUp, *LindexDn;
  LindexUp = luaL_checkarray_unit16_t(L, 2, &NShell);
  LindexDn = luaL_checkarray_unit16_t(L, 3, &Ntemp);
  if(NShell!=Ntemp)
  {
    return luaL_error(L, "Error in LuaOperatorNewLx, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
  }
  l=(NShell-1)/2;
  if(2*l+1!=NShell)
  {
    return luaL_error(L, "Error in LuaOperatorNewLx, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorLx(NFermion, l, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "ssp.operator");
  free(LindexUp);
  free(LindexDn);
  return 1;
}

static int LuaOperatorNewLy(lua_State * L)
{
  OperatorType *op;
  uint32_t Ntemp;
  uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t NShell,l;
  uint16_t *LindexUp, *LindexDn;
  LindexUp = luaL_checkarray_unit16_t(L, 2, &NShell);
  LindexDn = luaL_checkarray_unit16_t(L, 3, &Ntemp);
  if(NShell!=Ntemp)
  {
    return luaL_error(L, "Error in LuaOperatorNewLy, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
  }
  l=(NShell-1)/2;
  if(2*l+1!=NShell)
  {
    return luaL_error(L, "Error in LuaOperatorNewLy, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorLy(NFermion, l, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "ssp.operator");
  free(LindexUp);
  free(LindexDn);
  return 1;
}

static int LuaOperatorNewLz(lua_State * L)
{
  OperatorType *op;
  uint32_t Ntemp;
  uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t NShell,l;
  uint16_t *LindexUp, *LindexDn;
  LindexUp = luaL_checkarray_unit16_t(L, 2, &NShell);
  LindexDn = luaL_checkarray_unit16_t(L, 3, &Ntemp);
  if(NShell!=Ntemp)
  {
    return luaL_error(L, "Error in LuaOperatorNewLz, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
  }
  l=(NShell-1)/2;
  if(2*l+1!=NShell)
  {
    return luaL_error(L, "Error in LuaOperatorNewLz, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorLz(NFermion, l, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "ssp.operator");
  free(LindexUp);
  free(LindexDn);
  return 1;
}

static int LuaOperatorNewLplus(lua_State * L)
{
  OperatorType *op;
  uint32_t Ntemp;
  uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t NShell,l;
  uint16_t *LindexUp, *LindexDn;
  LindexUp = luaL_checkarray_unit16_t(L, 2, &NShell);
  LindexDn = luaL_checkarray_unit16_t(L, 3, &Ntemp);
  if(NShell!=Ntemp)
  {
    return luaL_error(L, "Error in LuaOperatorNewLplus, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
  }
  l=(NShell-1)/2;
  if(2*l+1!=NShell)
  {
    return luaL_error(L, "Error in LuaOperatorNewLplus, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorLplus(NFermion, l, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "ssp.operator");
  free(LindexUp);
  free(LindexDn);
  return 1;
}

static int LuaOperatorNewLmin(lua_State * L)
{
  OperatorType *op;
  uint32_t Ntemp;
  uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t NShell,l;
  uint16_t *LindexUp, *LindexDn;
  LindexUp = luaL_checkarray_unit16_t(L, 2, &NShell);
  LindexDn = luaL_checkarray_unit16_t(L, 3, &Ntemp);
  if(NShell!=Ntemp)
  {
    return luaL_error(L, "Error in LuaOperatorNewLmin, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
  }
  l=(NShell-1)/2;
  if(2*l+1!=NShell)
  {
    return luaL_error(L, "Error in LuaOperatorNewLmin, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorLmin(NFermion, l, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "ssp.operator");
  free(LindexUp);
  free(LindexDn);
  return 1;
}

static int LuaOperatorNewLsqr(lua_State * L)
{
  OperatorType *op;
  uint32_t Ntemp;
  uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t NShell,l;
  uint16_t *LindexUp, *LindexDn;
  LindexUp = luaL_checkarray_unit16_t(L, 2, &NShell);
  LindexDn = luaL_checkarray_unit16_t(L, 3, &Ntemp);
  if(NShell!=Ntemp)
  {
    return luaL_error(L, "Error in LuaOperatorNewLsqr, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
  }
  l=(NShell-1)/2;
  if(2*l+1!=NShell)
  {
    return luaL_error(L, "Error in LuaOperatorNewLsqr, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorLsqr(NFermion, l, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "ssp.operator");
  free(LindexUp);
  free(LindexDn);
  return 1;
}

static int LuaOperatorNewJx(lua_State * L)
{
  OperatorType *op;
  uint32_t Ntemp;
  uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t NShell,l;
  uint16_t *LindexUp, *LindexDn;
  LindexUp = luaL_checkarray_unit16_t(L, 2, &NShell);
  LindexDn = luaL_checkarray_unit16_t(L, 3, &Ntemp);
  if(NShell!=Ntemp)
  {
    return luaL_error(L, "Error in LuaOperatorNewJx, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
  }
  l=(NShell-1)/2;
  if(2*l+1!=NShell)
  {
    return luaL_error(L, "Error in LuaOperatorNewJx, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorJx(NFermion, l, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "ssp.operator");
  free(LindexUp);
  free(LindexDn);
  return 1;
}

static int LuaOperatorNewJy(lua_State * L)
{
  OperatorType *op;
  uint32_t Ntemp;
  uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t NShell,l;
  uint16_t *LindexUp, *LindexDn;
  LindexUp = luaL_checkarray_unit16_t(L, 2, &NShell);
  LindexDn = luaL_checkarray_unit16_t(L, 3, &Ntemp);
  if(NShell!=Ntemp)
  {
    return luaL_error(L, "Error in LuaOperatorNewJy, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
  }
  l=(NShell-1)/2;
  if(2*l+1!=NShell)
  {
    return luaL_error(L, "Error in LuaOperatorNewJy, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorJy(NFermion, l, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "ssp.operator");
  free(LindexUp);
  free(LindexDn);
  return 1;
}

static int LuaOperatorNewJz(lua_State * L)
{
  OperatorType *op;
  uint32_t Ntemp;
  uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t NShell,l;
  uint16_t *LindexUp, *LindexDn;
  LindexUp = luaL_checkarray_unit16_t(L, 2, &NShell);
  LindexDn = luaL_checkarray_unit16_t(L, 3, &Ntemp);
  if(NShell!=Ntemp)
  {
    return luaL_error(L, "Error in LuaOperatorNewJz, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
  }
  l=(NShell-1)/2;
  if(2*l+1!=NShell)
  {
    return luaL_error(L, "Error in LuaOperatorNewJz, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorJz(NFermion, l, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "ssp.operator");
  free(LindexUp);
  free(LindexDn);
  return 1;
}

static int LuaOperatorNewJplus(lua_State * L)
{
  OperatorType *op;
  uint32_t Ntemp;
  uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t NShell,l;
  uint16_t *LindexUp, *LindexDn;
  LindexUp = luaL_checkarray_unit16_t(L, 2, &NShell);
  LindexDn = luaL_checkarray_unit16_t(L, 3, &Ntemp);
  if(NShell!=Ntemp)
  {
    return luaL_error(L, "Error in LuaOperatorNewJplus, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
  }
  l=(NShell-1)/2;
  if(2*l+1!=NShell)
  {
    return luaL_error(L, "Error in LuaOperatorNewJplus, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorJplus(NFermion, l, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "ssp.operator");
  free(LindexUp);
  free(LindexDn);
  return 1;
}

static int LuaOperatorNewJmin(lua_State * L)
{
  OperatorType *op;
  uint32_t Ntemp;
  uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t NShell,l;
  uint16_t *LindexUp, *LindexDn;
  LindexUp = luaL_checkarray_unit16_t(L, 2, &NShell);
  LindexDn = luaL_checkarray_unit16_t(L, 3, &Ntemp);
  if(NShell!=Ntemp)
  {
    return luaL_error(L, "Error in LuaOperatorNewJmin, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
  }
  l=(NShell-1)/2;
  if(2*l+1!=NShell)
  {
    return luaL_error(L, "Error in LuaOperatorNewJmin, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorJmin(NFermion, l, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "ssp.operator");
  free(LindexUp);
  free(LindexDn);
  return 1;
}

static int LuaOperatorNewJsqr(lua_State * L)
{
  OperatorType *op;
  uint32_t Ntemp;
  uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t NShell,l;
  uint16_t *LindexUp, *LindexDn;
  LindexUp = luaL_checkarray_unit16_t(L, 2, &NShell);
  LindexDn = luaL_checkarray_unit16_t(L, 3, &Ntemp);
  if(NShell!=Ntemp)
  {
    return luaL_error(L, "Error in LuaOperatorNewJsqr, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
  }
  l=(NShell-1)/2;
  if(2*l+1!=NShell)
  {
    return luaL_error(L, "Error in LuaOperatorNewJsqr, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorJsqr(NFermion, l, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "ssp.operator");
  free(LindexUp);
  free(LindexDn);
  return 1;
}

static int LuaOperatorNewldots(lua_State * L)
{
  OperatorType *op;
  uint32_t Ntemp;
  uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t NShell,l;
  uint16_t *LindexUp, *LindexDn;
  LindexUp = luaL_checkarray_unit16_t(L, 2, &NShell);
  LindexDn = luaL_checkarray_unit16_t(L, 3, &Ntemp);
  if(NShell!=Ntemp)
  {
    return luaL_error(L, "Error in LuaOperatorNewldots, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
  }
  l=(NShell-1)/2;
  if(2*l+1!=NShell)
  {
    return luaL_error(L, "Error in LuaOperatorNewldots, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorldots(NFermion, l, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "ssp.operator");
  free(LindexUp);
  free(LindexDn);
  return 1;
}


static int LuaOperatorNewCF(lua_State * L)
{
  OperatorType *op;
  uint32_t i;
  uint32_t len = (uint32_t)lua_gettop(L);
  uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t NShell,Ntemp,lCr,lAn;
  uint16_t *LindexUpCr, *LindexUpAn, *LindexDnCr, *LindexDnAn;
  HarmonicsExpansion *Akm;
  if(len==3)
  {
    LindexUpCr = luaL_checkarray_unit16_t(L, 2, &NShell);
    LindexUpAn = LindexUpCr;
    LindexDnAn = LindexUpCr;
    LindexDnCr = LindexUpCr;
    lCr=(NShell-1)/2;
    lAn=lCr;
    if(2*lCr+1!=NShell)
    {
      return luaL_error(L, "Error in LuaOperatorNewCF, Lenght of Index (%d) is not equal to 2l+1 (%d)\n",NShell,lCr);
    }
  }
  else if(len==4)
  {
    LindexUpCr = luaL_checkarray_unit16_t(L, 2, &NShell);
    LindexUpAn = LindexUpCr;
    LindexDnCr = luaL_checkarray_unit16_t(L, 3, &Ntemp);
    LindexDnAn = LindexDnCr;
    lCr=(NShell-1)/2;
    lAn=lCr;
    if(NShell!=Ntemp)
    {
      return luaL_error(L, "Error in LuaOperatorNewCF, Index of spin up and spin down should relate to the same angular momentum (%d != %d)\n",NShell,Ntemp);
    }
    if(2*lCr+1!=NShell)
    {
      return luaL_error(L, "Error in LuaOperatorNewCF, Lenght of Index (%d) is not equal to 2l+1 (%d)\n",NShell,lCr);
    }
  }
  else if(len==6)
  {
    LindexUpCr = luaL_checkarray_unit16_t(L, 2, &NShell);
    LindexDnCr = luaL_checkarray_unit16_t(L, 3, &Ntemp);
    lCr=(NShell-1)/2;
    if(NShell!=Ntemp)
    {
      return luaL_error(L, "Error in LuaOperatorNewCF, Index of spin up and spin down should relate to the same angular momentum (%d != %d)\n",NShell,Ntemp);
    }
    if(2*lCr+1!=NShell)
    {
      return luaL_error(L, "Error in LuaOperatorNewCF, Lenght of Index (%d) is not equal to 2l+1 (%d)\n",NShell,lCr);
    }
    LindexUpAn = luaL_checkarray_unit16_t(L, 4, &NShell);
    LindexDnAn = luaL_checkarray_unit16_t(L, 5, &Ntemp);
    lAn=(NShell-1)/2;
    if(NShell!=Ntemp)
    {
      return luaL_error(L, "Error in LuaOperatorNewCF, Index of spin up and spin down should relate to the same angular momentum (%d != %d)\n",NShell,Ntemp);
    }
    if(2*lAn+1!=NShell)
    {
      return luaL_error(L, "Error in LuaOperatorNewCF, Lenght of Index (%d) is not equal to 2l+1 (%d)\n",NShell,lAn);
    }
  }
  else
  {
    return luaL_error(L, "Error in LuaOperatorNewCF, called with %d arguments, 3,4 or 6 expected\n",len);
  }
  Akm=luaL_checkarray_HarmonicsExpansion(L, -1);
  if(len==3)
  {
    for(i=0;i<Akm->N;i++)
    {
      Akm->Re[i]*=0.5;
      if(Akm->QComplex) Akm->Im[i]*=0.5;
    }
  }
  
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorCF(NFermion, lCr, lAn, LindexUpCr, LindexDnCr, LindexUpAn, LindexDnAn, *Akm, op);
  luaL_setmetatable(L, "ssp.operator");
  
  free(Akm->l);
  free(Akm->m);
  free(Akm->Re);
  if(Akm->QComplex)  free(Akm->Im);
  free(Akm);
  
  if(len==3)
  {
    free(LindexUpCr);
  }
  else if(len==4)
  {
    free(LindexUpCr);
    free(LindexDnCr);
  }
  else if(len==6)
  {
    free(LindexUpCr);
    free(LindexDnCr);
    free(LindexUpAn);
    free(LindexDnAn);
  }
  return 1;
}

static int LuaOperatorNewCreate(lua_State * L)
{
  int i;
  OperatorType *op;
  uint32_t len = (uint32_t)lua_gettop(L);
  if(len==2)
  {
    uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
    uint16_t Crindex = (uint16_t)luaL_checkunsigned(L, 2);
    op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    CreateOperatorCreate(NFermion, Crindex, op);
    luaL_setmetatable(L, "ssp.operator");
  }
  else if(len==3)
  {
    uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
    uint32_t NOpp,NNum;
    uint16_t *CrIndex = luaL_checkarray_unit16_t(L, 2, &NOpp);
    Complex16 *Prefactors = luaL_checkarray_complex16(L, 3, &NNum);
    if(NOpp!=NNum)
    {
      return luaL_error(L, "Error in LuaOperatorNewCreate, Lenght of Index (%d) is not equal to lenght of prefactors (%d)\n",NOpp,NNum);
    }
    double PrefactorsRe[NOpp],PrefactorsIm[NOpp];
    for(i=0;i<NOpp;i++)
    {
      PrefactorsRe[i]=Prefactors[i].Re;
      PrefactorsIm[i]=Prefactors[i].Im;
    }
    op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    CreateOperatorCreateList(NFermion, CrIndex, PrefactorsRe, PrefactorsIm, NOpp, op);
    luaL_setmetatable(L, "ssp.operator");
    free(CrIndex);
    free(Prefactors);
  }
  else
  {
    return luaL_error(L, "Error in LuaOperatorNewCreate, expected 2 or 3 arguments, not %d\n",len);
  }
  return 1;
}

static int LuaOperatorNewAnnihilate(lua_State * L)
{
  int i;
  OperatorType *op;
  uint32_t len = (uint32_t)lua_gettop(L);
  if(len==2)
  {
    uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
    uint16_t Crindex = (uint16_t)luaL_checkunsigned(L, 2);
    op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    CreateOperatorAnnihilate(NFermion, Crindex, op);
    luaL_setmetatable(L, "ssp.operator");
  }
  else if(len==3)
  {
    uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
    uint32_t NOpp,NNum;
    uint16_t *CrIndex = luaL_checkarray_unit16_t(L, 2, &NOpp);
    Complex16 *Prefactors = luaL_checkarray_complex16(L, 3, &NNum);
    if(NOpp!=NNum)
    {
      return luaL_error(L, "Error in LuaOperatorNewAnnihilate, Lenght of Index (%d) is not equal to lenght of prefactors (%d)\n",NOpp,NNum);
    }
    double PrefactorsRe[NOpp],PrefactorsIm[NOpp];
    for(i=0;i<NOpp;i++)
    {
      PrefactorsRe[i]=Prefactors[i].Re;
      PrefactorsIm[i]=Prefactors[i].Im;
    }
    op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    CreateOperatorAnnihilateList(NFermion, CrIndex, PrefactorsRe, PrefactorsIm, NOpp, op);
    luaL_setmetatable(L, "ssp.operator");
    free(CrIndex);
    free(Prefactors);
  }
  else
  {
    return luaL_error(L, "Error in LuaOperatorNewAnnihilate, expected 2 or 3 arguments, not %d\n",len);
  }
  return 1;
}


static int LuaOperatorNewNumber(lua_State * L)
{
  int i;
  OperatorType *op;
  uint32_t len = (uint32_t)lua_gettop(L);
  if(len==3)
  {
    uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
    uint16_t Crindex = (uint16_t)luaL_checkunsigned(L, 2);
    uint16_t Anindex = (uint16_t)luaL_checkunsigned(L, 3);
    op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    double one=1;
    double zero=0;
    CreateOperatorCreateAnnihilateList(NFermion, &Crindex, &Anindex, &one, &zero, 1, op);
    luaL_setmetatable(L, "ssp.operator");
  }
  else if(len==4)
  {
    uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
    uint32_t NOppCr,NOppAn,NNum;
    uint16_t *CrIndex = luaL_checkarray_unit16_t(L, 2, &NOppCr);
    uint16_t *AnIndex = luaL_checkarray_unit16_t(L, 3, &NOppAn);
    Complex16 *Prefactors = luaL_checkarray_complex16(L, 4, &NNum);
    if(NOppCr!=NOppAn)
    {
      return luaL_error(L, "Error in LuaOperatorNewNumber, Lenght of Index for Createion (%d) is not equal to lenght of index for Annihilation (%d)\n",NOppCr,NOppAn);
    }
    if(NOppCr!=NNum)
    {
      return luaL_error(L, "Error in LuaOperatorNewNumber, Lenght of Index (%d) is not equal to lenght of prefactors (%d)\n",NOppCr,NNum);
    }
    double PrefactorsRe[NOppCr],PrefactorsIm[NOppCr];
    for(i=0;i<NOppCr;i++)
    {
      PrefactorsRe[i]=Prefactors[i].Re;
      PrefactorsIm[i]=Prefactors[i].Im;
    }
    op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    CreateOperatorCreateAnnihilateList(NFermion, CrIndex, AnIndex, PrefactorsRe, PrefactorsIm, NOppCr, op);
    luaL_setmetatable(L, "ssp.operator");
    free(CrIndex);
    free(AnIndex);
    free(Prefactors);
  }
  else
  {
    return luaL_error(L, "Error in LuaOperatorNewAnnihilate, expected 3 or 4 arguments, not %d\n",len);
  }
  return 1;
}



/***********************/
/* Operator meta actions */
/***********************/

static int LuaOperatorAdd(lua_State * L)
{
  OperatorType *l = (OperatorType *) luaL_checkudata(L, 1, "ssp.operator");
  OperatorType *r = (OperatorType *) luaL_checkudata(L, 2, "ssp.operator");
  lua_settop(L, 0);
  OperatorType *res = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  luaL_setmetatable(L, "ssp.operator");
  OperatorCopy(res, *l);
  OperatorAdd(res, *r, 1., 0);
  return 1;
}

static int LuaOperatorSub(lua_State * L)
{
  OperatorType *l = (OperatorType *) luaL_checkudata(L, 1, "ssp.operator");
  OperatorType *r = (OperatorType *) luaL_checkudata(L, 2, "ssp.operator");
  OperatorType *res = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  luaL_setmetatable(L, "ssp.operator");
  OperatorCopy(res, *l);
  OperatorAdd(res, *r, -1., 0);
  return 1;
}


/*  multiplies two objects, whereby one is an operator
 *  a (complex) number times operator scales the operator - multiplication commutes
 *  a operator times operator is implemented as: A * B - multiplication does not commute
 *  a operator times wavefunction is implemented as O * psi and returns a wavefunction
 *  a wavefunction times operator is implemented as psi * O == ConjugateTranspose(O) * psi
 *  a table times operator is implemented as: {a,b,c} * O == {a * O, b * O, c * O}
 *  a operator times table is implemented as: O * {a,b,c} == {O * a, O * b, O * c}
 */
static int LuaOperatorMul(lua_State * L)
{
  uint32_t len = (uint32_t)lua_gettop(L);
  if(lua_type(L,len)==LUA_TTABLE)
  {
    /* psi * table */
    lua_newtable(L);
    lua_pushnil(L);                /* push a nil to the stack - needed for lua_next */
    while(lua_next(L, len) != 0)   /* push the next index and value to the end of the stack +1 on stack (last element replaced)*/
    {
      LuaOperatorMul(L);           /* multimply the value with the operator (first times last of the stack) +1 on stack*/
      lua_pushvalue(L, -3);        /* push the index to the end of the stack +1 on stack*/
      lua_insert(L, -2);           /* swap index and value */
      lua_settable(L, -5);         /* set new table value -2 on stack */
      lua_pop(L, 1);               /* remove one element from the stack -1 on stack */
    }
    return 1;
  }
  if(lua_type(L,1)==LUA_TTABLE)
  {
    /* psi * table */
    lua_newtable(L);
    lua_pushnil(L);                /* push a nil to the stack - needed for lua_next */
    while(lua_next(L, 1) != 0)     /* push the next index and value to the end of the stack +1 on stack (last element replaced)*/
    {
      lua_insert(L, 1);            /* set value from table to first element on stack */
      lua_pushvalue(L, len+1);     /* push the operator to the end of the stack +1 on stack*/
      LuaOperatorMul(L);           /* multimply the value with the operator (first times last of the stack) +1 on stack*/
      lua_pushvalue(L, -3);        /* push the index to the end of the stack +1 on stack*/
      lua_insert(L, -2);           /* swap index and value */
      lua_settable(L, -5);         /* set new table value -2 on stack */
      lua_pop(L, 1);               /* remove one element from the stack -1 on stack */
      lua_remove(L, 1);            /* remove first element from the stack -1 on stack */
    }
    return 1;
  }
  if (lua_isnumber(L, 1))
  {
    double b = (double) luaL_checknumber(L, 1);
    OperatorType *O = (OperatorType*) luaL_checkudata(L, len, "ssp.operator");
    OperatorType *ONew;
    ONew = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    OperatorCopy(ONew, *O);
    OperatorScale(ONew, b, 0);
    luaL_setmetatable(L, "ssp.operator");
    return 1;
  }
  if (lua_isnumber(L, len))
  {
    double b = (double) luaL_checknumber(L, len);
    OperatorType *O = (OperatorType*) luaL_checkudata(L, 1, "ssp.operator");
    OperatorType *ONew;
    ONew = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    OperatorCopy(ONew, *O);
    OperatorScale(ONew, b, 0);
    luaL_setmetatable(L, "ssp.operator");
    return 1;
  }
  if (luaL_testudata(L, 1, "Lua.complex") != NULL)
  {
    Complex16 *b = (Complex16 *) luaL_checkudata(L, 1, "Lua.complex");
    OperatorType *O = (OperatorType*) luaL_checkudata(L, len, "ssp.operator");
    OperatorType *ONew;
    ONew = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    OperatorCopy(ONew, *O);
    OperatorScale(ONew, b->Re, b->Im);
    luaL_setmetatable(L, "ssp.operator");
    return 1;
  }
  if (luaL_testudata(L, len, "Lua.complex") != NULL)
  {
    Complex16 *b = (Complex16 *) luaL_checkudata(L, len, "Lua.complex");
    OperatorType *O = (OperatorType*) luaL_checkudata(L, 1, "ssp.operator");
    OperatorType *ONew;
    ONew = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    OperatorCopy(ONew, *O);
    OperatorScale(ONew, b->Re, b->Im);
    luaL_setmetatable(L, "ssp.operator");
    return 1;
  }
  if (luaL_testudata(L, len, "ssp.wavefunction") != NULL)
  {
    /* O | psi> */
    WaveFunctionType *psi = (WaveFunctionType*) luaL_checkudata(L, len, "ssp.wavefunction");
    OperatorType *O = (OperatorType *) luaL_checkudata(L,1,"ssp.operator");
    WaveFunctionType *Opsi;
    Opsi = (WaveFunctionType *) lua_newuserdata(L, sizeof(WaveFunctionType));
    OperatorPsi(*O, *psi, Opsi);
    luaL_setmetatable(L, "ssp.wavefunction");
    return 1;
  }
  if (luaL_testudata(L, 1, "ssp.wavefunction") != NULL)
  {
    /* <psi | O */
    WaveFunctionType *psi = (WaveFunctionType*) luaL_checkudata(L, 1, "ssp.wavefunction");
    OperatorType *O = (OperatorType *) luaL_checkudata(L,len,"ssp.operator");
    OperatorType OT;
    OperatorConjugateTransposeNormalOrder(*O, &OT, &QDetNotOrderedFcFaBcBa);
    WaveFunctionType *Opsi;
    Opsi = (WaveFunctionType *) lua_newuserdata(L, sizeof(WaveFunctionType));
    OperatorPsi(OT, *psi, Opsi);
    FreeOperator(&OT);
    luaL_setmetatable(L, "ssp.wavefunction");
    return 1;
  }
  {
    /* O * O */
    OperatorType *l = (OperatorType *) luaL_checkudata(L, 1, "ssp.operator");
    OperatorType *r = (OperatorType *) luaL_checkudata(L, len, "ssp.operator");
    OperatorType *res = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    luaL_setmetatable(L, "ssp.operator");
    InitOperatorForMult(res, 2, l, r);
    OperatorMultiply2NormalOrder(res, *l, *r, 'N', 'N', &QDetNotOrderedFcFaBcBa);
  }
  return 1;
}



static int LuaOperatorPrint(lua_State * L)
{
  uint32_t len = (uint32_t)lua_gettop(L);
  if(len==0)
  {
    PrintOperatorReadable(*((OperatorType *) currentudata));
  }
  else if(lua_istable(L,-1))
  {
    printf("Options for print recognized but not implemented\n");
    PrintOperatorReadable(*((OperatorType *) currentudata));
  }
  else
  {
    return luaL_error(L, "Error in .Print(arg) unrecognized arguments in Print()\n");
  }
  return 0;
}

static int LuaOperatorChop(lua_State * L)
{
  uint32_t len = (uint32_t)lua_gettop(L);
  double border;
  switch (len)
  {
    case 0:
      border=10*DBL_EPSILON;
    break;
    case 1:
      border=luaL_checknumber(L, 1);
    break;
    default:
      return luaL_error(L, "Errorr, call to chop with %d arguments, expected 0 or 1\n",len);
    break;
  }
  OperatorChop((OperatorType *) currentudata, border);
  return 0;
}


static int LuaOperatorConjugate(lua_State * L)
{
  uint32_t len = (uint32_t)lua_gettop(L);
  if(len==0)
  {
    OperatorComplexConjugate((OperatorType *)currentudata);
  }
  else
  {
    return luaL_error(L, "Error in .Conjugate(arg) unrecognized arguments in Conjugate()\n");
  }
  return 0;
}

static int LuaOperatorConjugateTranspose(lua_State * L)
{
  uint32_t len = (uint32_t)lua_gettop(L);
  if(len==0)
  {
    OperatorType OppTemp;
    OperatorCopy(&OppTemp, *(OperatorType *)currentudata);
    FreeOperator((OperatorType *)currentudata);
    OperatorConjugateTransposeNormalOrder(OppTemp, (OperatorType *)currentudata, &QDetNotOrderedFcFaBaBc);
    FreeOperator(&OppTemp);
  }
  else
  {
    return luaL_error(L, "Error in .ConjugateTranspose(arg) unrecognized arguments in ConjugateTranspose()\n");
  }
  return 0;
}

static int LuaOperatorTranspose(lua_State * L)
{
  uint32_t len = (uint32_t)lua_gettop(L);
  if(len==0)
  {
    OperatorType OppTemp;
    OperatorCopy(&OppTemp, *(OperatorType *)currentudata);
    FreeOperator((OperatorType *)currentudata);
    OperatorTransposeNormalOrder(OppTemp, (OperatorType *)currentudata, &QDetNotOrderedFcFaBaBc);
    FreeOperator(&OppTemp);
  }
  else
  {
    return luaL_error(L, "Error in .ConjugateTranspose(arg) unrecognized arguments in ConjugateTranspose()\n");
  }
  return 0;
}

static int LuaOperatorIndex(lua_State * L)
{
  OperatorType *a = (OperatorType *) luaL_checkudata(L, 1, "ssp.operator");
  const char *str = lua_tostring(L, 2);
  currentudata = a;
  
  if (strcmp(str, "Print") == 0)
  {
    lua_pushcfunction(L, LuaOperatorPrint);
    return 1;
  }
  else if (strcmp(str, "Chop") == 0)
  {
    lua_pushcfunction(L, LuaOperatorChop);
    return 1;
  }
  if ((strcmp(str, "Conjugate") == 0)||(strcmp(str, "C") == 0)||(strcmp(str, "Dag") == 0))
  {
    lua_pushcfunction(L, LuaOperatorConjugate);
    return 1;
  }
  if ((strcmp(str, "ConjugateTranspose") == 0)||(strcmp(str, "CT")==0))
  {
    lua_pushcfunction(L, LuaOperatorConjugateTranspose);
    return 1;
  }
  if ((strcmp(str, "Transpose") == 0)||(strcmp(str, "T")==0))
  {
    lua_pushcfunction(L, LuaOperatorTranspose);
    return 1;
  }
  return luaL_error(L, "Error operator variable not found");;
}

static int LuaOperatorToString(lua_State * L)
{
  OperatorType *res = (OperatorType *) luaL_checkudata(L, 1, "ssp.operator");
  lua_pushstring(L, res->Name);
  return 1;
}

static int LuaOperatorDestroy(lua_State * L)
{
  OperatorType *a = (OperatorType *) luaL_checkudata(L, 1, "ssp.operator");
  FreeOperator(a);
  return 0;
}

static const struct luaL_Reg Operator_meta[] = {
  {"__add", LuaOperatorAdd},
  {"__sub", LuaOperatorSub},
  {"__mul", LuaOperatorMul},
  {"__index", LuaOperatorIndex},
  {"__tostring", LuaOperatorToString},
  {"__gc", LuaOperatorDestroy},
  {NULL, NULL}
};