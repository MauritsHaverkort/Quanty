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
      lua_pushfstring (L, "Error in LuaOperatorNewCoulomb, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
      lua_error(L);
      return 1;
    }
    l=(NShell-1)/2;
    if(2*l+1!=NShell)
    {
      lua_pushfstring (L, "Error in LuaOperatorNewCoulomb, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
      lua_error(L);
      return 1;
    }
    F = luaL_checkarray_double(L, 4, &Ntemp);
    if(Ntemp!=l+1)
    {
      lua_pushfstring (L, "Error in LuaOperatorNewCoulomb, Lenght of F (%d) is not equal to l+1\n",Ntemp);
      lua_error(L);
      return 1;
    }
    op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    CreateOperatorCoulombSingleShell(NFermion, l, LindexUp, LindexDn, F, op);
    luaL_setmetatable(L, "Lua.operator");
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
      lua_pushfstring (L, "Error in LuaOperatorNewCoulomb, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
      lua_error(L);
      return 1;
    }
    lA=(NShell-1)/2;
    if(2*lA+1!=NShell)
    {
      lua_pushfstring (L, "Error in LuaOperatorNewCoulomb, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
      lua_error(L);
      return 1;
    }
    LindexUpB = luaL_checkarray_unit16_t(L, 4, &NShell);
    LindexDnB = luaL_checkarray_unit16_t(L, 5, &Ntemp);
    if(NShell!=Ntemp)
    {
      lua_pushfstring (L, "Error in LuaOperatorNewCoulomb, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
      lua_error(L);
      return 1;
    }
    lB=(NShell-1)/2;
    if(2*lB+1!=NShell)
    {
      lua_pushfstring (L, "Error in LuaOperatorNewCoulomb, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
      lua_error(L);
      return 1;
    }
    F = luaL_checkarray_double(L, 6, &Ntemp);
    if(Ntemp!=(Min(lA,lB)+1))
    {
      lua_pushfstring (L, "Error in LuaOperatorNewCoulomb, Lenght of F (%d) is not equal to (Min(lA,lB)+1)1\n",Ntemp);
      lua_error(L);
      return 1;
    }
    G = luaL_checkarray_double(L, 7, &Ntemp);
    if(Ntemp!=(Min(lA,lB)+1))
    {
      lua_pushfstring (L, "Error in LuaOperatorNewCoulomb, Lenght of G (%d) is not equal to (Min(lA,lB)+1)1\n",Ntemp);
      lua_error(L);
      return 1;
    }
    op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    CreateOperatorCoulombTwoShellsConstantOccupation(NFermion, lA, lB, LindexUpA, LindexDnA, LindexUpB, LindexDnB, F, G, op);
    luaL_setmetatable(L, "Lua.operator");
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
    lua_pushfstring (L, "Error in LuaOperatorNewCoulomb, expected 4 or 7 arguments not %d",len);
    lua_error(L);
    return 1;
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
    lua_pushfstring (L, "Error in LuaOperatorNewSx, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
    lua_error(L);
    return 1;
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorSx(NFermion, NShell, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "Lua.operator");
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
    lua_pushfstring (L, "Error in LuaOperatorNewSy, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
    lua_error(L);
    return 1;
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorSy(NFermion, NShell, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "Lua.operator");
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
    lua_pushfstring (L, "Error in LuaOperatorNewSz, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
    lua_error(L);
    return 1;
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorSz(NFermion, NShell, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "Lua.operator");
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
    lua_pushfstring (L, "Error in LuaOperatorNewSsqr, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
    lua_error(L);
    return 1;
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorSsqr(NFermion, NShell, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "Lua.operator");
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
    lua_pushfstring (L, "Error in LuaOperatorNewSplus, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
    lua_error(L);
    return 1;
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorSplus(NFermion, NShell, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "Lua.operator");
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
    lua_pushfstring (L, "Error in LuaOperatorNewSmin, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
    lua_error(L);
    return 1;
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorSmin(NFermion, NShell, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "Lua.operator");
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
    lua_pushfstring (L, "Error in LuaOperatorNewLx, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
    lua_error(L);
    return 1;
  }
  l=(NShell-1)/2;
  if(2*l+1!=NShell)
  {
    lua_pushfstring (L, "Error in LuaOperatorNewLx, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
    lua_error(L);
    return 1;
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorLx(NFermion, l, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "Lua.operator");
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
    lua_pushfstring (L, "Error in LuaOperatorNewLy, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
    lua_error(L);
    return 1;
  }
  l=(NShell-1)/2;
  if(2*l+1!=NShell)
  {
    lua_pushfstring (L, "Error in LuaOperatorNewLy, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
    lua_error(L);
    return 1;
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorLy(NFermion, l, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "Lua.operator");
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
    lua_pushfstring (L, "Error in LuaOperatorNewLz, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
    lua_error(L);
    return 1;
  }
  l=(NShell-1)/2;
  if(2*l+1!=NShell)
  {
    lua_pushfstring (L, "Error in LuaOperatorNewLz, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
    lua_error(L);
    return 1;
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorLz(NFermion, l, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "Lua.operator");
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
    lua_pushfstring (L, "Error in LuaOperatorNewLplus, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
    lua_error(L);
    return 1;
  }
  l=(NShell-1)/2;
  if(2*l+1!=NShell)
  {
    lua_pushfstring (L, "Error in LuaOperatorNewLplus, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
    lua_error(L);
    return 1;
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorLplus(NFermion, l, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "Lua.operator");
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
    lua_pushfstring (L, "Error in LuaOperatorNewLmin, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
    lua_error(L);
    return 1;
  }
  l=(NShell-1)/2;
  if(2*l+1!=NShell)
  {
    lua_pushfstring (L, "Error in LuaOperatorNewLmin, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
    lua_error(L);
    return 1;
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorLmin(NFermion, l, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "Lua.operator");
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
    lua_pushfstring (L, "Error in LuaOperatorNewLsqr, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
    lua_error(L);
    return 1;
  }
  l=(NShell-1)/2;
  if(2*l+1!=NShell)
  {
    lua_pushfstring (L, "Error in LuaOperatorNewLsqr, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
    lua_error(L);
    return 1;
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorLsqr(NFermion, l, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "Lua.operator");
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
    lua_pushfstring (L, "Error in LuaOperatorNewJx, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
    lua_error(L);
    return 1;
  }
  l=(NShell-1)/2;
  if(2*l+1!=NShell)
  {
    lua_pushfstring (L, "Error in LuaOperatorNewJx, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
    lua_error(L);
    return 1;
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorJx(NFermion, l, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "Lua.operator");
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
    lua_pushfstring (L, "Error in LuaOperatorNewJy, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
    lua_error(L);
    return 1;
  }
  l=(NShell-1)/2;
  if(2*l+1!=NShell)
  {
    lua_pushfstring (L, "Error in LuaOperatorNewJy, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
    lua_error(L);
    return 1;
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorJy(NFermion, l, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "Lua.operator");
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
    lua_pushfstring (L, "Error in LuaOperatorNewJz, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
    lua_error(L);
    return 1;
  }
  l=(NShell-1)/2;
  if(2*l+1!=NShell)
  {
    lua_pushfstring (L, "Error in LuaOperatorNewJz, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
    lua_error(L);
    return 1;
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorJz(NFermion, l, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "Lua.operator");
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
    lua_pushfstring (L, "Error in LuaOperatorNewJplus, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
    lua_error(L);
    return 1;
  }
  l=(NShell-1)/2;
  if(2*l+1!=NShell)
  {
    lua_pushfstring (L, "Error in LuaOperatorNewJplus, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
    lua_error(L);
    return 1;
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorJplus(NFermion, l, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "Lua.operator");
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
    lua_pushfstring (L, "Error in LuaOperatorNewJmin, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
    lua_error(L);
    return 1;
  }
  l=(NShell-1)/2;
  if(2*l+1!=NShell)
  {
    lua_pushfstring (L, "Error in LuaOperatorNewJmin, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
    lua_error(L);
    return 1;
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorJmin(NFermion, l, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "Lua.operator");
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
    lua_pushfstring (L, "Error in LuaOperatorNewJsqr, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
    lua_error(L);
    return 1;
  }
  l=(NShell-1)/2;
  if(2*l+1!=NShell)
  {
    lua_pushfstring (L, "Error in LuaOperatorNewJsqr, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
    lua_error(L);
    return 0;
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorJsqr(NFermion, l, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "Lua.operator");
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
    lua_pushfstring (L, "Error in LuaOperatorNewldots, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
    lua_error(L);
    return 1;
  }
  l=(NShell-1)/2;
  if(2*l+1!=NShell)
  {
    lua_pushfstring (L, "Error in LuaOperatorNewldots, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
    lua_error(L);
    return 0;
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorldots(NFermion, l, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "Lua.operator");
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
      lua_pushfstring (L, "Error in LuaOperatorNewCF, Lenght of Index (%d) is not equal to 2l+1 (%d)\n",NShell,lCr);
      lua_error(L);
      return 1;
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
      lua_pushfstring (L, "Error in LuaOperatorNewCF, Index of spin up and spin down should relate to the same angular momentum (%d != %d)\n",NShell,Ntemp);
      lua_error(L);
      return 1;
    }
    if(2*lCr+1!=NShell)
    {
      lua_pushfstring (L, "Error in LuaOperatorNewCF, Lenght of Index (%d) is not equal to 2l+1 (%d)\n",NShell,lCr);
      lua_error(L);
      return 1;
    }
  }
  else if(len==6)
  {
    LindexUpCr = luaL_checkarray_unit16_t(L, 2, &NShell);
    LindexDnCr = luaL_checkarray_unit16_t(L, 3, &Ntemp);
    lCr=(NShell-1)/2;
    if(NShell!=Ntemp)
    {
      lua_pushfstring (L, "Error in LuaOperatorNewCF, Index of spin up and spin down should relate to the same angular momentum (%d != %d)\n",NShell,Ntemp);
      lua_error(L);
      return 1;
    }
    if(2*lCr+1!=NShell)
    {
      lua_pushfstring (L, "Error in LuaOperatorNewCF, Lenght of Index (%d) is not equal to 2l+1 (%d)\n",NShell,lCr);
      lua_error(L);
      return 1;
    }
    LindexUpAn = luaL_checkarray_unit16_t(L, 4, &NShell);
    LindexDnAn = luaL_checkarray_unit16_t(L, 5, &Ntemp);
    lAn=(NShell-1)/2;
    if(NShell!=Ntemp)
    {
      lua_pushfstring (L, "Error in LuaOperatorNewCF, Index of spin up and spin down should relate to the same angular momentum (%d != %d)\n",NShell,Ntemp);
      lua_error(L);
      return 1;
    }
    if(2*lAn+1!=NShell)
    {
      lua_pushfstring (L, "Error in LuaOperatorNewCF, Lenght of Index (%d) is not equal to 2l+1 (%d)\n",NShell,lAn);
      lua_error(L);
      return 1;
    }
  }
  else
  {
    lua_pushfstring (L, "Error in LuaOperatorNewCF, called with %d arguments, 3,4 or 6 expected\n",len);
    lua_error(L);
    return 1;
  }
  Akm=luaL_checkarray_HarmonicsExpansion(L, -1, &i);
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
  luaL_setmetatable(L, "Lua.operator");
  
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
  OperatorType *op;
  uint32_t len = (uint32_t)lua_gettop(L);
  if(len==2)
  {
    uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
    uint16_t Crindex = (uint16_t)luaL_checkunsigned(L, 2);
    op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    CreateOperatorCreate(NFermion, Crindex, op);
    luaL_setmetatable(L, "Lua.operator");
  }
  else if(len==3)
  {
    uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
    uint32_t NOpp,NNum;
    uint16_t *CrIndex = luaL_checkarray_unit16_t(L, 2, &NOpp);
    Complex16 *Prefactors = luaL_checkarray_complex16(L, 3, &NNum);
    if(NOpp!=NNum)
    {
      lua_pushfstring (L, "Error in LuaOperatorNewCreate, Lenght of Index (%d) is not equal to lenght of prefactors (%d)\n",NOpp,NNum);
      lua_error(L);
      return 1;
    }
    double PrefactorsRe[NOpp],PrefactorsIm[NOpp];
    for(int i=0;i<NOpp;i++)
    {
      PrefactorsRe[i]=Prefactors[i].Re;
      PrefactorsIm[i]=Prefactors[i].Im;
    }
    op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    CreateOperatorCreateList(NFermion, CrIndex, PrefactorsRe, PrefactorsIm, NOpp, op);
    luaL_setmetatable(L, "Lua.operator");
    free(CrIndex);
    free(Prefactors);
  }
  else
  {
    lua_pushfstring (L, "Error in LuaOperatorNewCreate, expected 2 or 3 arguments, not %d\n",len);
    lua_error(L);
    return 0;
  }
  return 1;
}

static int LuaOperatorNewAnnihilate(lua_State * L)
{
  OperatorType *op;
  uint32_t len = (uint32_t)lua_gettop(L);
  if(len==2)
  {
    uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
    uint16_t Crindex = (uint16_t)luaL_checkunsigned(L, 2);
    op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    CreateOperatorAnnihilate(NFermion, Crindex, op);
    luaL_setmetatable(L, "Lua.operator");
  }
  else if(len==3)
  {
    uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
    uint32_t NOpp,NNum;
    uint16_t *CrIndex = luaL_checkarray_unit16_t(L, 2, &NOpp);
    Complex16 *Prefactors = luaL_checkarray_complex16(L, 3, &NNum);
    if(NOpp!=NNum)
    {
      lua_pushfstring (L, "Error in LuaOperatorNewAnnihilate, Lenght of Index (%d) is not equal to lenght of prefactors (%d)\n",NOpp,NNum);
      lua_error(L);
      return 1;
    }
    double PrefactorsRe[NOpp],PrefactorsIm[NOpp];
    for(int i=0;i<NOpp;i++)
    {
      PrefactorsRe[i]=Prefactors[i].Re;
      PrefactorsIm[i]=Prefactors[i].Im;
    }
    op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    CreateOperatorAnnihilateList(NFermion, CrIndex, PrefactorsRe, PrefactorsIm, NOpp, op);
    luaL_setmetatable(L, "Lua.operator");
    free(CrIndex);
    free(Prefactors);
  }
  else
  {
    lua_pushfstring (L, "Error in LuaOperatorNewAnnihilate, expected 2 or 3 arguments, not %d\n",len);
    lua_error(L);
    return 0;
  }
  return 1;
}


static int LuaOperatorNewNumber(lua_State * L)
{
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
    luaL_setmetatable(L, "Lua.operator");
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
      lua_pushfstring (L, "Error in LuaOperatorNewNumber, Lenght of Index for Createion (%d) is not equal to lenght of index for Annihilation (%d)\n",NOppCr,NOppAn);
      lua_error(L);
      return 1;
    }
    if(NOppCr!=NNum)
    {
      lua_pushfstring (L, "Error in LuaOperatorNewNumber, Lenght of Index (%d) is not equal to lenght of prefactors (%d)\n",NOppCr,NNum);
      lua_error(L);
      return 1;
    }
    double PrefactorsRe[NOppCr],PrefactorsIm[NOppCr];
    for(int i=0;i<NOppCr;i++)
    {
      PrefactorsRe[i]=Prefactors[i].Re;
      PrefactorsIm[i]=Prefactors[i].Im;
    }
    op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    CreateOperatorCreateAnnihilateList(NFermion, CrIndex, AnIndex, PrefactorsRe, PrefactorsIm, NOppCr, op);
    luaL_setmetatable(L, "Lua.operator");
    free(CrIndex);
    free(AnIndex);
    free(Prefactors);
  }
  else
  {
    lua_pushfstring (L, "Error in LuaOperatorNewAnnihilate, expected 3 or 4 arguments, not %d\n",len);
    lua_error(L);
    return 1;
  }
  return 1;
}



/***********************/
/* operator arithmetic */
/***********************/

static int LuaOperatorMul(lua_State * L)
{
  if (luaL_testudata(L, 1, "Lua.operator") != NULL)
  {
    OperatorType *l = (OperatorType *) luaL_checkudata(L, 1, "Lua.operator");
    if (lua_isnumber(L, 2))
    {
      OperatorType *res = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
      luaL_setmetatable(L, "Lua.operator");
      double r = luaL_checknumber(L, 2);
      OperatorCopy(res, *l);
      OperatorScale(res, r, 0);
    }
    else if (luaL_testudata(L, 2, "Lua.complex") != NULL)
    {
      OperatorType *res = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
      luaL_setmetatable(L, "Lua.operator");
      Complex16 *r = (Complex16 *)luaL_checkudata(L, 2, "Lua.complex");
      OperatorCopy(res, *l);
      OperatorScale(res, r->Re, r->Im);
    }
    else if (luaL_testudata(L, 2, "Lua.wavefunction") != NULL)
    {
      WaveFunctionType *res = (WaveFunctionType *) lua_newuserdata(L, sizeof(WaveFunctionType));
      luaL_setmetatable(L, "Lua.wavefunction");
      WaveFunctionType *r = (WaveFunctionType *)luaL_checkudata(L,2,"Lua.wavefunction");
      OperatorPsi(*l, *r, res);
    }
    else
    {
      OperatorType *res = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
      luaL_setmetatable(L, "Lua.operator");
      OperatorType *r = (OperatorType *) luaL_checkudata(L, 2, "Lua.operator");
      InitOperatorForMult(res, 2, l, r);
      OperatorMultiply2NormalOrder(res, *l, *r, 'N', 'N', &QDetNotOrderedFcFaBcBa);
    }
  }
  else
  {
    OperatorType *l = (OperatorType *) luaL_checkudata(L, 2, "Lua.operator");
    if (lua_isnumber(L, 1))
    {
      OperatorType *res = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
      luaL_setmetatable(L, "Lua.operator");
      double r = luaL_checknumber(L, 1);
      OperatorCopy(res, *l);
      OperatorScale(res, r, 0);
    }
    else
    {
      lua_pushfstring (L, "Error in multiplication of operators probably a c error and not an error in the input script\n");
      lua_error(L);
    }
  }
  return 1;
}

static int LuaOperatorAdd(lua_State * L)
{
  OperatorType *l = (OperatorType *) luaL_checkudata(L, 1, "Lua.operator");
  OperatorType *r = (OperatorType *) luaL_checkudata(L, 2, "Lua.operator");
  lua_settop(L, 0);
  OperatorType *res = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  luaL_setmetatable(L, "Lua.operator");
  OperatorCopy(res, *l);
  OperatorAdd(res, *r, 1., 0);
  return 1;
}

static int LuaOperatorSub(lua_State * L)
{
  OperatorType *l = (OperatorType *) luaL_checkudata(L, 1, "Lua.operator");
  OperatorType *r = (OperatorType *) luaL_checkudata(L, 2, "Lua.operator");
  OperatorType *res = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  luaL_setmetatable(L, "Lua.operator");
  OperatorCopy(res, *l);
  OperatorAdd(res, *r, -1., 0);
  return 1;
}


static int LuaOperatorToString(lua_State * L)
{
  OperatorType *res = (OperatorType *) luaL_checkudata(L, 1, "Lua.operator");
  lua_pushstring(L, res->Name);
  return 1;
}

static int LuaOperatorPrintReadable(lua_State * L)
{
  PrintOperatorReadable(*((OperatorType *) currentudata));
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
      lua_pushfstring (L, "Errorr, call to chop with %d arguments, expected 0 or 1\n",len);
      lua_error(L);
    break;
  }
  OperatorChop((OperatorType *) currentudata, border);
  return 0;
}




/***********************/
/*   lua interfacing   */
/***********************/


static int LuaOperatorIndex(lua_State * L)
{
  OperatorType *a = (OperatorType *) luaL_checkudata(L, 1, "Lua.operator");
  const char *str = lua_tostring(L, 2);
  currentudata = a;
  
  if (strcmp(str, "Print") == 0)
  {
    lua_pushcfunction(L, LuaOperatorPrintReadable);
    return 1;
  }
  else if (strcmp(str, "Chop") == 0)
  {
    lua_pushcfunction(L, LuaOperatorChop);
    return 1;
  }
  return luaL_error(L, "Error operator variable not found");;
}

static int LuaOperatorDestroy(lua_State * L)
{
  OperatorType *a = (OperatorType *) luaL_checkudata(L, 1, "Lua.operator");
  FreeOperator(a);
  return 0;
}


static const struct luaL_Reg Operator_meta[] = {
  {"__index", LuaOperatorIndex},
  {"__add", LuaOperatorAdd},
  {"__sub", LuaOperatorSub},
  {"__mul", LuaOperatorMul},
  {"__tostring", LuaOperatorToString},
  {"__gc", LuaOperatorDestroy},
  {NULL, NULL}
};


static const struct luaL_Reg Operator_func[] = {
  {"NewCoulomb", LuaOperatorNewCoulomb},
  {"NewSx", LuaOperatorNewSx},
  {"NewSy", LuaOperatorNewSy},
  {"NewSz", LuaOperatorNewSz},
  {"NewSsqr", LuaOperatorNewSsqr},
  {"NewSplus", LuaOperatorNewSplus},
  {"NewSmin", LuaOperatorNewSmin},
  {"NewLx", LuaOperatorNewLx},
  {"NewLy", LuaOperatorNewLy},
  {"NewLz", LuaOperatorNewLz},
  {"NewLsqr", LuaOperatorNewLsqr},
  {"NewLplus", LuaOperatorNewLplus},
  {"NewLmin", LuaOperatorNewLmin},
  {"NewJx", LuaOperatorNewSx},
  {"NewJy", LuaOperatorNewSy},
  {"NewJz", LuaOperatorNewSz},
  {"NewJsqr", LuaOperatorNewSsqr},
  {"NewJplus", LuaOperatorNewSplus},
  {"NewJmin", LuaOperatorNewSmin},
  {"Newldots", LuaOperatorNewldots},
  {"NewCF", LuaOperatorNewCF},
  {"NewCreate", LuaOperatorNewCreate},
  {"NewAnnihilate", LuaOperatorNewAnnihilate},
  {"NewNumber", LuaOperatorNewNumber},
  {NULL, NULL}
};

int luaopen_Operator(lua_State * L)
{
  luaL_newmetatable(L, "Lua.operator");
  luaL_setfuncs(L, Operator_meta, 0);
  luaL_newlib(L, Operator_func);
  lua_setglobal(L, "Operator");
  return 1;
}