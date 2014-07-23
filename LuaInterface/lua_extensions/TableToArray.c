/*  returns a pointer to an array of lenght len_out
 *  of type _type
 *  found in the lua stack at position narg
 */

uint8_t* luaL_checkarray_unit8_t(lua_State *L, int narg, uint32_t *len_out)
{
  int i;
  luaL_checktype(L, narg, LUA_TTABLE);
  uint32_t len = (uint32_t)lua_objlen(L, narg);
  *len_out = len;
  uint8_t *buff = (uint8_t*)malloc(len*sizeof(uint8_t));
  for(i=0;i<len;i++)
  {
    lua_rawgeti(L, narg, i+1);
    buff[i]=(uint8_t)luaL_checkunsigned(L, -1);
    lua_pop(L, 1);
  }
  return buff;
}

uint16_t* luaL_checkarray_unit16_t(lua_State *L, int narg, uint32_t *len_out)
{
  int i;
  luaL_checktype(L, narg, LUA_TTABLE);
  uint32_t len = (uint32_t)lua_objlen(L, narg);
  *len_out = len;
  uint16_t *buff = (uint16_t*)malloc(len*sizeof(uint16_t));
  for(i=0;i<len;i++)
  {
    lua_rawgeti(L, narg, i+1);
    buff[i]=(uint16_t)luaL_checkunsigned(L, -1);
    lua_pop(L, 1);
  }
  return buff;
}

uint32_t* luaL_checkarray_unit32_t(lua_State *L, int narg, uint32_t *len_out)
{
  int i;
  luaL_checktype(L, narg, LUA_TTABLE);
  uint32_t len = (uint32_t)lua_objlen(L, narg);
  *len_out = len;
  uint32_t *buff = (uint32_t*)malloc(len*sizeof(uint32_t));
  for(i=0;i<len;i++)
  {
    lua_rawgeti(L, narg, i+1);
    buff[i]=(uint32_t)luaL_checkunsigned(L, -1);
    lua_pop(L, 1);
  }
  return buff;
}

double* luaL_checkarray_double(lua_State *L, int narg, uint32_t *len_out)
{
  int i;
  luaL_checktype(L, narg, LUA_TTABLE);
  uint32_t len = (uint32_t)lua_objlen(L, narg);
  *len_out = len;
  double *buff = (double*)malloc(len*sizeof(double));
  for(i=0;i<len;i++)
  {
    lua_rawgeti(L, narg, i+1);
    buff[i]=(double)luaL_checknumber(L, -1);
    lua_pop(L, 1);
  }
  return buff;
}


Complex16* luaL_checkarray_complex16(lua_State *L, int narg, uint32_t *len_out)
{
  int i;
  luaL_checktype(L, narg, LUA_TTABLE);
  uint32_t len = (uint32_t)lua_objlen(L, narg);
  *len_out = len;
  Complex16 *buff = (Complex16*)malloc(len*sizeof(Complex16));
  for(i=0;i<len;i++)
  {
    lua_rawgeti(L, narg, i+1);
    if(lua_isnumber(L, -1))
    {
      buff[i].Re=(double)luaL_checknumber(L, -1);
      lua_pop(L, 1);
      buff[i].Im=0;
    }
    else if (luaL_testudata(L, -1, "Lua.complex") != NULL)
    {
      buff[i] = *(Complex16 *)luaL_checkudata(L, -1, "Lua.complex");
      lua_pop(L, 1);
    }
    else
    {
      lua_pushfstring (L, "Error in reading complex table\n");
      lua_error(L);
    }
  }
  return buff;
}


MatrixType* luaL_checkMatrix(lua_State *L, int narg)
{
  uint32_t len_row;
  int i,j;
  MatrixType *buff;
  Complex16 * row;
  luaL_checktype(L, narg, LUA_TTABLE);
  uint32_t len = (uint32_t)lua_objlen(L, narg);
  buff=malloc(sizeof(MatrixType));
  buff->QComplex=1;
  buff->Nj=len;
  for(i=0;i<len;i++)
  {
    lua_rawgeti(L, narg, i+1);
    row = luaL_checkarray_complex16(L, -1, &len_row);
    if(i==0)
    {
      buff->Ni=len_row;
      InitMatrix(buff);
    }
    lua_pop(L, 1);
    if(len_row!=buff->Ni)
    {
      luaL_error(L, "Error while reading matrix. Dimensions of rows not constant. %d at first row %d at row %d",buff->Ni,len_row,i);
    }
    for(j=0;j<buff->Ni;j++)
    {
      buff->VRe[j][i]=row[j].Re;
      buff->VIm[j][i]=row[j].Im;
    }
    free(row);
  }
  MatrixChop(buff, 100*DBL_MIN);
  return buff;
}


WaveFunctionType* luaL_checkarray_WaveFunctionType(lua_State *L, int narg, uint32_t *len_out)
{
  int i;
  luaL_checktype(L, narg, LUA_TTABLE);
  uint32_t len = (uint32_t)lua_objlen(L, narg);
  *len_out = len;
  WaveFunctionType *buff = (WaveFunctionType*)malloc(len*sizeof(WaveFunctionType));
  for(i=0;i<len;i++)
  {
    lua_rawgeti(L, narg, i+1);
    buff[i] = *(WaveFunctionType *)luaL_checkudata(L, -1, "ssp.wavefunction");
    lua_pop(L, 1);
  }
  return buff;
}

OperatorType* luaL_checkarray_OperatorType(lua_State *L, int narg, uint32_t *len_out)
{
  int i;
  luaL_checktype(L, narg, LUA_TTABLE);
  uint32_t len = (uint32_t)lua_objlen(L, narg);
  *len_out = len;
  OperatorType *buff = (OperatorType*)malloc(len*sizeof(OperatorType));
  for(i=0;i<len;i++)
  {
    lua_rawgeti(L, narg, i+1);
    buff[i] = *(OperatorType *)luaL_checkudata(L, -1, "ssp.operator");
    lua_pop(L, 1);
  }
  return buff;
}

HarmonicsExpansion* luaL_checkarray_HarmonicsExpansion(lua_State *L, int narg, uint32_t *len_out)
{
  int i,lencheck;
  Complex16 temp;
  luaL_checktype(L, narg, LUA_TTABLE);
  uint32_t len = (uint32_t)lua_objlen(L, narg);
  *len_out = len;
  HarmonicsExpansion *buff = malloc(sizeof(HarmonicsExpansion));
  buff->N=len;
  buff->QComplex=0;
  buff->Im = malloc(len*sizeof(double));
  buff->Re = malloc(len*sizeof(double));
  buff->l = malloc(len*sizeof(int));
  buff->m = malloc(len*sizeof(int));
  for(i=0;i<len;i++)
  {
    lua_rawgeti(L, narg, i+1);
    luaL_checktype(L, -1, LUA_TTABLE);
    lencheck = (int)lua_objlen(L, -1);
    if(lencheck!=3)
    {
      lua_pushfstring (L, "Error Harmonics Expansion, found %d indices (expected 3 {l,m,val}) at level %d\n",lencheck,i+1);
      lua_error(L);
    }
    else
    {
      lua_rawgeti(L, -1, 1);
      buff->l[i]=(int)luaL_checkunsigned(L, -1);
      lua_pop(L, 1);
      lua_rawgeti(L, -1, 2);
      buff->m[i]=(int)luaL_checkinteger(L, -1);
      lua_pop(L, 1);
      if(Abs(buff->m[i])>buff->l[i])
      {
        luaL_error (L, "Error Harmonics Expansion, |m| (%d) > l (%d) at level %d\n",buff->m[i],buff->l[i],i+1);
      }
      lua_rawgeti(L, -1, 3);
      if(lua_isnumber(L, -1))
      {
        buff->Re[i] = (double)luaL_checknumber(L, -1);
        buff->Im[i]=0;
      }
      else
      {
        temp=*(Complex16 *)luaL_checkudata(L, -1, "Lua.complex");
        buff->Re[i]=temp.Re;
        buff->Im[i]=temp.Im;
        if(Abs(temp.Im)>10*(DBL_EPSILON*Abs(temp.Re)+DBL_MIN)) buff->QComplex=1;
      }
      lua_pop(L, 1);
    }
    lua_pop(L, 1);
  }
  return buff;
}

void LuaPushTable_HarmonicsExpansion(lua_State *L, HarmonicsExpansion *Akm)
{
  uint32_t i;
  Complex16 *res;
  lua_createtable(L, Akm->N, 0);
  for(i=0;i<Akm->N;i++)
  {
    lua_createtable(L, 3, 0);
    {
      lua_pushinteger(L, Akm->l[i]);
      lua_rawseti(L, -2, 1);
      lua_pushinteger(L, Akm->m[i]);
      lua_rawseti(L, -2, 2);
      if(Akm->QComplex)
      {
        res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
        res->Re=Akm->Re[i];
        res->Im=Akm->Im[i];
        luaL_setmetatable(L, "Lua.complex");
      }
      else
      {
        lua_pushnumber(L, Akm->Re[i]);
      }
      lua_rawseti(L, -2, 3);
    }
    lua_rawseti(L, -2, i+1);
  }
}

void shallowtablecopy(lua_State* L, int index)
{
  lua_newtable(L);
  lua_pushnil(L);
  while(lua_next(L, index) != 0)
  {
    lua_pushvalue(L, -2);
    lua_insert(L, -2);
    lua_settable(L, -4);
  }
}

void * LuaReturnVariableType(lua_State * L, int narg, char * type)
{
  int *buffint;
  double *buffdouble;
  void *buff;
  if (strcmp(type, "char*") == 0)
  {
    buff = (void*) luaL_checkstring(L, narg);
    return  buff;
  }
  if (strcmp(type, "int") == 0)
  {
    buffint = malloc(sizeof(int));
    *buffint = luaL_checkint(L, narg);
    return (void*) buffint;
  }
  if (strcmp(type, "double") == 0)
  {
    buffdouble = malloc(sizeof(double));
    *buffdouble = (double)luaL_checknumber(L, narg);
    return (void*) buffdouble;
  }
  
  
  
  luaL_error (L, "Error while reading variable type from stack. Type \"%s\" not recognized\n",type);
  return buff;
}


void * LuaFindOptionsInTable(lua_State * L, int narg, char* token, char* type)
{
  void * buff;
  char * tokenread;
  buff=NULL;
  luaL_checktype(L, narg, LUA_TTABLE);
  uint32_t len,lenelement;
  len = (uint32_t)lua_objlen(L, narg);
  uint32_t i;
  for(i=1;i<=len;i++)
  {
    lua_rawgeti(L, narg, i);
    luaL_checktype(L, -1, LUA_TTABLE);
    lenelement = (int)lua_objlen(L, -1);
    if(lenelement!=2) luaL_error (L, "Error in options, found %d indices (expected 2 {name,value} at level %d\n",lenelement,i);
    lua_rawgeti(L, -1, 1);
    tokenread = (char*) luaL_checkstring(L, -1);
    if (strcmp(tokenread, token) == 0)
    {
      lua_rawgeti(L, -2, 2);
      return LuaReturnVariableType(L, -1, type);
    }
    lua_pop(L,2);
  }
  return buff;
}
