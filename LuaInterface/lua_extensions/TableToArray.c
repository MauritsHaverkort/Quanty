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
        lua_pushfstring (L, "Error Harmonics Expansion, |m| (%d) > l (%d) at level %d\n",buff->m[i],buff->l[i],i+1);
        lua_error(L);
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