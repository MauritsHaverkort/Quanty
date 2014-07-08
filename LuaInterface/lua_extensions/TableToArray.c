/*  returns a pointer to an array of lenght len_out
 *  of type unit16_t
 *  found in the lua stack at position narg
 */
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