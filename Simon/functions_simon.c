/* Some makeshift functions */

static int LuaSimontest(lua_State * L)
{
  lua_pushstring(L, "Hello world");
  return 1;
}

static const struct luaL_Reg Simon_meta[] = {
  {NULL, NULL}
};

static const struct luaL_Reg Simon_func[] = {
  {"test", LuaSimontest},
  {NULL, NULL}
};
