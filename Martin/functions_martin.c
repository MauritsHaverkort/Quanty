/* Some makeshift functions */

static int Luatest(lua_State * L)
{
  lua_pushstring(L, "Hello world");
  return 1;
}

static const struct luaL_Reg Martin_meta[] = {
  {NULL, NULL}
};

static const struct luaL_Reg Martin_func[] = {
  {"test", Luatest},
  {NULL, NULL}
};
