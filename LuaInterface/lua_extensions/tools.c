static void luaMWH_StackDump (lua_State *L)
{
  int i,len;
  int top = lua_gettop(L);
  if(top<=0)
  {
    printf("empty stack (top=%d) at location %p\n",top,L);
  }
  printf("\nLua Stackdump");
  for (i = 1; i <= top; i++) {  /* repeat for each level */
    int t = lua_type(L, i);
    printf("\n%s - ",lua_typename(L, t));
    switch (t) {

      case LUA_TSTRING:  /* strings */
        printf("`%s'", lua_tostring(L, i));
        break;

      case LUA_TBOOLEAN:  /* booleans */
        printf(lua_toboolean(L, i) ? "true" : "false");
        break;

      case LUA_TNUMBER:  /* numbers */
        printf("%g", lua_tonumber(L, i));
        break;

      case LUA_TTABLE:  /* table */
        len = (int)lua_objlen(L, i);
        printf("%d entries", len);
        break;

      case LUA_TUSERDATA:  /* userdata */
        if (luaL_testudata(L, i, "Complex_Type") != NULL) printf("Complex");
        else if (luaL_testudata(L, i, "Operator_Type") != NULL) printf("Operator");
        else if (luaL_testudata(L, i, "Wavefunction_Type") != NULL) printf("Wavefunction");
        else if (luaL_testudata(L, i, "Spectra_Type") != NULL) printf("Spectra");
        else if (luaL_testudata(L, i, "TightBinding_Type") != NULL) printf("TightBinding");
        else if (luaL_testudata(L, i, "StuttgartCTRL_Type") != NULL) printf("StuttgartCTRL");
        else if (luaL_testudata(L, i, "StuttgartHAMR_Type") != NULL) printf("StuttgartHAMR");
        else if (luaL_testudata(L, i, "DresdenFPLO_Type") != NULL) printf("DresdenFPLO");
        else printf("Unrecognized");
        break;

      default:  /* other values */
        printf("%s", lua_typename(L, t));
        break;

    }
  }
  printf("\n");  /* end the listing */
  fflush(stdout);
}

static void luaMWH_RegisterGlobalFunctions (lua_State *L, const luaL_Reg *l)
{
  for (; l->name != NULL; l++)
  {  /* loop over all entries in regestry */
    lua_pushcfunction(L, l->func);
    lua_setglobal(L, l->name);
  }
}

uint32_t luaMWH_fulltablelen(lua_State *L, const int t)
{
  uint32_t len=0;
  lua_pushnil(L);
  while (lua_next(L, t) != 0)
  {
    lua_pop(L, 1);
    len++;
  }
  lua_pop(L, 1);
  return len;
}


