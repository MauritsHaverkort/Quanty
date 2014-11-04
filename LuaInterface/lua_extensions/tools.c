static void luaL_stackDump (lua_State *L) {
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
        if (luaL_testudata(L, i, "Lua.complex") != NULL) printf("complex");
        else if (luaL_testudata(L, i, "ssp.operator") != NULL) printf("operator");
        else if (luaL_testudata(L, i, "ssp.wavefunction") != NULL) printf("wavefunction");
        else if (luaL_testudata(L, i, "ssp.spectra") != NULL) printf("spectra");
        else printf("unrecognized");
        break;
        
      default:  /* other values */
        printf("%s", lua_typename(L, t));
        break;
        
    }
  }
  printf("\n");  /* end the listing */
  fflush(stdout);
}