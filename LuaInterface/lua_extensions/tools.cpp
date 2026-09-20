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
        if      (luaL_testudata(L, i, "Complex_Type")          != NULL) printf("Complex");
        else if (luaL_testudata(L, i, "Operator_Type")         != NULL) printf("Operator");
        else if (luaL_testudata(L, i, "Wavefunction_Type")     != NULL) printf("Wavefunction");
        else if (luaL_testudata(L, i, "Spectra_Type")          != NULL) printf("Spectra");
        else if (luaL_testudata(L, i, "TightBinding_Type")     != NULL) printf("TightBinding");
        else if (luaL_testudata(L, i, "StuttgartCTRL_Type")    != NULL) printf("StuttgartCTRL");
        else if (luaL_testudata(L, i, "StuttgartHAMR_Type")    != NULL) printf("StuttgartHAMR");
        else if (luaL_testudata(L, i, "DresdenFPLO_Type")      != NULL) printf("DresdenFPLO");
        else if (luaL_testudata(L, i, "Matrix_Type")           != NULL) printf("Matrix");
        else if (luaL_testudata(L, i, "ResponseFunction_Type") != NULL) printf("ResponseFunction");
        else if (luaL_testudata(L, i, "Graphics_Type")         != NULL) printf("Graphics_Type");
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



// PATH_MAX should come from limits.h, but is not guaranteed by the standard
#ifndef PATH_MAX
#define PATH_MAX 4096
#endif

/**
 * Return the absolute path of FileName in a newly allocated string.
 * The caller is responsible for freeing the string.
 * Returns NULL if the path can not be resolved, for example if the file does
 * not exist. Note that on unix like systems symbolic links are resolved,
 * while on windows they are not.
 */
char * MWHRealPath(const char * FileName)
{
  char * AbsName;
#if defined(_WIN32)
  if ((AbsName = (char*)malloc(_MAX_PATH * sizeof(char))) == NULL) return NULL;
  if (_fullpath(AbsName, FileName, _MAX_PATH) == NULL) { free(AbsName); return NULL; }
#else
  if ((AbsName = (char*)malloc(PATH_MAX * sizeof(char))) == NULL) return NULL;
  if (realpath(FileName, AbsName) == NULL) { free(AbsName); return NULL; }
#endif
  return AbsName;
}

/**
 * Return the directory that contains FileName in a newly allocated string.
 * The path is absolute and ends on a directory separator, such that the name
 * of a file next to FileName is found by simple concatenation.
 * The caller is responsible for freeing the string.
 * Returns NULL if the path can not be resolved.
 */
char * MWHFileDirectory(const char * FileName)
{
  char * AbsName;
  size_t i, last;
  if ((AbsName = MWHRealPath(FileName)) == NULL) return NULL;
  last = 0;
  for (i = 0; AbsName[i] != '\0'; i++)
  {
    if ((AbsName[i] == '/') || (AbsName[i] == '\\')) last = i + 1;
  }
  AbsName[last] = '\0';
  return AbsName;
}
