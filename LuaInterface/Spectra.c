static int LuaSpectrumPrint(lua_State * L)
{
  struct
  {
    char *stream;
  } Options;
  uint32_t len = (uint32_t)lua_gettop(L);
  if(len==0)
  {
    PrintSpectra(((SpectraTypeArray *)currentudata)->Spec, ((SpectraTypeArray *)currentudata)->NSpectra);
  }
  else if(lua_istable(L,-1))
  {
    Options.stream = (char*)LuaFindOptionsInTable(L, -1, "file", "char*");
    
    printf("%s\n",Options.stream);
    
    if(Options.stream==NULL)
    {
      PrintSpectra(((SpectraTypeArray *)currentudata)->Spec, ((SpectraTypeArray *)currentudata)->NSpectra);
    }
    else
    {
      FilePrintSpectra(((SpectraTypeArray *)currentudata)->Spec, ((SpectraTypeArray *)currentudata)->NSpectra, Options.stream, "w");
    }
  }
  else
  {
    return luaL_error(L, "Error in .Print(arg) unrecognized arguments in Print()\n");
  }
  return 0;
}

static int LuaSpectraAdd(lua_State * L)
{

  return 1;
}

static int LuaSpectraSub(lua_State * L)
{

  return 1;
}

static int LuaSpectraMul(lua_State * L)
{

  return 1;
}

static int LuaSpectraIndex(lua_State * L)
{
  SpectraTypeArray *a = (SpectraTypeArray*) luaL_checkudata(L, 1, "ssp.spectra");
  const char *str = lua_tostring(L, 2);
  currentudata = a;
  if (strcmp(str, "Print") == 0)
  {
    lua_pushcfunction(L, LuaSpectrumPrint);
    return 1;
  }
  return luaL_error(L, "Error operation \".%s()\" not defined on wavefunction\n",str);
}

static int LuaSpectraToString(lua_State * L)
{
  SpectraTypeArray *spec = (SpectraTypeArray*) luaL_checkudata(L, 1, "ssp.spectra");
  lua_pushstring(L, spec->Spec[0].Name);
  return 1;
}

static int LuaSpectraDestroy(lua_State * L)
{
  uint32_t i;
  SpectraTypeArray *spec = (SpectraTypeArray*) luaL_checkudata(L, 1, "ssp.spectra");
  for(i=0;i<spec->NSpectra;i++)
  {
    SpectrumFree(&spec->Spec[i]);
  }
  return 0;
}

static const struct luaL_Reg Spectra_meta[] = {
  {"__add", LuaSpectraAdd},
  {"__sub", LuaSpectraSub},
  {"__mul", LuaSpectraMul},
  {"__index", LuaSpectraIndex},
  {"__tostring", LuaSpectraToString},
  {"__gc", LuaSpectraDestroy},
  {NULL, NULL}
};