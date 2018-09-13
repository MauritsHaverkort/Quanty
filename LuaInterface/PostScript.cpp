/*
 *
 */
static int LuaPSOpen(lua_State * L)
{
  char* ret;
  int x,y;
  int len = (int)lua_gettop(L);
  if(len==0)
  {
    x=200;
    y=200;
  }
  else if(len==2)
  {
    x = (int)luaL_checknumber(L, 1);
    y = (int)luaL_checknumber(L, 2);
  }
  else
  {
    return luaL_error(L, "Error in PostScript.Open expected 0 or 2 arguments, not %d\n",len);
  }
  CreatePS(&ret,x,y);
  lua_pushstring(L, ret);
  free(ret);
  return 1;
}
/*
 *
 */
static int LuaPSClose(lua_State * L)
{
  char* ret;
  const char *in;
  in = luaL_checkstring(L, 1);
  ret=(char*)malloc(2024*((strlen(in)/2024)+2));
  strcpy(ret, in);
  ClosePS(&ret);
  lua_pushstring(L, ret);
  free(ret);
  return 1;
}
/*
 *
 */
static int LuaPSText(lua_State * L)
{
  int ilen;
  char* ret;
  const char* text;
  double x,y,size;
  char align;
  char font[]="Times-Roman";
  const char *in;
  const char *fontp=font;
  const char *alignp;
  int len = (int)lua_gettop(L);
  size=10;
  align='l';
  if(len<4) return luaL_error(L,"PostScript.Text requires minimal 4 arguments and not %d\n", len);
  in = luaL_checkstring(L, 1);
  text = luaL_checkstring(L, 2);
  x = (double)luaL_checknumber(L, 3);
  y = (double)luaL_checknumber(L, 4);
  for(ilen=5;ilen<=len;ilen++)
  {
    if(lua_isnumber(L, ilen))
    {
      size = (double)luaL_checknumber(L, ilen);
    }
    else
    {
      alignp = luaL_checkstring(L, ilen);
      if(strlen(alignp)==1)
      {
        align=alignp[0];
      }
      else
      {
        fontp=alignp;
      }
    }
  }
  
  
  
  ret=(char*)malloc(2024*((strlen(in)/2024)+2));
  strcpy(ret, in);
  
  
  
  
  ClosePS(&ret);
  lua_pushstring(L, ret);
  free(ret);
  return 1;
}
/*
 *
 */
static int LuaPSLine(lua_State * L)
{
  
  return 1;
}
/*
 *
 */
static int LuaPSPolygon(lua_State * L)
{
  
  return 1;
}





static const struct luaL_Reg PostScript_func[] = {
  {"Open", LuaPSOpen},
  {"Close", LuaPSClose},
  {"Text", LuaPSText},
  {"Line", LuaPSLine},
  {"Polygon", LuaPSPolygon},
  {NULL, NULL}
};
