// This file is not used in the Lua version of Quanty.
// Older C based routines can use functions in this file
// Do not add changes here, this file will be removed in the future

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
  // copy string found at the pointer on the lua stack to newly allocated memory
  ret=(char*)malloc(2024*((strlen(in)/2024)+2));
  strcpy(ret, in);
  if(strcmp(font, "Times-Roman")==0)
  {
    AddTextPS(&ret, text, x, y, size, align);
  }
  else if(strcmp(font, "Symbol")==0)
  {
    AddGreekTextPS(&ret, text, x, y, size, align);
  }
  else
  {
    luaL_error(L, "Unrecognized font in PostScript.Text (%s)\n",font);
  }
  lua_pushstring(L, ret);
  free(ret);
  return 1;
}
/*
 *
 */
static int LuaPSLine(lua_State * L)
{
  const char *in;
  char* ret;
  uint32_t RGBColor=0;
  uint32_t Nx,Ny,NRGB;
  double thickness=1;
  double *x,*y;
  uint32_t *RGBTable;
  int ilen;
  int len = (int)lua_gettop(L);
  if(len<3) return luaL_error(L,"PostScript.Line requires minimal 3 arguments and not %d\n", len);
  in = luaL_checkstring(L, 1);
  x = luaL_checkarray_double(L, 2, &Nx);
  y = luaL_checkarray_double(L, 3, &Ny);
  if(Nx!=Ny) return luaL_error(L,"PostScript.Line requires a table at position 2 and 3 of the same lenght (currently Nx=%d, Ny=%d)\n", Nx,Ny);
  for(ilen=4;ilen<=len;ilen++)
  {
    if(lua_isnumber(L, ilen)) thickness=(double)luaL_checknumber(L, ilen);
    if(lua_istable(L, ilen))
    {
      RGBTable = luaL_checkarray_unit32_t(L, ilen, &NRGB);
      if(NRGB!=3) return luaL_error(L,"PostScript.Line Color should be a table with 3 values (RGB) between 0 and 255\n");
      RGBColor=((((RGBTable[0]&255)<<8)+(RGBTable[1]&255))<<8)+(RGBTable[2]&255);
      free(RGBTable);
    }
  }
  // copy string found at the pointer on the lua stack to newly allocated memory
  ret=(char*)malloc(2024*((strlen(in)/2024)+2));
  strcpy(ret, in);
  AddLinePathPS(&ret, x, y, Nx, RGBColor, thickness);
  lua_pushstring(L, ret);
  free(ret);
  free(x);
  free(y);
  return 1;
}
/*
 *
 */
static int LuaPSPolygon(lua_State * L)
{
  const char *in;
  char* ret;
  uint32_t RGBColor=0;
  uint32_t Nx,Ny,NRGB;
  double thickness=1;
  double *x,*y;
  uint32_t *RGBTable;
  int ilen;
  int len = (int)lua_gettop(L);
  if(len<3) return luaL_error(L,"PostScript.Polygon requires minimal 3 arguments and not %d\n", len);
  in = luaL_checkstring(L, 1);
  x = luaL_checkarray_double(L, 2, &Nx);
  y = luaL_checkarray_double(L, 3, &Ny);
  if(Nx!=Ny) return luaL_error(L,"PostScript.Polygon requires a table at position 2 and 3 of the same lenght (currently Nx=%d, Ny=%d)\n", Nx,Ny);
  for(ilen=4;ilen<=len;ilen++)
  {
    if(lua_isnumber(L, ilen)) thickness=(double)luaL_checknumber(L, ilen);
    if(lua_istable(L, ilen))
    {
      RGBTable = luaL_checkarray_unit32_t(L, ilen, &NRGB);
      if(NRGB!=3) return luaL_error(L,"PostScript.Polygon Color should be a table with 3 values (RGB) between 0 and 255\n");
      RGBColor=((((RGBTable[0]&255)<<8)+(RGBTable[1]&255))<<8)+(RGBTable[2]&255);
      free(RGBTable);
    }
  }
  // copy string found at the pointer on the lua stack to newly allocated memory
  ret=(char*)malloc(2024*((strlen(in)/2024)+2));
  strcpy(ret, in);
  AddPolygonPS(&ret, x, y, Nx, RGBColor, thickness);
  lua_pushstring(L, ret);
  free(ret);
  free(x);
  free(y);
  return 1;
}

static int LuaPSBitmap(lua_State * L)
{
  const char *in;
  char* ret;
  double xmin, ymin, xmax, ymax;
  uint32_t *RGBColor;
  CompactMatrixType **mat=NULL;
  Complex16 *temp;
  int i,j,Nmat=0,Red,Green,Blue;
  // check intput types and length
  int len = (int)lua_gettop(L);
  if(len!=7) return luaL_error(L,"PostScript.Bitmap requires 7 arguments and not %d (String to add to, xmin, ymin, xmax, ymax, matrix or list of matrices, collorfunction\n", len);
  // read input
  in = luaL_checkstring(L, 1);
  xmin=(double)luaL_checknumber(L, 2);
  ymin=(double)luaL_checknumber(L, 3);
  xmax=(double)luaL_checknumber(L, 4);
  ymax=(double)luaL_checknumber(L, 5);
  // on position 6 we find either a compact matrix or a table of compact matrices
  if (luaL_testudata(L, 6, "Matrix_Type") != NULL)
  {
    Nmat=1;
    mat=(CompactMatrixType**)malloc(Nmat*sizeof(CompactMatrixType*));
    mat[0] = luaL_checkCompactMatrix(L, 6);
  }
  else if(lua_istable(L, 6))
  {
    lua_rawgeti(L, 6, 1);
    if(lua_istable(L, -1))
    {
      lua_rawgeti(L, -1, 1);
      if(lua_isnumber(L, -1))
      {
        Nmat=1;
        mat=(CompactMatrixType**)malloc(Nmat*sizeof(CompactMatrixType*));
        mat[0] = luaL_checkCompactMatrix(L, 6);
      }
      lua_pop(L, 1);
    }
    lua_pop(L, 1);
  }
  else
  {
    Nmat = (int)lua_objlen(L, 6);
    mat=(CompactMatrixType**)malloc(Nmat*sizeof(CompactMatrixType*));
    for(i=0;i<Nmat;i++)
    {
      lua_rawgeti(L, 6, i+1);
      mat[i] = luaL_checkCompactMatrix(L, -1);
      lua_pop(L, 1);
    }
  }
  // test dimensions
  for(i=1;i<Nmat;i++)
  {
    if((mat[i]->Ni!=mat[0]->Ni) || (mat[i]->Nj!=mat[0]->Nj)) return luaL_error(L, "Dimensions of matrices that specify the color must all be the same\n");
  }
  // copy string found at the pointer on the lua stack to newly allocated memory
  ret=(char*)malloc(2024*((strlen(in)/2024)+2));
  strcpy(ret, in);
  // use the function on position 7 to transfer the list of Compact matrices to RGB values
  RGBColor=(uint32_t*)malloc(mat[0]->Nj*mat[0]->Ni*sizeof(uint32_t));
  for(i=0;i<mat[0]->Nj*mat[0]->Ni;i++)
  {
    lua_pushvalue(L, 7);
    for(j=0;j<Nmat;j++)
    {
      if(mat[j]->QComplex==0)
      {
        lua_pushnumber(L, mat[j]->VRe[i]);
      }
      else
      {
        temp = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
        temp->Re = mat[j]->V[i].Re;
        temp->Im = mat[j]->V[i].Im;
        luaL_setmetatable(L, "Complex_Type");
      }
    }
    if(lua_pcall(L, Nmat, 3, 0))
    {
      return luaL_error(L, "Incorrect color function call\n");
    }
    Red=(int)luaL_checknumber(L, 8);
    Green=(int)luaL_checknumber(L, 9);
    Blue=(int)luaL_checknumber(L, 10);
    lua_pop(L, 3);
    if(Red<0) Red=0;
    if(Green<0) Green=0;
    if(Blue<0) Blue=0;
    if(Red>255) Red=255;
    if(Green>255) Green=255;
    if(Blue>255) Blue=255;
    RGBColor[i]=(((Red<<8)+Green)<<8)+Blue;
  }
  
  // create bitmap
  AddBitmapPS(&ret, xmin, ymin, xmax-xmin, ymax-ymin, RGBColor, mat[0]->Nj, mat[0]->Ni);
  lua_pushstring(L, ret);
  free(ret);
  for(i=0;i<Nmat;i++)
  {
    if(mat[i]->OnLuaStack==false)
    {
      FreeCompactMatrix(mat[i]);
      free(mat[i]);
    }
    
  }
  free(mat);
  free(RGBColor);
  return 1;
}

static const struct luaL_Reg PostScript_func[] = {
  {"Open", LuaPSOpen},
  {"Close", LuaPSClose},
  {"Text", LuaPSText},
  {"Line", LuaPSLine},
  {"Polygon", LuaPSPolygon},
  {"Bitmap", LuaPSBitmap},
  {NULL, NULL}
};
