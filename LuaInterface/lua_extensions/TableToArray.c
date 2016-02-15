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
    else if (luaL_testudata(L, -1, "Complex_Type") != NULL)
    {
      buff[i] = *(Complex16 *)luaL_checkudata(L, -1, "Complex_Type");
      lua_pop(L, 1);
    }
    else
    {
      luaL_error(L, "Error in reading complex table\n");
    }
  }
  return buff;
}

CompactVectorType luaL_checkVector(lua_State *L, int narg)
{
  uint32_t i;
  CompactVectorType buff;
  buff.V=luaL_checkarray_complex16(L, narg, &(buff.Ni));
  buff.QComplex=0;
  for(i=0;i<buff.Ni;i++)
  {
    if(Abs(buff.V[i].Im)>(10*DBL_MIN+10*Abs(buff.V[i].Re)*DBL_EPSILON)) buff.QComplex=1;
  }
  if(buff.QComplex==0)
  {
    double *VRe = malloc(buff.Ni*sizeof(double));
    for(i=0;i<buff.Ni;i++)
    {
      VRe[i]=buff.V[i].Re;
    }
    free(buff.V);
    buff.VRe=VRe;
  }
  return buff;
}

CompactMatrixType* luaL_checkCompactMatrix(lua_State *L, int narg)
{
  uint32_t len_row;
  int i,j;
  CompactMatrixType *buff;
  Complex16 * row;
  luaL_checktype(L, narg, LUA_TTABLE);
  uint32_t len = (uint32_t)lua_objlen(L, narg);
  buff=malloc(sizeof(CompactMatrixType));
  buff->QComplex=1;
  buff->Ni=len;
  for(i=0;i<len;i++)
  {
    lua_rawgeti(L, narg, i+1);
    row = luaL_checkarray_complex16(L, -1, &len_row);
    if(i==0)
    {
      buff->Nj=len_row;
      InitCompactMatrix(buff);
    }
    lua_pop(L, 1);
    if(len_row!=buff->Nj)
    {
      luaL_error(L, "Error while reading matrix. Dimensions of rows not constant. %d at first row %d at row %d",buff->Nj,len_row,i);
    }
    for(j=0;j<buff->Nj;j++)
    {
      buff->V[i*buff->Nj+j].Re=row[j].Re;
      buff->V[i*buff->Nj+j].Im=row[j].Im;
    }
    free(row);
  }
  CompactMatrixChop(buff, 100*DBL_MIN);
  return buff;
}


MatrixType* luaL_checkMatrix(lua_State *L, int narg)
{
  uint32_t len_row;
  int i,j;
  MatrixType *buff;
  Complex16 * row;
  luaL_checktype(L, narg, LUA_TTABLE);
  uint32_t len = (uint32_t)lua_objlen(L, narg);
  buff=malloc(sizeof(MatrixType));
  buff->QComplex=1;
  buff->Ni=len;
  for(i=0;i<len;i++)
  {
    lua_rawgeti(L, narg, i+1);
    row = luaL_checkarray_complex16(L, -1, &len_row);
    if(i==0)
    {
      buff->Nj=len_row;
      InitMatrix(buff);
    }
    lua_pop(L, 1);
    if(len_row!=buff->Nj)
    {
      luaL_error(L, "Error while reading matrix. Dimensions of rows not constant. %d at first row %d at row %d",buff->Nj,len_row,i);
    }
    for(j=0;j<buff->Nj;j++)
    {
      buff->VRe[i][j]=row[j].Re;
      buff->VIm[i][j]=row[j].Im;
    }
    free(row);
  }
  MatrixChop(buff, 100*DBL_MIN);
  return buff;
}


WaveFunctionType* luaL_checkarray_WaveFunctionType(lua_State *L, int narg, uint32_t *len_out)
{
  int i,j;
  luaL_checktype(L, narg, LUA_TTABLE);
  uint32_t len = (uint32_t)lua_objlen(L, narg);
  *len_out = len;
  WaveFunctionType *buff = (WaveFunctionType*)malloc(len*sizeof(WaveFunctionType));
  for(i=0;i<len;i++)
  {
    lua_rawgeti(L, narg, i+1);
    buff[i] = *(WaveFunctionType *)luaL_checkudata(L, -1, "Wavefunction_Type");
    lua_pop(L, 1);
    for(j=0;j<i;j++)
    {
      if(buff[i].ValueRe==buff[j].ValueRe)
      {
        WaveFunctionCopy(buff[j], &buff[i]);
      }
    }
  }
  return buff;
}

OperatorType* luaL_checkarray_OperatorType(lua_State *L, int narg, uint32_t *len_out)
{
  int i,j;
  luaL_checktype(L, narg, LUA_TTABLE);
  uint32_t len = (uint32_t)lua_objlen(L, narg);
  *len_out = len;
  OperatorType *buff = (OperatorType*)malloc(len*sizeof(OperatorType));
  for(i=0;i<len;i++)
  {
    lua_rawgeti(L, narg, i+1);
    buff[i] = *(OperatorType *)luaL_checkudata(L, -1, "Operator_Type");
    lua_pop(L, 1);
    for(j=0;j<i;j++)
    {
      if(buff[i].l==buff[j].l)
      {
        OperatorCopy(&buff[i], buff[j]);
      }
    }
  }
  return buff;
}

void luaL_setarray_DeterminatRestrictionType(lua_State *L, int narg, DeterminantRestrictionType * buff)
{
  uint32_t NFermion,NBoson,i;
  uint32_t len,len2,detleng,NWord;
  char * detchar;
  if(lua_isnil(L,narg))
  {
    buff->N=0;
    return;
  }
  luaL_checktype(L, narg, LUA_TTABLE);
  len = (uint32_t)lua_objlen(L, narg);
  if(len<3) luaL_error(L, "Restrictions should be given in the format {NFermion,NBoson,{\"determinant1\",minocc,maxocc},{\"determinant2\",minocc,maxocc},...}");
  buff->N=len-2;
  buff->Min=malloc(buff->N*sizeof(uint32_t));
  buff->Max=malloc(buff->N*sizeof(uint32_t));
  buff->Determinant=malloc(buff->N*sizeof(uint8_t*));
  lua_rawgeti(L, narg, 1);
  NFermion=(uint32_t)luaL_checkunsigned(L, -1);
  lua_pop(L, 1);
  lua_rawgeti(L, narg, 2);
  NBoson=(uint32_t)luaL_checkunsigned(L, -1);
  lua_pop(L, 1);
  NWord=NBoson + Min(NFermion, 1 + (NFermion - 1) / 8);
  for(i=0;i<buff->N;i++)
  {
    lua_rawgeti(L, narg, i+3);
    {
      luaL_checktype(L, -1, LUA_TTABLE);
      len2 = (uint32_t)lua_objlen(L, -1);
      if(len2!=3) luaL_error(L, "Restrictions should be given in the format {NFermion,NBoson,{\"determinant1\",minocc,maxocc},{\"determinant2\",minocc,maxocc},...}");
      lua_rawgeti(L, -1, 1);
      {
        if(lua_isstring(L, -1))
        {
          detchar = (char *) luaL_checkstring(L, -1);
          buff->Determinant[i] = StringToDet(NFermion, NBoson, detchar);
          detleng = NWord;
        }
        else
        {
          buff->Determinant[i]=luaL_checkarray_unit8_t(L, -1, &detleng);
        }
        if(detleng!=NWord)
        {
          luaL_error(L, "Error in Restrictions, determinant number %d should be a table of %d words\n",i,NBoson + Min(NFermion, 1 + (NFermion - 1) / 8));
        }
      }
      lua_pop(L, 1);
      lua_rawgeti(L, -1, 2);
      {
        buff->Min[i]=luaL_checkunsigned(L, -1);
      }
      lua_pop(L, 1);
      lua_rawgeti(L, -1, 3);
      {
        buff->Max[i]=luaL_checkunsigned(L, -1);
      }
      lua_pop(L, 1);
    }
    lua_pop(L, 1);
  }
  return;
}

DeterminantRestrictionType* luaL_checkarray_DeterminantRestrictionType(lua_State *L, int narg)
{
  DeterminantRestrictionType *buff;
  buff = malloc(sizeof(DeterminantRestrictionType));
  luaL_setarray_DeterminatRestrictionType(L, narg, buff);
  return buff;
}


HarmonicsExpansion* luaL_checkarray_HarmonicsExpansion(lua_State *L, int narg)
{
  int i,lencheck;
  Complex16 temp;
  luaL_checktype(L, narg, LUA_TTABLE);
  uint32_t len = (uint32_t)lua_objlen(L, narg);
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
      luaL_error(L, "Error Harmonics Expansion, found %d indices (expected 3 {l,m,val}) at level %d\n",lencheck,i+1);
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
        luaL_error (L, "Error Harmonics Expansion, |m| (%d) > l (%d) at level %d\n",buff->m[i],buff->l[i],i+1);
      }
      lua_rawgeti(L, -1, 3);
      if(lua_isnumber(L, -1))
      {
        buff->Re[i] = (double)luaL_checknumber(L, -1);
        buff->Im[i]=0;
      }
      else
      {
        temp=*(Complex16 *)luaL_checkudata(L, -1, "Complex_Type");
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

void LuaPushArray_double(lua_State *L, double const * const val, int const N)
{
  uint32_t i;
  lua_createtable(L, N, 0);
  for(i=0;i<N;i++)
  {
    lua_pushnumber(L, val[i]);
    lua_rawseti(L, -2, i+1);
  }
}


void LuaPushTable_HarmonicsExpansion(lua_State *L, HarmonicsExpansion *Akm)
{
  uint32_t i;
  Complex16 *res;
  lua_createtable(L, Akm->N, 0);
  for(i=0;i<Akm->N;i++)
  {
    lua_createtable(L, 3, 0);
    {
      lua_pushinteger(L, Akm->l[i]);
      lua_rawseti(L, -2, 1);
      lua_pushinteger(L, Akm->m[i]);
      lua_rawseti(L, -2, 2);
      if(Akm->QComplex)
      {
        res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
        res->Re=Akm->Re[i];
        res->Im=Akm->Im[i];
        luaL_setmetatable(L, "Complex_Type");
      }
      else
      {
        lua_pushnumber(L, Akm->Re[i]);
      }
      lua_rawseti(L, -2, 3);
    }
    lua_rawseti(L, -2, i+1);
  }
}

void LuaPushTable_WaveFunctionType(lua_State* L, WaveFunctionType const * const psi, uint32_t const Npsi )
{
  uint32_t i;
  WaveFunctionType *buff;
  lua_createtable(L, Npsi, 0);
  for(i=0;i<Npsi;i++)
  {
    buff = (WaveFunctionType*) lua_newuserdata(L, sizeof(WaveFunctionType));
    luaL_setmetatable(L, "Wavefunction_Type");
    *buff = psi[i];
    lua_rawseti(L, -2, i+1);
  }
}

void LuaPushTable_MatrixType(lua_State* L, MatrixType const * const M)
{
  uint32_t i,j;
  Complex16 *res;
  lua_createtable(L,M->Ni,0);
  for(i=0;i<M->Ni;i++)
  {
    lua_createtable(L,M->Nj,0);
    for(j=0;j<M->Nj;j++)
    {
      if(M->QComplex)
      {
        res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
        res->Re = M->VRe[i][j];
        res->Im = M->VIm[i][j];
        luaL_setmetatable(L, "Complex_Type");
      }
      else
      {
        lua_pushnumber(L, M->VRe[i][j]);
      }
      lua_rawseti(L, -2, j+1);
    }
    lua_rawseti(L, -2, i+1);
  }
}

void LuaPushTable_CompactMatrixType(lua_State* L, CompactMatrixType const * const M)
{
  uint32_t i,j;
  Complex16 *res;
  lua_createtable(L,M->Ni,0);
  for(i=0;i<M->Ni;i++)
  {
    lua_createtable(L,M->Nj,0);
    for(j=0;j<M->Nj;j++)
    {
      if(M->QComplex)
      {
        res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
        res->Re = M->V[i*M->Nj+j].Re;
        res->Im = M->V[i*M->Nj+j].Im;
        luaL_setmetatable(L, "Complex_Type");
      }
      else
      {
        lua_pushnumber(L, M->VRe[i*M->Nj+j]);
      }
      lua_rawseti(L, -2, j+1);
    }
    lua_rawseti(L, -2, i+1);
  }
}

void LuaPushTable_TriDiagonalMatrixType(lua_State* L, TriDiagonalMatrixType const * const M)
{
  uint32_t i,MNtrunc=M->N;
  for(i=1;i<M->N;i++)
  {
    if(Sqr(M->b[i])<0.01*DBL_EPSILON)
    {
      MNtrunc=i;
      break;
    }
  }
  lua_createtable(L, 3, 0);
  /* push a */
  LuaPushArray_double(L, M->a, MNtrunc+1);
  lua_rawseti(L, -2, 1);
  /* push b */
  LuaPushArray_double(L, M->b, MNtrunc);
  lua_rawseti(L, -2, 2);
  /* push mu */
  lua_pushnumber(L, M->mu);
  lua_rawseti(L, -2, 3);
}

/*  if at position len of the lua stack a table is found it is over writen by the table in matrix M
 *  otherwise an error is thrown
 */
void LuaInsertTable_CompactMatrixType(lua_State*L, int len, CompactMatrixType const * const M)
{
  luaL_checktype(L, len, LUA_TTABLE);
  LuaPushTable_CompactMatrixType(L, M);
  lua_replace(L, len);
}


void LuaPushTable_SpectraTypeArray(lua_State* L, SpectraTypeArray const * const M)
{
  uint32_t i,j;
  Complex16 *res;
  lua_createtable(L,M->NSpectra,0);
  for(i=0;i<M->NSpectra;i++)
  {
    lua_createtable(L,M->Spec[i].N+1,4);
    for(j=0;j<=M->Spec[i].N;j++)
    {
      lua_createtable(L,2,0);
      lua_pushnumber(L, M->Spec[i].Emin+j*M->Spec[i].dE);
      lua_rawseti(L, -2, 1);
      res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
      res->Re = M->Spec[i].Re[j];
      res->Im = M->Spec[i].Im[j];
      luaL_setmetatable(L, "Complex_Type");
      lua_rawseti(L, -2, 2);
      lua_rawseti(L, -2, j+1);
    }
    lua_pushnumber(L, M->Spec[i].EminPole);
    lua_setfield(L,-2,"EminPole");
    lua_pushnumber(L, M->Spec[i].EmaxPole);
    lua_setfield(L,-2,"EmaxPole");
    lua_pushnumber(L, M->Spec[i].Gamma);
    lua_setfield(L,-2,"Gamma");
    lua_pushnumber(L, Sqr(M->Spec[i].b[0]));
    lua_setfield(L,-2,"Integral");
    if(strcmp("", M->Spec[i].Name)!=0)
    {
      lua_pushstring(L, M->Spec[i].Name);
      lua_setfield(L, -2, "Name");
    }
    lua_rawseti(L, -2, i+1);
  }
  return;
}


void shallowtablecopy(lua_State* L, int index)
{
  lua_newtable(L);
  lua_pushnil(L);
  while(lua_next(L, index) != 0)
  {
    lua_pushvalue(L, -2);
    lua_insert(L, -2);
    lua_settable(L, -4);
  }
}

void * LuaReturnVariableType(lua_State * L, int narg, char * type)
{
  int *buffint;
  double *buffdouble;
  Bool *buffBoolean;
  DeterminantRestrictionType *buffDeterimantRestrictionType;
  void *buff;
  CompactVectorType *buffVector;
  if (strcmp(type, "char*") == 0)
  {
    buff = (void*) luaL_checkstring(L, narg);
    return  buff;
  }
  if (strcmp(type, "int") == 0)
  {
    buffint = malloc(sizeof(int));
    *buffint = (int)luaL_checkint(L, narg);
    return (void*) buffint;
  }
  if (strcmp(type, "uint32_t") == 0)
  {
    buffint = malloc(sizeof(uint32_t));
    *buffint = (uint32_t)luaL_checkint(L, narg);
    return (void*) buffint;
  }
  if (strcmp(type, "uint16_t") == 0)
  {
    buffint = malloc(sizeof(uint16_t));
    *buffint = (uint16_t)luaL_checkint(L, narg);
    return (void*) buffint;
  }
  if (strcmp(type, "uint8_t") == 0)
  {
    buffint = malloc(sizeof(uint8_t));
    *buffint = (uint8_t)luaL_checkint(L, narg);
    return (void*) buffint;
  }
  if (strcmp(type, "double") == 0)
  {
    buffdouble = malloc(sizeof(double));
    *buffdouble = (double)luaL_checknumber(L, narg);
    return (void*) buffdouble;
  }
  if (strcmp(type, "bool") == 0)
  {
    buffint = malloc(sizeof(int));
    buffBoolean = malloc(sizeof(Bool));
    *buffint = (int)lua_toboolean(L, narg);
    *buffBoolean= ((*buffint==1)? True: False);
    free(buffint);
    return (void*) buffBoolean;
  }
  if( strcmp(type, "Vector") ==0)
  {
    buffVector = malloc(sizeof(CompactVectorType));
    *buffVector = (CompactVectorType)luaL_checkVector(L,narg);
    return (void*) buffVector;
  }
  if (strcmp(type, "DeterminantRestrictionType") == 0)
  {
    buffDeterimantRestrictionType = (DeterminantRestrictionType*)luaL_checkarray_DeterminantRestrictionType(L, narg);
    return (void*) buffDeterimantRestrictionType;
  }
  luaL_error (L, "Error while reading variable type from stack. Type \"%s\" not recognized\n",type);
  return buff;
}

/*  if state narg in L  is a table
 *  and if this table has an element that is a table
 *  and if this table has two elements
 *  and if element one of that table is a string
 *  and if that string matches token
 *
 *  then read a variable of type type from the second position
 *
 *  if one of the conditions fails continue without error
 *  if the type of the variable at the second position is not of type type through an error.
 */
void * LuaFindOptionsInTable(lua_State * L, int narg, char* token, char* type)
{
  void * buff;
  char * tokenread;
  uint32_t i;
  uint32_t len,lenelement;
  buff=NULL;
  if(lua_istable(L, narg))
  {
    len = (uint32_t)lua_objlen(L, narg);
    for(i=1;i<=len;i++)
    {
      lua_rawgeti(L, narg, i);
      if(lua_istable(L, -1))
      {
        lenelement = (int)lua_objlen(L, -1);
        if(lenelement==2)
        {
          lua_rawgeti(L, -1, 1);
          if(lua_isstring(L, -1))
          {
            tokenread = (char*) luaL_checkstring(L, -1);
            if (strcmp(tokenread, token) == 0)
            {
              lua_rawgeti(L, -2, 2);
              buff = LuaReturnVariableType(L, -1, type);
              lua_pop(L,1);
            }
          }
          lua_pop(L,1);
        }
      }
      lua_pop(L,1);
    }
  }
  return buff;
}

/*  test if element at position narg of the stack could be a list of options
 *  returns 1 if:
 *    it is a table
 *    each element of that table is a table with two elements
 *      the first element of each table of length two is a string
 */
int32_t LuaTestOptionsList(lua_State *L, int narg)
{
  uint32_t i,len;
  if(!(lua_istable(L, narg))) return 0;
  len = (uint32_t)lua_objlen(L, narg);
  for(i=1;i<=len;i++)
  {
    lua_rawgeti(L, narg, i);
    if(!(lua_istable(L, -1)))
    {
      lua_pop(L,1);
      return 0;
    }
    if(!(((int)lua_objlen(L, -1))==2))
    {
      lua_pop(L,1);
      return 0;
    }
    lua_rawgeti(L,-1,1);
    if(!(lua_type(L, -1)==LUA_TSTRING))
    {
      lua_pop(L,1);
      lua_pop(L,1);
      return 0;
    }
    lua_pop(L,1);
    lua_pop(L,1);
  }
  return 1;
}

int32_t LuaCheckOptions(lua_State *L, int narg, LuaOptionsType *Options)
{
  char *tokenread, *Buffer;
  uint32_t i,j, len, tokenfound;
  /* check if L has a table at position narg and set len to the length of this table */
  luaL_checktype(L, narg, LUA_TTABLE);
  len = (uint32_t)lua_objlen(L, narg);
  /* check if there are options in the list at position narg that are not in the array Options */
  for(i=1;i<=len;i++)
  {
    lua_rawgeti(L, narg, i);
    if(!(lua_istable(L, -1)))
    {
      luaL_error(L, "Error in option %d, Options should be a LIST of two elements of the form {\"name\",value}\n",i);
    }
    if(!(((int)lua_objlen(L, -1))==2))
    {
      luaL_error(L, "Error in option %d, Options should be a list of TWO elements of the form {\"name\",value}\n",i);
    }
    lua_rawgeti(L, -1, 1);
    if(!(lua_type(L, -1)==LUA_TSTRING))
    {
      luaL_error(L, "Error in option %d, Options should be a list of two elements of the form {\"NAME\",value}\n",i);
    }
    tokenread = (char*) luaL_checkstring(L, -1);
    tokenfound=0;
    j=0;
    while (Options[j].name!=NULL)
    {
      if (strcmp(tokenread, Options[j].name) == 0)
      {
        tokenfound=1;
        lua_rawgeti(L, -2, 2);
        
        if (strcmp(Options[j].type, "char*") == 0)
        {
          *(char**)Options[j].variable=(char*)LuaReturnVariableType(L, -1, Options[j].type);
        }
        else if(strcmp(Options[j].type, "int") == 0)
        {
          *(int*)Options[j].variable=*(int*)LuaReturnVariableType(L, -1, Options[j].type);
        }
        else if(strcmp(Options[j].type, "uint32_t") == 0)
        {
          *(uint32_t*)Options[j].variable=*(uint32_t*)LuaReturnVariableType(L, -1, Options[j].type);
        }
        else if(strcmp(Options[j].type, "uint16_t") == 0)
        {
          *(uint16_t*)Options[j].variable=*(uint16_t*)LuaReturnVariableType(L, -1, Options[j].type);
        }
        else if(strcmp(Options[j].type, "uint8_t") == 0)
        {
          *(uint32_t*)Options[j].variable=*(uint32_t*)LuaReturnVariableType(L, -1, Options[j].type);
        }
        else if(strcmp(Options[j].type, "double") == 0)
        {
          *(double*)Options[j].variable=*(double*)LuaReturnVariableType(L, -1, Options[j].type);
        }
        else if(strcmp(Options[j].type, "bool") == 0)
        {
          *(Bool*)Options[j].variable=*(Bool*)LuaReturnVariableType(L, -1, Options[j].type);
        }
        else if(strcmp(Options[j].type, "Vector") ==0)
        {
          if(((CompactVectorType*)Options[j].variable)->Ni!=0)
          {
            if(((CompactVectorType*)Options[j].variable)->QComplex)
            {
              free(((CompactVectorType*)Options[j].variable)->V);
            }
            else
            {
              free(((CompactVectorType*)Options[j].variable)->VRe);
            }
          }
          *(CompactVectorType*)Options[j].variable=*(CompactVectorType*)LuaReturnVariableType(L, -1, Options[j].type);
        }
        else if(strcmp(Options[j].type, "DeterminantRestrictionType") == 0)
        {
          *(DeterminantRestrictionType*)Options[j].variable=*(DeterminantRestrictionType*)LuaReturnVariableType(L, -1, Options[j].type);
        }
        else
        {
          printf("ERROR: unrecognized type in LuaCheckOptions [%s]\n",Options[j].type);
        }
        lua_pop(L,1);
      }
      j++;
    }
    if(!tokenfound)
    {
      Buffer = malloc(4048);
      if(Buffer==NULL)
      {
        printf("ERROR: malloc failed in LuaCheckOptions\n");
        return 1;
      }
      Buffer[0]='\0';
      j=0;
      while (Options[j].name!=NULL)
      {
        sprintfAdd(&Buffer, "%s [%s]\n",Options[j].name,Options[j].type);
        j++;
      }
      luaL_error(L, "Error in option %d, Unrecognized option \"%s\"\n Possible values for options are:\n%s",i,tokenread,Buffer);
    }
    lua_pop(L,1);
    lua_pop(L,1);
  }
  return 0;
}

