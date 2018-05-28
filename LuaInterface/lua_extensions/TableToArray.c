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

char** luaL_checkarray_string(lua_State *L, int narg, uint32_t *len_out)
{
  int i;
  luaL_checktype(L, narg, LUA_TTABLE);
  uint32_t len = (uint32_t)lua_objlen(L, narg);
  *len_out = len;
  char**buff = (char**)malloc(len*sizeof(char*));
  for(i=0;i<len;i++)
  {
    lua_rawgeti(L, narg, i+1);
    buff[i] = (char*)luaL_checkstring(L,-1);
    lua_pop(L, 1);
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
  uint32_t len;
  // Matrix stored as pointer to structure on lua stack
  if (luaL_testudata(L, narg, "Matrix_Type") != NULL)
  {
    buff = (CompactMatrixType*)luaL_checkudata(L, narg, "Matrix_Type");
    if (buff!=NULL)
    {
      return buff;
    }
  }
  luaL_checktype(L, narg, LUA_TTABLE);
  len = (uint32_t)lua_objlen(L, narg);
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
  buff->a0=0;
  buff->b0=0;
  buff->mu=0;
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

/*  Check BlockGreensfunction reads a BlockGreensfunction from the stack
 *  it calls the one of the functions
 *     check (Block) TriDiagonalMatrixType
 *     check (Block) AndersonMatrixType
 *     check (Block) NaturalImpurityMatrixType
 *     check (Block) ListOfPoles
 *
 *  the later functions read from the lua stack a function of any of the types possible
 *  and transform this to the type requested
 *  If SuggestedGrid is set this will be used for the type conversion
 *  Note that the sugested grid is only used when needed, i.e. The Green's funciton can well be on a different Grid then sugested
 *
 *  The format of the different types on the lua stack are:
 *
 *  The format of TriDiagonalMatrixType should be
 *    { {a0,a1,a2,a3,...,an}, {b0,b1,b2,...,bn-1}, mu=0, type="Tri"}
 *    and defines the Green's function
 *    G0 = a0 + b0^2 . ( G1 )
 *       with Gi = 1/(w - ai - bi^2 ( Gi+1 ) )
 *       and Gn = 1/(w - an )
 *       and ai and bi real
 *
 *  The format of BlockTriDiagonalMatrixType should be
 *    { {A0,A1,A2,A3,..,An}, {B0,B1,B2,...,Bn-1}, mu=0, type="Tri"}
 *    and defines the Green's function
 *    G0 = A0 + B0^* . ( G1 ) . B0^T
 *       with Gi = 1/(w - Ai - Bi . ( Gi+1 ) . Bi^+)
 *       and Gn = 1/(w - An)
 *
 *  The format of AndersonMatrixType should be
 *    { {a0,a1,a2,a3,...,an}, {b0,b1,b2,...,bn-1}, mu=0, type="And"}
 *    and defines the Green's function
 *    G = a0 + b0^2 . ( 1/(w - a1 - sum_{i=1}^{n-1} bi^2 1/(w-ai+1) ) )
 *       and an and bn real
 *
 *  The format of BlockAndersonMatrixType should be
 *    { {A0,A1,A2,A3,...,An}, {B0,B1,B2,...,Bn-1}, mu=0, type="And"}
 *    and defines the Green's function
 *    G = A0 + B0^* . ( 1/(w - A1 - sum_{i=1}^{n-1} Bi . 1/(w-Ai+1) . Bi^+ ) ) . B0^T
 *
 *
 *  The format of NaturalImpurityMatrixType should be
 *    { {a0, b0}, val=Gval, con=Gcon, mu=0, type="Nat" }
 *    and defines the Green's function
 *    a0 + b0^2 ( Gval(w) + Gcon(w) )
 *
 *  The format of BlockNaturalImpurityMatrixType should be
 *    { {A0, B0}, val=Gval, con=Gcon, mu=0, type="Nat" }
 *    and defines the Green's function
 *    a0 + b0^2 ( Gval(w) + Gcon(w) )
 *
 *
 *  The format of ListOfPoles should be:
 *    { {a0,a1,a2,a3,..,an}, {b1,b2,b3,..,bn}, mu=0, type="ListOfPoles"}
 *    and defines the Green's function
 *    a0 + sum_{i=1^N} b_i * 1/(w-a_i)
 *
 *  The format of BlockListOfPoles should be:
 *    { {A0,a1,a2,a3,..,an}, {B1,B2,B3,..,Bn}, mu=0, type="ListOfPoles"}
 *    and defines the Green's function
 *    A0 + sum_{i=1^N} B_i * 1/(w-a_i)
 */

/*
 *  The format of TriDiagonalMatrixType should be
 *    { {a0,a1,a2,a3,...,an}, {b0,b1,b2,...,an-1}, mu=0, type="Tri"}
 *    and defines the Green's function
 *    G0 = a0 + b0^2 . ( G1 )
 *       with Gi = 1/(w - ai - bi^2 ( Gi+1 ) )
 *       stopping at  Gn = 1/(w - an )
 *       and ai and bi real
 */
TriDiagonalMatrixType* luaL_checkTriDiagonalMatrixType(lua_State *L, int narg)
{
  const char *type;
  uint32_t len;
  TriDiagonalMatrixType *buff;
  if(lua_istable(L, narg)==0)
  {
    luaL_error(L, "Error Green's function should be represented by a table");
  }
  lua_getfield(L, narg, "type");
  if(lua_isstring(L,-1)==0)
  {
    luaL_error(L, "Error Reading Green's function without type declaration.");
  }
  type = lua_tolstring(L, -1, NULL);
  if(strcmp(type, "Tri")==0)
  {
    lua_pop(L, 1);
    buff = malloc(sizeof(TriDiagonalMatrixType));
    if(lua_objlen(L,narg)!=2)
    {
      luaL_error(L, "Error a Green's function represented as a tri-diagonal matrix should be a list of two lists {{a0,a1,a2,...,an},{b0,b1,b2,...,bn-1}}");
    }
    lua_getfield(L, narg, "mu");
    if(lua_isnumber(L, -1)==1)
    {
      buff->mu=(double)luaL_checknumber(L, -1);
    }
    else
    {
      luaL_error(L, "Reading Green's function represented as a tri-diagonal matrix without mu defined");
    }
    lua_pop(L,1);
    lua_getfield(L, narg, "name");
    if(lua_isstring(L, -1)==1)
    {
      strncpy(buff->Name, luaL_checkstring(L, -1), MAXStrLength);
      buff->Name[MAXStrLength-1]='\0';
    }
    else if(lua_isnil(L, -1)==1)
    {
      buff->Name[0]='\0';
    }
    else
    {
      luaL_error(L, "Reading Green's function represented as a tri-diagonal matrix without name defined");
    }
    lua_pop(L,1);
    // read a
    lua_rawgeti(L, narg, 1);
    buff->a=luaL_checkarray_double(L, -1, &len);
    if(len==0) luaL_error(L, "Errorr reading Greens function with lenght of a equal to zero");
    buff->N=len-1;
    lua_pop(L, 1);
    // read b
    lua_rawgeti(L, narg, 2);
    buff->b=luaL_checkarray_double(L, -1, &len);
    if(len!=buff->N)
    {
      luaL_error(L, "Error Reading List of poles. Length of a (%d) should be one more than length of b (%d)",buff->N+1,len);
    }
    lua_pop(L, 1);
    return buff;
  }
  if(strcmp(type, "And")==0)
  {
    AndersonMatrixType *temp;
    lua_pop(L, 1);
    temp = luaL_checkAndersonMatrixType(L, narg);
    buff = malloc(sizeof(TriDiagonalMatrixType));
    AndersonMatrixToTriDiagonalMatrix(*temp, buff,temp->N);
    FreeAndersonMatrix(temp);
    free(temp);
    return buff;
  }
  if((strcmp(type, "NaturalImpurityOrbital")==0)||(strcmp(type, "Nat")==0))
  {
    BADoubleTriDiagonalMatrixType *temp;
    lua_pop(L, 1);
    temp = luaL_checkNatrualImpurityOrbitalMatrixType(L, narg);
    buff = malloc(sizeof(TriDiagonalMatrixType));
    BADoubleTriDiagonalMatrixToTriDiagonalMatrix(*temp, buff);
    FreeBADoubleTriDiagonalMatrix(temp);
    free(temp);
    return buff;
  }
  if(strcmp(type, "ListOfPoles")==0)
  {
    ListOfPoles *temp;
    lua_pop(L, 1);
    temp = luaL_checkListOfPoles(L, narg);
    buff = malloc(sizeof(TriDiagonalMatrixType));
    ListOfPolesToTridiagonalMatrix(*temp, buff,temp->N);
    FreeListOfPoles(temp);
    free(temp);
    return buff;
  }
  luaL_error(L, "Errorr reading Greens function of unknown type (%s)",type);
  return NULL;
}

/*
 *  The format of BlockTriDiagonalMatrixType should be
 *    { {A0,A1,A2,A3,...An}, {B0,B1,B2,...,Bn-1}, mu=0, type="Tri"}
 *    and defines the Green's function
 *    G0 = A0 + B0^* . ( G1 ) . B0^T
 *       with Gi = 1/(w - Ai - Bi . ( Gi+1 ) . Bi^+)
 *       with Gn = 1/(w - An)
 */
BlockTriDiagonalMatrixType* luaL_checkBlockTriDiagonalMatrixType(lua_State *L, int narg)
{
  const char *type;
  uint32_t len,i,j,k;
  Bool MakeReal=True;
  BlockTriDiagonalMatrixType *buff;
  if(lua_istable(L, narg)==0)
  {
    luaL_error(L, "Error Green's function should be represented by a table");
  }
  lua_getfield(L, narg, "type");
  if(lua_isstring(L,-1)==0)
  {
    luaL_error(L, "Error Reading Green's function without type declaration.");
  }
  type = lua_tolstring(L, -1, NULL);
  if(strcmp(type, "Tri")==0)
  {
    lua_pop(L, 1);
    buff = malloc(sizeof(BlockTriDiagonalMatrixType));
    if(lua_objlen(L,narg)!=2)
    {
      luaL_error(L, "Error a tridiagonal matrix should be a list of two lists {{A0,A1,...,An+1},{B0,B1,...,Bn}}, with A and B matrices and A Hermitian");
    }
    // read b
    lua_rawgeti(L, narg, 2);
    {
      len = (uint32_t)lua_objlen(L, -1);
      if(len==0) luaL_error(L, "Errorr reading Greens function with lenght of b equal to zero");
      buff->N=len;
      lua_rawgeti(L, -1, 1);
      buff->BlockSize = (uint32_t)lua_objlen(L, -1);
      lua_pop(L, 1);
      buff->QComplex=1;
      InitBlockTridiagonalMatrix(buff);
      for(i=0;i<buff->N;i++)
      {
        lua_rawgeti(L, -1, i+1); // get b_i onto the stack
        {
          if((uint32_t)lua_objlen(L, -1)!=buff->BlockSize)
          {
            luaL_error(L, "Errorr reading Greens function as block tridiagonal matrix, but all blocks must be equal in size. BlockSize (%d), size at i=%d is %d",buff->BlockSize,i+1,lua_objlen(L, -1));
          }
          for(j=0;j<buff->BlockSize;j++)
          {
            lua_rawgeti(L,-1,j+1);
            {
              if((uint32_t)lua_objlen(L, -1)!=buff->BlockSize)
              {
                luaL_error(L, "Errorr reading Greens function as block tridiagonal matrix, but all blocks must be equal in size. BlockSize (%d), size at i=%d j=%d is %d",buff->BlockSize,i+1,j+1,lua_objlen(L, -1));
              }
              for(k=0;k<buff->BlockSize;k++)
              {
                lua_rawgeti(L,-1,k+1);
                {
                  if(lua_isnumber(L, -1))
                  {
                    buff->b[i*(Sqr(buff->BlockSize))+k+j*buff->BlockSize].Re=(double)luaL_checknumber(L, -1);
                    buff->b[i*(Sqr(buff->BlockSize))+k+j*buff->BlockSize].Im=0;
                  }
                  else if (luaL_testudata(L, -1, "Complex_Type") != NULL)
                  {
                    MakeReal=False;
                    buff->b[i*(Sqr(buff->BlockSize))+k+j*buff->BlockSize] = *(Complex16 *)luaL_checkudata(L, -1, "Complex_Type");
                  }
                  else
                  {
                    luaL_error(L, "Error in reading complex table\n");
                  }
                }
                lua_pop(L, 1);
              }
            }
            lua_pop(L, 1);
          }
        }
        lua_pop(L, 1);
      }
    }
    lua_pop(L, 1);
    // read a
    lua_rawgeti(L, narg, 1);
    {
      len = (uint32_t)lua_objlen(L, -1);
      if(len!=buff->N+1) luaL_error(L, "Errorr reading Greens function in form of Block list of poles length of A(%d) should be one larger then length of B(%d)",len,buff->N);
      for(i=0;i<=buff->N;i++)
      {
        lua_rawgeti(L, -1, i+1); // get a_i onto the stack
        {
          if((uint32_t)lua_objlen(L, -1)!=buff->BlockSize)
          {
            luaL_error(L, "Errorr reading Greens function as block tridiagonal matrix, but all blocks must be equal in size. BlockSize (%d), size at i=%d is %d",buff->BlockSize,i+1,lua_objlen(L, -1));
          }
          for(j=0;j<buff->BlockSize;j++)
          {
            lua_rawgeti(L,-1,j+1);
            {
              if((uint32_t)lua_objlen(L, -1)!=buff->BlockSize)
              {
                luaL_error(L, "Errorr reading Greens function as block tridiagonal matrix, but all blocks must be equal in size. BlockSize (%d), size at i=%d j=%d is %d",buff->BlockSize,i+1,j+1,lua_objlen(L, -1));
              }
              for(k=0;k<buff->BlockSize;k++)
              {
                lua_rawgeti(L,-1,k+1);
                {
                  if(lua_isnumber(L, -1))
                  {
                    buff->a[i*(Sqr(buff->BlockSize))+k+j*buff->BlockSize].Re=(double)luaL_checknumber(L, -1);
                    buff->a[i*(Sqr(buff->BlockSize))+k+j*buff->BlockSize].Im=0;
                  }
                  else if (luaL_testudata(L, -1, "Complex_Type") != NULL)
                  {
                    MakeReal=False;
                    buff->a[i*(Sqr(buff->BlockSize))+k+j*buff->BlockSize] = *(Complex16 *)luaL_checkudata(L, -1, "Complex_Type");
                  }
                  else
                  {
                    luaL_error(L, "Error in reading complex table\n");
                  }
                }
                lua_pop(L, 1);
              }
            }
            lua_pop(L, 1);
          }
        }
        lua_pop(L, 1);
      }
    }
    lua_pop(L, 1);
    if(MakeReal==True)
    {
      ComplexBlockTridiagonalMatrixToReal(buff);
    }
    lua_getfield(L, narg, "mu");
    if(lua_isnumber(L, -1)==1)
    {
      buff->mu=(double)luaL_checknumber(L, -1);
    }
    else
    {
      luaL_error(L, "Reading Block tridiagonal matrix without mu defined");
    }
    lua_pop(L,1);
    lua_getfield(L, narg, "name");
    if(lua_isstring(L, -1)==1)
    {
      strncpy(buff->Name, luaL_checkstring(L, -1), MAXStrLength);
      buff->Name[MAXStrLength-1]='\0';
    }
    else if(lua_isnil(L, -1)==1)
    {
      buff->Name[0]='\0';
    }
    else
    {
      luaL_error(L, "Reading Block tridiagonal matrix without name correctly defined");
    }
    lua_pop(L,1);
    return buff;
  }
  if(strcmp(type, "And")==0)
  {
    BlockAndersonMatrixType *temp;
    lua_pop(L, 1);
    temp = luaL_checkBlockAndersonMatrixType(L, narg);
    buff = malloc(sizeof(BlockTriDiagonalMatrixType));
    BlockAndersonMatrixToBlockTridiagonalMatrix(*temp, buff);
    FreeBlockAndersonMatrix(temp);
    free(temp);
    return buff;
  }
  if((strcmp(type, "NaturalImpurityOrbital")==0)||(strcmp(type, "Nat")==0))
  {
    BlockNaturalImpurityOrbitalMatrixType *temp;
    lua_pop(L, 1);
    temp = luaL_checkBlockNatrualImpurityOrbitalMatrixType(L, narg);
    buff = malloc(sizeof(BlockTriDiagonalMatrixType));
    BlockNaturalImpurityOrbitalMatrixToBlockTridiagonalMatrix(*temp, buff);
    FreeBlockNaturalImpurityOrbitalMatrix(temp);
    free(temp);
    return buff;
  }
  if(strcmp(type, "ListOfPoles")==0)
  {
    BlockListOfPoles *temp;
    lua_pop(L, 1);
    temp = luaL_checkBlockListOfPoles(L, narg);
    buff = malloc(sizeof(BlockTriDiagonalMatrixType));
    BlockListOfPolesToBlockTridiagonalMatrix(*temp, buff, temp->N, 100*DBL_EPSILON);
    FreeBlockListOfPoles(temp);
    free(temp);
    return buff;
  }
  luaL_error(L, "Errorr reading Greens function of unknown type (%s)",type);
  return NULL;
}

/*
 *  The format of AndersonMatrixType should be
 *    { {a0,a1,a2,a3,...,an}, {b0,b1,b2,...,bn-1}, mu=0, type="And"}
 *    and defines the Green's function
 *    G0 = a0 + b0^2 . ( G1 )
 *       with Gn = 1/(w - an - bn^2 ( Gn+1 ) )
 *       and an and bn real
 */
AndersonMatrixType* luaL_checkAndersonMatrixType(lua_State *L, int narg)
{
  const char *type;
  uint32_t len;
  AndersonMatrixType *buff;
  if(lua_istable(L, narg)==0)
  {
    luaL_error(L, "Error Green's function should be represented by a table");
  }
  lua_getfield(L, narg, "type");
  if(lua_isstring(L,-1)==0)
  {
    luaL_error(L, "Error Reading Green's function without type declaration.");
  }
  type = lua_tolstring(L, -1, NULL);
  if(strcmp(type, "Tri")==0)
  {
    TriDiagonalMatrixType *temp;
    lua_pop(L, 1);
    temp = luaL_checkTriDiagonalMatrixType(L, narg);
    buff = malloc(sizeof(AndersonMatrixType));
    TriDiagonalMatrixToAndersonMatrix(*temp, buff);
    FreeTriDiagonalMatrix(temp);
    free(temp);
    return buff;
  }
  if(strcmp(type, "And")==0)
  {
    lua_pop(L, 1);
    buff = malloc(sizeof(AndersonMatrixType));
    if(lua_objlen(L,narg)!=2)
    {
      luaL_error(L, "Error a Green's function represented as a Anderson matrix should be a list of two lists {{a0,a1,a2,...,an},{b0,b1,b2,...,bn-1}}");
    }
    lua_getfield(L, narg, "mu");
    if(lua_isnumber(L, -1)==1)
    {
      buff->mu=(double)luaL_checknumber(L, -1);
    }
    else
    {
      luaL_error(L, "Reading Green's function represented as an Anderson matrix without mu defined");
    }
    lua_pop(L,1);
    lua_getfield(L, narg, "name");
    if(lua_isstring(L, -1)==1)
    {
      strncpy(buff->Name, luaL_checkstring(L, -1), MAXStrLength);
      buff->Name[MAXStrLength-1]='\0';
    }
    else if(lua_isnil(L, -1)==1)
    {
      buff->Name[0]='\0';
    }
    else
    {
      luaL_error(L, "Reading Green's function represented as an Anderson matrix without name defined");
    }
    lua_pop(L,1);
    // read a
    lua_rawgeti(L, narg, 1);
    buff->a=luaL_checkarray_double(L, -1, &len);
    if(len==0) luaL_error(L, "Errorr reading Greens function with lenght of a equal to zero");
    buff->N=len-1;
    lua_pop(L, 1);
    // read b
    lua_rawgeti(L, narg, 2);
    buff->b=luaL_checkarray_double(L, -1, &len);
    if(len!=buff->N)
    {
      luaL_error(L, "Error Reading Andeson matrix. Length of a (%d) should be one more than length of b (%d)",buff->N+1,len);
    }
    lua_pop(L, 1);
    return buff;
  }
  if((strcmp(type, "NaturalImpurityOrbital")==0)||(strcmp(type, "Nat")==0))
  {
    BADoubleTriDiagonalMatrixType *temp;
    lua_pop(L, 1);
    temp = luaL_checkNatrualImpurityOrbitalMatrixType(L, narg);
    buff = malloc(sizeof(AndersonMatrixType));
    BADoubleTriDiagonalMatrixToAndersonMatrix(*temp, buff);
    FreeBADoubleTriDiagonalMatrix(temp);
    free(temp);
    return buff;
  }
  if(strcmp(type, "ListOfPoles")==0)
  {
    ListOfPoles *temp;
    lua_pop(L, 1);
    temp = luaL_checkListOfPoles(L, narg);
    buff = malloc(sizeof(AndersonMatrixType));
    ListOfPolesToAndersonMatrix(*temp, buff);
    FreeListOfPoles(temp);
    free(temp);
    return buff;
  }
  luaL_error(L, "Errorr reading Greens function of unknown type (%s)",type);
  return NULL;
}

/*
 *  The format of BlockAndersonMatrixType should be
 *    { {A0,A1,A2,A3,...,An}, {B0,B1,B2,...,Bn-1}, mu=0, type="And"}
 *    and defines the Green's function
 *    G0 = A0 + B0^* . ( G1 ) . B0^T
 *       with Gn = 1/(w - An - Bn . ( Gn+1 ) . Bn^+ )
 */
BlockAndersonMatrixType* luaL_checkBlockAndersonMatrixType(lua_State *L, int narg)
{
  const char *type;
  uint32_t len,i,j,k;
  Bool MakeReal=True;
  BlockAndersonMatrixType *buff;
  if(lua_istable(L, narg)==0)
  {
    luaL_error(L, "Error Green's function should be represented by a table");
  }
  lua_getfield(L, narg, "type");
  if(lua_isstring(L,-1)==0)
  {
    luaL_error(L, "Error Reading Green's function without type declaration.");
  }
  type = lua_tolstring(L, -1, NULL);
  if(strcmp(type, "Tri")==0)
  {
    BlockTriDiagonalMatrixType *temp;
    lua_pop(L, 1);
    temp = luaL_checkBlockTriDiagonalMatrixType(L, narg);
    buff = malloc(sizeof(BlockAndersonMatrixType));
    BlockTriDiagonalMatrixToBlockAndersonMatrix(*temp, buff);
    FreeBlockTridiagonalMatrix(temp);
    free(temp);
    return buff;
  }
  if(strcmp(type, "And")==0)
  {
    lua_pop(L, 1);
    buff = malloc(sizeof(BlockAndersonMatrixType));
    if(lua_objlen(L,narg)!=2)
    {
      luaL_error(L, "Error an Anderson matrix should be a list of two lists {{A0,A1,...,An+1},{B0,B1,...,Bn}}, with A and B matrices and A Hermitian");
    }
    // read b
    lua_rawgeti(L, narg, 2);
    {
      len = (uint32_t)lua_objlen(L, -1);
      if(len==0) luaL_error(L, "Errorr reading Greens function with lenght of b equal to zero");
      buff->N=len;
      lua_rawgeti(L, -1, 1);
      buff->BlockSize = (uint32_t)lua_objlen(L, -1);
      lua_pop(L, 1);
      buff->QComplex=1;
      InitBlockAndersonMatrix(buff);
      for(i=0;i<buff->N;i++)
      {
        lua_rawgeti(L, -1, i+1); // get b_i onto the stack
        {
          if((uint32_t)lua_objlen(L, -1)!=buff->BlockSize)
          {
            luaL_error(L, "Errorr reading Greens function as block Anderson matrix, but all blocks must be equal in size. BlockSize (%d), size at i=%d is %d",buff->BlockSize,i+1,lua_objlen(L, -1));
          }
          for(j=0;j<buff->BlockSize;j++)
          {
            lua_rawgeti(L,-1,j+1);
            {
              if((uint32_t)lua_objlen(L, -1)!=buff->BlockSize)
              {
                luaL_error(L, "Errorr reading Greens function as block Anderson matrix, but all blocks must be equal in size. BlockSize (%d), size at i=%d j=%d is %d",buff->BlockSize,i+1,j+1,lua_objlen(L, -1));
              }
              for(k=0;k<buff->BlockSize;k++)
              {
                lua_rawgeti(L,-1,k+1);
                {
                  if(lua_isnumber(L, -1))
                  {
                    buff->b[i*(Sqr(buff->BlockSize))+k+j*buff->BlockSize].Re=(double)luaL_checknumber(L, -1);
                    buff->b[i*(Sqr(buff->BlockSize))+k+j*buff->BlockSize].Im=0;
                  }
                  else if (luaL_testudata(L, -1, "Complex_Type") != NULL)
                  {
                    MakeReal=False;
                    buff->b[i*(Sqr(buff->BlockSize))+k+j*buff->BlockSize] = *(Complex16 *)luaL_checkudata(L, -1, "Complex_Type");
                  }
                  else
                  {
                    luaL_error(L, "Error in reading complex table\n");
                  }
                }
                lua_pop(L, 1);
              }
            }
            lua_pop(L, 1);
          }
        }
        lua_pop(L, 1);
      }
    }
    lua_pop(L, 1);
    // read a
    lua_rawgeti(L, narg, 1);
    {
      len = (uint32_t)lua_objlen(L, -1);
      if(len!=buff->N+1) luaL_error(L, "Errorr reading Greens function in form of Block Anderson Matrix length of A(%d) should be one larger then length of B(%d)",len,buff->N);
      for(i=0;i<=buff->N;i++)
      {
        lua_rawgeti(L, -1, i+1); // get a_i onto the stack
        {
          if((uint32_t)lua_objlen(L, -1)!=buff->BlockSize)
          {
            luaL_error(L, "Errorr reading Greens function as block Anderson matrix, but all blocks must be equal in size. BlockSize (%d), size at i=%d is %d",buff->BlockSize,i+1,lua_objlen(L, -1));
          }
          for(j=0;j<buff->BlockSize;j++)
          {
            lua_rawgeti(L,-1,j+1);
            {
              if((uint32_t)lua_objlen(L, -1)!=buff->BlockSize)
              {
                luaL_error(L, "Errorr reading Greens function as block Anderson matrix, but all blocks must be equal in size. BlockSize (%d), size at i=%d j=%d is %d",buff->BlockSize,i+1,j+1,lua_objlen(L, -1));
              }
              for(k=0;k<buff->BlockSize;k++)
              {
                lua_rawgeti(L,-1,k+1);
                {
                  if(lua_isnumber(L, -1))
                  {
                    buff->a[i*(Sqr(buff->BlockSize))+k+j*buff->BlockSize].Re=(double)luaL_checknumber(L, -1);
                    buff->a[i*(Sqr(buff->BlockSize))+k+j*buff->BlockSize].Im=0;
                  }
                  else if (luaL_testudata(L, -1, "Complex_Type") != NULL)
                  {
                    MakeReal=False;
                    buff->a[i*(Sqr(buff->BlockSize))+k+j*buff->BlockSize] = *(Complex16 *)luaL_checkudata(L, -1, "Complex_Type");
                  }
                  else
                  {
                    luaL_error(L, "Error in reading complex table\n");
                  }
                }
                lua_pop(L, 1);
              }
            }
            lua_pop(L, 1);
          }
        }
        lua_pop(L, 1);
      }
    }
    lua_pop(L, 1);
    if(MakeReal==True)
    {
      ComplexBlockAndersonMatrixToReal(buff);
    }
    lua_getfield(L, narg, "mu");
    if(lua_isnumber(L, -1)==1)
    {
      buff->mu=(double)luaL_checknumber(L, -1);
    }
    else
    {
      luaL_error(L, "Reading Block Anderson matrix without mu defined");
    }
    lua_pop(L,1);
    lua_getfield(L, narg, "name");
    if(lua_isstring(L, -1)==1)
    {
      strncpy(buff->Name, luaL_checkstring(L, -1), MAXStrLength);
      buff->Name[MAXStrLength-1]='\0';
    }
    else if(lua_isnil(L, -1)==1)
    {
      buff->Name[0]='\0';
    }
    else
    {
      luaL_error(L, "Reading Block Anderson matrix without name correctly defined");
    }
    lua_pop(L,1);
    return buff;
  }
  if((strcmp(type, "NaturalImpurityOrbital")==0)||(strcmp(type, "Nat")==0))
  {
    BlockNaturalImpurityOrbitalMatrixType *temp;
    lua_pop(L, 1);
    temp = luaL_checkBlockNatrualImpurityOrbitalMatrixType(L, narg);
    buff = malloc(sizeof(BlockAndersonMatrixType));
    BlockNaturalImpurityOrbitalMatrixToBlockAndersonMatrix(*temp, buff);
    FreeBlockNaturalImpurityOrbitalMatrix(temp);
    free(temp);
    return buff;
  }
  if(strcmp(type, "ListOfPoles")==0)
  {
    BlockListOfPoles *temp;
    lua_pop(L, 1);
    temp = luaL_checkBlockListOfPoles(L, narg);
    buff = malloc(sizeof(BlockAndersonMatrixType));
    BlockListOfPolesToBlockAndersonMatrix(*temp, buff, temp->N);
    FreeBlockListOfPoles(temp);
    free(temp);
    return buff;
  }
  luaL_error(L, "Errorr reading Greens function of unknown type (%s)",type);
  return NULL;
}

/*
 *  The format of NaturalImpurityMatrixType should be
 *    { {a0, b0}, val=Gval, con=Gcon, mu=0, type="Nat" }
 *    and defines the Green's function
 *    a0 + b0^2 ( Gval(w) + Gcon(w) )
 */
BADoubleTriDiagonalMatrixType * luaL_checkNatrualImpurityOrbitalMatrixType(lua_State *L, int narg)
{
  const char *type;
  BADoubleTriDiagonalMatrixType *buff;
  if(lua_istable(L, narg)==0)
  {
    luaL_error(L, "Error Green's function should be represented by a table");
  }
  lua_getfield(L, narg, "type");
  if(lua_isstring(L,-1)==0)
  {
    luaL_error(L, "Error Reading Green's function without type declaration.");
  }
  type = lua_tolstring(L, -1, NULL);
  if(strcmp(type, "Tri")==0)
  {
    TriDiagonalMatrixType *temp;
    lua_pop(L, 1);
    temp = luaL_checkTriDiagonalMatrixType(L, narg);
    buff = malloc(sizeof(BADoubleTriDiagonalMatrixType));
    TriDiagonalMatrixToBADoubleTriDiagonalMatrix(*temp, buff, temp->N, temp->N, temp->mu, 1000*DBL_MIN);
    FreeTriDiagonalMatrix(temp);
    free(temp);
    return buff;
  }
  if(strcmp(type, "And")==0)
  {
    AndersonMatrixType *temp;
    lua_pop(L, 1);
    temp = luaL_checkAndersonMatrixType(L, narg);
    buff = malloc(sizeof(BADoubleTriDiagonalMatrixType));
    AndersonMatrixToBADoubleTriDiagonalMatrixSimple(*temp, buff);
    FreeAndersonMatrix(temp);
    free(temp);
    return buff;
  }
  if((strcmp(type, "NaturalImpurityOrbital")==0)||(strcmp(type, "Nat")==0))
  {
    TriDiagonalMatrixType *temp;
    buff = malloc(sizeof(BADoubleTriDiagonalMatrixType));
    lua_pop(L,1);
    /* read val */
    lua_getfield(L, narg, "val");
    if(lua_istable(L,-1)==0)
    {
      luaL_error(L, "Error Green's function of format NaturalImpurityOrbital should have a field val");
    }
    temp = luaL_checkTriDiagonalMatrixType(L, lua_gettop(L));
    buff->val = *temp;
    free(temp);
    lua_pop(L,1);
    /* read con */
    lua_getfield(L, narg, "con");
    if(lua_istable(L,-1)==0)
    {
      luaL_error(L, "Error Green's function of format NaturalImpurityOrbital should have a field con");
    }
    temp = luaL_checkTriDiagonalMatrixType(L, lua_gettop(L));
    buff->con = *temp;
    free(temp);
    lua_pop(L,1);
    /* read {a0,b0} */
    lua_rawgeti(L, narg, 1);
    if(lua_istable(L,-1)==0)
    {
      luaL_error(L, "Error Green's function of format NaturalImpurityOrbital should have the table {a0, b0} at position 1");
    }
    lua_rawgeti(L, -1, 1);
    buff->a0 = (double)luaL_checknumber(L, -1);
    lua_pop(L,1);
    lua_rawgeti(L, -1, 2);
    buff->b0 = (double)luaL_checknumber(L, -1);
    lua_pop(L,2);
    /* read mu */
    lua_getfield(L, narg, "mu");
    if(lua_isnumber(L, -1)==1)
    {
      buff->mu=(double)luaL_checknumber(L, -1);
    }
    else
    {
      luaL_error(L, "Reading NaturalImpurityOrbital matrix without mu defined");
    }
    lua_pop(L,1);
    /* read name */
    lua_getfield(L, narg, "name");
    if(lua_isstring(L, -1)==1)
    {
      strncpy(buff->Name, luaL_checkstring(L, -1), MAXStrLength);
      buff->Name[MAXStrLength-1]='\0';
    }
    else if(lua_isnil(L, -1)==1)
    {
      buff->Name[0]='\0';
    }
    else
    {
      luaL_error(L, "Reading NaturalImpurityOrbital matrix without name correctly defined");
    }
    lua_pop(L,1);
    /* read epsilon */
    lua_getfield(L, narg, "epsilon");
    if(lua_isnumber(L, -1)==1)
    {
      buff->Epsilon=(double)luaL_checknumber(L, -1);
    }
    else
    {
      buff->Epsilon=1000*DBL_EPSILON;
    }
    lua_pop(L,1);
    buff->acvB  = Sqr(buff->val.b[0]) * buff->val.a[1] + Sqr(buff->con.b[0]) * buff->con.a[1];
    buff->acvA  = Sqr(buff->val.b[0]) * buff->con.a[1] + Sqr(buff->con.b[0]) * buff->val.a[1];
    buff->acvBA =  buff->val.b[0] * buff->con.b[0] * (buff->val.a[1]-buff->con.a[1]);
    buff->bvB0  =  buff->val.b[0] * buff->val.b[1];
    buff->bcB0  =  buff->con.b[0] * buff->con.b[1];
    buff->bvA0  =  buff->con.b[0] * buff->val.b[1];
    buff->bcA0  = -buff->val.b[0] * buff->con.b[1];
    return buff;
  }
  if(strcmp(type, "ListOfPoles")==0)
  {
    ListOfPoles *temp;
    lua_pop(L, 1);
    temp = luaL_checkListOfPoles(L, narg);
    buff = malloc(sizeof(BADoubleTriDiagonalMatrixType));
    ListOfPolesToBADoubleTridiagonalMatrix(*temp, buff,temp->N,temp->N,temp->mu,1000*DBL_MIN);
    FreeListOfPoles(temp);
    free(temp);
    return buff;
  }
  luaL_error(L, "Errorr reading Greens function of unknown type (%s)",type);
  return NULL;
}

/*
 *
 */
BlockNaturalImpurityOrbitalMatrixType * luaL_checkBlockNatrualImpurityOrbitalMatrixType(lua_State *L, int narg)
{
  const char *type;
  BlockNaturalImpurityOrbitalMatrixType *buff;
  if(lua_istable(L, narg)==0)
  {
    luaL_error(L, "Error Green's function should be represented by a table");
  }
  lua_getfield(L, narg, "type");
  if(lua_isstring(L,-1)==0)
  {
    luaL_error(L, "Error Reading Green's function without type declaration.");
  }
  type = lua_tolstring(L, -1, NULL);
  if(strcmp(type, "Tri")==0)
  {
    BlockTriDiagonalMatrixType *temp;
    lua_pop(L, 1);
    temp = luaL_checkBlockTriDiagonalMatrixType(L, narg);
    buff = malloc(sizeof(BlockNaturalImpurityOrbitalMatrixType));
    BlockTriDiagonalMatrixToBlockNaturalImpurityMatrix(*temp, buff);
    FreeBlockTridiagonalMatrix(temp);
    free(temp);
    return buff;
  }
  if(strcmp(type, "And")==0)
  {
    BlockAndersonMatrixType *temp;
    lua_pop(L, 1);
    temp = luaL_checkBlockAndersonMatrixType(L, narg);
    buff = malloc(sizeof(BlockNaturalImpurityOrbitalMatrixType));
    BlockAndersonMatrixToBlockNaturalImpurityOrbitalMatrix(*temp, buff);
    FreeBlockAndersonMatrix(temp);
    free(temp);
    return buff;
  }
  if((strcmp(type, "NaturalImpurityOrbital")==0)||(strcmp(type, "Nat")==0))
  {
    BlockTriDiagonalMatrixType *temp;
    CompactMatrixType *tempM;
    lua_pop(L, 1);
    buff = malloc(sizeof(BlockNaturalImpurityOrbitalMatrixType));
    lua_getfield(L, narg, "val");
    if(lua_istable(L,-1)==0)
    {
      luaL_error(L, "Error Green's function of format NaturalImpurityOrbital should have a field val");
    }
    temp = luaL_checkBlockTriDiagonalMatrixType(L, lua_gettop(L));
    buff->val = *temp;
    free(temp);
    lua_pop(L,1);
    lua_getfield(L, narg, "con");
    if(lua_istable(L,-1)==0)
    {
      luaL_error(L, "Error Green's function of format NaturalImpurityOrbital should have a field con");
    }
    temp = luaL_checkBlockTriDiagonalMatrixType(L, lua_gettop(L));
    buff->con = *temp;
    free(temp);
    lua_pop(L,1);
    lua_rawgeti(L, narg, 1);
    if(lua_istable(L,-1)==0)
    {
      luaL_error(L, "Error Green's function of format NaturalImpurityOrbital should have the table {A0, B0} at position 1");
    }
    lua_rawgeti(L, -1, 1);
    tempM = luaL_checkCompactMatrix(L, lua_gettop(L));
    buff->a0 = *tempM;
    free(tempM);
    lua_pop(L,1);
    // read B0 from stack. if B0 is not the identity matrix scale val and con
    lua_rawgeti(L, -1, 2);
    tempM = luaL_checkCompactMatrix(L, lua_gettop(L));
    
    
    free(tempM);
    lua_pop(L,1);
    
    
    
    
    //TODO: set matrix TMT
    lua_pop(L,1);
    return buff;
  }
  if(strcmp(type, "ListOfPoles")==0)
  {
    BlockListOfPoles *temp;
    lua_pop(L, 1);
    temp = luaL_checkBlockListOfPoles(L, narg);
    buff = malloc(sizeof(BlockNaturalImpurityOrbitalMatrixType));
    BlockListOfPolesToBlockNaturalImpurityOrbitalMatrix(temp, buff, 1000 * DBL_MIN, temp->N, temp->N);
    FreeBlockListOfPoles(temp);
    free(temp);
    return buff;
  }
  luaL_error(L, "Errorr reading Greens function of unknown type (%s)",type);
  return NULL;
}

/*
 *  On the lua stack at position narg should be a response function of any type
 *  Returned is a list of poles
 *  transformation to list of poles depends on the setting SuggestedGrid
 *  if NSuggestedGrid=0 the transformation is made exact, thereby increasing the number of Blocks
 *  For a tri-diagonal representation with N Blocks of NBlockSize there will be N*NBlockSize poles
 *  each of NBlockSize.
 *
 *  The format of ListOfPoles should be:
 
 *    { {a0,a1,a2,a3}, {b1,b2,b3}, mu=0, type="ListOfPoles"}
 *    and defines the Green's function
 *    a0 + sum_{i=1^N} b_i * 1/(w-a_i)
 */
ListOfPoles * luaL_checkListOfPoles(lua_State *L, int narg)
{
  uint32_t i,len;
  const char *type;
  ListOfPoles *buff;
  if(lua_istable(L, narg)==0)
  {
    luaL_error(L, "Error Green's function should be represented by a table");
  }
  lua_getfield(L, narg, "type");
  if(lua_isstring(L,-1)==0)
  {
    luaL_error(L, "Error Reading Green's function without type declaration.");
  }
  type = lua_tolstring(L, -1, NULL);
  if(strcmp(type, "Tri")==0)
  {
    TriDiagonalMatrixType *temp;
    lua_pop(L, 1);
    temp = luaL_checkTriDiagonalMatrixType(L, narg);
    buff = malloc(sizeof(ListOfPoles));
    TriDiagonalMatrixToListOfPoles(*temp, buff);
    FreeTriDiagonalMatrix(temp);
    free(temp);
    return buff;
  }
  if(strcmp(type, "And")==0)
  {
    AndersonMatrixType *temp;
    lua_pop(L, 1);
    temp = luaL_checkAndersonMatrixType(L, narg);
    buff = malloc(sizeof(ListOfPoles));
    AndersonMatrixToListOfPoles(*temp, buff);
    FreeAndersonMatrix(temp);
    free(temp);
    return buff;
  }
  if((strcmp(type, "NaturalImpurityOrbital")==0)||(strcmp(type, "Nat")==0))
  {
    BADoubleTriDiagonalMatrixType *temp;
    lua_pop(L, 1);
    temp = luaL_checkNatrualImpurityOrbitalMatrixType(L, narg);
    buff = malloc(sizeof(ListOfPoles));
    BADoubleTriDiagonalMatrixToListOfPoles(*temp, buff);
    FreeBADoubleTriDiagonalMatrix(temp);
    free(temp);
    return buff;
  }
  if(strcmp(type, "ListOfPoles")==0)
  {
    buff=malloc(sizeof(ListOfPoles));
    lua_pop(L, 1);
    if(lua_objlen(L,narg)!=2)
    {
      luaL_error(L, "Error a list of pole should be a list of two lists {{a0,a1,a2,...,an},{b1^2,b2^2,...,bn^2}}");
    }
    lua_getfield(L, narg, "mu");
    if(lua_isnumber(L, -1)==1)
    {
      buff->mu=(double)luaL_checknumber(L, -1);
    }
    else
    {
      luaL_error(L, "Reading List of poles without mu defined");
    }
    lua_pop(L,1);
    lua_getfield(L, narg, "name");
    if(lua_isstring(L, -1)==1)
    {
      strncpy(buff->Name, luaL_checkstring(L, -1), MAXStrLength);
      buff->Name[MAXStrLength-1]='\0';
    }
    else if(lua_isnil(L, -1)==1)
    {
      buff->Name[0]='\0';
    }
    else
    {
      luaL_error(L, "Reading List of poles without name defined");
    }
    lua_pop(L,1);
    // read a
    lua_rawgeti(L, narg, 1);
    len = (uint32_t)lua_objlen(L, -1);
    if(len==0) luaL_error(L, "Errorr reading Greens function with lenght of a equal to zero");
    if(len==1) luaL_error(L, "Errorr reading Greens function with lenght of a equal to one");
    lua_rawgeti(L, -1, 1);
    buff->a0=(double)luaL_checknumber(L, -1);
    lua_pop(L, 1);
    buff->N=len-1;
    buff->a = (double*)malloc(buff->N*sizeof(double));
    for(i=0;i<buff->N;i++)
    {
      lua_rawgeti(L, -1, i+2);
      buff->a[i]=(double)luaL_checknumber(L, -1);
      lua_pop(L, 1);
    }
    lua_pop(L, 1);
    // read b
    lua_rawgeti(L, narg, 2);
    buff->b=luaL_checkarray_double(L, -1, &len);
    if(len!=buff->N)
    {
      luaL_error(L, "Error Reading List of poles. Length of a (%d) should be one more than length of b (%d)",buff->N+1,len);
    }
    lua_pop(L, 1);
    return buff;
  }
  luaL_error(L, "Errorr reading Greens function of unknown type (%s)",type);
  return NULL;
}


ListOfPoles * luaL_checkListOfPolesVariableGrid(lua_State *L, int narg, const double * const Grid, const uint32_t NGridPoints)
{
  //TODO: all
  return NULL;
}
/*
 *  On the lua stack at position narg should be a Block response function of any type
 *  Returned is a list of poles
 *  transformation to list of poles depends on the setting SuggestedGrid
 *  if NSuggestedGrid=0 the transformation is made exact, thereby increasing the number of Blocks
 *  For a tri-diagonal representation with N Blocks of NBlockSize there will be N*NBlockSize poles
 *  each of NBlockSize.
 *
 *  The format of BlockListOfPoles should be:
 *    { {A0,a1,a2,a3}, {B1,B2,B3}, mu=0, type="ListOfPoles"}
 *    and defines the Green's function
 *    A0 + sum_{i=1^N} B_i * 1/(w-a_i)
 *
 */
BlockListOfPoles * luaL_checkBlockListOfPoles(lua_State *L, int narg)
{
  Bool MakeReal=True;
  uint32_t i,j,k,len;
  const char *type;
  BlockListOfPoles *buff;
  if(lua_istable(L, narg)==0)
  {
    luaL_error(L, "Error Green's function should be represented by a table");
  }
  lua_getfield(L, narg, "type");
  if(lua_isstring(L,-1)==0)
  {
    luaL_error(L, "Error Reading Green's function without type declaration.");
  }
  type = lua_tolstring(L, -1, NULL);
  if(strcmp(type, "Tri")==0)
  {
    BlockTriDiagonalMatrixType *temp;
    lua_pop(L, 1);
    temp = luaL_checkBlockTriDiagonalMatrixType(L, narg);
    buff = malloc(sizeof(BlockListOfPoles));
    BlockTriDiagonalMatrixToBlockListOfPoles(temp, buff);
    FreeBlockTridiagonalMatrix(temp);
    free(temp);
    return buff;
  }
  if(strcmp(type, "And")==0)
  {
    BlockAndersonMatrixType *temp;
    lua_pop(L, 1);
    temp = luaL_checkBlockAndersonMatrixType(L, narg);
    buff = malloc(sizeof(BlockListOfPoles));
    BlockAndersonMatrixToBlockListOfPoles(*temp, buff);
    FreeBlockAndersonMatrix(temp);
    free(temp);
    return buff;
  }
  if((strcmp(type, "NaturalImpurityOrbital")==0)||(strcmp(type, "Nat")==0))
  {
    BlockNaturalImpurityOrbitalMatrixType *temp;
    lua_pop(L, 1);
    temp = luaL_checkBlockNatrualImpurityOrbitalMatrixType(L, narg);
    buff = malloc(sizeof(BlockListOfPoles));
    BlockNaturalImpurityOrbitalMatrixToBlockListOfPoles(*temp, buff);
    FreeBlockNaturalImpurityOrbitalMatrix(temp);
    free(temp);
    return buff;
  }
  if(strcmp(type, "ListOfPoles")==0)
  {
    lua_pop(L, 1);
    buff=malloc(sizeof(BlockListOfPoles));
    if(lua_objlen(L,narg)!=2)
    {
      luaL_error(L, "Error a list of poles should be a list of two lists {{a1,a2,...,an+1},{b1^2,b2^2,...,bn^2}}, with b matrices and a1 optionally a matrix");
    }
    // read b
    lua_rawgeti(L, narg, 2);
    {
      len = (uint32_t)lua_objlen(L, -1);
      if(len==0) luaL_error(L, "Errorr reading Greens function with lenght of b equal to zero");
      buff->N=len;
      lua_rawgeti(L, -1, 1);
      buff->BlockSize = (uint32_t)lua_objlen(L, -1);
      lua_pop(L, 1);
      buff->QComplex=1;
      InitBlockListOfPoles(buff);
      for(i=0;i<buff->N;i++)
      {
        lua_rawgeti(L, -1, i+1); // get b_i+1 onto the stack
        {
          if((uint32_t)lua_objlen(L, -1)!=buff->BlockSize)
          {
            luaL_error(L, "Errorr reading Greens function as block list of poles, but all blocks must be equal in size. BlockSize (%d), size at i=%d is %d",buff->BlockSize,i+1,lua_objlen(L, -1));
          }
          for(j=0;j<buff->BlockSize;j++)
          {
            lua_rawgeti(L,-1,j+1);
            {
              if((uint32_t)lua_objlen(L, -1)!=buff->BlockSize)
              {
                luaL_error(L, "Errorr reading Greens function as block list of poles, but all blocks must be equal in size. BlockSize (%d), size at i=%d j=%d is %d",buff->BlockSize,i+1,j+1,lua_objlen(L, -1));
              }
              for(k=0;k<buff->BlockSize;k++)
              {
                lua_rawgeti(L,-1,k+1);
                {
                  if(lua_isnumber(L, -1))
                  {
                    buff->V[i*(1+Sqr(buff->BlockSize))+1+k+j*buff->BlockSize].Re=(double)luaL_checknumber(L, -1);
                    buff->V[i*(1+Sqr(buff->BlockSize))+1+k+j*buff->BlockSize].Im=0;
                  }
                  else if (luaL_testudata(L, -1, "Complex_Type") != NULL)
                  {
                    MakeReal=False;
                    buff->V[i*(1+Sqr(buff->BlockSize))+1+k+j*buff->BlockSize] = *(Complex16 *)luaL_checkudata(L, -1, "Complex_Type");
                  }
                  else
                  {
                    luaL_error(L, "Error in reading complex table\n");
                  }
                }
                lua_pop(L, 1);
              }
            }
            lua_pop(L, 1);
          }
        }
        lua_pop(L, 1);
      }
    }
    lua_pop(L, 1);
    // read a
    lua_rawgeti(L, narg, 1);
    {
      len = (uint32_t)lua_objlen(L, -1);
      if(len!=buff->N+1)
      {
        luaL_error(L, "Error Reading List of poles. Length of a (%d) should be one more than length of b (%d)",len,buff->N);
      }
      lua_rawgeti(L, -1, 1);
      {
        if(lua_isnumber(L, -1))
        {
          for(i=0;i<Sqr(buff->BlockSize);i++)
          {
            buff->a0[i].Re=0;
            buff->a0[i].Im=0;
          }
          buff->a0[0].Re=(double)luaL_checknumber(L, -1);
          for(i=buff->BlockSize+1;i<Sqr(buff->BlockSize);i+=buff->BlockSize+1)
          {
            buff->a0[i].Re=buff->a0[0].Re;
          }
        }
        else if (luaL_testudata(L, -1, "Complex_Type") != NULL)
        {
          MakeReal=False;
          for(i=0;i<Sqr(buff->BlockSize);i++)
          {
            buff->a0[i].Re=0;
            buff->a0[i].Im=0;
          }
          buff->a0[0]=*(Complex16 *)luaL_checkudata(L, -1, "Complex_Type");
          for(i=buff->BlockSize+1;i<Sqr(buff->BlockSize);i+=buff->BlockSize+1)
          {
            buff->a0[i]=buff->a0[0];
          }
        }
        else
        {
          if((uint32_t)lua_objlen(L, -1)!=buff->BlockSize)
          {
            luaL_error(L, "Errorr reading Greens function as block list of poles, but all blocks must be equal in size. BlockSize (%d), size of a0 is %d",buff->BlockSize,lua_objlen(L, -1));
          }
          for(j=0;j<buff->BlockSize;j++)
          {
            lua_rawgeti(L,-1,j+1);
            {
              if((uint32_t)lua_objlen(L, -1)!=buff->BlockSize)
              {
                luaL_error(L, "Errorr reading Greens function as block list of poles, but all blocks must be equal in size. BlockSize (%d), size of a0(%d) is %d",buff->BlockSize,j+1,lua_objlen(L, -1));
              }
              for(k=0;k<buff->BlockSize;k++)
              {
                lua_rawgeti(L,-1,k+1);
                {
                  if(lua_isnumber(L, -1))
                  {
                    buff->a0[k+j*buff->BlockSize].Re=(double)luaL_checknumber(L, -1);
                    buff->a0[k+j*buff->BlockSize].Im=0;
                  }
                  else if (luaL_testudata(L, -1, "Complex_Type") != NULL)
                  {
                    MakeReal=False;
                    buff->a0[k+j*buff->BlockSize] = *(Complex16 *)luaL_checkudata(L, -1, "Complex_Type");
                  }
                  else
                  {
                    luaL_error(L, "Error in reading complex table\n");
                  }
                }
                lua_pop(L, 1);
              }
            }
            lua_pop(L, 1);
          }
        }
      }
      lua_pop(L, 1);
      for(i=0;i<buff->N;i++)
      {
        lua_rawgeti(L,-1,i+2);
        {
          buff->V[i*(1+Sqr(buff->BlockSize))].Re=(double)luaL_checknumber(L, -1);
          buff->V[i*(1+Sqr(buff->BlockSize))].Im=0;
        }
        lua_pop(L,1);
      }
    }
    lua_pop(L, 1);
    if(MakeReal==True)
    {
      MakeComplexBlockListOfPolesReal(buff);
    }
    lua_getfield(L, narg, "mu");
    if(lua_isnumber(L, -1)==1)
    {
      buff->mu=(double)luaL_checknumber(L, -1);
    }
    else
    {
      luaL_error(L, "Reading List of poles without mu defined");
    }
    lua_pop(L,1);
    lua_getfield(L, narg, "name");
    if(lua_isstring(L, -1)==1)
    {
      strncpy(buff->Name, luaL_checkstring(L, -1), MAXStrLength);
      buff->Name[MAXStrLength-1]='\0';
    }
    else if(lua_isnil(L, -1)==1)
    {
      buff->Name[0]='\0';
    }
    else
    {
      luaL_error(L, "Reading List of poles without name correctly defined");
    }
    lua_pop(L,1);
    return buff;
  }
  luaL_error(L, "Errorr reading Greens function of unknown type (%s)",type);
  return NULL;
}

BlockListOfPoles * luaL_checkBlockListOfPolesVariableGrid(lua_State *L, int narg, const double * const Grid, const uint32_t NGridPoints)
{
  const char *type;
  BlockListOfPoles *buff;
  if(lua_istable(L, narg)==0)
  {
    luaL_error(L, "Error Green's function should be represented by a table");
  }
  lua_getfield(L, narg, "type");
  if(lua_isstring(L,-1)==0)
  {
    luaL_error(L, "Error Reading Green's function without type declaration.");
  }
  type = lua_tolstring(L, -1, NULL);
  lua_pop(L, 1);
  if(strcmp(type, "Tri")==0)
  {
    BlockTriDiagonalMatrixType *temp;
    temp = luaL_checkBlockTriDiagonalMatrixType(L, narg);
    buff = malloc(sizeof(BlockListOfPoles));
    BlockTriDiagonalMatrixToBlockListOfPolesVariableGrid(temp, buff, Grid, NGridPoints);
    FreeBlockTridiagonalMatrix(temp);
    free(temp);
    return buff;
  }
  if(strcmp(type, "And")==0)
  {
    BlockAndersonMatrixType *temp;
    temp = luaL_checkBlockAndersonMatrixType(L, narg);
    buff = malloc(sizeof(BlockListOfPoles));
    BlockAndersonMatrixToBlockListOfPolesVariableGrid(*temp, buff, Grid, NGridPoints);
    FreeBlockAndersonMatrix(temp);
    free(temp);
    return buff;
  }
  if((strcmp(type, "NaturalImpurityOrbital")==0)||(strcmp(type, "Nat")==0))
  {
    BlockNaturalImpurityOrbitalMatrixType *temp;
    temp = luaL_checkBlockNatrualImpurityOrbitalMatrixType(L, narg);
    buff = malloc(sizeof(BlockListOfPoles));
    BlockNaturalImpurityOrbitalMatrixToBlockListOfPolesVariableGrid(*temp, buff, Grid, NGridPoints);
    FreeBlockNaturalImpurityOrbitalMatrix(temp);
    free(temp);
    return buff;
  }
  if(strcmp(type, "ListOfPoles")==0)
  {
    return luaL_checkBlockListOfPoles(L, narg);
  }
  luaL_error(L, "Errorr reading Greens function of unknown type (%s)",type);
  return NULL;
}



/*  reads a response function from the stack
 *  These functions can be in different formats (t, a, n, l) for tri-diagonal, anderson, natural or list of poles
 *  they can be in Block form or not
 *  If a pointer to an excisting structure is returned (user data on stack) a lowercase letter is returned
 *  If a newly created structure is returned (table on stack) an upper case letter is returned
 */
ResponsefunctionType * luaL_checkResponsefunction(lua_State *L, int narg)
{
  char const * type;
  ResponsefunctionType *buff;
  // Block matrix stored as pointer to structure on lua stack
  if (luaL_testudata(L, narg, "ResponseFunction_Type") != NULL)
  {
    buff = (ResponsefunctionType*)luaL_checkudata(L, narg, "ResponseFunction_Type");
    if (buff!=NULL)
    {
      return buff;
    }
  }
  // Block Matrix stored as lua table read sugested grid, type and allocate structure
  if(lua_istable(L, narg)==0)
  {
    luaL_error(L, "Error Green's function should be represented by a table");
  }
  buff=malloc(sizeof(ResponsefunctionType));
  // test if there is a suggested grid present. If present Stored in field SuggestedGrid
  lua_getfield(L, narg, "SuggestedGrid");
  if(lua_isnil(L, -1)==1)
  {
    // no sugested grid
    buff->NSuggestedGridPoints=0;
    buff->SuggestedGrid=NULL;
  }
  else if(lua_istable(L, -1)==1)
  {
    buff->SuggestedGrid = luaL_checkarray_double(L, -1, &(buff->NSuggestedGridPoints));
  }
  else
  {
    luaL_error(L, "Error reading Green's function: \"SuggestedGrid\" should be an table of doubles.");
  }
  lua_pop(L,1);
  // check if Response function is of Block type or not
  lua_rawgeti(L, narg, 1); // puts table a on the stack
  if(lua_objlen(L,-1)==0) luaL_error(L,"Error reading Green's function the length of the table is 0");
  lua_rawgeti(L,-1,1); // puts element A0 or a0 or the stack
  buff->Block = False;
  if(lua_istable(L,-1)) buff->Block = True;
  lua_pop(L,2);
  // read type
  lua_getfield(L, narg, "type");
  if(lua_isstring(L,-1)==0)
  {
    luaL_error(L, "Error Reading Green's function without valid type declaration.");
  }
  type = lua_tolstring(L, -1, NULL);
  if(strcmp(type, "Tri")==0)
  {
    lua_pop(L, 1);
    buff->type='T';
    if(buff->Block==True)
    {
      buff->data=(void*)luaL_checkBlockTriDiagonalMatrixType(L, narg);
    }
    else
    {
      buff->data=(void*)luaL_checkTriDiagonalMatrixType(L, narg);
    }
    return buff;
  }
  if(strcmp(type, "And")==0)
  {
    lua_pop(L, 1);
    buff->type='A';
    if(buff->Block==True)
    {
      buff->data = (void*)luaL_checkBlockAndersonMatrixType(L, narg);
    }
    else
    {
      buff->data = (void*)luaL_checkAndersonMatrixType(L, narg);
    }
    return buff;
  }
  if((strcmp(type, "NaturalImpurityOrbital")==0)||(strcmp(type, "Nat")==0))
  {
    lua_pop(L, 1);
    buff->type='N';
    if(buff->Block==True)
    {
      buff->data = (void*)luaL_checkBlockNatrualImpurityOrbitalMatrixType(L, narg);
    }
    else
    {
      buff->data = (void*)luaL_checkNatrualImpurityOrbitalMatrixType(L, narg);
    }
    return buff;
  }
  if(strcmp(type, "ListOfPoles")==0)
  {
    lua_pop(L, 1);
    buff->type='L';
    if(buff->Block==True)
    {
      if(buff->NSuggestedGridPoints==0)
      {
        buff->data = (void*)luaL_checkBlockListOfPoles(L, narg);
      }
      else
      {
        buff->data = (void*)luaL_checkBlockListOfPolesVariableGrid(L, narg, buff->SuggestedGrid, buff->NSuggestedGridPoints);
      }
    }
    else
    {
      if(buff->NSuggestedGridPoints==0)
      {
        buff->data = (void*)luaL_checkListOfPoles(L, narg);
      }
      else
      {
        buff->data = (void*)luaL_checkListOfPolesVariableGrid(L, narg, buff->SuggestedGrid, buff->NSuggestedGridPoints);
      }
    }
    return buff;
  }
  luaL_error(L, "Errorr reading Greens function of unknown type (%s)\nAllowed are:\nTri\nAnd\nNaturalImpurityOrbital\nListOfPoles",type);
  return NULL;
}


int lua_isBlockGreensfunction(lua_State *L, int narg)
{
  const char *type;
  if(luaL_checkudata (L, narg, "ResponseFunction_Type")!=NULL) return 1;
  if(lua_istable(L, narg)==0) return 0;
  lua_getfield(L, narg, "type");
  if(lua_isstring(L,-1)==0) return 0;
  type = lua_tolstring(L, -1, NULL);
  if(strcmp(type, "Tri")==0) {lua_pop(L, 1);return 1;}
  if(strcmp(type, "And")==0) {lua_pop(L, 1);return 1;}
  if((strcmp(type, "NaturalImpurityOrbital")==0)||(strcmp(type, "Nat")==0)) {lua_pop(L, 1);return 1;}
  if(strcmp(type, "ListOfPoles")==0) {lua_pop(L, 1);return 1;}
  lua_pop(L, 1);
  return 0;
}

int lua_isResponsefunction(lua_State *L, int narg)
{
  return lua_isBlockGreensfunction(L, narg);
}

WaveFunctionType* luaL_checkarray_WaveFunctionType(lua_State *L, int narg, uint32_t *len_out)
{
  int i;
  luaL_checktype(L, narg, LUA_TTABLE);
  *len_out = (uint32_t)lua_objlen(L, narg);
  WaveFunctionType *buff = (WaveFunctionType*)malloc((*len_out)*sizeof(WaveFunctionType));
  for(i=0;i<(*len_out);i++)
  {
    lua_rawgeti(L, narg, i+1);
    buff[i] = *(WaveFunctionType *)luaL_checkudata(L, -1, "Wavefunction_Type");
    lua_pop(L, 1);
  }
  return buff;
}

WaveFunctionType** luaL_checkarray_PointerWaveFunctionType(lua_State *L, int narg, uint32_t *len_out)
{
  int i;
  luaL_checktype(L, narg, LUA_TTABLE);
  *len_out = (uint32_t)lua_objlen(L, narg);
  WaveFunctionType **buff = (WaveFunctionType**)malloc((*len_out)*sizeof(WaveFunctionType*));
  for(i=0;i<(*len_out);i++)
  {
    lua_rawgeti(L, narg, i+1);
    buff[i] = (WaveFunctionType *)luaL_checkudata(L, -1, "Wavefunction_Type");
    lua_pop(L, 1);
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

void LuaPushArray_Complex16(lua_State *L, Complex16 const * const val, int const N)
{
  uint32_t i;
  Complex16 *res;
  lua_createtable(L, N, 0);
  for(i=0;i<N;i++)
  {
    if(val[i].Im==0)
    {
      lua_pushnumber(L, val[i].Re);
    }
    else
    {
      res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
      res->Re=val[i].Re;
      res->Im=val[i].Im;
      luaL_setmetatable(L, "Complex_Type");
    }
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


void LuaPushTable_WaveFunctionType_RemoveEmpty(lua_State* L, WaveFunctionType * const psi, uint32_t const Npsi )
{
  uint32_t i,j=0;
  WaveFunctionType *buff;
  lua_newtable(L);
  for(i=0;i<Npsi;i++)
  {
    if(psi[i].N==0)
    {
      WaveFunctionFree(&(psi[i]));
    }
    else
    {
      j++;
      buff = (WaveFunctionType*) lua_newuserdata(L, sizeof(WaveFunctionType));
      luaL_setmetatable(L, "Wavefunction_Type");
      *buff = psi[i];
      lua_rawseti(L, -2, j);
    }
  }
}

void LuaPushTable_IdentityMatrix(lua_State *L, uint32_t dim)
{
  uint32_t i,j;
  lua_createtable(L,dim,0);
  for(i=1;i<=dim;i++)
  {
    lua_createtable(L,dim,0);
    for(j=1;j<=dim;j++)
    {
      lua_pushnumber(L, ((i==j)?1:0));
      lua_rawseti(L,-2,j);
    }
    lua_rawseti(L,-2,i);
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
  /* set metatable */
  lua_getglobal(L, "MatrixMeta");
  lua_setmetatable(L, -2);
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
        if(M->V[i*M->Nj+j].Im==0)
        {
          lua_pushnumber(L, M->V[i*M->Nj+j].Re);
        }
        else
        {
          res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
          res->Re = M->V[i*M->Nj+j].Re;
          res->Im = M->V[i*M->Nj+j].Im;
          luaL_setmetatable(L, "Complex_Type");
        }
      }
      else
      {
        lua_pushnumber(L, M->VRe[i*M->Nj+j]);
      }
      lua_rawseti(L, -2, j+1);
    }
    lua_rawseti(L, -2, i+1);
  }
  /* set metatable */
  lua_getglobal(L, "MatrixMeta");
  lua_setmetatable(L, -2);
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
  lua_createtable(L, 2, 3);
  /* push a */
  LuaPushArray_double(L, M->a, MNtrunc+1);
  lua_rawseti(L, -2, 1);
  /* push b */
  LuaPushArray_double(L, M->b, MNtrunc);
  lua_rawseti(L, -2, 2);
  /* push mu */
  lua_pushnumber(L, M->mu);
  lua_setfield(L,-2,"mu");
  /* push type */
  lua_pushstring(L, "Tri");
  lua_setfield(L,-2,"type");
  /* push name */
  lua_pushstring(L, M->Name);
  lua_setfield(L,-2,"name");
  /* set metatable */
  lua_getglobal(L, "ResponseFunctionMeta");
  lua_setmetatable(L, -2);
}

void LuaPushTable_BlockTriDiagonalMatrixType(lua_State* L, BlockTriDiagonalMatrixType const * const M)
{
  uint32_t k,i,MNtrunc=M->N;
  double totb;
  CompactMatrixType block;
  //*   matrix k can be found at position [k*BlockSize*BlockSize]
  //*   element [i,j] of matrix k can be found at position [k*BlockSize*BlockSize+i*BlockSize+j]
  //*   stop the printing if the sum of the absolute value of all elements of b[i] is less then 0.01*DBL_EPSILON
  if(M->QComplex)
  {
    for(k=1;k<M->N;k++)
    {
      totb=0;
      for(i=k*M->BlockSize*M->BlockSize;i<(k+1)*M->BlockSize*M->BlockSize;i++)
      {
        totb+=M->b[i].Re*M->b[i].Re+M->b[i].Im*M->b[i].Im;
      }
      if(totb<0.01*DBL_EPSILON)
      {
        MNtrunc=k;
        break;
      }
    }
  }
  else
  {
    for(k=1;k<M->N;k++)
    {
      totb=0;
      for(i=k*M->BlockSize*M->BlockSize;i<(k+1)*M->BlockSize*M->BlockSize;i++)
      {
        totb+=M->bRe[i]*M->bRe[i];
      }
      if(totb<0.01*DBL_EPSILON)
      {
        MNtrunc=k;
        break;
      }
    }
  }
  block.QComplex=M->QComplex;
  block.Ni=M->BlockSize;
  block.Nj=M->BlockSize;
  block.a0=0;
  block.b0=0;
  block.mu=0;
  block.Name=NULL;
  lua_createtable(L, 2, 2);
  /* push a */
  if(block.QComplex)
  {
    block.V = M->a;
  }
  else
  {
    block.VRe = M->aRe;
  }
  lua_createtable(L, MNtrunc+1, 0);
  for(k=0;k<=MNtrunc;k++)
  {
    LuaPushTable_CompactMatrixType(L, &block);
    lua_rawseti(L, -2, k+1);
    if(block.QComplex)
    {
      block.V += Sqr(M->BlockSize);
    }
    else
    {
      block.VRe += Sqr(M->BlockSize);
    }
  }
  lua_rawseti(L, -2, 1);
  /* push b */
  if(block.QComplex)
  {
    block.V = M->b;
  }
  else
  {
    block.VRe = M->bRe;
  }
  lua_createtable(L, MNtrunc, 0);
  for(k=0;k<MNtrunc;k++)
  {
    LuaPushTable_CompactMatrixType(L, &block);
    lua_rawseti(L, -2, k+1);
    if(block.QComplex)
    {
      block.V += Sqr(M->BlockSize);
    }
    else
    {
      block.VRe += Sqr(M->BlockSize);
    }
  }
  lua_rawseti(L, -2, 2);
  /* push mu */
  lua_pushnumber(L, M->mu);
  lua_setfield(L,-2,"mu");
  /* push type */
  lua_pushstring(L, "Tri");
  lua_setfield(L,-2,"type");
  /* push name */
  lua_pushstring(L, M->Name);
  lua_setfield(L,-2,"name");
  /* set metatable */
  lua_getglobal(L, "ResponseFunctionMeta");
  lua_setmetatable(L, -2);
}

void LuaPushTable_AndersonMatrixType(lua_State* L, AndersonMatrixType const * const M)
{
  lua_createtable(L, 2, 3);
  /* push a */
  LuaPushArray_double(L, M->a, M->N+1);
  lua_rawseti(L, -2, 1);
  /* push b */
  LuaPushArray_double(L, M->b, M->N);
  lua_rawseti(L, -2, 2);
  /* push mu */
  lua_pushnumber(L, M->mu);
  lua_setfield(L,-2,"mu");
  /* push type */
  lua_pushstring(L, "And");
  lua_setfield(L,-2,"type");
  /* push name */
  lua_pushstring(L, M->Name);
  lua_setfield(L,-2,"name");
  /* set metatable */
  lua_getglobal(L, "ResponseFunctionMeta");
  lua_setmetatable(L, -2);
}

void LuaPushTable_BlockAndersonMatrixType(lua_State* L, BlockAndersonMatrixType const * const M)
{
  uint32_t k;
  CompactMatrixType block;
  //*   matrix k can be found at position [k*BlockSize*BlockSize]
  //*   element [i,j] of matrix k can be found at position [k*BlockSize*BlockSize+i*BlockSize+j]
  //*   stop the printing if the sum of the absolute value of all elements of b[i] is less than 0.01*DBL_EPSILON
  block.QComplex=M->QComplex;
  block.Ni=M->BlockSize;
  block.Nj=M->BlockSize;
  block.a0=0;
  block.b0=0;
  block.mu=0;
  block.Name=NULL;
  lua_createtable(L, 2, 3);
  /* push a */
  if(block.QComplex)
  {
    block.V = M->a;
  }
  else
  {
    block.VRe = M->aRe;
  }
  lua_createtable(L, M->N+1, 0);
  for(k=0;k<=M->N;k++)
  {
    LuaPushTable_CompactMatrixType(L, &block);
    lua_rawseti(L, -2, k+1);
    if(block.QComplex)
    {
      block.V += Sqr(M->BlockSize);
    }
    else
    {
      block.VRe += Sqr(M->BlockSize);
    }
  }
  lua_rawseti(L, -2, 1);
  /* push b */
  if(block.QComplex)
  {
    block.V = M->b;
  }
  else
  {
    block.VRe = M->bRe;
  }
  lua_createtable(L, M->N, 0);
  for(k=0;k<M->N;k++)
  {
    LuaPushTable_CompactMatrixType(L, &block);
    lua_rawseti(L, -2, k+1);
    if(block.QComplex)
    {
      block.V += Sqr(M->BlockSize);
    }
    else
    {
      block.VRe += Sqr(M->BlockSize);
    }
  }
  lua_rawseti(L, -2, 2);
  /* push mu */
  lua_pushnumber(L, M->mu);
  lua_setfield(L,-2,"mu");
  /* push type */
  lua_pushstring(L, "And");
  lua_setfield(L,-2,"type");
  /* push name */
  lua_pushstring(L, M->Name);
  lua_setfield(L,-2,"name");
  /* set metatable */
  lua_getglobal(L, "ResponseFunctionMeta");
  lua_setmetatable(L, -2);
}

void LuaPushTable_ListOfPoles(lua_State* L, ListOfPoles const * const M)
{
  uint32_t i;
  lua_createtable(L, 2, 3);
  /* push a */
  lua_createtable(L, M->N+1, 0);
  lua_pushnumber(L, M->a0);
  lua_rawseti(L, -2, 1);
  for(i=0;i<M->N;i++)
  {
    lua_pushnumber(L, M->a[i]);
    lua_rawseti(L, -2, i+2);
  }
  lua_rawseti(L, -2, 1);
  /* push b */
  LuaPushArray_double(L, M->b, M->N);
  lua_rawseti(L, -2, 2);
  /* push mu */
  lua_pushnumber(L, M->mu);
  lua_setfield(L,-2,"mu");
  /* push type */
  lua_pushstring(L, "ListOfPoles");
  lua_setfield(L,-2,"type");
  /* push name */
  lua_pushstring(L, M->Name);
  lua_setfield(L,-2,"name");
  /* set metatable */
  lua_getglobal(L, "ResponseFunctionMeta");
  lua_setmetatable(L, -2);
}

void LuaPushTable_BlockListOfPoles(lua_State* L, BlockListOfPoles const * const M)
{
  uint32_t k;
  CompactMatrixType block;
  block.QComplex=M->QComplex;
  block.Ni=M->BlockSize;
  block.Nj=M->BlockSize;
  block.a0=0;
  block.b0=0;
  block.mu=0;
  block.Name=NULL;
  lua_createtable(L, 2, 3);
  /* push a */
  if(block.QComplex)
  {
    block.V = M->a0;
  }
  else
  {
    block.VRe = M->a0Re;
  }
  lua_createtable(L, M->N+1, 0);
  LuaPushTable_CompactMatrixType(L, &block);
  lua_rawseti(L, -2, 1);
  if(M->QComplex)
  {
    for(k=0;k<M->N;k++)
    {
      lua_pushnumber(L, M->V[k*(1+Sqr(M->BlockSize))].Re);
      lua_rawseti(L, -2, k+2);
    }
  }
  else
  {
    for(k=0;k<M->N;k++)
    {
      lua_pushnumber(L, M->VRe[k*(1+Sqr(M->BlockSize))]);
      lua_rawseti(L, -2, k+2);
    }
  }
  lua_rawseti(L, -2, 1);
  /* push b */
  if(block.QComplex)
  {
    block.V = M->V+1;
  }
  else
  {
    block.VRe = M->VRe+1;
  }
  lua_createtable(L, M->N, 0);
  for(k=0;k<M->N;k++)
  {
    LuaPushTable_CompactMatrixType(L, &block);
    lua_rawseti(L, -2, k+1);
    if(block.QComplex)
    {
      block.V += Sqr(M->BlockSize)+1;
    }
    else
    {
      block.VRe += Sqr(M->BlockSize)+1;
    }
  }
  lua_rawseti(L, -2, 2);
  /* push mu */
  lua_pushnumber(L, M->mu);
  lua_setfield(L,-2,"mu");
  /* push type */
  lua_pushstring(L, "ListOfPoles");
  lua_setfield(L,-2,"type");
  /* push name */
  lua_pushstring(L, M->Name);
  lua_setfield(L,-2,"name");
  /* set metatable */
  lua_getglobal(L, "ResponseFunctionMeta");
  lua_setmetatable(L, -2);
}

void LuaPushTable_NaturalImpurityOrbital(lua_State* L, BADoubleTriDiagonalMatrixType * const M)
{
  lua_createtable(L, 1, 3);
  /* push {a0, b0} */
  lua_createtable(L,2,0);
  lua_pushnumber(L,M->a0);
  lua_rawseti(L, -2, 1);
  lua_pushnumber(L,M->b0);
  lua_rawseti(L, -2, 2);
  lua_rawseti(L, -2, 1);
  /* push val */
  LuaPushTable_TriDiagonalMatrixType(L, &(M->val));
  lua_setfield(L, -2, "val");
  /* push con */
  LuaPushTable_TriDiagonalMatrixType(L, &(M->con));
  lua_setfield(L, -2, "con");
  /* push mu */
  lua_pushnumber(L, M->mu);
  lua_setfield(L,-2,"mu");
  /* push Epsilon */
  lua_pushnumber(L, M->Epsilon);
  lua_setfield(L,-2,"epsilon");
  /* push type */
  lua_pushstring(L, "NaturalImpurityOrbital");
  lua_setfield(L,-2,"type");
  /* push name */
  lua_pushstring(L, M->Name);
  lua_setfield(L,-2,"name");
  /* set metatable */
  lua_getglobal(L, "ResponseFunctionMeta");
  lua_setmetatable(L, -2);
}


void LuaPushTable_BlockNaturalImpurityOrbital(lua_State* L, BlockNaturalImpurityOrbitalMatrixType * const M)
{
  lua_createtable(L, 1, 3);
  /* push {a0, b0} */
  lua_createtable(L,2,0);
  LuaPushTable_CompactMatrixType(L, &(M->a0));
  lua_rawseti(L, -2, 1);
  LuaPushTable_IdentityMatrix(L,M->BlockSize);
  lua_rawseti(L, -2, 2);
  lua_rawseti(L, -2, 1);
  /* push val */
  LuaPushTable_BlockTriDiagonalMatrixType(L, &(M->val));
  lua_setfield(L, -2, "val");
  /* push con */
  LuaPushTable_BlockTriDiagonalMatrixType(L, &(M->con));
  lua_setfield(L, -2, "con");
  /* push mu */
  lua_pushnumber(L, M->mu);
  lua_setfield(L,-2,"mu");
  /* push Epsilon */
  lua_pushnumber(L, M->Epsilon);
  lua_setfield(L,-2,"epsilon");
  /* push type */
  lua_pushstring(L, "NaturalImpurityOrbital");
  lua_setfield(L,-2,"type");
  /* push name */
  lua_pushstring(L, M->Name);
  lua_setfield(L,-2,"name");
  /* set metatable */
  lua_getglobal(L, "ResponseFunctionMeta");
  lua_setmetatable(L, -2);
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
    if(M->Spec[i].ResonantEnergy!=NAN)
    {
      lua_pushnumber(L, M->Spec[i].ResonantEnergy);
      lua_setfield(L,-2,"ResonantEnergy");
    }
    lua_rawseti(L, -2, i+1);
  }
  return;
}


/*  push a Block Greens function to the stack
 *  These functions can be in different formats (t, a, n, l) for tri-diagonal, anderson, natural or list of poles
 *  If a pointer to an excisting structure is returned (user data on stack) a lowercase letter is returned
 *  If a newly created structure is returned (table on stack) an upper case letter is returned
 */
void LuaPush_ResponseFunction(lua_State *L, ResponsefunctionType const * const M)
{
  ResponsefunctionType *res;
  if((M->type|32)==M->type) //lowercase
  {
    /* push a user data to the stack */
    res = (ResponsefunctionType *) lua_newuserdata(L, sizeof(ResponsefunctionType));
    *res = *M;
    luaL_setmetatable(L, "ResponseFunction_Type");
    return;
  }
  switch (M->type)
  {
    case 'T':
    {
      if(M->Block==True)
      {
        LuaPushTable_BlockTriDiagonalMatrixType(L, M->data);
      }
      else
      {
        LuaPushTable_TriDiagonalMatrixType(L, M->data);
      }
      return;
    }
    case 'A':
    {
      if(M->Block==True)
      {
        LuaPushTable_BlockAndersonMatrixType(L, M->data);
      }
      else
      {
        LuaPushTable_AndersonMatrixType(L, M->data);
      }
      return;
    }
    case 'L':
    {
      if(M->Block==True)
      {
        LuaPushTable_BlockListOfPoles(L, M->data);
      }
      else
      {
        LuaPushTable_ListOfPoles(L, M->data);
      }
      return;
    }
    case 'N':
    {
      if(M->Block==True)
      {
        LuaPushTable_BlockNaturalImpurityOrbital(L, M->data);
      }
      else
      {
        LuaPushTable_NaturalImpurityOrbital(L, M->data);
      }
      return;
    }
    default:
    {
      printf("Error in LuaPush_ResponseFunction: Type of response function should be either t/T,a/A,l/L, or n/N\n");
      return;
    }
  }
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
  uint32_t *buffuint, len;
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
  if (strcmp(type, "array-uint32_t-3") ==0)
  {
    buffuint = luaL_checkarray_unit32_t(L, narg, &len);
    if(len!=3)
      buffuint=NULL;
    return buffuint;
  }
  if (strcmp(type, "Vector") ==0)
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
  buff=NULL;
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

int strcicmp(char const *a, char const *b)
{
  int d;
  if(!(*a&&*b))
    return 1;
  for (;; a++, b++)
  {
    d = tolower(*a) - tolower(*b);
    if (d != 0 || !*a)
      return d;
  }
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
      if (strcicmp(tokenread, Options[j].name) == 0)
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
        else if(strcmp(Options[j].type, "array-uint32_t-3") == 0)
        {
          *(uint32_t**)Options[j].variable=(uint32_t*)LuaReturnVariableType(L, -1, Options[j].type);
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


static int LuaReadSingleSelfenergy(lua_State *L, int narg, OperatorType *V, SpectraTypeArray *SonGrid, ResponsefunctionType *SasMatrix)
{
  uint32_t n;
  /* read self energy */
  /* the self energy can be */
  // 1 a single object of spectra type
  // 2 a single object of ResponseFunction type
  // 3 A list of length 2 containing an operator or a list of operators and a self energy in format 1 or 2
  if ((lua_type(L,narg)==LUA_TUSERDATA)&&(luaL_testudata(L, narg, "Spectra_Type") != NULL))
  { // single self energy in spectra format
    V=NULL;
    SasMatrix=NULL;
    SonGrid = (SpectraTypeArray *) luaL_checkudata(L, narg, "Spectra_Type");
    if(SonGrid->NSpectra!=1)
    {
      return luaL_error(L, "Error in self energy. A spectrum object with more than one spectrum is given\n");
    }
  }
  else if(lua_isResponsefunction(L, narg))
  { // single self energy as matrix type (either tri, and, lp, or nai
    V=NULL;
    SonGrid=NULL;
    SasMatrix=luaL_checkResponsefunction(L, narg);
  }
  else if((lua_type(L,narg)==LUA_TTABLE)&&(lua_objlen(L, narg)==2))
  { // single or list of operators followed by a spectral function or response function as matrix type
    lua_rawgeti(L,narg,1);
    if (lua_type(L,narg)==LUA_TTABLE)
    {
      V = luaL_checkarray_OperatorType(L, -1, &n);
    }
    else
    {
      n=1;
      V = (OperatorType*)malloc(sizeof(OperatorType));
      V[0] = *(OperatorType *)luaL_checkudata(L, -1, "Operator_Type");
    }
    lua_pop(L,1);
    lua_rawgeti(L,narg,2);

    if ((lua_type(L,-1)==LUA_TUSERDATA)&&(luaL_testudata(L, -1, "Spectra_Type") != NULL))
    { // self energy in spectra format
      SasMatrix=NULL;
      SonGrid = (SpectraTypeArray *) luaL_checkudata(L, narg, "Spectra_Type");
      if(SonGrid->NSpectra!=n*n)
      {
        return luaL_error(L, "Error in self energy. A spectrum object for %d operators should have %d*%d elements instead of %d\n",n,n,n,SonGrid->NSpectra);
      }
    }
    else if(lua_isResponsefunction(L, narg))
    { // single self energy as matrix type (either tri, and, lp, or nai
      SonGrid=NULL;
      SasMatrix=luaL_checkResponsefunction(L, narg);
      if(((SasMatrix->type|32)=='t')&&(((BlockTriDiagonalMatrixType*)(SasMatrix->data))->BlockSize!=n)) return luaL_error(L, "Error dimension of self energy representation in tridiagonal form is not the same as number of operators\n");
      if(((SasMatrix->type|32)=='a')&&(((BlockAndersonMatrixType*)(SasMatrix->data))->BlockSize!=n)) return luaL_error(L, "Error dimension of self energy representation in anderson impurity matrix form is not the same as number of operators\n");
      if(((SasMatrix->type|32)=='l')&&(((BlockListOfPoles*)(SasMatrix->data))->BlockSize!=n)) return luaL_error(L, "Error dimension of self energy representation as list of poles form is not the same as number of operators\n");
      if(((SasMatrix->type|32)=='n')&&(((BlockNaturalImpurityOrbitalMatrixType*)(SasMatrix->data))->BlockSize!=n)) return luaL_error(L, "Error dimension of self energy representation in natural impurity orbital form is not the same as number of operators\n");
    }
    else
    {
      return luaL_error(L,"Error in format of self energy with operators, should be of Spectra_Type or ResponsFunction_Type\n");
    }
    lua_pop(L,1);
    return 0;
  }
  return luaL_error(L, "Error reading self energy, wrong format. Self energy should be either of Spectra_Type, of ResponsFunction_Type or a table containing an operator or table of operators and one of the previous two\n");
}

uint32_t FreeSelfEnergyArrayFromLua(uint32_t NSigma, OperatorType **V, SpectraTypeArray **SonGrid, ResponsefunctionType **SasMatrix)
{
  uint32_t i;
  for(i=0;i<NSigma;i++)
  {
    // Operators are copied into an array by dereferencing the OperatorType structure. The pointers do not need to be freed, the operator needs to be.
    if(V[i]!=NULL)
    {
      free(V[i]);
    }
    // Spectra type is all copied by pointers, Lua controls them
    // Response functions need care capital and small letters for type are allocated differently
    if(SasMatrix[i]!=NULL)
    {
      if((SasMatrix[i]->type&32)==0)
      {
        // capital letter for type, read as lua table and thus allocated the structure
        FreeBlockMatrix(SasMatrix[i]);
        free(SasMatrix[i]);
      }
    }
  }
  free(V);
  free(SonGrid);
  free(SasMatrix);
  return 0;
}
