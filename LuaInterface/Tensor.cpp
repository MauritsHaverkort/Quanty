/*
 *        Functions in library Tensor.*****
 */

/*  On stack 
 *  rank positive integer
 *  dimension either positive integer or list of length rank of positive integers
 *  optional variant list of rank bools
 *  Possible last input options {{"Userdata",true}}
 *
 *  output tensor object on stack
 */
static int LuaTensorZero(lua_State *L)
{
  int i;
  uint32_t rank;
  CompactTensorType *res;
  bool Userdata=false;
  uint32_t len = (uint32_t)lua_gettop(L);
  LuaOptionsType LuaFunctionOptions[]={
    {"Userdata","bool",&(Userdata)},
    {NULL, NULL, NULL}};
  if(LuaTestOptionsList(L,-1)==1)
  {
    LuaCheckOptions(L, -1, LuaFunctionOptions);
    len=len-1;
  }
  if(!((len==2)||(len==3)))
  {
    return luaL_error(L,"Tensor Zero called with %d arguments, 2 or 3 arguments expected\n", len);
  }
  if(Userdata==true)
  {
    res = (CompactTensorType *) lua_newuserdata(L, sizeof(CompactTensorType));
    luaL_setmetatable(L, "Tensor_Type");
    res->OnLuaStack=true;
  }
  else
  {
    res = (CompactTensorType*)malloc(sizeof(CompactTensorType));
  }
  res->QComplex=0;
  res->rank = (uint32_t)luaL_checkinteger(L, 1);
  if(res->rank==0)  luaL_error(L,"Trying to create a rank 0 tensor\n");
  if(lua_isnumber(L,2)==1)
  {
    res->dimensions = (uint32_t*)malloc(res->rank*sizeof(uint32_t));
    res->dimensions[0] = (uint32_t)luaL_checkunsigned(L, 2);
    for(i=1;i<res->rank;i++) res->dimensions[i]=res->dimensions[0];
  }
  else if(lua_istable(L,2)==1)
  {
    res->dimensions = luaL_checkarray_unit32_t(L, 2, &rank);
    if(rank!=res->rank) luaL_error(L,"Lengt of table defining the dimensions at position 2 should be equal to the rank of the tensor\n");
  }
  else return luaL_error(L,"second element in call to Tensor Zero should either be an integer spcifying the dimension of all ranks or a list of length rank");
  if(len==3)
  {
    res->variant=luaL_checkarray_bool(L, 2, &rank);
    if(rank!=res->rank) luaL_error(L,"Lengt of table defining the variant at position 3 should be equal to the rank of the tensor\n");
  }
  else
  {
    res->variant = (bool*)malloc(res->rank*sizeof(bool));
    for(i=0;i<res->rank;i++) res->variant[i]=false;
  }
  InitCompactTensor(res);
  if(Userdata==false)
  {
    res->OnLuaStack=false;
    LuaPushTable_CompactTensorType(L, res);
    FreeCompactTensor(res);
    free(res);
  }
  return 1;
}
static int LuaTensorIdentity(lua_State *L)
{
  int i,j,di;
  uint32_t rank;
  CompactTensorType *res;
  bool Userdata=false;
  uint32_t len = (uint32_t)lua_gettop(L);
  LuaOptionsType LuaFunctionOptions[]={
    {"Userdata","bool",&(Userdata)},
    {NULL, NULL, NULL}};
  if(LuaTestOptionsList(L,-1)==1)
  {
    LuaCheckOptions(L, -1, LuaFunctionOptions);
    len=len-1;
  }
  if(!((len==2)||(len==3)))
  {
    return luaL_error(L,"Tensor Zero called with %d arguments, 2 or 3 arguments expected\n", len);
  }
  if(Userdata==true)
  {
    res = (CompactTensorType *) lua_newuserdata(L, sizeof(CompactTensorType));
    luaL_setmetatable(L, "Tensor_Type");
    res->OnLuaStack=true;
  }
  else
  {
    res = (CompactTensorType*)malloc(sizeof(CompactTensorType));
  }
  res->QComplex=0;
  res->rank = (uint32_t)luaL_checkinteger(L, 1);
  if(res->rank==0)  luaL_error(L,"Trying to create a rank 0 tensor\n");
  if(lua_isnumber(L,2)==1)
  {
    res->dimensions = (uint32_t*)malloc(res->rank*sizeof(uint32_t));
    res->dimensions[0] = (uint32_t)luaL_checkunsigned(L, 2);
    for(i=1;i<res->rank;i++) res->dimensions[i]=res->dimensions[0];
  }
  else if(lua_istable(L,2)==1)
  {
    res->dimensions = luaL_checkarray_unit32_t(L, 2, &rank);
    if(rank!=res->rank) luaL_error(L,"Lengt of table defining the dimensions at position 2 should be equal to the rank of the tensor\n");
  }
  else return luaL_error(L,"second element in call to Tensor Zero should either be an integer spcifying the dimension of all ranks or a list of length rank");
  if(len==3)
  {
    res->variant=luaL_checkarray_bool(L, 2, &rank);
    if(rank!=res->rank) luaL_error(L,"Lengt of table defining the variant at position 3 should be equal to the rank of the tensor\n");
  }
  else
  {
    res->variant = (bool*)malloc(res->rank*sizeof(bool));
    for(i=0;i<res->rank;i++) res->variant[i]=false;
  }
  InitCompactTensor(res);
  di=1;
  for(j=1;j<res->rank;j++)
  {
    di *= res->dimensions[j];
    di += 1;
  }
  for(i=0;i<res->dimensions[0];i++)
  {
    for(j=1;j<res->rank;j++)
    {
      if(i==res->dimensions[j]) goto finishedFillingOnes;
    }
    if(res->QComplex==1)
    {
      res->V[i*di].Re=1;
    }
    else
    {
      res->VRe[i*di]=1;
    }
  }
  finishedFillingOnes:
  if(Userdata==false)
  {
    res->OnLuaStack=false;
    LuaPushTable_CompactTensorType(L, res);
    FreeCompactTensor(res);
    free(res);
  }
  return 1;
}

static int LuaTensorToTable(lua_State *L)
{
  CompactTensorType *A;
  A = (CompactTensorType *) luaL_checkCompactTensor(L, 1);
  LuaPushTable_CompactTensorType(L, A);
  FreeCompactTensor(A);
  return 1;}

static int LuaTensorToUserdata(lua_State *L)
{
  CompactTensorType *A,*B;
  A = (CompactTensorType *) luaL_checkCompactTensor(L, 1);
  B = (CompactTensorType *) lua_newuserdata(L, sizeof(CompactTensorType));
  luaL_setmetatable(L, "Tensor_Type");
  if (luaL_testudata(L, 1, "Tensor_Type") == NULL)
  {
    *B=*A;
    free(A);
  }
  else
  {
    CompactTensorCopy(A, B);
  }
  B->OnLuaStack=true;
  return 1;

}



/*
 *       Meta functions on tensors
 */

static int LuaTensorUnm(lua_State *L)
{
  //TODO:
  return luaL_error(L,"TODO");
  return 0;
}

static int LuaTensorAdd(lua_State *L)
{
  //TODO:
  return luaL_error(L,"TODO");
  return 0;
}

static int LuaTensorSub(lua_State *L)
{
  //TODO:
  return luaL_error(L,"TODO");
  return 0;
}

static int LuaTensorMul(lua_State *L)
{
  //TODO:
  return luaL_error(L,"TODO");
  return 0;
}

static int LuaTensorDiv(lua_State *L)
{
  //TODO:
  return luaL_error(L,"TODO");
  return 0;
}

static int LuaTensorPow(lua_State *L)
{
  //TODO:
  return luaL_error(L,"TODO");
  return 0;
}

static int LuaTensorEq(lua_State *L)
{
  //TODO:
  return luaL_error(L,"TODO");
  return 0;
}

static int LuaTensorLen(lua_State *L)
{
  //TODO:
  return luaL_error(L,"TODO");
  return 0;
}

static int LuaTensorToString(lua_State *L)
{
  CompactTensorType *M;
  DynamicString Mchar;
  M = (CompactTensorType*)luaL_checkudata(L, 1, "Tensor_Type");
  CompactTensorToString(M, &Mchar);
  lua_pushstring(L,Mchar.str);
  FreeDynamicString(&Mchar);
  FreeCompactTensor(M);
  return 1;
}

static int LuaTensorIndex(lua_State *L)
{
  //TODO:
  return luaL_error(L,"TODO");
  return 0;
}

static int LuaTensorDestroy(lua_State *L)
{
  CompactTensorType *M;
  M = (CompactTensorType*)luaL_checkudata(L, 1, "Tensor_Type");
  M->OnLuaStack=false; //give back the data
  FreeCompactTensor(M);
  return 0;
}





static const struct luaL_Reg Tensor_func[] = {
  {"Zero", LuaTensorZero},
  {"Identity", LuaTensorIdentity},
  {"ToTable", LuaTensorToTable},
  {"ToUserdata", LuaTensorToUserdata},
  {NULL, NULL}
};

static const struct luaL_Reg Tensor_meta[] = {
  {"__unm",      LuaTensorUnm},
  {"__add",      LuaTensorAdd},
  {"__sub",      LuaTensorSub},
  {"__mul",      LuaTensorMul},
  {"__div",      LuaTensorDiv},
  {"__pow",      LuaTensorPow},
  {"__eq",       LuaTensorEq},
  {"__len",      LuaTensorLen},
  {"__tostring", LuaTensorToString},
  {"__index",    LuaTensorIndex},
  {"__gc",       LuaTensorDestroy},
  {NULL, NULL}
};

static const struct luaL_Reg TensorTable_meta[] = {
  {"__unm", LuaTensorUnm},
  {"__add", LuaTensorAdd},
  {"__sub", LuaTensorSub},
  {"__mul", LuaTensorMul},
  {"__div", LuaTensorDiv},
  {"__pow", LuaTensorPow},
  {"__eq",  LuaTensorEq},
  {NULL, NULL}
};

