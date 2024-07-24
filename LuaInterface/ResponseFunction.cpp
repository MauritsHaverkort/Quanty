/*  Response functions are stored either as matrix representations and reside fully on the lua stack as tables
 *                                or as ResponsefunctionType
 *
 *  There are different types of response functions
 *    Tri
 *    And
 *    ListOfPoles
 *    NaturalImpurityOrbital
 *  each type can either be single or matrix represented
 *
 *  In order for the meta methods to function for table representations one needs to set the meta table "ResponseFunctionMeta" in lua
 *  There is no type-check, (i.e. setting the meta table is optional) for the functions acting on these tables.
 *
 *
 *
 */


/*   change between l,t,a, and n type
 *   change between representation as table or usertype (table capital, user data small)
 *   parameters on lua stack (responsefunction, typedef)
 *   typedef being a string of the form:
 *   Tri, And, ListOfPoles, NaturalImpurityOrbital, Userdata, Table, Dense
 *
 *   SingularValue is an option used to determine if the block size needs to be reduced in a tridiagonalization routine
 *   EffectivekBT is used to determine if a state belongs to the valence or conduction bands (occupied or empty using Fermi Dirac statistics)
 *   Deflate is a boolean option to determine if the final dimension should be the same (adding zero coupling) or is allowed to be reduced
 */
static int LuaChangeTypeResponseFunction(lua_State *L)
{
  ResponsefunctionType *Min, Mout;
  char const * convert;
  char type;
  double T=sqrt(DBL_MIN);
  double SingularValue=1000*DBL_EPSILON;
  bool deflate=true;
  //TODO: implement deflation when possible and deflate=true
  uint32_t len = (uint32_t)lua_gettop(L);
  LuaOptionsType LuaChangeTypeResponseFunctionOptions[]={
    {"SingularValue","double",&(SingularValue)},
    {"EffectivekBT","double",&(T)},
    {"Deflate","bool",&deflate},
    {NULL, NULL, NULL}};
  if(!(len==2)||(len==3)) return luaL_error(L, "Error in ChangeType expected 2 arguments and possible a third arguments specifing options, not %d\n",len);
  if(LuaTestOptionsList(L,-1)==1)
  {
    LuaCheckOptions(L, -1, LuaChangeTypeResponseFunctionOptions);
    len=len-1;
  }
  Min = luaL_checkResponsefunction(L,1);
  convert = luaL_checkstring(L, 2);
  if(strcmp(convert, "Dense")==0)
  {
    // Dense versions of response functions return 3 matrices, The response function, A0 and B0
    // check if we push table or userdata type
    if(Min->type<94)
    { // upercase push table on stack
      CompactMatrixType MCompactOut,A0,B0;
      ResponseFunctionToCompactMatrix(Min, &MCompactOut, deflate, SingularValue);
      LuaPushTable_CompactMatrixType(L, &MCompactOut);
      FreeCompactMatrix(&MCompactOut);
      if(Min->Block==true)
      {
        BlockResponseFunctionGetA0(Min, &A0);
        LuaPushTable_CompactMatrixType(L, &A0);
        FreeCompactMatrix(&A0);
        BlockResponseFunctionGetB0(Min, &B0);
        LuaPushTable_CompactMatrixType(L, &B0);
        FreeCompactMatrix(&B0);
        FreeResponsefunction(Min);
      }
      else
      {
        double a0,b0;
        ResponseFunctionGetA0(Min, &a0);
        ResponseFunctionGetB0(Min, &b0);
        lua_pushnumber(L, a0);
        lua_pushnumber(L, b0);
      }
    }
    else
    { // lowercase push userdata on stack
      CompactMatrixType *MCompactOut,*A0,*B0;
      MCompactOut = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
      luaL_setmetatable(L, "Matrix_Type");
      ResponseFunctionToCompactMatrix(Min, MCompactOut, deflate, SingularValue);
      MCompactOut->OnLuaStack=true;
      if(Min->Block==true)
      {
        A0 = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
        luaL_setmetatable(L, "Matrix_Type");
        BlockResponseFunctionGetA0(Min, A0);
        A0->OnLuaStack=true;
        B0 = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
        luaL_setmetatable(L, "Matrix_Type");
        BlockResponseFunctionGetB0(Min, B0);
        B0->OnLuaStack=true;
      }
      else
      {
        double a0,b0;
        ResponseFunctionGetA0(Min, &a0);
        ResponseFunctionGetB0(Min, &b0);
        lua_pushnumber(L, a0);
        lua_pushnumber(L, b0);
      }
    }
    return 3;
  }
  else if(strcmp(convert, "Tri")==0)
  {
    type='T'|(32&Min->type);  // t = T & 32, i.e. keep capital or not from input, but change letter to t/T
  }
  else if(strcmp(convert, "And")==0)
  {
    type='A'|(32&Min->type);
  }
  else if((strcmp(convert, "ListOfPoles")==0)||(strcmp(convert, "LP")==0))
  {
    type='L'|(32&Min->type);
  }
  else if((strcmp(convert, "NaturalImpurityOrbital")==0)||(strcmp(convert, "Nat")==0))
  {
    type='N'|(32&Min->type);
  }
  else if(strcmp(convert, "Userdata")==0)
  {
    type=(32)|Min->type; // to lower case
  }
  else if(strcmp(convert, "Table")==0)
  {
    type=(255-32)&Min->type; // to upper case
  }
  else
  {
    return luaL_error(L, "Error in ChangeType, second argument needs to be a string containgin either: \"And\", \"Tri\", \"ListOfPoles\", \"NaturalImpurityOrbital\", \"Nat\", \"Userdata\", \"Table\" and not %s\n",convert);
  }
  ConvertResponsefunction(Min, &Mout, type, T, SingularValue);
  LuaPush_ResponseFunction(L, &Mout);
  FreeResponsefunction(&Mout);
  FreeResponsefunction(Min);
  return 1;
}

static int LuaResponseFunctionToUserdata(lua_State *L)
{
  ResponsefunctionType *Min,Mout;
  uint32_t len = (uint32_t)lua_gettop(L);
  if(len!=1) return luaL_error(L, "Error in ToUserdata expected a single argument, not %d\n",len);
  Min = luaL_checkResponsefunction(L,1);
  if((Min->type|32)==Min->type) //lowercase, copy the response function and push it to the stack
  {
    ResponsefunctionCopy(Min,&Mout);
    LuaPush_ResponseFunction(L, &Mout);
    FreeResponsefunction(&Mout);
    FreeResponsefunction(Min);
  }
  else // upercase change to lower and push a table to the stack
  {
    Min->type+=32;
    LuaPush_ResponseFunction(L, Min);
    free(Min);
  }
  return 1;
}


static int LuaResponseFunctionToTable(lua_State *L)
{
  ResponsefunctionType *M;
  uint32_t len = (uint32_t)lua_gettop(L);
  if(len!=1) return luaL_error(L, "Error in ToTable expected a single argument, not %d\n",len);
  M = luaL_checkResponsefunction(L,1);
  switch(M->type)
  {
    case 'L':
    case 'l':
    {
      if(M->Block==true)
      {
        LuaPushTable_BlockListOfPoles(L, (BlockListOfPoles*)M->data);
      }
      else
      {
        LuaPushTable_ListOfPoles(L, (ListOfPoles*)(M->data));
      }
      break;
    }
    case 'A':
    case 'a':
    {
      if(M->Block==true)
      {
        LuaPushTable_BlockAndersonMatrixType(L, M->data);
      }
      else
      {
        LuaPushTable_AndersonMatrixType(L, M->data);
      }
      break;
    }
    case 'T':
    case 't':
    {
      if(M->Block==true)
      {
        LuaPushTable_BlockTriDiagonalMatrixType(L, M->data);
      }
      else
      {
        LuaPushTable_TriDiagonalMatrixType(L, M->data);
      }
      break;
    }
    case 'N':
    case 'n':
    {
      if(M->Block==true)
      {
        LuaPushTable_BlockNaturalImpurityOrbital(L, M->data);
      }
      else
      {
        LuaPushTable_NaturalImpurityOrbital(L, M->data);
      }
      break;
    }
    default:
    {
      printf("Error in ResponseFunctionToTable: unrecognized type %c",M->type);
      return 0;
    }
  }
  FreeResponsefunction(M);
  return 1;
}

/*
 *
 */
static int LuaReducePolesResponseFunction(lua_State *L)
{
  ResponsefunctionType *M, MRed;
  int Nred, Nredcon=0;

  uint32_t len = (uint32_t)lua_gettop(L);
  if(len>3) return luaL_error(L, "Error in ReducePoles: %d argument(s), 2 or 3 expected (response function, desired number of poles or ResponseFunction, desiered number of poles for val, desired number of poles for con).\n", len);

  M = luaL_checkResponsefunction(L,1);
  if(!lua_isnumber(L, 2)) return luaL_error(L, "Error in ReducePoles: 2nd argument should be the desired number of poles.\n");
  Nred = (int)lua_tointeger(L, 2);
  if (len==3)
  {
    if(!lua_isnumber(L, 3)) return luaL_error(L, "Error in ReducePoles: 3rd argument should be the desired number of poles in the conduction band.\n");
    Nredcon = (int)lua_tointeger(L, 3);
  }
  ConvertResponsefunction(M, &MRed, 'l', 0, 0);
  if(len==2)
  {
    ReduceBlockListOfPoles((BlockListOfPoles*)(MRed.data), Nred);
  }
  else
  {
    ReduceBlockListOfPoles((BlockListOfPoles*)(MRed.data), Nred, Nredcon);
  }
  LuaPush_ResponseFunction(L, &MRed);
  FreeResponsefunction(M);
  FreeResponsefunction(&MRed);
  return 1;
}

/*
 *
 */
static int LuaPlaceOnMeshResponseFunction(lua_State *L)
{
  ResponsefunctionType *M, MRed;
  uint32_t NE;
  double *mesh;

  uint32_t len = (uint32_t)lua_gettop(L);
  if(len!=2) return luaL_error(L, "Error in PlaceOnMesh: %d argument(s), 2 expected (ListOfPoles, desired list of grid points).\n", len);
  if(lua_type(L,2)!=LUA_TTABLE) return luaL_error(L, "Error in PlaceOnMesh: the 2nd argument should be a list of energy grid points.\n");

  M = luaL_checkResponsefunction(L,1);
  if((M->type|32)!='l')
  {
    return luaL_error(L, "Error in PlaceOnMesh: The first argument should be of type 'L'/'l' (ListOfPoles).\n");
  }
  mesh = luaL_checkarray_double(L, 2, &NE);
  ResponsefunctionCopy(M,&MRed);
  BlockListOfPolesPlaceOnMesh((BlockListOfPoles*)MRed.data, mesh, NE);
  LuaPush_ResponseFunction(L, &MRed);
  FreeResponsefunction(M);
  FreeResponsefunction(&MRed);
  return 1;
}

/*
 *
 */
static int LuaResponseFunctionEigensystem(lua_State *L)
{
  ResponsefunctionType *M;
  double *val;
  CompactMatrixType fun;
  M = luaL_checkResponsefunction(L,-1);
  
  switch(M->type)
  {
    case 'L':
    case 'l':
    {
      if(M->Block==true)
      {
        val = (double*)malloc(((BlockListOfPoles*)M->data)->N*((BlockListOfPoles*)M->data)->BlockSize*sizeof(double));
        BlockListOfPolesEigenSystem((BlockListOfPoles*)M->data, val, &fun);
      }
      else
      {
        val = (double*)malloc(((ListOfPoles*)M->data)->N*sizeof(double));
        ListOfPolesToEigensystemCompact((ListOfPoles*)(M->data), val, &fun);
      }
      break;
    }
    case 'A':
    case 'a':
    {
      if(M->Block==true)
      {
        val = (double*)malloc(((BlockAndersonMatrixType*)M->data)->N*((BlockAndersonMatrixType*)M->data)->BlockSize*sizeof(double));
        BlockAndersonMatrixEigenSystem((BlockAndersonMatrixType*)M->data, val, &fun);
      }
      else
      {
        val = (double*)malloc(((AndersonMatrixType*)M->data)->N*sizeof(double));
        AndersonMatrixToEigenSystemCompact((AndersonMatrixType*)M->data, val, &fun);
      }
      break;
    }
    case 'T':
    case 't':
    {
      if(M->Block==true)
      {
        val = (double*)malloc(((BlockTriDiagonalMatrixType*)M->data)->N*((BlockTriDiagonalMatrixType*)M->data)->BlockSize*sizeof(double));
        BlockTriDiagonalMatrixEigenSystem((BlockTriDiagonalMatrixType*)M->data, val, &fun);
      }
      else
      {
        val = (double*)malloc(((TriDiagonalMatrixType*)M->data)->N*sizeof(double));
        TriDiagonalMatrixToEigensystemCompact((TriDiagonalMatrixType*)M->data, val, &fun);
      }
      break;
    }
    case 'N':
    case 'n':
    {
      if(M->Block==true)
      {
        val = (double*)malloc((((BlockNaturalImpurityOrbitalMatrixType*)M->data)->val.N+ ((BlockNaturalImpurityOrbitalMatrixType*)M->data)->con.N )*((BlockNaturalImpurityOrbitalMatrixType*)M->data)->BlockSize*sizeof(double));
        BlockNaturalImpurityOrbitalMatrixEigenSystem((BlockNaturalImpurityOrbitalMatrixType*)M->data, val, &fun);
      }
      else
      {
        val = (double*)malloc((((BADoubleTriDiagonalMatrixType*)M->data)->val.N+((BADoubleTriDiagonalMatrixType*)M->data)->con.N+1)*sizeof(double));
        BADoubleTriDiagonalMatrixTypeToEigenSystemCompact((BADoubleTriDiagonalMatrixType*)M->data, val, &fun);
      }
      break;
    }
    default:
    {
      printf("Error in ResponseFunctionToTable: unrecognized type %c",M->type);
      return 0;
    }
  }
  LuaPushArray_double(L, val, fun.Ni);
  LuaPushTable_CompactMatrixType(L, &fun);
  FreeResponsefunction(M);
  FreeCompactMatrix(&fun);
  free(val);
  return 2;
}

/*  on stack two response functions and possible options
 *  The response functions must have the same dimension
 *  the First is G0, the second G
 *  calculated is G0^-1 - G^-1
 *  We transform G0 and G to an anderson form
 *  Inverting gives a list of poles
 *  we subtract list of poles and return sigma as a list of poles
 */
static int LuaCalculateSelfEnergy(lua_State *L)
{
  ResponsefunctionType *G0,*G,*Sigma;
  double SingularValue,T;
  LuaOptionsType CalculateSelfEnergyOptions[]={
    {"SingularValue","double",&(SingularValue)},
    {"EffectivekBT","double",&(T)},
    {NULL, NULL, NULL}};
  uint32_t len = (uint32_t)lua_gettop(L);
  if( LuaTestOptionsList(L,-1) )
  {
    LuaCheckOptions(L,-1,CalculateSelfEnergyOptions);
    len--;
  }
  if(len!=2)
  {
    return luaL_error(L, "CalculateSelfEnergy should be called with G0 and G as the argument (and an optional 3rd argument with options)\n Called with %d arguments\n",len);
  }
  G0 = luaL_checkResponsefunction(L,1);
  G  = luaL_checkResponsefunction(L,2);
  Sigma = (ResponsefunctionType *) lua_newuserdata(L, sizeof(ResponsefunctionType));
  CalculateSelfenergy(G0, G, Sigma);
  luaL_setmetatable(L, "ResponseFunction_Type");
  FreeResponsefunction(G0);
  FreeResponsefunction(G);
  return 1;
}

/*   Given Delta_0(omega) as a response function
 *   on a list of poles form it represents a tight binding Hamiltonian with:
 *      a_0 the onsite energy of the impurity
 *      b_i the hopping to a bath site with energy a_i
 *   Sigma(omega)
 *
 *   Calculate Delta(oemga) = Delta_0( omega - Sigma(omega) )
 */
static int LuaCalculateHybridizationFunction(lua_State *L)
{
  ResponsefunctionType *D0,*D,*Sigma;
  CompactVectorType Grid={0,0,NULL};
  double SingularValue=0,T=0;
  LuaOptionsType CalculateSelfEnergyOptions[]={
    {"SingularValue","double",&(SingularValue)},
    {"EffectivekBT","double",&(T)},
    {"EnergyGrid","Vector",&Grid},
    {NULL, NULL, NULL}};
  uint32_t len = (uint32_t)lua_gettop(L);
  if( LuaTestOptionsList(L,-1) )
  {
    LuaCheckOptions(L,-1,CalculateSelfEnergyOptions);
    len--;
  }
  if(len!=2)
  {
    return luaL_error(L, "CalculateHybridizationFunction should be called with D0 and Sigma as the argument (and an optional 3rd argument with options)\n Called with %d arguments\n",len);
  }
  D0    = luaL_checkResponsefunction(L,1);
  Sigma = luaL_checkResponsefunction(L,2);
  D     = (ResponsefunctionType *) lua_newuserdata(L, sizeof(ResponsefunctionType));
  if(Grid.Ni==0)
  {
    CalculateGBath(D0, Sigma, D);
  }
  else
  {
    CalculateGBathOnGrid(D0, Sigma, &Grid, D);
    free(Grid.VRe);
    Grid.Ni=0;
    Grid.QComplex=0;
  }
  luaL_setmetatable(L, "ResponseFunction_Type");
  FreeResponsefunction(D0);
  FreeResponsefunction(Sigma);
  return 1;
}

/*
 *
 */
static int LuaInvertEnergy(lua_State *L)
{
  ResponsefunctionType *Gin,Gout;
  Gin = luaL_checkResponsefunction(L,1);
  ResponseFunctionInvertEnergy(Gin,&Gout);
  LuaPush_ResponseFunction(L, &Gout);
  FreeResponsefunction(Gin);
  FreeResponsefunction(&Gout);
  return 1;
}


/*
 *
 */
static int LuaResponseFunctionUnm(lua_State * L)
{
  return luaL_error(L, "Responsefunctions can only store physical response functions. We know that Im[G]<=0 for all omega. For -G one would have Im[-G]>=0 and this would thus correspond to an unphisical response function\n");
}

/*
 *
 */
static int LuaResponseFunctionAdd(lua_State * L)
{
  ResponsefunctionType *Ma, *Mb, res;
  uint32_t len = (uint32_t)lua_gettop(L);
  if(len!=2) return luaL_error(L, "Error in ResponseFunctionAdd: %d argument(s), 2 expected.\n", len);
  Ma = luaL_checkResponsefunction(L,1);
  Mb = luaL_checkResponsefunction(L,2);
  if(Ma->Block != Mb->Block) return luaL_error(L, "Error in ResponseFunctionAdd: Block type of both functions needs to be the same\n");
  //TODO: if the functions are of Block type test dimensions
  
  ResponseFunctionAdd(Ma, Mb, &res);
  LuaPush_ResponseFunction(L, &res);
  FreeResponsefunction(Ma);
  FreeResponsefunction(Mb);
  FreeResponsefunction(&res);
  return 1;
}

/*
 *
 */
static int LuaResponseFunctionSub(lua_State * L)
{
  ResponsefunctionType *Ma, *Mb, res;
  uint32_t len = (uint32_t)lua_gettop(L);
  if(len!=2) return luaL_error(L, "Error in ResponseFunctionSub: %d argument(s), 2 expected.\n", len);
  Ma = luaL_checkResponsefunction(L,1);
  Mb = luaL_checkResponsefunction(L,2);
  if(Ma->Block != Mb->Block) return luaL_error(L, "Error in ResponseFunctionAdd: Block type of both functions needs to be the same\n");
  //TODO: if the functions are of Block type test dimensions
  
  ResponseFunctionSubtract(Ma, Mb, &res);
  LuaPush_ResponseFunction(L, &res);
  FreeResponsefunction(Ma);
  FreeResponsefunction(Mb);
  FreeResponsefunction(&res);
  return 1;
}

/*
 *
 */
static int LuaResponseFunctionMul(lua_State * L)
{
  //TODO: all
  printf("LuaResponseFunctionMul\n");
  return 0;
}

/*
 *
 */
static int LuaResponseFunctionDiv(lua_State * L)
{
  //TODO: all
  printf("LuaResponseFunctionDiv\n");
  return 0;
}


/*
 *
 */
static int LuaResponseFunctionCall(lua_State * L)
{
  ResponsefunctionType *res;
  double omega, gamma;
  if (lua_gettop(L)!=3)
  {
    return luaL_error(L, "Call to a response-function should contain two elements, (omega, gamma) with omega the energy at which the response function is evaluated and gamma an imaginary ofset from the real axis (FWMH lorenzian broadening)\n");
  }
  res = luaL_checkResponsefunction(L,1);
  omega = (double)luaL_checknumber(L, 2);
  gamma = (double)luaL_checknumber(L, 3);
  if(res->Block==true)
  {
    CompactMatrixType out;
    EvaluateBlockResponsefunctionType(res, &out, omega, gamma);
    LuaPushTable_CompactMatrixType(L, &out);
    FreeCompactMatrix(&out);  // added by Martin - otherwise heap is lost.
  }
  else
  {
    Complex16 *out;
    out = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
    *out = EvaluateResponsefunctionType(res, omega, gamma);
    luaL_setmetatable(L, "Complex_Type");
  }
  FreeResponsefunctionIncludingOriginalPointer(res);
  return 1;
}

/*
 *
 */
static int LuaResponseFunctionToString(lua_State * L)
{
  lua_pushfstring(L, "ResponseFunction in userdata format use ToTable() in order to get a table form");
  return 1;
}

/*  Given a Response function M
 *  push a table to the Lua stack containing 
 *     bitmasks for the different sections in the Hilbert space
 *     NF and NB as integers
 *     type as string
 *       bitmasks included are
 *         Core       CoreUp     CoreDn           -- a 0 string if FermionShift = 0
 *         Impurity   ImpurityUp ImpurityDn       -- Up and Dn only present if AddSpin = True
 *         BathVal    BathValUp  BathValDn        -- Not present if the format is Tri + Up and Dn only present if AddSpin = True
 *         BathCon    BathConUp  BathConDn        -- Not present if the format is Tri + Up and Dn only present if AddSpin = True
 *         BathExt    BathExtUp  BathExtDn        -- Only present if the format is Nat + Up and Dn only present if AddSpin = True
 *         Bath       BathUp     BathDn           -- Up and Dn only present if AddSpin = True
 *
 */
static int ResponseFunctionToLuaBitMask(lua_State *L, ResponsefunctionType *M, bool AddSpin, uint32_t FermionShift)
{
  int i,j,k;
  uint32_t NF,NB=0;
  uint32_t NCore, NImpurity, NBath, NBathVal, NBathCon, NBathExt;
  char *Core,*Impurity,*ImpurityUp,*ImpurityDn,*Bath,*BathUp,*BathDn;
  char *BathVal,*BathValUp,*BathValDn,*BathCon,*BathConDn,*BathConUp,*BathExt,*BathExtUp,*BathExtDn;
  
  lua_createtable(L, 0, 0);
  lua_pushboolean(L, M->Block);
  lua_setfield(L,-2,"Block");
  lua_pushboolean(L, AddSpin);
  lua_setfield(L,-2,"AddSpin");
  lua_pushinteger(L, FermionShift);
  lua_setfield(L,-2,"FermionShift");
  
  switch(M->type)
  {
    case 'L':
    case 'l':
    {
      return luaL_error(L, "Error in Responsefunction.CreateBitMasks() Can not create BitMaskFields for a response function stored as a list of poles\n");
      break;
    }
    case 'A':
    case 'a':
    {
      lua_pushstring(L, "And");
      lua_setfield(L,-2,"type");
      if(M->Block==true)
      {
        if(AddSpin==true)
        {
          NCore     = FermionShift;
          NImpurity = 2*((BlockAndersonMatrixType*)(M->data))->BlockSize;
          NBath     = NImpurity * (((BlockAndersonMatrixType*)(M->data))->N-1);
          NBathVal  = 0;
          NBathCon  = 0;
          NBathExt  = 0;
          NF = NCore + NImpurity + NBath;
          Core       = (char*)malloc((NF+1)*sizeof(char));
          Impurity   = (char*)malloc((NF+1)*sizeof(char));
          Bath       = (char*)malloc((NF+1)*sizeof(char));
          BathVal    = (char*)malloc((NF+1)*sizeof(char));
          BathCon    = (char*)malloc((NF+1)*sizeof(char));
          BathExt    = (char*)malloc((NF+1)*sizeof(char));
          ImpurityUp = (char*)malloc((NF+1)*sizeof(char));
          BathUp     = (char*)malloc((NF+1)*sizeof(char));
          BathValUp  = (char*)malloc((NF+1)*sizeof(char));
          BathConUp  = (char*)malloc((NF+1)*sizeof(char));
          BathExtUp  = (char*)malloc((NF+1)*sizeof(char));
          ImpurityDn = (char*)malloc((NF+1)*sizeof(char));
          BathDn     = (char*)malloc((NF+1)*sizeof(char));
          BathValDn  = (char*)malloc((NF+1)*sizeof(char));
          BathConDn  = (char*)malloc((NF+1)*sizeof(char));
          BathExtDn  = (char*)malloc((NF+1)*sizeof(char));
          for(k=0;k<NF;k++)
          {
            Core[k]       ='0';
            Impurity[k]   ='0';
            Bath[k]       ='0';
            BathVal[k]    ='0';
            BathCon[k]    ='0';
            BathExt[k]    ='0';
            ImpurityUp[k] ='0';
            BathUp[k]     ='0';
            BathValUp[k]  ='0';
            BathConUp[k]  ='0';
            BathExtUp[k]  ='0';
            ImpurityDn[k] ='0';
            BathDn[k]     ='0';
            BathValDn[k]  ='0';
            BathConDn[k]  ='0';
            BathExtDn[k]  ='0';
          }
          Core[NF]       ='\0';
          Impurity[NF]   ='\0';
          Bath[NF]       ='\0';
          BathVal[NF]    ='\0';
          BathCon[NF]    ='\0';
          BathExt[NF]    ='\0';
          ImpurityUp[NF] ='\0';
          BathUp[NF]     ='\0';
          BathValUp[NF]  ='\0';
          BathConUp[NF]  ='\0';
          BathExtUp[NF]  ='\0';
          ImpurityDn[NF] ='\0';
          BathDn[NF]     ='\0';
          BathValDn[NF]  ='\0';
          BathConDn[NF]  ='\0';
          BathExtDn[NF]  ='\0';
          for(k=0;k<NCore;k++)
          {
            Core[k]='1';
          }
          for(k=NCore;k<NCore+NImpurity;k++)
          {
            Impurity[k]='1'; ImpurityDn[k]='1'; ImpurityUp[k]='0';
            k++;
            Impurity[k]='1'; ImpurityDn[k]='0'; ImpurityUp[k]='1';
          }
          for(k=NCore+NImpurity;k<NF;k++)
          {
            Bath[k]='1'; BathDn[k]='1'; BathUp[k]='0';
            k++;
            Bath[k]='1'; BathDn[k]='0'; BathUp[k]='1';
          }
          k=NCore+NImpurity;
          if(((BlockAndersonMatrixType*)(M->data))->QComplex)
          {
            for(i=2;i<=((BlockAndersonMatrixType*)(M->data))->N;i++)
            {
              for(j=0;j<((BlockAndersonMatrixType*)(M->data))->BlockSize;j++)
              {
                if(((BlockAndersonMatrixType*)(M->data))->a[ i*((BlockAndersonMatrixType*)(M->data))->BlockSize*((BlockAndersonMatrixType*)(M->data))->BlockSize+j*((BlockAndersonMatrixType*)(M->data))->BlockSize+j ].Re < ((BlockAndersonMatrixType*)(M->data))->mu)
                {
                  NBathVal++;
                  BathVal[k]='1'; BathValDn[k]='1'; BathValUp[k]='0';
                  BathCon[k]='0'; BathConDn[k]='0'; BathConUp[k]='0';
                  BathExt[k]='0'; BathExtDn[k]='0'; BathExtUp[k]='0';
                  k++;
                  NBathVal++;
                  BathVal[k]='1'; BathValDn[k]='0'; BathValUp[k]='1';
                  BathCon[k]='0'; BathConDn[k]='0'; BathConUp[k]='0';
                  BathExt[k]='0'; BathExtDn[k]='0'; BathExtUp[k]='0';
                }
                else if(((BlockAndersonMatrixType*)(M->data))->a[ i*((BlockAndersonMatrixType*)(M->data))->BlockSize*((BlockAndersonMatrixType*)(M->data))->BlockSize+j*((BlockAndersonMatrixType*)(M->data))->BlockSize+j ].Re > ((BlockAndersonMatrixType*)(M->data))->mu)
                {
                  NBathCon++;
                  BathVal[k]='0'; BathValDn[k]='0'; BathValUp[k]='0';
                  BathCon[k]='1'; BathConDn[k]='1'; BathConUp[k]='0';
                  BathExt[k]='0'; BathExtDn[k]='0'; BathExtUp[k]='0';
                  k++;
                  NBathCon++;
                  BathVal[k]='0'; BathValDn[k]='0'; BathValUp[k]='0';
                  BathCon[k]='1'; BathConDn[k]='0'; BathConUp[k]='1';
                  BathExt[k]='0'; BathExtDn[k]='0'; BathExtUp[k]='0';
                }
                else
                {
                  NBathExt++;
                  BathVal[k]='0'; BathValDn[k]='0'; BathValUp[k]='0';
                  BathCon[k]='0'; BathConDn[k]='0'; BathConUp[k]='0';
                  BathExt[k]='1'; BathExtDn[k]='1'; BathExtUp[k]='0';
                  k++;
                  NBathExt++;
                  BathVal[k]='0'; BathValDn[k]='0'; BathValUp[k]='0';
                  BathCon[k]='0'; BathConDn[k]='0'; BathConUp[k]='0';
                  BathExt[k]='1'; BathExtDn[k]='0'; BathExtUp[k]='1';
                }
                k++;
              }
            }
          }
          else
          {
            for(i=2;i<=((BlockAndersonMatrixType*)(M->data))->N;i++)
            {
              for(j=0;j<((BlockAndersonMatrixType*)(M->data))->BlockSize;j++)
              {
                if(((BlockAndersonMatrixType*)(M->data))->aRe[ i*((BlockAndersonMatrixType*)(M->data))->BlockSize*((BlockAndersonMatrixType*)(M->data))->BlockSize+j*((BlockAndersonMatrixType*)(M->data))->BlockSize+j ] < ((BlockAndersonMatrixType*)(M->data))->mu)
                {
                  NBathVal++;
                  BathVal[k]='1'; BathValDn[k]='1'; BathValUp[k]='0';
                  BathCon[k]='0'; BathConDn[k]='0'; BathConUp[k]='0';
                  BathExt[k]='0'; BathExtDn[k]='0'; BathExtUp[k]='0';
                  k++;
                  NBathVal++;
                  BathVal[k]='1'; BathValDn[k]='0'; BathValUp[k]='1';
                  BathCon[k]='0'; BathConDn[k]='0'; BathConUp[k]='0';
                  BathExt[k]='0'; BathExtDn[k]='0'; BathExtUp[k]='0';
                }
                else if(((BlockAndersonMatrixType*)(M->data))->aRe[ i*((BlockAndersonMatrixType*)(M->data))->BlockSize*((BlockAndersonMatrixType*)(M->data))->BlockSize+j*((BlockAndersonMatrixType*)(M->data))->BlockSize+j ] > ((BlockAndersonMatrixType*)(M->data))->mu)
                {
                  NBathCon++;
                  BathVal[k]='0'; BathValDn[k]='0'; BathValUp[k]='0';
                  BathCon[k]='1'; BathConDn[k]='1'; BathConUp[k]='0';
                  BathExt[k]='0'; BathExtDn[k]='0'; BathExtUp[k]='0';
                  k++;
                  NBathCon++;
                  BathVal[k]='0'; BathValDn[k]='0'; BathValUp[k]='0';
                  BathCon[k]='1'; BathConDn[k]='0'; BathConUp[k]='1';
                  BathExt[k]='0'; BathExtDn[k]='0'; BathExtUp[k]='0';
                }
                else
                {
                  NBathExt++;
                  BathVal[k]='0'; BathValDn[k]='0'; BathValUp[k]='0';
                  BathCon[k]='0'; BathConDn[k]='0'; BathConUp[k]='0';
                  BathExt[k]='1'; BathExtDn[k]='1'; BathExtUp[k]='0';
                  k++;
                  NBathExt++;
                  BathVal[k]='0'; BathValDn[k]='0'; BathValUp[k]='0';
                  BathCon[k]='0'; BathConDn[k]='0'; BathConUp[k]='0';
                  BathExt[k]='1'; BathExtDn[k]='0'; BathExtUp[k]='1';
                }
                k++;
              }
            }
          }
          lua_pushstring(L, ImpurityUp);
          lua_setfield(L,-2,"ImpurityUp");
          free(ImpurityUp);
          lua_pushstring(L, ImpurityDn);
          lua_setfield(L,-2,"ImpurityDn");
          free(ImpurityDn);
          lua_pushstring(L, BathUp);
          lua_setfield(L,-2,"BathUp");
          free(BathUp);
          lua_pushstring(L, BathDn);
          lua_setfield(L,-2,"BathDn");
          free(BathDn);
          lua_pushstring(L, BathConUp);
          lua_setfield(L,-2,"BathConUp");
          free(BathConUp);
          lua_pushstring(L, BathConDn);
          lua_setfield(L,-2,"BathConDn");
          free(BathConDn);
          lua_pushstring(L, BathValUp);
          lua_setfield(L,-2,"BathValUp");
          free(BathValUp);
          lua_pushstring(L, BathValDn);
          lua_setfield(L,-2,"BathValDn");
          free(BathValDn);
          lua_pushstring(L, BathExtUp);
          lua_setfield(L,-2,"BathExtUp");
          free(BathExtUp);
          lua_pushstring(L, BathExtDn);
          lua_setfield(L,-2,"BathExtDn");
          free(BathExtDn);
        }
        else // AddSpin == false
        {
          NCore     = FermionShift;
          NImpurity = ((BlockAndersonMatrixType*)(M->data))->BlockSize;
          NBath     = NImpurity * (((BlockAndersonMatrixType*)(M->data))->N-1);
          NBathVal  = 0;
          NBathCon  = 0;
          NBathExt  = 0;
          NF = NCore + NImpurity + NBath;
          Core     = (char*)malloc((NF+1)*sizeof(char));
          Impurity = (char*)malloc((NF+1)*sizeof(char));
          Bath     = (char*)malloc((NF+1)*sizeof(char));
          BathVal  = (char*)malloc((NF+1)*sizeof(char));
          BathCon  = (char*)malloc((NF+1)*sizeof(char));
          BathExt  = (char*)malloc((NF+1)*sizeof(char));
          for(k=0;k<NF;k++)
          {
            Core[k]     ='0';
            Impurity[k] ='0';
            Bath[k]     ='0';
            BathVal[k]  ='0';
            BathCon[k]  ='0';
            BathExt[k]  ='0';
          }
          Core[NF]     ='\0';
          Impurity[NF] ='\0';
          Bath[NF]     ='\0';
          BathVal[NF]  ='\0';
          BathCon[NF]  ='\0';
          BathExt[NF]  ='\0';
          for(k=0;k<NCore;k++)
          {
            Core[k]='1';
          }
          for(k=NCore;k<NCore+NImpurity;k++)
          {
            Impurity[k]='1';
          }
          for(k=NCore+NImpurity;k<NF;k++)
          {
            Bath[k]='1';
          }
          k=NCore+NImpurity;
          if(((BlockAndersonMatrixType*)(M->data))->QComplex)
          {
            for(i=2;i<=((BlockAndersonMatrixType*)(M->data))->N;i++)
            {
              for(j=0;j<((BlockAndersonMatrixType*)(M->data))->BlockSize;j++)
              {
                if(((BlockAndersonMatrixType*)(M->data))->a[ i*((BlockAndersonMatrixType*)(M->data))->BlockSize*((BlockAndersonMatrixType*)(M->data))->BlockSize+j*((BlockAndersonMatrixType*)(M->data))->BlockSize+j ].Re < ((BlockAndersonMatrixType*)(M->data))->mu)
                {
                  NBathVal++;
                  BathVal[k]='1';
                  BathCon[k]='0';
                  BathExt[k]='0';
                }
                else if(((BlockAndersonMatrixType*)(M->data))->a[ i*((BlockAndersonMatrixType*)(M->data))->BlockSize*((BlockAndersonMatrixType*)(M->data))->BlockSize+j*((BlockAndersonMatrixType*)(M->data))->BlockSize+j ].Re > ((BlockAndersonMatrixType*)(M->data))->mu)
                {
                  NBathCon++;
                  BathVal[k]='0';
                  BathCon[k]='1';
                  BathExt[k]='0';
                }
                else
                {
                  NBathExt++;
                  BathVal[k]='0';
                  BathCon[k]='0';
                  BathExt[k]='1';
                }
                k++;
              }
            }
          }
          else
          {
            for(i=2;i<=((BlockAndersonMatrixType*)(M->data))->N;i++)
            {
              for(j=0;j<((BlockAndersonMatrixType*)(M->data))->BlockSize;j++)
              {
                if(((BlockAndersonMatrixType*)(M->data))->aRe[ i*((BlockAndersonMatrixType*)(M->data))->BlockSize*((BlockAndersonMatrixType*)(M->data))->BlockSize+j*((BlockAndersonMatrixType*)(M->data))->BlockSize+j ] < ((BlockAndersonMatrixType*)(M->data))->mu)
                {
                  NBathVal++;
                  BathVal[k]='1';
                  BathCon[k]='0';
                  BathExt[k]='0';
                }
                else if(((BlockAndersonMatrixType*)(M->data))->aRe[ i*((BlockAndersonMatrixType*)(M->data))->BlockSize*((BlockAndersonMatrixType*)(M->data))->BlockSize+j*((BlockAndersonMatrixType*)(M->data))->BlockSize+j ] > ((BlockAndersonMatrixType*)(M->data))->mu)
                {
                  NBathCon++;
                  BathVal[k]='0';
                  BathCon[k]='1';
                  BathExt[k]='0';
                }
                else
                {
                  NBathExt++;
                  BathVal[k]='0';
                  BathCon[k]='0';
                  BathExt[k]='1';
                }
                k++;
              }
            }
          }
        }
        lua_pushinteger(L, NCore);
        lua_setfield(L,-2,"NCore");
        lua_pushstring(L, Core);
        lua_setfield(L,-2,"Core");
        free(Core);
        lua_pushinteger(L, NImpurity);
        lua_setfield(L,-2,"NImpurity");
        lua_pushstring(L, Impurity);
        lua_setfield(L,-2,"Impurity");
        free(Impurity);
        lua_pushinteger(L, NBath);
        lua_setfield(L,-2,"NBath");
        lua_pushstring(L, Bath);
        lua_setfield(L,-2,"Bath");
        free(Bath);
        lua_pushinteger(L, NBathVal);
        lua_setfield(L,-2,"NBathVal");
        lua_pushstring(L, BathVal);
        lua_setfield(L,-2,"BathVal");
        free(BathVal);
        lua_pushinteger(L, NBathCon);
        lua_setfield(L,-2,"NBathCon");
        lua_pushstring(L, BathCon);
        lua_setfield(L,-2,"BathCon");
        free(BathCon);
        lua_pushinteger(L, NBathExt);
        lua_setfield(L,-2,"NBathExt");
        lua_pushstring(L, BathExt);
        lua_setfield(L,-2,"BathExt");
        free(BathExt);
        lua_pushinteger(L, NF);
        lua_setfield(L,-2,"NF");
        lua_pushinteger(L, NB);
        lua_setfield(L,-2,"NB");
      }
      else
      {
        if(AddSpin==true)
        {
          NCore     = FermionShift;
          NImpurity = 2;
          NBath     = NImpurity * (((AndersonMatrixType*)(M->data))->N-1);
          NBathVal  = 0;
          NBathCon  = 0;
          NBathExt  = 0;
          NF = NCore + NImpurity + NBath;
          Core       = (char*)malloc((NF+1)*sizeof(char));
          Impurity   = (char*)malloc((NF+1)*sizeof(char));
          Bath       = (char*)malloc((NF+1)*sizeof(char));
          BathVal    = (char*)malloc((NF+1)*sizeof(char));
          BathCon    = (char*)malloc((NF+1)*sizeof(char));
          BathExt    = (char*)malloc((NF+1)*sizeof(char));
          ImpurityUp = (char*)malloc((NF+1)*sizeof(char));
          BathUp     = (char*)malloc((NF+1)*sizeof(char));
          BathValUp  = (char*)malloc((NF+1)*sizeof(char));
          BathConUp  = (char*)malloc((NF+1)*sizeof(char));
          BathExtUp  = (char*)malloc((NF+1)*sizeof(char));
          ImpurityDn = (char*)malloc((NF+1)*sizeof(char));
          BathDn     = (char*)malloc((NF+1)*sizeof(char));
          BathValDn  = (char*)malloc((NF+1)*sizeof(char));
          BathConDn  = (char*)malloc((NF+1)*sizeof(char));
          BathExtDn  = (char*)malloc((NF+1)*sizeof(char));
          for(k=0;k<NF;k++)
          {
            Core[k]       ='0';
            Impurity[k]   ='0';
            Bath[k]       ='0';
            BathVal[k]    ='0';
            BathCon[k]    ='0';
            BathExt[k]    ='0';
            ImpurityUp[k] ='0';
            BathUp[k]     ='0';
            BathValUp[k]  ='0';
            BathConUp[k]  ='0';
            BathExtUp[k]  ='0';
            ImpurityDn[k] ='0';
            BathDn[k]     ='0';
            BathValDn[k]  ='0';
            BathConDn[k]  ='0';
            BathExtDn[k]  ='0';
          }
          Core[NF]       ='\0';
          Impurity[NF]   ='\0';
          Bath[NF]       ='\0';
          BathVal[NF]    ='\0';
          BathCon[NF]    ='\0';
          BathExt[NF]    ='\0';
          ImpurityUp[NF] ='\0';
          BathUp[NF]     ='\0';
          BathValUp[NF]  ='\0';
          BathConUp[NF]  ='\0';
          BathExtUp[NF]  ='\0';
          ImpurityDn[NF] ='\0';
          BathDn[NF]     ='\0';
          BathValDn[NF]  ='\0';
          BathConDn[NF]  ='\0';
          BathExtDn[NF]  ='\0';
          for(k=0;k<NCore;k++)
          {
            Core[k]='1';
          }
          for(k=NCore;k<NCore+NImpurity;k++)
          {
            Impurity[k]='1'; ImpurityDn[k]='1'; ImpurityUp[k]='0';
            k++;
            Impurity[k]='1'; ImpurityDn[k]='0'; ImpurityUp[k]='1';
          }
          for(k=NCore+NImpurity;k<NF;k++)
          {
            Bath[k]='1'; BathDn[k]='1'; BathUp[k]='0';
            k++;
            Bath[k]='1'; BathDn[k]='0'; BathUp[k]='1';
          }
          k=NCore+NImpurity;
          for(i=2;i<=((AndersonMatrixType*)(M->data))->N;i++)
          {
            if(((AndersonMatrixType*)(M->data))->a[ i ] < ((AndersonMatrixType*)(M->data))->mu)
            {
              NBathVal++;
              BathVal[k]='1'; BathValDn[k]='1'; BathValUp[k]='0';
              BathCon[k]='0'; BathConDn[k]='0'; BathConUp[k]='0';
              BathExt[k]='0'; BathExtDn[k]='0'; BathExtUp[k]='0';
              k++;
              NBathVal++;
              BathVal[k]='1'; BathValDn[k]='0'; BathValUp[k]='1';
              BathCon[k]='0'; BathConDn[k]='0'; BathConUp[k]='0';
              BathExt[k]='0'; BathExtDn[k]='0'; BathExtUp[k]='0';
            }
            else if(((AndersonMatrixType*)(M->data))->a[ i ] > ((AndersonMatrixType*)(M->data))->mu)
            {
              NBathCon++;
              BathVal[k]='0'; BathValDn[k]='0'; BathValUp[k]='0';
              BathCon[k]='1'; BathConDn[k]='1'; BathConUp[k]='0';
              BathExt[k]='0'; BathExtDn[k]='0'; BathExtUp[k]='0';
              k++;
              NBathCon++;
              BathVal[k]='0'; BathValDn[k]='0'; BathValUp[k]='0';
              BathCon[k]='1'; BathConDn[k]='0'; BathConUp[k]='1';
              BathExt[k]='0'; BathExtDn[k]='0'; BathExtUp[k]='0';
            }
            else
            {
              NBathExt++;
              BathVal[k]='0'; BathValDn[k]='0'; BathValUp[k]='0';
              BathCon[k]='0'; BathConDn[k]='0'; BathConUp[k]='0';
              BathExt[k]='1'; BathExtDn[k]='1'; BathExtUp[k]='0';
              k++;
              NBathExt++;
              BathVal[k]='0'; BathValDn[k]='0'; BathValUp[k]='0';
              BathCon[k]='0'; BathConDn[k]='0'; BathConUp[k]='0';
              BathExt[k]='1'; BathExtDn[k]='0'; BathExtUp[k]='1';
            }
            k++;
          }
          lua_pushstring(L, ImpurityUp);
          lua_setfield(L,-2,"ImpurityUp");
          free(ImpurityUp);
          lua_pushstring(L, ImpurityDn);
          lua_setfield(L,-2,"ImpurityDn");
          free(ImpurityDn);
          lua_pushstring(L, BathUp);
          lua_setfield(L,-2,"BathUp");
          free(BathUp);
          lua_pushstring(L, BathDn);
          lua_setfield(L,-2,"BathDn");
          free(BathDn);
          lua_pushstring(L, BathConUp);
          lua_setfield(L,-2,"BathConUp");
          free(BathConUp);
          lua_pushstring(L, BathConDn);
          lua_setfield(L,-2,"BathConDn");
          free(BathConDn);
          lua_pushstring(L, BathValUp);
          lua_setfield(L,-2,"BathValUp");
          free(BathValUp);
          lua_pushstring(L, BathValDn);
          lua_setfield(L,-2,"BathValDn");
          free(BathValDn);
          lua_pushstring(L, BathExtUp);
          lua_setfield(L,-2,"BathExtUp");
          free(BathExtUp);
          lua_pushstring(L, BathExtDn);
          lua_setfield(L,-2,"BathExtDn");
          free(BathExtDn);
        }
        else // AddSpin == false
        {
          NCore     = FermionShift;
          NImpurity = 1;
          NBath     = NImpurity * (((AndersonMatrixType*)(M->data))->N-1);
          NBathVal  = 0;
          NBathCon  = 0;
          NBathExt  = 0;
          NF = NCore + NImpurity + NBath;
          Core     = (char*)malloc((NF+1)*sizeof(char));
          Impurity = (char*)malloc((NF+1)*sizeof(char));
          Bath     = (char*)malloc((NF+1)*sizeof(char));
          BathVal  = (char*)malloc((NF+1)*sizeof(char));
          BathCon  = (char*)malloc((NF+1)*sizeof(char));
          BathExt  = (char*)malloc((NF+1)*sizeof(char));
          for(k=0;k<NF;k++)
          {
            Core[k]     ='0';
            Impurity[k] ='0';
            Bath[k]     ='0';
            BathVal[k]  ='0';
            BathCon[k]  ='0';
            BathExt[k]  ='0';
          }
          Core[NF]     ='\0';
          Impurity[NF] ='\0';
          Bath[NF]     ='\0';
          BathVal[NF]  ='\0';
          BathCon[NF]  ='\0';
          BathExt[NF]  ='\0';
          for(k=0;k<NCore;k++)
          {
            Core[k]='1';
          }
          for(k=NCore;k<NCore+NImpurity;k++)
          {
            Impurity[k]='1';
          }
          for(k=NCore+NImpurity;k<NF;k++)
          {
            Bath[k]='1';
          }
          k=NCore+NImpurity;
          for(i=2;i<=((AndersonMatrixType*)(M->data))->N;i++)
          {
            if(((AndersonMatrixType*)(M->data))->a[ i ] < ((AndersonMatrixType*)(M->data))->mu)
            {
              NBathVal++;
              BathVal[k]='1';
              BathCon[k]='0';
              BathExt[k]='0';
            }
            else if(((AndersonMatrixType*)(M->data))->a[ i ] > ((AndersonMatrixType*)(M->data))->mu)
            {
              NBathCon++;
              BathVal[k]='0';
              BathCon[k]='1';
              BathExt[k]='0';
            }
            else
            {
              NBathExt++;
              BathVal[k]='0';
              BathCon[k]='0';
              BathExt[k]='1';
            }
            k++;
          }
        }
        lua_pushinteger(L, NCore);
        lua_setfield(L,-2,"NCore");
        lua_pushstring(L, Core);
        lua_setfield(L,-2,"Core");
        free(Core);
        lua_pushinteger(L, NImpurity);
        lua_setfield(L,-2,"NImpurity");
        lua_pushstring(L, Impurity);
        lua_setfield(L,-2,"Impurity");
        free(Impurity);
        lua_pushinteger(L, NBath);
        lua_setfield(L,-2,"NBath");
        lua_pushstring(L, Bath);
        lua_setfield(L,-2,"Bath");
        free(Bath);
        lua_pushinteger(L, NBathVal);
        lua_setfield(L,-2,"NBathVal");
        lua_pushstring(L, BathVal);
        lua_setfield(L,-2,"BathVal");
        free(BathVal);
        lua_pushinteger(L, NBathCon);
        lua_setfield(L,-2,"NBathCon");
        lua_pushstring(L, BathCon);
        lua_setfield(L,-2,"BathCon");
        free(BathCon);
        lua_pushinteger(L, NBathExt);
        lua_setfield(L,-2,"NBathExt");
        lua_pushstring(L, BathExt);
        lua_setfield(L,-2,"BathExt");
        free(BathExt);
        lua_pushinteger(L, NF);
        lua_setfield(L,-2,"NF");
        lua_pushinteger(L, NB);
        lua_setfield(L,-2,"NB");
      }
      break;
    }
    case 'T':
    case 't':
    {
      lua_pushstring(L, "Tri");
      lua_setfield(L,-2,"type");
      if(M->Block==true)
      {
        if(AddSpin==true)
        {
          NCore     = FermionShift;
          NImpurity = 2*((BlockTriDiagonalMatrixType*)(M->data))->BlockSize;
          NBath     = 0;
          for(k=2;k<((BlockTriDiagonalMatrixType*)(M->data))->N+1;k++)
          {
            NBath+=2*((BlockTriDiagonalMatrixType*)(M->data))->BlockSizeList[k];
          }
          NF = NCore + NImpurity + NBath;
          Core       = (char*)malloc((NF+1)*sizeof(char));
          Impurity   = (char*)malloc((NF+1)*sizeof(char));
          Bath       = (char*)malloc((NF+1)*sizeof(char));
          ImpurityUp = (char*)malloc((NF+1)*sizeof(char));
          BathUp     = (char*)malloc((NF+1)*sizeof(char));
          ImpurityDn = (char*)malloc((NF+1)*sizeof(char));
          BathDn     = (char*)malloc((NF+1)*sizeof(char));
          for(k=0;k<NF;k++)
          {
            Core[k]       ='0';
            Impurity[k]   ='0';
            Bath[k]       ='0';
            ImpurityUp[k] ='0';
            BathUp[k]     ='0';
            ImpurityDn[k] ='0';
            BathDn[k]     ='0';
          }
          Core[NF]       ='\0';
          Impurity[NF]   ='\0';
          Bath[NF]       ='\0';
          ImpurityUp[NF] ='\0';
          BathUp[NF]     ='\0';
          ImpurityDn[NF] ='\0';
          BathDn[NF]     ='\0';
          for(k=0;k<NCore;k++)
          {
            Core[k]='1';
          }
          for(k=NCore;k<NCore+NImpurity;k++)
          {
            Impurity[k]='1'; ImpurityDn[k]='1'; ImpurityUp[k]='0';
            k++;
            Impurity[k]='1'; ImpurityDn[k]='0'; ImpurityUp[k]='1';
          }
          for(k=NCore+NImpurity;k<NF;k++)
          {
            Bath[k]='1'; BathDn[k]='1'; BathUp[k]='0';
            k++;
            Bath[k]='1'; BathDn[k]='0'; BathUp[k]='1';
          }
          k=NCore+NImpurity;
          lua_pushstring(L, ImpurityUp);
          lua_setfield(L,-2,"ImpurityUp");
          free(ImpurityUp);
          lua_pushstring(L, ImpurityDn);
          lua_setfield(L,-2,"ImpurityDn");
          free(ImpurityDn);
          lua_pushstring(L, BathUp);
          lua_setfield(L,-2,"BathUp");
          free(BathUp);
          lua_pushstring(L, BathDn);
          lua_setfield(L,-2,"BathDn");
          free(BathDn);
        }
        else // AddSpin == false
        {
          NCore     = FermionShift;
          NImpurity = ((BlockTriDiagonalMatrixType*)(M->data))->BlockSize;
          NBath     = 0;
          for(k=2;k<((BlockTriDiagonalMatrixType*)(M->data))->N+1;k++)
          {
            NBath+=((BlockTriDiagonalMatrixType*)(M->data))->BlockSizeList[k];
          }
          NBathVal  = 0;
          NBathCon  = 0;
          NBathExt  = 0;
          NF = NCore + NImpurity + NBath;
          Core     = (char*)malloc((NF+1)*sizeof(char));
          Impurity = (char*)malloc((NF+1)*sizeof(char));
          Bath     = (char*)malloc((NF+1)*sizeof(char));
          for(k=0;k<NF;k++)
          {
            Core[k]     ='0';
            Impurity[k] ='0';
            Bath[k]     ='0';
          }
          Core[NF]     ='\0';
          Impurity[NF] ='\0';
          Bath[NF]     ='\0';
          for(k=0;k<NCore;k++)
          {
            Core[k]='1';
          }
          for(k=NCore;k<NCore+NImpurity;k++)
          {
            Impurity[k]='1';
          }
          for(k=NCore+NImpurity;k<NF;k++)
          {
            Bath[k]='1';
          }
          k=NCore+NImpurity;
        }
        lua_pushinteger(L, NCore);
        lua_setfield(L,-2,"NCore");
        lua_pushstring(L, Core);
        lua_setfield(L,-2,"Core");
        free(Core);
        lua_pushinteger(L, NImpurity);
        lua_setfield(L,-2,"NImpurity");
        lua_pushstring(L, Impurity);
        lua_setfield(L,-2,"Impurity");
        free(Impurity);
        lua_pushinteger(L, NBath);
        lua_setfield(L,-2,"NBath");
        lua_pushstring(L, Bath);
        lua_setfield(L,-2,"Bath");
        free(Bath);
        lua_pushinteger(L, NF);
        lua_setfield(L,-2,"NF");
        lua_pushinteger(L, NB);
        lua_setfield(L,-2,"NB");
      }
      else
      {
        if(AddSpin==true)
        {
          NCore     = FermionShift;
          NImpurity = 2;
          NBath     = NImpurity * (((TriDiagonalMatrixType*)(M->data))->N-1);
          NBathVal  = 0;
          NBathCon  = 0;
          NBathExt  = 0;
          NF = NCore + NImpurity + NBath;
          Core       = (char*)malloc((NF+1)*sizeof(char));
          Impurity   = (char*)malloc((NF+1)*sizeof(char));
          Bath       = (char*)malloc((NF+1)*sizeof(char));
          ImpurityUp = (char*)malloc((NF+1)*sizeof(char));
          BathUp     = (char*)malloc((NF+1)*sizeof(char));
          ImpurityDn = (char*)malloc((NF+1)*sizeof(char));
          BathDn     = (char*)malloc((NF+1)*sizeof(char));
          for(k=0;k<NF;k++)
          {
            Core[k]       ='0';
            Impurity[k]   ='0';
            Bath[k]       ='0';
            ImpurityUp[k] ='0';
            BathUp[k]     ='0';
            ImpurityDn[k] ='0';
            BathDn[k]     ='0';
          }
          Core[NF]       ='\0';
          Impurity[NF]   ='\0';
          Bath[NF]       ='\0';
          ImpurityUp[NF] ='\0';
          BathUp[NF]     ='\0';
          ImpurityDn[NF] ='\0';
          BathDn[NF]     ='\0';
          for(k=0;k<NCore;k++)
          {
            Core[k]='1';
          }
          for(k=NCore;k<NCore+NImpurity;k++)
          {
            Impurity[k]='1'; ImpurityDn[k]='1'; ImpurityUp[k]='0';
            k++;
            Impurity[k]='1'; ImpurityDn[k]='0'; ImpurityUp[k]='1';
          }
          for(k=NCore+NImpurity;k<NF;k++)
          {
            Bath[k]='1'; BathDn[k]='1'; BathUp[k]='0';
            k++;
            Bath[k]='1'; BathDn[k]='0'; BathUp[k]='1';
          }
          k=NCore+NImpurity;
          lua_pushstring(L, ImpurityUp);
          lua_setfield(L,-2,"ImpurityUp");
          free(ImpurityUp);
          lua_pushstring(L, ImpurityDn);
          lua_setfield(L,-2,"ImpurityDn");
          free(ImpurityDn);
          lua_pushstring(L, BathUp);
          lua_setfield(L,-2,"BathUp");
          free(BathUp);
          lua_pushstring(L, BathDn);
          lua_setfield(L,-2,"BathDn");
          free(BathDn);
        }
        else // AddSpin == false
        {
          NCore     = FermionShift;
          NImpurity = 1;
          NBath     = NImpurity * (((TriDiagonalMatrixType*)(M->data))->N-1);
          NBathVal  = 0;
          NBathCon  = 0;
          NBathExt  = 0;
          NF = NCore + NImpurity + NBath;
          Core     = (char*)malloc((NF+1)*sizeof(char));
          Impurity = (char*)malloc((NF+1)*sizeof(char));
          Bath     = (char*)malloc((NF+1)*sizeof(char));
          for(k=0;k<NF;k++)
          {
            Core[k]     ='0';
            Impurity[k] ='0';
            Bath[k]     ='0';
          }
          Core[NF]     ='\0';
          Impurity[NF] ='\0';
          Bath[NF]     ='\0';
          for(k=0;k<NCore;k++)
          {
            Core[k]='1';
          }
          for(k=NCore;k<NCore+NImpurity;k++)
          {
            Impurity[k]='1';
          }
          for(k=NCore+NImpurity;k<NF;k++)
          {
            Bath[k]='1';
          }
          k=NCore+NImpurity;
        }
        lua_pushinteger(L, NCore);
        lua_setfield(L,-2,"NCore");
        lua_pushstring(L, Core);
        lua_setfield(L,-2,"Core");
        free(Core);
        lua_pushinteger(L, NImpurity);
        lua_setfield(L,-2,"NImpurity");
        lua_pushstring(L, Impurity);
        lua_setfield(L,-2,"Impurity");
        free(Impurity);
        lua_pushinteger(L, NBath);
        lua_setfield(L,-2,"NBath");
        lua_pushstring(L, Bath);
        lua_setfield(L,-2,"Bath");
        free(Bath);
        lua_pushinteger(L, NF);
        lua_setfield(L,-2,"NF");
        lua_pushinteger(L, NB);
        lua_setfield(L,-2,"NB");
      }
      break;
    }
    case 'N':
    case 'n':
    {
      lua_pushstring(L, "Nat");
      lua_setfield(L,-2,"type");
      if(M->Block==true)
      {
        if(AddSpin==true)
        {
          
        }
        else
        {
          
        }
      }
      else
      {
        if(AddSpin==true)
        {
          
        }
        else
        {
          
        }
      }
      break;
    }
    default:
    {
      return luaL_error(L, "Error in Responsefunction.CreateBitMasks() unrecognized type in response function (%s)\n",M->type);
    }
  }
  return 0;
}
static int LuaResponseFunctionToBitMask(lua_State *L)
{
  ResponsefunctionType *M;
  bool AddSpin=false;
  uint32_t FermionShift=0;
  LuaOptionsType LuaOptions[]={
    {"AddSpin","bool",&AddSpin},
    {"FermionShift","uint32_t",&FermionShift},
    {"NCore","uint32_t",&FermionShift},
    {NULL, NULL, NULL}};
  uint32_t len = (uint32_t) lua_gettop(L);
  if(LuaTestOptionsList(L, -1))
  {
    LuaCheckOptions(L, -1, LuaOptions);
    len--;
  }
  if(len!=1) return luaL_error(L, "Error in ResponseFunction.CreateBitMasks() It should be called with 1 argument, not with %d\n",len);
  M = luaL_checkResponsefunction(L,1);
  ResponseFunctionToLuaBitMask(L, M, AddSpin, FermionShift);
  FreeResponsefunctionIncludingOriginalPointer(M);
  return 1;
}

/*
 *
 */
static int LuaResponseFunctionToTightBinding(lua_State * L)
{
  TightBindingDef *TB;
  ResponsefunctionType *M;
  TB = (TightBindingDef*) lua_newuserdata(L, sizeof(TightBindingDef));
  luaL_setmetatable(L, "TightBinding_Type");
  M = luaL_checkResponsefunction(L,1);
  switch(M->type)
  {
    case 'L':
    case 'l':
    {
      return luaL_error(L, "Error in Responsefunction.ToTightBinding() Can not create Tight Binding representations for a response function stored as a list of poles\n");
      break;
    }
    case 'A':
    case 'a':
    {
      if(M->Block==true)
      {
        BlockAndersonMatrixToTightBinding((BlockAndersonMatrixType *)(M->data),TB);
      }
      else
      {
        //TODO: andersonmatrixToTightbinding
      }
      break;
    }
    case 'T':
    case 't':
    {
      if(M->Block==true)
      {
        BlockTriDiagonalMatrixToTightBinding((BlockTriDiagonalMatrixType *)(M->data),TB);
      }
      else
      {
        //TODO: tridiagonal matrix to TB
      }
      break;
    }
    case 'N':
    case 'n':
    {
      if(M->Block==true)
      {
        BlockNaturalImpurityOrbitalMatrixToTightBinding((BlockNaturalImpurityOrbitalMatrixType *)(M->data),TB);
      }
      else
      {
        //TODO: NatrualImpuritiOrbitalMatrixType to TB
      }
      break;
    }
    default:
    {
      return luaL_error(L, "Error in Responsefunction.ToTightbinding() unrecognized type for conversion (%s)\n",M->type);
      return 0;
    }
  }
  FreeResponsefunctionIncludingOriginalPointer(M);
  return 1;
}


/*
 *
 */
static int LuaResponseFunctionToOperator(lua_State * L)
{
  OperatorType *O;
  ResponsefunctionType *M;
  struct
  {
    bool AddSpin;
    uint32_t FermionShift;
  } Options;
  Options.AddSpin=false;
  Options.FermionShift=0;
  LuaOptionsType LuaOptions[]={
    {"AddSpin","bool",&(Options.AddSpin)},
    {"FermionShift","uint32_t",&(Options.FermionShift)},
    {NULL, NULL, NULL}};
  uint32_t len = (uint32_t) lua_gettop(L);
  if(LuaTestOptionsList(L, -1))
  {
    LuaCheckOptions(L, -1, LuaOptions);
    len--;
  }
  if(len!=1) return luaL_error(L, "Error in ResponseFunction.ToOperator() It should be called with 1 argument, not with %d\n",len);
  M = luaL_checkResponsefunction(L,1);
  O = (OperatorType*) lua_newuserdata(L, sizeof(OperatorType));
  luaL_setmetatable(L, "Operator_Type");

  
  switch(M->type)
  {
    case 'L':
    case 'l':
    {
      return luaL_error(L, "Error in Responsefunction.ToOperator() Can not create Operators for a response function stored as a list of poles\n");
      break;
    }
    case 'A':
    case 'a':
    {
      if(M->Block==true)
      {
        BlockAndersonMatrixToImpurityHamiltonian(O, (BlockAndersonMatrixType*)(M->data), Options.AddSpin, Options.FermionShift, 0);
      }
      else
      {
        
      }
      break;
    }
    case 'T':
    case 't':
    {
      if(M->Block==true)
      {
        BlockTridiagonalMatrixToImpurityHamiltonian(O, (BlockTriDiagonalMatrixType*)(M->data), Options.AddSpin, Options.FermionShift, 0);
      }
      else
      {
        
      }
      break;
    }
    case 'N':
    case 'n':
    {
      if(M->Block==true)
      {
        BlockNaturalImpurityOrbitalMatrixToImpurityHamiltonian(O, (BlockNaturalImpurityOrbitalMatrixType*)(M->data), Options.AddSpin, Options.FermionShift, 0);
      }
      else
      {
        
      }
      break;
    }
    default:
    {
      return luaL_error(L, "Error in Responsefunction.ToOperator() unrecognized type for conversion (%s)\n",M->type);
      return 0;
    }
  }
  FreeResponsefunctionIncludingOriginalPointer(M);
  return 1;
}


static int LuaResponseFunctionToCompactMatrix(lua_State *L)
{
  ResponsefunctionType *A;
  CompactMatrixType *B;
  struct {
    bool   Deflate;
    double SingularValue;
  } Options;
  Options.Deflate = false;
  Options.SingularValue = DBL_MIN;
  LuaOptionsType LuaOptions[]={
    {"Deflate",      "bool",  &(Options.Deflate)},
    {"SingularValue","double",&(Options.SingularValue)},
    {NULL, NULL, NULL}};
  if(LuaTestOptionsList(L, -1))
  {
    LuaCheckOptions(L, -1, LuaOptions);
  }
  B = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
  luaL_setmetatable(L, "Matrix_Type");
  B->OnLuaStack=true;
  A = luaL_checkResponsefunction(L,1);
  ResponseFunctionToCompactMatrix(A, B, Options.Deflate, Options.SingularValue);
  return 1;
}

static int LuaResponseFunctionToSpectra(lua_State *L)
{
  SpectraTypeArray *Spec;
  ResponsefunctionType *A;
  struct {
    double EnergyShift;
    double Emin;
    double Emax;
    double GammaG;
    double GammaL;
    int    NE;
  } Options;
  LuaOptionsType LuaOptions[]={
    {"EnergyShift","double",&(Options.EnergyShift)},
    {"Emin","double",&(Options.Emin)},
    {"Emax","double",&(Options.Emax)},
    {"NE","uint32_t",&(Options.NE)},
    {"GammaG","double",&(Options.GammaG)},
    {"GammaL","double",&(Options.GammaL)},
    {NULL, NULL, NULL}};
  Options.EnergyShift = 0;
  Options.Emin   =  DBL_MAX;
  Options.Emax   = -DBL_MAX;
  Options.GammaG = 0;
  Options.GammaL = 0;
  Options.NE     = 0;
  if(LuaTestOptionsList(L, -1))
  {
    LuaCheckOptions(L, -1, LuaOptions);
  }
  A = luaL_checkResponsefunction(L,1);
  Spec = (SpectraTypeArray *) lua_newuserdata(L, sizeof(SpectraTypeArray));
  luaL_setmetatable(L, "Spectra_Type");
  ResponseFunctionToSpectra(A,Spec, Options.EnergyShift, Options.Emin,Options.Emax,Options.GammaG,Options.GammaL,Options.NE);
  return 1;
}


/*
 *
 */
static int LuaResponseFunctionDestroy(lua_State * L)
{
  ResponsefunctionType *a = (ResponsefunctionType *) luaL_checkudata(L, 1, "ResponseFunction_Type");
  // convert type to captial letter, i.e. claim controle over the allocation
  a->type^=32;
  FreeResponsefunction(a);  
  return 0;
}

static const struct luaL_Reg ResponseFunction_func[] = {
  {"ChangeType", LuaChangeTypeResponseFunction},
  {"ReducePoles", LuaReducePolesResponseFunction},
  {"PlaceOnMesh", LuaPlaceOnMeshResponseFunction},
  {"ToTable", LuaResponseFunctionToTable},
  {"ToMatrix", LuaResponseFunctionToCompactMatrix},
  {"ToSpectra", LuaResponseFunctionToSpectra},
  {"ToUserdata", LuaResponseFunctionToUserdata},
  {"ToTightBinding", LuaResponseFunctionToTightBinding},
  {"ToOperator", LuaResponseFunctionToOperator},
  {"CalculateSelfEnergy", LuaCalculateSelfEnergy},
  {"CalculateHybridizationFunction", LuaCalculateHybridizationFunction},
  {"InvertEnergy", LuaInvertEnergy},
  {"CreateBitMasks", LuaResponseFunctionToBitMask},
  {NULL, NULL}
};

static const struct luaL_Reg ResponseFunction_meta[] = {
  {"__unm", LuaResponseFunctionUnm},
  {"__add", LuaResponseFunctionAdd},
  {"__sub", LuaResponseFunctionSub},
  {"__mul", LuaResponseFunctionMul},
  {"__div", LuaResponseFunctionDiv},
  {"__call", LuaResponseFunctionCall},
  {"__tostring", LuaResponseFunctionToString},
  {"__gc", LuaResponseFunctionDestroy},
  {NULL, NULL}
};

static const struct luaL_Reg ResponseFunctionTable_meta[] = {
  {"__unm", LuaResponseFunctionUnm},
  {"__add", LuaResponseFunctionAdd},
  {"__sub", LuaResponseFunctionSub},
  {"__mul", LuaResponseFunctionMul},
  {"__div", LuaResponseFunctionDiv},
  {"__call", LuaResponseFunctionCall},
  {NULL, NULL}
};


