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
 *   Tri, And, ListOfPoles, NaturalImpurityOrbital, DoubleTri, Userdata, Table, Dense
 *
 *   SingularValue is an option used to determine if the block size needs to be reduced in a tridiagonalization routine
 *   EffectivekBT is used to determine if a state belongs to the valence or conduction bands (occupied or empty using Fermi Dirac statistics)
 *   Deflate is a boolean option to determine if the final dimension should be the same (adding zero coupling) or is allowed to be reduced
 *   NTriMax is the largest chain length the result may have, NTriMaxVal and NTriMaxCon the same for the valence
 *     and the conduction chain of a NaturalImpurityOrbital or a DoubleTri matrix. Tri takes NTriMax, Nat and
 *     DoubleTri take all three, And and ListOfPoles have no chain and take none of them.
 */
static int LuaChangeTypeResponseFunction(lua_State *L)
{
  ResponsefunctionType *Min, Mout;
  char const * convert;
  char type;
  double ValueNotSet = -DBL_MAX;
  struct
  {
    char * Method;
    double SingularValue;
    double MinimalPoleWeigth;
    double MinimalPoleDistance;
    double T;
    bool deflate;
    bool relativeDistance;
    bool relativeWeigth;
    uint32_t NTriMax;
    uint32_t NTriMaxVal;
    uint32_t NTriMaxCon;
  } Options;
  // Standard values
  Options.Method = (char*)"Default";
  Options.SingularValue       = ValueNotSet;
  Options.MinimalPoleWeigth   = ValueNotSet;
  Options.MinimalPoleDistance = ValueNotSet;
  Options.T = ValueNotSet;
  Options.deflate=false;
  Options.relativeDistance=false;
  Options.relativeWeigth=false;
  Options.NTriMax    = NTriMaxNoLimit;
  Options.NTriMaxVal = NTriMaxNoLimit;
  Options.NTriMaxCon = NTriMaxNoLimit;
  uint32_t len = (uint32_t)lua_gettop(L);
  LuaOptionsType LuaChangeTypeResponseFunctionOptions[]={
    {"Method","char*",&(Options.Method)},
    {"SingularValue","double",&(Options.SingularValue)},
    {"MinimalPoleWeigth","double",&(Options.MinimalPoleWeigth)},
    {"MinimalPoleDistance","double",&(Options.MinimalPoleDistance)},
    {"EffectivekBT","double",&(Options.T)},
    {"Deflate","bool",&(Options.deflate)},
    {"relativeDistance","bool",&(Options.relativeDistance)},
    {"relativeWeigth","bool",&(Options.relativeWeigth)},
    {"NTriMax","uint32_t",&(Options.NTriMax)},
    {"NTriMaxVal","uint32_t",&(Options.NTriMaxVal)},
    {"NTriMaxCon","uint32_t",&(Options.NTriMaxCon)},
    {NULL, NULL, NULL}};
  if(LuaTestOptionsList(L,-1)==1)
  {
    LuaCheckOptions(L, -1, LuaChangeTypeResponseFunctionOptions);
    len=len-1;
  }
  if(len!=2) return luaL_error(L, "Error in ChangeType expected 2 arguments and possible a third arguments specifing options, not %d\n",len);
  /*  A relative pole distance is not implemented in BlockListOfPolesMergeDegeneratePoles
   *  yet, so refuse it here rather than let a conversion fail somewhere deep inside. */
  if(Options.relativeDistance==true) return luaL_error(L, "Error in ChangeType: relativeDistance is not implemented\n");
  if(Options.relativeWeigth==true) return luaL_error(L, "Error in ChangeType: relativeWeigth is not implemented\n");
  convert = luaL_checkstring(L, 2);
  /*  NTriMax, NTriMaxVal and NTriMaxCon put a ceiling on the length of the chain(s) of the
   *  result, a chain of length N being a[0..N] with b[0..N-1]. Only the representations
   *  that have a chain take them: Tri has one and takes NTriMax, Nat and DoubleTri have a
   *  valence and a conduction chain and take NTriMaxVal and NTriMaxCon, with NTriMax
   *  standing for both. A chain that comes out shorter on its own is left as it is, a
   *  longer one is cut without looking at the values further down. NTriMaxVal is what the
   *  single chain of a Tri travels on from here.
   */
  if(strcmp(convert, "Tri")==0)
  {
    if(Options.NTriMaxVal!=NTriMaxNoLimit) return luaL_error(L, "Error in ChangeType: a tri-diagonal matrix has one chain, its length is set with NTriMax and not with NTriMaxVal\n");
    if(Options.NTriMaxCon!=NTriMaxNoLimit) return luaL_error(L, "Error in ChangeType: a tri-diagonal matrix has one chain, its length is set with NTriMax and not with NTriMaxCon\n");
    Options.NTriMaxVal = Options.NTriMax;
  }
  else if((strcmp(convert, "NaturalImpurityOrbital")==0)||(strcmp(convert, "Nat")==0)||(strcmp(convert, "DoubleTri")==0)||(strcmp(convert, "DTri")==0))
  {
    if(Options.NTriMaxVal==NTriMaxNoLimit) Options.NTriMaxVal = Options.NTriMax;
    if(Options.NTriMaxCon==NTriMaxNoLimit) Options.NTriMaxCon = Options.NTriMax;
  }
  else if((Options.NTriMax!=NTriMaxNoLimit)||(Options.NTriMaxVal!=NTriMaxNoLimit)||(Options.NTriMaxCon!=NTriMaxNoLimit))
  {
    return luaL_error(L, "Error in ChangeType: a conversion to \"%s\" produces no chain and does not take NTriMax, NTriMaxVal or NTriMaxCon\n", convert);
  }
  Min = luaL_checkResponsefunction(L,1);
  // select the type and do the conversion
  if(strcmp(convert, "Dense")==0)
  {
    if(strcmp(Options.Method, "Default")==0)
    {
      // set standard values for options
      if(Options.SingularValue==ValueNotSet) Options.SingularValue=1000*DBL_EPSILON;
      // Dense versions of response functions return 3 matrices, The response function, A0 and B0
      // check if we push table or userdata type
      if(Min->type<94)
      { // upercase push table on stack
        CompactMatrixType MCompactOut,A0,B0;
        if(ResponseFunctionToCompactMatrix(Min, &MCompactOut, Options.deflate, Options.SingularValue))
        {
          return luaL_error(L, "Error in ChangeType(\"Dense\"): can not convert a response function of type %c to a dense matrix\n",Min->type);
        }
        LuaPushTable_CompactMatrixType(L, &MCompactOut);
        FreeCompactMatrix(&MCompactOut);
        if(Min->Block==true)
        {
          if(BlockResponseFunctionGetA0(Min, &A0))
          {
            return luaL_error(L, "Error in ChangeType(\"Dense\"): can not read A0 of a response function of type %c\n",Min->type);
          }
          LuaPushTable_CompactMatrixType(L, &A0);
          FreeCompactMatrix(&A0);
          if(BlockResponseFunctionGetB0(Min, &B0))
          {
            return luaL_error(L, "Error in ChangeType(\"Dense\"): can not read B0 of a response function of type %c\n",Min->type);
          }
          LuaPushTable_CompactMatrixType(L, &B0);
          FreeCompactMatrix(&B0);
          FreeResponsefunction(Min);
        }
        else
        {
          double a0,b0;
          if(ResponseFunctionGetA0(Min, &a0) || ResponseFunctionGetB0(Min, &b0))
          {
            return luaL_error(L, "Error in ChangeType(\"Dense\"): can not read a0 and b0 of a response function of type %c\n",Min->type);
          }
          lua_pushnumber(L, a0);
          lua_pushnumber(L, b0);
        }
      }
      else
      { // lowercase push userdata on stack
        CompactMatrixType *MCompactOut,*A0,*B0;
        MCompactOut = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
        luaL_setmetatable(L, "Matrix_Type");
        /*  Set the fresh user data to an empty matrix before it is filled. If the call below
         *  fails we leave through a lua error and the garbage collector will still run
         *  FreeCompactMatrix on this user data, which frees Name and V. With both NULL that is
         *  a no-op, on the uninitiallized struct it would free two arbitrary pointers.
         */
        SetCompactMatrixEmpty(MCompactOut);
        if(ResponseFunctionToCompactMatrix(Min, MCompactOut, Options.deflate, Options.SingularValue))
        {
          return luaL_error(L, "Error in ChangeType(\"Dense\"): can not convert a response function of type %c to a dense matrix\n",Min->type);
        }
        MCompactOut->OnLuaStack=true;
        if(Min->Block==true)
        {
          A0 = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
          luaL_setmetatable(L, "Matrix_Type");
          SetCompactMatrixEmpty(A0);
          if(BlockResponseFunctionGetA0(Min, A0))
          {
            return luaL_error(L, "Error in ChangeType(\"Dense\"): can not read A0 of a response function of type %c\n",Min->type);
          }
          A0->OnLuaStack=true;
          B0 = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
          luaL_setmetatable(L, "Matrix_Type");
          SetCompactMatrixEmpty(B0);
          if(BlockResponseFunctionGetB0(Min, B0))
          {
            return luaL_error(L, "Error in ChangeType(\"Dense\"): can not read B0 of a response function of type %c\n",Min->type);
          }
          B0->OnLuaStack=true;
        }
        else
        {
          double a0,b0;
          if(ResponseFunctionGetA0(Min, &a0) || ResponseFunctionGetB0(Min, &b0))
          {
            return luaL_error(L, "Error in ChangeType(\"Dense\"): can not read a0 and b0 of a response function of type %c\n",Min->type);
          }
          lua_pushnumber(L, a0);
          lua_pushnumber(L, b0);
        }
      }
      return 3;
    }
    else // options selection
    {
      return luaL_error(L, "Error in ChangeType Method \" %s \" not implemented for transformation to Dense",Options.Method);
    }
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
  else if((strcmp(convert, "DoubleTri")==0)||(strcmp(convert, "DTri")==0))
  {
    type='D'|(32&Min->type);
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
    return luaL_error(L, "Error in ChangeType, second argument needs to be a string containgin either: \"And\", \"Tri\", \"ListOfPoles\", \"NaturalImpurityOrbital\", \"Nat\", \"DoubleTri\", \"DTri\", \"Userdata\", \"Table\" and not %s\n",convert);
  }
  // test if options are valid and set standard values
  if(strcmp(convert, "Tri")==0)
  {
    if(strcmp(Options.Method, "Default")==0)
    {
      if(Options.SingularValue==ValueNotSet) Options.SingularValue=1000*DBL_EPSILON;;
      if(Options.MinimalPoleWeigth!=ValueNotSet) return luaL_error(L, "Error in ChangeType Method \" %s \" Does not suport MinimalPoleWeigth for conversion to Tri\n",Options.Method);

      if(Options.T!=ValueNotSet) return luaL_error(L, "Error in ChangeType Method \" %s \" Does not suport EffectivekBT for conversion to Tri\n",Options.Method);
    }
    else
    {
      return luaL_error(L, "Error in ChangeType Method \" %s \" not implemented for transformation to Tri\n",Options.Method);
    }
  }
  else if(strcmp(convert, "And")==0)
  {
    if(strcmp(Options.Method, "Default")==0)
    {
      if(Options.SingularValue==ValueNotSet) Options.SingularValue=1000*DBL_EPSILON;;
      if(Options.MinimalPoleWeigth!=ValueNotSet) return luaL_error(L, "Error in ChangeType Method \" %s \" Does not suport MinimalPoleWeigth for conversion to And\n",Options.Method);

      if(Options.T!=ValueNotSet) return luaL_error(L, "Error in ChangeType Method \" %s \" Does not suport EffectivekBT for conversion to And\n",Options.Method);
    }
    else
    {
      return luaL_error(L, "Error in ChangeType Method \" %s \" not implemented for transformation to And\n",Options.Method);
    }
  }
  else if((strcmp(convert, "ListOfPoles")==0)||(strcmp(convert, "LP")==0))
  {
    if(strcmp(Options.Method, "Default")==0)
    {
      if(Options.SingularValue==ValueNotSet) Options.SingularValue=1000*DBL_EPSILON;;
      if(Options.MinimalPoleWeigth!=ValueNotSet) return luaL_error(L, "Error in ChangeType Method \" %s \" Does not suport MinimalPoleWeigth for conversion to ListOfPoles\n",Options.Method);

      if(Options.T!=ValueNotSet) return luaL_error(L, "Error in ChangeType Method \" %s \" Does not suport EffectivekBT for conversion to ListOfPoles\n",Options.Method);
    }
    else
    {
      return luaL_error(L, "Error in ChangeType Method \" %s \" not implemented for transformation to ListOfPoles\n",Options.Method);
    }
  }
  else if((strcmp(convert, "NaturalImpurityOrbital")==0)||(strcmp(convert, "Nat")==0))
  {
    if(strcmp(Options.Method, "Default")==0)
    {
      if(Options.SingularValue==ValueNotSet) Options.SingularValue=1000*DBL_EPSILON;;
      if(Options.MinimalPoleWeigth!=ValueNotSet) return luaL_error(L, "Error in ChangeType Method \" %s \" Does not suport MinimalPoleWeigth\n",Options.Method);
      if(Options.T==ValueNotSet) Options.T = sqrt(DBL_MIN);
    }
    else
    {
      return luaL_error(L, "Error in ChangeType Method \" %s \" not implemented for transformation to NaturalImpurityOrbital\n",Options.Method);
    }
  }
  else if((strcmp(convert, "DoubleTri")==0)||(strcmp(convert, "DTri")==0))
  {
    if(strcmp(Options.Method, "Default")==0)
    {
      if(Options.SingularValue==ValueNotSet) Options.SingularValue=1000*DBL_EPSILON;;
      if(Options.MinimalPoleWeigth!=ValueNotSet) return luaL_error(L, "Error in ChangeType Method \" %s \" Does not suport MinimalPoleWeigth\n",Options.Method);
      if(Options.T==ValueNotSet) Options.T = sqrt(DBL_MIN);
    }
    else
    {
      return luaL_error(L, "Error in ChangeType Method \" %s \" not implemented for transformation to DoubleTri\n",Options.Method);
    }
  }
  else if(strcmp(convert, "Userdata")==0)
  {
    if(strcmp(Options.Method, "Default")==0)
    {
      if(Options.SingularValue!=ValueNotSet) return luaL_error(L, "Error in ChangeType Method \" %s \" Does not suport SingularValue\n",Options.Method);
      if(Options.MinimalPoleWeigth!=ValueNotSet) return luaL_error(L, "Error in ChangeType Method \" %s \" Does not suport MinimalPoleWeigth\n",Options.Method);
      if(Options.MinimalPoleDistance!=ValueNotSet) return luaL_error(L, "Error in ChangeType Method \" %s \" Does not suport MinimalPoleDistance\n",Options.Method);
      if(Options.T!=ValueNotSet) return luaL_error(L, "Error in ChangeType Method \" %s \" Does not suport EffectivekBT\n",Options.Method);
    }
    else
    {
      return luaL_error(L, "Error in ChangeType Method \" %s \" not implemented for transformation to Userdata\n",Options.Method);
    }
  }
  else if(strcmp(convert, "Table")==0)
  {
    if(strcmp(Options.Method, "Default")==0)
    {
      if(Options.SingularValue!=ValueNotSet) return luaL_error(L, "Error in ChangeType Method \" %s \" Does not suport SingularValue\n",Options.Method);
      if(Options.MinimalPoleWeigth!=ValueNotSet) return luaL_error(L, "Error in ChangeType Method \" %s \" Does not suport MinimalPoleWeigth\n",Options.Method);
      if(Options.MinimalPoleDistance!=ValueNotSet) return luaL_error(L, "Error in ChangeType Method \" %s \" Does not suport MinimalPoleDistance\n",Options.Method);
      if(Options.T!=ValueNotSet) return luaL_error(L, "Error in ChangeType Method \" %s \" Does not suport EffectivekBT\n",Options.Method);
    }
    else
    {
      return luaL_error(L, "Error in ChangeType Method \" %s \" not implemented for transformation to Table\n",Options.Method);
    }
  }
  
  if(ConvertResponsefunction(Min, &Mout, type, Options.Method, Options.T, Options.SingularValue, Options.MinimalPoleWeigth, Options.MinimalPoleDistance, Options.deflate, Options.relativeDistance, Options.relativeWeigth, Options.NTriMaxVal, Options.NTriMaxCon))
  {
    char typein=Min->type;
    FreeResponsefunction(Min);
    return luaL_error(L, "Error in ChangeType, the conversion of a response function of type %c to type %c is not implemented\n",typein,type);
  }
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
  else // upercase change to lower and push a pointer to the stack
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
    case 'D':
    case 'd':
    {
      if(M->Block==true)
      {
        LuaPushTable_BlockDoubleTriDiagonal(L, M->data);
      }
      else
      {
        LuaPushTable_DoubleTriDiagonal(L, M->data);
      }
      break;
    }
    default:
    {
      printf("Error in ResponseFunctionToTable: unrecognized type %c\n",M->type);
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
  ConvertResponsefunction(M, &MRed, 'l', "Default", 0, 0, 0, 0, false, false, false);
  if(M->Block==true)
  {
    if(len==2)
    {
      ReduceBlockListOfPoles((BlockListOfPoles*)(MRed.data), Nred);
    }
    else
    {
      ReduceBlockListOfPoles((BlockListOfPoles*)(MRed.data), Nred, Nredcon);
    }
  }
  else
  {
    if(len==2)
    {
      ReduceListOfPolesNumberFermi((ListOfPoles*)(MRed.data), Nred, ((ListOfPoles*)(MRed.data))->mu );
    }
    else
    {
      ReduceListOfPolesNumberFermi((ListOfPoles*)(MRed.data), Nred, Nredcon, ((ListOfPoles*)(MRed.data))->mu );
    }
  }
  LuaPush_ResponseFunction(L, &MRed);
  FreeResponsefunction(M);
  FreeResponsefunction(&MRed);
  return 1;
}

/*  On stack ResponseFunction, Emin, Emax
 *  Returns a ResponseFunction Whereby all poles in the window Emin to Emax are merged in a single pole
 */
static int LuaMergePolesInIntervalResponseFunction(lua_State *L)
{
  ResponsefunctionType *M, MRed;
  double Emin,Emax;
  
  
  uint32_t len = (uint32_t)lua_gettop(L);
  if(len!=3) return luaL_error(L, "Error in MergePolesInInterval: %d argument(s), 3 expected (response function, Emin, Emax.\n", len);

  M = luaL_checkResponsefunction(L,1);
  if(!lua_isnumber(L, 2)) return luaL_error(L, "Error in MergePolesInInterval: 2nd argument should be the Emin.\n");
  Emin=(double)luaL_checknumber(L, 2);
  if(!lua_isnumber(L, 3)) return luaL_error(L, "Error in MergePolesInInterval: 3rd argument should be the Emax.\n");
  Emax=(double)luaL_checknumber(L, 3);
  
  
  ConvertResponsefunction(M, &MRed, 'l', "Default", 0, 0, 0, 0, false, false, false);
  if(MRed.Block==true)
  {
    MergeBlockListOfPolesInInterval((BlockListOfPoles *)MRed.data, Emin,  Emax);
  }
  else
  {
    MergeListOfPolesInInterval((ListOfPoles *)MRed.data, Emin,  Emax);
  }
  LuaPush_ResponseFunction(L, &MRed);
  FreeResponsefunction(M);
  FreeResponsefunction(&MRed);
  return 1;
}


/*  On stack ResponseFunction, Emin, Emax
 *  Returns a ResponseFunction Whereby all poles in the window Emin to Emax are removed, i.e. their residue is set to zero
 */
static int LuaRemovePolesInIntervalResponseFunction(lua_State *L)
{
  ResponsefunctionType *M, MRed;
  double Emin,Emax;
  uint32_t len = (uint32_t)lua_gettop(L);
  if(len!=3) return luaL_error(L, "Error in RemovePolesInInterval: %d argument(s), 3 expected (response function, Emin, Emax.\n", len);

  M = luaL_checkResponsefunction(L,1);
  if(!lua_isnumber(L, 2)) return luaL_error(L, "Error in RemovePolesInInterval: 2nd argument should be the Emin.\n");
  Emin=(double)luaL_checknumber(L, 2);
  if(!lua_isnumber(L, 3)) return luaL_error(L, "Error in RemovePolesInInterval: 3rd argument should be the Emax.\n");
  Emax=(double)luaL_checknumber(L, 3);
  ConvertResponsefunction(M, &MRed, 'l', "Default", 0, 0, 0, 0, false, false, false);
  if(MRed.Block==true)
  {
    RemoveBlockListOfPolesInInterval((BlockListOfPoles *)MRed.data, Emin,  Emax);
  }
  else
  {
    RemoveListOfPolesInInterval((ListOfPoles *)MRed.data, Emin,  Emax);
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
  ResponsefunctionType *M, *MRed;
  uint32_t NE;
  double *mesh;

  uint32_t len = (uint32_t)lua_gettop(L);
  if(len!=2) return luaL_error(L, "Error in PlaceOnMesh: %d argument(s), 2 expected (ListOfPoles, desired list of grid points).\n", len);
  if(lua_type(L,2)!=LUA_TTABLE) return luaL_error(L, "Error in PlaceOnMesh: the 2nd argument should be a list of energy grid points.\n");

  M = luaL_checkResponsefunction(L,1);
  mesh = luaL_checkarray_double(L, 2, &NE);
  MRed = (ResponsefunctionType *) lua_newuserdata(L, sizeof(ResponsefunctionType));
  ConvertResponsefunction(M, MRed, 'l', "Default", 0, 0, 0, 0, false, false, false);
  if(MRed->Block==true)
  {
    BlockListOfPolesPlaceOnMesh((BlockListOfPoles*)MRed->data, mesh, NE);
  }
  else
  {
    ListOfPoles temp;
    ListOfPolesPlaceOnMesh((ListOfPoles*)MRed->data, &temp, mesh, NE);
    FreeListOfPoles((ListOfPoles*)MRed->data);
    *(ListOfPoles*)MRed->data = temp;
  }
  // give Lua the controle over the allocation
  MRed->type = MRed->type|32;
  luaL_setmetatable(L, "ResponseFunction_Type");
  FreeResponsefunction(M);
  return 1;
}

/*
 *
 */
/*  Pushes the eigenvalues and the eigenvectors of the matrix that stands in the denominator
 *  of G, which is the matrix ToMatrix returns: the poles for a list of poles, a[1..N] with
 *  its couplings for a tri-diagonal matrix, the observed site with both chains for a double
 *  tri-diagonal one. Its dimension is the number of levels the format stores, so a format
 *  that carries redundant levels -- a natural impurity orbital matrix one, a double
 *  tri-diagonal one two per bath block -- returns more eigenvalues than the number of poles
 *  of G, the extra ones belonging to levels that decouple.
 *
 *  A tri-diagonal matrix has its own LAPACK driver, banded in the block case, and keeps it:
 *  it costs O(N) and O(N BlockSize^2) memory against the O(N^2) of a dense matrix, and bath
 *  chains of a thousand sites are normal. Everything else goes through the dense matrix,
 *  which is the route the natural impurity orbital case already took. Four of the routines
 *  that used to be called here -- for a block list of poles, a block and a single valued
 *  Anderson matrix and a block natural impurity orbital matrix -- were empty stubs that
 *  returned without touching val or fun, and the uninitialised fun was then read and freed.
 */
static int LuaResponseFunctionEigensystem(lua_State *L)
{
  ResponsefunctionType *M;
  double *val;
  CompactMatrixType Mc, fun;
  M = luaL_checkResponsefunction(L,-1);
  if((M->type|32)=='t')
  {
    if(M->Block==true)
    {
      val = (double*)malloc(((BlockTriDiagonalMatrixType*)M->data)->N*((BlockTriDiagonalMatrixType*)M->data)->BlockSize*sizeof(double));
      if(val==NULL || BlockTriDiagonalMatrixEigenSystem((BlockTriDiagonalMatrixType*)M->data, val, &fun))
      {
        free(val);
        FreeResponsefunction(M);
        return luaL_error(L, "Error in Eigensystem: BlockTriDiagonalMatrixEigenSystem failed\n");
      }
    }
    else
    {
      val = (double*)malloc(((TriDiagonalMatrixType*)M->data)->N*sizeof(double));
      if(val==NULL || TriDiagonalMatrixToEigensystemCompact((TriDiagonalMatrixType*)M->data, val, &fun))
      {
        free(val);
        FreeResponsefunction(M);
        return luaL_error(L, "Error in Eigensystem: TriDiagonalMatrixToEigensystemCompact failed\n");
      }
    }
  }
  else
  {
    /*  the same Deflate and SingularValue that ToMatrix defaults to, so that the
     *  eigenvalues are those of the matrix ToMatrix hands out. Without deflation the
     *  matrix keeps the levels that a conversion left with zero coupling and their
     *  eigenvalues come along -- a block list of poles of three poles by two orbitals
     *  gives twelve of them where six are the poles of G. */
    if(ResponseFunctionToCompactMatrix(M, &Mc, true, DBL_EPSILON))
    {
      FreeResponsefunction(M);
      return luaL_error(L, "Error in Eigensystem: no matrix representation of a response function of type %c\n", M->type);
    }
    val = (double*)malloc(Mc.Ni*sizeof(double));
    if(val==NULL || CompactMatrixEigensystem(&Mc, &fun, val))
    {
      free(val);
      FreeCompactMatrix(&Mc);
      FreeResponsefunction(M);
      return luaL_error(L, "Error in Eigensystem: CompactMatrixEigensystem failed\n");
    }
    FreeCompactMatrix(&Mc);
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
 *
 *  epsilon is the minimal distance between two poles that belong to the same block in the anderson matrix representation to be considered different in energy
 */
static int LuaCalculateSelfEnergy(lua_State *L)
{
  ResponsefunctionType *G0,*G,*Sigma;
  double epsilon = 1000*DBL_EPSILON;
  LuaOptionsType CalculateSelfEnergyOptions[]={
    {"epsilon","double",&(epsilon)},
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
  CalculateSelfenergy(G0, G, Sigma, epsilon);
  luaL_setmetatable(L, "ResponseFunction_Type");
  FreeResponsefunction(G0);
  FreeResponsefunction(G);
  return 1;
}

/*   Given G_0(omega) and Sigma(omega) on the stack as response functons *
 *   Calculate 1 / ( G_0(w-Sigma(w))^{-1} + Sigma(w) )
 *   This yields the reponse function of a tight binding model whereby all sites get the self energy Sigma(omega)
 *   except for the site where we act on
 *
 *   Options allow to specify an energy grid for the final response function in Anderson representation
 */
static int LuaCalculateHybridizationFunction(lua_State *L)
{
  ResponsefunctionType *D0,*D,*Sigma;
  CompactVectorType Grid={0,0,NULL};
  double Epsilon=1000*DBL_EPSILON;
  LuaOptionsType CalculateSelfEnergyOptions[]={
    {"Epsilon","double",&(Epsilon)},
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
  if(D0->Block!=Sigma->Block) return luaL_error(L, "CalculateHybridizationFunction should be called with D0 and Sigma either both as block response functions or both as single valued response functions\n");
  if(BlockResponseFunctionGetBlockSize(D0) != BlockResponseFunctionGetBlockSize(Sigma)) return luaL_error(L, "CalculateHybridizationFunction should be called with D0 and Sigma with same dimension of BlockSize\n");
  D     = (ResponsefunctionType *) lua_newuserdata(L, sizeof(ResponsefunctionType));
  if(Grid.Ni==0)
  {
    CalculateGBath(D0, Sigma, D, Epsilon);
  }
  else
  {
    CalculateGBathOnGrid(D0, Sigma, &Grid, D, Epsilon);
    free(Grid.VRe);
    Grid.Ni=0;
    Grid.QComplex=0;
  }
  luaL_setmetatable(L, "ResponseFunction_Type");
  FreeResponsefunction(D0);
  FreeResponsefunction(Sigma);
  return 1;
}

static int LuaCalculateGBathFromGandSigma(lua_State *L)
{
  ResponsefunctionType *G,*Sigma,*GBath,*Gv=NULL,*Gc=NULL;
  CompactVectorType Grid={0,0,NULL};
  uint32_t Nv=UINT32_MAX, Nc=UINT32_MAX;
  double Epsilon=1000*DBL_EPSILON;
  char * GBathFormat=NULL;
  char * Method=NULL;
  LuaOptionsType CalculateybridizationFunctionOptions[]={
    {"Epsilon","double",&(Epsilon)},
    {"EnergyGrid","Vector",&Grid},
    {"Format","char*",&(GBathFormat)},
    {"Method","char*",&(Method)},
    {"Nc","uint32_t",&Nc},
    {"Nv","uint32_t",&Nv},
    {NULL, NULL, NULL}};
  uint32_t len = (uint32_t)lua_gettop(L);
  if( LuaTestOptionsList(L,-1) )
  {
    LuaCheckOptions(L,-1,CalculateybridizationFunctionOptions);
    len--;
  }
  G    = luaL_checkResponsefunction(L,1);
  Sigma = luaL_checkResponsefunction(L,2);
  if(BlockResponseFunctionGetBlockSize(G) != BlockResponseFunctionGetBlockSize(Sigma)) return luaL_error(L, "CalculateHybridizationFunction should be called with G and Sigma with same dimension of BlockSize\n");
  GBath = (ResponsefunctionType *) lua_newuserdata(L, sizeof(ResponsefunctionType));
  
  CalculateGBathFromGandSigma(G, Sigma, GBath, Epsilon, &Grid, Nv, Nc, Gv, Gc, Method, GBathFormat);
  // we give the controle of the allocation of the memory of Gout to Lua
  // turn the type to a lower case letter
  GBath->type = GBath->type | 32;
  luaL_setmetatable(L, "ResponseFunction_Type");
  FreeResponsefunction(G);
  FreeResponsefunction(Sigma);
  return 1;
}



/*   Given G_0(omega) and Sigma(omega) on the stack as response functons *
 *   Calculate G_0(omega-Sigma(omega))
 *   This yields the reponse function of a tight binding model whereby all sites get the self energy Sigma(omega)
 *
 *   Options allow to specify an energy grid for the final response function in ListOfPoles representation
 *
 *   We have implemented several methods using different approximations
 *     Not all methods are valid for all systems
 *     Some methods migth treat systems with off diagonal elements incorrectly
 *     Some mehtods might treat cases where we deflate the blocksize in a tri-diagonal form incorrectly
 */
static int LuaCalculateGfromG0andSigma(lua_State *L)
{
  ResponsefunctionType *G0,*G,*Sigma;
  CompactVectorType Grid={0,0,NULL};
  double Epsilon=1000*DBL_EPSILON;
  char* Method;
  Method = (char*)"Automatic";
  LuaOptionsType CalculateSelfEnergyOptions[]={
    {"Method","char*",&Method},
    {"Epsilon","double",&(Epsilon)},
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
  G0    = luaL_checkResponsefunction(L,1);
  Sigma = luaL_checkResponsefunction(L,2);
  if(G0->Block!=Sigma->Block) return luaL_error(L, "CalculateGfromG0andSigma should be called with G0 and Sigma either both as block response functions or both as single valued response functions\n");
  if(BlockResponseFunctionGetBlockSize(G0) != BlockResponseFunctionGetBlockSize(Sigma)) return luaL_error(L, "CalculateGfromG0andSigma should be called with G0 and Sigma with same dimension of BlockSize\n");
  G     = (ResponsefunctionType *) lua_newuserdata(L, sizeof(ResponsefunctionType));
  if(Grid.Ni==0)
  { // not on a grid
    if ((strcmp(Method,"TriDiagonalG0")==0) || ((Sigma->Block==true)&&(strcmp(Method,"Automatic")==0)))
    {
      // note that this method is not gennerally correct
      ResponseFunctionAddSelfEnergyDenseViaTri(G0, Sigma, G, Epsilon);
    }
    else if ((strcmp(Method,"ListOfPolesG0")==0) || ((Sigma->Block==false)&&(strcmp(Method,"Automatic")==0)))
    {
      // note that this method is not gennerally correct
      CalculateGfromG0andSigma(G0, Sigma, G, Epsilon);
    }
    else if (strcmp(Method,"AndG0")==0)
    {
      // note that this method is not gennerally correct
      ResponseFunctionAddSelfEnergyViaAndLP(G0, Sigma, G,Epsilon, NULL);
    }
    else
    {
      return luaL_error(L, "CalculateHybridizationFunction called with unkown method.\n");
    }
  }
  else
  { // on a grid
    if ((strcmp(Method,"TriDiagonalG0")==0) )
    {
      //TODO: memory problem
      // note that this method is not gennerally correct
      ResponseFunctionAddSelfEnergyOnGridViaG0Tri(G0, Sigma, &Grid, G, Epsilon);
    }
    else if ((strcmp(Method,"ListOfPolesG0")==0) || ((Sigma->Block==false)&&(strcmp(Method,"Automatic")==0)) || ((Sigma->Block==true)&&(strcmp(Method,"Automatic")==0)) )
    {
      // note that this method is not gennerally correct
      CalculateGfromG0andSigmaOnGrid(G0, Sigma, &Grid, G, Epsilon);
    }
    else if (strcmp(Method,"AndG0")==0)
    {
      // note that this method is not gennerally correct
      ResponseFunctionAddSelfEnergyViaAndLP(G0, Sigma, G,Epsilon, &Grid);
    }
    else
    {
      return luaL_error(L, "CalculateHybridizationFunction called with unkown method.\n");
    }
    free(Grid.VRe);
    Grid.Ni=0;
    Grid.QComplex=0;
  }
  luaL_setmetatable(L, "ResponseFunction_Type");
  FreeResponsefunction(G0);
  FreeResponsefunction(Sigma);
  return 1;
}
                                       





/*  G(mu+w) -> -G(mu-w)^T, the transpose alone and mu - w the plain reflection of the
 *  argument, nothing conjugated. On the stored numbers that is A0 -> -A0^T for the constant
 *  term, e_n -> 2 mu - e_n for the poles and w_n -> w_n^T for the residues. Why it is ^T
 *  and not ^+, and why the whole of G and not only its pole part follows, stands in the
 *  header of ResponseFunctionInvertEnergy.
 */
static int LuaResponseFunctionInvertEnergy(lua_State *L)
{
  ResponsefunctionType *Gin,Gout;
  Gin = luaL_checkResponsefunction(L,1);
  if(ResponseFunctionInvertEnergy(Gin,&Gout))
  {
    FreeResponsefunction(Gin);
    FreeResponsefunction(&Gout);
    return luaL_error(L, "Error inverting the energy of a response function\n");
  }
  LuaPush_ResponseFunction(L, &Gout);
  FreeResponsefunction(Gin);
  FreeResponsefunction(&Gout);
  return 1;
}



/*  G(w) -> G(w)^T
 *
 */
static int LuaResponseFunctionTranspose(lua_State *L)
{
  ResponsefunctionType *Gin,Gout;
  Gin = luaL_checkResponsefunction(L,1);
  if(ResponseFunctionTranspose(Gin,&Gout))
  {
    FreeResponsefunction(Gin);
    FreeResponsefunction(&Gout);
    return luaL_error(L, "Error transposing a response function\n");
  }
  LuaPush_ResponseFunction(L, &Gout);
  FreeResponsefunction(Gin);
  FreeResponsefunction(&Gout);
  return 1;
}


/*  When A0=0 and B0 = 1 this returns w - G(w)^{-1}
 *  G = A0 + B0^* 1 / (w - A1 - B1 1/(...) B1^+ ) B0^T
 *  Return (B0^T)^{-1} ( A1 - B1 1/(...) B1^+ ) (B0^*)^{-1}
 */
static int LuaResponseFunctionEffectiveHamiltonian(lua_State *L)
{
  CompactMatrixType *A0,*B0;
  ResponsefunctionType *Gin,*Heff;
  double ValueNotSet = -DBL_MAX;
  struct
  {
    double SingularValue;
    double MinimalPoleWeigth;
    double MinimalPoleDistance;
    double T;
  } Options;
  Options.SingularValue       = ValueNotSet;
  Options.MinimalPoleWeigth   = ValueNotSet;
  Options.MinimalPoleDistance = ValueNotSet;
  Options.T                   = ValueNotSet;
  /*  only SingularValue and MinimalPoleDistance reach anything: l -> a deflates on the
   *  first, d -> a and the merge of the bath poles use the second. T and MinimalPoleWeigth
   *  are read by nothing, so they are not offered. */
  LuaOptionsType LuaOptions[]={
    {"SingularValue","double",&(Options.SingularValue)},
//    {"MinimalPoleWeigth","double",&(Options.MinimalPoleWeigth)},
    {"MinimalPoleDistance","double",&(Options.MinimalPoleDistance)},
//    {"Temperature","double",&(Options.T)},
//    {"EffectivekBT","double",&(Options.T)},
    {NULL, NULL, NULL}};
  uint32_t len = (uint32_t) lua_gettop(L);
  if(LuaTestOptionsList(L, -1))
  {
    LuaCheckOptions(L, -1, LuaOptions);
    len--;
  }
  if(len!=1) return luaL_error(L, "ResponseFunction.EffectiveHamiltonian called with %d arguments, one argument expected and possibly a second specifying options",len);
  Gin = luaL_checkResponsefunction(L,1);
  Heff = (ResponsefunctionType *) lua_newuserdata(L, sizeof(ResponsefunctionType));
  A0 = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
  B0 = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
  if(ResponseFunctionEffectiveHamiltonian(Gin,Heff,A0,B0,Options.T,Options.SingularValue,Options.MinimalPoleWeigth,Options.MinimalPoleDistance))
  {
    return luaL_error(L, "Error in ResponseFunction.EffectiveHamiltonian()\n");
  }
  // we give the controle of the allocation of the memory of Heff to Lua
  // turn the type to a lower case letter
  Heff->type = Heff->type | 32;
  luaL_getmetatable(L, "ResponseFunction_Type");
  lua_setmetatable(L, -4);
  A0->OnLuaStack=true;
  luaL_getmetatable(L, "Matrix_Type");
  lua_setmetatable(L, -3);
  B0->OnLuaStack=true;
  luaL_getmetatable(L, "Matrix_Type");
  lua_setmetatable(L, -2);
  FreeResponsefunction(Gin);
  return 3;
}

static int LuaResponseFunctionFromEffectiveHamiltonian(lua_State *L)
{
 ResponsefunctionType *Gout,*Heff;
 double ValueNotSet = -DBL_MAX;
 struct
 {
   double SingularValue;
   double MinimalPoleWeigth;
   double MinimalPoleDistance;
   double T;
 } Options;
 Options.SingularValue       = ValueNotSet;
 Options.MinimalPoleWeigth   = ValueNotSet;
 Options.MinimalPoleDistance = ValueNotSet;
 Options.T                   = ValueNotSet;
 /*  the a -> l cell merges on MinimalPoleDistance and the n -> l one deflates on
  *  SingularValue; T and MinimalPoleWeigth are read by nothing on this route. */
 LuaOptionsType LuaOptions[]={
   {"SingularValue","double",&(Options.SingularValue)},
//   {"MinimalPoleWeigth","double",&(Options.MinimalPoleWeigth)},
   {"MinimalPoleDistance","double",&(Options.MinimalPoleDistance)},
//   {"Temperature","double",&(Options.T)},
//   {"EffectivekBT","double",&(Options.T)},
   {NULL, NULL, NULL}};
 uint32_t len = (uint32_t) lua_gettop(L);
 if(LuaTestOptionsList(L, -1))
 {
   LuaCheckOptions(L, -1, LuaOptions);
   len--;
 }
 if(len!=1) return luaL_error(L, "ResponseFunction.FromEffectiveHamiltonian called with %d arguments, one argument expected and possibly a second specifying options",len);
 Heff = luaL_checkResponsefunction(L,1);
 Gout = (ResponsefunctionType *) lua_newuserdata(L, sizeof(ResponsefunctionType));
 if(ResponseFunctionFromEffectiveHamiltonian(Heff, Gout, Options.T, Options.SingularValue, Options.MinimalPoleWeigth, Options.MinimalPoleDistance))
 {
   return luaL_error(L, "Error in ResponseFunction.FromEffectiveHamiltonian()\n");
 }
 // we give the controle of the allocation of the memory of Gout to Lua
 // turn the type to a lower case letter
 Gout->type = Gout->type | 32;
 luaL_setmetatable(L, "ResponseFunction_Type");
 FreeResponsefunction(Heff);
 return 1;
}

/*  Sigma, S = ResponseFunction.Downfold(G, R, {options})
 *
 *  The sub block of the inverse of G, inverted again,
 *
 *    Sigma(w) = ( R G^-1(w) R^+ )^-1 ,
 *
 *  with R a matrix of Nsub rows and as many columns as G has a block size. R = ( 0 | 1 )
 *  takes the lower right block of a two by two block structure, any other R of full row
 *  rank takes the subspace its rows span; the rows do not have to be orthonormal.
 *
 *  For G the augmented propagator of the operators a_m and q_m = [a_m,H1] and R selecting
 *  the q half, this is the symmetric improved estimator Sigma(w) - Sigma^HF of the self
 *  energy. Sigma^HF is not part of it and is read off the spectral weight of G itself,
 *  B0^* B0^T, which ResponseFunction.ChangeType(G,"Dense") returns.
 *
 *  S comes back as the spectral weight of the result, sum_i W_i, and is the sum rule the
 *  answer can be checked against.
 *
 *  A tri-diagonal G keeps its shape: its effective Hamiltonian is a chain, the projection
 *  is then a block Lanczos and nothing is diagonalised, and Sigma comes back a chain too.
 *  Every other type reaches the projection as a list of poles and Sigma comes back an
 *  Anderson matrix. ResponseFunction.ChangeType(Sigma,"ListOfPoles") turns either into the
 *  pole expansion Sigma(w) - Sigma^HF = sum_i W_i / ( w - a_i ), and pays a dense
 *  diagonalisation for it.
 *
 *  SingularValue deflates the directions of the spectral weight of G that carry no weight
 *  and MinimalPoleDistance merges poles that coincide, both as EffectiveHamiltonian takes
 *  them. NTriMax is a ceiling on the length of the chain the tri-diagonal route produces;
 *  leaving it out asks for no ceiling, which keeps the projection exact, and every other
 *  type refuses it, there being no chain to cut. See ResponseFunctionDownfold for what the
 *  routine does with them.
 */
static int LuaResponseFunctionDownfold(lua_State *L)
{
  ResponsefunctionType *G,*Gout;
  CompactMatrixType *R,*S;
  double ValueNotSet = -DBL_MAX;
  struct
  {
    double SingularValue;
    double MinimalPoleWeigth;
    double MinimalPoleDistance;
    double T;
    uint32_t NTriMax;
  } Options;
  Options.SingularValue       = ValueNotSet;
  Options.MinimalPoleWeigth   = ValueNotSet;
  Options.MinimalPoleDistance = ValueNotSet;
  Options.T                   = ValueNotSet;
  Options.NTriMax             = NTriMaxNoLimit;
  LuaOptionsType LuaOptions[]={
    {"SingularValue","double",&(Options.SingularValue)},
    {"MinimalPoleDistance","double",&(Options.MinimalPoleDistance)},
    {"NTriMax","uint32_t",&(Options.NTriMax)},
    {NULL, NULL, NULL}};
  uint32_t len = (uint32_t) lua_gettop(L);
  if(LuaTestOptionsList(L, -1))
  {
    LuaCheckOptions(L, -1, LuaOptions);
    len--;
  }
  if(len!=2) return luaL_error(L, "ResponseFunction.Downfold called with %d arguments, two arguments expected, a response function and the matrix selecting the subspace, and possibly a third specifying options",len);
  G = luaL_checkResponsefunction(L,1);
  R = luaL_checkCompactMatrix(L,2);
  if(G->Block==false)
  {
    return luaL_error(L, "Error in ResponseFunction.Downfold(). A single valued response function has no subspace to downfold onto\n");
  }
  if(ResponseFunctionGetBlockSize(G)!=R->Nj)
  {
    return luaL_error(L, "Error in ResponseFunction.Downfold(). The matrix selecting the subspace has %d columns and the response function has a block size of %d\n", R->Nj, ResponseFunctionGetBlockSize(G));
  }
  if(R->Ni>R->Nj)
  {
    return luaL_error(L, "Error in ResponseFunction.Downfold(). The matrix selecting the subspace has %d rows and %d columns, a subspace is not larger than the space it sits in\n", R->Ni, R->Nj);
  }
  Gout = (ResponsefunctionType *) lua_newuserdata(L, sizeof(ResponsefunctionType));
  S    = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
  if((Options.NTriMax!=NTriMaxNoLimit)&&((G->type|32)!='t'))
  {
    return luaL_error(L, "Error in ResponseFunction.Downfold(). NTriMax is a ceiling on the length of the chain the projection produces, and only a tri-diagonal response function is downfolded into a chain. A response function of type %c is downfolded one pole at a time and the result is a star\n", G->type|32);
  }
  if(ResponseFunctionDownfold(G,R,Gout,S,Options.T,Options.SingularValue,Options.MinimalPoleWeigth,Options.MinimalPoleDistance,Options.NTriMax))
  {
    return luaL_error(L, "Error in ResponseFunction.Downfold()\n");
  }
  // we give the controle of the allocation of the memory of Gout and S to Lua
  // turn the type to a lower case letter
  Gout->type = Gout->type | 32;
  luaL_getmetatable(L, "ResponseFunction_Type");
  lua_setmetatable(L, -3);
  S->OnLuaStack=true;
  luaL_getmetatable(L, "Matrix_Type");
  lua_setmetatable(L, -2);
  FreeResponsefunction(G);
  FreeCompactMatrix(R);
  return 2;
}

/*  G(w) -> G(w)[i1,i2,...in]
 *
 */
static int LuaResponseFunctionElement(lua_State *L)
{
  
  return 0;
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
  if(Ma->Block != Mb->Block) return luaL_error(L, "Error in ResponseFunctionSub: Block type of both functions needs to be the same\n");
  if(ResponseFunctionSubtract(Ma, Mb, &res)) return luaL_error(L, "Error in ResponseFunctionSub\n");
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
  double x;
  ResponsefunctionType *M,xM;
  if(lua_type(L,1)==LUA_TNUMBER)
  {
    x = luaL_checknumber(L, 1);
    M = luaL_checkResponsefunction(L,2);
  }
  else if(lua_type(L,2)==LUA_TNUMBER)
  {
    x = luaL_checknumber(L, 2);
    M = luaL_checkResponsefunction(L,1);
  }
  else
  {
    return luaL_error(L, "Error: You can only multiply a response function with a number\n");
  }
  if(x<0)  return luaL_error(L, "Error: response function multiplied by negative number\n");
  // if response function on lua stack copy otherwise set *xM to M
  xM = *M;
  if(M->type>94)
  { // lower case, data on lua stack, need to copy the data
    ResponsefunctionCopy(M, &xM);
  }
  if(ResponsefunctionMul(&xM,x))
  {
    FreeResponsefunction(M);
    return luaL_error(L, "Error scaling a response function by a number\n");
  }
  LuaPush_ResponseFunction(L, &xM);
  FreeResponsefunction(M);
  return 1;
}

/*
 *
 */
static int LuaResponseFunctionDiv(lua_State * L)
{
  double x;
  ResponsefunctionType *M,xM;
  if(lua_type(L,2)==LUA_TNUMBER)
  {
    x = luaL_checknumber(L, 2);
    M = luaL_checkResponsefunction(L,1);
    if(x==0) return luaL_error(L, "Error: response function divided by zero\n");
    if(x<0)  return luaL_error(L, "Error: response function divided by negative number\n");
    x=1/x;
  }
  else
  {
    return luaL_error(L, "Error: You can only divide a response function with a number\n");
  }
  // if response function on lua stack copy otherwise set *xM to M
  xM = *M;
  if(M->type>94)
  { // lower case, data on lua stack, need to copy the data
    ResponsefunctionCopy(M, &xM);
  }
  if(ResponsefunctionMul(&xM,x))
  {
    FreeResponsefunction(M);
    return luaL_error(L, "Error scaling a response function by a number\n");
  }
  LuaPush_ResponseFunction(L, &xM);
  FreeResponsefunction(M);
  return 1;
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
  ResponsefunctionType *res;
  res = luaL_checkResponsefunction(L,1);
    
  lua_pushcfunction(L, LuaResponseFunctionToTable);
  lua_pushvalue(L,1);
  lua_pcall(L, 1, 1, 0);
    
  if(lua_istable(L, -1)==0)
  {
    luaL_error(L, "Error in LuaResponseFunctionToString: Failed to convert ResponseFunction to a table");
  }

  lua_setglobal(L, "__temptable");
  if (luaL_dostring(L, "local result = tostring(__temptable); __temptable = nil; return result;") != LUA_OK)
  {
    luaL_error(L, "Error in LuaResponseFunctionToString: Successfuly converted userdata to table, but failed to convert table to string");
  }
  FreeResponsefunctionIncludingOriginalPointer(res);
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
static int ResponseFunctionToLuaBitMask(lua_State *L, ResponsefunctionType *M, bool AddSpin, uint32_t FermionShift, bool Deflate)
{
  int i,j,k;
  uint32_t NF,NB=0;
  uint32_t NCore, NImpurity, NBath, NBathVal, NBathCon, NBathExt, NImpurityExt;
  CompactMatrixType *T;
  char *Core,*Impurity,*ImpurityUp,*ImpurityDn,*Bath,*BathUp,*BathDn, *ImpurityExt, *ImpurityExtUp, *ImpurityExtDn;
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
          if(((BlockAndersonMatrixType*)(M->data))->N==0)
          {
            NImpurity=0;
            NBath=0;
          }
          else
          {
            NImpurity = 2*((BlockAndersonMatrixType*)(M->data))->BlockSize;
            NBath     = NImpurity * (((BlockAndersonMatrixType*)(M->data))->N-1);
          }
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
          if(((BlockAndersonMatrixType*)(M->data))->N==0)
          {
            NImpurity=0;
            NBath=0;
          }
          else
          {
            NImpurity = ((BlockAndersonMatrixType*)(M->data))->BlockSize;
            NBath     = NImpurity * (((BlockAndersonMatrixType*)(M->data))->N-1);
          }
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
          if(((AndersonMatrixType*)(M->data))->N==0)
          {
            NImpurity=0;
            NBath = 0;
          }
          else
          {
            NImpurity = 2;
            NBath     = NImpurity * (((AndersonMatrixType*)(M->data))->N-1);
          }
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
          if(((AndersonMatrixType*)(M->data))->N==0)
          {
            NImpurity=0;
            NBath = 0;
          }
          else
          {
            NImpurity = 1;
            NBath     = NImpurity * (((AndersonMatrixType*)(M->data))->N-1);
          }
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
          if(((BlockTriDiagonalMatrixType*)(M->data))->N==0)
          {
            NImpurity=0;
          }
          else
          {
            NImpurity = 2*((BlockTriDiagonalMatrixType*)(M->data))->BlockSize;
          }
          NCore     = FermionShift;
          NBath     = 0;
          if(Deflate==true)
          {
            for(k=2;k<((BlockTriDiagonalMatrixType*)(M->data))->N+1;k++)
            {
              NBath+=2*((BlockTriDiagonalMatrixType*)(M->data))->BlockSizeList[k];
            }
          }
          else
          {
            NBath = 2 * (((BlockTriDiagonalMatrixType*)(M->data))->N-1) * ((BlockTriDiagonalMatrixType*)(M->data))->BlockSize;
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
          if(((BlockTriDiagonalMatrixType*)(M->data))->N==0)
          {
            NImpurity=0;
          }
          else
          {
            NImpurity = ((BlockTriDiagonalMatrixType*)(M->data))->BlockSize;
          }
          NCore     = FermionShift;
          NBath     = 0;
          if(Deflate==true)
          {
            for(k=2;k<((BlockTriDiagonalMatrixType*)(M->data))->N+1;k++)
            {
              NBath+=((BlockTriDiagonalMatrixType*)(M->data))->BlockSizeList[k];
            }
          }
          else
          {
            NBath = (((BlockTriDiagonalMatrixType*)(M->data))->N-1) * ((BlockTriDiagonalMatrixType*)(M->data))->BlockSize;
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
          if(((TriDiagonalMatrixType*)(M->data))->N==0)
          {
            NImpurity = 0;
            NBath = 0;
          }
          else
          {
            NImpurity = 2;
            NBath     = NImpurity * (((TriDiagonalMatrixType*)(M->data))->N-1);
          }
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
          if(((TriDiagonalMatrixType*)(M->data))->N==0)
          {
            NImpurity = 0;
            NBath = 0;
          }
          else
          {
            NImpurity = 1;
            NBath     = NImpurity * (((TriDiagonalMatrixType*)(M->data))->N-1);
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
          NCore     = FermionShift;
          if(((BlockNaturalImpurityOrbitalMatrixType*)(M->data))->val.N + ((BlockNaturalImpurityOrbitalMatrixType*)(M->data))->con.N == 0)
          {
            NImpurity=0;
            NBathExt=0;
            NBathVal=0;
            NBathCon=0;
            NImpurityExt=0;
          }
          else
          {
            NImpurity = 2*((BlockNaturalImpurityOrbitalMatrixType*)(M->data))->BlockSize;
            NImpurityExt = 2*NImpurity;
            NBathExt  = NImpurity;
            if(Deflate==true)
            {
              NBathVal  = 0;
              for(k=2;k<=((BlockNaturalImpurityOrbitalMatrixType*)(M->data))->val.N;k++)
              {
                NBathVal+=2*((BlockNaturalImpurityOrbitalMatrixType*)(M->data))->val.BlockSizeList[k];
              }
              NBathCon  = 0;
              for(k=2;k<=((BlockNaturalImpurityOrbitalMatrixType*)(M->data))->con.N;k++)
              {
                NBathCon+=2*((BlockNaturalImpurityOrbitalMatrixType*)(M->data))->con.BlockSizeList[k];
              }
            }
            else
            {
              if(((BlockNaturalImpurityOrbitalMatrixType*)(M->data))->val.N==0)
              {
                NBathVal=0;
              }
              else
              {
                NBathVal = (((BlockNaturalImpurityOrbitalMatrixType*)(M->data))->val.N-1) * 2 * ((BlockNaturalImpurityOrbitalMatrixType*)(M->data))->BlockSize;
              }
              if(((BlockNaturalImpurityOrbitalMatrixType*)(M->data))->con.N==0)
              {
                NBathCon=0;
              }
              else
              {
                NBathCon = (((BlockNaturalImpurityOrbitalMatrixType*)(M->data))->con.N-1) * 2 * ((BlockNaturalImpurityOrbitalMatrixType*)(M->data))->BlockSize;
              }
            }
          }
          NBath     = NBathVal + NBathCon + NBathExt;
          NF = NCore + NImpurity + NBath;
          Core       = (char*)malloc((NF+1)*sizeof(char));
          Impurity   = (char*)malloc((NF+1)*sizeof(char));
          Bath       = (char*)malloc((NF+1)*sizeof(char));
          BathVal    = (char*)malloc((NF+1)*sizeof(char));
          BathCon    = (char*)malloc((NF+1)*sizeof(char));
          BathExt    = (char*)malloc((NF+1)*sizeof(char));
          ImpurityExt= (char*)malloc((NF+1)*sizeof(char));
          ImpurityUp = (char*)malloc((NF+1)*sizeof(char));
          BathUp     = (char*)malloc((NF+1)*sizeof(char));
          BathValUp  = (char*)malloc((NF+1)*sizeof(char));
          BathConUp  = (char*)malloc((NF+1)*sizeof(char));
          BathExtUp  = (char*)malloc((NF+1)*sizeof(char));
          ImpurityExtUp = (char*)malloc((NF+1)*sizeof(char));
          ImpurityDn = (char*)malloc((NF+1)*sizeof(char));
          BathDn     = (char*)malloc((NF+1)*sizeof(char));
          BathValDn  = (char*)malloc((NF+1)*sizeof(char));
          BathConDn  = (char*)malloc((NF+1)*sizeof(char));
          BathExtDn  = (char*)malloc((NF+1)*sizeof(char));
          ImpurityExtDn = (char*)malloc((NF+1)*sizeof(char));
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
            ImpurityExt[k]   ='0';
            ImpurityExtUp[k]   ='0';
            ImpurityExtDn[k]   ='0';

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
          ImpurityExt[NF]   ='\0';
          ImpurityExtUp[NF]   ='\0';
          ImpurityExtDn[NF]   ='\0';
          for(k=0;k<NCore;k++)
          {
            Core[k]='1';
          }
          for(k=NCore;k<NCore+NImpurity;k++)
          {
            Impurity[k]='1'; ImpurityDn[k]='1'; ImpurityUp[k]='0'; ImpurityExt[k]='1'; ImpurityExtDn[k]='1'; ImpurityExtUp[k]='0';
            k++;
            Impurity[k]='1'; ImpurityDn[k]='0'; ImpurityUp[k]='1'; ImpurityExt[k]='1'; ImpurityExtDn[k]='0'; ImpurityExtUp[k]='1';
          }
          for(k=NCore+NImpurity;k<NCore+NImpurity+NBathExt;k++)
          {
            BathExt[k]='1'; BathExtDn[k]='1'; BathExtUp[k]='0'; Bath[k]='1'; BathDn[k]='1'; ImpurityExt[k]='1'; ImpurityExtDn[k]='1'; ImpurityExtUp[k]='0';
            k++;
            BathExt[k]='1'; BathExtDn[k]='0'; BathExtUp[k]='1'; Bath[k]='1'; BathUp[k]='1'; ImpurityExt[k]='1'; ImpurityExtDn[k]='0'; ImpurityExtUp[k]='1';
          }
          for(k=NCore+NImpurity+NBathExt;k<NCore+NImpurity+NBathExt+NBathVal;k++)
          {
            BathVal[k]='1'; BathValDn[k]='1'; BathValUp[k]='0'; Bath[k]='1'; BathDn[k]='1';
            k++;
            BathVal[k]='1'; BathValDn[k]='0'; BathValUp[k]='1'; Bath[k]='1'; BathUp[k]='1';
          }
          for(k=NCore+NImpurity+NBathExt+NBathVal;k<NCore+NImpurity+NBathExt+NBathVal+NBathCon;k++)
          {
            BathCon[k]='1'; BathConDn[k]='1'; BathConUp[k]='0'; Bath[k]='1'; BathDn[k]='1';
            k++;
            BathCon[k]='1'; BathConDn[k]='0'; BathConUp[k]='1'; Bath[k]='1'; BathUp[k]='1';
          }
          lua_pushstring(L, ImpurityUp);
          lua_setfield(L,-2,"ImpurityUp");
          free(ImpurityUp);
          lua_pushstring(L, ImpurityDn);
          lua_setfield(L,-2,"ImpurityDn");
          free(ImpurityDn);
          lua_pushstring(L, ImpurityExtUp);
          lua_setfield(L,-2,"ImpurityExtUp");
          free(ImpurityExtUp);
          lua_pushstring(L, ImpurityExtDn);
          lua_setfield(L,-2,"ImpurityExtDn");
          free(ImpurityExtDn);
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
          if(((BlockNaturalImpurityOrbitalMatrixType*)(M->data))->val.N + ((BlockNaturalImpurityOrbitalMatrixType*)(M->data))->con.N == 0)
          {
            NImpurity=0;
            NBathExt=0;
            NBathVal=0;
            NBathCon=0;
            NImpurityExt=0;
          }
          else
          {
            NImpurity = ((BlockNaturalImpurityOrbitalMatrixType*)(M->data))->BlockSize;
            NImpurityExt=2*NImpurity;
            NBathExt  = NImpurity;
            if(Deflate==true)
            {
              NBathVal  = 0;
              for(k=2;k<=((BlockNaturalImpurityOrbitalMatrixType*)(M->data))->val.N;k++)
              {
                NBathVal+=((BlockNaturalImpurityOrbitalMatrixType*)(M->data))->val.BlockSizeList[k];
              }
              NBathCon  = 0;
              for(k=2;k<=((BlockNaturalImpurityOrbitalMatrixType*)(M->data))->con.N;k++)
              {
                NBathCon+=((BlockNaturalImpurityOrbitalMatrixType*)(M->data))->con.BlockSizeList[k];
              }
            }
            else
            {
              if(((BlockNaturalImpurityOrbitalMatrixType*)(M->data))->val.N==0)
              {
                NBathVal=0;
              }
              else
              {
                NBathVal = (((BlockNaturalImpurityOrbitalMatrixType*)(M->data))->val.N-1) * ((BlockNaturalImpurityOrbitalMatrixType*)(M->data))->BlockSize;
              }
              if(((BlockNaturalImpurityOrbitalMatrixType*)(M->data))->con.N==0)
              {
                NBathCon=0;
              }
              else
              {
                NBathCon = (((BlockNaturalImpurityOrbitalMatrixType*)(M->data))->con.N-1) * ((BlockNaturalImpurityOrbitalMatrixType*)(M->data))->BlockSize;
              }
            }
          }
          NBath     = NBathVal + NBathCon + NBathExt;
          NF = NCore + NImpurity + NBath;
          Core       = (char*)malloc((NF+1)*sizeof(char));
          Impurity   = (char*)malloc((NF+1)*sizeof(char));
          Bath       = (char*)malloc((NF+1)*sizeof(char));
          BathVal    = (char*)malloc((NF+1)*sizeof(char));
          BathCon    = (char*)malloc((NF+1)*sizeof(char));
          BathExt    = (char*)malloc((NF+1)*sizeof(char));
          ImpurityExt= (char*)malloc((NF+1)*sizeof(char));
          for(k=0;k<NF;k++)
          {
            Core[k]       ='0';
            Impurity[k]   ='0';
            Bath[k]       ='0';
            BathVal[k]    ='0';
            BathCon[k]    ='0';
            BathExt[k]    ='0';
            ImpurityExt[k]='0';
          }
          Core[NF]       ='\0';
          Impurity[NF]   ='\0';
          Bath[NF]       ='\0';
          BathVal[NF]    ='\0';
          BathCon[NF]    ='\0';
          BathExt[NF]    ='\0';
          ImpurityExt[NF]='\0';
          for(k=0;k<NCore;k++)
          {
            Core[k]='1';
          }
          for(k=NCore;k<NCore+NImpurity;k++)
          {
            Impurity[k]='1';
            ImpurityExt[k]='1';
          }
          for(k=NCore+NImpurity;k<NCore+NImpurity+NBathExt;k++)
          {
            BathExt[k]='1';
            Bath[k]='1';
            ImpurityExt[k]='1';
          }
          for(k=NCore+NImpurity+NBathExt;k<NCore+NImpurity+NBathExt+NBathVal;k++)
          {
            BathVal[k]='1';
            Bath[k]='1';
          }
          for(k=NCore+NImpurity+NBathExt+NBathVal;k<NCore+NImpurity+NBathExt+NBathVal+NBathCon;k++)
          {
            BathCon[k]='1';
            Bath[k]='1';
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
        lua_pushinteger(L, NImpurityExt);
        lua_setfield(L,-2,"NImpurityExt");
        lua_pushstring(L, ImpurityExt);
        lua_setfield(L,-2,"ImpurityExt");
        free(ImpurityExt);
        lua_pushinteger(L, NF);
        lua_setfield(L,-2,"NF");
        lua_pushinteger(L, NB);
        lua_setfield(L,-2,"NB");
        // push transformation matrix T
        T = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
        luaL_setmetatable(L, "Matrix_Type");
        T->OnLuaStack=true;
        lua_setfield(L,-2,"T");
        /*  the normalised T: the operator is built from TMT, which is the Hamiltonian in
         *  the orthonormal basis this transformation leads to */
        GetNormalizedTFromValandConInBlockNaturalImpurityOrbitalMatrixType((BlockNaturalImpurityOrbitalMatrixType*)(M->data), T);
        if(AddSpin==true)
        {
          CompactMatrixType Tspin;
          Tspin = MatrixEnlarge(*T, 2);
          T->OnLuaStack=false;
          FreeCompactMatrix(T);
          *T=Tspin;
          T->OnLuaStack=true;
        }
      }
      else
      {
        if(AddSpin==true)
        {
          NCore     = FermionShift;
          if(((BADoubleTriDiagonalMatrixType*)(M->data))->val.N + ((BADoubleTriDiagonalMatrixType*)(M->data))->con.N == 0)
          {
            NImpurity=0;
            NBathExt=0;
            NBathVal=0;
            NBathCon=0;
          }
          else
          {
            NImpurity = 2;
            NBathExt  = NImpurity;
            if(((BADoubleTriDiagonalMatrixType*)(M->data))->val.N==0)
            {
              NBathVal=0;
            }
            else
            {
              NBathVal = 2*((BADoubleTriDiagonalMatrixType*)(M->data))->val.N-2;
            }
            if(((BADoubleTriDiagonalMatrixType*)(M->data))->con.N==0)
            {
              NBathCon=0;
            }
            else
            {
              NBathCon = 2*((BADoubleTriDiagonalMatrixType*)(M->data))->con.N-2;
            }
          }
          NBath     = NBathVal + NBathCon + NBathExt;
          NF = NCore + NImpurity + NBath;
          Core          = (char*)malloc((NF+1)*sizeof(char));
          Impurity      = (char*)malloc((NF+1)*sizeof(char));
          ImpurityExt   = (char*)malloc((NF+1)*sizeof(char));
          Bath          = (char*)malloc((NF+1)*sizeof(char));
          BathVal       = (char*)malloc((NF+1)*sizeof(char));
          BathCon       = (char*)malloc((NF+1)*sizeof(char));
          BathExt       = (char*)malloc((NF+1)*sizeof(char));
          ImpurityUp    = (char*)malloc((NF+1)*sizeof(char));
          ImpurityExtUp = (char*)malloc((NF+1)*sizeof(char));
          BathUp        = (char*)malloc((NF+1)*sizeof(char));
          BathValUp     = (char*)malloc((NF+1)*sizeof(char));
          BathConUp     = (char*)malloc((NF+1)*sizeof(char));
          BathExtUp     = (char*)malloc((NF+1)*sizeof(char));
          ImpurityDn    = (char*)malloc((NF+1)*sizeof(char));
          ImpurityExtDn = (char*)malloc((NF+1)*sizeof(char));
          BathDn        = (char*)malloc((NF+1)*sizeof(char));
          BathValDn     = (char*)malloc((NF+1)*sizeof(char));
          BathConDn     = (char*)malloc((NF+1)*sizeof(char));
          BathExtDn     = (char*)malloc((NF+1)*sizeof(char));
          for(k=0;k<NF;k++)
          {
            Core[k]          ='0';
            Impurity[k]      ='0';
            ImpurityExt[k]   ='0';
            Bath[k]          ='0';
            BathVal[k]       ='0';
            BathCon[k]       ='0';
            BathExt[k]       ='0';
            ImpurityUp[k]    ='0';
            ImpurityExtUp[k] ='0';
            BathUp[k]        ='0';
            BathValUp[k]     ='0';
            BathConUp[k]     ='0';
            BathExtUp[k]     ='0';
            ImpurityDn[k]    ='0';
            ImpurityExtDn[k] ='0';
            BathDn[k]        ='0';
            BathValDn[k]     ='0';
            BathConDn[k]     ='0';
            BathExtDn[k]     ='0';
          }
          Core[NF]          ='\0';
          Impurity[NF]      ='\0';
          ImpurityExt[NF]   ='\0';
          Bath[NF]          ='\0';
          BathVal[NF]       ='\0';
          BathCon[NF]       ='\0';
          BathExt[NF]       ='\0';
          ImpurityUp[NF]    ='\0';
          ImpurityExtUp[NF] ='\0';
          BathUp[NF]        ='\0';
          BathValUp[NF]     ='\0';
          BathConUp[NF]     ='\0';
          BathExtUp[NF]     ='\0';
          ImpurityDn[NF]    ='\0';
          ImpurityExtDn[NF] ='\0';
          BathDn[NF]        ='\0';
          BathValDn[NF]     ='\0';
          BathConDn[NF]     ='\0';
          BathExtDn[NF]     ='\0';
          for(k=0;k<NCore;k++)
          {
            Core[k]='1';
          }
          for(k=NCore;k<NCore+NImpurity;k++)
          {
            Impurity[k]='1'; ImpurityDn[k]='1'; ImpurityUp[k]='0'; ImpurityExt[k]='1'; ImpurityExtDn[k]='1'; ImpurityExtUp[k]='0';
            k++;
            Impurity[k]='1'; ImpurityDn[k]='0'; ImpurityUp[k]='1'; ImpurityExt[k]='1'; ImpurityExtDn[k]='0'; ImpurityExtUp[k]='1';
          }
          for(k=NCore+NImpurity;k<NCore+NImpurity+NBathExt;k++)
          {
            BathExt[k]='1'; BathExtDn[k]='1'; BathExtUp[k]='0'; Bath[k]='1'; BathDn[k]='1'; BathUp[k]='0'; ImpurityExt[k]='1'; ImpurityExtDn[k]='1'; ImpurityExtUp[k]='0';
            k++;
            BathExt[k]='1'; BathExtDn[k]='0'; BathExtUp[k]='1'; Bath[k]='1'; BathDn[k]='0'; BathUp[k]='1'; ImpurityExt[k]='1'; ImpurityExtDn[k]='0'; ImpurityExtUp[k]='1';
          }
          for(k=NCore+NImpurity+NBathExt;k<NCore+NImpurity+NBathExt+NBathVal;k++)
          {
            BathVal[k]='1'; BathValDn[k]='1'; BathValUp[k]='0'; Bath[k]='1'; BathDn[k]='1'; BathUp[k]='0';
            k++;
            BathVal[k]='1'; BathValDn[k]='0'; BathValUp[k]='1'; Bath[k]='1'; BathDn[k]='0'; BathUp[k]='1';
          }
          for(k=NCore+NImpurity+NBathExt+NBathVal;k<NCore+NImpurity+NBathExt+NBathVal+NBathCon;k++)
          {
            BathCon[k]='1'; BathConDn[k]='1'; BathConUp[k]='0'; Bath[k]='1'; BathDn[k]='1'; BathUp[k]='0';
            k++;
            BathCon[k]='1'; BathConDn[k]='0'; BathConUp[k]='1'; Bath[k]='1'; BathDn[k]='0'; BathUp[k]='1';
          }
          lua_pushstring(L, ImpurityUp);
          lua_setfield(L,-2,"ImpurityUp");
          free(ImpurityUp);
          lua_pushstring(L, ImpurityDn);
          lua_setfield(L,-2,"ImpurityDn");
          free(ImpurityDn);
          lua_pushstring(L, ImpurityExtUp);
          lua_setfield(L,-2,"ImpurityExtUp");
          free(ImpurityExtUp);
          lua_pushstring(L, ImpurityExtDn);
          lua_setfield(L,-2,"ImpurityExtDn");
          free(ImpurityExtDn);
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
          if(((BADoubleTriDiagonalMatrixType*)(M->data))->val.N + ((BADoubleTriDiagonalMatrixType*)(M->data))->con.N == 0)
          {
            NImpurity=0;
            NBathExt=0;
            NBathVal=0;
            NBathCon=0;
          }
          else
          {
            NImpurity = 1;
            NBathExt  = NImpurity;
            if(((BADoubleTriDiagonalMatrixType*)(M->data))->val.N==0)
            {
              NBathVal=0;
            }
            else
            {
              NBathVal = ((BADoubleTriDiagonalMatrixType*)(M->data))->val.N-1;
            }
            if(((BADoubleTriDiagonalMatrixType*)(M->data))->con.N==0)
            {
              NBathCon=0;
            }
            else
            {
              NBathCon = ((BADoubleTriDiagonalMatrixType*)(M->data))->con.N-1;
            }
          }
          NBath     = NBathVal + NBathCon + NBathExt;
          NF = NCore + NImpurity + NBath;
          Core       = (char*)malloc((NF+1)*sizeof(char));
          Impurity   = (char*)malloc((NF+1)*sizeof(char));
          ImpurityExt= (char*)malloc((NF+1)*sizeof(char));
          Bath       = (char*)malloc((NF+1)*sizeof(char));
          BathVal    = (char*)malloc((NF+1)*sizeof(char));
          BathCon    = (char*)malloc((NF+1)*sizeof(char));
          BathExt    = (char*)malloc((NF+1)*sizeof(char));
          for(k=0;k<NF;k++)
          {
            Core[k]       ='0';
            Impurity[k]   ='0';
            ImpurityExt[k]='0';
            Bath[k]       ='0';
            BathVal[k]    ='0';
            BathCon[k]    ='0';
            BathExt[k]    ='0';
          }
          Core[NF]       ='\0';
          Impurity[NF]   ='\0';
          ImpurityExt[NF]='\0';
          Bath[NF]       ='\0';
          BathVal[NF]    ='\0';
          BathCon[NF]    ='\0';
          BathExt[NF]    ='\0';
          for(k=0;k<NCore;k++)
          {
            Core[k]='1';
          }
          for(k=NCore;k<NCore+NImpurity;k++)
          {
            Impurity[k]='1';
            ImpurityExt[k]='1';
          }
          for(k=NCore+NImpurity;k<NCore+NImpurity+NBathExt;k++)
          {
            BathExt[k]='1';
            Bath[k]='1';
            ImpurityExt[k]='1';
          }
          for(k=NCore+NImpurity+NBathExt;k<NCore+NImpurity+NBathExt+NBathVal;k++)
          {
            BathVal[k]='1';
            Bath[k]='1';
          }
          for(k=NCore+NImpurity+NBathExt+NBathVal;k<NCore+NImpurity+NBathExt+NBathVal+NBathCon;k++)
          {
            BathCon[k]='1';
            Bath[k]='1';
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
        lua_pushstring(L, ImpurityExt);
        lua_setfield(L,-2,"ImpurityExt");
        free(ImpurityExt);
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
        // push transformation matrix T
        T = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
        luaL_setmetatable(L, "Matrix_Type");
        lua_setfield(L,-2,"T");
        T->QComplex=0;
        T->Ni=2;
        T->Nj=2;
        InitCompactMatrix(T);
        T->OnLuaStack=true;
        if( ((BADoubleTriDiagonalMatrixType*)(M->data))->val.N==0  )
        {
          T->VRe[0]=0;
          T->VRe[3]=0;
        }
        else
        {
          T->VRe[0]= ( ((BADoubleTriDiagonalMatrixType*)(M->data))->val.b[0] );
          T->VRe[3]=-( ((BADoubleTriDiagonalMatrixType*)(M->data))->val.b[0] );
        }
        if( ((BADoubleTriDiagonalMatrixType*)(M->data))->con.N == 0)
        {
          T->VRe[1]=0;
          T->VRe[2]=0;
        }
        else
        {
          T->VRe[1]= ( ((BADoubleTriDiagonalMatrixType*)(M->data))->con.b[0] );
          T->VRe[2]= ( ((BADoubleTriDiagonalMatrixType*)(M->data))->con.b[0] );
        }
        if(AddSpin==true)
        {
          CompactMatrixType Tspin;
          Tspin = MatrixEnlarge(*T, 2);
          T->OnLuaStack=false;
          FreeCompactMatrix(T);
          *T=Tspin;
          T->OnLuaStack=true;
        }
      }
      break;
    }
    case 'D':
    case 'd':
    {
      lua_pushstring(L, "DoubleTri");
      lua_setfield(L,-2,"type");
      if(M->Block==true)
      {
        if(AddSpin==true)
        {
          NCore     = FermionShift;
          /*  a[0] of a chain is the impurity energy and a[1] to a[N] are all bath, so a
           *  chain contributes N bath blocks and the deflation sums run from 1. The
           *  natural orbital form this was copied from counts one block less per chain,
           *  a[1] being the natural orbital there and counted with the impurity, and its
           *  sums start at 2.  A chain of length zero needs no special case: N * BlockSize
           *  is then zero, and con.N == 0 is a common case for this type. */
          NImpurity = 2*((BlockDoubleTriDiagonalMatrixType*)(M->data))->BlockSize;
          if(Deflate==true)
          {
            NBathVal  = 0;
            for(k=1;k<=((BlockDoubleTriDiagonalMatrixType*)(M->data))->val.N;k++)
            {
              NBathVal+=2*((BlockDoubleTriDiagonalMatrixType*)(M->data))->val.BlockSizeList[k];
            }
            NBathCon  = 0;
            for(k=1;k<=((BlockDoubleTriDiagonalMatrixType*)(M->data))->con.N;k++)
            {
              NBathCon+=2*((BlockDoubleTriDiagonalMatrixType*)(M->data))->con.BlockSizeList[k];
            }
          }
          else
          {
            NBathVal = ((BlockDoubleTriDiagonalMatrixType*)(M->data))->val.N * 2 * ((BlockDoubleTriDiagonalMatrixType*)(M->data))->BlockSize;
            NBathCon = ((BlockDoubleTriDiagonalMatrixType*)(M->data))->con.N * 2 * ((BlockDoubleTriDiagonalMatrixType*)(M->data))->BlockSize;
          }
          NBath     = NBathVal + NBathCon;
          NF = NCore + NImpurity + NBath;
          Core       = (char*)malloc((NF+1)*sizeof(char));
          Impurity   = (char*)malloc((NF+1)*sizeof(char));
          Bath       = (char*)malloc((NF+1)*sizeof(char));
          BathVal    = (char*)malloc((NF+1)*sizeof(char));
          BathCon    = (char*)malloc((NF+1)*sizeof(char));
          ImpurityUp = (char*)malloc((NF+1)*sizeof(char));
          BathUp     = (char*)malloc((NF+1)*sizeof(char));
          BathValUp  = (char*)malloc((NF+1)*sizeof(char));
          BathConUp  = (char*)malloc((NF+1)*sizeof(char));
          ImpurityDn = (char*)malloc((NF+1)*sizeof(char));
          BathDn     = (char*)malloc((NF+1)*sizeof(char));
          BathValDn  = (char*)malloc((NF+1)*sizeof(char));
          BathConDn  = (char*)malloc((NF+1)*sizeof(char));
          for(k=0;k<NF;k++)
          {
            Core[k]       ='0';
            Impurity[k]   ='0';
            Bath[k]       ='0';
            BathVal[k]    ='0';
            BathCon[k]    ='0';
            ImpurityUp[k] ='0';
            BathUp[k]     ='0';
            BathValUp[k]  ='0';
            BathConUp[k]  ='0';
            ImpurityDn[k] ='0';
            BathDn[k]     ='0';
            BathValDn[k]  ='0';
            BathConDn[k]  ='0';

          }
          Core[NF]       ='\0';
          Impurity[NF]   ='\0';
          Bath[NF]       ='\0';
          BathVal[NF]    ='\0';
          BathCon[NF]    ='\0';
          ImpurityUp[NF] ='\0';
          BathUp[NF]     ='\0';
          BathValUp[NF]  ='\0';
          BathConUp[NF]  ='\0';
          ImpurityDn[NF] ='\0';
          BathDn[NF]     ='\0';
          BathValDn[NF]  ='\0';
          BathConDn[NF]  ='\0';
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
          for(k=NCore+NImpurity;k<NCore+NImpurity+NBathVal;k++)
          {
            BathVal[k]='1'; BathValDn[k]='1'; BathValUp[k]='0'; Bath[k]='1'; BathDn[k]='1';
            k++;
            BathVal[k]='1'; BathValDn[k]='0'; BathValUp[k]='1'; Bath[k]='1'; BathUp[k]='1';
          }
          for(k=NCore+NImpurity+NBathVal;k<NCore+NImpurity+NBathVal+NBathCon;k++)
          {
            BathCon[k]='1'; BathConDn[k]='1'; BathConUp[k]='0'; Bath[k]='1'; BathDn[k]='1';
            k++;
            BathCon[k]='1'; BathConDn[k]='0'; BathConUp[k]='1'; Bath[k]='1'; BathUp[k]='1';
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
        }
        else // AddSpin == false
        {
          NCore     = FermionShift;
          /*  a[0] of a chain is the impurity energy and a[1] to a[N] are all bath, so a
           *  chain contributes N bath blocks and the deflation sums run from 1. The
           *  natural orbital form this was copied from counts one block less per chain,
           *  a[1] being the natural orbital there and counted with the impurity, and its
           *  sums start at 2.  A chain of length zero needs no special case: N * BlockSize
           *  is then zero, and con.N == 0 is a common case for this type. */
          NImpurity = ((BlockDoubleTriDiagonalMatrixType*)(M->data))->BlockSize;
          if(Deflate==true)
          {
            NBathVal  = 0;
            for(k=1;k<=((BlockDoubleTriDiagonalMatrixType*)(M->data))->val.N;k++)
            {
              NBathVal+=((BlockDoubleTriDiagonalMatrixType*)(M->data))->val.BlockSizeList[k];
            }
            NBathCon  = 0;
            for(k=1;k<=((BlockDoubleTriDiagonalMatrixType*)(M->data))->con.N;k++)
            {
              NBathCon+=((BlockDoubleTriDiagonalMatrixType*)(M->data))->con.BlockSizeList[k];
            }
          }
          else
          {
            NBathVal = ((BlockDoubleTriDiagonalMatrixType*)(M->data))->val.N * ((BlockDoubleTriDiagonalMatrixType*)(M->data))->BlockSize;
            NBathCon = ((BlockDoubleTriDiagonalMatrixType*)(M->data))->con.N * ((BlockDoubleTriDiagonalMatrixType*)(M->data))->BlockSize;
          }
          NBath     = NBathVal + NBathCon;
          NF = NCore + NImpurity + NBath;
          Core       = (char*)malloc((NF+1)*sizeof(char));
          Impurity   = (char*)malloc((NF+1)*sizeof(char));
          Bath       = (char*)malloc((NF+1)*sizeof(char));
          BathVal    = (char*)malloc((NF+1)*sizeof(char));
          BathCon    = (char*)malloc((NF+1)*sizeof(char));
          for(k=0;k<NF;k++)
          {
            Core[k]       ='0';
            Impurity[k]   ='0';
            Bath[k]       ='0';
            BathVal[k]    ='0';
            BathCon[k]    ='0';
          }
          Core[NF]       ='\0';
          Impurity[NF]   ='\0';
          Bath[NF]       ='\0';
          BathVal[NF]    ='\0';
          BathCon[NF]    ='\0';
          for(k=0;k<NCore;k++)
          {
            Core[k]='1';
          }
          for(k=NCore;k<NCore+NImpurity;k++)
          {
            Impurity[k]='1';
          }
          for(k=NCore+NImpurity;k<NCore+NImpurity+NBathVal;k++)
          {
            BathVal[k]='1';
            Bath[k]='1';
          }
          for(k=NCore+NImpurity+NBathVal;k<NCore+NImpurity+NBathVal+NBathCon;k++)
          {
            BathCon[k]='1';
            Bath[k]='1';
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
          /*  a[0] of a chain is the impurity energy and a[1] to a[N] are all bath, so a
           *  chain contributes N bath blocks and the deflation sums run from 1. The
           *  natural orbital form this was copied from counts one block less per chain,
           *  a[1] being the natural orbital there and counted with the impurity, and its
           *  sums start at 2.  A chain of length zero needs no special case: N * BlockSize
           *  is then zero, and con.N == 0 is a common case for this type. */
          NImpurity = 2;
          NBathVal = 2 * ((DoubleTriDiagonalMatrixType*)(M->data))->val.N;
          NBathCon = 2 * ((DoubleTriDiagonalMatrixType*)(M->data))->con.N;
          NBath     = NBathVal + NBathCon;
          NF = NCore + NImpurity + NBath;
          Core          = (char*)malloc((NF+1)*sizeof(char));
          Impurity      = (char*)malloc((NF+1)*sizeof(char));
          Bath          = (char*)malloc((NF+1)*sizeof(char));
          BathVal       = (char*)malloc((NF+1)*sizeof(char));
          BathCon       = (char*)malloc((NF+1)*sizeof(char));
          ImpurityUp    = (char*)malloc((NF+1)*sizeof(char));
          BathUp        = (char*)malloc((NF+1)*sizeof(char));
          BathValUp     = (char*)malloc((NF+1)*sizeof(char));
          BathConUp     = (char*)malloc((NF+1)*sizeof(char));
          ImpurityDn    = (char*)malloc((NF+1)*sizeof(char));
          BathDn        = (char*)malloc((NF+1)*sizeof(char));
          BathValDn     = (char*)malloc((NF+1)*sizeof(char));
          BathConDn     = (char*)malloc((NF+1)*sizeof(char));
          for(k=0;k<NF;k++)
          {
            Core[k]          ='0';
            Impurity[k]      ='0';
            Bath[k]          ='0';
            BathVal[k]       ='0';
            BathCon[k]       ='0';
            ImpurityUp[k]    ='0';
            BathUp[k]        ='0';
            BathValUp[k]     ='0';
            BathConUp[k]     ='0';
            ImpurityDn[k]    ='0';
            BathDn[k]        ='0';
            BathValDn[k]     ='0';
            BathConDn[k]     ='0';
          }
          Core[NF]          ='\0';
          Impurity[NF]      ='\0';
          Bath[NF]          ='\0';
          BathVal[NF]       ='\0';
          BathCon[NF]       ='\0';
          ImpurityUp[NF]    ='\0';
          BathUp[NF]        ='\0';
          BathValUp[NF]     ='\0';
          BathConUp[NF]     ='\0';
          ImpurityDn[NF]    ='\0';
          BathDn[NF]        ='\0';
          BathValDn[NF]     ='\0';
          BathConDn[NF]     ='\0';
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
          for(k=NCore+NImpurity;k<NCore+NImpurity+NBathVal;k++)
          {
            BathVal[k]='1'; BathValDn[k]='1'; BathValUp[k]='0'; Bath[k]='1'; BathDn[k]='1'; BathUp[k]='0';
            k++;
            BathVal[k]='1'; BathValDn[k]='0'; BathValUp[k]='1'; Bath[k]='1'; BathDn[k]='0'; BathUp[k]='1';
          }
          for(k=NCore+NImpurity+NBathVal;k<NCore+NImpurity+NBathVal+NBathCon;k++)
          {
            BathCon[k]='1'; BathConDn[k]='1'; BathConUp[k]='0'; Bath[k]='1'; BathDn[k]='1'; BathUp[k]='0';
            k++;
            BathCon[k]='1'; BathConDn[k]='0'; BathConUp[k]='1'; Bath[k]='1'; BathDn[k]='0'; BathUp[k]='1';
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
        }
        else // AddSpin == false
        {
          NCore     = FermionShift;
          /*  a[0] of a chain is the impurity energy and a[1] to a[N] are all bath, so a
           *  chain contributes N bath blocks and the deflation sums run from 1. The
           *  natural orbital form this was copied from counts one block less per chain,
           *  a[1] being the natural orbital there and counted with the impurity, and its
           *  sums start at 2.  A chain of length zero needs no special case: N * BlockSize
           *  is then zero, and con.N == 0 is a common case for this type. */
          NImpurity = 1;
          NBathVal = ((DoubleTriDiagonalMatrixType*)(M->data))->val.N;
          NBathCon = ((DoubleTriDiagonalMatrixType*)(M->data))->con.N;
          NBath     = NBathVal + NBathCon;
          NF = NCore + NImpurity + NBath;
          Core       = (char*)malloc((NF+1)*sizeof(char));
          Impurity   = (char*)malloc((NF+1)*sizeof(char));
          Bath       = (char*)malloc((NF+1)*sizeof(char));
          BathVal    = (char*)malloc((NF+1)*sizeof(char));
          BathCon    = (char*)malloc((NF+1)*sizeof(char));
          for(k=0;k<NF;k++)
          {
            Core[k]       ='0';
            Impurity[k]   ='0';
            Bath[k]       ='0';
            BathVal[k]    ='0';
            BathCon[k]    ='0';
          }
          Core[NF]       ='\0';
          Impurity[NF]   ='\0';
          Bath[NF]       ='\0';
          BathVal[NF]    ='\0';
          BathCon[NF]    ='\0';
          for(k=0;k<NCore;k++)
          {
            Core[k]='1';
          }
          for(k=NCore;k<NCore+NImpurity;k++)
          {
            Impurity[k]='1';
          }
          for(k=NCore+NImpurity;k<NCore+NImpurity+NBathVal;k++)
          {
            BathVal[k]='1';
            Bath[k]='1';
          }
          for(k=NCore+NImpurity+NBathVal;k<NCore+NImpurity+NBathVal+NBathCon;k++)
          {
            BathCon[k]='1';
            Bath[k]='1';
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
        lua_pushinteger(L, NF);
        lua_setfield(L,-2,"NF");
        lua_pushinteger(L, NB);
        lua_setfield(L,-2,"NB");
      }
      break;
    }
    default:
    {
      return luaL_error(L, "Error in Responsefunction.CreateBitMasks() unrecognized type in response function (%c)\n",M->type);
    }
  }
  return 0;
}
static int LuaResponseFunctionToBitMask(lua_State *L)
{
  ResponsefunctionType *M;
  bool AddSpin=false;
  bool Deflate=false;
  uint32_t FermionShift=0;
  LuaOptionsType LuaOptions[]={
    {"Deflate","bool",&Deflate},
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
  ResponseFunctionToLuaBitMask(L, M, AddSpin, FermionShift, Deflate);
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
  M = luaL_checkResponsefunction(L,1);
  TB = (TightBindingDef*) lua_newuserdata(L, sizeof(TightBindingDef));
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
        return luaL_error(L, "Error in Responsefunction.ToTightBinding() Not yet implemented for this case\n");
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
        return luaL_error(L, "Error in Responsefunction.ToTightBinding() Not yet implemented for this case\n");
      }
      break;
    }
    case 'N':
    case 'n':
    {
      if(M->Block==true)
      {
        return luaL_error(L, "Error in Responsefunction.ToTightBinding() Not yet implemented for this case\n");
        BlockNaturalImpurityOrbitalMatrixToTightBinding((BlockNaturalImpurityOrbitalMatrixType *)(M->data),TB);
      }
      else
      {
        //TODO: NatrualImpuritiOrbitalMatrixType to TB
        return luaL_error(L, "Error in Responsefunction.ToTightBinding() Not yet implemented for this case\n");
      }
      break;
    }
    case 'D':
    case 'd':
    {
      //TODO: (Block)DoubleTriDiagonalMatrixType to TB
      return luaL_error(L, "Error in Responsefunction.ToTightBinding() Not yet implemented for this case\n");
      break;
    }
    default:
    {
      return luaL_error(L, "Error in Responsefunction.ToTightbinding() unrecognized type for conversion (%c)\n",M->type);
      return 0;
    }
  }
  FreeResponsefunctionIncludingOriginalPointer(M);
  luaL_setmetatable(L, "TightBinding_Type");
  return 1;
}


/*   On the lua stack we have a response function and optional a set of options
 *   pushes an operator on the stack
 *   push information about the operator to the stack
 *   if Ncon and or Nval is set and G is in And or Nat format aditionally pushes a table onto the stack that contians Gv, and or Gc
 */
static int LuaResponseFunctionToOperator(lua_State * L)
{
  OperatorType *O;
  ResponsefunctionType *M;
  uint32_t MNv=0, MNc=0;
  struct
  {
    bool AddSpin;
    bool Deflation;
    bool TransformImpurity;
    uint32_t Nv;
    uint32_t Nc;
    uint32_t FermionShift;
  } Options;
  Options.Deflation=false;
  Options.AddSpin=false;
  Options.FermionShift=0;
  Options.TransformImpurity = false;
  Options.Nv = UINT32_MAX;
  Options.Nc = UINT32_MAX;
  LuaOptionsType LuaOptions[]={
    {"AddSpin","bool",&(Options.AddSpin)},
 // Option not implemented  {"TransformImpurity","bool",&(Options.TransformImpurity)},
    {"FermionShift","uint32_t",&(Options.FermionShift)},
    {"Nc","uint32_t",&(Options.Nc)},
    {"Nv","uint32_t",&(Options.Nv)},
    {"Ncon","uint32_t",&(Options.Nc)},
    {"Nval","uint32_t",&(Options.Nv)},
    {NULL, NULL, NULL}};
  uint32_t len = (uint32_t) lua_gettop(L);
  if(LuaTestOptionsList(L, -1))
  {
    LuaCheckOptions(L, -1, LuaOptions);
    len--;
  }
  if(len!=1) return luaL_error(L, "Error in ResponseFunction.ToOperator() It should be called with 1 argument, not with %d\n",len);
  M = luaL_checkResponsefunction(L,1);
  if((Options.TransformImpurity==true)&&(!((M->type=='N')||(M->type=='n')))) return luaL_error(L, "Option TransformImpurity only availalbe for Response Functions in natural impurity orbital form");
  /*  Nv and Nc truncate the two chains, so they mean something for both types that have
   *  two: the natural orbital form and the double tri-diagonal one */
  if((Options.Nc!=UINT32_MAX)&&(!((M->type=='N')||(M->type=='n')||(M->type=='D')||(M->type=='d')))) return luaL_error(L, "Option Nc only availalbe for Response Functions in natural impurity orbital or double tri-diagonal form");
  if((Options.Nv!=UINT32_MAX)&&(!((M->type=='N')||(M->type=='n')||(M->type=='D')||(M->type=='d')))) return luaL_error(L, "Option Nv only availalbe for Response Functions in natural impurity orbital or double tri-diagonal form");
  O = (OperatorType*) lua_newuserdata(L, sizeof(OperatorType));
  luaL_setmetatable(L, "Operator_Type");
  switch(M->type)
  {
    case 'L':
    case 'l':
    {
      // make sure to close Lua and trhough an error without an aditional free of unallocated space
      InitOperator(O);
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
        AndersonMatrixToImpurityHamiltonian(O, (AndersonMatrixType*)M->data, Options.AddSpin, Options.FermionShift, 0);
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
        TriDiagonalMatrixToImpurityHamiltonian(O, (TriDiagonalMatrixType*)M->data, Options.AddSpin, Options.FermionShift, 0);
      }
      break;
    }
    case 'N':
    case 'n':
    {
      // returns the operator, and the table of bitmasks, but aditionally, when Nv or Nc are set also
      // a table containing 2 elements:
      //    GvalBath
      //    GconBath
      //    These contain the tails of the vanlence and conduction chain
      //
      if((Options.Nv==0)||(Options.Nc==0))
      {
        //TODO: make error go out gracious without crashing
        return luaL_error(L, "Error in ResponseFunction.ToOperator()\nFor Natural impurity orbital geometry the 1-st valence and conduction bath side are part of the impurity, setting Nv=0 or Nc=0 therefor is not an option\n");
      }
      if(M->Block==true)
      {
        // copy and shorten M->data if Nv and Nc set
        MNv = ((BlockNaturalImpurityOrbitalMatrixType*)(M->data))->val.N;
        MNc = ((BlockNaturalImpurityOrbitalMatrixType*)(M->data))->con.N;
        if(MNv>Options.Nv) ((BlockNaturalImpurityOrbitalMatrixType*)(M->data))->val.N=Options.Nv;
        if(MNc>Options.Nc) ((BlockNaturalImpurityOrbitalMatrixType*)(M->data))->con.N=Options.Nc;
        if(Options.TransformImpurity == true)
        {
          // TODO -- currently works as Options.TransformImpurity can not be true
        }
        else // impurity is not transformed
        {
          BlockNaturalImpurityOrbitalMatrixToImpurityHamiltonian(O, (BlockNaturalImpurityOrbitalMatrixType*)(M->data), Options.AddSpin, Options.FermionShift, 0);
        }
      }
      else
      {
        // copy and shorten M->data if Nv and Nc set
        MNv = ((BADoubleTriDiagonalMatrixType*)M->data)->val.N;
        MNc = ((BADoubleTriDiagonalMatrixType*)M->data)->con.N;
        if(MNv>Options.Nv) ((BADoubleTriDiagonalMatrixType*)M->data)->val.N=Options.Nv;
        if(MNc>Options.Nc) ((BADoubleTriDiagonalMatrixType*)M->data)->con.N=Options.Nc;
        if(Options.TransformImpurity == true)
        {
          // TODO -- currently works as Options.TransformImpurity can not be true
        }
        else // impurity is not transformed
        {
          BADoubleTriDiagonalMatrixToImpurityHamiltonian(O, (BADoubleTriDiagonalMatrixType*)M->data, Options.AddSpin, Options.FermionShift, 0);
        }
      }
      if((MNv>Options.Nv)||(MNc>Options.Nc))
      {
        ResponseFunctionToLuaBitMask(L, M, Options.AddSpin, Options.FermionShift,false);
        // undo changes to M
        if(M->Block==true)
        {
          // expand M->data if Nv and Nc set
          if(MNv>Options.Nv) ((BlockNaturalImpurityOrbitalMatrixType*)(M->data))->val.N=MNv;
          if(MNc>Options.Nc) ((BlockNaturalImpurityOrbitalMatrixType*)(M->data))->con.N=MNc;
        }
        else
        {
          // expand M->data if Nv and Nc set
          if(MNv>Options.Nv) ((BADoubleTriDiagonalMatrixType*)M->data)->val.N=MNv;
          if(MNc>Options.Nc) ((BADoubleTriDiagonalMatrixType*)M->data)->con.N=MNc;
        }
        // push Gv, and, Gc to stack
        lua_createtable(L,0,4);
        // nested if then else statements to select cases for AddSpin, Block and Nv,Nc set
        if(Options.AddSpin==true)
        {
          if(M->Block==true) // add spin and block true
          {
            if(MNv>Options.Nv)
            {
              ResponsefunctionType *GvalBath;
              lua_pushstring(L, "GvalBath");
              GvalBath = (ResponsefunctionType *) lua_newuserdata(L, sizeof(ResponsefunctionType));
              luaL_setmetatable(L, "ResponseFunction_Type");
              GvalBath->type = 't';
              GvalBath->Block=true;
              GvalBath->NSuggestedGridPoints=0;
              GvalBath->SuggestedGrid=NULL;
              GvalBath->data = malloc(sizeof(BlockTriDiagonalMatrixType));
              BlockTriDiagonalMatrixGetTailAndAddSpin(&(((BlockNaturalImpurityOrbitalMatrixType*)(M->data))->val), (BlockTriDiagonalMatrixType*)(GvalBath->data), MNv-Options.Nv);
              BlockTriDiagonalMatrixSetA0to0((BlockTriDiagonalMatrixType*)(GvalBath->data));
              lua_settable(L, -3);
              lua_pushstring(L, "GvalBathIndex");
              lua_createtable(L,((BlockTriDiagonalMatrixType*)(GvalBath->data))->BlockSize,0);
              if(Options.Nv==1)
              {
                for(int i=0;i<(((BlockTriDiagonalMatrixType*)(GvalBath->data))->BlockSize);i++)
                {
                  lua_pushnumber(L, i);
                  lua_rawseti(L, -2, i+1);
                }
              }
              else
              {
                for(int i=0;i<(((BlockTriDiagonalMatrixType*)(GvalBath->data))->BlockSize);i++)
                {
                  lua_pushnumber(L, i+(((BlockTriDiagonalMatrixType*)(GvalBath->data))->BlockSize)*(Options.Nv));
                  lua_rawseti(L, -2, i+1);
                }
              }
              lua_settable(L, -3);
            }
            else // push nil for G valence Bath
            {
              lua_pushstring(L, "GvalBath");
              lua_pushnil(L);
              lua_settable(L, -3);
              lua_pushstring(L, "GvalBathIndex");
              lua_createtable(L,0,0);
              lua_settable(L, -3);
            }
            if(MNc>Options.Nc)
            {
              ResponsefunctionType *GconBath;
              lua_pushstring(L, "GconBath");
              GconBath = (ResponsefunctionType *) lua_newuserdata(L, sizeof(ResponsefunctionType));
              luaL_setmetatable(L, "ResponseFunction_Type");
              GconBath->type = 't';
              GconBath->Block=true;
              GconBath->NSuggestedGridPoints=0;
              GconBath->SuggestedGrid=NULL;
              GconBath->data = malloc(sizeof(BlockTriDiagonalMatrixType));
              BlockTriDiagonalMatrixGetTailAndAddSpin(&(((BlockNaturalImpurityOrbitalMatrixType*)(M->data))->con), (BlockTriDiagonalMatrixType*)(GconBath->data), MNc-Options.Nc);
              BlockTriDiagonalMatrixSetA0to0((BlockTriDiagonalMatrixType*)(GconBath->data));
              lua_settable(L, -3);
              lua_pushstring(L, "GconBathIndex");
              lua_createtable(L,((BlockTriDiagonalMatrixType*)(GconBath->data))->BlockSize,0);
              if(Options.Nc==1)
              {
                for(int i=0;i<(((BlockTriDiagonalMatrixType*)(GconBath->data))->BlockSize);i++)
                {
                  lua_pushnumber(L, i+((BlockTriDiagonalMatrixType*)(GconBath->data))->BlockSize);
                  lua_rawseti(L, -2, i+1);
                }
              }
              else
              {
                if(Min(Options.Nv,MNv)==1)
                {
                  for(int i=0;i<(((BlockTriDiagonalMatrixType*)(GconBath->data))->BlockSize);i++)
                  {
                    lua_pushnumber(L, i+(((BlockTriDiagonalMatrixType*)(GconBath->data))->BlockSize)*(Options.Nc));
                    lua_rawseti(L, -2, i+1);
                  }
                }
                else
                {
                  for(int i=0;i<(((BlockTriDiagonalMatrixType*)(GconBath->data))->BlockSize);i++)
                  {
                    lua_pushnumber(L, i+(((BlockTriDiagonalMatrixType*)(GconBath->data))->BlockSize)*(Max(1,Min(Options.Nv,MNv))-1+(Options.Nc)));
                    lua_rawseti(L, -2, i+1);
                  }
                }
              }
              lua_settable(L, -3);
            }
            else // push nil for G conduction Bath
            {
              lua_pushstring(L, "GconBath");
              lua_pushnil(L);
              lua_settable(L, -3);
              lua_pushstring(L, "GconBathIndex");
              lua_createtable(L,0,0);
              lua_settable(L, -3);
            }
          }
          else // add spin true and block false
          {
            if(MNv>Options.Nv)
            {
              ResponsefunctionType *GvalBath;
              lua_pushstring(L, "GvalBath");
              GvalBath = (ResponsefunctionType *) lua_newuserdata(L, sizeof(ResponsefunctionType));
              luaL_setmetatable(L, "ResponseFunction_Type");
              GvalBath->type = 't';
              GvalBath->Block=true;
              GvalBath->NSuggestedGridPoints=0;
              GvalBath->SuggestedGrid=NULL;
              GvalBath->data = malloc(sizeof(BlockTriDiagonalMatrixType));
              TriDiagonalMatrixPartAddSpin((BlockTriDiagonalMatrixType*)(GvalBath->data), ((BADoubleTriDiagonalMatrixType*)(M->data))->val, Options.Nv+1, MNv);
              BlockTriDiagonalMatrixSetA0to0((BlockTriDiagonalMatrixType*)(GvalBath->data));
              lua_settable(L, -3);
              lua_pushstring(L, "GvalBathIndex");
              lua_createtable(L,((BlockTriDiagonalMatrixType*)(GvalBath->data))->BlockSize,0);
              if(Options.Nv==1)
              {
                for(int i=0;i<(((BlockTriDiagonalMatrixType*)(GvalBath->data))->BlockSize);i++)
                {
                  lua_pushnumber(L, i);
                  lua_rawseti(L, -2, i+1);
                }
              }
              else
              {
                for(int i=0;i<(((BlockTriDiagonalMatrixType*)(GvalBath->data))->BlockSize);i++)
                {
                  lua_pushnumber(L, i+(((BlockTriDiagonalMatrixType*)(GvalBath->data))->BlockSize)*(Options.Nv));
                  lua_rawseti(L, -2, i+1);
                }
              }
              lua_settable(L, -3);
            }
            else // push nil for G valence Bath
            {
              lua_pushstring(L, "GvalBath");
              lua_pushnil(L);
              lua_settable(L, -3);
              lua_pushstring(L, "GvalBathIndex");
              lua_createtable(L,0,0);
              lua_settable(L, -3);
            }
            if(MNc>Options.Nc)
            {
              ResponsefunctionType *GconBath;
              lua_pushstring(L, "GconBath");
              GconBath = (ResponsefunctionType *) lua_newuserdata(L, sizeof(ResponsefunctionType));
              luaL_setmetatable(L, "ResponseFunction_Type");
              GconBath->type = 't';
              GconBath->Block=true;
              GconBath->NSuggestedGridPoints=0;
              GconBath->SuggestedGrid=NULL;
              GconBath->data = malloc(sizeof(BlockTriDiagonalMatrixType));
              TriDiagonalMatrixPartAddSpin((BlockTriDiagonalMatrixType*)(GconBath->data), ((BADoubleTriDiagonalMatrixType*)(M->data))->con, Options.Nc+1, MNc);
              BlockTriDiagonalMatrixSetA0to0((BlockTriDiagonalMatrixType*)(GconBath->data));
              lua_settable(L, -3);
              lua_pushstring(L, "GconBathIndex");
              lua_createtable(L,((BlockTriDiagonalMatrixType*)(GconBath->data))->BlockSize,0);
              if(Options.Nc==1)
              {
                for(int i=0;i<(((BlockTriDiagonalMatrixType*)(GconBath->data))->BlockSize);i++)
                {
                  lua_pushnumber(L, i+((BlockTriDiagonalMatrixType*)(GconBath->data))->BlockSize);
                  lua_rawseti(L, -2, i+1);
                }
              }
              else
              {
                if(Min(Options.Nv,MNv)==1)
                {
                  for(int i=0;i<(((BlockTriDiagonalMatrixType*)(GconBath->data))->BlockSize);i++)
                  {
                    lua_pushnumber(L, i+(((BlockTriDiagonalMatrixType*)(GconBath->data))->BlockSize)*(Options.Nc));
                    lua_rawseti(L, -2, i+1);
                  }
                }
                else
                {
                  for(int i=0;i<(((BlockTriDiagonalMatrixType*)(GconBath->data))->BlockSize);i++)
                  {
                    lua_pushnumber(L, i+(((BlockTriDiagonalMatrixType*)(GconBath->data))->BlockSize)*(Max(1,Min(Options.Nv,MNv))-1+(Options.Nc)));
                    lua_rawseti(L, -2, i+1);
                  }
                }
              }
              lua_settable(L, -3);
            }
            else // push nil for G conduction Bath
            {
              lua_pushstring(L, "GconBath");
              lua_pushnil(L);
              lua_settable(L, -3);
              lua_pushstring(L, "GconBathIndex");
              lua_createtable(L,0,0);
              lua_settable(L, -3);
            }
          }
        }
        else // no spin added
        {
          if(M->Block==true) // no spin and block true
          {
            if(MNv>Options.Nv)
            {
              ResponsefunctionType *GvalBath;
              lua_pushstring(L, "GvalBath");
              GvalBath = (ResponsefunctionType *) lua_newuserdata(L, sizeof(ResponsefunctionType));
              luaL_setmetatable(L, "ResponseFunction_Type");
              GvalBath->type = 't';
              GvalBath->Block=true;
              GvalBath->NSuggestedGridPoints=0;
              GvalBath->SuggestedGrid=NULL;
              GvalBath->data = malloc(sizeof(BlockTriDiagonalMatrixType));
              BlockTriDiagonalMatrixGetTail(&(((BlockNaturalImpurityOrbitalMatrixType*)(M->data))->val), (BlockTriDiagonalMatrixType*)(GvalBath->data), MNv-Options.Nv);
              BlockTriDiagonalMatrixSetA0to0((BlockTriDiagonalMatrixType*)(GvalBath->data));
              lua_settable(L, -3);
              lua_pushstring(L, "GvalBathIndex");
              lua_createtable(L,((BlockTriDiagonalMatrixType*)(GvalBath->data))->BlockSize,0);
              if(Options.Nv==1)
              {
                for(int i=0;i<(((BlockTriDiagonalMatrixType*)(GvalBath->data))->BlockSize);i++)
                {
                  lua_pushnumber(L, i);
                  lua_rawseti(L, -2, i+1);
                }
              }
              else
              {
                for(int i=0;i<(((BlockTriDiagonalMatrixType*)(GvalBath->data))->BlockSize);i++)
                {
                  lua_pushnumber(L, i+(((BlockTriDiagonalMatrixType*)(GvalBath->data))->BlockSize)*(Options.Nv));
                  lua_rawseti(L, -2, i+1);
                }
              }
              lua_settable(L, -3);
            }
            else // push nil for G valence Bath
            {
              lua_pushstring(L, "GvalBath");
              lua_pushnil(L);
              lua_settable(L, -3);
              lua_pushstring(L, "GvalBathIndex");
              lua_createtable(L,0,0);
              lua_settable(L, -3);
            }
            if(MNc>Options.Nc)
            {
              ResponsefunctionType *GconBath;
              lua_pushstring(L, "GconBath");
              GconBath = (ResponsefunctionType *) lua_newuserdata(L, sizeof(ResponsefunctionType));
              luaL_setmetatable(L, "ResponseFunction_Type");
              GconBath->type = 't';
              GconBath->Block=true;
              GconBath->NSuggestedGridPoints=0;
              GconBath->SuggestedGrid=NULL;
              GconBath->data = malloc(sizeof(BlockTriDiagonalMatrixType));
              BlockTriDiagonalMatrixGetTail(&(((BlockNaturalImpurityOrbitalMatrixType*)(M->data))->con), (BlockTriDiagonalMatrixType*)(GconBath->data), MNc-Options.Nc);
              BlockTriDiagonalMatrixSetA0to0((BlockTriDiagonalMatrixType*)(GconBath->data));
              lua_settable(L, -3);
              lua_pushstring(L, "GconBathIndex");
              lua_createtable(L,((BlockTriDiagonalMatrixType*)(GconBath->data))->BlockSize,0);
              if(Options.Nc==1)
              {
                for(int i=0;i<(((BlockTriDiagonalMatrixType*)(GconBath->data))->BlockSize);i++)
                {
                  lua_pushnumber(L, i+((BlockTriDiagonalMatrixType*)(GconBath->data))->BlockSize);
                  lua_rawseti(L, -2, i+1);
                }
              }
              else
              {
                if(Min(Options.Nv,MNv)==1)
                {
                  for(int i=0;i<(((BlockTriDiagonalMatrixType*)(GconBath->data))->BlockSize);i++)
                  {
                    lua_pushnumber(L, i+(((BlockTriDiagonalMatrixType*)(GconBath->data))->BlockSize)*(Options.Nc));
                    lua_rawseti(L, -2, i+1);
                  }
                }
                else
                {
                  for(int i=0;i<(((BlockTriDiagonalMatrixType*)(GconBath->data))->BlockSize);i++)
                  {
                    lua_pushnumber(L, i+(((BlockTriDiagonalMatrixType*)(GconBath->data))->BlockSize)*(Max(1,Min(Options.Nv,MNv))-1+(Options.Nc)));
                    lua_rawseti(L, -2, i+1);
                  }
                }
              }
              lua_settable(L, -3);
            }
            else // push nil for G conduction Bath
            {
              lua_pushstring(L, "GconBath");
              lua_pushnil(L);
              lua_settable(L, -3);
              lua_pushstring(L, "GconBathIndex");
              lua_createtable(L,0,0);
              lua_settable(L, -3);
            }
          }
          else // no spin and block false
          {
            if(MNv>Options.Nv)
            {
              ResponsefunctionType *GvalBath;
              lua_pushstring(L, "GvalBath");
              GvalBath = (ResponsefunctionType *) lua_newuserdata(L, sizeof(ResponsefunctionType));
              luaL_setmetatable(L, "ResponseFunction_Type");
              GvalBath->type = 't';
              GvalBath->Block=false;
              GvalBath->NSuggestedGridPoints=0;
              GvalBath->SuggestedGrid=NULL;
              GvalBath->data = malloc(sizeof(TriDiagonalMatrixType));
              TriDiagonalMatrixPart((TriDiagonalMatrixType*)(GvalBath->data), ((BADoubleTriDiagonalMatrixType*)(M->data))->val, Options.Nv+1, MNv);
              ((TriDiagonalMatrixType*)(GvalBath->data))->a[0]=0;
              lua_settable(L, -3);
              lua_pushstring(L, "GvalBathIndex");
              lua_createtable(L,1,0);
              if(Options.Nv==1)
              {
                lua_pushnumber(L, 0);
                lua_rawseti(L, -2, 1);
              }
              else
              {
                lua_pushnumber(L, (Options.Nv));
                lua_rawseti(L, -2, 1);
              }
              lua_settable(L, -3);            }
            else // push nil for G valence Bath
            {
              lua_pushstring(L, "GvalBath");
              lua_pushnil(L);
              lua_settable(L, -3);
              lua_pushstring(L, "GvalBathIndex");
              lua_createtable(L,0,0);
              lua_settable(L, -3);
            }
            if(MNc>Options.Nc)
            {
              ResponsefunctionType *GconBath;
              lua_pushstring(L, "GconBath");
              GconBath = (ResponsefunctionType *) lua_newuserdata(L, sizeof(ResponsefunctionType));
              luaL_setmetatable(L, "ResponseFunction_Type");
              GconBath->type = 't';
              GconBath->Block=false;
              GconBath->NSuggestedGridPoints=0;
              GconBath->SuggestedGrid=NULL;
              GconBath->data = malloc(sizeof(TriDiagonalMatrixType));
              TriDiagonalMatrixPart((TriDiagonalMatrixType*)(GconBath->data), ((BADoubleTriDiagonalMatrixType*)(M->data))->con, Options.Nc+1, MNc);
              ((TriDiagonalMatrixType*)(GconBath->data))->a[0]=0;
              lua_settable(L, -3);
              lua_pushstring(L, "GconBathIndex");
              lua_createtable(L,1,0);
              if(Options.Nc==1)
              {
                lua_pushnumber(L, 1);
                lua_rawseti(L, -2, 1);
              }
              else
              {
                if(Min(Options.Nv,MNv)==1)
                {
                  lua_pushnumber(L, Options.Nc);
                }
                else
                {
                  lua_pushnumber(L, Max(1,Min(Options.Nv,MNv))-1+(Options.Nc));
                }
                lua_rawseti(L, -2, 1);
              }
              lua_settable(L, -3);
            }
            else // push nil for G conduction Bath
            {
              lua_pushstring(L, "GconBath");
              lua_pushnil(L);
              lua_settable(L, -3);
              lua_pushstring(L, "GconBathIndex");
              lua_createtable(L,0,0);
              lua_settable(L, -3);
            }
          }
        }
        FreeResponsefunctionIncludingOriginalPointer(M);
        return 3;
      }
      else
      {
        break;
      }
    }
    case 'D':
    case 'd':
    {
      // returns the operator, and the table of bitmasks, but aditionally, when Nv or Nc are set also
      // a table containing 2 elements:
      //    GvalBath
      //    GconBath
      //    These contain the tails of the vanlence and conduction chain
      //
      if(M->Block==true)
      {
        // copy and shorten M->data if Nv and Nc set
        MNv = ((BlockDoubleTriDiagonalMatrixType*)(M->data))->val.N;
        MNc = ((BlockDoubleTriDiagonalMatrixType*)(M->data))->con.N;
        if(MNv>Options.Nv) ((BlockDoubleTriDiagonalMatrixType*)(M->data))->val.N=Options.Nv;
        if(MNc>Options.Nc) ((BlockDoubleTriDiagonalMatrixType*)(M->data))->con.N=Options.Nc;
        if(BlockDoubleTriDiagonalMatrixToImpurityHamiltonian(O, (BlockDoubleTriDiagonalMatrixType*)(M->data), Options.AddSpin, Options.FermionShift, 0))
        {
          if(MNv>Options.Nv) ((BlockDoubleTriDiagonalMatrixType*)(M->data))->val.N=MNv;
          if(MNc>Options.Nc) ((BlockDoubleTriDiagonalMatrixType*)(M->data))->con.N=MNc;
          return luaL_error(L, "Error in Responsefunction.ToOperator() BlockDoubleTriDiagonalMatrixToImpurityHamiltonian is not implemented\n");
        }
      }
      else
      {
        // copy and shorten M->data if Nv and Nc set
        MNv = ((DoubleTriDiagonalMatrixType*)M->data)->val.N;
        MNc = ((DoubleTriDiagonalMatrixType*)M->data)->con.N;
        if(MNv>Options.Nv) ((DoubleTriDiagonalMatrixType*)M->data)->val.N=Options.Nv;
        if(MNc>Options.Nc) ((DoubleTriDiagonalMatrixType*)M->data)->con.N=Options.Nc;
        if(DoubleTriDiagonalMatrixToImpurityHamiltonian(O, (DoubleTriDiagonalMatrixType*)M->data, Options.AddSpin, Options.FermionShift, 0))
        {
          if(MNv>Options.Nv) ((DoubleTriDiagonalMatrixType*)M->data)->val.N=MNv;
          if(MNc>Options.Nc) ((DoubleTriDiagonalMatrixType*)M->data)->con.N=MNc;
          return luaL_error(L, "Error in Responsefunction.ToOperator() DoubleTriDiagonalMatrixToImpurityHamiltonian is not implemented\n");
        }
      }
      if((MNv>Options.Nv)||(MNc>Options.Nc))
      {
        ResponseFunctionToLuaBitMask(L, M, Options.AddSpin, Options.FermionShift,false);
        // undo changes to M
        if(M->Block==true)
        {
          // expand M->data if Nv and Nc set
          if(MNv>Options.Nv) ((BlockDoubleTriDiagonalMatrixType*)(M->data))->val.N=MNv;
          if(MNc>Options.Nc) ((BlockDoubleTriDiagonalMatrixType*)(M->data))->con.N=MNc;
        }
        else
        {
          // expand M->data if Nv and Nc set
          if(MNv>Options.Nv) ((DoubleTriDiagonalMatrixType*)M->data)->val.N=MNv;
          if(MNc>Options.Nc) ((DoubleTriDiagonalMatrixType*)M->data)->con.N=MNc;
        }
        // push Gv, and, Gc to stack
        lua_createtable(L,0,4);
        // nested if then else statements to select cases for AddSpin, Block and Nv,Nc set
        if(Options.AddSpin==true)
        {
          if(M->Block==true) // add spin and block true
          {
            if(MNv>Options.Nv)
            {
              ResponsefunctionType *GvalBath;
              lua_pushstring(L, "GvalBath");
              GvalBath = (ResponsefunctionType *) lua_newuserdata(L, sizeof(ResponsefunctionType));
              luaL_setmetatable(L, "ResponseFunction_Type");
              GvalBath->type = 't';
              GvalBath->Block=true;
              GvalBath->NSuggestedGridPoints=0;
              GvalBath->SuggestedGrid=NULL;
              GvalBath->data = malloc(sizeof(BlockTriDiagonalMatrixType));
              BlockTriDiagonalMatrixGetTailAndAddSpin(&(((BlockDoubleTriDiagonalMatrixType*)(M->data))->val), (BlockTriDiagonalMatrixType*)(GvalBath->data), MNv-Options.Nv);
              BlockTriDiagonalMatrixSetA0to0((BlockTriDiagonalMatrixType*)(GvalBath->data));
              lua_settable(L, -3);
              lua_pushstring(L, "GvalBathIndex");
              lua_createtable(L,((BlockTriDiagonalMatrixType*)(GvalBath->data))->BlockSize,0);
              for(int i=0;i<(((BlockTriDiagonalMatrixType*)(GvalBath->data))->BlockSize);i++)
              {
                lua_pushnumber(L, i+(((BlockTriDiagonalMatrixType*)(GvalBath->data))->BlockSize)*(Options.Nv));
                lua_rawseti(L, -2, i+1);
              }
              lua_settable(L, -3);
            }
            else // push nil for G valence Bath
            {
              lua_pushstring(L, "GvalBath");
              lua_pushnil(L);
              lua_settable(L, -3);
              lua_pushstring(L, "GvalBathIndex");
              lua_createtable(L,0,0);
              lua_settable(L, -3);
            }
            if(MNc>Options.Nc)
            {
              ResponsefunctionType *GconBath;
              lua_pushstring(L, "GconBath");
              GconBath = (ResponsefunctionType *) lua_newuserdata(L, sizeof(ResponsefunctionType));
              luaL_setmetatable(L, "ResponseFunction_Type");
              GconBath->type = 't';
              GconBath->Block=true;
              GconBath->NSuggestedGridPoints=0;
              GconBath->SuggestedGrid=NULL;
              GconBath->data = malloc(sizeof(BlockTriDiagonalMatrixType));
              BlockTriDiagonalMatrixGetTailAndAddSpin(&(((BlockDoubleTriDiagonalMatrixType*)(M->data))->con), (BlockTriDiagonalMatrixType*)(GconBath->data), MNc-Options.Nc);
              BlockTriDiagonalMatrixSetA0to0((BlockTriDiagonalMatrixType*)(GconBath->data));
              lua_settable(L, -3);
              lua_pushstring(L, "GconBathIndex");
              lua_createtable(L,((BlockTriDiagonalMatrixType*)(GconBath->data))->BlockSize,0);
              if(Options.Nc==0)
              {
                for(int i=0;i<(((BlockTriDiagonalMatrixType*)(GconBath->data))->BlockSize);i++)
                {
                  lua_pushnumber(L, i);
                  lua_rawseti(L, -2, i+1);
                }
              }
              else
              {
                for(int i=0;i<(((BlockTriDiagonalMatrixType*)(GconBath->data))->BlockSize);i++)
                {
                  lua_pushnumber(L, i+(((BlockTriDiagonalMatrixType*)(GconBath->data))->BlockSize)*(Min(MNv,Options.Nv)+Options.Nc));
                  lua_rawseti(L, -2, i+1);
                }
              }
              lua_settable(L, -3);
            }
            else // push nil for G conduction Bath
            {
              lua_pushstring(L, "GconBath");
              lua_pushnil(L);
              lua_settable(L, -3);
              lua_pushstring(L, "GconBathIndex");
              lua_createtable(L,0,0);
              lua_settable(L, -3);
            }
          }
          else // add spin true and block false
          {
            if(MNv>Options.Nv)
            {
              ResponsefunctionType *GvalBath;
              lua_pushstring(L, "GvalBath");
              GvalBath = (ResponsefunctionType *) lua_newuserdata(L, sizeof(ResponsefunctionType));
              luaL_setmetatable(L, "ResponseFunction_Type");
              GvalBath->type = 't';
              GvalBath->Block=true;
              GvalBath->NSuggestedGridPoints=0;
              GvalBath->SuggestedGrid=NULL;
              GvalBath->data = malloc(sizeof(BlockTriDiagonalMatrixType));
              TriDiagonalMatrixPartAddSpin((BlockTriDiagonalMatrixType*)(GvalBath->data), ((DoubleTriDiagonalMatrixType*)(M->data))->val, Options.Nv+1, MNv);
              BlockTriDiagonalMatrixSetA0to0((BlockTriDiagonalMatrixType*)(GvalBath->data));
              lua_settable(L, -3);
              lua_pushstring(L, "GvalBathIndex");
              lua_createtable(L,((BlockTriDiagonalMatrixType*)(GvalBath->data))->BlockSize,0);
              for(int i=0;i<(((BlockTriDiagonalMatrixType*)(GvalBath->data))->BlockSize);i++)
              {
                lua_pushnumber(L, i+(((BlockTriDiagonalMatrixType*)(GvalBath->data))->BlockSize)*(Options.Nv));
                lua_rawseti(L, -2, i+1);
              }
              lua_settable(L, -3);
            }
            else // push nil for G valence Bath
            {
              lua_pushstring(L, "GvalBath");
              lua_pushnil(L);
              lua_settable(L, -3);
              lua_pushstring(L, "GvalBathIndex");
              lua_createtable(L,0,0);
              lua_settable(L, -3);
            }
            if(MNc>Options.Nc)
            {
              ResponsefunctionType *GconBath;
              lua_pushstring(L, "GconBath");
              GconBath = (ResponsefunctionType *) lua_newuserdata(L, sizeof(ResponsefunctionType));
              luaL_setmetatable(L, "ResponseFunction_Type");
              GconBath->type = 't';
              GconBath->Block=true;
              GconBath->NSuggestedGridPoints=0;
              GconBath->SuggestedGrid=NULL;
              GconBath->data = malloc(sizeof(BlockTriDiagonalMatrixType));
              TriDiagonalMatrixPartAddSpin((BlockTriDiagonalMatrixType*)(GconBath->data), ((DoubleTriDiagonalMatrixType*)(M->data))->con, Options.Nc+1, MNc);
              BlockTriDiagonalMatrixSetA0to0((BlockTriDiagonalMatrixType*)(GconBath->data));
              lua_settable(L, -3);
              lua_pushstring(L, "GconBathIndex");
              lua_createtable(L,((BlockTriDiagonalMatrixType*)(GconBath->data))->BlockSize,0);
              if(Options.Nc==0)
              {
                for(int i=0;i<(((BlockTriDiagonalMatrixType*)(GconBath->data))->BlockSize);i++)
                {
                  lua_pushnumber(L, i);
                  lua_rawseti(L, -2, i+1);
                }
              }
              else
              {
                for(int i=0;i<(((BlockTriDiagonalMatrixType*)(GconBath->data))->BlockSize);i++)
                {
                  lua_pushnumber(L, i+(((BlockTriDiagonalMatrixType*)(GconBath->data))->BlockSize)*(Min(MNv,Options.Nv)+Options.Nc));
                  lua_rawseti(L, -2, i+1);
                }
              }
              lua_settable(L, -3);
            }
            else // push nil for G conduction Bath
            {
              lua_pushstring(L, "GconBath");
              lua_pushnil(L);
              lua_settable(L, -3);
              lua_pushstring(L, "GconBathIndex");
              lua_createtable(L,0,0);
              lua_settable(L, -3);
            }
          }
        }
        else // no spin added
        {
          if(M->Block==true) // no spin and block true
          {
            if(MNv>Options.Nv)
            {
              ResponsefunctionType *GvalBath;
              lua_pushstring(L, "GvalBath");
              GvalBath = (ResponsefunctionType *) lua_newuserdata(L, sizeof(ResponsefunctionType));
              luaL_setmetatable(L, "ResponseFunction_Type");
              GvalBath->type = 't';
              GvalBath->Block=true;
              GvalBath->NSuggestedGridPoints=0;
              GvalBath->SuggestedGrid=NULL;
              GvalBath->data = malloc(sizeof(BlockTriDiagonalMatrixType));
              BlockTriDiagonalMatrixGetTail(&(((BlockDoubleTriDiagonalMatrixType*)(M->data))->val), (BlockTriDiagonalMatrixType*)(GvalBath->data), MNv-Options.Nv);
              BlockTriDiagonalMatrixSetA0to0((BlockTriDiagonalMatrixType*)(GvalBath->data));
              lua_settable(L, -3);
              lua_pushstring(L, "GvalBathIndex");
              lua_createtable(L,((BlockTriDiagonalMatrixType*)(GvalBath->data))->BlockSize,0);
              for(int i=0;i<(((BlockTriDiagonalMatrixType*)(GvalBath->data))->BlockSize);i++)
              {
                lua_pushnumber(L, i+(((BlockTriDiagonalMatrixType*)(GvalBath->data))->BlockSize)*(Options.Nv));
                lua_rawseti(L, -2, i+1);
              }
              lua_settable(L, -3);
            }
            else // push nil for G valence Bath
            {
              lua_pushstring(L, "GvalBath");
              lua_pushnil(L);
              lua_settable(L, -3);
              lua_pushstring(L, "GvalBathIndex");
              lua_createtable(L,0,0);
              lua_settable(L, -3);
            }
            if(MNc>Options.Nc)
            {
              ResponsefunctionType *GconBath;
              lua_pushstring(L, "GconBath");
              GconBath = (ResponsefunctionType *) lua_newuserdata(L, sizeof(ResponsefunctionType));
              luaL_setmetatable(L, "ResponseFunction_Type");
              GconBath->type = 't';
              GconBath->Block=true;
              GconBath->NSuggestedGridPoints=0;
              GconBath->SuggestedGrid=NULL;
              GconBath->data = malloc(sizeof(BlockTriDiagonalMatrixType));
              BlockTriDiagonalMatrixGetTail(&(((BlockDoubleTriDiagonalMatrixType*)(M->data))->con), (BlockTriDiagonalMatrixType*)(GconBath->data), MNc-Options.Nc);
              BlockTriDiagonalMatrixSetA0to0((BlockTriDiagonalMatrixType*)(GconBath->data));
              lua_settable(L, -3);
              lua_pushstring(L, "GconBathIndex");
              lua_createtable(L,((BlockTriDiagonalMatrixType*)(GconBath->data))->BlockSize,0);
              if(Options.Nc==0)
              {
                for(int i=0;i<(((BlockTriDiagonalMatrixType*)(GconBath->data))->BlockSize);i++)
                {
                  lua_pushnumber(L, i);
                  lua_rawseti(L, -2, i+1);
                }
              }
              else
              {
                for(int i=0;i<(((BlockTriDiagonalMatrixType*)(GconBath->data))->BlockSize);i++)
                {
                  lua_pushnumber(L, i+(((BlockTriDiagonalMatrixType*)(GconBath->data))->BlockSize)*(Min(MNv,Options.Nv)+Options.Nc));
                  lua_rawseti(L, -2, i+1);
                }
              }
              lua_settable(L, -3);
            }
            else // push nil for G conduction Bath
            {
              lua_pushstring(L, "GconBath");
              lua_pushnil(L);
              lua_settable(L, -3);
              lua_pushstring(L, "GconBathIndex");
              lua_createtable(L,0,0);
              lua_settable(L, -3);
            }
          }
          else // no spin and block false
          {
            if(MNv>Options.Nv)
            {
              ResponsefunctionType *GvalBath;
              lua_pushstring(L, "GvalBath");
              GvalBath = (ResponsefunctionType *) lua_newuserdata(L, sizeof(ResponsefunctionType));
              luaL_setmetatable(L, "ResponseFunction_Type");
              GvalBath->type = 't';
              GvalBath->Block=false;
              GvalBath->NSuggestedGridPoints=0;
              GvalBath->SuggestedGrid=NULL;
              GvalBath->data = malloc(sizeof(TriDiagonalMatrixType));
              TriDiagonalMatrixPart((TriDiagonalMatrixType*)(GvalBath->data), ((DoubleTriDiagonalMatrixType*)(M->data))->val, Options.Nv+1, MNv);
              ((TriDiagonalMatrixType*)(GvalBath->data))->a[0]=0;
              lua_settable(L, -3);
              lua_pushstring(L, "GvalBathIndex");
              lua_createtable(L,((BlockTriDiagonalMatrixType*)(GvalBath->data))->BlockSize,0);
              lua_pushnumber(L, Options.Nv);
              lua_rawseti(L, -2, 1);
              lua_settable(L, -3);            }
            else // push nil for G valence Bath
            {
              lua_pushstring(L, "GvalBath");
              lua_pushnil(L);
              lua_settable(L, -3);
              lua_pushstring(L, "GvalBathIndex");
              lua_createtable(L,0,0);
              lua_settable(L, -3);
            }
            if(MNc>Options.Nc)
            {
              ResponsefunctionType *GconBath;
              lua_pushstring(L, "GconBath");
              GconBath = (ResponsefunctionType *) lua_newuserdata(L, sizeof(ResponsefunctionType));
              luaL_setmetatable(L, "ResponseFunction_Type");
              GconBath->type = 't';
              GconBath->Block=false;
              GconBath->NSuggestedGridPoints=0;
              GconBath->SuggestedGrid=NULL;
              GconBath->data = malloc(sizeof(TriDiagonalMatrixType));
              TriDiagonalMatrixPart((TriDiagonalMatrixType*)(GconBath->data), ((DoubleTriDiagonalMatrixType*)(M->data))->con, Options.Nc+1, MNc);
              ((TriDiagonalMatrixType*)(GconBath->data))->a[0]=0;
              lua_settable(L, -3);
              lua_pushstring(L, "GconBathIndex");
              lua_createtable(L,((BlockTriDiagonalMatrixType*)(GconBath->data))->BlockSize,0);
              if(Options.Nc==0)
              {
                lua_pushnumber(L, Options.Nc);
              }
              else
              {
                lua_pushnumber(L, Min(Options.Nv,MNv)+Options.Nc);
              }
              lua_rawseti(L, -2, 1);
              lua_settable(L, -3);
            }
            else // push nil for G conduction Bath
            {
              lua_pushstring(L, "GconBath");
              lua_pushnil(L);
              lua_settable(L, -3);
              lua_pushstring(L, "GconBathIndex");
              lua_createtable(L,0,0);
              lua_settable(L, -3);
            }
          }
        }
        FreeResponsefunctionIncludingOriginalPointer(M);
        return 3;
      }
      else
      {
        break;
      }
    }
    default:
    {
      return luaL_error(L, "Error in Responsefunction.ToOperator() unrecognized type for conversion (%c)\n",M->type);
      return 0;
    }
  }
  ResponseFunctionToLuaBitMask(L, M, Options.AddSpin, Options.FermionShift,false);
  FreeResponsefunctionIncludingOriginalPointer(M);
  return 2;
}


/* Calculate the density matrix of a Green's function given as a response function.
 * Input on stack:
 *   G : ResponseFunction_Type
 *   Options
 *       mu  -> chemical potential                      standard the chemical potential stored in G
 *       kBT -> temperature times the Boltzman constant  standard 0
 * Returns a BlockSize by BlockSize matrix (1 by 1 for a scalar Green's function)
 */
static int LuaResponseFunctionCalculateRho(lua_State *L)
{
  uint32_t len;
  ResponsefunctionType *A;
  CompactMatrixType *rho;
  struct
  {
    double mu;
    double kbT;
  } Options;

  len = (uint32_t)lua_gettop(L);
  if(!(len==1||len==2))
  {
    return luaL_error(L, "Error CalculateRho: called with %d arguments, 1 or 2 expected (a Green's function and options).\n",len);
  }
  A = luaL_checkResponsefunction(L,1);
  Options.mu = ResponseFunctionGetMu(A);
  Options.kbT= 0;
  LuaOptionsType LuaCalculateRhoOptions[]={
    {"mu", "double", &(Options.mu)},
    {"kbT","double", &(Options.kbT)},
    {"kBT","double", &(Options.kbT)},
    {NULL, NULL, NULL}};
  if(len==2)
  {
    LuaCheckOptions(L, 2, LuaCalculateRhoOptions);
  }
  /*  the density matrix first and the userdata to hold it after: a failed
   *  ResponseFunctionToDensityMatrix leaves nothing behind for the garbage collector to
   *  free, and the matrix does not reach lua at all */
  CompactMatrixType r;
  if(ResponseFunctionToDensityMatrix(A, &r, Options.mu, Options.kbT))
  {
    FreeResponsefunctionIncludingOriginalPointer(A);
    return luaL_error(L, "Error in CalculateRho: could not calculate the density matrix\n");
  }
  rho = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
  luaL_setmetatable(L, "Matrix_Type");
  *rho = r;
  rho->OnLuaStack=true;
  FreeResponsefunctionIncludingOriginalPointer(A);
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
  Options.Deflate = true;
  Options.SingularValue = DBL_EPSILON;
  LuaOptionsType LuaOptions[]={
    {"Deflate",      "bool",  &(Options.Deflate)},
    {"SingularValue","double",&(Options.SingularValue)},
    {NULL, NULL, NULL}};
  if(LuaTestOptionsList(L, -1))
  {
    LuaCheckOptions(L, -1, LuaOptions);
  }
  B = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
  A = luaL_checkResponsefunction(L,1);
  ResponseFunctionToCompactMatrix(A, B, Options.Deflate, Options.SingularValue);
  luaL_setmetatable(L, "Matrix_Type");
  B->OnLuaStack=true;
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
  A = luaL_checkResponsefunction(L,1);
  ResponseFunctionEstimateEminAndEmaxPole(A, &(Options.Emin), &(Options.Emax));
  Options.EnergyShift = 0;
  Options.GammaG = 0;
  Options.GammaL = 0;
  Options.NE     = 600;
  if(LuaTestOptionsList(L, -1))
  {
    LuaCheckOptions(L, -1, LuaOptions);
  }
  Spec = (SpectraTypeArray *) lua_newuserdata(L, sizeof(SpectraTypeArray));
  luaL_setmetatable(L, "Spectra_Type");
  /*  Set the fresh user data to an empty array before it is filled. If the call below fails we
   *  leave through a lua error and the garbage collector still runs FreeSpectraTypeArray on
   *  this user data.
   */
  Spec->NSpectra=0;
  Spec->Spec=NULL;
  Spec->Name=NULL;
  if(ResponseFunctionToSpectra(A,Spec, Options.EnergyShift, Options.Emin,Options.Emax,Options.GammaG,Options.GammaL,Options.NE))
  {
    return luaL_error(L, "Error in ToSpectra: can not create a spectrum for a response function of type %c\n",A->type);
  }
  return 1;
}


/*
 *
 */
static int LuaResponseFunctionRotate(lua_State * L)
{
  ResponsefunctionType *G, GRot;
  CompactMatrixType *U;
  double ValueNotSet = -DBL_MAX;
  struct
  {
    double SingularValue;
    double MinimalPoleWeigth;
    double MinimalPoleDistance;
    double T;
  } Options;
  Options.SingularValue       = ValueNotSet;
  Options.MinimalPoleWeigth   = ValueNotSet;
  Options.MinimalPoleDistance = ValueNotSet;
  Options.T                   = ValueNotSet;
  LuaOptionsType LuaOptions[]={
//    {"SingularValue","double",&(Options.SingularValue)},
//    {"MinimalPoleWeigth","double",&(Options.MinimalPoleWeigth)},
    {"MinimalPoleDistance","double",&(Options.MinimalPoleDistance)},
//    {"Temperature","double",&(Options.T)},
//    {"EffectivekBT","double",&(Options.T)},
    {NULL, NULL, NULL}};
  uint32_t len = (uint32_t)lua_gettop(L);
  if(LuaTestOptionsList(L, -1))
  {
    LuaCheckOptions(L, -1, LuaOptions);
    len--;
  }
  if(len!=2) return luaL_error(L, "ResponseFunction.Rotate called with %d arguments, two arguments expected and possibly a third specifying options",len);
  G = luaL_checkResponsefunction(L,1);
  U=luaL_checkCompactMatrix(L, 2);
  if(G->Block==false)
  {
    return luaL_error(L, "Error in ResponseFunction.Rotate(). Rotate rotates the elements in the matrix or Block representation of the response function. Can not rotate singel valued functions\n", len);
  }
  if(ResponseFunctionGetBlockSize(G)!=U->Nj) return luaL_error(L, "Error in ResponseFunction.Rotate(). Number of columns in roation matrix (Nj = %d) should be equal to the blocksize of the response function (%d)\n", U->Nj, ResponseFunctionGetBlockSize(G));
  if(BlockResponseFunctionRotate(G,U,&GRot,Options.T,Options.SingularValue,Options.MinimalPoleWeigth,Options.MinimalPoleDistance)!=0) return luaL_error(L, "Error in ResponseFunction.Rotate()\n", len);
  LuaPush_ResponseFunction(L, &GRot);
  FreeResponsefunction(G);
  FreeResponsefunction(&GRot);
  FreeCompactMatrix(U);
  return 1;
}

/*  compare two response functions
 *  two response functions are considered equal when evaluating them returns the same answer for all values of omega and gamma
 */
static int LuaResponseFunctionEq(lua_State *L)
{
  ResponsefunctionType *Ma, *Mb;
  uint32_t len = (uint32_t)lua_gettop(L);
  if(len!=2) return luaL_error(L, "Error in ResponseFunctionSub: %d argument(s), 2 expected.\n", len);
  Ma = luaL_checkResponsefunction(L,1);
  Mb = luaL_checkResponsefunction(L,2);
  if(ResponseFunctionEq(Ma, Mb, 100000*DBL_EPSILON)==true)
  {
    lua_pushboolean(L,1);
  }
  else
  {
    lua_pushboolean(L,0);
  }
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


/* on stack a table
 *
 */
static int LuaResponseFunctionNew(lua_State *L)
{
  ResponsefunctionType *A;
  if(lua_istable(L, 1)==0)
  {
    luaL_error(L, "Error NewResponseFunction expect as input a table");
  }
  A = luaL_checkResponsefunction(L,1);
  A->type+=32;
  LuaPush_ResponseFunction(L, A);
  FreeResponsefunction(A);
  return 1;
}



/*  Returns the length of a responsefunction object if it was transformed to a table
 *  This is quite a useless object. For natural orbitals the returned value is 1, for tri, and and list of poles it returns 2.
 */
static int LuaResponseFunctionLen(lua_State *L)
{
  ResponsefunctionType *a = (ResponsefunctionType *) luaL_checkudata(L, 1, "ResponseFunction_Type");
  switch(a->type)
  {
    case 'L':
    case 'l':
    {
      lua_pushinteger(L, 2);
      return 1;
      break;
    }
    case 'A':
    case 'a':
    {
      lua_pushinteger(L, 2);
      return 1;
      break;
    }
    case 'T':
    case 't':
    {
      lua_pushinteger(L, 2);
      return 1;
      break;
    }
    case 'N':
    case 'n':
    case 'D':
    case 'd':
    {
      /*  a natural impurity orbital and a double tri-diagonal matrix have one positional
       *  element in their table form, the { a0, b0 } pair, the two chains being the named
       *  val and con. The other three have two, a list of a and a list of b. */
      lua_pushinteger(L, 1);
      return 1;
      break;
    }
    default:
    {
      return luaL_error(L, "ResponseFunction of undifined type %c\n", a->type);
      return 0;
    }
  }
}

static int LuaResponseFunctionIndexMu(lua_State *L)
{
  ResponsefunctionType *a = (ResponsefunctionType *) luaL_checkudata(L, 1, "ResponseFunction_Type");
  switch(a->type)
  {
    case 'L':
    case 'l':
    {
      if(a->Block==true)
      {
        lua_pushnumber(L, ((BlockListOfPoles*)(a->data))->mu);
        return 1;
      }
      else
      {
        lua_pushnumber(L, ((ListOfPoles*)(a->data))->mu);
        return 1;
      }
    }
    case 'A':
    case 'a':
    {
      if(a->Block==true)
      {
        lua_pushnumber(L, ((BlockAndersonMatrixType*)(a->data))->mu);
        return 1;
      }
      else
      {
        lua_pushnumber(L, ((AndersonMatrixType*)(a->data))->mu);
        return 1;
      }
    }
    case 'T':
    case 't':
    {
      if(a->Block==true)
      {
        lua_pushnumber(L, ((BlockTriDiagonalMatrixType*)(a->data))->mu);
        return 1;
      }
      else
      {
        lua_pushnumber(L, ((TriDiagonalMatrixType*)(a->data))->mu);
        return 1;
      }
    }
    case 'N':
    case 'n':
    {
      if(a->Block==true)
      {
        lua_pushnumber(L, ((BlockNaturalImpurityOrbitalMatrixType*)(a->data))->mu);
        return 1;
      }
      else
      {
        lua_pushnumber(L, ((BADoubleTriDiagonalMatrixType*)(a->data))->mu);
        return 1;
      }
    }
    case 'D':
    case 'd':
    {
      if(a->Block==true)
      {
        lua_pushnumber(L, ((BlockDoubleTriDiagonalMatrixType*)(a->data))->mu);
        return 1;
      }
      else
      {
        lua_pushnumber(L, ((DoubleTriDiagonalMatrixType*)(a->data))->mu);
        return 1;
      }
    }
    default:
    {
      return luaL_error(L, "ResponseFunction of undifined type %c\n", a->type);
    }
  }
}
static int LuaResponseFunctionIndexType(lua_State *L)
{
  ResponsefunctionType *a = (ResponsefunctionType *) luaL_checkudata(L, 1, "ResponseFunction_Type");
  switch(a->type)
  {
    case 'L':
    case 'l':
    {
      lua_pushstring(L, "ListOfPoles");
      return 1;
      break;
    }
    case 'A':
    case 'a':
    {
      lua_pushstring(L, "And");
      return 1;
      break;
    }
    case 'T':
    case 't':
    {
      lua_pushstring(L, "Tri");
      return 1;
      break;
    }
    case 'N':
    case 'n':
    {
      lua_pushstring(L, "Nat");
      return 1;
      break;
    }
    case 'D':
    case 'd':
    {
      lua_pushstring(L, "DoubleTri");
      return 1;
      break;
    }
    default:
    {
      return luaL_error(L, "ResponseFunction of undifined type %c\n", a->type);
      return 0;
    }
  }
}

static int LuaResponseFunctionIndexName(lua_State *L)
{
  ResponsefunctionType *a = (ResponsefunctionType *) luaL_checkudata(L, 1, "ResponseFunction_Type");
  switch(a->type)
  {
    case 'L':
    case 'l':
    {
      if(a->Block==true)
      {
        lua_pushstring(L, ((BlockListOfPoles*)(a->data))->Name);
        return 1;
      }
      else
      {
        lua_pushstring(L, ((ListOfPoles*)(a->data))->Name);
        return 1;
      }
    }
    case 'A':
    case 'a':
    {
      if(a->Block==true)
      {
        lua_pushstring(L, ((BlockAndersonMatrixType*)(a->data))->Name);
        return 1;
      }
      else
      {
        lua_pushstring(L, ((AndersonMatrixType*)(a->data))->Name);
        return 1;
      }
    }
    case 'T':
    case 't':
    {
      if(a->Block==true)
      {
        lua_pushstring(L, ((BlockTriDiagonalMatrixType*)(a->data))->Name);
        return 1;
      }
      else
      {
        lua_pushstring(L, ((TriDiagonalMatrixType*)(a->data))->Name);
        return 1;
      }
    }
    case 'N':
    case 'n':
    {
      if(a->Block==true)
      {
        lua_pushstring(L, ((BlockNaturalImpurityOrbitalMatrixType*)(a->data))->Name);
        return 1;
      }
      else
      {
        lua_pushstring(L, ((BADoubleTriDiagonalMatrixType*)(a->data))->Name);
        return 1;
      }
    }
    case 'D':
    case 'd':
    {
      if(a->Block==true)
      {
        lua_pushstring(L, ((BlockDoubleTriDiagonalMatrixType*)(a->data))->Name);
        return 1;
      }
      else
      {
        lua_pushstring(L, ((DoubleTriDiagonalMatrixType*)(a->data))->Name);
        return 1;
      }
    }
    default:
    {
      return luaL_error(L, "ResponseFunction of undifined type %c\n", a->type);
    }
  }
}

static int LuaResponseFunctionIndexChop(lua_State *L)
{
  ResponsefunctionType *a = (ResponsefunctionType *) luaL_checkudata(L, lua_upvalueindex(1), "ResponseFunction_Type");
  uint32_t len = (uint32_t)lua_gettop(L);
  double epsilon=10*DBL_EPSILON;
  bool Deflate = true;
  LuaOptionsType LuaOptions[]={
    {"deflate","bool",&Deflate},
    {NULL, NULL, NULL}};
  if(len!=0)
  { // check for options if options found reduce len by one
    if(LuaTestOptionsList(L, -1))
    {
      LuaCheckOptions(L, -1, LuaOptions);
      len--;
    }
  }
  if(len==1)
  {
    epsilon=luaL_checknumber(L, 1);
  }
  if(len>1)
  {
    return luaL_error(L, "Errorr, call to Chop with %d arguments, expected 0 or 1 with possible options\n",len);
  }
  ResponsefunctionChop(a, epsilon, Deflate);
  return 0;
}

/*  on stack we have a response function represented as a table
 *  we need to chop this table inline
 *  we read the table into a struct, chop the struct and then overwrite the table with the struct
 */
static int LuaResponseFunctionTableIndexChop(lua_State *L)
{
  ResponsefunctionType *a = luaL_checkResponsefunction(L, lua_upvalueindex(1));
  
  
  
  uint32_t len = (uint32_t)lua_gettop(L);
  double epsilon=10*DBL_EPSILON;
  bool Deflate = true;
  LuaOptionsType LuaOptions[]={
    {"deflate","bool",&Deflate},
    {NULL, NULL, NULL}};
  if(len!=0)
  { // check for options if options found reduce len by one
    if(LuaTestOptionsList(L, -1))
    {
      LuaCheckOptions(L, -1, LuaOptions);
      len--;
    }
  }
  if(len==1)
  {
    epsilon=luaL_checknumber(L, 1);
  }
  if(len>1)
  {
    return luaL_error(L, "Errorr, call to Chop with %d arguments, expected 0 or 1 with possible options\n",len);
  }
  ResponsefunctionChop(a, epsilon, Deflate);
  LuaPush_ResponseFunction(L, a);
  // loop over all keys and values of the new table and set the values of the old table accordingly
  lua_pushnil(L);
  while (lua_next(L, -2) != 0)
  {
    lua_pushvalue(L, -2);
    lua_insert(L, -2);
    lua_settable(L, lua_upvalueindex(1));
  }
  lua_pop(L, 1);
  FreeResponsefunctionIncludingOriginalPointer(a);
  return 0;
}

static const struct luaL_Reg ResponseFunction_Index[] = {
  {"mu",     LuaResponseFunctionIndexMu},
  {"Mu",     LuaResponseFunctionIndexMu},
  {"type",   LuaResponseFunctionIndexType},
  {"Type",   LuaResponseFunctionIndexType},
  {"name",   LuaResponseFunctionIndexName},
  {"Name",   LuaResponseFunctionIndexName},
  {NULL, NULL}
};

static const struct luaL_Reg ResponseFunctionTable_Index[] = {
  {NULL, NULL}
};

static const struct luaL_Reg ResponseFunction_Method[] = {
  {"Chop",   LuaResponseFunctionIndexChop},
  {NULL, NULL}
};

static const struct luaL_Reg ResponseFunctionTable_Method[] = {
  {"Chop",   LuaResponseFunctionTableIndexChop},
  {NULL, NULL}
};

static int LuaResponseFunctionNewIndexMu(lua_State *L)
{
  double mu = luaL_checknumber(L, 3);
  ResponsefunctionType *a = (ResponsefunctionType *) luaL_checkudata(L, 1, "ResponseFunction_Type");
  switch(a->type)
  {
    case 'L':
    case 'l':
    {
      if(a->Block==true)
      {
        ((BlockListOfPoles*)(a->data))->mu = mu;
        return 0;
      }
      else
      {
        ((ListOfPoles*)(a->data))->mu = mu;
        return 0;
      }
    }
    case 'A':
    case 'a':
    {
      if(a->Block==true)
      {
        ((BlockAndersonMatrixType*)(a->data))->mu = mu;
        return 0;
      }
      else
      {
        ((AndersonMatrixType*)(a->data))->mu = mu;
        return 0;
      }
    }
    case 'T':
    case 't':
    {
      if(a->Block==true)
      {
        ((BlockTriDiagonalMatrixType*)(a->data))->mu = mu;
        return 0;
      }
      else
      {
        ((TriDiagonalMatrixType*)(a->data))->mu = mu;
        return 0;
      }
    }
    case 'N':
    case 'n':
    {
      if(a->Block==true)
      {
        ((BlockNaturalImpurityOrbitalMatrixType*)(a->data))->mu = mu;
        return 0;
      }
      else
      {
        ((BADoubleTriDiagonalMatrixType*)(a->data))->mu = mu;
        return 0;
      }
    }
    case 'D':
    case 'd':
    {
      if(a->Block==true)
      {
        ((BlockDoubleTriDiagonalMatrixType*)(a->data))->mu = mu;
        return 0;
      }
      else
      {
        ((DoubleTriDiagonalMatrixType*)(a->data))->mu = mu;
        return 0;
      }
    }
    default:
    {
      return luaL_error(L, "ResponseFunction of undifined type %c\n", a->type);
    }
  }
}
static int LuaResponseFunctionNewIndexName(lua_State *L)
{
  const char* name = luaL_checkstring(L, 3);
  ResponsefunctionType *a = (ResponsefunctionType *) luaL_checkudata(L, 1, "ResponseFunction_Type");
  switch(a->type)
  {
    case 'L':
    case 'l':
    {
      if(a->Block==true)
      {
        strncpy(((BlockListOfPoles*)(a->data))->Name, name, MAXStrLength);
        ((BlockListOfPoles*)(a->data))->Name[MAXStrLength-1]='\0';
        return 0;
      }
      else
      {
        strncpy(((ListOfPoles*)(a->data))->Name, name, MAXStrLength);
        ((ListOfPoles*)(a->data))->Name[MAXStrLength-1]='\0';
        return 0;
      }
    }
    case 'A':
    case 'a':
    {
      if(a->Block==true)
      {
        strncpy(((BlockAndersonMatrixType*)(a->data))->Name, name, MAXStrLength);
        ((BlockAndersonMatrixType*)(a->data))->Name[MAXStrLength-1]='\0';
        return 0;
      }
      else
      {
        strncpy(((AndersonMatrixType*)(a->data))->Name, name, MAXStrLength);
        ((AndersonMatrixType*)(a->data))->Name[MAXStrLength-1]='\0';
        return 0;
      }
    }
    case 'T':
    case 't':
    {
      if(a->Block==true)
      {
        strncpy(((BlockTriDiagonalMatrixType*)(a->data))->Name, name, MAXStrLength);
        ((BlockTriDiagonalMatrixType*)(a->data))->Name[MAXStrLength-1]='\0';
        return 0;
      }
      else
      {
        strncpy(((TriDiagonalMatrixType*)(a->data))->Name, name, MAXStrLength);
        ((TriDiagonalMatrixType*)(a->data))->Name[MAXStrLength-1]='\0';
        return 0;
      }
    }
    case 'N':
    case 'n':
    {
      if(a->Block==true)
      {
        strncpy(((BlockNaturalImpurityOrbitalMatrixType*)(a->data))->Name, name, MAXStrLength);
        ((BlockNaturalImpurityOrbitalMatrixType*)(a->data))->Name[MAXStrLength-1]='\0';
        return 0;
      }
      else
      {
        strncpy(((BADoubleTriDiagonalMatrixType*)(a->data))->Name, name, MAXStrLength);
        ((BADoubleTriDiagonalMatrixType*)(a->data))->Name[MAXStrLength-1]='\0';
        return 0;
      }
    }
    case 'D':
    case 'd':
    {
      if(a->Block==true)
      {
        strncpy(((BlockDoubleTriDiagonalMatrixType*)(a->data))->Name, name, MAXStrLength);
        ((BlockDoubleTriDiagonalMatrixType*)(a->data))->Name[MAXStrLength-1]='\0';
        return 0;
      }
      else
      {
        strncpy(((DoubleTriDiagonalMatrixType*)(a->data))->Name, name, MAXStrLength);
        ((DoubleTriDiagonalMatrixType*)(a->data))->Name[MAXStrLength-1]='\0';
        return 0;
      }
    }
    default:
    {
      return luaL_error(L, "ResponseFunction of undifined type %c\n", a->type);
    }
  }
}

static const struct luaL_Reg ResponseFunction_NewIndex[] = {
  {"mu",     LuaResponseFunctionNewIndexMu},
  {"Mu",     LuaResponseFunctionNewIndexMu},
  {"name",   LuaResponseFunctionNewIndexName},
  {"Name",   LuaResponseFunctionNewIndexName},
  {NULL, NULL}
};

static int LuaResponseFunctionIndex(lua_State *L)
{
  // the index key is at position 2 in the lua stack
  const char *key = lua_tostring(L, 2);
  if (key != NULL)
  {
    // Search the mapping table for index
    for (int i = 0; ResponseFunction_Index[i].name != NULL; ++i)
    {
      if (strcmp(key, ResponseFunction_Index[i].name) == 0)
      {
        // Match found, evaluate function
        return ResponseFunction_Index[i].func(L);
      }
    }
    // Search the mapping table for methods
    for (int i = 0; ResponseFunction_Method[i].name != NULL; ++i)
    {
      if (strcmp(key, ResponseFunction_Method[i].name) == 0)
      {
        lua_pushvalue(L, 1);
        lua_pushcclosure(L, ResponseFunction_Method[i].func,1);
        return 1;
      }
    }
  }
  return luaL_error(L, "ResponseFunction has no valid index %s\n", key);
}


static int LuaResponseFunctionTableIndex(lua_State *L)
{
  // the index key is at position 2 in the lua stack
  const char *key = lua_tostring(L, 2);
  if (key != NULL)
  {
    // Search the mapping table for index
    for (int i = 0; ResponseFunctionTable_Index[i].name != NULL; ++i)
    {
      if (strcmp(key, ResponseFunctionTable_Index[i].name) == 0)
      {
        // Match found, evaluate function
        return ResponseFunctionTable_Index[i].func(L);
      }
    }
    // Search the mapping table for methods
    for (int i = 0; ResponseFunctionTable_Method[i].name != NULL; ++i)
    {
      if (strcmp(key, ResponseFunctionTable_Method[i].name) == 0)
      {
        lua_pushvalue(L, 1);
        lua_pushcclosure(L, ResponseFunctionTable_Method[i].func,1);
        return 1;
      }
    }
  }
  // Fallback to original table
  lua_rawget(L, 1);
  return 1;
}

static int LuaResponseFunctionNewIndex(lua_State *L)
{
  // the index key is at position 2 in the lua stack
  const char *key = lua_tostring(L, 2);
  if (key != NULL)
  {
    // Search the mapping table
    for (int i = 0; ResponseFunction_NewIndex[i].name != NULL; ++i)
    {
      if (strcmp(key, ResponseFunction_NewIndex[i].name) == 0)
      {
        // Match found, evaluate function
        return ResponseFunction_NewIndex[i].func(L);
      }
    }
  }
  return luaL_error(L, "ResponseFunction has no valid index %s that can be changed\n", key);
}

static const struct luaL_Reg ResponseFunction_func[] = {
  {"ChangeType", LuaChangeTypeResponseFunction},
  {"Downfold", LuaResponseFunctionDownfold},
  {"EffectiveHamiltonian",LuaResponseFunctionEffectiveHamiltonian},
  {"FromEffectiveHamiltonian",LuaResponseFunctionFromEffectiveHamiltonian},
  {"New", LuaResponseFunctionNew},
  {"ReducePoles", LuaReducePolesResponseFunction},
  {"RemovePolesInInterval",LuaRemovePolesInIntervalResponseFunction},
  {"Rotate", LuaResponseFunctionRotate},
  {"MergePolesInInterval", LuaMergePolesInIntervalResponseFunction},
  {"PlaceOnMesh", LuaPlaceOnMeshResponseFunction},
  {"ToTable", LuaResponseFunctionToTable},
  {"ToMatrix", LuaResponseFunctionToCompactMatrix},
  {"ToSpectra", LuaResponseFunctionToSpectra},
  {"ToUserdata", LuaResponseFunctionToUserdata},
  {"ToTightBinding", LuaResponseFunctionToTightBinding},
  {"ToOperator", LuaResponseFunctionToOperator},
  {"CalculateRho", LuaResponseFunctionCalculateRho},
  {"CalculateSelfEnergy", LuaCalculateSelfEnergy},
  {"CalculateHybridizationFunction", LuaCalculateHybridizationFunction},
  {"CalculateGfromG0andSigma",LuaCalculateGfromG0andSigma},
  {"CalculateGbathfromGandSigma",LuaCalculateGBathFromGandSigma},
  {"InvertEnergy", LuaResponseFunctionInvertEnergy},
  {"Transpose", LuaResponseFunctionTranspose},
  {"Elements", LuaResponseFunctionElement},
  {"CreateBitMasks", LuaResponseFunctionToBitMask},
  {"Chop", LuaChop},
  {NULL, NULL}
};

static const struct luaL_Reg ResponseFunction_meta[] = {
  {"__unm", LuaResponseFunctionUnm},
  {"__add", LuaResponseFunctionAdd},
  {"__sub", LuaResponseFunctionSub},
  {"__mul", LuaResponseFunctionMul},
  {"__div", LuaResponseFunctionDiv},
  {"__call", LuaResponseFunctionCall},
  {"__eq",   LuaResponseFunctionEq},
  {"__tostring", LuaResponseFunctionToString},
  {"__gc", LuaResponseFunctionDestroy},
  {"__index", LuaResponseFunctionIndex},
  {"__newindex", LuaResponseFunctionNewIndex},
  {"__len", LuaResponseFunctionLen},
  {NULL, NULL}
};

static const struct luaL_Reg ResponseFunctionTable_meta[] = {
  {"__unm", LuaResponseFunctionUnm},
  {"__add", LuaResponseFunctionAdd},
  {"__sub", LuaResponseFunctionSub},
  {"__mul", LuaResponseFunctionMul},
  {"__div", LuaResponseFunctionDiv},
  {"__call", LuaResponseFunctionCall},
  {"__index", LuaResponseFunctionTableIndex},
  {"__eq",   LuaResponseFunctionEq},
  {NULL, NULL}
};

