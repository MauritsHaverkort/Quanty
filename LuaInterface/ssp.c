static int LuaHelp(lua_State *L)
{
  uint32_t len = (uint32_t)lua_gettop(L);
  if(len>1) return luaL_error(L, "Error in .Help(arg) expected 0 (general help) or 1 (help on function) arguments\n");
  if(len==0)
  {
    printf("-------------------------------------------------------------------------\n");
    printf("-          solid state physics model library  - lua interface           -\n");
    printf("-          (c) Maurits W. Haverkort                                     -\n");
    printf("-          http://www.solstatphys.org                                   -\n");
    printf("-------------------------------------------------------------------------\n");
    printf("\n\n");
    printf("Available functions\n");
    printf("\n");
    printf("ssp.Help()\n");
    printf("  - shows this file\n");
    printf("  - excepts as input one of the other functions in order to show help on that function\n");
    printf("  - example: ssp.Help(\"ssp.NewOperator\");\n");
    printf("\n");
    printf("ssp.NewOperator()\n");
    printf("  - creates a new operator\n");
    printf("  - expects as input on the first position a string describing the type of operator\n");
    printf("  - other parameters depend on the type of operator defined in the first string\n");
    printf("  - example: ssp.NewOperator(\"Sz\", NFermion, IndexUp, IndexDn);\n");
    printf("\n");
    printf("ssp.NewWavefunction()\n");
    printf("  - creates a new wavefunction\n");
    printf("  - expects as input the number of fermions, number of bosons, a table of determinants with the prefactor\n");
    printf("  - example: ssp.NewWavefunction(10,0,{{{0x80,0x00},1},{{0x04,0x00},1}});\n");
    printf("\n");
    printf("ssp.Copy() or ssp.Clone()\n");
    printf("  - creats a copy of a wavefunction or operator\n");
    printf("  - expects as input an object, returns a copy of that object\n");
    printf("  - example: ssp.Copy(psi);\n");
    printf("\n");
    printf("ssp.Conjugate() or ssp.C()\n");
    printf("  - creats a copy of a wavefunction or operator and conjugates this copy\n");
    printf("  - expects as input an object, returns a copy of that object\n");
    printf("  - example: ssp.Conjugate(psi);\n");
    printf("\n");
    printf("ssp.ConjugateTranspose() or ssp.CT()\n");
    printf("  - creats a copy of a wavefunction or operator and conjugates transposes this copy\n");
    printf("  - expects as input an object, returns a copy of that object\n");
    printf("  - example: ssp.ConjugateTranspose(psi);\n");
    printf("\n");
    printf("ssp.Transpose() or ssp.T()\n");
    printf("  - creats a copy of a wavefunction or operator and transposes this copy\n");
    printf("  - expects as input an object, returns a copy of that object\n");
    printf("  - example: ssp.Transpose(psi);\n");
    printf("\n");
    printf("ssp.EvolveToGroundstate()\n");
    printf("  - calculates eigenfunctions of the Hamiltonian, starting from the specified state\n");
    printf("  - expects as input an operator and a (list of) wavefunction(s), a list of options on the last position is allowed\n");
    printf("  - example ssp.EvolveToGroundstate(OppSy,{psi1,psi2});\n");
    printf("\n");
    printf("ssp.CreateSpectra()\n");
    printf("  - calculates spectra of a giving Hamiltonian, starting from the specified state and using a specific transition operator\n");
    printf("  - expects as input an operator a (list of) transition operator(s) and a (list of) wavefunction(s), a list of options on the last position is allowed\n");
    printf("  - example ssp.CreateSpectra(OppH,{T1,T2},{psi1,psi2});\n");
    printf("\n");
    
    printf("\n\n");
    printf("Available objects\n");
    printf("\n");
    printf("operator\n");
    printf("  - stores a quantum mechanical operator based on second quantization\n");
    printf("\n");
    printf("wavefunction\n");
    printf("  - stores a wavefunction\n");
    printf("\n");
    printf("complex\n");
    printf("  - stores a complex number\n");
    printf("\n");
    printf("spectra\n");
    printf("  - stores a spectral function\n");
    printf("\n");
    
  }
  if(len==1)
  {
    char const *Help =lua_tostring(L, 1);
    if (strcmp(Help, "ssp.NewOperator") == 0)
    {
      printf("ssp.NewOperator()\n");
      printf("  - creates a new operator\n");
      printf("  - expects as input on the first position a string describing the type of operator\n");
      printf("  - other parameters depend on the type of operator defined in the first string\n");
      printf("\n");
      printf("\n");
      printf("\n");
    }
    if (strcmp(Help, "ssp.NewWavefunction") == 0)
    {
      printf("ssp.NewWavefunction()\n");
      printf("  - creates a new wavefunction\n");
      printf("  - expects as input the number of fermions, number of bosons, a table of determinants with the prefactor\n");
      printf("\n");
      printf("\n");
      printf("\n");
    }
    if ((strcmp(Help, "ssp.Copy") == 0)||(strcmp(Help, "ssp.Clone") == 0))
    {
      printf("ssp.Copy() or ssp.Clone()\n");
      printf("  - creats a copy of a wavefunction or operator\n");
      printf("  - expects as input an object, returns a copy of that object\n");
      printf("\n");
      printf("\n");
      printf("\n");
    }
    if ((strcmp(Help, "ssp.Conjugate") == 0)||(strcmp(Help, "ssp.C") == 0))
    {
      printf("ssp.Conjugate() or ssp.C()\n");
      printf("  - creats a copy of a wavefunction or operator and conjugates this copy\n");
      printf("  - expects as input an object, returns a copy of that object\n");
      printf("\n");
      printf("\n");
      printf("\n");
    }
    if ((strcmp(Help, "ssp.ConjugateTranspose") == 0)||(strcmp(Help, "ssp.CT")==0))
    {
      printf("ssp.ConjugateTranspose() or ssp.CT()\n");
      printf("  - creats a copy of a wavefunction or operator and conjugates transposes this copy\n");
      printf("  - expects as input an object, returns a copy of that object\n");
      printf("\n");
      printf("\n");
      printf("\n");
    }
    if ((strcmp(Help, "ssp.Transpose") == 0)||(strcmp(Help, "ssp.T") == 0))
    {
      printf("ssp.Transpose() or ssp.T()\n");
      printf("  - creats a copy of a wavefunction or operator and transposes this copy\n");
      printf("  - expects as input an object, returns a copy of that object\n");
      printf("\n");
      printf("\n");
      printf("\n");
    }
    if (strcmp(Help, "ssp.EvolveToGroundstate") == 0)
    {
      printf("ssp.EvolveToGroundstate()\n");
      printf("  - calculates eigenfunctions of the Hamiltonian, starting from the specified state\n");
      printf("  - expects as input an operator and a (list of) wavefunction(s), a list of options on the last position is allowed\n");
      printf("\n");
      printf("\n");
      printf("\n");
    }
    if (strcmp(Help, "ssp.CreateSpectra") == 0)
    {
      printf("ssp.CreateSpectra()\n");
      printf("  - calculates spectra of a giving Hamiltonian, starting from the specified state and using a specific transition operator\n");
      printf("  - expects as input an operator a (list of) transition operator(s) and a (list of) wavefunction(s), a list of options on the last position is allowed\n");
      printf("\n");
      printf("\n");
      printf("\n");
    }
    
    
    
    
    if (strcmp(Help, "operator") == 0)
    {
      printf("operator\n");
      printf("  - stores a quantum mechanical operator based on second quantization\n");
      printf("\n");
      printf("\n");
      printf("\n");
    }
    if (strcmp(Help, "wavefunction") == 0)
    {
      printf("wavefunction\n");
      printf("  - stores a wavefunction\n");
      printf("\n");
      printf("\n");
      printf("\n");
    }
    if (strcmp(Help, "complex") == 0)
    {
      printf("complex\n");
      printf("  - stores a complex number\n");
      printf("\n");
      printf("\n");
      printf("\n");
    }
    if (strcmp(Help, "spectra") == 0)
    {
      printf("spectra\n");
      printf("  - stores a spectral function\n");
      printf("\n");
      printf("\n");
      printf("\n");
    }
    
    
    if (strcmp(Help, "") == 0)
    {
      printf("\n");
      printf("\n");
      printf("\n");
    }
  }
  return 0;
}



static int LuaNewOperator(lua_State * L)
{
  if(lua_gettop(L)==0)
  {
    lua_pushfstring (L, "ERROR: LuaNewOperator called with 0 arguments, expect at least one\n");
    lua_error(L);
    return 1;
  }
  if(lua_type(L, 1)!=LUA_TSTRING)
  {
    lua_pushfstring (L, "ERROR: first argument in LuaNewOperator should be a string (use \"Oppname\" as first argument)\n");
    lua_error(L);
    return 1;
  }
  char const *OT =lua_tostring(L, 1);
  if ((strcmp(OT, "Coulomb") == 0)||(strcmp(OT, "U") == 0))
  {
    lua_remove(L, 1);
    return LuaOperatorNewCoulomb(L);
  }
  if (strcmp(OT, "Sx") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewSx(L);
  }
  if (strcmp(OT, "Sy") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewSy(L);
  }
  if (strcmp(OT, "Sz") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewSz(L);
  }
  if (strcmp(OT, "Splus") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewSplus(L);
  }
  if (strcmp(OT, "Smin") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewSmin(L);
  }
  if ((strcmp(OT, "Ssqr") == 0)||(strcmp(OT, "S^2") == 0)||(strcmp(OT, "S2") == 0))
  {
    lua_remove(L, 1);
    return LuaOperatorNewSsqr(L);
  }
  if (strcmp(OT, "Lx") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewLx(L);
  }
  if (strcmp(OT, "Ly") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewLy(L);
  }
  if (strcmp(OT, "Lz") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewLz(L);
  }
  if (strcmp(OT, "Lplus") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewLplus(L);
  }
  if (strcmp(OT, "Lmin") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewLmin(L);
  }
  if ((strcmp(OT, "Lsqr") == 0)||(strcmp(OT, "L^2") == 0)||(strcmp(OT, "L2") == 0))
  {
    lua_remove(L, 1);
    return LuaOperatorNewLsqr(L);
  }
  if (strcmp(OT, "Jx") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewJx(L);
  }
  if (strcmp(OT, "Jy") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewJy(L);
  }
  if (strcmp(OT, "Jz") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewJz(L);
  }
  if (strcmp(OT, "Jplus") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewJplus(L);
  }
  if (strcmp(OT, "Jmin") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewJmin(L);
  }
  if ((strcmp(OT, "Jsqr") == 0)||(strcmp(OT, "J^2") == 0)||(strcmp(OT, "J2") == 0))
  {
    lua_remove(L, 1);
    return LuaOperatorNewJsqr(L);
  }
  if (strcmp(OT, "ldots") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewldots(L);
  }
  if (strcmp(OT, "CF") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewCF(L);
  }
  if ((strcmp(OT, "Create") == 0)||(strcmp(OT, "Cr") == 0))
  {
    lua_remove(L, 1);
    return LuaOperatorNewCreate(L);
  }
  if ((strcmp(OT, "Annihilate") == 0)||(strcmp(OT, "An") == 0))
  {
    lua_remove(L, 1);
    return LuaOperatorNewAnnihilate(L);
  }
  if ((strcmp(OT, "Number") == 0)||(strcmp(OT, "CrAn") == 0)||(strcmp(OT, "n") == 0))
  {
    lua_remove(L, 1);
    return LuaOperatorNewNumber(L);
  }
  return luaL_error(L, "ERROR: first argument in ssp.NewOperator (%s) is not recognized as a known operatortype\n See ssp.Help(\"ssp.NewOperator\") for more information\n",OT);
}

/*  on stack: string describing the symmetry
 *  integer describing the angular momentum
 *  table with eigen energies of type double
 *
 *  output table of Harmonics Expansion
 */
static int LuaPotentialExpandedOnYlm(lua_State *L)
{
  char *s;
  double * eigenvalues;
  HarmonicsExpansion Akm;
  uint32_t l,len;
  s=(char*)luaL_checkstring(L, 1);
  l=luaL_checkint(L,2);
  eigenvalues=luaL_checkarray_double(L, 3, &len);
  if(PotentialExpandedOnY(s, l, eigenvalues, &Akm))
  {
    luaL_error(L, "Potential of \"%s\" symmetry for angular momentum %d not found\n",s,l);
  }
  LuaPushTable_HarmonicsExpansion(L, &Akm);
  FreeHarmonicsExpansion(&Akm);free(eigenvalues);
  return 1;
}


static int LuaNewWaveFunction(lua_State * L)
{
  char *detchar;
  Complex16 temp;
  uint32_t QComplex=0,NDet,i,detleng;
  uint8_t *det;
	WaveFunctionType *psi;
  uint32_t len = (uint32_t)lua_gettop(L);
	psi = (WaveFunctionType*) lua_newuserdata(L, sizeof(WaveFunctionType));
  luaL_setmetatable(L, "ssp.wavefunction");
  psi->NFermion=(uint32_t)luaL_checkunsigned(L, 1);
  psi->NBoson  =(uint32_t)luaL_checkunsigned(L, 2);
  psi->QComplex=1;
  if(len==3)
  {
    psi->NBitsKey=16;
  }
  else if(len==4)
  {
    psi->NBitsKey  =(uint32_t)luaL_checkunsigned(L, 4);
  }
  else
  {
    lua_pushfstring (L, "Error in LuaWavefunctionNew, expected 3 or 4 arguments, not %d\n",len);
    lua_error(L);
    return 1;
  }
  WaveFunctionInit(psi);
  luaL_checktype(L, 3, LUA_TTABLE);
  NDet = (uint32_t)lua_objlen(L, 3);
  for(i=1;i<=NDet;i++)
  {
    lua_rawgeti(L, 3, i);
    luaL_checktype(L, -1, LUA_TTABLE);
    if(lua_objlen(L, -1)!=2)
    {
      lua_pushfstring (L, "Error in LuaWavefunctionNew, determinant number %d is not of the form {Determinant,prefactor}\n",i);
      lua_error(L);
      return 1;
    }
    else
    {
      lua_rawgeti(L, -1, 1);
      if(lua_isstring(L, -1))
      {
        detchar = (char *) luaL_checkstring(L, -1);
        det = StringToDet(psi->NFermion, psi->NBoson, detchar);
        detleng = psi->NBoson + Min(psi->NFermion, 1 + (psi->NFermion - 1) / 8);
      }
      else
      {
        det=luaL_checkarray_unit8_t(L, -1, &detleng);
      }
      if(detleng!=psi->NWord)
      {
        lua_pushfstring (L, "Error in LuaWavefunctionNew, determinant number %d should be a table of %d words\n",i,psi->NWord);
        lua_error(L);
        return 1;
      }
      lua_pop(L, 1);
      lua_rawgeti(L, -1, 2);
      if(lua_isnumber(L, -1))
      {
        temp.Re=(double)luaL_checknumber(L, -1);
        temp.Im=0;
      }
      else
      {
        QComplex=1;
        temp=*(Complex16 *)luaL_checkudata(L, -1, "Lua.complex");
      }
      lua_pop(L, 1);
      CheckDetPhysical(det, psi->NFermion, psi->NBoson);
      ComplexWaveFunctionAddElement(psi, det, temp.Re, temp.Im);
      free(det);
    }
    lua_pop(L, 1);
  }
  if(QComplex==0) MakeRealWaveFunctionFromComplex(*psi, psi);
  return 1;
}

/*  copy or clone the last element on the lua stack and place it on the stack
 *  works recursively on tables
 *  allocates new memory for user types
 */
static int Luacopy(lua_State * L)
{
  uint32_t len = (uint32_t)lua_gettop(L);
  switch (lua_type(L,len))
  {
    case LUA_TNIL:
      lua_pushnil(L);
      break;
    case LUA_TBOOLEAN:
      lua_pushboolean(L,lua_toboolean(L,len));
      break;
    case LUA_TNUMBER:
      lua_pushnumber(L,lua_tonumber(L,len));
      break;
    case LUA_TSTRING:
      lua_pushlstring(L,lua_tostring(L,len),lua_strlen(L,len));
      break;
    case LUA_TLIGHTUSERDATA:
      lua_pushlightuserdata(L,(void*)lua_touserdata(L,len));
      break;
    case
      LUA_TUSERDATA:
      if (luaL_testudata(L, len, "Lua.complex") != NULL)
      {
        Complex16 *a = (Complex16 *) luaL_checkudata(L, len, "Lua.complex");
        Complex16 *res;
        res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
        res->Re = a->Re;
        res->Im = a->Im;
        luaL_setmetatable(L, "Lua.complex");
        return 1;
      }
      if (luaL_testudata(L, len, "ssp.wavefunction") != NULL)
      {
        WaveFunctionType *psi = (WaveFunctionType*) luaL_checkudata(L, len, "ssp.wavefunction");
        WaveFunctionType *psiNew;
        psiNew = (WaveFunctionType *) lua_newuserdata(L, sizeof(WaveFunctionType));
        WaveFunctionCopy(*psi, psiNew);
        luaL_setmetatable(L, "ssp.wavefunction");
        return 1;
      }
      if (luaL_testudata(L, len, "ssp.operator") != NULL)
      {
        OperatorType *O = (OperatorType*) luaL_checkudata(L, len, "ssp.operator");
        OperatorType *ONew;
        ONew = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
        OperatorCopy(ONew, *O);
        luaL_setmetatable(L, "ssp.operator");
        return 1;
      }
      lua_pushfstring (L, "Error: copying or cloning an object of unrecognised type\n");
      lua_error(L);
      break;
    case LUA_TTABLE:
      lua_newtable(L);
      lua_pushnil(L);
      while(lua_next(L, len) != 0)
      {
        Luacopy(L);
        lua_pushvalue(L, -3);
        lua_insert(L, -2);
        lua_settable(L, -5);
        lua_pop(L, 1);
      }
      break;
    default:
      lua_pushfstring (L, "Error: copying or cloning an object of unrecognised type\n");
      lua_error(L);
      break;
  }
  return 1;
}

/*  copy or clone the last element on the lua stack and place it on the stack
 *  works recursively on tables
 *  allocates new memory for user types
 *  complex conjugates user types that can be conjugated
 */
static int Luaconjugate(lua_State * L)
{
  uint32_t len = (uint32_t)lua_gettop(L);
  switch (lua_type(L,len))
  {
    case LUA_TNIL:
      lua_pushnil(L);
      break;
    case LUA_TBOOLEAN:
      lua_pushboolean(L,lua_toboolean(L,len));
      break;
    case LUA_TNUMBER:
      lua_pushnumber(L,lua_tonumber(L,len));
      break;
    case LUA_TSTRING:
      lua_pushlstring(L,lua_tostring(L,len),lua_strlen(L,len));
      break;
    case LUA_TLIGHTUSERDATA:
      lua_pushlightuserdata(L,(void*)lua_touserdata(L,len));
      break;
    case
    LUA_TUSERDATA:
      if (luaL_testudata(L, len, "Lua.complex") != NULL)
      {
        Complex16 *a = (Complex16 *) luaL_checkudata(L, len, "Lua.complex");
        Complex16 *res;
        res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
        res->Re = a->Re;
        res->Im = -a->Im;
        luaL_setmetatable(L, "Lua.complex");
        return 1;
      }
      if (luaL_testudata(L, len, "ssp.wavefunction") != NULL)
      {
        WaveFunctionType *psi = (WaveFunctionType*) luaL_checkudata(L, len, "ssp.wavefunction");
        WaveFunctionType *psiNew;
        psiNew = (WaveFunctionType *) lua_newuserdata(L, sizeof(WaveFunctionType));
        WaveFunctionCopy(*psi, psiNew);
        WaveFunctionConjugate(*psiNew);
        luaL_setmetatable(L, "ssp.wavefunction");
        return 1;
      }
      if (luaL_testudata(L, len, "ssp.operator") != NULL)
      {
        OperatorType *O = (OperatorType*) luaL_checkudata(L, len, "ssp.operator");
        OperatorType *ONew;
        ONew = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
        OperatorCopy(ONew, *O);
        OperatorComplexConjugate(ONew);
        luaL_setmetatable(L, "ssp.operator");
        return 1;
      }
      lua_pushfstring (L, "Error: copying or cloning an object of unrecognised type\n");
      lua_error(L);
      break;
    case LUA_TTABLE:
      lua_newtable(L);
      lua_pushnil(L);
      while(lua_next(L, len) != 0)
      {
        Luacopy(L);
        lua_pushvalue(L, -3);
        lua_insert(L, -2);
        lua_settable(L, -5);
        lua_pop(L, 1);
      }
      break;
    default:
      lua_pushfstring (L, "Error: copying or cloning an object of unrecognised type\n");
      lua_error(L);
      break;
  }
  return 1;
}

static int Luaconjugatetranspose(lua_State * L)
{
  uint32_t len = (uint32_t)lua_gettop(L);
  switch (lua_type(L,len))
  {
    case LUA_TNIL:
      lua_pushnil(L);
      break;
    case LUA_TBOOLEAN:
      lua_pushboolean(L,lua_toboolean(L,len));
      break;
    case LUA_TNUMBER:
      lua_pushnumber(L,lua_tonumber(L,len));
      break;
    case LUA_TSTRING:
      lua_pushlstring(L,lua_tostring(L,len),lua_strlen(L,len));
      break;
    case LUA_TLIGHTUSERDATA:
      lua_pushlightuserdata(L,(void*)lua_touserdata(L,len));
      break;
    case
    LUA_TUSERDATA:
      if (luaL_testudata(L, len, "Lua.complex") != NULL)
      {
        Complex16 *a = (Complex16 *) luaL_checkudata(L, len, "Lua.complex");
        Complex16 *res;
        res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
        res->Re = a->Re;
        res->Im = -a->Im;
        luaL_setmetatable(L, "Lua.complex");
        return 1;
      }
      if (luaL_testudata(L, len, "ssp.wavefunction") != NULL)
      {
        WaveFunctionType *psi = (WaveFunctionType*) luaL_checkudata(L, len, "ssp.wavefunction");
        WaveFunctionType *psiNew;
        psiNew = (WaveFunctionType *) lua_newuserdata(L, sizeof(WaveFunctionType));
        WaveFunctionCopy(*psi, psiNew);
        WaveFunctionConjugate(*psiNew);
        luaL_setmetatable(L, "ssp.wavefunction");
        return 1;
      }
      if (luaL_testudata(L, len, "ssp.operator") != NULL)
      {
        OperatorType *O = (OperatorType*) luaL_checkudata(L, len, "ssp.operator");
        OperatorType *ONew;
        ONew = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
        OperatorConjugateTranspose(*O, ONew);
        luaL_setmetatable(L, "ssp.operator");
        return 1;
      }
      lua_pushfstring (L, "Error: copying or cloning an object of unrecognised type\n");
      lua_error(L);
      break;
    case LUA_TTABLE:
      lua_newtable(L);
      lua_pushnil(L);
      while(lua_next(L, len) != 0)
      {
        Luacopy(L);
        lua_pushvalue(L, -3);
        lua_insert(L, -2);
        lua_settable(L, -5);
        lua_pop(L, 1);
      }
      break;
    default:
      lua_pushfstring (L, "Error: copying or cloning an object of unrecognised type\n");
      lua_error(L);
      break;
  }
  return 1;
}

static int Luatranspose(lua_State * L)
{
  uint32_t len = (uint32_t)lua_gettop(L);
  switch (lua_type(L,len))
  {
    case LUA_TNIL:
      lua_pushnil(L);
      break;
    case LUA_TBOOLEAN:
      lua_pushboolean(L,lua_toboolean(L,len));
      break;
    case LUA_TNUMBER:
      lua_pushnumber(L,lua_tonumber(L,len));
      break;
    case LUA_TSTRING:
      lua_pushlstring(L,lua_tostring(L,len),lua_strlen(L,len));
      break;
    case LUA_TLIGHTUSERDATA:
      lua_pushlightuserdata(L,(void*)lua_touserdata(L,len));
      break;
    case
    LUA_TUSERDATA:
      if (luaL_testudata(L, len, "Lua.complex") != NULL)
      {
        Complex16 *a = (Complex16 *) luaL_checkudata(L, len, "Lua.complex");
        Complex16 *res;
        res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
        res->Re = a->Re;
        res->Im = a->Im;
        luaL_setmetatable(L, "Lua.complex");
        return 1;
      }
      if (luaL_testudata(L, len, "ssp.wavefunction") != NULL)
      {
        WaveFunctionType *psi = (WaveFunctionType*) luaL_checkudata(L, len, "ssp.wavefunction");
        WaveFunctionType *psiNew;
        psiNew = (WaveFunctionType *) lua_newuserdata(L, sizeof(WaveFunctionType));
        WaveFunctionCopy(*psi, psiNew);
        luaL_setmetatable(L, "ssp.wavefunction");
        return 1;
      }
      if (luaL_testudata(L, len, "ssp.operator") != NULL)
      {
        OperatorType *O = (OperatorType*) luaL_checkudata(L, len, "ssp.operator");
        OperatorType *ONew;
        ONew = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
        OperatorTranspose(*O, ONew);
        luaL_setmetatable(L, "ssp.operator");
        return 1;
      }
      lua_pushfstring (L, "Error: copying or cloning an object of unrecognised type\n");
      lua_error(L);
      break;
    case LUA_TTABLE:
      lua_newtable(L);
      lua_pushnil(L);
      while(lua_next(L, len) != 0)
      {
        Luacopy(L);
        lua_pushvalue(L, -3);
        lua_insert(L, -2);
        lua_settable(L, -5);
        lua_pop(L, 1);
      }
      break;
    default:
      return luaL_error(L, "Error: copying or cloning an object of unrecognised type\n");
      break;
  }
  return 1;
}


static int Luarotate(lua_State *L)
{
  struct
  {
    double zero;
  } Options;
  Options.zero=100*DBL_EPSILON;
  MatrixType *U;
  uint32_t len = (uint32_t)lua_gettop(L);
  if(len<2)
  {
    return luaL_error(L, "ssp.Rotate called with %d arguments, 2 expected\n",len);
  }
  if (luaL_testudata(L, len-1, "ssp.wavefunction") != NULL)
  {
    
    return 0;
  }
  if (luaL_testudata(L, len-1, "ssp.operator") != NULL)
  {
    U = luaL_checkMatrix(L, len);
    OperatorType *O = (OperatorType*) luaL_checkudata(L, len-1, "ssp.operator");
    if(U->Ni!=U->Nj) luaL_error(L, "Error while rotating operator. Rotation matrix should be square (Ni=%d Nj=%d)",U->Ni,U->Nj);
    if(U->Ni!=O->NFermion+O->NBoson) luaL_error(L, "Error while rotating operator. Rotation matrix should have dimension (%d) of number of fermions (%d) + bosons (%d) of operator",U->Ni,O->NFermion,O->NBoson);
    OperatorType *ONew;
    ONew = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    RotateOperatorNormalOrder(ONew, *O, *U,&QDetNotOrderedFcFaBaBc);
    luaL_setmetatable(L, "ssp.operator");
    OperatorChop(ONew, Options.zero);
    return 1;
  }
  if(lua_type(L,len-1)==LUA_TTABLE)
  {
    
    
    return 0;
  }
  
  
  
  return luaL_error(L, "ssp.Rotate called unrecognized type. Can only rotate wavefunctions operators and tables of these\n");
}



static int LuaCreateSpectra(lua_State *L)
{
  struct
  {
    uint32_t NTri;
    double epsilon;
    Bool restrictions;
    DeterminantRestrictionType DeterminantRestrictions;
    double Emin;
    double Emax;
    uint32_t NE;
    double Gamma;
    uint32_t Npsi;
    uint32_t NT;
  } Options;
  Options.NTri=200;
  Options.epsilon=10000*DBL_EPSILON;
  Options.restrictions=False;
  Options.Emin=DBL_MAX;
  Options.Emax=-DBL_MAX;
  Options.NE=1000;
  Options.Gamma=0;
  SpectraTypeArray *Spec;
  TriDiagonalMatrixType *GTri;
  OperatorType *H,*T;
  WaveFunctionType *psi,Tpsi;
  uint32_t i,iT,ipsi,len;
  double Eminpole,Emaxpole, *val, *psiEn, psiEntemp;
  
  /* check number of input parameters */
  len = (uint32_t)lua_gettop(L);
  if(!((len==3)||(len==4)))
  {
    return luaL_error(L, "Error CreateSpectra called with %d arguments, 3 or 4 expected\n",len);
  }
  
  /* read operator and wavefunctions */
  H = (OperatorType *) luaL_checkudata(L, 1, "ssp.operator");
  if(lua_type(L,2)==LUA_TTABLE)
  {
    T=luaL_checkarray_OperatorType(L, 2, &Options.NT);
  }
  else
  {
    Options.NT = 1;
    T = (OperatorType *) luaL_checkudata(L, 2, "ssp.operator");
  }
  if(lua_type(L,3)==LUA_TTABLE)
  {
    psi=luaL_checkarray_WaveFunctionType(L, 3, &Options.Npsi);
  }
  else
  {
    Options.Npsi = 1;
    psi = (WaveFunctionType *) luaL_checkudata(L, 3, "ssp.wavefunction");
  }
  /* read options */
  
  
  
  /* calculate spectra in Tridiagonal form */
  GTri = malloc(Options.NT*Options.Npsi*sizeof(TriDiagonalMatrixType));
  for(i=0;i<Options.NT*Options.Npsi;i++)
  {
    GTri[i].N=Options.NTri;
    InitTriDiagonalMatrix(&GTri[i]);
  }
  /* calculate groundstate energy */
  psiEn=malloc(Options.Npsi*Options.NT*sizeof(double));
  for(ipsi=0;ipsi<Options.Npsi;ipsi++)
  {
    PsiHermitianOperatorPsi(*H, psi[ipsi], &psiEntemp);
    for(iT=0;iT<Options.NT;iT++)
    {
      psiEn[iT*Options.Npsi+ipsi]=psiEntemp;
    }
  }
  /* loop over all Transition operators */
  for(iT=0;iT<Options.NT;iT++)
  {
    /* loop over all wavefunctions */
    for(ipsi=0;ipsi<Options.Npsi;ipsi++)
    {
      OperatorPsi(T[iT], psi[ipsi], &Tpsi);
      LanczosTriDiagonalize(*H, Tpsi, Options.epsilon, &GTri[iT*Options.Npsi+ipsi], Max(Options.Gamma,Options.epsilon));
    }
  }
  /* if Emin=DBL_MAX or Emax=DBL_MAX or Gamma==0 determine these */
  Eminpole=Options.Emin;
  Emaxpole=Options.Emax;
  if((Options.Emin==DBL_MAX)||(Options.Emax==-DBL_MAX))
  {
    /* get eigenvalues of all Tri diagonal matrices and set smallest and largest eigenvalue */
    /* range is 20% larger */
    val=malloc((Options.NTri)*sizeof(double));
    for(i=0;i<Options.NT*Options.Npsi;i++)
    {
      TriDiagonalMatrixToEigenvalues(GTri[i], val);
      if((val[0]-psiEn[i])<Eminpole)Eminpole=val[0]-psiEn[i];
      if((val[Options.NTri-1]-psiEn[i])>Emaxpole)Emaxpole=val[Options.NTri-1]-psiEn[i];
    }
    if(Options.Emin== DBL_MAX) Options.Emin = Eminpole - 0.25 * (Emaxpole-Eminpole);
    if(Options.Emax==-DBL_MAX) Options.Emax = Emaxpole + 0.25 * (Emaxpole-Eminpole);
  }
  if(Options.Gamma==0) Options.Gamma = 10*(Options.Emax-Options.Emin)/((double)Options.NE);
  /* transform tridiagonal form to spectral form */
  Spec = (SpectraTypeArray *) lua_newuserdata(L, sizeof(SpectraTypeArray));
  luaL_setmetatable(L, "ssp.spectra");
  Spec->NSpectra=Options.NT*Options.Npsi;
  Spec->Spec=malloc(Spec->NSpectra*sizeof(SpectraType));
  Spec->Spec[0].N=Options.NE;
  Spec->Spec[0].Emin=Options.Emin;
  Spec->Spec[0].Emax=Options.Emax;
  Spec->Spec[0].dE=(Options.Emax-Options.Emin)/((double)Options.NE);
  Spec->Spec[0].EminPole=Options.Emin;
  Spec->Spec[0].EmaxPole=Options.Emax;
  Spec->Spec[0].Gamma=Options.Gamma;
  for(i=0;i<Spec->NSpectra;i++)
  {
    SpectrumAllocate(Spec->Spec[0], &(Spec->Spec[i]));
    TriDiagonalMatrixToSpectrum(psiEn[i], GTri[i], &(Spec->Spec[i]));
  }
  for(i=0;i<Spec->NSpectra;i++)
  {
    FreeTriDiagonalMatrix(&GTri[i]);
  }
  free(GTri);
  free(psiEn);
  return 1;
}


/* given one operator and a (table of) wavefunctions calculate the ground-state (and excited states) starting from those wavefunctions.
 *
 */
static int LuaEvolveToGroundstate(lua_State * L)
{
  struct
  {
    uint32_t Npsi;
    double zero;
    double epsilon;
    Bool restrictions;
    DeterminantRestrictionType DeterminantRestrictions;
  } Options;
  Options.zero=100*DBL_EPSILON;
  Options.epsilon=10000*DBL_EPSILON;
  Options.restrictions=False;
  OperatorType *Opp;
  WaveFunctionType *psi,*temp;
  uint32_t i,len;
  /* check number of input parameters */
  len = (uint32_t)lua_gettop(L);
  if(!((len==2)||(len==3)))
  {
    return luaL_error(L, "Error EvolveToGroundstate called with %d arguments, 2 or 3 expected\n",len);
  }
  /* read operator and wavefunctions */
  Opp = (OperatorType *) luaL_checkudata(L, 1, "ssp.operator");
  if(lua_type(L,2)==LUA_TTABLE)
  {
    psi=luaL_checkarray_WaveFunctionType(L, 2, &Options.Npsi);
  }
  else
  {
    Options.Npsi = 1;
    psi = (WaveFunctionType *) luaL_checkudata(L, 2, "ssp.wavefunction");
  }
  /* read options */
  
  
  
  /* calculate groundstate */
  
  for (i = 0; i < Options.Npsi; i++)
  {
    printf("start of state #%d/%d",i+1,Options.Npsi);
    LanczosGroundState(*Opp, &(psi[i]), &Options.epsilon, Options.zero, 0, psi, i);
  }
  
  
  /* copy changed wavefunctions back to lua stack */
  if(lua_type(L,2)==LUA_TTABLE)
  {
    for(i=0;i<Options.Npsi;i++)
    {
      lua_rawgeti(L, 2, i+1);
      temp = (WaveFunctionType *)luaL_checkudata(L, -1, "ssp.wavefunction");
      *temp = psi[i];
      lua_pop(L, 1);
    }

  }
  return 0;
}


static const struct luaL_Reg ssp_functions[] = {
  {"Help", LuaHelp},
  {"NewOperator", LuaNewOperator},
  {"PotentialExpandedOnYlm", LuaPotentialExpandedOnYlm},
  {"NewWavefunction", LuaNewWaveFunction},
  {"Copy", Luacopy},
  {"Clone", Luacopy},
  {"Conjugate", Luaconjugate},
  {"C", Luaconjugate},
  {"ConjugateTranspose", Luaconjugatetranspose},
  {"CT", Luaconjugatetranspose},
  {"Transpose", Luatranspose},
  {"T", Luatranspose},
  {"Rotate", Luarotate},
  {"EvolveToGroundstate", LuaEvolveToGroundstate},
  {"CreateSpectra", LuaCreateSpectra},
  {"BraKet", LuaBraKet},           /* function in file BraKet.c */
  {"Eval", LuaEval},               /* function in file eval.c */
  {NULL, NULL}
};

void luaconstants(lua_State *L)
{
  Complex16 *luaI;
  luaI = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
  luaI->Re = 0.;
  luaI->Im = 1.;
  luaL_setmetatable(L, "Lua.complex");
  lua_setglobal(L, "I");
  lua_pushnumber(L, 3.14159265358979323846264338328);
  lua_setglobal(L, "Pi");
  lua_pushnumber(L, 2.71828182845904523536028747135);
  lua_setglobal(L, "E");
}

int luaopen_ssp(lua_State * L)
{
  luaL_newmetatable(L, "ssp.operator");
  luaL_setfuncs(L, Operator_meta, 0);
  
  luaL_newmetatable(L, "ssp.wavefunction");
  luaL_setfuncs(L, Wavefunction_meta, 0);
  
  luaL_newmetatable(L, "ssp.spectra");
  luaL_setfuncs(L, Spectra_meta, 0);
  
  luaL_newlib(L, ssp_functions);
  lua_setglobal(L, "ssp");
  
  
  return 1;
}