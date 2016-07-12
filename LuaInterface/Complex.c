static int LuaComplexNew(lua_State * L)
{
  Complex16 *res;
  double re = luaL_checknumber(L, 1);
  double im = luaL_checknumber(L, 2);
  res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
  res->Re = re;
  res->Im = im;
  luaL_setmetatable(L, "Complex_Type");
  return 1;
}

static int LuaComplexI(lua_State * L)
{
  Complex16 *res;
  res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
  res->Re = 0.;
  res->Im = 1.;
  luaL_setmetatable(L, "Complex_Type");
  return 1;
}

static int LuaComplexIndex(lua_State * L)
{
  Complex16 *a = (Complex16 *) luaL_checkudata(L, 1, "Complex_Type");
  const char *str = lua_tostring(L, 2);
  
  switch (str[0])
  {
    case 'r':
    case 'R':
      lua_pushnumber(L, a->Re);
      return 1;
      
    case 'i':
    case 'I':
      lua_pushnumber(L, a->Im);
      return 1;
  }
  return luaL_error(L, "Error Complex variable not found");;
}

static int LuaComplexNewIndex(lua_State * L)
{
  Complex16 *a = (Complex16 *) luaL_checkudata(L, 1, "Complex_Type");
  const char *str = lua_tostring(L, 2);
  double v = lua_tonumber(L, 3);
  
  switch (str[0])
  {
    case 'r':
    case 'R':
      a->Re = v;
      return 0;
      
    case 'i':
    case 'I':
      a->Im = v;
      return 0;
  }
  return luaL_error(L, "Error index %s not found in Complex Object\n",str);
}

static int LuaComplexAdd(lua_State * L)
{
  uint32_t len = (uint32_t)lua_gettop(L);
  switch (lua_type(L,1))
  {
    case LUA_TNIL:
    {
      // (nil) operator (original variable)
      return luaL_error(L, "Error: nil plus complex number\n");
    }
    case LUA_TBOOLEAN:
    {
      // (boolean) operator (original variable)
      return luaL_error(L, "Error: boolean plus complex number\n");
    }
    case LUA_TLIGHTUSERDATA:
    {
      // (ligth user data) operator (original variable)
      return luaL_error(L, "Error: lightuserdata plus complex number\n");
    }
    case LUA_TNUMBER:
    {
      // (number) operator (original variable)
      double a;
      Complex16 *b, *res;
      a = luaL_checknumber(L, 1);
      b = (Complex16 *) luaL_checkudata(L, len, "Complex_Type");
      res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
      res->Re =  a + b->Re;
      res->Im =      b->Im;
      luaL_setmetatable(L, "Complex_Type");
      return 1;
    }
    case LUA_TSTRING:
    {
      // (string) operator (original variable)
      return luaL_error(L, "Error: character string plus complex number\n");
    }
    case LUA_TTABLE:
    {
      lua_newtable(L);
      lua_pushnil(L);
      while(lua_next(L, 1) != 0)
      {
        lua_insert(L, 1);
        lua_pushvalue(L, -3);
        LuaComplexAdd(L);
        lua_pushvalue(L, -3);
        lua_replace(L, -3);
        lua_settable(L, -4);
        lua_remove(L, 1);
      }
      return 1;
    }
    default:
    {
      // first operator is of type (original variable)
      // read this variable and check what the second variable is
      Complex16 *a = (Complex16 *) luaL_checkudata(L, 1, "Complex_Type");
      switch (lua_type(L,len))
      {
        case LUA_TNIL:
        {
          // (original variable) operator (nil)
          return luaL_error(L, "Error: complex number plus nil\n");
        }
        case LUA_TBOOLEAN:
        {
          // (original variable) operator (boolean)
          return luaL_error(L, "Error: complex number plus boolean\n");
        }
        case LUA_TLIGHTUSERDATA:
        {
          // (original variable) operator (ligth user data)
          return luaL_error(L, "Error: complex number plus lightuserdata\n");
        }
        case LUA_TNUMBER:
        {
          // (original variable) operator (number)
          double b;
          Complex16 *res;
          b = luaL_checknumber(L, len);
          res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
          res->Re =  a->Re+b;
          res->Im =  a->Im;
          luaL_setmetatable(L, "Complex_Type");
          return 1;
        }
        case LUA_TSTRING:
        {
          // (original variable) operator (string)
          return luaL_error(L, "Error: complex number plus character string\n");
        }
        case LUA_TUSERDATA:
        {
          // (original variable) operator (user data)
          if (luaL_testudata(L, len, "Complex_Type") != NULL)
          {
            // (original variable) operator (user data - complex)
            Complex16 *b, *res;
            b = (Complex16 *) luaL_checkudata(L, len, "Complex_Type");
            res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
            res->Re = a->Re + b->Re;
            res->Im = a->Im + b->Im;
            luaL_setmetatable(L, "Complex_Type");
            return 1;
          }
          if (luaL_testudata(L, len, "Wavefunction_Type") != NULL)
          {
            // (original variable) operator (user data - wavefunction)
            return luaL_error(L, "Error: complex number plus wavefunction\n");
          }
          if (luaL_testudata(L, len, "Operator_Type") != NULL)
          {
            // (original variable) operator (user data - operator)
            OperatorType *res,*b;
            b = (OperatorType *) luaL_checkudata(L, 2, "Operator_Type");
            res = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
            luaL_setmetatable(L, "Operator_Type");
            OperatorCopy(res, *b);
            OperatorAddLader(res, NULL, 0, a->Re, a->Im);
            return 1;
          }
          if (luaL_testudata(L, len, "Spectra_Type") != NULL)
          {
            // (original variable) operator (user data - spectra)
            uint32_t i,j;
            SpectraTypeArray *res, *b;
            b = (SpectraTypeArray *) luaL_checkudata(L, len, "Spectra_Type");
            res = (SpectraTypeArray *) lua_newuserdata(L, sizeof(Complex16));
            res->NSpectra = b->NSpectra;
            res->Spec=malloc(res->NSpectra*sizeof(SpectraType));
            for(i=0;i<res->NSpectra;i++)
            {
              if(SpectrumAllocate(b->Spec[i], &(res->Spec[i])))
              {
                return luaL_error(L, "Error: not enough memory to allocate new spectra\n");
              }
              for(j=0;j<=res->Spec[i].N;j++)
              {
                res->Spec[i].Re[j]= b->Spec[i].Re[j]+a->Re;
                res->Spec[i].Im[j]= b->Spec[i].Im[j]+a->Im;
              }
            }
            luaL_setmetatable(L, "Spectra_Type");
            return 1;
          }
          return luaL_error(L, "Error: complex number plus unrecognized userdata\n");
        }
        case LUA_TTABLE:
        {
          lua_newtable(L);
          lua_pushnil(L);
          while(lua_next(L, len) != 0)
          {
            LuaComplexAdd(L);
            lua_pushvalue(L, -3);
            lua_insert(L, -2);
            lua_settable(L, -5);
            lua_pop(L, 1);
          }
          return 1;
        }
        default:
        {
          return luaL_error(L, "Error: complex number plus unrecognized type\n");
        }
      }
    }
  }
}

static int LuaComplexSub(lua_State * L)
{
  uint32_t len = (uint32_t)lua_gettop(L);
  switch (lua_type(L,1))
  {
    case LUA_TNIL:
    {
      // (nil) operator (original variable)
      return luaL_error(L, "Error: nil minus complex number\n");
    }
    case LUA_TBOOLEAN:
    {
      // (boolean) operator (original variable)
      return luaL_error(L, "Error: boolean minus complex number\n");
    }
    case LUA_TLIGHTUSERDATA:
    {
      // (ligth user data) operator (original variable)
      return luaL_error(L, "Error: lightuserdata minus complex number\n");
    }
    case LUA_TNUMBER:
    {
      // (number) operator (original variable)
      double a;
      Complex16 *b, *res;
      a = luaL_checknumber(L, 1);
      b = (Complex16 *) luaL_checkudata(L, len, "Complex_Type");
      res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
      res->Re =  a - b->Re;
      res->Im =    - b->Im;
      luaL_setmetatable(L, "Complex_Type");
      return 1;
    }
    case LUA_TSTRING:
    {
      // (string) operator (original variable)
      return luaL_error(L, "Error: character string minus complex number\n");
    }
    case LUA_TTABLE:
    {
      lua_newtable(L);
      lua_pushnil(L);
      while(lua_next(L, 1) != 0)
      {
        lua_insert(L, 1);
        lua_pushvalue(L, -3);
        LuaComplexSub(L);
        lua_pushvalue(L, -3);
        lua_replace(L, -3);
        lua_settable(L, -4);
        lua_remove(L, 1);
      }
      return 1;
    }
    default:
    {
      // first operator is of type (original variable)
      // read this variable and check what the second variable is
      Complex16 *a = (Complex16 *) luaL_checkudata(L, 1, "Complex_Type");
      switch (lua_type(L,len))
      {
        case LUA_TNIL:
        {
          // (original variable) operator (nil)
          return luaL_error(L, "Error: complex number minus nil\n");
        }
        case LUA_TBOOLEAN:
        {
          // (original variable) operator (boolean)
          return luaL_error(L, "Error: complex number minus boolean\n");
        }
        case LUA_TLIGHTUSERDATA:
        {
          // (original variable) operator (ligth user data)
          return luaL_error(L, "Error: complex number minus lightuserdata\n");
        }
        case LUA_TNUMBER:
        {
          // (original variable) operator (number)
          double b;
          Complex16 *res;
          b = luaL_checknumber(L, len);
          res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
          res->Re =  a->Re-b;
          res->Im =  a->Im;
          luaL_setmetatable(L, "Complex_Type");
          return 1;
        }
        case LUA_TSTRING:
        {
          // (original variable) operator (string)
          return luaL_error(L, "Error: complex number minus character string\n");
        }
        case LUA_TUSERDATA:
        {
          // (original variable) operator (user data)
          if (luaL_testudata(L, len, "Complex_Type") != NULL)
          {
            // (original variable) operator (user data - complex)
            Complex16 *b, *res;
            b = (Complex16 *) luaL_checkudata(L, len, "Complex_Type");
            res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
            res->Re = a->Re - b->Re;
            res->Im = a->Im - b->Im;
            luaL_setmetatable(L, "Complex_Type");
            return 1;
          }
          if (luaL_testudata(L, len, "Wavefunction_Type") != NULL)
          {
            // (original variable) operator (user data - wavefunction)
            return luaL_error(L, "Error: complex number minus wavefunction\n");
          }
          if (luaL_testudata(L, len, "Operator_Type") != NULL)
          {
            // (original variable) operator (user data - operator)
            OperatorType *res,*b;
            b = (OperatorType *) luaL_checkudata(L, 2, "Operator_Type");
            res = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
            luaL_setmetatable(L, "Operator_Type");
            OperatorCopy(res, *b);
            OperatorScale(res, -1, 0);
            OperatorAddLader(res, NULL, 0, a->Re, a->Im);
            return 1;
          }
          if (luaL_testudata(L, len, "Spectra_Type") != NULL)
          {
            // (original variable) operator (user data - spectra)
            uint32_t i,j;
            SpectraTypeArray *res, *b;
            b = (SpectraTypeArray *) luaL_checkudata(L, len, "Spectra_Type");
            res = (SpectraTypeArray *) lua_newuserdata(L, sizeof(Complex16));
            res->NSpectra = b->NSpectra;
            res->Spec=malloc(res->NSpectra*sizeof(SpectraType));
            for(i=0;i<res->NSpectra;i++)
            {
              if(SpectrumAllocate(b->Spec[i], &(res->Spec[i])))
              {
                return luaL_error(L, "Error: not enough memory to allocate new spectra\n");
              }
              for(j=0;j<=res->Spec[i].N;j++)
              {
                res->Spec[i].Re[j]= b->Spec[i].Re[j]-a->Re;
                res->Spec[i].Im[j]= b->Spec[i].Im[j]-a->Im;
              }
            }
            luaL_setmetatable(L, "Spectra_Type");
            return 1;
          }
          return luaL_error(L, "Error: complex number minus unrecognized userdata\n");
        }
        case LUA_TTABLE:
        {
          lua_newtable(L);
          lua_pushnil(L);
          while(lua_next(L, len) != 0)
          {
            LuaComplexSub(L);
            lua_pushvalue(L, -3);
            lua_insert(L, -2);
            lua_settable(L, -5);
            lua_pop(L, 1);
          }
          return 1;
        }
        default:
        {
          return luaL_error(L, "Error: complex number minus unrecognized type\n");
        }
      }
    }
  }
}

static int LuaComplexMul(lua_State * L)
{
  uint32_t len = (uint32_t)lua_gettop(L);
  switch (lua_type(L,1))
  {
    case LUA_TNIL:
    {
      // (nil) operator (original variable)
      return luaL_error(L, "Error: nil multiplied with complex number\n");
    }
    case LUA_TBOOLEAN:
    {
      // (boolean) operator (original variable)
      return luaL_error(L, "Error: boolean multiplied with complex number\n");
    }
    case LUA_TLIGHTUSERDATA:
    {
      // (ligth user data) operator (original variable)
      return luaL_error(L, "Error: lightuserdata multiplied with complex number\n");
    }
    case LUA_TNUMBER:
    {
      // (number) operator (original variable)
      double a;
      Complex16 *b, *res;
      a = luaL_checknumber(L, 1);
      b = (Complex16 *) luaL_checkudata(L, len, "Complex_Type");
      res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
      res->Re = a * b->Re;
      res->Im = a * b->Im;
      luaL_setmetatable(L, "Complex_Type");
      return 1;
    }
    case LUA_TSTRING:
    {
      // (string) operator (original variable)
      return luaL_error(L, "Error: character string multiplied with complex number\n");
    }
    case LUA_TTABLE:
    {
      lua_newtable(L);
      lua_pushnil(L);
      while(lua_next(L, 1) != 0)
      {
        lua_insert(L, 1);
        lua_pushvalue(L, -3);
        LuaComplexMul(L);
        lua_pushvalue(L, -3);
        lua_replace(L, -3);
        lua_settable(L, -4);
        lua_remove(L, 1);
      }
      return 1;
    }
    default:
    {
      // first operator is of type (original variable)
      // read this variable and check what the second variable is
      Complex16 *a = (Complex16 *) luaL_checkudata(L, 1, "Complex_Type");
      switch (lua_type(L,len))
      {
        case LUA_TNIL:
        {
          // (original variable) operator (nil)
          return luaL_error(L, "Error: complex number multiplied with nil\n");
        }
        case LUA_TBOOLEAN:
        {
          // (original variable) operator (boolean)
          return luaL_error(L, "Error: complex number multiplied with boolean\n");
        }
        case LUA_TLIGHTUSERDATA:
        {
          // (original variable) operator (ligth user data)
          return luaL_error(L, "Error: complex number multiplied with lightuserdata\n");
        }
        case LUA_TNUMBER:
        {
          // (original variable) operator (number)
          double b;
          Complex16 *res;
          b = luaL_checknumber(L, len);
          res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
          res->Re = a->Re * b;
          res->Im = a->Im * b;
          luaL_setmetatable(L, "Complex_Type");
          return 1;
        }
        case LUA_TSTRING:
        {
          // (original variable) operator (string)
          return luaL_error(L, "Error: complex number multiplied with character string\n");
        }
        case LUA_TUSERDATA:
        {
          // (original variable) operator (user data)
          if (luaL_testudata(L, len, "Complex_Type") != NULL)
          {
            // (original variable) operator (user data - complex)
            Complex16 *b, *res;
            b = (Complex16 *) luaL_checkudata(L, len, "Complex_Type");
            res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
            res->Re = (a->Re * b->Re - a->Im * b->Im);
            res->Im = (a->Im * b->Re + a->Re * b->Im);
            luaL_setmetatable(L, "Complex_Type");
            return 1;
          }
          if (luaL_testudata(L, len, "Wavefunction_Type") != NULL)
          {
            // (original variable) operator (user data - wavefunction)
            WaveFunctionType *b, *res;
            res = (WaveFunctionType *) lua_newuserdata(L, sizeof(WaveFunctionType));
            b = (WaveFunctionType *) luaL_checkudata(L, 2, "Wavefunction_Type");
            WaveFunctionCopy(*b, res);
            WaveFunctionMultiply(res, a->Re, a->Im);
            luaL_setmetatable(L, "Wavefunction_Type");
            return 1;
          }
          if (luaL_testudata(L, len, "Operator_Type") != NULL)
          {
            // (original variable) operator (user data - operator)
            OperatorType *b, *res;
            res = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
            b = (OperatorType *) luaL_checkudata(L, 2, "Operator_Type");
            OperatorCopy(res, *b);
            OperatorScale(res, a->Re, a->Im);
            luaL_setmetatable(L, "Operator_Type");
            return 1;
          }
          if (luaL_testudata(L, len, "Spectra_Type") != NULL)
          {
            // (original variable) operator (user data - spectra)
            uint32_t i;
            SpectraTypeArray *res, *b;
            b = (SpectraTypeArray *) luaL_checkudata(L, len, "Spectra_Type");
            res = (SpectraTypeArray *) lua_newuserdata(L, sizeof(SpectraTypeArray));
            
            if(CopySpectraArray(b, res)) return luaL_error(L, "Error: not enough memory to allocate new spectra\n");
            for(i=0;i<res->NSpectra;i++)
            {
              SpectraScalarMultiply(&(res->Spec[i]), a->Re, a->Im);
            }
            luaL_setmetatable(L, "Spectra_Type");
            return 1;
          }
          return luaL_error(L, "Error: complex number multiplied with unrecognized userdata\n");
        }
        case LUA_TTABLE:
        {
          lua_newtable(L);
          lua_pushnil(L);
          while(lua_next(L, len) != 0)
          {
            LuaComplexMul(L);
            lua_pushvalue(L, -3);
            lua_insert(L, -2);
            lua_settable(L, -5);
            lua_pop(L, 1);
          }
          return 1;
        }
        default:
        {
          return luaL_error(L, "Error: complex number multiplied with unrecognized type\n");
        }
      }
    }
  }
}

static int LuaComplexDiv(lua_State * L)
{
  uint32_t len = (uint32_t)lua_gettop(L);
  switch (lua_type(L,1))
  {
    case LUA_TNIL:
    {
      // (nil) operator (original variable)
      return luaL_error(L, "Error: nil divided by complex number\n");
    }
    case LUA_TBOOLEAN:
    {
      // (boolean) operator (original variable)
      return luaL_error(L, "Error: boolean divided by complex number\n");
    }
    case LUA_TLIGHTUSERDATA:
    {
      // (ligth user data) operator (original variable)
      return luaL_error(L, "Error: lightuserdata divided by complex number\n");
    }
    case LUA_TNUMBER:
    {
      // (number) operator (original variable)
      double a;
      Complex16 *b, *res;
      a = luaL_checknumber(L, 1);
      b = (Complex16 *) luaL_checkudata(L, len, "Complex_Type");
      res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
      res->Re =  a * b->Re / (Sqr(b->Im) + Sqr(b->Re));
      res->Im = -a * b->Im / (Sqr(b->Im) + Sqr(b->Re));
      luaL_setmetatable(L, "Complex_Type");
      return 1;
    }
    case LUA_TSTRING:
    {
      // (string) operator (original variable)
      return luaL_error(L, "Error: character string divided by complex number\n");
    }
    case LUA_TTABLE:
    {
      lua_newtable(L);
      lua_pushnil(L);
      while(lua_next(L, 1) != 0)
      {
        lua_insert(L, 1);
        lua_pushvalue(L, -3);
        LuaComplexDiv(L);
        lua_pushvalue(L, -3);
        lua_replace(L, -3);
        lua_settable(L, -4);
        lua_remove(L, 1);
      }
      return 1;
    }
    default:
    {
      // first operator is of type (original variable)
      // read this variable and check what the second variable is
      Complex16 *a = (Complex16 *) luaL_checkudata(L, 1, "Complex_Type");
      switch (lua_type(L,len))
      {
        case LUA_TNIL:
        {
          // (original variable) operator (nil)
          return luaL_error(L, "Error: complex number divided by nil\n");
        }
        case LUA_TBOOLEAN:
        {
          // (original variable) operator (boolean)
          return luaL_error(L, "Error: complex number divided by boolean\n");
        }
        case LUA_TLIGHTUSERDATA:
        {
          // (original variable) operator (ligth user data)
          return luaL_error(L, "Error: complex number divided by lightuserdata\n");
        }
        case LUA_TNUMBER:
        {
          // (original variable) operator (number)
          double b;
          Complex16 *res;
          b = luaL_checknumber(L, len);
          res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
          res->Re =  a->Re/b;
          res->Im =  a->Im/b;
          luaL_setmetatable(L, "Complex_Type");
          return 1;
        }
        case LUA_TSTRING:
        {
          // (original variable) operator (string)
          return luaL_error(L, "Error: complex number divided by character string\n");
        }
        case LUA_TUSERDATA:
        {
          // (original variable) operator (user data)
          if (luaL_testudata(L, len, "Complex_Type") != NULL)
          {
            // (original variable) operator (user data - complex)
            Complex16 *b, *res;
            b = (Complex16 *) luaL_checkudata(L, len, "Complex_Type");
            res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
            res->Re = (a->Re * b->Re + a->Im * b->Im) / (Sqr(b->Im) + Sqr(b->Re));
            res->Im = (a->Im * b->Re - a->Re * b->Im) / (Sqr(b->Im) + Sqr(b->Re));
            luaL_setmetatable(L, "Complex_Type");
            return 1;
          }
          if (luaL_testudata(L, len, "Wavefunction_Type") != NULL)
          {
            // (original variable) operator (user data - wavefunction)
            return luaL_error(L, "Error: complex number divided by wavefunction\n");
          }
          if (luaL_testudata(L, len, "Operator_Type") != NULL)
          {
            // (original variable) operator (user data - operator)
            return luaL_error(L, "Error: complex number divided by operator\n");
          }
          if (luaL_testudata(L, len, "Spectra_Type") != NULL)
          {
            // (original variable) operator (user data - spectra)
            uint32_t i,j;
            SpectraTypeArray *res, *b;
            b = (SpectraTypeArray *) luaL_checkudata(L, len, "Spectra_Type");
            res = (SpectraTypeArray *) lua_newuserdata(L, sizeof(Complex16));
            res->NSpectra = b->NSpectra;
            res->Spec = malloc(res->NSpectra * sizeof(SpectraType));
            for(i=0;i<res->NSpectra;i++)
            {
              if(SpectrumAllocate(b->Spec[i], &(res->Spec[i])))
              {
                return luaL_error(L, "Error: not enough memory to allocate new spectra\n");
              }
              for(j=0;j<=res->Spec[i].N;j++)
              {
                res->Spec[i].Re[j]= (a->Re*b->Spec[i].Re[j] + a->Im*b->Spec[i].Im[j])/(Sqr(b->Spec[i].Re[j])+Sqr(b->Spec[i].Im[j]));
                res->Spec[i].Im[j]= (a->Im*b->Spec[i].Re[j] - a->Re*b->Spec[i].Im[j])/(Sqr(b->Spec[i].Re[j])+Sqr(b->Spec[i].Im[j]));
              }
            }
            luaL_setmetatable(L, "Spectra_Type");
            return 1;
          }
          return luaL_error(L, "Error: complex number divided by unrecognized userdata\n");
        }
        case LUA_TTABLE:
        {
          lua_newtable(L);
          lua_pushnil(L);
          while(lua_next(L, len) != 0)
          {
            LuaComplexDiv(L);
            lua_pushvalue(L, -3);
            lua_insert(L, -2);
            lua_settable(L, -5);
            lua_pop(L, 1);
          }
          return 1;
        }
        default:
        {
          return luaL_error(L, "Error: complex number divided by unrecognized type\n");
        }
      }
    }
  }
}

static int LuaComplexPow(lua_State * L)
{
  uint32_t len = (uint32_t)lua_gettop(L);
  switch (lua_type(L,1))
  {
    case LUA_TNIL:
    {
      // (nil) operator (original variable)
      return luaL_error(L, "Error: nil to the power complex number\n");
    }
    case LUA_TBOOLEAN:
    {
      // (boolean) operator (original variable)
      return luaL_error(L, "Error: boolean to the power complex number\n");
    }
    case LUA_TLIGHTUSERDATA:
    {
      // (ligth user data) operator (original variable)
      return luaL_error(L, "Error: lightuserdata to the power complex number\n");
    }
    case LUA_TNUMBER:
    {
      // (number) operator (original variable)
      double a;
      Complex16 *b, *res;
      a = luaL_checknumber(L, 1);
      b = (Complex16 *) luaL_checkudata(L, len, "Complex_Type");
      res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
      res->Re = exp(log(a) * b->Re) * cos(log(a) * b->Im);
      res->Im = exp(log(a) * b->Re) * sin(log(a) * b->Im);
      luaL_setmetatable(L, "Complex_Type");
      return 1;
    }
    case LUA_TSTRING:
    {
      // (string) operator (original variable)
      return luaL_error(L, "Error: character string to the power complex number\n");
    }
    case LUA_TTABLE:
    {
      lua_newtable(L);
      lua_pushnil(L);
      while(lua_next(L, 1) != 0)
      {
        lua_insert(L, 1);
        lua_pushvalue(L, -3);
        LuaComplexPow(L);
        lua_pushvalue(L, -3);
        lua_replace(L, -3);
        lua_settable(L, -4);
        lua_remove(L, 1);
      }
      return 1;
    }
    default:
    {
      // first operator is of type (original variable)
      // read this variable and check what the second variable is
      Complex16 *a = (Complex16 *) luaL_checkudata(L, 1, "Complex_Type");
      switch (lua_type(L,len))
      {
        case LUA_TNIL:
        {
          // (original variable) operator (nil)
          return luaL_error(L, "Error: complex number to the power nil\n");
        }
        case LUA_TBOOLEAN:
        {
          // (original variable) operator (boolean)
          return luaL_error(L, "Error: complex number to the power boolean\n");
        }
        case LUA_TLIGHTUSERDATA:
        {
          // (original variable) operator (ligth user data)
          return luaL_error(L, "Error: complex number to the power lightuserdata\n");
        }
        case LUA_TNUMBER:
        {
          // (original variable) operator (number)
          double b;
          Complex16 *res;
          double complex za, zb, zc;
          b = luaL_checknumber(L, len);
          res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
          za = a->Re + a->Im * I;
          zb = b;
          zc = cpow(za,zb);
          res->Re=creal(zc);
          res->Im=cimag(zc);
          luaL_setmetatable(L, "Complex_Type");
          return 1;
        }
        case LUA_TSTRING:
        {
          // (original variable) operator (string)
          return luaL_error(L, "Error: complex number to the power string\n");
        }
        case LUA_TUSERDATA:
        {
          // (original variable) operator (user data)
          if (luaL_testudata(L, len, "Complex_Type") != NULL)
          {
            // (original variable) operator (user data - complex)
            Complex16 *res, *b;
            double complex za, zb, zc;
            b = (Complex16 *) luaL_checkudata(L, len, "Complex_Type");
            res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
            za = a->Re + a->Im * I;
            zb = b->Re + b->Im * I;
            zc = cpow(za,zb);
            res->Re=creal(zc);
            res->Im=cimag(zc);
            luaL_setmetatable(L, "Complex_Type");
            return 1;
          }
          if (luaL_testudata(L, len, "Wavefunction_Type") != NULL)
          {
            // (original variable) operator (user data - wavefunction)
            return luaL_error(L, "Error: complex number to the power wavefunction\n");
          }
          if (luaL_testudata(L, len, "Operator_Type") != NULL)
          {
            // (original variable) operator (user data - operator)
            return luaL_error(L, "Error: complex number to the power operator\n");
          }
          if (luaL_testudata(L, len, "Spectra_Type") != NULL)
          {
            // (original variable) operator (user data - spectra)
            return luaL_error(L, "Error: complex number to the power spectra\n");
          }
          return luaL_error(L, "Error: complex number to the power unrecognized userdata\n");
        }
        case LUA_TTABLE:
        {
          lua_newtable(L);
          lua_pushnil(L);
          while(lua_next(L, len) != 0)
          {
            LuaComplexPow(L);
            lua_pushvalue(L, -3);
            lua_insert(L, -2);
            lua_settable(L, -5);
            lua_pop(L, 1);
          }
          return 1;
        }
        default:
        {
          return luaL_error(L, "Error: complex number to the power unrecognized type\n");
        }
      }
    }
  }
}

static int LuaComplexEq(lua_State * L)
{
  Complex16 *l,*r;
  l = (Complex16 *) luaL_checkudata(L, 1, "Complex_Type");
  r = (Complex16 *) luaL_checkudata(L, 2, "Complex_Type");
  if((l->Re==r->Re)&&(l->Im==r->Im))
  {
    lua_pushboolean (L, 1);
  }
  else
  {
    lua_pushboolean (L, 0);
  }
  return 1;
}

static int LuaComplexUnm(lua_State * L)
{
  Complex16 *res;
  Complex16 *a = (Complex16 *) luaL_checkudata(L, 1, "Complex_Type");
  res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
  res->Re = -a->Re;
  res->Im = -a->Im;
  luaL_setmetatable(L, "Complex_Type");
  return 1;
}

static int LuaComplexExp(lua_State * L)
{
  Complex16 *res;
  Complex16 *a = (Complex16 *) luaL_checkudata(L, 1, "Complex_Type");
  res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
  res->Re = exp(a->Re) * cos(a->Im);
  res->Im = exp(a->Re) * sin(a->Im);
  luaL_setmetatable(L, "Complex_Type");
  return 1;
}

static int LuaComplexRe(lua_State * L)
{
  if(lua_isnumber(L, 1))
  {
    double a = luaL_checknumber(L, 1);
    lua_pushnumber(L, a);
  }
  else if (luaL_testudata(L, 1, "Complex_Type") != NULL)
  {
    Complex16 *a = (Complex16 *) luaL_checkudata(L, 1, "Complex_Type");
    lua_pushnumber(L, a->Re);
  }
  else
  {
    return luaL_error(L, "Can only take complex part of (complex) numbers");
  }
  return 1;
}

static int LuaComplexIm(lua_State * L)
{
  if(lua_isnumber(L, 1))
  {
    lua_pushnumber(L, 0);
  }
  else if (luaL_testudata(L, 1, "Complex_Type") != NULL)
  {
    Complex16 *a = (Complex16 *) luaL_checkudata(L, 1, "Complex_Type");
    lua_pushnumber(L, a->Im);
  }
  else
  {
    return luaL_error(L, "Can only take complex part of (complex) numbers");
  }
  return 1;
}

static int LuaComplexToString(lua_State * L)
{
  Complex16 *v = (Complex16 *) luaL_checkudata(L, 1, "Complex_Type");
  if(v->Im==0)
  {
    lua_pushfstring(L, "%f", v->Re);
  }
  else if(v->Im<0)
  {
    lua_pushfstring(L, "(%f - %f I)", v->Re, -v->Im);
  }
  else
  {
    lua_pushfstring(L, "(%f + %f I)", v->Re, v->Im);
  }
  return 1;
}

static const struct luaL_Reg Complex_meta[] = {
  {"__index", LuaComplexIndex},
  {"__add", LuaComplexAdd},
  {"__sub", LuaComplexSub},
  {"__mul", LuaComplexMul},
  {"__div", LuaComplexDiv},
  {"__pow", LuaComplexPow},
  {"__unm", LuaComplexUnm},
  {"__eq", LuaComplexEq},
  {"__newindex", LuaComplexNewIndex},
  {"__tostring", LuaComplexToString},
  {NULL, NULL}
};


static const struct luaL_Reg Complex_func[] = {
  {"New", LuaComplexNew},
  {"I", LuaComplexI},
  {"Exp", LuaComplexExp},
  {"Re", LuaComplexRe},
  {"re", LuaComplexRe},
  {"Im", LuaComplexIm},
  {"im", LuaComplexIm},
  {NULL, NULL}
};



int luaopen_Complex(lua_State * L)
{
  luaL_newmetatable(L, "Complex_Type");
  luaL_setfuncs(L, Complex_meta, 0);
  luaL_newlib(L, Complex_func);
  lua_setglobal(L, "Complex");
  return 1;
}
