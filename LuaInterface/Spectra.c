/*  Print a set of spectra
 *
 */
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
    
    VerbosityPrintf(0x2,"Spectra printed to file: %s\n",Options.stream);
    
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

/*  Return a subset of the spectra
 *  on stack a SpectrumTypeArray and either a number or a list of numbers (integers)
 */
static int LuaSpectrumElement(lua_State * L)
{
  SpectraTypeArray *a = (SpectraTypeArray *)currentudata;
  uint32_t len = (uint32_t)lua_gettop(L);
  if(len==1)
  {
    switch (lua_type(L,len))
    {
      case LUA_TNIL:
      {
        // (original variable) operator (nil)
        return luaL_error(L, "Error: can not find nill element in spectrum, expected a number or list of numbers\n");
      }
      case LUA_TBOOLEAN:
      {
        // (original variable) operator (boolean)
        return luaL_error(L, "Error: can not find boolean element in spectrum, expected a number or list of numbers\n");
      }
      case LUA_TLIGHTUSERDATA:
      {
        // (original variable) operator (ligth user data)
        return luaL_error(L, "Error: can not find lightuserdata element in spectrum, expected a number or list of numbers\n");
      }
      case LUA_TNUMBER:
      {
        // return element number
        uint32_t i;
        SpectraTypeArray *res;
        i = (uint32_t)luaL_checknumber(L, len);
        if(i>a->NSpectra)
        {
          return luaL_error(L, "Error: no element %d in spectra array with %d elements\n",i,a->NSpectra);
        }
        res = (SpectraTypeArray *) lua_newuserdata(L, sizeof(Complex16));
        res->NSpectra = 1;
        res->Spec=malloc(res->NSpectra*sizeof(SpectraType));
        if(SpectrumAllocate(a->Spec[i-1], &(res->Spec[0])))
        {
          return luaL_error(L, "Error: not enough memory to allocate new spectra\n");
        }
        CopySpectra(a->Spec[i-1], &(res->Spec[0]));
        luaL_setmetatable(L, "Spectra_Type");
        return 1;
      }
      case LUA_TSTRING:
      {
        // (original variable) operator (string)
        return luaL_error(L, "Error: can not find string element in spectrum, expected a number or list of numbers\n");
      }
      case LUA_TUSERDATA:
      {
        // (original variable) operator (user data)
        if (luaL_testudata(L, len, "Complex_Type") != NULL)
        {
          // (original variable) operator (user data - complex)
          return luaL_error(L, "Error: can not find complex valued number element in spectrum, expected a number or list of numbers\n");
        }
        if (luaL_testudata(L, len, "Wavefunction_Type") != NULL)
        {
          // (original variable) operator (user data - wavefunction)
          return luaL_error(L, "Error: can not find wavefunction element in spectrum, expected a number or list of numbers\n");
        }
        if (luaL_testudata(L, len, "Operator_Type") != NULL)
        {
          // (original variable) operator (user data - operator)
          return luaL_error(L, "Error: can not find operator element in spectrum, expected a number or list of numbers\n");
        }
        if (luaL_testudata(L, len, "Spectra_Type") != NULL)
        {
          // (original variable) operator (user data - spectra)
          return luaL_error(L, "Error: can not find spectrum element in spectrum, expected a number or list of numbers\n");
        }
        return luaL_error(L, "Error: unrecognized type of element\n");
      }
      case LUA_TTABLE:
      {
        uint32_t i,j;
        SpectraTypeArray *res;
        res = (SpectraTypeArray *) lua_newuserdata(L, sizeof(Complex16));
        res->NSpectra = (uint32_t) lua_objlen(L, len);
        res->Spec=malloc(res->NSpectra*sizeof(SpectraType));
        for(j=0;j<res->NSpectra;j++)
        {
          lua_rawgeti(L, len, j+1);
          i = (uint32_t)luaL_checknumber(L, -1);
          lua_pop(L,1);
          if(i>a->NSpectra)
          {
            return luaL_error(L, "Error: no element %d in spectra array with %d elements\n",i,a->NSpectra);
          }
          if(SpectrumAllocate(a->Spec[i-1], &(res->Spec[j])))
          {
            return luaL_error(L, "Error: not enough memory to allocate new spectra\n");
          }
          CopySpectra(a->Spec[i-1], &(res->Spec[j]));
        }
        luaL_setmetatable(L, "Spectra_Type");
        return 1;
      }
      default:
      {
        return luaL_error(L, "Error: unrecognized type of element\n");
      }
    }
  }
  else
  {
    return luaL_error(L, "Error in .Element(arg) unrecognized arguments in Element()\n");
  }
  return 0;
}


/*  Return a sum of the spectra
 *  on stack a SpectrumTypeArray and a list of numbers (doubles or complex) with the length of the list equal to the number of elements in the spectra
 *  list of lists are allowed
 */
static int LuaSpectrumSum(lua_State * L)
{
  uint32_t i,j,k;
  SpectraTypeArray *res;
  Complex16 alpha;
  SpectraTypeArray *a = (SpectraTypeArray *)currentudata;
  uint32_t len = (uint32_t)lua_gettop(L);
  res = (SpectraTypeArray *) lua_newuserdata(L, sizeof(Complex16));
  luaL_setmetatable(L, "Spectra_Type");
  res->NSpectra = len;
  res->Spec=malloc(res->NSpectra*sizeof(SpectraType));
  res->Spec[0].Emin    = a->Spec[0].Emin;
  res->Spec[0].Emax    = a->Spec[0].Emax;
  res->Spec[0].EminPole= a->Spec[0].EminPole;
  res->Spec[0].EmaxPole= a->Spec[0].EmaxPole;
  res->Spec[0].dE      = a->Spec[0].dE;
  res->Spec[0].Gamma   = a->Spec[0].Gamma;
  res->Spec[0].ZerothMoment = 0;
  res->Spec[0].FirstMoment = 0;
  res->Spec[0].a[0] = 0;
  res->Spec[0].a[1] = 0;
  res->Spec[0].a[2] = 0;
  res->Spec[0].a[3] = 0;
  res->Spec[0].a[4] = 0;
  res->Spec[0].b[0] = 0;
  res->Spec[0].b[1] = 0;
  res->Spec[0].b[2] = 0;
  res->Spec[0].b[3] = 0;
  res->Spec[0].b[4] = 0;
  for(i=1;i<a->NSpectra;i++)
  {
    if(a->Spec[i].Emin<res->Spec[0].Emin) res->Spec[0].Emin=a->Spec[i].Emin;
    if(a->Spec[i].Emax>res->Spec[0].Emax) res->Spec[0].Emax=a->Spec[i].Emax;
    if(a->Spec[i].EminPole<res->Spec[0].EminPole) res->Spec[0].EminPole=a->Spec[i].EminPole;
    if(a->Spec[i].EmaxPole>res->Spec[0].EmaxPole) res->Spec[0].EmaxPole=a->Spec[i].EmaxPole;
    if(a->Spec[i].dE<res->Spec[0].dE) res->Spec[0].dE=a->Spec[i].dE;
  }
  res->Spec[0].N = (uint32_t) (0.5 + (res->Spec[0].Emax-res->Spec[0].Emin)/res->Spec[0].dE);
  for(j=0;j<res->NSpectra;j++)
  {
    if(lua_type(L, j+1)!=LUA_TTABLE)
    {
      return luaL_error(L, "Error in .Sum(arg_1,arg_2,...) arg_%d should be a table\n",j+1);
    }
    if(((uint32_t) lua_objlen(L, j+1))!=a->NSpectra)
    {
      return luaL_error(L, "Error in .Sum(arg_1,arg_2,...) arg_%d should be a table of length %d\n",j+1,a->NSpectra);
    }
    if(SpectrumAllocate(res->Spec[0], &(res->Spec[j])))
    {
      return luaL_error(L, "Error: not enough memory to allocate new spectra\n");
    }
    for(k=0;k<a->NSpectra;k++)
    {
      lua_rawgeti(L, j+1, k+1);
      // get either a number or complex number from last position on stack
      if(lua_type(L,-1)==LUA_TNUMBER)
      {
        alpha.Im=0;
        alpha.Re = luaL_checknumber(L, -1);
      }
      else if((lua_type(L,-1)==LUA_TUSERDATA)&&(luaL_testudata(L, -1, "Complex_Type") != NULL))
      {
        alpha = *(Complex16 *) luaL_checkudata(L, -1, "Complex_Type");
      }
      else
      {
        return luaL_error(L, "Error in .Sum(arg_1,arg_2,...) elemnt %d in arg_%d should be a real or complex number\n",k+1,j+1);
      }
      lua_pop(L,1);
      SafeAddMultipliedSpectra(&(res->Spec[j]), a->Spec[k], alpha.Re, alpha.Im);
    }
  }
  return 1;
}

static int LuaSpectraAdd(lua_State * L)
{
  uint32_t len = (uint32_t)lua_gettop(L);
  switch (lua_type(L,1))
  {
    case LUA_TNIL:
    {
      // (nil) operator (original variable)
      return luaL_error(L, "Error: nil plus spectrum\n");
    }
    case LUA_TBOOLEAN:
    {
      // (boolean) operator (original variable)
      return luaL_error(L, "Error: boolean plus spectrum\n");
    }
    case LUA_TLIGHTUSERDATA:
    {
      // (ligth user data) operator (original variable)
      return luaL_error(L, "Error: lightuserdata plus spectrum\n");
    }
    case LUA_TNUMBER:
    {
      // (number) operator (original variable)
      uint32_t i,j;
      double a;
      SpectraTypeArray *res, *b;
      a = luaL_checknumber(L, 1);
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
          res->Spec[i].Re[j]= b->Spec[i].Re[j]+a;
          res->Spec[i].Im[j]= b->Spec[i].Im[j];
        }
      }
      luaL_setmetatable(L, "Spectra_Type");
      return 1;
    }
    case LUA_TSTRING:
    {
      // (string) operator (original variable)
      return luaL_error(L, "Error: character string plus spectrum\n");
    }
    case LUA_TTABLE:
    {
      lua_newtable(L);
      lua_pushnil(L);
      while(lua_next(L, 1) != 0)
      {
        lua_insert(L, 1);
        lua_pushvalue(L, -3);
        LuaSpectraAdd(L);
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
      SpectraTypeArray *a = (SpectraTypeArray *) luaL_checkudata(L, 1, "Spectra_Type");
      switch (lua_type(L,len))
      {
        case LUA_TNIL:
        {
          // (original variable) operator (nil)
          return luaL_error(L, "Error: spectra plus nil\n");
        }
        case LUA_TBOOLEAN:
        {
          // (original variable) operator (boolean)
          return luaL_error(L, "Error: spectra plus boolean\n");
        }
        case LUA_TLIGHTUSERDATA:
        {
          // (original variable) operator (ligth user data)
          return luaL_error(L, "Error: spectra plus lightuserdata\n");
        }
        case LUA_TNUMBER:
        {
          // (original variable) operator (number)
          uint32_t i,j;
          double b;
          SpectraTypeArray *res;
          b = luaL_checknumber(L, len);
          res = (SpectraTypeArray *) lua_newuserdata(L, sizeof(Complex16));
          res->NSpectra = a->NSpectra;
          res->Spec = malloc(res->NSpectra * sizeof(SpectraType));
          for(i=0;i<res->NSpectra;i++)
          {
            if(SpectrumAllocate(a->Spec[i], &(res->Spec[i])))
            {
              return luaL_error(L, "Error: not enough memory to allocate new spectra\n");
            }
            for(j=0;j<=res->Spec[i].N;j++)
            {
              res->Spec[i].Re[j]= a->Spec[i].Re[j]+b;
              res->Spec[i].Im[j]= a->Spec[i].Im[j];
            }
          }
          luaL_setmetatable(L, "Spectra_Type");
          return 1;
        }
        case LUA_TSTRING:
        {
          // (original variable) operator (string)
          return luaL_error(L, "Error: spectra plus string\n");
        }
        case LUA_TUSERDATA:
        {
          // (original variable) operator (user data)
          if (luaL_testudata(L, len, "Complex_Type") != NULL)
          {
            // (original variable) operator (user data - complex)
            uint32_t i,j;
            Complex16 *b;
            SpectraTypeArray *res;
            b = (Complex16 *) luaL_checkudata(L, len, "Complex_Type");
            res = (SpectraTypeArray *) lua_newuserdata(L, sizeof(Complex16));
            res->NSpectra = a->NSpectra;
            res->Spec = malloc(res->NSpectra * sizeof(SpectraType));
            for(i=0;i<res->NSpectra;i++)
            {
              if(SpectrumAllocate(a->Spec[i], &(res->Spec[i])))
              {
                return luaL_error(L, "Error: not enough memory to allocate new spectra\n");
              }
              for(j=0;j<=res->Spec[i].N;j++)
              {
                res->Spec[i].Re[j]= a->Spec[i].Re[j]+b->Re;
                res->Spec[i].Im[j]= a->Spec[i].Im[j]+b->Im;
              }
            }
            luaL_setmetatable(L, "Spectra_Type");
            return 1;
          }
          if (luaL_testudata(L, len, "Wavefunction_Type") != NULL)
          {
            // (original variable) operator (user data - wavefunction)
            return luaL_error(L, "Error: spectra plus wavefunction\n");
          }
          if (luaL_testudata(L, len, "Operator_Type") != NULL)
          {
            // (original variable) operator (user data - operator)
            return luaL_error(L, "Error: spectra plus operator\n");
          }
          if (luaL_testudata(L, len, "Spectra_Type") != NULL)
          {
            // (original variable) operator (user data - spectra)
            uint32_t i,j;
            SpectraTypeArray *res, *b;
            b = (SpectraTypeArray *) luaL_checkudata(L, len, "Spectra_Type");
            if(a->NSpectra!=b->NSpectra) return luaL_error(L, "Error: Adding to sets of spectra with different number of elements (%d) + (%d)\n",a->NSpectra,b->NSpectra);
            res = (SpectraTypeArray *) lua_newuserdata(L, sizeof(Complex16));
            res->NSpectra = a->NSpectra;
            res->Spec=malloc(res->NSpectra*sizeof(SpectraType));
            for(i=0;i<res->NSpectra;i++)
            {
              if((a->Spec[i].N==b->Spec[i].N)&&(a->Spec[i].Emin==b->Spec[i].Emin)&&(a->Spec[i].Emax==b->Spec[i].Emax))
              {
                if(SpectrumAllocate(b->Spec[i], &(res->Spec[i])))
                {
                  return luaL_error(L, "Error: not enough memory to allocate new spectra\n");
                }
                res->Spec[i].EminPole = Min(a->Spec[i].EminPole,b->Spec[i].EminPole);
                res->Spec[i].EmaxPole = Max(a->Spec[i].EmaxPole,b->Spec[i].EmaxPole);
                for(j=0;j<=res->Spec[i].N;j++)
                {
                  res->Spec[i].Re[j]= a->Spec[i].Re[j] + b->Spec[i].Re[j];
                  res->Spec[i].Im[j]= a->Spec[i].Im[j] + b->Spec[i].Im[j];
                }
              }
              else
              {
                return luaL_error(L, "Still to do: add two spectra with different energy ranges, spectra %d\n",i);
              }
            }
            luaL_setmetatable(L, "Spectra_Type");
            return 1;
          }
          return luaL_error(L, "Error: spectra plus unrecognized userdata\n");
        }
        case LUA_TTABLE:
        {
          lua_newtable(L);
          lua_pushnil(L);
          while(lua_next(L, len) != 0)
          {
            LuaSpectraAdd(L);
            lua_pushvalue(L, -3);
            lua_insert(L, -2);
            lua_settable(L, -5);
            lua_pop(L, 1);
          }
          return 1;
        }
        default:
        {
          return luaL_error(L, "Error: spectra plus unrecognized type\n");
        }
      }
    }
  }
}

static int LuaSpectraSub(lua_State * L)
{
  uint32_t len = (uint32_t)lua_gettop(L);
  switch (lua_type(L,1))
  {
    case LUA_TNIL:
    {
      // (nil) operator (original variable)
      return luaL_error(L, "Error: nil minus spectrum\n");
    }
    case LUA_TBOOLEAN:
    {
      // (boolean) operator (original variable)
      return luaL_error(L, "Error: boolean minus spectrum\n");
    }
    case LUA_TLIGHTUSERDATA:
    {
      // (ligth user data) operator (original variable)
      return luaL_error(L, "Error: lightuserdata minus spectrum\n");
    }
    case LUA_TNUMBER:
    {
      // (number) operator (original variable)
      uint32_t i,j;
      double a;
      SpectraTypeArray *res, *b;
      a = luaL_checknumber(L, 1);
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
          res->Spec[i].Re[j]= -b->Spec[i].Re[j]+a;
          res->Spec[i].Im[j]= -b->Spec[i].Im[j];
        }
      }
      luaL_setmetatable(L, "Spectra_Type");
      return 1;
    }
    case LUA_TSTRING:
    {
      // (string) operator (original variable)
      return luaL_error(L, "Error: character string minus spectrum\n");
    }
    case LUA_TTABLE:
    {
      lua_newtable(L);
      lua_pushnil(L);
      while(lua_next(L, 1) != 0)
      {
        lua_insert(L, 1);
        lua_pushvalue(L, -3);
        LuaSpectraSub(L);
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
      SpectraTypeArray *a = (SpectraTypeArray *) luaL_checkudata(L, 1, "Spectra_Type");
      switch (lua_type(L,len))
      {
        case LUA_TNIL:
        {
          // (original variable) operator (nil)
          return luaL_error(L, "Error: spectra minus nil\n");
        }
        case LUA_TBOOLEAN:
        {
          // (original variable) operator (boolean)
          return luaL_error(L, "Error: spectra minus boolean\n");
        }
        case LUA_TLIGHTUSERDATA:
        {
          // (original variable) operator (ligth user data)
          return luaL_error(L, "Error: spectra minus lightuserdata\n");
        }
        case LUA_TNUMBER:
        {
          // (original variable) operator (number)
          uint32_t i,j;
          double b;
          SpectraTypeArray *res;
          b = luaL_checknumber(L, len);
          res = (SpectraTypeArray *) lua_newuserdata(L, sizeof(Complex16));
          res->NSpectra = a->NSpectra;
          res->Spec = malloc(res->NSpectra * sizeof(SpectraType));
          for(i=0;i<res->NSpectra;i++)
          {
            if(SpectrumAllocate(a->Spec[i], &(res->Spec[i])))
            {
              return luaL_error(L, "Error: not enough memory to allocate new spectra\n");
            }
            for(j=0;j<=res->Spec[i].N;j++)
            {
              res->Spec[i].Re[j]= a->Spec[i].Re[j]-b;
              res->Spec[i].Im[j]= a->Spec[i].Im[j];
            }
          }
          luaL_setmetatable(L, "Spectra_Type");
          return 1;
        }
        case LUA_TSTRING:
        {
          // (original variable) operator (string)
          return luaL_error(L, "Error: spectra minus string\n");
        }
        case LUA_TUSERDATA:
        {
          // (original variable) operator (user data)
          if (luaL_testudata(L, len, "Complex_Type") != NULL)
          {
            // (original variable) operator (user data - complex)
            uint32_t i,j;
            Complex16 *b;
            SpectraTypeArray *res;
            b = (Complex16 *) luaL_checkudata(L, len, "Complex_Type");
            res = (SpectraTypeArray *) lua_newuserdata(L, sizeof(Complex16));
            res->NSpectra = a->NSpectra;
            res->Spec = malloc(res->NSpectra * sizeof(SpectraType));
            for(i=0;i<res->NSpectra;i++)
            {
              if(SpectrumAllocate(a->Spec[i], &(res->Spec[i])))
              {
                return luaL_error(L, "Error: not enough memory to allocate new spectra\n");
              }
              for(j=0;j<=res->Spec[i].N;j++)
              {
                res->Spec[i].Re[j]= a->Spec[i].Re[j]-b->Re;
                res->Spec[i].Im[j]= a->Spec[i].Im[j]-b->Im;
              }
            }
            luaL_setmetatable(L, "Spectra_Type");
            return 1;
          }
          if (luaL_testudata(L, len, "Wavefunction_Type") != NULL)
          {
            // (original variable) operator (user data - wavefunction)
            return luaL_error(L, "Error: spectra minus wavefunction\n");
          }
          if (luaL_testudata(L, len, "Operator_Type") != NULL)
          {
            // (original variable) operator (user data - operator)
            return luaL_error(L, "Error: spectra minus operator\n");
          }
          if (luaL_testudata(L, len, "Spectra_Type") != NULL)
          {
            // (original variable) operator (user data - spectra)
            uint32_t i,j;
            SpectraTypeArray *res, *b;
            b = (SpectraTypeArray *) luaL_checkudata(L, len, "Spectra_Type");
            if(a->NSpectra!=b->NSpectra) return luaL_error(L, "Error: Subtracting to sets of spectra with different number of elements (%d) + (%d)\n",a->NSpectra,b->NSpectra);
            res = (SpectraTypeArray *) lua_newuserdata(L, sizeof(Complex16));
            res->NSpectra = a->NSpectra;
            res->Spec = malloc(res->NSpectra * sizeof(SpectraType));
            for(i=0;i<res->NSpectra;i++)
            {
              if((a->Spec[i].N==b->Spec[i].N)&&(a->Spec[i].Emin==b->Spec[i].Emin)&&(a->Spec[i].Emax==b->Spec[i].Emax))
              {
                if(SpectrumAllocate(b->Spec[i], &(res->Spec[i])))
                {
                  return luaL_error(L, "Error: not enough memory to allocate new spectra\n");
                }
                res->Spec[i].EminPole = Min(a->Spec[i].EminPole,b->Spec[i].EminPole);
                res->Spec[i].EmaxPole = Max(a->Spec[i].EmaxPole,b->Spec[i].EmaxPole);
                for(j=0;j<=res->Spec[i].N;j++)
                {
                  res->Spec[i].Re[j]= a->Spec[i].Re[j] - b->Spec[i].Re[j];
                  res->Spec[i].Im[j]= a->Spec[i].Im[j] - b->Spec[i].Im[j];
                }
              }
              else
              {
                return luaL_error(L, "Still to do: subtract two spectra with different energy ranges, spectra %d\n",i);
              }
            }
            luaL_setmetatable(L, "Spectra_Type");
            return 1;
          }
          return luaL_error(L, "Error: spectra minus unrecognized userdata\n");
        }
        case LUA_TTABLE:
        {
          lua_newtable(L);
          lua_pushnil(L);
          while(lua_next(L, len) != 0)
          {
            LuaSpectraSub(L);
            lua_pushvalue(L, -3);
            lua_insert(L, -2);
            lua_settable(L, -5);
            lua_pop(L, 1);
          }
          return 1;
        }
        default:
        {
          return luaL_error(L, "Error: spectra minus unrecognized type\n");
        }
      }
    }
  }
}

static int LuaSpectraMul(lua_State * L)
{
  uint32_t len = (uint32_t)lua_gettop(L);
  switch (lua_type(L,1))
  {
    case LUA_TNIL:
    {
      // (nil) operator (original variable)
      return luaL_error(L, "Error: nil times spectrum\n");
    }
    case LUA_TBOOLEAN:
    {
      // (boolean) operator (original variable)
      return luaL_error(L, "Error: boolean times spectrum\n");
    }
    case LUA_TLIGHTUSERDATA:
    {
      // (ligth user data) operator (original variable)
      return luaL_error(L, "Error: lightuserdata times spectrum\n");
    }
    case LUA_TNUMBER:
    {
      // (number) operator (original variable)
      uint32_t i,j;
      double a;
      SpectraTypeArray *res, *b;
      a = luaL_checknumber(L, 1);
      b = (SpectraTypeArray *) luaL_checkudata(L, len, "Spectra_Type");
      res = (SpectraTypeArray *) lua_newuserdata(L, sizeof(Complex16));
      res->NSpectra = b->NSpectra;
      res->Spec = malloc(res->NSpectra*sizeof(SpectraType));
      for(i=0;i<res->NSpectra;i++)
      {
        if(SpectrumAllocate(b->Spec[i], &(res->Spec[i])))
        {
          return luaL_error(L, "Error: not enough memory to allocate new spectra\n");
        }
        for(j=0;j<=res->Spec[i].N;j++)
        {
          res->Spec[i].Re[j]= b->Spec[i].Re[j]*a;
          res->Spec[i].Im[j]= b->Spec[i].Im[j]*a;
        }
      }
      luaL_setmetatable(L, "Spectra_Type");
      return 1;
    }
    case LUA_TSTRING:
    {
      // (string) operator (original variable)
      return luaL_error(L, "Error: character string times spectrum\n");
    }
    case LUA_TTABLE:
    {
      lua_newtable(L);
      lua_pushnil(L);
      while(lua_next(L, 1) != 0)
      {
        lua_insert(L, 1);
        lua_pushvalue(L, -3);
        LuaSpectraMul(L);
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
      SpectraTypeArray *a = (SpectraTypeArray *) luaL_checkudata(L, 1, "Spectra_Type");
      switch (lua_type(L,len))
      {
        case LUA_TNIL:
        {
          // (original variable) operator (nil)
          return luaL_error(L, "Error: spectra times nil\n");
        }
        case LUA_TBOOLEAN:
        {
          // (original variable) operator (boolean)
          return luaL_error(L, "Error: spectra times boolean\n");
        }
        case LUA_TLIGHTUSERDATA:
        {
          // (original variable) operator (ligth user data)
          return luaL_error(L, "Error: spectra times lightuserdata\n");
        }
        case LUA_TNUMBER:
        {
          // (original variable) operator (number)
          uint32_t i,j;
          double b;
          SpectraTypeArray *res;
          b = luaL_checknumber(L, len);
          res = (SpectraTypeArray *) lua_newuserdata(L, sizeof(Complex16));
          res->NSpectra = a->NSpectra;
          res->Spec = malloc(res->NSpectra * sizeof(SpectraType));
          for(i=0;i<res->NSpectra;i++)
          {
            if(SpectrumAllocate(a->Spec[i], &(res->Spec[i])))
            {
              return luaL_error(L, "Error: not enough memory to allocate new spectra\n");
            }
            for(j=0;j<=res->Spec[i].N;j++)
            {
              res->Spec[i].Re[j]= a->Spec[i].Re[j]*b;
              res->Spec[i].Im[j]= a->Spec[i].Im[j]*b;
            }
          }
          luaL_setmetatable(L, "Spectra_Type");
          return 1;
        }
        case LUA_TSTRING:
        {
          // (original variable) operator (string)
          return luaL_error(L, "Error: spectra times string\n");
        }
        case LUA_TUSERDATA:
        {
          // (original variable) operator (user data)
          if (luaL_testudata(L, len, "Complex_Type") != NULL)
          {
            // (original variable) operator (user data - complex)
            uint32_t i,j;
            Complex16 *b;
            SpectraTypeArray *res;
            b = (Complex16 *) luaL_checkudata(L, len, "Complex_Type");
            res = (SpectraTypeArray *) lua_newuserdata(L, sizeof(Complex16));
            res->NSpectra = a->NSpectra;
            res->Spec = malloc(res->NSpectra * sizeof(SpectraType));
            for(i=0;i<res->NSpectra;i++)
            {
              if(SpectrumAllocate(a->Spec[i], &(res->Spec[i])))
              {
                return luaL_error(L, "Error: not enough memory to allocate new spectra\n");
              }
              for(j=0;j<=res->Spec[i].N;j++)
              {
                res->Spec[i].Re[j]= a->Spec[i].Re[j]*b->Re - a->Spec[i].Im[j]*b->Im;
                res->Spec[i].Im[j]= a->Spec[i].Re[j]*b->Im + a->Spec[i].Im[j]*b->Re;
              }
            }
            luaL_setmetatable(L, "Spectra_Type");
            return 1;
          }
          if (luaL_testudata(L, len, "Wavefunction_Type") != NULL)
          {
            // (original variable) operator (user data - wavefunction)
            return luaL_error(L, "Error: spectra times wavefunction\n");
          }
          if (luaL_testudata(L, len, "Operator_Type") != NULL)
          {
            // (original variable) operator (user data - operator)
            return luaL_error(L, "Error: spectra times operator\n");
          }
          if (luaL_testudata(L, len, "Spectra_Type") != NULL)
          {
            // (original variable) operator (user data - spectra)
            uint32_t i,j;
            SpectraTypeArray *res, *b;
            b = (SpectraTypeArray *) luaL_checkudata(L, len, "Spectra_Type");
            if(a->NSpectra!=b->NSpectra) return luaL_error(L, "Error: Multiplying to sets of spectra with different number of elements (%d) + (%d)\n",a->NSpectra,b->NSpectra);
            res = (SpectraTypeArray *) lua_newuserdata(L, sizeof(Complex16));
            res->NSpectra = a->NSpectra;
            res->Spec = malloc(res->NSpectra * sizeof(SpectraType));
            for(i=0;i<res->NSpectra;i++)
            {
              if((a->Spec[i].N==b->Spec[i].N)&&(a->Spec[i].Emin==b->Spec[i].Emin)&&(a->Spec[i].Emax==b->Spec[i].Emax))
              {
                if(SpectrumAllocate(b->Spec[i], &(res->Spec[i])))
                {
                  return luaL_error(L, "Error: not enough memory to allocate new spectra\n");
                }
                res->Spec[i].EminPole = Min(a->Spec[i].EminPole,b->Spec[i].EminPole);
                res->Spec[i].EmaxPole = Max(a->Spec[i].EmaxPole,b->Spec[i].EmaxPole);
                for(j=0;j<=res->Spec[i].N;j++)
                {
                  res->Spec[i].Re[j]= a->Spec[i].Re[j] * b->Spec[i].Re[j] - a->Spec[i].Im[j] * b->Spec[i].Im[j];
                  res->Spec[i].Im[j]= a->Spec[i].Re[j] * b->Spec[i].Im[j] - a->Spec[i].Im[j] * b->Spec[i].Re[j];
                }
              }
              else
              {
                return luaL_error(L, "Still to do: multiply two spectra with different energy ranges, spectra %d\n",i);
              }
            }
            luaL_setmetatable(L, "Spectra_Type");
            return 1;
          }
          return luaL_error(L, "Error: spectra times unrecognized userdata\n");
        }
        case LUA_TTABLE:
        {
          lua_newtable(L);
          lua_pushnil(L);
          while(lua_next(L, len) != 0)
          {
            LuaSpectraMul(L);
            lua_pushvalue(L, -3);
            lua_insert(L, -2);
            lua_settable(L, -5);
            lua_pop(L, 1);
          }
          return 1;
        }
        default:
        {
          return luaL_error(L, "Error: spectra times unrecognized type\n");
        }
      }
    }
  }
}

static int LuaSpectraDiv(lua_State * L)
{
  uint32_t len = (uint32_t)lua_gettop(L);
  switch (lua_type(L,1))
  {
    case LUA_TNIL:
    {
      // (nil) operator (original variable)
      return luaL_error(L, "Error: nil divided by spectrum\n");
    }
    case LUA_TBOOLEAN:
    {
      // (boolean) operator (original variable)
      return luaL_error(L, "Error: boolean divided by spectrum\n");
    }
    case LUA_TLIGHTUSERDATA:
    {
      // (ligth user data) operator (original variable)
      return luaL_error(L, "Error: lightuserdata divided by spectrum\n");
    }
    case LUA_TNUMBER:
    {
      // (number) operator (original variable)
      uint32_t i,j;
      double a;
      SpectraTypeArray *res, *b;
      a = luaL_checknumber(L, 1);
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
          res->Spec[i].Re[j]=  a*b->Spec[i].Re[j]/(Sqr(b->Spec[i].Re[j])+Sqr(b->Spec[i].Im[j]));
          res->Spec[i].Im[j]= -a*b->Spec[i].Im[j]/(Sqr(b->Spec[i].Re[j])+Sqr(b->Spec[i].Im[j]));
        }
      }
      luaL_setmetatable(L, "Spectra_Type");
      return 1;
    }
    case LUA_TSTRING:
    {
      // (string) operator (original variable)
      return luaL_error(L, "Error: character string divided by spectrum\n");
    }
    case LUA_TTABLE:
    {
      lua_newtable(L);
      lua_pushnil(L);
      while(lua_next(L, 1) != 0)
      {
        lua_insert(L, 1);
        lua_pushvalue(L, -3);
        LuaSpectraDiv(L);
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
      SpectraTypeArray *a = (SpectraTypeArray *) luaL_checkudata(L, 1, "Spectra_Type");
      switch (lua_type(L,len))
      {
        case LUA_TNIL:
        {
          // (original variable) operator (nil)
          return luaL_error(L, "Error: spectra divided by nil\n");
        }
        case LUA_TBOOLEAN:
        {
          // (original variable) operator (boolean)
          return luaL_error(L, "Error: spectra divided by boolean\n");
        }
        case LUA_TLIGHTUSERDATA:
        {
          // (original variable) operator (ligth user data)
          return luaL_error(L, "Error: spectra divided by lightuserdata\n");
        }
        case LUA_TNUMBER:
        {
          // (original variable) operator (number)
          uint32_t i,j;
          double b;
          SpectraTypeArray *res;
          b = luaL_checknumber(L, len);
          res = (SpectraTypeArray *) lua_newuserdata(L, sizeof(Complex16));
          res->NSpectra = a->NSpectra;
          res->Spec = malloc(res->NSpectra * sizeof(SpectraType));
          for(i=0;i<res->NSpectra;i++)
          {
            if(SpectrumAllocate(a->Spec[i], &(res->Spec[i])))
            {
              return luaL_error(L, "Error: not enough memory to allocate new spectra\n");
            }
            for(j=0;j<=res->Spec[i].N;j++)
            {
              res->Spec[i].Re[j]= a->Spec[i].Re[j]/b;
              res->Spec[i].Im[j]= a->Spec[i].Im[j]/b;
            }
          }
          luaL_setmetatable(L, "Spectra_Type");
          return 1;
        }
        case LUA_TSTRING:
        {
          // (original variable) operator (string)
          return luaL_error(L, "Error: spectra divided by string\n");
        }
        case LUA_TUSERDATA:
        {
          // (original variable) operator (user data)
          if (luaL_testudata(L, len, "Complex_Type") != NULL)
          {
            // (original variable) operator (user data - complex)
            uint32_t i,j;
            Complex16 *b;
            SpectraTypeArray *res;
            b = (Complex16 *) luaL_checkudata(L, len, "Complex_Type");
            res = (SpectraTypeArray *) lua_newuserdata(L, sizeof(Complex16));
            res->NSpectra = a->NSpectra;
            res->Spec = malloc(res->NSpectra * sizeof(SpectraType));
            for(i=0;i<res->NSpectra;i++)
            {
              if(SpectrumAllocate(a->Spec[i], &(res->Spec[i])))
              {
                return luaL_error(L, "Error: not enough memory to allocate new spectra\n");
              }
              for(j=0;j<=res->Spec[i].N;j++)
              {
                res->Spec[i].Re[j]= (a->Spec[i].Re[j]*b->Re + a->Spec[i].Im[j]*b->Im)/(Sqr(b->Re) + Sqr(b->Im));
                res->Spec[i].Im[j]= (a->Spec[i].Im[j]*b->Re - a->Spec[i].Re[j]*b->Im)/(Sqr(b->Re) + Sqr(b->Im));
              }
            }
            luaL_setmetatable(L, "Spectra_Type");
            return 1;
          }
          if (luaL_testudata(L, len, "Wavefunction_Type") != NULL)
          {
            // (original variable) operator (user data - wavefunction)
            return luaL_error(L, "Error: spectra divided by wavefunction\n");
          }
          if (luaL_testudata(L, len, "Operator_Type") != NULL)
          {
            // (original variable) operator (user data - operator)
            return luaL_error(L, "Error: spectra divided by operator\n");
          }
          if (luaL_testudata(L, len, "Spectra_Type") != NULL)
          {
            // (original variable) operator (user data - spectra)
            uint32_t i,j;
            SpectraTypeArray *res, *b;
            b = (SpectraTypeArray *) luaL_checkudata(L, len, "Spectra_Type");
            if(a->NSpectra!=b->NSpectra) return luaL_error(L, "Error: Dividing to sets of spectra with different number of elements (%d) + (%d)\n",a->NSpectra,b->NSpectra);
            res = (SpectraTypeArray *) lua_newuserdata(L, sizeof(Complex16));
            res->NSpectra = a->NSpectra;
            res->Spec = malloc(res->NSpectra * sizeof(SpectraType));
            for(i=0;i<res->NSpectra;i++)
            {
              if((a->Spec[i].N==b->Spec[i].N)&&(a->Spec[i].Emin==b->Spec[i].Emin)&&(a->Spec[i].Emax==b->Spec[i].Emax))
              {
                if(SpectrumAllocate(b->Spec[i], &(res->Spec[i])))
                {
                  return luaL_error(L, "Error: not enough memory to allocate new spectra\n");
                }
                res->Spec[i].EminPole = Min(a->Spec[i].EminPole,b->Spec[i].EminPole);
                res->Spec[i].EmaxPole = Max(a->Spec[i].EmaxPole,b->Spec[i].EmaxPole);
                for(j=0;j<=res->Spec[i].N;j++)
                {
                  res->Spec[i].Re[j]= (a->Spec[i].Re[j] * b->Spec[i].Re[j] + a->Spec[i].Im[j] * b->Spec[i].Im[j])/(Sqr(b->Spec[i].Re[j])+Sqr(b->Spec[i].Im[j]));
                  res->Spec[i].Im[j]= (a->Spec[i].Im[j] * b->Spec[i].Re[j] - a->Spec[i].Re[j] * b->Spec[i].Im[j])/(Sqr(b->Spec[i].Re[j])+Sqr(b->Spec[i].Im[j]));
                }
              }
              else
              {
                return luaL_error(L, "Still to do: divide two spectra with different energy ranges, spectra %d\n",i);
              }
            }
            luaL_setmetatable(L, "Spectra_Type");
            return 1;
          }
          return luaL_error(L, "Error: spectra divided by unrecognized userdata\n");
        }
        case LUA_TTABLE:
        {
          lua_newtable(L);
          lua_pushnil(L);
          while(lua_next(L, len) != 0)
          {
            LuaSpectraDiv(L);
            lua_pushvalue(L, -3);
            lua_insert(L, -2);
            lua_settable(L, -5);
            lua_pop(L, 1);
          }
          return 1;
        }
        default:
        {
          return luaL_error(L, "Error: spectra divided by unrecognized type\n");
        }
      }
    }
  }
}

static int LuaSpectraPow(lua_State * L)
{
  uint32_t len = (uint32_t)lua_gettop(L);
  switch (lua_type(L,1))
  {
    case LUA_TNIL:
    {
      // (nil) operator (original variable)
      return luaL_error(L, "Error: nil to the power spectrum\n");
    }
    case LUA_TBOOLEAN:
    {
      // (boolean) operator (original variable)
      return luaL_error(L, "Error: boolean to the power spectrum\n");
    }
    case LUA_TLIGHTUSERDATA:
    {
      // (ligth user data) operator (original variable)
      return luaL_error(L, "Error: lightuserdata to the power spectrum\n");
    }
    case LUA_TNUMBER:
    {
      // (number) operator (original variable)
      return luaL_error(L, "Error: number to the power spectrum\n");
    }
    case LUA_TSTRING:
    {
      // (string) operator (original variable)
      return luaL_error(L, "Error: character string to the power spectrum\n");
    }
    case LUA_TTABLE:
    {
      lua_newtable(L);
      lua_pushnil(L);
      while(lua_next(L, 1) != 0)
      {
        lua_insert(L, 1);
        lua_pushvalue(L, -3);
        LuaSpectraPow(L);
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
      SpectraTypeArray *a = (SpectraTypeArray *) luaL_checkudata(L, 1, "Spectra_Type");
      switch (lua_type(L,len))
      {
        case LUA_TNIL:
        {
          // (original variable) operator (nil)
          return luaL_error(L, "Error: spectra to the power nil\n");
        }
        case LUA_TBOOLEAN:
        {
          // (original variable) operator (boolean)
          return luaL_error(L, "Error: spectra to the power boolean\n");
        }
        case LUA_TLIGHTUSERDATA:
        {
          // (original variable) operator (ligth user data)
          return luaL_error(L, "Error: spectra to the power lightuserdata\n");
        }
        case LUA_TNUMBER:
        {
          // (original variable) operator (number)
          uint32_t i,j;
          double b;
          double complex za, zb, zr;
          SpectraTypeArray *res;
          b = luaL_checknumber(L, len);
          res = (SpectraTypeArray *) lua_newuserdata(L, sizeof(Complex16));
          res->NSpectra = a->NSpectra;
          res->Spec = malloc(res->NSpectra * sizeof(SpectraType));
          for(i=0;i<res->NSpectra;i++)
          {
            if(SpectrumAllocate(a->Spec[i], &(res->Spec[i])))
            {
              return luaL_error(L, "Error: not enough memory to allocate new spectra\n");
            }
            for(j=0;j<=res->Spec[i].N;j++)
            {
              za = a->Spec[i].Re[j] + I * a->Spec[i].Im[j];
              zb = b;
              zr = cpow(za, zb);
              res->Spec[i].Re[j]= creal(zr);
              res->Spec[i].Im[j]= cimag(zr);
            }
          }
          luaL_setmetatable(L, "Spectra_Type");
          return 1;
        }
        case LUA_TSTRING:
        {
          // (original variable) operator (string)
          return luaL_error(L, "Error: spectra to the power string\n");
        }
        case LUA_TUSERDATA:
        {
          // (original variable) operator (user data)
          if (luaL_testudata(L, len, "Complex_Type") != NULL)
          {
            // (original variable) operator (user data - complex)
            uint32_t i,j;
            Complex16 *b;
            double complex za, zb, zr;
            SpectraTypeArray *res;
            b = (Complex16 *) luaL_checkudata(L, len, "Complex_Type");
            res = (SpectraTypeArray *) lua_newuserdata(L, sizeof(Complex16));
            res->NSpectra = a->NSpectra;
            res->Spec = malloc(res->NSpectra * sizeof(SpectraType));
            for(i=0;i<res->NSpectra;i++)
            {
              if(SpectrumAllocate(a->Spec[i], &(res->Spec[i])))
              {
                return luaL_error(L, "Error: not enough memory to allocate new spectra\n");
              }
              for(j=0;j<=res->Spec[i].N;j++)
              {
                za = a->Spec[i].Re[j] + I * a->Spec[i].Im[j];
                zb = b->Re + I * b->Im;
                zr = cpow(za, zb);
                res->Spec[i].Re[j]= creal(zr);
                res->Spec[i].Im[j]= cimag(zr);
              }
            }
            luaL_setmetatable(L, "Spectra_Type");
            return 1;
          }
          if (luaL_testudata(L, len, "Wavefunction_Type") != NULL)
          {
            // (original variable) operator (user data - wavefunction)
            return luaL_error(L, "Error: spectra to the power wavefunction\n");
          }
          if (luaL_testudata(L, len, "Operator_Type") != NULL)
          {
            // (original variable) operator (user data - operator)
            return luaL_error(L, "Error: spectra to the power operator\n");
          }
          if (luaL_testudata(L, len, "Spectra_Type") != NULL)
          {
            // (original variable) operator (user data - spectra)
            return luaL_error(L, "Error: spectra to the power spectra\n");
          }
          return luaL_error(L, "Error: spectra to the power unrecognized userdata\n");
        }
        case LUA_TTABLE:
        {
          lua_newtable(L);
          lua_pushnil(L);
          while(lua_next(L, len) != 0)
          {
            LuaSpectraPow(L);
            lua_pushvalue(L, -3);
            lua_insert(L, -2);
            lua_settable(L, -5);
            lua_pop(L, 1);
          }
          return 1;
        }
        default:
        {
          return luaL_error(L, "Error: spectra to the power unrecognized type\n");
        }
      }
    }
  }
}

static int LuaSpectraIndex(lua_State * L)
{
  SpectraTypeArray *a = (SpectraTypeArray*) luaL_checkudata(L, 1, "Spectra_Type");
  const char *str = lua_tostring(L, 2);
  currentudata = a;
  if (strcmp(str, "Print") == 0)
  {
    lua_pushcfunction(L, LuaSpectrumPrint);
    return 1;
  }
  if (strcmp(str, "Element") == 0)
  {
    lua_pushcfunction(L, LuaSpectrumElement);
    return 1;
  }
  if (strcmp(str, "Sum") == 0)
  {
    lua_pushcfunction(L, LuaSpectrumSum);
    return 1;
  }
  return luaL_error(L, "Error operation \".%s()\" not defined on spectra\n",str);
}

static int LuaSpectraToString(lua_State * L)
{
  SpectraTypeArray *spec = (SpectraTypeArray*) luaL_checkudata(L, 1, "Spectra_Type");
  lua_pushstring(L, spec->Spec[0].Name);
  return 1;
}

static int LuaSpectraDestroy(lua_State * L)
{
  uint32_t i;
  SpectraTypeArray *spec = (SpectraTypeArray*) luaL_checkudata(L, 1, "Spectra_Type");
  for(i=0;i<spec->NSpectra;i++)
  {
    SpectrumFree(&spec->Spec[i]);
  }
  free(spec->Spec);
  return 0;
}

static const struct luaL_Reg Spectra_meta[] = {
  {"__add", LuaSpectraAdd},
  {"__sub", LuaSpectraSub},
  {"__mul", LuaSpectraMul},
  {"__div", LuaSpectraDiv},
  {"__pow", LuaSpectraPow},
  {"__index", LuaSpectraIndex},
  {"__tostring", LuaSpectraToString},
  {"__gc", LuaSpectraDestroy},
  {NULL, NULL}
};