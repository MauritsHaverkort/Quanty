/*  test if a string is a valid lua variable name
 *
 */
uint32_t QStringIsValidLuaVariableName(const char * const c)
{
  uint32_t i;
  if(strlen(c)==0) return 1;
  if( c[0] >= '0' && c[0] <= '9' ) return 1;
  for(i=0;i<strlen(c);i++)
  {
    if(!(  ((c[i] >= '0' ) && (c[i] <= '9'))
         ||((c[i] >= 'A' ) && (c[i] <= 'Z'))
         ||((c[i] >= 'a' ) && (c[i] <= 'z'))
         ||( c[i] == '_' )
         )) return 1;
  }
  return 0;
}



/*
 *  gets the variable with name VariableName from 'file' FileName
 *  and adds it to the global definition in L
 *
 *  The file at FileName can be of different type, currently implemented
 *    - ASCII text (FileName is a directory)
 *
 *  returns 0 on succes, returns 1 on failure
 */
uint32_t LuaGet(lua_State *L, const char * const FileName, const char * const VariableName)
{
  struct stat info;
  // test if VariableName is lua safe
  if(QStringIsValidLuaVariableName(VariableName))
  {
    return luaL_error(L, "%s is not a valid token for a variable name in Quanty (only alpha numerical characters and '_' are allowed)\n",VariableName);
  }
  // determine type of FileName
  if( stat( FileName, &info ) != 0 )
  {
    // FileName does not seem to point to anything
    return luaL_error(L, "Tried to open file %s for reading, but file does not excists\n",FileName);
  }
  else if( info.st_mode & S_IFDIR )
  {
    // FileName points to a directory
    return LuaGetASCII(L, FileName, VariableName);
  }
  else
  {
    // FileName points to a file
    return luaL_error(L, "Tried to open file %s for reading with get, but filetype not recognized\n",FileName);
  }
  return 0;
}

/*
 *  puts the global variable with name VariableName
 *  and adds it to the file at FileName
 *
 *  The file at FileName can be of different type, currently implemented
 *    - ASCII text (FileName is a directory)
 *
 *  returns 0 on succes, returns 1 on failure
 */
uint32_t LuaPut(lua_State *L, const char * const FileName, const char * const VariableName, const char * const FileType)
{
  struct stat info;
  // test if VariableName is lua safe
  if(QStringIsValidLuaVariableName(VariableName))
  {
    return luaL_error(L, "%s is not a valid token for a variable name in Quanty (only alpha numerical characters and '_' are allowed)\n",VariableName);
  }
  // test if VariableName is defined and put on on the Lua stack
  lua_getglobal(L, VariableName);
  if(lua_isnil(L, -1))
  {
    luaL_warning(L, "try to save variable %s to file %s, but variable is undefined (nil)\n",VariableName,FileName);
    return 0;
  }
  // determine type of FileName
  if( stat( FileName, &info ) != 0 )
  {
    // FileName does not seem to point to anything
    // we need to create a new file of type FileType
    if((FileType==NULL)||(strcmp(FileType, "ASCII")==0))
    {
      // create a file of type ASCII
#if defined(_WIN32)
      mkdir(FileName);
#else
      mkdir(FileName, 0700);
#endif
      if(!(FileType==NULL))
      {
        lua_remove(L, 3);
      }
      return LuaPutASCII(L, FileName, VariableName);
    }
    if(strcmp(FileType, "hdf5")==0)
    {
      luaL_error(L, "hdf5 file format not finished\n");
      //TODO: create hdf5 file and put variable in it
    }
    luaL_error(L, "Tried to put variable to new file of unrecognized type %s\n",FileType);
  }
  else if( info.st_mode & S_IFDIR )
  {
    // FileName points to a directory -- type should be ASCII
    if(!(FileType==NULL))
    {
      if(!(strcmp(FileType, "ASCII")==0))
      {
        return luaL_error(L, "tried to put variable %s to file %s of indecated type %s, but file is of type ASCII on disk\n",VariableName,FileName,FileType);
      }
      lua_remove(L, 3);
    }
    return LuaPutASCII(L, FileName, VariableName);
  }
  else
  {
    // FileName points to a file
    return luaL_error(L, "Tried to open file %s for writing with put, but filetype not recognized\n",FileName);
    //TODO: recognize filetype and put variable in it
  }
  return 0;
}

/*  ============================================
 *  ============================================
 */

/* ASCII implementation of get and put
 *
 * database structure is given as folder with files
 * each variable name corresponds to file name
 * each file contains of excecutable Quanty code that creates the original variable back in memory
 *
 */


/*  Quanty variables contain lower case / upper case letters, nubmers or the character '_'
 *  In order to turn them into save system independent file names we change upper case letters to 'A' -> '^a'
 *  The maximum file length is given by MaxLength=200
 */
char * QuantyVariableToValidFilename(const char * str)
{
  const uint32_t MaxLength=200;
  uint32_t len,i,ofset=0,upper=0;
  char * retchar;
  len=(uint32_t)strlen(str);
  for(i=0;i<len;i++)
  {
    if((str[i]>='A')&&(str[i]<='Z')) upper++;
  }
  retchar = (char*)malloc(sizeof(char)*Min(MaxLength+1,1+len+upper));
  retchar[Min(MaxLength,len+upper)]='\0';
  for(i=0;i<Min(MaxLength,len+upper);i++)
  {
    if((str[i-ofset]>='A')&&(str[i-ofset]<='Z'))
    {
      retchar[i]='^';
      ofset++;
      i++;
      retchar[i]=str[i-ofset]+32;
    }
    else
    {
      retchar[i]=str[i-ofset];
    }
  }
  return retchar;
}

uint32_t LuaGetASCII(lua_State *L, const char * const DirName, const char * const VariableName)
{
  char *FileName;
  char *FullName;
  FileName = QuantyVariableToValidFilename(VariableName);
  FullName = (char*)malloc(sizeof(char)*(strlen(FileName)+strlen(DirName)+2+7));
  sprintf(FullName, "%s/%s.Quanty",DirName,FileName);
  lua_pop(L, 2);
  lua_pushstring(L, FullName);
  free(FileName);
  free(FullName);
  return luaB_dofile(L);
}

/*  Writes a double such that lua reads back the same double. Whole numbers are written without
 *  a decimal point, all others with the fewest of 15, 16 or 17 significant digits that give the
 *  number back (17 always do). Lua has no literal for inf and nan, they are written as expressions.
 */
void LuaASCIIPrintDouble(FILE * stream, double const n)
{
  char buf[32];
  int digits;
  if(isnan(n))
  {
    fprintf(stream,"(0/0)");
    return;
  }
  if(isinf(n))
  {
    fprintf(stream,(n>0)?"(1/0)":"(-1/0)");
    return;
  }
  if((n==floor(n))&&(n<INT32_MAX)&&(n>INT32_MIN))
  {
    fprintf(stream,"%d",(int32_t)n);
    return;
  }
  for(digits=15;digits<=17;digits++)
  {
    snprintf(buf, sizeof(buf), "%.*G", digits, n);
    if(strtod(buf, NULL)==n) break;
  }
  fputs(buf, stream);
}

/*  Writes the value at stack position narg as lua code that recreates it. Nothing is left on the
 *  stack; values are found by their absolute position, so everything pushed here is popped again
 *  before returning.
 *  Raw is true while the content of a ResponseFunction.New(...) is written. Tables in there that
 *  carry the ResponseFunctionMeta metatable (the val and con parts of a NaturalImpurityOrbital)
 *  are read back as plain tables and must not be wrapped in a ResponseFunction.New(...) of their own.
 */
static uint32_t LuaASCIIStreamPrintValue(lua_State *L, int narg, FILE * stream, bool const Raw)
{
  narg = lua_absindex(L, narg);
  if(lua_type(L,narg)==LUA_TNIL)
  {
    fprintf(stream, "nil");
    return 0;
  }
  if(lua_type(L,narg)==LUA_TNUMBER)
  {
    LuaASCIIPrintDouble(stream, lua_tonumber(L, narg));
    return 0;
  }
  if(lua_type(L,narg)==LUA_TBOOLEAN)
  {
    int n;
    n = lua_toboolean(L, narg);
    if(n==0)
    {
      fprintf(stream,"false");
    }
    else
    {
      fprintf(stream,"true");
    }
    return 0;
  }
  if(lua_type(L,narg)==LUA_TSTRING)
  {
    const char *n;
    uint32_t i;
    size_t len;
    n = lua_tolstring(L, narg, &len);
    fputc('\"', stream);
    for(i=0;i<len;i++)
    {
      unsigned char const c = (unsigned char) n[i];
      if(c=='\\')
      {
        fputc('\\', stream);
        fputc('\\', stream);
      }
      else if(c=='\"')
      {
        fputc('\\', stream);
        fputc('\"', stream);
      }
      else if(c=='\n')
      {
        fputc('\\', stream);
        fputc('n', stream);
      }
      else if((c<32)||(c==127))
      {
        // other control characters as a three digit decimal escape. With fewer digits the
        // characters that follow could be read as part of it ("\0" followed by "1" is "\01")
        fprintf(stream, "\\%03u", (unsigned int) c);
      }
      else
      {
        fputc(c, stream);
      }
    }
    fputc('\"', stream);
    return 0;
  }
  if(lua_type(L,narg)==LUA_TTABLE)
  {
    // check for specific meta data and then print differently
    if ((!Raw) && lua_getmetatable(L, narg))
    {
      lua_getglobal(L, "ResponseFunctionMeta");
      if(lua_rawequal(L, -1, -2)==1)
      {
        lua_pop(L, 2);
        goto ExportResponseFunction;
      }
      lua_pop(L, 2);
    }
    // print a table without meta data
    bool FirstElement=false;
    size_t i,tablelen = lua_objlen(L, narg);
    fprintf(stream, "{ ");
    /* array elements */
    for(i=0;i<tablelen;i++)
    {
      lua_rawgeti(L, narg, (int)i+1);                                    // put next element on stack
      if((i!=0)&&(lua_type(L, -1)==LUA_TTABLE)) fprintf(stream, "\n  "); // if next element is a table start on a new line
      LuaASCIIStreamPrintValue(L, lua_gettop(L), stream, Raw);           // turn next element into string
      lua_pop(L, 1);                                                     // remove element from stack
      if(i!=(tablelen-1)) fprintf(stream, " , ");                        // unless it is the last element print a ,
    }
    /* hash ellements */
    if(tablelen==0) FirstElement=true;
    lua_pushnil(L);
    while (lua_next(L,narg)!=0)
    {
      if((lua_type(L,-2)==LUA_TNUMBER) && (lua_tonumber(L, -2)==floor(lua_tonumber(L, -2))) && (0<lua_tonumber(L, -2)) && (lua_tonumber(L, -2)<=tablelen))
      {
        lua_pop(L, 1);                        // if key value is in the range [1,tablen] the value is already printed in the array part
      }
      else
      {
        if(FirstElement)
        {
          FirstElement=false;
        }
        else
        {
          fprintf(stream," ,\n  ");                      // print a , unless it is the first element
        }
        fprintf(stream,"[ ");
        LuaASCIIStreamPrintValue(L, lua_gettop(L)-1, stream, Raw); // print the key
        fprintf(stream," ] = ");
        LuaASCIIStreamPrintValue(L, lua_gettop(L), stream, Raw);   // print the value
        lua_pop(L, 1);                                   // pop the value
      }
    }
    fprintf(stream," }");                                // close the table
    return 0;
  }
  if(lua_type(L,narg)==LUA_TFUNCTION)
  {
    return luaL_error(L, "Error in LuaASCIIStreamPrint unrecognized type of function\n");
  }
  if(lua_type(L,narg)==LUA_TUSERDATA)
  {
    if (luaL_testudata(L, narg, "Complex_Type") != NULL)
    {
      Complex16 *n;
      n = (Complex16 *) luaL_checkudata(L, narg, "Complex_Type");
      if(Abs(n->Im)<10*DBL_MIN)
      {
        LuaASCIIPrintDouble(stream, n->Re);
      }
      else
      {
        fprintf(stream,"( ");
        LuaASCIIPrintDouble(stream, n->Re);
        if(n->Im>0)
        {
          fprintf(stream," + I * ");
        }
        else
        {
          fprintf(stream," - I * ");
        }
        LuaASCIIPrintDouble(stream, Abs(n->Im));
        fprintf(stream," )");
      }
      return 0;
    }
    if (luaL_testudata(L, narg, "Operator_Type") != NULL)
    {
      bool komma;
      int32_t i,j,l;
      OperatorType *a = (OperatorType *) luaL_checkudata(L, narg, "Operator_Type");
      fprintf(stream, "NewOperator( %d, %d,{\n",a->NFermion,a->NBoson);
      komma=false;
      for(l=0;l<=a->MaxLength;l++)
      {
        if(a->l[l].QComplex)
        {
          for(i=0;i<a->l[l].N;i++)
          {
            if(komma)
            {
              fprintf(stream,",\n {");
            }
            else
            {
              komma=true;
              fprintf(stream," {");
            }
            for(j=l-1;j>=0;j--)
            {
              fprintf(stream," %c%4u,",(((a->l[l].o[i * l + j]) & CrAnMask) ? ' ' : '-'),(a->l[l].o[i*l+j]&NotCrAnMask));
            }
            fprintf(stream,"%23.16E ",a->l[l].V[i].Re);
            if(a->l[l].V[i].Im<0)
            {
              fprintf(stream,"- I *");
            }
            else
            {
              fprintf(stream,"+ I *");
            }
            fprintf(stream,"%23.16E }",Abs(a->l[l].V[i].Im));
          }
        }
        else
        {
          for(i=0;i<a->l[l].N;i++)
          {
            if(komma)
            {
              fprintf(stream,",\n {");
            }
            else
            {
              komma=true;
              fprintf(stream," {");
            }
            for(j=l-1;j>=0;j--)
            {
              fprintf(stream,"%c%4u,",(((a->l[l].o[i * l + j]) & CrAnMask) ? ' ' : '-'),(a->l[l].o[i*l+j]&NotCrAnMask));
            }
            fprintf(stream,"%23.16E }",a->l[l].VRe[i]);
          }
        }
      }
      fprintf(stream, "},\n{{\"Name\",");
      lua_pushstring(L, a->Name);
      LuaASCIIStreamPrintValue(L, -1, stream, Raw);
      lua_pop(L, 1);
      fprintf(stream,"},{\"NBitsKey\",{");
      if(a->MaxLength>=1) fprintf(stream, "%d",a->l[1].NBitsKey);
      for(i=2;i<=a->MaxLength;i++)
      {
        fprintf(stream, ",%d",a->l[i].NBitsKey);
      }
      fprintf(stream, "}}");
      if(a->Res.N!=0)
      {
        fprintf(stream, ", {\"Restrictions\",{ %d, %d",a->NFermion,a->NBoson);
        for(i=0;i<a->Res.N;i++)
        {
          char *det;
          det = DetToString(a->NFermion, a->NBoson, a->Res.Determinant[i]);
          fprintf(stream,", {\"%s\",%u,%u}",det,a->Res.Min[i],a->Res.Max[i]);
          free(det);
        }
        fprintf(stream,"}}");
      }
      fprintf(stream,"})");
      return 0;
    }
    if (luaL_testudata(L, narg, "Wavefunction_Type") != NULL)
    {
      uint32_t i,j;
      WaveFunctionType *a = (WaveFunctionType *) luaL_checkudata(L, narg, "Wavefunction_Type");
      fprintf(stream, "NewWavefunction( %d, %d,{\n",a->NFermion,a->NBoson);
      if(a->QComplex)
      {
        for(i=1;i<a->N+1;i++)
        {
          fprintf(stream,"{\"");
          for(j=0;j<a->NFermion;j++)
          {
            fprintf(stream, "%c",(a->Determinant[BlokElementDet(i, a->NWord, (j >> 3))] & (128 >> (j & 7)))?'1':'0');
          }
          for(j=a->NWord-a->NBoson;j<a->NWord;j++)
          {
            fprintf(stream, " %2.2X", a->Determinant[BlokElementDet(i, a->NWord, j)]);
          }
          fprintf(stream, "\",%23.16E",a->ValueRe[BlokElement(i)]);
          if(a->ValueIm[BlokElement(i)]<0)
          {
            fprintf(stream," - I * ");
          }
          else
          {
            fprintf(stream," + I * ");
          }
          fprintf(stream, "%23.16E}",Abs(a->ValueIm[BlokElement(i)]));
          if(i!=a->N) fprintf(stream, ",\n");
        }
      }
      else
      {
        for(i=1;i<a->N+1;i++)
        {
          fprintf(stream,"{\"");
          for(j=0;j<a->NFermion;j++)
          {
            fprintf(stream, "%c",(a->Determinant[BlokElementDet(i, a->NWord, (j >> 3))] & (128 >> (j & 7)))?'1':'0');
          }
          for(j=a->NWord-a->NBoson;j<a->NWord;j++)
          {
            fprintf(stream, " %2.2X", a->Determinant[BlokElementDet(i, a->NWord, j)]);
          }
          fprintf(stream, "\",%23.16E}",a->ValueRe[BlokElement(i)]);
          if(i!=a->N) fprintf(stream, ",\n");
        }
      }
      fprintf(stream, "},\n{{\"Name\",");
      lua_pushstring(L, a->Name);
      LuaASCIIStreamPrintValue(L, -1, stream, Raw);
      lua_pop(L, 1);
      fprintf(stream,"},{\"NBitsKey\",%d}})",a->NBitsKey);
      return 0;
    }
    if (luaL_testudata(L, narg, "TightBinding_Type") != NULL)
    {
      //TODO: Tight binding. Name, Units, Cell and Atoms can be read from lua and written as a
      //      table, but reading Hopping gives a string and not the list of hoppings
      return luaL_error(L, "Error in io.put: variables of type TightBinding can not be saved (not implemented)\n");
    }
    if (luaL_testudata(L, narg, "Spectra_Type") != NULL)
    {
      SpectraTypeArray *n;
      uint32_t i,j;
      n = (SpectraTypeArray *) luaL_checkudata(L, narg, "Spectra_Type");
      fprintf(stream,"Spectra.FromTable({\n");
      for(i=0;i<n->NSpectra;i++)
      {
        if(i!=0)fprintf(stream,",\n");
        fprintf(stream,"{");
        for(j=0;j<=n->Spec[i].N;j++)
        {
          fprintf(stream,"{%23.16E,%23.16E %c I * %23.16E},",n->Spec[i].Emin+j*n->Spec[i].dE,n->Spec[i].Re[j],(n->Spec[i].Im[j]<0)?'-':'+',Abs(n->Spec[i].Im[j]));
        }
        if(strcmp("", n->Spec[i].Name)!=0)
        {
          fprintf(stream,"[\"Name\"]=\"%s\",",n->Spec[i].Name);
        }
        fprintf(stream,"[\"EminPole\"]=%23.16E,",n->Spec[i].EminPole);
        fprintf(stream,"[\"EmaxPole\"]=%23.16E,",n->Spec[i].EmaxPole);
        fprintf(stream,"[\"Gamma\"]=%23.16E,",n->Spec[i].Gamma);
        fprintf(stream,"[\"Integral\"]=%23.16E}",Sqr(n->Spec[i].b[0]));
      }
      fprintf(stream,"})");
      return 0;
    }
    if (luaL_testudata(L, narg, "ResponseFunction_Type") != NULL)
    {
      ExportResponseFunction:
      fprintf(stream,"ResponseFunction.New(");
      // turn the response function in a table, which ends up on top of the stack. It is called
      // with a copy of the value as its only argument, as ToTable reads position 1 and wants
      // nothing else on the stack; the value can be deep inside a table
      lua_pushcfunction(L, LuaResponseFunctionToTable);
      lua_pushvalue(L, narg);
      lua_call(L, 1, 1);
      lua_pushnil(L);
      lua_setmetatable(L, -2);
      LuaASCIIStreamPrintValue(L, lua_gettop(L), stream, true);
      lua_pop(L, 1);
      fprintf(stream,")\n");
      return 0;
    }
    if (luaL_testudata(L, narg, "Matrix_Type") != NULL)
    {
      fprintf(stream,"Matrix.New(");
      // turn the matrix in a table, on top of the stack, as for the response function above
      lua_pushcfunction(L, LuaMatrixToTable);
      lua_pushvalue(L, narg);
      lua_call(L, 1, 1);
      lua_pushnil(L);
      lua_setmetatable(L, -2);
      LuaASCIIStreamPrintValue(L, lua_gettop(L), stream, Raw);
      lua_pop(L, 1);
      fprintf(stream,")\n");
      return 0;
    }
    {
      // types that exist, but have no ASCII representation yet
      const char * const NotSaved[] = {"Graphics_Type", "BitMap_Type", "Tensor_Type", "GaussianRadialFunction_Type", "InterpolatingFunction_Type"};
      uint32_t i;
      for(i=0;i<sizeof(NotSaved)/sizeof(NotSaved[0]);i++)
      {
        if (luaL_testudata(L, narg, NotSaved[i]) != NULL)
        {
          return luaL_error(L, "Error in io.put: variables of type %s can not be saved (not implemented)\n", NotSaved[i]);
        }
      }
    }
    return luaL_error(L, "Error in LuaASCIIStreamPrint unrecognized type of user data\n");
  }
  if(lua_type(L,narg)==LUA_TTHREAD)
  {
    return luaL_error(L, "Error in LuaASCIIStreamPrint unrecognized type of thread\n");
  }
  if(lua_type(L,narg)==LUA_TLIGHTUSERDATA)
  {
    return luaL_error(L, "Error in LuaASCIIStreamPrint unrecognized type of light user data\n");
  }
  return luaL_error(L, "Error in LuaASCIIStreamPrint unrecognized type of data\n");
}

uint32_t LuaASCIIStreamPrint(lua_State *L, int narg, FILE * stream)
{
  return LuaASCIIStreamPrintValue(L, narg, stream, false);
}

uint32_t LuaPutASCII(lua_State *L, const char * const DirName, const char * const VariableName)
{
  FILE *fp;
  char *FileName;
  char *FullName;
  FileName = QuantyVariableToValidFilename(VariableName);
  FullName = (char*)malloc(sizeof(char)*(strlen(FileName)+strlen(DirName)+2+7));
  sprintf(FullName, "%s/%s.Quanty",DirName,FileName);
  free(FileName);
  fp=fopen(FullName, "w");
  if(fp==NULL)
  {
    printf("Open of file \"%s\" failed\n",FullName);
    printf("Could not save variable \"%s\"\n",VariableName);
    free(FullName);
    return 1;
  }
  free(FullName);
  fprintf(fp, "%s = ",VariableName);
  LuaASCIIStreamPrint(L,3,fp);
  fclose(fp);
  return 0;
}
