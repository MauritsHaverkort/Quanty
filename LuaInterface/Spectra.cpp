/*  copy a spectrum which is the last element on the stack
 *
 */
static int LuaSpectraCopy(lua_State * L)
{
  uint32_t len;
  len = (uint32_t)lua_gettop(L);
  if(len!=1) return luaL_error(L, "Error Spectra.Copy called with %d arguments, 1 expected\n",len);
  if(!((lua_type(L,1)==LUA_TUSERDATA)&&(luaL_testudata(L, len, "Spectra_Type") != NULL))) return luaL_error(L, "Error Spectra.Copy first element should be of spectra type\n");
  return Luacopy(L);
}

static int LuaSpectraBroaden(lua_State * L)
{
  uint32_t len,i,j, jList;
  MatrixType *FWHMList;
  SpectraTypeArray *a,*res;
  double GGaus,GLorentz,*GGausVar,*GLorentzVar;
  len = (uint32_t)lua_gettop(L);
  if(len!=3) return luaL_error(L, "Error Spectra.Broaden called with %d arguments, 3 arguments expected\n",len);
  if(!((lua_type(L,1)==LUA_TUSERDATA)&&(luaL_testudata(L, 1, "Spectra_Type") != NULL))) return luaL_error(L, "Error Spectra.Broaden first element should be of spectra type\n");
  a = (SpectraTypeArray *) luaL_checkudata(L, 1, "Spectra_Type");
  res = (SpectraTypeArray *) lua_newuserdata(L, sizeof(SpectraTypeArray));
  CopySpectraArray(a, res);
  luaL_setmetatable(L, "Spectra_Type");
  if(lua_type(L,2)==LUA_TNUMBER)
  {
    GGaus=luaL_checknumber(L, 2);
    for(i=0;i<res->NSpectra;i++)
    {
      BroadenG(&(res->Spec[i]), GGaus);
    }
  }
  else if (lua_type(L,2)==LUA_TTABLE)
  {
    FWHMList = luaL_checkMatrix(L, 2);
    if((FWHMList->QComplex)||(FWHMList->Nj!=2)) return luaL_error(L, "List of variable Gaussian broadening at position 2 should be of the form {{Energy_1, FWHM_1},...} with Energy and FWHM of type real\n");
    for(i=1;i<FWHMList->Ni;i++)
    {
      if(FWHMList->VRe[i-1][0]>FWHMList->VRe[i][0]) return luaL_error(L, "List of variable Gaussian broadening at position 2 should be of the form {{Energy_1, FWHM_1},{Energy_2, FWHM_2},...,{Energy_n, FWHM_n}},\nwith Energy_1<Energh_2<...<Energh_n");
    }
    GGausVar= (double*)malloc((res->Spec[0].N+1)*sizeof(double));
    for(i=0;i<res->NSpectra;i++)
    {
      GGausVar= (double*)realloc(GGausVar, (res->Spec[i].N+1)*sizeof(double));
      if(GGausVar==NULL) return luaL_error(L, "Fatal out of memory, intermediate allocation failed\n");
      jList=0;
      for(j=0;j<=res->Spec[i].N;j++)
      {
        while((jList<FWHMList->Ni)&&(res->Spec[i].Emin+j*res->Spec[i].dE>FWHMList->VRe[jList][0])) jList++;
        if(jList==0)
        {
          GGausVar[j]=FWHMList->VRe[jList][1];
        }
        else if(jList==FWHMList->Ni)
        {
          GGausVar[j]=FWHMList->VRe[jList-1][1];
        }
        else
        {
          GGausVar[j]=FWHMList->VRe[jList-1][1] + (FWHMList->VRe[jList][1]-FWHMList->VRe[jList-1][1]) * (res->Spec[i].Emin+j*res->Spec[i].dE-FWHMList->VRe[jList-1][0])/(FWHMList->VRe[jList][0]-FWHMList->VRe[jList-1][0]);
        }
      }
      BroadenGvar(&(res->Spec[i]), GGausVar);
    }
    free(GGausVar);
    FreeMatrix(FWHMList);
  }
  else return luaL_error(L, "Error Spectra.Broaden 2-nd argument should be either a number defining the Gaussian FWHM,\nor a table of the form {{Energy_1, FWHM_1},...}\n",len);
  if(lua_type(L,3)==LUA_TNUMBER)
  {
    GLorentz=luaL_checknumber(L, 3);
    for(i=0;i<res->NSpectra;i++)
    {
      BroadenL(&(res->Spec[i]), GLorentz);
    }
  }
  else if (lua_type(L,3)==LUA_TTABLE)
  {
    FWHMList = luaL_checkMatrix(L, 3);
    if((FWHMList->QComplex)||(FWHMList->Nj!=2)) return luaL_error(L, "List of variable Lorenzian broadening at position 3 should be of the form {{Energy_1, FWHM_1},...} with Energy and FWHM of type real\n");
    for(i=1;i<FWHMList->Ni;i++)
    {
      if(FWHMList->VRe[i-1][0]>FWHMList->VRe[i][0]) return luaL_error(L, "List of variable Lorenzian broadening at position 3 should be of the form {{Energy_1, FWHM_1},{Energy_2, FWHM_2},...,{Energy_n, FWHM_n}},\nwith Energy_1<Energh_2<...<Energh_n");
    }
    GLorentzVar= (double*)malloc((res->Spec[0].N+1)*sizeof(double));
    for(i=0;i<res->NSpectra;i++)
    {
      GLorentzVar= (double*)realloc(GLorentzVar, (res->Spec[i].N+1)*sizeof(double));
      if(GLorentzVar==NULL) return luaL_error(L, "Fatal out of memory, intermediate allocation failed\n");
      jList=0;
      for(j=0;j<=res->Spec[i].N;j++)
      {
        while((jList<FWHMList->Ni)&&(res->Spec[i].Emin+j*res->Spec[i].dE>FWHMList->VRe[jList][0])) jList++;
        if(jList==0)
        {
          GLorentzVar[j]=FWHMList->VRe[jList][1];
        }
        else if(jList==FWHMList->Ni)
        {
          GLorentzVar[j]=FWHMList->VRe[jList-1][1];
        }
        else
        {
          GLorentzVar[j]=FWHMList->VRe[jList-1][1] + (FWHMList->VRe[jList][1]-FWHMList->VRe[jList-1][1]) * (res->Spec[i].Emin+j*res->Spec[i].dE-FWHMList->VRe[jList-1][0])/(FWHMList->VRe[jList][0]-FWHMList->VRe[jList-1][0]);
        }
      }
      BroadenLvar(&(res->Spec[i]), GLorentzVar);
    }
    free(GLorentzVar);
    FreeMatrix(FWHMList);
  }
  else return luaL_error(L, "Error Spectra.Broaden 3-rd argument should be either a number defining the Lorenzian FWHM,\nor a table of the form {{Energy_1, FWHM_1},...}\n",len);
  return 1;
}

static int LuaSpectraBroadenMethod(lua_State * L)
{
  uint32_t len,i,j, jList;
  MatrixType *FWHMList;
  SpectraTypeArray *res;
  double GGaus,GLorentz,*GGausVar,*GLorentzVar;
  len = (uint32_t)lua_gettop(L);
  if(len!=2) return luaL_error(L, "Error Spectra.Broaden called with %d arguments, 2 arguments expected\n",len);
  res = (SpectraTypeArray*) currentudata;
  if(lua_type(L,1)==LUA_TNUMBER)
  {
    GGaus=luaL_checknumber(L, 1);
    for(i=0;i<res->NSpectra;i++)
    {
      BroadenG(&(res->Spec[i]), GGaus);
    }
  }
  else if (lua_type(L,1)==LUA_TTABLE)
  {
    FWHMList = luaL_checkMatrix(L, 1);
    if((FWHMList->QComplex)||(FWHMList->Nj!=2)) return luaL_error(L, "List of variable Gaussian broadening at position 1 should be of the form {{Energy_1, FWHM_1},...} with Energy and FWHM of type real\n");
    for(i=1;i<FWHMList->Ni;i++)
    {
      if(FWHMList->VRe[i-1][0]>FWHMList->VRe[i][0]) return luaL_error(L, "List of variable Gaussian broadening at position 1 should be of the form {{Energy_1, FWHM_1},{Energy_2, FWHM_2},...,{Energy_n, FWHM_n}},\nwith Energy_1<Energh_2<...<Energh_n");
    }
    GGausVar= (double*)malloc((res->Spec[0].N+1)*sizeof(double));
    for(i=0;i<res->NSpectra;i++)
    {
      GGausVar= (double*)realloc(GGausVar, (res->Spec[i].N+1)*sizeof(double));
      if(GGausVar==NULL) return luaL_error(L, "Fatal out of memory, intermediate allocation failed\n");
      jList=0;
      for(j=0;j<=res->Spec[i].N;j++)
      {
        while((jList<FWHMList->Ni)&&(res->Spec[i].Emin+j*res->Spec[i].dE>FWHMList->VRe[jList][0])) jList++;
        if(jList==0)
        {
          GGausVar[j]=FWHMList->VRe[jList][1];
        }
        else if(jList==FWHMList->Ni)
        {
          GGausVar[j]=FWHMList->VRe[jList-1][1];
        }
        else
        {
          GGausVar[j]=FWHMList->VRe[jList-1][1] + (FWHMList->VRe[jList][1]-FWHMList->VRe[jList-1][1]) * (res->Spec[i].Emin+j*res->Spec[i].dE-FWHMList->VRe[jList-1][0])/(FWHMList->VRe[jList][0]-FWHMList->VRe[jList-1][0]);
        }
      }
      BroadenGvar(&(res->Spec[i]), GGausVar);
    }
    free(GGausVar);
    FreeMatrix(FWHMList);
  }
  else return luaL_error(L, "Error Spectra.Broaden 1-st argument should be either a number defining the Gaussian FWHM,\nor a table of the form {{Energy_1, FWHM_1},...}\n",len);
  if(lua_type(L,2)==LUA_TNUMBER)
  {
    GLorentz=luaL_checknumber(L, 2);
    for(i=0;i<res->NSpectra;i++)
    {
      BroadenL(&(res->Spec[i]), GLorentz);
    }
  }
  else if (lua_type(L,2)==LUA_TTABLE)
  {
    FWHMList = luaL_checkMatrix(L, 2);
    if((FWHMList->QComplex)||(FWHMList->Nj!=2)) return luaL_error(L, "List of variable Lorenzian broadening at position 2 should be of the form {{Energy_1, FWHM_1},...} with Energy and FWHM of type real\n");
    for(i=1;i<FWHMList->Ni;i++)
    {
      if(FWHMList->VRe[i-1][0]>FWHMList->VRe[i][0]) return luaL_error(L, "List of variable Lorenzian broadening at position 2 should be of the form {{Energy_1, FWHM_1},{Energy_2, FWHM_2},...,{Energy_n, FWHM_n}},\nwith Energy_1<Energh_2<...<Energh_n");
    }
    GLorentzVar= (double*)malloc((res->Spec[0].N+1)*sizeof(double));
    for(i=0;i<res->NSpectra;i++)
    {
      GLorentzVar= (double*)realloc(GLorentzVar, (res->Spec[i].N+1)*sizeof(double));
      if(GLorentzVar==NULL) return luaL_error(L, "Fatal out of memory, intermediate allocation failed\n");
      jList=0;
      for(j=0;j<=res->Spec[i].N;j++)
      {
        while((jList<FWHMList->Ni)&&(res->Spec[i].Emin+j*res->Spec[i].dE>FWHMList->VRe[jList][0])) jList++;
        if(jList==0)
        {
          GLorentzVar[j]=FWHMList->VRe[jList][1];
        }
        else if(jList==FWHMList->Ni)
        {
          GLorentzVar[j]=FWHMList->VRe[jList-1][1];
        }
        else
        {
          GLorentzVar[j]=FWHMList->VRe[jList-1][1] + (FWHMList->VRe[jList][1]-FWHMList->VRe[jList-1][1]) * (res->Spec[i].Emin+j*res->Spec[i].dE-FWHMList->VRe[jList-1][0])/(FWHMList->VRe[jList][0]-FWHMList->VRe[jList-1][0]);
        }
      }
      BroadenLvar(&(res->Spec[i]), GLorentzVar);
    }
    free(GLorentzVar);
    FreeMatrix(FWHMList);
  }
  else return luaL_error(L, "Error Spectra.Broaden 2-nd argument should be either a number defining the Lorenzian FWHM,\nor a table of the form {{Energy_1, FWHM_1},...}\n",len);
  return 0;
}

static int LuaSpectraShift(lua_State * L)
{
  uint32_t len,i;
  SpectraTypeArray *a,*res;
  double shift;
  len = (uint32_t)lua_gettop(L);
  if(len!=2) return luaL_error(L, "Error Spectra.Shift called with %d arguments, 2 arguments expected\n",len);
  if(!((lua_type(L,1)==LUA_TUSERDATA)&&(luaL_testudata(L, 1, "Spectra_Type") != NULL))) return luaL_error(L, "Error Spectra.Shift first element should be of spectra type\n");
  a = (SpectraTypeArray *) luaL_checkudata(L, 1, "Spectra_Type");
  if(lua_type(L,2)==LUA_TNUMBER)
  {
    shift=luaL_checknumber(L, 2);
    res = (SpectraTypeArray *) lua_newuserdata(L, sizeof(SpectraTypeArray));
    CopySpectraArray(a, res);
    for(i=0;i<res->NSpectra;i++)
    {
      res->Spec[i].Emax+=shift;
      res->Spec[i].Emin+=shift;
      res->Spec[i].EmaxPole+=shift;
      res->Spec[i].EminPole+=shift;
      res->Spec[i].FirstMoment+=shift;
      res->Spec[i].a[0]+=shift;
      res->Spec[i].a[1]+=shift;
      res->Spec[i].a[2]+=shift;
      res->Spec[i].a[3]+=shift;
      res->Spec[i].a[4]+=shift;
    }
    luaL_setmetatable(L, "Spectra_Type");
    return 1;
  }
  else if(lua_type(L,2)==LUA_TTABLE)
  {
    if(a->NSpectra!=lua_objlen(L,2)) return luaL_error(L, "Error Spectra.Shift Table with shifts is of length %d, but there are %d spectra\n",lua_objlen(L,2),a->NSpectra);
    res = (SpectraTypeArray *) lua_newuserdata(L, sizeof(SpectraTypeArray));
    CopySpectraArray(a, res);
    for(i=0;i<res->NSpectra;i++)
    {
      lua_rawgeti(L, 2, i+1);
      shift=luaL_checknumber(L, -1);
      lua_pop(L, 1);
      res->Spec[i].Emax+=shift;
      res->Spec[i].Emin+=shift;
      res->Spec[i].EmaxPole+=shift;
      res->Spec[i].EminPole+=shift;
      res->Spec[i].FirstMoment+=shift;
      res->Spec[i].a[0]+=shift;
      res->Spec[i].a[1]+=shift;
      res->Spec[i].a[2]+=shift;
      res->Spec[i].a[3]+=shift;
      res->Spec[i].a[4]+=shift;
    }
    luaL_setmetatable(L, "Spectra_Type");
    return 1;
  }
  return luaL_error(L, "Error Spectra.Shift second argument should be a number or a table of numbers\n");
}

static int LuaSpectraShiftMethod(lua_State * L)
{
  uint32_t len,i;
  SpectraTypeArray *res;
  double shift;
  len = (uint32_t)lua_gettop(L);
  if(len!=1) return luaL_error(L, "Error Spectra.Shift called with %d arguments, 1 argument expected\n",len);
  res = (SpectraTypeArray*) currentudata;
  if(lua_type(L,1)==LUA_TNUMBER)
  {
    shift=luaL_checknumber(L, 1);
    for(i=0;i<res->NSpectra;i++)
    {
      res->Spec[i].Emax+=shift;
      res->Spec[i].Emin+=shift;
      res->Spec[i].EmaxPole+=shift;
      res->Spec[i].EminPole+=shift;
      res->Spec[i].FirstMoment+=shift;
      res->Spec[i].a[0]+=shift;
      res->Spec[i].a[1]+=shift;
      res->Spec[i].a[2]+=shift;
      res->Spec[i].a[3]+=shift;
      res->Spec[i].a[4]+=shift;
    }
    return 0;
  }
  else if(lua_type(L,1)==LUA_TTABLE)
  {
    if(res->NSpectra!=lua_objlen(L,1)) return luaL_error(L, "Error Spectra.Shift Table with shifts is of length %d, but there are %d spectra\n",lua_objlen(L,1),res->NSpectra);
    for(i=0;i<res->NSpectra;i++)
    {
      lua_rawgeti(L, 1, i+1);
      shift=luaL_checknumber(L, -1);
      lua_pop(L, 1);
      res->Spec[i].Emax+=shift;
      res->Spec[i].Emin+=shift;
      res->Spec[i].EmaxPole+=shift;
      res->Spec[i].EminPole+=shift;
      res->Spec[i].FirstMoment+=shift;
      res->Spec[i].a[0]+=shift;
      res->Spec[i].a[1]+=shift;
      res->Spec[i].a[2]+=shift;
      res->Spec[i].a[3]+=shift;
      res->Spec[i].a[4]+=shift;
    }
    return 0;
  }
  return luaL_error(L, "Error Spectra.Shift second argument should be a number or a table of numbers\n");
}

static int LuaSpectraToTable(lua_State * L)
{
  uint32_t len;
  SpectraTypeArray *a;
  len = (uint32_t)lua_gettop(L);
  if(len!=1) return luaL_error(L, "Error Spectra.ToTable called with %d arguments, 1 argument expected\n",len);
  if(!((lua_type(L,1)==LUA_TUSERDATA)&&(luaL_testudata(L, 1, "Spectra_Type") != NULL))) return luaL_error(L, "Error Spectra.ToTable first element should be of spectra type\n");
  a = (SpectraTypeArray *) luaL_checkudata(L, 1, "Spectra_Type");
  LuaPushTable_SpectraTypeArray(L, a);
  return 1;
}


static int LuaTableToSpectra(lua_State * L)
{
  SpectraTypeArray *res;
  Complex16 * tempval;
  uint32_t i,j;
  double Enow,epsilon;
  if(lua_gettop(L)!=1) return luaL_error(L, "Error Spectra.FromTable called with %d arguments, 1 argument (a table) expected\n",lua_gettop(L));
  luaL_checktype(L, 1, LUA_TTABLE);
  res = (SpectraTypeArray *) lua_newuserdata(L, sizeof(SpectraTypeArray));
  luaL_setmetatable(L, "Spectra_Type");
  res->NSpectra = 0;
  res->Spec= (SpectraType*)malloc(lua_objlen(L, 1)*sizeof(SpectraType));
  for(i=0;i<(uint32_t)lua_objlen(L, 1);i++)
  {
    lua_rawgeti(L, 1, i+1);
    if(lua_type(L,-1)!=LUA_TTABLE) return luaL_error(L, "Error Spectra.FromTable the %d-th spectrum is not a table\n",i+1);
    res->Spec[i].N=(uint32_t)lua_objlen(L, -1);
    if(res->Spec[i].N==0) return luaL_error(L, "Error Spectra.FromTable the %d-th spectrum has no elements\n",i+1);
    res->Spec[i].N--;
    res->Spec[i].Re= (double*)malloc((res->Spec[i].N+1)*sizeof(double));
    res->Spec[i].Im= (double*)malloc((res->Spec[i].N+1)*sizeof(double));
    res->NSpectra++;
    // get Emin
    lua_rawgeti(L, -1, 1);
    if(lua_type(L,-1)!=LUA_TTABLE) return luaL_error(L, "Error Spectra.FromTable the %d-th spectrum should be a table of elements consisting of {energy,Intensity}\n",i+1);
    lua_rawgeti(L, -1, 1);
    res->Spec[i].Emin=(double)lua_tonumber(L, -1);
    lua_pop(L, 2);
    // get Emax
    lua_rawgeti(L, -1, res->Spec[i].N+1);
    if(lua_type(L,-1)!=LUA_TTABLE) return luaL_error(L, "Error Spectra.FromTable the %d-th spectrum should be a table of elements consisting of {energy,Intensity}\n",i+1);
    lua_rawgeti(L, -1, 1);
    res->Spec[i].Emax=(double)lua_tonumber(L, -1);
    lua_pop(L, 2);
    // set meta information
    res->Spec[i].dE=(res->Spec[i].Emax-res->Spec[i].Emin)/res->Spec[i].N;
    lua_getfield(L, -1, "Name");
    if(lua_type(L, -1)==LUA_TNIL)
    {
      lua_pop(L, 1);
      strcpy(res->Spec[i].Name, "");
    }
    else if(lua_type(L,-1)==LUA_TSTRING)
    {
      strncpy(res->Spec[i].Name, lua_tostring(L, -1), MAXStrLength);
      lua_pop(L,1);
    }
    else return luaL_error(L, "Error Spectra.FromTable \"Name\" of spectrum %d should be a string\n",lua_gettop(L));
    lua_getfield(L, -1,"ResonantEnergy");
    if(lua_type(L, -1)==LUA_TNIL)
    {
      lua_pop(L,1);
      res->Spec[i].ResonantEnergy=NAN;
    }
    else if(lua_type(L,-1)==LUA_TNUMBER)
    {
      res->Spec[i].ResonantEnergy=(double)lua_tonumber(L, -1);
      lua_pop(L,1);
    }
    else return luaL_error(L, "Error Spectra.FromTable \"ResonantEnergy\" of spectrum %d should be a number\n",lua_gettop(L));
    lua_getfield(L, -1, "EminPole");
    if(lua_type(L, -1)==LUA_TNIL)
    {
      lua_pop(L, 1);
      res->Spec[i].EminPole=res->Spec[i].Emin;
    }
    else if(lua_type(L,-1)==LUA_TNUMBER)
    {
      res->Spec[i].EminPole=(double)lua_tonumber(L, -1);
      lua_pop(L,1);
    }
    else return luaL_error(L, "Error Spectra.FromTable \"EminPole\" of spectrum %d should be a number\n",lua_gettop(L));
    lua_getfield(L, -1, "EmaxPole");
    if(lua_type(L, -1)==LUA_TNIL)
    {
      lua_pop(L, 1);
      res->Spec[i].EmaxPole=res->Spec[i].Emax;
    }
    else if(lua_type(L,-1)==LUA_TNUMBER)
    {
      res->Spec[i].EmaxPole=(double)lua_tonumber(L, -1);
      lua_pop(L,1);
    }
    else return luaL_error(L, "Error Spectra.FromTable \"EmaxPole\" of spectrum %d should be a number\n",lua_gettop(L));
    lua_getfield(L, -1, "Gamma");
    if(lua_type(L, -1)==LUA_TNIL)
    {
      lua_pop(L, 1);
      res->Spec[i].Gamma=0;
    }
    else if(lua_type(L,-1)==LUA_TNUMBER)
    {
      res->Spec[i].Gamma=(double)lua_tonumber(L, -1);
      lua_pop(L,1);
    }
    else return luaL_error(L, "Error Spectra.FromTable \"Gamma\" of spectrum %d should be a number\n",lua_gettop(L));
    // read the intensity
    epsilon = Abs(res->Spec[i].Emax-res->Spec[i].Emin)*DBL_EPSILON*res->Spec[i].N;
    Enow=res->Spec[i].Emin;
    for(j=0;j<=res->Spec[i].N;j++)
    {
      lua_rawgeti(L, -1, j+1);
      if(lua_type(L,-1)!=LUA_TTABLE) return luaL_error(L, "Error Spectra.FromTable the %d-th spectrum should be a table of elements consisting of {energy,Intensity}\n",i+1);
      lua_rawgeti(L, -1, 2);
      if (lua_isnumber(L, -1))
      {
        res->Spec[i].Re[j]=lua_tonumber(L, -1);
        res->Spec[i].Im[j]=0;
      }
      else
      {
        tempval = (Complex16 *) luaL_checkudata(L, -1, "Complex_Type");
        res->Spec[i].Re[j] = tempval->Re;
        res->Spec[i].Im[j] = tempval->Im;
      }
      lua_rawgeti(L, -2, 1);
      if(Abs((double)lua_tonumber(L, -1)-Enow)>epsilon) return luaL_error(L, "Error Spectra.FromTable the %d-th spectrum should be a table of elements consisting of {energy,Intensity}, whereby the energy needs to be a list of increasing energies with constant difference.\n(Energy steps not constant.)",i+1);
      lua_pop(L, 3);
      Enow+=res->Spec[i].dE;
    }
    lua_getfield(L, -1, "Integral");
    if(lua_type(L, -1)==LUA_TNIL)
    {
      lua_pop(L, 1);
      res->Spec[i].b[0]=0;
      for(j=0;j<=res->Spec[i].N;j++)
      {
        res->Spec[i].b[0]-=res->Spec[i].Im[j];
      }
      res->Spec[i].b[0]*=res->Spec[i].dE/Pi;
      res->Spec[i].b[0]=sqrt(res->Spec[i].b[0]);
    }
    else if(lua_type(L,-1)==LUA_TNUMBER)
    {
      res->Spec[i].b[0]=sqrt((double)lua_tonumber(L, -1));
      lua_pop(L,1);
    }
    else return luaL_error(L, "Error Spectra.FromTable \"Integral\" of spectrum %d should be a number\n",lua_gettop(L));
    lua_pop(L,1);
  }
  return 1;
}

static int LuaSpectraMultiplyByEnergy(lua_State * L)
{
  uint32_t len,i,j;
  SpectraTypeArray *a,*res;
  len = (uint32_t)lua_gettop(L);
  if(len!=1) return luaL_error(L, "Error Spectra.MultiplyByEnergy called with %d arguments, 1 argument expected\n",len);
  if(!((lua_type(L,1)==LUA_TUSERDATA)&&(luaL_testudata(L, 1, "Spectra_Type") != NULL))) return luaL_error(L, "Error Spectra.MultiplyByEnergy first element should be of spectra type\n");
  a = (SpectraTypeArray *) luaL_checkudata(L, 1, "Spectra_Type");
  res = (SpectraTypeArray *) lua_newuserdata(L, sizeof(SpectraTypeArray));
  CopySpectraArray(a, res);
  for(i=0;i<res->NSpectra;i++)
  {
    for(j=0;j<=res->Spec[i].N;j++)
    {
      res->Spec[i].Re[j]*=(res->Spec[i].Emin+j*res->Spec[i].dE);
      res->Spec[i].Im[j]*=(res->Spec[i].Emin+j*res->Spec[i].dE);
    }
  }
  luaL_setmetatable(L, "Spectra_Type");
  return 1;
}

static int LuaSpectraMultiplyByEnergyMethod(lua_State * L)
{
  uint32_t len,i,j;
  SpectraTypeArray *res = (SpectraTypeArray*) currentudata;
  len = (uint32_t)lua_gettop(L);
  if(len!=0) return luaL_error(L, "Error Spectra.MultiplyByEnergy called with %d arguments, 0 arguments expected\n",len);
  for(i=0;i<res->NSpectra;i++)
  {
    for(j=0;j<=res->Spec[i].N;j++)
    {
      res->Spec[i].Re[j]*=(res->Spec[i].Emin+j*res->Spec[i].dE);
      res->Spec[i].Im[j]*=(res->Spec[i].Emin+j*res->Spec[i].dE);
    }
  }
  return 0;
}

static int LuaSpectraSum(lua_State * L)
{
  uint32_t len,i,j,k;
  SpectraTypeArray *a,*res;
  Complex16 alpha;
  len = (uint32_t)lua_gettop(L);
  if(len<2) return luaL_error(L, "Error Spectra.Sum called with %d arguments, 2 or more arguments expected\n",len);
  if(!((lua_type(L,1)==LUA_TUSERDATA)&&(luaL_testudata(L, 1, "Spectra_Type") != NULL))) return luaL_error(L, "Error Spectra.Sum first element should be of spectra type\n");
  a = (SpectraTypeArray *) luaL_checkudata(L, 1, "Spectra_Type");
  res = (SpectraTypeArray *) lua_newuserdata(L, sizeof(Complex16));
  luaL_setmetatable(L, "Spectra_Type");
  res->NSpectra = len-1;
  res->Spec= (SpectraType*)malloc(res->NSpectra*sizeof(SpectraType));
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
  strcpy(res->Spec[0].Name, "Spectra");
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
    if(lua_type(L, j+2)!=LUA_TTABLE)
    {
      return luaL_error(L, "Error in Spectra.Sum(Spectra,arg_1,arg_2,...) arg_%d should be a table\n",j+1);
    }
    if(((uint32_t) lua_objlen(L, j+2))!=a->NSpectra)
    {
      return luaL_error(L, "Error in Spectra.Sum(Spectra,arg_1,arg_2,...) arg_%d should be a table of length %d\n",j+1,a->NSpectra);
    }
    if(SpectrumAllocate(res->Spec[0], &(res->Spec[j])))
    {
      return luaL_error(L, "Error: not enough memory to allocate new spectra\n");
    }
    for(k=0;k<a->NSpectra;k++)
    {
      lua_rawgeti(L, j+2, k+1);
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
        return luaL_error(L, "Error in Spectra.Sum(Spectra,arg_1,arg_2,...) elemnt %d in arg_%d should be a real or complex number\n",k+1,j+1);
      }
      lua_pop(L,1);
      SafeAddMultipliedSpectra(&(res->Spec[j]), a->Spec[k], alpha.Re, alpha.Im);
    }
  }
  return 1;
}

/*  Return a sum of the spectra
 *  on stack a SpectrumTypeArray and a list of numbers (doubles or complex) with the length of the list equal to the number of elements in the spectra
 *  list of lists are allowed
 */
static int LuaSpectraSumMethod(lua_State * L)
{
  uint32_t i,j,k;
  SpectraTypeArray *res;
  Complex16 alpha;
  SpectraTypeArray *a = (SpectraTypeArray *)currentudata;
  uint32_t len = (uint32_t)lua_gettop(L);
  res = (SpectraTypeArray *) lua_newuserdata(L, sizeof(Complex16));
  luaL_setmetatable(L, "Spectra_Type");
  res->NSpectra = len;
  res->Spec= (SpectraType*)malloc(res->NSpectra*sizeof(SpectraType));
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

static int LuaSpectraElement(lua_State * L)
{
  SpectraTypeArray *a;
  uint32_t len = (uint32_t)lua_gettop(L);
  if(len!=2) return luaL_error(L, "Error Spectra.Element called with %d arguments, 2 arguments expected\n",len);
  if(!((lua_type(L,1)==LUA_TUSERDATA)&&(luaL_testudata(L, 1, "Spectra_Type") != NULL))) return luaL_error(L, "Error Spectra.Element first argument should be of spectra type\n");
  a = (SpectraTypeArray *) luaL_checkudata(L, 1, "Spectra_Type");
  switch (lua_type(L,2))
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
      res->Spec= (SpectraType*)malloc(res->NSpectra*sizeof(SpectraType));
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
      res->Spec= (SpectraType*)malloc(res->NSpectra*sizeof(SpectraType));
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

/*  Return a subset of the spectra
 *  on stack a SpectrumTypeArray and either a number or a list of numbers (integers)
 */
static int LuaSpectraElementMethod(lua_State * L)
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
        res->Spec= (SpectraType*)malloc(res->NSpectra*sizeof(SpectraType));
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
        res->Spec= (SpectraType*)malloc(res->NSpectra*sizeof(SpectraType));
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

static int LuaSpectraToFunction(lua_State * L)
{
  return 1;
}

/*  Print a set of spectra
 *
 */
static int LuaSpectraPrintMethod(lua_State * L)
{
  struct
  {
    char *stream=NULL;
  } Options;
  uint32_t len = (uint32_t)lua_gettop(L);
  LuaOptionsType LuaSpectraPrintOptions[]={
    {"file","char*",&(Options.stream)},
    {NULL, NULL, NULL}};
  if(!((len==0)||(len==1)))
  {
    return luaL_error(L, "Error in .Print(arg) caled with %d arguments 0 or 1 expected\n",len);
  }
  if(len==1)
  {
    LuaCheckOptions(L, -1, LuaSpectraPrintOptions);
  }
  if(Options.stream==NULL)
  {
    PrintSpectra(((SpectraTypeArray *)currentudata)->Spec, ((SpectraTypeArray *)currentudata)->NSpectra);
  }
  else
  {
    FilePrintSpectra(((SpectraTypeArray *)currentudata)->Spec, ((SpectraTypeArray *)currentudata)->NSpectra, Options.stream, "w");
    VerbosityPrintf(0x2,"Spectra printed to file: %s\n",Options.stream);
  }
  return 0;
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
      res = (SpectraTypeArray *) lua_newuserdata(L, sizeof(SpectraTypeArray));
      if(CopySpectraArray(b, res)) return luaL_error(L, "Error: not enough memory to allocate new spectra\n");
      for(i=0;i<res->NSpectra;i++)
      {
        for(j=0;j<=res->Spec[i].N;j++)
        {
          res->Spec[i].Re[j]+=a;
        }
        res->Spec[i].a[0]+=a;
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
          res = (SpectraTypeArray *) lua_newuserdata(L, sizeof(SpectraTypeArray));
          if(CopySpectraArray(a, res)) return luaL_error(L, "Error: not enough memory to allocate new spectra\n");
          for(i=0;i<res->NSpectra;i++)
          {
            for(j=0;j<=res->Spec[i].N;j++)
            {
              res->Spec[i].Re[j]+=b;
            }
            res->Spec[i].a[0]+=b;
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
            res = (SpectraTypeArray *) lua_newuserdata(L, sizeof(SpectraTypeArray));
            if(CopySpectraArray(a, res)) return luaL_error(L, "Error: not enough memory to allocate new spectra\n");
            for(i=0;i<res->NSpectra;i++)
            {
              for(j=0;j<=res->Spec[i].N;j++)
              {
                res->Spec[i].Re[j]+=b->Re;
                res->Spec[i].Re[j]+=b->Im;
              }
              //TODO: recalculate phase, a and b
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
            uint32_t i;
            SpectraTypeArray *res, *b;
            b = (SpectraTypeArray *) luaL_checkudata(L, len, "Spectra_Type");
            if(a->NSpectra!=b->NSpectra) return luaL_error(L, "Error: Adding two sets of spectra with different number of elements (%d) + (%d)\n",a->NSpectra,b->NSpectra);
            res = (SpectraTypeArray *) lua_newuserdata(L, sizeof(SpectraTypeArray));
            if(CopySpectraArray(a, res)) return luaL_error(L, "Error: not enough memory to allocate new spectra\n");
            for(i=0;i<res->NSpectra;i++)
            {
              SafeAddMultipliedSpectra(&(res->Spec[i]), b->Spec[i], 1, 0);
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
      res = (SpectraTypeArray *) lua_newuserdata(L, sizeof(SpectraTypeArray));
      if(CopySpectraArray(b, res)) return luaL_error(L, "Error: not enough memory to allocate new spectra\n");
      for(i=0;i<res->NSpectra;i++)
      {
        for(j=0;j<=res->Spec[i].N;j++)
        {
          res->Spec[i].Re[j]*=-1;
          res->Spec[i].Im[j]*=-1;
          res->Spec[i].Re[j]+=a;
        }
        res->Spec[i].a[0]-=a;
        res->Spec[i].phase+=Pi;
        if(res->Spec[i].phase>2*Pi) res->Spec[i].phase-=2*Pi;
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
          res = (SpectraTypeArray *) lua_newuserdata(L, sizeof(SpectraTypeArray));
          if(CopySpectraArray(a, res)) return luaL_error(L, "Error: not enough memory to allocate new spectra\n");
          for(i=0;i<res->NSpectra;i++)
          {
            for(j=0;j<=res->Spec[i].N;j++)
            {
              res->Spec[i].Re[j]-=b;
            }
            res->Spec[i].a[0]-=b;
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
            res = (SpectraTypeArray *) lua_newuserdata(L, sizeof(SpectraTypeArray));
            if(CopySpectraArray(a, res)) return luaL_error(L, "Error: not enough memory to allocate new spectra\n");
            for(i=0;i<res->NSpectra;i++)
            {
              for(j=0;j<=res->Spec[i].N;j++)
              {
                res->Spec[i].Re[j]-=b->Re;
                res->Spec[i].Re[j]-=b->Im;
              }
              //TODO: recalculate phase, a and b
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
            uint32_t i;
            SpectraTypeArray *res, *b;
            b = (SpectraTypeArray *) luaL_checkudata(L, len, "Spectra_Type");
            if(a->NSpectra!=b->NSpectra) return luaL_error(L, "Error: Subtracting to sets of spectra with different number of elements (%d) + (%d)\n",a->NSpectra,b->NSpectra);
            res = (SpectraTypeArray *) lua_newuserdata(L, sizeof(SpectraTypeArray));
            if(CopySpectraArray(a, res)) return luaL_error(L, "Error: not enough memory to allocate new spectra\n");
            for(i=0;i<res->NSpectra;i++)
            {
              SafeAddMultipliedSpectra(&(res->Spec[i]), b->Spec[i], -1, 0);
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
      uint32_t i;
      double a;
      SpectraTypeArray *res, *b;
      a = luaL_checknumber(L, 1);
      b = (SpectraTypeArray *) luaL_checkudata(L, len, "Spectra_Type");
      res = (SpectraTypeArray *) lua_newuserdata(L, sizeof(SpectraTypeArray));
      if(CopySpectraArray(b, res)) return luaL_error(L, "Error: not enough memory to allocate new spectra\n");
      for(i=0;i<res->NSpectra;i++)
      {
        SpectraScalarMultiply(&(res->Spec[i]), a, 0);
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
          uint32_t i;
          double b;
          SpectraTypeArray *res;
          b = luaL_checknumber(L, len);
          res = (SpectraTypeArray *) lua_newuserdata(L, sizeof(SpectraTypeArray));
          if(CopySpectraArray(a, res)) return luaL_error(L, "Error: not enough memory to allocate new spectra\n");
          for(i=0;i<res->NSpectra;i++)
          {
            SpectraScalarMultiply(&(res->Spec[i]), b, 0);
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
            uint32_t i;
            Complex16 *b;
            SpectraTypeArray *res;
            b = (Complex16 *) luaL_checkudata(L, len, "Complex_Type");
            res = (SpectraTypeArray *) lua_newuserdata(L, sizeof(SpectraTypeArray));
            if(CopySpectraArray(a, res)) return luaL_error(L, "Error: not enough memory to allocate new spectra\n");
            for(i=0;i<res->NSpectra;i++)
            {
              SpectraScalarMultiply(&(res->Spec[i]), b->Re, b->Im);
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
            res->Spec = (SpectraType*)malloc(res->NSpectra * sizeof(SpectraType));
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
                  res->Spec[i].Im[j]= a->Spec[i].Re[j] * b->Spec[i].Im[j] + a->Spec[i].Im[j] * b->Spec[i].Re[j];
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
      res->Spec = (SpectraType*)malloc(res->NSpectra * sizeof(SpectraType));
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
          uint32_t i;
          double b;
          SpectraTypeArray *res;
          b = luaL_checknumber(L, len);
          if(b==0) return luaL_error(L, "Error: spectra divided by zero\n");
          b=1/b;
          res = (SpectraTypeArray *) lua_newuserdata(L, sizeof(SpectraTypeArray));
          if(CopySpectraArray(a, res)) return luaL_error(L, "Error: not enough memory to allocate new spectra\n");
          for(i=0;i<res->NSpectra;i++)
          {
            SpectraScalarMultiply(&(res->Spec[i]), b, 0);
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
            uint32_t i;
            Complex16 *b;
            SpectraTypeArray *res;
            b = (Complex16 *) luaL_checkudata(L, len, "Complex_Type");
            res = (SpectraTypeArray *) lua_newuserdata(L, sizeof(SpectraTypeArray));
            if(CopySpectraArray(a, res)) return luaL_error(L, "Error: not enough memory to allocate new spectra\n");
            for(i=0;i<res->NSpectra;i++)
            {
              SpectraScalarMultiply(&(res->Spec[i]), b->Re/(Sqr(b->Re) + Sqr(b->Im)), b->Im/(Sqr(b->Re) + Sqr(b->Im)));
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
            res->Spec = (SpectraType*)malloc(res->NSpectra * sizeof(SpectraType));
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
          Complex16 c, ac;
          SpectraTypeArray *res;
          b = luaL_checknumber(L, len);
          res = (SpectraTypeArray *) lua_newuserdata(L, sizeof(Complex16));
          res->NSpectra = a->NSpectra;
          res->Spec = (SpectraType*)malloc(res->NSpectra * sizeof(SpectraType));
          for(i=0;i<res->NSpectra;i++)
          {
            if(SpectrumAllocate(a->Spec[i], &(res->Spec[i])))
            {
              return luaL_error(L, "Error: not enough memory to allocate new spectra\n");
            }
            for(j=0;j<=res->Spec[i].N;j++)
            {
              c.Re = a->Spec[i].Re[j];
              c.Im = a->Spec[i].Im[j];
              ac = Cpowd(c, b);
              res->Spec[i].Re[j]= ac.Re;
              res->Spec[i].Im[j]= ac.Im;
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
            Complex16 *b,c,ac;
            SpectraTypeArray *res;
            b = (Complex16 *) luaL_checkudata(L, len, "Complex_Type");
            res = (SpectraTypeArray *) lua_newuserdata(L, sizeof(Complex16));
            res->NSpectra = a->NSpectra;
            res->Spec = (SpectraType*)malloc(res->NSpectra * sizeof(SpectraType));
            for(i=0;i<res->NSpectra;i++)
            {
              if(SpectrumAllocate(a->Spec[i], &(res->Spec[i])))
              {
                return luaL_error(L, "Error: not enough memory to allocate new spectra\n");
              }
              for(j=0;j<=res->Spec[i].N;j++)
              {
                c.Re = a->Spec[i].Re[j];
                c.Im = a->Spec[i].Im[j];
                ac = Cpow(c, *b);
                res->Spec[i].Re[j]= ac.Re;
                res->Spec[i].Im[j]= ac.Im;
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
  if (strcmp(str, "Broaden") == 0)
  {
    lua_pushcfunction(L, LuaSpectraBroadenMethod);
    return 1;
  }
  if (strcmp(str, "Element") == 0)
  {
    printf("WARNING: G.Element is deprecated. Please use Spectra.Element instead");
    lua_pushcfunction(L, LuaSpectraElementMethod);
    return 1;
  }
  if (strcmp(str, "Print") == 0)
  {
    lua_pushcfunction(L, LuaSpectraPrintMethod);
    return 1;
  }
  if (strcmp(str, "Shift") == 0)
  {
    lua_pushcfunction(L, LuaSpectraShiftMethod);
    return 1;
  }
  if (strcmp(str, "Sum") == 0)
  {
    printf("WARNING: G.Sum is deprecated. Please use Spectra.Sum instead");
    lua_pushcfunction(L, LuaSpectraSumMethod);
    return 1;
  }
  if (strcmp(str, "MultiplyByEnergy") == 0)
  {
    lua_pushcfunction(L, LuaSpectraMultiplyByEnergyMethod);
    return 1;
  }
  if (strcmp(str, "N") == 0)
  {
    lua_pushinteger(L, a->NSpectra);
    return 1;
  }
  if (strcmp(str, "NE") == 0)
  {
    uint32_t i;
    lua_createtable(L, a->NSpectra, 0);
    for(i=0;i<a->NSpectra;i++)
    {
      lua_pushinteger(L, a->Spec[i].N);
      lua_rawseti(L, -2, i+1);
    }
    return 1;
  }
  if (strcmp(str, "dE") == 0)
  {
    uint32_t i;
    lua_createtable(L, a->NSpectra, 0);
    for(i=0;i<a->NSpectra;i++)
    {
      lua_pushnumber(L, a->Spec[i].dE);
      lua_rawseti(L, -2, i+1);
    }
    return 1;
  }
  if (strcmp(str, "Emin") == 0)
  {
    uint32_t i;
    lua_createtable(L, a->NSpectra, 0);
    for(i=0;i<a->NSpectra;i++)
    {
      lua_pushnumber(L, a->Spec[i].Emin);
      lua_rawseti(L, -2, i+1);
    }
    return 1;
  }
  if (strcmp(str, "Emax") == 0)
  {
    uint32_t i;
    lua_createtable(L, a->NSpectra, 0);
    for(i=0;i<a->NSpectra;i++)
    {
      lua_pushnumber(L, a->Spec[i].Emax);
      lua_rawseti(L, -2, i+1);
    }
    return 1;
  }
  if (strcmp(str, "EminPole") == 0)
  {
    uint32_t i;
    lua_createtable(L, a->NSpectra, 0);
    for(i=0;i<a->NSpectra;i++)
    {
      lua_pushnumber(L, a->Spec[i].EminPole);
      lua_rawseti(L, -2, i+1);
    }
    return 1;
  }
  if (strcmp(str, "EmaxPole") == 0)
  {
    uint32_t i;
    lua_createtable(L, a->NSpectra, 0);
    for(i=0;i<a->NSpectra;i++)
    {
      lua_pushnumber(L, a->Spec[i].EmaxPole);
      lua_rawseti(L, -2, i+1);
    }
    return 1;
  }
  if (strcmp(str, "Gamma") == 0)
  {
    uint32_t i;
    lua_createtable(L, a->NSpectra, 0);
    for(i=0;i<a->NSpectra;i++)
    {
      lua_pushnumber(L, a->Spec[i].Gamma);
      lua_rawseti(L, -2, i+1);
    }
    return 1;
  }
  if (strcmp(str, "Integral") == 0)
  {
    uint32_t i;
    lua_createtable(L, a->NSpectra, 0);
    for(i=0;i<a->NSpectra;i++)
    {
      lua_pushnumber(L, Sqr(a->Spec[i].b[0]));
      lua_rawseti(L, -2, i+1);
    }
    return 1;
  }
  if (strcmp(str, "ResonantEnergy") == 0)
  {
    uint32_t i;
    lua_createtable(L, a->NSpectra, 0);
    for(i=0;i<a->NSpectra;i++)
    {
      lua_pushnumber(L, a->Spec[i].ResonantEnergy);
      lua_rawseti(L, -2, i+1);
    }
    return 1;
  }
  if (strcmp(str, "Name") == 0)
  {
    uint32_t i;
    lua_createtable(L, a->NSpectra, 0);
    for(i=0;i<a->NSpectra;i++)
    {
      lua_pushstring(L, a->Spec[i].Name);
      lua_rawseti(L, -2, i+1);
    }
    return 1;
  }
  return luaL_error(L, "Error operation \".%s()\" not defined on spectra\n",str);
}

static int LuaSpectraToString(lua_State * L)
{
  SpectraTypeArray *spec = (SpectraTypeArray*) luaL_checkudata(L, 1, "Spectra_Type");
  char *Buffer;
  StringPrintSpectra(spec->Spec, spec->NSpectra, &Buffer);
  lua_pushstring(L, Buffer);
  free(Buffer);
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

static const struct luaL_Reg Spectra_func[] = {
  {"Broaden", LuaSpectraBroaden},
  {"Copy", LuaSpectraCopy},
  {"Clone", LuaSpectraCopy},
  {"Element", LuaSpectraElement},
  {"MultiplyByEnergy", LuaSpectraMultiplyByEnergy},
  {"Shift", LuaSpectraShift},
  {"Sum", LuaSpectraSum},
  //  {"ToFunction", LuaSpectraToFunction},
  {"ToTable", LuaSpectraToTable},
  {"FromTable", LuaTableToSpectra}, // deprecated, use ToUserdata
  {"ToUserdata", LuaTableToSpectra},
  {NULL, NULL}
};
