/*  push a new TightBinding object to the lua stack
 *  the tight binding object is empty as its properties need to be set via its ellements
 */
static int LuaNewTightBinding(lua_State * L)
{
  TightBindingDef *TB;
  uint32_t len = (uint32_t)lua_gettop(L);
  if(len==0)
  {
    TB = (TightBindingDef*) lua_newuserdata(L, sizeof(TightBindingDef));
    InitTigthBindingDef(TB);
    luaL_setmetatable(L, "TightBinding_Type");
  }
  else
  {
    return luaL_error(L, "Error in NewTightBinding(arg) unrecognized arguments arg, should be empty\n");
  }
  return 1;
}

static int LuaTightBindingCopy(lua_State * L)
{
  TightBindingDef *TB = (TightBindingDef*) luaL_checkudata(L, -1, "TightBinding_Type");
  TightBindingDef *TBNew;
  TBNew = (TightBindingDef *) lua_newuserdata(L, sizeof(TightBindingDef));
  CopyTightBindingDef(TB, TBNew);
  luaL_setmetatable(L, "TightBinding_Type");
  return 1;
}


static int LuaTightBindingChop(lua_State * L)
{
  double epsilon=10*DBL_EPSILON;
  uint32_t len = (uint32_t)lua_gettop(L);
  struct
  {
    bool Distance;
    bool Deflate;
  } Options;
  LuaOptionsType LuaChopOptions[]={
    {"Distance" ,        "bool",                      &(Options.Distance)},
    {"Deflate"  ,        "bool",                      &(Options.Deflate)},
    {NULL, NULL, NULL}};
  Options.Distance=false;
  Options.Deflate=true;
  if(LuaTestOptionsList(L, -1))
  {
    LuaCheckOptions(L, len, LuaChopOptions);
    len=len-1;
  }
  if(len!=1)
  {
    epsilon=luaL_checknumber(L,len);
    len=len-1;
  }
  TightBindingDef *TB = (TightBindingDef *) luaL_checkudata(L, 1, "TightBinding_Type");
  TightBindingDef *TBnew;
  TBnew = (TightBindingDef*) lua_newuserdata(L, sizeof(TightBindingDef));
  CopyTightBindingDef(TB, TBnew);
  if (Options.Distance)
  {
    ChopTightBindingDefDistance(TBnew, epsilon);
  }
  else
  {
    ChopTightBindingDefSize(TBnew, epsilon);
  }
  luaL_setmetatable(L, "TightBinding_Type");
  return 1;
}

static int LuaTightBindingAddSpin(lua_State * L)
{
  TightBindingDef *TB = (TightBindingDef*) luaL_checkudata(L, -1, "TightBinding_Type");
  TightBindingDef *TBNew;
  TBNew = (TightBindingDef *) lua_newuserdata(L, sizeof(TightBindingDef));
  CopyTightBindingDefAddSpin(TB, TBNew);
  luaL_setmetatable(L, "TightBinding_Type");
  return 1;
}

/*  TightBinding.ShiftOnsiteEnergy(TB, dE) leaves TB untouched and returns a shifted copy
 */
static int LuaTightBindingShiftOnsiteEnergy(lua_State * L)
{
  TightBindingDef *TB = (TightBindingDef*) luaL_checkudata(L, 1, "TightBinding_Type");
  double shift = (double)luaL_checknumber(L, 2);
  TightBindingDef *TBNew;
  TBNew = (TightBindingDef *) lua_newuserdata(L, sizeof(TightBindingDef));
  CopyTightBindingDef(TB, TBNew);
  luaL_setmetatable(L, "TightBinding_Type");
  TightBindingShiftOnsiteEnergy(TBNew, shift);
  return 1;
}

/*  TB.ShiftOnsiteEnergy(dE) shifts TB itself, TB is the upvalue set in LuaTightBindingIndex
 */
static int LuaTightBindingShiftOnsiteEnergyMethod(lua_State * L)
{
  TightBindingDef *TB = (TightBindingDef*) luaL_checkudata(L, lua_upvalueindex(1), "TightBinding_Type");
  double shift = (double)luaL_checknumber(L, 1);
  TightBindingShiftOnsiteEnergy(TB, shift);
  return 0;
}

static int LuaTightBindingDerivative(lua_State * L)
{
  TightBindingDef *TBDx, *TBDy, *TBDz;
  TightBindingDef *TB = (TightBindingDef*) luaL_checkudata(L, -1, "TightBinding_Type");
  TBDx = (TightBindingDef *) lua_newuserdata(L, sizeof(TightBindingDef));
  luaL_setmetatable(L, "TightBinding_Type");
  TBDy = (TightBindingDef *) lua_newuserdata(L, sizeof(TightBindingDef));
  luaL_setmetatable(L, "TightBinding_Type");
  TBDz = (TightBindingDef *) lua_newuserdata(L, sizeof(TightBindingDef));
  luaL_setmetatable(L, "TightBinding_Type");
  TightBindingDerivative(TB, TBDx, TBDy , TBDz);
  return 3;
}

/*  On stack l1, l2, {x,y,z}, {l1_l2_sigma, ...}
 *  Return the Slater Koster hopping matrix elements as table on a basis of Kubic Harmonics
 */
static int LuaSlaterKosterTable(lua_State * L)
{
  double x,y,z;
  double params[5];
  
  uint32_t len = (uint32_t)lua_gettop(L);
  if(len!=4)
  {
    luaL_error(L, "Error in SlaterKosterTable, expected on stack l1, l2, {x,y,z} ,{hopping parameters}\nFound only %d parameters",len);
  }
  uint32_t l1 = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t l2 = (uint32_t)luaL_checkunsigned(L, 2);
  if((l1>4)||(l2>4))
  {
    luaL_error(L,"Slater Koster Table only implemented for s,p,d,f and g orbitals");
  }
  luaL_checktype(L, 3, LUA_TTABLE);
  len = (uint32_t)lua_objlen(L, 3);
  if(len !=3)
  {
    luaL_error(L, "Error in SlaterKosterTable, expected the directional cosins {x,y,z} on position 3 on the stack. Found a table with length %d\n",len);
  }
  lua_rawgeti(L, 3, 1);
  x = (double)luaL_checknumber(L, -1);
  lua_pop(L, 1);
  lua_rawgeti(L, 3, 2);
  y = (double)luaL_checknumber(L, -1);
  lua_pop(L, 1);
  lua_rawgeti(L, 3, 3);
  z = (double)luaL_checknumber(L, -1);
  lua_pop(L, 1);
  if(Abs(x)+Abs(y)+Abs(z)==0)
  {
    luaL_error(L, "Error in SlaterKosterTable, Hopping vector should not be {0,0,0}\n",len);
  }
  luaL_checktype(L,4, LUA_TTABLE);
  len = (uint32_t)lua_objlen(L, 4);
  if(len !=Min(l1,l2)+1)
  {
    luaL_error(L, "Error in SlaterKosterTable, expected the hopping parameters on position 4 on the stack.\nThe number of hopping parameters should be Min(l1,l2)+1.\nFound a table with length %d\n",len);
  }
  for (int i=0;i<=Min(l1,l2);i++)
  {
    lua_rawgeti(L, 4, i+1);
    params[i] = (double)luaL_checknumber(L, -1);
    lua_pop(L, 1);
  }
  CompactMatrixType *H;
  H = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
  if( SlaterKosterTable(l1, l2, x, y, z, params, H) )
  {
    luaL_error(L, "SlaterKosterTable failed\nFatal Error\n");
  }
  luaL_setmetatable(L, "Matrix_Type");
  H->OnLuaStack=true;
  return 1;
}



static int LuaFileReadStuttgartCTRL(lua_State * L)
{
  StuttgartCTRL *CTRL = NULL;
  uint32_t len = (uint32_t)lua_gettop(L);
  if(len!=1) return luaL_error(L, "Error in FileReadStuttgartCTRL(arg).\n The argument should be the file name.");

  const char *flname = luaL_checkstring(L, 1);
  CTRL = (StuttgartCTRL*) lua_newuserdata(L, sizeof(StuttgartCTRL));
  FileReadStuttgartCTRL(CTRL, flname);
  luaL_setmetatable(L, "StuttgartCTRL_Type");

  return 1;
}

static int LuaFileReadStuttgartHAMR(lua_State * L)
{
  StuttgartHAMR *HAMR = NULL;
  uint32_t len = (uint32_t)lua_gettop(L);
  if(len!=1) return luaL_error(L, "Error in FileReadStuttgartHAMR(arg).\n The argument should be the file name.");

  const char *flname = luaL_checkstring(L, 1);
  HAMR = (StuttgartHAMR*) lua_newuserdata(L, sizeof(StuttgartHAMR));
  FileReadStuttgartHAMR(HAMR, flname);
  luaL_setmetatable(L, "StuttgartHAMR_Type");

  return 1;
}

static int LuaFileReadDresdenFPLO(lua_State * L)
{
  DresdenFPLO *FPLO = NULL;
  uint32_t len = (uint32_t)lua_gettop(L);
  if(len!=1) return luaL_error(L, "Error in FileReadDresdenFPLO(arg).\n The argument should be the file name of the wan output.");

  const char *flname = luaL_checkstring(L, 1);
  FPLO = (DresdenFPLO*) lua_newuserdata(L, sizeof(DresdenFPLO));
  if(FileReadDresdenFPLO(FPLO, flname)) return luaL_error(L, "Error while reading or opening file %s in FileReadDresdenFPLO(arg).",flname);
  luaL_setmetatable(L, "DresdenFPLO_Type");
  return 1;
}

static int LuaReadFPLO(lua_State * L)
{
  DresdenFPLO *FPLO = NULL;
  uint32_t len = (uint32_t)lua_gettop(L);
  if(len!=1) return luaL_error(L, "Error in ReadFPLO(arg).\n The argument should be the folder where the FPLO calculation is done");
  const char *dirname = luaL_checkstring(L, 1);
  FPLO = (DresdenFPLO*) lua_newuserdata(L, sizeof(DresdenFPLO));
  if(ReadFPLO(FPLO, dirname)) return luaL_error(L, "Error while reading or opening files in directory %s in ReadFPLO(arg).",dirname);
  luaL_setmetatable(L, "DresdenFPLO_Type");
  return 1;
}

static int LuaTightBindingDefFromStuttgartCTRLHAMR(lua_State * L)
{
  TightBindingDef *TB;
  StuttgartCTRL *CTRL = NULL;
  StuttgartHAMR *HAMR = NULL;
  uint32_t len = (uint32_t)lua_gettop(L);
  if (len != 2)
  {
    return luaL_error(L, "Error in TightBindingDefFromStuttgartCTRLHAMR(arg).\nIt takes two arguments (CTRL and HAMR).");
  }

  CTRL = (StuttgartCTRL *) luaL_checkudata(L, 1, "StuttgartCTRL_Type");
  HAMR = (StuttgartHAMR *) luaL_checkudata(L, 2, "StuttgartHAMR_Type");

  TB = (TightBindingDef*) lua_newuserdata(L, sizeof(TightBindingDef));
  TightBindingDefFromStuttgartCTRLHAMR(*CTRL, *HAMR, TB);
  luaL_setmetatable(L, "TightBinding_Type");

  return 1;
}

static int LuaTightBindingDefFromDresdenFPLO(lua_State * L)
{
  TightBindingDef *TB;
  DresdenFPLO *FPLO = NULL;
  uint32_t len = (uint32_t)lua_gettop(L);
  if (len != 1)
  {
    return luaL_error(L, "Error in TightBindingDefFromDresdenFPLO(arg).\nIt takes a DresdenFPLO type as the single argument.");
  }

  FPLO = (DresdenFPLO *) luaL_checkudata(L, 1, "DresdenFPLO_Type");
  TB = (TightBindingDef*) lua_newuserdata(L, sizeof(TightBindingDef));

  TightBindingDefFromDresdenFPLO(*FPLO, TB);
  luaL_setmetatable(L, "TightBinding_Type");

  return 1;
}


int ReturnIndex(DresdenFPLO *FPLO, int sort, int n, int l, int m ,int spin)
{
  int i;
  for(i=0;i<FPLO->lDOS.N;i++)
  {
    if( (FPLO->lDOS.sort[i] == sort ) &&
        (FPLO->lDOS.n[i]    == n    ) &&
        (FPLO->lDOS.l[i]    == l    ) &&
        (FPLO->lDOS.m[i]    == m    ) &&
        (FPLO->lDOS.spin[i] == spin )    ) return i;
  }
  printf("Error: sort nlmsigma Index not found\n %d %d %d %d %d\n",sort,n,l,m,spin);
  return -1;
}


/*  On stack is a Dresden FPLO type and a list or string of orbitals to include
 *  The list can be a single string of the form "3d"
 *  Or a list giving the orbitals more selectively
 */
static int LuaHybridizationFunctionFromFPLO(lua_State * L)
{
  int i;
  int site;
  int *OrbitalList=NULL;
  int NOrbitals;
  const char *nlms;
  int n,l,m,s;
  bool hasm,hass;
  DresdenFPLO *FPLO = NULL;
  ResponsefunctionType *Hyb;
  BlockListOfPoles *LP;
  uint32_t len = (uint32_t)lua_gettop(L);
  // possible read options from the stack at position 3
  struct
  {
    bool Normalize=true;
    bool SetNegativeToZero=true;
    bool MergeNegativeToNeighbour=false;
    CompactVectorType EnergyWindow={0,0,NULL};
  }Options;
  LuaOptionsType LuaOptions[]={
    {"Normalize",                "bool",   &(Options.Normalize)},
    {"SetNegativeToZero",        "bool",   &(Options.SetNegativeToZero)},
    {"MergeNegativeToNeighbour", "bool",   &(Options.MergeNegativeToNeighbour)},
    {"EnergyWindow",             "Vector", &(Options.EnergyWindow)},
    {NULL, NULL, NULL}};
  if(len==4)
  {
    LuaCheckOptions(L, 4, LuaOptions);
    len--;
  }
  if(len!=3) luaL_error(L, "Error in HybridizationFunctionFromFPLO, expected on stack an FPLO input, the site number and a definition of the orbitals to include and possible a list of options\n");
  // check if energy window is a valid window
  if(Options.EnergyWindow.Ni!=0)
  {
    if(Options.EnergyWindow.Ni!=2) return luaL_error(L, "ERROR: Energy window should be a list of two real numbers related to the minimal and maximum energy used of the DOS to generate the hybridization function\n");
    if(Options.EnergyWindow.QComplex!=0) return luaL_error(L, "ERROR: Energy window should be a list of two real numbers related to the minimal and maximum energy used of the DOS to generate the hybridization function. Numbers are not real\n");
  }
  // read the FPLO input from the first position of the stack
  FPLO = (DresdenFPLO *) luaL_checkudata(L, 1, "DresdenFPLO_Type");
  // read the site from the stack
  site = (uint32_t)luaL_checkinteger(L,2);
  // read the list of orbitals to include from the second position of the stack
  // Options are nl in the form of "1s" or "3d"
  // A specific nlm value in the form of "1s0" "2p-1", "2p1" or "2p+1"
  // A specific nlsigma value in the form of "1s_1" or "1s_2"
  // A specific nlmsigna falue in the form of "1s0_1" or "2p-1_2"
  // Or a table of these strings
  if( lua_istable(L,3) )
  {
    return luaL_error(L, "ERROR: currently no tables as input for orbitals suported\n");
  }
  else
  {
    nlms = lua_tostring(L,3);
    ReadNonRelativisticOrbital(nlms, &n, &l, &m, &s, &hasm, &hass);
    if(hasm==true)
    {
      NOrbitals = 2;
      if( (FPLO->Options.MSpin==1) || (hass==true) ) NOrbitals=1;
    }
    else
    {
      NOrbitals = 2 * l + 1;
      if( (FPLO->Options.MSpin==2) && (hass==false) ) NOrbitals*=2;
    }
    OrbitalList=(int*)malloc(NOrbitals*sizeof(int));
    if(hasm==true)
    {
      if (hass==true)
      {
        OrbitalList[0] = ReturnIndex(FPLO, site, n, l, m, s);
      }
      else if (FPLO->Options.MSpin==1)
      {
        OrbitalList[0] = ReturnIndex(FPLO, site, n, l, m, 1);
      }
      else
      {
        OrbitalList[0] = ReturnIndex(FPLO, site, n, l, m, 1);
        OrbitalList[1] = ReturnIndex(FPLO, site, n, l, m, 2);
      }
    }
    else
    {
      if (hass==true)
      {
        for(m=-l;m<=l;m++)
        {
          OrbitalList[l+m] = ReturnIndex(FPLO, site, n, l, m, s);
        }
      }
      else if (FPLO->Options.MSpin==1)
      {
        for(m=-l;m<=l;m++)
        {
          OrbitalList[l+m] = ReturnIndex(FPLO, site, n, l, m, 1);
        }
      }
      else
      {
        for(m=-l;m<=l;m++)
        {
          OrbitalList[2*l+0+2*m] = ReturnIndex(FPLO, site, n, l, m, 1);
          OrbitalList[2*l+1+2*m] = ReturnIndex(FPLO, site, n, l, m, 2);
        }
      }
    }
  }
  for(i=0;i<NOrbitals;i++)
  {
    if (OrbitalList[i]==-1) return luaL_error(L, "ERROR: orbital index not found\n");
  }
  // push Hybridizaton function to the stack
  Hyb = (ResponsefunctionType*) lua_newuserdata(L, sizeof(ResponsefunctionType));
  // fill the hybridization function with data
  Hyb->type = 'l';
  Hyb->Block=true;
  LP = (BlockListOfPoles*)malloc(sizeof(BlockListOfPoles));
  Hyb->data = LP;
  Hyb->SuggestedGrid=NULL;
  Hyb->NSuggestedGridPoints=0;
  // set the dimension in LP to the number of orbitals;
  LP->BlockSize=NOrbitals;
  HybridizationFunctionFromFPLO(FPLO, LP, OrbitalList, Options.MergeNegativeToNeighbour, Options.SetNegativeToZero, Options.Normalize, Options.EnergyWindow);
  luaL_setmetatable(L, "ResponseFunction_Type");
  free(OrbitalList);
  return 1;
}




/* on stack:
 *   TB : TightBinding_Type
 *
 *   Options
 *       Emin
 *       Emax
 *       NE
 *       "Nk" -> {Nkx,Nky,Nkz} number of k-points       standard {40,40,40}
 *       Method -> "mesh", ...                          standard "mesh"
 *       Character -> {matrix} or "full" or "total" standard "full" rho is psi^* (x) psi, rho2 = M^T . psi^* (x) psi . M^*
 *       Broadening -> "None", {"Gauss",FMHM}, {"Lorentz, FWHM}
 *   returns a Spectra_Type
 */
static int LuaCalculateDOS(lua_State *L)
{
  //TODO: character
  uint32_t i, len;
  TightBindingDef *a;
  SpectraTypeArray *DOS;

  struct
  {
    int *Nk;
    //char* Method;   //only method implemented is the mesh/grid.
    const char* Character;
    const char* Broadening;
    double FWHM;
    double Emin;
    double Emax;
    uint32_t NE;
  } Options;
  Options.Nk = (int*)malloc(3*sizeof(int));
  Options.Nk[0]=40;
  Options.Nk[1]=40;
  Options.Nk[2]=40;
  //Options.Method    = "Mesh";
  Options.Character = "Full";
  Options.Broadening= "None";
  Options.FWHM=0;
  Options.Emin=-10;
  Options.Emax=10;
  Options.NE=2000;
  LuaOptionsType LuaCalculateDOSOptions[]={
    {"Nk", "array-uint32_t-3", &(Options.Nk)},
    //{"Method",     "char*",    &(Options.Method)},
    {"Character",  "char*",    &(Options.Character)},
    {"Broadening", "char*",    &(Options.Broadening)},
    {"FWHM",       "double",   &(Options.FWHM)},
    {"Emin",       "double",   &(Options.Emin)},
    {"Emax",       "double",   &(Options.Emax)},
    {"NE",         "uint32_t", &(Options.NE)},
    {NULL, NULL, NULL}};

  len = (uint32_t)lua_gettop(L);
  if(!(len==1||len==2))
  {
    return luaL_error(L, "Error CalculateDOS: called with %d arguments, 1 or 2 expected (TightBindingDef, and options).\n",len);
  }
  a = (TightBindingDef*) luaL_checkudata(L, 1, "TightBinding_Type");
  if(len==2)
  {
    LuaCheckOptions(L, -1, LuaCalculateDOSOptions);
  }

  DOS = (SpectraTypeArray *) lua_newuserdata(L, sizeof(SpectraTypeArray));
  luaL_setmetatable(L, "Spectra_Type");
  DOS->Name=(char*)malloc(4*sizeof(char));
  DOS->Name[0]='D';
  DOS->Name[1]='O';
  DOS->Name[2]='S';
  DOS->Name[3]='\0';
  if(strcicmp(Options.Character, "Full")==0)
  {
    DOS->NSpectra = 1;
    DOS->Spec= (SpectraType*)malloc(DOS->NSpectra*sizeof(SpectraType));
    DOS->Spec[0].N    = Options.NE;
    DOS->Spec[0].Emin = Options.Emin;
    DOS->Spec[0].Emax = Options.Emax;
    DOS->Spec[0].dE   =(DOS->Spec[0].Emax-DOS->Spec[0].Emin)/DOS->Spec[0].N;
    SpectrumAllocate(DOS->Spec[0], &(DOS->Spec[0]));
    strcpy(DOS->Spec[0].Name, "Total DOS");

    TightBindingDOS(*a, Options.Nk, &(DOS->Spec[0]));
    // Note there are some inconsistencies regarding how the DOS is normalized in the original functions.
    SpectraScalarMultiply(&(DOS->Spec[0]), DOS->Spec[0].dE, 0.);
  }
  else
  {
    printf("Warning CalculateDOS: character-resolved DOS only outputs the -Im part of the spectra to the real part.");
    DOS->NSpectra = (*a).NFermion+1;
    printf("%d\n", DOS->NSpectra);
    DOS->Spec= (SpectraType*)malloc(DOS->NSpectra*sizeof(SpectraType));
    DOS->Spec[0].N    = Options.NE;
    DOS->Spec[0].Emin = Options.Emin;
    DOS->Spec[0].Emax = Options.Emax;
    DOS->Spec[0].dE   =(DOS->Spec[0].Emax-DOS->Spec[0].Emin)/DOS->Spec[0].N;
    for(i=0;i<DOS->NSpectra-1;i++)
    {
      SpectrumAllocate(DOS->Spec[0], &(DOS->Spec[i]));
      sprintf(DOS->Spec[i].Name, "DOS %i/%i", i, DOS->NSpectra-1);
    }
    SpectrumAllocate(DOS->Spec[0], &(DOS->Spec[DOS->NSpectra-1]));
    strcpy(DOS->Spec[DOS->NSpectra-1].Name, "Total DOS");

    TightBindingDOSCharacter(*a, Options.Nk, DOS->Spec, a->NFermion+1, &EigenFunctionToCharactersStandard);
  }

  if(strcicmp(Options.Broadening, "Gaussian")==0||
     strcicmp(Options.Broadening, "Gauss")==0||
     strcicmp(Options.Broadening, "G")==0)
  {
    for(i=0;i<DOS->NSpectra;i++)
    {
      BroadenG(&(DOS->Spec[i]), Options.FWHM);
    }
  }
  else if(strcicmp(Options.Broadening, "Lorentzian")==0||
          strcicmp(Options.Broadening, "Lorentz")==0||
          strcicmp(Options.Broadening, "L")==0)
  {
    for(i=0;i<DOS->NSpectra;i++)
    {
      BroadenL(&(DOS->Spec[i]), Options.FWHM);
    }
  }
  else if(strcicmp(Options.Broadening, "None")!=0)
  {
    printf("\nWarning: Unrecognized broadening method, using Gaussian.\n");
    for(i=0;i<DOS->NSpectra;i++)
    {
      BroadenL(&(DOS->Spec[i]), Options.FWHM);
    }
  }

  return 1;
}

/* On stack we have
 *     1) a tight binding Hamiltonian, HTB
 *     2) a path in k-space of the form {{Label, k-point}, N1 , {Label, k-point} , N2, {Label, k-point}}
 *     3) a function that assigns a double value (Character) between 0 and 1 to a vector of length HTB.NFermion  (used as Black-White band-color or band-thickness)
 *     3,4)     same as 3 (used as Black-White band-color and band-thickness)
 *     3,4,5)   same as 3 (used as RGB band-color)
 *     3,4,5,6) same as 3 (used as RGB band-color and band-thickness)
 *     -1) possible options
 *  In place of a function a matrix P of dimension HTB.NFermion can be given, the character of eigenstate
 *  psi is then Re(sum_ij psi_i^* P_ij psi_j).
 *
 *     we return
 *     1) a matrix of dimensions HTB.NFermion by sum_i N_i + (#(kpath)+1)/2 with the band eigenbvalues
 *     2) (3,4,5,6) a matrix of dimensions HTB.NFermion by sum_i N_i + (#(kpath)+1)/2 with the band character
 *
 *  With a self energy on position 3, in any of the forms CalculateG accepts, followed by 0, 1 or 3
 *  characters, see LuaCalculateBandsSelfEnergy.
 */
static bool LuaIsResponseFunctionObject(lua_State *L, int idx);
/*  Re(sum_ij u_i^* P_ij u_j) for u = ur + I ui of length n and P an n by n matrix, real or complex
 */
static double CharacterQuadraticForm(CompactMatrixType const * const P, double const * const ur, double const * const ui, int const n)
{
  double ret=0, pr, pi;
  for(int i=0;i<n;i++)
  {
    for(int j=0;j<n;j++)
    {
      pr = (P->QComplex) ? (P->V[i*P->Nj+j].Re) : (P->VRe[i*P->Nj+j]);
      pi = (P->QComplex) ? (P->V[i*P->Nj+j].Im) : (0);
      // Re( (ur_i - I ui_i) (pr + I pi) (ur_j + I ui_j) )
      ret += pr * (ur[i]*ur[j] + ui[i]*ui[j]) - pi * (ur[i]*ui[j] - ui[i]*ur[j]);
    }
  }
  return ret;
}
static int LuaCalculateBandsSelfEnergy(lua_State *L, TightBindingDef const * const HTB, BandstructurePathType const * const kList, int const Nktot, int const iSigma, int const iChar, int const NCharacter);
static int LuaCalculateBands(lua_State *L)
{
  // check (and read) the input
  LuaOptionsType LuaCalculateBandOptions[]={
  {NULL, NULL, NULL}};
  int len = (int)lua_gettop(L);
  if (LuaTestOptionsList(L,-1))
  {
    LuaCheckOptions(L,-1,LuaCalculateBandOptions);
    len--;
  }
  if(len<2)
  {
    return luaL_error(L, "Error CalculateBands: called with %d arguments (excluding options), minimal 2 expected (TightBindingDef, and k-path).\n",len);
  }
  if(len>6)
  {
    return luaL_error(L, "Error CalculateBands: called with %d arguments (excluding options), maximal 6 expected (TightBindingDef, k-path, (0 to 4 times) a function that transforms an eigenstate to a character (value between 0 and 1) ).\n",len);
  }
  int NCharacter = len-2;
  TightBindingDef *HTB = (TightBindingDef*) luaL_checkudata(L, 1, "TightBinding_Type");
  // a self energy on position 3 is a response function or a table per atom, a character is a function or a matrix
  bool WithSigma = (len>=3) && (LuaIsResponseFunctionObject(L, 3) || (lua_istable(L, 3) && !lua_isfunction(L, 3)));
  if(WithSigma)
  {
    NCharacter = len-3;
    if(!(NCharacter==0 || NCharacter==1 || NCharacter==3))
    {
      return luaL_error(L, "Error CalculateBands: with a self energy 0, 1 or 3 characters are possible, not %d (a band-thickness is not defined for a spectral function)\n", NCharacter);
    }
  }
  // check (and read) the k-list
  BandstructurePathType kList;
  if(!(lua_istable(L, 2)))
  {
    return luaL_error(L, "Error CalculateBands: input at position 2 should be a path in k-space of the form {{Label, k-point}, N1 , {Label, k-point} , N2, {Label, k-point}}).\n");
  }
  int lenklist = (int)lua_objlen(L, 2);
  int Nktot = (lenklist+1)/2;
  if((lenklist&1)==0)
  {
    return luaL_error(L, "Error CalculateBands: input at position 2 should be a path in k-space of the form {{Label, k-point}, N1 , {Label, k-point} , N2, {Label, k-point}}).\nThe length should be odd and not %d\n",lenklist);
  }
  kList.NQ = (lenklist+1)/2;
  kList.Q  = (BandstructurePointType*)malloc(kList.NQ*sizeof(BandstructurePointType));
  kList.Units = HTB->Cell.Units;
  for(int i=1;i<=lenklist;i+=2)
  {
    // test k-point
    lua_rawgeti(L, 2, i); // put the table {Label, k-vec} on the stack
    if(!(lua_istable(L, -1)))
    {
      lua_pop(L,1);
      return luaL_error(L, "Wrong format for k-list at position %d\n",i);
    }
    if(!(((int)lua_objlen(L, -1))==2))
    {
      lua_pop(L,1);
      return luaL_error(L, "Wrong format for k-list at position %d\n",i);
    }
    lua_rawgeti(L,-1,1); // put the string Label on the stack
    if(!(lua_type(L, -1)==LUA_TSTRING))
    {
      lua_pop(L,1);
      lua_pop(L,1);
      return luaL_error(L, "Wrong format for k-list at position %d\n",i);
    }
    kList.Q[(i-1)/2].label = lua_tolstring(L, -1, NULL)[0];
    lua_pop(L,1); // remove the string from the stack
    lua_rawgeti(L,-1,2); // put the k-point on the stack
    if(!(lua_istable(L, -1)))
    {
      lua_pop(L,1);
      return luaL_error(L, "Wrong format for k-list at position %d\n",i);
    }
    if(!(((int)lua_objlen(L, -1))==3))
    {
      lua_pop(L,1);
      return luaL_error(L, "Wrong format for k-list at position %d\n",i);
    }
    for(int j=0;j<3;j++)
    {
      lua_rawgeti(L, -1, j+1);
      kList.Q[(i-1)/2].q[j] = (double)luaL_checknumber(L, -1);
      lua_pop(L,1);
    }
    lua_pop(L,1); // remove the k-point from the stack
    lua_pop(L,1); // remove the table {Lable, k-vec} from the stack
    // test integer N
    if(i==lenklist)
    {
      kList.Q[(i-1)/2].Nq = 0;
    }
    else
    {
      lua_rawgeti(L, 2, i+1); // put the integer Ni on the stack
      if (!lua_isnumber(L, -1))
      {
        return luaL_error(L, "Wrong format for k-list at position %d\nExpected an integer\n",i+1);
      }
      kList.Q[(i-1)/2].Nq = (int) lua_tointeger(L, -1);
      if(kList.Q[(i-1)/2].Nq<0)
      {
        return luaL_error(L, "Number of k-points in the path can not be negative at position %d\n",i+1);
      }
      Nktot += kList.Q[(i-1)/2].Nq;
      lua_pop(L,1);
    }
  }
  // finished reading the k-list
  if(WithSigma)
  {
    int ret = LuaCalculateBandsSelfEnergy(L, HTB, &kList, Nktot, 3, 4, NCharacter);
    free(kList.Q);
    return ret;
  }
  // produce output
  // we will place these on the lua stack later by copying the pointers
  CompactMatrixType val;
  CompactMatrixType *Character;
  Character = (CompactMatrixType *)malloc(NCharacter * sizeof(CompactMatrixType));
  val.Ni = Nktot;
  val.Nj = HTB->NFermion;
  val.QComplex = 0;
  InitCompactMatrix(&val);
  strcpy(val.Name, "Band Energies");
  for(int i=0;i<NCharacter;i++)
  {
    Character[i].Ni = Nktot;
    Character[i].Nj = HTB->NFermion;
    Character[i].QComplex = 0;
    InitCompactMatrix(&(Character[i]));
    sprintf(Character[i].Name, "Band Character %d",i+1);
  }
  // check that position 3 to 2+NCharacterare indeed functions or matrices of dimension NFermion
  for(int i=0;i<NCharacter;i++)
  {
    CompactMatrixType *P = (CompactMatrixType*)luaL_testudata(L, 3+i, "Matrix_Type");
    if (P!=NULL)
    {
      if(P->Ni!=HTB->NFermion || P->Nj!=HTB->NFermion)
      {
        return luaL_error(L, "Error CalculateBands: the character matrix on position %d is %d by %d, it should be %d by %d\n", 3+i, P->Ni, P->Nj, HTB->NFermion, HTB->NFermion);
      }
    }
    else if (!lua_isfunction(L, 3+i))
    {
      return luaL_error(L, "Expected a function that takes an eigenstate and returns a character (number of type double bewteen 0 and 1), or a matrix (Matrix.ToUserdata) on position %d\n",3+i);
    }
  }
  // loop over all k-points
  {
    CompactMatrixType H, *T;
    T = (CompactMatrixType*)lua_newuserdata(L,sizeof(CompactMatrixType));
    T->OnLuaStack=false; // warning danger -- we should not loose track of T on the lua stack before we make T save to be garbage collected
    luaL_setmetatable(L, "Matrix_Type");
    int iT = lua_gettop(L);
    int ik=0;
    double k[3];
    for(int i=0; i<kList.NQ; i++)
    {
      for(int j=0; j<=kList.Q[i].Nq; j++)
      {
        if(j==0)
        {
          k[0] = kList.Q[i].q[0];
          k[1] = kList.Q[i].q[1];
          k[2] = kList.Q[i].q[2];
        }
        else
        {
          k[0] = kList.Q[i].q[0] * (kList.Q[i].Nq + 1 - j) / (kList.Q[i].Nq + 1) + kList.Q[i+1].q[0] * j / (kList.Q[i].Nq + 1) ;
          k[1] = kList.Q[i].q[1] * (kList.Q[i].Nq + 1 - j) / (kList.Q[i].Nq + 1) + kList.Q[i+1].q[1] * j / (kList.Q[i].Nq + 1) ;
          k[2] = kList.Q[i].q[2] * (kList.Q[i].Nq + 1 - j) / (kList.Q[i].Nq + 1) + kList.Q[i+1].q[2] * j / (kList.Q[i].Nq + 1) ;
        }
        CreateHamiltonianFromTightBindingDef(*HTB, k, &H);
        CompactMatrixEigensystem(&H, T, &(val.VRe[ik*val.Nj]));
        // transform eigenstates to character
        for(int c=0;c<NCharacter;c++)
        {
          CompactMatrixType *P = (CompactMatrixType*)luaL_testudata(L, 3+c, "Matrix_Type");
          if(P!=NULL)
          {
            double ur[T->Nj], ui[T->Nj];
            for(int l=0;l<T->Nj;l++)
            {
              for(int m=0;m<T->Nj;m++)
              {
                ur[m] = (T->QComplex) ? (T->V[l*T->Nj+m].Re) : (T->VRe[l*T->Nj+m]);
                ui[m] = (T->QComplex) ? (T->V[l*T->Nj+m].Im) : (0);
              }
              Character[c].VRe[ik*T->Nj+l] = CharacterQuadraticForm(P, ur, ui, T->Nj);
            }
            continue;
          }
          lua_pushvalue(L, 3+c);             // Lua function
          lua_createtable(L, 0, 2);
          lua_pushvalue(L, iT);              // pointer to T
          lua_setfield(L,-2,"__Matrix");
          lua_pushnumber(L, 0);
          lua_setfield(L,-2,"__Row");
          lua_getglobal(L, "MatrixRowMeta");
          lua_setmetatable(L, -2);
          for(int l=0;l<T->Nj;l++)
          {
            lua_pushvalue(L, -2); // copy function to stack as pcall removes it
            lua_pushvalue(L, -2); // copy matrix row to stack as pcal removes it
            lua_pushnumber(L, l+1);
            lua_setfield(L,-2,"__Row");
            if (lua_pcall(L, 1, 1, 0) != LUA_OK)
            {
              const char *err = lua_tostring(L, -1);
              luaL_error(L, "Error in character function %d: %s", c+1, err);
            }
            Character[c].VRe[ik*T->Nj+l] = (double)luaL_checknumber(L, -1);
            lua_pop(L,1);
          }
          lua_pop(L,1);
          lua_pop(L,1);
        }
        FreeCompactMatrix(T);
        FreeCompactMatrix(&H);
        ik++;
      }
    }
    // T is a matrix on the Lua Stack, but the last action was to free T. We need to allocate T and give it to the Lua stack such that Lua can Garbage collect T. :-)
    T->Ni=0;
    T->Nj=0;
    T->QComplex=0;
    InitCompactMatrix(T);
    T->OnLuaStack=true;
  }
  // place result on lua stack
  CompactMatrixType* M;
  M=(CompactMatrixType*)lua_newuserdata(L,sizeof(CompactMatrixType));
  luaL_setmetatable(L, "Matrix_Type");
  *M=val;
  M->OnLuaStack=true;
  for(int i=0;i<NCharacter;i++)
  {
    M=(CompactMatrixType*)lua_newuserdata(L,sizeof(CompactMatrixType));
    luaL_setmetatable(L, "Matrix_Type");
    *M=Character[i];
    M->OnLuaStack=true;
  }
  // free matrices
  free(Character);
  free(kList.Q);
  return 1+NCharacter;
}

/*
 *
 */
static int LuaCalculateFermi(lua_State *L)
{
  //TODO: all

  return 1;
}

/*  true if the value at stack index idx is a response function, either as user data or as a lua table
 *  with a string field type. A table of response functions has no field type.
 */
static bool LuaIsResponseFunctionObject(lua_State *L, int idx)
{
  bool ret;
  if(luaL_testudata(L, idx, "ResponseFunction_Type")!=NULL) return true;
  if(lua_istable(L, idx)==0) return false;
  lua_getfield(L, idx, "type");
  ret = (lua_type(L, -1)==LUA_TSTRING);
  lua_pop(L, 1);
  return ret;
}

static char const * ResponseFunctionTypeName(char const type)
{
  switch(type|32)
  {
    case 'l': return "ListOfPoles";
    case 't': return "Tri";
    case 'a': return "And";
    case 'n': return "Nat";
    case 'd': return "DoubleTri";
  }
  return "unknown";
}

/*  the self energies per atom or per shell collected by LuaReadSelfEnergyPerAtom */
typedef struct SelfEnergyEntriesType
{
  uint32_t N;                     /* number of entries */
  ResponsefunctionType **Sigma;   /* the self energy of each entry */
  uint32_t *Offset;               /* index of the first orbital of each entry */
  int *Atom;                      /* atom of each entry */
  int *Shell;                     /* shell of each entry, -1 for all shells of the atom */
  bool *Owned;                    /* true if Sigma was allocated here and has to be freed */
  double mu;                      /* mu of the first entry */
}
SelfEnergyEntriesType;

/*  writes "atom i (name)" or "atom i (name) shell j (name)" into where, which has room for 2*MAXStrLength+64 */
static void SelfEnergyEntryName(char * const where, TightBindingDef const * const TB, int const iatom, int const ishell)
{
  AtomType const * const Atom = &(TB->Cell.Atom[iatom]);
  if(ishell<0) snprintf(where, 2*MAXStrLength+64, "atom %d (%s)", iatom+1, Atom->Name);
  else         snprintf(where, 2*MAXStrLength+64, "atom %d (%s) shell %d (%s)", iatom+1, Atom->Name, ishell+1, Atom->ShellName[ishell]);
}

/*  Read the value at stack index idx as the self energy of atom iatom and shell ishell (ishell==-1 for all
 *  shells of the atom) and add it to E. The number 0 is a zero self energy and adds nothing.
 */
static void LuaReadSelfEnergyEntry(lua_State *L, int const idx, TightBindingDef const * const TB, int const iatom, int const ishell, SelfEnergyEntriesType * const E)
{
  AtomType const * const Atom = &(TB->Cell.Atom[iatom]);
  ResponsefunctionType *rf;
  char where[2*MAXStrLength+64], first[2*MAXStrLength+64];
  int j, dim=0;
  SelfEnergyEntryName(where, TB, iatom, ishell);
  if(lua_type(L, idx)==LUA_TNUMBER)
  {
    if(lua_tonumber(L, idx)!=0) luaL_error(L, "Error CalculateG: the self energy of %s is a number, only 0 is allowed, for a zero self energy\n", where);
    return;
  }
  if(!LuaIsResponseFunctionObject(L, idx))
  {
    if(ishell<0) luaL_error(L, "Error CalculateG: the self energy of %s should be a response function, a table of response functions per shell or 0\n", where);
    else         luaL_error(L, "Error CalculateG: the self energy of %s should be a response function or 0\n", where);
  }
  rf = luaL_checkResponsefunction(L, idx);
  if(ishell<0) for(j=0;j<Atom->Nshell;j++) dim += Atom->shelldegeneracy[j];
  else         dim = Atom->shelldegeneracy[ishell];
  if((int)ResponseFunctionGetBlockSize(rf)!=dim)
  {
    luaL_error(L, "Error CalculateG: the self energy of %s has dimension %d, but %s has %d orbitals\n", where, (int)ResponseFunctionGetBlockSize(rf), where, dim);
  }
  if(E->N==0)
  {
    E->mu = ResponseFunctionGetMu(rf);
  }
  else
  {
    SelfEnergyEntryName(first, TB, E->Atom[0], E->Shell[0]);
    if((rf->type|32)!=(E->Sigma[0]->type|32))
    {
      luaL_error(L, "Error CalculateG: the self energies per atom or shell must all be of the same type, but %s has type %s and %s has type %s. Use ResponseFunction.ChangeType to give them one type\n",
                 first, ResponseFunctionTypeName(E->Sigma[0]->type), where, ResponseFunctionTypeName(rf->type));
    }
    if(fabs(ResponseFunctionGetMu(rf)-E->mu) > 10*DBL_EPSILON*Max(1.0,fabs(E->mu)))
    {
      luaL_error(L, "Error CalculateG: the self energies per atom or shell must all have the same mu, but %s has mu = %f and %s has mu = %f\n",
                 first, E->mu, where, ResponseFunctionGetMu(rf));
    }
  }
  E->Sigma[E->N]  = rf;
  E->Offset[E->N] = (Atom->Nshell>0) ? (TB->FermionIndex[iatom][(ishell<0)?0:ishell]) : (0);
  E->Atom[E->N]   = iatom;
  E->Shell[E->N]  = ishell;
  E->Owned[E->N]  = ((rf->type&32)==0); // an upper case type was read from a lua table and allocated by luaL_checkResponsefunction
  E->N++;
}

/*  Read a self energy given per atom, or per atom and shell, from the table at stack index narg and return
 *  it as one block response function of dimension TB->NFermion (malloc'd, C owns it), or NULL if all entries
 *  are zero. The table has one entry per atom of the tight binding model, which is either
 *    -- a response function of the dimension of the atom (all its shells)
 *    -- a table with one entry per shell, each a response function of the dimension of the shell, or 0
 *    -- the number 0, a zero self energy
 *  All response functions must be of the same type and have the same mu. Single valued ones are copied
 *  into block form first. The result is the direct sum, see ResponseFunctionDirectSum.
 */
static ResponsefunctionType * LuaReadSelfEnergyPerAtom(lua_State *L, int const narg, TightBindingDef const * const TB)
{
  uint32_t i, iatom, maxentry=0, NAtoms=TB->Cell.NAtoms;
  int ishell, idx;
  ResponsefunctionType *rf, *Sigma=NULL;
  SelfEnergyEntriesType E;
  if(lua_objlen(L, narg)!=NAtoms)
  {
    luaL_error(L, "Error CalculateG: the self energy per atom is a table of length %d, but the tight binding model has %d atoms\n", (int)lua_objlen(L, narg), (int)NAtoms);
  }
  for(iatom=0;iatom<NAtoms;iatom++)
  {
    maxentry += (TB->Cell.Atom[iatom].Nshell>0) ? (TB->Cell.Atom[iatom].Nshell) : (1);
  }
  E.N      = 0;
  E.mu     = 0;
  E.Sigma  = (ResponsefunctionType**)malloc(maxentry*sizeof(ResponsefunctionType*));
  E.Offset = (uint32_t*)malloc(maxentry*sizeof(uint32_t));
  E.Atom   = (int*)malloc(maxentry*sizeof(int));
  E.Shell  = (int*)malloc(maxentry*sizeof(int));
  E.Owned  = (bool*)malloc(maxentry*sizeof(bool));
  for(iatom=0;iatom<NAtoms;iatom++)
  {
    lua_rawgeti(L, narg, iatom+1);
    idx = lua_gettop(L);
    if(lua_istable(L, idx) && !LuaIsResponseFunctionObject(L, idx))
    { // a table with one entry per shell
      if(lua_objlen(L, idx)!=(size_t)TB->Cell.Atom[iatom].Nshell)
      {
        luaL_error(L, "Error CalculateG: the self energy of atom %d (%s) is a table of length %d, but the atom has %d shells\n", iatom+1, TB->Cell.Atom[iatom].Name, (int)lua_objlen(L, idx), TB->Cell.Atom[iatom].Nshell);
      }
      for(ishell=0;ishell<TB->Cell.Atom[iatom].Nshell;ishell++)
      {
        lua_rawgeti(L, idx, ishell+1);
        LuaReadSelfEnergyEntry(L, lua_gettop(L), TB, iatom, ishell, &E);
        lua_pop(L, 1);
      }
    }
    else
    {
      LuaReadSelfEnergyEntry(L, idx, TB, iatom, -1, &E);
    }
    lua_pop(L, 1);
  }
  if(E.N!=0)
  {
    // single valued response functions are copied into block form of BlockSize 1
    for(i=0;i<E.N;i++)
    {
      if(E.Sigma[i]->Block==true) continue;
      rf = (ResponsefunctionType*)malloc(sizeof(ResponsefunctionType));
      if(ResponseFunctionToBlockResponseFunction(E.Sigma[i], rf))
      {
        luaL_error(L, "Error CalculateG: could not copy the self energy of atom %d (%s) into block form\n", E.Atom[i]+1, TB->Cell.Atom[E.Atom[i]].Name);
      }
      if(E.Owned[i]) FreeResponsefunctionIncludingOriginalPointer(E.Sigma[i]);
      E.Sigma[i] = rf;
      E.Owned[i] = true;
    }
    Sigma = (ResponsefunctionType*)malloc(sizeof(ResponsefunctionType));
    if(ResponseFunctionDirectSum(E.Sigma, E.Offset, E.N, TB->NFermion, Sigma))
    {
      luaL_error(L, "Error CalculateG: could not merge the self energies per atom or shell into one\n");
    }
    for(i=0;i<E.N;i++)
    {
      if(E.Owned[i]) FreeResponsefunctionIncludingOriginalPointer(E.Sigma[i]);
    }
  }
  free(E.Sigma);
  free(E.Offset);
  free(E.Atom);
  free(E.Shell);
  free(E.Owned);
  return Sigma;
}

/*  Read the self energy at stack index idx in any of the forms CalculateG accepts (a response function,
 *  or a table per atom or per atom and shell, see LuaReadSelfEnergyPerAtom) and return it embedded as
 *  the hermitian matrix
 *       ( a0    C     )
 *   E = ( C^+   Hbath )
 *  whose first block has dimension TB->NFermion (see BlockTriDiagonalMatrixInversedToCompactMatrix and
 *  its siblings). A zero self energy gives the zero matrix of dimension TB->NFermion.
 */
static void LuaReadSelfEnergyEmbedded(lua_State *L, int const idx, TightBindingDef const * const TB, CompactMatrixType * const E)
{
  ResponsefunctionType *Sigma, *SigmaBlock;
  bool FreeSigma=false;
  uint32_t BlockSize;
  int bad=0;
  if(lua_istable(L, idx) && !LuaIsResponseFunctionObject(L, idx))
  {
    Sigma = LuaReadSelfEnergyPerAtom(L, idx, TB);
    FreeSigma = (Sigma!=NULL);
  }
  else
  {
    Sigma = luaL_checkResponsefunction(L, idx);
    FreeSigma = ((Sigma->type&32)==0);
    if(Sigma->Block!=true)
    {
      SigmaBlock = (ResponsefunctionType*)malloc(sizeof(ResponsefunctionType));
      if(ResponseFunctionToBlockResponseFunction(Sigma, SigmaBlock))
      {
        luaL_error(L, "Error: could not copy the self energy into block form\n");
      }
      if(FreeSigma) FreeResponsefunctionIncludingOriginalPointer(Sigma);
      Sigma = SigmaBlock;
      FreeSigma = true;
    }
    BlockSize = ResponseFunctionGetBlockSize(Sigma);
    if(BlockSize != (uint32_t)TB->NFermion)
    {
      luaL_error(L, "Error: the BlockSize of the self energy (%d) and NFermion of the tight binding Hamiltonian (%d) need to be equal\n", (int)BlockSize, TB->NFermion);
    }
  }
  if(Sigma==NULL)
  {
    E->Ni=TB->NFermion; E->Nj=TB->NFermion; E->QComplex=0; E->mu=0;
    InitCompactMatrix(E);
    for(int i=0;i<TB->NFermion*TB->NFermion;i++) E->VRe[i]=0;
    return;
  }
  switch(Sigma->type)
  {
    case 'l': case 'L': bad=BlockListOfPolesInversedToCompactMatrix((BlockListOfPoles*)Sigma->data, E); break;
    case 't': case 'T': bad=BlockTriDiagonalMatrixInversedToCompactMatrix((BlockTriDiagonalMatrixType*)Sigma->data, E); break;
    case 'a': case 'A': bad=BlockAndersonMatrixInversedToCompactMatrix((BlockAndersonMatrixType*)Sigma->data, E); break;
    case 'n': case 'N': bad=BlockNaturalImpurityOrbitalMatrixInversedToCompactMatrix((BlockNaturalImpurityOrbitalMatrixType*)Sigma->data, E); break;
    case 'd': case 'D': bad=BlockDoubleTriDiagonalMatrixInversedToCompactMatrix((BlockDoubleTriDiagonalMatrixType*)Sigma->data, E); break;
    default: bad=1;
  }
  if(FreeSigma) FreeResponsefunctionIncludingOriginalPointer(Sigma);
  if(bad) luaL_error(L, "Error: could not embed the self energy as a matrix\n");
}

/*  Return the character at stack index iChar, a function that takes an eigenstate or a matrix
 *  (Matrix_Type), as the complex n by n matrix P with character(psi) = Re(sum_ij psi_i^* P_ij psi_j).
 *  A function is evaluated on e_i, (e_i+e_j)/sqrt(2) and (e_i+I e_j)/sqrt(2), which gives P by
 *  polarisation, and then on two more vectors to check that it is such a quadratic form. Only then
 *  has the character of a spectral function, Tr[P A(k,w)], a meaning.
 */
static void LuaCharacterToMatrix(lua_State *L, int const iChar, int const n, int const c, CompactMatrixType * const P)
{
  CompactMatrixType *Pin = (CompactMatrixType*)luaL_testudata(L, iChar, "Matrix_Type");
  int i, j, r, iR, NRow = n*n+2;
  double f[NRow], ur[n], ui[n], norm, diff, scale=1;
  P->Ni=n; P->Nj=n; P->QComplex=1;
  InitCompactMatrix(P);
  sprintf(P->Name, "Character %d", c+1);
  if(Pin!=NULL)
  {
    if(Pin->Ni!=n || Pin->Nj!=n)
    {
      luaL_error(L, "Error CalculateBands: character %d is a %d by %d matrix, it should be %d by %d\n", c+1, Pin->Ni, Pin->Nj, n, n);
    }
    for(i=0;i<n*n;i++)
    {
      P->V[i].Re = (Pin->QComplex) ? (Pin->V[i].Re) : (Pin->VRe[i]);
      P->V[i].Im = (Pin->QComplex) ? (Pin->V[i].Im) : (0);
    }
    return;
  }
  if(!lua_isfunction(L, iChar))
  {
    luaL_error(L, "Error CalculateBands: character %d should be a function that takes an eigenstate and returns a number, or a matrix (Matrix.ToUserdata)\n", c+1);
  }
  // the test vectors, as rows of a matrix on the Lua stack
  CompactMatrixType *R = (CompactMatrixType*)lua_newuserdata(L, sizeof(CompactMatrixType));
  R->Ni=NRow; R->Nj=n; R->QComplex=1;
  InitCompactMatrix(R);
  strcpy(R->Name, "Character test vectors");
  R->OnLuaStack=true;
  luaL_setmetatable(L, "Matrix_Type");
  iR = lua_gettop(L);
  for(i=0;i<NRow*n;i++) { R->V[i].Re=0; R->V[i].Im=0; }
  r = n;
  for(i=0;i<n;i++)
  {
    R->V[i*n+i].Re = 1;
    for(j=i+1;j<n;j++)
    {
      R->V[r*n+i].Re = 1/sqrt(2.0); R->V[r*n+j].Re = 1/sqrt(2.0); r++;
      R->V[r*n+i].Re = 1/sqrt(2.0); R->V[r*n+j].Im = 1/sqrt(2.0); r++;
    }
  }
  for(r=n*n;r<NRow;r++)
  {
    norm=0;
    for(i=0;i<n;i++)
    {
      R->V[r*n+i].Re = (r==n*n) ? (cos(1.1*i+0.3)+0.5) : (fmod(i,3.0)-0.8);
      R->V[r*n+i].Im = (r==n*n) ? (sin(0.7*i*i+0.9))   : (0.3*cos(2.3*i+0.1));
      norm += Sqr(R->V[r*n+i].Re) + Sqr(R->V[r*n+i].Im);
    }
    for(i=0;i<n;i++) { R->V[r*n+i].Re/=sqrt(norm); R->V[r*n+i].Im/=sqrt(norm); }
  }
  for(r=0;r<NRow;r++)
  {
    lua_pushvalue(L, iChar);
    lua_createtable(L, 0, 2);
    lua_pushvalue(L, iR);
    lua_setfield(L,-2,"__Matrix");
    lua_pushnumber(L, r+1);
    lua_setfield(L,-2,"__Row");
    lua_getglobal(L, "MatrixRowMeta");
    lua_setmetatable(L, -2);
    if (lua_pcall(L, 1, 1, 0) != LUA_OK)
    {
      const char *err = lua_tostring(L, -1);
      luaL_error(L, "Error in character function %d: %s", c+1, err);
    }
    f[r] = (double)luaL_checknumber(L, -1);
    lua_pop(L,1);
  }
  // polarisation: f((e_i+e_j)/sqrt(2)) = (P_ii+P_jj)/2 + Re P_ij, f((e_i+I e_j)/sqrt(2)) = (P_ii+P_jj)/2 - Im P_ij
  for(i=0;i<n;i++) { P->V[i*n+i].Re = f[i]; P->V[i*n+i].Im = 0; }
  r = n;
  for(i=0;i<n;i++)
  {
    for(j=i+1;j<n;j++)
    {
      P->V[i*n+j].Re =   f[r]   - (f[i]+f[j])/2;
      P->V[i*n+j].Im = -(f[r+1] - (f[i]+f[j])/2);
      P->V[j*n+i].Re =  P->V[i*n+j].Re;
      P->V[j*n+i].Im = -P->V[i*n+j].Im;
      r+=2;
    }
  }
  for(i=0;i<n*n;i++) scale += fabs(P->V[i].Re) + fabs(P->V[i].Im);
  for(r=n*n;r<NRow;r++)
  {
    for(i=0;i<n;i++) { ur[i]=R->V[r*n+i].Re; ui[i]=R->V[r*n+i].Im; }
    diff = f[r] - CharacterQuadraticForm(P, ur, ui, n);
    if(fabs(diff) > 1E-10*scale)
    {
      luaL_error(L, "Error CalculateBands: with a self energy character function %d must be a quadratic form, character(psi) = sum_ij psi_i^* P_ij psi_j, "
                    "as only then Tr[P A(k,w)] is its spectral function. It differs from one by %g on a test vector.\n", c+1, diff);
    }
  }
  lua_pop(L,1);
}

/*  How much memory the poles of the k-points may take before CalculateBands walks the path in
 *  blocks.  One k-point costs (1+NCharacter) N doubles, N the dimension of H(k) plus the
 *  embedded self energy, so 64 MB is 20000 k-points of a one band model with 100 poles of
 *  Sigma, or 2000 of a ten band model with 200.                                            */
static size_t const BandsBufferBytes = 64*1024*1024;

/*  CalculateBands with a self energy. On stack
 *     iSigma) a self energy in any of the forms CalculateG accepts
 *     iChar .. iChar+NCharacter-1) the characters (0, 1 or 3), functions of a quadratic form or matrices, see LuaCharacterToMatrix
 *  At every k of the path the self energy embedded as a hermitian matrix E (see LuaReadSelfEnergyEmbedded)
 *  with H(k) added to its first block is diagonalised. The first HTB->NFermion components u of eigenvector m
 *  give the poles of the band spectral functions:
 *     Tr[A(k,w)]     = sum_m |u_m|^2       delta(w-e_m)
 *     Tr[P_c A(k,w)] = sum_m u_m^dag P_c u_m delta(w-e_m)
 *  H(k) and E are added in the gauge of CalculateG, constant phase within the unit cell, in which a k
 *  independent self energy couples orbitals within the same unit cell. u is then rotated to the gauge of
 *  CalculateBands without self energy, u_i -> e^{I k.r_i} u_i, such that a character means the same with
 *  and without self energy. Poles with a total weight below 1000 DBL_EPSILON are dropped.
 *
 *  We return 1 + NCharacter tables with one scalar list of poles per k-point: the total, and one per character.
 */
static int LuaCalculateBandsSelfEnergy(lua_State *L, TightBindingDef const * const HTB, BandstructurePathType const * const kList, int const Nktot, int const iSigma, int const iChar, int const NCharacter)
{
  int const n = HTB->NFermion;
  int const NSet = 1+NCharacter;
  int i, m, s, a, sh, o, ik, Np, N, iTab;
  CompactMatrixType E, P[3];
  ResponsefunctionType *G;
  ListOfPoles *LOP;
  for(i=0;i<NCharacter;i++) LuaCharacterToMatrix(L, iChar+i, n, i, &(P[i]));
  LuaReadSelfEnergyEmbedded(L, iSigma, HTB, &E);
  N = E.Ni;
  double *rpos = (double*)malloc(3*n*sizeof(double));
  if(rpos==NULL) luaL_error(L, "Error CalculateBands: out of memory\n");
  // position of the atom of every orbital
  for(a=0;a<HTB->Cell.NAtoms;a++)
  {
    for(sh=0;sh<HTB->Cell.Atom[a].Nshell;sh++)
    {
      for(o=0;o<HTB->Cell.Atom[a].shelldegeneracy[sh];o++)
      {
        i = HTB->FermionIndex[a][sh]+o;
        rpos[3*i+0] = HTB->Cell.Atom[a].Position[0];
        rpos[3*i+1] = HTB->Cell.Atom[a].Position[1];
        rpos[3*i+2] = HTB->Cell.Atom[a].Position[2];
      }
    }
  }
  iTab = lua_gettop(L)+1;
  for(s=0;s<NSet;s++) lua_createtable(L, Nktot, 0);
  /*  the k-points of the path, so that the loop over them can be a flat one */
  double *kall = (double*)malloc(3*(size_t)Nktot*sizeof(double));
  if(kall==NULL) luaL_error(L, "Error CalculateBands: out of memory for %d k-points\n", Nktot);
  ik=0;
  for(int iq=0; iq<kList->NQ; iq++)
  {
    for(int jq=0; jq<=kList->Q[iq].Nq; jq++)
    {
      for(i=0;i<3;i++)
      {
        if(jq==0) kall[3*ik+i] = kList->Q[iq].q[i];
        else      kall[3*ik+i] = kList->Q[iq].q[i] * (kList->Q[iq].Nq + 1 - jq) / (kList->Q[iq].Nq + 1) + kList->Q[iq+1].q[i] * jq / (kList->Q[iq].Nq + 1) ;
      }
      ik++;
    }
  }
  /*  The k-points are independent and are diagonalized in parallel into a plain buffer; the
   *  Lua objects are built from that buffer afterwards, as Lua may only be touched by one
   *  thread.  The buffer holds KBlock k-points at a time, so that a long path with a self
   *  energy of many poles does not ask for more memory than the calculation itself.       */
  size_t const PerK = (size_t)(1+NSet)*(size_t)N*sizeof(double);
  int KBlock = (int)(BandsBufferBytes/PerK);
  if(KBlock<1) KBlock=1;
  if(KBlock>Nktot) KBlock=Nktot;
  double *eAll  = (double*)malloc((size_t)KBlock*(size_t)N*sizeof(double));
  double *wAll  = (double*)malloc((size_t)NSet*(size_t)KBlock*(size_t)N*sizeof(double));
  int    *NpAll = (int*)malloc((size_t)KBlock*sizeof(int));
  if((eAll==NULL)||(wAll==NULL)||(NpAll==NULL))
  {
    free(kall); free(eAll); free(wAll); free(NpAll);
    luaL_error(L, "Error CalculateBands: out of memory for the poles of %d k-points of dimension %d\n", KBlock, N);
  }
  for(int kstart=0; kstart<Nktot; kstart+=KBlock)
  {
    int const kend = Min(kstart+KBlock, Nktot);
    int failed=0;
#pragma omp parallel
    {
      double *valt = (double*)malloc(N*sizeof(double));
      double *urt  = (double*)malloc(n*sizeof(double));
      double *uit  = (double*)malloc(n*sizeof(double));
      int ikk;
      if((valt==NULL)||(urt==NULL)||(uit==NULL)) failed=1;
#pragma omp for schedule(static)
      for(ikk=kstart; ikk<kend; ikk++)
      {
        int ii, jj, mm, ss, Np1, lk=ikk-kstart;
        double tr1, ti1, phase1, wt1;
        double const * const kk = kall+3*(size_t)ikk;
        double * const el = eAll+(size_t)lk*(size_t)N;
        CompactMatrixType HS1, T1;
        MatrixType Hcp1;
        if(failed) continue;
        CreateHamiltonianFromTightBindingDefConstantPhaseInUnitCell(*HTB, kk, &Hcp1);
        HS1.Ni=N; HS1.Nj=N; HS1.QComplex=(E.QComplex || Hcp1.QComplex);
        InitCompactMatrix(&HS1);
        if(HS1.QComplex)
        {
          for(ii=0;ii<N*N;ii++)
          {
            HS1.V[ii].Re = (E.QComplex) ? (E.V[ii].Re) : (E.VRe[ii]);
            HS1.V[ii].Im = (E.QComplex) ? (E.V[ii].Im) : (0);
          }
          for(ii=0;ii<n;ii++)
          {
            for(jj=0;jj<n;jj++)
            {
              HS1.V[ii*N+jj].Re += Hcp1.VRe[ii][jj];
              if(Hcp1.QComplex) HS1.V[ii*N+jj].Im += Hcp1.VIm[ii][jj];
            }
          }
        }
        else
        {
          for(ii=0;ii<N*N;ii++) HS1.VRe[ii] = E.VRe[ii];
          for(ii=0;ii<n;ii++) for(jj=0;jj<n;jj++) HS1.VRe[ii*N+jj] += Hcp1.VRe[ii][jj];
        }
        FreeMatrix(&Hcp1);
        CompactMatrixEigensystem(&HS1, &T1, valt);
        Np1=0;
        for(mm=0;mm<N;mm++)
        {
          wt1=0;
          for(ii=0;ii<n;ii++)
          {
            tr1 = (T1.QComplex) ? (T1.V[mm*T1.Nj+ii].Re) : (T1.VRe[mm*T1.Nj+ii]);
            ti1 = (T1.QComplex) ? (T1.V[mm*T1.Nj+ii].Im) : (0);
            double *ri = rpos+3*ii;
            phase1 = Dot(kk, ri);
            urt[ii] = tr1*cos(phase1) - ti1*sin(phase1);
            uit[ii] = tr1*sin(phase1) + ti1*cos(phase1);
            wt1 += Sqr(urt[ii]) + Sqr(uit[ii]);
          }
          if(wt1 < 1000*DBL_EPSILON) continue;
          el[Np1] = valt[mm];
          wAll[(size_t)lk*(size_t)N+Np1] = wt1;
          for(ss=1;ss<NSet;ss++) wAll[((size_t)ss*(size_t)KBlock+(size_t)lk)*(size_t)N+Np1] = CharacterQuadraticForm(&(P[ss-1]), urt, uit, n);
          Np1++;
        }
        NpAll[lk] = Np1;
        FreeCompactMatrix(&T1);
        FreeCompactMatrix(&HS1);
      }
      free(valt); free(urt); free(uit);
    }
    if(failed)
    {
      free(kall); free(eAll); free(wAll); free(NpAll);
      luaL_error(L, "Error CalculateBands: out of memory while diagonalizing the k-points\n");
    }
    /*  and now, one thread, the response function of every k-point of this block */
    for(ik=kstart; ik<kend; ik++)
    {
      int const lk = ik-kstart;
      Np = NpAll[lk];
      for(s=0;s<NSet;s++)
      {
        LOP = (ListOfPoles*)malloc(sizeof(ListOfPoles));
        LOP->N = Np;
        InitListOfPoles(LOP);
        if(s==0) strcpy(LOP->Name, "Band spectral function");
        else     sprintf(LOP->Name, "Band character %d", s);
        for(m=0;m<Np;m++) { LOP->a[m] = eAll[(size_t)lk*(size_t)N+m]; LOP->b[m] = wAll[((size_t)s*(size_t)KBlock+(size_t)lk)*(size_t)N+m]; }
        G = (ResponsefunctionType*)lua_newuserdata(L, sizeof(ResponsefunctionType));
        G->type = 'l';
        G->Block = false;
        G->data = (void*)LOP;
        G->SuggestedGrid = NULL;
        G->NSuggestedGridPoints = 0;
        luaL_setmetatable(L, "ResponseFunction_Type");
        lua_rawseti(L, iTab+s, ik+1);
      }
    }
  }
  free(kall); free(eAll); free(wAll); free(NpAll);
  for(s=0;s<NSet;s++) lua_pushvalue(L, iTab+s);
  FreeCompactMatrix(&E);
  for(i=0;i<NCharacter;i++) FreeCompactMatrix(&(P[i]));
  free(rpos);
  return NSet;
}

/* Calculate the Green's function from a given tight-binding Hamiltonian.
 * Input on stack:
 *   TB    : TightBinding_Type
 *   Sigma : Optional a self energy, either as a ResponseFunctionObject or as a table of Responsefunction objects or as a table of tables of response function objects
 *           A table has one entry per atom: a response function for the atom, a table with a response function per shell, or 0 for a zero self energy (also per shell).
 *           All must be of the same type and have the same mu, see LuaReadSelfEnergyPerAtom.
 *   Options
 *       Emin
 *       Emax
 *       NE
 *       EnergyGrid -> arbirary increasing energy grid for the poles of the response function (in place of Emin, Emax, NE)
 *       Nk -> {Nkx,Nky,Nkz} number of k-points       standard {40,40,40}
 *       Type
 *       HybridizationFunction -> also return the augmented propagator of every atom
 *       Method -> left out: the default route, which diagonalizes H(k) plus the embedded self energy.
 *                 "Mesh": G is evaluated off the real axis and its residues on the grid are taken
 *                 without a broadening (TightBindingHamiltonianPlusSelfEnergyOnMesh). The cost does
 *                 not grow with the number of poles of Sigma. The grid is extended until it holds
 *                 all poles of G, and the result has poles on the grid points only.
 *       Lambda -> Method "Mesh" only: the distances of the points added to extend the grid grow by
 *                 this factor, from the end spacing on (default 2, at least 1; 1 extends uniformly)
 */
static int LuaCalculateG(lua_State *L)
{
  uint32_t len, i;
  TightBindingDef *a;
  ResponsefunctionType *G, *Sigma=NULL;
  uint32_t SigmaBlockSize;
  bool FreeSigma=false;
  struct
  {
    int *Nk;
    char* Method;
    double Lambda;
    double Emin;
    double Emax;
    uint32_t NE;
    const char* Type;
    bool HybridizationFunction;
  } Options;
  Options.Nk = NULL;
  Options.Method = NULL;
  Options.Lambda = NAN;
  bool Mesh = false;
  Options.Emin=NAN;
  Options.Emax=NAN;
  Options.NE=0;
  Options.Type="ListOfPoles";
  Options.HybridizationFunction=false;
  CompactVectorType Grid={0,0,NULL};
  LuaOptionsType LuaCalculateGOptions[]={
    {"Nk", "array-uint32_t-3", &(Options.Nk)},
    {"Method",     "char*",    &(Options.Method)},
    {"Lambda",     "double",   &(Options.Lambda)},
    {"Emin",       "double",   &(Options.Emin)},
    {"Emax",       "double",   &(Options.Emax)},
    {"NE",         "uint32_t", &(Options.NE)},
    {"EnergyGrid", "Vector",   &Grid},
    {"Type",       "char*",    &(Options.Type)},
    {"HybridizationFunction", "bool", &(Options.HybridizationFunction)},
    {NULL, NULL, NULL}};
  len = (uint32_t)lua_gettop(L);
  if (LuaTestOptionsList(L,-1))
  {
    LuaCheckOptions(L,-1,LuaCalculateGOptions);
    len--;
  }
  if(Options.Method!=NULL)
  {
    if(strcicmp(Options.Method,"Mesh")!=0)
    {
      return luaL_error(L, "Error CalculateG: unknown Method \"%s\"; leave it out for the default, or ask for \"Mesh\"\n", Options.Method);
    }
    Mesh = true;
  }
  if(!isnan(Options.Lambda) && !Mesh)
  {
    return luaL_error(L, "Error CalculateG: Lambda belongs to Method \"Mesh\"\n");
  }
  if(isnan(Options.Lambda)) Options.Lambda = 2;
  if(!(Options.Lambda >= 1))
  {
    return luaL_error(L, "Error CalculateG: Lambda is %g, it has to be at least 1\n", Options.Lambda);
  }
  if (Grid.Ni != 0)
  {
    if (!isnan(Options.Emin) || !isnan(Options.Emax) || Options.NE != 0)
      return luaL_error(L, "Error CalculateG: option EnergyGrid cannot be combined with Emin, Emax or NE.\n");
    if (Grid.Ni < 2) return luaL_error(L, "Error CalculateG: option EnergyGrid needs at least 2 energies.\n");
    if (Grid.QComplex) return luaL_error(L, "Error CalculateG: option EnergyGrid needs a real vector.\n");
    for (i=0; (i+1)<Grid.Ni; i++)
    {
      if (!(Grid.VRe[i+1] > Grid.VRe[i]))
      {
        return luaL_error(L, "Error CalculateG: EnergyGrid must be strictly increasing, but element %d (%f) is not larger than element %d (%f).\n",
                          (int)i+2, Grid.VRe[i+1], (int)i+1, Grid.VRe[i]);
      }
    }
  }
  else
  {
    if (isnan(Options.Emin)) Options.Emin = -10;
    if (isnan(Options.Emax)) Options.Emax =  10;
    if (Options.NE == 0)     Options.NE   =  2000;
    if (Options.Emax <= Options.Emin)
      return luaL_error(L, "Error CalculateG: Emax needs to be larger than Emin.\n");
  }
  if(!(len==1||len==2))
  {
    return luaL_error(L, "Error CalculateG: called with %d arguments and possible options, 1 or 2 expected (TightBindingDef, and optional a self energy).\n",len);
  }
  a = (TightBindingDef*) luaL_checkudata(L, 1, "TightBinding_Type");
  if(len==2)
  {
    if(lua_istable(L, 2) && !LuaIsResponseFunctionObject(L, 2))
    { // a self energy per atom or per atom and shell, merged into one
      Sigma = LuaReadSelfEnergyPerAtom(L, 2, a);
      FreeSigma = (Sigma!=NULL);
    }
    else
    {
      Sigma = luaL_checkResponsefunction(L, 2);
      if(Sigma->Block!=true)
      {
        // a single valued self energy is copied into the block form of the same type with BlockSize 1
        ResponsefunctionType *SigmaBlock = (ResponsefunctionType*)malloc(sizeof(ResponsefunctionType));
        if(ResponseFunctionToBlockResponseFunction(Sigma, SigmaBlock))
        {
          return luaL_error(L, "Error CalculateG: could not copy the self energy into block form\n");
        }
        if(FreeSigma==true) FreeResponsefunctionIncludingOriginalPointer(Sigma);
        Sigma = SigmaBlock;
        FreeSigma = true;
      }
      switch(Sigma->type)
      {
        case 'l': case 'L': SigmaBlockSize = ((BlockListOfPoles*)Sigma->data)->BlockSize; break;
        case 't': case 'T': SigmaBlockSize = ((BlockTriDiagonalMatrixType*)Sigma->data)->BlockSize; break;
        case 'a': case 'A': SigmaBlockSize = ((BlockAndersonMatrixType*)Sigma->data)->BlockSize; break;
        case 'n': case 'N': SigmaBlockSize = ((BlockNaturalImpurityOrbitalMatrixType*)Sigma->data)->BlockSize; break;
        case 'd': case 'D': SigmaBlockSize = ((BlockDoubleTriDiagonalMatrixType*)Sigma->data)->BlockSize; break;
        default: return luaL_error(L, "Error CalculateG: self energy of unknown type %c\n", Sigma->type);
      }
      if(SigmaBlockSize != a->NFermion)
      {
        return luaL_error(L, "Error BlockSize of Sigma (%d) and NFermion of TB Hamiltonian (%d) need to be equivalent ", SigmaBlockSize, a->NFermion);
      }
    }
  }
  if(Options.Nk==NULL)
  {
    Options.Nk = (int*)malloc(3*sizeof(int));
    Options.Nk[0]=40;
    Options.Nk[1]=40;
    Options.Nk[2]=40;
  }
  if(Options.Nk[0]<1) return luaL_error(L, "Error first number of k-points needs to be larger than 0");
  if(Options.Nk[1]<1) return luaL_error(L, "Error second number of k-points needs to be larger than 0");
  if(Options.Nk[2]<1) return luaL_error(L, "Error third number of k-points needs to be larger than 0");
  /*  The arguments stay on the stack until everything that points into them has been
   *  read: Sigma is the memory of a response function userdata, and Options.Type and
   *  Grid point into the option table. Lua takes the top nresults values as the
   *  result, so the output is simply pushed on top of them. */
  /*  Asking for the hybridisation routes everything through the grid form, which fills the
   *  augmented propagators from the same eigenvectors that build G. A uniform grid there
   *  gives the same G as Emin, Emax and NE do, to 1E-13. */
  BlockListOfPoles *Gaug=NULL;
  uint32_t NGaug=0;
  double *hybgrid=NULL; uint32_t hybngrid=0;
  CompactMatrixType Ehyb;
  if(Options.HybridizationFunction||Mesh)
  {
    if(Grid.Ni!=0) { hybgrid=Grid.VRe; hybngrid=Grid.Ni; }
    else
    {
      hybngrid = Options.NE+1;
      hybgrid  = (double*)malloc(hybngrid*sizeof(double));
      for(i=0;i<hybngrid;i++) hybgrid[i] = Options.Emin + i*(Options.Emax-Options.Emin)/((double)Options.NE);
    }
  }
  if(Options.HybridizationFunction&&!Mesh)
  {
    if(Sigma==NULL)
    {
      Ehyb.Ni=a->NFermion; Ehyb.Nj=a->NFermion; Ehyb.QComplex=0; Ehyb.mu=0;
      InitCompactMatrix(&Ehyb);
      for(i=0;i<(uint32_t)(a->NFermion*a->NFermion);i++) Ehyb.VRe[i]=0;
    }
    else
    {
      int bad=0;
      switch(Sigma->type)
      {
        case 'l': case 'L': bad=BlockListOfPolesInversedToCompactMatrix((BlockListOfPoles*)Sigma->data, &Ehyb); break;
        case 't': case 'T': bad=BlockTriDiagonalMatrixInversedToCompactMatrix((BlockTriDiagonalMatrixType*)Sigma->data, &Ehyb); break;
        case 'a': case 'A': bad=BlockAndersonMatrixInversedToCompactMatrix((BlockAndersonMatrixType*)Sigma->data, &Ehyb); break;
        case 'n': case 'N': bad=BlockNaturalImpurityOrbitalMatrixInversedToCompactMatrix((BlockNaturalImpurityOrbitalMatrixType*)Sigma->data, &Ehyb); break;
        case 'd': case 'D': bad=BlockDoubleTriDiagonalMatrixInversedToCompactMatrix((BlockDoubleTriDiagonalMatrixType*)Sigma->data, &Ehyb); break;
        default: bad=1;
      }
      if(bad) return luaL_error(L, "Error CalculateG: could not embed the self energy for the hybridization function\n");
    }
  }
  if(Options.HybridizationFunction)
  {
    NGaug = a->Cell.NAtoms;
    Gaug  = (BlockListOfPoles*)malloc(NGaug*sizeof(BlockListOfPoles));
  }
  G = (ResponsefunctionType *) lua_newuserdata(L, sizeof(ResponsefunctionType));
  luaL_setmetatable(L, "ResponseFunction_Type");
  G->NSuggestedGridPoints = 0;
  G->type = 'l';
  G->Block = true;
  G->data = (void*) malloc(sizeof(BlockListOfPoles));
  if(Mesh)
  {
    /*  the same grid as the default route, the grid points are the poles     */
    if(TightBindingHamiltonianPlusSelfEnergyOnMesh(a, Options.Nk, Sigma, hybgrid, hybngrid, Options.Lambda, (BlockListOfPoles*)G->data, Gaug))
    {
      return luaL_error(L, "Error CalculateG with Method \"Mesh\"\n");
    }
    if(Grid.Ni==0) free(hybgrid);
  }
  else if(Options.HybridizationFunction)
  {
    if(TightBindingHamiltonianPlusEmbeddedSelfEnergyToBlockListOfPolesGrid(*a, Options.Nk, &Ehyb, (BlockListOfPoles*)G->data, hybgrid, hybngrid, Gaug))
    {
      return luaL_error(L, "Error CalculateG: could not calculate the hybridization function\n");
    }
    FreeCompactMatrix(&Ehyb);
    if(Grid.Ni==0) free(hybgrid);
  }
  else if(Sigma==NULL)
  {
    if (Grid.Ni==0)
    {
      TightBindingHamiltonianToBlockListOfPoles(*a, Options.Nk, (BlockListOfPoles*)G->data, Options.Emin, Options.Emax, Options.NE);
    }
    else
    {
      TightBindingHamiltonianToBlockListOfPolesGrid(*a, Options.Nk, (BlockListOfPoles*)G->data, Grid.VRe, Grid.Ni);
    }
  }
  else
  {
    BlockListOfPoles *GBLP = (BlockListOfPoles*)G->data;
    switch(Sigma->type)
    {
      case 'l': case 'L':
        if (Grid.Ni==0) TightBindingHamiltonianPlusSelfEnergyToBlockListOfPoles(*a, Options.Nk, *(BlockListOfPoles*)Sigma->data, GBLP, Options.Emin, Options.Emax, Options.NE);
        else            TightBindingHamiltonianPlusSelfEnergyToBlockListOfPolesGrid(*a, Options.Nk, *(BlockListOfPoles*)Sigma->data, GBLP, Grid.VRe, Grid.Ni);
        break;
      case 't': case 'T':
        if (Grid.Ni==0) TightBindingHamiltonianPlusSelfEnergyToBlockListOfPoles(*a, Options.Nk, *(BlockTriDiagonalMatrixType*)Sigma->data, GBLP, Options.Emin, Options.Emax, Options.NE);
        else            TightBindingHamiltonianPlusSelfEnergyToBlockListOfPolesGrid(*a, Options.Nk, *(BlockTriDiagonalMatrixType*)Sigma->data, GBLP, Grid.VRe, Grid.Ni);
        break;
      case 'a': case 'A':
        if (Grid.Ni==0) TightBindingHamiltonianPlusSelfEnergyToBlockListOfPoles(*a, Options.Nk, *(BlockAndersonMatrixType*)Sigma->data, GBLP, Options.Emin, Options.Emax, Options.NE);
        else            TightBindingHamiltonianPlusSelfEnergyToBlockListOfPolesGrid(*a, Options.Nk, *(BlockAndersonMatrixType*)Sigma->data, GBLP, Grid.VRe, Grid.Ni);
        break;
      case 'n': case 'N':
        if (Grid.Ni==0) TightBindingHamiltonianPlusSelfEnergyToBlockListOfPoles(*a, Options.Nk, *(BlockNaturalImpurityOrbitalMatrixType*)Sigma->data, GBLP, Options.Emin, Options.Emax, Options.NE);
        else            TightBindingHamiltonianPlusSelfEnergyToBlockListOfPolesGrid(*a, Options.Nk, *(BlockNaturalImpurityOrbitalMatrixType*)Sigma->data, GBLP, Grid.VRe, Grid.Ni);
        break;
      case 'd': case 'D':
        if (Grid.Ni==0) TightBindingHamiltonianPlusSelfEnergyToBlockListOfPoles(*a, Options.Nk, *(BlockDoubleTriDiagonalMatrixType*)Sigma->data, GBLP, Options.Emin, Options.Emax, Options.NE);
        else            TightBindingHamiltonianPlusSelfEnergyToBlockListOfPolesGrid(*a, Options.Nk, *(BlockDoubleTriDiagonalMatrixType*)Sigma->data, GBLP, Grid.VRe, Grid.Ni);
        break;
    }
  }
  ChopBlockListOfPoles((BlockListOfPoles*)G->data, 1000*DBL_EPSILON);
  free(Options.Nk);
  if(FreeSigma==true)
  {
    FreeResponsefunctionIncludingOriginalPointer(Sigma);
  }
  if(strcmp(Options.Type, "ListOfPoles")!=0)
  {
    /*  ChangeType reads a response function at index 1 and a type at index 2 and
     *  refuses any other number of arguments, so it is left those two and nothing
     *  else. Both are pushed before the arguments under them are dropped. */
    int iG = lua_gettop(L);
    lua_pushstring(L, Options.Type);
    lua_pushvalue(L, iG);
    lua_replace(L, 1);
    lua_replace(L, 2);
    lua_settop(L, 2);
    LuaChangeTypeResponseFunction(L);
  }
  if(Gaug!=NULL)
  {
    lua_newtable(L);
    for(i=0;i<NGaug;i++)
    {
      ResponsefunctionType *R = (ResponsefunctionType*) lua_newuserdata(L, sizeof(ResponsefunctionType));
      luaL_setmetatable(L, "ResponseFunction_Type");
      R->type  = 'l';
      R->Block = true;
      R->NSuggestedGridPoints = 0;
      R->SuggestedGrid = NULL;
      R->data  = (void*)malloc(sizeof(BlockListOfPoles));
      *((BlockListOfPoles*)R->data) = Gaug[i];
      lua_rawseti(L, -2, i+1);
    }
    free(Gaug);
    return 2;
  }
  return 1;
}

static int LuaCalculateMu(lua_State *L)
{
  uint32_t len;
  TightBindingDef *A;
  len = (uint32_t)lua_gettop(L);
  A = (TightBindingDef*) luaL_checkudata(L, 1, "TightBinding_Type");
  int NStateEF, NoccEF;
  double mu, Emin, Emax, gap;
  struct
  {
    int *Nk;
    double nd;
    double kBT;
    //char* Method;   //only method implemented is the mesh/grid.
  } Options;
  Options.Nk = (int*)malloc(3*sizeof(int));
  Options.Nk[0]=40;
  Options.Nk[1]=40;
  Options.Nk[2]=40;
  Options.nd =0.5 * A->NFermion;
  Options.kBT=1E-10;
  LuaOptionsType LuaCalculateMuOptions[]={
    {"Nk", "array-uint32_t-3", &(Options.Nk)},
    {"n", "double"           , &(Options.nd)},
    {"kbT","double"          , &(Options.kBT)},
    {"kBT","double"          , &(Options.kBT)},
    {NULL, NULL, NULL}};
  if (LuaTestOptionsList(L,len))
  {
    LuaCheckOptions(L,len,LuaCalculateMuOptions);
    len--;
  }
  TightBindingDetermineMu(A, Options.Nk, Options.nd, Options.kBT, &mu, &Emin, &Emax, &gap, &NStateEF, &NoccEF);
  lua_pushnumber(L,mu);
  lua_pushnumber(L,Emin);
  lua_pushnumber(L,Emax);
  lua_pushnumber(L,gap);
  lua_pushnumber(L,NStateEF);
  lua_pushnumber(L,NoccEF);
  return 6;
}



/* Calculate the Density matrix from a given tight-binding Hamiltonian or from a given
 * Green's function represented as a response function.
 * Input on stack:
 *   TB : TightBinding_Type
 *   Options
 *       Nk -> {Nkx,Nky,Nkz} number of k-points       standard {40,40,40}
 *  or
 *   G  : ResponseFunction_Type
 *   Options (see LuaResponseFunctionCalculateRho)
 */
static int LuaCalculateRho(lua_State *L)
{
  uint32_t len;
  TightBindingDef *a;
  CompactMatrixType *rho;

  /* a Green's function is handled by the response function version of CalculateRho */
  if(luaL_testudata(L, 1, "TightBinding_Type") == NULL)
  {
    return LuaResponseFunctionCalculateRho(L);
  }

  struct
  {
    int *Nk;
    double mu;
    double kbT;
    //char* Method;   //only method implemented is the mesh/grid.
  } Options;
  Options.Nk = (int*)malloc(3*sizeof(int));
  Options.Nk[0]=40;
  Options.Nk[1]=40;
  Options.Nk[2]=40;
  Options.mu =0;
  Options.kbT=DBL_EPSILON;
  //Options.Method    = "Mesh";
  LuaOptionsType LuaCalculateGOptions[]={
    {"Nk", "array-uint32_t-3", &(Options.Nk)},
    {"mu", "double"          , &(Options.mu)},
    {"kbT","double"          , &(Options.kbT)},
    {"kBT","double"          , &(Options.kbT)},
    //{"Method",     "char*",    &(Options.Method)},
    {NULL, NULL, NULL}};
  
  len = (uint32_t)lua_gettop(L);
  if(!(len==1||len==2))
  {
    return luaL_error(L, "Error CalculateRho: called with %d arguments, 1 or 2 expected (TightBindingDef, and options).\n",len);
  }
  a = (TightBindingDef*) luaL_checkudata(L, 1, "TightBinding_Type");
  // if the reciprocal lattice is zero set Nk to zero
  if((Abs(a->Cell.gUnitCell[0][0])+Abs(a->Cell.gUnitCell[0][1])+Abs(a->Cell.gUnitCell[0][2]))==0) Options.Nk[0]=1;
  if((Abs(a->Cell.gUnitCell[1][0])+Abs(a->Cell.gUnitCell[1][1])+Abs(a->Cell.gUnitCell[1][2]))==0) Options.Nk[1]=1;
  if((Abs(a->Cell.gUnitCell[2][0])+Abs(a->Cell.gUnitCell[2][1])+Abs(a->Cell.gUnitCell[2][2]))==0) Options.Nk[2]=1;
  if(len==2)
  {
    LuaCheckOptions(L, 2, LuaCalculateGOptions);
  }
  
  lua_settop(L, 0);
  rho = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
  luaL_setmetatable(L, "Matrix_Type");
  rho->OnLuaStack=true;
  TightBindingHamiltonianToDensityMatrix(a, Options.Nk, rho, Options.mu, Options.kbT);
  free(Options.Nk);
  
  return 1;
}

static int LuaPushIndicesDictFromTBWithoutSpin(lua_State *L, TightBindingDef *TB)
{
  uint32_t i, j,iAtom, iShell, iOrb, index[TB->NFermion], indexAtom[TB->NFermion], indexShell[TB->NFermion],AtomLength=0,ndigits=1;
  char tempchar[4*MAXStrLength];
  // set print formatting
  if(TB->Cell.NAtoms>999999)
  {
    ndigits=10;
  }
  else if(TB->Cell.NAtoms>99999)
  {
    ndigits=6;
  }
  else if(TB->Cell.NAtoms>9999)
  {
    ndigits=5;
  }
  else if(TB->Cell.NAtoms>999)
  {
    ndigits=4;
  }
  else if(TB->Cell.NAtoms>99)
  {
    ndigits=3;
  }
  else if(TB->Cell.NAtoms>9)
  {
    ndigits=2;
  }
  // All
  j=0;
  for(i=0;i<TB->NFermion;i++)
  {
    index[i]=i;
  }
  lua_createtable(L, TB->NFermion, TB->NFermion);
  lua_pushstring(L, "All");
  lua_rawseti(L,-2,++j);
  lua_pushstring(L, "All");
  LuaPushArray_uint32_t(L, index, TB->NFermion);
  lua_settable(L, -3);
  i=0;
  for(iAtom=0;iAtom<TB->Cell.NAtoms;iAtom++)
  {
    for(iShell=0;iShell<TB->Cell.Atom[iAtom].Nshell;iShell++)
    {
      for(iOrb=0;iOrb<TB->Cell.Atom[iAtom].shelldegeneracy[iShell];iOrb++)
      {
        //AtomName_ShellName_OrbitalName index
        snprintf(tempchar, sizeof(tempchar), "%0*d%s%s%s%s%s%s",ndigits,iAtom,"_",TB->Cell.Atom[iAtom].Name,"_",TB->Cell.Atom[iAtom].ShellName[iShell],"_",TB->Cell.Atom[iAtom].OrbitalName[iShell][iOrb]);
        lua_pushstring(L,tempchar);
        lua_rawseti(L,-2,++j);
        lua_pushstring(L,tempchar);
        LuaPushArray_uint32_t(L,&i,1);
        lua_settable(L, -3);
        indexShell[iOrb]=i;
        indexAtom[AtomLength]=i;
        AtomLength++;
        i++;
      }
      //AtomName_ShellName index
      snprintf(tempchar, sizeof(tempchar), "%0*d%s%s%s%s",ndigits,iAtom,"_",TB->Cell.Atom[iAtom].Name,"_",TB->Cell.Atom[iAtom].ShellName[iShell]);
      lua_pushstring(L,tempchar);
      lua_rawseti(L,-2,++j);
      lua_pushstring(L,tempchar);
      LuaPushArray_uint32_t(L, indexShell, TB->Cell.Atom[iAtom].shelldegeneracy[iShell]);
      lua_settable(L, -3);
    }
    // AtomName index
    snprintf(tempchar, sizeof(tempchar), "%0*d%s%s",ndigits,iAtom,"_",TB->Cell.Atom[iAtom].Name);
    lua_pushstring(L,tempchar);
    lua_rawseti(L,-2,++j);
    lua_pushstring(L,tempchar);
    LuaPushArray_uint32_t(L, indexAtom, AtomLength);
    lua_settable(L, -3);
    AtomLength=0;
  }
  return 1;
}


static int LuaPushIndicesDictFromTBWithSpin(lua_State *L, TightBindingDef *TB)
{
  uint32_t i, j,iAtom, iShell, iOrb, index[2*TB->NFermion], indexAtom[2*TB->NFermion], indexShell[2*TB->NFermion], indexOrb[2], AtomLength=0,ndigits=1;
  char tempchar[4*MAXStrLength];
  // set print formatting
  if(TB->Cell.NAtoms>999999)
  {
    ndigits=10;
  }
  else if(TB->Cell.NAtoms>99999)
  {
    ndigits=6;
  }
  else if(TB->Cell.NAtoms>9999)
  {
    ndigits=5;
  }
  else if(TB->Cell.NAtoms>999)
  {
    ndigits=4;
  }
  else if(TB->Cell.NAtoms>99)
  {
    ndigits=3;
  }
  else if(TB->Cell.NAtoms>9)
  {
    ndigits=2;
  }
  lua_createtable(L, 2*TB->NFermion, 2*TB->NFermion);
  j=0;
  // All
  for(i=0;i<2*TB->NFermion;i++)
  {
    index[i]=i;
  }
  lua_pushstring(L, "All");
  lua_rawseti(L,-2,++j);
  lua_pushstring(L, "All");
  LuaPushArray_uint32_t(L, index, 2*TB->NFermion);
  lua_settable(L, -3);
  for(i=0;i<TB->NFermion;i++)
  {
    index[i]=2*i;
  }
  lua_pushstring(L, "All_Dn");
  lua_rawseti(L,-2,++j);
  lua_pushstring(L, "All_Dn");
  LuaPushArray_uint32_t(L, index, TB->NFermion);
  lua_settable(L, -3);
  for(i=0;i<TB->NFermion;i++)
  {
    index[i]=2*i+1;
  }
  lua_pushstring(L, "All_Up");
  lua_rawseti(L,-2,++j);
  lua_pushstring(L, "All_Up");
  LuaPushArray_uint32_t(L, index, TB->NFermion);
  lua_settable(L, -3);
  i=0;
  for(iAtom=0;iAtom<TB->Cell.NAtoms;iAtom++)
  {
    for(iShell=0;iShell<TB->Cell.Atom[iAtom].Nshell;iShell++)
    {
      for(iOrb=0;iOrb<TB->Cell.Atom[iAtom].shelldegeneracy[iShell];iOrb++)
      {
        //AtomName_ShellName_OrbitalName index
        snprintf(tempchar, sizeof(tempchar), "%0*d%s%s%s%s%s%s",ndigits,iAtom,"_",TB->Cell.Atom[iAtom].Name,"_",TB->Cell.Atom[iAtom].ShellName[iShell],"_",TB->Cell.Atom[iAtom].OrbitalName[iShell][iOrb]);
        lua_pushstring(L,tempchar);
        lua_rawseti(L,-2,++j);
        lua_pushstring(L,tempchar);
        indexOrb[0]=2*i;
        indexOrb[1]=2*i+1;
        LuaPushArray_uint32_t(L,indexOrb,2);
        lua_settable(L, -3);
        indexShell[2*iOrb]=2*i;
        indexShell[2*iOrb+1]=2*i+1;
        indexAtom[2*AtomLength]=2*i;
        indexAtom[2*AtomLength+1]=2*i+1;
        AtomLength++;
        i++;
      }
      //AtomName_ShellName index
      snprintf(tempchar, sizeof(tempchar), "%0*d%s%s%s%s",ndigits,iAtom,"_",TB->Cell.Atom[iAtom].Name,"_",TB->Cell.Atom[iAtom].ShellName[iShell]);
      lua_pushstring(L,tempchar);
      lua_rawseti(L,-2,++j);
      lua_pushstring(L,tempchar);
      LuaPushArray_uint32_t(L, indexShell, 2*TB->Cell.Atom[iAtom].shelldegeneracy[iShell]);
      lua_settable(L, -3);
    }
    // AtomName index
    snprintf(tempchar, sizeof(tempchar), "%0*d%s%s",ndigits,iAtom,"_",TB->Cell.Atom[iAtom].Name);
    lua_pushstring(L,tempchar);
    lua_rawseti(L,-2,++j);
    lua_pushstring(L,tempchar);
    LuaPushArray_uint32_t(L, indexAtom, 2*AtomLength);
    lua_settable(L, -3);
    AtomLength=0;
  }
  
  i=0;
  for(iAtom=0;iAtom<TB->Cell.NAtoms;iAtom++)
  {
    for(iShell=0;iShell<TB->Cell.Atom[iAtom].Nshell;iShell++)
    {
      for(iOrb=0;iOrb<TB->Cell.Atom[iAtom].shelldegeneracy[iShell];iOrb++)
      {
        //AtomName_ShellName_OrbitalName index
        snprintf(tempchar, sizeof(tempchar), "%0*d%s%s%s%s%s%s_Dn",ndigits,iAtom,"_",TB->Cell.Atom[iAtom].Name,"_",TB->Cell.Atom[iAtom].ShellName[iShell],"_",TB->Cell.Atom[iAtom].OrbitalName[iShell][iOrb]);
        lua_pushstring(L,tempchar);
        lua_rawseti(L,-2,++j);
        lua_pushstring(L,tempchar);
        indexOrb[0]=2*i;
        LuaPushArray_uint32_t(L,indexOrb,1);
        lua_settable(L, -3);
        indexShell[iOrb]=2*i;
        indexAtom[AtomLength]=2*i;
        AtomLength++;
        i++;
      }
      //AtomName_ShellName index
      snprintf(tempchar, sizeof(tempchar), "%0*d%s%s%s%s_Dn",ndigits,iAtom,"_",TB->Cell.Atom[iAtom].Name,"_",TB->Cell.Atom[iAtom].ShellName[iShell]);
      lua_pushstring(L,tempchar);
      lua_rawseti(L,-2,++j);
      lua_pushstring(L,tempchar);
      LuaPushArray_uint32_t(L, indexShell, TB->Cell.Atom[iAtom].shelldegeneracy[iShell]);
      lua_settable(L, -3);
    }
    // AtomName index
    snprintf(tempchar, sizeof(tempchar), "%0*d%s%s_Dn",ndigits,iAtom,"_",TB->Cell.Atom[iAtom].Name);
    lua_pushstring(L,tempchar);
    lua_rawseti(L,-2,++j);
    lua_pushstring(L,tempchar);
    LuaPushArray_uint32_t(L, indexAtom, AtomLength);
    lua_settable(L, -3);
    AtomLength=0;
  }
  
  i=0;
  for(iAtom=0;iAtom<TB->Cell.NAtoms;iAtom++)
  {
    for(iShell=0;iShell<TB->Cell.Atom[iAtom].Nshell;iShell++)
    {
      for(iOrb=0;iOrb<TB->Cell.Atom[iAtom].shelldegeneracy[iShell];iOrb++)
      {
        //AtomName_ShellName_OrbitalName index
        snprintf(tempchar, sizeof(tempchar), "%0*d%s%s%s%s%s%s_Up",ndigits,iAtom,"_",TB->Cell.Atom[iAtom].Name,"_",TB->Cell.Atom[iAtom].ShellName[iShell],"_",TB->Cell.Atom[iAtom].OrbitalName[iShell][iOrb]);
        lua_pushstring(L,tempchar);
        lua_rawseti(L,-2,++j);
        lua_pushstring(L,tempchar);
        indexOrb[0]=2*i+1;
        LuaPushArray_uint32_t(L,indexOrb,1);
        lua_settable(L, -3);
        indexShell[iOrb]=2*i+1;
        indexAtom[AtomLength]=2*i+1;
        AtomLength++;
        i++;
      }
      //AtomName_ShellName index
      snprintf(tempchar, sizeof(tempchar), "%0*d%s%s%s%s_Up",ndigits,iAtom,"_",TB->Cell.Atom[iAtom].Name,"_",TB->Cell.Atom[iAtom].ShellName[iShell]);
      lua_pushstring(L,tempchar);
      lua_rawseti(L,-2,++j);
      lua_pushstring(L,tempchar);
      LuaPushArray_uint32_t(L, indexShell, TB->Cell.Atom[iAtom].shelldegeneracy[iShell]);
      lua_settable(L, -3);
    }
    // AtomName index
    snprintf(tempchar, sizeof(tempchar), "%0*d%s%s_Up",ndigits,iAtom,"_",TB->Cell.Atom[iAtom].Name);
    lua_pushstring(L,tempchar);
    lua_rawseti(L,-2,++j);
    lua_pushstring(L,tempchar);
    LuaPushArray_uint32_t(L, indexAtom, AtomLength);
    lua_settable(L, -3);
    AtomLength=0;
  }
  return 1;
}


static int LuaPushIndicesDictFromTB(lua_State *L, TightBindingDef *TB, bool AddSpin)
{
  if(AddSpin==true)
  {
    return LuaPushIndicesDictFromTBWithSpin(L, TB);
  }
  else
  {
    return LuaPushIndicesDictFromTBWithoutSpin(L, TB);
  }
}


/*  Creates an Hamiltonian representing a real space cluster
 *  Input on stack
 *    A tight binding Definition
 *    A cluster definition
 *       Possible formats:
 *         {"periodic", {a,b,c}}
 *         {"open", atoms_to_include}
 *            atoms_to_include can have the following format
 *              {Atom.Name, Atom.Position, {Atom.Shells}}   with Atom.Shells is an option. If not given all shells of this atom are included. (same for orbitals per shell)
 *    Aditional Operator to add (optional)
 *         Fermion index should follow the definition in the original tight binding Hamiltonian, when spin is added it is added in the form down/up, down/up, ...
 *    Options (optional)
 */
static int LuaCreateClusterHamiltonian(lua_State *L)
{
  MatrixType **AtomShellOrbitalMap;
  uint32_t len, lenclusterdef,i,j,ishell,jshell,iorb,jorb,l,nreturn=1;
  char *tokenread;
  TightBindingDef *a, aSuper, *aSuperReturn;
  CrystalStructureSuperCellType SuperCell;
  uint32_t periodic;
  OperatorType *H;
  CompactVectorType k;
  bool AddSpin=false;
  bool QOperators=false, QOptions=false, ReturnTBSuperCell=false, ReturnTBCellPositions=false, ReturnTBKVectors=false, ReturnTBIndicesDict=false;
  OperatorType *AditionalOperator=NULL;
  MatrixType U;
  double AtomPositionAccuracy=0.001;
  k.Ni=3;
  k.QComplex=0;
  k.VRe= (double*)calloc(3,sizeof(double));
  LuaOptionsType LuaCreateClusterHamiltonianOptions[]={
    {"AddSpin", "bool", &AddSpin},
    {"ReturnTBSuperCell", "bool", &ReturnTBSuperCell},
    {"ReturnTBCellPositions", "bool", &ReturnTBCellPositions},
    {"ReturnTBKVectors", "bool", &ReturnTBKVectors},
    {"ReturnTBIndicesDict", "bool", &ReturnTBIndicesDict},
    {"k", "Vector", &k},
    {NULL, NULL, NULL}};
  /* check number of input parameters */
  len = (uint32_t)lua_gettop(L);
  if(!((len==2)||(len==3)||(len==4)))
  {
    return luaL_error(L, "Error CreateClusterHamiltonian called with %d arguments, 2, 3 or 4 expected\n",len);
  }
  /* read the tight binding defintion from the first place on the lua stack */
  a = (TightBindingDef *) luaL_checkudata(L, 1, "TightBinding_Type");
  /* check the cluster defintion from the second position on the stack */
  if(lua_type(L,2)!=LUA_TTABLE)
  {
    return luaL_error(L, "Error CreateClusterHamiltonian second argument should be a table of the form {\"periodic\", {a,b,c}, (optional)atoms_to_include } or {\"open\", atoms_to_include}\n");
  }
  lenclusterdef = (uint32_t)lua_objlen(L, 2);
  if(!((lenclusterdef==2)||(lenclusterdef==3)))
  {
    return luaL_error(L, "Error CreateClusterHamiltonian second argument should be a table of the form {\"periodic\", {a,b,c}, (optional)atoms_to_include } or {\"open\", atoms_to_include}\n");
  }
  /* read if we have open or periodic boundary conditions */
  lua_rawgeti(L, 2, 1);
  if(lua_type(L, -1)!=LUA_TSTRING)
  {
    return luaL_error(L, "Error in CreateClusterHamiltonian second argument should be a table with the first element either the string \"open\" or \"periodic\"\n");
  }
  tokenread = (char*) luaL_checkstring(L, -1);
  if (strcmp(tokenread, "periodic") == 0)
  {
    periodic=1;
  }
  else if (strcmp(tokenread, "open") == 0)
  {
    periodic=0;
  }
  else
  {
    return luaL_error(L, "Error in CreateClusterHamiltonian second argument should be a table with the first element either the string \"open\" or \"periodic\"\n");
  }
  lua_pop(L,1);
  /* if peridic read the super-cell parameters */
  if(periodic)
  {
    lua_rawgeti(L, 2, 2);
    if((lua_type(L, -1)!=LUA_TTABLE)||((uint32_t)lua_objlen(L, -1)!=3))
    {
      return luaL_error(L, "Error in CreateClusterHamiltonian second argument should be a table with the second element a table of the form {a,b,c}\n");
    }
    for(i=1;i<=3;i++)
    {
      lua_rawgeti(L, -1, i);
      if((lua_type(L, -1)!=LUA_TTABLE)||((uint32_t)lua_objlen(L, -1)!=3))
      {
        return luaL_error(L, "Error in CreateClusterHamiltonian second argument should be a table with the second element a table of the form {a,b,c} and element %d (a,b,c) should be a list of 3 integers\n",i);
      }
      for(j=1;j<=3;j++)
      {
        lua_rawgeti(L, -1, j);
        SuperCell.rSuperCell[i-1][j-1]=(int32_t)luaL_checkinteger(L, -1);
        lua_pop(L,1);
      }
      lua_pop(L,1);
    }
    lua_pop(L,1);
    /* read a list of atoms to include / exclude */
    if(lenclusterdef==3)
    {
      //TODO:
      return luaL_error(L, "Error in CreateClusterHamiltonian including only some of the atoms / shells / orbitals not yet implemented \n");
    }
    else
    {
      SuperCell.Crystal=a->Cell;
    }
    AddAditionalInformationToSuperCell(&SuperCell);
    /* deterimine if there are aditional operators or aditional options */
    if(len==2)
    {
      QOperators=false;
      QOptions=false;
    }
    else
    {
      if(len==4)
      {
        QOperators=true;
        QOptions=true;
      }
      /* read additional operators from the stack */
      if (luaL_testudata(L, 3, "Operator_Type") != NULL)
      {
        QOperators=true;
        AditionalOperator = (OperatorType *) luaL_checkudata(L, 3, "Operator_Type");
        if(len==3)
        {
          QOptions=false;
        }
      }
      else
      {
        QOperators=false;
        QOptions=true;
        if(len==4)
        {
          return luaL_error(L, "Error CreateClusterHamiltonian third argument should be an operator defining aditional interactions within the unit cell\n");
        }
      }
    }
    /* read options */
    if(QOptions==true)
    {
      LuaCheckOptions(L, -1, LuaCreateClusterHamiltonianOptions);
      if(k.Ni!=3)
      {
        return luaL_error(L, "Error CreateClusterHamiltonian optional k-vector should have a length of 3\n");
      }
      if(k.QComplex)
      {
        return luaL_error(L, "Error CreateClusterHamiltonian optional k-vector should be real\n");
      }
    }
    /* create the cluster hamiltonian */
    CreateTightBindingSuperCell(*a, &aSuper, SuperCell);
    H=(OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    luaL_setmetatable(L, "Operator_Type");
    if(AddSpin==true)
    {
      CreateHamiltonianOperatorFromTightBindingDefPlusSpin(aSuper, k.VRe, H);
    }
    else
    {
      CreateHamiltonianOperatorFromTightBindingDef(aSuper, k.VRe, H);
    }
    /* add aditional operators */
    if(QOperators)
    {
      if(AddSpin)
      {
        if(AditionalOperator->NFermion!=(2*a->NFermion)) return luaL_error(L, "Error CreateClusterHamiltonian aditional operator should have the same number of spin-orbitals (%f) as the cluster has (%f)\n",(double)AditionalOperator->NFermion,(double)a->NFermion);
      }
      else
      {
        if(AditionalOperator->NFermion!=a->NFermion) return luaL_error(L, "Error CreateClusterHamiltonian aditional operator should have the same number of spin-orbitals (%f) as the cluster has (%f)\n",(double)AditionalOperator->NFermion,(double)a->NFermion);
      }
      U.Ni=H->NFermion;
      U.Nj=AditionalOperator->NFermion;
      U.QComplex=0;
      InitMatrix(&U);
      for(i=0;i<SuperCell.NCell;i++)
      {
        MatrixScalarMultiply(&U, 0, 0);
        for(j=0;j<AditionalOperator->NFermion;j++)
        {
          U.VRe[j+i*AditionalOperator->NFermion][j]=1;
        }
        //TODO: test if function servived changes made in order in mappedOperatorAddNormalOrder. (i.e. row/column major)
        MappedOperatorAddNormalOrder(H, *AditionalOperator, U, &QDetNotOrderedFcFaBcBa);
      }
    }
    SetQComplexOperator(H);
    /* return TBsupercell if asked for */
    if(ReturnTBSuperCell)
    {
      aSuperReturn=(TightBindingDef*) lua_newuserdata(L, sizeof(TightBindingDef));
      luaL_setmetatable(L, "TightBinding_Type");
      *aSuperReturn=aSuper;
      nreturn++;
    }
    if(ReturnTBCellPositions)
    {
      lua_createtable(L, SuperCell.NCell, 0);
      for(i=0;i<SuperCell.NCell;i++)
      {
        lua_createtable(L, 3, 0);
        lua_pushnumber(L, SuperCell.CellPos[i][0]);
        lua_rawseti(L, -2, 1);
        lua_pushnumber(L, SuperCell.CellPos[i][1]);
        lua_rawseti(L, -2, 2);
        lua_pushnumber(L, SuperCell.CellPos[i][2]);
        lua_rawseti(L, -2, 3);
        lua_rawseti(L, -2, i+1);
      }
      nreturn++;
    }
    if(ReturnTBKVectors)
    {
      lua_createtable(L, SuperCell.NCell, 0);
      for(i=0;i<SuperCell.NCell;i++)
      {
        lua_createtable(L, 3, 0);
        lua_pushnumber(L, SuperCell.kVector[i][0]);
        lua_rawseti(L, -2, 1);
        lua_pushnumber(L, SuperCell.kVector[i][1]);
        lua_rawseti(L, -2, 2);
        lua_pushnumber(L, SuperCell.kVector[i][2]);
        lua_rawseti(L, -2, 3);
        lua_rawseti(L, -2, i+1);
      }
      nreturn++;
    }
    if(ReturnTBIndicesDict)
    {
      LuaPushIndicesDictFromTB(L, &aSuper, AddSpin);
      nreturn++;
    }
    if(ReturnTBSuperCell==false)
    {
      FreeTightBindingDef(&aSuper);
    }
    return nreturn;
  }
  else
  { // open systems, i.e. finite size cluster
    TightBindingDef* ClusterDef;
    int *ClusterAtomPosInTB;
    nreturn=2;
    // place return values on stack
    H=(OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    luaL_setmetatable(L, "Operator_Type");
    //
    H->MaxLength=0;
    H->l=NULL;
    H->Res.N=0;
    ClusterDef=(TightBindingDef*) lua_newuserdata(L, sizeof(TightBindingDef));
    luaL_setmetatable(L, "TightBinding_Type");
    //
    ClusterDef->NFermion=0;
    ClusterDef->FermionIndex=NULL;
    ClusterDef->nH=0;
    ClusterDef->H=NULL;
    // set name
    strcpy(ClusterDef->Name, a->Name);
    // set cell
    ClusterDef->Cell.rUnitCell[0][0]=INFINITY;
    ClusterDef->Cell.rUnitCell[0][1]=0;
    ClusterDef->Cell.rUnitCell[0][2]=0;
    ClusterDef->Cell.rUnitCell[1][0]=0;
    ClusterDef->Cell.rUnitCell[1][1]=INFINITY;
    ClusterDef->Cell.rUnitCell[1][2]=0;
    ClusterDef->Cell.rUnitCell[2][0]=0;
    ClusterDef->Cell.rUnitCell[2][1]=0;
    ClusterDef->Cell.rUnitCell[2][2]=INFINITY;
    ClusterDef->Cell.gUnitCell[0][0]=0;
    ClusterDef->Cell.gUnitCell[0][1]=0;
    ClusterDef->Cell.gUnitCell[0][2]=0;
    ClusterDef->Cell.gUnitCell[1][0]=0;
    ClusterDef->Cell.gUnitCell[1][1]=0;
    ClusterDef->Cell.gUnitCell[1][2]=0;
    ClusterDef->Cell.gUnitCell[2][0]=0;
    ClusterDef->Cell.gUnitCell[2][1]=0;
    ClusterDef->Cell.gUnitCell[2][2]=0;
    ClusterDef->Cell.Units=a->Cell.Units;
    // set atomic units to absolute, as there is no unit cell
    ClusterDef->Cell.Units &= ~(((uint32_t)1)<<3);
    // read moleclar structure
    lua_rawgeti(L, 2, 2);
    if(lua_type(L,-1)!=LUA_TTABLE)
    {
      return luaL_error(L, "Error CreateClusterHamiltonian for open clusters expects a list of atoms to include after the token \"open\"\n");
    }
    // read the specific atoms
    ClusterDef->Cell.NAtoms=(uint32_t)lua_objlen(L, -1);
    ClusterDef->Cell.Atom=(AtomType*)malloc(ClusterDef->Cell.NAtoms*sizeof(AtomType));
    ClusterAtomPosInTB=(int*)malloc(ClusterDef->Cell.NAtoms*sizeof(int));
    for(i=0;i<ClusterDef->Cell.NAtoms;i++)
    {
      lua_rawgeti(L, -1, i+1);
      LuaReadAtomWithDefault(L, -1, &(ClusterDef->Cell.Atom[i]), &(a->Cell),&(ClusterAtomPosInTB[i]));
      lua_pop(L, 1);
    }
    lua_pop(L, 1);
    SetNFermionIndexTightBindingDef(ClusterDef);
    
    
    // read options
    /* deterimine if there are aditional operators or aditional options */
    if(len==2)
    {
      QOperators=false;
      QOptions=false;
    }
    else
    {
      if(len==4)
      {
        QOperators=true;
        QOptions=true;
      }
      /* read additional operators from the stack */
      if (luaL_testudata(L, 3, "Operator_Type") != NULL)
      {
        QOperators=true;
        AditionalOperator = (OperatorType *) luaL_checkudata(L, 3, "Operator_Type");
        if(len==3)
        {
          QOptions=false;
        }
      }
      else
      {
        QOperators=false;
        QOptions=true;
        if(len==4)
        {
          return luaL_error(L, "Error CreateClusterHamiltonian third argument should be an operator defining aditional interactions within the unit cell\n");
        }
      }
    }
    /* read options */
    if(QOptions==true)
    {
      LuaCheckOptions(L, len, LuaCreateClusterHamiltonianOptions);
      if(k.Ni!=3)
      {
        return luaL_error(L, "Error CreateClusterHamiltonian optional k-vector should have a length of 3\n");
      }
    }
    // loop over all atoms and shells for creation and annihilation site
    // search for connecting vector in TD definition
    // set Operator as matrix
    ClusterDef->nH = 0;
		AtomShellOrbitalMap = (MatrixType**) malloc(ClusterDef->Cell.NAtoms * sizeof(MatrixType*));
		for (i = 0; i < ClusterDef->Cell.NAtoms; i++)
		{
			ClusterDef->nH += ClusterDef->Cell.Atom[i].Nshell;
			AtomShellOrbitalMap[i] = (MatrixType*) malloc(ClusterDef->Cell.Atom[i].Nshell * sizeof(MatrixType));
			for (ishell = 0; ishell < ClusterDef->Cell.Atom[i].Nshell; ishell++)
			{
				AtomShellOrbitalMap[i][ishell].QComplex = 0;
				AtomShellOrbitalMap[i][ishell].Ni = ClusterDef->Cell.Atom[i].shelldegeneracy[ishell];
				j = FindAtomPosModulusUnitCell(&(ClusterDef->Cell.Atom[i]), &(a->Cell), AtomPositionAccuracy);
				if (j > a->Cell.NAtoms - 1)
					return luaL_error(L, "Atom %s at position (%f,%f,%f) not found in tight binding definition!\n", ClusterDef->Cell.Atom[i].Name, ClusterDef->Cell.Atom[i].Position[0],
							ClusterDef->Cell.Atom[i].Position[1], ClusterDef->Cell.Atom[i].Position[2]);
				for (jshell = 0; jshell < a->Cell.Atom[j].Nshell; jshell++)
				{
					bool sameShellName = strcmp(ClusterDef->Cell.Atom[i].ShellName[ishell], a->Cell.Atom[j].ShellName[jshell]) == 0;
					if (sameShellName)
					{
						AtomShellOrbitalMap[i][ishell].Nj = a->Cell.Atom[j].shelldegeneracy[jshell];
						InitMatrix(&(AtomShellOrbitalMap[i][ishell]));
						for (iorb = 0; iorb < AtomShellOrbitalMap[i][ishell].Ni; iorb++)
						{
							for (jorb = 0; jorb < AtomShellOrbitalMap[i][ishell].Nj; jorb++)
							{
								if (strcmp(ClusterDef->Cell.Atom[i].OrbitalName[ishell][iorb], a->Cell.Atom[j].OrbitalName[jshell][jorb]) == 0)
								{
									AtomShellOrbitalMap[i][ishell].VRe[iorb][jorb] = 1;
									break;
								}
								else if (jorb == AtomShellOrbitalMap[i][ishell].Nj - 1)
								{
									return luaL_error(L, "Orbital %s of shell %s of atom %s not found in tight binding definition!\n", ClusterDef->Cell.Atom[i].OrbitalName[ishell][iorb],
											ClusterDef->Cell.Atom[i].ShellName[ishell], ClusterDef->Cell.Atom[i].Name);
								}
							}
						}
						break;
					}
					else if (jshell == a->Cell.Atom[j].Nshell - 1)
					{
						return luaL_error(L, "Shell %s of atom %s not found in tight binding definition!\n", ClusterDef->Cell.Atom[i].ShellName[ishell], ClusterDef->Cell.Atom[i].Name);
					}
				}
			}
		}
    ClusterDef->nH*=ClusterDef->nH;
    ClusterDef->H=(TightBindingDefHopping*)malloc(ClusterDef->nH * sizeof(TightBindingDefHopping));
    ClusterDef->nH=0;
    for(i=0;i<ClusterDef->Cell.NAtoms;i++)
    {
      for(ishell=0;ishell<ClusterDef->Cell.Atom[i].Nshell;ishell++)
      {
        for(j=0;j<ClusterDef->Cell.NAtoms;j++)
        {
          for(jshell=0;jshell<ClusterDef->Cell.Atom[j].Nshell;jshell++)
          {
            for(l=0;l<a->nH;l++)
            {
              if(  (a->H[l].Atompos1==ClusterAtomPosInTB[i])
                 &&(a->H[l].Atompos2==ClusterAtomPosInTB[j])
                 &&(strcmp(a->Cell.Atom[a->H[l].Atompos1].Name,ClusterDef->Cell.Atom[i].Name)==0)
                 &&(strcmp(a->Cell.Atom[a->H[l].Atompos1].ShellName[a->H[l].Shellpos1],ClusterDef->Cell.Atom[i].ShellName[ishell])==0)
                 &&(strcmp(a->Cell.Atom[a->H[l].Atompos2].Name,ClusterDef->Cell.Atom[j].Name)==0)
                 &&(strcmp(a->Cell.Atom[a->H[l].Atompos2].ShellName[a->H[l].Shellpos2],ClusterDef->Cell.Atom[j].ShellName[jshell])==0)
                 &&(EqualFloat(a->H[l].ConnectingVector[0], ClusterDef->Cell.Atom[j].Position[0]-ClusterDef->Cell.Atom[i].Position[0],AtomPositionAccuracy)==1)
                 &&(EqualFloat(a->H[l].ConnectingVector[1], ClusterDef->Cell.Atom[j].Position[1]-ClusterDef->Cell.Atom[i].Position[1],AtomPositionAccuracy)==1)
                 &&(EqualFloat(a->H[l].ConnectingVector[2], ClusterDef->Cell.Atom[j].Position[2]-ClusterDef->Cell.Atom[i].Position[2],AtomPositionAccuracy)==1) )
              {
                ClusterDef->H[ClusterDef->nH].Atompos1=i;
                ClusterDef->H[ClusterDef->nH].Shellpos1=ishell;
                ClusterDef->H[ClusterDef->nH].Atompos2=j;
                ClusterDef->H[ClusterDef->nH].Shellpos2=jshell;
                ClusterDef->H[ClusterDef->nH].TransVecs[0]=0;
                ClusterDef->H[ClusterDef->nH].TransVecs[1]=0;
                ClusterDef->H[ClusterDef->nH].TransVecs[2]=0;
                ClusterDef->H[ClusterDef->nH].ConnectingVector[0]=ClusterDef->Cell.Atom[j].Position[0]-ClusterDef->Cell.Atom[i].Position[0];
                ClusterDef->H[ClusterDef->nH].ConnectingVector[1]=ClusterDef->Cell.Atom[j].Position[1]-ClusterDef->Cell.Atom[i].Position[1];
                ClusterDef->H[ClusterDef->nH].ConnectingVector[2]=ClusterDef->Cell.Atom[j].Position[2]-ClusterDef->Cell.Atom[i].Position[2];
                MatrixTransformLeftRight(&(ClusterDef->H[ClusterDef->nH].Hop), a->H[l].Hop, AtomShellOrbitalMap[i][ishell], AtomShellOrbitalMap[j][jshell]);
                ClusterDef->nH++;
              }
            }
          }
        }
      }
    }
    // build Hamiltonian from cluster tight binding definition
    if(AddSpin==true)
    {
      CreateHamiltonianOperatorFromTightBindingDefPlusSpin(*ClusterDef, k.VRe, H);
    }
    else
    {
      CreateHamiltonianOperatorFromTightBindingDef(*ClusterDef, k.VRe, H);
    }
    // add additional operator
    if(QOperators==true)
    {
      if(H->NFermion==AditionalOperator->NFermion)
      {
        OperatorAddNormalOrder(H, *AditionalOperator, 1, 0, &QDetNotOrderedFcFaBcBa);
      }
      else
      {
        return luaL_error(L, "Number of Fermions in tight binding Hamiltonian (%i) not equal to number of Fermions in additional operator (%i)\n",H->NFermion,AditionalOperator->NFermion);
      }
    }
    for(i=0;i<ClusterDef->Cell.NAtoms;i++)
    {
      if(AtomShellOrbitalMap[i]!=NULL) free(AtomShellOrbitalMap[i]);
    }
    free(AtomShellOrbitalMap);
    free(ClusterAtomPosInTB);
    
    if(ReturnTBIndicesDict)
    {
      LuaPushIndicesDictFromTB(L, ClusterDef, AddSpin);
      nreturn++;
    }
    
    return nreturn;
  }
}

static int LuaFindAllAtomsInsideSphere(lua_State *L)
{
  char StdUnits[]="Absolute";
  char* Units=StdUnits;
  CrystalStructureType A, B;
  double c[3],r;
  A.Units=8;
  uint32_t len,i,j;
  LuaOptionsType LuaFindAllAtomsInsideSphereOptions[]={
    {"Units","char*",Units},
    {NULL, NULL, NULL}};
  len = (uint32_t)lua_gettop(L);
  if((len!=4)&&(len!=5))
  {
    return luaL_error(L, "Error FindAllAtomsInsideSphere expect 4 or 5 arguments and not %d\n",len);
  }
  if(lua_type(L,1)!=LUA_TTABLE)
  {
    return luaL_error(L, "Error FindAllAtomsInsideSphere expect a list of atoms on position 1\n");
  }
  if(lua_type(L,2)!=LUA_TTABLE)
  {
    return luaL_error(L, "Error FindAllAtomsInsideSphere expects a table (3 by 3) defining the unit cell on position 2\n");
  }
  if(lua_type(L,3)!=LUA_TTABLE)
  {
    return luaL_error(L, "Error FindAllAtomsInsideSphere expects a table length 3 defining the center of the sphere on position 3\n");
  }
  if(lua_type(L,4)!=LUA_TNUMBER)
  {
    return luaL_error(L, "Error FindAllAtomsInsideSphere expects a number on position 4 defining the sphere radius\n");
  }
  len = (uint32_t)lua_objlen(L, 1);
  A.NAtoms=len;
  A.Atom= (AtomType*)malloc(len*sizeof(AtomType));
  for(i=0;i<len;i++)
  {
    lua_rawgeti(L, 1, i+1);
    LuaReadAtom(L, -1, &(A.Atom[i]));
    lua_pop(L, 1);
  }
  if(lua_objlen(L, 2)!=3)
  {
    return luaL_error(L, "Error FindAllAtomsInsideSphere expects unitcell definition on position 2\n");
  }
  for(i=0;i<3;i++)
  {
    lua_rawgeti(L, 2, i+1);
    if(lua_type(L,-1)!=LUA_TTABLE)
    {
      return luaL_error(L, "Error FindAllAtomsInsideSphere expects unitcell definition on position 2\n");
    }
    for(j=0;j<3;j++)
    {
      lua_rawgeti(L, -1, j+1);
      A.rUnitCell[i][j]=(double)luaL_checknumber(L, -1);
      lua_pop(L,1);
    }
    lua_pop(L, 1);
  }
  len = (uint32_t)lua_gettop(L);
  for(j=0;j<3;j++)
  {
    lua_rawgeti(L, 3, j+1);
    c[j]=(double)luaL_checknumber(L, -1);
    lua_pop(L,1);
  }
  r=(double)luaL_checknumber(L, 4);
  if(len==5)
  {
    LuaCheckOptions(L, 5, LuaFindAllAtomsInsideSphereOptions);
  }
  if(strcmp(Units, "Absolute")==0) A.Units&=~8;
  if(strcmp(Units, "Relative")==0) A.Units|=8;
  FindAllAtomsInsideSphere(&A, c, r, &B);
  lua_createtable(L, B.NAtoms, 0);
  for(i=0;i<B.NAtoms;i++)
  {
    LuaPushAtom(L, &(B.Atom[i]));
    lua_rawseti(L, -2, i+1);
  }
  FreeCrystalStructureType(&A);
  FreeCrystalStructureType(&B);
  return 1;
}


/*
 * reads a atom definition from the Lua stack
 * format:
 *           Atom   = {Name, Position, {Shell1, Shell2, ...}}
 *           Shelli = {Name, {orb1, orb2, ...}}
 */
static int LuaReadAtom(lua_State *L, int narg, AtomType *a)
{
  uint32_t i,j,len;
  double *arr;
  if(lua_type(L,narg)==LUA_TTABLE)
  {
    len = (uint32_t)lua_objlen(L, narg);
    if(len!=3) return luaL_error(L, "Error ReadAtom expects a list of the form {Name, Position, {Shell1, Shell2, ...}}\n");
    // get the name of the atom
    lua_rawgeti(L, narg, 1);
    strncpy(a->Name, luaL_checkstring(L, -1), MAXStrLength-1);
    lua_pop(L, 1);
    // set the atom number
    a->AtomNumber=AtomNameToNumber(a->Name);
    // get the position of the atom
    lua_rawgeti(L, narg, 2);
    arr=luaL_checkarray_double(L,-1,&len);
    if(len!=3) return luaL_error(L, "Error ReadAtom position of atom %s should be a list of 3 numbers\n",a->Name);
    a->Position[0]=arr[0];a->Position[1]=arr[1];a->Position[2]=arr[2];
    free(arr);
    lua_pop(L, 1);
    // get the list with shells
    lua_rawgeti(L, narg, 3);
    if(lua_type(L,-1)!=LUA_TTABLE) return luaL_error(L, "Error ReadAtom shells of atom %s should be a list of shells\n",a->Name);
    a->Nshell = (int)lua_objlen(L, -1);
    a->ShellName = (char**)malloc(a->Nshell*sizeof(char*));
    a->shelldegeneracy = (int*)malloc(a->Nshell*sizeof(int));
    a->OrbitalName = (char***)malloc(a->Nshell*sizeof(char**));
    for(i=0;i<a->Nshell;i++)
    {
      lua_rawgeti(L, -1, i+1);
      if(lua_type(L,-1)!=LUA_TTABLE) return luaL_error(L, "Error ReadAtom shell %d of atom %s should be a list containing the name and a list of orbitals\n",i+1,a->Name);
      len = (uint32_t)lua_objlen(L, -1);
      if(len!=2) return luaL_error(L, "Error ReadAtom shell %d of atom %s should be a list containing the name and a list of orbitals\n",i+1,a->Name);
      lua_rawgeti(L, -1, 1);
      a->ShellName[i]= (char*)malloc(MAXStrLength*sizeof(char));
      strncpy(a->ShellName[i], luaL_checkstring(L, -1), MAXStrLength-1);
      lua_pop(L, 1);
      lua_rawgeti(L, -1, 2);
      if(lua_type(L,-1)!=LUA_TTABLE) return luaL_error(L, "Error ReadAtom shell %s(%d) of atom %s should be a list containing the name and a list of orbitals\n",a->ShellName[i],i+1,a->Name);
      a->shelldegeneracy[i] = (uint32_t)lua_objlen(L, -1);
      a->OrbitalName[i] = (char**)malloc(a->shelldegeneracy[i]*sizeof(char*));
      for(j=0;j<a->shelldegeneracy[i];j++)
      {
        lua_rawgeti(L, -1, j+1);
        a->OrbitalName[i][j]= (char*)malloc(MAXStrLength*sizeof(char));
        strncpy(a->OrbitalName[i][j], luaL_checkstring(L, -1), MAXStrLength-1);
        lua_pop(L, 1);
      }
      lua_pop(L, 1);
      lua_pop(L, 1);
    }
    lua_pop(L, 1);
  }
  else
  {
    return luaL_error(L, "Error ReadAtom expects a list of the form {Name, Position, {Shell1, Shell2, ...}}\n");
  }
  return 0;
}

/*  Read atom name and position from lua stack
 *  From name and position determine which atom one refers to in the CrystalStructure
 *     Name must be equal
 *     if G= 0
 *     (Position_atom - Position_atom_in_structure) = 0
 *     else
 *     (Position_atom - Position_atom_in_structure) . G = n or n pi depending on units
 */
static int LuaReadAtomWithDefault(lua_State *L, int narg, AtomType *a, CrystalStructureType *def, int* AtomIndexInCrystalStructure)
{
  uint32_t i,j,k,l,len;
  double *arr;
  a->Name[0]='\0';
  a->Nshell=0;
  //a->OrbitalName=NULL;
  //a->shelldegeneracy=NULL;
  //a->ShellName=NULL;
  if(lua_type(L,narg)==LUA_TTABLE)
  {
    len = (uint32_t)lua_objlen(L, narg);
    if(!((len==3)||(len==2))) return luaL_error(L, "Error ReadAtom expects a list of the form {Name, Position, {Shell1, Shell2, ...}} or {Name, Position}\n");
    // get the name of the atom
    lua_rawgeti(L, narg, 1);
    strncpy(a->Name, luaL_checkstring(L, -1), MAXStrLength-1);
    lua_pop(L, 1);
    // set the atom number
    a->AtomNumber=AtomNameToNumber(a->Name);
    // get the position of the atom
    lua_rawgeti(L, narg, 2);
    arr=luaL_checkarray_double(L,-1,&len);
    if(len!=3) return luaL_error(L, "Error ReadAtom position of atom %s should be a list of 3 numbers\n",a->Name);
    a->Position[0]=arr[0];a->Position[1]=arr[1];a->Position[2]=arr[2];
    free(arr);
    lua_pop(L, 1);
    // determine AtomIndexInCrystalStructure
    *AtomIndexInCrystalStructure = FindAtomPosModulusUnitCell(a, def, 0.0001);
    if((len==2)&&(*AtomIndexInCrystalStructure==INT_MAX)) luaL_error(L, "Error Atom %s not found at position %f %f %f in tight binding object\n",a->Name,a->Position[0],a->Position[1],a->Position[2]);
    // get the list with shells
    len = (uint32_t)lua_objlen(L, narg);
    if(len==3)
    {
      lua_rawgeti(L, narg, 3);
      if(lua_type(L,-1)!=LUA_TTABLE) return luaL_error(L, "Error ReadAtom shells of atom %s at position %f %f %f, should be a list of shells\n",a->Name,a->Position[0],a->Position[1],a->Position[2]);
      a->Nshell = (int)lua_objlen(L, -1);
      a->ShellName = (char**)malloc(a->Nshell*sizeof(char*));
      a->shelldegeneracy = (int*)malloc(a->Nshell*sizeof(int));
      a->OrbitalName = (char***)malloc(a->Nshell*sizeof(char**));
      for(i=0;i<a->Nshell;i++)
      {
        lua_rawgeti(L, -1, i+1);
        if(lua_type(L,-1)==LUA_TTABLE)
        {
          // shell given by a table of the form {"shellname", {"orb1name", "orb2name", ...}}
          len = (uint32_t)lua_objlen(L, -1);
          if(len!=2) return luaL_error(L, "Error ReadAtom shell %d of atom %s at position %f %f %f should be a list containing the name and a list of orbitals\n",i+1,a->Name,a->Position[0],a->Position[1],a->Position[2]);
          lua_rawgeti(L, -1, 1);
          a->ShellName[i]= (char*)malloc(MAXStrLength*sizeof(char));
          strncpy(a->ShellName[i], luaL_checkstring(L, -1), MAXStrLength-1);
          lua_pop(L, 1);
          lua_rawgeti(L, -1, 2);
          if(lua_type(L,-1)!=LUA_TTABLE) return luaL_error(L, "Error ReadAtom shell %s(%d) of atom %s at position %f %f %f should be a list containing the name and a list of orbitals\n",a->ShellName[i],i+1,a->Name,a->Position[0],a->Position[1],a->Position[2]);
          a->shelldegeneracy[i] = (uint32_t)lua_objlen(L, -1);
          a->OrbitalName[i] = (char**)malloc(a->shelldegeneracy[i]*sizeof(char*));
          for(j=0;j<a->shelldegeneracy[i];j++)
          {
            lua_rawgeti(L, -1, j+1);
            a->OrbitalName[i][j]= (char*)malloc(MAXStrLength*sizeof(char));
            strncpy(a->OrbitalName[i][j], luaL_checkstring(L, -1), MAXStrLength-1);
            lua_pop(L, 1);
          }
          lua_pop(L, 1);
        }
        else
        {
          // shell given by a string specifing its name, all orbitals taken from def
          a->ShellName[i]= (char*)malloc(MAXStrLength*sizeof(char));
          strncpy(a->ShellName[i], luaL_checkstring(L, -1), MAXStrLength-1);
          // find atom and shell in def and set the list of orbitals
          if(*AtomIndexInCrystalStructure!=INT_MAX)
          {
            j=*AtomIndexInCrystalStructure;
            {
              for(k=0;k<def->Atom[j].Nshell;k++)
              {
                if(strcmp(a->ShellName[i], def->Atom[j].ShellName[k])==0)
                {
                  a->shelldegeneracy[i] = def->Atom[j].shelldegeneracy[k];
                  a->OrbitalName[i] = (char**)malloc(a->shelldegeneracy[k]*sizeof(char*));
                  for(l=0;l<a->shelldegeneracy[k];l++)
                  {
                    a->OrbitalName[k][l]= (char*)malloc(MAXStrLength*sizeof(char));
                    strncpy(a->OrbitalName[k][l], def->Atom[j].OrbitalName[k][l], MAXStrLength-1);
                  }
                  goto FoundShell;
                }
              }
              return luaL_error(L, "Found atom with label %s and shell with label %s but without a definition of orbitals\n Shell of atom is not in the refference structure so I can't set the shells and orbitals to default\n", a->Name, a->ShellName[i]);
            }
          }
          return luaL_error(L, "Found atom with label %s but without a definition of shells and orbitals\n Atom is not in the refference structure so I can't set the shells and orbitals to default\n", a->Name);
          FoundShell:;
        }
        lua_pop(L, 1);
      }
      lua_pop(L, 1);
    }
    else
    {
      // find atom in def and set the list with shells
      for(i=0;i<def->NAtoms;i++)
      {
        if(strcmp(a->Name, def->Atom[i].Name)==0)
        {
          a->Nshell = def->Atom[i].Nshell;
          a->ShellName = (char**)malloc(a->Nshell*sizeof(char*));
          a->shelldegeneracy = (int*)malloc(a->Nshell*sizeof(int));
          a->OrbitalName = (char***)malloc(a->Nshell*sizeof(char**));
          for(j=0;j<a->Nshell;j++)
          {
            a->ShellName[j]= (char*)malloc(MAXStrLength*sizeof(char));
            strncpy(a->ShellName[j], def->Atom[i].ShellName[j], MAXStrLength-1);
            a->shelldegeneracy[j] = def->Atom[i].shelldegeneracy[j];
            a->OrbitalName[j] = (char**)malloc(a->shelldegeneracy[j]*sizeof(char*));
            for(k=0;k<a->shelldegeneracy[j];k++)
            {
              a->OrbitalName[j][k]= (char*)malloc(MAXStrLength*sizeof(char));
              strncpy(a->OrbitalName[j][k], def->Atom[i].OrbitalName[j][k], MAXStrLength-1);
            }
          }
          goto FoundAtom;
        }
      }
      return luaL_error(L, "Found atom with label %s at position %f %f %f, but without a definition of shells and orbitals\n Atom is not in the refference structure so I can't set the shells and orbitals to default\n", a->Name,a->Position[0],a->Position[1],a->Position[2]);
      FoundAtom:;
    }
  }
  else
  {
    return luaL_error(L, "Error ReadAtom expects a list of the form {Name, Position, {Shell1, Shell2, ...}} or {Name, Position}\n");
  }
  return 0;
}

/*
 * push a definition of an atom to the stack as list
 * format:
 *           Atom   = {Name, Position, {Shell1, Shell2, ...}}
 *           Shelli = {Name, {orb1, orb2, ...}}
 */
static int LuaPushAtom(lua_State *L, AtomType const *a)
{
  uint32_t i,j;
  lua_createtable(L, 3, 0);
  {
    lua_pushstring(L, a->Name);
    lua_rawseti(L, -2, 1); // set name of atom to first position
    lua_createtable(L, 3, 0);
    {
      for(i=0;i<3;i++)
      {
        lua_pushnumber(L, a->Position[i]);
        lua_rawseti(L, -2, i+1);
      }
    }
    lua_rawseti(L, -2, 2); // set table with atom coordinates to second position
    lua_createtable(L, a->Nshell, 0);
    {
      for(i=0;i<a->Nshell;i++)
      {
        lua_createtable(L,2,0);
        {
          lua_pushstring(L, a->ShellName[i]);
          lua_rawseti(L, -2, 1); // set name of the shell to first position
          lua_createtable(L,a->shelldegeneracy[i],0);
          for(j=0;j<a->shelldegeneracy[i];j++)
          {
            lua_pushstring(L, a->OrbitalName[i][j]);
            lua_rawseti(L, -2, j+1);
          }
          lua_rawseti(L, -2, 2); // set a table with a list of the orbitals to the second position
        }
        lua_rawseti(L, -2, i+1); // set a table with the shell name and orbital definitions to position i+1
      }
    }
    lua_rawseti(L, -2, 3); // set a table with the shell definitions to the third position
  }
  return 1;
}

/*
 * reads a Hopping definition from the Lua stack
 * format
 *        Hopping = {"AtomName1.ShellName1","AtomName2.ShellName2",Vec,H}
 *        with H a table of tables of real or complex numbers / matrix
 *
 */
static int LuaReadHopping(lua_State *L, int narg, CrystalStructureType *Cell, TightBindingDefHopping *H)
{
  double buff,CellVolume,r[3],idouble[3];
  uint32_t i,j;
  char AtomShell1[2*MAXStrLength+1];
  char AtomShell2[2*MAXStrLength+1];
  char *Atom1,*Atom2,*Shell1,*Shell2;
  bool comparepos=false;
  double AtomPositionAccuracy=0.001;
  // test input
  if((lua_type(L,narg)!=LUA_TTABLE)) return luaL_error(L, "ReadHopping expects a table of the form {\"AtomName1.ShellName1\",\"AtomName2.ShellName2\",Vec,H}\n");
  if (lua_objlen(L, narg)!=4) return luaL_error(L, "ReadHopping expects a table with 4 elements of the form {\"AtomName1.ShellName1\",\"AtomName2.ShellName2\",Vec,H}\n");
  // read atom and shell 1
  lua_rawgeti(L, narg, 1);
  if(lua_istable(L, -1))
  {
    comparepos=true;
    lua_rawgeti(L, -1, 1);
    strncpy(AtomShell1, luaL_checkstring(L, -1), 2*MAXStrLength);
    lua_pop(L,1);
    lua_rawgeti(L, -1, 2);
    for(i=0;i<3;i++)
    {
      lua_rawgeti(L, -1, i+1);
      idouble[i]=luaL_checknumber(L, -1);
      lua_pop(L,1);
    }
    lua_pop(L,1);
  }
  else
  {
    strncpy(AtomShell1, luaL_checkstring(L, -1), 2*MAXStrLength);
  }
  lua_pop(L, 1);
  Atom1=AtomShell1;
  Shell1=strchr(AtomShell1 , '.');
  if(Shell1==NULL) return luaL_error(L, "ReadHopping expects a string of the form \"AtomName.ShellName\" defining the first atom and shell in the hopping got %s\n",AtomShell1);
  *Shell1='\0';
  Shell1++;
  H->Atompos1=Cell->NAtoms;
  for(i=0;i<Cell->NAtoms;i++)
  {
    if(strcmp(Cell->Atom[i].Name, Atom1)==0)
    {
      if((comparepos==false)||(EqualFloat(Cell->Atom[i].Position[0], idouble[0],AtomPositionAccuracy)&&EqualFloat(Cell->Atom[i].Position[1], idouble[1],AtomPositionAccuracy)&&EqualFloat(Cell->Atom[i].Position[2], idouble[2],AtomPositionAccuracy)))
      {
        H->Atompos1=i;
        break;
      }
    }
  }
  if(H->Atompos1==Cell->NAtoms) return luaL_error(L, "Atom name %s in ReadHopping not found in Cell definition\n",Atom1);
  H->Shellpos1=Cell->Atom[H->Atompos1].Nshell;
  for(i=0;i<Cell->Atom[H->Atompos1].Nshell;i++)
  {
    if(strcmp(Cell->Atom[H->Atompos1].ShellName[i], Shell1)==0)
    {
      H->Shellpos1=i;
      break;
    }
  }
  if(H->Shellpos1==Cell->Atom[H->Atompos1].Nshell) return luaL_error(L, "Shell name %s of Atom %s in ReadHopping not found in Cell definition\n",Shell1,Atom1);
  // read atom and shell 2
  comparepos=false;
  lua_rawgeti(L, narg, 2);
  if(lua_istable(L, -1))
  {
    comparepos=true;
    lua_rawgeti(L, -1, 1);
    strncpy(AtomShell2, luaL_checkstring(L, -1), 2*MAXStrLength);
    lua_pop(L,1);
    lua_rawgeti(L, -1, 2);
    for(i=0;i<3;i++)
    {
      lua_rawgeti(L, -1, i+1);
      idouble[i]=luaL_checknumber(L, -1);
      lua_pop(L,1);
    }
    lua_pop(L,1);
  }
  else
  {
    strncpy(AtomShell2, luaL_checkstring(L, -1), 2*MAXStrLength);
  }
  lua_pop(L, 1);
  Atom2=AtomShell2;
  Shell2=strchr(AtomShell2 , '.');
  if(Shell2==NULL) return luaL_error(L, "ReadHopping expects a string of the form \"AtomName.ShellName\" defining the second atom and shell in the hopping got %s\n",AtomShell2);
  *Shell2='\0';
  Shell2++;
  H->Atompos2=Cell->NAtoms;
  for(i=0;i<Cell->NAtoms;i++)
  {
    if(strcmp(Cell->Atom[i].Name, Atom2)==0)
    {
      if((comparepos==false)||(EqualFloat(Cell->Atom[i].Position[0], idouble[0],AtomPositionAccuracy)&&EqualFloat(Cell->Atom[i].Position[1], idouble[1],AtomPositionAccuracy)&&EqualFloat(Cell->Atom[i].Position[2], idouble[2],AtomPositionAccuracy)))
      {
        H->Atompos2=i;
        break;
      }
    }
  }
  if(H->Atompos2==Cell->NAtoms) return luaL_error(L, "Atom name %s in ReadHopping not found in Cell definition\n",Atom2);
  H->Shellpos2=Cell->Atom[H->Atompos2].Nshell;
  for(i=0;i<Cell->Atom[H->Atompos2].Nshell;i++)
  {
    if(strcmp(Cell->Atom[H->Atompos2].ShellName[i], Shell2)==0)
    {
      H->Shellpos2=i;
      break;
    }
  }
  if(H->Shellpos2==Cell->Atom[H->Atompos2].Nshell) return luaL_error(L, "Shell name %s of Atom %s in ReadHopping not found in Cell definition\n",Shell2,Atom2);
  // read trans vector
  lua_rawgeti(L, narg, 3);
  if((lua_type(L,-1)!=LUA_TTABLE)) return luaL_error(L, "ReadHopping expects a table of the form {x,y,z} or {i,j,k} for the vector defining the connection between atom %s and atom %s\n",AtomShell1,AtomShell2);
  if (lua_objlen(L, -1)!=3) return luaL_error(L, "ReadHopping expects a table of the form {x,y,z} or {i,j,k} for the vector defining the connection between atom %s and atom %s\n",AtomShell1,AtomShell2);
  if((Cell->Units&8)==0) // absolute units
  {
    for(i=0;i<3;i++)
    {
      lua_rawgeti(L, -1, i+1);
      H->ConnectingVector[i]=(double)luaL_checknumber(L, -1);
      lua_pop(L,1);
    }
    lua_pop(L,1);
    // create TransVec from connecting Vector
    if(NormSqr(Cell->rUnitCell[0])==INFINITY)
    {
      if(NormSqr(Cell->rUnitCell[1])==INFINITY)
      {
        if(NormSqr(Cell->rUnitCell[2])==INFINITY)
        {  // a, b and c infinite
          H->TransVecs[0]=0;
          H->TransVecs[1]=0;
          H->TransVecs[2]=0;
        }
        else
        {  // a and b infinite c finite
          return luaL_error(L,"Sorry input for slabs not completed (setting of translation vector needs to be done");
        }
      }
      else
      {
        if(NormSqr(Cell->rUnitCell[1])==INFINITY)
        {  // a and c infinte b finite
          return luaL_error(L,"Sorry input for slabs not completed (setting of translation vector needs to be done");
        }
        else
        {  // a infinite b and c finite
          return luaL_error(L,"Sorry input for slabs not completed (setting of translation vector needs to be done");
        }
      }
    }
    else
    {
      if(NormSqr(Cell->rUnitCell[1])==INFINITY)
      {
        if(NormSqr(Cell->rUnitCell[2])==INFINITY)
        {  // a finite b and c infinite
          return luaL_error(L,"Sorry input for slabs not completed (setting of translation vector needs to be done");
        }
        else
        {  // a and c finite b infinite
          return luaL_error(L,"Sorry input for slabs not completed (setting of translation vector needs to be done");
        }
      }
      else
      {
        if(NormSqr(Cell->rUnitCell[2])==INFINITY)
        {  // a and b finite c infinite
          return luaL_error(L,"Sorry input for slabs not completed (setting of translation vector needs to be done");
        }
        else
        {  // a, b and c finite
          // {x,y,z} = i a + j b + k c + AtomPos2 - AtomPos1
          // r = {x,y,z} - AtomPos2 + AtomPos1
          // i = Cross[b,c] . r / (Cross[a,b].c)
          // j = Cross[c,a] . r / (Cross[a,b].c)
          // k = Cross[a,b] . r / (Cross[a,b].c)
          CellVolume=ParallelalPipedVolume(Cell->rUnitCell[0], Cell->rUnitCell[1], Cell->rUnitCell[2]);
          if(fabs(CellVolume)<100*sqrt(DBL_MIN)) return luaL_error(L, "ReadHopping called with a cell defined that has a zero volume:\n    a={%d %d %d}\n    b={%d %d %d}\n    c={%d %d %2d}\n",Cell->rUnitCell[0][0],Cell->rUnitCell[0][1],Cell->rUnitCell[0][2],Cell->rUnitCell[1][0],Cell->rUnitCell[1][1],Cell->rUnitCell[1][2],Cell->rUnitCell[2][0],Cell->rUnitCell[2][1],Cell->rUnitCell[2][2]);
          r[0]=H->ConnectingVector[0]+Cell->Atom[H->Atompos1].Position[0]-Cell->Atom[H->Atompos2].Position[0];
          r[1]=H->ConnectingVector[1]+Cell->Atom[H->Atompos1].Position[1]-Cell->Atom[H->Atompos2].Position[1];
          r[2]=H->ConnectingVector[2]+Cell->Atom[H->Atompos1].Position[2]-Cell->Atom[H->Atompos2].Position[2];
          idouble[0]=ParallelalPipedVolume(Cell->rUnitCell[1], Cell->rUnitCell[2], r)/CellVolume;
          idouble[1]=ParallelalPipedVolume(Cell->rUnitCell[2], Cell->rUnitCell[0], r)/CellVolume;
          idouble[2]=ParallelalPipedVolume(Cell->rUnitCell[0], Cell->rUnitCell[1], r)/CellVolume;
          H->TransVecs[0]=(int)(idouble[0]+(idouble[0]>0?1:-1)*0.5);
          H->TransVecs[1]=(int)(idouble[1]+(idouble[1]>0?1:-1)*0.5);
          H->TransVecs[2]=(int)(idouble[2]+(idouble[2]>0?1:-1)*0.5);
          if((fabs(idouble[0]-H->TransVecs[0])>0.0000001)||(fabs(idouble[1]-H->TransVecs[1])>0.0000001)||(fabs(idouble[2]-H->TransVecs[2])>0.0000001)) return luaL_error(L, "ReadHopping connecting vector given in absolute coordinates, but not an integer number of unit cells plus the distance in the unit cell between the atoms:\n   {x,y,z} = {%f %f %f} (%f %f %f)\n   {i,j,k} = {%f %f %f}\n",(double)H->ConnectingVector[0],(double)H->ConnectingVector[1],(double)H->ConnectingVector[2],(double)H->TransVecs[0],(double)H->TransVecs[1],(double)H->TransVecs[2],idouble[0],idouble[1],idouble[2]);
        }
      }
    }
  }
  else // relative units
  {
    for(i=0;i<3;i++)
    {
      lua_rawgeti(L, -1, i+1);
      buff=(double)luaL_checknumber(L, -1);
      H->TransVecs[i]=(int)buff;
      if(fabs(buff-H->TransVecs[i])>100*DBL_EPSILON) return luaL_error(L, "Units for atom positions is set to Relative, i.e. all vectors must be given in units of a,b and c.\nReadHopping expects a table of the form {ai,bi,ci} for the vector defining the connection between atom %s and atom %s\nWith ai, bi and ci being integers, found %d\n",AtomShell1,AtomShell2,buff);
      lua_pop(L,1);
    }
    lua_pop(L,1);
    // create Connecting Vector from Trans Vec
    for(i=0;i<3;i++) // i loops over x,y,z
    {
      H->ConnectingVector[i]=0;
      for(j=0;j<3;j++) // j loops over a,b,c
      {
        H->ConnectingVector[i] += (H->TransVecs[j]+Cell->Atom[H->Atompos2].Position[j]-Cell->Atom[H->Atompos2].Position[j]) * Cell->rUnitCell[j][i];
      }
    }
  }
  // read H
  lua_rawgeti(L, narg, 4);
  if((lua_type(L,-1)!=LUA_TTABLE) && (luaL_testudata(L, -1, "Matrix_Type") == NULL) ) return luaL_error(L, "ReadHopping expects a matrix for the definition of H between atom %s and %s\n",AtomShell1,AtomShell2);
  H->Hop = *luaL_checkMatrix(L, -1);
  if( (H->Hop.Ni!=Cell->Atom[H->Atompos1].shelldegeneracy[H->Shellpos1]) || (H->Hop.Nj!=Cell->Atom[H->Atompos2].shelldegeneracy[H->Shellpos2]) )
  {
    printf("Matrix found\n");
    PrintMatrix(H->Hop);
    printf("Expected dimension, %d %d\n",Cell->Atom[H->Atompos1].shelldegeneracy[H->Shellpos1],Cell->Atom[H->Atompos2].shelldegeneracy[H->Shellpos2]);
    return luaL_error(L, "ReadHopping found inconsistency in matrix defining the hopping");
  }
  lua_pop(L,1);
  return 0;
}

/*
 * push a Hopping definition to the Lua stack
 * format
 *          Hopping = {"AtomName1.ShellName1","AtomName2.ShellName2",Vec,H}
 *        with H a table of tables of real or complex numbers / matrix
 *
 */
static int LuaPushHopping(lua_State *L, CrystalStructureType const *Cell, TightBindingDefHopping const *H)
{
  //TODO:
  return 1;
}


static int LuaTightBindingIndex(lua_State * L)
{
  uint32_t i,j;
  TightBindingDef *a = (TightBindingDef *) luaL_checkudata(L, 1, "TightBinding_Type");
  const char *str = lua_tostring(L, 2);
  if (strcmp(str, "Name") == 0)
  {
    lua_pushstring(L, a->Name);
    return 1;
  }
  else if (strcmp(str, "Cell")==0)
  {
    lua_createtable(L,3,0);
    for(i=0;i<3;i++)
    {
      lua_createtable(L,3,0);
      for(j=0;j<3;j++)
      {
        lua_pushnumber(L, a->Cell.rUnitCell[i][j]);
        lua_rawseti(L, -2, j+1);
      }
      lua_rawseti(L, -2, i+1);
    }
    return 1;
  }
  else if (strcmp(str, "ReciprocalCell")==0)
  {
    lua_createtable(L,3,0);
    for(i=0;i<3;i++)
    {
      lua_createtable(L,3,0);
      for(j=0;j<3;j++)
      {
        lua_pushnumber(L, a->Cell.gUnitCell[i][j]);
        lua_rawseti(L, -2, j+1);
      }
      lua_rawseti(L, -2, i+1);
    }
    return 1;
  }
  else if (strcmp(str, "NAtoms") == 0)
  {
    lua_pushunsigned(L, a->Cell.NAtoms);
    return 1;
  }
  else if (strcmp(str, "Atoms") == 0)
  {
    lua_createtable(L, a->Cell.NAtoms, 0);
    for(i=0;i<a->Cell.NAtoms;i++)
    {
      LuaPushAtom(L, &(a->Cell.Atom[i]));
      lua_rawseti(L, -2, i+1);
    }
    return 1;
  }
  else if (strcmp(str, "Units")==0)
  {
    lua_createtable(L,3,0);
    if((a->Cell.Units&1)==1)
    {
      lua_pushfstring(L, "2Pi");
      lua_rawseti(L, -2, 1);
    }
    if((a->Cell.Units&1)==0)
    {
      lua_pushfstring(L, "NoPi");
      lua_rawseti(L, -2, 1);
    }
    if((a->Cell.Units&6)==0)
    {
      lua_pushfstring(L, "Angstrom");
      lua_rawseti(L, -2, 2);
    }
    if((a->Cell.Units&6)==2)
    {
      lua_pushfstring(L, "Bohr");
      lua_rawseti(L, -2, 2);
    }
    if((a->Cell.Units&6)==4)
    {
      lua_pushfstring(L, "nm");
      lua_rawseti(L, -2, 2);
    }
    if((a->Cell.Units&8)==0)
    {
      lua_pushfstring(L, "Absolute");
      lua_rawseti(L, -2, 3);
    }
    if((a->Cell.Units&8)==8)
    {
      lua_pushfstring(L, "Relative");
      lua_rawseti(L, -2, 3);
    }
    return 1;
  }
  else if (strcmp(str, "Hopping") == 0)
  {
    //TODO: push a table to the stack with the hopping definitions as defined when reading the variable
    return 1;
  }
  else if (strcmp(str, "Hk")==0)
  {
    //TODO: push a matrix to the stack as table representing Hk(kx,ky,kz) read options to define either constant or kx phase in the unit cell
    return 1;
  }
  else if (strcmp(str, "NF")==0)
  {
    lua_pushinteger(L, a->NFermion);
    return 1;
  }
  else if (strcmp(str, "ShiftOnsiteEnergy")==0)
  {
    lua_pushvalue(L, 1);
    lua_pushcclosure(L, LuaTightBindingShiftOnsiteEnergyMethod, 1);
    return 1;
  }
  return luaL_error(L, "Error index %s not found in TightBinding Object\n",str);
}


static int LuaTightBindingNewIndex(lua_State * L)
{
  uint32_t i,j,len;
  TightBindingDef *a = (TightBindingDef *) luaL_checkudata(L, 1, "TightBinding_Type");
  const char *str = lua_tostring(L, 2);

  if(strcmp(str, "Name")==0)
  {
    strncpy(a->Name, luaL_checkstring(L, 3), MAXStrLength-1);
  }
  else if(strcmp(str, "Cell")==0)
  {
    luaL_checktype(L, 3, LUA_TTABLE);
    len = (uint32_t)lua_objlen(L, 3);
    if(len!=3)
    {
      return luaL_error(L, "Error TightBindingObject.Cell should be of the form { a, b, c }\n");
    }
    for(i=0;i<3;i++)
    {
      lua_rawgeti(L, 3, i+1);
      luaL_checktype(L, -1, LUA_TTABLE);
      len = (uint32_t)lua_objlen(L, -1);
      if(len!=3)
      {
        return luaL_error(L, "Error TightBindingObject.Cell should be of the form { {a1,a2,a3}, {b1,b2,b3}, {c1,c2,c3} }\n");
      }
      for(j=0;j<3;j++)
      {
        lua_rawgeti(L, -1, j+1);
        a->Cell.rUnitCell[i][j]=(double)luaL_checknumber(L, -1);
        lua_pop(L, 1);
      }
      lua_pop(L, 1);
    }
    if(fabs(ParallelalPipedVolume(a->Cell.rUnitCell[0], a->Cell.rUnitCell[1], a->Cell.rUnitCell[2]))<100*sqrt(DBL_MIN)) return luaL_error(L, "Error TightBindingObject.Cell should be of the form { {a1,a2,a3}, {b1,b2,b3}, {c1,c2,c3} } and span a non-zero volume.\n       i.e. a, b and c must be independent vectors and they are not\n");
    ReciprocalLatticeVector(a->Cell.rUnitCell, a->Cell.gUnitCell,a->Cell.Units);
    //TODO: set ConnectionVector anew
  }
  else if(strcmp(str, "Atoms")==0)
  {
    if(a->nH!=0)
    {
      luaL_warning(L, "Warning Changed the settings of atoms for a tight binding hamiltonian (%s)\n Previous definition of all hoppings will be arrased\n",a->Name);
      for(i=0;i<a->nH;i++)
      {
        FreeMatrix(&(a->H[i].Hop));
      }
      free(a->H);
      a->nH=0;
      if(a->FermionIndex!=NULL)
      {
        for (i = 0; i < a->Cell.NAtoms; i++)
        {
          free(a->FermionIndex[i]);
        }
        if(a->Cell.NAtoms!=0) free(a->FermionIndex);
      }
      a->Cell.NAtoms=0;
    }
    // format of atoms  = {Atom1, Atom2, ...}
    //           Atomi  = {Name, Position, {Shell1, Shell2, ...}}
    //           Shelli = {Name, {orb1, orb2, ...}}
    if(lua_type(L,3)==LUA_TTABLE)
    {
      len = (uint32_t)lua_objlen(L, 3);
      a->Cell.NAtoms=len;
      a->Cell.Atom= (AtomType*)malloc(len*sizeof(AtomType));
      for(i=0;i<len;i++)
      {
        lua_rawgeti(L, 3, i+1);
        LuaReadAtom(L, -1, &(a->Cell.Atom[i]));
        lua_pop(L, 1);
      }
    }
    else
    {
      return luaL_error(L, "Error TightBindingObject.Atoms expect a list of atoms\n");
    }
    SetNFermionIndexTightBindingDef(a);
  }
  else if(strcmp(str, "Units")==0)
  {
    if(lua_type(L,3)==LUA_TNUMBER)
    {
      a->Cell.Units=(uint32_t)luaL_checkunsigned(L, 3);
    }
    else if((lua_type(L,3)==LUA_TSTRING)||(lua_type(L,3)==LUA_TTABLE))
    {
      if(lua_type(L,3)==LUA_TSTRING)
      {
        str = luaL_checkstring(L, 3);
        len=0;
      }
      else
      {
        len = (uint32_t)lua_objlen(L, 3);
      }
      for(i=0;i<=len;i++)
      {
        if(i==0) i++;
        if(len!=0)
        {
          lua_rawgeti(L, 3, i);
          str = luaL_checkstring(L, -1);
        }
        if(strcmp(str, "2Pi")==0)
        {
          a->Cell.Units |= (uint32_t)1;       // set bit 0 to 1
        }
        else if(strcmp(str, "NoPi")==0)
        {
          a->Cell.Units &= ~((uint32_t)1);    // set bit 0 to 0
        }
        else if((strcmp(str, "Angstrom")==0)||(strcmp(str, "A")==0))
        {
          a->Cell.Units &= ~((uint32_t)1<<1); // set bit 1 to 0
          a->Cell.Units &= ~((uint32_t)1<<2); // set bit 2 to 0
        }
        else if((strcmp(str, "Bohr")==0)||(strcmp(str, "a0")==0))
        {
          a->Cell.Units |= (uint32_t)1<<1;    // set bit 1 to 1
          a->Cell.Units &= ~((uint32_t)1<<2); // set bit 2 to 0
        }
        else if((strcmp(str, "nanometer")==0)||(strcmp(str, "nm")==0))
        {
          a->Cell.Units &= ~((uint32_t)1<<1); // set bit 1 to 0
          a->Cell.Units |= (uint32_t)1<<2;    // set bit 2 to 1
        }
        else if(strcmp(str, "Absolute")==0)
        {
          a->Cell.Units &= ~((uint32_t)1<<3); // set bit 3 to 0
        }
        else if(strcmp(str, "Relative")==0)
        {
          a->Cell.Units |= (uint32_t)1<<3;    // set bit 3 to 1
        }
        else
        {
          return luaL_error(L, "Error TightBindingObject.Unit %s is not defined\n",str);
        }
        if(len!=0)
        {
          lua_pop(L, 1);
        }
      }
    }
    else
    {
      return luaL_error(L, "Error not a valid type assignment for TightBinding Object (%s) should be integer, string or table of strings\n",lua_typename(L, lua_type(L, 3)));
    }
    ReciprocalLatticeVector(a->Cell.rUnitCell, a->Cell.gUnitCell,a->Cell.Units);
    //TODO: set ConnectionVector anew
  }
  else if(strcmp(str, "Hopping")==0)
  {
    //
    if((lua_type(L,3)!=LUA_TTABLE)) return luaL_error(L, "Error TightBinding.Hopping should be a list of hopping definitions {{def1}, {def2}, ...}\n");
    a->nH = (uint32_t)lua_objlen(L, 3);
    a->H= (TightBindingDefHopping*)malloc(a->nH*sizeof(TightBindingDefHopping));
    for(i=0;i<a->nH;i++) // make sure that whne LuaReadHopping fails we can free a without addressing unassigned / random pointers
    {
      a->H[i].Hop.Ni=0;
      a->H[i].Hop.Nj=0;
      a->H[i].Hop.QComplex=0;
      a->H[i].Hop.VRe=NULL;
    }
    for(i=0;i<a->nH;i++)
    {
      lua_rawgeti(L, 3, i+1);
      LuaReadHopping(L, -1, &(a->Cell), &(a->H[i]));
      lua_pop(L, 1);
    }
  }
  else if(strcmp(str, "Hk")==0)
  {
    return luaL_error(L, "Error can not set Hk in a TightBinding Object\n",str);
  }
  else
  {
    return luaL_error(L, "Error index %s not found in TightBinding Object\n",str);
  }
  return 0;
}

static int LuaTightBindingToString(lua_State * L)
{
  TightBindingDef *res = (TightBindingDef *) luaL_checkudata(L, 1, "TightBinding_Type");
  char *Buffer;
  StringPrintTightBindingDef(*res, &Buffer);
  lua_pushstring(L, Buffer);
  free(Buffer);
  return 1;
}

static int LuaCrystalStructureToCIF(lua_State *L)
{
  char *Buffer;
  TightBindingDef *res = (TightBindingDef *) luaL_checkudata(L, 1, "TightBinding_Type");
  CrystalStructureTypeToCif(&(res->Cell), &Buffer);
  lua_pushstring(L, Buffer);
  free(Buffer);
  return 1;
}

static int LuaTightBindingCall(lua_State * L)
{
  uint32_t len = (uint32_t)lua_gettop(L);
  double k[3] = {0.0,0.0,0.0};
  bool AsOperator=false;
  bool ConstantPhaseInUnitCell=false;
  bool AddSpin=false;
  LuaOptionsType LuaOptions[]={
    {"AsOperator","bool"             ,&AsOperator},
    {"ConstantPhaseInUnitCell","bool",&ConstantPhaseInUnitCell},
    {"AddSpin","bool",&AddSpin},
    {NULL, NULL, NULL}};
  if(len!=0)
  { // check for options if options found reduce len by one
    if(LuaTestOptionsList(L, -1))
    {
      LuaCheckOptions(L, -1, LuaOptions);
      len--;
    }
  }
  if((len<2)||(len>4)) return luaL_error(L, "Error TightBindingObject(...) should be called with a k-vector of length 1,2, or 3 and possible options\n");
  TightBindingDef *a = (TightBindingDef *) luaL_checkudata(L, 1, "TightBinding_Type");
  if((len==2)&&(lua_type(L,2)==LUA_TTABLE))
  {
    len = (int)lua_objlen(L, 2);
    if((len<1)||(len>3)) luaL_error(L, "Error TightBindingObject(...) should be called with a k-vector of length 1,2, or 3 and possible options\n");
    for(int i=0;i<len;i++)
    {
      lua_rawgeti(L, 2, i+1);
      k[i]=luaL_checknumber(L, -1);
      lua_pop(L, 1);
    }
  }
  else
  {
    k[0]=luaL_checknumber(L, 2);
    if(len>2) k[1]=luaL_checknumber(L, 3);
    if(len>3) k[2]=luaL_checknumber(L, 4);
  }
  if((AsOperator==false)&&(ConstantPhaseInUnitCell==false)&&(AddSpin==false))
  {
    MatrixType H;
    CreateHamiltonianFromTightBindingDef(*a, k, &H);
    lua_push(L,H);
    FreeMatrix(&H);
    return 1;
  }
  if((AsOperator==false )&&(ConstantPhaseInUnitCell==false)&&(AddSpin==true ))
  {
    MatrixType H,Hp;
    CreateHamiltonianFromTightBindingDef(*a, k, &Hp);
    MatrixCopyAddSpin(&H, Hp);
    lua_push(L,H);
    FreeMatrix(&H);
    FreeMatrix(&Hp);
    return 1;
  }
  if((AsOperator==false )&&(ConstantPhaseInUnitCell==true )&&(AddSpin==false))
  {
    MatrixType H;
    CreateHamiltonianFromTightBindingDefConstantPhaseInUnitCell(*a, k, &H);
    lua_push(L,H);
    FreeMatrix(&H);
    return 1;
  }
  if((AsOperator==false )&&(ConstantPhaseInUnitCell==true )&&(AddSpin==true ))
  {
    MatrixType H,Hp;
    CreateHamiltonianFromTightBindingDefConstantPhaseInUnitCell(*a, k, &Hp);
    MatrixCopyAddSpin(&H, Hp);
    lua_push(L,H);
    FreeMatrix(&H);
    FreeMatrix(&Hp);
    return 1;

  }
  if((AsOperator==true )&&(ConstantPhaseInUnitCell==false)&&(AddSpin==false))
  {
    
  }
  if((AsOperator==true )&&(ConstantPhaseInUnitCell==false)&&(AddSpin==true ))
  {

  }
  if((AsOperator==true )&&(ConstantPhaseInUnitCell==true )&&(AddSpin==false))
  {
    
  }
  if((AsOperator==true )&&(ConstantPhaseInUnitCell==true )&&(AddSpin==true ))
  {

  }
  return luaL_error(L, "Error TightBindingObject(...) incompattible set of options");
}

static int LuaTightBindingDestroy(lua_State * L)
{
  TightBindingDef *a = (TightBindingDef *) luaL_checkudata(L, 1, "TightBinding_Type");
  FreeTightBindingDef(a);
  return 0;
}

static int LuaDresdenFPLOToString(lua_State * L)
{
  lua_pushstring(L, "FPLO input");
  return 1;
}

static int LuaDresdenFPLODestroy(lua_State * L)
{
	DresdenFPLO *a = (DresdenFPLO *) luaL_checkudata(L, 1, "DresdenFPLO_Type");
  FreeDresdenFPLO(a);
  return 0;
}

static const struct luaL_Reg TightBinding_func[] = {
  {"New",               LuaNewTightBinding},
  {"Chop",              LuaTightBindingChop},
  {"Copy",              LuaTightBindingCopy},
  {"AddSpin",           LuaTightBindingAddSpin},
  {"Derivative",        LuaTightBindingDerivative},
  {"ShiftOnsiteEnergy", LuaTightBindingShiftOnsiteEnergy},
  {"CalculateG",        LuaCalculateG},
  {"CalculateMu",       LuaCalculateMu},
  {"CalculateBands",    LuaCalculateBands},
  {"SlaterKosterTable", LuaSlaterKosterTable},
  {"ExportCIF",         LuaCrystalStructureToCIF},
  {NULL, NULL}
};
  

static const struct luaL_Reg TightBinding_meta[] = {
  {"__index",    LuaTightBindingIndex},
  {"__newindex", LuaTightBindingNewIndex},
  {"__tostring", LuaTightBindingToString},
  {"__call",     LuaTightBindingCall},
  {"__gc",       LuaTightBindingDestroy},
  {NULL, NULL}
};

static const struct luaL_Reg StuttgartCTRL_meta[] = {
  {NULL, NULL}
};

static const struct luaL_Reg StuttgartHAMR_meta[] = {
  {NULL, NULL}
};

static const struct luaL_Reg DresdenFPLO_meta[] = {
	{"__gc", LuaDresdenFPLODestroy},
  {"__tostring", LuaDresdenFPLOToString},
  {NULL, NULL}
};
