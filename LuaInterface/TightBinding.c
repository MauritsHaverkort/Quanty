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
    luaL_setmetatable(L, "TightBinding_Type");
    TB->nH=0;
    TB->NFermion=0;
    TB->Cell.NAtoms=0;
    TB->Cell.Units=1;
    TB->FermionIndex=NULL;
  }
  else
  {
    return luaL_error(L, "Error in NewTightBinding(arg) unrecognized arguments arg, should be empty\n");
  }
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


/* on stack: 
 *   TB : TightBinding_Type
 *   {Emin,Emax,NE}
 *   Options
 *       "Nk" -> {Nkx,Nky,Nkz} number of k-points       standard {40,40,40}
 *       Method -> "mesh", ...                          standard "mesh"
 *       Character -> {matrix} or "full" or "total" standard "full" rho is psi^* (x) psi, rho2 = M^T . psi^* (x) psi . M^*
 *       Broadening -> "None", {"Gauss",FMHM}, {"Lorentz, FWHM}
 *   returns a Spectra_Type
 */
static int LuaCalculateG(lua_State *L)
{
  //TODO: all
  struct
  {
    double Nk[3];
    char Method[64];
    char Character[64];
    char Broadening[64];
    double FWHM;
    double Emin;
    double Emax;
    uint32_t NE;
  } Options;
  Options.Nk[0]=40;
  Options.Nk[1]=40;
  Options.Nk[2]=40;
  strcpy(Options.Method,     "Mesh");
  strcpy(Options.Character,  "Full");
  strcpy(Options.Broadening, "None");
  Options.FWHM=0;
  Options.Emin=-10;
  Options.Emax=10;
  Options.NE=2000;
  LuaOptionsType LuaCalculateGOptions[]={
    {"Nk","array-uint32_t-3",&(Options.Nk[0])},
    {"Method","char*",&(Options.Method[0])},
    {"Character","char*",&(Options.Character[0])},
    {"Broadening", "char*", &(Options.Broadening[0])},
    {"FWHM","double",&(Options.FWHM)},
    {"Emin","double",&(Options.Emin)},
    {"Emax","double",&(Options.Emax)},
    {"NE","uint32_t",&(Options.NE)},
    {NULL, NULL, NULL}};
  
  return 1;
}

/*
 *
 */
static int LuaCalculateBands(lua_State *L)
{
  //TODO: all
  
  return 1;
}

/*
 *
 */
static int LuaCalculateFermi(lua_State *L)
{
  //TODO: all
  
  return 1;
}

/*
 *
 */
static int LuaCalculateDOS(lua_State *L)
{
  return LuaCalculateG(L);
}

/*  Creates an Hamiltonian representing a real space cluster
 *  Input on stack
 *    A tight binding Definition
 *    A cluster definition
 *       Possible formats:
 *         {"periodic", {a,b,c}, atoms_to_include_optional}
 *         {"open", atoms_to_include}
 *    Aditional Operator to add (optional)
 *         Fermion index should follow the definition in the original tight binding Hamiltonian, when spin is added it is added in the form down/up, down/up, ...
 *    Options (optional)
 */
static int LuaCreateClusterHamiltonian(lua_State *L)
{
  uint32_t len, lenclusterdef,i,j,nreturn=1;
  char *tokenread;
  TightBindingDef *a, aSuper, *aSuperReturn;
  CrystalStructureSuperCellType SuperCell;
  uint32_t periodic;
  OperatorType *H;
  CompactVectorType k;
  Bool AddSpin=False;
  Bool QOperators=False, QOptions=False, ReturnTBSuperCell=False, ReturnTBCellPositions=False, ReturnTBKVectors=False;
  OperatorType *AditionalOperator=NULL;
  MatrixType U;
  k.Ni=3;
  k.QComplex=0;
  k.VRe=calloc(3,sizeof(double));
  LuaOptionsType LuaCreateClusterHamiltonianOptions[]={
    {"AddSpin", "bool", &AddSpin},
    {"ReturnTBSuperCell", "bool", &ReturnTBSuperCell},
    {"ReturnTBCellPositions", "bool", &ReturnTBCellPositions},
    {"ReturnTBKVectors", "bool", &ReturnTBKVectors},
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
  /* read the cluster defintion from the second position on the stack */
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
  }
  else
  {
    ReturnTBSuperCell=False;
    ReturnTBKVectors=False;
    ReturnTBCellPositions=False;
    //TODO:
    return luaL_error(L, "Error in CreateClusterHamiltonian open boundary conditions not yet implemented \n");
  }
  /* deterimine if there are aditional operators or aditional options */
  if(len==2)
  {
    QOperators=False;
    QOptions=False;
  }
  else
  {
    if(len==4)
    {
      QOperators=True;
      QOptions=True;
    }
    /* read additional operators from the stack */
    if (luaL_testudata(L, 3, "Operator_Type") != NULL)
    {
      QOperators=True;
      AditionalOperator = (OperatorType *) luaL_checkudata(L, 3, "Operator_Type");
      if(len==3)
      {
        QOptions=False;
      }
    }
    else
    {
      QOperators=False;
      QOptions=True;
      if(len==4)
      {
        return luaL_error(L, "Error CreateClusterHamiltonian third argument should be an operator defining aditional interactions within the unit cell\n");
      }
    }
  }
  /* read options */
  if(QOptions==True)
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
  if(AddSpin==True)
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
  else
  {
    FreeTightBindingDef(&aSuper);
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
  return nreturn;
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
    a->ShellName = malloc(a->Nshell*sizeof(char*));
    a->shelldegeneracy = malloc(a->Nshell*sizeof(int));
    a->OrbitalName = malloc(a->Nshell*sizeof(char**));
    for(i=0;i<a->Nshell;i++)
    {
      lua_rawgeti(L, -1, i+1);
      if(lua_type(L,-1)!=LUA_TTABLE) return luaL_error(L, "Error ReadAtom shell %d of atom %s should be a list containing the name and a list of orbitals\n",i+1,a->Name);
      len = (uint32_t)lua_objlen(L, -1);
      if(len!=2) return luaL_error(L, "Error ReadAtom shell %d of atom %s should be a list containing the name and a list of orbitals\n",i+1,a->Name);
      lua_rawgeti(L, -1, 1);
      a->ShellName[i]=malloc(MAXStrLength*sizeof(char));
      strncpy(a->ShellName[i], luaL_checkstring(L, -1), MAXStrLength-1);
      lua_pop(L, 1);
      lua_rawgeti(L, -1, 2);
      if(lua_type(L,-1)!=LUA_TTABLE) return luaL_error(L, "Error ReadAtom shell %s(%d) of atom %s should be a list containing the name and a list of orbitals\n",a->ShellName[i],i+1,a->Name);
      a->shelldegeneracy[i] = (uint32_t)lua_objlen(L, -1);
      a->OrbitalName[i] = malloc(a->shelldegeneracy[i]*sizeof(char*));
      for(j=0;j<a->shelldegeneracy[i];j++)
      {
        lua_rawgeti(L, -1, j+1);
        a->OrbitalName[i][j]=malloc(MAXStrLength*sizeof(char));
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
            lua_pushstring(L, &(a->ShellName[i][j]));
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
 *          Hopping = {"AtomName1.ShellName1","AtomName2.ShellName2",Vec,H}
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
  // test input
  if((lua_type(L,narg)!=LUA_TTABLE)) return luaL_error(L, "ReadHopping expects a table of the form {\"AtomName1.ShellName1\",\"AtomName2.ShellName2\",Vec,H}\n");
  if (lua_objlen(L, narg)!=4) return luaL_error(L, "ReadHopping expects a table with 4 elements of the form {\"AtomName1.ShellName1\",\"AtomName2.ShellName2\",Vec,H}\n");
  // read atom and shell 1
  lua_rawgeti(L, narg, 1);
  strncpy(AtomShell1, luaL_checkstring(L, -1), 2*MAXStrLength);
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
      H->Atompos1=i;
      break;
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
  lua_rawgeti(L, narg, 2);
  strncpy(AtomShell2, luaL_checkstring(L, -1), 2*MAXStrLength);
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
      H->Atompos2=i;
      break;
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
  if((lua_type(L,-1)!=LUA_TTABLE)) return luaL_error(L, "ReadHopping expects a matrix for the definition of H between atom %s and %s\n",AtomShell1,AtomShell2);
  H->Hop = *luaL_checkMatrix(L, -1);
  
  
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
    ResiprocalLatticeVector(a->Cell.rUnitCell, a->Cell.gUnitCell,a->Cell.Units);
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
      a->Cell.Atom=malloc(len*sizeof(AtomType));
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
    ResiprocalLatticeVector(a->Cell.rUnitCell, a->Cell.gUnitCell,a->Cell.Units);
    //TODO: set ConnectionVector anew
  }
  else if(strcmp(str, "Hopping")==0)
  {
    //
    if((lua_type(L,3)!=LUA_TTABLE)) return luaL_error(L, "Error TightBinding.Hopping should be a list of hopping definitions {{def1}, {def2}, ...}\n");
    a->nH = (uint32_t)lua_objlen(L, 3);
    a->H=malloc(a->nH*sizeof(TightBindingDefHopping));
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
  return 1;
}

static int LuaTightBindingDestroy(lua_State * L)
{
  TightBindingDef *a = (TightBindingDef *) luaL_checkudata(L, 1, "TightBinding_Type");
  FreeTightBindingDef(a);
  return 0;
}




static const struct luaL_Reg TightBinding_meta[] = {
  {"__index", LuaTightBindingIndex},
  {"__newindex", LuaTightBindingNewIndex},
  {"__tostring", LuaTightBindingToString},
  {"__gc", LuaTightBindingDestroy},
  {NULL, NULL}
};

static const struct luaL_Reg StuttgartCTRL_meta[] = {
  {NULL, NULL}
};

static const struct luaL_Reg StuttgartHAMR_meta[] = {
  {NULL, NULL}
};
