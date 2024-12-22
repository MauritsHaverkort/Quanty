/*
 *
 */
static int LuaMatrixTranspose(lua_State * L)
{
  CompactMatrixType *A;
  uint32_t len = (uint32_t)lua_gettop(L);
  A = (CompactMatrixType *) luaL_checkCompactMatrix(L, 1);
  if(len!=1) luaL_error(L, "Matrix.Transpose called with %d arguments, only one expected",len);
  if (luaL_testudata(L, 1, "Matrix_Type") == NULL)
  {
    CompactMatrixTranspose(A);
    LuaPushTable_CompactMatrixType(L, A);
    FreeCompactMatrix(A);
    free(A);
  }
  else
  {
    CompactMatrixType *B;
    B = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
    luaL_setmetatable(L, "Matrix_Type");
    B->OnLuaStack=true;
    CompactMatrixCopy(*A, B);
    CompactMatrixTranspose(B);
  }
  return 1;
}

/*
 *
 */
static int LuaMatrixConjugateTranspose(lua_State * L)
{
  CompactMatrixType *A;
  uint32_t len = (uint32_t)lua_gettop(L);
  A = (CompactMatrixType *) luaL_checkCompactMatrix(L, 1);
  if(len!=1) luaL_error(L, "Matrix.ConjugateTranspose called with %d arguments, only one expected",len);
  if (luaL_testudata(L, 1, "Matrix_Type") == NULL)
  {
    CompactMatrixConjugateTranspose(A);
    LuaPushTable_CompactMatrixType(L, A);
    FreeCompactMatrix(A);
    free(A);
  }
  else
  {
    CompactMatrixType *B;
    B = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
    luaL_setmetatable(L, "Matrix_Type");
    B->OnLuaStack=true;
    CompactMatrixCopy(*A, B);
    CompactMatrixConjugateTranspose(B);
  }
  return 1;
}

/*
 *
 */
static int LuaMatrixConjugate(lua_State * L)
{
  CompactMatrixType *A;
  uint32_t len = (uint32_t)lua_gettop(L);
  A = (CompactMatrixType *) luaL_checkCompactMatrix(L, 1);
  if(len!=1) luaL_error(L, "Matrix.Conjugate called with %d arguments, only one expected",len);
  if (luaL_testudata(L, 1, "Matrix_Type") == NULL)
  {
    CompactMatrixConjugate(A);
    LuaPushTable_CompactMatrixType(L, A);
    FreeCompactMatrix(A);
    free(A);
  }
  else
  {
    CompactMatrixType *B;
    B = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
    luaL_setmetatable(L, "Matrix_Type");
    B->OnLuaStack=true;
    CompactMatrixCopy(*A, B);
    CompactMatrixConjugate(B);
  }
  return 1;
}

/*
 *
 */
static int LuaMatrixInverse(lua_State * L)
{
  CompactMatrixType *A;
  uint32_t len = (uint32_t)lua_gettop(L);
  A = (CompactMatrixType *) luaL_checkCompactMatrix(L, 1);
  if(len!=1) luaL_error(L, "Matrix.Inverse called with %d arguments, only one expected",len);
  if(A->Ni!=A->Nj) return luaL_error(L,"Matrix Inverse only works on square matrices\n");
  if (luaL_testudata(L, 1, "Matrix_Type") == NULL)
  {
    CompactMatrixInverse(A);
    LuaPushTable_CompactMatrixType(L, A);
    FreeCompactMatrix(A);
    free(A);
  }
  else
  {
    CompactMatrixType *B;
    B = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
    luaL_setmetatable(L, "Matrix_Type");
    B->OnLuaStack=true;
    CompactMatrixCopy(*A, B);
    CompactMatrixInverse(B);
  }
  return 1;
}

/*
 *
 */
static int LuaMatrixSqrt(lua_State * L)
{
  CompactMatrixType *A,*B;
  uint32_t len = (uint32_t)lua_gettop(L);
  A = (CompactMatrixType *) luaL_checkCompactMatrix(L, 1);
  if(A->Ni!=A->Nj) return luaL_error(L,"Matrix sqrt only works on square matrices\n");
  if(len!=1) luaL_error(L, "Matrix.Sqrt called with %d arguments, only one expected",len);
  if(QCompactMatrixHermitian(*A, DBL_EPSILON)==0)
  {
    return luaL_error(L,"Matrix sqrt only works on Hermitian matrices\n");
    if (luaL_testudata(L, 1, "Matrix_Type") == NULL)
    { // matrix on stack as table
      B= (CompactMatrixType*)malloc(sizeof(CompactMatrixType));
      if(NonHermitianCompactMatrixSqrt(A,B)) return luaL_error(L,"Matrix Sqrt failed\n");
      LuaPushTable_CompactMatrixType(L, B);
      FreeCompactMatrix(A);
      free(A);
      FreeCompactMatrix(B);
      free(B);
    }
    else
    { // matrix on stack as userdata
      B = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
      luaL_setmetatable(L, "Matrix_Type");
      B->OnLuaStack=true;
      if(NonHermitianCompactMatrixSqrt(A,B)) return luaL_error(L,"Matrix Sqrt failed\n");
    }
  }
  else
  {
    if (luaL_testudata(L, 1, "Matrix_Type") == NULL)
    { // matrix on stack as table
      B= (CompactMatrixType*)malloc(sizeof(CompactMatrixType));
      if(CompactMatrixSqrt(A,B)) return luaL_error(L,"Matrix Sqrt failed\n");
      LuaPushTable_CompactMatrixType(L, B);
      FreeCompactMatrix(A);
      free(A);
      FreeCompactMatrix(B);
      free(B);
    }
    else
    { // matrix on stack as userdata
      B = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
      luaL_setmetatable(L, "Matrix_Type");
      B->OnLuaStack=true;
      if(CompactMatrixSqrt(A,B)) return luaL_error(L,"Matrix Sqrt failed\n");
    }
    return 1;
  }
}

/*
 *
 */
static int LuaMatrixEigensystem(lua_State * L)
{
  CompactMatrixType *H;
  uint32_t len = (uint32_t)lua_gettop(L);
  if(len!=1) luaL_error(L, "Matrix.Eigensystem called with %d arguments, only one expected",len);
  H=luaL_checkCompactMatrix(L, 1);
  if(H->Ni!=H->Nj) return luaL_error(L,"Matrix Eigensystem only works on square matrices\n");
  if(QCompactMatrixHermitian(*H, 10*DBL_EPSILON)==0)
  {
    printf("WARNING: non hermitian matrix found\n");
    printf("Using left and right handed eigensystem\n");
    printf("With potential complex eigenvalues\n");
    if (luaL_testudata(L, 1, "Matrix_Type") == NULL)
    {
      CompactMatrixType funL, funR;
      Complex16 *val;
      val= (Complex16*)malloc(H->Ni*sizeof(Complex16));
      NonHermitianCompactMatrixEigensystem(H,&funL, &funR,val);
      LuaPushArray_Complex16(L, val, funL.Ni);
      free(val);
      FreeCompactMatrix(H);
      LuaPushTable_CompactMatrixType(L, &funL);
      FreeCompactMatrix(&funL);
      LuaPushTable_CompactMatrixType(L, &funR);
      FreeCompactMatrix(&funR);
    }
    else
    {
      CompactMatrixType *funL, *funR;
      Complex16 *val;
      val= (Complex16*)malloc(H->Ni*sizeof(Complex16));
      funL = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
      luaL_setmetatable(L, "Matrix_Type");
      funL->OnLuaStack=true;
      funR = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
      luaL_setmetatable(L, "Matrix_Type");
      funR->OnLuaStack=true;
      NonHermitianCompactMatrixEigensystem(H,funL, funR,val);
      LuaPushArray_Complex16(L, val, funL->Ni);
      lua_insert (L, -3);
      free(val);
    }
    return 3;
  }
  else
  {
    if (luaL_testudata(L, 1, "Matrix_Type") == NULL)
    {
      CompactMatrixType fun;
      double *val;
      val= (double*)malloc(H->Ni*sizeof(double));
      CompactMatrixEigensystem(H,&fun,val);
      LuaPushArray_double(L, val, fun.Ni);
      free(val);
      FreeCompactMatrix(H);
      LuaPushTable_CompactMatrixType(L, &fun);
      FreeCompactMatrix(&fun);
    }
    else
    {
      CompactMatrixType *fun,funtemp;
      double *val;
      val= (double*)malloc(H->Ni*sizeof(double));
      CompactMatrixEigensystem(H,&funtemp,val);
      LuaPushArray_double(L, val, funtemp.Ni);
      free(val);
      fun = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
      *fun=funtemp;
      luaL_setmetatable(L, "Matrix_Type");
      fun->OnLuaStack=true;
    }
    return 2;
  }
}



/*
 *
 */
static int LuaMatrixTridiagonalize(lua_State * L)
{
  CompactMatrixType *H,*T;
  ResponsefunctionType *HTri;
  uint32_t len = (uint32_t)lua_gettop(L);
  struct
  {
    CompactMatrixType StartingVectors;
    uint32_t BlockSize;
    uint32_t NTri;
    double   SingularValue;
    bool     Reorthogonalize;
    bool     ReturnU;
  } Options;
  Options.StartingVectors.Ni=0;
  Options.BlockSize=1;
  Options.NTri=UINT32_MAX;
  Options.SingularValue=sqrt(DBL_EPSILON);
  Options.Reorthogonalize=true;
  Options.ReturnU=true;
  LuaOptionsType LuaFunctionOptions[]={
 //   {"StartingVectors"}, "CompactMatrix", &(Options.StartingVectors)},
    {"BlockSize",        "uint32_t",      &(Options.BlockSize)},
    {"NTri",             "uint32_t",      &(Options.NTri)},
    {"SingularValue",    "double",        &(Options.SingularValue)},
    {"ReturnU",          "bool",          &(Options.ReturnU)},
    {NULL, NULL, NULL}};
  if(LuaTestOptionsList(L,-1)==1)
  {
    LuaCheckOptions(L, -1, LuaFunctionOptions);
    len=len-1;
  }
  if(len!=1) luaL_error(L, "Matrix.TriDiagonalize called with %d arguments, only one expected (with possible options on the second position)\n",len);
  H=luaL_checkCompactMatrix(L, 1);
  if(H->Ni!=H->Nj) return luaL_error(L,"Matrix.TriDiagonalize only works on square matrices\n");
  if(QCompactMatrixHermitian(*H, 10*DBL_EPSILON)==0)
  {
    return luaL_error(L,"Matrix.TriDiagonalize only works on Hermitian matrices\nIf the matrix should be Hermitian, but only has a small deviation you can use Chop( (M + ConjugateTranspose(M) )/2 ) as input\n");
  }
  if((Options.StartingVectors.Ni==0)&&(Options.BlockSize==1))
  { // use dense method from Lapack
    HTri = (ResponsefunctionType *) lua_newuserdata(L, sizeof(ResponsefunctionType));
    luaL_setmetatable(L, "ResponseFunction_Type");
    HTri->type='t';
    HTri->Block=false;
    HTri->data=(void*)malloc(sizeof(TriDiagonalMatrixType));
    HTri->NSuggestedGridPoints=0;
    HTri->SuggestedGrid=NULL;
    if(Options.ReturnU)
    {
      T = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
      luaL_setmetatable(L, "Matrix_Type");
      T->OnLuaStack=true;
    }
    else
    {
      T=NULL;
    }
    if(CompactMatrixTriDiagonalize(H, T, (TriDiagonalMatrixType*)(HTri->data)))
    {
      return luaL_error(L,"CompactMatrixTriDiagonalize failed\n");
    }
    if(Options.ReturnU)
    {
      return 2;
    }
    else
    {
      return 1;
    }
  }
  else
  { // use sparse method, i.e. Lanczos
    HTri = (ResponsefunctionType *) lua_newuserdata(L, sizeof(ResponsefunctionType));
    luaL_setmetatable(L, "ResponseFunction_Type");
    HTri->type='t';
    HTri->Block=true;
    HTri->data=(void*)malloc(sizeof(BlockTriDiagonalMatrixType));
    HTri->NSuggestedGridPoints=0;
    HTri->SuggestedGrid=NULL;
    if(Options.NTri==UINT32_MAX) Options.NTri=H->Ni;
    if(Options.StartingVectors.Ni==0)
    {
      Options.StartingVectors.Ni=Options.BlockSize;
      Options.StartingVectors.Nj=H->Nj;
      Options.StartingVectors.QComplex=0;
      InitCompactIdentityMatrix(&(Options.StartingVectors));
    }
    if(CompactMatrixBlockBandDiagonalize(H, &(Options.StartingVectors), (BlockTriDiagonalMatrixType *)(HTri->data), Options.SingularValue, Options.NTri, Options.NTri, Options.Reorthogonalize, Options.ReturnU))
    {
      return luaL_error(L,"CompactMatrixBlockBandDiagonalize failed\n");
    }
    if(Options.ReturnU)
    {
      T = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
      luaL_setmetatable(L, "Matrix_Type");
      *T=Options.StartingVectors;
      T->OnLuaStack=true;
      return 2;
    }
    else
    {
      return 1;
    }
  }
}








/*  on stack Matrix A and Matrix B
 *  computes B^* . A . B^T
 *
 */
static int LuaMatrixRotate(lua_State * L)
{
  CompactMatrixType *A,*B;
  uint32_t len = (uint32_t)lua_gettop(L);
  if(len!=2) luaL_error(L, "Matrix.Rotate called with %d arguments, two arguments expected",len);
  A=luaL_checkCompactMatrix(L, 1);
  B=luaL_checkCompactMatrix(L, 2);
  if(A->Ni!=A->Nj) return luaL_error(L, "Matrix.Rotate needs the first matrix to be square, at the moment Ni=%d and Nj=%d",A->Ni,A->Nj);
  if(A->Ni!=B->Nj) return luaL_error(L, "Matrix.Rotate needs the first matrix to be of the same dimension (%d) as the number of Columns in the second matrix (%d)",A->Ni,B->Nj);
  if ((luaL_testudata(L, 1, "Matrix_Type") == NULL)&&(luaL_testudata(L, 2, "Matrix_Type") == NULL))
  { // matrices on stack as table
    CompactMatrixType C;
    if(TransformCompactMatrixcT(A, B, &C)==1) return luaL_error(L, "Error in Matrix.Rotate");
    LuaPushTable_CompactMatrixType(L, &C);
    FreeCompactMatrix(&C);
  }
  else
  {
    CompactMatrixType *C;
    C = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
    luaL_setmetatable(L, "Matrix_Type");
    if(TransformCompactMatrixcT(A, B, C)==1) return luaL_error(L, "Error in Matrix.Rotate");
    C->OnLuaStack=true;
  }
  if(A->OnLuaStack==false)
  {
    FreeCompactMatrix(A);
    free(A);
  }
  if(B->OnLuaStack==false)
  {
    FreeCompactMatrix(B);
    free(B);
  }
  return 1;
}


/*
 *
 */
static int LuaMatrixExp(lua_State * L)
{
  CompactMatrixType *H;
  uint32_t len = (uint32_t)lua_gettop(L);
  if(len!=1) luaL_error(L, "Matrix.Exp called with %d arguments, only one expected",len);
  H=luaL_checkCompactMatrix(L, 1);
  if(H->Ni!=H->Nj) return luaL_error(L,"Matrix Exp only works on square matrices\n");
  if(QCompactMatrixHermitian(*H, 2*DBL_EPSILON)==0)
  {
    if (luaL_testudata(L, 1, "Matrix_Type") == NULL)
    { // non hermitian matrix stored as table
      CompactMatrixType ExpA;
      if(NonHermitianCompactMatrixExp(H, &ExpA)) return luaL_error(L,"Matrix Exp failed (probably singular needs aditional algorithm to be implemented)\n");
      FreeCompactMatrix(H);
      LuaPushTable_CompactMatrixType(L, &ExpA);
      FreeCompactMatrix(&ExpA);
    }
    else
    { // non hermitian matrix stored as userdata
      CompactMatrixType *ExpA;
      ExpA = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
      luaL_setmetatable(L, "Matrix_Type");
      ExpA->OnLuaStack=true;
      if(NonHermitianCompactMatrixExp(H, ExpA)) return luaL_error(L,"Matrix Exp failed (probably singular needs aditional algorithm to be implemented)\n");
    }
  }
  else
  {
    if (luaL_testudata(L, 1, "Matrix_Type") == NULL)
    { // Hermitian matrix stored as table
      CompactMatrixType ExpA;
      if(HermitianCompactMatrixExp(H, &ExpA)) return luaL_error(L,"Matrix Exp failed (should not happen for Hermitian matrices)\n");
      FreeCompactMatrix(H);
      LuaPushTable_CompactMatrixType(L, &ExpA);
      FreeCompactMatrix(&ExpA);
    }
    else
    {
      CompactMatrixType *ExpA;
      ExpA = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
      luaL_setmetatable(L, "Matrix_Type");
      ExpA->OnLuaStack=true;
      if(HermitianCompactMatrixExp(H, ExpA)) return luaL_error(L,"Matrix Exp failed (should not happen for Hermitian matrices)\n");
    }
  }
  return 1;
}

/*  function Matrix.Random[{min,max},{Ni,Nj}]
 *  returns an Ni by Nj matrix with random numbers in the interval min,max
 *  Possible third input options {{"Userdata",true}}
 */
static int LuaMatrixRandom(lua_State * L)
{
  uint32_t i,j,k,*Nij;
  CompactMatrixType *res;
  Complex16 *minmax;
  bool Userdata=true;
  uint32_t len = (uint32_t)lua_gettop(L);
  LuaOptionsType LuaFunctionOptions[]={
    {"Userdata","bool",&(Userdata)},
    {NULL, NULL, NULL}};
  if (len==3)
  {
    LuaCheckOptions(L, -1, LuaFunctionOptions);
  }
  if(!((len==2)||(len==3)))
  {
    return luaL_error(L,"Matrix Random called with %d arguments, 2 or 3 arguments expected\n",len);
  }
  luaL_checktype(L, 1, LUA_TTABLE);
  if(lua_objlen(L, 1)!=2) luaL_error(L,"Matrix Random called with table of length %d at posittion 1, expected {min,max}\n",lua_objlen(L, 1));
  luaL_checktype(L, 2, LUA_TTABLE);
  if(lua_objlen(L, 2)!=2) luaL_error(L,"Matrix Random called with table of length %d at posittion 2, expected {Ni,Nj}\n",lua_objlen(L, 2));
  if(Userdata==true)
  {
    res = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
    luaL_setmetatable(L, "Matrix_Type");
    res->OnLuaStack=true;
  }
  else
  {
    res = (CompactMatrixType*)malloc(sizeof(CompactMatrixType));
  }
  minmax = luaL_checkarray_complex16(L, 1, &len);
  res->QComplex=0;
  if(minmax[0].Im!=0) res->QComplex=1;
  if(minmax[1].Im!=0) res->QComplex=1;
  Nij    = luaL_checkarray_unit32_t(L, 2, &len);
  res->Ni = Nij[0];
  res->Nj = Nij[1];
  free(Nij);
  InitCompactMatrix(res);
  k=0;
  if(res->QComplex==0)
  {
    for(i=0;i<res->Ni;i++)
    {
      for(j=0;j<res->Nj;j++)
      {
        res->VRe[k]=minmax[0].Re+(minmax[1].Re-minmax[0].Re)*((double)rand()/(double)(RAND_MAX));
        k++;
      }
    }
  }
  else
  {
    for(i=0;i<res->Ni;i++)
    {
      for(j=0;j<res->Nj;j++)
      {
        res->V[k].Re=minmax[0].Re+(minmax[1].Re-minmax[0].Re)*((double)rand()/(double)(RAND_MAX));
        res->V[k].Im=minmax[0].Im+(minmax[1].Im-minmax[0].Im)*((double)rand()/(double)(RAND_MAX));
        k++;
      }
    }
  }
  free(minmax);
  if(Userdata==false)
  {
    LuaPushTable_CompactMatrixType(L, res);
    FreeCompactMatrix(res);
    free(res);
  }
  return 1;
}


/*  functoin Matrix.Zero(dim)
 *  returns an dim by dim Matrix filed with zeros
 *  Possible third input options {{"Userdata",true}}
 */
static int LuaMatrixZero(lua_State * L)
{
  CompactMatrixType *res;
  bool Userdata=true;
  uint32_t len = (uint32_t)lua_gettop(L);
  LuaOptionsType LuaFunctionOptions[]={
    {"Userdata","bool",&(Userdata)},
    {NULL, NULL, NULL}};
  if (len==2)
  {
    LuaCheckOptions(L, -1, LuaFunctionOptions);
  }
  if(!((len==1)||(len==2)))
  {
    return luaL_error(L,"Matrix Zero called with %d arguments, 1 or 2 arguments expected\n", len);
  }
  if(Userdata==true)
  {
    res = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
    luaL_setmetatable(L, "Matrix_Type");
    res->OnLuaStack=true;
  }
  else
  {
    res = (CompactMatrixType*)malloc(sizeof(CompactMatrixType));
  }
  res->QComplex=0;
  res->Ni = (uint32_t)luaL_checkinteger(L, 1);
  res->Nj = res->Ni;
  InitCompactMatrix(res);
  if(Userdata==false)
  {
    LuaPushTable_CompactMatrixType(L, res);
    FreeCompactMatrix(res);
    free(res);
  }
  return 1;
}

/*  functoin Matrix.Identity(dim)
 *  returns an dim by dim Identity Matrix
 *  Possible third input options {{"Userdata",true}}
 */
static int LuaMatrixIdentity(lua_State * L)
{
  CompactMatrixType *res;
  bool Userdata=true;
  uint32_t len = (uint32_t)lua_gettop(L);
  LuaOptionsType LuaFunctionOptions[]={
    {"Userdata","bool",&(Userdata)},
    {NULL, NULL, NULL}};
  if (len==2)
  {
    LuaCheckOptions(L, -1, LuaFunctionOptions);
  }
  if(!((len==1)||(len==2)))
  {
    return luaL_error(L,"Matrix Identity called with %d arguments, 1 or 2 arguments expected\n");
  }
  if(Userdata==true)
  {
    res = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
    luaL_setmetatable(L, "Matrix_Type");
    res->OnLuaStack=true;
  }
  else
  {
    res = (CompactMatrixType*)malloc(sizeof(CompactMatrixType));
  }
  res->QComplex=0;
  res->Ni = (uint32_t)luaL_checkinteger(L, 1);
  res->Nj = res->Ni;
	if (res->Ni != 0)
		InitCompactIdentityMatrix(res);
	else
		InitCompactMatrix(res);
  if(Userdata==false)
  {
    LuaPushTable_CompactMatrixType(L, res);
    FreeCompactMatrix(res);
    free(res);
  }
  return 1;
}

/*  functoin Matrix.Diagonal({x1,x2,x3,...,xn})
 *  returns an n by n Matrix with x1 to xn on the diagonal
 *  Possible third input options {{"Userdata",true}}
 */
static int LuaMatrixDiagonal(lua_State * L)
{
  Complex16 *diag;
  CompactMatrixType *res;
  bool Userdata=true;
  uint32_t i,k,len = (uint32_t)lua_gettop(L);
  LuaOptionsType LuaFunctionOptions[]={
    {"Userdata","bool",&(Userdata)},
    {NULL, NULL, NULL}};
  if (len==2)
  {
    LuaCheckOptions(L, -1, LuaFunctionOptions);
  }
  if(!((len==1)||(len==2)))
  {
    return luaL_error(L,"Matrix Diagonal called with %d arguments, 1 or 2 arguments expected\n");
  }
  if(Userdata==true)
  {
    res = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
    luaL_setmetatable(L, "Matrix_Type");
    res->OnLuaStack=true;
  }
  else
  {
    res = (CompactMatrixType*)malloc(sizeof(CompactMatrixType));
  }
  diag = luaL_checkarray_complex16(L, 1, &len);
  res->QComplex=0;
  for(i=0;i<len;i++)
  {
    if(diag[i].Im!=0)
    {
      res->QComplex=1;
      break;
    }
  }
  res->Ni = len;
  res->Nj = len;
  InitCompactIdentityMatrix(res);
  k=0;
  if(res->QComplex==0)
  {
    for(i=0;i<res->Ni;i++)
    {
      res->VRe[k]=diag[i].Re;
      k+=res->Ni+1;
    }
  }
  else
  {
    for(i=0;i<res->Ni;i++)
    {
      res->V[k].Re=diag[i].Re;
      res->V[k].Im=diag[i].Im;
      k+=res->Ni+1;
    }
  }
  if(Userdata==false)
  {
    LuaPushTable_CompactMatrixType(L, res);
    FreeCompactMatrix(res);
    free(res);
  }
  free(diag);	//Added by Martin, otherwise heap is lost
  return 1;
}


/*
 *
 */
static int LuaMatrixToTable(lua_State * L)
{
  CompactMatrixType *A;
  A = (CompactMatrixType *) luaL_checkCompactMatrix(L, 1);
  LuaPushTable_CompactMatrixType(L, A);
  FreeCompactMatrix(A);
  return 1;
}

/*
 *
 */
static int LuaMatrixToUserdata(lua_State * L)
{
  CompactMatrixType *A,*B;
  A = (CompactMatrixType *) luaL_checkCompactMatrix(L, 1);
  B = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
  luaL_setmetatable(L, "Matrix_Type");
  if (luaL_testudata(L, 1, "Matrix_Type") == NULL)
  {
    *B=*A;
    free(A);
  }
  else
  {
    CompactMatrixCopy(*A, B);
  }
  B->OnLuaStack=true;
  return 1;
}

/* pushes -M to stack
 *
 */
static int LuaMatrixUnm(lua_State * L)
{
  CompactMatrixType *A;
  A = (CompactMatrixType *) luaL_checkCompactMatrix(L, 1);
  if (luaL_testudata(L, 1, "Matrix_Type") == NULL)
  {
    CompactMatrixScale(A, -1);
    LuaPushTable_CompactMatrixType(L, A);
    FreeCompactMatrix(A);
    free(A);
  }
  else
  {
    CompactMatrixType *B;
    B = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
    luaL_setmetatable(L, "Matrix_Type");
    B->OnLuaStack=true;
    CompactMatrixCopy(*A, B);
    CompactMatrixScale(B, -1);
  }
  return 1;
}

/*
 *
 */
static int LuaMatrixAdd(lua_State * L)
{
  CompactMatrixType *A,*B,*C;
  A = (CompactMatrixType *) luaL_checkCompactMatrix(L, 1);
  B = (CompactMatrixType *) luaL_checkCompactMatrix(L, 2);
  if((A->Ni != B->Ni) || (A->Nj != B->Nj))
  {
    return luaL_error(L, "Error: dimensions of matrix A (%d %d) and B (%d %d) in A+B not equivalent\n",A->Ni,A->Nj,B->Ni,B->Nj);
  }
  if (luaL_testudata(L, 1, "Matrix_Type") == NULL)
  { // A is read from a table
    CompactMatrixAdd(A, B);
    if (luaL_testudata(L, 2, "Matrix_Type") == NULL)
    { // B is read from a table
      LuaPushTable_CompactMatrixType(L, A);
      FreeCompactMatrix(A);
      free(A);
      FreeCompactMatrix(B);
      free(B);
    }
    else
    {
      C = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
      luaL_setmetatable(L, "Matrix_Type");
      *C=*A;
      C->OnLuaStack=true;
      free(A);
    }
  }
  else if (luaL_testudata(L, 2, "Matrix_Type") == NULL)
  { // B is read from a table
    CompactMatrixAdd(B, A);
    C = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
    luaL_setmetatable(L, "Matrix_Type");
    *C=*B;
    C->OnLuaStack=true;
    free(B);
  }
  else
  {
    C = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
    luaL_setmetatable(L, "Matrix_Type");
    C->OnLuaStack=true;
    CompactMatrixCopy(*A, C);
    CompactMatrixAdd(C, B);
  }
  return 1;
}

/*
 *
 */
static int LuaMatrixSub(lua_State * L)
{
  CompactMatrixType *A,*B,*C;
  A = (CompactMatrixType *) luaL_checkCompactMatrix(L, 1);
  B = (CompactMatrixType *) luaL_checkCompactMatrix(L, 2);
  if((A->Ni != B->Ni) || (A->Nj != B->Nj))
  {
    return luaL_error(L, "Error: dimensions of matrix A (%d %d) and B (%d %d) in A-B not equivalent\n",A->Ni,A->Nj,B->Ni,B->Nj);
  }
  if (luaL_testudata(L, 1, "Matrix_Type") == NULL)
  { // A is read from a table
    CompactMatrixSub(A, B);
    if (luaL_testudata(L, 2, "Matrix_Type") == NULL)
    { // B is read from a table
      LuaPushTable_CompactMatrixType(L, A);
      FreeCompactMatrix(A);
      free(A);
      FreeCompactMatrix(B);
      free(B);
    }
    else
    {
      C = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
      luaL_setmetatable(L, "Matrix_Type");
      *C=*A;
      C->OnLuaStack=true;
      free(A);
    }
  }
  else if (luaL_testudata(L, 2, "Matrix_Type") == NULL)
  { // B is read from a table
    CompactMatrixSub(B, A);
    CompactMatrixScale(B, -1);
    C = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
    luaL_setmetatable(L, "Matrix_Type");
    *C=*B;
    C->OnLuaStack=true;
    free(B);
  }
  else
  {
    C = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
    luaL_setmetatable(L, "Matrix_Type");
    C->OnLuaStack=true;
    CompactMatrixCopy(*A, C);
    CompactMatrixSub(C, B);
  }
  return 1;
}

/*
 *
 */
static int LuaMatrixMul(lua_State * L)
{
  switch (lua_type(L,1))
  {
    case LUA_TNUMBER:
    {
      // (number) operator (matrix)
      double a;
      CompactMatrixType *B,*C;
      a = luaL_checknumber(L, 1);
      B = (CompactMatrixType *) luaL_checkCompactMatrix(L, 2);
      if (luaL_testudata(L, 2, "Matrix_Type") == NULL)
      { // Matrix on stack as table
        CompactMatrixScale(B, a);
        LuaPushTable_CompactMatrixType(L, B);
        FreeCompactMatrix(B);
        free(B);
      }
      else
      { // Matrix on stack as userdata
        C = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
        luaL_setmetatable(L, "Matrix_Type");
        C->OnLuaStack=true;
        CompactMatrixCopy(*B, C);
        CompactMatrixScale(C, a);
      }
      return 1;
    }
    case LUA_TUSERDATA:
    {
      if (luaL_testudata(L, 1, "Complex_Type") != NULL)
      {
        // (complex number) operator (matrix)
        Complex16 *a;
        CompactMatrixType *B,*C;
        a = (Complex16 *) luaL_checkudata(L, 1, "Complex_Type");
        B = (CompactMatrixType *) luaL_checkCompactMatrix(L, 2);
        if (luaL_testudata(L, 2, "Matrix_Type") == NULL)
        { // Matrix on stack as table
          CompactMatrixScaleComplex(B, *a);
          LuaPushTable_CompactMatrixType(L, B);
          FreeCompactMatrix(B);
          free(B);
        }
        else
        { // Matrix on stack as userdata
          C = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
          luaL_setmetatable(L, "Matrix_Type");
          C->OnLuaStack=true;
          CompactMatrixCopy(*B, C);
          CompactMatrixScaleComplex(C, *a);
        }
        return 1;
      }
      if (luaL_testudata(L, 1, "Matrix_Type") != NULL)
      {
        // (matrix as user data) operator (???)
        switch (lua_type(L,2))
        {
          case LUA_TNUMBER:
          {
            // (matrix as userdata) operator (number)
            double b;
            CompactMatrixType *A,*C;
            A = (CompactMatrixType*)luaL_checkudata(L, 1, "Matrix_Type");
            b = luaL_checknumber(L, 2);
            C = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
            luaL_setmetatable(L, "Matrix_Type");
            C->OnLuaStack=true;
            CompactMatrixCopy(*A, C);
            CompactMatrixScale(C, b);
            return 1;
          }
          case LUA_TUSERDATA:
          {
            if (luaL_testudata(L, 2, "Complex_Type") != NULL)
            {
              // (matrix as user data) operator (complex number)
              Complex16 *b;
              CompactMatrixType *A,*C;
              A = (CompactMatrixType*)luaL_checkudata(L, 1, "Matrix_Type");
              b = (Complex16 *) luaL_checkudata(L, 2, "Complex_Type");
              C = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
              luaL_setmetatable(L, "Matrix_Type");
              C->OnLuaStack=true;
              CompactMatrixCopy(*A, C);
              CompactMatrixScaleComplex(C, *b);
              return 1;
            }
            if (luaL_testudata(L, 2, "Matrix_Type") != NULL)
            {
              // (matrix as use data) operator (matrix as user data)
              CompactMatrixType *A,*B,*C;
              A = (CompactMatrixType*)luaL_checkudata(L, 1, "Matrix_Type");
              B = (CompactMatrixType*)luaL_checkudata(L, 2, "Matrix_Type");
              C = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
              luaL_setmetatable(L, "Matrix_Type");
              C->OnLuaStack=true;
              if(A->Nj!=B->Ni) return luaL_error(L,"Matrix multiplication requirs the number of columns of the left matrix (%d) to be equal to the number of rows (%d) of the right matrix\n",A->Nj,B->Ni);
              CompactMatrixMultiply('N', 'N', A, B, C);
              return 1;
            }
          }
          case LUA_TTABLE:
          {
            // (matrix as user data) operator (matrix as table)
            CompactMatrixType *A,*B,*C;
            A = (CompactMatrixType*) luaL_checkudata(L, 1, "Matrix_Type");
            B = (CompactMatrixType *) luaL_checkCompactMatrix(L, 2);
            C = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
            luaL_setmetatable(L, "Matrix_Type");
            C->OnLuaStack=true;
            if(A->Nj!=B->Ni) return luaL_error(L,"Matrix multiplication requirs the number of columns of the left matrix (%d) to be equal to the number of rows (%d) of the right matrix\n",A->Nj,B->Ni);
            CompactMatrixMultiply('N', 'N', A, B, C);
            FreeCompactMatrix(B);
            free(B);
            return 1;
          }
          default:
          {
            return luaL_error(L,"Unrecognized type (%s) at position two of multiplication \n",lua_type(L,1));
          }
        }
      }
    }
    case LUA_TTABLE:
    {
      // (matrix as table) operator (???)
      switch (lua_type(L,2))
      {
        case LUA_TNUMBER:
        {
          // (matrix as table) operator (number)
          double b;
          CompactMatrixType *A;
          A = (CompactMatrixType *) luaL_checkCompactMatrix(L, 1);
          b = luaL_checknumber(L, 2);
          CompactMatrixScale(A, b);
          LuaPushTable_CompactMatrixType(L, A);
          FreeCompactMatrix(A);
          free(A);
          return 1;
        }
        case LUA_TUSERDATA:
        {
          if (luaL_testudata(L, 2, "Complex_Type") != NULL)
          {
            // (matrix as table) operator (complex number)
            Complex16 *b;
            CompactMatrixType *A;
            A = (CompactMatrixType *) luaL_checkCompactMatrix(L, 1);
            b = (Complex16 *) luaL_checkudata(L, 2, "Complex_Type");
            if (luaL_testudata(L, 2, "Matrix_Type") == NULL)
            CompactMatrixScaleComplex(A, *b);
            LuaPushTable_CompactMatrixType(L, A);
            FreeCompactMatrix(A);
            free(A);
            return 1;
          }
          if (luaL_testudata(L, 2, "Matrix_Type") != NULL)
          {
            // (matrix as table) operator (matrix as user data)
            CompactMatrixType *A,*B,*C;
            A = (CompactMatrixType *) luaL_checkCompactMatrix(L, 1);
            B = (CompactMatrixType *) (CompactMatrixType*)luaL_checkudata(L, 2, "Matrix_Type");
            C = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
            luaL_setmetatable(L, "Matrix_Type");
            C->OnLuaStack=true;
            if(A->Nj!=B->Ni) return luaL_error(L,"Matrix multiplication requirs the number of columns of the left matrix (%d) to be equal to the number of rows (%d) of the right matrix\n",A->Nj,B->Ni);
            CompactMatrixMultiply('N', 'N', A, B, C);
            FreeCompactMatrix(A);
            free(A);
            return 1;
          }
        }
        case LUA_TTABLE:
        {
          // (matrix as table) operator (matrix as table)
          CompactMatrixType *A,*B,C;
          A = (CompactMatrixType *) luaL_checkCompactMatrix(L, 1);
          B = (CompactMatrixType *) luaL_checkCompactMatrix(L, 2);
          if(A->Nj!=B->Ni) return luaL_error(L,"Matrix multiplication requirs the number of columns of the left matrix (%d) to be equal to the number of rows (%d) of the right matrix\n",A->Nj,B->Ni);
          CompactMatrixMultiply('N', 'N', A, B, &C);
          LuaPushTable_CompactMatrixType(L, &C);
          FreeCompactMatrix(A);
          free(A);
          FreeCompactMatrix(B);
          free(B);
          FreeCompactMatrix(&C);
          return 1;
        }
        default:
        {
          return luaL_error(L,"Unrecognized type (%s) at position two of multiplication \n",lua_type(L,1));
        }
      }
    }
    default:
    {
      return luaL_error(L,"Unrecognized type (%d) at position one of multiplication \n",lua_type(L,1));
    }
  }
}

/*
 *
 */
static int LuaMatrixDiv(lua_State * L)
{
  switch (lua_type(L,1))
  {
    case LUA_TNUMBER:
    {
      // (number) operator (matrix)
      double a;
      CompactMatrixType *B,*C;
      a = luaL_checknumber(L, 1);
      B = (CompactMatrixType *) luaL_checkCompactMatrix(L, 2);
      if (luaL_testudata(L, 2, "Matrix_Type") == NULL)
      { // Matrix on stack as table
        if(B->Nj!=B->Ni) return luaL_error(L,"Matrix inversion only works for square matrices dimensions are (%d %d)\n",B->Ni,B->Nj);
        if(CompactMatrixInverse(B)!=0) return luaL_error(L,"Matrix inversion failed, probably the matrix is singular\n");
        CompactMatrixScale(B, a);
        LuaPushTable_CompactMatrixType(L, B);
        FreeCompactMatrix(B);
        free(B);
      }
      else
      { // Matrix on stack as userdata
        C = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
        luaL_setmetatable(L, "Matrix_Type");
        C->OnLuaStack=true;
        CompactMatrixCopy(*B, C);
        if(C->Nj!=C->Ni) return luaL_error(L,"Matrix inversion only works for square matrices dimensions are (%d %d)\n",C->Ni,C->Nj);
        if(CompactMatrixInverse(C)!=0) return luaL_error(L,"Matrix inversion failed, probably the matrix is singular\n");
        CompactMatrixScale(C, a);
      }
      return 1;
    }
    case LUA_TUSERDATA:
    {
      if (luaL_testudata(L, 1, "Complex_Type") != NULL)
      {
        // (complex number) operator (matrix)
        Complex16 *a;
        CompactMatrixType *B,*C;
        a = (Complex16 *) luaL_checkudata(L, 1, "Complex_Type");
        B = (CompactMatrixType *) luaL_checkCompactMatrix(L, 2);
        if (luaL_testudata(L, 2, "Matrix_Type") == NULL)
        { // Matrix on stack as table
          if(B->Nj!=B->Ni) return luaL_error(L,"Matrix inversion only works for square matrices dimensions are (%d %d)\n",B->Ni,B->Nj);
          if(CompactMatrixInverse(B)!=0) return luaL_error(L,"Matrix inversion failed, probably the matrix is singular\n");
          CompactMatrixScaleComplex(B, *a);
          LuaPushTable_CompactMatrixType(L, B);
          FreeCompactMatrix(B);
          free(B);
        }
        else
        { // Matrix on stack as userdata
          C = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
          luaL_setmetatable(L, "Matrix_Type");
          C->OnLuaStack=true;
          CompactMatrixCopy(*B, C);
          if(C->Nj!=C->Ni) return luaL_error(L,"Matrix inversion only works for square matrices dimensions are (%d %d)\n",C->Ni,C->Nj);
          if(CompactMatrixInverse(C)!=0) return luaL_error(L,"Matrix inversion failed, probably the matrix is singular\n");
          CompactMatrixScaleComplex(C, *a);
        }
        return 1;
      }
      if (luaL_testudata(L, 1, "Matrix_Type") != NULL)
      {
        // (matrix as user data) operator (???)
        switch (lua_type(L,2))
        {
          case LUA_TNUMBER:
          {
            // (matrix as userdata) operator (number)
            double b;
            CompactMatrixType *A,*C;
            A = (CompactMatrixType*)luaL_checkudata(L, 1, "Matrix_Type");
            b = luaL_checknumber(L, 2);
            C = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
            luaL_setmetatable(L, "Matrix_Type");
            C->OnLuaStack=true;
            CompactMatrixCopy(*A, C);
            CompactMatrixScale(C, 1/b);
            return 1;
          }
          case LUA_TUSERDATA:
          {
            if (luaL_testudata(L, 2, "Complex_Type") != NULL)
            {
              // (matrix as user data) operator (complex number)
              Complex16 *b,binv;
              CompactMatrixType *A,*C;
              A = (CompactMatrixType*)luaL_checkudata(L, 1, "Matrix_Type");
              b = (Complex16 *) luaL_checkudata(L, 2, "Complex_Type");
              C = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
              luaL_setmetatable(L, "Matrix_Type");
              C->OnLuaStack=true;
              CompactMatrixCopy(*A, C);
              binv.Re = b->Re / (Sqr(b->Re) + Sqr(b->Im));
              binv.Im = b->Im / (Sqr(b->Re) + Sqr(b->Im));
              CompactMatrixScaleComplex(C, binv);
              return 1;
            }
            if (luaL_testudata(L, 2, "Matrix_Type") != NULL)
            {
              // (matrix as use data) operator (matrix as user data)
              CompactMatrixType *A,*B,*C,Binv;
              A = (CompactMatrixType*)luaL_checkudata(L, 1, "Matrix_Type");
              B = (CompactMatrixType*)luaL_checkudata(L, 2, "Matrix_Type");
              C = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
              luaL_setmetatable(L, "Matrix_Type");
              C->OnLuaStack=true;
              CompactMatrixCopy(*B, &Binv);
              if(Binv.Nj!=Binv.Ni) return luaL_error(L,"Matrix inversion only works for square matrices dimensions are (%d %d)\n",Binv.Ni,Binv.Nj);
              if(CompactMatrixInverse(&Binv)!=0) return luaL_error(L,"Matrix inversion failed, probably the matrix is singular\n");
              if(A->Nj!=B->Ni) return luaL_error(L,"Matrix multiplication requirs the number of columns of the left matrix (%d) to be equal to the number of rows (%d) of the right matrix\n",A->Nj,B->Ni);
              CompactMatrixMultiply('N', 'N', A, &Binv, C);
              FreeCompactMatrix(&Binv);
              return 1;
            }
          }
          case LUA_TTABLE:
          {
            // (matrix as user data) operator (matrix as table)
            CompactMatrixType *A,*B,*C;
            A = (CompactMatrixType*) luaL_checkudata(L, 1, "Matrix_Type");
            B = (CompactMatrixType *) luaL_checkCompactMatrix(L, 2);
            C = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
            luaL_setmetatable(L, "Matrix_Type");
            C->OnLuaStack=true;
            if(B->Nj!=B->Ni) return luaL_error(L,"Matrix inversion only works for square matrices dimensions are (%d %d)\n",B->Ni,B->Nj);
            if(CompactMatrixInverse(B)!=0) return luaL_error(L,"Matrix inversion failed, probably the matrix is singular\n");
            if(A->Nj!=B->Ni) return luaL_error(L,"Matrix multiplication requirs the number of columns of the left matrix (%d) to be equal to the number of rows (%d) of the right matrix\n",A->Nj,B->Ni);
            CompactMatrixMultiply('N', 'N', A, B, C);
            FreeCompactMatrix(B);
            free(B);
            return 1;
          }
          default:
          {
            return luaL_error(L,"Unrecognized type (%s) at position two of multiplication \n",lua_type(L,1));
          }
        }
      }
    }
    case LUA_TTABLE:
    {
      // (matrix as table) operator (???)
      switch (lua_type(L,2))
      {
        case LUA_TNUMBER:
        {
          // (matrix as table) operator (number)
          double b;
          CompactMatrixType *A;
          A = (CompactMatrixType *) luaL_checkCompactMatrix(L, 1);
          b = luaL_checknumber(L, 2);
          CompactMatrixScale(A, 1/b);
          LuaPushTable_CompactMatrixType(L, A);
          FreeCompactMatrix(A);
          free(A);
          return 1;
        }
        case LUA_TUSERDATA:
        {
          if (luaL_testudata(L, 2, "Complex_Type") != NULL)
          {
            // (matrix as table) operator (complex number)
            Complex16 *b,binv;
            CompactMatrixType *A;
            A = (CompactMatrixType *) luaL_checkCompactMatrix(L, 1);
            b = (Complex16 *) luaL_checkudata(L, 2, "Complex_Type");
            binv.Re = b->Re / (Sqr(b->Re) + Sqr(b->Im));
            binv.Im = b->Im / (Sqr(b->Re) + Sqr(b->Im));
            CompactMatrixScaleComplex(A, binv);
            LuaPushTable_CompactMatrixType(L, A);
            FreeCompactMatrix(A);
            free(A);
            return 1;
          }
          if (luaL_testudata(L, 2, "Matrix_Type") != NULL)
          {
            // (matrix as table) operator (matrix as user data)
            CompactMatrixType *A,*B,*C,Binv;
            A = (CompactMatrixType *) luaL_checkCompactMatrix(L, 1);
            B = (CompactMatrixType *) (CompactMatrixType*)luaL_checkudata(L, 2, "Matrix_Type");
            C = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
            luaL_setmetatable(L, "Matrix_Type");
            C->OnLuaStack=true;
            CompactMatrixCopy(*B, &Binv);
            if(Binv.Nj!=Binv.Ni) return luaL_error(L,"Matrix inversion only works for square matrices dimensions are (%d %d)\n",Binv.Ni,Binv.Nj);
            if(CompactMatrixInverse(&Binv)!=0) return luaL_error(L,"Matrix inversion failed, probably the matrix is singular\n");
            if(A->Nj!=B->Ni) return luaL_error(L,"Matrix multiplication requirs the number of columns of the left matrix (%d) to be equal to the number of rows (%d) of the right matrix\n",A->Nj,B->Ni);
            CompactMatrixMultiply('N', 'N', A, B, C);
            FreeCompactMatrix(A);
            free(A);
            FreeCompactMatrix(&Binv);
            return 1;
          }
        }
        case LUA_TTABLE:
        {
          // (matrix as table) operator (matrix as table)
          CompactMatrixType *A,*B,C;
          A = (CompactMatrixType *) luaL_checkCompactMatrix(L, 1);
          B = (CompactMatrixType *) luaL_checkCompactMatrix(L, 2);
          if(B->Nj!=B->Ni) return luaL_error(L,"Matrix inversion only works for square matrices dimensions are (%d %d)\n",B->Ni,B->Nj);
          if(CompactMatrixInverse(B)!=0) return luaL_error(L,"Matrix inversion failed, probably the matrix is singular\n");
          if(A->Nj!=B->Ni) return luaL_error(L,"Matrix multiplication requirs the number of columns of the left matrix (%d) to be equal to the number of rows (%d) of the right matrix\n",A->Nj,B->Ni);
          CompactMatrixMultiply('N', 'N', A, B, &C);
          LuaPushTable_CompactMatrixType(L, &C);
          FreeCompactMatrix(A);
          free(A);
          FreeCompactMatrix(B);
          free(B);
          FreeCompactMatrix(&C);
          return 1;
        }
        default:
        {
          return luaL_error(L,"Unrecognized type (%s) at position two of multiplication \n",lua_type(L,1));
        }
      }
    }
    default:
    {
      return luaL_error(L,"Unrecognized type (%d) at position one of multiplication \n",lua_type(L,1));
    }
  }
}

/*
 *
 */
static int LuaMatrixPow(lua_State * L)
{
  double n;
  CompactMatrixType *A,*B;
  n = luaL_checknumber(L, 2);
  if(Abs(n-rint(n))>(10*DBL_EPSILON*Abs(n)+10*DBL_MIN)) return luaL_error(L,"Matrix power only works on integers\n");
  A = (CompactMatrixType *) luaL_checkCompactMatrix(L, 1);
  if(A->Ni!=A->Nj) return luaL_error(L,"Matrix power only works on square matrices\n");
  if (luaL_testudata(L, 1, "Matrix_Type") == NULL)
  { // matrix on stack as table
    B= (CompactMatrixType*)malloc(sizeof(CompactMatrixType));
    if(CompactMatrixPow(A,B,(int)n)) return luaL_error(L,"Matrix power failed\n");
    LuaPushTable_CompactMatrixType(L, B);
    FreeCompactMatrix(A);
    free(A);
    FreeCompactMatrix(B);
    free(B);
  }
  else
  { // matrix on stack as userdata
    B = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
    luaL_setmetatable(L, "Matrix_Type");
    B->OnLuaStack=true;
    if(CompactMatrixPow(A,B,(int)n)) return luaL_error(L,"Matrix power failed\n");
  }
  return 1;
}


/*
 *
 */
static int LuaMatrixToString(lua_State * L)
{
	CompactMatrixType *A;
	A = (CompactMatrixType *) luaL_checkCompactMatrix(L, 1);
	std::stringstream sstream;
	int colWidth;
	for (int i = 0; i < A->Ni; ++i)
	{
		if (i == 0)
			sstream << "{ { ";
		else
			sstream << "  { ";
		for (int j = 0; j < A->Nj; ++j)
		{
			if (j != 0)
				sstream << ", ";
			if (A->QComplex)
			{
				colWidth = 16;
				string tmpStr = numberToStringMaxWidth(A->V[j + i * (A->Nj)].Re, colWidth / 2 - 1);
        if (A->V[j + i * (A->Nj)].Im == 0)
        {
          sstream << tmpStr << "  ";
        }
        else
        {
          if (A->V[j + i * (A->Nj)].Im < 0)
            tmpStr += " -";
          else
            tmpStr += " +";
          tmpStr += numberToStringMaxWidth(Abs(A->V[j + i * (A->Nj)].Im), colWidth - (int)tmpStr.length() - 1);
          sstream << tmpStr << " I";
        }
        sstream << blankFill(colWidth, tmpStr.length());
			}
      else
			{
				colWidth = 8;
				string tmpStr = numberToStringMaxWidth(A->VRe[j + i * (A->Nj)], colWidth - 1);
				sstream << tmpStr;
				sstream << blankFill(colWidth, tmpStr.length());
			}
		}
		if (i == A->Ni - 1)
			sstream << "} }" << endl;
		else
			sstream << "} ," << endl;
	}
	lua_pushstring(L, sstream.str());
	return 1;
}

static int LuaMatrixElem(lua_State *L)
{
	CompactMatrixType M;
	int i=0, j=0;
	int len = lua_gettop(L);
	if (len != 2) return luaL_error(L, "Error: Matrix.Elem called with %d arguments, 2 expected: i and j.", len);
	if (!luaL_isType(L, 1, i)) return luaL_error(L, "Error: First element in Matrix.Elem should be an integer!");
	if (!luaL_isType(L, 2, j)) return luaL_error(L, "Error: Second element in Matrix.Elem should be an integer!");
	lua_check(L, lua_upvalueindex(1), M);
	lua_check(L, 1, i);
	lua_check(L, 2, j);
	if (abs(i) > M.Ni || i == 0)
	{
		return luaL_error(L, "Error in Matrix.Elem: Index i=%d lies outside of range +-[1,%d].", i, M.Ni);
	}
	if (abs(j) > M.Nj || j == 0)
	{
		return luaL_error(L, "Error in Matrix.Elem: Index j=%d lies outside of range +-[1,%d].", j, M.Nj);
	}
	if (i < 0)
		i = M.Ni + i;
	else
		--i;
	if (j < 0)
		j = M.Nj + j;
	else
		--j;
	if ((bool) M.QComplex)
		lua_push(L, M.V[i * M.Nj + j]);
	else
		lua_push(L, M.VRe[i * M.Nj + j]);
	return 1;
}

/*  Index can be either: an integer in the range 1 to Ni or from -1 to -Ni.
 *                       a table of integers in the set [-Ni, ... ,-1,1,..., Ni]
 *                          -- these two options return an MatrixRow object
 *                       the string "Ni", "Nj", "NRow", "NCol", "NColumn", "dim", "Dim", "elem", "Elem"
 *                          -- these return integers ("Ni", "Nj", "NRow", "NCol", "NColumn")
 *                                          a table of two integers ("dim", "Dim")
 *                                          a double or complex value ("elem", "Elem")
 *
 *   On stack we have at position 1 the matrix (for sure)
 *   On position two the index --  which needs to be checked.
 *                                 wrong indices thorugh an error
 */
static int LuaMatrixIndex(lua_State *L)
{
  CompactMatrixType* A = (CompactMatrixType*)luaL_checkudata(L, 1, "Matrix_Type");
  switch (lua_type(L,2))
  {
    case LUA_TNUMBER:
    {
      int i = luaL_checkint(L,2);
      if(i==0) return luaL_error(L, "Error: Matrix indexed with zero. Indices must be in the range [-Nrow(%d), ..., -1,1,...,Nrow] in the matrix\n",A->Ni);
      if(i> (int)A->Ni) return luaL_error(L, "Error: Matrix indexed with number %d. Indices can not be larger than the number of rows (%d) in the matrix\n",i,A->Ni);
      if(i<-((int)A->Ni)) return luaL_error(L, "Error: Matrix indexed with number %d. Indices can not be smaller than minus the number of rows (%d) in the matrix\n",i,A->Ni);
      // we now return a table with the elements __Matrix that contains a pointer to the original matrix and the elment __Row that contains the number i
      //                       this table has the meta method MatrixRowMeta set
      lua_createtable(L, 0, 2);
      lua_pushvalue(L, 1);
      lua_setfield(L,-2,"__Matrix");
      lua_pushnumber(L, i);
      lua_setfield(L,-2,"__Row");
      lua_getglobal(L, "MatrixRowMeta");
      lua_setmetatable(L, -2);
      return 1;
    }
    break;
    case LUA_TTABLE:
    {
      int i,j;
      int len = (int)lua_objlen(L, 2);
      lua_createtable(L, 0, 2);
      lua_pushvalue(L, 1);
      lua_setfield(L,-2,"__Matrix");
      lua_createtable(L, len, 0);
      for(i=1;i<=len;i++)
      {
        // get number 1 from the table at position 2 of the stack.
        // check if the element in the table containgin the indices are in the range [-Ni,...,-1,1,...Ni]
        // push the elment to the table at position -1 on the stack
        lua_rawgeti(L,2, i);
        j=(int)luaL_checkint(L, -1);
        if(j ==   0) return luaL_error(L, "Error: Matrix indexed with zero. Indices must be in the range [-Nrow(%d), ..., -1,1,...,Nrow] in the matrix\n",A->Ni);
        if(j<-((int)A->Ni)) return luaL_error(L, "Error: Matrix indexed with number %d. Indices can not be smaller than minus the number of rows (%d) in the matrix\n",j,A->Ni);
        if(j> (int)A->Ni) return luaL_error(L, "Error: Matrix indexed with number %d. Indices can not be larger than the number of rows (%d) in the matrix\n",j,A->Ni);
        lua_rawseti(L, -2, i);
      }
      lua_setfield(L,-2,"__Row");
      lua_getglobal(L, "MatrixRowMeta");
      lua_setmetatable(L, -2);
      return 1;
    }
    break;
    case LUA_TSTRING:
    {
      // strings allowed: "Ni"=="NRow", "Nj"=="Ncol"=="Ncolumn", "dim"=="Dim", "elem"=="Elem",
      const char *lua_string = lua_tostring(L, 2);
      if ((strcmp(lua_string, "Ni") == 0)||(strcmp(lua_string, "NRow") == 0)||(strcmp(lua_string, "Nrow") == 0))
      {
        lua_pushnumber(L, A->Ni);
        return 1;
      }
      if ((strcmp(lua_string, "Nj") == 0)||(strcmp(lua_string, "NCol") == 0)||(strcmp(lua_string, "Ncol") == 0)||(strcmp(lua_string, "NColumn") == 0)||(strcmp(lua_string, "Ncolumn") == 0))
      {
        lua_pushnumber(L, A->Nj);
        return 1;
      }
      if ((strcmp(lua_string, "Dim") == 0)||(strcmp(lua_string, "dim") == 0))
      {
        lua_createtable(L, 2, 0);
        lua_pushnumber(L, A->Ni);
        lua_rawseti(L, -2, 1);
        lua_pushnumber(L, A->Nj);
        lua_rawseti(L, -2, 2);
        return 1;
      }
      if  ((strcmp(lua_string, "Elem") == 0)||(strcmp(lua_string, "elem") == 0))
      {
        lua_pushvalue(L, 1);
        lua_pushcclosure(L, LuaMatrixElem, 1);
        return 1;
      }
      return luaL_error(L, "Error: Matrix indexed with unrecognised string %s",lua_string);
    }
    break;
    default:
      return luaL_error(L, "Error: Matrix indexed with unrecognised type!");
    break;
  }
}

/* the first position on the Lua stack is a userdata of Matrix_Type
 *   index must be an integer or a table of integers
 *   the integers must be unique and in the range [-Ni, ... ,-1,1,...,Ni]
 *   or a string Ni, Nrow, Nj, Ncol
 *
 *   the right hand side must be a table of numbers (real or complex) with the length equal to Nj
 *      or a table of tables of numbers with the  number of rows equal to the nubmers in the index and the length equal to Nj.
 *
 */
static int LuaMatrixNewIndex(lua_State *L)
{
  CompactMatrixType* A = (CompactMatrixType*)luaL_checkudata(L, 1, "Matrix_Type");
  switch (lua_type(L,2))
  {
    case LUA_TNUMBER:
    {
      int j;
      int i = luaL_checkint(L,2);
      if(i==0) return luaL_error(L, "Error: Matrix indexed with zero. Indices must be in the range [-Nrow(%d), ..., -1,1,...,Nrow] in the matrix\n",A->Ni);
      if(i> (int)A->Ni) return luaL_error(L, "Error: Matrix indexed with number %d. Indices can not be larger than the number of rows (%d) in the matrix\n",i,A->Ni);
      if(i<-((int)A->Ni)) return luaL_error(L, "Error: Matrix indexed with number %d. Indices can not be smaller than minus the number of rows (%d) in the matrix\n",i,A->Ni);
      // we now read a table of length Nj
      if(lua_type(L,3)!=LUA_TTABLE) return luaL_error(L, "Error: Matrix you try to set row %d equal to somehting that is not a table\n",i);
      int len = (int)lua_objlen(L, 3);
      if(len!=A->Nj) return luaL_error(L, "Error: Matrix you try to set row %d equal to a table of length %d whereas the matrix has %d columns\n",i,len,A->Nj);
      for(j=1;j<=len;j++)
      {
        lua_rawgeti(L, 3, j);
        if(lua_isnumber(L, -1))
        {
          if(A->QComplex==0)
          {
            A->VRe[(i-1)*A->Nj +(j-1)] = (double)luaL_checknumber(L, -1);
          }
          else
          {
            A->V[(i-1)*A->Nj +(j-1)].Re = (double)luaL_checknumber(L, -1);
            A->V[(i-1)*A->Nj +(j-1)].Im = 0;
          }
          lua_pop(L, 1);
        }
        else if (luaL_testudata(L, -1, "Complex_Type") != NULL)
        {
          if(A->QComplex==0)
          {
            MakeCompactRealMatrixComplex(A);
          }
          A->V[(i-1)*A->Nj +(j-1)] = *(Complex16 *)luaL_checkudata(L, -1, "Complex_Type");
          lua_pop(L, 1);
        }
        else
        {
          luaL_error(L, "Error in reading table at position %d real or complex number expected \n",j);
        }
      }
      return 0;
    }
    break;
    case LUA_TTABLE:
    {
      //TODO: Matrix indexed with a table set equal to a table of tables or matrix.
      return luaL_error(L, "Error: Matrix indexed with a table. Not implemented");
      return 0;
    }
    break;
    case LUA_TSTRING:
    {
      return luaL_error(L, "Error: Matrix indexed with a string assigned to a value. That is not an option.");
    }
    break;
    default:
      return luaL_error(L, "Error: Matrix indexed with unrecognised type!");
    break;
  }
}
/*   the index of MatrixRow can be a number or table
 *     the indices in the index must be in range of 1 to Nj at the moment this call is evaluated.
 *     the indices stored in Matrix row must be in range of 1 to Ni at the moment his call is evaluated.
 *     in particularly the matrix stored in __Matrix might have changed between the call to Matrix index and MatrixRow index
 *
 *     returned is a table or a table of tables with meta method matrix set
 */
static int LuaMatrixRowIndex(lua_State *L)
{
  int i,j,k,l;
  CompactMatrixType* A;
  lua_getfield(L, 1, "__Matrix");
  A = (CompactMatrixType*)luaL_checkudata(L, -1, "Matrix_Type");
  lua_pop(L, 1);
  lua_getfield(L, 1, "__Row");
  switch (lua_type(L,-1))
  {
    case LUA_TNUMBER:
    { // the rows are indexed with a single number
      i = luaL_checkint(L,-1);
      lua_pop(L, 1);
      switch (lua_type(L,2))
      {
        case LUA_TNUMBER:
        { // the columns and rows are indexed with a single number
          j = luaL_checkint(L,2);
          // check if both i and j are in range
          if(i==0) return luaL_error(L, "Error: Matrix indexed with zero. Indices must be in the range [-Nrow(%d), ..., -1,1,...,Nrow] in the matrix\n",A->Ni);
          if(j==0) return luaL_error(L, "Error: Matrix indexed with zero. Indices must be in the range [-Ncol(%d), ..., -1,1,...,Ncol] in the matrix\n",A->Nj);
          if(Abs(i)>A->Ni) return luaL_error(L, "Error: Matrix index out of range (%d). Indices must be in the range [-Nrow(%d), ..., -1,1,...,Nrow] in the matrix\n",i,A->Ni);
          if(Abs(j)>A->Nj) return luaL_error(L, "Error: Matrix index out of range (%d). Indices must be in the range [-Ncol(%d), ..., -1,1,...,Ncol] in the matrix\n",j,A->Nj);
          if(i<0) i = (int)(A->Ni) +1 + i;
          if(j<0) j = (int)(A->Nj) +1 + j;
          // push a single number to the stack
          if(A->QComplex==0)
          {
            lua_pushnumber(L, A->VRe[j-1 + (i-1) * A->Nj]);
          }
          else
          {
            Complex16 * res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
            res->Re=A->V[j-1 + (i-1) * A->Nj].Re;
            res->Im=A->V[j-1 + (i-1) * A->Nj].Im;
            luaL_setmetatable(L, "Complex_Type");
          }
          return 1;
        }
        break;
        case LUA_TTABLE:
        { // the columns are indexed with a table (at position 2 on the stack) the rows with a number (i)
          int len = (int)lua_objlen(L, 2);
          lua_createtable(L, len, 0);
          for(k=1;k<=len;k++)
          {
            lua_rawgeti(L,2, k);
            j=(int)luaL_checkint(L, -1);
            lua_pop(L, 1);
            if(j ==   0) return luaL_error(L, "Error: Matrix indexed with zero. Indices must be in the range [-Ncol(%d), ..., -1,1,...,Ncol] in the matrix\n",A->Nj);
            if(Abs(j)>A->Nj) return luaL_error(L, "Error: Matrix index out of range (%d). Indices must be in the range [-Ncol(%d), ..., -1,1,...,Ncol] in the matrix\n",j,A->Nj);
            if(j<0) j = (int)(A->Nj) +1 + j;
            if(A->QComplex==0)
            {
              lua_pushnumber(L, A->VRe[j-1 + (i-1) * A->Nj]);
            }
            else
            {
              Complex16 * res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
              res->Re=A->V[j-1 + (i-1) * A->Nj].Re;
              res->Im=A->V[j-1 + (i-1) * A->Nj].Im;
              luaL_setmetatable(L, "Complex_Type");
            }
            lua_rawseti(L, -2, k);
          }
          return 1;
        }
        break;
        default:
          return luaL_error(L, "Error: Matrix indexed with unrecognised type!");
        break;
      }
    }
    break;
    case LUA_TTABLE:
    { // the rows are indexed with a table
      switch (lua_type(L,2))
      {
        case LUA_TNUMBER:
        { // the columns are indexed with a single number the rows with a table (at positon 3==-1 on the stack)
          j = luaL_checkint(L,2);
          if(j ==   0) return luaL_error(L, "Error: Matrix indexed with zero. Indices must be in the range [-Ncol(%d), ..., -1,1,...,Ncol] in the matrix\n",A->Nj);
          if(Abs(j)>A->Nj) return luaL_error(L, "Error: Matrix index out of range (%d). Indices must be in the range [-Ncol(%d), ..., -1,1,...,Ncol] in the matrix\n",j,A->Nj);
          
          int len = (int)lua_objlen(L, 3);
          lua_createtable(L, len, 0);
          for(k=1;k<=len;k++)
          {
            lua_rawgeti(L,3, k);
            i=(int)luaL_checkint(L, -1);
            lua_pop(L, 1);
            if(i ==   0) return luaL_error(L, "Error: Matrix indexed with zero. Indices must be in the range [-Nrow(%d), ..., -1,1,...,Nrow] in the matrix\n",A->Ni);
            if(Abs(i)>A->Ni) return luaL_error(L, "Error: Matrix index out of range (%d). Indices must be in the range [-Nrow(%d), ..., -1,1,...,Nrow] in the matrix\n",i,A->Ni);
            if(i<0) i = (int)(A->Ni) +1 + i;
            if(A->QComplex==0)
            {
              lua_pushnumber(L, A->VRe[j-1 + (i-1) * A->Nj]);
            }
            else
            {
              Complex16 * res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
              res->Re=A->V[j-1 + (i-1) * A->Nj].Re;
              res->Im=A->V[j-1 + (i-1) * A->Nj].Im;
              luaL_setmetatable(L, "Complex_Type");
            }
            lua_rawseti(L, -2, k);
          }
          return 1;
        }
        break;
        case LUA_TTABLE:
        { // the columns (at position 2 on the stack) and the rows (at positon 3==-1 on the stack) are indexed with a table
          int lenj = (int)lua_objlen(L, 2);
          int leni = (int)lua_objlen(L, 3);
          lua_createtable(L, leni, 0);
          lua_getglobal(L, "MatrixMeta");
          lua_setmetatable(L, -2);
          for(k=1;k<=leni;k++)
          {
            lua_rawgeti(L,3, k);
            i=(int)luaL_checkint(L, -1);
            lua_pop(L, 1);
            if(i ==   0) return luaL_error(L, "Error: Matrix indexed with zero. Indices must be in the range [-Nrow(%d), ..., -1,1,...,Nrow] in the matrix\n",A->Ni);
            if(Abs(i)>A->Ni) return luaL_error(L, "Error: Matrix index out of range (%d). Indices must be in the range [-Nrow(%d), ..., -1,1,...,Nrow] in the matrix\n",i,A->Ni);
            if(i<0) i = (int)(A->Ni) +1 + i;
            lua_createtable(L, lenj, 0);
            for(l=1;l<=lenj;l++)
            {
              lua_rawgeti(L,2, l);
              j=(int)luaL_checkint(L, -1);
              lua_pop(L, 1);
              if(j ==   0) return luaL_error(L, "Error: Matrix indexed with zero. Indices must be in the range [-Ncol(%d), ..., -1,1,...,Ncol] in the matrix\n",A->Nj);
              if(Abs(j)>A->Nj) return luaL_error(L, "Error: Matrix index out of range (%d). Indices must be in the range [-Ncol(%d), ..., -1,1,...,Ncol] in the matrix\n",j,A->Nj);
              if(j<0) j = (int)(A->Nj) +1 +j;
              if(A->QComplex==0)
              {
                lua_pushnumber(L, A->VRe[j-1 + (i-1) * A->Nj]);
              }
              else
              {
                Complex16 * res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
                res->Re=A->V[j-1 + (i-1) * A->Nj].Re;
                res->Im=A->V[j-1 + (i-1) * A->Nj].Im;
                luaL_setmetatable(L, "Complex_Type");
              }
              lua_rawseti(L, -2, l);
            }
            lua_rawseti(L, -2, k);
          }
          return 1;
        }
        break;
        default:
          return luaL_error(L, "Error: Matrix indexed with unrecognised type!");
        break;
      }
    }
    break;
    default:
      return luaL_error(L, "Error: Matrix indexed with unrecognised type!");
    break;
  }
}

/*   On stack a MatrixRow object, an index and a new value to assign to the MatrixRow
 *   The index can either be a number or a table
 *   The value to assign must be either a number a table, a table of tables or a matrix
 *   The dimenstion of the object to assign is given by the index and __Row in the MatrixRow
 */
static int LuaMatrixRowNewIndex(lua_State *L)
{
  int i,j,k,l;
  int len,leni,lenj;
  CompactMatrixType* A;
  lua_getfield(L, 1, "__Matrix");
  A = (CompactMatrixType*)luaL_checkudata(L, -1, "Matrix_Type");
  lua_pop(L, 1);
  lua_getfield(L, 1, "__Row");
  switch (lua_type(L,4))
  {
    case LUA_TNUMBER:
    { // the rows are indexed with a single number
      i = luaL_checkint(L,4);
      lua_pop(L, 1);
      switch (lua_type(L,2))
      {
        case LUA_TNUMBER:
        { // the columns and rows are indexed with a single number
          j = luaL_checkint(L,2);
          // check if both i and j are in range
          if(i==0) return luaL_error(L, "Error: Matrix indexed with zero. Indices must be in the range [-Nrow(%d), ..., -1,1,...,Nrow] in the matrix\n",A->Ni);
          if(j==0) return luaL_error(L, "Error: Matrix indexed with zero. Indices must be in the range [-Ncol(%d), ..., -1,1,...,Ncol] in the matrix\n",A->Nj);
          if(Abs(i)>A->Ni) return luaL_error(L, "Error: Matrix index out of range (%d). Indices must be in the range [-Nrow(%d), ..., -1,1,...,Nrow] in the matrix\n",i,A->Ni);
          if(Abs(j)>A->Nj) return luaL_error(L, "Error: Matrix index out of range (%d). Indices must be in the range [-Ncol(%d), ..., -1,1,...,Ncol] in the matrix\n",j,A->Nj);
          if(i<0) i = (int)(A->Ni) +1 + i;
          if(j<0) j = (int)(A->Nj) +1 + j;
          // read a single number from the stack position 3 and set element i,j
          if(lua_isnumber(L, 3))
          {
            if(A->QComplex==0)
            {
              A->VRe[j-1 + (i-1) * A->Nj]=(double)luaL_checknumber(L, 3);
            }
            else
            {
              A->V[j-1 + (i-1) * A->Nj].Re=(double)luaL_checknumber(L, 3);
              A->V[j-1 + (i-1) * A->Nj].Im=0;
            }
          }
          else if (luaL_testudata(L, -1, "Complex_Type") != NULL)
          {
            if(A->QComplex==0)
            {
              MakeCompactRealMatrixComplex(A);
            }
            A->V[j-1 + (i-1) * A->Nj] = *(Complex16 *)luaL_checkudata(L, 3, "Complex_Type");
          }
          else
          {
            return luaL_error(L, "Error trying to set matrix element %d,%d to a value of unknown type\n",i,j);
          }
          return 0;
        }
        break;
        case LUA_TTABLE:
        { // the columns are indexed with a table (at position 2 on the stack) the rows with a number (i)
          lenj = (int)lua_objlen(L, 2);
          if(lua_type(L,3)!=LUA_TTABLE) return luaL_error(L, "Error: Matrix indexed with a number (%d) and a table. You try to set this to be equal to something that is not a table\n",i);
          int  len = (int)lua_objlen(L, 3);
          if(len!=lenj) return luaL_error(L, "Error: Matrix indexed with a number (%d) and a table of length (%d). You try to set this to be equal to a table of length (%d)\n",i,lenj,len);
          for(k=1;k<=lenj;k++)
          {
            lua_rawgeti(L,2, k);
            j=(int)luaL_checkint(L, -1);
            lua_pop(L, 1);
            if(j ==   0) return luaL_error(L, "Error: Matrix indexed with zero. Indices must be in the range [-Ncol(%d), ..., -1,1,...,Ncol] in the matrix\n",A->Nj);
            if(Abs(j)>A->Nj) return luaL_error(L, "Error: Matrix index out of range (%d). Indices must be in the range [-Ncol(%d), ..., -1,1,...,Ncol] in the matrix\n",j,A->Nj);
            if(j<0) j = (int)(A->Nj) +1 + j;
            lua_rawgeti(L,3, k);
            if(lua_isnumber(L, -1))
            {
              if(A->QComplex==0)
              {
                A->VRe[j-1 + (i-1) * A->Nj]=(double)luaL_checknumber(L, -1);
              }
              else
              {
                A->V[j-1 + (i-1) * A->Nj].Re=(double)luaL_checknumber(L, -1);
                A->V[j-1 + (i-1) * A->Nj].Im=0;
              }
            }
            else if (luaL_testudata(L, -1, "Complex_Type") != NULL)
            {
              if(A->QComplex==0)
              {
                MakeCompactRealMatrixComplex(A);
              }
              A->V[j-1 + (i-1) * A->Nj] = *(Complex16 *)luaL_checkudata(L, -1, "Complex_Type");
            }
            else
            {
              return luaL_error(L, "Error trying to set matrix element %d,%d to a value of unknown type\n",i,j);
            }
            lua_pop(L,1);
          }
          return 0;
        }
        break;
        default:
          return luaL_error(L, "Error: Matrix indexed with unrecognised type!");
        break;
      }
    }
    break;
    case LUA_TTABLE:
    { // the rows are indexed with a table
      switch (lua_type(L,2))
      {
        case LUA_TNUMBER:
        { // the columns are indexed with a single number the rows with a table (at positon 4 on the stack)
          j = luaL_checkint(L,2);
          if(j ==   0) return luaL_error(L, "Error: Matrix indexed with zero. Indices must be in the range [-Ncol(%d), ..., -1,1,...,Ncol] in the matrix\n",A->Nj);
          if(Abs(j)>A->Nj) return luaL_error(L, "Error: Matrix index out of range (%d). Indices must be in the range [-Ncol(%d), ..., -1,1,...,Ncol] in the matrix\n",j,A->Nj);
          leni = (int)lua_objlen(L, 4);
          if(lua_type(L,3)!=LUA_TTABLE) return luaL_error(L, "Error: Matrix indexed with a table and a number (%d). You try to set this to be equal to something that is not a table\n",j);
          len = (int)lua_objlen(L, 3);
          if(len!=leni) return luaL_error(L, "Error: Matrix indexed with a table of length (%d) and a number (%d). You try to set this to be equal to a table of length (%d)\n",leni,j,len);
          for(k=1;k<=leni;k++)
          {
            lua_rawgeti(L,4, k);
            i=(int)luaL_checkint(L, -1);
            lua_pop(L, 1);
            if(i ==   0) return luaL_error(L, "Error: Matrix indexed with zero. Indices must be in the range [-Nrow(%d), ..., -1,1,...,Nrow] in the matrix\n",A->Ni);
            if(Abs(i)>A->Ni) return luaL_error(L, "Error: Matrix index out of range (%d). Indices must be in the range [-Nrow(%d), ..., -1,1,...,Nrow] in the matrix\n",i,A->Ni);
            if(i<0) i = (int)(A->Ni) +1 + i;
            lua_rawgeti(L,3, k);
            if(lua_isnumber(L, -1))
            {
              if(A->QComplex==0)
              {
                A->VRe[j-1 + (i-1) * A->Nj]=(double)luaL_checknumber(L, -1);
              }
              else
              {
                A->V[j-1 + (i-1) * A->Nj].Re=(double)luaL_checknumber(L, -1);
                A->V[j-1 + (i-1) * A->Nj].Im=0;
              }
            }
            else if (luaL_testudata(L, -1, "Complex_Type") != NULL)
            {
              if(A->QComplex==0)
              {
                MakeCompactRealMatrixComplex(A);
              }
              A->V[j-1 + (i-1) * A->Nj] = *(Complex16 *)luaL_checkudata(L, -1, "Complex_Type");
            }
            else
            {
              return luaL_error(L, "Error trying to set matrix element %d,%d to a value of unknown type\n",i,j);
            }
            lua_pop(L,1);
          }
          return 0;
        }
        break;
        case LUA_TTABLE:
        { // the columns (at position 2 on the stack) and the rows (at positon 4 on the stack) are indexed with a table
          lenj = (int)lua_objlen(L, 2);
          leni = (int)lua_objlen(L, 4);
          if (luaL_testudata(L, 3, "Matrix_Type") != NULL)
          {
            CompactMatrixType *M = (CompactMatrixType*)luaL_checkudata(L, 3, "Matrix_Type");
            if(M->Ni!=leni) luaL_error(L, "Error trying to set matrix elements for %d rows by %d numbers\n",leni,M->Ni);
            if(M->Nj!=lenj) luaL_error(L, "Error trying to set matrix elements for %d columns by %d numbers\n",lenj,M->Nj);
            if((M->QComplex!=0)&&(A->QComplex==0)) MakeCompactRealMatrixComplex(A);
            for(k=1;k<=leni;k++)
            {
              lua_rawgeti(L,4, k);
              i=(int)luaL_checkint(L, -1);
              lua_pop(L, 1);
              if(i ==   0) return luaL_error(L, "Error: Matrix indexed with zero. Indices must be in the range [-Nrow(%d), ..., -1,1,...,Nrow] in the matrix\n",A->Ni);
              if(Abs(i)>A->Ni) return luaL_error(L, "Error: Matrix index out of range (%d). Indices must be in the range [-Nrow(%d), ..., -1,1,...,Nrow] in the matrix\n",i,A->Ni);
              if(i<0) i = (int)(A->Ni) +1 + i;
              for(l=1;l<=lenj;l++)
              {
                lua_rawgeti(L,2, l);
                j=(int)luaL_checkint(L, -1);
                lua_pop(L, 1);
                if(j ==   0) return luaL_error(L, "Error: Matrix indexed with zero. Indices must be in the range [-Ncol(%d), ..., -1,1,...,Ncol] in the matrix\n",A->Nj);
                if(Abs(j)>A->Nj) return luaL_error(L, "Error: Matrix index out of range (%d). Indices must be in the range [-Ncol(%d), ..., -1,1,...,Ncol] in the matrix\n",j,A->Nj);
                if(j<0) j = (int)(A->Nj) +1 + j;
                if(A->QComplex==0)
                {
                  A->VRe[j-1 + (i-1) * A->Nj]=M->VRe[l-1 + (k-1) * M->Nj];
                }
                else
                {
                  if(M->QComplex==0)
                  {
                    A->V[j-1 + (i-1) * A->Nj].Re=M->VRe[l-1 + (k-1) * M->Nj];
                    A->V[j-1 + (i-1) * A->Nj].Re=0;
                  }
                  else
                  {
                    A->V[j-1 + (i-1) * A->Nj]=M->V[l-1 + (k-1) * M->Nj];
                  }
                }
              }
            }
          }
          if(lua_type(L,3)!=LUA_TTABLE) return luaL_error(L, "Error: Matrix indexed with a table and a table. You try to set this to be equal to something that is not a table of tables or a matrix\n");
          
          int  len = (int)lua_objlen(L, 3);
          if(len!=leni) return luaL_error(L, "Error: Matrix indexed with tables of length %d and %d. You try to set this to be equal to a table of length (%d)\n",leni,lenj,len);
          for(k=1;k<=leni;k++)
          {
            lua_rawgeti(L,4, k);
            i=(int)luaL_checkint(L, -1);
            lua_pop(L, 1);
            if(i ==   0) return luaL_error(L, "Error: Matrix indexed with zero. Indices must be in the range [-Nrow(%d), ..., -1,1,...,Nrow] in the matrix\n",A->Ni);
            if(Abs(i)>A->Ni) return luaL_error(L, "Error: Matrix index out of range (%d). Indices must be in the range [-Nrow(%d), ..., -1,1,...,Nrow] in the matrix\n",i,A->Ni);
            if(i<0) i = (int)(A->Ni) +1 + i;
            lua_rawgeti(L,3, k);
            if(lua_type(L,-1)!=LUA_TTABLE) return luaL_error(L, "Error: Matrix indexed with a table and a table. You try to set this to be equal to something that is not a table of tables or a matrix at row %d\n",k);
            len = (int)lua_objlen(L, -1);
            if(len!=lenj) return luaL_error(L, "Error: Matrix indexed with tables of length %d and %d. You try to set row %d to be equal to a table of length (%d)\n",leni,lenj,k,len);
            for(l=1;l<=lenj;l++)
            {
              lua_rawgeti(L,2, l);
              j=(int)luaL_checkint(L, -1);
              lua_pop(L, 1);
              if(j ==   0) return luaL_error(L, "Error: Matrix indexed with zero. Indices must be in the range [-Ncol(%d), ..., -1,1,...,Ncol] in the matrix\n",A->Nj);
              if(Abs(j)>A->Nj) return luaL_error(L, "Error: Matrix index out of range (%d). Indices must be in the range [-Ncol(%d), ..., -1,1,...,Ncol] in the matrix\n",j,A->Nj);
              if(j<0) j = (int)(A->Nj) +1 + j;
              lua_rawgeti(L,-1, l);
              if(lua_isnumber(L, -1))
              {
                if(A->QComplex==0)
                {
                  A->VRe[j-1 + (i-1) * A->Nj]=(double)luaL_checknumber(L, -1);
                }
                else
                {
                  A->V[j-1 + (i-1) * A->Nj].Re=(double)luaL_checknumber(L, -1);
                  A->V[j-1 + (i-1) * A->Nj].Im=0;
                }
              }
              else if (luaL_testudata(L, -1, "Complex_Type") != NULL)
              {
                if(A->QComplex==0)
                {
                  MakeCompactRealMatrixComplex(A);
                }
                A->V[j-1 + (i-1) * A->Nj] = *(Complex16 *)luaL_checkudata(L, -1, "Complex_Type");
              }
              else
              {
                return luaL_error(L, "Error trying to set matrix element %d,%d to a value of unknown type\n",i,j);
              }
              lua_pop(L, 1);
            }
            lua_pop(L, 1);
          }
          return 0;
        }
        break;
        default:
          return luaL_error(L, "Error: Matrix indexed with unrecognised type!");
        break;
      }
    }
    break;
    default:
      return luaL_error(L, "Error: Matrix indexed with unrecognised type!");
    break;
  }
}


static int LuaMatrixRowLen(lua_State *L)
{
  int len;
  lua_getfield(L, 1, "__Row");
  switch (lua_type(L,-1))
  {
    case LUA_TNUMBER:
    {
      CompactMatrixType* A;
      lua_getfield(L, 1, "__Matrix");
      A = (CompactMatrixType*)luaL_checkudata(L, -1, "Matrix_Type");
      len=(int)(A->Nj);
      lua_pop(L, 1);
    }
    break;
    case LUA_TTABLE:
    {
      len = (int)lua_objlen(L, -1);
    }
    break;
    default:
      return luaL_error(L, "Error: Matrix indexed with unrecognised type!");
    break;
  }
  lua_pop(L, 1);
  lua_pushinteger(L,len);
  return 1;
}

static int LuaMatrixRowToString(lua_State *L)
{
  std::stringstream sstream;
  int colWidth, i, Ni;
  CompactMatrixType* A;
  lua_getfield(L, 1, "__Matrix");
  A = (CompactMatrixType*)luaL_checkudata(L, -1, "Matrix_Type");
  lua_pop(L, 1);
  lua_getfield(L, 1, "__Row");
  switch (lua_type(L,-1))
  {
    case LUA_TNUMBER:
    {
      i = luaL_checkint(L,-1);
      i=i-1;
      lua_pop(L, 1);
      sstream << " { ";
      for (int j = 0; j < A->Nj; ++j)
      {
        if (j != 0) sstream << ", ";
        if (A->QComplex)
        {
          colWidth = 16;
          string tmpStr = numberToStringMaxWidth(A->V[j + i * (A->Nj)].Re, colWidth / 2 - 1);
          if (A->V[j + i * (A->Nj)].Im < 0)
            tmpStr += " -";
          else
            tmpStr += " +";
          tmpStr += numberToStringMaxWidth(Abs(A->V[j + i * (A->Nj)].Im), colWidth - (int)tmpStr.length() - 1);
          sstream << tmpStr << " I";
          sstream << blankFill(colWidth, tmpStr.length());
        }
        else
        {
          colWidth = 8;
          string tmpStr = numberToStringMaxWidth(A->VRe[j + i * (A->Nj)], colWidth - 1);
          sstream << tmpStr;
          sstream << blankFill(colWidth, tmpStr.length());
        }
      }
      sstream << "}" << endl;
      lua_pushstring(L, sstream.str());
      return 1;
    }
    break;
    case LUA_TTABLE:
    {
      Ni = (int)lua_objlen(L, -1);
      for (int ii = 0; ii < Ni; ++ii)
      {
        lua_rawgeti(L,-1, ii+1);
        i=(int)luaL_checkint(L, -1);
        i=i-1;
        lua_pop(L, 1);
        
        if (ii == 0)
          sstream << "{ { ";
        else
          sstream << "  { ";
        for (int j = 0; j < A->Nj; ++j)
        {
          if (j != 0)
            sstream << ", ";
          if (A->QComplex)
          {
            colWidth = 16;
            string tmpStr = numberToStringMaxWidth(A->V[j + i * (A->Nj)].Re, colWidth / 2 - 1);
            if (A->V[j + i * (A->Nj)].Im < 0)
              tmpStr += " -";
            else
              tmpStr += " +";
            tmpStr += numberToStringMaxWidth(Abs(A->V[j + i * (A->Nj)].Im), colWidth - (int)tmpStr.length() - 1);
            sstream << tmpStr << " I";
            sstream << blankFill(colWidth, tmpStr.length());
          } else
          {
            colWidth = 8;
            string tmpStr = numberToStringMaxWidth(A->VRe[j + i * (A->Nj)], colWidth - 1);
            sstream << tmpStr;
            sstream << blankFill(colWidth, tmpStr.length());
          }
        }
        if (ii == Ni - 1)
          sstream << "} }" << endl;
        else
          sstream << "} ," << endl;
      }
      lua_pushstring(L, sstream.str());
      return 1;
    }
    break;
    default:
      return luaL_error(L, "Error: Matrix indexed with unrecognised type!");
    break;
  }
}


/*
 *
 */
static int LuaMatrixDestroy(lua_State * L)
{
  CompactMatrixType *A;
  A = (CompactMatrixType*)luaL_checkudata(L, 1, "Matrix_Type");
  A->OnLuaStack=false; //give back the data
  FreeCompactMatrix(A);
  return 0;
}

//Feed it a table of tables of tables of tables of numbers, and it will return a table of tables.
static int LuaMatrixFlatten(lua_State * L)
{
	int len = lua_gettop(L);
	if (len != 1) return luaL_error(L, "MatrixFlatten called with %d arguments. 1 expected, Table of tables of tables of tables of numbers.\n", len);
	luaL_checktype(L, 1, LUA_TTABLE);
	int numOfRowsOuter, numOfColsOuter;
  int i,j;
  bool ReturnAsUserdata = false;
	numOfRowsOuter = (int) lua_objlen(L, 1);
	lua_rawgeti(L, 1, 1);
	luaL_checktype(L, -1, LUA_TTABLE);
	numOfColsOuter = (int) lua_objlen(L, -1);
	lua_pop(L, 1);
	vector<vector<CompactMatrixType*> > allMs(numOfRowsOuter);
	for (i = 0; i < numOfRowsOuter; ++i)
	{
		lua_rawgeti(L, 1, i + 1);
		luaL_checktype(L, -1, LUA_TTABLE);
		vector<CompactMatrixType*> col (numOfColsOuter);
		allMs[i] = col;
		for (j = 0; j < numOfColsOuter; ++j)
		{
			lua_rawgeti(L, -1, j + 1);
			if (lua_isnumber(L, -1))
			{
				if (luaL_checknumber(L, -1) == 0)
        {
          allMs[i][j] = (CompactMatrixType*)malloc(sizeof(CompactMatrixType));
					*allMs[i][j] = CompactnullMatrix();
        }
				else
        {
					return luaL_error(L, "MatrixFlatten called with entry %d. Matrix or 0 expected.\n", luaL_checknumber(L, -1));
        }
			}
      else
			{
				lua_check(L, -1, allMs[i][j]);
        if(allMs[i][j]->OnLuaStack==true) ReturnAsUserdata=true;
			}
			lua_pop(L, 1);
		}
		lua_pop(L, 1);
	}
  if (ReturnAsUserdata==false)
  {
    CompactMatrixType fullM;
    fullM = MatrixFlatten(allMs, true);
    LuaPushTable_CompactMatrixType(L, &fullM);
    FreeCompactMatrix(&fullM);
  }
  else
  {
    CompactMatrixType *fullM;
    fullM = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
    luaL_setmetatable(L, "Matrix_Type");
    *fullM = MatrixFlatten(allMs, true);
    fullM->OnLuaStack=true;
  }
  for (i = 0; i < numOfRowsOuter; ++i)
  {
    for (j = 0; j < numOfColsOuter; ++j)
    {
      if(allMs[i][j]->OnLuaStack==false)
      {
        FreeCompactMatrix(allMs[i][j]);
        free(allMs[i][j]);
      }
    }
  }
  return 1;
}

static int LuaMatrixAddSpin(lua_State * L)
{
	CompactMatrixType *M;
	uint32_t len = (uint32_t) lua_gettop(L);
	if (len != 1)
	{
		return luaL_error(L, "MatrixAddSpin called with %d arguments. 1 expected, a matrix.\n", len);
	}
  M = luaL_checkCompactMatrix(L, 1);
  if (luaL_testudata(L, 1, "Matrix_Type") == NULL)
  {
    CompactMatrixType fullM;
    fullM = MatrixEnlarge(*M, 2);
    LuaPushTable_CompactMatrixType(L, &fullM);
    FreeCompactMatrix(&fullM);
  }
  else
  {
    CompactMatrixType *fullM;
    fullM = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
    luaL_setmetatable(L, "Matrix_Type");
    *fullM = MatrixEnlarge(*M, 2);
    fullM->OnLuaStack=true;
  }
  FreeCompactMatrix(M);
	return 1;
}

static int LuaMatrixToOperator(lua_State * L)
{
  CompactMatrixType M;
  uint len = (uint) lua_gettop(L);
  if (len != 1)
    return luaL_error(L, "MatrixToOperator called with %d arguments, only 1 expected.", len);
  if (!luaL_isType(L, 1, M))
    return luaL_error(L, "MatrixToOperator called with an argument that is not a matrix.");
  lua_check(L, 1, M);
  lua_push(L, CompactMatrixToOperator(M));
  FreeCompactMatrix(&M);
  return 1;
}

static int LuaMatrixEnlarge(lua_State * L)
{
	CompactMatrixType *M;
	int factor;
	uint32_t len = (uint32_t) lua_gettop(L);
	if (len != 2)
	{
		return luaL_error(L, "MatrixEnlarge called with %d arguments. 2 expected, a matrix and the factor by which it is to be expanded.\n", len);
	}
	M = luaL_checkCompactMatrix(L, 1);
  factor = (int)luaL_checkinteger(L, 2);
  if (luaL_testudata(L, 1, "Matrix_Type") == NULL)
  {
    CompactMatrixType fullM;
    fullM = MatrixEnlarge(*M, factor);
    LuaPushTable_CompactMatrixType(L, &fullM);
    FreeCompactMatrix(&fullM);
  }
  else
  {
    CompactMatrixType *fullM;
    fullM = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
    luaL_setmetatable(L, "Matrix_Type");
    *fullM = MatrixEnlarge(*M, factor);
    fullM->OnLuaStack=true;
  }
  return 1;
}

static int LuaMatrixExpand(lua_State * L)
{
	CompactMatrixType *M;
	int newSize;
  Complex16 newEntries={0,0};
	uint32_t len = (uint32_t) lua_gettop(L);
	if (len < 2 || len > 3)
	{
		return luaL_error(L, "MatrixExpand called with %d arguments. 2 or 3 expected, a matrix, the new matrix size, and possibly a value for the entries on the new diagonal.\n", len);
	}
	M = luaL_checkCompactMatrix(L, 1);
	newSize = (int) luaL_checkinteger(L, 2);
	if (len == 3)
	{
		if (luaL_testudata(L, 3, "Complex_Type") != NULL)
		{
			newEntries = *(Complex16 *) luaL_checkudata(L, 3, "Complex_Type");
		}
    else
		{
			newEntries.Re = luaL_checknumber(L, 3);
			newEntries.Im = 0;
		}
	}
  if (luaL_testudata(L, 1, "Matrix_Type") == NULL)
  {
    CompactMatrixType fullM;
    fullM = MatrixExpand(*M, newSize, newEntries);
    LuaPushTable_CompactMatrixType(L, &fullM);
    FreeCompactMatrix(&fullM);
  }
  else
  {
    CompactMatrixType *fullM;
    fullM = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
    luaL_setmetatable(L, "Matrix_Type");
    *fullM = MatrixExpand(*M, newSize, newEntries);
    fullM->OnLuaStack=true;
  }
	return 1;
}

static int LuaSubmatrix(lua_State * L)
{
	CompactMatrixType *M, *Msub;
	uint32_t len = (uint32_t) lua_gettop(L);
	if (len != 2 && len != 3 && len != 5)
	{
		return luaL_error(L,
				"Submatrix called with %d arguments. Expected either 2 (matrix and size of submatrix, or alternatively matrix and list of indices), 3 (matrix, number of rows and number of columns of submatrix, or alternatively matrix, table with first and last index of row, table with first and last index of column) or 5 (matrix, first row index, number of rows, first column index, number of columns).\n",
				len);
	}
  // read matrix from stack
	M = luaL_checkCompactMatrix(L, 1);
  // push sub matrix to stack as a usertype or allocate memory
  if(M->OnLuaStack==true)
  {
    Msub = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
    luaL_setmetatable(L, "Matrix_Type");
    Msub->OnLuaStack=true;
  }
  else
  {
    Msub = (CompactMatrixType*)malloc(sizeof(CompactMatrixType));
    Msub->OnLuaStack=false;
  }
	if (len == 2)
	{
		int size=0;
		vector < uint32_t > inds;
		if (luaL_isType(L, 2, size))
		{
			size = (int) luaL_checkinteger(L, 2);
			*Msub = Submatrix(*M, size);
		}
		else if (luaL_isType(L, 2, inds))
		{
			lua_check(L, 2, inds);
      for (auto &val : inds)
      {
        if(val==0)
        {
          Msub->QComplex=0;
          Msub->Ni=0;
          Msub->Nj=0;
          InitCompactMatrix(Msub);
          return luaL_error(L,"Matrix elements are indexed starting at 1. Index 0 found\n");
        }
        else if(val>0)
        {
          val -= 1;
        }
        else
        {
          if(M->Ni!=M->Nj)
          {
            Msub->QComplex=0;
            Msub->Ni=0;
            Msub->Nj=0;
            InitCompactMatrix(Msub);
            return luaL_error(L,"Negative index found for submatrix. Counting from the end can only be used for square matrices\n");
          }
          val = M->Ni + val;
          if(val<1)
          {
            Msub->QComplex=0;
            Msub->Ni=0;
            Msub->Nj=0;
            InitCompactMatrix(Msub);
            return luaL_error(L,"Negative index found for submatrix with absolute value larger than dimension of matrix\n");
          }
        }
      }
			*Msub = Submatrix(*M, inds);
		}
		else
		{
      Msub->QComplex=0;
      Msub->Ni=0;
      Msub->Nj=0;
      InitCompactMatrix(Msub);
			return luaL_error(L, "Error: Submatrix called with unknown type for second argument!");
		}
	}
  else if (len == 3)
	{
		if (lua_istable(L, 2) && lua_istable(L, 3))
		{
			if (lua_objlen(L, 2) != 2 || lua_objlen(L, 3) != 2)
      {
        Msub->QComplex=0;
        Msub->Ni=0;
        Msub->Nj=0;
        InitCompactMatrix(Msub);
        return luaL_error(L, "Error, Submatrix expected input (M, {r0, rn}, {c0, cn}).");
      }
			lua_rawgeti(L, 2, 1);
			int rowStart = (int) luaL_checkinteger(L, -1);
			--rowStart;
			lua_pop(L, 1);
			lua_rawgeti(L, 2, 2);
			int rows = (int) luaL_checkinteger(L, -1) - rowStart;
			lua_pop(L, 1);
			lua_rawgeti(L, 3, 1);
			int colStart = (int) luaL_checkinteger(L, -1);
			--colStart;
			lua_pop(L, 1);
			lua_rawgeti(L, 3, 2);
			int cols = (int) luaL_checkinteger(L, -1) - colStart;
			lua_pop(L, 1);
			*Msub = Submatrix(*M, rowStart, rows, colStart, cols);
		} else
		{
			int rows = (int) luaL_checkinteger(L, 2);
			int cols = (int) luaL_checkinteger(L, 3);
			*Msub = Submatrix(*M, rows, cols);
		}
	}
  else if (len == 5)
	{
		int rowStart = (int) luaL_checkinteger(L, 2);
		--rowStart;
		int rows = (int) luaL_checkinteger(L, 3);
		int colStart = (int) luaL_checkinteger(L, 4);
		--colStart;
		int cols = (int) luaL_checkinteger(L, 5);
		*Msub = Submatrix(*M, rowStart, rows, colStart, cols);
	}
	if (isNullMatrix(*Msub))
  {
    Msub->QComplex=0;
    Msub->Ni=0;
    Msub->Nj=0;
    InitCompactMatrix(Msub);
    return luaL_error(L, "Submatrix failed.");
  }
  FreeCompactMatrix(M);
  if(Msub->OnLuaStack==false)
  {
    LuaPushTable_CompactMatrixType(L, Msub);
    FreeCompactMatrix(Msub);
    free(Msub);
  }
  return 1;
}

static int LuaMatrixEq(lua_State * L)
{
	CompactMatrixType *l, *r;
	uint32_t elem,i,j;
	// read matrix from stack
	l = luaL_checkCompactMatrix(L, 1);
	r = luaL_checkCompactMatrix(L, 2);
	if ((l->Ni != r->Ni) || (l->Nj != r->Nj) || (l->QComplex != r->QComplex))
	{
		lua_push(L, false);
    FreeCompactMatrixIncludingPointer(l);
    FreeCompactMatrixIncludingPointer(r);
		return 1;
	}
	if (l->QComplex)
	{
		for (i = 0; i < l->Ni; ++i)
		{
			for (j = 0; j < l->Nj; ++j)
			{
				elem = i * l->Nj + j;
				if ((l->V[elem].Re != r->V[elem].Re) || (l->V[elem].Im != r->V[elem].Im))
				{
					lua_push(L, false);
          FreeCompactMatrixIncludingPointer(l);
          FreeCompactMatrixIncludingPointer(r);
					return 1;
				}
			}
		}
	}
  else
	{
		for (i = 0; i < l->Ni; ++i)
		{
			for (j = 0; j < l->Nj; ++j)
			{
				elem = i * l->Nj + j;
				if (l->VRe[elem] != r->VRe[elem])
				{
					lua_push(L, false);
          FreeCompactMatrixIncludingPointer(l);
          FreeCompactMatrixIncludingPointer(r);
					return 1;
				}
			}
		}
	}
	lua_push(L, true);
  FreeCompactMatrixIncludingPointer(l);
  FreeCompactMatrixIncludingPointer(r);
	return 1;
}

static int LuaMatrixLen(lua_State * L)
{
	CompactMatrixType M;
	lua_check(L, 1, M);
	lua_push(L, M.Ni);
	return 1;
}

static int LauMatrixNew(lua_State *L)
{
  bool Userdata=true;
  uint32_t len = (uint32_t)lua_gettop(L);
  LuaOptionsType LuaFunctionOptions[]={
    {"Userdata","bool",&(Userdata)},
    {NULL, NULL, NULL}};
  if (len==2)
  {
    LuaCheckOptions(L, -1, LuaFunctionOptions);
    len--;
  }
  if(len!=1)
  {
    luaL_error(L, "Error NewMatrix expect as input a table and or a possible set of options.");
  }
  CompactMatrixType *A,*B;
  if(lua_istable(L, 1)==0)
  {
    luaL_error(L, "Error NewMatrix expect as input a table");
  }
  A = (CompactMatrixType *) luaL_checkCompactMatrix(L, 1);
  if(Userdata==true)
  {
    if(A->OnLuaStack==true)
    {
      B = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
      luaL_setmetatable(L, "Matrix_Type");
      CompactMatrixCopy(*A, B);
    }
    else
    {
      B = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
      luaL_setmetatable(L, "Matrix_Type");
      A->OnLuaStack=true;
      *B = *A;
    }
  }
  else
  {
    LuaPushTable_CompactMatrixType(L, A);
    FreeCompactMatrix(A);
  }
  return 1;
}


static const struct luaL_Reg Matrix_func[] = {
	{"AddSpin", LuaMatrixAddSpin},
  {"Conjugate", LuaMatrixConjugate},
  {"ConjugateTranspose", LuaMatrixConjugateTranspose},
  {"Diagonal", LuaMatrixDiagonal},
	{"Enlarge", LuaMatrixEnlarge},
  {"Exp", LuaMatrixExp},
	{"Expand", LuaMatrixExpand},
	{"Flatten",LuaMatrixFlatten},
  {"Identity", LuaMatrixIdentity},
  {"Inverse", LuaMatrixInverse},
  {"New", LauMatrixNew},
  {"Random", LuaMatrixRandom},
  {"RemoveSpin", LuaRemoveSpin},
  {"Rotate", LuaMatrixRotate},
  {"Sqrt", LuaMatrixSqrt},
	{"Sub",LuaSubmatrix},
  {"ToOperator", LuaMatrixToOperator},
  {"ToTable", LuaMatrixToTable},
  {"ToUserdata", LuaMatrixToUserdata},
  {"Transpose", LuaMatrixTranspose},
  {"Zero", LuaMatrixZero},
  {"Eigensystem", LuaMatrixEigensystem},
  {"Tridiagonalize", LuaMatrixTridiagonalize},
  {NULL, NULL}
};

static const struct luaL_Reg Matrix_meta[] = {
  {"__unm", LuaMatrixUnm},
  {"__add", LuaMatrixAdd},
  {"__sub", LuaMatrixSub},
  {"__mul", LuaMatrixMul},
  {"__div", LuaMatrixDiv},
  {"__pow", LuaMatrixPow},
	{"__eq", LuaMatrixEq},
	{"__len", LuaMatrixLen},
  {"__tostring", LuaMatrixToString},
  {"__index", LuaMatrixIndex},
  {"__newindex", LuaMatrixNewIndex},
  {"__gc", LuaMatrixDestroy},
  {NULL, NULL}
};

static const struct luaL_Reg Matrix_Row_meta[] = {
  {"__len", LuaMatrixRowLen},
  {"__tostring", LuaMatrixRowToString},
  {"__index", LuaMatrixRowIndex},
  {"__newindex", LuaMatrixRowNewIndex},
  {NULL, NULL}
};


//TODO: Matrices stored as tables do not behave the same way as matrices stored as userdata when indexed.........

static const struct luaL_Reg MatrixTable_meta[] = {
  {"__unm", LuaMatrixUnm},
  {"__add", LuaMatrixAdd},
  {"__sub", LuaMatrixSub},
  {"__mul", LuaMatrixMul},
  {"__div", LuaMatrixDiv},
  {"__pow", LuaMatrixPow},
	{"__eq", LuaMatrixEq},
  {NULL, NULL}
};
