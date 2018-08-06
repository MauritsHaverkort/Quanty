/*
 *
 */
static int LuaMatrixTranspose(lua_State * L)
{
  CompactMatrixType *A;
  A = (CompactMatrixType *) luaL_checkCompactMatrix(L, 1);
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
  A = (CompactMatrixType *) luaL_checkCompactMatrix(L, 1);
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
  A = (CompactMatrixType *) luaL_checkCompactMatrix(L, 1);
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
  A = (CompactMatrixType *) luaL_checkCompactMatrix(L, 1);
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
  A = (CompactMatrixType *) luaL_checkCompactMatrix(L, 1);
  if(A->Ni!=A->Nj) return luaL_error(L,"Matrix sqrt only works on square matrices\n");
  if(QCompactMatrixHermitian(*A, DBL_EPSILON)==0)
  {
    //TODO: matrix sqrt on non hermitian matrices
    return luaL_error(L,"Matrix sqrt only works on Hermitian matrices\n");
    //remove the line above
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
  H=luaL_checkCompactMatrix(L, 1);
  if(H->Ni!=H->Nj) return luaL_error(L,"Matrix Eigensystem only works on square matrices\n");
  if(QCompactMatrixHermitian(*H, 10*DBL_EPSILON)==0)
  {
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
      funR = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
      luaL_setmetatable(L, "Matrix_Type");
      NonHermitianCompactMatrixEigensystem(H,funL, funR,val);
      LuaPushArray_Complex16(L, val, funL->Ni);
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
    }
    return 2;
  }
}

/*
 *
 */
static int LuaMatrixExp(lua_State * L)
{
  CompactMatrixType *H;
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
  bool Userdata=false;
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
    return luaL_error(L,"Matrix Random called with %d arguments, 2 or 3 arguments expected\n");
  }
  luaL_checktype(L, 1, LUA_TTABLE);
  if(lua_objlen(L, 1)!=2) luaL_error(L,"Matrix Random called with table of length %d at posittion 1, expected {min,max}\n",lua_objlen(L, 1));
  luaL_checktype(L, 2, LUA_TTABLE);
  if(lua_objlen(L, 2)!=2) luaL_error(L,"Matrix Random called with table of length %d at posittion 2, expected {Ni,Nj}\n",lua_objlen(L, 2));
  if(Userdata==true)
  {
    res = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
    luaL_setmetatable(L, "Matrix_Type");
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
  bool Userdata=false;
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
  bool Userdata=false;
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
  }
  else
  {
    res = (CompactMatrixType*)malloc(sizeof(CompactMatrixType));
  }
  res->QComplex=0;
  res->Ni = (uint32_t)luaL_checkinteger(L, 1);
  res->Nj = res->Ni;
  InitCompactIdentityMatrix(res);
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
  bool Userdata=false;
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
  if (luaL_testudata(L, 1, "Matrix_Type") == NULL)
  {
    FreeCompactMatrix(A);
    free(A);
  }
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
      free(A);
    }
  }
  else if (luaL_testudata(L, 2, "Matrix_Type") == NULL)
  { // B is read from a table
    CompactMatrixAdd(B, A);
    C = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
    luaL_setmetatable(L, "Matrix_Type");
    *C=*B;
    free(B);
  }
  else
  {
    C = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
    luaL_setmetatable(L, "Matrix_Type");
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
    free(B);
  }
  else
  {
    C = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
    luaL_setmetatable(L, "Matrix_Type");
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
    if(CompactMatrixPow(A,B,(int)n)) return luaL_error(L,"Matrix power failed\n");
  }
  return 1;
}


/*
 *
 */
static int LuaMatrixToString(lua_State * L)
{
  lua_pushfstring(L, "Matrix in userdata format use ToTable() in order to get a table form");
  return 1;
}


/*
 *
 */
static int LuaMatrixIndex(lua_State * L)
{
  uint32_t i,j;
  double n;
  CompactMatrixType *A, *B;
  Complex16 * res;
  A = (CompactMatrixType *) luaL_checkCompactMatrix(L, 1);
  n =  luaL_checknumber(L, 2);
  i = (uint32_t) n;
  if(Abs(n-(double)i)>(1+Abs(n))*DBL_EPSILON*10) return luaL_error(L, "Index of an Array should be an integer\n");
  if(i==0) return luaL_error(L, "Index of an Array should be positive (non-zero)\n");
  if(A->Ni==1)
  {
    if(i>A->Nj) luaL_error(L, "Index %i larger then size of Array %i\n",i,A->Nj);
    if(A->QComplex==0)
    {
      lua_pushnumber(L, A->VRe[i-1]);
    }
    else
    {
      res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
      res->Re=A->V[i-1].Re;
      res->Im=A->V[i-1].Im;
      luaL_setmetatable(L, "Complex_Type");
    }
  }
  else
  {
    if(i>A->Ni) luaL_error(L, "Index %i larger then size of Array %i\n",i,A->Ni);
    B = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
    luaL_setmetatable(L, "Matrix_Type");
    B->Ni=1;
    B->Nj=A->Nj;
    B->QComplex=A->QComplex;
    InitCompactMatrix(B);
    if(A->QComplex==0)
    {
      for(j=0;j<A->Nj;j++)
      {
        B->VRe[j]=A->VRe[(i-1)*A->Nj+j];
      }
    }
    else
    {
      for(j=0;j<A->Nj;j++)
      {
        B->V[j]=A->V[(i-1)*A->Nj+j];
      }
    }
  }
  return 1;
}

/*
 *
 */
static int LuaMatrixDestroy(lua_State * L)
{
  CompactMatrixType *A;
  A = (CompactMatrixType*)luaL_checkudata(L, 1, "Matrix_Type");
  FreeCompactMatrix(A);
  return 0;
}

//Feed it a table of tables of tables of tables of numbers, and it will return a table of tables.
static int LuaMatrixFlatten(lua_State * L)
{
	int len = lua_gettop(L);
	if (len != 1)
		return luaL_error(L, "MatrixFlatten called with %d arguments. 1 expected, Table of tables of tables of tables of numbers.\n", len);
	luaL_checktype(L, 1, LUA_TTABLE);
	int numOfRowsOuter, numOfColsOuter;
	numOfRowsOuter = (int) lua_objlen(L, 1);
	lua_rawgeti(L, 1, 1);
	luaL_checktype(L, -1, LUA_TTABLE);
	numOfColsOuter = (int) lua_objlen(L, -1);
	lua_pop(L, 1);
	vector<vector<MatrixType> > allMs(numOfRowsOuter);
	for (int i = 0; i < numOfRowsOuter; ++i)
	{
		lua_rawgeti(L, 1, i + 1);
		luaL_checktype(L, -1, LUA_TTABLE);
		vector<MatrixType> col (numOfColsOuter);
		allMs[i] = col;
		for (int j = 0; j < numOfColsOuter; ++j)
		{
			lua_rawgeti(L, -1, j + 1);
			if (lua_isnumber(L, -1))
			{
				if (luaL_checknumber(L, -1) == 0)
					allMs[i][j] = nullMatrix();
				else
					return luaL_error(L, "MatrixFlatten called with entry %d. Matrix or 0 expected.\n", luaL_checknumber(L, -1));
			} else
			{
				allMs[i][j]  = *luaL_checkMatrix(L, -1);
			}
			lua_pop(L, 1);
		}
		lua_pop(L, 1);
	}
	MatrixType fullM = MatrixFlatten(allMs);
	LuaPushTable_MatrixType(L, &fullM);
	return 1;
}

static int LuaMatrixAddSpin(lua_State * L)
{
	MatrixType M;
	uint32_t len = (uint32_t) lua_gettop(L);
	if (len != 1)
	{
		return luaL_error(L, "MatrixAddSpin called with %d arguments. 1 expected, a matrix.\n", len);
	}
	M = *luaL_checkMatrix(L, 1);
	MatrixType fullM = MatrixExpand(M, 2);
	LuaPushTable_MatrixType(L, &fullM);
	return 1;
}

static int LuaMatrixExpand(lua_State * L)
{
	MatrixType M;
	int factor;
	uint32_t len = (uint32_t) lua_gettop(L);
	if (len != 2)
	{
		return luaL_error(L, "MatrixExpand called with %d arguments. 2 expected, a matrix and the factor by which it is to be expanded.\n", len);
	}
	M = *luaL_checkMatrix(L, 1);
	factor = (int)luaL_checkinteger(L, 2);
	MatrixType fullM = MatrixExpand(M, factor);
	LuaPushTable_MatrixType(L, &fullM);
	return 1;
}


static const struct luaL_Reg Matrix_func[] = {
	{"AddSpin", LuaMatrixAddSpin},
  {"Conjugate", LuaMatrixConjugate},
  {"ConjugateTranspose", LuaMatrixConjugateTranspose},
  {"Diagonal", LuaMatrixDiagonal},
  {"Exp", LuaMatrixExp},
	{"Expand", LuaMatrixExpand},
	{"Flatten",LuaMatrixFlatten},
  {"Identity", LuaMatrixIdentity},
  {"Inverse", LuaMatrixInverse},
  {"Random", LuaMatrixRandom},
  {"Sqrt", LuaMatrixSqrt},
  {"ToTable", LuaMatrixToTable},
  {"ToUserdata", LuaMatrixToUserdata},
  {"Transpose", LuaMatrixTranspose},
  {"Zero", LuaMatrixZero},
  {NULL, NULL}
};

static const struct luaL_Reg Matrix_meta[] = {
  {"__unm", LuaMatrixUnm},
  {"__add", LuaMatrixAdd},
  {"__sub", LuaMatrixSub},
  {"__mul", LuaMatrixMul},
  {"__div", LuaMatrixDiv},
  {"__pow", LuaMatrixPow},
  {"__tostring", LuaMatrixToString},
  {"__index", LuaMatrixIndex},
  {"__gc", LuaMatrixDestroy},
  {NULL, NULL}
};

static const struct luaL_Reg MatrixTable_meta[] = {
  {"__unm", LuaMatrixUnm},
  {"__add", LuaMatrixAdd},
  {"__sub", LuaMatrixSub},
  {"__mul", LuaMatrixMul},
  {"__div", LuaMatrixDiv},
  {"__pow", LuaMatrixPow},
  {NULL, NULL}
};
