/*   calculates < psi | psi > or < psi | O | psi >
 
 */
static int LuaBraKet(lua_State * L)
{
  uint32_t len = (uint32_t)lua_gettop(L);
  if (len==2)
  {
    if(lua_type(L,1)==LUA_TTABLE)
    {
      if(lua_type(L,2)==LUA_TTABLE)
      { // table times table
        if(lua_rawequal(L, 1, 2))
        { // same table left and right
          WaveFunctionType *psiA;
          CompactMatrixType res;
          uint32_t i,j;
          psiA = luaL_checkarray_WaveFunctionType(L, 1, &(res.Ni));
          res.Nj=res.Ni;
          res.QComplex=1;
          InitCompactMatrix(&res);
          for(i=0;i<res.Ni;i++)
          {
            WaveFunctionCalculateNormSqr(psiA[i], &(res.V[i+i*res.Nj].Re));
            res.V[i+i*res.Nj].Im=0;
            for(j=i+1;j<res.Nj;j++)
            {
              WaveFunctionDotProduct(psiA[i], psiA[j], &(res.V[j+i*res.Nj].Re), &(res.V[j+i*res.Nj].Im));
              res.V[i+j*res.Nj].Re = res.V[j+i*res.Nj].Re;
              res.V[i+j*res.Nj].Im =-res.V[j+i*res.Nj].Im;
            }
          }
          LuaPushTable_CompactMatrixType(L, &res);
          FreeCompactMatrix(&res);
          free(psiA);
          return 1;
        }
        else
        {
          WaveFunctionType *psiA, *psiB;
          uint32_t i,j;
          CompactMatrixType res;
          psiA = luaL_checkarray_WaveFunctionType(L, 1, &(res.Ni));
          psiB = luaL_checkarray_WaveFunctionType(L, 2, &(res.Nj));
          res.QComplex=1;
          InitCompactMatrix(&res);
          for(i=0;i<res.Ni;i++)
          {
            for(j=0;j<res.Nj;j++)
            {
              if(QWaveFunctionSameBasis(&(psiA[i]), &(psiB[j])))
              {
                WaveFunctionDotProductSameBasis(psiA[i], psiB[j], &(res.V[j+i*res.Nj].Re), &(res.V[j+i*res.Nj].Im));
              }
              else
              {
                WaveFunctionDotProduct(psiA[i], psiB[j], &(res.V[j+i*res.Nj].Re), &(res.V[j+i*res.Nj].Im));
              }
            }
          }
          LuaPushTable_CompactMatrixType(L, &res);
          FreeCompactMatrix(&res);
          free(psiA);
          free(psiB);
          return 1;
        }
      }
      else
      { // table times function
        WaveFunctionType *psiA, *psiB;
        uint32_t i;
        CompactVectorType res;
        psiA = luaL_checkarray_WaveFunctionType(L, 1, &(res.Ni));
        psiB = (WaveFunctionType*) luaL_checkudata(L, 2, "Wavefunction_Type");
        res.QComplex=1;
        res.V= (Complex16*)malloc(res.Ni*sizeof(Complex16));
        for(i=0;i<res.Ni;i++)
        {
          if(QWaveFunctionSame(&(psiA[i]), psiB))
          {
            WaveFunctionCalculateNorm(*psiA, &(res.V[i].Re));
            res.V[i].Im=0;
          }
          else if(QWaveFunctionSameBasis(&(psiA[i]), psiB))
          {
            WaveFunctionDotProductSameBasis(psiA[i], *psiB, &(res.V[i].Re), &(res.V[i].Im));
          }
          else
          {
            WaveFunctionDotProduct(psiA[i], *psiB, &(res.V[i].Re), &(res.V[i].Im));
          }
        }
        LuaPushArray_Complex16(L, res.V, res.Ni);
        free(res.V);
        free(psiA);
        return 1;
      }
    }
    else
    {
      if(lua_type(L,2)==LUA_TTABLE)
      { // function times table
        WaveFunctionType *psiA, *psiB;
        uint32_t i;
        CompactVectorType res;
        psiA = (WaveFunctionType*) luaL_checkudata(L, 1, "Wavefunction_Type");
        psiB = luaL_checkarray_WaveFunctionType(L, 2, &(res.Ni));
        res.QComplex=1;
        res.V= (Complex16*)malloc(res.Ni*sizeof(Complex16));
        for(i=0;i<res.Ni;i++)
        {
          if(QWaveFunctionSame(psiA, &(psiB[i])))
          {
            WaveFunctionCalculateNorm(*psiA, &(res.V[i].Re));
            res.V[i].Im=0;
          }
          else if(QWaveFunctionSameBasis(psiA, &(psiB[i])))
          {
            WaveFunctionDotProductSameBasis(*psiA, psiB[i], &(res.V[i].Re), &(res.V[i].Im));
          }
          else
          {
            WaveFunctionDotProduct(*psiA, psiB[i], &(res.V[i].Re), &(res.V[i].Im));
          }
        }
        LuaPushArray_Complex16(L, res.V, res.Ni);
        free(res.V);
        free(psiB);
        return 1;
      }
      else
      { // function times function
        if(lua_rawequal(L, 1, 2))
        {
          WaveFunctionType *psiA;
          double res;
          psiA = (WaveFunctionType*) luaL_checkudata(L, 1, "Wavefunction_Type");
          WaveFunctionCalculateNormSqr(*psiA, &res);
          lua_pushnumber(L, res);
          return 1;
        }
        else
        {
          WaveFunctionType *psiA, *psiB;
          Complex16 res;
          psiA = (WaveFunctionType*) luaL_checkudata(L, 1, "Wavefunction_Type");
          psiB = (WaveFunctionType*) luaL_checkudata(L, 2, "Wavefunction_Type");
          if(QWaveFunctionSameBasis(psiA, psiB))
          {
            WaveFunctionDotProductSameBasis(*psiA, *psiB, &(res.Re), &(res.Im));
          }
          else
          {
            WaveFunctionDotProduct(*psiA, *psiB, &(res.Re), &(res.Im));
          }
          if(res.Im==0)
          {
            lua_pushnumber(L, res.Re);
          }
          else
          {
            Complex16 *temp;
            temp = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
            temp->Re = res.Re;
            temp->Im = res.Im;
            luaL_setmetatable(L, "Complex_Type");
          }
          return 1;
        }
      }
    }
  }
  else if (len==3)
  {
    if(lua_type(L,1)==LUA_TTABLE)
    {
      if(lua_type(L,3)==LUA_TTABLE)
      { // table times table
        if(lua_rawequal(L, 1, 3))
        { // same table left and right
          WaveFunctionType *psiA;
          OperatorType *O;
          MatrixType res;
          psiA = luaL_checkarray_WaveFunctionType(L, 1, &(res.Ni));
          O    = (OperatorType*)     luaL_checkudata(L, 2, "Operator_Type");
          res.Nj=res.Ni;
          res.QComplex=1;
          InitMatrix(&res);
          PsiOperatorPsiMatrix(*O, psiA, &res);
          LuaPushTable_MatrixType(L, &res);
          FreeMatrix(&res);
          free(psiA);
          return 1;
        }
        else
        {
          WaveFunctionType *psiA,*psiB;
          OperatorType *O;
          MatrixType res;
          psiA = luaL_checkarray_WaveFunctionType(L, 1, &(res.Ni));
          O    = (OperatorType*)     luaL_checkudata(L, 2, "Operator_Type");
          psiB = luaL_checkarray_WaveFunctionType(L, 3, &(res.Nj));
          res.QComplex=1;
          InitMatrix(&res);
          PsiPrimeOperatorPsiMatrix(*O, psiA, psiB, &res);
          LuaPushTable_MatrixType(L, &res);
          FreeMatrix(&res);
          free(psiA);
          free(psiB);
          return 1;
        }
      }
      else
      { // table times function
        WaveFunctionType *psiA,*psiB;
        OperatorType *O;
        uint32_t i;
        CompactVectorType res;
        psiA = luaL_checkarray_WaveFunctionType(L, 1, &(res.Ni));
        O    = (OperatorType*)     luaL_checkudata(L, 2, "Operator_Type");
        psiB = (WaveFunctionType*) luaL_checkudata(L, 3, "Wavefunction_Type");
        res.QComplex=1;
        res.V= (Complex16*)malloc(res.Ni*sizeof(Complex16));
        for(i=0;i<res.Ni;i++)
        {
          if(QWaveFunctionSame(&(psiA[i]),psiB))
          {
            PsiOperatorPsi(*O, *psiB, &(res.V[i].Re), &(res.V[i].Im));
          }
          else
          {
            PsiPrimeOperatorPsi(*O, psiA[i], *psiB, &(res.V[i].Re), &(res.V[i].Im));
          }
        }
        LuaPushArray_Complex16(L, res.V, res.Ni);
        free(res.V);
        free(psiB);
        return 1;
      }
    }
    else
    {
      if(lua_type(L,3)==LUA_TTABLE)
      { // function times table
        WaveFunctionType *psiA,*psiB;
        OperatorType *O;
        uint32_t i;
        CompactVectorType res;
        psiA = (WaveFunctionType*) luaL_checkudata(L, 1, "Wavefunction_Type");
        O    = (OperatorType*)     luaL_checkudata(L, 2, "Operator_Type");
        psiB = luaL_checkarray_WaveFunctionType(L, 3, &(res.Ni));
        res.QComplex=1;
        res.V= (Complex16*)malloc(res.Ni*sizeof(Complex16));
        for(i=0;i<res.Ni;i++)
        {
          if(QWaveFunctionSame(psiA, &(psiB[i])))
          {
            PsiOperatorPsi(*O, *psiA, &(res.V[i].Re), &(res.V[i].Im));
          }
          else
          {
            PsiPrimeOperatorPsi(*O, *psiA, psiB[i], &(res.V[i].Re), &(res.V[i].Im));
          }
        }
        LuaPushArray_Complex16(L, res.V, res.Ni);
        free(res.V);
        free(psiB);
        return 1;
      }
      else
      { // function times function
        WaveFunctionType *psiA, *psiB;
        OperatorType *O;
        Complex16 res;
        psiA = (WaveFunctionType*) luaL_checkudata(L, 1, "Wavefunction_Type");
        O    = (OperatorType*)     luaL_checkudata(L, 2, "Operator_Type");
        psiB = (WaveFunctionType*) luaL_checkudata(L, 3, "Wavefunction_Type");
        if(QWaveFunctionSame(psiA, psiB))
        {
          PsiOperatorPsi(*O, *psiA, &(res.Re), &(res.Im));
        }
        else
        {
          PsiPrimeOperatorPsi(*O, *psiA, *psiB, &(res.Re), &(res.Im));
        }
        if(Abs(res.Im)<10*(DBL_EPSILON*Abs(res.Re)+DBL_MIN))
        {
          lua_pushnumber(L, res.Re);
        }
        else
        {
          Complex16 *temp;
          temp = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
          temp->Re = res.Re;
          temp->Im = res.Im;
          luaL_setmetatable(L, "Complex_Type");
        }
        return 1;
      }
    }
  }
  return luaL_error(L, "Error in BraKet, called with %d arguments, 2 or 3 expected\n",len);
}
