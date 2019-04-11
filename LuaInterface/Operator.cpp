/*  define the NewOperator functions
 *  depending on the first argument of NewOperator different types of standard operators are created. Here the implementation
 */

static int LuaOperatorNewCoulomb(lua_State * L)
{
  OperatorType *op;
  uint32_t Ntemp,NShell,l,lA,lB,lC,lD;
  bool jbasis=false;
  uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t len = (uint32_t)lua_gettop(L);
  if(lua_isnumber(L,3)) jbasis=true;
  if(NFermion==0) return luaL_error(L, "Coulomb operator acting on zero fermions is not defined");
  if ((len==4)&&(jbasis==false)) /* creates the coulomb operator in a single shell of angular momentum L */
  {
    uint16_t *LindexUp, *LindexDn;
    double *F;
    LindexDn = luaL_checkarray_unit16_t(L, 2, &NShell);
    LindexUp = luaL_checkarray_unit16_t(L, 3, &Ntemp);
    if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexUp, NShell)==false) return luaL_error(L, "Index up contains a refference outside the range [0,%d]",NFermion-1);
    if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexDn, Ntemp )==false) return luaL_error(L, "Index dn contains a refference outside the range [0,%d]",NFermion-1);
    if(NShell!=Ntemp)
    {
      return luaL_error(L, "Error in LuaOperatorNewCoulomb, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
    }
    l=(NShell-1)/2;
    if(2*l+1!=NShell)
    {
      return luaL_error(L, "Error in LuaOperatorNewCoulomb, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
    }
    F = luaL_checkarray_double(L, 4, &Ntemp);
    if(Ntemp!=l+1)
    {
      return luaL_error(L, "Error in LuaOperatorNewCoulomb, Lenght of F (%d) is not equal to l+1\n",Ntemp);
    }
    op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    CreateOperatorCoulombSingleShell(NFermion, l, LindexUp, LindexDn, F, op);
    luaL_setmetatable(L, "Operator_Type");
    free(LindexUp);
    free(LindexDn);
    free(F);
    return 1;
  }
  else if ((len==7)&&(jbasis==false)) /* creates the coulomb operator between two shells with LA and LB conserving the occupation of the A and B shell */
  {
    uint16_t *LindexUpA, *LindexDnA, *LindexUpB, *LindexDnB;
    double *F,*G;
    LindexDnA = luaL_checkarray_unit16_t(L, 2, &NShell);
    LindexUpA = luaL_checkarray_unit16_t(L, 3, &Ntemp);
    if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexUpA, NShell)==false) return luaL_error(L, "Index up contains a refference outside the range [0,%d]",NFermion-1);
    if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexDnA, Ntemp )==false) return luaL_error(L, "Index dn contains a refference outside the range [0,%d]",NFermion-1);
    if(NShell!=Ntemp)
    {
      return luaL_error(L, "Error in LuaOperatorNewCoulomb, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
    }
    lA=(NShell-1)/2;
    if(2*lA+1!=NShell)
    {
      return luaL_error(L, "Error in LuaOperatorNewCoulomb, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
    }
    LindexDnB = luaL_checkarray_unit16_t(L, 4, &NShell);
    LindexUpB = luaL_checkarray_unit16_t(L, 5, &Ntemp);
    if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexUpB, NShell)==false) return luaL_error(L, "Index up contains a refference outside the range [0,%d]",NFermion-1);
    if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexDnB, Ntemp )==false) return luaL_error(L, "Index dn contains a refference outside the range [0,%d]",NFermion-1);
    if(NShell!=Ntemp)
    {
      return luaL_error(L, "Error in LuaOperatorNewCoulomb, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
    }
    lB=(NShell-1)/2;
    if(2*lB+1!=NShell)
    {
      return luaL_error(L, "Error in LuaOperatorNewCoulomb, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
    }
    F = luaL_checkarray_double(L, 6, &Ntemp);
    if(Ntemp!=(Min(lA,lB)+1))
    {
      return luaL_error(L, "Error in LuaOperatorNewCoulomb, Lenght of F (%d) is not equal to (Min(lA,lB)+1)\n",Ntemp);
    }
    G = luaL_checkarray_double(L, 7, &Ntemp);
    if(Ntemp!=(Min(lA,lB)+1))
    {
      return luaL_error(L, "Error in LuaOperatorNewCoulomb, Lenght of G (%d) is not equal to (Min(lA,lB)+1)\n",Ntemp);
    }
    op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    CreateOperatorCoulombTwoShellsConstantOccupation(NFermion, lA, lB, LindexUpA, LindexDnA, LindexUpB, LindexDnB, F, G, op);
    luaL_setmetatable(L, "Operator_Type");
    free(LindexUpA);
    free(LindexDnA);
    free(LindexUpB);
    free(LindexDnB);
    free(F);
    free(G);
    return 1;
  }
  else if ((len==10)&&(jbasis==false)) /* creates the general coulomb operator between four shels with LA, LB, LC and LD */
  {
    uint16_t *LindexUpA, *LindexDnA, *LindexUpB, *LindexDnB, *LindexUpC, *LindexDnC, *LindexUpD, *LindexDnD;
    double *R;
    LindexDnA = luaL_checkarray_unit16_t(L, 2, &NShell);
    LindexUpA = luaL_checkarray_unit16_t(L, 3, &Ntemp);
    if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexUpA, NShell)==false) return luaL_error(L, "Index up of first shell contains a refference outside the range [0,%d]",NFermion-1);
    if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexDnA, Ntemp )==false) return luaL_error(L, "Index dn of first shell contains a refference outside the range [0,%d]",NFermion-1);
    if(NShell!=Ntemp)
    {
      return luaL_error(L, "Error in LuaOperatorNewCoulomb, For first shell lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
    }
    lA=(NShell-1)/2;
    if(2*lA+1!=NShell)
    {
      return luaL_error(L, "Error in LuaOperatorNewCoulomb, For first shell lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
    }
    LindexDnB = luaL_checkarray_unit16_t(L, 4, &NShell);
    LindexUpB = luaL_checkarray_unit16_t(L, 5, &Ntemp);
    if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexUpB, NShell)==false) return luaL_error(L, "Index up of second shell contains a refference outside the range [0,%d]",NFermion-1);
    if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexDnB, Ntemp )==false) return luaL_error(L, "Index dn of second shell contains a refference outside the range [0,%d]",NFermion-1);
    if(NShell!=Ntemp)
    {
      return luaL_error(L, "Error in LuaOperatorNewCoulomb, For second shell lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
    }
    lB=(NShell-1)/2;
    if(2*lB+1!=NShell)
    {
      return luaL_error(L, "Error in LuaOperatorNewCoulomb, For second shell lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
    }
    LindexDnC = luaL_checkarray_unit16_t(L, 6, &NShell);
    LindexUpC = luaL_checkarray_unit16_t(L, 7, &Ntemp);
    if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexUpC, NShell)==false) return luaL_error(L, "Index up of third shell contains a refference outside the range [0,%d]",NFermion-1);
    if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexDnC, Ntemp )==false) return luaL_error(L, "Index dn of third shell contains a refference outside the range [0,%d]",NFermion-1);
    if(NShell!=Ntemp)
    {
      return luaL_error(L, "Error in LuaOperatorNewCoulomb, For third shell lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
    }
    lC=(NShell-1)/2;
    if(2*lC+1!=NShell)
    {
      return luaL_error(L, "Error in LuaOperatorNewCoulomb, For third shell lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
    }
    LindexDnD = luaL_checkarray_unit16_t(L, 8, &NShell);
    LindexUpD = luaL_checkarray_unit16_t(L, 9, &Ntemp);
    if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexUpD, NShell)==false) return luaL_error(L, "Index up of forth shell contains a refference outside the range [0,%d]",NFermion-1);
    if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexDnD, Ntemp )==false) return luaL_error(L, "Index dn of forth shell contains a refference outside the range [0,%d]",NFermion-1);
    if(NShell!=Ntemp)
    {
      return luaL_error(L, "Error in LuaOperatorNewCoulomb, For forth shell lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
    }
    lD=(NShell-1)/2;
    if(2*lD+1!=NShell)
    {
      return luaL_error(L, "Error in LuaOperatorNewCoulomb, For forth shell lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
    }
    if((lA+lB+lC+lD)&1)
    {
      return luaL_error(L, "Error in LuaOperatorNewCoulomb, lA(%d)+lB(%d)+lC(%d)+lD(%d) should be even.",lA,lB,lC,lD);
    }
    R = luaL_checkarray_double(L, 10, &Ntemp);
    if(Ntemp!=(Min(lA+lC,lB+lD)-Max(Max(lA,lC)-Min(lA,lC),Max(lB,lD)-Min(lB,lD)))/2+1)
    {
      return luaL_error(L, "Error in LuaOperatorNewCoulomb, Lenght of R (%d) is not equal to (Min(lA+lC,lB+lD)-Max(Abs(lA-lC),Abs(lB-lD)))/2+1\n",Ntemp);
    }
    op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    CreateOperatorCoulombFourShells(NFermion, lA, lB, lC, lD, LindexUpA, LindexDnA, LindexUpB, LindexDnB, LindexUpC, LindexDnC, LindexUpD, LindexDnD, R, op);
    luaL_setmetatable(L, "Operator_Type");
    free(LindexUpA);
    free(LindexDnA);
    free(LindexUpB);
    free(LindexDnB);
    free(LindexUpC);
    free(LindexDnC);
    free(LindexUpD);
    free(LindexDnD);
    free(R);
    return 1;
  }
  else if (len==4) /* creates the coulomb operator in a single shell of total angular momentum J */
  {
    // Althought the large and small componten of the wave function have a differnt l, the angular part works out such that they have the same form of the operator
    // The resulting slater integral is
    // Int dr1 dr2 Min[r1,r2]^k / Max[r1,r2]^(k+1)  ( G(r1)^2 G(r2)^2 + 2 G(r1)^2 F(r2)^2 + F(r1)^2 F(r2)^2 )
    // with k 0, 2, 4, ... and k <= 2j
    //
    uint32_t twoj;
    uint16_t *Jindex;
    int kappa;
    double *F;
    Jindex = luaL_checkarray_unit16_t(L, 2, &NShell);
    if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), Jindex, NShell)==false) return luaL_error(L, "Index of Jshell contains a refference outside the range [0,%d]",NFermion-1);
    twoj=(NShell-1);
    if(NShell&1)
    {
      return luaL_error(L, "Error in LuaOperatorNewCoulomb, Lenght of Index for J shell (%d) is not equal to 2j+1 for half integer j\n",NShell);
    }
    kappa = luaL_checknumber(L, 3);
    if(twoj!=2*Abs(kappa)-1)
    {
      return luaL_error(L, "Error in LuaOperatorNewCoulomb, j (%d/2) is not equal to |kappa|-1/2 (kappa = %d)\n",twoj,kappa);
    }
    F = luaL_checkarray_double(L, 4, &Ntemp);
    if(Ntemp!=(twoj+1)/2)
    {
      return luaL_error(L, "Error in LuaOperatorNewCoulomb, Lenght of Flarge (%d) is not equal to (2j+1)/2 (%d+1)/2\n",Ntemp,twoj);
    }
    op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    CreateOperatorCoulombSingleRelativisticShell(NFermion, kappa, Jindex, F, op);
    luaL_setmetatable(L, "Operator_Type");
    free(Jindex);
    free(F);
    return 1;
  }
  else if (len==7) /* creates the coulomb operator between two shells with JA and JB conserving the occupation of the A and B shell */
  {
    // Althought the large and small componten of the wave function have a differnt l, the angular part works out such that they have the same form of the operator
    // The resulting slater integral is
    // F[k] = Int dr1 dr2 Min[r1,r2]^k / Max[r1,r2]^(k+1)  ( G1(r1)^2 G2(r2)^2 + G1(r1)^2 F2(r2)^2 + G2(r1)^2 F1(r2)^2 + F1(r1)^2 F2(r2)^2 )
    // with k 0, 2, 4, ... and k <= 2 Min[j1,j2]
    // G[k] = Int dr1 dr2 Min[r1,r2]^k / Max[r1,r2]^(k+1)  ( G1(r1) G2(r1) G1(r2) G2(r2) + F1(r1) F2(r1) F1(r2) F2(r2)  )
    // with k |j1-j2|, |j1-j2|+2, |j1-j2|+4, ... and k <= j1+j2
    //
    int kappaA, kappaB;
    uint32_t twojA,twojB;
    uint16_t *JindexA,*JindexB;
    double *F, *G;
    JindexA = luaL_checkarray_unit16_t(L, 2, &NShell);
    if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), JindexA, NShell)==false) return luaL_error(L, "Index of first J shell contains a refference outside the range [0,%d]",NFermion-1);
    twojA   = NShell-1;
    if(NShell&1)
    {
      return luaL_error(L, "Error in LuaOperatorNewCoulomb, Lenght of Index for first J shell (%d) is not equal to 2j+1 for half integer j\n",NShell);
    }
    kappaA = luaL_checknumber(L, 3);
    if(twojA!=2*Abs(kappaA)-1)
    {
      return luaL_error(L, "Error in LuaOperatorNewCoulomb, for first shell: j (%d/2) is not equal to |kappa|-1/2 (kappa = %d)\n",twojA,kappaA);
    }
    JindexB = luaL_checkarray_unit16_t(L, 4, &NShell);
    if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), JindexB, NShell)==false) return luaL_error(L, "Index of second J shell contains a refference outside the range [0,%d]",NFermion-1);
    twojB   = NShell-1;
    if(NShell&1)
    {
      return luaL_error(L, "Error in LuaOperatorNewCoulomb, Lenght of Index for second J shell (%d) is not equal to 2j+1 for half integer j\n",NShell);
    }
    kappaB = luaL_checknumber(L, 5);
    if(twojB!=2*Abs(kappaB)-1)
    {
      return luaL_error(L, "Error in LuaOperatorNewCoulomb, for second shell: j (%d/2) is not equal to |kappa|-1/2 (kappa = %d)\n",twojB,kappaB);
    }
    F = luaL_checkarray_double(L, 6, &Ntemp);
    if(Ntemp!=(Min(twojA,twojB)+1)/2)
    {
      return luaL_error(L, "Error in LuaOperatorNewCoulomb, Lenght of F (%d) is not equal to (Min(2 j1 (%d),2 j2 (%d))+1)/2\n",Ntemp,twojA,twojB);
    }
    G = luaL_checkarray_double(L, 7, &Ntemp);
    if(Ntemp!=(Min(twojA,twojB)+1)/2)
    {
      return luaL_error(L, "Error in LuaOperatorNewCoulomb, Lenght of G (%d) is not equal to (Min(2 j1 (%d),2 j2 (%d))+1)/2\n",Ntemp,twojA,twojB);
    }
    op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    CreateOperatorCoulombTwoRelativisticShellsConstantOccupation(NFermion, kappaA, kappaB, JindexA, JindexB, F, G, op);
    luaL_setmetatable(L, "Operator_Type");
    free(JindexA);
    free(JindexB);
    free(F);
    free(G);
    return 1;
  }
  else if (len==10) /* creates the general coulomb operator between four shells with JA, JB, JC and JD */
  {
    // Althought the large and small componten of the wave function have a differnt l, the angular part works out such that they have the same form of the operator
    // The resulting slater integral is
    // R[k] = Int dr1 dr2 Min[r1,r2]^k / Max[r1,r2]^(k+1)  ( G3(r1)G1(r1) + F3(r1)F1(r1) ) ( G4(r2)G2(r2) + F4(r2)F2(r2) )
    // with k Max[|j1-j3|,|j2-j4|], 2+Max[|j1-j3|,|j2-j4|], ... and k <= Min[j1+j3,j2+j4]
    //
    int kappaA, kappaB, kappaC, kappaD;
    uint32_t twojA,twojB,twojC,twojD;
    uint16_t *JindexA,*JindexB,*JindexC,*JindexD;
    double *R;
    JindexA = luaL_checkarray_unit16_t(L, 2, &NShell);
    if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), JindexA, NShell)==false) return luaL_error(L, "Index of first J shell contains a refference outside the range [0,%d]",NFermion-1);
    twojA   = NShell-1;
    if(NShell&1)
    {
      return luaL_error(L, "Error in LuaOperatorNewCoulomb, Lenght of Index for first J shell (%d) is not equal to 2j+1 for half integer j\n",NShell);
    }
    kappaA = luaL_checknumber(L, 3);
    if(twojA!=2*Abs(kappaA)-1)
    {
      return luaL_error(L, "Error in LuaOperatorNewCoulomb, for first shell: j (%d/2) is not equal to |kappa|-1/2 (kappa = %d)\n",twojA,kappaA);
    }
    JindexB = luaL_checkarray_unit16_t(L, 4, &NShell);
    if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), JindexB, NShell)==false) return luaL_error(L, "Index of second J shell contains a refference outside the range [0,%d]",NFermion-1);
    twojB   = NShell-1;
    if(NShell&1)
    {
      return luaL_error(L, "Error in LuaOperatorNewCoulomb, Lenght of Index for second J shell (%d) is not equal to 2j+1 for half integer j\n",NShell);
    }
    kappaB = luaL_checknumber(L, 5);
    if(twojB!=2*Abs(kappaB)-1)
    {
      return luaL_error(L, "Error in LuaOperatorNewCoulomb, for second shell: j (%d/2) is not equal to |kappa|-1/2 (kappa = %d)\n",twojB,kappaB);
    }
    JindexC = luaL_checkarray_unit16_t(L, 6, &NShell);
    if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), JindexC, NShell)==false) return luaL_error(L, "Index of third J shell contains a refference outside the range [0,%d]",NFermion-1);
    twojC   = NShell-1;
    if(NShell&1)
    {
      return luaL_error(L, "Error in LuaOperatorNewCoulomb, Lenght of Index for third J shell (%d) is not equal to 2j+1 for half integer j\n",NShell);
    }
    kappaC = luaL_checknumber(L, 7);
    if(twojC!=2*Abs(kappaC)-1)
    {
      return luaL_error(L, "Error in LuaOperatorNewCoulomb, for third shell: j (%d/2) is not equal to |kappa|-1/2 (kappa = %d)\n",twojC,kappaC);
    }
    JindexD = luaL_checkarray_unit16_t(L, 8, &NShell);
    if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), JindexD, NShell)==false) return luaL_error(L, "Index of forth J shell contains a refference outside the range [0,%d]",NFermion-1);
    twojD   = NShell-1;
    if(NShell&1)
    {
      return luaL_error(L, "Error in LuaOperatorNewCoulomb, Lenght of Index for forth J shell (%d) is not equal to 2j+1 for half integer j\n",NShell);
    }
    kappaD = luaL_checknumber(L, 9);
    if(twojD!=2*Abs(kappaD)-1)
    {
      return luaL_error(L, "Error in LuaOperatorNewCoulomb, for forth shell: j (%d/2) is not equal to |kappa|-1/2 (kappa = %d)\n",twojD,kappaD);
    }
    R = luaL_checkarray_double(L, 10, &Ntemp);
    if(Ntemp!=(1+(Min(twojA+twojC,twojB+twojD)-Max(Max(twojA,twojC)-Min(twojA,twojC),Max(twojB,twojD)-Min(twojB,twojD)))/4))
    {
      return luaL_error(L, "Error in LuaOperatorNewCoulomb (10), Lenght of F (%d) is not equal to (Max(2 j1 (%d) + 2 j3 (%d), 2 j2 (%d) + 2 j4 (%d)) - Min( |2j1-2j3|,|2j2-2j4| ))/2)+1\n",Ntemp,twojA,twojC,twojB,twojD);
    }
    op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    CreateOperatorCoulombFourRelativisticShells(NFermion, kappaA, kappaB, kappaC, kappaD, JindexA, JindexB, JindexC, JindexD, R, op);
    luaL_setmetatable(L, "Operator_Type");
    free(JindexA);
    free(JindexB);
    free(JindexC);
    free(JindexD);
    free(R);
    return 1;
  }
  else
  {
    return luaL_error(L, "Error in LuaOperatorNewCoulomb, expected 4, 7, or 10 arguments for Coulomb repulsion in ls coupling basis or 3, 5, or 6 arguments for Coulomb repulsion in jj coupling basis, not %d",len);
  }
  return 1;
}

static int LuaOperatorNewSx(lua_State * L)
{
  OperatorType *op;
  uint32_t Ntemp;
  uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t NShell;
  uint16_t *LindexUp, *LindexDn;
  if(NFermion==0) return luaL_error(L, "Spin operator acting on zero fermions is not defined");
  LindexUp = luaL_checkarray_unit16_t(L, 2, &NShell);
  LindexDn = luaL_checkarray_unit16_t(L, 3, &Ntemp);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexUp, NShell)==false) return luaL_error(L, "Index up contains a refference outside the range [0,%d]",NFermion-1);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexDn, Ntemp )==false) return luaL_error(L, "Index dn contains a refference outside the range [0,%d]",NFermion-1);
  if(NShell!=Ntemp)
  {
    return luaL_error(L, "Error in LuaOperatorNewSx, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorSx(NFermion, NShell, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "Operator_Type");
  free(LindexUp);
  free(LindexDn);
  return 1;
}

static int LuaOperatorNewSy(lua_State * L)
{
  OperatorType *op;
  uint32_t Ntemp;
  uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t NShell;
  uint16_t *LindexUp, *LindexDn;
  if(NFermion==0) return luaL_error(L, "Spin operator acting on zero fermions is not defined");
  LindexUp = luaL_checkarray_unit16_t(L, 2, &NShell);
  LindexDn = luaL_checkarray_unit16_t(L, 3, &Ntemp);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexUp, NShell)==false) return luaL_error(L, "Index up contains a refference outside the range [0,%d]",NFermion-1);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexDn, Ntemp )==false) return luaL_error(L, "Index dn contains a refference outside the range [0,%d]",NFermion-1);
  if(NShell!=Ntemp)
  {
    return luaL_error(L, "Error in LuaOperatorNewSy, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorSy(NFermion, NShell, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "Operator_Type");
  free(LindexUp);
  free(LindexDn);
  return 1;
}

static int LuaOperatorNewSz(lua_State * L)
{
  OperatorType *op;
  uint32_t Ntemp;
  uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t NShell;
  uint16_t *LindexUp, *LindexDn;
  if(NFermion==0) return luaL_error(L, "Spin operator acting on zero fermions is not defined");
  LindexUp = luaL_checkarray_unit16_t(L, 2, &NShell);
  LindexDn = luaL_checkarray_unit16_t(L, 3, &Ntemp);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexUp, NShell)==false) return luaL_error(L, "Index up contains a refference outside the range [0,%d]",NFermion-1);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexDn, Ntemp )==false) return luaL_error(L, "Index dn contains a refference outside the range [0,%d]",NFermion-1);
  if(NShell!=Ntemp)
  {
    return luaL_error(L, "Error in LuaOperatorNewSz, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorSz(NFermion, NShell, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "Operator_Type");
  free(LindexUp);
  free(LindexDn);
  return 1;
}

static int LuaOperatorNewSsqr(lua_State * L)
{
  OperatorType *op;
  uint32_t Ntemp;
  uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t NShell;
  uint16_t *LindexUp, *LindexDn;
  if(NFermion==0) return luaL_error(L, "Spin operator acting on zero fermions is not defined");
  LindexUp = luaL_checkarray_unit16_t(L, 2, &NShell);
  LindexDn = luaL_checkarray_unit16_t(L, 3, &Ntemp);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexUp, NShell)==false) return luaL_error(L, "Index up contains a refference outside the range [0,%d]",NFermion-1);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexDn, Ntemp )==false) return luaL_error(L, "Index dn contains a refference outside the range [0,%d]",NFermion-1);
  if(NShell!=Ntemp)
  {
    return luaL_error(L, "Error in LuaOperatorNewSsqr, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorSsqr(NFermion, NShell, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "Operator_Type");
  free(LindexUp);
  free(LindexDn);
  return 1;
}

static int LuaOperatorNewSplus(lua_State * L)
{
  OperatorType *op;
  uint32_t Ntemp;
  uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t NShell;
  uint16_t *LindexUp, *LindexDn;
  if(NFermion==0) return luaL_error(L, "Spin operator acting on zero fermions is not defined");
  LindexUp = luaL_checkarray_unit16_t(L, 2, &NShell);
  LindexDn = luaL_checkarray_unit16_t(L, 3, &Ntemp);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexUp, NShell)==false) return luaL_error(L, "Index up contains a refference outside the range [0,%d]",NFermion-1);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexDn, Ntemp )==false) return luaL_error(L, "Index dn contains a refference outside the range [0,%d]",NFermion-1);
  if(NShell!=Ntemp)
  {
    return luaL_error(L, "Error in LuaOperatorNewSplus, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorSplus(NFermion, NShell, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "Operator_Type");
  free(LindexUp);
  free(LindexDn);
  return 1;
}

static int LuaOperatorNewSmin(lua_State * L)
{
  OperatorType *op;
  uint32_t Ntemp;
  uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t NShell;
  uint16_t *LindexUp, *LindexDn;
  if(NFermion==0) return luaL_error(L, "Spin operator acting on zero fermions is not defined");
  LindexUp = luaL_checkarray_unit16_t(L, 2, &NShell);
  LindexDn = luaL_checkarray_unit16_t(L, 3, &Ntemp);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexUp, NShell)==false) return luaL_error(L, "Index up contains a refference outside the range [0,%d]",NFermion-1);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexDn, Ntemp )==false) return luaL_error(L, "Index dn contains a refference outside the range [0,%d]",NFermion-1);
  if(NShell!=Ntemp)
  {
    return luaL_error(L, "Error in LuaOperatorNewSmin, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorSmin(NFermion, NShell, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "Operator_Type");
  free(LindexUp);
  free(LindexDn);
  return 1;
}

static int LuaOperatorNewLx(lua_State * L)
{
  struct
  {
    const char *basis;
    uint32_t ibasis;   //  0  000 SphericalHarmonicY, alternating dn,up from m=-l to m=l
                       //  1  001 TesseralHarmonicZ,  alternating dn,up sin first, than cos
                       //  2  010 KubicHarmonicK,     alternatind dn,up order according to irep a,b,e,t
                       //  3  011 jjz basis,          j=l-1/2 first, than j=l+1/2 from m=-j to m=j
  } Options;
  Options.ibasis=0;
  OperatorType *op;
  uint32_t Ntemp;
  uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t NShell,l;
  uint16_t *LindexUp, *LindexDn;
  LuaOptionsType LuaOperatorOptions[]={
    {"basis","char*",&(Options.basis)},
    {NULL, NULL, NULL}};
  Options.basis="SphericalHarmonicsY";
  if(NFermion==0) return luaL_error(L, "Angular momentum operator acting on zero fermions is not defined");
  if(LuaTestOptionsList(L, -1))
  {
    LuaCheckOptions(L, -1, LuaOperatorOptions);
    if((strcmp(Options.basis, "SphericalHarmonicsY")==0)||(strcmp(Options.basis, "SphericalHarmonics")==0)||(strcmp(Options.basis, "BasisY")==0)||(strcmp(Options.basis, "Y")==0))
    {
      Options.ibasis=0;
    }
    else if((strcmp(Options.basis, "TesseralHarmonicsZ")==0)||(strcmp(Options.basis, "TesseralHarmonics")==0)||(strcmp(Options.basis, "BasisZ")==0)||(strcmp(Options.basis, "Z")==0))
    {
      Options.ibasis=1;
    }
    else if((strcmp(Options.basis, "KubicHarmonicsK")==0)||(strcmp(Options.basis, "KubicHarmonics")==0)||(strcmp(Options.basis, "BasisK")==0)||(strcmp(Options.basis, "K")==0)||(strcmp(Options.basis, "CubicHarmonicsK")==0)||(strcmp(Options.basis, "CubicHarmonics")==0))
    {
      Options.ibasis=2;
    }
    else if((strcmp(Options.basis, "jjz")==0)||(strcmp(Options.basis, "Basisjjz")==0)||(strcmp(Options.basis, "Basisj")==0)||(strcmp(Options.basis, "j")==0))
    {
      Options.ibasis=3;
    }
    else
    {
      return luaL_error(L, "Unrecognized value for basis %s, use SphericalHarmonicsY, TesseralHarmonicsZ, KubicHarmonicsK, or jjz\n",Options.basis);
    }
  }
  LindexUp = luaL_checkarray_unit16_t(L, 2, &NShell);
  LindexDn = luaL_checkarray_unit16_t(L, 3, &Ntemp);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexUp, NShell)==false) return luaL_error(L, "Index up contains a refference outside the range [0,%d]",NFermion-1);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexDn, Ntemp )==false) return luaL_error(L, "Index dn contains a refference outside the range [0,%d]",NFermion-1);
  if(Options.ibasis==3)
  {
    if((NShell+2)!=Ntemp)
    {
      return luaL_error(L, "Error in LuaOperatorNewLx, Lenght of IndexJmin (%d) is not equal to lenght of IndexJplus-2 (%d)\n",NShell,Ntemp);
    }
    l=NShell/2;
    if((2*l)!=NShell)
    {
      return luaL_error(L, "Error in LuaOperatorNewLx, Lenght of IndexJmin (%d) is not equal to 2l\n",NShell);
    }
  }
  else
  {
    if(NShell!=Ntemp)
    {
      return luaL_error(L, "Error in LuaOperatorNewLx, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
    }
    l=(NShell-1)/2;
    if(2*l+1!=NShell)
    {
      return luaL_error(L, "Error in LuaOperatorNewLx, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
    }
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  switch (Options.ibasis)
  {
    case 0:
      CreateOperatorLx(NFermion, l, LindexUp, LindexDn, op);
      break;
    case 1:
      CreateOperatorLxZ(NFermion, l, LindexUp, LindexDn, op);
      break;
    case 2:
      CreateOperatorLxK(NFermion, l, LindexUp, LindexDn, op);
      break;
    case 3:
      CreateOperatorLxjjz(NFermion, l, LindexUp, LindexDn, op);
      break;
  }
  luaL_setmetatable(L, "Operator_Type");
  free(LindexUp);
  free(LindexDn);
  return 1;
}

static int LuaOperatorNewLy(lua_State * L)
{
  OperatorType *op;
  uint32_t Ntemp;
  uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t NShell,l;
  uint16_t *LindexUp, *LindexDn;
  if(NFermion==0) return luaL_error(L, "Angular momentum operator acting on zero fermions is not defined");
  LindexUp = luaL_checkarray_unit16_t(L, 2, &NShell);
  LindexDn = luaL_checkarray_unit16_t(L, 3, &Ntemp);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexUp, NShell)==false) return luaL_error(L, "Index up contains a refference outside the range [0,%d]",NFermion-1);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexDn, Ntemp )==false) return luaL_error(L, "Index dn contains a refference outside the range [0,%d]",NFermion-1);
  if(NShell!=Ntemp)
  {
    return luaL_error(L, "Error in LuaOperatorNewLy, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
  }
  l=(NShell-1)/2;
  if(2*l+1!=NShell)
  {
    return luaL_error(L, "Error in LuaOperatorNewLy, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorLy(NFermion, l, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "Operator_Type");
  free(LindexUp);
  free(LindexDn);
  return 1;
}

static int LuaOperatorNewLz(lua_State * L)
{
  OperatorType *op;
  uint32_t Ntemp;
  uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t NShell,l;
  uint16_t *LindexUp, *LindexDn;
  if(NFermion==0) return luaL_error(L, "Angular momentum operator acting on zero fermions is not defined");
  LindexUp = luaL_checkarray_unit16_t(L, 2, &NShell);
  LindexDn = luaL_checkarray_unit16_t(L, 3, &Ntemp);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexUp, NShell)==false) return luaL_error(L, "Index up contains a refference outside the range [0,%d]",NFermion-1);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexDn, Ntemp )==false) return luaL_error(L, "Index dn contains a refference outside the range [0,%d]",NFermion-1);
  if(NShell!=Ntemp)
  {
    return luaL_error(L, "Error in LuaOperatorNewLz, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
  }
  l=(NShell-1)/2;
  if(2*l+1!=NShell)
  {
    return luaL_error(L, "Error in LuaOperatorNewLz, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorLz(NFermion, l, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "Operator_Type");
  free(LindexUp);
  free(LindexDn);
  return 1;
}

static int LuaOperatorNewLplus(lua_State * L)
{
  OperatorType *op;
  uint32_t Ntemp;
  uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t NShell,l;
  uint16_t *LindexUp, *LindexDn;
  if(NFermion==0) return luaL_error(L, "Angular momentum operator acting on zero fermions is not defined");
  LindexUp = luaL_checkarray_unit16_t(L, 2, &NShell);
  LindexDn = luaL_checkarray_unit16_t(L, 3, &Ntemp);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexUp, NShell)==false) return luaL_error(L, "Index up contains a refference outside the range [0,%d]",NFermion-1);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexDn, Ntemp )==false) return luaL_error(L, "Index dn contains a refference outside the range [0,%d]",NFermion-1);
  if(NShell!=Ntemp)
  {
    return luaL_error(L, "Error in LuaOperatorNewLplus, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
  }
  l=(NShell-1)/2;
  if(2*l+1!=NShell)
  {
    return luaL_error(L, "Error in LuaOperatorNewLplus, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorLplus(NFermion, l, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "Operator_Type");
  free(LindexUp);
  free(LindexDn);
  return 1;
}

static int LuaOperatorNewLmin(lua_State * L)
{
  OperatorType *op;
  uint32_t Ntemp;
  uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t NShell,l;
  uint16_t *LindexUp, *LindexDn;
  if(NFermion==0) return luaL_error(L, "Angular momentum operator acting on zero fermions is not defined");
  LindexUp = luaL_checkarray_unit16_t(L, 2, &NShell);
  LindexDn = luaL_checkarray_unit16_t(L, 3, &Ntemp);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexUp, NShell)==false) return luaL_error(L, "Index up contains a refference outside the range [0,%d]",NFermion-1);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexDn, Ntemp )==false) return luaL_error(L, "Index dn contains a refference outside the range [0,%d]",NFermion-1);
  if(NShell!=Ntemp)
  {
    return luaL_error(L, "Error in LuaOperatorNewLmin, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
  }
  l=(NShell-1)/2;
  if(2*l+1!=NShell)
  {
    return luaL_error(L, "Error in LuaOperatorNewLmin, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorLmin(NFermion, l, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "Operator_Type");
  free(LindexUp);
  free(LindexDn);
  return 1;
}

static int LuaOperatorNewLsqr(lua_State * L)
{
  OperatorType *op;
  uint32_t Ntemp;
  uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t NShell,l;
  uint16_t *LindexUp, *LindexDn;
  if(NFermion==0) return luaL_error(L, "Angular momentum operator acting on zero fermions is not defined");
  LindexUp = luaL_checkarray_unit16_t(L, 2, &NShell);
  LindexDn = luaL_checkarray_unit16_t(L, 3, &Ntemp);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexUp, NShell)==false) return luaL_error(L, "Index up contains a refference outside the range [0,%d]",NFermion-1);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexDn, Ntemp )==false) return luaL_error(L, "Index dn contains a refference outside the range [0,%d]",NFermion-1);
  if(NShell!=Ntemp)
  {
    return luaL_error(L, "Error in LuaOperatorNewLsqr, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
  }
  l=(NShell-1)/2;
  if(2*l+1!=NShell)
  {
    return luaL_error(L, "Error in LuaOperatorNewLsqr, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorLsqr(NFermion, l, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "Operator_Type");
  free(LindexUp);
  free(LindexDn);
  return 1;
}

static int LuaOperatorNewJx(lua_State * L)
{
  OperatorType *op;
  uint32_t Ntemp;
  uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t NShell,l;
  uint16_t *LindexUp, *LindexDn;
  if(NFermion==0) return luaL_error(L, "Total Angular momentum operator acting on zero fermions is not defined");
  LindexUp = luaL_checkarray_unit16_t(L, 2, &NShell);
  LindexDn = luaL_checkarray_unit16_t(L, 3, &Ntemp);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexUp, NShell)==false) return luaL_error(L, "Index up contains a refference outside the range [0,%d]",NFermion-1);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexDn, Ntemp )==false) return luaL_error(L, "Index dn contains a refference outside the range [0,%d]",NFermion-1);
  if(NShell!=Ntemp)
  {
    return luaL_error(L, "Error in LuaOperatorNewJx, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
  }
  l=(NShell-1)/2;
  if(2*l+1!=NShell)
  {
    return luaL_error(L, "Error in LuaOperatorNewJx, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorJx(NFermion, l, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "Operator_Type");
  free(LindexUp);
  free(LindexDn);
  return 1;
}

static int LuaOperatorNewJy(lua_State * L)
{
  OperatorType *op;
  uint32_t Ntemp;
  uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t NShell,l;
  uint16_t *LindexUp, *LindexDn;
  if(NFermion==0) return luaL_error(L, "Total Angular momentum operator acting on zero fermions is not defined");
  LindexUp = luaL_checkarray_unit16_t(L, 2, &NShell);
  LindexDn = luaL_checkarray_unit16_t(L, 3, &Ntemp);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexUp, NShell)==false) return luaL_error(L, "Index up contains a refference outside the range [0,%d]",NFermion-1);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexDn, Ntemp )==false) return luaL_error(L, "Index dn contains a refference outside the range [0,%d]",NFermion-1);
  if(NShell!=Ntemp)
  {
    return luaL_error(L, "Error in LuaOperatorNewJy, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
  }
  l=(NShell-1)/2;
  if(2*l+1!=NShell)
  {
    return luaL_error(L, "Error in LuaOperatorNewJy, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorJy(NFermion, l, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "Operator_Type");
  free(LindexUp);
  free(LindexDn);
  return 1;
}

static int LuaOperatorNewJz(lua_State * L)
{
  OperatorType *op;
  uint32_t Ntemp;
  uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t NShell,l;
  uint16_t *LindexUp, *LindexDn;
  if(NFermion==0) return luaL_error(L, "Total Angular momentum operator acting on zero fermions is not defined");
  LindexUp = luaL_checkarray_unit16_t(L, 2, &NShell);
  LindexDn = luaL_checkarray_unit16_t(L, 3, &Ntemp);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexUp, NShell)==false) return luaL_error(L, "Index up contains a refference outside the range [0,%d]",NFermion-1);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexDn, Ntemp )==false) return luaL_error(L, "Index dn contains a refference outside the range [0,%d]",NFermion-1);
  if(NShell!=Ntemp)
  {
    return luaL_error(L, "Error in LuaOperatorNewJz, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
  }
  l=(NShell-1)/2;
  if(2*l+1!=NShell)
  {
    return luaL_error(L, "Error in LuaOperatorNewJz, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorJz(NFermion, l, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "Operator_Type");
  free(LindexUp);
  free(LindexDn);
  return 1;
}

static int LuaOperatorNewJplus(lua_State * L)
{
  OperatorType *op;
  uint32_t Ntemp;
  uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t NShell,l;
  uint16_t *LindexUp, *LindexDn;
  if(NFermion==0) return luaL_error(L, "Total Angular momentum operator acting on zero fermions is not defined");
  LindexUp = luaL_checkarray_unit16_t(L, 2, &NShell);
  LindexDn = luaL_checkarray_unit16_t(L, 3, &Ntemp);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexUp, NShell)==false) return luaL_error(L, "Index up contains a refference outside the range [0,%d]",NFermion-1);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexDn, Ntemp )==false) return luaL_error(L, "Index dn contains a refference outside the range [0,%d]",NFermion-1);
  if(NShell!=Ntemp)
  {
    return luaL_error(L, "Error in LuaOperatorNewJplus, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
  }
  l=(NShell-1)/2;
  if(2*l+1!=NShell)
  {
    return luaL_error(L, "Error in LuaOperatorNewJplus, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorJplus(NFermion, l, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "Operator_Type");
  free(LindexUp);
  free(LindexDn);
  return 1;
}

static int LuaOperatorNewJmin(lua_State * L)
{
  OperatorType *op;
  uint32_t Ntemp;
  uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t NShell,l;
  uint16_t *LindexUp, *LindexDn;
  if(NFermion==0) return luaL_error(L, "Total Angular momentum operator acting on zero fermions is not defined");
  LindexUp = luaL_checkarray_unit16_t(L, 2, &NShell);
  LindexDn = luaL_checkarray_unit16_t(L, 3, &Ntemp);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexUp, NShell)==false) return luaL_error(L, "Index up contains a refference outside the range [0,%d]",NFermion-1);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexDn, Ntemp )==false) return luaL_error(L, "Index dn contains a refference outside the range [0,%d]",NFermion-1);
  if(NShell!=Ntemp)
  {
    return luaL_error(L, "Error in LuaOperatorNewJmin, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
  }
  l=(NShell-1)/2;
  if(2*l+1!=NShell)
  {
    return luaL_error(L, "Error in LuaOperatorNewJmin, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorJmin(NFermion, l, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "Operator_Type");
  free(LindexUp);
  free(LindexDn);
  return 1;
}

static int LuaOperatorNewJsqr(lua_State * L)
{
  OperatorType *op;
  uint32_t Ntemp;
  uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t NShell,l;
  uint16_t *LindexUp, *LindexDn;
  if(NFermion==0) return luaL_error(L, "Total Angular momentum operator acting on zero fermions is not defined");
  LindexUp = luaL_checkarray_unit16_t(L, 2, &NShell);
  LindexDn = luaL_checkarray_unit16_t(L, 3, &Ntemp);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexUp, NShell)==false) return luaL_error(L, "Index up contains a refference outside the range [0,%d]",NFermion-1);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexDn, Ntemp )==false) return luaL_error(L, "Index dn contains a refference outside the range [0,%d]",NFermion-1);
  if(NShell!=Ntemp)
  {
    return luaL_error(L, "Error in LuaOperatorNewJsqr, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
  }
  l=(NShell-1)/2;
  if(2*l+1!=NShell)
  {
    return luaL_error(L, "Error in LuaOperatorNewJsqr, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorJsqr(NFermion, l, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "Operator_Type");
  free(LindexUp);
  free(LindexDn);
  return 1;
}

static int LuaOperatorNewTx(lua_State * L)
{
  OperatorType *op;
  uint32_t Ntemp;
  uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t NShell,l;
  uint16_t *LindexUp, *LindexDn;
  if(NFermion==0) return luaL_error(L, "Magnetic dipole operator acting on zero fermions is not defined");
  LindexUp = luaL_checkarray_unit16_t(L, 2, &NShell);
  LindexDn = luaL_checkarray_unit16_t(L, 3, &Ntemp);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexUp, NShell)==false) return luaL_error(L, "Index up contains a refference outside the range [0,%d]",NFermion-1);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexDn, Ntemp )==false) return luaL_error(L, "Index dn contains a refference outside the range [0,%d]",NFermion-1);
  if(NShell!=Ntemp)
  {
    return luaL_error(L, "Error in LuaOperatorNewTx, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
  }
  l=(NShell-1)/2;
  if(2*l+1!=NShell)
  {
    return luaL_error(L, "Error in LuaOperatorNewTx, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorTx(NFermion, l, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "Operator_Type");
  free(LindexUp);
  free(LindexDn);
  return 1;
}

static int LuaOperatorNewTy(lua_State * L)
{
  OperatorType *op;
  uint32_t Ntemp;
  uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t NShell,l;
  uint16_t *LindexUp, *LindexDn;
  if(NFermion==0) return luaL_error(L, "Magnetic dipole operator acting on zero fermions is not defined");
  LindexUp = luaL_checkarray_unit16_t(L, 2, &NShell);
  LindexDn = luaL_checkarray_unit16_t(L, 3, &Ntemp);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexUp, NShell)==false) return luaL_error(L, "Index up contains a refference outside the range [0,%d]",NFermion-1);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexDn, Ntemp )==false) return luaL_error(L, "Index dn contains a refference outside the range [0,%d]",NFermion-1);
  if(NShell!=Ntemp)
  {
    return luaL_error(L, "Error in LuaOperatorNewTy, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
  }
  l=(NShell-1)/2;
  if(2*l+1!=NShell)
  {
    return luaL_error(L, "Error in LuaOperatorNewTy, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorTy(NFermion, l, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "Operator_Type");
  free(LindexUp);
  free(LindexDn);
  return 1;
}

static int LuaOperatorNewTz(lua_State * L)
{
  OperatorType *op;
  uint32_t Ntemp;
  uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t NShell,l;
  uint16_t *LindexUp, *LindexDn;
  if(NFermion==0) return luaL_error(L, "Magnetic dipole operator acting on zero fermions is not defined");
  LindexUp = luaL_checkarray_unit16_t(L, 2, &NShell);
  LindexDn = luaL_checkarray_unit16_t(L, 3, &Ntemp);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexUp, NShell)==false) return luaL_error(L, "Index up contains a refference outside the range [0,%d]",NFermion-1);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexDn, Ntemp )==false) return luaL_error(L, "Index dn contains a refference outside the range [0,%d]",NFermion-1);
  if(NShell!=Ntemp)
  {
    return luaL_error(L, "Error in LuaOperatorNewTz, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
  }
  l=(NShell-1)/2;
  if(2*l+1!=NShell)
  {
    return luaL_error(L, "Error in LuaOperatorNewTz, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorTz(NFermion, l, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "Operator_Type");
  free(LindexUp);
  free(LindexDn);
  return 1;
}

static int LuaOperatorNewTplus(lua_State * L)
{
  OperatorType *op;
  uint32_t Ntemp;
  uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t NShell,l;
  uint16_t *LindexUp, *LindexDn;
  if(NFermion==0) return luaL_error(L, "Magnetic dipole operator acting on zero fermions is not defined");
  LindexUp = luaL_checkarray_unit16_t(L, 2, &NShell);
  LindexDn = luaL_checkarray_unit16_t(L, 3, &Ntemp);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexUp, NShell)==false) return luaL_error(L, "Index up contains a refference outside the range [0,%d]",NFermion-1);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexDn, Ntemp )==false) return luaL_error(L, "Index dn contains a refference outside the range [0,%d]",NFermion-1);
  if(NShell!=Ntemp)
  {
    return luaL_error(L, "Error in LuaOperatorNewTplus, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
  }
  l=(NShell-1)/2;
  if(2*l+1!=NShell)
  {
    return luaL_error(L, "Error in LuaOperatorNewTplus, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorTplus(NFermion, l, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "Operator_Type");
  free(LindexUp);
  free(LindexDn);
  return 1;
}

static int LuaOperatorNewTmin(lua_State * L)
{
  OperatorType *op;
  uint32_t Ntemp;
  uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t NShell,l;
  uint16_t *LindexUp, *LindexDn;
  if(NFermion==0) return luaL_error(L, "Magnetic dipole operator acting on zero fermions is not defined");
  LindexUp = luaL_checkarray_unit16_t(L, 2, &NShell);
  LindexDn = luaL_checkarray_unit16_t(L, 3, &Ntemp);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexUp, NShell)==false) return luaL_error(L, "Index up contains a refference outside the range [0,%d]",NFermion-1);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexDn, Ntemp )==false) return luaL_error(L, "Index dn contains a refference outside the range [0,%d]",NFermion-1);
  if(NShell!=Ntemp)
  {
    return luaL_error(L, "Error in LuaOperatorNewTmin, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
  }
  l=(NShell-1)/2;
  if(2*l+1!=NShell)
  {
    return luaL_error(L, "Error in LuaOperatorNewTmin, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorTmin(NFermion, l, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "Operator_Type");
  free(LindexUp);
  free(LindexDn);
  return 1;
}

static int LuaOperatorNewTsqr(lua_State * L)
{
  OperatorType *op;
  uint32_t Ntemp;
  uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t NShell,l;
  uint16_t *LindexUp, *LindexDn;
  if(NFermion==0) return luaL_error(L, "Magnetic dipole operator acting on zero fermions is not defined");
  LindexUp = luaL_checkarray_unit16_t(L, 2, &NShell);
  LindexDn = luaL_checkarray_unit16_t(L, 3, &Ntemp);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexUp, NShell)==false) return luaL_error(L, "Index up contains a refference outside the range [0,%d]",NFermion-1);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexDn, Ntemp )==false) return luaL_error(L, "Index dn contains a refference outside the range [0,%d]",NFermion-1);
  if(NShell!=Ntemp)
  {
    return luaL_error(L, "Error in LuaOperatorNewTsqr, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
  }
  l=(NShell-1)/2;
  if(2*l+1!=NShell)
  {
    return luaL_error(L, "Error in LuaOperatorNewTsqr, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorTsqr(NFermion, l, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "Operator_Type");
  free(LindexUp);
  free(LindexDn);
  return 1;
}

static int LuaOperatorNewQxx(lua_State * L)
{
  OperatorType *op;
  uint32_t Ntemp;
  uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t NShell,l;
  uint16_t *LindexUp, *LindexDn;
  if(NFermion==0) return luaL_error(L, "Quadrupole operator acting on zero fermions is not defined");
  LindexUp = luaL_checkarray_unit16_t(L, 2, &NShell);
  LindexDn = luaL_checkarray_unit16_t(L, 3, &Ntemp);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexUp, NShell)==false) return luaL_error(L, "Index up contains a refference outside the range [0,%d]",NFermion-1);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexDn, Ntemp )==false) return luaL_error(L, "Index dn contains a refference outside the range [0,%d]",NFermion-1);
  if(NShell!=Ntemp)
  {
    return luaL_error(L, "Error in LuaOperatorNewQxx, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
  }
  l=(NShell-1)/2;
  if(2*l+1!=NShell)
  {
    return luaL_error(L, "Error in LuaOperatorNewQxx, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorQxx(NFermion, l, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "Operator_Type");
  free(LindexUp);
  free(LindexDn);
  return 1;
}

static int LuaOperatorNewQyy(lua_State * L)
{
  OperatorType *op;
  uint32_t Ntemp;
  uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t NShell,l;
  uint16_t *LindexUp, *LindexDn;
  if(NFermion==0) return luaL_error(L, "Quadrupole operator acting on zero fermions is not defined");
  LindexUp = luaL_checkarray_unit16_t(L, 2, &NShell);
  LindexDn = luaL_checkarray_unit16_t(L, 3, &Ntemp);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexUp, NShell)==false) return luaL_error(L, "Index up contains a refference outside the range [0,%d]",NFermion-1);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexDn, Ntemp )==false) return luaL_error(L, "Index dn contains a refference outside the range [0,%d]",NFermion-1);
  if(NShell!=Ntemp)
  {
    return luaL_error(L, "Error in LuaOperatorNewQyy, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
  }
  l=(NShell-1)/2;
  if(2*l+1!=NShell)
  {
    return luaL_error(L, "Error in LuaOperatorNewQyy, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorQyy(NFermion, l, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "Operator_Type");
  free(LindexUp);
  free(LindexDn);
  return 1;
}

static int LuaOperatorNewQzz(lua_State * L)
{
  OperatorType *op;
  uint32_t Ntemp;
  uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t NShell,l;
  uint16_t *LindexUp, *LindexDn;
  if(NFermion==0) return luaL_error(L, "Quadrupole operator acting on zero fermions is not defined");
  LindexUp = luaL_checkarray_unit16_t(L, 2, &NShell);
  LindexDn = luaL_checkarray_unit16_t(L, 3, &Ntemp);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexUp, NShell)==false) return luaL_error(L, "Index up contains a refference outside the range [0,%d]",NFermion-1);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexDn, Ntemp )==false) return luaL_error(L, "Index dn contains a refference outside the range [0,%d]",NFermion-1);
  if(NShell!=Ntemp)
  {
    return luaL_error(L, "Error in LuaOperatorNewQzz, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
  }
  l=(NShell-1)/2;
  if(2*l+1!=NShell)
  {
    return luaL_error(L, "Error in LuaOperatorNewQzz, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorQzz(NFermion, l, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "Operator_Type");
  free(LindexUp);
  free(LindexDn);
  return 1;
}

static int LuaOperatorNewQxy(lua_State * L)
{
  OperatorType *op;
  uint32_t Ntemp;
  uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t NShell,l;
  uint16_t *LindexUp, *LindexDn;
  if(NFermion==0) return luaL_error(L, "Quadrupole operator acting on zero fermions is not defined");
  LindexUp = luaL_checkarray_unit16_t(L, 2, &NShell);
  LindexDn = luaL_checkarray_unit16_t(L, 3, &Ntemp);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexUp, NShell)==false) return luaL_error(L, "Index up contains a refference outside the range [0,%d]",NFermion-1);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexDn, Ntemp )==false) return luaL_error(L, "Index dn contains a refference outside the range [0,%d]",NFermion-1);
  if(NShell!=Ntemp)
  {
    return luaL_error(L, "Error in LuaOperatorNewQxy, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
  }
  l=(NShell-1)/2;
  if(2*l+1!=NShell)
  {
    return luaL_error(L, "Error in LuaOperatorNewQxy, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorQxy(NFermion, l, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "Operator_Type");
  free(LindexUp);
  free(LindexDn);
  return 1;
}

static int LuaOperatorNewQxz(lua_State * L)
{
  OperatorType *op;
  uint32_t Ntemp;
  uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t NShell,l;
  uint16_t *LindexUp, *LindexDn;
  if(NFermion==0) return luaL_error(L, "Quadrupole operator acting on zero fermions is not defined");
  LindexUp = luaL_checkarray_unit16_t(L, 2, &NShell);
  LindexDn = luaL_checkarray_unit16_t(L, 3, &Ntemp);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexUp, NShell)==false) return luaL_error(L, "Index up contains a refference outside the range [0,%d]",NFermion-1);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexDn, Ntemp )==false) return luaL_error(L, "Index dn contains a refference outside the range [0,%d]",NFermion-1);
  if(NShell!=Ntemp)
  {
    return luaL_error(L, "Error in LuaOperatorNewQxz, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
  }
  l=(NShell-1)/2;
  if(2*l+1!=NShell)
  {
    return luaL_error(L, "Error in LuaOperatorNewQxz, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorQxz(NFermion, l, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "Operator_Type");
  free(LindexUp);
  free(LindexDn);
  return 1;
}

static int LuaOperatorNewQyz(lua_State * L)
{
  OperatorType *op;
  uint32_t Ntemp;
  uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t NShell,l;
  uint16_t *LindexUp, *LindexDn;
  if(NFermion==0) return luaL_error(L, "Quadrupole operator acting on zero fermions is not defined");
  LindexUp = luaL_checkarray_unit16_t(L, 2, &NShell);
  LindexDn = luaL_checkarray_unit16_t(L, 3, &Ntemp);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexUp, NShell)==false) return luaL_error(L, "Index up contains a refference outside the range [0,%d]",NFermion-1);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexDn, Ntemp )==false) return luaL_error(L, "Index dn contains a refference outside the range [0,%d]",NFermion-1);
  if(NShell!=Ntemp)
  {
    return luaL_error(L, "Error in LuaOperatorNewQyz, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
  }
  l=(NShell-1)/2;
  if(2*l+1!=NShell)
  {
    return luaL_error(L, "Error in LuaOperatorNewQyz, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorQyz(NFermion, l, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "Operator_Type");
  free(LindexUp);
  free(LindexDn);
  return 1;
}

static int LuaOperatorNewldots(lua_State * L)
{
  OperatorType *op;
  uint32_t Ntemp;
  uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t NShell,l;
  uint16_t *LindexUp, *LindexDn;
  if(NFermion==0) return luaL_error(L, "Spin-orbit operator acting on zero fermions is not defined");
  LindexUp = luaL_checkarray_unit16_t(L, 2, &NShell);
  LindexDn = luaL_checkarray_unit16_t(L, 3, &Ntemp);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexUp, NShell)==false) return luaL_error(L, "Index up contains a refference outside the range [0,%d]",NFermion-1);
  if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexDn, Ntemp )==false) return luaL_error(L, "Index dn contains a refference outside the range [0,%d]",NFermion-1);
  if(NShell!=Ntemp)
  {
    return luaL_error(L, "Error in LuaOperatorNewldots, Lenght of IndexUp (%d) is not equal to lenght of IndexDn (%d)\n",NShell,Ntemp);
  }
  l=(NShell-1)/2;
  if(2*l+1!=NShell)
  {
    return luaL_error(L, "Error in LuaOperatorNewldots, Lenght of IndexUp (%d) is not equal to 2l+1\n",NShell);
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorldots(NFermion, l, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "Operator_Type");
  free(LindexUp);
  free(LindexDn);
  return 1;
}

static int LuaOperatorNewCF(lua_State * L)
{
  OperatorType *op;
  uint32_t i;
  uint32_t len = (uint32_t)lua_gettop(L);
  uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
  uint32_t NShell,Ntemp,lCr,lAn;
  uint16_t *LindexUpCr, *LindexUpAn, *LindexDnCr, *LindexDnAn;
  HarmonicsExpansion *Akm;
  if(NFermion==0) return luaL_error(L, "Crystal-field operator acting on zero fermions is not defined");
  if(len==3)
  {
    LindexUpCr = luaL_checkarray_unit16_t(L, 2, &NShell);
    if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexUpCr, NShell)==false) return luaL_error(L, "Index up contains a refference outside the range [0,%d]",NFermion-1);
    LindexUpAn = LindexUpCr;
    LindexDnAn = LindexUpCr;
    LindexDnCr = LindexUpCr;
    lCr=(NShell-1)/2;
    lAn=lCr;
    if(2*lCr+1!=NShell)
    {
      return luaL_error(L, "Error in LuaOperatorNewCF, Lenght of Index (%d) is not equal to 2l+1 (%d)\n",NShell,lCr);
    }
  }
  else if(len==4)
  {
    LindexUpCr = luaL_checkarray_unit16_t(L, 2, &NShell);
    LindexUpAn = LindexUpCr;
    LindexDnCr = luaL_checkarray_unit16_t(L, 3, &Ntemp);
    LindexDnAn = LindexDnCr;
    if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexUpCr, NShell)==false) return luaL_error(L, "Index up contains a refference outside the range [0,%d]",NFermion-1);
    if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexDnCr, Ntemp )==false) return luaL_error(L, "Index dn contains a refference outside the range [0,%d]",NFermion-1);
    lCr=(NShell-1)/2;
    lAn=lCr;
    if(NShell!=Ntemp)
    {
      return luaL_error(L, "Error in LuaOperatorNewCF, Index of spin up and spin down should relate to the same angular momentum (%d != %d)\n",NShell,Ntemp);
    }
    if(2*lCr+1!=NShell)
    {
      return luaL_error(L, "Error in LuaOperatorNewCF, Lenght of Index (%d) is not equal to 2l+1 (%d)\n",NShell,lCr);
    }
  }
  else if(len==6)
  {
    LindexUpCr = luaL_checkarray_unit16_t(L, 2, &NShell);
    LindexDnCr = luaL_checkarray_unit16_t(L, 3, &Ntemp);
    if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexUpCr, NShell)==false) return luaL_error(L, "Index up contains a refference outside the range [0,%d]",NFermion-1);
    if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexDnCr, Ntemp )==false) return luaL_error(L, "Index dn contains a refference outside the range [0,%d]",NFermion-1);
    lCr=(NShell-1)/2;
    if(NShell!=Ntemp)
    {
      return luaL_error(L, "Error in LuaOperatorNewCF, Index of spin up and spin down should relate to the same angular momentum (%d != %d)\n",NShell,Ntemp);
    }
    if(2*lCr+1!=NShell)
    {
      return luaL_error(L, "Error in LuaOperatorNewCF, Lenght of Index (%d) is not equal to 2l+1 (%d)\n",NShell,lCr);
    }
    LindexUpAn = luaL_checkarray_unit16_t(L, 4, &NShell);
    LindexDnAn = luaL_checkarray_unit16_t(L, 5, &Ntemp);
    if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexUpAn, NShell)==false) return luaL_error(L, "Index up contains a refference outside the range [0,%d]",NFermion-1);
    if(Test_uint16_t_array_in_range(0, (uint16_t)(NFermion-1), LindexDnAn, Ntemp )==false) return luaL_error(L, "Index dn contains a refference outside the range [0,%d]",NFermion-1);
    lAn=(NShell-1)/2;
    if(NShell!=Ntemp)
    {
      return luaL_error(L, "Error in LuaOperatorNewCF, Index of spin up and spin down should relate to the same angular momentum (%d != %d)\n",NShell,Ntemp);
    }
    if(2*lAn+1!=NShell)
    {
      return luaL_error(L, "Error in LuaOperatorNewCF, Lenght of Index (%d) is not equal to 2l+1 (%d)\n",NShell,lAn);
    }
  }
  else
  {
    return luaL_error(L, "Error in LuaOperatorNewCF, called with %d arguments, 3,4 or 6 expected\n",len);
  }
  Akm=luaL_checkarray_HarmonicsExpansion(L, -1);
  if(len==3)
  {
    for(i=0;i<Akm->N;i++)
    {
      Akm->Re[i]*=0.5;
      if(Akm->QComplex) Akm->Im[i]*=0.5;
    }
  }

  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorCF(NFermion, lCr, lAn, LindexUpCr, LindexDnCr, LindexUpAn, LindexDnAn, *Akm, op);
  luaL_setmetatable(L, "Operator_Type");

  free(Akm->l);
  free(Akm->m);
  free(Akm->Re);
  if(Akm->QComplex)  free(Akm->Im);
  free(Akm);

  if(len==3)
  {
    free(LindexUpCr);
  }
  else if(len==4)
  {
    free(LindexUpCr);
    free(LindexDnCr);
  }
  else if(len==6)
  {
    free(LindexUpCr);
    free(LindexDnCr);
    free(LindexUpAn);
    free(LindexDnAn);
  }
  return 1;
}

static int LuaOperatorNewCreate(lua_State * L)
{
  int i;
  OperatorType *op;
  uint32_t len = (uint32_t)lua_gettop(L);
  if(len==2)
  {
    uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
    uint16_t Crindex = (uint16_t)luaL_checkunsigned(L, 2);
    op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    CreateOperatorCreate(NFermion, Crindex, op);
    luaL_setmetatable(L, "Operator_Type");
  }
  else if(len==3)
  {
    uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
    uint32_t NOpp,NNum;
    uint16_t *CrIndex = luaL_checkarray_unit16_t(L, 2, &NOpp);
    Complex16 *Prefactors = luaL_checkarray_complex16(L, 3, &NNum);
    if(NOpp!=NNum)
    {
      return luaL_error(L, "Error in LuaOperatorNewCreate, Lenght of Index (%d) is not equal to lenght of prefactors (%d)\n",NOpp,NNum);
    }
    double PrefactorsRe[NOpp],PrefactorsIm[NOpp];
    for(i=0;i<NOpp;i++)
    {
      PrefactorsRe[i]=Prefactors[i].Re;
      PrefactorsIm[i]=Prefactors[i].Im;
    }
    op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    CreateOperatorCreateList(NFermion, CrIndex, PrefactorsRe, PrefactorsIm, NOpp, op);
    luaL_setmetatable(L, "Operator_Type");
    free(CrIndex);
    free(Prefactors);
  }
  else
  {
    return luaL_error(L, "Error in LuaOperatorNewCreate, expected 2 or 3 arguments, not %d\n",len);
  }
  return 1;
}

static int LuaOperatorNewAnnihilate(lua_State * L)
{
  int i;
  OperatorType *op;
  uint32_t len = (uint32_t)lua_gettop(L);
  if(len==2)
  {
    uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
    uint16_t Crindex = (uint16_t)luaL_checkunsigned(L, 2);
    op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    CreateOperatorAnnihilate(NFermion, Crindex, op);
    luaL_setmetatable(L, "Operator_Type");
  }
  else if(len==3)
  {
    uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
    uint32_t NOpp,NNum;
    uint16_t *CrIndex = luaL_checkarray_unit16_t(L, 2, &NOpp);
    Complex16 *Prefactors = luaL_checkarray_complex16(L, 3, &NNum);
    if(NOpp!=NNum)
    {
      return luaL_error(L, "Error in LuaOperatorNewAnnihilate, Lenght of Index (%d) is not equal to lenght of prefactors (%d)\n",NOpp,NNum);
    }
    double PrefactorsRe[NOpp],PrefactorsIm[NOpp];
    for(i=0;i<NOpp;i++)
    {
      PrefactorsRe[i]=Prefactors[i].Re;
      PrefactorsIm[i]=Prefactors[i].Im;
    }
    op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    CreateOperatorAnnihilateList(NFermion, CrIndex, PrefactorsRe, PrefactorsIm, NOpp, op);
    luaL_setmetatable(L, "Operator_Type");
    free(CrIndex);
    free(Prefactors);
  }
  else
  {
    return luaL_error(L, "Error in LuaOperatorNewAnnihilate, expected 2 or 3 arguments, not %d\n",len);
  }
  return 1;
}

static int LuaOperatorNewIdentity(lua_State * L)
{
  OperatorType *op;
  uint32_t len = (uint32_t)lua_gettop(L);
  if(len==1)
  {
    uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
    op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    CreateOperatorIdentity(NFermion, op);
    luaL_setmetatable(L, "Operator_Type");
    return 1;
  }
  else
  {
    return luaL_error(L, "Error in NewOperator(\"Identity\",...), expected 1 aditional argument, not %d\n",len);
  }
}

static int LuaOperatorNewNumber(lua_State * L)
{
  int i;
  OperatorType *op;
  uint32_t len = (uint32_t)lua_gettop(L);
  if(len==3)
  {
    uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
    uint16_t Crindex = (uint16_t)luaL_checkunsigned(L, 2);
    uint16_t Anindex = (uint16_t)luaL_checkunsigned(L, 3);
    op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    double one=1;
    double zero=0;
    CreateOperatorCreateAnnihilateList(NFermion, &Crindex, &Anindex, &one, &zero, 1, op);
    luaL_setmetatable(L, "Operator_Type");
  }
  else if(len==4)
  {
    uint32_t NFermion = (uint32_t)luaL_checkunsigned(L, 1);
    uint32_t NOppCr,NOppAn,NNum;
    uint16_t *CrIndex = luaL_checkarray_unit16_t(L, 2, &NOppCr);
    uint16_t *AnIndex = luaL_checkarray_unit16_t(L, 3, &NOppAn);
    Complex16 *Prefactors = luaL_checkarray_complex16(L, 4, &NNum);
    if(NOppCr!=NOppAn)
    {
      return luaL_error(L, "Error in LuaOperatorNewNumber, Lenght of Index for Createion (%d) is not equal to lenght of index for Annihilation (%d)\n",NOppCr,NOppAn);
    }
    if(NOppCr!=NNum)
    {
      return luaL_error(L, "Error in LuaOperatorNewNumber, Lenght of Index (%d) is not equal to lenght of prefactors (%d)\n",NOppCr,NNum);
    }
    double PrefactorsRe[NOppCr],PrefactorsIm[NOppCr];
    for(i=0;i<NOppCr;i++)
    {
      PrefactorsRe[i]=Prefactors[i].Re;
      PrefactorsIm[i]=Prefactors[i].Im;
    }
    op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    CreateOperatorCreateAnnihilateList(NFermion, CrIndex, AnIndex, PrefactorsRe, PrefactorsIm, NOppCr, op);
    luaL_setmetatable(L, "Operator_Type");
    free(CrIndex);
    free(AnIndex);
    free(Prefactors);
  }
  else
  {
    return luaL_error(L, "Error in LuaOperatorNewAnnihilate, expected 3 or 4 arguments, not %d\n",len);
  }
  return 1;
}

static int LuaOperatorExp(lua_State * L)
{
  uint32_t len = (uint32_t)lua_gettop(L);
  OperatorType *l = (OperatorType *) luaL_checkudata(L, 1, "Operator_Type");
  OperatorType *res = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  Complex16 alpha={1.0,0.0};
  luaL_setmetatable(L, "Operator_Type");
  if(len==1)
  {
    if(HermitianOperatorNormalOrderExponent(*l, alpha, res, &QDetNotOrderedFcFaBcBa))
    {
      return luaL_error(L, "Error, call to OperatorExp failed\n");
    }
  }
  else if(len==2)
  {
    if (lua_isnumber(L, 2))
    {
      alpha.Re = (double) luaL_checknumber(L, len);
      alpha.Im=0.0;
    }
    else if(luaL_testudata(L, 2, "Complex_Type") != NULL)
    {
      alpha = * (Complex16 *) luaL_checkudata(L, 2, "Complex_Type");
    }
    else
    {
      return luaL_error(L, "Error, call to OperatorExp with 2 arguments, first should be an operator, second is an optional real or complex number. Type found at position two is %s\n",lua_typename(L,lua_type(L,2)));
    }
    if(HermitianOperatorNormalOrderExponent(*l, alpha, res, &QDetNotOrderedFcFaBcBa))
    {
      return luaL_error(L, "Error, call to OperatorExp failed\n");
    }
  }
  else
  {
    return luaL_error(L, "Error, call to OperatorExp with %d arguments, expected 1 or 2\n",len);
  }
  return 1;
}

/***********************/
/* Operator meta actions */
/***********************/

static int LuaOperatorAdd(lua_State * L)
{
  OperatorType *l,*r,*res;
  double a;
  Complex16 *b;
  if (lua_isnumber(L, 1))
  {
    a = (double) luaL_checknumber(L, 1);
    r = (OperatorType *) luaL_checkudata(L, 2, "Operator_Type");
    res = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    luaL_setmetatable(L, "Operator_Type");
    OperatorCopy(res, *r);
    OperatorAddLader(res, NULL, 0, a, 0);
    return 1;
  }
  if (lua_isnumber(L, 2))
  {
    l = (OperatorType *) luaL_checkudata(L, 1, "Operator_Type");
    a = (double) luaL_checknumber(L, 2);
    res = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    luaL_setmetatable(L, "Operator_Type");
    OperatorCopy(res, *l);
    OperatorAddLader(res, NULL, 0, a, 0);
    return 1;
  }
  if (luaL_testudata(L, 1, "Complex_Type") != NULL)
  {
    b = (Complex16 *) luaL_checkudata(L, 1, "Complex_Type");
    r = (OperatorType *) luaL_checkudata(L, 2, "Operator_Type");
    res = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    luaL_setmetatable(L, "Operator_Type");
    OperatorCopy(res, *r);
    OperatorAddLader(res, NULL, 0, b->Re, b->Im);
    return 1;
  }
  if (luaL_testudata(L, 2, "Complex_Type") != NULL)
  {
    l = (OperatorType *) luaL_checkudata(L, 1, "Operator_Type");
    b = (Complex16 *) luaL_checkudata(L, 2, "Complex_Type");
    res = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    luaL_setmetatable(L, "Operator_Type");
    OperatorCopy(res, *l);
    OperatorAddLader(res, NULL, 0, b->Re, b->Im);
    return 1;
  }
  {
    l = (OperatorType *) luaL_checkudata(L, 1, "Operator_Type");
    r = (OperatorType *) luaL_checkudata(L, 2, "Operator_Type");
    if(l->NFermion!=r->NFermion) return luaL_error(L, "Adding two operators with different number of Fermions (%d and %d)\n",l->NFermion,r->NFermion);
    if(l->NBoson!=r->NBoson) return luaL_error(L, "Adding two operators with different number of Bosons (%d and %d)\n",l->NBoson,r->NBoson);
    res = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    luaL_setmetatable(L, "Operator_Type");
    OperatorCopy(res, *l);
    OperatorAdd(res, *r, 1., 0);
    return 1;
  }
}

static int LuaOperatorSub(lua_State * L)
{
  OperatorType *l,*r,*res;
  double a;
  Complex16 *b;
  if (lua_isnumber(L, 1))
  {
    a = (double) luaL_checknumber(L, 1);
    r = (OperatorType *) luaL_checkudata(L, 2, "Operator_Type");
    res = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    luaL_setmetatable(L, "Operator_Type");
    OperatorCopy(res, *r);
    OperatorScale(res, -1, 0);
    OperatorAddLader(res, NULL, 0, a, 0);
    return 1;
  }
  if (lua_isnumber(L, 2))
  {
    l = (OperatorType *) luaL_checkudata(L, 1, "Operator_Type");
    a = (double) luaL_checknumber(L, 2);
    res = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    luaL_setmetatable(L, "Operator_Type");
    OperatorCopy(res, *l);
    OperatorAddLader(res, NULL, 0, -a, 0);
    return 1;
  }
  if (luaL_testudata(L, 1, "Complex_Type") != NULL)
  {
    b = (Complex16 *) luaL_checkudata(L, 1, "Complex_Type");
    r = (OperatorType *) luaL_checkudata(L, 2, "Operator_Type");
    res = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    luaL_setmetatable(L, "Operator_Type");
    OperatorCopy(res, *r);
    OperatorScale(res, -1, 0);
    OperatorAddLader(res, NULL, 0, b->Re, b->Im);
    return 1;
  }
  if (luaL_testudata(L, 2, "Complex_Type") != NULL)
  {
    l = (OperatorType *) luaL_checkudata(L, 1, "Operator_Type");
    b = (Complex16 *) luaL_checkudata(L, 2, "Complex_Type");
    res = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    luaL_setmetatable(L, "Operator_Type");
    OperatorCopy(res, *l);
    OperatorAddLader(res, NULL, 0, -b->Re, -b->Im);
    return 1;
  }
  {
    l = (OperatorType *) luaL_checkudata(L, 1, "Operator_Type");
    r = (OperatorType *) luaL_checkudata(L, 2, "Operator_Type");
    if(l->NFermion!=r->NFermion) return luaL_error(L, "Adding two operators with different number of Fermions (%d and %d)\n",l->NFermion,r->NFermion);
    if(l->NBoson!=r->NBoson) return luaL_error(L, "Adding two operators with different number of Bosons (%d and %d)\n",l->NBoson,r->NBoson);
    res = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    luaL_setmetatable(L, "Operator_Type");
    OperatorCopy(res, *l);
    OperatorAdd(res, *r, -1., 0);
    return 1;
  }
}

/* divide an operator by a number or complex number
 *
 */
static int LuaOperatorDiv(lua_State * L)
{
  if (lua_isnumber(L, 1))
  {
    return luaL_error(L, "Error, number devided by an operator\n");
  }
  if(lua_type(L,1)==LUA_TTABLE)
  {
    return luaL_error(L, "Error, table devided by an operator\n");
  }
  if (lua_isnumber(L, 2))
  {
    double b = (double) luaL_checknumber(L, 2);
    OperatorType *O = (OperatorType*) luaL_checkudata(L, 1, "Operator_Type");
    OperatorType *ONew;
    ONew = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    OperatorCopy(ONew, *O);
    OperatorScale(ONew, 1/b, 0);
    luaL_setmetatable(L, "Operator_Type");
    return 1;
  }
  if (luaL_testudata(L, 2, "Complex_Type") != NULL)
  {
    Complex16 *b = (Complex16 *) luaL_checkudata(L, 2, "Complex_Type");
    OperatorType *O = (OperatorType*) luaL_checkudata(L, 1, "Operator_Type");
    OperatorType *ONew;
    ONew = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    OperatorCopy(ONew, *O);
    OperatorScale(ONew, b->Re/(Sqr(b->Re)+Sqr(b->Im)), -b->Im/(Sqr(b->Re)+Sqr(b->Im)));
    luaL_setmetatable(L, "Operator_Type");
    return 1;
  }
  return luaL_error(L, "Error, operator devided by an object of unrecognized type\n");
}

/*  multiplies two objects, whereby one is an operator
 *  a (complex) number times operator scales the operator - multiplication commutes
 *  a operator times operator is implemented as: A * B - multiplication does not commute
 *  a operator times wavefunction is implemented as O * psi and returns a wavefunction
 *  a wavefunction times operator is implemented as psi * O == ConjugateTranspose(O) * psi
 *  a table times operator is implemented as: {a,b,c} * O == {a * O, b * O, c * O}
 *  a operator times table is implemented as: O * {a,b,c} == {O * a, O * b, O * c}
 */
static int LuaOperatorMul(lua_State * L)
{
  uint32_t len = (uint32_t)lua_gettop(L);
  if(lua_type(L,len)==LUA_TTABLE)
  {
    /* psi * table */
    lua_newtable(L);
    lua_pushnil(L);                /* push a nil to the stack - needed for lua_next */
    while(lua_next(L, len) != 0)   /* push the next index and value to the end of the stack +1 on stack (last element replaced)*/
    {
      LuaOperatorMul(L);           /* multimply the value with the operator (first times last of the stack) +1 on stack*/
      lua_pushvalue(L, -3);        /* push the index to the end of the stack +1 on stack*/
      lua_insert(L, -2);           /* swap index and value */
      lua_settable(L, -5);         /* set new table value -2 on stack */
      lua_pop(L, 1);               /* remove one element from the stack -1 on stack */
    }
    return 1;
  }
  if(lua_type(L,1)==LUA_TTABLE)
  {
    /* psi * table */
    lua_newtable(L);
    lua_pushnil(L);                /* push a nil to the stack - needed for lua_next */
    while(lua_next(L, 1) != 0)     /* push the next index and value to the end of the stack +1 on stack (last element replaced)*/
    {
      lua_insert(L, 1);            /* set value from table to first element on stack */
      lua_pushvalue(L, len+1);     /* push the operator to the end of the stack +1 on stack*/
      LuaOperatorMul(L);           /* multimply the value with the operator (first times last of the stack) +1 on stack*/
      lua_pushvalue(L, -3);        /* push the index to the end of the stack +1 on stack*/
      lua_insert(L, -2);           /* swap index and value */
      lua_settable(L, -5);         /* set new table value -2 on stack */
      lua_pop(L, 1);               /* remove one element from the stack -1 on stack */
      lua_remove(L, 1);            /* remove first element from the stack -1 on stack */
    }
    return 1;
  }
  if (lua_isnumber(L, 1))
  {
    double b = (double) luaL_checknumber(L, 1);
    OperatorType *O = (OperatorType*) luaL_checkudata(L, len, "Operator_Type");
    OperatorType *ONew;
    ONew = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    OperatorCopy(ONew, *O);
    OperatorScale(ONew, b, 0);
    luaL_setmetatable(L, "Operator_Type");
    return 1;
  }
  if (lua_isnumber(L, len))
  {
    double b = (double) luaL_checknumber(L, len);
    OperatorType *O = (OperatorType*) luaL_checkudata(L, 1, "Operator_Type");
    OperatorType *ONew;
    ONew = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    OperatorCopy(ONew, *O);
    OperatorScale(ONew, b, 0);
    luaL_setmetatable(L, "Operator_Type");
    return 1;
  }
  if (luaL_testudata(L, 1, "Complex_Type") != NULL)
  {
    Complex16 *b = (Complex16 *) luaL_checkudata(L, 1, "Complex_Type");
    OperatorType *O = (OperatorType*) luaL_checkudata(L, len, "Operator_Type");
    OperatorType *ONew;
    ONew = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    OperatorCopy(ONew, *O);
    OperatorScale(ONew, b->Re, b->Im);
    luaL_setmetatable(L, "Operator_Type");
    return 1;
  }
  if (luaL_testudata(L, len, "Complex_Type") != NULL)
  {
    Complex16 *b = (Complex16 *) luaL_checkudata(L, len, "Complex_Type");
    OperatorType *O = (OperatorType*) luaL_checkudata(L, 1, "Operator_Type");
    OperatorType *ONew;
    ONew = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    OperatorCopy(ONew, *O);
    OperatorScale(ONew, b->Re, b->Im);
    luaL_setmetatable(L, "Operator_Type");
    return 1;
  }
  if (luaL_testudata(L, len, "Wavefunction_Type") != NULL)
  {
    /* O | psi> */
    WaveFunctionType *psi = (WaveFunctionType*) luaL_checkudata(L, len, "Wavefunction_Type");
    OperatorType *O = (OperatorType *) luaL_checkudata(L,1,"Operator_Type");
    WaveFunctionType *Opsi;
    Opsi = (WaveFunctionType *) lua_newuserdata(L, sizeof(WaveFunctionType));
    OperatorPsi(*O, *psi, Opsi);
    luaL_setmetatable(L, "Wavefunction_Type");
    return 1;
  }
  if (luaL_testudata(L, 1, "Wavefunction_Type") != NULL)
  {
    // this part never excecutes while it is handled by the wave-function multiplication rules
    /* <psi | O */
    WaveFunctionType *psi = (WaveFunctionType*) luaL_checkudata(L, 1, "Wavefunction_Type");
    OperatorType *O = (OperatorType *) luaL_checkudata(L,len,"Operator_Type");
    OperatorType OT;
    OperatorConjugateTransposeNormalOrder(*O, &OT, &QDetNotOrderedFcFaBcBa);
    WaveFunctionType *Opsi;
    Opsi = (WaveFunctionType *) lua_newuserdata(L, sizeof(WaveFunctionType));
    OperatorPsi(OT, *psi, Opsi);
    FreeOperator(&OT);
    luaL_setmetatable(L, "Wavefunction_Type");
    return 1;
  }
  {
    /* O * O */
    OperatorType *l = (OperatorType *) luaL_checkudata(L, 1, "Operator_Type");
    OperatorType *r = (OperatorType *) luaL_checkudata(L, len, "Operator_Type");
    OperatorType *res = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    luaL_setmetatable(L, "Operator_Type");
    InitOperatorForMult(res, 2, l, r);
    OperatorMultiply2NormalOrder(res, *l, *r, 'N', 'N', &QDetNotOrderedFcFaBcBa);
  }
  return 1;
}

static int LuaOperatorUnm(lua_State * L)
{
  OperatorType *res, *a;
  a = (OperatorType *) luaL_checkudata(L, 1, "Operator_Type");
  res = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  OperatorCopy(res, *a);
  OperatorScale(res, -1, 0);
  luaL_setmetatable(L, "Operator_Type");
  return 1;
}

static int LuaOperatorPrint(lua_State * L)
{
  struct
  {
    char *stream;
    bool Full;
  } Options;
  uint32_t len = (uint32_t)lua_gettop(L);
  LuaOptionsType LuaOperatorPrintOptions[]={
    {"file", "char*", &(Options.stream)},
    {"Full", "bool", &(Options.Full)},
    {NULL, NULL, NULL}};
  if(!((len==0)||(len==1)))
  {
    return luaL_error(L, "Error in .Print(arg) unrecognized arguments in Print()\n");
  }
  if(len==0)
  {
    PrintOperatorReadable(*((OperatorType *) currentudata));
    return 0;
  }
  Options.Full=false;
  Options.stream=NULL;
  LuaCheckOptions(L, -1, LuaOperatorPrintOptions);
  if(Options.stream==NULL)
  {
    if(Options.Full==false)
    {
      PrintOperatorReadable(*((OperatorType *) currentudata));
    }
    else
    {
      PrintOperatorFullReadable(*((OperatorType *) currentudata));
    }
  }
  else
  {
    if(Options.Full==false)
    {
      FilePrintOperatorReadable(*((OperatorType *) currentudata), Options.stream, "w");
    }
    else
    {
      FilePrintOperatorFullReadable(*((OperatorType *) currentudata), Options.stream, "w");
    }
  }
  return 0;
}

static int LuaOperatorChop(lua_State * L)
{
  uint32_t len = (uint32_t)lua_gettop(L);
  double border;
  switch (len)
  {
    case 0:
      border=10*DBL_EPSILON;
    break;
    case 1:
      border=luaL_checknumber(L, 1);
    break;
    default:
      return luaL_error(L, "Errorr, call to chop with %d arguments, expected 0 or 1\n",len);
    break;
  }
  OperatorChop((OperatorType *) currentudata, border);
  OperatorSetMaxLength((OperatorType *) currentudata);
  return 0;;
}

static int LuaOperatorSetName(lua_State * L)
{
  uint32_t len = (uint32_t)lua_gettop(L);
  switch (len)
  {
    case 0:
      OperatorSetName((OperatorType *) currentudata, "");
    break;
    case 1:
      OperatorSetName((OperatorType *) currentudata, luaL_checkstring(L,1));
    break;
    default:
      return luaL_error(L, "Errorr, call to SetName with %d arguments, expected 0 or 1\n",len);
    break;
  }
  return 0;
}

static int LuaOperatorConjugate(lua_State * L)
{
  uint32_t len = (uint32_t)lua_gettop(L);
  if(len==0)
  {
    OperatorComplexConjugate((OperatorType *)currentudata);
  }
  else
  {
    return luaL_error(L, "Error in .Conjugate(arg) unrecognized arguments in Conjugate()\n");
  }
  return 0;
}

static int LuaOperatorConjugateTranspose(lua_State * L)
{
  //TODO: write the function OperatorTransposeInline
  uint32_t len = (uint32_t)lua_gettop(L);
  if(len==0)
  {
    OperatorType OppTemp;
    OperatorCopy(&OppTemp, *(OperatorType *)currentudata);
    FreeOperator((OperatorType *)currentudata);
    OperatorConjugateTransposeNormalOrder(OppTemp, (OperatorType *)currentudata, &QDetNotOrderedFcFaBaBc);
    FreeOperator(&OppTemp);
  }
  else
  {
    return luaL_error(L, "Error in .ConjugateTranspose(arg) unrecognized arguments in ConjugateTranspose()\n");
  }
  return 0;
}

static int LuaOperatorTranspose(lua_State * L)
{
  //TODO: write the function OperatorTransposeInline
  uint32_t len = (uint32_t)lua_gettop(L);
  if(len==0)
  {
    OperatorType OppTemp;
    OperatorCopy(&OppTemp, *(OperatorType *)currentudata);
    FreeOperator((OperatorType *)currentudata);
    OperatorTransposeNormalOrder(OppTemp, (OperatorType *)currentudata, &QDetNotOrderedFcFaBaBc);
    FreeOperator(&OppTemp);
  }
  else
  {
    return luaL_error(L, "Error in .ConjugateTranspose(arg) unrecognized arguments in ConjugateTranspose()\n");
  }
  return 0;
}

static int LuaOperatorIndex(lua_State * L)
{
  OperatorType *a = (OperatorType *) luaL_checkudata(L, 1, "Operator_Type");
  const char *str = lua_tostring(L, 2);
  currentudata = a;

  if (strcmp(str, "Print") == 0)
  {
    lua_pushcfunction(L, LuaOperatorPrint);
    return 1;
  }
  if (strcmp(str, "Chop") == 0)
  {
    lua_pushcfunction(L, LuaOperatorChop);
    return 1;
  }
  if ((strcmp(str, "Conjugate") == 0)||(strcmp(str, "C") == 0)||(strcmp(str, "Dag") == 0))
  {
    lua_pushcfunction(L, LuaOperatorConjugate);
    return 1;
  }
  if ((strcmp(str, "ConjugateTranspose") == 0)||(strcmp(str, "CT")==0))
  {
    lua_pushcfunction(L, LuaOperatorConjugateTranspose);
    return 1;
  }
  if ((strcmp(str, "Transpose") == 0)||(strcmp(str, "T")==0))
  {
    lua_pushcfunction(L, LuaOperatorTranspose);
    return 1;
  }
  if (strcmp(str, "Name") == 0)
  {
    lua_pushstring(L, a->Name);
    return 1;
  }
  if (strcmp(str, "NF") == 0)
  {
    lua_pushinteger(L, a->NFermion);
    return 1;
  }
  if (strcmp(str, "NB") == 0)
  {
    lua_pushinteger(L, a->NBoson);
    return 1;
  }
  if (strcmp(str, "Complex") == 0)
  {
    lua_pushboolean(L, a->QComplex);
    return 1;
  }
  if (strcmp(str, "N") == 0)
  {
    uint32_t i=0;
    lua_createtable(L, a->MaxLength, 1);
    lua_pushnumber(L, i);
    lua_pushnumber(L, a->l[i].N);
    lua_settable(L, -3);
    for(i=1;i<=a->MaxLength;i++)
    {
      lua_pushnumber(L, a->l[i].N);
      lua_rawseti(L, -2, i);
    }
    return 1;
  }
  if (strcmp(str, "Restrictions") == 0)
  {
    uint32_t i;
    char *s;
    if (a->Res.N==0)
    {
      lua_pushnil(L);
      return 1;
    }
    lua_createtable(L, a->Res.N+2, 0);
    lua_pushnumber(L, a->NFermion);
    lua_rawseti(L, -2, 1);
    lua_pushnumber(L, a->NBoson);
    lua_rawseti(L, -2, 2);
    for(i=0;i<a->Res.N;i++)
    {
      lua_createtable(L, 3, 0);
      s=DetToString(a->NFermion, a->NBoson, a->Res.Determinant[i]);
      lua_pushstring(L,s);
      free(s);
      lua_rawseti(L, -2, 1);
      lua_pushnumber(L, a->Res.Min[i]);
      lua_rawseti(L, -2, 2);
      lua_pushnumber(L, a->Res.Max[i]);
      lua_rawseti(L, -2, 3);
      lua_rawseti(L, -2, i+3);
    }
    return 1;
  }
  return luaL_error(L, "Error operator variable not found");;
}

static int LuaOperatorNewIndex(lua_State * L)
{
  OperatorType *a = (OperatorType *) luaL_checkudata(L, 1, "Operator_Type");
  const char *str = lua_tostring(L, 2);
  if(strcmp(str, "Name")==0)
  {
    OperatorSetName(a, lua_tostring(L, 3));
    return 0;
  }
  if(strcmp(str, "Complex")==0)
  {
    return luaL_error(L, "Errorr, you can not change the value Complex of an operator %s.\nIn order to change you have to modify the operator itselves, the value of Complex will follow automatically.",a->Name);
  }
  if(strcmp(str, "N")==0)
  {
    return luaL_error(L, "Errorr, you can not change the value representing the number of operators in %s\nIf you add or subtract something from the operator the number will follow automatically.",a->Name);
  }
  if(strcmp(str, "NB")==0)
  {
    uint32_t NBnew;
    NBnew = (uint32_t)luaL_checknumber(L, 3);
    if(NBnew+a->NFermion>32768) return luaL_error(L, "Maximum number of fermions + bosons allowed is 2^15 = 32768\nTrying to change the number of bosons from %d to %d with %d fermions\n",a->NBoson,NBnew,a->NFermion);
    OperatorChangeNFandNB(a, a->NFermion, NBnew);
    //TODO: allow bosons - remove this line when debuged
    return luaL_error(L, "Current version does not suport Bosons\n");
    return 0;
  }
  if(strcmp(str, "NF")==0)
  {
    uint32_t NFnew;
    NFnew = (uint32_t)luaL_checknumber(L, 3);
    if(NFnew+a->NBoson>32768) return luaL_error(L, "Maximum number of fermions + bosons allowed is 2^15 = 32768\nTrying to change the number of fermions from %d to %d with %d bosons\n",a->NFermion,NFnew,a->NBoson);
    OperatorChangeNFandNB(a, NFnew, a->NBoson);
    return 0;
  }
  if (strcmp(str, "Restrictions") == 0)
  {
    uint32_t i;
    for(i=0;i<a->Res.N;i++)
    {
      free(a->Res.Determinant[i]);
    }
    if(a->Res.N!=0)
    {
      free(a->Res.Min);
      free(a->Res.Max);
      free(a->Res.Determinant);
    }
    a->Res.N=0;
    if(!lua_isnil(L,3))
    {
      // check if NF and NB agree with size of operator and through an error if not

      lua_rawgeti(L, 3, 1);
      if(a->NFermion!=(uint32_t)luaL_checkunsigned(L, -1))
      {
        return luaL_error(L, "Number of Fermions in definition of restrictions (%d) not the same as number of Fermions in operator (%d)\n",(uint32_t)luaL_checkunsigned(L, -1),a->NFermion);
      }
      lua_pop(L, 1);
      lua_rawgeti(L, 3, 2);
      if(a->NBoson!=(uint32_t)luaL_checkunsigned(L, -1))
      {
        return luaL_error(L, "Number of Bosons in definition of restrictions (%d) not the same as number of Bosons in operator (%d)\n",(uint32_t)luaL_checkunsigned(L, -1),a->NBoson);
      }
      lua_pop(L, 1);
      luaL_setarray_DeterminatRestrictionType(L,3,&(a->Res));
    }
    return 0;
  }
  return luaL_error(L, "Error index %s not found in Operator Object\n",str);
}

static int LuaOperatorToString(lua_State * L)
{
  OperatorType *res = (OperatorType *) luaL_checkudata(L, 1, "Operator_Type");
  char *Buffer;
  StringPrintOperatorReadable(*res, &Buffer);
  lua_pushstring(L, Buffer);
  return 1;
}

static int LuaOperatorDestroy(lua_State * L)
{
  OperatorType *a = (OperatorType *) luaL_checkudata(L, 1, "Operator_Type");
  FreeOperator(a);
  return 0;
}

static const struct luaL_Reg Operator_meta[] = {
  {"__add", LuaOperatorAdd},
  {"__sub", LuaOperatorSub},
  {"__mul", LuaOperatorMul},
  {"__div", LuaOperatorDiv},
  {"__unm", LuaOperatorUnm},
  {"__index", LuaOperatorIndex},
  {"__newindex", LuaOperatorNewIndex},
  {"__tostring", LuaOperatorToString},
  {"__gc", LuaOperatorDestroy},
  {NULL, NULL}
};
