/*  These functions transform matrices to operators for some standard models
 *
 */

/*   There are a few similar versions. The funcitons XXXToAndersonImpurityOperator are "old" and transform the matrix
 *   to an impurity matrix in the order where we first have all spin up and then all spin down
 *
 */
 
/*
 *
 */
int TriDiagonalMatrixToImpurityHamiltonian(OperatorType * O, TriDiagonalMatrixType * M, bool AddSpin, uint32_t FermionShift, uint32_t NF)
{
  int i,k;
  uint16_t oadd[2];
  InitNamedOperator(O, "Single Site Tridiagonal Impurity Hamiltonian");
  O->NFermion = FermionShift + M->N;
  if(AddSpin==true) O->NFermion = FermionShift + 2*M->N;
  if(NF!=0)
  {
    if(NF<O->NFermion)
    {
      printf("ERROR in TriDiagonalMatrixToImpurityHamiltonian. NF (%d) < NFermion (%d)\n",NF,O->NFermion);
      return 1;
    }
    else
    {
      O->NFermion = NF;
    }
  }
  if(M->N==0) return 0;
  O->NWord = Min(O->NFermion, 1 + (O->NFermion - 1) / 8);
  InitOperatorFixedLength(O, 2, 8);
  if(AddSpin==true)
  {
    k=FermionShift;
    for(i=1;i<M->N;i++)
    {
      oadd[1] = CrAnMask + k;
      oadd[0] = k;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, M->a[i], &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + k+2;
      oadd[0] = k;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, M->b[i], &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + k;
      oadd[0] = k+2;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, M->b[i], &QDetNotOrderedFcFaBcBa);
      k++;
      oadd[1] = CrAnMask + k;
      oadd[0] = k;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, M->a[i], &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + k+2;
      oadd[0] = k;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, M->b[i], &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + k;
      oadd[0] = k+2;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, M->b[i], &QDetNotOrderedFcFaBcBa);
      k++;
    }
    oadd[1] = CrAnMask + k;
    oadd[0] = k;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, M->a[M->N], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + k +1;
    oadd[0] = k +1;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, M->a[M->N], &QDetNotOrderedFcFaBcBa);
  }
  else
  {
    k=FermionShift;
    for(i=1;i<M->N;i++)
    {
      oadd[1] = CrAnMask + k;
      oadd[0] = k;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, M->a[i], &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + k+1;
      oadd[0] = k;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, M->b[i], &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + k;
      oadd[0] = k+1;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, M->b[i], &QDetNotOrderedFcFaBcBa);
      k++;
    }
    oadd[1] = CrAnMask + k;
    oadd[0] = k;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, M->a[M->N], &QDetNotOrderedFcFaBcBa);
  }
  return 0;
}
int AndersonMatrixToImpurityHamiltonian(OperatorType * O, AndersonMatrixType * M, bool AddSpin, uint32_t FermionShift, uint32_t NF)
{
  int i,k;
  uint16_t oadd[2];
  InitNamedOperator(O, "Single Site Anderson Impurity Hamiltonian");
  O->NFermion = FermionShift + M->N;
  if(AddSpin==true) O->NFermion = FermionShift + 2*M->N;
  if(NF!=0)
  {
    if(NF<O->NFermion)
    {
      printf("ERROR in AndersonMatrixToImpurityHamiltonian. NF (%d) < NFermion (%d)\n",NF,O->NFermion);
      return 1;
    }
    else
    {
      O->NFermion = NF;
    }
  }
  if(M->N==0) return 0;
  O->NWord = Min(O->NFermion, 1 + (O->NFermion - 1) / 8);
  InitOperatorFixedLength(O, 2, 8);
  if(AddSpin==true)
  {
    k=FermionShift;
    for(i=1;i<M->N;i++)
    {
      oadd[1] = CrAnMask + k;
      oadd[0] = k;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, M->a[i], &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + k+2;
      oadd[0] = FermionShift;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, M->b[i], &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + FermionShift;
      oadd[0] = k+2;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, M->b[i], &QDetNotOrderedFcFaBcBa);
      k++;
      oadd[1] = CrAnMask + k;
      oadd[0] = k;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, M->a[i], &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + k+2;
      oadd[0] = FermionShift+1;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, M->b[i], &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + FermionShift+1;
      oadd[0] = k+2;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, M->b[i], &QDetNotOrderedFcFaBcBa);
      k++;
    }
    oadd[1] = CrAnMask + k;
    oadd[0] = k;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, M->a[M->N], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + k +1;
    oadd[0] = k +1;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, M->a[M->N], &QDetNotOrderedFcFaBcBa);
  }
  else
  {
    k=FermionShift;
    for(i=1;i<M->N;i++)
    {
      oadd[1] = CrAnMask + k;
      oadd[0] = k;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, M->a[i], &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + k+1;
      oadd[0] = FermionShift;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, M->b[i], &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + FermionShift;
      oadd[0] = k+1;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, M->b[i], &QDetNotOrderedFcFaBcBa);
      k++;
    }
    oadd[1] = CrAnMask + k;
    oadd[0] = k;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, M->a[M->N], &QDetNotOrderedFcFaBcBa);
  }
  return 0;
}
int BADoubleTriDiagonalMatrixToImpurityHamiltonian(OperatorType * O, BADoubleTriDiagonalMatrixType * M, bool AddSpin, uint32_t FermionShift, uint32_t NF)
{
  int i,k;
  uint16_t oadd[2];
  InitNamedOperator(O, "Single Site Natural Impurity Hamiltonian");
  O->NFermion = FermionShift + M->val.N + M->con.N;
  if(AddSpin==true) O->NFermion = FermionShift + 2*(M->val.N+M->con.N);
  if(NF!=0)
  {
    if(NF<O->NFermion)
    {
      printf("ERROR in TriDiagonalMatrixToImpurityHamiltonian. NF (%d) < NFermion (%d)\n",NF,O->NFermion);
      return 1;
    }
    else
    {
      O->NFermion = NF;
    }
  }
  if((M->val.N==0) || (M->con.N)==0) return 0;
  O->NWord = Min(O->NFermion, 1 + (O->NFermion - 1) / 8);
  InitOperatorFixedLength(O, 2, 8);
  if(AddSpin==true)
  {
    oadd[1] = CrAnMask + 0 + FermionShift;
    oadd[0] = 0 + FermionShift;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, M->acvB,  &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 1 + FermionShift;
    oadd[0] = 1 + FermionShift;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, M->acvB,  &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 2 + FermionShift;
    oadd[0] = 2 + FermionShift;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, M->acvA,  &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 3 + FermionShift;
    oadd[0] = 3 + FermionShift;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, M->acvA,  &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 0 + FermionShift;
    oadd[0] = 2 + FermionShift;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, M->acvBA, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 2 + FermionShift;
    oadd[0] = 0 + FermionShift;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, M->acvBA, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 1 + FermionShift;
    oadd[0] = 3 + FermionShift;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, M->acvBA, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 3 + FermionShift;
    oadd[0] = 1 + FermionShift;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, M->acvBA, &QDetNotOrderedFcFaBcBa);
    
    k=FermionShift+4;
    
    if(M->val.N>1)
    {
      oadd[1] = CrAnMask + 4 + FermionShift;
      oadd[0] = 0 + FermionShift;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, M->bvB0, &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + 0 + FermionShift;
      oadd[0] = 4 + FermionShift;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, M->bvB0, &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + 5 + FermionShift;
      oadd[0] = 1 + FermionShift;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, M->bvB0, &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + 1 + FermionShift;
      oadd[0] = 5 + FermionShift;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, M->bvB0, &QDetNotOrderedFcFaBcBa);
      
      oadd[1] = CrAnMask + 4 + FermionShift;
      oadd[0] = 2 + FermionShift;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, M->bvA0, &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + 2 + FermionShift;
      oadd[0] = 4 + FermionShift;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, M->bvA0, &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + 5 + FermionShift;
      oadd[0] = 3 + FermionShift;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, M->bvA0, &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + 3 + FermionShift;
      oadd[0] = 5 + FermionShift;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, M->bvA0, &QDetNotOrderedFcFaBcBa);
      for(i=2;i<M->val.N;i++)
      {
        oadd[1] = CrAnMask + k;
        oadd[0] = k;
        RealOperatorAddLaderNormalOrder(O, oadd, 2, M->val.a[i], &QDetNotOrderedFcFaBcBa);
        oadd[1] = CrAnMask + k+2;
        oadd[0] = k;
        RealOperatorAddLaderNormalOrder(O, oadd, 2, M->val.b[i], &QDetNotOrderedFcFaBcBa);
        oadd[1] = CrAnMask + k;
        oadd[0] = k+2;
        RealOperatorAddLaderNormalOrder(O, oadd, 2, M->val.b[i], &QDetNotOrderedFcFaBcBa);
        k++;
        oadd[1] = CrAnMask + k;
        oadd[0] = k;
        RealOperatorAddLaderNormalOrder(O, oadd, 2, M->val.a[i], &QDetNotOrderedFcFaBcBa);
        oadd[1] = CrAnMask + k+2;
        oadd[0] = k;
        RealOperatorAddLaderNormalOrder(O, oadd, 2, M->val.b[i], &QDetNotOrderedFcFaBcBa);
        oadd[1] = CrAnMask + k;
        oadd[0] = k+2;
        RealOperatorAddLaderNormalOrder(O, oadd, 2, M->val.b[i], &QDetNotOrderedFcFaBcBa);
        k++;
      }
      oadd[1] = CrAnMask + k;
      oadd[0] = k;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, M->val.a[M->val.N], &QDetNotOrderedFcFaBcBa);
      k++;
      oadd[1] = CrAnMask + k;
      oadd[0] = k;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, M->val.a[M->val.N], &QDetNotOrderedFcFaBcBa);
      k++;
    }
    
    if(M->con.N>1)
    {
      oadd[1] = CrAnMask + 2*M->val.N + 2 + FermionShift;
      oadd[0] = 0 + FermionShift;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, M->bcB0, &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + 0 + FermionShift;
      oadd[0] = 2*M->val.N + 2 + FermionShift;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, M->bcB0, &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + 2*M->val.N + 3 + FermionShift;
      oadd[0] = 1 + FermionShift;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, M->bcB0, &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + 1 + FermionShift;
      oadd[0] = 2*M->val.N + 3 + FermionShift;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, M->bcB0, &QDetNotOrderedFcFaBcBa);
      
      oadd[1] = CrAnMask + 2*M->val.N + 2 + FermionShift;
      oadd[0] = 2 + FermionShift;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, M->bcA0, &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + 2 + FermionShift;
      oadd[0] = 2*M->val.N + 2 + FermionShift;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, M->bcA0, &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + 2*M->val.N + 3 + FermionShift;
      oadd[0] = 3 + FermionShift;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, M->bcA0, &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + 3 + FermionShift;
      oadd[0] = 2*M->val.N + 3 + FermionShift;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, M->bcA0, &QDetNotOrderedFcFaBcBa);
      for(i=2;i<M->con.N;i++)
      {
        oadd[1] = CrAnMask + k;
        oadd[0] = k;
        RealOperatorAddLaderNormalOrder(O, oadd, 2, M->con.a[i], &QDetNotOrderedFcFaBcBa);
        oadd[1] = CrAnMask + k+2;
        oadd[0] = k;
        RealOperatorAddLaderNormalOrder(O, oadd, 2, M->con.b[i], &QDetNotOrderedFcFaBcBa);
        oadd[1] = CrAnMask + k;
        oadd[0] = k+2;
        RealOperatorAddLaderNormalOrder(O, oadd, 2, M->con.b[i], &QDetNotOrderedFcFaBcBa);
        k++;
        oadd[1] = CrAnMask + k;
        oadd[0] = k;
        RealOperatorAddLaderNormalOrder(O, oadd, 2, M->con.a[i], &QDetNotOrderedFcFaBcBa);
        oadd[1] = CrAnMask + k+2;
        oadd[0] = k;
        RealOperatorAddLaderNormalOrder(O, oadd, 2, M->con.b[i], &QDetNotOrderedFcFaBcBa);
        oadd[1] = CrAnMask + k;
        oadd[0] = k+2;
        RealOperatorAddLaderNormalOrder(O, oadd, 2, M->con.b[i], &QDetNotOrderedFcFaBcBa);
        k++;
      }
      oadd[1] = CrAnMask + k;
      oadd[0] = k;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, M->con.a[M->con.N], &QDetNotOrderedFcFaBcBa);
      k++;
      oadd[1] = CrAnMask + k;
      oadd[0] = k;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, M->con.a[M->con.N], &QDetNotOrderedFcFaBcBa);
    }
  }
  else
  {
    oadd[1] = CrAnMask + 0 + FermionShift;
    oadd[0] = 0 + FermionShift;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, M->acvB,  &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 1 + FermionShift;
    oadd[0] = 1 + FermionShift;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, M->acvA,  &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 0 + FermionShift;
    oadd[0] = 1 + FermionShift;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, M->acvBA, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 1 + FermionShift;
    oadd[0] = 0 + FermionShift;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, M->acvBA, &QDetNotOrderedFcFaBcBa);
    k=FermionShift+2;
    if(M->val.N>1)
    {
      oadd[1] = CrAnMask + 2 + FermionShift;
      oadd[0] = 0 + FermionShift;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, M->bvB0, &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + 0 + FermionShift;
      oadd[0] = 2 + FermionShift;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, M->bvB0, &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + 2 + FermionShift;
      oadd[0] = 1 + FermionShift;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, M->bvA0, &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + 1 + FermionShift;
      oadd[0] = 2 + FermionShift;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, M->bvA0, &QDetNotOrderedFcFaBcBa);
      
      for(i=2;i<M->val.N;i++)
      {
        oadd[1] = CrAnMask + k;
        oadd[0] = k;
        RealOperatorAddLaderNormalOrder(O, oadd, 2, M->val.a[i], &QDetNotOrderedFcFaBcBa);
        oadd[1] = CrAnMask + k+1;
        oadd[0] = k;
        RealOperatorAddLaderNormalOrder(O, oadd, 2, M->val.b[i], &QDetNotOrderedFcFaBcBa);
        oadd[1] = CrAnMask + k;
        oadd[0] = k+1;
        RealOperatorAddLaderNormalOrder(O, oadd, 2, M->val.b[i], &QDetNotOrderedFcFaBcBa);
        k++;
      }
      oadd[1] = CrAnMask + k;
      oadd[0] = k;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, M->val.a[M->val.N], &QDetNotOrderedFcFaBcBa);
      k++;
    }
    if(M->con.N>1)
    {
      oadd[1] = CrAnMask + M->val.N + 1 + FermionShift;
      oadd[0] = 0 + FermionShift;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, M->bcB0, &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + 0 + FermionShift;
      oadd[0] = M->val.N + 1 + FermionShift;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, M->bcB0, &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + M->val.N + 1 + FermionShift;
      oadd[0] = 1 + FermionShift;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, M->bcA0, &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + 1 + FermionShift;
      oadd[0] = M->val.N + 1 + FermionShift;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, M->bcA0, &QDetNotOrderedFcFaBcBa);
      for(i=2;i<M->con.N;i++)
      {
        oadd[1] = CrAnMask + k;
        oadd[0] = k;
        RealOperatorAddLaderNormalOrder(O, oadd, 2, M->con.a[i], &QDetNotOrderedFcFaBcBa);
        oadd[1] = CrAnMask + k+1;
        oadd[0] = k;
        RealOperatorAddLaderNormalOrder(O, oadd, 2, M->con.b[i], &QDetNotOrderedFcFaBcBa);
        oadd[1] = CrAnMask + k;
        oadd[0] = k+1;
        RealOperatorAddLaderNormalOrder(O, oadd, 2, M->con.b[i], &QDetNotOrderedFcFaBcBa);
        k++;
      }
      oadd[1] = CrAnMask + k;
      oadd[0] = k;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, M->con.a[M->con.N], &QDetNotOrderedFcFaBcBa);
    }
  }
  return 0;
}


/*  One impurity site with a valence and a conduction chain hanging off it, as an operator.
 *  The single valued twin of BlockDoubleTriDiagonalMatrixToImpurityHamiltonian, and this
 *  type is real throughout, so none of the conjugates of the block version arise.
 *
 *  The orbitals are
 *
 *      [ core | impurity | val_1 ... val_Nv | con_1 ... con_Nc ]
 *
 *  which is the layout ResponseFunction.CreateBitMasks reports and the one
 *  DoubleTriDiagonalMatrixToCompactMatrix builds. a[0] of a chain is the impurity energy --
 *  the two chains share the site, so it is val.a[0] + con.a[0] -- and b[0] is the hopping
 *  from the impurity into the first site of that chain, with a[1] to a[N] and b[1] to
 *  b[N-1] the chain itself. A chain therefore contributes N sites, one more per chain than
 *  the natural orbital form puts in its chains, where a[1] is the natural orbital itself.
 *
 *  That is what makes the loop below read b[i-1] for the bond into site i: the site before
 *  val_i is val_{i-1} for i > 1 and the impurity for i = 1, and the bond between them is
 *  b[i-1] in both cases.
 *
 *  M->a0 and M->b0 are the offset and the scaling of the Green's function,
 *  G = a0 + b0^2 [1/(w-H)]_00, and are not part of H, as in every other
 *  …ToImpurityHamiltonian.
 *
 *  An empty chain is a normal case for this type rather than an edge case -- a chemical
 *  potential outside the bath gives one -- and needs no code of its own: the loop over that
 *  chain does not run and its a[0] still contributes its share of the impurity energy. Note
 *  this cannot be handed to TriDiagonalMatrixToImpurityHamiltonian when one chain is empty,
 *  the way the natural orbital form can, because that builds a[1] to a[N] and treats a[0]
 *  and b[0] as the offset and the scaling: it would drop the impurity site and lose the
 *  other chain's share of its energy.
 *
 *  With AddSpin each orbital becomes two neighbouring fermions and the shift is not
 *  doubled: orbital n spin s is fermion FermionShift + 2*n + s.
 */
int DoubleTriDiagonalMatrixToImpurityHamiltonian(OperatorType * O, DoubleTriDiagonalMatrixType * M, bool AddSpin, uint32_t FermionShift, uint32_t NF)
{
  uint32_t i,s,ns,N,ip,iq;
  uint16_t oadd[2];
  N  = 1 + M->val.N + M->con.N;
  ns = (AddSpin==true) ? 2 : 1;
  InitNamedOperator(O, "Double Tridiagonal Impurity Hamiltonian");
  O->NBoson = 0;
  O->NFermion = FermionShift + ns*N;
  if(NF!=0)
  {
    if(NF<O->NFermion)
    {
      printf("ERROR in DoubleTriDiagonalMatrixToImpurityHamiltonian. NF (%d) < NFermion (%d)\n",NF,O->NFermion);
      return 1;
    }
    else
    {
      O->NFermion = NF;
    }
  }
  O->NWord = Min(O->NFermion, 1 + (O->NFermion - 1) / 8);
  InitOperatorFixedLength(O, 2, 8);
  for(s=0;s<ns;s++)
  {
    /*  the impurity site, its energy shared by the two chains */
    ip = FermionShift + s;
    oadd[1] = CrAnMask + ip;
    oadd[0] = ip;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, M->val.a[0] + M->con.a[0], &QDetNotOrderedFcFaBcBa);
    /*  the valence chain, site i at orbital i and bonded to orbital i-1 by b[i-1] */
    for(i=1;i<=M->val.N;i++)
    {
      iq = FermionShift + ns*i + s;
      ip = FermionShift + ns*(i-1) + s;
      oadd[1] = CrAnMask + iq;
      oadd[0] = iq;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, M->val.a[i], &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + ip;
      oadd[0] = iq;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, M->val.b[i-1], &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + iq;
      oadd[0] = ip;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, M->val.b[i-1], &QDetNotOrderedFcFaBcBa);
    }
    /*  and the conduction chain, site i at orbital val.N+i, the impurity being the site
     *  before con_1 */
    for(i=1;i<=M->con.N;i++)
    {
      iq = FermionShift + ns*(M->val.N+i) + s;
      if(i==1) ip = FermionShift + s;
      else     ip = FermionShift + ns*(M->val.N+i-1) + s;
      oadd[1] = CrAnMask + iq;
      oadd[0] = iq;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, M->con.a[i], &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + ip;
      oadd[0] = iq;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, M->con.b[i-1], &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + iq;
      oadd[0] = ip;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, M->con.b[i-1], &QDetNotOrderedFcFaBcBa);
    }
  }
  return 0;
}

/*  OLD - use TriDiagonalMatrixToImpurityHamiltonian instead
 *
 * Creates an anderson impurity Hamiltonian
 *
 *  The impurity is at index i=0 for spin up and i=1 for spin down
 *  The spin up bath is from i=2       to i=Mup.N
 *  The spin dn bath is from i=Mup.N+1 to i=Mup.N+Mdn.N-1
 */
uint32_t TriDiagonalMatrixToAndersonImpurityOperator(OperatorType * O, TriDiagonalMatrixType Mup, TriDiagonalMatrixType Mdn, double U, double mu)
{
  uint32_t i;
  uint16_t oadd[4];
  InitNamedOperator(O, "Single Site Anderson Impurity Hamiltonian");
  O->NFermion = Mup.N + Mdn.N;
  O->NWord = Min(O->NFermion, 1 + (O->NFermion - 1) / 8);
  InitOperatorFixedLength(O, 2, 8);
  InitOperatorFixedLength(O, 4, 8);
  oadd[3] = CrAnMask + 0;
  oadd[2] = CrAnMask + 1;
  oadd[1] = 1;
  oadd[0] = 0;
  RealOperatorAddLaderNormalOrder(O, oadd, 4, U, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 0;
  oadd[0] = 0;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, mu, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 1;
  oadd[0] = 1;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, mu, &QDetNotOrderedFcFaBcBa);
  /* spin up */
  oadd[1] = CrAnMask + 0;
  oadd[0] = 0;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.a[1], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 0;
  oadd[0] = 2;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.b[1], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 2;
  oadd[0] = 0;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.b[1], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 2;
  oadd[0] = 2;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.a[2], &QDetNotOrderedFcFaBcBa);
  for (i = 3; i <= Mup.N; i++)
  {
    oadd[1] = CrAnMask + i - 1;
    oadd[0] = i;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.b[i - 1], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i;
    oadd[0] = i - 1;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.b[i - 1], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i;
    oadd[0] = i;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.a[i], &QDetNotOrderedFcFaBcBa);
  }
  /* spin dn */
  oadd[1] = CrAnMask + 1;
  oadd[0] = 1;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.a[1], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 1;
  oadd[0] = Mup.N + 1;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.b[1], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + Mup.N + 1;
  oadd[0] = 1;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.b[1], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + Mup.N + 1;
  oadd[0] = Mup.N + 1;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.a[2], &QDetNotOrderedFcFaBcBa);
  for (i = 2; i < Mdn.N; i++)
  {
    oadd[1] = CrAnMask + Mup.N + i - 1;
    oadd[0] = Mup.N + i;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.b[i], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + Mup.N + i;
    oadd[0] = Mup.N + i - 1;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.b[i], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + Mup.N + i;
    oadd[0] = Mup.N + i;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.a[i + 1], &QDetNotOrderedFcFaBcBa);
  }
  return 0;
}

/*  OLD - use AndersonMatrixToImpurityHamiltonian instead
 *
 * Creates an anderson impurity Hamiltonian
 *
 *  The impurity is at index i=0 for spin up and i=1 for spin down
 *  The spin up bath is from i=2       to i=Mup.N
 *  The spin dn bath is from i=Mup.N+1 to i=Mup.N+Mdn.N-1
 */
uint32_t AndersonMatrixToAndersonImpurityOperator(OperatorType * O, AndersonMatrixType Mup, AndersonMatrixType Mdn, double U, double mu)
{
  uint32_t i;
  uint16_t oadd[4];
  InitNamedOperator(O, "Single Site Anderson Impurity Hamiltonian");
  O->NFermion = Mup.N + Mdn.N;
  O->NWord = Min(O->NFermion, 1 + (O->NFermion - 1) / 8);
  InitOperatorFixedLength(O, 2, 8);
  InitOperatorFixedLength(O, 4, 8);
  oadd[3] = CrAnMask + 0;
  oadd[2] = CrAnMask + 1;
  oadd[1] = 1;
  oadd[0] = 0;
  RealOperatorAddLaderNormalOrder(O, oadd, 4, U, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 0;
  oadd[0] = 0;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, mu, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 1;
  oadd[0] = 1;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, mu, &QDetNotOrderedFcFaBcBa);
  /* spin up */
  oadd[1] = CrAnMask + 0;
  oadd[0] = 0;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.a[1], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 0;
  oadd[0] = 2;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.b[1], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 2;
  oadd[0] = 0;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.b[1], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 2;
  oadd[0] = 2;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.a[2], &QDetNotOrderedFcFaBcBa);
  for (i = 3; i <= Mup.N; i++)
  {
    oadd[1] = CrAnMask + 0;
    oadd[0] = i;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.b[i - 1], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i;
    oadd[0] = 0;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.b[i - 1], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i;
    oadd[0] = i;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.a[i], &QDetNotOrderedFcFaBcBa);
  }
  /* spin dn */
  oadd[1] = CrAnMask + 1;
  oadd[0] = 1;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.a[1], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 1;
  oadd[0] = Mup.N + 1;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.b[1], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + Mup.N + 1;
  oadd[0] = 1;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.b[1], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + Mup.N + 1;
  oadd[0] = Mup.N + 1;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.a[2], &QDetNotOrderedFcFaBcBa);
  for (i = 2; i < Mdn.N; i++)
  {
    oadd[1] = CrAnMask + 1;
    oadd[0] = Mup.N + i;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.b[i], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + Mup.N + i;
    oadd[0] = 1;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.b[i], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + Mup.N + i;
    oadd[0] = Mup.N + i;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.a[i + 1], &QDetNotOrderedFcFaBcBa);
  }
  return 0;
}


/*  Creates an anderson impurity Hamiltonian
 *
 *  The impurity is at index i=0 for spin up and i=1 for spin down
 *  The first bath site is at index i=2 for spin up and i=3 for spin down
 *  The spin up bath is from i=4                     to i=Mup.val.N+Mup.con.N+3
 *  The spin dn bath is from i=Mup.val.N+Mup.con.N+4 to i=Mup.val.N+Mup.con.N+Mdn.val.N+Mdn.con.N+3
 */
uint32_t SecondDoubleTriDiagonalMatrixToAndersonImpurityOperator(OperatorType * O,SecondDoubleTriDiagonalMatrixType Mup,SecondDoubleTriDiagonalMatrixType Mdn, double U,double mu)
{
  uint32_t i;
  uint16_t oadd[4];
  InitNamedOperator(O, "Single Site Anderson Impurity Hamiltonian");
  O->NFermion = Mup.Dt.val.N + Mup.Dt.con.N + Mdn.Dt.val.N + Mdn.Dt.con.N + 4;
  O->NWord = Min(O->NFermion, 1 + (O->NFermion - 1) / 8);
  InitOperatorFixedLength(O, 2, 8);
  InitOperatorFixedLength(O, 4, 8);
  /* impurity */
  oadd[3] = CrAnMask + 0;
  oadd[2] = CrAnMask + 1;
  oadd[1] = 1;
  oadd[0] = 0;
  RealOperatorAddLaderNormalOrder(O, oadd, 4, U, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 0;
  oadd[0] = 0;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, mu, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 1;
  oadd[0] = 1;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, mu, &QDetNotOrderedFcFaBcBa);
  /* first bath site spin up */
  oadd[1] = CrAnMask + 0;
  oadd[0] = 0;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.Dt.a0, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 0;
  oadd[0] = 2;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.Dt.b0, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 2;
  oadd[0] = 0;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.Dt.b0, &QDetNotOrderedFcFaBcBa);
  /* first bath site spin dn */
  oadd[1] = CrAnMask + 1;
  oadd[0] = 1;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.Dt.a0, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 1;
  oadd[0] = 3;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.Dt.b0, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 3;
  oadd[0] = 1;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.Dt.b0, &QDetNotOrderedFcFaBcBa);
  /* spin up valence */
  oadd[1] = CrAnMask + 2;
  oadd[0] = 2;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.Dt.val.a[0], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 2;
  oadd[0] = 4;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.Dt.val.b[0], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 4;
  oadd[0] = 2;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.Dt.val.b[0], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 4;
  oadd[0] = 4;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.Dt.val.a[1], &QDetNotOrderedFcFaBcBa);
  for (i = 2; i <= Mup.Dt.val.N; i++)
  {
    oadd[1] = CrAnMask + i + 2;
    oadd[0] = i + 3;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.Dt.val.b[i - 1], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i + 3;
    oadd[0] = i + 2;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.Dt.val.b[i - 1], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i + 3;
    oadd[0] = i + 3;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.Dt.val.a[i], &QDetNotOrderedFcFaBcBa);
  }
  /* spin up conductance */
  oadd[1] = CrAnMask + 2;
  oadd[0] = 2;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.Dt.con.a[0], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 2;
  oadd[0] = Mup.Dt.val.N + 4;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.Dt.con.b[0], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + Mup.Dt.val.N + 4;
  oadd[0] = 2;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.Dt.con.b[0], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + Mup.Dt.val.N + 4;
  oadd[0] = Mup.Dt.val.N + 4;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.Dt.con.a[1], &QDetNotOrderedFcFaBcBa);
  for (i = 2; i <= Mup.Dt.con.N; i++)
  {
    oadd[1] = CrAnMask + i + 2 + Mup.Dt.val.N;
    oadd[0] = i + 3 + Mup.Dt.val.N;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.Dt.con.b[i - 1], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i + 3 + Mup.Dt.val.N;
    oadd[0] = i + 2 + Mup.Dt.val.N;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.Dt.con.b[i - 1], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i + 3 + Mup.Dt.val.N;
    oadd[0] = i + 3 + Mup.Dt.val.N;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.Dt.con.a[i], &QDetNotOrderedFcFaBcBa);
  }
  /* spin dn valence */
  oadd[1] = CrAnMask + 3;
  oadd[0] = 3;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.Dt.val.a[0], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 3;
  oadd[0] = Mup.Dt.val.N + Mup.Dt.con.N + 4;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.Dt.val.b[0], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + Mup.Dt.val.N + Mup.Dt.con.N + 4;
  oadd[0] = 3;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.Dt.val.b[0], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + Mup.Dt.val.N + Mup.Dt.con.N + 4;
  oadd[0] = Mup.Dt.val.N + Mup.Dt.con.N + 4;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.Dt.val.a[1], &QDetNotOrderedFcFaBcBa);
  for (i = 2; i <= Mdn.Dt.val.N; i++)
  {
    oadd[1] = CrAnMask + i + Mup.Dt.val.N + Mup.Dt.con.N + 2;
    oadd[0] = i + Mup.Dt.val.N + Mup.Dt.con.N + 3;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.Dt.val.b[i - 1], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i + Mup.Dt.val.N + Mup.Dt.con.N + 3;
    oadd[0] = i + Mup.Dt.val.N + Mup.Dt.con.N + 2;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.Dt.val.b[i - 1], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i + Mup.Dt.val.N + Mup.Dt.con.N + 3;
    oadd[0] = i + Mup.Dt.val.N + Mup.Dt.con.N + 3;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.Dt.val.a[i], &QDetNotOrderedFcFaBcBa);
  }
  /* spin dn conductance */
  oadd[1] = CrAnMask + 3;
  oadd[0] = 3;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.Dt.con.a[0], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 3;
  oadd[0] = Mup.Dt.val.N + Mup.Dt.con.N + Mdn.Dt.val.N + 4;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.Dt.con.b[0], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + Mup.Dt.val.N + Mup.Dt.con.N + Mdn.Dt.val.N + 4;
  oadd[0] = 3;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.Dt.con.b[0], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + Mup.Dt.val.N + Mup.Dt.con.N + Mdn.Dt.val.N + 4;
  oadd[0] = Mup.Dt.val.N + Mup.Dt.con.N + Mdn.Dt.val.N + 4;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.Dt.con.a[1], &QDetNotOrderedFcFaBcBa);
  for (i = 2; i <= Mdn.Dt.con.N; i++)
  {
    oadd[1] = CrAnMask + i + Mup.Dt.val.N + Mup.Dt.con.N + Mdn.Dt.val.N + 2;
    oadd[0] = i + Mup.Dt.val.N + Mup.Dt.con.N + Mdn.Dt.val.N + 3;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.Dt.con.b[i - 1], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i + Mup.Dt.val.N + Mup.Dt.con.N + Mdn.Dt.val.N + 3;
    oadd[0] = i + Mup.Dt.val.N + Mup.Dt.con.N + Mdn.Dt.val.N + 2;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.Dt.con.b[i - 1], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i + Mup.Dt.val.N + Mup.Dt.con.N + Mdn.Dt.val.N + 3;
    oadd[0] = i + Mup.Dt.val.N + Mup.Dt.con.N + Mdn.Dt.val.N + 3;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.Dt.con.a[i], &QDetNotOrderedFcFaBcBa);
  }
  return 0;
}

/*  Creates an anderson impurity Hamiltonian
 *
 *  The impurity is at index i=0 for spin up and i=1 for spin down
 *  The spin up bath is from i=2                     to i=Mup.val.N+Mup.con.N+1
 *  The spin dn bath is from i=Mup.val.N+Mup.con.N+2 to i=Mup.val.N+Mup.con.N+Mdn.val.N+Mdn.con.N+1
 */
uint32_t DoubleTriDiagonalMatrixToAndersonImpurityOperator(OperatorType * O, DoubleTriDiagonalMatrixType Mup,
                                                           DoubleTriDiagonalMatrixType Mdn, double U, double mu)
{
  uint32_t i;
  uint16_t oadd[4];
  InitNamedOperator(O, "Single Site Anderson Impurity Hamiltonian");
  O->NFermion = Mup.val.N + Mup.con.N + Mdn.val.N + Mdn.con.N + 2;
  O->NWord = Min(O->NFermion, 1 + (O->NFermion - 1) / 8);
  InitOperatorFixedLength(O, 2, 8);
  InitOperatorFixedLength(O, 4, 8);
  oadd[3] = CrAnMask + 0;
  oadd[2] = CrAnMask + 1;
  oadd[1] = 1;
  oadd[0] = 0;
  RealOperatorAddLaderNormalOrder(O, oadd, 4, U, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 0;
  oadd[0] = 0;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, mu, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 1;
  oadd[0] = 1;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, mu, &QDetNotOrderedFcFaBcBa);
  /* spin up valence */
  oadd[1] = CrAnMask + 0;
  oadd[0] = 0;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.val.a[0], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 0;
  oadd[0] = 2;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.val.b[0], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 2;
  oadd[0] = 0;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.val.b[0], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 2;
  oadd[0] = 2;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.val.a[1], &QDetNotOrderedFcFaBcBa);
  for (i = 2; i <= Mup.val.N; i++)
  {
    oadd[1] = CrAnMask + i;
    oadd[0] = i + 1;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.val.b[i - 1], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i + 1;
    oadd[0] = i;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.val.b[i - 1], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i + 1;
    oadd[0] = i + 1;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.val.a[i], &QDetNotOrderedFcFaBcBa);
  }
  /* spin up conductance */
  oadd[1] = CrAnMask + 0;
  oadd[0] = 0;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.con.a[0], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 0;
  oadd[0] = Mup.val.N + 2;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.con.b[0], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + Mup.val.N + 2;
  oadd[0] = 0;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.con.b[0], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + Mup.val.N + 2;
  oadd[0] = Mup.val.N + 2;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.con.a[1], &QDetNotOrderedFcFaBcBa);
  for (i = 2; i <= Mup.con.N; i++)
  {
    oadd[1] = CrAnMask + i + Mup.val.N;
    oadd[0] = i + 1 + Mup.val.N;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.con.b[i - 1], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i + 1 + Mup.val.N;
    oadd[0] = i + Mup.val.N;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.con.b[i - 1], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i + 1 + Mup.val.N;
    oadd[0] = i + 1 + Mup.val.N;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.con.a[i], &QDetNotOrderedFcFaBcBa);
  }
  /* spin dn valence */
  oadd[1] = CrAnMask + 1;
  oadd[0] = 1;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.val.a[0], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 1;
  oadd[0] = Mup.val.N + Mup.con.N + 2;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.val.b[0], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + Mup.val.N + Mup.con.N + 2;
  oadd[0] = 1;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.val.b[0], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + Mup.val.N + Mup.con.N + 2;
  oadd[0] = Mup.val.N + Mup.con.N + 2;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.val.a[1], &QDetNotOrderedFcFaBcBa);
  for (i = 2; i <= Mdn.val.N; i++)
  {
    oadd[1] = CrAnMask + i + Mup.val.N + Mup.con.N;
    oadd[0] = i + Mup.val.N + Mup.con.N + 1;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.val.b[i - 1], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i + Mup.val.N + Mup.con.N + 1;
    oadd[0] = i + Mup.val.N + Mup.con.N;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.val.b[i - 1], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i + Mup.val.N + Mup.con.N + 1;
    oadd[0] = i + Mup.val.N + Mup.con.N + 1;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.val.a[i], &QDetNotOrderedFcFaBcBa);
  }
  /* spin dn conductance */
  oadd[1] = CrAnMask + 1;
  oadd[0] = 1;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.con.a[0], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 1;
  oadd[0] = Mup.val.N + Mup.con.N + Mdn.val.N + 2;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.con.b[0], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + Mup.val.N + Mup.con.N + Mdn.val.N + 2;
  oadd[0] = 1;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.con.b[0], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + Mup.val.N + Mup.con.N + Mdn.val.N + 2;
  oadd[0] = Mup.val.N + Mup.con.N + Mdn.val.N + 2;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.con.a[1], &QDetNotOrderedFcFaBcBa);
  for (i = 2; i <= Mdn.con.N; i++)
  {
    oadd[1] = CrAnMask + i + Mup.val.N + Mup.con.N + Mdn.val.N;
    oadd[0] = i + Mup.val.N + Mup.con.N + Mdn.val.N + 1;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.con.b[i - 1], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i + Mup.val.N + Mup.con.N + Mdn.val.N + 1;
    oadd[0] = i + Mup.val.N + Mup.con.N + Mdn.val.N;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.con.b[i - 1], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i + Mup.val.N + Mup.con.N + Mdn.val.N + 1;
    oadd[0] = i + Mup.val.N + Mup.con.N + Mdn.val.N + 1;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.con.a[i], &QDetNotOrderedFcFaBcBa);
  }
  return 0;
}



/*  Creates an anderson impurity Hamiltonian
 *
 *  The impurity is at index i=0 for spin up and i=1 for spin down
 *  The first bath site is at index i=2 for spin up and i=3 for spin down
 *  The spin up bath is from i=4                     to i=Mup.val.N+Mup.con.N+3
 *  The spin dn bath is from i=Mup.val.N+Mup.con.N+4 to i=Mup.val.N+Mup.con.N+Mdn.val.N+Mdn.con.N+3
 */
uint32_t BADoubleTriDiagonalMatrixToAndersonImpurityOperator(OperatorType * O, BADoubleTriDiagonalMatrixType Mup, BADoubleTriDiagonalMatrixType Mdn, double U, double mu)
{
  uint32_t i;
  uint16_t oadd[4];
  InitNamedOperator(O, "Single Site Anderson Impurity Hamiltonian");
  O->NFermion = Mup.val.N + Mup.con.N + Mdn.val.N + Mdn.con.N;
  O->NWord = Min(O->NFermion, 1 + (O->NFermion - 1) / 8);
  InitOperatorFixedLength(O, 2, 8);
  InitOperatorFixedLength(O, 4, 8);
  /* impurity */
  oadd[3] = CrAnMask + 0;
  oadd[2] = CrAnMask + 1;
  oadd[1] = 1;
  oadd[0] = 0;
  RealOperatorAddLaderNormalOrder(O, oadd, 4, U, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 0;
  oadd[0] = 0;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.acvB - mu, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 1;
  oadd[0] = 1;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.acvB - mu, &QDetNotOrderedFcFaBcBa);
  /* first bath site spin up */
  oadd[1] = CrAnMask + 2;
  oadd[0] = 2;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.acvA, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 0;
  oadd[0] = 2;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.acvBA, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 2;
  oadd[0] = 0;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.acvBA, &QDetNotOrderedFcFaBcBa);
  /* first bath site spin dn */
  oadd[1] = CrAnMask + 3;
  oadd[0] = 3;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.acvA, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 1;
  oadd[0] = 3;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.acvBA, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 3;
  oadd[0] = 1;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.acvBA, &QDetNotOrderedFcFaBcBa);
  /* spin up valence */
  if(Mup.val.N>1)
  {
    oadd[1] = CrAnMask + 0;
    oadd[0] = 4;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.bvB0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4;
    oadd[0] = 0;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.bvB0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 2;
    oadd[0] = 4;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.bvA0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4;
    oadd[0] = 2;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.bvA0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4;
    oadd[0] = 4;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.val.a[2], &QDetNotOrderedFcFaBcBa);
  }
  for (i = 2; i < Mup.val.N; i++)
  {
    oadd[1] = CrAnMask + i + 2;
    oadd[0] = i + 3;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.val.b[i], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i + 3;
    oadd[0] = i + 2;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.val.b[i], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i + 3;
    oadd[0] = i + 3;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.val.a[i + 1], &QDetNotOrderedFcFaBcBa);
  }
  /* spin up conductance */
  if(Mup.con.N>1)
  {
    oadd[1] = CrAnMask + 0;
    oadd[0] = 4 + Mup.val.N - 1;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.bcB0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4 + Mup.val.N - 1;
    oadd[0] = 0;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.bcB0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 2;
    oadd[0] = 4 + Mup.val.N - 1;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.bcA0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4 + Mup.val.N - 1;
    oadd[0] = 2;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.bcA0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4 + Mup.val.N - 1;
    oadd[0] = 4 + Mup.val.N - 1;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.con.a[2], &QDetNotOrderedFcFaBcBa);
  }
  for (i = 2; i < Mup.con.N; i++)
  {
    oadd[1] = CrAnMask + i + 2 + Mup.val.N - 1;
    oadd[0] = i + 3 + Mup.val.N - 1;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.con.b[i], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i + 3 + Mup.val.N - 1;
    oadd[0] = i + 2 + Mup.val.N - 1;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.con.b[i], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i + 3 + Mup.val.N - 1;
    oadd[0] = i + 3 + Mup.val.N - 1;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.con.a[i + 1], &QDetNotOrderedFcFaBcBa);
  }
  /* spin dn valence */
  if(Mdn.val.N>1)
  {
    oadd[1] = CrAnMask + 1;
    oadd[0] = 4 + Mup.val.N + Mup.con.N - 2;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.bvB0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4 + Mup.val.N + Mup.con.N - 2;
    oadd[0] = 1;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.bvB0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 3;
    oadd[0] = 4 + Mup.val.N + Mup.con.N - 2;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.bvA0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4 + Mup.val.N + Mup.con.N - 2;
    oadd[0] = 3;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.bvA0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4 + Mup.val.N + Mup.con.N - 2;
    oadd[0] = 4 + Mup.val.N + Mup.con.N - 2;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.val.a[2], &QDetNotOrderedFcFaBcBa);
  }
  for (i = 2; i < Mdn.val.N; i++)
  {
    oadd[1] = CrAnMask + i + 2 + Mup.val.N + Mup.con.N - 2;
    oadd[0] = i + 3 + Mup.val.N + Mup.con.N - 2;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.val.b[i], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i + 3 + Mup.val.N + Mup.con.N - 2;
    oadd[0] = i + 2 + Mup.val.N + Mup.con.N - 2;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.val.b[i], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i + 3 + Mup.val.N + Mup.con.N - 2;
    oadd[0] = i + 3 + Mup.val.N + Mup.con.N - 2;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.val.a[i + 1], &QDetNotOrderedFcFaBcBa);
  }
  /* spin dn conductance */
  if(Mdn.con.N>1)
  {
    oadd[1] = CrAnMask + 1;
    oadd[0] = 4 + Mup.val.N + Mup.con.N + Mdn.val.N - 3;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.bcB0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4 + Mup.val.N + Mup.con.N + Mdn.val.N - 3;
    oadd[0] = 1;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.bcB0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 3;
    oadd[0] = 4 + Mup.val.N + Mup.con.N + Mdn.val.N - 3;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.bcA0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4 + Mup.val.N + Mup.con.N + Mdn.val.N - 3;
    oadd[0] = 3;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.bcA0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4 + Mup.val.N + Mup.con.N + Mdn.val.N - 3;
    oadd[0] = 4 + Mup.val.N + Mup.con.N + Mdn.val.N - 3;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.con.a[2], &QDetNotOrderedFcFaBcBa);
  }
  for (i = 2; i < Mdn.con.N; i++)
  {
    oadd[1] = CrAnMask + i + 2 + Mup.val.N + Mup.con.N + Mdn.val.N - 3;
    oadd[0] = i + 3 + Mup.val.N + Mup.con.N + Mdn.val.N - 3;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.con.b[i], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i + 3 + Mup.val.N + Mup.con.N + Mdn.val.N - 3;
    oadd[0] = i + 2 + Mup.val.N + Mup.con.N + Mdn.val.N - 3;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.con.b[i], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i + 3 + Mup.val.N + Mup.con.N + Mdn.val.N - 3;
    oadd[0] = i + 3 + Mup.val.N + Mup.con.N + Mdn.val.N - 3;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.con.a[i + 1], &QDetNotOrderedFcFaBcBa);
  }
  return 0;
}




/*  transform an operator to a matrix on a basis of single slater determinants
 *  the basis is given by the determinants in psi, in order they apear in psi
 *  The prefactors in psi are ignored
 */
uint32_t OperatorToDeterminantMatrix(const OperatorType * const O, const WaveFunctionType * const psi, CompactMatrixType *M)
{
  uint32_t i, j, k, l, m, ii, ij, ik;
  double phase;
  uint8_t *Det;
  M->Ni=psi->N;
  M->Nj=psi->N;
  if((O->QComplex==0))
  {
    M->QComplex=0;
  }
  else
  {
    M->QComplex=1;
  }
  if(InitCompactMatrix(M)) return 1;
  strcpy(M->Name, O->Name);
  Det = (uint8_t*)malloc(psi->NWord * sizeof(uint8_t));
  if (Det == NULL) AllocFailed();
  for(i=1;i<=psi->N;i++)
  {
    ii = i >> NBitsBlok;
    ik = i & BlokMask;
    ij = psi->NWord * ik;
    for (l = 0; l <= O->MaxLength; l++)
    {
      if (O->l[l].QComplex==0)
      {
        if(M->QComplex==0)
        { // all real
          for (j = 0; j < O->l[l].N; j++)
          {
            for (k = 0; k < psi->NWord; k++) Det[k] = psi->Determinant[ii][ij + k];
            phase = OperatorLader(Det, &(O->l[l].o[l * j]), l, O->NFermion);
            if (phase)
            {
              m=WaveFunctionIndexDet(psi, Det);
              if(m!=0)
              {
                M->VRe[(m-1)*M->Ni+(i-1)] += phase * O->l[l].VRe[j];
              }
            }
          }
        }
        else
        { // M complex and O real
          for (j = 0; j < O->l[l].N; j++)
          {
            for (k = 0; k < psi->NWord; k++) Det[k] = psi->Determinant[ii][ij + k];
            phase = OperatorLader(Det, &(O->l[l].o[l * j]), l, O->NFermion);
            if (phase)
            {
              m=WaveFunctionIndexDet(psi, Det);
              if(m!=0)
              {
                M->V[(m-1)*M->Ni+(i-1)].Re += phase * O->l[l].VRe[j];
              }
            }
          }
        }
      }
      else
      { // O complex
        for (j = 0; j < O->l[l].N; j++)
        {
          for (k = 0; k < psi->NWord; k++) Det[k] = psi->Determinant[ii][ij + k];
          phase = OperatorLader(Det, &(O->l[l].o[l * j]), l, O->NFermion);
          if (phase)
          {
            m=WaveFunctionIndexDet(psi, Det);
            if(m!=0)
            {
              M->V[(m-1)*M->Ni+(i-1)].Re += phase * O->l[l].V[j].Re;
              M->V[(m-1)*M->Ni+(i-1)].Im += phase * O->l[l].V[j].Im;
            }
          }
        }
      }
    }
  }
  free(Det);
  return 0;
}

/*  transform an operator to a matrix on a basis of single slater determinants
 *  the basis is given by the determinants in psi, in order they apear in psi
 *  The prefactors in psi are ignored
 *  The matrix brings a state on the basis of psi_j to a state on the basis of psi_i
 *  M->Ni = psi_i.N and M->Nj = psi_j.N
 */
uint32_t OperatorToDeterminantMatrix2(const OperatorType * const O, const WaveFunctionType * const psi_i, const WaveFunctionType * const psi_j, CompactMatrixType *M)
{
  uint32_t i, j, k, l, m, ii, ij, ik;
  double phase;
  uint8_t *Det;
  M->Ni=psi_i->N;
  M->Nj=psi_j->N;
  if((O->QComplex==0))
  {
    M->QComplex=0;
  }
  else
  {
    M->QComplex=1;
  }
  if(InitCompactMatrix(M)) return 1;
  strcpy(M->Name, O->Name);
  Det = (uint8_t*)malloc(Max(psi_i->NWord,psi_j->NWord) * sizeof(uint8_t));
  if (Det == NULL) AllocFailed();
  for(i=1;i<=psi_j->N;i++)
  {
    ii = i >> NBitsBlok;
    ik = i & BlokMask;
    ij = psi_j->NWord * ik;
    for (l = 0; l <= O->MaxLength; l++)
    {
      if (O->l[l].QComplex==0)
      {
        if(M->QComplex==0)
        { // all real
          for (j = 0; j < O->l[l].N; j++)
          {
            for (k = 0; k < psi_j->NWord; k++) Det[k] = psi_j->Determinant[ii][ij + k];
            phase = OperatorLader(Det, &(O->l[l].o[l * j]), l, O->NFermion);
            if (phase)
            {
              m=WaveFunctionIndexDet(psi_i, Det);
              if(m!=0)
              {
                M->VRe[(m-1)*M->Nj+(i-1)] += phase * O->l[l].VRe[j];
              }
            }
          }
        }
        else
        { // M complex and O real
          for (j = 0; j < O->l[l].N; j++)
          {
            for (k = 0; k < psi_j->NWord; k++) Det[k] = psi_j->Determinant[ii][ij + k];
            phase = OperatorLader(Det, &(O->l[l].o[l * j]), l, O->NFermion);
            if (phase)
            {
              m=WaveFunctionIndexDet(psi_i, Det);
              if(m!=0)
              {
                M->V[(m-1)*M->Nj+(i-1)].Re += phase * O->l[l].VRe[j];
              }
            }
          }
        }
      }
      else
      { // O complex
        for (j = 0; j < O->l[l].N; j++)
        {
          for (k = 0; k < psi_j->NWord; k++) Det[k] = psi_j->Determinant[ii][ij + k];
          phase = OperatorLader(Det, &(O->l[l].o[l * j]), l, O->NFermion);
          if (phase)
          {
            m=WaveFunctionIndexDet(psi_i, Det);
            if(m!=0)
            {
              M->V[(m-1)*M->Nj+(i-1)].Re += phase * O->l[l].V[j].Re;
              M->V[(m-1)*M->Nj+(i-1)].Im += phase * O->l[l].V[j].Im;
            }
          }
        }
      }
    }
  }
  free(Det);
  return 0;
}
uint32_t OperatorToDeterminantMatrixWithPrefactorsInPsi(const OperatorType * const O, const WaveFunctionType * const psi, CompactMatrixType *M)
{
  uint32_t i, j, k, l, m, ii, ij, ik;
  double phase;
  uint8_t *Det;
  M->Ni=psi->N;
  M->Nj=psi->N;
  if((O->QComplex==0)&&(psi->QComplex==0))
  {
    M->QComplex=0;
  }
  else
  {
    M->QComplex=1;
  }
  InitCompactMatrix(M);
  strcpy(M->Name, O->Name);
  Det = (uint8_t*)malloc(psi->NWord * sizeof(uint8_t));
  if (Det == NULL) AllocFailed();
  if(psi->QComplex==0)
  {
    for(i=1;i<=psi->N;i++)
    {
      ii = i >> NBitsBlok;
      ik = i & BlokMask;
      ij = psi->NWord * ik;
      for (l = 0; l <= O->MaxLength; l++)
      {
        if (O->l[l].QComplex==0)
        {
          if(M->QComplex==0)
          { // all real
            for (j = 0; j < O->l[l].N; j++)
            {
              for (k = 0; k < psi->NWord; k++) Det[k] = psi->Determinant[ii][ij + k];
              phase = OperatorLader(Det, &(O->l[l].o[l * j]), l, O->NFermion);
              if (phase)
              {
                m=WaveFunctionIndexDet(psi, Det);
                if(m!=0)
                {
                  M->VRe[(m-1)*M->Ni+(i-1)] += phase * psi->ValueRe[BlokElement(m)] * O->l[l].VRe[j] * psi->ValueRe[ii][ik];
                }
              }
            }
          }
          else
          { // M complex psi and O real
            for (j = 0; j < O->l[l].N; j++)
            {
              for (k = 0; k < psi->NWord; k++) Det[k] = psi->Determinant[ii][ij + k];
              phase = OperatorLader(Det, &(O->l[l].o[l * j]), l, O->NFermion);
              if (phase)
              {
                m=WaveFunctionIndexDet(psi, Det);
                if(m!=0)
                {
                  M->V[(m-1)*M->Ni+(i-1)].Re += phase * psi->ValueRe[BlokElement(m)] * O->l[l].VRe[j] * psi->ValueRe[ii][ik];
                }
              }
            }
          }
        }
        else
        { // psi real O complex
          for (j = 0; j < O->l[l].N; j++)
          {
            for (k = 0; k < psi->NWord; k++) Det[k] = psi->Determinant[ii][ij + k];
            phase = OperatorLader(Det, &(O->l[l].o[l * j]), l, O->NFermion);
            if (phase)
            {
              m=WaveFunctionIndexDet(psi, Det);
              if(m!=0)
              {
                M->V[(m-1)*M->Ni+(i-1)].Re += phase * psi->ValueRe[BlokElement(m)] * O->l[l].V[j].Re * psi->ValueRe[ii][ik];
                M->V[(m-1)*M->Ni+(i-1)].Im += phase * psi->ValueRe[BlokElement(m)] * O->l[l].V[j].Im * psi->ValueRe[ii][ik];
              }
            }
          }
        }
      }
    }
  }
  else
  {
    for(i=1;i<=psi->N;i++)
    {
      ii = i >> NBitsBlok;
      ik = i & BlokMask;
      ij = psi->NWord * ik;
      for (l = 0; l <= O->MaxLength; l++)
      {
        if (O->l[l].QComplex==0)
        { // psi complex O real
          for (j = 0; j < O->l[l].N; j++)
          {
            for (k = 0; k < psi->NWord; k++) Det[k] = psi->Determinant[ii][ij + k];
            phase = OperatorLader(Det, &(O->l[l].o[l * j]), l, O->NFermion);
            if (phase)
            {
              m=WaveFunctionIndexDet(psi, Det);
              if(m!=0)
              {
                M->V[(m-1)*M->Ni+(i-1)].Re += phase * ( psi->ValueRe[BlokElement(m)] * O->l[l].VRe[j] * psi->ValueRe[ii][ik]
                                                       + psi->ValueIm[BlokElement(m)] * O->l[l].VRe[j] * psi->ValueIm[ii][ik] );
                M->V[(m-1)*M->Ni+(i-1)].Im += phase * ( psi->ValueRe[BlokElement(m)] * O->l[l].VRe[j] * psi->ValueIm[ii][ik]
                                                       - psi->ValueIm[BlokElement(m)] * O->l[l].VRe[j] * psi->ValueRe[ii][ik] );
              }
            }
          }
        }
        else
        { // psi complex O complex
          for (j = 0; j < O->l[l].N; j++)
          {
            for (k = 0; k < psi->NWord; k++) Det[k] = psi->Determinant[ii][ij + k];
            phase = OperatorLader(Det, &(O->l[l].o[l * j]), l, O->NFermion);
            if (phase)
            {
              m=WaveFunctionIndexDet(psi, Det);
              if(m!=0)
              {
                M->V[(m-1)*M->Ni+(i-1)].Re += phase * ( psi->ValueRe[BlokElement(m)] * O->l[l].V[j].Re * psi->ValueRe[ii][ik]
                                                       + psi->ValueIm[BlokElement(m)] * O->l[l].V[j].Re * psi->ValueIm[ii][ik]
                                                       - psi->ValueRe[BlokElement(m)] * O->l[l].V[j].Im * psi->ValueIm[ii][ik]
                                                       + psi->ValueIm[BlokElement(m)] * O->l[l].V[j].Im * psi->ValueRe[ii][ik]);
                M->V[(m-1)*M->Ni+(i-1)].Im += phase * ( psi->ValueRe[BlokElement(m)] * O->l[l].V[j].Re * psi->ValueIm[ii][ik]
                                                       - psi->ValueIm[BlokElement(m)] * O->l[l].V[j].Re * psi->ValueRe[ii][ik]
                                                       + psi->ValueRe[BlokElement(m)] * O->l[l].V[j].Im * psi->ValueRe[ii][ik]
                                                       + psi->ValueIm[BlokElement(m)] * O->l[l].V[j].Im * psi->ValueIm[ii][ik]);
              }
            }
          }
        }
      }
    }
  }
  free(Det);
  return 0;
}


