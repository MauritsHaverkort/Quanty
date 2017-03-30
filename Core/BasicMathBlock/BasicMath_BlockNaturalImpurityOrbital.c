/*
 *     ( A[0][0]=Ai    | A[0][1]=Bib   | A[0][2]=Biv   | 0             | 0           | A[0][5]=Bic   | 0             | 0           )
 *     ( A[1][0]=Bib^+ | A[1][1]=Ab    | A[1][2]=Bbv   | 0             | 0           | A[0][5]=Bbc   | 0             | 0           )
 *     ( A[2][0]=Biv^+ | A[2][1]=Bbv^+ | A[2][2]=av2   | A[2][3]=bv2   | 0           | 0             | 0             | 0           )
 *   A=( 0             | 0             | A[3][2]=bv2^+ | A[3][3]=av3   | A[3][4]=bv3 | 0             | 0             | 0           )
 *     ( 0             | 0             | 0             | A[4][3]=bv3^+ | A[4][4]=av4 | 0             | 0             | 0           )
 *     ( A[5][0]=Bic^+ | A[5][1]=Bbc^+ | 0             | 0             | 0           | A[5][5]=ac2   | A[5][6]=bc2   | 0           )
 *     ( 0             | 0             | 0             | 0             | 0           | A[6][5]=bc2^+ | A[6][6]=ac3   | A[6][7]=bc3 )
 *     ( 0             | 0             | 0             | 0             | 0           | 0             | A[7][6]=bc3^+ | A[7][7]=ac4 )
 *
 *  A = T Ap T^C
 *
 *     ( A[0][0]=av1   | 0             | A[0][2]=bv1   | 0             | 0           | 0             | 0             | 0           )
 *     ( 0             | A[1][1]=ac1   | 0             | 0             | 0           | A[0][5]=bc1   | 0             | 0           )
 *     ( A[2][0]=bv1^+ | 0             | A[2][2]=av2   | A[2][3]=bv2   | 0           | 0             | 0             | 0           )
 *  Ap=( 0             | 0             | A[3][2]=bv2^+ | A[3][3]=av3   | A[3][4]=bv3 | 0             | 0             | 0           )
 *     ( 0             | 0             | 0             | A[4][3]=bv3^+ | A[4][4]=av4 | 0             | 0             | 0           )
 *     ( 0             | A[5][1]=bc1^+ | 0             | 0             | 0           | A[5][5]=ac2   | A[5][6]=bc2   | 0           )
 *     ( 0             | 0             | 0             | 0             | 0           | A[6][5]=bc2^+ | A[6][6]=ac3   | A[6][7]=bc3 )
 *     ( 0             | 0             | 0             | 0             | 0           | 0             | A[7][6]=bc3^+ | A[7][7]=ac4 )
 *
 *
 *     ( A[0][0]=bv0   | A[0][1]=bc0   | 0             | 0             | 0           | 0             | 0             | 0           )
 *     ( A[1][0]=~bc0  | A[1][1]=~-bv0 | 0             | 0             | 0           | 0             | 0             | 0           )
 *     ( 0             | 0             | 1             | 0             | 0           | 0             | 0             | 0           )
 *   T=( 0             | 0             | 0             | 1             | 0           | 0             | 0             | 0           )
 *     ( 0             | 0             | 0             | 0             | 1           | 0             | 0             | 0           )
 *     ( 0             | 0             | 0             | 0             | 0           | 1             | 0             | 0           )
 *     ( 0             | 0             | 0             | 0             | 0           | 0             | 1             | 0           )
 *     ( 0             | 0             | 0             | 0             | 0           | 0             | 0             | 1           )
 *
 *  with,
 *   T T^C = 1
 *
 *   (A B)(i,j) = \sum_k A(i,k) B(k,j)
 *
 *   y[i] = \sum_j A(i,j) x[j]
 *
 *   element (i,j) is found at position [i*Nj + j]
 *
 *   Note that the block A[2;;4][2;;4] and the block A[5;;7][5;;7] are normal tri-diagonal matrices
 *   The states av1 and ac1 in Ap are linear combinations of the impurity site and a special linear combination of the other sites
 *   such that the impurity has an occupation n, the other site an occupation 1-n, av1 an occupation 1 and ac1 an occupaiton of 0.
 */

/*
 *   Initialize a struct of BlockListOfPoles type
 *   Init val and con as Tridiagonal matrices
 *   Init a0 as a   BlockSize by   BlockSize matrix
 *   Init T  as a 4*BlockSize by 4*BlockSize matrix
 *   QComplex, Blocksize, val.N and con.N need to be set
 */
uint32_t InitBlockNaturalImpurityOrbitalMatrix(BlockNaturalImpurityOrbitalMatrixType *const A)
{
  A->val.BlockSize=A->BlockSize;
  A->val.QComplex=A->QComplex;
  InitBlockTridiagonalMatrix(&(A->val));
  A->con.BlockSize=A->BlockSize;
  A->con.QComplex=A->QComplex;
  InitBlockTridiagonalMatrix(&(A->con));
  A->a0.QComplex=A->QComplex;
  A->a0.Ni=A->BlockSize;
  A->a0.Nj=A->BlockSize;
  InitCompactMatrix(&(A->a0));
  A->TMT.QComplex=A->QComplex;
  A->TMT.Ni=4*A->BlockSize;
  A->TMT.Nj=4*A->BlockSize;
  InitCompactIdentityMatrix(&(A->TMT));
  return 0;
}

uint32_t FreeBlockNaturalImpurityOrbitalMatrix(BlockNaturalImpurityOrbitalMatrixType *const A)
{
  FreeBlockTridiagonalMatrix(&(A->val));
  FreeBlockTridiagonalMatrix(&(A->con));
  FreeCompactMatrix(&(A->a0));
  FreeCompactMatrix(&(A->TMT));
  return 0;
}

uint32_t BlockNaturalImpurityOrbitalMatrixChop(BlockNaturalImpurityOrbitalMatrixType *const A, double const Epsilon)
{
  uint32_t i, KeepComplex=0;
  if(A->QComplex)
  {
    for(i=0;i<A->BlockSize*A->BlockSize*(A->val.N+1);i++)
    {
      if(Abs(A->val.a[i].Re)<Epsilon)A->val.a[i].Re=0;
      if(Abs(A->val.a[i].Im)<Epsilon)
      {
        A->val.a[i].Im=0;
      }
      else
      {
        KeepComplex=1;
      }
    }
    for(i=0;i<A->BlockSize*A->BlockSize*(A->val.N);i++)
    {
      if(Abs(A->val.b[i].Re)<Epsilon)A->val.b[i].Re=0;
      if(Abs(A->val.b[i].Im)<Epsilon)
      {
        A->val.b[i].Im=0;
      }
      else
      {
        KeepComplex=1;
      }
    }
    for(i=0;i<A->BlockSize*A->BlockSize*(A->con.N+1);i++)
    {
      if(Abs(A->con.a[i].Re)<Epsilon)A->con.a[i].Re=0;
      if(Abs(A->con.a[i].Im)<Epsilon)
      {
        A->con.a[i].Im=0;
      }
      else
      {
        KeepComplex=1;
      }
    }
    for(i=0;i<A->BlockSize*A->BlockSize*(A->con.N);i++)
    {
      if(Abs(A->con.b[i].Re)<Epsilon)A->con.b[i].Re=0;
      if(Abs(A->con.b[i].Im)<Epsilon)
      {
        A->con.b[i].Im=0;
      }
      else
      {
        KeepComplex=1;
      }
    }
    for(i=0;i<A->BlockSize*A->BlockSize;i++)
    {
      if(Abs(A->a0.V[i].Re)<Epsilon)A->a0.V[i].Re=0;
      if(Abs(A->a0.V[i].Im)<Epsilon)
      {
        A->a0.V[i].Im=0;
      }
      else
      {
        KeepComplex=1;
      }
    }
    for(i=0;i<4*A->BlockSize*A->BlockSize;i++)
    {
      if(Abs(A->TMT.V[i].Re)<Epsilon)A->TMT.V[i].Re=0;
      if(Abs(A->TMT.V[i].Im)<Epsilon)
      {
        A->TMT.V[i].Im=0;
      }
      else
      {
        KeepComplex=1;
      }
    }
    if(!KeepComplex)
    {
      ComplexBlockTridiagonalMatrixToReal(&A->val);
      ComplexBlockTridiagonalMatrixToReal(&A->con);
      MakeCompactComplexMatrixReal(&A->a0);
      MakeCompactComplexMatrixReal(&A->TMT);
      A->QComplex=0;
    }
  }
  else
  {
    for(i=0;i<A->BlockSize*A->BlockSize*(A->val.N+1);i++)
    {
      if(Abs(A->val.aRe[i])<Epsilon)A->val.aRe[i]=0;
    }
    for(i=0;i<A->BlockSize*A->BlockSize*(A->val.N);i++)
    {
      if(Abs(A->val.bRe[i])<Epsilon)A->val.bRe[i]=0;
    }
    for(i=0;i<A->BlockSize*A->BlockSize*(A->con.N+1);i++)
    {
      if(Abs(A->con.aRe[i])<Epsilon)A->con.aRe[i]=0;
    }
    for(i=0;i<A->BlockSize*A->BlockSize*(A->con.N);i++)
    {
      if(Abs(A->con.bRe[i])<Epsilon)A->con.bRe[i]=0;
    }
    for(i=0;i<A->BlockSize*A->BlockSize;i++)
    {
      if(Abs(A->a0.VRe[i])<Epsilon)A->a0.VRe[i]=0;
    }
    for(i=0;i<4*A->BlockSize*A->BlockSize;i++)
    {
      if(Abs(A->TMT.VRe[i])<Epsilon)A->TMT.VRe[i]=0;
    }
  }
  return 0;
}


/*
 *  B <- A
 *  B not initialized on entry
 */
uint32_t CopyBlockNaturalImpurityOrbitalMatrix(BlockNaturalImpurityOrbitalMatrixType A, BlockNaturalImpurityOrbitalMatrixType *B)
{
  //TODO: all
  return 0;
}


