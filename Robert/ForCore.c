/* Functions that can be merged into core once tested to Maurits' standards... */




/* Should go in BasicMath_CompactMatrixVector.c
 * 
 * M<-A
 */
uint32_t CompactMatrixCopyToMatrix(MatrixType * const M, CompactMatrixTypeOld const A)
{
  uint32_t i;
  M->Ni = A.Ni;
  M->Nj = A.Nj;
  M->QComplex = A.QComplex;
  if (InitMatrix(M)) return 1;
  strcpy(M->Name, A.Name);
  if (M->QComplex)
  {
    for (i = 0; i < M->Ni*M->Nj; i++)
    {
      M->VRe[0][i] = A.V[i].Re;
      M->VIm[0][i] = A.V[i].Im;
    }
  }
  else
  {
    for (i = 0; i < M->Ni*M->Nj; i++)
    {
      M->VRe[0][i] = A.VRe[i];
    }
  }
  return 0;
}


/* BasicMath_MatrixVector.c
 * 
 * Issue?:  I've found that I have to calculate T^{trans}AT to get the same results as the non-BLAS version which
 * calculates TAT^{trans}  (I haven't had a complex case to test the conjugation yet)
 * This must be an issue with the conversion between matrix types or something?
 */
uint32_t MatrixTransformBLAS(MatrixType * const M, MatrixType const A, MatrixType const T)
{
  CompactMatrixTypeOld AComp, TComp,TComp2,BComp;
  if ((A.Ni != A.Nj) || (A.Ni != T.Nj))
  {
    printf("Invallid dimensions for matrix transformation A.Ni should be T.Nj and A should be square\n");
    printf("A.Ni=%i A.Nj=%i T.Ni=%i T.Nj=%i \n", A.Ni, A.Nj, T.Ni, T.Nj);
    return 1;
  }
  
  printf("MatrixTransormBLAS\nA.Ni=%i A.Nj=%i T.Ni=%i T.Nj=%i \n", A.Ni, A.Nj, T.Ni, T.Nj);fflush(stdout);
  
  if(MatrixCopyToCompact(&AComp,A))
  {
    printf("Matrix conversion to compact failed in MatrixTransformBLAS\n");
    return 1;
  }
  if(MatrixCopyToCompact(&TComp,T))
  {
    printf("Matrix conversion to compact failed in MatrixTransformBLAS\n");
    return 1;
  }
  if(MatrixCopyToCompact(&TComp2,T))
  {
    printf("Matrix conversion to compact failed in MatrixTransformBLAS\n");
    return 1;
  }
  if(CompactMatrixMultiply3('C','N','N',&TComp,&AComp,&TComp2,&BComp))
  {
    printf("CompactMatrixMultiply3 failed in MatrixTransformBLAS\n");
    return 1;
  }
  if(CompactMatrixCopyToMatrix(M,BComp))
  {
    printf("CompactMatrixCopyToMatrix failed in MatrixTransformBLAS\n");
    return 1;
  }
    
  FreeCompactMatrix(&AComp);
  FreeCompactMatrix(&TComp);
  FreeCompactMatrix(&TComp2);
  FreeCompactMatrix(&BComp);
  
  return 0;
}



/* Should go in BasicMath_CompactMatrixVector.c
 * 
 * !!Not all cases thoroughly tested yet!!
 * 
 * Converts the square section enclosed by iStart and iEnd to a CompactMatrix
 * 
 * Useful for returning the bath part of a block tridiagonal matrix
 *
 */
uint32_t BlockTriDiagonalMatrixPartToCompactMatrix(BlockTriDiagonalMatrixType const * A, CompactMatrixTypeOld * const B, const uint32_t iStart, const uint32_t iEnd)
{
  int i,j,k,apos,bpos,bposL,bposR;
  B->QComplex=A->QComplex;
  B->Ni=(iEnd - iStart +1)* A->BlockSize;
  B->Nj=(iEnd - iStart +1) * A->BlockSize;
  if(InitCompactMatrix(B))
  {
    printf("InitCompactMatrix failed in BlockTriDiagonalMatrixToCompactMatrix\n");
    return 1;
  }
  bpos=0;
  apos=(1+iStart)*A->BlockSize*A->BlockSize;
  if(A->QComplex)
  {
    for(k=iStart;k<=iEnd;k++)
    {
      for(i=0;i<A->BlockSize;i++)
      {
        for(j=0;j<A->BlockSize;j++)
        {
          B->V[bpos].Re=A->a[apos].Re;
          B->V[bpos].Im=A->a[apos].Im;
          bpos++;
          apos++;
        }
        bpos+=B->Nj-A->BlockSize;
      }
      bpos+=A->BlockSize;
    }
    apos=(1+iStart)*A->BlockSize*A->BlockSize;
    bposL=A->BlockSize;
    bposR=A->BlockSize*B->Nj;
    for(k=iStart;k<iEnd;k++)
    {
      for(i=0;i<A->BlockSize;i++)
      {
        for(j=0;j<A->BlockSize;j++)
        {
          B->V[bposL].Re= A->b[apos].Re;
          B->V[bposR].Re= A->b[apos].Re;
          B->V[bposL].Im=-A->b[apos].Im;
          B->V[bposR].Im= A->b[apos].Im;
          bposR++;
          bposL+=B->Nj;
          apos++;
        }
        bposR+=B->Nj-A->BlockSize;
        bposL+=(1-A->BlockSize*B->Nj);
      }
      bposR+=A->BlockSize;
      bposL+=(A->BlockSize*B->Nj);
    }
  }
  else
  {
    for(k=iStart;k<=iEnd;k++)
    {
      for(i=0;i<A->BlockSize;i++)
      {
        for(j=0;j<A->BlockSize;j++)
        {
          B->VRe[bpos]=A->aRe[apos];
          bpos++;
          apos++;
        }
        bpos+=B->Nj-A->BlockSize;
      }
      bpos+=A->BlockSize;
    }
    apos=(1+iStart)*A->BlockSize*A->BlockSize;
    bposL=A->BlockSize;
    bposR=A->BlockSize*B->Nj;
    for(k=iStart;k<iEnd;k++)
    {
      for(i=0;i<A->BlockSize;i++)
      {
        for(j=0;j<A->BlockSize;j++)
        {
          B->VRe[bposL]=A->bRe[apos];
          B->VRe[bposR]=A->bRe[apos];
          bposR++;
          bposL+=B->Nj;
          apos++;
        }
        bposR+=B->Nj-A->BlockSize;
        bposL+=(1-A->BlockSize*B->Nj);
      }
      bposR+=A->BlockSize;
      bposL+=(A->BlockSize*B->Nj);
    }
  }
  return 0;
}

/*  Should go in BasicMath_BlockTypeConversion.c
 *
 * Takes a blocktridiagonal matrix as input
 * Diagonalizes the bath part of this matrix, and uses the eigenvectors to transform the
 * block coupling the impurity to the bath into a row of blocks coupling the impurity
 * to the new Anderson sites with energies obtained from eigenvalues.
 * 
 * !!Does not work as intended yet!!
 * Two issues:
 * 1. I get mixed hopping in the anderson matrix: xz hops to both xz/yz, and same for yz. 
 * 2. The bath orbitals become reordered because the eigenvalues from the diagonalization
 *    are sorted. Not sure if this is an issue.
 * 
 * End result is a spectrum calculated after using this is slightly different from the original
 * spectrum using the BlockTridiagonalMatrix, and I think it is due to issue #1.
 * 
 */
uint32_t BlockTriDiagonalMatrixToBlockAndersonMatrix(BlockTriDiagonalMatrixType const ATri, BlockAndersonMatrixType * const B)
{
  int i,j,bpos,cpos,epos,k;
  CompactMatrixTypeOld Afull, Afun, C,D,E;
  double *Aval;
  
  /* Convert the portion after the first block to compact matrix */
  if(BlockTriDiagonalMatrixPartToCompactMatrix(&ATri,&Afull,1,ATri.N-1))
  {
    printf("BlockTriDiagonalMatrixPartToCompactMatrix failed in BlockTriDiagonalMatrixToBlockAndersonMatrix\n");
    return 1;
  }
  
  /* Diagonalize compact matrix */
  Aval=malloc(ATri.BlockSize*(ATri.N-1)*sizeof(double));
  if(CompactMatrixEigensystem(&Afull, &Afun, Aval))
  {
    printf("CompactMatrixEigensystem failed in BlockTriDiagonalMatrixToBlockAndersonMatrix\n");
    return 1;
  }
  
  B->QComplex = ATri.QComplex;
  B->N = ATri.N;
  B->BlockSize = ATri.BlockSize;
  InitBlockAndersonMatrix(B);
  
  if(ATri.QComplex)
  {
    /* First diagonal block is same as ATri first diagonal block */
    bpos = B->BlockSize*B->BlockSize;
    for(i = 0; i < B->BlockSize; i++)
    {
      for(j = 0; j < B->BlockSize; j++)
      {
        B->a[bpos].Re = ATri.a[bpos].Re;
        bpos++;
      }
    }
    
    /* Set diagonals to eigenvalues */
    cpos = 0;
    for(k = 1; k < ATri.N; k++)
    {
      for(i = 0; i < B->BlockSize; i++)
      {
        B->a[(k+1)*B->BlockSize*B->BlockSize + i*B->BlockSize + i].Re = Aval[cpos++];
      }
    }
    
    /* Transform first block of ATri using eigenfunctions from diagonalization */
    C.QComplex = B->QComplex;
    C.Ni = B->BlockSize;
    C.Nj = B->BlockSize;
    InitCompactMatrix(&C);
    bpos = B->BlockSize*B->BlockSize;
    for(i = 0; i < B->BlockSize*B->BlockSize; i++)
    {
      C.V[i].Re = ATri.b[bpos++].Re;
      C.V[i].Im = ATri.b[bpos++].Im;
    }
    
    D.QComplex = B->QComplex;
    D.Nj = B->BlockSize;
    D.Ni = B->BlockSize*(B->N-1);
    InitCompactMatrix(&D);
    for(i = 0; i < B->BlockSize; i++)
    {
      for(j = 0; j < B->BlockSize*(B->N-1); j++)
      {
        D.V[i*B->BlockSize*(B->N-1)+ j].Re = Afun.V[i*B->BlockSize*(B->N-1)+j].Re; 
        D.V[i*B->BlockSize*(B->N-1) + j].Im = Afun.V[i*B->BlockSize*(B->N-1)+j].Im; 
      }
    }
    
    /* Transform the block that coupled impurity to chain into a row of blocks
       coupling impurity to anderson bath
     */
    CompactMatrixMultiply('N', 'N', &C, &D, &E);
    /* Now copy E into the beta values of block matrix*/
    bpos = B->BlockSize*B->BlockSize;
    for(k = 1; k < B->BlockSize; k++)
    {
      for(i = 0; i < B->BlockSize; i++)
      {
        for(j = 0; j < B->BlockSize; j++)
        {
          B->b[bpos].Re = E.V[j*B->BlockSize*(B->N-1)+(k-1)*B->BlockSize+i].Re;
          B->b[bpos++].Im = E.V[j*B->BlockSize*(B->N-1)+(k-1)*B->BlockSize+i].Im;
        }
      }
    }
  }
  else
  {
    /* Real */
    /* First diagonal block is same as ATri first diagonal block */
    bpos = B->BlockSize*B->BlockSize;
    for(i = 0; i < B->BlockSize; i++)
    {
      for(j = 0; j < B->BlockSize; j++)
      {
        B->aRe[bpos] = ATri.aRe[bpos];
        bpos++;
      }
    }
    
    /* Set diagonals to eigenvalues */
    cpos = 0;
    for(k = 1; k < ATri.N; k++)
    {
      for(i = 0; i < B->BlockSize; i++)
      {
        B->aRe[(k+1)*B->BlockSize*B->BlockSize + i*B->BlockSize + i] = Aval[cpos++];
      }
    }
    
    /* Transform first block of ATri using eigenfunctions from diagonalization */
    C.QComplex = B->QComplex;
    C.Ni = B->BlockSize;
    C.Nj = B->BlockSize;
    InitCompactMatrix(&C);
    bpos = B->BlockSize*B->BlockSize;
    for(i = 0; i < B->BlockSize*B->BlockSize; i++)
    {
      C.VRe[i] = ATri.bRe[bpos++];
    }
    
    D.QComplex = B->QComplex;
    D.Nj = B->BlockSize;
    D.Ni = B->BlockSize*(B->N-1);
    InitCompactMatrix(&D);
    for(i = 0; i < B->BlockSize; i++)
    {
      for(j = 0; j < B->BlockSize*(B->N-1); j++)
      {
        D.VRe[i*B->BlockSize*(B->N-1) + j] = Afun.VRe[i*B->BlockSize*(B->N-1)+j]; 
       }
    }
    
    /* Transform the block that coupled impurity to chain into a row of blocks
       coupling impurity to anderson bath */
    CompactMatrixMultiply('N', 'N', &C, &D, &E);
    
    /* Now copy E into the beta values of anderson block matrix*/
    bpos = B->BlockSize*B->BlockSize;
    for(k = 1; k < B->N; k++)
    {
      for(i = 0; i < B->BlockSize; i++)
      {
        for(j = 0; j < B->BlockSize; j++)
        {
          B->bRe[bpos++] = E.VRe[j*B->BlockSize*(B->N-1)+(k-1)*B->BlockSize+i];
        }
      }
    }
  }
  
  FreeCompactMatrix(&C);
  FreeCompactMatrix(&D);
  FreeCompactMatrix(&E);
  FreeCompactMatrix(&Afun);
  FreeCompactMatrix(&Afull);
  free(Aval);
  
  return 0;
}


/*  BasicMath_BlockTypeConversion.c
 *
 * Converts a BlockListofpoles to blockandersonmatrix, using a blocktridiagonalmatrix as intermediate
 * (I was not sure how to do the conversion directly)
 * 
 * 
 */
uint32_t BlockListOfPolesToBlockAndersonMatrix(BlockListOfPoles const A, BlockAndersonMatrixType *const B, int const NBlockMax)
{
  BlockTriDiagonalMatrixType ATri;
  
  if(BlockListOfPolesToBlockTridiagonalMatrix(A, &ATri, NBlockMax, DBL_EPSILON))
  {
    printf("BlockListOfPolesToBlockTridiagonalMatrix failed in BlockListOfPolesToBlockAndersonMatrix\n");
    return 1;
  }
  
  BlockTridiagonalMatrixChop(&ATri,10*DBL_EPSILON);
  
  if(BlockTriDiagonalMatrixToBlockAndersonMatrix(ATri,B))
  {
    printf("BlockTriDiagonalMatrixToBlockAndersonMatrix failed in BlockListOfPolesToBlockAndersonMatrix\n");
    return 1;
  }

  return 0;
}

/*  on entry A contains the Green's function in a BlockListOfPoles format.
 *
 * Converts a BlockListofpoles to blockandersonmatrix, using a blocktridiagonalmatrix as intermediate
 * (I was not sure how to do the conversion directly)
 * 
 * 
 *  AT specifies how to reduce the dimension
 *
 */
uint32_t BlockListOfPolesToBlockAndersonMatrixReduceDimension(BlockListOfPoles const A, BlockAndersonMatrixType *B, CompactMatrixTypeOld const *AT, int const NBlockMax)
{
  BlockTriDiagonalMatrixType ATri;
  
  if(BlockListOfPolesToBlockTridiagonalMatrixReduceDimension(A, &ATri, AT, NBlockMax, DBL_EPSILON))
  {
    printf("BlockListOfPolesToBlockTridiagonalMatrixReduceDimension failed in BlockListOfPolesToBlockAndersonMatrixReduceDimension\n");
    return 1;
  }
  
  BlockTridiagonalMatrixChop(&ATri,1000*DBL_EPSILON);
  
  if(BlockTriDiagonalMatrixToBlockAndersonMatrix(ATri,B))
  {
    printf("BlockTriDiagonalMatrixToBlockAndersonMatrix failed in BlockListOfPolesToBlockAndersonMatrixReduceDimension\n");
    return 1;
  }
  
  return 0;
}







/* BasicMath_BlockMatrixToOperator.c 
 * 
 * !!Complex case is not quite completed. Everything is there, but gPert is TriDiagonalMatrixType which does
 *   not accept complex values.  So I always just give it the real part of the b value
 * 
 * 
 * Starting from a BlockTridiagonalMatrix creates an impurity Hamiltonian and represents the bath as
 *  perturbation operators and tridiagonal matrices. NCell specifies how many sites to include in the
 * "impurity" part (the rest will be in the perturbation). NCell=1 is just the actual impurity
 * 
 *  H on input not initiallized
 *  NBlockMax the maximal number of blocks in the Hamiltonian, this includes the impurity.
 *   NBlockMax=2 would thus correspond to a triditional ligand field calculation
 *  if Spin=0 each block in A is taken once into the Hamiltonian
 *  If Spin =1 (!0) each block is copied twice, once for spin up, once for spin down. The order is block up, block down
 *  NFermionShift is the numbers of Fermions shifted. (i.e. if you need a 2p core before a d impurity block NFermionShift=6 and the p-shell is from 0 to 5, the d-shell form 6 to 15.
 *  NFermion is the total number of Fermions in the basis. 
 *  If NFermion=0 on input the number will be calculated such that H fits.
 */
uint32_t BlockTridiagonalMatrixToImpurityHamiltonianPerturbation(OperatorType *const H, OperatorType * const VPert, TriDiagonalMatrixType * const gPert, BlockTriDiagonalMatrixType const *A, uint32_t NCell, uint32_t const NBlockMax, uint32_t const Spin, uint32_t const NFermionShift, uint32_t const NFermion)
{
  uint32_t i,j,k,ia,ifer,jfer,ib;
  uint16_t oadd[2];
  
  
  if(!(NCell < NBlockMax)||!(NCell < A->N)||(NBlockMax<2)||(A->N<2))
  {
    printf("Not enough bath sites for perturbation method\n");
    return 1;
  }
  
  InitNamedOperator(H, "Impurity Hamiltonian");
  H->NBoson=0;
  if(Spin==0)
  {
    H->NFermion=NFermionShift+NCell*A->BlockSize;
  }
  else
  {
    H->NFermion=NFermionShift+2*NCell*A->BlockSize;
  }
  if(NFermion!=0 )
  {
    if(NFermion>=H->NFermion)
    {
      H->NFermion=NFermion;
    }
    else
    {
      printf("Warning in BlockTridiagonalMatrixToImpurityHamiltonian\n");
      printf("Function called with NFermion=%i, but total number of Fermions must be larger than %i\n",NFermion,H->NFermion);
      printf("Number of Fermions enlarged\n");
    }
  }
  
  InitOperatorFixedLength(H, 2, 8);
  
  if(Spin)
    k = 2;
  else
    k = 1;
  
  /* Initialize perturbation operators */
  for(i = 0; i < k*A->BlockSize; i++)
  {
    InitOperator(&VPert[i]);
    VPert[i].NFermion = H->NFermion;
    VPert[i].NBoson = 0;
    VPert[i].NWord = Min(VPert[i].NFermion, 1 + (VPert[i].NFermion - 1) / 8);
    VPert[i].QComplex = 0;
    InitOperatorFixedLength(&VPert[i],1,8);
  }
  
  /* Initialize Perturbation tridiagonal matrices */
  for(i = 0; i < k*A->BlockSize; i++)
  {
    gPert[i].N = NBlockMax - NCell;
    InitTriDiagonalMatrix(&gPert[i]);
  }

  
  
  
  if(A->QComplex)
  {
    /* !!Complex not completed: TridiagonalMatrix doesn't have a complex option, but it seems the b values would
     * be complex in this case? */
    
    MakeRealOperatorFixedLengthComplex(H, 2);
    for(i = 0; i < k*A->BlockSize; i++)
    {
      MakeRealOperatorFixedLengthComplex(&VPert[i],1);

    }
    
    if(Spin==0)
    {
      /* a values first*/
      ia=A->BlockSize*A->BlockSize;
      for(k=0;k<Min(NBlockMax,A->N);k++)
      {
        if(k < NCell)
        {
          /* Part up to NCell goes in H */
          ifer=NFermionShift+k*A->BlockSize;
          for(i=0;i<A->BlockSize;i++)
          {
            jfer=NFermionShift+k*A->BlockSize;
            for(j=0;j<A->BlockSize;j++)
            {
              oadd[1]=CrAnMask+ifer;
              oadd[0]=jfer;
              ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->a[ia].Re, A->a[ia].Im, &QDetNotOrderedFcFaBcBa);
              ia++;
              jfer++;
            }
            ifer++;
          }
        }
        else
        {
          /* Blocks after NCell go in gPert*/
          for(i = 0; i < A->BlockSize; i++)
          {
            gPert[i].a[k-NCell+1] = A->a[(k+1)*A->BlockSize*A->BlockSize + i*A->BlockSize + i].Re;
          }
        }
      }
      
      /* now b values */
      ib=A->BlockSize*A->BlockSize;
      for(k=0;k<Min(NBlockMax-1,A->N-1);k++)
      {
        if(k < (NCell-1))
        {
          ifer=NFermionShift+(k+1)*A->BlockSize;
          for(i=0;i<A->BlockSize;i++)
          {
            jfer=NFermionShift+k*A->BlockSize;
            for(j=0;j<A->BlockSize;j++)
            {
              oadd[1]=CrAnMask+ifer;
              oadd[0]=jfer;
              ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->b[ib].Re, A->b[ib].Im, &QDetNotOrderedFcFaBcBa);
              oadd[1]=CrAnMask+jfer;
              oadd[0]=ifer;
              ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->b[ib].Re, A->b[ib].Im, &QDetNotOrderedFcFaBcBa);
              ib++;
              jfer++;
            }
            ifer++;
          }
        }
        else if(k == (NCell-1))
        {
          /* Hopping from NCell to NCell+1 goes in VPert */
          for(i = 0; i < A->BlockSize; i++)
          {
            oadd[0] = CrAnMask + 2*k*A->BlockSize + i;
            ComplexOperatorAddLaderNormalOrder(&VPert[i],oadd,1,A->b[(k+1)*A->BlockSize*A->BlockSize+i*A->BlockSize+i].Re,A->b[(k+1)*A->BlockSize*A->BlockSize+i*A->BlockSize+i].Im, &QDetNotOrderedFcFaBcBa);
          }
        }
        else
        {
          /* The rest goes in gPert */
          for(i = 0; i < A->BlockSize; i++)
          {
            gPert[i].b[(k-NCell+1)] = A->b[(k+1)*A->BlockSize*A->BlockSize+i*A->BlockSize+i].Re;
          }
        }
      }
    }
    else /* real with spin */
    {
      /* a values first*/
      ia=A->BlockSize*A->BlockSize;
      for(k=0;k<Min(NBlockMax,A->N);k++)
      {
        if(k < NCell)
        {
          /* Part up to NCell goes in H */
          ifer=NFermionShift+2*k*A->BlockSize;
          for(i=0;i<A->BlockSize;i++)
          {
            jfer=NFermionShift+2*k*A->BlockSize;
            for(j=0;j<A->BlockSize;j++)
            {
              oadd[1]=CrAnMask+ifer;
              oadd[0]=jfer;
              ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->a[ia].Re, A->a[ia].Im, &QDetNotOrderedFcFaBcBa);
              oadd[1]=CrAnMask+ifer+A->BlockSize;
              oadd[0]=jfer+A->BlockSize;
              ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->a[ia].Re, A->a[ia].Im, &QDetNotOrderedFcFaBcBa);
              ia++;
              jfer++;
            }
            ifer++;
          }
        }
        else
        {
          /* Blocks after NCell go in gPert*/
          for(i = 0; i < A->BlockSize; i++)
          {
            gPert[i].a[k-NCell+1] = A->a[(k+1)*A->BlockSize*A->BlockSize + i*A->BlockSize + i].Re;
            gPert[i+A->BlockSize].a[k-NCell+1] = A->a[(k+1)*A->BlockSize*A->BlockSize + i*A->BlockSize + i].Re;
          }
        }
      }
      
      /* now b values */
      ib=A->BlockSize*A->BlockSize;
      for(k=0;k<Min(NBlockMax-1,A->N-1);k++)
      {
        if(k < (NCell-1))
        {
          ifer=NFermionShift+2*(k+1)*A->BlockSize;
          for(i=0;i<A->BlockSize;i++)
          {
            jfer=NFermionShift+2*k*A->BlockSize;
            for(j=0;j<A->BlockSize;j++)
            {
              oadd[1]=CrAnMask+ifer;
              oadd[0]=jfer;
              ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->b[ib].Re, A->b[ib].Im, &QDetNotOrderedFcFaBcBa);
              oadd[1]=CrAnMask+jfer;
              oadd[0]=ifer;
              ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->b[ib].Re, A->b[ib].Im, &QDetNotOrderedFcFaBcBa);
              oadd[1]=CrAnMask+ifer+A->BlockSize;
              oadd[0]=jfer+A->BlockSize;
              ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->b[ib].Re, A->b[ib].Im, &QDetNotOrderedFcFaBcBa);
              oadd[1]=CrAnMask+jfer+A->BlockSize;
              oadd[0]=ifer+A->BlockSize;
              ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->b[ib].Re, A->b[ib].Im, &QDetNotOrderedFcFaBcBa);
              ib++;
              jfer++;
            }
            ifer++;
          }
        }
        else if(k == (NCell-1))
        {
          /* Hopping from NCell to NCell+1 goes in VPert */
          for(i = 0; i < A->BlockSize; i++)
          {
            oadd[0] = CrAnMask + 2*k*A->BlockSize + i;
            ComplexOperatorAddLaderNormalOrder(&VPert[i],oadd,1,A->b[(k+1)*A->BlockSize*A->BlockSize+i*A->BlockSize+i].Re,A->b[(k+1)*A->BlockSize*A->BlockSize+i*A->BlockSize+i].Im, &QDetNotOrderedFcFaBcBa);
            oadd[0] = CrAnMask + 2*k*A->BlockSize + i + A->BlockSize;
            ComplexOperatorAddLaderNormalOrder(&VPert[i+A->BlockSize],oadd,1,A->b[(k+1)*A->BlockSize*A->BlockSize+i*A->BlockSize+i].Re,A->b[(k+1)*A->BlockSize*A->BlockSize+i*A->BlockSize+i].Im, &QDetNotOrderedFcFaBcBa);
          }
        }
        else
        {
          /* The rest goes in gPert */
          for(i = 0; i < A->BlockSize; i++)
          {
            gPert[i].b[(k-NCell+1)] = A->b[(k+1)*A->BlockSize*A->BlockSize+i*A->BlockSize+i].Re;
            gPert[i+A->BlockSize].b[(k-NCell+1)] = A->b[(k+1)*A->BlockSize*A->BlockSize+i*A->BlockSize+i].Re;
          }
        }
      }
    }
  }
  else /* H is real */
  {
    if(Spin==0)
    {
      /* a values first*/
      ia=A->BlockSize*A->BlockSize;
      for(k=0;k<Min(NBlockMax,A->N);k++)
      {
        if(k < NCell)
        {
          /* Part up to NCell goes in H */
          ifer=NFermionShift+k*A->BlockSize;
          for(i=0;i<A->BlockSize;i++)
          {
            jfer=NFermionShift+k*A->BlockSize;
            for(j=0;j<A->BlockSize;j++)
            {
              oadd[1]=CrAnMask+ifer;
              oadd[0]=jfer;
              RealOperatorAddLaderNormalOrder(H, oadd, 2, A->aRe[ia], &QDetNotOrderedFcFaBcBa);
              ia++;
              jfer++;
            }
            ifer++;
          }
        }
        else
        {
          /* Blocks after NCell go in gPert*/
          for(i = 0; i < A->BlockSize; i++)
          {
            gPert[i].a[k-NCell+1] = A->aRe[(k+1)*A->BlockSize*A->BlockSize + i*A->BlockSize + i];
          }
        }
      }
      
      /* now b values */
      ib=A->BlockSize*A->BlockSize;
      for(k=0;k<Min(NBlockMax-1,A->N-1);k++)
      {
        if(k < (NCell-1))
        {
          ifer=NFermionShift+(k+1)*A->BlockSize;
          for(i=0;i<A->BlockSize;i++)
          {
            jfer=NFermionShift+k*A->BlockSize;
            for(j=0;j<A->BlockSize;j++)
            {
              oadd[1]=CrAnMask+ifer;
              oadd[0]=jfer;
              RealOperatorAddLaderNormalOrder(H, oadd, 2, A->bRe[ib], &QDetNotOrderedFcFaBcBa);
              oadd[1]=CrAnMask+jfer;
              oadd[0]=ifer;
              RealOperatorAddLaderNormalOrder(H, oadd, 2, A->bRe[ib], &QDetNotOrderedFcFaBcBa);
              ib++;
              jfer++;
            }
            ifer++;
          }
        }
        else if(k == (NCell-1))
        {
          /* Hopping from NCell to NCell+1 goes in VPert */
          for(i = 0; i < A->BlockSize; i++)
          {
            oadd[0] = CrAnMask + 2*k*A->BlockSize + i;
            RealOperatorAddLaderNormalOrder(&VPert[i],oadd,1,A->bRe[(k+1)*A->BlockSize*A->BlockSize+i*A->BlockSize+i], &QDetNotOrderedFcFaBcBa);
          }
        }
        else
        {
          /* The rest goes in gPert */
          for(i = 0; i < A->BlockSize; i++)
          {
            gPert[i].b[(k-NCell+1)] = A->bRe[(k+1)*A->BlockSize*A->BlockSize+i*A->BlockSize+i];
          }
        }
      }
    }
    else /* real with spin */
    {
      /* a values first*/
      ia=A->BlockSize*A->BlockSize;
      for(k=0;k<Min(NBlockMax,A->N);k++)
      {
        if(k < NCell)
        {
          /* Part up to NCell goes in H */
          ifer=NFermionShift+2*k*A->BlockSize;
          for(i=0;i<A->BlockSize;i++)
          {
            jfer=NFermionShift+2*k*A->BlockSize;
            for(j=0;j<A->BlockSize;j++)
            {
              oadd[1]=CrAnMask+ifer;
              oadd[0]=jfer;
              RealOperatorAddLaderNormalOrder(H, oadd, 2, A->aRe[ia], &QDetNotOrderedFcFaBcBa);
              oadd[1]=CrAnMask+ifer+A->BlockSize;
              oadd[0]=jfer+A->BlockSize;
              RealOperatorAddLaderNormalOrder(H, oadd, 2, A->aRe[ia], &QDetNotOrderedFcFaBcBa);
              ia++;
              jfer++;
            }
            ifer++;
          }
        }
        else
        {
          /* Blocks after NCell go in gPert*/
          for(i = 0; i < A->BlockSize; i++)
          {
            gPert[i].a[k-NCell+1] = A->aRe[(k+1)*A->BlockSize*A->BlockSize + i*A->BlockSize + i];
            gPert[i+A->BlockSize].a[k-NCell+1] = A->aRe[(k+1)*A->BlockSize*A->BlockSize + i*A->BlockSize + i];
          }
        }
      }
      
      /* now b values */
      ib=A->BlockSize*A->BlockSize;
      for(k=0;k<Min(NBlockMax-1,A->N-1);k++)
      {
        if(k < (NCell-1))
        {
          ifer=NFermionShift+2*(k+1)*A->BlockSize;
          for(i=0;i<A->BlockSize;i++)
          {
            jfer=NFermionShift+2*k*A->BlockSize;
            for(j=0;j<A->BlockSize;j++)
            {
              oadd[1]=CrAnMask+ifer;
              oadd[0]=jfer;
              RealOperatorAddLaderNormalOrder(H, oadd, 2, A->bRe[ib], &QDetNotOrderedFcFaBcBa);
              oadd[1]=CrAnMask+jfer;
              oadd[0]=ifer;
              RealOperatorAddLaderNormalOrder(H, oadd, 2, A->bRe[ib], &QDetNotOrderedFcFaBcBa);
              oadd[1]=CrAnMask+ifer+A->BlockSize;
              oadd[0]=jfer+A->BlockSize;
              RealOperatorAddLaderNormalOrder(H, oadd, 2, A->bRe[ib], &QDetNotOrderedFcFaBcBa);
              oadd[1]=CrAnMask+jfer+A->BlockSize;
              oadd[0]=ifer+A->BlockSize;
              RealOperatorAddLaderNormalOrder(H, oadd, 2, A->bRe[ib], &QDetNotOrderedFcFaBcBa);
              ib++;
              jfer++;
            }
            ifer++;
          }
        }
        else if(k == (NCell-1))
        {
          /* Hopping from NCell to NCell+1 goes in VPert */
          for(i = 0; i < A->BlockSize; i++)
          {
            oadd[0] = CrAnMask + 2*k*A->BlockSize + i;
            RealOperatorAddLaderNormalOrder(&VPert[i],oadd,1,A->bRe[(k+1)*A->BlockSize*A->BlockSize+i*A->BlockSize+i], &QDetNotOrderedFcFaBcBa);
            oadd[0] = CrAnMask + 2*k*A->BlockSize + i + A->BlockSize;
            RealOperatorAddLaderNormalOrder(&VPert[i+A->BlockSize],oadd,1,A->bRe[(k+1)*A->BlockSize*A->BlockSize+i*A->BlockSize+i], &QDetNotOrderedFcFaBcBa);
          }
        }
        else
        {
          /* The rest goes in gPert */
          for(i = 0; i < A->BlockSize; i++)
          {
            gPert[i].b[(k-NCell+1)] = A->bRe[(k+1)*A->BlockSize*A->BlockSize+i*A->BlockSize+i];
            gPert[i+A->BlockSize].b[(k-NCell+1)] = A->bRe[(k+1)*A->BlockSize*A->BlockSize+i*A->BlockSize+i];
          }
        }
      }
    }
  }
  return 0;
}

/*  Starting from a BlockAndersonMatrix creates an impurity Hamiltonian
 *  H on input not initiallized
 *  NBlockMax the maximal number of blocks in the Hamiltonian, this includes the impurity.
 *   NBlockMax=2 would thus correspond to a triditional ligand field calculation
 *  if Spin=0 each block in A is taken once into the Hamiltonian
 *  If Spin =1 (!0) each block is copied twice, once for spin up, once for spin down. The order is block up, block down
 *  NFermionShift is the numbers of Fermions shifted. (i.e. if you need a 2p core before a d impurity block NFermionShift=6 and the p-shell is from 0 to 5, the d-shell form 6 to 15.
 *  NFermion is the total number of Fermions in the basis. 
 *  If NFermion=0 on input the number will be calculated such that H fits.
 */
uint32_t BlockAndersonMatrixToImpurityHamiltonian(OperatorType *const H, BlockAndersonMatrixType const *A, uint32_t const NBlockMax, uint32_t const Spin, uint32_t const NFermionShift, uint32_t const NFermion)
{
  uint32_t i,j,k,ia,ifer,jfer,ib;
  uint16_t oadd[2];
  InitNamedOperator(H, "Impurity Hamiltonian");
  H->NBoson=0;
  if(Spin==0)
  {
    H->NFermion=NFermionShift+NBlockMax*A->BlockSize;
  }
  else
  {
    H->NFermion=NFermionShift+2*NBlockMax*A->BlockSize;
  }
  if(NFermion!=0 )
  {
    if(NFermion >= H->NFermion)
    {
      H->NFermion=NFermion;
    }
    else
    {
      printf("Warning in BlockTridiagonalMatrixToImpurityHamiltonian\n");
      printf("Function called with NFermion=%i, but total number of Fermions must be larger than %i\n",NFermion,H->NFermion);
      printf("Number of Fermions enlarged\n");
    }
  }
  InitOperatorFixedLength(H, 2, 8);
  if(A->QComplex)
  {
    MakeRealOperatorFixedLengthComplex(H, 2);
    if(Spin==0)
    {
      ia=A->BlockSize*A->BlockSize;
      for(k=0;k<Min(NBlockMax,A->N);k++)
      {
        ifer=NFermionShift+k*A->BlockSize;
        for(i=0;i<A->BlockSize;i++)
        {
          jfer=NFermionShift+k*A->BlockSize;
          for(j=0;j<A->BlockSize;j++)
          {
            oadd[1]=CrAnMask+ifer;
            oadd[0]=jfer;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->a[ia].Re, A->a[ia].Im, &QDetNotOrderedFcFaBcBa);
            ia++;
            jfer++;
          }
          ifer++;
        }
      }
      ib=A->BlockSize*A->BlockSize;
      for(k=0;k<Min(NBlockMax-1,A->N-1);k++)
      {
        ifer=NFermionShift;//+(k+1)*A->BlockSize;
        for(i=0;i<A->BlockSize;i++)
        {
          jfer=NFermionShift+(k+1)*A->BlockSize;
          for(j=0;j<A->BlockSize;j++)
          {
            oadd[1]=CrAnMask+ifer;
            oadd[0]=jfer;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->b[ib].Re, A->b[ib].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+jfer;
            oadd[0]=ifer;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->b[ib].Re,-A->b[ib].Im, &QDetNotOrderedFcFaBcBa);
            ib++;
            jfer++;
          }
          ifer++;
        }
      }
    }
    else /* Complex with spin */
    {
      ia=A->BlockSize*A->BlockSize;
      for(k=0;k<Min(NBlockMax,A->N);k++)
      {
        ifer=NFermionShift+2*k*A->BlockSize;
        for(i=0;i<A->BlockSize;i++)
        {
          jfer=NFermionShift+2*k*A->BlockSize;
          for(j=0;j<A->BlockSize;j++)
          {
            oadd[1]=CrAnMask+ifer;
            oadd[0]=jfer;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->a[ia].Re, A->a[ia].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+ifer+A->BlockSize;
            oadd[0]=jfer+A->BlockSize;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->a[ia].Re, A->a[ia].Im, &QDetNotOrderedFcFaBcBa);
            ia++;
            jfer++;
          }
          ifer++;
        }
      }
      ib=A->BlockSize*A->BlockSize;
      for(k=0;k<Min(NBlockMax-1,A->N-1);k++)
      {
        ifer=NFermionShift;//+2*(k+1)*A->BlockSize;
        for(i=0;i<A->BlockSize;i++)
        {
          jfer=NFermionShift+2*(k+1)*A->BlockSize;
          for(j=0;j<A->BlockSize;j++)
          {
            oadd[1]=CrAnMask+ifer;
            oadd[0]=jfer;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->b[ib].Re, A->b[ib].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+jfer;
            oadd[0]=ifer;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->b[ib].Re,-A->b[ib].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+ifer+A->BlockSize;
            oadd[0]=jfer+A->BlockSize;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->b[ib].Re, A->b[ib].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+jfer+A->BlockSize;
            oadd[0]=ifer+A->BlockSize;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->b[ib].Re,-A->b[ib].Im, &QDetNotOrderedFcFaBcBa);
            ib++;
            jfer++;
          }
          ifer++;
        }
      }
    }
  }
  else /* H is real */
  {
    if(Spin==0)
    {
      ia=A->BlockSize*A->BlockSize;
      for(k=0;k<Min(NBlockMax,A->N);k++)
      {
        ifer=NFermionShift+k*A->BlockSize;
        for(i=0;i<A->BlockSize;i++)
        {
          jfer=NFermionShift+k*A->BlockSize;
          for(j=0;j<A->BlockSize;j++)
          {
            oadd[1]=CrAnMask+ifer;
            oadd[0]=jfer;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->aRe[ia], &QDetNotOrderedFcFaBcBa);
            ia++;
            jfer++;
          }
          ifer++;
        }
      }
      ib=A->BlockSize*A->BlockSize;
      for(k=0;k<Min(NBlockMax-1,A->N-1);k++)
      {
        ifer=NFermionShift;//+(k+1)*A->BlockSize;
        for(i=0;i<A->BlockSize;i++)
        {
          jfer=NFermionShift+(k+1)*A->BlockSize;
          for(j=0;j<A->BlockSize;j++)
          {
            oadd[1]=CrAnMask+ifer;
            oadd[0]=jfer;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->bRe[ib], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+jfer;
            oadd[0]=ifer;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->bRe[ib], &QDetNotOrderedFcFaBcBa);
            ib++;
            jfer++;
          }
          ifer++;
        }
      }
    }
    else /* real with spin */
    {
      ia=A->BlockSize*A->BlockSize;
      for(k=0;k<Min(NBlockMax,A->N);k++)
      {
        ifer=NFermionShift+2*k*A->BlockSize;
        for(i=0;i<A->BlockSize;i++)
        {
          jfer=NFermionShift+2*k*A->BlockSize;
          for(j=0;j<A->BlockSize;j++)
          {
            oadd[1]=CrAnMask+ifer;
            oadd[0]=jfer;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->aRe[ia], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+ifer+A->BlockSize;
            oadd[0]=jfer+A->BlockSize;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->aRe[ia], &QDetNotOrderedFcFaBcBa);
            ia++;
            jfer++;
          }
          ifer++;
        }
      }
      ib=A->BlockSize*A->BlockSize;
      for(k=0;k<Min(NBlockMax-1,A->N-1);k++)
      {
        ifer=NFermionShift;//+2*(k+1)*A->BlockSize;
        for(i=0;i<A->BlockSize;i++)
        {
          jfer=NFermionShift+2*(k+1)*A->BlockSize;
          for(j=0;j<A->BlockSize;j++)
          {
            oadd[1]=CrAnMask+ifer;
            oadd[0]=jfer;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->bRe[ib], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+jfer;
            oadd[0]=ifer;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->bRe[ib], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+ifer+A->BlockSize;
            oadd[0]=jfer+A->BlockSize;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->bRe[ib], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+jfer+A->BlockSize;
            oadd[0]=ifer+A->BlockSize;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->bRe[ib], &QDetNotOrderedFcFaBcBa);
            ib++;
            jfer++;
          }
          ifer++;
        }
      }
    }
  }
  return 0;
}


/*  From BasicMath_MatrixTypeConversion. Only difference is it calls my MatrixTransformBLAS instead of MatrixTransform
 * 
 * 
 * Returns the Generalized Eigensystem of a non-orthnormal matrix.
 *  Given a basis psi_i
 *  Matrix M contains the values psi_i^* . H . psi_j
 *  Matrix S contains the values psi_i^* . psi_j
 *
 *  on exit val defines a diagonal matrix H and
 *  H = T . M . T^{\dag}
 *  I = T . S . T^{\dag}
 *  S = Sroot . Sroot^{\dag}
 *
 */
uint32_t GeneralizedEigensystemSingularBLAS(MatrixType const M, MatrixType const S, MatrixType * const T,
                                        MatrixType * const Sroot, double *val)
{
  uint32_t i, j, NewDim, retval = 0;
  double Sval[S.Ni], bound;
  MatrixType EigenSysS, TMT, TST, Ttemp;
  MatrixCopy(&TST, S);
  MatrixScalarMultiply(&TST, -1, 0);
  Eigensystem(TST, &EigenSysS, Sval);
  MatrixScalarMultiply(&TST, 0, 0);
  MatrixTransformBLAS(&TMT, M, EigenSysS);
  NewDim = 1;
  bound = 0;
  for (i = 0; i < M.Ni; i++)
  {
    bound -= Sval[i];
  }
  bound *= DBL_EPSILON * 100;
  bound += DBL_MIN * 100;
  for (i = 0; i < M.Ni; i++)
  {
    if (-Sval[i] > bound)
    {
      TST.VRe[i][i] = -Sval[i];
      NewDim = i + 1;
    }
    else
    {
      for (j = 0; j < M.Nj; j++)
      {
        TMT.VRe[i][j] = 0;
        TMT.VRe[j][i] = 0;
      }
      TST.VRe[i][i] = 1;
    }
  }
  printf("Singularity border =%10.3E : Number of states kept = %i\n", bound, NewDim);
  printf("Norm of basis ( 0...7): ");
  for (i = 0; i < Min(8, S.Ni); i++)
  {
    printf("%10.3E ", -Sval[i]);
  }
  printf("..........\n");
  if (GeneralizedEigenvalues(TMT, TST, val))
  {
    printf("GeneralizedEigenvalues failed in GeneralizedEigensystemSingular\n");
    retval = 1;
  }
  for (i = NewDim; i < M.Ni; i++)
  {
    TMT.VRe[i][i] = val[M.Ni - 1] + 100 * DBL_EPSILON * Abs(val[M.Ni - 1]);
  }
  if (GeneralizedEigensystem(TMT, TST, &Ttemp, Sroot, val))
  {
    printf("GeneralizedEigensystem failed in GeneralizedEigensystemSingular\n");
    retval = 1;
  }
  for (i = NewDim; i < M.Ni; i++)
  {
    for (j = 0; j < M.Ni; j++)
    {
      Ttemp.VRe[j][i] = 0;
    }
  }
  MatrixMultiply(T, Ttemp, EigenSysS);
  FreeMatrix(&EigenSysS);
  FreeMatrix(&TMT);
  FreeMatrix(&TST);
  FreeMatrix(&Ttemp);
  return retval;
}





/* From Lanczos.c  Only difference is it calls my GeneralizedEigensystemSingularBLAS instead of the existing version*/
uint32_t LanczosTriDiagonalizeDysonLowMemBLAS(OperatorType const O, WaveFunctionType psi, double epsilon,
                                          uint32_t const NKrylov, OperatorType const *V, TriDiagonalMatrixType const *g,
                                          uint32_t const NV, TriDiagonalMatrixType * Gout)
{
  uint32_t i, j;
  TriDiagonalMatrixType OKrylov;
  MatrixType Stemp, Htemp, HtempEigenFun, StempRoot;
  MatrixType Hval, *VMat;
  WaveFunctionType *PsiKrylov;
  OperatorType Otemp;
  TimeType Time;
  InitializeTiming(&Time);
  BeginTiming("Total", &Time);
  printf("\nStart of LanczosTriDiagonalizeDysonLowMem\n\n");
  /* Initiallize matrices */
  Htemp.Ni = NKrylov;
  Htemp.Nj = NKrylov;
  Htemp.QComplex = (O.QComplex || psi.QComplex);
  for (i = 0; i < NV; i++)
  {
    Htemp.QComplex = (Htemp.QComplex || V[i].QComplex);
  }
  if (InitMatrix(&Htemp))
  {
    printf("InitMatrix failed in LanczosTriDiagonalizeDyson\n");
    fflush(stdout);
    return 1;
  };
  Stemp.Ni = NKrylov;
  Stemp.Nj = NKrylov;
  Stemp.QComplex = Htemp.QComplex;
  if (InitMatrix(&Stemp))
  {
    printf("InitMatrix failed in LanczosTriDiagonalizeDyson\n");
    fflush(stdout);
    return 1;
  };
  Hval.Ni = NV;
  Hval.Nj = NKrylov;
  Hval.QComplex = 0;
  if (InitMatrix(&Hval))
  {
    printf("InitMatrix failed in LanczosTriDiagonalizeDyson\n");
    fflush(stdout);
    return 1;
  };
  if ((VMat = calloc(NV, sizeof(MatrixType))) == NULL)
  {
    printf("calloc failed in LanczosTriDiagonalizeDyson\n");
    fflush(stdout);
    return 1;
  };
  /* Create Krylov basis of O */
  OKrylov.N = NKrylov;
  InitTriDiagonalMatrix(&OKrylov);
  if ((PsiKrylov = calloc(NKrylov, sizeof(WaveFunctionType))) == NULL)
  {
    printf("calloc failed in LanczosTriDiagonalizeDyson\n");
    fflush(stdout);
    return 1;
  };
  BeginTiming("Krylov", &Time);
  LanczosTriDiagonalizeKrylov(O, psi, epsilon, &OKrylov, PsiKrylov);
  EndTiming("Krylov", &Time);
  /* loop over all V[i] */
  BeginTiming("Creation of matrix elements", &Time);
  for (i = 0; i < NV; i++)
  {
    printf("Create perturbation %i from %i\n", i + 1, NV);
    /* create extended non orthonormal matrices on original basis */
    /* initiallize Otemp for V^{\dag} V */
    InitOperatorForMult(&Otemp, 2, &V[i], &V[i]);
    /* create Stemp */
    OperatorMultiply2NormalOrder(&Otemp, V[i], V[i], 'D', 'N', &QDetNotOrderedFcFaBcBa);
    PsiOperatorPsiMatrix(Otemp, PsiKrylov, &Stemp);
    FreeOperator(&Otemp);
    /* initiallize Otemp for V^{\dag} O V */
    InitOperatorForMult(&Otemp, 3, &V[i], &O, &V[i]);
    /* create Htemp */
    OperatorMultiply3NormalOrder(&Otemp, V[i], O, V[i], 'D', 'N', 'N', &QDetNotOrderedFcFaBcBa);
    PsiOperatorPsiMatrix(Otemp, PsiKrylov, &Htemp);
    FreeOperator(&Otemp);
    /* diagonalize Htemp with repsect to Stemp */
    GeneralizedEigensystemSingularBLAS(Htemp, Stemp, &HtempEigenFun, &StempRoot, Hval.VRe[i]);
    MatrixMultiply(&(VMat[i]), HtempEigenFun, Stemp);
    /* Free HtempEigenFun */
    FreeMatrix(&HtempEigenFun);
  }
  /* Free PsiKrylov */
  for (j = 0; j < NKrylov; j++)
  {
    WaveFunctionFree(&(PsiKrylov[j]));
  }
  free(PsiKrylov);
  FreeMatrix(&Htemp);
  FreeMatrix(&Stemp);
  EndTiming("Creation of matrix elements", &Time);
  /* create Gout */
  BeginTiming("Final Tridiagonalize", &Time);
  if (TridiagonalizeImpurityCoupledToBathDyson(OKrylov, VMat, Hval, g, NV, Gout))
  {
    printf("TridiagonalizeImpurityCoupledToBathDyson failed in LanczosTriDiagonalizeDyson\n");
    FreeTriDiagonalMatrix(&OKrylov);
    FreeMatrix(&Hval);
    for (i = 0; i < NV; i++)
    {
      FreeMatrix(&(VMat[i]));
    }
    free(VMat);
    return 1;
  }
  EndTiming("Final Tridiagonalize", &Time);
  /* Free things allocated */
  FreeTriDiagonalMatrix(&OKrylov);
  FreeMatrix(&Hval);
  for (i = 0; i < NV; i++)
  {
    FreeMatrix(&(VMat[i]));
  }
  free(VMat);
  EndTiming("Total", &Time);
  PrintTiming(Time);
  FreeTiming(&Time);
  return 0;
}