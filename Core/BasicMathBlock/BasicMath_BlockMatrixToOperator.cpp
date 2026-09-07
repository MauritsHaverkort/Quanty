
/*  Starting from a BlockTridiagonalMatrix creates an impurity Hamiltonian
 *  H on input not initiallized
 *  if Spin=0 each block in A is taken once into the Hamiltonian
 *  If Spin =1 (!0) each block is copied twice, once for spin up, once for spin down. The order is block up, block down
 *  NFermionShift is the numbers of Fermions shifted. (i.e. if you need a 2p core before a d impurity block NFermionShift=6 and the p-shell is from 0 to 5, the d-shell form 6 to 15.
 *  NFermion is the total number of Fermions in the basis. 
 *  If NFermion=0 on input the number will be calculated such that H fits.
 */
uint32_t BlockTridiagonalMatrixToImpurityHamiltonian(OperatorType *const H, BlockTriDiagonalMatrixType const *A, bool const Spin, uint32_t const NFermionShift, uint32_t const NFermion)
{
  uint32_t i,j,k,ia,ifer,jfer,ib;
  uint16_t oadd[2];
  InitNamedOperator(H, "Impurity Hamiltonian");
  H->NBoson=0;
  if(Spin==false)
  {
    H->NFermion=NFermionShift+A->N*A->BlockSize;
  }
  else
  {
    H->NFermion=NFermionShift+2*A->N*A->BlockSize;
  }
  if(NFermion!=0 )
  {
    if(NFermion>H->NFermion)
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
  H->NWord=Min(H->NFermion,1+(H->NFermion-1)/8);
  if(A->N==0) return 0;
  InitOperatorFixedLength(H, 2, 8);
  if(A->QComplex)
  {
    MakeRealOperatorFixedLengthComplex(H, 2);
    if(Spin==0)
    {
      ia=A->BlockSize*A->BlockSize;
      for(k=0;k<A->N;k++)
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
      for(k=0;k<A->N-1;k++)
      {
        ifer=NFermionShift+k*A->BlockSize;
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
      for(k=0;k<A->N;k++)
      {
        ifer=NFermionShift+k*A->BlockSize;
        for(i=0;i<A->BlockSize;i++)
        {
          jfer=NFermionShift+k*A->BlockSize;
          for(j=0;j<A->BlockSize;j++)
          {
            oadd[1]=CrAnMask+2*ifer;
            oadd[0]=2*jfer;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->a[ia].Re, A->a[ia].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+2*ifer+1;
            oadd[0]=2*jfer+1;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->a[ia].Re, A->a[ia].Im, &QDetNotOrderedFcFaBcBa);
            ia++;
            jfer++;
          }
          ifer++;
        }
      }
      ib=A->BlockSize*A->BlockSize;
      for(k=0;k<A->N-1;k++)
      {
        ifer=NFermionShift+k*A->BlockSize;
        for(i=0;i<A->BlockSize;i++)
        {
          jfer=NFermionShift+(k+1)*A->BlockSize;
          for(j=0;j<A->BlockSize;j++)
          {
            oadd[1]=CrAnMask+2*ifer;
            oadd[0]=2*jfer;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->b[ib].Re, A->b[ib].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+2*jfer;
            oadd[0]=2*ifer;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->b[ib].Re,-A->b[ib].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+2*ifer+1;
            oadd[0]=2*jfer+1;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->b[ib].Re, A->b[ib].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+2*jfer+1;
            oadd[0]=2*ifer+1;
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
      for(k=0;k<A->N;k++)
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
      for(k=0;k<A->N-1;k++)
      {
        ifer=NFermionShift+k*A->BlockSize;
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
      for(k=0;k<A->N;k++)
      {
        ifer=NFermionShift+k*A->BlockSize;
        for(i=0;i<A->BlockSize;i++)
        {
          jfer=NFermionShift+k*A->BlockSize;
          for(j=0;j<A->BlockSize;j++)
          {
            oadd[1]=CrAnMask+2*ifer;
            oadd[0]=2*jfer;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->aRe[ia], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+2*ifer+1;
            oadd[0]=2*jfer+1;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->aRe[ia], &QDetNotOrderedFcFaBcBa);
            ia++;
            jfer++;
          }
          ifer++;
        }
      }
      ib=A->BlockSize*A->BlockSize;
      for(k=0;k<A->N-1;k++)
      {
        ifer=NFermionShift+k*A->BlockSize;
        for(i=0;i<A->BlockSize;i++)
        {
          jfer=NFermionShift+(k+1)*A->BlockSize;
          for(j=0;j<A->BlockSize;j++)
          {
            oadd[1]=CrAnMask+2*ifer;
            oadd[0]=2*jfer;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->bRe[ib], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+2*ifer+1;
            oadd[0]=2*jfer+1;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->bRe[ib], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+2*jfer;
            oadd[0]=2*ifer;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->bRe[ib], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+2*jfer+1;
            oadd[0]=2*ifer+1;
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


/* 
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
  //TODO: test
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
  
  H->NWord = Min(H->NFermion, 1 + (H->NFermion - 1) / 8);
  
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
uint32_t BlockAndersonMatrixToImpurityHamiltonian(OperatorType *const H, BlockAndersonMatrixType const *A, bool const Spin, uint32_t const NFermionShift, uint32_t const NFermion)
{
  uint32_t i,j,k,ia,ifer,jfer,ib;
  uint16_t oadd[2];
  InitNamedOperator(H, "Impurity Hamiltonian");
  H->NBoson=0;
  if(Spin==false)
  {
    H->NFermion=NFermionShift+A->N*A->BlockSize;
  }
  else
  {
    H->NFermion=NFermionShift+2*A->N*A->BlockSize;
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
  H->NWord=Min(H->NFermion,1+(H->NFermion-1)/8);
  if(A->N==0) return 0;
  InitOperatorFixedLength(H, 2, 8);
  if(A->QComplex)
  {
    MakeRealOperatorFixedLengthComplex(H, 2);
    if(Spin==0)
    {
      ia=A->BlockSize*A->BlockSize;
      for(k=0;k<A->N;k++)
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
      for(k=0;k<A->N-1;k++)
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
      for(k=0;k<A->N;k++)
      {
        ifer=NFermionShift+k*A->BlockSize;
        for(i=0;i<A->BlockSize;i++)
        {
          jfer=NFermionShift+k*A->BlockSize;
          for(j=0;j<A->BlockSize;j++)
          {
            oadd[1]=CrAnMask+2*ifer;
            oadd[0]=2*jfer;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->a[ia].Re, A->a[ia].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+2*ifer+1;
            oadd[0]=2*jfer+1;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->a[ia].Re, A->a[ia].Im, &QDetNotOrderedFcFaBcBa);
            ia++;
            jfer++;
          }
          ifer++;
        }
      }
      ib=A->BlockSize*A->BlockSize;
      for(k=0;k<A->N-1;k++)
      {
        ifer=NFermionShift;//+(k+1)*A->BlockSize;
        for(i=0;i<A->BlockSize;i++)
        {
          jfer=NFermionShift+(k+1)*A->BlockSize;
          for(j=0;j<A->BlockSize;j++)
          {
            oadd[1]=CrAnMask+2*ifer;
            oadd[0]=2*jfer;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->b[ib].Re, A->b[ib].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+2*jfer;
            oadd[0]=2*ifer;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->b[ib].Re,-A->b[ib].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+2*ifer+1;
            oadd[0]=2*jfer+1;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->b[ib].Re, A->b[ib].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+2*jfer+1;
            oadd[0]=2*ifer+1;
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
      for(k=0;k<A->N;k++)
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
      for(k=0;k<A->N-1;k++)
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
      for(k=0;k<A->N;k++)
      {
        ifer=NFermionShift+k*A->BlockSize;
        for(i=0;i<A->BlockSize;i++)
        {
          jfer=NFermionShift+k*A->BlockSize;
          for(j=0;j<A->BlockSize;j++)
          {
            oadd[1]=CrAnMask+2*ifer;
            oadd[0]=2*jfer;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->aRe[ia], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+2*ifer+1;
            oadd[0]=2*jfer+1;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->aRe[ia], &QDetNotOrderedFcFaBcBa);
            ia++;
            jfer++;
          }
          ifer++;
        }
      }
      ib=A->BlockSize*A->BlockSize;
      for(k=0;k<A->N-1;k++)
      {
        ifer=NFermionShift;//+(k+1)*A->BlockSize;
        for(i=0;i<A->BlockSize;i++)
        {
          jfer=NFermionShift+(k+1)*A->BlockSize;
          for(j=0;j<A->BlockSize;j++)
          {
            oadd[1]=CrAnMask+2*ifer;
            oadd[0]=2*jfer;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->bRe[ib], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+2*ifer+1;
            oadd[0]=2*jfer+1;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->bRe[ib], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+2*jfer;
            oadd[0]=2*ifer;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->bRe[ib], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+2*jfer+1;
            oadd[0]=2*ifer+1;
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


/*  Starting from a BlockTridiagonalMatrix creates an impurity Hamiltonian
 *  H on input not initiallized
 *
 *  The impurity can be found as the first block
 *  The second block is the partially filled bath
 *  Block 3 to 2+NBlockMaxVal are the occupied valence bath orbitals
 *  Block 3+NBlockMaxVal to 2+NBlockMaxVal+NBlockMaxCon are the unoccupied conduction bath orbitals
 *
 *  if Spin=0 each block in A is taken once into the Hamiltonian
 *  If Spin =1 (!0) each block is copied twice, once for spin up, once for spin down. The order is block up, block down
 *  NFermionShift is the numbers of Fermions shifted. (i.e. if you need a 2p core before a d impurity block NFermionShift=6 and the p-shell is from 0 to 5, the d-shell form 6 to 15.
 *  NFermion is the total number of Fermions in the basis.
 *  If NFermion=0 on input the number will be calculated such that H fits.
 */
uint32_t BlockNaturalImpurityOrbitalMatrixToImpurityHamiltonian(OperatorType *const H, BlockNaturalImpurityOrbitalMatrixType const *A, bool const Spin, uint32_t const NFermionShift, uint32_t const NFermion)
{
  bool Deflation=false; // this should be a parameter in the function, but then we should also implement it when the operator is created
  uint32_t i,j,k,ia,ifer,jfer,ib;
  uint16_t oadd[2];
  InitNamedOperator(H, "Impurity Hamiltonian");
  H->NBoson=0;
  H->NFermion=NFermionShift;
  if (Deflation==true)
  {
    for(i=1;i<=A->val.N;i++)
    {
      H->NFermion+=A->val.BlockSizeList[i];
    }
    for(i=1;i<=A->con.N;i++)
    {
      H->NFermion+=A->con.BlockSizeList[i];
    }
  }
  else
  {
    H->NFermion += A->val.N * A->BlockSize;
    H->NFermion += A->con.N * A->BlockSize;
  }
  if(Spin==true)
  {
    H->NFermion=2*(H->NFermion-NFermionShift)+NFermionShift;
  }
  if(NFermion!=0 )
  {
    if(NFermion>H->NFermion)
    {
      H->NFermion=NFermion;
    }
    else
    {
      printf("Warning in BlockNaturalImpurityOrbitalMatrixToImpurityHamiltonian\n");
      printf("Function called with NFermion=%i, but total number of Fermions must be larger than %i\n",NFermion,H->NFermion);
      printf("Number of Fermions enlarged\n");
    }
  }
  H->NWord=Min(H->NFermion,1+(H->NFermion-1)/8);
  InitOperatorFixedLength(H, 2, 8);
  if((A->val.N==0)&&(A->con.N==0))
  {
    return 0;
  }
  if(A->val.N==0)
  {
    FreeOperator(H);
    return BlockTridiagonalMatrixToImpurityHamiltonian(H, &(A->con), Spin, NFermionShift, NFermion);
  }
  if(A->con.N==0)
  {
    FreeOperator(H);
    return BlockTridiagonalMatrixToImpurityHamiltonian(H, &(A->val), Spin, NFermionShift, NFermion);
  }
  if((A->val.N==1)&&(A->con.N==1))
  {
    if(A->QComplex)
    {
      if(Spin==0)
      {
        for(i=0;i<2*A->BlockSize;i++)
        {
          for(j=0;j<2*A->BlockSize;j++)
          {
            oadd[1]=CrAnMask+NFermionShift+i;
            oadd[0]=NFermionShift+j;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[i*A->TMT.Nj+j].Re, A->TMT.V[i*A->TMT.Nj+j].Im, &QDetNotOrderedFcFaBcBa);
          }
        }
      }
      else
      {
        for(i=0;i<2*A->BlockSize;i++)
        {
          for(j=0;j<2*A->BlockSize;j++)
          {
            oadd[1]=CrAnMask+NFermionShift+2*i;
            oadd[0]=NFermionShift+2*j;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[i*A->TMT.Nj+j].Re, A->TMT.V[i*A->TMT.Nj+j].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+2*i+1;
            oadd[0]=NFermionShift+2*j+1;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[i*A->TMT.Nj+j].Re, A->TMT.V[i*A->TMT.Nj+j].Im, &QDetNotOrderedFcFaBcBa);
          }
        }
      }
    }
    else
    {
      if(Spin==0)
      {
        for(i=0;i<2*A->BlockSize;i++)
        {
          for(j=0;j<2*A->BlockSize;j++)
          {
            oadd[1]=CrAnMask+NFermionShift+i;
            oadd[0]=NFermionShift+j;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[i*A->TMT.Nj+j], &QDetNotOrderedFcFaBcBa);
          }
        }
      }
      else
      {
        for(i=0;i<2*A->BlockSize;i++)
        {
          for(j=0;j<2*A->BlockSize;j++)
          {
            oadd[1]=CrAnMask+NFermionShift+2*i;
            oadd[0]=NFermionShift+2*j;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[i*A->TMT.Nj+j], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+2*i+1;
            oadd[0]=NFermionShift+2*j+1;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[i*A->TMT.Nj+j], &QDetNotOrderedFcFaBcBa);
          }
        }
      }
    }
    return 0;
  }
  if(A->QComplex)
  {
    MakeRealOperatorFixedLengthComplex(H, 2);
    if(Spin==0)
    {
      /* transform matrix TMT to Operator */
      for(i=0;i<A->BlockSize;i++)
      {
        for(j=0;j<A->BlockSize;j++)
        {
          oadd[1]=CrAnMask+NFermionShift+i;
          oadd[0]=NFermionShift+j;
          ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[i*A->TMT.Nj+j].Re, A->TMT.V[i*A->TMT.Nj+j].Im, &QDetNotOrderedFcFaBcBa);
          oadd[1]=CrAnMask+NFermionShift+i+A->BlockSize;
          oadd[0]=NFermionShift+j;
          ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[(i+A->BlockSize)*A->TMT.Nj+j].Re, A->TMT.V[(i+A->BlockSize)*A->TMT.Nj+j].Im, &QDetNotOrderedFcFaBcBa);
          oadd[1]=CrAnMask+NFermionShift+i;
          oadd[0]=NFermionShift+j+A->BlockSize;
          ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[i*A->TMT.Nj+j+A->BlockSize].Re, A->TMT.V[i*A->TMT.Nj+j+A->BlockSize].Im, &QDetNotOrderedFcFaBcBa);
          oadd[1]=CrAnMask+NFermionShift+i+A->BlockSize;
          oadd[0]=NFermionShift+j+A->BlockSize;
          ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[(i+A->BlockSize)*A->TMT.Nj+j+A->BlockSize].Re, A->TMT.V[(i+A->BlockSize)*A->TMT.Nj+j+A->BlockSize].Im, &QDetNotOrderedFcFaBcBa);
          if(A->val.N>1)
          {
            oadd[1]=CrAnMask+NFermionShift+i+2*A->BlockSize;
            oadd[0]=NFermionShift+j;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[(i+2*A->BlockSize)*A->TMT.Nj+j].Re, A->TMT.V[(i+2*A->BlockSize)*A->TMT.Nj+j].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+i+2*A->BlockSize;
            oadd[0]=NFermionShift+j+A->BlockSize;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[(i+2*A->BlockSize)*A->TMT.Nj+j+A->BlockSize].Re, A->TMT.V[(i+2*A->BlockSize)*A->TMT.Nj+j+A->BlockSize].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+i;
            oadd[0]=NFermionShift+j+2*A->BlockSize;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[i*A->TMT.Nj+j+2*A->BlockSize].Re, A->TMT.V[i*A->TMT.Nj+j+2*A->BlockSize].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+i+A->BlockSize;
            oadd[0]=NFermionShift+j+2*A->BlockSize;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[(i+A->BlockSize)*A->TMT.Nj+j+2*A->BlockSize].Re, A->TMT.V[(i+A->BlockSize)*A->TMT.Nj+j+2*A->BlockSize].Im, &QDetNotOrderedFcFaBcBa);
          }
          if(A->con.N>1)
          {
            oadd[1]=CrAnMask+NFermionShift+i+(1+A->val.N)*A->BlockSize;
            oadd[0]=NFermionShift+j;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[(i+3*A->BlockSize)*A->TMT.Nj+j].Re, A->TMT.V[(i+3*A->BlockSize)*A->TMT.Nj+j].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+i+(1+A->val.N)*A->BlockSize;
            oadd[0]=NFermionShift+j+A->BlockSize;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[(i+3*A->BlockSize)*A->TMT.Nj+j+A->BlockSize].Re, A->TMT.V[(i+3*A->BlockSize)*A->TMT.Nj+j+A->BlockSize].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+i;
            oadd[0]=NFermionShift+j+(1+A->val.N)*A->BlockSize;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[i*A->TMT.Nj+j+3*A->BlockSize].Re, A->TMT.V[i*A->TMT.Nj+j+3*A->BlockSize].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+i+A->BlockSize;
            oadd[0]=NFermionShift+j+(1+A->val.N)*A->BlockSize;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[(i+A->BlockSize)*A->TMT.Nj+j+3*A->BlockSize].Re, A->TMT.V[(i+A->BlockSize)*A->TMT.Nj+j+3*A->BlockSize].Im, &QDetNotOrderedFcFaBcBa);
          }
        }
      }
      /* tridiagonal matrix val to operator */
      ia=2*A->BlockSize*A->BlockSize;
      for(k=2;k<=A->val.N;k++)
      {
        ifer=NFermionShift+k*A->BlockSize;
        for(i=0;i<A->BlockSize;i++)
        {
          jfer=NFermionShift+k*A->BlockSize;
          for(j=0;j<A->BlockSize;j++)
          {
            oadd[1]=CrAnMask+ifer;
            oadd[0]=jfer;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->val.a[ia].Re, A->val.a[ia].Im, &QDetNotOrderedFcFaBcBa);
            ia++;
            jfer++;
          }
          ifer++;
        }
      }
      ib=2*A->BlockSize*A->BlockSize;
      for(k=2;k<A->val.N;k++)
      {
        ifer=NFermionShift+(k+1)*A->BlockSize;
        for(i=0;i<A->BlockSize;i++)
        {
          jfer=NFermionShift+k*A->BlockSize;
          for(j=0;j<A->BlockSize;j++)
          {
            oadd[1]=CrAnMask+ifer;
            oadd[0]=jfer;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->val.b[ib].Re, A->val.b[ib].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+jfer;
            oadd[0]=ifer;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->val.b[ib].Re,-A->val.b[ib].Im, &QDetNotOrderedFcFaBcBa);
            ib++;
            jfer++;
          }
          ifer++;
        }
      }
      /* tridiagonal matrix con to operator */
      ia=2*A->BlockSize*A->BlockSize;
      for(k=2+A->val.N;k<=A->val.N+A->con.N;k++)
      {
        ifer=NFermionShift+(k-1)*A->BlockSize;
        for(i=0;i<A->BlockSize;i++)
        {
          jfer=NFermionShift+(k-1)*A->BlockSize;
          for(j=0;j<A->BlockSize;j++)
          {
            oadd[1]=CrAnMask+ifer;
            oadd[0]=jfer;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->con.a[ia].Re, A->con.a[ia].Im, &QDetNotOrderedFcFaBcBa);
            ia++;
            jfer++;
          }
          ifer++;
        }
      }
      ib=2*A->BlockSize*A->BlockSize;
      for(k=2+A->val.N;k<A->val.N+A->con.N;k++)
      {
        ifer=NFermionShift+(k+0)*A->BlockSize;
        for(i=0;i<A->BlockSize;i++)
        {
          jfer=NFermionShift+(k-1)*A->BlockSize;
          for(j=0;j<A->BlockSize;j++)
          {
            oadd[1]=CrAnMask+ifer;
            oadd[0]=jfer;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->con.b[ib].Re, A->con.b[ib].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+jfer;
            oadd[0]=ifer;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->con.b[ib].Re,-A->con.b[ib].Im, &QDetNotOrderedFcFaBcBa);
            ib++;
            jfer++;
          }
          ifer++;
        }
      }
    }
    else /* complex with spin */
    {
      /* transform matrix TMT to Operator */
      for(i=0;i<A->BlockSize;i++)
      {
        for(j=0;j<A->BlockSize;j++)
        {
          oadd[1]=CrAnMask+NFermionShift+2*i;
          oadd[0]=NFermionShift+2*j;
          ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[i*A->TMT.Nj+j].Re, A->TMT.V[i*A->TMT.Nj+j].Im, &QDetNotOrderedFcFaBcBa);
          oadd[1]=CrAnMask+NFermionShift+2*(i+A->BlockSize);
          oadd[0]=NFermionShift+2*j;
          ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[(i+A->BlockSize)*A->TMT.Nj+j].Re, A->TMT.V[(i+A->BlockSize)*A->TMT.Nj+j].Im, &QDetNotOrderedFcFaBcBa);
          oadd[1]=CrAnMask+NFermionShift+2*i;
          oadd[0]=NFermionShift+2*(j+A->BlockSize);
          ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[i*A->TMT.Nj+j+A->BlockSize].Re, A->TMT.V[i*A->TMT.Nj+j+A->BlockSize].Im, &QDetNotOrderedFcFaBcBa);
          oadd[1]=CrAnMask+NFermionShift+2*(i+A->BlockSize);
          oadd[0]=NFermionShift+2*(j+A->BlockSize);
          ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[(i+A->BlockSize)*A->TMT.Nj+j+A->BlockSize].Re, A->TMT.V[(i+A->BlockSize)*A->TMT.Nj+j+A->BlockSize].Im, &QDetNotOrderedFcFaBcBa);
          
          oadd[1]=CrAnMask+NFermionShift+2*i+1;
          oadd[0]=NFermionShift+2*j+1;
          ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[i*A->TMT.Nj+j].Re, A->TMT.V[i*A->TMT.Nj+j].Im, &QDetNotOrderedFcFaBcBa);
          oadd[1]=CrAnMask+NFermionShift+2*(i+A->BlockSize)+1;
          oadd[0]=NFermionShift+2*j+1;
          ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[(i+A->BlockSize)*A->TMT.Nj+j].Re, A->TMT.V[(i+A->BlockSize)*A->TMT.Nj+j].Im, &QDetNotOrderedFcFaBcBa);
          oadd[1]=CrAnMask+NFermionShift+2*i+1;
          oadd[0]=NFermionShift+2*(j+A->BlockSize)+1;
          ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[i*A->TMT.Nj+j+A->BlockSize].Re, A->TMT.V[i*A->TMT.Nj+j+A->BlockSize].Im, &QDetNotOrderedFcFaBcBa);
          oadd[1]=CrAnMask+NFermionShift+2*(i+A->BlockSize)+1;
          oadd[0]=NFermionShift+2*(j+A->BlockSize)+1;
          ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[(i+A->BlockSize)*A->TMT.Nj+j+A->BlockSize].Re, A->TMT.V[(i+A->BlockSize)*A->TMT.Nj+j+A->BlockSize].Im, &QDetNotOrderedFcFaBcBa);
          
          if(A->val.N>1)
          {
            oadd[1]=CrAnMask+NFermionShift+2*(i+2*A->BlockSize);
            oadd[0]=NFermionShift+2*j;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[(i+2*A->BlockSize)*A->TMT.Nj+j].Re, A->TMT.V[(i+2*A->BlockSize)*A->TMT.Nj+j].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+2*(i+2*A->BlockSize);
            oadd[0]=NFermionShift+2*(j+A->BlockSize);
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[(i+2*A->BlockSize)*A->TMT.Nj+j+A->BlockSize].Re, A->TMT.V[(i+2*A->BlockSize)*A->TMT.Nj+j+A->BlockSize].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+2*i;
            oadd[0]=NFermionShift+2*(j+2*A->BlockSize);
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[i*A->TMT.Nj+j+2*A->BlockSize].Re, A->TMT.V[i*A->TMT.Nj+j+2*A->BlockSize].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+2*(i+A->BlockSize);
            oadd[0]=NFermionShift+2*(j+2*A->BlockSize);
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[(i+A->BlockSize)*A->TMT.Nj+j+2*A->BlockSize].Re, A->TMT.V[(i+A->BlockSize)*A->TMT.Nj+j+2*A->BlockSize].Im, &QDetNotOrderedFcFaBcBa);
            
            oadd[1]=CrAnMask+NFermionShift+2*(i+2*A->BlockSize)+1;
            oadd[0]=NFermionShift+2*j+1;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[(i+2*A->BlockSize)*A->TMT.Nj+j].Re, A->TMT.V[(i+2*A->BlockSize)*A->TMT.Nj+j].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+2*(i+2*A->BlockSize)+1;
            oadd[0]=NFermionShift+2*(j+A->BlockSize)+1;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[(i+2*A->BlockSize)*A->TMT.Nj+j+A->BlockSize].Re, A->TMT.V[(i+2*A->BlockSize)*A->TMT.Nj+j+A->BlockSize].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+2*i+1;
            oadd[0]=NFermionShift+2*(j+2*A->BlockSize)+1;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[i*A->TMT.Nj+j+2*A->BlockSize].Re, A->TMT.V[i*A->TMT.Nj+j+2*A->BlockSize].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+2*(i+A->BlockSize)+1;
            oadd[0]=NFermionShift+2*(j+2*A->BlockSize)+1;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[(i+A->BlockSize)*A->TMT.Nj+j+2*A->BlockSize].Re, A->TMT.V[(i+A->BlockSize)*A->TMT.Nj+j+2*A->BlockSize].Im, &QDetNotOrderedFcFaBcBa);
          }
          if(A->con.N>1)
          {
            oadd[1]=CrAnMask+NFermionShift+2*(i+(1+A->val.N)*A->BlockSize);
            oadd[0]=NFermionShift+2*j;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[(i+3*A->BlockSize)*A->TMT.Nj+j].Re, A->TMT.V[(i+3*A->BlockSize)*A->TMT.Nj+j].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+2*(i+(1+A->val.N)*A->BlockSize);
            oadd[0]=NFermionShift+2*(j+A->BlockSize);
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[(i+3*A->BlockSize)*A->TMT.Nj+j+A->BlockSize].Re, A->TMT.V[(i+3*A->BlockSize)*A->TMT.Nj+j+A->BlockSize].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+2*i;
            oadd[0]=NFermionShift+2*(j+(1+A->val.N)*A->BlockSize);
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[i*A->TMT.Nj+j+3*A->BlockSize].Re, A->TMT.V[i*A->TMT.Nj+j+3*A->BlockSize].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+2*(i+A->BlockSize);
            oadd[0]=NFermionShift+2*(j+(1+A->val.N)*A->BlockSize);
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[(i+A->BlockSize)*A->TMT.Nj+j+3*A->BlockSize].Re, A->TMT.V[(i+A->BlockSize)*A->TMT.Nj+j+3*A->BlockSize].Im, &QDetNotOrderedFcFaBcBa);
            
            oadd[1]=CrAnMask+NFermionShift+2*(i+(1+A->val.N)*A->BlockSize)+1;
            oadd[0]=NFermionShift+2*j+1;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[(i+3*A->BlockSize)*A->TMT.Nj+j].Re, A->TMT.V[(i+3*A->BlockSize)*A->TMT.Nj+j].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+2*(i+(1+A->val.N)*A->BlockSize)+1;
            oadd[0]=NFermionShift+2*(j+A->BlockSize)+1;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[(i+3*A->BlockSize)*A->TMT.Nj+j+A->BlockSize].Re, A->TMT.V[(i+3*A->BlockSize)*A->TMT.Nj+j+A->BlockSize].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+2*i+1;
            oadd[0]=NFermionShift+2*(j+(1+A->val.N)*A->BlockSize)+1;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[i*A->TMT.Nj+j+3*A->BlockSize].Re, A->TMT.V[i*A->TMT.Nj+j+3*A->BlockSize].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+2*(i+A->BlockSize)+1;
            oadd[0]=NFermionShift+2*(j+(1+A->val.N)*A->BlockSize)+1;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[(i+A->BlockSize)*A->TMT.Nj+j+3*A->BlockSize].Re, A->TMT.V[(i+A->BlockSize)*A->TMT.Nj+j+3*A->BlockSize].Im, &QDetNotOrderedFcFaBcBa);
          }
        }
      }
      /* tridiagonal matrix val to operator */
      ia=2*A->BlockSize*A->BlockSize;
      for(k=2;k<=A->val.N;k++)
      {
        ifer=k*A->BlockSize;
        for(i=0;i<A->BlockSize;i++)
        {
          jfer=k*A->BlockSize;
          for(j=0;j<A->BlockSize;j++)
          {
            oadd[1]=CrAnMask+2*ifer+NFermionShift;
            oadd[0]=2*jfer+NFermionShift;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->val.a[ia].Re, A->val.a[ia].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+2*ifer+1+NFermionShift;
            oadd[0]=2*jfer+1+NFermionShift;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->val.a[ia].Re, A->val.a[ia].Im, &QDetNotOrderedFcFaBcBa);
            ia++;
            jfer++;
          }
          ifer++;
        }
      }
      ib=2*A->BlockSize*A->BlockSize;
      for(k=2;k<A->val.N;k++)
      {
        ifer=(k+1)*A->BlockSize;
        for(i=0;i<A->BlockSize;i++)
        {
          jfer=k*A->BlockSize;
          for(j=0;j<A->BlockSize;j++)
          {
            oadd[1]=CrAnMask+2*ifer+NFermionShift;
            oadd[0]=2*jfer+NFermionShift;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->val.b[ib].Re, A->val.b[ib].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+2*jfer+NFermionShift;
            oadd[0]=2*ifer+NFermionShift;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->val.b[ib].Re,-A->val.b[ib].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+2*ifer+1+NFermionShift;
            oadd[0]=2*jfer+1+NFermionShift;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->val.b[ib].Re, A->val.b[ib].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+2*jfer+1+NFermionShift;
            oadd[0]=2*ifer+1+NFermionShift;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->val.b[ib].Re,-A->val.b[ib].Im, &QDetNotOrderedFcFaBcBa);
            ib++;
            jfer++;
          }
          ifer++;
        }
      }
      /* tridiagonal matrix con to operator */
      ia=2*A->BlockSize*A->BlockSize;
      for(k=2+A->val.N;k<=A->val.N+A->con.N;k++)
      {
        ifer=(k-1)*A->BlockSize;
        for(i=0;i<A->BlockSize;i++)
        {
          jfer=(k-1)*A->BlockSize;
          for(j=0;j<A->BlockSize;j++)
          {
            oadd[1]=CrAnMask+2*ifer+NFermionShift;
            oadd[0]=2*jfer+NFermionShift;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->con.a[ia].Re, A->con.a[ia].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+2*ifer+1+NFermionShift;
            oadd[0]=2*jfer+1+NFermionShift;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->con.a[ia].Re, A->con.a[ia].Im, &QDetNotOrderedFcFaBcBa);
            ia++;
            jfer++;
          }
          ifer++;
        }
      }
      ib=2*A->BlockSize*A->BlockSize;
      for(k=2+A->val.N;k<A->val.N+A->con.N;k++)
      {
        ifer=(k+0)*A->BlockSize;
        for(i=0;i<A->BlockSize;i++)
        {
          jfer=(k-1)*A->BlockSize;
          for(j=0;j<A->BlockSize;j++)
          {
            oadd[1]=CrAnMask+2*ifer+NFermionShift;
            oadd[0]=2*jfer+NFermionShift;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->con.b[ib].Re, A->con.b[ib].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+2*jfer+NFermionShift;
            oadd[0]=2*ifer+NFermionShift;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->con.b[ib].Re,-A->con.b[ib].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+2*ifer+1+NFermionShift;
            oadd[0]=2*jfer+1+NFermionShift;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->con.b[ib].Re, A->con.b[ib].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+2*jfer+1+NFermionShift;
            oadd[0]=2*ifer+1+NFermionShift;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->con.b[ib].Re,-A->con.b[ib].Im, &QDetNotOrderedFcFaBcBa);
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
      /* transform matrix TMT to Operator */
      for(i=0;i<A->BlockSize;i++)
      {
        for(j=0;j<A->BlockSize;j++)
        {
          oadd[1]=CrAnMask+NFermionShift+i;
          oadd[0]=NFermionShift+j;
          RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[i*A->TMT.Nj+j], &QDetNotOrderedFcFaBcBa);
          oadd[1]=CrAnMask+NFermionShift+i+A->BlockSize;
          oadd[0]=NFermionShift+j;
          RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[(i+A->BlockSize)*A->TMT.Nj+j], &QDetNotOrderedFcFaBcBa);
          oadd[1]=CrAnMask+NFermionShift+i;
          oadd[0]=NFermionShift+j+A->BlockSize;
          RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[i*A->TMT.Nj+j+A->BlockSize], &QDetNotOrderedFcFaBcBa);
          oadd[1]=CrAnMask+NFermionShift+i+A->BlockSize;
          oadd[0]=NFermionShift+j+A->BlockSize;
          RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[(i+A->BlockSize)*A->TMT.Nj+j+A->BlockSize], &QDetNotOrderedFcFaBcBa);
          if(A->val.N>1)
          {
            oadd[1]=CrAnMask+NFermionShift+i+2*A->BlockSize;
            oadd[0]=NFermionShift+j;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[(i+2*A->BlockSize)*A->TMT.Nj+j], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+i+2*A->BlockSize;
            oadd[0]=NFermionShift+j+A->BlockSize;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[(i+2*A->BlockSize)*A->TMT.Nj+j+A->BlockSize], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+i;
            oadd[0]=NFermionShift+j+2*A->BlockSize;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[i*A->TMT.Nj+j+2*A->BlockSize], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+i+A->BlockSize;
            oadd[0]=NFermionShift+j+2*A->BlockSize;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[(i+A->BlockSize)*A->TMT.Nj+j+2*A->BlockSize], &QDetNotOrderedFcFaBcBa);
          }
          if(A->con.N>1)
          {
            oadd[1]=CrAnMask+NFermionShift+i+(1+A->val.N)*A->BlockSize;
            oadd[0]=NFermionShift+j;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[(i+3*A->BlockSize)*A->TMT.Nj+j], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+i+(1+A->val.N)*A->BlockSize;
            oadd[0]=NFermionShift+j+A->BlockSize;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[(i+3*A->BlockSize)*A->TMT.Nj+j+A->BlockSize], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+i;
            oadd[0]=NFermionShift+j+(1+A->val.N)*A->BlockSize;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[i*A->TMT.Nj+j+3*A->BlockSize], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+i+A->BlockSize;
            oadd[0]=NFermionShift+j+(1+A->val.N)*A->BlockSize;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[(i+A->BlockSize)*A->TMT.Nj+j+3*A->BlockSize], &QDetNotOrderedFcFaBcBa);
          }
        }
      }
      /* tridiagonal matrix val to operator */
      ia=2*A->BlockSize*A->BlockSize;
      for(k=2;k<=A->val.N;k++)
      {
        ifer=NFermionShift+k*A->BlockSize;
        for(i=0;i<A->BlockSize;i++)
        {
          jfer=NFermionShift+k*A->BlockSize;
          for(j=0;j<A->BlockSize;j++)
          {
            oadd[1]=CrAnMask+ifer;
            oadd[0]=jfer;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->val.aRe[ia], &QDetNotOrderedFcFaBcBa);
            ia++;
            jfer++;
          }
          ifer++;
        }
      }
      ib=2*A->BlockSize*A->BlockSize;
      for(k=2;k<A->val.N;k++)
      {
        ifer=NFermionShift+(k+1)*A->BlockSize;
        for(i=0;i<A->BlockSize;i++)
        {
          jfer=NFermionShift+k*A->BlockSize;
          for(j=0;j<A->BlockSize;j++)
          {
            oadd[1]=CrAnMask+ifer;
            oadd[0]=jfer;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->val.bRe[ib], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+jfer;
            oadd[0]=ifer;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->val.bRe[ib], &QDetNotOrderedFcFaBcBa);
            ib++;
            jfer++;
          }
          ifer++;
        }
      }
      /* tridiagonal matrix con to operator */
      ia=2*A->BlockSize*A->BlockSize;
      for(k=2+A->val.N;k<=A->val.N+A->con.N;k++)
      {
        ifer=NFermionShift+(k-1)*A->BlockSize;
        for(i=0;i<A->BlockSize;i++)
        {
          jfer=NFermionShift+(k-1)*A->BlockSize;
          for(j=0;j<A->BlockSize;j++)
          {
            oadd[1]=CrAnMask+ifer;
            oadd[0]=jfer;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->con.aRe[ia], &QDetNotOrderedFcFaBcBa);
            ia++;
            jfer++;
          }
          ifer++;
        }
      }
      ib=2*A->BlockSize*A->BlockSize;
      for(k=2+A->val.N;k<A->val.N+A->con.N;k++)
      {
        ifer=NFermionShift+(k+0)*A->BlockSize;
        for(i=0;i<A->BlockSize;i++)
        {
          jfer=NFermionShift+(k-1)*A->BlockSize;
          for(j=0;j<A->BlockSize;j++)
          {
            oadd[1]=CrAnMask+ifer;
            oadd[0]=jfer;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->con.bRe[ib], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+jfer;
            oadd[0]=ifer;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->con.bRe[ib], &QDetNotOrderedFcFaBcBa);
            ib++;
            jfer++;
          }
          ifer++;
        }
      }
    }
    else /* real with spin */
    {
      /* transform matrix TMT to Operator */
      for(i=0;i<A->BlockSize;i++)
      {
        for(j=0;j<A->BlockSize;j++)
        {
          oadd[1]=CrAnMask+NFermionShift+2*i;
          oadd[0]=NFermionShift+2*j;
          RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[i*A->TMT.Nj+j], &QDetNotOrderedFcFaBcBa);
          oadd[1]=CrAnMask+NFermionShift+2*(i+A->BlockSize);
          oadd[0]=NFermionShift+2*j;
          RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[(i+A->BlockSize)*A->TMT.Nj+j], &QDetNotOrderedFcFaBcBa);
          oadd[1]=CrAnMask+NFermionShift+2*i;
          oadd[0]=NFermionShift+2*(j+A->BlockSize);
          RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[i*A->TMT.Nj+j+A->BlockSize], &QDetNotOrderedFcFaBcBa);
          oadd[1]=CrAnMask+NFermionShift+2*(i+A->BlockSize);
          oadd[0]=NFermionShift+2*(j+A->BlockSize);
          RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[(i+A->BlockSize)*A->TMT.Nj+j+A->BlockSize], &QDetNotOrderedFcFaBcBa);
          
          oadd[1]=CrAnMask+NFermionShift+2*i+1;
          oadd[0]=NFermionShift+2*j+1;
          RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[i*A->TMT.Nj+j], &QDetNotOrderedFcFaBcBa);
          oadd[1]=CrAnMask+NFermionShift+2*(i+A->BlockSize)+1;
          oadd[0]=NFermionShift+2*j+1;
          RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[(i+A->BlockSize)*A->TMT.Nj+j], &QDetNotOrderedFcFaBcBa);
          oadd[1]=CrAnMask+NFermionShift+2*i+1;
          oadd[0]=NFermionShift+2*(j+A->BlockSize)+1;
          RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[i*A->TMT.Nj+j+A->BlockSize], &QDetNotOrderedFcFaBcBa);
          oadd[1]=CrAnMask+NFermionShift+2*(i+A->BlockSize)+1;
          oadd[0]=NFermionShift+2*(j+A->BlockSize)+1;
          RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[(i+A->BlockSize)*A->TMT.Nj+j+A->BlockSize], &QDetNotOrderedFcFaBcBa);
     
          if(A->val.N>1)
          {
            oadd[1]=CrAnMask+NFermionShift+2*(i+2*A->BlockSize);
            oadd[0]=NFermionShift+2*j;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[(i+2*A->BlockSize)*A->TMT.Nj+j], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+2*(i+2*A->BlockSize);
            oadd[0]=NFermionShift+2*(j+A->BlockSize);
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[(i+2*A->BlockSize)*A->TMT.Nj+j+A->BlockSize], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+2*i;
            oadd[0]=NFermionShift+2*(j+2*A->BlockSize);
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[i*A->TMT.Nj+j+2*A->BlockSize], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+2*(i+A->BlockSize);
            oadd[0]=NFermionShift+2*(j+2*A->BlockSize);
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[(i+A->BlockSize)*A->TMT.Nj+j+2*A->BlockSize], &QDetNotOrderedFcFaBcBa);
            
            oadd[1]=CrAnMask+NFermionShift+2*(i+2*A->BlockSize)+1;
            oadd[0]=NFermionShift+2*j+1;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[(i+2*A->BlockSize)*A->TMT.Nj+j], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+2*(i+2*A->BlockSize)+1;
            oadd[0]=NFermionShift+2*(j+A->BlockSize)+1;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[(i+2*A->BlockSize)*A->TMT.Nj+j+A->BlockSize], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+2*i+1;
            oadd[0]=NFermionShift+2*(j+2*A->BlockSize)+1;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[i*A->TMT.Nj+j+2*A->BlockSize], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+2*(i+A->BlockSize)+1;
            oadd[0]=NFermionShift+2*(j+2*A->BlockSize)+1;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[(i+A->BlockSize)*A->TMT.Nj+j+2*A->BlockSize], &QDetNotOrderedFcFaBcBa);
          }
          if(A->con.N>1)
          {
            oadd[1]=CrAnMask+NFermionShift+2*(i+(1+A->val.N)*A->BlockSize);
            oadd[0]=NFermionShift+2*j;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[(i+3*A->BlockSize)*A->TMT.Nj+j], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+2*(i+(1+A->val.N)*A->BlockSize);
            oadd[0]=NFermionShift+2*(j+A->BlockSize);
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[(i+3*A->BlockSize)*A->TMT.Nj+j+A->BlockSize], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+2*i;
            oadd[0]=NFermionShift+2*(j+(1+A->val.N)*A->BlockSize);
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[i*A->TMT.Nj+j+3*A->BlockSize], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+2*(i+A->BlockSize);
            oadd[0]=NFermionShift+2*(j+(1+A->val.N)*A->BlockSize);
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[(i+A->BlockSize)*A->TMT.Nj+j+3*A->BlockSize], &QDetNotOrderedFcFaBcBa);
            
            oadd[1]=CrAnMask+NFermionShift+2*(i+(1+A->val.N)*A->BlockSize)+1;
            oadd[0]=NFermionShift+2*j+1;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[(i+3*A->BlockSize)*A->TMT.Nj+j], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+2*(i+(1+A->val.N)*A->BlockSize)+1;
            oadd[0]=NFermionShift+2*(j+A->BlockSize)+1;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[(i+3*A->BlockSize)*A->TMT.Nj+j+A->BlockSize], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+2*i+1;
            oadd[0]=NFermionShift+2*(j+(1+A->val.N)*A->BlockSize)+1;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[i*A->TMT.Nj+j+3*A->BlockSize], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+2*(i+A->BlockSize)+1;
            oadd[0]=NFermionShift+2*(j+(1+A->val.N)*A->BlockSize)+1;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[(i+A->BlockSize)*A->TMT.Nj+j+3*A->BlockSize], &QDetNotOrderedFcFaBcBa);
            
          }
        }
      }
      /* tridiagonal matrix val to operator */
      ia=2*A->BlockSize*A->BlockSize;
      for(k=2;k<=A->val.N;k++)
      {
        ifer=k*A->BlockSize;
        for(i=0;i<A->BlockSize;i++)
        {
          jfer=k*A->BlockSize;
          for(j=0;j<A->BlockSize;j++)
          {
            oadd[1]=CrAnMask+2*ifer+NFermionShift;
            oadd[0]=2*jfer+NFermionShift;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->val.aRe[ia], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+2*ifer+1+NFermionShift;
            oadd[0]=2*jfer+1+NFermionShift;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->val.aRe[ia], &QDetNotOrderedFcFaBcBa);
            ia++;
            jfer++;
          }
          ifer++;
        }
      }
      ib=2*A->BlockSize*A->BlockSize;
      for(k=2;k<A->val.N;k++)
      {
        ifer=(k+1)*A->BlockSize;
        for(i=0;i<A->BlockSize;i++)
        {
          jfer=k*A->BlockSize;
          for(j=0;j<A->BlockSize;j++)
          {
            oadd[1]=CrAnMask+2*ifer+NFermionShift;
            oadd[0]=2*jfer+NFermionShift;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->val.bRe[ib], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+2*jfer+NFermionShift;
            oadd[0]=2*ifer+NFermionShift;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->val.bRe[ib], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+2*ifer+1+NFermionShift;
            oadd[0]=2*jfer+1+NFermionShift;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->val.bRe[ib], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+2*jfer+1+NFermionShift;
            oadd[0]=2*ifer+1+NFermionShift;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->val.bRe[ib], &QDetNotOrderedFcFaBcBa);
            ib++;
            jfer++;
          }
          ifer++;
        }
      }
      /* tridiagonal matrix con to operator */
      ia=2*A->BlockSize*A->BlockSize;
      for(k=2+A->val.N;k<=A->val.N+A->con.N;k++)
      {
        ifer=(k-1)*A->BlockSize;
        for(i=0;i<A->BlockSize;i++)
        {
          jfer=(k-1)*A->BlockSize;
          for(j=0;j<A->BlockSize;j++)
          {
            oadd[1]=CrAnMask+2*ifer+NFermionShift;
            oadd[0]=2*jfer+NFermionShift;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->con.aRe[ia], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+2*ifer+1+NFermionShift;
            oadd[0]=2*jfer+1+NFermionShift;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->con.aRe[ia], &QDetNotOrderedFcFaBcBa);
            ia++;
            jfer++;
          }
          ifer++;
        }
      }
      ib=2*A->BlockSize*A->BlockSize;
      for(k=2+A->val.N;k<A->val.N+A->con.N;k++)
      {
        ifer=(k+0)*A->BlockSize;
        for(i=0;i<A->BlockSize;i++)
        {
          jfer=(k-1)*A->BlockSize;
          for(j=0;j<A->BlockSize;j++)
          {
            oadd[1]=CrAnMask+2*ifer+NFermionShift;
            oadd[0]=2*jfer+NFermionShift;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->con.bRe[ib], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+2*jfer+NFermionShift;
            oadd[0]=2*ifer+NFermionShift;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->con.bRe[ib], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+2*ifer+1+NFermionShift;
            oadd[0]=2*jfer+1+NFermionShift;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->con.bRe[ib], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+2*jfer+1+NFermionShift;
            oadd[0]=2*ifer+1+NFermionShift;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->con.bRe[ib], &QDetNotOrderedFcFaBcBa);
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


uint32_t BlockNaturalImpurityOrbitalMatrixToImpurityHamiltonianPerturpation(OperatorType *const H, OperatorType * const VPerturbation, TriDiagonalMatrixType * const gPerturbation, BlockNaturalImpurityOrbitalMatrixType const *A, uint32_t const NBlockMaxVal, uint32_t const NBlockMaxCon, uint32_t const Spin, uint32_t const NFermionShift, uint32_t const NFermion)
{
  //TODO: implement
  printf("To DO\n");
  return 0;
}

/*  Add one block of a chain to an impurity Hamiltonian.
 *
 *    aNotB == true    the diagonal block C->a[iblock], added at orbitals (ip,ip). It is
 *                     hermitian and covers both of its own triangles, so nothing else
 *                     is added for it.
 *    aNotB == false   the coupling C->b[iblock], added from orbitals ip to iq together
 *                     with its hermitian conjugate from iq back to ip.
 *
 *  Conj adds the conjugate of the block instead of the block. That is what the head
 *  coupling b[0] of a chain needs and what a coupling inside a chain must not have: b[0]
 *  enters the Green's function as b0^* . G . b0^T where the couplings below it enter as
 *  b . G . b^+ (see BlockTriDiagonalMatrixType in StructDef.cpp), and a hermitian matrix
 *  with X in the upper right block gives X . G . X^+, so X is conj(b0) at the head and b
 *  everywhere else. BlockDoubleTriDiagonalMatrixToCompactMatrixDeflate builds the same
 *  matrix and says the same thing.
 *
 *  ip and iq count orbitals from the start of the impurity problem, not fermions. With
 *  Spin each orbital becomes two neighbouring fermions and the block is added once for
 *  each, so orbital n spin s is fermion NFermionShift + 2*n + s. The shift itself is not
 *  doubled: the core orbitals in front of it are not part of the impurity problem. That is
 *  what the single valued family does and what CreateBitMasks reports.
 *
 *  QComplex says whether H holds complex numbers, which it does whenever either chain
 *  does; the block itself is read according to its own chain, since a mixed matrix can
 *  have one real chain and one complex one.
 */
static void BlockDoubleTriAddBlockToH(OperatorType *const H, BlockTriDiagonalMatrixType const * const C, bool const aNotB, uint32_t const iblock, uint32_t const ip, uint32_t const iq, uint32_t const NFermionShift, bool const Spin, bool const QComplex, bool const Conj)
{
  uint32_t i,j,s,ns,bs,base,fp,fq;
  uint16_t oadd[2];
  double re,im;
  bs   = C->BlockSize;
  base = iblock * bs * bs;
  ns   = (Spin==true) ? 2 : 1;
  for(i=0;i<bs;i++)
  {
    for(j=0;j<bs;j++)
    {
      if(C->QComplex)
      {
        if(aNotB==true) { re = C->a[base+i*bs+j].Re; im = C->a[base+i*bs+j].Im; }
        else            { re = C->b[base+i*bs+j].Re; im = C->b[base+i*bs+j].Im; }
      }
      else
      {
        if(aNotB==true) re = C->aRe[base+i*bs+j];
        else            re = C->bRe[base+i*bs+j];
        im = 0;
      }
      if(Conj==true) im = -im;
      for(s=0;s<ns;s++)
      {
        fp = NFermionShift + ns*(ip+i) + s;
        fq = NFermionShift + ns*(iq+j) + s;
        oadd[1] = CrAnMask + fp;
        oadd[0] = fq;
        if(QComplex) ComplexOperatorAddLaderNormalOrder(H, oadd, 2, re, im, &QDetNotOrderedFcFaBcBa);
        else         RealOperatorAddLaderNormalOrder(H, oadd, 2, re, &QDetNotOrderedFcFaBcBa);
        if(aNotB==false)
        {
          oadd[1] = CrAnMask + fq;
          oadd[0] = fp;
          if(QComplex) ComplexOperatorAddLaderNormalOrder(H, oadd, 2, re, -im, &QDetNotOrderedFcFaBcBa);
          else         RealOperatorAddLaderNormalOrder(H, oadd, 2, re, &QDetNotOrderedFcFaBcBa);
        }
      }
    }
  }
}

/*  One impurity block with a valence and a conduction chain hanging off it, as an operator.
 *
 *  The orbitals are
 *
 *      [ core | impurity | val_1 ... val_Nv | con_1 ... con_Nc ]
 *
 *  one block each, which is the layout ResponseFunction.CreateBitMasks reports and the one
 *  BlockDoubleTriDiagonalMatrixToCompactMatrixDeflate builds. a[0] of a chain is the
 *  impurity block -- the two chains share that block, so it is val.a[0] + con.a[0] -- and
 *  b[0] is the hopping from the impurity into the first block of that chain, with a[1] to
 *  a[N] and b[1] to b[N-1] the chain itself. So a chain contributes N blocks here, one more
 *  per chain than the natural orbital form puts in its chains, where a[1] is the natural
 *  orbital itself and belongs to the impurity side.
 *
 *  A->a0 and A->b0 are the offset and the scaling of the Green's function,
 *  G = A0 + B0^* [1/(w-H)]_00 B0^T, and are not part of H, as in every other
 *  …ToImpurityHamiltonian. Read them with BlockResponseFunctionGetA0 and ...GetB0.
 *
 *  An empty chain is a normal case for this type rather than an edge case -- it is what a
 *  chemical potential outside the bath gives -- and needs none of its own code: the loops
 *  over that chain do not run and its a[0] still contributes its share of the impurity
 *  energy. Note this cannot be handed to BlockTridiagonalMatrixToImpurityHamiltonian when
 *  one chain is empty, the way the natural orbital form can, because that builds a[1] to
 *  a[N] and treats a[0] and b[0] as the offset and the scaling: it would drop the impurity
 *  block and lose the other chain's share of its energy.
 *
 *  Deflation is not implemented, as in the rest of this family; every block is BlockSize
 *  wide. The chains of a d can carry a BlockSizeList, so this is the place to start when
 *  deflation is wanted here.
 */
uint32_t BlockDoubleTriDiagonalMatrixToImpurityHamiltonian(OperatorType *const H, BlockDoubleTriDiagonalMatrixType const *A, bool const Spin, uint32_t const NFermionShift, uint32_t const NFermion)
{
  uint32_t bs, k, NBlock, Nv, Nc;
  bool QComplex;
  bs = A->BlockSize;
  Nv = A->val.N;
  Nc = A->con.N;
  NBlock = 1 + Nv + Nc;
  /*  an operator has no mixed representation, so it is complex when either chain is. a0 and
   *  b0 are not part of H, so their being complex does not make H complex. */
  QComplex = ((A->val.QComplex!=0)||(A->con.QComplex!=0));
  InitNamedOperator(H, "Impurity Hamiltonian");
  H->NBoson=0;
  if(Spin==false) H->NFermion = NFermionShift + NBlock*bs;
  else            H->NFermion = NFermionShift + 2*NBlock*bs;
  if(NFermion!=0)
  {
    if(NFermion>H->NFermion)
    {
      H->NFermion=NFermion;
    }
    else
    {
      printf("Warning in BlockDoubleTriDiagonalMatrixToImpurityHamiltonian\n");
      printf("Function called with NFermion=%i, but total number of Fermions must be larger than %i\n",NFermion,H->NFermion);
      printf("Number of Fermions enlarged\n");
    }
  }
  H->NWord=Min(H->NFermion,1+(H->NFermion-1)/8);
  InitOperatorFixedLength(H, 2, 8);
  if(QComplex==true) MakeRealOperatorFixedLengthComplex(H, 2);
  /*  the impurity block, val.a[0] + con.a[0]. The two shares are added as two blocks;
   *  operator terms accumulate, and this way each is read with its own chain's QComplex */
  BlockDoubleTriAddBlockToH(H, &(A->val), true, 0, 0, 0, NFermionShift, Spin, QComplex, false);
  BlockDoubleTriAddBlockToH(H, &(A->con), true, 0, 0, 0, NFermionShift, Spin, QComplex, false);
  /*  the site energies of the two chains, a[k] of block k */
  for(k=1;k<=Nv;k++)
  {
    BlockDoubleTriAddBlockToH(H, &(A->val), true, k, k*bs, k*bs, NFermionShift, Spin, QComplex, false);
  }
  for(k=1;k<=Nc;k++)
  {
    BlockDoubleTriAddBlockToH(H, &(A->con), true, k, (Nv+k)*bs, (Nv+k)*bs, NFermionShift, Spin, QComplex, false);
  }
  /*  the impurity to first block hoppings, which are the two heads and take the conjugate */
  if(Nv!=0)
  {
    BlockDoubleTriAddBlockToH(H, &(A->val), false, 0, 0, bs, NFermionShift, Spin, QComplex, true);
  }
  if(Nc!=0)
  {
    BlockDoubleTriAddBlockToH(H, &(A->con), false, 0, 0, (Nv+1)*bs, NFermionShift, Spin, QComplex, true);
  }
  /*  and the hoppings inside the chains, b[k] between block k and block k+1, as stored */
  for(k=1;k<Nv;k++)
  {
    BlockDoubleTriAddBlockToH(H, &(A->val), false, k, k*bs, (k+1)*bs, NFermionShift, Spin, QComplex, false);
  }
  for(k=1;k<Nc;k++)
  {
    BlockDoubleTriAddBlockToH(H, &(A->con), false, k, (Nv+k)*bs, (Nv+k+1)*bs, NFermionShift, Spin, QComplex, false);
  }
  return 0;
}
