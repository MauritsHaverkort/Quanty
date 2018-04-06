#include "BasicMath_BlockMatrixToOperator.h"

/*  Starting from a BlockTridiagonalMatrix creates an impurity Hamiltonian
 *  H on input not initiallized
 *  NBlockMax the maximal number of blocks in the Hamiltonian, this includes the impurity.
 *   NBlockMax=2 would thus correspond to a triditional ligand field calculation
 *  if Spin=0 each block in A is taken once into the Hamiltonian
 *  If Spin =1 (!0) each block is copied twice, once for spin up, once for spin down. The order is block up, block down
 *  NFermionShift is the numbers of Fermions shifted. (i.e. if you need a 2p core before a d impurity block NFermionShift=6 and the p-shell is from 0 to 5, the d-shell form 6 to 15.
 *  NFermion is the total number of Fermions in the basis. 
 *  If NFermion=0 on input the number will be calculated such that H fits.
 */
uint32_t BlockTridiagonalMatrixToImpurityHamiltonian(OperatorType *const H, BlockTriDiagonalMatrixType const *A, uint32_t const NBlockMax, uint32_t const Spin, uint32_t const NFermionShift, uint32_t const NFermion)
{
  //TODO: test add H to Block tridiag using block lanczos
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
        ifer=NFermionShift+2*k*A->BlockSize;
        for(i=0;i<A->BlockSize;i++)
        {
          jfer=NFermionShift+2*(k+1)*A->BlockSize;
          for(j=0;j<A->BlockSize;j++)
          {
            oadd[1]=CrAnMask+ifer;
            oadd[0]=jfer;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->b[ib].Re,-A->b[ib].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+jfer;
            oadd[0]=ifer;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->b[ib].Re, A->b[ib].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+ifer+A->BlockSize;
            oadd[0]=jfer+A->BlockSize;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->b[ib].Re,-A->b[ib].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+jfer+A->BlockSize;
            oadd[0]=ifer+A->BlockSize;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->b[ib].Re, A->b[ib].Im, &QDetNotOrderedFcFaBcBa);
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
        ifer=NFermionShift+2*k*A->BlockSize;
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
  //TODO: test
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
uint32_t BlockNaturalImpurityOrbitalMatrixToImpurityHamiltonian(OperatorType *const H, BlockNaturalImpurityOrbitalMatrixType const *A, uint32_t const NBlockMaxVal, uint32_t const NBlockMaxCon, uint32_t const Spin, uint32_t const NFermionShift, uint32_t const NFermion)
{
  //TODO: test
  uint32_t i,j,k,ia,ifer,jfer,ib;
  uint16_t oadd[2];
  InitNamedOperator(H, "Impurity Hamiltonian");
  H->NBoson=0;
  if(Spin==0)
  {
    H->NFermion=NFermionShift+(2+NBlockMaxVal+NBlockMaxCon)*A->BlockSize;
  }
  else
  {
    H->NFermion=NFermionShift+2*(2+NBlockMaxVal+NBlockMaxCon)*A->BlockSize;
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
  InitOperatorFixedLength(H, 4, 8);
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
          if(NBlockMaxVal!=0)
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
          if(NBlockMaxCon!=0)
          {
            oadd[1]=CrAnMask+NFermionShift+i+(2+NBlockMaxVal)*A->BlockSize;
            oadd[0]=NFermionShift+j;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[(i+3*A->BlockSize)*A->TMT.Nj+j].Re, A->TMT.V[(i+3*A->BlockSize)*A->TMT.Nj+j].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+i+(2+NBlockMaxVal)*A->BlockSize;
            oadd[0]=NFermionShift+j+A->BlockSize;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[(i+3*A->BlockSize)*A->TMT.Nj+j+A->BlockSize].Re, A->TMT.V[(i+3*A->BlockSize)*A->TMT.Nj+j+A->BlockSize].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+i;
            oadd[0]=NFermionShift+j+(2+NBlockMaxVal)*A->BlockSize;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[i*A->TMT.Nj+j+3*A->BlockSize].Re, A->TMT.V[i*A->TMT.Nj+j+3*A->BlockSize].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+i+A->BlockSize;
            oadd[0]=NFermionShift+j+(2+NBlockMaxVal)*A->BlockSize;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[(i+A->BlockSize)*A->TMT.Nj+j+3*A->BlockSize].Re, A->TMT.V[(i+A->BlockSize)*A->TMT.Nj+j+3*A->BlockSize].Im, &QDetNotOrderedFcFaBcBa);
          }
        }
      }
      /* tridiagonal matrix val to operator */
      ia=2*A->BlockSize*A->BlockSize;
      for(k=2;k<=Min(NBlockMaxVal+1,A->val.N);k++)
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
      for(k=2;k<Min(NBlockMaxVal+1,A->val.N);k++)
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
      for(k=2+NBlockMaxVal;k<=NBlockMaxVal+Min(NBlockMaxCon+1,A->con.N);k++)
      {
        ifer=NFermionShift+k*A->BlockSize;
        for(i=0;i<A->BlockSize;i++)
        {
          jfer=NFermionShift+k*A->BlockSize;
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
      for(k=2+NBlockMaxVal;k<NBlockMaxVal+Min(NBlockMaxCon+1,A->con.N);k++)
      {
        ifer=NFermionShift+(k+1)*A->BlockSize;
        for(i=0;i<A->BlockSize;i++)
        {
          jfer=NFermionShift+k*A->BlockSize;
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
      /* transform matrix TMT to Operator spin up*/
      for(i=0;i<A->BlockSize;i++)
      {
        for(j=0;j<A->BlockSize;j++)
        {
          oadd[1]=CrAnMask+NFermionShift+i;
          oadd[0]=NFermionShift+j;
          ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[i*A->TMT.Nj+j].Re, A->TMT.V[i*A->TMT.Nj+j].Im, &QDetNotOrderedFcFaBcBa);
          oadd[1]=CrAnMask+NFermionShift+i+2*A->BlockSize;
          oadd[0]=NFermionShift+j;
          ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[(i+A->BlockSize)*A->TMT.Nj+j].Re, A->TMT.V[(i+A->BlockSize)*A->TMT.Nj+j].Im, &QDetNotOrderedFcFaBcBa);
          oadd[1]=CrAnMask+NFermionShift+i;
          oadd[0]=NFermionShift+j+2*A->BlockSize;
          ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[i*A->TMT.Nj+j+A->BlockSize].Re, A->TMT.V[i*A->TMT.Nj+j+A->BlockSize].Im, &QDetNotOrderedFcFaBcBa);
          oadd[1]=CrAnMask+NFermionShift+i+2*A->BlockSize;
          oadd[0]=NFermionShift+j+2*A->BlockSize;
          ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[(i+A->BlockSize)*A->TMT.Nj+j+A->BlockSize].Re, A->TMT.V[(i+A->BlockSize)*A->TMT.Nj+j+A->BlockSize].Im, &QDetNotOrderedFcFaBcBa);
          if(NBlockMaxVal!=0)
          {
            oadd[1]=CrAnMask+NFermionShift+i+4*A->BlockSize;
            oadd[0]=NFermionShift+j;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[(i+2*A->BlockSize)*A->TMT.Nj+j].Re, A->TMT.V[(i+2*A->BlockSize)*A->TMT.Nj+j].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+i+4*A->BlockSize;
            oadd[0]=NFermionShift+j+2*A->BlockSize;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[(i+2*A->BlockSize)*A->TMT.Nj+j+A->BlockSize].Re, A->TMT.V[(i+2*A->BlockSize)*A->TMT.Nj+j+A->BlockSize].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+i;
            oadd[0]=NFermionShift+j+4*A->BlockSize;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[i*A->TMT.Nj+j+2*A->BlockSize].Re, A->TMT.V[i*A->TMT.Nj+j+2*A->BlockSize].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+i+2*A->BlockSize;
            oadd[0]=NFermionShift+j+4*A->BlockSize;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[(i+A->BlockSize)*A->TMT.Nj+j+2*A->BlockSize].Re, A->TMT.V[(i+A->BlockSize)*A->TMT.Nj+j+2*A->BlockSize].Im, &QDetNotOrderedFcFaBcBa);
          }
          if(NBlockMaxCon!=0)
          {
            oadd[1]=CrAnMask+NFermionShift+i+2*(2+NBlockMaxVal)*A->BlockSize;
            oadd[0]=NFermionShift+j;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[(i+3*A->BlockSize)*A->TMT.Nj+j].Re, A->TMT.V[(i+3*A->BlockSize)*A->TMT.Nj+j].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+i+2*(2+NBlockMaxVal)*A->BlockSize;
            oadd[0]=NFermionShift+j+2*A->BlockSize;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[(i+3*A->BlockSize)*A->TMT.Nj+j+A->BlockSize].Re, A->TMT.V[(i+3*A->BlockSize)*A->TMT.Nj+j+A->BlockSize].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+i;
            oadd[0]=NFermionShift+j+2*(2+NBlockMaxVal)*A->BlockSize;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[i*A->TMT.Nj+j+3*A->BlockSize].Re, A->TMT.V[i*A->TMT.Nj+j+3*A->BlockSize].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+i+2*A->BlockSize;
            oadd[0]=NFermionShift+j+2*(2+NBlockMaxVal)*A->BlockSize;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[(i+A->BlockSize)*A->TMT.Nj+j+3*A->BlockSize].Re, A->TMT.V[(i+A->BlockSize)*A->TMT.Nj+j+3*A->BlockSize].Im, &QDetNotOrderedFcFaBcBa);
          }
        }
      }
      /* tridiagonal matrix val to operator spin up*/
      ia=2*A->BlockSize*A->BlockSize;
      for(k=2;k<=Min(NBlockMaxVal+1,A->val.N);k++)
      {
        ifer=NFermionShift+2*k*A->BlockSize;
        for(i=0;i<A->BlockSize;i++)
        {
          jfer=NFermionShift+2*k*A->BlockSize;
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
      for(k=2;k<Min(NBlockMaxVal+1,A->val.N);k++)
      {
        ifer=NFermionShift+2*(k+1)*A->BlockSize;
        for(i=0;i<A->BlockSize;i++)
        {
          jfer=NFermionShift+2*k*A->BlockSize;
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
      /* tridiagonal matrix con to operator spin up*/
      ia=2*A->BlockSize*A->BlockSize;
      for(k=2+NBlockMaxVal;k<=NBlockMaxVal+Min(NBlockMaxCon+1,A->con.N);k++)
      {
        ifer=NFermionShift+2*k*A->BlockSize;
        for(i=0;i<A->BlockSize;i++)
        {
          jfer=NFermionShift+2*k*A->BlockSize;
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
      for(k=2+NBlockMaxVal;k<NBlockMaxVal+Min(NBlockMaxCon+1,A->con.N);k++)
      {
        ifer=NFermionShift+2*(k+1)*A->BlockSize;
        for(i=0;i<A->BlockSize;i++)
        {
          jfer=NFermionShift+2*k*A->BlockSize;
          for(j=0;j<A->BlockSize;j++)
          {
            oadd[1]=CrAnMask+ifer;
            oadd[0]=jfer;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->con.b[ib].Re, A->con.b[ib].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+jfer;
            oadd[0]=ifer;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->con.b[ib].Re, A->con.b[ib].Im, &QDetNotOrderedFcFaBcBa);
            ib++;
            jfer++;
          }
          ifer++;
        }
      }
      /* transform matrix TMT to Operator spin dn*/
      for(i=0;i<A->BlockSize;i++)
      {
        for(j=0;j<A->BlockSize;j++)
        {
          oadd[1]=CrAnMask+NFermionShift+i+A->BlockSize;
          oadd[0]=NFermionShift+j+A->BlockSize;
          ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[i*A->TMT.Nj+j].Re, A->TMT.V[i*A->TMT.Nj+j].Im, &QDetNotOrderedFcFaBcBa);
          oadd[1]=CrAnMask+NFermionShift+i+3*A->BlockSize;
          oadd[0]=NFermionShift+j+A->BlockSize;
          ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[(i+A->BlockSize)*A->TMT.Nj+j].Re, A->TMT.V[(i+A->BlockSize)*A->TMT.Nj+j].Im, &QDetNotOrderedFcFaBcBa);
          oadd[1]=CrAnMask+NFermionShift+i+A->BlockSize;
          oadd[0]=NFermionShift+j+3*A->BlockSize;
          ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[i*A->TMT.Nj+j+A->BlockSize].Re, A->TMT.V[i*A->TMT.Nj+j+A->BlockSize].Im, &QDetNotOrderedFcFaBcBa);
          oadd[1]=CrAnMask+NFermionShift+i+3*A->BlockSize;
          oadd[0]=NFermionShift+j+3*A->BlockSize;
          ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[(i+A->BlockSize)*A->TMT.Nj+j+A->BlockSize].Re, A->TMT.V[(i+A->BlockSize)*A->TMT.Nj+j+A->BlockSize].Im, &QDetNotOrderedFcFaBcBa);
          if(NBlockMaxVal!=0)
          {
            oadd[1]=CrAnMask+NFermionShift+i+5*A->BlockSize;
            oadd[0]=NFermionShift+j+A->BlockSize;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[(i+2*A->BlockSize)*A->TMT.Nj+j].Re, A->TMT.V[(i+2*A->BlockSize)*A->TMT.Nj+j].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+i+5*A->BlockSize;
            oadd[0]=NFermionShift+j+3*A->BlockSize;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[(i+2*A->BlockSize)*A->TMT.Nj+j+A->BlockSize].Re, A->TMT.V[(i+2*A->BlockSize)*A->TMT.Nj+j+A->BlockSize].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+i+A->BlockSize;
            oadd[0]=NFermionShift+j+5*A->BlockSize;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[i*A->TMT.Nj+j+2*A->BlockSize].Re, A->TMT.V[i*A->TMT.Nj+j+2*A->BlockSize].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+i+3*A->BlockSize;
            oadd[0]=NFermionShift+j+5*A->BlockSize;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[(i+A->BlockSize)*A->TMT.Nj+j+2*A->BlockSize].Re, A->TMT.V[(i+A->BlockSize)*A->TMT.Nj+j+2*A->BlockSize].Im, &QDetNotOrderedFcFaBcBa);
          }
          if(NBlockMaxCon!=0)
          {
            oadd[1]=CrAnMask+NFermionShift+i+2*(2+NBlockMaxVal)*A->BlockSize+A->BlockSize;
            oadd[0]=NFermionShift+j+A->BlockSize;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[(i+3*A->BlockSize)*A->TMT.Nj+j].Re, A->TMT.V[(i+3*A->BlockSize)*A->TMT.Nj+j].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+i+2*(2+NBlockMaxVal)*A->BlockSize+A->BlockSize;
            oadd[0]=NFermionShift+j+3*A->BlockSize;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[(i+3*A->BlockSize)*A->TMT.Nj+j+A->BlockSize].Re, A->TMT.V[(i+3*A->BlockSize)*A->TMT.Nj+j+A->BlockSize].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+i+A->BlockSize;
            oadd[0]=NFermionShift+j+2*(2+NBlockMaxVal)*A->BlockSize+A->BlockSize;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[i*A->TMT.Nj+j+3*A->BlockSize].Re, A->TMT.V[i*A->TMT.Nj+j+3*A->BlockSize].Im, &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+i+3*A->BlockSize;
            oadd[0]=NFermionShift+j+2*(2+NBlockMaxVal)*A->BlockSize+A->BlockSize;
            ComplexOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.V[(i+A->BlockSize)*A->TMT.Nj+j+3*A->BlockSize].Re, A->TMT.V[(i+A->BlockSize)*A->TMT.Nj+j+3*A->BlockSize].Im, &QDetNotOrderedFcFaBcBa);
          }
        }
      }
      /* tridiagonal matrix val to operator spin dn*/
      ia=2*A->BlockSize*A->BlockSize;
      for(k=2;k<=Min(NBlockMaxVal+1,A->val.N);k++)
      {
        ifer=NFermionShift+2*k*A->BlockSize+A->BlockSize;
        for(i=0;i<A->BlockSize;i++)
        {
          jfer=NFermionShift+2*k*A->BlockSize+A->BlockSize;
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
      for(k=2;k<Min(NBlockMaxVal+1,A->val.N);k++)
      {
        ifer=NFermionShift+2*(k+1)*A->BlockSize+A->BlockSize;
        for(i=0;i<A->BlockSize;i++)
        {
          jfer=NFermionShift+2*k*A->BlockSize+A->BlockSize;
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
      /* tridiagonal matrix con to operator spin dn*/
      ia=2*A->BlockSize*A->BlockSize;
      for(k=2+NBlockMaxVal;k<=NBlockMaxVal+Min(NBlockMaxCon+1,A->con.N);k++)
      {
        ifer=NFermionShift+2*k*A->BlockSize+A->BlockSize;
        for(i=0;i<A->BlockSize;i++)
        {
          jfer=NFermionShift+2*k*A->BlockSize+A->BlockSize;
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
      for(k=2+NBlockMaxVal;k<NBlockMaxVal+Min(NBlockMaxCon+1,A->con.N);k++)
      {
        ifer=NFermionShift+2*(k+1)*A->BlockSize+A->BlockSize;
        for(i=0;i<A->BlockSize;i++)
        {
          jfer=NFermionShift+2*k*A->BlockSize+A->BlockSize;
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
          if(NBlockMaxVal!=0)
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
          if(NBlockMaxCon!=0)
          {
            oadd[1]=CrAnMask+NFermionShift+i+(2+NBlockMaxVal)*A->BlockSize;
            oadd[0]=NFermionShift+j;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[(i+3*A->BlockSize)*A->TMT.Nj+j], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+i+(2+NBlockMaxVal)*A->BlockSize;
            oadd[0]=NFermionShift+j+A->BlockSize;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[(i+3*A->BlockSize)*A->TMT.Nj+j+A->BlockSize], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+i;
            oadd[0]=NFermionShift+j+(2+NBlockMaxVal)*A->BlockSize;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[i*A->TMT.Nj+j+3*A->BlockSize], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+i+A->BlockSize;
            oadd[0]=NFermionShift+j+(2+NBlockMaxVal)*A->BlockSize;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[(i+A->BlockSize)*A->TMT.Nj+j+3*A->BlockSize], &QDetNotOrderedFcFaBcBa);
          }
          
        }
      }
      /* tridiagonal matrix val to operator */
      ia=2*A->BlockSize*A->BlockSize;
      for(k=2;k<=Min(NBlockMaxVal+1,A->val.N);k++)
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
      for(k=2;k<Min(NBlockMaxVal+1,A->val.N);k++)
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
      for(k=2+NBlockMaxVal;k<=NBlockMaxVal+Min(NBlockMaxCon+1,A->con.N);k++)
      {
        ifer=NFermionShift+k*A->BlockSize;
        for(i=0;i<A->BlockSize;i++)
        {
          jfer=NFermionShift+k*A->BlockSize;
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
      for(k=2+NBlockMaxVal;k<NBlockMaxVal+Min(NBlockMaxCon+1,A->con.N);k++)
      {
        ifer=NFermionShift+(k+1)*A->BlockSize;
        for(i=0;i<A->BlockSize;i++)
        {
          jfer=NFermionShift+k*A->BlockSize;
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
      /* transform matrix TMT to Operator spin up*/
      for(i=0;i<A->BlockSize;i++)
      {
        for(j=0;j<A->BlockSize;j++)
        {
          oadd[1]=CrAnMask+NFermionShift+i;
          oadd[0]=NFermionShift+j;
          RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[i*A->TMT.Nj+j], &QDetNotOrderedFcFaBcBa);
          oadd[1]=CrAnMask+NFermionShift+i+2*A->BlockSize;
          oadd[0]=NFermionShift+j;
          RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[(i+A->BlockSize)*A->TMT.Nj+j], &QDetNotOrderedFcFaBcBa);
          oadd[1]=CrAnMask+NFermionShift+i;
          oadd[0]=NFermionShift+j+2*A->BlockSize;
          RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[i*A->TMT.Nj+j+A->BlockSize], &QDetNotOrderedFcFaBcBa);
          oadd[1]=CrAnMask+NFermionShift+i+2*A->BlockSize;
          oadd[0]=NFermionShift+j+2*A->BlockSize;
          RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[(i+A->BlockSize)*A->TMT.Nj+j+A->BlockSize], &QDetNotOrderedFcFaBcBa);
          if(NBlockMaxVal!=0)
          {
            oadd[1]=CrAnMask+NFermionShift+i+4*A->BlockSize;
            oadd[0]=NFermionShift+j;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[(i+2*A->BlockSize)*A->TMT.Nj+j], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+i+4*A->BlockSize;
            oadd[0]=NFermionShift+j+2*A->BlockSize;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[(i+2*A->BlockSize)*A->TMT.Nj+j+A->BlockSize], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+i;
            oadd[0]=NFermionShift+j+4*A->BlockSize;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[i*A->TMT.Nj+j+2*A->BlockSize], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+i+2*A->BlockSize;
            oadd[0]=NFermionShift+j+4*A->BlockSize;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[(i+A->BlockSize)*A->TMT.Nj+j+2*A->BlockSize], &QDetNotOrderedFcFaBcBa);
          }
          if(NBlockMaxCon!=0)
          {
            oadd[1]=CrAnMask+NFermionShift+i+2*(2+NBlockMaxVal)*A->BlockSize;
            oadd[0]=NFermionShift+j;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[(i+3*A->BlockSize)*A->TMT.Nj+j], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+i+2*(2+NBlockMaxVal)*A->BlockSize;
            oadd[0]=NFermionShift+j+2*A->BlockSize;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[(i+3*A->BlockSize)*A->TMT.Nj+j+A->BlockSize], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+i;
            oadd[0]=NFermionShift+j+2*(2+NBlockMaxVal)*A->BlockSize;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[i*A->TMT.Nj+j+3*A->BlockSize], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+i+2*A->BlockSize;
            oadd[0]=NFermionShift+j+2*(2+NBlockMaxVal)*A->BlockSize;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[(i+A->BlockSize)*A->TMT.Nj+j+3*A->BlockSize], &QDetNotOrderedFcFaBcBa);
          }
        }
      }
      /* tridiagonal matrix val to operator spin up*/
      ia=2*A->BlockSize*A->BlockSize;
      for(k=2;k<=Min(NBlockMaxVal+1,A->val.N);k++)
      {
        ifer=NFermionShift+2*k*A->BlockSize;
        for(i=0;i<A->BlockSize;i++)
        {
          jfer=NFermionShift+2*k*A->BlockSize;
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
      for(k=2;k<Min(NBlockMaxVal+1,A->val.N);k++)
      {
        ifer=NFermionShift+2*(k+1)*A->BlockSize;
        for(i=0;i<A->BlockSize;i++)
        {
          jfer=NFermionShift+2*k*A->BlockSize;
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
      /* tridiagonal matrix con to operator spin up*/
      ia=2*A->BlockSize*A->BlockSize;
      for(k=2+NBlockMaxVal;k<=NBlockMaxVal+Min(NBlockMaxCon+1,A->con.N);k++)
      {
        ifer=NFermionShift+2*k*A->BlockSize;
        for(i=0;i<A->BlockSize;i++)
        {
          jfer=NFermionShift+2*k*A->BlockSize;
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
      for(k=2+NBlockMaxVal;k<NBlockMaxVal+Min(NBlockMaxCon+1,A->con.N);k++)
      {
        ifer=NFermionShift+2*(k+1)*A->BlockSize;
        for(i=0;i<A->BlockSize;i++)
        {
          jfer=NFermionShift+2*k*A->BlockSize;
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
      /* transform matrix TMT to Operator spin dn*/
      for(i=0;i<A->BlockSize;i++)
      {
        for(j=0;j<A->BlockSize;j++)
        {
          oadd[1]=CrAnMask+NFermionShift+i+A->BlockSize;
          oadd[0]=NFermionShift+j+A->BlockSize;
          RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[i*A->TMT.Nj+j], &QDetNotOrderedFcFaBcBa);
          oadd[1]=CrAnMask+NFermionShift+i+3*A->BlockSize;
          oadd[0]=NFermionShift+j+A->BlockSize;
          RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[(i+A->BlockSize)*A->TMT.Nj+j], &QDetNotOrderedFcFaBcBa);
          oadd[1]=CrAnMask+NFermionShift+i+A->BlockSize;
          oadd[0]=NFermionShift+j+3*A->BlockSize;
          RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[i*A->TMT.Nj+j+A->BlockSize], &QDetNotOrderedFcFaBcBa);
          oadd[1]=CrAnMask+NFermionShift+i+3*A->BlockSize;
          oadd[0]=NFermionShift+j+3*A->BlockSize;
          RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[(i+A->BlockSize)*A->TMT.Nj+j+A->BlockSize], &QDetNotOrderedFcFaBcBa);
          if(NBlockMaxVal!=0)
          {
            oadd[1]=CrAnMask+NFermionShift+i+5*A->BlockSize;
            oadd[0]=NFermionShift+j+A->BlockSize;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[(i+2*A->BlockSize)*A->TMT.Nj+j], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+i+5*A->BlockSize;
            oadd[0]=NFermionShift+j+3*A->BlockSize;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[(i+2*A->BlockSize)*A->TMT.Nj+j+A->BlockSize], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+i+A->BlockSize;
            oadd[0]=NFermionShift+j+5*A->BlockSize;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[i*A->TMT.Nj+j+2*A->BlockSize], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+i+3*A->BlockSize;
            oadd[0]=NFermionShift+j+5*A->BlockSize;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[(i+A->BlockSize)*A->TMT.Nj+j+2*A->BlockSize], &QDetNotOrderedFcFaBcBa);
          }
          if(NBlockMaxCon!=0)
          {
            oadd[1]=CrAnMask+NFermionShift+i+2*(2+NBlockMaxVal)*A->BlockSize+A->BlockSize;
            oadd[0]=NFermionShift+j+A->BlockSize;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[(i+3*A->BlockSize)*A->TMT.Nj+j], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+i+2*(2+NBlockMaxVal)*A->BlockSize+A->BlockSize;
            oadd[0]=NFermionShift+j+3*A->BlockSize;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[(i+3*A->BlockSize)*A->TMT.Nj+j+A->BlockSize], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+i+A->BlockSize;
            oadd[0]=NFermionShift+j+2*(2+NBlockMaxVal)*A->BlockSize+A->BlockSize;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[i*A->TMT.Nj+j+3*A->BlockSize], &QDetNotOrderedFcFaBcBa);
            oadd[1]=CrAnMask+NFermionShift+i+3*A->BlockSize;
            oadd[0]=NFermionShift+j+2*(2+NBlockMaxVal)*A->BlockSize+A->BlockSize;
            RealOperatorAddLaderNormalOrder(H, oadd, 2, A->TMT.VRe[(i+A->BlockSize)*A->TMT.Nj+j+3*A->BlockSize], &QDetNotOrderedFcFaBcBa);
          }
        }
      }
      /* tridiagonal matrix val to operator spin dn*/
      ia=2*A->BlockSize*A->BlockSize;
      for(k=2;k<=Min(NBlockMaxVal+1,A->val.N);k++)
      {
        ifer=NFermionShift+2*k*A->BlockSize+A->BlockSize;
        for(i=0;i<A->BlockSize;i++)
        {
          jfer=NFermionShift+2*k*A->BlockSize+A->BlockSize;
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
      for(k=2;k<Min(NBlockMaxVal+1,A->val.N);k++)
      {
        ifer=NFermionShift+2*(k+1)*A->BlockSize+A->BlockSize;
        for(i=0;i<A->BlockSize;i++)
        {
          jfer=NFermionShift+2*k*A->BlockSize+A->BlockSize;
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
      /* tridiagonal matrix con to operator spin dn*/
      ia=2*A->BlockSize*A->BlockSize;
      for(k=2+NBlockMaxVal;k<=NBlockMaxVal+Min(NBlockMaxCon+1,A->con.N);k++)
      {
        ifer=NFermionShift+2*k*A->BlockSize+A->BlockSize;
        for(i=0;i<A->BlockSize;i++)
        {
          jfer=NFermionShift+2*k*A->BlockSize+A->BlockSize;
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
      for(k=2+NBlockMaxVal;k<NBlockMaxVal+Min(NBlockMaxCon+1,A->con.N);k++)
      {
        ifer=NFermionShift+2*(k+1)*A->BlockSize+A->BlockSize;
        for(i=0;i<A->BlockSize;i++)
        {
          jfer=NFermionShift+2*k*A->BlockSize+A->BlockSize;
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
  }
  return 0;
}


uint32_t BlockNaturalImpurityOrbitalMatrixToImpurityHamiltonianPerturpation(OperatorType *const H, OperatorType * const VPerturbation, TriDiagonalMatrixType * const gPerturbation, BlockNaturalImpurityOrbitalMatrixType const *A, uint32_t const NBlockMaxVal, uint32_t const NBlockMaxCon, uint32_t const Spin, uint32_t const NFermionShift, uint32_t const NFermion)
{
  //TODO: implement
  printf("To DO\n");
  return 0;
}
