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
  InitOperatorFixedLength(H, 2, 8);
  if(A->QComplex)
  {
  
  }
  else /* H is real */
  {
    if(Spin==0)
    {
      /* copy matrix TMT to Operator */
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
      
      ia=2*A->BlockSize*A->BlockSize;
      for(k=2;k<=Min(NBlockMaxVal,A->val.N);k++)
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
      for(k=2;k<Min(NBlockMaxVal,A->val.N);k++)
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
      
      ia=2*A->BlockSize*A->BlockSize;
      for(k=2+NBlockMaxVal;k<=Min(NBlockMaxCon+NBlockMaxVal,A->con.N+NBlockMaxVal);k++)
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
      for(k=2+NBlockMaxVal;k<Min(NBlockMaxCon+NBlockMaxVal,A->con.N+NBlockMaxVal);k++)
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
    
    }
  }
  return 0;
}