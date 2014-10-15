/*  Starting from a set of Npsi wavefunctions calculate the lowest N states with the same basis
 *
 *  On input:
 *    H operator whose eigenstates are determined
 *    psi a vector of Npsi functions.
 *       Next, Determiantn and KeyList are the same for all functions and point to the same location.
 *    zero smallest positive number considered non-zero for convergence tests
 */
uint32_t BlockLanczosGroundStateConserveBasisRR(const OperatorType H, WaveFunctionType * const psi, uint32_t const Npsi, const double zero)
{
  printf("To Do\n");
  
  WaveFunctionType v[Npsi][3];
  WaveFunctionType * vflat=&(v[0][0]);
  const uint32_t NKrylovMax = 100, NKrylovStart = 10;
  uint32_t converged=0, NKrylov=NKrylovStart;
  uint32_t i,j,k,vim,vi,vip;
  double Norm;
  BlockTriDiagonalMatrixType Htri;
  // init Htri
  Htri.QComplex = H.QComplex;
  Htri.BlockSize=Npsi;
  Htri.N=NKrylovMax;
  InitBlockTridiagonalMatrix(&Htri);
  // copy psi[i] into v[i][0], only copy pointers to allocated memory
  for(i=0;i<Npsi;i++)
  {
    v[i][0]=psi[i];
  }
  // allocate memory for v[i][1] and v[i][2]
  for(i=0;i<Npsi;i++)
  {
    if(OperatorPsiConserveBasis(H, v[i][0], &v[i][1]))
    {
      printf("OperatorPsiConserveBasis failed in BlockLanczosGroundStateConserveBasis\n");
      for(j=0;j<i;j++)
      {
        WaveFunctionFreeValues(&v[j][1]);
        WaveFunctionFreeValues(&v[j][2]);
      }
      return 1;
    }
    if(OperatorPsiConserveBasis(H, v[i][1], &v[i][2]))
    {
      printf("OperatorPsiConserveBasis failed in BlockLanczosGroundStateConserveBasis\n");
      WaveFunctionFreeValues(&v[i][1]);
      for(j=0;j<i;j++)
      {
        WaveFunctionFreeValues(&v[j][1]);
        WaveFunctionFreeValues(&v[j][2]);
      }
      return 1;
    }
  }
  for(i=1;i<3*Npsi;i++)
  {
    WaveFunctionOrthogonalizeConserveBasis(&(vflat[i]), vflat, i);
    WaveFunctionCalculateNorm(vflat[i], &Norm);
    
  }
  
  
  while(!converged)
  {
    // start an other Block Lanczos loop
    for(i=1;i<NKrylov;i++)
    {
      // loop NKrylov times to create a tri-diagonal representation of H
      vip = (i + 1) % 3;
      vi = (i) % 3;
      vim = (i - 1) % 3;
      for(j=0;j<Npsi;j++)
      {
        OperatorPsiConserveBasisPreAllocated(&H, &(v[j][vi]), &(v[j][vip]));
      }
      // a_i = <vi | H | vi > == <vi | vip > (at the moment)
      for(j=0;j<Npsi;j++)
      {
        for(k=0;k<Npsi;k++)
        {
          //WaveFunctionDotProductSameBasis(v[j][vi], v[k][vip], <#double *DotRe#>, <#double *DotIm#>);
        }
      }
      
      
      
      
      
    }
      
    
    
    // Calculate the eigen-states in the Krylov basis
    
    // Given the eigen-states in the Krylov basis reconstruct the ground-state in the original basis
    
    // test for convergence
    // variance of all states (basis conserving) must be smaller than zero
  
  }
  return 0;
}