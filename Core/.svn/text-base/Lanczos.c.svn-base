/*
*   Lanczos Subroutines
*   Functions Specified:
*     LanczosGroundState returns the lowest eigenvector of a given hamiltonian
*     LanczosTriDiagonalize returns a list containing a and b of a tridiagonal matrix given an hamiltonian and initial wave function
*
*/

uint32_t LanczosGroundStateRR(const OperatorType H, WaveFunctionType *const psi, double *epsilon, const double zero, double shiftin, WaveFunctionType *OtherPsi, uint32_t NOtherPsi)
{
  printf("\nStart of LanczosGroundStateRR\n\n");
  fflush(stdout);
  const uint32_t NKrylovMax=100;
  uint32_t NKrylov;
  WaveFunctionType v[3];
  double energy=4*zero,energy1=2*zero,energy2;
  double shift;
  uint32_t cycles,i,vi,vip,vim;
  MatrixType Eigen;
  TriDiagonalMatrixType M;
  double val[NKrylovMax];
  M.N=NKrylovMax;
  InitTriDiagonalMatrix(&M);
  while((energy-energy1)>(zero*(1.0+fabs(energy))))
  {
    cycles=0;
    energy=energy1;
    shift=-shiftin-energy*0.49;
    Repeat:
      if(OperatorPsiRR(H,*psi,&v[0]))                           goto OutOffMemv0;
      RealWaveFunctionMultiply(*psi,shift);
      if(RealWaveFunctionAddUnionBasis(&v[0],*psi))                 goto OutOffMemv0add;
      if(RealWaveFunctionOrthogonalize(&v[0], OtherPsi, NOtherPsi)) goto OutOffMemv0Ortho;
      RealWaveFunctionNormalize(&v[0]);
      printf("BasisSetSize After Expansion: %12i Truncation: ",v[0].N);
      RealWaveFunctionTruncateBasis(&v[0], *epsilon);
      printf("%12i\n",v[0].N);
      fflush(stdout);
      if(OperatorPsiConserveBasisRR(H,v[0],&v[1]))              goto OutOffMemv1;
      RealWaveFunctionMultiply(v[0], shift);
      RealWaveFunctionAddSameBasis(v[1],v[0]);
      RealWaveFunctionNormalize(&v[1]);
      if(OperatorPsiConserveBasisRR(H,v[1],&v[2]))              goto OutOffMemv2;
      RealWaveFunctionMultiply(v[1], shift);
      RealWaveFunctionAddSameBasis(v[2],v[1]);
      RealWaveFunctionNormalize(&v[2]);
      RealWaveFunctionFree(psi);
      if(OperatorPsiConserveBasisRR(H,v[2],psi))                goto OutOffMempsi;
      RealWaveFunctionMultiply(v[2], shift);
      RealWaveFunctionAddSameBasis(*psi,v[2]);
      RealWaveFunctionOrthogonalizeConserveBasis(psi, OtherPsi, NOtherPsi);
      RealWaveFunctionNormalize(psi);
    goto EnoughMem;
      OutOffMempsi:
        RealWaveFunctionFreeValues(&v[2]);
      OutOffMemv2:
        RealWaveFunctionFreeValues(&v[1]);
      OutOffMemv1:
      OutOffMemv0Ortho:
      OutOffMemv0add:
        RealWaveFunctionFreeValues(&v[0]);
        RealWaveFunctionMultiply(*psi,1.0/shift);
      OutOffMemv0:
        printf("Out off memory increase epsilon from %22.15E to %22.15E\n",*epsilon,*epsilon*2);
        *epsilon*=2.0;
        if(psi->NBlok==1) {printf("Outof memory in LanczosGroundStateRR,\npsi->NBLok=1 so no way to recover.\n");AllocFailed();}
        printf("BasisSetSize Before Truncation: %12i ",psi->N);
        RealWaveFunctionTruncateBasis(psi,*epsilon);
        printf("After Truncation: %12i\n",psi->N);
        fflush(stdout);
    goto Repeat;
    EnoughMem:

    energy1=2*zero;
    energy2=4*zero;
    while((cycles<3)||((energy2-energy1)>(zero*(1.0+fabs(energy1)))))
    {
      NKrylov=5+cycles%(NKrylovMax-5);
      M.N=NKrylov;
      
      cycles++;
      if((cycles&15)==15)
      {
        printf("Reached %i-th cycle with energy %22.15E\n",cycles,energy1);
        fflush(stdout);
      }
      energy2=energy1;
      RealWaveFunctionMultiply(v[0],0);
      RealWaveFunctionAddSameBasis(v[0],*psi);
      OperatorPsiConserveBasisPreAllocatedRR(H,v[0],&(v[1]));
      M.a[1]=RealWaveFunctionDotProduct(v[1],v[0]);
      for(i=1;i<M.N;i++)
      {
        vip=(i+1)%3;
        vi =(i  )%3;
        vim=(i-1)%3;
        RealWaveFunctionAddMultipliedConserveBasis(&(v[vi]),v[vim],-M.a[i]);
        M.b[i]=RealWaveFunctionCalculateNorm(v[vi]);
        if(M.b[i]>DBL_MIN*10) RealWaveFunctionMultiply(v[vi],(double)1.0/M.b[i]);
        else RealWaveFunctionMultiply(v[vi],0.0);
        RealWaveFunctionOrthogonalizeConserveBasis(&v[vi], OtherPsi, NOtherPsi);
        OperatorPsiConserveBasisPreAllocatedRR(H,v[vi],&(v[vip]));
        M.a[i+1]=RealWaveFunctionDotProduct(v[vip],v[vi]);
        RealWaveFunctionAddMultipliedConserveBasis(&(v[vip]),v[vim],-M.b[i]);
      }
      TriDiagonalMatrixToEigensystem(M, &Eigen, val);
      RealWaveFunctionMultiply(v[0],0);
      RealWaveFunctionAddSameBasis(v[0],*psi);
      RealWaveFunctionMultiply(*psi,Eigen.VRe[0][0]);
      OperatorPsiConserveBasisPreAllocatedRR(H,v[0],&(v[1]));
      for(i=1;i<M.N;i++)
      {
        vip=(i+1)%3;
        vi =(i  )%3;
        vim=(i-1)%3;
        RealWaveFunctionAddMultipliedSameBasis(v[vi],v[vim],-M.a[i]);
        if(M.b[i]>DBL_MIN*10) RealWaveFunctionMultiply(v[vi],(double)1.0/M.b[i]);
        else RealWaveFunctionMultiply(v[vi],0.0);
        RealWaveFunctionOrthogonalizeConserveBasis(&v[vi], OtherPsi, NOtherPsi);
        OperatorPsiConserveBasisPreAllocatedRR(H,v[vi],&(v[vip]));
        RealWaveFunctionAddMultipliedSameBasis(v[vip],v[vim],-M.b[i]);
        RealWaveFunctionAddMultipliedSameBasis(*psi, v[vi], Eigen.VRe[0][i]);
      }
      RealWaveFunctionOrthogonalizeConserveBasis(psi, OtherPsi, NOtherPsi);
      RealWaveFunctionNormalize(psi);
      PsiOperatorPsiRR(H,*psi,&energy1);
    }
    RealWaveFunctionFreeValues(&v[0]);
    RealWaveFunctionFreeValues(&v[1]);
    RealWaveFunctionFreeValues(&v[2]);
    FreeMatrix(&Eigen);
    printf("Energy: %22.15E  ",energy1);
    printf("NumberOfLaczosCycles: %i\n",cycles);
    fflush(stdout);
  }
  return 0;
}

uint32_t LanczosGroundStateRC(const OperatorType H, WaveFunctionType *const psi, double *epsilon, const double zero, double shiftin, WaveFunctionType *OtherPsi, uint32_t NOtherPsi)
{
  printf("\nStart of LanczosGroundStateRC\n\n");
  fflush(stdout);
  const uint32_t NKrylovMax=100;
  uint32_t NKrylov;
  WaveFunctionType v[3];
  double energy=4*zero,energy1=2*zero,energy2;
  double Imzero;
  double shift;
  uint32_t cycles,i,vi,vip,vim;
  MatrixType Eigen;
  TriDiagonalMatrixType M;
  double val[NKrylovMax];
  M.N=NKrylovMax;
  InitTriDiagonalMatrix(&M);
  while((energy-energy1)>(zero*(1.0+fabs(energy))))
  {
    cycles=0;
    energy=energy1;
    shift=-shiftin-energy*0.49;
    Repeat:
      if(OperatorPsiRC(H,*psi,&v[0]))                           goto OutOffMemv0;
      ComplexWaveFunctionMultiply(*psi,shift);
      if(ComplexWaveFunctionAddUnionBasis(&v[0],*psi))                 goto OutOffMemv0add;
      if(ComplexWaveFunctionOrthogonalize(&v[0], OtherPsi, NOtherPsi)) goto OutOffMemv0Ortho;
      ComplexWaveFunctionNormalize(&v[0]);
      printf("BasisSetSize After Expansion: %12i Truncation: ",v[0].N);
      ComplexWaveFunctionTruncateBasis(&v[0], *epsilon);
      printf("%12i\n",v[0].N);
      fflush(stdout);
      if(OperatorPsiConserveBasisRC(H,v[0],&v[1]))              goto OutOffMemv1;
      ComplexWaveFunctionMultiply(v[0], shift);
      ComplexWaveFunctionAddSameBasis(v[1],v[0]);
      ComplexWaveFunctionNormalize(&v[1]);
      if(OperatorPsiConserveBasisRC(H,v[1],&v[2]))              goto OutOffMemv2;
      ComplexWaveFunctionMultiply(v[1], shift);
      ComplexWaveFunctionAddSameBasis(v[2],v[1]);
      ComplexWaveFunctionNormalize(&v[2]);
      ComplexWaveFunctionFree(psi);
      if(OperatorPsiConserveBasisRC(H,v[2],psi))                goto OutOffMempsi;
      ComplexWaveFunctionMultiply(v[2], shift);
      ComplexWaveFunctionAddSameBasis(*psi,v[2]);
      ComplexWaveFunctionOrthogonalizeConserveBasis(psi, OtherPsi, NOtherPsi);
      ComplexWaveFunctionNormalize(psi);
    goto EnoughMem;
      OutOffMempsi:
        ComplexWaveFunctionFreeValues(&v[2]);
      OutOffMemv2:
        ComplexWaveFunctionFreeValues(&v[1]);
      OutOffMemv1:
      OutOffMemv0Ortho:
      OutOffMemv0add:
        ComplexWaveFunctionFreeValues(&v[0]);
        ComplexWaveFunctionMultiply(*psi,1.0/shift);
      OutOffMemv0:
        printf("Out off memory increase epsilon from %22.15E to %22.15E\n",*epsilon,*epsilon*2);
        *epsilon*=2.0;
        if(psi->NBlok==1) {printf("Outof memory in LanczosGroundStateRC,\npsi->NBLok=1 so no way to recover.\n");AllocFailed();}
        printf("BasisSetSize Before Truncation: %12i ",psi->N);
        ComplexWaveFunctionTruncateBasis(psi,*epsilon);
        printf("After Truncation: %12i\n",psi->N);
        fflush(stdout);
    goto Repeat;
    EnoughMem:

    energy1=2*zero;
    energy2=4*zero;
    while((cycles<3)||((energy2-energy1)>(zero*(1.0+fabs(energy1)))))
    {
      NKrylov=5+cycles%(NKrylovMax-5);
      M.N=NKrylov;
      
      cycles++;
      if((cycles&15)==15)
      {
        printf("Reached %i-th cycle with energy %22.15E\n",cycles,energy1);
        fflush(stdout);
      }
      energy2=energy1;
      ComplexWaveFunctionMultiply(v[0],0);
      ComplexWaveFunctionAddSameBasis(v[0],*psi);
      OperatorPsiConserveBasisPreAllocatedRC(H,v[0],&(v[1]));
      ComplexWaveFunctionDotProduct(v[1],v[0],&(M.a[1]),&Imzero);
      for(i=1;i<M.N;i++)
      {
        vip=(i+1)%3;
        vi =(i  )%3;
        vim=(i-1)%3;
        ComplexWaveFunctionAddMultipliedConserveBasis(&(v[vi]),v[vim],-M.a[i]);
        M.b[i]=ComplexWaveFunctionCalculateNorm(v[vi]);
        if(M.b[i]>DBL_MIN*10) ComplexWaveFunctionMultiply(v[vi],(double)1.0/M.b[i]);
        else ComplexWaveFunctionMultiply(v[vi],0.0);
        ComplexWaveFunctionOrthogonalizeConserveBasis(&v[vi], OtherPsi, NOtherPsi);
        OperatorPsiConserveBasisPreAllocatedRC(H,v[vi],&(v[vip]));
        ComplexWaveFunctionDotProduct(v[vip],v[vi],&(M.a[i+1]),&Imzero);
        ComplexWaveFunctionAddMultipliedConserveBasis(&(v[vip]),v[vim],-M.b[i]);
      }
      TriDiagonalMatrixToEigensystem(M, &Eigen, val);
      ComplexWaveFunctionMultiply(v[0],0);
      ComplexWaveFunctionAddSameBasis(v[0],*psi);
      ComplexWaveFunctionMultiply(*psi,Eigen.VRe[0][0]);
      OperatorPsiConserveBasisPreAllocatedRC(H,v[0],&(v[1]));
      for(i=1;i<M.N;i++)
      {
        vip=(i+1)%3;
        vi =(i  )%3;
        vim=(i-1)%3;
        ComplexWaveFunctionAddMultipliedSameBasis(v[vi],v[vim],-M.a[i]);
        if(M.b[i]>DBL_MIN*10) RealWaveFunctionMultiply(v[vi],(double)1.0/M.b[i]);
        else ComplexWaveFunctionMultiply(v[vi],0.0);
        ComplexWaveFunctionOrthogonalizeConserveBasis(&v[vi], OtherPsi, NOtherPsi);
        OperatorPsiConserveBasisPreAllocated(H,v[vi],&(v[vip]));
        ComplexWaveFunctionAddMultipliedSameBasis(v[vip],v[vim],-M.b[i]);
        ComplexWaveFunctionAddMultipliedSameBasis(*psi, v[vi], Eigen.VRe[0][i]);
      }
      ComplexWaveFunctionOrthogonalizeConserveBasis(psi, OtherPsi, NOtherPsi);
      ComplexWaveFunctionNormalize(psi);
      PsiOperatorPsiRHC(H,*psi,&energy1);
    }
    ComplexWaveFunctionFreeValues(&v[0]);
    ComplexWaveFunctionFreeValues(&v[1]);
    ComplexWaveFunctionFreeValues(&v[2]);
    FreeMatrix(&Eigen);
    printf("Energy: %22.15E  ",energy1);
    printf("NumberOfLaczosCycles: %i\n",cycles);
    fflush(stdout);
  }
  return 0;
}

uint32_t LanczosGroundStateMC(const OperatorType H, WaveFunctionType *const psi, double *epsilon, const double zero, double shiftin, WaveFunctionType *OtherPsi, uint32_t NOtherPsi)
{
  printf("\nStart of LanczosGroundStateMC\n\n");
  fflush(stdout);
  const uint32_t NKrylovMax=100;
  uint32_t NKrylov;
  WaveFunctionType v[3];
  double energy=4*zero,energy1=2*zero,energy2;
  double Imzero;
  double shift;
  uint32_t cycles,i,vi,vip,vim;
  MatrixType Eigen;
  TriDiagonalMatrixType M;
  double val[NKrylovMax];
  M.N=NKrylovMax;
  InitTriDiagonalMatrix(&M);
  while((energy-energy1)>(zero*(1.0+fabs(energy))))
  {
    cycles=0;
    energy=energy1;
    shift=-shiftin-energy*0.49;
    Repeat:
      if(OperatorPsiMC(H,*psi,&v[0]))                           goto OutOffMemv0;
      ComplexWaveFunctionMultiply(*psi,shift);
      if(ComplexWaveFunctionAddUnionBasis(&v[0],*psi))                 goto OutOffMemv0add;
      if(ComplexWaveFunctionOrthogonalize(&v[0], OtherPsi, NOtherPsi)) goto OutOffMemv0Ortho;
      ComplexWaveFunctionNormalize(&v[0]);
      printf("BasisSetSize After Expansion: %12i Truncation: ",v[0].N);
      ComplexWaveFunctionTruncateBasis(&v[0], *epsilon);
      printf("%12i\n",v[0].N);
      fflush(stdout);
      if(OperatorPsiConserveBasisMC(H,v[0],&v[1]))              goto OutOffMemv1;
      ComplexWaveFunctionMultiply(v[0], shift);
      ComplexWaveFunctionAddSameBasis(v[1],v[0]);
      ComplexWaveFunctionNormalize(&v[1]);
      if(OperatorPsiConserveBasisMC(H,v[1],&v[2]))              goto OutOffMemv2;
      ComplexWaveFunctionMultiply(v[1], shift);
      ComplexWaveFunctionAddSameBasis(v[2],v[1]);
      ComplexWaveFunctionNormalize(&v[2]);
      ComplexWaveFunctionFree(psi);
      if(OperatorPsiConserveBasisMC(H,v[2],psi))                goto OutOffMempsi;
      ComplexWaveFunctionMultiply(v[2], shift);
      ComplexWaveFunctionAddSameBasis(*psi,v[2]);
      ComplexWaveFunctionOrthogonalizeConserveBasis(psi, OtherPsi, NOtherPsi);
      ComplexWaveFunctionNormalize(psi);
    goto EnoughMem;
      OutOffMempsi:
        ComplexWaveFunctionFreeValues(&v[2]);
      OutOffMemv2:
        ComplexWaveFunctionFreeValues(&v[1]);
      OutOffMemv1:
      OutOffMemv0Ortho:
      OutOffMemv0add:
        ComplexWaveFunctionFreeValues(&v[0]);
        ComplexWaveFunctionMultiply(*psi,1.0/shift);
      OutOffMemv0:
        printf("Out off memory increase epsilon from %22.15E to %22.15E\n",*epsilon,*epsilon*2);
        *epsilon*=2.0;
        if(psi->NBlok==1) {printf("Outof memory in LanczosGroundStateMC,\npsi->NBLok=1 so no way to recover.\n");AllocFailed();}
        printf("BasisSetSize Before Truncation: %12i ",psi->N);
        ComplexWaveFunctionTruncateBasis(psi,*epsilon);
        printf("After Truncation: %12i\n",psi->N);
        fflush(stdout);
    goto Repeat;
    EnoughMem:

    energy1=2*zero;
    energy2=4*zero;
    while((cycles<3)||((energy2-energy1)>(zero*(1.0+fabs(energy1)))))
    {
      NKrylov=5+cycles%(NKrylovMax-5);
      M.N=NKrylov;
      
      cycles++;
      if((cycles&15)==15)
      {
        printf("Reached %i-th cycle with energy %22.15E\n",cycles,energy1);
        fflush(stdout);
      }
      energy2=energy1;
      ComplexWaveFunctionMultiply(v[0],0);
      ComplexWaveFunctionAddSameBasis(v[0],*psi);
      OperatorPsiConserveBasisPreAllocatedMC(H,v[0],&(v[1]));
      ComplexWaveFunctionDotProduct(v[1],v[0],&(M.a[1]),&Imzero);
      for(i=1;i<M.N;i++)
      {
        vip=(i+1)%3;
        vi =(i  )%3;
        vim=(i-1)%3;
        ComplexWaveFunctionAddMultipliedConserveBasis(&(v[vi]),v[vim],-M.a[i]);
        M.b[i]=ComplexWaveFunctionCalculateNorm(v[vi]);
        if(M.b[i]>DBL_MIN*10) ComplexWaveFunctionMultiply(v[vi],(double)1.0/M.b[i]);
        else ComplexWaveFunctionMultiply(v[vi],0.0);
        ComplexWaveFunctionOrthogonalizeConserveBasis(&v[vi], OtherPsi, NOtherPsi);
        OperatorPsiConserveBasisPreAllocatedMC(H,v[vi],&(v[vip]));
        ComplexWaveFunctionDotProduct(v[vip],v[vi],&(M.a[i+1]),&Imzero);
        ComplexWaveFunctionAddMultipliedConserveBasis(&(v[vip]),v[vim],-M.b[i]);
      }
      TriDiagonalMatrixToEigensystem(M, &Eigen, val);
      ComplexWaveFunctionMultiply(v[0],0);
      ComplexWaveFunctionAddSameBasis(v[0],*psi);
      ComplexWaveFunctionMultiply(*psi,Eigen.VRe[0][0]);
      OperatorPsiConserveBasisPreAllocatedMC(H,v[0],&(v[1]));
      for(i=1;i<M.N;i++)
      {
        vip=(i+1)%3;
        vi =(i  )%3;
        vim=(i-1)%3;
        ComplexWaveFunctionAddMultipliedSameBasis(v[vi],v[vim],-M.a[i]);
        if(M.b[i]>DBL_MIN*10) ComplexWaveFunctionMultiply(v[vi],(double)1.0/M.b[i]);
        else ComplexWaveFunctionMultiply(v[vi],0.0);
        ComplexWaveFunctionOrthogonalizeConserveBasis(&v[vi], OtherPsi, NOtherPsi);
        OperatorPsiConserveBasisPreAllocatedMC(H,v[vi],&(v[vip]));
        ComplexWaveFunctionAddMultipliedSameBasis(v[vip],v[vim],-M.b[i]);
        ComplexWaveFunctionAddMultipliedSameBasis(*psi, v[vi], Eigen.VRe[0][i]);
      }
      ComplexWaveFunctionOrthogonalizeConserveBasis(psi, OtherPsi, NOtherPsi);
      ComplexWaveFunctionNormalize(psi);
      PsiOperatorPsiMHC(H,*psi,&energy1);
    }
    ComplexWaveFunctionFreeValues(&v[0]);
    ComplexWaveFunctionFreeValues(&v[1]);
    ComplexWaveFunctionFreeValues(&v[2]);
    FreeMatrix(&Eigen);
    printf("Energy: %22.15E  ",energy1);
    printf("NumberOfLaczosCycles: %i\n",cycles);
    fflush(stdout);
  }
  return 0;
}

uint32_t LanczosGroundStateCC(const OperatorType H, WaveFunctionType *const psi, double *epsilon, const double zero, double shiftin, WaveFunctionType *OtherPsi, uint32_t NOtherPsi)
{
  printf("\nStart of LanczosGroundStateCC\n\n");
  fflush(stdout);
  const uint32_t NKrylovMax=100;
  uint32_t NKrylov;
  WaveFunctionType v[3];
  double energy=4*zero,energy1=2*zero,energy2;
  double Imzero;
  double shift;
  uint32_t cycles,i,vi,vip,vim;
  MatrixType Eigen;
  TriDiagonalMatrixType M;
  double val[NKrylovMax];
  M.N=NKrylovMax;
  InitTriDiagonalMatrix(&M);
  while((energy-energy1)>(zero*(1.0+fabs(energy))))
  {
    cycles=0;
    energy=energy1;
    shift=-shiftin-energy*0.49;
    Repeat:
      if(OperatorPsiCC(H,*psi,&v[0]))                           goto OutOffMemv0;
      ComplexWaveFunctionMultiply(*psi,shift);
      if(ComplexWaveFunctionAddUnionBasis(&v[0],*psi))                 goto OutOffMemv0add;
      if(ComplexWaveFunctionOrthogonalize(&v[0], OtherPsi, NOtherPsi)) goto OutOffMemv0Ortho;
      ComplexWaveFunctionNormalize(&v[0]);
      printf("BasisSetSize After Expansion: %12i Truncation: ",v[0].N);
      ComplexWaveFunctionTruncateBasis(&v[0], *epsilon);
      printf("%12i\n",v[0].N);
      fflush(stdout);
      if(OperatorPsiConserveBasisCC(H,v[0],&v[1]))              goto OutOffMemv1;
      ComplexWaveFunctionMultiply(v[0], shift);
      ComplexWaveFunctionAddSameBasis(v[1],v[0]);
      ComplexWaveFunctionNormalize(&v[1]);
      if(OperatorPsiConserveBasisCC(H,v[1],&v[2]))              goto OutOffMemv2;
      ComplexWaveFunctionMultiply(v[1], shift);
      ComplexWaveFunctionAddSameBasis(v[2],v[1]);
      ComplexWaveFunctionNormalize(&v[2]);
      ComplexWaveFunctionFree(psi);
      if(OperatorPsiConserveBasisCC(H,v[2],psi))                goto OutOffMempsi;
      ComplexWaveFunctionMultiply(v[2], shift);
      ComplexWaveFunctionAddSameBasis(*psi,v[2]);
      ComplexWaveFunctionOrthogonalizeConserveBasis(psi, OtherPsi, NOtherPsi);
      ComplexWaveFunctionNormalize(psi);
    goto EnoughMem;
      OutOffMempsi:
        ComplexWaveFunctionFreeValues(&v[2]);
      OutOffMemv2:
        ComplexWaveFunctionFreeValues(&v[1]);
      OutOffMemv1:
      OutOffMemv0Ortho:
      OutOffMemv0add:
        ComplexWaveFunctionFreeValues(&v[0]);
        ComplexWaveFunctionMultiply(*psi,1.0/shift);
      OutOffMemv0:
        printf("Out off memory increase epsilon from %22.15E to %22.15E\n",*epsilon,*epsilon*2);
        *epsilon*=2.0;
        if(psi->NBlok==1) {printf("Outof memory in LanczosGroundStateCC,\npsi->NBLok=1 so no way to recover.\n");AllocFailed();}
        printf("BasisSetSize Before Truncation: %12i ",psi->N);
        ComplexWaveFunctionTruncateBasis(psi,*epsilon);
        printf("After Truncation: %12i\n",psi->N);
        fflush(stdout);
    goto Repeat;
    EnoughMem:

    energy1=2*zero;
    energy2=4*zero;
    while((cycles<3)||((energy2-energy1)>(zero*(1.0+fabs(energy1)))))
    {
      NKrylov=5+cycles%(NKrylovMax-5);
      M.N=NKrylov;
      
      cycles++;
      if((cycles&15)==15)
      {
        printf("Reached %i-th cycle with energy %22.15E\n",cycles,energy1);
        fflush(stdout);
      }
      energy2=energy1;
      ComplexWaveFunctionMultiply(v[0],0);
      ComplexWaveFunctionAddSameBasis(v[0],*psi);
      OperatorPsiConserveBasisPreAllocatedCC(H,v[0],&(v[1]));
      ComplexWaveFunctionDotProduct(v[1],v[0],&(M.a[0]),&Imzero);
      for(i=1;i<M.N;i++)
      {
        vip=(i+1)%3;
        vi =(i  )%3;
        vim=(i-1)%3;
        ComplexWaveFunctionAddMultipliedConserveBasis(&(v[vi]),v[vim],-M.a[i]);
        M.b[i]=ComplexWaveFunctionCalculateNorm(v[vi]);
        if(M.b[i]>DBL_MIN*10) ComplexWaveFunctionMultiply(v[vi],(double)1.0/M.b[i]);
        else ComplexWaveFunctionMultiply(v[vi],0.0);
        ComplexWaveFunctionOrthogonalizeConserveBasis(&v[vi], OtherPsi, NOtherPsi);
        OperatorPsiConserveBasisPreAllocatedCC(H,v[vi],&(v[vip]));
        ComplexWaveFunctionDotProduct(v[vip],v[vi],&(M.a[i+1]),&Imzero);
        ComplexWaveFunctionAddMultipliedConserveBasis(&(v[vip]),v[vim],-M.b[i]);
      }
      TriDiagonalMatrixToEigensystem(M, &Eigen, val);
      ComplexWaveFunctionMultiply(v[0],0);
      ComplexWaveFunctionAddSameBasis(v[0],*psi);
      ComplexWaveFunctionMultiply(*psi,Eigen.VRe[0][0]);
      OperatorPsiConserveBasisPreAllocatedCC(H,v[0],&(v[1]));
      for(i=1;i<M.N;i++)
      {
        vip=(i+1)%3;
        vi =(i  )%3;
        vim=(i-1)%3;
        ComplexWaveFunctionAddMultipliedSameBasis(v[vi],v[vim],-M.a[i]);
        if(M.b[i]>DBL_MIN*10) ComplexWaveFunctionMultiply(v[vi],(double)1.0/M.b[i]);
        else ComplexWaveFunctionMultiply(v[vi],0.0);
        ComplexWaveFunctionOrthogonalizeConserveBasis(&v[vi], OtherPsi, NOtherPsi);
        OperatorPsiConserveBasisPreAllocatedCC(H,v[vi],&(v[vip]));
        ComplexWaveFunctionAddMultipliedSameBasis(v[vip],v[vim],-M.b[i]);
        ComplexWaveFunctionAddMultipliedSameBasis(*psi, v[vi], Eigen.VRe[0][i]);
      }
      ComplexWaveFunctionOrthogonalizeConserveBasis(psi, OtherPsi, NOtherPsi);
      ComplexWaveFunctionNormalize(psi);
      PsiOperatorPsiCHC(H,*psi,&energy1);
    }
    ComplexWaveFunctionFreeValues(&v[0]);
    ComplexWaveFunctionFreeValues(&v[1]);
    ComplexWaveFunctionFreeValues(&v[2]);
    FreeMatrix(&Eigen);
    printf("Energy: %22.15E  ",energy1);
    printf("NumberOfLaczosCycles: %i\n",cycles);
    fflush(stdout);
  }
  return 0;
}

uint32_t LanczosGroundStateCR(const OperatorType H, WaveFunctionType *const psi, double *epsilon, const double zero, double shiftin, WaveFunctionType *OtherPsi, uint32_t NOtherPsi)
{
  if(MakeComplexWaveFunctionFromReal(psi,*psi)){printf("Lanczos called with a real wave-function and a complex Hamiltonian.\nTurned the wave-function complex in step 1.\nBut failed.\n");fflush(stdout);return 1;}  
  return LanczosGroundStateCC(H, psi, epsilon, zero, shiftin, OtherPsi, NOtherPsi);
}

uint32_t LanczosGroundStateMR(const OperatorType H, WaveFunctionType *const psi, double *epsilon, const double zero, double shiftin, WaveFunctionType *OtherPsi, uint32_t NOtherPsi)
{
  if(MakeComplexWaveFunctionFromReal(psi,*psi)){printf("Lanczos called with a real wave-function and a mixed complex Hamiltonian.\nTurned the wave-function complex in step 1.\nBut failed.\n");fflush(stdout);return 1;}  
  return LanczosGroundStateCC(H, psi, epsilon, zero, shiftin, OtherPsi, NOtherPsi);
}

uint32_t LanczosTriDiagonalizeRR(const OperatorType H, WaveFunctionType psi, double epsilon, TriDiagonalMatrixType *const M, double Gamma)
{
  printf("\nStart of LanczosTriDiagonalizeRR\n\n");
  fflush(stdout);
  WaveFunctionType v[3];
  uint32_t i,vi,vip,vim;
  double temp,errorestimate=1.0;
  double QRe[3],QIm[3],PRe[3],PIm[3],ValEst[3];
  QRe[0]=1;
  QRe[1]=0;
  QRe[2]=0;
  QIm[0]=0;
  QIm[1]=0;
  QIm[2]=0;
  PRe[0]=0;
  PRe[1]=0;
  PRe[2]=1;
  PIm[0]=0;
  PIm[1]=0;
  PIm[2]=0;
  ValEst[0]=0.0;
  ValEst[1]=0.0;
  ValEst[2]=0.0;
  /* copy the name of psi to v[0] in order to do the loop more cleanely  */
  v[0].N=psi.N;
  v[0].NBlok=psi.NBlok;
  v[0].NFermion=psi.NFermion;
  v[0].NBoson=psi.NBoson;
  v[0].NWord=psi.NWord;
  v[0].NBitsKey=psi.NBitsKey;
  v[0].ValueRe=psi.ValueRe;
  v[0].Determinant=psi.Determinant;
  v[0].Next=psi.Next;
  v[0].KeyList=psi.KeyList;
  v[0].QComplex=psi.QComplex;
  /* finished copying */
  M->a[0]=0;
  M->b[0]=RealWaveFunctionCalculateNorm(v[0]);
  if(M->b[0]==0)
  {
    RealWaveFunctionMultiply(v[0],(double)0.0);
  }
  else
  {
    RealWaveFunctionMultiply(v[0],(double)1.0/M->b[0]);
  }
  printf("BasisSetSize: %12i ",v[0].N);
  RealWaveFunctionTruncateBasis(&(v[0]),epsilon);
  printf("%12i\n",v[0].N);
  while(OperatorPsiRR(H,v[0],&(v[1])))
  { /* Allocation of Memory failed */
    printf("Allocation of new memory failed in LanczosTriDiagonalize for the function OperatorPsi\n");
    printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",epsilon,epsilon*2.0);
    fflush(stdout);
    epsilon*=2.0;
    if(v[0].NBlok==1)AllocFailed();
    RealWaveFunctionTruncateBasis(&(v[0]),epsilon);
  }
  M->a[1]=RealWaveFunctionDotProduct(v[1],v[0]);
  printf("BasisSetSize: %12i ",v[1].N);
  for(i=1;i<M->N;i++)
  {
    vip=(i+1)%3;
    vi =(i  )%3;
    vim=(i-1)%3;
    while(RealWaveFunctionAddMultipliedUnionBasis(&(v[vi]),v[vim],-M->a[i]))
    { /*Allocation of Memory failed */
      printf("Allocation of new memory failed in LanczosTriDiagonalize for the function RealWaveFunctionAddMultipliedUnionBasis(&(v[vi]),v[vim],-M->a[i])\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",epsilon,epsilon*2.0);
      fflush(stdout);
      epsilon*=2.0;
      temp=RealWaveFunctionCalculateNorm(v[vi]);
      if(v[vi].NBlok==1)AllocFailed();
      RealWaveFunctionTruncateBasis(&(v[vi]),epsilon*temp/errorestimate);
    }
    M->b[i]=RealWaveFunctionCalculateNorm(v[vi]);
    if(M->b[i]==0)
    {
      RealWaveFunctionMultiply(v[vi],(double)0.0);
    }
    else
    {
      RealWaveFunctionMultiply(v[vi],(double)1.0/M->b[i]);
    }
    RealWaveFunctionTruncateBasis(&(v[vi]),epsilon/errorestimate);
    printf("%12i ",v[vi].N);
    printf("| M.a[%3i] M.b[%3i]: %12.5E %12.5E | Err: %12.5E\n",i,i-1,M->a[i],M->b[i-1],errorestimate);
    /*
    QRe[vi]=M->b[i-1]*QRe[vip]-0.5*Gamma*QIm[vim];
    QIm[vi]=M->b[i-1]*QIm[vip]+0.5*Gamma*QRe[vim];
    PRe[vi]=M->b[i-1]*PRe[vip]-0.5*Gamma*PIm[vim];
    PIm[vi]=M->b[i-1]*PIm[vip]+0.5*Gamma*PRe[vim];
    ValEst[vi]=fabs((PIm[vi]*QRe[vi]-PRe[vi]*QIm[vi])/(QIm[vi]*QIm[vi]+QRe[vi]*QRe[vi]));
    if(i>5)errorestimate=fabs(2*(ValEst[vi]-ValEst[vim])/(Gamma));
    */
    errorestimate*=M->b[i]*M->b[i]/(M->b[i]*M->b[i]+4*Gamma*Gamma);
    while(OperatorPsiRR(H,v[vi],&(v[vip])))
    { /*Allocation of Memory failed */
      printf("Allocation of new memory failed in LanczosTriDiagonalize for the function OperatorPsi\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",epsilon,epsilon*2.0);
      fflush(stdout);
      epsilon*=2.0;
      if(v[vi].NBlok==1)AllocFailed();
      RealWaveFunctionTruncateBasis(&(v[vi]),epsilon/errorestimate);
    }
    printf("BasisSetSize: %12i ",v[vip].N);
    M->a[i+1]=RealWaveFunctionDotProduct(v[vip],v[vi]);
    while(RealWaveFunctionAddMultipliedUnionBasis(&(v[vip]),v[vim],-M->b[i]))
    { /*Allocation of Memory failed */
      printf("Allocation of new memory failed in LanczosTriDiagonalize for the function RealWaveFunctionAddMultipliedUnionBasis(&(v[vip]),v[vim],-M->b[i])\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",epsilon,epsilon*2.0);
      fflush(stdout);
      epsilon*=2.0;
      temp=RealWaveFunctionCalculateNorm(v[vip]);
      if(v[vip].NBlok==1)AllocFailed();
      RealWaveFunctionTruncateBasis(&(v[vip]),epsilon*temp/errorestimate);
    }
    RealWaveFunctionFree(&v[vim]);
  }
  printf("             | M.a[%3i] M.b[%3i]: %12.5E %12.5E | Err: %12.5E\n\n",i,i-1,M->a[i],M->b[i-1],errorestimate);
  if(M->N>1)
  {
    RealWaveFunctionFree(&v[vi]);
    RealWaveFunctionFree(&v[vip]);
  }
  fflush(stdout);
  return 0;
}

uint32_t LanczosTriDiagonalizeCR(const OperatorType H, WaveFunctionType psi, double epsilon, TriDiagonalMatrixType *const M, double Gamma)
{
  uint32_t retval;
  WaveFunctionType Cpsi;
  if(MakeComplexWaveFunctionFromReal(&Cpsi,psi)){printf("LanczosTriDiagonalizeCR called with a real wave-function and a complex Hamiltonian.\nTurned the wave-function complex in step 1.\nBut failed.\n");fflush(stdout);return 1;}  
  retval=LanczosTriDiagonalizeCC(H, Cpsi, epsilon, M, Gamma);
  MakeRealWaveFunctionFromComplex(Cpsi,&Cpsi);
  return retval;
}

uint32_t LanczosTriDiagonalizeMR(const OperatorType H, WaveFunctionType psi, double epsilon, TriDiagonalMatrixType *const M, double Gamma)
{
  uint32_t retval;
  WaveFunctionType Cpsi;
  if(MakeComplexWaveFunctionFromReal(&Cpsi,psi)){printf("LanczosTriDiagonalizeMR called with a real wave-function and a mixed complex Hamiltonian.\nTurned the wave-function complex in step 1.\nBut failed.\n");fflush(stdout);return 1;}  
  retval=LanczosTriDiagonalizeMC(H, Cpsi, epsilon, M, Gamma);
  MakeRealWaveFunctionFromComplex(Cpsi,&Cpsi);
  return retval;
}

uint32_t LanczosTriDiagonalizeRC(const OperatorType H, WaveFunctionType psi, double epsilon, TriDiagonalMatrixType *const M, double Gamma)
{
  printf("\nStart of LanczosTriDiagonalizeRC\n\n");
  fflush(stdout);
  WaveFunctionType v[3];
  uint32_t i,vi,vip,vim;
  double temp,errorestimate=1.0,junk;
  /* copy the name of psi to v[0] in order to do the loop more cleanely  */
  v[0].N=psi.N;
  v[0].NBlok=psi.NBlok;
  v[0].NFermion=psi.NFermion;
  v[0].NBoson=psi.NBoson;
  v[0].NWord=psi.NWord;
  v[0].NBitsKey=psi.NBitsKey;
  v[0].ValueRe=psi.ValueRe;
  v[0].ValueIm=psi.ValueIm;
  v[0].Determinant=psi.Determinant;
  v[0].Next=psi.Next;
  v[0].KeyList=psi.KeyList;
  v[0].QComplex=psi.QComplex;
  /* finished copying */
  M->a[0]=0;
  M->b[0]=ComplexWaveFunctionCalculateNorm(v[0]);
  ComplexWaveFunctionMultiply(v[0],(double)1.0/M->b[0]);
  printf("BasisSetSize: %12i ",v[0].N);
  ComplexWaveFunctionTruncateBasis(&(v[0]),epsilon);
  printf("%12i\n",v[0].N);
  while(OperatorPsiRC(H,v[0],&(v[1])))
  { /* Allocation of Memory failed */
    printf("Allocation of new memory failed in LanczosTriDiagonalize for the function OperatorPsi\n");
    printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",epsilon,epsilon*2.0);
    fflush(stdout);
    epsilon*=2.0;
    if(v[0].NBlok==1)AllocFailed();
    ComplexWaveFunctionTruncateBasis(&(v[0]),epsilon);
  }
  ComplexWaveFunctionDotProduct(v[1],v[0],&M->a[1],&junk);
  printf("BasisSetSize: %12i ",v[1].N);
  for(i=1;i<M->N;i++)
  {
    vip=(i+1)%3;
    vi =(i  )%3;
    vim=(i-1)%3;
    while(ComplexWaveFunctionAddMultipliedUnionBasis(&(v[vi]),v[vim],-M->a[i]))
    { /*Allocation of Memory failed */
      printf("Allocation of new memory failed in LanczosTriDiagonalize for the function RealWaveFunctionAddMultipliedUnionBasis(&(v[vi]),v[vim],-M->a[i])\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",epsilon,epsilon*2.0);
      fflush(stdout);
      epsilon*=2.0;
      temp=ComplexWaveFunctionCalculateNorm(v[vi]);
      if(v[vi].NBlok==1)AllocFailed();
      ComplexWaveFunctionTruncateBasis(&(v[vi]),epsilon*temp/errorestimate);
    }
    M->b[i]=ComplexWaveFunctionCalculateNorm(v[vi]);
    ComplexWaveFunctionMultiply(v[vi],(double)1.0/M->b[i]);
    ComplexWaveFunctionTruncateBasis(&(v[vi]),epsilon/errorestimate);
    printf("%12i ",v[vi].N);
    printf("| M.a[%3i] M.b[%3i]: %12.5E %12.5E | Err: %12.5E\n",i,i-1,M->a[i],M->b[i-1],errorestimate);
    errorestimate*=M->b[i]*M->b[i]/(M->b[i]*M->b[i]+4*Gamma*Gamma);
    while(OperatorPsiRC(H,v[vi],&(v[vip])))
    { /*Allocation of Memory failed */
      printf("Allocation of new memory failed in LanczosTriDiagonalize for the function OperatorPsi\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",epsilon,epsilon*2.0);
      fflush(stdout);
      epsilon*=2.0;
      if(v[vi].NBlok==1)AllocFailed();
      ComplexWaveFunctionTruncateBasis(&(v[vi]),epsilon/errorestimate);
    }
    printf("BasisSetSize: %12i ",v[vip].N);
    ComplexWaveFunctionDotProduct(v[vip],v[vi],&M->a[i+1],&junk);
    while(ComplexWaveFunctionAddMultipliedUnionBasis(&(v[vip]),v[vim],-M->b[i]))
    { /*Allocation of Memory failed */
      printf("Allocation of new memory failed in LanczosTriDiagonalize for the function RealWaveFunctionAddMultipliedUnionBasis(&(v[vip]),v[vim],-M->b[i])\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",epsilon,epsilon*2.0);
      fflush(stdout);
      epsilon*=2.0;
      temp=ComplexWaveFunctionCalculateNorm(v[vip]);
      if(v[vip].NBlok==1)AllocFailed();
      ComplexWaveFunctionTruncateBasis(&(v[vip]),epsilon*temp/errorestimate);
    }
    ComplexWaveFunctionFree(&v[vim]);
  }
  printf("             | M.a[%3i] M.b[%3i]: %12.5E %12.5E | Err: %12.5E\n\n",i,i-1,M->a[i],M->b[i-1],errorestimate);
  if(M->N>1)
  {
    ComplexWaveFunctionFree(&v[vi]);
    ComplexWaveFunctionFree(&v[vip]);
  }
  fflush(stdout);
  return 0;
}

uint32_t LanczosTriDiagonalizeCC(const OperatorType H, WaveFunctionType psi, double epsilon, TriDiagonalMatrixType *const M, double Gamma)
{
  printf("\nStart of LanczosTriDiagonalizeCC\n\n");
  fflush(stdout);
  WaveFunctionType v[3];
  uint32_t i,vi,vip,vim;
  double temp,errorestimate=1.0,junk;
  /* copy the name of psi to v[0] in order to do the loop more cleanely  */
  v[0].N=psi.N;
  v[0].NBlok=psi.NBlok;
  v[0].NFermion=psi.NFermion;
  v[0].NBoson=psi.NBoson;
  v[0].NWord=psi.NWord;
  v[0].NBitsKey=psi.NBitsKey;
  v[0].ValueRe=psi.ValueRe;
  v[0].ValueIm=psi.ValueIm;
  v[0].Determinant=psi.Determinant;
  v[0].Next=psi.Next;
  v[0].KeyList=psi.KeyList;
  v[0].QComplex=psi.QComplex;
  /* finished copying */
  M->a[0]=0;
  M->b[0]=ComplexWaveFunctionCalculateNorm(v[0]);
  ComplexWaveFunctionMultiply(v[0],(double)1.0/M->b[0]);
  printf("BasisSetSize: %12i ",v[0].N);
  ComplexWaveFunctionTruncateBasis(&(v[0]),epsilon);
  printf("%12i\n",v[0].N);
  while(OperatorPsiCC(H,v[0],&(v[1])))
  { /* Allocation of Memory failed */
    printf("Allocation of new memory failed in LanczosTriDiagonalize for the function OperatorPsi\n");
    printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",epsilon,epsilon*2.0);
    fflush(stdout);
    epsilon*=2.0;
    if(v[0].NBlok==1)AllocFailed();
    ComplexWaveFunctionTruncateBasis(&(v[0]),epsilon);
  }
  ComplexWaveFunctionDotProduct(v[1],v[0],&M->a[1],&junk);
  printf("BasisSetSize: %12i ",v[1].N);
  for(i=1;i<M->N;i++)
  {
    vip=(i+1)%3;
    vi =(i  )%3;
    vim=(i-1)%3;
    while(ComplexWaveFunctionAddMultipliedUnionBasis(&(v[vi]),v[vim],-M->a[i]))
    { /*Allocation of Memory failed */
      printf("Allocation of new memory failed in LanczosTriDiagonalize for the function RealWaveFunctionAddMultipliedUnionBasis(&(v[vi]),v[vim],-M->a[i])\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",epsilon,epsilon*2.0);
      fflush(stdout);
      epsilon*=2.0;
      temp=ComplexWaveFunctionCalculateNorm(v[vi]);
      if(v[vi].NBlok==1)AllocFailed();
      ComplexWaveFunctionTruncateBasis(&(v[vi]),epsilon*temp/errorestimate);
    }
    M->b[i]=ComplexWaveFunctionCalculateNorm(v[vi]);
    ComplexWaveFunctionMultiply(v[vi],(double)1.0/M->b[i]);
    ComplexWaveFunctionTruncateBasis(&(v[vi]),epsilon/errorestimate);
    printf("%12i ",v[vi].N);
    printf("| M.a[%3i] M.b[%3i]: %12.5E %12.5E | Err: %12.5E\n",i,i-1,M->a[i],M->b[i-1],errorestimate);
    errorestimate*=M->b[i]*M->b[i]/(M->b[i]*M->b[i]+4*Gamma*Gamma);
    while(OperatorPsiCC(H,v[vi],&(v[vip])))
    { /*Allocation of Memory failed */
      printf("Allocation of new memory failed in LanczosTriDiagonalize for the function OperatorPsi\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",epsilon,epsilon*2.0);
      fflush(stdout);
      epsilon*=2.0;
      if(v[vi].NBlok==1)AllocFailed();
      ComplexWaveFunctionTruncateBasis(&(v[vi]),epsilon/errorestimate);
    }
    printf("BasisSetSize: %12i ",v[vip].N);
    ComplexWaveFunctionDotProduct(v[vip],v[vi],&M->a[i+1],&junk);
    while(ComplexWaveFunctionAddMultipliedUnionBasis(&(v[vip]),v[vim],-M->b[i]))
    { /*Allocation of Memory failed */
      printf("Allocation of new memory failed in LanczosTriDiagonalize for the function RealWaveFunctionAddMultipliedUnionBasis(&(v[vip]),v[vim],-M->b[i])\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",epsilon,epsilon*2.0);
      fflush(stdout);
      epsilon*=2.0;
      temp=ComplexWaveFunctionCalculateNorm(v[vip]);
      if(v[vip].NBlok==1)AllocFailed();
      ComplexWaveFunctionTruncateBasis(&(v[vip]),epsilon*temp/errorestimate);
    }
    ComplexWaveFunctionFree(&v[vim]);
  }
  printf("             | M.a[%3i] M.b[%3i]: %12.5E %12.5E | Err: %12.5E\n\n",i,i-1,M->a[i],M->b[i-1],errorestimate);
  if(M->N>1)
  {
    ComplexWaveFunctionFree(&v[vi]);
    ComplexWaveFunctionFree(&v[vip]);
  }
  fflush(stdout);
  return 0;
}

uint32_t LanczosTriDiagonalizeMC(const OperatorType H, WaveFunctionType psi, double epsilon, TriDiagonalMatrixType *const M, double Gamma)
{
  printf("\nStart of LanczosTriDiagonalizeMC\n\n");
  fflush(stdout);
  WaveFunctionType v[3];
  uint32_t i,vi,vip,vim;
  double temp,errorestimate=1.0,junk;
  /* copy the name of psi to v[0] in order to do the loop more cleanely  */
  v[0].N=psi.N;
  v[0].NBlok=psi.NBlok;
  v[0].NFermion=psi.NFermion;
  v[0].NBoson=psi.NBoson;
  v[0].NWord=psi.NWord;
  v[0].NBitsKey=psi.NBitsKey;
  v[0].ValueRe=psi.ValueRe;
  v[0].ValueIm=psi.ValueIm;
  v[0].Determinant=psi.Determinant;
  v[0].Next=psi.Next;
  v[0].KeyList=psi.KeyList;
  v[0].QComplex=psi.QComplex;
  /* finished copying */
  M->a[0]=0;
  M->b[0]=ComplexWaveFunctionCalculateNorm(v[0]);
  ComplexWaveFunctionMultiply(v[0],(double)1.0/M->b[0]);
  printf("BasisSetSize: %12i ",v[0].N);
  ComplexWaveFunctionTruncateBasis(&(v[0]),epsilon);
  printf("%12i\n",v[0].N);
  while(OperatorPsiMC(H,v[0],&(v[1])))
  { /* Allocation of Memory failed */
    printf("Allocation of new memory failed in LanczosTriDiagonalize for the function OperatorPsi\n");
    printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",epsilon,epsilon*2.0);
    fflush(stdout);
    epsilon*=2.0;
    if(v[0].NBlok==1)AllocFailed();
    ComplexWaveFunctionTruncateBasis(&(v[0]),epsilon);
  }
  ComplexWaveFunctionDotProduct(v[1],v[0],&M->a[1],&junk);
  printf("BasisSetSize: %12i ",v[1].N);
  for(i=1;i<M->N;i++)
  {
    vip=(i+1)%3;
    vi =(i  )%3;
    vim=(i-1)%3;
    while(ComplexWaveFunctionAddMultipliedUnionBasis(&(v[vi]),v[vim],-M->a[i]))
    { /*Allocation of Memory failed */
      printf("Allocation of new memory failed in LanczosTriDiagonalize for the function RealWaveFunctionAddMultipliedUnionBasis(&(v[vi]),v[vim],-M->a[i])\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",epsilon,epsilon*2.0);
      fflush(stdout);
      epsilon*=2.0;
      temp=ComplexWaveFunctionCalculateNorm(v[vi]);
      if(v[vi].NBlok==1)AllocFailed();
      ComplexWaveFunctionTruncateBasis(&(v[vi]),epsilon*temp/errorestimate);
    }
    M->b[i]=ComplexWaveFunctionCalculateNorm(v[vi]);
    ComplexWaveFunctionMultiply(v[vi],(double)1.0/M->b[i]);
    ComplexWaveFunctionTruncateBasis(&(v[vi]),epsilon/errorestimate);
    printf("%12i ",v[vi].N);
    printf("| M.a[%3i] M.b[%3i]: %12.5E %12.5E | Err: %12.5E\n",i,i-1,M->a[i],M->b[i-1],errorestimate);
    errorestimate*=M->b[i]*M->b[i]/(M->b[i]*M->b[i]+4*Gamma*Gamma);
    while(OperatorPsiMC(H,v[vi],&(v[vip])))
    { /*Allocation of Memory failed */
      printf("Allocation of new memory failed in LanczosTriDiagonalize for the function OperatorPsi\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",epsilon,epsilon*2.0);
      fflush(stdout);
      epsilon*=2.0;
      if(v[vi].NBlok==1)AllocFailed();
      ComplexWaveFunctionTruncateBasis(&(v[vi]),epsilon/errorestimate);
    }
    printf("BasisSetSize: %12i ",v[vip].N);
    ComplexWaveFunctionDotProduct(v[vip],v[vi],&M->a[i+1],&junk);
    while(ComplexWaveFunctionAddMultipliedUnionBasis(&(v[vip]),v[vim],-M->b[i]))
    { /*Allocation of Memory failed */
      printf("Allocation of new memory failed in LanczosTriDiagonalize for the function RealWaveFunctionAddMultipliedUnionBasis(&(v[vip]),v[vim],-M->b[i])\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",epsilon,epsilon*2.0);
      fflush(stdout);
      epsilon*=2.0;
      temp=ComplexWaveFunctionCalculateNorm(v[vip]);
      if(v[vip].NBlok==1)AllocFailed();
      ComplexWaveFunctionTruncateBasis(&(v[vip]),epsilon*temp/errorestimate);
    }
    ComplexWaveFunctionFree(&v[vim]);
  }
  printf("             | M.a[%3i] M.b[%3i]: %12.5E %12.5E | Err: %12.5E\n\n",i,i-1,M->a[i],M->b[i-1],errorestimate);
  if(M->N>1)
  {
    ComplexWaveFunctionFree(&v[vi]);
    ComplexWaveFunctionFree(&v[vip]);
  }
  fflush(stdout);
  return 0;
}



uint32_t LanczosTriDiagonalizePlusDensityRR(const OperatorType H, WaveFunctionType psi, double epsilon, TriDiagonalMatrixType *const M, MatrixType *const Rho, double Gamma)
{
  printf("\nStart of LanczosTriDiagonalizePlusDensityRR\n\n");
  fflush(stdout);
  WaveFunctionType v[3];
  uint32_t i,vi,vip,vim;
  double temp,errorestimate=1.0;
  double QRe[3],QIm[3],PRe[3],PIm[3],ValEst[3];
  QRe[0]=1;
  QRe[1]=0;
  QRe[2]=0;
  QIm[0]=0;
  QIm[1]=0;
  QIm[2]=0;
  PRe[0]=0;
  PRe[1]=0;
  PRe[2]=1;
  PIm[0]=0;
  PIm[1]=0;
  PIm[2]=0;
  ValEst[0]=0.0;
  ValEst[1]=0.0;
  ValEst[2]=0.0;
  /* copy the name of psi to v[0] in order to do the loop more cleanely  */
  v[0].N=psi.N;
  v[0].NBlok=psi.NBlok;
  v[0].NFermion=psi.NFermion;
  v[0].NBoson=psi.NBoson;
  v[0].NWord=psi.NWord;
  v[0].NBitsKey=psi.NBitsKey;
  v[0].ValueRe=psi.ValueRe;
  v[0].Determinant=psi.Determinant;
  v[0].Next=psi.Next;
  v[0].KeyList=psi.KeyList;
  v[0].QComplex=psi.QComplex;
  /* finished copying */
  M->a[0]=0;
  M->b[0]=RealWaveFunctionCalculateNorm(v[0]);
  if(M->b[0]==0)
  {
    RealWaveFunctionMultiply(v[0],(double)0.0);
  }
  else
  {
    RealWaveFunctionMultiply(v[0],(double)1.0/M->b[0]);
  }
  CalculateDensity(Rho->VRe[0],v[0]);
  printf("BasisSetSize: %12i ",v[0].N);
  RealWaveFunctionTruncateBasis(&(v[0]),epsilon);
  printf("%12i\n",v[0].N);
  while(OperatorPsiRR(H,v[0],&(v[1])))
  { /* Allocation of Memory failed */
    printf("Allocation of new memory failed in LanczosTriDiagonalize for the function OperatorPsi\n");
    printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",epsilon,epsilon*2.0);
    fflush(stdout);
    epsilon*=2.0;
    if(v[0].NBlok==1)AllocFailed();
    RealWaveFunctionTruncateBasis(&(v[0]),epsilon);
  }
  M->a[1]=RealWaveFunctionDotProduct(v[1],v[0]);
  printf("BasisSetSize: %12i ",v[1].N);
  for(i=1;i<M->N;i++)
  {
    vip=(i+1)%3;
    vi =(i  )%3;
    vim=(i-1)%3;
    while(RealWaveFunctionAddMultipliedUnionBasis(&(v[vi]),v[vim],-M->a[i]))
    { /*Allocation of Memory failed */
      printf("Allocation of new memory failed in LanczosTriDiagonalize for the function RealWaveFunctionAddMultipliedUnionBasis(&(v[vi]),v[vim],-M->a[i])\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",epsilon,epsilon*2.0);
      fflush(stdout);
      epsilon*=2.0;
      temp=RealWaveFunctionCalculateNorm(v[vi]);
      if(v[vi].NBlok==1)AllocFailed();
      RealWaveFunctionTruncateBasis(&(v[vi]),epsilon*temp/errorestimate);
    }
    M->b[i]=RealWaveFunctionCalculateNorm(v[vi]);
    if(M->b[i]==0)
    {
      RealWaveFunctionMultiply(v[vi],(double)0.0);
    }
    else
    {
      RealWaveFunctionMultiply(v[vi],(double)1.0/M->b[i]);
    }
    CalculateDensity(Rho->VRe[i],v[vi]);
    RealWaveFunctionTruncateBasis(&(v[vi]),epsilon/errorestimate);
    printf("%12i ",v[vi].N);
    printf("| M.a[%3i] M.b[%3i]: %12.5E %12.5E | Err: %12.5E\n",i,i-1,M->a[i],M->b[i-1],errorestimate);
    /*
    QRe[vi]=M->b[i-1]*QRe[vip]-0.5*Gamma*QIm[vim];
    QIm[vi]=M->b[i-1]*QIm[vip]+0.5*Gamma*QRe[vim];
    PRe[vi]=M->b[i-1]*PRe[vip]-0.5*Gamma*PIm[vim];
    PIm[vi]=M->b[i-1]*PIm[vip]+0.5*Gamma*PRe[vim];
    ValEst[vi]=fabs((PIm[vi]*QRe[vi]-PRe[vi]*QIm[vi])/(QIm[vi]*QIm[vi]+QRe[vi]*QRe[vi]));
    if(i>5)errorestimate=fabs(2*(ValEst[vi]-ValEst[vim])/(Gamma));
    */
    errorestimate*=M->b[i]*M->b[i]/(M->b[i]*M->b[i]+4*Gamma*Gamma);
    while(OperatorPsiRR(H,v[vi],&(v[vip])))
    { /*Allocation of Memory failed */
      printf("Allocation of new memory failed in LanczosTriDiagonalize for the function OperatorPsi\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",epsilon,epsilon*2.0);
      fflush(stdout);
      epsilon*=2.0;
      if(v[vi].NBlok==1)AllocFailed();
      RealWaveFunctionTruncateBasis(&(v[vi]),epsilon/errorestimate);
    }
    printf("BasisSetSize: %12i ",v[vip].N);
    M->a[i+1]=RealWaveFunctionDotProduct(v[vip],v[vi]);
    while(RealWaveFunctionAddMultipliedUnionBasis(&(v[vip]),v[vim],-M->b[i]))
    { /*Allocation of Memory failed */
      printf("Allocation of new memory failed in LanczosTriDiagonalize for the function RealWaveFunctionAddMultipliedUnionBasis(&(v[vip]),v[vim],-M->b[i])\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",epsilon,epsilon*2.0);
      fflush(stdout);
      epsilon*=2.0;
      temp=RealWaveFunctionCalculateNorm(v[vip]);
      if(v[vip].NBlok==1)AllocFailed();
      RealWaveFunctionTruncateBasis(&(v[vip]),epsilon*temp/errorestimate);
    }
    RealWaveFunctionFree(&v[vim]);
  }
  printf("             | M.a[%3i] M.b[%3i]: %12.5E %12.5E | Err: %12.5E\n\n",i,i-1,M->a[i],M->b[i-1],errorestimate);
  if(M->N>1)
  {
    RealWaveFunctionFree(&v[vi]);
    RealWaveFunctionFree(&v[vip]);
  }
  fflush(stdout);
  return 0;
}

uint32_t LanczosTriDiagonalizePlusDensityCR(const OperatorType H, WaveFunctionType psi, double epsilon, TriDiagonalMatrixType *const M, MatrixType *const Rho, double Gamma)
{
  uint32_t retval;
  WaveFunctionType Cpsi;
  if(MakeComplexWaveFunctionFromReal(&Cpsi,psi)){printf("LanczosTriDiagonalizePlusDensityCR called with a real wave-function and a complex Hamiltonian.\nTurned the wave-function complex in step 1.\nBut failed.\n");return 1;}  
  retval=LanczosTriDiagonalizePlusDensityCC(H, Cpsi, epsilon, M, Rho, Gamma);
  MakeRealWaveFunctionFromComplex(Cpsi,&Cpsi);
  return retval;
}

uint32_t LanczosTriDiagonalizePlusDensityMR(const OperatorType H, WaveFunctionType psi, double epsilon, TriDiagonalMatrixType *const M, MatrixType *const Rho, double Gamma)
{
  uint32_t retval;
  WaveFunctionType Cpsi;
  if(MakeComplexWaveFunctionFromReal(&Cpsi,psi)){printf("LanczosTriDiagonalizePlusDensityMR called with a real wave-function and a mixed complex Hamiltonian.\nTurned the wave-function complex in step 1.\nBut failed.\n");return 1;}  
  retval=LanczosTriDiagonalizePlusDensityMC(H, Cpsi, epsilon, M, Rho, Gamma);
  MakeRealWaveFunctionFromComplex(Cpsi,&Cpsi);
  return retval;
}

uint32_t LanczosTriDiagonalizePlusDensityRC(const OperatorType H, WaveFunctionType psi, double epsilon, TriDiagonalMatrixType *const M, MatrixType *const Rho, double Gamma)
{
  printf("\nStart of LanczosTriDiagonalizePlusDensityRC\n\n");
  fflush(stdout);
  WaveFunctionType v[3];
  uint32_t i,vi,vip,vim;
  double temp,errorestimate=1.0,junk;
  /* copy the name of psi to v[0] in order to do the loop more cleanely  */
  v[0].N=psi.N;
  v[0].NBlok=psi.NBlok;
  v[0].NFermion=psi.NFermion;
  v[0].NBoson=psi.NBoson;
  v[0].NWord=psi.NWord;
  v[0].NBitsKey=psi.NBitsKey;
  v[0].ValueRe=psi.ValueRe;
  v[0].ValueIm=psi.ValueIm;
  v[0].Determinant=psi.Determinant;
  v[0].Next=psi.Next;
  v[0].KeyList=psi.KeyList;
  v[0].QComplex=psi.QComplex;
  /* finished copying */
  M->a[0]=0;
  M->b[0]=ComplexWaveFunctionCalculateNorm(v[0]);
  ComplexWaveFunctionMultiply(v[0],(double)1.0/M->b[0]);
  printf("BasisSetSize: %12i ",v[0].N);
  ComplexWaveFunctionTruncateBasis(&(v[0]),epsilon);
  printf("%12i\n",v[0].N);
  while(OperatorPsiRC(H,v[0],&(v[1])))
  { /* Allocation of Memory failed */
    printf("Allocation of new memory failed in LanczosTriDiagonalize for the function OperatorPsi\n");
    printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",epsilon,epsilon*2.0);
    fflush(stdout);
    epsilon*=2.0;
    if(v[0].NBlok==1)AllocFailed();
    ComplexWaveFunctionTruncateBasis(&(v[0]),epsilon);
  }
  ComplexWaveFunctionDotProduct(v[1],v[0],&M->a[1],&junk);
  printf("BasisSetSize: %12i ",v[1].N);
  for(i=1;i<M->N;i++)
  {
    vip=(i+1)%3;
    vi =(i  )%3;
    vim=(i-1)%3;
    while(ComplexWaveFunctionAddMultipliedUnionBasis(&(v[vi]),v[vim],-M->a[i]))
    { /*Allocation of Memory failed */
      printf("Allocation of new memory failed in LanczosTriDiagonalize for the function RealWaveFunctionAddMultipliedUnionBasis(&(v[vi]),v[vim],-M->a[i])\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",epsilon,epsilon*2.0);
      fflush(stdout);
      epsilon*=2.0;
      temp=ComplexWaveFunctionCalculateNorm(v[vi]);
      if(v[vi].NBlok==1)AllocFailed();
      ComplexWaveFunctionTruncateBasis(&(v[vi]),epsilon*temp/errorestimate);
    }
    M->b[i]=ComplexWaveFunctionCalculateNorm(v[vi]);
    ComplexWaveFunctionMultiply(v[vi],(double)1.0/M->b[i]);
    ComplexWaveFunctionTruncateBasis(&(v[vi]),epsilon/errorestimate);
    printf("%12i ",v[vi].N);
    printf("| M.a[%3i] M.b[%3i]: %12.5E %12.5E | Err: %12.5E\n",i,i-1,M->a[i],M->b[i-1],errorestimate);
    errorestimate*=M->b[i]*M->b[i]/(M->b[i]*M->b[i]+4*Gamma*Gamma);
    while(OperatorPsiRC(H,v[vi],&(v[vip])))
    { /*Allocation of Memory failed */
      printf("Allocation of new memory failed in LanczosTriDiagonalize for the function OperatorPsi\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",epsilon,epsilon*2.0);
      fflush(stdout);
      epsilon*=2.0;
      if(v[vi].NBlok==1)AllocFailed();
      ComplexWaveFunctionTruncateBasis(&(v[vi]),epsilon/errorestimate);
    }
    printf("BasisSetSize: %12i ",v[vip].N);
    ComplexWaveFunctionDotProduct(v[vip],v[vi],&M->a[i+1],&junk);
    while(ComplexWaveFunctionAddMultipliedUnionBasis(&(v[vip]),v[vim],-M->b[i]))
    { /*Allocation of Memory failed */
      printf("Allocation of new memory failed in LanczosTriDiagonalize for the function RealWaveFunctionAddMultipliedUnionBasis(&(v[vip]),v[vim],-M->b[i])\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",epsilon,epsilon*2.0);
      fflush(stdout);
      epsilon*=2.0;
      temp=ComplexWaveFunctionCalculateNorm(v[vip]);
      if(v[vip].NBlok==1)AllocFailed();
      ComplexWaveFunctionTruncateBasis(&(v[vip]),epsilon*temp/errorestimate);
    }
    ComplexWaveFunctionFree(&v[vim]);
  }
  printf("             | M.a[%3i] M.b[%3i]: %12.5E %12.5E | Err: %12.5E\n\n",i,i-1,M->a[i],M->b[i-1],errorestimate);
  if(M->N>1)
  {
    ComplexWaveFunctionFree(&v[vi]);
    ComplexWaveFunctionFree(&v[vip]);
  }
  fflush(stdout);
  return 0;
}

uint32_t LanczosTriDiagonalizePlusDensityCC(const OperatorType H, WaveFunctionType psi, double epsilon, TriDiagonalMatrixType *const M, MatrixType *const Rho, double Gamma)
{
  printf("\nStart of LanczosTriDiagonalizePlusDensityCC\n\n");
  fflush(stdout);
  WaveFunctionType v[3];
  uint32_t i,vi,vip,vim;
  double temp,errorestimate=1.0,junk;
  /* copy the name of psi to v[0] in order to do the loop more cleanely  */
  v[0].N=psi.N;
  v[0].NBlok=psi.NBlok;
  v[0].NFermion=psi.NFermion;
  v[0].NBoson=psi.NBoson;
  v[0].NWord=psi.NWord;
  v[0].NBitsKey=psi.NBitsKey;
  v[0].ValueRe=psi.ValueRe;
  v[0].ValueIm=psi.ValueIm;
  v[0].Determinant=psi.Determinant;
  v[0].Next=psi.Next;
  v[0].KeyList=psi.KeyList;
  v[0].QComplex=psi.QComplex;
  /* finished copying */
  M->a[0]=0;
  M->b[0]=ComplexWaveFunctionCalculateNorm(v[0]);
  ComplexWaveFunctionMultiply(v[0],(double)1.0/M->b[0]);
  printf("BasisSetSize: %12i ",v[0].N);
  ComplexWaveFunctionTruncateBasis(&(v[0]),epsilon);
  printf("%12i\n",v[0].N);
  while(OperatorPsiCC(H,v[0],&(v[1])))
  { /* Allocation of Memory failed */
    printf("Allocation of new memory failed in LanczosTriDiagonalize for the function OperatorPsi\n");
    printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",epsilon,epsilon*2.0);
    fflush(stdout);
    epsilon*=2.0;
    if(v[0].NBlok==1)AllocFailed();
    ComplexWaveFunctionTruncateBasis(&(v[0]),epsilon);
  }
  ComplexWaveFunctionDotProduct(v[1],v[0],&M->a[1],&junk);
  printf("BasisSetSize: %12i ",v[1].N);
  for(i=1;i<M->N;i++)
  {
    vip=(i+1)%3;
    vi =(i  )%3;
    vim=(i-1)%3;
    while(ComplexWaveFunctionAddMultipliedUnionBasis(&(v[vi]),v[vim],-M->a[i]))
    { /*Allocation of Memory failed */
      printf("Allocation of new memory failed in LanczosTriDiagonalize for the function RealWaveFunctionAddMultipliedUnionBasis(&(v[vi]),v[vim],-M->a[i])\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",epsilon,epsilon*2.0);
      fflush(stdout);
      epsilon*=2.0;
      temp=ComplexWaveFunctionCalculateNorm(v[vi]);
      if(v[vi].NBlok==1)AllocFailed();
      ComplexWaveFunctionTruncateBasis(&(v[vi]),epsilon*temp/errorestimate);
    }
    M->b[i]=ComplexWaveFunctionCalculateNorm(v[vi]);
    ComplexWaveFunctionMultiply(v[vi],(double)1.0/M->b[i]);
    ComplexWaveFunctionTruncateBasis(&(v[vi]),epsilon/errorestimate);
    printf("%12i ",v[vi].N);
    printf("| M.a[%3i] M.b[%3i]: %12.5E %12.5E | Err: %12.5E\n",i,i-1,M->a[i],M->b[i-1],errorestimate);
    errorestimate*=M->b[i]*M->b[i]/(M->b[i]*M->b[i]+4*Gamma*Gamma);
    while(OperatorPsiCC(H,v[vi],&(v[vip])))
    { /*Allocation of Memory failed */
      printf("Allocation of new memory failed in LanczosTriDiagonalize for the function OperatorPsi\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",epsilon,epsilon*2.0);
      fflush(stdout);
      epsilon*=2.0;
      if(v[vi].NBlok==1)AllocFailed();
      ComplexWaveFunctionTruncateBasis(&(v[vi]),epsilon/errorestimate);
    }
    printf("BasisSetSize: %12i ",v[vip].N);
    ComplexWaveFunctionDotProduct(v[vip],v[vi],&M->a[i+1],&junk);
    while(ComplexWaveFunctionAddMultipliedUnionBasis(&(v[vip]),v[vim],-M->b[i]))
    { /*Allocation of Memory failed */
      printf("Allocation of new memory failed in LanczosTriDiagonalize for the function RealWaveFunctionAddMultipliedUnionBasis(&(v[vip]),v[vim],-M->b[i])\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",epsilon,epsilon*2.0);
      fflush(stdout);
      epsilon*=2.0;
      temp=ComplexWaveFunctionCalculateNorm(v[vip]);
      if(v[vip].NBlok==1)AllocFailed();
      ComplexWaveFunctionTruncateBasis(&(v[vip]),epsilon*temp/errorestimate);
    }
    ComplexWaveFunctionFree(&v[vim]);
  }
  printf("             | M.a[%3i] M.b[%3i]: %12.5E %12.5E | Err: %12.5E\n\n",i,i-1,M->a[i],M->b[i-1],errorestimate);
  if(M->N>1)
  {
    ComplexWaveFunctionFree(&v[vi]);
    ComplexWaveFunctionFree(&v[vip]);
  }
  fflush(stdout);
  return 0;
}

uint32_t LanczosTriDiagonalizePlusDensityMC(const OperatorType H, WaveFunctionType psi, double epsilon, TriDiagonalMatrixType *const M, MatrixType *const Rho, double Gamma)
{
  printf("\nStart of LanczosTriDiagonalizePlusDensityMC\n\n");
  fflush(stdout);
  WaveFunctionType v[3];
  uint32_t i,vi,vip,vim;
  double temp,errorestimate=1.0,junk;
  /* copy the name of psi to v[0] in order to do the loop more cleanely  */
  v[0].N=psi.N;
  v[0].NBlok=psi.NBlok;
  v[0].NFermion=psi.NFermion;
  v[0].NBoson=psi.NBoson;
  v[0].NWord=psi.NWord;
  v[0].NBitsKey=psi.NBitsKey;
  v[0].ValueRe=psi.ValueRe;
  v[0].ValueIm=psi.ValueIm;
  v[0].Determinant=psi.Determinant;
  v[0].Next=psi.Next;
  v[0].KeyList=psi.KeyList;
  v[0].QComplex=psi.QComplex;
  /* finished copying */
  M->a[0]=0;
  M->b[0]=ComplexWaveFunctionCalculateNorm(v[0]);
  ComplexWaveFunctionMultiply(v[0],(double)1.0/M->b[0]);
  printf("BasisSetSize: %12i ",v[0].N);
  ComplexWaveFunctionTruncateBasis(&(v[0]),epsilon);
  printf("%12i\n",v[0].N);
  while(OperatorPsiMC(H,v[0],&(v[1])))
  { /* Allocation of Memory failed */
    printf("Allocation of new memory failed in LanczosTriDiagonalize for the function OperatorPsi\n");
    printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",epsilon,epsilon*2.0);
    fflush(stdout);
    epsilon*=2.0;
    if(v[0].NBlok==1)AllocFailed();
    ComplexWaveFunctionTruncateBasis(&(v[0]),epsilon);
  }
  ComplexWaveFunctionDotProduct(v[1],v[0],&M->a[1],&junk);
  printf("BasisSetSize: %12i ",v[1].N);
  for(i=1;i<M->N;i++)
  {
    vip=(i+1)%3;
    vi =(i  )%3;
    vim=(i-1)%3;
    while(ComplexWaveFunctionAddMultipliedUnionBasis(&(v[vi]),v[vim],-M->a[i]))
    { /*Allocation of Memory failed */
      printf("Allocation of new memory failed in LanczosTriDiagonalize for the function RealWaveFunctionAddMultipliedUnionBasis(&(v[vi]),v[vim],-M->a[i])\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",epsilon,epsilon*2.0);
      fflush(stdout);
      epsilon*=2.0;
      temp=ComplexWaveFunctionCalculateNorm(v[vi]);
      if(v[vi].NBlok==1)AllocFailed();
      ComplexWaveFunctionTruncateBasis(&(v[vi]),epsilon*temp/errorestimate);
    }
    M->b[i]=ComplexWaveFunctionCalculateNorm(v[vi]);
    ComplexWaveFunctionMultiply(v[vi],(double)1.0/M->b[i]);
    ComplexWaveFunctionTruncateBasis(&(v[vi]),epsilon/errorestimate);
    printf("%12i ",v[vi].N);
    printf("| M.a[%3i] M.b[%3i]: %12.5E %12.5E | Err: %12.5E\n",i,i-1,M->a[i],M->b[i-1],errorestimate);
    errorestimate*=M->b[i]*M->b[i]/(M->b[i]*M->b[i]+4*Gamma*Gamma);
    while(OperatorPsiMC(H,v[vi],&(v[vip])))
    { /*Allocation of Memory failed */
      printf("Allocation of new memory failed in LanczosTriDiagonalize for the function OperatorPsi\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",epsilon,epsilon*2.0);
      fflush(stdout);
      epsilon*=2.0;
      if(v[vi].NBlok==1)AllocFailed();
      ComplexWaveFunctionTruncateBasis(&(v[vi]),epsilon/errorestimate);
    }
    printf("BasisSetSize: %12i ",v[vip].N);
    ComplexWaveFunctionDotProduct(v[vip],v[vi],&M->a[i+1],&junk);
    while(ComplexWaveFunctionAddMultipliedUnionBasis(&(v[vip]),v[vim],-M->b[i]))
    { /*Allocation of Memory failed */
      printf("Allocation of new memory failed in LanczosTriDiagonalize for the function RealWaveFunctionAddMultipliedUnionBasis(&(v[vip]),v[vim],-M->b[i])\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",epsilon,epsilon*2.0);
      fflush(stdout);
      epsilon*=2.0;
      temp=ComplexWaveFunctionCalculateNorm(v[vip]);
      if(v[vip].NBlok==1)AllocFailed();
      ComplexWaveFunctionTruncateBasis(&(v[vip]),epsilon*temp/errorestimate);
    }
    ComplexWaveFunctionFree(&v[vim]);
  }
  printf("             | M.a[%3i] M.b[%3i]: %12.5E %12.5E | Err: %12.5E\n\n",i,i-1,M->a[i],M->b[i-1],errorestimate);
  if(M->N>1)
  {
    ComplexWaveFunctionFree(&v[vi]);
    ComplexWaveFunctionFree(&v[vip]);
  }
  fflush(stdout);
  return 0;
}






uint32_t LanczosGroundState(const OperatorType H, WaveFunctionType *const psi, double *epsilon, const double zero, double shiftin, WaveFunctionType *OtherPsi, uint32_t NOtherPsi)
{
  if(psi->QComplex)
  {
    if(H.QComplex==0)
    {
      return LanczosGroundStateRC(H, psi, epsilon, zero, shiftin, OtherPsi, NOtherPsi);
    }
    else if(H.QComplex==1)
    {
      return LanczosGroundStateCC(H, psi, epsilon, zero, shiftin, OtherPsi, NOtherPsi);
    }
    else
    {
      return LanczosGroundStateMC(H, psi, epsilon, zero, shiftin, OtherPsi, NOtherPsi);
    }
  }
  else
  {
    if(H.QComplex==0)
    {
      return LanczosGroundStateRR(H, psi, epsilon, zero, shiftin, OtherPsi, NOtherPsi);
    }
    else if(H.QComplex==1)
    {
      if(MakeComplexWaveFunctionFromReal(psi, *psi)) {printf("MakeComplexWaveFunctionFromReal in LanczosGroundState failed\n");return 1;}
      return LanczosGroundStateCC(H, psi, epsilon, zero, shiftin, OtherPsi, NOtherPsi);
    }
    else
    {
      if(MakeComplexWaveFunctionFromReal(psi, *psi)) {printf("MakeComplexWaveFunctionFromReal in LanczosGroundState failed\n");return 1;}
      return LanczosGroundStateMC(H, psi, epsilon, zero, shiftin, OtherPsi, NOtherPsi);
    }
  }
}


uint32_t LanczosTriDiagonalize(const OperatorType H, WaveFunctionType psi, double epsilon, TriDiagonalMatrixType *const M, double Gamma)
{
  WaveFunctionType Cpsi;
  if(psi.QComplex)
  {
    if(H.QComplex==0)
    {
      return LanczosTriDiagonalizeRC(H, psi, epsilon, M, Gamma);
    }
    else if(H.QComplex==1)
    {
      return LanczosTriDiagonalizeCC(H, psi, epsilon, M, Gamma);
    }
    else
    {
      return LanczosTriDiagonalizeMC(H, psi, epsilon, M, Gamma);
    }
  }
  else
  {
    if(H.QComplex==0)
    {
      return LanczosTriDiagonalizeRR(H, psi, epsilon, M, Gamma);
    }
    else if(H.QComplex==1)    
    {
      if(MakeComplexWaveFunctionFromReal(&Cpsi, psi)){printf("MakeComplexWaveFunctionFromReal failed in LanczosTriDiagonalize\n");return 1;}
      return LanczosTriDiagonalizeCC(H, Cpsi, epsilon, M, Gamma);
    }
    else
    {
      if(MakeComplexWaveFunctionFromReal(&Cpsi, psi)){printf("MakeComplexWaveFunctionFromReal failed in LanczosTriDiagonalize\n");return 1;}
      return LanczosTriDiagonalizeMC(H, Cpsi, epsilon, M, Gamma);
    }
  }
}


uint32_t LanczosTriDiagonalizePlusDensity(const OperatorType H, WaveFunctionType psi, double epsilon, TriDiagonalMatrixType *const M, MatrixType *const Rho, double Gamma)
{
  WaveFunctionType Cpsi;
  if(psi.QComplex)
  {
    if(H.QComplex==0)
    {
      return LanczosTriDiagonalizePlusDensityRC(H, psi, epsilon, M, Rho, Gamma);
    }
    else if(H.QComplex==1)
    {
      return LanczosTriDiagonalizePlusDensityCC(H, psi, epsilon, M, Rho, Gamma);
    }
    else
    {
      return LanczosTriDiagonalizePlusDensityMC(H, psi, epsilon, M, Rho, Gamma);
    }
  }
  else
  {
    if(H.QComplex==0)
    {
      return LanczosTriDiagonalizePlusDensityRR(H, psi, epsilon, M, Rho, Gamma);
    }
    else if(H.QComplex==1)    
    {
      if(MakeComplexWaveFunctionFromReal(&Cpsi, psi)){printf("MakeComplexWaveFunctionFromReal failed in LanczosTriDiagonalize\n");return 1;}
      return LanczosTriDiagonalizePlusDensityCC(H, Cpsi, epsilon, M, Rho, Gamma);
    }
    else
    {
      if(MakeComplexWaveFunctionFromReal(&Cpsi, psi)){printf("MakeComplexWaveFunctionFromReal failed in LanczosTriDiagonalize\n");return 1;}
      return LanczosTriDiagonalizePlusDensityMC(H, Cpsi, epsilon, M, Rho, Gamma);
    }
  }
}








