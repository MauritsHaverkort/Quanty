









/* =========================================== old Lanczos functions ============================================= */




/*
*   Eigenstates of a 2 by 2 matrix
*   The lowest eigen vector and energy can be found in Energy[0] and {Functions[0][0],Functions[0][1]}
*   Exept when the matrix on entry is already diagonal or close to diagonal |b1|<(DBL_EPSILON*|a0-a1|)
*     in that case the returned eigen-system are the entered basis functions
*/
void Eigensystem2(const double ain0, const double ain1, const double b1, double Energy[], double Functions[][2])
{
  double shift,a0,a1;
  shift=(ain0+ain1)/2.0;
  a0=ain0-shift;
  a1=ain1-shift;
  if(fabs(b1)<(10*(DBL_EPSILON*fabs(a0-a1)+DBL_MIN)))
  {
    if(a0<a1)
    {
      Energy[0]=a0;
      Energy[1]=a1;
      Functions[0][0]=1;
      Functions[0][1]=0;
      Functions[1][0]=0;
      Functions[1][1]=1;
    }
    else
    {
      Energy[0]=a1;
      Energy[1]=a0;
      Functions[0][0]=0;
      Functions[0][1]=1;
      Functions[1][0]=1;
      Functions[1][1]=0;
    }
  }
  else
  {
    Energy[0]=0.5*(a0+a1-sqrt(Sqr(a0-a1)+4*Sqr(b1)));
    Energy[1]=0.5*(a0+a1+sqrt(Sqr(a0-a1)+4*Sqr(b1)));
    /* use only the vector for wich |E-a| is largest */
    if(fabs((Energy[0]-a1))>fabs(Energy[1]-a1))
    {/* fun 0 stable */
      Functions[0][0]=(Energy[0]-a1)/sqrt(Sqr(Energy[0]-a1)+Sqr(b1));
      Functions[0][1]=(b1)/sqrt(Sqr(Energy[0]-a1)+Sqr(b1));
      Functions[1][1]=Functions[0][0];
      Functions[1][0]=-Functions[0][1];
    }
    else
    {/* fun 1 stable */
      Functions[1][0]=(Energy[1]-a1)/sqrt(Sqr(Energy[1]-a1)+Sqr(b1));
      Functions[1][1]=(b1)/sqrt(Sqr(Energy[1]-a1)+Sqr(b1));
      Functions[0][1]=Functions[1][0];
      Functions[0][0]=-Functions[1][1];
    }
  }
  Energy[0]+=shift;
  Energy[1]+=shift;
  return;
}

uint32_t Eigensystem3(const double a0, const double a1, const double a2, const double b1, const double b2, double Energy[3], double Functions[3][3])
{
  double e[3]={b1,b2,0};
  double g, r, p, f, b, s, c, t, shift;
  int nIter;
  int m,l,i,k;

  shift=(a0+a1+a2)/3.0;
  Energy[0]=a0-shift;
  Energy[1]=a1-shift;
  Energy[2]=a2-shift;
  Functions[0][0]=1; Functions[0][1]=0; Functions[0][2]=0;
  Functions[1][0]=0; Functions[1][1]=1; Functions[1][2]=0;
  Functions[2][0]=0; Functions[2][1]=0; Functions[2][2]=1;

  for (l=0; l < 2; l++)
  {
    nIter = 0;
    while (1)
    {
      // Check for convergence and exit iteration loop if off-diagonal
      // element e(l) is zero
      for (m=l; m <= 1; m++)
      {
        g = fabs(Energy[m])+fabs(Energy[m+1]);
        if (fabs(e[m]) + g == g)
          break;
      }
      if (m == l)
        break;
      
      if (nIter++ >= 30)
        return 1;

      // Calculate g = d_m - k
      g = (Energy[l+1] - Energy[l]) / (e[l] + e[l]);
      r = sqrt(Sqr(g) + 1.0);
      if (g > 0)
        g = Energy[m] - Energy[l] + e[l]/(g + r);
      else
        g = Energy[m] - Energy[l] + e[l]/(g - r);

      s = c = 1.0;
      p = 0.0;
      for (i=m-1; i >= l; i--)
      {
        f = s * e[i];
        b = c * e[i];
        if (fabs(f) > fabs(g))
        {
          c      = g / f;
          r      = sqrt(Sqr(c) + 1.0);
          e[i+1] = f * r;
          c     *= (s = 1.0/r);
        }
        else
        {
          s      = f / g;
          r      = sqrt(Sqr(s) + 1.0);
          e[i+1] = g * r;
          s     *= (c = 1.0/r);
        }
        
        g = Energy[i+1] - p;
        r = (Energy[i] - g)*s + 2.0*c*b;
        p = s * r;
        Energy[i+1] = g + p;
        g = c*r - b;

        for (k=0; k < 3; k++)
        {
          t = Functions[i+1][k];
          Functions[i+1][k] = s*Functions[i][k] + c*t;
          Functions[i][k]   = c*Functions[i][k] - s*t;
        }

      }
      Energy[l] -= p;
      e[l]  = g;
      e[m]  = 0.0;
    }
  }
  Energy[0]+=shift;
  Energy[1]+=shift;
  Energy[2]+=shift;
  if((Energy[0]<=Energy[1])&&(Energy[1]<=Energy[2]))
  {
    return 0;
  }
  if((Energy[0]<=Energy[2])&&(Energy[2]<=Energy[1]))
  {
    r=Energy[1];
    Energy[1]=Energy[2];
    Energy[2]=r;
    e[0]=Functions[1][0]; e[1]=Functions[1][1]; e[2]=Functions[1][2];
    Functions[1][0]=Functions[2][0]; Functions[1][1]=Functions[2][1]; Functions[1][2]=Functions[2][2];
    Functions[2][0]=e[0]; Functions[2][1]=e[1]; Functions[2][2]=e[2];
    return 0;
  }
  if((Energy[1]<=Energy[0])&&(Energy[0]<=Energy[2]))
  {
    r=Energy[1];
    Energy[1]=Energy[0];
    Energy[0]=r;
    e[0]=Functions[1][0]; e[1]=Functions[1][1]; e[2]=Functions[1][2];
    Functions[1][0]=Functions[0][0]; Functions[1][1]=Functions[0][1]; Functions[1][2]=Functions[0][2];
    Functions[0][0]=e[0]; Functions[0][1]=e[1]; Functions[0][2]=e[2];
    return 0;
  }
  if((Energy[2]<=Energy[1])&&(Energy[1]<=Energy[0]))
  {
    r=Energy[2];
    Energy[2]=Energy[0];
    Energy[0]=r;
    e[0]=Functions[2][0]; e[1]=Functions[2][1]; e[2]=Functions[2][2];
    Functions[2][0]=Functions[0][0]; Functions[2][1]=Functions[0][1]; Functions[2][2]=Functions[0][2];
    Functions[0][0]=e[0]; Functions[0][1]=e[1]; Functions[0][2]=e[2];
    return 0;
  }
  if((Energy[1]<=Energy[2])&&(Energy[2]<=Energy[0]))
  {
    r=Energy[1];
    Energy[1]=Energy[2];
    Energy[2]=Energy[0];
    Energy[0]=r;
    e[0]=Functions[1][0]; e[1]=Functions[1][1]; e[2]=Functions[1][2];
    Functions[1][0]=Functions[2][0]; Functions[1][1]=Functions[2][1]; Functions[1][2]=Functions[2][2];
    Functions[2][0]=Functions[0][0]; Functions[2][1]=Functions[0][1]; Functions[2][2]=Functions[0][2];
    Functions[0][0]=e[0]; Functions[0][1]=e[1]; Functions[0][2]=e[2];
    return 0;
  }
  if((Energy[1]<=Energy[2])&&(Energy[2]<=Energy[0]))
  {
    r=Energy[2];
    Energy[2]=Energy[1];
    Energy[1]=Energy[0];
    Energy[0]=r;
    e[0]=Functions[2][0]; e[1]=Functions[2][1]; e[2]=Functions[2][2];
    Functions[2][0]=Functions[1][0]; Functions[2][1]=Functions[1][1]; Functions[2][2]=Functions[1][2];
    Functions[1][0]=Functions[0][0]; Functions[1][1]=Functions[0][1]; Functions[1][2]=Functions[0][2];
    Functions[0][0]=e[0]; Functions[0][1]=e[1]; Functions[0][2]=e[2];
    return 0;
  }
  return 2;
}



/*
*   Returns a single wave function orthogonal to the functions given in OhterPsi with the lowest energy, given
*     H         Hamiltonian
*     psi       Starting wave function
*     epsilon   Value to wich psi is truncated before a single excitation is added
*     zero      Smallest (positive) value that should be distinguised from zero
*     shiftin   Positive value larger then the largest eigenvalue
*     OtherPsi  Wave-functions the solution has to be orthogonal to
*     NOtherPsi Number of wave-functions to orthogonalize to
*
*
*   The values in psi are distroyed and the returned vector has the "same" storage location as the input wave function
*
*   algorithm:
*   while(notconverged)
*     psi = H psi
*     for(i=0;i<NOtherPsi;i++)
*     {
*       psi-=(OtherPsi_i . psi) OtherPsi_i
*     }
*     psi = psi/|psi|
*   end while
*   return psi
*/
void LanczosGroundStateFullRestarted(const OperatorType H, WaveFunctionType *const psi, double *epsilon, const double zero, double shiftin, WaveFunctionType *OtherPsi, uint32_t NOtherPsi)
{
  printf("Start of LanczosGroundState\n\n");
  fflush(stdout);
  WaveFunctionType Opsi;
  double energy=4*zero,energy1=2*zero,energy2;
  double shift;
  uint32_t cycles;
  while((energy-energy1)>(zero*(1.0+fabs(energy))))
  {
    cycles=0;
    energy=energy1;
    shift=-shiftin-energy*0.49;
    while(OperatorPsi(H,*psi,&Opsi))                           /* Newly allocate memory for Opsi both in Value, Determinant and Next */
    {/* Allocation Failed */
      printf("Allocation of new memory failed in LanczosGroundState for the function OperatorPsi\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",*epsilon,*epsilon*2.0);
      fflush(stdout);
      *epsilon*=2.0;
      if(psi->NBlok==1)AllocFailed();
      RealWaveFunctionTruncateBasis(psi,*epsilon);
     }
    printf("BasisSetSizeAfter Expansion/Truncation : %12i ",Opsi.N);
    RealWaveFunctionMultiply(*psi,shift);
    while(RealWaveFunctionAddUnionBasis(&Opsi,*psi))               /* Possibly expand allocated memory of Opsi */
    {/* Allocation Failed */
      printf("Allocation of new memory failed in LanczosGroundState for the function RealWaveFunctionAddUnionBasis\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",*epsilon,*epsilon*2.0);
      fflush(stdout);
      *epsilon*=2.0;
      RealWaveFunctionFreeValues(&Opsi);  
      RealWaveFunctionMultiply(*psi,1.0/shift);
      if(psi->NBlok==1)AllocFailed();
      RealWaveFunctionTruncateBasis(psi,*epsilon);
      OperatorPsi(H,*psi,&Opsi);
      RealWaveFunctionMultiply(*psi,shift);
    }
    RealWaveFunctionFree(psi);                                     /* Remove allocated memory for psi from memory */
    while(RealWaveFunctionOrthogonalize(&Opsi, OtherPsi, NOtherPsi))
    {/* Allocation Failed */
      printf("Allocation of new memory failed in LanczosGroundState for the function RealWaveFunctionOrthogonalize\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",*epsilon,*epsilon*2.0);
      fflush(stdout);
      *epsilon*=2.0;
      if(Opsi.NBlok==1)AllocFailed();
      RealWaveFunctionTruncateBasis(&Opsi,*epsilon);
    }
    RealWaveFunctionNormalize(&Opsi);
    RealWaveFunctionTruncateBasis(&Opsi,*epsilon);                 /* Possibly remove allocated memory for Opsi */
    printf("%12i\n",Opsi.N);
    while(OperatorPsiConserveBasis(H,Opsi,psi))                /* Newly allocate memory for psi in Value. Determinant and Next point to same positions as Opsi*/
    {/* Allocation Failed */
      printf("Allocation of new memory failed in LanczosGroundState for the function OperatorPsiConserveBasis\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",*epsilon,*epsilon*2.0);
      fflush(stdout);
      *epsilon*=2.0;
      if(Opsi.NBlok==1)AllocFailed();
      RealWaveFunctionTruncateBasis(&Opsi,*epsilon);
    }
    RealWaveFunctionMultiply(Opsi,shift);
    RealWaveFunctionAddSameBasis(*psi,Opsi);
    RealWaveFunctionNormalize(psi);
    energy1=2*zero;
    energy2=4*zero;
    while((energy2-energy1)>(zero*(1.0+fabs(energy1))))
    {
      cycles++;
      if((cycles&255)==255) {printf("Reached %i-th cycle with energy %22.15E\n",cycles,energy1); fflush(stdout);}
      OperatorPsiConserveBasisPreAllocated(H,*psi,&Opsi);
      RealWaveFunctionOrthogonalizeConserveBasis(&Opsi, OtherPsi, NOtherPsi);
      energy2=RealWaveFunctionDotProductSameBasis(*psi,Opsi);
      RealWaveFunctionMultiply(*psi,shift);
      RealWaveFunctionAddSameBasis(Opsi,*psi);
      RealWaveFunctionNormalize(&Opsi);
      OperatorPsiConserveBasisPreAllocated(H,Opsi,psi);
      RealWaveFunctionOrthogonalizeConserveBasis(psi, OtherPsi, NOtherPsi);
      energy1=RealWaveFunctionDotProductSameBasis(*psi,Opsi);
      RealWaveFunctionMultiply(Opsi,shift);
      RealWaveFunctionAddSameBasis(*psi,Opsi);
      RealWaveFunctionNormalize(psi);
    }
    RealWaveFunctionFreeValues(&Opsi);                              /* Free Opsi.ValueRe, but don't free Opsi.Determinant as this is the same space as where psi.Determinant points to */
    printf("Energy: %22.15E  ",energy1);
    printf("NumberOfLaczosCycles: %i\n",cycles);
  }
  printf("\n");
  fflush(stdout);
}

void LanczosGroundStateDouble(const OperatorType H, WaveFunctionType *const psi, double *epsilon, const double zero, double shiftin, WaveFunctionType *OtherPsi, uint32_t NOtherPsi)
{
  printf("Start of LanczosGroundState\n\n");
  fflush(stdout);
  WaveFunctionType Opsi;
  double energy=4*zero,energy1=2*zero,energy2;
  double shift,a0,a1,b1;
  uint32_t cycles;
  while((energy-energy1)>(0.5*zero*(1.0+fabs(energy))))
  {
    cycles=0;
    energy=energy1;
    shift=-shiftin-energy*0.49;
    while(OperatorPsi(H,*psi,&Opsi))                           /* Newly allocate memory for Opsi both in Value, Determinant and Next */
    {/* Allocation Failed */
      printf("Allocation of new memory failed in LanczosGroundState for the function OperatorPsi\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",*epsilon,*epsilon*2.0);
      fflush(stdout);
      *epsilon*=2.0;
      if(psi->NBlok==1)AllocFailed();
      RealWaveFunctionTruncateBasis(psi,*epsilon);
     }
    printf("BasisSetSizeAfter Expansion/Truncation : %12i ",Opsi.N);
    RealWaveFunctionMultiply(*psi,shift);
    while(RealWaveFunctionAddUnionBasis(&Opsi,*psi))               /* Possibly expand allocated memory of Opsi */
    {/* Allocation Failed */
      printf("Allocation of new memory failed in LanczosGroundState for the function RealWaveFunctionAddUnionBasis\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",*epsilon,*epsilon*2.0);
      fflush(stdout);
      *epsilon*=2.0;
      RealWaveFunctionFreeValues(&Opsi);  
      RealWaveFunctionMultiply(*psi,1.0/shift);
      if(psi->NBlok==1)AllocFailed();
      RealWaveFunctionTruncateBasis(psi,*epsilon);
      OperatorPsi(H,*psi,&Opsi);
      RealWaveFunctionMultiply(*psi,shift);
    }
    RealWaveFunctionFree(psi);                                     /* Remove allocated memory for psi from memory */
    while(RealWaveFunctionOrthogonalize(&Opsi, OtherPsi, NOtherPsi))
    {/* Allocation Failed */
      printf("Allocation of new memory failed in LanczosGroundState for the function RealWaveFunctionOrthogonalize\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",*epsilon,*epsilon*2.0);
      fflush(stdout);
      *epsilon*=2.0;
      if(Opsi.NBlok==1)AllocFailed();
      RealWaveFunctionTruncateBasis(&Opsi,*epsilon);
    }
    RealWaveFunctionNormalize(&Opsi);
    RealWaveFunctionTruncateBasis(&Opsi,*epsilon);                 /* Possibly remove allocated memory for Opsi */
    printf("%12i\n",Opsi.N);
    while(OperatorPsiConserveBasis(H,Opsi,psi))                /* Newly allocate memory for psi in Value. Determinant and Next point to same positions as Opsi*/
    {/* Allocation Failed */
      printf("Allocation of new memory failed in LanczosGroundState for the function OperatorPsiConserveBasis\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",*epsilon,*epsilon*2.0);
      fflush(stdout);
      *epsilon*=2.0;
      if(Opsi.NBlok==1)AllocFailed();
      RealWaveFunctionTruncateBasis(&Opsi,*epsilon);
    }
    RealWaveFunctionMultiply(Opsi,shift);
    RealWaveFunctionAddSameBasis(*psi,Opsi);
    RealWaveFunctionNormalize(psi);
    energy1=2*zero;
    energy2=4*zero;
    while((energy2-energy1)>(zero*(1.0+fabs(energy1))))
    {
      cycles++;
      if((cycles&255)==255) {printf("Reached %i-th cycle with energy %22.15E\n",cycles,energy1); fflush(stdout);}
      energy2=energy1;
      OperatorPsiConserveBasisPreAllocated(H,*psi,&Opsi);
      a0=RealWaveFunctionDotProductSameBasis(*psi,Opsi);
      RealWaveFunctionAddMultipliedSameBasis(Opsi,*psi,-a0);
      b1=RealWaveFunctionCalculateNorm(Opsi);
      if(b1==0)
      {
        energy1=a0;
      }
      else
      {
        RealWaveFunctionMultiply(Opsi,(double)1.0/b1);
        PsiOperatorPsiRR(H,Opsi,&a1);
        energy1=0.5*(a0+a1-sqrt(a0*a0-2*a0*a1+a1*a1+4*b1*b1));
        RealWaveFunctionMultiply(*psi,(energy1-a1)*2);
        RealWaveFunctionAddMultipliedSameBasis(*psi,Opsi,2*b1);
        RealWaveFunctionOrthogonalizeConserveBasis(psi, OtherPsi, NOtherPsi);
        RealWaveFunctionNormalize(psi);
      }
    }
    RealWaveFunctionFreeValues(&Opsi);                              /* Free  Opsi.ValueRe, but don't free  Opsi.Determinant as this is the same space as where psi.Determinant points to */
    printf("Energy: %22.15E  ",energy1);
    printf("NumberOfLaczosCycles: %i\n",cycles);
  }
  printf("\n");
  fflush(stdout);
}

uint32_t LanczosGroundStateRRTripple(const OperatorType H, WaveFunctionType *const psi, double *epsilon, const double zero, double shiftin, WaveFunctionType *OtherPsi, uint32_t NOtherPsi)
{
  printf("Start of LanczosGroundStateRR\n\n");
  fflush(stdout);
  WaveFunctionType Opsi,OOpsi;
  double energy=4*zero,energy1=2*zero,energy2;
  double shift,a0,a1,a2,b1,b2;
  double EnT[3],FunT[3][3];
  uint32_t cycles;
  
  while(RealWaveFunctionOrthogonalize(psi, OtherPsi, NOtherPsi))
  {/* Allocation Failed */
    printf("Allocation of new memory failed in LanczosGroundState for the function RealWaveFunctionOrthogonalize\n");
    printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",*epsilon,*epsilon*2.0);
    fflush(stdout);
    *epsilon*=2.0;
    if(Opsi.NBlok==1)AllocFailed();
    RealWaveFunctionTruncateBasis(psi,*epsilon);
  }
  RealWaveFunctionNormalize(psi);
  RealWaveFunctionTruncateBasis(psi,*epsilon);
  
  while((energy-energy1)>(zero*(1.0+fabs(energy))))
  {
    cycles=0;
    energy=energy1;
    shift=-shiftin-energy*0.49;
    while(OperatorPsiRR(H,*psi,&Opsi))                           /* Newly allocate memory for Opsi both in Value, Determinant and Next */
    {/* Allocation Failed */
      printf("Allocation of new memory failed in LanczosGroundState for the function OperatorPsi\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",*epsilon,*epsilon*2.0);
      fflush(stdout);
      *epsilon*=2.0;
      if(psi->NBlok==1)AllocFailed();
      RealWaveFunctionTruncateBasis(psi,*epsilon);
     }
    printf("BasisSetSizeAfter Expansion/Truncation : %12i ",Opsi.N);
    RealWaveFunctionMultiply(*psi,shift);
    while(RealWaveFunctionAddUnionBasis(&Opsi,*psi))               /* Possibly expand allocated memory of Opsi */
    {/* Allocation Failed */
      printf("Allocation of new memory failed in LanczosGroundState for the function RealWaveFunctionAddUnionBasis\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",*epsilon,*epsilon*2.0);
      fflush(stdout);
      *epsilon*=2.0;
      RealWaveFunctionFreeValues(&Opsi);  
      RealWaveFunctionMultiply(*psi,1.0/shift);
      if(psi->NBlok==1)AllocFailed();
      RealWaveFunctionTruncateBasis(psi,*epsilon);
      OperatorPsiRR(H,*psi,&Opsi);
      RealWaveFunctionMultiply(*psi,shift);
    }
    RealWaveFunctionFree(psi);                                     /* Remove allocated memory for psi from memory */
    while(RealWaveFunctionOrthogonalize(&Opsi, OtherPsi, NOtherPsi))
    {/* Allocation Failed */
      printf("Allocation of new memory failed in LanczosGroundState for the function RealWaveFunctionOrthogonalize\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",*epsilon,*epsilon*2.0);
      fflush(stdout);
      *epsilon*=2.0;
      if(Opsi.NBlok==1)AllocFailed();
      RealWaveFunctionTruncateBasis(&Opsi,*epsilon);
    }
    RealWaveFunctionNormalize(&Opsi);
    RealWaveFunctionTruncateBasis(&Opsi,*epsilon);                 /* Possibly remove allocated memory for Opsi */
    printf("%12i\n",Opsi.N);
    while(OperatorPsiConserveBasisRR(H,Opsi,&OOpsi))              /* Newly allocate memory for OOpsi in Value. Determinant and Next point to same positions as Opsi*/
    {/* Allocation Failed */
      printf("Allocation of new memory failed in LanczosGroundState for the function OperatorPsiConserveBasis\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",*epsilon,*epsilon*2.0);
      fflush(stdout);
      *epsilon*=2.0;
      if(Opsi.NBlok==1)AllocFailed();
      RealWaveFunctionTruncateBasis(&Opsi,*epsilon);
    }
    RealWaveFunctionMultiply(Opsi,shift);
    RealWaveFunctionAddSameBasis(OOpsi,Opsi);
    RealWaveFunctionNormalize(&OOpsi);
    while(OperatorPsiConserveBasisRR(H,OOpsi,psi))              /* Newly allocate memory for psi in Value. Determinant and Next point to same positions as Opsi*/
    {/* Allocation Failed */
      printf("Allocation of new memory failed in LanczosGroundState for the function OperatorPsiConserveBasis\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",*epsilon,*epsilon*2.0);
      fflush(stdout);
      *epsilon*=2.0;
      if(Opsi.NBlok==1)AllocFailed();
      RealWaveFunctionTruncateBasis(&Opsi,*epsilon);
      RealWaveFunctionFreeValues(&OOpsi);
      OperatorPsiConserveBasisRR(H,Opsi,&OOpsi);
    }
    RealWaveFunctionMultiply(OOpsi,shift);
    RealWaveFunctionAddSameBasis(*psi,OOpsi);
    RealWaveFunctionNormalize(psi);  
    energy1=2*zero;
    energy2=4*zero;
    while((cycles<4)||((energy2-energy1)>(zero*(1.0+fabs(energy1)))))
    {
      cycles++;
      if((cycles&255)==255)
      {
        printf("Reached %i-th cycle with energy %22.15E\n",cycles,energy1); fflush(stdout);
      }
      energy2=energy1;
      OperatorPsiConserveBasisPreAllocatedRR(H,*psi,&Opsi);
      a0=RealWaveFunctionDotProductSameBasis(*psi,Opsi);
      RealWaveFunctionAddMultipliedSameBasis(Opsi,*psi,-a0);
      b1=RealWaveFunctionCalculateNorm(Opsi);
      if(b1==0)
      {
        a1=fabs(a0)+2.0;
        a2=fabs(a0)+fabs(a1)+2.0;
        b2=0;
      }
      else
      {
        RealWaveFunctionMultiply(Opsi,(double)1.0/b1);
        OperatorPsiConserveBasisPreAllocatedRR(H,Opsi,&OOpsi);
        a1=RealWaveFunctionDotProductSameBasis(Opsi,OOpsi);
        RealWaveFunctionAddMultipliedSameBasis(OOpsi,Opsi,-a1);
        b2=RealWaveFunctionCalculateNorm(OOpsi);
        if(b2==0)
        {
          a2=fabs(a0)+fabs(a1)+2.0;
        }
        else
        {
          RealWaveFunctionMultiply(OOpsi,(double)1.0/b2);
          PsiOperatorPsiRR(H,OOpsi,&a2);
        }
      }
      Eigensystem3(a0,a1,a2,b1,b2,EnT,FunT);
      energy1=EnT[0];
      RealWaveFunctionMultiply(*psi,FunT[0][0]);
      RealWaveFunctionAddMultipliedSameBasis(*psi,Opsi,FunT[0][1]);
      RealWaveFunctionAddMultipliedSameBasis(*psi,OOpsi,FunT[0][2]);
      RealWaveFunctionOrthogonalizeConserveBasis(psi, OtherPsi, NOtherPsi);
      RealWaveFunctionNormalize(psi);
    }
    RealWaveFunctionFreeValues(&Opsi);                              /* Free  Opsi.ValueRe, but don't free  Opsi.Determinant as this is the same space as where psi.Determinant points to */
    RealWaveFunctionFreeValues(&OOpsi);                             /* Free OOpsi.ValueRe, but don't free OOpsi.Determinant as this is the same space as where psi.Determinant points to */
    printf("Energy: %22.15E  ",energy1);
    printf("NumberOfLaczosCycles: %i\n",cycles);
  }
  printf("\n");
  fflush(stdout);
  return 0;
}


uint32_t LanczosGroundStateRCTripple(const OperatorType H, WaveFunctionType *const psi, double *epsilon, const double zero, double shiftin, WaveFunctionType *OtherPsi, uint32_t NOtherPsi)
{
  printf("Start of LanczosGroundStateRC\n\n");
  fflush(stdout);
  WaveFunctionType Opsi,OOpsi;
  double energy=4*zero,energy1=2*zero,energy2;
  double shift,a0,a1,a2,b1,b2,a0Im,a1Im;
  double EnT[3],FunT[3][3];
  uint32_t cycles;
  while((energy-energy1)>(zero*(1.0+fabs(energy))))
  {
    cycles=0;
    energy=energy1;
    shift=-shiftin-energy*0.49;
    while(OperatorPsiRC(H,*psi,&Opsi))                           /* Newly allocate memory for Opsi both in Value, Determinant and Next */
    {/* Allocation Failed */
      printf("Allocation of new memory failed in LanczosGroundState for the function OperatorPsi\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",*epsilon,*epsilon*2.0);
      fflush(stdout);
      *epsilon*=2.0;
      if(psi->NBlok==1)AllocFailed();
      ComplexWaveFunctionTruncateBasis(psi,*epsilon);
     }
    printf("BasisSetSizeAfter Expansion/Truncation : %12i ",Opsi.N);
    ComplexWaveFunctionMultiply(*psi,shift);
    while(ComplexWaveFunctionAddUnionBasis(&Opsi,*psi))               /* Possibly expand allocated memory of Opsi */
    {/* Allocation Failed */
      printf("Allocation of new memory failed in LanczosGroundState for the function RealWaveFunctionAddUnionBasis\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",*epsilon,*epsilon*2.0);
      fflush(stdout);
      *epsilon*=2.0;
      ComplexWaveFunctionFreeValues(&Opsi);  
      ComplexWaveFunctionMultiply(*psi,1.0/shift);
      if(psi->NBlok==1)AllocFailed();
      ComplexWaveFunctionTruncateBasis(psi,*epsilon);
      OperatorPsiRC(H,*psi,&Opsi);
      ComplexWaveFunctionMultiply(*psi,shift);
    }
    ComplexWaveFunctionFree(psi);                                     /* Remove allocated memory for psi from memory */
    while(ComplexWaveFunctionOrthogonalize(&Opsi, OtherPsi, NOtherPsi))
    {/* Allocation Failed */
      printf("Allocation of new memory failed in LanczosGroundState for the function RealWaveFunctionOrthogonalize\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",*epsilon,*epsilon*2.0);
      fflush(stdout);
      *epsilon*=2.0;
      if(Opsi.NBlok==1)AllocFailed();
      ComplexWaveFunctionTruncateBasis(&Opsi,*epsilon);
    }
    ComplexWaveFunctionNormalize(&Opsi);
    ComplexWaveFunctionTruncateBasis(&Opsi,*epsilon);                 /* Possibly remove allocated memory for Opsi */
    printf("%12i\n",Opsi.N);
    while(OperatorPsiConserveBasisRC(H,Opsi,&OOpsi))              /* Newly allocate memory for OOpsi in Value. Determinant and Next point to same positions as Opsi*/
    {/* Allocation Failed */
      printf("Allocation of new memory failed in LanczosGroundState for the function OperatorPsiConserveBasis\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",*epsilon,*epsilon*2.0);
      fflush(stdout);
      *epsilon*=2.0;
      if(Opsi.NBlok==1)AllocFailed();
      ComplexWaveFunctionTruncateBasis(&Opsi,*epsilon);
    }
    ComplexWaveFunctionMultiply(Opsi,shift);
    ComplexWaveFunctionAddSameBasis(OOpsi,Opsi);
    ComplexWaveFunctionNormalize(&OOpsi);
    while(OperatorPsiConserveBasisRC(H,OOpsi,psi))              /* Newly allocate memory for psi in Value. Determinant and Next point to same positions as Opsi*/
    {/* Allocation Failed */
      printf("Allocation of new memory failed in LanczosGroundState for the function OperatorPsiConserveBasis\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",*epsilon,*epsilon*2.0);
      fflush(stdout);
      *epsilon*=2.0;
      if(Opsi.NBlok==1)AllocFailed();
      ComplexWaveFunctionTruncateBasis(&Opsi,*epsilon);
      ComplexWaveFunctionFreeValues(&OOpsi);
      OperatorPsiConserveBasisRC(H,Opsi,&OOpsi);
    }
    ComplexWaveFunctionMultiply(OOpsi,shift);
    ComplexWaveFunctionAddSameBasis(*psi,OOpsi);
    ComplexWaveFunctionNormalize(psi);  
    energy1=2*zero;
    energy2=4*zero;
    while((cycles<4)||((energy2-energy1)>(zero*(1.0+fabs(energy1)))))
    {
      cycles++;
      if((cycles&255)==255) {printf("Reached %i-th cycle with energy %22.15E\n",cycles,energy1); fflush(stdout);}
      energy2=energy1;
      OperatorPsiConserveBasisPreAllocatedRC(H,*psi,&Opsi);
      ComplexWaveFunctionDotProductSameBasis(*psi,Opsi,&a0,&a0Im);
      ComplexWaveFunctionAddMultipliedSameBasis(Opsi,*psi,-a0);
      b1=ComplexWaveFunctionCalculateNorm(Opsi);
      if(b1==0)
      {
        a1=fabs(a0)+2.0;
        a2=fabs(a0)+fabs(a1)+2.0;
        b2=0;
      }
      else
      {
        ComplexWaveFunctionMultiply(Opsi,(double)1.0/b1);
        OperatorPsiConserveBasisPreAllocatedRC(H,Opsi,&OOpsi);
        ComplexWaveFunctionDotProductSameBasis(Opsi,OOpsi,&a1,&a1Im);
        ComplexWaveFunctionAddMultipliedSameBasis(OOpsi,Opsi,-a1);
        b2=ComplexWaveFunctionCalculateNorm(OOpsi);
        if(b2==0)
        {
          a2=fabs(a0)+fabs(a1)+2.0;
        }
        else
        {
          ComplexWaveFunctionMultiply(OOpsi,(double)1.0/b2);
          PsiOperatorPsiRHC(H,OOpsi,&a2);
        }
      }
      Eigensystem3(a0,a1,a2,b1,b2,EnT,FunT);
      energy1=EnT[0];
      ComplexWaveFunctionMultiply(*psi,FunT[0][0]);
      ComplexWaveFunctionAddMultipliedSameBasis(*psi,Opsi,FunT[0][1]);
      ComplexWaveFunctionAddMultipliedSameBasis(*psi,OOpsi,FunT[0][2]);
      ComplexWaveFunctionOrthogonalizeConserveBasis(psi, OtherPsi, NOtherPsi);
      ComplexWaveFunctionNormalize(psi);
    }
    ComplexWaveFunctionFreeValues(&Opsi);                              /* Free  Opsi.ValueRe, but don't free  Opsi.Determinant as this is the same space as where psi.Determinant points to */
    ComplexWaveFunctionFreeValues(&OOpsi);                             /* Free OOpsi.ValueRe, but don't free OOpsi.Determinant as this is the same space as where psi.Determinant points to */
    printf("Energy: %22.15E  ",energy1);
    printf("NumberOfLaczosCycles: %i\n",cycles);
  }
  printf("\n");
  fflush(stdout);
  return 0;
}

uint32_t LanczosGroundStateCCTripple(const OperatorType H, WaveFunctionType *const psi, double *epsilon, const double zero, double shiftin, WaveFunctionType *OtherPsi, uint32_t NOtherPsi)
{
  printf("Start of LanczosGroundStateCC\n\n");
  fflush(stdout);
  WaveFunctionType Opsi,OOpsi;
  double energy=4*zero,energy1=2*zero,energy2;
  double shift,a0,a1,a2,b1,b2,a0Im,a1Im;
  double EnT[3],FunT[3][3];
  uint32_t cycles;
  while((energy-energy1)>(zero*(1.0+fabs(energy))))
  {
    cycles=0;
    energy=energy1;
    shift=-shiftin-energy*0.49;
    while(OperatorPsiCC(H,*psi,&Opsi))                           /* Newly allocate memory for Opsi both in Value, Determinant and Next */
    {/* Allocation Failed */
      printf("Allocation of new memory failed in LanczosGroundState for the function OperatorPsi\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",*epsilon,*epsilon*2.0);
      fflush(stdout);
      *epsilon*=2.0;
      if(psi->NBlok==1)AllocFailed();
      ComplexWaveFunctionTruncateBasis(psi,*epsilon);
     }
    printf("BasisSetSizeAfter Expansion/Truncation : %12i ",Opsi.N);
    ComplexWaveFunctionMultiply(*psi,shift);
    while(ComplexWaveFunctionAddUnionBasis(&Opsi,*psi))               /* Possibly expand allocated memory of Opsi */
    {/* Allocation Failed */
      printf("Allocation of new memory failed in LanczosGroundState for the function RealWaveFunctionAddUnionBasis\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",*epsilon,*epsilon*2.0);
      fflush(stdout);
      *epsilon*=2.0;
      ComplexWaveFunctionFreeValues(&Opsi);  
      ComplexWaveFunctionMultiply(*psi,1.0/shift);
      if(psi->NBlok==1)AllocFailed();
      ComplexWaveFunctionTruncateBasis(psi,*epsilon);
      OperatorPsiCC(H,*psi,&Opsi);
      ComplexWaveFunctionMultiply(*psi,shift);
    }
    ComplexWaveFunctionFree(psi);                                     /* Remove allocated memory for psi from memory */
    while(ComplexWaveFunctionOrthogonalize(&Opsi, OtherPsi, NOtherPsi))
    {/* Allocation Failed */
      printf("Allocation of new memory failed in LanczosGroundState for the function RealWaveFunctionOrthogonalize\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",*epsilon,*epsilon*2.0);
      fflush(stdout);
      *epsilon*=2.0;
      if(Opsi.NBlok==1)AllocFailed();
      ComplexWaveFunctionTruncateBasis(&Opsi,*epsilon);
    }
    ComplexWaveFunctionNormalize(&Opsi);
    ComplexWaveFunctionTruncateBasis(&Opsi,*epsilon);                 /* Possibly remove allocated memory for Opsi */
    printf("%12i\n",Opsi.N);
    while(OperatorPsiConserveBasisCC(H,Opsi,&OOpsi))              /* Newly allocate memory for OOpsi in Value. Determinant and Next point to same positions as Opsi*/
    {/* Allocation Failed */
      printf("Allocation of new memory failed in LanczosGroundState for the function OperatorPsiConserveBasis\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",*epsilon,*epsilon*2.0);
      fflush(stdout);
      *epsilon*=2.0;
      if(Opsi.NBlok==1)AllocFailed();
      ComplexWaveFunctionTruncateBasis(&Opsi,*epsilon);
    }
    ComplexWaveFunctionMultiply(Opsi,shift);
    ComplexWaveFunctionAddSameBasis(OOpsi,Opsi);
    ComplexWaveFunctionNormalize(&OOpsi);
    while(OperatorPsiConserveBasisCC(H,OOpsi,psi))              /* Newly allocate memory for psi in Value. Determinant and Next point to same positions as Opsi*/
    {/* Allocation Failed */
      printf("Allocation of new memory failed in LanczosGroundState for the function OperatorPsiConserveBasis\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",*epsilon,*epsilon*2.0);
      fflush(stdout);
      *epsilon*=2.0;
      if(Opsi.NBlok==1)AllocFailed();
      ComplexWaveFunctionTruncateBasis(&Opsi,*epsilon);
      ComplexWaveFunctionFreeValues(&OOpsi);
      OperatorPsiConserveBasisCC(H,Opsi,&OOpsi);
    }
    ComplexWaveFunctionMultiply(OOpsi,shift);
    ComplexWaveFunctionAddSameBasis(*psi,OOpsi);
    ComplexWaveFunctionNormalize(psi);  
    energy1=2*zero;
    energy2=4*zero;
    while((cycles<4)||((energy2-energy1)>(zero*(1.0+fabs(energy1)))))
    {
      cycles++;
      if((cycles&255)==255) {printf("Reached %i-th cycle with energy %22.15E\n",cycles,energy1); fflush(stdout);}
      energy2=energy1;
      OperatorPsiConserveBasisPreAllocatedCC(H,*psi,&Opsi);
      ComplexWaveFunctionDotProductSameBasis(*psi,Opsi,&a0,&a0Im);
      ComplexWaveFunctionAddMultipliedSameBasis(Opsi,*psi,-a0);
      b1=ComplexWaveFunctionCalculateNorm(Opsi);
      if(b1==0)
      {
        a1=fabs(a0)+2.0;
        a2=fabs(a0)+fabs(a1)+2.0;
        b2=0;
      }
      else
      {
        ComplexWaveFunctionMultiply(Opsi,(double)1.0/b1);
        OperatorPsiConserveBasisPreAllocatedCC(H,Opsi,&OOpsi);
        ComplexWaveFunctionDotProductSameBasis(Opsi,OOpsi,&a1,&a1Im);
        ComplexWaveFunctionAddMultipliedSameBasis(OOpsi,Opsi,-a1);
        b2=ComplexWaveFunctionCalculateNorm(OOpsi);
        if(b2==0)
        {
          a2=fabs(a0)+fabs(a1)+2.0;
        }
        else
        {
          ComplexWaveFunctionMultiply(OOpsi,(double)1.0/b2);
          PsiOperatorPsiCHC(H,OOpsi,&a2);
        }
      }
      Eigensystem3(a0,a1,a2,b1,b2,EnT,FunT);
      energy1=EnT[0];
      ComplexWaveFunctionMultiply(*psi,FunT[0][0]);
      ComplexWaveFunctionAddMultipliedSameBasis(*psi,Opsi,FunT[0][1]);
      ComplexWaveFunctionAddMultipliedSameBasis(*psi,OOpsi,FunT[0][2]);
      ComplexWaveFunctionOrthogonalizeConserveBasis(psi, OtherPsi, NOtherPsi);
      ComplexWaveFunctionNormalize(psi);
    }
    ComplexWaveFunctionFreeValues(&Opsi);                              /* Free  Opsi.ValueRe, but don't free  Opsi.Determinant as this is the same space as where psi.Determinant points to */
    ComplexWaveFunctionFreeValues(&OOpsi);                             /* Free OOpsi.ValueRe, but don't free OOpsi.Determinant as this is the same space as where psi.Determinant points to */
    printf("Energy: %22.15E  ",energy1);
    printf("NumberOfLaczosCycles: %i\n",cycles);
  }
  printf("\n");
  fflush(stdout);
  return 0;
}


void LanczosTriDiagonalizeConserveBasisPreAllocated(const OperatorType H, double shift, WaveFunctionType v[], TriDiagonalMatrixType *const M)
{
  uint32_t i,vi,vip,vim;
  M->b[0]=1.0;
  OperatorPsiConserveBasisPreAllocated(H,v[0],&(v[1]));
  RealWaveFunctionAddMultipliedConserveBasis(&(v[1]),v[0],shift);
  M->a[0]=RealWaveFunctionDotProduct(v[1],v[0]);
  for(i=1;i<M->N;i++)
  {
    vip=(i+1)%3;
    vi =(i  )%3;
    vim=(i-1)%3;
    RealWaveFunctionAddMultipliedConserveBasis(&(v[vi]),v[vim],-M->a[i-1]);
    M->b[i]=RealWaveFunctionCalculateNorm(v[vi]);
    RealWaveFunctionMultiply(v[vi],(double)1.0/M->b[i]);
    OperatorPsiConserveBasisPreAllocated(H,v[vi],&(v[vip]));
    RealWaveFunctionAddMultipliedConserveBasis(&(v[vip]),v[vim],shift);
    M->a[i]=RealWaveFunctionDotProduct(v[vip],v[vi]);
    RealWaveFunctionAddMultipliedConserveBasis(&(v[vip]),v[vim],-M->b[i]);
  }
}

