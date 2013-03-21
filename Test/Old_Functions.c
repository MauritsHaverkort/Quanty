









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


/* =========================================== old DMFT functions ============================================= */



int oldmainDyson(DMFTSettings HSet, TimeType *Time)
{
  WaveFunctionType Psi,OPsi,PsiTemp;
  TriDiagonalMatrixType G0,GBath,GclPup,GclPdn,GclMup,GclMdn,Gclup,Gcldn,Gcl,Sigma,GBathNew,GBathOld,gPerturbation[4];
  DoubleTriDiagonalMatrixType SigmaDouble;
  SecondDoubleTriDiagonalMatrixType GAnd,GAndTrunk;
  OperatorType HAnd, OPESdn,OPESup,OIPESdn,OIPESup,VPerturbation[4];
  MatrixType T;
  SpectraType Spectra[10];
  char GBathFileName[255],PsiFileName[255],SpecFileName[25],Buffer[255];
  uint32_t iOuter,i,PsiNBathV,PsiNBathC;
  uint16_t oadd[1];
  double Ef=0.0,epsilonAnd=0.00001,mu;
  double PsiEnergy,PsiQuality,PsiHist[20],*PsiRho;
  /* read G0 */
  if(FileReadPlainTriDiagonalMatrixType(&G0,HSet.G0FileName)){printf("Fatal Error in main\n Failed to read G0 from file %s\n",HSet.G0FileName);return 1;};
  strcpy(G0.Name,"G0");
  /* read G_Bath_iStart and psi_iStart*/
  if(HSet.iStart==0)
  {
    sprintf(GBathFileName,"%s000",HSet.GRootFileName);
    sprintf(PsiFileName,"%s000",HSet.WaveFunctionRootFileName);
  }
  else if(HSet.iStart<10)
  {
    sprintf(GBathFileName,"%s00%1i",HSet.GRootFileName,HSet.iStart);
    sprintf(PsiFileName,"%s00%1i",HSet.WaveFunctionRootFileName,HSet.iStart);
  }
  else if(HSet.iStart<100)
  {
    sprintf(GBathFileName,"%s0%2i",HSet.GRootFileName,HSet.iStart);
    sprintf(PsiFileName,"%s0%2i",HSet.WaveFunctionRootFileName,HSet.iStart);
  }
  else
  {
    sprintf(GBathFileName,"%s%i",HSet.GRootFileName,HSet.iStart);
    sprintf(PsiFileName,"%s%i",HSet.WaveFunctionRootFileName,HSet.iStart);
  }
  if(FileReadPlainTriDiagonalMatrixType(&GBath,GBathFileName)){printf("Fatal Error in main\n Failed to read GBath from file %s !\n",GBathFileName);return 1;};
  if(file_exists(PsiFileName))
  {
    Psi.NBitsKey=HSet.NBitsKey;
    if(FileReadWaveFunctionPlainBinary(&Psi, PsiFileName)){printf("Fatal Error in main\n Failed to read Psi from file %s !\n",PsiFileName);return 1;};
    printf("\nRead the starting wave function from the file %s\n",PsiFileName);
    sscanf(Psi.Name,"%s %i %i",Buffer,&PsiNBathC,&PsiNBathV);
    printf("Wave-function \"%s\" found has %i conduction bath sites, %i valence bath sites and %i Fermions in basis\n",Psi.Name,PsiNBathC,PsiNBathV,Psi.NFermion);
    if((2*(PsiNBathC+PsiNBathV+2))!=Psi.NFermion)
    {
      printf("The total number of Fermions is incompatible witht the bath sizes\n");
      printf("This wave-function can not be used to start from\n");
      return 1;
    }
  }
  else
  {
    printf("\nWave-function file: '%s' does not excist\n Started from new function\n",PsiFileName);
    CreateWaveFunctionDMFTSecond(&Psi, HSet.NBath, HSet.NBath, HSet.NBitsKey);
  }
  strcpy(GBath.Name,"Gbath");
  /* do the DMFT loops */
  for(iOuter=HSet.iStart;iOuter<HSet.iStart+HSet.Ni;iOuter++)
  {
    BeginTiming("DMFTLoop",Time);
    printf("\n\n***********************************************************************\n");
    printf("Start of loop %i from %i to %i\n",iOuter+1,HSet.iStart,HSet.iStart+HSet.Ni);
    printf("***********************************************************************\n\n\n");
    /* Create splitted anderson matrix */
    TriDiagonalMatrixToSecondDoubleTriDiagonalMatrix(GBath, &T, &GAnd, Ef, epsilonAnd);
    /* Create GAnd and GAndTrunk as a copy pointing to the same memory */
    GAndTrunk.Ef=GAnd.Ef;
    GAndTrunk.Epsilon=GAnd.Epsilon;
    GAndTrunk.a0=GAnd.a0;
    GAndTrunk.b0=GAnd.b0;
    GAndTrunk.Dt.Ef=GAnd.Dt.Ef;
    GAndTrunk.Dt.Epsilon=GAnd.Dt.Epsilon;
    GAndTrunk.Dt.a0=GAnd.Dt.a0;
    GAndTrunk.Dt.b0=GAnd.Dt.b0;
    GAndTrunk.Dt.val.N=GAnd.Dt.val.N;
    GAndTrunk.Dt.val.a=GAnd.Dt.val.a;
    GAndTrunk.Dt.val.b=GAnd.Dt.val.b;
    GAndTrunk.Dt.con.N=GAnd.Dt.con.N;
    GAndTrunk.Dt.con.a=GAnd.Dt.con.a;
    GAndTrunk.Dt.con.b=GAnd.Dt.con.b;
    /* Create VPerturbation */
    for(i=0;i<4;i++)
    {
      InitOperator(&(VPerturbation[i]));
      VPerturbation[i].NFermion=4+2*(Min(GAnd.Dt.val.N,HSet.NBath)+Min(GAnd.Dt.con.N,HSet.NBath));
      VPerturbation[i].NBoson=0;
      VPerturbation[i].NWord=Min(VPerturbation[i].NFermion,1+(VPerturbation[i].NFermion-1)/8);
      VPerturbation[i].QComplex=0;
      InitOperatorFixedLength(&(VPerturbation[i]),1,1);
    }
    if(GAnd.Dt.val.N>HSet.NBath)
    {
      oadd[0]=CrAnMask+3+HSet.NBath;
      RealOperatorAddLader(&(VPerturbation[0]), oadd, 1, GAnd.Dt.val.b[HSet.NBath]);
      oadd[0]=CrAnMask+3+Min(GAnd.Dt.con.N,HSet.NBath)+2*HSet.NBath;
      RealOperatorAddLader(&(VPerturbation[1]), oadd, 1, GAnd.Dt.val.b[HSet.NBath]);
    }
    if(GAnd.Dt.con.N>HSet.NBath)
    {
      oadd[0]=3+Min(GAnd.Dt.val.N,HSet.NBath)+HSet.NBath;
      RealOperatorAddLader(&(VPerturbation[2]), oadd, 1, GAnd.Dt.con.b[HSet.NBath]);
      oadd[0]=3+2*Min(GAnd.Dt.val.N,HSet.NBath)+2*HSet.NBath;
      RealOperatorAddLader(&(VPerturbation[3]), oadd, 1, GAnd.Dt.con.b[HSet.NBath]);
    }
    /* Create gPerturbation */
    if(GAnd.Dt.val.N>HSet.NBath)
    {
      gPerturbation[0].N=GAnd.Dt.val.N-HSet.NBath;
      gPerturbation[1].N=GAnd.Dt.val.N-HSet.NBath;
      InitTriDiagonalMatrix(&(gPerturbation[0]));
      InitTriDiagonalMatrix(&(gPerturbation[1]));
      for(i=1;i<GAnd.Dt.val.N-HSet.NBath;i++)
      {
        gPerturbation[0].a[i]=-GAnd.Dt.val.a[i+HSet.NBath];
        gPerturbation[0].b[i]=GAnd.Dt.val.b[i+HSet.NBath];
        gPerturbation[1].a[i]=-GAnd.Dt.val.a[i+HSet.NBath];
        gPerturbation[1].b[i]=GAnd.Dt.val.b[i+HSet.NBath];
      }
      gPerturbation[0].a[GAnd.Dt.val.N-HSet.NBath]=-GAnd.Dt.val.a[GAnd.Dt.val.N];
      gPerturbation[1].a[GAnd.Dt.val.N-HSet.NBath]=-GAnd.Dt.val.a[GAnd.Dt.val.N];
    }
    else
    {
      gPerturbation[0].N=1;
      gPerturbation[1].N=1;
      InitTriDiagonalMatrix(&(gPerturbation[0]));
      InitTriDiagonalMatrix(&(gPerturbation[1]));
      gPerturbation[0].a[1]=0;
      gPerturbation[1].a[1]=0;
    }
    if(GAnd.Dt.con.N>HSet.NBath)
    {
      gPerturbation[2].N=GAnd.Dt.con.N-HSet.NBath;
      gPerturbation[3].N=GAnd.Dt.con.N-HSet.NBath;
      InitTriDiagonalMatrix(&(gPerturbation[2]));
      InitTriDiagonalMatrix(&(gPerturbation[3]));
      for(i=1;i<GAnd.Dt.con.N-HSet.NBath;i++)
      {
        gPerturbation[2].a[i]=GAnd.Dt.con.a[i+HSet.NBath];
        gPerturbation[2].b[i]=GAnd.Dt.con.b[i+HSet.NBath];
        gPerturbation[3].a[i]=GAnd.Dt.con.a[i+HSet.NBath];
        gPerturbation[3].b[i]=GAnd.Dt.con.b[i+HSet.NBath];
      }
      gPerturbation[2].a[GAnd.Dt.con.N-HSet.NBath]=GAnd.Dt.con.a[GAnd.Dt.con.N];
      gPerturbation[3].a[GAnd.Dt.con.N-HSet.NBath]=GAnd.Dt.con.a[GAnd.Dt.con.N];
    }
    else
    {
      gPerturbation[2].N=1;
      gPerturbation[3].N=1;
      InitTriDiagonalMatrix(&(gPerturbation[2]));
      InitTriDiagonalMatrix(&(gPerturbation[3]));
      gPerturbation[2].a[1]=0;
      gPerturbation[3].a[1]=0;
    };
    /* Create HAnd */
    FreeMatrix(&T);
    if(GAnd.Dt.con.N<HSet.NBath)
    {
      printf("Warning: GAnd.Dt.con.N(%i)<HSet.NBath(%i)\n",GAnd.Dt.con.N,HSet.NBath);
    }
    else
    {
      GAndTrunk.Dt.con.N=HSet.NBath;
    }
    if(GAnd.Dt.val.N<HSet.NBath)
    {
      printf("Warning: GAnd.Dt.val.N(%i)<HSet.NBath(%i)\n",GAnd.Dt.val.N,HSet.NBath);
    }
    else
    {
      GAndTrunk.Dt.val.N=HSet.NBath;
    }
    mu=-HSet.U/2;
    SecondDoubleTriDiagonalMatrixToAndersonImpurityOperator(&HAnd, GAndTrunk, GAndTrunk, HSet.U, mu);
    /* test if psi has correct dimensions and correct if needed */
    sscanf(Psi.Name,"%s %i %i",Buffer,&PsiNBathC,&PsiNBathV);
    if(PsiNBathC!=GAndTrunk.Dt.con.N)
    {
      printf("\nBefore diagonalizing:\nThe wavefunction has %i conduction bath sites\nThe Anderson Hamiltonian %i\n",PsiNBathC,GAndTrunk.Dt.con.N);
      if(PsiNBathC>GAndTrunk.Dt.con.N)
      {
        printf("Removing %i unoccupied fermions from the wave-function\n",PsiNBathC-GAndTrunk.Dt.con.N);
        for(i=0;i<PsiNBathC-GAndTrunk.Dt.con.N;i++)
        {
          WaveFunctionRemoveFermion(Psi, &PsiTemp, Psi.NFermion-1);
          WaveFunctionFree(&Psi);
          WaveFunctionRemoveFermion(PsiTemp, &Psi, 4+PsiNBathV+PsiNBathC-i-1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      else
      {
        printf("Adding %i unoccupied fermions to the wave-function\n",GAndTrunk.Dt.con.N-PsiNBathC);
        for(i=0;i<GAndTrunk.Dt.con.N-PsiNBathC;i++)
        {
          WaveFunctionAddFermion(Psi, &PsiTemp, Psi.NFermion,0);
          WaveFunctionFree(&Psi);
          WaveFunctionAddFermion(PsiTemp, &Psi, 4+PsiNBathV+PsiNBathC+i,0);
          WaveFunctionFree(&PsiTemp);
        }
      }
    }
    if(PsiNBathV!=GAndTrunk.Dt.val.N)
    {
      printf("\nBefore diagonalizing:\nThe wavefunction has %i valence bath sites\nThe Anderson Hamiltonian %i\n",PsiNBathV,GAndTrunk.Dt.val.N);
      if(PsiNBathV>GAndTrunk.Dt.val.N)
      {
        printf("Removing %i occupied fermions from the wave-function\n",PsiNBathV-GAndTrunk.Dt.val.N);
        for(i=0;i<PsiNBathV-GAndTrunk.Dt.val.N;i++)
        {
          WaveFunctionRemoveFermion(Psi, &PsiTemp, Psi.NFermion-GAndTrunk.Dt.con.N-1);
          WaveFunctionFree(&Psi);
          WaveFunctionRemoveFermion(PsiTemp, &Psi, 4+PsiNBathV-i-1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      else
      {
        printf("Adding %i occupied fermions to the wave-function\n",GAndTrunk.Dt.val.N-PsiNBathV);
        for(i=0;i<GAndTrunk.Dt.val.N-PsiNBathV;i++)
        {
          WaveFunctionAddFermion(Psi, &PsiTemp, Psi.NFermion-GAndTrunk.Dt.val.N,1);
          WaveFunctionFree(&Psi);
          WaveFunctionAddFermion(PsiTemp, &Psi, 4+PsiNBathV+i,1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      printf("Corrected basis\n");
      PrintWaveFunctionReadableTruncated(Psi, 0.1);
    }
    sprintf(Psi.Name,"DMFT %i %i",GAndTrunk.Dt.con.N,GAndTrunk.Dt.val.N);
    /* create ground-state wave-function */
    BeginTiming("Lanczos GroundState",Time);
    LanczosGroundState(HAnd,&Psi,&HSet.zeroLancWavefunction,HSet.zeroLancEnergy,HSet.U,&Psi,0);
    LanczosGroundState(HAnd,&Psi,&HSet.zeroLancWavefunction,HSet.zeroLancEnergy,HSet.U,&Psi,0);
    EndTiming("Lanczos GroundState",Time);
    PsiHermitianOperatorPsi(HAnd, Psi, &PsiEnergy);
    PsiHermitianOperatorSqrPsi(HAnd,Psi,&PsiQuality);
    PsiQuality/=Sqr(PsiEnergy);
    PrintWaveFunctionReadableTruncated(Psi, 0.1);
    CalculatePrefactorHistogram(PsiHist, Psi, 10.0, 20);
    printf("\nHistogram of prefactors in Wave-function\n");
    for(i=0;i<19;i++)
    {
      printf("10^(-%2i) to 10^(-%2i) | %22.15E\n",i,i+1,PsiHist[i]);
    }
    printf("10^(-%2i) to 0        | %22.15E\n",i,PsiHist[19]);
    printf("\n");
    printf("Energy                  <psi|H^2|psi>/(<psi|H|psi>^2)\n");
    printf("%22.15E  %22.15E\n",PsiEnergy,PsiQuality);
    PsiRho=calloc(Psi.NFermion,sizeof(double));
    if(PsiRho==NULL){printf("calloc failed in main for PsiRho\n Fatal Error\n");fflush(stdout);return 1;}
    CalculateDensity(PsiRho,Psi);
    printf("Occupation \n");
    for(i=0;i<Psi.NFermion;i++)
    {
      printf("%5i %22.15E\n",i,PsiRho[i]);
    }
    free(PsiRho);
    /*                                           */
    /* create spin up and spin down PES and IPES */
    /*                                           */
    CreatePESdnOperatorDMFT(&OPESdn, HAnd.NFermion);
    CreatePESupOperatorDMFT(&OPESup, HAnd.NFermion);
    CreateIPESdnOperatorDMFT(&OIPESdn, HAnd.NFermion);
    CreateIPESupOperatorDMFT(&OIPESup, HAnd.NFermion);
    GclMup.N=Max(HSet.NTriDiag,G0.N);
    InitTriDiagonalMatrix(&GclMup);
    GclMdn.N=Max(HSet.NTriDiag,G0.N);
    InitTriDiagonalMatrix(&GclMdn);
    GclPup.N=Max(HSet.NTriDiag,G0.N);
    InitTriDiagonalMatrix(&GclPup);
    GclPdn.N=Max(HSet.NTriDiag,G0.N);
    InitTriDiagonalMatrix(&GclPdn);
    BeginTiming("Lanczos Tri Diagonalization Dyson",Time);
    OperatorPsi(OPESup,Psi,&OPsi);
    LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation, gPerturbation, 4, &GclMup);
    EndTiming("Lanczos Tri Diagonalization Dyson",Time);
    BeginTiming("Lanczos Tri Diagonalization Dyson",Time);
    OperatorPsi(OPESdn,Psi,&OPsi);
    LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation, gPerturbation, 4, &GclMdn);
    EndTiming("Lanczos Tri Diagonalization Dyson",Time);
    BeginTiming("Lanczos Tri Diagonalization Dyson",Time);
    OperatorPsi(OIPESup,Psi,&OPsi);
    LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation, gPerturbation, 4, &GclPup);
    EndTiming("Lanczos Tri Diagonalization Dyson",Time);
    BeginTiming("Lanczos Tri Diagonalization Dyson",Time);
    OperatorPsi(OIPESdn,Psi,&OPsi);
    LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation, gPerturbation, 4, &GclPdn);
    EndTiming("Lanczos Tri Diagonalization Dyson",Time);
    /* G math - create Gcl as the sum of PES - IPES */
    printf("Start G math\n");
    BeginTiming("G math",Time);
    ShiftTriDiagonalMatrix(&GclMdn, -PsiEnergy);
    ShiftTriDiagonalMatrix(&GclMup, -PsiEnergy);
    ShiftTriDiagonalMatrix(&GclPdn, -PsiEnergy);
    ShiftTriDiagonalMatrix(&GclPup, -PsiEnergy);
    for(i=1;i<=GclMdn.N;i++)
    {
      GclMdn.a[i]*=-1;
    }
    for(i=1;i<=GclMup.N;i++)
    {
      GclMup.a[i]*=-1;
    }
    strcpy(GclPup.Name,"Gcluster-Plus-Up");
    strcpy(GclPdn.Name,"Gcluster-Plus-Dn");
    strcpy(GclMup.Name,"Gcluster-Min-Up");
    strcpy(GclMdn.Name,"Gcluster-Min-Dn");
    printf("  Start SumGreensfunctionsTriDiagonalRepresentation\n");fflush(stdout);
    SumGreensfunctionsTriDiagonalRepresentation(GclMup, GclPup, &Gclup);
    strcpy(Gclup.Name,"Gcluster-Up");
    SumGreensfunctionsTriDiagonalRepresentation(GclMdn, GclPdn, &Gcldn);
    strcpy(Gcldn.Name,"Gcluster-Dn");
    SumGreensfunctionsTriDiagonalRepresentation(Gcldn, Gclup, &Gcl);
    strcpy(Gcl.Name,"Gcluster");
    Gcl.b[0]=1.0;
    printf("Gcl.N = %i ",Gcl.N);
    SymmetrizeTridiagonalMatrix(&Gcl);
    printf("=> %i\n",Gcl.N);
    /* G math - create Self energy From Gcl and GBath */
    BeginTiming("G math Selfenergy",Time);
    printf("  Start SelfEnergyTriDiagonalRepresentation\n");
    printf("G0.N = %i ",G0.N);
    SelfEnergyTriDiagonalRepresentation(GBath, Gcl, &Sigma);
    printf("Sigma.N = %i ",Sigma.N);
    TriDiagonalMatrixToDoubleTriDiagonalMatrix(Sigma, &T, &SigmaDouble, Ef, epsilonAnd);
    FreeMatrix(&T);
    Sigma.N=Min(Sigma.N,Max(4*HSet.NTriDiag,2*G0.N));
    printf("=> %i ",Sigma.N);
    strcpy(Sigma.Name,"Sigma");
    SymmetrizeTridiagonalMatrix(&Sigma);
    printf("=> %i ",Sigma.N);
    Sigma.N=Min(Sigma.N,Max(2*HSet.NTriDiag,2*G0.N));
    printf("=> %i\n",Sigma.N);
    EndTiming("G math Selfenergy",Time);
    
    
    printf("  Start GBathNewTriDiagonalRepresentation\n");
    GBathTriDiagonalRepresentation(G0, Sigma, &GBathNew,Max(2*HSet.NTriDiag,2*G0.N));
    GBathNew.N=Min(GBathNew.N,Max(4*HSet.NTriDiag,2*G0.N));
    strcpy(GBathNew.Name,"Gbath-New");
    SymmetrizeTridiagonalMatrix(&GBathNew);
    GBathNew.N=Min(GBathNew.N,Max(4*HSet.NTriDiag,2*G0.N));
    
    
    /* calculate GBath for the next loop by mixing GBath with GBathNew */
    printf("Start mixing GBath with GBathNew\n");
    TriDiagonalMatrixCopy(&GBathOld, GBath);
    strcpy(GBathOld.Name,"Gbath-Old");
    FreeTriDiagonalMatrix(&GBath);
    GBathOld.b[0]=sqrt(HSet.Mix);
    GBathNew.b[0]=sqrt(1-HSet.Mix);
    SumGreensfunctionsTriDiagonalRepresentation(GBathNew, GBathOld, &GBath);
    strcpy(GBath.Name,"Gbath");
    GBathOld.b[0]=1;
    GBathNew.b[0]=1;
    SymmetrizeTridiagonalMatrix(&GBath);
    GBath.N=Min(GBath.N,Max(4*HSet.NTriDiag,2*G0.N));
    EndTiming("G math",Time);
    /* create spectra */
    printf("Start create spectra\n");
    Spectra[0].N       =HSet.NE;
    Spectra[0].Emin    =HSet.Emin;
    Spectra[0].Emax    =HSet.Emax;
    Spectra[0].dE      =(HSet.Emax-HSet.Emin)/HSet.NE;
    Spectra[0].EminPole=HSet.Emin;
    Spectra[0].EmaxPole=HSet.Emax;
    Spectra[0].Gamma   =HSet.Gamma;
    for(i=0;i<10;i++)
    {
      SpectrumAllocate(Spectra[0], &Spectra[i]);
    }
    TriDiagonalMatrixToSpectrum(0,GBath,&(Spectra[0]));              /*  3 */
    TriDiagonalMatrixToSpectrum(0,GBathOld,&(Spectra[1]));           /*  5 */
    TriDiagonalMatrixToSpectrum(0,GBathNew,&(Spectra[2]));           /*  7 */
    SecondDoubleTriDiagonalMatrixToSpectrum(0,GAnd,&(Spectra[3]));   /*  9 */
    TriDiagonalMatrixToSpectrum(0,G0,&(Spectra[4]));                 /* 11 */
    TriDiagonalMatrixToSpectrum(0,Sigma,&(Spectra[5]));              /* 13 */
    DoubleTriDiagonalMatrixToSpectrum(0,SigmaDouble,&(Spectra[6]));  /* 15 */
    TriDiagonalMatrixToSpectrum(0,GclPup,&(Spectra[7]));             /* 17 */
    TriDiagonalMatrixToSpectrum(0,GclMup,&(Spectra[8]));             /* 19 */
    TriDiagonalMatrixToSpectrum(0,Gcl,&(Spectra[9]));                /* 21 */
    /* write to disk */
    printf("Start write to dist\n");
    if(iOuter==0)
    {
      sprintf(GBathFileName,"%s001",HSet.GRootFileName);
      sprintf(PsiFileName,"%s001",HSet.WaveFunctionRootFileName);
      sprintf(SpecFileName,"%s001",HSet.SpecRootFileName);
    }
    else if(iOuter<9)
    {
      sprintf(GBathFileName,"%s00%1i",HSet.GRootFileName,iOuter+1);
      sprintf(PsiFileName,"%s00%1i",HSet.WaveFunctionRootFileName,iOuter+1);
      sprintf(SpecFileName,"%s00%1i",HSet.SpecRootFileName,iOuter+1);
    }
    else if(iOuter<99)
    {
      sprintf(GBathFileName,"%s0%2i",HSet.GRootFileName,iOuter+1);
      sprintf(PsiFileName,"%s0%2i",HSet.WaveFunctionRootFileName,iOuter+1);
      sprintf(SpecFileName,"%s0%2i",HSet.SpecRootFileName,iOuter+1);
    }
    else if(iOuter<999)
    {
      sprintf(GBathFileName,"%s%3i",HSet.GRootFileName,iOuter+1);
      sprintf(PsiFileName,"%s%3i",HSet.WaveFunctionRootFileName,iOuter+1);
      sprintf(SpecFileName,"%s%3i",HSet.SpecRootFileName,iOuter+1);
    }
    FilePrintSpectra(Spectra, 10, SpecFileName, "w");
    FilePrintPlainTriDiagonalMatrix(GBath,GBathFileName,"w");
    FilePrintWaveFunctionPlainBinary(Psi, PsiFileName, "w");
    /* free variables */
    printf("Start free variables\n");
    for(i=0;i<4;i++)
    {
      FreeOperator(&(VPerturbation[i]));
      FreeTriDiagonalMatrix(&(gPerturbation[i]));
    }
    FreeOperator(&HAnd);
    FreeOperator(&OPESdn);
    FreeOperator(&OPESup);
    FreeOperator(&OIPESdn);
    FreeOperator(&OIPESup);
    FreeTriDiagonalMatrix(&GclMup);
    FreeTriDiagonalMatrix(&GclMdn);
    FreeTriDiagonalMatrix(&GclPup);
    FreeTriDiagonalMatrix(&GclPdn);
    FreeTriDiagonalMatrix(&Gcl);
    FreeTriDiagonalMatrix(&Sigma);
    FreeTriDiagonalMatrix(&GBathNew);
    FreeTriDiagonalMatrix(&GBathOld);
    FreeDoubleTriDiagonalMatrix(&(GAnd.Dt));
    EndTiming("DMFTLoop",Time);
    
    PrintTiming(*Time);
  }
  return 0;
};


int oldmainNoDyson(DMFTSettings HSet, TimeType *Time)
{
  WaveFunctionType Psi,OPsi,PsiTemp;
  TriDiagonalMatrixType G0,GBath,GclPup,GclPdn,GclMup,GclMdn,Gclup,Gcldn,Gcl,Sigma,GBathNew,GBathOld;
  SecondDoubleTriDiagonalMatrixType GAnd;
  OperatorType HAnd, OPESdn,OPESup,OIPESdn,OIPESup;
  MatrixType T;
  SpectraType Spectra[10];
  char GBathFileName[255],PsiFileName[255],SpecFileName[25],Buffer[255];
  uint32_t iOuter,i,PsiNBathV,PsiNBathC;
  double Ef=0.0,epsilonAnd=0.0000001,mu;
  double PsiEnergy,PsiQuality,PsiHist[20],*PsiRho;
  /* read G0 */
  if(FileReadPlainTriDiagonalMatrixType(&G0,HSet.G0FileName)){printf("Fatal Error in main\n Failed to read G0 from file %s\n",HSet.G0FileName);return 1;};
  strcpy(G0.Name,"G0");
  /* read G_Bath_iStart and psi_iStart*/
  if(HSet.iStart==0)
  {
    sprintf(GBathFileName,"%s000",HSet.GRootFileName);
    sprintf(PsiFileName,"%s000",HSet.WaveFunctionRootFileName);
  }
  else if(HSet.iStart<10)
  {
    sprintf(GBathFileName,"%s00%1i",HSet.GRootFileName,HSet.iStart);
    sprintf(PsiFileName,"%s00%1i",HSet.WaveFunctionRootFileName,HSet.iStart);
  }
  else if(HSet.iStart<100)
  {
    sprintf(GBathFileName,"%s0%2i",HSet.GRootFileName,HSet.iStart);
    sprintf(PsiFileName,"%s0%2i",HSet.WaveFunctionRootFileName,HSet.iStart);
  }
  else
  {
    sprintf(GBathFileName,"%s%i",HSet.GRootFileName,HSet.iStart);
    sprintf(PsiFileName,"%s%i",HSet.WaveFunctionRootFileName,HSet.iStart);
  }
  if(FileReadPlainTriDiagonalMatrixType(&GBath,GBathFileName)){printf("Fatal Error in main\n Failed to read GBath from file %s !\n",GBathFileName);return 1;};
  if(file_exists(PsiFileName))
  {
    Psi.NBitsKey=HSet.NBitsKey;
    if(FileReadWaveFunctionPlainBinary(&Psi, PsiFileName)){printf("Fatal Error in main\n Failed to read Psi from file %s !\n",PsiFileName);return 1;};
    printf("\nRead the starting wave function from the file %s\n",PsiFileName);
    sscanf(Psi.Name,"%s %i %i",Buffer,&PsiNBathC,&PsiNBathV);
    printf("Wave-function \"%s\" found has %i conduction bath sites, %i valence bath sites and %i Fermions in basis\n",Psi.Name,PsiNBathC,PsiNBathV,Psi.NFermion);
    if((2*(PsiNBathC+PsiNBathV+2))!=Psi.NFermion)
    {
      printf("The total number of Fermions is incompatible witht the bath sizes\n");
      printf("This wave-function can not be used to start from\n");
      return 1;
    }
  }
  else
  {
    printf("\nWave-function file: '%s' does not excist\n Started from new function\n",PsiFileName);
    CreateWaveFunctionDMFTSecond(&Psi, HSet.NBath, HSet.NBath, HSet.NBitsKey);
  }
  strcpy(GBath.Name,"Gbath");
  /* do the DMFT loops */
  for(iOuter=HSet.iStart;iOuter<HSet.iStart+HSet.Ni;iOuter++)
  {
    BeginTiming("DMFTLoop",Time);
    printf("\n\n***********************************************************************\n");
    printf("Start of loop %i from %i to %i\n",iOuter+1,HSet.iStart,HSet.iStart+HSet.Ni);
    printf("***********************************************************************\n\n\n");
    TriDiagonalMatrixToSecondDoubleTriDiagonalMatrix(GBath, &T, &GAnd, Ef, epsilonAnd);
    
    PrintTriDiagonalMatrix(GBath);
    PrintDoubleTriDiagonalMatrix(GAnd.Dt);
    
    FreeMatrix(&T);
    if(GAnd.Dt.con.N<HSet.NBath)
    {
      printf("Warning: GAnd.Dt.con.N(%i)<HSet.NBath(%i)\n",GAnd.Dt.con.N,HSet.NBath);
    }
    else
    {
      GAnd.Dt.con.N=HSet.NBath;
    }
    if(GAnd.Dt.val.N<HSet.NBath)
    {
      printf("Warning: GAnd.Dt.val.N(%i)<HSet.NBath(%i)\n",GAnd.Dt.val.N,HSet.NBath);
    }
    else
    {
      GAnd.Dt.val.N=HSet.NBath;
    }
    mu=-HSet.U/2;
    SecondDoubleTriDiagonalMatrixToAndersonImpurityOperator(&HAnd, GAnd, GAnd, HSet.U, mu);
    PrintOperatorReadable(HAnd);
    /* test if psi has correct dimensions and correct if needed */
    sscanf(Psi.Name,"%s %i %i",Buffer,&PsiNBathC,&PsiNBathV);
    if(PsiNBathC!=GAnd.Dt.con.N)
    {
      printf("\nBefore diagonalizing:\nThe wavefunction has %i conduction bath sites\nThe Anderson Hamiltonian %i\n",PsiNBathC,GAnd.Dt.con.N);
      if(PsiNBathC>GAnd.Dt.con.N)
      {
        printf("Removing %i unoccupied fermions from the wave-function\n",PsiNBathC-GAnd.Dt.con.N);
        for(i=0;i<PsiNBathC-GAnd.Dt.con.N;i++)
        {
          WaveFunctionRemoveFermion(Psi, &PsiTemp, Psi.NFermion-1);
          WaveFunctionFree(&Psi);
          WaveFunctionRemoveFermion(PsiTemp, &Psi, 4+PsiNBathV+PsiNBathC-i-1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      else
      {
        printf("Adding %i unoccupied fermions to the wave-function\n",GAnd.Dt.con.N-PsiNBathC);
        for(i=0;i<GAnd.Dt.con.N-PsiNBathC;i++)
        {
          WaveFunctionAddFermion(Psi, &PsiTemp, Psi.NFermion,0);
          WaveFunctionFree(&Psi);
          WaveFunctionAddFermion(PsiTemp, &Psi, 4+PsiNBathV+PsiNBathC+i,0);
          WaveFunctionFree(&PsiTemp);
        }
      }
    }
    if(PsiNBathV!=GAnd.Dt.val.N)
    {
      printf("\nBefore diagonalizing:\nThe wavefunction has %i valence bath sites\nThe Anderson Hamiltonian %i\n",PsiNBathV,GAnd.Dt.val.N);
      if(PsiNBathV>GAnd.Dt.val.N)
      {
        printf("Removing %i occupied fermions from the wave-function\n",PsiNBathV-GAnd.Dt.val.N);
        for(i=0;i<PsiNBathV-GAnd.Dt.val.N;i++)
        {
          WaveFunctionRemoveFermion(Psi, &PsiTemp, Psi.NFermion-GAnd.Dt.con.N-1);
          WaveFunctionFree(&Psi);
          WaveFunctionRemoveFermion(PsiTemp, &Psi, 4+PsiNBathV-i-1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      else
      {
        printf("Adding %i occupied fermions to the wave-function\n",GAnd.Dt.val.N-PsiNBathV);
        for(i=0;i<GAnd.Dt.val.N-PsiNBathV;i++)
        {
          WaveFunctionAddFermion(Psi, &PsiTemp, Psi.NFermion-GAnd.Dt.val.N,1);
          WaveFunctionFree(&Psi);
          WaveFunctionAddFermion(PsiTemp, &Psi, 4+PsiNBathV+i,1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      printf("Corrected basis\n");
      PrintWaveFunctionReadableTruncated(Psi, 0.1);
    }
    sprintf(Psi.Name,"DMFT %i %i",GAnd.Dt.con.N,GAnd.Dt.val.N);
    /* create ground-state wave-funciton */
    BeginTiming("Lanczos GroundState",Time);
    LanczosGroundState(HAnd,&Psi,&HSet.zeroLancWavefunction,HSet.zeroLancEnergy,HSet.U,&Psi,0);
    LanczosGroundState(HAnd,&Psi,&HSet.zeroLancWavefunction,HSet.zeroLancEnergy,HSet.U,&Psi,0);
    EndTiming("Lanczos GroundState",Time);
    PsiHermitianOperatorPsi(HAnd, Psi, &PsiEnergy);
    PsiHermitianOperatorSqrPsi(HAnd,Psi,&PsiQuality);
    PsiQuality/=Sqr(PsiEnergy);
    PrintWaveFunctionReadableTruncated(Psi, 0.1);
    CalculatePrefactorHistogram(PsiHist, Psi, 10.0, 20);
    printf("\nHistogram of prefactors in Wave-function\n");
    for(i=0;i<19;i++)
    {
      printf("10^(-%2i) to 10^(-%2i) | %22.15E\n",i,i+1,PsiHist[i]);
    }
    printf("10^(-%2i) to 0        | %22.15E\n",i,PsiHist[19]);
    printf("\n");
    printf("Energy                  <psi|H^2|psi>/(<psi|H|psi>^2)\n");
    printf("%22.15E  %22.15E\n",PsiEnergy,PsiQuality);
    PsiRho=calloc(Psi.NFermion,sizeof(double));
    if(PsiRho==NULL){printf("calloc failed in main for PsiRho\n Fatal Error\n");fflush(stdout);return 1;}
    CalculateDensity(PsiRho,Psi);
    printf("Occupation \n");
    for(i=0;i<Psi.NFermion;i++)
    {
      printf("%5i %22.15E\n",i,PsiRho[i]);
    }
    free(PsiRho);
    /* create spin up and spin down PES and IPES */
    CreatePESdnOperatorDMFT(&OPESdn, HAnd.NFermion);
    CreatePESupOperatorDMFT(&OPESup, HAnd.NFermion);
    CreateIPESdnOperatorDMFT(&OIPESdn, HAnd.NFermion);
    CreateIPESupOperatorDMFT(&OIPESup, HAnd.NFermion);
    GclMup.N=HSet.NTriDiag;
    InitTriDiagonalMatrix(&GclMup);
    GclMdn.N=HSet.NTriDiag;
    InitTriDiagonalMatrix(&GclMdn);
    GclPup.N=HSet.NTriDiag;
    InitTriDiagonalMatrix(&GclPup);
    GclPdn.N=HSet.NTriDiag;
    InitTriDiagonalMatrix(&GclPdn);
    BeginTiming("Lanczos Tri Diagonalization PES up",Time);
    OperatorPsi(OPESup,Psi,&OPsi);
    LanczosTriDiagonalize(HAnd, OPsi, HSet.zeroLancWavefunction,&GclMup,HSet.Gamma);
    EndTiming("Lanczos Tri Diagonalization PES up",Time);
    BeginTiming("Lanczos Tri Diagonalization PES dn",Time);
    OperatorPsi(OPESdn,Psi,&OPsi);
    LanczosTriDiagonalize(HAnd, OPsi, HSet.zeroLancWavefunction,&GclMdn,HSet.Gamma);
    EndTiming("Lanczos Tri Diagonalization PES dn",Time);
    BeginTiming("Lanczos Tri Diagonalization IPES up",Time);
    OperatorPsi(OIPESup,Psi,&OPsi);
    LanczosTriDiagonalize(HAnd, OPsi, HSet.zeroLancWavefunction,&GclPup,HSet.Gamma);
    EndTiming("Lanczos Tri Diagonalization IPES up",Time);
    BeginTiming("Lanczos Tri Diagonalization IPES dn",Time);
    OperatorPsi(OIPESdn,Psi,&OPsi);
    LanczosTriDiagonalize(HAnd, OPsi, HSet.zeroLancWavefunction,&GclPdn,HSet.Gamma);
    EndTiming("Lanczos Tri Diagonalization IPES dn",Time);
    /* G math */
    printf("Start G math\n");
    ShiftTriDiagonalMatrix(&GclMdn, -PsiEnergy);
    ShiftTriDiagonalMatrix(&GclMup, -PsiEnergy);
    ShiftTriDiagonalMatrix(&GclPdn, -PsiEnergy);
    ShiftTriDiagonalMatrix(&GclPup, -PsiEnergy);
    for(i=1;i<=GclMdn.N;i++)
    {
      GclMdn.a[i]*=-1;
    }
    for(i=1;i<=GclMup.N;i++)
    {
      GclMup.a[i]*=-1;
    }
    strcpy(GclPup.Name,"Gcluster-Plus-Up");
    strcpy(GclPdn.Name,"Gcluster-Plus-Dn");
    strcpy(GclMup.Name,"Gcluster-Min-Up");
    strcpy(GclMdn.Name,"Gcluster-Min-Dn");
    printf("  Start SumGreensfunctionsTriDiagonalRepresentation\n");fflush(stdout);
    SumGreensfunctionsTriDiagonalRepresentation(GclPup, GclMup, &Gclup);
    strcpy(Gclup.Name,"Gcluster-Up");
    SumGreensfunctionsTriDiagonalRepresentation(GclPdn, GclMdn, &Gcldn);
    strcpy(Gcldn.Name,"Gcluster-Dn");
    SumGreensfunctionsTriDiagonalRepresentation(Gclup, Gcldn, &Gcl);
    strcpy(Gcl.Name,"Gcluster");
    Gcl.b[0]=1.0;
    // SymmetrizeTridiagonalMatrix(&Gcl);
    Gcl.N=Min(Gcl.N,Max(4*HSet.NTriDiag,G0.N));
    printf("  Start SelfEnergyTriDiagonalRepresentation\n");
    SelfEnergyTriDiagonalRepresentation(G0, Gcl, &Sigma);
    Sigma.N=Min(Sigma.N,Max(4*HSet.NTriDiag,G0.N));
    strcpy(Sigma.Name,"Sigma");
    // SymmetrizeTridiagonalMatrix(&Sigma);
    Sigma.N=Min(Sigma.N,Max(4*HSet.NTriDiag,G0.N));
    
    
    
    printf("  Start GBathNewTriDiagonalRepresentation\n");
    // SumGreensfunctionsTriDiagonalRepresentation(Gcl, Sigma, &GBathNew);
    GBathTriDiagonalRepresentation(G0, Sigma, &GBathNew,Max(4*HSet.NTriDiag,4*G0.N));
    
    
    
    GBathNew.N=Min(GBathNew.N,Max(4*HSet.NTriDiag,G0.N));
    strcpy(GBathNew.Name,"Gbath-New");
    // SymmetrizeTridiagonalMatrix(&GBathNew);
    GBathNew.N=Min(GBathNew.N,Max(4*HSet.NTriDiag,G0.N));
    /* calculate GBath for the next loop by mixing GBath with GBathNew */
    printf("Start mixing GBath with GBathNew\n");
    TriDiagonalMatrixCopy(&GBathOld, GBath);
    strcpy(GBathOld.Name,"Gbath-Old");
    FreeTriDiagonalMatrix(&GBath);
    GBathOld.b[0]=sqrt(HSet.Mix);
    GBathNew.b[0]=sqrt(1-HSet.Mix);
    SumGreensfunctionsTriDiagonalRepresentation(GBathNew, GBathOld, &GBath);
    strcpy(GBath.Name,"Gbath");
    GBathOld.b[0]=1;
    GBathNew.b[0]=1;
    //  SymmetrizeTridiagonalMatrix(&GBath);
    GBath.N=Min(GBath.N,Max(4*HSet.NTriDiag,G0.N));
    /* create spectra */
    printf("Start create spectra\n");
    Spectra[0].N       =HSet.NE;
    Spectra[0].Emin    =HSet.Emin;
    Spectra[0].Emax    =HSet.Emax;
    Spectra[0].dE      =(HSet.Emax-HSet.Emin)/HSet.NE;
    Spectra[0].EminPole=HSet.Emin;
    Spectra[0].EmaxPole=HSet.Emax;
    Spectra[0].Gamma   =HSet.Gamma;
    for(i=0;i<10;i++)
    {
      SpectrumAllocate(Spectra[0], &Spectra[i]);
    }
    TriDiagonalMatrixToSpectrum(0,GBath,&(Spectra[0]));
    TriDiagonalMatrixToSpectrum(0,GBathOld,&(Spectra[1]));
    TriDiagonalMatrixToSpectrum(0,GBathNew,&(Spectra[2]));
    SecondDoubleTriDiagonalMatrixToSpectrum(0,GAnd,&(Spectra[3]));
    TriDiagonalMatrixToSpectrum(0,G0,&(Spectra[4]));
    TriDiagonalMatrixToSpectrum(0,Sigma,&(Spectra[5]));
    //    TriDiagonalMatrixToSpectrum(0,GLocal,&(Spectra[6]));
    TriDiagonalMatrixToSpectrum(0,GclPup,&(Spectra[7]));
    TriDiagonalMatrixToSpectrum(0,GclMup,&(Spectra[8]));
    TriDiagonalMatrixToSpectrum(0,Gcl,&(Spectra[9]));
    /* write to disk */
    printf("Start write to dist\n");
    if(iOuter==0)
    {
      sprintf(GBathFileName,"%s001",HSet.GRootFileName);
      sprintf(PsiFileName,"%s001",HSet.WaveFunctionRootFileName);
      sprintf(SpecFileName,"%s001",HSet.SpecRootFileName);
    }
    else if(iOuter<9)
    {
      sprintf(GBathFileName,"%s00%1i",HSet.GRootFileName,iOuter+1);
      sprintf(PsiFileName,"%s00%1i",HSet.WaveFunctionRootFileName,iOuter+1);
      sprintf(SpecFileName,"%s00%1i",HSet.SpecRootFileName,iOuter+1);
    }
    else if(iOuter<99)
    {
      sprintf(GBathFileName,"%s0%2i",HSet.GRootFileName,iOuter+1);
      sprintf(PsiFileName,"%s0%2i",HSet.WaveFunctionRootFileName,iOuter+1);
      sprintf(SpecFileName,"%s0%2i",HSet.SpecRootFileName,iOuter+1);
    }
    else if(iOuter<999)
    {
      sprintf(GBathFileName,"%s%3i",HSet.GRootFileName,iOuter+1);
      sprintf(PsiFileName,"%s%3i",HSet.WaveFunctionRootFileName,iOuter+1);
      sprintf(SpecFileName,"%s%3i",HSet.SpecRootFileName,iOuter+1);
    }
    FilePrintSpectra(Spectra, 10, SpecFileName, "w");
    FilePrintPlainTriDiagonalMatrix(GBath,GBathFileName,"w");
    FilePrintWaveFunctionPlainBinary(Psi, PsiFileName, "w");
    /* free variables */
    printf("Start free variables\n");
    FreeOperator(&HAnd);
    FreeOperator(&OPESdn);
    FreeOperator(&OPESup);
    FreeOperator(&OIPESdn);
    FreeOperator(&OIPESup);
    FreeTriDiagonalMatrix(&GclMup);
    FreeTriDiagonalMatrix(&GclMdn);
    FreeTriDiagonalMatrix(&GclPup);
    FreeTriDiagonalMatrix(&GclPdn);
    FreeTriDiagonalMatrix(&Gcl);
    //   FreeTriDiagonalMatrix(&GLocal);
    FreeTriDiagonalMatrix(&Sigma);
    FreeTriDiagonalMatrix(&GBathNew);
    FreeTriDiagonalMatrix(&GBathOld);
    FreeDoubleTriDiagonalMatrix(&(GAnd.Dt));
    EndTiming("DMFTLoop",Time);
    PrintTiming(*Time);
  }
  return 0;
};


int old2mainDyson(DMFTSettings HSet, TimeType *Time)
{
  WaveFunctionType Psi,OPsi,PsiTemp;
  TriDiagonalMatrixType G0,GBath,GclPup,GclPdn,GclMup,GclMdn,Gclup,Gcldn,Gcl,Sigma,GLocal,GBathNew,GBathOld,gPerturbation[4];
  DoubleTriDiagonalMatrixType GAnd,GAndTrunk;
  OperatorType HAnd, OPESdn,OPESup,OIPESdn,OIPESup,VPerturbation[4];
  MatrixType T;
  SpectraType Spectra[10];
  char GBathFileName[255],PsiFileName[255],SpecFileName[25],Buffer[255];
  uint32_t iOuter,i,PsiNBathV,PsiNBathC;
  uint16_t oadd[1];
  double Ef=0.0,epsilonAnd=0.1,mu;
  double PsiEnergy,PsiQuality,PsiHist[20],*PsiRho;
  /* read G0 */
  if(FileReadPlainTriDiagonalMatrixType(&G0,HSet.G0FileName)){printf("Fatal Error in main\n Failed to read G0 from file %s\n",HSet.G0FileName);return 1;};
  strcpy(G0.Name,"G0");
  /* read G_Bath_iStart and psi_iStart*/
  if(HSet.iStart==0)
  {
    sprintf(GBathFileName,"%s000",HSet.GRootFileName);
    sprintf(PsiFileName,"%s000",HSet.WaveFunctionRootFileName);
  }
  else if(HSet.iStart<10)
  {
    sprintf(GBathFileName,"%s00%1i",HSet.GRootFileName,HSet.iStart);
    sprintf(PsiFileName,"%s00%1i",HSet.WaveFunctionRootFileName,HSet.iStart);
  }
  else if(HSet.iStart<100)
  {
    sprintf(GBathFileName,"%s0%2i",HSet.GRootFileName,HSet.iStart);
    sprintf(PsiFileName,"%s0%2i",HSet.WaveFunctionRootFileName,HSet.iStart);
  }
  else
  {
    sprintf(GBathFileName,"%s%i",HSet.GRootFileName,HSet.iStart);
    sprintf(PsiFileName,"%s%i",HSet.WaveFunctionRootFileName,HSet.iStart);
  }
  if(FileReadPlainTriDiagonalMatrixType(&GBath,GBathFileName)){printf("Fatal Error in main\n Failed to read GBath from file %s !\n",GBathFileName);return 1;};
  if(file_exists(PsiFileName))
  {
    Psi.NBitsKey=HSet.NBitsKey;
    if(FileReadWaveFunctionPlainBinary(&Psi, PsiFileName)){printf("Fatal Error in main\n Failed to read Psi from file %s !\n",PsiFileName);return 1;};
    printf("\nRead the starting wave function from the file %s\n",PsiFileName);
    sscanf(Psi.Name,"%s %i %i",Buffer,&PsiNBathC,&PsiNBathV);
    printf("Wave-function \"%s\" found has %i conduction bath sites, %i valence bath sites and %i Fermions in basis\n",Psi.Name,PsiNBathC,PsiNBathV,Psi.NFermion);
    if((2*(PsiNBathC+PsiNBathV+1))!=Psi.NFermion)
    {
      printf("The total number of Fermions is incompatible witht the bath sizes\n");
      printf("This wave-function can not be used to start from\n");
      return 1;
    }
  }
  else
  {
    printf("\nWave-function file: '%s' does not excist\n Started from new function\n",PsiFileName);
    CreateWaveFunctionDMFT(&Psi, HSet.NBath, HSet.NBath, HSet.NBitsKey);
  }
  strcpy(GBath.Name,"Gbath");
  /* do the DMFT loops */
  for(iOuter=HSet.iStart;iOuter<HSet.iStart+HSet.Ni;iOuter++)
  {
    BeginTiming("DMFTLoop",Time);
    printf("\n\n***********************************************************************\n");
    printf("Start of loop %i from %i to %i\n",iOuter+1,HSet.iStart,HSet.iStart+HSet.Ni);
    printf("***********************************************************************\n\n\n");
    /* Create splitted anderson matrix */
    TriDiagonalMatrixToDoubleTriDiagonalMatrix(GBath, &T, &GAnd, Ef, epsilonAnd);
    /* Create GAnd and GAndTrunk as a copy pointing to the same memory */
    GAndTrunk.Ef=GAnd.Ef;
    GAndTrunk.Epsilon=GAnd.Epsilon;
    GAndTrunk.a0=GAnd.a0;
    GAndTrunk.b0=GAnd.b0;
    GAndTrunk.val.N=GAnd.val.N;
    GAndTrunk.val.a=GAnd.val.a;
    GAndTrunk.val.b=GAnd.val.b;
    GAndTrunk.con.N=GAnd.con.N;
    GAndTrunk.con.a=GAnd.con.a;
    GAndTrunk.con.b=GAnd.con.b;
    /* Create VPerturbation */
    for(i=0;i<4;i++)
    {
      InitOperator(&(VPerturbation[i]));
      VPerturbation[i].NFermion=2+2*(Min(GAnd.val.N,HSet.NBath)+Min(GAnd.con.N,HSet.NBath));
      VPerturbation[i].NBoson=0;
      VPerturbation[i].NWord=Min(VPerturbation[i].NFermion,1+(VPerturbation[i].NFermion-1)/8);
      VPerturbation[i].QComplex=0;
      InitOperatorFixedLength(&(VPerturbation[i]),1,1);
    }
    if(GAnd.val.N>HSet.NBath)
    {
      oadd[0]=CrAnMask+1+HSet.NBath;
      RealOperatorAddLader(&(VPerturbation[0]), oadd, 1, GAnd.val.b[HSet.NBath]);
      oadd[0]=CrAnMask+1+Min(GAnd.con.N,HSet.NBath)+2*HSet.NBath;
      RealOperatorAddLader(&(VPerturbation[1]), oadd, 1, GAnd.val.b[HSet.NBath]);
    }
    if(GAnd.con.N>HSet.NBath)
    {
      oadd[0]=1+Min(GAnd.val.N,HSet.NBath)+HSet.NBath;
      RealOperatorAddLader(&(VPerturbation[2]), oadd, 1, GAnd.con.b[HSet.NBath]);
      oadd[0]=1+2*Min(GAnd.val.N,HSet.NBath)+2*HSet.NBath;
      RealOperatorAddLader(&(VPerturbation[3]), oadd, 1, GAnd.con.b[HSet.NBath]);
    }
    /* Create gPerturbation */
    if(GAnd.val.N>HSet.NBath)
    {
      gPerturbation[0].N=GAnd.val.N-HSet.NBath;
      gPerturbation[1].N=GAnd.val.N-HSet.NBath;
      InitTriDiagonalMatrix(&(gPerturbation[0]));
      InitTriDiagonalMatrix(&(gPerturbation[1]));
      for(i=1;i<GAnd.val.N-HSet.NBath;i++)
      {
        gPerturbation[0].a[i]=-GAnd.val.a[i+HSet.NBath];
        gPerturbation[0].b[i]=GAnd.val.b[i+HSet.NBath];
        gPerturbation[1].a[i]=-GAnd.val.a[i+HSet.NBath];
        gPerturbation[1].b[i]=GAnd.val.b[i+HSet.NBath];
      }
      gPerturbation[0].a[GAnd.val.N-HSet.NBath]=-GAnd.val.a[GAnd.val.N];
      gPerturbation[1].a[GAnd.val.N-HSet.NBath]=-GAnd.val.a[GAnd.val.N];
    }
    else
    {
      gPerturbation[0].N=1;
      gPerturbation[1].N=1;
      InitTriDiagonalMatrix(&(gPerturbation[0]));
      InitTriDiagonalMatrix(&(gPerturbation[1]));
      gPerturbation[0].a[1]=0;
      gPerturbation[1].a[1]=0;
    }
    if(GAnd.con.N>HSet.NBath)
    {
      gPerturbation[2].N=GAnd.con.N-HSet.NBath;
      gPerturbation[3].N=GAnd.con.N-HSet.NBath;
      InitTriDiagonalMatrix(&(gPerturbation[2]));
      InitTriDiagonalMatrix(&(gPerturbation[3]));
      for(i=1;i<GAnd.con.N-HSet.NBath;i++)
      {
        gPerturbation[2].a[i]=GAnd.con.a[i+HSet.NBath];
        gPerturbation[2].b[i]=GAnd.con.b[i+HSet.NBath];
        gPerturbation[3].a[i]=GAnd.con.a[i+HSet.NBath];
        gPerturbation[3].b[i]=GAnd.con.b[i+HSet.NBath];
      }
      gPerturbation[2].a[GAnd.con.N-HSet.NBath]=GAnd.con.a[GAnd.con.N];
      gPerturbation[3].a[GAnd.con.N-HSet.NBath]=GAnd.con.a[GAnd.con.N];
    }
    else
    {
      gPerturbation[2].N=1;
      gPerturbation[3].N=1;
      InitTriDiagonalMatrix(&(gPerturbation[2]));
      InitTriDiagonalMatrix(&(gPerturbation[3]));
      gPerturbation[2].a[1]=0;
      gPerturbation[3].a[1]=0;
    };
    /* Create HAnd */
    FreeMatrix(&T);
    if(GAnd.con.N<HSet.NBath)
    {
      printf("Warning: GAnd.con.N(%i)<HSet.NBath(%i)\n",GAnd.con.N,HSet.NBath);
    }
    else
    {
      GAndTrunk.con.N=HSet.NBath;
    }
    if(GAnd.val.N<HSet.NBath)
    {
      printf("Warning: GAnd.val.N(%i)<HSet.NBath(%i)\n",GAnd.val.N,HSet.NBath);
    }
    else
    {
      GAndTrunk.val.N=HSet.NBath;
    }
    mu=-HSet.U/2;
    DoubleTriDiagonalMatrixToAndersonImpurityOperator(&HAnd, GAndTrunk, GAndTrunk, HSet.U, mu);
    /* test if psi has correct dimensions and correct if needed */
    sscanf(Psi.Name,"%s %i %i",Buffer,&PsiNBathC,&PsiNBathV);
    if(PsiNBathC!=GAndTrunk.con.N)
    {
      printf("\nBefore diagonalizing:\nThe wavefunction has %i conduction bath sites\nThe Anderson Hamiltonian %i\n",PsiNBathC,GAndTrunk.con.N);
      if(PsiNBathC>GAndTrunk.con.N)
      {
        printf("Removing %i unoccupied fermions from the wave-function\n",PsiNBathC-GAndTrunk.con.N);
        for(i=0;i<PsiNBathC-GAndTrunk.con.N;i++)
        {
          WaveFunctionRemoveFermion(Psi, &PsiTemp, Psi.NFermion-1);
          WaveFunctionFree(&Psi);
          WaveFunctionRemoveFermion(PsiTemp, &Psi, 2+PsiNBathV+PsiNBathC-i-1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      else
      {
        printf("Adding %i unoccupied fermions to the wave-function\n",GAndTrunk.con.N-PsiNBathC);
        for(i=0;i<GAndTrunk.con.N-PsiNBathC;i++)
        {
          WaveFunctionAddFermion(Psi, &PsiTemp, Psi.NFermion,0);
          WaveFunctionFree(&Psi);
          WaveFunctionAddFermion(PsiTemp, &Psi, 2+PsiNBathV+PsiNBathC+i,0);
          WaveFunctionFree(&PsiTemp);
        }
      }
    }
    if(PsiNBathV!=GAndTrunk.val.N)
    {
      printf("\nBefore diagonalizing:\nThe wavefunction has %i valence bath sites\nThe Anderson Hamiltonian %i\n",PsiNBathV,GAndTrunk.val.N);
      if(PsiNBathV>GAndTrunk.val.N)
      {
        printf("Removing %i occupied fermions from the wave-function\n",PsiNBathV-GAndTrunk.val.N);
        for(i=0;i<PsiNBathV-GAndTrunk.val.N;i++)
        {
          WaveFunctionRemoveFermion(Psi, &PsiTemp, Psi.NFermion-GAndTrunk.con.N-1);
          WaveFunctionFree(&Psi);
          WaveFunctionRemoveFermion(PsiTemp, &Psi, 2+PsiNBathV-i-1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      else
      {
        printf("Adding %i occupied fermions to the wave-function\n",GAndTrunk.val.N-PsiNBathV);
        for(i=0;i<GAndTrunk.val.N-PsiNBathV;i++)
        {
          WaveFunctionAddFermion(Psi, &PsiTemp, Psi.NFermion-GAndTrunk.val.N,1);
          WaveFunctionFree(&Psi);
          WaveFunctionAddFermion(PsiTemp, &Psi, 2+PsiNBathV+i,1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      printf("Corrected basis\n");
      PrintWaveFunctionReadableTruncated(Psi, 0.1);
    }
    sprintf(Psi.Name,"DMFT %i %i",GAndTrunk.con.N,GAndTrunk.val.N);
    /* create ground-state wave-function */
    BeginTiming("Lanczos GroundState",Time);
    LanczosGroundState(HAnd,&Psi,&HSet.zeroLancWavefunction,HSet.zeroLancEnergy,HSet.U,&Psi,0);
    LanczosGroundState(HAnd,&Psi,&HSet.zeroLancWavefunction,HSet.zeroLancEnergy,HSet.U,&Psi,0);
    EndTiming("Lanczos GroundState",Time);
    PsiHermitianOperatorPsi(HAnd, Psi, &PsiEnergy);
    PsiHermitianOperatorSqrPsi(HAnd,Psi,&PsiQuality);
    PsiQuality/=Sqr(PsiEnergy);
    PrintWaveFunctionReadableTruncated(Psi, 0.1);
    CalculatePrefactorHistogram(PsiHist, Psi, 10.0, 20);
    printf("\nHistogram of prefactors in Wave-function\n");
    for(i=0;i<19;i++)
    {
      printf("10^(-%2i) to 10^(-%2i) | %22.15E\n",i,i+1,PsiHist[i]);
    }
    printf("10^(-%2i) to 0        | %22.15E\n",i,PsiHist[19]);
    printf("\n");
    printf("Energy                  <psi|H^2|psi>/(<psi|H|psi>^2)\n");
    printf("%22.15E  %22.15E\n",PsiEnergy,PsiQuality);
    PsiRho=calloc(Psi.NFermion,sizeof(double));
    if(PsiRho==NULL){printf("calloc failed in main for PsiRho\n Fatal Error\n");fflush(stdout);return 1;}
    CalculateDensity(PsiRho,Psi);
    printf("Occupation \n");
    for(i=0;i<Psi.NFermion;i++)
    {
      printf("%5i %22.15E\n",i,PsiRho[i]);
    }
    free(PsiRho);
    /*                                           */
    /* create spin up and spin down PES and IPES */
    /*                                           */
    CreatePESdnOperatorDMFT(&OPESdn, HAnd.NFermion);
    CreatePESupOperatorDMFT(&OPESup, HAnd.NFermion);
    CreateIPESdnOperatorDMFT(&OIPESdn, HAnd.NFermion);
    CreateIPESupOperatorDMFT(&OIPESup, HAnd.NFermion);
    GclMup.N=Max(HSet.NTriDiag,G0.N);
    InitTriDiagonalMatrix(&GclMup);
    GclMdn.N=Max(HSet.NTriDiag,G0.N);
    InitTriDiagonalMatrix(&GclMdn);
    GclPup.N=Max(HSet.NTriDiag,G0.N);
    InitTriDiagonalMatrix(&GclPup);
    GclPdn.N=Max(HSet.NTriDiag,G0.N);
    InitTriDiagonalMatrix(&GclPdn);
    BeginTiming("Lanczos Tri Diagonalization Dyson",Time);
    OperatorPsi(OPESup,Psi,&OPsi);
    LanczosTriDiagonalizeDyson(HAnd, OPsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation, gPerturbation, 4, &GclMup);
    EndTiming("Lanczos Tri Diagonalization Dyson",Time);
    BeginTiming("Lanczos Tri Diagonalization Dyson",Time);
    OperatorPsi(OPESdn,Psi,&OPsi);
    LanczosTriDiagonalizeDyson(HAnd, OPsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation, gPerturbation, 4, &GclMdn);
    EndTiming("Lanczos Tri Diagonalization Dyson",Time);
    BeginTiming("Lanczos Tri Diagonalization Dyson",Time);
    OperatorPsi(OIPESup,Psi,&OPsi);
    LanczosTriDiagonalizeDyson(HAnd, OPsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation, gPerturbation, 4, &GclPup);
    EndTiming("Lanczos Tri Diagonalization Dyson",Time);
    BeginTiming("Lanczos Tri Diagonalization Dyson",Time);
    OperatorPsi(OIPESdn,Psi,&OPsi);
    LanczosTriDiagonalizeDyson(HAnd, OPsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation, gPerturbation, 4, &GclPdn);
    EndTiming("Lanczos Tri Diagonalization Dyson",Time);
    /* G math */
    printf("Start G math\n");
    ShiftTriDiagonalMatrix(&GclMdn, -PsiEnergy);
    ShiftTriDiagonalMatrix(&GclMup, -PsiEnergy);
    ShiftTriDiagonalMatrix(&GclPdn, -PsiEnergy);
    ShiftTriDiagonalMatrix(&GclPup, -PsiEnergy);
    for(i=1;i<=GclMdn.N;i++)
    {
      GclMdn.a[i]*=-1;
    }
    for(i=1;i<=GclMup.N;i++)
    {
      GclMup.a[i]*=-1;
    }
    strcpy(GclPup.Name,"Gcluster-Plus-Up");
    strcpy(GclPdn.Name,"Gcluster-Plus-Dn");
    strcpy(GclMup.Name,"Gcluster-Min-Up");
    strcpy(GclMdn.Name,"Gcluster-Min-Dn");
    printf("  Start SumGreensfunctionsTriDiagonalRepresentation\n");fflush(stdout);
    SumGreensfunctionsTriDiagonalRepresentation(GclPup, GclMup, &Gclup);
    strcpy(Gclup.Name,"Gcluster-Up");
    SumGreensfunctionsTriDiagonalRepresentation(GclPdn, GclMdn, &Gcldn);
    strcpy(Gcldn.Name,"Gcluster-Dn");
    SumGreensfunctionsTriDiagonalRepresentation(Gclup, Gcldn, &Gcl);
    strcpy(Gcl.Name,"Gcluster");
    Gcl.b[0]=1.0;
    printf("Gcl.N = %i ",Gcl.N);
    SymmetrizeTridiagonalMatrix(&Gcl);
    Gcl.N=Min(Gcl.N,Max(4*HSet.NTriDiag,2*G0.N));
    printf("=> %i\n",Gcl.N);
    printf("  Start SelfEnergyTriDiagonalRepresentation\n");
    printf("GBath.N = %i ",GBath.N);
    SelfEnergyTriDiagonalRepresentation(GBath, Gcl, &Sigma);
    printf("Sigma.N = %i ",Sigma.N);
    Sigma.N=Min(Sigma.N,Max(4*HSet.NTriDiag,2*G0.N));
    printf("=> %i ",Sigma.N);
    strcpy(Sigma.Name,"Sigma");
    SymmetrizeTridiagonalMatrix(&Sigma);
    printf("=> %i ",Sigma.N);
    Sigma.N=Min(Sigma.N,Max(4*HSet.NTriDiag,2*G0.N));
    printf("=> %i\n",Sigma.N);
    printf("  Start GLocalTriDiagonalRepresentation\n");
    GLocalTriDiagonalRepresentation(G0, Sigma, &GLocal, Sigma.N+G0.N-1);
    printf("Glocal.N = %i ",GLocal.N);
    GLocal.N=Min(GLocal.N,Max(4*HSet.NTriDiag,2*G0.N));
    printf("=> %i ",GLocal.N);
    strcpy(GLocal.Name,"Glocal");
    SymmetrizeTridiagonalMatrix(&GLocal);
    printf("=> %i ",GLocal.N);
    GLocal.N=Min(GLocal.N,Max(4*HSet.NTriDiag,2*G0.N));
    printf("=> %i\n",GLocal.N);
    printf("  Start GBathNewTriDiagonalRepresentation\n");
    GBathNewTriDiagonalRepresentation(GLocal, Sigma, &GBathNew);
    GBathNew.N=Min(GBathNew.N,Max(4*HSet.NTriDiag,2*G0.N));
    strcpy(GBathNew.Name,"Gbath-New");
    SymmetrizeTridiagonalMatrix(&GBathNew);
    GBathNew.N=Min(GBathNew.N,Max(4*HSet.NTriDiag,2*G0.N));
    /* calculate GBath for the next loop by mixing GBath with GBathNew */
    printf("Start mixing GBath with GBathNew\n");
    TriDiagonalMatrixCopy(&GBathOld, GBath);
    strcpy(GBathOld.Name,"Gbath-Old");
    FreeTriDiagonalMatrix(&GBath);
    GBathOld.b[0]=sqrt(HSet.Mix);
    GBathNew.b[0]=sqrt(1-HSet.Mix);
    SumGreensfunctionsTriDiagonalRepresentation(GBathNew, GBathOld, &GBath);
    strcpy(GBath.Name,"Gbath");
    GBathOld.b[0]=1;
    GBathNew.b[0]=1;
    SymmetrizeTridiagonalMatrix(&GBath);
    GBath.N=Min(GBath.N,Max(4*HSet.NTriDiag,2*G0.N));
    /* create spectra */
    printf("Start create spectra\n");
    Spectra[0].N       =HSet.NE;
    Spectra[0].Emin    =HSet.Emin;
    Spectra[0].Emax    =HSet.Emax;
    Spectra[0].dE      =(HSet.Emax-HSet.Emin)/HSet.NE;
    Spectra[0].EminPole=HSet.Emin;
    Spectra[0].EmaxPole=HSet.Emax;
    Spectra[0].Gamma   =HSet.Gamma;
    for(i=0;i<10;i++)
    {
      SpectrumAllocate(Spectra[0], &Spectra[i]);
    }
    TriDiagonalMatrixToSpectrum(0,GBath,&(Spectra[0]));        /*  3 */
    TriDiagonalMatrixToSpectrum(0,GBathOld,&(Spectra[1]));     /*  5 */
    TriDiagonalMatrixToSpectrum(0,GBathNew,&(Spectra[2]));     /*  7 */
    DoubleTriDiagonalMatrixToSpectrum(0,GAnd,&(Spectra[3]));   /*  9 */
    TriDiagonalMatrixToSpectrum(0,G0,&(Spectra[4]));           /* 11 */
    TriDiagonalMatrixToSpectrum(0,Sigma,&(Spectra[5]));        /* 13 */
    TriDiagonalMatrixToSpectrum(0,GLocal,&(Spectra[6]));       /* 15 */
    TriDiagonalMatrixToSpectrum(0,Gclup,&(Spectra[7]));        /* 17 */
    TriDiagonalMatrixToSpectrum(0,Gcldn,&(Spectra[8]));        /* 19 */
    TriDiagonalMatrixToSpectrum(0,Gcl,&(Spectra[9]));          /* 21 */
    /* write to disk */
    printf("Start write to dist\n");
    if(iOuter==0)
    {
      sprintf(GBathFileName,"%s001",HSet.GRootFileName);
      sprintf(PsiFileName,"%s001",HSet.WaveFunctionRootFileName);
      sprintf(SpecFileName,"%s001",HSet.SpecRootFileName);
    }
    else if(iOuter<9)
    {
      sprintf(GBathFileName,"%s00%1i",HSet.GRootFileName,iOuter+1);
      sprintf(PsiFileName,"%s00%1i",HSet.WaveFunctionRootFileName,iOuter+1);
      sprintf(SpecFileName,"%s00%1i",HSet.SpecRootFileName,iOuter+1);
    }
    else if(iOuter<99)
    {
      sprintf(GBathFileName,"%s0%2i",HSet.GRootFileName,iOuter+1);
      sprintf(PsiFileName,"%s0%2i",HSet.WaveFunctionRootFileName,iOuter+1);
      sprintf(SpecFileName,"%s0%2i",HSet.SpecRootFileName,iOuter+1);
    }
    else if(iOuter<999)
    {
      sprintf(GBathFileName,"%s%3i",HSet.GRootFileName,iOuter+1);
      sprintf(PsiFileName,"%s%3i",HSet.WaveFunctionRootFileName,iOuter+1);
      sprintf(SpecFileName,"%s%3i",HSet.SpecRootFileName,iOuter+1);
    }
    FilePrintSpectra(Spectra, 10, SpecFileName, "w");
    FilePrintPlainTriDiagonalMatrix(GBath,GBathFileName,"w");
    FilePrintWaveFunctionPlainBinary(Psi, PsiFileName, "w");
    /* free variables */
    printf("Start free variables\n");
    for(i=0;i<4;i++)
    {
      FreeOperator(&(VPerturbation[i]));
      FreeTriDiagonalMatrix(&(gPerturbation[i]));
    }
    FreeOperator(&HAnd);
    FreeOperator(&OPESdn);
    FreeOperator(&OPESup);
    FreeOperator(&OIPESdn);
    FreeOperator(&OIPESup);
    FreeTriDiagonalMatrix(&GclMup);
    FreeTriDiagonalMatrix(&GclMdn);
    FreeTriDiagonalMatrix(&GclPup);
    FreeTriDiagonalMatrix(&GclPdn);
    FreeTriDiagonalMatrix(&Gcl);
    FreeTriDiagonalMatrix(&GLocal);
    FreeTriDiagonalMatrix(&Sigma);
    FreeTriDiagonalMatrix(&GBathNew);
    FreeTriDiagonalMatrix(&GBathOld);
    FreeDoubleTriDiagonalMatrix(&GAnd);
    EndTiming("DMFTLoop",Time);
    
    PrintTiming(*Time);
  }
  return 0;
};

int oldmainDysonRho(int argc, const char * argv[])
{
  DMFTSettings HSet;
  WaveFunctionType Psi,OPsi,PsiTemp;
  TriDiagonalMatrixType G0,GBath,GclPup,GclPdn,GclMup,GclMdn,Gclup,Gcldn,Gcl,Sigma,GLocal,GBathNew,GBathOld,Gtemp,Gperturb[4];
  DoubleTriDiagonalMatrixType GAnd,GAndOrig;
  OperatorType HAnd, OPESdn,OPESup,OIPESdn,OIPESup;
  TimeType Time;
  MatrixType T,Rho,RhoBound;
  SpectraType Spectra[10];
  char GBathFileName[255],PsiFileName[255],SpecFileName[25],Buffer[255];
  uint32_t iOuter,i,PsiNBathV,PsiNBathC;
  double Ef=0.0,epsilonAnd=0.0000001,mu;
  double PsiEnergy,PsiQuality,PsiHist[20],*PsiRho;
  InitializeTiming(&Time);
  BeginTiming("Total",&Time);
  TimeStampAll();
  argc-=1;
  /* read the DMFT settings */
  if(FileReadDMFTSettings(&HSet,argv[argc])){printf("Fatal Error in main\n Failed to read DMFT Settings from file %s\n",argv[argc]);return 1;};
  PrintDMFTSettings(HSet);
  /* read G0 */
  if(FileReadPlainTriDiagonalMatrixType(&G0,HSet.G0FileName)){printf("Fatal Error in main\n Failed to read G0 from file %s\n",HSet.G0FileName);return 1;};
  strcpy(G0.Name,"G0");
  /* read G_Bath_iStart and psi_iStart*/
  if(HSet.iStart==0)
  {
    sprintf(GBathFileName,"%s000",HSet.GRootFileName);
    sprintf(PsiFileName,"%s000",HSet.WaveFunctionRootFileName);
  }
  else if(HSet.iStart<10)
  {
    sprintf(GBathFileName,"%s00%1i",HSet.GRootFileName,HSet.iStart);
    sprintf(PsiFileName,"%s00%1i",HSet.WaveFunctionRootFileName,HSet.iStart);
  }
  else if(HSet.iStart<100)
  {
    sprintf(GBathFileName,"%s0%2i",HSet.GRootFileName,HSet.iStart);
    sprintf(PsiFileName,"%s0%2i",HSet.WaveFunctionRootFileName,HSet.iStart);
  }
  else
  {
    sprintf(GBathFileName,"%s%i",HSet.GRootFileName,HSet.iStart);
    sprintf(PsiFileName,"%s%i",HSet.WaveFunctionRootFileName,HSet.iStart);
  }
  if(FileReadPlainTriDiagonalMatrixType(&GBath,GBathFileName)){printf("Fatal Error in main\n Failed to read GBath from file %s !\n",GBathFileName);return 1;};
  if(file_exists(PsiFileName))
  {
    Psi.NBitsKey=HSet.NBitsKey;
    if(FileReadWaveFunctionPlainBinary(&Psi, PsiFileName)){printf("Fatal Error in main\n Failed to read Psi from file %s !\n",PsiFileName);return 1;};
    printf("\nRead the starting wave function from the file %s\n",PsiFileName);
    sscanf(Psi.Name,"%s %i %i",Buffer,&PsiNBathC,&PsiNBathV);
    printf("Wave-function \"%s\" found has %i conduction bath sites, %i valence bath sites and %i Fermions in basis\n",Psi.Name,PsiNBathC,PsiNBathV,Psi.NFermion);
    if((2*(PsiNBathC+PsiNBathV+1))!=Psi.NFermion)
    {
      printf("The total number of Fermions is incompatible witht the bath sizes\n");
      printf("This wave-function can not be used to start from\n");
      return 1;
    }
  }
  else
  {
    printf("\nWave-function file: '%s' does not excist\n Started from new function\n",PsiFileName);
    CreateWaveFunctionDMFT(&Psi, HSet.NBath, HSet.NBath, HSet.NBitsKey);
  }
  strcpy(GBath.Name,"Gbath");
  /* allocate Rho matrices */
  Rho.QComplex=0;
  Rho.Nj=2*(HSet.NBath+1);
  Rho.Ni=HSet.NTriDiag;
  if(InitMatrix(&Rho)){printf("InitMatrix(&Rho) failed in main\n");fflush(stdout);return 1;};
  RhoBound.QComplex=0;
  RhoBound.Nj=4;
  RhoBound.Ni=HSet.NTriDiag;
  if(InitMatrix(&RhoBound)){printf("InitMatrix(&RhoBound) failed in main\n");fflush(stdout);return 1;};
  /* do the DMFT loops */
  for(iOuter=HSet.iStart;iOuter<HSet.iStart+HSet.Ni;iOuter++)
  {
    BeginTiming("DMFTLoop",&Time);
    printf("\n\n***********************************************************************\n");
    printf("Start of loop %i from %i to %i\n",iOuter+1,HSet.iStart,HSet.iStart+HSet.Ni);
    printf("***********************************************************************\n\n\n");
    /* create GBath as a double tri-diagonal matrix */
    TriDiagonalMatrixToDoubleTriDiagonalMatrix(GBath, &T, &GAnd, Ef, epsilonAnd);
    FreeMatrix(&T);
    /* make GAndOrig point to GAnd with original size */
    GAndOrig=GAnd;
    /* reduce size of GBath for exact diagonalization */
    if(GAnd.con.N<HSet.NBath)
    {
      printf("Warning: GAnd.con.N(%i)<HSet.NBath(%i)\n",GAnd.con.N,HSet.NBath);
    }
    else
    {
      GAnd.con.N=HSet.NBath;
    }
    if(GAnd.val.N<HSet.NBath)
    {
      printf("Warning: GAnd.val.N(%i)<HSet.NBath(%i)\n",GAnd.val.N,HSet.NBath);
    }
    else
    {
      GAnd.val.N=HSet.NBath;
    }
    mu=-HSet.U/2;
    DoubleTriDiagonalMatrixToAndersonImpurityOperator(&HAnd, GAnd, GAnd, HSet.U, mu);
    /* test if psi has correct dimensions and correct if needed */
    sscanf(Psi.Name,"%s %i %i",Buffer,&PsiNBathC,&PsiNBathV);
    if(PsiNBathC!=GAnd.con.N)
    {
      printf("\nBefore diagonalizing:\nThe wavefunction has %i conduction bath sites\nThe Anderson Hamiltonian %i\n",PsiNBathC,GAnd.con.N);
      if(PsiNBathC>GAnd.con.N)
      {
        printf("Removing %i unoccupied fermions from the wave-function\n",PsiNBathC-GAnd.con.N);
        for(i=0;i<PsiNBathC-GAnd.con.N;i++)
        {
          WaveFunctionRemoveFermion(Psi, &PsiTemp, Psi.NFermion-1);
          WaveFunctionFree(&Psi);
          WaveFunctionRemoveFermion(PsiTemp, &Psi, 2+PsiNBathV+PsiNBathC-i-1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      else
      {
        printf("Adding %i unoccupied fermions to the wave-function\n",GAnd.con.N-PsiNBathC);
        for(i=0;i<GAnd.con.N-PsiNBathC;i++)
        {
          WaveFunctionAddFermion(Psi, &PsiTemp, Psi.NFermion,0);
          WaveFunctionFree(&Psi);
          WaveFunctionAddFermion(PsiTemp, &Psi, 2+PsiNBathV+PsiNBathC+i,0);
          WaveFunctionFree(&PsiTemp);
        }
      }
    }
    if(PsiNBathV!=GAnd.val.N)
    {
      printf("\nBefore diagonalizing:\nThe wavefunction has %i valence bath sites\nThe Anderson Hamiltonian %i\n",PsiNBathV,GAnd.val.N);
      if(PsiNBathV>GAnd.val.N)
      {
        printf("Removing %i occupied fermions from the wave-function\n",PsiNBathV-GAnd.val.N);
        for(i=0;i<PsiNBathV-GAnd.val.N;i++)
        {
          WaveFunctionRemoveFermion(Psi, &PsiTemp, Psi.NFermion-GAnd.con.N-1);
          WaveFunctionFree(&Psi);
          WaveFunctionRemoveFermion(PsiTemp, &Psi, 2+PsiNBathV-i-1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      else
      {
        printf("Adding %i occupied fermions to the wave-function\n",GAnd.val.N-PsiNBathV);
        for(i=0;i<GAnd.val.N-PsiNBathV;i++)
        {
          WaveFunctionAddFermion(Psi, &PsiTemp, Psi.NFermion-GAnd.val.N,1);
          WaveFunctionFree(&Psi);
          WaveFunctionAddFermion(PsiTemp, &Psi, 2+PsiNBathV+i,1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      printf("Corrected basis\n");
      PrintWaveFunctionReadableTruncated(Psi, 0.1);
    }
    sprintf(Psi.Name,"DMFT %i %i",GAnd.con.N,GAnd.val.N);
    /* create ground-state wave-function */
    BeginTiming("Lanczos GroundState",&Time);
    LanczosGroundState(HAnd,&Psi,&HSet.zeroLancWavefunction,HSet.zeroLancEnergy,HSet.U,&Psi,0);
    LanczosGroundState(HAnd,&Psi,&HSet.zeroLancWavefunction,HSet.zeroLancEnergy,HSet.U,&Psi,0);
    EndTiming("Lanczos GroundState",&Time);
    PsiHermitianOperatorPsi(HAnd, Psi, &PsiEnergy);
    PsiHermitianOperatorSqrPsi(HAnd,Psi,&PsiQuality);
    PsiQuality/=Sqr(PsiEnergy);
    PrintWaveFunctionReadableTruncated(Psi, 0.1);
    CalculatePrefactorHistogram(PsiHist, Psi, 10.0, 20);
    printf("\nHistogram of prefactors in Wave-function\n");
    for(i=0;i<19;i++)
    {
      printf("10^(-%2i) to 10^(-%2i) | %22.15E\n",i,i+1,PsiHist[i]);
    }
    printf("10^(-%2i) to 0        | %22.15E\n",i,PsiHist[19]);
    printf("\n");
    printf("Energy                  <psi|H^2|psi>/(<psi|H|psi>^2)\n");
    printf("%22.15E  %22.15E\n",PsiEnergy,PsiQuality);
    PsiRho=calloc(Psi.NFermion,sizeof(double));
    if(PsiRho==NULL){printf("calloc failed in main for PsiRho\n Fatal Error\n");fflush(stdout);return 1;}
    CalculateDensity(PsiRho,Psi);
    printf("Occupation \n");
    for(i=0;i<Psi.NFermion;i++)
    {
      printf("%5i %22.15E\n",i,PsiRho[i]);
    }
    free(PsiRho);
    /* create Gperturb */
    if(GAndOrig.val.N<=HSet.NBath)
    {
      Gperturb[0].N=0;
      Gperturb[2].N=0;
    }
    else
    {
      Gperturb[0].N=GAndOrig.val.N-HSet.NBath;
      Gperturb[0].a=&GAndOrig.val.a[HSet.NBath];
      Gperturb[0].b=&GAndOrig.val.b[HSet.NBath];
      Gperturb[2].N=GAndOrig.val.N-HSet.NBath;
      Gperturb[2].a=&GAndOrig.val.a[HSet.NBath];
      Gperturb[2].b=&GAndOrig.val.b[HSet.NBath];
    }
    if(GAndOrig.con.N<=HSet.NBath)
    {
      Gperturb[1].N=0;
      Gperturb[3].N=0;
    }
    else
    {
      Gperturb[1].N=GAndOrig.con.N-HSet.NBath;
      Gperturb[1].a=&GAndOrig.con.a[HSet.NBath];
      Gperturb[1].b=&GAndOrig.con.b[HSet.NBath];
      Gperturb[3].N=GAndOrig.con.N-HSet.NBath;
      Gperturb[3].a=&GAndOrig.con.a[HSet.NBath];
      Gperturb[3].b=&GAndOrig.con.b[HSet.NBath];
    }
    
    
    
    /* create spin up and spin down PES and IPES */
    CreatePESdnOperatorDMFT(&OPESdn, HAnd.NFermion);
    CreatePESupOperatorDMFT(&OPESup, HAnd.NFermion);
    CreateIPESdnOperatorDMFT(&OIPESdn, HAnd.NFermion);
    CreateIPESupOperatorDMFT(&OIPESup, HAnd.NFermion);
    GclMup.N=HSet.NTriDiag;
    InitTriDiagonalMatrix(&GclMup);
    GclMdn.N=HSet.NTriDiag;
    InitTriDiagonalMatrix(&GclMdn);
    GclPup.N=HSet.NTriDiag;
    InitTriDiagonalMatrix(&GclPup);
    GclPdn.N=HSet.NTriDiag;
    InitTriDiagonalMatrix(&GclPdn);
    Gtemp.N=HSet.NTriDiag;
    InitTriDiagonalMatrix(&Gtemp);
    BeginTiming("Lanczos Tri Diagonalization PES up",&Time);
    OperatorPsi(OPESup,Psi,&OPsi);
    LanczosTriDiagonalizePlusDensity(HAnd, OPsi, HSet.zeroLancWavefunction,&Gtemp,&Rho,HSet.Gamma);
    /* copy the relavant elements from rho to rhobound */
    for(i=0;i<Rho.Ni;i++)
    {
      RhoBound.VRe[i][0]=1-Rho.VRe[i][GAnd.val.N+1];
      RhoBound.VRe[i][1]=Rho.VRe[i][GAnd.val.N+GAnd.con.N+1];
      RhoBound.VRe[i][2]=1-Rho.VRe[i][GAnd.val.N+GAnd.con.N+GAnd.val.N+1];
      RhoBound.VRe[i][3]=Rho.VRe[i][GAnd.val.N+GAnd.con.N+GAnd.val.N+GAnd.con.N+1];
    }
    PrintMatrix(RhoBound);
    /* use Dyson to add more bath sites to spectrum than ground-state */
    TriDiagonalMatrixDysonScalar(Gtemp, &GclMup, RhoBound, 4, Gperturb);
    EndTiming("Lanczos Tri Diagonalization PES up",&Time);
    BeginTiming("Lanczos Tri Diagonalization PES dn",&Time);
    OperatorPsi(OPESdn,Psi,&OPsi);
    LanczosTriDiagonalizePlusDensity(HAnd, OPsi, HSet.zeroLancWavefunction,&Gtemp,&Rho,HSet.Gamma);
    /* copy the relavant elements from rho to rhobound */
    for(i=0;i<Rho.Ni;i++)
    {
      RhoBound.VRe[i][0]=1-Rho.VRe[i][GAnd.val.N+1];
      RhoBound.VRe[i][1]=Rho.VRe[i][GAnd.val.N+GAnd.con.N+1];
      RhoBound.VRe[i][2]=1-Rho.VRe[i][GAnd.val.N+GAnd.con.N+GAnd.val.N+1];
      RhoBound.VRe[i][3]=Rho.VRe[i][GAnd.val.N+GAnd.con.N+GAnd.val.N+GAnd.con.N+1];
    }
    PrintMatrix(RhoBound);
    /* use Dyson to add more bath sites to spectrum than ground-state */
    TriDiagonalMatrixDysonScalar(Gtemp, &GclMdn, RhoBound, 4, Gperturb);
    PrintTriDiagonalMatrix(Gtemp);
    PrintTriDiagonalMatrix(GclMdn);
    EndTiming("Lanczos Tri Diagonalization PES dn",&Time);
    BeginTiming("Lanczos Tri Diagonalization IPES up",&Time);
    OperatorPsi(OIPESup,Psi,&OPsi);
    LanczosTriDiagonalizePlusDensity(HAnd, OPsi, HSet.zeroLancWavefunction,&Gtemp,&Rho,HSet.Gamma);
    /* copy the relavant elements from rho to rhobound */
    for(i=0;i<Rho.Ni;i++)
    {
      RhoBound.VRe[i][0]=1-Rho.VRe[i][GAnd.val.N+1];
      RhoBound.VRe[i][1]=Rho.VRe[i][GAnd.val.N+GAnd.con.N+1];
      RhoBound.VRe[i][2]=1-Rho.VRe[i][GAnd.val.N+GAnd.con.N+GAnd.val.N+1];
      RhoBound.VRe[i][3]=Rho.VRe[i][GAnd.val.N+GAnd.con.N+GAnd.val.N+GAnd.con.N+1];
    }
    PrintMatrix(RhoBound);
    /* use Dyson to add more bath sites to spectrum than ground-state */
    TriDiagonalMatrixDysonScalar(Gtemp, &GclPup, RhoBound, 4, Gperturb);
    EndTiming("Lanczos Tri Diagonalization IPES up",&Time);
    BeginTiming("Lanczos Tri Diagonalization IPES dn",&Time);
    OperatorPsi(OIPESdn,Psi,&OPsi);
    LanczosTriDiagonalizePlusDensity(HAnd, OPsi, HSet.zeroLancWavefunction,&Gtemp,&Rho,HSet.Gamma);
    /* copy the relavant elements from rho to rhobound */
    for(i=0;i<Rho.Ni;i++)
    {
      RhoBound.VRe[i][0]=1-Rho.VRe[i][GAnd.val.N+1];
      RhoBound.VRe[i][1]=Rho.VRe[i][GAnd.val.N+GAnd.con.N+1];
      RhoBound.VRe[i][2]=1-Rho.VRe[i][GAnd.val.N+GAnd.con.N+GAnd.val.N+1];
      RhoBound.VRe[i][3]=Rho.VRe[i][GAnd.val.N+GAnd.con.N+GAnd.val.N+GAnd.con.N+1];
    }
    PrintMatrix(RhoBound);
    /* use Dyson to add more bath sites to spectrum than ground-state */
    TriDiagonalMatrixDysonScalar(Gtemp, &GclPdn, RhoBound, 4, Gperturb);
    PrintTriDiagonalMatrix(Gtemp);
    PrintTriDiagonalMatrix(GclPdn);
    EndTiming("Lanczos Tri Diagonalization IPES dn",&Time);
    /* G math */
    printf("Start G math\n");
    ShiftTriDiagonalMatrix(&GclMdn, -PsiEnergy);
    ShiftTriDiagonalMatrix(&GclMup, -PsiEnergy);
    ShiftTriDiagonalMatrix(&GclPdn, -PsiEnergy);
    ShiftTriDiagonalMatrix(&GclPup, -PsiEnergy);
    for(i=1;i<=GclMdn.N;i++)
    {
      GclMdn.a[i]*=-1;
    }
    for(i=1;i<=GclMup.N;i++)
    {
      GclMup.a[i]*=-1;
    }
    strcpy(GclPup.Name,"Gcluster-Plus-Up");
    strcpy(GclPdn.Name,"Gcluster-Plus-Dn");
    strcpy(GclMup.Name,"Gcluster-Min-Up");
    strcpy(GclMdn.Name,"Gcluster-Min-Dn");
    printf("  Start SumGreensfunctionsTriDiagonalRepresentation\n");fflush(stdout);
    SumGreensfunctionsTriDiagonalRepresentation(GclPup, GclMup, &Gclup);
    strcpy(Gclup.Name,"Gcluster-Up");
    SumGreensfunctionsTriDiagonalRepresentation(GclPdn, GclMdn, &Gcldn);
    strcpy(Gcldn.Name,"Gcluster-Dn");
    SumGreensfunctionsTriDiagonalRepresentation(Gclup, Gcldn, &Gcl);
    strcpy(Gcl.Name,"Gcluster");
    Gcl.b[0]=1.0;
    SymmetrizeTridiagonalMatrix(&Gcl);
    Gcl.N=Min(Gcl.N,4*HSet.NTriDiag);
    printf("  Start SelfEnergyTriDiagonalRepresentation\n");
    SelfEnergyTriDiagonalRepresentation(GBath, Gcl, &Sigma);
    Sigma.N=Min(Sigma.N,4*HSet.NTriDiag);
    strcpy(Sigma.Name,"Sigma");
    SymmetrizeTridiagonalMatrix(&Sigma);
    Sigma.N=Min(Sigma.N,4*HSet.NTriDiag);
    printf("  Start GLocalTriDiagonalRepresentation\n");
    GLocalTriDiagonalRepresentation(G0, Sigma, &GLocal, HSet.NTriDiag);
    GLocal.N=Min(GLocal.N,4*HSet.NTriDiag);
    strcpy(GLocal.Name,"Glocal");
    SymmetrizeTridiagonalMatrix(&GLocal);
    GLocal.N=Min(GLocal.N,4*HSet.NTriDiag);
    printf("  Start GBathNewTriDiagonalRepresentation\n");
    GBathNewTriDiagonalRepresentation(GLocal, Sigma, &GBathNew);
    GBathNew.N=Min(GBathNew.N,4*HSet.NTriDiag);
    strcpy(GBathNew.Name,"Gbath-New");
    SymmetrizeTridiagonalMatrix(&GBathNew);
    GBathNew.N=Min(GBathNew.N,4*HSet.NTriDiag);
    /* calculate GBath for the next loop by mixing GBath with GBathNew */
    printf("Start mixing GBath with GBathNew\n");
    TriDiagonalMatrixCopy(&GBathOld, GBath);
    strcpy(GBathOld.Name,"Gbath-Old");
    FreeTriDiagonalMatrix(&GBath);
    GBathOld.b[0]=sqrt(HSet.Mix);
    GBathNew.b[0]=sqrt(1-HSet.Mix);
    SumGreensfunctionsTriDiagonalRepresentation(GBathNew, GBathOld, &GBath);
    strcpy(GBath.Name,"Gbath");
    GBathOld.b[0]=1;
    GBathNew.b[0]=1;
    SymmetrizeTridiagonalMatrix(&GBath);
    GBath.N=Min(GBath.N,4*HSet.NTriDiag);
    /* create spectra */
    printf("Start create spectra\n");
    Spectra[0].N       =HSet.NE;
    Spectra[0].Emin    =HSet.Emin;
    Spectra[0].Emax    =HSet.Emax;
    Spectra[0].dE      =(HSet.Emax-HSet.Emin)/HSet.NE;
    Spectra[0].EminPole=HSet.Emin;
    Spectra[0].EmaxPole=HSet.Emax;
    Spectra[0].Gamma   =HSet.Gamma;
    for(i=0;i<10;i++)
    {
      SpectrumAllocate(Spectra[0], &Spectra[i]);
    }
    TriDiagonalMatrixToSpectrum(0,GBath,&(Spectra[0]));
    TriDiagonalMatrixToSpectrum(0,GBathOld,&(Spectra[1]));
    TriDiagonalMatrixToSpectrum(0,GBathNew,&(Spectra[2]));
    DoubleTriDiagonalMatrixToSpectrum(0,GAnd,&(Spectra[3]));
    TriDiagonalMatrixToSpectrum(0,G0,&(Spectra[4]));
    TriDiagonalMatrixToSpectrum(0,Sigma,&(Spectra[5]));
    TriDiagonalMatrixToSpectrum(0,GLocal,&(Spectra[6]));
    TriDiagonalMatrixToSpectrum(0,Gclup,&(Spectra[7]));
    TriDiagonalMatrixToSpectrum(0,Gcldn,&(Spectra[8]));
    TriDiagonalMatrixToSpectrum(0,Gcl,&(Spectra[9]));
    /* write to disk */
    printf("Start write to dist\n");
    if(iOuter==0)
    {
      sprintf(GBathFileName,"%s001",HSet.GRootFileName);
      sprintf(PsiFileName,"%s001",HSet.WaveFunctionRootFileName);
      sprintf(SpecFileName,"%s001",HSet.SpecRootFileName);
    }
    else if(iOuter<9)
    {
      sprintf(GBathFileName,"%s00%1i",HSet.GRootFileName,iOuter+1);
      sprintf(PsiFileName,"%s00%1i",HSet.WaveFunctionRootFileName,iOuter+1);
      sprintf(SpecFileName,"%s00%1i",HSet.SpecRootFileName,iOuter+1);
    }
    else if(iOuter<99)
    {
      sprintf(GBathFileName,"%s0%2i",HSet.GRootFileName,iOuter+1);
      sprintf(PsiFileName,"%s0%2i",HSet.WaveFunctionRootFileName,iOuter+1);
      sprintf(SpecFileName,"%s0%2i",HSet.SpecRootFileName,iOuter+1);
    }
    else if(iOuter<999)
    {
      sprintf(GBathFileName,"%s%3i",HSet.GRootFileName,iOuter+1);
      sprintf(PsiFileName,"%s%3i",HSet.WaveFunctionRootFileName,iOuter+1);
      sprintf(SpecFileName,"%s%3i",HSet.SpecRootFileName,iOuter+1);
    }
    FilePrintSpectra(Spectra, 10, SpecFileName, "w");
    FilePrintPlainTriDiagonalMatrix(GBath,GBathFileName,"w");
    FilePrintWaveFunctionPlainBinary(Psi, PsiFileName, "w");
    /* free variables */
    printf("Start free variables\n");
    FreeOperator(&HAnd);
    FreeOperator(&OPESdn);
    FreeOperator(&OPESup);
    FreeOperator(&OIPESdn);
    FreeOperator(&OIPESup);
    FreeTriDiagonalMatrix(&GclMup);
    FreeTriDiagonalMatrix(&GclMdn);
    FreeTriDiagonalMatrix(&GclPup);
    FreeTriDiagonalMatrix(&GclPdn);
    FreeTriDiagonalMatrix(&Gtemp);
    FreeTriDiagonalMatrix(&Gcl);
    FreeTriDiagonalMatrix(&GLocal);
    FreeTriDiagonalMatrix(&Sigma);
    FreeTriDiagonalMatrix(&GBathNew);
    FreeTriDiagonalMatrix(&GBathOld);
    FreeDoubleTriDiagonalMatrix(&GAnd);
    EndTiming("DMFTLoop",&Time);
  }
  EndTiming("Total",&Time);
  PrintTiming(Time);
  return 0;
};

int old2mainNoDyson(DMFTSettings HSet, TimeType *Time)
{
  WaveFunctionType Psi,OPsi,PsiTemp;
  TriDiagonalMatrixType G0,GBath,GclPup,GclPdn,GclMup,GclMdn,Gclup,Gcldn,Gcl,Sigma,GLocal,GBathNew,GBathOld;
  DoubleTriDiagonalMatrixType GAnd;
  OperatorType HAnd, OPESdn,OPESup,OIPESdn,OIPESup;
  MatrixType T;
  SpectraType Spectra[10];
  char GBathFileName[255],PsiFileName[255],SpecFileName[25],Buffer[255];
  uint32_t iOuter,i,PsiNBathV,PsiNBathC;
  double Ef=0.0,epsilonAnd=0.0000001,mu;
  double PsiEnergy,PsiQuality,PsiHist[20],*PsiRho;
  /* read G0 */
  if(FileReadPlainTriDiagonalMatrixType(&G0,HSet.G0FileName)){printf("Fatal Error in main\n Failed to read G0 from file %s\n",HSet.G0FileName);return 1;};
  strcpy(G0.Name,"G0");
  /* read G_Bath_iStart and psi_iStart*/
  if(HSet.iStart==0)
  {
    sprintf(GBathFileName,"%s000",HSet.GRootFileName);
    sprintf(PsiFileName,"%s000",HSet.WaveFunctionRootFileName);
  }
  else if(HSet.iStart<10)
  {
    sprintf(GBathFileName,"%s00%1i",HSet.GRootFileName,HSet.iStart);
    sprintf(PsiFileName,"%s00%1i",HSet.WaveFunctionRootFileName,HSet.iStart);
  }
  else if(HSet.iStart<100)
  {
    sprintf(GBathFileName,"%s0%2i",HSet.GRootFileName,HSet.iStart);
    sprintf(PsiFileName,"%s0%2i",HSet.WaveFunctionRootFileName,HSet.iStart);
  }
  else
  {
    sprintf(GBathFileName,"%s%i",HSet.GRootFileName,HSet.iStart);
    sprintf(PsiFileName,"%s%i",HSet.WaveFunctionRootFileName,HSet.iStart);
  }
  if(FileReadPlainTriDiagonalMatrixType(&GBath,GBathFileName)){printf("Fatal Error in main\n Failed to read GBath from file %s !\n",GBathFileName);return 1;};
  if(file_exists(PsiFileName))
  {
    Psi.NBitsKey=HSet.NBitsKey;
    if(FileReadWaveFunctionPlainBinary(&Psi, PsiFileName)){printf("Fatal Error in main\n Failed to read Psi from file %s !\n",PsiFileName);return 1;};
    printf("\nRead the starting wave function from the file %s\n",PsiFileName);
    sscanf(Psi.Name,"%s %i %i",Buffer,&PsiNBathC,&PsiNBathV);
    printf("Wave-function \"%s\" found has %i conduction bath sites, %i valence bath sites and %i Fermions in basis\n",Psi.Name,PsiNBathC,PsiNBathV,Psi.NFermion);
    if((2*(PsiNBathC+PsiNBathV+1))!=Psi.NFermion)
    {
      printf("The total number of Fermions is incompatible witht the bath sizes\n");
      printf("This wave-function can not be used to start from\n");
      return 1;
    }
  }
  else
  {
    printf("\nWave-function file: '%s' does not excist\n Started from new function\n",PsiFileName);
    CreateWaveFunctionDMFT(&Psi, HSet.NBath, HSet.NBath, HSet.NBitsKey);
  }
  strcpy(GBath.Name,"Gbath");
  /* do the DMFT loops */
  for(iOuter=HSet.iStart;iOuter<HSet.iStart+HSet.Ni;iOuter++)
  {
    BeginTiming("DMFTLoop",Time);
    printf("\n\n***********************************************************************\n");
    printf("Start of loop %i from %i to %i\n",iOuter+1,HSet.iStart,HSet.iStart+HSet.Ni);
    printf("***********************************************************************\n\n\n");
    TriDiagonalMatrixToDoubleTriDiagonalMatrix(GBath, &T, &GAnd, Ef, epsilonAnd);
    FreeMatrix(&T);
    if(GAnd.con.N<HSet.NBath)
    {
      printf("Warning: GAnd.con.N(%i)<HSet.NBath(%i)\n",GAnd.con.N,HSet.NBath);
    }
    else
    {
      GAnd.con.N=HSet.NBath;
    }
    if(GAnd.val.N<HSet.NBath)
    {
      printf("Warning: GAnd.val.N(%i)<HSet.NBath(%i)\n",GAnd.val.N,HSet.NBath);
    }
    else
    {
      GAnd.val.N=HSet.NBath;
    }
    mu=-HSet.U/2;
    DoubleTriDiagonalMatrixToAndersonImpurityOperator(&HAnd, GAnd, GAnd, HSet.U, mu);
    /* test if psi has correct dimensions and correct if needed */
    sscanf(Psi.Name,"%s %i %i",Buffer,&PsiNBathC,&PsiNBathV);
    if(PsiNBathC!=GAnd.con.N)
    {
      printf("\nBefore diagonalizing:\nThe wavefunction has %i conduction bath sites\nThe Anderson Hamiltonian %i\n",PsiNBathC,GAnd.con.N);
      if(PsiNBathC>GAnd.con.N)
      {
        printf("Removing %i unoccupied fermions from the wave-function\n",PsiNBathC-GAnd.con.N);
        for(i=0;i<PsiNBathC-GAnd.con.N;i++)
        {
          WaveFunctionRemoveFermion(Psi, &PsiTemp, Psi.NFermion-1);
          WaveFunctionFree(&Psi);
          WaveFunctionRemoveFermion(PsiTemp, &Psi, 2+PsiNBathV+PsiNBathC-i-1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      else
      {
        printf("Adding %i unoccupied fermions to the wave-function\n",GAnd.con.N-PsiNBathC);
        for(i=0;i<GAnd.con.N-PsiNBathC;i++)
        {
          WaveFunctionAddFermion(Psi, &PsiTemp, Psi.NFermion,0);
          WaveFunctionFree(&Psi);
          WaveFunctionAddFermion(PsiTemp, &Psi, 2+PsiNBathV+PsiNBathC+i,0);
          WaveFunctionFree(&PsiTemp);
        }
      }
    }
    if(PsiNBathV!=GAnd.val.N)
    {
      printf("\nBefore diagonalizing:\nThe wavefunction has %i valence bath sites\nThe Anderson Hamiltonian %i\n",PsiNBathV,GAnd.val.N);
      if(PsiNBathV>GAnd.val.N)
      {
        printf("Removing %i occupied fermions from the wave-function\n",PsiNBathV-GAnd.val.N);
        for(i=0;i<PsiNBathV-GAnd.val.N;i++)
        {
          WaveFunctionRemoveFermion(Psi, &PsiTemp, Psi.NFermion-GAnd.con.N-1);
          WaveFunctionFree(&Psi);
          WaveFunctionRemoveFermion(PsiTemp, &Psi, 2+PsiNBathV-i-1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      else
      {
        printf("Adding %i occupied fermions to the wave-function\n",GAnd.val.N-PsiNBathV);
        for(i=0;i<GAnd.val.N-PsiNBathV;i++)
        {
          WaveFunctionAddFermion(Psi, &PsiTemp, Psi.NFermion-GAnd.val.N,1);
          WaveFunctionFree(&Psi);
          WaveFunctionAddFermion(PsiTemp, &Psi, 2+PsiNBathV+i,1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      printf("Corrected basis\n");
      PrintWaveFunctionReadableTruncated(Psi, 0.1);
    }
    sprintf(Psi.Name,"DMFT %i %i",GAnd.con.N,GAnd.val.N);
    /* create ground-state wave-funciton */
    BeginTiming("Lanczos GroundState",Time);
    LanczosGroundState(HAnd,&Psi,&HSet.zeroLancWavefunction,HSet.zeroLancEnergy,HSet.U,&Psi,0);
    LanczosGroundState(HAnd,&Psi,&HSet.zeroLancWavefunction,HSet.zeroLancEnergy,HSet.U,&Psi,0);
    EndTiming("Lanczos GroundState",Time);
    PsiHermitianOperatorPsi(HAnd, Psi, &PsiEnergy);
    PsiHermitianOperatorSqrPsi(HAnd,Psi,&PsiQuality);
    PsiQuality/=Sqr(PsiEnergy);
    PrintWaveFunctionReadableTruncated(Psi, 0.1);
    CalculatePrefactorHistogram(PsiHist, Psi, 10.0, 20);
    printf("\nHistogram of prefactors in Wave-function\n");
    for(i=0;i<19;i++)
    {
      printf("10^(-%2i) to 10^(-%2i) | %22.15E\n",i,i+1,PsiHist[i]);
    }
    printf("10^(-%2i) to 0        | %22.15E\n",i,PsiHist[19]);
    printf("\n");
    printf("Energy                  <psi|H^2|psi>/(<psi|H|psi>^2)\n");
    printf("%22.15E  %22.15E\n",PsiEnergy,PsiQuality);
    PsiRho=calloc(Psi.NFermion,sizeof(double));
    if(PsiRho==NULL){printf("calloc failed in main for PsiRho\n Fatal Error\n");fflush(stdout);return 1;}
    CalculateDensity(PsiRho,Psi);
    printf("Occupation \n");
    for(i=0;i<Psi.NFermion;i++)
    {
      printf("%5i %22.15E\n",i,PsiRho[i]);
    }
    free(PsiRho);
    /* create spin up and spin down PES and IPES */
    CreatePESdnOperatorDMFT(&OPESdn, HAnd.NFermion);
    CreatePESupOperatorDMFT(&OPESup, HAnd.NFermion);
    CreateIPESdnOperatorDMFT(&OIPESdn, HAnd.NFermion);
    CreateIPESupOperatorDMFT(&OIPESup, HAnd.NFermion);
    GclMup.N=HSet.NTriDiag;
    InitTriDiagonalMatrix(&GclMup);
    GclMdn.N=HSet.NTriDiag;
    InitTriDiagonalMatrix(&GclMdn);
    GclPup.N=HSet.NTriDiag;
    InitTriDiagonalMatrix(&GclPup);
    GclPdn.N=HSet.NTriDiag;
    InitTriDiagonalMatrix(&GclPdn);
    BeginTiming("Lanczos Tri Diagonalization PES up",Time);
    OperatorPsi(OPESup,Psi,&OPsi);
    LanczosTriDiagonalize(HAnd, OPsi, HSet.zeroLancWavefunction,&GclMup,HSet.Gamma);
    EndTiming("Lanczos Tri Diagonalization PES up",Time);
    BeginTiming("Lanczos Tri Diagonalization PES dn",Time);
    OperatorPsi(OPESdn,Psi,&OPsi);
    LanczosTriDiagonalize(HAnd, OPsi, HSet.zeroLancWavefunction,&GclMdn,HSet.Gamma);
    EndTiming("Lanczos Tri Diagonalization PES dn",Time);
    BeginTiming("Lanczos Tri Diagonalization IPES up",Time);
    OperatorPsi(OIPESup,Psi,&OPsi);
    LanczosTriDiagonalize(HAnd, OPsi, HSet.zeroLancWavefunction,&GclPup,HSet.Gamma);
    EndTiming("Lanczos Tri Diagonalization IPES up",Time);
    BeginTiming("Lanczos Tri Diagonalization IPES dn",Time);
    OperatorPsi(OIPESdn,Psi,&OPsi);
    LanczosTriDiagonalize(HAnd, OPsi, HSet.zeroLancWavefunction,&GclPdn,HSet.Gamma);
    EndTiming("Lanczos Tri Diagonalization IPES dn",Time);
    /* G math */
    printf("Start G math\n");
    ShiftTriDiagonalMatrix(&GclMdn, -PsiEnergy);
    ShiftTriDiagonalMatrix(&GclMup, -PsiEnergy);
    ShiftTriDiagonalMatrix(&GclPdn, -PsiEnergy);
    ShiftTriDiagonalMatrix(&GclPup, -PsiEnergy);
    for(i=1;i<=GclMdn.N;i++)
    {
      GclMdn.a[i]*=-1;
    }
    for(i=1;i<=GclMup.N;i++)
    {
      GclMup.a[i]*=-1;
    }
    strcpy(GclPup.Name,"Gcluster-Plus-Up");
    strcpy(GclPdn.Name,"Gcluster-Plus-Dn");
    strcpy(GclMup.Name,"Gcluster-Min-Up");
    strcpy(GclMdn.Name,"Gcluster-Min-Dn");
    printf("  Start SumGreensfunctionsTriDiagonalRepresentation\n");fflush(stdout);
    SumGreensfunctionsTriDiagonalRepresentation(GclPup, GclMup, &Gclup);
    strcpy(Gclup.Name,"Gcluster-Up");
    SumGreensfunctionsTriDiagonalRepresentation(GclPdn, GclMdn, &Gcldn);
    strcpy(Gcldn.Name,"Gcluster-Dn");
    SumGreensfunctionsTriDiagonalRepresentation(Gclup, Gcldn, &Gcl);
    strcpy(Gcl.Name,"Gcluster");
    Gcl.b[0]=1.0;
    SymmetrizeTridiagonalMatrix(&Gcl);
    Gcl.N=Min(Gcl.N,Max(4*HSet.NTriDiag,G0.N));
    printf("  Start SelfEnergyTriDiagonalRepresentation\n");
    SelfEnergyTriDiagonalRepresentation(GBath, Gcl, &Sigma);
    Sigma.N=Min(Sigma.N,Max(4*HSet.NTriDiag,G0.N));
    strcpy(Sigma.Name,"Sigma");
    SymmetrizeTridiagonalMatrix(&Sigma);
    Sigma.N=Min(Sigma.N,Max(4*HSet.NTriDiag,G0.N));
    printf("  Start GLocalTriDiagonalRepresentation\n");
    GLocalTriDiagonalRepresentation(G0, Sigma, &GLocal, HSet.NTriDiag);
    GLocal.N=Min(GLocal.N,Max(4*HSet.NTriDiag,G0.N));
    strcpy(GLocal.Name,"Glocal");
    SymmetrizeTridiagonalMatrix(&GLocal);
    GLocal.N=Min(GLocal.N,Max(4*HSet.NTriDiag,G0.N));
    printf("  Start GBathNewTriDiagonalRepresentation\n");
    GBathNewTriDiagonalRepresentation(GLocal, Sigma, &GBathNew);
    GBathNew.N=Min(GBathNew.N,Max(4*HSet.NTriDiag,G0.N));
    strcpy(GBathNew.Name,"Gbath-New");
    SymmetrizeTridiagonalMatrix(&GBathNew);
    GBathNew.N=Min(GBathNew.N,Max(4*HSet.NTriDiag,G0.N));
    /* calculate GBath for the next loop by mixing GBath with GBathNew */
    printf("Start mixing GBath with GBathNew\n");
    TriDiagonalMatrixCopy(&GBathOld, GBath);
    strcpy(GBathOld.Name,"Gbath-Old");
    FreeTriDiagonalMatrix(&GBath);
    GBathOld.b[0]=sqrt(HSet.Mix);
    GBathNew.b[0]=sqrt(1-HSet.Mix);
    SumGreensfunctionsTriDiagonalRepresentation(GBathNew, GBathOld, &GBath);
    strcpy(GBath.Name,"Gbath");
    GBathOld.b[0]=1;
    GBathNew.b[0]=1;
    SymmetrizeTridiagonalMatrix(&GBath);
    GBath.N=Min(GBath.N,Max(4*HSet.NTriDiag,G0.N));
    /* create spectra */
    printf("Start create spectra\n");
    Spectra[0].N       =HSet.NE;
    Spectra[0].Emin    =HSet.Emin;
    Spectra[0].Emax    =HSet.Emax;
    Spectra[0].dE      =(HSet.Emax-HSet.Emin)/HSet.NE;
    Spectra[0].EminPole=HSet.Emin;
    Spectra[0].EmaxPole=HSet.Emax;
    Spectra[0].Gamma   =HSet.Gamma;
    for(i=0;i<10;i++)
    {
      SpectrumAllocate(Spectra[0], &Spectra[i]);
    }
    TriDiagonalMatrixToSpectrum(0,GBath,&(Spectra[0]));
    TriDiagonalMatrixToSpectrum(0,GBathOld,&(Spectra[1]));
    TriDiagonalMatrixToSpectrum(0,GBathNew,&(Spectra[2]));
    DoubleTriDiagonalMatrixToSpectrum(0,GAnd,&(Spectra[3]));
    TriDiagonalMatrixToSpectrum(0,G0,&(Spectra[4]));
    TriDiagonalMatrixToSpectrum(0,Sigma,&(Spectra[5]));
    TriDiagonalMatrixToSpectrum(0,GLocal,&(Spectra[6]));
    TriDiagonalMatrixToSpectrum(0,Gclup,&(Spectra[7]));
    TriDiagonalMatrixToSpectrum(0,Gcldn,&(Spectra[8]));
    TriDiagonalMatrixToSpectrum(0,Gcl,&(Spectra[9]));
    /* write to disk */
    printf("Start write to dist\n");
    if(iOuter==0)
    {
      sprintf(GBathFileName,"%s001",HSet.GRootFileName);
      sprintf(PsiFileName,"%s001",HSet.WaveFunctionRootFileName);
      sprintf(SpecFileName,"%s001",HSet.SpecRootFileName);
    }
    else if(iOuter<9)
    {
      sprintf(GBathFileName,"%s00%1i",HSet.GRootFileName,iOuter+1);
      sprintf(PsiFileName,"%s00%1i",HSet.WaveFunctionRootFileName,iOuter+1);
      sprintf(SpecFileName,"%s00%1i",HSet.SpecRootFileName,iOuter+1);
    }
    else if(iOuter<99)
    {
      sprintf(GBathFileName,"%s0%2i",HSet.GRootFileName,iOuter+1);
      sprintf(PsiFileName,"%s0%2i",HSet.WaveFunctionRootFileName,iOuter+1);
      sprintf(SpecFileName,"%s0%2i",HSet.SpecRootFileName,iOuter+1);
    }
    else if(iOuter<999)
    {
      sprintf(GBathFileName,"%s%3i",HSet.GRootFileName,iOuter+1);
      sprintf(PsiFileName,"%s%3i",HSet.WaveFunctionRootFileName,iOuter+1);
      sprintf(SpecFileName,"%s%3i",HSet.SpecRootFileName,iOuter+1);
    }
    FilePrintSpectra(Spectra, 10, SpecFileName, "w");
    FilePrintPlainTriDiagonalMatrix(GBath,GBathFileName,"w");
    FilePrintWaveFunctionPlainBinary(Psi, PsiFileName, "w");
    /* free variables */
    printf("Start free variables\n");
    FreeOperator(&HAnd);
    FreeOperator(&OPESdn);
    FreeOperator(&OPESup);
    FreeOperator(&OIPESdn);
    FreeOperator(&OIPESup);
    FreeTriDiagonalMatrix(&GclMup);
    FreeTriDiagonalMatrix(&GclMdn);
    FreeTriDiagonalMatrix(&GclPup);
    FreeTriDiagonalMatrix(&GclPdn);
    FreeTriDiagonalMatrix(&Gcl);
    FreeTriDiagonalMatrix(&GLocal);
    FreeTriDiagonalMatrix(&Sigma);
    FreeTriDiagonalMatrix(&GBathNew);
    FreeTriDiagonalMatrix(&GBathOld);
    FreeDoubleTriDiagonalMatrix(&GAnd);
    EndTiming("DMFTLoop",Time);
    PrintTiming(*Time);
  }
  return 0;
};

int oldmainNoDysonSecond(DMFTSettings HSet, TimeType *Time)
{
  WaveFunctionType Psi,OPsi,PsiTemp;
  TriDiagonalMatrixType G0,GBath,GclPup,GclPdn,GclMup,GclMdn,Gclup,Gcldn,Gcl,Sigma,GLocal,GBathNew,GBathOld;
  SecondDoubleTriDiagonalMatrixType GAnd;
  OperatorType HAnd, OPESdn,OPESup,OIPESdn,OIPESup;
  MatrixType T;
  SpectraType Spectra[10];
  char GBathFileName[255],PsiFileName[255],SpecFileName[25],Buffer[255];
  uint32_t iOuter,i,PsiNBathV,PsiNBathC;
  double Ef=0.0,epsilonAnd=0.0000001,mu;
  double PsiEnergy,PsiQuality,PsiHist[20],*PsiRho;
  /* read G0 */
  if(FileReadPlainTriDiagonalMatrixType(&G0,HSet.G0FileName)){printf("Fatal Error in main\n Failed to read G0 from file %s\n",HSet.G0FileName);return 1;};
  strcpy(G0.Name,"G0");
  /* read G_Bath_iStart and psi_iStart*/
  if(HSet.iStart==0)
  {
    sprintf(GBathFileName,"%s000",HSet.GRootFileName);
    sprintf(PsiFileName,"%s000",HSet.WaveFunctionRootFileName);
  }
  else if(HSet.iStart<10)
  {
    sprintf(GBathFileName,"%s00%1i",HSet.GRootFileName,HSet.iStart);
    sprintf(PsiFileName,"%s00%1i",HSet.WaveFunctionRootFileName,HSet.iStart);
  }
  else if(HSet.iStart<100)
  {
    sprintf(GBathFileName,"%s0%2i",HSet.GRootFileName,HSet.iStart);
    sprintf(PsiFileName,"%s0%2i",HSet.WaveFunctionRootFileName,HSet.iStart);
  }
  else
  {
    sprintf(GBathFileName,"%s%i",HSet.GRootFileName,HSet.iStart);
    sprintf(PsiFileName,"%s%i",HSet.WaveFunctionRootFileName,HSet.iStart);
  }
  if(FileReadPlainTriDiagonalMatrixType(&GBath,GBathFileName)){printf("Fatal Error in main\n Failed to read GBath from file %s !\n",GBathFileName);return 1;};
  if(file_exists(PsiFileName))
  {
    Psi.NBitsKey=HSet.NBitsKey;
    if(FileReadWaveFunctionPlainBinary(&Psi, PsiFileName)){printf("Fatal Error in main\n Failed to read Psi from file %s !\n",PsiFileName);return 1;};
    printf("\nRead the starting wave function from the file %s\n",PsiFileName);
    sscanf(Psi.Name,"%s %i %i",Buffer,&PsiNBathC,&PsiNBathV);
    printf("Wave-function \"%s\" found has %i conduction bath sites, %i valence bath sites and %i Fermions in basis\n",Psi.Name,PsiNBathC,PsiNBathV,Psi.NFermion);
    if((2*(PsiNBathC+PsiNBathV+2))!=Psi.NFermion)
    {
      printf("The total number of Fermions is incompatible witht the bath sizes\n");
      printf("This wave-function can not be used to start from\n");
      return 1;
    }
  }
  else
  {
    printf("\nWave-function file: '%s' does not excist\n Started from new function\n",PsiFileName);
    CreateWaveFunctionDMFTSecond(&Psi, HSet.NBath, HSet.NBath, HSet.NBitsKey);
  }
  strcpy(GBath.Name,"Gbath");
  /* do the DMFT loops */
  for(iOuter=HSet.iStart;iOuter<HSet.iStart+HSet.Ni;iOuter++)
  {
    BeginTiming("DMFTLoop",Time);
    printf("\n\n***********************************************************************\n");
    printf("Start of loop %i from %i to %i\n",iOuter+1,HSet.iStart,HSet.iStart+HSet.Ni);
    printf("***********************************************************************\n\n\n");
    TriDiagonalMatrixToSecondDoubleTriDiagonalMatrix(GBath, &T, &GAnd, Ef, epsilonAnd);
    
    PrintTriDiagonalMatrix(GBath);
    PrintDoubleTriDiagonalMatrix(GAnd.Dt);
    
    FreeMatrix(&T);
    if(GAnd.Dt.con.N<HSet.NBath)
    {
      printf("Warning: GAnd.Dt.con.N(%i)<HSet.NBath(%i)\n",GAnd.Dt.con.N,HSet.NBath);
    }
    else
    {
      GAnd.Dt.con.N=HSet.NBath;
    }
    if(GAnd.Dt.val.N<HSet.NBath)
    {
      printf("Warning: GAnd.Dt.val.N(%i)<HSet.NBath(%i)\n",GAnd.Dt.val.N,HSet.NBath);
    }
    else
    {
      GAnd.Dt.val.N=HSet.NBath;
    }
    mu=-HSet.U/2;
    SecondDoubleTriDiagonalMatrixToAndersonImpurityOperator(&HAnd, GAnd, GAnd, HSet.U, mu);
    PrintOperatorReadable(HAnd);
    /* test if psi has correct dimensions and correct if needed */
    sscanf(Psi.Name,"%s %i %i",Buffer,&PsiNBathC,&PsiNBathV);
    if(PsiNBathC!=GAnd.Dt.con.N)
    {
      printf("\nBefore diagonalizing:\nThe wavefunction has %i conduction bath sites\nThe Anderson Hamiltonian %i\n",PsiNBathC,GAnd.Dt.con.N);
      if(PsiNBathC>GAnd.Dt.con.N)
      {
        printf("Removing %i unoccupied fermions from the wave-function\n",PsiNBathC-GAnd.Dt.con.N);
        for(i=0;i<PsiNBathC-GAnd.Dt.con.N;i++)
        {
          WaveFunctionRemoveFermion(Psi, &PsiTemp, Psi.NFermion-1);
          WaveFunctionFree(&Psi);
          WaveFunctionRemoveFermion(PsiTemp, &Psi, 4+PsiNBathV+PsiNBathC-i-1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      else
      {
        printf("Adding %i unoccupied fermions to the wave-function\n",GAnd.Dt.con.N-PsiNBathC);
        for(i=0;i<GAnd.Dt.con.N-PsiNBathC;i++)
        {
          WaveFunctionAddFermion(Psi, &PsiTemp, Psi.NFermion,0);
          WaveFunctionFree(&Psi);
          WaveFunctionAddFermion(PsiTemp, &Psi, 4+PsiNBathV+PsiNBathC+i,0);
          WaveFunctionFree(&PsiTemp);
        }
      }
    }
    if(PsiNBathV!=GAnd.Dt.val.N)
    {
      printf("\nBefore diagonalizing:\nThe wavefunction has %i valence bath sites\nThe Anderson Hamiltonian %i\n",PsiNBathV,GAnd.Dt.val.N);
      if(PsiNBathV>GAnd.Dt.val.N)
      {
        printf("Removing %i occupied fermions from the wave-function\n",PsiNBathV-GAnd.Dt.val.N);
        for(i=0;i<PsiNBathV-GAnd.Dt.val.N;i++)
        {
          WaveFunctionRemoveFermion(Psi, &PsiTemp, Psi.NFermion-GAnd.Dt.con.N-1);
          WaveFunctionFree(&Psi);
          WaveFunctionRemoveFermion(PsiTemp, &Psi, 4+PsiNBathV-i-1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      else
      {
        printf("Adding %i occupied fermions to the wave-function\n",GAnd.Dt.val.N-PsiNBathV);
        for(i=0;i<GAnd.Dt.val.N-PsiNBathV;i++)
        {
          WaveFunctionAddFermion(Psi, &PsiTemp, Psi.NFermion-GAnd.Dt.val.N,1);
          WaveFunctionFree(&Psi);
          WaveFunctionAddFermion(PsiTemp, &Psi, 4+PsiNBathV+i,1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      printf("Corrected basis\n");
      PrintWaveFunctionReadableTruncated(Psi, 0.1);
    }
    sprintf(Psi.Name,"DMFT %i %i",GAnd.Dt.con.N,GAnd.Dt.val.N);
    /* create ground-state wave-funciton */
    BeginTiming("Lanczos GroundState",Time);
    LanczosGroundState(HAnd,&Psi,&HSet.zeroLancWavefunction,HSet.zeroLancEnergy,HSet.U,&Psi,0);
    LanczosGroundState(HAnd,&Psi,&HSet.zeroLancWavefunction,HSet.zeroLancEnergy,HSet.U,&Psi,0);
    EndTiming("Lanczos GroundState",Time);
    PsiHermitianOperatorPsi(HAnd, Psi, &PsiEnergy);
    PsiHermitianOperatorSqrPsi(HAnd,Psi,&PsiQuality);
    PsiQuality/=Sqr(PsiEnergy);
    PrintWaveFunctionReadableTruncated(Psi, 0.1);
    CalculatePrefactorHistogram(PsiHist, Psi, 10.0, 20);
    printf("\nHistogram of prefactors in Wave-function\n");
    for(i=0;i<19;i++)
    {
      printf("10^(-%2i) to 10^(-%2i) | %22.15E\n",i,i+1,PsiHist[i]);
    }
    printf("10^(-%2i) to 0        | %22.15E\n",i,PsiHist[19]);
    printf("\n");
    printf("Energy                  <psi|H^2|psi>/(<psi|H|psi>^2)\n");
    printf("%22.15E  %22.15E\n",PsiEnergy,PsiQuality);
    PsiRho=calloc(Psi.NFermion,sizeof(double));
    if(PsiRho==NULL){printf("calloc failed in main for PsiRho\n Fatal Error\n");fflush(stdout);return 1;}
    CalculateDensity(PsiRho,Psi);
    printf("Occupation \n");
    for(i=0;i<Psi.NFermion;i++)
    {
      printf("%5i %22.15E\n",i,PsiRho[i]);
    }
    free(PsiRho);
    /* create spin up and spin down PES and IPES */
    CreatePESdnOperatorDMFT(&OPESdn, HAnd.NFermion);
    CreatePESupOperatorDMFT(&OPESup, HAnd.NFermion);
    CreateIPESdnOperatorDMFT(&OIPESdn, HAnd.NFermion);
    CreateIPESupOperatorDMFT(&OIPESup, HAnd.NFermion);
    GclMup.N=HSet.NTriDiag;
    InitTriDiagonalMatrix(&GclMup);
    GclMdn.N=HSet.NTriDiag;
    InitTriDiagonalMatrix(&GclMdn);
    GclPup.N=HSet.NTriDiag;
    InitTriDiagonalMatrix(&GclPup);
    GclPdn.N=HSet.NTriDiag;
    InitTriDiagonalMatrix(&GclPdn);
    BeginTiming("Lanczos Tri Diagonalization PES up",Time);
    OperatorPsi(OPESup,Psi,&OPsi);
    LanczosTriDiagonalize(HAnd, OPsi, HSet.zeroLancWavefunction,&GclMup,HSet.Gamma);
    EndTiming("Lanczos Tri Diagonalization PES up",Time);
    BeginTiming("Lanczos Tri Diagonalization PES dn",Time);
    OperatorPsi(OPESdn,Psi,&OPsi);
    LanczosTriDiagonalize(HAnd, OPsi, HSet.zeroLancWavefunction,&GclMdn,HSet.Gamma);
    EndTiming("Lanczos Tri Diagonalization PES dn",Time);
    BeginTiming("Lanczos Tri Diagonalization IPES up",Time);
    OperatorPsi(OIPESup,Psi,&OPsi);
    LanczosTriDiagonalize(HAnd, OPsi, HSet.zeroLancWavefunction,&GclPup,HSet.Gamma);
    EndTiming("Lanczos Tri Diagonalization IPES up",Time);
    BeginTiming("Lanczos Tri Diagonalization IPES dn",Time);
    OperatorPsi(OIPESdn,Psi,&OPsi);
    LanczosTriDiagonalize(HAnd, OPsi, HSet.zeroLancWavefunction,&GclPdn,HSet.Gamma);
    EndTiming("Lanczos Tri Diagonalization IPES dn",Time);
    /* G math */
    printf("Start G math\n");
    ShiftTriDiagonalMatrix(&GclMdn, -PsiEnergy);
    ShiftTriDiagonalMatrix(&GclMup, -PsiEnergy);
    ShiftTriDiagonalMatrix(&GclPdn, -PsiEnergy);
    ShiftTriDiagonalMatrix(&GclPup, -PsiEnergy);
    for(i=1;i<=GclMdn.N;i++)
    {
      GclMdn.a[i]*=-1;
    }
    for(i=1;i<=GclMup.N;i++)
    {
      GclMup.a[i]*=-1;
    }
    strcpy(GclPup.Name,"Gcluster-Plus-Up");
    strcpy(GclPdn.Name,"Gcluster-Plus-Dn");
    strcpy(GclMup.Name,"Gcluster-Min-Up");
    strcpy(GclMdn.Name,"Gcluster-Min-Dn");
    printf("  Start SumGreensfunctionsTriDiagonalRepresentation\n");fflush(stdout);
    SumGreensfunctionsTriDiagonalRepresentation(GclPup, GclMup, &Gclup);
    strcpy(Gclup.Name,"Gcluster-Up");
    SumGreensfunctionsTriDiagonalRepresentation(GclPdn, GclMdn, &Gcldn);
    strcpy(Gcldn.Name,"Gcluster-Dn");
    SumGreensfunctionsTriDiagonalRepresentation(Gclup, Gcldn, &Gcl);
    strcpy(Gcl.Name,"Gcluster");
    Gcl.b[0]=1.0;
    // SymmetrizeTridiagonalMatrix(&Gcl);
    Gcl.N=Min(Gcl.N,Max(4*HSet.NTriDiag,G0.N));
    printf("  Start SelfEnergyTriDiagonalRepresentation\n");
    SelfEnergyTriDiagonalRepresentation(GBath, Gcl, &Sigma);
    Sigma.N=Min(Sigma.N,Max(4*HSet.NTriDiag,G0.N));
    strcpy(Sigma.Name,"Sigma");
    // SymmetrizeTridiagonalMatrix(&Sigma);
    Sigma.N=Min(Sigma.N,Max(4*HSet.NTriDiag,G0.N));
    printf("  Start GLocalTriDiagonalRepresentation\n");
    GLocalTriDiagonalRepresentation(G0, Sigma, &GLocal, HSet.NTriDiag);
    GLocal.N=Min(GLocal.N,Max(4*HSet.NTriDiag,G0.N));
    strcpy(GLocal.Name,"Glocal");
    // SymmetrizeTridiagonalMatrix(&GLocal);
    GLocal.N=Min(GLocal.N,Max(4*HSet.NTriDiag,G0.N));
    printf("  Start GBathNewTriDiagonalRepresentation\n");
    GBathNewTriDiagonalRepresentation(GLocal, Sigma, &GBathNew);
    GBathNew.N=Min(GBathNew.N,Max(4*HSet.NTriDiag,G0.N));
    strcpy(GBathNew.Name,"Gbath-New");
    // SymmetrizeTridiagonalMatrix(&GBathNew);
    GBathNew.N=Min(GBathNew.N,Max(4*HSet.NTriDiag,G0.N));
    /* calculate GBath for the next loop by mixing GBath with GBathNew */
    printf("Start mixing GBath with GBathNew\n");
    TriDiagonalMatrixCopy(&GBathOld, GBath);
    strcpy(GBathOld.Name,"Gbath-Old");
    FreeTriDiagonalMatrix(&GBath);
    GBathOld.b[0]=sqrt(HSet.Mix);
    GBathNew.b[0]=sqrt(1-HSet.Mix);
    SumGreensfunctionsTriDiagonalRepresentation(GBathNew, GBathOld, &GBath);
    strcpy(GBath.Name,"Gbath");
    GBathOld.b[0]=1;
    GBathNew.b[0]=1;
    //  SymmetrizeTridiagonalMatrix(&GBath);
    GBath.N=Min(GBath.N,Max(4*HSet.NTriDiag,G0.N));
    /* create spectra */
    printf("Start create spectra\n");
    Spectra[0].N       =HSet.NE;
    Spectra[0].Emin    =HSet.Emin;
    Spectra[0].Emax    =HSet.Emax;
    Spectra[0].dE      =(HSet.Emax-HSet.Emin)/HSet.NE;
    Spectra[0].EminPole=HSet.Emin;
    Spectra[0].EmaxPole=HSet.Emax;
    Spectra[0].Gamma   =HSet.Gamma;
    for(i=0;i<10;i++)
    {
      SpectrumAllocate(Spectra[0], &Spectra[i]);
    }
    TriDiagonalMatrixToSpectrum(0,GBath,&(Spectra[0]));
    TriDiagonalMatrixToSpectrum(0,GBathOld,&(Spectra[1]));
    TriDiagonalMatrixToSpectrum(0,GBathNew,&(Spectra[2]));
    SecondDoubleTriDiagonalMatrixToSpectrum(0,GAnd,&(Spectra[3]));
    TriDiagonalMatrixToSpectrum(0,G0,&(Spectra[4]));
    TriDiagonalMatrixToSpectrum(0,Sigma,&(Spectra[5]));
    TriDiagonalMatrixToSpectrum(0,GLocal,&(Spectra[6]));
    TriDiagonalMatrixToSpectrum(0,GclPup,&(Spectra[7]));
    TriDiagonalMatrixToSpectrum(0,GclMup,&(Spectra[8]));
    TriDiagonalMatrixToSpectrum(0,Gcl,&(Spectra[9]));
    /* write to disk */
    printf("Start write to dist\n");
    if(iOuter==0)
    {
      sprintf(GBathFileName,"%s001",HSet.GRootFileName);
      sprintf(PsiFileName,"%s001",HSet.WaveFunctionRootFileName);
      sprintf(SpecFileName,"%s001",HSet.SpecRootFileName);
    }
    else if(iOuter<9)
    {
      sprintf(GBathFileName,"%s00%1i",HSet.GRootFileName,iOuter+1);
      sprintf(PsiFileName,"%s00%1i",HSet.WaveFunctionRootFileName,iOuter+1);
      sprintf(SpecFileName,"%s00%1i",HSet.SpecRootFileName,iOuter+1);
    }
    else if(iOuter<99)
    {
      sprintf(GBathFileName,"%s0%2i",HSet.GRootFileName,iOuter+1);
      sprintf(PsiFileName,"%s0%2i",HSet.WaveFunctionRootFileName,iOuter+1);
      sprintf(SpecFileName,"%s0%2i",HSet.SpecRootFileName,iOuter+1);
    }
    else if(iOuter<999)
    {
      sprintf(GBathFileName,"%s%3i",HSet.GRootFileName,iOuter+1);
      sprintf(PsiFileName,"%s%3i",HSet.WaveFunctionRootFileName,iOuter+1);
      sprintf(SpecFileName,"%s%3i",HSet.SpecRootFileName,iOuter+1);
    }
    FilePrintSpectra(Spectra, 10, SpecFileName, "w");
    FilePrintPlainTriDiagonalMatrix(GBath,GBathFileName,"w");
    FilePrintWaveFunctionPlainBinary(Psi, PsiFileName, "w");
    /* free variables */
    printf("Start free variables\n");
    FreeOperator(&HAnd);
    FreeOperator(&OPESdn);
    FreeOperator(&OPESup);
    FreeOperator(&OIPESdn);
    FreeOperator(&OIPESup);
    FreeTriDiagonalMatrix(&GclMup);
    FreeTriDiagonalMatrix(&GclMdn);
    FreeTriDiagonalMatrix(&GclPup);
    FreeTriDiagonalMatrix(&GclPdn);
    FreeTriDiagonalMatrix(&Gcl);
    FreeTriDiagonalMatrix(&GLocal);
    FreeTriDiagonalMatrix(&Sigma);
    FreeTriDiagonalMatrix(&GBathNew);
    FreeTriDiagonalMatrix(&GBathOld);
    FreeDoubleTriDiagonalMatrix(&(GAnd.Dt));
    EndTiming("DMFTLoop",Time);
    PrintTiming(*Time);
  }
  return 0;
};

int oldmainDysonSecond(DMFTSettings HSet, TimeType *Time)
{
  WaveFunctionType Psi,OPsi,PsiTemp;
  TriDiagonalMatrixType G0,GBath,GclPup,GclPdn,GclMup,GclMdn,Gclup,Gcldn,Gcl,Sigma,GLocal,GBathNew,GBathOld,gPerturbation[4];
  SecondDoubleTriDiagonalMatrixType GAnd,GAndTrunk;
  OperatorType HAnd, OPESdn,OPESup,OIPESdn,OIPESup,VPerturbation[4];
  MatrixType T;
  SpectraType Spectra[10];
  char GBathFileName[255],PsiFileName[255],SpecFileName[25],Buffer[255];
  uint32_t iOuter,i,PsiNBathV,PsiNBathC;
  uint16_t oadd[1];
  double Ef=0.0,epsilonAnd=0.00001,mu;
  double PsiEnergy,PsiQuality,PsiHist[20],*PsiRho;
  /* read G0 */
  if(FileReadPlainTriDiagonalMatrixType(&G0,HSet.G0FileName)){printf("Fatal Error in main\n Failed to read G0 from file %s\n",HSet.G0FileName);return 1;};
  strcpy(G0.Name,"G0");
  /* read G_Bath_iStart and psi_iStart*/
  if(HSet.iStart==0)
  {
    sprintf(GBathFileName,"%s000",HSet.GRootFileName);
    sprintf(PsiFileName,"%s000",HSet.WaveFunctionRootFileName);
  }
  else if(HSet.iStart<10)
  {
    sprintf(GBathFileName,"%s00%1i",HSet.GRootFileName,HSet.iStart);
    sprintf(PsiFileName,"%s00%1i",HSet.WaveFunctionRootFileName,HSet.iStart);
  }
  else if(HSet.iStart<100)
  {
    sprintf(GBathFileName,"%s0%2i",HSet.GRootFileName,HSet.iStart);
    sprintf(PsiFileName,"%s0%2i",HSet.WaveFunctionRootFileName,HSet.iStart);
  }
  else
  {
    sprintf(GBathFileName,"%s%i",HSet.GRootFileName,HSet.iStart);
    sprintf(PsiFileName,"%s%i",HSet.WaveFunctionRootFileName,HSet.iStart);
  }
  if(FileReadPlainTriDiagonalMatrixType(&GBath,GBathFileName)){printf("Fatal Error in main\n Failed to read GBath from file %s !\n",GBathFileName);return 1;};
  if(file_exists(PsiFileName))
  {
    Psi.NBitsKey=HSet.NBitsKey;
    if(FileReadWaveFunctionPlainBinary(&Psi, PsiFileName)){printf("Fatal Error in main\n Failed to read Psi from file %s !\n",PsiFileName);return 1;};
    printf("\nRead the starting wave function from the file %s\n",PsiFileName);
    sscanf(Psi.Name,"%s %i %i",Buffer,&PsiNBathC,&PsiNBathV);
    printf("Wave-function \"%s\" found has %i conduction bath sites, %i valence bath sites and %i Fermions in basis\n",Psi.Name,PsiNBathC,PsiNBathV,Psi.NFermion);
    if((2*(PsiNBathC+PsiNBathV+2))!=Psi.NFermion)
    {
      printf("The total number of Fermions is incompatible witht the bath sizes\n");
      printf("This wave-function can not be used to start from\n");
      return 1;
    }
  }
  else
  {
    printf("\nWave-function file: '%s' does not excist\n Started from new function\n",PsiFileName);
    CreateWaveFunctionDMFTSecond(&Psi, HSet.NBath, HSet.NBath, HSet.NBitsKey);
  }
  strcpy(GBath.Name,"Gbath");
  /* do the DMFT loops */
  for(iOuter=HSet.iStart;iOuter<HSet.iStart+HSet.Ni;iOuter++)
  {
    BeginTiming("DMFTLoop",Time);
    printf("\n\n***********************************************************************\n");
    printf("Start of loop %i from %i to %i\n",iOuter+1,HSet.iStart,HSet.iStart+HSet.Ni);
    printf("***********************************************************************\n\n\n");
    /* Create splitted anderson matrix */
    TriDiagonalMatrixToSecondDoubleTriDiagonalMatrix(GBath, &T, &GAnd, Ef, epsilonAnd);
    /* Create GAnd and GAndTrunk as a copy pointing to the same memory */
    GAndTrunk.Ef=GAnd.Ef;
    GAndTrunk.Epsilon=GAnd.Epsilon;
    GAndTrunk.a0=GAnd.a0;
    GAndTrunk.b0=GAnd.b0;
    GAndTrunk.Dt.Ef=GAnd.Dt.Ef;
    GAndTrunk.Dt.Epsilon=GAnd.Dt.Epsilon;
    GAndTrunk.Dt.a0=GAnd.Dt.a0;
    GAndTrunk.Dt.b0=GAnd.Dt.b0;
    GAndTrunk.Dt.val.N=GAnd.Dt.val.N;
    GAndTrunk.Dt.val.a=GAnd.Dt.val.a;
    GAndTrunk.Dt.val.b=GAnd.Dt.val.b;
    GAndTrunk.Dt.con.N=GAnd.Dt.con.N;
    GAndTrunk.Dt.con.a=GAnd.Dt.con.a;
    GAndTrunk.Dt.con.b=GAnd.Dt.con.b;
    /* Create VPerturbation */
    for(i=0;i<4;i++)
    {
      InitOperator(&(VPerturbation[i]));
      VPerturbation[i].NFermion=4+2*(Min(GAnd.Dt.val.N,HSet.NBath)+Min(GAnd.Dt.con.N,HSet.NBath));
      VPerturbation[i].NBoson=0;
      VPerturbation[i].NWord=Min(VPerturbation[i].NFermion,1+(VPerturbation[i].NFermion-1)/8);
      VPerturbation[i].QComplex=0;
      InitOperatorFixedLength(&(VPerturbation[i]),1,1);
    }
    if(GAnd.Dt.val.N>HSet.NBath)
    {
      oadd[0]=CrAnMask+3+HSet.NBath;
      RealOperatorAddLader(&(VPerturbation[0]), oadd, 1, GAnd.Dt.val.b[HSet.NBath]);
      oadd[0]=CrAnMask+3+Min(GAnd.Dt.con.N,HSet.NBath)+2*HSet.NBath;
      RealOperatorAddLader(&(VPerturbation[1]), oadd, 1, GAnd.Dt.val.b[HSet.NBath]);
    }
    if(GAnd.Dt.con.N>HSet.NBath)
    {
      oadd[0]=3+Min(GAnd.Dt.val.N,HSet.NBath)+HSet.NBath;
      RealOperatorAddLader(&(VPerturbation[2]), oadd, 1, GAnd.Dt.con.b[HSet.NBath]);
      oadd[0]=3+2*Min(GAnd.Dt.val.N,HSet.NBath)+2*HSet.NBath;
      RealOperatorAddLader(&(VPerturbation[3]), oadd, 1, GAnd.Dt.con.b[HSet.NBath]);
    }
    /* Create gPerturbation */
    if(GAnd.Dt.val.N>HSet.NBath)
    {
      gPerturbation[0].N=GAnd.Dt.val.N-HSet.NBath;
      gPerturbation[1].N=GAnd.Dt.val.N-HSet.NBath;
      InitTriDiagonalMatrix(&(gPerturbation[0]));
      InitTriDiagonalMatrix(&(gPerturbation[1]));
      for(i=1;i<GAnd.Dt.val.N-HSet.NBath;i++)
      {
        gPerturbation[0].a[i]=-GAnd.Dt.val.a[i+HSet.NBath];
        gPerturbation[0].b[i]=GAnd.Dt.val.b[i+HSet.NBath];
        gPerturbation[1].a[i]=-GAnd.Dt.val.a[i+HSet.NBath];
        gPerturbation[1].b[i]=GAnd.Dt.val.b[i+HSet.NBath];
      }
      gPerturbation[0].a[GAnd.Dt.val.N-HSet.NBath]=-GAnd.Dt.val.a[GAnd.Dt.val.N];
      gPerturbation[1].a[GAnd.Dt.val.N-HSet.NBath]=-GAnd.Dt.val.a[GAnd.Dt.val.N];
    }
    else
    {
      gPerturbation[0].N=1;
      gPerturbation[1].N=1;
      InitTriDiagonalMatrix(&(gPerturbation[0]));
      InitTriDiagonalMatrix(&(gPerturbation[1]));
      gPerturbation[0].a[1]=0;
      gPerturbation[1].a[1]=0;
    }
    if(GAnd.Dt.con.N>HSet.NBath)
    {
      gPerturbation[2].N=GAnd.Dt.con.N-HSet.NBath;
      gPerturbation[3].N=GAnd.Dt.con.N-HSet.NBath;
      InitTriDiagonalMatrix(&(gPerturbation[2]));
      InitTriDiagonalMatrix(&(gPerturbation[3]));
      for(i=1;i<GAnd.Dt.con.N-HSet.NBath;i++)
      {
        gPerturbation[2].a[i]=GAnd.Dt.con.a[i+HSet.NBath];
        gPerturbation[2].b[i]=GAnd.Dt.con.b[i+HSet.NBath];
        gPerturbation[3].a[i]=GAnd.Dt.con.a[i+HSet.NBath];
        gPerturbation[3].b[i]=GAnd.Dt.con.b[i+HSet.NBath];
      }
      gPerturbation[2].a[GAnd.Dt.con.N-HSet.NBath]=GAnd.Dt.con.a[GAnd.Dt.con.N];
      gPerturbation[3].a[GAnd.Dt.con.N-HSet.NBath]=GAnd.Dt.con.a[GAnd.Dt.con.N];
    }
    else
    {
      gPerturbation[2].N=1;
      gPerturbation[3].N=1;
      InitTriDiagonalMatrix(&(gPerturbation[2]));
      InitTriDiagonalMatrix(&(gPerturbation[3]));
      gPerturbation[2].a[1]=0;
      gPerturbation[3].a[1]=0;
    };
    /* Create HAnd */
    FreeMatrix(&T);
    if(GAnd.Dt.con.N<HSet.NBath)
    {
      printf("Warning: GAnd.Dt.con.N(%i)<HSet.NBath(%i)\n",GAnd.Dt.con.N,HSet.NBath);
    }
    else
    {
      GAndTrunk.Dt.con.N=HSet.NBath;
    }
    if(GAnd.Dt.val.N<HSet.NBath)
    {
      printf("Warning: GAnd.Dt.val.N(%i)<HSet.NBath(%i)\n",GAnd.Dt.val.N,HSet.NBath);
    }
    else
    {
      GAndTrunk.Dt.val.N=HSet.NBath;
    }
    mu=-HSet.U/2;
    SecondDoubleTriDiagonalMatrixToAndersonImpurityOperator(&HAnd, GAndTrunk, GAndTrunk, HSet.U, mu);
    /* test if psi has correct dimensions and correct if needed */
    sscanf(Psi.Name,"%s %i %i",Buffer,&PsiNBathC,&PsiNBathV);
    if(PsiNBathC!=GAndTrunk.Dt.con.N)
    {
      printf("\nBefore diagonalizing:\nThe wavefunction has %i conduction bath sites\nThe Anderson Hamiltonian %i\n",PsiNBathC,GAndTrunk.Dt.con.N);
      if(PsiNBathC>GAndTrunk.Dt.con.N)
      {
        printf("Removing %i unoccupied fermions from the wave-function\n",PsiNBathC-GAndTrunk.Dt.con.N);
        for(i=0;i<PsiNBathC-GAndTrunk.Dt.con.N;i++)
        {
          WaveFunctionRemoveFermion(Psi, &PsiTemp, Psi.NFermion-1);
          WaveFunctionFree(&Psi);
          WaveFunctionRemoveFermion(PsiTemp, &Psi, 4+PsiNBathV+PsiNBathC-i-1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      else
      {
        printf("Adding %i unoccupied fermions to the wave-function\n",GAndTrunk.Dt.con.N-PsiNBathC);
        for(i=0;i<GAndTrunk.Dt.con.N-PsiNBathC;i++)
        {
          WaveFunctionAddFermion(Psi, &PsiTemp, Psi.NFermion,0);
          WaveFunctionFree(&Psi);
          WaveFunctionAddFermion(PsiTemp, &Psi, 4+PsiNBathV+PsiNBathC+i,0);
          WaveFunctionFree(&PsiTemp);
        }
      }
    }
    if(PsiNBathV!=GAndTrunk.Dt.val.N)
    {
      printf("\nBefore diagonalizing:\nThe wavefunction has %i valence bath sites\nThe Anderson Hamiltonian %i\n",PsiNBathV,GAndTrunk.Dt.val.N);
      if(PsiNBathV>GAndTrunk.Dt.val.N)
      {
        printf("Removing %i occupied fermions from the wave-function\n",PsiNBathV-GAndTrunk.Dt.val.N);
        for(i=0;i<PsiNBathV-GAndTrunk.Dt.val.N;i++)
        {
          WaveFunctionRemoveFermion(Psi, &PsiTemp, Psi.NFermion-GAndTrunk.Dt.con.N-1);
          WaveFunctionFree(&Psi);
          WaveFunctionRemoveFermion(PsiTemp, &Psi, 4+PsiNBathV-i-1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      else
      {
        printf("Adding %i occupied fermions to the wave-function\n",GAndTrunk.Dt.val.N-PsiNBathV);
        for(i=0;i<GAndTrunk.Dt.val.N-PsiNBathV;i++)
        {
          WaveFunctionAddFermion(Psi, &PsiTemp, Psi.NFermion-GAndTrunk.Dt.val.N,1);
          WaveFunctionFree(&Psi);
          WaveFunctionAddFermion(PsiTemp, &Psi, 4+PsiNBathV+i,1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      printf("Corrected basis\n");
      PrintWaveFunctionReadableTruncated(Psi, 0.1);
    }
    sprintf(Psi.Name,"DMFT %i %i",GAndTrunk.Dt.con.N,GAndTrunk.Dt.val.N);
    /* create ground-state wave-function */
    BeginTiming("Lanczos GroundState",Time);
    LanczosGroundState(HAnd,&Psi,&HSet.zeroLancWavefunction,HSet.zeroLancEnergy,HSet.U,&Psi,0);
    LanczosGroundState(HAnd,&Psi,&HSet.zeroLancWavefunction,HSet.zeroLancEnergy,HSet.U,&Psi,0);
    EndTiming("Lanczos GroundState",Time);
    PsiHermitianOperatorPsi(HAnd, Psi, &PsiEnergy);
    PsiHermitianOperatorSqrPsi(HAnd,Psi,&PsiQuality);
    PsiQuality/=Sqr(PsiEnergy);
    PrintWaveFunctionReadableTruncated(Psi, 0.1);
    CalculatePrefactorHistogram(PsiHist, Psi, 10.0, 20);
    printf("\nHistogram of prefactors in Wave-function\n");
    for(i=0;i<19;i++)
    {
      printf("10^(-%2i) to 10^(-%2i) | %22.15E\n",i,i+1,PsiHist[i]);
    }
    printf("10^(-%2i) to 0        | %22.15E\n",i,PsiHist[19]);
    printf("\n");
    printf("Energy                  <psi|H^2|psi>/(<psi|H|psi>^2)\n");
    printf("%22.15E  %22.15E\n",PsiEnergy,PsiQuality);
    PsiRho=calloc(Psi.NFermion,sizeof(double));
    if(PsiRho==NULL){printf("calloc failed in main for PsiRho\n Fatal Error\n");fflush(stdout);return 1;}
    CalculateDensity(PsiRho,Psi);
    printf("Occupation \n");
    for(i=0;i<Psi.NFermion;i++)
    {
      printf("%5i %22.15E\n",i,PsiRho[i]);
    }
    free(PsiRho);
    /*                                           */
    /* create spin up and spin down PES and IPES */
    /*                                           */
    CreatePESdnOperatorDMFT(&OPESdn, HAnd.NFermion);
    CreatePESupOperatorDMFT(&OPESup, HAnd.NFermion);
    CreateIPESdnOperatorDMFT(&OIPESdn, HAnd.NFermion);
    CreateIPESupOperatorDMFT(&OIPESup, HAnd.NFermion);
    GclMup.N=Max(HSet.NTriDiag,G0.N);
    InitTriDiagonalMatrix(&GclMup);
    GclMdn.N=Max(HSet.NTriDiag,G0.N);
    InitTriDiagonalMatrix(&GclMdn);
    GclPup.N=Max(HSet.NTriDiag,G0.N);
    InitTriDiagonalMatrix(&GclPup);
    GclPdn.N=Max(HSet.NTriDiag,G0.N);
    InitTriDiagonalMatrix(&GclPdn);
    BeginTiming("Lanczos Tri Diagonalization Dyson",Time);
    OperatorPsi(OPESup,Psi,&OPsi);
    LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation, gPerturbation, 4, &GclMup);
    EndTiming("Lanczos Tri Diagonalization Dyson",Time);
    BeginTiming("Lanczos Tri Diagonalization Dyson",Time);
    OperatorPsi(OPESdn,Psi,&OPsi);
    LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation, gPerturbation, 4, &GclMdn);
    EndTiming("Lanczos Tri Diagonalization Dyson",Time);
    BeginTiming("Lanczos Tri Diagonalization Dyson",Time);
    OperatorPsi(OIPESup,Psi,&OPsi);
    LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation, gPerturbation, 4, &GclPup);
    EndTiming("Lanczos Tri Diagonalization Dyson",Time);
    BeginTiming("Lanczos Tri Diagonalization Dyson",Time);
    OperatorPsi(OIPESdn,Psi,&OPsi);
    LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, HSet.zeroLancWavefunction, HSet.NTriDiag, VPerturbation, gPerturbation, 4, &GclPdn);
    EndTiming("Lanczos Tri Diagonalization Dyson",Time);
    /* G math */
    printf("Start G math\n");
    BeginTiming("G math",Time);
    ShiftTriDiagonalMatrix(&GclMdn, -PsiEnergy);
    ShiftTriDiagonalMatrix(&GclMup, -PsiEnergy);
    ShiftTriDiagonalMatrix(&GclPdn, -PsiEnergy);
    ShiftTriDiagonalMatrix(&GclPup, -PsiEnergy);
    for(i=1;i<=GclMdn.N;i++)
    {
      GclMdn.a[i]*=-1;
    }
    for(i=1;i<=GclMup.N;i++)
    {
      GclMup.a[i]*=-1;
    }
    strcpy(GclPup.Name,"Gcluster-Plus-Up");
    strcpy(GclPdn.Name,"Gcluster-Plus-Dn");
    strcpy(GclMup.Name,"Gcluster-Min-Up");
    strcpy(GclMdn.Name,"Gcluster-Min-Dn");
    printf("  Start SumGreensfunctionsTriDiagonalRepresentation\n");fflush(stdout);
    SumGreensfunctionsTriDiagonalRepresentation(GclPup, GclMup, &Gclup);
    strcpy(Gclup.Name,"Gcluster-Up");
    SumGreensfunctionsTriDiagonalRepresentation(GclPdn, GclMdn, &Gcldn);
    strcpy(Gcldn.Name,"Gcluster-Dn");
    SumGreensfunctionsTriDiagonalRepresentation(Gclup, Gcldn, &Gcl);
    strcpy(Gcl.Name,"Gcluster");
    Gcl.b[0]=1.0;
    printf("Gcl.N = %i ",Gcl.N);
    // SymmetrizeTridiagonalMatrix(&Gcl);
    Gcl.N=Min(Gcl.N,Max(4*HSet.NTriDiag,2*G0.N));
    printf("=> %i\n",Gcl.N);
    BeginTiming("G math Selfenergy",Time);
    printf("  Start SelfEnergyTriDiagonalRepresentation\n");
    printf("GBath.N = %i ",GBath.N);
    SelfEnergyTriDiagonalRepresentation(GBath, Gcl, &Sigma);
    printf("Sigma.N = %i ",Sigma.N);
    Sigma.N=Min(Sigma.N,Max(4*HSet.NTriDiag,2*G0.N));
    printf("=> %i ",Sigma.N);
    strcpy(Sigma.Name,"Sigma");
    // SymmetrizeTridiagonalMatrix(&Sigma);
    printf("=> %i ",Sigma.N);
    Sigma.N=Min(Sigma.N,Max(4*HSet.NTriDiag,2*G0.N));
    printf("=> %i\n",Sigma.N);
    EndTiming("G math Selfenergy",Time);
    printf("  Start GLocalTriDiagonalRepresentation\n");
    GLocalTriDiagonalRepresentation(G0, Sigma, &GLocal, Sigma.N+G0.N-1);
    printf("Glocal.N = %i ",GLocal.N);
    GLocal.N=Min(GLocal.N,Max(4*HSet.NTriDiag,2*G0.N));
    printf("=> %i ",GLocal.N);
    strcpy(GLocal.Name,"Glocal");
    // SymmetrizeTridiagonalMatrix(&GLocal);
    printf("=> %i ",GLocal.N);
    GLocal.N=Min(GLocal.N,Max(4*HSet.NTriDiag,2*G0.N));
    printf("=> %i\n",GLocal.N);
    printf("  Start GBathNewTriDiagonalRepresentation\n");
    GBathNewTriDiagonalRepresentation(GLocal, Sigma, &GBathNew);
    GBathNew.N=Min(GBathNew.N,Max(4*HSet.NTriDiag,2*G0.N));
    strcpy(GBathNew.Name,"Gbath-New");
    // SymmetrizeTridiagonalMatrix(&GBathNew);
    GBathNew.N=Min(GBathNew.N,Max(4*HSet.NTriDiag,2*G0.N));
    /* calculate GBath for the next loop by mixing GBath with GBathNew */
    printf("Start mixing GBath with GBathNew\n");
    TriDiagonalMatrixCopy(&GBathOld, GBath);
    strcpy(GBathOld.Name,"Gbath-Old");
    FreeTriDiagonalMatrix(&GBath);
    GBathOld.b[0]=sqrt(HSet.Mix);
    GBathNew.b[0]=sqrt(1-HSet.Mix);
    SumGreensfunctionsTriDiagonalRepresentation(GBathNew, GBathOld, &GBath);
    strcpy(GBath.Name,"Gbath");
    GBathOld.b[0]=1;
    GBathNew.b[0]=1;
    // SymmetrizeTridiagonalMatrix(&GBath);
    GBath.N=Min(GBath.N,Max(4*HSet.NTriDiag,2*G0.N));
    EndTiming("G math",Time);
    /* create spectra */
    printf("Start create spectra\n");
    Spectra[0].N       =HSet.NE;
    Spectra[0].Emin    =HSet.Emin;
    Spectra[0].Emax    =HSet.Emax;
    Spectra[0].dE      =(HSet.Emax-HSet.Emin)/HSet.NE;
    Spectra[0].EminPole=HSet.Emin;
    Spectra[0].EmaxPole=HSet.Emax;
    Spectra[0].Gamma   =HSet.Gamma;
    for(i=0;i<10;i++)
    {
      SpectrumAllocate(Spectra[0], &Spectra[i]);
    }
    TriDiagonalMatrixToSpectrum(0,GBath,&(Spectra[0]));              /*  3 */
    TriDiagonalMatrixToSpectrum(0,GBathOld,&(Spectra[1]));           /*  5 */
    TriDiagonalMatrixToSpectrum(0,GBathNew,&(Spectra[2]));           /*  7 */
    SecondDoubleTriDiagonalMatrixToSpectrum(0,GAnd,&(Spectra[3]));   /*  9 */
    TriDiagonalMatrixToSpectrum(0,G0,&(Spectra[4]));                 /* 11 */
    TriDiagonalMatrixToSpectrum(0,Sigma,&(Spectra[5]));              /* 13 */
    TriDiagonalMatrixToSpectrum(0,GLocal,&(Spectra[6]));             /* 15 */
    TriDiagonalMatrixToSpectrum(0,GclPup,&(Spectra[7]));             /* 17 */
    TriDiagonalMatrixToSpectrum(0,GclMup,&(Spectra[8]));             /* 19 */
    TriDiagonalMatrixToSpectrum(0,Gcl,&(Spectra[9]));                /* 21 */
    /* write to disk */
    printf("Start write to dist\n");
    if(iOuter==0)
    {
      sprintf(GBathFileName,"%s001",HSet.GRootFileName);
      sprintf(PsiFileName,"%s001",HSet.WaveFunctionRootFileName);
      sprintf(SpecFileName,"%s001",HSet.SpecRootFileName);
    }
    else if(iOuter<9)
    {
      sprintf(GBathFileName,"%s00%1i",HSet.GRootFileName,iOuter+1);
      sprintf(PsiFileName,"%s00%1i",HSet.WaveFunctionRootFileName,iOuter+1);
      sprintf(SpecFileName,"%s00%1i",HSet.SpecRootFileName,iOuter+1);
    }
    else if(iOuter<99)
    {
      sprintf(GBathFileName,"%s0%2i",HSet.GRootFileName,iOuter+1);
      sprintf(PsiFileName,"%s0%2i",HSet.WaveFunctionRootFileName,iOuter+1);
      sprintf(SpecFileName,"%s0%2i",HSet.SpecRootFileName,iOuter+1);
    }
    else if(iOuter<999)
    {
      sprintf(GBathFileName,"%s%3i",HSet.GRootFileName,iOuter+1);
      sprintf(PsiFileName,"%s%3i",HSet.WaveFunctionRootFileName,iOuter+1);
      sprintf(SpecFileName,"%s%3i",HSet.SpecRootFileName,iOuter+1);
    }
    FilePrintSpectra(Spectra, 10, SpecFileName, "w");
    FilePrintPlainTriDiagonalMatrix(GBath,GBathFileName,"w");
    FilePrintWaveFunctionPlainBinary(Psi, PsiFileName, "w");
    /* free variables */
    printf("Start free variables\n");
    for(i=0;i<4;i++)
    {
      FreeOperator(&(VPerturbation[i]));
      FreeTriDiagonalMatrix(&(gPerturbation[i]));
    }
    FreeOperator(&HAnd);
    FreeOperator(&OPESdn);
    FreeOperator(&OPESup);
    FreeOperator(&OIPESdn);
    FreeOperator(&OIPESup);
    FreeTriDiagonalMatrix(&GclMup);
    FreeTriDiagonalMatrix(&GclMdn);
    FreeTriDiagonalMatrix(&GclPup);
    FreeTriDiagonalMatrix(&GclPdn);
    FreeTriDiagonalMatrix(&Gcl);
    FreeTriDiagonalMatrix(&GLocal);
    FreeTriDiagonalMatrix(&Sigma);
    FreeTriDiagonalMatrix(&GBathNew);
    FreeTriDiagonalMatrix(&GBathOld);
    FreeDoubleTriDiagonalMatrix(&(GAnd.Dt));
    EndTiming("DMFTLoop",Time);
    
    PrintTiming(*Time);
  }
  return 0;};