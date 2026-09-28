//  Functions originated by Robert Green
//

/*  Calculate < psi | T1^+ G1^+ T2^+ G2 T2 G1 T1 | psi >
 */
uint32_t ResonantSpectra(const OperatorType OppH1, const OperatorType OppH2, const OperatorType T1, const OperatorType *T2, uint32_t NT2, const WaveFunctionType psi, TriDiagonalMatrixType **SpecTri2, SpectraType **Spec, double energy1, double energy2, ConvergenceSettingsType *settings)
{
  /*Given wavefunction psi (usually ground state, but could be excited)*/
  
  /*Calculates resonant spectra using <i|T*GTGT*GT|i> in the following steps:
   
   1. Calculate T|i>
   2. Use T|i> to transform H to tridiagonal form in Krylov basis
   3. Invert (w - H_tri + iGamma) to get G(w)
   4. Calculate G(w)T|i>
   5. Calculate T*G(w)T|i>
   6. Use the vector from step 5 to calculate RIXS spectrum using Lanczos
   */
  
  /* settings that need to be set:
   * NTri1, NTri2
   * Emin1, Emax1, NE1, Gamma1 - for automatic determination use DBL_MAX, -DBL_MAX, 0, 0
   * Emin2, Emax2, NE2, Gamma2 - for automatic determination use DBL_MAX, -DBL_MAX, 0, 0
   * epsilon
   */
  
  WaveFunctionType Tpsi,GTpsi,TGTpsi, PsiBasis[settings->NTri1];
  TriDiagonalMatrixType SpecTri1;
  MatrixType EFun;
  double EVal[settings->NTri1], EValInvRe[settings->NTri1], EValInvIm[settings->NTri1],tempDoubleRe,tempDoubleIm, EVal2[settings->NTri2], Norm;
  int i,p,k;
  double wRes;
  
  VerbosityPrintf(0x0001,"Start of ResonantSpectra using Lanczos\n");
  
  /* calcualte ground state energy */
  if(energy1==DBL_MAX)
  {
    PsiHermitianOperatorPsi(OppH1, psi, &energy1);
    WaveFunctionCalculateNormSqr(psi, &Norm);
    if(Norm > DBL_EPSILON) energy1/=Norm;
  }
  if(energy2==DBL_MAX)
  {
    PsiHermitianOperatorPsi(OppH2, psi, &energy2);
    WaveFunctionCalculateNormSqr(psi, &Norm);
    if(Norm > DBL_EPSILON) energy2/=Norm;
  }
  
  /* Calculate T|i> */
  if(settings->restrictions1.N==0)
  {
    OperatorPsi(T1,psi,&Tpsi);
  }
  else
  {
    OperatorPsiRestricted(T1, psi, &Tpsi, settings->restrictions1);
  }
  
  /* Use T|i> to transform H to tridiagonal form in Krylov basis */
  /* Using the version which gives back the Krylov basis vectors */
  SpecTri1.N = settings->NTri1;
  InitTriDiagonalMatrix(&SpecTri1);
  if(settings->restrictions1.N==0)
  {
    LanczosTriDiagonalizeKrylov(OppH1,Tpsi,settings->epsilon,&SpecTri1,PsiBasis);
  }
  else
  {
    LanczosTriDiagonalizeKrylovRestricted(OppH1, Tpsi, settings->restrictions1, settings->epsilon, &SpecTri1, PsiBasis);
  }
  
  /*Compute eigenvalues and eigenvectors, then S^-1 = QAQ^T where
   Q are the eigenvectors and A is a matrix with the eigenvalue inverses on the diagonal.
   */
  if(TriDiagonalMatrixToEigensystem(SpecTri1, &EFun, EVal))
  {
    printf("Problem with TridiagonalMatrixToEigenSystem\n");
  }
  
  /* determine Emin1, Emax1, NE1 and Gamma1 if not set */
  if(settings->Emin1>( 0.001*DBL_MAX) && (settings->Emax1<(-0.001*DBL_MAX)))
  {
    settings->Emin1=1.1*EVal[0]-0.1*EVal[settings->NTri1-1];
    settings->Emax1=1.1*EVal[settings->NTri1-1]-0.1*EVal[0];
    VerbosityPrintf(0x0002, "Automatically set energy range for intermediate states Emin=%22.15E Emax=%22.15E\n",settings->Emin1,settings->Emax1);
  }
  if(settings->Emin1>( 0.001*DBL_MAX))
  {
    settings->Emin1=1.1*EVal[0]-0.1*settings->Emax1;
    VerbosityPrintf(0x0002, "Automatically set energy range for intermediate states Emin=%22.15E Emax=%22.15E\n",settings->Emin1,settings->Emax1);
  }
  if(settings->Emax1<(-0.001*DBL_MAX))
  {
    settings->Emax1=1.1*EVal[settings->NTri1-1]-0.1*settings->Emin1;
    VerbosityPrintf(0x0002, "Automatically set energy range for intermediate states Emin=%22.15E Emax=%22.15E\n",settings->Emin1,settings->Emax1);
  }
  if((settings->Gamma1==0)&&(settings->NE1==0))
  {
    settings->NE1=1000;
    settings->Gamma1=5*(settings->Emax1-settings->Emin1)/settings->NE1;
    VerbosityPrintf(0x0002, "Automatically set number of energy points (NE=%d) and lifetime (Gamma=%22.15E) for intermediate state\n",settings->NE1,settings->Gamma1);
  }
  if(settings->NE1==0)
  {
    settings->NE1    =(uint32_t)(5*(settings->Emax1-settings->Emin1)/settings->Gamma1);
    VerbosityPrintf(0x0002, "Automatically set number of energy points (NE=%d) for intermediate state\n",settings->NE1);
  }
  if(settings->Gamma1==0)
  {
    settings->Gamma1 =(5*(settings->Emax1-settings->Emin1)/settings->NE1);
    VerbosityPrintf(0x0002, "Automatically set lifetime (Gamma=%22.15E) for intermediate state\n",settings->Gamma1);
  }
  
  /* allocate SpecTri2 */
  *SpecTri2 = (TriDiagonalMatrixType*)malloc((settings->NE1+1)*NT2*sizeof(TriDiagonalMatrixType));
  if(*SpecTri2==NULL)
  {
    printf("malloc of SpecTri2 failed in ResonantSpectra\n");
    return 1;
  }
  
  for(p = 0; p <= settings->NE1; p++)
  {
    /* New method for GT|i>:
     1. Tridiagonalize H using T|i> (outside wres loop)
     2. Calculate eigenvalues (EVal) and eigenfunctions (EFun) of Htri (outside)
     3. For each wres
     a. EValInv[i] = 1/(wres - EVal[i] + iGammaRes)
     b. GT|i> = sum_j EFun[i][j]*EValIn[j]*K[j], where K is Krylov basis vectors
     */
    
    /* Calculate 1/(E_0 + wres - EVal[i] + iGammaRes), split into real and imaginary */
    
    wRes = settings->Emin1 + (settings->Emax1-settings->Emin1)*p/((double)settings->NE1) + energy1;
    for(i = 0; i < settings->NTri1; i++)
    {
      EValInvRe[i] =  (wRes - EVal[i])       / ((wRes - EVal[i])*(wRes - EVal[i])+(settings->Gamma1/2.0)*(settings->Gamma1/2.0));
      EValInvIm[i] = -(settings->Gamma1/2.0) / ((wRes - EVal[i])*(wRes - EVal[i])+(settings->Gamma1/2.0)*(settings->Gamma1/2.0));
    }
    
    /* We are trying to calculate Q[0,:]A^(-1)Q^(-1).  Because A is diagonal, we were
     able to calculate B = Q[0,:]*I*Q^(-1) outside the loop. Now we calculate diag(A)*B */
    /* Calculate Q*A^(-1)*Q^(-1)[:,0] */
    GTpsi.NFermion = psi.NFermion;
    GTpsi.NBoson   = psi.NBoson;
    GTpsi.QComplex = psi.QComplex;
    GTpsi.NBitsKey = psi.NBitsKey;
    if(WaveFunctionInit(&GTpsi))
    {
      printf("WaveFunctionInit failed in ResonantSpectra 00");
      return 1;
    }
    for(i = 0; i < settings->NTri1; i++)
    {
      tempDoubleRe = 0;
      tempDoubleIm = 0;
      
      if(EFun.QComplex)
      {
        printf("Complex eigenfunctions for G1 in ResonantSpectra\n");
        for(k = 0; k < settings->NTri1; k++)
        {
          /*suggested correction*/
          tempDoubleRe += EValInvRe[k]*EFun.VRe[k][i]*EFun.VRe[k][0] -
          EValInvIm[k]*EFun.VIm[k][i]*EFun.VRe[k][0] +
          EValInvIm[k]*EFun.VRe[k][i]*EFun.VIm[k][0] +
          EValInvRe[k]*EFun.VIm[k][i]*EFun.VIm[k][0];
          tempDoubleIm += EValInvIm[k]*EFun.VRe[k][i]*EFun.VRe[k][0] +
          EValInvRe[k]*EFun.VIm[k][i]*EFun.VRe[k][0] -
          EValInvRe[k]*EFun.VRe[k][i]*EFun.VIm[k][0] +
          EValInvIm[k]*EFun.VIm[k][i]*EFun.VIm[k][0];
          
          /* original version:
          tempDoubleRe += EValInvRe[k]*EFun.VRe[k][i]*EFun.VRe[k][0] -
          EValInvIm[k]*EFun.VIm[k][i]*EFun.VRe[k][0] -
          EValInvIm[k]*EFun.VRe[k][i]*EFun.VIm[k][0] -
          EValInvRe[k]*EFun.VIm[k][i]*EFun.VIm[k][0];
          tempDoubleIm += EValInvIm[k]*EFun.VRe[k][i]*EFun.VRe[k][0] +
          EValInvRe[k]*EFun.VIm[k][i]*EFun.VRe[k][0] +
          EValInvRe[k]*EFun.VRe[k][i]*EFun.VIm[k][0] -
          EValInvIm[k]*EFun.VIm[k][i]*EFun.VIm[k][0];
          */
        }
      }
      else
      {
        for(k = 0; k < settings->NTri1; k++)
        {
          tempDoubleRe += EValInvRe[k]*EFun.VRe[k][i]*EFun.VRe[k][0];
          tempDoubleIm += EValInvIm[k]*EFun.VRe[k][i]*EFun.VRe[k][0];
        }
      }
      WaveFunctionAddComplexMultipliedUnionBasis(&GTpsi, PsiBasis[i],SpecTri1.b[0]*tempDoubleRe,SpecTri1.b[0]*tempDoubleIm);
    }
    for(i = 0; i < NT2; i++)
    {
      /* Calculate T*G(w)T|i> */
      if(settings->restrictions2.N==0)
      {
        OperatorPsi(T2[i],GTpsi,&TGTpsi);
      }
      else
      {
        OperatorPsiRestricted(T2[i], GTpsi, &TGTpsi, settings->restrictions2);
      }
      /* Calculate the RIXS spectrum */
      (*SpecTri2)[p + i*(settings->NE1+1)].N=settings->NTri2;
      InitTriDiagonalMatrix(&((*SpecTri2)[p + i*(settings->NE1+1)]));
      if(settings->restrictions2.N==0)
      {
        LanczosTriDiagonalize(OppH2, TGTpsi, settings->epsilon, &((*SpecTri2)[p + i*(settings->NE1+1)]), (sqrt(DBL_EPSILON)+settings->Gamma2));
      }
      else
      {
        LanczosTriDiagonalizeRestricted(OppH2, TGTpsi, settings->restrictions2, settings->epsilon, &((*SpecTri2)[p + i*(settings->NE1+1)]), (sqrt(DBL_EPSILON)+settings->Gamma2));
      }
    }
    WaveFunctionFree(&GTpsi);
  }
  /* create spectra from tridiagonal matrices */
  *Spec = (SpectraType*)malloc(NT2 * (settings->NE1+1) * sizeof(SpectraType));
  for(i=0; i<NT2 * (settings->NE1+1);i++)
  {
    /* determine min and max poles */
    TriDiagonalMatrixToEigenvalues((*SpecTri2)[i], EVal2);
    (*Spec)[i].EminPole = EVal2[0]-energy2;
    (*Spec)[i].EmaxPole = EVal2[settings->NTri2-1]-energy2;
  }
  /* determine the energy range of the spectrum if not give */
  EVal2[0]=(*Spec)[0].EminPole;
  EVal2[1]=(*Spec)[0].EmaxPole;
  for(i=1; i<NT2 * (settings->NE1+1);i++)
  {
    if((*Spec)[i].EminPole<EVal2[0]) EVal2[0]=(*Spec)[i].EminPole;
    if((*Spec)[i].EmaxPole>EVal2[1]) EVal2[1]=(*Spec)[i].EmaxPole;
  }
  if((settings->Emin2>( 0.001*DBL_MAX)) || (settings->Emax2<(-0.001*DBL_MAX)))
  {
    settings->Emin2 = 1.1 * EVal2[0] - 0.1*EVal2[1];
    settings->Emax2 = 1.1 * EVal2[1] - 0.1*EVal2[0];
    VerbosityPrintf(0x0002, "Automatically set energy range for final states Emin=%22.15E Emax=%22.15E\n",settings->Emin2,settings->Emax2);
  }
  if(settings->Emin2>( 0.001*DBL_MAX))
  {
    settings->Emin2=1.1*EVal2[0]-0.1*settings->Emax2;
    VerbosityPrintf(0x0002, "Automatically set energy range for final states Emin=%22.15E Emax=%22.15E\n",settings->Emin2,settings->Emax2);
  }
  if(settings->Emax2<(-0.001*DBL_MAX))
  {
    settings->Emax2=1.1*EVal2[1]-0.1*settings->Emin2;
    VerbosityPrintf(0x0002, "Automatically set energy range for final states Emin=%22.15E Emax=%22.15E\n",settings->Emin2,settings->Emax2);
  }
  if((settings->Gamma2==0)&&(settings->NE2==0))
  {
    settings->NE2=1000;
    settings->Gamma2=5*(settings->Emax2-settings->Emin2)/settings->NE2;
    VerbosityPrintf(0x0002, "Automatically set number of energy points (NE=%d) and lifetime (Gamma=%22.15E) for final state\n",settings->NE2,settings->Gamma2);
  }
  if(settings->NE2==0)
  {
    settings->NE2    =(uint32_t)(5*(settings->Emax2-settings->Emin2)/settings->Gamma2);
    VerbosityPrintf(0x0002, "Automatically set number of energy points (NE=%d) for final state\n",settings->NE2);
  }
  if(settings->Gamma2==0)
  {
    settings->Gamma2 =(5*(settings->Emax2-settings->Emin2)/settings->NE2);
    VerbosityPrintf(0x0002, "Automatically set lifetime (Gamma=%22.15E) for final state\n",settings->Gamma2);
  }
  for(i=0;i<NT2*(settings->NE1+1);i++)
  {
    (*Spec)[i].N    =settings->NE2;
    (*Spec)[i].Emin =settings->Emin2;
    (*Spec)[i].Emax =settings->Emax2;
    (*Spec)[i].dE   =(settings->Emax2-settings->Emin2)/settings->NE2;
    (*Spec)[i].Gamma=settings->Gamma2;
    (*Spec)[i].ResonantEnergy = settings->Emin1 + ((double)(i % (settings->NE1+1))/(settings->NE1)) * (settings->Emax1-settings->Emin1);
    sprintf((*Spec)[i].Name, "Spectrum wres %i/%i T2 %i/%i",(int)(i % (settings->NE1+1)),settings->NE1,(int)(i/(settings->NE1+1)),NT2);
    SpectrumAllocate((*Spec)[i], &(*Spec)[i]);
    TriDiagonalMatrixToSpectrum(energy2, (*SpecTri2)[i], &(*Spec)[i]);
  }
  return 0;
}

/*  Calculate < psi | T1^+ G1^+ T2^+ G2 T2 G1 T1 | psi >
 *  psi, T1 and T2 are vectors of Length Np, NT1 and NT2
 *  psi and T1 are contract to a single vector of length Np*NT1
 *  returned are (Np*NT1) * NT2 spectra
 */
uint32_t ResonantSpectraDense(WaveFunctionType * H1psiBasis, WaveFunctionType * H2psiBasis,const OperatorType OppH1, const OperatorType OppH2, const OperatorType *T1, uint32_t NT1, const OperatorType *T2, uint32_t NT2, const WaveFunctionType *psi, uint32_t Npsi, SpectraType **Spec, double energy1, double energy2, ConvergenceSettingsType *settings)
{
  /*Given wavefunction psi (usually ground state, but could be excited)*/
  
  /*Calculates resonant spectra using <i|T*GTGT*GT|i> in the following steps:
   
   1. Calculate MT1i  = T1|i> for all elments of T1 and |i> on a basis of H1psiBasis as a dense matrix
      The dimensions of MT1i are Ni = Npsi*NT1 Nj=H1psiBasis.N
      the order is {T1[0] * psi[0], T1[1] * psi[0], ..., T1[NT1-1] * psi[0], ..., ..., T1[0] * psi[Npsi-1], T1[1] * psi[Npsi-1], ..., T1[NT1-1] * psi[Npsi-1]}
   2. Calcualte MfT2m[k] = a dense matrix rerpesenting T2 on a basis of H1psiBasis and H2psiBasis for each k from 0 to NT2
      The dimensions of MfT2m[k] are Ni = N2psiBasis.N by Nj = H1psiBasis.N
   3. Calculate Hm the dense matrix representing H1 on a basis of H1psiBasis
      The dimensions of Hm are Ni = H1psiBasis.N by Nj = H1psiBasis.N
   4. Calculate Hf the dense matrix representing H2 on a basis of H2psiBasis
      The dimensions of Hf are Ni = H2psiBasis.N by Nj = H2psiBasis.N
   
   5. calculate Valm and Funm as the eigenvalues and functions of Hm
      The dimensions of Funm are Ni = H1psiBasis.N by Nj = H1psiBasis.N
   6. calculate A = Funm^* * MT1i^T
      The dimensions of A are Ni = H1psiBasis.N by Nj = Npis*NT1
   7. calculate Valf and Funf as the eigenvalues and functions of Hf
      The dimensions of Funf are Ni = H2psiBasis.N by Nj = H2psiBasis.N
   8. calcualte B[k] = Funf^* * MfT2m[k] * Funm^T
      The dimensions of B[k] are Ni = H2psiBasis.N by Nj = H1psiBasis.N
   
   For each inicdent photon energy wi[l] do (l from 0 to Nwi)
  
     9. Calcualte the matrix C as a copy of A with each element scaled: Cij = Aij * 1/(wi[l] - Em[i] + Ei[j] + i Gm)
        The dimensions of C are Ni = H1psiBasis.N by Nj = Npsi*NT1
   
     For each decay operator T2[k] do
   
       10. Calculate the matrix D = B[k] * C
           The dimensions of D are Ni = H2psiBasis.N by Nj = Npsi*NT1
           Dij contains the scattering amplitude from initial T1psi[j] to final state i with resonant energy wi[l] and decay operator T2[k]
       11. Square all elements of matrix D to get the scattering intensity
           Dij = Dij^* * Dij
   
       For each transfer energy wr[m] do (m from 0 to Nwr)
       
         12. Spec[T1psi[j],T2[k],wi[l]][wr[m]] = sum_i 1 / (wr[m] - Ef[i] + Ei[j] + i Gf)
             The index of the spectrum is:  l + Nwr * ( k + NT2 * ( jT1 + NT1 * ( jpsi ) ) )  where j = ( jT1 + NT1 * ( jpsi ) )
   
   */
  CompactMatrixType MT1i, *MfT2m, Hm, Hf, Funm, Funf, A, B[NT2], C, D;
  WaveFunctionType Opsi;
  double *Valm, *Valf, energy1List[Npsi], energy2List[Npsi], wRes, wLoss, Norm, energy1ListMin, energy1ListMax, energy2ListMin, energy2ListMax;
  MfT2m= (CompactMatrixType*)malloc(NT2 * sizeof(CompactMatrixType));
  Valm = (double*)malloc(H1psiBasis->N*sizeof(double));
  Valf = (double*)malloc(H2psiBasis->N*sizeof(double));
  uint32_t i,j,k,l,m,n,o;
  Complex16 Lorentz;
  
  VerbosityPrintf(0x0001,"Start of ResonantSpectra using Dense Methods\n");
  
  //
  *Spec = (SpectraType*)malloc(NT1 * Npsi * NT2 * (settings->NE1+1) * sizeof(SpectraType));
  if(Spec==NULL)
  {
    printf("ERROR: malloc failed in ResonantSpectraDense\n");
    return 1;
  }
  // Calculate the energy of the initial states used as a shift of the spectra
  if(energy1==DBL_MAX)
  {
    for(i=0;i<Npsi;i++)
    {
      PsiHermitianOperatorPsi(OppH1, psi[i], &energy1List[i]);
      WaveFunctionCalculateNormSqr(psi[i], &Norm);
      if(Norm > DBL_EPSILON) energy1List[i]/=Norm;
    }
  }
  else
  {
    for(i=0;i<Npsi;i++)
    {
      energy1List[i]=energy1;
    }
  }
  if(energy2==DBL_MAX)
  {
    for(i=0;i<Npsi;i++)
    {
      PsiHermitianOperatorPsi(OppH2, psi[i], &energy2List[i]);
      WaveFunctionCalculateNormSqr(psi[i], &Norm);
      if(Norm > DBL_EPSILON) energy2List[i]/=Norm;
    }
  }
  else
  {
    for(i=0;i<Npsi;i++)
    {
      energy2List[i]=energy2;
    }
  }
  // get extrema of shifts
  energy1ListMin=energy1List[0];
  energy1ListMax=energy1List[0];
  energy2ListMin=energy2List[0];
  energy2ListMax=energy2List[0];
  for(i=1;i<Npsi;i++)
  {
    if(energy1List[i]<energy1ListMin) energy1ListMin=energy1List[i];
    if(energy1List[i]>energy1ListMax) energy1ListMax=energy1List[i];
    if(energy2List[i]<energy2ListMin) energy2ListMin=energy2List[i];
    if(energy2List[i]>energy2ListMax) energy2ListMax=energy2List[i];
  }
  //1. Calculate MT1i  = T1|i> for all elments of T1 and |i> on a basis of H1psiBasis as a dense matrix
  //   The dimensions of MT1i are Ni = Npsi*NT1 Nj=H1psiBasis.N
  //   the order is {T1[0] * psi[0], T1[1] * psi[0], ..., T1[NT1-1] * psi[0], ..., ..., T1[0] * psi[Npsi-1], T1[1] * psi[Npsi-1], ..., T1[NT1-1] * psi[Npsi-1]}
  MT1i.Ni = Npsi * NT1;
  MT1i.Nj = H1psiBasis->N;
  MT1i.QComplex = H1psiBasis->QComplex;
  InitCompactMatrix(&MT1i);
  k=0;
  for(i=0;i<Npsi;i++)
  {
    for(j=0;j<NT1;j++)
    {
      if(settings->restrictions1.N==0)
      {
        OperatorPsi(T1[j], psi[i], &Opsi);
      }
      else
      {
        OperatorPsiRestricted(T1[j], psi[i], &Opsi, settings->restrictions1);
      }
      WaveFunctionMultiply(H1psiBasis, 0, 0);
      WaveFunctionAddConserveBasis(H1psiBasis, Opsi);
      WaveFunctionFree(&Opsi);
      if(H1psiBasis->QComplex==0)
      {
        for(l=1;l<=H1psiBasis->N;l++)
        {
          MT1i.VRe[k] = H1psiBasis->ValueRe[BlokElement(l)];
          k++;
        }
      }
      else
      {
        for(l=1;l<=H1psiBasis->N;l++)
        {
          MT1i.V[k].Re = H1psiBasis->ValueRe[BlokElement(l)];
          MT1i.V[k].Im = H1psiBasis->ValueIm[BlokElement(l)];
          k++;
        }
      }
    }
  }
  //2. Calcualte MfT2m[k] = a dense matrix rerpesenting T2 on a basis of H1psiBasis and H2psiBasis for each k from 0 to NT2
  //   The dimensions of MfT2m[k] are Ni = H2psiBasis.N by Nj = H1psiBasis.N
  for(k=0;k<NT2;k++)
  {
    OperatorToDeterminantMatrix2(&(T2[k]), H2psiBasis, H1psiBasis, &(MfT2m[k]));
  }
  //3. Calculate Hm the dense matrix representing H1 on a basis of H1psiBasis
  //   The dimensions of Hm are Ni = H1psiBasis.N by Nj = H1psiBasis.N
  OperatorToDeterminantMatrix(&OppH1, H1psiBasis, &Hm);
  //4. Calculate Hf the dense matrix representing H2 on a basis of H2psiBasis
  //   The dimensions of Hf are Ni = H2psiBasis.N by Nj = H2psiBasis.N
  OperatorToDeterminantMatrix(&OppH2, H2psiBasis, &Hf);
  //5. calculate Valm and Funm as the eigenvalues and functions of Hm
  //   The dimensions of Funm are Ni = H1psiBasis.N by Nj = H1psiBasis.N
  CompactMatrixEigensystem(&Hm, &Funm, Valm);
  //6. calculate A = Funm^* * MT1i^T
  //   The dimensions of A are Ni = H1psiBasis.N by Nj = Npis*NT1
  CompactMatrixMultiply('*', 'T', &Funm, &MT1i, &A);
  //7. calculate Valf and Funf as the eigenvalues and functions of Hf
  //   The dimensions of Funf are Ni = H2psiBasis.N by Nj = H2psiBasis.N
  CompactMatrixEigensystem(&Hf, &Funf, Valf);
  //8. calcualte B[k] = Funf^* * MfT2m[k] * Funm^T
  //   The dimensions of B[k] are Ni = H2psiBasis.N by Nj = H1psiBasis.N
  for(k=0;k<NT2;k++)
  {
    CompactMatrixMultiply('N', 'T', &Funm, &(MfT2m[k]), &C);
    CompactMatrixMultiply('*', 'T', &Funf, &C, &(B[k]));
    FreeCompactMatrix(&C);
  }
  // free eigenfunctions
  FreeCompactMatrix(&Funm);
  FreeCompactMatrix(&Funf);
  //
  /* determine Emin1, Emax1, NE1 and Gamma1 if not set */
  if(settings->Emin1>( 0.001*DBL_MAX) && (settings->Emax1<(-0.001*DBL_MAX)))
  {
    settings->Emin1=1.1*Valm[0]-0.1*Valm[H1psiBasis->N-1]-energy1ListMax;
    settings->Emax1=1.1*Valm[H1psiBasis->N-1]-0.1*Valm[0]-energy1ListMin;
    VerbosityPrintf(0x0002, "Automatically set energy range for intermediate states Emin=%22.15E Emax=%22.15E\n",settings->Emin1,settings->Emax1);
  }
  if(settings->Emin1>( 0.001*DBL_MAX))
  {
    settings->Emin1=1.1*Valm[0]-0.1*settings->Emax1-energy1ListMax;
    VerbosityPrintf(0x0002, "Automatically set energy range for intermediate states Emin=%22.15E Emax=%22.15E\n",settings->Emin1,settings->Emax1);
  }
  if(settings->Emax1<(-0.001*DBL_MAX))
  {
    settings->Emax1=1.1*Valm[H1psiBasis->N-1]-0.1*settings->Emin1-energy1ListMin;
    VerbosityPrintf(0x0002, "Automatically set energy range for intermediate states Emin=%22.15E Emax=%22.15E\n",settings->Emin1,settings->Emax1);
  }
  if((settings->Gamma1==0)&&(settings->NE1==0))
  {
    settings->NE1=1000;
    settings->Gamma1=5*(settings->Emax1-settings->Emin1)/settings->NE1;
    VerbosityPrintf(0x0002, "Automatically set number of energy points (NE=%d) and lifetime (Gamma=%22.15E) for intermediate state\n",settings->NE1,settings->Gamma1);
  }
  if(settings->NE1==0)
  {
    settings->NE1    =(uint32_t)(5*(settings->Emax1-settings->Emin1)/settings->Gamma1);
    VerbosityPrintf(0x0002, "Automatically set number of energy points (NE=%d) for intermediate state\n",settings->NE1);
  }
  if(settings->Gamma1==0)
  {
    settings->Gamma1 =(5*(settings->Emax1-settings->Emin1)/settings->NE1);
    VerbosityPrintf(0x0002, "Automatically set lifetime (Gamma=%22.15E) for intermediate state\n",settings->Gamma1);
  }
  /* determine Emin2, Emax2, NE2 and Gamma2 if not set */
  if(settings->Emin2>( 0.001*DBL_MAX) && (settings->Emax2<(-0.001*DBL_MAX)))
  {
    settings->Emin2=1.1*Valf[0]-0.1*Valf[H2psiBasis->N-1]-energy2ListMax;
    settings->Emax2=1.1*Valf[H2psiBasis->N-1]-0.1*Valf[0]-energy2ListMin;
    VerbosityPrintf(0x0002, "Automatically set energy range for final states Emin=%22.15E Emax=%22.15E\n",settings->Emin2,settings->Emax2);
  }
  if(settings->Emin2>( 0.001*DBL_MAX))
  {
    settings->Emin2=1.1*Valf[0]-0.1*settings->Emax2-energy2ListMax;
    VerbosityPrintf(0x0002, "Automatically set energy range for final states Emin=%22.15E Emax=%22.15E\n",settings->Emin2,settings->Emax2);
  }
  if(settings->Emax2<(-0.001*DBL_MAX))
  {
    settings->Emax2=1.1*Valf[H2psiBasis->N-1]-0.1*settings->Emin2-energy2ListMin;
    VerbosityPrintf(0x0002, "Automatically set energy range for final states Emin=%22.15E Emax=%22.15E\n",settings->Emin2,settings->Emax2);
  }
  if((settings->Gamma2==0)&&(settings->NE2==0))
  {
    settings->NE2=1000;
    settings->Gamma2=5*(settings->Emax2-settings->Emin2)/settings->NE2;
    VerbosityPrintf(0x0002, "Automatically set number of energy points (NE=%d) and lifetime (Gamma=%22.15E) for final state\n",settings->NE2,settings->Gamma2);
  }
  if(settings->NE2==0)
  {
    settings->NE2    =(uint32_t)(5*(settings->Emax2-settings->Emin2)/settings->Gamma2);
    VerbosityPrintf(0x0002, "Automatically set number of energy points (NE=%d) for final state\n",settings->NE2);
  }
  if(settings->Gamma2==0)
  {
    settings->Gamma2 =(5*(settings->Emax2-settings->Emin2)/settings->NE2);
    VerbosityPrintf(0x0002, "Automatically set lifetime (Gamma=%22.15E) for final state\n",settings->Gamma2);
  }
  //   For each inicdent photon energy wi[l] do (l from 0 to Nwi)
  C.Ni=A.Ni;
  C.Nj=A.Nj;
  C.QComplex=1;
  InitCompactMatrix(&C);
  for(l = 0; l <= settings->NE1; l++)
  {
    //9. Calcualte the matrix C as a copy of A with each element scaled: Cij = Aij * 1/(wi[l] - Em[i] + Ei[j] + i Gm)
    //   The dimensions of C are Ni = H1psiBasis.N by Nj = Npsi*NT1
    /* Calculate 1/(E_0 + wres - EVal[i] + iGammaRes), split into real and imaginary */
    wRes = settings->Emin1 + (settings->Emax1-settings->Emin1)*l/((double)settings->NE1);


    o=0;
    if(A.QComplex==0)
    {
      for(i=0;i<H1psiBasis->N;i++)
      {
        for(j=0;j<Npsi;j++)
        {
          Lorentz.Re = 1/(Sqr(wRes - Valm[i] + energy1List[j]) + Sqr(settings->Gamma1/2));
          Lorentz.Im = -Lorentz.Re * settings->Gamma1/2;
          Lorentz.Re *= wRes - Valm[i] + energy1List[j];
          for(k=0;k<NT1;k++)
          {
            C.V[o].Re = A.VRe[o] * Lorentz.Re;
            C.V[o].Im = A.VRe[o] * Lorentz.Im;
            o++;
          }
        }
      }
    }
    else
    {
      for(i=0;i<H1psiBasis->N;i++)
      {
        for(j=0;j<Npsi;j++)
        {
          Lorentz.Re = 1/(Sqr(wRes - Valm[i] + energy1List[j]) + Sqr(settings->Gamma1*0.5));
          Lorentz.Im = -Lorentz.Re * settings->Gamma1 * 0.5;
          Lorentz.Re *= (wRes - Valm[i] + energy1List[j]);
          for(k=0;k<NT1;k++)
          {
            C.V[o].Re = A.V[o].Re * Lorentz.Re - A.V[o].Im * Lorentz.Im;
            C.V[o].Im = A.V[o].Re * Lorentz.Im + A.V[o].Im * Lorentz.Re;
            o++;
          }
        }
      }
    }
    //   For each decay operator T2[k] do
    for(k=0;k<NT2;k++)
    {
      //10. Calculate the matrix D = B[k] * C
      //    The dimensions of D are Ni = H2psiBasis.N by Nj = Npsi*NT1
      //    Dij contains the scattering amplitude from initial T1psi[j] to final state i with resonant energy wi[l] and decay operator T2[k]
      CompactMatrixMultiply('N', 'N', &(B[k]), &C, &D);
      //11. Square all elements of matrix D to get the scattering intensity
      //    Dij = Dij^* * Dij
      o=0;
      for(i=0;i<D.Ni;i++)
      {
        for(j=0;j<D.Nj;j++)
        {
          D.V[o].Re*=D.V[o].Re;
          D.V[o].Re+=Sqr(D.V[o].Im);
          o++;
        }
      }
      // loop over all T1psi
      for(j=0;j<NT1*Npsi;j++)
      {
        n = l + (settings->NE1+1) * (k + NT2 * j);
        // allocate the spectra
        (*Spec)[n].N    =settings->NE2;
        (*Spec)[n].Emin =settings->Emin2;
        (*Spec)[n].Emax =settings->Emax2;
        (*Spec)[n].dE   =(settings->Emax2-settings->Emin2)/settings->NE2;
        (*Spec)[n].Gamma=settings->Gamma2;
        (*Spec)[n].ResonantEnergy = settings->Emin1 + ((double)(l)/(settings->NE1)) * (settings->Emax1-settings->Emin1);
        sprintf((*Spec)[n].Name, "Spectrum wres %i/%i T2 %i/%i T1 %i/%i psi %i/%i",l,settings->NE1,k,NT2,j%NT1,NT1,(int)(j/NT1),Npsi);
        SpectrumAllocate((*Spec)[n], &(*Spec)[n]);
        //    calculate the integral
        for(i=0;i<H2psiBasis->N;i++)
        {
          (*Spec)[n].b[0]  += D.V[i*D.Nj+j].Re;
        }
        //    For each transfer energy wr[m] do (m from 0 to Nwr)
        for(m=0;m<=settings->NE2;m++)
        {
          wLoss = settings->Emin2 + ((double)m/(double)(settings->NE2))*(settings->Emax2 - settings->Emin2);
          //12. Spec[T1psi[j],T2[k],wi[l]][wr[m]] = sum_i Dij / (wr[m] - Ef[i] + Ei[j] + i Gf)
          //    The index of the spectrum is:  l + Nwr * ( k + NT2 * ( jT1 + NT1 * ( jpsi ) ) )  where j = ( jT1 + NT1 * ( jpsi ) )
          for(i=0;i<H2psiBasis->N;i++)
          {
            (*Spec)[n].Im[m] -= D.V[i*D.Nj+j].Re * (settings->Gamma2*0.5)                 / (Sqr(wLoss - Valf[i] + energy2List[j/NT1]) + Sqr(settings->Gamma2*0.5));
            (*Spec)[n].Re[m] += D.V[i*D.Nj+j].Re * (wLoss - Valf[i] + energy2List[j/NT1]) / (Sqr(wLoss - Valf[i] + energy2List[j/NT1]) + Sqr(settings->Gamma2*0.5));
          } // end loop over all final states
        } // end loop over energy loss photon energies index m
        (*Spec)[n].b[0] = sqrt((*Spec)[n].b[0]);
        if(H2psiBasis->N==0)
        {
          (*Spec)[n].EmaxPole=0;
          (*Spec)[n].EminPole=0;
        }
        else
        {
          (*Spec)[n].EmaxPole = Valf[H2psiBasis->N-1] -  energy2List[j/NT1];
          (*Spec)[n].EminPole = Valf[0] -  energy2List[j/NT1];
        }
      } // end loop over T1psi indexed by j
      FreeCompactMatrix(&D);
    } // end loop over decay operators index k
  } // end loop over incomming photon energies index l
  FreeCompactMatrix(&A);
  for(k=0;k<NT2;k++)
  {
    FreeCompactMatrix(&(B[k]));
  }
  FreeCompactMatrix(&C);
  free(Valm);
  free(Valf);
  free(MfT2m);
  return 0;
}

// TODO:
uint32_t ResonantSpectraDenseTensor(WaveFunctionType * H1psiBasis, WaveFunctionType * H2psiBasis,const OperatorType OppH1, const OperatorType OppH2, const OperatorType *T1, uint32_t NT1, const OperatorType *T2, uint32_t NT2, const WaveFunctionType *psi, uint32_t Npsi, SpectraType **Spec, double energy1, double energy2, ConvergenceSettingsType *settings)
{
 /*Given wavefunction psi (usually ground state, but could be excited)*/
 
 /*Calculates resonant spectra using <i|T*GTGT*GT|i> in the following steps:
  
  1. Calculate MT1i  = T1|i> for all elments of T1 and |i> on a basis of H1psiBasis as a dense matrix
     The dimensions of MT1i are Ni = Npsi*NT1 Nj=H1psiBasis.N
     the order is {T1[0] * psi[0], T1[1] * psi[0], ..., T1[NT1-1] * psi[0], ..., ..., T1[0] * psi[Npsi-1], T1[1] * psi[Npsi-1], ..., T1[NT1-1] * psi[Npsi-1]}
  2. Calcualte MfT2m[k] = a dense matrix rerpesenting T2 on a basis of H1psiBasis and H2psiBasis for each k from 0 to NT2
     The dimensions of MfT2m[k] are Ni = N2psiBasis.N by Nj = H1psiBasis.N
  3. Calculate Hm the dense matrix representing H1 on a basis of H1psiBasis
     The dimensions of Hm are Ni = H1psiBasis.N by Nj = H1psiBasis.N
  4. Calculate Hf the dense matrix representing H2 on a basis of H2psiBasis
     The dimensions of Hf are Ni = H2psiBasis.N by Nj = H2psiBasis.N
  
  5. calculate Valm and Funm as the eigenvalues and functions of Hm
     The dimensions of Funm are Ni = H1psiBasis.N by Nj = H1psiBasis.N
  6. calculate A = Funm^* * MT1i^T
     The dimensions of A are Ni = H1psiBasis.N by Nj = Npsi*NT1
  7. calculate Valf and Funf as the eigenvalues and functions of Hf
     The dimensions of Funf are Ni = H2psiBasis.N by Nj = H2psiBasis.N
  8. calcualte B[k] = Funf^* * MfT2m[k] * Funm^T
     The dimensions of B[k] are Ni = H2psiBasis.N by Nj = H1psiBasis.N
  
  For each inicdent photon energy wi[l] do (l from 0 to Nwi)
 
    9. Calcualte the matrix C as a copy of A with each element scaled: Cij = Aij * 1/(wi[l] - Em[i] + Ei[j] + i Gm)
       The dimensions of C are Ni = H1psiBasis.N by Nj = Npsi*NT1
  
    //Here is now where The dense Tensor version differs from the non tensor one. we can apply Create spectra in tensor mode. 
  
    For each decay operator T2[k] do
  
      //10. Calculate the matrix R[k] = B[k] * C
      //    The dimensions of R[k] are Ni = H2psiBasis.N by Nj = Npsi*NT1
      //    R[k]_ij contains the scattering amplitude from initial T1psi[j] to final state i with resonant energy wi[l] and decay operator T2[k]
  
    For each transfer energy wr[m] do (m from 0 to Nwr)

    //11. Compute the spectra for each combinations of Psi and T1,T2 
    //    each tensor element is spec[t]= <Psi[s]|T1_dag[m] G1_dag(w1) T2_dag[n] * G2(w2) * T2[j] * G1(w1) * T1[i]|Psi[r]>
    //    where t= l + (settings->NE1+1) * (n + NT2 *(m+NT1*(j+NT2*(i+NT1*(s+Npsi*r)))));   


      
  
  */
 CompactMatrixType MT1i, *MfT2m, Hm, Hf, Funm, Funf, A, B[NT2], C, R[NT2];
 WaveFunctionType Opsi;
 double *Valm, *Valf, energy1List[Npsi], energy2List[Npsi], wRes, wLoss, Norm, energy1ListMin, energy1ListMax, energy2ListMin, energy2ListMax;
 MfT2m= (CompactMatrixType*)malloc(NT2 * sizeof(CompactMatrixType));
 Valm = (double*)malloc(H1psiBasis->N*sizeof(double));
 Valf = (double*)malloc(H2psiBasis->N*sizeof(double));
 uint32_t i,j,k,l,m,n,o,r,s,t,a,b,kRr,kRs;
 Complex16 Lorentz;
 
 VerbosityPrintf(0x0001,"Start of ResonantSpectra Tensor using Dense Methods\n");
 
 //
 *Spec = (SpectraType*)malloc(Sqr(NT1 * Npsi * NT2) * (settings->NE1+1) * sizeof(SpectraType));
 if(Spec==NULL)
 {
   printf("ERROR: malloc failed in ResonantSpectraDense\n");
   return 1;
 }
 // Calculate the energy of the initial states used as a shift of the spectra
 if(energy1==DBL_MAX)
 {
   for(i=0;i<Npsi;i++)
   {
     PsiHermitianOperatorPsi(OppH1, psi[i], &energy1List[i]);
     WaveFunctionCalculateNormSqr(psi[i], &Norm);
     if(Norm > DBL_EPSILON) energy1List[i]/=Norm;
   }
 }
 else
 {
   for(i=0;i<Npsi;i++)
   {
     energy1List[i]=energy1;
   }
 }
 if(energy2==DBL_MAX)
 {
   for(i=0;i<Npsi;i++)
   {
     PsiHermitianOperatorPsi(OppH2, psi[i], &energy2List[i]);
     WaveFunctionCalculateNormSqr(psi[i], &Norm);
     if(Norm > DBL_EPSILON) energy2List[i]/=Norm;
   }
 }
 else
 {
   for(i=0;i<Npsi;i++)
   {
     energy2List[i]=energy2;
   }
 }
 // get extrema of shifts
 energy1ListMin=energy1List[0];
 energy1ListMax=energy1List[0];
 energy2ListMin=energy2List[0];
 energy2ListMax=energy2List[0];
 for(i=1;i<Npsi;i++)
 {
   if(energy1List[i]<energy1ListMin) energy1ListMin=energy1List[i];
   if(energy1List[i]>energy1ListMax) energy1ListMax=energy1List[i];
   if(energy2List[i]<energy2ListMin) energy2ListMin=energy2List[i];
   if(energy2List[i]>energy2ListMax) energy2ListMax=energy2List[i];
 }
 //1. Calculate MT1i  = T1|i> for all elments of T1 and |i> on a basis of H1psiBasis as a dense matrix
 //   The dimensions of MT1i are Ni = Npsi*NT1 Nj=H1psiBasis.N
 //   the order is {T1[0] * psi[0], T1[1] * psi[0], ..., T1[NT1-1] * psi[0], ..., ..., T1[0] * psi[Npsi-1], T1[1] * psi[Npsi-1], ..., T1[NT1-1] * psi[Npsi-1]}
 MT1i.Ni = Npsi * NT1;
 MT1i.Nj = H1psiBasis->N;
 MT1i.QComplex = H1psiBasis->QComplex;
 InitCompactMatrix(&MT1i);
 k=0;
 for(i=0;i<Npsi;i++)
 {
   for(j=0;j<NT1;j++)
   {
     if(settings->restrictions1.N==0)
     {
       OperatorPsi(T1[j], psi[i], &Opsi);
     }
     else
     {
       OperatorPsiRestricted(T1[j], psi[i], &Opsi, settings->restrictions1);
     }
     WaveFunctionMultiply(H1psiBasis, 0, 0);
     //why add????
     WaveFunctionAddConserveBasis(H1psiBasis, Opsi);
     WaveFunctionFree(&Opsi);
     if(H1psiBasis->QComplex==0)
     {
       for(l=1;l<=H1psiBasis->N;l++)
       {
         MT1i.VRe[k] = H1psiBasis->ValueRe[BlokElement(l)];
         k++;
       }
     }
     else
     {
       for(l=1;l<=H1psiBasis->N;l++)
       {
         MT1i.V[k].Re = H1psiBasis->ValueRe[BlokElement(l)];
         MT1i.V[k].Im = H1psiBasis->ValueIm[BlokElement(l)];
         k++;
       }
     }
   }
 }
 //2. Calcualte MfT2m[k] = a dense matrix rerpesenting T2 on a basis of H1psiBasis and H2psiBasis for each k from 0 to NT2
 //   The dimensions of MfT2m[k] are Ni = H2psiBasis.N by Nj = H1psiBasis.N
 for(k=0;k<NT2;k++)
 {
   OperatorToDeterminantMatrix2(&(T2[k]), H2psiBasis, H1psiBasis, &(MfT2m[k]));
 }
 //3. Calculate Hm the dense matrix representing H1 on a basis of H1psiBasis
 //   The dimensions of Hm are Ni = H1psiBasis.N by Nj = H1psiBasis.N
 OperatorToDeterminantMatrix(&OppH1, H1psiBasis, &Hm);
 //4. Calculate Hf the dense matrix representing H2 on a basis of H2psiBasis
 //   The dimensions of Hf are Ni = H2psiBasis.N by Nj = H2psiBasis.N
 OperatorToDeterminantMatrix(&OppH2, H2psiBasis, &Hf);
 //5. calculate Valm and Funm as the eigenvalues and functions of Hm
 //   The dimensions of Funm are Ni = H1psiBasis.N by Nj = H1psiBasis.N
 CompactMatrixEigensystem(&Hm, &Funm, Valm);
 //6. calculate A = Funm^* * MT1i^T
 //   The dimensions of A are Ni = H1psiBasis.N by Nj = Npis*NT1
 CompactMatrixMultiply('*', 'T', &Funm, &MT1i, &A);
 //7. calculate Valf and Funf as the eigenvalues and functions of Hf
 //   The dimensions of Funf are Ni = H2psiBasis.N by Nj = H2psiBasis.N
 CompactMatrixEigensystem(&Hf, &Funf, Valf);
 //8. calcualte B[k] = Funf^* * MfT2m[k] * Funm^T
 //   The dimensions of B[k] are Ni = H2psiBasis.N by Nj = H1psiBasis.N
 for(k=0;k<NT2;k++)
 {
   CompactMatrixMultiply('N', 'T', &Funm, &(MfT2m[k]), &C);
   CompactMatrixMultiply('*', 'T', &Funf, &C, &(B[k]));
   FreeCompactMatrix(&C);
 }
 // free eigenfunctions
 FreeCompactMatrix(&Funm);
 FreeCompactMatrix(&Funf);
 //
 /* determine Emin1, Emax1, NE1 and Gamma1 if not set */
 if(settings->Emin1>( 0.001*DBL_MAX) && (settings->Emax1<(-0.001*DBL_MAX)))
 {
   settings->Emin1=1.1*Valm[0]-0.1*Valm[H1psiBasis->N-1]-energy1ListMax;
   settings->Emax1=1.1*Valm[H1psiBasis->N-1]-0.1*Valm[0]-energy1ListMin;
   VerbosityPrintf(0x0002, "Automatically set energy range for intermediate states Emin=%22.15E Emax=%22.15E\n",settings->Emin1,settings->Emax1);
 }
 if(settings->Emin1>( 0.001*DBL_MAX))
 {
   settings->Emin1=1.1*Valm[0]-0.1*settings->Emax1-energy1ListMax;
   VerbosityPrintf(0x0002, "Automatically set energy range for intermediate states Emin=%22.15E Emax=%22.15E\n",settings->Emin1,settings->Emax1);
 }
 if(settings->Emax1<(-0.001*DBL_MAX))
 {
   settings->Emax1=1.1*Valm[H1psiBasis->N-1]-0.1*settings->Emin1-energy1ListMin;
   VerbosityPrintf(0x0002, "Automatically set energy range for intermediate states Emin=%22.15E Emax=%22.15E\n",settings->Emin1,settings->Emax1);
 }
 if((settings->Gamma1==0)&&(settings->NE1==0))
 {
   settings->NE1=1000;
   settings->Gamma1=5*(settings->Emax1-settings->Emin1)/settings->NE1;
   VerbosityPrintf(0x0002, "Automatically set number of energy points (NE=%d) and lifetime (Gamma=%22.15E) for intermediate state\n",settings->NE1,settings->Gamma1);
 }
 if(settings->NE1==0)
 {
   settings->NE1    =(uint32_t)(5*(settings->Emax1-settings->Emin1)/settings->Gamma1);
   VerbosityPrintf(0x0002, "Automatically set number of energy points (NE=%d) for intermediate state\n",settings->NE1);
 }
 if(settings->Gamma1==0)
 {
   settings->Gamma1 =(5*(settings->Emax1-settings->Emin1)/settings->NE1);
   VerbosityPrintf(0x0002, "Automatically set lifetime (Gamma=%22.15E) for intermediate state\n",settings->Gamma1);
 }
 /* determine Emin2, Emax2, NE2 and Gamma2 if not set */
 if(settings->Emin2>( 0.001*DBL_MAX) && (settings->Emax2<(-0.001*DBL_MAX)))
 {
   settings->Emin2=1.1*Valf[0]-0.1*Valf[H2psiBasis->N-1]-energy2ListMax;
   settings->Emax2=1.1*Valf[H2psiBasis->N-1]-0.1*Valf[0]-energy2ListMin;
   VerbosityPrintf(0x0002, "Automatically set energy range for final states Emin=%22.15E Emax=%22.15E\n",settings->Emin2,settings->Emax2);
 }
 if(settings->Emin2>( 0.001*DBL_MAX))
 {
   settings->Emin2=1.1*Valf[0]-0.1*settings->Emax2-energy2ListMax;
   VerbosityPrintf(0x0002, "Automatically set energy range for final states Emin=%22.15E Emax=%22.15E\n",settings->Emin2,settings->Emax2);
 }
 if(settings->Emax2<(-0.001*DBL_MAX))
 {
   settings->Emax2=1.1*Valf[H2psiBasis->N-1]-0.1*settings->Emin2-energy2ListMin;
   VerbosityPrintf(0x0002, "Automatically set energy range for final states Emin=%22.15E Emax=%22.15E\n",settings->Emin2,settings->Emax2);
 }
 if((settings->Gamma2==0)&&(settings->NE2==0))
 {
   settings->NE2=1000;
   settings->Gamma2=5*(settings->Emax2-settings->Emin2)/settings->NE2;
   VerbosityPrintf(0x0002, "Automatically set number of energy points (NE=%d) and lifetime (Gamma=%22.15E) for final state\n",settings->NE2,settings->Gamma2);
 }
 if(settings->NE2==0)
 {
   settings->NE2    =(uint32_t)(5*(settings->Emax2-settings->Emin2)/settings->Gamma2);
   VerbosityPrintf(0x0002, "Automatically set number of energy points (NE=%d) for final state\n",settings->NE2);
 }
 if(settings->Gamma2==0)
 {
   settings->Gamma2 =(5*(settings->Emax2-settings->Emin2)/settings->NE2);
   VerbosityPrintf(0x0002, "Automatically set lifetime (Gamma=%22.15E) for final state\n",settings->Gamma2);
 }
 //   For each inicdent photon energy wi[l] do (l from 0 to Nwi)
 C.Ni=A.Ni;
 C.Nj=A.Nj;
 C.QComplex=1;
 InitCompactMatrix(&C);
 for(l = 0; l <= settings->NE1; l++)
 {
   //9. Calcualte the matrix C as a copy of A with each element scaled: Cij = Aij * 1/(wi[l] - Em[i] + Ei[j] + i Gm)
   //   The dimensions of C are Ni = H1psiBasis.N by Nj = Npsi*NT1
   /* Calculate 1/(E_0 + wres - EVal[i] + iGammaRes), split into real and imaginary */
   wRes = settings->Emin1 + (settings->Emax1-settings->Emin1)*l/((double)settings->NE1);
   o=0;
   if(A.QComplex==0)
   {
     for(i=0;i<H1psiBasis->N;i++)
     {
       for(j=0;j<Npsi;j++)
       {
         Lorentz.Re = (wRes - Valm[i] + energy1List[j])/(Sqr(wRes - Valm[i] + energy1List[j]) + Sqr(settings->Gamma1*0.5));
         Lorentz.Im = -(settings->Gamma1 * 0.5)/(Sqr(wRes - Valm[i] + energy1List[j]) + Sqr(settings->Gamma1*0.5));
         for(k=0;k<NT1;k++)
         {
           C.V[o].Re = A.VRe[o] * Lorentz.Re;
           C.V[o].Im = A.VRe[o] * Lorentz.Im;
           o++;
         }
       }
     }
   }
   else
   {
     for(i=0;i<H1psiBasis->N;i++)
     {
       for(j=0;j<Npsi;j++)
       {
         Lorentz.Re = (wRes - Valm[i] + energy1List[j])/(Sqr(wRes - Valm[i] + energy1List[j]) + Sqr(settings->Gamma1*0.5));
         Lorentz.Im = -(settings->Gamma1 * 0.5)/(Sqr(wRes - Valm[i] + energy1List[j]) + Sqr(settings->Gamma1*0.5));
         for(k=0;k<NT1;k++)
         {
           C.V[o].Re = A.V[o].Re * Lorentz.Re - A.V[o].Im * Lorentz.Im;
           C.V[o].Im = A.V[o].Re * Lorentz.Im + A.V[o].Im * Lorentz.Re;
           o++;
         }
       }
     }
   }
   //   For each decay operator T2[k] do
   for(k=0;k<NT2;k++)
   {
     //10. Calculate the matrix R[k] = B[k] * C
     //    The dimensions of R[k] are Ni = H2psiBasis.N by Nj = Npsi*NT1
     //    R[k]_ij contains the scattering amplitude from initial T1psi[j] to final state i with resonant energy wi[l] and decay operator T2[k]
     CompactMatrixMultiply('N', 'N', &(B[k]), &C, &(R[k]));
   }

  //11. Compute the spectra for each combinations of Psi and T1,T2 
  //    each tensor element is spec[t]= <Psi[s]|T1_dag[m] G1_dag(w1) T2_dag[n] * G2(w2) * T2[j] * G1(w1) * T1[i]|Psi[r]>
  //    where t= l + (settings->NE1+1) * (n + NT2 *(m+NT1*(j+NT2*(i+NT1*(s+Npsi*r)))));
  
  for(r=0;r<Npsi;r++)
  {
    for(s=0;s<Npsi;s++)
    {
      for(i=0;i<NT1;i++)
      {
        for(j=0;j<NT2;j++)
        {
          for(m=0;m<NT1;m++)
          {
            for(n=0;n<NT2;n++)
            {

              // OLD: define indexing output, from outest to innermost: Psiright[r],psileft[s],T1[i],T2[j],T1_dag[m],T2_dag[n],omega1[l]
              // NEW: define indexing output, from outest to innermost: psileft[s],Psiright[r],T1_dag[m],T2_dag[n],T1[i],T2[j],omega1[l]

              t = l + (settings->NE1+1) * (j + NT2 *(i+NT1*(n+NT2*(m+NT1*(r+Npsi*s)))));
              
              //select right and ledt T|psi>
              kRr = i+r*NT1;
              kRs = m+s*NT1;

              // allocate the spectra
              (*Spec)[t].N    =settings->NE2;
              (*Spec)[t].Emin =settings->Emin2;
              (*Spec)[t].Emax =settings->Emax2;
              (*Spec)[t].dE   =(settings->Emax2-settings->Emin2)/settings->NE2;
              (*Spec)[t].Gamma=settings->Gamma2;
              (*Spec)[t].ResonantEnergy = settings->Emin1 + ((double)(l)/(settings->NE1)) * (settings->Emax1-settings->Emin1);
              sprintf((*Spec)[t].Name, "Spectrum wres %i/%i T2_dag %i/%i T1_dag %i/%i T2 %i/%i T1 %i/%i psi_left %i/%i psi_right %i/%i",l,settings->NE1,n,NT2,m,NT1,j,NT2,i,NT1,s,Npsi,r,Npsi);
              SpectrumAllocate((*Spec)[t], &(*Spec)[t]);
              //    calculate the integral 
              for(a=0;a<H2psiBasis->N;a++)
              {
                //TODO: does it make sense for the Tensor off diag elements to have integral defined?
               // (*Spec)[t].b[0]  += R[j].V[kRr+a*R[j].Nj].Re*R[n].V[kRs+a*R[n].Nj].Re;
                //(*Spec)[t].b[0]  += R[j].V[kRr+a*R[j].Nj].Im*R[n].V[kRs+a*R[n].Nj].Im;
              }
              //    For each transfer energy wr[m] do (m from 0 to Nwr)
              for(b=0;b<settings->NE2+1;b++)
              {
                wLoss = settings->Emin2 + ((double)b/(double)(settings->NE2))*(settings->Emax2 - settings->Emin2);
                //12. Spec[T1psi[j],T2[k],wi[l]][wr[m]] = sum_a R^{*}_a * R_a/ (wr[m] - Ef[a] + Ei[j] + i Gf)
                //    The index of the spectrum is:  l + Nwr * ( k + NT2 * ( jT1 + NT1 * ( jpsi ) ) )  where j = ( jT1 + NT1 * ( jpsi ) )
                for(a=0;a<H2psiBasis->N;a++)
                {
                  //imagionary and real parts of the lorentzian for the final state
                  double tempDoubleRe = (wLoss - Valf[a] + energy2List[0]) / (Sqr(wLoss - Valf[a] + energy2List[0]) + Sqr(settings->Gamma2*0.5));
                  double tempDoubleIm = -(settings->Gamma2*0.5)             / (Sqr(wLoss - Valf[a] + energy2List[0]) + Sqr(settings->Gamma2*0.5));
                  //TODO: energy2List[0] how to deal with the offset for cross components between wavefunctions?
                  
                  //The dimensions of R[k] are Ni = H2psiBasis.N by Nj = Npsi*NT1
                  //sum_a R^*_{a,kRs}*R^*_{a,kRr}*Gf_a
                  (*Spec)[t].Im[b] += R[n].V[kRs+a*R[n].Nj].Im* R[j].V[kRr+a*R[j].Nj].Im * tempDoubleIm +
                                      R[n].V[kRs+a*R[n].Nj].Re* R[j].V[kRr+a*R[j].Nj].Re * tempDoubleIm +
                                      R[n].V[kRs+a*R[n].Nj].Re* R[j].V[kRr+a*R[j].Nj].Im * tempDoubleRe -
                                      R[n].V[kRs+a*R[n].Nj].Im* R[j].V[kRr+a*R[j].Nj].Re * tempDoubleRe;

                  (*Spec)[t].Re[b] += R[n].V[kRs+a*R[n].Nj].Im* R[j].V[kRr+a*R[j].Nj].Im * tempDoubleRe +
                                      R[n].V[kRs+a*R[n].Nj].Re* R[j].V[kRr+a*R[j].Nj].Re * tempDoubleRe -
                                      R[n].V[kRs+a*R[n].Nj].Re* R[j].V[kRr+a*R[j].Nj].Im * tempDoubleIm +
                                      R[n].V[kRs+a*R[n].Nj].Im* R[j].V[kRr+a*R[j].Nj].Re * tempDoubleIm;
                } // end loop over all final states
              } // end loop over energy loss photon energies index b
              //(*Spec)[t].b[0] = sqrt((*Spec)[t].b[0]);
              if(H2psiBasis->N==0)
              {
                (*Spec)[t].EmaxPole=0;
                (*Spec)[t].EminPole=0;
              }
              else
              {
                //TODO: energy2List[0] how to deal with the offset for cross components between wavefunctions?
                (*Spec)[t].EmaxPole = Valf[H2psiBasis->N-1] -  energy2List[0];
                (*Spec)[t].EminPole = Valf[0] -  energy2List[0];
              }
            } //end loop 
          } // end loop 
        } // end loop over T1psi indexed by j
      } // end loop over decay operators index k
    }//end loop over psi
  }//end loop over psi 
 } // end loop over incomming photon energies index l
 FreeCompactMatrix(&A);
 for(k=0;k<NT2;k++)
 {
   FreeCompactMatrix(&(B[k]));
   FreeCompactMatrix(&(R[k]));
 }
 FreeCompactMatrix(&C);
 free(Valm);
 free(Valf);
 free(MfT2m);
 return 0;
}



uint32_t ResonantSpectraPartialDense(const OperatorType OppH1, const OperatorType OppH2, const OperatorType *T1, uint32_t NT1, const OperatorType *T2, uint32_t NT2, const WaveFunctionType *psi, uint32_t Npsi, TriDiagonalMatrixType **SpecTri2, SpectraType **Spec, double energy1, double energy2, ConvergenceSettingsType *settings)
{
  /*Calculates resonant spectra using <i|T*GTGT*GT|i> = sum_f |sum_m <f|T2|m><m|T1|i> 1/(W-Em+energy1+iG1/2)|^2 1/(w-Ef+energy2 + iG2/2) in the following steps:
   
   1.  Calculate T1|i> for all T1 and |i>
   2.  Use all T1|i>'s to transform H1 to tridiagonal form in Krylov basis using Block Lanczos
   3.  Calculate eigen-states of Krylov basis (|m>) and eigen values Em
   4.  Calculate matrix <m|T1|i>
   5.  Calculate T2 |m> for all |m>'s
   6.  Use all T2 |m>'s to transform H2 to tridiagonal form in Krylov basis using Block Lanczos
   7.  Calculate eigen-states of Krylov basis |f> and eigen values Ef
   8.  Calculate matrix <f|T2|m> for all T2's
   9.  for each W calculate M_{fT2, T1i}(Em) = |sum_m <f|T2|m><m|T1|i> 1/(W-Em+energy1List[i]+iG1/2)|^2
   10. for each w calcualte sum_f M_{fT2, T1i}(Em) 1/(w-Ef+energy2List[i] + iG2/2)
   */
  
  /* settings that need to be set:
   * NTri1, NTri2
   * Emin1, Emax1, NE1, Gamma1 - for automatic determination use DBL_MAX, -DBL_MAX, 0, 0
   * Emin2, Emax2, NE2, Gamma2 - for automatic determination use DBL_MAX, -DBL_MAX, 0, 0
   * epsilon
   */
  
  WaveFunctionType *Tpsi, *PsiBasis[settings->NTri1];
  uint32_t ipsi,iT1,iT1psi,i,j,k;
  double energy1List[Npsi],energy2List[Npsi], *H1T1psival;
  BlockTriDiagonalMatrixType H1T1psi;
  CompactMatrixType H1T1psifun, mKTpsi;
  uint32_t BlockSize[settings->NTri1+1];
  BlockSize[0]=NT1*Npsi;
  
  // calculate energy?List
  for(ipsi=0;ipsi<Npsi;ipsi++)
  {
    if(energy1==DBL_MAX)
    {
      PsiHermitianOperatorPsi(OppH1, psi[ipsi], &energy1List[ipsi]);
    }
    else
    {
      energy1List[ipsi]=energy1;
    }
    if(energy2==DBL_MAX)
    {
      PsiHermitianOperatorPsi(OppH2, psi[ipsi], &energy2List[ipsi]);
    }
    else
    {
      energy2List[ipsi]=energy2;
    }
  }
  // 1. calculate T1|i>
  Tpsi= (WaveFunctionType*)malloc(NT1*Npsi*sizeof(WaveFunctionType));
  for(ipsi=0;ipsi<Npsi;ipsi++)
  {
    for(iT1=0;iT1<NT1;iT1++)
    {
      if(settings->restrictions1.N==0)
      {
        OperatorPsi(T1[iT1], psi[ipsi], &Tpsi[ipsi*NT1+iT1]);
      }
      else
      {
        OperatorPsiRestricted(T1[iT1], psi[ipsi], &Tpsi[ipsi*NT1+iT1], settings->restrictions1);
      }
    }
  }
  // 2. Block diagonalize H1 starting from Tpsi
  /* Use T|i>==Tpsi to transform H to Block tridiagonal form in Krylov basis */
  /* Using the version which gives back the Krylov basis vectors */
  H1T1psi.N = settings->NTri1;
  H1T1psi.BlockSize = Npsi*NT1;
  H1T1psi.QComplex = 0;
  if(OppH1.QComplex!=0) H1T1psi.QComplex = 1;
  for(iT1psi=0;iT1psi<H1T1psi.BlockSize;iT1psi++)
  {
    if(Tpsi[iT1psi].QComplex!=0) H1T1psi.QComplex = 1;
  }
  InitBlockTridiagonalMatrix(&H1T1psi);
  if(settings->restrictions1.N==0)
  {
    LanczosBlockTriDiagonalizeKrylov(OppH1, Tpsi, settings->epsilon, settings->SingularValue, &H1T1psi, PsiBasis, &BlockSize[1]);
  }
  else
  {
    LanczosBlockTriDiagonalizeKrylovRestricted(OppH1, Tpsi, settings->restrictions1, settings->epsilon, settings->SingularValue, &H1T1psi, PsiBasis, &BlockSize[1]);
  }
  // 3.  Calculate eigen-states of Krylov basis (|m>) and eigen values Em
  H1T1psival= (double*)malloc(H1T1psi.N*H1T1psi.BlockSize*sizeof(double));
  BlockTriDiagonalMatrixEigenSystemWithDeflation(&H1T1psi, H1T1psival, &H1T1psifun, BlockSize);
  // 4.  Calculate matrix <m|T1|i>
  mKTpsi.Ni=Npsi*NT1;
  mKTpsi.Nj=H1T1psifun.Ni;
  mKTpsi.QComplex=H1T1psi.QComplex;
  InitCompactMatrix(&mKTpsi);
  k=0;
  for(i=0;i<mKTpsi.Ni;i++)
  {
    for(j=0;i<mKTpsi.Nj;j++)
    {
      mKTpsi.VRe[k]=0;
      //RealWaveFunctionDotProduct(PsiBasis[j], <#const WaveFunctionType *b#>);
      k++;
    }
  }
  
  
  
  // 5.  Calculate T2 |m> for all |m>'s
  
  // 6.  Use all T2 |m>'s to transform H2 to tridiagonal form in Krylov basis using Block Lanczos
  
  // 7.  Calculate eigen-states of Krylov basis |f> and eigen values Ef
  
  // 8.  Calculate matrix <f|T2|m> for all T2's
  
  // 9.  for each W calculate M_{fT2, T1i}(Em) = |sum_m <f|T2|m><m|T1|i> 1/(W-Em+energy1List[i]+iG1/2)|^2
  
  // 10. for each w calcualte sum_f M_{fT2, T1i}(Em) 1/(w-Ef+energy2List[i] + iG2/2)
  
  free(Tpsi);
  return 0;
}


/*  Calculate < psi | T1^+ G1^+ T2^+ G2 T2 G1 T1 | psi >
 *  psi, T1 and T2 are vectors of Length Np, NT1 and NT2
 *  psi and T1 are contract to a single vector of length Np*NT1
 *  returned is a matrix of dimensions Sqr(Npsi* NT1 * NT2) in of flattened array of spetra, where the innest index runs over the excitation energy wres.*/
uint32_t ResonantSpectraTensor(const OperatorType OppH1, const OperatorType OppH2, const OperatorType *T1, uint32_t NT1, const OperatorType *T2, uint32_t NT2, const WaveFunctionType *psi, uint32_t Npsi, TriDiagonalMatrixType **SpecTri2, SpectraType **Spec, double energy1, double energy2, ConvergenceSettingsType *settings)
{
  /*Given wavefunction psi (usually ground state, but could be excited)*/
  
  /*Calculates resonant spectra using <i|T*GTGT*GT|i> in the following steps:
   
   1. Calculate T1|i> for all elments of T1 and |i> giving a vector of length Npsi*NT1 states and put all on the same basis.
   2. Use T1|i> to transform H1 to block tridiagonal form in Krylov basis
   3. Compute the Block tridiagonal eigenvalues and eigenfunctions of H1 in the Krylov basis
   for each wRes do:
     4. Invert (wRes - H_tri + iGamma) to get G(wRes) as a matrix of dimension (Npsi*NT1*Ntri1)^2.
     5. Calculate G(wRes)T|i> as a vector of length Npsi*NT1
     6. Calculate T2*G(wRes)T1|i> as a vector of length Npsi*NT1* NT2
     7. Use T2*G(wRes)T1|i> to transform H2 to Block tridiagonal and get GTri2 for the same wRes
     8. Convert to spectra. Here, the block algorithm automatically combines into Tensor form, just as in CreateSpectra.
     9. Save in *Spec with the desired indexing.
   */
  
  /* settings that need to be set:
   * NTri1, NTri2
   * Emin1, Emax1, NE1, Gamma1 - for automatic determination use DBL_MAX, -DBL_MAX, 0, 0
   * Emin2, Emax2, NE2, Gamma2 - for automatic determination use DBL_MAX, -DBL_MAX, 0, 0
   * epsilon
   */

  /*  The greensfunctions are "shifted" by energy1 and energy2.
   *  If energy1=DBL_MAX or energy2=DBL_MAX on input the shift is calculated as energy1 = <psi|H1|psi> and energy2 = <psi|H2|psi> this is well defined if Npsi=1
   *  If Npsi>1 and energy? is not DBL_MAX energy? of the input is used for all states
   *  If energy?=DBL_MAX we set energy? = <psi[0] | H? | psi[0]> and print a warning
   */

  

  //Possibly TODO:
  //      - Krylov true Option.
  //      - Also save response function  Gtri2

  WaveFunctionType Tpsi[Npsi*NT1],GTpsi[Npsi*NT1],TGTpsi[Npsi*NT1*NT2],**PsiBasis=NULL;
  BlockTriDiagonalMatrixType SpecTri1;
  CompactMatrixType EFun,G1Kryl;
  uint32_t *BlockSizeN=NULL;
  Complex16 *Gdiag=NULL;
  double EVal[Npsi*NT1*settings->NTri1], EValInvRe[Npsi*NT1*settings->NTri1], EValInvIm[Npsi*NT1*settings->NTri1], tempDoubleRe,tempDoubleIm,wRes,*Eval2,Emin2pole,Emax2pole;
  int i,ipsi,iT1,p,m,n,j,o,r,s;
  

  /* calcualte ground state energy */
  // TODO: add this warning to the dense version as well

  if(energy1==DBL_MAX)
  {
    if(Npsi!=1)
    {
      printf("Warning in function ResonantSpectraTensor\nFunction called with %i wavefunctions but no fixed energy shift set.\n",Npsi);
      printf("For the calculation of G1 we assume G1 = 1/(w-H1+energy1 + i Gamma/2) with energy1 = <psi[0] | H1 | psi[0]> for all psi[i].\n");
      printf("This does not seem to make physical sence to me.\n");
      printf("If you are convinced that what you do is ok, please contact me to discuss and I will change this warning.\n");
    }
    PsiHermitianOperatorPsi(OppH1, psi[0], &energy1);
  }
  if(energy2==DBL_MAX)
  {
    if(Npsi!=1)
    {
      printf("Warning in function ResonantSpectraTensor\nFunction called with %i wavefunctions but no fixed energy shift set.\n",Npsi);
      printf("For the calculation of G2 we assume G2 = 1/(w-H2+energy2 + i Gamma/2) with energy2 = <psi[0] | H2 | psi[0]> for all psi[i].\n");
      printf("This does not seem to make physical sence to me.\n");
      printf("If you are convinced that what you do is ok, please contact me to discuss and I will change this warning.\n");
    }
    PsiHermitianOperatorPsi(OppH2, psi[0], &energy2);
  }
  
  
  //1. Calculate T1|i> for all elments of T1 and |i> giving a vector of length Npsi*NT1 
  i=0;

  for(ipsi=0;ipsi<Npsi;ipsi++)
  {
    for(iT1=0;iT1<NT1;iT1++)
    {
      {
        if(settings->restrictions1.N==0)
        {
          OperatorPsi(T1[iT1],psi[ipsi],&Tpsi[i]);
        }
        else
        {
          OperatorPsiRestricted(T1[iT1],psi[ipsi],&Tpsi[i], settings->restrictions1);
        }
      }
      if(Tpsi[i].QComplex)
      {
        ComplexWaveFunctionTruncateBasis(&Tpsi[i],settings->epsilon);
      }
      else
      {
        RealWaveFunctionTruncateBasis(&Tpsi[i],settings->epsilon);
      }
      i++;
    }
  }




  /* Copy the created states to put them on the same basis */
  WaveFunctionType TpsiSameBasis[NT1 * Npsi];
  
  WaveFunctionCopy(Tpsi[0], &TpsiSameBasis[0]);
  WaveFunctionFree(&(Tpsi[0]));
  for(i=1;i< NT1 * Npsi;i++)
  {
    if(WaveFunctionAddMultipliedUnionBasis(&TpsiSameBasis[0], Tpsi[i], 0)) printf("malloc failed. Out of memory during initialization in CalculateSpectra\n");
  }
  for(i=1;i< NT1 * Npsi;i++)
  {
    if(WaveFunctionInitCopyBasis(&(TpsiSameBasis[i]), TpsiSameBasis[0])) printf("malloc failed. Out of memory during initialization in CalculateSpectra\n");
    WaveFunctionAddConserveBasis(&(TpsiSameBasis[i]), Tpsi[i]);
    WaveFunctionFree(&(Tpsi[i]));
  }

 
  
  //2. Use T|i>==Tpsi to transform H1 to Block tridiagonal form in Krylov basis */
  //   Using the version which gives back the Krylov basis vectors */
  SpecTri1.N = settings->NTri1;
  SpecTri1.BlockSize = Npsi*NT1;
  SpecTri1.QComplex = 0;

  if(OppH1.QComplex!=0) SpecTri1.QComplex = 1;
  for(i=0;i<SpecTri1.BlockSize;i++)
  {
    if(TpsiSameBasis[i].QComplex!=0) SpecTri1.QComplex = 1;
  }
  InitBlockTridiagonalMatrix(&SpecTri1);
  SpecTri1.mu=0;//psiEn[0];


  BlockSizeN= (uint32_t*)malloc(settings->NTri1*sizeof(uint32_t));
  PsiBasis= (WaveFunctionType**)malloc(settings->NTri1*sizeof(WaveFunctionType*));
  
  if(settings->restrictions1.N==0)
  {
    LanczosBlockTriDiagonalizeKrylov(OppH1, TpsiSameBasis, settings->epsilon, settings->SingularValue, &SpecTri1, PsiBasis, BlockSizeN);
  }
  else
  {
    LanczosBlockTriDiagonalizeKrylovRestricted(OppH1, TpsiSameBasis, settings->restrictions1, settings->epsilon, settings->SingularValue, &SpecTri1, PsiBasis, BlockSizeN);
  }
  BlockTridiagonalMatrixChop(&SpecTri1, settings->epsilon);

  printf("SpecTri1.b[0] first 3 x 3 block:\n");
  if(SpecTri1.N==0 || SpecTri1.BlockSize==0)
  {
    printf("  (empty)\n");
  }
  else
  {
    for(i=0;i<Min(3,SpecTri1.BlockSize);i++)
    {
      for(j=0;j<Min(3,SpecTri1.BlockSize);j++)
      {
        if(SpecTri1.QComplex)
        {
          printf("(% .15E,% .15E) ", SpecTri1.b[i*SpecTri1.BlockSize+j].Re, SpecTri1.b[i*SpecTri1.BlockSize+j].Im);
        }
        else
        {
          printf("% .15E ", SpecTri1.bRe[i*SpecTri1.BlockSize+j]);
        }
      }
      printf("\n");
    }
  }


  /*3. Compute the Block tridiagonal eigenvalues and eigenfunctions of H1 in the Krylov basis.
  The standard eygensistem function computes SpecTri1 = Efun^T * D(Eval) * Efun, where D is the diagonal matrix of eigenvalues.
   */
  /* Eval has length (settings->NTri1)*(Npsi*NT1), with inner loop over (Npsi*NT1)*/
  //if(BlockTriDiagonalMatrixEigenSystemWithDeflation(&SpecTri1, EVal, &EFun, &SpecTri1.BlockSize))
  if(BlockTriDiagonalMatrixEigenSystem(&SpecTri1, EVal, &EFun))
  {
    printf("Problem with BlockTridiagonalMatrixToEigenSystem\n");
  }
  //we take the conjugate transpose of EFun for later uses
  CompactMatrixConjugateTranspose(&EFun);


  /* determine Eminpole and Emaxpole */
  /* get eigenvalues Block Tridiagonal matrix and set smallest and largest eigenvalue */
  {
  double Emin1pole = DBL_MAX;
  double Emax1pole =-DBL_MAX;
  Emin1pole=EVal[0]-energy1;
  Emax1pole=EVal[SpecTri1.N*SpecTri1.BlockSize-1]-energy1;

  /* if Emin=DBL_MAX or Emax=-DBL_MAX or Gamma==0 determine these */
  if(Emin1pole==Emax1pole)
  {
    if(Emin1pole==0)
    {
      if(settings->Emin1== DBL_MAX) settings->Emin1 =-1;
      if(settings->Emax1==-DBL_MAX) settings->Emax1 = 1;
    }
    if(Emin1pole<0)
    {
      if(settings->Emin1== DBL_MAX) settings->Emin1 = 2.125*Emin1pole;
      if(settings->Emax1==-DBL_MAX) settings->Emax1 =-0.125*Emin1pole;
    }
    if(Emin1pole>0)
    {
      if(settings->Emin1== DBL_MAX) settings->Emin1 =-0.125*Emin1pole;
      if(settings->Emax1==-DBL_MAX) settings->Emax1 = 2.125*Emin1pole;
    }
  }
  if(settings->Emin1== DBL_MAX) settings->Emin1 = Emin1pole - 0.125 * (Emax1pole-Emin1pole);
  if(settings->Emax1==-DBL_MAX) settings->Emax1 = Emax1pole + 0.125 * (Emax1pole-Emin1pole);
  if(settings->Gamma1==0) settings->Gamma1 = 10*(settings->Emax1-settings->Emin1)/((double)settings->NE1);
  }

  
  /* allocate memory for the output */
  *Spec = (SpectraType*)malloc((Sqr(Npsi * NT1 * NT2)*(settings->NE1+1)) * sizeof(SpectraType));

  
  /*Start loop over NE1+1 excitation energies*/
  for(p = 0; p <= settings->NE1; p++)
  {
    /*
    4. Invert (w - H_tri + iGamma) to get G(w) as a matrix of dimension (Npsi*NT1*Ntri1)^2:
     Calculate 1/(E_0 + wres - EVal[i] + iGammaRes), split into real and imaginary 
    */
  
    Gdiag = (Complex16*)malloc((settings->NTri1)*(Npsi)*(NT1) * sizeof(Complex16));
    wRes = settings->Emin1 + (settings->Emax1-settings->Emin1)*p/((double)settings->NE1) + energy1;

    for(i = 0; i <(settings->NTri1)*(Npsi)*(NT1); i++)
    {
      Gdiag[i].Re =  (wRes - EVal[i])       / ((wRes - EVal[i])*(wRes - EVal[i])+(settings->Gamma1/2.0)*(settings->Gamma1/2.0));
      Gdiag[i].Im = -(settings->Gamma1/2.0) / ((wRes - EVal[i])*(wRes - EVal[i])+(settings->Gamma1/2.0)*(settings->Gamma1/2.0));
    }

    /*5. Calculate G(w)T|i> as a vector of length Npsi*NT1
    
    The goal is to compute a list of vectors GTpsi[j] = G(wres)T|i> for all {T1|i>} combination indexed by j=0,...,Npsi*NT1-1.
    We want tocompute:
      GTpsi[j] = sum_{m}^{NTri1}sum_{n}^{Blocksize[m]} alpha^{j}_{m,n} |PsiBasis[m][n]>
    where
      alpha^{j}_{m,n} = sum_{i}^{Npsi*NT1} sum_{r}^{Npsi*NT1*NTri1} Efun_{r,n+m*Npsi*NT1} * Gdiag[r] * Efun^{*}_{r,i} * B0_{j,i}.
    */
    
    /*First we rotate the diagonal propagator into the basis of Krylov states */
    /* we use GKryl = Efun^T * Gdiag * Efun^* */
    /* Note that before the p loop we take the ComplexTranspose of Efun*/
    TransformComplexDiagonalMatrixByCompactMatrixTconATtran(Gdiag, &EFun, &G1Kryl);

    o = 0; //index for TGT|i>
    for (j=0;j< NT1 * Npsi;j++)
    {
      GTpsi[j].NFermion = TpsiSameBasis[j].NFermion;
      GTpsi[j].NBoson   = TpsiSameBasis[j].NBoson;
      GTpsi[j].QComplex = TpsiSameBasis[j].QComplex;
      GTpsi[j].NBitsKey = TpsiSameBasis[j].NBitsKey;
      if(WaveFunctionInit(&GTpsi[j]))
      {
        printf("WaveFunctionInit failed in ResonantSpectra 00");
        return 1;
      }
 
      for(m = 0; m < settings->NTri1 ; m++)
      {
        fflush(stdout);
        if(BlockSizeN[m]==0) break;
        for(n = 0; n < BlockSizeN[m]; n++)
        {
          uint32_t idx;
          double alphaRe = 0;
          double alphaIm = 0;
          for (i=0;i< NT1 * Npsi;i++)
          { 

            /*
            alpha^{j}_{m,n} = sum_{i}^{Npsi*NT1} GKryl[n+m*Npsi*NT1,i] * B0_[j,i]
            */

            idx = i+(NT1*Npsi*(0+settings->NTri1*(n+(NT1*Npsi*m)))); 

            if(SpecTri1.QComplex)
            {
              alphaRe += SpecTri1.b[j*(NT1*Npsi)+i].Re *  G1Kryl.V[idx].Re - SpecTri1.b[j*(NT1*Npsi)+i].Im *  G1Kryl.V[idx].Im;
              alphaIm += SpecTri1.b[j*(NT1*Npsi)+i].Re *  G1Kryl.V[idx].Im + SpecTri1.b[j*(NT1*Npsi)+i].Im *  G1Kryl.V[idx].Re;
            }
            else
            {
              alphaRe += SpecTri1.bRe[j*(NT1*Npsi)+i] *  G1Kryl.V[idx].Re;
              alphaIm += SpecTri1.bRe[j*(NT1*Npsi)+i] *  G1Kryl.V[idx].Im;
            }             
          }//close loop over i

          WaveFunctionAddComplexMultipliedUnionBasis(&GTpsi[j], PsiBasis[m][n], alphaRe, alphaIm);//check the minus alpha here
        }//close loop over n
      }//close loop over m
   
      // 6. Calculate T2*G(w)T1|i> as a vector of length Npsi*NT1* NT2
      for(s = 0; s < NT2; s++)
      {
        TGTpsi[o].NFermion = GTpsi[j].NFermion;
        TGTpsi[o].NBoson = GTpsi[j].NBoson;
        TGTpsi[o].QComplex = GTpsi[j].QComplex;
        TGTpsi[o].NBitsKey = GTpsi[j].NBitsKey;

        if(WaveFunctionInit(&TGTpsi[o]))
        {
          printf("WaveFunctionInit failed in ResonantSpectra 00");
          return 1;
        }

        /* Calculate T*G(w)T|i> */
        if(settings->restrictions2.N==0)
        {
          OperatorPsi(T2[s],GTpsi[j],&TGTpsi[o]);
        }
        else
        {
          OperatorPsiRestricted(T2[s], GTpsi[j], &TGTpsi[o], settings->restrictions2);
        }
        
        if(TGTpsi[o].QComplex)
        {
          ComplexWaveFunctionTruncateBasis(&TGTpsi[o],settings->epsilon);
        }
        else
        {
          RealWaveFunctionTruncateBasis(&TGTpsi[o],settings->epsilon);
        }
        o++;
      }//close loop over s (T2 operators)
    WaveFunctionFree(&(GTpsi[j]));
    } //close loop over j T|i>
    

    /* calculate spectra*/
    WaveFunctionType TGTpsiSameBasis[NT2*NT1*Npsi];
    BlockTriDiagonalMatrixType GTri2;

    /* place starting TGTpsi functions on same basis */
    WaveFunctionCopy(TGTpsi[0], &TGTpsiSameBasis[0]);
    WaveFunctionFree(&(TGTpsi[0]));
    for(i=1;i<NT2*NT1*Npsi;i++)
    {
      if(WaveFunctionAddMultipliedUnionBasis(&TGTpsiSameBasis[0], TGTpsi[i], 0))
      {
        printf("malloc failed. Out of memory during initialization in CalculateSpectra\n");
      }
    }
    for(i=1;i<NT2*NT1*Npsi;i++)
    {
      if(WaveFunctionInitCopyBasis(&(TGTpsiSameBasis[i]), TGTpsiSameBasis[0]))
      {
        printf("malloc failed. Out of memory during initialization in CalculateSpectra\n");
      }
      WaveFunctionAddConserveBasis(&(TGTpsiSameBasis[i]), TGTpsi[i]);
      WaveFunctionFree(&(TGTpsi[i]));
    }

    //7. Use T2*G(w)T1|i> to transform H2 to Block tridiagonal 
    /* init BlockTriDiagonal Matrix */
    GTri2.BlockSize = NT2*NT1*Npsi;
    GTri2.N = settings->NTri2;
    /* create block tridiagonal representation */
    GTri2.QComplex=0;
    if(OppH2.QComplex) GTri2.QComplex=1;
    for(i=0;i<GTri2.BlockSize;i++)
    {
      if(TGTpsiSameBasis[i].QComplex) GTri2.QComplex=1;
    }
    InitBlockTridiagonalMatrix(&GTri2);
    GTri2.mu=0;//psiEn[0];

  
    if(settings->restrictions2.N==0)
    {
      LanczosBlockTriDiagonalize(OppH2, TGTpsiSameBasis, settings->epsilon, settings->SingularValue, &GTri2, Max(settings->Gamma2,settings->epsilon));
    }
    else
    {
      LanczosBlockTriDiagonalizeRestricted(OppH2, TGTpsiSameBasis, settings->restrictions2, settings->epsilon, settings->SingularValue, &GTri2, Max(settings->Gamma2,settings->epsilon));
    }
    BlockTridiagonalMatrixChop(&GTri2, settings->epsilon);


    /* determine Eminpole and Emaxpole */
    /* get eigenvalues Block Tridiagonal matrix and set smallest and largest eigenvalue */
    
    //set Eminpole and Emaxpole (they are in the output and I need them for auto range determination
    Eval2= (double*)malloc(GTri2.N*GTri2.BlockSize*sizeof(double));
    BlockTriDiagonalMatrixEigenValues(&GTri2,Eval2);

    Emin2pole=Eval2[0]-energy2;
    Emax2pole=Eval2[GTri2.N*GTri2.BlockSize-1]-energy2;
    free(Eval2);
    /* if Emin=DBL_MAX or Emax=-DBL_MAX or Gamma==0 determine these */
    if(Emin2pole==Emax2pole)
    {
      if(Emin2pole==0)
      {
        if(settings->Emin2== DBL_MAX) settings->Emin2 =-1;
        if(settings->Emax2==-DBL_MAX) settings->Emax2 = 1;
      }
      if(Emin2pole<0)
      {
        if(settings->Emin2== DBL_MAX) settings->Emin2 = 2.125*Emin2pole;
        if(settings->Emax2==-DBL_MAX) settings->Emax2 =-0.125*Emin2pole ;
      }
      if(Emin2pole>0)
      {
        if(settings->Emin2== DBL_MAX) settings->Emin2 =-0.125*Emin2pole;
        if(settings->Emax2==-DBL_MAX) settings->Emax2 = 2.125*Emin2pole;
      }
    }
    if(settings->Emin2== DBL_MAX) settings->Emin2 = Emin2pole - 0.125 * (Emax2pole-Emin2pole);
    if(settings->Emax2==-DBL_MAX) settings->Emax2 = Emax2pole + 0.125 * (Emax2pole-Emin2pole);
    if(settings->Gamma2==0) settings->Gamma2 = 10*(settings->Emax2-settings->Emin2)/((double)settings->NE2);
    

    /*8. Convert GTri2 to spectra. 
    The block algorithm should automatically combine into tensor, just as in CreateSpectra */
    SpectraTypeArray Specwp; 
    Specwp.Name=(char*)malloc(sizeof(char));
    Specwp.Name[0]='\0';
    Specwp.NSpectra=Sqr(NT2 * NT1 * Npsi);
    Specwp.Spec= (SpectraType*)malloc(Specwp.NSpectra*sizeof(SpectraType));
    printf("Energy2: %f\n", energy2);
    BlockTriDiagonalMatrixToSpectraWithBroadening(energy2, &GTri2, Specwp.Spec, settings->Emin2, settings->Emax2, settings->NE2, settings->Gamma2);
    //ShiftBlockTridiagonalMatrix(&GTri2,energy2);

    for(r=0;r<Npsi;r++)
    {
      for(s=0;s<Npsi;s++)
      {
        for(i=0;i<NT1;i++)
        {
          for(j=0;j<NT2;j++)
          {
            for(m=0;m<NT1;m++)
            {
              for(n=0;n<NT2;n++)
              {

              // the expercted index from the Block Algorithm
              uint32_t idxOld = j + NT2*(i+NT1*(s+Npsi*(n+NT2*(m+NT1*r))));
              // put in the same indexing as the Dense method
              uint32_t idxNew = j + NT2*(i+NT1*(n+NT2*(m+NT1*(r+Npsi*s))));

              Specwp.Spec[idxOld].EminPole=Emin2pole;
              Specwp.Spec[idxOld].EmaxPole=Emax2pole;
            
              //ChopSpectra(&Specwp.Spec[i], settings->epsilon);
              Specwp.Spec[idxOld].ResonantEnergy = settings->Emin1 + ((double)(p % (settings->NE1+1))/(settings->NE1)) * (settings->Emax1-settings->Emin1);
              //sprintf(Specwp.Spec[idxOld].Name, "Spectrum wres %i/%i Psi %i/%i T1 %i/%i T2 %i/%i",p,settings->NE1,s,Npsi,r,Npsi,i,NT1,j,NT2);  
              SpectrumAllocate(Specwp.Spec[idxOld], &(*Spec)[p + idxNew*(settings->NE1+1)]);
              
              //copying to avoid changing after freeing Specwp.Spec[i]
              CopySpectra(Specwp.Spec[idxOld],  &(*Spec)[p + idxNew*(settings->NE1+1)]);
              SpectrumFree(&Specwp.Spec[idxOld]);
              }
            }
          }
        }
      }
    }    
    free(Specwp.Name);
    Specwp.NSpectra=0;
    FreeBlockTridiagonalMatrix(&GTri2);
  }
  free(BlockSizeN);
  WaveFunctionFree(*PsiBasis);
  FreeBlockTridiagonalMatrix(&SpecTri1);
  return 0;
}


/*  Calculate < psi | Texc^+ G TDec[i] G2 TDec[i]^+ G Texc | psi >
 *  Input expects TDec^+
 */
uint32_t ResonantSpectraOrig(const OperatorType OppH, const OperatorType OppH2, const OperatorType Texc, const OperatorType *TDec, uint32_t NDec,  const WaveFunctionType psi, double energy, double epsilon, uint32_t numwRes, SpectraType **spectra, double *wRes, double gammaRes, double gamma, double Gfwhm, uint32_t NKrylov1, uint32_t NKrylov2)
{
  /*Given wavefunction psi (usually ground state, but could be excited)*/
  
  /*
  Calculates resonant spectra using <i|T*GTGT*GT|i> in the following steps:
   
   1. Calculate T|i>
   2. Use T|i> to transform H to tridiagonal form in Krylov basis
   3. Invert (w - H_tri + iGamma) to get G(w)
   4. Calculate G(w)T|i>
   5. Calculate T*G(w)T|i>
   6. Use the vector from step 5 to calculate RIXS spectrum using Lanczos
   */
  
  OperatorType OppH2Rot;
  WaveFunctionType Tpsi,GTpsi,TGTpsi, TGTpsiRot, PsiBasis[NKrylov1];
  TriDiagonalMatrixType SpecTri1, SpecTri2;
  MatrixType tempMat,EFun,Density,DensityPart,TRot;
  double DensityValues[psi.NFermion];
  double EVal[NKrylov1], EValInvRe[NKrylov1], EValInvIm[NKrylov1],tempDoubleRe,tempDoubleIm;
  int i,p,j,k;
  int const rotate = 0;
  
  /* Calculate T|i> */
  OperatorPsi(Texc,psi,&Tpsi);
  
  /* Use T|i> to transform H to tridiagonal form in Krylov basis */
  /* Using the version which gives back the Krylov basis vectors */
  SpecTri1.N = NKrylov1;
  InitTriDiagonalMatrix(&SpecTri1);
  LanczosTriDiagonalizeKrylov(OppH,Tpsi,epsilon,&SpecTri1,PsiBasis);
  
  
  /*Compute eigenvalues and eigenvectors, then S^-1 = QAQ^T where
   Q are the eigenvectors and A is a matrix with the eigenvalue inverses on the diagonal.
   */
  if(TriDiagonalMatrixToEigensystem(SpecTri1, &EFun, EVal))
  {
    printf("Problem with TridiagonalMatrixToEigenSystem\n");
  }
  
  for(p = 0; p < numwRes; p++)
  {
    /* New method for GT|i>:
     1. Tridiagonalize H using T|i> (outside wres loop)
     2. Calculate eigenvalues (EVal) and eigenfunctions (EFun) of Htri (outside)
     3. For each wres
     a. EValInv[i] = 1/(wres - EVal[i] + iGammaRes)
     b. GT|i> = sum_j EFun[i][j]*EValIn[j]*K[j], where K is Krylov basis vectors
     */
    
    /* Calculate 1/(E_0 + wres - EVal[i] + iGammaRes), split into real and imaginary */
    for(i = 0; i < NKrylov1; i++)
    {
      EValInvRe[i] = (wRes[p] + energy - EVal[i]) / ((wRes[p] + energy - EVal[i])*(wRes[p] + energy - EVal[i])+(gammaRes/2.0)*(gammaRes/2.0));
      EValInvIm[i] = -(gammaRes/2.0) / ((wRes[p] + energy - EVal[i])*(wRes[p] + energy - EVal[i])+(gammaRes/2.0)*(gammaRes/2.0));
    }
    
    /* We are trying to calculate Q[0,:]A^(-1)Q^(-1).  Because A is diagonal, we were
     able to calculate B = Q[0,:]*I*Q^(-1) outside the loop. Now we calculate diag(A)*B */
    /* Calculate Q*A^(-1)*Q^(-1)[:,0] */
    GTpsi.NFermion = psi.NFermion;
    GTpsi.NBoson = psi.NBoson;
    GTpsi.QComplex = psi.QComplex;
    GTpsi.NBitsKey = psi.NBitsKey;
    if(WaveFunctionInit(&GTpsi))
    {
      printf("WaveFunctionInit failed in ResonantSpectraOrig 00");
      return 1;
    }
    for(i = 0; i < NKrylov1; i++)
    {
      tempDoubleRe = 0;
      tempDoubleIm = 0;
      
      if(EFun.QComplex)
      {
        for(k = 0; k < NKrylov1; k++)
        {
          tempDoubleRe += EValInvRe[k]*EFun.VRe[k][i]*EFun.VRe[k][0] -
          EValInvIm[k]*EFun.VIm[k][i]*EFun.VRe[k][0] -
          EValInvIm[k]*EFun.VRe[k][i]*EFun.VIm[k][0] -
          EValInvRe[k]*EFun.VIm[k][i]*EFun.VIm[k][0];
          tempDoubleIm += EValInvIm[k]*EFun.VRe[k][i]*EFun.VRe[k][0] +
          EValInvRe[k]*EFun.VIm[k][i]*EFun.VRe[k][0] +
          EValInvRe[k]*EFun.VRe[k][i]*EFun.VIm[k][0] -
          EValInvIm[k]*EFun.VIm[k][i]*EFun.VIm[k][0];
        }
      }
      else
      {
        for(k = 0; k < NKrylov1; k++)
        {
          tempDoubleRe += EValInvRe[k]*EFun.VRe[k][i]*EFun.VRe[k][0];
          tempDoubleIm += EValInvIm[k]*EFun.VRe[k][i]*EFun.VRe[k][0];
        }
      }
      WaveFunctionAddComplexMultipliedUnionBasis(&GTpsi, PsiBasis[i],SpecTri1.b[0]*tempDoubleRe,SpecTri1.b[0]*tempDoubleIm);
    }
    for(i = 0; i < NDec; i++)
    {
      /* Calculate T*G(w)T|i> */
      OperatorPsi(TDec[i],GTpsi,&TGTpsi);
      
      if(rotate)
      {
        printf("\nStarting density matrix calculations and rotations\n"); fflush(stdout);
        /* Calculate density matrix and perform rotations */
        WaveFunctionCalculateNorm(TGTpsi, &tempDoubleRe);
        WaveFunctionNormalize(&TGTpsi);
        
        CalculateDensityMatrix(&Density, TGTpsi);
        MatrixPart(&DensityPart, Density, 16, psi.NFermion-1,16,psi.NFermion-1);
        Eigensystem(DensityPart, &tempMat, DensityValues);
        FreeMatrix(&tempMat);
        
        printf("\n***********************************************\n");
        printf("\n\nEigenvalues of bath part of density matrix:\n\n");
        printf("\n***********************************************\n");
        
        for(k = 0; k < psi.NFermion-16; k++)
          printf("%f\n",DensityValues[k]);
        
        fflush(stdout);
        /* Rotate orbitals to diagonalize bath part of density matrix */
        /* First zero out all off diagonals outside bath block, and shift
         correlated diagonals to -ve so the order is maintained     */
        for(j = 0; j < psi.NFermion; j++)
        {
          for(k = 0; k < psi.NFermion; k++)
          {
            if((j < 16) || (k < 16))
            {
              if(j == k)
              {
                Density.VRe[j][j] = j - 17;
                Density.VIm[j][j] = 0.0;
              }
              else
              {
                Density.VRe[j][k] = 0.0;
                Density.VIm[j][k] = 0.0;
              }
            }
          }
        }
        
        Eigensystem(Density,&TRot,DensityValues);
        MatrixConjugateTranspose(&TRot);
        printf("\nWavefunction size before rotation: %d\n", TGTpsi.N); fflush(stdout);
        //PrintMatrix(TRot);
        
        //rescale by the norm to keep correct intensity
        WaveFunctionMultiply(&TGTpsi, tempDoubleRe, 0.0);
        RotateWaveFunction(&TGTpsiRot, TGTpsi, TRot);
        printf("\nWavefunction size after rotation, before truncation: %d\n", TGTpsi.N); fflush(stdout);
        
        RotateOperatorNormalOrder(&OppH2Rot, OppH2, TRot, &QDetNotOrderedFcFaBaBc);//&QDetNotOrderedFcFaBcBa);
        
        OperatorChop(&OppH2Rot, epsilon);
        WaveFunctionFree(&TGTpsi);
        
        WaveFunctionTruncateBasis(&TGTpsiRot, epsilon);
        printf("\nWavefunction size after rotation, after truncation: %d\n", TGTpsiRot.N); fflush(stdout);
        
        FreeMatrix(&Density);
        FreeMatrix(&DensityPart);
        /* End density matrix calculations and rotations */
      }
      /* Calculate the RIXS spectrum */
      sprintf(spectra[i][p].Name, "%d",i);
      spectra[i][p].N = 2000;
      spectra[i][p].Emin = -5.0;
      spectra[i][p].Emax = 15.0;
      spectra[i][p].dE = (spectra[i][p].Emax - spectra[i][p].Emin) / spectra[i][p].N;
      spectra[i][p].Gamma = gamma;
      SpectrumAllocate(spectra[i][p], &spectra[i][p]);
      SpecTri2.N = NKrylov2;
      InitTriDiagonalMatrix(&SpecTri2);
      
      if(rotate)
        LanczosTriDiagonalize(OppH2Rot, TGTpsiRot, epsilon, &SpecTri2, gamma);
      else
        LanczosTriDiagonalize(OppH2, TGTpsi, epsilon, &SpecTri2, gamma);
      
      TriDiagonalMatrixToSpectrum(energy, SpecTri2, &spectra[i][p]);
      BroadenG(&spectra[i][p], Gfwhm);
      FreeTriDiagonalMatrix(&SpecTri2);
    }
    WaveFunctionFree(&GTpsi);
  }
  return 0;
}

/*  Calculates <i|T_1^dag G^dag T_2 T_2^dag G T_1 |i>
 *  T2 is allowed to be a vector of NDec different operators
 *  input:
 *    OppH opperator defining G
 *    Texc Operator T1
 *    TDec pointer to a list of operators T2
 *    NDec number of elements in TDec
 *    energy of type double, ofset energy of ground-state
 *    epsilon of type double, smallest number kept in determinant
 *    NKrylov1 of type uint32_t number of states in Krylov basis
 *    spectra pointer to NDec spectra (must be allocated, not initialized)
 *      spectra[0].N number of points in spectrum [minus one such that dE = (Emax-Emin)/N]
 *      spectra[0].Emin minimal energy
 *      spectra[0].Emax maximal energy
 *      spectra[0].Gamma lifetime of final state excitation
 *        if spectra[0].N is 0 it is set to 1000
 *        if Emin = - DBL_MAX it is set to Eminpole - 0.125 (Emaxpole-Eminpole)
 *        if Emax =   DBL_MAX it is set to Emaxpole + 0.125 (Emaxpole-Eminpole)
 *        if spectra[0].Gamma is 0 it is set to 10*(Emax-Emin)/N
 *
 *  output:
 *    spectra[0] to spectra[NDec-1] allocated spectra with the FY spectra for TDec[0] to TDec[NDec-1]
 */
uint32_t FluorescenceYieldRestricted(const OperatorType OppH, DeterminantRestrictionType const Restrictions, const OperatorType Texc, const OperatorType *TDec, uint32_t NDec,  const WaveFunctionType psi, double energy, double epsilon, uint32_t NKrylov1, SpectraType *spectra)
{
  /*Given wavefunction psi (usually ground state, but could be excited)*/
  
  /*Calculates resonant spectra using <i|T*GTT*GT|i> in the following steps:
   
   1. Calculate T|i>
   2. Use T|i> to transform H to tridiagonal form in Krylov basis
   3. Invert (w - H_tri + iGamma) to get G(w)
   4. Calculate G(w)T|i>
   5. Calculate T*G(w)T|i>
   6. Calculate the norm squared of T*G(w)T|i>
   */
  
  WaveFunctionType Tpsi,GTpsi,TGTpsi, PsiBasis[NKrylov1];
  TriDiagonalMatrixType SpecTri1;
  MatrixType EFun;
  double EVal[NKrylov1], EValInvRe[NKrylov1], EValInvIm[NKrylov1],tempDoubleRe,tempDoubleIm,wRes;
  int i,p,k;
  /* Calculate T|i> */
  OperatorPsi(Texc,psi,&Tpsi);
  
  /* Use T|i> to transform H to tridiagonal form in Krylov basis */
  /* Using the version which gives back the Krylov basis vectors */
  SpecTri1.N = NKrylov1;
  InitTriDiagonalMatrix(&SpecTri1);
  if(Restrictions.N==0)
  {
    LanczosTriDiagonalizeKrylov(OppH,Tpsi,epsilon,&SpecTri1,PsiBasis);
  }
  else
  {
    LanczosTriDiagonalizeKrylovRestricted(OppH, Tpsi, Restrictions, epsilon, &SpecTri1, PsiBasis);
  }
  /*Compute eigenvalues and eigenvectors, then S^-1 = QAQ^T where
   Q are the eigenvectors and A is a matrix with the eigenvalue
   inverses on the diagonal. */
  if(TriDiagonalMatrixToEigensystem(SpecTri1, &EFun, EVal))
  {
    printf("Problem with TridiagonalMatrixToEigenSystem in FluorescenceYieldRestricted\n");
  }
  /* set spectra values when auto initialization is wanted */
  if(spectra[0].N==0) spectra[0].N=1000;
  spectra[0].EminPole=EVal[0]-energy;
  spectra[0].EmaxPole=EVal[NKrylov1-1]-energy;
  if(spectra[0].Emin== DBL_MAX) spectra[0].Emin = spectra[0].EminPole - 0.125 * (spectra[0].EmaxPole-spectra[0].EminPole);
  if(spectra[0].Emax==-DBL_MAX) spectra[0].Emax = spectra[0].EmaxPole + 0.125 * (spectra[0].EmaxPole-spectra[0].EminPole);
  if(spectra[0].Gamma==0) spectra[0].Gamma=10*(spectra[0].Emax-spectra[0].Emin)/((double)spectra[0].N);
  /* Initialize Spectra */
  spectra[0].dE=(spectra[0].Emax-spectra[0].Emin)/((double)spectra[0].N);
  for(i = 0; i < NDec; i++)
  {
    sprintf(spectra[i].Name, "Spectrum %i/%i",i,NDec);
    SpectrumAllocate(spectra[0], &spectra[i]);
  }
  
  for(p = 0; p <= spectra[0].N; p++)
  {
    /* New method for GT|i>:
     1. Tridiagonalize H using T|i> (outside wres loop)
     2. Calculate eigenvalues (EVal) and eigenfunctions (EFun) of Htri (outside)
     3. For each wres
     a. EValInv[i] = 1/(wres - EVal[i] + iGammaRes)
     b. GT|i> = sum_j EFun[i][j]*EValIn[j]*K[j], where K is Krylov basis vectors
     */
    
    /* Calculate 1/(E_0 + wres - EVal[i] + iGammaRes), split into real and imaginary */
    wRes=spectra[0].Emin+p*spectra[0].dE;
    for(i = 0; i < NKrylov1; i++)
    {
      EValInvRe[i] =  (wRes + energy - EVal[i]) / (Sqr(wRes + energy - EVal[i])+Sqr(spectra[0].Gamma/2.0));
      EValInvIm[i] = -(spectra[0].Gamma/2.0)    / (Sqr(wRes + energy - EVal[i])+Sqr(spectra[0].Gamma/2.0));
    }
    
    /* Calculate Q*A^(-1)*Q^(-1)[:,0] */
    GTpsi.NFermion = psi.NFermion;
    GTpsi.NBoson = psi.NBoson;
    GTpsi.QComplex = psi.QComplex;
    GTpsi.NBitsKey = psi.NBitsKey;
    if(WaveFunctionInit(&GTpsi))
    {
      printf("WaveFunctionInit failed in FluorescenceYieldRestricted 00");
      return 1;
    }
    for(i = 0; i < NKrylov1; i++)
    {
      tempDoubleRe = 0;
      tempDoubleIm = 0;
      
      if(EFun.QComplex)
      {
        for(k = 0; k < NKrylov1; k++)
        {
          tempDoubleRe += EValInvRe[k]*EFun.VRe[k][i]*EFun.VRe[k][0] -
          EValInvIm[k]*EFun.VIm[k][i]*EFun.VRe[k][0] -
          EValInvIm[k]*EFun.VRe[k][i]*EFun.VIm[k][0] -
          EValInvRe[k]*EFun.VIm[k][i]*EFun.VIm[k][0];
          tempDoubleIm += EValInvIm[k]*EFun.VRe[k][i]*EFun.VRe[k][0] +
          EValInvRe[k]*EFun.VIm[k][i]*EFun.VRe[k][0] +
          EValInvRe[k]*EFun.VRe[k][i]*EFun.VIm[k][0] -
          EValInvIm[k]*EFun.VIm[k][i]*EFun.VIm[k][0];
        }
      }
      else
      {
        for(k = 0; k < NKrylov1; k++)
        {
          tempDoubleRe += EValInvRe[k]*EFun.VRe[k][i]*EFun.VRe[k][0];
          tempDoubleIm += EValInvIm[k]*EFun.VRe[k][i]*EFun.VRe[k][0];
        }
      }
      WaveFunctionAddComplexMultipliedUnionBasis(&GTpsi, PsiBasis[i],SpecTri1.b[0]*tempDoubleRe,SpecTri1.b[0]*tempDoubleIm);
    }
    for(i = 0; i < NDec; i++)
    {
      /* Calculate T*G(w)T|i> */
      OperatorPsi(TDec[i],GTpsi,&TGTpsi);
      
      /* Calculate Spectrum */
      WaveFunctionCalculateNorm(TGTpsi, &spectra[i].Re[p]);
      WaveFunctionFree(&TGTpsi);
      /* Want norm squared, only real part is returned, im part is zero. FY is NOT a response function, but actually the spectra correspond to the expectation value of an Hermitian operator */
      spectra[i].Re[p] = Sqr(spectra[i].Re[p]);
      spectra[i].Im[p]=0;
    }
    WaveFunctionFree(&GTpsi);
  }
  // free PsiBasis
  for(i=0;i<NKrylov1;i++)
  {
    WaveFunctionFree(&PsiBasis[i]);
  }
  return 0;
}

uint32_t FluorescenceYield(const OperatorType OppH, const OperatorType Texc, const OperatorType *TDec, uint32_t NDec,  const WaveFunctionType psi, double energy, double epsilon, uint32_t NKrylov1, SpectraType *spectra)
{
  DeterminantRestrictionType Restrictions;
  Restrictions.N=0;
  return FluorescenceYieldRestricted(OppH, Restrictions, Texc, TDec, NDec, psi, energy, epsilon, NKrylov1, spectra);
}


