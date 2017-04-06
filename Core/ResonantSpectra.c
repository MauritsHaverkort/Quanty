//  Functions originated by Robert Green
//

/*  Calculate < psi | T1^+ G1^+ T2^+ G2 T2 G1 T1 | psi >
 */
uint32_t ResonantSpectra(const OperatorType OppH1, const OperatorType OppH2, const OperatorType T1, const OperatorType *T2, uint32_t NT2, const WaveFunctionType psi, TriDiagonalMatrixType **SpecTri2, SpectraType **Spec, ConvergenceSettingsType *settings)
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
  double EVal[settings->NTri1], EValInvRe[settings->NTri1], EValInvIm[settings->NTri1],tempDoubleRe,tempDoubleIm, EVal2[settings->NTri2];
  int i,p,k;
  double energy1, energy2, wRes;
  
  /* calcualte ground state energy */
  PsiHermitianOperatorPsi(OppH1, psi, &energy1);
  PsiHermitianOperatorPsi(OppH2, psi, &energy2);
  
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
    VerbosityPrintf(0x0002, "Automatically set number of energy points (NE=%d) and for intermediate state\n",settings->NE1);
  }
  if(settings->Gamma1==0)
  {
    settings->Gamma1 =(5*(settings->Emax1-settings->Emin1)/settings->NE1);
    VerbosityPrintf(0x0002, "Automatically set lifetime (Gamma=%22.15E) for intermediate state\n",settings->Gamma1);
  }
  
  /* allocate SpecTri2 */
  *SpecTri2 = malloc((settings->NE1+1)*NT2*sizeof(TriDiagonalMatrixType));
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
        for(k = 0; k < settings->NTri1; k++)
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
  *Spec = malloc(NT2 * (settings->NE1+1) * sizeof(SpectraType));
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
    sprintf((*Spec)[i].Name, "Spectrum %i/%i",i,NT2);
    SpectrumAllocate((*Spec)[i], &(*Spec)[i]);
    TriDiagonalMatrixToSpectrum(energy2, (*SpecTri2)[i], &(*Spec)[i]);
  }
  return 0;
}



/*  Calculate < psi | T1^+ G1^+ T2^+ G2 T2 G1 T1 | psi >
 *  psi, T1 and T2 are vectors of Length Np, NT1 and NT2
 *  psi and T1 are contract to a single vector of length Np*NT1
 *  returned is a matrix of dimensions (Np*NT1) * NT2
 */
uint32_t ResonantSpectraTensor(const OperatorType OppH1, const OperatorType OppH2, const OperatorType *T1, uint32_t NT1, const OperatorType *T2, uint32_t NT2, const WaveFunctionType *psi, uint32_t Npsi, TriDiagonalMatrixType **SpecTri2, SpectraType **Spec, double energy1, double energy2, ConvergenceSettingsType *settings)
{
  /*Given wavefunction psi (usually ground state, but could be excited)*/
  
  /*Calculates resonant spectra using <i|T*GTGT*GT|i> in the following steps:
   
   1. Calculate T1|i> for all elments of T1 and |i> giving a vector of length Npsi*NT1
   2. Use T|i> to transform H to block tridiagonal form in Krylov basis
   3. Invert (w - H_tri + iGamma) to get G(w) as a matrix of dimension Npsi*NT1
   4. Calculate G(w)T|i> as a vector of length Npsi*NT1
   5. Calculate T2*G(w)T1|i> as a vector of length (Npsi*NT1) * NT2
   6. Use the vector from step 5 to calculate RIXS spectrum using Block Lanczos
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
  
  //TODO: All
  
  WaveFunctionType Tpsi[Npsi*NT1],GTpsi,TGTpsi, *PsiBasis[settings->NTri1];
  BlockTriDiagonalMatrixType SpecTri1;
  MatrixType EFun;
  double BlockSize[settings->NTri];
  double EVal[settings->NTri1], EValInvRe[settings->NTri1], EValInvIm[settings->NTri1],tempDoubleRe,tempDoubleIm, EVal2[settings->NTri2];
  int i,ipsi,iT1,iT2,p,k;
  double wRes;
  
  /* calcualte ground state energy */
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
  
  /* Calculate T|i> */
  i=0;
  if(settings->restrictions1.N==0)
  {
    for(ipsi=0;ipsi<Npsi;ipsi++)
    {
      for(iT1=0;iT1<NT1;iT1++)
      {
        OperatorPsi(T1[iT1],psi[ipsi],&Tpsi[i]);
        i++;
      }
    }
  }
  else
  {
    for(ipsi=0;ipsi<Npsi;ipsi++)
    {
      for(iT1=0;iT1<NT1;iT1++)
      {
        OperatorPsiRestricted(T1[iT1],psi[ipsi],&Tpsi[i], settings->restrictions1);
        i++;
      }
    }
  }
  
  /* Use T|i>==Tpsi to transform H to BLock tridiagonal form in Krylov basis */
  /* Using the version which gives back the Krylov basis vectors */
  SpecTri1.N = settings->NTri1;
  SpecTri1.BlockSize = Npsi*NT1;
  SpecTri1.QComplex = 0;
  if(OppH1.QComplex!=0) SpecTri1.QComplex = 1;
  for(i=0;i<SpecTri1.BlockSize;i++)
  {
    if(Tpsi[i].QComplex!=0) SpecTri1.QComplex = 1;
  }
  InitBlockTridiagonalMatrix(&SpecTri1);
  if(settings->restrictions1.N==0)
  {
    LanczosBlockTriDiagonalizeKrylov(OppH1, Tpsi, settings->epsilon, &SpecTri1, PsiBasis, BlockSize);
  }
  else
  {
    LanczosBlockTriDiagonalizeKrylovRestricted(OppH1, Tpsi, settings->restrictions1, settings->epsilon, &SpecTri1, PsiBasis, BlockSize);
  }
  
  /*Compute eigenvalues and eigenvectors, then S^-1 = QAQ^T where
   Q are the eigenvectors and A is a matrix with the eigenvalue inverses on the diagonal.
   */
  //TODO: eigensystem of BlockTridiagonal matrix with deflation
  //if(TriDiagonalMatrixToEigensystem(SpecTri1, &EFun, EVal))
  //{
  //  printf("Problem with TridiagonalMatrixToEigenSystem\n");
  //}
  
  
  
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
    VerbosityPrintf(0x0002, "Automatically set number of energy points (NE=%d) and for intermediate state\n",settings->NE1);
  }
  if(settings->Gamma1==0)
  {
    settings->Gamma1 =(5*(settings->Emax1-settings->Emin1)/settings->NE1);
    VerbosityPrintf(0x0002, "Automatically set lifetime (Gamma=%22.15E) for intermediate state\n",settings->Gamma1);
  }
  
  /* allocate SpecTri2 */
  *SpecTri2 = malloc((settings->NE1+1)*NT2*sizeof(TriDiagonalMatrixType));
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
    //GTpsi.NFermion = psi.NFermion;
    //GTpsi.NBoson   = psi.NBoson;
    //GTpsi.QComplex = psi.QComplex;
    //GTpsi.NBitsKey = psi.NBitsKey;
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
        for(k = 0; k < settings->NTri1; k++)
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
        for(k = 0; k < settings->NTri1; k++)
        {
          tempDoubleRe += EValInvRe[k]*EFun.VRe[k][i]*EFun.VRe[k][0];
          tempDoubleIm += EValInvIm[k]*EFun.VRe[k][i]*EFun.VRe[k][0];
        }
      }
      //WaveFunctionAddComplexMultipliedUnionBasis(&GTpsi, PsiBasis[i],SpecTri1.b[0]*tempDoubleRe,SpecTri1.b[0]*tempDoubleIm);
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
  *Spec = malloc(NT2 * (settings->NE1+1) * sizeof(SpectraType));
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
    sprintf((*Spec)[i].Name, "Spectrum %i/%i",i,NT2);
    SpectrumAllocate((*Spec)[i], &(*Spec)[i]);
    TriDiagonalMatrixToSpectrum(energy2, (*SpecTri2)[i], &(*Spec)[i]);
  }
  return 0;
}


/*  Calculate < psi | Texc^+ G TDec[i] G2 TDec[i]^+ G Texc | psi >
 *  Input expects TDec^+
 */
uint32_t ResonantSpectraOrig(const OperatorType OppH, const OperatorType OppH2, const OperatorType Texc, const OperatorType *TDec, uint32_t NDec,  const WaveFunctionType psi, double energy, double epsilon, uint32_t numwRes, SpectraType (*spectra)[numwRes], double *wRes, double gammaRes, double gamma, double Gfwhm, uint32_t NKrylov1, uint32_t NKrylov2)
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


