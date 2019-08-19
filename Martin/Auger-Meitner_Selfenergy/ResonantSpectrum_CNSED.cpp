//  Functions originated by Robert Green
//

/*  Calculate < psi | T1^+ G1^+ T2^+ G2 T2 G1 T1 | psi >
 */
uint32_t ResonantSpectrumCNSED(const OperatorType OppH1, const OperatorType OppH2, const OperatorType T1, const OperatorType *T2,
                               uint32_t NT2, const WaveFunctionType psi, TriDiagonalMatrixType **SpecTri2, SpectraType **Spec,
                               double energy1, double energy2, ConvergenceSettingsType *settings,
                               const vector<vector<OperatorType> > &Uk, const vector<vector<vector<double> > > &Fk, const vector<double> Econt)
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
  
// Martin start	---------------------------------------------------------------------------------------

    double gammaSE = 1;
    uint32_t Npoles1 = settings->NTri1;
    uint32_t Nw1 = settings->NE1+1;
    uint32_t Nw2 = settings->NE2+1;
    uint32_t Npoles2 = settings->NTri2;
    uint32_t Nkappa = (uint32_t) Uk.size();
    vector<vector<double> > polesSE1 = vector<vector<double> >(Npoles1 * Nkappa);
    vector<vector<Complex16> > selfen = vector<vector<Complex16> >(Npoles1);
    vector<WaveFunctionType> EigenBasis1;
    WaveFunctionType UEb[Npoles1];
    RotateWaveFunctions(PsiBasis,Npoles1,EFun,EigenBasis1);
    OperatorType Usum;
    BlockTriDiagonalMatrixType M; 
    M.BlockSize = Npoles1;
    M.N = settings->NTri2; //TODO: implement NTri3
    M.QComplex = EigenBasis1[0].QComplex + OppH1.QComplex; 
    InitBlockTridiagonalMatrix(&M);
    M.mu=0;
    uint32_t BlockSize[M.N];
    WaveFunctionType** KrylBasCont = (WaveFunctionType**) malloc(M.N * sizeof(WaveFunctionType*)); //TODO: free
    vector<vector<Complex16> > SelfEnergy, SEtemp;
    for (uint32_t i = 0; i < (uint32_t) Uk.size(); i++) // loop over kappa
    {
        InitOperator(&Usum);
        for (uint32_t k = 0; k < (uint32_t) Uk[i].size(); k++) // sum over slater integrals k
        {
                OperatorAdd(&Usum,Uk[i][k]);
        }
        for (uint32_t j = 0; j < Npoles1; j++)
        {
            UEb[j].NFermion = EigenBasis1[0].NFermion;
            UEb[j].NBoson = EigenBasis1[0].NBoson;
            UEb[j].QComplex = EigenBasis1[0].QComplex + Usum.QComplex; 
            UEb[j].NBitsKey = EigenBasis1[0].NBitsKey;
            WaveFunctionInit(&(UEb[j]));
            OperatorPsi(Usum, EigenBasis1[j], &(UEb[j]) );
        }
        LanczosBlockTriDiagonalizeKrylov(OppH1, UEb, settings->epsilon,&M, KrylBasCont, BlockSize);
        vector<double> EAuger = vector<double>(M.N*M.BlockSize);
        CompactMatrixType fun;  // Should not be allocated acording to Documentation of CompactMatrixMultiply called in BlockTriDiagonalMatrixEigenSystem
        BlockTriDiagonalMatrixEigenSystem(&M, &(EAuger[0]), &fun);
        double wMin = settings->Emin1 + energy1;                    //TODO: What is this shift? is it needed here or does it have to be removed?
        double wMax = settings->Emax1 + energy1;
        double dw   = (settings->Emax1-settings->Emin1)/((double)settings->NE1);
        vector<WaveFunctionType> psiCont;
        RotateWaveFunctions(KrylBasCont,M.N,M.BlockSize,fun,psiCont);
        CalcCNSE(psiCont, EigenBasis1, Uk[i], wMin, wMax, dw, gammaSE, Fk[i], Econt, EAuger, SEtemp);  
        
        /*for (uint32_t j = 0; j < Npoles; j++)
        {
            WaveFunctionFree(&UEb[j]);  // Documentation of LanczosBlockTriDiagonalizeKrylov says that these are already freeed
        }*/
        FreeOperator(&Usum);
    }
// Martin end	---------------------------------------------------------------------------------------  
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



