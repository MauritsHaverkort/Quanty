//  Functions originated by Robert Green
//

/*  Calculate < psi | T1^+ G1^+ T2^+ G2 T2 G1 T1 | psi >
 */
uint32_t ResonantSpectrumCNSED(const OperatorType OppH1, const OperatorType OppH2, const OperatorType T1, const OperatorType *T2,
                               uint32_t NT2, const WaveFunctionType psi, TriDiagonalMatrixType **SpecTri2, SpectraType **Spec,
                               double energy1, double energy2, ConvergenceSettingsType *settings, const vector<vector<OperatorType> > &Uk,
                               const vector<vector<vector<double> > > &Fk, const vector<vector<double> > &Econt, double gammaSE)
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
  //PrintWaveFunctionReadable(Tpsi);
  
  /* Use T|i> to transform H to tridiagonal form in Krylov basis */
  /* Using the version which gives back the Krylov basis vectors */
  SpecTri1.N = settings->NTri1;
  InitTriDiagonalMatrix(&SpecTri1);
  if(settings->restrictions1.N==0)
  {
    LanczosTriDiagonalizeKrylov(OppH1,Tpsi,settings->epsilon,&SpecTri1,PsiBasis);
    //PrintWaveFunctionReadable(PsiBasis[0]);
    //printf("settings->epsilon = %E\n",settings->epsilon);
  }
  else
  {
    LanczosTriDiagonalizeKrylovRestricted(OppH1, Tpsi, settings->restrictions1, settings->epsilon, &SpecTri1, PsiBasis);
  }
  /*
  printf("a\tb\n");
  for (uint32_t j = 0; j < SpecTri1.N; j++)
  {
      printf("%E\t%E\n",SpecTri1.a[j],SpecTri1.b[j]);
  }*/
  //PrintWaveFunctionReadable(psi);
  //PrintOperatorReadable(T1);
    
  /*Compute eigenvalues and eigenvectors, then S^-1 = QAQ^T where
   Q are the eigenvectors and A is a matrix with the eigenvalue inverses on the diagonal.
   */
  if(TriDiagonalMatrixToEigensystem(SpecTri1, &EFun, EVal))
  {
    printf("Problem with TridiagonalMatrixToEigenSystem\n");
  }
  //for (uint32_t j = 0; j < settings->NTri1; j++) printf("%E\t",EVal[j]);printf("\n#Eval = %u\n",settings->NTri1);
  
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

    uint32_t Npoles1 = settings->NTri1;
    uint32_t Nw1 = settings->NE1+1;
    // mwh not used uint32_t Nw2 = settings->NE2+1;
    uint32_t Nkappa = (uint32_t) Uk.size();
    vector<vector<double> > polesSE1 = vector<vector<double> >(Npoles1 * Nkappa);
    vector<vector<Complex16> > selfen = vector<vector<Complex16> >(Npoles1);
    vector<WaveFunctionType> EigenBasis1;
    WaveFunctionType UEb[Npoles1];
    RotateWaveFunctions(PsiBasis,Npoles1,EFun,EigenBasis1);
    //PrintWaveFunctionReadable(EigenBasis1[0]);
    
    OperatorType Usum;
    
    vector<vector<Complex16> > SelfEnergy, SEtemp;
    SelfEnergy = vector<vector<Complex16> >(Npoles1);
    printf("energy1 = %E\n",energy1);
    double wMin = settings->Emin1 ;//+ energy1;                    //TODO: What is this shift? is it needed here or does it have to be removed?
    double wMax = settings->Emax1 ;//+ energy1;
    double dw   = (settings->Emax1-settings->Emin1)/((double)settings->NE1);
    if (gammaSE < 0 ) 
    {
        if (Econt[0].size() > 1)
        {
            gammaSE = 3 * (Econt[0][1] - Econt[0][0]);
        }
        else
        {
            gammaSE = settings->Gamma2;
        }
    }
        
    for (uint32_t i = 0; i < Npoles1; i++)
    {
        SelfEnergy[i] = vector<Complex16>(Nw1);
    }
    for (uint32_t i = 0; i < (uint32_t) Uk.size(); i++) // loop over kappa
    {
        OperatorCopy(&Usum,Uk[i][0]);
        for (uint32_t k = 1; k < (uint32_t) Uk[i].size(); k++) // sum over slater integrals k
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
            WaveFunctionTruncateBasis(&(UEb[j]),1E-14);
            //PrintWaveFunctionReadable(UEb[j]);
        }
        BlockTriDiagonalMatrixType M; 
        M.BlockSize = Npoles1; // -- ************ WAAAAAARRRRRRRRNIIIIIIIIIIIIIING ***-++++++++++++++++++++++---/+-
        M.N = settings->NTri+1; 
        M.QComplex = EigenBasis1[0].QComplex + OppH1.QComplex; 
        InitBlockTridiagonalMatrix(&M);
        M.mu=0;
        vector<uint32_t> BlockSize = vector<uint32_t>(M.N);
        WaveFunctionType** KrylBasCont = (WaveFunctionType**) malloc(M.N * sizeof(WaveFunctionType*)); //TODO: free
        LanczosBlockTriDiagonalizeKrylov(OppH1, UEb, settings->epsilon,&M, KrylBasCont, &(BlockSize[0]));
        M.N--;
        BlockSize.pop_back();
        //PrintWaveFunctionReadable(KrylBasCont[0][0]);
/*        printf("a:\n");
        for (uint32_t k = 1; k < M.N; k++)
        {
            for (uint32_t m = 0; m < M.BlockSize; m++)
            {
                for (uint32_t n = 0; n < M.BlockSize; n++)
                {
                    printf("%E\t",M.aRe[k*M.BlockSize*M.BlockSize + m*M.BlockSize +n]);
                }
                printf("\n");
            }
            printf("\n");
        }
        printf("b:\n");
        for (uint32_t k = 1; k < M.N-1; k++)
        {
            for (uint32_t m = 0; m < M.BlockSize; m++)
            {
                for (uint32_t n = 0; n < M.BlockSize; n++)
                {
                    printf("%E\t",M.bRe[k*M.BlockSize*M.BlockSize + m*M.BlockSize +n]);
                }
                printf("\n");
            }
            printf("\n");
        }*/
        CompactMatrixType CM;
      printf("\n\n\nmwh test if the chagne I made is ok\n\n\n");
      // mwh this function is double
      // BlockTriMatrixToCompact(M,BlockSize,CM);
      // mwh please use
      // and test if BlockSize is set correctly in M
      BlockTriDiagonalMatrixToCompactMatrixWithDeflation(&M, &CM);
        /*printf("*************************************************************************\n");
        PrintMartix(CM.VRe,CM.Ni);
        printf("*************************************************************************\n");*/
        vector<double> EAuger = vector<double>(CM.Ni);
        CompactMatrixType fun;  // Should not be allocated acording to Documentation of CompactMatrixMultiply called in BlockTriDiagonalMatrixEigenSystem
        BlockTriDiagonalMatrixEigenSystem(&M, &(EAuger[0]), &fun);
        for (uint32_t j = 0; j < (uint32_t) EAuger.size(); j++) printf("%E\t",EAuger[j]);printf("\n#Eauger = %u\n", (uint32_t)EAuger.size());
        EAuger[0] = 0;
        CompactMatrixEigensystem(&CM,&fun,&(EAuger[0]));
        for (uint32_t j = 0; j < (uint32_t) EAuger.size(); j++) printf("%E\t",EAuger[j]);printf("\n#Eauger = %u\n",(uint32_t)EAuger.size());
        FreeCompactMatrix(&CM);
        vector<WaveFunctionType> psiCont;
        RotateWaveFunctions(KrylBasCont,M.N,M.BlockSize,BlockSize,fun,psiCont);
        
        //PrintWaveFunctionReadable(psiCont[0]);
 /*       for (uint32_t k = 0; k < M.N; k++)
        {
            for (uint32_t j = 0; j < BlockSize[k]; j++)
                WaveFunctionFree(&(KrylBasCont[k])[j]);
            //free(KrylBasCont[k];
        }
        free(KrylBasCont);*/
        FreeBlockTridiagonalMatrix(&M);
        /*         printf("psiCont.size = %u, fun.Ni = %u\n",psiCont.size(),fun.Ni);
        for (uint32_t fu = 0; fu < EigenBasis1.size(); fu++)
        {
            for ( uint32_t ba = 0; ba < psiCont.size(); ba++)
                printf("%E\t",fun.VRe[fu*fun.Nj +ba]);
            printf("\n");
            PrintWaveFunctionReadable(EigenBasis1[fu]);
        }
        printf("\n");
 
        // Debug below here -------------------****************************---------------------------------------********************************-----------------<<<<<<<<<<<<<<<<<<<<
*/      printf("wMin = %E\twMax = %E\tdw = %E\t gamma = %E\n",wMin,wMax,dw,gammaSE);
        CalcCNSE(psiCont, EigenBasis1, Uk[i], wMin, wMax, dw, gammaSE, Fk[i], Econt[i], EAuger, SEtemp); 
        for (uint32_t j = 0; j < Npoles1; j++)
        {
            //printf("%E\t%E\n",SEtemp[j][0].Re,SEtemp[j][0].Im);
            AddVectors(SelfEnergy[j],SEtemp[j]);
        }        
        /*for (uint32_t j = 0; j < Npoles; j++)
        {
            WaveFunctionFree(&UEb[j]);  // Documentation of LanczosBlockTriDiagonalizeKrylov says that these are already freeed
        }*/
        FreeOperator(&Usum);
    }
    printf("%E\t%E\n",SelfEnergy[0][0].Re,SelfEnergy[0][0].Im);
    /*for (uint32_t i = 0; i < Npoles1; i++)
    {
        printf("%E\t%E\n",SelfEnergy[i][0].Re,SelfEnergy[i][0].Im);
    }*/
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
      //EValInvRe[i] =  (wRes - EVal[i])       / ((wRes - EVal[i])*(wRes - EVal[i])+(settings->Gamma1/2.0)*(settings->Gamma1/2.0));
      //EValInvIm[i] = -(settings->Gamma1/2.0) / ((wRes - EVal[i])*(wRes - EVal[i])+(settings->Gamma1/2.0)*(settings->Gamma1/2.0));
      EValInvRe[i] =  (wRes - EVal[i])       / ((wRes - EVal[i])*(wRes - EVal[i]) + SelfEnergy[i][p].Im*SelfEnergy[i][p].Im);
      //EValInvIm[i] = (SelfEnergy[i][p].Im - settings->Gamma1/2.0)/ ((wRes - EVal[i])*(wRes - EVal[i]) + (SelfEnergy[i][p].Im - settings->Gamma1/2.0)*(SelfEnergy[i][p].Im - settings->Gamma1/2.0));
      EValInvIm[i] = (SelfEnergy[i][p].Im )/ ((wRes - EVal[i])*(wRes - EVal[i]) + (SelfEnergy[i][p].Im )*(SelfEnergy[i][p].Im));
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

uint32_t ResonantSpectrumCNSED(const OperatorType OppH1, const OperatorType OppH2, const OperatorType T1, const OperatorType *T2,
                               uint32_t NT2, const WaveFunctionType psi, TriDiagonalMatrixType **SpecTri2, SpectraType **Spec,
                               double energy1, double energy2, ConvergenceSettingsType *settings, const vector<vector<Complex16> > &SelfEnergy)
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
  //PrintWaveFunctionReadable(Tpsi);
  
  /* Use T|i> to transform H to tridiagonal form in Krylov basis */
  /* Using the version which gives back the Krylov basis vectors */
  SpecTri1.N = settings->NTri1;
  InitTriDiagonalMatrix(&SpecTri1);
  if(settings->restrictions1.N==0)
  {
    LanczosTriDiagonalizeKrylov(OppH1,Tpsi,settings->epsilon,&SpecTri1,PsiBasis);
    //PrintWaveFunctionReadable(PsiBasis[0]);
    //printf("settings->epsilon = %E\n",settings->epsilon);
  }
  else
  {
    LanczosTriDiagonalizeKrylovRestricted(OppH1, Tpsi, settings->restrictions1, settings->epsilon, &SpecTri1, PsiBasis);
  }
  /*
  printf("a\tb\n");
  for (uint32_t j = 0; j < SpecTri1.N; j++)
  {
      printf("%E\t%E\n",SpecTri1.a[j],SpecTri1.b[j]);
  }*/
  //PrintWaveFunctionReadable(psi);
  //PrintOperatorReadable(T1);
    
  /*Compute eigenvalues and eigenvectors, then S^-1 = QAQ^T where
   Q are the eigenvectors and A is a matrix with the eigenvalue inverses on the diagonal.
   */
  if(TriDiagonalMatrixToEigensystem(SpecTri1, &EFun, EVal))
  {
    printf("Problem with TridiagonalMatrixToEigenSystem\n");
  }
  //for (uint32_t j = 0; j < settings->NTri1; j++) printf("%E\t",EVal[j]);printf("\n#Eval = %u\n",settings->NTri1);
  
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
      //EValInvRe[i] =  (wRes - EVal[i])       / ((wRes - EVal[i])*(wRes - EVal[i])+(settings->Gamma1/2.0)*(settings->Gamma1/2.0));
      //EValInvIm[i] = -(settings->Gamma1/2.0) / ((wRes - EVal[i])*(wRes - EVal[i])+(settings->Gamma1/2.0)*(settings->Gamma1/2.0));
      EValInvRe[i] =  (wRes - EVal[i] - SelfEnergy[i][p].Re)       / ((wRes - EVal[i]- SelfEnergy[i][p].Re)*(wRes - EVal[i]- SelfEnergy[i][p].Re) + SelfEnergy[i][p].Im*SelfEnergy[i][p].Im);
      //EValInvIm[i] = (SelfEnergy[i][p].Im - settings->Gamma1/2.0)/ ((wRes - EVal[i])*(wRes - EVal[i]) + (SelfEnergy[i][p].Im - settings->Gamma1/2.0)*(SelfEnergy[i][p].Im - settings->Gamma1/2.0));
      EValInvIm[i] = (SelfEnergy[i][p].Im )/ ((wRes - EVal[i]- SelfEnergy[i][p].Re)*(wRes - EVal[i]- SelfEnergy[i][p].Re) + (SelfEnergy[i][p].Im )*(SelfEnergy[i][p].Im));
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

static int LuaResSpectCNSED(lua_State *L)
{
    struct
  {
    ConvergenceSettingsType Settings;
    double E0, E01, E02;
    uint32_t Npsi;
    uint32_t NT1;
    uint32_t NT2;
    double GammaSE;
    bool Tensor;
    bool DenseMethod;
  } Options;
  Options.Settings=InitConvergenceSettings();
  Options.E0=DBL_MAX;
  Options.E01=DBL_MAX;
  Options.E02=DBL_MAX;
  Options.Tensor=false;
  Options.DenseMethod=false;
  Options.Settings.NTri = 5;
  Options.GammaSE = -1;
  ConvergenceSettingsType SettingsNow;
  SpectraTypeArray *Spec;
  SpectraType *SpecT1psi;
  TriDiagonalMatrixType *SpecTri;
  OperatorType *H1,*H2,*T1, *T2;
  WaveFunctionType *psi; //mwh not used HTpsiDenseBasis,HTHTpsiDenseBasis,Opsi,TpsiSum;
  uint32_t iT,ipsi,len,i;
  LuaOptionsType LuaCreateResonantSpectraOptions[]={
    {"NTri1","uint32_t",&(Options.Settings.NTri1)},
    {"NTri2","uint32_t",&(Options.Settings.NTri2)},
    {"Tensor","bool",&(Options.Tensor)},
    {"restrictions1", "DeterminantRestrictionType", &(Options.Settings.restrictions1)},
    {"restrictions2", "DeterminantRestrictionType", &(Options.Settings.restrictions2)},
    {"Restrictions1", "DeterminantRestrictionType", &(Options.Settings.restrictions1)},
    {"Restrictions2", "DeterminantRestrictionType", &(Options.Settings.restrictions2)},
    {"epsilon","double",&(Options.Settings.epsilon)},
    {"Epsilon","double",&(Options.Settings.epsilon)},
    {"E0","double",&(Options.E0)},
    {"E0_1","double",&(Options.E01)},
    {"E0_2","double",&(Options.E02)},
    {"Emin1","double",&(Options.Settings.Emin1)},
    {"Emax1","double",&(Options.Settings.Emax1)},
    {"NE1","uint32_t",&(Options.Settings.NE1)},
    {"Gamma1","double",&(Options.Settings.Gamma1)},
    {"Emin2","double",&(Options.Settings.Emin2)},
    {"Emax2","double",&(Options.Settings.Emax2)},
    {"NE2","uint32_t",&(Options.Settings.NE2)},
    {"Gamma2","double",&(Options.Settings.Gamma2)},
    {"DenseBorder","uint32_t",&(Options.Settings.DenseBorder)},
    {"ForceDenseMethod","uint32_t",&(Options.DenseMethod)},
	{"NTriSE","uint32_t",&(Options.Settings.NTri)},
	{"GammaSE","double",&(Options.GammaSE)},
    {NULL, NULL, NULL}};
  printf("Start ResonantSpectrum\n");
  /* free garbage full cycle*/
  if(L->l_G->gcrunning) lua_gc(L, LUA_GCCOLLECT, 0);
  /* check number of input parameters */
  len = (uint32_t)lua_gettop(L);
  if(!((len==9)||(len==8)))
  {
    return luaL_error(L, "Error LuaCreateResonantSpectra called with %d arguments,8 or 9 expected (H1, H2, T1, T2^+, psi,Uk,Fk,EnergyContinuum, {options})\n",len);
  }
  /* read operator and wavefunctions */
  H1 = (OperatorType *) luaL_checkudata(L, 1, "Operator_Type");
  H2 = (OperatorType *) luaL_checkudata(L, 2, "Operator_Type");
  if(lua_type(L,3)==LUA_TTABLE)
  {
    T1=luaL_checkarray_OperatorType(L, 3, &Options.NT1);
    if(Options.NT1==0) return luaL_error(L, "Errof CreateResonantSpectra called with a table of transition operators (1) of length 0\n");
  }
  else
  {
    Options.NT1 = 1;
    T1 = (OperatorType *) luaL_checkudata(L, 3, "Operator_Type");
  }
  if(lua_type(L,4)==LUA_TTABLE)
  {
    T2=luaL_checkarray_OperatorType(L, 4, &Options.NT2);
    if(Options.NT2==0) return luaL_error(L, "Errof CreateResonantSpectra called with a table of transition operators (2) of length 0\n");
  }
  else
  {
    Options.NT2 = 1;
    T2 = (OperatorType *) luaL_checkudata(L, 4, "Operator_Type");
  }
  if(lua_type(L,5)==LUA_TTABLE)
  {
    psi=luaL_checkarray_WaveFunctionType(L, 5, &Options.Npsi);
    if(Options.Npsi==0) return luaL_error(L, "Errof CreateResonantSpectra called with a table of wave functions of length 0\n");
  }
  else
  {
    Options.Npsi = 1;
    psi = (WaveFunctionType *) luaL_checkudata(L, 5, "Wavefunction_Type");
  }
  vector<vector<OperatorType> > Uk;
  vector<vector<vector<double> > > Fk;
  lua_check(L,6,Uk);
  vector<vector<double> > EnergyContinuum;
  lua_check(L,7,Fk);
  lua_check(L,8,EnergyContinuum);

  /* read options */
  if (len==9)
  {
    LuaCheckOptions(L, -1, LuaCreateResonantSpectraOptions);
  }
  if(Options.E01==DBL_MAX)Options.E01=Options.E0;
  if(Options.E02==DBL_MAX)Options.E02=Options.E0;
  /* test if operator is Hermitian */
  switch(QOperatorHermitian(*H1, 0.0001*Options.Settings.zero))
  {
    case 3 :
    {
      printf("WARNING: QOperatorHermitian failed in LuaCreateResonantSpectra, this should not happen\n");
      printf("my guess is that you are out of memory and will crash soon\n");
      break;
    }
    case 2 :
    {
      luaL_warning(L, "CreateResonantSpectra called with an operator (H1) that is not fully Hermitian\nStrongly sugested to stop the calculation and repair this\n\n\n\nDo it\n\n\nReally Check this\n\n\nDid you?\n\n\n");
      break;
    }
    case 1 :
    {
      // all fine H is Hermitian upto zero
      break;
    }
    case 0 :
    {
      luaL_error(L, "CreateResonantSpectra called with an operator (H1) that is not Hermitian\n");
      break;
    }
  }
  switch(QOperatorHermitian(*H2, 0.0001*Options.Settings.zero))
  {
    case 3 :
    {
      printf("WARNING: QOperatorHermitian failed in LuaCreateResonantSpectra, this should not happen\n");
      printf("my guess is that you are out of memory and will crash soon\n");
      break;
    }
    case 2 :
    {
      luaL_warning(L, "CreateResonantSpectra called with an operator (H2) that is not fully Hermitian\n");
      break;
    }
    case 1 :
    {
      // all fine H is Hermitian upto zero
      break;
    }
    case 0 :
    {
      luaL_error(L, "CreateResonantSpectra called with an operator (H2) that is not Hermitian\n");
      break;
    }
  }
  /* allocate spectra and a pointer to return matrix on lua stack */
  Spec = (SpectraTypeArray *) lua_newuserdata(L, sizeof(SpectraTypeArray));
  luaL_setmetatable(L, "Spectra_Type");
  Spec->NSpectra=0;
  Spec->Spec=NULL;

    if(Options.Tensor==false)
    {
        for(ipsi=0;ipsi<Options.Npsi;ipsi++)
        {
        for(iT=0;iT<Options.NT1;iT++)
        {
            SettingsNow = Options.Settings;
            ResonantSpectrumCNSED(*H1, *H2, T1[iT], T2, Options.NT2, psi[ipsi], &SpecTri, &SpecT1psi, Options.E01, Options.E02, &(SettingsNow), Uk, Fk, EnergyContinuum,Options.GammaSE);
            Spec->NSpectra+=(SettingsNow.NE1+1)*Options.NT2;
            Spec->Spec= (SpectraType*)realloc(Spec->Spec, Spec->NSpectra*sizeof(SpectraType));
            for(i=0;i<(SettingsNow.NE1+1)*Options.NT2;i++)
            {
            FreeTriDiagonalMatrix(&(SpecTri[i]));
            }
            free(SpecTri);
            for(i=0;i<(SettingsNow.NE1+1)*Options.NT2;i++)
            {
            sprintf(SpecT1psi[i].Name, "%s T1 %i/%i psi %i/%i",SpecT1psi[i].Name,iT,Options.NT1,ipsi,Options.Npsi);
            Spec->Spec[i+Spec->NSpectra-(SettingsNow.NE1+1)*Options.NT2] = SpecT1psi[i];
            }
            free(SpecT1psi);
        }
        }
    }
    else
    {

        //ResonantSpectraTensor(*H1, *H2, T1, Options.NT1, T2, Options.NT2, psi, Options.Npsi, &SpecTri, &(Spec->Spec[0]), Options.E01, Options.E02, &Options.Settings);
        return luaL_error(L,"Tensor form for Resonant spectra not yet available");
    }

  return 1;
    
}


static int LuaResSpectCNSEDInput(lua_State *L)
{
    struct
  {
    ConvergenceSettingsType Settings;
    double E0, E01, E02;
    uint32_t Npsi;
    uint32_t NT1;
    uint32_t NT2;
    double GammaSE;
    bool Tensor;
    bool DenseMethod;
  } Options;
  Options.Settings=InitConvergenceSettings();
  Options.E0=DBL_MAX;
  Options.E01=DBL_MAX;
  Options.E02=DBL_MAX;
  Options.Tensor=false;
  Options.DenseMethod=false;
  Options.Settings.NTri = 5;
  Options.GammaSE = -1;
  ConvergenceSettingsType SettingsNow;
  SpectraTypeArray *Spec;
  SpectraType *SpecT1psi;
  TriDiagonalMatrixType *SpecTri;
  OperatorType *H1,*H2,*T1, *T2;
  WaveFunctionType *psi; // mwh not used, HTpsiDenseBasis,HTHTpsiDenseBasis,Opsi,TpsiSum;
  uint32_t iT,ipsi,len,i;
  LuaOptionsType LuaCreateResonantSpectraOptions[]={
    {"NTri1","uint32_t",&(Options.Settings.NTri1)},
    {"NTri2","uint32_t",&(Options.Settings.NTri2)},
    {"Tensor","bool",&(Options.Tensor)},
    {"restrictions1", "DeterminantRestrictionType", &(Options.Settings.restrictions1)},
    {"restrictions2", "DeterminantRestrictionType", &(Options.Settings.restrictions2)},
    {"Restrictions1", "DeterminantRestrictionType", &(Options.Settings.restrictions1)},
    {"Restrictions2", "DeterminantRestrictionType", &(Options.Settings.restrictions2)},
    {"epsilon","double",&(Options.Settings.epsilon)},
    {"Epsilon","double",&(Options.Settings.epsilon)},
    {"E0","double",&(Options.E0)},
    {"E0_1","double",&(Options.E01)},
    {"E0_2","double",&(Options.E02)},
    {"Emin1","double",&(Options.Settings.Emin1)},
    {"Emax1","double",&(Options.Settings.Emax1)},
    {"NE1","uint32_t",&(Options.Settings.NE1)},
    {"Gamma1","double",&(Options.Settings.Gamma1)},
    {"Emin2","double",&(Options.Settings.Emin2)},
    {"Emax2","double",&(Options.Settings.Emax2)},
    {"NE2","uint32_t",&(Options.Settings.NE2)},
    {"Gamma2","double",&(Options.Settings.Gamma2)},
    {"DenseBorder","uint32_t",&(Options.Settings.DenseBorder)},
    {"ForceDenseMethod","uint32_t",&(Options.DenseMethod)},
	{"NTriSE","uint32_t",&(Options.Settings.NTri)},
	{"GammaSE","double",&(Options.GammaSE)},
    {NULL, NULL, NULL}};
  printf("Start ResonantSpectrum\n");
  /* free garbage full cycle*/
  if(L->l_G->gcrunning) lua_gc(L, LUA_GCCOLLECT, 0);
  /* check number of input parameters */
  len = (uint32_t)lua_gettop(L);
  if(!((len==6)||(len==7)))
  {
    return luaL_error(L, "Error LuaCreateResonantSpectra called with %d arguments,6 or 7 expected (H1, H2, T1, T2^+, psi,Uk,Fk,EnergyContinuum, {options})\n",len);
  }
  /* read operator and wavefunctions */
  H1 = (OperatorType *) luaL_checkudata(L, 1, "Operator_Type");
  H2 = (OperatorType *) luaL_checkudata(L, 2, "Operator_Type");
  if(lua_type(L,3)==LUA_TTABLE)
  {
    T1=luaL_checkarray_OperatorType(L, 3, &Options.NT1);
    if(Options.NT1==0) return luaL_error(L, "Errof CreateResonantSpectra called with a table of transition operators (1) of length 0\n");
  }
  else
  {
    Options.NT1 = 1;
    T1 = (OperatorType *) luaL_checkudata(L, 3, "Operator_Type");
  }
  if(lua_type(L,4)==LUA_TTABLE)
  {
    T2=luaL_checkarray_OperatorType(L, 4, &Options.NT2);
    if(Options.NT2==0) return luaL_error(L, "Errof CreateResonantSpectra called with a table of transition operators (2) of length 0\n");
  }
  else
  {
    Options.NT2 = 1;
    T2 = (OperatorType *) luaL_checkudata(L, 4, "Operator_Type");
  }
  if(lua_type(L,5)==LUA_TTABLE)
  {
    psi=luaL_checkarray_WaveFunctionType(L, 5, &Options.Npsi);
    if(Options.Npsi==0) return luaL_error(L, "Errof CreateResonantSpectra called with a table of wave functions of length 0\n");
  }
  else
  {
    Options.Npsi = 1;
    psi = (WaveFunctionType *) luaL_checkudata(L, 5, "Wavefunction_Type");
  }
  vector<vector<Complex16> > selfen;
  lua_check(L,6,selfen);
  /* read options */
  if (len==7)
  {
    LuaCheckOptions(L, -1, LuaCreateResonantSpectraOptions);
  }
  if(Options.E01==DBL_MAX)Options.E01=Options.E0;
  if(Options.E02==DBL_MAX)Options.E02=Options.E0;
  /* test if operator is Hermitian */
  switch(QOperatorHermitian(*H1, 0.0001*Options.Settings.zero))
  {
    case 3 :
    {
      printf("WARNING: QOperatorHermitian failed in LuaCreateResonantSpectra, this should not happen\n");
      printf("my guess is that you are out of memory and will crash soon\n");
      break;
    }
    case 2 :
    {
      luaL_warning(L, "CreateResonantSpectra called with an operator (H1) that is not fully Hermitian\nStrongly sugested to stop the calculation and repair this\n\n\n\nDo it\n\n\nReally Check this\n\n\nDid you?\n\n\n");
      break;
    }
    case 1 :
    {
      // all fine H is Hermitian upto zero
      break;
    }
    case 0 :
    {
      luaL_error(L, "CreateResonantSpectra called with an operator (H1) that is not Hermitian\n");
      break;
    }
  }
  switch(QOperatorHermitian(*H2, 0.0001*Options.Settings.zero))
  {
    case 3 :
    {
      printf("WARNING: QOperatorHermitian failed in LuaCreateResonantSpectra, this should not happen\n");
      printf("my guess is that you are out of memory and will crash soon\n");
      break;
    }
    case 2 :
    {
      luaL_warning(L, "CreateResonantSpectra called with an operator (H2) that is not fully Hermitian\n");
      break;
    }
    case 1 :
    {
      // all fine H is Hermitian upto zero
      break;
    }
    case 0 :
    {
      luaL_error(L, "CreateResonantSpectra called with an operator (H2) that is not Hermitian\n");
      break;
    }
  }
  /* allocate spectra and a pointer to return matrix on lua stack */
  Spec = (SpectraTypeArray *) lua_newuserdata(L, sizeof(SpectraTypeArray));
  luaL_setmetatable(L, "Spectra_Type");
  Spec->NSpectra=0;
  Spec->Spec=NULL;

    if(Options.Tensor==false)
    {
        for(ipsi=0;ipsi<Options.Npsi;ipsi++)
        {
        for(iT=0;iT<Options.NT1;iT++)
        {
            SettingsNow = Options.Settings;
            ResonantSpectrumCNSED(*H1, *H2, T1[iT], T2, Options.NT2, psi[ipsi], &SpecTri, &SpecT1psi, Options.E01, Options.E02, &(SettingsNow), selfen);
            Spec->NSpectra+=(SettingsNow.NE1+1)*Options.NT2;
            Spec->Spec= (SpectraType*)realloc(Spec->Spec, Spec->NSpectra*sizeof(SpectraType));
            for(i=0;i<(SettingsNow.NE1+1)*Options.NT2;i++)
            {
            FreeTriDiagonalMatrix(&(SpecTri[i]));
            }
            free(SpecTri);
            for(i=0;i<(SettingsNow.NE1+1)*Options.NT2;i++)
            {
            sprintf(SpecT1psi[i].Name, "%s T1 %i/%i psi %i/%i",SpecT1psi[i].Name,iT,Options.NT1,ipsi,Options.Npsi);
            Spec->Spec[i+Spec->NSpectra-(SettingsNow.NE1+1)*Options.NT2] = SpecT1psi[i];
            }
            free(SpecT1psi);
        }
        }
    }
    else
    {

        //ResonantSpectraTensor(*H1, *H2, T1, Options.NT1, T2, Options.NT2, psi, Options.Npsi, &SpecTri, &(Spec->Spec[0]), Options.E01, Options.E02, &Options.Settings);
        return luaL_error(L,"Tensor form for Resonant spectra not yet available");
    }

  return 1;
    
}





































