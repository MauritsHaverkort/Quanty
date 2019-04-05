// Function copy pasetd from REsonantSpectra.cpp and self energy added

uint32_t ResonantSpectraSelfEnergy(const OperatorType OppH1, const OperatorType OppH2, const OperatorType T1, const OperatorType *T2, const vector<OperatorType*> &U, uint32_t NT2, const WaveFunctionType psi, TriDiagonalMatrixType **SpecTri2, SpectraType **Spec, double energy1, double energy2, ConvergenceSettingsType *settings, const vector<double> &EnergyContinuum, uint32_t NTriSE)
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

//settings->NTri1 = 4; // TODO Remove that ---------------------------------------------------------------------------------------------------
  WaveFunctionType Tpsi,GTpsi,TGTpsi, PsiBasis[settings->NTri1],*EigenBasis;
  TriDiagonalMatrixType SpecTri1;
  MatrixType EFun;
  double EVal[settings->NTri1], EValInvRe[settings->NTri1], EValInvIm[settings->NTri1],tempDoubleRe,tempDoubleIm, EVal2[settings->NTri2], Norm;
  int i,p,k;
  double wRes;
  uint32_t NSelfEnergy = (uint32_t) U.size();

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

  /* 1) Rotate PsiBasis to EigenBasis
   * 2) Calculate SelfEnergy SE on EigenBasis in tri-diagonal form
   */

  printf("NumberOf stuff : %d\nEV[0] = %f\n",settings->NTri1, EVal[0]);
  //PrintMatrix(EFun);
  RotateWaveFunctions(PsiBasis,settings->NTri1,EFun,EigenBasis);
  BlockTriDiagonalMatrixType SE[NSelfEnergy];
  for (i = 0; i < NSelfEnergy; i++)
  {
	  CreateSpectrum(OppH1,*U[i],EigenBasis,settings->NTri1,settings->epsilon,SE[i],NTriSE);
  }
  //PrintBlockTriDiagonalMatrix(SE[0]);
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
    /*
     * 1) Evaluate SE(w,Gamma1)
     * 2) Add Eigenvalues to w
     * 3) Invert
     * 4) U^T.Gd
     * 5) (4).U_{;0}
     */

    CompactMatrixType SEw, temp;
    EvaluateBlockTriDiagonalMatrix(&SE[0],&SEw,wRes-EnergyContinuum[0],settings->Gamma1);
/*    if (p == 0)
    {
    	printf("wRes = %.16f\ngamma = %.16f\nSE[0].N = %d\n",wRes, settings->Gamma1,SE[0].N);
    	PrintCompactMatrixMinimal(SEw);
    }
*/
    for(i = 1; i < NSelfEnergy; i++)
    {
    	EvaluateBlockTriDiagonalMatrix(&SE[i],&temp,wRes-EnergyContinuum[i],settings->Gamma1);
    	CompactMatrixAdd(&SEw,&temp);
    	FreeCompactMatrix(&temp);
    }
    for(i = 0; i < settings->NTri1; i++)
    {
    	// Note that SE - (w - H) is easier, so later we nedd to add a minus sign
    	SEw.V[i+SEw.Ni*i] -= (wRes - EVal[i]); // TODO Imaginary Part includes broadening,too (i0) - shall I use it????????????????????????????????????????????????????????????
    	SEw.V[i+SEw.Ni*i].Im -= settings->Gamma1 /2.;

    	//EValInvIm[i] = -(settings->Gamma1/2.0) / ((wRes - EVal[i])*(wRes - EVal[i])+(settings->Gamma1/2.0)*(settings->Gamma1/2.0));
    }
    CompactMatrixInverse(&SEw);
    CompactMatrixType UTGdiag, Ueigen;
    Ueigen = MatrixTypeToCompactMatrixType(EFun);
    CompactMatrixMultiply('C','N',&Ueigen,&SEw,&UTGdiag);
/*    if (p == 0)
	{
		PrintCompactMatrixMinimal(UTGdiag);
	}-------------------------------------------------------------------------------------------------------------------- Debug below here
*/
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
          tempDoubleRe += UTGdiag.V[i*UTGdiag.Nj + k].Re * Ueigen.V[k*Ueigen.Nj].Re
        		  	    - UTGdiag.V[i*UTGdiag.Nj + k].Im * Ueigen.V[k*Ueigen.Nj].Im;
          tempDoubleIm += UTGdiag.V[i*UTGdiag.Nj + k].Re * Ueigen.V[k*Ueigen.Nj].Im
        		  	    + UTGdiag.V[i*UTGdiag.Nj + k].Im * Ueigen.V[k*Ueigen.Nj].Re;
        }
      }
      else
      {
        for(k = 0; k < settings->NTri1; k++)
        {
          tempDoubleRe += UTGdiag.V[i*UTGdiag.Nj + k].Re * Ueigen.VRe[k*Ueigen.Nj];
          tempDoubleIm += UTGdiag.V[i*UTGdiag.Nj + k].Im * Ueigen.VRe[k*Ueigen.Nj];
        }
      }
      // Here we correct for the missing minus sign (see above) by multiplying with -b0
      WaveFunctionAddComplexMultipliedUnionBasis(&GTpsi, PsiBasis[i],-SpecTri1.b[0]*tempDoubleRe,-SpecTri1.b[0]*tempDoubleIm);
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


