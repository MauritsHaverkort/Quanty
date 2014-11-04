//  Functions by Robert Green
//


uint32_t ResonantSpectra(const OperatorType OppH, const OperatorType Texc, const OperatorType *TDec, uint32_t NDec,  const WaveFunctionType psi, double energy, double epsilon, uint32_t numwRes, SpectraType (*spectra)[numwRes], double *wRes, double gammaRes, double gamma, double Gfwhm, uint32_t NKrylov1, uint32_t NKrylov2)
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
  
  OperatorType OppHRot;
  WaveFunctionType Tpsi,GTpsi,TGTpsi, TGTpsiRot, PsiBasis[NKrylov1];
  TriDiagonalMatrixType SpecTri1, SpecTri2;
  MatrixType tempMat,EFun,Density,DensityPart,TRot;
  double DensityValues[psi.NFermion];
  double EVal[NKrylov1], EValInvRe[NKrylov1], EValInvIm[NKrylov1],tempDoubleRe,tempDoubleIm;
  int i,p,j,k;
  int rotate = 0;
  
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
    WaveFunctionInit(&GTpsi);
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
        
        RotateOperatorNormalOrder(&OppHRot, OppH, TRot, &QDetNotOrderedFcFaBaBc);//&QDetNotOrderedFcFaBcBa);
        
        OperatorChop(&OppHRot, epsilon);
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
        LanczosTriDiagonalize(OppHRot, TGTpsiRot, epsilon, &SpecTri2, gamma);
      else
        LanczosTriDiagonalize(OppH, TGTpsi, epsilon, &SpecTri2, gamma);
      
      TriDiagonalMatrixToSpectrum(energy, SpecTri2, &spectra[i][p]);
      BroadenG(&spectra[i][p], Gfwhm);
      FreeTriDiagonalMatrix(&SpecTri2);
    }
    WaveFunctionFree(&GTpsi);
  }
  return 0;
}

uint32_t FluorescenceYield(const OperatorType OppH, const OperatorType Texc, const OperatorType *TDec, uint32_t NDec,  const WaveFunctionType psi, double energy, double epsilon, uint32_t numwRes, SpectraType spectra[NDec], double *wRes, double gammaRes, double Gfwhm, uint32_t NKrylov1)
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
  double EVal[NKrylov1], EValInvRe[NKrylov1], EValInvIm[NKrylov1],tempDoubleRe,tempDoubleIm;
  int i,p,k;
  
  /* Calculate T|i> */
  OperatorPsi(Texc,psi,&Tpsi);
  
  /* Use T|i> to transform H to tridiagonal form in Krylov basis */
  /* Using the version which gives back the Krylov basis vectors */
  SpecTri1.N = NKrylov1;
  InitTriDiagonalMatrix(&SpecTri1);
  LanczosTriDiagonalizeKrylov(OppH,Tpsi,epsilon,&SpecTri1,PsiBasis);
  
  /*Compute eigenvalues and eigenvectors, then S^-1 = QAQ^T where
   Q are the eigenvectors and A is a matrix with the eigenvalue
   inverses on the diagonal. */
  if(TriDiagonalMatrixToEigensystem(SpecTri1, &EFun, EVal))
  {
    printf("Problem with TridiagonalMatrixToEigenSystem\n");
  }
  
  /* Initialize Spectra */
  for(i = 0; i < NDec; i++)
  {
    sprintf(spectra[i].Name, "%d",i);
    spectra[i].N = numwRes;
    spectra[i].Emin = wRes[0];
    spectra[i].Emax = wRes[numwRes-1];
    spectra[i].dE = (wRes[1] - wRes[0]);
    spectra[i].Gamma = gammaRes;
    SpectrumAllocate(spectra[i], &spectra[i]);
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
    
    /* Calculate Q*A^(-1)*Q^(-1)[:,0] */
    GTpsi.NFermion = psi.NFermion;
    GTpsi.NBoson = psi.NBoson;
    GTpsi.QComplex = psi.QComplex;
    GTpsi.NBitsKey = psi.NBitsKey;
    WaveFunctionInit(&GTpsi);
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
      WaveFunctionCalculateNorm(TGTpsi, &spectra[i].Im[p]);
      WaveFunctionFree(&TGTpsi);
      /* Want norm squared, set negative for convention */
      spectra[i].Im[p] = - spectra[i].Im[p] * spectra[i].Im[p];
    }
    WaveFunctionFree(&GTpsi);
  }
  for(i = 0; i < NDec; i++)
  {
    /* Construct real part and broaden */
    KramersKronigImToRe(&spectra[i]);
    BroadenG(&spectra[i], Gfwhm);
  }
  return 0;
}



