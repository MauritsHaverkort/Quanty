// These functions have been merged into the main core
// changes to RobertTest have been made accordingly - easy to undo


uint32_t TriDiagonalMatrixBProduct(TriDiagonalMatrixType *SpecTri, uint32_t blow, uint32_t bhigh, double *bProduct)
{
  int i;

  if(blow > 0 && bhigh < SpecTri->N)
  {
    if(bhigh < blow)
      *bProduct = 1.0;//SpecTri->b[blow];
    else
    {
      *bProduct = 1.0;
      for(i = blow; i <= bhigh; i++)
        *bProduct *= SpecTri->b[i];
    }
    return 0;
  }
  else
    return 1;
}


uint32_t TriDiagonalMatrixToGreensFunction(TriDiagonalMatrixType * const SpecTri, double wRes, double energy, double gamma, MatrixType * const invMat)
{
  double thetaRe[SpecTri->N+1], thetaIm[SpecTri->N+1], phiRe[SpecTri->N + 1], phiIm[SpecTri->N + 1], bProduct, denom, shift;
  int i, j;

  /* flip sign of tridiagonal matrix (because we want inv(E+w-H+igamma))*/
//  TriDiagonalMatrixScalarMultiply(SpecTri, -1.0);
//  /* shift by E+w if necessary*/
//  if((energy + wRes) != 0)
//    ShiftTriDiagonalMatrix(SpecTri, energy + wRes);
  shift = energy + wRes;
  /* Initialize theta and phi at start of recurrence. Note that at

     http://en.wikipedia.org/wiki/Tridiagonal_matrix#Inversion

     theta indices are defined from 0...N (N+1 values), and that is how they
     are here as well. But the matrix indices in the algorithm start at 1,1 so
     the theta indices are shifted w.r.t. the matrix indices during the main
     loop later.  The used range of Phi is defined from 2 ... N+1 (N values),
     so its indices are shifted as well below.
  */
  thetaRe[0] = 1; thetaIm[0] = 0;
  thetaRe[1] = (shift - SpecTri->a[1]); thetaIm[1] = gamma;

  phiRe[SpecTri->N-1] = (shift - SpecTri->a[SpecTri->N]); phiIm[SpecTri->N-1] = gamma;
  phiRe[SpecTri->N] = 1.0; phiIm[SpecTri->N] = 0.0;
  /* Calculate the theta and phi recurrences */
  for(i = 2; i < SpecTri->N+1; i++)
  {
  /*Spectri->b should be multiplied by -1 as well, but it cancels out because always squared*/
    thetaRe[i] = (shift - SpecTri->a[i])*thetaRe[i-1] - gamma*thetaIm[i-1] - SpecTri->b[i-1]*SpecTri->b[i-1]*thetaRe[i-2];
    thetaIm[i] = (shift - SpecTri->a[i])*thetaIm[i-1] + gamma * thetaRe[i-1] - SpecTri->b[i-1]*SpecTri->b[i-1]*thetaIm[i-2];
  }

  for(i = SpecTri->N - 2; i >= 0; i--)
  {
    phiRe[i] = (shift - SpecTri->a[i+1])*phiRe[i+1] - gamma*phiIm[i+1] - SpecTri->b[i+1]*SpecTri->b[i+1]*phiRe[i+2];
    phiIm[i] = (shift - SpecTri->a[i+1])*phiIm[i+1] + gamma*phiRe[i+1] - SpecTri->b[i+1]*SpecTri->b[i+1]*phiIm[i+2];
  }

  /* Compute the matrix inverse */
  invMat->QComplex = 1;
  invMat->Ni = SpecTri->N;
  invMat->Nj = SpecTri->N;
  InitMatrix(invMat);
  for(i = 0; i < SpecTri->N; i++)
  {
    for(j = i; j < SpecTri->N; j++)
    {
      if(TriDiagonalMatrixBProduct(SpecTri, i+1, j - 1+1, &bProduct))
      {
        printf("TriDiagonalMatrixBProduct failed in TriDiagonalMatrixToGreensFunction\n");
        return 1;
      }
      invMat->VRe[i][j] = 0;
      invMat->VIm[i][j] = 0;
      /*Watch for divide by zeros. Shouldn't be a problem for nonsingular matrices.*/
      denom = thetaRe[SpecTri->N]*thetaRe[SpecTri->N] + thetaIm[SpecTri->N]*thetaIm[SpecTri->N];
      if(denom != 0)
      {
        /*Doing real and imaginary parts separately for:
        
                        bProduct*theta[i]*phi[j+1]/theta[N]
        */
        invMat->VRe[i][j] = bProduct*(thetaRe[i]*phiRe[j+1]*thetaRe[SpecTri->N] + 
                                      thetaRe[i]*phiIm[j+1]*thetaIm[SpecTri->N] -
                                      thetaIm[i]*phiIm[j+1]*thetaRe[SpecTri->N] +
                                      thetaIm[i]*phiRe[j+1]*thetaIm[SpecTri->N])/denom;
        invMat->VIm[i][j] = bProduct*(-thetaRe[i]*phiRe[j+1]*thetaIm[SpecTri->N] +
                                       thetaRe[i]*phiIm[j+1]*thetaRe[SpecTri->N] +
                                       thetaIm[i]*phiRe[j+1]*thetaRe[SpecTri->N] +
                                       thetaIm[i]*phiIm[j+1]*thetaIm[SpecTri->N])/denom; 
      }
      /* Handle the (-1)^(i+j) part. I'm assuming this way would
         be faster than pow(i+j,-1), but I'm not positive.  */
      if((i+j)%2 != 0)
      {
        invMat->VRe[i][j] *= -1;
        invMat->VIm[i][j] *= -1;
      }
      /* Symmetric part */
      invMat->VRe[j][i] = invMat->VRe[i][j];
      invMat->VIm[j][i] = invMat->VIm[i][j];
    }
  }
  return 0;
}


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

uint32_t StreamPrintRIXSSpectra(const uint32_t numwRes, SpectraType (*Spec)[numwRes], double *wRes, const uint32_t NSpec, FILE * stream)
{
  uint32_t i, j, p, NE;
  fprintf(stream, "Energy           ");
  for(i = 0; i < NSpec; i++)
  {
    /*Header for real part*/
    for(j = 0; j < numwRes; j++)
    {
      fprintf(stream, "%f          ", wRes[j]);
    }
    /*Header for imag part*/
    for(j = 0; j < numwRes; j++)
    {
      fprintf(stream, "%f          ", wRes[j]);
    }

  }
  fprintf(stream, "\n");

  NE = Spec[0][0].N;
  for (i = 1; i < NSpec; i++)
  {
    for(j = 0; j < numwRes; j++)
      NE = Max(NE, Spec[i][j].N);
  }
  for (j = 0; j <= NE; j++)
  {
    fprintf(stream, "%19.12E", (double)(Spec[0][0].dE * j + Spec[0][0].Emin));
    for (i = 0; i < NSpec; i++)
    {

      for(p = 0; p < numwRes; p++)
      {
        if(j <= Spec[i][p].N)
          fprintf(stream, " %22.15E ", Spec[i][p].Re[j]);
        else
          fprintf(stream, " %22.15E ", 0.0);
      }
      for(p = 0; p < numwRes; p++)
      {
        if(j <= Spec[i][p].N)
          fprintf(stream, " %22.15E ", Spec[i][p].Im[j]);
        else
          fprintf(stream, " %22.15E ", 0.0);
      }

    }
    fprintf(stream, "\n");
  }
  return 0;
}

uint32_t FilePrintRIXSSpectra(const uint32_t numwRes, SpectraType (*Spec)[numwRes], double *wRes, const uint32_t NSpec, const char *const flname, const char *const mode)
{
  FILE  *stream;
  uint32_t retval;
  if ((strcmp(mode, "w") && strcmp(mode, "a")) != 0)
  {
    printf
      ("ERROR in FilePrintRIXSSpectra.\nmode must be either \"w\" for overwrite or \"a\" for append\n %s as input is not excepted\n",
       mode);
    return 1;
  }
  if ((stream = fopen(flname, mode)) == NULL)
  {
    printf("ERROR\nCan't open %s for writing with mode %s.\n", flname, mode);
    return 1;
  }
  retval = StreamPrintRIXSSpectra(numwRes, Spec, wRes, NSpec, stream);
  fclose(stream);
  return retval;
}

uint32_t StreamPrintRIXSSpectrum(const uint32_t numwRes, SpectraType (*Spec)[numwRes], double *wRes, int printReal, const uint32_t iSpec, FILE * stream)
{
  //printReal == 1 : print the real part to the file
  //printReal == other : print the imaginary part
  
  uint32_t i, j, p, NE;
  fprintf(stream, "Energy           ");
  /*Header */
  for(j = 0; j < numwRes; j++)
  {
    fprintf(stream, "%f          ", wRes[j]);
  }

  fprintf(stream, "\n");

  NE = Spec[iSpec][0].N;

  for (j = 0; j <= NE; j++)
  {
    fprintf(stream, "%19.12E", (double)(Spec[iSpec][0].dE * j + Spec[iSpec][0].Emin));
    if(printReal == 1)
    {
      for(p = 0; p < numwRes; p++)
      {
        fprintf(stream, " %22.15E ", Spec[iSpec][p].Re[j]);
      }
    }
    else
    {
      for(p = 0; p < numwRes; p++)
      {
        fprintf(stream, " %22.15E ", -Spec[iSpec][p].Im[j]);
      }
    }
    fprintf(stream, "\n");
  }
  return 0;
}

uint32_t FilePrintRIXSSpectrum(const uint32_t numwRes, SpectraType (*Spec)[numwRes], double *wRes, int printReal, const uint32_t iSpec, const char *const flname, const char *const mode)
{
  //iSpec is index of the spectrum you want to print (e.g xx might 0, xy 1, xz 2, etc)
  FILE  *stream;
  uint32_t retval;
  
  if ((strcmp(mode, "w") && strcmp(mode, "a")) != 0)
  {
    printf
      ("ERROR in FilePrintRIXSSpectra.\nmode must be either \"w\" for overwrite or \"a\" for append\n %s as input is not excepted\n",
       mode);
    return 1;
  }
  if ((stream = fopen(flname, mode)) == NULL)
  {
    printf("ERROR\nCan't open %s for writing with mode %s.\n", flname, mode);
    return 1;
  }
  retval = StreamPrintRIXSSpectrum(numwRes, Spec, wRes, printReal, iSpec, stream);
  fclose(stream);
  return retval;
}