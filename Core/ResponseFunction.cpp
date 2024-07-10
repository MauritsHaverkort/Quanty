/*  This function takes as imput a pointer to an empty SpectraTypeArray
 *  If GammaL is zero we calcualte the pole energy and add spectral weight to the grid points left and right of this
 *  Otherwise we broaden with a Lorenzian
 *  We broaden afterwards with a Gaussian
 */
uint32_t ResponseFunctionToSpectra(const ResponsefunctionType *M, SpectraTypeArray *A, double GroundStateEnergy, double Emin, double Emax, double GammaG, double GammaL, uint32_t NE)
{
  if(M->Block==true) return BlockResponseFunctionToSpectra(M, A, GroundStateEnergy, Emin, Emax, GammaG, GammaL, NE);
  A->NSpectra=1;
  A->Spec = (SpectraType*)malloc(sizeof(SpectraType));
  A->Name=(char*)malloc(MAXStrLength*sizeof(char));
  switch (M->type|32)
  {
    case 'l':
    {
      strncpy(A->Name, ((ListOfPoles*)M->data)->Name, MAXStrLength);
    }
    break;
    case 't':
    {
      strncpy(A->Name, ((TriDiagonalMatrixType*)M->data)->Name, MAXStrLength);
    }
    break;
    case 'a':
    {
      strncpy(A->Name, ((AndersonMatrixType*)M->data)->Name, MAXStrLength);
    }
    break;
    case 'n':
    {
      strncpy(A->Name, ((DoubleTriDiagonalMatrixType*)M->data)->Name, MAXStrLength);
    }
    break;
  }
  return NoneBlockResponseFunctionToSpectra(M, A->Spec, GroundStateEnergy, Emin, Emax, GammaG, GammaL, NE);
}


uint32_t BlockResponseFunctionToSpectra(const ResponsefunctionType *M, SpectraTypeArray *A, double GroundStateEnergy, double Emin, double Emax, double GammaG, double GammaL, uint32_t NE)
{
  int i;
  int ExtraPoints = 0;
  double EboarderShift =0;
  const double dE = (Emax - Emin) / NE;
  if(GammaG!=0)
  {
    ExtraPoints = ( 4 * GammaG ) / dE ;
    EboarderShift = ExtraPoints * dE;
    Emin -= EboarderShift;
    Emax += EboarderShift;
    NE   += 2*ExtraPoints;
  }
  if(GammaL==0) // transform to spectrum with zero Gamma
  {
    switch (M->type|32)
    {
      case 'l':
      {
        A->Name=(char*)malloc(MAXStrLength*sizeof(char));
        A->NSpectra = Sqr((((BlockListOfPoles*)(M->data))->BlockSize));
        A->Spec = (SpectraType*)malloc( A->NSpectra*sizeof(SpectraType));
        strncpy(A->Name, ((BlockListOfPoles*)M->data)->Name, MAXStrLength);
        BlockListOfPolesToSpectra(GroundStateEnergy, (BlockListOfPoles*)(M->data), A->Spec, Emin, Emax, NE);
      }
      break;
      case 't':
      {
        A->Name=(char*)malloc(MAXStrLength*sizeof(char));
        A->NSpectra = Sqr((((BlockTriDiagonalMatrixType*)(M->data))->BlockSize));
        A->Spec = (SpectraType*)malloc( A->NSpectra*sizeof(SpectraType));
        strncpy(A->Name, ((BlockTriDiagonalMatrixType*)M->data)->Name, MAXStrLength);
        BlockTriDiagonalMatrixToSpectra(GroundStateEnergy,(BlockTriDiagonalMatrixType*)(M->data), A->Spec, Emin, Emax, NE);
      }
      break;
      case 'a':
      {
        A->Name=(char*)malloc(MAXStrLength*sizeof(char));
        A->NSpectra = Sqr((((BlockAndersonMatrixType*)(M->data))->BlockSize));
        A->Spec = (SpectraType*)malloc( A->NSpectra*sizeof(SpectraType));
        strncpy(A->Name, ((BlockAndersonMatrixType*)M->data)->Name, MAXStrLength);
        BlockAndersonMatrixToSpectra(GroundStateEnergy,(BlockAndersonMatrixType*)(M->data), A->Spec, Emin, Emax, NE);
      }
      break;
      case 'n':
      {
        A->Name=(char*)malloc(MAXStrLength*sizeof(char));
        A->NSpectra = Sqr((((BlockNaturalImpurityOrbitalMatrixType*)(M->data))->BlockSize));
        A->Spec = (SpectraType*)malloc( A->NSpectra*sizeof(SpectraType));
        strncpy(A->Name, ((BlockNaturalImpurityOrbitalMatrixType*)M->data)->Name, MAXStrLength);
        BlockNaturalImpurityOrbitalMatrixToSpectra(GroundStateEnergy,(BlockNaturalImpurityOrbitalMatrixType*)(M->data), A->Spec, Emin, Emax, NE);
      }
      break;
    }
  }
  else // transform to spectra with non-zero Gamma
  {
    switch (M->type|32)
    {
      case 'l':
      {
        A->Name=(char*)malloc(MAXStrLength*sizeof(char));
        A->NSpectra = Sqr((((BlockListOfPoles*)(M->data))->BlockSize));
        A->Spec = (SpectraType*)malloc( A->NSpectra*sizeof(SpectraType));
        strncpy(A->Name, ((BlockListOfPoles*)M->data)->Name, MAXStrLength);
        BlockListOfPolesToSpectraWithBroadening(GroundStateEnergy, (BlockListOfPoles*)(M->data), A->Spec, Emin, Emax, NE, GammaL);
      }
      break;
      case 't':
      {
        A->Name=(char*)malloc(MAXStrLength*sizeof(char));
        A->NSpectra = Sqr((((BlockTriDiagonalMatrixType*)(M->data))->BlockSize));
        A->Spec = (SpectraType*)malloc( A->NSpectra*sizeof(SpectraType));
        strncpy(A->Name, ((BlockTriDiagonalMatrixType*)M->data)->Name, MAXStrLength);
        BlockTriDiagonalMatrixToSpectraWithBroadening(GroundStateEnergy,(BlockTriDiagonalMatrixType*)(M->data), A->Spec, Emin, Emax, NE, GammaL);
      }
      break;
      case 'a':
      {
        A->Name=(char*)malloc(MAXStrLength*sizeof(char));
        A->NSpectra = Sqr((((BlockAndersonMatrixType*)(M->data))->BlockSize));
        A->Spec = (SpectraType*)malloc( A->NSpectra*sizeof(SpectraType));
        strncpy(A->Name, ((BlockAndersonMatrixType*)M->data)->Name, MAXStrLength);
        BlockAndersonMatrixToSpectraWithBroadening(GroundStateEnergy,(BlockAndersonMatrixType*)(M->data), A->Spec, Emin, Emax, NE, GammaL);
      }
      break;
      case 'n':
      {
        A->Name=(char*)malloc(MAXStrLength*sizeof(char));
        A->NSpectra = Sqr((((BlockNaturalImpurityOrbitalMatrixType*)(M->data))->BlockSize));
        A->Spec = (SpectraType*)malloc( A->NSpectra*sizeof(SpectraType));
        strncpy(A->Name, ((BlockNaturalImpurityOrbitalMatrixType*)M->data)->Name, MAXStrLength);
        BlockNaturalImpurityOrbitalMatrixToSpectraWithBroadening(GroundStateEnergy,(BlockNaturalImpurityOrbitalMatrixType*)(M->data), A->Spec, Emin, Emax, NE, GammaL);
      }
      break;
    }
  }
  if(GammaG!=0)
  {
    for(i=0;i<A->NSpectra;i++)
    {
      BroadenG(&(A->Spec[i]), GammaG);
    }
    RemovePointsSpectraArray(A, ExtraPoints, ExtraPoints);
  }
  return 0;
}

uint32_t NoneBlockResponseFunctionToSpectra(const ResponsefunctionType *M, SpectraType *A, double GroundStateEnergy, double Emin, double Emax, double GammaG, double GammaL, uint32_t NE)
{
  int i;
  int ExtraPoints = 0;
  double EboarderShift =0;
  const double dE = (Emax - Emin) / NE;
  if(GammaG!=0)
  {
    ExtraPoints = ( 4 * GammaG ) / dE ;
    EboarderShift = ExtraPoints * dE;
    Emin -= EboarderShift;
    Emax += EboarderShift;
    NE   += 2*ExtraPoints;
  }
  A->phase=0;
  A->Emax=Emax;
  A->Emin=Emin;
  A->EmaxPole=Emax;
  A->EminPole=Emin;
  A->Gamma=GammaL;
  A->N=NE;
  A->dE=dE;
  A->ResonantEnergy=NAN;
  SpectrumAllocate(*A, A);
  if(GammaL==0) // transform to spectrum with zero Gamma
  {
    switch (M->type|32)
    {
      case 'l':
      {
        ListOfPolesToSpectrum(GroundStateEnergy, *(ListOfPoles*)M->data, A);
      }
        break;
      case 't':
      {
        TriDiagonalMatrixToSpectrum(GroundStateEnergy, *(TriDiagonalMatrixType*)M->data, A);
      }
        break;
      case 'a':
      {
        AndersonImpurityMatrixToSpectrum(GroundStateEnergy, *(AndersonMatrixType*)M->data, A);
      }
        break;
      case 'n':
      {
        BADoubleTriDiagonalMatrixToSpectrum(GroundStateEnergy, *(BADoubleTriDiagonalMatrixType*)M->data, A);
      }
        break;
    }
  }
  else // transform to spectra with non-zero Gamma
  {
    
    printf("TODO: NoneBlockResponseFunctionToSpectra with Lorentzian broadening\n");

    switch (M->type|32)
    {
     case 'l':
     {
       
     }
     break;
     case 't':
     {
     
     }
     break;
     case 'a':
     {
     
     }
     break;
     case 'n':
     {
     
     }
     break;
   }
 }
 if(GammaG!=0)
 {
   BroadenG(A, GammaG);
   RemovePointsSpectra(A, ExtraPoints, ExtraPoints);
 }
 return 0;
}
