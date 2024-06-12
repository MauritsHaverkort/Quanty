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
  if(GammaL==0) // transform to spectrum with zero Gamma
  {
    switch (M->type|32)
    {
      case 'l':
      {
        A->Name=(char*)malloc(MAXStrLength*sizeof(char));
        strncpy(A->Name, ((BlockListOfPoles*)M->data)->Name, MAXStrLength);
        A->Spec = (SpectraType*)malloc( Sqr((((BlockListOfPoles*)(M->data))->BlockSize))*sizeof(SpectraType));
        BlockListOfPolesToSpectra(GroundStateEnergy, (BlockListOfPoles*)(M->data), A->Spec, Emin, Emax, NE);
      }
      break;
      case 't':
      {
        printf("TODO: Block Tridiagonal matrix to spectrum without broadening\n");
      }
      break;
      case 'a':
      {
        printf("TODO: Block Anderson matrix to spectrum without broadening\n");
      }
      break;
      case 'n':
      {
        printf("TODO: Block Natural Impurity matrix to spectrum without broadening\n");
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
        strncpy(A->Name, ((BlockListOfPoles*)M->data)->Name, MAXStrLength);
        A->Spec = (SpectraType*)malloc( Sqr((((BlockListOfPoles*)(M->data))->BlockSize))*sizeof(SpectraType));
        BlockListOfPolesToSpectraWithBroadening(GroundStateEnergy, (BlockListOfPoles*)(M->data), A->Spec, Emin, Emax, NE, GammaL);
      }
      break;
      case 't':
      {
        A->Name=(char*)malloc(MAXStrLength*sizeof(char));
        strncpy(A->Name, ((BlockTriDiagonalMatrixType*)M->data)->Name, MAXStrLength);
        A->Spec = (SpectraType*)malloc( Sqr((((BlockTriDiagonalMatrixType*)(M->data))->BlockSize))*sizeof(SpectraType));
        BlockTriDiagonalMatrixToSpectraWithBroadening(GroundStateEnergy,(BlockTriDiagonalMatrixType*)(M->data), A->Spec, Emin, Emax, NE, GammaL);
      }
      break;
      case 'a':
      {
        A->Name=(char*)malloc(MAXStrLength*sizeof(char));
        strncpy(A->Name, ((BlockAndersonMatrixType*)M->data)->Name, MAXStrLength);
        A->Spec = (SpectraType*)malloc( Sqr((((BlockAndersonMatrixType*)(M->data))->BlockSize))*sizeof(SpectraType));
        BlockAndersonMatrixToSpectraWithBroadening(GroundStateEnergy,(BlockAndersonMatrixType*)(M->data), A->Spec, Emin, Emax, NE, GammaL);
      }
      break;
      case 'n':
      {
        printf("TODO: Block Natural Impurity matrix to spectrum\n");
      }
      break;
    }
  }
  if(GammaG!=0)
  {
    for(i=0;i<=A->NSpectra;i++)
    {
      BroadenG(&(A->Spec[i]), GammaG);
    }
  }
  return 0;
}

uint32_t NoneBlockResponseFunctionToSpectra(const ResponsefunctionType *M, SpectraType *A, double GroundStateEnergy, double Emin, double Emax, double GammaG, double GammaL, uint32_t NE)
{
  printf("TODO: NoneBlockResponseFunctionToSpectra\n");
  if(GammaL==0) // transform to spectrum with zero Gamma
  {
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
  else // transform to spectra with non-zero Gamma
  {
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
 }
 return 0;
}

ResponsefunctionType *ConvertResponseFunction(ResponsefunctionType *M, const char type, double T, double SingularValue)
{
  if(M->Block==true) return ConvertBlockMatrix(M, type, T, SingularValue);
  //TODO: convert non-block response functions
  // see BasicMatch_MatrixTypeConverson.cpp
  printf("//TODO: convert non-block response functions \n");
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
  printf("ERROR: unrecognized type in response function\n");
  return NULL;
}
