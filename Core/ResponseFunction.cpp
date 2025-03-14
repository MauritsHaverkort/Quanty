/*  Using the Gershgoring Circle Theorem we find bounds for the eigenvalues of a matrix
 *
 *
 */
uint32_t ResponseFunctionEstimateEminAndEmaxPole(const ResponsefunctionType *const M, double *EminPole, double *EmaxPole)
{
  int i,j,k;
  double Rrow;
  *EminPole= DBL_MAX;
  *EmaxPole=-DBL_MAX;
  if(M->Block==true)
  { // Block matrices
    switch (M->type|32)
    {
      case 'l':
      {
        BlockListOfPoles *A=(BlockListOfPoles*)(M->data);
        if(A->N==0)
        {
          *EminPole=-1;
          *EmaxPole= 1;
          return 0;
        }
        if(A->QComplex==0)
        {
          for(i=0;i<A->N;i++)
          {
            if(A->N==1)
            {
              *EminPole=A->VRe[0]-1;
              *EmaxPole=A->VRe[0]+1;
              return 0;
            }
            if( A->VRe[i*(1+Sqr(A->N))] < *EminPole ) *EminPole=A->VRe[i*(1+Sqr(A->N))];
            if( A->VRe[i*(1+Sqr(A->N))] > *EmaxPole ) *EmaxPole=A->VRe[i*(1+Sqr(A->N))];
          }
        }
        else
        {
          if(A->N==1)
          {
            *EminPole=A->V[0].Re-1;
            *EmaxPole=A->V[0].Re+1;
            return 0;
          }
          for(i=0;i<A->N;i++)
          {
            if( A->V[i*(1+Sqr(A->N))].Re < *EminPole ) *EminPole=A->V[i*(1+Sqr(A->N))].Re;
            if( A->V[i*(1+Sqr(A->N))].Re > *EmaxPole ) *EmaxPole=A->V[i*(1+Sqr(A->N))].Re;
          }
        }
        return 0;
      }
      break;
      case 't':
      {
        BlockTriDiagonalMatrixType *A=(BlockTriDiagonalMatrixType*)(M->data);
        if(A->N==0)
        {
          *EminPole=-1;
          *EmaxPole= 1;
          return 0;
        }
        if(A->QComplex==0)
        { //real
          for(i=0;i<A->N;i++)
          {
            for(j=0;j<A->BlockSize;j++)
            {
              Rrow=0;
              if(i!=0)
              {
                for(k=0;k<A->BlockSize;k++)
                {
                  Rrow += fabs(A->bRe[i * Sqr(A->BlockSize) + j * A->BlockSize + k]);
                }
              }
              for(k=0;k<j;k++)
              {
                Rrow += fabs(A->aRe[(i+1) * Sqr(A->BlockSize) + j * A->BlockSize + k]);
              }
              for(k=j+1;k<A->BlockSize;k++)
              {
                Rrow += fabs(A->aRe[(i+1) * Sqr(A->BlockSize) + j * A->BlockSize + k]);
              }
              if(i!=A->N-1)
              {
                for(k=0;k<A->BlockSize;k++)
                {
                  Rrow += fabs(A->bRe[(i+1) * Sqr(A->BlockSize) + j * A->BlockSize + k]);
                }
              }
              if(A->aRe[(i+1)*Sqr(A->BlockSize)+j+j*A->BlockSize] - Rrow < *EminPole) *EminPole = A->aRe[(i+1)*Sqr(A->BlockSize)+j+j*A->BlockSize] - Rrow;
              if(A->aRe[(i+1)*Sqr(A->BlockSize)+j+j*A->BlockSize] + Rrow > *EmaxPole) *EmaxPole = A->aRe[(i+1)*Sqr(A->BlockSize)+j+j*A->BlockSize] + Rrow;
            }
          }
        }
        else
        { //complex
          for(i=0;i<A->N;i++)
          {
            for(j=0;j<A->BlockSize;j++)
            {
              Rrow=0;
              if(i!=0)
              {
                for(k=0;k<A->BlockSize;k++)
                {
                  Rrow += sqrt(Sqr(A->b[i * Sqr(A->BlockSize) + j * A->BlockSize + k].Re)+Sqr(A->b[i * Sqr(A->BlockSize) + j * A->BlockSize + k].Im));
                }
              }
              for(k=0;k<j;k++)
              {
                Rrow += sqrt(Sqr(A->a[(i+1) * Sqr(A->BlockSize) + j * A->BlockSize + k].Re) + Sqr(A->a[(i+1) * Sqr(A->BlockSize) + j * A->BlockSize + k].Im));
              }
              for(k=j+1;k<A->BlockSize;k++)
              {
                Rrow += sqrt(Sqr(A->a[(i+1) * Sqr(A->BlockSize) + j * A->BlockSize + k].Re) + Sqr(A->a[(i+1) * Sqr(A->BlockSize) + j * A->BlockSize + k].Im));
              }
              if(i!=A->N-1)
              {
                for(k=0;k<A->BlockSize;k++)
                {
                  Rrow += sqrt(Sqr(A->b[(i+1) * Sqr(A->BlockSize) + j * A->BlockSize + k].Re) + Sqr(A->b[(i+1) * Sqr(A->BlockSize) + j * A->BlockSize + k].Im));
                }
              }
              if(A->a[(i+1)*Sqr(A->BlockSize)+j+j*A->BlockSize].Re - Rrow < *EminPole) *EminPole = A->a[(i+1)*Sqr(A->BlockSize)+j+j*A->BlockSize].Re - Rrow;
              if(A->a[(i+1)*Sqr(A->BlockSize)+j+j*A->BlockSize].Re + Rrow > *EmaxPole) *EmaxPole = A->a[(i+1)*Sqr(A->BlockSize)+j+j*A->BlockSize].Re + Rrow;
            }
          }
        }
      }
      break;
      case 'a':
      {
        BlockAndersonMatrixType *A=(BlockAndersonMatrixType*)(M->data);
        if(A->N==0)
        {
          *EminPole=-1;
          *EmaxPole= 1;
          return 0;
        }
        if(A->QComplex==0)
        { //real
          for(j=0;j<A->BlockSize;j++) // loop over all rows in A0
          {
            Rrow=0;
            for(i=1;i<j;i++)
            {
              Rrow += fabs(A->aRe[(1)*Sqr(A->BlockSize)+i+j*A->BlockSize]);
            }
            for(i=j+1;i<A->BlockSize;i++)
            {
              Rrow += fabs(A->aRe[(1)*Sqr(A->BlockSize)+i+j*A->BlockSize]);
            }
            for(i=1;i<A->N;i++) // loop over all blocks in B
            {
              for(k=0;k<A->BlockSize;k++) // loop over all Columns in Bi
              {
                Rrow += fabs(A->bRe[(i)*Sqr(A->BlockSize)+k+j*A->BlockSize]);
              }
            }
            if(A->aRe[(1)*Sqr(A->BlockSize)+j+j*A->BlockSize] - Rrow < *EminPole) *EminPole = A->aRe[(1)*Sqr(A->BlockSize)+j+j*A->BlockSize] - Rrow;
            if(A->aRe[(1)*Sqr(A->BlockSize)+j+j*A->BlockSize] + Rrow > *EmaxPole) *EmaxPole = A->aRe[(1)*Sqr(A->BlockSize)+j+j*A->BlockSize] + Rrow;
          }
          for(i=1;i<A->N;i++) // loop over all blocks in B
          {
            for(j=0;j<A->BlockSize;j++) // loop over all rows in Ai+1
            {
              Rrow=0;
              for(k=0;k<j;k++)
              {
                Rrow += fabs(A->aRe[(i+1)*Sqr(A->BlockSize)+k+j*A->BlockSize]);
              }
              for(k=j+1;k<A->BlockSize;k++)
              {
                Rrow += fabs(A->aRe[(i+1)*Sqr(A->BlockSize)+k+j*A->BlockSize]);
              }
              for(k=0;k<A->BlockSize;k++)
              {
                Rrow += fabs(A->bRe[(i)*Sqr(A->BlockSize)+k+j*A->BlockSize]);
              }
              if(A->aRe[(i+1)*Sqr(A->BlockSize)+j+j*A->BlockSize] - Rrow < *EminPole) *EminPole = A->aRe[(i+1)*Sqr(A->BlockSize)+j+j*A->BlockSize] - Rrow;
              if(A->aRe[(i+1)*Sqr(A->BlockSize)+j+j*A->BlockSize] + Rrow > *EmaxPole) *EmaxPole = A->aRe[(i+1)*Sqr(A->BlockSize)+j+j*A->BlockSize] + Rrow;
            }
          }
        }
        else
        { // complex
          for(j=0;j<A->BlockSize;j++) // loop over all rows in A0
          {
            Rrow=0;
            for(i=1;i<j;i++)
            {
              Rrow += sqrt(Sqr(A->a[(1)*Sqr(A->BlockSize)+i+j*A->BlockSize].Re) + Sqr(A->a[(1)*Sqr(A->BlockSize)+i+j*A->BlockSize].Im));
            }
            for(i=j+1;i<A->BlockSize;i++)
            {
              Rrow += sqrt(Sqr(A->a[(1)*Sqr(A->BlockSize)+i+j*A->BlockSize].Re) + Sqr(A->a[(1)*Sqr(A->BlockSize)+i+j*A->BlockSize].Im));
            }
            for(i=1;i<A->N;i++) // loop over all blocks in B
            {
              for(k=0;k<A->BlockSize;k++) // loop over all Columns in Bi
              {
                Rrow += sqrt(Sqr(A->b[(i)*Sqr(A->BlockSize)+k+j*A->BlockSize].Re) + Sqr(A->b[(i)*Sqr(A->BlockSize)+k+j*A->BlockSize].Im));
              }
            }
            if(A->a[(1)*Sqr(A->BlockSize)+j+j*A->BlockSize].Re - Rrow < *EminPole) *EminPole = A->a[(1)*Sqr(A->BlockSize)+j+j*A->BlockSize].Re - Rrow;
            if(A->a[(1)*Sqr(A->BlockSize)+j+j*A->BlockSize].Re + Rrow > *EmaxPole) *EmaxPole = A->a[(1)*Sqr(A->BlockSize)+j+j*A->BlockSize].Re + Rrow;
          }
          for(i=1;i<A->N;i++) // loop over all blocks in B
          {
            for(j=0;j<A->BlockSize;j++) // loop over all rows in Ai+1
            {
              Rrow=0;
              for(k=0;k<j;k++)
              {
                Rrow += sqrt(Sqr(A->a[(i+1)*Sqr(A->BlockSize)+k+j*A->BlockSize].Re) + Sqr(A->a[(i+1)*Sqr(A->BlockSize)+k+j*A->BlockSize].Im));
              }
              for(k=j+1;k<A->BlockSize;k++)
              {
                Rrow += sqrt(Sqr(A->a[(i+1)*Sqr(A->BlockSize)+k+j*A->BlockSize].Re) + Sqr(A->a[(i+1)*Sqr(A->BlockSize)+k+j*A->BlockSize].Im));
              }
              for(k=0;k<A->BlockSize;k++)
              {
                Rrow += sqrt(Sqr(A->b[(i)*Sqr(A->BlockSize)+k+j*A->BlockSize].Re) + Sqr(A->b[(i)*Sqr(A->BlockSize)+k+j*A->BlockSize].Im));
              }
              if(A->a[(i+1)*Sqr(A->BlockSize)+j+j*A->BlockSize].Re - Rrow < *EminPole) *EminPole = A->a[(i+1)*Sqr(A->BlockSize)+j+j*A->BlockSize].Re - Rrow;
              if(A->a[(i+1)*Sqr(A->BlockSize)+j+j*A->BlockSize].Re + Rrow > *EmaxPole) *EmaxPole = A->a[(i+1)*Sqr(A->BlockSize)+j+j*A->BlockSize].Re + Rrow;
            }
          }
        }
      }
      break;
      case 'n':
      {
        double temp;
        ResponsefunctionType R;
        BlockNaturalImpurityOrbitalMatrixType *A=(BlockNaturalImpurityOrbitalMatrixType*)(M->data);
        R.type='t';
        R.Block=true;
        R.SuggestedGrid=NULL;
        R.NSuggestedGridPoints=0;
        R.data = (void*)(&(A->con));
        ResponseFunctionEstimateEminAndEmaxPole(&R, &temp, EmaxPole);
        R.data = (void*)(&(A->val));
        ResponseFunctionEstimateEminAndEmaxPole(&R, EminPole, &temp);
        return 0;
      }
      break;
    }
  }
  else
  { // single valued entries
    switch (M->type|32)
    {
      case 'l':
      {
        ListOfPoles *A=(ListOfPoles*)(M->data);
        if(A->N==0)
        {
          *EminPole=-1;
          *EmaxPole= 1;
          return 0;
        }
        for(i=0;i<A->N;i++)
        {
          if(A->a[i]<*EminPole) *EminPole=A->a[i];
          if(A->a[i]>*EmaxPole) *EmaxPole=A->a[i];
        }
        return 0;
      }
      break;
      case 't':
      {
        TriDiagonalMatrixType *A=(TriDiagonalMatrixType*)(M->data);
        if(A->N==0)
        {
          *EminPole=-1;
          *EmaxPole= 1;
          return 0;
        }
        if(A->N==1)
        {
          *EminPole=A->a[1]-1;
          *EmaxPole=A->a[1]+1;
          return 0;
        }
        *EminPole = A->a[1] - fabs(A->b[1]);
        *EmaxPole = A->a[1] + fabs(A->b[1]);
        for(i=2;i<A->N;i++)
        {
          if(*EminPole > A->a[i] - fabs(A->b[i-1]) - fabs(A->b[i])) *EminPole = A->a[i] - fabs(A->b[i-1]) - fabs(A->b[i]);
          if(*EmaxPole < A->a[i] + fabs(A->b[i-1]) + fabs(A->b[i])) *EmaxPole = A->a[i] + fabs(A->b[i-1]) + fabs(A->b[i]);
        }
        if(*EminPole > A->a[A->N] - fabs(A->b[A->N-1])) *EminPole = A->a[A->N] - fabs(A->b[A->N-1]);
        if(*EmaxPole < A->a[A->N] + fabs(A->b[A->N-1])) *EmaxPole = A->a[A->N] + fabs(A->b[A->N-1]);
        return 0;
      }
      break;
      case 'a':
      {
        AndersonMatrixType *A=(AndersonMatrixType*)(M->data);
        if(A->N==0)
        {
          *EminPole=-1;
          *EmaxPole= 1;
          return 0;
        }
        if(A->N==1)
        {
          *EminPole=A->a[1]-1;
          *EmaxPole=A->a[1]+1;
          return 0;
        }
        *EminPole=A->a[1];
        *EmaxPole=A->a[1];
        for(i=1;i<A->N;i++)
        {
          *EminPole -= fabs(A->b[i]);
          *EmaxPole += fabs(A->b[i]);
        }
        for(i=2;i<=A->N;i++)
        {
          if(*EminPole > A->a[i] - 2*fabs(A->b[i-1]) ) *EminPole = A->a[i] - 2*fabs(A->b[i-1]) ;
          if(*EmaxPole < A->a[i] + 2*fabs(A->b[i-1]) ) *EmaxPole = A->a[i] + 2*fabs(A->b[i-1]) ;
        }
        return 0;
      }
      break;
      case 'n':
      {
        double temp;
        ResponsefunctionType R;
        BADoubleTriDiagonalMatrixType *A=(BADoubleTriDiagonalMatrixType*)(M->data);
        R.type='t';
        R.Block=false;
        R.SuggestedGrid=NULL;
        R.NSuggestedGridPoints=0;
        R.data = (void*)(&(A->con));
        ResponseFunctionEstimateEminAndEmaxPole(&R, &temp, EmaxPole);
        R.data = (void*)(&(A->val));
        ResponseFunctionEstimateEminAndEmaxPole(&R, EminPole, &temp);
        return 0;
      }
      break;
    }
  }
  
  return 0;
}

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
 if(GammaG!=0)
 {
   BroadenG(A, GammaG);
   RemovePointsSpectra(A, ExtraPoints, ExtraPoints);
 }
 return 0;
}

uint32_t ResponseFunctionGetBlockSize(const ResponsefunctionType *const a)
{
  if(a->Block==false) return 1;
  switch(a->type)
  {
    case 'L':
    case 'l':
    {
      return ((BlockListOfPoles*)a->data)->BlockSize;
      break;
    }
    case 'A':
    case 'a':
    {
      return ((BlockAndersonMatrixType*)a->data)->BlockSize;
      break;
    }
    case 'T':
    case 't':
    {
      return ((BlockTriDiagonalMatrixType*)a->data)->BlockSize;
      break;
    }
    case 'N':
    case 'n':
    {
      return ((BlockNaturalImpurityOrbitalMatrixType*)a->data)->BlockSize;
      break;
    }
    default:
    {
      printf("Error in ResponseFunctionGetBlockSize: unrecognized type %c",a->type);
      return 0;
    }
  }
}

uint32_t ResponseFunctionGetQComplex(const ResponsefunctionType *const a)
{
  if(a->Block==false) return 0;
  switch(a->type)
  {
    case 'L':
    case 'l':
    {
      return ((BlockListOfPoles*)a->data)->QComplex;
      break;
    }
    case 'A':
    case 'a':
    {
      return ((BlockAndersonMatrixType*)a->data)->QComplex;
      break;
    }
    case 'T':
    case 't':
    {
      return ((BlockTriDiagonalMatrixType*)a->data)->QComplex;
      break;
    }
    case 'N':
    case 'n':
    {
      return ((BlockNaturalImpurityOrbitalMatrixType*)a->data)->QComplex;
      break;
    }
    default:
    {
      printf("Error in ResponseFunctionGetBlockSize: unrecognized type %c",a->type);
      return 0;
    }
  }
}

uint32_t ResponseFunctionGetA0(const ResponsefunctionType *const a, CompactMatrixType *M)
{
  M->Ni = ResponseFunctionGetBlockSize(a);
  M->Nj = M->Ni;
  M->QComplex = ResponseFunctionGetQComplex(a);
  InitCompactMatrix(M);
  if(a->Block==false)
  {
    switch(a->type)
    {
      case 'L':
      case 'l':
      {
        M->VRe[0] = ((ListOfPoles*)a->data)->a0;
        break;
      }
      case 'A':
      case 'a':
      {
        M->VRe[0] = ((AndersonMatrixType*)a->data)->a[0];
        break;
      }
      case 'T':
      case 't':
      {
        M->VRe[0] = ((TriDiagonalMatrixType*)a->data)->a[0];
        break;
      }
      case 'N':
      case 'n':
      {
        M->VRe[0] = ((BADoubleTriDiagonalMatrixType*)a->data)->a0;
        break;
      }
      default:
      {
        printf("Error in ResponseFunctionGetBlockSize: unrecognized type %c",a->type);
        return 1;
      }
    }
  }
  else
  {
    switch(a->type)
    {
      case 'L':
      case 'l':
      {
        if(M->QComplex==0)
        {
          for(int i=0; i<Sqr(M->Ni);i++)
          {
            M->VRe[i] = ((BlockListOfPoles*)a->data)->a0Re[i];
          }
        }
        else
        {
          for(int i=0; i<Sqr(M->Ni);i++)
          {
            M->V[i].Re = ((BlockListOfPoles*)a->data)->a0[i].Re;
            M->V[i].Im = ((BlockListOfPoles*)a->data)->a0[i].Im;
          }
        }
        break;
      }
      case 'A':
      case 'a':
      {
        if(M->QComplex==0)
        {
          for(int i=0; i<Sqr(M->Ni);i++)
          {
            M->VRe[i] = ((BlockAndersonMatrixType*)a->data)->aRe[i];
          }
        }
        else
        {
          for(int i=0; i<Sqr(M->Ni);i++)
          {
            M->V[i].Re = ((BlockAndersonMatrixType*)a->data)->a[i].Re;
            M->V[i].Im = ((BlockAndersonMatrixType*)a->data)->a[i].Im;
          }
        }
        break;
      }
      case 'T':
      case 't':
      {
        if(M->QComplex==0)
        {
          for(int i=0; i<Sqr(M->Ni);i++)
          {
            M->VRe[i] = ((BlockTriDiagonalMatrixType*)a->data)->aRe[i];
          }
        }
        else
        {
          for(int i=0; i<Sqr(M->Ni);i++)
          {
            M->V[i].Re = ((BlockTriDiagonalMatrixType*)a->data)->a[i].Re;
            M->V[i].Im = ((BlockTriDiagonalMatrixType*)a->data)->a[i].Im;
          }
        }
        break;
      }
      case 'N':
      case 'n':
      {
        if(M->QComplex==0)
        {
          for(int i=0; i<Sqr(M->Ni);i++)
          {
            M->VRe[i] = ((BlockNaturalImpurityOrbitalMatrixType*)a->data)->a0.VRe[i];
          }
        }
        else
        {
          for(int i=0; i<Sqr(M->Ni);i++)
          {
            M->V[i].Re = ((BlockNaturalImpurityOrbitalMatrixType*)a->data)->a0.V[i].Re;
            M->V[i].Im = ((BlockNaturalImpurityOrbitalMatrixType*)a->data)->a0.V[i].Im;
          }
        }
        break;
      }
      default:
      {
        printf("Error in ResponseFunctionGetBlockSize: unrecognized type %c",a->type);
        return 0;
      }
    }
  }
  return 0;
}

uint32_t EvaluateResponsefunction(ResponsefunctionType const *A, CompactMatrixType * Gw, const double omega, const double Gamma)
{
  if(A->Block==true) return EvaluateBlockResponsefunctionType(A, Gw, omega, Gamma);
  Gw->Ni=1;
  Gw->Nj=1;
  Gw->QComplex=1;
  InitCompactMatrix(Gw);
  Gw->V[0] = EvaluateResponsefunctionType(A, omega, Gamma);
  return 0;
}

bool ResponseFunctionEq(const ResponsefunctionType *const a, const ResponsefunctionType *const b, double zero)
{
  int i,j;
  double omega0, omega1, omega0shift;
  CompactMatrixType Ma,Mb,A0;
  // test dimension
  if(ResponseFunctionGetBlockSize(a)!=ResponseFunctionGetBlockSize(b)) return false;
  // compare a0
  ResponseFunctionGetA0(a, &Ma);
  ResponseFunctionGetA0(b, &A0);
  CompactMatrixSub(&Ma, &A0);
  if(Ma.QComplex==0)
  {
    for(i=0;i<Ma.Ni*Ma.Nj;i++)
    {
      if(Abs(Ma.VRe[i])>zero)
      {
        FreeCompactMatrix(&Ma);
        FreeCompactMatrix(&A0);
        return false;
      }
    }
  }
  else
  {
    for(i=0;i<Ma.Ni*Ma.Nj;i++)
    {
      if((Abs(Ma.V[i].Re)>zero)||(Abs(Ma.V[i].Im)>zero))
      {
        FreeCompactMatrix(&Ma);
        FreeCompactMatrix(&A0);
        return false;
      }
    }
  }
  // determine gridpoints and gamma for evaluation and compare values after evaluation
  omega0=0;
  omega1 = 1E10;
  zero=zero*(Ma.Ni*Ma.Nj);
  FreeCompactMatrix(&Ma);
  if(A0.QComplex==0) MakeCompactRealMatrixComplex(&A0);;
  for(j=1;j<13;j++)
  {
    omega0shift=0;
    EvaluateResponsefunction(a, &Ma, omega0-omega1, omega1/10.0);
    EvaluateResponsefunction(b, &Mb, omega0-omega1, omega1/10.0);
    for(i=0;i<Ma.Ni*Ma.Nj;i+=Ma.Ni+1)
    {
      omega0shift+= (1.0/(Ma.V[i].Re-A0.V[i].Re));
    }
    CompactMatrixSub(&Ma, &Mb);
    FreeCompactMatrix(&Mb);
    for(i=0;i<Ma.Ni*Ma.Nj;i++)
    {
      if((Abs(Ma.V[i].Re)>zero)||(Abs(Ma.V[i].Im)>zero))
      {
        FreeCompactMatrix(&Ma);
        FreeCompactMatrix(&A0);
        return false;
      }
    }
    FreeCompactMatrix(&Ma);
    EvaluateResponsefunction(a, &Ma, omega0+omega1, omega1/10.0);
    EvaluateResponsefunction(b, &Mb, omega0+omega1, omega1/10.0);
    for(i=0;i<Ma.Ni*Ma.Nj;i+=Ma.Ni+1)
    {
      omega0shift+= 1.0/(Ma.V[i].Re-A0.V[i].Re);
    }
    CompactMatrixSub(&Ma, &Mb);
    FreeCompactMatrix(&Mb);
    for(i=0;i<Ma.Ni*Ma.Nj;i++)
    {
      if((Abs(Ma.V[i].Re)>zero)||(Abs(Ma.V[i].Im)>zero))
      {
        FreeCompactMatrix(&Ma);
        FreeCompactMatrix(&A0);
        return false;
      }
    }
    omega0shift=omega0shift*(-0.5/(double)Ma.Ni);
    omega0 += omega0shift;
    omega1 *=0.1;
    FreeCompactMatrix(&Ma);
  }
  FreeCompactMatrix(&A0);
  return true;
}


/*   Transform G0 and Sigma to Tri-diagonal matrices (if they are not already in that form
 *   Create one big matrix of dimension G->N * (Sigma->N+1)
 *
 *   Epsilon is needed when creating a tri-diagonal matrix to determine if one needs to deflate
 */
int ResponseFunctionAddSelfEnergyDenseViaTri(const ResponsefunctionType *const G0, const ResponsefunctionType *const Sigma, ResponsefunctionType * G, double Epsilon)
{
  if(G0->Block!=Sigma->Block)
  {
    printf("Error in ResponseFunctionAddSelfEnergyDenseViaTri\n G0 and Sigma must either be both of Block type, or both of single value type\n");
    return 1;
  }
  if(G0->Block==true)
  {
    G->Block=true;
    G->data=malloc(sizeof(BlockListOfPoles));
    G->type='l';
    G->NSuggestedGridPoints=0;
    G->SuggestedGrid=NULL;
    BlockTriDiagonalMatrixType *G0Data, *SigmaData;
    if((G0->type=='t')||(G0->type=='T'))
    {
      G0Data = (BlockTriDiagonalMatrixType*)G0->data;
    }
    else
    {
      G0Data = (BlockTriDiagonalMatrixType*)malloc(sizeof(BlockTriDiagonalMatrixType));
      switch (G0->type|32)
      {
        case 'l':
        {
          BlockListOfPolesToBlockTridiagonalMatrix(*(BlockListOfPoles*)G0->data, G0Data, 1+((BlockListOfPoles*)G0->data)->N, Epsilon);
        }
        break;
        case 'a':
        {
          BlockAndersonMatrixToBlockTridiagonalMatrix((BlockAndersonMatrixType*)G0->data, G0Data, Epsilon);
        }
        break;
        case 'n':
        {
          BlockNaturalImpurityOrbitalMatrixToBlockTridiagonalMatrix((BlockNaturalImpurityOrbitalMatrixType*)G0->data, G0Data, Epsilon);
        }
        break;
        default:
        {
          printf("Error in ResponseFunctionAddSelfEnergyDenseViaTri unrecognized type to convert to: %c",G0->type);
        }
        break;
      }
    }
    if((Sigma->type=='t')||(Sigma->type=='T'))
    {
      SigmaData = (BlockTriDiagonalMatrixType*)Sigma->data;
    }
    else
    {
      SigmaData = (BlockTriDiagonalMatrixType*)malloc(sizeof(BlockTriDiagonalMatrixType));
      switch (Sigma->type|32)
      {
        case 'l':
        {
          BlockListOfPolesToBlockTridiagonalMatrix(*(BlockListOfPoles*)Sigma->data, SigmaData, 1+((BlockListOfPoles*)G0->data)->N, Epsilon);
        }
        break;
        case 'a':
        {
          BlockAndersonMatrixToBlockTridiagonalMatrix((BlockAndersonMatrixType*)Sigma->data, SigmaData, Epsilon);
        }
        break;
        case 'n':
        {
          BlockNaturalImpurityOrbitalMatrixToBlockTridiagonalMatrix((BlockNaturalImpurityOrbitalMatrixType*)Sigma->data, SigmaData, Epsilon);
        }
        break;
        default:
        {
          printf("Error in ResponseFunctionAddSelfEnergyDenseViaTri unrecognized type to convert to: %c",Sigma->type);
        }
        break;
      }
    }
    if(ResponseFunctionAddSelfEnergyDenseBlockTriTriLP(G0Data, SigmaData, (BlockListOfPoles*)G->data))
    {
      printf("Error in ResponseFunctionAddSelfEnergyDenseViaTri\n function ResponseFunctionAddSelfEnergyDenseBlockTriTriLP failed\n");
    }
    if(!((G0->type=='t')||(G0->type=='T')))
    {
      FreeBlockTridiagonalMatrix(G0Data);
      free(G0Data);
    }
    if(!((Sigma->type=='t')||(Sigma->type=='T')))
    {
      FreeBlockTridiagonalMatrix(SigmaData);
      free(SigmaData);
    }
  }
  else // non block version
  {
    G->Block=false;
    G->data=malloc(sizeof(ListOfPoles));
    G->type='l';
    G->NSuggestedGridPoints=0;
    G->SuggestedGrid=NULL;
    TriDiagonalMatrixType *G0Data, *SigmaData;
    if((G0->type=='t')||(G0->type=='T'))
    {
      G0Data = (TriDiagonalMatrixType*)G0->data;
    }
    else
    {
      G0Data = (TriDiagonalMatrixType*)malloc(sizeof(TriDiagonalMatrixType));
      switch (G0->type|32)
      {
        case 'l':
        {
          ListOfPolesToTridiagonalMatrix(*(ListOfPoles*)G0->data, G0Data, 1+((ListOfPoles*)G0->data)->N);
        }
        break;
        case 'a':
        {
          AndersonMatrixToTriDiagonalMatrix(*(AndersonMatrixType*)G0->data, G0Data, ((AndersonMatrixType*)G0->data)->N );
        }
        break;
        case 'n':
        {
          BADoubleTriDiagonalMatrixToTriDiagonalMatrix((BADoubleTriDiagonalMatrixType*)G0->data, G0Data);
        }
        break;
        default:
        {
          printf("Error in ResponseFunctionAddSelfEnergyDenseViaTri unrecognized type to convert to: %c",G0->type);
        }
        break;
      }
    }
    if((Sigma->type=='t')||(Sigma->type=='T'))
    {
      SigmaData = (TriDiagonalMatrixType*)Sigma->data;
    }
    else
    {
      SigmaData = (TriDiagonalMatrixType*)malloc(sizeof(TriDiagonalMatrixType));
      switch (Sigma->type|32)
      {
        case 'l':
        {
          ListOfPolesToTridiagonalMatrix(*(ListOfPoles*)Sigma->data, SigmaData, 1+((ListOfPoles*)G0->data)->N);
        }
        break;
        case 'a':
        {
          AndersonMatrixToTriDiagonalMatrix(*(AndersonMatrixType*)Sigma->data, SigmaData, ((AndersonMatrixType*)Sigma->data)->N);
        }
        break;
        case 'n':
        {
          BADoubleTriDiagonalMatrixToTriDiagonalMatrix((BADoubleTriDiagonalMatrixType*)Sigma->data, SigmaData);
        }
        break;
        default:
        {
          printf("Error in ResponseFunctionAddSelfEnergyDenseViaTri unrecognized type to convert to: %c",Sigma->type);
        }
        break;
      }
    }
    if(ResponseFunctionAddSelfEnergyDenseTriTriLP(G0Data, SigmaData, (ListOfPoles*)G->data))
    {
      printf("Error in ResponseFunctionAddSelfEnergyDenseViaTri\n function ResponseFunctionAddSelfEnergyDenseBlockTriTriLP failed\n");
    }
    if(!((G0->type=='t')||(G0->type=='T')))
    {
      FreeTriDiagonalMatrix(G0Data);
      free(G0Data);
    }
    if(!((Sigma->type=='t')||(Sigma->type=='T')))
    {
      FreeTriDiagonalMatrix(SigmaData);
      free(SigmaData);
    }
  }
  return 0;
}
int ResponseFunctionAddSelfEnergyDenseTriTriLP(const TriDiagonalMatrixType *const G0, const TriDiagonalMatrixType *const Sigma, ListOfPoles * G)
{
  BlockTriDiagonalMatrixType G0Block, SigmaBlock;
  BlockListOfPoles GBlock;
  G0Block.N=G0->N;
  G0Block.BlockSize=1;
  G0Block.QComplex=0;
  G0Block.aRe = G0->a;
  G0Block.bRe = G0->b;
  SigmaBlock.N=Sigma->N;
  SigmaBlock.BlockSize=1;
  SigmaBlock.QComplex=0;
  SigmaBlock.aRe = Sigma->a;
  SigmaBlock.bRe = Sigma->b;
  ResponseFunctionAddSelfEnergyDenseBlockTriTriLP(&G0Block,&SigmaBlock,&GBlock);
  G->N = GBlock.N;
  InitListOfPoles(G);
  G->mu = GBlock.mu;
  G->a0 = GBlock.a0Re[0];
  for(int i=0;i<G->N;i++)
  {
    G->a[i] = GBlock.VRe[2*i];
    G->b[i] = GBlock.VRe[2*i+1];
  }
  FreeBlockListOfPoles(&GBlock);
  return 0;
}
int ResponseFunctionAddSelfEnergyDenseBlockTriTriLP(const BlockTriDiagonalMatrixType *const G0, const BlockTriDiagonalMatrixType *const Sigma, BlockListOfPoles * G)
{
  CompactMatrixType GMat;
  GMat.Ni=G0->BlockSize * G0->N * (Sigma->N+1);
  GMat.Nj=G0->BlockSize * G0->N * (Sigma->N+1);
  GMat.QComplex = G0->QComplex | Sigma->QComplex;
  if(InitCompactMatrix(&GMat))
  {
    printf("Error: InitCompactMatrix failed in ResponseFunctionAddSelfEnergyDenseBlockTriTriLP\n");
    return 1;
  }
  // add Sigma to GMat
  if(GMat.QComplex==0)
  {  // GMat and Sigma real
    for(int k=0;k<Sigma->N;k++)
    {
      for(int l=0;l<G0->N;l++)
      {
        for(int i=0;i<G0->BlockSize;i++)
        {
          for(int j=0;j<G0->BlockSize;j++)
          {
            GMat.VRe[ i * GMat.Nj + j + (GMat.Nj+1) * ((k + l * (Sigma->N+1)) * Sigma->BlockSize )                              ] = Sigma->aRe[j+i*Sigma->BlockSize + k*Sqr(Sigma->BlockSize)];
            GMat.VRe[ i * GMat.Nj + j + (GMat.Nj+1) * ((k + l * (Sigma->N+1)) * Sigma->BlockSize ) + Sigma->BlockSize           ] = Sigma->bRe[j+i*Sigma->BlockSize + k*Sqr(Sigma->BlockSize)];
            GMat.VRe[ i * GMat.Nj + j + (GMat.Nj+1) * ((k + l * (Sigma->N+1)) * Sigma->BlockSize ) + Sigma->BlockSize * GMat.Nj ] = Sigma->bRe[i+j*Sigma->BlockSize + k*Sqr(Sigma->BlockSize)];
          }
        }
      }
    }
    {
      int k=Sigma->N;
      for(int l=0;l<G0->N;l++)
      {
        for(int i=0;i<G0->BlockSize;i++)
        {
          for(int j=0;j<G0->BlockSize;j++)
          {
            GMat.VRe[ i * GMat.Nj + j + (GMat.Nj+1) * ((k + l * (Sigma->N+1)) * Sigma->BlockSize ) ] = Sigma->aRe[j+i*Sigma->BlockSize + k*Sqr(Sigma->BlockSize)];
          }
        }
      }
    }
  }
  else if(Sigma->QComplex==0)
  { // Sigma Real but GMat complex
    for(int k=0;k<Sigma->N;k++)
    {
      for(int l=0;l<G0->N;l++)
      {
        for(int i=0;i<G0->BlockSize;i++)
        {
          for(int j=0;j<G0->BlockSize;j++)
          {
            GMat.V[ i * GMat.Nj + j + (GMat.Nj+1) * ((k + l * (Sigma->N+1)) * Sigma->BlockSize )                              ].Re = Sigma->aRe[j+i*Sigma->BlockSize + k*Sqr(Sigma->BlockSize)];
            GMat.V[ i * GMat.Nj + j + (GMat.Nj+1) * ((k + l * (Sigma->N+1)) * Sigma->BlockSize ) + Sigma->BlockSize           ].Re = Sigma->bRe[j+i*Sigma->BlockSize + k*Sqr(Sigma->BlockSize)];
            GMat.V[ i * GMat.Nj + j + (GMat.Nj+1) * ((k + l * (Sigma->N+1)) * Sigma->BlockSize ) + Sigma->BlockSize * GMat.Nj ].Re = Sigma->bRe[i+j*Sigma->BlockSize + k*Sqr(Sigma->BlockSize)];
          }
        }
      }
    }
    {
      int k=Sigma->N;
      for(int l=0;l<G0->N;l++)
      {
        for(int i=0;i<G0->BlockSize;i++)
        {
          for(int j=0;j<G0->BlockSize;j++)
          {
            GMat.V[ i * GMat.Nj + j + (GMat.Nj+1) * ((k + l * (Sigma->N+1)) * Sigma->BlockSize ) ].Re = Sigma->aRe[j+i*Sigma->BlockSize + k*Sqr(Sigma->BlockSize)];
          }
        }
      }
    }
  }
  else
  {  // Sigma and GMat complex
    for(int k=0;k<Sigma->N;k++)
    {
      for(int l=0;l<G0->N;l++)
      {
        for(int i=0;i<G0->BlockSize;i++)
        {
          for(int j=0;j<G0->BlockSize;j++)
          {
            GMat.V[ i * GMat.Nj + j + (GMat.Nj+1) * ((k + l * (Sigma->N+1)) * Sigma->BlockSize )                              ].Re = Sigma->a[j+i*Sigma->BlockSize + k*Sqr(Sigma->BlockSize)].Re;
            GMat.V[ i * GMat.Nj + j + (GMat.Nj+1) * ((k + l * (Sigma->N+1)) * Sigma->BlockSize ) + Sigma->BlockSize           ].Re = Sigma->b[j+i*Sigma->BlockSize + k*Sqr(Sigma->BlockSize)].Re;
            GMat.V[ i * GMat.Nj + j + (GMat.Nj+1) * ((k + l * (Sigma->N+1)) * Sigma->BlockSize ) + Sigma->BlockSize * GMat.Nj ].Re = Sigma->b[i+j*Sigma->BlockSize + k*Sqr(Sigma->BlockSize)].Re;
            GMat.V[ i * GMat.Nj + j + (GMat.Nj+1) * ((k + l * (Sigma->N+1)) * Sigma->BlockSize )                              ].Im = Sigma->a[j+i*Sigma->BlockSize + k*Sqr(Sigma->BlockSize)].Im;
            GMat.V[ i * GMat.Nj + j + (GMat.Nj+1) * ((k + l * (Sigma->N+1)) * Sigma->BlockSize ) + Sigma->BlockSize           ].Im = Sigma->b[j+i*Sigma->BlockSize + k*Sqr(Sigma->BlockSize)].Im;
            GMat.V[ i * GMat.Nj + j + (GMat.Nj+1) * ((k + l * (Sigma->N+1)) * Sigma->BlockSize ) + Sigma->BlockSize * GMat.Nj ].Im =-Sigma->b[i+j*Sigma->BlockSize + k*Sqr(Sigma->BlockSize)].Im;
          }
        }
      }
    }
    {
      int k=Sigma->N;
      for(int l=0;l<G0->N;l++)
      {
        for(int i=0;i<G0->BlockSize;i++)
        {
          for(int j=0;j<G0->BlockSize;j++)
          {
            GMat.V[ i * GMat.Nj + j + (GMat.Nj+1) * ((k + l * (Sigma->N+1)) * Sigma->BlockSize ) ].Re = Sigma->a[j+i*Sigma->BlockSize + k*Sqr(Sigma->BlockSize)].Re;
            GMat.V[ i * GMat.Nj + j + (GMat.Nj+1) * ((k + l * (Sigma->N+1)) * Sigma->BlockSize ) ].Im = Sigma->a[j+i*Sigma->BlockSize + k*Sqr(Sigma->BlockSize)].Im;
          }
        }
      }
    }
  }
  // add G0 to GMat
  if(GMat.QComplex==0)
  { // all real
    for(int l=0;l<G0->N-1;l++)
    {
      for(int i=0;i<G0->BlockSize;i++)
      {
        for(int j=0;j<G0->BlockSize;j++)
        {
          GMat.VRe[ i * GMat.Nj + j + l * (GMat.Nj+1) * (Sigma->N+1) * Sigma->BlockSize                                             ] += G0->aRe[j+i*G0->BlockSize + (l+1)*Sqr(G0->BlockSize)];
          GMat.VRe[ i * GMat.Nj + j + l * (GMat.Nj+1) * (Sigma->N+1) * Sigma->BlockSize +           (Sigma->N+1) * Sigma->BlockSize ]  = G0->bRe[j+i*G0->BlockSize + (l+1)*Sqr(G0->BlockSize)];
          GMat.VRe[ i * GMat.Nj + j + l * (GMat.Nj+1) * (Sigma->N+1) * Sigma->BlockSize + GMat.Nj * (Sigma->N+1) * Sigma->BlockSize ]  = G0->bRe[i+j*G0->BlockSize + (l+1)*Sqr(G0->BlockSize)];
        }
      }
    }
    {
      int l=G0->N-1;
      for(int i=0;i<G0->BlockSize;i++)
      {
        for(int j=0;j<G0->BlockSize;j++)
        {
          GMat.VRe[ i * GMat.Nj + j + l * (GMat.Nj+1) * (Sigma->N+1) * Sigma->BlockSize                                             ] += G0->aRe[j+i*G0->BlockSize + (l+1)*Sqr(G0->BlockSize)];
        }
      }
    }
  }
  else if(G0->QComplex==0)
  { // G0 real but GMat complex
    for(int l=0;l<G0->N-1;l++)
    {
      for(int i=0;i<G0->BlockSize;i++)
      {
        for(int j=0;j<G0->BlockSize;j++)
        {
          GMat.V[ i * GMat.Nj + j + l * (GMat.Nj+1) * (Sigma->N+1) * Sigma->BlockSize                                             ].Re += G0->aRe[j+i*G0->BlockSize + (l+1)*Sqr(G0->BlockSize)];
          GMat.V[ i * GMat.Nj + j + l * (GMat.Nj+1) * (Sigma->N+1) * Sigma->BlockSize +           (Sigma->N+1) * Sigma->BlockSize ].Re  = G0->bRe[j+i*G0->BlockSize + (l+1)*Sqr(G0->BlockSize)];
          GMat.V[ i * GMat.Nj + j + l * (GMat.Nj+1) * (Sigma->N+1) * Sigma->BlockSize + GMat.Nj * (Sigma->N+1) * Sigma->BlockSize ].Re  = G0->bRe[i+j*G0->BlockSize + (l+1)*Sqr(G0->BlockSize)];
        }
      }
    }
    {
      int l=G0->N-1;
      for(int i=0;i<G0->BlockSize;i++)
      {
        for(int j=0;j<G0->BlockSize;j++)
        {
          GMat.V[ i * GMat.Nj + j + l * (GMat.Nj+1) * (Sigma->N+1) * Sigma->BlockSize                                             ].Re += G0->aRe[j+i*G0->BlockSize + (l+1)*Sqr(G0->BlockSize)];
        }
      }
    }
    
  }
  else
  {  // G0 and GMat complex
    for(int l=0;l<G0->N-1;l++)
    {
      for(int i=0;i<G0->BlockSize;i++)
      {
        for(int j=0;j<G0->BlockSize;j++)
        {
          GMat.V[ i * GMat.Nj + j + l * (GMat.Nj+1) * (Sigma->N+1) * Sigma->BlockSize                                             ].Re += G0->a[j+i*G0->BlockSize + (l+1)*Sqr(G0->BlockSize)].Re;
          GMat.V[ i * GMat.Nj + j + l * (GMat.Nj+1) * (Sigma->N+1) * Sigma->BlockSize +           (Sigma->N+1) * Sigma->BlockSize ].Re  = G0->b[j+i*G0->BlockSize + (l+1)*Sqr(G0->BlockSize)].Re;
          GMat.V[ i * GMat.Nj + j + l * (GMat.Nj+1) * (Sigma->N+1) * Sigma->BlockSize + GMat.Nj * (Sigma->N+1) * Sigma->BlockSize ].Re  = G0->b[i+j*G0->BlockSize + (l+1)*Sqr(G0->BlockSize)].Re;
          GMat.V[ i * GMat.Nj + j + l * (GMat.Nj+1) * (Sigma->N+1) * Sigma->BlockSize                                             ].Im += G0->a[j+i*G0->BlockSize + (l+1)*Sqr(G0->BlockSize)].Im;
          GMat.V[ i * GMat.Nj + j + l * (GMat.Nj+1) * (Sigma->N+1) * Sigma->BlockSize +           (Sigma->N+1) * Sigma->BlockSize ].Im  = G0->b[j+i*G0->BlockSize + (l+1)*Sqr(G0->BlockSize)].Im;
          GMat.V[ i * GMat.Nj + j + l * (GMat.Nj+1) * (Sigma->N+1) * Sigma->BlockSize + GMat.Nj * (Sigma->N+1) * Sigma->BlockSize ].Im  =-G0->b[i+j*G0->BlockSize + (l+1)*Sqr(G0->BlockSize)].Im;
        }
      }
    }
    {
      int l=G0->N-1;
      for(int i=0;i<G0->BlockSize;i++)
      {
        for(int j=0;j<G0->BlockSize;j++)
        {
          GMat.VRe[ i * GMat.Nj + j + l * (GMat.Nj+1) * (Sigma->N+1) * Sigma->BlockSize                                             ] += G0->aRe[j+i*G0->BlockSize + (l+1)*Sqr(G0->BlockSize)];
        }
      }
    }
    
  }
  HermitianCompactMatrixToBlockListOfPoles(&GMat, G, G0->BlockSize);
  FreeCompactMatrix(&GMat);
  G->mu = G0->mu;
  return 0;
}

/*  we represent G0 as a TriDiagonal matrix
 *  We represent Sigma "inverse" as an anderson matrix within a compact matrix.
 *      the compact matrix has an aditional dimension inlcuding BlockSize * NGridPoints space
 *
 *  Starting at the botom of the tri-diagonal matrix of G0 we have
 *
 *
 *   SigmaNow = ( G0.An + S.A0 , S.B1, S.B2, S.B3, ... )
 *              (         S.B1 , S.A1, 0   , 0   , ... )
 *              (         S.B2 , 0   , S.A2, 0   , ... )
 *              (         S.B3 , 0   , 0   , S.A3, ... )
 *              (         ...  , ... , ... , ... , ... )
 *   We trun this matrix into a BlockListOfPoles and that into an Anderson matrix to get
 *   the next iteration of SigmaNow
 *   SigmaNow = ( G0.An-1 + S.A0 , S.B1, S.B2, S.B3, ..., G0.Bn , )
 *              (           S.B1 , S.A1, 0   , 0   , ..., 0     , )
 *              (           S.B2 , 0   , S.A2, 0   , ..., 0     , )
 *              (           S.B3 , 0   , 0   , S.A3, ..., 0     , )
 *              (           ...  , ... , ... , ... , ..., 0     , )
 *              (           G0.Bn, 0   , 0   , 0   , 0  , A1    , B1, B2, B3, ... )
 *              (           0    , 0   , 0   , 0   , 0  , B1    , A2, 0 , 0 , ... )
 *              (           0    , 0   , 0   , 0   , 0  , B2    , 0 , A3, 0 , ... )
 *              (           0    , 0   , 0   , 0   , 0  , B3    , 0 , 0 , A4, ... )
 *              (           0    , 0   , 0   , 0   , 0  , ...   ,...,...,..., ... )
 */
int ResponseFunctionAddSelfEnergyOnGridViaG0Tri(const ResponsefunctionType *const G0, const ResponsefunctionType *const Sigma, CompactVectorType *const Grid, ResponsefunctionType * G, double Epsilon)
{
  CompactMatrixType GMat;
  uint32_t SigmaDim;
  if(G0->Block!=Sigma->Block)
  {
    printf("Error in ResponseFunctionAddSelfEnergyOnGridViaG0Tri\n Sigma and G0 not both Block or single valued");
    return 1;
  }
  if(G0->Block==true)
  {
    // set G0Tri to a pointer representing G0 as a tri diagonal matrix, either within G0 or newly allcated and transformed
    BlockTriDiagonalMatrixType *G0Data;
    if((G0->type=='t')||(G0->type=='T'))
    {
      G0Data = (BlockTriDiagonalMatrixType*)G0->data;
    }
    else
    {
      G0Data = (BlockTriDiagonalMatrixType*)malloc(sizeof(BlockTriDiagonalMatrixType));
      switch (G0->type|32)
      {
        case 'l':
        {
          BlockListOfPolesToBlockTridiagonalMatrix(*(BlockListOfPoles*)G0->data, G0Data, 1+((BlockListOfPoles*)G0->data)->N, Epsilon);
        }
        break;
        case 'a':
        {
          BlockAndersonMatrixToBlockTridiagonalMatrix((BlockAndersonMatrixType*)G0->data, G0Data, Epsilon);
        }
        break;
        case 'n':
        {
          BlockNaturalImpurityOrbitalMatrixToBlockTridiagonalMatrix((BlockNaturalImpurityOrbitalMatrixType*)G0->data, G0Data, Epsilon);
        }
        break;
        default:
        {
          printf("Error in ResponseFunctionAddSelfEnergyDenseViaTri unrecognized type to convert to: %c",G0->type);
        }
        break;
      }
    }
    // set GMat to the correct dimension and fill the upper left corner with 1/Sigma
    switch (Sigma->type|32)
    {
      case 't':
      {
        BlockTriDiagonalMatrixType* SData = (BlockTriDiagonalMatrixType*)Sigma->data;
        if(SData->BlockSize != G0Data->BlockSize)
        {
          printf("Error in ResponseFunctionAddSelfEnergyOnGridViaG0Tri\n BlockSize of G0 and Sigma not equal");
        }
        SigmaDim = SData->BlockSize * (SData->N+1);
        GMat.Ni = G0Data->BlockSize * Grid->Ni + SigmaDim;
        GMat.Nj = GMat.Ni;
        GMat.QComplex = G0Data->QComplex | SData->QComplex;
        if(InitCompactMatrix(&GMat))
        {
          printf("Error InitCompactMatrix failed in ResponseFunctionAddSelfEnergyOnGridViaG0Tri\n");
          return 1;
        }
        // fill the matrix
        if(GMat.QComplex==0)
        {
          for(int k=0;k<SData->N;k++)
          {
            for(int i=0;i<SData->BlockSize;i++)
            {
              for(int j=0;j<SData->BlockSize;j++)
              {
                GMat.VRe[j+i*GMat.Nj + k*(GMat.Nj+1)*(SData->BlockSize)                           ] = SData->aRe[j+i*SData->BlockSize + k*Sqr(SData->BlockSize)];
                GMat.VRe[j+i*GMat.Nj + k*(GMat.Nj+1)*(SData->BlockSize) + SData->BlockSize        ] = SData->bRe[j+i*SData->BlockSize + k*Sqr(SData->BlockSize)];
                GMat.VRe[j+i*GMat.Nj + k*(GMat.Nj+1)*(SData->BlockSize) + SData->BlockSize*GMat.Nj] = SData->bRe[i+j*SData->BlockSize + k*Sqr(SData->BlockSize)];
              }
            }
          }
          {
            int k=SData->N;
            for(int i=0;i<SData->BlockSize;i++)
            {
              for(int j=0;j<SData->BlockSize;j++)
              {
                GMat.VRe[j+i*GMat.Nj + k*(GMat.Nj+1)*(SData->BlockSize)                           ] = SData->aRe[j+i*SData->BlockSize + k*Sqr(SData->BlockSize)];
              }
            }
          }
        }
        else if (SData->QComplex==0)
        { // GMat complex but SData real
          printf("TODO");
        }
        else
        { // GMat and SData both complex
          printf("TODO");
        }
      }
      break;
      case 'a':
      {
        printf("TODO");
      }
      break;
      case 'n':
      {
        printf("TODO");
      }
      break;
      case 'l':
      {
        // transform Sigma to tridiagonal type and call the function again.
        ResponsefunctionType SigmaTri;
        ConvertResponsefunction(Sigma, &SigmaTri, 't', 0, Epsilon);
        // We already transformed G0ResTri to tridiagonal format so lets reuse this
        ResponsefunctionType G0ResTri;
        G0ResTri.type='t';
        G0ResTri.data=G0Data;
        G0ResTri.Block=true;
        G0ResTri.NSuggestedGridPoints=0;
        G0ResTri.SuggestedGrid=NULL;
        // and now restart the function
        ResponseFunctionAddSelfEnergyOnGridViaG0Tri(&G0ResTri, &SigmaTri, Grid, G, Epsilon);
        FreeResponsefunction(&SigmaTri);
        // we need to free G0 if this was not called as a tridiagonal format
        if(!((G0->type=='t')||(G0->type=='T')))
        {
          FreeBlockTridiagonalMatrix(G0Data);
          free(G0Data);
        }
        return 0;
      }
      break;
    }
    G->Block=true;
    G->type='l';
    G->NSuggestedGridPoints=0;
    G->SuggestedGrid=NULL;
    G->data = malloc(sizeof(BlockListOfPoles));
    ResponseFunctionAddSelfEnergyOnGridBlockTriCompactLP(G0Data, &GMat, SigmaDim, Grid, (BlockListOfPoles * )G->data);
    if(!((G0->type=='t')||(G0->type=='T')))
    {
      FreeBlockTridiagonalMatrix(G0Data);
      free(G0Data);
    }
    FreeCompactMatrix(&GMat);
  }
  else // non block function
  {
    printf("TODO");
  }
  return 0;
}
int ResponseFunctionAddSelfEnergyOnGridTriCompactLP(const TriDiagonalMatrixType *const G0, CompactMatrixType * Sigma, uint32_t SigmaDim, CompactVectorType *const Grid, ListOfPoles * G)
{
  //TODO
  return 0;
}
/*   The top left SigmaDim by SigmaDim block of Sigma contains the self energy
 *   G0 is a tri-diagonal matrix
 *   Grid the energy grid we want to use
 *   We now need to calculate G iteratively from botom to top of the Tri-Diagonal matrix
 */
int ResponseFunctionAddSelfEnergyOnGridBlockTriCompactLP(const BlockTriDiagonalMatrixType *const G0, CompactMatrixType * Sigma, uint32_t SigmaDim, CompactVectorType *const Grid, BlockListOfPoles * G)
{
  if(Grid->QComplex!=0)
  {
    printf("Error in ResponseFunctionAddSelfEnergyOnGridBlockTriCompactLP Grid is complex");
    return 1;
  }
  if(Sigma->QComplex==0)
  {
    return ResponseFunctionAddSelfEnergyOnGridBlockTriCompactLPRealReal(G0, Sigma, SigmaDim, Grid, G);
  }
  else
  {
    if (G0->QComplex==0)
    {
      return ResponseFunctionAddSelfEnergyOnGridBlockTriCompactLPRealComplex(G0, Sigma, SigmaDim, Grid, G);
    }
    else
    {
      return ResponseFunctionAddSelfEnergyOnGridBlockTriCompactLPComplexComplex(G0, Sigma, SigmaDim, Grid, G);
    }
  }
}
int ResponseFunctionAddSelfEnergyOnGridBlockTriCompactLPRealReal(const BlockTriDiagonalMatrixType *const G0, CompactMatrixType * Sigma, uint32_t SigmaDim, CompactVectorType *const Grid, BlockListOfPoles * G)
{
  BlockListOfPoles SigmaTemp;
  BlockTriDiagonalMatrixType Gminor;
  if(G0->N==0)
  {
    G->BlockSize=G0->BlockSize;
    G->N = G0->N;
    G->QComplex = G0->QComplex;
    InitBlockListOfPoles(G);
    memcpy(G->a0Re, G0->aRe, Sqr(G->BlockSize)*sizeof(double));
    return 0;
  }
  // we will overwrite the top left BlockSize by BlockSize part of Sigma
  // we copy it to a new matrix and then can refill it whenever needed
  CompactMatrixType SigmaA0;
  SigmaA0.Ni       = G0->BlockSize;
  SigmaA0.Nj       = G0->BlockSize;
  SigmaA0.QComplex = 0;
  InitCompactMatrix(&SigmaA0);
  for(int i=0;i<SigmaA0.Ni;i++)
  {
    for(int j=0;j<SigmaA0.Nj;j++)
    {
      SigmaA0.VRe[j+i*SigmaA0.Nj] = Sigma->VRe[j+i*Sigma->Nj];
    }
  }
  // The first step can be done different from the others
  // We have a matrix Sigma of size Sigma->Nj, this can fit Nstart self energies of size SigmaDim,
  // with Nstart = floor(Sigma->Nj / SigmaDim)
  int Nstart = Min(Sigma->Nj / SigmaDim,G0->N);
  { // first step is different
    // copy sigma
    for(int k=1;k<Nstart;k++)
    {
      for(int i=0;i<SigmaDim;i++)
      {
        for(int j=0;j<SigmaDim;j++)
        {
          Sigma->VRe[j+i*Sigma->Nj + k*(SigmaDim*(1+Sigma->Nj))]=Sigma->VRe[j+i*Sigma->Nj];
        }
      }
    }
    // add G0 to Sigma
    // We can add the last Nstart elements of G0 to the matrix
    for(int k=0;k<Nstart-1;k++)
    {
      for(int i=0;i<G0->BlockSize;i++)
      {
        for(int j=0;j<G0->BlockSize;j++)
        {
          Sigma->VRe[j+i*Sigma->Nj + k*(SigmaDim*(1+Sigma->Nj))]                      += G0->aRe[j+i*G0->BlockSize + (G0->N + k - Nstart + 1) * Sqr(G0->BlockSize)];
          Sigma->VRe[j+i*Sigma->Nj + k*(SigmaDim*(1+Sigma->Nj)) + SigmaDim]           += G0->bRe[j+i*G0->BlockSize + (G0->N + k - Nstart + 1) * Sqr(G0->BlockSize)];
          Sigma->VRe[j+i*Sigma->Nj + k*(SigmaDim*(1+Sigma->Nj)) + SigmaDim*Sigma->Nj] += G0->bRe[i+j*G0->BlockSize + (G0->N + k - Nstart + 1) * Sqr(G0->BlockSize)];
        }
      }
    }
    {
      int k = Nstart-1;
      for(int i=0;i<G0->BlockSize;i++)
      {
        for(int j=0;j<G0->BlockSize;j++)
        {
          Sigma->VRe[j+i*Sigma->Nj + k*(SigmaDim*(1+Sigma->Nj))] += G0->aRe[j+i*G0->BlockSize + (G0->N + k - Nstart + 1) * Sqr(G0->BlockSize)];
        }
      }
    }
    // Nstart blocks of sigma set diagonalize the matrix of dimension Nstart * SigmaDim
    HermitianCompactMatrixToBlockListOfPoles(Sigma, &SigmaTemp, G0->BlockSize);
    BlockListOfPolesPlaceOnMesh(&SigmaTemp, Grid->VRe, Grid->Ni);
  }
  // set botom block of Sigma to zero
  for(int i=0;i<SigmaDim;i++)
  {
    for(int j=SigmaDim;j<Sigma->Nj;j++)
    {
      Sigma->VRe[j+i*Sigma->Nj]=0;
    }
  }
  for(int i=SigmaDim;i<Sigma->Ni;i++)
  {
    for(int j=0;j<Sigma->Nj;j++)
    {
      Sigma->VRe[j+i*Sigma->Nj]=0;
    }
  }
  for(int k=G0->N - Nstart; k>0;k--)
  {
    // we placed the botom Nstart block of G0 with self energy on a mesh into SigmaTemp
    // we now need to add layer i
    // set left top block of Sigma to G0->Ak - SigmaA0
    // set hopping from level k to all below to G0->Bk
    for(int i=0;i<G0->BlockSize;i++)
    {
      for(int j=0;j<G0->BlockSize;j++)
      {
        Sigma->VRe[j+i*Sigma->Nj] = G0->aRe[j+i*G0->BlockSize + k * Sqr(G0->BlockSize)]  + SigmaA0.VRe[j+i*SigmaA0.Nj];
        Sigma->VRe[j+i*Sigma->Nj + SigmaDim] = G0->bRe[j+i*G0->BlockSize + k * Sqr(G0->BlockSize)];
        Sigma->VRe[j+(i+SigmaDim)*Sigma->Nj] = G0->bRe[i+j*G0->BlockSize + k * Sqr(G0->BlockSize)];
      }
    }
    // set the botom part to SigmaTemp. A0 = 0 and B0=1 by construction
    BlockListOfPolesToBlockTridiagonalMatrix(SigmaTemp, &Gminor, SigmaTemp.N, DBL_EPSILON);
    for (int l=0;l<Gminor.N-1;l++)
    {
      for(int i=0;i<G0->BlockSize;i++)
      {
        for(int j=0;j<G0->BlockSize;j++)
        {
          Sigma->VRe[j+i*Sigma->Nj+ (SigmaDim + G0->BlockSize * l) * (1+Sigma->Nj)                            ] = Gminor.aRe[(l+1)*Sqr(G0->BlockSize) + i*G0->BlockSize + j];
          Sigma->VRe[j+i*Sigma->Nj+ (SigmaDim + G0->BlockSize * l) * (1+Sigma->Nj) + G0->BlockSize            ] = Gminor.bRe[(l+1)*Sqr(G0->BlockSize) + i*G0->BlockSize + j];
          Sigma->VRe[j+i*Sigma->Nj+ (SigmaDim + G0->BlockSize * l) * (1+Sigma->Nj) + G0->BlockSize * Sigma->Nj] = Gminor.bRe[(l+1)*Sqr(G0->BlockSize) + j*G0->BlockSize + i];
        }
      }
    }
    {
      int l=Gminor.N-1;
      for(int i=0;i<G0->BlockSize;i++)
      {
        for(int j=0;j<G0->BlockSize;j++)
        {
          Sigma->VRe[j+i*Sigma->Nj+ (SigmaDim + G0->BlockSize * l) * (1+Sigma->Nj)]=Gminor.aRe[(l+1)*Sqr(G0->BlockSize) + i*G0->BlockSize + j];
        }
      }
    }
    for (int i= (SigmaDim + G0->BlockSize * Gminor.N); i < Sigma->Ni; i++)
    {
      for (int j= (SigmaDim + G0->BlockSize * Gminor.N); j < Sigma->Nj; j++)
      {
        Sigma->VRe[j+i*Sigma->Nj]=0;
      }
    }
    FreeBlockTridiagonalMatrix(&Gminor);
    FreeBlockListOfPoles(&SigmaTemp);
    HermitianCompactMatrixToBlockListOfPoles(Sigma, &SigmaTemp, G0->BlockSize);
    BlockListOfPolesPlaceOnMesh(&SigmaTemp, Grid->VRe, Grid->Ni);
  }
  // We placed G0 with self energy on a mesh, We now need to add G->A0 and multiply with G0->B0
  CompactMatrixType B0;
  B0.QComplex = G0->QComplex;
  B0.Ni = G0->BlockSize;
  B0.Nj = G0->BlockSize;
  B0.Name = (char*)G0->Name;
  if(G0->QComplex==0)
  {
    B0.VRe = G0->bRe;
  }
  else
  {
    B0.V = G0->b;
  }
  BlockListOfPolesTransformcT(&SigmaTemp, G, &B0);
  if(SigmaTemp.QComplex==0)
  {
    for(int i=0;i<G0->BlockSize;i++)
    {
      for(int j=0;j<G0->BlockSize;j++)
      {
        SigmaTemp.a0Re[j+i*SigmaTemp.BlockSize] = G0->aRe[j+i*G0->BlockSize];
      }
    }
  }
  else if (G0->QComplex==0)
  {
    for(int i=0;i<G0->BlockSize;i++)
    {
      for(int j=0;j<G0->BlockSize;j++)
      {
        SigmaTemp.a0[j+i*SigmaTemp.BlockSize].Re = G0->aRe[j+i*G0->BlockSize];
        SigmaTemp.a0[j+i*SigmaTemp.BlockSize].Im = 0;
      }
    }
  }
  else
  {
    for(int i=0;i<G0->BlockSize;i++)
    {
      for(int j=0;j<G0->BlockSize;j++)
      {
        SigmaTemp.a0[j+i*SigmaTemp.BlockSize].Re = G0->a[j+i*G0->BlockSize].Re;
        SigmaTemp.a0[j+i*SigmaTemp.BlockSize].Im = G0->a[j+i*G0->BlockSize].Im;
      }
    }
  }
  FreeBlockListOfPoles(&SigmaTemp);
  FreeCompactMatrix(&SigmaA0);
  return 0;
}
int ResponseFunctionAddSelfEnergyOnGridBlockTriCompactLPRealComplex(const BlockTriDiagonalMatrixType *const G0, CompactMatrixType * Sigma, uint32_t SigmaDim, CompactVectorType *const Grid, BlockListOfPoles * G)
{
  printf("TODO\n");
  return 1;
}
int ResponseFunctionAddSelfEnergyOnGridBlockTriCompactLPComplexComplex(const BlockTriDiagonalMatrixType *const G0, CompactMatrixType * Sigma, uint32_t SigmaDim, CompactVectorType *const Grid, BlockListOfPoles * G)
{
  printf("TODO\n");
  return 1;
}


/*   We have G0 = A_0 + B_0^* ( 1 / (w - A1 - sum_k  B_k  1/(w-A_{k+1}) B_k^+ ) B_0^T
 *           S  = As_0 + sum_i Bs_i (1/as_i) Bs_i^+  // note Bs_i Bs_i^+ is stored
 *
 *   We need G = A_0 + B_0^* ( 1 / (w - A1 - As_0 - sum_i Bs_i (1/as_i) Bs_i^+   - sum_k  B_k  1/(w-A_{k+1} - As_0 - sum_i Bs_i (1/as_i) Bs_i^+ ) B_k^+ ) B_0^T
 *
 *   for each k we generate the Anderson matrix
 *      B_k  1/(w-A_{k+1} - As_0 - sum_i Bs_i (1/as_i) Bs_i^+ ) B_k^+
 *   and transform it to a list of poles that is added to G in the denominator
 *
 *  Note that the Ak's might not be on the original basis that is used for sigma. If this is the case the funcitons return the incorrect result
 *
 */
int ResponseFunctionAddSelfEnergyViaAndLP(const ResponsefunctionType *const G0, const ResponsefunctionType *const Sigma, ResponsefunctionType * G, double Epsilon, CompactVectorType *const Grid)
{
  if(G0->Block!=Sigma->Block)
  {
    printf("Error in ResponseFunctionAddSelfEnergyViaAndLP\n G0 and Sigma must either be both of Block type, or both of single value type\n");
    return 1;
  }
  if(G0->Block==true)
  {
    G->Block=true;
    G->data=malloc(sizeof(BlockAndersonMatrixType));
    G->type='a';
    G->NSuggestedGridPoints=0;
    G->SuggestedGrid=NULL;
    BlockAndersonMatrixType *G0Data;
    BlockListOfPoles        *SigmaData;
    if((G0->type=='a')||(G0->type=='A'))
    {
      G0Data = (BlockAndersonMatrixType*)G0->data;
    }
    else
    {
      G0Data = (BlockAndersonMatrixType*)malloc(sizeof(BlockAndersonMatrixType));
      switch (G0->type|32)
      {
        case 'l':
        {
          BlockListOfPolesToBlockAndersonMatrix(*(BlockListOfPoles*)G0->data, G0Data, 1+((BlockListOfPoles*)G0->data)->N);
        }
        break;
        case 't':
        {
          BlockTriDiagonalMatrixToBlockAndersonMatrix(*(BlockTriDiagonalMatrixType*)G0->data, G0Data);
        }
        break;
        case 'n':
        {
          BlockNaturalImpurityOrbitalMatrixToBlockAndersonMatrix((BlockNaturalImpurityOrbitalMatrixType*)G0->data, G0Data);
        }
        break;
        default:
        {
          printf("Error in ResponseFunctionAddSelfEnergyViaAndLP unrecognized type to convert to: %c",G0->type);
        }
        break;
      }
    }
    if((Sigma->type=='l')||(Sigma->type=='L'))
    {
      SigmaData = (BlockListOfPoles*)Sigma->data;
    }
    else
    {
      SigmaData = (BlockListOfPoles*)malloc(sizeof(BlockListOfPoles));
      switch (Sigma->type|32)
      {
        case 'a':
        {
          BlockAndersonMatrixToBlockListOfPoles((BlockAndersonMatrixType*)Sigma->data, SigmaData);
        }
        break;
        case 't':
        {
          BlockTriDiagonalMatrixToBlockListOfPoles((BlockTriDiagonalMatrixType*)Sigma->data, SigmaData);
        }
        break;
        case 'n':
        {
          BlockNaturalImpurityOrbitalMatrixToBlockListOfPoles((BlockNaturalImpurityOrbitalMatrixType*)Sigma->data, SigmaData);
        }
        break;
        default:
        {
          printf("Error in ResponseFunctionAddSelfEnergyViaAndLP unrecognized type to convert to: %c",Sigma->type);
        }
        break;
      }
    }
    if((Grid==NULL)||Grid->Ni==0)
    {
      if(ResponseFunctionAddSelfEnergyBlockAndLPAnd(G0Data, SigmaData, (BlockAndersonMatrixType*)G->data, Epsilon))
      {
        printf("Error in ResponseFunctionAddSelfEnergyViaAndLP\n function ResponseFunctionAddSelfEnergyBlockAndLPAnd failed\n");
      }
    }
    else
    {
      if(ResponseFunctionAddSelfEnergyOnGridBlockAndLPAnd(G0Data, SigmaData, (BlockAndersonMatrixType*)G->data, Epsilon, Grid))
      {
        printf("Error in ResponseFunctionAddSelfEnergyViaAndLP\n function ResponseFunctionAddSelfEnergyOnGridBlockAndLPAnd failed\n");
      }
    }
    if(!((G0->type=='a')||(G0->type=='A')))
    {
      FreeBlockAndersonMatrix(G0Data);
      free(G0Data);
    }
    if(!((Sigma->type=='l')||(Sigma->type=='L')))
    {
      FreeBlockListOfPoles(SigmaData);
      free(SigmaData);
    }
  }
  else // non block version
  {
    G->Block=false;
    G->data=malloc(sizeof(AndersonMatrixType));
    G->type='a';
    G->NSuggestedGridPoints=0;
    G->SuggestedGrid=NULL;
    AndersonMatrixType *G0Data;
    ListOfPoles        *SigmaData;
    if((G0->type=='a')||(G0->type=='A'))
    {
      G0Data = (AndersonMatrixType*)G0->data;
    }
    else
    {
      G0Data = (AndersonMatrixType*)malloc(sizeof(AndersonMatrixType));
      switch (G0->type|32)
      {
        case 'l':
        {
          ListOfPolesToAndersonMatrix(*(ListOfPoles*)G0->data, G0Data);
        }
        break;
        case 't':
        {
          TriDiagonalMatrixToAndersonMatrix(*(TriDiagonalMatrixType*)G0->data, G0Data);
        }
        break;
        case 'n':
        {
          BADoubleTriDiagonalMatrixToAndersonMatrix((BADoubleTriDiagonalMatrixType*)G0->data, G0Data);
        }
        break;
        default:
        {
          printf("Error in ResponseFunctionAddSelfEnergyViaAndLP unrecognized type to convert to: %c",G0->type);
        }
        break;
      }
    }
    if((Sigma->type=='l')||(Sigma->type=='L'))
    {
      SigmaData = (ListOfPoles*)Sigma->data;
    }
    else
    {
      SigmaData = (ListOfPoles*)malloc(sizeof(BlockListOfPoles));
      switch (Sigma->type|32)
      {
        case 'a':
        {
          AndersonMatrixToListOfPoles(*(AndersonMatrixType*)Sigma->data, SigmaData);
        }
        break;
        case 't':
        {
          TriDiagonalMatrixToListOfPoles(*(TriDiagonalMatrixType*)Sigma->data, SigmaData);
        }
        break;
        case 'n':
        {
          BADoubleTriDiagonalMatrixToListOfPoles((BADoubleTriDiagonalMatrixType*)Sigma->data, SigmaData);
        }
        break;
        default:
        {
          printf("Error in ResponseFunctionAddSelfEnergyViaAndLP unrecognized type to convert to: %c",Sigma->type);
        }
        break;
      }
    }
    if((Grid==NULL)||Grid->Ni==0)
    {
      if(ResponseFunctionAddSelfEnergyAndLPAnd(G0Data, SigmaData, (AndersonMatrixType*)G->data, Epsilon))
      {
        printf("Error in ResponseFunctionAddSelfEnergyViaAndLP\n function ResponseFunctionAddSelfEnergyAndLPAnd failed\n");
      }
    }
    else
    {
      if(ResponseFunctionAddSelfEnergyOnGridAndLPAnd(G0Data, SigmaData, (AndersonMatrixType*)G->data, Epsilon, Grid))
      {
        printf("Error in ResponseFunctionAddSelfEnergyViaAndLP\n function ResponseFunctionAddSelfEnergyOnGridAndLPAnd failed\n");
      }
    }
    if(!((G0->type=='a')||(G0->type=='A')))
    {
      FreeAndersonMatrix(G0Data);
      free(G0Data);
    }
    if(!((Sigma->type=='l')||(Sigma->type=='L')))
    {
      FreeListOfPoles(SigmaData);
      free(SigmaData);
    }
  }
  return 0;
}
int ResponseFunctionAddSelfEnergyAndLPAnd(const AndersonMatrixType *const G0, const ListOfPoles *const Sigma, AndersonMatrixType * G, double Epsilon)
{
  if(G0->N==0)
  {
    AndersonMatrixCopy(G, *G0);
    return 0;
  }
  uint32_t i, j;
  ListOfPoles bSigmaPole;
  AndersonMatrixType bSigma;
  G->N = G0->N * (Sigma->N + 1) + 1;
  if (InitAndersonMatrix(G))
  {
    printf("InitAndersonMatrix failed in ResponseFunctionAddSelfEnergyAndLPAnd\n");
    return 1;
  }
  strcpy(G->Name, "GBath");
  bSigma.N = Sigma->N + 1;
  if (InitAndersonMatrix(&bSigma))
  {
    printf("InitAndersonMatrix failed in ResponseFunctionAddSelfEnergyAndLPAnd\n");
    return 1;
  }
  bSigma.a[0] = 0;
  for (j = 0; j < Sigma->N; j++)
  {
    bSigma.a[j + 2] = Sigma->a[j];
    bSigma.b[j + 1] = sqrt(Sigma->b[j]);
    G->a[j+2] = bSigma.a[j + 2];
    G->b[j+1] = bSigma.b[j + 1];
  }
  for (i = 1; i < G0->N; i++)
  {
    bSigma.b[0] = G0->b[i];
    bSigma.a[1] = G0->a[i + 1] + Sigma->a0;    if (AndersonMatrixToListOfPoles(bSigma, &bSigmaPole))
    {
      printf("AndersonMatrixToListOfPoles failed in ResponseFunctionAddSelfEnergyAndLPAnd\n");
      return 1;
    }
    for (j = 0; j < bSigmaPole.N; j++)
    {
      G->a[j + i * (Sigma->N + 1) + 1] = bSigmaPole.a[j];
      G->b[j + i * (Sigma->N + 1) + 0] = sqrt(bSigmaPole.b[j]);
    }
    for (j = bSigmaPole.N; j < bSigma.N; j++)
    {
      G->a[j + i * (Sigma->N + 1) + 1] = 0;
      G->b[j + i * (Sigma->N + 1) + 0] = 0;
    }
    FreeListOfPoles(&bSigmaPole);
  }
  G->a[0] = G0->a[0];
  G->b[0] = G0->b[0];
  G->a[1] = G0->a[1] + Sigma->a0;
  G->mu = G0->mu;
  FreeAndersonMatrix(&bSigma);
  OrderAndersonMatrix(G);
  return 0;
}
int ResponseFunctionAddSelfEnergyBlockAndLPAnd(const BlockAndersonMatrixType *const G0, const BlockListOfPoles *const Sigma, BlockAndersonMatrixType * G, double Epsilon)
{
  if(G0->N==0)
  {
    CopyBlockAndersonMatrix(*G0, G);
    return 0;
  }
  int i,j,k,l;
  CompactMatrixType SigmaInv,G0Bi,GBiSqr,GBi;
  BlockListOfPoles G0iWithSigma;
  // write the inverse of Sigma as a matrix
  BlockListOfPolesInversedToCompactMatrix(Sigma, &SigmaInv);
  if((G0->QComplex!=0)&&(Sigma->QComplex==0)) MakeCompactRealMatrixComplex(&SigmaInv);
  // we use G0Bi and GBiSqr as a dummie pointer Compact matrix to elements in G0
  G0Bi.Ni=G0->BlockSize;
  G0Bi.Nj=G0->BlockSize;
  G0Bi.QComplex=G0->QComplex;
  G0Bi.Name=SigmaInv.Name;
  GBiSqr.Ni=G0->BlockSize;
  GBiSqr.Nj=G0->BlockSize;
  GBiSqr.QComplex=SigmaInv.QComplex;
  GBiSqr.Name=SigmaInv.Name;
  // initialize G -- this can be very large
  G->BlockSize = G0->BlockSize;
  G->QComplex = SigmaInv.QComplex;
  G->N = 1 + (G0->N-1) * (Sigma->N+1) * G0->BlockSize + Sigma->N;
  InitBlockAndersonMatrix(G);
  // set G.A0 and G.B0
  if(G->QComplex==0)
  {
    memcpy(G->aRe, G0->aRe, 2*Sqr(G0->BlockSize)*sizeof(double));
    memcpy(G->bRe, G0->bRe, Sqr(G0->BlockSize)*sizeof(double));
  }
  else
  {
    if(G0->QComplex==0)
    {
      for(i=0;i<Sqr(G0->BlockSize);i++)
      {
        G->b[i].Re = G0->bRe[i];
        G->b[i].Im = 0;
      }
      for(i=0;i<2*Sqr(G0->BlockSize);i++)
      {
        G->a[i].Re = G0->aRe[i];
        G->a[i].Im = 0;
      }
    }
    else
    {
      memcpy(G->a, G0->a, 2*Sqr(G0->BlockSize)*sizeof(Complex16));
      memcpy(G->b, G0->b, Sqr(G0->BlockSize)*sizeof(Complex16));
    }
  }
  G->mu=G0->mu;
  G->N=1;
  // add self energy to G->A1
  if(G->QComplex==0)
  {
    for(i=0;i<Sqr(G0->BlockSize);i++)
    {
      G->aRe[i+Sqr(G0->BlockSize)]+=Sigma->a0Re[i];
    }
  }
  else
  {
    if(Sigma->QComplex==0)
    {
      for(i=0;i<Sqr(G0->BlockSize);i++)
      {
        G->a[i+Sqr(G0->BlockSize)].Re+=Sigma->a0Re[i];
      }
    }
    else
    {
      for(i=0;i<Sqr(G0->BlockSize);i++)
      {
        G->a[i+Sqr(G0->BlockSize)].Re+=Sigma->a0[i].Re;
        G->a[i+Sqr(G0->BlockSize)].Im+=Sigma->a0[i].Im;
      }
    }
  }
  // add self energy poles to G
  if(G->QComplex==0)
  {
    for(j=0;j<Sigma->N;j++)
    {
      GBiSqr.VRe = &(Sigma->VRe[1 + j * (Sqr(Sigma->BlockSize)+1)]);
      CompactMatrixSqrt(&GBiSqr, &GBi);
      memcpy(&(G->bRe[G->N * Sqr(G->BlockSize)]), GBi.VRe, Sqr(G->BlockSize)*sizeof(double));
      FreeCompactMatrix(&GBi);
      G->N++;
      for(k=0;k<G->BlockSize;k++)
      {
        G->aRe[G->N*Sqr(G->BlockSize) + k * (1+G->BlockSize)] = Sigma->VRe[j * (Sqr(Sigma->BlockSize)+1)];
      }
    }
  }
  else
  {
    if(Sigma->QComplex==0)
    {
      for(j=0;j<Sigma->N;j++)
      {
        GBiSqr.VRe = &(Sigma->VRe[1 + j * (Sqr(Sigma->BlockSize)+1)]);
        CompactMatrixSqrt(&GBiSqr, &GBi);
        memcpy(&(G->bRe[G->N * Sqr(G->BlockSize)]), GBi.VRe, Sqr(G->BlockSize)*sizeof(double));
        for(k=0;k<Sqr(G->BlockSize);k++)
        {
          G->b[G->N*Sqr(G->BlockSize) + k].Re = GBi.VRe[k];
          G->b[G->N*Sqr(G->BlockSize) + k].Re = 0;
        }
        FreeCompactMatrix(&GBi);
        G->N++;
        for(k=0;k<G->BlockSize;k++)
        {
          G->a[G->N*Sqr(G->BlockSize) + k * (1+G->BlockSize)].Re = Sigma->VRe[j * (Sqr(Sigma->BlockSize)+1)];
          G->a[G->N*Sqr(G->BlockSize) + k * (1+G->BlockSize)].Re = 0;
        }
      }
    }
    else
    {
      for(j=0;j<Sigma->N;j++)
      {
        GBiSqr.V = &(Sigma->V[1 + j * (Sqr(Sigma->BlockSize)+1)]);
        CompactMatrixSqrt(&GBiSqr, &GBi);
        memcpy(&(G->b[G->N * Sqr(G->BlockSize)]), GBi.V, Sqr(G->BlockSize)*sizeof(Complex16));
        FreeCompactMatrix(&GBi);
        G->N++;
        for(k=0;k<G->BlockSize;k++)
        {
          G->a[G->N*Sqr(G->BlockSize) + k * (1+G->BlockSize)].Re = Sigma->V[j * (Sqr(Sigma->BlockSize)+1)].Re;
          G->a[G->N*Sqr(G->BlockSize) + k * (1+G->BlockSize)].Im = Sigma->V[j * (Sqr(Sigma->BlockSize)+1)].Im;
        }
      }
    }
  }
  // loop over all blocks in G0 (i=1;i<G0.N;i++ )
  for(i=1;i<G0->N;i++)
  {
    // Set the left top block of SigmaInv equal to G0.A[i+1] + Sigma.A[0]
    if(SigmaInv.QComplex==0)
    {
      for(k=0;k<G0->BlockSize;k++)
      {
        for(l=0;l<G0->BlockSize;l++)
        {
          SigmaInv.VRe[l + k*SigmaInv.Nj] = G0->aRe[(i+1) * Sqr(G0->BlockSize) + l + k*(G0->BlockSize)] + Sigma->a0Re[l + k*(G0->BlockSize)];
        }
      }
    }
    else
    {
      if(Sigma->QComplex==0)
      {
        for(k=0;k<G0->BlockSize;k++)
        {
          for(l=0;l<G0->BlockSize;l++)
          {
            SigmaInv.V[l + k*SigmaInv.Nj].Re = G0->a[(i+1) * Sqr(G0->BlockSize) + l + k*(G0->BlockSize)].Re + Sigma->a0Re[l + k*(G0->BlockSize)];
            SigmaInv.V[l + k*SigmaInv.Nj].Im = G0->a[(i+1) * Sqr(G0->BlockSize) + l + k*(G0->BlockSize)].Im;
          }
        }
      }
      else if(G0->QComplex==0)
      {
        for(k=0;k<G0->BlockSize;k++)
        {
          for(l=0;l<G0->BlockSize;l++)
          {
            SigmaInv.V[l + k*SigmaInv.Nj].Re = G0->aRe[(i+1) * Sqr(G0->BlockSize) + l + k*(G0->BlockSize)] + Sigma->a0[l + k*(G0->BlockSize)].Re;
            SigmaInv.V[l + k*SigmaInv.Nj].Im =                                                               Sigma->a0[l + k*(G0->BlockSize)].Im;
          }
        }
      }
      else
      { // all complex
        for(k=0;k<G0->BlockSize;k++)
        {
          for(l=0;l<G0->BlockSize;l++)
          {
            SigmaInv.V[l + k*SigmaInv.Nj].Re = G0->a[(i+1) * Sqr(G0->BlockSize) + l + k*(G0->BlockSize)].Re + Sigma->a0[l + k*(G0->BlockSize)].Re;
            SigmaInv.V[l + k*SigmaInv.Nj].Im = G0->a[(i+1) * Sqr(G0->BlockSize) + l + k*(G0->BlockSize)].Im + Sigma->a0[l + k*(G0->BlockSize)].Im;
          }
        }
      }
    }
    // SigmaInv to BlockListOfPoles
    HermitianCompactMatrixToBlockListOfPoles(&SigmaInv, &G0iWithSigma, G0->BlockSize);
    // Merge poles with energy differences smaller than Epsilon
    BlockListOfPolesMergeDegeneratePoles(&G0iWithSigma, Epsilon);
    // rotate BlockListOpPoles by G0.B[i]
    if(G0->QComplex==0)
    {
      G0Bi.VRe = &(G0->bRe[i*Sqr(G0->BlockSize)]);
    }
    else
    {
      G0Bi.V = &(G0->b[i*Sqr(G0->BlockSize)]);
    }
    BlockListOfPolesTransformInlineNC(&G0iWithSigma, &G0Bi);
    // List of Poles to AndersonImpurity matrix part
    if(G->QComplex==0)
    {
      for(j=0;j<G0iWithSigma.N;j++)
      {
        GBiSqr.VRe = &(G0iWithSigma.VRe[1 + j * (Sqr(G0iWithSigma.BlockSize)+1)]);
        CompactMatrixSqrt(&GBiSqr, &GBi);
        memcpy(&(G->bRe[G->N * Sqr(G->BlockSize)]), GBi.VRe, Sqr(G->BlockSize)*sizeof(double));
        FreeCompactMatrix(&GBi);
        G->N++;
        for(k=0;k<G->BlockSize;k++)
        {
          G->aRe[G->N*Sqr(G->BlockSize) + k * (1+G->BlockSize)] = G0iWithSigma.VRe[j * (Sqr(G0iWithSigma.BlockSize)+1)];
        }
      }
    }
    else
    {
      for(j=0;j<G0iWithSigma.N;j++)
      {
        GBiSqr.V = &(G0iWithSigma.V[1 + j * (Sqr(G0iWithSigma.BlockSize)+1)]);
        CompactMatrixSqrt(&GBiSqr, &GBi);
        memcpy(&(G->b[G->N * Sqr(G->BlockSize)]), GBi.V, Sqr(G->BlockSize)*sizeof(double));
        FreeCompactMatrix(&GBi);
        G->N++;
        for(k=0;k<G->BlockSize;k++)
        {
          G->a[G->N*Sqr(G->BlockSize) + k * (1+G->BlockSize)] = G0iWithSigma.V[j * (Sqr(G0iWithSigma.BlockSize)+1)];
        }
      }
    }
    FreeBlockListOfPoles(&G0iWithSigma);
  }
  FreeCompactMatrix(&SigmaInv);
  return 0;
}
int ResponseFunctionAddSelfEnergyOnGridAndLPAnd(const AndersonMatrixType *const G0, const ListOfPoles *const Sigma, AndersonMatrixType * G, double Epsilon, CompactVectorType *const Grid)
{
  if(G0->N==0)
  {
    AndersonMatrixCopy(G, *G0);
    return 0;
  }
  if(Grid->Ni==0)
  {
    printf("Zero gridpoints in ResponseFunctionAddSelfEnergyOnGridAndLPAnd. That does not make sence.\n");
    return 1;
  }
  uint32_t i, j;
  ListOfPoles bSigmaPole,Ghyb;
  AndersonMatrixType bSigma;
  // we take the hybridization part as a list of poles on a grid to which we add
  // once all added we bring it to an anderson form
  Ghyb.N=Grid->Ni;
  InitListOfPoles(&Ghyb);
  for(i=0;i<Ghyb.N;i++)
  {
    Ghyb.a[i]=Grid->VRe[i];
    Ghyb.b[i]=0;
  }
  AddToListOfPolesOnGridOrdered(&Ghyb, Sigma, Grid->VRe[0], Grid->VRe[Grid->Ni-1]);
  bSigma.N = Sigma->N + 1;
  if (InitAndersonMatrix(&bSigma))
  {
    printf("InitAndersonMatrix failed in ResponseFunctionAddSelfEnergyOnGridAndLPAnd\n");
    return 1;
  }
  bSigma.a[0] = 0;
  for (i = 1; i < G0->N; i++)
  {
    bSigma.b[0] = G0->b[i];
    bSigma.a[1] = G0->a[i + 1] + Sigma->a0;
    for (j = 0; j < Sigma->N; j++)
    {
      bSigma.a[j + 2] =  Sigma->a[j];
      bSigma.b[j + 1] = sqrt(Sigma->b[j]);
    }
    if (AndersonMatrixToListOfPoles(bSigma, &bSigmaPole))
    {
      printf("AndersonMatrixToListOfPoles failed in ResponseFunctionAddSelfEnergyOnGridAndLPAnd\n");
      return 1;
    }
    // add bSigmaPole to the Ghyb
    // both are sorted
    AddToListOfPolesOnGridOrdered(&Ghyb, &bSigmaPole, Grid->VRe[0], Grid->VRe[Grid->Ni-1]);
    FreeListOfPoles(&bSigmaPole);
  }
  G->N = Grid->Ni+1;
  if (InitAndersonMatrix(G))
  {
    printf("InitAndersonMatrix failed in GBathAndListAndRepresentation\n");
    return 1;
  }
  for(i=0;i<Grid->Ni;i++)
  {
    G->a[i+2] = Ghyb.a[i];
    G->b[i+1] = sqrt(Ghyb.b[i]);
  }
  G->a[0] = G0->a[0];
  G->b[0] = G0->b[0];
  G->a[1] = G0->a[1] + Sigma->a[0];
  G->mu = G0->mu;
  strcpy(G->Name, "G");
  FreeListOfPoles(&Ghyb);
  FreeAndersonMatrix(&bSigma);
  OrderAndersonMatrix(G);
  return 0;
}
int ResponseFunctionAddSelfEnergyOnGridBlockAndLPAnd(const BlockAndersonMatrixType *const G0, const BlockListOfPoles *const Sigma, BlockAndersonMatrixType * G, double Epsilon, CompactVectorType *const Grid)
{
  int i,k,l;
  BlockListOfPoles G0iWithSigma,GInv;
  CompactMatrixType G0Bi,GBiSqr,GBi;
  if(G0->N==0)
  {
    CopyBlockAndersonMatrix(*G0, G);
    return 0;
  }
  CompactMatrixType SigmaInv;
  // write the inverse of Sigma as a matrix
  BlockListOfPolesInversedToCompactMatrix(Sigma, &SigmaInv);
  if((G0->QComplex!=0)&&(Sigma->QComplex==0)) MakeCompactRealMatrixComplex(&SigmaInv);
  // we use G0Bi and GBiSqr as a dummie pointer Compact matrix to elements in G0 and G
  G0Bi.Ni=G0->BlockSize;
  G0Bi.Nj=G0->BlockSize;
  G0Bi.QComplex=G0->QComplex;
  G0Bi.Name=SigmaInv.Name;
  // initialize G
  G->BlockSize = G0->BlockSize;
  G->QComplex = SigmaInv.QComplex;
  G->N = 1 + Grid->Ni;
  InitBlockAndersonMatrix(G);
  GBiSqr.Ni=G->BlockSize;
  GBiSqr.Nj=G->BlockSize;
  GBiSqr.QComplex=G->QComplex;
  GBiSqr.Name=G->Name;
  // initialize Ginv
  GInv.N = Grid->Ni;
  GInv.QComplex = G->QComplex;
  GInv.BlockSize = G->BlockSize;
  InitBlockListOfPoles(&GInv);
  if(GInv.QComplex==0)
  {
    for(i=0;i<Grid->Ni;i++)
    {
      GInv.VRe[i*(1+Sqr(GInv.BlockSize))] = Grid->VRe[i];
    }
  }
  else
  {
    for(i=0;i<Grid->Ni;i++)
    {
      GInv.V[i*(1+Sqr(GInv.BlockSize))].Re = Grid->VRe[i];
    }
  }
  // set G.A0 and G.B0
  if(G->QComplex==0)
  {
    memcpy(G->aRe, G0->aRe, 2*Sqr(G0->BlockSize)*sizeof(double));
    memcpy(G->bRe, G0->bRe, Sqr(G0->BlockSize)*sizeof(double));
  }
  else if(G0->QComplex==0)
  {
    for(i=0;i<Sqr(G0->BlockSize);i++)
    {
      G->b[i].Re = G0->bRe[i];
      G->b[i].Im = 0;
    }
    for(i=0;i<2*Sqr(G0->BlockSize);i++)
    {
      G->a[i].Re = G0->aRe[i];
      G->a[i].Im = 0;
    }
  }
  else
  {
    memcpy(G->a, G0->a, 2*Sqr(G0->BlockSize)*sizeof(Complex16));
    memcpy(G->b, G0->b, Sqr(G0->BlockSize)*sizeof(Complex16));
  }
  G->mu=G0->mu;
  // Add sigma to Ginv (we neglect Sigma.A0 in Ginv)
  AddBlockListOfPolesOnMeshOfA(&GInv, Sigma);
  // add Sigma.A0 to G.A1
  if(G->QComplex==0)
  {
    for(i=0;i<Sqr(G0->BlockSize);i++)
    {
      G->aRe[i + Sqr(G0->BlockSize)] += Sigma->a0Re[i];
    }
  }
  else if(Sigma->QComplex==0)
  {
    for(i=0;i<Sqr(G0->BlockSize);i++)
    {
      G->a[i + Sqr(G0->BlockSize)].Re += Sigma->a0Re[i];
    }
  }
  else
  {
    for(i=0;i<Sqr(G0->BlockSize);i++)
    {
      G->a[i + Sqr(G0->BlockSize)].Re += Sigma->a0[i].Re;
      G->a[i + Sqr(G0->BlockSize)].Im += Sigma->a0[i].Im;
    }
  }
  // loop over all blocks in G0 (i=1;i<G0.N;i++ )
  for(i=1;i<G0->N;i++)
  {
    // Set the left top block of SigmaInv equal to G0.A[i+1] + Sigma.A[0]
    if(SigmaInv.QComplex==0)
    {
      for(k=0;k<G0->BlockSize;k++)
      {
        for(l=0;l<G0->BlockSize;l++)
        {
          SigmaInv.VRe[l + k*SigmaInv.Nj] = G0->aRe[(i+1) * Sqr(G0->BlockSize) + l + k*(G0->BlockSize)] + Sigma->a0Re[l + k*(G0->BlockSize)];
        }
      }
    }
    else
    {
      if(Sigma->QComplex==0)
      {
        for(k=0;k<G0->BlockSize;k++)
        {
          for(l=0;l<G0->BlockSize;l++)
          {
            SigmaInv.V[l + k*SigmaInv.Nj].Re = G0->a[(i+1) * Sqr(G0->BlockSize) + l + k*(G0->BlockSize)].Re + Sigma->a0Re[l + k*(G0->BlockSize)];
            SigmaInv.V[l + k*SigmaInv.Nj].Im = G0->a[(i+1) * Sqr(G0->BlockSize) + l + k*(G0->BlockSize)].Im;
          }
        }
      }
      else if(G0->QComplex==0)
      {
        for(k=0;k<G0->BlockSize;k++)
        {
          for(l=0;l<G0->BlockSize;l++)
          {
            SigmaInv.V[l + k*SigmaInv.Nj].Re = G0->aRe[(i+1) * Sqr(G0->BlockSize) + l + k*(G0->BlockSize)] + Sigma->a0[l + k*(G0->BlockSize)].Re;
            SigmaInv.V[l + k*SigmaInv.Nj].Im =                                                               Sigma->a0[l + k*(G0->BlockSize)].Im;
          }
        }
      }
      else
      { // all complex
        for(k=0;k<G0->BlockSize;k++)
        {
          for(l=0;l<G0->BlockSize;l++)
          {
            SigmaInv.V[l + k*SigmaInv.Nj].Re = G0->a[(i+1) * Sqr(G0->BlockSize) + l + k*(G0->BlockSize)].Re + Sigma->a0[l + k*(G0->BlockSize)].Re;
            SigmaInv.V[l + k*SigmaInv.Nj].Im = G0->a[(i+1) * Sqr(G0->BlockSize) + l + k*(G0->BlockSize)].Im + Sigma->a0[l + k*(G0->BlockSize)].Im;
          }
        }
      }
    }
    // SigmaInv to BlockListOfPoles
    HermitianCompactMatrixToBlockListOfPoles(&SigmaInv, &G0iWithSigma, G0->BlockSize);
    // rotate BlockListOpPoles by G0.B[i]
    if(G0->QComplex==0)
    {
      G0Bi.VRe = &(G0->bRe[i*Sqr(G0->BlockSize)]);
    }
    else
    {
      G0Bi.V = &(G0->b[i*Sqr(G0->BlockSize)]);
    }
    BlockListOfPolesTransformInlineNC(&G0iWithSigma, &G0Bi);
    // Add ListofPoles on Grid
    AddBlockListOfPolesOnMeshOfA(&GInv, &G0iWithSigma);
    FreeBlockListOfPoles(&G0iWithSigma);
  }
  FreeCompactMatrix(&SigmaInv);
  // add GInv into G (GInv.A0 = 0 by construction)
  if(G->QComplex==0)
  {
    for(i=0;i<Grid->Ni;i++)
    {
      GBiSqr.VRe = &(GInv.VRe[i * (1+Sqr(G->BlockSize)) + 1]);
      CompactMatrixSqrt(&GBiSqr, &GBi);
      for(k=0;k<G->BlockSize;k++) G->aRe[k*(1+G->BlockSize) + (i+2)*Sqr(G->BlockSize)] = GInv.VRe[i * (1+Sqr(G->BlockSize))];
      memcpy(&(G->bRe[(i+1)*Sqr(G->BlockSize)]), GBi.VRe, Sqr(G->BlockSize)*sizeof(double));
      FreeCompactMatrix(&GBi);
    }
  }
  else
  {
    for(i=0;i<Grid->Ni;i++)
    {
      GBiSqr.V = &(GInv.V[i * (1+Sqr(G->BlockSize)) + 1]);
      CompactMatrixSqrt(&GBiSqr, &GBi);
      for(k=0;k<G->BlockSize;k++) G->a[k*(1+G->BlockSize) + (i+2)*Sqr(G->BlockSize)].Re = GInv.V[i * (1+Sqr(G->BlockSize))].Re;
      memcpy(&(G->b[(i+1)*Sqr(G->BlockSize)]), GBi.V, Sqr(G->BlockSize)*sizeof(double));
      FreeCompactMatrix(&GBi);
    }
  }
  return 0;
}
