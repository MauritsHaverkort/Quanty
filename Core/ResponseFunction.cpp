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
