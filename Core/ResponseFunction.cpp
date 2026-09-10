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
            if( A->VRe[i*(1+Sqr(A->BlockSize))] < *EminPole ) *EminPole=A->VRe[i*(1+Sqr(A->BlockSize))];
            if( A->VRe[i*(1+Sqr(A->BlockSize))] > *EmaxPole ) *EmaxPole=A->VRe[i*(1+Sqr(A->BlockSize))];
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
            if( A->V[i*(1+Sqr(A->BlockSize))].Re < *EminPole ) *EminPole=A->V[i*(1+Sqr(A->BlockSize))].Re;
            if( A->V[i*(1+Sqr(A->BlockSize))].Re > *EmaxPole ) *EmaxPole=A->V[i*(1+Sqr(A->BlockSize))].Re;
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
      case 'd':
      {
        double Emin2, Emax2, diag, Rrow, av, ac;
        ResponsefunctionType R;
        BlockDoubleTriDiagonalMatrixType *A=(BlockDoubleTriDiagonalMatrixType*)(M->data);
        /*  the two chains, each as a block tri-diagonal matrix. Unlike the natural impurity
         *  orbital case above both chains are ordinary baths here, so take the widest range
         *  of the two rather than the valence minimum and the conduction maximum.
         */
        R.type='t';
        R.Block=true;
        R.SuggestedGrid=NULL;
        R.NSuggestedGridPoints=0;
        R.data = (void*)(&(A->val));
        ResponseFunctionEstimateEminAndEmaxPole(&R, EminPole, EmaxPole);
        R.data = (void*)(&(A->con));
        ResponseFunctionEstimateEminAndEmaxPole(&R, &Emin2, &Emax2);
        if(*EminPole > Emin2) *EminPole = Emin2;
        if(*EmaxPole < Emax2) *EmaxPole = Emax2;
        /*  A block tri-diagonal matrix starts at a[1], so the impurity block that this type
         *  puts in front of the two chains is not covered above. Add a Gershgorin bound for
         *  its rows. Summing the magnitudes of the val and the con contribution separately
         *  only widens the bound, which is safe. Each chain is read by its own QComplex, it
         *  need not be the one of A for a mixed object.
         */
        for(j=0;j<A->BlockSize;j++)
        {
          diag = 0;
          Rrow = 0;
          for(k=0;k<A->BlockSize;k++)
          {
            uint32_t idx = j*A->BlockSize + k;
            av = A->val.QComplex ? sqrt(Sqr(A->val.a[idx].Re)+Sqr(A->val.a[idx].Im)) : fabs(A->val.aRe[idx]);
            ac = A->con.QComplex ? sqrt(Sqr(A->con.a[idx].Re)+Sqr(A->con.a[idx].Im)) : fabs(A->con.aRe[idx]);
            if(k==j)
            {
              diag += A->val.QComplex ? A->val.a[idx].Re : A->val.aRe[idx];
              diag += A->con.QComplex ? A->con.a[idx].Re : A->con.aRe[idx];
            }
            else
            {
              Rrow += av + ac;
            }
            if(A->val.N!=0) Rrow += A->val.QComplex ? sqrt(Sqr(A->val.b[idx].Re)+Sqr(A->val.b[idx].Im)) : fabs(A->val.bRe[idx]);
            if(A->con.N!=0) Rrow += A->con.QComplex ? sqrt(Sqr(A->con.b[idx].Re)+Sqr(A->con.b[idx].Im)) : fabs(A->con.bRe[idx]);
          }
          if(diag - Rrow < *EminPole) *EminPole = diag - Rrow;
          if(diag + Rrow > *EmaxPole) *EmaxPole = diag + Rrow;
        }
        if(*EminPole >= *EmaxPole)
        {
          *EminPole -= 1;
          *EmaxPole += 1;
        }
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
      case 'd':
      {
        double eimp, Rrow;
        DoubleTriDiagonalMatrixType *A=(DoubleTriDiagonalMatrixType*)(M->data);
        /*  Gershgorin bound over the rows of the matrix this type represents: one impurity
         *  site with energy val.a[0]+con.a[0], coupled to the two chains by val.b[0] and
         *  con.b[0]. The chains themselves are val.a[1..N] and con.a[1..N]. val.b has
         *  length val.N, so b[0] only exists when the chain is not empty.
         */
        eimp = A->val.a[0] + A->con.a[0];
        Rrow = 0;
        if(A->val.N!=0) Rrow += fabs(A->val.b[0]);
        if(A->con.N!=0) Rrow += fabs(A->con.b[0]);
        *EminPole = eimp - Rrow;
        *EmaxPole = eimp + Rrow;
        for(i=1;i<=A->val.N;i++)
        {
          Rrow = fabs(A->val.b[i-1]);
          if(i<A->val.N) Rrow += fabs(A->val.b[i]);
          if(*EminPole > A->val.a[i] - Rrow) *EminPole = A->val.a[i] - Rrow;
          if(*EmaxPole < A->val.a[i] + Rrow) *EmaxPole = A->val.a[i] + Rrow;
        }
        for(i=1;i<=A->con.N;i++)
        {
          Rrow = fabs(A->con.b[i-1]);
          if(i<A->con.N) Rrow += fabs(A->con.b[i]);
          if(*EminPole > A->con.a[i] - Rrow) *EminPole = A->con.a[i] - Rrow;
          if(*EmaxPole < A->con.a[i] + Rrow) *EmaxPole = A->con.a[i] + Rrow;
        }
        if(*EminPole >= *EmaxPole)
        {
          *EminPole -= 1;
          *EmaxPole += 1;
        }
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
      strncpy(A->Name, ((BADoubleTriDiagonalMatrixType*)M->data)->Name, MAXStrLength);
    }
    break;
    case 'd':
    {
      strncpy(A->Name, ((DoubleTriDiagonalMatrixType*)M->data)->Name, MAXStrLength);
    }
    break;
    default:
    {
      printf("Unrecognized type %c in ResponseFunctionToSpectra\n",M->type);
      return 1;
    }
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
      case 'd':
      {
        A->Name=(char*)malloc(MAXStrLength*sizeof(char));
        A->NSpectra = Sqr((((BlockDoubleTriDiagonalMatrixType*)(M->data))->BlockSize));
        A->Spec = (SpectraType*)malloc( A->NSpectra*sizeof(SpectraType));
        strncpy(A->Name, ((BlockDoubleTriDiagonalMatrixType*)M->data)->Name, MAXStrLength);
        if(BlockDoubleTriDiagonalMatrixToSpectra(GroundStateEnergy,(BlockDoubleTriDiagonalMatrixType*)(M->data), A->Spec, Emin, Emax, NE))
        {
          /*  Nothing inside A->Spec was allocated, so hand back an empty array. Leaving the
           *  half built one would have FreeSpectraTypeArray free spectra that do not exist.
           */
          free(A->Spec);
          free(A->Name);
          A->Spec=NULL;
          A->Name=NULL;
          A->NSpectra=0;
          return 1;
        }
      }
      break;
      default:
      {
        printf("Unrecognized type %c in BlockResponseFunctionToSpectra\n",M->type);
        return 1;
      }
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
      case 'd':
      {
        A->Name=(char*)malloc(MAXStrLength*sizeof(char));
        A->NSpectra = Sqr((((BlockDoubleTriDiagonalMatrixType*)(M->data))->BlockSize));
        A->Spec = (SpectraType*)malloc( A->NSpectra*sizeof(SpectraType));
        strncpy(A->Name, ((BlockDoubleTriDiagonalMatrixType*)M->data)->Name, MAXStrLength);
        BlockDoubleTriDiagonalMatrixToSpectraWithBroadening(GroundStateEnergy,(BlockDoubleTriDiagonalMatrixType*)(M->data), A->Spec, Emin, Emax, NE, GammaL);
      }
      break;
      default:
      {
        printf("Unrecognized type %c in BlockResponseFunctionToSpectra\n",M->type);
        return 1;
      }
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
      case 'd':
      {
        DoubleTriDiagonalMatrixToSpectrum(GroundStateEnergy, *(DoubleTriDiagonalMatrixType*)M->data, A);
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
     case 'd':
     {
       DoubleTriDiagonalMatrixToSpectrum(GroundStateEnergy, *(DoubleTriDiagonalMatrixType*)M->data, A);
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

/*  Given a Green's function in any of the representations a response function can take
 *  the one particle density matrix of the impurity block is returned as a matrix of
 *  size BlockSize by BlockSize (1 by 1 for a scalar response function)
 *    rho = sum_i f(a_i, mu, kbT) B_i
 *  All representations that are not a list of poles are first converted to a list of poles.
 *  Note that this conversion is a diagonalization of the full matrix and thus not cheap.
 *  Note furthermore that the conversion projects onto the first block, the occupation of the
 *  bath states is not part of rho.
 */
/*  The distance below which two poles may be merged when all that is asked of the list is a
 *  density matrix.
 *
 *  rho = sum_i f(a_i,mu,kbT) W_i, and a merge puts the surviving pole at the trace weighted
 *  mean of the two, so the trace of rho -- the occupation -- is right to second order in the
 *  distance delta, while rho itself is right to first order only: the first order term
 *  cancels in the trace and not off it. That error is bounded by |f'|max delta ||W||, which
 *  is delta ||W|| / (4 kbT), so delta = 4 kbT DBL_EPSILON keeps it at the level of the
 *  arithmetic itself. It scales with kbT, which is the point: a merge is only invisible as
 *  far as the Fermi function is flat.
 *
 *  At kbT <= 0 no positive delta qualifies -- f is a step, and merging two poles that
 *  straddle mu moves a whole weight matrix from filled to half filled -- so zero is passed,
 *  which merges nothing at all, the test inside BlockListOfPolesMergeDegeneratePoles being a
 *  strict <.
 *
 *  Merging is never needed for a density matrix. It exists to remove the arbitrary
 *  eigenbasis inside a group of degenerate poles before something else is built out of the
 *  list, and rho does not see that basis: a unitary mixing inside the group redistributes
 *  the individual W_i and leaves their sum, and f is the same for all of them. So the rule
 *  above is a bound on the harm a merge may do and not a requirement.
 */
static double DensityMatrixMinimalPoleDistance(double const kbT)
{
  if(kbT<=0) return 0;
  return 4 * kbT * DBL_EPSILON;
}

uint32_t ResponseFunctionToDensityMatrix(const ResponsefunctionType *M, CompactMatrixType *rho, double mu, double kbT)
{
  if(M->Block)
  {
    BlockListOfPoles G;
    switch (M->type|32)
    {
      case 'l':
      {
        return BlockListOfPolesToDensityMatrix((BlockListOfPoles*)(M->data), rho, mu, kbT);
      }
      case 't':
      {
        if(BlockTriDiagonalMatrixToBlockListOfPoles((BlockTriDiagonalMatrixType*)(M->data), &G, DensityMatrixMinimalPoleDistance(kbT), false))
        {
          printf("BlockTriDiagonalMatrixToBlockListOfPoles failed in ResponseFunctionToDensityMatrix\n");
          return 1;
        }
      }
      break;
      case 'a':
      {
        if(BlockAndersonMatrixToBlockListOfPoles((BlockAndersonMatrixType*)(M->data), &G, DensityMatrixMinimalPoleDistance(kbT), false))
        {
          printf("BlockAndersonMatrixToBlockListOfPoles failed in ResponseFunctionToDensityMatrix\n");
          return 1;
        }
      }
      break;
      case 'n':
      {
        if(BlockNaturalImpurityOrbitalMatrixToBlockListOfPoles((BlockNaturalImpurityOrbitalMatrixType*)(M->data), &G, DensityMatrixMinimalPoleDistance(kbT), false))
        {
          printf("BlockNaturalImpurityOrbitalMatrixToBlockListOfPoles failed in ResponseFunctionToDensityMatrix\n");
          return 1;
        }
      }
      break;
      case 'd':
      {
        if(BlockDoubleTriDiagonalMatrixToBlockListOfPoles((BlockDoubleTriDiagonalMatrixType*)(M->data), &G, DensityMatrixMinimalPoleDistance(kbT), false))
        {
          printf("BlockDoubleTriDiagonalMatrixToBlockListOfPoles failed in ResponseFunctionToDensityMatrix\n");
          return 1;
        }
      }
      break;
      default:
      {
        printf("Error in ResponseFunctionToDensityMatrix: unrecognized type %c\n",M->type);
        return 1;
      }
    }
    BlockListOfPolesToDensityMatrix(&G, rho, mu, kbT);
    FreeBlockListOfPoles(&G);
    return 0;
  }
  else
  {
    ListOfPoles L;
    switch (M->type|32)
    {
      case 'l':
      {
        return ListOfPolesToDensityMatrix(*(ListOfPoles*)(M->data), rho, mu, kbT);
      }
      case 't':
      {
        if(TriDiagonalMatrixToListOfPoles(*(TriDiagonalMatrixType*)(M->data), &L))
        {
          printf("TriDiagonalMatrixToListOfPoles failed in ResponseFunctionToDensityMatrix\n");
          return 1;
        }
      }
      break;
      case 'a':
      {
        if(AndersonMatrixToListOfPoles(*(AndersonMatrixType*)(M->data), &L))
        {
          printf("AndersonMatrixToListOfPoles failed in ResponseFunctionToDensityMatrix\n");
          return 1;
        }
      }
      break;
      case 'n':
      {
        if(BADoubleTriDiagonalMatrixToListOfPoles((BADoubleTriDiagonalMatrixType*)(M->data), &L))
        {
          printf("BADoubleTriDiagonalMatrixToListOfPoles failed in ResponseFunctionToDensityMatrix\n");
          return 1;
        }
      }
      break;
      case 'd':
      {
        /*  no merge distance to pass: the single valued conversions take none, a list of
         *  poles of a scalar G having no weight matrices whose basis could be arbitrary */
        if(DoubleTriDiagonalMatrixToListOfPoles(*(DoubleTriDiagonalMatrixType*)(M->data), &L))
        {
          printf("DoubleTriDiagonalMatrixToListOfPoles failed in ResponseFunctionToDensityMatrix\n");
          return 1;
        }
      }
      break;
      default:
      {
        printf("Error in ResponseFunctionToDensityMatrix: unrecognized type %c\n",M->type);
        return 1;
      }
    }
    ListOfPolesToDensityMatrix(L, rho, mu, kbT);
    FreeListOfPoles(&L);
    return 0;
  }
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
    case 'D':
    case 'd':
    {
      return ((BlockDoubleTriDiagonalMatrixType*)a->data)->BlockSize;
      break;
    }
    default:
    {
      printf("Error in ResponseFunctionGetBlockSize: unrecognized type %c\n",a->type);
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
    case 'D':
    case 'd':
    {
      return ((BlockDoubleTriDiagonalMatrixType*)a->data)->QComplex;
      break;
    }
    default:
    {
      printf("Error in ResponseFunctionGetQComplex: unrecognized type %c\n",a->type);
      return 0;
    }
  }
}

/*  returns the chemical potential stored in the response function
 */
double ResponseFunctionGetMu(const ResponsefunctionType *const a)
{
  switch(a->type)
  {
    case 'L':
    case 'l':
    {
      if(a->Block) return ((BlockListOfPoles*)a->data)->mu;
      return ((ListOfPoles*)a->data)->mu;
    }
    case 'A':
    case 'a':
    {
      if(a->Block) return ((BlockAndersonMatrixType*)a->data)->mu;
      return ((AndersonMatrixType*)a->data)->mu;
    }
    case 'T':
    case 't':
    {
      if(a->Block) return ((BlockTriDiagonalMatrixType*)a->data)->mu;
      return ((TriDiagonalMatrixType*)a->data)->mu;
    }
    case 'N':
    case 'n':
    {
      if(a->Block) return ((BlockNaturalImpurityOrbitalMatrixType*)a->data)->mu;
      return ((BADoubleTriDiagonalMatrixType*)a->data)->mu;
    }
    case 'D':
    case 'd':
    {
      if(a->Block) return ((BlockDoubleTriDiagonalMatrixType*)a->data)->mu;
      return ((DoubleTriDiagonalMatrixType*)a->data)->mu;
    }
    default:
    {
      printf("Error in ResponseFunctionGetMu: unrecognized type %c\n",a->type);
      return 0;
    }
  }
}

uint32_t ResponseFunctionGetA0(const ResponsefunctionType *const a, CompactMatrixType *M)
{
  M->Ni = ResponseFunctionGetBlockSize(a);
  M->Nj = M->Ni;
  /*  A CompactMatrixType is not a compound object and has no mixed state, so the 0/1/2
   *  QComplex of the response function is narrowed to 0 or 1 here.
   */
  M->QComplex = ResponseFunctionGetQComplex(a) ? 1 : 0;
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
      case 'D':
      case 'd':
      {
        M->VRe[0] = ((DoubleTriDiagonalMatrixType*)a->data)->a0;
        break;
      }
      default:
      {
        printf("Error in ResponseFunctionGetA0: unrecognized type %c\n",a->type);
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
        /*  a0 carries its own QComplex, which for a mixed response function (QComplex 2)
         *  need not be the one of the response function, so read a0 by its own flag and
         *  write M by M's. a0 complex forces M complex, so nothing is ever dropped.
         */
        CompactMatrixType const * const a0 = &(((BlockNaturalImpurityOrbitalMatrixType*)a->data)->a0);
        for(int i=0; i<Sqr(M->Ni);i++)
        {
          double re = a0->QComplex ? a0->V[i].Re : a0->VRe[i];
          double im = a0->QComplex ? a0->V[i].Im : 0;
          if(M->QComplex==0)
          {
            M->VRe[i] = re;
          }
          else
          {
            M->V[i].Re = re;
            M->V[i].Im = im;
          }
        }
        break;
      }
      case 'D':
      case 'd':
      {
        /*  a0 carries its own QComplex, which for a mixed response function (QComplex 2)
         *  need not be the one of the response function, so read a0 by its own flag and
         *  write M by M's. a0 complex forces M complex, so nothing is ever dropped.
         */
        CompactMatrixType const * const a0 = &(((BlockDoubleTriDiagonalMatrixType*)a->data)->a0);
        for(int i=0; i<Sqr(M->Ni);i++)
        {
          double re = a0->QComplex ? a0->V[i].Re : a0->VRe[i];
          double im = a0->QComplex ? a0->V[i].Im : 0;
          if(M->QComplex==0)
          {
            M->VRe[i] = re;
          }
          else
          {
            M->V[i].Re = re;
            M->V[i].Im = im;
          }
        }
        break;
      }
      default:
      {
        printf("Error in ResponseFunctionGetA0: unrecognized type %c\n",a->type);
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
          BlockAndersonMatrixToBlockTridiagonalMatrix((BlockAndersonMatrixType*)G0->data, G0Data, Epsilon, -1, false);
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
          BlockAndersonMatrixToBlockTridiagonalMatrix((BlockAndersonMatrixType*)Sigma->data, SigmaData, Epsilon, -1, false);
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
          BlockAndersonMatrixToBlockTridiagonalMatrix((BlockAndersonMatrixType*)G0->data, G0Data, Epsilon, -1, false);
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
        GMat.QComplex = G0Data->QComplex || SData->QComplex;
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
          for(int k=0;k<SData->N;k++)
          {
            for(int i=0;i<SData->BlockSize;i++)
            {
              for(int j=0;j<SData->BlockSize;j++)
              {
                GMat.V[j+i*GMat.Nj + k*(GMat.Nj+1)*(SData->BlockSize)                           ].Re = SData->aRe[j+i*SData->BlockSize + k*Sqr(SData->BlockSize)];
                GMat.V[j+i*GMat.Nj + k*(GMat.Nj+1)*(SData->BlockSize)                           ].Im = 0;
                GMat.V[j+i*GMat.Nj + k*(GMat.Nj+1)*(SData->BlockSize) + SData->BlockSize        ].Re = SData->bRe[j+i*SData->BlockSize + k*Sqr(SData->BlockSize)];
                GMat.V[j+i*GMat.Nj + k*(GMat.Nj+1)*(SData->BlockSize) + SData->BlockSize        ].Im = 0;
                GMat.V[j+i*GMat.Nj + k*(GMat.Nj+1)*(SData->BlockSize) + SData->BlockSize*GMat.Nj].Re = SData->bRe[i+j*SData->BlockSize + k*Sqr(SData->BlockSize)];
                GMat.V[j+i*GMat.Nj + k*(GMat.Nj+1)*(SData->BlockSize) + SData->BlockSize*GMat.Nj].Im = 0;
              }
            }
          }
          {
            int k=SData->N;
            for(int i=0;i<SData->BlockSize;i++)
            {
              for(int j=0;j<SData->BlockSize;j++)
              {
                GMat.V[j+i*GMat.Nj + k*(GMat.Nj+1)*(SData->BlockSize)                           ].Re = SData->aRe[j+i*SData->BlockSize + k*Sqr(SData->BlockSize)];
                GMat.V[j+i*GMat.Nj + k*(GMat.Nj+1)*(SData->BlockSize)                           ].Im = 0;
              }
            }
          }
        }
        else
        { // GMat and SData both complex
          for(int k=0;k<SData->N;k++)
          {
            for(int i=0;i<SData->BlockSize;i++)
            {
              for(int j=0;j<SData->BlockSize;j++)
              {
                GMat.V[j+i*GMat.Nj + k*(GMat.Nj+1)*(SData->BlockSize)                           ].Re = SData->a[j+i*SData->BlockSize + k*Sqr(SData->BlockSize)].Re;
                GMat.V[j+i*GMat.Nj + k*(GMat.Nj+1)*(SData->BlockSize)                           ].Im = SData->a[j+i*SData->BlockSize + k*Sqr(SData->BlockSize)].Im;
                GMat.V[j+i*GMat.Nj + k*(GMat.Nj+1)*(SData->BlockSize) + SData->BlockSize        ].Re = SData->b[j+i*SData->BlockSize + k*Sqr(SData->BlockSize)].Re;
                GMat.V[j+i*GMat.Nj + k*(GMat.Nj+1)*(SData->BlockSize) + SData->BlockSize        ].Im = SData->b[j+i*SData->BlockSize + k*Sqr(SData->BlockSize)].Im;
                GMat.V[j+i*GMat.Nj + k*(GMat.Nj+1)*(SData->BlockSize) + SData->BlockSize*GMat.Nj].Re = SData->b[i+j*SData->BlockSize + k*Sqr(SData->BlockSize)].Re;
                GMat.V[j+i*GMat.Nj + k*(GMat.Nj+1)*(SData->BlockSize) + SData->BlockSize*GMat.Nj].Im =-SData->b[i+j*SData->BlockSize + k*Sqr(SData->BlockSize)].Im;
              }
            }
          }
          {
            int k=SData->N;
            for(int i=0;i<SData->BlockSize;i++)
            {
              for(int j=0;j<SData->BlockSize;j++)
              {
                GMat.V[j+i*GMat.Nj + k*(GMat.Nj+1)*(SData->BlockSize)                           ].Re = SData->a[j+i*SData->BlockSize + k*Sqr(SData->BlockSize)].Re;
                GMat.V[j+i*GMat.Nj + k*(GMat.Nj+1)*(SData->BlockSize)                           ].Im = SData->a[j+i*SData->BlockSize + k*Sqr(SData->BlockSize)].Im;
              }
            }
          }
        }
      }
      break;
      case 'a':
      {
        // transform Sigma to tridiagonal type and call the function again.
        // this part can be replaced by filling the matrix with the anderson form without the aditional transformation
        ResponsefunctionType SigmaTri;
        ConvertResponsefunction(Sigma, &SigmaTri, 't', "Default", 0, Epsilon, 0, 0, false, false, false);
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
      case 'n':
      {
        // transform Sigma to tridiagonal type and call the function again.
        // this part can be replaced by filling the matrix with the natural orbital form without the aditional transformation
        ResponsefunctionType SigmaTri;
        ConvertResponsefunction(Sigma, &SigmaTri, 't', "Default", 0, Epsilon, 0, 0, false, false, false);
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
      case 'l':
      {
        // transform Sigma to tridiagonal type and call the function again.
        ResponsefunctionType SigmaTri;
        ConvertResponsefunction(Sigma, &SigmaTri, 't', "Default", 0, Epsilon, 0, 0, false, false, false);
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
      default:
      {
        printf("ERORR in ResponseFunctionAddSelfEnergyOnGridViaG0Tri\n Unrecognized type for Sigma\n");
        return 1;
      }
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
    printf("ERORR in ResponseFunctionAddSelfEnergyOnGridViaG0Tri\n Method only implemented for Block response functions\n");
    return 1;
  }
  return 0;
}
int ResponseFunctionAddSelfEnergyOnGridTriCompactLP(const TriDiagonalMatrixType *const G0, CompactMatrixType * Sigma, uint32_t SigmaDim, CompactVectorType *const Grid, ListOfPoles * G)
{
  printf("TODO\n");
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
          BlockListOfPolesToBlockAndersonMatrix(*(BlockListOfPoles*)G0->data, G0Data, 1+((BlockListOfPoles*)G0->data)->N, 1000*DBL_EPSILON);
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
          BlockAndersonMatrixToBlockListOfPoles((BlockAndersonMatrixType*)Sigma->data, SigmaData, 1000*DBL_EPSILON,false);
        }
        break;
        case 't':
        {
          BlockTriDiagonalMatrixToBlockListOfPoles((BlockTriDiagonalMatrixType*)Sigma->data, SigmaData, 1000*DBL_EPSILON,false);
        }
        break;
        case 'n':
        {
          BlockNaturalImpurityOrbitalMatrixToBlockListOfPoles((BlockNaturalImpurityOrbitalMatrixType*)Sigma->data, SigmaData, -1, false);
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
    BlockListOfPolesMergeDegeneratePoles(&G0iWithSigma, Epsilon, false);
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

/* ===========================================================================
 *
 *   G(w) for a Hamiltonian H in a space P plus a self energy given as a set of
 *   operators V and block response functions Sigma
 *
 *   The mathematics is written out in front of CalculateGforHwithSigmaAsOperator
 *   below.  The helpers here are shared with nothing else yet, but they are
 *   deliberately not static and carry generic names: should the same thing be
 *   implemented a second time somewhere in the amalgamated build the linker
 *   complains instead of the two copies silently drifting apart.
 *
 * =========================================================================== */

/*  relative cut on the eigenvalues of the overlap matrix of PsiBasis            */
static double const SigmaAsOperatorEpsOrtho = 1e-10;
/*  relative cut on the eigenvalues of the extended overlap matrix Sext          */
static double const SigmaAsOperatorEpsDefl = 1e-10;
/*  relative cut when a grid residue or a pole residue is factorized             */
static double const SigmaAsOperatorEpsChol = 1e-10;
/*  relative tolerance on negative spectral weight in a residue of Sigma         */
static double const SigmaAsOperatorEpsPSD = 1e-8;
/*  relative cut in the re-orthogonalization of the outer block Lanczos          */
static double const SigmaAsOperatorEpsLanczos = 1e-12;
/*  above this many bath states the poles are binned on a grid instead of kept   */
static uint32_t const SigmaAsOperatorMaxBathExact = 20000;
/*  memory (in Complex16) that the grid accumulators are allowed to take         */
static double const SigmaAsOperatorMaxGridMemory = 1.25e8;
/*  hard cap on the dimension of the extended overlap matrix Sext                */
static uint32_t const SigmaAsOperatorMaxLi = 20000;

/*  Marks a compact matrix as holding nothing, so that a call to
 *  FreeCompactMatrix on an error path before it was ever initialized is a no-op
 */
static void CompactMatrixSetUnallocated(CompactMatrixType *const M)
{
  M->Ni = 0;
  M->Nj = 0;
  M->QComplex = 1;
  M->V = NULL;
  M->Name = NULL;
  M->a0 = 0;
  M->b0 = 1;
  M->mu = 0;
  M->OnLuaStack = false;
}

/*  Eigen decomposition of a Hermitian matrix with the convention
 *
 *    A_{jk} = sum_n val[n] fun[n*N+j] conj( fun[n*N+k] )
 *
 *  i.e. row n of fun holds eigenvector n and the eigenvalues come out in
 *  ascending order.  This is what CompactMatrixEigensystem already returns: it
 *  hands the row major array to LAPACK as if it were column major, so for a
 *  complex matrix LAPACK diagonalizes the transpose, and the conjugation at the
 *  end of that routine undoes it.  The wrapper exists to state the convention
 *  in one place: getting it wrong is silent for a real symmetric matrix and
 *  transposes the answer for a complex Hermitian one.
 */
int HermitianCompactMatrixEigensystemRowVectors(CompactMatrixType const *const A, CompactMatrixType *const fun, double *const val)
{
  if (CompactMatrixEigensystem(A, fun, val))
  {
    printf("CompactMatrixEigensystem failed in HermitianCompactMatrixEigensystemRowVectors\n");
    return 1;
  }
  return 0;
}

/*  Out_{ab} = < PsiBasis[a] | PsiBasis[b] >
 *  Out is allocated NPsi x NPsi complex on entry
 */
int PsiListPsiListCompactMatrix(WaveFunctionType const *const PsiBasis, uint32_t const NPsi, CompactMatrixType *const Out)
{
  uint32_t a, b;
  double Re, Im;
  for (b = 0; b < NPsi; b++)
  {
    for (a = 0; a <= b; a++)
    {
      WaveFunctionDotProduct(PsiBasis[a], PsiBasis[b], &Re, &Im);
      Out->V[a * NPsi + b].Re = Re;
      Out->V[a * NPsi + b].Im = Im;
      Out->V[b * NPsi + a].Re = Re;
      Out->V[b * NPsi + a].Im = -Im;
    }
  }
  return 0;
}

/*  Out_{ab} = < PsiBasis[a] | O | PsiBasis[b] >
 *  Out is allocated NPsi x NPsi complex on entry
 *  O is applied NPsi times, one scratch wave function is resident
 */
int PsiListOperatorPsiListCompactMatrix(OperatorType const *const O, WaveFunctionType const *const PsiBasis, uint32_t const NPsi, CompactMatrixType *const Out)
{
  uint32_t a, b;
  double Re, Im;
  WaveFunctionType Opsi;
  for (b = 0; b < NPsi; b++)
  {
    if (OperatorPsi(*O, PsiBasis[b], &Opsi))
    {
      printf("OperatorPsi failed in PsiListOperatorPsiListCompactMatrix\n");
      return 1;
    }
    for (a = 0; a < NPsi; a++)
    {
      WaveFunctionDotProduct(PsiBasis[a], Opsi, &Re, &Im);
      Out->V[a * NPsi + b].Re = Re;
      Out->V[a * NPsi + b].Im = Im;
    }
    WaveFunctionFree(&Opsi);
  }
  return 0;
}

/*  A view on rows [row0,row0+nrow) of A.  The view shares memory with A and must
 *  never be freed, which OnLuaStack takes care of.
 */
CompactMatrixType CompactMatrixRowBlockView(CompactMatrixType const *const A, uint32_t const row0, uint32_t const nrow)
{
  CompactMatrixType v;
  v = *A;
  v.Ni = nrow;
  v.V = &(A->V[(size_t)row0 * A->Nj]);
  v.OnLuaStack = true;
  return v;
}

/*  Orthonormalizes the rows of the pair (R1,R2) -- one vector is the row of R1
 *  followed by the row of R2 -- by modified Gram Schmidt, done twice for
 *  stability, and returns the transformation L, allocated Ni x Ni on entry:
 *
 *    R_old[q] = sum_p L[q*Ni+p] R_new[p]
 *
 *  Rows that are linearly dependent on the earlier ones are set to zero and get
 *  a zero row in L.  Zero rows keep the block size uniform without adding
 *  anything to the matrix that is being built, which is what the block Lanczos
 *  below needs when the Krylov space runs out before order blocks are filled.
 *  R2 may have zero columns.
 */
int RowBlockPairOrthonormalizeWithTransform(CompactMatrixType *const R1, CompactMatrixType *const R2, CompactMatrixType *const L, double const zero)
{
  uint32_t p, q, r, pass, n1 = R1->Nj, n2 = R2->Nj;
  double cRe, cIm, norm;
  SetCompactMatrixToZero(L);
  for (q = 0; q < R1->Ni; q++)
  {
    for (pass = 0; pass < 2; pass++)
    {
      for (p = 0; p < q; p++)
      {
        cRe = 0;
        cIm = 0;
        for (r = 0; r < n1; r++)
        {
          cRe += R1->V[(size_t)p * n1 + r].Re * R1->V[(size_t)q * n1 + r].Re + R1->V[(size_t)p * n1 + r].Im * R1->V[(size_t)q * n1 + r].Im;
          cIm += R1->V[(size_t)p * n1 + r].Re * R1->V[(size_t)q * n1 + r].Im - R1->V[(size_t)p * n1 + r].Im * R1->V[(size_t)q * n1 + r].Re;
        }
        for (r = 0; r < n2; r++)
        {
          cRe += R2->V[(size_t)p * n2 + r].Re * R2->V[(size_t)q * n2 + r].Re + R2->V[(size_t)p * n2 + r].Im * R2->V[(size_t)q * n2 + r].Im;
          cIm += R2->V[(size_t)p * n2 + r].Re * R2->V[(size_t)q * n2 + r].Im - R2->V[(size_t)p * n2 + r].Im * R2->V[(size_t)q * n2 + r].Re;
        }
        L->V[q * L->Nj + p].Re += cRe;
        L->V[q * L->Nj + p].Im += cIm;
        for (r = 0; r < n1; r++)
        {
          R1->V[(size_t)q * n1 + r].Re -= cRe * R1->V[(size_t)p * n1 + r].Re - cIm * R1->V[(size_t)p * n1 + r].Im;
          R1->V[(size_t)q * n1 + r].Im -= cRe * R1->V[(size_t)p * n1 + r].Im + cIm * R1->V[(size_t)p * n1 + r].Re;
        }
        for (r = 0; r < n2; r++)
        {
          R2->V[(size_t)q * n2 + r].Re -= cRe * R2->V[(size_t)p * n2 + r].Re - cIm * R2->V[(size_t)p * n2 + r].Im;
          R2->V[(size_t)q * n2 + r].Im -= cRe * R2->V[(size_t)p * n2 + r].Im + cIm * R2->V[(size_t)p * n2 + r].Re;
        }
      }
    }
    norm = 0;
    for (r = 0; r < n1; r++) norm += Sqr(R1->V[(size_t)q * n1 + r].Re) + Sqr(R1->V[(size_t)q * n1 + r].Im);
    for (r = 0; r < n2; r++) norm += Sqr(R2->V[(size_t)q * n2 + r].Re) + Sqr(R2->V[(size_t)q * n2 + r].Im);
    norm = sqrt(norm);
    if (norm <= zero)
    {
      for (r = 0; r < n1; r++) { R1->V[(size_t)q * n1 + r].Re = 0; R1->V[(size_t)q * n1 + r].Im = 0; }
      for (r = 0; r < n2; r++) { R2->V[(size_t)q * n2 + r].Re = 0; R2->V[(size_t)q * n2 + r].Im = 0; }
      L->V[q * L->Nj + q].Re = 0;
      L->V[q * L->Nj + q].Im = 0;
    }
    else
    {
      L->V[q * L->Nj + q].Re = norm;
      L->V[q * L->Nj + q].Im = 0;
      for (r = 0; r < n1; r++) { R1->V[(size_t)q * n1 + r].Re /= norm; R1->V[(size_t)q * n1 + r].Im /= norm; }
      for (r = 0; r < n2; r++) { R2->V[(size_t)q * n2 + r].Re /= norm; R2->V[(size_t)q * n2 + r].Im /= norm; }
    }
  }
  return 0;
}

/*  The residue matrix and the energy of pole n of a block list of poles.
 *  W is allocated BlockSize x BlockSize complex on entry.
 */
int BlockListOfPolesGetPole(BlockListOfPoles const *const A, uint32_t const n, CompactMatrixType *const W, double *const energy)
{
  uint32_t k, bs = A->BlockSize, stride = 1 + bs * bs;
  if (A->QComplex)
  {
    *energy = A->V[n * stride].Re;
    for (k = 0; k < bs * bs; k++)
    {
      W->V[k].Re = A->V[n * stride + 1 + k].Re;
      W->V[k].Im = A->V[n * stride + 1 + k].Im;
    }
  }
  else
  {
    *energy = A->VRe[n * stride];
    for (k = 0; k < bs * bs; k++)
    {
      W->V[k].Re = A->VRe[n * stride + 1 + k];
      W->V[k].Im = 0;
    }
  }
  SquareCompactMatrixMakeHermitian(W);
  return 0;
}

/*  The constant term a0 of a block list of poles.
 *  A0 is allocated BlockSize x BlockSize complex on entry.
 */
int BlockListOfPolesGetA0(BlockListOfPoles const *const A, CompactMatrixType *const A0)
{
  uint32_t k;
  for (k = 0; k < A->BlockSize * A->BlockSize; k++)
  {
    if (A->QComplex)
    {
      A0->V[k].Re = A->a0[k].Re;
      A0->V[k].Im = A->a0[k].Im;
    }
    else
    {
      A0->V[k].Re = A->a0Re[k];
      A0->V[k].Im = 0;
    }
  }
  return 0;
}

/*  Factorizes a Hermitian positive semi definite matrix R (N x N) as
 *
 *    R = Y^dag Y ,   Y of size rho x N
 *
 *  keeping the eigenvalues above SingularValue * lambda_max.  Y is initialized
 *  here and has to be freed by the caller.  The discarded and the negative
 *  weight are added to DiscardedWeight and NegativeWeight, which lets the caller
 *  report over the whole run instead of once per matrix.
 */
int HermitianCompactMatrixToRectangularRoot(CompactMatrixType const *const R, CompactMatrixType *const Y, double const SingularValue, double *const DiscardedWeight, double *const NegativeWeight)
{
  uint32_t k, p, rho, N = R->Ni;
  double *val, lmax, s;
  CompactMatrixType fun;
  val = (double*)malloc(N * sizeof(double));
  if (val == NULL)
  {
    printf("malloc failed in HermitianCompactMatrixToRectangularRoot\n");
    return 1;
  }
  if (HermitianCompactMatrixEigensystemRowVectors(R, &fun, val))
  {
    free(val);
    return 1;
  }
  lmax = Max(val[N - 1], 0.0);
  rho = 0;
  for (k = 0; k < N; k++)
  {
    if (val[k] > SingularValue * lmax) rho++;
    else if (val[k] < -SingularValue * lmax) *NegativeWeight -= val[k];
    else if (val[k] > 0) *DiscardedWeight += val[k];
  }
  Y->Ni = Max(rho, 1);
  Y->Nj = N;
  Y->QComplex = 1;
  if (InitCompactMatrix(Y))
  {
    printf("InitCompactMatrix failed in HermitianCompactMatrixToRectangularRoot\n");
    FreeCompactMatrix(&fun);
    free(val);
    return 1;
  }
  Y->Ni = rho;
  rho = 0;
  for (k = 0; k < N; k++)
  {
    if (val[k] <= SingularValue * lmax) continue;
    s = sqrt(val[k]);
    for (p = 0; p < N; p++)
    {
      Y->V[rho * N + p].Re = s * fun.V[k * N + p].Re;
      Y->V[rho * N + p].Im = -s * fun.V[k * N + p].Im;
    }
    rho++;
  }
  FreeCompactMatrix(&fun);
  free(val);
  return 0;
}


/*
 *  Given a one particle basis | phi_{\tau} > with \tau in [0,Ntau] ( Ntau = H->NFermion = PsiBasis[i].NFermion )
 *  Given an operator H and basis PsiBasis[i] with i in [0,NPsiBasis)
 *  Given a self energy defined by the operator
 *    SigmaOpp = sum_{i=0}^{NSigma-1} SigmaOpp[i]
 *    SigmaOpp[i] = sum_{j,k} SigmaOpp[i]_{jk}
 *    SigmaOpp[i]_{jk} = V[i]_j^{\dag} Sigma[i]_{j,k}(omega) V[i]_k
 *      with j,k in [0,Sigma[i].BlockSize)
 *  calculate
 *    G_{lm}(w) = < PsiBasis[l] | 1 / ( w - H - SigmaOpp(w) ) | PsiBasis[m] >
 *    with l,m in [0,BlockSize)
 *  and return it in Block Anderson (arrowhead) form in G.
 *
 *  --------------------------------------------------------------------------
 *  What is calculated
 *  --------------------------------------------------------------------------
 *
 *  Sigma[i] is a block list of poles,
 *
 *    Sigma[i]_{jj'}(w) = a0[i]_{jj'} + sum_n W[i][n]_{jj'} / ( w - a[i][n] )
 *
 *  with W[i][n] Hermitian and positive semi definite.  Reading the poles as
 *  one particle states outside P, at energy a[i][n] and hopping to the channel
 *  operators V[i][j], the total Hamiltonian on P plus one excitation into the
 *  bath (the space P + Q1) is H plus those bath energies plus the hopping.  H
 *  itself contains neither, and must not.
 *
 *  Eliminating Q1 exactly (Feshbach) gives, without any approximation beyond
 *  the truncation at one bath excitation,
 *
 *    Sigma_{ab}(w) = sum_i sum_n sum_{jj'} W[i][n]_{jj'}
 *                    < a | V[i][j]^dag ( w - a[i][n] - H )^{-1} V[i][j'] | b >
 *
 *  The resolvent of H sits between the two channel operators, so every bath
 *  energy a[i][n] is split by the whole spectrum of H in the space reached by
 *  V[i][j] | PsiBasis[b] >.  Dropping that splitting would give the naive
 *  hybridization-only self energy, which is wrong.
 *
 *  H commutes with the bath operators, so the same intermediate space serves
 *  all poles of a block: it has to be built once per i, not once per pole.
 *
 *  --------------------------------------------------------------------------
 *  How it is calculated
 *  --------------------------------------------------------------------------
 *
 *  B  PsiBasis is a Krylov basis and is only locally orthogonal, so the overlap
 *     matrix SMat_{ab} = <a|b> is formed, Loewdin orthogonalized (dropping the
 *     null directions that the loss of global orthogonality produces) and
 *     everything is transformed into that orthonormal basis of dimension
 *     North <= NPsiBasis.  HMat_{ab} = <a|H|b> is formed with NPsiBasis
 *     applications of H.
 *
 *  C  Per i, with mu = (j,a) a composite index of length Li = bi * North,
 *
 *       Sext_{mu,nu} = < a | V[i][j]^dag       V[i][j'] | b >
 *       Hext_{mu,nu} = < a | V[i][j]^dag   H   V[i][j'] | b >
 *
 *     Both are matrix elements of operators that act inside P, so they are
 *     built from operator products and cost no wave function outside P.  Sext
 *     is diagonalized and deflated, Hext is diagonalized in what is left.
 *
 *     This is the one approximation in the routine.  It is Rayleigh-Ritz in
 *     the space spanned by V[i][j]|PsiBasis[a]>, and that span is not in
 *     general closed under H -- closing it would be a block Lanczos in the
 *     N-1 (or N+1) sector, which needs the whole deflated seed block resident
 *     as wave functions and is what makes the exact treatment expensive.  So:
 *
 *       -  whenever the span happens to be invariant under H, and in
 *          particular whenever the channels reach the whole sector, the result
 *          is exact within P + Q1 and agrees with an explicit bath to machine
 *          precision;
 *       -  otherwise the first two moments of the intermediate resolvent are
 *          still exact, which is what matters because Sigma is only ever
 *          evaluated inside a resolvent at w + i Gamma.  On a deliberately
 *          bad case (three modes in P of which only one is a channel, Gamma
 *          a tenth of the width of the intermediate spectrum) this leaves a
 *          10 percent error on a spectrum for which leaving the self energy
 *          out altogether is 70 percent wrong.
 *
 *  D  Z[i]_{n',(j,a)} = < n' | V[i][j] | PsiBasis[a] >, which follows from the
 *     two diagonalizations without touching a wave function again.  The zeroth
 *     moment sum rule Z^dag Z == Sext is checked.
 *
 *  E  The bath of the arrowhead.  Each triple (i,n,n') contributes a pole at
 *     a[i][n] + E[i][n'] with residue G^dag W[i][n] G, G = Z[i][n'] read as a
 *     bi x North matrix.  Below SigmaAsOperatorMaxBathExact bath states these
 *     are kept one by one and the result is exact.  Above it, or when the
 *     caller supplies a grid, the poles are binned on the grid conserving the
 *     zeroth and first moment of each bin, which caps the bath at Ng * North.
 *
 *  F  Block Lanczos on the arrowhead from the BlockSize seed vectors, order
 *     blocks, full re-orthogonalization, the arrowhead structure used in the
 *     matrix vector product.  This keeps 2*order moments of the BlockSize x
 *     BlockSize matrix spectral function exact.
 *
 *  G  Chain to star and the result padded up to the block size of B0.
 *
 *  --------------------------------------------------------------------------
 *  Arguments
 *  --------------------------------------------------------------------------
 *
 *  H                    the Hamiltonian in P, without the bath and without the
 *                       hybridization that Sigma stands for
 *  PsiBasis, NPsiBasis  a basis of the part of P that is kept.  The first
 *                       BlockSize entries are the seeds, the rest expands them
 *  BlockSize            number of seed vectors, i.e. BlockSizeN[0] of the block
 *                       Lanczos that produced PsiBasis, after its deflation
 *  V                    V[i][j] with j < Sigma[i].BlockSize, the channels
 *  Sigma, NSigma        the self energies.  NSigma == 0 is legal and gives the
 *                       plain resolvent of H
 *  order                number of blocks of the returned Anderson matrix
 *  SuggestedGrid        grid to bin the poles of Sigma on, NULL for automatic.
 *                       Binning is exact in the zeroth and first moment of a
 *                       bin and errs at O(dg^2) in the second, so a spacing of
 *                       about the intended broadening is enough.  Dense over
 *                       the spectral window and sparse in the tails is right
 *  B0                   BlockSizeFull x BlockSize, the transformation from the
 *                       deflated seed block back to the full set of transition
 *                       operators (b[0] of the block tridiagonal matrix that
 *                       produced PsiBasis).  Its number of rows sets the block
 *                       size of G.  NULL means the identity
 *  G                    the result, not initialized on entry
 *
 *  returns 0 on success, 1 on failure
 */
int CalculateGforHwithSigmaAsOperator(const OperatorType *const H, const WaveFunctionType *const PsiBasis, const uint32_t NPsiBasis, const uint32_t BlockSize,  OperatorType ** V, const BlockListOfPoles *const Sigma, const uint32_t NSigma, const uint32_t order, const double *const SuggestedGrid, const uint32_t NSuggestedGridPoints, const CompactMatrixType *const B0, BlockAndersonMatrixType * G)
{
  uint32_t i, j, jp, k, n, np, p, q, a;
  uint32_t Npsi = NPsiBasis, m = BlockSize, NFull, North = 0;
  uint32_t Dtot, Nbath = 0, NbathMax = 0, orderUse, Ng = 0;
  uint32_t QComplexOut;
  int ierr = 0;
  double *Sval = NULL, *Grid = NULL, *Ebath = NULL;
  double **Elist = NULL;
  double norm, Re, Emin, Emax;
  double DiscardedWeight = 0, NegativeWeight = 0;
  CompactMatrixType SMat, HMat, Sfun, T, TS, HO, Cseed, Ybath, Ablk, Bblk, Lblk;
  CompactMatrixType *Zlist = NULL, *Rgrid = NULL;
  uint32_t *Dlist = NULL;
  BlockTriDiagonalMatrixType chain;
  BlockAndersonMatrixType GSmall;
  OperatorType Heff, *Hptr;
  bool QHeffAllocated = false;

  CompactMatrixSetUnallocated(&SMat);  CompactMatrixSetUnallocated(&HMat);
  CompactMatrixSetUnallocated(&Sfun);  CompactMatrixSetUnallocated(&T);
  CompactMatrixSetUnallocated(&TS);    CompactMatrixSetUnallocated(&HO);
  CompactMatrixSetUnallocated(&Cseed); CompactMatrixSetUnallocated(&Ybath);
  CompactMatrixSetUnallocated(&Ablk);  CompactMatrixSetUnallocated(&Bblk);
  CompactMatrixSetUnallocated(&Lblk);
  chain.a = NULL;
  GSmall.a = NULL;

  /* ---------------------------------------------------------------------- */
  /*  A   input                                                             */
  /* ---------------------------------------------------------------------- */
  if ((H == NULL) || (PsiBasis == NULL) || (G == NULL))
  {
    printf("CalculateGforHwithSigmaAsOperator called with a NULL pointer for H, PsiBasis or G\n");
    return 1;
  }
  if ((Npsi == 0) || (m == 0) || (m > Npsi))
  {
    printf("CalculateGforHwithSigmaAsOperator called with NPsiBasis=%u and BlockSize=%u, needs 0 < BlockSize <= NPsiBasis\n", Npsi, m);
    return 1;
  }
  if (order == 0)
  {
    printf("CalculateGforHwithSigmaAsOperator called with order=0, needs at least one block\n");
    return 1;
  }
  if (B0 == NULL)
  {
    NFull = m;
  }
  else
  {
    NFull = B0->Ni;
    if (B0->Nj != m)
    {
      printf("CalculateGforHwithSigmaAsOperator called with B0 of size %u x %u, the number of columns should be BlockSize=%u\n", B0->Ni, B0->Nj, m);
      return 1;
    }
    if (NFull < m)
    {
      printf("CalculateGforHwithSigmaAsOperator called with B0 of size %u x %u, the number of rows can not be smaller than BlockSize\n", B0->Ni, B0->Nj);
      return 1;
    }
  }
  for (i = 0; i < NSigma; i++)
  {
    if (Sigma[i].BlockSize == 0)
    {
      printf("CalculateGforHwithSigmaAsOperator called with Sigma[%u].BlockSize=0\n", i);
      return 1;
    }
    if ((V == NULL) || (V[i] == NULL))
    {
      printf("CalculateGforHwithSigmaAsOperator called with no operators for Sigma[%u]\n", i);
      return 1;
    }
  }
  QComplexOut = (H->QComplex != 0);
  for (a = 0; a < Npsi; a++) if (PsiBasis[a].QComplex) QComplexOut = 1;
  for (i = 0; i < NSigma; i++)
  {
    if (Sigma[i].QComplex) QComplexOut = 1;
    for (j = 0; j < Sigma[i].BlockSize; j++) if (V[i][j].QComplex) QComplexOut = 1;
  }
  if ((B0 != NULL) && (B0->QComplex)) QComplexOut = 1;

  /* ---------------------------------------------------------------------- */
  /*  A2  the constant term of a self energy is a static potential in P       */
  /*                                                                          */
  /*      Sigma[i](w) = a0[i] + sum_n W[i][n] / ( w - a[i][n] ) has, next to   */
  /*      the poles that become bath states, a frequency independent piece.    */
  /*      That piece is an ordinary one particle operator in P and has to act  */
  /*      wherever H acts -- in the outer space and just as much inside the    */
  /*      resolvent of the intermediate space.  Adding it to H once here is    */
  /*      both cheaper and more consistent than adding it to the outer matrix  */
  /*      alone, which would leave the intermediate spectrum unshifted.        */
  /* ---------------------------------------------------------------------- */
  VerbosityPrintf(0x0001, "CalculateGforHwithSigmaAsOperator: NPsiBasis=%u BlockSize=%u NSigma=%u order=%u\n", Npsi, m, NSigma, order);
  Hptr = (OperatorType*)H;
  for (i = 0; i < NSigma; i++)
  {
    uint32_t bi = Sigma[i].BlockSize;
    CompactMatrixType A0;
    OperatorType Oprod;
    A0.Ni = bi; A0.Nj = bi; A0.QComplex = 1;
    if (InitCompactMatrix(&A0)) { printf("InitCompactMatrix failed in CalculateGforHwithSigmaAsOperator\n"); ierr = 1; goto cleanup; }
    BlockListOfPolesGetA0(&(Sigma[i]), &A0);
    for (j = 0; j < bi; j++)
    {
      for (jp = 0; jp < bi; jp++)
      {
        if ((fabs(A0.V[j * bi + jp].Re) < 100 * DBL_MIN) && (fabs(A0.V[j * bi + jp].Im) < 100 * DBL_MIN)) continue;
        if (QHeffAllocated == false)
        {
          if (OperatorCopy(&Heff, *H)) { printf("OperatorCopy failed in CalculateGforHwithSigmaAsOperator\n"); FreeCompactMatrix(&A0); ierr = 1; goto cleanup; }
          QHeffAllocated = true;
          Hptr = &Heff;
        }
        InitOperatorForMult(&Oprod, 2, &(V[i][j]), &(V[i][jp]));
        OperatorMultiply2NormalOrder(&Oprod, V[i][j], V[i][jp], 'D', 'N', &QDetNotOrderedFcFaBcBa);
        OperatorAddNormalOrder(&Heff, Oprod, A0.V[j * bi + jp].Re, A0.V[j * bi + jp].Im, &QDetNotOrderedFcFaBcBa);
        FreeOperator(&Oprod);
      }
    }
    FreeCompactMatrix(&A0);
  }

  /* ---------------------------------------------------------------------- */
  /*  B   overlap and Hamiltonian in the given basis, then orthonormalize    */
  /* ---------------------------------------------------------------------- */
  SMat.Ni = Npsi; SMat.Nj = Npsi; SMat.QComplex = 1;
  HMat.Ni = Npsi; HMat.Nj = Npsi; HMat.QComplex = 1;
  if (InitCompactMatrix(&SMat) || InitCompactMatrix(&HMat))
  {
    printf("InitCompactMatrix failed in CalculateGforHwithSigmaAsOperator for the %u x %u outer matrices\n", Npsi, Npsi);
    ierr = 1; goto cleanup;
  }
  if (PsiListPsiListCompactMatrix(PsiBasis, Npsi, &SMat)) { ierr = 1; goto cleanup; }
  if (PsiListOperatorPsiListCompactMatrix(Hptr, PsiBasis, Npsi, &HMat)) { ierr = 1; goto cleanup; }
  SquareCompactMatrixMakeHermitian(&HMat);
  /*  Loewdin: |e_p> = sum_a T_{a,p} |PsiBasis[a]>, T_{a,p} = x_{p,a}/sqrt(s_p) */
  Sval = (double*)malloc(Npsi * sizeof(double));
  if (Sval == NULL) { printf("malloc failed in CalculateGforHwithSigmaAsOperator\n"); ierr = 1; goto cleanup; }
  if (HermitianCompactMatrixEigensystemRowVectors(&SMat, &Sfun, Sval)) { ierr = 1; goto cleanup; }
  North = 0;
  for (p = 0; p < Npsi; p++) if (Sval[p] > SigmaAsOperatorEpsOrtho * Sval[Npsi - 1]) North++;
  if (North < m)
  {
    printf("CalculateGforHwithSigmaAsOperator: the %u basis functions span only %u directions, fewer than the %u seeds\n", Npsi, North, m);
    ierr = 1; goto cleanup;
  }
  if (North != Npsi)
  {
    VerbosityPrintf(0x0001, "CalculateGforHwithSigmaAsOperator: PsiBasis lost orthogonality, %u of %u directions kept\n", North, Npsi);
  }
  T.Ni = Npsi; T.Nj = North; T.QComplex = 1;
  if (InitCompactMatrix(&T)) { printf("InitCompactMatrix failed in CalculateGforHwithSigmaAsOperator\n"); ierr = 1; goto cleanup; }
  q = 0;
  for (p = 0; p < Npsi; p++)
  {
    if (Sval[p] <= SigmaAsOperatorEpsOrtho * Sval[Npsi - 1]) continue;
    norm = 1.0 / sqrt(Sval[p]);
    for (a = 0; a < Npsi; a++)
    {
      T.V[a * North + q].Re = norm * Sfun.V[p * Npsi + a].Re;
      T.V[a * North + q].Im = norm * Sfun.V[p * Npsi + a].Im;
    }
    q++;
  }
  FreeCompactMatrix(&Sfun); CompactMatrixSetUnallocated(&Sfun);
  /*  HO = T^dag HMat T  */
  TS.Ni = North; TS.Nj = Npsi; TS.QComplex = 1;
  HO.Ni = North; HO.Nj = North; HO.QComplex = 1;
  if (InitCompactMatrix(&TS) || InitCompactMatrix(&HO)) { printf("InitCompactMatrix failed in CalculateGforHwithSigmaAsOperator\n"); ierr = 1; goto cleanup; }
  CompactMatrixAddMultiply('C', 'N', &T, &HMat, &TS, 1, 0, 0, 0);
  CompactMatrixAddMultiply('N', 'N', &TS, &T, &HO, 1, 0, 0, 0);
  SquareCompactMatrixMakeHermitian(&HO);
  /*  Cseed = ( T^dag SMat )[:,0:m], the seeds expressed in the new basis  */
  CompactMatrixAddMultiply('C', 'N', &T, &SMat, &TS, 1, 0, 0, 0);
  Cseed.Ni = North; Cseed.Nj = m; Cseed.QComplex = 1;
  if (InitCompactMatrix(&Cseed)) { printf("InitCompactMatrix failed in CalculateGforHwithSigmaAsOperator\n"); ierr = 1; goto cleanup; }
  for (p = 0; p < North; p++)
  {
    for (a = 0; a < m; a++)
    {
      Cseed.V[p * m + a] = TS.V[p * Npsi + a];
    }
  }
  FreeCompactMatrix(&HMat); CompactMatrixSetUnallocated(&HMat);
  FreeCompactMatrix(&SMat); CompactMatrixSetUnallocated(&SMat);
  FreeCompactMatrix(&TS); CompactMatrixSetUnallocated(&TS);

  /* ---------------------------------------------------------------------- */
  /*  C,D the intermediate space and the tensor Z, one per self energy       */
  /* ---------------------------------------------------------------------- */
  if (NSigma != 0)
  {
    Zlist = (CompactMatrixType*)calloc(NSigma, sizeof(CompactMatrixType));
    Elist = (double**)calloc(NSigma, sizeof(double*));
    Dlist = (uint32_t*)calloc(NSigma, sizeof(uint32_t));
    if ((Zlist == NULL) || (Elist == NULL) || (Dlist == NULL))
    {
      printf("malloc failed in CalculateGforHwithSigmaAsOperator\n");
      ierr = 1; goto cleanup;
    }
  }
  for (i = 0; i < NSigma; i++)
  {
    uint32_t bi = Sigma[i].BlockSize, Li = bi * North, di;
    CompactMatrixType SextO, HextO, Mtmp, MtmpO, TM, Xdefl, Hred, Sxfun, ufun, Ws, Chk;
    OperatorType Oprod;
    double *sxval = NULL;
    Li = bi * North;
    if (Li > SigmaAsOperatorMaxLi)
    {
      printf("CalculateGforHwithSigmaAsOperator: the extended space of self energy %u has dimension %u x %u = %u\n", i, bi, North, Li);
      printf("  which is above the compiled in limit SigmaAsOperatorMaxLi = %u.\n", SigmaAsOperatorMaxLi);
      printf("  Reduce NTri, the block size, or the number of channels of the self energy.\n");
      ierr = 1; goto cleanup;
    }
    SextO.Ni = Li;    SextO.Nj = Li;    SextO.QComplex = 1;
    HextO.Ni = Li;    HextO.Nj = Li;    HextO.QComplex = 1;
    Mtmp.Ni  = Npsi;  Mtmp.Nj  = Npsi;  Mtmp.QComplex  = 1;
    MtmpO.Ni = North; MtmpO.Nj = North; MtmpO.QComplex = 1;
    TM.Ni    = North; TM.Nj    = Npsi;  TM.QComplex    = 1;
    if (InitCompactMatrix(&SextO) || InitCompactMatrix(&HextO) || InitCompactMatrix(&Mtmp) ||
        InitCompactMatrix(&MtmpO) || InitCompactMatrix(&TM))
    {
      printf("InitCompactMatrix failed in CalculateGforHwithSigmaAsOperator for self energy %u\n", i);
      ierr = 1; goto cleanup;
    }
    /*  Sext and Hext.  V[i][j]^dag V[i][j'] and V[i][j]^dag H V[i][j'] act
     *  inside P, so no wave function outside P is ever created.               */
    for (j = 0; j < bi; j++)
    {
      for (jp = j; jp < bi; jp++)
      {
        InitOperatorForMult(&Oprod, 2, &(V[i][j]), &(V[i][jp]));
        OperatorMultiply2NormalOrder(&Oprod, V[i][j], V[i][jp], 'D', 'N', &QDetNotOrderedFcFaBcBa);
        if (PsiListOperatorPsiListCompactMatrix(&Oprod, PsiBasis, Npsi, &Mtmp)) { FreeOperator(&Oprod); ierr = 1; goto cleanup; }
        FreeOperator(&Oprod);
        CompactMatrixAddMultiply('C', 'N', &T, &Mtmp, &TM, 1, 0, 0, 0);
        CompactMatrixAddMultiply('N', 'N', &TM, &T, &MtmpO, 1, 0, 0, 0);
        for (p = 0; p < North; p++)
        {
          for (q = 0; q < North; q++)
          {
            SextO.V[(j * North + p) * Li + jp * North + q] = MtmpO.V[p * North + q];
            SextO.V[(jp * North + q) * Li + j * North + p].Re =  MtmpO.V[p * North + q].Re;
            SextO.V[(jp * North + q) * Li + j * North + p].Im = -MtmpO.V[p * North + q].Im;
          }
        }
        InitOperatorForMult(&Oprod, 3, &(V[i][j]), Hptr, &(V[i][jp]));
        OperatorMultiply3NormalOrder(&Oprod, V[i][j], *Hptr, V[i][jp], 'D', 'N', 'N', &QDetNotOrderedFcFaBcBa);
        if (PsiListOperatorPsiListCompactMatrix(&Oprod, PsiBasis, Npsi, &Mtmp)) { FreeOperator(&Oprod); ierr = 1; goto cleanup; }
        FreeOperator(&Oprod);
        CompactMatrixAddMultiply('C', 'N', &T, &Mtmp, &TM, 1, 0, 0, 0);
        CompactMatrixAddMultiply('N', 'N', &TM, &T, &MtmpO, 1, 0, 0, 0);
        for (p = 0; p < North; p++)
        {
          for (q = 0; q < North; q++)
          {
            HextO.V[(j * North + p) * Li + jp * North + q] = MtmpO.V[p * North + q];
            HextO.V[(jp * North + q) * Li + j * North + p].Re =  MtmpO.V[p * North + q].Re;
            HextO.V[(jp * North + q) * Li + j * North + p].Im = -MtmpO.V[p * North + q].Im;
          }
        }
      }
    }
    FreeCompactMatrix(&Mtmp);
    FreeCompactMatrix(&MtmpO);
    FreeCompactMatrix(&TM);
    /*  deflate Sext                                                          */
    sxval = (double*)malloc(Li * sizeof(double));
    if (sxval == NULL) { printf("malloc failed in CalculateGforHwithSigmaAsOperator\n"); ierr = 1; goto cleanup; }
    if (HermitianCompactMatrixEigensystemRowVectors(&SextO, &Sxfun, sxval)) { free(sxval); ierr = 1; goto cleanup; }
    di = 0;
    for (k = 0; k < Li; k++) if (sxval[k] > SigmaAsOperatorEpsDefl * sxval[Li - 1]) di++;
    if (di == 0)
    {
      VerbosityPrintf(0x0001, "CalculateGforHwithSigmaAsOperator: self energy %u does not couple to PsiBasis, skipped\n", i);
      Dlist[i] = 0;
      free(sxval);
      FreeCompactMatrix(&Sxfun);
      FreeCompactMatrix(&SextO);
      FreeCompactMatrix(&HextO);
      continue;
    }
    VerbosityPrintf(0x0001, "CalculateGforHwithSigmaAsOperator: self energy %u, extended space %u deflated to %u\n", i, Li, di);
    /*  Xdefl_{mu,alpha} = x_{alpha,mu} / sqrt(s_alpha)                        */
    Xdefl.Ni = Li; Xdefl.Nj = di; Xdefl.QComplex = 1;
    Ws.Ni = di; Ws.Nj = Li; Ws.QComplex = 1;
    Hred.Ni = di; Hred.Nj = di; Hred.QComplex = 1;
    TM.Ni = di; TM.Nj = Li; TM.QComplex = 1;
    if (InitCompactMatrix(&Xdefl) || InitCompactMatrix(&Ws) || InitCompactMatrix(&Hred) || InitCompactMatrix(&TM))
    {
      printf("InitCompactMatrix failed in CalculateGforHwithSigmaAsOperator for self energy %u\n", i);
      free(sxval); ierr = 1; goto cleanup;
    }
    q = 0;
    for (k = 0; k < Li; k++)
    {
      if (sxval[k] <= SigmaAsOperatorEpsDefl * sxval[Li - 1]) continue;
      norm = 1.0 / sqrt(sxval[k]);
      for (p = 0; p < Li; p++)
      {
        Xdefl.V[p * di + q].Re = norm * Sxfun.V[k * Li + p].Re;
        Xdefl.V[p * di + q].Im = norm * Sxfun.V[k * Li + p].Im;
        /*  Ws_{alpha,nu} = sqrt(s_alpha) conj( x_{alpha,nu} ) = <chi_alpha|w_nu> */
        Ws.V[q * Li + p].Re =  sxval[k] * norm * Sxfun.V[k * Li + p].Re;
        Ws.V[q * Li + p].Im = -sxval[k] * norm * Sxfun.V[k * Li + p].Im;
      }
      q++;
    }
    FreeCompactMatrix(&Sxfun);
    /*  Hred = Xdefl^dag Hext Xdefl, diagonalize                              */
    CompactMatrixAddMultiply('C', 'N', &Xdefl, &HextO, &TM, 1, 0, 0, 0);
    CompactMatrixAddMultiply('N', 'N', &TM, &Xdefl, &Hred, 1, 0, 0, 0);
    SquareCompactMatrixMakeHermitian(&Hred);
    FreeCompactMatrix(&TM);
    FreeCompactMatrix(&Xdefl);
    FreeCompactMatrix(&HextO);
    Elist[i] = (double*)malloc(di * sizeof(double));
    if (Elist[i] == NULL) { printf("malloc failed in CalculateGforHwithSigmaAsOperator\n"); free(sxval); ierr = 1; goto cleanup; }
    if (HermitianCompactMatrixEigensystemRowVectors(&Hred, &ufun, Elist[i])) { free(sxval); ierr = 1; goto cleanup; }
    FreeCompactMatrix(&Hred);
    free(sxval);
    /*  Z_{n',nu} = sum_alpha conj( u_{n',alpha} ) <chi_alpha|w_nu>            */
    Zlist[i].Ni = di; Zlist[i].Nj = Li; Zlist[i].QComplex = 1;
    if (InitCompactMatrix(&(Zlist[i]))) { printf("InitCompactMatrix failed in CalculateGforHwithSigmaAsOperator\n"); FreeCompactMatrix(&ufun); ierr = 1; goto cleanup; }
    CompactMatrixAddMultiply('*', 'N', &ufun, &Ws, &(Zlist[i]), 1, 0, 0, 0);
    FreeCompactMatrix(&ufun);
    FreeCompactMatrix(&Ws);
    Dlist[i] = di;
    /*  the zeroth moment of Sigma is exact by construction, so Z^dag Z has to
     *  reproduce Sext.  A failure here means the deflation above is wrong and
     *  everything downstream is meaningless, so it is always checked.         */
    Chk.Ni = Li; Chk.Nj = Li; Chk.QComplex = 1;
    if (InitCompactMatrix(&Chk)) { printf("InitCompactMatrix failed in CalculateGforHwithSigmaAsOperator\n"); ierr = 1; goto cleanup; }
    CompactMatrixAddMultiply('C', 'N', &(Zlist[i]), &(Zlist[i]), &Chk, 1, 0, 0, 0);
    norm = 0;
    Re = 0;
    for (k = 0; k < Li * Li; k++)
    {
      norm = Max(norm, fabs(Chk.V[k].Re - SextO.V[k].Re) + fabs(Chk.V[k].Im - SextO.V[k].Im));
      Re = Max(Re, fabs(SextO.V[k].Re) + fabs(SextO.V[k].Im));
    }
    if (norm > 1e-8 * Max(Re, 1.0))
    {
      printf("CalculateGforHwithSigmaAsOperator: the zeroth moment sum rule of self energy %u is violated by %g\n", i, norm);
      printf("  Z^dag Z should reproduce the extended overlap matrix.  The result is not to be trusted.\n");
    }
    FreeCompactMatrix(&Chk);
    FreeCompactMatrix(&SextO);
  }

  /* ---------------------------------------------------------------------- */
  /*  E   the bath of the arrowhead                                          */
  /* ---------------------------------------------------------------------- */
  /*  the range the poles cover, how many there are, and whether the residues
   *  are positive semi definite.  A residue that is not has no Hamiltonian
   *  representation: the arrowhead below would not be Hermitian and the block
   *  Lanczos on it would be invalid, not merely inaccurate.                    */
  Emin = DBL_MAX;
  Emax = -DBL_MAX;
  NbathMax = 0;
  for (i = 0; i < NSigma; i++)
  {
    uint32_t bi = Sigma[i].BlockSize;
    CompactMatrixType Wm, Wfun;
    double *wval, en, wmax;
    if (Dlist[i] == 0) continue;
    Wm.Ni = bi; Wm.Nj = bi; Wm.QComplex = 1;
    wval = (double*)malloc(bi * sizeof(double));
    if (InitCompactMatrix(&Wm) || (wval == NULL)) { printf("malloc failed in CalculateGforHwithSigmaAsOperator\n"); ierr = 1; goto cleanup; }
    for (n = 0; n < Sigma[i].N; n++)
    {
      BlockListOfPolesGetPole(&(Sigma[i]), n, &Wm, &en);
      for (np = 0; np < Dlist[i]; np++)
      {
        Emin = Min(Emin, en + Elist[i][np]);
        Emax = Max(Emax, en + Elist[i][np]);
      }
      NbathMax += Dlist[i] * bi;
      if (HermitianCompactMatrixEigensystemRowVectors(&Wm, &Wfun, wval)) { free(wval); FreeCompactMatrix(&Wm); ierr = 1; goto cleanup; }
      FreeCompactMatrix(&Wfun);
      wmax = Max(wval[bi - 1], 0.0);
      for (k = 0; k < bi; k++)
      {
        if (wval[k] < -SigmaAsOperatorEpsPSD * Max(wmax, DBL_MIN)) NegativeWeight -= wval[k];
      }
    }
    free(wval);
    FreeCompactMatrix(&Wm);
  }
  if (NbathMax != 0)
  {
    bool UseGrid = false;
    if ((SuggestedGrid != NULL) && (NSuggestedGridPoints >= 2)) UseGrid = true;
    if (NbathMax > SigmaAsOperatorMaxBathExact) UseGrid = true;
    if (UseGrid == false)
    {
      /*  E exact.  Every pole of every self energy, split by every state of the
       *  intermediate space, becomes its own bath state.                      */
      CompactMatrixType Wm, Wfun, Lt, Gm, Ysub;
      double *wval, wmax;
      uint32_t r;
      Ybath.Ni = NbathMax; Ybath.Nj = North; Ybath.QComplex = 1;
      if (InitCompactMatrix(&Ybath)) { printf("InitCompactMatrix failed in CalculateGforHwithSigmaAsOperator for the bath\n"); ierr = 1; goto cleanup; }
      Ebath = (double*)calloc(NbathMax, sizeof(double));
      if (Ebath == NULL) { printf("malloc failed in CalculateGforHwithSigmaAsOperator\n"); ierr = 1; goto cleanup; }
      Nbath = 0;
      for (i = 0; i < NSigma; i++)
      {
        uint32_t bi = Sigma[i].BlockSize, Li = bi * North;
        if (Dlist[i] == 0) continue;
        Wm.Ni = bi; Wm.Nj = bi; Wm.QComplex = 1;
        Lt.Ni = bi; Lt.Nj = bi; Lt.QComplex = 1;
        Ysub.Ni = bi; Ysub.Nj = North; Ysub.QComplex = 1;
        wval = (double*)malloc(bi * sizeof(double));
        if (InitCompactMatrix(&Wm) || InitCompactMatrix(&Lt) || InitCompactMatrix(&Ysub) || (wval == NULL))
        {
          printf("InitCompactMatrix failed in CalculateGforHwithSigmaAsOperator\n"); ierr = 1; goto cleanup;
        }
        for (n = 0; n < Sigma[i].N; n++)
        {
          double en;
          BlockListOfPolesGetPole(&(Sigma[i]), n, &Wm, &en);
          if (HermitianCompactMatrixEigensystemRowVectors(&Wm, &Wfun, wval)) { ierr = 1; goto cleanup; }
          wmax = Max(wval[bi - 1], 0.0);
          r = 0;
          for (k = 0; k < bi; k++)
          {
            if (wval[k] > SigmaAsOperatorEpsChol * wmax) r++;
            else if (wval[k] > 0) DiscardedWeight += wval[k];
          }
          if (r == 0) { FreeCompactMatrix(&Wfun); continue; }
          Lt.Ni = r;
          Ysub.Ni = r;
          q = 0;
          for (k = 0; k < bi; k++)
          {
            if (wval[k] <= SigmaAsOperatorEpsChol * wmax) continue;
            norm = sqrt(wval[k]);
            for (j = 0; j < bi; j++)
            {
              Lt.V[q * bi + j].Re =  norm * Wfun.V[k * bi + j].Re;
              Lt.V[q * bi + j].Im = -norm * Wfun.V[k * bi + j].Im;
            }
            q++;
          }
          FreeCompactMatrix(&Wfun);
          for (np = 0; np < Dlist[i]; np++)
          {
            Gm.Ni = bi; Gm.Nj = North; Gm.QComplex = 1; Gm.OnLuaStack = true;
            Gm.Name = Zlist[i].Name;
            Gm.a0 = 0; Gm.b0 = 1; Gm.mu = 0;
            Gm.V = &(Zlist[i].V[(size_t)np * Li]);
            CompactMatrixAddMultiply('N', 'N', &Lt, &Gm, &Ysub, 1, 0, 0, 0);
            for (k = 0; k < r; k++)
            {
              Ebath[Nbath] = en + Elist[i][np];
              for (p = 0; p < North; p++) Ybath.V[(size_t)Nbath * North + p] = Ysub.V[k * North + p];
              Nbath++;
            }
          }
        }
        free(wval);
        FreeCompactMatrix(&Wm);
        FreeCompactMatrix(&Lt);
        FreeCompactMatrix(&Ysub);
      }
      Ybath.Ni = Nbath;
      VerbosityPrintf(0x0001, "CalculateGforHwithSigmaAsOperator: %u bath states, poles kept exactly\n", Nbath);
    }
    else
    {
      /*  E binned.  The zeroth and the first moment of every bin are exact, the
       *  second is off by O(dg^2), and the bath is capped at Ng * North.       */
      CompactMatrixType Wm, Gm, WG, Yc;
      uint32_t c, clow, chigh;
      double uleft, uright, en;
      if ((SuggestedGrid != NULL) && (NSuggestedGridPoints >= 2))
      {
        for (k = 1; k < NSuggestedGridPoints; k++)
        {
          if (SuggestedGrid[k] <= SuggestedGrid[k - 1])
          {
            printf("CalculateGforHwithSigmaAsOperator: the grid is not strictly increasing at position %u (%g <= %g)\n", k, SuggestedGrid[k], SuggestedGrid[k - 1]);
            ierr = 1; goto cleanup;
          }
        }
        Ng = NSuggestedGridPoints + ((Emin < SuggestedGrid[0]) ? 1 : 0) + ((Emax > SuggestedGrid[NSuggestedGridPoints - 1]) ? 1 : 0);
        Grid = (double*)malloc(Ng * sizeof(double));
        if (Grid == NULL) { printf("malloc failed in CalculateGforHwithSigmaAsOperator\n"); ierr = 1; goto cleanup; }
        k = 0;
        if (Emin < SuggestedGrid[0])
        {
          Grid[k++] = Emin;
          VerbosityPrintf(0x0001, "CalculateGforHwithSigmaAsOperator: grid extended down to %g to hold all poles\n", Emin);
        }
        for (q = 0; q < NSuggestedGridPoints; q++) Grid[k++] = SuggestedGrid[q];
        if (Emax > SuggestedGrid[NSuggestedGridPoints - 1])
        {
          Grid[k++] = Emax;
          VerbosityPrintf(0x0001, "CalculateGforHwithSigmaAsOperator: grid extended up to %g to hold all poles\n", Emax);
        }
      }
      else
      {
        Ng = Max(2, 8 * order);
        if (Ng > NbathMax) Ng = Max(2, NbathMax);
        if ((double)Ng * (double)North * (double)North > SigmaAsOperatorMaxGridMemory)
        {
          Ng = (uint32_t)(SigmaAsOperatorMaxGridMemory / ((double)North * (double)North));
          Ng = Max(2, Ng);
          printf("CalculateGforHwithSigmaAsOperator: the automatic grid is cut back to %u points to stay inside the memory limit\n", Ng);
          printf("  Supply GridPoints if a different distribution is wanted.\n");
        }
        if (Emax - Emin < 100 * DBL_MIN) { Emin -= 0.5; Emax += 0.5; }
        Grid = (double*)malloc(Ng * sizeof(double));
        if (Grid == NULL) { printf("malloc failed in CalculateGforHwithSigmaAsOperator\n"); ierr = 1; goto cleanup; }
        for (k = 0; k < Ng; k++) Grid[k] = Emin + (Emax - Emin) * (double)k / (double)(Ng - 1);
      }
      VerbosityPrintf(0x0001, "CalculateGforHwithSigmaAsOperator: %u poles binned on %u grid points\n", NbathMax, Ng);
      Rgrid = (CompactMatrixType*)calloc(Ng, sizeof(CompactMatrixType));
      if (Rgrid == NULL) { printf("malloc failed in CalculateGforHwithSigmaAsOperator\n"); ierr = 1; goto cleanup; }
      for (c = 0; c < Ng; c++)
      {
        Rgrid[c].Ni = North; Rgrid[c].Nj = North; Rgrid[c].QComplex = 1;
        if (InitCompactMatrix(&(Rgrid[c])))
        {
          printf("InitCompactMatrix failed in CalculateGforHwithSigmaAsOperator for grid point %u of %u\n", c, Ng);
          printf("  %u accumulators of %u x %u do not fit in memory, use fewer grid points.\n", Ng, North, North);
          ierr = 1; goto cleanup;
        }
      }
      for (i = 0; i < NSigma; i++)
      {
        uint32_t bi = Sigma[i].BlockSize, Li = bi * North;
        if (Dlist[i] == 0) continue;
        Wm.Ni = bi; Wm.Nj = bi; Wm.QComplex = 1;
        WG.Ni = bi; WG.Nj = North; WG.QComplex = 1;
        if (InitCompactMatrix(&Wm) || InitCompactMatrix(&WG)) { printf("InitCompactMatrix failed in CalculateGforHwithSigmaAsOperator\n"); ierr = 1; goto cleanup; }
        for (np = 0; np < Dlist[i]; np++)
        {
          Gm.Ni = bi; Gm.Nj = North; Gm.QComplex = 1; Gm.OnLuaStack = true;
          Gm.Name = Zlist[i].Name;
          Gm.a0 = 0; Gm.b0 = 1; Gm.mu = 0;
          Gm.V = &(Zlist[i].V[(size_t)np * Li]);
          for (n = 0; n < Sigma[i].N; n++)
          {
            BlockListOfPolesGetPole(&(Sigma[i]), n, &Wm, &en);
            en += Elist[i][np];
            /*  bisection for the interval that holds en                        */
            clow = 0;
            chigh = Ng - 1;
            while (chigh - clow > 1)
            {
              c = (clow + chigh) / 2;
              if (Grid[c] <= en) clow = c; else chigh = c;
            }
            FindWeigth(&(Grid[clow]), &en, &(Grid[chigh]), &uleft, &uright);
            CompactMatrixAddMultiply('N', 'N', &Wm, &Gm, &WG, 1, 0, 0, 0);
            CompactMatrixAddMultiply('C', 'N', &Gm, &WG, &(Rgrid[clow]),  uleft, 0, 1, 0);
            CompactMatrixAddMultiply('C', 'N', &Gm, &WG, &(Rgrid[chigh]), uright, 0, 1, 0);
          }
        }
        FreeCompactMatrix(&Wm);
        FreeCompactMatrix(&WG);
      }
      /*  factorize every bin.  R_c is a positively weighted sum of G^dag W G
       *  and is positive semi definite by construction, so the truncation only
       *  ever throws weight away and never creates negative weight.            */
      Nbath = 0;
      for (c = 0; c < Ng; c++)
      {
        SquareCompactMatrixMakeHermitian(&(Rgrid[c]));
        if (HermitianCompactMatrixToRectangularRoot(&(Rgrid[c]), &Yc, SigmaAsOperatorEpsChol, &DiscardedWeight, &NegativeWeight)) { ierr = 1; goto cleanup; }
        FreeCompactMatrix(&(Rgrid[c]));
        Rgrid[c] = Yc;
        Nbath += Yc.Ni;
      }
      Ybath.Ni = Max(Nbath, 1); Ybath.Nj = North; Ybath.QComplex = 1;
      if (InitCompactMatrix(&Ybath)) { printf("InitCompactMatrix failed in CalculateGforHwithSigmaAsOperator for the bath\n"); ierr = 1; goto cleanup; }
      Ybath.Ni = Nbath;
      Ebath = (double*)calloc(Max(Nbath, 1), sizeof(double));
      if (Ebath == NULL) { printf("malloc failed in CalculateGforHwithSigmaAsOperator\n"); ierr = 1; goto cleanup; }
      k = 0;
      for (c = 0; c < Ng; c++)
      {
        for (q = 0; q < Rgrid[c].Ni; q++)
        {
          Ebath[k] = Grid[c];
          for (p = 0; p < North; p++) Ybath.V[(size_t)k * North + p] = Rgrid[c].V[q * North + p];
          k++;
        }
        FreeCompactMatrix(&(Rgrid[c]));
      }
      free(Rgrid); Rgrid = NULL;
      VerbosityPrintf(0x0001, "CalculateGforHwithSigmaAsOperator: %u bath states after binning\n", Nbath);
    }
  }
  if (NegativeWeight > 0)
  {
    printf("CalculateGforHwithSigmaAsOperator: %g negative spectral weight found in the residues of Sigma\n", NegativeWeight);
    printf("  A self energy with a residue that is not positive semi definite has no Hamiltonian representation.\n");
    printf("  The negative part is dropped, the result is an approximation of an unphysical input.\n");
  }
  if (DiscardedWeight > 0)
  {
    VerbosityPrintf(0x0001, "CalculateGforHwithSigmaAsOperator: %g spectral weight below the truncation threshold discarded\n", DiscardedWeight);
  }

  /* ---------------------------------------------------------------------- */
  /*  F   block Lanczos on the arrowhead                                     */
  /*                                                                         */
  /*      M = (  HO      Ybath^dag             )                             */
  /*          (  Ybath   diag( Ebath )         )                             */
  /*                                                                         */
  /*      the Schur complement of which on the outer block is                */
  /*      w - HO - Sigma(w), so the resolvent of M in the seed block is G.    */
  /*      The structure is used in the matrix vector product, M is never      */
  /*      formed.                                                             */
  /* ---------------------------------------------------------------------- */
  Dtot = North + Nbath;
  orderUse = order;
  if ((double)orderUse * (double)m > (double)Dtot) orderUse = (Dtot + m - 1) / m;
  if (orderUse == 0) orderUse = 1;
  if (orderUse != order)
  {
    VerbosityPrintf(0x0001, "CalculateGforHwithSigmaAsOperator: order reduced from %u to %u, the space has only %u dimensions\n", order, orderUse, Dtot);
  }
  {
    CompactMatrixType Qoall, Qball, Yo, Yb, Qok, Qbk, Qoj, Qbj, Qonext, Qbnext;
    uint32_t l, t;
    Qoall.Ni = orderUse * m; Qoall.Nj = North;        Qoall.QComplex = 1;
    Qball.Ni = orderUse * m; Qball.Nj = Max(Nbath, 1); Qball.QComplex = 1;
    Yo.Ni = m; Yo.Nj = North;         Yo.QComplex = 1;
    Yb.Ni = m; Yb.Nj = Max(Nbath, 1); Yb.QComplex = 1;
    Ablk.Ni = m; Ablk.Nj = m; Ablk.QComplex = 1;
    Bblk.Ni = m; Bblk.Nj = m; Bblk.QComplex = 1;
    Lblk.Ni = m; Lblk.Nj = m; Lblk.QComplex = 1;
    if (InitCompactMatrix(&Qoall) || InitCompactMatrix(&Qball) || InitCompactMatrix(&Yo) ||
        InitCompactMatrix(&Yb) || InitCompactMatrix(&Ablk) || InitCompactMatrix(&Bblk) || InitCompactMatrix(&Lblk))
    {
      printf("InitCompactMatrix failed in CalculateGforHwithSigmaAsOperator for the block Lanczos\n");
      printf("  %u vectors of length %u do not fit in memory, use a smaller order.\n", orderUse * m, Dtot);
      FreeCompactMatrix(&Qoall); FreeCompactMatrix(&Qball); FreeCompactMatrix(&Yo); FreeCompactMatrix(&Yb);
      ierr = 1; goto cleanup;
    }
    Qball.Nj = Nbath;
    Yb.Nj = Nbath;
    chain.N = orderUse;
    chain.BlockSize = m;
    chain.QComplex = 1;
    if (InitBlockTridiagonalMatrix(&chain))
    {
      printf("InitBlockTridiagonalMatrix failed in CalculateGforHwithSigmaAsOperator\n");
      FreeCompactMatrix(&Qoall); FreeCompactMatrix(&Qball); FreeCompactMatrix(&Yo); FreeCompactMatrix(&Yb);
      ierr = 1; goto cleanup;
    }
    chain.mu = 0;
    /*  the seed block, the m seeds in the orthonormalized outer basis         */
    for (l = 0; l < m; l++)
    {
      for (p = 0; p < North; p++) Qoall.V[l * North + p] = Cseed.V[p * m + l];
    }
    Qok = CompactMatrixRowBlockView(&Qoall, 0, m);
    Qbk = CompactMatrixRowBlockView(&Qball, 0, m);
    RowBlockPairOrthonormalizeWithTransform(&Qok, &Qbk, &Lblk, SigmaAsOperatorEpsLanczos);
    memcpy(chain.b, Lblk.V, Sqr(m) * sizeof(Complex16));
    for (k = 0; k < orderUse; k++)
    {
      Qok = CompactMatrixRowBlockView(&Qoall, k * m, m);
      Qbk = CompactMatrixRowBlockView(&Qball, k * m, m);
      /*  Y = M Q_k                                                           */
      CompactMatrixAddMultiply('N', 'T', &Qok, &HO, &Yo, 1, 0, 0, 0);
      if (Nbath)
      {
        CompactMatrixAddMultiply('N', '*', &Qbk, &Ybath, &Yo, 1, 0, 1, 0);
        CompactMatrixAddMultiply('N', 'T', &Qok, &Ybath, &Yb, 1, 0, 0, 0);
        for (l = 0; l < m; l++)
        {
          for (t = 0; t < Nbath; t++)
          {
            Yb.V[(size_t)l * Nbath + t].Re += Ebath[t] * Qbk.V[(size_t)l * Nbath + t].Re;
            Yb.V[(size_t)l * Nbath + t].Im += Ebath[t] * Qbk.V[(size_t)l * Nbath + t].Im;
          }
        }
      }
      /*  A_k = Q_k^dag M Q_k                                                 */
      CompactMatrixAddMultiply('*', 'T', &Qok, &Yo, &Ablk, 1, 0, 0, 0);
      if (Nbath) CompactMatrixAddMultiply('*', 'T', &Qbk, &Yb, &Ablk, 1, 0, 1, 0);
      SquareCompactMatrixMakeHermitian(&Ablk);
      memcpy(&(chain.a[(size_t)(k + 1) * Sqr(m)]), Ablk.V, Sqr(m) * sizeof(Complex16));
      if (k + 1 == orderUse) break;
      /*  the residual                                                        */
      CompactMatrixAddMultiply('T', 'N', &Ablk, &Qok, &Yo, -1, 0, 1, 0);
      if (Nbath) CompactMatrixAddMultiply('T', 'N', &Ablk, &Qbk, &Yb, -1, 0, 1, 0);
      if (k != 0)
      {
        Qoj = CompactMatrixRowBlockView(&Qoall, (k - 1) * m, m);
        Qbj = CompactMatrixRowBlockView(&Qball, (k - 1) * m, m);
        memcpy(Bblk.V, &(chain.b[(size_t)k * Sqr(m)]), Sqr(m) * sizeof(Complex16));
        CompactMatrixAddMultiply('T', 'N', &Bblk, &Qoj, &Yo, -1, 0, 1, 0);
        if (Nbath) CompactMatrixAddMultiply('T', 'N', &Bblk, &Qbj, &Yb, -1, 0, 1, 0);
      }
      /*  full re-orthogonalization, twice; the vectors are cheap here         */
      for (q = 0; q < 2; q++)
      {
        for (j = 0; j <= k; j++)
        {
          Qoj = CompactMatrixRowBlockView(&Qoall, j * m, m);
          Qbj = CompactMatrixRowBlockView(&Qball, j * m, m);
          CompactMatrixAddMultiply('*', 'T', &Qoj, &Yo, &Bblk, 1, 0, 0, 0);
          if (Nbath) CompactMatrixAddMultiply('*', 'T', &Qbj, &Yb, &Bblk, 1, 0, 1, 0);
          CompactMatrixAddMultiply('T', 'N', &Bblk, &Qoj, &Yo, -1, 0, 1, 0);
          if (Nbath) CompactMatrixAddMultiply('T', 'N', &Bblk, &Qbj, &Yb, -1, 0, 1, 0);
        }
      }
      /*  normalize what is left, zero rows keep the block size uniform        */
      norm = 0;
      for (l = 0; l < m; l++)
      {
        Re = 0;
        for (p = 0; p < North; p++) Re += Sqr(Yo.V[l * North + p].Re) + Sqr(Yo.V[l * North + p].Im);
        for (t = 0; t < Nbath; t++) Re += Sqr(Yb.V[(size_t)l * Nbath + t].Re) + Sqr(Yb.V[(size_t)l * Nbath + t].Im);
        norm = Max(norm, sqrt(Re));
      }
      RowBlockPairOrthonormalizeWithTransform(&Yo, &Yb, &Lblk, SigmaAsOperatorEpsLanczos * Max(norm, DBL_MIN));
      Qonext = CompactMatrixRowBlockView(&Qoall, (k + 1) * m, m);
      Qbnext = CompactMatrixRowBlockView(&Qball, (k + 1) * m, m);
      memcpy(Qonext.V, Yo.V, (size_t)m * North * sizeof(Complex16));
      if (Nbath) memcpy(Qbnext.V, Yb.V, (size_t)m * Nbath * sizeof(Complex16));
      /*  R_old[i] = sum_j Q_{k+1}[j] L[i][j] and B_k[i][j] = <v_k^i|M|v_{k+1}^j>
       *  = conj( L[i][j] ), the block that sits at position (k,k+1)            */
      for (p = 0; p < m; p++)
      {
        for (q = 0; q < m; q++)
        {
          chain.b[(size_t)(k + 1) * Sqr(m) + p * m + q].Re =  Lblk.V[p * m + q].Re;
          chain.b[(size_t)(k + 1) * Sqr(m) + p * m + q].Im = -Lblk.V[p * m + q].Im;
        }
      }
    }
    FreeCompactMatrix(&Qoall);
    FreeCompactMatrix(&Qball);
    FreeCompactMatrix(&Yo);
    FreeCompactMatrix(&Yb);
  }

  /* ---------------------------------------------------------------------- */
  /*  G   chain to star and pad back up to the full block size               */
  /* ---------------------------------------------------------------------- */
  if (BlockTriDiagonalMatrixToBlockAndersonMatrix(chain, &GSmall))
  {
    printf("BlockTriDiagonalMatrixToBlockAndersonMatrix failed in CalculateGforHwithSigmaAsOperator\n");
    ierr = 1; goto cleanup;
  }
  G->N = GSmall.N;
  G->BlockSize = NFull;
  G->QComplex = 1;
  if (InitBlockAndersonMatrix(G))
  {
    printf("InitBlockAndersonMatrix failed in CalculateGforHwithSigmaAsOperator\n");
    ierr = 1; goto cleanup;
  }
  G->mu = 0;
  strcpy(G->Name, "G with self energy");
  for (k = 0; k <= G->N; k++)
  {
    for (p = 0; p < m; p++)
    {
      for (q = 0; q < m; q++)
      {
        G->a[(size_t)k * Sqr(NFull) + p * NFull + q] = GSmall.a[(size_t)k * Sqr(m) + p * m + q];
      }
    }
  }
  for (k = 1; k < G->N; k++)
  {
    for (p = 0; p < m; p++)
    {
      for (q = 0; q < m; q++)
      {
        G->b[(size_t)k * Sqr(NFull) + p * NFull + q] = GSmall.b[(size_t)k * Sqr(m) + p * m + q];
      }
    }
  }
  /*  b[0] carries the deflated seed block back to the full set of transition
   *  operators.  B0 has the rank deficiency of the seeds in it, so its columns
   *  beyond BlockSize are zero and the padded part of the chain never shows.   */
  for (p = 0; p < NFull; p++)
  {
    for (q = 0; q < NFull; q++)
    {
      G->b[p * NFull + q].Re = 0;
      G->b[p * NFull + q].Im = 0;
    }
  }
  for (p = 0; p < NFull; p++)
  {
    for (q = 0; q < m; q++)
    {
      for (k = 0; k < m; k++)
      {
        Complex16 b0pk;
        if (B0 == NULL)
        {
          b0pk.Re = (p == k) ? 1 : 0;
          b0pk.Im = 0;
        }
        else if (B0->QComplex)
        {
          b0pk = B0->V[p * m + k];
        }
        else
        {
          b0pk.Re = B0->VRe[p * m + k];
          b0pk.Im = 0;
        }
        G->b[p * NFull + q].Re += b0pk.Re * GSmall.b[k * m + q].Re - b0pk.Im * GSmall.b[k * m + q].Im;
        G->b[p * NFull + q].Im += b0pk.Re * GSmall.b[k * m + q].Im + b0pk.Im * GSmall.b[k * m + q].Re;
      }
    }
  }
  if (QComplexOut == 0)
  {
    norm = 0;
    Re = 0;
    for (k = 0; k < (G->N + 1) * Sqr(NFull); k++) { norm = Max(norm, fabs(G->a[k].Im)); Re = Max(Re, fabs(G->a[k].Re)); }
    for (k = 0; k < G->N * Sqr(NFull); k++)       { norm = Max(norm, fabs(G->b[k].Im)); Re = Max(Re, fabs(G->b[k].Re)); }
    if (norm > 1e-10 * Max(Re, 1.0))
    {
      printf("CalculateGforHwithSigmaAsOperator: all input is real but the result has an imaginary part of %g, kept complex\n", norm);
    }
    else
    {
      ComplexBlockAndersonMatrixToReal(G);
    }
  }

cleanup:
  FreeCompactMatrix(&SMat);
  FreeCompactMatrix(&HMat);
  FreeCompactMatrix(&Sfun);
  FreeCompactMatrix(&T);
  FreeCompactMatrix(&TS);
  FreeCompactMatrix(&HO);
  FreeCompactMatrix(&Cseed);
  FreeCompactMatrix(&Ybath);
  FreeCompactMatrix(&Ablk);
  FreeCompactMatrix(&Bblk);
  FreeCompactMatrix(&Lblk);
  if (chain.a != NULL) FreeBlockTridiagonalMatrix(&chain);
  if (GSmall.a != NULL) FreeBlockAndersonMatrix(&GSmall);
  if (Rgrid != NULL)
  {
    for (k = 0; k < Ng; k++) FreeCompactMatrix(&(Rgrid[k]));
    free(Rgrid);
  }
  if (Zlist != NULL)
  {
    for (i = 0; i < NSigma; i++) FreeCompactMatrix(&(Zlist[i]));
    free(Zlist);
  }
  if (Elist != NULL)
  {
    for (i = 0; i < NSigma; i++) free(Elist[i]);
    free(Elist);
  }
  free(Dlist);
  if (QHeffAllocated) FreeOperator(&Heff);
  free(Sval);
  free(Grid);
  free(Ebath);
  return ierr;
}
