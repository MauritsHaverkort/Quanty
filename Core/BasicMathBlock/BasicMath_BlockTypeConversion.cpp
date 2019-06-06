/*  These functions convert one type of representation to an other
 *
 *  Conversions are:
 *    Evaluate__ (at one specific point)
 *    __ToSpectra (return a spectrum object)
 */


/*
 *  Evaluate__ calcualtes the Green's function at one specific point
 */

/*
 *
 */
uint32_t EvaluateBlockListOfPoles(BlockListOfPoles const *A, CompactMatrixType * Gw, const double omega, const double Gamma)
{
  int i,iblock,itot;
  Complex16 Lorentzian;
  Gw->QComplex=1;
  Gw->Ni=A->BlockSize;
  Gw->Nj=A->BlockSize;
  Gw->a0=0;
  Gw->b0=0;
  Gw->mu=A->mu;
  InitCompactMatrix(Gw);
  if(A->QComplex==0)
  {// real
    // add A0
    for(i=0;i<A->BlockSize*A->BlockSize;i++)
    {
      Gw->V[i].Re=A->VRe[A->N*(1+A->BlockSize*A->BlockSize)+i];
      Gw->V[i].Im=0;
    }
    // loop over all poles
    itot=0;
    for(iblock=0;iblock<A->N;iblock++)
    {
      Lorentzian.Re = (omega-A->VRe[itot]) / (Sqr(Gamma/2.0) + Sqr(omega-A->VRe[itot]));
      Lorentzian.Im = (-Gamma/2.0        ) / (Sqr(Gamma/2.0) + Sqr(omega-A->VRe[itot]));
      itot++;
      for(i=0;i<A->BlockSize*A->BlockSize;i++)
      {
        Gw->V[i].Re += A->VRe[itot] * Lorentzian.Re;
        Gw->V[i].Im += A->VRe[itot] * Lorentzian.Im;
        itot++;
      }
    }
  }
  else
  {// complex
    // add A0
    for(i=0;i<A->BlockSize*A->BlockSize;i++)
    {
      Gw->V[i].Re=A->V[A->N*(1+A->BlockSize*A->BlockSize)+i].Re;
      Gw->V[i].Im=A->V[A->N*(1+A->BlockSize*A->BlockSize)+i].Im;
    }
    // loop over all poles
    itot=0;
    for(iblock=0;iblock<A->N;iblock++)
    {
      Lorentzian.Re = (omega-A->V[itot].Re) / (Sqr(Gamma) + Sqr(omega-A->V[itot].Re));
      Lorentzian.Im = (-Gamma             ) / (Sqr(Gamma) + Sqr(omega-A->V[itot].Re));
      itot++;
      for(i=0;i<A->BlockSize*A->BlockSize;i++)
      {
        Gw->V[i].Re += A->V[itot].Re * Lorentzian.Re - A->V[itot].Im * Lorentzian.Im;
        Gw->V[i].Im += A->V[itot].Re * Lorentzian.Im + A->V[itot].Im * Lorentzian.Re;
        itot++;
      }
    }
  }
  return 0;
}

uint32_t EvaluateBlockListOfPoles(void const *A, CompactMatrixType * Gw, const double omega, const double Gamma)
{
	BlockListOfPoles *ATypecast = (BlockListOfPoles*) A;
	return EvaluateBlockListOfPoles(ATypecast, Gw, omega, Gamma);
}


/*
 *  G(w) = A0 + B0^* ( 1 / ( w - A1 - \sum_i ( B_i 1 / ( w - A_{i+1} ) . B_i^+ ) ) ) . B0^T
 */
uint32_t EvaluateBlockAndersonMatrix(BlockAndersonMatrixType const *A, CompactMatrixType * Gw, const double omega, const double Gamma)
{
  uint32_t i,iblock;
  CompactMatrixType a,b,Gwi;
  a.Ni=A->BlockSize;
  a.Nj=A->BlockSize;
  b.Ni=A->BlockSize;
  b.Nj=A->BlockSize;
  a.QComplex=A->QComplex;
  b.QComplex=A->QComplex;
  a.a0=0;
  a.b0=0;
  a.mu=0;
  b.a0=0;
  b.b0=0;
  b.mu=0;
  Gw->QComplex=1;
  Gw->Ni=A->BlockSize;
  Gw->Nj=A->BlockSize;
  Gw->a0=0;
  Gw->b0=0;
  Gw->mu=A->mu;
  Gwi.QComplex=1;
  Gwi.Ni=A->BlockSize;
  Gwi.Nj=A->BlockSize;
  Gwi.a0=0;
  Gwi.b0=0;
  Gwi.mu=A->mu;
  InitCompactMatrix(Gw);
  InitCompactMatrix(&Gwi);
  for(i=0;i<A->BlockSize*A->BlockSize;i++)
  {
    Gw->V[i].Re=0;
    Gw->V[i].Im=0;
  }
  // calculate Gi = b_{i-1} . 1/(w-a_i) . b_{i-1}^+ and add Gi to Gw
  for(iblock=A->N;iblock>1;iblock--)
  {
    /* set a and b to point to right matrices */
    if(A->QComplex)
    {
      a.V=&(A->a[iblock    *A->BlockSize*A->BlockSize]);
      b.V=&(A->b[(iblock-1)*A->BlockSize*A->BlockSize]);
    }
    else
    {
      a.VRe=&(A->aRe[iblock    *A->BlockSize*A->BlockSize]);
      b.VRe=&(A->bRe[(iblock-1)*A->BlockSize*A->BlockSize]);
    }
    // Gi = b_{i-1} . 1/(w-a_i) . b_{i-1}^+
    SetCompactMatrixToZero(&Gwi);
    CompactMatrixAdd(&Gwi, &a);
    for(i=0;i<A->BlockSize*A->BlockSize;i++)
    {
      Gwi.V[i].Re*=-1;
      Gwi.V[i].Im*=-1;
    }
    for(i=0;i<A->BlockSize*A->BlockSize;i+=A->BlockSize+1)
    {
      Gwi.V[i].Re+=omega;
      Gwi.V[i].Im+=Gamma/2.0;
    }
    CompactMatrixInverse(&Gwi);
    TransformCompactMatrixInlineNC(&Gwi, &b); /* Gwi <- b Gw b^dag */
    CompactMatrixAdd(Gw, &Gwi);
  }
  // add a1
  iblock = 1;
  if(A->QComplex)
  {
    a.V=&(A->a[iblock    *A->BlockSize*A->BlockSize]);
  }
  else
  {
    a.VRe=&(A->aRe[iblock    *A->BlockSize*A->BlockSize]);
  }
  CompactMatrixAdd(Gw, &a);
  // multipli by -1
  for(i=0;i<A->BlockSize*A->BlockSize;i++)
  {
    Gw->V[i].Re*=-1;
    Gw->V[i].Im*=-1;
  }
  // add w
  for(i=0;i<A->BlockSize*A->BlockSize;i+=A->BlockSize+1)
  {
    Gw->V[i].Re+=omega;
    Gw->V[i].Im+=Gamma/2.0;
  }
  // inverse
  CompactMatrixInverse(Gw);
  // multiply by b0
  if(A->QComplex)
  {
    b.V=A->b;
  }
  else
  {
    b.VRe=A->bRe;
  }
  TransformCompactMatrixInlinecT(Gw, &b); /* Gw <- b^* Gw b^T */
  // add a0
  if(A->QComplex)
  {
    a.V=A->a;
  }
  else
  {
    a.VRe=A->aRe;
  }
  CompactMatrixAdd(Gw, &a);
  return 0;

}

uint32_t EvaluateBlockAndersonMatrix(void const *A, CompactMatrixType * Gw, const double omega, const double Gamma)
{
	BlockAndersonMatrixType *ATypecast = (BlockAndersonMatrixType*) A;
	return EvaluateBlockAndersonMatrix(ATypecast, Gw, omega, Gamma);
}

/*  Transform a BlockTriDiagonalMatrix representation of G to a spectral representation of G at one specific energy.
 *
 *  On entry:
 *    A contains the BlockTriDiagonal matrix to transform
 *    G is a pointer to SpectraType wich allocats A->BlockSize^2 spectra
 *    Emin and Emax set the min and max energy in the spectrum
 *    N the number of energies.
 *    Gamma sets the lorenzian broadening Gamma > 3*(Emax-Emin)/N
 *
 *    Gw = a0 + b0 (1/( w + i Gamma/2 - a1 - b1  ( w + i Gamma/2 - a2 - b2 ( ... ) b2^+ ) b1^+ )) b0^+
 *
 *    Gw is not initialized on entry but is when function finishes
 *
 */
uint32_t EvaluateBlockTriDiagonalMatrix(BlockTriDiagonalMatrixType const * A, CompactMatrixType * Gw, const double omega, const double Gamma)
{
  int i,iblock;
  CompactMatrixType a,b;
  a.Ni=A->BlockSize;
  a.Nj=A->BlockSize;
  b.Ni=A->BlockSize;
  b.Nj=A->BlockSize;
  a.QComplex=A->QComplex;
  b.QComplex=A->QComplex;
  a.a0=0;
  a.b0=0;
  a.mu=0;
  b.a0=0;
  b.b0=0;
  b.mu=0;
  Gw->QComplex=1;
  Gw->Ni=A->BlockSize;
  Gw->Nj=A->BlockSize;
  Gw->a0=0;
  Gw->b0=0;
  Gw->mu=A->mu;
  InitCompactMatrix(Gw);
  for(i=0;i<A->BlockSize*A->BlockSize;i++)
  {
    Gw->V[i].Re=0;
    Gw->V[i].Im=0;
  }
  for(iblock=A->N;iblock>0;iblock--)
  {
    /* set a and b to point to right matrices */
    if(A->QComplex)
    {
      a.V=&(A->a[iblock    *A->BlockSize*A->BlockSize]);
      b.V=&(A->b[(iblock-1)*A->BlockSize*A->BlockSize]);
    }
    else
    {
      a.VRe=&(A->aRe[iblock    *A->BlockSize*A->BlockSize]);
      b.VRe=&(A->bRe[(iblock-1)*A->BlockSize*A->BlockSize]);
    }
    CompactMatrixAdd(Gw, &a);
    for(i=0;i<A->BlockSize*A->BlockSize;i++)
    {
      Gw->V[i].Re*=-1;
      Gw->V[i].Im*=-1;
    }
    for(i=0;i<A->BlockSize*A->BlockSize;i+=A->BlockSize+1)
    {
      Gw->V[i].Re+=omega;
      Gw->V[i].Im+=Gamma/2.0;
    }
    CompactMatrixInverse(Gw);
    if(iblock==1)
    {
      TransformCompactMatrixInlinecT(Gw, &b); /* Gw <- b^* Gw b^T */
    }
    else
    {
      TransformCompactMatrixInlineNC(Gw, &b); /* Gw <- b Gw b^dag */
    }
  }
  // add a0
  if(A->QComplex)
  {
    a.V=A->a;
  }
  else
  {
    a.VRe=A->aRe;
  }
  CompactMatrixAdd(Gw, &a);
  return 0;
}

uint32_t EvaluateBlockTriDiagonalMatrix(void const * A, CompactMatrixType * Gw, const double omega, const double Gamma)
{
	BlockTriDiagonalMatrixType *ATypecast = (BlockTriDiagonalMatrixType*) A;
	return EvaluateBlockTriDiagonalMatrix(ATypecast, Gw, omega, Gamma);
}

/*
 *
 */
uint32_t EvaluateBlockNaturalImpurityOrbital(BlockNaturalImpurityOrbitalMatrixType const *A, CompactMatrixType * Gw, const double omega, const double Gamma)
{
  //TODO: all
  printf("Not implemented EvaluateBlockaturalImpurityOrbital\n");
  return 0;
}

uint32_t EvaluateBlockNaturalImpurityOrbital(void const *A, CompactMatrixType * Gw, const double omega, const double Gamma)
{
	BlockNaturalImpurityOrbitalMatrixType *ATypecast = (BlockNaturalImpurityOrbitalMatrixType*) A;
	return EvaluateBlockNaturalImpurityOrbital(ATypecast, Gw, omega, Gamma);
}

uint32_t EvaluateBlockResponsefunctionType(ResponsefunctionType const *A, CompactMatrixType * Gw, const double omega, const double Gamma)
{
  switch (A->type)
  {
    case 't':
    case 'T':
    {
      return EvaluateBlockTriDiagonalMatrix(A->data, Gw, omega, Gamma);
    }
    case 'a':
    case 'A':
    {
      return EvaluateBlockAndersonMatrix(A->data, Gw, omega, Gamma);
    }
    case 'l':
    case 'L':
    {
      return EvaluateBlockListOfPoles(A->data, Gw, omega, Gamma);
    }
    case 'n':
    case 'N':
    {
      return EvaluateBlockNaturalImpurityOrbital(A->data, Gw, omega, Gamma);
    }
    default:
    {
      printf("Error in LuaPush_ResponseFunction: Type of Response function should be either T/t,A/a,L/l, or N/n and not %c\n",A->type);
      return 0;
    }
  }
  return 0;
}


/*  =========================================================================
 *  Block list of poles to __
 *  ========================================================================= */


/*  Transform a Block List of Poles representation of G to a spectral representation of G.
 *
 *  On entry:
 *    A contains the BlockLisOfPoles to transform
 *    G is a pointer to SpectraType wich is not allocated
 *    Emin and Emax set the min and max energy in the spectrum
 *    N the number of energies.
 *
 *    G = sum_i lim_{n->0^+} b_i^* 1/(w - a_i + I n) b_i^T
 *    G = sum_i lim_{n->0^+} b^2 / (w - a_i + I n)
 *
 *    Implementation first the im part of lim_{n->0^+} 1/(w - a_i + I n) is added as a "delta" function moved to the energies E_i and E_{i+1} such that
 *      E_i < a and E_{i+1} > a
 *
 *      after this is done for all poles the real part of lim_{n->0^+} 1/(w - a_i + I n) is added.
 */
uint32_t BlockListOfPolesToSpectra(BlockListOfPoles const A, SpectraType *G, double Emin, double Emax, uint32_t N)
{
  int i,j,k;
  uint32_t iEmin,iEplus;
  double wmin,wplus;
  SpectraType *Gtemp;
  G[0].N=N;
  G[0].Emin=Emin;
  G[0].Emax=Emax;
  G[0].EminPole=Emin;
  G[0].EmaxPole=Emax;
  G[0].dE=(Emax-Emin)/((double)N);
  G[0].Gamma=0;
  G[0].ResonantEnergy=NAN;
  strcpy(G[0].Name, A.Name);
  for(i=0;i<A.BlockSize*A.BlockSize;i++)
  {
    SpectrumAllocate(G[0], &(G[i]));
  }
  if(A.QComplex)
  {
    Gtemp= (SpectraType*)malloc(A.BlockSize*A.BlockSize*sizeof(SpectraType));
    for(i=0;i<A.BlockSize*A.BlockSize;i++)
    {
      SpectrumAllocate(G[0], &(Gtemp[i]));
    }
    /* j loops over all blocks in A */
    for(j=0;j<(A.BlockSize*A.BlockSize+1)*A.N;j+=A.BlockSize*A.BlockSize+1)
    {
      /* the poles have energy A.V[j].Re */
      if(A.V[j].Re<=Emin)
      {
        iEmin=0;
        iEplus=0;
        wmin=1;
        wplus=0;
      }
      else if(A.V[j].Re>=Emax)
      {
        iEmin=N;
        iEplus=N;
        wmin=1;
        wplus=0;
      }
      else
      {
        iEmin=(uint32_t)((A.V[j].Re-Emin)/G[0].dE);
        iEplus=iEmin+1;
        wplus=((A.V[j].Re-Emin)/G[0].dE)-iEmin;
        wmin=1-wplus;
      }
      /* add the im part of the lorenzian to energy iEmin and iEplus */
      for(i=0;i<A.BlockSize*A.BlockSize;i++)
      {
        Gtemp[i].Im[iEmin]-=wmin*A.V[j+i+1].Re;
        Gtemp[i].Re[iEmin]+=wmin*A.V[j+i+1].Im;
        Gtemp[i].Im[iEplus]-=wplus*A.V[j+i+1].Re;
        Gtemp[i].Re[iEplus]+=wplus*A.V[j+i+1].Im;
      }
    }
  }
  else
  {
    Gtemp=G;
    for(j=0;j<(A.BlockSize*A.BlockSize+1)*A.N;j+=A.BlockSize*A.BlockSize+1)
    {
      /* the poles have energy A.V[j].Re */
      if(A.VRe[j]<=Emin)
      {
        iEmin=0;
        iEplus=0;
        wmin=1;
        wplus=0;
      }
      else if(A.VRe[j]>=Emax)
      {
        iEmin=N;
        iEplus=N;
        wmin=1;
        wplus=0;
      }
      else
      {
        iEmin=(uint32_t)((A.VRe[j]-Emin)/G[0].dE);
        iEplus=iEmin+1;
        wplus=((A.VRe[j]-Emin)/G[0].dE)-iEmin;
        wmin=1-wplus;
      }
      /* add the im part of the lorenzian to energy iEmin and iEplus */
      for(i=0;i<A.BlockSize*A.BlockSize;i++)
      {
        Gtemp[i].Im[iEmin]-=wmin*A.VRe[j+i+1];
        Gtemp[i].Im[iEplus]-=wplus*A.VRe[j+i+1];
      }
    }
  }
  /* the im part of the lorenzians is there, we now need to add the real part */
  /* note that in the complex case some lorenzians have a complex pre-factor */
  /* weight w_i at energy E_i contributes w_i / (w - E_i) to all energies w, except at w=E_i */
  for(i=0;i<A.BlockSize*A.BlockSize;i++)
  {
    for(j=0;j<=N;j++)
    {
      for(k=0;k<j;k++)
      {
        G[i].Re[k]-=Gtemp[i].Im[j]/(Pi*(k-j));
      }
      for(k=j+1;k<=N;k++)
      {
        G[i].Re[k]-=Gtemp[i].Im[j]/(Pi*(k-j));
      }
    }
  }
  if(A.QComplex)
  {
    for(i=0;i<A.BlockSize*A.BlockSize;i++)
    {
      for(j=0;j<=N;j++)
      {
        for(k=0;k<j;k++)
        {
          G[i].Im[k]+=Gtemp[i].Re[j]/(Pi*(k-j));
        }
        for(k=j+1;k<=N;k++)
        {
          G[i].Im[k]+=Gtemp[i].Re[j]/(Pi*(k-j));
        }
      }
      for(j=0;j<=N;j++)
      {
        G[i].Re[j]+=Gtemp[i].Re[j];
        G[i].Im[j]+=Gtemp[i].Im[j];
      }
      SpectrumFree(&(Gtemp[i]));
    }
    free(Gtemp);
  }
  for(i=0;i<A.BlockSize*A.BlockSize;i++)
  {
    for(j=0;j<=N;j++)
    {
      G[i].Im[j]*=Pi/G[i].dE;
      G[i].Re[j]*=Pi/G[i].dE;
    }
  }
  /* add a0 */
  if(A.QComplex)
  {
    for(i=0;i<A.BlockSize*A.BlockSize;i++)
    {
      for(j=0;j<=N;j++)
      {
        G[i].Re[j]+=A.a0[i].Re;
        G[i].Im[j]+=A.a0[i].Im;
      }
    }
  }
  else
  {
    for(i=0;i<A.BlockSize*A.BlockSize;i++)
    {
      for(j=0;j<=N;j++)
      {
        G[i].Re[j]+=A.a0Re[i];
      }
    }
  }
  return 0;
}
uint32_t BlockListOfPolesToSpectraWithBroadening(BlockListOfPoles const A, SpectraType *G, double Emin, double Emax, uint32_t N, double Gamma)
{
  int iE;
  uint32_t i,iblock,itot;
  Complex16 Lorentzian;
  double omega;
  G[0].N=N;
  G[0].Emin=Emin;
  G[0].Emax=Emax;
  G[0].EminPole=Emin;
  G[0].EmaxPole=Emax;
  G[0].dE=(Emax-Emin)/((double)N);
  G[0].Gamma=Gamma;
  G[0].ResonantEnergy=NAN;
  strcpy(G[0].Name, A.Name);
  for(i=0;i<A.BlockSize*A.BlockSize;i++)
  {
    SpectrumAllocate(G[0], &(G[i]));
  }
  if(A.QComplex==0)
  {// real
#pragma omp parallel for private(omega,i,itot,iblock,Lorentzian)
    for(iE=0;iE<=N;iE++)
    {
      omega = ((double)iE/(double)N)*(Emax-Emin)+Emin;
      // add A0
      for(i=0;i<A.BlockSize*A.BlockSize;i++)
      {
        G[i].Re[iE]=A.VRe[A.N*(1+A.BlockSize*A.BlockSize)+i];
        G[i].Im[iE]=0;
      }
      // loop over all poles
      itot=0;
      for(iblock=0;iblock<A.N;iblock++)
      {
        Lorentzian.Re = (omega-A.VRe[itot]) / (Sqr(Gamma/2.0) + Sqr(omega-A.VRe[itot]));
        Lorentzian.Im = (-Gamma/2.0       ) / (Sqr(Gamma/2.0) + Sqr(omega-A.VRe[itot]));
        itot++;
        for(i=0;i<A.BlockSize*A.BlockSize;i++)
        {
          G[i].Re[iE] += A.VRe[itot] * Lorentzian.Re;
          G[i].Im[iE] += A.VRe[itot] * Lorentzian.Im;
          itot++;
        }
      }
    }
  }
  else
  {// complex
#pragma omp parallel for private(omega,i,itot,iblock,Lorentzian)
    for(iE=0;iE<=N;iE++)
    {
      omega = ((double)iE/(double)N)*(Emax-Emin)+Emin;
      // add A0
      for(i=0;i<A.BlockSize*A.BlockSize;i++)
      {
        G[i].Re[iE]=A.V[A.N*(1+A.BlockSize*A.BlockSize)+i].Re;
        G[i].Im[iE]=A.V[A.N*(1+A.BlockSize*A.BlockSize)+i].Im;
      }
      // loop over all poles
      itot=0;
      for(iblock=0;iblock<A.N;iblock++)
      {
        Lorentzian.Re = (omega-A.V[itot].Re) / (Sqr(Gamma/2.0) + Sqr(omega-A.V[itot].Re));
        Lorentzian.Im = (-Gamma/2.         ) / (Sqr(Gamma/2.0) + Sqr(omega-A.V[itot].Re));
        itot++;
        for(i=0;i<A.BlockSize*A.BlockSize;i++)
        {
          G[i].Re[iE] += A.V[itot].Re * Lorentzian.Re - A.V[itot].Im * Lorentzian.Im;
          G[i].Im[iE] += A.V[itot].Re * Lorentzian.Im + A.V[itot].Im * Lorentzian.Re;
          itot++;
        }
      }
    }
  }
  return 0;
}


uint32_t BlockListOfPolesEigenValues(const BlockListOfPoles * const A, double * val)
{
  //TODO all
  return 0;
}
uint32_t BlockListOfPolesEigenSystem(const BlockListOfPoles * const A, double * val, CompactMatrixType * fun)
{
  //TODO all
  return 0;
}




/*  on entry A contains the Green's function in a BlockListOfPoles format.
 *  on excit B contains the Green's function in a BlockTridiagonal format. (note that a0 and b0 are aditional prefactors)
 *  The returned matrix is always of the dimension NBLockMax by A.BlockSize, if couplings are small (SingularValue) they are set to zero.
 */
uint32_t BlockListOfPolesToBlockTridiagonalMatrix(BlockListOfPoles const A, BlockTriDiagonalMatrixType *const B, int const NBlockMax, double const SingularValue)
{
  double *alpha;
  int k,i,ialpha,iN;
  char name[]="temp";
  CompactMatrixType betaisqr,betai,beta,T;
  
  betaisqr.QComplex=A.QComplex;
  betaisqr.Ni=A.BlockSize;
  betaisqr.Nj=A.BlockSize;
  betaisqr.Name=name;
  
  beta.QComplex=A.QComplex;
  beta.Ni=A.BlockSize*A.N;
  beta.Nj=A.BlockSize;
  InitCompactMatrix(&beta);
  
  alpha= (double*)malloc(A.N*A.BlockSize*sizeof(double));
  ialpha=0;
  iN=0;
  /* create the starting vectors beta as the square root of the spectral weight (DOS) in list of poles */
  if(A.QComplex)
  {
    for(k=0;k<A.N;k++)
    {
      betaisqr.V = &(A.V[iN+1]);
      CompactMatrixSqrt(&betaisqr, &betai);
      memcpy(&(beta.V[iN-k]), betai.V, sizeof(Complex16)*A.BlockSize*A.BlockSize);
      FreeCompactMatrix(&betai);
      for(i=0;i<A.BlockSize;i++)
      {
        alpha[ialpha]=A.V[iN].Re;
        ialpha++;
      }
      iN+=1+A.BlockSize*A.BlockSize;
    }
  }
  else /* A is real */
  {
    for(k=0;k<A.N;k++)
    {
      betaisqr.VRe = &(A.VRe[iN+1]);
      CompactMatrixSqrt(&betaisqr, &betai);
      memcpy(&(beta.VRe[iN-k]), betai.VRe, sizeof(double)*A.BlockSize*A.BlockSize);
      FreeCompactMatrix(&betai);
      for(i=0;i<A.BlockSize;i++)
      {
        alpha[ialpha]=A.VRe[iN];
        ialpha++;
      }
      iN+=1+A.BlockSize*A.BlockSize;
    }
  }
  CompactMatrixConjugateTranspose(&beta);
  DiagonalMatrixBlockTridiagonalize(alpha, &beta, &T, B, NBlockMax, SingularValue);
  B->mu=A.mu;
  free(alpha);
  FreeCompactMatrix(&beta);
  if(B->QComplex)
  {
    for(i=0;i<A.BlockSize*A.BlockSize;i++)
    {
      B->a[i]=A.a0[i];
    }
  }
  else
  {
    for(i=0;i<A.BlockSize*A.BlockSize;i++)
    {
      B->aRe[i]=A.a0Re[i];
    }
  }
  return 0;
}

/*  on entry A contains the Green's function in a BlockListOfPoles format.
 *  on excit B contains the Green's function in a BlockTridiagonal format. (note that a0 and b0 are aditional prefactors)
 *  The returned matrix is always of the dimension NBLockMax by A.BlockSize, if couplings are small (SingularValue) they are set to zero.
 *
 *  AT is a matrix of dimension Ni times Nj=A.BlockSize.
 *  The actual tridiagonalization is done for the starting vectors AT * beta
 *
 */
uint32_t BlockListOfPolesToBlockTridiagonalMatrixReduceDimension(BlockListOfPoles const A, BlockTriDiagonalMatrixType *B, CompactMatrixType const *AT,  int const NBlockMax, double const SingularValue)
{
  double *alpha;
  int k,i,ialpha,iN;
  char name[]="temp";
  CompactMatrixType betaisqr,betai,beta,T,ATbeta;
  
  betaisqr.QComplex=A.QComplex;
  betaisqr.Ni=A.BlockSize;
  betaisqr.Nj=A.BlockSize;
  betaisqr.Name=name;
  
  beta.QComplex=A.QComplex;
  beta.Ni=A.BlockSize*A.N;
  beta.Nj=A.BlockSize;
  InitCompactMatrix(&beta);
  
  alpha= (double*)malloc(A.N*A.BlockSize*sizeof(double));
  ialpha=0;
  iN=0;
  if(A.QComplex)
  {
    for(k=0;k<A.N;k++)
    {
      betaisqr.V = &(A.V[iN+1]);
      CompactMatrixSqrt(&betaisqr, &betai);
      memcpy(&(beta.V[iN-k]), betai.V, sizeof(Complex16)*A.BlockSize*A.BlockSize);
      FreeCompactMatrix(&betai);
      for(i=0;i<A.BlockSize;i++)
      {
        alpha[ialpha]=A.V[iN].Re;
        ialpha++;
      }
      iN+=1+A.BlockSize*A.BlockSize;
    }
  }
  else /* A is real */
  {
    for(k=0;k<A.N;k++)
    {
      betaisqr.VRe = &(A.VRe[iN+1]);
      CompactMatrixSqrt(&betaisqr, &betai);
      memcpy(&(beta.VRe[iN-k]), betai.VRe, sizeof(double)*A.BlockSize*A.BlockSize);
      FreeCompactMatrix(&betai);
      for(i=0;i<A.BlockSize;i++)
      {
        alpha[ialpha]=A.VRe[iN];
        ialpha++;
      }
      iN+=1+A.BlockSize*A.BlockSize;
    }
  }
  // note that we need to complex transpose beta to get a set of NBlock row vectors
  CompactMatrixMultiply('N', 'C', AT, &beta, &ATbeta);
  FreeCompactMatrix(&beta);
  DiagonalMatrixBlockTridiagonalize(alpha, &ATbeta, &T, B, NBlockMax, SingularValue);
  B->mu=A.mu;
  free(alpha);
  FreeCompactMatrix(&ATbeta);
  return 0;
}


/*  on entry A contains the Green's function in a BlockListOfPoles format.
 *  on exit B contains the Green's function in a NaturalImpurityOrbital format.
 *
 */
uint32_t BlockListOfPolesToBlockNaturalImpurityOrbitalMatrix(BlockListOfPoles const *A, BlockNaturalImpurityOrbitalMatrixType * const B, double const kbT, uint32_t const NBlockMaxVal, uint32_t const NBlockMaxCon)
{
  BlockListOfPoles Acon, Aval;
  B->mu=A->mu;
  B->QComplex=A->QComplex;
  B->BlockSize = A->BlockSize;
  strcpy(B->Name, A->Name);
  B->a0.QComplex=A->QComplex;
  B->a0.Ni=A->BlockSize;
  B->a0.Nj=A->BlockSize;
  InitCompactMatrix(&B->a0);
  if(A->QComplex)
  {
    memcpy(B->a0.V, A->a0, A->BlockSize*A->BlockSize*sizeof(Complex16));
  }
  else
  {
    memcpy(B->a0.VRe, A->a0Re, A->BlockSize*A->BlockSize*sizeof(double));
  }
  CopyBlockListOfPolesCon(*A, &Acon, A->mu, kbT);
  CopyBlockListOfPolesVal(*A, &Aval, A->mu, kbT);
  BlockListOfPolesToBlockTridiagonalMatrix(Acon, &(B->con), NBlockMaxCon, 1000*DBL_EPSILON);
  BlockListOfPolesToBlockTridiagonalMatrix(Aval, &(B->val), NBlockMaxVal, 1000*DBL_EPSILON);
  SetTMTFromValandConInBlockNaturalImpurityOrbitalMatrixType(B);
  return 0;
}

/*  Input:
 *    A         BlockListOfPoles to converge to Nat. Imp. orbitlals
 *    AT        Matrix defining the reduction of the impurity size
 *              AT->Nj original dimension of A
 *              AT->Ni new dimension
 *    kBT       Effective artificial temperature used to determine the density matrix
 *    NBlockVal Lenght of chain in result
 *    NBlockCon Lenght of chain in result
 *  Output
 *    B resutling G on Nat. Imp. Orb. basis
 */
uint32_t BlockListOfPolesToBlockNaturalImpurityOrbitalMatrixReduceDimension(BlockListOfPoles const *A, BlockNaturalImpurityOrbitalMatrixType * const B, CompactMatrixType const *AT, double const kbT, uint32_t const NBlockMaxVal, uint32_t const NBlockMaxCon)
{
  CompactMatrixType Aa0;
  BlockListOfPoles Acon, Aval;
  B->mu=A->mu;
  B->QComplex=A->QComplex;
  B->BlockSize = AT->Ni;
  strcpy(B->Name, A->Name);
  // define a0 as a compact matrix such that I can call transform
  Aa0.QComplex=A->QComplex;
  Aa0.Ni=A->BlockSize;
  Aa0.Nj=A->BlockSize;
  Aa0.a0=0;
  Aa0.b0=1;
  if(Aa0.QComplex)
  {
    Aa0.V=A->a0;
  }
  else
  {
    Aa0.VRe=A->a0Re;
  }
  TransformCompactMatrixNC(&Aa0, AT, &(B->a0));
  CopyBlockListOfPolesCon(*A, &Acon, A->mu, kbT);
  CopyBlockListOfPolesVal(*A, &Aval, A->mu, kbT);
  BlockListOfPolesToBlockTridiagonalMatrixReduceDimension(Acon, &(B->con), AT, NBlockMaxCon+1, 1000*DBL_EPSILON);
  BlockListOfPolesToBlockTridiagonalMatrixReduceDimension(Aval, &(B->val), AT, NBlockMaxVal+1, 1000*DBL_EPSILON);
  SetTMTFromValandConInBlockNaturalImpurityOrbitalMatrixType(B);
  return 0;
}


/*
 *
 * Converts a BlockListofpoles to blockandersonmatrix, using a blocktridiagonalmatrix as intermediate
 *
 *
 */
uint32_t BlockListOfPolesToBlockAndersonMatrix(BlockListOfPoles const A, BlockAndersonMatrixType *const B, int const NBlockMax)
{
  BlockTriDiagonalMatrixType ATri;
  if(BlockListOfPolesToBlockTridiagonalMatrix(A, &ATri, NBlockMax, DBL_EPSILON))
  {
    printf("BlockListOfPolesToBlockTridiagonalMatrix failed in BlockListOfPolesToBlockAndersonMatrix\n");
    return 1;
  }
  BlockTridiagonalMatrixChop(&ATri,10*DBL_EPSILON);
  if(BlockTriDiagonalMatrixToBlockAndersonMatrix(ATri,B))
  {
    printf("BlockTriDiagonalMatrixToBlockAndersonMatrix failed in BlockListOfPolesToBlockAndersonMatrix\n");
    return 1;
  }
  FreeBlockTridiagonalMatrix(&ATri);
  return 0;
}

/*  on entry A contains the Green's function in a BlockListOfPoles format.
 *
 * Converts a BlockListofpoles to blockandersonmatrix, using a blocktridiagonalmatrix as intermediate
 *
 *
 *  AT specifies how to reduce the dimension
 *
 */
uint32_t BlockListOfPolesToBlockAndersonMatrixReduceDimension(BlockListOfPoles const A, BlockAndersonMatrixType *B, CompactMatrixType const *AT, int const NBlockMax)
{
  BlockTriDiagonalMatrixType ATri;
  if(BlockListOfPolesToBlockTridiagonalMatrixReduceDimension(A, &ATri, AT, NBlockMax, DBL_EPSILON))
  {
    printf("BlockListOfPolesToBlockTridiagonalMatrixReduceDimension failed in BlockListOfPolesToBlockAndersonMatrixReduceDimension\n");
    return 1;
  }
  BlockTridiagonalMatrixChop(&ATri,100*DBL_EPSILON);
  if(BlockTriDiagonalMatrixToBlockAndersonMatrix(ATri,B))
  {
    printf("BlockTriDiagonalMatrixToBlockAndersonMatrix failed in BlockListOfPolesToBlockAndersonMatrixReduceDimension\n");
    return 1;
  }
  FreeBlockTridiagonalMatrix(&ATri);
  return 0;
}



/*  =========================================================================
 *  BlockAndersonMatrixType to __
 *  ========================================================================= */



/*
 *
 */
uint32_t BlockAndersonMatrixToSpectra(BlockAndersonMatrixType const A, SpectraType *G, double Emin, double Emax, uint32_t N)
{
  //TODO: all
  printf("Not implemented BlockAndersonMatrixToSpectra\n");
  return 1;
}


uint32_t BlockAndersonMatrixEigenValues(const BlockAndersonMatrixType * const A, double * val)
{
  //TODO all
  return 0;
}
uint32_t BlockAndersonMatrixEigenSystem(const BlockAndersonMatrixType * const A, double * val, CompactMatrixType * fun)
{
  //TODO all
  return 0;
}



uint32_t BlockAndersonMatrixToBlockListOfPoles(BlockAndersonMatrixType const AAnd, BlockListOfPoles * const B)
{
  //TODO: all
  printf("Not implemented BlockAndersonMatrixToBlockListOfPoles\n");
  return 1;
}
uint32_t BlockAndersonMatrixToBlockListOfPolesVariableGrid(BlockAndersonMatrixType const AAnd, BlockListOfPoles * const B, const double * const Grid, const uint32_t NGridPoints)
{
  //TODO: all
  printf("Not implemented BlockAndersonMatrixToBlockListOfPoles\n");
  return 1;
}
uint32_t BlockAndersonMatrixToBlockTridiagonalMatrix(BlockAndersonMatrixType const AAnd, BlockTriDiagonalMatrixType * const B)
{
  //TODO: all
  printf("Not implemented BlockAndersonMatrixToBlockTridiagonalMatrix\n");
  return 1;
}
uint32_t BlockAndersonMatrixToBlockNaturalImpurityOrbitalMatrix(BlockAndersonMatrixType const AAnd, BlockNaturalImpurityOrbitalMatrixType * const B)
{
  //TODO: all
  printf("Not implemented BlockAndersonMatrixToBlockNaturalImpurityOrbitalMatrix\n");
  return 1;
}



/*  =========================================================================
 *  BlockTriDiagonalMatrixType to __
 *  ========================================================================= */



/*  Transform a BlockTriDiagonalMatrix representation of G to a spectral representation of G.
 *
 *  On entry:
 *    A contains the BlockTriDiagonal matrix to transform
 *    G is a pointer to SpectraType wich allocats A->BlockSize^2 spectra
 *    Emin and Emax set the min and max energy in the spectrum
 *    N the number of energies.
 *    Gamma sets the lorenzian broadening Gamma > 3*(Emax-Emin)/N
 *
 *    G = a0 + b0 (1/( w + i Gamma/2 - a1 - b1  ( w + i Gamma/2 - a2 - b2 ( ... ) b2^+ ) b1^+ )) b0^+
 *
 */
uint32_t BlockTriDiagonalMatrixToSpectra(const double Energy, BlockTriDiagonalMatrixType const * A, SpectraType *const G, double const Emin, double const Emax, uint32_t const N, double const Gamma)
{
  int i,iw,iblock;
  double omega;
  CompactMatrixType Gw,a,b;
  G[0].N=N;
  G[0].Emin=Emin;
  G[0].Emax=Emax;
  G[0].EminPole=Emin;
  G[0].EmaxPole=Emax;
  G[0].dE=(Emax-Emin)/((double)N);
  G[0].Gamma=Gamma;
  G[0].ResonantEnergy=NAN;
  strcpy(G[0].Name, A->Name);
  for(i=0;i<A->BlockSize*A->BlockSize;i++)
  {
    SpectrumAllocate(G[0], &(G[i]));
  }
  a.Ni=A->BlockSize;
  a.Nj=A->BlockSize;
  b.Ni=A->BlockSize;
  b.Nj=A->BlockSize;
  a.QComplex=A->QComplex;
  b.QComplex=A->QComplex;
  Gw.QComplex=1;
  Gw.Ni=A->BlockSize;
  Gw.Nj=A->BlockSize;
  InitCompactMatrix(&Gw);
  for(iw=0;iw<=N;iw++)
  {
    omega = Emin + iw * G[0].dE + Energy;
    for(i=0;i<A->BlockSize*A->BlockSize;i++)
    {
      Gw.V[i].Re=0;
      Gw.V[i].Im=0;
    }
    for(iblock=A->N;iblock>0;iblock--)
    {
      /* set a and b to point to right matrices */
      if(A->QComplex)
      {
        a.V=&(A->a[iblock    *A->BlockSize*A->BlockSize]);
        b.V=&(A->b[(iblock-1)*A->BlockSize*A->BlockSize]);
      }
      else
      {
        a.VRe=&(A->aRe[iblock    *A->BlockSize*A->BlockSize]);
        b.VRe=&(A->bRe[(iblock-1)*A->BlockSize*A->BlockSize]);
      }
      CompactMatrixAdd(&Gw, &a);
      for(i=0;i<A->BlockSize*A->BlockSize;i++)
      {
        Gw.V[i].Re*=-1;
        Gw.V[i].Im*=-1;
      }
      for(i=0;i<A->BlockSize*A->BlockSize;i+=A->BlockSize+1)
      {
        Gw.V[i].Re+=omega;
        Gw.V[i].Im+=Gamma/2.0;
      }
      CompactMatrixInverse(&Gw);
      if(iblock==1)
      {
        TransformCompactMatrixInlinecT(&Gw, &b); /* Gw <- b Gw b^T */
      }
      else
      {
        TransformCompactMatrixInlineNC(&Gw, &b); /* Gw <- b Gw b^dag */
      }
    }
    // add a0
    if(A->QComplex)
    {
      a.V=A->a;
    }
    else
    {
      a.VRe=A->aRe;
    }
    CompactMatrixAdd(&Gw, &a);
    // add Gw to G
    for(i=0;i<A->BlockSize*A->BlockSize;i++)
    {
      G[i].Re[iw]=Gw.V[i].Re;
      G[i].Im[iw]=Gw.V[i].Im;
    }
  }
  FreeCompactMatrix(&Gw);
  return 0;
}



// Given a matrix in BlockTriDiagonal form calculate the Eigenvalues
// val is allocated with a minimal length of A->N * A->BlockSize
uint32_t BlockTriDiagonalMatrixEigenValues(const BlockTriDiagonalMatrixType * const A, double * val)
{
  uint32_t i,j,k,iBlock,jBlock;
  char const *UPLO="U";
  char const *VECT="N";
  const uint32_t N=A->N*A->BlockSize;
  const uint32_t KD = A->BlockSize-1;
  const uint32_t LDAB = A->BlockSize;
  TriDiagonalMatrixType M;
  int INFO=0;
  M.N=N;
  InitTriDiagonalMatrix(&M);
  if(A->QComplex==0)
  {
    double * AB;
    double *WORK;
    AB= (double*)calloc(LDAB*N,sizeof(double));
    if(AB==NULL)
    {
      printf("ERROR: malloc failed in BlockTriDiagonalMatrixEigenValues\n");
      return 1;
    }
    WORK= (double*)malloc((N+1)*sizeof(double));
    if(WORK==NULL)
    {
      free(AB);
      printf("ERROR: malloc failed in BlockTriDiagonalMatrixEigenValues\n");
      return 1;
    }
    // use dsbtrd to transform A into a tri diagonal form
    //fill matrix AB
    // if UPLO = 'U', AB(i,j-i) = A(i,j); (fortran and c column v.s. row major)
    j=0;
    for(k=0;k<A->N;k++)
    {
      i=0;
      for(iBlock=0;iBlock<A->BlockSize;iBlock++)
      {
        for(jBlock=iBlock;jBlock<A->BlockSize;jBlock++)
        {
          AB[j]=A->aRe[i];
          i++;
          j++;
        }
        j+=iBlock; // set the aditional values of the row to zero
      }
    }
    // Tridiagonalize the Block tridiagonal matrix
    dsbtrd(VECT, UPLO, &N, &KD, AB, &LDAB, &(M.a[1]), &(M.b[1]), NULL, &N, WORK, &INFO);
    free(AB);
    free(WORK);
    if(INFO!=0)
    {
      printf("ERROR: dsbtrd failed in BlockTriDiagonalMatrixEigenValues with INFO %i\n",INFO);
      return 1;
    }
  }
  else
  {
    Complex16 * AB;
    Complex16 *WORK;
    AB= (Complex16*)calloc(LDAB*N,sizeof(Complex16));	// changed by Martin: initalize to zero bugfix 4.6.2019
    if(AB==NULL)
    {
      printf("ERROR: malloc failed in BlockTriDiagonalMatrixEigenValues\n");
      return 1;
    }
    WORK= (Complex16*)malloc(N*sizeof(Complex16));
    if(WORK==NULL)
    {
      free(AB);
      printf("ERROR: malloc failed in BlockTriDiagonalMatrixEigenValues\n");
      return 1;
    }
    //fill matrix AB
    // if UPLO = 'U', AB(i,j-i) = A(i,j); (fortran and c column v.s. row major)
    j=0;
    for(k=0;k<A->N;k++)
    {
      i=0;
      for(iBlock=0;iBlock<A->BlockSize;iBlock++)
      {
        for(jBlock=iBlock;jBlock<A->BlockSize;jBlock++)
        {
          AB[j].Re=A->a[i].Re;
          AB[j].Im=A->a[i].Im;
          i++;
          j++;
        }
        j+=iBlock; // set the aditional values of the row to zero
      }
    }
    // Tridiagonalize the Block tridiagonal matrix
    zhbtrd(VECT, UPLO, &N, &KD, AB, &LDAB, &(M.a[1]), &(M.b[1]), NULL, &N, WORK, &INFO);
    free(AB);
    free(WORK);
    if(INFO!=0)
    {
      printf("ERROR: zhbtrd failed in BlockTriDiagonalMatrixEigenValues with INFO %i\n",INFO);
      return 1;
    }
  }
  TriDiagonalMatrixToEigenvalues(M, val);
  FreeTriDiagonalMatrix(&M);
  return 0;
}

uint32_t BlockTriDiagonalMatrixEigenSystem(const BlockTriDiagonalMatrixType * const A, double * val, CompactMatrixType * fun)
{
  uint32_t j,k,iBlock,jBlock;
  char const *UPLO="L";
  //TODO: rewrite for UPLO=U such that we can use A*B^T aftwerwards (faster) Needs new filling of matrix AB
  char const *VECT="V";
  CompactMatrixType Q, funQ;
  const uint32_t N=A->N*A->BlockSize;
  const uint32_t KD = 2*A->BlockSize-1;
  const uint32_t LDAB = 2*A->BlockSize;
  TriDiagonalMatrixType M;
  int INFO=0;
  M.N=N;
  InitTriDiagonalMatrix(&M);
  Q.Ni=N;
  Q.Nj=N;
  Q.QComplex=A->QComplex;
  if(InitCompactMatrix(&Q))
  {
    printf("ERROR: InitCompactMatrix failed in BlockTriDiagonalMatrixEigenSystem\n");
    return 1;
  }
  if(A->QComplex==0)
  {
    double * AB;
    double *WORK;
    AB= (double*)calloc(LDAB*N,sizeof(double));
    if(AB==NULL)
    {
      printf("ERROR: malloc failed in BlockTriDiagonalMatrixEigenSystem\n");
      return 1;
    }
    WORK= (double*)malloc(N*(N+1)*sizeof(double));
    if(WORK==NULL)
    {
      free(AB);
      printf("ERROR: malloc failed in BlockTriDiagonalMatrixEigenSystem\n");
      return 1;
    }
    // use dsbtrd to transform A into a tri diagonal form
    //fill matrix AB
    j=0;
    for(k=1;k<A->N;k++)
    {
      for(iBlock=0;iBlock<A->BlockSize;iBlock++)
      {
        for(jBlock=iBlock;jBlock<A->BlockSize;jBlock++)
        {
          AB[j]=A->aRe[k*A->BlockSize*A->BlockSize + iBlock*A->BlockSize + jBlock];
          j++;
        }
        for(jBlock=0;jBlock<A->BlockSize;jBlock++)
        {
          AB[j]=A->bRe[k*A->BlockSize*A->BlockSize + iBlock*A->BlockSize + jBlock];
          j++;
        }
        j+=iBlock;
      }
    }
    k=A->N;
    for(iBlock=0;iBlock<A->BlockSize;iBlock++)
    {
      for(jBlock=iBlock;jBlock<A->BlockSize;jBlock++)
      {
        AB[j]=A->aRe[k*A->BlockSize*A->BlockSize + iBlock*A->BlockSize + jBlock];
        j++;
      }
      j+=A->BlockSize+iBlock;
    }
    // Tridiagonalize the Block tridiagonal matrix
    dsbtrd(VECT, UPLO, &N, &KD, AB, &LDAB, &(M.a[1]), &(M.b[1]), Q.VRe, &N, WORK, &INFO);
    free(AB);
    free(WORK);
    if(INFO!=0)
    {
      printf("ERROR: dsbtrd failed in BlockTriDiagonalMatrixEigenSystem with INFO %i\n",INFO);
      return 1;
    }
  }
  else
  {
    Complex16 * AB;
    Complex16 *WORK;
    AB= (Complex16*)calloc(LDAB*N,sizeof(Complex16));
    if(AB==NULL)
    {
      printf("ERROR: malloc failed in BlockTriDiagonalMatrixEigenSystem\n");
      return 1;
    }
    WORK= (Complex16*)malloc(N*sizeof(Complex16));
    if(WORK==NULL)
    {
      free(AB);
      printf("ERROR: malloc failed in BlockTriDiagonalMatrixEigenSystem\n");
      return 1;
    }
    //fill matrix AB
    j=0;
    for(k=1;k<A->N;k++)
    {
      for(iBlock=0;iBlock<A->BlockSize;iBlock++)
      {
        for(jBlock=iBlock;jBlock<A->BlockSize;jBlock++)
        {
          AB[j]=A->a[k*A->BlockSize*A->BlockSize + iBlock*A->BlockSize + jBlock];
          j++;
        }
        for(jBlock=0;jBlock<A->BlockSize;jBlock++)
        {
          AB[j]=A->b[k*A->BlockSize*A->BlockSize + iBlock*A->BlockSize + jBlock];
          j++;
        }
        j+=iBlock;
      }
    }
    k=A->N;
    for(iBlock=0;iBlock<A->BlockSize;iBlock++)
    {
      for(jBlock=iBlock;jBlock<A->BlockSize;jBlock++)
      {
        AB[j]=A->a[k*A->BlockSize*A->BlockSize + iBlock*A->BlockSize + jBlock];
        j++;
      }
      j+=A->BlockSize+iBlock;
    }
    // Tridiagonalize the Block tridiagonal matrix
    zhbtrd(VECT, UPLO, &N, &KD, AB, &LDAB, &(M.a[1]), &(M.b[1]), Q.V, &N, WORK, &INFO);
    free(AB);
    free(WORK);
    if(INFO!=0)
    {
      printf("ERROR: zhbtrd failed in BlockTriDiagonalMatrixEigenSystem with INFO %i\n",INFO);
      return 1;
    }
  }
  
  
  
  
  TriDiagonalMatrixToEigensystemCompact(&M, val, &funQ);
  CompactMatrixMultiply('N', 'N', &funQ, &Q, fun);
  CompactMatrixConjugate(fun);
  FreeCompactMatrix(&Q);
  FreeCompactMatrix(&funQ);
  FreeTriDiagonalMatrix(&M);
  return 0;
}

uint32_t BlockTriDiagonalMatrixEigenSystemWithDeflation(const BlockTriDiagonalMatrixType * const A, double * val, CompactMatrixType * fun, const uint32_t * const BlockSize)
{
  uint32_t j,k,iBlock,jBlock;
  char const *UPLO="L";
  //TODO: rewrite for UPLO=U such that we can use A*B^T aftwerwards (faster) Needs new filling of matrix AB
  char const *VECT="V";
  CompactMatrixType Q, funQ;
  uint32_t N=0;
  const uint32_t KD = 2*A->BlockSize-1;
  const uint32_t LDAB = 2*A->BlockSize;
  TriDiagonalMatrixType M;
  int INFO=0;
  for(j=1;j<=A->N;j++)
  {
    N+=BlockSize[j];
  }
  M.N=N;
  InitTriDiagonalMatrix(&M);
  Q.Ni=N;
  Q.Nj=N;
  Q.QComplex=A->QComplex;
  if(InitCompactMatrix(&Q))
  {
    printf("ERROR: InitCompactMatrix failed in BlockTriDiagonalMatrixEigenSystem\n");
    return 1;
  }
  if(A->QComplex==0)
  {
    double * AB;
    double *WORK;
    AB= (double*)calloc(LDAB*N,sizeof(double));
    if(AB==NULL)
    {
      printf("ERROR: malloc failed in BlockTriDiagonalMatrixEigenSystem\n");
      return 1;
    }
    WORK= (double*)malloc(N*(N+1)*sizeof(double));
    if(WORK==NULL)
    {
      free(AB);
      printf("ERROR: malloc failed in BlockTriDiagonalMatrixEigenSystem\n");
      return 1;
    }
    // use dsbtrd to transform A into a tri diagonal form
    //fill matrix AB
    j=0;
    for(k=1;k<A->N;k++)
    {
      for(iBlock=0;iBlock<BlockSize[k];iBlock++)
      {
        for(jBlock=iBlock;jBlock<BlockSize[k];jBlock++)
        {
          AB[j]=A->aRe[k*A->BlockSize*A->BlockSize + iBlock*A->BlockSize + jBlock];
          j++;
        }
        j+=A->BlockSize-BlockSize[k];
        for(jBlock=0;jBlock<BlockSize[k+1];jBlock++)
        {
          AB[j]=A->bRe[k*A->BlockSize*A->BlockSize + iBlock*A->BlockSize + jBlock];
          j++;
        }
        j+=A->BlockSize-BlockSize[k+1];
        j+=iBlock;
      }
    }
    k=A->N;
    for(iBlock=0;iBlock<BlockSize[A->N];iBlock++)
    {
      for(jBlock=iBlock;jBlock<BlockSize[A->N];jBlock++)
      {
        AB[j]=A->aRe[k*A->BlockSize*A->BlockSize + iBlock*A->BlockSize + jBlock];
        j++;
      }
      j+=A->BlockSize-BlockSize[A->N];
      j+=A->BlockSize+iBlock;
    }
    // Tridiagonalize the Block tridiagonal matrix
    dsbtrd(VECT, UPLO, &N, &KD, AB, &LDAB, &(M.a[1]), &(M.b[1]), Q.VRe, &N, WORK, &INFO);
    free(AB);
    free(WORK);
    if(INFO!=0)
    {
      printf("ERROR: dsbtrd failed in BlockTriDiagonalMatrixEigenSystem with INFO %i\n",INFO);
      return 1;
    }
  }
  else
  {
    Complex16 * AB;
    Complex16 *WORK;
    AB= (Complex16*)calloc(LDAB*N,sizeof(Complex16));
    if(AB==NULL)
    {
      printf("ERROR: malloc failed in BlockTriDiagonalMatrixEigenSystem\n");
      return 1;
    }
    WORK= (Complex16*)malloc(N*sizeof(Complex16));
    if(WORK==NULL)
    {
      free(AB);
      printf("ERROR: malloc failed in BlockTriDiagonalMatrixEigenSystem\n");
      return 1;
    }
    //fill matrix AB
    j=0;
    for(k=1;k<A->N;k++)
    {
      for(iBlock=0;iBlock<BlockSize[k];iBlock++)
      {
        for(jBlock=iBlock;jBlock<BlockSize[k];jBlock++)
        {
          AB[j]=A->a[k*A->BlockSize*A->BlockSize + iBlock*A->BlockSize + jBlock];
          j++;
        }
        j+=A->BlockSize-BlockSize[k];
        for(jBlock=0;jBlock<BlockSize[k+1];jBlock++)
        {
          AB[j]=A->b[k*A->BlockSize*A->BlockSize + iBlock*A->BlockSize + jBlock];
          j++;
        }
        j+=A->BlockSize-BlockSize[k+1];
        j+=iBlock;
      }
    }
    k=A->N;
    for(iBlock=0;iBlock<BlockSize[A->N];iBlock++)
    {
      for(jBlock=iBlock;jBlock<BlockSize[A->N];jBlock++)
      {
        AB[j]=A->a[k*A->BlockSize*A->BlockSize + iBlock*A->BlockSize + jBlock];
        j++;
      }
      j+=A->BlockSize-BlockSize[A->N];
      j+=A->BlockSize+iBlock;
    }
    // Tridiagonalize the Block tridiagonal matrix
    zhbtrd(VECT, UPLO, &N, &KD, AB, &LDAB, &(M.a[1]), &(M.b[1]), Q.V, &N, WORK, &INFO);
    free(AB);
    free(WORK);
    if(INFO!=0)
    {
      printf("ERROR: zhbtrd failed in BlockTriDiagonalMatrixEigenSystem with INFO %i\n",INFO);
      return 1;
    }
  }
  TriDiagonalMatrixToEigensystemCompact(&M, val, &funQ);
  CompactMatrixMultiply('N', 'N', &funQ, &Q, fun);
  CompactMatrixConjugate(fun);
  FreeCompactMatrix(&Q);
  FreeCompactMatrix(&funQ);
  FreeTriDiagonalMatrix(&M);
  return 0;
}

uint32_t BlockTriDiagonalMatrixToBlockListOfPoles(BlockTriDiagonalMatrixType const * A, BlockListOfPoles * const G)
{
  //TODO: all
  printf("not yet implemented BlockTriDiagonalMatrixToBlockListOfPoles\n");
  return 1;
}

uint32_t BlockTriDiagonalMatrixToBlockListOfPolesVariableGrid(BlockTriDiagonalMatrixType const * A, BlockListOfPoles * const G, const double * const Grid, const uint32_t NGridPoints)
{
  //TODO: all
  printf("not yet implemented BlockTriDiagonalMatrixToBlockListOfPolesVariableGrid\n");
  return 1;
}

/*  On entry A represents a Green's function in a block tridiagonal form
 *  On excit G represents the same Green's function by a list of poles (energy weigth)
 *  The energy mesh of G is given by Emin to Emax in steps of NE
 *
 *  This routine goes via a dense representation of the matrix and should be avoided, or replaced
 */
uint32_t BlockTriDiagonalMatrixToBlockListOfPolesFixedGrid(BlockTriDiagonalMatrixType const * A, BlockListOfPoles * const G, double const Emin, double const Emax, uint32_t const NE)
{
  int i,j,k,irho,kfun,iE;
  CompactMatrixType Afull, Afun, Rho, b0;
  double energy, dE=(Emax-Emin)/((double)NE), rest, totmin=0, totmax=0, totnow;
  double *Aval;
  /* first turn the matrix into compact form */
  if(BlockTriDiagonalMatrixToCompactMatrix(A,&Afull))
  {
    printf("BlockTriDiagonalMatrixToCompactMatrix failed in BlockTriDiagonalMatrixToBlockListOfPolesFixedGrid\n");
    return 1;
  }
  /* calcualte the eigen vectors and eigen values */
  Aval= (double*)malloc(A->BlockSize*A->N*sizeof(double));
  if(CompactMatrixEigensystem(&Afull, &Afun, Aval))
  {
    printf("CompactMatrixEigensystem failed in BlockTriDiagonalMatrixToBlockListOfPolesFixedGrid\n");
    return 1;
  }
  /* initiallize the memory for the block list of poles */
  G->BlockSize=A->BlockSize;
  G->N=NE+3;
  G->QComplex=A->QComplex;
  strcpy(G->Name, A->Name);
  if(InitBlockListOfPoles(G))
  {
    printf("InitBlockListOfPoles failed in BlockTriDiagonalMatrixToBlockListOfPolesFixedGrid\n");
    return 1;
  }
  energy=Emin-dE;
  /* initiallize the memory for rho of a specific eigen vector */
  Rho.QComplex=Afun.QComplex;
  Rho.Ni=A->BlockSize;
  Rho.Nj=A->BlockSize;
  if(InitCompactMatrix(&Rho))
  {
    printf("InitCompactMatrix failed in BlockTriDiagonalMatrixToBlockListOfPolesFixedGrid\n");
    return 1;
  }
  if(A->QComplex)
  {
    /* set the energy of the poles */
    for(i=0;i<G->N*(1+G->BlockSize*G->BlockSize);i+=(1+G->BlockSize*G->BlockSize))
    {
      G->V[i].Re=energy;
      G->V[i].Im=0;
      energy+=dE;
    }
    /* set the weight of the poles */
    kfun=0;
    for(k=0;k<Afun.Ni;k++)
    {
      /* create beta^2 */
      irho=0;
      for(i=0;i<A->BlockSize;i++)
      {
        for(j=0;j<A->BlockSize;j++)
        {
          Rho.V[irho].Re = Afun.V[kfun+i].Re*Afun.V[kfun+j].Re + Afun.V[kfun+i].Im*Afun.V[kfun+j].Im;
          Rho.V[irho].Im =-Afun.V[kfun+i].Re*Afun.V[kfun+j].Im + Afun.V[kfun+i].Im*Afun.V[kfun+j].Re;
          irho++;
        }
      }
      kfun+=Afun.Nj;
      /* add beta^2 to the correct position on the list of poles */
      rest = (Aval[k] - Emin) / dE;
      iE = (int)(floor(rest));
      if (iE >= (int)NE)
      {
        totnow=0;
        for(i=0;i<A->BlockSize;i++)
        {
          totnow+=Rho.V[i+i*A->BlockSize].Re;
        }
        if(totmax+totnow<10*DBL_MIN)
        {
          totmax+=10*DBL_MIN;
          totnow+=10*DBL_MIN;
        }
        G->V[(NE+2)*(1+G->BlockSize*G->BlockSize)].Re = (totmax * G->V[(NE+2)*(1+G->BlockSize*G->BlockSize)].Re + totnow * Aval[k]) / (totmax + totnow);
        irho=0;
        for(i=0;i<A->BlockSize;i++)
        {
          for(j=0;j<A->BlockSize;j++)
          {
            G->V[(NE+2)*(1+G->BlockSize*G->BlockSize)+1+j+i*A->BlockSize].Re+=Rho.V[irho].Re;
            G->V[(NE+2)*(1+G->BlockSize*G->BlockSize)+1+j+i*A->BlockSize].Im+=Rho.V[irho].Im;
            irho++;
          }
        }
        totmax+=totnow;
      }
      else if (iE < 0)
      {
        totnow=0;
        irho=0;
        for(i=0;i<A->BlockSize;i++)
        {
          totnow+=Rho.V[irho].Re;
          irho+=A->BlockSize+1;
        }
        if(totmin+totnow<10*DBL_MIN)
        {
          totmin+=10*DBL_MIN;
          totnow+=10*DBL_MIN;
        }
        G->V[0].Re = (totmin * G->V[0].Re + totnow * Aval[k]) / (totmin + totnow);
        irho=0;
        for(i=0;i<A->BlockSize;i++)
        {
          for(j=0;j<A->BlockSize;j++)
          {
            G->V[1+j+i*A->BlockSize].Re+=Rho.V[irho].Re;
            G->V[1+j+i*A->BlockSize].Im+=Rho.V[irho].Im;
            irho++;
          }
        }
        totmin+=totnow;
      }
      else
      {
        rest -= iE;  /* rest is a number between 0 and 1 pole at iE+1 gets 1-rest pole at iE+2 gest rest */
        irho=0;
        for(i=0;i<A->BlockSize;i++)
        {
          for(j=0;j<A->BlockSize;j++)
          {
            G->V[(iE+1)*(1+G->BlockSize*G->BlockSize)+1+j+i*A->BlockSize].Re+=(1-rest)*Rho.V[irho].Re;
            G->V[(iE+2)*(1+G->BlockSize*G->BlockSize)+1+j+i*A->BlockSize].Re+=(  rest)*Rho.V[irho].Re;
            G->V[(iE+1)*(1+G->BlockSize*G->BlockSize)+1+j+i*A->BlockSize].Im+=(1-rest)*Rho.V[irho].Im;
            G->V[(iE+2)*(1+G->BlockSize*G->BlockSize)+1+j+i*A->BlockSize].Im+=(  rest)*Rho.V[irho].Im;
            irho++;
          }
        }
      }
    }
  }
  else /*real */
  {
    /* set the energy of the poles */
    for(i=0;i<G->N*(1+G->BlockSize*G->BlockSize);i+=(1+G->BlockSize*G->BlockSize))
    {
      G->VRe[i]=energy;
      energy+=dE;
    }
    /* set the weight of the poles */
    kfun=0;
    for(k=0;k<Afun.Ni;k++)
    {
      /* create beta^2 */
      irho=0;
      for(i=0;i<A->BlockSize;i++)
      {
        for(j=0;j<A->BlockSize;j++)
        {
          Rho.VRe[irho] = Afun.VRe[kfun+i]*Afun.VRe[kfun+j];
          irho++;
        }
      }
      kfun+=Afun.Nj;
      /* add beta^2 to the correct position on the list of poles */
      rest = (Aval[k] - Emin) / dE;
      iE = (int)(floor(rest));
      if (iE >= (int)NE)
      {
        totnow=0;
        for(i=0;i<A->BlockSize;i++)
        {
          totnow+=Rho.VRe[i+i*A->BlockSize];
        }
        if(totmax+totnow<10*DBL_MIN)
        {
          totmax+=10*DBL_MIN;
          totnow+=10*DBL_MIN;
        }
        G->VRe[(NE+2)*(1+G->BlockSize*G->BlockSize)] = (totmax * G->VRe[(NE+2)*(1+G->BlockSize*G->BlockSize)] + totnow * Aval[k]) / (totmax + totnow);
        irho=0;
        for(i=0;i<A->BlockSize;i++)
        {
          for(j=0;j<A->BlockSize;j++)
          {
            G->VRe[(NE+2)*(1+G->BlockSize*G->BlockSize)+1+j+i*A->BlockSize]+=Rho.VRe[irho];
            irho++;
          }
        }
        totmax+=totnow;
      }
      else if (iE < 0)
      {
        totnow=0;
        irho=0;
        for(i=0;i<A->BlockSize;i++)
        {
          totnow+=Rho.VRe[irho];
          irho+=A->BlockSize+1;
        }
        if(totmin+totnow<10*DBL_MIN)
        {
          totmin+=10*DBL_MIN;
          totnow+=10*DBL_MIN;
        }
        G->VRe[0] = (totmin * G->VRe[0] + totnow * Aval[k]) / (totmin + totnow);
        irho=0;
        for(i=0;i<A->BlockSize;i++)
        {
          for(j=0;j<A->BlockSize;j++)
          {
            G->VRe[1+j+i*A->BlockSize]+=Rho.VRe[irho];
            irho++;
          }
        }
        totmin+=totnow;
      }
      else
      {
        rest -= iE;  /* rest is a number between 0 and 1 pole at iE+1 gets 1-rest pole at iE+2 gest rest */
        irho=0;
        for(i=0;i<A->BlockSize;i++)
        {
          for(j=0;j<A->BlockSize;j++)
          {
            G->VRe[(iE+1)*(1+G->BlockSize*G->BlockSize)+1+j+i*A->BlockSize]+=(1-rest)*Rho.VRe[irho];
            G->VRe[(iE+2)*(1+G->BlockSize*G->BlockSize)+1+j+i*A->BlockSize]+=(  rest)*Rho.VRe[irho];
            irho++;
          }
        }
      }
    }
  }
  b0.QComplex=A->QComplex;
  b0.Ni=A->BlockSize;
  b0.Nj=A->BlockSize;
  if(b0.QComplex)
  {
    b0.V=A->b;
  }
  else
  {
    b0.VRe=A->bRe;
  }
  BlockListOfPolesTransformInlineNC(G, &b0);//
  ChopBlockListOfPoles(G, DBL_EPSILON);
  G->mu = A->mu;
  FreeCompactMatrix(&Rho);
  FreeCompactMatrix(&Afull);
  FreeCompactMatrix(&Afun);
  return 0;
}


/*
 *
 * Takes a blocktridiagonal matrix as input
 * Diagonalizes the bath part of this matrix, and uses the eigenvectors to transform the
 * block coupling the impurity to the bath into a row of blocks coupling the impurity
 * to the new Anderson sites with energies obtained from eigenvalues.
 *
 * !!Does not work as intended yet!!
 * Two issues:
 * 1. I get mixed hopping in the anderson matrix: xz hops to both xz/yz, and same for yz.
 * 2. The bath orbitals become reordered because the eigenvalues from the diagonalization
 *    are sorted. Not sure if this is an issue.
 *
 * End result is a spectrum calculated after using this is slightly different from the original
 * spectrum using the BlockTridiagonalMatrix, and I think it is due to issue #1.
 *
 */
uint32_t BlockTriDiagonalMatrixToBlockAndersonMatrix(BlockTriDiagonalMatrixType const ATri, BlockAndersonMatrixType * const B)
{
  int i,j,bpos,cpos,k;
  CompactMatrixType Afull, Afun, C,D,E;
  double *Aval;
  
  /* Convert the portion after the first block to compact matrix */
  if(BlockTriDiagonalMatrixPartToCompactMatrix(&ATri,&Afull,1,ATri.N-1))
  {
    printf("BlockTriDiagonalMatrixPartToCompactMatrix failed in BlockTriDiagonalMatrixToBlockAndersonMatrix\n");
    return 1;
  }
  /* Diagonalize compact matrix */
  Aval= (double*)malloc(ATri.BlockSize*(ATri.N-1)*sizeof(double));
  if(CompactMatrixEigensystem(&Afull, &Afun, Aval))
  {
    printf("CompactMatrixEigensystem failed in BlockTriDiagonalMatrixToBlockAndersonMatrix\n");
    return 1;
  }
  B->QComplex = ATri.QComplex;
  B->N = ATri.N;
  B->BlockSize = ATri.BlockSize;
  InitBlockAndersonMatrix(B);
  if(ATri.QComplex)
  {
    /* copy A0 and B0 */
    bpos = 0;
    for(i = 0; i < B->BlockSize; i++)
    {
      for(j = 0; j < B->BlockSize; j++)
      {
        B->a[bpos].Re = ATri.a[bpos].Re;
        B->a[bpos].Im = ATri.a[bpos].Im;
        B->b[bpos].Re = ATri.b[bpos].Re;
        B->b[bpos].Im = ATri.b[bpos].Im;
        bpos++;
      }
    }
    /* First diagonal block is same as ATri first diagonal block */
    bpos = B->BlockSize*B->BlockSize;
    for(i = 0; i < B->BlockSize; i++)
    {
      for(j = 0; j < B->BlockSize; j++)
      {
        B->a[bpos].Re = ATri.a[bpos].Re;
        B->a[bpos].Im = ATri.a[bpos].Im;
        bpos++;
      }
    }
    /* Set diagonals to eigenvalues */
    cpos = 0;
    for(k = 1; k < ATri.N; k++)
    {
      for(i = 0; i < B->BlockSize; i++)
      {
        B->a[(k+1)*B->BlockSize*B->BlockSize + i*B->BlockSize + i].Re = Aval[cpos++];
      }
    }
    /* Transform first block of ATri using eigenfunctions from diagonalization */
    C.QComplex = B->QComplex;
    C.Ni = B->BlockSize;
    C.Nj = B->BlockSize;
    InitCompactMatrix(&C);
    bpos = B->BlockSize*B->BlockSize;
    for(i = 0; i < B->BlockSize*B->BlockSize; i++)
    {
      C.V[i].Re = ATri.b[bpos++].Re;
      C.V[i].Im = ATri.b[bpos++].Im;
    }
    D.QComplex = B->QComplex;
    D.Ni = B->BlockSize;
    D.Nj = B->BlockSize*(B->N-1);
    InitCompactMatrix(&D);
    for(i = 0; i < B->BlockSize; i++)
    {
      for(j = 0; j < B->BlockSize*(B->N-1); j++)
      {
        D.V[i*B->BlockSize*(B->N-1)+ j].Re = Afun.V[j*B->BlockSize*(B->N-1)+i].Re;
        D.V[i*B->BlockSize*(B->N-1) + j].Im = Afun.V[j*B->BlockSize*(B->N-1)+i].Im;
      }
    }
    /* Transform the block that coupled impurity to chain into a row of blocks
     coupling impurity to anderson bath */
    CompactMatrixMultiply('N', 'N', &C, &D, &E);
    /* Now copy E into the beta values of block matrix*/
    bpos = B->BlockSize*B->BlockSize;
    for(k = 0; k < B->N-1; k++)
    {
      for(i = 0; i < B->BlockSize; i++)
      {
        for(j = 0; j < B->BlockSize; j++)
        {
          B->b[bpos].Re   = E.V[i*B->BlockSize*(B->N-1)+k*B->BlockSize+j].Re;
          B->b[bpos++].Im = E.V[i*B->BlockSize*(B->N-1)+k*B->BlockSize+j].Im;
        }
      }
    }
  }
  else
  {
    /* Real */
    /* copy A0 and B0 */
    bpos = 0;
    for(i = 0; i < B->BlockSize; i++)
    {
      for(j = 0; j < B->BlockSize; j++)
      {
        B->aRe[bpos] = ATri.aRe[bpos];
        B->bRe[bpos] = ATri.bRe[bpos];
        bpos++;
      }
    }
    /* First diagonal block is same as ATri first diagonal block */
    bpos = B->BlockSize*B->BlockSize;
    for(i = 0; i < B->BlockSize; i++)
    {
      for(j = 0; j < B->BlockSize; j++)
      {
        B->aRe[bpos] = ATri.aRe[bpos];
        bpos++;
      }
    }
    /* Set diagonals to eigenvalues */
    cpos = 0;
    for(k = 1; k < ATri.N; k++)
    {
      for(i = 0; i < B->BlockSize; i++)
      {
        B->aRe[(k+1)*B->BlockSize*B->BlockSize + i*B->BlockSize + i] = Aval[cpos++];
      }
    }
    /* Transform first block of ATri using eigenfunctions from diagonalization */
    C.QComplex = B->QComplex;
    C.Ni = B->BlockSize;
    C.Nj = B->BlockSize;
    InitCompactMatrix(&C);
    bpos = B->BlockSize*B->BlockSize;
    for(i = 0; i < B->BlockSize*B->BlockSize; i++)
    {
      C.VRe[i] = ATri.bRe[bpos++];
    }
    D.QComplex = B->QComplex;
    D.Ni = B->BlockSize;
    D.Nj = B->BlockSize*(B->N-1);
    InitCompactMatrix(&D);
    for(i = 0; i < B->BlockSize; i++)
    {
      for(j = 0; j < B->BlockSize*(B->N-1); j++)
      {
        D.VRe[i*B->BlockSize*(B->N-1) + j] = Afun.VRe[j*B->BlockSize*(B->N-1)+i];
      }
    }
    /* Transform the block that coupled impurity to chain into a row of blocks
     coupling impurity to anderson bath */
    CompactMatrixMultiply('N', 'N', &C, &D, &E);
    /* Now copy E into the beta values of anderson block matrix*/
    bpos = B->BlockSize*B->BlockSize;
    for(k = 0; k < B->N-1; k++)
    {
      for(i = 0; i < B->BlockSize; i++)
      {
        for(j = 0; j < B->BlockSize; j++)
        {
          B->bRe[bpos++] = E.VRe[i*(B->N-1)*B->BlockSize+k*B->BlockSize+j];
        }
      }
    }
  }
  FreeCompactMatrix(&C);
  FreeCompactMatrix(&D);
  FreeCompactMatrix(&E);
  FreeCompactMatrix(&Afun);
  FreeCompactMatrix(&Afull);
  free(Aval);
  return 0;
}

uint32_t BlockTriDiagonalMatrixToBlockNaturalImpurityMatrix(BlockTriDiagonalMatrixType const ATri, BlockNaturalImpurityOrbitalMatrixType * const B)
{
  //TODO: all
  printf("Not implemented BlockTriDiagonalMatrixToBlockNaturalImpurityMatrix\n");
  return 1;
}


/*  =========================================================================
 *  BlockNaturalImpurityOrbitalMatrixType to __
 *  ========================================================================= */


/*  G = (T Gp T^+)[1,1]
 *
 *       (Gv 0  )
 *  Gp = (0  Gc )
 *
 *       ( bv0    bc0 )
 *  T  = (~bc0  ~-bv0 )
 *
 *  G = bc0 gc bc0^+ + bv0 gv bv0^+
 */
uint32_t BlockNaturalImpurityOrbitalToSpectra(BlockNaturalImpurityOrbitalMatrixType const * A, SpectraType *const G, double const Emin, double const Emax, uint32_t const N, double const Gamma)
{
  SpectraType Gval[A->BlockSize*A->BlockSize], Gcon[A->BlockSize*A->BlockSize];
  uint32_t i,j;

  BlockTriDiagonalMatrixToSpectra(0, &(A->val), Gval, Emin, Emax, N, Gamma);
  BlockTriDiagonalMatrixToSpectra(0, &(A->con), Gcon, Emin, Emax, N, Gamma);

  G[0].N=N;
  G[0].Emin=Emin;
  G[0].Emax=Emax;
  G[0].EminPole=Emin;
  G[0].EmaxPole=Emax;
  G[0].dE=(Emax-Emin)/((double)N);
  G[0].Gamma=Gamma;
  for(i=0;i<A->BlockSize*A->BlockSize;i++)
  {
    SpectrumAllocate(G[0], &(G[i]));
  }

  if(A->QComplex)
  {
    for(i=0;i<A->BlockSize*A->BlockSize;i++)
    {
      for(j=0;j<=N;j++)
      {
        G[i].Re[j] = Gval[i].Re[j] + Gcon[i].Re[j];
        G[i].Im[j] = Gval[i].Im[j] + Gcon[i].Im[j];
      }
    }
  }
  else
  {
    for(i=0;i<A->BlockSize*A->BlockSize;i++)
    {
      for(j=0;j<=N;j++)
      {
        G[i].Re[j] = Gval[i].Re[j] + Gcon[i].Re[j];
      }
    }
  }

  for(i=0;i<A->BlockSize*A->BlockSize;i++)
  {
    SpectrumFree(&Gval[i]);
    SpectrumFree(&Gcon[i]);
  }
  return 0;
}

uint32_t BlockNaturalImpurityOrbitalMatrixEigenValues(const BlockNaturalImpurityOrbitalMatrixType * const A, double * val)
{
  //TODO all
  return 0;
}
uint32_t BlockNaturalImpurityOrbitalMatrixEigenSystem(const BlockNaturalImpurityOrbitalMatrixType * const A, double * val, CompactMatrixType * fun)
{
  //TODO all
  return 0;
}

uint32_t BlockNaturalImpurityOrbitalMatrixToBlockAndersonMatrix(BlockNaturalImpurityOrbitalMatrixType const AAnd, BlockAndersonMatrixType * const B)
{
  //TODO: all
  printf("Not implemented BlockNaturalImpurityOrbitalMatrixToBlockAndersonMatrix\n");
  return 1;
}
uint32_t BlockNaturalImpurityOrbitalMatrixToBlockListOfPoles(BlockNaturalImpurityOrbitalMatrixType const AAnd, BlockListOfPoles * const B)
{
  //TODO: all
  printf("Not implemented BlockNaturalImpurityOrbitalMatrixToBlockListOfPoles\n");
  return 1;
}
uint32_t BlockNaturalImpurityOrbitalMatrixToBlockListOfPolesVariableGrid(BlockNaturalImpurityOrbitalMatrixType const AAnd, BlockListOfPoles * const B, const double * const Grid, const uint32_t NGridPoints)
{
  //TODO: all
  printf("Not implemented BlockNaturalImpurityOrbitalMatrixToBlockListOfPolesVariableGrid\n");
  return 1;
}

uint32_t BlockNaturalImpurityOrbitalMatrixToBlockTridiagonalMatrix(BlockNaturalImpurityOrbitalMatrixType const AAnd, BlockTriDiagonalMatrixType * const B)
{
  //TODO: all
  printf("Not implemented BlockNaturalImpurityOrbitalMatrixToBlockTridiagonalMatrix\n");
  return 1;
}



/*  =========================================================================
 *  Inversed BlockListOfPoles to __
 *  ========================================================================= */


/*  Turn a block list of poles into a matrix
 *       ( a0   b1  b2  b3  ..  bn  )
 *       ( b1^+ a1  0   0   ..  0   )
 *   B = ( b2^+ 0   a2  0   ..  0   )
 *       ( ..   ..  ..  ..  ..  ..  )
 *       ( bn^+ 0   0   0   ..  an  )
 *
 *   The greensfunction (w-B)^-1_[blocksize,blocksize] is equal to the inverse of the greensfunction defined by BlockListOfPoles
 *
 */
uint32_t BlockListOfPolesInversedToCompactMatrix(BlockListOfPoles const * const A, CompactMatrixType * const B)
{
  uint32_t i,j,k,kA,kBa;
  B->QComplex = A->QComplex;
  B->Ni = (A->N+1) * A->BlockSize;
  B->Nj = (A->N+1) * A->BlockSize;
  B->mu = A->mu;
  InitCompactMatrix(B);
  CompactMatrixType Anow,Asqrt;
  Anow.Ni=A->BlockSize;
  Anow.Nj=A->BlockSize;
  Anow.Name=(char*)malloc(MAXStrLength*sizeof(char));
  strcpy(Anow.Name, "CompactMatrix");
  if(A->QComplex)
  {
    Anow.QComplex=1;
    for(i=0;i<A->BlockSize;i++)
    {
      for(j=0;j<A->BlockSize;j++)
      {
        B->V[i*B->Nj+j].Re=A->a0[i*A->BlockSize+j].Re;
        B->V[i*B->Nj+j].Im=A->a0[i*A->BlockSize+j].Im;
      }
    }
    kA=0;
    kBa=A->BlockSize*(B->Nj+1);
    for(k=0;k<A->N;k++)
    {
      Anow.V=&(A->V[kA+1]);
      CompactMatrixSqrt(&Anow, &Asqrt);
      for(i=0;i<A->BlockSize;i++)
      {
        B->V[kBa].Re=A->V[kA].Re;
        kBa+=B->Nj+1;
        for(j=0;j<A->BlockSize;j++)
        {
          B->V[(i+(k+1)*A->BlockSize)*B->Nj+j].Re = Asqrt.V[i*A->BlockSize+j].Re;
          B->V[(i+(k+1)*A->BlockSize)+j*B->Nj].Re = Asqrt.V[i*A->BlockSize+j].Re;
          B->V[(i+(k+1)*A->BlockSize)*B->Nj+j].Im = Asqrt.V[i*A->BlockSize+j].Im;
          B->V[(i+(k+1)*A->BlockSize)+j*B->Nj].Im =-Asqrt.V[i*A->BlockSize+j].Im;
        }
      }
      FreeCompactMatrix(&Asqrt);
      kA+=A->BlockSize*A->BlockSize+1;
    }
  }
  else
  {
    Anow.QComplex=0;
    for(i=0;i<A->BlockSize;i++)
    {
      for(j=0;j<A->BlockSize;j++)
      {
        B->VRe[i*B->Nj+j]=A->a0Re[i*A->BlockSize+j];
      }
    }
    kA=0;
    kBa=A->BlockSize*(B->Nj+1);
    for(k=0;k<A->N;k++)
    {
      Anow.VRe=&(A->VRe[kA+1]);
      CompactMatrixSqrt(&Anow, &Asqrt);
      for(i=0;i<A->BlockSize;i++)
      {
        B->VRe[kBa]=A->VRe[kA];
        kBa+=B->Nj+1;
        for(j=0;j<A->BlockSize;j++)
        {
          B->VRe[(i+(k+1)*A->BlockSize)*B->Nj+j]=Asqrt.VRe[i*A->BlockSize+j];
          B->VRe[(i+(k+1)*A->BlockSize)+j*B->Nj]=Asqrt.VRe[i*A->BlockSize+j];
        }
      }
      FreeCompactMatrix(&Asqrt);
      kA+=A->BlockSize*A->BlockSize+1;
    }
  }
  return 0;
}

/*  =========================================================================
 *  TightBindingDef to __
 *  ========================================================================= */


/*  Given an tight binding Hamiltonian and a k-mesh defined by the number of divisions of a, b and c (Nk)
 *  returns the green function matrix of size A.NFermion by A.NFermion
 *          on an energy mesh with N+3 points. N+1 points from Emin to Emax and one point below Emin and one above Emax such that all weight and first moment is ok
 *  Calculates the eigen values and functions for each k
 *  for each k and eigen value i
 *  weight is added to the density matrix according to Fermi-dirac statistics at temperature kBT and chemical potential mu
 */
uint32_t TightBindingHamiltonianToDensityMatrix(TightBindingDef const * const A, int Nk[3], CompactMatrixType *rho, double mu, double kbT)
{
  double q[3],diq[3],val[A->NFermion],weigth;
  int iq[3],ival,k,Nktot=Nk[0]*Nk[1]*Nk[2];
  MatrixType H;
  CompactMatrixType T;
  
  rho->Ni=A->NFermion;
  rho->Nj=A->NFermion;
  rho->QComplex=1;
  InitCompactMatrix(rho);
  rho->mu = mu;
  rho->a0=0;
  rho->b0=1;
  strcpy(rho->Name, A->Name);
  for (iq[0] = 0; iq[0] < Nk[0]; iq[0]++)
  {
    for (iq[1] = 0; iq[1] < Nk[1]; iq[1]++)
    {
      for (iq[2] = 0; iq[2] < Nk[2]; iq[2]++)
      {
        diq[0] = iq[0] / ((double)(Nk[0]));
        diq[1] = iq[1] / ((double)(Nk[1]));
        diq[2] = iq[2] / ((double)(Nk[2]));
        q[0] = diq[0] * A->Cell.gUnitCell[0][0] + diq[1] * A->Cell.gUnitCell[1][0] + diq[2] * A->Cell.gUnitCell[2][0];
        q[1] = diq[0] * A->Cell.gUnitCell[0][1] + diq[1] * A->Cell.gUnitCell[1][1] + diq[2] * A->Cell.gUnitCell[2][1];
        q[2] = diq[0] * A->Cell.gUnitCell[0][2] + diq[1] * A->Cell.gUnitCell[1][2] + diq[2] * A->Cell.gUnitCell[2][2];
        CreateHamiltonianFromTightBindingDefConstantPhaseInUnitCell(*A, q, &H);
        EigensystemToCompact(H, &T, val);
        for(ival=0;ival<A->NFermion;ival++)
        {
          /* calculate contribution for eigen-value ival */
          weigth = FermiFunction(val[ival], mu, kbT);
          if(H.QComplex)
          {
            ComplexEigenFunctionToCharacter(&(T.V[ival*A->NFermion]), &H);
            for(k=0;k<A->NFermion*A->NFermion;k++)
            {
              rho->V[k].Re += weigth * H.VRe[0][k];
              rho->V[k].Im += weigth * H.VIm[0][k];
            }
          }
          else
          {
            RealEigenFunctionToCharacter(&(T.VRe[ival*A->NFermion]), &H);
            for(k=0;k<A->NFermion*A->NFermion;k++)
            {
              rho->V[k].Re += weigth * H.VRe[0][k];
            }
          }
        }
        FreeMatrix(&H);
      }
    }
  }
  
  for(k=0;k<rho->Ni*rho->Nj;k++)
  {
    rho->V[k].Re/=Nktot;
    rho->V[k].Im/=Nktot;
  }
  
  return 0;
}



/*  Given an tight binding Hamiltonian and a k-mesh defined by the number of divisions of a, b and c (Nk)
 *  returns the green function matrix of size A.NFermion by A.NFermion
 *          on an energy mesh with N+3 points. N+1 points from Emin to Emax and one point below Emin and one above Emax such that all weight and first moment is ok
 *  Calculates the eigen values and functions for each k
 *  for each k and eigen value i
 *  weight is added to the pole at iE and iE+1 such that Energy[iE] < val < Energy[iE+1]
 *  The matrix h_{jk} is fun[i,j] * fun[i,k]^* is added
 */
uint32_t TightBindingHamiltonianToBlockListOfPoles(TightBindingDef const A, int Nk[3], BlockListOfPoles * const G, double const Emin, double const Emax, uint32_t const NE)
{
  MatrixType H;
  CompactMatrixType T;
  double q[3], val[A.NFermion], rest, totmin=0, totmax=0, totnow;
  int    iq[3],ival,iE,i,j,ii;
  double diq[3],dE=(Emax-Emin)/((double)NE),Nktot;
  G->BlockSize=A.NFermion;
  G->N=NE+3;
  G->QComplex=1;
  strcpy(G->Name, A.Name);
  if(InitBlockListOfPoles(G))
  {
    printf("InitBlockListOfPoles failed in TightBindingHamiltonianToBlockListOfPoles\n");
    return 1;
  }
  ii=-1;
  for(i=0;i<G->N*(1+G->BlockSize*G->BlockSize);i+=(1+G->BlockSize*G->BlockSize))
  {
    G->V[i].Re=Emin + ii * dE;
    ii++;
  }
  for (iq[0] = 0; iq[0] < Nk[0]; iq[0]++)
  {
    for (iq[1] = 0; iq[1] < Nk[1]; iq[1]++)
    {
      for (iq[2] = 0; iq[2] < Nk[2]; iq[2]++)
      {
        diq[0] = iq[0] / ((double)(Nk[0]));
        diq[1] = iq[1] / ((double)(Nk[1]));
        diq[2] = iq[2] / ((double)(Nk[2]));
        q[0] = diq[0] * A.Cell.gUnitCell[0][0] + diq[1] * A.Cell.gUnitCell[1][0] + diq[2] * A.Cell.gUnitCell[2][0];
        q[1] = diq[0] * A.Cell.gUnitCell[0][1] + diq[1] * A.Cell.gUnitCell[1][1] + diq[2] * A.Cell.gUnitCell[2][1];
        q[2] = diq[0] * A.Cell.gUnitCell[0][2] + diq[1] * A.Cell.gUnitCell[1][2] + diq[2] * A.Cell.gUnitCell[2][2];
        CreateHamiltonianFromTightBindingDefConstantPhaseInUnitCell(A, q, &H);
        EigensystemToCompact(H, &T, val);
        //printf("\n");
        for(ival=0;ival<A.NFermion;ival++)
        {
          //printf("%f ",val[ival]);
          /* calculate matrix B for eigen-value ival */
          if(H.QComplex)
          {
            ComplexEigenFunctionToCharacter(&(T.V[ival*A.NFermion]), &H);
          }
          else
          {
            RealEigenFunctionToCharacter(&(T.VRe[ival*A.NFermion]), &H);
          }
          /* determine the position of the pole in the list of poles */
          rest = (val[ival] - Emin) / dE;
          iE = (int)(floor(rest));
          if (iE >= (int)NE)
          {
            totnow=0;
            for(i=0;i<A.NFermion;i++)
            {
              totnow+=H.VRe[i][i];
            }
            if(totmax+totnow<10*DBL_MIN)
            {
              totmax+=10*DBL_MIN;
              totnow+=10*DBL_MIN;
            }
            G->V[(NE+2)*(1+G->BlockSize*G->BlockSize)].Re = (totmax * G->V[(NE+2)*(1+G->BlockSize*G->BlockSize)].Re + totnow * val[ival]) / (totmax + totnow);
            for(i=0;i<A.NFermion;i++)
            {
              for(j=0;j<A.NFermion;j++)
              {
                G->V[(NE+2)*(1+G->BlockSize*G->BlockSize)+1+j+i*A.NFermion].Re+=H.VRe[i][j];
                if(H.QComplex)
                {
                  G->V[(NE+2)*(1+G->BlockSize*G->BlockSize)+1+j+i*A.NFermion].Im+=H.VIm[i][j];
                }
              }
            }
            totmax+=totnow;
          }
          else if (iE < 0)
          {
            totnow=0;
            for(i=0;i<A.NFermion;i++)
            {
              totnow+=H.VRe[i][i];
            }
            if(totmin+totnow<10*DBL_MIN)
            {
              totmin+=10*DBL_MIN;
              totnow+=10*DBL_MIN;
            }
            G->V[0].Re = (totmin * G->V[0].Re + totnow * val[ival]) / (totmin + totnow);
            for(i=0;i<A.NFermion;i++)
            {
              for(j=0;j<A.NFermion;j++)
              {
                G->V[1+j+i*A.NFermion].Re+=H.VRe[i][j];
                if(H.QComplex)
                {
                  G->V[1+j+i*A.NFermion].Im+=H.VIm[i][j];
                }
              }
            }
            totmin+=totnow;
          }
          else
          {
            rest -= iE;  /* rest is a number between 0 and 1 pole at iE+1 gets 1-rest pole at iE+2 gest rest */
            for(i=0;i<A.NFermion;i++)
            {
              for(j=0;j<A.NFermion;j++)
              {
                G->V[(iE+1)*(1+G->BlockSize*G->BlockSize)+1+j+i*A.NFermion].Re+=(1-rest)*H.VRe[i][j];
                G->V[(iE+2)*(1+G->BlockSize*G->BlockSize)+1+j+i*A.NFermion].Re+=(  rest)*H.VRe[i][j];
                if(H.QComplex)
                {
                  G->V[(iE+1)*(1+G->BlockSize*G->BlockSize)+1+j+i*A.NFermion].Im+=(1-rest)*H.VIm[i][j];
                  G->V[(iE+2)*(1+G->BlockSize*G->BlockSize)+1+j+i*A.NFermion].Im+=(  rest)*H.VIm[i][j];
                }
              }
            }
          }
        }
      }
      FreeMatrix(&H);
      FreeCompactMatrix(&T);
    }
  }
  Nktot=Nk[0]*Nk[1]*Nk[2];
  for(i=0;i<G->N*(1+G->BlockSize*G->BlockSize);i+=(1+G->BlockSize*G->BlockSize))
  {
    for(j=1;j<1+G->BlockSize*G->BlockSize;j++)
    {
      G->V[i+j].Re/=Nktot;
      G->V[i+j].Im/=Nktot;
    }
  }
  ChopBlockListOfPoles(G, DBL_EPSILON);
  return 0;
}




/*  Given an tight binding Hamiltonian and a k-mesh defined by the number of divisions of a, b and c (Nk)
 *  furhtermore given a list of unit cell translations given as vectors in units of a, b and c of length Nvec.
 *  returns the green function matrix of size Nvec * A.NFermion by Nvec * A.NFermion
 *          on an energy mesh with N+3 points. N+1 points from Emin to Emax and one point below Emin and one above Emax such that all weight and first moment is ok
 *  Calculates the eigen values and functions for each k
 *  for each k and eigen value i
 *  weight is added to the pole at iE and iE+1 such that Energy[iE] < val < Energy[iE+1]
 *  The matrix h_{jk} is fun[i,j] * fun[i,k]^* is added
 *  a Function at a shifted unit cell is multiplied by e^{i kvec . x }
 */
uint32_t TightBindingHamiltonianMultiSiteToBlockListOfPoles(TightBindingDef const A, int Nk[3], uint32_t const Nvec, int const *vec, BlockListOfPoles * const G, double const Emin, double const Emax, uint32_t const NE)
{
  MatrixType H;
  CompactMatrixType T;
  double q[3], val[A.NFermion], rest, totmin=0, totmax=0, totnow, dvec[3];
  double const piornot = ((A.Cell.Units&1) ? (1) : (Pi));  // if A.Cell.Units last bit is 1 then r.g = pi, otherwise it evaluates to 1, piornot * r.g now eveluates always to pi.
  int    iq[3],ival,iE,i,j,ivec,jvec;
  double diq[3],dE=(Emax-Emin)/((double)NE),energy,Nktot;
  Complex16 phase;
  G->BlockSize=A.NFermion * Nvec;
  G->N=NE+3;
  G->QComplex=1;
  strcpy(G->Name, A.Name);
  if(InitBlockListOfPoles(G))
  {
    printf("InitBlockListOfPoles failed in TightBindingHamiltonianToBlockListOfPoles\n");
    return 1;
  }
  energy=Emin-dE;
  for(i=0;i<G->N*(1+G->BlockSize*G->BlockSize);i+=(1+G->BlockSize*G->BlockSize))
  {
    G->V[i].Re=energy;
    energy+=dE;
  }
  for (iq[0] = 0; iq[0] < Nk[0]; iq[0]++)
  {
    for (iq[1] = 0; iq[1] < Nk[1]; iq[1]++)
    {
      for (iq[2] = 0; iq[2] < Nk[2]; iq[2]++)
      {
        diq[0] = iq[0] / ((double)(Nk[0]));
        diq[1] = iq[1] / ((double)(Nk[1]));
        diq[2] = iq[2] / ((double)(Nk[2]));
        q[0] = diq[0] * A.Cell.gUnitCell[0][0] + diq[1] * A.Cell.gUnitCell[1][0] + diq[2] * A.Cell.gUnitCell[2][0];
        q[1] = diq[0] * A.Cell.gUnitCell[0][1] + diq[1] * A.Cell.gUnitCell[1][1] + diq[2] * A.Cell.gUnitCell[2][1];
        q[2] = diq[0] * A.Cell.gUnitCell[0][2] + diq[1] * A.Cell.gUnitCell[1][2] + diq[2] * A.Cell.gUnitCell[2][2];
        CreateHamiltonianFromTightBindingDefConstantPhaseInUnitCell(A, q, &H);
        EigensystemToCompact(H, &T, val);
        for(ival=0;ival<A.NFermion;ival++)
        {
          /* calculate character matrix for eigen-value ival - note that H changes meaning from Hamiltonian to Character matrix */
          if(H.QComplex)
          {
            ComplexEigenFunctionToCharacterMultiSite(&(T.V[ival*A.NFermion]), &H, Nvec, vec);
          }
          else
          {
            RealEigenFunctionToCharacterMultiSite(&(T.VRe[ival*A.NFermion]), &H, Nvec, vec);
          }
          /* determine the position of the pole in the list of poles */
          rest = (val[ival] - Emin) / dE;
          iE = (int)(floor(rest));
          if (iE >= NE)
          {
            totnow=0;
            for(i=0;i<A.NFermion;i++)
            {
              totnow+=H.VRe[i][i];
            }
            if(totmax+totnow<10*DBL_MIN)
            {
              totmax+=10*DBL_MIN;
              totnow+=10*DBL_MIN;
            }
            G->V[(NE+2)*(1+G->BlockSize*G->BlockSize)].Re = (totmax * G->V[(NE+2)*(1+G->BlockSize*G->BlockSize)].Re + totnow * val[ival]) / (totmax + totnow);
            for(i=0;i<A.NFermion;i++)
            {
              for(ivec=0;ivec<Nvec;ivec++)
              {
                for(j=0;j<A.NFermion;j++)
                {
                  for(jvec=0;jvec<Nvec;jvec++)
                  {
                    dvec[0] = (vec[3*ivec+0]-vec[3*jvec+0]) * A.Cell.rUnitCell[0][0] + (vec[3*ivec+1]-vec[3*jvec+1]) * A.Cell.rUnitCell[1][0] + (vec[3*ivec+2]-vec[3*jvec+2]) * A.Cell.rUnitCell[2][0];
                    dvec[1] = (vec[3*ivec+0]-vec[3*jvec+0]) * A.Cell.rUnitCell[0][1] + (vec[3*ivec+1]-vec[3*jvec+1]) * A.Cell.rUnitCell[1][1] + (vec[3*ivec+2]-vec[3*jvec+2]) * A.Cell.rUnitCell[2][1];
                    dvec[2] = (vec[3*ivec+0]-vec[3*jvec+0]) * A.Cell.rUnitCell[0][2] + (vec[3*ivec+1]-vec[3*jvec+1]) * A.Cell.rUnitCell[1][2] + (vec[3*ivec+2]-vec[3*jvec+2]) * A.Cell.rUnitCell[2][2];
                    phase.Re=cos(piornot * Dot(dvec,q));
                    phase.Im=sin(piornot * Dot(dvec,q));
                    G->V[(NE+2)*(1+G->BlockSize*G->BlockSize)+1 + j+jvec*A.NFermion + (i+ivec*A.NFermion)*G->BlockSize].Re += phase.Re*H.VRe[i][j];
                    G->V[(NE+2)*(1+G->BlockSize*G->BlockSize)+1 + j+jvec*A.NFermion + (i+ivec*A.NFermion)*G->BlockSize].Im += phase.Im*H.VRe[i][j];
                    if(H.QComplex)
                    {
                      G->V[(NE+2)*(1+G->BlockSize*G->BlockSize)+1 + j+jvec*A.NFermion + (i+ivec*A.NFermion)*G->BlockSize].Re -= phase.Im*H.VIm[i][j];
                      G->V[(NE+2)*(1+G->BlockSize*G->BlockSize)+1 + j+jvec*A.NFermion + (i+ivec*A.NFermion)*G->BlockSize].Im += phase.Re*H.VIm[i][j];
                    }
                  }
                }
              }
            }
            totmax+=totnow;
          }
          else if (iE < 0)
          {
            totnow=0;
            for(i=0;i<A.NFermion;i++)
            {
              totnow+=H.VRe[i][i];
            }
            if(totmin+totnow<10*DBL_MIN)
            {
              totmin+=10*DBL_MIN;
              totnow+=10*DBL_MIN;
            }
            G->V[0].Re = (totmin * G->V[0].Re + totnow * val[ival]) / (totmin + totnow);
            for(i=0;i<A.NFermion;i++)
            {
              for(ivec=0;ivec<Nvec;ivec++)
              {
                for(j=0;j<A.NFermion;j++)
                {
                  for(jvec=0;jvec<Nvec;jvec++)
                  {
                    dvec[0] = (vec[3*ivec+0]-vec[3*jvec+0]) * A.Cell.rUnitCell[0][0] + (vec[3*ivec+1]-vec[3*jvec+1]) * A.Cell.rUnitCell[1][0] + (vec[3*ivec+2]-vec[3*jvec+2]) * A.Cell.rUnitCell[2][0];
                    dvec[1] = (vec[3*ivec+0]-vec[3*jvec+0]) * A.Cell.rUnitCell[0][1] + (vec[3*ivec+1]-vec[3*jvec+1]) * A.Cell.rUnitCell[1][1] + (vec[3*ivec+2]-vec[3*jvec+2]) * A.Cell.rUnitCell[2][1];
                    dvec[2] = (vec[3*ivec+0]-vec[3*jvec+0]) * A.Cell.rUnitCell[0][2] + (vec[3*ivec+1]-vec[3*jvec+1]) * A.Cell.rUnitCell[1][2] + (vec[3*ivec+2]-vec[3*jvec+2]) * A.Cell.rUnitCell[2][2];
                    phase.Re=cos(piornot * Dot(dvec,q));
                    phase.Im=sin(piornot * Dot(dvec,q));
                    G->V[1 + j+jvec*A.NFermion + (i+ivec*A.NFermion)*G->BlockSize].Re += phase.Re*H.VRe[i][j];
                    G->V[1 + j+jvec*A.NFermion + (i+ivec*A.NFermion)*G->BlockSize].Im += phase.Im*H.VRe[i][j];
                    if(H.QComplex)
                    {
                      G->V[1 + j+jvec*A.NFermion + (i+ivec*A.NFermion)*G->BlockSize].Re -= phase.Im*H.VIm[i][j];
                      G->V[1 + j+jvec*A.NFermion + (i+ivec*A.NFermion)*G->BlockSize].Im += phase.Re*H.VIm[i][j];
                    }
                  }
                }
              }
            }
            totmin+=totnow;
          }
          else
          {
            rest -= iE;  /* rest is a number between 0 and 1 pole at iE+1 gets 1-rest pole at iE+2 gest rest */
            for(ivec=0;ivec<Nvec;ivec++)
            {
              for(jvec=0;jvec<Nvec;jvec++)
              {
                dvec[0] = (vec[3*ivec+0]-vec[3*jvec+0]) * A.Cell.rUnitCell[0][0] + (vec[3*ivec+1]-vec[3*jvec+1]) * A.Cell.rUnitCell[1][0] + (vec[3*ivec+2]-vec[3*jvec+2]) * A.Cell.rUnitCell[2][0];
                dvec[1] = (vec[3*ivec+0]-vec[3*jvec+0]) * A.Cell.rUnitCell[0][1] + (vec[3*ivec+1]-vec[3*jvec+1]) * A.Cell.rUnitCell[1][1] + (vec[3*ivec+2]-vec[3*jvec+2]) * A.Cell.rUnitCell[2][1];
                dvec[2] = (vec[3*ivec+0]-vec[3*jvec+0]) * A.Cell.rUnitCell[0][2] + (vec[3*ivec+1]-vec[3*jvec+1]) * A.Cell.rUnitCell[1][2] + (vec[3*ivec+2]-vec[3*jvec+2]) * A.Cell.rUnitCell[2][2];
                phase.Re=cos(piornot * Dot(dvec,q));
                phase.Im=sin(piornot * Dot(dvec,q));

               //printf("phase for ivec %i jvec %i is %12.5E + I %12.5E at q is %12.5E %12.5E %12.5E and vec is  %12.5E %12.5E %12.5E q.dvec = %12.5E piornot = %12.5E\n",ivec, jvec, phase.Re, phase.Im, q[0], q[1], q[2], dvec[0], dvec[1], dvec[2], Dot(dvec,q),piornot);

                for(i=0;i<A.NFermion;i++)
                {
                  for(j=0;j<A.NFermion;j++)
                  {
                    G->V[(iE+1)*(1+G->BlockSize*G->BlockSize)+1 + j+jvec*A.NFermion + (i+ivec*A.NFermion)*G->BlockSize].Re += phase.Re*(1-rest)*H.VRe[i][j];
                    G->V[(iE+1)*(1+G->BlockSize*G->BlockSize)+1 + j+jvec*A.NFermion + (i+ivec*A.NFermion)*G->BlockSize].Im += phase.Im*(1-rest)*H.VRe[i][j];
                    G->V[(iE+2)*(1+G->BlockSize*G->BlockSize)+1 + j+jvec*A.NFermion + (i+ivec*A.NFermion)*G->BlockSize].Re += phase.Re*(  rest)*H.VRe[i][j];
                    G->V[(iE+2)*(1+G->BlockSize*G->BlockSize)+1 + j+jvec*A.NFermion + (i+ivec*A.NFermion)*G->BlockSize].Im += phase.Im*(  rest)*H.VRe[i][j];
                    if(H.QComplex)
                    {
                      G->V[(iE+1)*(1+G->BlockSize*G->BlockSize)+1 + j+jvec*A.NFermion + (i+ivec*A.NFermion)*G->BlockSize].Re -= phase.Im*(1-rest)*H.VIm[i][j];
                      G->V[(iE+1)*(1+G->BlockSize*G->BlockSize)+1 + j+jvec*A.NFermion + (i+ivec*A.NFermion)*G->BlockSize].Im += phase.Re*(1-rest)*H.VIm[i][j];
                      G->V[(iE+2)*(1+G->BlockSize*G->BlockSize)+1 + j+jvec*A.NFermion + (i+ivec*A.NFermion)*G->BlockSize].Re -= phase.Im*(  rest)*H.VIm[i][j];
                      G->V[(iE+2)*(1+G->BlockSize*G->BlockSize)+1 + j+jvec*A.NFermion + (i+ivec*A.NFermion)*G->BlockSize].Im += phase.Re*(  rest)*H.VIm[i][j];
                    }
                  }
                }
              }
            }
          }
        }
      }
      FreeMatrix(&H);
      FreeCompactMatrix(&T);
    }
  }
  Nktot=Nk[0]*Nk[1]*Nk[2];
  for(i=0;i<G->N*(1+G->BlockSize*G->BlockSize);i+=(1+G->BlockSize*G->BlockSize))
  {
    for(j=1;j<1+G->BlockSize*G->BlockSize;j++)
    {
      G->V[i+j].Re/=Nktot;
      G->V[i+j].Im/=Nktot;
    }
  }
  ChopBlockListOfPoles(G, DBL_EPSILON);
  return 0;
}


/*  Given an tight binding Hamiltonian and a k-mesh defined by the number of divisions of a, b and c (Nk)
 *  returns the green function matrix of size A.NFermion by A.NFermion
 *          on an energy mesh with N+3 points. N+1 points from Emin to Emax and one point below Emin and one above Emax such that all weight and first moment is ok
 *  Calculates the eigen values and functions for each k including a self energy input as a Block List of poles
 *  for each k and eigen value i
 *  weight is added to the pole at iE and iE+1 such that Energy[iE] < val < Energy[iE+1]
 *  The matrix h_{jk} is fun[i,j] * fun[i,k]^* is added
 *
 *  Sigma must have dimensions BlockSize equal to A.NFermion
 */
uint32_t TightBindingHamiltonianPlusSelfEnergyToBlockListOfPoles(TightBindingDef const A, int Nk[3], BlockListOfPoles const Sigma, BlockListOfPoles * const G, double const Emin, double const Emax, uint32_t const NE)
{
  MatrixType H;
  CompactMatrixType T, HSigma, HSigmaComplex;
  double q[3], val[A.NFermion*(Sigma.N+1)], rest, totmin=0, totmax=0, totnow;
  int    iq[3],ival,iE,i,j;
  double diq[3],dE=(Emax-Emin)/((double)NE),energy,Nktot;
  // test if sigma size is ok
  if(Sigma.BlockSize!=A.NFermion)
  {
    printf("in TightBindingHamiltonianPlusSelfEnergyToBlockListOfPoles\nThe dimension of the self energy (%i) is not equal to the dimension of the tight binding Hamiltonian (%i)\n",Sigma.BlockSize,A.NFermion);
    return 1;
  }
  // init G
  G->BlockSize=A.NFermion;
  G->N=NE+3;
  G->QComplex=1;
  strcpy(G->Name, A.Name);
  if(InitBlockListOfPoles(G))
  {
    printf("InitBlockListOfPoles failed in TightBindingHamiltonianToBlockListOfPoles\n");
    return 1;
  }
  // set energy of G
  energy=Emin-dE;
  for(i=0;i<G->N*(1+G->BlockSize*G->BlockSize);i+=(1+G->BlockSize*G->BlockSize))
  {
    G->V[i].Re=energy;
    energy+=dE;
  }
  // create HSigma
  BlockListOfPolesInversedToCompactMatrix(&Sigma, &HSigma);
  if(Sigma.QComplex==0)
  {
    CompactMatrixCopy(HSigma, &HSigmaComplex);
    MakeCompactRealMatrixComplex(&HSigmaComplex);
  }
  for (iq[0] = 0; iq[0] < Nk[0]; iq[0]++)
  {
    for (iq[1] = 0; iq[1] < Nk[1]; iq[1]++)
    {
      for (iq[2] = 0; iq[2] < Nk[2]; iq[2]++)
      {
        diq[0] = iq[0] / ((double)(Nk[0]));
        diq[1] = iq[1] / ((double)(Nk[1]));
        diq[2] = iq[2] / ((double)(Nk[2]));
        q[0] = diq[0] * A.Cell.gUnitCell[0][0] + diq[1] * A.Cell.gUnitCell[1][0] + diq[2] * A.Cell.gUnitCell[2][0];
        q[1] = diq[0] * A.Cell.gUnitCell[0][1] + diq[1] * A.Cell.gUnitCell[1][1] + diq[2] * A.Cell.gUnitCell[2][1];
        q[2] = diq[0] * A.Cell.gUnitCell[0][2] + diq[1] * A.Cell.gUnitCell[1][2] + diq[2] * A.Cell.gUnitCell[2][2];
        CreateHamiltonianFromTightBindingDefConstantPhaseInUnitCell(A, q, &H);
        if(H.QComplex)
        {
          if(Sigma.QComplex)
          {
            for(i=0;i<A.NFermion;i++)
            {
              for(j=0;j<A.NFermion;j++)
              {
                HSigma.V[i*HSigma.Nj+j].Re=Sigma.a0[i*Sigma.BlockSize+j].Re + H.VRe[i][j];
                HSigma.V[i*HSigma.Nj+j].Im=Sigma.a0[i*Sigma.BlockSize+j].Im + H.VIm[i][j];
              }
            }
            CompactMatrixEigensystem(&HSigma, &T, val);
          }
          else
          {
            for(i=0;i<A.NFermion;i++)
            {
              for(j=0;j<A.NFermion;j++)
              {
                HSigmaComplex.V[i*HSigma.Nj+j].Re=Sigma.a0Re[i*Sigma.BlockSize+j] + H.VRe[i][j];
                HSigmaComplex.V[i*HSigma.Nj+j].Im=                                  H.VIm[i][j];
              }
            }
            CompactMatrixEigensystem(&HSigmaComplex, &T, val);
          }
        }
        else
        {
          if(Sigma.QComplex)
          {
            for(i=0;i<A.NFermion;i++)
            {
              for(j=0;j<A.NFermion;j++)
              {
                HSigma.V[i*HSigma.Nj+j].Re=Sigma.a0[i*Sigma.BlockSize+j].Re + H.VRe[i][j];
                HSigma.V[i*HSigma.Nj+j].Im=Sigma.a0[i*Sigma.BlockSize+j].Im;
              }
            }
            CompactMatrixEigensystem(&HSigma, &T, val);
          }
          else
          {
            for(i=0;i<A.NFermion;i++)
            {
              for(j=0;j<A.NFermion;j++)
              {
                HSigma.VRe[i*HSigma.Nj+j]=Sigma.a0Re[i*Sigma.BlockSize+j] + H.VRe[i][j];
              }
            }
            CompactMatrixEigensystem(&HSigma, &T, val);
          }
        }
        if(T.QComplex)
        {
          if(H.QComplex==0) MakeRealMatrixComplex(&H);
        }
        else
        {
          if(H.QComplex) MakeComplexMatrixReal(&H);
        }
        for(ival=0;ival<T.Nj;ival++)
        {
          /* calculate matrix B for eigen-value ival */
          if(H.QComplex)
          {
            ComplexEigenFunctionToCharacter(&(T.V[ival*T.Nj]), &H);
          }
          else
          {
            RealEigenFunctionToCharacter(&(T.VRe[ival*T.Nj]), &H);
          }
          /* determine the position of the pole in the list of poles */
          rest = (val[ival] - Emin) / dE;
          iE = (int)(floor(rest));
          if (iE >= NE)
          {
            totnow=0;
            for(i=0;i<A.NFermion;i++)
            {
              totnow+=H.VRe[i][i];
            }
            if(totmax+totnow<10*DBL_MIN)
            {
              totmax+=10*DBL_MIN;
              totnow+=10*DBL_MIN;
            }
            G->V[(NE+2)*(1+G->BlockSize*G->BlockSize)].Re = (totmax * G->V[(NE+2)*(1+G->BlockSize*G->BlockSize)].Re + totnow * val[ival]) / (totmax + totnow);
            for(i=0;i<A.NFermion;i++)
            {
              for(j=0;j<A.NFermion;j++)
              {
                G->V[(NE+2)*(1+G->BlockSize*G->BlockSize)+1+j+i*A.NFermion].Re+=H.VRe[i][j];
                if(H.QComplex)
                {
                  G->V[(NE+2)*(1+G->BlockSize*G->BlockSize)+1+j+i*A.NFermion].Im+=H.VIm[i][j];
                }
              }
            }
            totmax+=totnow;
          }
          else if (iE < 0)
          {
            totnow=0;
            for(i=0;i<A.NFermion;i++)
            {
              totnow+=H.VRe[i][i];
            }
            if(totmin+totnow<10*DBL_MIN)
            {
              totmin+=10*DBL_MIN;
              totnow+=10*DBL_MIN;
            }
            G->V[0].Re = (totmin * G->V[0].Re + totnow * val[ival]) / (totmin + totnow);
            for(i=0;i<A.NFermion;i++)
            {
              for(j=0;j<A.NFermion;j++)
              {
                G->V[1+j+i*A.NFermion].Re+=H.VRe[i][j];
                if(H.QComplex)
                {
                  G->V[1+j+i*A.NFermion].Im+=H.VIm[i][j];
                }
              }
            }
            totmin+=totnow;
          }
          else
          {
            rest -= iE;  /* rest is a number between 0 and 1 pole at iE+1 gets 1-rest pole at iE+2 gest rest */
            for(i=0;i<A.NFermion;i++)
            {
              for(j=0;j<A.NFermion;j++)
              {
                G->V[(iE+1)*(1+G->BlockSize*G->BlockSize)+1+j+i*A.NFermion].Re+=(1-rest)*H.VRe[i][j];
                G->V[(iE+2)*(1+G->BlockSize*G->BlockSize)+1+j+i*A.NFermion].Re+=(  rest)*H.VRe[i][j];
                if(H.QComplex)
                {
                  G->V[(iE+1)*(1+G->BlockSize*G->BlockSize)+1+j+i*A.NFermion].Im+=(1-rest)*H.VIm[i][j];
                  G->V[(iE+2)*(1+G->BlockSize*G->BlockSize)+1+j+i*A.NFermion].Im+=(  rest)*H.VIm[i][j];
                }
              }
            }
          }
        }
      }
      FreeMatrix(&H);
      FreeCompactMatrix(&T);
    }
  }
  Nktot=Nk[0]*Nk[1]*Nk[2];
  for(i=0;i<G->N*(1+G->BlockSize*G->BlockSize);i+=(1+G->BlockSize*G->BlockSize))
  {
    for(j=1;j<1+G->BlockSize*G->BlockSize;j++)
    {
      G->V[i+j].Re/=Nktot;
      G->V[i+j].Im/=Nktot;
    }
  }
  ChopBlockListOfPoles(G, DBL_EPSILON);
  return 0;
}


/*  Given an tight binding Hamiltonian and a k-mesh defined by the number of divisions of a, b and c (Nk)
 *  returns the green function matrix of size Nvec*A.NFermion by Nvec*A.NFermion
 *          on an energy mesh with N+3 points. N+1 points from Emin to Emax and one point below Emin and one above Emax such that all weight and first moment is ok
 *  Calculates the eigen values and functions for each k including a self energy input as a Block List of poles
 *  for each k and eigen value i
 *  weight is added to the pole at iE and iE+1 such that Energy[iE] < val < Energy[iE+1]
 *  The matrix h_{jk} is fun[i,j] * fun[i,k]^* is added
 *
 *  Sigma.Blocksize must be equal to A.NFermion
 *    All sites have the same self energy !!!
 */
uint32_t TightBindingHamiltonianPlusSelfEnergyMultiSiteToBlockListOfPoles(TightBindingDef const A, int Nk[3], uint32_t const Nvec, int const *vec, BlockListOfPoles const Sigma, BlockListOfPoles * const G, double const Emin, double const Emax, uint32_t const NE)
{
  MatrixType H;
  CompactMatrixType T, HSigma, HSigmaComplex;
  double q[3], val[A.NFermion*(Sigma.N+1)], rest, totmin=0, totmax=0, totnow, dvec[3];
  double const piornot = ((A.Cell.Units&1) ? (1) : (Pi));  // if A.Cell.Units last bit is 1 then r.g = pi, otherwise it evaluates to 1, piornot * r.g now eveluates always to pi.
  int    iq[3],ival,iE,i,j,ivec,jvec;
  double diq[3],dE=(Emax-Emin)/((double)NE),energy,Nktot;
  Complex16 phase;
  // test if sigma size is ok
  if(Sigma.BlockSize!=A.NFermion)
  {
    printf("in TightBindingHamiltonianPlusSelfEnergyToBlockListOfPoles\nThe dimension of the self energy (%i) is not equal to the dimension of the tight binding Hamiltonian (%i) \n",Sigma.BlockSize,A.NFermion);
    return 1;
  }
  // init G
  G->BlockSize=A.NFermion * Nvec;
  G->N=NE+3;
  G->QComplex=1;
  strcpy(G->Name, A.Name);
  if(InitBlockListOfPoles(G))
  {
    printf("InitBlockListOfPoles failed in TightBindingHamiltonianToBlockListOfPoles\n");
    return 1;
  }
  // set energy of G
  energy=Emin-dE;
  for(i=0;i<G->N*(1+G->BlockSize*G->BlockSize);i+=(1+G->BlockSize*G->BlockSize))
  {
    G->V[i].Re=energy;
    energy+=dE;
  }
  // create HSigma
  BlockListOfPolesInversedToCompactMatrix(&Sigma, &HSigma);
  if(Sigma.QComplex==0)
  {
    CompactMatrixCopy(HSigma, &HSigmaComplex);
    MakeCompactRealMatrixComplex(&HSigmaComplex);
  }
  for (iq[0] = 0; iq[0] < Nk[0]; iq[0]++)
  {
    for (iq[1] = 0; iq[1] < Nk[1]; iq[1]++)
    {
      for (iq[2] = 0; iq[2] < Nk[2]; iq[2]++)
      {
        diq[0] = iq[0] / ((double)(Nk[0]));
        diq[1] = iq[1] / ((double)(Nk[1]));
        diq[2] = iq[2] / ((double)(Nk[2]));
        q[0] = diq[0] * A.Cell.gUnitCell[0][0] + diq[1] * A.Cell.gUnitCell[1][0] + diq[2] * A.Cell.gUnitCell[2][0];
        q[1] = diq[0] * A.Cell.gUnitCell[0][1] + diq[1] * A.Cell.gUnitCell[1][1] + diq[2] * A.Cell.gUnitCell[2][1];
        q[2] = diq[0] * A.Cell.gUnitCell[0][2] + diq[1] * A.Cell.gUnitCell[1][2] + diq[2] * A.Cell.gUnitCell[2][2];
        CreateHamiltonianFromTightBindingDefConstantPhaseInUnitCell(A, q, &H);
        if(H.QComplex)
        {
          if(Sigma.QComplex)
          {
            for(i=0;i<A.NFermion;i++)
            {
              for(j=0;j<A.NFermion;j++)
              {
                HSigma.V[i*HSigma.Nj+j].Re=Sigma.a0[i*Sigma.BlockSize+j].Re + H.VRe[i][j];
                HSigma.V[i*HSigma.Nj+j].Im=Sigma.a0[i*Sigma.BlockSize+j].Im + H.VIm[i][j];
              }
            }
            CompactMatrixEigensystem(&HSigma, &T, val);
          }
          else
          {
            for(i=0;i<A.NFermion;i++)
            {
              for(j=0;j<A.NFermion;j++)
              {
                HSigmaComplex.V[i*HSigma.Nj+j].Re=Sigma.a0Re[i*Sigma.BlockSize+j] + H.VRe[i][j];
                HSigmaComplex.V[i*HSigma.Nj+j].Im=                                  H.VIm[i][j];
              }
            }
            CompactMatrixEigensystem(&HSigmaComplex, &T, val);
          }
        }
        else
        {
          if(Sigma.QComplex)
          {
            for(i=0;i<A.NFermion;i++)
            {
              for(j=0;j<A.NFermion;j++)
              {
                HSigma.V[i*HSigma.Nj+j].Re=Sigma.a0[i*Sigma.BlockSize+j].Re + H.VRe[i][j];
                HSigma.V[i*HSigma.Nj+j].Im=Sigma.a0[i*Sigma.BlockSize+j].Im;
              }
            }
            CompactMatrixEigensystem(&HSigma, &T, val);
          }
          else
          {
            for(i=0;i<A.NFermion;i++)
            {
              for(j=0;j<A.NFermion;j++)
              {
                HSigma.VRe[i*HSigma.Nj+j]=Sigma.a0Re[i*Sigma.BlockSize+j] + H.VRe[i][j];
              }
            }
            CompactMatrixEigensystem(&HSigma, &T, val);
          }
        }
        if(T.QComplex)
        {
          if(H.QComplex==0) MakeRealMatrixComplex(&H);
        }
        else
        {
          if(H.QComplex) MakeComplexMatrixReal(&H);
        }
        for(ival=0;ival<T.Nj;ival++)
        {
          /* calculate character matrix for eigen-value ival - note that H changes meaning from Hamiltonian to Character matrix */
          if(H.QComplex)
          {
            ComplexEigenFunctionToCharacter(&(T.V[ival*T.Nj]), &H);
          }
          else
          {
            RealEigenFunctionToCharacter(&(T.VRe[ival*T.Nj]), &H);
          }
          /* determine the position of the pole in the list of poles */
          rest = (val[ival] - Emin) / dE;
          iE = (int)(floor(rest));
          if (iE >= NE)
          {
            totnow=0;
            for(i=0;i<A.NFermion;i++)
            {
              totnow+=H.VRe[i][i];
            }
            if(totmax+totnow<10*DBL_MIN)
            {
              totmax+=10*DBL_MIN;
              totnow+=10*DBL_MIN;
            }
            G->V[(NE+2)*(1+G->BlockSize*G->BlockSize)].Re = (totmax * G->V[(NE+2)*(1+G->BlockSize*G->BlockSize)].Re + totnow * val[ival]) / (totmax + totnow);
            for(i=0;i<A.NFermion;i++)
            {
              for(ivec=0;ivec<Nvec;ivec++)
              {
                for(j=0;j<A.NFermion;j++)
                {
                  for(jvec=0;jvec<Nvec;jvec++)
                  {
                    dvec[0] = (vec[3*ivec+0]-vec[3*jvec+0]) * A.Cell.rUnitCell[0][0] + (vec[3*ivec+1]-vec[3*jvec+1]) * A.Cell.rUnitCell[1][0] + (vec[3*ivec+2]-vec[3*jvec+2]) * A.Cell.rUnitCell[2][0];
                    dvec[1] = (vec[3*ivec+0]-vec[3*jvec+0]) * A.Cell.rUnitCell[0][1] + (vec[3*ivec+1]-vec[3*jvec+1]) * A.Cell.rUnitCell[1][1] + (vec[3*ivec+2]-vec[3*jvec+2]) * A.Cell.rUnitCell[2][1];
                    dvec[2] = (vec[3*ivec+0]-vec[3*jvec+0]) * A.Cell.rUnitCell[0][2] + (vec[3*ivec+1]-vec[3*jvec+1]) * A.Cell.rUnitCell[1][2] + (vec[3*ivec+2]-vec[3*jvec+2]) * A.Cell.rUnitCell[2][2];
                    phase.Re=cos(piornot * Dot(dvec,q));
                    phase.Im=sin(piornot * Dot(dvec,q));
                    G->V[(NE+2)*(1+G->BlockSize*G->BlockSize)+1 + j+jvec*A.NFermion + (i+ivec*A.NFermion)*G->BlockSize].Re += phase.Re*H.VRe[i][j];
                    G->V[(NE+2)*(1+G->BlockSize*G->BlockSize)+1 + j+jvec*A.NFermion + (i+ivec*A.NFermion)*G->BlockSize].Im += phase.Im*H.VRe[i][j];
                    if(H.QComplex)
                    {
                      G->V[(NE+2)*(1+G->BlockSize*G->BlockSize)+1 + j+jvec*A.NFermion + (i+ivec*A.NFermion)*G->BlockSize].Re -= phase.Im*H.VIm[i][j];
                      G->V[(NE+2)*(1+G->BlockSize*G->BlockSize)+1 + j+jvec*A.NFermion + (i+ivec*A.NFermion)*G->BlockSize].Im += phase.Re*H.VIm[i][j];
                    }
                  }
                }
              }
            }
            totmax+=totnow;
          }
          else if (iE < 0)
          {
            totnow=0;
            for(i=0;i<A.NFermion;i++)
            {
              totnow+=H.VRe[i][i];
            }
            if(totmin+totnow<10*DBL_MIN)
            {
              totmin+=10*DBL_MIN;
              totnow+=10*DBL_MIN;
            }
            G->V[0].Re = (totmin * G->V[0].Re + totnow * val[ival]) / (totmin + totnow);
            for(i=0;i<A.NFermion;i++)
            {
              for(ivec=0;ivec<Nvec;ivec++)
              {
                for(j=0;j<A.NFermion;j++)
                {
                  for(jvec=0;jvec<Nvec;jvec++)
                  {
                    dvec[0] = (vec[3*ivec+0]-vec[3*jvec+0]) * A.Cell.rUnitCell[0][0] + (vec[3*ivec+1]-vec[3*jvec+1]) * A.Cell.rUnitCell[1][0] + (vec[3*ivec+2]-vec[3*jvec+2]) * A.Cell.rUnitCell[2][0];
                    dvec[1] = (vec[3*ivec+0]-vec[3*jvec+0]) * A.Cell.rUnitCell[0][1] + (vec[3*ivec+1]-vec[3*jvec+1]) * A.Cell.rUnitCell[1][1] + (vec[3*ivec+2]-vec[3*jvec+2]) * A.Cell.rUnitCell[2][1];
                    dvec[2] = (vec[3*ivec+0]-vec[3*jvec+0]) * A.Cell.rUnitCell[0][2] + (vec[3*ivec+1]-vec[3*jvec+1]) * A.Cell.rUnitCell[1][2] + (vec[3*ivec+2]-vec[3*jvec+2]) * A.Cell.rUnitCell[2][2];
                    phase.Re=cos(piornot * Dot(dvec,q));
                    phase.Im=sin(piornot * Dot(dvec,q));
                    G->V[1 + j+jvec*A.NFermion + (i+ivec*A.NFermion)*G->BlockSize].Re += phase.Re*H.VRe[i][j];
                    G->V[1 + j+jvec*A.NFermion + (i+ivec*A.NFermion)*G->BlockSize].Im += phase.Im*H.VRe[i][j];
                    if(H.QComplex)
                    {
                      G->V[1 + j+jvec*A.NFermion + (i+ivec*A.NFermion)*G->BlockSize].Re -= phase.Im*H.VIm[i][j];
                      G->V[1 + j+jvec*A.NFermion + (i+ivec*A.NFermion)*G->BlockSize].Im += phase.Re*H.VIm[i][j];
                    }
                  }
                }
              }
            }
            totmin+=totnow;
          }
          else
          {
            rest -= iE;  /* rest is a number between 0 and 1 pole at iE+1 gets 1-rest pole at iE+2 gest rest */
            for(ivec=0;ivec<Nvec;ivec++)
            {
              for(jvec=0;jvec<Nvec;jvec++)
              {
                dvec[0] = (vec[3*ivec+0]-vec[3*jvec+0]) * A.Cell.rUnitCell[0][0] + (vec[3*ivec+1]-vec[3*jvec+1]) * A.Cell.rUnitCell[1][0] + (vec[3*ivec+2]-vec[3*jvec+2]) * A.Cell.rUnitCell[2][0];
                dvec[1] = (vec[3*ivec+0]-vec[3*jvec+0]) * A.Cell.rUnitCell[0][1] + (vec[3*ivec+1]-vec[3*jvec+1]) * A.Cell.rUnitCell[1][1] + (vec[3*ivec+2]-vec[3*jvec+2]) * A.Cell.rUnitCell[2][1];
                dvec[2] = (vec[3*ivec+0]-vec[3*jvec+0]) * A.Cell.rUnitCell[0][2] + (vec[3*ivec+1]-vec[3*jvec+1]) * A.Cell.rUnitCell[1][2] + (vec[3*ivec+2]-vec[3*jvec+2]) * A.Cell.rUnitCell[2][2];
                phase.Re=cos(piornot * Dot(dvec,q));
                phase.Im=sin(piornot * Dot(dvec,q));
                for(i=0;i<A.NFermion;i++)
                {
                  for(j=0;j<A.NFermion;j++)
                  {
                    G->V[(iE+1)*(1+G->BlockSize*G->BlockSize)+1 + j+jvec*A.NFermion + (i+ivec*A.NFermion)*G->BlockSize].Re += phase.Re*(1-rest)*H.VRe[i][j];
                    G->V[(iE+1)*(1+G->BlockSize*G->BlockSize)+1 + j+jvec*A.NFermion + (i+ivec*A.NFermion)*G->BlockSize].Im += phase.Im*(1-rest)*H.VRe[i][j];
                    G->V[(iE+2)*(1+G->BlockSize*G->BlockSize)+1 + j+jvec*A.NFermion + (i+ivec*A.NFermion)*G->BlockSize].Re += phase.Re*(  rest)*H.VRe[i][j];
                    G->V[(iE+2)*(1+G->BlockSize*G->BlockSize)+1 + j+jvec*A.NFermion + (i+ivec*A.NFermion)*G->BlockSize].Im += phase.Im*(  rest)*H.VRe[i][j];
                    if(H.QComplex)
                    {
                      G->V[(iE+1)*(1+G->BlockSize*G->BlockSize)+1 + j+jvec*A.NFermion + (i+ivec*A.NFermion)*G->BlockSize].Re -= phase.Im*(1-rest)*H.VIm[i][j];
                      G->V[(iE+1)*(1+G->BlockSize*G->BlockSize)+1 + j+jvec*A.NFermion + (i+ivec*A.NFermion)*G->BlockSize].Im += phase.Re*(1-rest)*H.VIm[i][j];
                      G->V[(iE+2)*(1+G->BlockSize*G->BlockSize)+1 + j+jvec*A.NFermion + (i+ivec*A.NFermion)*G->BlockSize].Re -= phase.Im*(  rest)*H.VIm[i][j];
                      G->V[(iE+2)*(1+G->BlockSize*G->BlockSize)+1 + j+jvec*A.NFermion + (i+ivec*A.NFermion)*G->BlockSize].Im += phase.Re*(  rest)*H.VIm[i][j];
                    }
                  }
                }
              }
            }
          }
        }
      }
      FreeMatrix(&H);
      FreeCompactMatrix(&T);
    }
  }
  Nktot=Nk[0]*Nk[1]*Nk[2];
  for(i=0;i<G->N*(1+G->BlockSize*G->BlockSize);i+=(1+G->BlockSize*G->BlockSize))
  {
    for(j=1;j<1+G->BlockSize*G->BlockSize;j++)
    {
      G->V[i+j].Re/=Nktot;
      G->V[i+j].Im/=Nktot;
    }
  }
  ChopBlockListOfPoles(G, DBL_EPSILON);
  return 0;
}




/*  Given an tight binding Hamiltonian and a k-mesh defined by the number of divisions of a, b and c (Nk)
 *  returns the green function matrix of size 2*A.NFermion by 2*A.NFermion
 *          on an energy mesh with N+3 points. N+1 points from Emin to Emax and one point below Emin and one above Emax such that all weight and first moment is ok
 *  Calculates the eigen values and functions for each k including a self energy input as a Block List of poles
 *  for each k and eigen value i
 *  weight is added to the pole at iE and iE+1 such that Energy[iE] < val < Energy[iE+1]
 *  The matrix h_{jk} is fun[i,j] * fun[i,k]^* is added
 *
 *  Sigma must have dimensions BlockSize equal to 2 * A.NFermion
 */
uint32_t TightBindingHamiltonianAddSpinPlusSelfEnergyToBlockListOfPoles(TightBindingDef const A, int Nk[3], BlockListOfPoles const Sigma, BlockListOfPoles * const G, double const Emin, double const Emax, uint32_t const NE)
{
  MatrixType H, PsiChar;
  CompactMatrixType T, HSigma, HSigmaComplex;
  double q[3], val[2*A.NFermion*(Sigma.N+1)], rest, totmin=0, totmax=0, totnow;
  int    iq[3],ival,iE,i,j;
  double diq[3],dE=(Emax-Emin)/((double)NE),energy,Nktot;
  // initiallize PsiChar, a matrix storing the density of a wave-function
  PsiChar.Ni=2*A.NFermion;
  PsiChar.Nj=2*A.NFermion;
  PsiChar.QComplex=1;
  InitMatrix(&PsiChar);
  // test if sigma size is ok
  if(Sigma.BlockSize!= 2 * A.NFermion)
  {
    printf("in TightBindingHamiltonianAddSpinPlusSelfEnergyToBlockListOfPoles\nThe dimension of the self energy (%i) is not equal to the dimension of the tight binding Hamiltonian (%i) times 2 for the spin\n",Sigma.BlockSize,A.NFermion);
    return 1;
  }
  // init G
  G->BlockSize=2*A.NFermion;
  G->N=NE+3;
  G->QComplex=1;
  strcpy(G->Name, A.Name);
  if(InitBlockListOfPoles(G))
  {
    printf("InitBlockListOfPoles failed in TightBindingHamiltonianAddSpinPlusSelfEnergyToBlockListOfPoles\n");
    return 1;
  }
  // set energy of G
  energy=Emin-dE;
  for(i=0;i<G->N*(1+G->BlockSize*G->BlockSize);i+=(1+G->BlockSize*G->BlockSize))
  {
    G->V[i].Re=energy;
    energy+=dE;
  }
  // create HSigma
  BlockListOfPolesInversedToCompactMatrix(&Sigma, &HSigma);
  if(Sigma.QComplex==0)
  {
    CompactMatrixCopy(HSigma, &HSigmaComplex);
    MakeCompactRealMatrixComplex(&HSigmaComplex);
  }
  for (iq[0] = 0; iq[0] < Nk[0]; iq[0]++)
  {
    for (iq[1] = 0; iq[1] < Nk[1]; iq[1]++)
    {
      for (iq[2] = 0; iq[2] < Nk[2]; iq[2]++)
      {
        diq[0] = iq[0] / ((double)(Nk[0]));
        diq[1] = iq[1] / ((double)(Nk[1]));
        diq[2] = iq[2] / ((double)(Nk[2]));
        q[0] = diq[0] * A.Cell.gUnitCell[0][0] + diq[1] * A.Cell.gUnitCell[1][0] + diq[2] * A.Cell.gUnitCell[2][0];
        q[1] = diq[0] * A.Cell.gUnitCell[0][1] + diq[1] * A.Cell.gUnitCell[1][1] + diq[2] * A.Cell.gUnitCell[2][1];
        q[2] = diq[0] * A.Cell.gUnitCell[0][2] + diq[1] * A.Cell.gUnitCell[1][2] + diq[2] * A.Cell.gUnitCell[2][2];
        CreateHamiltonianFromTightBindingDefConstantPhaseInUnitCell(A, q, &H);
        if(H.QComplex)
        {
          if(Sigma.QComplex)
          {
            for(i=0;i<A.NFermion;i++)
            {
              for(j=0;j<A.NFermion;j++)
              {
                HSigma.V[i*HSigma.Nj+j].Re=Sigma.a0[i*Sigma.BlockSize+j].Re + H.VRe[i][j];
                HSigma.V[i*HSigma.Nj+j].Im=Sigma.a0[i*Sigma.BlockSize+j].Im + H.VIm[i][j];
                HSigma.V[(i+A.NFermion)*HSigma.Nj+(j+A.NFermion)].Re=Sigma.a0[(i+A.NFermion)*Sigma.BlockSize+(j+A.NFermion)].Re + H.VRe[i][j];
                HSigma.V[(i+A.NFermion)*HSigma.Nj+(j+A.NFermion)].Im=Sigma.a0[(i+A.NFermion)*Sigma.BlockSize+(j+A.NFermion)].Im + H.VIm[i][j];
              }
            }
            CompactMatrixEigensystem(&HSigma, &T, val);
          }
          else
          {
            for(i=0;i<A.NFermion;i++)
            {
              for(j=0;j<A.NFermion;j++)
              {
                HSigmaComplex.V[i*HSigma.Nj+j].Re=Sigma.a0Re[i*Sigma.BlockSize+j] + H.VRe[i][j];
                HSigmaComplex.V[i*HSigma.Nj+j].Im=                                  H.VIm[i][j];
                HSigmaComplex.V[(i+A.NFermion)*HSigma.Nj+(j+A.NFermion)].Re=Sigma.a0Re[(i+A.NFermion)*Sigma.BlockSize+(j+A.NFermion)] + H.VRe[i][j];
                HSigmaComplex.V[(i+A.NFermion)*HSigma.Nj+(j+A.NFermion)].Im=                                                            H.VIm[i][j];
              }
            }
            CompactMatrixEigensystem(&HSigmaComplex, &T, val);
          }
        }
        else
        {
          if(Sigma.QComplex)
          {
            for(i=0;i<A.NFermion;i++)
            {
              for(j=0;j<A.NFermion;j++)
              {
                HSigma.V[i*HSigma.Nj+j].Re=Sigma.a0[i*Sigma.BlockSize+j].Re + H.VRe[i][j];
                HSigma.V[i*HSigma.Nj+j].Im=Sigma.a0[i*Sigma.BlockSize+j].Im;
                HSigma.V[(i+A.NFermion)*HSigma.Nj+(j+A.NFermion)].Re=Sigma.a0[(i+A.NFermion)*Sigma.BlockSize+(j+A.NFermion)].Re + H.VRe[i][j];
                HSigma.V[(i+A.NFermion)*HSigma.Nj+(j+A.NFermion)].Im=Sigma.a0[(i+A.NFermion)*Sigma.BlockSize+(j+A.NFermion)].Im;
              }
            }
            CompactMatrixEigensystem(&HSigma, &T, val);
          }
          else
          {
            for(i=0;i<A.NFermion;i++)
            {
              for(j=0;j<A.NFermion;j++)
              {
                HSigma.VRe[i*HSigma.Nj+j]=Sigma.a0Re[i*Sigma.BlockSize+j] + H.VRe[i][j];
                HSigma.VRe[(i+A.NFermion)*HSigma.Nj+(j+A.NFermion)]=Sigma.a0Re[(i+A.NFermion)*Sigma.BlockSize+(j+A.NFermion)] + H.VRe[i][j];
              }
            }
            CompactMatrixEigensystem(&HSigma, &T, val);
          }
        }
        FreeMatrix(&H);
        for(ival=0;ival<T.Nj;ival++)
        {
          /* calculate matrix B for eigen-value ival */
          if(T.QComplex)
          {
            ComplexEigenFunctionToCharacter(&(T.V[ival*T.Nj]), &PsiChar);
          }
          else
          {
            RealEigenFunctionToCharacter(&(T.VRe[ival*T.Nj]), &PsiChar);
          }
          /* determine the position of the pole in the list of poles */
          rest = (val[ival] - Emin) / dE;
          iE = (int)(floor(rest));
          if (iE >= NE)
          {
            totnow=0;
            for(i=0;i<2*A.NFermion;i++)
            {
              totnow+=PsiChar.VRe[i][i];
            }
            if(totmax+totnow<10*DBL_MIN)
            {
              totmax+=10*DBL_MIN;
              totnow+=10*DBL_MIN;
            }
            G->V[(NE+2)*(1+G->BlockSize*G->BlockSize)].Re = (totmax * G->V[(NE+2)*(1+G->BlockSize*G->BlockSize)].Re + totnow * val[ival]) / (totmax + totnow);
            for(i=0;i<2*A.NFermion;i++)
            {
              for(j=0;j<2*A.NFermion;j++)
              {
                G->V[(NE+2)*(1+G->BlockSize*G->BlockSize)+1+j+i*2*A.NFermion].Re+=PsiChar.VRe[i][j];
                if(T.QComplex)
                {
                  G->V[(NE+2)*(1+G->BlockSize*G->BlockSize)+1+j+i*2*A.NFermion].Im+=PsiChar.VIm[i][j];
                }
              }
            }
            totmax+=totnow;
          }
          else if (iE < 0)
          {
            totnow=0;
            for(i=0;i<2*A.NFermion;i++)
            {
              totnow+=PsiChar.VRe[i][i];
            }
            if(totmin+totnow<10*DBL_MIN)
            {
              totmin+=10*DBL_MIN;
              totnow+=10*DBL_MIN;
            }
            G->V[0].Re = (totmin * G->V[0].Re + totnow * val[ival]) / (totmin + totnow);
            for(i=0;i<2*A.NFermion;i++)
            {
              for(j=0;j<2*A.NFermion;j++)
              {
                G->V[1+j+i*2*A.NFermion].Re+=PsiChar.VRe[i][j];
                if(T.QComplex)
                {
                  G->V[1+j+i*2*A.NFermion].Im+=PsiChar.VIm[i][j];
                }
              }
            }
            totmin+=totnow;
          }
          else
          {
            rest -= iE;  /* rest is a number between 0 and 1 pole at iE+1 gets 1-rest pole at iE+2 gest rest */
            for(i=0;i<2*A.NFermion;i++)
            {
              for(j=0;j<2*A.NFermion;j++)
              {
                G->V[(iE+1)*(1+G->BlockSize*G->BlockSize)+1+j+i*2*A.NFermion].Re+=(1-rest)*PsiChar.VRe[i][j];
                G->V[(iE+2)*(1+G->BlockSize*G->BlockSize)+1+j+i*2*A.NFermion].Re+=(  rest)*PsiChar.VRe[i][j];
                if(T.QComplex)
                {
                  G->V[(iE+1)*(1+G->BlockSize*G->BlockSize)+1+j+i*2*A.NFermion].Im+=(1-rest)*PsiChar.VIm[i][j];
                  G->V[(iE+2)*(1+G->BlockSize*G->BlockSize)+1+j+i*2*A.NFermion].Im+=(  rest)*PsiChar.VIm[i][j];
                }
              }
            }
          }
        }
      }
      FreeCompactMatrix(&T);
    }
  }
  FreeMatrix(&PsiChar);
  Nktot=Nk[0]*Nk[1]*Nk[2];
  for(i=0;i<G->N*(1+G->BlockSize*G->BlockSize);i+=(1+G->BlockSize*G->BlockSize))
  {
    for(j=1;j<1+G->BlockSize*G->BlockSize;j++)
    {
      G->V[i+j].Re/=Nktot;
      G->V[i+j].Im/=Nktot;
    }
  }
  ChopBlockListOfPoles(G, DBL_EPSILON);
  return 0;
}


/*  Given an tight binding Hamiltonian and a k-mesh defined by the number of divisions of a, b and c (Nk)
 *  returns the green function matrix of size 2 * Nvec * A.NFermion by 2 * Nvec * A.NFermion
 *          on an energy mesh with N+3 points. N+1 points from Emin to Emax and one point below Emin and one above Emax such that all weight and first moment is ok
 *  Calculates the eigen values and functions for each k including a self energy input as a Block List of poles
 *  for each k and eigen value i
 *  weight is added to the pole at iE and iE+1 such that Energy[iE] < val < Energy[iE+1]
 *  The matrix h_{jk} is fun[i,j] * fun[i,k]^* is added
 *
 *  Sigma.Blocksize must be equal to 2 * A.NFermion
 */
uint32_t TightBindingHamiltonianAddSpinPlusSelfEnergyMultiSiteToBlockListOfPoles(TightBindingDef const A, int Nk[3], uint32_t const Nvec, int const *vec, BlockListOfPoles const Sigma, BlockListOfPoles * const G, double const Emin, double const Emax, uint32_t const NE)
{
  MatrixType H, PsiChar;
  CompactMatrixType T, HSigma, HSigmaComplex;
  double q[3], val[A.NFermion*(Sigma.N+1)], rest, totmin=0, totmax=0, totnow, dvec[3];
  double const piornot = ((A.Cell.Units&1) ? (1) : (Pi));  // if A.Cell.Units last bit is 1 then r.g = pi, otherwise it evaluates to 1, piornot * r.g now eveluates always to pi.
  int    iq[3],ival,iE,i,j,ivec,jvec;
  double diq[3],dE=(Emax-Emin)/((double)NE),energy,Nktot;
  Complex16 phase;
  // initiallize PsiChar, a matrix storing the density of a wave-function
  PsiChar.Ni=2*A.NFermion;
  PsiChar.Nj=2*A.NFermion;
  PsiChar.QComplex=1;
  InitMatrix(&PsiChar);
  // test if sigma size is ok
  if(Sigma.BlockSize!=2*A.NFermion)
  {
    printf("in TightBindingHamiltonianAddSpinPlusSelfEnergyMultiSiteToBlockListOfPoles\nThe dimension of the self energy (%i) is not equal to the dimension of the tight binding Hamiltonian (%i) times 2 for the spin\n",Sigma.BlockSize,A.NFermion);
    return 1;
  }
  // init G
  G->BlockSize=2 * A.NFermion * Nvec;
  G->N=NE+3;
  G->QComplex=1;
  strcpy(G->Name, A.Name);
  if(InitBlockListOfPoles(G))
  {
    printf("InitBlockListOfPoles failed in TightBindingHamiltonianAddSpinPlusSelfEnergyMultiSiteToBlockListOfPoles\n");
    return 1;
  }
  // set energy of G
  energy=Emin-dE;
  for(i=0;i<G->N*(1+G->BlockSize*G->BlockSize);i+=(1+G->BlockSize*G->BlockSize))
  {
    G->V[i].Re=energy;
    energy+=dE;
  }
  // create HSigma
  BlockListOfPolesInversedToCompactMatrix(&Sigma, &HSigma);
  if(Sigma.QComplex==0)
  {
    CompactMatrixCopy(HSigma, &HSigmaComplex);
    MakeCompactRealMatrixComplex(&HSigmaComplex);
  }
  for (iq[0] = 0; iq[0] < Nk[0]; iq[0]++)
  {
    for (iq[1] = 0; iq[1] < Nk[1]; iq[1]++)
    {
      for (iq[2] = 0; iq[2] < Nk[2]; iq[2]++)
      {
        diq[0] = iq[0] / ((double)(Nk[0]));
        diq[1] = iq[1] / ((double)(Nk[1]));
        diq[2] = iq[2] / ((double)(Nk[2]));
        q[0] = diq[0] * A.Cell.gUnitCell[0][0] + diq[1] * A.Cell.gUnitCell[1][0] + diq[2] * A.Cell.gUnitCell[2][0];
        q[1] = diq[0] * A.Cell.gUnitCell[0][1] + diq[1] * A.Cell.gUnitCell[1][1] + diq[2] * A.Cell.gUnitCell[2][1];
        q[2] = diq[0] * A.Cell.gUnitCell[0][2] + diq[1] * A.Cell.gUnitCell[1][2] + diq[2] * A.Cell.gUnitCell[2][2];
        CreateHamiltonianFromTightBindingDefConstantPhaseInUnitCell(A, q, &H);
        if(H.QComplex)
        {
          if(Sigma.QComplex)
          {
            for(i=0;i<A.NFermion;i++)
            {
              for(j=0;j<A.NFermion;j++)
              {
                HSigma.V[i*HSigma.Nj+j].Re=Sigma.a0[i*Sigma.BlockSize+j].Re + H.VRe[i][j];
                HSigma.V[i*HSigma.Nj+j].Im=Sigma.a0[i*Sigma.BlockSize+j].Im + H.VIm[i][j];
                HSigma.V[(i+A.NFermion)*HSigma.Nj+(j+A.NFermion)].Re=Sigma.a0[(i+A.NFermion)*Sigma.BlockSize+(j+A.NFermion)].Re + H.VRe[i][j];
                HSigma.V[(i+A.NFermion)*HSigma.Nj+(j+A.NFermion)].Im=Sigma.a0[(i+A.NFermion)*Sigma.BlockSize+(j+A.NFermion)].Im + H.VIm[i][j];
              }
            }
            CompactMatrixEigensystem(&HSigma, &T, val);
          }
          else
          {
            for(i=0;i<A.NFermion;i++)
            {
              for(j=0;j<A.NFermion;j++)
              {
                HSigmaComplex.V[i*HSigma.Nj+j].Re=Sigma.a0Re[i*Sigma.BlockSize+j] + H.VRe[i][j];
                HSigmaComplex.V[i*HSigma.Nj+j].Im=                                  H.VIm[i][j];
                HSigmaComplex.V[(i+A.NFermion)*HSigma.Nj+(j+A.NFermion)].Re=Sigma.a0Re[(i+A.NFermion)*Sigma.BlockSize+(j+A.NFermion)] + H.VRe[i][j];
                HSigmaComplex.V[(i+A.NFermion)*HSigma.Nj+(j+A.NFermion)].Im=                                                            H.VIm[i][j];
              }
            }
            CompactMatrixEigensystem(&HSigmaComplex, &T, val);
          }
        }
        else
        {
          if(Sigma.QComplex)
          {
            for(i=0;i<A.NFermion;i++)
            {
              for(j=0;j<A.NFermion;j++)
              {
                HSigma.V[i*HSigma.Nj+j].Re=Sigma.a0[i*Sigma.BlockSize+j].Re + H.VRe[i][j];
                HSigma.V[i*HSigma.Nj+j].Im=Sigma.a0[i*Sigma.BlockSize+j].Im;
                HSigma.V[(i+A.NFermion)*HSigma.Nj+(j+A.NFermion)].Re=Sigma.a0[(i+A.NFermion)*Sigma.BlockSize+(j+A.NFermion)].Re + H.VRe[i][j];
                HSigma.V[(i+A.NFermion)*HSigma.Nj+(j+A.NFermion)].Im=Sigma.a0[(i+A.NFermion)*Sigma.BlockSize+(j+A.NFermion)].Im;
              }
            }
            CompactMatrixEigensystem(&HSigma, &T, val);
          }
          else
          {
            for(i=0;i<A.NFermion;i++)
            {
              for(j=0;j<A.NFermion;j++)
              {
                HSigma.VRe[i*HSigma.Nj+j]=Sigma.a0Re[i*Sigma.BlockSize+j] + H.VRe[i][j];
                HSigma.VRe[(i+A.NFermion)*HSigma.Nj+(j+A.NFermion)]=Sigma.a0Re[(i+A.NFermion)*Sigma.BlockSize+(j+A.NFermion)] + H.VRe[i][j];
              }
            }
            CompactMatrixEigensystem(&HSigma, &T, val);
          }
        }
        FreeMatrix(&H);
        for(ival=0;ival<T.Nj;ival++)
        {
          /* calculate character matrix for eigen-value ival - note that H changes meaning from Hamiltonian to Character matrix */
          if(T.QComplex)
          {
            ComplexEigenFunctionToCharacter(&(T.V[ival*T.Nj]), &PsiChar);
          }
          else
          {
            RealEigenFunctionToCharacter(&(T.VRe[ival*T.Nj]), &PsiChar);
          }
          /* determine the position of the pole in the list of poles */
          rest = (val[ival] - Emin) / dE;
          iE = (int)(floor(rest));
          if (iE >= NE)
          {
            totnow=0;
            for(i=0;i<2*A.NFermion;i++)
            {
              totnow+=PsiChar.VRe[i][i];
            }
            if(totmax+totnow<10*DBL_MIN)
            {
              totmax+=10*DBL_MIN;
              totnow+=10*DBL_MIN;
            }
            G->V[(NE+2)*(1+G->BlockSize*G->BlockSize)].Re = (totmax * G->V[(NE+2)*(1+G->BlockSize*G->BlockSize)].Re + totnow * val[ival]) / (totmax + totnow);
            for(i=0;i<2*A.NFermion;i++)
            {
              for(ivec=0;ivec<Nvec;ivec++)
              {
                for(j=0;j<2*A.NFermion;j++)
                {
                  for(jvec=0;jvec<Nvec;jvec++)
                  {
                    dvec[0] = (vec[3*ivec+0]-vec[3*jvec+0]) * A.Cell.rUnitCell[0][0] + (vec[3*ivec+1]-vec[3*jvec+1]) * A.Cell.rUnitCell[1][0] + (vec[3*ivec+2]-vec[3*jvec+2]) * A.Cell.rUnitCell[2][0];
                    dvec[1] = (vec[3*ivec+0]-vec[3*jvec+0]) * A.Cell.rUnitCell[0][1] + (vec[3*ivec+1]-vec[3*jvec+1]) * A.Cell.rUnitCell[1][1] + (vec[3*ivec+2]-vec[3*jvec+2]) * A.Cell.rUnitCell[2][1];
                    dvec[2] = (vec[3*ivec+0]-vec[3*jvec+0]) * A.Cell.rUnitCell[0][2] + (vec[3*ivec+1]-vec[3*jvec+1]) * A.Cell.rUnitCell[1][2] + (vec[3*ivec+2]-vec[3*jvec+2]) * A.Cell.rUnitCell[2][2];
                    phase.Re=cos(piornot * Dot(dvec,q));
                    phase.Im=sin(piornot * Dot(dvec,q));
                    G->V[(NE+2)*(1+G->BlockSize*G->BlockSize)+1 + j+jvec*2*A.NFermion + (i+ivec*2*A.NFermion)*G->BlockSize].Re += phase.Re*PsiChar.VRe[i][j];
                    G->V[(NE+2)*(1+G->BlockSize*G->BlockSize)+1 + j+jvec*2*A.NFermion + (i+ivec*2*A.NFermion)*G->BlockSize].Im += phase.Im*PsiChar.VRe[i][j];
                    if(T.QComplex)
                    {
                      G->V[(NE+2)*(1+G->BlockSize*G->BlockSize)+1 + j+jvec*2*A.NFermion + (i+ivec*2*A.NFermion)*G->BlockSize].Re -= phase.Im*PsiChar.VIm[i][j];
                      G->V[(NE+2)*(1+G->BlockSize*G->BlockSize)+1 + j+jvec*2*A.NFermion + (i+ivec*2*A.NFermion)*G->BlockSize].Im += phase.Re*PsiChar.VIm[i][j];
                    }
                  }
                }
              }
            }
            totmax+=totnow;
          }
          else if (iE < 0)
          {
            totnow=0;
            for(i=0;i<2*A.NFermion;i++)
            {
              totnow+=PsiChar.VRe[i][i];
            }
            if(totmin+totnow<10*DBL_MIN)
            {
              totmin+=10*DBL_MIN;
              totnow+=10*DBL_MIN;
            }
            G->V[0].Re = (totmin * G->V[0].Re + totnow * val[ival]) / (totmin + totnow);
            for(i=0;i<2*A.NFermion;i++)
            {
              for(ivec=0;ivec<Nvec;ivec++)
              {
                for(j=0;j<2*A.NFermion;j++)
                {
                  for(jvec=0;jvec<Nvec;jvec++)
                  {
                    dvec[0] = (vec[3*ivec+0]-vec[3*jvec+0]) * A.Cell.rUnitCell[0][0] + (vec[3*ivec+1]-vec[3*jvec+1]) * A.Cell.rUnitCell[1][0] + (vec[3*ivec+2]-vec[3*jvec+2]) * A.Cell.rUnitCell[2][0];
                    dvec[1] = (vec[3*ivec+0]-vec[3*jvec+0]) * A.Cell.rUnitCell[0][1] + (vec[3*ivec+1]-vec[3*jvec+1]) * A.Cell.rUnitCell[1][1] + (vec[3*ivec+2]-vec[3*jvec+2]) * A.Cell.rUnitCell[2][1];
                    dvec[2] = (vec[3*ivec+0]-vec[3*jvec+0]) * A.Cell.rUnitCell[0][2] + (vec[3*ivec+1]-vec[3*jvec+1]) * A.Cell.rUnitCell[1][2] + (vec[3*ivec+2]-vec[3*jvec+2]) * A.Cell.rUnitCell[2][2];
                    phase.Re=cos(piornot * Dot(dvec,q));
                    phase.Im=sin(piornot * Dot(dvec,q));
                    G->V[1 + j+jvec*2*A.NFermion + (i+ivec*2*A.NFermion)*G->BlockSize].Re += phase.Re*PsiChar.VRe[i][j];
                    G->V[1 + j+jvec*2*A.NFermion + (i+ivec*2*A.NFermion)*G->BlockSize].Im += phase.Im*PsiChar.VRe[i][j];
                    if(T.QComplex)
                    {
                      G->V[1 + j+jvec*2*A.NFermion + (i+ivec*2*A.NFermion)*G->BlockSize].Re -= phase.Im*PsiChar.VIm[i][j];
                      G->V[1 + j+jvec*2*A.NFermion + (i+ivec*2*A.NFermion)*G->BlockSize].Im += phase.Re*PsiChar.VIm[i][j];
                    }
                  }
                }
              }
            }
            totmin+=totnow;
          }
          else
          {
            rest -= iE;  /* rest is a number between 0 and 1 pole at iE+1 gets 1-rest pole at iE+2 gest rest */
            for(ivec=0;ivec<Nvec;ivec++)
            {
              for(jvec=0;jvec<Nvec;jvec++)
              {
                dvec[0] = (vec[3*ivec+0]-vec[3*jvec+0]) * A.Cell.rUnitCell[0][0] + (vec[3*ivec+1]-vec[3*jvec+1]) * A.Cell.rUnitCell[1][0] + (vec[3*ivec+2]-vec[3*jvec+2]) * A.Cell.rUnitCell[2][0];
                dvec[1] = (vec[3*ivec+0]-vec[3*jvec+0]) * A.Cell.rUnitCell[0][1] + (vec[3*ivec+1]-vec[3*jvec+1]) * A.Cell.rUnitCell[1][1] + (vec[3*ivec+2]-vec[3*jvec+2]) * A.Cell.rUnitCell[2][1];
                dvec[2] = (vec[3*ivec+0]-vec[3*jvec+0]) * A.Cell.rUnitCell[0][2] + (vec[3*ivec+1]-vec[3*jvec+1]) * A.Cell.rUnitCell[1][2] + (vec[3*ivec+2]-vec[3*jvec+2]) * A.Cell.rUnitCell[2][2];
                phase.Re=cos(piornot * Dot(dvec,q));
                phase.Im=sin(piornot * Dot(dvec,q));
                for(i=0;i<2*A.NFermion;i++)
                {
                  for(j=0;j<2*A.NFermion;j++)
                  {
                    G->V[(iE+1)*(1+G->BlockSize*G->BlockSize)+1 + j+jvec*2*A.NFermion + (i+ivec*2*A.NFermion)*G->BlockSize].Re += phase.Re*(1-rest)*PsiChar.VRe[i][j];
                    G->V[(iE+1)*(1+G->BlockSize*G->BlockSize)+1 + j+jvec*2*A.NFermion + (i+ivec*2*A.NFermion)*G->BlockSize].Im += phase.Im*(1-rest)*PsiChar.VRe[i][j];
                    G->V[(iE+2)*(1+G->BlockSize*G->BlockSize)+1 + j+jvec*2*A.NFermion + (i+ivec*2*A.NFermion)*G->BlockSize].Re += phase.Re*(  rest)*PsiChar.VRe[i][j];
                    G->V[(iE+2)*(1+G->BlockSize*G->BlockSize)+1 + j+jvec*2*A.NFermion + (i+ivec*2*A.NFermion)*G->BlockSize].Im += phase.Im*(  rest)*PsiChar.VRe[i][j];
                    if(T.QComplex)
                    {
                      G->V[(iE+1)*(1+G->BlockSize*G->BlockSize)+1 + j+jvec*2*A.NFermion + (i+ivec*2*A.NFermion)*G->BlockSize].Re -= phase.Im*(1-rest)*PsiChar.VIm[i][j];
                      G->V[(iE+1)*(1+G->BlockSize*G->BlockSize)+1 + j+jvec*2*A.NFermion + (i+ivec*2*A.NFermion)*G->BlockSize].Im += phase.Re*(1-rest)*PsiChar.VIm[i][j];
                      G->V[(iE+2)*(1+G->BlockSize*G->BlockSize)+1 + j+jvec*2*A.NFermion + (i+ivec*2*A.NFermion)*G->BlockSize].Re -= phase.Im*(  rest)*PsiChar.VIm[i][j];
                      G->V[(iE+2)*(1+G->BlockSize*G->BlockSize)+1 + j+jvec*2*A.NFermion + (i+ivec*2*A.NFermion)*G->BlockSize].Im += phase.Re*(  rest)*PsiChar.VIm[i][j];
                    }
                  }
                }
              }
            }
          }
        }
      }
      FreeCompactMatrix(&T);
    }
  }
  FreeMatrix(&PsiChar);
  Nktot=Nk[0]*Nk[1]*Nk[2];
  for(i=0;i<G->N*(1+G->BlockSize*G->BlockSize);i+=(1+G->BlockSize*G->BlockSize))
  {
    for(j=1;j<1+G->BlockSize*G->BlockSize;j++)
    {
      G->V[i+j].Re/=Nktot;
      G->V[i+j].Im/=Nktot;
    }
  }
  ChopBlockListOfPoles(G, DBL_EPSILON);
  return 0;
}





/*  On input
 *    matrix A is diagonal and real and of dimension psi->Nj by psi->Nj. Values are stored as a vector of doubles
 *    matrix psi is of dimension psi->Ni * psi->Nj
 *
 *  Calculated is a block tridiagonal matrix ABlock that has the same Green's function as psi . 1/(w-A) . psi^{C}
 *  aka [1/(w-ABlock)]_[[1;;psi->Ni,1;;psi->Ni]] = psi . 1/(w-A) . psi^{C}
 *
 *  on output matrix T is of dimension psi->Ni*NBlock by psi->Nj such that T^{dag} A T is block tridiagonal of dimension Ni*NBlock with psi on element 1 to psi->Ni
 *    ABlock = T^dag A T
 *
 * if psi is non orthogonal then an orthogonal form of psi is stored as the first block in T. b0 is psi psi_{ortho}^+
 */
uint32_t DiagonalMatrixBlockTridiagonalize(double const *A, CompactMatrixType const *psi, CompactMatrixType *const T, BlockTriDiagonalMatrixType *const ABlock, double NBlock, double SingularValue)
{
  char const *charC="C";
  char const *charN="N";
  int i,ii,j,k,oldindex,newindex,ipos,iipos,ia,ib,jcol,icol,allzero;
  double norm,dot,dotRe,dotIm;
  CompactMatrixType Tnow;

  
  printf("In DiagonalMatrixBlockTridiagonalize\n");
  PrintCompactMatrixMinimal(*psi);
  for(i=0;i<psi->Nj;i++)
  {
    printf("A[%2i] = %22.15E\n",i,A[i]);
  }
  printf("\n\n\n %g \n",NBlock);
  
  
  
  ABlock->N=NBlock;
  ABlock->BlockSize=psi->Ni;
  ABlock->QComplex=psi->QComplex;
  InitBlockTridiagonalMatrix(ABlock);
  T->Ni=psi->Ni*NBlock;
  T->Nj=psi->Nj;
  T->QComplex=psi->QComplex;
  if(InitCompactMatrix(T))
  {
    printf("InitCompactMatrix failed in DiagonalMatrixBlockTridiagonalize\n");
    return 1;
  }
  Tnow.Name=T->Name;
  Tnow.QComplex=T->QComplex;
  if(T->QComplex)
  {/* calculate T */
    Complex16 one, zero;
    one.Re=1.0;  one.Im=0.0;
    zero.Re=0.0; zero.Im=0.0;
    /* set first vector block */
    memcpy(T->V, psi->V, psi->Nj*psi->Ni*sizeof(Complex16));
    Tnow.V=T->V;
    Tnow.Ni=psi->Ni;
    Tnow.Nj=psi->Nj;
    /* orthonormalize Tnow */
    CompactMatrixGramSchmidtOrthogonalizeRow(&Tnow, SingularValue, true);
    zgemm(charC,charN,&ABlock->BlockSize,&ABlock->BlockSize,&psi->Nj,&one,Tnow.V,&psi->Nj,psi->V,&psi->Nj,&zero,ABlock->b,&ABlock->BlockSize);
    /* loop over psi_{n+1} = A psi_{n} and orthonormalization */
    oldindex=0;
    newindex=psi->Ni*psi->Nj;
    for(k=1;k<ABlock->N;k++)
    {
      Tnow.Ni+=psi->Ni;
      for(i=0;i<psi->Ni;i++)
      {
        for(j=0;j<psi->Nj;j++)
        {
          Tnow.V[newindex].Re = A[j] * Tnow.V[oldindex].Re;
          Tnow.V[newindex].Im = A[j] * Tnow.V[oldindex].Im;
          oldindex++;
          newindex++;
        }
      }
      allzero=1;
      for(i=Tnow.Ni-psi->Ni;i<Tnow.Ni;i++)
      {
        ipos=i*Tnow.Nj;
        /* orthogonalize i to ii=0 to i */
        for(ii=0;ii<i;ii++)
        {
          iipos=ii*Tnow.Nj;
          dotRe=0;
          dotIm=0;
          for(j=ipos;j<ipos+Tnow.Nj;j++)
          {
            dotRe+=Tnow.V[iipos].Re*Tnow.V[j].Re + Tnow.V[iipos].Im*Tnow.V[j].Im;
            dotIm+=Tnow.V[iipos].Re*Tnow.V[j].Im - Tnow.V[iipos].Im*Tnow.V[j].Re;
            iipos++;
          }
          iipos=ii*Tnow.Nj;
          for(j=ipos;j<ipos+Tnow.Nj;j++)
          {
            Tnow.V[j].Re -=dotRe*Tnow.V[iipos].Re - dotIm*Tnow.V[iipos].Im;
            Tnow.V[j].Im -=dotRe*Tnow.V[iipos].Im + dotIm*Tnow.V[iipos].Re;
            iipos++;
          }
        }
        /* normalize vector i */
        norm=0;
        for(j=ipos;j<ipos+Tnow.Nj;j++)
        {
          norm+=Sqr(Tnow.V[j].Re)+Sqr(Tnow.V[j].Im);
        }
        if(norm<SingularValue)
        {
          norm=0;
        }
        else
        {
          allzero=0;
          norm=1.0/sqrt(norm);
        }
        for(j=ipos;j<ipos+Tnow.Nj;j++)
        {
          Tnow.V[j].Re*=norm;
          Tnow.V[j].Im*=norm;
        }
      }
      if(allzero)
      {
        ABlock->N=k;
      }
    }
    
    
    printf("transformation matrix T\n");
    PrintCompactMatrixMinimal(*T);
    CompactMatrixType test1, test2;
    CompactMatrixMultiply('N', 'C', T, T, &test1);
    CompactMatrixChop(&test1, 1000*DBL_EPSILON);
    printf("matrix T T^C = test1 should be identity matrix\n");
    PrintCompactMatrixMinimal(test1);
    TransformDiagonalMatrixByCompactMatrixTconATtran(A, T, &test2);
    CompactMatrixChop(&test2, 1000*DBL_EPSILON);
    printf("matrix T^C A T = test2\n");
    PrintCompactMatrixMinimal(test2);
    

    /*  calculate ABlock = T^dag A T
     *  ABlock(i,j) = sum_k T(k,i)^* A(k,k) T(k,j)
     */
    iipos=ABlock->BlockSize*ABlock->BlockSize; /* block position in tridiagonal matrix */
    for(ia=0;ia<ABlock->N;ia++)
    {
      for(i=0;i<ABlock->BlockSize;i++)
      {
        for(j=0;j<ABlock->BlockSize;j++)
        {
          icol=(i+ia*ABlock->BlockSize)*T->Nj;
          jcol=(j+ia*ABlock->BlockSize)*T->Nj;
          ABlock->a[iipos].Re=0;
          ABlock->a[iipos].Im=0;
          for(k=0;k<T->Nj;k++)
          {
            ABlock->a[iipos].Re+=T->V[jcol].Re*A[k]*T->V[icol].Re + T->V[jcol].Im*A[k]*T->V[icol].Im;
            ABlock->a[iipos].Im+=T->V[jcol].Re*A[k]*T->V[icol].Im - T->V[jcol].Im*A[k]*T->V[icol].Re;
            icol++;
            jcol++;
          }
          iipos++;
        }
      }
    }
    iipos=ABlock->BlockSize*ABlock->BlockSize; /* position in tridiagonal matrix */
    for(ib=0;ib<(ABlock->N-1);ib++)
    {
      for(i=0;i<ABlock->BlockSize;i++)
      {
        for(j=0;j<ABlock->BlockSize;j++)
        {
          icol=(i+ib*ABlock->BlockSize)*T->Nj;
          jcol=(j+(ib+1)*ABlock->BlockSize)*T->Nj;
          ABlock->b[iipos].Re=0;
          ABlock->b[iipos].Im=0;
          for(k=0;k<T->Nj;k++)
          {
            ABlock->b[iipos].Re+=T->V[jcol].Re*A[k]*T->V[icol].Re + T->V[jcol].Im*A[k]*T->V[icol].Im;
            ABlock->b[iipos].Im+=T->V[jcol].Re*A[k]*T->V[icol].Im - T->V[jcol].Im*A[k]*T->V[icol].Re;
            icol++;
            jcol++;
          }
          iipos++;
        }
      }
    }
  }
  else /* real */
  {/* calculate T */
    double one=1.0, zero=0.0;
    /* set first vector block */
    memcpy(T->VRe, psi->VRe, psi->Nj*psi->Ni*sizeof(double));
    Tnow.VRe=T->VRe;
    Tnow.Ni=psi->Ni;
    Tnow.Nj=psi->Nj;
    /* orthonormalize Tnow */
    CompactMatrixGramSchmidtOrthogonalizeRow(&Tnow, SingularValue, true);
    dgemm(charC,charN,&ABlock->BlockSize,&ABlock->BlockSize,&psi->Nj,&one,Tnow.VRe,&psi->Nj,psi->VRe,&psi->Nj,&zero,ABlock->bRe,&ABlock->BlockSize);
    /* loop over psi_{n+1} = A psi_{n} and orthonormalization */
    oldindex=0;
    newindex=psi->Ni*psi->Nj;
    for(k=1;k<ABlock->N;k++)
    {
      Tnow.Ni+=psi->Ni;
      for(i=0;i<psi->Ni;i++)
      {
        for(j=0;j<psi->Nj;j++)
        {
          Tnow.VRe[newindex] = A[j] * Tnow.VRe[oldindex];
          oldindex++;
          newindex++;
        }
      }
      allzero=1;
      for(i=Tnow.Ni-psi->Ni;i<Tnow.Ni;i++)
      {
        ipos=i*Tnow.Nj;
        /* orthogonalize i to j=0 to i */
        for(ii=0;ii<i;ii++)
        {
          iipos=ii*Tnow.Nj;
          dot=0;
          for(j=ipos;j<ipos+Tnow.Nj;j++)
          {
            dot+=Tnow.VRe[j]*Tnow.VRe[iipos];
            iipos++;
          }
          iipos=ii*Tnow.Nj;
          for(j=ipos;j<ipos+Tnow.Nj;j++)
          {
            Tnow.VRe[j] -=dot*Tnow.VRe[iipos];
            iipos++;
          }
        }
        /* normalize vector i */
        norm=0;
        for(j=ipos;j<ipos+Tnow.Nj;j++)
        {
          norm+=Sqr(Tnow.VRe[j]);
        }
        if(norm<SingularValue)
        {
          norm=0;
        }
        else
        {
          allzero=0;
          norm=1.0/sqrt(norm);
        }
        for(j=ipos;j<ipos+Tnow.Nj;j++)
        {
          Tnow.VRe[j]*=norm;
        }
      }
      if(allzero)
      {
        ABlock->N=k;
      }
    }
    /* calculate T^+ A T */
    iipos=ABlock->BlockSize*ABlock->BlockSize; /* position in tridiagonal matrix */
    for(ia=0;ia<ABlock->N;ia++)
    {
      for(i=0;i<ABlock->BlockSize;i++)
      {
        for(j=0;j<ABlock->BlockSize;j++)
        {
          icol=(i+ia*ABlock->BlockSize)*T->Nj;
          jcol=(j+ia*ABlock->BlockSize)*T->Nj;
          ABlock->aRe[iipos]=0;
          for(k=0;k<T->Nj;k++)
          {
            ABlock->aRe[iipos]+=T->VRe[jcol]*A[k]*T->VRe[icol];
            icol++;
            jcol++;
          }
          iipos++;
        }
      }
    }
    iipos=ABlock->BlockSize*ABlock->BlockSize; /* position in tridiagonal matrix */
    for(ib=0;ib<(ABlock->N-1);ib++)
    {
      for(i=0;i<ABlock->BlockSize;i++)
      {
        for(j=0;j<ABlock->BlockSize;j++)
        {
          icol=(i+ib*ABlock->BlockSize)*T->Nj;
          jcol=(j+(ib+1)*ABlock->BlockSize)*T->Nj;
          ABlock->bRe[iipos]=0;
          for(k=0;k<T->Nj;k++)
          {
            ABlock->bRe[iipos]+=T->VRe[jcol]*A[k]*T->VRe[icol];
            icol++;
            jcol++;
          }
          iipos++;
        }
      }
    }
  }
  return 0;
}



/*  =========================================================================
 *  ... to TigthBinding
 *  ========================================================================= */

/* transform a Block matrix into a tight binding Hamiltonian, such that the one particle Green's function 
 * is equal to the response function represented by the Block matrix
 * This is only possible when A0=0 and B0=1
 *
 * The tight binding contains N Atoms labeled from 1 to N
 * Each atom has M->BlockSize orbitals
 */

uint32_t BlockMatrixToTightBinding(BlockListOfPoles const * const M, TightBindingDef * TB)
{
  //TODO: BlockListOfPoles to TB
  return 0;
}
/* The onsite energy of each atom is given by A_i
 * The hopping between atom i and i+1 is given by B_i
 */
uint32_t BlockMatrixToTightBinding(BlockTriDiagonalMatrixType const * const M, TightBindingDef * TB)
{
  uint32_t i,iH,j,iorb,jorb;
  strcpy(TB->Name,M->Name);
  TB->Cell.rUnitCell[0][0]=INFINITY; TB->Cell.rUnitCell[0][1]=0;        TB->Cell.rUnitCell[0][2]=0;
  TB->Cell.rUnitCell[1][0]=0;        TB->Cell.rUnitCell[1][1]=INFINITY; TB->Cell.rUnitCell[1][2]=0;
  TB->Cell.rUnitCell[2][0]=0;        TB->Cell.rUnitCell[2][1]=0;        TB->Cell.rUnitCell[2][2]=INFINITY;
  TB->Cell.gUnitCell[0][0]=0;        TB->Cell.gUnitCell[0][1]=0;        TB->Cell.gUnitCell[0][2]=0;
  TB->Cell.gUnitCell[1][0]=0;        TB->Cell.gUnitCell[1][1]=0;        TB->Cell.gUnitCell[1][2]=0;
  TB->Cell.gUnitCell[2][0]=0;        TB->Cell.gUnitCell[2][1]=0;        TB->Cell.gUnitCell[2][2]=0;
  TB->Cell.NAtoms=M->N;
  TB->Cell.Atom=(AtomType*)malloc(TB->Cell.NAtoms*sizeof(AtomType));
  j=1;
  if(M->N>9)
  {
    j=2;
    if(M->N>99)
    {
      j=3;
      if(M->N>999)
      {
        j=4;
        if(M->N>9999)
        {
          j=5;
          if(M->N>99999)
          {
            j=6;
            if(M->N>999999)
            {
              j=10;
            }
          }
        }
      }
    }
  }
  jorb=1;
  if(M->BlockSize>9)
  {
    jorb=2;
    if(M->BlockSize>99)
    {
      jorb=3;
      if(M->BlockSize>999)
      {
        jorb=4;
        if(M->BlockSize>9999)
        {
          jorb=5;
          if(M->BlockSize>99999)
          {
            jorb=6;
            if(M->BlockSize>999999)
            {
              jorb=10;
            }
          }
        }
      }
    }
  }
  for(i=0;i<TB->Cell.NAtoms;i++)
  {
    sprintf(TB->Cell.Atom[i].Name, "%0*d",j,i);
    TB->Cell.Atom[i].AtomNumber=0;
    TB->Cell.Atom[i].Position[0]=0;  TB->Cell.Atom[i].Position[1]=0;  TB->Cell.Atom[i].Position[2]=0;
    TB->Cell.Atom[i].Nshell=1;
    TB->Cell.Atom[i].ShellName=(char**)(malloc(sizeof(char*)));
    TB->Cell.Atom[i].ShellName[0]=(char*)(malloc(2*sizeof(char)));
    TB->Cell.Atom[i].ShellName[0][0]='x';
    TB->Cell.Atom[i].ShellName[0][1]='\0';
    TB->Cell.Atom[i].shelldegeneracy=(int*)(malloc(sizeof(int)));
    TB->Cell.Atom[i].shelldegeneracy[0]=M->BlockSize;
    TB->Cell.Atom[i].OrbitalName=(char***)(malloc(sizeof(char**)));
    TB->Cell.Atom[i].OrbitalName[0]=(char**)(malloc(TB->Cell.Atom[i].shelldegeneracy[0]*sizeof(char*)));
    for(iorb=0;iorb<M->BlockSize;iorb++)
    {
      TB->Cell.Atom[i].OrbitalName[0][iorb]=(char*)malloc((j+1)*sizeof(char));
      sprintf(TB->Cell.Atom[i].OrbitalName[0][iorb], "%0*d",jorb,iorb);
    }
  }
  TB->Cell.Units=0;
  
  if(M->N==0)
  {
    TB->nH=0;
  }
  else
  {
    TB->nH=2*M->N-1;
  }
  TB->H=(TightBindingDefHopping*)malloc(TB->nH*sizeof(TightBindingDefHopping));
  TB->NFermion=M->N*M->BlockSize;
  SetNFermionIndexTightBindingDef(TB);
  iH=0;
  // onsite energy
  for(i=1;i<=M->N;i++)
  {
    TB->H[iH].Atompos1=i-1;
    TB->H[iH].Shellpos1=0;
    TB->H[iH].Atompos2=i-1;
    TB->H[iH].Shellpos2=0;
    TB->H[iH].TransVecs[0]=0;
    TB->H[iH].TransVecs[1]=0;
    TB->H[iH].TransVecs[2]=0;
    TB->H[iH].ConnectingVector[0]=0;
    TB->H[iH].ConnectingVector[1]=0;
    TB->H[iH].ConnectingVector[2]=0;
    TB->H[iH].Hop.Ni=M->BlockSize;
    TB->H[iH].Hop.Nj=M->BlockSize;
    TB->H[iH].Hop.QComplex=M->QComplex;
    InitMatrix(&(TB->H[iH].Hop));
    if(M->QComplex==0)
    {
      for(j=0;j<Sqr(M->BlockSize);j++)
      {
        TB->H[iH].Hop.VRe[0][j]=M->aRe[(i)*Sqr(M->BlockSize)+j];
      }
    }
    else
    {
      for(j=0;j<Sqr(M->BlockSize);j++)
      {
        TB->H[iH].Hop.VRe[0][j]=M->a[(i)*Sqr(M->BlockSize)+j].Re;
        TB->H[iH].Hop.VIm[0][j]=M->a[(i)*Sqr(M->BlockSize)+j].Im;
      }
    }
    iH++;
  }
  // hopping
  for(i=1;i<M->N;i++)
  {
    TB->H[iH].Atompos1=i-1;
    TB->H[iH].Shellpos1=0;
    TB->H[iH].Atompos2=i;
    TB->H[iH].Shellpos2=0;
    TB->H[iH].TransVecs[0]=0;
    TB->H[iH].TransVecs[1]=0;
    TB->H[iH].TransVecs[2]=0;
    TB->H[iH].ConnectingVector[0]=0;
    TB->H[iH].ConnectingVector[1]=0;
    TB->H[iH].ConnectingVector[2]=0;
    TB->H[iH].Hop.Ni=M->BlockSize;
    TB->H[iH].Hop.Nj=M->BlockSize;
    TB->H[iH].Hop.QComplex=M->QComplex;
    InitMatrix(&(TB->H[iH].Hop));
    if(M->QComplex==0)
    {
      for(j=0;j<Sqr(M->BlockSize);j++)
      {
        TB->H[iH].Hop.VRe[0][j]=M->bRe[(i)*Sqr(M->BlockSize)+j];
      }
    }
    else
    {
      for(j=0;j<Sqr(M->BlockSize);j++)
      {
        TB->H[iH].Hop.VRe[0][j]=M->b[(i)*Sqr(M->BlockSize)+j].Re;
        TB->H[iH].Hop.VIm[0][j]=M->b[(i)*Sqr(M->BlockSize)+j].Im;
      }
    }
    iH++;
  }
  return 0;
}

uint32_t BlockMatrixToTightBinding(BlockAndersonMatrixType const * const M, TightBindingDef * TB)
{
  uint32_t i,iH,j,iorb,jorb;
  strcpy(TB->Name,M->Name);
  TB->Cell.rUnitCell[0][0]=INFINITY; TB->Cell.rUnitCell[0][1]=0;        TB->Cell.rUnitCell[0][2]=0;
  TB->Cell.rUnitCell[1][0]=0;        TB->Cell.rUnitCell[1][1]=INFINITY; TB->Cell.rUnitCell[1][2]=0;
  TB->Cell.rUnitCell[2][0]=0;        TB->Cell.rUnitCell[2][1]=0;        TB->Cell.rUnitCell[2][2]=INFINITY;
  TB->Cell.gUnitCell[0][0]=0;        TB->Cell.gUnitCell[0][1]=0;        TB->Cell.gUnitCell[0][2]=0;
  TB->Cell.gUnitCell[1][0]=0;        TB->Cell.gUnitCell[1][1]=0;        TB->Cell.gUnitCell[1][2]=0;
  TB->Cell.gUnitCell[2][0]=0;        TB->Cell.gUnitCell[2][1]=0;        TB->Cell.gUnitCell[2][2]=0;
  TB->Cell.NAtoms=M->N;
  TB->Cell.Atom=(AtomType*)malloc(TB->Cell.NAtoms*sizeof(AtomType));
  j=1;
  if(M->N>9)
  {
    j=2;
    if(M->N>99)
    {
      j=3;
      if(M->N>999)
      {
        j=4;
        if(M->N>9999)
        {
          j=5;
          if(M->N>99999)
          {
            j=6;
            if(M->N>999999)
            {
              j=10;
            }
          }
        }
      }
    }
  }
  jorb=1;
  if(M->BlockSize>9)
  {
    jorb=2;
    if(M->BlockSize>99)
    {
      jorb=3;
      if(M->BlockSize>999)
      {
        jorb=4;
        if(M->BlockSize>9999)
        {
          jorb=5;
          if(M->BlockSize>99999)
          {
            jorb=6;
            if(M->BlockSize>999999)
            {
              jorb=10;
            }
          }
        }
      }
    }
  }
  for(i=0;i<TB->Cell.NAtoms;i++)
  {
    sprintf(TB->Cell.Atom[i].Name, "%0*d",j,i);
    TB->Cell.Atom[i].AtomNumber=0;
    TB->Cell.Atom[i].Position[0]=0;  TB->Cell.Atom[i].Position[1]=0;  TB->Cell.Atom[i].Position[2]=0;
    TB->Cell.Atom[i].Nshell=1;
    TB->Cell.Atom[i].ShellName=(char**)(malloc(sizeof(char*)));
    TB->Cell.Atom[i].ShellName[0]=(char*)(malloc(2*sizeof(char)));
    TB->Cell.Atom[i].ShellName[0][0]='x';
    TB->Cell.Atom[i].ShellName[0][1]='\0';
    TB->Cell.Atom[i].shelldegeneracy=(int*)(malloc(sizeof(int)));
    TB->Cell.Atom[i].shelldegeneracy[0]=M->BlockSize;
    TB->Cell.Atom[i].OrbitalName=(char***)(malloc(sizeof(char**)));
    TB->Cell.Atom[i].OrbitalName[0]=(char**)(malloc(TB->Cell.Atom[i].shelldegeneracy[0]*sizeof(char*)));
    for(iorb=0;iorb<M->BlockSize;iorb++)
    {
      TB->Cell.Atom[i].OrbitalName[0][iorb]=(char*)malloc((j+1)*sizeof(char));
      sprintf(TB->Cell.Atom[i].OrbitalName[0][iorb], "%0*d",jorb,iorb);
    }
  }
  TB->Cell.Units=0;
  
  if(M->N==0)
  {
    TB->nH=0;
  }
  else
  {
    TB->nH=2*M->N-1;
  }
  TB->H=(TightBindingDefHopping*)malloc(TB->nH*sizeof(TightBindingDefHopping));
  TB->NFermion=M->N*M->BlockSize;
  SetNFermionIndexTightBindingDef(TB);
  iH=0;
  // onsite energy
  for(i=1;i<=M->N;i++)
  {
    TB->H[iH].Atompos1=i-1;
    TB->H[iH].Shellpos1=0;
    TB->H[iH].Atompos2=i-1;
    TB->H[iH].Shellpos2=0;
    TB->H[iH].TransVecs[0]=0;
    TB->H[iH].TransVecs[1]=0;
    TB->H[iH].TransVecs[2]=0;
    TB->H[iH].ConnectingVector[0]=0;
    TB->H[iH].ConnectingVector[1]=0;
    TB->H[iH].ConnectingVector[2]=0;
    TB->H[iH].Hop.Ni=M->BlockSize;
    TB->H[iH].Hop.Nj=M->BlockSize;
    TB->H[iH].Hop.QComplex=M->QComplex;
    InitMatrix(&(TB->H[iH].Hop));
    if(M->QComplex==0)
    {
      for(j=0;j<Sqr(M->BlockSize);j++)
      {
        TB->H[iH].Hop.VRe[0][j]=M->aRe[(i)*Sqr(M->BlockSize)+j];
      }
    }
    else
    {
      for(j=0;j<Sqr(M->BlockSize);j++)
      {
        TB->H[iH].Hop.VRe[0][j]=M->a[(i)*Sqr(M->BlockSize)+j].Re;
        TB->H[iH].Hop.VIm[0][j]=M->a[(i)*Sqr(M->BlockSize)+j].Im;
      }
    }
    iH++;
  }
  // hopping
  for(i=1;i<M->N;i++)
  {
    TB->H[iH].Atompos1=0;
    TB->H[iH].Shellpos1=0;
    TB->H[iH].Atompos2=i;
    TB->H[iH].Shellpos2=0;
    TB->H[iH].TransVecs[0]=0;
    TB->H[iH].TransVecs[1]=0;
    TB->H[iH].TransVecs[2]=0;
    TB->H[iH].ConnectingVector[0]=0;
    TB->H[iH].ConnectingVector[1]=0;
    TB->H[iH].ConnectingVector[2]=0;
    TB->H[iH].Hop.Ni=M->BlockSize;
    TB->H[iH].Hop.Nj=M->BlockSize;
    TB->H[iH].Hop.QComplex=M->QComplex;
    InitMatrix(&(TB->H[iH].Hop));
    if(M->QComplex==0)
    {
      for(j=0;j<Sqr(M->BlockSize);j++)
      {
        TB->H[iH].Hop.VRe[0][j]=M->bRe[(i)*Sqr(M->BlockSize)+j];
      }
    }
    else
    {
      for(j=0;j<Sqr(M->BlockSize);j++)
      {
        TB->H[iH].Hop.VRe[0][j]=M->b[(i)*Sqr(M->BlockSize)+j].Re;
        TB->H[iH].Hop.VIm[0][j]=M->b[(i)*Sqr(M->BlockSize)+j].Im;
      }
    }
    iH++;
  }
  return 0;
}

uint32_t BlockMatrixToTightBinding(BlockNaturalImpurityOrbitalMatrixType const * const M, TightBindingDef * TB)
{
  //TODO: BlockNaturalImpurityOrbitalMatrixType to TB
  return 0;
}
