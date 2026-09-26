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
      Lorentzian.Re = (omega-A->V[itot].Re) / (Sqr(Gamma/2.0) + Sqr(omega-A->V[itot].Re));
      Lorentzian.Im = (-Gamma/2.0         ) / (Sqr(Gamma/2.0) + Sqr(omega-A->V[itot].Re));
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
  if(A->N!=0)
  {
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
  FreeCompactMatrix(&Gwi);
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

uint32_t EvaluateBlockDoubleTriDiagonalMatrix(BlockDoubleTriDiagonalMatrixType const * A, CompactMatrixType * Gw, const double omega, const double Gamma)
{
  Complex16 z = {omega,Gamma/2};
  CompactMatrixType val,con, det;
  EvaluateBlockTriDiagonalMatrix(&(A->val), &val, omega, Gamma);
  EvaluateBlockTriDiagonalMatrix(&(A->con), &con, omega, Gamma);
  det.Ni=val.Ni;
  det.Nj=val.Nj;
  det.QComplex=1;
  InitCompactIdentityMatrix(&det);
  CompactMatrixScaleComplex(&det, z);
  CompactMatrixSub(&det, &val);
  CompactMatrixSub(&det, &con);
  CompactMatrixInverse(&det);
  CompactMatrixMultiply3('*', 'N', 'T', &(A->b0), &det, &(A->b0), Gw);
  CompactMatrixAdd(Gw,&(A->a0));
  FreeCompactMatrix(&val);
  FreeCompactMatrix(&con);
  FreeCompactMatrix(&det);
  return 0;
}
uint32_t EvaluateBlockDoubleTriDiagonalMatrix(void const * A, CompactMatrixType * Gw, const double omega, const double Gamma)
{
  BlockDoubleTriDiagonalMatrixType *ATypecast = (BlockDoubleTriDiagonalMatrixType*) A;
  return EvaluateBlockDoubleTriDiagonalMatrix(ATypecast, Gw, omega, Gamma);
}


/*
 *
 */
uint32_t EvaluateBlockNaturalImpurityOrbital(BlockNaturalImpurityOrbitalMatrixType const *A, CompactMatrixType * Gw, const double omega, const double Gamma)
{
  CompactMatrixType con,val;
  if(CompactMatrixCopy(A->a0, Gw))
  {
    printf("CompactMatrixCopy failed in EvaluateBlockNaturalImpurityOrbital\n");
    return 1;
  }
  EvaluateBlockTriDiagonalMatrix(&(A->con), &con, omega, Gamma);
  CompactMatrixAdd(Gw, &con);
  FreeCompactMatrix(&con);
  EvaluateBlockTriDiagonalMatrix(&(A->val), &val, omega, Gamma);
  CompactMatrixAdd(Gw, &val);
  FreeCompactMatrix(&val);
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
    case 'd':
    case 'D':
    {
      return EvaluateBlockDoubleTriDiagonalMatrix(A->data, Gw, omega, Gamma);
    }
    default:
    {
      printf("Error in EvaluateBlockResponsefunctionType: Type of Response function should be either T/t,A/a,L/l, N/n, or D/d and not %c\n",A->type);
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
uint32_t BlockListOfPolesToSpectra(double GroundStateEnergy, BlockListOfPoles const *A, SpectraType *G, double Emin, double Emax, uint32_t N)
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
  strcpy(G[0].Name, A->Name);
  for(i=0;i<A->BlockSize*A->BlockSize;i++)
  {
    SpectrumAllocate(G[0], &(G[i]));
  }
  if(A->QComplex)
  {
    Gtemp= (SpectraType*)malloc(A->BlockSize*A->BlockSize*sizeof(SpectraType));
    for(i=0;i<A->BlockSize*A->BlockSize;i++)
    {
      SpectrumAllocate(G[0], &(Gtemp[i]));
    }
    /* j loops over all blocks in A */
    for(j=0;j<(A->BlockSize*A->BlockSize+1)*A->N;j+=A->BlockSize*A->BlockSize+1)
    {
      /* the poles have energy A.V[j].Re */
      if(A->V[j].Re<=Emin)
      {
        iEmin=0;
        iEplus=0;
        wmin=1;
        wplus=0;
      }
      else if(A->V[j].Re>=Emax)
      {
        iEmin=N;
        iEplus=N;
        wmin=1;
        wplus=0;
      }
      else
      {
        iEmin=(uint32_t)((A->V[j].Re-Emin)/G[0].dE);
        iEplus=iEmin+1;
        wplus=((A->V[j].Re-Emin)/G[0].dE)-iEmin;
        wmin=1-wplus;
      }
      /* add the im part of the lorenzian to energy iEmin and iEplus */
      for(i=0;i<A->BlockSize*A->BlockSize;i++)
      {
        Gtemp[i].Im[iEmin]-=wmin*A->V[j+i+1].Re;
        Gtemp[i].Re[iEmin]+=wmin*A->V[j+i+1].Im;
        Gtemp[i].Im[iEplus]-=wplus*A->V[j+i+1].Re;
        Gtemp[i].Re[iEplus]+=wplus*A->V[j+i+1].Im;
      }
    }
  }
  else
  {
    Gtemp=G;
    for(j=0;j<(A->BlockSize*A->BlockSize+1)*A->N;j+=A->BlockSize*A->BlockSize+1)
    {
      /* the poles have energy A.V[j].Re */
      if(A->VRe[j]<=Emin)
      {
        iEmin=0;
        iEplus=0;
        wmin=1;
        wplus=0;
      }
      else if(A->VRe[j]>=Emax)
      {
        iEmin=N;
        iEplus=N;
        wmin=1;
        wplus=0;
      }
      else
      {
        iEmin=(uint32_t)((A->VRe[j]-Emin)/G[0].dE);
        iEplus=iEmin+1;
        wplus=((A->VRe[j]-Emin)/G[0].dE)-iEmin;
        wmin=1-wplus;
      }
      /* add the im part of the lorenzian to energy iEmin and iEplus */
      for(i=0;i<A->BlockSize*A->BlockSize;i++)
      {
        Gtemp[i].Im[iEmin]-=wmin*A->VRe[j+i+1];
        Gtemp[i].Im[iEplus]-=wplus*A->VRe[j+i+1];
      }
    }
  }
  /* the im part of the lorenzians is there, we now need to add the real part */
  /* note that in the complex case some lorenzians have a complex pre-factor */
  /* weight w_i at energy E_i contributes w_i / (w - E_i) to all energies w, except at w=E_i */
  for(i=0;i<A->BlockSize*A->BlockSize;i++)
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
  if(A->QComplex)
  {
    for(i=0;i<A->BlockSize*A->BlockSize;i++)
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
  for(i=0;i<A->BlockSize*A->BlockSize;i++)
  {
    for(j=0;j<=N;j++)
    {
      G[i].Im[j]*=Pi/G[i].dE;
      G[i].Re[j]*=Pi/G[i].dE;
    }
  }
  /* add a0 */
  if(A->QComplex)
  {
    for(i=0;i<A->BlockSize*A->BlockSize;i++)
    {
      for(j=0;j<=N;j++)
      {
        G[i].Re[j]+=A->a0[i].Re;
        G[i].Im[j]+=A->a0[i].Im;
      }
    }
  }
  else
  {
    for(i=0;i<A->BlockSize*A->BlockSize;i++)
    {
      for(j=0;j<=N;j++)
      {
        G[i].Re[j]+=A->a0Re[i];
      }
    }
  }
  return 0;
}
uint32_t BlockListOfPolesToSpectraWithBroadening(double GroundStateEnergy, BlockListOfPoles const *A, SpectraType *G, double Emin, double Emax, uint32_t N, double Gamma)
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
  strcpy(G[0].Name, A->Name);
  for(i=0;i<A->BlockSize*A->BlockSize;i++)
  {
    SpectrumAllocate(G[0], &(G[i]));
  }
  if(A->QComplex==0)
  {// real
#pragma omp parallel for private(omega,i,itot,iblock,Lorentzian)
    for(iE=0;iE<=N;iE++)
    {
      omega = ((double)iE/(double)N)*(Emax-Emin)+Emin;
      // add A0
      for(i=0;i<A->BlockSize*A->BlockSize;i++)
      {
        G[i].Re[iE]=A->VRe[A->N*(1+A->BlockSize*A->BlockSize)+i];
        G[i].Im[iE]=0;
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
          G[i].Re[iE] += A->VRe[itot] * Lorentzian.Re;
          G[i].Im[iE] += A->VRe[itot] * Lorentzian.Im;
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
      for(i=0;i<A->BlockSize*A->BlockSize;i++)
      {
        G[i].Re[iE]=A->V[A->N*(1+A->BlockSize*A->BlockSize)+i].Re;
        G[i].Im[iE]=A->V[A->N*(1+A->BlockSize*A->BlockSize)+i].Im;
      }
      // loop over all poles
      itot=0;
      for(iblock=0;iblock<A->N;iblock++)
      {
        Lorentzian.Re = (omega-A->V[itot].Re) / (Sqr(Gamma/2.0) + Sqr(omega-A->V[itot].Re));
        Lorentzian.Im = (-Gamma/2.          ) / (Sqr(Gamma/2.0) + Sqr(omega-A->V[itot].Re));
        itot++;
        for(i=0;i<A->BlockSize*A->BlockSize;i++)
        {
          G[i].Re[iE] += A->V[itot].Re * Lorentzian.Re - A->V[itot].Im * Lorentzian.Im;
          G[i].Im[iE] += A->V[itot].Re * Lorentzian.Im + A->V[itot].Im * Lorentzian.Re;
          itot++;
        }
      }
    }
  }
  return 0;
}

/*  val needs to be allocated befor the function is called
 *
 */
uint32_t BlockListOfPolesEigenValues(const BlockListOfPoles * const A, double * val)
{
  int i;
  if(A->QComplex==0)
  {
    for(i=0;i<A->N;i++)
    {
      val[i] = A->VRe[i*(1+A->BlockSize*A->BlockSize)];
    }
  }
  else
  {
    for(i=0;i<A->N;i++)
    {
      val[i] = A->V[i*(1+A->BlockSize*A->BlockSize)].Re;
    }
  }
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
  return BlockListOfPolesToBlockTridiagonalMatrix(A, B, NBlockMax, SingularValue, NTriMaxNoLimit);
}

/*  NTriMax is a ceiling on the number of blocks of B, honoured by stopping the block
 *  Lanczos there rather than by building the whole chain and cutting it. See
 *  DiagonalMatrixBlockTridiagonalize for why this is not what NBlockMax does.
 */
uint32_t BlockListOfPolesToBlockTridiagonalMatrix(BlockListOfPoles const A, BlockTriDiagonalMatrixType *const B, int const NBlockMax, double const SingularValue, uint32_t const NTriMax)
{
  double *alpha;
  int k,i,ialpha,iN;
  char name[]="temp";
  CompactMatrixType betaisqr,betai,beta,T;
  if(A.N==0)
  { /*  No poles, so there is no chain and nothing to tri-diagonalize: what is left is
     *  the offset a0. The block Lanczos below must not be handed a starting space of
     *  dimension zero -- CompactMatrixBlockBandDiagonalize would allocate for no
     *  blocks and then report one. */
    B->N=0;
    B->BlockSize=A.BlockSize;
    B->QComplex=A.QComplex;
    if(InitBlockTridiagonalMatrix(B))
    {
      printf("InitBlockTridiagonalMatrix failed in BlockListOfPolesToBlockTridiagonalMatrix\n");
      fflush(stdout);
      return 1;
    }
    B->mu=A.mu;
    strncpy(B->Name, A.Name, MAXStrLength);
    B->Name[MAXStrLength-1]='\0';
    if(A.QComplex) memcpy(B->a,   A.a0,   Sqr(A.BlockSize)*sizeof(Complex16));
    else           memcpy(B->aRe, A.a0Re, Sqr(A.BlockSize)*sizeof(double));
    return 0;
  }
  
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
  CompactMatrixTranspose(&beta);
  DiagonalMatrixBlockTridiagonalize(alpha, &beta, &T, B, NBlockMax, SingularValue, NTriMax);
  B->mu=A.mu;
  free(alpha);
  FreeCompactMatrix(&beta);
  FreeCompactMatrix(&T);
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

/*  A block list of poles to the block double tri-diagonal form, through the Anderson
 *  form, as the scalar type does.
 *
 *  The detour is not avoidable: a list of poles holds the poles of G, while the two
 *  chains of the double tri-diagonal form are the two halves of the hybridization and
 *  sit in the denominator. Getting from one to the other means inverting, which is
 *  exactly what the Anderson form is: BlockListOfPolesToBlockAndersonMatrix
 *  tri-diagonalizes the poles and then reads the bath off the chain.
 *
 *  Note this is not the route to the natural orbital form. That one splits the poles of
 *  G directly and adds the two resulting Green's functions, so it never needs the
 *  inversion -- see BlockListOfPolesToBlockNaturalImpurityOrbitalMatrix, which looks
 *  deceptively similar and is a different calculation.
 *
 *  B is not allocated on entry.
 */
uint32_t BlockListOfPolesToBlockDoubleTriDiagonalMatrix(BlockListOfPoles const *A, BlockDoubleTriDiagonalMatrixType * const B, double const kbT, double const SingularValue, double const MinimalPoleDistance)
{
  BlockAndersonMatrixType AAnd;
  /*  the +1 mirrors the l -> a cell of ConvertBlockMatrix, which carries a TODO saying
   *  it should not be needed */
  if (BlockListOfPolesToBlockAndersonMatrix(*A, &AAnd, A->N + 1, SingularValue))
  {
    printf("BlockListOfPolesToBlockAndersonMatrix failed in BlockListOfPolesToBlockDoubleTriDiagonalMatrix\n");
    fflush(stdout);
    return 1;
  }
  if (BlockAndersonMatrixToBlockDoubleTriDiagonalMatrix(&AAnd, B, kbT, SingularValue, MinimalPoleDistance))
  {
    printf("BlockAndersonMatrixToBlockDoubleTriDiagonalMatrix failed in BlockListOfPolesToBlockDoubleTriDiagonalMatrix\n");
    fflush(stdout);
    FreeBlockAndersonMatrix(&AAnd);
    return 1;
  }
  FreeBlockAndersonMatrix(&AAnd);
  return 0;
}

/*  on entry A contains the Green's function in a BlockListOfPoles format.
 *  on exit B contains the Green's function in a NaturalImpurityOrbital format.
 *
 */
uint32_t BlockListOfPolesToBlockNaturalImpurityOrbitalMatrix(BlockListOfPoles const *A, BlockNaturalImpurityOrbitalMatrixType * const B, double const kbT, double const SingularValue, uint32_t const NBlockMaxVal, uint32_t const NBlockMaxCon)
{
  return BlockListOfPolesToBlockNaturalImpurityOrbitalMatrix(A, B, kbT, SingularValue, NBlockMaxVal, NBlockMaxCon, NTriMaxNoLimit, NTriMaxNoLimit);
}

/*  NTriMaxVal and NTriMaxCon are ceilings on the number of blocks of the two chains, each
 *  honoured by stopping that chain's block Lanczos there. TMT is built from the chains
 *  afterwards and so comes out right for the chains that are kept.
 */
uint32_t BlockListOfPolesToBlockNaturalImpurityOrbitalMatrix(BlockListOfPoles const *A, BlockNaturalImpurityOrbitalMatrixType * const B, double const kbT, double const SingularValue, uint32_t const NBlockMaxVal, uint32_t const NBlockMaxCon, uint32_t const NTriMaxVal, uint32_t const NTriMaxCon)
{
  BlockListOfPoles Acon, Aval;
  B->mu=A->mu;
  B->Epsilon=kbT;
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
  if(Acon.QComplex)
  {
    memset(Acon.a0, 0, Acon.BlockSize*Acon.BlockSize*sizeof(Complex16));
  }
  else
  {
    memset(Acon.a0Re, 0, Acon.BlockSize*Acon.BlockSize*sizeof(double));
  }
  if(Aval.QComplex)
  {
    memset(Aval.a0, 0, Aval.BlockSize*Aval.BlockSize*sizeof(Complex16));
  }
  else
  {
    memset(Aval.a0Re, 0, Aval.BlockSize*Aval.BlockSize*sizeof(double));
  }
  BlockListOfPolesToBlockTridiagonalMatrix(Acon, &(B->con), NBlockMaxCon, SingularValue, NTriMaxCon);
  BlockListOfPolesToBlockTridiagonalMatrix(Aval, &(B->val), NBlockMaxVal, SingularValue, NTriMaxVal);
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
uint32_t BlockListOfPolesToBlockNaturalImpurityOrbitalMatrixReduceDimension(BlockListOfPoles const *A, BlockNaturalImpurityOrbitalMatrixType * const B, CompactMatrixType const *AT, double const kbT, double const SingularValue, uint32_t const NBlockMaxVal, uint32_t const NBlockMaxCon)
{
  CompactMatrixType Aa0;
  BlockListOfPoles Acon, Aval;
  B->mu=A->mu;
  B->Epsilon=kbT;
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
  BlockListOfPolesToBlockTridiagonalMatrixReduceDimension(Acon, &(B->con), AT, NBlockMaxCon+1, SingularValue);
  BlockListOfPolesToBlockTridiagonalMatrixReduceDimension(Aval, &(B->val), AT, NBlockMaxVal+1, SingularValue);
  SetTMTFromValandConInBlockNaturalImpurityOrbitalMatrixType(B);
  return 0;
}


/*
 *
 * Converts a BlockListofpoles to blockandersonmatrix, using a blocktridiagonalmatrix as intermediate
 *
 *
 */
uint32_t BlockListOfPolesToBlockAndersonMatrix(BlockListOfPoles const A, BlockAndersonMatrixType *const B, int const NBlockMax, double const SingularValue)
{
  if(A.N==0)
  {
    B->N=0;
    B->BlockSize = A.BlockSize;
    B->QComplex = A.QComplex;
    InitBlockAndersonMatrix(B);
    B->mu = A.mu;
    strncpy(B->Name, A.Name, MAXStrLength);
    B->Name[MAXStrLength-1]='\0';
    /*  a0Re and a0 are the pointers to the offset block, not the block itself, and A is
     *  passed by value: &(A.a0Re) is the address of that pointer on the stack, so this
     *  used to copy BlockSize^2 doubles from an eight byte object and leave G = a0 as
     *  whatever happened to follow it */
    if(A.QComplex==0)
    {
      memcpy(B->aRe, A.a0Re, Sqr(A.BlockSize)*sizeof(double));
    }
    else
    {
      memcpy(B->a, A.a0, Sqr(A.BlockSize)*sizeof(Complex16));
    }
    return 0;
  }
  BlockTriDiagonalMatrixType ATri;
  /*  this used to hard code 1000*DBL_EPSILON and ignore its own SingularValue, so the
   *  option reached the function and stopped there */
  if(BlockListOfPolesToBlockTridiagonalMatrix(A, &ATri, NBlockMax, SingularValue))
  {
    printf("BlockListOfPolesToBlockTridiagonalMatrix failed in BlockListOfPolesToBlockAndersonMatrix\n");
    return 1;
  }
  BlockTridiagonalMatrixChop(&ATri,10*DBL_EPSILON);
  if(BlockTriDiagonalMatrixToBlockAndersonMatrix(ATri,B))
  {
    printf("BlockTriDiagonalMatrixToBlockAndersonMatrix failed in BlockListOfPolesToBlockAndersonMatrix\n");
    FreeBlockTridiagonalMatrix(&ATri);
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


/*  Given a Green's function in a block list of poles representation
 *    G = A0 + sum_i B_i / (w - a_i)
 *  the one particle density matrix of the block is
 *    rho = sum_i f(a_i, mu, kbT) B_i
 *  with f the Fermi function. A pole exactly at mu is counted half.
 *  The constant offset A0 carries no spectral weight and does not contribute to rho.
 *  rho is a matrix of size BlockSize by BlockSize
 */
uint32_t BlockListOfPolesToDensityMatrix(BlockListOfPoles const * const A, CompactMatrixType * const rho, double const mu, double const kbT)
{
  uint32_t i, k, itot=0, NBlock=Sqr(A->BlockSize);
  double weigth;
  rho->Ni=A->BlockSize;
  rho->Nj=A->BlockSize;
  rho->QComplex=A->QComplex;
  InitCompactMatrix(rho);
  rho->mu=mu;
  strncpy(rho->Name, A->Name, MAXStrLength);
  if(A->QComplex)
  {
    for(i=0;i<A->N;i++)
    {
      weigth = FermiFunction(A->V[itot].Re, mu, kbT);
      itot++;
      for(k=0;k<NBlock;k++)
      {
        rho->V[k].Re += weigth * A->V[itot].Re;
        rho->V[k].Im += weigth * A->V[itot].Im;
        itot++;
      }
    }
  }
  else
  {
    for(i=0;i<A->N;i++)
    {
      weigth = FermiFunction(A->VRe[itot], mu, kbT);
      itot++;
      for(k=0;k<NBlock;k++)
      {
        rho->VRe[k] += weigth * A->VRe[itot];
        itot++;
      }
    }
  }
  return 0;
}


/*  =========================================================================
 *  BlockAndersonMatrixType to __
 *  ========================================================================= */



/*
 *
 */
uint32_t BlockAndersonMatrixToSpectra(double Energy, const BlockAndersonMatrixType * const A, SpectraType *G, double Emin, double Emax, uint32_t N)
{
  BlockListOfPoles ALP;
  BlockAndersonMatrixToBlockListOfPoles(A, &ALP, 1000*DBL_EPSILON, false);
  BlockListOfPolesToSpectra(Energy, &ALP, G, Emin, Emax, N);
  FreeBlockListOfPoles(&ALP);
  return 0;
}


/*
 *
 */
uint32_t BlockAndersonMatrixToSpectraWithBroadening(double Energy, const BlockAndersonMatrixType * const A, SpectraType *G, double Emin, double Emax, uint32_t N, double const Gamma)
{
  int i,j,iE;
  double omega;
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
  #pragma omp parallel
  {
    CompactMatrixType Gw;
    #pragma omp for private(omega,i,j)
    for(iE=0;iE<=N;iE++)
    {
      omega = ((double)iE/(double)N)*(Emax-Emin)+Emin-Energy;
      EvaluateBlockAndersonMatrix(A, &Gw, omega, Gamma);
      for(i=0;i<A->BlockSize;i++)
      {
        for(j=0;j<A->BlockSize;j++)
        {
          G[i*A->BlockSize+j].Re[iE] = Gw.V[i*A->BlockSize+j].Re;
          G[i*A->BlockSize+j].Im[iE] = Gw.V[i*A->BlockSize+j].Im;
        }
      }
      FreeCompactMatrix(&Gw);
    }
  }
  return 0;
}


uint32_t BlockAndersonMatrixEigenValues(const BlockAndersonMatrixType * const A, double * val)
{
  //TODO all
  printf("Not implemented BlockAndersonMatrixEigenValues\n");
  return 1;
}
uint32_t BlockAndersonMatrixEigenSystem(const BlockAndersonMatrixType * const A, double * val, CompactMatrixType * fun)
{
  
  printf("Not implemented BlockAndersonMatrixEigenSystem\n");
  return 1;
}



uint32_t BlockAndersonMatrixToBlockListOfPoles(BlockAndersonMatrixType const * A, BlockListOfPoles * const G, double const MinimalPoleDistance, bool relativeMinimumPoleDistance)
{
  if(A->N==0)
  {
    G->N=0;
    G->BlockSize=A->BlockSize;
    G->QComplex=A->QComplex;
    InitBlockListOfPoles(G);
    G->mu=A->mu;
    if(A->QComplex==0)
    {
      memcpy(G->a0Re, A->aRe, Sqr(A->BlockSize)*sizeof(double));
    }
    else
    {
      memcpy(G->a0, A->a, Sqr(A->BlockSize)*sizeof(Complex16));
    }
    return 0;
  }
  uint32_t i,j;
  CompactMatrixType Afull, vec;
  if(BlockAndersonMatrixToCompactMatrix(A,&Afull))
  {
    printf("BlockTriDiagonalMatrixToCompactMatrix failed in BlockTriDiagonalMatrixToBlockListOfPoles\n");
    return 1;
  }
  vec.Ni=A->BlockSize;
  vec.Nj=Afull.Nj;
  vec.QComplex=A->QComplex;
  InitCompactMatrix(&vec);
  if(A->QComplex)
  {
    for(i=0;i<A->BlockSize;i++)
    {
      for(j=0;j<A->BlockSize;j++)
      {
        vec.V[i*vec.Nj+j] = A->b[i*A->BlockSize+j];
      }
    }
  }
  else
  {
    for(i=0;i<A->BlockSize;i++)
    {
      for(j=0;j<A->BlockSize;j++)
      {
        vec.VRe[i*vec.Nj+j] = A->bRe[i*A->BlockSize+j];
      }
    }
  }
  if(HermitianCompactMatrixToBlockListOfPoles(&Afull, G, &vec))
  {
    printf("HermitianCompactMatrixToListOfPoles failed in BlockTriDiagonalMatrixToBlockListOfPoles\n");
    return 1;
  }
  FreeCompactMatrix(&Afull);
  FreeCompactMatrix(&vec);
  // set A0
  if(A->QComplex)
  {
    memcpy(G->a0,A->a,G->BlockSize*G->BlockSize*sizeof(Complex16));
  }
  else
  {
    memcpy(G->a0Re,A->aRe,G->BlockSize*G->BlockSize*sizeof(double));
  }
  //
  G->mu=A->mu;
  // merge degenerate poles
  BlockListOfPolesMergeDegeneratePoles(G, MinimalPoleDistance, relativeMinimumPoleDistance);
  return 0;
}
/*  On entry A holds a Green's function as a BlockAndersonMatrix,
 *
 *    G(w)     = A0 + B0^* . ( w - A1 - Delta(w) )^-1 . B0^T
 *    Delta(w) = sum_k B_k . ( w - A_{k+1} )^-1 . B_k^+          k = 1 ... A->N-1
 *
 *  On exit Delta holds everything that sits in that denominator, A1 + Delta(w), as a
 *  block list of poles: a0 = A1, and one pole per bath state.
 *
 *  This is the block twin of AndersonMatrixInverseToListOfPoles, which for the scalar
 *  type sets a0 = a[1] and lists the bath as b_i^2 at a_{i+2}. Here a bath block is a
 *  matrix rather than a number, so each is diagonalized on its own,
 *  A_{k+1} = U_k L_k U_k^+, and contributes BlockSize poles at its eigenvalues:
 *
 *    Delta(w) = sum_{k,m} v_{k,m} v_{k,m}^+ / ( w - l_{k,m} ),   v_{k,m} = (B_k U_k)[:,m]
 *
 *  Every weight is therefore a rank one Hermitian matrix, and that is the point rather
 *  than a shortcoming. Summing the rank one weights over a degenerate group of levels
 *  gives B_k P B_k^+ with P the spectral projector of that group, which does not depend
 *  on the arbitrary basis LAPACK picks inside it, while each term separately does. So a
 *  caller that needs the symmetry of a degenerate multiplet back merges the degenerate
 *  poles, which adds their weight matrices -- spectral weight being the only thing that
 *  may be added. See BlockListOfPolesMergeDegeneratePoles.
 *
 *  The bath blocks are uncoupled, so they are diagonalized one at a time, BlockSize by
 *  BlockSize: O(NBath * BlockSize^3), and the full matrix is never assembled. Do not
 *  confuse this with BlockAndersonMatrixToBlockListOfPoles, which builds the whole
 *  matrix dense, diagonalizes it at O((N*BlockSize)^3) and returns the poles of G. The
 *  poles of G are the right source for the natural orbital form, which adds two Green's
 *  functions; the poles of the hybridization are the right source for the double
 *  tri-diagonal form, whose two chains sit in the denominator.
 *
 *  The poles come out ordered by energy, which is what merging degenerate ones expects.
 *
 *  Delta is not allocated on entry.
 */
uint32_t BlockAndersonMatrixInverseToBlockListOfPoles(BlockAndersonMatrixType const * A, BlockListOfPoles * const Delta)
{
  uint32_t k, m, i, j, BS, BS2, NBath, ipole, base;
  char name[] = "bath block";
  CompactMatrixType Abath, Bk, fun, Vk;
  double *val;

  BS    = A->BlockSize;
  BS2   = BS * BS;
  NBath = (A->N == 0) ? 0 : A->N - 1;

  Delta->N         = NBath * BS;
  Delta->BlockSize = BS;
  Delta->QComplex  = A->QComplex;
  if (InitBlockListOfPoles(Delta))
  {
    printf("InitBlockListOfPoles failed in BlockAndersonMatrixInverseToBlockListOfPoles\n");
    fflush(stdout);
    return 1;
  }
  strncpy(Delta->Name, A->Name, MAXStrLength);
  Delta->Name[MAXStrLength-1] = '\0';
  Delta->mu = A->mu;

  /*  a0 is the impurity block A1. A matrix of dimension zero does not have one, and
   *  InitBlockListOfPoles has already left a0 at zero for that case */
  if (A->N != 0)
  {
    if (A->QComplex) memcpy(Delta->a0,   &(A->a[BS2]),   BS2 * sizeof(Complex16));
    else             memcpy(Delta->a0Re, &(A->aRe[BS2]), BS2 * sizeof(double));
  }
  if (NBath == 0) return 0;

  if ((val = (double*)malloc(BS * sizeof(double))) == NULL)
  {
    printf("malloc for val failed in BlockAndersonMatrixInverseToBlockListOfPoles\n");
    fflush(stdout);
    FreeBlockListOfPoles(Delta);
    return 1;
  }
  /*  Abath and Bk are views on A's own storage and are never freed */
  Abath.Name = name;  Abath.Ni = BS;  Abath.Nj = BS;  Abath.QComplex = A->QComplex;
  Abath.a0 = 0;  Abath.b0 = 1;  Abath.mu = 0;  Abath.OnLuaStack = false;
  Bk = Abath;

  ipole = 0;
  for (k = 0; k < NBath; k++)
  {
    /*  the bath block A_{k+1} is A->a[k+2], the hopping B_k into it is A->b[k+1] */
    if (A->QComplex)
    {
      Abath.V = &(A->a[(k+2)*BS2]);
      Bk.V    = &(A->b[(k+1)*BS2]);
    }
    else
    {
      Abath.VRe = &(A->aRe[(k+2)*BS2]);
      Bk.VRe    = &(A->bRe[(k+1)*BS2]);
    }
    /*  CompactMatrixEigensystem returns fun with fun^* . A . fun^T diagonal, so the
     *  eigenvectors are the rows of fun and U_k is fun^T */
    if (CompactMatrixEigensystem(&Abath, &fun, val))
    {
      printf("CompactMatrixEigensystem failed for bath block %u in BlockAndersonMatrixInverseToBlockListOfPoles\n", k);
      fflush(stdout);
      free(val);
      FreeBlockListOfPoles(Delta);
      return 1;
    }
    if (CompactMatrixMultiply('N', 'T', &Bk, &fun, &Vk))
    {
      printf("CompactMatrixMultiply failed for bath block %u in BlockAndersonMatrixInverseToBlockListOfPoles\n", k);
      fflush(stdout);
      FreeCompactMatrix(&fun);
      free(val);
      FreeBlockListOfPoles(Delta);
      return 1;
    }
    /*  column m of Vk is v_{k,m}, and the pole at l_{k,m} carries v_{k,m} v_{k,m}^+ */
    for (m = 0; m < BS; m++)
    {
      base = ipole * (1 + BS2);
      if (A->QComplex)
      {
        Delta->V[base].Re = val[m];
        Delta->V[base].Im = 0;
        for (i = 0; i < BS; i++)
        {
          for (j = 0; j < BS; j++)
          {
            Delta->V[base + 1 + i*BS + j].Re =  Vk.V[i*BS+m].Re * Vk.V[j*BS+m].Re + Vk.V[i*BS+m].Im * Vk.V[j*BS+m].Im;
            Delta->V[base + 1 + i*BS + j].Im = -Vk.V[i*BS+m].Re * Vk.V[j*BS+m].Im + Vk.V[i*BS+m].Im * Vk.V[j*BS+m].Re;
          }
        }
      }
      else
      {
        Delta->VRe[base] = val[m];
        for (i = 0; i < BS; i++)
        {
          for (j = 0; j < BS; j++)
          {
            Delta->VRe[base + 1 + i*BS + j] = Vk.VRe[i*BS+m] * Vk.VRe[j*BS+m];
          }
        }
      }
      ipole++;
    }
    FreeCompactMatrix(&fun);
    FreeCompactMatrix(&Vk);
  }
  free(val);
  OrderBlockListOfPoles(Delta);
  return 0;
}

uint32_t BlockAndersonMatrixToBlockListOfPolesVariableGrid(BlockAndersonMatrixType const AAnd, BlockListOfPoles * const B, const double * const Grid, const uint32_t NGridPoints)
{
  //TODO: all
  printf("Not implemented BlockAndersonMatrixToBlockListOfPoles\n");
  return 1;
}
uint32_t BlockAndersonMatrixToBlockTridiagonalMatrix(BlockAndersonMatrixType const * A, BlockTriDiagonalMatrixType * const B, const double SingularValue, double const MinimalPoleDistance, bool const relativeDistance)
{
  BlockListOfPoles G;
  BlockAndersonMatrixToBlockListOfPoles(A, &G, MinimalPoleDistance, relativeDistance);
  //TODO: here there should not be a +1 for the dimension
  //TODO: test what is wrong here
  BlockListOfPolesToBlockTridiagonalMatrix(G, B, A->N+1, SingularValue);
  FreeBlockListOfPoles(&G);
  return 0;
}
/*  Split the bath of a block Anderson matrix at the chemical potential and put the two
 *  halves in the two chains of a block double tri-diagonal matrix,
 *
 *    G(w) = A0 + B0^* . ( w - Gval(w) - Gcon(w) )^-1 . B0^T
 *
 *  The chains are the two halves of the hybridization as seen from the observed block,
 *  so they sit in the denominator; that is what makes this different from the natural
 *  orbital form, which adds two Green's functions and is therefore built from the poles
 *  of G. Here the source is the pole list of the hybridization,
 *  BlockAndersonMatrixInverseToBlockListOfPoles.
 *
 *  Four steps, and the order matters:
 *
 *   1. the bath as poles, one rank one weight per bath state;
 *   2. merge the degenerate ones. This restores the symmetry of a degenerate multiplet:
 *      the rank one weights of a degenerate group add to B P B^+ with P the spectral
 *      projector, which does not depend on the basis the eigensolver picked inside the
 *      group, while each weight separately does;
 *   3. the Fermi split, f to the valence chain and 1-f to the conduction chain. It comes
 *      after the merge so that degenerate partners get an identical f by construction
 *      and the split cannot undo the symmetry the merge just restored;
 *   4. a block Lanczos on each half.
 *
 *  MinimalPoleDistance has to stay far below any physical splitting: micro eV Zeeman
 *  splittings on a band of tens of eV are the case this type exists for. The error it
 *  is meant to hide is the O(eps * |A_k|) of the individual BlockSize by BlockSize
 *  diagonalizations, which does not grow with the number of poles, so the default is a
 *  few eps times the band width. Pass a value <= 0 to get that default.
 *
 *  The impurity block is shared by the two branches: val.a[0] + con.a[0] is what enters
 *  the denominator and each is given half of it. How it is split is a gauge with no
 *  canonical value, so equal halves is a convention and not an invariant.
 *
 *  A chain can come out empty, with mu below or above the whole bath, and both can with
 *  no bath at all. B is not allocated on entry.
 */
uint32_t BlockAndersonMatrixToBlockDoubleTriDiagonalMatrix(BlockAndersonMatrixType const * A, BlockDoubleTriDiagonalMatrixType * const B, const double kbT, double const SingularValue, double const MinimalPoleDistance)
{
  uint32_t BS, BS2, NBath, i, NBlockMax;
  BlockListOfPoles Delta, Dval, Dcon;

  BS    = A->BlockSize;
  BS2   = BS * BS;
  NBath = (A->N == 0) ? 0 : A->N - 1;

  if (BlockAndersonMatrixInverseToBlockListOfPoles(A, &Delta))
  {
    printf("BlockAndersonMatrixInverseToBlockListOfPoles failed in BlockAndersonMatrixToBlockDoubleTriDiagonalMatrix\n");
    fflush(stdout);
    return 1;
  }
  if (Delta.N != 0)
  {
    /*  a negative MinimalPoleDistance asks BlockListOfPolesMergeDegeneratePoles for its
     *  standard value, which is the only place that decides it */
    if (BlockListOfPolesMergeDegeneratePoles(&Delta, MinimalPoleDistance, false))
    {
      printf("BlockListOfPolesMergeDegeneratePoles failed in BlockAndersonMatrixToBlockDoubleTriDiagonalMatrix\n");
      fflush(stdout);
      FreeBlockListOfPoles(&Delta);
      return 1;
    }
  }
  /*  CopyBlockListOfPoles* copies a0 along and scales only the weights */
  if (CopyBlockListOfPolesVal(Delta, &Dval, A->mu, kbT) || CopyBlockListOfPolesCon(Delta, &Dcon, A->mu, kbT))
  {
    printf("CopyBlockListOfPolesVal/Con failed in BlockAndersonMatrixToBlockDoubleTriDiagonalMatrix\n");
    fflush(stdout);
    FreeBlockListOfPoles(&Delta);
    return 1;
  }
  FreeBlockListOfPoles(&Delta);
  /*  half of the impurity block into each head. BlockListOfPolesToBlockTridiagonalMatrix
   *  writes a pole list's a0 into a[0] of the chain, so this needs no special casing */
  for (i = 0; i < BS2; i++)
  {
    if (A->QComplex)
    {
      Dval.a0[i].Re *= 0.5;  Dval.a0[i].Im *= 0.5;
      Dcon.a0[i].Re *= 0.5;  Dcon.a0[i].Im *= 0.5;
    }
    else
    {
      Dval.a0Re[i] *= 0.5;
      Dcon.a0Re[i] *= 0.5;
    }
  }
  /*  a level outside the +-40 kbT window has an exactly zero weight in one of the two
   *  halves; dropping those is exact, and it is what makes an absent branch a pole list
   *  with no poles */
  BlockListOfPolesRemoveZeroWeightPoles(&Dval);
  BlockListOfPolesRemoveZeroWeightPoles(&Dcon);
  /*  a chain spans at most the whole bath, NBath*BlockSize states. Blocks shrink where
   *  the block Lanczos deflates, so in the limit that is NBath*BlockSize blocks of
   *  dimension one; the routine stops of its own accord once the Krylov space is
   *  exhausted, so this bounds an allocation and does not truncate anything. */
  NBlockMax = NBath * BS;
  if (BlockListOfPolesToBlockTridiagonalMatrix(Dval, &(B->val), NBlockMax, SingularValue))
  {
    printf("BlockListOfPolesToBlockTridiagonalMatrix failed for val in BlockAndersonMatrixToBlockDoubleTriDiagonalMatrix\n");
    fflush(stdout);
    FreeBlockListOfPoles(&Dval);
    FreeBlockListOfPoles(&Dcon);
    return 1;
  }
  if (BlockListOfPolesToBlockTridiagonalMatrix(Dcon, &(B->con), NBlockMax, SingularValue))
  {
    printf("BlockListOfPolesToBlockTridiagonalMatrix failed for con in BlockAndersonMatrixToBlockDoubleTriDiagonalMatrix\n");
    fflush(stdout);
    FreeBlockTridiagonalMatrix(&(B->val));
    FreeBlockListOfPoles(&Dval);
    FreeBlockListOfPoles(&Dcon);
    return 1;
  }
  FreeBlockListOfPoles(&Dval);
  FreeBlockListOfPoles(&Dcon);

  B->BlockSize = BS;
  B->QComplex  = A->QComplex;
  B->mu        = A->mu;
  B->kBT       = kbT;
  strncpy(B->Name, A->Name, MAXStrLength);
  B->Name[MAXStrLength-1] = '\0';
  /*  a0 and b0 of the matrix itself. A matrix of dimension zero has an a[0] but no
   *  b[0]; b0 then stays zero, which is what says nothing couples to it. */
  B->a0.QComplex = A->QComplex;  B->a0.Ni = BS;  B->a0.Nj = BS;
  B->b0.QComplex = A->QComplex;  B->b0.Ni = BS;  B->b0.Nj = BS;
  if (InitCompactMatrix(&(B->a0)) || InitCompactMatrix(&(B->b0)))
  {
    printf("InitCompactMatrix failed in BlockAndersonMatrixToBlockDoubleTriDiagonalMatrix\n");
    fflush(stdout);
    return 1;
  }
  if (A->QComplex)
  {
    memcpy(B->a0.V, A->a, BS2 * sizeof(Complex16));
    if (A->N != 0) memcpy(B->b0.V, A->b, BS2 * sizeof(Complex16));
  }
  else
  {
    memcpy(B->a0.VRe, A->aRe, BS2 * sizeof(double));
    if (A->N != 0) memcpy(B->b0.VRe, A->bRe, BS2 * sizeof(double));
  }
  return 0;
}
uint32_t BlockAndersonMatrixToBlockNaturalImpurityOrbitalMatrix(BlockAndersonMatrixType const * A, BlockNaturalImpurityOrbitalMatrixType * const B, const double kbT, double const SingularValue, double const MinimalPoleDistance, bool const relativeDistance)
{
  BlockListOfPoles G;
  BlockAndersonMatrixToBlockListOfPoles(A, &G, MinimalPoleDistance, relativeDistance);
  //TODO: here there should not be a +1 for the dimension
  //TODO: test what is wrong here
  
  BlockListOfPolesToBlockNaturalImpurityOrbitalMatrix(&G, B, kbT, SingularValue, A->N+1, A->N+1);
  
  FreeBlockListOfPoles(&G);
  return 0;
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
uint32_t BlockTriDiagonalMatrixToSpectraWithBroadening(const double Energy, BlockTriDiagonalMatrixType const * A, SpectraType *const G, double const Emin, double const Emax, uint32_t const N, double const Gamma)
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

uint32_t BlockTriDiagonalMatrixToSpectra(const double Energy, BlockTriDiagonalMatrixType const * A, SpectraType *const G, double const Emin, double const Emax, uint32_t const N)
{
  BlockListOfPoles ALP;
  /*  a negative MinimalPoleDistance asks for the standard value */
  BlockTriDiagonalMatrixToBlockListOfPoles(A, &ALP, -1, false);
  BlockListOfPolesToSpectra(Energy, &ALP, G, Emin, Emax, N);
  FreeBlockListOfPoles(&ALP);
  return 0;
}


// Given a matrix in BlockTriDiagonal form calculate the Eigenvalues
// val is allocated with a minimal length of A->N * A->BlockSize
uint32_t BlockTriDiagonalMatrixEigenValues(const BlockTriDiagonalMatrixType * const A, double * val)
{
  uint32_t j,k,iBlock,jBlock;
  char const *UPLO="L";
  char const *VECT="N";
  const uint32_t N=A->N*A->BlockSize;
  const uint32_t KD = 2*A->BlockSize-1;
  const uint32_t LDAB = 2*A->BlockSize;
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

uint32_t BlockTriDiagonalMatrixToBlockListOfPoles(BlockTriDiagonalMatrixType const * A, BlockListOfPoles * const G, double const MinimalPoleDistance, bool RelativeDistance)
{
  if(A->N==0)
  {
    G->N=0;
    G->QComplex=A->QComplex;
    G->BlockSize=A->BlockSize;
    InitBlockListOfPoles(G);
    G->mu = A->mu;
    if(A->QComplex==0)
    {
      memcpy(G->a0Re, A->aRe, Sqr(A->BlockSize)*sizeof(double));
    }
    else
    {
      memcpy(G->a0, A->a, Sqr(A->BlockSize)*sizeof(Complex16));
    }
    return 0;
  }
  uint32_t i,j;
  CompactMatrixType Afull, vec;
  if(BlockTriDiagonalMatrixToCompactMatrix(A,&Afull))
  {
    printf("BlockTriDiagonalMatrixToCompactMatrix failed in BlockTriDiagonalMatrixToBlockListOfPoles\n");
    return 1;
  }
  vec.Ni=A->BlockSize;
  vec.Nj=Afull.Nj;
  vec.QComplex=A->QComplex;
  InitCompactMatrix(&vec);
  if(A->QComplex)
  {
    for(i=0;i<A->BlockSize;i++)
    {
      for(j=0;j<A->BlockSize;j++)
      {
        vec.V[i*vec.Nj+j] = A->b[i*A->BlockSize+j];
      }
    }
  }
  else
  {
    for(i=0;i<A->BlockSize;i++)
    {
      for(j=0;j<A->BlockSize;j++)
      {
        vec.VRe[i*vec.Nj+j] = A->bRe[i*A->BlockSize+j];
      }
    }
  }
  if(HermitianCompactMatrixToBlockListOfPoles(&Afull, G, &vec))
  {
    printf("HermitianCompactMatrixToListOfPoles failed in BlockTriDiagonalMatrixToBlockListOfPoles\n");
    return 1;
  }
  FreeCompactMatrix(&Afull);
  FreeCompactMatrix(&vec);
  // set A0
  if(A->QComplex)
  {
    memcpy(G->a0,A->a,G->BlockSize*G->BlockSize*sizeof(Complex16));
  }
  else
  {
    memcpy(G->a0Re,A->aRe,G->BlockSize*G->BlockSize*sizeof(double));
  }
  //
  G->mu=A->mu;
  /*  merge degenerate poles. A negative MinimalPoleDistance asks
   *  BlockListOfPolesMergeDegeneratePoles for its standard value; this used to compute a
   *  band width here and ignore both of its own parameters */
  return BlockListOfPolesMergeDegeneratePoles(G, MinimalPoleDistance, RelativeDistance);
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
 */
uint32_t BlockTriDiagonalMatrixToBlockAndersonMatrix(BlockTriDiagonalMatrixType const ATri, BlockAndersonMatrixType * const B)
{
  if(ATri.N==0)
  {
    B->N=0;
    B->QComplex=ATri.QComplex;
    B->BlockSize=ATri.BlockSize;
    InitBlockAndersonMatrix(B);
    B->mu = ATri.mu;
    if(ATri.QComplex==0)
    {
      memcpy(B->aRe, ATri.aRe, Sqr(ATri.BlockSize)*sizeof(double));
    }
    else
    {
      memcpy(B->a, ATri.a, Sqr(ATri.BlockSize)*sizeof(Complex16));
    }
    return 0;
  }
  
  if(ATri.N==1)
  {
    B->N=1;
    B->QComplex=ATri.QComplex;
    B->BlockSize=ATri.BlockSize;
    InitBlockAndersonMatrix(B);
    B->mu = ATri.mu;
    if(ATri.QComplex==0)
    {
      memcpy(B->aRe, ATri.aRe, 2*Sqr(ATri.BlockSize)*sizeof(double));
      memcpy(B->bRe, ATri.bRe,   Sqr(ATri.BlockSize)*sizeof(double));
    }
    else
    {
      memcpy(B->a, ATri.a, 2*Sqr(ATri.BlockSize)*sizeof(Complex16));
      memcpy(B->b, ATri.b,   Sqr(ATri.BlockSize)*sizeof(Complex16));
    }
    return 0;
  }

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
  /*  the N == 0 and N == 1 branches above carry mu across; this one used to drop it,
   *  which puts the chemical potential at zero for anything converted through here --
   *  and a Fermi split further downstream then happens in the wrong place */
  B->mu = ATri.mu;
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
      C.V[i].Re = ATri.b[bpos].Re;
      C.V[i].Im = ATri.b[bpos].Im;
      bpos++;
    }
    D.QComplex = B->QComplex;
    D.Ni = B->BlockSize;
    D.Nj = B->BlockSize*(B->N-1);
    InitCompactMatrix(&D);
    for(i = 0; i < B->BlockSize; i++)
    {
      for(j = 0; j < B->BlockSize*(B->N-1); j++)
      {
        D.V[i*B->BlockSize*(B->N-1) + j].Re = Afun.V[j*B->BlockSize*(B->N-1)+i].Re;
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
          B->b[bpos].Re = E.V[i*B->BlockSize*(B->N-1)+k*B->BlockSize+j].Re;
          B->b[bpos].Im = E.V[i*B->BlockSize*(B->N-1)+k*B->BlockSize+j].Im;
          bpos++;
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
          B->bRe[bpos] = E.VRe[i*(B->N-1)*B->BlockSize+k*B->BlockSize+j];
          bpos++;
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

/*  A block tri-diagonal matrix to the block double tri-diagonal form, through the
 *  Anderson form -- the same route the scalar type takes.
 *
 *  Both steps need the bath diagonalized: BlockTriDiagonalMatrixToBlockAndersonMatrix
 *  diagonalizes the chain from the second block onwards and transforms the coupling,
 *  and there is no cheaper exact way to get the poles of a chain's hybridization. So
 *  composing costs one dense O((N*BlockSize)^3) diagonalization, which is what a
 *  direct route would cost too, plus the block-local ones of the second step.
 *
 *  The Anderson form in between has mixed hoppings and a reordered bath -- both a
 *  matter of the basis the eigensolver picks, see ChangeTypeAndersonBlock.Quanty --
 *  and neither matters here: BlockAndersonMatrixToBlockDoubleTriDiagonalMatrix
 *  diagonalizes every bath block again on its way to the hybridization poles.
 */
uint32_t BlockTriDiagonalMatrixToBlockDoubleTriDiagonalMatrix(BlockTriDiagonalMatrixType const ATri, BlockDoubleTriDiagonalMatrixType * const B, const double kBT, double const SingularValue, double const MinimalPoleDistance)
{
  BlockAndersonMatrixType AAnd;
  if (BlockTriDiagonalMatrixToBlockAndersonMatrix(ATri, &AAnd))
  {
    printf("BlockTriDiagonalMatrixToBlockAndersonMatrix failed in BlockTriDiagonalMatrixToBlockDoubleTriDiagonalMatrix\n");
    fflush(stdout);
    return 1;
  }
  if (BlockAndersonMatrixToBlockDoubleTriDiagonalMatrix(&AAnd, B, kBT, SingularValue, MinimalPoleDistance))
  {
    printf("BlockAndersonMatrixToBlockDoubleTriDiagonalMatrix failed in BlockTriDiagonalMatrixToBlockDoubleTriDiagonalMatrix\n");
    fflush(stdout);
    FreeBlockAndersonMatrix(&AAnd);
    return 1;
  }
  FreeBlockAndersonMatrix(&AAnd);
  return 0;
}

uint32_t BlockTriDiagonalMatrixToBlockNaturalImpurityMatrix(BlockTriDiagonalMatrixType const ATri, BlockNaturalImpurityOrbitalMatrixType * const B, const double kBT, double const SingularValue, double const MinimalPoleDistance, bool const relativeDistance)
{
  BlockListOfPoles ABlock;
  if(BlockTriDiagonalMatrixToBlockListOfPoles(&ATri, &ABlock, MinimalPoleDistance, relativeDistance))
  {
    printf("BlockTriDiagonalMatrixToBlockListOfPoles failed in BlockTriDiagonalMatrixToBlockNaturalImpurityMatrix\n");
    return 1;
  }
  BlockListOfPolesToBlockNaturalImpurityOrbitalMatrix(&ABlock, B, kBT, SingularValue, ATri.N, ATri.N);
  FreeBlockListOfPoles(&ABlock);
  return 0;
}




/*  =========================================================================
 *  BlockDoubleTriDiagonalMatrixType to __
 *  ========================================================================= */

/*  A block double tri-diagonal matrix to a block tri-diagonal one, through the Anderson
 *  form, which is what the four inbound cells of this row do in the other direction.
 *
 *  The two forms split G differently and the only thing they agree on is G itself, so
 *  something in between has to hold G plainly. The Anderson form does: the two chains
 *  here are the two halves of the hybridization and sit in the denominator, while a
 *  chain of the tri-diagonal form is the whole of G continued-fraction expanded from the
 *  observed block. Going d -> a puts the hybridization back into one bath, and
 *  BlockAndersonMatrixToBlockTridiagonalMatrix then takes the poles of G and runs a block
 *  Lanczos from the observed block.
 *
 *  Both parameters are needed and neither is decided here. SingularValue is what
 *  BlockListOfPolesToBlockTridiagonalMatrix deflates on, so it belongs to the second
 *  step; MinimalPoleDistance is the merge tolerance of the first, and it is forwarded to
 *  the second as well, which merges the poles of G on its way.
 *
 *  Cost and dimension. The first step is O(N*BlockSize) small diagonalizations, the
 *  second one dense O((N*BlockSize)^3) on the whole Anderson matrix. That matrix is the
 *  larger of the two, since d -> a stores a rank one weight as a full BlockSize by
 *  BlockSize sqrt(W) and so comes out a factor BlockSize wider than the d it came from --
 *  it is the intermediate that sets the cost, not the answer. The chain that comes back
 *  is not the interleaved band reduction the scalar type uses, and is not expected to
 *  have the scalar version's dimension 1 + val.N + con.N; it is as long as the block
 *  Lanczos needs before the Krylov space of the observed block is exhausted.
 *
 *  H is not allocated on entry.
 */
uint32_t BlockDoubleTriDiagonalMatrixToBlockTriDiagonalMatrix(BlockDoubleTriDiagonalMatrixType const * M, BlockTriDiagonalMatrixType * const H, double const SingularValue, double const MinimalPoleDistance, bool const relativeDistance)
{
  BlockAndersonMatrixType AAnd;
  if(BlockDoubleTriDiagonalMatrixToBlockAndersonMatrix(M, &AAnd, MinimalPoleDistance, relativeDistance))
  {
    printf("BlockDoubleTriDiagonalMatrixToBlockAndersonMatrix failed in BlockDoubleTriDiagonalMatrixToBlockTriDiagonalMatrix\n");
    fflush(stdout);
    return 1;
  }
  if(BlockAndersonMatrixToBlockTridiagonalMatrix(&AAnd, H, SingularValue, MinimalPoleDistance, relativeDistance))
  {
    printf("BlockAndersonMatrixToBlockTridiagonalMatrix failed in BlockDoubleTriDiagonalMatrixToBlockTriDiagonalMatrix\n");
    fflush(stdout);
    FreeBlockAndersonMatrix(&AAnd);
    return 1;
  }
  FreeBlockAndersonMatrix(&AAnd);
  /*  the name is the one thing the route loses: the intermediate pole list is created
   *  with a name of its own */
  strncpy(H->Name, M->Name, MAXStrLength);
  H->Name[MAXStrLength-1] = '\0';
  return 0;
}

/*  The bath of a block double tri-diagonal matrix is its two chains side by side, and one
 *  Anderson block is made per merged pole.
 *
 *  The two chains sit in the same denominator,
 *
 *    G(w) = A0 + B0^* . ( w - Gval(w) - Gcon(w) )^-1 . B0^T
 *
 *  and each of them is an ordinary block tri-diagonal response function, so
 *  Gval(w) = val.a[0] + sum_k W_k / (w - e_k) over the poles of its own minor. Adding the
 *  two pole lists gives everything that stands in that denominator: a0 = val.a[0] +
 *  con.a[0], which is the observed block the two branches share, and one pole per bath
 *  state. That is the Anderson denominator, w - A1 - sum_k B_k (w - A_{k+1})^-1 B_k^+, with
 *
 *    A_{k+1} = e_k * I         (BlockSize by BlockSize)
 *    B_k     = sqrt(W_k)       (full, dense, no columns dropped)
 *
 *  since sum_k sqrt(W_k) (w - e_k I)^-1 sqrt(W_k)^+ = sum_k W_k / (w - e_k), exactly.
 *
 *  Why this form and not a smaller one. The dimension of the matrix representing a response
 *  function is an efficiency question; breaking a symmetry is not acceptable at all.
 *  Shifting every bath level of a t2g irep by the same 100 meV is a poor but usable
 *  approximation, shifting one of the three by 1 meV changes magnetic susceptibilities. So
 *  dimension loses. A_{k+1} = e_k * I is what makes a degeneracy structural rather than
 *  numerical: all BlockSize states of a block share one stored number, so nothing
 *  downstream can split them -- not roundoff, not a later shift, not a chop. Packing states
 *  of different energies onto one diagonal would give every state its own number and its
 *  own opportunity to drift.
 *
 *  sqrt(W) is canonical although CompactMatrixSqrt computes it from an arbitrary
 *  eigenbasis: any f(W) = U f(L) U^+ depends only on the spectral projectors of W, so the
 *  arbitrary basis inside a degenerate group cancels exactly. Individual eigenvectors of W
 *  are not canonical, and a factorization into exactly rank(W) columns is defined only up
 *  to a unitary U(rank) -- minimal dimension and canonicality are in genuine conflict.
 *  Zero couplings are therefore kept and the bath comes out N1 * BlockSize states for N1
 *  merged poles; BlockDoubleTriDiagonalMatrixChop removes what carries nothing. Compression
 *  belongs earlier, in placing the poles on a grid, which merges them onto shared energies
 *  before this conversion and makes each grid point one block.
 *
 *  Where the conjugates go. A block matrix couples its head as b0^* ( ) b0^T and every
 *  level below it as b ( ) b^+, see the struct comments; both forms here have that same
 *  head, so a0 and b0 are copied as they are. One level down the two differ, and the pole
 *  list of a chain comes out in the head convention -- HermitianCompactMatrixToBlockListOfPoles
 *  computes b0^* (w-A)^-1 b0^T and stores W = x^* x^T -- while a bath block of the Anderson
 *  form sits below and wants B B^+. The two read the same Hermitian matrix with the index
 *  order swapped, and the Hermitian square root is what crosses between them:
 *  CompactMatrixSqrt returns B = B^+ with B B^+ = W, so sqrt(W) is already the object the
 *  bath wants. Note that the same B in a head slot would give W^T instead.
 *
 *  A block double tri-diagonal matrix carries four independent QComplex flags, one per
 *  part, and its own may be 2 for a mixed one, while an Anderson matrix has a single flag.
 *  So the result is complex when any part is, and the pole list is made complex before the
 *  two halves are added -- which also keeps AddBlockListOfPoles off its mixed branch.
 *
 *  MinimalPoleDistance is forwarded and never decided here: a negative value reaches
 *  BlockListOfPolesMergeDegeneratePoles, the one place that turns it into the standard
 *  10 * DBL_EPSILON * (BandWidth + DBL_EPSILON). Merging is not optional. The rank one
 *  weights of a degenerate group add to a spectral projector, which does not depend on the
 *  basis the eigensolver picked inside the group, while each weight separately does.
 *
 *  H is not allocated on entry.
 */
uint32_t BlockDoubleTriDiagonalMatrixToBlockAndersonMatrix(BlockDoubleTriDiagonalMatrixType const * M, BlockAndersonMatrixType * const H, double const MinimalPoleDistance, bool const relativeDistance)
{
  uint32_t BS, BS2, k, i, stride;
  bool QZero;
  double e;
  char name[] = "block";
  BlockListOfPoles L, Lcon;
  CompactMatrixType W, SqrtW, view;

  BS    = M->BlockSize;
  BS2   = BS * BS;
  /*  a matrix of dimension zero: nothing couples to the observed block and G = a0. The
   *  scalar type has the same early return on b0 == 0 */
  QZero = BlockDoubleTriDiagonalMatrixB0IsZero(M);
  H->BlockSize = BS;
  H->QComplex  = (M->a0.QComplex!=0);
  if(QZero==false)
  {
    if((M->b0.QComplex!=0)||(M->val.QComplex!=0)||(M->con.QComplex!=0)) H->QComplex = 1;
    if(BlockTriDiagonalMatrixToBlockListOfPoles(&(M->val), &L, MinimalPoleDistance, relativeDistance))
    {
      printf("BlockTriDiagonalMatrixToBlockListOfPoles failed for val in BlockDoubleTriDiagonalMatrixToBlockAndersonMatrix\n");
      fflush(stdout);
      return 1;
    }
    if(BlockTriDiagonalMatrixToBlockListOfPoles(&(M->con), &Lcon, MinimalPoleDistance, relativeDistance))
    {
      printf("BlockTriDiagonalMatrixToBlockListOfPoles failed for con in BlockDoubleTriDiagonalMatrixToBlockAndersonMatrix\n");
      fflush(stdout);
      FreeBlockListOfPoles(&L);
      return 1;
    }
    /*  everything the Anderson matrix stores has to be of one type, and a0 or b0 alone
     *  being complex is enough to make it so -- the chains may still be real */
    if((H->QComplex!=0)&&(L.QComplex==0))
    {
      if(MakeRealBlockListOfPolesComplex(&L))
      {
        printf("MakeRealBlockListOfPolesComplex failed in BlockDoubleTriDiagonalMatrixToBlockAndersonMatrix\n");
        fflush(stdout);
        FreeBlockListOfPoles(&L);
        FreeBlockListOfPoles(&Lcon);
        return 1;
      }
    }
    /*  the two a0 blocks add, which is the observed block val.a[0] + con.a[0], and the
     *  poles concatenate */
    if(AddBlockListOfPoles(&L, Lcon))
    {
      printf("AddBlockListOfPoles failed in BlockDoubleTriDiagonalMatrixToBlockAndersonMatrix\n");
      fflush(stdout);
      FreeBlockListOfPoles(&L);
      FreeBlockListOfPoles(&Lcon);
      return 1;
    }
    FreeBlockListOfPoles(&Lcon);
    /*  a level that the Fermi split left in both chains arrives twice and has to be put
     *  back together here. BlockListOfPolesMergeDegeneratePoles orders the list itself */
    if(BlockListOfPolesMergeDegeneratePoles(&L, MinimalPoleDistance, relativeDistance))
    {
      printf("BlockListOfPolesMergeDegeneratePoles failed in BlockDoubleTriDiagonalMatrixToBlockAndersonMatrix\n");
      fflush(stdout);
      FreeBlockListOfPoles(&L);
      return 1;
    }
  }
  H->N = (QZero==true) ? 0 : 1 + L.N;
  if(InitBlockAndersonMatrix(H))
  {
    printf("InitBlockAndersonMatrix failed in BlockDoubleTriDiagonalMatrixToBlockAndersonMatrix\n");
    fflush(stdout);
    if(QZero==false) FreeBlockListOfPoles(&L);
    return 1;
  }
  strncpy(H->Name, M->Name, MAXStrLength);
  H->Name[MAXStrLength-1] = '\0';
  H->mu = M->mu;
  /*  view is a window on H's own storage, so what is added lands in H. H is complex
   *  whenever a0 or b0 is, which is what keeps CompactMatrixAdd from reallocating the
   *  window instead of writing through it */
  view.Name = name;  view.Ni = BS;  view.Nj = BS;  view.QComplex = H->QComplex;
  view.a0 = 0;  view.b0 = 1;  view.mu = 0;  view.OnLuaStack = false;
  if(H->QComplex) view.V = H->a;  else view.VRe = H->aRe;
  CompactMatrixAdd(&view, &(M->a0));            /* InitBlockAndersonMatrix left a[0] zero */
  if(H->N==0)
  {
    return 0;
  }
  if(H->QComplex) view.V = H->b;  else view.VRe = H->bRe;
  SetCompactMatrixToZero(&view);                /* and b[0] the identity */
  CompactMatrixAdd(&view, &(M->b0));
  /*  W is a window on the pole list, one weight at a time, and is never freed */
  W.Name = name;  W.Ni = BS;  W.Nj = BS;  W.QComplex = L.QComplex;
  W.a0 = 0;  W.b0 = 1;  W.mu = 0;  W.OnLuaStack = false;
  /*  A1 is the observed block shared by the two branches */
  if(H->QComplex) view.V = &(H->a[BS2]);  else view.VRe = &(H->aRe[BS2]);
  if(L.QComplex)  W.V    = L.a0;          else W.VRe    = L.a0Re;
  CompactMatrixAdd(&view, &W);
  /*  one bath block per merged pole */
  stride = 1 + BS2;
  for(k=0;k<L.N;k++)
  {
    if(L.QComplex)
    {
      e     = L.V[k*stride].Re;
      W.V   = &(L.V[k*stride+1]);
    }
    else
    {
      e     = L.VRe[k*stride];
      W.VRe = &(L.VRe[k*stride+1]);
    }
    /*  a weight that is not quite positive is reported by CompactMatrixSqrt and clamped;
     *  that is a warning about the input and not a reason to fail the conversion */
    CompactMatrixSqrt(&W, &SqrtW);
    if(H->QComplex)
    {
      for(i=0;i<BS;i++) H->a[(k+2)*BS2 + i*BS + i].Re = e;
      memcpy(&(H->b[(k+1)*BS2]), SqrtW.V, BS2*sizeof(Complex16));
    }
    else
    {
      for(i=0;i<BS;i++) H->aRe[(k+2)*BS2 + i*BS + i] = e;
      memcpy(&(H->bRe[(k+1)*BS2]), SqrtW.VRe, BS2*sizeof(double));
    }
    FreeCompactMatrix(&SqrtW);
  }
  FreeBlockListOfPoles(&L);
  return 0;
}

/*  The natural orbital form adds its two chains where the double tri-diagonal form nests
 *  them -- see the two struct comments -- so there is no rotation from one to the other:
 *  the two split the same G differently. What they do agree on is G, and the form that
 *  holds G plainly is the list of poles, so that is the route. It is the one the scalar
 *  DoubleTriDiagonalMatrixToBADoubleTriDiagonalMatrix takes as well.
 *
 *  Going through the poles also means d -> n splits at mu exactly as l -> n does, both
 *  through CopyBlockListOfPolesVal and CopyBlockListOfPolesCon, so the two agree by
 *  construction rather than by accident.
 *
 *  Note this is the poles of G and not the poles of the hybridization, which is what
 *  d -> a wants: the observed level is a pole of G and not one of the bath. That is the
 *  whole difference between the two readings of the same matrix.
 *
 *  Three parameters, none of them decided here. MinimalPoleDistance is the merge tolerance
 *  of the first step, kBT the temperature of the Fermi split of the second, and
 *  SingularValue what the two block Lanczos runs of the second step deflate on. The chain
 *  lengths are bounded by the number of poles, which is what the scalar version passes
 *  too; it bounds an allocation rather than truncating anything.
 *
 *  H is not allocated on entry.
 */
uint32_t BlockDoubleTriDiagonalMatrixToBlockNaturalImpurityMatrix(BlockDoubleTriDiagonalMatrixType const * M, BlockNaturalImpurityOrbitalMatrixType * const H, double const kBT, double const SingularValue, double const MinimalPoleDistance, bool const relativeDistance)
{
  BlockListOfPoles L;
  if(BlockDoubleTriDiagonalMatrixToBlockListOfPoles(M, &L, MinimalPoleDistance, relativeDistance))
  {
    printf("BlockDoubleTriDiagonalMatrixToBlockListOfPoles failed in BlockDoubleTriDiagonalMatrixToBlockNaturalImpurityMatrix\n");
    fflush(stdout);
    return 1;
  }
  if(BlockListOfPolesToBlockNaturalImpurityOrbitalMatrix(&L, H, kBT, SingularValue, L.N, L.N))
  {
    printf("BlockListOfPolesToBlockNaturalImpurityOrbitalMatrix failed in BlockDoubleTriDiagonalMatrixToBlockNaturalImpurityMatrix\n");
    fflush(stdout);
    FreeBlockListOfPoles(&L);
    return 1;
  }
  FreeBlockListOfPoles(&L);
  strncpy(H->Name, M->Name, MAXStrLength);
  H->Name[MAXStrLength-1] = '\0';
  return 0;
}

/*  Copy one BlockSize by BlockSize block, applying op on the way: 'N' as it is, 'T'
 *  transposed, 'C' conjugated, 'D' conjugate transposed. dst may be complex where src is
 *  real, which is what lets a real chain be written into a complex matrix.
 *
 *  Both are views on somebody else's storage and neither is allocated or freed here.
 */
static uint32_t BlockCopyWithOp(CompactMatrixType * const dst, CompactMatrixType const * const src, char const op)
{
  uint32_t i, j, n;
  double re, im;
  n = src->Ni;
  for(i=0;i<n;i++)
  {
    for(j=0;j<n;j++)
    {
      if(src->QComplex) { re = src->V[i*n+j].Re; im = src->V[i*n+j].Im; }
      else              { re = src->VRe[i*n+j];  im = 0;                }
      if((op=='C')||(op=='D')) im = -im;
      if((op=='T')||(op=='D'))
      {
        if(dst->QComplex) { dst->V[j*n+i].Re = re; dst->V[j*n+i].Im = im; }
        else                dst->VRe[j*n+i]  = re;
      }
      else
      {
        if(dst->QComplex) { dst->V[i*n+j].Re = re; dst->V[i*n+j].Im = im; }
        else                dst->VRe[i*n+j]  = re;
      }
    }
  }
  return 0;
}

/*  The whole matrix is block tri-diagonal, with the observed block in the middle.
 *
 *  The poles of G are the eigenvalues of the whole matrix, with weights what the observed
 *  block sees of each eigenvector -- not the two chain spectra, which are the poles of the
 *  hybridization and are what BlockDoubleTriDiagonalMatrixToBlockAndersonMatrix wants.
 *
 *  So this needs the whole matrix diagonalized, and it is cheap. Read it as
 *
 *      val_N ... val_2, val_1, observed, con_1, con_2, ... con_N
 *
 *  and it is an ordinary block tri-diagonal matrix of 1 + val.N + con.N blocks, half
 *  bandwidth 2*BlockSize-1, which BlockTriDiagonalMatrixEigenSystem reduces with
 *  dsbtrd / zhbtrd and diagonalizes. Only the observed block sits in the middle rather
 *  than at an end, which costs nothing here: the weights are read from the middle block
 *  of the eigenvectors instead of the first. (The Anderson route that d -> t takes is
 *  needed only when the observed block has to come first, which a Lanczos chain requires
 *  and an eigenvalue decomposition does not.)
 *
 *  Where the conjugates go. A block matrix couples its head as b0^* ( ) b0^T and every
 *  level below it as b ( ) b^+, see the struct comments, and the two chains of a d enter
 *  the denominator through their heads. So the block coupling the observed block to val_1
 *  is conj(val.b[0]) and not val.b[0]: eliminating the chain then gives
 *  conj(b0) ( ) conj(b0)^+ = b0^* ( ) b0^T, which is what EvaluateBlockTriDiagonalMatrix
 *  computes. Reversing the valence half exchanges the upper and lower blocks on top of
 *  that, so the bonds inside it are transposed as well:
 *
 *    T.b[j]        j = 1 ... val.N-1     val.b[val.N-j]^+     inner val bonds, reversed
 *    T.b[val.N]                          val.b[0]^T           val_1 to the observed block
 *    T.b[val.N+1]                        conj(con.b[0])       observed block to con_1
 *    T.b[val.N+1+k] k = 1 ... con.N-1    con.b[k]             inner con bonds, as stored
 *
 *  Two of those four differ from a plain copy already for a real matrix, and the other two
 *  only for a complex one.
 *
 *  The weights. With the observed block at block index val.N,
 *
 *    G(w) = a0 + b0^* [ (w - H)^-1 ]_{observed} b0^T = a0 + sum_k W_k / (w - e_k)
 *    W_k  = b0^* u_k u_k^+ b0^T   with u_k the part of eigenvector k on that block
 *
 *  which is what a BlockListOfPoles holds: the weight itself, B_i = B'_i B'_i^C, and not
 *  a coupling. sqrt(W) is what turns a pole back into a matrix -- it is what d -> a puts
 *  in the Anderson bath and what BlockListOfPolesToBlockTridiagonalMatrix forms as its
 *  starting block -- and has no business here.
 *
 *  BlockTriDiagonalMatrixEigenSystem returns fun with fun^* H fun^T diagonal, so the rows
 *  of fun are the eigenvectors, the same convention CompactMatrixEigensystem uses. The two
 *  steps below are then the ones HermitianCompactMatrixToBlockListOfPoles takes, so this
 *  list comes out with the same index convention as every other one; that function is not
 *  called itself only because it would diagonalize densely and throw the band away.
 *
 *  Poles of zero weight are kept, as everywhere in the d row; Chop and ReducePoles remove
 *  them. The eigenvalues come out ascending, and the merge orders in any case.
 *
 *  H is not allocated on entry.
 */
uint32_t BlockDoubleTriDiagonalMatrixToBlockListOfPoles(BlockDoubleTriDiagonalMatrixType const * M, BlockListOfPoles * const H, double const MinimalPoleDistance, bool const relativeDistance)
{
  uint32_t BS, BS2, N, Nv, Nc, imp, n, k, i, j, base, stride;
  char name[] = "block";
  BlockTriDiagonalMatrixType T;
  CompactMatrixType fun, vec, funvec, dst, src;
  double *val;

  BS  = M->BlockSize;
  BS2 = BS * BS;
  Nv  = M->val.N;
  Nc  = M->con.N;
  N   = 1 + Nv + Nc;               /* blocks of the matrix whose resolvent is G */
  imp = Nv;                        /* 0 based block index of the observed block */
  n   = N * BS;

  H->BlockSize = BS;
  H->QComplex  = (M->a0.QComplex!=0);
  /*  a matrix of dimension zero: nothing couples to the observed block and G = a0 */
  if(BlockDoubleTriDiagonalMatrixB0IsZero(M))
  {
    H->N = 0;
    if(InitBlockListOfPoles(H))
    {
      printf("InitBlockListOfPoles failed in BlockDoubleTriDiagonalMatrixToBlockListOfPoles\n");
      fflush(stdout);
      return 1;
    }
    if(H->QComplex) memcpy(H->a0,   M->a0.V,   BS2*sizeof(Complex16));
    else            memcpy(H->a0Re, M->a0.VRe, BS2*sizeof(double));
    H->mu = M->mu;
    strncpy(H->Name, M->Name, MAXStrLength);
    H->Name[MAXStrLength-1] = '\0';
    return 0;
  }
  if((M->b0.QComplex!=0)||(M->val.QComplex!=0)||(M->con.QComplex!=0)) H->QComplex = 1;

  /*  the valence chain backwards, then the observed block, then the conduction chain */
  T.N        = N;
  T.BlockSize= BS;
  T.QComplex = H->QComplex;
  if(InitBlockTridiagonalMatrix(&T))
  {
    printf("InitBlockTridiagonalMatrix failed in BlockDoubleTriDiagonalMatrixToBlockListOfPoles\n");
    fflush(stdout);
    return 1;
  }
  T.mu = M->mu;
  dst.Name = name;  dst.Ni = BS;  dst.Nj = BS;  dst.QComplex = T.QComplex;
  dst.a0 = 0;  dst.b0 = 1;  dst.mu = 0;  dst.OnLuaStack = false;
  src = dst;
  /*  the diagonal, a[1] to a[N]. InitBlockTridiagonalMatrix left everything zero */
  src.QComplex = M->val.QComplex;
  for(k=1;k<=Nv;k++)
  {
    if(T.QComplex) dst.V   = &(T.a[k*BS2]);        else dst.VRe = &(T.aRe[k*BS2]);
    if(src.QComplex) src.V = &(M->val.a[(Nv+1-k)*BS2]); else src.VRe = &(M->val.aRe[(Nv+1-k)*BS2]);
    BlockCopyWithOp(&dst, &src, 'N');
  }
  src.QComplex = M->con.QComplex;
  for(k=1;k<=Nc;k++)
  {
    if(T.QComplex) dst.V   = &(T.a[(Nv+1+k)*BS2]); else dst.VRe = &(T.aRe[(Nv+1+k)*BS2]);
    if(src.QComplex) src.V = &(M->con.a[k*BS2]);   else src.VRe = &(M->con.aRe[k*BS2]);
    BlockCopyWithOp(&dst, &src, 'N');
  }
  /*  the observed block is shared: its energy is the plain sum of the two heads.
   *  CompactMatrixAdd adds a real block into a complex one and does not reallocate the
   *  window, T being complex whenever either chain is */
  if(T.QComplex) dst.V = &(T.a[(Nv+1)*BS2]); else dst.VRe = &(T.aRe[(Nv+1)*BS2]);
  src.QComplex = M->val.QComplex;
  if(src.QComplex) src.V = M->val.a; else src.VRe = M->val.aRe;
  CompactMatrixAdd(&dst, &src);
  src.QComplex = M->con.QComplex;
  if(src.QComplex) src.V = M->con.a; else src.VRe = M->con.aRe;
  CompactMatrixAdd(&dst, &src);
  /*  and the bonds, with the conjugates of the table above */
  src.QComplex = M->val.QComplex;
  for(k=1;k<Nv;k++)
  {
    if(T.QComplex) dst.V   = &(T.b[k*BS2]);            else dst.VRe = &(T.bRe[k*BS2]);
    if(src.QComplex) src.V = &(M->val.b[(Nv-k)*BS2]);  else src.VRe = &(M->val.bRe[(Nv-k)*BS2]);
    BlockCopyWithOp(&dst, &src, 'D');
  }
  if(Nv!=0)
  {
    if(T.QComplex) dst.V   = &(T.b[Nv*BS2]);           else dst.VRe = &(T.bRe[Nv*BS2]);
    if(src.QComplex) src.V = M->val.b;                 else src.VRe = M->val.bRe;
    BlockCopyWithOp(&dst, &src, 'T');
  }
  src.QComplex = M->con.QComplex;
  if(Nc!=0)
  {
    if(T.QComplex) dst.V   = &(T.b[(Nv+1)*BS2]);       else dst.VRe = &(T.bRe[(Nv+1)*BS2]);
    if(src.QComplex) src.V = M->con.b;                 else src.VRe = M->con.bRe;
    BlockCopyWithOp(&dst, &src, 'C');
  }
  for(k=1;k<Nc;k++)
  {
    if(T.QComplex) dst.V   = &(T.b[(Nv+1+k)*BS2]);     else dst.VRe = &(T.bRe[(Nv+1+k)*BS2]);
    if(src.QComplex) src.V = &(M->con.b[k*BS2]);       else src.VRe = &(M->con.bRe[k*BS2]);
    BlockCopyWithOp(&dst, &src, 'N');
  }

  if((val = (double*)malloc(n*sizeof(double)))==NULL)
  {
    printf("malloc for val failed in BlockDoubleTriDiagonalMatrixToBlockListOfPoles\n");
    fflush(stdout);
    FreeBlockTridiagonalMatrix(&T);
    return 1;
  }
  if(BlockTriDiagonalMatrixEigenSystem(&T, val, &fun))
  {
    printf("BlockTriDiagonalMatrixEigenSystem failed in BlockDoubleTriDiagonalMatrixToBlockListOfPoles\n");
    fflush(stdout);
    FreeBlockTridiagonalMatrix(&T);
    free(val);
    return 1;
  }
  FreeBlockTridiagonalMatrix(&T);

  /*  b0 in the columns of the observed block, so that funvec = conj(fun) . vec^T holds
   *  b0^* u_k for every eigenvector k */
  vec.Ni = BS;  vec.Nj = n;  vec.QComplex = H->QComplex;
  if(InitCompactMatrix(&vec))
  {
    printf("InitCompactMatrix failed in BlockDoubleTriDiagonalMatrixToBlockListOfPoles\n");
    fflush(stdout);
    FreeCompactMatrix(&fun);
    free(val);
    return 1;
  }
  for(i=0;i<BS;i++)
  {
    for(j=0;j<BS;j++)
    {
      if(vec.QComplex)
      {
        if(M->b0.QComplex) vec.V[i*n + imp*BS + j]    = M->b0.V[i*BS+j];
        else               vec.V[i*n + imp*BS + j].Re = M->b0.VRe[i*BS+j];
      }
      else                 vec.VRe[i*n + imp*BS + j]  = M->b0.VRe[i*BS+j];
    }
  }
  if(CompactMatrixMultiply('*', 'T', &fun, &vec, &funvec))
  {
    printf("CompactMatrixMultiply failed in BlockDoubleTriDiagonalMatrixToBlockListOfPoles\n");
    fflush(stdout);
    FreeCompactMatrix(&fun);
    FreeCompactMatrix(&vec);
    free(val);
    return 1;
  }
  FreeCompactMatrix(&fun);
  FreeCompactMatrix(&vec);

  H->N = n;
  if(InitBlockListOfPoles(H))
  {
    printf("InitBlockListOfPoles failed in BlockDoubleTriDiagonalMatrixToBlockListOfPoles\n");
    fflush(stdout);
    FreeCompactMatrix(&funvec);
    free(val);
    return 1;
  }
  strncpy(H->Name, M->Name, MAXStrLength);
  H->Name[MAXStrLength-1] = '\0';
  H->mu = M->mu;
  stride = 1 + BS2;
  for(k=0;k<n;k++)
  {
    base = k * stride;
    if(H->QComplex)
    {
      H->V[base].Re = val[k];
      H->V[base].Im = 0;
      for(i=0;i<BS;i++)
      {
        for(j=0;j<BS;j++)
        {
          H->V[base+1+i*BS+j].Re =  funvec.V[k*BS+j].Re*funvec.V[k*BS+i].Re + funvec.V[k*BS+j].Im*funvec.V[k*BS+i].Im;
          H->V[base+1+i*BS+j].Im = -funvec.V[k*BS+j].Re*funvec.V[k*BS+i].Im + funvec.V[k*BS+j].Im*funvec.V[k*BS+i].Re;
        }
      }
    }
    else
    {
      H->VRe[base] = val[k];
      for(i=0;i<BS;i++)
      {
        for(j=0;j<BS;j++)
        {
          H->VRe[base+1+i*BS+j] = funvec.VRe[k*BS+j]*funvec.VRe[k*BS+i];
        }
      }
    }
  }
  FreeCompactMatrix(&funvec);
  free(val);
  /*  a0 is the offset of the whole function and is carried along untouched */
  if(H->QComplex)
  {
    if(M->a0.QComplex) memcpy(H->a0, M->a0.V, BS2*sizeof(Complex16));
    else for(i=0;i<BS2;i++) { H->a0[i].Re = M->a0.VRe[i]; H->a0[i].Im = 0; }
  }
  else
  {
    memcpy(H->a0Re, M->a0.VRe, BS2*sizeof(double));
  }
  /*  a negative MinimalPoleDistance asks BlockListOfPolesMergeDegeneratePoles for its
   *  standard value, and that is the one place deciding it. Merging is what restores the
   *  symmetry of a degenerate multiplet: the weights of a degenerate group add to a
   *  spectral projector, which does not depend on the basis the eigensolver picked
   *  inside it, while each weight separately does */
  return BlockListOfPolesMergeDegeneratePoles(H, MinimalPoleDistance, relativeDistance);
}
uint32_t BlockDoubleTriDiagonalMatrixToCompactMatrix(BlockDoubleTriDiagonalMatrixType const M, CompactMatrixType * const H)
{
  //TODO
  printf("TODO: BlockDoubleTriDiagonalMatrixToCompactMatrix\n");
  return 1;
}


uint32_t BlockDoubleTriDiagonalMatrixToSpectraWithBroadening(const double Energy, BlockDoubleTriDiagonalMatrixType const * A, SpectraType *const G, double const Emin, double const Emax, uint32_t const N, double const Gamma)
{
  int iw,i;
  double omega;
  CompactMatrixType Gw;
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
  for(iw=0;iw<=N;iw++)
  {
    omega = Emin + iw * G[0].dE + Energy;
    /*  Evaluate allocates Gw on every call, so it has to be freed inside the loop and not
     *  once at the end.
     */
    EvaluateBlockDoubleTriDiagonalMatrix(A, &Gw, omega, Gamma);
    // add Gw to G, Gw is always complex as Evaluate inverts w + i Gamma/2 - ...
    for(i=0;i<A->BlockSize*A->BlockSize;i++)
    {
      G[i].Re[iw]=Gw.V[i].Re;
      G[i].Im[iw]=Gw.V[i].Im;
    }
    FreeCompactMatrix(&Gw);
  }
  return 0;
}

uint32_t BlockDoubleTriDiagonalMatrixToSpectra(const double Energy, BlockDoubleTriDiagonalMatrixType const * A, SpectraType *const G, double const Emin, double const Emax, uint32_t const N)
{
 BlockListOfPoles ALP;
 /*  A spectrum without broadening is a set of delta functions, which needs the poles. Do not
  *  continue when that conversion is not there: ALP would be left uninitiallized.
  */
 if(BlockDoubleTriDiagonalMatrixToBlockListOfPoles(A, &ALP, -1, false))
 {
   printf("Error in BlockDoubleTriDiagonalMatrixToSpectra: can not create a spectrum without broadening without BlockDoubleTriDiagonalMatrixToBlockListOfPoles\n");
   return 1;
 }
 BlockListOfPolesToSpectra(Energy, &ALP, G, Emin, Emax, N);
 FreeBlockListOfPoles(&ALP);
 return 0;
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
uint32_t BlockNaturalImpurityOrbitalMatrixToSpectraWithBroadening(double GroundStateEnergy, BlockNaturalImpurityOrbitalMatrixType const * A, SpectraType *const G, double const Emin, double const Emax, uint32_t const N, double const Gamma)
{
  SpectraType Gval[A->BlockSize*A->BlockSize], Gcon[A->BlockSize*A->BlockSize];
  uint32_t i,j;

  BlockTriDiagonalMatrixToSpectraWithBroadening(GroundStateEnergy, &(A->val), Gval, Emin, Emax, N, Gamma);
  BlockTriDiagonalMatrixToSpectraWithBroadening(GroundStateEnergy, &(A->con), Gcon, Emin, Emax, N, Gamma);

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

  for(i=0;i<A->BlockSize*A->BlockSize;i++)
  {
    for(j=0;j<=N;j++)
    {
      G[i].Re[j] = Gval[i].Re[j] + Gcon[i].Re[j];
      G[i].Im[j] = Gval[i].Im[j] + Gcon[i].Im[j];
    }
  }

  for(i=0;i<A->BlockSize*A->BlockSize;i++)
  {
    SpectrumFree(&Gval[i]);
    SpectrumFree(&Gcon[i]);
  }
  return 0;
}


uint32_t BlockNaturalImpurityOrbitalMatrixToSpectra(double GroundStateEnergy, BlockNaturalImpurityOrbitalMatrixType const * A, SpectraType *const G, double const Emin, double const Emax, uint32_t const N)
{
 SpectraType Gval[A->BlockSize*A->BlockSize], Gcon[A->BlockSize*A->BlockSize];
 uint32_t i,j;

 BlockTriDiagonalMatrixToSpectra(GroundStateEnergy, &(A->val), Gval, Emin, Emax, N);
 BlockTriDiagonalMatrixToSpectra(GroundStateEnergy, &(A->con), Gcon, Emin, Emax, N);

 G[0].N=N;
 G[0].Emin=Emin;
 G[0].Emax=Emax;
 G[0].EminPole=Emin;
 G[0].EmaxPole=Emax;
 G[0].dE=(Emax-Emin)/((double)N);
 G[0].Gamma=0;
 for(i=0;i<A->BlockSize*A->BlockSize;i++)
 {
   SpectrumAllocate(G[0], &(G[i]));
 }
  if(A->a0.QComplex==0)
  {
    for(i=0;i<A->BlockSize*A->BlockSize;i++)
    {
      for(j=0;j<=N;j++)
      {
        G[i].Re[j] = Gval[i].Re[j] + Gcon[i].Re[j] + A->a0.VRe[i];
        G[i].Im[j] = Gval[i].Im[j] + Gcon[i].Im[j] ;
      }
    }
  }
  else
  {
    for(i=0;i<A->BlockSize*A->BlockSize;i++)
    {
      for(j=0;j<=N;j++)
      {
        G[i].Re[j] = Gval[i].Re[j] + Gcon[i].Re[j] + A->a0.V[i].Re;
        G[i].Im[j] = Gval[i].Im[j] + Gcon[i].Im[j] + A->a0.V[i].Im;
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
  printf("ToDo BlockNaturalImpurityOrbitalMatrixEigenSystem\n");
  //TODO all
  return 0;
}

uint32_t BlockNaturalImpurityOrbitalMatrixToBlockAndersonMatrix(BlockNaturalImpurityOrbitalMatrixType const * A, BlockAndersonMatrixType * const B)
{
  //TODO: algorithm can be optimized
  BlockListOfPoles ALP;
  BlockNaturalImpurityOrbitalMatrixToBlockListOfPoles(A,&ALP, -1, false);
  //TODO: the plus one should not be here
  BlockListOfPolesToBlockAndersonMatrix(ALP, B, ALP.N+1, 1000*DBL_EPSILON);
  FreeBlockListOfPoles(&ALP);
  return 0;
}


uint32_t BlockNaturalImpurityOrbitalMatrixToBlockListOfPoles(BlockNaturalImpurityOrbitalMatrixType const * A, BlockListOfPoles * const B, double const MinimalPoleDistance, bool const relativeDistance)
{
  uint32_t i;
  BlockListOfPoles Bval;
  BlockTriDiagonalMatrixToBlockListOfPoles(&(A->con), B, MinimalPoleDistance, relativeDistance);
  BlockTriDiagonalMatrixToBlockListOfPoles(&(A->val), &Bval, MinimalPoleDistance, relativeDistance);
  AddBlockListOfPoles(B, Bval);
  FreeBlockListOfPoles(&Bval);
  if((B->QComplex==0)&&(A->a0.QComplex==0))
  { // all real
    for(i=0;i<Sqr(B->BlockSize);i++)
    {
      B->a0Re[i]+=A->a0.VRe[i];
    }
  }
  else if(A->a0.QComplex==0)
  { // B Complex but A->a0 real
    for(i=0;i<Sqr(B->BlockSize);i++)
    {
      B->a0[i].Re+=A->a0.VRe[i];
    }
  }
  else if(B->QComplex==0)
  { // B real, but A->a0 complex
    MakeRealBlockListOfPolesComplex(B);
    for(i=0;i<Sqr(B->BlockSize);i++)
    {
      B->a0[i]+=A->a0.V[i];
    }
  }
  else
  { // all complex
    for(i=0;i<Sqr(B->BlockSize);i++)
    {
      B->a0[i]+=A->a0.V[i];
    }
  }
  return 0;
}
uint32_t BlockNaturalImpurityOrbitalMatrixToBlockListOfPolesVariableGrid(BlockNaturalImpurityOrbitalMatrixType const AAnd, BlockListOfPoles * const B, const double * const Grid, const uint32_t NGridPoints)
{
  //TODO: all
  printf("Not implemented BlockNaturalImpurityOrbitalMatrixToBlockListOfPolesVariableGrid\n");
  return 1;
}

/*  A block natural impurity orbital matrix to the block double tri-diagonal form,
 *  through the Anderson form.
 *
 *  These two types are not a change of basis on each other. The natural orbital form
 *  adds two Green's functions, G = A0 + B0^* (gval + gcon) B0^T, while the double
 *  tri-diagonal form nests its two chains in one denominator. The only thing they agree
 *  on is G, so the route is through a form that holds G plainly, and the Anderson form
 *  is the one that then hands over a bath to split.
 *
 *  BlockNaturalImpurityOrbitalMatrixToBlockAndersonMatrix goes by way of the poles of
 *  G, so this composition diagonalizes twice. Its own comment says the algorithm can
 *  be optimized; it reproduces G to 5e-13, which is what matters here, and a direct
 *  route would still have to invert somewhere.
 *
 *  B is not allocated on entry.
 */
uint32_t BlockNaturalImpurityOrbitalMatrixToBlockDoubleTridiagonalMatrix(BlockNaturalImpurityOrbitalMatrixType const * A, BlockDoubleTriDiagonalMatrixType * const B, const double kBT, const double SingularValue, double const MinimalPoleDistance)
{
  BlockAndersonMatrixType AAnd;
  if (BlockNaturalImpurityOrbitalMatrixToBlockAndersonMatrix(A, &AAnd))
  {
    printf("BlockNaturalImpurityOrbitalMatrixToBlockAndersonMatrix failed in BlockNaturalImpurityOrbitalMatrixToBlockDoubleTridiagonalMatrix\n");
    fflush(stdout);
    return 1;
  }
  if (BlockAndersonMatrixToBlockDoubleTriDiagonalMatrix(&AAnd, B, kBT, SingularValue, MinimalPoleDistance))
  {
    printf("BlockAndersonMatrixToBlockDoubleTriDiagonalMatrix failed in BlockNaturalImpurityOrbitalMatrixToBlockDoubleTridiagonalMatrix\n");
    fflush(stdout);
    FreeBlockAndersonMatrix(&AAnd);
    return 1;
  }
  FreeBlockAndersonMatrix(&AAnd);
  return 0;
}


uint32_t BlockNaturalImpurityOrbitalMatrixToBlockTridiagonalMatrix(BlockNaturalImpurityOrbitalMatrixType const * A, BlockTriDiagonalMatrixType * const B, const double SingularValue)
{
  //TODO: algorithm can be optimized use that one needs to tridiagonalize an already banded matrix
  BlockListOfPoles ALP;
  BlockNaturalImpurityOrbitalMatrixToBlockListOfPoles(A,&ALP, -1, false);
  //TODO: the plus one should not be here
  BlockListOfPolesToBlockTridiagonalMatrix(ALP, B, ALP.N+1, SingularValue);
  FreeBlockListOfPoles(&ALP);
  return 0;
}



/*  =========================================================================
 *  Inversed BlockListOfPoles to __
 *  ========================================================================= */


/*  Turn a block list of poles into a matrix
 *       ( A0   B1  B2  B3  ..  Bn  )
 *       ( B1^+ a1  0   0   ..  0   )
 *   B = ( B2^+ 0   a2  0   ..  0   )
 *       ( ..   ..  ..  ..  ..  ..  )
 *       ( Bn^+ 0   0   0   ..  an  )
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
  free(Anow.Name);
  return 0;
}

/*  =========================================================================
 *  Self energies embedded as a matrix
 *
 *  Each representation of a self energy has the form
 *     Sigma(w) = a0 + C . 1/(w - Hbath) . C^+
 *  and is written as a matrix
 *       ( a0    C     )
 *   B = ( C^+   Hbath )
 *  such that (w-B)^-1_[blocksize,blocksize] = 1/(w-Sigma). The first block is the site the self energy
 *  acts on, the other blocks are the bath of the representation, block by block as stored.
 *  ========================================================================= */

/*  Add the n by n block M (row major) to B at block position (brow,bcol). mode 'N' adds M, mode '*' adds
 *  M^*. If AddAdjoint the adjoint of what is added is also added at block position (bcol,brow).
 *  M is real (MRe) or complex (M) according to QComplex, B can only be real if M is real.
 */
static void AddBlockToCompactMatrix(CompactMatrixType * const B, uint32_t const brow, uint32_t const bcol, uint32_t const n, uint32_t const QComplex, double const * const MRe, Complex16 const * const M, char const mode, bool const AddAdjoint)
{
  uint32_t i,j,row,col;
  double re,im;
  for(i=0;i<n;i++)
  {
    for(j=0;j<n;j++)
    {
      row = brow*n+i;
      col = bcol*n+j;
      if(QComplex)
      {
        re = M[i*n+j].Re;
        im = (mode=='*') ? (-M[i*n+j].Im) : (M[i*n+j].Im);
      }
      else
      {
        re = MRe[i*n+j];
        im = 0;
      }
      if(B->QComplex)
      {
        B->V[row*B->Nj+col].Re += re;
        B->V[row*B->Nj+col].Im += im;
        if(AddAdjoint)
        {
          B->V[col*B->Nj+row].Re += re;
          B->V[col*B->Nj+row].Im -= im;
        }
      }
      else
      {
        B->VRe[row*B->Nj+col] += re;
        if(AddAdjoint)
        {
          B->VRe[col*B->Nj+row] += re;
        }
      }
    }
  }
}

/*  Attach the block tridiagonal matrix
 *     T(w) = a0 + b0^* . 1/(w - a1 - b1 . ( ... ) . b1^+) . b0^T
 *  to block host of B: a0 is added to the host block, a1 ... aN are placed at the block positions
 *  first ... first+N-1, b0^* couples the host to a1 and bk couples ak to ak+1.
 */
static void AttachBlockTriDiagonalMatrixToCompactMatrix(CompactMatrixType * const B, BlockTriDiagonalMatrixType const * const T, uint32_t const host, uint32_t const first)
{
  uint32_t k,n=T->BlockSize;
  AddBlockToCompactMatrix(B, host, host, n, T->QComplex, T->aRe, T->a, 'N', false);
  for(k=1;k<=T->N;k++)
  {
    AddBlockToCompactMatrix(B, first+k-1, first+k-1, n, T->QComplex, &(T->aRe[k*n*n]), &(T->a[k*n*n]), 'N', false);
  }
  if(T->N>0)
  {
    AddBlockToCompactMatrix(B, host, first, n, T->QComplex, T->bRe, T->b, '*', true);
  }
  for(k=1;k<T->N;k++)
  {
    AddBlockToCompactMatrix(B, first+k-1, first+k, n, T->QComplex, &(T->bRe[k*n*n]), &(T->b[k*n*n]), 'N', true);
  }
}

static uint32_t InitEmbeddedSelfEnergy(CompactMatrixType * const B, uint32_t const NBlocks, uint32_t const BlockSize, uint32_t const QComplex, double const mu, char const * const caller)
{
  B->QComplex = QComplex;
  B->Ni = NBlocks * BlockSize;
  B->Nj = NBlocks * BlockSize;
  if(InitCompactMatrix(B))
  {
    printf("InitCompactMatrix failed in %s\n", caller);
    return 1;
  }
  B->mu = mu;
  return 0;
}

/*  Turn a block tridiagonal matrix, representing a self energy
 *     Sigma = a0 + b0^* . 1/(w - a1 - b1 . ( ... ) . b1^+) . b0^T
 *  into a matrix
 *       ( a0     b0^*   0      0   ..  )
 *       ( b0^T   a1     b1     0   ..  )
 *   B = ( 0      b1^+   a2     b2  ..  )
 *       ( ..     ..     ..     ..  ..  )
 *       ( 0      0      ..   bN-1^+ aN )
 *
 *   The greensfunction (w-B)^-1_[blocksize,blocksize] is equal to 1/(w-Sigma), i.e. the inverse of the
 *   greensfunction defined by the BlockTriDiagonalMatrix. B has dimension (A->N+1)*A->BlockSize
 *
 */
uint32_t BlockTriDiagonalMatrixInversedToCompactMatrix(BlockTriDiagonalMatrixType const * const A, CompactMatrixType * const B)
{
  if(InitEmbeddedSelfEnergy(B, A->N+1, A->BlockSize, A->QComplex, A->mu, "BlockTriDiagonalMatrixInversedToCompactMatrix")) return 1;
  AttachBlockTriDiagonalMatrixToCompactMatrix(B, A, 0, 1);
  return 0;
}

/*  Turn a block Anderson matrix, representing a self energy
 *     Sigma = A0 + B0^* ( 1 / ( w - A1 - \sum_i ( B_i 1 / ( w - A_{i+1} ) . B_i^+ ) ) ) . B0^T
 *  into a matrix
 *       ( A0     B0^*   0      0   ..  )
 *       ( B0^T   A1     B1     B2  ..  )
 *   B = ( 0      B1^+   A2     0   ..  )
 *       ( 0      B2^+   0      A3  ..  )
 *       ( ..     ..     ..     ..  ..  )
 *
 *   such that (w-B)^-1_[blocksize,blocksize] is equal to 1/(w-Sigma). B has dimension (A->N+1)*A->BlockSize
 */
uint32_t BlockAndersonMatrixInversedToCompactMatrix(BlockAndersonMatrixType const * const A, CompactMatrixType * const B)
{
  uint32_t k,n=A->BlockSize;
  if(InitEmbeddedSelfEnergy(B, A->N+1, A->BlockSize, A->QComplex, A->mu, "BlockAndersonMatrixInversedToCompactMatrix")) return 1;
  for(k=0;k<=A->N;k++)
  {
    AddBlockToCompactMatrix(B, k, k, n, A->QComplex, &(A->aRe[k*n*n]), &(A->a[k*n*n]), 'N', false);
  }
  if(A->N>0)
  {
    AddBlockToCompactMatrix(B, 0, 1, n, A->QComplex, A->bRe, A->b, '*', true);
  }
  for(k=1;k<A->N;k++)
  {
    AddBlockToCompactMatrix(B, 1, k+1, n, A->QComplex, &(A->bRe[k*n*n]), &(A->b[k*n*n]), 'N', true);
  }
  return 0;
}

/*  Turn a block natural impurity orbital matrix, representing a self energy
 *     Sigma = a0 + Gval(w) + Gcon(w)
 *  with Gval and Gcon block tridiagonal matrices, into a matrix
 *       ( a0+av0+ac0  bv0^*  ..  bc0^*  ..  )
 *       ( bv0^T       av1    bv1  0     ..  )
 *   B = ( ..          bv1^+  ..   0     ..  )
 *       ( bc0^T       0      0    ac1   bc1 )
 *       ( ..          0      0    bc1^+ ..  )
 *
 *   such that (w-B)^-1_[blocksize,blocksize] is equal to 1/(w-Sigma).
 *   B has dimension (1+A->val.N+A->con.N)*A->BlockSize. TMT is not needed.
 */
uint32_t BlockNaturalImpurityOrbitalMatrixInversedToCompactMatrix(BlockNaturalImpurityOrbitalMatrixType const * const A, CompactMatrixType * const B)
{
  uint32_t n=A->BlockSize;
  if((A->val.BlockSize!=n)||(A->con.BlockSize!=n)||(A->a0.Ni!=n)||(A->a0.Nj!=n))
  {
    printf("Block sizes of a0, val and con are not equal to BlockSize in BlockNaturalImpurityOrbitalMatrixInversedToCompactMatrix\n");
    return 1;
  }
  if(InitEmbeddedSelfEnergy(B, 1+A->val.N+A->con.N, n, A->QComplex|A->a0.QComplex|A->val.QComplex|A->con.QComplex, A->mu, "BlockNaturalImpurityOrbitalMatrixInversedToCompactMatrix")) return 1;
  AddBlockToCompactMatrix(B, 0, 0, n, A->a0.QComplex, A->a0.VRe, A->a0.V, 'N', false);
  AttachBlockTriDiagonalMatrixToCompactMatrix(B, &(A->val), 0, 1);
  AttachBlockTriDiagonalMatrixToCompactMatrix(B, &(A->con), 0, 1+A->val.N);
  return 0;
}

/*  Turn a block double tridiagonal matrix, representing a self energy
 *     Sigma = a0 + b0^* . 1/(w - Gval(w) - Gcon(w)) . b0^T
 *  with Gval and Gcon block tridiagonal matrices, into a matrix
 *       ( a0     b0^*     0      ..  0      ..  )
 *       ( b0^T   av0+ac0  bv0^*  ..  bc0^*  ..  )
 *   B = ( 0      bv0^T    av1    ..  0      ..  )
 *       ( ..     ..       ..     ..  0      ..  )
 *       ( 0      bc0^T    0      0   ac1    ..  )
 *       ( ..     ..       0      0   ..     ..  )
 *
 *   such that (w-B)^-1_[blocksize,blocksize] is equal to 1/(w-Sigma).
 *   B has dimension (2+A->val.N+A->con.N)*A->BlockSize.
 */
uint32_t BlockDoubleTriDiagonalMatrixInversedToCompactMatrix(BlockDoubleTriDiagonalMatrixType const * const A, CompactMatrixType * const B)
{
  uint32_t n=A->BlockSize;
  if((A->val.BlockSize!=n)||(A->con.BlockSize!=n)||(A->a0.Ni!=n)||(A->a0.Nj!=n)||(A->b0.Ni!=n)||(A->b0.Nj!=n))
  {
    printf("Block sizes of a0, b0, val and con are not equal to BlockSize in BlockDoubleTriDiagonalMatrixInversedToCompactMatrix\n");
    return 1;
  }
  if(InitEmbeddedSelfEnergy(B, 2+A->val.N+A->con.N, n, A->QComplex|A->a0.QComplex|A->b0.QComplex|A->val.QComplex|A->con.QComplex, A->mu, "BlockDoubleTriDiagonalMatrixInversedToCompactMatrix")) return 1;
  AddBlockToCompactMatrix(B, 0, 0, n, A->a0.QComplex, A->a0.VRe, A->a0.V, 'N', false);
  AddBlockToCompactMatrix(B, 0, 1, n, A->b0.QComplex, A->b0.VRe, A->b0.V, '*', true);
  AttachBlockTriDiagonalMatrixToCompactMatrix(B, &(A->val), 1, 2);
  AttachBlockTriDiagonalMatrixToCompactMatrix(B, &(A->con), 1, 2+A->val.N);
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
        FreeCompactMatrix(&T);
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
            ComplexEigenFunctionToResidue(&(T.V[ival*A.NFermion]), &H);
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
        FreeMatrix(&H);
        FreeCompactMatrix(&T);
      }
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

/* Given a strictly increasing array E[0..n-1], return the index i of the interval
 * that contains val, i.e.  E[i] <= val < E[i+1].
 *   returns -1   if val <  E[0]     (below the grid)
 *   returns n-1  if val >= E[n-1]   (at or above the last grid point)
 * The second convention is taken from TightBindingHamiltonianToBlockListOfPoles(), where val == Emax is also
 * counted as "outside" and pushed into the pole above the grid; see also the notes for TightBindingHamiltonianToBlockListOfPolesGrid().
 */
static int FindEnergyInterval(double const * const E, int const n, double const val)
{
  int lo, hi, mid;
  if (val <  E[0])   return -1;
  if (val >= E[n-1]) return n-1;
  lo = 0;
  hi = n-1;
  while (hi - lo > 1)
  {
    mid = (lo + hi) / 2;
    if (val < E[mid]) hi = mid;
    else              lo = mid;
  }
  return lo;
}

/*  Given a tight binding Hamiltonian, a k-mesh defined by the number of divisions of
 *  a, b and c (Nk), and an arbitrary strictly increasing energy grid Egrid[0..NEgrid-1],
 *  returns the Green's function matrix of size A.NFermion by A.NFermion as a list of
 *  NEgrid+2 poles: the NEgrid grid points, plus one pole below and one above the grid
 *  that collect everything outside, positioned such that weight and first moment are exact.
 *  For each k and eigenvalue i the weight h_{jk} = fun[i,j] * fun[i,k]^* is distributed
 *  linearly over the two poles bracketing the eigenvalue.
 */
uint32_t TightBindingHamiltonianToBlockListOfPolesGrid(TightBindingDef const A, int Nk[3], BlockListOfPoles * const G, double const * const Egrid, uint32_t const NEgrid)
{
  MatrixType H;
  CompactMatrixType T;
  double q[3], val[A.NFermion], rest, totmin=0, totmax=0, totnow;
  int    iq[3],ival,iE,i,j,ii;
  double diq[3],Nktot;
  if(NEgrid<2)
  {
    printf("Error in TightBindingHamiltonianToBlockListOfPolesGrid: the energy grid needs at least 2 points\n");
    return 1;
  }
  G->BlockSize=A.NFermion;
  G->N=NEgrid+2;
  G->QComplex=1;
  strcpy(G->Name, A.Name);
  if(InitBlockListOfPoles(G))
  {
    printf("InitBlockListOfPoles failed in TightBindingHamiltonianToBlockListOfPolesGrid\n");
    return 1;
  }
  ii=-1;
  for(i=0;i<G->N*(1+G->BlockSize*G->BlockSize);i+=(1+G->BlockSize*G->BlockSize))
  {
    if(ii<0)                  G->V[i].Re=Egrid[0] - (Egrid[1] - Egrid[0]);
    else if(ii<(int)NEgrid)   G->V[i].Re=Egrid[ii];
    else                      G->V[i].Re=Egrid[NEgrid-1] + (Egrid[NEgrid-1] - Egrid[NEgrid-2]);
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
            ComplexEigenFunctionToResidue(&(T.V[ival*A.NFermion]), &H);
          }
          else
          {
            RealEigenFunctionToCharacter(&(T.VRe[ival*A.NFermion]), &H);
          }
          /* determine the position of the pole in the list of poles */
          iE = FindEnergyInterval(Egrid, (int)NEgrid, val[ival]);
          if (iE >= (int)NEgrid-1)
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
            G->V[(NEgrid+1)*(1+G->BlockSize*G->BlockSize)].Re = (totmax * G->V[(NEgrid+1)*(1+G->BlockSize*G->BlockSize)].Re + totnow * val[ival]) / (totmax + totnow);
            for(i=0;i<A.NFermion;i++)
            {
              for(j=0;j<A.NFermion;j++)
              {
                G->V[(NEgrid+1)*(1+G->BlockSize*G->BlockSize)+1+j+i*A.NFermion].Re+=H.VRe[i][j];
                if(H.QComplex)
                {
                  G->V[(NEgrid+1)*(1+G->BlockSize*G->BlockSize)+1+j+i*A.NFermion].Im+=H.VIm[i][j];
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
            rest = (val[ival] - Egrid[iE]) / (Egrid[iE+1] - Egrid[iE]);  /* rest is a number between 0 and 1 pole at iE+1 gets 1-rest pole at iE+2 gest rest */
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
        FreeMatrix(&H);
        FreeCompactMatrix(&T);
      }
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
        FreeMatrix(&H);
        FreeCompactMatrix(&T);
      }
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
            ComplexEigenFunctionToResidue(&(T.V[ival*T.Nj]), &H);
          }
          else
          {
            RealEigenFunctionToCharacter(&(T.VRe[ival*T.Nj]), &H);
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
        FreeMatrix(&H);
        FreeCompactMatrix(&T);
      }
    }
  }
  FreeCompactMatrix(&HSigma);
  if(Sigma.QComplex==0) FreeCompactMatrix(&HSigmaComplex);
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


/*  Given a tight binding Hamiltonian, a k-mesh defined by the number of divisions of
 *  a, b and c (Nk), and an arbitrary strictly increasing energy grid Egrid[0..NEgrid-1],
 *  returns the Green's function matrix of size A.NFermion by A.NFermion as a list of
 *  NEgrid+2 poles: the NEgrid grid points, plus one pole below and one above the grid
 *  that collect everything outside, positioned such that weight and first moment are exact.
 *  For each k and eigenvalue i the weight h_{jk} = fun[i,j] * fun[i,k]^* is distributed
 *  linearly over the two poles bracketing the eigenvalue.
 *  Sigma must have dimensions BlockSize equal to A.NFermion
 */
uint32_t TightBindingHamiltonianPlusSelfEnergyToBlockListOfPolesGrid(TightBindingDef const A, int Nk[3], BlockListOfPoles const Sigma, BlockListOfPoles * const G, double const * const Egrid, uint32_t const NEgrid)
{
  MatrixType H;
  CompactMatrixType T, HSigma, HSigmaComplex;
  double q[3], val[A.NFermion*(Sigma.N+1)], rest, totmin=0, totmax=0, totnow;
  int    iq[3],ival,iE,i,j,ii;
  double diq[3],Nktot;
  // test if sigma size is ok
  if(Sigma.BlockSize!=A.NFermion)
  {
    printf("in TightBindingHamiltonianPlusSelfEnergyToBlockListOfPoles\nThe dimension of the self energy (%i) is not equal to the dimension of the tight binding Hamiltonian (%i)\n",Sigma.BlockSize,A.NFermion);
    return 1;
  }
  // init G
  G->BlockSize=A.NFermion;
  G->N=NEgrid+2;
  G->QComplex=1;
  strcpy(G->Name, A.Name);
  if(InitBlockListOfPoles(G))
  {
    printf("InitBlockListOfPoles failed in TightBindingHamiltonianToBlockListOfPoles\n");
    return 1;
  }
  // set energy of G
  ii=-1;
  for(i=0;i<G->N*(1+G->BlockSize*G->BlockSize);i+=(1+G->BlockSize*G->BlockSize))
  {
    if(ii<0)                  G->V[i].Re=Egrid[0] - (Egrid[1] - Egrid[0]);
    else if(ii<(int)NEgrid)   G->V[i].Re=Egrid[ii];
    else                      G->V[i].Re=Egrid[NEgrid-1] + (Egrid[NEgrid-1] - Egrid[NEgrid-2]);
    ii++;
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
            ComplexEigenFunctionToResidue(&(T.V[ival*T.Nj]), &H);
          }
          else
          {
            RealEigenFunctionToCharacter(&(T.VRe[ival*T.Nj]), &H);
          }
          /* determine the position of the pole in the list of poles */
          iE = FindEnergyInterval(Egrid, (int)NEgrid, val[ival]);
          if (iE >= (int)NEgrid-1)
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
            G->V[(NEgrid+1)*(1+G->BlockSize*G->BlockSize)].Re = (totmax * G->V[(NEgrid+1)*(1+G->BlockSize*G->BlockSize)].Re + totnow * val[ival]) / (totmax + totnow);
            for(i=0;i<A.NFermion;i++)
            {
              for(j=0;j<A.NFermion;j++)
              {
                G->V[(NEgrid+1)*(1+G->BlockSize*G->BlockSize)+1+j+i*A.NFermion].Re+=H.VRe[i][j];
                if(H.QComplex)
                {
                  G->V[(NEgrid+1)*(1+G->BlockSize*G->BlockSize)+1+j+i*A.NFermion].Im+=H.VIm[i][j];
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
            rest = (val[ival] - Egrid[iE]) / (Egrid[iE+1] - Egrid[iE]);  /* rest is a number between 0 and 1 pole at iE+1 gets 1-rest pole at iE+2 gest rest */
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
        FreeMatrix(&H);
        FreeCompactMatrix(&T);
      }
    }
  }
  FreeCompactMatrix(&HSigma);
  if(Sigma.QComplex==0) FreeCompactMatrix(&HSigmaComplex);
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


/*  Per atom augmented propagators, accumulated in the same k-loop as the local Green
 *  function by TightBindingHamiltonianPlusEmbeddedSelfEnergyToBlockListOfPolesGrid.
 *
 *  The augmented propagator of atom a is the 2 n_a by 2 n_a block response function of the
 *  operator pair ( c_a , H c_a ), where c_a are the orbitals of atom a and H is the tight
 *  binding Hamiltonian. For a one body H the second operator is the commutator [c_a, H], so
 *  this is the same construction as the augmented propagator (a, [a,H1]) that gives the self
 *  energy, with the one body Hamiltonian in place of the interaction:
 *
 *      Gaug = ( <c|Q|c>      <c|Q H|c>   )          Q = ( w - E - H(k) )^-1
 *             ( <c|H Q|c>    <c|H Q H|c> )
 *
 *  Its Schur complement onto the second block,
 *
 *      S = <c|HQH|c> - <c|HQ|c> <c|Q|c>^-1 <c|QH|c> = zeta_a - (G_aa)^-1 - <H>_aa ,
 *
 *  is the hybridisation function of atom a, with zeta = w - Sigma(w), so
 *  ResponseFunction.Downfold on the second block returns it. The identity needs Sigma to be
 *  block diagonal in the atom index -- the site local self energy of single site DMFT -- and
 *  it then survives the projection onto one atom, because the projector commutes with zeta.
 *
 *  What this buys is conditioning. zeta appears only inside the resolvent and is never formed
 *  on its own, so nothing large is ever subtracted. Taking the hybridisation as zeta - G^-1
 *  instead loses all significance where zeta is large, which is what happens near the pole a
 *  Mott insulating self energy has at the chemical potential.
 *
 *  Every eigenvector contributes the outer product of v = ( u_a , (H u)_a ), where u_a are
 *  the components of its first A.NFermion entries that belong to atom a, distributed over the
 *  two poles bracketing the eigenvalue exactly as the weight of G is.
 */
typedef struct
{
  int NAtoms;
  int *off;                  /* first orbital of each atom */
  int *nor;                  /* number of orbitals of each atom */
  double *totmin;            /* weight collected in the pole below the grid, per atom */
  double *totmax;            /* and in the pole above it */
  BlockListOfPoles *Gaug;
  Complex16 *hu;             /* work space of A.NFermion */
  Complex16 *v;              /* work space of max 2 n_a */
}
AugmentedPerAtomType;

static uint32_t AugmentedPerAtomInit(AugmentedPerAtomType * const S, TightBindingDef const A, BlockListOfPoles * const Gaug, double const * const Egrid, uint32_t const NEgrid)
{
  int a,i,j,ii,nb=0;
  S->NAtoms = (int)A.Cell.NAtoms;
  S->Gaug   = Gaug;
  S->off    = (int*)malloc(S->NAtoms*sizeof(int));
  S->nor    = (int*)malloc(S->NAtoms*sizeof(int));
  S->totmin = (double*)calloc(S->NAtoms,sizeof(double));
  S->totmax = (double*)calloc(S->NAtoms,sizeof(double));
  S->hu     = (Complex16*)malloc(A.NFermion*sizeof(Complex16));
  S->v      = NULL;
  if(S->off==NULL||S->nor==NULL||S->totmin==NULL||S->totmax==NULL||S->hu==NULL) return 1;
  for(a=0;a<S->NAtoms;a++)
  {
    S->off[a] = A.FermionIndex[a][0];
    S->nor[a] = 0;
    for(j=0;j<A.Cell.Atom[a].Nshell;j++) S->nor[a] += A.Cell.Atom[a].shelldegeneracy[j];
    if(2*S->nor[a]>nb) nb = 2*S->nor[a];
  }
  S->v = (Complex16*)malloc((nb>0?nb:1)*sizeof(Complex16));
  if(S->v==NULL) return 1;
  for(a=0;a<S->NAtoms;a++)
  {
    Gaug[a].BlockSize = 2*S->nor[a];
    Gaug[a].N         = NEgrid+2;
    Gaug[a].QComplex  = 1;
    Gaug[a].mu        = 0;
    strncpy(Gaug[a].Name, A.Cell.Atom[a].Name, MAXStrLength-1);
    Gaug[a].Name[MAXStrLength-1]='\0';
    if(InitBlockListOfPoles(&(Gaug[a]))) return 1;
    ii=-1;
    for(i=0;i<(int)(Gaug[a].N*(1+Gaug[a].BlockSize*Gaug[a].BlockSize));i+=(1+Gaug[a].BlockSize*Gaug[a].BlockSize))
    {
      if(ii<0)                Gaug[a].V[i].Re = Egrid[0] - (Egrid[1] - Egrid[0]);
      else if(ii<(int)NEgrid) Gaug[a].V[i].Re = Egrid[ii];
      else                    Gaug[a].V[i].Re = Egrid[NEgrid-1] + (Egrid[NEgrid-1] - Egrid[NEgrid-2]);
      ii++;
    }
  }
  return 0;
}

/*  Must be called before ComplexEigenFunctionToResidue overwrites H with the residue. */
static void AugmentedPerAtomAdd(AugmentedPerAtomType * const S, MatrixType const * const H, CompactMatrixType const * const T, int const ival, int const NFermion, double const val, int const iE, double const * const Egrid, uint32_t const NEgrid)
{
  int a,i,j,p,r,bs,ib;
  double rest, totnow, wRe, wIm, uRe, uIm;
  for(i=0;i<NFermion;i++)
  {
    S->hu[i].Re=0; S->hu[i].Im=0;
    for(j=0;j<NFermion;j++)
    {
      if(T->QComplex) { uRe=T->V[ival*T->Nj+j].Re; uIm=T->V[ival*T->Nj+j].Im; }
      else            { uRe=T->VRe[ival*T->Nj+j];  uIm=0; }
      S->hu[i].Re += H->VRe[i][j]*uRe - ((H->QComplex)?(H->VIm[i][j]*uIm):(0));
      S->hu[i].Im += H->VRe[i][j]*uIm + ((H->QComplex)?(H->VIm[i][j]*uRe):(0));
    }
  }
  for(a=0;a<S->NAtoms;a++)
  {
    bs = 2*S->nor[a];
    if(bs==0) continue;
    for(p=0;p<S->nor[a];p++)
    {
      if(T->QComplex) { S->v[p] = T->V[ival*T->Nj + S->off[a]+p]; }
      else            { S->v[p].Re = T->VRe[ival*T->Nj + S->off[a]+p]; S->v[p].Im = 0; }
      S->v[S->nor[a]+p] = S->hu[S->off[a]+p];
    }
    totnow=0;
    for(p=0;p<bs;p++) totnow += S->v[p].Re*S->v[p].Re + S->v[p].Im*S->v[p].Im;
    if (iE >= (int)NEgrid-1)
    {
      ib = (int)(NEgrid+1);
      if(S->totmax[a]+totnow<10*DBL_MIN) { S->totmax[a]+=10*DBL_MIN; totnow+=10*DBL_MIN; }
      S->Gaug[a].V[ib*(1+bs*bs)].Re = (S->totmax[a]*S->Gaug[a].V[ib*(1+bs*bs)].Re + totnow*val)/(S->totmax[a]+totnow);
      S->totmax[a]+=totnow;
      rest=0;
    }
    else if (iE < 0)
    {
      ib = 0;
      if(S->totmin[a]+totnow<10*DBL_MIN) { S->totmin[a]+=10*DBL_MIN; totnow+=10*DBL_MIN; }
      S->Gaug[a].V[0].Re = (S->totmin[a]*S->Gaug[a].V[0].Re + totnow*val)/(S->totmin[a]+totnow);
      S->totmin[a]+=totnow;
      rest=0;
    }
    else
    {
      ib   = iE+1;
      rest = (val - Egrid[iE]) / (Egrid[iE+1] - Egrid[iE]);
    }
    /*  residue v v^dag, the same convention as ComplexEigenFunctionToResidue */
    for(p=0;p<bs;p++)
    {
      for(r=0;r<bs;r++)
      {
        wRe = S->v[p].Re*S->v[r].Re + S->v[p].Im*S->v[r].Im;
        wIm = S->v[p].Im*S->v[r].Re - S->v[p].Re*S->v[r].Im;
        S->Gaug[a].V[ib*(1+bs*bs)+1+r+p*bs].Re += (1-rest)*wRe;
        S->Gaug[a].V[ib*(1+bs*bs)+1+r+p*bs].Im += (1-rest)*wIm;
        if(rest!=0)
        {
          S->Gaug[a].V[(ib+1)*(1+bs*bs)+1+r+p*bs].Re += rest*wRe;
          S->Gaug[a].V[(ib+1)*(1+bs*bs)+1+r+p*bs].Im += rest*wIm;
        }
      }
    }
  }
}

static void AugmentedPerAtomFinish(AugmentedPerAtomType * const S, double const Nktot)
{
  int a,i,j,bs;
  for(a=0;a<S->NAtoms;a++)
  {
    bs = 2*S->nor[a];
    if(bs==0) continue;
    for(i=0;i<(int)(S->Gaug[a].N*(1+bs*bs));i+=(1+bs*bs))
    {
      for(j=1;j<1+bs*bs;j++)
      {
        S->Gaug[a].V[i+j].Re/=Nktot;
        S->Gaug[a].V[i+j].Im/=Nktot;
      }
    }
    ChopBlockListOfPoles(&(S->Gaug[a]), DBL_EPSILON);
  }
}

static void AugmentedPerAtomFree(AugmentedPerAtomType * const S)
{
  free(S->off); free(S->nor); free(S->totmin); free(S->totmax); free(S->hu); free(S->v);
}

/*  Given a tight binding Hamiltonian, a k-mesh defined by the number of divisions of a, b and c (Nk),
 *  a self energy embedded as a matrix
 *       ( a0    C     )
 *   E = ( C^+   Hbath )
 *  (see BlockTriDiagonalMatrixInversedToCompactMatrix and its siblings) and an arbitrary strictly
 *  increasing energy grid Egrid[0..NEgrid-1], returns the Green's function as a list of NEgrid+2 poles,
 *  as TightBindingHamiltonianPlusSelfEnergyToBlockListOfPolesGrid.
 *  For each k the matrix E with H(k) added to its first block is diagonalised (dense). The first
 *  A.NFermion components of each eigenvector give the weight, which is distributed linearly over the
 *  two poles bracketing the eigenvalue.
 *
 *  The first block of E must have dimension A.NFermion
 */
uint32_t TightBindingHamiltonianPlusEmbeddedSelfEnergyToBlockListOfPolesGrid(TightBindingDef const A, int Nk[3], CompactMatrixType const * const E, BlockListOfPoles * const G, double const * const Egrid, uint32_t const NEgrid, BlockListOfPoles * const Gaug)
{
  MatrixType H, Hk;
  AugmentedPerAtomType Aug;
  CompactMatrixType T, HSigma, HSigmaComplex;
  double q[3], *val, rest, totmin=0, totmax=0, totnow;
  int    iq[3],ival,iE,i,j,ii;
  double diq[3],Nktot;
  // init G
  G->BlockSize=A.NFermion;
  G->N=NEgrid+2;
  G->QComplex=1;
  strcpy(G->Name, A.Name);
  if(InitBlockListOfPoles(G))
  {
    printf("InitBlockListOfPoles failed in TightBindingHamiltonianPlusEmbeddedSelfEnergyToBlockListOfPolesGrid\n");
    return 1;
  }
  // set energy of G
  ii=-1;
  for(i=0;i<G->N*(1+G->BlockSize*G->BlockSize);i+=(1+G->BlockSize*G->BlockSize))
  {
    if(ii<0)                  G->V[i].Re=Egrid[0] - (Egrid[1] - Egrid[0]);
    else if(ii<(int)NEgrid)   G->V[i].Re=Egrid[ii];
    else                      G->V[i].Re=Egrid[NEgrid-1] + (Egrid[NEgrid-1] - Egrid[NEgrid-2]);
    ii++;
  }
  /*  one augmented propagator per atom, filled from the same eigenvectors */
  if(Gaug!=NULL)
  {
    if(AugmentedPerAtomInit(&Aug, A, Gaug, Egrid, NEgrid))
    {
      printf("AugmentedPerAtomInit failed in TightBindingHamiltonianPlusEmbeddedSelfEnergyToBlockListOfPolesGrid\n");
      return 1;
    }
  }
  // create HSigma, the first block (a0) gets H(k) added for each k
  CompactMatrixCopy(*E, &HSigma);
  if(E->QComplex==0)
  {
    CompactMatrixCopy(*E, &HSigmaComplex);
    MakeCompactRealMatrixComplex(&HSigmaComplex);
  }
  val = (double*)malloc(HSigma.Nj*sizeof(double));
  if(val==NULL)
  {
    printf("malloc failed in TightBindingHamiltonianPlusEmbeddedSelfEnergyToBlockListOfPolesGrid\n");
    return 1;
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
        if(E->QComplex)
        {
          for(i=0;i<A.NFermion;i++)
          {
            for(j=0;j<A.NFermion;j++)
            {
              HSigma.V[i*HSigma.Nj+j].Re=E->V[i*E->Nj+j].Re + H.VRe[i][j];
              HSigma.V[i*HSigma.Nj+j].Im=E->V[i*E->Nj+j].Im + ((H.QComplex) ? (H.VIm[i][j]) : (0));
            }
          }
          CompactMatrixEigensystem(&HSigma, &T, val);
        }
        else if(H.QComplex)
        {
          for(i=0;i<A.NFermion;i++)
          {
            for(j=0;j<A.NFermion;j++)
            {
              HSigmaComplex.V[i*HSigma.Nj+j].Re=E->VRe[i*E->Nj+j] + H.VRe[i][j];
              HSigmaComplex.V[i*HSigma.Nj+j].Im=                    H.VIm[i][j];
            }
          }
          CompactMatrixEigensystem(&HSigmaComplex, &T, val);
        }
        else
        {
          for(i=0;i<A.NFermion;i++)
          {
            for(j=0;j<A.NFermion;j++)
            {
              HSigma.VRe[i*HSigma.Nj+j]=E->VRe[i*E->Nj+j] + H.VRe[i][j];
            }
          }
          CompactMatrixEigensystem(&HSigma, &T, val);
        }
        if(T.QComplex)
        {
          if(H.QComplex==0) MakeRealMatrixComplex(&H);
        }
        else
        {
          if(H.QComplex) MakeComplexMatrixReal(&H);
        }
        /*  ComplexEigenFunctionToResidue writes the residue into H, so H holds the
         *  Hamiltonian only until the first eigenvector. Keep a copy of it for H u. */
        if(Gaug!=NULL) MatrixCopy(&Hk, H);
        for(ival=0;ival<T.Nj;ival++)
        {
          /* determine the position of the pole in the list of poles */
          iE = FindEnergyInterval(Egrid, (int)NEgrid, val[ival]);
          if(Gaug!=NULL) AugmentedPerAtomAdd(&Aug, &Hk, &T, ival, A.NFermion, val[ival], iE, Egrid, NEgrid);
          /* calculate matrix B for eigen-value ival */
          if(H.QComplex)
          {
            ComplexEigenFunctionToResidue(&(T.V[ival*T.Nj]), &H);
          }
          else
          {
            RealEigenFunctionToCharacter(&(T.VRe[ival*T.Nj]), &H);
          }
          if (iE >= (int)NEgrid-1)
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
            G->V[(NEgrid+1)*(1+G->BlockSize*G->BlockSize)].Re = (totmax * G->V[(NEgrid+1)*(1+G->BlockSize*G->BlockSize)].Re + totnow * val[ival]) / (totmax + totnow);
            for(i=0;i<A.NFermion;i++)
            {
              for(j=0;j<A.NFermion;j++)
              {
                G->V[(NEgrid+1)*(1+G->BlockSize*G->BlockSize)+1+j+i*A.NFermion].Re+=H.VRe[i][j];
                if(H.QComplex)
                {
                  G->V[(NEgrid+1)*(1+G->BlockSize*G->BlockSize)+1+j+i*A.NFermion].Im+=H.VIm[i][j];
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
            rest = (val[ival] - Egrid[iE]) / (Egrid[iE+1] - Egrid[iE]);  /* rest is a number between 0 and 1 pole at iE+1 gets 1-rest pole at iE+2 gest rest */
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
        if(Gaug!=NULL) FreeMatrix(&Hk);
        FreeMatrix(&H);
        FreeCompactMatrix(&T);
      }
    }
  }
  free(val);
  FreeCompactMatrix(&HSigma);
  if(E->QComplex==0) FreeCompactMatrix(&HSigmaComplex);
  Nktot=Nk[0]*Nk[1]*Nk[2];
  for(i=0;i<G->N*(1+G->BlockSize*G->BlockSize);i+=(1+G->BlockSize*G->BlockSize))
  {
    for(j=1;j<1+G->BlockSize*G->BlockSize;j++)
    {
      G->V[i+j].Re/=Nktot;
      G->V[i+j].Im/=Nktot;
    }
  }
  if(Gaug!=NULL)
  {
    AugmentedPerAtomFinish(&Aug, Nktot);
    AugmentedPerAtomFree(&Aug);
  }
  ChopBlockListOfPoles(G, DBL_EPSILON);
  return 0;
}

/*  As above, on the energy mesh Emin to Emax in NE steps, which gives NE+3 poles as
 *  TightBindingHamiltonianPlusSelfEnergyToBlockListOfPoles
 */
static uint32_t TightBindingHamiltonianPlusEmbeddedSelfEnergyToBlockListOfPoles(TightBindingDef const A, int Nk[3], CompactMatrixType const * const E, BlockListOfPoles * const G, double const Emin, double const Emax, uint32_t const NE)
{
  uint32_t i, ret;
  double dE=(Emax-Emin)/((double)NE);
  double *Egrid = (double*)malloc((NE+1)*sizeof(double));
  if(Egrid==NULL)
  {
    printf("malloc failed in TightBindingHamiltonianPlusEmbeddedSelfEnergyToBlockListOfPoles\n");
    return 1;
  }
  for(i=0;i<=NE;i++)
  {
    Egrid[i]=Emin+i*dE;
  }
  ret = TightBindingHamiltonianPlusEmbeddedSelfEnergyToBlockListOfPolesGrid(A, Nk, E, G, Egrid, NE+1, NULL);
  free(Egrid);
  return ret;
}


/*  As TightBindingHamiltonianPlusSelfEnergyToBlockListOfPoles(Grid) above, but with the self energy given
 *  as a block tridiagonal matrix, a block Anderson matrix, a block natural impurity orbital matrix or a
 *  block double tridiagonal matrix. The self energy is embedded as a matrix of dimension
 *     (1 + number of bath blocks) * A.NFermion
 *  (2 + ... for the double tridiagonal matrix), which is diagonalised for each k. A list of poles with N
 *  poles gives dimension (N+1)*A.NFermion, and the exact conversion of a tridiagonal matrix with N blocks
 *  gives up to N*A.NFermion poles, each with a rank one weight.
 *
 *  Sigma must have dimensions BlockSize equal to A.NFermion
 */
static uint32_t TestSelfEnergyBlockSize(uint32_t const BlockSize, uint32_t const NFermion, char const * const caller)
{
  if(BlockSize!=NFermion)
  {
    printf("in %s\nThe dimension of the self energy (%i) is not equal to the dimension of the tight binding Hamiltonian (%i)\n",caller,BlockSize,NFermion);
    return 1;
  }
  return 0;
}

uint32_t TightBindingHamiltonianPlusSelfEnergyToBlockListOfPoles(TightBindingDef const A, int Nk[3], BlockTriDiagonalMatrixType const Sigma, BlockListOfPoles * const G, double const Emin, double const Emax, uint32_t const NE)
{
  uint32_t ret;
  CompactMatrixType E;
  if(TestSelfEnergyBlockSize(Sigma.BlockSize, A.NFermion, "TightBindingHamiltonianPlusSelfEnergyToBlockListOfPoles")) return 1;
  if(BlockTriDiagonalMatrixInversedToCompactMatrix(&Sigma, &E)) return 1;
  ret = TightBindingHamiltonianPlusEmbeddedSelfEnergyToBlockListOfPoles(A, Nk, &E, G, Emin, Emax, NE);
  FreeCompactMatrix(&E);
  return ret;
}

uint32_t TightBindingHamiltonianPlusSelfEnergyToBlockListOfPolesGrid(TightBindingDef const A, int Nk[3], BlockTriDiagonalMatrixType const Sigma, BlockListOfPoles * const G, double const * const Egrid, uint32_t const NEgrid)
{
  uint32_t ret;
  CompactMatrixType E;
  if(TestSelfEnergyBlockSize(Sigma.BlockSize, A.NFermion, "TightBindingHamiltonianPlusSelfEnergyToBlockListOfPolesGrid")) return 1;
  if(BlockTriDiagonalMatrixInversedToCompactMatrix(&Sigma, &E)) return 1;
  ret = TightBindingHamiltonianPlusEmbeddedSelfEnergyToBlockListOfPolesGrid(A, Nk, &E, G, Egrid, NEgrid, NULL);
  FreeCompactMatrix(&E);
  return ret;
}

uint32_t TightBindingHamiltonianPlusSelfEnergyToBlockListOfPoles(TightBindingDef const A, int Nk[3], BlockAndersonMatrixType const Sigma, BlockListOfPoles * const G, double const Emin, double const Emax, uint32_t const NE)
{
  uint32_t ret;
  CompactMatrixType E;
  if(TestSelfEnergyBlockSize(Sigma.BlockSize, A.NFermion, "TightBindingHamiltonianPlusSelfEnergyToBlockListOfPoles")) return 1;
  if(BlockAndersonMatrixInversedToCompactMatrix(&Sigma, &E)) return 1;
  ret = TightBindingHamiltonianPlusEmbeddedSelfEnergyToBlockListOfPoles(A, Nk, &E, G, Emin, Emax, NE);
  FreeCompactMatrix(&E);
  return ret;
}

uint32_t TightBindingHamiltonianPlusSelfEnergyToBlockListOfPolesGrid(TightBindingDef const A, int Nk[3], BlockAndersonMatrixType const Sigma, BlockListOfPoles * const G, double const * const Egrid, uint32_t const NEgrid)
{
  uint32_t ret;
  CompactMatrixType E;
  if(TestSelfEnergyBlockSize(Sigma.BlockSize, A.NFermion, "TightBindingHamiltonianPlusSelfEnergyToBlockListOfPolesGrid")) return 1;
  if(BlockAndersonMatrixInversedToCompactMatrix(&Sigma, &E)) return 1;
  ret = TightBindingHamiltonianPlusEmbeddedSelfEnergyToBlockListOfPolesGrid(A, Nk, &E, G, Egrid, NEgrid, NULL);
  FreeCompactMatrix(&E);
  return ret;
}

uint32_t TightBindingHamiltonianPlusSelfEnergyToBlockListOfPoles(TightBindingDef const A, int Nk[3], BlockNaturalImpurityOrbitalMatrixType const Sigma, BlockListOfPoles * const G, double const Emin, double const Emax, uint32_t const NE)
{
  uint32_t ret;
  CompactMatrixType E;
  if(TestSelfEnergyBlockSize(Sigma.BlockSize, A.NFermion, "TightBindingHamiltonianPlusSelfEnergyToBlockListOfPoles")) return 1;
  if(BlockNaturalImpurityOrbitalMatrixInversedToCompactMatrix(&Sigma, &E)) return 1;
  ret = TightBindingHamiltonianPlusEmbeddedSelfEnergyToBlockListOfPoles(A, Nk, &E, G, Emin, Emax, NE);
  FreeCompactMatrix(&E);
  return ret;
}

uint32_t TightBindingHamiltonianPlusSelfEnergyToBlockListOfPolesGrid(TightBindingDef const A, int Nk[3], BlockNaturalImpurityOrbitalMatrixType const Sigma, BlockListOfPoles * const G, double const * const Egrid, uint32_t const NEgrid)
{
  uint32_t ret;
  CompactMatrixType E;
  if(TestSelfEnergyBlockSize(Sigma.BlockSize, A.NFermion, "TightBindingHamiltonianPlusSelfEnergyToBlockListOfPolesGrid")) return 1;
  if(BlockNaturalImpurityOrbitalMatrixInversedToCompactMatrix(&Sigma, &E)) return 1;
  ret = TightBindingHamiltonianPlusEmbeddedSelfEnergyToBlockListOfPolesGrid(A, Nk, &E, G, Egrid, NEgrid, NULL);
  FreeCompactMatrix(&E);
  return ret;
}

uint32_t TightBindingHamiltonianPlusSelfEnergyToBlockListOfPoles(TightBindingDef const A, int Nk[3], BlockDoubleTriDiagonalMatrixType const Sigma, BlockListOfPoles * const G, double const Emin, double const Emax, uint32_t const NE)
{
  uint32_t ret;
  CompactMatrixType E;
  if(TestSelfEnergyBlockSize(Sigma.BlockSize, A.NFermion, "TightBindingHamiltonianPlusSelfEnergyToBlockListOfPoles")) return 1;
  if(BlockDoubleTriDiagonalMatrixInversedToCompactMatrix(&Sigma, &E)) return 1;
  ret = TightBindingHamiltonianPlusEmbeddedSelfEnergyToBlockListOfPoles(A, Nk, &E, G, Emin, Emax, NE);
  FreeCompactMatrix(&E);
  return ret;
}

uint32_t TightBindingHamiltonianPlusSelfEnergyToBlockListOfPolesGrid(TightBindingDef const A, int Nk[3], BlockDoubleTriDiagonalMatrixType const Sigma, BlockListOfPoles * const G, double const * const Egrid, uint32_t const NEgrid)
{
  uint32_t ret;
  CompactMatrixType E;
  if(TestSelfEnergyBlockSize(Sigma.BlockSize, A.NFermion, "TightBindingHamiltonianPlusSelfEnergyToBlockListOfPolesGrid")) return 1;
  if(BlockDoubleTriDiagonalMatrixInversedToCompactMatrix(&Sigma, &E)) return 1;
  ret = TightBindingHamiltonianPlusEmbeddedSelfEnergyToBlockListOfPolesGrid(A, Nk, &E, G, Egrid, NEgrid, NULL);
  FreeCompactMatrix(&E);
  return ret;
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
        FreeMatrix(&H);
        FreeCompactMatrix(&T);
      }
    }
  }
  FreeCompactMatrix(&HSigma);
  if(Sigma.QComplex==0) FreeCompactMatrix(&HSigmaComplex);
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
          if (iE >= (int)NE)
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
        FreeCompactMatrix(&T);
      }
    }
  }
  FreeMatrix(&PsiChar);
  FreeCompactMatrix(&HSigma);
  if(Sigma.QComplex==0) FreeCompactMatrix(&HSigmaComplex);
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
          if (iE >= (int)NE)
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
        FreeCompactMatrix(&T);
      }
    }
  }
  FreeMatrix(&PsiChar);
  FreeCompactMatrix(&HSigma);
  if(Sigma.QComplex==0) FreeCompactMatrix(&HSigmaComplex);
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
uint32_t DiagonalMatrixBlockTridiagonalize(double const *A, CompactMatrixType *psi, CompactMatrixType *const T, BlockTriDiagonalMatrixType *const ABlock, double NBlock, double SingularValue)
{
  return DiagonalMatrixBlockTridiagonalize(A, psi, T, ABlock, NBlock, SingularValue, NTriMaxNoLimit);
}

/*  NTriMax is a ceiling on the number of blocks of ABlock. NBlock is not that ceiling and
 *  never has been: it is handed in by every caller but was dropped here, psi->Nj going to
 *  CompactMatrixBlockBandDiagonalize in its place, and the callers pass a number (the
 *  count of poles) that is smaller than the chain a deflating block Lanczos can reach, so
 *  starting to honour it would shorten chains that other parts of the code rely on. It is
 *  left alone and NTriMax carries the ceiling instead.
 */
uint32_t DiagonalMatrixBlockTridiagonalize(double const *A, CompactMatrixType *psi, CompactMatrixType *const T, BlockTriDiagonalMatrixType *const ABlock, double NBlock, double SingularValue, uint32_t const NTriMax)
{
  int i;
  uint32_t NTri;
  CompactMatrixType M;
  M.QComplex=0;
  M.Ni=psi->Nj;
  M.Nj=psi->Nj;
  InitCompactMatrix(&M);
  for(i=0;i<psi->Nj;i++)
  {
    M.VRe[i + i * M.Ni] = A[i];
  }
  CompactMatrixCopy(*psi, T);
  /*  The block Lanczos cannot be asked for zero blocks, it writes b0 before its loop, so
   *  one block is the floor; a caller that wants none cuts that block away afterwards. */
  NTri = psi->Nj;
  if(NTriMax<NTri) NTri=NTriMax;
  if(NTri<1) NTri=1;
  CompactMatrixBlockBandDiagonalize(&M, T, ABlock, SingularValue, NTri, psi->Nj, true, true);
  FreeCompactMatrix(&M);
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

/*   G(w) = A0 + B0 ( sum_i B_i 1/(w-A_i) B_i^T ) B0^T
 *   only accept: A0 = 0 && B0 = 1
 *
 */
uint32_t BlockListOfPolesToTightBinding(BlockListOfPoles const * const M, TightBindingDef * TB)
{
  //TODO: BlockListOfPoles to TB
  
  return 0;
}

/* The onsite energy of each atom is given by A_i
 * The hopping between atom i and i+1 is given by B_i
 */
uint32_t BlockTriDiagonalMatrixToTightBinding(BlockTriDiagonalMatrixType const * const M, TightBindingDef * TB)
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

uint32_t BlockAndersonMatrixToTightBinding(BlockAndersonMatrixType const * const M, TightBindingDef * TB)
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

uint32_t BlockNaturalImpurityOrbitalMatrixToTightBinding(BlockNaturalImpurityOrbitalMatrixType const * const M, TightBindingDef * TB)
{
  //TODO: BlockNaturalImpurityOrbitalMatrixType to TB
  return 0;
}


/*  =========================================================================
 *  Direct sums
 *
 *  Given m response functions A[i] of the same type with block sizes n_i and an offset Offset[i] each, return
 *  the response function B of block size BlockSize that is block diagonal, with A[i] in the rows and columns
 *  Offset[i] ... Offset[i]+n_i-1. A[i] may be NULL, which is a zero response function. The blocks may not
 *  overlap; rows and columns that no A[i] covers are zero. Nothing is transformed, the numbers of A[i] are
 *  copied into place. Chains and baths of different length are padded with zero blocks, which couple to
 *  nothing and so do not change the function.
 *  The result is complex if any A[i] is complex. mu is taken from the first A[i] that is not NULL, the
 *  caller is expected to have checked that they are equal.
 *  ========================================================================= */

/*  Copy the n by n block (row major) that starts at element Astart of A into the Nj wide matrix B, with its
 *  first element at element Bstart of B. A is real (ARe) or complex (A) according to AQComplex and likewise
 *  B, B can only be real if A is real.
 */
static void CopyBlockIntoDirectSum(double * const BRe, Complex16 * const B, uint32_t const BQComplex, uint32_t const Bstart, uint32_t const Nj, double const * const ARe, Complex16 const * const A, uint32_t const AQComplex, uint32_t const Astart, uint32_t const n)
{
  uint32_t i,j;
  for(i=0;i<n;i++)
  {
    for(j=0;j<n;j++)
    {
      if(BQComplex)
      {
        if(AQComplex)
        {
          B[Bstart+i*Nj+j] = A[Astart+i*n+j];
        }
        else
        {
          B[Bstart+i*Nj+j].Re = ARe[Astart+i*n+j];
          B[Bstart+i*Nj+j].Im = 0;
        }
      }
      else
      {
        BRe[Bstart+i*Nj+j] = ARe[Astart+i*n+j];
      }
    }
  }
}

/*  For a list of poles the poles of all A[i] are put in one list, each with its weight in block i and zeros
 *  elsewhere. Poles of different A[i] at the same energy are then merged into one pole with a block diagonal
 *  weight by BlockListOfPolesMergeDegeneratePoles with its standard distance (a few DBL_EPSILON times the
 *  band width), which keeps the number of poles from growing with m when all A[i] are on one energy grid.
 */
uint32_t BlockListOfPolesDirectSum(BlockListOfPoles const * const * const A, uint32_t const * const Offset, uint32_t const m, uint32_t const BlockSize, BlockListOfPoles * const B)
{
  uint32_t i,p,q,n=BlockSize,stride=1+BlockSize*BlockSize;
  bool first=true;
  B->BlockSize = BlockSize;
  B->N = 0;
  B->QComplex = 0;
  for(i=0;i<m;i++)
  {
    if(A[i]==NULL) continue;
    B->N += A[i]->N;
    B->QComplex |= A[i]->QComplex;
  }
  if(InitBlockListOfPoles(B))
  {
    printf("InitBlockListOfPoles failed in BlockListOfPolesDirectSum\n");
    return 1;
  }
  strcpy(B->Name, "Direct sum");
  p=0;
  for(i=0;i<m;i++)
  {
    if(A[i]==NULL) continue;
    if(first) {B->mu = A[i]->mu; first=false;}
    for(q=0;q<A[i]->N;q++)
    {
      if(B->QComplex) B->V[p*stride].Re = (A[i]->QComplex) ? (A[i]->V[q*(1+Sqr(A[i]->BlockSize))].Re) : (A[i]->VRe[q*(1+Sqr(A[i]->BlockSize))]);
      else            B->VRe[p*stride]  = A[i]->VRe[q*(1+Sqr(A[i]->BlockSize))];
      CopyBlockIntoDirectSum(B->VRe, B->V, B->QComplex, p*stride+1+Offset[i]*n+Offset[i], n, A[i]->VRe, A[i]->V, A[i]->QComplex, q*(1+Sqr(A[i]->BlockSize))+1, A[i]->BlockSize);
      p++;
    }
    CopyBlockIntoDirectSum(B->a0Re, B->a0, B->QComplex, Offset[i]*n+Offset[i], n, A[i]->a0Re, A[i]->a0, A[i]->QComplex, 0, A[i]->BlockSize);
  }
  return BlockListOfPolesMergeDegeneratePoles(B, -1, false);
}

uint32_t BlockTriDiagonalMatrixDirectSum(BlockTriDiagonalMatrixType const * const * const A, uint32_t const * const Offset, uint32_t const m, uint32_t const BlockSize, BlockTriDiagonalMatrixType * const B)
{
  uint32_t i,k,n=BlockSize;
  bool first=true;
  B->BlockSize = BlockSize;
  B->N = 0;
  B->QComplex = 0;
  B->mu = 0;
  for(i=0;i<m;i++)
  {
    if(A[i]==NULL) continue;
    if(A[i]->N > B->N) B->N = A[i]->N;
    B->QComplex |= A[i]->QComplex;
  }
  if(InitBlockTridiagonalMatrix(B))
  {
    printf("InitBlockTridiagonalMatrix failed in BlockTriDiagonalMatrixDirectSum\n");
    return 1;
  }
  strcpy(B->Name, "Direct sum");
  // InitBlockTridiagonalMatrix sets b0 to the identity, here it is built from the blocks of A
  for(k=0;k<n*n;k++)
  {
    if(B->N==0) break;
    if(B->QComplex) {B->b[k].Re=0; B->b[k].Im=0;}
    else            B->bRe[k]=0;
  }
  for(i=0;i<m;i++)
  {
    if(A[i]==NULL) continue;
    if(first) {B->mu = A[i]->mu; first=false;}
    for(k=0;k<=A[i]->N;k++)
    {
      CopyBlockIntoDirectSum(B->aRe, B->a, B->QComplex, k*n*n+Offset[i]*n+Offset[i], n, A[i]->aRe, A[i]->a, A[i]->QComplex, k*Sqr(A[i]->BlockSize), A[i]->BlockSize);
    }
    for(k=0;k<A[i]->N;k++)
    {
      CopyBlockIntoDirectSum(B->bRe, B->b, B->QComplex, k*n*n+Offset[i]*n+Offset[i], n, A[i]->bRe, A[i]->b, A[i]->QComplex, k*Sqr(A[i]->BlockSize), A[i]->BlockSize);
    }
  }
  return 0;
}

uint32_t BlockAndersonMatrixDirectSum(BlockAndersonMatrixType const * const * const A, uint32_t const * const Offset, uint32_t const m, uint32_t const BlockSize, BlockAndersonMatrixType * const B)
{
  uint32_t i,k,n=BlockSize;
  bool first=true;
  B->BlockSize = BlockSize;
  B->N = 0;
  B->QComplex = 0;
  B->mu = 0;
  for(i=0;i<m;i++)
  {
    if(A[i]==NULL) continue;
    if(A[i]->N > B->N) B->N = A[i]->N;
    B->QComplex |= A[i]->QComplex;
  }
  if(InitBlockAndersonMatrix(B))
  {
    printf("InitBlockAndersonMatrix failed in BlockAndersonMatrixDirectSum\n");
    return 1;
  }
  strcpy(B->Name, "Direct sum");
  // InitBlockAndersonMatrix sets b0 to the identity, here it is built from the blocks of A
  for(k=0;k<n*n;k++)
  {
    if(B->N==0) break;
    if(B->QComplex) {B->b[k].Re=0; B->b[k].Im=0;}
    else            B->bRe[k]=0;
  }
  for(i=0;i<m;i++)
  {
    if(A[i]==NULL) continue;
    if(first) {B->mu = A[i]->mu; first=false;}
    for(k=0;k<=A[i]->N;k++)
    {
      CopyBlockIntoDirectSum(B->aRe, B->a, B->QComplex, k*n*n+Offset[i]*n+Offset[i], n, A[i]->aRe, A[i]->a, A[i]->QComplex, k*Sqr(A[i]->BlockSize), A[i]->BlockSize);
    }
    for(k=0;k<A[i]->N;k++)
    {
      CopyBlockIntoDirectSum(B->bRe, B->b, B->QComplex, k*n*n+Offset[i]*n+Offset[i], n, A[i]->bRe, A[i]->b, A[i]->QComplex, k*Sqr(A[i]->BlockSize), A[i]->BlockSize);
    }
  }
  return 0;
}

/*  direct sum of the n_i by n_i compact matrices A[i] (NULL is zero) into the BlockSize by BlockSize matrix B */
static uint32_t CompactMatrixDirectSum(CompactMatrixType const * const * const A, uint32_t const * const Offset, uint32_t const m, uint32_t const BlockSize, CompactMatrixType * const B)
{
  uint32_t i;
  B->Ni = BlockSize;
  B->Nj = BlockSize;
  B->QComplex = 0;
  for(i=0;i<m;i++)
  {
    if(A[i]!=NULL) B->QComplex |= A[i]->QComplex;
  }
  if(InitCompactMatrix(B))
  {
    printf("InitCompactMatrix failed in CompactMatrixDirectSum\n");
    return 1;
  }
  for(i=0;i<m;i++)
  {
    if(A[i]==NULL) continue;
    CopyBlockIntoDirectSum(B->VRe, B->V, B->QComplex, Offset[i]*BlockSize+Offset[i], BlockSize, A[i]->VRe, A[i]->V, A[i]->QComplex, 0, A[i]->Ni);
  }
  return 0;
}

/*  a0, val and con are summed, TMT is calculated from the summed val and con */
uint32_t BlockNaturalImpurityOrbitalMatrixDirectSum(BlockNaturalImpurityOrbitalMatrixType const * const * const A, uint32_t const * const Offset, uint32_t const m, uint32_t const BlockSize, BlockNaturalImpurityOrbitalMatrixType * const B)
{
  uint32_t i, ret=0;
  bool first=true;
  BlockTriDiagonalMatrixType const **val = (BlockTriDiagonalMatrixType const **)malloc(m*sizeof(BlockTriDiagonalMatrixType*));
  BlockTriDiagonalMatrixType const **con = (BlockTriDiagonalMatrixType const **)malloc(m*sizeof(BlockTriDiagonalMatrixType*));
  CompactMatrixType const **a0 = (CompactMatrixType const **)malloc(m*sizeof(CompactMatrixType*));
  if((val==NULL)||(con==NULL)||(a0==NULL))
  {
    printf("malloc failed in BlockNaturalImpurityOrbitalMatrixDirectSum\n");
    return 1;
  }
  B->BlockSize = BlockSize;
  B->QComplex = 0;
  B->mu = 0;
  B->Epsilon = 0;
  for(i=0;i<m;i++)
  {
    val[i] = (A[i]==NULL) ? (NULL) : (&(A[i]->val));
    con[i] = (A[i]==NULL) ? (NULL) : (&(A[i]->con));
    a0[i]  = (A[i]==NULL) ? (NULL) : (&(A[i]->a0));
    if(A[i]==NULL) continue;
    if(first) {B->mu = A[i]->mu; B->Epsilon = A[i]->Epsilon; first=false;}
    B->QComplex |= A[i]->QComplex;
  }
  strcpy(B->Name, "Direct sum");
  if(BlockTriDiagonalMatrixDirectSum(val, Offset, m, BlockSize, &(B->val))) ret=1;
  if(BlockTriDiagonalMatrixDirectSum(con, Offset, m, BlockSize, &(B->con))) ret=1;
  if(CompactMatrixDirectSum(a0, Offset, m, BlockSize, &(B->a0))) ret=1;
  if(ret==0) ret = SetTMTFromValandConInBlockNaturalImpurityOrbitalMatrixType(B);
  free(val);
  free(con);
  free(a0);
  return ret;
}

/*  a0, b0, val and con are summed. kBT, only used to generate the representation, is taken from the first */
uint32_t BlockDoubleTriDiagonalMatrixDirectSum(BlockDoubleTriDiagonalMatrixType const * const * const A, uint32_t const * const Offset, uint32_t const m, uint32_t const BlockSize, BlockDoubleTriDiagonalMatrixType * const B)
{
  uint32_t i, ret=0;
  bool first=true;
  BlockTriDiagonalMatrixType const **val = (BlockTriDiagonalMatrixType const **)malloc(m*sizeof(BlockTriDiagonalMatrixType*));
  BlockTriDiagonalMatrixType const **con = (BlockTriDiagonalMatrixType const **)malloc(m*sizeof(BlockTriDiagonalMatrixType*));
  CompactMatrixType const **a0 = (CompactMatrixType const **)malloc(m*sizeof(CompactMatrixType*));
  CompactMatrixType const **b0 = (CompactMatrixType const **)malloc(m*sizeof(CompactMatrixType*));
  if((val==NULL)||(con==NULL)||(a0==NULL)||(b0==NULL))
  {
    printf("malloc failed in BlockDoubleTriDiagonalMatrixDirectSum\n");
    return 1;
  }
  B->BlockSize = BlockSize;
  B->mu = 0;
  B->kBT = 0;
  for(i=0;i<m;i++)
  {
    val[i] = (A[i]==NULL) ? (NULL) : (&(A[i]->val));
    con[i] = (A[i]==NULL) ? (NULL) : (&(A[i]->con));
    a0[i]  = (A[i]==NULL) ? (NULL) : (&(A[i]->a0));
    b0[i]  = (A[i]==NULL) ? (NULL) : (&(A[i]->b0));
    if(A[i]==NULL) continue;
    if(first) {B->mu = A[i]->mu; B->kBT = A[i]->kBT; first=false;}
  }
  strcpy(B->Name, "Direct sum");
  if(BlockTriDiagonalMatrixDirectSum(val, Offset, m, BlockSize, &(B->val))) ret=1;
  if(BlockTriDiagonalMatrixDirectSum(con, Offset, m, BlockSize, &(B->con))) ret=1;
  if(CompactMatrixDirectSum(a0, Offset, m, BlockSize, &(B->a0))) ret=1;
  if(CompactMatrixDirectSum(b0, Offset, m, BlockSize, &(B->b0))) ret=1;
  SetBlockDoubleTriDiagonalMatrixQComplex(B);
  free(val);
  free(con);
  free(a0);
  free(b0);
  return ret;
}
