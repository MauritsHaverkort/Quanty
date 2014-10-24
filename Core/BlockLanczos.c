/*  given an Determinant and a set of Restrictions
 *  If I leave iFermion empty in the Determinant and I will be allowed to make changes to all fermions at higher index positions will I
 *  A) still be able to fullfill the restrictions - return 0
 *  B) not be able to fullfill the restrictions - return 1
 */
uint32_t QMustOccupy(uint32_t const iFermion, uint8_t *Det, DeterminantRestrictionType const Restrictions)
{
  printf("to do\n");
  return 0;
}

/*  given an Determinant and a set of Restrictions
 *  If I add an electron to iFermion in the Determinant and I will be allowed to make changes to all fermions at higher index positions will I
 *  A) still be able to fullfill the restrictions - return 0
 *  B) not be able to fullfill the restrictions - return 1
 */
uint32_t QCanNotOccupy(uint32_t const iFermion, uint8_t *Det, DeterminantRestrictionType const Restrictions)
{
  printf("to do\n");
  return 0;
}

/*  given a set of restrictions create Npsi orthonormal functions on a reasonably large basis that can serve as a starting point for
 *  finding the ground-state of H
 */
uint32_t CreateStartingFunctions(const OperatorType H, uint32_t const Npsi, DeterminantRestrictionType const StartingPoint, WaveFunctionType * psi, uint32_t const NBitsKey)
{
  double Norm;
  WaveFunctionType psi0;
  uint8_t *Determinant;
  uint32_t i,j,iFermion,MustOccupy,CannotOccupy;
  // if nothing has to be done return
  if(Npsi==0) return 0;
  // allocate psi0
  psi0.NFermion=H.NFermion;
  psi0.NBoson=H.NBoson;
  psi0.QComplex=(H.QComplex==0) ? 0 : 1;
  psi0.NBitsKey=NBitsKey;
  if(WaveFunctionInit(&(psi0)))
  {
    printf("WaveFunctionInit failed in CreateStartingFunctions\n");
    return 1;
  }
  // create some determinants fullfilling the StartingPoint conditions
  // "some" is equal to 10+3*Npsi
  Determinant = calloc(psi0.NWord, sizeof(uint8_t));
  for(i=0;i<10+3*Npsi;i++)
  {
    for(iFermion=0;iFermion<H.NFermion;iFermion++)
    {
      MustOccupy=QMustOccupy(iFermion, Determinant, StartingPoint);
      CannotOccupy=QCanNotOccupy(iFermion, Determinant, StartingPoint);
      if(MustOccupy&&CannotOccupy)
      {
        printf("Cannot fullfill all restrictions in CreateStartingFunctions\n");
        free(Determinant);
        return 1;
      }
      if(MustOccupy)
      {
        OperatorCr(Determinant, i, psi->NFermion);
      }
      else if(CannotOccupy)
      {
        // do nothing
      }
      else
      {
        if(rand()&1) OperatorCr(Determinant, i, psi->NFermion);
      }
    }
    WaveFunctionAddElement(&psi0, Determinant, 1, 0);
    for(j=0;j<psi0.NWord;j++)
    {
      Determinant[j]=0;
    }
  }
  // randomize psi0
  WaveFunctionRandomize(&(psi0));
  // allocate psi
  psi=malloc(Npsi*sizeof(WaveFunctionType));
  // add states by making psi = O psi0
  if(OperatorPsiRestricted(H, psi0, &(psi[0]), StartingPoint))
  {
    printf("OperatorPsiRestricted failed in CreateStartingFunctions\n");
    return 1;
  }
  // init other functions
  for(i=1;i<Npsi;i++)
  {
    if(WaveFunctionInitCopyBasis(&(psi[i]), (psi[0])))
    {
      printf("WaveFunctionInitCopyBasis failed in CreateStartingFunctions\n");
      return 1;
    }
  }
  // randomize and orthonormalize psi
  for(i=1;i<Npsi;i++)
  {
    WaveFunctionRandomize(&(psi[i]));
    WaveFunctionOrthogonalizeSameBasis(&(psi[i]), &(psi[0]), i);
    WaveFunctionCalculateNorm(psi[i], &Norm);
    if(Norm<1000*DBL_EPSILON)
    {
      WaveFunctionMultiply(&(psi[i]), 0, 0);
    }
    else
    {
      WaveFunctionMultiply(&(psi[i]), 1/Norm, 0);
    }
  }
  return 0;
}


/*  Calculate Opsi[i] = H psi[i] for i=0 to Npsi
 *  on input all psi[i]'s are on the same basis
 *  on excit all Opsi[i]'s are on the same basis
 */
uint32_t BlockOperatorPsiRestricted(const OperatorType H, WaveFunctionType const * const psi, WaveFunctionType * const Opsi, uint32_t const Npsi, DeterminantRestrictionType const Restrictions)
{
  uint32_t i,j;
  if(OperatorPsiRestricted(H, psi[0], &(Opsi[0]), Restrictions))
  {
    printf("OperatorPsi failed in BlockOperatorPsi\n");
    return 1;
  }
  for(i=1;i<Npsi;i++)
  {
    if(WaveFunctionInitCopyBasis(&(Opsi[i]), (Opsi[0])))
    {
      for(j=1;j<i;j++)
      {
        WaveFunctionFreeValues(&Opsi[j]);
      }
      WaveFunctionFree(&Opsi[0]);
      free(Opsi);
      printf("WaveFunctionInitCopyBasis failed in BlockOperatorPsi\n");
      return 1;
    }
  }
  for(i=1;i<Npsi;i++)
  {
    OperatorPsiConserveBasisPreAllocated(&H, &(psi[i]), &(Opsi[i]));
  }
  return 0;
}

/*  Calculate Opsi[i] = H psi[i] for i=0 to Npsi
 *  on input all psi[i]'s are on the same basis
 *  on excit all Opsi[i]'s are on the same basis
 */
uint32_t BlockOperatorPsi(const OperatorType H, WaveFunctionType const * const psi, WaveFunctionType * const Opsi, uint32_t const Npsi)
{
  uint32_t i,j;
  if(OperatorPsi(H, psi[0], &(Opsi[0])))
  {
    printf("OperatorPsi failed in BlockOperatorPsi\n");
    return 1;
  }
  for(i=1;i<Npsi;i++)
  {
    if(WaveFunctionInitCopyBasis(&(Opsi[i]), (Opsi[0])))
    {
      for(j=1;j<i;j++)
      {
        WaveFunctionFreeValues(&Opsi[j]);
      }
      WaveFunctionFree(&Opsi[0]);
      free(Opsi);
      printf("WaveFunctionInitCopyBasis failed in BlockOperatorPsi\n");
      return 1;
    }
  }
  for(i=1;i<Npsi;i++)
  {
    OperatorPsiConserveBasisPreAllocated(&H, &(psi[i]), &(Opsi[i]));
  }
  return 0;
}

void BlockWaveFunctionFree(WaveFunctionType *psi, const uint32_t Npsi)
{
  uint32_t i;
  for(i=1;i<Npsi;i++)
  {
    WaveFunctionFreeValues(&(psi[i]));
  }
  WaveFunctionFree(&(psi[0]));
}

void RealBlockWaveFunctionTruncateBasis(WaveFunctionType *const psi, uint32_t const Npsi, const double epsilon)
{
  uint32_t startN, iBlok, iElement, ipsi, KeepDeterminant;
  for(ipsi=0;ipsi<Npsi;ipsi++)
  {
    psi[ipsi].ValueRe[0][0] = 0;
  }
  startN = psi[0].N;
  psi[0].N = 0;
  for (iElement = 0; iElement < ((uint32_t) 1 << (psi[0].NBitsKey)); ++iElement)
  {                             /* Whipe out the entire KeyList */
    psi[0].KeyList[iElement] = 0;
  }
  for (iBlok = 0; iBlok < psi[0].NBlok - 1; iBlok++)
  {
    for (iElement = 0; iElement < BlokSize; iElement++)
    {
      psi[0].Next[iBlok][iElement] = 0;
      KeepDeterminant=0;
      for(ipsi=0;ipsi<Npsi;ipsi++)
      {
        if(Abs(psi[ipsi].ValueRe[iBlok][iElement])<epsilon)
        {
          psi[ipsi].ValueRe[iBlok][iElement]=0;
        }
        else
        {
          KeepDeterminant=1;
        }
      }
      if (KeepDeterminant)
      {
        RealWaveFunctionAddElement(&(psi[0]), &(psi[0].Determinant[iBlok][iElement * psi[0].NWord]), psi[0].ValueRe[iBlok][iElement]);
        for(ipsi=1;ipsi<Npsi;ipsi++)
        {
          psi[ipsi].ValueRe[BlokElement(psi[0].N)] = psi[ipsi].ValueRe[iBlok][iElement];
        }
      }
    }
  }
  iBlok = psi[0].NBlok - 1;
  for (iElement = 0; iElement < 1 + startN - iBlok * BlokSize; ++iElement)
  {
    psi[0].Next[iBlok][iElement] = 0;
    KeepDeterminant=0;
    for(ipsi=0;ipsi<Npsi;ipsi++)
    {
      if(Abs(psi[ipsi].ValueRe[iBlok][iElement])<epsilon)
      {
        psi[ipsi].ValueRe[iBlok][iElement]=0;
      }
      else
      {
        KeepDeterminant=1;
      }
    }
    if (KeepDeterminant)
    {
      RealWaveFunctionAddElement(&(psi[0]), &(psi[0].Determinant[iBlok][iElement * psi[0].NWord]), psi[0].ValueRe[iBlok][iElement]);
      for(ipsi=1;ipsi<Npsi;ipsi++)
      {
        psi[ipsi].ValueRe[BlokElement(psi[0].N)] = psi[ipsi].ValueRe[iBlok][iElement];
      }
    }
  }
  for (iBlok = 1 + (psi[0].N >> NBitsBlok); iBlok < psi[0].NBlok; iBlok++)
  {
    free(psi[0].Next[iBlok]);
    free(psi[0].Determinant[iBlok]);
    for(ipsi=0;ipsi<Npsi;ipsi++)
    {
      free(psi[ipsi].ValueRe[iBlok]);
    }
  }
  psi[0].NBlok = 1 + (psi[0].N >> NBitsBlok);
  psi[0].Next = realloc(psi[0].Next, psi[0].NBlok * sizeof(uint32_t *));
  psi[0].Determinant = realloc(psi[0].Determinant, psi[0].NBlok * sizeof(uint8_t *));
  for(ipsi=0;ipsi<Npsi;ipsi++)
  {
    psi[ipsi].ValueRe = realloc(psi[ipsi].ValueRe, psi[0].NBlok * sizeof(double *));
  }
  for(ipsi=1;ipsi<Npsi;ipsi++)
  {
    psi[ipsi].N           = psi[0].N;
    psi[ipsi].NBlok       = psi[0].NBlok;
    psi[ipsi].Next        = psi[0].Next;
    psi[ipsi].Determinant = psi[0].Determinant;
  }
}

void ComplexBlockWaveFunctionTruncateBasis(WaveFunctionType *const psi, uint32_t const Npsi, const double epsilon)
{
  double const epsilonSqr=Sqr(epsilon);
  uint32_t startN, iBlok, iElement, ipsi, KeepDeterminant;
  for(ipsi=0;ipsi<Npsi;ipsi++)
  {
    psi[ipsi].ValueRe[0][0] = 0;
    psi[ipsi].ValueIm[0][0] = 0;
  }
  startN = psi[0].N;
  psi[0].N = 0;
  for (iElement = 0; iElement < ((uint32_t) 1 << (psi[0].NBitsKey)); ++iElement)
  {                             /* Whipe out the entire KeyList */
    psi[0].KeyList[iElement] = 0;
  }
  for (iBlok = 0; iBlok < psi[0].NBlok - 1; iBlok++)
  {
    for (iElement = 0; iElement < BlokSize; iElement++)
    {
      psi[0].Next[iBlok][iElement] = 0;
      KeepDeterminant=0;
      for(ipsi=0;ipsi<Npsi;ipsi++)
      {
        if((Sqr(psi[ipsi].ValueRe[iBlok][iElement])+Sqr(psi[ipsi].ValueIm[iBlok][iElement]))<epsilonSqr)
        {
          psi[ipsi].ValueRe[iBlok][iElement]=0;
          psi[ipsi].ValueIm[iBlok][iElement]=0;
        }
        else
        {
          KeepDeterminant=1;
        }
      }
      if (KeepDeterminant)
      {
        ComplexWaveFunctionAddElement(&(psi[0]), &(psi[0].Determinant[iBlok][iElement * psi[0].NWord]), psi[0].ValueRe[iBlok][iElement], psi[0].ValueIm[iBlok][iElement]);
        for(ipsi=1;ipsi<Npsi;ipsi++)
        {
          psi[ipsi].ValueRe[BlokElement(psi[0].N)] = psi[ipsi].ValueRe[iBlok][iElement];
          psi[ipsi].ValueIm[BlokElement(psi[0].N)] = psi[ipsi].ValueIm[iBlok][iElement];
        }
      }
    }
  }
  iBlok = psi[0].NBlok - 1;
  for (iElement = 0; iElement < 1 + startN - iBlok * BlokSize; ++iElement)
  {
    psi[0].Next[iBlok][iElement] = 0;
    KeepDeterminant=0;
    for(ipsi=0;ipsi<Npsi;ipsi++)
    {
      if((Sqr(psi[ipsi].ValueRe[iBlok][iElement])+Sqr(psi[ipsi].ValueIm[iBlok][iElement]))<epsilonSqr)
      {
        psi[ipsi].ValueRe[iBlok][iElement]=0;
        psi[ipsi].ValueIm[iBlok][iElement]=0;
      }
      else
      {
        KeepDeterminant=1;
      }
    }
    if (KeepDeterminant)
    {
      ComplexWaveFunctionAddElement(&(psi[0]), &(psi[0].Determinant[iBlok][iElement * psi[0].NWord]), psi[0].ValueRe[iBlok][iElement], psi[0].ValueIm[iBlok][iElement]);
      for(ipsi=1;ipsi<Npsi;ipsi++)
      {
        psi[ipsi].ValueRe[BlokElement(psi[0].N)] = psi[ipsi].ValueRe[iBlok][iElement];
        psi[ipsi].ValueIm[BlokElement(psi[0].N)] = psi[ipsi].ValueIm[iBlok][iElement];
      }
    }
  }
  for (iBlok = 1 + (psi[0].N >> NBitsBlok); iBlok < psi[0].NBlok; iBlok++)
  {
    free(psi[0].Next[iBlok]);
    free(psi[0].Determinant[iBlok]);
    for(ipsi=0;ipsi<Npsi;ipsi++)
    {
      free(psi[ipsi].ValueRe[iBlok]);
      free(psi[ipsi].ValueIm[iBlok]);
    }
  }
  psi[0].NBlok = 1 + (psi[0].N >> NBitsBlok);
  psi[0].Next = realloc(psi[0].Next, psi[0].NBlok * sizeof(uint32_t *));
  psi[0].Determinant = realloc(psi[0].Determinant, psi[0].NBlok * sizeof(uint8_t *));
  for(ipsi=0;ipsi<Npsi;ipsi++)
  {
    psi[ipsi].ValueRe = realloc(psi[ipsi].ValueRe, psi[0].NBlok * sizeof(double *));
    psi[ipsi].ValueIm = realloc(psi[ipsi].ValueIm, psi[0].NBlok * sizeof(double *));
  }
  for(ipsi=1;ipsi<Npsi;ipsi++)
  {
    psi[ipsi].N           = psi[0].N;
    psi[ipsi].NBlok       = psi[0].NBlok;
    psi[ipsi].Next        = psi[0].Next;
    psi[ipsi].Determinant = psi[0].Determinant;
  }
}

void BlockWaveFunctionTruncateBasis(WaveFunctionType *const psi, uint32_t const Npsi, const double epsilon)
{
  if(psi[0].QComplex)
  {
    ComplexBlockWaveFunctionTruncateBasis(psi, Npsi, epsilon);
  }
  else
  {
    RealBlockWaveFunctionTruncateBasis(psi, Npsi, epsilon);
  }
}

/*  use a dense mehtod to find the Npsi lowest states of Operator H
 *  The basis is given by the psi[0].N determinants
 *  psi[0] to psi[Npsi-1] are allocated and live on the same basis
 *  on output psi[i].V is changed to represent the Npsi lowest eigen states
 */
uint32_t BlockLanczosGroundStateConserveBasisDense(const OperatorType H, WaveFunctionType * const psi, uint32_t const Npsi)
{
  uint32_t i;
  CompactMatrixType Hdens;
  Hdens.QComplex = ((H.QComplex)||(psi[0].QComplex));
  Hdens.Ni=psi[0].N;
  Hdens.Nj=psi[0].N;
  InitCompactMatrix(&Hdens);
  // fill the matrix Hdens
  
  WaveFunctionType Hpsi;
  WaveFunctionInitCopyBasis(&Hpsi, *psi);
  for(i=0;i<psi[0].N;i++)
  {
    
  }
  WaveFunctionFreeValues(&Hpsi);
  
  
  return 0;
}

/*  Use Krylov Lanczos method to find lowest Npsi eigenstates of H
 *  Starting vectors are in v[0]
 *  v[1] to v[NKrylovStart] are allocated
 *  after a restart increase KrylovBasis by NKrylovStep if one is able to allocate the needed v's for this.
 *  all states on same basis
 *
 *  on excess free memory of v[i], except v[0] which contains the Npsi eigen states
 */
uint32_t BlockLanczosGroundStateConserveBasisKrylov(const OperatorType H, WaveFunctionType ** const v, uint32_t const Npsi, const double zero)
{
  printf("to do\n");
  return 0;
}

/*  Use Krylov Lanczos method to find lowest Npsi eigenstates of H
 *  Starting vectors are in v[0]
 *  v[1] and v[2] are allocated
 *  all states on same basis
 *
 *  on excess free memory of v[1] and v[2], except v[0] which contains the Npsi eigen states
 */
uint32_t BlockLanczosGroundStateConserveBasisKrylov3(const OperatorType H, WaveFunctionType ** const v, uint32_t const Npsi, const double zero)
{
  printf("to do\n");
  return 0;
}


/*   Finds in the basis given by psi the lowest Npsi eigenstates of H
 *   Input:
 *     Npsi:     number of eigenstates to be found
 *     psi[0].N: number of determinants in basis
 *     zero:     smallest positive real to be considered different from 0 in convergent tests
 *
 *   Internal constants:
 *     DenseBoarder:  (the number of basis states below which a dense calculation is done)
 *     NKrlylovStart: (the number of KrylovVectors to start with)
 *
 *   if psi[0].N < DenseBoarder do a dense calculation
 *   if Npsi > psi[0]/2 do a dense calculation
 *
 */
uint32_t BlockLanczosGroundStateConserveBasis(const OperatorType H, WaveFunctionType * const psi, uint32_t const Npsi, const double zero)
{
  uint32_t DenseBoarder=100;
  uint32_t NKrylovStart=15;
  uint32_t NKrylovMax=100;
  
  uint32_t i,j,k,l;
  WaveFunctionType *v[NKrylovMax];
  // for small basises do a dense calculation
  if(psi[0].N<=DenseBoarder) return BlockLanczosGroundStateConserveBasisDense(H,psi, Npsi);
  // create the first 3 Block vectors in memory. If this does not work return 1
  v[0]=psi;
  v[1]=malloc(Npsi*sizeof(WaveFunctionType));
  v[2]=malloc(Npsi*sizeof(WaveFunctionType));
  for(i=0;i<Npsi;i++)
  {
    if(WaveFunctionInitCopyBasis(&(v[0][i]), (v[1][i])))
    {
      for(j=0;j<i;j++)
      {
        WaveFunctionFreeValues(&v[1][j]);
      }
      free(v[1]);
      free(v[2]);
      printf("WaveFunctionInitCopyBasis failed in BlockLanczosGroundStateConserveBasis for first block (%i-th function)\n",i);
      return 1;
    }
  }
  for(i=0;i<Npsi;i++)
  {
    if(WaveFunctionInitCopyBasis(&(v[0][i]), (v[2][i])))
    {
      for(j=0;j<i;j++)
      {
        WaveFunctionFreeValues(&v[2][j]);
      }
      for(j=0;j<Npsi;j++)
      {
        WaveFunctionFreeValues(&v[1][j]);
      }
      free(v[1]);
      free(v[2]);
      printf("WaveFunctionInitCopyBasis failed in BlockLanczosGroundStateConserveBasis for second block (%i-th function)\n",i);
      return 1;
    }
  }
  // create block 3 to block NKrylovStart if this does not work use Krylov recreating the states if it works use Krylov remembering the states
  for(k=3;k<NKrylovStart;k++)
  {
    for(i=0;i<Npsi;i++)
    {
      if(WaveFunctionInitCopyBasis(&(v[0][i]), (v[k][i])))
      {
        for(j=0;j<i;j++)
        {
          WaveFunctionFreeValues(&v[k][j]);
        }
        for(l=3;l<k;l++)
        {
          for(j=0;j<Npsi;j++)
          {
            WaveFunctionFreeValues(&v[l][j]);
          }
        }
        return BlockLanczosGroundStateConserveBasisKrylov3(H, v, Npsi, zero);
      }
    }
  }
  return BlockLanczosGroundStateConserveBasisKrylov(H, v, Npsi, zero);
}

/*  Find the lowest Npsi eigen-states of the operator H
 *
 *  Input:
 *    H:             Operator
 *    Npsi:          number of states to converge
 *    zero:          smallest positive real larger than zero considered to be equal to zero in convergence tests
 *    epsilon:       smallest norm of prefactor of a determinant kept in a wave function
 *    StartingPoint: a set of restrictions defining the starting point wave function
 *    psi:           a pointer to Npsi uninitialized functions
 *    NBitsKey:      The number of bits in the key
 *
 *  Output:
 *    psi:           a list of Npsi wave-functions representing the Npsi lowest eigen-states
 *    epsilon:       will change if not enough memory is available.
 *
 *  The restricted version uses the DeterminantRestrictions always, the full version first calls the restricted version to find a starting point
 */
uint32_t BlockGroundStateRestricted(const OperatorType H, uint32_t const Npsi, const double zero, double * const epsilon, DeterminantRestrictionType const Restrictions, WaveFunctionType * const psi, const uint32_t NBitsKey)
{
  uint32_t i,Qconverged=0, AtLeastAnOtherLoop=0;
  WaveFunctionType Opsi[Npsi];
  double EnergyOld[Npsi],Energy[Npsi],EnergySqr[Npsi];
  // create Npsi independent functions fullfiling the starting point conditions
  // all functions must have the same basis
  // if H is complex or mixed psi will be complex, if H is real psi is real
  if(CreateStartingFunctions(H, Npsi, Restrictions, psi, NBitsKey))
  {
    printf("CreateStartingFunctions failed in BlockGroundState\n");
    return 1;
  }
  // create starting energy and variance, and check for convergence. Do not check for convergence here, but only in the unrestricted version
  Qconverged=0;
  for(i=0;i<Npsi;i++)
  {
    PsiHermitianOperatorPsi(H, psi[i], &(Energy[i]));
    PsiHermitianOperatorSqrPsi(H, psi[i], &(EnergySqr[i]));
    if((EnergySqr[i]-Sqr(Energy[i]))>zero) Qconverged=0;
  }
  // loop until converged
  while(!Qconverged)
  {
    // expand basis
    while(BlockOperatorPsiRestricted(H, psi, Opsi, Npsi, Restrictions))
    {
      // was not able to evaluate O psi, need to truncate psi further
      if(psi[0].NBlok==1)
      {
        printf("Run outof memory need to truncate but NBlok=1\n You're screwed\n");
        return 1;
      }
      *epsilon/=2;
      BlockWaveFunctionTruncateBasis(psi,Npsi,*epsilon);
    }
    BlockWaveFunctionFree(psi,Npsi);
    for(i=0;i<Npsi;i++)
    {
      psi[i]=Opsi[i];
    }
    // find ground-state in basis
    while(BlockLanczosGroundStateConserveBasis(H, psi, Npsi, zero))
    {
      AtLeastAnOtherLoop=1;
      printf("Run outof memory in BlockLanczosGroundStateConserveBasis, that's not good an no good solution is implemented\n");
      // was not able to evaluate ground-state in basis
      if(psi[0].NBlok==1)
      {
        printf("Run outof memory need to truncate but NBlok=1\n You're screwed\n");
        return 1;
      }
      *epsilon*=2;
      BlockWaveFunctionTruncateBasis(psi,Npsi,*epsilon);
    }
    // truncate basis
    BlockWaveFunctionTruncateBasis(psi,Npsi,*epsilon);
    // you're converged when the variance of all states is small or when the energy gain is close to zero
    Qconverged=1;
    for(i=0;i<Npsi;i++)
    {
      EnergyOld[i]=Energy[i];
      PsiHermitianOperatorPsi(H, psi[i], &(Energy[i]));
      PsiHermitianOperatorSqrPsiRestricted(H, psi[i], &(EnergySqr[i]), Restrictions);
      if((EnergySqr[i]-Sqr(Energy[i]))>zero) Qconverged=0;
    }
    if(Qconverged==0)
    {
      Qconverged=1;
      for(i=0;i<Npsi;i++)
      {
        if((EnergyOld[i]-Energy[i])>zero) Qconverged=0;
      }
    }
    if(AtLeastAnOtherLoop)
    {
      Qconverged=0;
      AtLeastAnOtherLoop=0;
    }
  }
  return 0;
}
uint32_t BlockGroundState(const OperatorType H, uint32_t const Npsi, const double zero, double * const epsilon, DeterminantRestrictionType const StartingPoint, WaveFunctionType * const psi, const uint32_t NBitsKey)
{
  uint32_t i,Qconverged=0, AtLeastAnOtherLoop=0;
  WaveFunctionType Opsi[Npsi];
  double EnergyOld[Npsi],Energy[Npsi],EnergySqr[Npsi];
  // First create the eigen states using the StartingPoint restrictions
  if(BlockGroundStateRestricted(H, Npsi, zero, epsilon, StartingPoint, psi, NBitsKey))
  {
    printf("BlockGroundStateRestricted failed in BlockGroundState\n");
    return 1;
  }
  // then free those restrictions and find the eigen states in the full space
  // create starting energy and variance, and check for convergence. If converged do nothing, if not converged run loop till convergence is reached
  Qconverged=1;
  for(i=0;i<Npsi;i++)
  {
    PsiHermitianOperatorPsi(H, psi[i], &(Energy[i]));
    PsiHermitianOperatorSqrPsi(H, psi[i], &(EnergySqr[i]));
    if((EnergySqr[i]-Sqr(Energy[i]))>zero) Qconverged=0;
  }
  // loop until converged
  while(!Qconverged)
  {
    // expand basis
    while(BlockOperatorPsi(H, psi, Opsi, Npsi))
    {
      // was not able to evaluate O psi, need to truncate psi further
      if(psi[0].NBlok==1)
      {
        printf("Run outof memory need to truncate but NBlok=1\n You're screwed\n");
        return 1;
      }
      *epsilon/=2;
      BlockWaveFunctionTruncateBasis(psi,Npsi,*epsilon);
    }
    BlockWaveFunctionFree(psi,Npsi);
    for(i=0;i<Npsi;i++)
    {
      psi[i]=Opsi[i];
    }
    // find ground-state in basis
    while(BlockLanczosGroundStateConserveBasis(H, psi, Npsi, zero))
    {
      AtLeastAnOtherLoop=1;
      printf("Run outof memory in BlockLanczosGroundStateConserveBasis, that's not good an no good solution is implemented\n");
      // was not able to evaluate ground-state in basis
      if(psi[0].NBlok==1)
      {
        printf("Run outof memory need to truncate but NBlok=1\n You're screwed\n");
        return 1;
      }
      *epsilon*=2;
      BlockWaveFunctionTruncateBasis(psi,Npsi,*epsilon);
    }
    // truncate basis
    BlockWaveFunctionTruncateBasis(psi,Npsi,*epsilon);
    // you're converged when the variance of all states is small or when the energy gain is close to zero
    Qconverged=1;
    for(i=0;i<Npsi;i++)
    {
      EnergyOld[i]=Energy[i];
      PsiHermitianOperatorPsi(H, psi[i], &(Energy[i]));
      PsiHermitianOperatorSqrPsi(H, psi[i], &(EnergySqr[i]));
      if((EnergySqr[i]-Sqr(Energy[i]))>zero) Qconverged=0;
    }
    if(Qconverged==0)
    {
      Qconverged=1;
      for(i=0;i<Npsi;i++)
      {
        if((EnergyOld[i]-Energy[i])>zero) Qconverged=0;
      }
    }
    if(AtLeastAnOtherLoop)
    {
      Qconverged=0;
      AtLeastAnOtherLoop=0;
    }
  }
  return 0;
}
