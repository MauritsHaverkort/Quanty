/*  given an Determinant and a set of Restrictions
 *  If I leave iFermion empty in the Determinant and I will be allowed to make changes to all fermions at higher index positions will I
 *  A) still be able to fullfill the restrictions - return 0
 *  B) not be able to fullfill the restrictions - return 1
 */
uint32_t QMustOccupy(uint32_t const iFermion, uint32_t const NFermion, uint8_t *det, DeterminantRestrictionType const Restrictions)
{
  uint32_t i, j, k, Nelectron, NextraOrbitals;
  uint32_t NWord=Min(NFermion,1+(NFermion-1)/8);
  uint32_t iWord=(iFermion/8);
  uint8_t detfiltered;
  
  // singles - test each restriction seperately. If one restriction tells you the orbital must be occupied return 1
  for(i=0;i<Restrictions.N;i++)
  {
    // if this restriction does not apply to iFermion do nothing for this restriction
    if(Restrictions.Determinant[i][iWord]&(128>>(iFermion&7)))
    {
      // Nelectron - calculate how many electrons are in the orbitals (0 to iFermion-1) the restrictions apply to
      Nelectron=0;
      for(j=0;j<iWord;j++)
      {
        detfiltered=Restrictions.Determinant[i][j]&det[j];
        for(k=0;k<=7;k++)
        {
          Nelectron+=(detfiltered>>k)&1;
        }
      }
      detfiltered=Restrictions.Determinant[i][iWord]&det[iWord];
      for(k=0;k<(iFermion&7);k++)
      {
        Nelectron+=(detfiltered>>(7-k))&1;
      }
      // NextraOrbitals - calculate how many orbitals there are left the restrictions apply to (iFermion to NFermion-1)
      NextraOrbitals=0;
      detfiltered=Restrictions.Determinant[i][iWord];
      for(k=1+(iFermion&7);k<=7;k++)
      {
        NextraOrbitals+=(detfiltered>>(7-k))&1;
      }
      for(j=iWord+1;j+1<NWord;j++)
      {
        detfiltered=Restrictions.Determinant[i][j];
        for(k=0;k<=7;k++)
        {
          NextraOrbitals+=(detfiltered>>k)&1;
        }
      }
      if(iWord!=NWord-1)
      {
        detfiltered=Restrictions.Determinant[i][NWord-1];
        for(k=0;k<=((NFermion-1)&7);k++)
        {
          NextraOrbitals+=(detfiltered>>(7-k))&1;
        }
      }
      // test restrictions
      if((NextraOrbitals+Nelectron)<Restrictions.Min[i]) return 1;
    }
  }
  return 0;
}

/*  given an Determinant and a set of Restrictions
 *  If I add an electron to iFermion in the Determinant and I will be allowed to make changes to all fermions at higher index positions will I
 *  A) still be able to fullfill the restrictions - return 0
 *  B) not be able to fullfill the restrictions - return 1
 */
uint32_t QCanNotOccupy(uint32_t const iFermion, uint32_t const NFermion, uint8_t *det, DeterminantRestrictionType const Restrictions)
{
  uint32_t i, j, k, Nelectron;
  uint32_t iWord=(iFermion/8);
  uint8_t detfiltered;
  // singles - test each restriction seperately. If one restriction tells you the orbital must be occupied return 1
  for(i=0;i<Restrictions.N;i++)
  {
    // if this restriction does not apply to iFermion do nothing for this Restriction
    if(Restrictions.Determinant[i][iWord]&(128>>(iFermion&7)))
    {
      // Nelectron - calculate how many electrons are in the orbitals (0 to iFermion-1) the restrictions apply to
      Nelectron=0;
      for(j=0;j<iWord;j++)
      {
        detfiltered=Restrictions.Determinant[i][j]&det[j];
        for(k=0;k<=7;k++)
        {
          Nelectron+=(detfiltered>>k)&1;
        }
      }
      detfiltered=Restrictions.Determinant[i][iWord]&det[iWord];
      for(k=0;k<(iFermion&7);k++)
      {
        Nelectron+=(detfiltered>>(7-k))&1;
      }
      // test restrictions
      if((Nelectron)>=Restrictions.Max[i]) return 1;
    }
  }
  return 0;
}

uint32_t ExtendBasisUsingRestrictions(DeterminantRestrictionType const StartingPoint, WaveFunctionType * psi)
{
  uint16_t oadd[2];
  uint32_t iN, iFermion, jFermion;
  OperatorType O;
  WaveFunctionType psi0;
  psi0 = *psi;
  InitOperator(&O);
  O.NFermion=psi0.NFermion;
  O.NBoson=psi0.NBoson;
  O.NWord=psi0.NWord;
  InitOperatorFixedLength(&O, 2, 8);
  RealOperatorAddLader(&O, oadd, 0, 1);
  for(iN=0;iN<StartingPoint.N;iN++)
  {
    for(iFermion=0;iFermion<O.NFermion;iFermion++)
    {
      for(jFermion=iFermion+1;jFermion<O.NFermion;jFermion++)
      {
        if((StartingPoint.Determinant[iN][(iFermion)>>3]&(128>>(iFermion&7)))&&(StartingPoint.Determinant[iN][(jFermion)>>3]&(128>>(jFermion&7))))
        {
          oadd[0]=CrAnMask+(uint16_t)iFermion;
          oadd[1]=(uint16_t)jFermion;
          RealOperatorAddLader(&O, oadd, 2, 1);
          oadd[0]=CrAnMask+(uint16_t)jFermion;
          oadd[1]=(uint16_t)iFermion;
          RealOperatorAddLader(&O, oadd, 2, 1);
        }
      }
    }
  }
  if(OperatorPsiRestricted(O, psi0, psi, StartingPoint))
  {
    printf("Warning: OperatorPsiRestricted failed in ExtendBasisUsingRestrictions\n");
    *psi = psi0;
    return 1;
  }
  WaveFunctionFree(&psi0);
  FreeOperator(&O);
  return 0;
}


/*  given a set of restrictions create Npsi orthonormal functions on a reasonably large basis that can serve as a starting point for
 *  finding the ground-state of H
 */
uint32_t CreateStartingFunctions(const OperatorType H, uint32_t const Npsi, DeterminantRestrictionType const StartingPoint, WaveFunctionType * psi, uint32_t const NBitsKey)
{
  double Norm, occiFermion, avvOcc, totOcc, dice;
  WaveFunctionType psi0;
  uint8_t *Determinant;
  uint32_t i,j,iFermion,MustOccupy,CannotOccupy, idet;
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
  // "some" is equal to 200+3*Npsi
  Determinant = calloc(psi0.NWord, sizeof(uint8_t));
  avvOcc = 0.5;
  for(i=0;i<Min(300+100*Npsi,5000);i++)
  {
    totOcc=0;
    for(iFermion=0;iFermion<H.NFermion;iFermion++)
    {
      MustOccupy=QMustOccupy(iFermion, H.NFermion, Determinant, StartingPoint);
      CannotOccupy=QCanNotOccupy(iFermion, H.NFermion, Determinant, StartingPoint);
      if(MustOccupy&&CannotOccupy)
      {
        printf("Cannot fullfill all restrictions in CreateStartingFunctions\n");
        free(Determinant);
        return 1;
      }
      if(MustOccupy)
      {
        OperatorCr(Determinant, (uint16_t)iFermion, psi0.NFermion);
        totOcc++;
      }
      else if(CannotOccupy)
      {
        // do nothing
      }
      else
      {
        // if current spinorbital should be filled or not is random
        // in order to make sure all states get equal changes the change if an orbital gets filled depends on the
        // occupation of that orbital in psi0
        // calculate occupation of orbital iFermion in psi0
        occiFermion = 0;
        for (idet = 1; idet <= psi0.N; idet++)
        {
          if (psi0.Determinant[BlokElementDet(idet, psi0.NWord, (iFermion >> 3))] & (128 >> (iFermion & 7)))
          {
            occiFermion += Abs(psi0.ValueRe[BlokElement(idet)]);
          }
        }
        // create random value between 0 and 1
        dice=rand()/(double)RAND_MAX;
        // if occupation < average occupation times number of electrons in the system, add an electron, otherwise not
        // add randomness in this by multiplying with (2*dice)^2
        // due to restrictions the avverage occupation might be off we thus run the avverage occupation from 0 to 200%
        if(occiFermion < Sqr(2 * dice) * avvOcc * Sqr(i+1)/200.)
        {
          OperatorCr(Determinant, (uint16_t)iFermion, psi0.NFermion);
          totOcc++;
        }
      }
    }
    WaveFunctionAddElement(&psi0, Determinant, 1, 0);
    for(j=0;j<psi0.NWord;j++)
    {
      Determinant[j]=0;
    }
    avvOcc*=i;
    avvOcc+=(totOcc/(double)H.NFermion);
    avvOcc/=(i+1.0);
  }
#ifdef PrintDebugInfo
  VerbosityPrintf(0x20000,"Random state psi0 created in CreateStartingFunctions\n");
  if(GlobalVerbosity&0x20000)
  {
    PrintWaveFunctionReadable(psi0);
  }
#endif
  // randomize psi0
  WaveFunctionRandomize(&(psi0));
  // add states by making psi = O psi0
  if(OperatorPsiRestricted(H, psi0, &(psi[0]), StartingPoint))
  {
    printf("OperatorPsiRestricted failed in CreateStartingFunctions\n");
    return 1;
  }
  WaveFunctionAddUnionBasis(&(psi[0]), psi0);
  WaveFunctionFree(&psi0);
  if(psi[0].N<Npsi)
  {
    // did not find Npsi independent determinants, extand basis further
    ExtendBasisUsingRestrictions(StartingPoint, &(psi[0]));
  }
  if(psi[0].N<Npsi)
  {
    VerbosityPrintf(0x0002, "Warning: Number of states to find smaller than number of determinants in random starting vector\n");
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
  for(i=0;i<Npsi;i++)
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
  free(Determinant);
  return 0;
}

/*  Calculate Opsi[i] = H psi[i] for i=0 to Npsi
 *  on input all psi[i]'s are on the same basis
 *  on excit all Opsi[i]'s are on the same basis as psi[i]'s
 *  memory for values of Opsi is already allocated.
 */
uint32_t BlockOperatorPsiConserveBasisPreAllocated(const OperatorType H, WaveFunctionType const * const psi, WaveFunctionType * const Opsi, uint32_t const Npsi)
{
  uint32_t i;
  for(i=0;i<Npsi;i++)
  {
    OperatorPsiConserveBasisPreAllocated(&H, &(psi[i]), &(Opsi[i]));
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
  VerbosityPrintf(0x0002,"Start of BlockOperatorPsiRestricted\n");
  if(OperatorPsiRestricted(H, psi[0], &(Opsi[0]), Restrictions))
  {
    printf("OperatorPsi failed in BlockOperatorPsi\n");
    return 1;
  }
  WaveFunctionAddMultipliedUnionBasis(&(Opsi[0]), psi[0], 0);
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
  VerbosityPrintf(0x0002,"Start of BlockOperatorPsi\n");
  if(OperatorPsi(H, psi[0], &(Opsi[0])))
  {
    printf("OperatorPsi failed in BlockOperatorPsi\n");
    return 1;
  }
  WaveFunctionAddMultipliedUnionBasis(&(Opsi[0]), psi[0], 0);
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

uint32_t BlockWaveFunctionCalculatOverlapMatrix(WaveFunctionType *const psi, uint32_t const Npsi, CompactMatrixType *const S)
{
  //TODO: test after converging matrix format
  uint32_t i,j;
  S->Ni=Npsi;
  S->Nj=Npsi;
  S->QComplex=psi[0].QComplex;
  InitCompactMatrix(S);
  if(psi[0].QComplex)
  {
    for(i=0;i<Npsi;i++)
    {
      for(j=0;j<Npsi;j++)
      {
        ComplexWaveFunctionDotProductSameBasis(psi[i], psi[j], &(S->V[j+i*Npsi].Re), &(S->V[j+i*Npsi].Im));      }
    }
  }
  else
  {
    for(i=0;i<Npsi;i++)
    {
      for(j=0;j<Npsi;j++)
      {
        S->VRe[j+i*Npsi]=RealWaveFunctionDotProductSameBasis(psi[i], psi[j]);
      }
    }
  }
  return 0;
}

/*  use a dense mehtod to find the Npsi lowest states of Operator H
 *  The basis is given by the psi[0].N determinants
 *  psi[0] to psi[Npsi-1] are allocated and live on the same basis
 *  on output psi[i].V is changed to represent the Npsi lowest eigen states
 */
uint32_t BlockLanczosGroundStateConserveBasisDense(const OperatorType H, WaveFunctionType * const psi, uint32_t const Npsi)
{
  int iOMP;
  uint32_t i,j,k;
  CompactMatrixType Hdens,fun;
  double val[psi[0].N];
  VerbosityPrintf(0x0100,"Start of BlockLanczosGroundStateConserveBasisDense\n");
  // test if psi is complex when H is complex, if not turn it into a complex function
  if((H.QComplex)&&(psi[0].QComplex==0))
  {
    VerbosityPrintf(0x0100,"Complex Hamiltonian on real wave-function, make functions complex\n");
    for(i=0;i<Npsi;i++)
    {
      if(MakeComplexWaveFunctionFromReal(&(psi[i]), psi[i]))
      {
        for(j=0;j<i;j++)
        {
          psi[j].QComplex=0;
          for(k=0;k<psi[j].NBlok;k++)
          {
            free(psi[j].ValueIm[k]);
          }
          free(psi[j].ValueIm);
        }
        printf("MakeComplexWaveFunctionFromReal failed in BlockLanczosGroundStateConserveBasisDense\n");
        return 1;
      }
    }
  }
  // create Hdens
  if(GlobalTimeOn&0x40) BeginTiming("Create Hdens", &GlobalTime);
  Hdens.QComplex = psi[0].QComplex;
  Hdens.Ni=psi[0].N;
  Hdens.Nj=psi[0].N;
  if(InitCompactMatrix(&Hdens))
  {
    printf("InitCompactMatrix failed in BlockLanczosGroundStateConserveBasisDense\n");
    return 1;
  }
  // fill the matrix Hdens
#pragma omp parallel private(i,j)
  {
    WaveFunctionType Hpsi;
    if(WaveFunctionInitCopyBasis(&Hpsi, psi[0]))
    {
      printf("WaveFunctionInitCopyBasis failed in BlockLanczosGroundStateConserveBasisDense\n");
      // need to clean things in parallel session before return - to do
      //return 1;
    }
    if(Hdens.QComplex)
    {
#pragma omp for
      for(iOMP=0;iOMP<psi[0].N;iOMP++)
      {
        j=(uint32_t)iOMP;
        OperatorDeterminantPreAllocated(&H, &(psi[0].Determinant[BlokElementDet(j+1,psi[0].NWord,0)]), &Hpsi);
        for(i=0;i<psi[0].N;i++)
        {
          Hdens.V[j+i*Hdens.Nj].Re=Hpsi.ValueRe[BlokElement(i+1)];
          Hdens.V[j+i*Hdens.Nj].Im=Hpsi.ValueIm[BlokElement(i+1)];
        }
      }
    }
    else
    {
#pragma omp for
      for(iOMP=0;iOMP<psi[0].N;iOMP++)
      {
        j=(uint32_t)iOMP;
        OperatorDeterminantPreAllocated(&H, &(psi[0].Determinant[BlokElementDet(j+1,psi[0].NWord,0)]), &Hpsi);
        for(i=0;i<psi[0].N;i++)
        {
          Hdens.VRe[j+i*Hdens.Nj]=Hpsi.ValueRe[BlokElement(i+1)];
        }
      }
    }
    WaveFunctionFreeValues(&Hpsi);
  }
  if(GlobalVerbosity&0x1000000)
  {
    printf("matrix to diagonalize\n");
    PrintCompactMatrixMinimal(Hdens);
  }
  if(GlobalTimeOn&0x40) EndTiming("Create Hdens", &GlobalTime);
  if(GlobalTimeOn&0x40) BeginTiming("Diagonalize Hdens", &GlobalTime);
  // diagonalize dense matrix
  CompactMatrixEigensystem(&Hdens, &fun, val);
  FreeCompactMatrix(&Hdens);
  if(fun.QComplex)
  {
    for(i=0;i<Min(Npsi,psi[0].N);i++)
    {
      for(j=0;j<fun.Nj;j++)
      {
        psi[i].ValueRe[BlokElement(j+1)]=fun.V[j+i*Hdens.Nj].Re;
        psi[i].ValueIm[BlokElement(j+1)]=fun.V[j+i*Hdens.Nj].Im;
      }
    }
    for(i=psi[0].N;i<Npsi;i++)
    {
      for(j=0;j<fun.Nj;j++)
      {
        psi[i].ValueRe[BlokElement(j+1)]=0;
        psi[i].ValueIm[BlokElement(j+1)]=0;
      }
    }
  }
  else
  {
    for(i=0;i<Min(Npsi,psi[0].N);i++)
    {
      for(j=0;j<fun.Nj;j++)
      {
        psi[i].ValueRe[BlokElement(j+1)]=fun.VRe[j+i*Hdens.Nj];
      }
    }
    for(i=psi[0].N;i<Npsi;i++)
    {
      for(j=0;j<fun.Nj;j++)
      {
        psi[i].ValueRe[BlokElement(j+1)]=0;
      }
    }
  }
  FreeCompactMatrix(&fun);
  if(GlobalTimeOn&0x40) EndTiming("Diagonalize Hdens", &GlobalTime);
  for(i=0;i<psi[0].N;i++)
  {
    VerbosityPrintf(0x0100,"val[%3i] = %22.15E\n",i,val[i]);
  }
  return 0;
}

/*  orthonormalize a set of Block Wave functions
 *
 */
uint32_t ListOfBlockWaveFunctionsOrthonormalizeSameBasis(WaveFunctionType **const v, uint32_t const Npsi, uint32_t const NBlock)
{
  uint32_t iBlock, ipsi, jBlock, jpsi;
  double Norm;
  double DotRe, DotIm;
  const double MinNorm=sqrt(DBL_EPSILON);
  for(iBlock=0;iBlock<NBlock;iBlock++)
  {
    for(ipsi=0;ipsi<Npsi;ipsi++)
    {
      for(jBlock=0;jBlock<iBlock;jBlock++)
      {
        for(jpsi=0;jpsi<Npsi;jpsi++)
        {
          WaveFunctionDotProductSameBasis(v[jBlock][jpsi], v[iBlock][ipsi], &DotRe, &DotIm);
          WaveFunctionAddComplexMultipliedConserveBasis(&(v[iBlock][ipsi]), v[jBlock][jpsi], -DotRe, -DotIm);
        }
      }
      for(jpsi=0;jpsi<ipsi;jpsi++)
      {
        WaveFunctionDotProductSameBasis(v[iBlock][jpsi], v[iBlock][ipsi], &DotRe, &DotIm);
        WaveFunctionAddComplexMultipliedConserveBasis(&(v[iBlock][ipsi]), v[iBlock][jpsi], -DotRe, -DotIm);
      }
      WaveFunctionCalculateNorm(v[iBlock][ipsi], &Norm);
      if(Norm<MinNorm)
      {
        Norm=0;
      }
      else
      {
        Norm=1.0/Norm;
      }
      WaveFunctionMultiply(&(v[iBlock][ipsi]), Norm, 0);
    }
  }
  return 0;
}
uint32_t ListOfBlockWaveFunctionsOrthonormalizeLastBlockSameBasis(WaveFunctionType **const v, uint32_t const Npsi, uint32_t const NBlock)
{
  uint32_t iBlock=NBlock-1, ipsi, jBlock, jpsi;
  double Norm;
  double DotRe, DotIm;
  const double MinNorm=sqrt(DBL_EPSILON);
  for(ipsi=0;ipsi<Npsi;ipsi++)
  {
    for(jBlock=0;jBlock<iBlock;jBlock++)
    {
      for(jpsi=0;jpsi<Npsi;jpsi++)
      {
        WaveFunctionDotProductSameBasis(v[jBlock][jpsi], v[iBlock][ipsi], &DotRe, &DotIm);
        WaveFunctionAddComplexMultipliedConserveBasis(&(v[iBlock][ipsi]), v[jBlock][jpsi], -DotRe, -DotIm);
      }
    }
    for(jpsi=0;jpsi<ipsi;jpsi++)
    {
      WaveFunctionDotProductSameBasis(v[iBlock][jpsi], v[iBlock][ipsi], &DotRe, &DotIm);
      WaveFunctionAddComplexMultipliedConserveBasis(&(v[iBlock][ipsi]), v[iBlock][jpsi], -DotRe, -DotIm);
    }
    WaveFunctionCalculateNorm(v[iBlock][ipsi], &Norm);
    if(Norm<MinNorm)
    {
      Norm=0;
    }
    else
    {
      Norm=1.0/Norm;
    }
    WaveFunctionMultiply(&(v[iBlock][ipsi]), Norm, 0);
  }
  return 0;
}
uint32_t BlockWaveFunctionOrthonormalizeSameBasis(WaveFunctionType *const v, uint32_t const Npsi)
{
  uint32_t ipsi, jpsi;
  double Norm;
  double DotRe, DotIm;
  const double MinNorm=sqrt(DBL_EPSILON);
  for(ipsi=0;ipsi<Npsi;ipsi++)
  {
    for(jpsi=0;jpsi<ipsi;jpsi++)
    {
      WaveFunctionDotProductSameBasis(v[jpsi], v[ipsi], &DotRe, &DotIm);
      WaveFunctionAddComplexMultipliedConserveBasis(&(v[ipsi]), v[jpsi], -DotRe, -DotIm);
    }
    WaveFunctionCalculateNorm(v[ipsi], &Norm);
    if(Norm<MinNorm)
    {
      Norm=0;
    }
    else
    {
      Norm=1.0/Norm;
    }
    WaveFunctionMultiply(&(v[ipsi]), Norm, 0);
  }
  return 0;
}

/*  Use Band Lanczos method to find lowest Npsi eigenstates of H
 *  Starting vectors are in v[0] to v[Npis-1]
 *  v[0] to v[2*Npsi+NKrylovStart-1] are allocated
 *  after a restart increase KrylovBasis by NKrylovStep if one is able to allocate the needed v's for this.
 *  maximum allocate v[2*Npsi+NKrylovMax]
 *  all states on same basis
 *
 *  on excess free memory of v[i], except v[0] to v[Npsi-1] which contains the Npsi eigen states
 */
uint32_t BlockLanczosGroundStateConserveBasisKrylovReal(const OperatorType H, WaveFunctionType * const v, uint32_t const Npsi, ConvergenceSettingsType * Settings)
{
  VerbosityPrintf(0x0100,"In BlockLanczosGroundStateConserveBasisKrylovReal\n");
  //CompactMatrixType tempS; // needed for testing of norm during convergence loops
  int    info;
  const char *jobz = "V";
  const char *uplo = "U";
  int    lwork = 1 + 6 * (Npsi+Settings->NKrylovMax) + 2 * Sqr((Npsi+Settings->NKrylovMax));
  int    liwork = 3 + 5 * (Npsi+Settings->NKrylovMax);
  int    *iwork;
  double *work;
  double *val;
  double Norm, DotRe;
  uint32_t NKrylovNow=Npsi+Settings->NKrylovStart;
  uint32_t i,j, NpsiDeflated, NorthoExtra, Deflated[Npsi],irestart=0;
  CompactMatrixType HKrylov;
  uint32_t QConverged=0;
  double **temp;
  // init several matrices
  val=malloc((Npsi+Settings->NKrylovMax)*sizeof(double));
  iwork = malloc(liwork*sizeof(int));
  if(iwork==NULL)
  {
    printf("malloc for iwork failed in BlockLanczosGroundStateConserveBasisKrylovReal\n");
    free(val);
    return 1;
  }
  work = malloc(lwork*sizeof(double));
  if(work==NULL)
  {
    printf("malloc for work failed in BlockLanczosGroundStateConserveBasisKrylovReal\n");
    free(val);
    free(iwork);
    return 1;
  }
  HKrylov.Ni=Npsi+Settings->NKrylovMax;
  HKrylov.Nj=Npsi+Settings->NKrylovMax;
  HKrylov.QComplex=0;
  if(InitCompactMatrix(&HKrylov))
  {
    printf("InitCompactMatrix failed in BlockLanczosGroundStateConserveBasisKrylov\n");
    free(val);
    free(iwork);
    free(work);
    return 1;
  }
  // start timing
  if(GlobalTimeOn&0x40) BeginTiming("Restarted Block Lanczos Conserve Basis", &GlobalTime);
  // calculate ground-state using band lanczos until converged
  while(!QConverged)
  {
    irestart++;
    // initialize variables
    NpsiDeflated=Npsi;
    NorthoExtra=0;
    memset(HKrylov.VRe, 0, HKrylov.Ni*HKrylov.Nj*sizeof(double));
    // loop over NKrylovNow states
    for(i=0;i<NKrylovNow;i++)
    {
      if(i==Npsi+1)VerbosityPrintf(0x0002,"  Restart loop %i with a Krylov basis of %i and a full basis of %i\n",irestart,NKrylovNow,v[0].N);
      // normalize v[i] and set <v[i] | H | v[i-NpsiDeflated]>
      Norm=RealWaveFunctionCalculateNorm(v[i]);
      //printf("i=%i Norm is %22.15E\n",i,Norm);
      while(Sqr(Norm)<(DBL_EPSILON+0.1*Abs(Settings->zero))*((i<NpsiDeflated)?1:(Sqr(DBL_EPSILON)+Sqr(HKrylov.VRe[0])+Sqr(HKrylov.VRe[i-NpsiDeflated+(i-NpsiDeflated)*HKrylov.Nj]))))
      {
        VerbosityPrintf(0x0100,"      Deflating from %i to %i at level %i\n",NpsiDeflated,NpsiDeflated-1,i);
        // Need to deflate
        if(i>=NpsiDeflated)
        {
          Deflated[NorthoExtra]=i-NpsiDeflated;
          NorthoExtra++;
        }
        NpsiDeflated--;
        if(NpsiDeflated==0)
        {
          // basis complete at i basis functions
          // set NKrylovNow to i and free aditional states
          QConverged=1;
          for(j=i+Npsi;j<NKrylovNow+Npsi;j++)
          {
            RealWaveFunctionFreeValues(&v[j]);
          }
          NKrylovNow=i;
          break;
        }
        temp=v[i].ValueRe;
        for(j=i;j<i+NpsiDeflated;j++)
        {
          v[j].ValueRe=v[j+1].ValueRe;
        }
        v[i+NpsiDeflated].ValueRe=temp;
        Norm=RealWaveFunctionCalculateNorm(v[i]);
        //printf("i=%i Norm is %22.15E\n",i,Norm);
      }
      if(i>=NpsiDeflated) HKrylov.VRe[i+(i-NpsiDeflated)*HKrylov.Nj]=Norm;
      RealWaveFunctionMultiply(v[i], 1/Norm);
      // orthogonalize v[j] to v[i] the NpsiDeflated candiates agains the latest found state)
      // for j from i+1 to i+NpsiDeflated-1
      // set < v[i] | H | v[j-NpsiDeflated] >
      for(j=i+1;j<i+NpsiDeflated;j++)
      {
        DotRe = RealWaveFunctionDotProductSameBasis(v[i], v[j]);
        RealWaveFunctionAddMultipliedSameBasis(v[j], v[i], -DotRe);
        if(j>=NpsiDeflated) HKrylov.VRe[i+(j-NpsiDeflated)*HKrylov.Nj]=DotRe;
      }
      // calculate v[i+NpsiDeflated] = H v[i] and onsite term ( < v[i] | H | v[i] >
      OperatorPsiConserveBasisPreAllocatedRR(&H, &v[i], &v[i+NpsiDeflated]);
      HKrylov.VRe[i+i*HKrylov.Nj]=RealWaveFunctionDotProductSameBasis(v[i], v[i+NpsiDeflated]);
      // orthogonalize v[i+NpsiDeflated] to v[j] (the new candiate agains all found states)
      // for j from i-NpsiDeflated to i-1
      for(j=(i-Min(i, NpsiDeflated));j<=i;j++)
      {
        HKrylov.VRe[j+i*HKrylov.Nj]=HKrylov.VRe[i+j*HKrylov.Nj];
        RealWaveFunctionAddMultipliedSameBasis(v[i+NpsiDeflated], v[j], -HKrylov.VRe[j+i*HKrylov.Nj]);
      }
      /* // print aditional information on convergence and norms of v
      printf("i+NpsiDeflated = %i + %i\n",i,NpsiDeflated);
      BlockWaveFunctionCalculatOverlapMatrix(v, NpsiDeflated+Npsi, &tempS);
      PrintCompactMatrixMinimal(tempS);
      FreeCompactMatrix(&tempS);
      */
      // aditional orthogonalization due to deflation of v[i+NpsiDeflated]
      for(j=0;j<NorthoExtra;j++)
      {
        HKrylov.VRe[Deflated[j]+i*HKrylov.Nj]=RealWaveFunctionDotProduct(&v[Deflated[j]], &v[i+NpsiDeflated]);
        HKrylov.VRe[i+Deflated[j]*HKrylov.Nj]=HKrylov.VRe[Deflated[j]+i*HKrylov.Nj];
        RealWaveFunctionAddMultipliedSameBasis(v[i+NpsiDeflated], v[Deflated[j]], -HKrylov.VRe[Deflated[j]+i*HKrylov.Nj]);
      }
      // aditional reorthogonalization
      for(j=0;j<=i;j++)
      {
        DotRe = RealWaveFunctionDotProductSameBasis(v[j], v[i+NpsiDeflated]);
        RealWaveFunctionAddMultipliedSameBasis(v[i+NpsiDeflated], v[j], -DotRe);
      }
    }
    //Calculate Eigen system of HKrylov
    dsyevd(jobz, uplo, &NKrylovNow, HKrylov.VRe, &HKrylov.Ni, val, work, &lwork, iwork, &liwork, &info);
    if(info)
    {
      printf("\n\n\nWARNING\n\n");
      printf("dsyevd failed with info %i in BlockLanczosGroundStateConserveBasisKrylovReal\n",info);
      printf("Hopefully the next restart is better\n\n\n");
    }
    VerbosityPrintf(0x0100,"  The eigen-values of the current lowest %i states are:\n",Npsi);
    for(i=0;i<Npsi;i++)
    {
      VerbosityPrintf(0x0100,"    Eigenvalue %3i = %22.15E\n",i,val[i]);
    }
    fflush(stdout);
    //Calculate new v[0] to v[Npsi-1]
    for(i=NKrylovNow;i<NKrylovNow+Npsi;i++)
    {
      RealWaveFunctionMultiply(v[i], 0);
    }
    for(i=0;i<Npsi;i++)
    {
      for(j=0;j<NKrylovNow;j++)
      {
        RealWaveFunctionAddMultipliedSameBasis(v[i+NKrylovNow], v[j], HKrylov.VRe[j+i*HKrylov.Nj]);
      }
    }
    for(i=0;i<Npsi;i++)
    {
      temp=v[i].ValueRe;
      v[i].ValueRe=v[i+NKrylovNow].ValueRe;
      v[i+NKrylovNow].ValueRe=temp;
    }
    // Expand Krylov basis if possible by memory constraints
    if(NKrylovNow+Settings->NKrylovStep+Npsi<Settings->NKrylovMax)
    {
      for(i=NKrylovNow+Npsi;i<NKrylovNow+2*Npsi+Settings->NKrylovStep;i++)
      {
        if(WaveFunctionInitCopyBasis(&v[i], v[0]))
        {
          for(j=NKrylovNow+Npsi;j<i;j++)
          {
            RealWaveFunctionFreeValues(&v[j]);
          }
          goto OutOfMemNoKrylovIncrease;
        }
      }
      NKrylovNow+=Settings->NKrylovStep+Npsi;
    OutOfMemNoKrylovIncrease:;
    }
    VerbosityPrintf(0x0100,"\n");
  }
  if(GlobalTimeOn&0x40) EndTiming("Restarted Block Lanczos Conserve Basis", &GlobalTime);
  free(val);
  free(iwork);
  free(work);
  FreeCompactMatrix(&HKrylov);
  for(i=Npsi;i<Npsi+NKrylovNow;i++)
  {
    WaveFunctionFreeValues(&v[i]);
  }
  return 0;
}

uint32_t BlockLanczosGroundStateConserveBasisKrylovComplex(const OperatorType H, WaveFunctionType * const v, uint32_t const Npsi, ConvergenceSettingsType * Settings)
{
  VerbosityPrintf(0x0100,"In BlockLanczosGroundStateConserveBasisKrylovComplex\n");
  //CompactMatrixType tempS; // needed for testing of norm during convergence loops
  int    info;
  const char *jobz = "V";
  const char *uplo = "U";
  int    lwork  = 1 + 6 * (Npsi+Settings->NKrylovMax) + 2 * Sqr((Npsi+Settings->NKrylovMax));
  int    lcwork = 2 * (Npsi+Settings->NKrylovMax) + Sqr((Npsi+Settings->NKrylovMax));
  int    liwork = 3 + 5 * (Npsi+Settings->NKrylovMax);
  int    *iwork;
  double *work;
  Complex16 *cwork;
  double *val;
  double Norm, DotRe, DotIm;
  uint32_t NKrylovNow=Npsi+Settings->NKrylovStart;
  uint32_t i,j, NpsiDeflated, NorthoExtra, Deflated[Npsi],irestart=0;
  CompactMatrixType HKrylov;
  uint32_t QConverged=0;
  double **tempRe, **tempIm;
  // init several matrices
  val=malloc((Npsi+Settings->NKrylovMax)*sizeof(double));
  iwork = malloc(liwork*sizeof(int));
  if(iwork==NULL)
  {
    printf("malloc for iwork failed in BlockLanczosGroundStateConserveBasisKrylovComplex\n");
    free(val);
    return 1;
  }
  work = malloc(lwork*sizeof(double));
  if(work==NULL)
  {
    printf("malloc for work failed in BlockLanczosGroundStateConserveBasisKrylovComplex\n");
    free(val);
    free(iwork);
    return 1;
  }
  cwork = malloc(lcwork*sizeof(Complex16));
  if(cwork==NULL)
  {
    printf("malloc for cwork failed in BlockLanczosGroundStateConserveBasisKrylovComplex\n");
    free(val);
    free(iwork);
    free(work);
    return 1;
  }
  HKrylov.Ni=Npsi+Settings->NKrylovMax;
  HKrylov.Nj=Npsi+Settings->NKrylovMax;
  HKrylov.QComplex=1;
  if(InitCompactMatrix(&HKrylov))
  {
    printf("InitCompactMatrix failed in BlockLanczosGroundStateConserveBasisKrylovComplex\n");
    free(val);
    free(iwork);
    free(work);
    return 1;
  }
  // start timing
  if(GlobalTimeOn&0x40) BeginTiming("Restarted Block Lanczos Conserve Basis Complex", &GlobalTime);
  // calculate ground-state using band lanczos until converged
  while(!QConverged)
  {
    irestart++;
    // initialize variables
    NpsiDeflated=Npsi;
    NorthoExtra=0;
    memset(HKrylov.V, 0, HKrylov.Ni*HKrylov.Nj*sizeof(Complex16));
    // loop over NKrylovNow states
    for(i=0;i<NKrylovNow;i++)
    {
      if(i==Npsi+1)VerbosityPrintf(0x0002,"  Restart loop %i with a Krylov basis of %i and a full basis of %i\n",irestart,NKrylovNow,v[0].N);
      // normalize v[i] and set <v[i] | H | v[i-NpsiDeflated]>
      Norm=ComplexWaveFunctionCalculateNorm(v[i]);
      //printf("i=%i Norm is %22.15E\n",i,Norm);
      while(Sqr(Norm)<(DBL_EPSILON+0.1*Abs(Settings->zero))*((i<NpsiDeflated)?1:(Sqr(DBL_EPSILON)+Sqr(HKrylov.V[0].Re)+Sqr(HKrylov.V[i-NpsiDeflated+(i-NpsiDeflated)*HKrylov.Nj].Re))))
      {
        VerbosityPrintf(0x0100,"      Deflating from %i to %i at level %i\n",NpsiDeflated,NpsiDeflated-1,i);
        // Need to deflate
        if(i>=NpsiDeflated)
        {
          Deflated[NorthoExtra]=i-NpsiDeflated;
          NorthoExtra++;
        }
        NpsiDeflated--;
        if(NpsiDeflated==0)
        {
          // basis complete at i basis functions
          // set NKrylovNow to i and free aditional states
          QConverged=1;
          for(j=i+Npsi;j<NKrylovNow+Npsi;j++)
          {
            ComplexWaveFunctionFreeValues(&v[j]);
          }
          NKrylovNow=i;
          break;
        }
        tempRe=v[i].ValueRe;
        tempIm=v[i].ValueIm;
        for(j=i;j<i+NpsiDeflated;j++)
        {
          v[j].ValueRe=v[j+1].ValueRe;
          v[j].ValueIm=v[j+1].ValueIm;
        }
        v[i+NpsiDeflated].ValueRe=tempRe;
        v[i+NpsiDeflated].ValueIm=tempIm;
        Norm=ComplexWaveFunctionCalculateNorm(v[i]);
        //printf("i=%i Norm is %22.15E\n",i,Norm);
      }
      if(i>=NpsiDeflated)
      {
        HKrylov.V[i+(i-NpsiDeflated)*HKrylov.Nj].Re=Norm;
        HKrylov.V[i+(i-NpsiDeflated)*HKrylov.Nj].Im=0;
      }
      ComplexWaveFunctionMultiply(v[i], 1/Norm);
      // orthogonalize the NpsiDeflated candiates agains the latest found state
      // orthogonalize v[j] to v[i]
      // for j from i+1 to i+NpsiDeflated-1
      // set < v[i] | H | v[j-NpsiDeflated] >
      for(j=i+1;j<i+NpsiDeflated;j++)
      {
        ComplexWaveFunctionDotProductSameBasis(v[i], v[j], &DotRe, &DotIm);
        ComplexWaveFunctionAddComplexMultipliedConserveBasis(&v[j], v[i], -DotRe, -DotIm);
        if(j>=NpsiDeflated)
        {
          HKrylov.V[i+(j-NpsiDeflated)*HKrylov.Nj].Re= DotRe;
          HKrylov.V[i+(j-NpsiDeflated)*HKrylov.Nj].Im= DotIm;
        }
      }
      // calculate v[i+NpsiDeflated] = H v[i] and onsite term < v[i] | H | v[i] >
      if(H.QComplex==1)
      {
        OperatorPsiConserveBasisPreAllocatedCC(&H, &v[i], &v[i+NpsiDeflated]);
      }
      else
      {
        OperatorPsiConserveBasisPreAllocatedMC(&H, &v[i], &v[i+NpsiDeflated]);
      }
      ComplexWaveFunctionDotProductSameBasis(v[i], v[i+NpsiDeflated], &HKrylov.V[i+i*HKrylov.Nj].Re, &DotIm);
      // orthogonalize v[i+NpsiDeflated] to v[j] (the new candiate agains the latest found states)
      // for j from i-NpsiDeflated to i-1
      for(j=(i-Min(i, NpsiDeflated));j<=i;j++)
      {
        HKrylov.V[j+i*HKrylov.Nj].Re= HKrylov.V[i+j*HKrylov.Nj].Re;
        HKrylov.V[j+i*HKrylov.Nj].Im=-HKrylov.V[i+j*HKrylov.Nj].Im;
        ComplexWaveFunctionAddComplexMultipliedSameBasis(v[i+NpsiDeflated], v[j], -HKrylov.V[j+i*HKrylov.Nj].Re, -HKrylov.V[j+i*HKrylov.Nj].Im);
      }
      /*
       // print aditional information on convergence and norms of v
       printf("i+NpsiDeflated = %i + %i\n",i,NpsiDeflated);
       BlockWaveFunctionCalculatOverlapMatrix(v, NpsiDeflated+i, &tempS);
       CompactMatrixChop(&tempS, 10000*DBL_EPSILON);
       PrintCompactMatrixMinimal(tempS);
       FreeCompactMatrix(&tempS);
      */
      // aditional orthogonalization due to deflation of v[i+NpsiDeflated]
      for(j=0;j<NorthoExtra;j++)
      {
        ComplexWaveFunctionDotProductSameBasis(v[Deflated[j]], v[i+NpsiDeflated], &HKrylov.V[Deflated[j]+i*HKrylov.Nj].Re, &HKrylov.V[Deflated[j]+i*HKrylov.Nj].Im);
        HKrylov.V[i+Deflated[j]*HKrylov.Nj].Re= HKrylov.V[Deflated[j]+i*HKrylov.Nj].Re;
        HKrylov.V[i+Deflated[j]*HKrylov.Nj].Im=-HKrylov.V[Deflated[j]+i*HKrylov.Nj].Im;
        ComplexWaveFunctionAddComplexMultipliedSameBasis(v[i+NpsiDeflated], v[Deflated[j]], -HKrylov.V[Deflated[j]+i*HKrylov.Nj].Re, -HKrylov.V[Deflated[j]+i*HKrylov.Nj].Im);
      }
      // double orthogonalization - in principle only needed for some cases -- large spread in eigen values (10^2 to 10^-2) -- only for last loop ???
      for(j=0;j<=i;j++)
      {
        ComplexWaveFunctionDotProductSameBasis(v[j], v[i+NpsiDeflated], &DotRe, &DotIm);
        ComplexWaveFunctionAddComplexMultipliedSameBasis(v[i+NpsiDeflated], v[j], -DotRe, -DotIm);
      }
    }
    //Calculate Eigen system of HKrylov
    zheevd(jobz, uplo, &NKrylovNow, HKrylov.V, &HKrylov.Ni, val, cwork, &lcwork, work, &lwork, iwork, &liwork, &info);
    if(info)
    {
      printf("\n\n\nWARNING\n\n");
      printf("zheevd failed with info %i in BlockLanczosGroundStateConserveBasisKrylovComplex\n",info);
      printf("Hopefully the next restart is better\n\n\n");
    }
    VerbosityPrintf(0x0100,"  The eigen-values of the current lowest %i states are:\n",Npsi);
    for(i=0;i<Npsi;i++)
    {
      VerbosityPrintf(0x0100,"    Eigenvalue %3i = %22.15E\n",i,val[i]);
    }
    fflush(stdout);
    //Calculate new v[0] to v[Npsi-1]
    for(i=NKrylovNow;i<NKrylovNow+Npsi;i++)
    {
      ComplexWaveFunctionMultiply(v[i], 0);
    }
    for(i=0;i<Npsi;i++)
    {
      for(j=0;j<NKrylovNow;j++)
      {
        ComplexWaveFunctionAddComplexMultipliedConserveBasis(&v[i+NKrylovNow], v[j], HKrylov.V[j+i*HKrylov.Nj].Re, HKrylov.V[j+i*HKrylov.Nj].Im);
      }
    }
    for(i=0;i<Npsi;i++)
    {
      tempRe=v[i].ValueRe;
      tempIm=v[i].ValueIm;
      v[i].ValueRe=v[i+NKrylovNow].ValueRe;
      v[i].ValueIm=v[i+NKrylovNow].ValueIm;
      v[i+NKrylovNow].ValueRe=tempRe;
      v[i+NKrylovNow].ValueIm=tempIm;
    }
    // Expand Krylov basis if possible by memory constraints
    if(NKrylovNow+Settings->NKrylovStep+Npsi<Settings->NKrylovMax)
    {
      for(i=NKrylovNow+Npsi;i<NKrylovNow+2*Npsi+Settings->NKrylovStep;i++)
      {
        if(WaveFunctionInitCopyBasis(&v[i], v[0]))
        {
          for(j=NKrylovNow+Npsi;j<i;j++)
          {
            ComplexWaveFunctionFreeValues(&v[j]);
          }
          goto OutOfMemNoKrylovIncrease;
        }
      }
      NKrylovNow+=Settings->NKrylovStep+Npsi;
    OutOfMemNoKrylovIncrease:;
    }
    VerbosityPrintf(0x0100,"\n");
  }
  if(GlobalTimeOn&0x40) EndTiming("Restarted Block Lanczos Conserve Basis Complex", &GlobalTime);
  free(val);
  free(iwork);
  free(work);
  free(cwork);
  FreeCompactMatrix(&HKrylov);
  for(i=Npsi;i<Npsi+NKrylovNow;i++)
  {
    ComplexWaveFunctionFreeValues(&v[i]);
  }
  return 0;
}

uint32_t BlockLanczosGroundStateConserveBasisKrylov(const OperatorType H, WaveFunctionType * const v, uint32_t const Npsi, ConvergenceSettingsType * Settings)
{
  if(H.QComplex)
  {
    return BlockLanczosGroundStateConserveBasisKrylovComplex(H, v, Npsi, Settings);
  }
  else
  {
    return BlockLanczosGroundStateConserveBasisKrylovReal(H, v, Npsi, Settings);
  }
}


/*  Use Band Lanczos method to find lowest Npsi eigenstates of H
 *  Starting vectors are in v[0] to v[Npsi-1]
 *  v[0] to v[2*Npsi] are allocated
 *  all states on same basis
 *
 *  there are NpsiAlloc > 2 Npsi+1 states to work with
 *  on excess free memory of v[Npsi] to v[NpsiAlloc-1]
 */
uint32_t BlockLanczosGroundStateConserveBasisKrylovRecalculate(const OperatorType H, WaveFunctionType * const v, uint32_t const Npsi, uint32_t const NpsiAlloc, ConvergenceSettingsType * Settings)
{
  uint32_t temp;
  
  printf("to do BlockLanczosGroundStateConserveBasisKrylovRecalculate\n");
  printf("Cheap fix needs to be improved\n");
  
  temp=Settings->NKrylovStart;
  Settings->NKrylovStart=Min(Settings->NKrylovStart,NpsiAlloc-2*Npsi);
  
  if(BlockLanczosGroundStateConserveBasisKrylov(H, v, Npsi, Settings))
  {
    printf("BlockLanczosGroundStateConserveBasisKrylov failed in BlockLanczosGroundStateConserveBasis\n");
    return 1;
  }
  Settings->NKrylovStart=temp;
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
 *   if psi[0].N <= DenseBoarder do a dense calculation
 *
 */
uint32_t BlockLanczosGroundStateConserveBasis(const OperatorType H, WaveFunctionType * const psi, uint32_t const Npsi, ConvergenceSettingsType * Settings)
{
  uint32_t i,MaxAlloc=2*Npsi+Settings->NKrylovMax, StartAlloc=2*Npsi+Settings->NKrylovStart;
  WaveFunctionType *v;
  VerbosityPrintf(0x10000000,"Start of BlockLanczosGroundStateConserveBasis\n");
  // for small basises do a dense calculation
  if(psi[0].N<=Settings->DenseBoarder) return BlockLanczosGroundStateConserveBasisDense(H,psi, Npsi);
  // create the vectors v in memory.
  //   if all 2*Npsi+NKrylovMax are created return BlockLanczosGroundStateConserveBasisKrylov
  //   if only 2*Npsi+1 states can be created return BlockLanczosGroundStateConserveBasisKrylovRecalculate
  //   if less then 2*Npsi+1 states can be created return 1
  v=malloc(MaxAlloc*sizeof(WaveFunctionType));
  if(v==NULL)
  {
    printf("malloc failed in BlockLanczosGroundStateConserveBasis\n");
    return 1;
  }
  for(i=0;i<Npsi;i++)
  {
    v[i]=psi[i];
  }
  for(i=Npsi;i<StartAlloc;i++)
  {
    if(WaveFunctionInitCopyBasis(&(v[i]), (v[0])))
    {
      if(i<2*Npsi+1)
      {
        return 1;
      }
      StartAlloc=i-1;
      if(BlockLanczosGroundStateConserveBasisKrylovRecalculate(H, v, Npsi, StartAlloc, Settings))
      {
        printf("BlockLanczosGroundStateConserveBasisKrylovRecalculate failed in BlockLanczosGroundStateConserveBasis\n");
        return 1;
      }
      for(i=0;i<Npsi;i++)
      {
        psi[i]=v[i];
      }
      free(v);
      return 0;
    }
  }
  if(BlockLanczosGroundStateConserveBasisKrylov(H, v, Npsi, Settings))
  {
    printf("BlockLanczosGroundStateConserveBasisKrylov failed in BlockLanczosGroundStateConserveBasis\n");
    return 1;
  }
  for(i=0;i<Npsi;i++)
  {
    psi[i]=v[i];
  }
  free(v);
  return 0;
}

/*  Find the lowest Npsi eigen-states of the operator H
 *
 *  Input:
 *    H:             Operator
 *    Npsi:          number of states to converge
 *    Settings:      ConvergenceSettingsType* pointer to important values that need to be set:
 *                   Settings->zero    smallest possitive double that is considered different from zero when comparing convergence creteria
 *                   Settings->epsilon smallest positive double that is considered different from zero when truncating wavefunctions
 *    StartingPoint: a set of restrictions defining the starting point wave function
 *                   if StartingPoint.N = 0 it is assumed that the functions psi are intialized and one starts from there
 *    Restrictions:  a set of restrictions makeing the basis smaler
 *    psi:           a pointer to Npsi functions (If StartingPoint.N=0 they should be initialized as form a basis for the starting point of the Lanczos routine, otherwise they should be unititiallized)
 *    NBitsKey:      The number of bits in the key (only used when StartingPoint.N != 0
 *
 *  Output:
 *    psi:           a list of Npsi wave-functions representing the Npsi lowest eigen-states
 *    epsilon:       will change if not enough memory is available.
 *
 */
uint32_t BlockGroundState(const OperatorType H, uint32_t const Npsi, ConvergenceSettingsType * Settings, DeterminantRestrictionType const StartingPoint, DeterminantRestrictionType const Restrictions, WaveFunctionType * const psi, const uint32_t NBitsKey)
{
  uint32_t i,Qconverged=0, AtLeastAnOtherLoop=0,iloop=0;;
  WaveFunctionType Opsi[Npsi];
  double Energy[2][Npsi],EnergySqr[2][Npsi],temp;
  // if StartingPoint.N != 0
  // create Npsi independent functions fullfiling the starting point conditions
  // all functions must have the same basis
  // if H is complex or mixed psi will be complex, if H is real psi is real
  VerbosityPrintf(0x10000000,"Start of BlockGroundState\n");
  if((StartingPoint.N!=0)&&(!(NBitsKey==0)))
  {
    if(CreateStartingFunctions(H, Npsi, StartingPoint, psi, NBitsKey))
    {
      printf("CreateStartingFunctions failed in BlockGroundState\n");
      return 1;
    }
    BlockGroundState(H, Npsi, Settings, StartingPoint, StartingPoint, psi, 0);
  }
  VerbosityPrintf(0x20000,"\nStarting function 1 to Npsi=%i in BlockGroundState\n\n",Npsi);
  if(GlobalVerbosity&0x20000)
  {
      for(i=0;i<Npsi;i++) PrintWaveFunctionReadable(psi[i]);
  }
  // set energy[0]
  for(i=0;i<Npsi;i++)
  {
    Energy[0][i]=1E+99;
  }
  VerbosityPrintf(0x0002,"Start of BlockGroundState. Converge %i states to an energy with relative variance smaller than %22.15E\n\n",Npsi,Settings->zero);
  // loop until converged
  while((!Qconverged)&&(iloop<100))
  {
    iloop++;
    // expand basis depending on Restrictions.N use either restrictions or don't use them
    while( (Restrictions.N) ? BlockOperatorPsiRestricted(H, psi, Opsi, Npsi, Restrictions) : BlockOperatorPsi(H, psi, Opsi, Npsi))
    {
      // was not able to evaluate O psi, need to truncate psi further
      if(psi[0].NBlok==1)
      {
        printf("Run outof memory need to truncate but NBlok=1\n You're screwed\n");
        return 1;
      }
      Settings->epsilon*=2;
      BlockWaveFunctionTruncateBasis(psi,Npsi,Settings->epsilon);
    }
    // calculate Energy and variance of Energy - test for convergence
    Qconverged=1;
    VerbosityPrintf(0x0100,"Number of determinants in basis %i (and basis of H psi %i)\n",psi[0].N,Opsi[0].N);
    VerbosityPrintf(0x0100,"State # |   Energy               |  Variance of Energy    |  Energy change\n");
    for(i=0;i<Npsi;i++)
    {
      Energy[1][i]=Energy[0][i];
      EnergySqr[1][i]=EnergySqr[0][i];
      WaveFunctionDotProduct(psi[i], Opsi[i], &Energy[0][i], &temp);
      WaveFunctionDotProductSameBasis(Opsi[i], Opsi[i], &EnergySqr[0][i], &temp);
      VerbosityPrintf(0x0100,"%7i | %22.15E | %22.15E | %22.15E\n",i,Energy[0][i],EnergySqr[0][i]-Sqr(Energy[0][i]),Energy[0][i]-Energy[1][i]);
      if((EnergySqr[0][i]-Sqr(Energy[0][i])>=(EnergySqr[0][i]*Settings->zero))&&(((Energy[1][i]-Energy[0][i])>fabs(Energy[0][i]*10*DBL_EPSILON))||(Energy[0][i]-Energy[1][i]>fabs(Energy[0][i]*Settings->zero)))) Qconverged=0;
    }
    VerbosityPrintf(0x0100,"\n");
    VerbosityPrintf(0x10000000,"Converged %i\n",Qconverged);
    if(AtLeastAnOtherLoop)
    {
      Qconverged=0;
      AtLeastAnOtherLoop=0;
    }
    if(Qconverged)
    {
      // free Opsi
      BlockWaveFunctionFree(Opsi,Npsi);
    }
    else
    {
      VerbosityPrintf(0x0002,"Outer loop %3i, Number of Determinants: %9i %9i last variance %22.15E\n",iloop,psi[0].N,Opsi[0].N,EnergySqr[0][Npsi-1]-Sqr(Energy[0][Npsi-1]));
      // free psi and set psi=Opsi
      BlockWaveFunctionFree(psi,Npsi);
      for(i=0;i<Npsi;i++)
      {
        psi[i]=Opsi[i];
      }
      // find ground-state in basis
      while(BlockLanczosGroundStateConserveBasis(H, psi, Npsi, Settings))
      {
        AtLeastAnOtherLoop=1;
        printf("Run outof memory in BlockLanczosGroundStateConserveBasis, that's not good an no good solution is implemented\n");
        // was not able to evaluate ground-state in basis
        if(psi[0].NBlok==1)
        {
          printf("Run outof memory need to truncate but NBlok=1\n You're screwed\n");
          return 1;
        }
        Settings->epsilon*=2;
        BlockWaveFunctionTruncateBasis(psi,Npsi,Settings->epsilon);
      }
      // truncate basis
      BlockWaveFunctionTruncateBasis(psi,Npsi,Settings->epsilon);
    }
  }
  return 0;
}

