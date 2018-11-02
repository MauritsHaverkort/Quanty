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

uint32_t ExtendBasisCompletenessSmallSystems(DeterminantRestrictionType const StartingPoint, WaveFunctionType * psi)
{
  uint32_t i,iFull;
  uint8_t  newdet[4]={0,0,0,0};
  if(psi->NFermion>30)
  {
    printf("Error in ExtendBasisCompletenessSmallSystems max number of Fermions is 30\n");
    return 1;
  }
  iFull=0;
  for(i=0;i<psi->NFermion;i++)
  {
    iFull+=1<<(31-i);
  }  
  if(psi->QComplex)
  {
    if(QDetNotRestricted(newdet, psi->NFermion, psi->NBoson, psi->NWord, StartingPoint))
    {
      if(ComplexWaveFunctionAddElement(psi, newdet, 1,1))
      {
        return 1;
      }
    }
    for(i=iFull;i>(iFull>>(psi->NFermion));i-=1<<(32-psi->NFermion))
    {
      newdet[0]=(uint8_t)(i>>24);
      newdet[1]=(uint8_t)(i>>16);
      newdet[2]=(uint8_t)(i>>8);
      newdet[3]=(uint8_t)(i);
      if(QDetNotRestricted(newdet, psi->NFermion, psi->NBoson, psi->NWord, StartingPoint))
      {
        if(ComplexWaveFunctionAddElement(psi, newdet, 1,1))
        {
          return 1;
        }
      }
    }
  }
  else
  {
    if(QDetNotRestricted(newdet, psi->NFermion, psi->NBoson, psi->NWord, StartingPoint))
    {
      if(RealWaveFunctionAddElement(psi, newdet, 1))
      {
        return 1;
      }
    }
    for(i=iFull;i>(iFull>>(psi->NFermion));i-=1<<(32-psi->NFermion))
    {
      newdet[0]=(uint8_t)(i>>24);
      newdet[1]=(uint8_t)(i>>16);
      newdet[2]=(uint8_t)(i>>8);
      newdet[3]=(uint8_t)(i);
      if(QDetNotRestricted(newdet, psi->NFermion, psi->NBoson, psi->NWord, StartingPoint))
      {
        if(RealWaveFunctionAddElement(psi, newdet, 1))
        {
          return 1;
        }
      }
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
  if(psi0.NFermion<=24)
  {
    ExtendBasisCompletenessSmallSystems(StartingPoint, &(psi0));
#ifdef PrintDebugInfo
    VerbosityPrintf(0x20000,"All possible states psi0 created in CreateStartingFunctions\n");
    if(GlobalVerbosity&0x20000)
    {
      PrintWaveFunctionReadable(psi0);
    }
#endif
    // randomize psi0
    WaveFunctionRandomize(&(psi0));
    // set psi[0]
    psi[0] = psi0;
  }
  else
  {
    Determinant = (uint8_t*)calloc(psi0.NWord, sizeof(uint8_t));
    // strategy based on randomness and occupations
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
    free(Determinant);
  }
  // print warning if there are not enough starting determinant
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
  return 0;
}

/*  use a dense method to find the Npsi lowest states of Operator H
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


/*******************************************************************
 *********************** Block Tri-diagonalize *********************
 *******************************************************************/


/*  Calculates a Block tri diagonal representation of H starting from the wavefunctions psi
 *  On input
 *  H Hamiltonian
 *  psi[0] to psi[M->BlockSize-1] wavefunctions
 *  epsilon smallest prefactor for a determinant considered different from zero
 *  Gamma imaginary part of broadening. (used to determine convergence)
 *  M->N number of tri-diagonal steps made in the calculation. (when converged wave-functions truncate to lenght zero and everyghing is fast, M->N is there to set a maximum)
 *    Note that current implementations do not test enough on convergence and one normally does not reach a state where code thinks we are converged, (although we are)
 *
 *  Different versions for real and complex, Conserve basis and restrictions on the operator. 
 *  Specific versions are "unsafe" no checking of imput conditions is done.
 *
 *  For the noral version the wave-functions in the list psi are distroyed (and freed) after a succesfull calculation
 *  all wavefunctions must be on the same basis
 *  For the Krylov version the wave-funcitons in the list psi is expanded and contains the krylov basis after a succesfull calculation
 *  Each block is on the same basis and has length M.BlockSize
 *  There are M.N Blocks in total
 *
 */



uint32_t LanczosBlockTriDiagonalizeRestrictedRC(const OperatorType H, WaveFunctionType * const psi, DeterminantRestrictionType Restrictions, double epsilon, BlockTriDiagonalMatrixType * const M, const double Gamma)
{
  uint32_t i,j,k,vi,vip,vim;
  CompactMatrixType Mab;
  WaveFunctionType *v[4];
  uint32_t NSingular,Npsi=M->BlockSize,Npsivi=M->BlockSize,Npsivim=M->BlockSize;
  VerbosityPrintf(0x0001,"Start of LanczosBlockTriDiagonalizeRestrictedRC\n");
  fflush(stdout);
  for(i=0;i<M->BlockSize;i++)
  {
    M->b[i+i*M->BlockSize].Re=0;
  }
  v[0] = psi;
  v[1] = (WaveFunctionType*)malloc(Npsi*sizeof(WaveFunctionType));
  v[2] = (WaveFunctionType*)malloc(Npsi*sizeof(WaveFunctionType));
  v[3] = (WaveFunctionType*)malloc(Npsi*sizeof(WaveFunctionType));
  Mab.Ni=M->BlockSize;
  Mab.Nj=M->BlockSize;
  Mab.QComplex=1;
  ComplexWaveFunctionListSVDOrthogonalizeSameBasis(v[0], v[1], &Npsi, epsilon, true); // v[1] = orthonormal[v0]
  NSingular=M->BlockSize-Npsi;
  if(NSingular!=0)
  {
    VerbosityPrintf(0x0001,"Reduce dimension from %u to %u in block %u\n",M->BlockSize,M->BlockSize-NSingular,0);
  }
  // b[0] = v[0] . v[1]^+ or b[0]^T = v[1]^* . v[0]^T
  for(i=0;i<Npsi+NSingular;i++)
  {
    for(j=0;j<Npsi;j++)
    {
      ComplexWaveFunctionDotProductSameBasis(v[1][j], v[0][i], &(M->b[i * M->BlockSize + j].Re), &(M->b[i * M->BlockSize + j].Im));
    }
  }
  vi=1;
  vim=0;
  vip=2;
  Npsivi=Npsi;
  for(k=1;k<M->N;k++)
  {
    // v[vip] = H v[vi]^T
    while(BlockOperatorPsiRestrictedRC(H, v[vi], v[vip], Restrictions, Npsivi))
    {
      printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRestrictedRC for the function BlockOperatorPsiRestrictedRC\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      if (v[vi][0].NBlok == 1) AllocFailed();
      ComplexBlockWaveFunctionTruncateBasis(v[vi], Npsivi, epsilon);
    }
    // a[k] = v[vi]^* . v[vip]^T  (v[vi]^* . H . v[vi]^T)
    Mab.V = &(M->a[k*Sqr(M->BlockSize)]);
    ComplexBlockWaveFunctionDotProduct(v[vi], Npsi, v[vip], Npsi, &Mab); //Npsivi = Npsi - changes at Lowdin Orthonormalization
    if(k!=1)
    {
      // b[k-1] = v[vim]^* . v[vip]^T  (v[vim]^* . H . v[vi]^T)
      Mab.V = &(M->b[(k-1)*Sqr(M->BlockSize)]);
      ComplexBlockWaveFunctionDotProduct(v[vim], Npsivim, v[vip], Npsi, &Mab);
    }
    NSingular=1;
    while(NSingular)
    {
      // orthogonalize v[vip] to v[vi] and v[vim]
      // we use matrix multiplication between v[vip] and v[vi] which is a numerical unstable form of Gram Schmidt
      // we heal this by repeating the orthogonalization twice
      Mab.V=M->a;
      // w = v[vip]^* . v[vi]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vi], Npsivi, &Mab);
      // v[vip] -= w^* . v[vi]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vi], Npsivi, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRestrictedRC for the function ComplexBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vi],  Npsivi, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,   epsilon);
      }
      // w = v[vip]^* . v[vi]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vi], Npsivi, &Mab);
      // v[vip] -= w^* . v[vi]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vi], Npsivi, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRestrictedRC for the function ComplexBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vi],  Npsivi, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,   epsilon);
      }
      // w = v[vip]^* . v[vim]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vim], Npsivim, &Mab);
      // v[vip] -= w^* . v[vim]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vim], Npsivim, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRestrictedRC for the function ComplexBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vim], Npsivim, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,    epsilon);
      }
      // w = v[vip]^* . v[vim]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vim], Npsivim, &Mab);
      // v[vip] -= w^* . v[vim]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vim], Npsivim, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRestrictedRC for the function ComplexBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vim], Npsivim, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,    epsilon);
      }
      // w = v[vip]^* . v[vi]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vi], Npsivi, &Mab);
      // v[vip] -= w^* . v[vi]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vi], Npsivi, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRestrictedRC for the function ComplexBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vi],  Npsivi, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,   epsilon);
      }
      // w = v[vip]^* . v[vim]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vim], Npsivim, &Mab);
      // v[vip] -= w^* . v[vim]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vim], Npsivim, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRestrictedRC for the function ComplexBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vim], Npsivim, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,    epsilon);
      }
      // Lowdin orthonormalize v[vip] into v[3], free v[vip] and set v[vip] equal to v[3]
      NSingular=Npsi;
      ComplexWaveFunctionListSVDOrthogonalizeSameBasis(v[vip], v[3], &Npsi, epsilon, true);
      NSingular-=Npsi;
      if(NSingular!=0)
      {
        VerbosityPrintf(0x0001, "Reduce dimension from %u to %u in block %u\n",Npsi+NSingular,Npsi,k);
      }
      for(i=0;i<Npsi+NSingular;i++)
      {
        ComplexWaveFunctionFreeValues(&v[vip][i]);
      }
      for(i=0;i<Npsi;i++)
      {
        v[vip][i]=v[3][i];
      }
      // remove small values from v[vip]
      ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi, epsilon);
    }
    // remove vim from memory and set vim to vi and vi to vip
    if(k==1)
    {
      // nothing to free
    }
    else if(k==2)
    {
      for(i=0;i<Npsi;i++)
      {
        ComplexWaveFunctionFreeValues(&(v[vim][i]));
      }
    }
    else
    {
      BlockWaveFunctionFree(v[vim], Npsivim);
    }
    vim=vi;
    vi=vip;
    vip=(vip+1)%3;
    Npsivim=Npsivi;
    Npsivi=Npsi;
  }
  if(M->N<=1)
  {
    for(i=0;i<Npsivi;i++)
    {
      ComplexWaveFunctionFreeValues(&(v[1][i]));
    }
  }
  else if(M->N==2)
  {
    for(i=0;i<Npsivim;i++)
    {
      ComplexWaveFunctionFreeValues(&(v[vim][i]));
    }
    BlockWaveFunctionFree(v[vi],  Npsivi);
  }
  else
  {
    BlockWaveFunctionFree(v[vim], Npsivim);
    BlockWaveFunctionFree(v[vi],  Npsivi);
  }
  for(i=0;i<Sqr(M->BlockSize);i++)
  {
    M->a[i].Re=0;
    M->a[i].Im=0;
  }
  free(v[1]);
  free(v[2]);
  free(v[3]);
  return 0;
}
uint32_t LanczosBlockTriDiagonalizeRestrictedMC(const OperatorType H, WaveFunctionType * const psi, DeterminantRestrictionType Restrictions, double epsilon, BlockTriDiagonalMatrixType * const M, const double Gamma)
{
  uint32_t i,j,k,vi,vip,vim;
  CompactMatrixType Mab;
  WaveFunctionType *v[4];
  uint32_t NSingular,Npsi=M->BlockSize,Npsivi=M->BlockSize,Npsivim=M->BlockSize;
  VerbosityPrintf(0x0001,"Start of LanczosBlockTriDiagonalizeRestrictedMC\n");
  fflush(stdout);
  for(i=0;i<M->BlockSize;i++)
  {
    M->b[i+i*M->BlockSize].Re=0;
  }
  v[0] = psi;
  v[1] = (WaveFunctionType*)malloc(Npsi*sizeof(WaveFunctionType));
  v[2] = (WaveFunctionType*)malloc(Npsi*sizeof(WaveFunctionType));
  v[3] = (WaveFunctionType*)malloc(Npsi*sizeof(WaveFunctionType));
  Mab.Ni=M->BlockSize;
  Mab.Nj=M->BlockSize;
  Mab.QComplex=1;
  ComplexWaveFunctionListSVDOrthogonalizeSameBasis(v[0], v[1], &Npsi, epsilon, true); // v[1] = orthonormal[v0]
  NSingular=M->BlockSize-Npsi;
  if(NSingular!=0)
  {
    VerbosityPrintf(0x0001,"Reduce dimension from %u to %u in block %u\n",M->BlockSize,M->BlockSize-NSingular,0);
  }
  // b[0] = v[0] . v[1]^+ or b[0]^T = v[1]^* . v[0]^T
  for(i=0;i<Npsi+NSingular;i++)
  {
    for(j=0;j<Npsi;j++)
    {
      ComplexWaveFunctionDotProductSameBasis(v[1][j], v[0][i], &(M->b[i * M->BlockSize + j].Re), &(M->b[i * M->BlockSize + j].Im));
    }
  }
  vi=1;
  vim=0;
  vip=2;
  Npsivi=Npsi;
  for(k=1;k<M->N;k++)
  {
    // v[vip] = H v[vi]^T
    while(BlockOperatorPsiRestrictedMC(H, v[vi], v[vip], Restrictions, Npsivi))
    {
      printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRestrictedMC for the function BlockOperatorPsiRestrictedMC\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      if (v[vi][0].NBlok == 1) AllocFailed();
      ComplexBlockWaveFunctionTruncateBasis(v[vi], Npsivi, epsilon);
    }
    // a[k] = v[vi]^* . v[vip]^T  (v[vi]^* . H . v[vi]^T)
    Mab.V = &(M->a[k*Sqr(M->BlockSize)]);
    ComplexBlockWaveFunctionDotProduct(v[vi], Npsi, v[vip], Npsi, &Mab); //Npsivi = Npsi - changes at Lowdin Orthonormalization
    if(k!=1)
    {
      // b[k-1] = v[vim]^* . v[vip]^T  (v[vim]^* . H . v[vi]^T)
      Mab.V = &(M->b[(k-1)*Sqr(M->BlockSize)]);
      ComplexBlockWaveFunctionDotProduct(v[vim], Npsivim, v[vip], Npsi, &Mab);
    }
    NSingular=1;
    while(NSingular)
    {
      // orthogonalize v[vip] to v[vi] and v[vim]
      // we use matrix multiplication between v[vip] and v[vi] which is a numerical unstable form of Gram Schmidt
      // we heal this by repeating the orthogonalization twice
      Mab.V=M->a;
      // w = v[vip]^* . v[vi]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vi], Npsivi, &Mab);
      // v[vip] -= w^* . v[vi]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vi], Npsivi, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRestrictedMC for the function ComplexBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vi],  Npsivi, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,   epsilon);
      }
      // w = v[vip]^* . v[vi]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vi], Npsivi, &Mab);
      // v[vip] -= w^* . v[vi]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vi], Npsivi, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRestrictedMC for the function ComplexBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vi],  Npsivi, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,   epsilon);
      }
      // w = v[vip]^* . v[vim]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vim], Npsivim, &Mab);
      // v[vip] -= w^* . v[vim]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vim], Npsivim, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRestrictedMC for the function ComplexBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vim], Npsivim, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,    epsilon);
      }
      // w = v[vip]^* . v[vim]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vim], Npsivim, &Mab);
      // v[vip] -= w^* . v[vim]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vim], Npsivim, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRestrictedMC for the function ComplexBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vim], Npsivim, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,    epsilon);
      }
      // w = v[vip]^* . v[vi]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vi], Npsivi, &Mab);
      // v[vip] -= w^* . v[vi]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vi], Npsivi, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRestrictedMC for the function ComplexBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vi],  Npsivi, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,   epsilon);
      }
      // w = v[vip]^* . v[vim]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vim], Npsivim, &Mab);
      // v[vip] -= w^* . v[vim]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vim], Npsivim, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRestrictedMC for the function ComplexBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vim], Npsivim, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,    epsilon);
      }
      // Lowdin orthonormalize v[vip] into v[3], free v[vip] and set v[vip] equal to v[3]
      NSingular=Npsi;
      ComplexWaveFunctionListSVDOrthogonalizeSameBasis(v[vip], v[3], &Npsi, epsilon, true);
      NSingular-=Npsi;
      if(NSingular!=0)
      {
        VerbosityPrintf(0x0001, "Reduce dimension from %u to %u in block %u\n",Npsi+NSingular,Npsi,k);
      }
      for(i=0;i<Npsi+NSingular;i++)
      {
        ComplexWaveFunctionFreeValues(&v[vip][i]);
      }
      for(i=0;i<Npsi;i++)
      {
        v[vip][i]=v[3][i];
      }
      // remove small values from v[vip]
      ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi, epsilon);
    }
    // remove vim from memory and set vim to vi and vi to vip
    if(k==1)
    {
      // nothing to free
    }
    else if(k==2)
    {
      for(i=0;i<Npsi;i++)
      {
        ComplexWaveFunctionFreeValues(&(v[vim][i]));
      }
    }
    else
    {
      BlockWaveFunctionFree(v[vim], Npsivim);
    }
    vim=vi;
    vi=vip;
    vip=(vip+1)%3;
    Npsivim=Npsivi;
    Npsivi=Npsi;
  }
  if(M->N<=1)
  {
    for(i=0;i<Npsivi;i++)
    {
      ComplexWaveFunctionFreeValues(&(v[1][i]));
    }
  }
  else if(M->N==2)
  {
    for(i=0;i<Npsivim;i++)
    {
      ComplexWaveFunctionFreeValues(&(v[vim][i]));
    }
    BlockWaveFunctionFree(v[vi],  Npsivi);
  }
  else
  {
    BlockWaveFunctionFree(v[vim], Npsivim);
    BlockWaveFunctionFree(v[vi],  Npsivi);
  }
  for(i=0;i<Sqr(M->BlockSize);i++)
  {
    M->a[i].Re=0;
    M->a[i].Im=0;
  }
  free(v[1]);
  free(v[2]);
  free(v[3]);
  return 0;
}
uint32_t LanczosBlockTriDiagonalizeRestrictedCC(const OperatorType H, WaveFunctionType * const psi, DeterminantRestrictionType Restrictions, double epsilon, BlockTriDiagonalMatrixType * const M, const double Gamma)
{
  uint32_t i,j,k,vi,vip,vim;
  CompactMatrixType Mab;
  WaveFunctionType *v[4];
  uint32_t NSingular,Npsi=M->BlockSize,Npsivi=M->BlockSize,Npsivim=M->BlockSize;
  VerbosityPrintf(0x0001,"Start of LanczosBlockTriDiagonalizeRestrictedCC\n");
  fflush(stdout);
  for(i=0;i<M->BlockSize;i++)
  {
    M->b[i+i*M->BlockSize].Re=0;
  }
  v[0] = psi;
  v[1] = (WaveFunctionType*)malloc(Npsi*sizeof(WaveFunctionType));
  v[2] = (WaveFunctionType*)malloc(Npsi*sizeof(WaveFunctionType));
  v[3] = (WaveFunctionType*)malloc(Npsi*sizeof(WaveFunctionType));
  Mab.Ni=M->BlockSize;
  Mab.Nj=M->BlockSize;
  Mab.QComplex=1;
  ComplexWaveFunctionListSVDOrthogonalizeSameBasis(v[0], v[1], &Npsi, epsilon, true); // v[1] = orthonormal[v0]
  NSingular=M->BlockSize-Npsi;
  if(NSingular!=0)
  {
    VerbosityPrintf(0x0001,"Reduce dimension from %u to %u in block %u\n",M->BlockSize,M->BlockSize-NSingular,0);
  }
  // b[0] = v[0] . v[1]^+ or b[0]^T = v[1]^* . v[0]^T
  for(i=0;i<Npsi+NSingular;i++)
  {
    for(j=0;j<Npsi;j++)
    {
      ComplexWaveFunctionDotProductSameBasis(v[1][j], v[0][i], &(M->b[i * M->BlockSize + j].Re), &(M->b[i * M->BlockSize + j].Im));
    }
  }
  vi=1;
  vim=0;
  vip=2;
  Npsivi=Npsi;
  for(k=1;k<M->N;k++)
  {
    // v[vip] = H v[vi]^T
    while(BlockOperatorPsiRestrictedCC(H, v[vi], v[vip], Restrictions, Npsivi))
    {
      printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRestrictedCC for the function BlockOperatorPsiRestrictedCC\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      if (v[vi][0].NBlok == 1) AllocFailed();
      ComplexBlockWaveFunctionTruncateBasis(v[vi], Npsivi, epsilon);
    }
    // a[k] = v[vi]^* . v[vip]^T  (v[vi]^* . H . v[vi]^T)
    Mab.V = &(M->a[k*Sqr(M->BlockSize)]);
    ComplexBlockWaveFunctionDotProduct(v[vi], Npsi, v[vip], Npsi, &Mab); //Npsivi = Npsi - changes at Lowdin Orthonormalization
    if(k!=1)
    {
      // b[k-1] = v[vim]^* . v[vip]^T  (v[vim]^* . H . v[vi]^T)
      Mab.V = &(M->b[(k-1)*Sqr(M->BlockSize)]);
      ComplexBlockWaveFunctionDotProduct(v[vim], Npsivim, v[vip], Npsi, &Mab);
    }
    NSingular=1;
    while(NSingular)
    {
      // orthogonalize v[vip] to v[vi] and v[vim]
      // we use matrix multiplication between v[vip] and v[vi] which is a numerical unstable form of Gram Schmidt
      // we heal this by repeating the orthogonalization twice and once for vi and vim (factor of 3 reduction in speed but needed in tests)
      Mab.V=M->a;
      // w = v[vip]^* . v[vi]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vi], Npsivi, &Mab);
      // v[vip] -= w^* . v[vi]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vi], Npsivi, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRestrictedCC for the function ComplexBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vi],  Npsivi, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,   epsilon);
      }
      // w = v[vip]^* . v[vi]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vi], Npsivi, &Mab);
      // v[vip] -= w^* . v[vi]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vi], Npsivi, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRestrictedCC for the function ComplexBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vi],  Npsivi, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,   epsilon);
      }
      // w = v[vip]^* . v[vim]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vim], Npsivim, &Mab);
      // v[vip] -= w^* . v[vim]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vim], Npsivim, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRestrictedCC for the function ComplexBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vim][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vim],  Npsivim, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip],  Npsi,   epsilon);
      }
      // w = v[vip]^* . v[vim]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vim], Npsivim, &Mab);
      // v[vip] -= w^* . v[vim]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vim], Npsivim, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRestrictedCC for the function ComplexBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vim][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vim],  Npsivim, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip],  Npsi,   epsilon);
      }
      // w = v[vip]^* . v[vi]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vi], Npsivi, &Mab);
      // v[vip] -= w^* . v[vi]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vi], Npsivi, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRestrictedCC for the function ComplexBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vi],  Npsivi, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,   epsilon);
      }
      // w = v[vip]^* . v[vim]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vim], Npsivim, &Mab);
      // v[vip] -= w^* . v[vim]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vim], Npsivim, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRestrictedCC for the function ComplexBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vim][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vim],  Npsivim, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip],  Npsi,   epsilon);
      }
      // Lowdin orthonormalize v[vip] into v[3], free v[vip] and set v[vip] equal to v[3]
      NSingular=Npsi;
      ComplexWaveFunctionListSVDOrthogonalizeSameBasis(v[vip], v[3], &Npsi, epsilon, true);
      NSingular-=Npsi;
      if(NSingular!=0)
      {
        VerbosityPrintf(0x0001, "Reduce dimension from %u to %u in block %u\n",Npsi+NSingular,Npsi,k);
      }
      for(i=0;i<Npsi+NSingular;i++)
      {
        ComplexWaveFunctionFreeValues(&v[vip][i]);
      }
      for(i=0;i<Npsi;i++)
      {
        v[vip][i]=v[3][i];
      }
      // remove small values from v[vip]
      ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi, epsilon);
    }
    // remove vim from memory and set vim to vi and vi to vip
    if(k==1)
    {
      // nothing to free
    }
    else if(k==2)
    {
      for(i=0;i<Npsi;i++)
      {
        ComplexWaveFunctionFreeValues(&(v[vim][i]));
      }
    }
    else
    {
      BlockWaveFunctionFree(v[vim], Npsivim);
    }
    vim=vi;
    vi=vip;
    vip=(vip+1)%3;
    Npsivim=Npsivi;
    Npsivi=Npsi;
  }
  if(M->N<=1)
  {
    for(i=0;i<Npsivi;i++)
    {
      ComplexWaveFunctionFreeValues(&(v[1][i]));
    }
  }
  else if(M->N==2)
  {
    for(i=0;i<Npsivim;i++)
    {
      ComplexWaveFunctionFreeValues(&(v[vim][i]));
    }
    BlockWaveFunctionFree(v[vi],  Npsivi);
  }
  else
  {
    BlockWaveFunctionFree(v[vim], Npsivim);
    BlockWaveFunctionFree(v[vi],  Npsivi);
  }
  for(i=0;i<Sqr(M->BlockSize);i++)
  {
    M->a[i].Re=0;
    M->a[i].Im=0;
  }
  free(v[1]);
  free(v[2]);
  free(v[3]);
  return 0;
}
uint32_t LanczosBlockTriDiagonalizeRestrictedRR(const OperatorType H, WaveFunctionType * const psi, DeterminantRestrictionType Restrictions, double epsilon, BlockTriDiagonalMatrixType * const M, const double Gamma)
{
  uint32_t i,j,k,vi,vip,vim;
  CompactMatrixType Mab;
  WaveFunctionType *v[4];
  uint32_t NSingular,Npsi=M->BlockSize,Npsivi=M->BlockSize,Npsivim=M->BlockSize;
  VerbosityPrintf(0x0001,"Start of LanczosBlockTriDiagonalizeRestrictedRR\n");
  fflush(stdout);
  for(i=0;i<M->BlockSize;i++)
  {
    M->bRe[i+i*M->BlockSize]=0;
  }
  v[0] = psi;
  v[1] = (WaveFunctionType*)malloc(Npsi*sizeof(WaveFunctionType));
  v[2] = (WaveFunctionType*)malloc(Npsi*sizeof(WaveFunctionType));
  v[3] = (WaveFunctionType*)malloc(Npsi*sizeof(WaveFunctionType));
  Mab.Ni=M->BlockSize;
  Mab.Nj=M->BlockSize;
  Mab.QComplex=0;
  RealWaveFunctionListSVDOrthogonalizeSameBasis(v[0], v[1], &Npsi, epsilon, true); // v[1] = orthonormal[v0]
  NSingular=M->BlockSize-Npsi;
  if(NSingular!=0)
  {
    VerbosityPrintf(0x0001,"Reduce dimension from %u to %u in block %u\n",M->BlockSize,M->BlockSize-NSingular,0);
  }
  // b[0] = v[0] . v[1]^+ or b[0]^T = v[1]^* . v[0]^T
  for(i=0;i<Npsi+NSingular;i++)
  {
    for(j=0;j<Npsi;j++)
    {
      M->bRe[i * M->BlockSize + j] = RealWaveFunctionDotProductSameBasis(v[1][j], v[0][i]);
    }
  }
  vi=1;
  vim=0;
  vip=2;
  Npsivi=Npsi;
  for(k=1;k<M->N;k++)
  {
    // v[vip] = H v[vi]^T
    while(BlockOperatorPsiRestrictedRR(H, v[vi], v[vip], Restrictions, Npsivi))
    {
      printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRestrictedRR for the function BlockOperatorPsiRestrictedRR\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      if (v[vi][0].NBlok == 1) AllocFailed();
      RealBlockWaveFunctionTruncateBasis(v[vi], Npsivi, epsilon);
    }
    // a[k] = v[vi]^* . v[vip]^T  (v[vi]^* . H . v[vi]^T)
    Mab.VRe = &(M->aRe[k*Sqr(M->BlockSize)]);
    RealBlockWaveFunctionDotProduct(v[vi], Npsi, v[vip], Npsi, &Mab); //Npsivi = Npsi - changes at Lowdin Orthonormalization
    if(k!=1)
    {
      // b[k-1] = v[vim]^* . v[vip]^T  (v[vim]^* . H . v[vi]^T)
      Mab.VRe = &(M->bRe[(k-1)*Sqr(M->BlockSize)]);
      RealBlockWaveFunctionDotProduct(v[vim], Npsivim, v[vip], Npsi, &Mab);
    }
    NSingular=1;
    while(NSingular)
    {
      // orthogonalize v[vip] to v[vi] and v[vim]
      // we use matrix multiplication between v[vip] and v[vi] which is a numerical unstable form of Gram Schmidt
      // we heal this by repeating the orthogonalization twice for each matrix vector and once for vi and vim (factor of 3 reduction in speed but needed in tests)
      Mab.VRe=M->aRe;
      // w = v[vip]^* . v[vi]^T
      CompactMatrixScale(&Mab, 0);
      RealBlockWaveFunctionDotProduct(v[vip], Npsi, v[vi], Npsivi, &Mab);
      // v[vip] -= w^* . v[vi]
      while(RealBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vi], Npsivi, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRestrictedRR for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        RealBlockWaveFunctionTruncateBasis(v[vi],  Npsivi, epsilon);
        RealBlockWaveFunctionTruncateBasis(v[vip], Npsi,   epsilon);
      }
      // w = v[vip]^* . v[vi]^T
      CompactMatrixScale(&Mab, 0);
      RealBlockWaveFunctionDotProduct(v[vip], Npsi, v[vi], Npsivi, &Mab);
      // v[vip] -= w^* . v[vi]
      while(RealBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vi], Npsivi, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRestrictedRR for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        RealBlockWaveFunctionTruncateBasis(v[vi],  Npsivi, epsilon);
        RealBlockWaveFunctionTruncateBasis(v[vip], Npsi,   epsilon);
      }
      // w = v[vip]^* . v[vim]^T
      CompactMatrixScale(&Mab, 0);
      RealBlockWaveFunctionDotProduct(v[vip], Npsi, v[vim], Npsivim, &Mab);
      // v[vip] -= w^* . v[vim]
      while(RealBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vim], Npsivim, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRestrictedRR for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vim][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        RealBlockWaveFunctionTruncateBasis(v[vim], Npsivim, epsilon);
        RealBlockWaveFunctionTruncateBasis(v[vip], Npsi,    epsilon);
      }
      // w = v[vip]^* . v[vim]^T
      CompactMatrixScale(&Mab, 0);
      RealBlockWaveFunctionDotProduct(v[vip], Npsi, v[vim], Npsivim, &Mab);
      // v[vip] -= w^* . v[vim]
      while(RealBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vim], Npsivim, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRestrictedRR for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vim][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        RealBlockWaveFunctionTruncateBasis(v[vim], Npsivim, epsilon);
        RealBlockWaveFunctionTruncateBasis(v[vip], Npsi,    epsilon);
      }
      // w = v[vip]^* . v[vi]^T
      CompactMatrixScale(&Mab, 0);
      RealBlockWaveFunctionDotProduct(v[vip], Npsi, v[vi], Npsivi, &Mab);
      // v[vip] -= w^* . v[vi]
      while(RealBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vi], Npsivi, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRestrictedRR for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        RealBlockWaveFunctionTruncateBasis(v[vi],  Npsivi, epsilon);
        RealBlockWaveFunctionTruncateBasis(v[vip], Npsi,   epsilon);
      }
      // w = v[vip]^* . v[vim]^T
      CompactMatrixScale(&Mab, 0);
      RealBlockWaveFunctionDotProduct(v[vip], Npsi, v[vim], Npsivim, &Mab);
      // v[vip] -= w^* . v[vim]
      while(RealBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vim], Npsivim, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRestrictedRR for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vim][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        RealBlockWaveFunctionTruncateBasis(v[vim], Npsivim, epsilon);
        RealBlockWaveFunctionTruncateBasis(v[vip], Npsi,    epsilon);
      }
      // Lowdin orthonormalize v[vip] into v[3], free v[vip] and set v[vip] equal to v[3]
      NSingular=Npsi;
      RealWaveFunctionListSVDOrthogonalizeSameBasis(v[vip], v[3], &Npsi, epsilon, true);
      NSingular-=Npsi;
      if(NSingular!=0)
      {
        VerbosityPrintf(0x0001, "Reduce dimension from %u to %u in block %u\n",Npsi+NSingular,Npsi,k);
      }
      for(i=0;i<Npsi+NSingular;i++)
      {
        WaveFunctionFreeValues(&v[vip][i]);
      }
      for(i=0;i<Npsi;i++)
      {
        v[vip][i]=v[3][i];
      }
      // remove small values from v[vip]
      BlockWaveFunctionTruncateBasis(v[vip], Npsi, epsilon);
    }
    // remove vim from memory and set vim to vi and vi to vip
    if(k==1)
    {
      // nothing to free
    }
    else if(k==2)
    {
      for(i=0;i<Npsi;i++)
      {
        RealWaveFunctionFreeValues(&(v[vim][i]));
      }
    }
    else
    {
      BlockWaveFunctionFree(v[vim], Npsivim);
    }
    vim=vi;
    vi=vip;
    vip=(vip+1)%3;
    Npsivim=Npsivi;
    Npsivi=Npsi;
  }
  if(M->N<=1)
  {
    for(i=0;i<Npsivi;i++)
    {
      RealWaveFunctionFreeValues(&(v[1][i]));
    }
  }
  else if(M->N==2)
  {
    for(i=0;i<Npsivim;i++)
    {
      RealWaveFunctionFreeValues(&(v[vim][i]));
    }
    BlockWaveFunctionFree(v[vi],  Npsivi);
  }
  else
  {
    BlockWaveFunctionFree(v[vim], Npsivim);
    BlockWaveFunctionFree(v[vi],  Npsivi);
  }
  for(i=0;i<Sqr(M->BlockSize);i++)
  {
    M->aRe[i]=0;
  }
  free(v[1]);
  free(v[2]);
  free(v[3]);
  return 0;
}
uint32_t LanczosBlockTriDiagonalizeRestrictedMR(const OperatorType H, WaveFunctionType * const psi, DeterminantRestrictionType Restrictions, double epsilon, BlockTriDiagonalMatrixType * const M, const double Gamma)
{
  // make psi complex, call MC
  uint32_t i;
  VerbosityPrintf(0x0001,"Start of LanczosBlockTriDiagonalizeRestrictedMR\n");
  fflush(stdout);
  for(i=0;i<M->BlockSize;i++)
  {
    if (MakeComplexWaveFunctionFromReal(&psi[i], psi[i]))
    {
      printf("LanczosBlockTriDiagonalizeRestrictedMR called with a real wave-function and a complex Hamiltonian.\nTurned the wave-function complex in step 1.\nBut failed at function %i.\n",i);
      fflush(stdout);
      for(;i>0;i--)
      {
        MakeRealWaveFunctionFromComplex(psi[i-1], &psi[i-1]);
      }
      return 1;
    }
  }
  return LanczosBlockTriDiagonalizeRestrictedMC(H,psi,Restrictions,epsilon,M,Gamma);
}
uint32_t LanczosBlockTriDiagonalizeRestrictedCR(const OperatorType H, WaveFunctionType * const psi, DeterminantRestrictionType Restrictions, double epsilon, BlockTriDiagonalMatrixType * const M, const double Gamma)
{
  // make psi complex, and call CC
  uint32_t i;
  VerbosityPrintf(0x0001,"Start of LanczosBlockTriDiagonalizeRestrictedCR\n");
  fflush(stdout);
  for(i=0;i<M->BlockSize;i++)
  {
    if (MakeComplexWaveFunctionFromReal(&psi[i], psi[i]))
    {
      printf("LanczosBlockTriDiagonalizeRestrictedCR called with a real wave-function and a complex Hamiltonian.\nTurned the wave-function complex in step 1.\nBut failed at function %i.\n",i);
      fflush(stdout);
      for(;i>0;i--)
      {
        MakeRealWaveFunctionFromComplex(psi[i-1], &psi[i-1]);
      }
      return 1;
    }
  }
  return LanczosBlockTriDiagonalizeRestrictedCC(H,psi,Restrictions,epsilon,M,Gamma);
}



uint32_t LanczosBlockTriDiagonalizeRC(const OperatorType H, WaveFunctionType * const psi, double epsilon, BlockTriDiagonalMatrixType * const M, const double Gamma)
{
  uint32_t i,j,k,vi,vip,vim;
  CompactMatrixType Mab;
  WaveFunctionType *v[4];
  uint32_t NSingular,Npsi=M->BlockSize,Npsivi=M->BlockSize,Npsivim=M->BlockSize;
  VerbosityPrintf(0x0001,"Start of LanczosBlockTriDiagonalizeRC\n");
  fflush(stdout);
  for(i=0;i<M->BlockSize;i++)
  {
    M->b[i+i*M->BlockSize].Re=0;
  }
  v[0] = psi;
  v[1] = (WaveFunctionType*)malloc(Npsi*sizeof(WaveFunctionType));
  v[2] = (WaveFunctionType*)malloc(Npsi*sizeof(WaveFunctionType));
  v[3] = (WaveFunctionType*)malloc(Npsi*sizeof(WaveFunctionType));
  Mab.Ni=M->BlockSize;
  Mab.Nj=M->BlockSize;
  Mab.QComplex=1;
  ComplexWaveFunctionListSVDOrthogonalizeSameBasis(v[0], v[1], &Npsi, epsilon, true); // v[1] = orthonormal[v0]
  NSingular=M->BlockSize-Npsi;
  if(NSingular!=0)
  {
    VerbosityPrintf(0x0001,"Reduce dimension from %u to %u in block %u\n",M->BlockSize,M->BlockSize-NSingular,0);
  }
  // b[0] = v[0] . v[1]^+ or b[0]^T = v[1]^* . v[0]^T
  for(i=0;i<Npsi+NSingular;i++)
  {
    for(j=0;j<Npsi;j++)
    {
      ComplexWaveFunctionDotProductSameBasis(v[1][j], v[0][i], &(M->b[i * M->BlockSize + j].Re), &(M->b[i * M->BlockSize + j].Im));
    }
  }
  vi=1;
  vim=0;
  vip=2;
  Npsivi=Npsi;
  for(k=1;k<M->N;k++)
  {
    // v[vip] = H v[vi]^T
    while(BlockOperatorPsiRC(H, v[vi], v[vip], Npsivi))
    {
      printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRR for the function BlockOperatorPsiRR\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      if (v[vi][0].NBlok == 1) AllocFailed();
      ComplexBlockWaveFunctionTruncateBasis(v[vi], Npsivi, epsilon);
    }
    // a[k] = v[vi]^* . v[vip]^T  (v[vi]^* . H . v[vi]^T)
    Mab.V = &(M->a[k*Sqr(M->BlockSize)]);
    ComplexBlockWaveFunctionDotProduct(v[vi], Npsi, v[vip], Npsi, &Mab); //Npsivi = Npsi - changes at Lowdin Orthonormalization
    if(k!=1)
    {
      // b[k-1] = v[vim]^* . v[vip]^T  (v[vim]^* . H . v[vi]^T)
      Mab.V = &(M->b[(k-1)*Sqr(M->BlockSize)]);
      ComplexBlockWaveFunctionDotProduct(v[vim], Npsivim, v[vip], Npsi, &Mab);
    }
    NSingular=1;
    while(NSingular)
    {
      // orthogonalize v[vip] to v[vi] and v[vim]
      // we use matrix multiplication between v[vip] and v[vi] which is a numerical unstable form of Gram Schmidt
      // we heal this by repeating the orthogonalization twice
      Mab.V=M->a;
      // w = v[vip]^* . v[vi]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vi], Npsivi, &Mab);
      // v[vip] -= w^* . v[vi]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vi], Npsivi, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRC for the function ComplexBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vi],  Npsivi, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,   epsilon);
      }
      // w = v[vip]^* . v[vi]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vi], Npsivi, &Mab);
      // v[vip] -= w^* . v[vi]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vi], Npsivi, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRC for the function ComplexBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vi],  Npsivi, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,   epsilon);
      }
      // w = v[vip]^* . v[vim]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vim], Npsivim, &Mab);
      // v[vip] -= w^* . v[vim]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vim], Npsivim, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRC for the function ComplexBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vim], Npsivim, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,    epsilon);
      }
      // w = v[vip]^* . v[vim]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vim], Npsivim, &Mab);
      // v[vip] -= w^* . v[vim]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vim], Npsivim, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRC for the function ComplexBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vim], Npsivim, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,    epsilon);
      }
      // w = v[vip]^* . v[vi]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vi], Npsivi, &Mab);
      // v[vip] -= w^* . v[vi]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vi], Npsivi, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRC for the function ComplexBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vi],  Npsivi, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,   epsilon);
      }
      // w = v[vip]^* . v[vim]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vim], Npsivim, &Mab);
      // v[vip] -= w^* . v[vim]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vim], Npsivim, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRC for the function ComplexBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vim], Npsivim, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,    epsilon);
      }
      // Lowdin orthonormalize v[vip] into v[3], free v[vip] and set v[vip] equal to v[3]
      NSingular=Npsi;
      ComplexWaveFunctionListSVDOrthogonalizeSameBasis(v[vip], v[3], &Npsi, epsilon, true);
      NSingular-=Npsi;
      if(NSingular!=0)
      {
        VerbosityPrintf(0x0001, "Reduce dimension from %u to %u in block %u\n",Npsi+NSingular,Npsi,k);
      }
      for(i=0;i<Npsi+NSingular;i++)
      {
        ComplexWaveFunctionFreeValues(&v[vip][i]);
      }
      for(i=0;i<Npsi;i++)
      {
        v[vip][i]=v[3][i];
      }
      // remove small values from v[vip]
      ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi, epsilon);
    }
    // remove vim from memory and set vim to vi and vi to vip
    if(k==1)
    {
      // nothing to free
    }
    else if(k==2)
    {
      for(i=0;i<Npsi;i++)
      {
        ComplexWaveFunctionFreeValues(&(v[vim][i]));
      }
    }
    else
    {
      BlockWaveFunctionFree(v[vim], Npsivim);
    }
    vim=vi;
    vi=vip;
    vip=(vip+1)%3;
    Npsivim=Npsivi;
    Npsivi=Npsi;
  }
  if(M->N<=1)
  {
    for(i=0;i<Npsivi;i++)
    {
      ComplexWaveFunctionFreeValues(&(v[1][i]));
    }
  }
  else if(M->N==2)
  {
    for(i=0;i<Npsivim;i++)
    {
      ComplexWaveFunctionFreeValues(&(v[vim][i]));
    }
    BlockWaveFunctionFree(v[vi],  Npsivi);
  }
  else
  {
    BlockWaveFunctionFree(v[vim], Npsivim);
    BlockWaveFunctionFree(v[vi],  Npsivi);
  }
  for(i=0;i<Sqr(M->BlockSize);i++)
  {
    M->a[i].Re=0;
    M->a[i].Im=0;
  }
  free(v[1]);
  free(v[2]);
  free(v[3]);
  return 0;
}
uint32_t LanczosBlockTriDiagonalizeMC(const OperatorType H, WaveFunctionType * const psi, double epsilon, BlockTriDiagonalMatrixType * const M, const double Gamma)
{
  uint32_t i,j,k,vi,vip,vim;
  CompactMatrixType Mab;
  WaveFunctionType *v[4];
  uint32_t NSingular,Npsi=M->BlockSize,Npsivi=M->BlockSize,Npsivim=M->BlockSize;
  VerbosityPrintf(0x0001,"Start of LanczosBlockTriDiagonalizeMC\n");
  fflush(stdout);
  for(i=0;i<M->BlockSize;i++)
  {
    M->b[i+i*M->BlockSize].Re=0;
  }
  v[0] = psi;
  v[1] = (WaveFunctionType*)malloc(Npsi*sizeof(WaveFunctionType));
  v[2] = (WaveFunctionType*)malloc(Npsi*sizeof(WaveFunctionType));
  v[3] = (WaveFunctionType*)malloc(Npsi*sizeof(WaveFunctionType));
  Mab.Ni=M->BlockSize;
  Mab.Nj=M->BlockSize;
  Mab.QComplex=1;
  ComplexWaveFunctionListSVDOrthogonalizeSameBasis(v[0], v[1], &Npsi, epsilon, true); // v[1] = orthonormal[v0]
  NSingular=M->BlockSize-Npsi;
  if(NSingular!=0)
  {
    VerbosityPrintf(0x0001,"Reduce dimension from %u to %u in block %u\n",M->BlockSize,M->BlockSize-NSingular,0);
  }
  // b[0] = v[0] . v[1]^+ or b[0]^T = v[1]^* . v[0]^T
  for(i=0;i<Npsi+NSingular;i++)
  {
    for(j=0;j<Npsi;j++)
    {
      ComplexWaveFunctionDotProductSameBasis(v[1][j], v[0][i], &(M->b[i * M->BlockSize + j].Re), &(M->b[i * M->BlockSize + j].Im));
    }
  }
  vi=1;
  vim=0;
  vip=2;
  Npsivi=Npsi;
  for(k=1;k<M->N;k++)
  {
    // v[vip] = H v[vi]^T
    while(BlockOperatorPsiMC(H, v[vi], v[vip], Npsivi))
    {
      printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRR for the function BlockOperatorPsiRR\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      if (v[vi][0].NBlok == 1) AllocFailed();
      ComplexBlockWaveFunctionTruncateBasis(v[vi], Npsivi, epsilon);
    }
    // a[k] = v[vi]^* . v[vip]^T  (v[vi]^* . H . v[vi]^T)
    Mab.V = &(M->a[k*Sqr(M->BlockSize)]);
    ComplexBlockWaveFunctionDotProduct(v[vi], Npsi, v[vip], Npsi, &Mab); //Npsivi = Npsi - changes at Lowdin Orthonormalization
    if(k!=1)
    {
      // b[k-1] = v[vim]^* . v[vip]^T  (v[vim]^* . H . v[vi]^T)
      Mab.V = &(M->b[(k-1)*Sqr(M->BlockSize)]);
      ComplexBlockWaveFunctionDotProduct(v[vim], Npsivim, v[vip], Npsi, &Mab);
    }
    NSingular=1;
    while(NSingular)
    {
      // orthogonalize v[vip] to v[vi] and v[vim]
      // we use matrix multiplication between v[vip] and v[vi] which is a numerical unstable form of Gram Schmidt
      // we heal this by repeating the orthogonalization twice
      Mab.V=M->a;
      // w = v[vip]^* . v[vi]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vi], Npsivi, &Mab);
      // v[vip] -= w^* . v[vi]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vi], Npsivi, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeMC for the function ComplexBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vi],  Npsivi, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,   epsilon);
      }
      // w = v[vip]^* . v[vi]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vi], Npsivi, &Mab);
      // v[vip] -= w^* . v[vi]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vi], Npsivi, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeMC for the function ComplexBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vi],  Npsivi, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,   epsilon);
      }
      // w = v[vip]^* . v[vim]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vim], Npsivim, &Mab);
      // v[vip] -= w^* . v[vim]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vim], Npsivim, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeMC for the function ComplexBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vim], Npsivim, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,    epsilon);
      }
      // w = v[vip]^* . v[vim]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vim], Npsivim, &Mab);
      // v[vip] -= w^* . v[vim]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vim], Npsivim, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeMC for the function ComplexBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vim], Npsivim, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,    epsilon);
      }
      // w = v[vip]^* . v[vi]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vi], Npsivi, &Mab);
      // v[vip] -= w^* . v[vi]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vi], Npsivi, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeMC for the function ComplexBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vi],  Npsivi, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,   epsilon);
      }
      // w = v[vip]^* . v[vim]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vim], Npsivim, &Mab);
      // v[vip] -= w^* . v[vim]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vim], Npsivim, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeMC for the function ComplexBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vim], Npsivim, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,    epsilon);
      }
      // Lowdin orthonormalize v[vip] into v[3], free v[vip] and set v[vip] equal to v[3]
      NSingular=Npsi;
      ComplexWaveFunctionListSVDOrthogonalizeSameBasis(v[vip], v[3], &Npsi, epsilon, true);
      NSingular-=Npsi;
      if(NSingular!=0)
      {
        VerbosityPrintf(0x0001, "Reduce dimension from %u to %u in block %u\n",Npsi+NSingular,Npsi,k);
      }
      for(i=0;i<Npsi+NSingular;i++)
      {
        ComplexWaveFunctionFreeValues(&v[vip][i]);
      }
      for(i=0;i<Npsi;i++)
      {
        v[vip][i]=v[3][i];
      }
      // remove small values from v[vip]
      ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi, epsilon);
    }
    // remove vim from memory and set vim to vi and vi to vip
    if(k==1)
    {
      // nothing to free
    }
    else if(k==2)
    {
      for(i=0;i<Npsi;i++)
      {
        ComplexWaveFunctionFreeValues(&(v[vim][i]));
      }
    }
    else
    {
      BlockWaveFunctionFree(v[vim], Npsivim);
    }
    vim=vi;
    vi=vip;
    vip=(vip+1)%3;
    Npsivim=Npsivi;
    Npsivi=Npsi;
  }
  if(M->N<=1)
  {
    for(i=0;i<Npsivi;i++)
    {
      ComplexWaveFunctionFreeValues(&(v[1][i]));
    }
  }
  else if(M->N==2)
  {
    for(i=0;i<Npsivim;i++)
    {
      ComplexWaveFunctionFreeValues(&(v[vim][i]));
    }
    BlockWaveFunctionFree(v[vi],  Npsivi);
  }
  else
  {
    BlockWaveFunctionFree(v[vim], Npsivim);
    BlockWaveFunctionFree(v[vi],  Npsivi);
  }
  for(i=0;i<Sqr(M->BlockSize);i++)
  {
    M->a[i].Re=0;
    M->a[i].Im=0;
  }
  free(v[1]);
  free(v[2]);
  free(v[3]);
  return 0;
}
uint32_t LanczosBlockTriDiagonalizeCC(const OperatorType H, WaveFunctionType * const psi, double epsilon, BlockTriDiagonalMatrixType * const M, const double Gamma)
{
  uint32_t i,j,k,vi,vip,vim;
  CompactMatrixType Mab;
  WaveFunctionType *v[4];
  uint32_t NSingular,Npsi=M->BlockSize,Npsivi=M->BlockSize,Npsivim=M->BlockSize;
  VerbosityPrintf(0x0001,"Start of LanczosBlockTriDiagonalizeCC\n");
  fflush(stdout);
  for(i=0;i<M->BlockSize;i++)
  {
    M->b[i+i*M->BlockSize].Re=0;
  }
  v[0] = psi;
  v[1] = (WaveFunctionType*)malloc(Npsi*sizeof(WaveFunctionType));
  v[2] = (WaveFunctionType*)malloc(Npsi*sizeof(WaveFunctionType));
  v[3] = (WaveFunctionType*)malloc(Npsi*sizeof(WaveFunctionType));
  Mab.Ni=M->BlockSize;
  Mab.Nj=M->BlockSize;
  Mab.QComplex=1;
  ComplexWaveFunctionListSVDOrthogonalizeSameBasis(v[0], v[1], &Npsi, epsilon, true); // v[1] = orthonormal[v0]
  NSingular=M->BlockSize-Npsi;
  if(NSingular!=0)
  {
    VerbosityPrintf(0x0001,"Reduce dimension from %u to %u in block %u\n",M->BlockSize,M->BlockSize-NSingular,0);
  }
  // b[0] = v[0] . v[1]^+ or b[0]^T = v[1]^* . v[0]^T
  for(i=0;i<Npsi+NSingular;i++)
  {
    for(j=0;j<Npsi;j++)
    {
      ComplexWaveFunctionDotProductSameBasis(v[1][j], v[0][i], &(M->b[i * M->BlockSize + j].Re), &(M->b[i * M->BlockSize + j].Im));
    }
  }
  vi=1;
  vim=0;
  vip=2;
  Npsivi=Npsi;
  for(k=1;k<M->N;k++)
  {
    // v[vip] = H v[vi]^T
    while(BlockOperatorPsiCC(H, v[vi], v[vip], Npsivi))
    {
      printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeCC for the function BlockOperatorPsiCC\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      if (v[vi][0].NBlok == 1) AllocFailed();
      ComplexBlockWaveFunctionTruncateBasis(v[vi], Npsivi, epsilon);
    }
    // a[k] = v[vi]^* . v[vip]^T  (v[vi]^* . H . v[vi]^T)
    Mab.V = &(M->a[k*Sqr(M->BlockSize)]);
    ComplexBlockWaveFunctionDotProduct(v[vi], Npsi, v[vip], Npsi, &Mab); //Npsivi = Npsi - changes at Lowdin Orthonormalization
    if(k!=1)
    {
      // b[k-1] = v[vim]^* . v[vip]^T  (v[vim]^* . H . v[vi]^T)
      Mab.V = &(M->b[(k-1)*Sqr(M->BlockSize)]);
      ComplexBlockWaveFunctionDotProduct(v[vim], Npsivim, v[vip], Npsi, &Mab);
    }
    NSingular=1;
    while(NSingular)
    {
      // orthogonalize v[vip] to v[vi] and v[vim]
      // we use matrix multiplication between v[vip] and v[vi] which is a numerical unstable form of Gram Schmidt
      // we heal this by repeating the orthogonalization twice and once for vi and vim (factor of 3 reduction in speed but needed in tests)
      Mab.V=M->a;
      // w = v[vip]^* . v[vi]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vi], Npsivi, &Mab);
      // v[vip] -= w^* . v[vi]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vi], Npsivi, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeCC for the function ComplexBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vi],  Npsivi, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,   epsilon);
      }
      // w = v[vip]^* . v[vi]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vi], Npsivi, &Mab);
      // v[vip] -= w^* . v[vi]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vi], Npsivi, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeCC for the function ComplexBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vi],  Npsivi, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,   epsilon);
      }
      // w = v[vip]^* . v[vim]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vim], Npsivim, &Mab);
      // v[vip] -= w^* . v[vim]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vim], Npsivim, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeCC for the function ComplexBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vim][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vim],  Npsivim, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip],  Npsi,   epsilon);
      }
      // w = v[vip]^* . v[vim]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vim], Npsivim, &Mab);
      // v[vip] -= w^* . v[vim]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vim], Npsivim, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeCC for the function ComplexBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vim][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vim],  Npsivim, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip],  Npsi,   epsilon);
      }
      // w = v[vip]^* . v[vi]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vi], Npsivi, &Mab);
      // v[vip] -= w^* . v[vi]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vi], Npsivi, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeCC for the function ComplexBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vi],  Npsivi, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,   epsilon);
      }
      // w = v[vip]^* . v[vim]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vim], Npsivim, &Mab);
      // v[vip] -= w^* . v[vim]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vim], Npsivim, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeCC for the function ComplexBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vim][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vim],  Npsivim, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip],  Npsi,   epsilon);
      }
      // Lowdin orthonormalize v[vip] into v[3], free v[vip] and set v[vip] equal to v[3]
      NSingular=Npsi;
      ComplexWaveFunctionListSVDOrthogonalizeSameBasis(v[vip], v[3], &Npsi, epsilon, true);
      NSingular-=Npsi;
      if(NSingular!=0)
      {
        VerbosityPrintf(0x0001, "Reduce dimension from %u to %u in block %u\n",Npsi+NSingular,Npsi,k);
      }
      for(i=0;i<Npsi+NSingular;i++)
      {
        ComplexWaveFunctionFreeValues(&v[vip][i]);
      }
      for(i=0;i<Npsi;i++)
      {
        v[vip][i]=v[3][i];
      }
      // remove small values from v[vip]
      ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi, epsilon);
    }
    // remove vim from memory and set vim to vi and vi to vip
    if(k==1)
    {
      // nothing to free
    }
    else if(k==2)
    {
      for(i=0;i<Npsi;i++)
      {
        ComplexWaveFunctionFreeValues(&(v[vim][i]));
      }
    }
    else
    {
      BlockWaveFunctionFree(v[vim], Npsivim);
    }
    vim=vi;
    vi=vip;
    vip=(vip+1)%3;
    Npsivim=Npsivi;
    Npsivi=Npsi;
  }
  if(M->N<=1)
  {
    for(i=0;i<Npsivi;i++)
    {
      ComplexWaveFunctionFreeValues(&(v[1][i]));
    }
  }
  else if(M->N==2)
  {
    for(i=0;i<Npsivim;i++)
    {
      ComplexWaveFunctionFreeValues(&(v[vim][i]));
    }
    BlockWaveFunctionFree(v[vi],  Npsivi);
  }
  else
  {
    BlockWaveFunctionFree(v[vim], Npsivim);
    BlockWaveFunctionFree(v[vi],  Npsivi);
  }
  for(i=0;i<Sqr(M->BlockSize);i++)
  {
    M->a[i].Re=0;
    M->a[i].Im=0;
  }
  free(v[1]);
  free(v[2]);
  free(v[3]);
  return 0;
}
uint32_t LanczosBlockTriDiagonalizeRR(const OperatorType H, WaveFunctionType * const psi, double epsilon, BlockTriDiagonalMatrixType * const M, const double Gamma)
{
  uint32_t i,j,k,vi,vip,vim;
  CompactMatrixType Mab;
  WaveFunctionType *v[4];
  uint32_t NSingular,Npsi=M->BlockSize,Npsivi=M->BlockSize,Npsivim=M->BlockSize;
  VerbosityPrintf(0x0001,"Start of LanczosBlockTriDiagonalizeRR\n");
  fflush(stdout);
  for(i=0;i<M->BlockSize;i++)
  {
    M->bRe[i+i*M->BlockSize]=0;
  }
  v[0] = psi;
  v[1] = (WaveFunctionType*)malloc(Npsi*sizeof(WaveFunctionType));
  v[2] = (WaveFunctionType*)malloc(Npsi*sizeof(WaveFunctionType));
  v[3] = (WaveFunctionType*)malloc(Npsi*sizeof(WaveFunctionType));
  Mab.Ni=M->BlockSize;
  Mab.Nj=M->BlockSize;
  Mab.QComplex=0;
  RealWaveFunctionListSVDOrthogonalizeSameBasis(v[0], v[1], &Npsi, epsilon, true); // v[1] = orthonormal[v0]
  NSingular=M->BlockSize-Npsi;
  if(NSingular!=0)
  {
    VerbosityPrintf(0x0001,"Reduce dimension from %u to %u in block %u\n",M->BlockSize,M->BlockSize-NSingular,0);
  }
  // b[0] = v[0] . v[1]^+ or b[0]^T = v[1]^* . v[0]^T
  for(i=0;i<Npsi+NSingular;i++)
  {
    for(j=0;j<Npsi;j++)
    {
      M->bRe[i * M->BlockSize + j] = RealWaveFunctionDotProductSameBasis(v[1][j], v[0][i]);
    }
  }
  vi=1;
  vim=0;
  vip=2;
  Npsivi=Npsi;
  VerbosityPrintf(0x0100,"BasisSetSize: %12i ", v[1][0].N);
  for(k=1;k<M->N;k++)
  {
    // v[vip] = H v[vi]^T
    while(BlockOperatorPsiRR(H, v[vi], v[vip], Npsivi))
    {
      printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRR for the function BlockOperatorPsiRR\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      if (v[vi][0].NBlok == 1) AllocFailed();
      RealBlockWaveFunctionTruncateBasis(v[vi], Npsivi, epsilon);
    }
    VerbosityPrintf(0x0100,"%12i ", v[vi][0].N);
    VerbosityPrintf(0x0100,"| M.a[%3i] M.b[%3i]: %12.5E %12.5E | Err: %12.5E\n", i, i - 1, M->a[i], M->b[i - 1]);
    // a[k] = v[vi]^* . v[vip]^T  (v[vi]^* . H . v[vi]^T)
    Mab.VRe = &(M->aRe[k*Sqr(M->BlockSize)]);
    RealBlockWaveFunctionDotProduct(v[vi], Npsi, v[vip], Npsi, &Mab); //Npsivi = Npsi - changes at Lowdin Orthonormalization
    if(k!=1)
    {
      // b[k-1] = v[vim]^* . v[vip]^T  (v[vim]^* . H . v[vi]^T)
      Mab.VRe = &(M->bRe[(k-1)*Sqr(M->BlockSize)]);
      RealBlockWaveFunctionDotProduct(v[vim], Npsivim, v[vip], Npsi, &Mab);
    }
    NSingular=1;
    while(NSingular)
    {
      // orthogonalize v[vip] to v[vi] and v[vim]
      // we use matrix multiplication between v[vip] and v[vi] which is a numerical unstable form of Gram Schmidt
      // we heal this by repeating the orthogonalization twice for each matrix vector and once for vi and vim (factor of 3 reduction in speed but needed in tests)
      Mab.VRe=M->aRe;
      // w = v[vip]^* . v[vi]^T
      CompactMatrixScale(&Mab, 0);
      RealBlockWaveFunctionDotProduct(v[vip], Npsi, v[vi], Npsivi, &Mab);
      // v[vip] -= w^* . v[vi]
      while(RealBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vi], Npsivi, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRR for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        RealBlockWaveFunctionTruncateBasis(v[vi],  Npsivi, epsilon);
        RealBlockWaveFunctionTruncateBasis(v[vip], Npsi,   epsilon);
      }
      // w = v[vip]^* . v[vi]^T
      CompactMatrixScale(&Mab, 0);
      RealBlockWaveFunctionDotProduct(v[vip], Npsi, v[vi], Npsivi, &Mab);
      // v[vip] -= w^* . v[vi]
      while(RealBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vi], Npsivi, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRR for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        RealBlockWaveFunctionTruncateBasis(v[vi],  Npsivi, epsilon);
        RealBlockWaveFunctionTruncateBasis(v[vip], Npsi,   epsilon);
      }
      // w = v[vip]^* . v[vim]^T
      CompactMatrixScale(&Mab, 0);
      RealBlockWaveFunctionDotProduct(v[vip], Npsi, v[vim], Npsivim, &Mab);
      // v[vip] -= w^* . v[vim]
      while(RealBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vim], Npsivim, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRR for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vim][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        RealBlockWaveFunctionTruncateBasis(v[vim], Npsivim, epsilon);
        RealBlockWaveFunctionTruncateBasis(v[vip], Npsi,    epsilon);
      }
      // w = v[vip]^* . v[vim]^T
      CompactMatrixScale(&Mab, 0);
      RealBlockWaveFunctionDotProduct(v[vip], Npsi, v[vim], Npsivim, &Mab);
      // v[vip] -= w^* . v[vim]
      while(RealBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vim], Npsivim, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRR for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vim][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        RealBlockWaveFunctionTruncateBasis(v[vim], Npsivim, epsilon);
        RealBlockWaveFunctionTruncateBasis(v[vip], Npsi,    epsilon);
      }
      // w = v[vip]^* . v[vi]^T
      CompactMatrixScale(&Mab, 0);
      RealBlockWaveFunctionDotProduct(v[vip], Npsi, v[vi], Npsivi, &Mab);
      // v[vip] -= w^* . v[vi]
      while(RealBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vi], Npsivi, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRR for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        RealBlockWaveFunctionTruncateBasis(v[vi],  Npsivi, epsilon);
        RealBlockWaveFunctionTruncateBasis(v[vip], Npsi,   epsilon);
      }
      // w = v[vip]^* . v[vim]^T
      CompactMatrixScale(&Mab, 0);
      RealBlockWaveFunctionDotProduct(v[vip], Npsi, v[vim], Npsivim, &Mab);
      // v[vip] -= w^* . v[vim]
      while(RealBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vim], Npsivim, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRR for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vim][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        RealBlockWaveFunctionTruncateBasis(v[vim], Npsivim, epsilon);
        RealBlockWaveFunctionTruncateBasis(v[vip], Npsi,    epsilon);
      }
      // Lowdin orthonormalize v[vip] into v[3], free v[vip] and set v[vip] equal to v[3]
      NSingular=Npsi;
      RealWaveFunctionListSVDOrthogonalizeSameBasis(v[vip], v[3], &Npsi, epsilon, true);
      NSingular-=Npsi;
      if(NSingular!=0)
      {
        VerbosityPrintf(0x0001, "Reduce dimension from %u to %u in block %u\n",Npsi+NSingular,Npsi,k);
      }
      for(i=0;i<Npsi+NSingular;i++)
      {
        WaveFunctionFreeValues(&v[vip][i]);
      }
      for(i=0;i<Npsi;i++)
      {
        v[vip][i]=v[3][i];
      }
      // remove small values from v[vip]
      BlockWaveFunctionTruncateBasis(v[vip], Npsi, epsilon);
    }
    // remove vim from memory and set vim to vi and vi to vip
    if(k==1)
    {
      // nothing to free
    }
    else if(k==2)
    {
      for(i=0;i<Npsi;i++)
      {
        RealWaveFunctionFreeValues(&(v[vim][i]));
      }
    }
    else
    {
      BlockWaveFunctionFree(v[vim], Npsivim);
    }
    vim=vi;
    vi=vip;
    vip=(vip+1)%3;
    Npsivim=Npsivi;
    Npsivi=Npsi;
  }
  if(M->N<=1)
  {
    for(i=0;i<Npsivi;i++)
    {
      RealWaveFunctionFreeValues(&(v[1][i]));
    }
  }
  else if(M->N==2)
  {
    for(i=0;i<Npsivim;i++)
    {
      RealWaveFunctionFreeValues(&(v[vim][i]));
    }
    BlockWaveFunctionFree(v[vi],  Npsivi);
  }
  else
  {
    BlockWaveFunctionFree(v[vim], Npsivim);
    BlockWaveFunctionFree(v[vi],  Npsivi);
  }
  for(i=0;i<Sqr(M->BlockSize);i++)
  {
    M->aRe[i]=0;
  }
  free(v[1]);
  free(v[2]);
  free(v[3]);
  return 0;
}
uint32_t LanczosBlockTriDiagonalizeMR(const OperatorType H, WaveFunctionType * const psi, double epsilon, BlockTriDiagonalMatrixType * const M, const double Gamma)
{
  // make psi complex, call MC
  uint32_t i;
  VerbosityPrintf(0x0001,"Start of LanczosBlockTriDiagonalizeMR\n");
  fflush(stdout);
  for(i=0;i<M->BlockSize;i++)
  {
    if (MakeComplexWaveFunctionFromReal(&psi[i], psi[i]))
    {
      printf("LanczosBlockTriDiagonalizeMR called with a real wave-function and a complex Hamiltonian.\nTurned the wave-function complex in step 1.\nBut failed at function %i.\n",i);
      fflush(stdout);
      for(;i>0;i--)
      {
        MakeRealWaveFunctionFromComplex(psi[i-1], &psi[i-1]);
      }
      return 1;
    }
  }
  return LanczosBlockTriDiagonalizeMC(H,psi,epsilon,M,Gamma);
}
uint32_t LanczosBlockTriDiagonalizeCR(const OperatorType H, WaveFunctionType * const psi, double epsilon, BlockTriDiagonalMatrixType * const M, const double Gamma)
{
  // make psi complex, and call CC
  uint32_t i;
  VerbosityPrintf(0x0001,"Start of LanczosBlockTriDiagonalizeCR\n");
  fflush(stdout);
  for(i=0;i<M->BlockSize;i++)
  {
    if (MakeComplexWaveFunctionFromReal(&psi[i], psi[i]))
    {
      printf("LanczosBlockTriDiagonalizeCR called with a real wave-function and a complex Hamiltonian.\nTurned the wave-function complex in step 1.\nBut failed at function %i.\n",i);
      fflush(stdout);
      for(;i>0;i--)
      {
        MakeRealWaveFunctionFromComplex(psi[i-1], &psi[i-1]);
      }
      return 1;
    }
  }
  return LanczosBlockTriDiagonalizeCC(H,psi,epsilon,M,Gamma);
}



uint32_t LanczosBlockTriDiagonalizeKrylovRestrictedRC(const OperatorType H, WaveFunctionType * const psi, DeterminantRestrictionType Restrictions, double epsilon, BlockTriDiagonalMatrixType * const M, WaveFunctionType ** psiBasis, uint32_t *BlockSize)
{
  uint32_t i,j,k,vi,vip,vim;
  CompactMatrixType Mab;
  WaveFunctionType *v[4];
  uint32_t NSingular,Npsi=M->BlockSize,Npsivi=M->BlockSize,Npsivim=M->BlockSize;
  VerbosityPrintf(0x0001,"Start of LanczosBlockTriDiagonalizeRC\n");
  fflush(stdout);
  for(i=0;i<M->BlockSize;i++)
  {
    M->b[i+i*M->BlockSize].Re=0;
    M->b[i+i*M->BlockSize].Im=0;
  }
  v[0] = psi;
  v[3] = (WaveFunctionType*)malloc(Npsi*sizeof(WaveFunctionType));
  
  Mab.Ni=M->BlockSize;
  Mab.Nj=M->BlockSize;
  Mab.QComplex=1;
  v[1] = (WaveFunctionType*)malloc(Npsi*sizeof(WaveFunctionType));
  ComplexWaveFunctionListSVDOrthogonalize(v[0], v[1], &Npsi, epsilon, true); // v[1] = orthonormal[v0]
  BlockSize[0]=Npsi;
  psiBasis[0]=v[1];
  if(BlockSize[0]==0)
  {
    free(v[1]);
    psiBasis[0]=NULL;
    for(i=0;i<Sqr(M->BlockSize);i++)
    {
      M->a[i].Re=0;
      M->b[i].Re=0;
      M->a[i].Im=0;
      M->b[i].Im=0;
    }
    free(v[3]);
    return 0;
  }
  NSingular=M->BlockSize-Npsi;
  if(NSingular!=0)
  {
    VerbosityPrintf(0x0001,"Reduce dimension from %u to %u in block %u\n",M->BlockSize,M->BlockSize-NSingular,0);
  }
  // b[0] = v[0] . v[1]^+ or b[0]^T = v[1]^* . v[0]^T
  for(i=0;i<Npsi+NSingular;i++)
  {
    for(j=0;j<Npsi;j++)
    {
      ComplexWaveFunctionDotProduct(&(v[1][j]), &(v[0][i]), &(M->b[i * M->BlockSize + j].Re), &(M->b[i * M->BlockSize + j].Im));
    }
  }
  vi=1;
  vim=0;
  vip=2;
  Npsivi=Npsi;
  for(k=1;k<M->N;k++)
  {
    // v[vip] = H v[vi]^T
    v[vip] = (WaveFunctionType*)malloc(Npsivi*sizeof(WaveFunctionType));
    while(BlockOperatorPsiRestrictedRC(H, v[vi], v[vip], Restrictions, Npsivi))
    {
      printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRC for the function BlockOperatorPsiRR\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      if (v[vi][0].NBlok == 1) AllocFailed();
      ComplexBlockWaveFunctionTruncateBasis(v[vi], Npsivi, epsilon);
    }
    // a[k] = v[vi]^* . v[vip]^T  (v[vi]^* . H . v[vi]^T)
    Mab.V = &(M->a[k*Sqr(M->BlockSize)]);
    ComplexBlockWaveFunctionDotProduct(v[vi], Npsi, v[vip], Npsi, &Mab); //Npsivi = Npsi - changes at Lowdin Orthonormalization
    if(k!=1)
    {
      // b[k-1] = v[vim]^* . v[vip]^T  (v[vim]^* . H . v[vi]^T)
      Mab.V = &(M->b[(k-1)*Sqr(M->BlockSize)]);
      ComplexBlockWaveFunctionDotProduct(v[vim], Npsivim, v[vip], Npsi, &Mab);
    }
    NSingular=1;
    while(NSingular)
    {
      // orthogonalize v[vip] to v[vi] and v[vim]
      // we use matrix multiplication between v[vip] and v[vi] which is a numerical unstable form of Gram Schmidt
      // we heal this by repeating the orthogonalization twice for each matrix vector and once for vi and vim (factor of 3 reduction in speed but needed in tests)
      Mab.V=M->a;
      // w = v[vip]^* . v[vi]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vi], Npsivi, &Mab);
      // v[vip] -= w^* . v[vi]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vi], Npsivi, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRC for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vi],  Npsivi, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,   epsilon);
      }
      // w = v[vip]^* . v[vi]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vi], Npsivi, &Mab);
      // v[vip] -= w^* . v[vi]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vi], Npsivi, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRC for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vi],  Npsivi, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,   epsilon);
      }
      // w = v[vip]^* . v[vim]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vim], Npsivim, &Mab);
      // v[vip] -= w^* . v[vim]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vim], Npsivim, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRC for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vim][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vim], Npsivim, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,    epsilon);
      }
      // w = v[vip]^* . v[vim]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vim], Npsivim, &Mab);
      // v[vip] -= w^* . v[vim]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vim], Npsivim, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRC for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vim][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vim], Npsivim, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,    epsilon);
      }
      // w = v[vip]^* . v[vi]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vi], Npsivi, &Mab);
      // v[vip] -= w^* . v[vi]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vi], Npsivi, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRC for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vi],  Npsivi, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,   epsilon);
      }
      // w = v[vip]^* . v[vim]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vim], Npsivim, &Mab);
      // v[vip] -= w^* . v[vim]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vim], Npsivim, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRC for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vim][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vim], Npsivim, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,    epsilon);
      }
      // Lowdin orthonormalize v[vip] into v[3], free v[vip] and set v[vip] equal to v[3]
      NSingular=Npsi;
      ComplexWaveFunctionListSVDOrthogonalizeSameBasis(v[vip], v[3], &Npsi, epsilon, true);
      NSingular-=Npsi;
      if(NSingular!=0)
      {
        VerbosityPrintf(0x0001, "Reduce dimension from %u to %u in block %u\n",Npsi+NSingular,Npsi,k);
      }
      for(i=0;i<Npsi+NSingular;i++)
      {
        ComplexWaveFunctionFreeValues(&v[vip][i]);
      }
      for(i=0;i<Npsi;i++)
      {
        v[vip][i]=v[3][i];
      }
      // remove small values from v[vip]
      ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi, epsilon);
    }
    BlockSize[k]=Npsi;
    psiBasis[k]=v[vip];
    if(BlockSize[k]==0)
    {
      free(v[vip]);
      psiBasis[k]=NULL;
      for(i=0;i<Sqr(M->BlockSize);i++)
      {
        M->a[i].Re=0;
        M->a[i].Im=0;
      }
      free(v[3]);
      return 0;
    }
    // set pointers for next state
    vim=vi;
    vi=vip;
    vip=(vip+1)%3;
    Npsivim=Npsivi;
    Npsivi=Npsi;
  }
  for(i=0;i<Sqr(M->BlockSize);i++)
  {
    M->a[i].Re=0;
    M->a[i].Im=0;
  }
  free(v[3]);
  return 0;
}
uint32_t LanczosBlockTriDiagonalizeKrylovRestrictedMC(const OperatorType H, WaveFunctionType * const psi, DeterminantRestrictionType Restrictions, double epsilon, BlockTriDiagonalMatrixType * const M, WaveFunctionType ** psiBasis, uint32_t *BlockSize)
{
  uint32_t i,j,k,vi,vip,vim;
  CompactMatrixType Mab;
  WaveFunctionType *v[4];
  uint32_t NSingular,Npsi=M->BlockSize,Npsivi=M->BlockSize,Npsivim=M->BlockSize;
  VerbosityPrintf(0x0001,"Start of LanczosBlockTriDiagonalizeMC\n");
  fflush(stdout);
  for(i=0;i<M->BlockSize;i++)
  {
    M->b[i+i*M->BlockSize].Re=0;
    M->b[i+i*M->BlockSize].Im=0;
  }
  v[0] = psi;
  v[3] = (WaveFunctionType*)malloc(Npsi*sizeof(WaveFunctionType));
  
  Mab.Ni=M->BlockSize;
  Mab.Nj=M->BlockSize;
  Mab.QComplex=1;
  v[1] = (WaveFunctionType*)malloc(Npsi*sizeof(WaveFunctionType));
  ComplexWaveFunctionListSVDOrthogonalize(v[0], v[1], &Npsi, epsilon, true); // v[1] = orthonormal[v0]
  BlockSize[0]=Npsi;
  psiBasis[0]=v[1];
  if(BlockSize[0]==0)
  {
    free(v[1]);
    psiBasis[0]=NULL;
    for(i=0;i<Sqr(M->BlockSize);i++)
    {
      M->a[i].Re=0;
      M->b[i].Re=0;
      M->a[i].Im=0;
      M->b[i].Im=0;
    }
    free(v[3]);
    return 0;
  }
  NSingular=M->BlockSize-Npsi;
  if(NSingular!=0)
  {
    VerbosityPrintf(0x0001,"Reduce dimension from %u to %u in block %u\n",M->BlockSize,M->BlockSize-NSingular,0);
  }
  // b[0] = v[0] . v[1]^+ or b[0]^T = v[1]^* . v[0]^T
  for(i=0;i<Npsi+NSingular;i++)
  {
    for(j=0;j<Npsi;j++)
    {
      ComplexWaveFunctionDotProduct(&(v[1][j]), &(v[0][i]), &(M->b[i * M->BlockSize + j].Re), &(M->b[i * M->BlockSize + j].Im));
    }
  }
  vi=1;
  vim=0;
  vip=2;
  Npsivi=Npsi;
  for(k=1;k<M->N;k++)
  {
    // v[vip] = H v[vi]^T
    v[vip] = (WaveFunctionType*)malloc(Npsivi*sizeof(WaveFunctionType));
    while(BlockOperatorPsiRestrictedMC(H, v[vi], v[vip], Restrictions, Npsivi))
    {
      printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeMC for the function BlockOperatorPsiRR\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      if (v[vi][0].NBlok == 1) AllocFailed();
      ComplexBlockWaveFunctionTruncateBasis(v[vi], Npsivi, epsilon);
    }
    // a[k] = v[vi]^* . v[vip]^T  (v[vi]^* . H . v[vi]^T)
    Mab.V = &(M->a[k*Sqr(M->BlockSize)]);
    ComplexBlockWaveFunctionDotProduct(v[vi], Npsi, v[vip], Npsi, &Mab); //Npsivi = Npsi - changes at Lowdin Orthonormalization
    if(k!=1)
    {
      // b[k-1] = v[vim]^* . v[vip]^T  (v[vim]^* . H . v[vi]^T)
      Mab.V = &(M->b[(k-1)*Sqr(M->BlockSize)]);
      ComplexBlockWaveFunctionDotProduct(v[vim], Npsivim, v[vip], Npsi, &Mab);
    }
    NSingular=1;
    while(NSingular)
    {
      // orthogonalize v[vip] to v[vi] and v[vim]
      // we use matrix multiplication between v[vip] and v[vi] which is a numerical unstable form of Gram Schmidt
      // we heal this by repeating the orthogonalization twice for each matrix vector and once for vi and vim (factor of 3 reduction in speed but needed in tests)
      Mab.V=M->a;
      // w = v[vip]^* . v[vi]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vi], Npsivi, &Mab);
      // v[vip] -= w^* . v[vi]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vi], Npsivi, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeMC for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vi],  Npsivi, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,   epsilon);
      }
      // w = v[vip]^* . v[vi]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vi], Npsivi, &Mab);
      // v[vip] -= w^* . v[vi]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vi], Npsivi, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeMC for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vi],  Npsivi, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,   epsilon);
      }
      // w = v[vip]^* . v[vim]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vim], Npsivim, &Mab);
      // v[vip] -= w^* . v[vim]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vim], Npsivim, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeMC for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vim][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vim], Npsivim, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,    epsilon);
      }
      // w = v[vip]^* . v[vim]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vim], Npsivim, &Mab);
      // v[vip] -= w^* . v[vim]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vim], Npsivim, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeMC for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vim][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vim], Npsivim, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,    epsilon);
      }
      // w = v[vip]^* . v[vi]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vi], Npsivi, &Mab);
      // v[vip] -= w^* . v[vi]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vi], Npsivi, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeMC for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vi],  Npsivi, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,   epsilon);
      }
      // w = v[vip]^* . v[vim]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vim], Npsivim, &Mab);
      // v[vip] -= w^* . v[vim]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vim], Npsivim, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeMC for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vim][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vim], Npsivim, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,    epsilon);
      }
      // Lowdin orthonormalize v[vip] into v[3], free v[vip] and set v[vip] equal to v[3]
      NSingular=Npsi;
      ComplexWaveFunctionListSVDOrthogonalizeSameBasis(v[vip], v[3], &Npsi, epsilon, true);
      NSingular-=Npsi;
      if(NSingular!=0)
      {
        VerbosityPrintf(0x0001, "Reduce dimension from %u to %u in block %u\n",Npsi+NSingular,Npsi,k);
      }
      for(i=0;i<Npsi+NSingular;i++)
      {
        ComplexWaveFunctionFreeValues(&v[vip][i]);
      }
      for(i=0;i<Npsi;i++)
      {
        v[vip][i]=v[3][i];
      }
      // remove small values from v[vip]
      ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi, epsilon);
    }
    BlockSize[k]=Npsi;
    psiBasis[k]=v[vip];
    if(BlockSize[k]==0)
    {
      free(v[vip]);
      psiBasis[k]=NULL;
      for(i=0;i<Sqr(M->BlockSize);i++)
      {
        M->a[i].Re=0;
        M->a[i].Im=0;
      }
      free(v[3]);
      return 0;
    }
    // set pointers for next state
    vim=vi;
    vi=vip;
    vip=(vip+1)%3;
    Npsivim=Npsivi;
    Npsivi=Npsi;
  }
  for(i=0;i<Sqr(M->BlockSize);i++)
  {
    M->a[i].Re=0;
    M->a[i].Im=0;
  }
  free(v[3]);
  return 0;
}
uint32_t LanczosBlockTriDiagonalizeKrylovRestrictedCC(const OperatorType H, WaveFunctionType * const psi, DeterminantRestrictionType Restrictions, double epsilon, BlockTriDiagonalMatrixType * const M, WaveFunctionType ** psiBasis, uint32_t *BlockSize)
{
  uint32_t i,j,k,vi,vip,vim;
  CompactMatrixType Mab;
  WaveFunctionType *v[4];
  uint32_t NSingular,Npsi=M->BlockSize,Npsivi=M->BlockSize,Npsivim=M->BlockSize;
  VerbosityPrintf(0x0001,"Start of LanczosBlockTriDiagonalizeCC\n");
  fflush(stdout);
  for(i=0;i<M->BlockSize;i++)
  {
    M->b[i+i*M->BlockSize].Re=0;
    M->b[i+i*M->BlockSize].Im=0;
  }
  v[0] = psi;
  v[3] = (WaveFunctionType*)malloc(Npsi*sizeof(WaveFunctionType));
  
  Mab.Ni=M->BlockSize;
  Mab.Nj=M->BlockSize;
  Mab.QComplex=1;
  v[1] = (WaveFunctionType*)malloc(Npsi*sizeof(WaveFunctionType));
  ComplexWaveFunctionListSVDOrthogonalize(v[0], v[1], &Npsi, epsilon, true); // v[1] = orthonormal[v0]
  BlockSize[0]=Npsi;
  psiBasis[0]=v[1];
  if(BlockSize[0]==0)
  {
    free(v[1]);
    psiBasis[0]=NULL;
    for(i=0;i<Sqr(M->BlockSize);i++)
    {
      M->a[i].Re=0;
      M->b[i].Re=0;
      M->a[i].Im=0;
      M->b[i].Im=0;
    }
    free(v[3]);
    return 0;
  }
  NSingular=M->BlockSize-Npsi;
  if(NSingular!=0)
  {
    VerbosityPrintf(0x0001,"Reduce dimension from %u to %u in block %u\n",M->BlockSize,M->BlockSize-NSingular,0);
  }
  // b[0] = v[0] . v[1]^+ or b[0]^T = v[1]^* . v[0]^T
  for(i=0;i<Npsi+NSingular;i++)
  {
    for(j=0;j<Npsi;j++)
    {
      ComplexWaveFunctionDotProduct(&(v[1][j]), &(v[0][i]), &(M->b[i * M->BlockSize + j].Re), &(M->b[i * M->BlockSize + j].Im));
    }
  }
  vi=1;
  vim=0;
  vip=2;
  Npsivi=Npsi;
  for(k=1;k<M->N;k++)
  {
    // v[vip] = H v[vi]^T
    v[vip] = (WaveFunctionType*)malloc(Npsivi*sizeof(WaveFunctionType));
    while(BlockOperatorPsiRestrictedCC(H, v[vi], v[vip], Restrictions, Npsivi))
    {
      printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeCC for the function BlockOperatorPsiRR\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      if (v[vi][0].NBlok == 1) AllocFailed();
      ComplexBlockWaveFunctionTruncateBasis(v[vi], Npsivi, epsilon);
    }
    // a[k] = v[vi]^* . v[vip]^T  (v[vi]^* . H . v[vi]^T)
    Mab.V = &(M->a[k*Sqr(M->BlockSize)]);
    ComplexBlockWaveFunctionDotProduct(v[vi], Npsi, v[vip], Npsi, &Mab); //Npsivi = Npsi - changes at Lowdin Orthonormalization
    if(k!=1)
    {
      // b[k-1] = v[vim]^* . v[vip]^T  (v[vim]^* . H . v[vi]^T)
      Mab.V = &(M->b[(k-1)*Sqr(M->BlockSize)]);
      ComplexBlockWaveFunctionDotProduct(v[vim], Npsivim, v[vip], Npsi, &Mab);
    }
    NSingular=1;
    while(NSingular)
    {
      // orthogonalize v[vip] to v[vi] and v[vim]
      // we use matrix multiplication between v[vip] and v[vi] which is a numerical unstable form of Gram Schmidt
      // we heal this by repeating the orthogonalization twice for each matrix vector and once for vi and vim (factor of 3 reduction in speed but needed in tests)
      Mab.V=M->a;
      // w = v[vip]^* . v[vi]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vi], Npsivi, &Mab);
      // v[vip] -= w^* . v[vi]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vi], Npsivi, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeCC for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vi],  Npsivi, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,   epsilon);
      }
      // w = v[vip]^* . v[vi]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vi], Npsivi, &Mab);
      // v[vip] -= w^* . v[vi]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vi], Npsivi, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeCC for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vi],  Npsivi, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,   epsilon);
      }
      // w = v[vip]^* . v[vim]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vim], Npsivim, &Mab);
      // v[vip] -= w^* . v[vim]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vim], Npsivim, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeCC for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vim][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vim], Npsivim, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,    epsilon);
      }
      // w = v[vip]^* . v[vim]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vim], Npsivim, &Mab);
      // v[vip] -= w^* . v[vim]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vim], Npsivim, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeCC for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vim][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vim], Npsivim, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,    epsilon);
      }
      // w = v[vip]^* . v[vi]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vi], Npsivi, &Mab);
      // v[vip] -= w^* . v[vi]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vi], Npsivi, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeCC for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vi],  Npsivi, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,   epsilon);
      }
      // w = v[vip]^* . v[vim]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vim], Npsivim, &Mab);
      // v[vip] -= w^* . v[vim]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vim], Npsivim, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeCC for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vim][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vim], Npsivim, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,    epsilon);
      }
      // Lowdin orthonormalize v[vip] into v[3], free v[vip] and set v[vip] equal to v[3]
      NSingular=Npsi;
      ComplexWaveFunctionListSVDOrthogonalizeSameBasis(v[vip], v[3], &Npsi, epsilon, true);
      NSingular-=Npsi;
      if(NSingular!=0)
      {
        VerbosityPrintf(0x0001, "Reduce dimension from %u to %u in block %u\n",Npsi+NSingular,Npsi,k);
      }
      for(i=0;i<Npsi+NSingular;i++)
      {
        ComplexWaveFunctionFreeValues(&v[vip][i]);
      }
      for(i=0;i<Npsi;i++)
      {
        v[vip][i]=v[3][i];
      }
      // remove small values from v[vip]
      ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi, epsilon);
    }
    BlockSize[k]=Npsi;
    psiBasis[k]=v[vip];
    if(BlockSize[k]==0)
    {
      free(v[vip]);
      psiBasis[k]=NULL;
      for(i=0;i<Sqr(M->BlockSize);i++)
      {
        M->a[i].Re=0;
        M->a[i].Im=0;
      }
      free(v[3]);
      return 0;
    }
    // set pointers for next state
    vim=vi;
    vi=vip;
    vip=(vip+1)%3;
    Npsivim=Npsivi;
    Npsivi=Npsi;
  }
  for(i=0;i<Sqr(M->BlockSize);i++)
  {
    M->a[i].Re=0;
    M->a[i].Im=0;
  }
  free(v[3]);
  return 0;
}
uint32_t LanczosBlockTriDiagonalizeKrylovRestrictedRR(const OperatorType H, WaveFunctionType * const psi, DeterminantRestrictionType Restrictions, double epsilon, BlockTriDiagonalMatrixType * const M, WaveFunctionType ** psiBasis, uint32_t *BlockSize)
{
  uint32_t i,j,k,vi,vip,vim;
  CompactMatrixType Mab;
  WaveFunctionType *v[4];
  uint32_t NSingular,Npsi=M->BlockSize,Npsivi=M->BlockSize,Npsivim=M->BlockSize;
  VerbosityPrintf(0x0001,"Start of LanczosBlockTriDiagonalizeKrylovRestrictedRR\n");
  fflush(stdout);
  for(i=0;i<M->BlockSize;i++)
  {
    M->bRe[i+i*M->BlockSize]=0;
  }
  v[0] = psi;
  v[3] = (WaveFunctionType*)malloc(Npsi*sizeof(WaveFunctionType));
  
  Mab.Ni=M->BlockSize;
  Mab.Nj=M->BlockSize;
  Mab.QComplex=0;
  v[1] = (WaveFunctionType*)malloc(Npsi*sizeof(WaveFunctionType));
  RealWaveFunctionListSVDOrthogonalize(v[0], v[1], &Npsi, epsilon, true); // v[1] = orthonormal[v0]
  BlockSize[0]=Npsi;
  psiBasis[0]=v[1];
  if(BlockSize[0]==0)
  {
    free(v[1]);
    psiBasis[0]=NULL;
    for(i=0;i<Sqr(M->BlockSize);i++)
    {
      M->aRe[i]=0;
      M->bRe[i]=0;
    }
    free(v[3]);
    return 0;
  }
  NSingular=M->BlockSize-Npsi;
  if(NSingular!=0)
  {
    VerbosityPrintf(0x0001,"Reduce dimension from %u to %u in block %u\n",M->BlockSize,M->BlockSize-NSingular,0);
  }
  // b[0] = v[0] . v[1]^+ or b[0]^T = v[1]^* . v[0]^T
  for(i=0;i<Npsi+NSingular;i++)
  {
    for(j=0;j<Npsi;j++)
    {
      M->bRe[i * M->BlockSize + j] = RealWaveFunctionDotProduct(&(v[1][j]), &(v[0][i]));
    }
  }
  vi=1;
  vim=0;
  vip=2;
  Npsivi=Npsi;
  for(k=1;k<M->N;k++)
  {
    // v[vip] = H v[vi]^T
    v[vip] = (WaveFunctionType*)malloc(Npsivi*sizeof(WaveFunctionType));
    while(BlockOperatorPsiRestrictedRR(H, v[vi], v[vip], Restrictions, Npsivi))
    {
      printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRR for the function BlockOperatorPsiRR\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      if (v[vi][0].NBlok == 1) AllocFailed();
      RealBlockWaveFunctionTruncateBasis(v[vi], Npsivi, epsilon);
    }
    // a[k] = v[vi]^* . v[vip]^T  (v[vi]^* . H . v[vi]^T)
    Mab.VRe = &(M->aRe[k*Sqr(M->BlockSize)]);
    RealBlockWaveFunctionDotProduct(v[vi], Npsi, v[vip], Npsi, &Mab); //Npsivi = Npsi - changes at Lowdin Orthonormalization
    if(k!=1)
    {
      // b[k-1] = v[vim]^* . v[vip]^T  (v[vim]^* . H . v[vi]^T)
      Mab.VRe = &(M->bRe[(k-1)*Sqr(M->BlockSize)]);
      RealBlockWaveFunctionDotProduct(v[vim], Npsivim, v[vip], Npsi, &Mab);
    }
    NSingular=1;
    while(NSingular)
    {
      // orthogonalize v[vip] to v[vi] and v[vim]
      // we use matrix multiplication between v[vip] and v[vi] which is a numerical unstable form of Gram Schmidt
      // we heal this by repeating the orthogonalization twice for each matrix vector and once for vi and vim (factor of 3 reduction in speed but needed in tests)
      Mab.VRe=M->aRe;
      // w = v[vip]^* . v[vi]^T
      CompactMatrixScale(&Mab, 0);
      RealBlockWaveFunctionDotProduct(v[vip], Npsi, v[vi], Npsivi, &Mab);
      // v[vip] -= w^* . v[vi]
      while(RealBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vi], Npsivi, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRR for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        RealBlockWaveFunctionTruncateBasis(v[vi],  Npsivi, epsilon);
        RealBlockWaveFunctionTruncateBasis(v[vip], Npsi,   epsilon);
      }
      // w = v[vip]^* . v[vi]^T
      CompactMatrixScale(&Mab, 0);
      RealBlockWaveFunctionDotProduct(v[vip], Npsi, v[vi], Npsivi, &Mab);
      // v[vip] -= w^* . v[vi]
      while(RealBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vi], Npsivi, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRR for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        RealBlockWaveFunctionTruncateBasis(v[vi],  Npsivi, epsilon);
        RealBlockWaveFunctionTruncateBasis(v[vip], Npsi,   epsilon);
      }
      // w = v[vip]^* . v[vim]^T
      CompactMatrixScale(&Mab, 0);
      RealBlockWaveFunctionDotProduct(v[vip], Npsi, v[vim], Npsivim, &Mab);
      // v[vip] -= w^* . v[vim]
      while(RealBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vim], Npsivim, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRR for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vim][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        RealBlockWaveFunctionTruncateBasis(v[vim], Npsivim, epsilon);
        RealBlockWaveFunctionTruncateBasis(v[vip], Npsi,    epsilon);
      }
      // w = v[vip]^* . v[vim]^T
      CompactMatrixScale(&Mab, 0);
      RealBlockWaveFunctionDotProduct(v[vip], Npsi, v[vim], Npsivim, &Mab);
      // v[vip] -= w^* . v[vim]
      while(RealBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vim], Npsivim, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRR for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vim][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        RealBlockWaveFunctionTruncateBasis(v[vim], Npsivim, epsilon);
        RealBlockWaveFunctionTruncateBasis(v[vip], Npsi,    epsilon);
      }
      // w = v[vip]^* . v[vi]^T
      CompactMatrixScale(&Mab, 0);
      RealBlockWaveFunctionDotProduct(v[vip], Npsi, v[vi], Npsivi, &Mab);
      // v[vip] -= w^* . v[vi]
      while(RealBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vi], Npsivi, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRR for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        RealBlockWaveFunctionTruncateBasis(v[vi],  Npsivi, epsilon);
        RealBlockWaveFunctionTruncateBasis(v[vip], Npsi,   epsilon);
      }
      // w = v[vip]^* . v[vim]^T
      CompactMatrixScale(&Mab, 0);
      RealBlockWaveFunctionDotProduct(v[vip], Npsi, v[vim], Npsivim, &Mab);
      // v[vip] -= w^* . v[vim]
      while(RealBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vim], Npsivim, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRR for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vim][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        RealBlockWaveFunctionTruncateBasis(v[vim], Npsivim, epsilon);
        RealBlockWaveFunctionTruncateBasis(v[vip], Npsi,    epsilon);
      }
      // Lowdin orthonormalize v[vip] into v[3], free v[vip] and set v[vip] equal to v[3]
      NSingular=Npsi;
      RealWaveFunctionListSVDOrthogonalizeSameBasis(v[vip], v[3], &Npsi, epsilon, true);
      NSingular-=Npsi;
      if(NSingular!=0)
      {
        VerbosityPrintf(0x0001, "Reduce dimension from %u to %u in block %u\n",Npsi+NSingular,Npsi,k);
      }
      for(i=0;i<Npsi+NSingular;i++)
      {
        WaveFunctionFreeValues(&v[vip][i]);
      }
      for(i=0;i<Npsi;i++)
      {
        v[vip][i]=v[3][i];
      }
      // remove small values from v[vip]
      BlockWaveFunctionTruncateBasis(v[vip], Npsi, epsilon);
    }
    BlockSize[k]=Npsi;
    psiBasis[k]=v[vip];
    if(BlockSize[k]==0)
    {
      free(v[vip]);
      psiBasis[k]=NULL;
      for(i=0;i<Sqr(M->BlockSize);i++)
      {
        M->aRe[i]=0;
      }
      free(v[3]);
      return 0;
    }
    // set pointers for next state
    vim=vi;
    vi=vip;
    vip=(vip+1)%3;
    Npsivim=Npsivi;
    Npsivi=Npsi;
  }
  for(i=0;i<Sqr(M->BlockSize);i++)
  {
    M->aRe[i]=0;
  }
  free(v[3]);
  return 0;
}
uint32_t LanczosBlockTriDiagonalizeKrylovRestrictedMR(const OperatorType H, WaveFunctionType * const psi, DeterminantRestrictionType Restrictions, double epsilon, BlockTriDiagonalMatrixType * const M, WaveFunctionType ** psiBasis, uint32_t *BlockSize)
{
  // make psi complex, call MC
  uint32_t i;
  VerbosityPrintf(0x0001,"Start of LanczosBlockTriDiagonalizeKrylovRestrictedMR\n");
  fflush(stdout);
  for(i=0;i<M->BlockSize;i++)
  {
    if (MakeComplexWaveFunctionFromReal(&psi[i], psi[i]))
    {
      printf("LanczosBlockTriDiagonalizeKrylovRestrictedMR called with a real wave-function and a complex Hamiltonian.\nTurned the wave-function complex in step 1.\nBut failed at function %i.\n",i);
      fflush(stdout);
      for(;i>0;i--)
      {
        MakeRealWaveFunctionFromComplex(psi[i-1], &psi[i-1]);
      }
      return 1;
    }
  }
  return LanczosBlockTriDiagonalizeKrylovRestrictedMC(H,psi,Restrictions,epsilon,M,psiBasis,BlockSize);
}
uint32_t LanczosBlockTriDiagonalizeKrylovRestrictedCR(const OperatorType H, WaveFunctionType * const psi, DeterminantRestrictionType Restrictions, double epsilon, BlockTriDiagonalMatrixType * const M, WaveFunctionType ** psiBasis, uint32_t *BlockSize)
{
  // make psi complex, and call CC
  uint32_t i;
  VerbosityPrintf(0x0001,"Start of LanczosBlockTriDiagonalizeKrylovRestrictedCR\n");
  fflush(stdout);
  for(i=0;i<M->BlockSize;i++)
  {
    if (MakeComplexWaveFunctionFromReal(&psi[i], psi[i]))
    {
      printf("LanczosBlockTriDiagonalizeKrylovRestrictedCR called with a real wave-function and a complex Hamiltonian.\nTurned the wave-function complex in step 1.\nBut failed at function %i.\n",i);
      fflush(stdout);
      for(;i>0;i--)
      {
        MakeRealWaveFunctionFromComplex(psi[i-1], &psi[i-1]);
      }
      return 1;
    }
  }
  return LanczosBlockTriDiagonalizeKrylovRestrictedCC(H,psi,Restrictions,epsilon,M,psiBasis,BlockSize);
}

uint32_t LanczosBlockTriDiagonalizeKrylovRC(const OperatorType H, WaveFunctionType * const psi, double epsilon, BlockTriDiagonalMatrixType * const M, WaveFunctionType ** psiBasis, uint32_t *BlockSize)
{
  uint32_t i,j,k,vi,vip,vim;
  CompactMatrixType Mab;
  WaveFunctionType *v[4];
  uint32_t NSingular,Npsi=M->BlockSize,Npsivi=M->BlockSize,Npsivim=M->BlockSize;
  VerbosityPrintf(0x0001,"Start of LanczosBlockTriDiagonalizeRC\n");
  fflush(stdout);
  for(i=0;i<M->BlockSize;i++)
  {
    M->b[i+i*M->BlockSize].Re=0;
    M->b[i+i*M->BlockSize].Im=0;
  }
  v[0] = psi;
  v[3] = (WaveFunctionType*)malloc(Npsi*sizeof(WaveFunctionType));
  
  Mab.Ni=M->BlockSize;
  Mab.Nj=M->BlockSize;
  Mab.QComplex=1;
  v[1] = (WaveFunctionType*)malloc(Npsi*sizeof(WaveFunctionType));
  ComplexWaveFunctionListSVDOrthogonalize(v[0], v[1], &Npsi, epsilon, true); // v[1] = orthonormal[v0]
  BlockSize[0]=Npsi;
  psiBasis[0]=v[1];
  if(BlockSize[0]==0)
  {
    free(v[1]);
    psiBasis[0]=NULL;
    for(i=0;i<Sqr(M->BlockSize);i++)
    {
      M->a[i].Re=0;
      M->b[i].Re=0;
      M->a[i].Im=0;
      M->b[i].Im=0;
    }
    free(v[3]);
    return 0;
  }
  NSingular=M->BlockSize-Npsi;
  if(NSingular!=0)
  {
    VerbosityPrintf(0x0001,"Reduce dimension from %u to %u in block %u\n",M->BlockSize,M->BlockSize-NSingular,0);
  }
  // b[0] = v[0] . v[1]^+ or b[0]^T = v[1]^* . v[0]^T
  for(i=0;i<Npsi+NSingular;i++)
  {
    for(j=0;j<Npsi;j++)
    {
      ComplexWaveFunctionDotProduct(&(v[1][j]), &(v[0][i]), &(M->b[i * M->BlockSize + j].Re), &(M->b[i * M->BlockSize + j].Im));
    }
  }
  vi=1;
  vim=0;
  vip=2;
  Npsivi=Npsi;
  for(k=1;k<M->N;k++)
  {
    // v[vip] = H v[vi]^T
    v[vip] = (WaveFunctionType*)malloc(Npsivi*sizeof(WaveFunctionType));
    while(BlockOperatorPsiRC(H, v[vi], v[vip], Npsivi))
    {
      printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRC for the function BlockOperatorPsiRR\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      if (v[vi][0].NBlok == 1) AllocFailed();
      ComplexBlockWaveFunctionTruncateBasis(v[vi], Npsivi, epsilon);
    }
    // a[k] = v[vi]^* . v[vip]^T  (v[vi]^* . H . v[vi]^T)
    Mab.V = &(M->a[k*Sqr(M->BlockSize)]);
    ComplexBlockWaveFunctionDotProduct(v[vi], Npsi, v[vip], Npsi, &Mab); //Npsivi = Npsi - changes at Lowdin Orthonormalization
    if(k!=1)
    {
      // b[k-1] = v[vim]^* . v[vip]^T  (v[vim]^* . H . v[vi]^T)
      Mab.V = &(M->b[(k-1)*Sqr(M->BlockSize)]);
      ComplexBlockWaveFunctionDotProduct(v[vim], Npsivim, v[vip], Npsi, &Mab);
    }
    NSingular=1;
    while(NSingular)
    {
      // orthogonalize v[vip] to v[vi] and v[vim]
      // we use matrix multiplication between v[vip] and v[vi] which is a numerical unstable form of Gram Schmidt
      // we heal this by repeating the orthogonalization twice for each matrix vector and once for vi and vim (factor of 3 reduction in speed but needed in tests)
      Mab.V=M->a;
      // w = v[vip]^* . v[vi]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vi], Npsivi, &Mab);
      // v[vip] -= w^* . v[vi]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vi], Npsivi, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRC for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vi],  Npsivi, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,   epsilon);
      }
      // w = v[vip]^* . v[vi]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vi], Npsivi, &Mab);
      // v[vip] -= w^* . v[vi]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vi], Npsivi, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRC for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vi],  Npsivi, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,   epsilon);
      }
      // w = v[vip]^* . v[vim]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vim], Npsivim, &Mab);
      // v[vip] -= w^* . v[vim]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vim], Npsivim, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRC for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vim][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vim], Npsivim, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,    epsilon);
      }
      // w = v[vip]^* . v[vim]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vim], Npsivim, &Mab);
      // v[vip] -= w^* . v[vim]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vim], Npsivim, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRC for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vim][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vim], Npsivim, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,    epsilon);
      }
      // w = v[vip]^* . v[vi]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vi], Npsivi, &Mab);
      // v[vip] -= w^* . v[vi]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vi], Npsivi, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRC for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vi],  Npsivi, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,   epsilon);
      }
      // w = v[vip]^* . v[vim]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vim], Npsivim, &Mab);
      // v[vip] -= w^* . v[vim]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vim], Npsivim, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRC for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vim][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vim], Npsivim, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,    epsilon);
      }
      // Lowdin orthonormalize v[vip] into v[3], free v[vip] and set v[vip] equal to v[3]
      NSingular=Npsi;
      ComplexWaveFunctionListSVDOrthogonalizeSameBasis(v[vip], v[3], &Npsi, epsilon, true);
      NSingular-=Npsi;
      if(NSingular!=0)
      {
        VerbosityPrintf(0x0001, "Reduce dimension from %u to %u in block %u\n",Npsi+NSingular,Npsi,k);
      }
      for(i=0;i<Npsi+NSingular;i++)
      {
        ComplexWaveFunctionFreeValues(&v[vip][i]);
      }
      for(i=0;i<Npsi;i++)
      {
        v[vip][i]=v[3][i];
      }
      // remove small values from v[vip]
      ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi, epsilon);
    }
    BlockSize[k]=Npsi;
    psiBasis[k]=v[vip];
    if(BlockSize[k]==0)
    {
      free(v[vip]);
      psiBasis[k]=NULL;
      for(i=0;i<Sqr(M->BlockSize);i++)
      {
        M->a[i].Re=0;
        M->a[i].Im=0;
      }
      free(v[3]);
      return 0;
    }
    // set pointers for next state
    vim=vi;
    vi=vip;
    vip=(vip+1)%3;
    Npsivim=Npsivi;
    Npsivi=Npsi;
  }
  for(i=0;i<Sqr(M->BlockSize);i++)
  {
    M->a[i].Re=0;
    M->a[i].Im=0;
  }
  free(v[3]);
  return 0;
}
uint32_t LanczosBlockTriDiagonalizeKrylovMC(const OperatorType H, WaveFunctionType * const psi, double epsilon, BlockTriDiagonalMatrixType * const M, WaveFunctionType ** psiBasis, uint32_t *BlockSize)
{
  uint32_t i,j,k,vi,vip,vim;
  CompactMatrixType Mab;
  WaveFunctionType *v[4];
  uint32_t NSingular,Npsi=M->BlockSize,Npsivi=M->BlockSize,Npsivim=M->BlockSize;
  VerbosityPrintf(0x0001,"Start of LanczosBlockTriDiagonalizeMC\n");
  fflush(stdout);
  for(i=0;i<M->BlockSize;i++)
  {
    M->b[i+i*M->BlockSize].Re=0;
    M->b[i+i*M->BlockSize].Im=0;
  }
  v[0] = psi;
  v[3] = (WaveFunctionType*)malloc(Npsi*sizeof(WaveFunctionType));
  
  Mab.Ni=M->BlockSize;
  Mab.Nj=M->BlockSize;
  Mab.QComplex=1;
  v[1] = (WaveFunctionType*)malloc(Npsi*sizeof(WaveFunctionType));
  ComplexWaveFunctionListSVDOrthogonalize(v[0], v[1], &Npsi, epsilon, true); // v[1] = orthonormal[v0]
  BlockSize[0]=Npsi;
  psiBasis[0]=v[1];
  if(BlockSize[0]==0)
  {
    free(v[1]);
    psiBasis[0]=NULL;
    for(i=0;i<Sqr(M->BlockSize);i++)
    {
      M->a[i].Re=0;
      M->b[i].Re=0;
      M->a[i].Im=0;
      M->b[i].Im=0;
    }
    free(v[3]);
    return 0;
  }
  NSingular=M->BlockSize-Npsi;
  if(NSingular!=0)
  {
    VerbosityPrintf(0x0001,"Reduce dimension from %u to %u in block %u\n",M->BlockSize,M->BlockSize-NSingular,0);
  }
  // b[0] = v[0] . v[1]^+ or b[0]^T = v[1]^* . v[0]^T
  for(i=0;i<Npsi+NSingular;i++)
  {
    for(j=0;j<Npsi;j++)
    {
      ComplexWaveFunctionDotProduct(&(v[1][j]), &(v[0][i]), &(M->b[i * M->BlockSize + j].Re), &(M->b[i * M->BlockSize + j].Im));
    }
  }
  vi=1;
  vim=0;
  vip=2;
  Npsivi=Npsi;
  for(k=1;k<M->N;k++)
  {
    // v[vip] = H v[vi]^T
    v[vip] = (WaveFunctionType*)malloc(Npsivi*sizeof(WaveFunctionType));
    while(BlockOperatorPsiMC(H, v[vi], v[vip], Npsivi))
    {
      printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeMC for the function BlockOperatorPsiRR\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      if (v[vi][0].NBlok == 1) AllocFailed();
      ComplexBlockWaveFunctionTruncateBasis(v[vi], Npsivi, epsilon);
    }
    // a[k] = v[vi]^* . v[vip]^T  (v[vi]^* . H . v[vi]^T)
    Mab.V = &(M->a[k*Sqr(M->BlockSize)]);
    ComplexBlockWaveFunctionDotProduct(v[vi], Npsi, v[vip], Npsi, &Mab); //Npsivi = Npsi - changes at Lowdin Orthonormalization
    if(k!=1)
    {
      // b[k-1] = v[vim]^* . v[vip]^T  (v[vim]^* . H . v[vi]^T)
      Mab.V = &(M->b[(k-1)*Sqr(M->BlockSize)]);
      ComplexBlockWaveFunctionDotProduct(v[vim], Npsivim, v[vip], Npsi, &Mab);
    }
    NSingular=1;
    while(NSingular)
    {
      // orthogonalize v[vip] to v[vi] and v[vim]
      // we use matrix multiplication between v[vip] and v[vi] which is a numerical unstable form of Gram Schmidt
      // we heal this by repeating the orthogonalization twice for each matrix vector and once for vi and vim (factor of 3 reduction in speed but needed in tests)
      Mab.V=M->a;
      // w = v[vip]^* . v[vi]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vi], Npsivi, &Mab);
      // v[vip] -= w^* . v[vi]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vi], Npsivi, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeMC for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vi],  Npsivi, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,   epsilon);
      }
      // w = v[vip]^* . v[vi]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vi], Npsivi, &Mab);
      // v[vip] -= w^* . v[vi]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vi], Npsivi, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeMC for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vi],  Npsivi, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,   epsilon);
      }
      // w = v[vip]^* . v[vim]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vim], Npsivim, &Mab);
      // v[vip] -= w^* . v[vim]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vim], Npsivim, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeMC for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vim][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vim], Npsivim, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,    epsilon);
      }
      // w = v[vip]^* . v[vim]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vim], Npsivim, &Mab);
      // v[vip] -= w^* . v[vim]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vim], Npsivim, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeMC for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vim][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vim], Npsivim, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,    epsilon);
      }
      // w = v[vip]^* . v[vi]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vi], Npsivi, &Mab);
      // v[vip] -= w^* . v[vi]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vi], Npsivi, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeMC for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vi],  Npsivi, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,   epsilon);
      }
      // w = v[vip]^* . v[vim]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vim], Npsivim, &Mab);
      // v[vip] -= w^* . v[vim]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vim], Npsivim, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeMC for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vim][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vim], Npsivim, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,    epsilon);
      }
      // Lowdin orthonormalize v[vip] into v[3], free v[vip] and set v[vip] equal to v[3]
      NSingular=Npsi;
      ComplexWaveFunctionListSVDOrthogonalizeSameBasis(v[vip], v[3], &Npsi, epsilon, true);
      NSingular-=Npsi;
      if(NSingular!=0)
      {
        VerbosityPrintf(0x0001, "Reduce dimension from %u to %u in block %u\n",Npsi+NSingular,Npsi,k);
      }
      for(i=0;i<Npsi+NSingular;i++)
      {
        ComplexWaveFunctionFreeValues(&v[vip][i]);
      }
      for(i=0;i<Npsi;i++)
      {
        v[vip][i]=v[3][i];
      }
      // remove small values from v[vip]
      ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi, epsilon);
    }
    BlockSize[k]=Npsi;
    psiBasis[k]=v[vip];
    if(BlockSize[k]==0)
    {
      free(v[vip]);
      psiBasis[k]=NULL;
      for(i=0;i<Sqr(M->BlockSize);i++)
      {
        M->a[i].Re=0;
        M->a[i].Im=0;
      }
      free(v[3]);
      return 0;
    }
    // set pointers for next state
    vim=vi;
    vi=vip;
    vip=(vip+1)%3;
    Npsivim=Npsivi;
    Npsivi=Npsi;
  }
  for(i=0;i<Sqr(M->BlockSize);i++)
  {
    M->a[i].Re=0;
    M->a[i].Im=0;
  }
  free(v[3]);
  return 0;
}
uint32_t LanczosBlockTriDiagonalizeKrylovCC(const OperatorType H, WaveFunctionType * const psi, double epsilon, BlockTriDiagonalMatrixType * const M, WaveFunctionType ** psiBasis, uint32_t *BlockSize)
{
  uint32_t i,j,k,vi,vip,vim;
  CompactMatrixType Mab;
  WaveFunctionType *v[4];
  uint32_t NSingular,Npsi=M->BlockSize,Npsivi=M->BlockSize,Npsivim=M->BlockSize;
  VerbosityPrintf(0x0001,"Start of LanczosBlockTriDiagonalizeCC\n");
  fflush(stdout);
  for(i=0;i<M->BlockSize;i++)
  {
    M->b[i+i*M->BlockSize].Re=0;
    M->b[i+i*M->BlockSize].Im=0;
  }
  v[0] = psi;
  v[3] = (WaveFunctionType*)malloc(Npsi*sizeof(WaveFunctionType));
  
  Mab.Ni=M->BlockSize;
  Mab.Nj=M->BlockSize;
  Mab.QComplex=1;
  v[1] = (WaveFunctionType*)malloc(Npsi*sizeof(WaveFunctionType));
  ComplexWaveFunctionListSVDOrthogonalize(v[0], v[1], &Npsi, epsilon, true); // v[1] = orthonormal[v0]
  BlockSize[0]=Npsi;
  psiBasis[0]=v[1];
  if(BlockSize[0]==0)
  {
    free(v[1]);
    psiBasis[0]=NULL;
    for(i=0;i<Sqr(M->BlockSize);i++)
    {
      M->a[i].Re=0;
      M->b[i].Re=0;
      M->a[i].Im=0;
      M->b[i].Im=0;
    }
    free(v[3]);
    return 0;
  }
  NSingular=M->BlockSize-Npsi;
  if(NSingular!=0)
  {
    VerbosityPrintf(0x0001,"Reduce dimension from %u to %u in block %u\n",M->BlockSize,M->BlockSize-NSingular,0);
  }
  // b[0] = v[0] . v[1]^+ or b[0]^T = v[1]^* . v[0]^T
  for(i=0;i<Npsi+NSingular;i++)
  {
    for(j=0;j<Npsi;j++)
    {
      ComplexWaveFunctionDotProduct(&(v[1][j]), &(v[0][i]), &(M->b[i * M->BlockSize + j].Re), &(M->b[i * M->BlockSize + j].Im));
    }
  }
  vi=1;
  vim=0;
  vip=2;
  Npsivi=Npsi;
  for(k=1;k<M->N;k++)
  {
    // v[vip] = H v[vi]^T
    v[vip] = (WaveFunctionType*)malloc(Npsivi*sizeof(WaveFunctionType));
    while(BlockOperatorPsiCC(H, v[vi], v[vip], Npsivi))
    {
      printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeCC for the function BlockOperatorPsiRR\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      if (v[vi][0].NBlok == 1) AllocFailed();
      ComplexBlockWaveFunctionTruncateBasis(v[vi], Npsivi, epsilon);
    }
    // a[k] = v[vi]^* . v[vip]^T  (v[vi]^* . H . v[vi]^T)
    Mab.V = &(M->a[k*Sqr(M->BlockSize)]);
    ComplexBlockWaveFunctionDotProduct(v[vi], Npsi, v[vip], Npsi, &Mab); //Npsivi = Npsi - changes at Lowdin Orthonormalization
    if(k!=1)
    {
      // b[k-1] = v[vim]^* . v[vip]^T  (v[vim]^* . H . v[vi]^T)
      Mab.V = &(M->b[(k-1)*Sqr(M->BlockSize)]);
      ComplexBlockWaveFunctionDotProduct(v[vim], Npsivim, v[vip], Npsi, &Mab);
    }
    NSingular=1;
    while(NSingular)
    {
      // orthogonalize v[vip] to v[vi] and v[vim]
      // we use matrix multiplication between v[vip] and v[vi] which is a numerical unstable form of Gram Schmidt
      // we heal this by repeating the orthogonalization twice for each matrix vector and once for vi and vim (factor of 3 reduction in speed but needed in tests)
      Mab.V=M->a;
      // w = v[vip]^* . v[vi]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vi], Npsivi, &Mab);
      // v[vip] -= w^* . v[vi]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vi], Npsivi, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeCC for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vi],  Npsivi, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,   epsilon);
      }
      // w = v[vip]^* . v[vi]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vi], Npsivi, &Mab);
      // v[vip] -= w^* . v[vi]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vi], Npsivi, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeCC for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vi],  Npsivi, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,   epsilon);
      }
      // w = v[vip]^* . v[vim]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vim], Npsivim, &Mab);
      // v[vip] -= w^* . v[vim]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vim], Npsivim, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeCC for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vim][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vim], Npsivim, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,    epsilon);
      }
      // w = v[vip]^* . v[vim]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vim], Npsivim, &Mab);
      // v[vip] -= w^* . v[vim]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vim], Npsivim, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeCC for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vim][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vim], Npsivim, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,    epsilon);
      }
      // w = v[vip]^* . v[vi]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vi], Npsivi, &Mab);
      // v[vip] -= w^* . v[vi]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vi], Npsivi, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeCC for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vi],  Npsivi, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,   epsilon);
      }
      // w = v[vip]^* . v[vim]^T
      CompactMatrixScale(&Mab, 0);
      ComplexBlockWaveFunctionDotProduct(v[vip], Npsi, v[vim], Npsivim, &Mab);
      // v[vip] -= w^* . v[vim]
      while(ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vim], Npsivim, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeCC for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vim][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        ComplexBlockWaveFunctionTruncateBasis(v[vim], Npsivim, epsilon);
        ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi,    epsilon);
      }
      // Lowdin orthonormalize v[vip] into v[3], free v[vip] and set v[vip] equal to v[3]
      NSingular=Npsi;
      ComplexWaveFunctionListSVDOrthogonalizeSameBasis(v[vip], v[3], &Npsi, epsilon, true);
      NSingular-=Npsi;
      if(NSingular!=0)
      {
        VerbosityPrintf(0x0001, "Reduce dimension from %u to %u in block %u\n",Npsi+NSingular,Npsi,k);
      }
      for(i=0;i<Npsi+NSingular;i++)
      {
        ComplexWaveFunctionFreeValues(&v[vip][i]);
      }
      for(i=0;i<Npsi;i++)
      {
        v[vip][i]=v[3][i];
      }
      // remove small values from v[vip]
      ComplexBlockWaveFunctionTruncateBasis(v[vip], Npsi, epsilon);
    }
    BlockSize[k]=Npsi;
    psiBasis[k]=v[vip];
    if(BlockSize[k]==0)
    {
      free(v[vip]);
      psiBasis[k]=NULL;
      for(i=0;i<Sqr(M->BlockSize);i++)
      {
        M->a[i].Re=0;
        M->a[i].Im=0;
      }
      free(v[3]);
      return 0;
    }
    // set pointers for next state
    vim=vi;
    vi=vip;
    vip=(vip+1)%3;
    Npsivim=Npsivi;
    Npsivi=Npsi;
  }
  for(i=0;i<Sqr(M->BlockSize);i++)
  {
    M->a[i].Re=0;
    M->a[i].Im=0;
  }
  free(v[3]);
  return 0;
}
uint32_t LanczosBlockTriDiagonalizeKrylovRR(const OperatorType H, WaveFunctionType * const psi, double epsilon, BlockTriDiagonalMatrixType * const M, WaveFunctionType ** psiBasis, uint32_t *BlockSize)
{
  uint32_t i,j,k,vi,vip,vim;
  CompactMatrixType Mab;
  WaveFunctionType *v[4];
  uint32_t NSingular,Npsi=M->BlockSize,Npsivi=M->BlockSize,Npsivim=M->BlockSize;
  VerbosityPrintf(0x0001,"Start of LanczosBlockTriDiagonalizeRR\n");
  fflush(stdout);
  for(i=0;i<M->BlockSize;i++)
  {
    M->bRe[i+i*M->BlockSize]=0;
  }
  v[0] = psi;
  v[3] = (WaveFunctionType*)malloc(Npsi*sizeof(WaveFunctionType));
  
  Mab.Ni=M->BlockSize;
  Mab.Nj=M->BlockSize;
  Mab.QComplex=0;
  v[1] = (WaveFunctionType*)malloc(Npsi*sizeof(WaveFunctionType));
  RealWaveFunctionListSVDOrthogonalize(v[0], v[1], &Npsi, epsilon, true); // v[1] = orthonormal[v0]
  BlockSize[0]=Npsi;
  psiBasis[0]=v[1];
  if(BlockSize[0]==0)
  {
    free(v[1]);
    psiBasis[0]=NULL;
    for(i=0;i<Sqr(M->BlockSize);i++)
    {
      M->aRe[i]=0;
      M->bRe[i]=0;
    }
    free(v[3]);
    return 0;
  }
  NSingular=M->BlockSize-Npsi;
  if(NSingular!=0)
  {
    VerbosityPrintf(0x0001,"Reduce dimension from %u to %u in block %u\n",M->BlockSize,M->BlockSize-NSingular,0);
  }
  // b[0] = v[0] . v[1]^+ or b[0]^T = v[1]^* . v[0]^T
  for(i=0;i<Npsi+NSingular;i++)
  {
    for(j=0;j<Npsi;j++)
    {
      M->bRe[i * M->BlockSize + j] = RealWaveFunctionDotProduct(&(v[1][j]), &(v[0][i]));
    }
  }
  vi=1;
  vim=0;
  vip=2;
  Npsivi=Npsi;
  for(k=1;k<M->N;k++)
  {
    // v[vip] = H v[vi]^T
    v[vip] = (WaveFunctionType*)malloc(Npsivi*sizeof(WaveFunctionType));
    while(BlockOperatorPsiRR(H, v[vi], v[vip], Npsivi))
    {
      printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRR for the function BlockOperatorPsiRR\n");
      printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
      fflush(stdout);
      epsilon *= 2.0;
      if (v[vi][0].NBlok == 1) AllocFailed();
      RealBlockWaveFunctionTruncateBasis(v[vi], Npsivi, epsilon);
    }
    // a[k] = v[vi]^* . v[vip]^T  (v[vi]^* . H . v[vi]^T)
    Mab.VRe = &(M->aRe[k*Sqr(M->BlockSize)]);
    RealBlockWaveFunctionDotProduct(v[vi], Npsi, v[vip], Npsi, &Mab); //Npsivi = Npsi - changes at Lowdin Orthonormalization
    if(k!=1)
    {
      // b[k-1] = v[vim]^* . v[vip]^T  (v[vim]^* . H . v[vi]^T)
      Mab.VRe = &(M->bRe[(k-1)*Sqr(M->BlockSize)]);
      RealBlockWaveFunctionDotProduct(v[vim], Npsivim, v[vip], Npsi, &Mab);
    }
    NSingular=1;
    while(NSingular)
    {
      // orthogonalize v[vip] to v[vi] and v[vim]
      // we use matrix multiplication between v[vip] and v[vi] which is a numerical unstable form of Gram Schmidt
      // we heal this by repeating the orthogonalization twice for each matrix vector and once for vi and vim (factor of 3 reduction in speed but needed in tests)
      Mab.VRe=M->aRe;
      // w = v[vip]^* . v[vi]^T
      CompactMatrixScale(&Mab, 0);
      RealBlockWaveFunctionDotProduct(v[vip], Npsi, v[vi], Npsivi, &Mab);
      // v[vip] -= w^* . v[vi]
      while(RealBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vi], Npsivi, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRR for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        RealBlockWaveFunctionTruncateBasis(v[vi],  Npsivi, epsilon);
        RealBlockWaveFunctionTruncateBasis(v[vip], Npsi,   epsilon);
      }
      // w = v[vip]^* . v[vi]^T
      CompactMatrixScale(&Mab, 0);
      RealBlockWaveFunctionDotProduct(v[vip], Npsi, v[vi], Npsivi, &Mab);
      // v[vip] -= w^* . v[vi]
      while(RealBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vi], Npsivi, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRR for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        RealBlockWaveFunctionTruncateBasis(v[vi],  Npsivi, epsilon);
        RealBlockWaveFunctionTruncateBasis(v[vip], Npsi,   epsilon);
      }
      // w = v[vip]^* . v[vim]^T
      CompactMatrixScale(&Mab, 0);
      RealBlockWaveFunctionDotProduct(v[vip], Npsi, v[vim], Npsivim, &Mab);
      // v[vip] -= w^* . v[vim]
      while(RealBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vim], Npsivim, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRR for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vim][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        RealBlockWaveFunctionTruncateBasis(v[vim], Npsivim, epsilon);
        RealBlockWaveFunctionTruncateBasis(v[vip], Npsi,    epsilon);
      }
      // w = v[vip]^* . v[vim]^T
      CompactMatrixScale(&Mab, 0);
      RealBlockWaveFunctionDotProduct(v[vip], Npsi, v[vim], Npsivim, &Mab);
      // v[vip] -= w^* . v[vim]
      while(RealBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vim], Npsivim, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRR for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vim][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        RealBlockWaveFunctionTruncateBasis(v[vim], Npsivim, epsilon);
        RealBlockWaveFunctionTruncateBasis(v[vip], Npsi,    epsilon);
      }
      // w = v[vip]^* . v[vi]^T
      CompactMatrixScale(&Mab, 0);
      RealBlockWaveFunctionDotProduct(v[vip], Npsi, v[vi], Npsivi, &Mab);
      // v[vip] -= w^* . v[vi]
      while(RealBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vi], Npsivi, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRR for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vi][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        RealBlockWaveFunctionTruncateBasis(v[vi],  Npsivi, epsilon);
        RealBlockWaveFunctionTruncateBasis(v[vip], Npsi,   epsilon);
      }
      // w = v[vip]^* . v[vim]^T
      CompactMatrixScale(&Mab, 0);
      RealBlockWaveFunctionDotProduct(v[vip], Npsi, v[vim], Npsivim, &Mab);
      // v[vip] -= w^* . v[vim]
      while(RealBlockWaveFunctionSubtractMultipliedUnionBasis(v[vip], Npsi, v[vim], Npsivim, Mab))
      {
        printf("Allocation of new memory failed in LanczosBlockTriDiagonalizeRR for the function RealBlockWaveFunctionSubtractMultipliedUnionBasis\n");
        printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n", epsilon, epsilon * 2.0);
        fflush(stdout);
        epsilon *= 2.0;
        if ((v[vim][0].NBlok == 1)&&(v[vip][0].NBlok == 1)) AllocFailed();
        RealBlockWaveFunctionTruncateBasis(v[vim], Npsivim, epsilon);
        RealBlockWaveFunctionTruncateBasis(v[vip], Npsi,    epsilon);
      }
      // Lowdin orthonormalize v[vip] into v[3], free v[vip] and set v[vip] equal to v[3]
      NSingular=Npsi;
      RealWaveFunctionListSVDOrthogonalizeSameBasis(v[vip], v[3], &Npsi, epsilon, true);
      NSingular-=Npsi;
      if(NSingular!=0)
      {
        VerbosityPrintf(0x0001, "Reduce dimension from %u to %u in block %u\n",Npsi+NSingular,Npsi,k);
      }
      for(i=0;i<Npsi+NSingular;i++)
      {
        WaveFunctionFreeValues(&v[vip][i]);
      }
      for(i=0;i<Npsi;i++)
      {
        v[vip][i]=v[3][i];
      }
      // remove small values from v[vip]
      BlockWaveFunctionTruncateBasis(v[vip], Npsi, epsilon);
    }
    BlockSize[k]=Npsi;
    psiBasis[k]=v[vip];
    if(BlockSize[k]==0)
    {
      free(v[vip]);
      psiBasis[k]=NULL;
      for(i=0;i<Sqr(M->BlockSize);i++)
      {
        M->aRe[i]=0;
      }
      free(v[3]);
      return 0;
    }
    // set pointers for next state
    vim=vi;
    vi=vip;
    vip=(vip+1)%3;
    Npsivim=Npsivi;
    Npsivi=Npsi;
  }
  for(i=0;i<Sqr(M->BlockSize);i++)
  {
    M->aRe[i]=0;
  }
  free(v[3]);
  return 0;
}
uint32_t LanczosBlockTriDiagonalizeKrylovMR(const OperatorType H, WaveFunctionType * const psi, double epsilon, BlockTriDiagonalMatrixType * const M, WaveFunctionType ** psiBasis, uint32_t *BlockSize)
{
  // make psi complex, call MC
  uint32_t i;
  VerbosityPrintf(0x0001,"Start of LanczosBlockTriDiagonalizeKrylovMR\n");
  fflush(stdout);
  for(i=0;i<M->BlockSize;i++)
  {
    if (MakeComplexWaveFunctionFromReal(&psi[i], psi[i]))
    {
      printf("LanczosBlockTriDiagonalizeKrylovMR called with a real wave-function and a complex Hamiltonian.\nTurned the wave-function complex in step 1.\nBut failed at function %i.\n",i);
      fflush(stdout);
      for(;i>0;i--)
      {
        MakeRealWaveFunctionFromComplex(psi[i-1], &psi[i-1]);
      }
      return 1;
    }
  }
  return LanczosBlockTriDiagonalizeKrylovMC(H,psi,epsilon,M,psiBasis,BlockSize);
}
uint32_t LanczosBlockTriDiagonalizeKrylovCR(const OperatorType H, WaveFunctionType * const psi, double epsilon, BlockTriDiagonalMatrixType * const M, WaveFunctionType ** psiBasis, uint32_t *BlockSize)
{
  // make psi complex, and call CC
  uint32_t i;
  VerbosityPrintf(0x0001,"Start of LanczosBlockTriDiagonalizeKrylovCR\n");
  fflush(stdout);
  for(i=0;i<M->BlockSize;i++)
  {
    if (MakeComplexWaveFunctionFromReal(&psi[i], psi[i]))
    {
      printf("LanczosBlockTriDiagonalizeKrylovCR called with a real wave-function and a complex Hamiltonian.\nTurned the wave-function complex in step 1.\nBut failed at function %i.\n",i);
      fflush(stdout);
      for(;i>0;i--)
      {
        MakeRealWaveFunctionFromComplex(psi[i-1], &psi[i-1]);
      }
      return 1;
    }
  }
  return LanczosBlockTriDiagonalizeKrylovCC(H,psi,epsilon,M,psiBasis,BlockSize);
}


uint32_t LanczosBlockTriDiagonalizeRestricted(const OperatorType H, WaveFunctionType * const psi, DeterminantRestrictionType Restrictions, double epsilon, BlockTriDiagonalMatrixType * const M, const double Gamma)
{
  uint32_t i,psiComplex=0;
  VerbosityPrintf(0x0001,"Start of LanczosBlockTriDiagonalizeRestricted\n");
  // all psi complex ?
  for(i=0;i<M->BlockSize;i++)
  {
    if(psi[i].QComplex)psiComplex++;
  }
  if(!((psiComplex==M->BlockSize)||(psiComplex==0)))
  {
    for(i=0;i<M->BlockSize;i++)
    {
      if (MakeComplexWaveFunctionFromReal(&psi[i], psi[i]))
      {
        printf("LanczosBlockTriDiagonalizeRestricted called with a set of wave-functions that are mixed real and complex.\nTurned the wave-function complex in step 1.\nBut failed at function %i.\n",i);
        fflush(stdout);
        for(;i>0;i--)
        {
          MakeRealWaveFunctionFromComplex(psi[i-1], &psi[i-1]);
        }
        return 1;
      }
    }
  }
  if (psiComplex)
  {
    if (H.QComplex == 0)
    {
      return LanczosBlockTriDiagonalizeRestrictedRC(H, psi, Restrictions, epsilon, M, Gamma);
    }
    else if (H.QComplex == 1)
    {
      return LanczosBlockTriDiagonalizeRestrictedCC(H, psi, Restrictions, epsilon, M, Gamma);
    }
    else
    {
      return LanczosBlockTriDiagonalizeRestrictedMC(H, psi, Restrictions, epsilon, M, Gamma);
    }
  }
  else
  {
    if (H.QComplex == 0)
    {
      return LanczosBlockTriDiagonalizeRestrictedRR(H, psi, Restrictions, epsilon, M, Gamma);
    }
    else if (H.QComplex == 1)
    {
      return LanczosBlockTriDiagonalizeRestrictedCR(H, psi, Restrictions, epsilon, M, Gamma);
    }
    else
    {
      return LanczosBlockTriDiagonalizeRestrictedMR(H, psi, Restrictions, epsilon, M, Gamma);
    }
  }
}

uint32_t LanczosBlockTriDiagonalize(const OperatorType H, WaveFunctionType * const psi, double epsilon, BlockTriDiagonalMatrixType * const M, const double Gamma)
{
  uint32_t i,psiComplex=0;
  VerbosityPrintf(0x0001,"Start of LanczosBlockTriDiagonalize\n");
  // restrictions ?
  if(H.Res.N!=0)
  {
    BlockWaveFunctionRestrict(psi, M->BlockSize, H.Res);
    return LanczosBlockTriDiagonalizeRestricted(H, psi, H.Res, epsilon, M, Gamma);
  }
  // all psi complex ?
  for(i=0;i<M->BlockSize;i++)
  {
    if(psi[i].QComplex)psiComplex++;
  }
  if(!((psiComplex==M->BlockSize)||(psiComplex==0)))
  {
    for(i=0;i<M->BlockSize;i++)
    {
      if (MakeComplexWaveFunctionFromReal(&psi[i], psi[i]))
      {
        printf("LanczosBlockTriDiagonalize called with a set of wave-functions that are mixed real and complex.\nTurned the wave-function complex in step 1.\nBut failed at function %i.\n",i);
        fflush(stdout);
        for(;i>0;i--)
        {
          MakeRealWaveFunctionFromComplex(psi[i-1], &psi[i-1]);
        }
        return 1;
      }
    }
  }
  if (psiComplex)
  {
    if (H.QComplex == 0)
    {
      return LanczosBlockTriDiagonalizeRC(H, psi, epsilon, M, Gamma);
    }
    else if (H.QComplex == 1)
    {
      return LanczosBlockTriDiagonalizeCC(H, psi, epsilon, M, Gamma);
    }
    else
    {
      return LanczosBlockTriDiagonalizeMC(H, psi, epsilon, M, Gamma);
    }
  }
  else
  {
    if (H.QComplex == 0)
    {
      return LanczosBlockTriDiagonalizeRR(H, psi, epsilon, M, Gamma);
    }
    else if (H.QComplex == 1)
    {
      return LanczosBlockTriDiagonalizeCR(H, psi, epsilon, M, Gamma);
    }
    else
    {
      return LanczosBlockTriDiagonalizeMR(H, psi, epsilon, M, Gamma);
    }
  }
}

uint32_t LanczosBlockTriDiagonalizeKrylovRestricted(const OperatorType H, WaveFunctionType * const psi, DeterminantRestrictionType Restrictions, double epsilon, BlockTriDiagonalMatrixType * const M, WaveFunctionType ** psiBasis, uint32_t *BlockSize)
{
  uint32_t i,psiComplex=0;
  VerbosityPrintf(0x0001,"Start of LanczosBlockTriDiagonalize Krylov Restricted \n");
  // all psi complex ?
  for(i=0;i<M->BlockSize;i++)
  {
    if(psi[i].QComplex)psiComplex++;
  }
  if(!((psiComplex==M->BlockSize)||(psiComplex==0)))
  {
    for(i=0;i<M->BlockSize;i++)
    {
      if (MakeComplexWaveFunctionFromReal(&psi[i], psi[i]))
      {
        printf("LanczosBlockTriDiagonalizeKrylovRestricted called with a set of wave-functions that are mixed real and complex.\nTurned the wave-function complex in step 1.\nBut failed at function %i.\n",i);
        fflush(stdout);
        for(;i>0;i--)
        {
          MakeRealWaveFunctionFromComplex(psi[i-1], &psi[i-1]);
        }
        return 1;
      }
    }
  }
  if (psiComplex)
  {
    if (H.QComplex == 0)
    {
      return LanczosBlockTriDiagonalizeKrylovRestrictedRC(H, psi, Restrictions, epsilon, M, psiBasis, BlockSize);
    }
    else if (H.QComplex == 1)
    {
      return LanczosBlockTriDiagonalizeKrylovRestrictedCC(H, psi, Restrictions, epsilon, M, psiBasis, BlockSize);
    }
    else
    {
      return LanczosBlockTriDiagonalizeKrylovRestrictedMC(H, psi, Restrictions, epsilon, M, psiBasis, BlockSize);
    }
  }
  else
  {
    if (H.QComplex == 0)
    {
      return LanczosBlockTriDiagonalizeKrylovRestrictedRR(H, psi, Restrictions, epsilon, M, psiBasis, BlockSize);
    }
    else if (H.QComplex == 1)
    {
      return LanczosBlockTriDiagonalizeKrylovRestrictedCR(H, psi, Restrictions, epsilon, M, psiBasis, BlockSize);
    }
    else
    {
      return LanczosBlockTriDiagonalizeKrylovRestrictedMR(H, psi, Restrictions, epsilon, M, psiBasis, BlockSize);
    }
  }
}
 
/*  The Krylov basis is returned as a list of M->N pointers each poiter allocates a set of wave-functions that have the same basis
 *  The number of functions to allocate is M->BlockSize-Nsingular[i] = BlockSize[i]
 *  The wave functions psi are freed (or used in psiBasis) on succesfull exit
 */
uint32_t LanczosBlockTriDiagonalizeKrylov(const OperatorType H, WaveFunctionType * const psi, double epsilon, BlockTriDiagonalMatrixType * const M, WaveFunctionType ** psiBasis, uint32_t *BlockSize)
{
  uint32_t i,psiComplex=0;
  VerbosityPrintf(0x0001,"Start of LanczosBlockTriDiagonalize Krylov\n");
  // restrictions ?
  if(H.Res.N!=0)
  {
    BlockWaveFunctionRestrict(psi, M->BlockSize, H.Res);
    return LanczosBlockTriDiagonalizeKrylovRestricted(H, psi, H.Res, epsilon, M, psiBasis, BlockSize);
  }
  // all psi complex ?
  for(i=0;i<M->BlockSize;i++)
  {
    if(psi[i].QComplex)psiComplex++;
  }
  if(!((psiComplex==M->BlockSize)||(psiComplex==0)))
  {
    for(i=0;i<M->BlockSize;i++)
    {
      if (MakeComplexWaveFunctionFromReal(&psi[i], psi[i]))
      {
        printf("LanczosBlockTriDiagonalizeKrylov called with a set of wave-functions that are mixed real and complex.\nTurned the wave-function complex in step 1.\nBut failed at function %i.\n",i);
        fflush(stdout);
        for(;i>0;i--)
        {
          MakeRealWaveFunctionFromComplex(psi[i-1], &psi[i-1]);
        }
        return 1;
      }
    }
  }
  if (psiComplex)
  {
    if (H.QComplex == 0)
    {
      return LanczosBlockTriDiagonalizeKrylovRC(H, psi, epsilon, M, psiBasis, BlockSize);
    }
    else if (H.QComplex == 1)
    {
      return LanczosBlockTriDiagonalizeKrylovCC(H, psi, epsilon, M, psiBasis, BlockSize);
    }
    else
    {
      return LanczosBlockTriDiagonalizeKrylovMC(H, psi, epsilon, M, psiBasis, BlockSize);
    }
  }
  else
  {
    if (H.QComplex == 0)
    {
      return LanczosBlockTriDiagonalizeKrylovRR(H, psi, epsilon, M, psiBasis, BlockSize);
    }
    else if (H.QComplex == 1)
    {
      return LanczosBlockTriDiagonalizeKrylovCR(H, psi, epsilon, M, psiBasis, BlockSize);
    }
    else
    {
      return LanczosBlockTriDiagonalizeKrylovMR(H, psi, epsilon, M, psiBasis, BlockSize);
    }
  }
}

/*******************************************************************
 **************************** Ground-state *************************
 *******************************************************************/

uint32_t TimeToDeflateInBlockLanczosKrylov(double Norm, double Zero, uint32_t i, uint32_t NpsiDeflated, CompactMatrixType *HKrylov)
{
  uint32_t j;
  double tot=sqrt(DBL_MIN);
  if(HKrylov->QComplex==0)
  {
    for(j=0;j<i;j++)
    {
      tot+=Sqr(HKrylov->VRe[j+j*HKrylov->Nj]);
    }
  }
  else
  {
    for(j=0;j<i;j++)
    {
      tot+=Sqr(HKrylov->V[j+j*HKrylov->Nj].Re);
    }
  }
  
  if(Sqr(Norm)<tot*Abs(Zero)) return 1;
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
 *
 *  for restarts the energy is shifted to the previous largest energy in the krylov basis
 */
uint32_t BlockLanczosGroundStateConserveBasisKrylovReal(OperatorType * H, WaveFunctionType * const v, uint32_t const Npsi, ConvergenceSettingsType * Settings)
{
  VerbosityPrintf(0x0100,"In BlockLanczosGroundStateConserveBasisKrylovReal\n");
  //CompactMatrixType tempS; // needed for testing of norm during convergence loops
  int    info;
  const char *jobz = "V";
  const char *uplo = "U";
  uint32_t    lwork = 1 + 6 * (Npsi+Settings->NKrylovMax) + 2 * Sqr((Npsi+Settings->NKrylovMax));
  uint32_t    liwork = 3 + 5 * (Npsi+Settings->NKrylovMax);
  int    *iwork;
  double *work;
  double *val;
  double Norm, DotRe, HE0;
  uint32_t NKrylovNow=Npsi+Settings->NKrylovStart;
  uint32_t i,j, NpsiDeflated, NorthoExtra, Deflated[Npsi],irestart=0;
  CompactMatrixType HKrylov;
  uint32_t QConverged=0;
  double **temp;
  // get value of H->l[0] so we can shift and restore the operator
  if(H->l[0].NAloc==0) OperatorAddLader(H, NULL, 0, 0, 0); // allocate the operator with length 0
  if(H->l[0].QComplex==0)
  {
    HE0 = H->l[0].VRe[0];
  }
  else
  {
    HE0 = H->l[0].V[0].Re;
  }
  // init several matrices
  val= (double*)malloc((Npsi+Settings->NKrylovMax)*sizeof(double));
  iwork = (int*)malloc(liwork*sizeof(int));
  if(iwork==NULL)
  {
    printf("malloc for iwork failed in BlockLanczosGroundStateConserveBasisKrylovReal\n");
    free(val);
    return 1;
  }
  work = (double*)malloc(lwork*sizeof(double));
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
      while(TimeToDeflateInBlockLanczosKrylov(Norm, Settings->zero, i, NpsiDeflated, &HKrylov))
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
      OperatorPsiConserveBasisPreAllocatedRR(H, &v[i], &v[i+NpsiDeflated]);
      HKrylov.VRe[i+i*HKrylov.Nj]=RealWaveFunctionDotProductSameBasis(v[i], v[i+NpsiDeflated]);
      // The Hamiltonian should be negative definite. If the onsite term is positive we should shift H
      // we can do this by changing H and subtracting the energy from all previous diagonal elements
      // but that should be implemented stop the loop as an alternative
      //if(HKrylov.VRe[i+i*HKrylov.Nj]>0)
      //{
      //  NKrylovNow = Min(i+Npsi,NKrylovNow);
      //}
      // orthogonalize v[i+NpsiDeflated] to v[j] (the new candiate agains all found states)
      // for j from i-NpsiDeflated to i
      for(j=(i-Min(i, NpsiDeflated));j<=i;j++)
      {
        HKrylov.VRe[j+i*HKrylov.Nj]=HKrylov.VRe[i+j*HKrylov.Nj];
        RealWaveFunctionAddMultipliedSameBasis(v[i+NpsiDeflated], v[j], -HKrylov.VRe[j+i*HKrylov.Nj]);
      }
      /* // print aditional information on convergence and norms of v
      printf("i+NpsiDeflated = %i + %i\n",i,NpsiDeflated);
      BlockWaveFunctionCalculatOverlapMatrixSerial(v, NpsiDeflated+Npsi, &tempS);
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
    if(Npsi<(NKrylovNow-1))
    {
      VerbosityPrintf(0x0100,"    Eigenvalue ...\n");
      VerbosityPrintf(0x0100,"    Eigenvalue %3i = %22.15E\n",NKrylovNow-1,val[NKrylovNow-1]);
    }
    fflush(stdout);
    // shift the Hamiltonian
    if(strcmp("Automatic",Settings->ShiftScheme)==0)
    {
      OperatorAddLader(H, NULL, 0, -(val[NKrylovNow-1]+val[Npsi-1])/2, 0);
    }
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
  // undo shift
  if(H->l[0].QComplex==0)
  {
    H->l[0].VRe[0] = HE0;
  }
  else
  {
    H->l[0].V[0].Re = HE0;
  }
  return 0;
}

uint32_t BlockLanczosGroundStateConserveBasisKrylovComplex(OperatorType * H, WaveFunctionType * const v, uint32_t const Npsi, ConvergenceSettingsType * Settings)
{
  VerbosityPrintf(0x0100,"In BlockLanczosGroundStateConserveBasisKrylovComplex\n");
  //CompactMatrixType tempS; // needed for testing of norm during convergence loops
  int    info;
  const char *jobz = "V";
  const char *uplo = "U";
  uint32_t    lwork  = 1 + 6 * (Npsi+Settings->NKrylovMax) + 2 * Sqr((Npsi+Settings->NKrylovMax));
  uint32_t    lcwork = 2 * (Npsi+Settings->NKrylovMax) + Sqr((Npsi+Settings->NKrylovMax));
  uint32_t    liwork = 3 + 5 * (Npsi+Settings->NKrylovMax);
  int    *iwork;
  double *work;
  Complex16 *cwork;
  double *val;
  double Norm, DotRe, DotIm, HE0;
  uint32_t NKrylovNow=Npsi+Settings->NKrylovStart;
  uint32_t i,j, NpsiDeflated, NorthoExtra, Deflated[Npsi],irestart=0;
  CompactMatrixType HKrylov;
  uint32_t QConverged=0;
  double **tempRe, **tempIm;
  // get value of H->l[0] so we can shift and restore the operator
  if(H->l[0].NAloc==0) OperatorAddLader(H, NULL, 0, 0, 0); // allocate the operator with length 0
  if(H->l[0].QComplex==0)
  {
    HE0 = H->l[0].VRe[0];
  }
  else
  {
    HE0 = H->l[0].V[0].Re;
  }
  // init several matrices
  val= (double*)malloc((Npsi+Settings->NKrylovMax)*sizeof(double));
  iwork = (int*)malloc(liwork*sizeof(int));
  if(iwork==NULL)
  {
    printf("malloc for iwork failed in BlockLanczosGroundStateConserveBasisKrylovComplex\n");
    free(val);
    return 1;
  }
  work = (double*)malloc(lwork*sizeof(double));
  if(work==NULL)
  {
    printf("malloc for work failed in BlockLanczosGroundStateConserveBasisKrylovComplex\n");
    free(val);
    free(iwork);
    return 1;
  }
  cwork = (Complex16*)malloc(lcwork*sizeof(Complex16));
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
      while(TimeToDeflateInBlockLanczosKrylov(Norm, Settings->zero, i, NpsiDeflated, &HKrylov))
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
      if(H->QComplex==1)
      {
        OperatorPsiConserveBasisPreAllocatedCC(H, &v[i], &v[i+NpsiDeflated]);
      }
      else
      {
        OperatorPsiConserveBasisPreAllocatedMC(H, &v[i], &v[i+NpsiDeflated]);
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
       BlockWaveFunctionCalculatOverlapMatrixSerial(v, NpsiDeflated+i, &tempS);
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
    if(Npsi<(NKrylovNow-1))
    {
      VerbosityPrintf(0x0100,"    Eigenvalue ...\n");
      VerbosityPrintf(0x0100,"    Eigenvalue %3i = %22.15E\n",NKrylovNow-1,val[NKrylovNow-1]);
    }
    fflush(stdout);
    // shift the Hamiltonian
    if(strcmp("Automatic",Settings->ShiftScheme)==0)
    {
      OperatorAddLader(H, NULL, 0, -(val[NKrylovNow-1]+val[Npsi-1])/2, 0);
    }
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
  // undo shift
  if(H->l[0].QComplex==0)
  {
    H->l[0].VRe[0] = HE0;
  }
  else
  {
    H->l[0].V[0].Re = HE0;
  }
  return 0;
}

uint32_t BlockLanczosGroundStateConserveBasisKrylov(OperatorType *H, WaveFunctionType * const v, uint32_t const Npsi, ConvergenceSettingsType * Settings)
{
  if(H->QComplex)
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
uint32_t BlockLanczosGroundStateConserveBasisKrylovRecalculate(OperatorType * H, WaveFunctionType * const v, uint32_t const Npsi, uint32_t const NpsiAlloc, ConvergenceSettingsType * Settings)
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
 *     DenseBorder:  (the number of basis states below which a dense calculation is done)
 *     NKrlylovStart: (the number of KrylovVectors to start with)
 *
 *   if psi[0].N <= DenseBorder do a dense calculation
 *
 */
uint32_t BlockLanczosGroundStateConserveBasis(OperatorType * H, WaveFunctionType * const psi, uint32_t const Npsi, ConvergenceSettingsType * Settings)
{
  uint32_t i,MaxAlloc=2*Npsi+Settings->NKrylovMax, StartAlloc=2*Npsi+Settings->NKrylovStart;
  WaveFunctionType *v;
  VerbosityPrintf(0x10000000,"Start of BlockLanczosGroundStateConserveBasis\n");
  // for small basises do a dense calculation
  if(psi[0].N<=Settings->DenseBorder) return BlockLanczosGroundStateConserveBasisDense(*H,psi, Npsi);
  // create the vectors v in memory.
  //   if all 2*Npsi+NKrylovMax are created return BlockLanczosGroundStateConserveBasisKrylov
  //   if only 2*Npsi+1 states can be created return BlockLanczosGroundStateConserveBasisKrylovRecalculate
  //   if less then 2*Npsi+1 states can be created return 1
  v= (WaveFunctionType*)malloc(MaxAlloc*sizeof(WaveFunctionType));
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
 *    H:             Operator assumed to be Hermitian (not tested and function fails if it is not)
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
uint32_t BlockGroundState(OperatorType *H, uint32_t const Npsi, ConvergenceSettingsType * Settings, DeterminantRestrictionType const StartingPoint, DeterminantRestrictionType const Restrictions, WaveFunctionType * const psi, const uint32_t NBitsKey)
{
  uint32_t i,Qconverged=0, AtLeastAnOtherLoop=0,iloop=0,NBitsKeyNew;
  WaveFunctionType Opsi[Npsi];
  double Energy[2][Npsi],EnergySqr[2][Npsi],temp,shiftE,shiftVar,HE0,totESqr;
  // in order to allow shifts we need H.l[0].VRe[0]
  if(H->l[0].NAloc==0) OperatorAddLader(H, NULL, 0, 0, 0); // allocate the operator with length 0
  if(H->l[0].QComplex==0)
  {
    HE0 = H->l[0].VRe[0];
  }
  else
  {
    HE0 = H->l[0].V[0].Re;
  }
  // if StartingPoint.N != 0
  // create Npsi independent functions fullfiling the starting point conditions
  // all functions must have the same basis
  // if H is complex or mixed psi will be complex, if H is real psi is real
  VerbosityPrintf(0x10000000,"Start of BlockGroundState\n");
  if((StartingPoint.N!=0)&&(!(NBitsKey==0)))
  {
    if(CreateStartingFunctions(*H, Npsi, StartingPoint, psi, NBitsKey))
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
  // Change the size of the Hash lookuptable
  NBitsKeyNew = 8;
  while(psi[0].N>(1<<(NBitsKeyNew-1)))
  {
    NBitsKeyNew++;
  }
  NBitsKeyNew++;
  if(psi[0].NBitsKey!=NBitsKeyNew) BlockWaveFunctionResetKeylist(psi, NBitsKeyNew, Npsi);
  // create some feeling for the energy of the state
  WaveFunctionCopy(psi[0], &(Opsi[0]));
  WaveFunctionRandomize(&(Opsi[0]));
  WaveFunctionNormalize(&(Opsi[0]));
  PsiHermitianOperatorPsi(*H, Opsi[0], &totESqr);
  totESqr*=totESqr;
  WaveFunctionFree(&(Opsi[0]));
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
    while( (Restrictions.N) ? BlockOperatorPsiSerialRestricted(*H, psi, Opsi, Npsi, Restrictions) : BlockOperatorPsiSerial(*H, psi, Opsi, Npsi))
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
    // Change the size of the Hash lookuptable
    NBitsKeyNew = 8;
    while(Opsi[0].N>(1<<(NBitsKeyNew-1)))
    {
      NBitsKeyNew++;
    }
    if(Opsi[0].NBitsKey!=NBitsKeyNew) BlockWaveFunctionResetKeylist(Opsi, NBitsKeyNew, Npsi);
    // calculate Energy and variance of Energy - these values are needed to test for convergence
    Qconverged=1;
    VerbosityPrintf(0x0100,"Number of determinants in basis %i (and basis of H psi %i) with size of lookuptable %i\n",psi[0].N,Opsi[0].N,1<<Opsi[0].NBitsKey);
    VerbosityPrintf(0x0100,"State # |   Energy               |  Variance of Energy    |  Energy change\n");
    for(i=0;i<Npsi;i++)
    {
      Energy[1][i]=Energy[0][i];
      EnergySqr[1][i]=EnergySqr[0][i];
      WaveFunctionDotProduct(psi[i], Opsi[i], &Energy[0][i], &temp);
      WaveFunctionDotProductSameBasis(Opsi[i], Opsi[i], &EnergySqr[0][i], &temp);
      totESqr+=EnergySqr[0][i];
      VerbosityPrintf(0x0100,"%7i | %22.15E | %22.15E | %22.15E\n",i,Energy[0][i],EnergySqr[0][i]-Sqr(Energy[0][i]),Energy[0][i]-Energy[1][i]);
    }
    // set Qconverged
    for(i=0;i<Npsi;i++)
    {
      if(EnergySqr[0][i]-Sqr(Energy[0][i])>(DBL_MIN+totESqr*Settings->zero/Npsi)) Qconverged=0;
    }
    VerbosityPrintf(0x0100,"\n");
    VerbosityPrintf(0x10000000,"Converged %i\n",Qconverged);
    // converge to ground-state in current basis or stop if converged
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
      if(!(strcmp("None", Settings->ShiftScheme)==0))
      {
        // we shift the Hamiltonian by the largest energy of psi[i]
        // and by ten times the largest variance of psi[i]
        shiftE=-1E+99;
        shiftVar=-1E+99;
        for(i=0;i<Npsi;i++)
        {
          shiftE   = Max(shiftE,Energy[0][i]);
          shiftVar = Max(shiftVar, EnergySqr[0][i] - Sqr(Energy[0][i]));
        }
        shiftVar = sqrt(Abs(shiftVar)+1000*DBL_MIN);
        // shift H
        OperatorAddLader(H, NULL, 0, -shiftE-10*shiftVar, 0);
        // appli this shift in retrospect for the calculation of Opsi = H psi
        for(i=0;i<Npsi;i++)
        {
          WaveFunctionAddMultipliedConserveBasis(&(Opsi[i]), psi[i], -shiftE-10*shiftVar);
        }
      }
      // print info
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
      // undo shift
      if(H->l[0].QComplex==0)
      {
        H->l[0].VRe[0] = HE0;
      }
      else
      {
        H->l[0].V[0].Re = HE0;
      }
      // truncate basis
      BlockWaveFunctionTruncateBasis(psi,Npsi,Settings->epsilon);
    }
  }
  return 0;
}

