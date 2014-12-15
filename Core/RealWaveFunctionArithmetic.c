

/*  Remove Wave function from memory
*   On entry psi contains a wave function with psi.NBlok alocations
*   On excit all allocate memory is freed
*/
void RealWaveFunctionFree(WaveFunctionType * const psi)
{
  uint32_t i;
  for (i = 0; i < psi->NBlok; i++)
  {
    free(psi->ValueRe[i]);
    free(psi->Next[i]);
    free(psi->Determinant[i]);
  }
  free(psi->ValueRe);
  free(psi->Next);
  free(psi->Determinant);
  free(psi->KeyList);
}

void RealWaveFunctionReset(WaveFunctionType * const psi)
{
  uint32_t i;
  for (i = 1; i < psi->NBlok; i++)
  {
    free(psi->ValueRe[i]);
    free(psi->Next[i]);
    free(psi->Determinant[i]);
  }
  psi->N = 0;
  psi->NBlok=1;
  memset(psi->KeyList, 0, ( 1 << (psi->NBitsKey)) * sizeof(uint32_t));
  psi->ValueRe[0][0] = 0.0;
}


/*  Remove The values of a wave function from memory, but keep the determiants and the pointers to next determinants with the same Hash value
*   Usefull for more then one wave function defined on the same basis where psi.Next and psi.Determinant point to the same location as used by an other wave function
*   On entry psi contains a wave function with psi.NBlok alocations
*   On excit all allocate memory for psi.ValueRe is freed
*/
void RealWaveFunctionFreeValues(WaveFunctionType * const psi)
{
  uint32_t i;
  for (i = 0; i < psi->NBlok; i++)
  {
    free(psi->ValueRe[i]);
  }
  free(psi->ValueRe);
}

/*
*   add ellement with single slater determinant det and pre-factor a to wave function psi
*   returns 0 if succeded
*   returns 1 if no memory available in order to add det to psi
*   If succeded psi -> psi+det
*   If not succeded psi -> psi
*/
uint32_t RealWaveFunctionAddElement(WaveFunctionType * const psi, const uint8_t * const det, const double a)
{
  uint32_t i, NextIndex, ThereIsANext, PsiAIsPsiB;
  void  *temp;
  /* test if the key that belongs to psi already excists */
  if (psi->KeyList[hash(det, psi->NWord, psi->NBitsKey)] == 0)
  {                             /* add a new element to psi with unique new key */
    psi->N++;
    /* test if one has to add a new block to psi in order to store the new element of psi */
    if (psi->N >= BlokSize * psi->NBlok)
    {                           /* add a new block to psi */
      psi->NBlok++;
      temp = psi->ValueRe;
      psi->ValueRe = realloc(psi->ValueRe, psi->NBlok * sizeof(double *));
      if (psi->ValueRe == NULL)
      {                         /* Allocation failed */
        printf("alloc failed RealWaveFunctionAddElement 01\n");
        fflush(stdout);
        /* restore psi and return 0 */
        psi->ValueRe = temp;
        psi->NBlok--;
        psi->N--;
        return (1);
      }
      temp = psi->Next;
      psi->Next = realloc(psi->Next, psi->NBlok * sizeof(uint32_t *));
      if (psi->Next == NULL)
      {                         /* Allocation failed */
        printf("alloc failed RealWaveFunctionAddElement 02\n");
        fflush(stdout);
        /* restore psi and return 0 */
        psi->Next = temp;
        psi->NBlok--;
        psi->N--;
        psi->ValueRe = realloc(psi->ValueRe, psi->NBlok * sizeof(double *));
        return (1);
      }
      temp = psi->Determinant;
      psi->Determinant = realloc(psi->Determinant, psi->NBlok * sizeof(uint8_t *));
      if (psi->Determinant == NULL)
      {                         /* Allocation failed */
        printf("alloc failed RealWaveFunctionAddElement 03\n");
        fflush(stdout);
        /* restore psi and return 0 */
        psi->Determinant = temp;
        psi->NBlok--;
        psi->N--;
        psi->ValueRe = realloc(psi->ValueRe, psi->NBlok * sizeof(double *));
        psi->Next = realloc(psi->Next, psi->NBlok * sizeof(uint32_t *));
        return (1);
      }
      psi->ValueRe[psi->NBlok - 1] = malloc(BlokSize * sizeof(double));
      if (psi->ValueRe[psi->NBlok - 1] == NULL)
      {                         /* Allocation failed */
        printf("alloc failed RealWaveFunctionAddElement 04\n");
        fflush(stdout);
        /* restore psi and return 0 */
        psi->NBlok--;
        psi->N--;
        psi->ValueRe = realloc(psi->ValueRe, psi->NBlok * sizeof(double *));
        psi->Next = realloc(psi->Next, psi->NBlok * sizeof(uint32_t *));
        psi->Determinant = realloc(psi->Determinant, psi->NBlok * sizeof(uint8_t *));
        return (1);
      }
      psi->Next[psi->NBlok - 1] = malloc(BlokSize * sizeof(uint32_t));
      if (psi->Next[psi->NBlok - 1] == NULL)
      {                         /* Allocation failed */
        printf("alloc failed RealWaveFunctionAddElement 05\n");
        fflush(stdout);
        /* restore psi and return 0 */
        free(psi->ValueRe[psi->NBlok - 1]);
        psi->NBlok--;
        psi->N--;
        psi->ValueRe = realloc(psi->ValueRe, psi->NBlok * sizeof(double *));
        psi->Next = realloc(psi->Next, psi->NBlok * sizeof(uint32_t *));
        psi->Determinant = realloc(psi->Determinant, psi->NBlok * sizeof(uint8_t *));
        return (1);
      }
      psi->Determinant[psi->NBlok - 1] = malloc(BlokSize * psi->NWord * sizeof(uint8_t));
      if (psi->Determinant[psi->NBlok - 1] == NULL)
      {                         /* Allocation failed */
        printf("alloc failed RealWaveFunctionAddElement 06\n");
        fflush(stdout);
        /* restore psi and return 0 */
        free(psi->ValueRe[psi->NBlok - 1]);
        free(psi->Next[psi->NBlok - 1]);
        psi->NBlok--;
        psi->N--;
        psi->ValueRe = realloc(psi->ValueRe, psi->NBlok * sizeof(double *));
        psi->Next = realloc(psi->Next, psi->NBlok * sizeof(uint32_t *));
        psi->Determinant = realloc(psi->Determinant, psi->NBlok * sizeof(uint8_t *));
        return (1);
      }
    }                           /* Added a new block to psi */
    psi->ValueRe[BlokElement(psi->N)] = a;
    psi->Next[BlokElement(psi->N)] = 0;
    for (i = 0; i < psi->NWord; i++)
    {
      psi->Determinant[BlokElementDet(psi->N, psi->NWord, i)] = det[i];
    }
    psi->KeyList[hash(det, psi->NWord, psi->NBitsKey)] = psi->N;
  }
  else
  {                             /* The key is already taken, test if this is due to a collision, or if psi is already in the basis set */
    ThereIsANext = 1;
    NextIndex = psi->KeyList[hash(det, psi->NWord, psi->NBitsKey)];
    while (ThereIsANext)
    {
      PsiAIsPsiB = 1;
      for (i = 0; i < psi->NWord; i++)
      {
        if (psi->Determinant[BlokElementDet(NextIndex, psi->NWord, i)] != det[i])
          PsiAIsPsiB = 0;
      }
      if (PsiAIsPsiB)
      {                         /* psi already in basis set add and leave loop  */
        ThereIsANext = 0;
        psi->ValueRe[BlokElement(NextIndex)] += a;
      }
      else
      {                         /* psi not equal to this one, if there is a next continue, otherwise add */
        if (!psi->Next[BlokElement(NextIndex)])
        {                       /* add a new element to psi with already used key */
          psi->N++;
          /* test if one has to add a new block to psi in order to store the new element of psi */
          if (psi->N >= BlokSize * psi->NBlok)
          {                     /* add a new block to psi */
            psi->NBlok++;
            temp = psi->ValueRe;
            psi->ValueRe = realloc(psi->ValueRe, psi->NBlok * sizeof(double *));
            if (psi->ValueRe == NULL)
            {                   /* Allocation failed */
              printf("alloc failed RealWaveFunctionAddElement 11\n");
              fflush(stdout);
              /* restore psi and return 1 */
              psi->ValueRe = temp;
              psi->NBlok--;
              psi->N--;
              return (1);
            }
            temp = psi->Next;
            psi->Next = realloc(psi->Next, psi->NBlok * sizeof(uint32_t *));
            if (psi->Next == NULL)
            {                   /* Allocation failed */
              printf("alloc failed RealWaveFunctionAddElement 12\n");
              fflush(stdout);
              /* restore psi and return 1 */
              psi->Next = temp;
              psi->NBlok--;
              psi->N--;
              psi->ValueRe = realloc(psi->ValueRe, psi->NBlok * sizeof(double *));
              return (1);
            }
            temp = psi->Determinant;
            psi->Determinant = realloc(psi->Determinant, psi->NBlok * sizeof(uint8_t *));
            if (psi->Determinant == NULL)
            {                   /* Allocation failed */
              printf("alloc failed RealWaveFunctionAddElement 13\n");
              fflush(stdout);
              /* restore psi and return 1 */
              psi->Determinant = temp;
              psi->NBlok--;
              psi->N--;
              psi->ValueRe = realloc(psi->ValueRe, psi->NBlok * sizeof(double *));
              psi->Next = realloc(psi->Next, psi->NBlok * sizeof(uint32_t *));
              return (1);
            }
            psi->ValueRe[psi->NBlok - 1] = malloc(BlokSize * sizeof(double));
            if (psi->ValueRe[psi->NBlok - 1] == NULL)
            {                   /* Allocation failed */
              printf("alloc failed RealWaveFunctionAddElement 14\n");
              fflush(stdout);
              /* restore psi and return 1 */
              psi->NBlok--;
              psi->N--;
              psi->ValueRe = realloc(psi->ValueRe, psi->NBlok * sizeof(double *));
              psi->Next = realloc(psi->Next, psi->NBlok * sizeof(uint32_t *));
              psi->Determinant = realloc(psi->Determinant, psi->NBlok * sizeof(uint8_t *));
              return (1);
            }
            psi->Next[psi->NBlok - 1] = malloc(BlokSize * sizeof(uint32_t));
            if (psi->Next[psi->NBlok - 1] == NULL)
            {                   /* Allocation failed */
              printf("alloc failed RealWaveFunctionAddElement 15\n");
              fflush(stdout);
              /* restore psi and return 1 */
              free(psi->ValueRe[psi->NBlok - 1]);
              psi->NBlok--;
              psi->N--;
              psi->ValueRe = realloc(psi->ValueRe, psi->NBlok * sizeof(double *));
              psi->Next = realloc(psi->Next, psi->NBlok * sizeof(uint32_t *));
              psi->Determinant = realloc(psi->Determinant, psi->NBlok * sizeof(uint8_t *));
              return (1);
            }
            psi->Determinant[psi->NBlok - 1] = malloc(BlokSize * psi->NWord * sizeof(uint8_t));
            if (psi->Determinant[psi->NBlok - 1] == NULL)
            {                   /* Allocation failed */
              printf("alloc failed RealWaveFunctionAddElement 16\n");
              fflush(stdout);
              /* restore psi and return 1 */
              free(psi->ValueRe[psi->NBlok - 1]);
              free(psi->Next[psi->NBlok - 1]);
              psi->NBlok--;
              psi->N--;
              psi->ValueRe = realloc(psi->ValueRe, psi->NBlok * sizeof(double *));
              psi->Next = realloc(psi->Next, psi->NBlok * sizeof(uint32_t *));
              psi->Determinant = realloc(psi->Determinant, psi->NBlok * sizeof(uint8_t *));
              return (1);
            }
          }                     /* Added a new block to psi */
          /* add the element */
          psi->ValueRe[BlokElement(psi->N)] = a;
          psi->Next[BlokElement(psi->N)] = 0;
          for (i = 0; i < psi->NWord; i++)
          {
            psi->Determinant[BlokElementDet(psi->N, psi->NWord, i)] = det[i];
          }
          psi->Next[BlokElement(NextIndex)] = psi->N;
          ThereIsANext = 0;
        }
        else
        {                       /* the key is in use and there still is at least one more element with the same key, continu checking */
          NextIndex = psi->Next[BlokElement(NextIndex)];
        }
      }
    }
  }
  return (0);
}


/*
*   add element with single slater determinant det and prefactor a to wave function psi
*   But only add it iff det is already present in the basis
*   As psi only stores the number of non zero elements and pointers to the real values the value of psi is not changed.
*   The values psi points to are naturaly changed
*/
void RealWaveFunctionAddElementConserveBasis(WaveFunctionType const *psi, const uint8_t * const det, const double a)
{
  uint32_t i, NextIndex, PsiAIsPsiB;
  /* test if the key that belongs to psi already excists */
  if (a < 10*DBL_MIN) return;
  NextIndex=psi->KeyList[hash(det, psi->NWord, psi->NBitsKey)];
  if (!NextIndex)
  {                             /* det not in key, do nothing */
    return;
  }
  else
  {                             /* The key is already taken, test if this is due to a collision, or if psi is already in the basis set */
    while (1)
    {
      PsiAIsPsiB = 1;
      for (i = 0; i < psi->NWord; i++)
      {
        if (psi->Determinant[BlokElementDet(NextIndex, psi->NWord, i)] != det[i])
        {
          PsiAIsPsiB = 0;
        }
      }
      if (PsiAIsPsiB)
      {                         /* psi already in basis set add and leave loop  */
#pragma omp atomic
        psi->ValueRe[BlokElement(NextIndex)] += a;
        return;
      }
      else
      {                         /* psi not equal to this one, if there is a next continue, otherwise leave */
        if (!psi->Next[BlokElement(NextIndex)])
        {                       /* hash(psi) had a collision, but det is different. There are no other psi's with the same hash, so det not in basis lets leave */
          return;
        }
        else
        {                       /* the key is in use and there still is at least one more element with the same key, continu checking */
          NextIndex = psi->Next[BlokElement(NextIndex)];
        }
      }
    }
  }
}


/*  return the prefactor of det in wave function psi
 *
 */
double RealWaveFunctionPrefactorDet(const WaveFunctionType psi, const uint8_t * const det)
{
  uint32_t i, NextIndex, PsiAIsPsiB;
  /* test if the key that belongs to psi already excists */
  if (psi.KeyList[hash(det, psi.NWord, psi.NBitsKey)] == 0)
  {                             /* det not in key return 0 */
    return 0.0;
  }
  else
  {                             /* The key is already taken, test if this is due to a collision, or if psi is already in the basis set */
    NextIndex = psi.KeyList[hash(det, psi.NWord, psi.NBitsKey)];
    while (1)
    {
      PsiAIsPsiB = 1;
      for (i = 0; i < psi.NWord; i++)
      {
        if (psi.Determinant[BlokElementDet(NextIndex, psi.NWord, i)] != det[i]) PsiAIsPsiB = 0;
      }
      if (PsiAIsPsiB)
      {                         /* psi already in basis set return prefactor and leave loop  */
        return psi.ValueRe[BlokElement(NextIndex)];
      }
      else
      {                         /* psi not equal to this one, if there is a next continue, otherwise leave */
        if (!psi.Next[BlokElement(NextIndex)])
        {                       /* hash(psi) had a collision, but det is different. There are no other psi's with the same hash, so det not in basis lets leave */
          return 0.0;
        }
        else
        {                       /* the key is in use and there still is at least one more element with the same key, continu checking */
          NextIndex = psi.Next[BlokElement(NextIndex)];
        }
      }
    }
  }
}

/*  Parallel calculating CPU's can't write simultanious to the same wave function.
*   therefor they have a small buffer in wich elemets to add are stored and when this buffer is full the entire buffer is added to the wave function
*   This has to be done in a critical region, i.e. only one CPU at a time can do this routine
*   return 0 if everything worked
*   return 1 if out off memory
*/
uint32_t RealWaveFunctionAddElementOMPMiniFlush(WaveFunctionType * const psi, uint32_t * NMini,
                                                const uint8_t * const MiniDet, const double *const MiniValue)
{
  uint32_t i, NoMoreMem = 0;
#pragma omp critical (AddToPsi)
  {
    for (i = 0; i < *NMini; i++)
    {
      if (RealWaveFunctionAddElement(psi, &MiniDet[i * psi->NWord], MiniValue[i]))
      {
        printf("RealWaveFunctionAddElement failed in RealWaveFunctionAddElementOMPMiniFlush\n");
        fflush(stdout);
        NoMoreMem = 1;
        i = *NMini;
        *NMini = 0;
      }
    }
#pragma omp flush
  }
  *NMini = 0;
  return NoMoreMem;
}


/*
*   a->a+b
*   Look out this function does not alter the values stored in a nor b
*   As a and b only store pointers
*   This function does alter the values a points to
*/
void RealWaveFunctionAddSameBasis(WaveFunctionType a, const WaveFunctionType b)
{
  uint32_t iBlok, iElement;
  int    iOMP;
  iBlok = a.NBlok - 1;
#pragma omp parallel for private(iElement)
  for (iOMP = 0; iOMP <= (int)(a.N - iBlok * BlokSize); ++iOMP)
  {
    iElement = (uint32_t) iOMP;
    a.ValueRe[iBlok][iElement] += b.ValueRe[iBlok][iElement];
  }
#pragma omp parallel for private(iElement,iBlok)
  for (iOMP = 0; iOMP < (int)(a.NBlok - 1); iOMP++)
  {
    iBlok = (uint32_t) iOMP;
    for (iElement = 0; iElement < BlokSize; ++iElement)
    {
      a.ValueRe[iBlok][iElement] += b.ValueRe[iBlok][iElement];
    }
  }
}


/*
*   a->a+alpha b
*   Look out this function does not alter the values stored in a nor b
*   As a and b only store pointers
*   This function does alter the values a points to
*/
void RealWaveFunctionAddMultipliedSameBasis(WaveFunctionType a, const WaveFunctionType b, const double alpha)
{
  uint32_t iBlok, iElement;
  int    iOMP;
  iBlok = a.NBlok - 1;
#pragma omp parallel for private(iElement)
  for (iOMP = 0; iOMP <= (int)(a.N - iBlok * BlokSize); ++iOMP)
  {
    iElement = (uint32_t) iOMP;
    a.ValueRe[iBlok][iElement] += alpha * b.ValueRe[iBlok][iElement];
  }
#pragma omp parallel for private(iElement,iBlok)
  for (iOMP = 0; iOMP < (int)(a.NBlok - 1); iOMP++)
  {
    iBlok = (uint32_t) iOMP;
    for (iElement = 0; iElement < BlokSize; ++iElement)
    {
      a.ValueRe[iBlok][iElement] += alpha * b.ValueRe[iBlok][iElement];
    }
  }
}

/*
*   (a+new allocation) <- a+b
*   This function leaves b in the memory and adds values to a
*   When aditional space for a needs to be allocated it will
*   Therefor values of a might change and a is returned from the function
*   If a is an empty function b is copied into a
*   return 0 if succeded
*   return 1 if not enough memory in order to allocate a
*   If not enough memory a is not changed, but might have a larger basis
*/
uint32_t RealWaveFunctionAddUnionBasis(WaveFunctionType * const a, const WaveFunctionType b)
{
  uint32_t i, j;
  for (i = b.N; i != 0; --i)
  {
    if (RealWaveFunctionAddElement(a, &b.Determinant[BlokElementDet(i, b.NWord, 0)], b.ValueRe[BlokElement(i)]) != 0)
    {                           /* Allocation of new memory for a failed Remove all elements added so far and return 1 */
      for (j = b.N; j != i; --j)
      {
        RealWaveFunctionAddElementConserveBasis(a, &b.Determinant[BlokElementDet(i, b.NWord, 0)],-b.ValueRe[BlokElement(i)]);
      }
      return 1;
    }
  }
  return 0;
}

/*
*   (a+new allocation) <- a+alpha b
*   This function leaves b in the memory and adds values to a
*   When aditional space for a needs to be allocated it will
*   Therefor values of a might change and a is returned from the function
*   If a is an empty function b is copied into a (use of copy is strongly recommanded)
*/
uint32_t RealWaveFunctionAddMultipliedUnionBasis(WaveFunctionType * const a, const WaveFunctionType b, const double alpha)
{
  uint32_t i, j;
  for (i = b.N; i != 0; --i)
  {
    if (RealWaveFunctionAddElement(a, &b.Determinant[BlokElementDet(i, b.NWord, 0)], alpha * b.ValueRe[BlokElement(i)])
        != 0)
    {                           /* Allocation of new memory for a failed Remove all elements added so far and return 1 */
      for (j = b.N; j != i; --j)
      {
        RealWaveFunctionAddElementConserveBasis(a, &b.Determinant[BlokElementDet(i, b.NWord, 0)],-alpha * b.ValueRe[BlokElement(i)]);
      }
      return 1;
    }
  }
  return 0;
}


/*
 *   (a) <- a+ b
 *   This function leaves b in the memory and adds values to a
 *   if b_i is not in the basis of a it will be ignored
 *   Therefor values of a might change and a is returned from the function
 *   If a is an empty function b is copied into a (use of copy is strongly recommanded)
 */
void RealWaveFunctionAddConserveBasis(WaveFunctionType * const a, const WaveFunctionType b)
{
  uint32_t i;
  int    iOMP;
#pragma omp parallel for private(i)
  for (iOMP = 1; iOMP <= (int)b.N; iOMP++)
  {
    i = (uint32_t) iOMP;
    RealWaveFunctionAddElementConserveBasis(a, &b.Determinant[BlokElementDet(i, b.NWord, 0)],b.ValueRe[BlokElement(i)]);
  }
}

/*
*   (a) <- a+alpha b
*   This function leaves b in the memory and adds values to a
*   if b_i is not in the basis of a it will be ignored
*   Therefor values of a might change and a is returned from the function
*   If a is an empty function b is copied into a (use of copy is strongly recommanded)
*/
void RealWaveFunctionAddMultipliedConserveBasis(WaveFunctionType * const a, const WaveFunctionType b,const double alpha)
{
  uint32_t i;
  int    iOMP;
#pragma omp parallel for private(i)
  for (iOMP = 1; iOMP <= (int)b.N; iOMP++)
  {
    i = (uint32_t) iOMP;
    RealWaveFunctionAddElementConserveBasis(a, &b.Determinant[BlokElementDet(i, b.NWord, 0)],alpha * b.ValueRe[BlokElement(i)]);
  }
}


/*
*   a->a/|a|
*   Look out this function does not alter the values stored in a
*   As a only stores pointers
*   This function does alter the values a points to
*/
void RealWaveFunctionNormalize(WaveFunctionType * const a)
{
  uint32_t iBlok, iElement;
  int    iOMP;
  double NormA = 0, NormB = 0, Norm = 0;
  volatile double c, t, y;
  a->ValueRe[0][0] = 0;
  iBlok = a->NBlok - 1;
  c = 0;
#pragma omp parallel for firstprivate(c) private(iElement,t,y) reduction(+ : NormA)
  for (iOMP = 0; iOMP <= (int)(a->N - iBlok * BlokSize); iOMP++)
  {
    iElement = (uint32_t) iOMP;
    y = a->ValueRe[iBlok][iElement] * a->ValueRe[iBlok][iElement] - c;
    t = NormA + y;
    c = t - NormA;
    c -= y;
    NormA = t;
  }
  c = 0;
#pragma omp parallel for firstprivate(c) private(iElement,iBlok,t,y) reduction(+ : NormB)
  for (iOMP = 0; iOMP < (int)(a->NBlok - 1); iOMP++)
  {
    iBlok = (uint32_t) iOMP;
    for (iElement = 0; iElement < BlokSize; ++iElement)
    {
      y = a->ValueRe[iBlok][iElement] * a->ValueRe[iBlok][iElement] - c;
      t = NormB + y;
      c = t - NormB;
      c -= y;
      NormB = t;
    }
  }
  Norm = sqrt(NormA + NormB);
  iBlok = a->NBlok - 1;
  if(Norm<10*DBL_MIN) return;
#pragma omp parallel for private(iElement)
  for (iOMP = 0; iOMP <= (int)(a->N - iBlok * BlokSize); iOMP++)
  {
    iElement = (uint32_t) iOMP;
    a->ValueRe[iBlok][iElement] /= Norm;
  }
#pragma omp parallel for private(iElement,iBlok)
  for (iOMP = 0; iOMP < (int)(a->NBlok - 1); iOMP++)
  {
    iBlok = (uint32_t) iOMP;
    for (iElement = 0; iElement < BlokSize; ++iElement)
    {
      a->ValueRe[iBlok][iElement] /= Norm;
    }
  }
}

/*
*   a->a*b
*   Look out this function does not alter the values stored in a
*   As a only stores pointers
*   This function does alter the values a points to
*/
void RealWaveFunctionMultiply(WaveFunctionType a, const double b)
{
  uint32_t iBlok, iElement;
  int    iOMP;
  iBlok = a.NBlok - 1;
#pragma omp parallel for private(iElement)
  for (iOMP = 0; iOMP <= (int)(a.N - iBlok * BlokSize); ++iOMP)
  {
    iElement = (uint32_t) iOMP;
    a.ValueRe[iBlok][iElement] *= b;
  }
#pragma omp parallel for private(iElement,iBlok)
  for (iOMP = 0; iOMP < (int)(a.NBlok - 1); iOMP++)
  {
    iBlok = (uint32_t) iOMP;
    for (iElement = 0; iElement < BlokSize; ++iElement)
    {
      a.ValueRe[iBlok][iElement] *= b;
    }
  }
}


/*
*   <-|a|
*/
double RealWaveFunctionCalculateNorm(const WaveFunctionType a)
{
  uint32_t iBlok, iElement;
  volatile double c, t, y;
  int    iOMP;
  double NormA = 0, NormB = 0, Norm = 0;
  iBlok = a.NBlok - 1;
  c = 0;
#pragma omp parallel for firstprivate(c) private(iElement,t,y) reduction(+ : NormA)
  for (iOMP = 0; iOMP <= (int)(a.N - iBlok * BlokSize); ++iOMP)
  {
    iElement = (uint32_t) iOMP;
    y = a.ValueRe[iBlok][iElement] * a.ValueRe[iBlok][iElement] - c;
    t = NormA + y;
    c = t - NormA;
    c -= y;
    NormA = t;
  }
  c = 0;
#pragma omp parallel for firstprivate(c) private(iElement,iBlok,t,y) reduction(+ : NormB)
  for (iOMP = 0; iOMP < (int)(a.NBlok - 1); iOMP++)
  {
    iBlok = (uint32_t) iOMP;
    for (iElement = 0; iElement < BlokSize; ++iElement)
    {
      y = a.ValueRe[iBlok][iElement] * a.ValueRe[iBlok][iElement] - c;
      t = NormB + y;
      c = t - NormB;
      c -= y;
      NormB = t;
    }
  }
  Norm = sqrt(NormA + NormB);
  return Norm;
}

/*
*   calculate <a|b>
*/
double RealWaveFunctionDotProductSameBasis(const WaveFunctionType a, const WaveFunctionType b)
{
  uint32_t iBlok, iElement;
  int    iOMP;
  double DotA = 0, DotB = 0, Dot = 0;
  volatile double c, t, y;
  iBlok = a.NBlok - 1;
  c = 0;
#pragma omp parallel for firstprivate(c) private(iElement,t,y) reduction(+ : DotA)
  for (iOMP = 0; iOMP <= (int)(a.N - iBlok * BlokSize); ++iOMP)
  {
    iElement = (uint32_t) iOMP;
    y = a.ValueRe[iBlok][iElement] * b.ValueRe[iBlok][iElement] - c;
    t = DotA + y;
    c = t - DotA;
    c -= y;
    DotA = t;
  }
  c = 0;
#pragma omp parallel for firstprivate(c) private(iElement,iBlok,t,y) reduction(+ : DotB)
  for (iOMP = 0; iOMP < (int)(a.NBlok - 1); iOMP++)
  {
    iBlok = (uint32_t) iOMP;
    for (iElement = 0; iElement < BlokSize; ++iElement)
    {
      y = a.ValueRe[iBlok][iElement] * b.ValueRe[iBlok][iElement] - c;
      t = DotB + y;
      c = t - DotB;
      c -= y;
      DotB = t;
    }
  }
  Dot = DotA + DotB;
  return Dot;
}

/*
*   calculate <a|b>
*/
double RealWaveFunctionDotProduct(const WaveFunctionType *a, const WaveFunctionType *b)
{
  uint32_t i, j, NextIndex, PsiAIsPsiB;
  int    iOMP;
  double Dot = 0;
  volatile double c, t, y;
  c = 0;

  /* switch coefficients to reduce number of hash lookups */
  if (a->N > b->N)
  {
    const WaveFunctionType *temp = a;
    a = b;
    b = temp;
  }

#pragma omp parallel for firstprivate(c) private(i,j,NextIndex,PsiAIsPsiB,t,y) reduction(+ : Dot)
  for (iOMP = 1; iOMP <= (int)(a->N); ++iOMP)
  {
    i = (uint32_t) iOMP;
    NextIndex = b->KeyList[hash(&a->Determinant[BlokElementDet(i, a->NWord, 0)], a->NWord, a->NBitsKey)];
    while (NextIndex)
    {
      PsiAIsPsiB = 1;
      for (j = 0; j < b->NWord; j++)
      {
        if (b->Determinant[BlokElementDet(NextIndex, b->NWord, j)] != a->Determinant[BlokElementDet(i, a->NWord, j)])
        {
          PsiAIsPsiB = 0;
        }
      }
      if (PsiAIsPsiB)
      {                         /* a=b found add the product to Dot and leave loop  */
        y = a->ValueRe[BlokElement(i)] * b->ValueRe[BlokElement(NextIndex)] - c;
        t = Dot + y;
        c = t - Dot;
        c -= y;
        Dot = t;
        NextIndex = 0;
      }
      else
      {
        NextIndex = b->Next[BlokElement(NextIndex)];
      }
    }
  }
  return Dot;
}


/* Set the pre-factors of a to a random value
*/
void RealWaveFunctionRandomize(WaveFunctionType * const a)
{
  uint32_t iBlok, iElement;
  int    iOMP;
  iBlok = a->NBlok - 1;
#pragma omp parallel for private(iElement)
  for (iOMP = 0; iOMP <= (int)(a->N - iBlok * BlokSize); ++iOMP)
  {
    iElement = (uint32_t) iOMP;
    a->ValueRe[iBlok][iElement] = rand() - RAND_MAX / 2;
  }
#pragma omp parallel for private(iElement,iBlok)
  for (iOMP = 0; iOMP < (int)(a->NBlok - 1); iOMP++)
  {
    iBlok = (uint32_t) iOMP;
    for (iElement = 0; iElement < BlokSize; ++iElement)
    {
      a->ValueRe[iBlok][iElement] = rand() - RAND_MAX / 2;
    }
  }
  a->ValueRe[0][0] = 0;
  RealWaveFunctionNormalize(a);
}

/*
*   remove determinants from the basis that have a prefactor smaller then epsilon  
*/
void RealWaveFunctionTruncateBasis(WaveFunctionType * const a, const double epsilon)
{
  uint32_t startN, iBlok, iElement;
  a->ValueRe[0][0] = 0;
  startN = a->N;
  a->N = 0;
  for (iElement = 0; iElement < ((uint32_t) 1 << (a->NBitsKey)); ++iElement)
  {                             /* Whipe out the entire KeyList */
    a->KeyList[iElement] = 0;
  }
  for (iBlok = 0; iBlok < a->NBlok - 1; iBlok++)
  {
    for (iElement = 0; iElement < BlokSize; iElement++)
    {
      a->Next[iBlok][iElement] = 0;
      if ((a->ValueRe[iBlok][iElement] > epsilon) || (-a->ValueRe[iBlok][iElement] > epsilon))
        if (RealWaveFunctionAddElement(a, &a->Determinant[iBlok][iElement * a->NWord], a->ValueRe[iBlok][iElement]) != 0)
          AllocFailed();
    }
  }
  iBlok = a->NBlok - 1;
  for (iElement = 0; iElement < 1 + startN - iBlok * BlokSize; ++iElement)
  {
    a->Next[iBlok][iElement] = 0;
    if ((a->ValueRe[iBlok][iElement] > epsilon) || (-a->ValueRe[iBlok][iElement] > epsilon))
      if (RealWaveFunctionAddElement(a, &a->Determinant[iBlok][iElement * a->NWord], a->ValueRe[iBlok][iElement]) != 0)
        AllocFailed();
  }
  for (iBlok = 1 + (a->N >> NBitsBlok); iBlok < a->NBlok; iBlok++)
  {
    free(a->ValueRe[iBlok]);
    free(a->Next[iBlok]);
    free(a->Determinant[iBlok]);
  }
  a->NBlok = 1 + (a->N >> NBitsBlok);
  a->ValueRe = realloc(a->ValueRe, a->NBlok * sizeof(double *));
  a->Next = realloc(a->Next, a->NBlok * sizeof(uint32_t *));
  a->Determinant = realloc(a->Determinant, a->NBlok * sizeof(uint8_t *));
}

/*
 *   Make psi orthogonal to OtherPsi, psi and other psi are on the same basis
 */
void RealWaveFunctionOrthogonalizeSameBasis(WaveFunctionType * const psi, const WaveFunctionType * const OtherPsi, const uint32_t Npsi)
{
  uint32_t i;
  double a;
  for (i = 0; i < Npsi; i++)
  {
    a = RealWaveFunctionDotProductSameBasis(*psi, OtherPsi[i]);
    RealWaveFunctionAddMultipliedSameBasis(*psi, OtherPsi[i], -a);
  }
}

/*  
*   Make psi orthogonal to OtherPsi, but conserve the basis set of psi
*/
void RealWaveFunctionOrthogonalizeConserveBasis(WaveFunctionType * const psi, const WaveFunctionType * const OtherPsi, const uint32_t Npsi)
{
  uint32_t i;
  double a;
  for (i = 0; i < Npsi; i++)
  {
    a = RealWaveFunctionDotProduct(psi, &OtherPsi[i]);
    RealWaveFunctionAddMultipliedConserveBasis(psi, OtherPsi[i], -a);
  }
}

/*  
*   Make psi orthogonal to OtherPsi, if needed enlarge the basis set of psi
*/
uint32_t RealWaveFunctionOrthogonalize(WaveFunctionType * const psi, const WaveFunctionType * const OtherPsi, const uint32_t Npsi)
{
  uint32_t i;
  double a;
  for (i = 0; i < Npsi; i++)
  {
    a = RealWaveFunctionDotProduct(psi, &OtherPsi[i]);
    if (RealWaveFunctionAddMultipliedUnionBasis(psi, OtherPsi[i], -a))
    {
      printf("RealWaveFunctionAddMultipliedUnionBasis failed in RealWaveFunctionOrthogonalize\n");
      return 1;
    }
  }
  return 0;
}


uint32_t CreateRealRandomWaveFunction(WaveFunctionType * const psi, uint32_t const NFermions, uint32_t const nelectrons)
{
  if (nelectrons > NFermions)
  {
    printf("Error in CreateRealRandomWaveFunction\n");
    printf("You tried to make a funciton with %i electrons in %i places\n", nelectrons, NFermions);
    return 1;
  };
  uint32_t idet, ielec, iferm, iofset;
  uint8_t *Determinant;
  psi->NBitsKey = 16;
  psi->QComplex = 0;
  psi->NBoson = 0;
  psi->NFermion = NFermions;
  if (WaveFunctionInit(psi))
  {
    printf("WaveFunctionInit failed in CreateRealRandomWaveFunction\n");
    fflush(stdout);
    return 1;
  };
  sprintf(psi->Name, "WaveFunction");
  for (idet = 0; idet < 50; idet++)
  {
    Determinant = calloc(psi->NWord, sizeof(uint8_t));
    iofset = 0;
    for (ielec = 0; ielec < nelectrons; ielec++)
    {
      if (rand() & 1)
      {
        if (NFermions - nelectrons > iofset)
          iofset++;
        if (rand() & 1)
        {
          if (NFermions - nelectrons > iofset)
            iofset++;
          if (rand() & 1)
          {
            if (NFermions - nelectrons > iofset)
              iofset++;
            if (rand() & 1)
            {
              if (NFermions - nelectrons > iofset)
                iofset++;
              if (rand() & 1)
              {
                if (NFermions - nelectrons > iofset)
                  iofset++;
                if (rand() & 1)
                {
                  if (NFermions - nelectrons > iofset)
                    iofset++;
                }
              }
            }
          }
        }
      }
      iferm = ielec + iofset;
      OperatorCr(Determinant, iferm, psi->NFermion);
    }
    RealWaveFunctionAddElement(psi, Determinant, 1.0);
    free(Determinant);
  }
  RealWaveFunctionRandomize(psi);
  return 0;
}



/* on entry 
 *  A is an empty wave function with allocated memory (init A is done)
 *  Det contains the determinant to add
 *  U the rotation matrix
 *  Elements mixing bosons and fermions in U are neglected
 *  U is assumed to be square and unitary
 *  Det is assumed to have A->NFermion fermions and A->NBoson bosons
 *
 */
uint32_t RealMappedDetToWaveFunction(WaveFunctionType * const A, uint8_t *Det, MatrixType const U)
{
  WaveFunctionType Temp;
  uint8_t DetZero[A->NWord];
  uint16_t oadd[1];
  uint32_t i,j;
  OperatorType O;
  InitOperator(&O);
  O.NFermion=A->NFermion;
  O.NBoson=A->NBoson;
  O.NWord=A->NWord;
  
  InitOperatorFixedLength(&O, 1, 8);
  
  
  for(i=0;i<A->NWord;i++)
  {
    DetZero[i]=0;
  }
  
  RealWaveFunctionAddElement(A, DetZero, 1.0);
  
  for(i=0;i<A->NFermion;i++)
  {
    oadd[0]=CrAnMask | (uint16_t)i;
    RealOperatorAddLader(&O, oadd, 1, 0.0);
  }
  
  for(i=0;i<A->NFermion;i++)
  {
    /* test if Fermion i is occupied */
    if(Det[i >> 3] & (128 >> (i & 7)))
    {
      for(j=0;j<A->NFermion;j++)
      {
        O.l[1].VRe[j]=U.VRe[i][j];
      }
      OperatorPsiRR(O, *A, &Temp);
      RealWaveFunctionTruncateBasis(&Temp, 10*DBL_EPSILON);
      WaveFunctionFree(A);
      *A=Temp;
    }
    
  }
  for(i=A->NFermion;i<A->NFermion+A->NBoson;i++)
  {
    printf("Boson rotation not yet implemted\n");
    return 1;
  }
  
  FreeOperator(&O);
  return 0;
}


/* psiA <- psiA + U psiB
 *  U is a (NFermionB+NBosonB)*(NFermionA+NBosonA) demensional matrix that defines the transformation
 *  a_j        = \sum_i U_{i,j}   a_i
 *  a_j^{\dag} = \sum_i U_{i,j}^* a_i^{dag}
 *
 * assume A B and U are real (version RRR)
 */
uint32_t MappedWaveFunctionAddRRR(WaveFunctionType * const A, WaveFunctionType const B, MatrixType const U)
{
  uint8_t *Det;
  uint32_t idetB;
  WaveFunctionType Udet;
  
  
  Udet.NFermion=A->NFermion;
  Udet.NBoson=A->NBoson;
  Udet.QComplex=0;
  Udet.NBitsKey=A->NBitsKey;
  WaveFunctionInit(&Udet);
  
  /* loop over all determinants in B */
  for(idetB=1;idetB<=B.N;idetB++)
  {
    Det=&(B.Determinant[BlokElementDet(idetB, B.NWord, 0)]);
    RealWaveFunctionReset(&Udet);
    RealMappedDetToWaveFunction(&Udet, Det, U);
    RealWaveFunctionAddMultipliedUnionBasis(A, Udet, B.ValueRe[BlokElement(idetB)]);
  }
  RealWaveFunctionTruncateBasis(A, 10*DBL_EPSILON);
  RealWaveFunctionNormalize(A);
  WaveFunctionFree(&Udet);
  
  return 0;
}

