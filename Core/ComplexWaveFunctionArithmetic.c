
/*  Remove Wave function from memory
*   On entry psi contains a wave function with psi.NBlok alocations
*   On excit all allocate memory is freed
*/
void ComplexWaveFunctionFree(WaveFunctionType * const psi)
{
  uint32_t i;
  for (i = 0; i < psi->NBlok; i++)
  {
    free(psi->ValueRe[i]);
    free(psi->ValueIm[i]);
    free(psi->Next[i]);
    free(psi->Determinant[i]);
  }
  free(psi->ValueRe);
  free(psi->ValueIm);
  free(psi->Next);
  free(psi->Determinant);
  free(psi->KeyList);
}

void ComplexWaveFunctionReset(WaveFunctionType * const psi)
{
  uint32_t i;
  for (i = 1; i < psi->NBlok; i++)
  {
    free(psi->ValueRe[i]);
    free(psi->ValueIm[i]);
    free(psi->Next[i]);
    free(psi->Determinant[i]);
  }
  psi->N = 0;
  psi->NBlok=1;
  memset(psi->KeyList, 0, ( 1 << (psi->NBitsKey)) * sizeof(uint32_t));
  psi->ValueRe[0][0] = 0.0;
  psi->ValueIm[0][0] = 0.0;
}


/*  Remove The values of a wave function from memory, but keep the determiants and the pointers to next determinants with the same Hash value
*   Usefull for more then one wave function defined on the same basis where psi.Next and psi.Determinant point to the same location as used by an other wave function
*   On entry psi contains a wave function with psi.NBlok alocations
*   On excit all allocate memory for psi.ValueRe is freed
*/
void ComplexWaveFunctionFreeValues(WaveFunctionType * const psi)
{
  uint32_t i;
  for (i = 0; i < psi->NBlok; i++)
  {
    free(psi->ValueRe[i]);
    free(psi->ValueIm[i]);
  }
  free(psi->ValueRe);
  free(psi->ValueIm);
}

/*
*   add ellement with single slater determinant det and pre-factor a to wave function psi
*   returns 0 if succeded
*   returns 1 if no memory available in order to add det to psi
*   If succeded psi -> psi+det
*   If not succeded psi -> psi
*/
uint32_t ComplexWaveFunctionAddElement(WaveFunctionType * const psi, const uint8_t * const det, const double aRe, const double aIm)
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
        printf("alloc failed RealWaveFunctionAddElement 01 Re\n");
        fflush(stdout);
        /* restore psi and return 0 */
        psi->ValueRe = temp;
        psi->NBlok--;
        return (1);
      }
      temp = psi->ValueIm;
      psi->ValueIm = realloc(psi->ValueIm, psi->NBlok * sizeof(double *));
      if (psi->ValueIm == NULL)
      {                         /* Allocation failed */
        printf("alloc failed RealWaveFunctionAddElement 01 Im\n");
        fflush(stdout);
        /* restore psi and return 0 */
        psi->ValueIm = temp;
        psi->NBlok--;
        psi->ValueRe = realloc(psi->ValueRe, psi->NBlok * sizeof(double *));
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
        psi->ValueRe = realloc(psi->ValueRe, psi->NBlok * sizeof(double *));
        psi->ValueIm = realloc(psi->ValueIm, psi->NBlok * sizeof(double *));
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
        psi->ValueRe = realloc(psi->ValueRe, psi->NBlok * sizeof(double *));
        psi->ValueIm = realloc(psi->ValueIm, psi->NBlok * sizeof(double *));
        psi->Next = realloc(psi->Next, psi->NBlok * sizeof(uint32_t *));
        return (1);
      }
      psi->ValueRe[psi->NBlok - 1] = malloc(BlokSize * sizeof(double));
      if (psi->ValueRe[psi->NBlok - 1] == NULL)
      {                         /* Allocation failed */
        printf("alloc failed RealWaveFunctionAddElement 04 Re\n");
        fflush(stdout);
        /* restore psi and return 0 */
        psi->NBlok--;
        psi->ValueRe = realloc(psi->ValueRe, psi->NBlok * sizeof(double *));
        psi->ValueIm = realloc(psi->ValueIm, psi->NBlok * sizeof(double *));
        psi->Next = realloc(psi->Next, psi->NBlok * sizeof(uint32_t *));
        psi->Determinant = realloc(psi->Determinant, psi->NBlok * sizeof(uint8_t *));
        return (1);
      }
      psi->ValueIm[psi->NBlok - 1] = malloc(BlokSize * sizeof(double));
      if (psi->ValueIm[psi->NBlok - 1] == NULL)
      {                         /* Allocation failed */
        printf("alloc failed RealWaveFunctionAddElement 04 Re\n");
        fflush(stdout);
        /* restore psi and return 0 */
        free(psi->ValueRe[psi->NBlok - 1]);
        psi->NBlok--;
        psi->ValueRe = realloc(psi->ValueRe, psi->NBlok * sizeof(double *));
        psi->ValueIm = realloc(psi->ValueIm, psi->NBlok * sizeof(double *));
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
        free(psi->ValueIm[psi->NBlok - 1]);
        psi->NBlok--;
        psi->ValueRe = realloc(psi->ValueRe, psi->NBlok * sizeof(double *));
        psi->ValueIm = realloc(psi->ValueIm, psi->NBlok * sizeof(double *));
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
        free(psi->ValueIm[psi->NBlok - 1]);
        free(psi->Next[psi->NBlok - 1]);
        psi->NBlok--;
        psi->ValueRe = realloc(psi->ValueRe, psi->NBlok * sizeof(double *));
        psi->ValueIm = realloc(psi->ValueIm, psi->NBlok * sizeof(double *));
        psi->Next = realloc(psi->Next, psi->NBlok * sizeof(uint32_t *));
        psi->Determinant = realloc(psi->Determinant, psi->NBlok * sizeof(uint8_t *));
        return (1);
      }
    }                           /* Added a new block to psi */
    psi->ValueRe[BlokElement(psi->N)] = aRe;
    psi->ValueIm[BlokElement(psi->N)] = aIm;
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
        psi->ValueRe[BlokElement(NextIndex)] += aRe;
        psi->ValueIm[BlokElement(NextIndex)] += aIm;
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
              printf("alloc failed RealWaveFunctionAddElement 11 Re\n");
              fflush(stdout);
              /* restore psi and return 1 */
              psi->ValueRe = temp;
              psi->NBlok--;
              return (1);
            }
            temp = psi->ValueIm;
            psi->ValueIm = realloc(psi->ValueIm, psi->NBlok * sizeof(double *));
            if (psi->ValueIm == NULL)
            {                   /* Allocation failed */
              printf("alloc failed RealWaveFunctionAddElement 11 Im\n");
              fflush(stdout);
              /* restore psi and return 1 */
              psi->ValueIm = temp;
              psi->NBlok--;
              psi->ValueRe = realloc(psi->ValueRe, psi->NBlok * sizeof(double *));
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
              psi->ValueRe = realloc(psi->ValueRe, psi->NBlok * sizeof(double *));
              psi->ValueIm = realloc(psi->ValueIm, psi->NBlok * sizeof(double *));
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
              psi->ValueRe = realloc(psi->ValueRe, psi->NBlok * sizeof(double *));
              psi->ValueIm = realloc(psi->ValueIm, psi->NBlok * sizeof(double *));
              psi->Next = realloc(psi->Next, psi->NBlok * sizeof(uint32_t *));
              return (1);
            }
            psi->ValueRe[psi->NBlok - 1] = malloc(BlokSize * sizeof(double));
            if (psi->ValueRe[psi->NBlok - 1] == NULL)
            {                   /* Allocation failed */
              printf("alloc failed RealWaveFunctionAddElement 14 Re\n");
              fflush(stdout);
              /* restore psi and return 1 */
              psi->NBlok--;
              psi->ValueRe = realloc(psi->ValueRe, psi->NBlok * sizeof(double *));
              psi->Next = realloc(psi->Next, psi->NBlok * sizeof(uint32_t *));
              psi->Determinant = realloc(psi->Determinant, psi->NBlok * sizeof(uint8_t *));
              return (1);
            }
            psi->ValueIm[psi->NBlok - 1] = malloc(BlokSize * sizeof(double));
            if (psi->ValueIm[psi->NBlok - 1] == NULL)
            {                   /* Allocation failed */
              printf("alloc failed RealWaveFunctionAddElement 14 Re\n");
              fflush(stdout);
              /* restore psi and return 1 */
              free(psi->ValueRe[psi->NBlok - 1]);
              psi->NBlok--;
              psi->ValueIm = realloc(psi->ValueIm, psi->NBlok * sizeof(double *));
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
              free(psi->ValueIm[psi->NBlok - 1]);
              psi->NBlok--;
              psi->ValueRe = realloc(psi->ValueRe, psi->NBlok * sizeof(double *));
              psi->ValueIm = realloc(psi->ValueIm, psi->NBlok * sizeof(double *));
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
              free(psi->ValueIm[psi->NBlok - 1]);
              free(psi->Next[psi->NBlok - 1]);
              psi->NBlok--;
              psi->ValueRe = realloc(psi->ValueRe, psi->NBlok * sizeof(double *));
              psi->ValueIm = realloc(psi->ValueIm, psi->NBlok * sizeof(double *));
              psi->Next = realloc(psi->Next, psi->NBlok * sizeof(uint32_t *));
              psi->Determinant = realloc(psi->Determinant, psi->NBlok * sizeof(uint8_t *));
              return (1);
            }
          }                     /* Added a new block to psi */
          /* add the element */
          psi->ValueRe[BlokElement(psi->N)] = aRe;
          psi->ValueIm[BlokElement(psi->N)] = aIm;
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
*   add ellement with single slater determinant det and prefactor a to wave function psi
*   But only add it iff det is already present in the basis
*   As psi only stores the number of non zero elements and pointers to the real values the value of psi is not changed.
*   The values psi points to are natturally changed
*/
void ComplexWaveFunctionAddElementConserveBasis(WaveFunctionType psi, const uint8_t * const det, const double aRe, const double aIm)
{
  uint32_t i, NextIndex, PsiAIsPsiB;
  if ((Sqr(aRe)+Sqr(aIm))<100*DBL_MIN) return;
  /* test if the key that belongs to psi already excists */
  NextIndex = psi.KeyList[hash(det, psi.NWord, psi.NBitsKey)];
  if (NextIndex == 0)
  {                             /* det not in key, do nothing */
    return;
  }
  else
  {                             /* The key is already taken, test if this is due to a collision, or if psi is already in the basis set */
    while (1)
    {
      PsiAIsPsiB = 1;
      for (i = 0; i < psi.NWord; i++)
      {
        if (psi.Determinant[BlokElementDet(NextIndex, psi.NWord, i)] != det[i])
        {
          PsiAIsPsiB = 0;
        }
      }
      if (PsiAIsPsiB)
      {                         /* psi already in basis set add and leave loop  */
#pragma omp atomic
        psi.ValueRe[BlokElement(NextIndex)] += aRe;
#pragma omp atomic
        psi.ValueIm[BlokElement(NextIndex)] += aIm;
        return;
      }
      else
      {                         /* psi not equal to this one, if there is a next continue, otherwise leave */
        if (!psi.Next[BlokElement(NextIndex)])
        {                       /* hash(psi) had a collision, but det is different. There are no other psi's with the same hash, so det not in basis lets leave */
          return;
        }
        else
        {                       /* the key is in use and there still is at least one more element with the same key, continu checking */
          NextIndex = psi.Next[BlokElement(NextIndex)];
        }
      }
    }
  }
}

/*
*   return the prefactor of det in wave function psi
*/
void ComplexWaveFunctionPrefactorDet(const WaveFunctionType psi, const uint8_t * const det, double *aRe, double *aIm)
{
  uint32_t i, NextIndex, PsiAIsPsiB;
  /* test if the key that belongs to psi already excists */
  if (psi.KeyList[hash(det, psi.NWord, psi.NBitsKey)] == 0)
  {                             /* det not in key retrun 0.0 */
    *aRe = 0.0;
    *aIm = 0.0;
    return;
  }
  else
  {                             /* The key is already taken, test if this is due to a collision, or if psi is already in the basis set */
    NextIndex = psi.KeyList[hash(det, psi.NWord, psi.NBitsKey)];
    while (1)
    {
      PsiAIsPsiB = 1;
      for (i = 0; i < psi.NWord; i++)
      {
        if (psi.Determinant[BlokElementDet(NextIndex, psi.NWord, i)] != det[i])
          PsiAIsPsiB = 0;
      }
      if (PsiAIsPsiB)
      {                         /* psi already in basis set retunr prefactor and leave loop  */
        *aRe = psi.ValueRe[BlokElement(NextIndex)];
        *aIm = psi.ValueIm[BlokElement(NextIndex)];
        return;
      }
      else
      {                         /* psi not equal to this one, if there is a next continue, otherwise leave */
        if (!psi.Next[BlokElement(NextIndex)])
        {                       /* hash(psi) had a collision, but det is different. There are no other psi's with the same hash, so det not in basis lets leave */
          *aRe = 0.0;
          *aIm = 0.0;
          return;
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
uint32_t ComplexWaveFunctionAddElementOMPMiniFlush(WaveFunctionType * const psi, uint32_t * NMini,
                                                   const uint8_t * const MiniDet, const double *const MiniValueRe,
                                                   const double *const MiniValueIm)
{
  uint32_t i, NoMoreMem = 0;
#pragma omp critical (AddToPsi)
  {
    for (i = 0; i < *NMini; i++)
    {
      if (ComplexWaveFunctionAddElement(psi, &MiniDet[i * psi->NWord], MiniValueRe[i], MiniValueIm[i]))
      {
        printf("ComplexWaveFunctionAddElement failed in ComplexWaveFunctionAddElementOMPMiniFlush\n");
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
void ComplexWaveFunctionAddSameBasis(WaveFunctionType a, const WaveFunctionType b)
{
  uint32_t iBlok, iElement;
  int    iOMP;
  iBlok = a.NBlok - 1;
#pragma omp parallel for private(iElement)
  for (iOMP = 0; iOMP <= (int)(a.N - iBlok * BlokSize); ++iOMP)
  {
    iElement = (uint32_t) iOMP;
    a.ValueRe[iBlok][iElement] += b.ValueRe[iBlok][iElement];
    a.ValueIm[iBlok][iElement] += b.ValueIm[iBlok][iElement];
  }
#pragma omp parallel for private(iElement,iBlok)
  for (iOMP = 0; iOMP < (int)(a.NBlok - 1); iOMP++)
  {
    iBlok = (uint32_t) iOMP;
    for (iElement = 0; iElement < BlokSize; ++iElement)
    {
      a.ValueRe[iBlok][iElement] += b.ValueRe[iBlok][iElement];
      a.ValueIm[iBlok][iElement] += b.ValueIm[iBlok][iElement];
    }
  }
}

/*
*   a->a+alpha b
*   Look out this function does not alter the values stored in a nor b
*   As a and b only store pointers
*   This function does alter the values a points to
*/
void ComplexWaveFunctionAddMultipliedSameBasis(WaveFunctionType a, const WaveFunctionType b, const double alpha)
{
  uint32_t iBlok, iElement;
  int    iOMP;
  iBlok = a.NBlok - 1;
#pragma omp parallel for private(iElement)
  for (iOMP = 0; iOMP <= (int)(a.N - iBlok * BlokSize); ++iOMP)
  {
    iElement = (uint32_t) iOMP;
    a.ValueRe[iBlok][iElement] += alpha * b.ValueRe[iBlok][iElement];
    a.ValueIm[iBlok][iElement] += alpha * b.ValueIm[iBlok][iElement];
  }
#pragma omp parallel for private(iElement,iBlok)
  for (iOMP = 0; iOMP < (int)(a.NBlok - 1); iOMP++)
  {
    iBlok = (uint32_t) iOMP;
    for (iElement = 0; iElement < BlokSize; ++iElement)
    {
      a.ValueRe[iBlok][iElement] += alpha * b.ValueRe[iBlok][iElement];
      a.ValueIm[iBlok][iElement] += alpha * b.ValueIm[iBlok][iElement];
    }
  }
}

/*
*   a->a+alpha b
*   Look out this function does not alter the values stored in a nor b
*   As a and b only store pointers
*   This function does alter the values a points to
*/
void ComplexWaveFunctionAddComplexMultipliedSameBasis(WaveFunctionType a, const WaveFunctionType b, const double alphaRe, const double alphaIm)
{
  uint32_t iBlok, iElement;
  int    iOMP;
  iBlok = a.NBlok - 1;
#pragma omp parallel for private(iElement)
  for (iOMP = 0; iOMP <= (int)(a.N - iBlok * BlokSize); ++iOMP)
  {
    iElement = (uint32_t) iOMP;
    a.ValueRe[iBlok][iElement] += alphaRe * b.ValueRe[iBlok][iElement] - alphaIm * b.ValueIm[iBlok][iElement];
    a.ValueIm[iBlok][iElement] += alphaRe * b.ValueIm[iBlok][iElement] + alphaIm * b.ValueRe[iBlok][iElement];
  }
#pragma omp parallel for private(iElement,iBlok)
  for (iOMP = 0; iOMP < (int)(a.NBlok - 1); iOMP++)
  {
    iBlok = (uint32_t) iOMP;
    for (iElement = 0; iElement < BlokSize; ++iElement)
    {
      a.ValueRe[iBlok][iElement] += alphaRe * b.ValueRe[iBlok][iElement] - alphaIm * b.ValueIm[iBlok][iElement];
      a.ValueIm[iBlok][iElement] += alphaRe * b.ValueIm[iBlok][iElement] + alphaIm * b.ValueRe[iBlok][iElement];
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
uint32_t ComplexWaveFunctionAddUnionBasis(WaveFunctionType * const a, const WaveFunctionType b)
{
  uint32_t i, j;
  for (i = b.N; i != 0; --i)
  {
    if (ComplexWaveFunctionAddElement
        (a, &b.Determinant[BlokElementDet(i, b.NWord, 0)], b.ValueRe[BlokElement(i)], b.ValueIm[BlokElement(i)]) != 0)
    {                           /* Allocation of new memory for a failed Remove all elements added so far and return 1 */
      for (j = b.N; j != i; --j)
      {
        ComplexWaveFunctionAddElementConserveBasis(*a, &b.Determinant[BlokElementDet(i, b.NWord, 0)],
                                                   -b.ValueRe[BlokElement(i)], -b.ValueIm[BlokElement(i)]);
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
uint32_t ComplexWaveFunctionAddMultipliedUnionBasis(WaveFunctionType * const a, const WaveFunctionType b, const double alphaRe)
{
  uint32_t i, j;
  for (i = b.N; i != 0; --i)
  {
    if (ComplexWaveFunctionAddElement
        (a, &b.Determinant[BlokElementDet(i, b.NWord, 0)], alphaRe * b.ValueRe[BlokElement(i)],
         alphaRe * b.ValueIm[BlokElement(i)]) != 0)
    {                           /* Allocation of new memory for a failed Remove all elements added so far and return 1 */
      for (j = b.N; j != i; --j)
      {
        ComplexWaveFunctionAddElementConserveBasis(*a, &b.Determinant[BlokElementDet(i, b.NWord, 0)],
                                                   -alphaRe * b.ValueRe[BlokElement(i)],
                                                   -alphaRe * b.ValueIm[BlokElement(i)]);
      }
      return 1;
    }
  }
  return 0;
}

uint32_t ComplexWaveFunctionAddComplexMultipliedUnionBasis(WaveFunctionType * const a, const WaveFunctionType b,
                                                           const double alphaRe, const double alphaIm)
{
  uint32_t i, j;
  for (i = b.N; i != 0; --i)
  {
    if (ComplexWaveFunctionAddElement
        (a, &b.Determinant[BlokElementDet(i, b.NWord, 0)],
         alphaRe * b.ValueRe[BlokElement(i)] - alphaIm * b.ValueIm[BlokElement(i)],
         alphaRe * b.ValueIm[BlokElement(i)] + alphaIm * b.ValueRe[BlokElement(i)]) != 0)
    {                           /* Allocation of new memory for a failed Remove all elements added so far and return 1 */
      for (j = b.N; j != i; --j)
      {
        ComplexWaveFunctionAddElementConserveBasis(*a, &b.Determinant[BlokElementDet(i, b.NWord, 0)],
                                                   -alphaRe * b.ValueRe[BlokElement(i)] +
                                                   alphaIm * b.ValueIm[BlokElement(i)],
                                                   -alphaRe * b.ValueIm[BlokElement(i)] -
                                                   alphaIm * b.ValueRe[BlokElement(i)]);
      }
      return 1;
    }
  }
  return 0;
}

uint32_t RealWaveFunctionAddComplexMultipliedUnionBasis(WaveFunctionType * const a, const WaveFunctionType b,
                                                           const double alphaRe, const double alphaIm)
{
  uint32_t i, j;
  for (i = b.N; i != 0; --i)
  {
    if (ComplexWaveFunctionAddElement
        (a, &b.Determinant[BlokElementDet(i, b.NWord, 0)],
         alphaRe * b.ValueRe[BlokElement(i)],alphaIm * b.ValueRe[BlokElement(i)]) != 0)
    {                           /* Allocation of new memory for a failed Remove all elements added so far and return 1 */
      for (j = b.N; j != i; --j)
      {
        ComplexWaveFunctionAddElementConserveBasis(*a, &b.Determinant[BlokElementDet(i, b.NWord, 0)],
                                                   -alphaRe * b.ValueRe[BlokElement(i)] , -alphaIm * b.ValueRe[BlokElement(i)]);
      }
      return 1;
    }
  }
  return 0;
}

/*
 *   (a) <- a b
 *   This function leaves b in the memory and adds values to a
 *   if b_i is not in the basis of a it will be ignored
 *   Therefor values of a might change and a is returned from the function
 *   If a is an empty function b is copied into a (use of copy is strongly recommanded)
 */
void ComplexWaveFunctionAddConserveBasis(WaveFunctionType * const a, const WaveFunctionType b)
{
  uint32_t i;
  int    iOMP;
#pragma omp parallel for private(i)
  for (iOMP = 1; iOMP <= (int)b.N; iOMP++)
  {
    i = (uint32_t) iOMP;
    ComplexWaveFunctionAddElementConserveBasis(*a, &b.Determinant[BlokElementDet(i, b.NWord, 0)],
                                               b.ValueRe[BlokElement(i)], b.ValueIm[BlokElement(i)]);
  }
}

/*
*   (a) <- a+alpha b
*   This function leaves b in the memory and adds values to a
*   if b_i is not in the basis of a it will be ignored
*   Therefor values of a might change and a is returned from the function
*   If a is an empty function b is copied into a (use of copy is strongly recommanded)
*/
void ComplexWaveFunctionAddMultipliedConserveBasis(WaveFunctionType * const a, const WaveFunctionType b,
                                                   const double alphaRe)
{
  uint32_t i;
  int    iOMP;
#pragma omp parallel for private(i)
  for (iOMP = 1; iOMP <= (int)b.N; iOMP++)
  {
    i = (uint32_t) iOMP;
    ComplexWaveFunctionAddElementConserveBasis(*a, &b.Determinant[BlokElementDet(i, b.NWord, 0)], alphaRe * b.ValueRe[BlokElement(i)], alphaRe * b.ValueIm[BlokElement(i)]);
  }
}

void ComplexWaveFunctionAddComplexMultipliedConserveBasis(WaveFunctionType * const a, const WaveFunctionType b, const double alphaRe, const double alphaIm)
{
  uint32_t i;
  int    iOMP;
#pragma omp parallel for private(i)
  for (iOMP = 1; iOMP <= (int)b.N; iOMP++)
  {
    i = (uint32_t) iOMP;
    ComplexWaveFunctionAddElementConserveBasis(*a, &b.Determinant[BlokElementDet(i, b.NWord, 0)], alphaRe * b.ValueRe[BlokElement(i)] - alphaIm * b.ValueIm[BlokElement(i)], alphaRe * b.ValueIm[BlokElement(i)] + alphaIm * b.ValueRe[BlokElement(i)]);
  }
}

/*
*   a->a/|a|
*   Look out this function does not alter the values stored in a
*   As a only stores pointers
*   This function does alter the values a points to
*/
void ComplexWaveFunctionNormalize(WaveFunctionType * const a)
{
  uint32_t iBlok, iElement;
  int    iOMP;
  double Norm = 0, NormA = 0, NormB = 0;
  a->ValueRe[0][0] = 0;
  a->ValueIm[0][0] = 0;
  iBlok = a->NBlok - 1;
#pragma omp parallel for private(iElement) reduction(+ : NormA)
  for (iOMP = 0; iOMP <= (int)(a->N - iBlok * BlokSize); iOMP++)
  {
    iElement = (uint32_t) iOMP;
    NormA +=
      a->ValueRe[iBlok][iElement] * a->ValueRe[iBlok][iElement] +
      a->ValueIm[iBlok][iElement] * a->ValueIm[iBlok][iElement];
  }
#pragma omp parallel for private(iElement,iBlok) reduction(+ : NormB)
  for (iOMP = 0; iOMP < (int)(a->NBlok - 1); iOMP++)
  {
    iBlok = (uint32_t) iOMP;
    for (iElement = 0; iElement < BlokSize; ++iElement)
    {
      NormB +=
        a->ValueRe[iBlok][iElement] * a->ValueRe[iBlok][iElement] +
        a->ValueIm[iBlok][iElement] * a->ValueIm[iBlok][iElement];
    }
  }
  Norm = sqrt(NormA + NormB);
  iBlok = a->NBlok - 1;
#pragma omp parallel for private(iElement)
  for (iOMP = 0; iOMP <= (int)(a->N - iBlok * BlokSize); iOMP++)
  {
    iElement = (uint32_t) iOMP;
    a->ValueRe[iBlok][iElement] /= Norm;
    a->ValueIm[iBlok][iElement] /= Norm;
  }
#pragma omp parallel for private(iElement,iBlok)
  for (iOMP = 0; iOMP < (int)(a->NBlok - 1); iOMP++)
  {
    iBlok = (uint32_t) iOMP;
    for (iElement = 0; iElement < BlokSize; ++iElement)
    {
      a->ValueRe[iBlok][iElement] /= Norm;
      a->ValueIm[iBlok][iElement] /= Norm;
    }
  }
}

/*
*   a->a*b
*   Look out this function does not alter the values stored in a
*   As a only stores pointers
*   This function does alter the values a points to
*/
void ComplexWaveFunctionMultiply(WaveFunctionType a, const double b)
{
  uint32_t iBlok, iElement;
  int    iOMP;
  iBlok = a.NBlok - 1;
#pragma omp parallel for private(iElement)
  for (iOMP = 0; iOMP <= (int)(a.N - iBlok * BlokSize); ++iOMP)
  {
    iElement = (uint32_t) iOMP;
    a.ValueRe[iBlok][iElement] *= b;
    a.ValueIm[iBlok][iElement] *= b;
  }
#pragma omp parallel for private(iElement,iBlok)
  for (iOMP = 0; iOMP < (int)(a.NBlok - 1); iOMP++)
  {
    iBlok = (uint32_t) iOMP;
    for (iElement = 0; iElement < BlokSize; ++iElement)
    {
      a.ValueRe[iBlok][iElement] *= b;
      a.ValueIm[iBlok][iElement] *= b;
    }
  }
}


/*
 *   a->a^*
 *   Look out this function does not alter the values stored in a
 *   As a only stores pointers
 *   This function does alter the values a points to
 */
void ComplexWaveFunctionConjugate(WaveFunctionType a)
{
  uint32_t iBlok, iElement;
  int    iOMP;
  iBlok = a.NBlok - 1;
#pragma omp parallel for private(iElement)
  for (iOMP = 0; iOMP <= (int)(a.N - iBlok * BlokSize); ++iOMP)
  {
    iElement = (uint32_t) iOMP;
    a.ValueIm[iBlok][iElement] *= -1;
  }
#pragma omp parallel for private(iElement,iBlok)
  for (iOMP = 0; iOMP < (int)(a.NBlok - 1); iOMP++)
  {
    iBlok = (uint32_t) iOMP;
    for (iElement = 0; iElement < BlokSize; ++iElement)
    {
      a.ValueIm[iBlok][iElement] *= -1;
    }
  }
}


/*
*   a->a*b
*   Look out this function does not alter the values stored in a
*   As a only stores pointers
*   This function does alter the values a points to
*/
void ComplexWaveFunctionComplexMultiply(WaveFunctionType * a, const double bRe, const double bIm)
{
  uint32_t iBlok, iElement;
  double temp;
  int    iOMP;
  iBlok = a->NBlok - 1;
#pragma omp parallel for private(iElement,temp)
  for (iOMP = 0; iOMP <= (int)(a->N - iBlok * BlokSize); ++iOMP)
  {
    iElement = (uint32_t) iOMP;
    temp = bRe * a->ValueRe[iBlok][iElement] - bIm * a->ValueIm[iBlok][iElement];
    a->ValueIm[iBlok][iElement] = bIm * a->ValueRe[iBlok][iElement] + bRe * a->ValueIm[iBlok][iElement];
    a->ValueRe[iBlok][iElement] = temp;
  }
#pragma omp parallel for private(iElement,iBlok,temp)
  for (iOMP = 0; iOMP < (int)(a->NBlok - 1); iOMP++)
  {
    iBlok = (uint32_t) iOMP;
    for (iElement = 0; iElement < BlokSize; ++iElement)
    {
      temp = bRe * a->ValueRe[iBlok][iElement] - bIm * a->ValueIm[iBlok][iElement];
      a->ValueIm[iBlok][iElement] = bIm * a->ValueRe[iBlok][iElement] + bRe * a->ValueIm[iBlok][iElement];
      a->ValueRe[iBlok][iElement] = temp;
    }
  }
}

/*
*   <-|a|
*/
double ComplexWaveFunctionCalculateNorm(const WaveFunctionType a)
{
  uint32_t iBlok, iElement;
  int    iOMP;
  double NormA = 0, NormB = 0, Norm = 0;
  iBlok = a.NBlok - 1;
#pragma omp parallel for private(iElement) reduction(+ : NormA)
  for (iOMP = 0; iOMP <= (int)(a.N - iBlok * BlokSize); ++iOMP)
  {
    iElement = (uint32_t) iOMP;
    NormA += a.ValueRe[iBlok][iElement] * a.ValueRe[iBlok][iElement] + a.ValueIm[iBlok][iElement] * a.ValueIm[iBlok][iElement];
  }
#pragma omp parallel for private(iElement,iBlok) reduction(+ : NormB)
  for (iOMP = 0; iOMP < (int)(a.NBlok - 1); iOMP++)
  {
    iBlok = (uint32_t) iOMP;
    for (iElement = 0; iElement < BlokSize; ++iElement)
    {
      NormB += a.ValueRe[iBlok][iElement] * a.ValueRe[iBlok][iElement] + a.ValueIm[iBlok][iElement] * a.ValueIm[iBlok][iElement];
    }
  }
  Norm = sqrt(NormA + NormB);
  return Norm;
}

/*
*   calculate <a|b>
*   The first entry is complex conjugated
*/
void ComplexWaveFunctionDotProductSameBasis(const WaveFunctionType a, const WaveFunctionType b, double *DotRe, double *DotIm)
{
  uint32_t iBlok, iElement;
  int    iOMP;
  double DotReOMPA = 0, DotImOMPA = 0, DotReOMPB = 0, DotImOMPB = 0;
  iBlok = a.NBlok - 1;
#pragma omp parallel for private(iElement) reduction(+ : DotReOMPA,DotImOMPA)
  for (iOMP = 0; iOMP <= (int)(a.N - iBlok * BlokSize); ++iOMP)
  {
    iElement = (uint32_t) iOMP;
    DotReOMPA += a.ValueRe[iBlok][iElement] * b.ValueRe[iBlok][iElement] + a.ValueIm[iBlok][iElement] * b.ValueIm[iBlok][iElement];
    DotImOMPA += a.ValueRe[iBlok][iElement] * b.ValueIm[iBlok][iElement] - a.ValueIm[iBlok][iElement] * b.ValueRe[iBlok][iElement];
  }
#pragma omp parallel for private(iElement,iBlok) reduction(+ : DotReOMPB, DotImOMPB)
  for (iOMP = 0; iOMP < (int)(a.NBlok - 1); iOMP++)
  {
    iBlok = (uint32_t) iOMP;
    for (iElement = 0; iElement < BlokSize; ++iElement)
    {
      DotReOMPB += a.ValueRe[iBlok][iElement] * b.ValueRe[iBlok][iElement] + a.ValueIm[iBlok][iElement] * b.ValueIm[iBlok][iElement];
      DotImOMPB += a.ValueRe[iBlok][iElement] * b.ValueIm[iBlok][iElement] - a.ValueIm[iBlok][iElement] * b.ValueRe[iBlok][iElement];
    }
  }
  *DotRe = DotReOMPA + DotReOMPB;
  *DotIm = DotImOMPA + DotImOMPB;
}

/*
*   calculate <a|b>
*/
void ComplexWaveFunctionDotProduct(const WaveFunctionType *a, const WaveFunctionType *b, double *DotRe, double *DotIm)
{
  uint32_t i, j, NextIndex, PsiAIsPsiB;
  int    iOMP;
  double DotReOMP = 0, DotImOMP = 0;  
  int exchanged = 0;

  /* switch coefficients to reduce number of hash lookups */
  if (a->N > b->N)
  {
    const WaveFunctionType *temp = a;
    a = b;
    b = temp;
    exchanged = 1;
  }


#pragma omp parallel for private(i,j,NextIndex,PsiAIsPsiB) reduction(+ : DotReOMP,DotImOMP)
  for (iOMP = 1; iOMP <= (int)(a->N); ++iOMP)
  {
    i = (uint32_t) iOMP;
    NextIndex = b->KeyList[hash(&a->Determinant[BlokElementDet(i, a->NWord, 0)], b->NWord, b->NBitsKey)];
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
        DotReOMP +=
          a->ValueRe[BlokElement(i)] * b->ValueRe[BlokElement(NextIndex)] +
          a->ValueIm[BlokElement(i)] * b->ValueIm[BlokElement(NextIndex)];
        DotImOMP +=
          a->ValueRe[BlokElement(i)] * b->ValueIm[BlokElement(NextIndex)] -
          a->ValueIm[BlokElement(i)] * b->ValueRe[BlokElement(NextIndex)];
        NextIndex = 0;
      }
      else
      {
        NextIndex = b->Next[BlokElement(NextIndex)];
      }
    }
  }
  *DotRe = DotReOMP;
  *DotIm = DotImOMP*(exchanged?-1:1);
}

/* Set the pre-factors of a to a random value
*/
void ComplexWaveFunctionRandomize(WaveFunctionType * const a)
{
  uint32_t iBlok, iElement;
  int    iOMP;
  iBlok = a->NBlok - 1;
#pragma omp parallel for private(iElement)
  for (iOMP = 0; iOMP <= (int)(a->N - iBlok * BlokSize); ++iOMP)
  {
    iElement = (uint32_t) iOMP;
    a->ValueRe[iBlok][iElement] = rand() - RAND_MAX / 2;
    a->ValueIm[iBlok][iElement] = rand() - RAND_MAX / 2;
  }
#pragma omp parallel for private(iElement,iBlok)
  for (iOMP = 0; iOMP < (int)(a->NBlok - 1); iOMP++)
  {
    iBlok = (uint32_t) iOMP;
    for (iElement = 0; iElement < BlokSize; ++iElement)
    {
      a->ValueRe[iBlok][iElement] = rand() - RAND_MAX / 2;
      a->ValueIm[iBlok][iElement] = rand() - RAND_MAX / 2;
    }
  }
  ComplexWaveFunctionNormalize(a);
}

/*
*   remove determinants from the basis that have a prefactor smaller then epsilon  
*/
void ComplexWaveFunctionTruncateBasis(WaveFunctionType * const a, const double epsilon)
{
  const double Sqrepsilon = Sqr(epsilon);
  uint32_t startN, iBlok, iElement;
  a->ValueRe[0][0] = 0;
  a->ValueIm[0][0] = 0;
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
      if (Sqr(a->ValueRe[iBlok][iElement]) + Sqr(a->ValueIm[iBlok][iElement]) > Sqrepsilon)
        if (ComplexWaveFunctionAddElement(a, &a->Determinant[iBlok][iElement * a->NWord], a->ValueRe[iBlok][iElement], a->ValueIm[iBlok][iElement]) != 0) AllocFailed();
    }
  }
  iBlok = a->NBlok - 1;
  for (iElement = 0; iElement < 1 + startN - iBlok * BlokSize; ++iElement)
  {
    a->Next[iBlok][iElement] = 0;
    if (Sqr(a->ValueRe[iBlok][iElement]) + Sqr(a->ValueIm[iBlok][iElement]) > Sqrepsilon)
      if (ComplexWaveFunctionAddElement(a, &a->Determinant[iBlok][iElement * a->NWord], a->ValueRe[iBlok][iElement], a->ValueIm[iBlok][iElement]) != 0) AllocFailed();
  }
  for (iBlok = 1 + (a->N >> NBitsBlok); iBlok < a->NBlok; iBlok++)
  {
    free(a->ValueRe[iBlok]);
    free(a->ValueIm[iBlok]);
    free(a->Next[iBlok]);
    free(a->Determinant[iBlok]);
  }
  a->NBlok = 1 + (a->N >> NBitsBlok);
  a->ValueRe = realloc(a->ValueRe, a->NBlok * sizeof(double *));
  a->ValueIm = realloc(a->ValueIm, a->NBlok * sizeof(double *));
  a->Next = realloc(a->Next, a->NBlok * sizeof(uint32_t *));
  a->Determinant = realloc(a->Determinant, a->NBlok * sizeof(uint8_t *));
}

/*
 *  change the size of the keylist
 */
uint32_t ComplexWaveFunctionResetKeylist(WaveFunctionType * const a, const uint32_t NBitsKey)
{
  const double Sqrepsilon = 10*DBL_MIN;
  const uint32_t startN=a->N;
  uint32_t iBlok, iElement;
  uint32_t *tempKeyList;
  a->ValueRe[0][0] = 0;
  a->ValueIm[0][0] = 0;
  tempKeyList=a->KeyList;
  if((a->KeyList=realloc(a->KeyList, sizeof(uint32_t)*(1<<NBitsKey)))==NULL)
  {
    a->KeyList=tempKeyList;
    printf("realloc failed in RealWaveFunctionResetKeylist\n");
    return 1;
  }
  a->NBitsKey=NBitsKey;
  a->N = 0;
  for (iBlok = 0; iBlok < a->NBlok - 1; iBlok++)
  {
    for (iElement = 0; iElement < BlokSize; iElement++)
    {
      a->Next[iBlok][iElement] = 0;
      if (Sqr(a->ValueRe[iBlok][iElement]) + Sqr(a->ValueIm[iBlok][iElement]) > Sqrepsilon)
        if (ComplexWaveFunctionAddElement(a, &a->Determinant[iBlok][iElement * a->NWord], a->ValueRe[iBlok][iElement], a->ValueIm[iBlok][iElement]) != 0) AllocFailed();
    }
  }
  iBlok = a->NBlok - 1;
  for (iElement = 0; iElement < 1 + startN - iBlok * BlokSize; ++iElement)
  {
    a->Next[iBlok][iElement] = 0;
    if (Sqr(a->ValueRe[iBlok][iElement]) + Sqr(a->ValueIm[iBlok][iElement]) > Sqrepsilon)
      if (ComplexWaveFunctionAddElement(a, &a->Determinant[iBlok][iElement * a->NWord], a->ValueRe[iBlok][iElement], a->ValueIm[iBlok][iElement]) != 0) AllocFailed();
  }
  for (iBlok = 1 + (a->N >> NBitsBlok); iBlok < a->NBlok; iBlok++)
  {
    free(a->ValueRe[iBlok]);
    free(a->ValueIm[iBlok]);
    free(a->Next[iBlok]);
    free(a->Determinant[iBlok]);
  }
  a->NBlok = 1 + (a->N >> NBitsBlok);
  a->ValueRe = realloc(a->ValueRe, a->NBlok * sizeof(double *));
  a->ValueIm = realloc(a->ValueIm, a->NBlok * sizeof(double *));
  a->Next = realloc(a->Next, a->NBlok * sizeof(uint32_t *));
  a->Determinant = realloc(a->Determinant, a->NBlok * sizeof(uint8_t *));
  return 0;
}

/*
 *   remove determinants from the basis that do not fullfill the restrictions
 */
void ComplexWaveFunctionRestrict(WaveFunctionType * const a, const DeterminantRestrictionType res)
{
  uint32_t startN, iBlok, iElement;
  a->ValueRe[0][0] = 0;
  a->ValueIm[0][0] = 0;
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
      if (QDetNotRestricted(&(a->Determinant[iBlok][iElement*a->NWord]), a->NFermion, a->NBoson, a->NWord, res))
        if (ComplexWaveFunctionAddElement(a, &a->Determinant[iBlok][iElement * a->NWord], a->ValueRe[iBlok][iElement], a->ValueIm[iBlok][iElement]) != 0) AllocFailed();
    }
  }
  iBlok = a->NBlok - 1;
  for (iElement = 0; iElement < 1 + startN - iBlok * BlokSize; ++iElement)
  {
    a->Next[iBlok][iElement] = 0;
    if (QDetNotRestricted(&(a->Determinant[iBlok][iElement*a->NWord]), a->NFermion, a->NBoson, a->NWord, res))
      if (ComplexWaveFunctionAddElement(a, &a->Determinant[iBlok][iElement * a->NWord], a->ValueRe[iBlok][iElement], a->ValueIm[iBlok][iElement]) != 0) AllocFailed();
  }
  for (iBlok = 1 + (a->N >> NBitsBlok); iBlok < a->NBlok; iBlok++)
  {
    free(a->ValueRe[iBlok]);
    free(a->ValueIm[iBlok]);
    free(a->Next[iBlok]);
    free(a->Determinant[iBlok]);
  }
  a->NBlok = 1 + (a->N >> NBitsBlok);
  a->ValueRe = realloc(a->ValueRe, a->NBlok * sizeof(double *));
  a->ValueIm = realloc(a->ValueIm, a->NBlok * sizeof(double *));
  a->Next = realloc(a->Next, a->NBlok * sizeof(uint32_t *));
  a->Determinant = realloc(a->Determinant, a->NBlok * sizeof(uint8_t *));
}

/*  remove determinants from the basis that have a prefactor smaller then epsilon
 *  if all imaginary parts are in absolute value smaller than epsilon make complex function real
 */
void ComplexWaveFunctionTruncateBasisAllowReal(WaveFunctionType * const a, const double epsilon)
{
  const double Sqrepsilon = Sqr(epsilon);
  Bool keepcomplex=False;
  uint32_t startN, iBlok, iElement;
  a->ValueRe[0][0] = 0;
  a->ValueIm[0][0] = 0;
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
      if (Sqr(a->ValueRe[iBlok][iElement]) + Sqr(a->ValueIm[iBlok][iElement]) > Sqrepsilon)
      {
        if (ComplexWaveFunctionAddElement(a, &a->Determinant[iBlok][iElement * a->NWord], a->ValueRe[iBlok][iElement], a->ValueIm[iBlok][iElement]) != 0) AllocFailed();
        if (Sqr(a->ValueIm[iBlok][iElement]) > Sqrepsilon) keepcomplex = True;
      }
    }
  }
  iBlok = a->NBlok - 1;
  for (iElement = 0; iElement < 1 + startN - iBlok * BlokSize; ++iElement)
  {
    a->Next[iBlok][iElement] = 0;
    if (Sqr(a->ValueRe[iBlok][iElement]) + Sqr(a->ValueIm[iBlok][iElement]) > Sqrepsilon)
    {
      if (ComplexWaveFunctionAddElement(a, &a->Determinant[iBlok][iElement * a->NWord], a->ValueRe[iBlok][iElement], a->ValueIm[iBlok][iElement]) != 0) AllocFailed();
      if (Sqr(a->ValueIm[iBlok][iElement]) > Sqrepsilon) keepcomplex = True;
    }
  }
  for (iBlok = 1 + (a->N >> NBitsBlok); iBlok < a->NBlok; iBlok++)
  {
    free(a->ValueRe[iBlok]);
    free(a->ValueIm[iBlok]);
    free(a->Next[iBlok]);
    free(a->Determinant[iBlok]);
  }
  a->NBlok = 1 + (a->N >> NBitsBlok);
  a->ValueRe = realloc(a->ValueRe, a->NBlok * sizeof(double *));
  a->Next = realloc(a->Next, a->NBlok * sizeof(uint32_t *));
  a->Determinant = realloc(a->Determinant, a->NBlok * sizeof(uint8_t *));
  if(keepcomplex)
  {
    a->ValueIm = realloc(a->ValueIm, a->NBlok * sizeof(double *));
  }
  else
  {
    for (iBlok = 0; iBlok < a->NBlok; iBlok++)
    {
      free(a->ValueIm[iBlok]);
    }
    free(a->ValueIm);
    a->QComplex = 0;
  }
}


/*
 *   Make psi orthogonal to OtherPsi, psi and other psi are on the same basis
 */
void ComplexWaveFunctionOrthogonalizeSameBasis(WaveFunctionType * const psi, const WaveFunctionType * const OtherPsi, const uint32_t Npsi)
{
  uint32_t i;
  double aRe, aIm;
  for (i = 0; i < Npsi; i++)
  {
    ComplexWaveFunctionDotProductSameBasis(OtherPsi[i], *psi, &aRe, &aIm);
    ComplexWaveFunctionAddComplexMultipliedSameBasis(*psi, OtherPsi[i], -aRe, -aIm);
  }
}

/*
 *   Make psi orthogonal to OtherPsi, but conserve the basis set of psi
 */
void ComplexWaveFunctionOrthogonalizeConserveBasis(WaveFunctionType * const psi, const WaveFunctionType * const OtherPsi, const uint32_t Npsi)
{
  uint32_t i;
  double aRe, aIm;
  for (i = 0; i < Npsi; i++)
  {
    ComplexWaveFunctionDotProduct(&OtherPsi[i], psi, &aRe, &aIm);
    ComplexWaveFunctionAddComplexMultipliedConserveBasis(psi, OtherPsi[i], -aRe, -aIm);
  }
}

/*  
*   Make psi orthogonal to OtherPsi, if needed enlarge the basis set of psi
*/
uint32_t ComplexWaveFunctionOrthogonalize(WaveFunctionType * const psi, const WaveFunctionType * const OtherPsi,
                                          const uint32_t Npsi)
{
  uint32_t i;
  double aRe, aIm;
  for (i = 0; i < Npsi; i++)
  {
    ComplexWaveFunctionDotProduct(&OtherPsi[i], psi, &aRe, &aIm);
    if (ComplexWaveFunctionAddComplexMultipliedUnionBasis(psi, OtherPsi[i], -aRe, -aIm))
    {
      return 1;
    }
  }
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
uint32_t ComplexMappedDetToWaveFunction(WaveFunctionType * const A, uint8_t *Det, MatrixType const U)
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
  
  ComplexWaveFunctionAddElement(A, DetZero, 1.0, 0.0);
  
  for(i=0;i<A->NFermion;i++)
  {
    oadd[0]=CrAnMask | (uint16_t)i;
    ComplexOperatorAddLader(&O, oadd, 1, 0.0);
  }
  
  for(i=0;i<A->NFermion;i++)
  {
    /* test if Fermion i is occupied */
    if(Det[i >> 3] & (128 >> (i & 7)))
    {
      for(j=0;j<A->NFermion;j++)
      {
        O.l[1].V[j].Re=U.VRe[j][i];
        O.l[1].V[j].Im=-U.VIm[j][i];
      }
      OperatorPsiCC(O, *A, &Temp);
      ComplexWaveFunctionTruncateBasis(&Temp, 10*DBL_EPSILON);
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
 *  U is a (NFermionA+NBosonA)*(NFermionB+NBosonB) demensional matrix that defines the transformation
 *  For operators we have
 *  a_i        = \sum_j U_{i,j}   a_j
 *  a_i^{\dag} = \sum_j U_{i,j}^* a_j^{dag}
 *  writing:
 *    psi = \sum_j alpha_j a_j^{dag} |0>
 *    we have
 *    U psi = \sum_j U_{i,j}^* alpha_j a_j^{dag} |0>
 *
 * assume A B and U are complex (version CCC)
 */
uint32_t MappedWaveFunctionAddCCC(WaveFunctionType * const A, WaveFunctionType const B, MatrixType const U)
{
  uint8_t *Det;
  uint32_t idetB;
  WaveFunctionType Udet;
  Udet.NFermion=A->NFermion;
  Udet.NBoson=A->NBoson;
  Udet.QComplex=1;
  Udet.NBitsKey=A->NBitsKey;
  WaveFunctionInit(&Udet);
  /* loop over all determinants in B */
  for(idetB=1;idetB<=B.N;idetB++)
  {
    Det=&(B.Determinant[BlokElementDet(idetB, B.NWord, 0)]);
    ComplexWaveFunctionReset(&Udet);
    ComplexMappedDetToWaveFunction(&Udet, Det, U);
    ComplexWaveFunctionAddComplexMultipliedUnionBasis(A, Udet, B.ValueRe[BlokElement(idetB)], B.ValueIm[BlokElement(idetB)]);
  }
  ComplexWaveFunctionTruncateBasis(A, 10*DBL_EPSILON);
  ComplexWaveFunctionNormalize(A);
  WaveFunctionFree(&Udet);
  return 0;
}
uint32_t MappedWaveFunctionAddCCR(WaveFunctionType * const A, WaveFunctionType const B, MatrixType const U)
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
    RealWaveFunctionAddComplexMultipliedUnionBasis(A, Udet, B.ValueRe[BlokElement(idetB)], B.ValueIm[BlokElement(idetB)]);
  }
  ComplexWaveFunctionTruncateBasis(A, 10*DBL_EPSILON);
  ComplexWaveFunctionNormalize(A);
  WaveFunctionFree(&Udet);
  return 0;
}
uint32_t MappedWaveFunctionAddCRC(WaveFunctionType * const A, WaveFunctionType const B, MatrixType const U)
{
  uint8_t *Det;
  uint32_t idetB;
  WaveFunctionType Udet;
  Udet.NFermion=A->NFermion;
  Udet.NBoson=A->NBoson;
  Udet.QComplex=1;
  Udet.NBitsKey=A->NBitsKey;
  WaveFunctionInit(&Udet);
  /* loop over all determinants in B */
  for(idetB=1;idetB<=B.N;idetB++)
  {
    Det=&(B.Determinant[BlokElementDet(idetB, B.NWord, 0)]);
    ComplexWaveFunctionReset(&Udet);
    ComplexMappedDetToWaveFunction(&Udet, Det, U);
    ComplexWaveFunctionAddMultipliedUnionBasis(A, Udet, B.ValueRe[BlokElement(idetB)]);
  }
  ComplexWaveFunctionTruncateBasis(A, 10*DBL_EPSILON);
  ComplexWaveFunctionNormalize(A);
  WaveFunctionFree(&Udet);
  return 0;
}
