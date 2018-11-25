/*  Block Wavefunctions are arrays of wave-function types wich have the same basis.
 *  for the array psi[] the values
 *  psi[i].QComplex
 *        .N
 *        .NBlok
 *        .NFermion
 *        .NBoson
 *        .Nword
 *        .Next
 *        .Determinant
 *        .NBitskey
 *        .Keylist
 *  should all be the same for all [i]. It is the programers responsibility to make sure this is the case. Otherwise these functions will produce junk
 *  psi[i].Name
 *        .ValueRe
 *        .ValueIm
 *  define the Npsi functions
 */

/*  Free the wave-functions psi[i] from memory
 *
 */
void BlockWaveFunctionFree(WaveFunctionType *psi, const uint32_t Npsi)
{
  uint32_t i;
  for(i=1;i<Npsi;i++)
  {
    WaveFunctionFreeValues(&(psi[i]));
  }
  if(Npsi>0) WaveFunctionFree(&(psi[0]));
}

uint32_t MakeComplexBlockWaveFunctionFromReal(WaveFunctionType *psi, const uint32_t Npsi)
{
  uint32_t i,j,k;
  for(i=0;i<Npsi;i++)
  {
    psi[i].ValueIm = (double**)malloc(psi[0].NBlok*sizeof(double*));
    if(psi[i].ValueIm==NULL)
    {
      for(j=0;j<i;j++)
      {
        free(psi[j].ValueIm);
      }
      printf("malloc failed in MakeComplexBlockWaveFunctionFromReal\n");
      return 1;
    }
  }
  for(i=0;i<Npsi;i++)
  {
    for(j=0;j<psi[0].NBlok;j++)
    {
      psi[i].ValueIm[j] = (double*)calloc(BlokSize,sizeof(double));
      if(psi[i].ValueIm[j]==NULL)
      {
        for(k=0;k<j;k++)
        {
          free(psi[i].ValueIm[k]);
        }
        for(k=0;k<i;k++)
        {
          for(j=0;j<psi[0].NBlok;j++)
          {
            free(psi[k].ValueIm[j]);
          }
          free(psi[k].ValueIm);
        }
        printf("malloc failed in MakeComplexBlockWaveFunctionFromReal\n");
        return 1;
      }
    }
  }
  for(i=0;i<Npsi;i++)
  {
    psi[i].QComplex=1;
  }
  return 0;
}


/*
 *   add ellement with single slater determinant det and Npsi pre-factors a[i] to Npsi wave functions psi[i]
 *   all psi[i] must have the same basis (all values are the same execept for ValueRe
 *   returns 0 if succeded
 *   returns 1 if no memory available in order to add det to psi
 *   If succeded psi -> psi+det
 *   If not succeded psi -> psi
 */
uint32_t RealBlockWaveFunctionAddElement(WaveFunctionType * const psi, const uint8_t * const det, const double * const a, const uint32_t Npsi)
{
  uint32_t i, j, k, NextIndex, ThereIsANext, PsiAIsPsiB;
  void  *temp;
  if(Npsi==0) return 0;
  /* test if the key that belongs to psi already excists */
  if (psi[0].KeyList[hash(det, psi[0].NWord, psi[0].NBitsKey)] == 0)
  {                             /* add a new element to psi with unique new key */
    for(j=0;j<Npsi;j++) psi[j].N++;
    /* test if one has to add a new block to psi in order to store the new element of psi */
    if (psi[0].N >= BlokSize * psi[0].NBlok)
    {                           /* add a new block to psi */
      temp = (uint32_t **)realloc(psi[0].Next, (psi[0].NBlok+1) * sizeof(uint32_t *));
      if(temp==NULL)
      {                         /* Allocation failed */
        printf("alloc failed RealBlockWaveFunctionAddElement Next\n");
        fflush(stdout);
        /* restore psi and return 0 */
        for(j=0;j<Npsi;j++) psi[j].N--;
        return 1;
      }
      for(j=0;j<Npsi;j++) psi[j].Next = (uint32_t**)temp;
      temp = (uint32_t **)realloc(psi[0].Determinant, (psi[0].NBlok+1) * sizeof(uint32_t *));
      if(temp==NULL)
      {                         /* Allocation failed */
        printf("alloc failed RealBlockWaveFunctionAddElement Determinant\n");
        fflush(stdout);
        /* restore psi and return 0 */
        for(j=0;j<Npsi;j++) psi[j].N--;
        return 1;
      }
      for(j=0;j<Npsi;j++) psi[j].Determinant = (uint8_t**)temp;
      for(j=0;j<Npsi;j++)
      {
        temp = (double **)realloc(psi[j].ValueRe, (psi[j].NBlok+1) * sizeof(double *));
        if (temp == NULL)
        {                         /* Allocation failed */
          printf("alloc failed RealBlockWaveFunctionAddElement ValueRe %i\n",j);
          fflush(stdout);
          /* restore psi and return 0 */
          for(j=0;j<Npsi;j++) psi[j].N--;
          return (1);
        }
        psi[j].ValueRe = (double**)temp;
      }
      psi[0].Next[psi[0].NBlok] = (uint32_t*)malloc(BlokSize * sizeof(uint32_t));
      if (psi[0].Next[psi[0].NBlok] == NULL)
      {                         /* Allocation failed */
        printf("alloc failed RealBlockWaveFunctionAddElement Next[-1]\n");
        fflush(stdout);
        /* restore psi and return 0 */
        for(j=0;j<Npsi;j++) psi[j].N--;
        return (1);
      }
      psi[0].Determinant[psi[0].NBlok] = (uint8_t*)malloc(BlokSize * psi[0].NWord * sizeof(uint8_t));
      if (psi[0].Determinant[psi[0].NBlok] == NULL)
      {                         /* Allocation failed */
        printf("alloc failed RealBlockWaveFunctionAddElement Determinant[-1]\n");
        fflush(stdout);
        /* restore psi and return 0 */
        for(j=0;j<Npsi;j++) psi[j].N--;
        free(psi[0].Next[psi[0].NBlok]);
        return (1);
      }
      for(j=0;j<Npsi;j++)
      {
        psi[j].ValueRe[psi[0].NBlok] = (double*)malloc(BlokSize * sizeof(double));
        if(psi[j].ValueRe[psi[0].NBlok]==NULL)
        {
          printf("alloc failed RealBlockWaveFunctionAddElement ValueRe[%i]\n",j);
          fflush(stdout);
          /* restore psi and return 0 */
          for(k=0;k<Npsi;k++) psi[k].N--;
          free(psi[0].Next[psi[0].NBlok]);
          free(psi[0].Determinant[psi[0].NBlok]);
          for(k=0;k<j;k++) free(psi[k].ValueRe[psi[0].NBlok]);
          return (1);
        }
      }
      for(j=0;j<Npsi;j++) psi[j].NBlok++;
    }                           /* Added a new block to psi */
    for(i=0;i<Npsi;i++)
    {
      psi[i].ValueRe[BlokElement(psi[0].N)] = a[i];
    }
    psi[0].Next[BlokElement(psi[0].N)] = 0;
    for (i = 0; i < psi[0].NWord; i++)
    {
      psi[0].Determinant[BlokElementDet(psi[0].N, psi[0].NWord, i)] = det[i];
    }
    psi[0].KeyList[hash(det, psi[0].NWord, psi[0].NBitsKey)] = psi[0].N;
  }
  else
  {                             /* The key is already taken, test if this is due to a collision, or if psi is already in the basis set */
    ThereIsANext = 1;
    NextIndex = psi[0].KeyList[hash(det, psi[0].NWord, psi[0].NBitsKey)];
    while (ThereIsANext)
    {
      PsiAIsPsiB = 1;
      for (i = 0; i < psi[0].NWord; i++)
      {
        if (psi[0].Determinant[BlokElementDet(NextIndex, psi[0].NWord, i)] != det[i]) PsiAIsPsiB = 0;
      }
      if (PsiAIsPsiB)
      {                         /* psi already in basis set add and leave loop  */
        ThereIsANext = 0;
        for(i=0;i<Npsi;i++)
        {
          psi[i].ValueRe[BlokElement(NextIndex)] += a[i];
        }
      }
      else
      {                         /* psi not equal to this one, if there is a next continue, otherwise add */
        if (!psi[0].Next[BlokElement(NextIndex)])
        {                       /* add a new element to psi with already used key */
          for(j=0;j<Npsi;j++) psi[j].N++;
          /* test if one has to add a new block to psi in order to store the new element of psi */
          if (psi[0].N >= BlokSize * psi[0].NBlok)
          {                           /* add a new block to psi */
            temp = (uint32_t **)realloc(psi[0].Next, (psi[0].NBlok+1) * sizeof(uint32_t *));
            if(temp==NULL)
            {                         /* Allocation failed */
              printf("alloc failed RealBlockWaveFunctionAddElement Next\n");
              fflush(stdout);
              /* restore psi and return 0 */
              for(j=0;j<Npsi;j++) psi[j].N--;
              return 1;
            }
            for(j=0;j<Npsi;j++) psi[j].Next = (uint32_t**)temp;
            temp = (uint32_t **)realloc(psi[0].Determinant, (psi[0].NBlok+1) * sizeof(uint32_t *));
            if(temp==NULL)
            {                         /* Allocation failed */
              printf("alloc failed RealBlockWaveFunctionAddElement Determinant\n");
              fflush(stdout);
              /* restore psi and return 0 */
              for(j=0;j<Npsi;j++) psi[j].N--;
              return 1;
            }
            for(j=0;j<Npsi;j++) psi[j].Determinant = (uint8_t**)temp;
            for(j=0;j<Npsi;j++)
            {
              temp = (double **)realloc(psi[j].ValueRe, (psi[j].NBlok+1) * sizeof(double *));
              if (temp == NULL)
              {                         /* Allocation failed */
                printf("alloc failed RealBlockWaveFunctionAddElement ValueRe %i\n",j);
                fflush(stdout);
                /* restore psi and return 0 */
                for(j=0;j<Npsi;j++) psi[j].N--;
                return (1);
              }
              psi[j].ValueRe = (double**)temp;
            }
            psi[0].Next[psi[0].NBlok] = (uint32_t*)malloc(BlokSize * sizeof(uint32_t));
            if (psi[0].Next[psi[0].NBlok] == NULL)
            {                         /* Allocation failed */
              printf("alloc failed RealBlockWaveFunctionAddElement Next[-1]\n");
              fflush(stdout);
              /* restore psi and return 0 */
              for(j=0;j<Npsi;j++) psi[j].N--;
              return (1);
            }
            psi[0].Determinant[psi[0].NBlok] = (uint8_t*)malloc(BlokSize * psi[0].NWord * sizeof(uint8_t));
            if (psi[0].Determinant[psi[0].NBlok] == NULL)
            {                         /* Allocation failed */
              printf("alloc failed RealBlockWaveFunctionAddElement Determinant[-1]\n");
              fflush(stdout);
              /* restore psi and return 0 */
              for(j=0;j<Npsi;j++) psi[j].N--;
              free(psi[0].Next[psi[0].NBlok]);
              return (1);
            }
            for(j=0;j<Npsi;j++)
            {
              psi[j].ValueRe[psi[0].NBlok] = (double*)malloc(BlokSize * sizeof(double));
              if(psi[j].ValueRe[psi[0].NBlok]==NULL)
              {
                printf("alloc failed RealBlockWaveFunctionAddElement ValueRe[%i]\n",j);
                fflush(stdout);
                /* restore psi and return 0 */
                for(k=0;k<Npsi;k++) psi[k].N--;
                free(psi[0].Next[psi[0].NBlok]);
                free(psi[0].Determinant[psi[0].NBlok]);
                for(k=0;k<j;k++) free(psi[k].ValueRe[psi[0].NBlok]);
                return (1);
              }
            }
            for(j=0;j<Npsi;j++) psi[j].NBlok++;
          }          /* Added a new block to psi */
          /* add the element */
          for(i=0;i<Npsi;i++) psi[i].ValueRe[BlokElement(psi[0].N)] = a[i];
          psi[0].Next[BlokElement(psi[0].N)] = 0;
          for (i = 0; i < psi[0].NWord; i++)
          {
            psi[0].Determinant[BlokElementDet(psi[0].N, psi[0].NWord, i)] = det[i];
          }
          psi[0].Next[BlokElement(NextIndex)] = psi[0].N;
          ThereIsANext = 0;
        }
        else
        {                       /* the key is in use and there still is at least one more element with the same key, continu checking */
          NextIndex = psi[0].Next[BlokElement(NextIndex)];
        }
      }
    }
  }
  return (0);
}
uint32_t ComplexBlockWaveFunctionAddElement(WaveFunctionType * const psi, const uint8_t * const det, const double * aRe, const double * aIm, const uint32_t Npsi)
{
  uint32_t i, j, k, NextIndex, ThereIsANext, PsiAIsPsiB;
  void  *temp;
  if(Npsi==0) return 0;
  /* test if the key that belongs to psi already excists */
  if (psi[0].KeyList[hash(det, psi[0].NWord, psi[0].NBitsKey)] == 0)
  {                             /* add a new element to psi with unique new key */
    for(j=0;j<Npsi;j++) psi[j].N++;
    /* test if one has to add a new block to psi in order to store the new element of psi */
    if (psi[0].N >= BlokSize * psi[0].NBlok)
    {                           /* add a new block to psi */
      temp = (uint32_t **)realloc(psi[0].Next, (psi[0].NBlok+1) * sizeof(uint32_t *));
      if(temp==NULL)
      {                         /* Allocation failed */
        printf("alloc failed ComplexBlockWaveFunctionAddElement Next\n");
        fflush(stdout);
        /* restore psi and return 0 */
        for(j=0;j<Npsi;j++) psi[j].N--;
        return 1;
      }
      for(j=0;j<Npsi;j++) psi[j].Next = (uint32_t**)temp;
      temp = (uint32_t **)realloc(psi[0].Determinant, (psi[0].NBlok+1) * sizeof(uint32_t *));
      if(temp==NULL)
      {                         /* Allocation failed */
        printf("alloc failed ComplexBlockWaveFunctionAddElement Determinant\n");
        fflush(stdout);
        /* restore psi and return 0 */
        for(j=0;j<Npsi;j++) psi[j].N--;
        return 1;
      }
      for(j=0;j<Npsi;j++) psi[j].Determinant = (uint8_t**)temp;
      for(j=0;j<Npsi;j++)
      {
        temp = (double **)realloc(psi[j].ValueRe, (psi[j].NBlok+1) * sizeof(double *));
        if (temp == NULL)
        {                         /* Allocation failed */
          printf("alloc failed ComplexBlockWaveFunctionAddElement ValueRe %i\n",j);
          fflush(stdout);
          /* restore psi and return 0 */
          for(j=0;j<Npsi;j++) psi[j].N--;
          return (1);
        }
        psi[j].ValueRe = (double**)temp;
      }
      for(j=0;j<Npsi;j++)
      {
        temp = (double **)realloc(psi[j].ValueIm, (psi[j].NBlok+1) * sizeof(double *));
        if (temp == NULL)
        {                         /* Allocation failed */
          printf("alloc failed ComplexBlockWaveFunctionAddElement ValueIm %i\n",j);
          fflush(stdout);
          /* restore psi and return 0 */
          for(j=0;j<Npsi;j++) psi[j].N--;
          return (1);
        }
        psi[j].ValueIm = (double**)temp;
      }
      psi[0].Next[psi[0].NBlok] = (uint32_t*)malloc(BlokSize * sizeof(uint32_t));
      if (psi[0].Next[psi[0].NBlok] == NULL)
      {                         /* Allocation failed */
        printf("alloc failed ComplexBlockWaveFunctionAddElement Next[-1]\n");
        fflush(stdout);
        /* restore psi and return 0 */
        for(j=0;j<Npsi;j++) psi[j].N--;
        return (1);
      }
      psi[0].Determinant[psi[0].NBlok] = (uint8_t*)malloc(BlokSize * psi[0].NWord * sizeof(uint8_t));
      if (psi[0].Determinant[psi[0].NBlok] == NULL)
      {                         /* Allocation failed */
        printf("alloc failed ComplexBlockWaveFunctionAddElement Determinant[-1]\n");
        fflush(stdout);
        /* restore psi and return 0 */
        for(j=0;j<Npsi;j++) psi[j].N--;
        free(psi[0].Next[psi[0].NBlok]);
        return (1);
      }
      for(j=0;j<Npsi;j++)
      {
        psi[j].ValueRe[psi[0].NBlok] = (double*)malloc(BlokSize * sizeof(double));
        if(psi[j].ValueRe[psi[0].NBlok]==NULL)
        {
          printf("alloc failed ComplexBlockWaveFunctionAddElement ValueRe[%i]\n",j);
          fflush(stdout);
          /* restore psi and return 0 */
          for(k=0;k<Npsi;k++) psi[k].N--;
          free(psi[0].Next[psi[0].NBlok]);
          free(psi[0].Determinant[psi[0].NBlok]);
          for(k=0;k<j;k++) free(psi[k].ValueRe[psi[0].NBlok]);
          return (1);
        }
      }
      for(j=0;j<Npsi;j++)
      {
        psi[j].ValueIm[psi[0].NBlok] = (double*)malloc(BlokSize * sizeof(double));
        if(psi[j].ValueIm[psi[0].NBlok]==NULL)
        {
          printf("alloc failed ComplexBlockWaveFunctionAddElement ValueRe[%i]\n",j);
          fflush(stdout);
          /* restore psi and return 0 */
          for(k=0;k<Npsi;k++) psi[k].N--;
          free(psi[0].Next[psi[0].NBlok]);
          free(psi[0].Determinant[psi[0].NBlok]);
          for(k=0;k<Npsi;k++) free(psi[k].ValueRe[psi[0].NBlok]);
          for(k=0;k<j;k++) free(psi[k].ValueIm[psi[0].NBlok]);
          return (1);
        }
      }
      for(j=0;j<Npsi;j++) psi[j].NBlok++;
    }                           /* Added a new block to psi */
    for(i=0;i<Npsi;i++)
    {
      psi[i].ValueRe[BlokElement(psi[0].N)] = aRe[i];
      psi[i].ValueIm[BlokElement(psi[0].N)] = aIm[i];
    }
    psi[0].Next[BlokElement(psi[0].N)] = 0;
    for (i = 0; i < psi[0].NWord; i++)
    {
      psi[0].Determinant[BlokElementDet(psi[0].N, psi[0].NWord, i)] = det[i];
    }
    psi[0].KeyList[hash(det, psi[0].NWord, psi[0].NBitsKey)] = psi[0].N;
  }
  else
  {                             /* The key is already taken, test if this is due to a collision, or if psi is already in the basis set */
    ThereIsANext = 1;
    NextIndex = psi[0].KeyList[hash(det, psi[0].NWord, psi[0].NBitsKey)];
    while (ThereIsANext)
    {
      PsiAIsPsiB = 1;
      for (i = 0; i < psi[0].NWord; i++)
      {
        if (psi[0].Determinant[BlokElementDet(NextIndex, psi[0].NWord, i)] != det[i]) PsiAIsPsiB = 0;
      }
      if (PsiAIsPsiB)
      {                         /* psi already in basis set add and leave loop  */
        ThereIsANext = 0;
        for(i=0;i<Npsi;i++)
        {
          psi[i].ValueRe[BlokElement(NextIndex)] += aRe[i];
          psi[i].ValueIm[BlokElement(NextIndex)] += aIm[i];
        }
      }
      else
      {                         /* psi not equal to this one, if there is a next continue, otherwise add */
        if (!psi[0].Next[BlokElement(NextIndex)])
        {                       /* add a new element to psi with already used key */
          for(j=0;j<Npsi;j++) psi[j].N++;
          /* test if one has to add a new block to psi in order to store the new element of psi */
          if (psi[0].N >= BlokSize * psi[0].NBlok)
          {                           /* add a new block to psi */
            temp = (uint32_t **)realloc(psi[0].Next, (psi[0].NBlok+1) * sizeof(uint32_t *));
            if(temp==NULL)
            {                         /* Allocation failed */
              printf("alloc failed ComplexBlockWaveFunctionAddElement Next\n");
              fflush(stdout);
              /* restore psi and return 0 */
              for(j=0;j<Npsi;j++) psi[j].N--;
              return 1;
            }
            for(j=0;j<Npsi;j++) psi[j].Next = (uint32_t**)temp;
            temp = (uint32_t **)realloc(psi[0].Determinant, (psi[0].NBlok+1) * sizeof(uint32_t *));
            if(temp==NULL)
            {                         /* Allocation failed */
              printf("alloc failed ComplexBlockWaveFunctionAddElement Determinant\n");
              fflush(stdout);
              /* restore psi and return 0 */
              for(j=0;j<Npsi;j++) psi[j].N--;
              return 1;
            }
            for(j=0;j<Npsi;j++) psi[j].Determinant = (uint8_t**)temp;
            for(j=0;j<Npsi;j++)
            {
              temp = (double **)realloc(psi[j].ValueRe, (psi[j].NBlok+1) * sizeof(double *));
              if (temp == NULL)
              {                         /* Allocation failed */
                printf("alloc failed ComplexBlockWaveFunctionAddElement ValueRe %i\n",j);
                fflush(stdout);
                /* restore psi and return 0 */
                for(j=0;j<Npsi;j++) psi[j].N--;
                return (1);
              }
              psi[j].ValueRe = (double**)temp;
            }
            for(j=0;j<Npsi;j++)
            {
              temp = (double **)realloc(psi[j].ValueIm, (psi[j].NBlok+1) * sizeof(double *));
              if (temp == NULL)
              {                         /* Allocation failed */
                printf("alloc failed ComplexBlockWaveFunctionAddElement ValueRe %i\n",j);
                fflush(stdout);
                /* restore psi and return 0 */
                for(j=0;j<Npsi;j++) psi[j].N--;
                return (1);
              }
              psi[j].ValueIm = (double**)temp;
            }
            psi[0].Next[psi[0].NBlok] = (uint32_t*)malloc(BlokSize * sizeof(uint32_t));
            if (psi[0].Next[psi[0].NBlok] == NULL)
            {                         /* Allocation failed */
              printf("alloc failed ComplexBlockWaveFunctionAddElement Next[-1]\n");
              fflush(stdout);
              /* restore psi and return 0 */
              for(j=0;j<Npsi;j++) psi[j].N--;
              return (1);
            }
            psi[0].Determinant[psi[0].NBlok] = (uint8_t*)malloc(BlokSize * psi[0].NWord * sizeof(uint8_t));
            if (psi[0].Determinant[psi[0].NBlok] == NULL)
            {                         /* Allocation failed */
              printf("alloc failed ComplexBlockWaveFunctionAddElement Determinant[-1]\n");
              fflush(stdout);
              /* restore psi and return 0 */
              for(j=0;j<Npsi;j++) psi[j].N--;
              free(psi[0].Next[psi[0].NBlok]);
              return (1);
            }
            for(j=0;j<Npsi;j++)
            {
              psi[j].ValueRe[psi[0].NBlok] = (double*)malloc(BlokSize * sizeof(double));
              if(psi[j].ValueRe[psi[0].NBlok]==NULL)
              {
                printf("alloc failed ComplexBlockWaveFunctionAddElement ValueRe[%i]\n",j);
                fflush(stdout);
                /* restore psi and return 0 */
                for(k=0;k<Npsi;k++) psi[k].N--;
                free(psi[0].Next[psi[0].NBlok]);
                free(psi[0].Determinant[psi[0].NBlok]);
                for(k=0;k<j;k++) free(psi[k].ValueRe[psi[0].NBlok]);
                return (1);
              }
            }
            for(j=0;j<Npsi;j++)
            {
              psi[j].ValueIm[psi[0].NBlok] = (double*)malloc(BlokSize * sizeof(double));
              if(psi[j].ValueIm[psi[0].NBlok]==NULL)
              {
                printf("alloc failed ComplexBlockWaveFunctionAddElement ValueRe[%i]\n",j);
                fflush(stdout);
                /* restore psi and return 0 */
                for(k=0;k<Npsi;k++) psi[k].N--;
                free(psi[0].Next[psi[0].NBlok]);
                free(psi[0].Determinant[psi[0].NBlok]);
                for(k=0;k<Npsi;k++) free(psi[k].ValueRe[psi[0].NBlok]);
                for(k=0;k<j;k++) free(psi[k].ValueIm[psi[0].NBlok]);
                return (1);
              }
            }
            for(j=0;j<Npsi;j++) psi[j].NBlok++;
          }          /* Added a new block to psi */
          /* add the element */
          for(i=0;i<Npsi;i++)
          {
            psi[i].ValueRe[BlokElement(psi[0].N)] = aRe[i];
            psi[i].ValueIm[BlokElement(psi[0].N)] = aIm[i];
          }
          psi[0].Next[BlokElement(psi[0].N)] = 0;
          for (i = 0; i < psi[0].NWord; i++)
          {
            psi[0].Determinant[BlokElementDet(psi[0].N, psi[0].NWord, i)] = det[i];
          }
          psi[0].Next[BlokElement(NextIndex)] = psi[0].N;
          ThereIsANext = 0;
        }
        else
        {                       /* the key is in use and there still is at least one more element with the same key, continu checking */
          NextIndex = psi[0].Next[BlokElement(NextIndex)];
        }
      }
    }
  }
  return (0);
}
uint32_t BlockWaveFunctionAddElement(WaveFunctionType * const psi, const uint8_t * const det, const double * aRe, const double * aIm, const uint32_t Npsi)
{
  double imsum,resum;
  uint32_t i;
  if(Npsi==0) return 0;
  if(psi[0].QComplex==0)
  {
    if(aIm==NULL)
    {
      return RealBlockWaveFunctionAddElement(psi, det, aRe, Npsi);
    }
    imsum=0;
    resum=0;
    for(i=0;i<Npsi;i++)
    {
      imsum+=Sqr(aIm[i]);
      resum+=Sqr(aRe[i]);
    }
    if(imsum < 100 * Sqr(DBL_EPSILON) * resum)
    {
      return RealBlockWaveFunctionAddElement(psi, det, aRe, Npsi);
    }
    if(MakeComplexBlockWaveFunctionFromReal(psi, Npsi))
    {
      printf("MakeComplexBlockWaveFunctionFromReal failed in BlockWaveFunctionAddElement\n");
      return 1;
    }
    return ComplexBlockWaveFunctionAddElement(psi, det, aRe, aIm, Npsi);
  }
  if(aIm==NULL)
  {
    double aIm0[Npsi];
    memset(aIm0, 0.0, Npsi);
    return ComplexBlockWaveFunctionAddElement(psi, det, aRe, aIm0, Npsi);
  }
  return ComplexBlockWaveFunctionAddElement(psi, det, aRe, aIm, Npsi);
}


/*
 *   add ellement with single slater determinant det and Npsi pre-factors a[i] to Npsi wave functions psi[i]
 *   but only add it if the determinant is already in the basis
 */
void RealBlockWaveFunctionAddElementConserveBasis(WaveFunctionType const *psi, const uint8_t * const det, const double * a, const uint32_t Npsi)
{
  uint32_t i, NextIndex, PsiAIsPsiB;
  if(Npsi==0) return;
  /* test if the key that belongs to psi already excists */
  NextIndex=psi[0].KeyList[hash(det, psi[0].NWord, psi[0].NBitsKey)];
  if (!NextIndex)
  {                            /* det not in key, do nothing */
    return;
  }
  else
  {                             /* The key is already taken, test if this is due to a collision, or if psi is already in the basis set */
	while (1)
    {
      PsiAIsPsiB = 1;
      for (i = 0; i < psi[0].NWord; i++)
      {
        if (psi[0].Determinant[BlokElementDet(NextIndex, psi[0].NWord, i)] != det[i])
        {
        	PsiAIsPsiB = 0;
        }
      }
      if (PsiAIsPsiB)
      {                         /* psi already in basis set add and leave loop  */
    	for(i=0;i<Npsi;i++)
        {
#pragma omp atomic
          psi[i].ValueRe[BlokElement(NextIndex)] += a[i];
        }
        return;
      }
      else
      {                         /* psi not equal to this one, if there is a next continue, otherwise leave */
        if (!psi[0].Next[BlokElement(NextIndex)])
        {                       /* hash(psi) had a collision, but det is different. There are no other psi's with the same hash, so det not in basis lets leave */
          return;
        }
        else
        {                       /* the key is in use and there still is at least one more element with the same key, continu checking */
          NextIndex = psi[0].Next[BlokElement(NextIndex)];
        }
      }
    }
  }
}
void ComplexBlockWaveFunctionAddElementConserveBasis(WaveFunctionType const *psi, const uint8_t * const det, const double * aRe, const double * aIm, const uint32_t Npsi)
{
  uint32_t i, NextIndex, PsiAIsPsiB;
  if(Npsi==0) return;
  /* test if the key that belongs to psi already excists */
  NextIndex=psi[0].KeyList[hash(det, psi[0].NWord, psi[0].NBitsKey)];
  if (!NextIndex)
  {                             /* det not in key, do nothing */
    return;
  }
  else
  {                             /* The key is already taken, test if this is due to a collision, or if psi is already in the basis set */
    while (1)
    {
      PsiAIsPsiB = 1;
      for (i = 0; i < psi[0].NWord; i++)
      {
        if (psi[0].Determinant[BlokElementDet(NextIndex, psi[0].NWord, i)] != det[i])
        {
          PsiAIsPsiB = 0;
        }
      }
      if (PsiAIsPsiB)
      {                         /* psi already in basis set add and leave loop  */
        for(i=0;i<Npsi;i++)
        {
#pragma omp atomic
          psi[i].ValueRe[BlokElement(NextIndex)] += aRe[i];
#pragma omp atomic
          psi[i].ValueIm[BlokElement(NextIndex)] += aIm[i];
        }
        return;
      }
      else
      {                         /* psi not equal to this one, if there is a next continue, otherwise leave */
        if (!psi[0].Next[BlokElement(NextIndex)])
        {                       /* hash(psi) had a collision, but det is different. There are no other psi's with the same hash, so det not in basis lets leave */
          return;
        }
        else
        {                       /* the key is in use and there still is at least one more element with the same key, continu checking */
          NextIndex = psi[0].Next[BlokElement(NextIndex)];
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
uint32_t RealBlockWaveFunctionAddElementOMPMiniFlush(WaveFunctionType * const psi, uint32_t * NMini, const uint8_t * const MiniDet, const double *const MiniValue, const uint32_t Npsi, const bool ExpandBasis)
{
  uint32_t i, NoMoreMem = 0;
#pragma omp critical (AddToPsi)
  {
    for (i = 0; i < *NMini; i++)
    {
      if (ExpandBasis)
      {
    	  if (RealBlockWaveFunctionAddElement(psi, &MiniDet[i * psi->NWord], &(MiniValue[i*Npsi]), Npsi))
		  {
			printf("RealBlockWaveFunctionAddElement failed in RealBlockWaveFunctionAddElementOMPMiniFlush\n");
			fflush(stdout);
			NoMoreMem = 1;
			i = *NMini;
			*NMini = 0;
		  }
      }
      else
      {
    	  RealBlockWaveFunctionAddElementConserveBasis(psi, &MiniDet[i * psi->NWord], &(MiniValue[i*Npsi]), Npsi);
      }
    }
#pragma omp flush
  }
  *NMini = 0;
  return NoMoreMem;
}
uint32_t ComplexBlockWaveFunctionAddElementOMPMiniFlush(WaveFunctionType * const psi, uint32_t * NMini, const uint8_t * const MiniDet, const double *const MiniValueRe, const double *const MiniValueIm, const uint32_t Npsi)
{
  uint32_t i, NoMoreMem = 0;
#pragma omp critical (AddToPsi)
  {
    for (i = 0; i < *NMini; i++)
    {
      if (ComplexBlockWaveFunctionAddElement(psi, &MiniDet[i * psi->NWord], &(MiniValueRe[i*Npsi]), &(MiniValueIm[i*Npsi]), Npsi))
      {
        printf("ComplexBlockWaveFunctionAddElement failed in ComplexBlockWaveFunctionAddElementOMPMiniFlush\n");
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





/*   a and b are block Wavefunctions of the same length Npsi
 *   (a+new allocation) <- a+b
 *   a[i] <- a[i] + b[i]
 *   This function leaves b in the memory and adds values to a
 *   When aditional space for a needs to be allocated it will
 *   Therefor values of a might change and a is returned from the function
 *   If a is an empty function b is copied into a
 *   return 0 if succeded
 *   return 1 if not enough memory in order to allocate a
 *   If not enough memory a is not changed, but might have a larger basis
 */
uint32_t RealBlockWaveFunctionAddUnionBasis(WaveFunctionType * const a, const WaveFunctionType * const b, const uint32_t Npsi)
{
  uint32_t i, j, k;
  double bVRe[Npsi];
  for (i = b[0].N; i != 0; --i)
  {
    for(j=0;j<Npsi;j++)
    {
      bVRe[j] = b[j].ValueRe[BlokElement(i)];
    }
    if (RealBlockWaveFunctionAddElement(a, &b[0].Determinant[BlokElementDet(i, b[0].NWord, 0)], bVRe, Npsi) != 0)
    {                           /* Allocation of new memory for a failed Remove all elements added so far and return 1 */
      for (j = b[0].N; j != i; --j)
      {
        for(k=0;k<Npsi;k++)
        {
          bVRe[j] = -b[k].ValueRe[BlokElement(i)];
        }
        RealBlockWaveFunctionAddElementConserveBasis(a, &b[0].Determinant[BlokElementDet(i, b[0].NWord, 0)],bVRe, Npsi);
      }
      return 1;
    }
  }
  return 0;
}
uint32_t ComplexBlockWaveFunctionAddUnionBasis(WaveFunctionType * const a, const WaveFunctionType * const b, const uint32_t Npsi)
{
  uint32_t i, j, k;
  double bVRe[Npsi];
  double bVIm[Npsi];
  for (i = b[0].N; i != 0; --i)
  {
    for(j=0;j<Npsi;j++)
    {
      bVRe[j] = b[j].ValueRe[BlokElement(i)];
      bVIm[j] = b[j].ValueIm[BlokElement(i)];
    }
    if (ComplexBlockWaveFunctionAddElement(a, &b[0].Determinant[BlokElementDet(i, b[0].NWord, 0)], bVRe, bVIm, Npsi) != 0)
    {                           /* Allocation of new memory for a failed Remove all elements added so far and return 1 */
      for (j = b[0].N; j != i; --j)
      {
        for(k=0;k<Npsi;k++)
        {
          bVRe[j] = -b[k].ValueRe[BlokElement(i)];
          bVIm[j] = -b[k].ValueIm[BlokElement(i)];
        }
        ComplexBlockWaveFunctionAddElementConserveBasis(a, &b[0].Determinant[BlokElementDet(i, b[0].NWord, 0)],bVRe,bVIm, Npsi);
      }
      return 1;
    }
  }
  return 0;
}

/*   a and b are block Wavefunctions that might have different length, alpha is a matrix
 *   (a+new allocation) <- a+alpha^* b
 *   a[i] = a[i] + sum_j alpha[i][j]^* b[j]
 *   This function leaves b in the memory and adds values to a
 *   When aditional space for a needs to be allocated it will
 *   Therefor values of a might change and a is returned from the function
 *   If a is an empty function b is copied into a (use of copy is strongly recommanded)
 */
uint32_t RealBlockWaveFunctionAddMultipliedUnionBasis(WaveFunctionType * const a, const uint32_t Npsia, const WaveFunctionType * const b, const uint32_t Npsib, const CompactMatrixType alpha)
{
  uint32_t i, j, k;
  double bVRe[Npsia];
  for(k=b[0].N;k!=0;--k)
  {
    for(i=0;i<Npsia;i++)
    {
      bVRe[i]=0;
      for(j=0;j<Npsib;j++)
      {
        bVRe[i] += alpha.VRe[i*alpha.Nj+j] * b[j].ValueRe[BlokElement(k)];
      }
    }
    if(RealBlockWaveFunctionAddElement(a, &b[0].Determinant[BlokElementDet(k, b[0].NWord, 0)], bVRe, Npsia))
    {                          /* Allocation of new memory for a failed Remove all elements added so far and return 1 */
      for(k=k+1;k<=b[0].N;k++)
      {
        for(i=0;i<Npsia;i++)
        {
          bVRe[i]=0;
          for(j=0;j<Npsib;j++)
          {
            bVRe[i] -= alpha.VRe[i*alpha.Nj+j] * b[j].ValueRe[BlokElement(k)];
          }
        }
        RealBlockWaveFunctionAddElementConserveBasis(a, &b[0].Determinant[BlokElementDet(k, b[0].NWord, 0)], bVRe, Npsia);
      }
      return 1;
    }
  }
  return 0;
}
uint32_t ComplexBlockWaveFunctionAddMultipliedUnionBasis(WaveFunctionType * const a, const uint32_t Npsia, const WaveFunctionType * const b, const uint32_t Npsib, const CompactMatrixType alpha)
{
  uint32_t i, j, k;
  double bVRe[Npsia];
  double bVIm[Npsia];
  for(k=b[0].N;k!=0;--k)
  {
    for(i=0;i<Npsia;i++)
    {
      bVRe[i]=0;
      bVIm[i]=0;
      for(j=0;j<Npsib;j++)
      {
        bVRe[i] += alpha.V[i*alpha.Nj+j].Re * b[j].ValueRe[BlokElement(k)] + alpha.V[i*alpha.Nj+j].Im * b[j].ValueIm[BlokElement(k)];
        bVIm[i] += alpha.V[i*alpha.Nj+j].Re * b[j].ValueIm[BlokElement(k)] - alpha.V[i*alpha.Nj+j].Im * b[j].ValueRe[BlokElement(k)];
      }
    }
    if(ComplexBlockWaveFunctionAddElement(a, &b[0].Determinant[BlokElementDet(k, b[0].NWord, 0)], bVRe, bVIm, Npsia))
    {                          /* Allocation of new memory for a failed Remove all elements added so far and return 1 */
      for(k=k+1;k<=b[0].N;k++)
      {
        for(i=0;i<Npsia;i++)
        {
          bVRe[i]=0;
          bVIm[i]=0;
          for(j=0;j<Npsib;j++)
          {
            bVRe[i] -= alpha.V[i*alpha.Nj+j].Re * b[j].ValueRe[BlokElement(k)] + alpha.V[i*alpha.Nj+j].Im * b[j].ValueIm[BlokElement(k)];
            bVIm[i] -= alpha.V[i*alpha.Nj+j].Re * b[j].ValueIm[BlokElement(k)] - alpha.V[i*alpha.Nj+j].Im * b[j].ValueRe[BlokElement(k)];
          }
        }
        ComplexBlockWaveFunctionAddElementConserveBasis(a, &b[0].Determinant[BlokElementDet(k, b[0].NWord, 0)], bVRe, bVIm, Npsia);
      }
      return 1;
    }
  }
  return 0;
}

/*   a and b are block Wavefunctions that might have different length, alpha is a matrix
 *   (a+new allocation) <- a-alpha^* b
 *   a[i] = a[i] - sum_j alpha[i][j]^* b[j]
 *   This function leaves b in the memory and adds values to a
 *   When aditional space for a needs to be allocated it will
 *   Therefor values of a might change and a is returned from the function
 *   If a is an empty function b is copied into a (use of copy is strongly recommanded)
 */
uint32_t RealBlockWaveFunctionSubtractMultipliedUnionBasis(WaveFunctionType * const a, const uint32_t Npsia, const WaveFunctionType * const b, const uint32_t Npsib, const CompactMatrixType alpha)
{
  uint32_t i, j, k;
  double bVRe[Npsia];
  if(Npsib==0) return 0;
  for(k=b[0].N;k!=0;--k)
  {
    for(i=0;i<Npsia;i++)
    {
      bVRe[i]=0;
      for(j=0;j<Npsib;j++)
      {
        bVRe[i] -= alpha.VRe[i*alpha.Nj+j] * b[j].ValueRe[BlokElement(k)];
      }
    }
    if(RealBlockWaveFunctionAddElement(a, &b[0].Determinant[BlokElementDet(k, b[0].NWord, 0)], bVRe, Npsia))
    {                          /* Allocation of new memory for a failed Remove all elements added so far and return 1 */
      for(k=k+1;k<=b[0].N;k++)
      {
        for(i=0;i<Npsia;i++)
        {
          bVRe[i]=0;
          for(j=0;j<Npsib;j++)
          {
            bVRe[i] += alpha.VRe[i*alpha.Nj+j] * b[j].ValueRe[BlokElement(k)];
          }
        }
        RealBlockWaveFunctionAddElementConserveBasis(a, &b[0].Determinant[BlokElementDet(k, b[0].NWord, 0)], bVRe, Npsia);
      }
      return 1;
    }
  }
  return 0;
}
uint32_t ComplexBlockWaveFunctionSubtractMultipliedUnionBasis(WaveFunctionType * const a, const uint32_t Npsia, const WaveFunctionType * const b, const uint32_t Npsib, const CompactMatrixType alpha)
{
  uint32_t i, j, k;
  double bVRe[Npsia];
  double bVIm[Npsia];
  if(Npsib==0) return 0;
  for(k=b[0].N;k!=0;--k)
  {
    for(i=0;i<Npsia;i++)
    {
      bVRe[i]=0;
      bVIm[i]=0;
      for(j=0;j<Npsib;j++)
      {
        bVRe[i] -= alpha.V[i*alpha.Nj+j].Re * b[j].ValueRe[BlokElement(k)] + alpha.V[i*alpha.Nj+j].Im * b[j].ValueIm[BlokElement(k)];
        bVIm[i] -= alpha.V[i*alpha.Nj+j].Re * b[j].ValueIm[BlokElement(k)] - alpha.V[i*alpha.Nj+j].Im * b[j].ValueRe[BlokElement(k)];
      }
    }
    if(ComplexBlockWaveFunctionAddElement(a, &b[0].Determinant[BlokElementDet(k, b[0].NWord, 0)], bVRe, bVIm, Npsia))
    {                          /* Allocation of new memory for a failed Remove all elements added so far and return 1 */
      for(k=k+1;k<=b[0].N;k++)
      {
        for(i=0;i<Npsia;i++)
        {
          bVRe[i]=0;
          bVIm[i]=0;
          for(j=0;j<Npsib;j++)
          {
            bVRe[i] += alpha.V[i*alpha.Nj+j].Re * b[j].ValueRe[BlokElement(k)] + alpha.V[i*alpha.Nj+j].Im * b[j].ValueIm[BlokElement(k)];
            bVIm[i] += alpha.V[i*alpha.Nj+j].Re * b[j].ValueIm[BlokElement(k)] - alpha.V[i*alpha.Nj+j].Im * b[j].ValueRe[BlokElement(k)];
          }
        }
        ComplexBlockWaveFunctionAddElementConserveBasis(a, &b[0].Determinant[BlokElementDet(k, b[0].NWord, 0)], bVRe, bVIm, Npsia);
      }
      return 1;
    }
  }
  return 0;
}



/*
 *  remove determinants from the basis that do not fullfill the restrictions
 */
void RealBlockWaveFunctionRestrict(WaveFunctionType * const a, const uint32_t Npsi, const DeterminantRestrictionType res)
{
  uint32_t i,startN, iBlok, iElement;
  double aval[Npsi];
  if(Npsi==0) return ;
  startN = a[0].N;
  for(i=0;i<Npsi;i++)
  {
    a[i].ValueRe[0][0] = 0;
    a[i].N=0;
  }
  for (iElement = 0; iElement < ((uint32_t) 1 << (a[0].NBitsKey)); ++iElement)
  {                             /* Whipe out the entire KeyList */
    a[0].KeyList[iElement] = 0;
  }
  for (iBlok = 0; iBlok < a[0].NBlok - 1; iBlok++)
  {
    for (iElement = 0; iElement < BlokSize; iElement++)
    {
      a[0].Next[iBlok][iElement] = 0;
      if (QDetNotRestricted(&(a[0].Determinant[iBlok][iElement*a[0].NWord]), a[0].NFermion, a[0].NBoson, a[0].NWord, res))
      {
        for(i=0;i<Npsi;i++)
        {
          aval[i] = a[i].ValueRe[iBlok][iElement];
        }
        if (RealBlockWaveFunctionAddElement(a, &a[0].Determinant[iBlok][iElement * a[0].NWord], aval, Npsi)) AllocFailed();
      }
    }
  }
  iBlok = a->NBlok - 1;
  for (iElement = 0; iElement < 1 + startN - iBlok * BlokSize; ++iElement)
  {
    a[0].Next[iBlok][iElement] = 0;
    if (QDetNotRestricted(&(a[0].Determinant[iBlok][iElement*a->NWord]), a[0].NFermion, a[0].NBoson, a[0].NWord, res))
    {
      for(i=0;i<Npsi;i++)
      {
        aval[i] = a[i].ValueRe[iBlok][iElement];
      }
      if (RealBlockWaveFunctionAddElement(a, &a[0].Determinant[iBlok][iElement * a[0].NWord], aval, Npsi)) AllocFailed();
    }
  }
  for (iBlok = 1 + (a[0].N >> NBitsBlok); iBlok < a[0].NBlok; iBlok++)
  {
    for(i=0;i<Npsi;i++) free(a[i].ValueRe[iBlok]);
    free(a[0].Next[iBlok]);
    free(a[0].Determinant[iBlok]);
  }
  for(i=0;i<Npsi;i++)
  {
    a[i].NBlok = 1 + (a[0].N >> NBitsBlok);
    a[i].ValueRe = (double **)realloc(a[i].ValueRe, a[0].NBlok * sizeof(double *));
  }
  a[0].Next = (uint32_t **)realloc(a[0].Next, a[0].NBlok * sizeof(uint32_t *));
  a[0].Determinant = (uint8_t **)realloc(a[0].Determinant, a[0].NBlok * sizeof(uint8_t *));
}
void ComplexBlockWaveFunctionRestrict(WaveFunctionType * const a, const uint32_t Npsi, const DeterminantRestrictionType res)
{
  uint32_t i,startN, iBlok, iElement;
  double avalRe[Npsi];
  double avalIm[Npsi];
  if(Npsi==0) return;
  startN = a[0].N;
  for(i=0;i<Npsi;i++)
  {
    a[i].ValueRe[0][0] = 0;
    a[i].ValueIm[0][0] = 0;
    a[i].N=0;
  }
  for (iElement = 0; iElement < ((uint32_t) 1 << (a[0].NBitsKey)); ++iElement)
  {                             /* Whipe out the entire KeyList */
    a[0].KeyList[iElement] = 0;
  }
  for (iBlok = 0; iBlok < a[0].NBlok - 1; iBlok++)
  {
    for (iElement = 0; iElement < BlokSize; iElement++)
    {
      a[0].Next[iBlok][iElement] = 0;
      if (QDetNotRestricted(&(a[0].Determinant[iBlok][iElement*a[0].NWord]), a[0].NFermion, a[0].NBoson, a[0].NWord, res))
      {
        for(i=0;i<Npsi;i++)
        {
          avalRe[i] = a[i].ValueRe[iBlok][iElement];
          avalIm[i] = a[i].ValueIm[iBlok][iElement];
        }
        if (ComplexBlockWaveFunctionAddElement(a, &a[0].Determinant[iBlok][iElement * a[0].NWord], avalRe, avalIm, Npsi)) AllocFailed();
      }
    }
  }
  iBlok = a->NBlok - 1;
  for (iElement = 0; iElement < 1 + startN - iBlok * BlokSize; ++iElement)
  {
    a[0].Next[iBlok][iElement] = 0;
    if (QDetNotRestricted(&(a[0].Determinant[iBlok][iElement*a->NWord]), a[0].NFermion, a[0].NBoson, a[0].NWord, res))
    {
      for(i=0;i<Npsi;i++)
      {
        avalRe[i] = a[i].ValueRe[iBlok][iElement];
        avalIm[i] = a[i].ValueIm[iBlok][iElement];
      }
      if (ComplexBlockWaveFunctionAddElement(a, &a[0].Determinant[iBlok][iElement * a[0].NWord], avalRe, avalIm, Npsi)) AllocFailed();
    }
  }
  for (iBlok = 1 + (a[0].N >> NBitsBlok); iBlok < a[0].NBlok; iBlok++)
  {
    for(i=0;i<Npsi;i++)
    {
      free(a[i].ValueRe[iBlok]);
      free(a[i].ValueIm[iBlok]);
    }
    free(a[0].Next[iBlok]);
    free(a[0].Determinant[iBlok]);
  }
  for(i=0;i<Npsi;i++)
  {
    a[i].NBlok = 1 + (a[0].N >> NBitsBlok);
    a[i].ValueRe = (double **)realloc(a[i].ValueRe, a[0].NBlok * sizeof(double *));
  }
  a[0].Next = (uint32_t **)realloc(a[0].Next, a[0].NBlok * sizeof(uint32_t *));
  a[0].Determinant = (uint8_t **)realloc(a[0].Determinant, a[0].NBlok * sizeof(uint8_t *));
}
void BlockWaveFunctionRestrict(WaveFunctionType * const a, const uint32_t Npsi, const DeterminantRestrictionType res)
{
  if(Npsi==0) return;
  if(a[0].QComplex==0)
  {
    RealBlockWaveFunctionRestrict(a, Npsi, res);
  }
  else
  {
    ComplexBlockWaveFunctionRestrict(a, Npsi, res);
  }
}




/*  Calculate Opsi[i] = H psi[i] for i=0 to Npsi
 *  on input all psi[i]'s are on the same basis
 *  on excit all Opsi[i]'s are on the same basis as psi[i]'s
 *  memory for values of Opsi is already allocated.
 *
 *  The serial version multiplies each wave-function of the Block seperately and only creates the final basis of each state after the first multiplication
 *  The version below multiplies all blocks simultanious
 *  The second variant needs less flops, but much more data in memory (larger cash, larger ) which one is faster depends (needs system and problem dependent testing)
 */
uint32_t BlockOperatorPsiSerialConserveBasisPreAllocated(const OperatorType H, WaveFunctionType const * const psi, WaveFunctionType * const Opsi, uint32_t const Npsi)
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
uint32_t BlockOperatorPsiSerialRestricted(const OperatorType H, WaveFunctionType const * const psi, WaveFunctionType * const Opsi, uint32_t const Npsi, DeterminantRestrictionType const Restrictions)
{
  uint32_t i,j;
  VerbosityPrintf(0x0002,"Start of BlockOperatorPsiSerialRestricted\n");
  if(Npsi==0) return 0;
  if(OperatorPsiRestricted(H, psi[0], &(Opsi[0]), Restrictions))
  {
    printf("OperatorPsi failed in BlockOperatorPsiSerial\n");
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
      printf("WaveFunctionInitCopyBasis failed in BlockOperatorPsiSerial\n");
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
uint32_t BlockOperatorPsiSerial(const OperatorType H, WaveFunctionType const * const psi, WaveFunctionType * const Opsi, uint32_t const Npsi)
{
  uint32_t i,j;
  if(H.Res.N!=0) return BlockOperatorPsiSerialRestricted(H, psi, Opsi, Npsi, H.Res);
  VerbosityPrintf(0x0002,"Start of BlockOperatorPsiSerial\n");
  if(Npsi==0) return 0;
  if(OperatorPsi(H, psi[0], &(Opsi[0])))
  {
    printf("OperatorPsi failed in BlockOperatorPsiSerial\n");
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
      printf("WaveFunctionInitCopyBasis failed in BlockOperatorPsiSerial\n");
      return 1;
    }
  }
  for(i=1;i<Npsi;i++)
  {
    OperatorPsiConserveBasisPreAllocated(&H, &(psi[i]), &(Opsi[i]));
  }
  return 0;
}

/*  Calculate Opsi[i] = H psi[i] for i=0 to Npsi using restrictions
 *  on input all psi[i]'s are on the same basis
 *  on excit all Opsi[i]'s are on the same basis
 */
uint32_t BlockOperatorPsiRestrictedRR(const OperatorType O, WaveFunctionType const * const psi, WaveFunctionType * const Opsi, DeterminantRestrictionType const Restrictions, uint32_t const Npsi)
{
  //TODO: replace the serial version by a version that runs all states in a block simultanious
  return BlockOperatorPsiSerialRestricted(O, psi, Opsi, Npsi, Restrictions);
}
uint32_t BlockOperatorPsiRestrictedRC(const OperatorType O, WaveFunctionType const * const psi, WaveFunctionType * const Opsi, DeterminantRestrictionType const Restrictions, uint32_t const Npsi)
{
  //TODO: replace the serial version by a version that runs all states in a block simultanious
  return BlockOperatorPsiSerialRestricted(O, psi, Opsi, Npsi, Restrictions);
}
uint32_t BlockOperatorPsiRestrictedCR(const OperatorType O, WaveFunctionType const * const psi, WaveFunctionType * const Opsi, DeterminantRestrictionType const Restrictions, uint32_t const Npsi)
{
  //TODO: replace the serial version by a version that runs all states in a block simultanious
  return BlockOperatorPsiSerialRestricted(O, psi, Opsi, Npsi, Restrictions);
}
uint32_t BlockOperatorPsiRestrictedCC(const OperatorType O, WaveFunctionType const * const psi, WaveFunctionType * const Opsi, DeterminantRestrictionType const Restrictions, uint32_t const Npsi)
{
  //TODO: replace the serial version by a version that runs all states in a block simultanious
  return BlockOperatorPsiSerialRestricted(O, psi, Opsi, Npsi, Restrictions);
}
uint32_t BlockOperatorPsiRestrictedMR(const OperatorType O, WaveFunctionType const * const psi, WaveFunctionType * const Opsi, DeterminantRestrictionType const Restrictions, uint32_t const Npsi)
{
  //TODO: replace the serial version by a version that runs all states in a block simultanious
  return BlockOperatorPsiSerialRestricted(O, psi, Opsi, Npsi, Restrictions);
}
uint32_t BlockOperatorPsiRestrictedMC(const OperatorType O, WaveFunctionType const * const psi, WaveFunctionType * const Opsi, DeterminantRestrictionType const Restrictions, uint32_t const Npsi)
{
  //TODO: replace the serial version by a version that runs all states in a block simultanious
  return BlockOperatorPsiSerialRestricted(O, psi, Opsi, Npsi, Restrictions);
}
uint32_t BlockOperatorPsiRestricted(const OperatorType O, WaveFunctionType const * const psi, WaveFunctionType * const Opsi, DeterminantRestrictionType const Restrictions, uint32_t const Npsi)
{
  if(Npsi==0) return 0;
  if (O.QComplex == 0)
  {
    if (psi[0].QComplex)
    {
      return BlockOperatorPsiRestrictedRC(O, psi, Opsi, Restrictions, Npsi);
    }
    return BlockOperatorPsiRestrictedRR(O, psi, Opsi, Restrictions, Npsi);
  }
  if (O.QComplex == 1)
  {
    if (psi[0].QComplex)
    {
      return BlockOperatorPsiRestrictedCC(O, psi, Opsi, Restrictions, Npsi);
    }
    return BlockOperatorPsiRestrictedCR(O, psi, Opsi, Restrictions, Npsi);
  }
  if (psi[0].QComplex)
  {
    return BlockOperatorPsiRestrictedMC(O, psi, Opsi, Restrictions, Npsi);
  }
  return BlockOperatorPsiRestrictedMR(O, psi, Opsi, Restrictions, Npsi);
}

/*  Calculate Opsi[i] = H psi[i] for i=0 to Npsi
 *  on input all psi[i]'s are on the same basis
 *  on excit all Opsi[i]'s are on the same basis
 */
uint32_t BlockOperatorPsiRR(const OperatorType O, WaveFunctionType const * const psi, WaveFunctionType * const Opsi, uint32_t const Npsi, const bool ExpandBasis)
{
  uint32_t i, j, k, l, m, ii, ij, ik, OutOffMem;
  int    iOMP;
  double phase;
  if(Npsi==0) return 0;
  Opsi[0].QComplex = 0;
  Opsi[0].NFermion = psi[0].NFermion;
  Opsi[0].NBoson = psi[0].NBoson;
  Opsi[0].NBitsKey = psi[0].NBitsKey;
  if (!ExpandBasis)
  {
	  WaveFunctionInitCopyBasis(&(Opsi[0]), psi[0]);
  }
  else
  {
	  if (WaveFunctionInit(&(Opsi[0])))
	  {
		printf("Could not initialize Opsi[0] in BlockOperatorPsiRR\n");
		return 1;
	  };
  }
  for(i=1;i<Npsi;i++)
  {
    if(WaveFunctionInitCopyBasis(&(Opsi[i]), Opsi[0]))
    {
      printf("Could not initialize Opsi[%i] in BlockOperatorPsiRR\n",i);
      for(j=1;j<i;j++)
      {
        WaveFunctionFreeValues(&(Opsi[j]));
      }
      WaveFunctionFree(&(Opsi[0]));
      return 1;
    }
  }
  OutOffMem = 0;
#pragma omp parallel private(i,j,k,l,m,ii,ij,ik,phase)
  {
    uint8_t Det[Opsi->NWord];
    uint8_t MiniDet[MiniBlokSize * Opsi->NWord];
    double MiniValue[MiniBlokSize][Npsi];
    uint32_t NMini;
    NMini = 0;
#pragma omp for
    for (iOMP = 1; iOMP <= psi[0].N; iOMP++)
    {
      if (!OutOffMem)
      {
        i = (uint32_t) iOMP;
        ii = i >> NBitsBlok;
        ik = i & BlokMask;
        ij = psi[0].NWord * ik;
        for (l = 0; l <= O.MaxLength; l++)
        {
          for (j = 0; j < O.l[l].N; j++)
          {
            for (k = 0; k < Opsi->NWord; k++)
              Det[k] = psi[0].Determinant[ii][ij + k];
            phase = OperatorLader(Det, &(O.l[l].o[j * l]), l, O.NFermion);
            if (fabs(phase) > (10 * DBL_MIN))
            {
              for(m=0;m<Npsi;m++) MiniValue[NMini][m] = phase * O.l[l].VRe[j] * psi[m].ValueRe[ii][ik];
              for (k = 0; k < Opsi->NWord; k++)
                MiniDet[NMini * Opsi->NWord + k] = Det[k];
              NMini++;
              if (NMini == MiniBlokSize)
              {                 /* flush Mini buffer */
                if (RealBlockWaveFunctionAddElementOMPMiniFlush(Opsi, &NMini, MiniDet, MiniValue[0], Npsi, ExpandBasis))
                {               /* Allocation failed */
                  printf("RealWaveFunctionAddElementOMPMiniFlush failed in BlockOperatorPsiRR\n 01");
                  fflush(stdout);
                  OutOffMem = 1;
#pragma omp flush(OutOffMem)
                  j = O.l[l].N;
                  l = O.MaxLength;
                }
              }
            }
          }                     /* end of loop over different operators of lenght l */
        }                       /* end of loop over operators of different length */
      }                         /* Out Off Mem test */
    }                           /* end of loop over elements of psi this loop runs parallel */
    if (RealBlockWaveFunctionAddElementOMPMiniFlush(Opsi, &NMini, MiniDet, MiniValue[0], Npsi, ExpandBasis))
    {                           /* Allocation failed */
      printf("RealBlockWaveFunctionAddElementOMPMiniFlush failed in BlockOperatorPsiRR 02\n");
      fflush(stdout);
      OutOffMem = 1;
#pragma omp flush(OutOffMem)
    }
  }                             /* end of parallel region */
  if (OutOffMem)
  {
    for(i=1;i<Npsi;i++)
    {
      RealWaveFunctionFreeValues(&(Opsi[i]));
    }
    RealWaveFunctionFree(&(Opsi[0]));
    return 1;
  }
  return 0;
}
uint32_t BlockOperatorPsiRC(const OperatorType O, WaveFunctionType const * const psi, WaveFunctionType * const Opsi, uint32_t const Npsi)
{
  uint32_t i, j, k, l, m, ii, ij, ik, OutOffMem;
  int    iOMP;
  double phase;
  if(Npsi==0) return 0;
  Opsi[0].QComplex = 1;
  Opsi[0].NFermion = psi[0].NFermion;
  Opsi[0].NBoson = psi[0].NBoson;
  Opsi[0].NBitsKey = psi[0].NBitsKey;
  if (WaveFunctionInit(&(Opsi[0])))
  {
    printf("Could not initialize Opsi[0] in BlockOperatorPsiRC\n");
    return 1;
  };
  for(i=1;i<Npsi;i++)
  {
    if(WaveFunctionInitCopyBasis(&(Opsi[i]), Opsi[0]))
    {
      printf("Could not initialize Opsi[%i] in BlockOperatorPsiRC\n",i);
      for(j=1;j<i;j++)
      {
        WaveFunctionFreeValues(&(Opsi[j]));
      }
      WaveFunctionFree(&(Opsi[0]));
      return 1;
    }
  }
  OutOffMem = 0;
#pragma omp parallel private(i,j,k,l,m,ii,ij,ik,phase)
  {
    uint8_t Det[Opsi->NWord];
    uint8_t MiniDet[MiniBlokSize * Opsi->NWord];
    double MiniValueRe[MiniBlokSize][Npsi];
    double MiniValueIm[MiniBlokSize][Npsi];
    uint32_t NMini = 0;
#pragma omp for
    for (iOMP = 1; iOMP <= psi[0].N; iOMP++)
    {
      if (!OutOffMem)
      {
        i = (uint32_t) iOMP;
        ii = i >> NBitsBlok;
        ik = i & BlokMask;
        ij = psi[0].NWord * ik;
        for (l = 0; l <= O.MaxLength; l++)
        {
          for (j = 0; j < O.l[l].N; j++)
          {
            for (k = 0; k < Opsi->NWord; k++)
              Det[k] = psi[0].Determinant[ii][ij + k];
            phase = OperatorLader(Det, &(O.l[l].o[j * l]), l, O.NFermion);
            if (fabs(phase) > (10 * DBL_MIN))
            {
              for(m=0;m<Npsi;m++)
              {
                MiniValueRe[NMini][m] = phase * O.l[l].VRe[j] * psi[m].ValueRe[ii][ik];
                MiniValueIm[NMini][m] = phase * O.l[l].VRe[j] * psi[m].ValueIm[ii][ik];
              }
              for (k = 0; k < Opsi->NWord; k++) MiniDet[NMini * Opsi->NWord + k] = Det[k];
              NMini++;
              if (NMini == MiniBlokSize)
              {                 /* flush Mini buffer */
                if (ComplexBlockWaveFunctionAddElementOMPMiniFlush(Opsi, &NMini, MiniDet, MiniValueRe[0], MiniValueIm[0], Npsi))
                {               /* Allocation failed */
                  printf("ComplexBlockWaveFunctionAddElementOMPMiniFlush failed in BlockOperatorPsiRC\n 01");
                  fflush(stdout);
                  OutOffMem = 1;
#pragma omp flush(OutOffMem)
                  j = O.l[l].N;
                  l = O.MaxLength;
                }
              }
            }
          }                     /* end of loop over different operators of lenght l */
        }                       /* end of loop over operators of different length */
      }                         /* Out Off Mem test */
    }                           /* end of loop over elements of psi this loop runs parallel */
    if (ComplexBlockWaveFunctionAddElementOMPMiniFlush(Opsi, &NMini, MiniDet, MiniValueRe[0], MiniValueIm[0], Npsi))
    {                           /* Allocation failed */
      printf("ComplexBlockWaveFunctionAddElementOMPMiniFlush failed in BlockOperatorPsiRC 02\n");
      fflush(stdout);
      OutOffMem = 1;
#pragma omp flush(OutOffMem)
    }
  }                             /* end of parallel region */
  if (OutOffMem)
  {
    for(i=1;i<Npsi;i++)
    {
      ComplexWaveFunctionFreeValues(&(Opsi[i]));
    }
    ComplexWaveFunctionFree(&(Opsi[0]));
    return 1;
  }
  return 0;
}
uint32_t BlockOperatorPsiCR(const OperatorType O, WaveFunctionType const * const psi, WaveFunctionType * const Opsi, uint32_t const Npsi)
{
  uint32_t i, j, k, l, m, ii, ij, ik, OutOffMem;
  int    iOMP;
  double phase;
  if(Npsi==0) return 0;
  Opsi[0].QComplex = 1;
  Opsi[0].NFermion = psi[0].NFermion;
  Opsi[0].NBoson = psi[0].NBoson;
  Opsi[0].NBitsKey = psi[0].NBitsKey;
  if (WaveFunctionInit(&(Opsi[0])))
  {
    printf("Could not initialize Opsi[0] in BlockOperatorPsiCR\n");
    return 1;
  };
  for(i=1;i<Npsi;i++)
  {
    if(WaveFunctionInitCopyBasis(&(Opsi[i]), Opsi[0]))
    {
      printf("Could not initialize Opsi[%i] in BlockOperatorPsiCR\n",i);
      for(j=1;j<i;j++)
      {
        WaveFunctionFreeValues(&(Opsi[j]));
      }
      WaveFunctionFree(&(Opsi[0]));
      return 1;
    }
  }
  OutOffMem = 0;
#pragma omp parallel private(i,j,k,l,m,ii,ij,ik,phase)
  {
    uint8_t Det[Opsi->NWord];
    uint8_t MiniDet[MiniBlokSize * Opsi->NWord];
    double MiniValueRe[MiniBlokSize][Npsi];
    double MiniValueIm[MiniBlokSize][Npsi];
    uint32_t NMini;
    NMini = 0;
#pragma omp for
    for (iOMP = 1; iOMP <= psi[0].N; iOMP++)
    {
      if (!OutOffMem)
      {
        i = (uint32_t) iOMP;
        ii = i >> NBitsBlok;
        ik = i & BlokMask;
        ij = psi[0].NWord * ik;
        for (l = 0; l <= O.MaxLength; l++)
        {
          for (j = 0; j < O.l[l].N; j++)
          {
            for (k = 0; k < Opsi->NWord; k++) Det[k] = psi[0].Determinant[ii][ij + k];
            phase = OperatorLader(Det, &(O.l[l].o[j * l]), l, O.NFermion);
            if (fabs(phase) > (10 * DBL_MIN))
            {
              for(m=0;m<Npsi;m++)
              {
                MiniValueRe[NMini][m] = phase * O.l[l].V[j].Re * psi[m].ValueRe[ii][ik];
                MiniValueIm[NMini][m] = phase * O.l[l].V[j].Im * psi[m].ValueRe[ii][ik];
              }
              for (k = 0; k < Opsi->NWord; k++) MiniDet[NMini * Opsi->NWord + k] = Det[k];
              NMini++;
              if (NMini == MiniBlokSize)
              {                 /* flush Mini buffer */
                if (ComplexBlockWaveFunctionAddElementOMPMiniFlush(Opsi, &NMini, MiniDet, MiniValueRe[0], MiniValueIm[0], Npsi))
                {               /* Allocation failed */
                  printf("ComplexBlockWaveFunctionAddElementOMPMiniFlush failed in BlockOperatorPsiCR\n 01");
                  fflush(stdout);
                  OutOffMem = 1;
#pragma omp flush(OutOffMem)
                  j = O.l[l].N;
                  l = O.MaxLength;
                }
              }
            }
          }                     /* end of loop over different operators of lenght l */
        }                       /* end of loop over operators of different length */
      }                         /* Out Off Mem test */
    }                           /* end of loop over elements of psi this loop runs parallel */
    if (ComplexBlockWaveFunctionAddElementOMPMiniFlush(Opsi, &NMini, MiniDet, MiniValueRe[0], MiniValueIm[0], Npsi))
    {                           /* Allocation failed */
      printf("ComplexBlockWaveFunctionAddElementOMPMiniFlush failed in BlockOperatorPsiCR 02\n");
      fflush(stdout);
      OutOffMem = 1;
#pragma omp flush(OutOffMem)
    }
  }                             /* end of parallel region */
  if (OutOffMem)
  {
    for(i=1;i<Npsi;i++)
    {
      ComplexWaveFunctionFreeValues(&(Opsi[i]));
    }
    ComplexWaveFunctionFree(&(Opsi[0]));
    return 1;
  }
  return 0;
}
uint32_t BlockOperatorPsiCC(const OperatorType O, WaveFunctionType const * const psi, WaveFunctionType * const Opsi, uint32_t const Npsi)
{
  uint32_t i, j, k, l, m, ii, ij, ik, OutOffMem;
  int    iOMP;
  double phase;
  if(Npsi==0) return 0;
  Opsi[0].QComplex = 1;
  Opsi[0].NFermion = psi[0].NFermion;
  Opsi[0].NBoson = psi[0].NBoson;
  Opsi[0].NBitsKey = psi[0].NBitsKey;
  if (WaveFunctionInit(&(Opsi[0])))
  {
    printf("Could not initialize Opsi[0] in BlockOperatorPsiCC\n");
    return 1;
  };
  for(i=1;i<Npsi;i++)
  {
    if(WaveFunctionInitCopyBasis(&(Opsi[i]), Opsi[0]))
    {
      printf("Could not initialize Opsi[%i] in BlockOperatorPsiCC\n",i);
      for(j=1;j<i;j++)
      {
        WaveFunctionFreeValues(&(Opsi[j]));
      }
      WaveFunctionFree(&(Opsi[0]));
      return 1;
    }
  }
  OutOffMem = 0;
#pragma omp parallel private(i,j,k,l,m,ii,ij,ik,phase)
  {
    uint8_t Det[Opsi->NWord];
    uint8_t MiniDet[MiniBlokSize * Opsi->NWord];
    double MiniValueRe[MiniBlokSize][Npsi];
    double MiniValueIm[MiniBlokSize][Npsi];
    uint32_t NMini;
    NMini = 0;
#pragma omp for
    for (iOMP = 1; iOMP <= psi[0].N; iOMP++)
    {
      if (!OutOffMem)
      {
        i = (uint32_t) iOMP;
        ii = i >> NBitsBlok;
        ik = i & BlokMask;
        ij = psi[0].NWord * ik;
        for (l = 0; l <= O.MaxLength; l++)
        {
          for (j = 0; j < O.l[l].N; j++)
          {
            for (k = 0; k < Opsi->NWord; k++) Det[k] = psi[0].Determinant[ii][ij + k];
            phase = OperatorLader(Det, &(O.l[l].o[j * l]), l, O.NFermion);
            if (fabs(phase) > (10 * DBL_MIN))
            {
              for(m=0;m<Npsi;m++)
              {
                MiniValueRe[NMini][m] = phase * (O.l[l].V[j].Re * psi[m].ValueRe[ii][ik] - O.l[l].V[j].Im * psi[m].ValueIm[ii][ik]);
                MiniValueIm[NMini][m] = phase * (O.l[l].V[j].Re * psi[m].ValueIm[ii][ik] + O.l[l].V[j].Im * psi[m].ValueRe[ii][ik]);
              }
              for (k = 0; k < Opsi->NWord; k++) MiniDet[NMini * Opsi->NWord + k] = Det[k];
              NMini++;
              if (NMini == MiniBlokSize)
              {                 /* flush Mini buffer */
                if (ComplexBlockWaveFunctionAddElementOMPMiniFlush(Opsi, &NMini, MiniDet, MiniValueRe[0], MiniValueIm[0], Npsi))
                {               /* Allocation failed */
                  printf("ComplexBlockWaveFunctionAddElementOMPMiniFlush failed in BlockOperatorPsiCC\n 01");
                  fflush(stdout);
                  OutOffMem = 1;
#pragma omp flush(OutOffMem)
                  j = O.l[l].N;
                  l = O.MaxLength;
                }
              }
            }
          }                     /* end of loop over different operators of lenght l */
        }                       /* end of loop over operators of different length */
      }                         /* Out Off Mem test */
    }                           /* end of loop over elements of psi this loop runs parallel */
    if (ComplexBlockWaveFunctionAddElementOMPMiniFlush(Opsi, &NMini, MiniDet, MiniValueRe[0], MiniValueIm[0], Npsi))
    {                           /* Allocation failed */
      printf("ComplexBlockWaveFunctionAddElementOMPMiniFlush failed in BlockOperatorPsiCC 02\n");
      fflush(stdout);
      OutOffMem = 1;
#pragma omp flush(OutOffMem)
    }
  }                             /* end of parallel region */
  if (OutOffMem)
  {
    for(i=1;i<Npsi;i++)
    {
      ComplexWaveFunctionFreeValues(&(Opsi[i]));
    }
    ComplexWaveFunctionFree(&(Opsi[0]));
    return 1;
  }
  return 0;
}
uint32_t BlockOperatorPsiMR(const OperatorType O, WaveFunctionType const * const psi, WaveFunctionType * const Opsi, uint32_t const Npsi)
{
  uint32_t i, j, k, l, m, ii, ij, ik, OutOffMem;
  int    iOMP;
  double phase;
  if(Npsi==0) return 0;
  Opsi[0].QComplex = 1;
  Opsi[0].NFermion = psi[0].NFermion;
  Opsi[0].NBoson = psi[0].NBoson;
  Opsi[0].NBitsKey = psi[0].NBitsKey;
  if (WaveFunctionInit(&(Opsi[0])))
  {
    printf("Could not initialize Opsi[0] in BlockOperatorPsiRR\n");
    return 1;
  };
  for(i=1;i<Npsi;i++)
  {
    if(WaveFunctionInitCopyBasis(&(Opsi[i]), Opsi[0]))
    {
      printf("Could not initialize Opsi[%i] in BlockOperatorPsiRR\n",i);
      for(j=1;j<i;j++)
      {
        WaveFunctionFreeValues(&(Opsi[j]));
      }
      WaveFunctionFree(&(Opsi[0]));
      return 1;
    }
  }
  OutOffMem = 0;
#pragma omp parallel private(i,j,k,l,m,ii,ij,ik,phase)
  {
    uint8_t Det[Opsi->NWord];
    uint8_t MiniDet[MiniBlokSize * Opsi->NWord];
    double MiniValueRe[MiniBlokSize][Npsi];
    double MiniValueIm[MiniBlokSize][Npsi];
    uint32_t NMini;
    NMini = 0;
#pragma omp for
    for (iOMP = 1; iOMP <= psi[0].N; iOMP++)
    {
      if (!OutOffMem)
      {
        i = (uint32_t) iOMP;
        ii = i >> NBitsBlok;
        ik = i & BlokMask;
        ij = psi[0].NWord * ik;
        for (l = 0; l <= O.MaxLength; l++)
        {
          if(O.l[l].QComplex==0)
          { /* operator real at this length */
            for (j = 0; j < O.l[l].N; j++)
            {
              for (k = 0; k < Opsi->NWord; k++) Det[k] = psi[0].Determinant[ii][ij + k];
              phase = OperatorLader(Det, &(O.l[l].o[j * l]), l, O.NFermion);
              if (fabs(phase) > (10 * DBL_MIN))
              {
                for(m=0;m<Npsi;m++)
                {
                  MiniValueRe[NMini][m] = phase * O.l[l].VRe[j] * psi[m].ValueRe[ii][ik];
                  MiniValueIm[NMini][m] = 0;
                }
                for (k = 0; k < Opsi->NWord; k++) MiniDet[NMini * Opsi->NWord + k] = Det[k];
                NMini++;
                if (NMini == MiniBlokSize)
                {                 /* flush Mini buffer */
                  if (ComplexBlockWaveFunctionAddElementOMPMiniFlush(Opsi, &NMini, MiniDet, MiniValueRe[0], MiniValueIm[0], Npsi))
                  {               /* Allocation failed */
                    printf("ComplexBlockWaveFunctionAddElementOMPMiniFlush failed in BlockOperatorPsiMR\n 01");
                    fflush(stdout);
                    OutOffMem = 1;
#pragma omp flush(OutOffMem)
                    j = O.l[l].N;
                    l = O.MaxLength;
                  }
                }
              }
            }                     /* end of loop over different operators of lenght l */
          }
          else
          { /* operator complex at this length */
            for (j = 0; j < O.l[l].N; j++)
            {
              for (k = 0; k < Opsi->NWord; k++) Det[k] = psi[0].Determinant[ii][ij + k];
              phase = OperatorLader(Det, &(O.l[l].o[j * l]), l, O.NFermion);
              if (fabs(phase) > (10 * DBL_MIN))
              {
                for(m=0;m<Npsi;m++)
                {
                  MiniValueRe[NMini][m] = phase * O.l[l].V[j].Re * psi[m].ValueRe[ii][ik];
                  MiniValueIm[NMini][m] = phase * O.l[l].V[j].Im * psi[m].ValueRe[ii][ik];
                }
                for (k = 0; k < Opsi->NWord; k++) MiniDet[NMini * Opsi->NWord + k] = Det[k];
                NMini++;
                if (NMini == MiniBlokSize)
                {                 /* flush Mini buffer */
                  if (ComplexBlockWaveFunctionAddElementOMPMiniFlush(Opsi, &NMini, MiniDet, MiniValueRe[0], MiniValueIm[0], Npsi))
                  {               /* Allocation failed */
                    printf("ComplexBlockWaveFunctionAddElementOMPMiniFlush failed in BlockOperatorPsiMR\n 01");
                    fflush(stdout);
                    OutOffMem = 1;
#pragma omp flush(OutOffMem)
                    j = O.l[l].N;
                    l = O.MaxLength;
                  }
                }
              }
            }                     /* end of loop over different operators of lenght l */
          }
        }                       /* end of loop over operators of different length */
      }                         /* Out Off Mem test */
    }                           /* end of loop over elements of psi this loop runs parallel */
    if (ComplexBlockWaveFunctionAddElementOMPMiniFlush(Opsi, &NMini, MiniDet, MiniValueRe[0], MiniValueIm[0], Npsi))
    {                           /* Allocation failed */
      printf("ComplexBlockWaveFunctionAddElementOMPMiniFlush failed in BlockOperatorPsiMR 02\n");
      fflush(stdout);
      OutOffMem = 1;
#pragma omp flush(OutOffMem)
    }
  }                             /* end of parallel region */
  if (OutOffMem)
  {
    for(i=1;i<Npsi;i++)
    {
      ComplexWaveFunctionFreeValues(&(Opsi[i]));
    }
    ComplexWaveFunctionFree(&(Opsi[0]));
    return 1;
  }
  return 0;
}
uint32_t BlockOperatorPsiMC(const OperatorType O, WaveFunctionType const * const psi, WaveFunctionType * const Opsi, uint32_t const Npsi)
{
  uint32_t i, j, k, l, m, ii, ij, ik, OutOffMem;
  int    iOMP;
  double phase;
  if(Npsi==0) return 0;
  Opsi[0].QComplex = 1;
  Opsi[0].NFermion = psi[0].NFermion;
  Opsi[0].NBoson = psi[0].NBoson;
  Opsi[0].NBitsKey = psi[0].NBitsKey;
  if (WaveFunctionInit(&(Opsi[0])))
  {
    printf("Could not initialize Opsi[0] in BlockOperatorPsiMC\n");
    return 1;
  };
  for(i=1;i<Npsi;i++)
  {
    if(WaveFunctionInitCopyBasis(&(Opsi[i]), Opsi[0]))
    {
      printf("Could not initialize Opsi[%i] in BlockOperatorPsiMC\n",i);
      for(j=1;j<i;j++)
      {
        WaveFunctionFreeValues(&(Opsi[j]));
      }
      WaveFunctionFree(&(Opsi[0]));
      return 1;
    }
  }
  OutOffMem = 0;
#pragma omp parallel private(i,j,k,l,m,ii,ij,ik,phase)
  {
    uint8_t Det[Opsi->NWord];
    uint8_t MiniDet[MiniBlokSize * Opsi->NWord];
    double MiniValueRe[MiniBlokSize][Npsi];
    double MiniValueIm[MiniBlokSize][Npsi];
    uint32_t NMini;
    NMini = 0;
#pragma omp for
    for (iOMP = 1; iOMP <= psi[0].N; iOMP++)
    {
      if (!OutOffMem)
      {
        i = (uint32_t) iOMP;
        ii = i >> NBitsBlok;
        ik = i & BlokMask;
        ij = psi[0].NWord * ik;
        for (l = 0; l <= O.MaxLength; l++)
        {
          if(O.l[l].QComplex==0)
          { // operator real at this length
            for (j = 0; j < O.l[l].N; j++)
            {
              for (k = 0; k < Opsi->NWord; k++) Det[k] = psi[0].Determinant[ii][ij + k];
              phase = OperatorLader(Det, &(O.l[l].o[j * l]), l, O.NFermion);
              if (fabs(phase) > (10 * DBL_MIN))
              {
                for(m=0;m<Npsi;m++)
                {
                  MiniValueRe[NMini][m] = phase * O.l[l].VRe[j] * psi[m].ValueRe[ii][ik];
                  MiniValueIm[NMini][m] = phase * O.l[l].VRe[j] * psi[m].ValueIm[ii][ik];
                }
                for (k = 0; k < Opsi->NWord; k++) MiniDet[NMini * Opsi->NWord + k] = Det[k];
                NMini++;
                if (NMini == MiniBlokSize)
                {                 /* flush Mini buffer */
                  if (ComplexBlockWaveFunctionAddElementOMPMiniFlush(Opsi, &NMini, MiniDet, MiniValueRe[0], MiniValueIm[0], Npsi))
                  {               /* Allocation failed */
                    printf("ComplexBlockWaveFunctionAddElementOMPMiniFlush failed in BlockOperatorPsiMC\n 01");
                    fflush(stdout);
                    OutOffMem = 1;
#pragma omp flush(OutOffMem)
                    j = O.l[l].N;
                    l = O.MaxLength;
                  }
                }
              }
            }                   /* end of loop over different operators of lenght l */
          }
          else
          { // operator complex at this length
            for (j = 0; j < O.l[l].N; j++)
            {
              for (k = 0; k < Opsi->NWord; k++) Det[k] = psi[0].Determinant[ii][ij + k];
              phase = OperatorLader(Det, &(O.l[l].o[j * l]), l, O.NFermion);
              if (fabs(phase) > (10 * DBL_MIN))
              {
                for(m=0;m<Npsi;m++)
                {
                  MiniValueRe[NMini][m] = phase * (O.l[l].V[j].Re * psi[m].ValueRe[ii][ik] - O.l[l].V[j].Im * psi[m].ValueIm[ii][ik] );
                  MiniValueIm[NMini][m] = phase * (O.l[l].V[j].Re * psi[m].ValueIm[ii][ik] + O.l[l].V[j].Im * psi[m].ValueRe[ii][ik] );
                }
                for (k = 0; k < Opsi->NWord; k++) MiniDet[NMini * Opsi->NWord + k] = Det[k];
                NMini++;
                if (NMini == MiniBlokSize)
                {                 /* flush Mini buffer */
                  if (ComplexBlockWaveFunctionAddElementOMPMiniFlush(Opsi, &NMini, MiniDet, MiniValueRe[0], MiniValueIm[0], Npsi))
                  {               /* Allocation failed */
                    printf("ComplexBlockWaveFunctionAddElementOMPMiniFlush failed in BlockOperatorPsiMC\n 01");
                    fflush(stdout);
                    OutOffMem = 1;
#pragma omp flush(OutOffMem)
                    j = O.l[l].N;
                    l = O.MaxLength;
                  }
                }
              }
            }                   /* end of loop over different operators of lenght l */
          }
        }                       /* end of loop over operators of different length */
      }                         /* Out Off Mem test */
    }                           /* end of loop over elements of psi this loop runs parallel */
    if (ComplexBlockWaveFunctionAddElementOMPMiniFlush(Opsi, &NMini, MiniDet, MiniValueRe[0], MiniValueIm[0], Npsi))
    {                           /* Allocation failed */
      printf("ComplexBlockWaveFunctionAddElementOMPMiniFlush failed in BlockOperatorPsiMC 02\n");
      fflush(stdout);
      OutOffMem = 1;
#pragma omp flush(OutOffMem)
    }
  }                             /* end of parallel region */
  if (OutOffMem)
  {
    for(i=1;i<Npsi;i++)
    {
      ComplexWaveFunctionFreeValues(&(Opsi[i]));
    }
    ComplexWaveFunctionFree(&(Opsi[0]));
    return 1;
  }
  return 0;
}
uint32_t BlockOperatorPsi(const OperatorType O, WaveFunctionType const * const psi, WaveFunctionType * const Opsi, uint32_t const Npsi)
{
  if(Npsi==0) return 0;
  if(O.Res.N!=0) return BlockOperatorPsiRestricted(O, psi, Opsi, O.Res, Npsi);
  if (O.QComplex == 0)
  {
    if (psi[0].QComplex)
    {
      return BlockOperatorPsiRC(O, psi, Opsi, Npsi);
    }
    return BlockOperatorPsiRR(O, psi, Opsi, Npsi);
  }
  if (O.QComplex == 1)
  {
    if (psi[0].QComplex)
    {
      return BlockOperatorPsiCC(O, psi, Opsi, Npsi);
    }
    return BlockOperatorPsiCR(O, psi, Opsi, Npsi);
  }
  if (psi[0].QComplex)
  {
    return BlockOperatorPsiMC(O, psi, Opsi, Npsi);
  }
  return BlockOperatorPsiMR(O, psi, Opsi, Npsi);
}


/*  Remove determinants from memory for which all blocks have a prefactor smaller than epsilon
 *  set prefactor of determinats to zero that are smaller than epsilon
 */
void RealBlockWaveFunctionTruncateBasis(WaveFunctionType *const psi, uint32_t const Npsi, const double epsilon)
{
  uint32_t startN, iBlok, iElement, ipsi, KeepDeterminant;
  if(Npsi==0) return;
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
  psi[0].Next = (uint32_t **)realloc(psi[0].Next, psi[0].NBlok * sizeof(uint32_t *));
  psi[0].Determinant = (uint8_t **)realloc(psi[0].Determinant, psi[0].NBlok * sizeof(uint8_t *));
  for(ipsi=0;ipsi<Npsi;ipsi++)
  {
    psi[ipsi].ValueRe = (double **)realloc(psi[ipsi].ValueRe, psi[0].NBlok * sizeof(double *));
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
  if(Npsi==0) return;
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
  psi[0].Next = (uint32_t **)realloc(psi[0].Next, psi[0].NBlok * sizeof(uint32_t *));
  psi[0].Determinant = (uint8_t **)realloc(psi[0].Determinant, psi[0].NBlok * sizeof(uint8_t *));
  for(ipsi=0;ipsi<Npsi;ipsi++)
  {
    psi[ipsi].ValueRe = (double **)realloc(psi[ipsi].ValueRe, psi[0].NBlok * sizeof(double *));
    psi[ipsi].ValueIm = (double **)realloc(psi[ipsi].ValueIm, psi[0].NBlok * sizeof(double *));
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
  if(Npsi==0) return;
  if(psi[0].QComplex)
  {
    ComplexBlockWaveFunctionTruncateBasis(psi, Npsi, epsilon);
  }
  else
  {
    RealBlockWaveFunctionTruncateBasis(psi, Npsi, epsilon);
  }
}

uint32_t RealBlockWaveFunctionResetKeylist(WaveFunctionType * const a, const uint32_t NBitsKey, uint32_t const Npsi)
{
  uint32_t *tempKeyList,i,iBlok,iElement;
  const uint32_t startN=a[0].N;
  double aval[Npsi];
  // wipeout old keylist and replace it by a new one
  if((tempKeyList= (uint32_t*)calloc((1<<NBitsKey), sizeof(uint32_t)))==NULL)
  {
    printf("realloc failed in RealWaveFunctionResetKeylist\n");
    return 1;
  }
  free(a[0].KeyList);
  for(i=0;i<Npsi;i++)
  {
    a[i].ValueRe[0][0] = 0;
    a[i].NBitsKey=NBitsKey;
    a[i].N = 0;
    a[i].KeyList=tempKeyList;
  }
  // add all values back (and rest a.next while doing so)
  // first loop over all Blokcs
  for (iBlok = 0; iBlok < a[0].NBlok - 1; iBlok++)
  {
    for (iElement = 0; iElement < BlokSize; iElement++)
    {
      a[0].Next[iBlok][iElement] = 0;
      if((iBlok==0)&&(iElement==0))continue;
      for(i=0;i<Npsi;i++)
      {
        aval[i] = a[i].ValueRe[iBlok][iElement];
      }
      if (RealBlockWaveFunctionAddElement(a, &(a[0].Determinant[iBlok][iElement * a[0].NWord]), aval, Npsi)!= 0) AllocFailed();
    }
  }
  // for the last block loop over all elements
  iBlok = a[0].NBlok - 1;
  for (iElement = 0; iElement < 1 + startN - iBlok * BlokSize; ++iElement)
  {
    a[0].Next[iBlok][iElement] = 0;
    if((iBlok==0)&&(iElement==0))continue;
    for(i=0;i<Npsi;i++)
    {
      aval[i] = a[i].ValueRe[iBlok][iElement];
    }
    if (RealBlockWaveFunctionAddElement(a, &(a[0].Determinant[iBlok][iElement * a[0].NWord]), aval, Npsi)!= 0) AllocFailed();
  }
  return 0;
}
uint32_t ComplexBlockWaveFunctionResetKeylist(WaveFunctionType * const a, const uint32_t NBitsKey, uint32_t const Npsi)
{
  uint32_t *tempKeyList,i,iBlok,iElement;
  const uint32_t startN=a[0].N;
  double avalRe[Npsi],avalIm[Npsi];
  // wipeout old keylist and replace it by a new one
  if((tempKeyList= (uint32_t*)calloc((1<<NBitsKey), sizeof(uint32_t)))==NULL)
  {
    printf("realloc failed in RealWaveFunctionResetKeylist\n");
    return 1;
  }
  free(a[0].KeyList);
  for(i=0;i<Npsi;i++)
  {
    a[i].ValueRe[0][0] = 0;
    a[i].ValueIm[0][0] = 0;
    a[i].NBitsKey=NBitsKey;
    a[i].N = 0;
    a[i].KeyList=tempKeyList;
  }
  // add all values back (and rest a.next while doing so)
  // first loop over all Blokcs
  for (iBlok = 0; iBlok < a[0].NBlok - 1; iBlok++)
  {
    for (iElement = 0; iElement < BlokSize; iElement++)
    {
      a[0].Next[iBlok][iElement] = 0;
      if((iBlok==0)&&(iElement==0))continue;
      for(i=0;i<Npsi;i++)
      {
        avalRe[i] = a[i].ValueRe[iBlok][iElement];
        avalIm[i] = a[i].ValueRe[iBlok][iElement];
      }
      if (ComplexBlockWaveFunctionAddElement(a, &(a[0].Determinant[iBlok][iElement * a[0].NWord]), avalRe, avalIm, Npsi)!= 0) AllocFailed();
    }
  }
  // for the last block loop over all elements
  iBlok = a[0].NBlok - 1;
  for (iElement = 0; iElement < 1 + startN - iBlok * BlokSize; ++iElement)
  {
    a[0].Next[iBlok][iElement] = 0;
    if((iBlok==0)&&(iElement==0))continue;
    for(i=0;i<Npsi;i++)
    {
      avalRe[i] = a[i].ValueRe[iBlok][iElement];
      avalIm[i] = a[i].ValueIm[iBlok][iElement];
    }
    if (ComplexBlockWaveFunctionAddElement(a, &(a[0].Determinant[iBlok][iElement * a[0].NWord]), avalRe, avalIm, Npsi)!= 0) AllocFailed();
  }
  return 0;
}

uint32_t BlockWaveFunctionResetKeylist(WaveFunctionType * const a, const uint32_t NBitsKey, uint32_t const Npsi)
{
  if(Npsi==0) return 0;
  if (a->QComplex == 0)
  {
    return RealBlockWaveFunctionResetKeylist(a, NBitsKey, Npsi);
  }
  else
  {
    return ComplexBlockWaveFunctionResetKeylist(a, NBitsKey, Npsi);
  }
}





/*  calculate S[i][j] = < psi[i] | psi[j] >
 *  S not initialized on input
 */
uint32_t BlockWaveFunctionCalculatOverlapMatrixSerial(WaveFunctionType *const psi, uint32_t const Npsi, CompactMatrixType *const S)
{
  uint32_t i,j;
  S->Ni=Npsi;
  S->Nj=Npsi;
  if(Npsi==0)
  {
    S->QComplex=0;
  }
  else
  {
    S->QComplex=psi[0].QComplex;
  }
  InitCompactMatrix(S);
  if(S->QComplex)
  {
    for(i=0;i<Npsi;i++)
    {
      for(j=0;j<Npsi;j++)
      {
        ComplexWaveFunctionDotProductSameBasis(psi[i], psi[j], &(S->V[j+i*Npsi].Re), &(S->V[j+i*Npsi].Im));
      }
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

/*
 *   calculate S[i,j] = < a[i] | b[j] > with a and b lists on the same basis, but a and b on a different basis
 *   S should be initiallized before the function call
 */
uint32_t RealBlockWaveFunctionDotProduct( const WaveFunctionType * a, const uint32_t Npsia, const WaveFunctionType * b, const uint32_t Npsib, CompactMatrixType *S)
{
  uint32_t failed=0;
  if((Npsia==0)||(Npsib==0)) return 0;
  /* switch coefficients to reduce number of hash lookups */
  if (a[0].N > b[0].N)
  {
    CompactMatrixType ST;
    uint32_t i,j;
    ST.Ni=S->Nj;
    ST.Nj=S->Ni;
    ST.QComplex=S->QComplex;
    InitCompactMatrix(&ST);
    if(RealBlockWaveFunctionDotProduct(b, Npsib, a, Npsia, &ST))
    {
      printf("RealBlockWaveFunctionDotProduct failed in RealBlockWaveFunctionDotProduct\n");
      return 1;
    }
    for(i=0;i<S->Ni;i++)
    {
      for(j=0;j<S->Nj;j++)
      {
        S->VRe[i*S->Nj+j] = ST.VRe[j*ST.Nj+i];
      }
    }
    FreeCompactMatrix(&ST);
    return 0;
  }
  CompactMatrixScale(S, 0);
#pragma omp parallel
  {
    uint32_t i, j, k, index, NextIndex, PsiAIsPsiB;
    int    iOMP;
    volatile double t,y;
    CompactMatrixType Sprivate,c;
    Sprivate.Ni=S->Ni;
    Sprivate.Nj=S->Nj;
    Sprivate.QComplex=0;
    if(InitCompactMatrix(&Sprivate)) failed =1;
    c.Ni=S->Ni;
    c.Nj=S->Nj;
    c.QComplex=0;
    if(InitCompactMatrix(&c)) failed =1;
    if(failed==0)
    {
#pragma omp for
      for (iOMP = 1; iOMP <= (int)(a[0].N); ++iOMP)
      {
        k = (uint32_t) iOMP;
        NextIndex = b[0].KeyList[hash(&a[0].Determinant[BlokElementDet(k, a[0].NWord, 0)], b[0].NWord, b[0].NBitsKey)];
        while (NextIndex)
        {
          PsiAIsPsiB = 1;
          for (j = 0; j < b[0].NWord; j++)
          {
            if (b[0].Determinant[BlokElementDet(NextIndex, b[0].NWord, j)] != a[0].Determinant[BlokElementDet(k, a[0].NWord, j)])
            {
              PsiAIsPsiB = 0;
            }
          }
          if (PsiAIsPsiB)
          {                         /* a=b found add the product to Dot and leave loop  */
            for(i=0;i<Npsia;i++)
            {
              index = i * S->Ni;
              for(j=0;j<Npsib;j++)
              {
                y = a[i].ValueRe[BlokElement(k)] * b[j].ValueRe[BlokElement(NextIndex)] - c.VRe[index];
                t = Sprivate.VRe[index] + y;
                c.VRe[index] = t - Sprivate.VRe[index];
                c.VRe[index] -= y;
                Sprivate.VRe[index] = t;
                index++;
              }
            }
            NextIndex = 0;
          }
          else
          {
            NextIndex = b[0].Next[BlokElement(NextIndex)];
          }
        }
      }
#pragma omp critical
      for(i=0;i<Sqr(Sprivate.Ni);i++)
      {
        S->VRe[i]+=Sprivate.VRe[i];
      }
      FreeCompactMatrix(&c);
      FreeCompactMatrix(&Sprivate);
    }
  }
  return failed;
}
uint32_t ComplexBlockWaveFunctionDotProduct( const WaveFunctionType * a, const uint32_t Npsia, const WaveFunctionType * b, const uint32_t Npsib, CompactMatrixType *S)
{
  uint32_t failed=0;
  if((Npsia==0)||(Npsib==0)) return 0;
  /* switch coefficients to reduce number of hash lookups */
  if (a[0].N > b[0].N)
  {
    CompactMatrixType ST;
    uint32_t i,j;
    ST.Ni=S->Nj;
    ST.Nj=S->Ni;
    ST.QComplex=S->QComplex;
    InitCompactMatrix(&ST);
    if(ComplexBlockWaveFunctionDotProduct(b, Npsib, a, Npsia, &ST))
    {
      printf("ComplexBlockWaveFunctionDotProduct failed in ComplexBlockWaveFunctionDotProduct\n");
      return 1;
    }
    for(i=0;i<S->Ni;i++)
    {
      for(j=0;j<S->Nj;j++)
      {
        S->V[i*S->Nj+j].Re = ST.V[j*ST.Nj+i].Re;
        S->V[i*S->Nj+j].Im =-ST.V[j*ST.Nj+i].Im;
      }
    }
    FreeCompactMatrix(&ST);
    return 0;
  }
  CompactMatrixScale(S, 0);
#pragma omp parallel
  {
    uint32_t i, j, k, index, NextIndex, PsiAIsPsiB;
    int    iOMP;
    volatile Complex16 t,y;
    CompactMatrixType Sprivate,c;
    Sprivate.Ni=S->Ni;
    Sprivate.Nj=S->Nj;
    Sprivate.QComplex=1;
    if(InitCompactMatrix(&Sprivate)) failed =1;
    c.Ni=S->Ni;
    c.Nj=S->Nj;
    c.QComplex=1;
    if(InitCompactMatrix(&c)) failed =1;
    if(failed==0)
    {
#pragma omp for
      for (iOMP = 1; iOMP <= (int)(a[0].N); ++iOMP)
      {
        k = (uint32_t) iOMP;
        NextIndex = b[0].KeyList[hash(&a[0].Determinant[BlokElementDet(k, a[0].NWord, 0)], b[0].NWord, b[0].NBitsKey)];
        while (NextIndex)
        {
          PsiAIsPsiB = 1;
          for (j = 0; j < b[0].NWord; j++)
          {
            if (b[0].Determinant[BlokElementDet(NextIndex, b[0].NWord, j)] != a[0].Determinant[BlokElementDet(k, a[0].NWord, j)])
            {
              PsiAIsPsiB = 0;
            }
          }
          if (PsiAIsPsiB)
          {                         /* a=b found add the product to Dot and leave loop  */
            for(i=0;i<Npsia;i++)
            {
              index = i * S->Ni;
              for(j=0;j<Npsib;j++)
              {
                y.Re = a[i].ValueRe[BlokElement(k)] * b[j].ValueRe[BlokElement(NextIndex)] + a[i].ValueIm[BlokElement(k)] * b[j].ValueIm[BlokElement(NextIndex)]  - c.V[index].Re;
                y.Im = a[i].ValueRe[BlokElement(k)] * b[j].ValueIm[BlokElement(NextIndex)] - a[i].ValueIm[BlokElement(k)] * b[j].ValueRe[BlokElement(NextIndex)]  - c.V[index].Im;
                t.Re = Sprivate.V[index].Re + y.Re;
                t.Im = Sprivate.V[index].Im + y.Im;
                c.V[index].Re = t.Re - Sprivate.V[index].Re;
                c.V[index].Im = t.Im - Sprivate.V[index].Im;
                c.V[index].Re -= y.Re;
                c.V[index].Im -= y.Im;
                Sprivate.V[index].Re = t.Re;
                Sprivate.V[index].Im = t.Im;
                index++;
              }
            }
            NextIndex = 0;
          }
          else
          {
            NextIndex = b[0].Next[BlokElement(NextIndex)];
          }
        }
      }
#pragma omp critical
      for(i=0;i<Sqr(Sprivate.Ni);i++)
      {
        S->V[i].Re+=Sprivate.V[i].Re;
        S->V[i].Im+=Sprivate.V[i].Im;
      }
      FreeCompactMatrix(&c);
      FreeCompactMatrix(&Sprivate);
    }
  }
  return failed;
}
uint32_t BlockWaveFunctionDotProduct( const WaveFunctionType * a, const uint32_t Npsia, const WaveFunctionType * b, const uint32_t Npsib, CompactMatrixType *S)
{
  if((Npsia==0)||(Npsib==0)) return 0;
  if(a[0].QComplex==0)
  {
    if(b[0].QComplex==0)
    {
      if(S->QComplex==0)
      {
        // all real
        return RealBlockWaveFunctionDotProduct(a, Npsia, b, Npsib, S);
      }
      else
      {
        // psia and psib real S complex
        CompactMatrixType Sreal;
        uint32_t i;
        Sreal = *S;
        Sreal.QComplex=0;
        Sreal.VRe=(double*)calloc(Sreal.Ni*Sreal.Nj,sizeof(double));
        if(Sreal.VRe==NULL)
        {
          printf("malloc failed in BlockWaveFunctionDotProduct for real allocation of matrix S\n");
          fflush(stdout);
          return 1;
        }
        if(RealBlockWaveFunctionDotProduct(a, Npsia, b, Npsib, &Sreal))
        {
          printf("RealBlockWaveFunctionDotProduct failed in BlockWaveFunctionDotProduct\n");
          FreeCompactMatrix(&Sreal);
          return 1;
        }
        for(i=0;i<Sreal.Ni*Sreal.Nj;i++)
        {
          S->V[i].Re=Sreal.VRe[i];
          S->V[i].Im=0;
        }
        FreeCompactMatrix(&Sreal);
        return 0;
      }
    }
    else
    {
      // psia real psib complex
      //TODO: real.complex dot
      printf("Dot product of real times complex Block function not implemented\n");
      return 1;
    }
  }
  else
  {
    if(b[0].QComplex==0)
    {
      // psia complex psib real
      //TODO: complex.real dot
      printf("Dot product of complex times real Block function not implemented\n");
      return 1;
    }
    else
    {
      if(S->QComplex==0)
      {
        // a and b complex S real, i.e. only save the real part of the dot product
        VerbosityPrintf(0xFFFF, "Warning in BlockWaveFunctionDotProduct\nMultiplication of two complex Block Wave-functions (%s and %s) stored as a real dot-product\n",a[0].Name,b[0].Name);
        CompactMatrixType Sim;
        uint32_t i;
        Sim=*S;
        Sim.QComplex=1;
        Sim.V=(Complex16*)calloc(Sim.Ni*Sim.Nj,sizeof(Complex16));
        if(Sim.V==NULL)
        {
          printf("malloc failed in BlockWaveFunctionDotProduct for complex allocation of matrix S\n");
          fflush(stdout);
          return 1;
        }
        if(ComplexBlockWaveFunctionDotProduct(a, Npsia, b, Npsib, &Sim))
        {
          printf("ComplexBlockWaveFunctionDotProduct failed in BlockWaveFunctionDotProduct\n");
          FreeCompactMatrix(&Sim);
          return 1;
        }
        for(i=0;i<Sim.Ni*Sim.Nj;i++)
        {
          S->V[i].Re=Sim.VRe[i];
          S->V[i].Im=0;
        }
        FreeCompactMatrix(&Sim);
        return 0;
      }
      else
      {
        // all complex
        return ComplexBlockWaveFunctionDotProduct(a, Npsia, b, Npsib, S);
      }
    }
  }
}



