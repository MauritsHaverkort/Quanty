
/*
*   Input pointer to a determinant and an integer that specifies the position that has to be created
*   Output integer containing the phase (-1,1, or 0)
*   Effects on the determinant pointed to
*     If element (bit) i of the determinant was zero it now will be one
*/
double inline OperatorCr(uint8_t * const Determinant, const uint16_t i, const uint32_t NFermion)
{
  uint8_t detblok;
  uint32_t j;
  if (i >= NFermion)
  {
    j = i - NFermion + ((NFermion + 7) >> 3);
    if (Determinant[j] == 255)
    {
      return 0;
    }
    else
    {
      return sqrt(++Determinant[j]);  /* add one to Determinant[j] and then return sqrt(Determinant[j]) */
    }
  }
  if (Determinant[i >> 3] & (128 >> (i & 7)))
  {
    return 0;
  }
  Determinant[i >> 3] += (128 >> (i & 7));
  detblok = Determinant[i >> 3] >> (8 - (i & 7));
  for (j = i >> 3; j != 0;)
  {
    j--;
    detblok ^= Determinant[j];
  }
  detblok ^= detblok >> 4;
  detblok &= 15;
  detblok ^= detblok >> 2;
  detblok &= 3;
  detblok ^= detblok >> 1;
  detblok &= 1;
  if (detblok)
  {
    return -1;
  }
  return 1;
}

/*
*   Input pointer to a determinant and an integer that specifies the position that has to be annihilated
*   Output integer containing the phase (-1,1, or 0)
*   Effects on the determinant pointed to
*     If element (bit) i of the determinant was one it now will be zero
*/
double inline OperatorAn(uint8_t * const Determinant, const uint16_t i, const uint32_t NFermion)
{
  uint8_t detblok;
  uint32_t j;
  if (i >= NFermion)
  {
    j = i - NFermion + ((NFermion + 7) >> 3);
    if (Determinant[j] == 0)
    {
      return 0;
    }
    else
    {
      return sqrt(Determinant[j]--);  /* return sqrt(Determinant[j]) and subtract one from Determinant[j] */
    }
  }
  if (!(Determinant[i >> 3] & (128 >> (i & 7))))
  {
    return 0;
  }
  Determinant[i >> 3] -= (128 >> (i & 7));
  detblok = Determinant[i >> 3] >> (8 - (i & 7));
  for (j = i >> 3; j != 0;)
  {
    j--;
    detblok ^= Determinant[j];
  }
  detblok ^= detblok >> 4;
  detblok &= 15;
  detblok ^= detblok >> 2;
  detblok &= 3;
  detblok ^= detblok >> 1;
  detblok &= 1;
  if (detblok)
  {
    return -1;
  }
  return 1;
}

/*
 *  Input a pointer to a determinant, a pointer to a list of positions that have to be created or annihilated, the lenght of the list and the number of Ferminons
 *  and return the phase and new Determinant
 *
 *  Slight warning: there is no check on the maximum value of o[i]; if one tries to annihilate a particle that is not in the basis set then a bit somewhere someplace
 *  in memory will be set..... o[i] should be smaller then NFermion+NBoson-1 for all i.
 */
double inline OperatorLader(uint8_t * const Determinant, const uint16_t * const o, const uint32_t lenght,
                     const uint32_t NFermion)
{
  double phase = 1;
  uint32_t i;
  for (i = 0; i < lenght; i++)
  {
    if (o[i] & CrAnMask)
    {
      phase *= OperatorCr(Determinant, (o[i] & NotCrAnMask), NFermion);
    }
    else
    {
      phase *= OperatorAn(Determinant, o[i], NFermion);
    }
  }
  return phase;
}


uint32_t inline PsiAddOperatorPsiConserveBasisCC(const OperatorType *O, const WaveFunctionType *psi,
                                          WaveFunctionType * const Opsi)
{
  return 0;
}

uint32_t inline PsiAddOperatorPsiConserveBasisCR(const OperatorType *O, const WaveFunctionType * psi,
                                          WaveFunctionType * const Opsi)
{
  return 0;
}

uint32_t inline PsiAddOperatorPsiConserveBasisRC(const OperatorType *O, const WaveFunctionType *psi,
                                          WaveFunctionType * const Opsi)
{
  return 0;
}

uint32_t inline PsiAddOperatorPsiConserveBasisRR(const OperatorType *O, const WaveFunctionType *psi,
                                          WaveFunctionType * const Opsi)
{
  return 0;
}

uint32_t inline PsiAddOperatorPsiConserveBasisRRC(const OperatorType *O, const WaveFunctionType *psi,
                                           WaveFunctionType * const Opsi)
{
  return 0;
}

uint32_t inline PsiAddOperatorPsiConserveBasisMC(const OperatorType *O, const WaveFunctionType *psi,
                                          WaveFunctionType * const Opsi)
{
  return 0;
}

uint32_t inline PsiAddOperatorPsiConserveBasisMR(const OperatorType *O, const WaveFunctionType *psi,
                                          WaveFunctionType * const Opsi)
{
  return 0;
}

uint32_t OperatorPsiConserveBasisCC(const OperatorType O, const WaveFunctionType psi, WaveFunctionType * const Opsi)
{
  uint32_t i, ii, ij, ik, j, k, l;
  int    iOMP;
  double phase;
  uint8_t *Det;
  Opsi->N = psi.N;
  Opsi->NBlok = psi.NBlok;
  Opsi->NFermion = psi.NFermion;
  Opsi->NBoson = psi.NBoson;
  Opsi->NWord = psi.NWord;
  Opsi->Next = psi.Next;
  Opsi->Determinant = psi.Determinant;
  Opsi->NBitsKey = psi.NBitsKey;
  Opsi->KeyList = psi.KeyList;
  Opsi->QComplex = psi.QComplex;
  Opsi->ValueRe = malloc(Opsi->NBlok * sizeof(double *));
  if (Opsi->ValueRe == NULL)
  {                             /* Allocation Failed */
    printf("alloc failed OperatorPsiConserveBasis 01Re\n");
    fflush(stdout);
    Opsi->N = 0;
    Opsi->NBlok = 0;
    return 1;
  }
  Opsi->ValueIm = malloc(Opsi->NBlok * sizeof(double *));
  if (Opsi->ValueIm == NULL)
  {                             /* Allocation Failed */
    printf("alloc failed OperatorPsiConserveBasis 01Im\n");
    fflush(stdout);
    Opsi->N = 0;
    Opsi->NBlok = 0;
    free(Opsi->ValueRe);
    return 1;
  }
  for (i = 0; i < Opsi->NBlok; i++)
  {
    Opsi->ValueRe[i] = calloc(BlokSize, sizeof(double));
    Opsi->ValueIm[i] = calloc(BlokSize, sizeof(double));
    if (Opsi->ValueIm[i] == NULL)
    {                           /* Allocation Failed */
      printf("alloc failed OperatorPsiConserveBasis 02\n");
      fflush(stdout);
      if (Opsi->ValueRe[i] != NULL)
        free(Opsi->ValueRe[i]);
      for (j = 0; j < i; j++)
      {
        free(Opsi->ValueRe[j]);
        free(Opsi->ValueIm[j]);
      }
      Opsi->N = 0;
      Opsi->NBlok = 0;
      free(Opsi->ValueRe);
      free(Opsi->ValueIm);
      return 1;
    }
  }
#pragma omp parallel private(Det,i,j,k,l,ii,ij,ik,phase,iOMP)
  {
    Det = malloc(Opsi->NWord * sizeof(uint8_t));
    if (Det == NULL)
      AllocFailed();
#pragma omp for
    for (iOMP = 1; iOMP <= (int)psi.N; iOMP++)
    {
      i = (uint32_t) iOMP;
      ii = i >> NBitsBlok;
      ik = i & BlokMask;
      ij = psi.NWord * ik;
      for (l = 0; l <= O.MaxLength; l++)
      {                         /* act with a lader operator of length l on Determinant */
        for (j = 0; j < O.l[l].N; j++)
        {                       /* take opperator number k of length l and act on Determinant */
          for (k = 0; k < Opsi->NWord; k++)
            Det[k] = psi.Determinant[ii][ij + k];
          phase = OperatorLader(Det, &(O.l[l].o[j * l]), l, O.NFermion);
          ComplexWaveFunctionAddElementConserveBasis(*Opsi, Det,
                                                     phase * (O.l[l].V[j].Re * psi.ValueRe[ii][ik] -
                                                              O.l[l].V[j].Im * psi.ValueIm[ii][ik]),
                                                     phase * (O.l[l].V[j].Re * psi.ValueIm[ii][ik] +
                                                              O.l[l].V[j].Im * psi.ValueRe[ii][ik]));
        }
      }
    }
    free(Det);
  }
  return 0;
}

uint32_t OperatorPsiConserveBasisCR(const OperatorType O, const WaveFunctionType psi, WaveFunctionType * const Opsi)
{
  uint32_t i, ii, ij, ik, j, k, l;
  int    iOMP;
  double phase;
  uint8_t *Det;
  Opsi->N = psi.N;
  Opsi->NBlok = psi.NBlok;
  Opsi->NFermion = psi.NFermion;
  Opsi->NBoson = psi.NBoson;
  Opsi->NWord = psi.NWord;
  Opsi->Next = psi.Next;
  Opsi->Determinant = psi.Determinant;
  Opsi->NBitsKey = psi.NBitsKey;
  Opsi->KeyList = psi.KeyList;
  Opsi->QComplex = 1;
  Opsi->ValueRe = malloc(Opsi->NBlok * sizeof(double *));
  if (Opsi->ValueRe == NULL)
  {                             /* Allocation Failed */
    printf("alloc failed OperatorPsiConserveBasis 01Re\n");
    fflush(stdout);
    Opsi->N = 0;
    Opsi->NBlok = 0;
    return 1;
  }
  Opsi->ValueIm = malloc(Opsi->NBlok * sizeof(double *));
  {                             /* Allocation Failed */
    printf("alloc failed OperatorPsiConserveBasis 01Im\n");
    fflush(stdout);
    Opsi->N = 0;
    Opsi->NBlok = 0;
    free(Opsi->ValueRe);
    return 1;
  }
  for (i = 0; i < Opsi->NBlok; i++)
  {
    Opsi->ValueRe[i] = calloc(BlokSize, sizeof(double));
    Opsi->ValueIm[i] = calloc(BlokSize, sizeof(double));
    if (Opsi->ValueIm[i] == NULL)
    {                           /* Allocation Failed */
      printf("alloc failed OperatorPsiConserveBasis 02\n");
      fflush(stdout);
      if (Opsi->ValueRe[i] != NULL)
        free(Opsi->ValueRe[i]);
      for (j = 0; j < i; j++)
      {
        free(Opsi->ValueRe[j]);
        free(Opsi->ValueIm[j]);
      }
      Opsi->N = 0;
      Opsi->NBlok = 0;
      free(Opsi->ValueRe);
      free(Opsi->ValueIm);
      return 1;
    }
  }
#pragma omp parallel private(Det,i,j,k,l,ii,ij,ik,phase,iOMP)
  {
    Det = malloc(Opsi->NWord * sizeof(uint8_t));
    if (Det == NULL)
      AllocFailed();
#pragma omp for
    for (iOMP = 1; iOMP <= (int)psi.N; iOMP++)
    {
      i = (uint32_t) iOMP;
      ii = i >> NBitsBlok;
      ik = i & BlokMask;
      ij = psi.NWord * ik;
      for (l = 0; l <= O.MaxLength; l++)
      {                         /* act with a lader operator of length l on Determinant */
        for (j = 0; j < O.l[l].N; j++)
        {                       /* take opperator number k of length l and act on Determinant */
          for (k = 0; k < Opsi->NWord; k++)
            Det[k] = psi.Determinant[ii][ij + k];
          phase = OperatorLader(Det, &(O.l[l].o[j * l]), l, O.NFermion);
          ComplexWaveFunctionAddElementConserveBasis(*Opsi, Det, phase * (O.l[l].V[j].Re * psi.ValueRe[ii][ik]),
                                                     phase * (O.l[l].V[j].Im * psi.ValueRe[ii][ik]));
        }
      }
    }
    free(Det);
  }
  return 0;
}

uint32_t OperatorPsiConserveBasisRC(const OperatorType O, const WaveFunctionType psi, WaveFunctionType * const Opsi)
{
  uint32_t i, ii, ij, ik, j, k, l;
  int    iOMP;
  double phase;
  uint8_t *Det;
  Opsi->N = psi.N;
  Opsi->NBlok = psi.NBlok;
  Opsi->NFermion = psi.NFermion;
  Opsi->NBoson = psi.NBoson;
  Opsi->NWord = psi.NWord;
  Opsi->Next = psi.Next;
  Opsi->Determinant = psi.Determinant;
  Opsi->NBitsKey = psi.NBitsKey;
  Opsi->KeyList = psi.KeyList;
  Opsi->QComplex = psi.QComplex;
  Opsi->ValueRe = malloc(Opsi->NBlok * sizeof(double *));
  if (Opsi->ValueRe == NULL)
  {                             /* Allocation Failed */
    printf("alloc failed OperatorPsiConserveBasis 01Re\n");
    fflush(stdout);
    Opsi->N = 0;
    Opsi->NBlok = 0;
    return 1;
  }
  Opsi->ValueIm = malloc(Opsi->NBlok * sizeof(double *));
  {                             /* Allocation Failed */
    printf("alloc failed OperatorPsiConserveBasis 01Im\n");
    fflush(stdout);
    Opsi->N = 0;
    Opsi->NBlok = 0;
    free(Opsi->ValueRe);
    return 1;
  }
  for (i = 0; i < Opsi->NBlok; i++)
  {
    Opsi->ValueRe[i] = calloc(BlokSize, sizeof(double));
    Opsi->ValueIm[i] = calloc(BlokSize, sizeof(double));
    if (Opsi->ValueIm[i] == NULL)
    {                           /* Allocation Failed */
      printf("alloc failed OperatorPsiConserveBasis 02\n");
      fflush(stdout);
      if (Opsi->ValueRe[i] != NULL)
        free(Opsi->ValueRe[i]);
      for (j = 0; j < i; j++)
      {
        free(Opsi->ValueRe[j]);
        free(Opsi->ValueIm[j]);
      }
      Opsi->N = 0;
      Opsi->NBlok = 0;
      free(Opsi->ValueRe);
      free(Opsi->ValueIm);
      return 1;
    }
  }
#pragma omp parallel private(Det,i,j,k,l,ii,ij,ik,phase,iOMP)
  {
    Det = malloc(Opsi->NWord * sizeof(uint8_t));
    if (Det == NULL)
      AllocFailed();
#pragma omp for
    for (iOMP = 1; iOMP <= (int)psi.N; iOMP++)
    {
      i = (uint32_t) iOMP;
      ii = i >> NBitsBlok;
      ik = i & BlokMask;
      ij = psi.NWord * ik;
      for (l = 0; l <= O.MaxLength; l++)
      {                         /* act with a lader operator of length l on Determinant */
        for (j = 0; j < O.l[l].N; j++)
        {                       /* take opperator number k of length l and act on Determinant */
          for (k = 0; k < Opsi->NWord; k++)
            Det[k] = psi.Determinant[ii][ij + k];
          phase = OperatorLader(Det, &(O.l[l].o[j * l]), l, O.NFermion);
          ComplexWaveFunctionAddElementConserveBasis(*Opsi, Det, phase * O.l[l].VRe[j] * psi.ValueRe[ii][ik],
                                                     phase * O.l[l].VRe[j] * psi.ValueIm[ii][ik]);
        }
      }
    }
    free(Det);
  }
  return 0;
}

uint32_t OperatorPsiConserveBasisRR(const OperatorType O, const WaveFunctionType psi, WaveFunctionType * const Opsi)
{
  uint32_t i, ii, ij, ik, j, k, l;
  int    iOMP;
  double phase;
  uint8_t *Det;
  Opsi->N = psi.N;
  Opsi->NBlok = psi.NBlok;
  Opsi->NFermion = psi.NFermion;
  Opsi->NBoson = psi.NBoson;
  Opsi->NWord = psi.NWord;
  Opsi->Next = psi.Next;
  Opsi->Determinant = psi.Determinant;
  Opsi->NBitsKey = psi.NBitsKey;
  Opsi->KeyList = psi.KeyList;
  Opsi->QComplex = psi.QComplex;
  Opsi->ValueRe = malloc(Opsi->NBlok * sizeof(double *));
  if (Opsi->ValueRe == NULL)
  {                             /* Allocation Failed */
    printf("alloc failed OperatorPsiConserveBasis 01\n");
    fflush(stdout);
    Opsi->N = 0;
    Opsi->NBlok = 0;
    return 1;
  }
  for (i = 0; i < Opsi->NBlok; i++)
  {
    Opsi->ValueRe[i] = calloc(BlokSize, sizeof(double));
    if (Opsi->ValueRe[i] == NULL)
    {                           /* Allocation Failed */
      printf("alloc failed OperatorPsiConserveBasis 02\n");
      fflush(stdout);
      for (j = 0; j < i; j++)
      {
        free(Opsi->ValueRe[j]);
      }
      Opsi->N = 0;
      Opsi->NBlok = 0;
      return 1;
    }
  }
#pragma omp parallel private(Det,i,j,k,l,ii,ij,ik,phase,iOMP)
  {
    Det = malloc(Opsi->NWord * sizeof(uint8_t));
    if (Det == NULL)
      AllocFailed();
#pragma omp for
    for (iOMP = 1; iOMP <= (int)psi.N; iOMP++)
    {
      i = (uint32_t) iOMP;
      ii = i >> NBitsBlok;
      ik = i & BlokMask;
      ij = psi.NWord * ik;
      for (l = 0; l <= O.MaxLength; l++)
      {                         /* act with a lader operator of length l on Determinant */
        for (j = 0; j < O.l[l].N; j++)
        {                       /* take opperator number k of length l and act on Determinant */
          for (k = 0; k < Opsi->NWord; k++)
            Det[k] = psi.Determinant[ii][ij + k];
          phase = OperatorLader(Det, &(O.l[l].o[j * l]), l, O.NFermion);
          RealWaveFunctionAddElementConserveBasis(Opsi, Det, phase * O.l[l].VRe[j] * psi.ValueRe[ii][ik]);
        }
      }
    }
    free(Det);
  }
  return 0;
}

uint32_t OperatorPsiConserveBasisMC(const OperatorType O, const WaveFunctionType psi, WaveFunctionType * const Opsi)
{
  uint32_t i, ii, ij, ik, j, k, l;
  int    iOMP;
  double phase;
  uint8_t *Det;
  Opsi->N = psi.N;
  Opsi->NBlok = psi.NBlok;
  Opsi->NFermion = psi.NFermion;
  Opsi->NBoson = psi.NBoson;
  Opsi->NWord = psi.NWord;
  Opsi->Next = psi.Next;
  Opsi->Determinant = psi.Determinant;
  Opsi->NBitsKey = psi.NBitsKey;
  Opsi->KeyList = psi.KeyList;
  Opsi->QComplex = psi.QComplex;
  Opsi->ValueRe = malloc(Opsi->NBlok * sizeof(double *));
  if (Opsi->ValueRe == NULL)
  {                             /* Allocation Failed */
    printf("alloc failed OperatorPsiConserveBasis 01Re\n");
    fflush(stdout);
    Opsi->N = 0;
    Opsi->NBlok = 0;
    return 1;
  }
  Opsi->ValueIm = malloc(Opsi->NBlok * sizeof(double *));
  if (Opsi->ValueIm == NULL)
  {                             /* Allocation Failed */
    printf("alloc failed OperatorPsiConserveBasis 01Im\n");
    fflush(stdout);
    Opsi->N = 0;
    Opsi->NBlok = 0;
    free(Opsi->ValueRe);
    return 1;
  }
  for (i = 0; i < Opsi->NBlok; i++)
  {
    Opsi->ValueRe[i] = calloc(BlokSize, sizeof(double));
    Opsi->ValueIm[i] = calloc(BlokSize, sizeof(double));
    if (Opsi->ValueIm[i] == NULL)
    {                           /* Allocation Failed */
      printf("alloc failed OperatorPsiConserveBasis 02\n");
      fflush(stdout);
      if (Opsi->ValueRe[i] != NULL)
        free(Opsi->ValueRe[i]);
      for (j = 0; j < i; j++)
      {
        free(Opsi->ValueRe[j]);
        free(Opsi->ValueIm[j]);
      }
      Opsi->N = 0;
      Opsi->NBlok = 0;
      free(Opsi->ValueRe);
      free(Opsi->ValueIm);
      return 1;
    }
  }
#pragma omp parallel private(Det,i,j,k,l,ii,ij,ik,phase,iOMP)
  {
    Det = malloc(Opsi->NWord * sizeof(uint8_t));
    if (Det == NULL)
      AllocFailed();
#pragma omp for
    for (iOMP = 1; iOMP <= (int)psi.N; iOMP++)
    {
      i = (uint32_t) iOMP;
      ii = i >> NBitsBlok;
      ik = i & BlokMask;
      ij = psi.NWord * ik;
      if (O.l[0].N != 0)
        for (l = 0; l <= O.MaxLength; l++)
        {                       /* act with a lader operator of length l on Determinant */
          if (O.l[l].QComplex)
          {
            for (j = 0; j < O.l[l].N; j++)
            {                   /* take opperator number k of length l and act on Determinant */
              for (k = 0; k < Opsi->NWord; k++)
                Det[k] = psi.Determinant[ii][ij + k];
              phase = OperatorLader(Det, &(O.l[l].o[j * l]), l, O.NFermion);
              ComplexWaveFunctionAddElementConserveBasis(*Opsi, Det,
                                                         phase * (O.l[l].V[j].Re * psi.ValueRe[ii][ik] -
                                                                  O.l[l].V[j].Im * psi.ValueIm[ii][ik]),
                                                         phase * (O.l[l].V[j].Re * psi.ValueIm[ii][ik] +
                                                                  O.l[l].V[j].Im * psi.ValueRe[ii][ik]));
            }
          }
          else
          {
            for (j = 0; j < O.l[l].N; j++)
            {                   /* take opperator number k of length l and act on Determinant */
              for (k = 0; k < Opsi->NWord; k++)
                Det[k] = psi.Determinant[ii][ij + k];
              phase = OperatorLader(Det, &(O.l[l].o[j * l]), l, O.NFermion);
              ComplexWaveFunctionAddElementConserveBasis(*Opsi, Det, phase * O.l[l].VRe[j] * psi.ValueRe[ii][ik],
                                                         phase * O.l[l].VRe[j] * psi.ValueIm[ii][ik]);
            }
          }
        }
    }
    free(Det);
  }
  return 0;
}

uint32_t OperatorPsiConserveBasisMR(const OperatorType O, const WaveFunctionType psi, WaveFunctionType * const Opsi)
{
  uint32_t i, ii, ij, ik, j, k, l;
  int    iOMP;
  double phase;
  uint8_t *Det;
  Opsi->N = psi.N;
  Opsi->NBlok = psi.NBlok;
  Opsi->NFermion = psi.NFermion;
  Opsi->NBoson = psi.NBoson;
  Opsi->NWord = psi.NWord;
  Opsi->Next = psi.Next;
  Opsi->Determinant = psi.Determinant;
  Opsi->NBitsKey = psi.NBitsKey;
  Opsi->KeyList = psi.KeyList;
  Opsi->QComplex = 1;
  Opsi->ValueRe = malloc(Opsi->NBlok * sizeof(double *));
  if (Opsi->ValueRe == NULL)
  {                             /* Allocation Failed */
    printf("alloc failed OperatorPsiConserveBasis 01Re\n");
    fflush(stdout);
    Opsi->N = 0;
    Opsi->NBlok = 0;
    return 1;
  }
  Opsi->ValueIm = malloc(Opsi->NBlok * sizeof(double *));
  {                             /* Allocation Failed */
    printf("alloc failed OperatorPsiConserveBasis 01Im\n");
    fflush(stdout);
    Opsi->N = 0;
    Opsi->NBlok = 0;
    free(Opsi->ValueRe);
    return 1;
  }
  for (i = 0; i < Opsi->NBlok; i++)
  {
    Opsi->ValueRe[i] = calloc(BlokSize, sizeof(double));
    Opsi->ValueIm[i] = calloc(BlokSize, sizeof(double));
    if (Opsi->ValueIm[i] == NULL)
    {                           /* Allocation Failed */
      printf("alloc failed OperatorPsiConserveBasis 02\n");
      fflush(stdout);
      if (Opsi->ValueRe[i] != NULL)
        free(Opsi->ValueRe[i]);
      for (j = 0; j < i; j++)
      {
        free(Opsi->ValueRe[j]);
        free(Opsi->ValueIm[j]);
      }
      Opsi->N = 0;
      Opsi->NBlok = 0;
      free(Opsi->ValueRe);
      free(Opsi->ValueIm);
      return 1;
    }
  }
#pragma omp parallel private(Det,i,j,k,l,ii,ij,ik,phase,iOMP)
  {
    Det = malloc(Opsi->NWord * sizeof(uint8_t));
    if (Det == NULL)
      AllocFailed();
#pragma omp for
    for (iOMP = 1; iOMP <= (int)psi.N; iOMP++)
    {
      i = (uint32_t) iOMP;
      ii = i >> NBitsBlok;
      ik = i & BlokMask;
      ij = psi.NWord * ik;
      if (O.l[0].N != 0)
        for (l = 0; l <= O.MaxLength; l++)
        {                       /* act with a lader operator of length l on Determinant */
          if (O.l[l].QComplex)
          {
            for (j = 0; j < O.l[l].N; j++)
            {                   /* take opperator number k of length l and act on Determinant */
              for (k = 0; k < Opsi->NWord; k++)
                Det[k] = psi.Determinant[ii][ij + k];
              phase = OperatorLader(Det, &(O.l[l].o[j * l]), l, O.NFermion);
              ComplexWaveFunctionAddElementConserveBasis(*Opsi, Det, phase * (O.l[l].V[j].Re * psi.ValueRe[ii][ik]),
                                                         phase * (O.l[l].V[j].Im * psi.ValueRe[ii][ik]));
            }
          }
          else
          {
            for (j = 0; j < O.l[l].N; j++)
            {                   /* take opperator number k of length l and act on Determinant */
              for (k = 0; k < Opsi->NWord; k++)
                Det[k] = psi.Determinant[ii][ij + k];
              phase = OperatorLader(Det, &(O.l[l].o[j * l]), l, O.NFermion);
              ComplexWaveFunctionAddElementConserveBasis(*Opsi, Det, phase * O.l[l].VRe[j] * psi.ValueRe[ii][ik], 0.0);
            }
          }
        }
    }
    free(Det);
  }
  return 0;
}

uint32_t PsiAddOperatorPsiConserveBasisPreAllocatedCC(const OperatorType *O, const WaveFunctionType *psi,
                                                      WaveFunctionType * const Opsi)
{
  return 0;
}

uint32_t PsiAddOperatorPsiConserveBasisPreAllocatedCR(const OperatorType *O, const WaveFunctionType *psi,
                                                      WaveFunctionType * const Opsi)
{
  return 0;
}

uint32_t PsiAddOperatorPsiConserveBasisPreAllocatedRC(const OperatorType *O, const WaveFunctionType *psi,
                                                      WaveFunctionType * const Opsi)
{
  return 0;
}

uint32_t PsiAddOperatorPsiConserveBasisPreAllocatedRR(const OperatorType *O, const WaveFunctionType *psi,
                                                      WaveFunctionType * const Opsi)
{
  return 0;
}

uint32_t PsiAddOperatorPsiConserveBasisPreAllocatedRRC(const OperatorType *O, const WaveFunctionType *psi,
                                                       WaveFunctionType * const Opsi)
{
  return 0;
}

uint32_t PsiAddOperatorPsiConserveBasisPreAllocatedMC(const OperatorType *O, const WaveFunctionType *psi,
                                                      WaveFunctionType * const Opsi)
{
  return 0;
}

uint32_t PsiAddOperatorPsiConserveBasisPreAllocatedMR(const OperatorType *O, const WaveFunctionType *psi,
                                                      WaveFunctionType * const Opsi)
{
  return 0;
}

uint32_t OperatorPsiConserveBasisPreAllocatedCC(const OperatorType *O, const WaveFunctionType *psi,
                                                WaveFunctionType * const Opsi)
{
  uint32_t i, j, k, l, ii, ij, ik;
  double phase;
  uint8_t *Det;
  int    iOMP;
  if (!(Opsi->QComplex))
    if (MakeComplexWaveFunctionFromReal(Opsi, *Opsi))
    {
      printf
        ("OperatorPsiConserveBasisPreAllocatedCC called with allocation of Opsi being real.\nAllocation faild when tried to expand real to complex\n");
      return 1;
    }
#pragma omp parallel for private(i,j)
  for (iOMP = 0; iOMP < Opsi->NBlok; iOMP++)
  {
    i = (uint32_t) iOMP;
    for (j = 0; j < BlokSize; j++)
    {
      Opsi->ValueRe[i][j] = 0;
      Opsi->ValueIm[i][j] = 0;
    }
  }
#pragma omp parallel private(Det,i,j,k,l,ii,ij,ik,phase,iOMP)
  {
    Det = malloc(Opsi->NWord * sizeof(uint8_t));
    if (Det == NULL)
      AllocFailed();
#pragma omp for
    for (iOMP = 1; iOMP <= (int)psi->N; iOMP++)
    {
      i = (uint32_t) iOMP;
      ii = i >> NBitsBlok;
      ik = i & BlokMask;
      ij = psi->NWord * ik;
      for (l = 0; l <= O->MaxLength; l++)
      {                         /* act with a lader operator of length l on Determinant */
        for (j = 0; j < O->l[l].N; j++)
        {                       /* take opperator number k of length l and act on Determinant */
          for (k = 0; k < Opsi->NWord; k++)
            Det[k] = psi->Determinant[ii][ij + k];
          phase = OperatorLader(Det, &(O->l[l].o[j * l]), l, O->NFermion);
          ComplexWaveFunctionAddElementConserveBasis(*Opsi, Det,
                                                     phase * (O->l[l].V[j].Re * psi->ValueRe[ii][ik] -
                                                              O->l[l].V[j].Im * psi->ValueIm[ii][ik]),
                                                     phase * (O->l[l].V[j].Re * psi->ValueIm[ii][ik] +
                                                              O->l[l].V[j].Im * psi->ValueRe[ii][ik]));
        }
      }
    }
    free(Det);
  }
  return 0;
}

uint32_t OperatorPsiConserveBasisPreAllocatedCR(const OperatorType *O, const WaveFunctionType *psi,
                                                WaveFunctionType * const Opsi)
{
  uint32_t i, j, k, l, ii, ij, ik;
  double phase;
  uint8_t *Det;
  int    iOMP;
  if (!(Opsi->QComplex))
    if (MakeComplexWaveFunctionFromReal(Opsi, *Opsi))
    {
      printf
        ("OperatorPsiConserveBasisPreAllocatedCR called with allocation of Opsi being real.\nAllocation faild when tried to expand real to complex\n");
      return 1;
    }
#pragma omp parallel for private(i,j)
  for (iOMP = 0; iOMP < Opsi->NBlok; iOMP++)
  {
    i = (uint32_t) iOMP;
    for (j = 0; j < BlokSize; j++)
    {
      Opsi->ValueRe[i][j] = 0;
      Opsi->ValueIm[i][j] = 0;
    }
  }
#pragma omp parallel private(Det,i,j,k,l,ii,ij,ik,phase,iOMP)
  {
    Det = malloc(Opsi->NWord * sizeof(uint8_t));
    if (Det == NULL)
      AllocFailed();
#pragma omp for
    for (iOMP = 1; iOMP <= (int)psi->N; iOMP++)
    {
      i = (uint32_t) iOMP;
      ii = i >> NBitsBlok;
      ik = i & BlokMask;
      ij = psi->NWord * ik;
      for (l = 0; l <= O->MaxLength; l++)
      {                         /* act with a lader operator of length l on Determinant */
        for (j = 0; j < O->l[l].N; j++)
        {                       /* take opperator number k of length l and act on Determinant */
          for (k = 0; k < Opsi->NWord; k++)
            Det[k] = psi->Determinant[ii][ij + k];
          phase = OperatorLader(Det, &(O->l[l].o[j * l]), l, O->NFermion);
          ComplexWaveFunctionAddElementConserveBasis(*Opsi, Det, phase * (O->l[l].V[j].Re * psi->ValueRe[ii][ik]),
                                                     phase * (O->l[l].V[j].Im * psi->ValueRe[ii][ik]));
        }
      }
    }
    free(Det);
  }
  return 0;
}

uint32_t OperatorPsiConserveBasisPreAllocatedRC(const OperatorType *O, const WaveFunctionType *psi,
                                                WaveFunctionType * const Opsi)
{
  uint32_t i, j, k, l, ii, ij, ik;
  double phase;
  uint8_t *Det;
  int    iOMP;
  if (!(Opsi->QComplex))
    if (MakeComplexWaveFunctionFromReal(Opsi, *Opsi))
    {
      printf
        ("OperatorPsiConserveBasisPreAllocatedCC called with allocation of Opsi being real.\nAllocation faild when tried to expand real to complex\n");
      return 1;
    }
#pragma omp parallel for private(i,j)
  for (iOMP = 0; iOMP < Opsi->NBlok; iOMP++)
  {
    i = (uint32_t) iOMP;
    for (j = 0; j < BlokSize; j++)
    {
      Opsi->ValueRe[i][j] = 0;
      Opsi->ValueIm[i][j] = 0;
    }
  }
#pragma omp parallel private(Det,i,j,k,l,ii,ij,ik,phase,iOMP)
  {
    Det = malloc(Opsi->NWord * sizeof(uint8_t));
    if (Det == NULL)
      AllocFailed();
#pragma omp for
    for (iOMP = 1; iOMP <= (int)psi->N; iOMP++)
    {
      i = (uint32_t) iOMP;
      ii = i >> NBitsBlok;
      ik = i & BlokMask;
      ij = psi->NWord * ik;
      for (l = 0; l <= O->MaxLength; l++)
      {                         /* act with a lader operator of length l on Determinant */
        for (j = 0; j < O->l[l].N; j++)
        {                       /* take opperator number k of length l and act on Determinant */
          for (k = 0; k < Opsi->NWord; k++)
            Det[k] = psi->Determinant[ii][ij + k];
          phase = OperatorLader(Det, &(O->l[l].o[j * l]), l, O->NFermion);
          ComplexWaveFunctionAddElementConserveBasis(*Opsi, Det, phase * O->l[l].VRe[j] * psi->ValueRe[ii][ik],
                                                     phase * O->l[l].VRe[j] * psi->ValueIm[ii][ik]);
        }
      }
    }
    free(Det);
  }
  return 0;
}

uint32_t OperatorPsiConserveBasisPreAllocatedRR(const OperatorType *O, const WaveFunctionType *psi,
                                                WaveFunctionType * const Opsi)
{
  uint32_t i, j, k, l, ii, ij, ik;
  double phase;
  uint8_t *Det;
  int    iOMP;
  if (Opsi->QComplex)
    MakeRealWaveFunctionFromComplex(*Opsi, Opsi);
#pragma omp parallel for private(i,j)
  for (iOMP = 0; iOMP < Opsi->NBlok; iOMP++)
  {
    i = (uint32_t) iOMP;
    for (j = 0; j < BlokSize; j++)
    {
      Opsi->ValueRe[i][j] = 0;
    }
  }
#pragma omp parallel private(Det,i,j,k,l,ii,ij,ik,phase,iOMP)
  {
    Det = malloc(Opsi->NWord * sizeof(uint8_t));
    if (Det == NULL)
      AllocFailed();
#pragma omp for
    for (iOMP = 1; iOMP <= (int)psi->N; iOMP++)
    {
      i = (uint32_t) iOMP;
      ii = i >> NBitsBlok;
      ik = i & BlokMask;
      ij = psi->NWord * ik;
      for (l = 0; l <= O->MaxLength; l++)
      {                         /* act with a lader operator of length l on Determinant */
        for (j = 0; j < O->l[l].N; j++)
        {                       /* take opperator number k of length l and act on Determinant */
          for (k = 0; k < Opsi->NWord; k++)
            Det[k] = psi->Determinant[ii][ij + k];
          phase = OperatorLader(Det, &(O->l[l].o[j * l]), l, O->NFermion);
          RealWaveFunctionAddElementConserveBasis(Opsi, Det, phase * O->l[l].VRe[j] * psi->ValueRe[ii][ik]);
        }
      }
    }
    free(Det);
  }
  return 0;
}

uint32_t OperatorPsiConserveBasisPreAllocatedMC(const OperatorType *O, const WaveFunctionType *psi,
                                                WaveFunctionType * const Opsi)
{
  uint32_t i, j, k, l, ii, ij, ik;
  double phase;
  uint8_t *Det;
  int    iOMP;
  if (!(Opsi->QComplex))
    if (MakeComplexWaveFunctionFromReal(Opsi, *Opsi))
    {
      printf("OperatorPsiConserveBasisPreAllocatedCC called with allocation of Opsi being real.\nAllocation faild when tried to expand real to complex\n");
      return 1;
    }
#pragma omp parallel for private(i,j)
  for (iOMP = 0; iOMP < Opsi->NBlok; iOMP++)
  {
    i = (uint32_t) iOMP;
    for (j = 0; j < BlokSize; j++)
    {
      Opsi->ValueRe[i][j] = 0;
      Opsi->ValueIm[i][j] = 0;
    }
  }
#pragma omp parallel private(Det,i,j,k,l,ii,ij,ik,phase,iOMP)
  {
    Det = malloc(Opsi->NWord * sizeof(uint8_t));
    if (Det == NULL)
      AllocFailed();
#pragma omp for
    for (iOMP = 1; iOMP <= (int)psi->N; iOMP++)
    {
      i = (uint32_t) iOMP;
      ii = i >> NBitsBlok;
      ik = i & BlokMask;
      ij = psi->NWord * ik;
      for (l = 0; l <= O->MaxLength; l++)
      {                         /* act with a lader operator of length l on Determinant */
        if (O->l[l].QComplex)
        {
          for (j = 0; j < O->l[l].N; j++)
          {                     /* take opperator number k of length l and act on Determinant */
            for (k = 0; k < Opsi->NWord; k++)
              Det[k] = psi->Determinant[ii][ij + k];
            phase = OperatorLader(Det, &(O->l[l].o[j * l]), l, O->NFermion);
            ComplexWaveFunctionAddElementConserveBasis(*Opsi, Det,
                                                       phase * (O->l[l].V[j].Re * psi->ValueRe[ii][ik] -
                                                                O->l[l].V[j].Im * psi->ValueIm[ii][ik]),
                                                       phase * (O->l[l].V[j].Re * psi->ValueIm[ii][ik] +
                                                                O->l[l].V[j].Im * psi->ValueRe[ii][ik]));
          }
        }
        else
        {
          for (j = 0; j < O->l[l].N; j++)
          {                     /* take opperator number k of length l and act on Determinant */
            for (k = 0; k < Opsi->NWord; k++)
              Det[k] = psi->Determinant[ii][ij + k];
            phase = OperatorLader(Det, &(O->l[l].o[j * l]), l, O->NFermion);
            ComplexWaveFunctionAddElementConserveBasis(*Opsi, Det, phase * O->l[l].VRe[j] * psi->ValueRe[ii][ik],
                                                       phase * O->l[l].VRe[j] * psi->ValueIm[ii][ik]);
          }
        }
      }
    }
    free(Det);
  }
  return 0;
}

uint32_t OperatorPsiConserveBasisPreAllocatedMR(const OperatorType *O, const WaveFunctionType *psi,
                                                WaveFunctionType * const Opsi)
{
  uint32_t i, j, k, l, ii, ij, ik;
  double phase;
  uint8_t *Det;
  int    iOMP;
  if (!(Opsi->QComplex))
    if (MakeComplexWaveFunctionFromReal(Opsi, *Opsi))
    {
      printf("OperatorPsiConserveBasisPreAllocatedCR called with allocation of Opsi being real.\nAllocation faild when tried to expand real to complex\n");
      return 1;
    }
#pragma omp parallel for private(i,j)
  for (iOMP = 0; iOMP < Opsi->NBlok; iOMP++)
  {
    i = (uint32_t) iOMP;
    for (j = 0; j < BlokSize; j++)
    {
      Opsi->ValueRe[i][j] = 0;
      Opsi->ValueIm[i][j] = 0;
    }
  }
#pragma omp parallel private(Det,i,j,k,l,ii,ij,ik,phase,iOMP)
  {
    Det = malloc(Opsi->NWord * sizeof(uint8_t));
    if (Det == NULL)
      AllocFailed();
#pragma omp for
    for (iOMP = 1; iOMP <= (int)psi->N; iOMP++)
    {
      i = (uint32_t) iOMP;
      ii = i >> NBitsBlok;
      ik = i & BlokMask;
      ij = psi->NWord * ik;
      for (l = 0; l <= O->MaxLength; l++)
      {                         /* act with a lader operator of length l on Determinant */
        if (O->l[l].QComplex)
        {
          for (j = 0; j < O->l[l].N; j++)
          {                     /* take opperator number k of length l and act on Determinant */
            for (k = 0; k < Opsi->NWord; k++)
              Det[k] = psi->Determinant[ii][ij + k];
            phase = OperatorLader(Det, &(O->l[l].o[j * l]), l, O->NFermion);
            ComplexWaveFunctionAddElementConserveBasis(*Opsi, Det, phase * (O->l[l].V[j].Re * psi->ValueRe[ii][ik]),
                                                       phase * (O->l[l].V[j].Im * psi->ValueRe[ii][ik]));
          }
        }
        else
        {
          for (j = 0; j < O->l[l].N; j++)
          {                     /* take opperator number k of length l and act on Determinant */
            for (k = 0; k < Opsi->NWord; k++)
              Det[k] = psi->Determinant[ii][ij + k];
            phase = OperatorLader(Det, &(O->l[l].o[j * l]), l, O->NFermion);
            ComplexWaveFunctionAddElementConserveBasis(*Opsi, Det, phase * O->l[l].VRe[j] * psi->ValueRe[ii][ik], 0.0);
          }
        }
      }
    }
    free(Det);
  }
  return 0;
}

uint32_t PsiAddOperatorPsiCC(const OperatorType O, const WaveFunctionType psi, WaveFunctionType * const Opsi)
{
  return 0;
}

uint32_t PsiAddOperatorPsiCR(const OperatorType O, const WaveFunctionType psi, WaveFunctionType * const Opsi)
{
  return 0;
}

uint32_t PsiAddOperatorPsiRC(const OperatorType O, const WaveFunctionType psi, WaveFunctionType * const Opsi)
{
  return 0;
}

uint32_t PsiAddOperatorPsiRR(const OperatorType O, const WaveFunctionType psi, WaveFunctionType * const Opsi)
{
  return 0;
}

uint32_t PsiAddOperatorPsiRRC(const OperatorType O, const WaveFunctionType psi, WaveFunctionType * const Opsi)
{
  return 0;
}

uint32_t PsiAddOperatorPsiMC(const OperatorType O, const WaveFunctionType psi, WaveFunctionType * const Opsi)
{
  return 0;
}

uint32_t PsiAddOperatorPsiMR(const OperatorType O, const WaveFunctionType psi, WaveFunctionType * const Opsi)
{
  return 0;
}

uint32_t OperatorPsiCC(const OperatorType O, const WaveFunctionType psi, WaveFunctionType * const Opsi)
{
  uint32_t i, j, k, l, ii, ij, ik, OutOffMem;
  int    iOMP;
  double phase;
  Opsi->QComplex = 1;
  Opsi->NFermion = psi.NFermion;
  Opsi->NBoson = psi.NBoson;
  Opsi->NBitsKey = psi.NBitsKey;
  if (WaveFunctionInit(Opsi))
  {
    printf("Could not initialize Opsi in OperatorPsiCC\n");
    return 1;
  };
  OutOffMem = 0;
#pragma omp parallel private(i,j,k,l,ii,ij,ik,phase)
  {
    uint8_t Det[Opsi->NWord];
    uint8_t MiniDet[MiniBlokSize * Opsi->NWord];
    double MiniValueRe[MiniBlokSize];
    double MiniValueIm[MiniBlokSize];
    uint32_t NMini;
    NMini = 0;
#pragma omp for
    for (iOMP = 1; iOMP <= psi.N; iOMP++)
    {
      if (!OutOffMem)
      {
        i = (uint32_t) iOMP;
        ii = i >> NBitsBlok;
        ik = i & BlokMask;
        ij = psi.NWord * ik;
        for (l = 0; l <= O.MaxLength; l++)
        {
          for (j = 0; j < O.l[l].N; j++)
          {
            for (k = 0; k < Opsi->NWord; k++)
              Det[k] = psi.Determinant[ii][ij + k];
            phase = OperatorLader(Det, &(O.l[l].o[j * l]), l, O.NFermion);
            if (fabs(phase) > (10 * DBL_MIN))
            {
              MiniValueRe[NMini] =
                phase * (O.l[l].V[j].Re * psi.ValueRe[ii][ik] - O.l[l].V[j].Im * psi.ValueIm[ii][ik]);
              MiniValueIm[NMini] =
                phase * (O.l[l].V[j].Re * psi.ValueIm[ii][ik] + O.l[l].V[j].Im * psi.ValueRe[ii][ik]);
              for (k = 0; k < Opsi->NWord; k++)
                MiniDet[NMini * Opsi->NWord + k] = Det[k];
              NMini++;
              if (NMini == MiniBlokSize)
              {                 /* flush Mini buffer */
                if (ComplexWaveFunctionAddElementOMPMiniFlush(Opsi, &NMini, MiniDet, MiniValueRe, MiniValueIm))
                {               /* Allocation failed */
                  printf("ComplexWaveFunctionAddElementOMPMiniFlush failed in OperatorPsiCC\n 01");
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
    if (ComplexWaveFunctionAddElementOMPMiniFlush(Opsi, &NMini, MiniDet, MiniValueRe, MiniValueIm))
    {                           /* Allocation failed */
      printf("ComplexWaveFunctionAddElementOMPMiniFlush failed in OperatorPsiCC 02\n");
      fflush(stdout);
      OutOffMem = 1;
#pragma omp flush(OutOffMem)
    }
  }                             /* end of parallel region */
  if (OutOffMem)
  {
    RealWaveFunctionFree(Opsi);
    return 1;
  }
  return 0;
}

uint32_t OperatorPsiCR(const OperatorType O, const WaveFunctionType psi, WaveFunctionType * const Opsi)
{
  uint32_t i, j, k, l, ii, ij, ik, OutOffMem;
  int    iOMP;
  double phase;
  Opsi->QComplex = 1;
  Opsi->NFermion = psi.NFermion;
  Opsi->NBoson = psi.NBoson;
  Opsi->NBitsKey = psi.NBitsKey;
  if (WaveFunctionInit(Opsi))
  {
    printf("Could not initialize Opsi in OperatorPsiCR\n");
    return 1;
  };
  OutOffMem = 0;
#pragma omp parallel private(i,j,k,l,ii,ij,ik,phase)
  {
    uint8_t Det[Opsi->NWord];
    uint8_t MiniDet[MiniBlokSize * Opsi->NWord];
    double MiniValueRe[MiniBlokSize];
    double MiniValueIm[MiniBlokSize];
    uint32_t NMini;
    NMini = 0;
#pragma omp for
    for (iOMP = 1; iOMP <= psi.N; iOMP++)
    {
      if (!OutOffMem)
      {
        i = (uint32_t) iOMP;
        ii = i >> NBitsBlok;
        ik = i & BlokMask;
        ij = psi.NWord * ik;
        for (l = 0; l <= O.MaxLength; l++)
        {
          for (j = 0; j < O.l[l].N; j++)
          {
            for (k = 0; k < Opsi->NWord; k++)
              Det[k] = psi.Determinant[ii][ij + k];
            phase = OperatorLader(Det, &(O.l[l].o[j * l]), l, O.NFermion);
            if (fabs(phase) > (10 * DBL_MIN))
            {
              MiniValueRe[NMini] = phase * O.l[l].V[j].Re * psi.ValueRe[ii][ik];
              MiniValueIm[NMini] = phase * O.l[l].V[j].Im * psi.ValueRe[ii][ik];
              for (k = 0; k < Opsi->NWord; k++)
                MiniDet[NMini * Opsi->NWord + k] = Det[k];
              NMini++;
              if (NMini == MiniBlokSize)
              {                 /* flush Mini buffer */
                if (ComplexWaveFunctionAddElementOMPMiniFlush(Opsi, &NMini, MiniDet, MiniValueRe, MiniValueIm))
                {               /* Allocation failed */
                  printf("ComplexWaveFunctionAddElementOMPMiniFlush failed in OperatorPsiCR\n 01");
                  fflush(stdout);
                  OutOffMem = 1;
#pragma omp flush(OutOffMem)
                  j = O.l[l].N;
                  l = O.MaxLength;
                }               /* end of allocation failed */
              }                 /* end of flush Mini buffer */
            }                   /* end of adding a new element */
          }                     /* end of loop over different operators of lenght l */
        }                       /* end of loop over operators of different length */
      }                         /* Out Off Mem test */
    }                           /* end of loop over elements of psi this loop runs parallel */
    if (ComplexWaveFunctionAddElementOMPMiniFlush(Opsi, &NMini, MiniDet, MiniValueRe, MiniValueIm))
    {                           /* Allocation failed */
      printf("ComplexWaveFunctionAddElementOMPMiniFlush failed in OperatorPsiCR 02\n");
      fflush(stdout);
      OutOffMem = 1;
#pragma omp flush(OutOffMem)
    }
  }                             /* end of parallel region */
  if (OutOffMem)
  {
    RealWaveFunctionFree(Opsi);
    return 1;
  }
  return 0;
}

uint32_t OperatorPsiRC(const OperatorType O, const WaveFunctionType psi, WaveFunctionType * const Opsi)
{
  uint32_t i, j, k, l, ii, ij, ik, OutOffMem;
  int    iOMP;
  double phase;
  Opsi->QComplex = 1;
  Opsi->NFermion = psi.NFermion;
  Opsi->NBoson = psi.NBoson;
  Opsi->NBitsKey = psi.NBitsKey;
  if (WaveFunctionInit(Opsi))
  {
    printf("Could not initialize Opsi in OperatorPsiRC\n");
    return 1;
  };
  OutOffMem = 0;
#pragma omp parallel private(i,j,k,l,ii,ij,ik,phase)
  {
    uint8_t Det[Opsi->NWord];
    uint8_t MiniDet[MiniBlokSize * Opsi->NWord];
    double MiniValueRe[MiniBlokSize];
    double MiniValueIm[MiniBlokSize];
    uint32_t NMini;
    NMini = 0;
#pragma omp for
    for (iOMP = 1; iOMP <= psi.N; iOMP++)
    {
      if (!OutOffMem)
      {
        i = (uint32_t) iOMP;
        ii = i >> NBitsBlok;
        ik = i & BlokMask;
        ij = psi.NWord * ik;
        for (l = 0; l <= O.MaxLength; l++)
        {
          for (j = 0; j < O.l[l].N; j++)
          {
            for (k = 0; k < Opsi->NWord; k++)
              Det[k] = psi.Determinant[ii][ij + k];
            phase = OperatorLader(Det, &(O.l[l].o[j * l]), l, O.NFermion);
            if (fabs(phase) > (10 * DBL_MIN))
            {
              MiniValueRe[NMini] = phase * O.l[l].VRe[j] * psi.ValueRe[ii][ik];
              MiniValueIm[NMini] = phase * O.l[l].VRe[j] * psi.ValueIm[ii][ik];
              for (k = 0; k < Opsi->NWord; k++)
                MiniDet[NMini * Opsi->NWord + k] = Det[k];
              NMini++;
              if (NMini == MiniBlokSize)
              {                 /* flush Mini buffer */
                if (ComplexWaveFunctionAddElementOMPMiniFlush(Opsi, &NMini, MiniDet, MiniValueRe, MiniValueIm))
                {               /* Allocation failed */
                  printf("ComplexWaveFunctionAddElementOMPMiniFlush failed in OperatorPsiRC\n 01");
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
    if (ComplexWaveFunctionAddElementOMPMiniFlush(Opsi, &NMini, MiniDet, MiniValueRe, MiniValueIm))
    {                           /* Allocation failed */
      printf("ComplexWaveFunctionAddElementOMPMiniFlush failed in OperatorPsiRC 02\n");
      fflush(stdout);
      OutOffMem = 1;
#pragma omp flush(OutOffMem)
    }
  }                             /* end of parallel region */
  if (OutOffMem)
  {
    RealWaveFunctionFree(Opsi);
    return 1;
  }
  return 0;
}

uint32_t OperatorPsiRR(const OperatorType O, const WaveFunctionType psi, WaveFunctionType * const Opsi)
{
  uint32_t i, j, k, l, ii, ij, ik, OutOffMem;
  int    iOMP;
  double phase;
  Opsi->QComplex = 0;
  Opsi->NFermion = psi.NFermion;
  Opsi->NBoson = psi.NBoson;
  Opsi->NBitsKey = psi.NBitsKey;
  if (WaveFunctionInit(Opsi))
  {
    printf("Could not initialize Opsi in OperatorPsiRR\n");
    return 1;
  };
  OutOffMem = 0;
#pragma omp parallel private(i,j,k,l,ii,ij,ik,phase)
  {
    uint8_t Det[Opsi->NWord];
    uint8_t MiniDet[MiniBlokSize * Opsi->NWord];
    double MiniValue[MiniBlokSize];
    uint32_t NMini;
    NMini = 0;
#pragma omp for
    for (iOMP = 1; iOMP <= psi.N; iOMP++)
    {
      if (!OutOffMem)
      {
        i = (uint32_t) iOMP;
        ii = i >> NBitsBlok;
        ik = i & BlokMask;
        ij = psi.NWord * ik;
        for (l = 0; l <= O.MaxLength; l++)
        {
          for (j = 0; j < O.l[l].N; j++)
          {
            for (k = 0; k < Opsi->NWord; k++)
              Det[k] = psi.Determinant[ii][ij + k];
            phase = OperatorLader(Det, &(O.l[l].o[j * l]), l, O.NFermion);
            if (fabs(phase) > (10 * DBL_MIN))
            {
              MiniValue[NMini] = phase * O.l[l].VRe[j] * psi.ValueRe[ii][ik];
              for (k = 0; k < Opsi->NWord; k++)
                MiniDet[NMini * Opsi->NWord + k] = Det[k];
              NMini++;
              if (NMini == MiniBlokSize)
              {                 /* flush Mini buffer */
                if (RealWaveFunctionAddElementOMPMiniFlush(Opsi, &NMini, MiniDet, MiniValue))
                {               /* Allocation failed */
                  printf("RealWaveFunctionAddElementOMPMiniFlush failed in OperatorPsiRR\n 01");
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
    if (RealWaveFunctionAddElementOMPMiniFlush(Opsi, &NMini, MiniDet, MiniValue))
    {                           /* Allocation failed */
      printf("RealWaveFunctionAddElementOMPMiniFlush failed in OperatorPsiRR 02\n");
      fflush(stdout);
      OutOffMem = 1;
#pragma omp flush(OutOffMem)
    }
  }                             /* end of parallel region */
  if (OutOffMem)
  {
    RealWaveFunctionFree(Opsi);
    return 1;
  }
  return 0;
}

uint32_t OperatorPsiMC(const OperatorType O, const WaveFunctionType psi, WaveFunctionType * const Opsi)
{
  uint32_t i, j, k, l, ii, ij, ik, OutOffMem;
  int    iOMP;
  double phase;
  Opsi->QComplex = 1;
  Opsi->NFermion = psi.NFermion;
  Opsi->NBoson = psi.NBoson;
  Opsi->NBitsKey = psi.NBitsKey;
  if (WaveFunctionInit(Opsi))
  {
    printf("Could not initialize Opsi in OperatorPsiCC\n");
    return 1;
  };
  OutOffMem = 0;
#pragma omp parallel private(i,j,k,l,ii,ij,ik,phase)
  {
    uint8_t Det[Opsi->NWord];
    uint8_t MiniDet[MiniBlokSize * Opsi->NWord];
    double MiniValueRe[MiniBlokSize];
    double MiniValueIm[MiniBlokSize];
    uint32_t NMini;
    NMini = 0;
#pragma omp for
    for (iOMP = 1; iOMP <= psi.N; iOMP++)
    {
      if (!OutOffMem)
      {
        i = (uint32_t) iOMP;
        ii = i >> NBitsBlok;
        ik = i & BlokMask;
        ij = psi.NWord * ik;
        for (l = 0; l <= O.MaxLength; l++)
        {
          if (O.l[l].QComplex == 0)
          {                     /* Operator is real at this length */
            for (j = 0; j < O.l[l].N; j++)
            {
              for (k = 0; k < Opsi->NWord; k++)
                Det[k] = psi.Determinant[ii][ij + k];
              phase = OperatorLader(Det, &(O.l[l].o[j * l]), l, O.NFermion);
              if (fabs(phase) > (10 * DBL_MIN))
              {
                MiniValueRe[NMini] = phase * O.l[l].VRe[j] * psi.ValueRe[ii][ik];
                MiniValueIm[NMini] = phase * O.l[l].VRe[j] * psi.ValueIm[ii][ik];
                for (k = 0; k < Opsi->NWord; k++)
                  MiniDet[NMini * Opsi->NWord + k] = Det[k];
                NMini++;
                if (NMini == MiniBlokSize)
                {               /* flush Mini buffer */
                  if (ComplexWaveFunctionAddElementOMPMiniFlush(Opsi, &NMini, MiniDet, MiniValueRe, MiniValueIm))
                  {             /* Allocation failed */
                    printf("ComplexWaveFunctionAddElementOMPMiniFlush failed in OperatorPsiRC\n 01");
                    fflush(stdout);
                    OutOffMem = 1;
#pragma omp flush(OutOffMem)
                    j = O.l[l].N;
                    l = O.MaxLength;
                  }
                }               /* end of condition that Mini buffer is full and needs to be flushed */
              }                 /* end of coditional that phase is above numerical zero */
            }                   /* end of loop over operators of length l (counter j) */
          }
          else
          {                     /* Operator is complex at this length */
            for (j = 0; j < O.l[l].N; j++)
            {
              for (k = 0; k < Opsi->NWord; k++)
                Det[k] = psi.Determinant[ii][ij + k];
              phase = OperatorLader(Det, &(O.l[l].o[j * l]), l, O.NFermion);
              if (fabs(phase) > (10 * DBL_MIN))
              {
                MiniValueRe[NMini] =
                  phase * (O.l[l].V[j].Re * psi.ValueRe[ii][ik] - O.l[l].V[j].Im * psi.ValueIm[ii][ik]);
                MiniValueIm[NMini] =
                  phase * (O.l[l].V[j].Re * psi.ValueIm[ii][ik] + O.l[l].V[j].Im * psi.ValueRe[ii][ik]);
                for (k = 0; k < Opsi->NWord; k++)
                  MiniDet[NMini * Opsi->NWord + k] = Det[k];
                NMini++;
                if (NMini == MiniBlokSize)
                {               /* flush Mini buffer */
                  if (ComplexWaveFunctionAddElementOMPMiniFlush(Opsi, &NMini, MiniDet, MiniValueRe, MiniValueIm))
                  {             /* Allocation failed */
                    printf("ComplexWaveFunctionAddElementOMPMiniFlush failed in OperatorPsiCC\n 01");
                    fflush(stdout);
                    OutOffMem = 1;
#pragma omp flush(OutOffMem)
                    j = O.l[l].N;
                    l = O.MaxLength;
                  }
                }
              }
            }
          }                     /* end of loop over different operators of lenght l */
        }                       /* end of loop over operators of different length */
      }                         /* Out Off Mem test */
    }                           /* end of loop over elements of psi this loop runs parallel */
    if (ComplexWaveFunctionAddElementOMPMiniFlush(Opsi, &NMini, MiniDet, MiniValueRe, MiniValueIm))
    {                           /* Allocation failed */
      printf("ComplexWaveFunctionAddElementOMPMiniFlush failed in OperatorPsiCC 02\n");
      fflush(stdout);
      OutOffMem = 1;
#pragma omp flush(OutOffMem)
    }
  }                             /* end of parallel region */
  if (OutOffMem)
  {
    RealWaveFunctionFree(Opsi);
    return 1;
  }
  return 0;
}

uint32_t OperatorPsiMR(const OperatorType O, const WaveFunctionType psi, WaveFunctionType * const Opsi)
{
  uint32_t i, j, k, l, ii, ij, ik, OutOffMem;
  int    iOMP;
  double phase;
  Opsi->QComplex = 1;
  Opsi->NFermion = psi.NFermion;
  Opsi->NBoson = psi.NBoson;
  Opsi->NBitsKey = psi.NBitsKey;
  if (WaveFunctionInit(Opsi))
  {
    printf("Could not initialize Opsi in OperatorPsiCR\n");
    return 1;
  };
  OutOffMem = 0;
#pragma omp parallel private(i,j,k,l,ii,ij,ik,phase)
  {
    uint8_t Det[Opsi->NWord];
    uint8_t MiniDet[MiniBlokSize * Opsi->NWord];
    double MiniValueRe[MiniBlokSize];
    double MiniValueIm[MiniBlokSize];
    uint32_t NMini;
    NMini = 0;
#pragma omp for
    for (iOMP = 1; iOMP <= psi.N; iOMP++)
    {
      if (!OutOffMem)
      {
        i = (uint32_t) iOMP;
        ii = i >> NBitsBlok;
        ik = i & BlokMask;
        ij = psi.NWord * ik;
        for (l = 0; l <= O.MaxLength; l++)
        {
          if (O.l[l].QComplex == 0)
          {                     /* Operator real at this length */
            for (j = 0; j < O.l[l].N; j++)
            {
              for (k = 0; k < Opsi->NWord; k++)
                Det[k] = psi.Determinant[ii][ij + k];
              phase = OperatorLader(Det, &(O.l[l].o[j * l]), l, O.NFermion);
              if (fabs(phase) > (10 * DBL_MIN))
              {
                MiniValueRe[NMini] = phase * O.l[l].VRe[j] * psi.ValueRe[ii][ik];
                for (k = 0; k < Opsi->NWord; k++)
                  MiniDet[NMini * Opsi->NWord + k] = Det[k];
                NMini++;
                if (NMini == MiniBlokSize)
                {               /* flush Mini buffer */
                  if (RealWaveFunctionAddElementOMPMiniFlush(Opsi, &NMini, MiniDet, MiniValueRe))
                  {             /* Allocation failed */
                    printf("RealWaveFunctionAddElementOMPMiniFlush failed in OperatorPsiRR\n 01");
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
          {                     /* Operator complex at this length */
            for (j = 0; j < O.l[l].N; j++)
            {
              for (k = 0; k < Opsi->NWord; k++)
                Det[k] = psi.Determinant[ii][ij + k];
              phase = OperatorLader(Det, &(O.l[l].o[j * l]), l, O.NFermion);
              if (fabs(phase) > (10 * DBL_MIN))
              {
                MiniValueRe[NMini] = phase * O.l[l].V[j].Re * psi.ValueRe[ii][ik];
                MiniValueIm[NMini] = phase * O.l[l].V[j].Im * psi.ValueRe[ii][ik];
                for (k = 0; k < Opsi->NWord; k++)
                  MiniDet[NMini * Opsi->NWord + k] = Det[k];
                NMini++;
                if (NMini == MiniBlokSize)
                {               /* flush Mini buffer */
                  if (ComplexWaveFunctionAddElementOMPMiniFlush(Opsi, &NMini, MiniDet, MiniValueRe, MiniValueIm))
                  {             /* Allocation failed */
                    printf("ComplexWaveFunctionAddElementOMPMiniFlush failed in OperatorPsiCR\n 01");
                    fflush(stdout);
                    OutOffMem = 1;
#pragma omp flush(OutOffMem)
                    j = O.l[l].N;
                    l = O.MaxLength;
                  }             /* end of allocation failed */
                }               /* end of flush Mini buffer */
              }                 /* end of adding a new element */
            }                   /* end of choice between real or complex at this length */
          }                     /* end of loop over different operators of lenght l */
        }                       /* end of loop over operators of different length */
      }                         /* Out Off Mem test */
    }                           /* end of loop over elements of psi this loop runs parallel */
    if (ComplexWaveFunctionAddElementOMPMiniFlush(Opsi, &NMini, MiniDet, MiniValueRe, MiniValueIm))
    {                           /* Allocation failed */
      printf("ComplexWaveFunctionAddElementOMPMiniFlush failed in OperatorPsiCR 02\n");
      fflush(stdout);
      OutOffMem = 1;
#pragma omp flush(OutOffMem)
    }
  }                             /* end of parallel region */
  if (OutOffMem)
  {
    RealWaveFunctionFree(Opsi);
    return 1;
  }
  return 0;
}

uint32_t PsiAddOperatorPsiSerialCC(const OperatorType O, const WaveFunctionType psi, WaveFunctionType * const Opsi)
{
  return 0;
}

uint32_t PsiAddOperatorPsiSerialCR(const OperatorType O, const WaveFunctionType psi, WaveFunctionType * const Opsi)
{
  return 0;
}

uint32_t PsiAddOperatorPsiSerialRC(const OperatorType O, const WaveFunctionType psi, WaveFunctionType * const Opsi)
{
  return 0;
}

uint32_t PsiAddOperatorPsiSerialRR(const OperatorType O, const WaveFunctionType psi, WaveFunctionType * const Opsi)
{
  return 0;
}

uint32_t PsiAddOperatorPsiSerialRRC(const OperatorType O, const WaveFunctionType psi, WaveFunctionType * const Opsi)
{
  return 0;
}

uint32_t PsiAddOperatorPsiSerialMC(const OperatorType O, const WaveFunctionType psi, WaveFunctionType * const Opsi)
{
  return 0;
}

uint32_t PsiAddOperatorPsiSerialMR(const OperatorType O, const WaveFunctionType psi, WaveFunctionType * const Opsi)
{
  return 0;
}

uint32_t OperatorPsiSerialCC(const OperatorType O, const WaveFunctionType psi, WaveFunctionType * const Opsi)
{
  return 0;
}

uint32_t OperatorPsiSerialCR(const OperatorType O, const WaveFunctionType psi, WaveFunctionType * const Opsi)
{
  return 0;
}

uint32_t OperatorPsiSerialRC(const OperatorType O, const WaveFunctionType psi, WaveFunctionType * const Opsi)
{
  return 0;
}

uint32_t OperatorPsiSerialRR(const OperatorType O, const WaveFunctionType psi, WaveFunctionType * const Opsi)
{
  return 0;
}

uint32_t OperatorPsiSerialMC(const OperatorType O, const WaveFunctionType psi, WaveFunctionType * const Opsi)
{
  return 0;
}

uint32_t OperatorPsiSerialMR(const OperatorType O, const WaveFunctionType psi, WaveFunctionType * const Opsi)
{
  return 0;
}

uint32_t PsiOperatorPsiCC(const OperatorType O, const WaveFunctionType psi, double *VRe, double *VIm)
{
  uint32_t i, j, k, l, ii, ij, ik;
  double phase;
  uint8_t *Det;
  int    iOMP;
  double ExpecValueRe = 0;
  double ExpecValueIm = 0;
  double aRe, aIm;
#pragma omp parallel private(aRe,aIm,Det,i,j,k,l,ii,ij,ik,phase,iOMP) reduction(+ : ExpecValueRe, ExpecValueIm)
  {
    Det = malloc(psi.NWord * sizeof(uint8_t));
    if (Det == NULL)
      AllocFailed();
#pragma omp for
    for (iOMP = 1; iOMP <= (int)psi.N; iOMP++)
    {
      i = (uint32_t) iOMP;
      ii = i >> NBitsBlok;
      ik = i & BlokMask;
      ij = psi.NWord * ik;
      for (l = 0; l <= O.MaxLength; l++)
      {
        for (j = 0; j < O.l[l].N; j++)
        {
          for (k = 0; k < psi.NWord; k++)
            Det[k] = psi.Determinant[ii][ij + k];
          phase = OperatorLader(Det, &(O.l[l].o[l * j]), l, O.NFermion);
          if (phase)
          {
            ComplexWaveFunctionPrefactorDet(psi, Det, &aRe, &aIm);
            ExpecValueRe +=
              phase * (aRe * O.l[l].V[j].Re * psi.ValueRe[ii][ik] - aRe * O.l[l].V[j].Im * psi.ValueIm[ii][ik] +
                       aIm * O.l[l].V[j].Re * psi.ValueIm[ii][ik] + aIm * O.l[l].V[j].Im * psi.ValueRe[ii][ik]);
            ExpecValueIm +=
              phase * (aIm * O.l[l].V[j].Im * psi.ValueIm[ii][ik] - aIm * O.l[l].V[j].Re * psi.ValueRe[ii][ik] +
                       aRe * O.l[l].V[j].Im * psi.ValueRe[ii][ik] + aRe * O.l[l].V[j].Re * psi.ValueIm[ii][ik]);
          }
        }
      }
    }
    free(Det);
  }
  *VRe = ExpecValueRe;
  *VIm = ExpecValueIm;
  return 0;
}

uint32_t PsiOperatorPsiCR(const OperatorType O, const WaveFunctionType psi, double *VRe, double *VIm)
{
  uint32_t i, j, k, l, ii, ij, ik;
  double phase;
  uint8_t *Det;
  int    iOMP;
  double ExpecValueRe = 0;
  double ExpecValueIm = 0;
#pragma omp parallel private(Det,i,j,k,l,ii,ij,ik,phase,iOMP) reduction(+ : ExpecValueRe, ExpecValueIm)
  {
    Det = malloc(psi.NWord * sizeof(uint8_t));
    if (Det == NULL)
      AllocFailed();
#pragma omp for
    for (iOMP = 1; iOMP <= (int)psi.N; iOMP++)
    {
      i = (uint32_t) iOMP;
      ii = i >> NBitsBlok;
      ik = i & BlokMask;
      ij = psi.NWord * ik;
      for (l = 0; l <= O.MaxLength; l++)
      {
        for (j = 0; j < O.l[l].N; j++)
        {
          for (k = 0; k < psi.NWord; k++)
            Det[k] = psi.Determinant[ii][ij + k];
          phase = OperatorLader(Det, &(O.l[l].o[l * j]), l, O.NFermion);
          if (phase)
          {
            ExpecValueRe += phase * RealWaveFunctionPrefactorDet(psi, Det) * O.l[l].V[j].Re * psi.ValueRe[ii][ik];
            ExpecValueIm += phase * RealWaveFunctionPrefactorDet(psi, Det) * O.l[l].V[j].Im * psi.ValueRe[ii][ik];
          }
        }
      }
    }
    free(Det);
  }
  *VRe = ExpecValueRe;
  *VIm = ExpecValueIm;
  return 0;
}

uint32_t PsiOperatorPsiRC(const OperatorType O, const WaveFunctionType psi, double *VRe, double *VIm)
{
  uint32_t i, j, k, l, ii, ij, ik;
  double phase;
  uint8_t *Det;
  int    iOMP;
  double ExpecValueRe = 0;
  double ExpecValueIm = 0;
  double aRe, aIm;
#pragma omp parallel private(aRe,aIm,Det,i,j,k,l,ii,ij,ik,phase,iOMP) reduction(+ : ExpecValueRe, ExpecValueIm)
  {
    Det = malloc(psi.NWord * sizeof(uint8_t));
    if (Det == NULL)
      AllocFailed();
#pragma omp for
    for (iOMP = 1; iOMP <= (int)psi.N; iOMP++)
    {
      i = (uint32_t) iOMP;
      ii = i >> NBitsBlok;
      ik = i & BlokMask;
      ij = psi.NWord * ik;
      for (l = 0; l <= O.MaxLength; l++)
      {
        for (j = 0; j < O.l[l].N; j++)
        {
          for (k = 0; k < psi.NWord; k++)
            Det[k] = psi.Determinant[ii][ij + k];
          phase = OperatorLader(Det, &(O.l[l].o[l * j]), l, O.NFermion);
          if (phase)
          {
            ComplexWaveFunctionPrefactorDet(psi, Det, &aRe, &aIm);
            ExpecValueRe +=
              phase * (aRe * O.l[l].VRe[j] * psi.ValueRe[ii][ik] + aIm * O.l[l].VRe[j] * psi.ValueIm[ii][ik]);
            ExpecValueIm +=
              phase * (-aIm * O.l[l].VRe[j] * psi.ValueRe[ii][ik] + aRe * O.l[l].VRe[j] * psi.ValueIm[ii][ik]);
          }
        }
      }
    }
    free(Det);
  }
  *VRe = ExpecValueRe;
  *VIm = ExpecValueIm;
  return 0;
}

uint32_t PsiOperatorPsiRR(const OperatorType O, const WaveFunctionType psi, double *VRe)
{
  uint32_t i, j, k, l, ii, ij, ik;
  double phase;
  uint8_t *Det;
  int    iOMP;
  double ExpecValue = 0;
  volatile double c, t, y;
  c = 0;
#pragma omp parallel firstprivate(c) private(Det,i,j,k,l,ii,ij,ik,phase,iOMP,t,y) reduction(+ : ExpecValue)
  {
    Det = malloc(psi.NWord * sizeof(uint8_t));
    if (Det == NULL)
      AllocFailed();
#pragma omp for
    for (iOMP = 1; iOMP <= (int)psi.N; iOMP++)
    {
      i = (uint32_t) iOMP;
      ii = i >> NBitsBlok;
      ik = i & BlokMask;
      ij = psi.NWord * ik;
      for (l = 0; l <= O.MaxLength; l++)
      {
        for (j = 0; j < O.l[l].N; j++)
        {
          for (k = 0; k < psi.NWord; k++)
            Det[k] = psi.Determinant[ii][ij + k];
          phase = OperatorLader(Det, &(O.l[l].o[l * j]), l, O.NFermion);
          if (phase)
          {
            y = (phase * RealWaveFunctionPrefactorDet(psi, Det) * O.l[l].VRe[j] * psi.ValueRe[ii][ik]) - c;
            t = ExpecValue + y;
            c = t - ExpecValue;
            c -= y;
            ExpecValue = t;
          }
        }
      }
    }
    free(Det);
  }
  *VRe = ExpecValue;
  return 0;
}

uint32_t PsiOperatorPsiMC(const OperatorType O, const WaveFunctionType psi, double *VRe, double *VIm)
{
  uint32_t i, j, k, l, ii, ij, ik;
  double phase;
  uint8_t *Det;
  int    iOMP;
  double ExpecValueRe = 0;
  double ExpecValueIm = 0;
  double aRe, aIm;
#pragma omp parallel private(aRe,aIm,Det,i,j,k,l,ii,ij,ik,phase,iOMP) reduction(+ : ExpecValueRe, ExpecValueIm)
  {
    Det = malloc(psi.NWord * sizeof(uint8_t));
    if (Det == NULL)
      AllocFailed();
#pragma omp for
    for (iOMP = 1; iOMP <= (int)psi.N; iOMP++)
    {
      i = (uint32_t) iOMP;
      ii = i >> NBitsBlok;
      ik = i & BlokMask;
      ij = psi.NWord * ik;
      for (l = 0; l <= O.MaxLength; l++)
      {
        if (O.l[l].QComplex == 0)
        {
          for (j = 0; j < O.l[l].N; j++)
          {
            for (k = 0; k < psi.NWord; k++)
              Det[k] = psi.Determinant[ii][ij + k];
            phase = OperatorLader(Det, &(O.l[l].o[l * j]), l, O.NFermion);
            if (phase)
            {
              ComplexWaveFunctionPrefactorDet(psi, Det, &aRe, &aIm);
              ExpecValueRe +=
                phase * (aRe * O.l[l].VRe[j] * psi.ValueRe[ii][ik] + aIm * O.l[l].VRe[j] * psi.ValueIm[ii][ik]);
              ExpecValueIm +=
                phase * (-aIm * O.l[l].VRe[j] * psi.ValueRe[ii][ik] + aRe * O.l[l].VRe[j] * psi.ValueIm[ii][ik]);
            }
          }
        }
        else
        {
          for (j = 0; j < O.l[l].N; j++)
          {
            for (k = 0; k < psi.NWord; k++)
              Det[k] = psi.Determinant[ii][ij + k];
            phase = OperatorLader(Det, &(O.l[l].o[l * j]), l, O.NFermion);
            if (phase)
            {
              ComplexWaveFunctionPrefactorDet(psi, Det, &aRe, &aIm);
              ExpecValueRe +=
                phase * (aRe * O.l[l].V[j].Re * psi.ValueRe[ii][ik] - aRe * O.l[l].V[j].Im * psi.ValueIm[ii][ik] +
                         aIm * O.l[l].V[j].Re * psi.ValueIm[ii][ik] + aIm * O.l[l].V[j].Im * psi.ValueRe[ii][ik]);
              ExpecValueIm +=
                phase * (aIm * O.l[l].V[j].Im * psi.ValueIm[ii][ik] - aIm * O.l[l].V[j].Re * psi.ValueRe[ii][ik] +
                         aRe * O.l[l].V[j].Im * psi.ValueRe[ii][ik] + aRe * O.l[l].V[j].Re * psi.ValueIm[ii][ik]);
            }
          }
        }
      }
    }
    free(Det);
  }
  *VRe = ExpecValueRe;
  *VIm = ExpecValueIm;
  return 0;
}

uint32_t PsiOperatorPsiMR(const OperatorType O, const WaveFunctionType psi, double *VRe, double *VIm)
{
  uint32_t i, j, k, l, ii, ij, ik;
  double phase;
  uint8_t *Det;
  int    iOMP;
  double ExpecValueRe = 0;
  double ExpecValueIm = 0;
#pragma omp parallel private(Det,i,j,k,l,ii,ij,ik,phase,iOMP) reduction(+ : ExpecValueRe, ExpecValueIm)
  {
    Det = malloc(psi.NWord * sizeof(uint8_t));
    if (Det == NULL)
      AllocFailed();
#pragma omp for
    for (iOMP = 1; iOMP <= (int)psi.N; iOMP++)
    {
      i = (uint32_t) iOMP;
      ii = i >> NBitsBlok;
      ik = i & BlokMask;
      ij = psi.NWord * ik;
      for (l = 0; l <= O.MaxLength; l++)
      {
        if (O.l[l].QComplex == 0)
        {
          for (j = 0; j < O.l[l].N; j++)
          {
            for (k = 0; k < psi.NWord; k++)
              Det[k] = psi.Determinant[ii][ij + k];
            phase = OperatorLader(Det, &(O.l[l].o[l * j]), l, O.NFermion);
            if (phase)
              ExpecValueRe += phase * RealWaveFunctionPrefactorDet(psi, Det) * O.l[l].VRe[j] * psi.ValueRe[ii][ik];
          }
        }
        else
        {
          for (j = 0; j < O.l[l].N; j++)
          {
            for (k = 0; k < psi.NWord; k++)
              Det[k] = psi.Determinant[ii][ij + k];
            phase = OperatorLader(Det, &(O.l[l].o[l * j]), l, O.NFermion);
            if (phase)
            {
              ExpecValueRe += phase * RealWaveFunctionPrefactorDet(psi, Det) * O.l[l].V[j].Re * psi.ValueRe[ii][ik];
              ExpecValueIm += phase * RealWaveFunctionPrefactorDet(psi, Det) * O.l[l].V[j].Im * psi.ValueRe[ii][ik];
            }
          }
        }
      }
    }
    free(Det);
  }
  *VRe = ExpecValueRe;
  *VIm = ExpecValueIm;
  return 0;
}

uint32_t PsiOperatorPsiCHC(const OperatorType O, const WaveFunctionType psi, double *VRe)
{
  uint32_t i, j, k, l, ii, ij, ik;
  double phase;
  uint8_t *Det;
  int    iOMP;
  double ExpecValueRe = 0;
  double aRe, aIm;
#pragma omp parallel private(aRe,aIm,Det,i,j,k,l,ii,ij,ik,phase,iOMP) reduction(+ : ExpecValueRe)
  {
    Det = malloc(psi.NWord * sizeof(uint8_t));
    if (Det == NULL)
      AllocFailed();
#pragma omp for
    for (iOMP = 1; iOMP <= (int)psi.N; iOMP++)
    {
      i = (uint32_t) iOMP;
      ii = i >> NBitsBlok;
      ik = i & BlokMask;
      ij = psi.NWord * ik;
      for (l = 0; l <= O.MaxLength; l++)
      {
        for (j = 0; j < O.l[l].N; j++)
        {
          for (k = 0; k < psi.NWord; k++)
            Det[k] = psi.Determinant[ii][ij + k];
          phase = OperatorLader(Det, &(O.l[l].o[l * j]), l, O.NFermion);
          if (phase)
          {
            ComplexWaveFunctionPrefactorDet(psi, Det, &aRe, &aIm);
            ExpecValueRe +=
              phase * (aRe * O.l[l].V[j].Re * psi.ValueRe[ii][ik] - aRe * O.l[l].V[j].Im * psi.ValueIm[ii][ik] +
                       aIm * O.l[l].V[j].Re * psi.ValueIm[ii][ik] + aIm * O.l[l].V[j].Im * psi.ValueRe[ii][ik]);
          }
        }
      }
    }
    free(Det);
  }
  *VRe = ExpecValueRe;
  return 0;
}

uint32_t PsiOperatorPsiCHR(const OperatorType O, const WaveFunctionType psi, double *VRe)
{
  uint32_t i, j, k, l, ii, ij, ik;
  double phase;
  uint8_t *Det;
  int    iOMP;
  double ExpecValueRe = 0;
#pragma omp parallel private(Det,i,j,k,l,ii,ij,ik,phase,iOMP) reduction(+ : ExpecValueRe)
  {
    Det = malloc(psi.NWord * sizeof(uint8_t));
    if (Det == NULL)
      AllocFailed();
#pragma omp for
    for (iOMP = 1; iOMP <= (int)psi.N; iOMP++)
    {
      i = (uint32_t) iOMP;
      ii = i >> NBitsBlok;
      ik = i & BlokMask;
      ij = psi.NWord * ik;
      for (l = 0; l <= O.MaxLength; l++)
      {
        for (j = 0; j < O.l[l].N; j++)
        {
          for (k = 0; k < psi.NWord; k++)
            Det[k] = psi.Determinant[ii][ij + k];
          phase = OperatorLader(Det, &(O.l[l].o[l * j]), l, O.NFermion);
          if (phase)
          {
            ExpecValueRe += phase * RealWaveFunctionPrefactorDet(psi, Det) * O.l[l].V[j].Re * psi.ValueRe[ii][ik];
          }
        }
      }
    }
    free(Det);
  }
  *VRe = ExpecValueRe;
  return 0;
}

uint32_t PsiOperatorPsiRHC(const OperatorType O, const WaveFunctionType psi, double *VRe)
{
  uint32_t i, j, k, l, ii, ij, ik;
  double phase;
  uint8_t *Det;
  int    iOMP;
  double ExpecValueRe = 0;
  double aRe, aIm;
#pragma omp parallel private(aRe,aIm,Det,i,j,k,l,ii,ij,ik,phase,iOMP) reduction(+ : ExpecValueRe)
  {
    Det = malloc(psi.NWord * sizeof(uint8_t));
    if (Det == NULL)
      AllocFailed();
#pragma omp for
    for (iOMP = 1; iOMP <= (int)psi.N; iOMP++)
    {
      i = (uint32_t) iOMP;
      ii = i >> NBitsBlok;
      ik = i & BlokMask;
      ij = psi.NWord * ik;
      for (l = 0; l <= O.MaxLength; l++)
      {
        for (j = 0; j < O.l[l].N; j++)
        {
          for (k = 0; k < psi.NWord; k++)
            Det[k] = psi.Determinant[ii][ij + k];
          phase = OperatorLader(Det, &(O.l[l].o[l * j]), l, O.NFermion);
          if (phase)
          {
            ComplexWaveFunctionPrefactorDet(psi, Det, &aRe, &aIm);
            ExpecValueRe +=
              phase * (aRe * O.l[l].VRe[j] * psi.ValueRe[ii][ik] + aIm * O.l[l].VRe[j] * psi.ValueIm[ii][ik]);
          }
        }
      }
    }
    free(Det);
  }
  *VRe = ExpecValueRe;
  return 0;
}

uint32_t PsiOperatorPsiRHR(const OperatorType O, const WaveFunctionType psi, double *VRe)
{
  uint32_t i, j, k, l, ii, ij, ik;
  double phase;
  uint8_t *Det;
  int    iOMP;
  double ExpecValue = 0;
#pragma omp parallel private(Det,i,j,k,l,ii,ij,ik,phase,iOMP) reduction(+ : ExpecValue)
  {
    Det = malloc(psi.NWord * sizeof(uint8_t));
    if (Det == NULL)
      AllocFailed();
#pragma omp for
    for (iOMP = 1; iOMP <= (int)psi.N; iOMP++)
    {
      i = (uint32_t) iOMP;
      ii = i >> NBitsBlok;
      ik = i & BlokMask;
      ij = psi.NWord * ik;
      for (l = 0; l <= O.MaxLength; l++)
      {
        for (j = 0; j < O.l[l].N; j++)
        {
          for (k = 0; k < psi.NWord; k++)
            Det[k] = psi.Determinant[ii][ij + k];
          phase = OperatorLader(Det, &(O.l[l].o[l * j]), l, O.NFermion);
          if (phase)
            ExpecValue += phase * RealWaveFunctionPrefactorDet(psi, Det) * O.l[l].VRe[j] * psi.ValueRe[ii][ik];
        }
      }
    }
    free(Det);
  }
  *VRe = ExpecValue;
  return 0;
}

uint32_t PsiOperatorPsiMHC(const OperatorType O, const WaveFunctionType psi, double *VRe)
{
  uint32_t i, j, k, l, ii, ij, ik;
  double phase;
  uint8_t *Det;
  int    iOMP;
  double ExpecValueRe = 0;
  double aRe, aIm;
#pragma omp parallel private(aRe,aIm,Det,i,j,k,l,ii,ij,ik,phase,iOMP) reduction(+ : ExpecValueRe)
  {
    Det = malloc(psi.NWord * sizeof(uint8_t));
    if (Det == NULL)
      AllocFailed();
#pragma omp for
    for (iOMP = 1; iOMP <= (int)psi.N; iOMP++)
    {
      i = (uint32_t) iOMP;
      ii = i >> NBitsBlok;
      ik = i & BlokMask;
      ij = psi.NWord * ik;
      for (l = 0; l <= O.MaxLength; l++)
      {
        if (O.l[l].QComplex == 0)
        {
          for (j = 0; j < O.l[l].N; j++)
          {
            for (k = 0; k < psi.NWord; k++)
              Det[k] = psi.Determinant[ii][ij + k];
            phase = OperatorLader(Det, &(O.l[l].o[l * j]), l, O.NFermion);
            if (phase)
            {
              ComplexWaveFunctionPrefactorDet(psi, Det, &aRe, &aIm);
              ExpecValueRe +=
                phase * (aRe * O.l[l].VRe[j] * psi.ValueRe[ii][ik] + aIm * O.l[l].VRe[j] * psi.ValueIm[ii][ik]);
            }
          }
        }
        else
        {
          for (j = 0; j < O.l[l].N; j++)
          {
            for (k = 0; k < psi.NWord; k++)
              Det[k] = psi.Determinant[ii][ij + k];
            phase = OperatorLader(Det, &(O.l[l].o[l * j]), l, O.NFermion);
            if (phase)
            {
              ComplexWaveFunctionPrefactorDet(psi, Det, &aRe, &aIm);
              ExpecValueRe +=
                phase * (aRe * O.l[l].V[j].Re * psi.ValueRe[ii][ik] - aRe * O.l[l].V[j].Im * psi.ValueIm[ii][ik] +
                         aIm * O.l[l].V[j].Re * psi.ValueIm[ii][ik] + aIm * O.l[l].V[j].Im * psi.ValueRe[ii][ik]);
            }
          }
        }
      }
    }
    free(Det);
  }
  *VRe = ExpecValueRe;
  return 0;
}

uint32_t PsiOperatorPsiMHR(const OperatorType O, const WaveFunctionType psi, double *VRe)
{
  uint32_t i, j, k, l, ii, ij, ik;
  double phase;
  uint8_t *Det;
  int    iOMP;
  double ExpecValueRe = 0;
#pragma omp parallel private(Det,i,j,k,l,ii,ij,ik,phase,iOMP) reduction(+ : ExpecValueRe)
  {
    Det = malloc(psi.NWord * sizeof(uint8_t));
    if (Det == NULL)
      AllocFailed();
#pragma omp for
    for (iOMP = 1; iOMP <= (int)psi.N; iOMP++)
    {
      i = (uint32_t) iOMP;
      ii = i >> NBitsBlok;
      ik = i & BlokMask;
      ij = psi.NWord * ik;
      for (l = 0; l <= O.MaxLength; l++)
      {
        if (O.l[l].QComplex == 0)
        {
          for (j = 0; j < O.l[l].N; j++)
          {
            for (k = 0; k < psi.NWord; k++)
              Det[k] = psi.Determinant[ii][ij + k];
            phase = OperatorLader(Det, &(O.l[l].o[l * j]), l, O.NFermion);
            if (phase)
              ExpecValueRe += phase * RealWaveFunctionPrefactorDet(psi, Det) * O.l[l].VRe[j] * psi.ValueRe[ii][ik];
          }
        }
        else
        {
          for (j = 0; j < O.l[l].N; j++)
          {
            for (k = 0; k < psi.NWord; k++)
              Det[k] = psi.Determinant[ii][ij + k];
            phase = OperatorLader(Det, &(O.l[l].o[l * j]), l, O.NFermion);
            if (phase)
            {
              ExpecValueRe += phase * RealWaveFunctionPrefactorDet(psi, Det) * O.l[l].V[j].Re * psi.ValueRe[ii][ik];
            }
          }
        }
      }
    }
    free(Det);
  }
  *VRe = ExpecValueRe;
  return 0;
}

uint32_t PsiOperatorSqrPsiCC(const OperatorType O, const WaveFunctionType psi, double *VRe, double *VIm)
{
  uint32_t i, j, k, l, ii, ij, ik, lS, jS;
  double phase, phaseS;
  uint8_t *Det, *DetSqr;
  int    iOMP;
  double ExpecValueRe = 0;
  double ExpecValueIm = 0;
  double aRe, aIm, temp;
#pragma omp parallel private(temp,aRe,aIm,Det,DetSqr,i,j,k,l,lS,jS,ii,ij,ik,phase,phaseS,iOMP) reduction(+ : ExpecValueRe, ExpecValueIm)
  {
    Det = malloc(psi.NWord * sizeof(uint8_t));
    if (Det == NULL)
      AllocFailed();
    DetSqr = malloc(psi.NWord * sizeof(uint8_t));
    if (DetSqr == NULL)
      AllocFailed();
#pragma omp for
    for (iOMP = 1; iOMP <= (int)psi.N; iOMP++)
    {
      i = (uint32_t) iOMP;
      ii = i >> NBitsBlok;
      ik = i & BlokMask;
      ij = psi.NWord * ik;
      for (l = 0; l <= O.MaxLength; l++)
      {
        for (j = 0; j < O.l[l].N; j++)
        {
          for (k = 0; k < psi.NWord; k++)
            Det[k] = psi.Determinant[ii][ij + k];
          phase = OperatorLader(Det, &(O.l[l].o[l * j]), l, O.NFermion);
          if (phase)
            for (lS = 0; lS <= O.MaxLength; lS++)
            {
              for (jS = 0; jS < O.l[lS].N; jS++)
              {
                for (k = 0; k < psi.NWord; k++)
                  DetSqr[k] = Det[k];
                phaseS = OperatorLader(DetSqr, &(O.l[lS].o[lS * jS]), lS, O.NFermion);
                if (phaseS)
                {
                  ComplexWaveFunctionPrefactorDet(psi, DetSqr, &aRe, &aIm);
                  aIm *= -1;
                  temp = aRe;
                  aRe = temp * O.l[l].V[j].Re - aIm * O.l[l].V[j].Im;
                  aIm = temp * O.l[l].V[j].Im + aIm * O.l[l].V[j].Re;
                  temp = aRe;
                  aRe = temp * O.l[lS].V[jS].Re - aIm * O.l[lS].V[jS].Im;
                  aIm = temp * O.l[lS].V[jS].Im + aIm * O.l[lS].V[jS].Re;
                  temp = aRe;
                  aRe = temp * psi.ValueRe[ii][ik] - aIm * psi.ValueIm[ii][ik];
                  aIm = temp * psi.ValueIm[ii][ik] + aIm * psi.ValueRe[ii][ik];
                  ExpecValueRe += phase * phaseS * aRe;
                  ExpecValueIm += phase * phaseS * aIm;
                }
              }
            }
        }
      }
    }
    free(Det);
    free(DetSqr);
  }
  *VRe = ExpecValueRe;
  *VIm = ExpecValueIm;
  return 0;
}

uint32_t PsiOperatorSqrPsiCR(const OperatorType O, const WaveFunctionType psi, double *VRe, double *VIm)
{
  uint32_t i, j, k, l, ii, ij, ik, lS, jS;
  double phase, phaseS;
  uint8_t *Det, *DetSqr;
  int    iOMP;
  double ExpecValueRe = 0;
  double ExpecValueIm = 0;
#pragma omp parallel private(Det,DetSqr,i,j,k,l,lS,jS,ii,ij,ik,phase,phaseS,iOMP) reduction(+ : ExpecValueRe, ExpecValueIm)
  {
    Det = malloc(psi.NWord * sizeof(uint8_t));
    if (Det == NULL)
      AllocFailed();
    DetSqr = malloc(psi.NWord * sizeof(uint8_t));
    if (DetSqr == NULL)
      AllocFailed();
#pragma omp for
    for (iOMP = 1; iOMP <= (int)psi.N; iOMP++)
    {
      i = (uint32_t) iOMP;
      ii = i >> NBitsBlok;
      ik = i & BlokMask;
      ij = psi.NWord * ik;
      for (l = 0; l <= O.MaxLength; l++)
      {
        for (j = 0; j < O.l[l].N; j++)
        {
          for (k = 0; k < psi.NWord; k++)
            Det[k] = psi.Determinant[ii][ij + k];
          phase = OperatorLader(Det, &(O.l[l].o[l * j]), l, O.NFermion);
          if (phase)
            for (lS = 0; lS <= O.MaxLength; lS++)
            {
              for (jS = 0; jS < O.l[lS].N; jS++)
              {
                for (k = 0; k < psi.NWord; k++)
                  DetSqr[k] = Det[k];
                phaseS = OperatorLader(DetSqr, &(O.l[lS].o[lS * jS]), lS, O.NFermion);
                if (phaseS)
                {
                  ExpecValueRe +=
                    phase * phaseS * RealWaveFunctionPrefactorDet(psi,
                                                                  Det) * psi.ValueRe[ii][ik] * (O.l[l].V[j].Re *
                                                                                                O.l[lS].V[jS].Re -
                                                                                                O.l[l].V[j].Im *
                                                                                                O.l[lS].V[jS].Im);
                  ExpecValueIm +=
                    phase * phaseS * RealWaveFunctionPrefactorDet(psi,
                                                                  Det) * psi.ValueRe[ii][ik] * (O.l[l].V[j].Re *
                                                                                                O.l[lS].V[jS].Im +
                                                                                                O.l[l].V[j].Im *
                                                                                                O.l[lS].V[jS].Re);
                }
              }
            }
        }
      }
    }
    free(Det);
    free(DetSqr);
  }
  *VRe = ExpecValueRe;
  *VIm = ExpecValueIm;
  return 0;
}

uint32_t PsiOperatorSqrPsiRC(const OperatorType O, const WaveFunctionType psi, double *VRe, double *VIm)
{
  uint32_t i, j, k, l, ii, ij, ik, lS, jS;
  double phase, phaseS;
  uint8_t *Det, *DetSqr;
  int    iOMP;
  double ExpecValueRe = 0;
  double ExpecValueIm = 0;
  double aRe, aIm;
#pragma omp parallel private(aRe,aIm,Det,DetSqr,i,j,k,l,lS,jS,ii,ij,ik,phase,phaseS,iOMP) reduction(+ : ExpecValueRe, ExpecValueIm)
  {
    Det = malloc(psi.NWord * sizeof(uint8_t));
    if (Det == NULL)
      AllocFailed();
    DetSqr = malloc(psi.NWord * sizeof(uint8_t));
    if (DetSqr == NULL)
      AllocFailed();
#pragma omp for
    for (iOMP = 1; iOMP <= (int)psi.N; iOMP++)
    {
      i = (uint32_t) iOMP;
      ii = i >> NBitsBlok;
      ik = i & BlokMask;
      ij = psi.NWord * ik;
      for (l = 0; l <= O.MaxLength; l++)
      {
        for (j = 0; j < O.l[l].N; j++)
        {
          for (k = 0; k < psi.NWord; k++)
            Det[k] = psi.Determinant[ii][ij + k];
          phase = OperatorLader(Det, &(O.l[l].o[l * j]), l, O.NFermion);
          if (phase)
            for (lS = 0; lS <= O.MaxLength; lS++)
            {
              for (jS = 0; jS < O.l[lS].N; jS++)
              {
                for (k = 0; k < psi.NWord; k++)
                  DetSqr[k] = Det[k];
                phaseS = OperatorLader(DetSqr, &(O.l[lS].o[lS * jS]), lS, O.NFermion);
                if (phaseS)
                {
                  ComplexWaveFunctionPrefactorDet(psi, DetSqr, &aRe, &aIm);
                  ExpecValueRe +=
                    phase * phaseS * (aRe * O.l[l].VRe[j] * O.l[lS].VRe[jS] * psi.ValueRe[ii][ik] +
                                      aIm * O.l[l].VRe[j] * O.l[lS].VRe[jS] * psi.ValueIm[ii][ik]);
                  ExpecValueIm +=
                    phase * phaseS * (-aIm * O.l[l].VRe[j] * O.l[lS].VRe[jS] * psi.ValueRe[ii][ik] +
                                      aRe * O.l[l].VRe[j] * O.l[lS].VRe[jS] * psi.ValueIm[ii][ik]);
                }
              }
            }
        }
      }
    }
    free(Det);
    free(DetSqr);
  }
  *VRe = ExpecValueRe;
  *VIm = ExpecValueIm;
  return 0;
}

uint32_t PsiOperatorSqrPsiRR(const OperatorType O, const WaveFunctionType psi, double *VRe)
{
  uint32_t i, j, k, l, ii, ij, ik, lS, jS;
  double phase, phaseS;
  uint8_t *Det, *DetSqr;
  int    iOMP;
  double ExpecValue = 0;
#pragma omp parallel private(Det,DetSqr,i,j,k,l,lS,jS,ii,ij,ik,phase,phaseS,iOMP) reduction(+ : ExpecValue)
  {
    Det = malloc(psi.NWord * sizeof(uint8_t));
    if (Det == NULL)
      AllocFailed();
    DetSqr = malloc(psi.NWord * sizeof(uint8_t));
    if (DetSqr == NULL)
      AllocFailed();
#pragma omp for
    for (iOMP = 1; iOMP <= (int)psi.N; iOMP++)
    {
      i = (uint32_t) iOMP;
      ii = i >> NBitsBlok;
      ik = i & BlokMask;
      ij = psi.NWord * ik;
      for (l = 0; l <= O.MaxLength; l++)
      {
        for (j = 0; j < O.l[l].N; j++)
        {
          for (k = 0; k < psi.NWord; k++)
            Det[k] = psi.Determinant[ii][ij + k];
          phase = OperatorLader(Det, &(O.l[l].o[l * j]), l, O.NFermion);
          if (phase)
            for (lS = 0; lS <= O.MaxLength; lS++)
            {
              for (jS = 0; jS < O.l[lS].N; jS++)
              {
                for (k = 0; k < psi.NWord; k++)
                  DetSqr[k] = Det[k];
                phaseS = OperatorLader(DetSqr, &(O.l[lS].o[lS * jS]), lS, O.NFermion);
                if (phaseS)
                  ExpecValue +=
                    phase * phaseS * RealWaveFunctionPrefactorDet(psi,
                                                                  DetSqr) * O.l[l].VRe[j] * O.l[lS].VRe[jS] *
                    psi.ValueRe[ii][ik];
              }
            }
        }
      }
    }
    free(Det);
    free(DetSqr);
  }
  *VRe = ExpecValue;
  return 0;
}

uint32_t PsiOperatorSqrPsiMC(const OperatorType O, const WaveFunctionType psi, double *VRe, double *VIm)
{
  uint32_t i, j, k, l, ii, ij, ik, lS, jS;
  double phase, phaseS;
  uint8_t *Det, *DetSqr;
  int    iOMP;
  double ExpecValueRe = 0;
  double ExpecValueIm = 0;
  double aRe, aIm, temp;
#pragma omp parallel private(temp,aRe,aIm,Det,DetSqr,i,j,k,l,lS,jS,ii,ij,ik,phase,phaseS,iOMP) reduction(+ : ExpecValueRe, ExpecValueIm)
  {
    Det = malloc(psi.NWord * sizeof(uint8_t));
    if (Det == NULL)
      AllocFailed();
    DetSqr = malloc(psi.NWord * sizeof(uint8_t));
    if (DetSqr == NULL)
      AllocFailed();
#pragma omp for
    for (iOMP = 1; iOMP <= (int)psi.N; iOMP++)
    {
      i = (uint32_t) iOMP;
      ii = i >> NBitsBlok;
      ik = i & BlokMask;
      ij = psi.NWord * ik;
      for (l = 0; l <= O.MaxLength; l++)
      {
        if (O.l[l].QComplex == 0)
        {                       /* Operator for length l is real */
          for (j = 0; j < O.l[l].N; j++)
          {
            for (k = 0; k < psi.NWord; k++)
              Det[k] = psi.Determinant[ii][ij + k];
            phase = OperatorLader(Det, &(O.l[l].o[l * j]), l, O.NFermion);
            if (phase)
              for (lS = 0; lS <= O.MaxLength; lS++)
              {
                if (O.l[l].QComplex == 0)
                {               /* Operator for length lS is real */
                  for (jS = 0; jS < O.l[lS].N; jS++)
                  {
                    for (k = 0; k < psi.NWord; k++)
                      DetSqr[k] = Det[k];
                    phaseS = OperatorLader(DetSqr, &(O.l[lS].o[lS * jS]), lS, O.NFermion);
                    if (phaseS)
                    {
                      ComplexWaveFunctionPrefactorDet(psi, DetSqr, &aRe, &aIm);
                      ExpecValueRe +=
                        phase * phaseS * (aRe * O.l[l].VRe[j] * O.l[lS].VRe[jS] * psi.ValueRe[ii][ik] +
                                          aIm * O.l[l].VRe[j] * O.l[lS].VRe[jS] * psi.ValueIm[ii][ik]);
                      ExpecValueIm +=
                        phase * phaseS * (-aIm * O.l[l].VRe[j] * O.l[lS].VRe[jS] * psi.ValueRe[ii][ik] +
                                          aRe * O.l[l].VRe[j] * O.l[lS].VRe[jS] * psi.ValueIm[ii][ik]);
                    }
                  }
                }
                else
                {               /* Operator for length lS is complex */
                  for (jS = 0; jS < O.l[lS].N; jS++)
                  {
                    for (k = 0; k < psi.NWord; k++)
                      DetSqr[k] = Det[k];
                    phaseS = OperatorLader(DetSqr, &(O.l[lS].o[lS * jS]), lS, O.NFermion);
                    if (phaseS)
                    {
                      ComplexWaveFunctionPrefactorDet(psi, DetSqr, &aRe, &aIm);
                      aIm *= -1;
                      aRe = aRe * O.l[l].VRe[j];
                      aIm = aIm * O.l[l].VRe[j];
                      temp = aRe;
                      aRe = temp * O.l[lS].V[jS].Re - aIm * O.l[lS].V[jS].Im;
                      aIm = temp * O.l[lS].V[jS].Im + aIm * O.l[lS].V[jS].Re;
                      temp = aRe;
                      aRe = temp * psi.ValueRe[ii][ik] - aIm * psi.ValueIm[ii][ik];
                      aIm = temp * psi.ValueIm[ii][ik] + aIm * psi.ValueRe[ii][ik];
                      ExpecValueRe += phase * phaseS * aRe;
                      ExpecValueIm += phase * phaseS * aIm;
                    }
                  }
                }
              }
          }
        }
        else
        {                       /* Operator for length l is complex */
          for (j = 0; j < O.l[l].N; j++)
          {
            for (k = 0; k < psi.NWord; k++)
              Det[k] = psi.Determinant[ii][ij + k];
            phase = OperatorLader(Det, &(O.l[l].o[l * j]), l, O.NFermion);
            if (phase)
              for (lS = 0; lS <= O.MaxLength; lS++)
              {
                if (O.l[l].QComplex == 0)
                {               /* Operator for length lS is real */
                  for (jS = 0; jS < O.l[lS].N; jS++)
                  {
                    for (k = 0; k < psi.NWord; k++)
                      DetSqr[k] = Det[k];
                    phaseS = OperatorLader(DetSqr, &(O.l[lS].o[lS * jS]), lS, O.NFermion);
                    if (phaseS)
                    {
                      ComplexWaveFunctionPrefactorDet(psi, DetSqr, &aRe, &aIm);
                      aIm *= -1;
                      temp = aRe;
                      aRe = temp * O.l[l].V[j].Re - aIm * O.l[l].V[j].Im;
                      aIm = temp * O.l[l].V[j].Im + aIm * O.l[l].V[j].Re;
                      aRe = aRe * O.l[lS].VRe[jS];
                      aIm = aIm * O.l[lS].VRe[jS];
                      temp = aRe;
                      aRe = temp * psi.ValueRe[ii][ik] - aIm * psi.ValueIm[ii][ik];
                      aIm = temp * psi.ValueIm[ii][ik] + aIm * psi.ValueRe[ii][ik];
                      ExpecValueRe += phase * phaseS * aRe;
                      ExpecValueIm += phase * phaseS * aIm;
                    }
                  }
                }
                else
                {               /* Operator for length lS is complex */
                  for (jS = 0; jS < O.l[lS].N; jS++)
                  {
                    for (k = 0; k < psi.NWord; k++)
                      DetSqr[k] = Det[k];
                    phaseS = OperatorLader(DetSqr, &(O.l[lS].o[lS * jS]), lS, O.NFermion);
                    if (phaseS)
                    {
                      ComplexWaveFunctionPrefactorDet(psi, DetSqr, &aRe, &aIm);
                      aIm *= -1;
                      temp = aRe;
                      aRe = temp * O.l[l].V[j].Re - aIm * O.l[l].V[j].Im;
                      aIm = temp * O.l[l].V[j].Im + aIm * O.l[l].V[j].Re;
                      temp = aRe;
                      aRe = temp * O.l[lS].V[jS].Re - aIm * O.l[lS].V[jS].Im;
                      aIm = temp * O.l[lS].V[jS].Im + aIm * O.l[lS].V[jS].Re;
                      temp = aRe;
                      aRe = temp * psi.ValueRe[ii][ik] - aIm * psi.ValueIm[ii][ik];
                      aIm = temp * psi.ValueIm[ii][ik] + aIm * psi.ValueRe[ii][ik];
                      ExpecValueRe += phase * phaseS * aRe;
                      ExpecValueIm += phase * phaseS * aIm;
                    }
                  }
                }
              }
          }
        }
      }
    }
    free(Det);
    free(DetSqr);
  }
  *VRe = ExpecValueRe;
  *VIm = ExpecValueIm;
  return 0;
}

uint32_t PsiOperatorSqrPsiMR(const OperatorType O, const WaveFunctionType psi, double *VRe, double *VIm)
{
  uint32_t i, j, k, l, ii, ij, ik, lS, jS;
  double phase, phaseS;
  uint8_t *Det, *DetSqr;
  int    iOMP;
  double ExpecValueRe = 0;
  double ExpecValueIm = 0;
#pragma omp parallel private(Det,DetSqr,i,j,k,l,lS,jS,ii,ij,ik,phase,phaseS,iOMP) reduction(+ : ExpecValueRe, ExpecValueIm)
  {
    Det = malloc(psi.NWord * sizeof(uint8_t));
    if (Det == NULL)
      AllocFailed();
    DetSqr = malloc(psi.NWord * sizeof(uint8_t));
    if (DetSqr == NULL)
      AllocFailed();
#pragma omp for
    for (iOMP = 1; iOMP <= (int)psi.N; iOMP++)
    {
      i = (uint32_t) iOMP;
      ii = i >> NBitsBlok;
      ik = i & BlokMask;
      ij = psi.NWord * ik;
      for (l = 0; l <= O.MaxLength; l++)
      {
        if (O.l[l].QComplex == 0)
        {                       /* Operator at length l real */
          for (j = 0; j < O.l[l].N; j++)
          {
            for (k = 0; k < psi.NWord; k++)
              Det[k] = psi.Determinant[ii][ij + k];
            phase = OperatorLader(Det, &(O.l[l].o[l * j]), l, O.NFermion);
            if (phase)
              for (lS = 0; lS <= O.MaxLength; lS++)
              {
                if (O.l[lS].QComplex == 0)
                {               /* Operator at length lS real */
                  for (jS = 0; jS < O.l[lS].N; jS++)
                  {
                    for (k = 0; k < psi.NWord; k++)
                      DetSqr[k] = Det[k];
                    phaseS = OperatorLader(DetSqr, &(O.l[lS].o[lS * jS]), lS, O.NFermion);
                    if (phaseS)
                    {
                      ExpecValueRe +=
                        phase * phaseS * RealWaveFunctionPrefactorDet(psi,
                                                                      Det) * psi.ValueRe[ii][ik] * (O.l[l].VRe[j] *
                                                                                                    O.l[lS].VRe[jS]);
                    }
                  }
                }
                else
                {               /* Operator at length lS complex */
                  for (jS = 0; jS < O.l[lS].N; jS++)
                  {
                    for (k = 0; k < psi.NWord; k++)
                      DetSqr[k] = Det[k];
                    phaseS = OperatorLader(DetSqr, &(O.l[lS].o[lS * jS]), lS, O.NFermion);
                    if (phaseS)
                    {
                      ExpecValueRe +=
                        phase * phaseS * RealWaveFunctionPrefactorDet(psi,
                                                                      Det) * psi.ValueRe[ii][ik] * (O.l[l].VRe[j] *
                                                                                                    O.l[lS].V[jS].Re);
                      ExpecValueIm +=
                        phase * phaseS * RealWaveFunctionPrefactorDet(psi,
                                                                      Det) * psi.ValueRe[ii][ik] * (O.l[l].VRe[j] *
                                                                                                    O.l[lS].V[jS].Im);
                    }
                  }
                }
              }
          }
        }
        else
        {                       /* Operator at length l complex */
          for (j = 0; j < O.l[l].N; j++)
          {
            for (k = 0; k < psi.NWord; k++)
              Det[k] = psi.Determinant[ii][ij + k];
            phase = OperatorLader(Det, &(O.l[l].o[l * j]), l, O.NFermion);
            if (phase)
              for (lS = 0; lS <= O.MaxLength; lS++)
              {
                if (O.l[lS].QComplex == 0)
                {               /* Operator at length lS real */
                  for (jS = 0; jS < O.l[lS].N; jS++)
                  {
                    for (k = 0; k < psi.NWord; k++)
                      DetSqr[k] = Det[k];
                    phaseS = OperatorLader(DetSqr, &(O.l[lS].o[lS * jS]), lS, O.NFermion);
                    if (phaseS)
                    {
                      ExpecValueRe +=
                        phase * phaseS * RealWaveFunctionPrefactorDet(psi,
                                                                      Det) * psi.ValueRe[ii][ik] * (O.l[l].V[j].Re *
                                                                                                    O.l[lS].VRe[jS]);
                      ExpecValueIm +=
                        phase * phaseS * RealWaveFunctionPrefactorDet(psi,
                                                                      Det) * psi.ValueRe[ii][ik] * (O.l[l].V[j].Im *
                                                                                                    O.l[lS].VRe[jS]);
                    }
                  }
                }
                else
                {               /* Operator at length lS complex */
                  for (jS = 0; jS < O.l[lS].N; jS++)
                  {
                    for (k = 0; k < psi.NWord; k++)
                      DetSqr[k] = Det[k];
                    phaseS = OperatorLader(DetSqr, &(O.l[lS].o[lS * jS]), lS, O.NFermion);
                    if (phaseS)
                    {
                      ExpecValueRe +=
                        phase * phaseS * RealWaveFunctionPrefactorDet(psi,
                                                                      Det) * psi.ValueRe[ii][ik] * (O.l[l].V[j].Re *
                                                                                                    O.l[lS].V[jS].Re -
                                                                                                    O.l[l].V[j].Im *
                                                                                                    O.l[lS].V[jS].Im);
                      ExpecValueIm +=
                        phase * phaseS * RealWaveFunctionPrefactorDet(psi,
                                                                      Det) * psi.ValueRe[ii][ik] * (O.l[l].V[j].Re *
                                                                                                    O.l[lS].V[jS].Im +
                                                                                                    O.l[l].V[j].Im *
                                                                                                    O.l[lS].V[jS].Re);
                    }
                  }
                }
              }
          }
        }
      }
    }
    free(Det);
    free(DetSqr);
  }
  *VRe = ExpecValueRe;
  *VIm = ExpecValueIm;
  return 0;
}

uint32_t PsiOperatorSqrPsiCHC(const OperatorType O, const WaveFunctionType psi, double *VRe)
{
  uint32_t i, j, k, l, ii, ij, ik, lS, jS;
  double phase, phaseS;
  uint8_t *Det, *DetSqr;
  int    iOMP;
  double ExpecValueRe = 0;
  double aRe, aIm, temp;
#pragma omp parallel private(temp,aRe,aIm,Det,DetSqr,i,j,k,l,lS,jS,ii,ij,ik,phase,phaseS,iOMP) reduction(+ : ExpecValueRe)
  {
    Det = malloc(psi.NWord * sizeof(uint8_t));
    if (Det == NULL)
      AllocFailed();
    DetSqr = malloc(psi.NWord * sizeof(uint8_t));
    if (DetSqr == NULL)
      AllocFailed();
#pragma omp for
    for (iOMP = 1; iOMP <= (int)psi.N; iOMP++)
    {
      i = (uint32_t) iOMP;
      ii = i >> NBitsBlok;
      ik = i & BlokMask;
      ij = psi.NWord * ik;
      for (l = 0; l <= O.MaxLength; l++)
      {
        for (j = 0; j < O.l[l].N; j++)
        {
          for (k = 0; k < psi.NWord; k++)
            Det[k] = psi.Determinant[ii][ij + k];
          phase = OperatorLader(Det, &(O.l[l].o[l * j]), l, O.NFermion);
          if (phase)
            for (lS = 0; lS <= O.MaxLength; lS++)
            {
              for (jS = 0; jS < O.l[lS].N; jS++)
              {
                for (k = 0; k < psi.NWord; k++)
                  DetSqr[k] = Det[k];
                phaseS = OperatorLader(DetSqr, &(O.l[lS].o[lS * jS]), lS, O.NFermion);
                if (phaseS)
                {
                  ComplexWaveFunctionPrefactorDet(psi, DetSqr, &aRe, &aIm);
                  aIm *= -1;
                  temp = aRe;
                  aRe = temp * O.l[l].V[j].Re - aIm * O.l[l].V[j].Im;
                  aIm = temp * O.l[l].V[j].Im + aIm * O.l[l].V[j].Re;
                  temp = aRe;
                  aRe = temp * O.l[lS].V[jS].Re - aIm * O.l[lS].V[jS].Im;
                  aIm = temp * O.l[lS].V[jS].Im + aIm * O.l[lS].V[jS].Re;
                  temp = aRe;
                  aRe = temp * psi.ValueRe[ii][ik] - aIm * psi.ValueIm[ii][ik];
                  ExpecValueRe += phase * phaseS * aRe;
                }
              }
            }
        }
      }
    }
    free(Det);
    free(DetSqr);
  }
  *VRe = ExpecValueRe;
  return 0;
}

uint32_t PsiOperatorSqrPsiCHR(const OperatorType O, const WaveFunctionType psi, double *VRe)
{
  uint32_t i, j, k, l, ii, ij, ik, lS, jS;
  double phase, phaseS;
  uint8_t *Det, *DetSqr;
  int    iOMP;
  double ExpecValueRe = 0;
#pragma omp parallel private(Det,DetSqr,i,j,k,l,lS,jS,ii,ij,ik,phase,phaseS,iOMP) reduction(+ : ExpecValueRe)
  {
    Det = malloc(psi.NWord * sizeof(uint8_t));
    if (Det == NULL)
      AllocFailed();
    DetSqr = malloc(psi.NWord * sizeof(uint8_t));
    if (DetSqr == NULL)
      AllocFailed();
#pragma omp for
    for (iOMP = 1; iOMP <= (int)psi.N; iOMP++)
    {
      i = (uint32_t) iOMP;
      ii = i >> NBitsBlok;
      ik = i & BlokMask;
      ij = psi.NWord * ik;
      for (l = 0; l <= O.MaxLength; l++)
      {
        for (j = 0; j < O.l[l].N; j++)
        {
          for (k = 0; k < psi.NWord; k++)
            Det[k] = psi.Determinant[ii][ij + k];
          phase = OperatorLader(Det, &(O.l[l].o[l * j]), l, O.NFermion);
          if (phase)
            for (lS = 0; lS <= O.MaxLength; lS++)
            {
              for (jS = 0; jS < O.l[lS].N; jS++)
              {
                for (k = 0; k < psi.NWord; k++)
                  DetSqr[k] = Det[k];
                phaseS = OperatorLader(DetSqr, &(O.l[lS].o[lS * jS]), lS, O.NFermion);
                if (phaseS)
                {
                  ExpecValueRe +=
                    phase * phaseS * RealWaveFunctionPrefactorDet(psi,
                                                                  Det) * psi.ValueRe[ii][ik] * (O.l[l].V[j].Re *
                                                                                                O.l[lS].V[jS].Re -
                                                                                                O.l[l].V[j].Im *
                                                                                                O.l[lS].V[jS].Im);
                }
              }
            }
        }
      }
    }
    free(Det);
    free(DetSqr);
  }
  *VRe = ExpecValueRe;
  return 0;
}

uint32_t PsiOperatorSqrPsiRHC(const OperatorType O, const WaveFunctionType psi, double *VRe)
{
  uint32_t i, j, k, l, ii, ij, ik, lS, jS;
  double phase, phaseS;
  uint8_t *Det, *DetSqr;
  int    iOMP;
  double ExpecValueRe = 0;
  double aRe, aIm;
#pragma omp parallel private(aRe,aIm,Det,DetSqr,i,j,k,l,lS,jS,ii,ij,ik,phase,phaseS,iOMP) reduction(+ : ExpecValueRe)
  {
    Det = malloc(psi.NWord * sizeof(uint8_t));
    if (Det == NULL)
      AllocFailed();
    DetSqr = malloc(psi.NWord * sizeof(uint8_t));
    if (DetSqr == NULL)
      AllocFailed();
#pragma omp for
    for (iOMP = 1; iOMP <= (int)psi.N; iOMP++)
    {
      i = (uint32_t) iOMP;
      ii = i >> NBitsBlok;
      ik = i & BlokMask;
      ij = psi.NWord * ik;
      for (l = 0; l <= O.MaxLength; l++)
      {
        for (j = 0; j < O.l[l].N; j++)
        {
          for (k = 0; k < psi.NWord; k++)
            Det[k] = psi.Determinant[ii][ij + k];
          phase = OperatorLader(Det, &(O.l[l].o[l * j]), l, O.NFermion);
          if (phase)
            for (lS = 0; lS <= O.MaxLength; lS++)
            {
              for (jS = 0; jS < O.l[lS].N; jS++)
              {
                for (k = 0; k < psi.NWord; k++)
                  DetSqr[k] = Det[k];
                phaseS = OperatorLader(DetSqr, &(O.l[lS].o[lS * jS]), lS, O.NFermion);
                if (phaseS)
                {
                  ComplexWaveFunctionPrefactorDet(psi, DetSqr, &aRe, &aIm);
                  ExpecValueRe +=
                    phase * phaseS * (aRe * O.l[l].VRe[j] * O.l[lS].VRe[jS] * psi.ValueRe[ii][ik] +
                                      aIm * O.l[l].VRe[j] * O.l[lS].VRe[jS] * psi.ValueIm[ii][ik]);
                }
              }
            }
        }
      }
    }
    free(Det);
    free(DetSqr);
  }
  *VRe = ExpecValueRe;
  return 0;
}

uint32_t PsiOperatorSqrPsiRHR(const OperatorType O, const WaveFunctionType psi, double *VRe)
{
  uint32_t i, j, k, l, ii, ij, ik, lS, jS;
  double phase, phaseS;
  uint8_t *Det, *DetSqr;
  int    iOMP;
  double ExpecValue = 0;
#pragma omp parallel private(Det,DetSqr,i,j,k,l,lS,jS,ii,ij,ik,phase,phaseS,iOMP) reduction(+ : ExpecValue)
  {
    Det = malloc(psi.NWord * sizeof(uint8_t));
    if (Det == NULL)
      AllocFailed();
    DetSqr = malloc(psi.NWord * sizeof(uint8_t));
    if (DetSqr == NULL)
      AllocFailed();
#pragma omp for
    for (iOMP = 1; iOMP <= (int)psi.N; iOMP++)
    {
      i = (uint32_t) iOMP;
      ii = i >> NBitsBlok;
      ik = i & BlokMask;
      ij = psi.NWord * ik;
      for (l = 0; l <= O.MaxLength; l++)
      {
        for (j = 0; j < O.l[l].N; j++)
        {
          for (k = 0; k < psi.NWord; k++)
            Det[k] = psi.Determinant[ii][ij + k];
          phase = OperatorLader(Det, &(O.l[l].o[l * j]), l, O.NFermion);
          if (phase)
            for (lS = 0; lS <= O.MaxLength; lS++)
            {
              for (jS = 0; jS < O.l[lS].N; jS++)
              {
                for (k = 0; k < psi.NWord; k++)
                  DetSqr[k] = Det[k];
                phaseS = OperatorLader(DetSqr, &(O.l[lS].o[lS * jS]), lS, O.NFermion);
                if (phaseS)
                  ExpecValue +=
                    phase * phaseS * RealWaveFunctionPrefactorDet(psi,
                                                                  DetSqr) * O.l[l].VRe[j] * O.l[lS].VRe[jS] *
                    psi.ValueRe[ii][ik];
              }
            }
        }
      }
    }
    free(Det);
    free(DetSqr);
  }
  *VRe = ExpecValue;
  return 0;
}

uint32_t PsiOperatorSqrPsiMHC(const OperatorType O, const WaveFunctionType psi, double *VRe)
{
  uint32_t i, j, k, l, ii, ij, ik, lS, jS;
  double phase, phaseS;
  uint8_t *Det, *DetSqr;
  int    iOMP;
  double ExpecValueRe = 0;
  double aRe, aIm, temp;
#pragma omp parallel private(temp,aRe,aIm,Det,DetSqr,i,j,k,l,lS,jS,ii,ij,ik,phase,phaseS,iOMP) reduction(+ : ExpecValueRe)
  {
    Det = malloc(psi.NWord * sizeof(uint8_t));
    if (Det == NULL)
      AllocFailed();
    DetSqr = malloc(psi.NWord * sizeof(uint8_t));
    if (DetSqr == NULL)
      AllocFailed();
#pragma omp for
    for (iOMP = 1; iOMP <= (int)psi.N; iOMP++)
    {
      i = (uint32_t) iOMP;
      ii = i >> NBitsBlok;
      ik = i & BlokMask;
      ij = psi.NWord * ik;
      for (l = 0; l <= O.MaxLength; l++)
      {
        if (O.l[l].QComplex == 0)
        {                       /* Operator with length l real */
          for (j = 0; j < O.l[l].N; j++)
          {
            for (k = 0; k < psi.NWord; k++)
              Det[k] = psi.Determinant[ii][ij + k];
            phase = OperatorLader(Det, &(O.l[l].o[l * j]), l, O.NFermion);
            if (phase)
              for (lS = 0; lS <= O.MaxLength; lS++)
              {
                if (O.l[lS].QComplex == 0)
                {               /* Operator with length lS real */
                  for (jS = 0; jS < O.l[lS].N; jS++)
                  {
                    for (k = 0; k < psi.NWord; k++)
                      DetSqr[k] = Det[k];
                    phaseS = OperatorLader(DetSqr, &(O.l[lS].o[lS * jS]), lS, O.NFermion);
                    if (phaseS)
                    {
                      ComplexWaveFunctionPrefactorDet(psi, DetSqr, &aRe, &aIm);
                      ExpecValueRe +=
                        phase * phaseS * (aRe * O.l[l].VRe[j] * O.l[lS].VRe[jS] * psi.ValueRe[ii][ik] +
                                          aIm * O.l[l].VRe[j] * O.l[lS].VRe[jS] * psi.ValueIm[ii][ik]);
                    }
                  }
                }
                else
                {               /* Operator with length lS complex */
                  for (jS = 0; jS < O.l[lS].N; jS++)
                  {
                    for (k = 0; k < psi.NWord; k++)
                      DetSqr[k] = Det[k];
                    phaseS = OperatorLader(DetSqr, &(O.l[lS].o[lS * jS]), lS, O.NFermion);
                    if (phaseS)
                    {
                      ComplexWaveFunctionPrefactorDet(psi, DetSqr, &aRe, &aIm);
                      aIm *= -1;
                      aRe = aRe * O.l[l].VRe[j];
                      aIm = aIm * O.l[l].VRe[j];
                      temp = aRe;
                      aRe = temp * O.l[lS].V[jS].Re - aIm * O.l[lS].V[jS].Im;
                      aIm = temp * O.l[lS].V[jS].Im + aIm * O.l[lS].V[jS].Re;
                      temp = aRe;
                      aRe = temp * psi.ValueRe[ii][ik] - aIm * psi.ValueIm[ii][ik];
                      ExpecValueRe += phase * phaseS * aRe;
                    }
                  }
                }
              }
          }
        }
        else
        {                       /* Operator with length l complex */
          for (j = 0; j < O.l[l].N; j++)
          {
            for (k = 0; k < psi.NWord; k++)
              Det[k] = psi.Determinant[ii][ij + k];
            phase = OperatorLader(Det, &(O.l[l].o[l * j]), l, O.NFermion);
            if (phase)
              for (lS = 0; lS <= O.MaxLength; lS++)
              {
                if (O.l[lS].QComplex == 0)
                {               /* Operator with length lS real */
                  for (jS = 0; jS < O.l[lS].N; jS++)
                  {
                    for (k = 0; k < psi.NWord; k++)
                      DetSqr[k] = Det[k];
                    phaseS = OperatorLader(DetSqr, &(O.l[lS].o[lS * jS]), lS, O.NFermion);
                    if (phaseS)
                    {
                      ComplexWaveFunctionPrefactorDet(psi, DetSqr, &aRe, &aIm);
                      aIm *= -1;
                      temp = aRe;
                      aRe = temp * O.l[l].V[j].Re - aIm * O.l[l].V[j].Im;
                      aIm = temp * O.l[l].V[j].Im + aIm * O.l[l].V[j].Re;
                      aRe = aRe * O.l[lS].VRe[jS];
                      aIm = aIm * O.l[lS].VRe[jS];
                      temp = aRe;
                      aRe = temp * psi.ValueRe[ii][ik] - aIm * psi.ValueIm[ii][ik];
                      ExpecValueRe += phase * phaseS * aRe;
                    }
                  }
                }
                else
                {               /* Operator with length lS complex */
                  for (jS = 0; jS < O.l[lS].N; jS++)
                  {
                    for (k = 0; k < psi.NWord; k++)
                      DetSqr[k] = Det[k];
                    phaseS = OperatorLader(DetSqr, &(O.l[lS].o[lS * jS]), lS, O.NFermion);
                    if (phaseS)
                    {
                      ComplexWaveFunctionPrefactorDet(psi, DetSqr, &aRe, &aIm);
                      aIm *= -1;
                      temp = aRe;
                      aRe = temp * O.l[l].V[j].Re - aIm * O.l[l].V[j].Im;
                      aIm = temp * O.l[l].V[j].Im + aIm * O.l[l].V[j].Re;
                      temp = aRe;
                      aRe = temp * O.l[lS].V[jS].Re - aIm * O.l[lS].V[jS].Im;
                      aIm = temp * O.l[lS].V[jS].Im + aIm * O.l[lS].V[jS].Re;
                      temp = aRe;
                      aRe = temp * psi.ValueRe[ii][ik] - aIm * psi.ValueIm[ii][ik];
                      ExpecValueRe += phase * phaseS * aRe;
                    }
                  }
                }
              }
          }
        }
      }
    }
    free(Det);
    free(DetSqr);
  }
  *VRe = ExpecValueRe;
  return 0;
}

uint32_t PsiOperatorSqrPsiMHR(const OperatorType O, const WaveFunctionType psi, double *VRe)
{
  uint32_t i, j, k, l, ii, ij, ik, lS, jS;
  double phase, phaseS;
  uint8_t *Det, *DetSqr;
  int    iOMP;
  double ExpecValueRe = 0;
#pragma omp parallel private(Det,DetSqr,i,j,k,l,lS,jS,ii,ij,ik,phase,phaseS,iOMP) reduction(+ : ExpecValueRe)
  {
    Det = malloc(psi.NWord * sizeof(uint8_t));
    if (Det == NULL)
      AllocFailed();
    DetSqr = malloc(psi.NWord * sizeof(uint8_t));
    if (DetSqr == NULL)
      AllocFailed();
#pragma omp for
    for (iOMP = 1; iOMP <= (int)psi.N; iOMP++)
    {
      i = (uint32_t) iOMP;
      ii = i >> NBitsBlok;
      ik = i & BlokMask;
      ij = psi.NWord * ik;
      for (l = 0; l <= O.MaxLength; l++)
      {
        if (O.l[l].QComplex == 0)
        {
          for (j = 0; j < O.l[l].N; j++)
          {
            for (k = 0; k < psi.NWord; k++)
              Det[k] = psi.Determinant[ii][ij + k];
            phase = OperatorLader(Det, &(O.l[l].o[l * j]), l, O.NFermion);
            if (phase)
              for (lS = 0; lS <= O.MaxLength; lS++)
              {
                if (O.l[lS].QComplex == 0)
                {
                  for (jS = 0; jS < O.l[lS].N; jS++)
                  {
                    for (k = 0; k < psi.NWord; k++)
                      DetSqr[k] = Det[k];
                    phaseS = OperatorLader(DetSqr, &(O.l[lS].o[lS * jS]), lS, O.NFermion);
                    if (phaseS)
                    {
                      ExpecValueRe +=
                        phase * phaseS * RealWaveFunctionPrefactorDet(psi,
                                                                      Det) * psi.ValueRe[ii][ik] * (O.l[l].VRe[j] *
                                                                                                    O.l[lS].VRe[jS]);
                    }
                  }
                }
                else
                {
                  for (jS = 0; jS < O.l[lS].N; jS++)
                  {
                    for (k = 0; k < psi.NWord; k++)
                      DetSqr[k] = Det[k];
                    phaseS = OperatorLader(DetSqr, &(O.l[lS].o[lS * jS]), lS, O.NFermion);
                    if (phaseS)
                    {
                      ExpecValueRe +=
                        phase * phaseS * RealWaveFunctionPrefactorDet(psi,
                                                                      Det) * psi.ValueRe[ii][ik] * (O.l[l].VRe[j] *
                                                                                                    O.l[lS].V[jS].Re);
                    }
                  }
                }
              }
          }
        }
        else
        {
          for (j = 0; j < O.l[l].N; j++)
          {
            for (k = 0; k < psi.NWord; k++)
              Det[k] = psi.Determinant[ii][ij + k];
            phase = OperatorLader(Det, &(O.l[l].o[l * j]), l, O.NFermion);
            if (phase)
              for (lS = 0; lS <= O.MaxLength; lS++)
              {
                if (O.l[lS].QComplex == 0)
                {
                  for (jS = 0; jS < O.l[lS].N; jS++)
                  {
                    for (k = 0; k < psi.NWord; k++)
                      DetSqr[k] = Det[k];
                    phaseS = OperatorLader(DetSqr, &(O.l[lS].o[lS * jS]), lS, O.NFermion);
                    if (phaseS)
                    {
                      ExpecValueRe +=
                        phase * phaseS * RealWaveFunctionPrefactorDet(psi,
                                                                      Det) * psi.ValueRe[ii][ik] * (O.l[l].V[j].Re *
                                                                                                    O.l[lS].VRe[jS]);
                    }
                  }
                }
                else
                {
                  for (jS = 0; jS < O.l[lS].N; jS++)
                  {
                    for (k = 0; k < psi.NWord; k++)
                      DetSqr[k] = Det[k];
                    phaseS = OperatorLader(DetSqr, &(O.l[lS].o[lS * jS]), lS, O.NFermion);
                    if (phaseS)
                    {
                      ExpecValueRe +=
                        phase * phaseS * RealWaveFunctionPrefactorDet(psi,
                                                                      Det) * psi.ValueRe[ii][ik] * (O.l[l].V[j].Re *
                                                                                                    O.l[lS].V[jS].Re -
                                                                                                    O.l[l].V[j].Im *
                                                                                                    O.l[lS].V[jS].Im);
                    }
                  }
                }
              }
          }
        }
      }
    }
    free(Det);
    free(DetSqr);
  }
  *VRe = ExpecValueRe;
  return 0;
}

/*
*   Calculate O psi whereby the basis set used is not changed
*   input O psi and Opsi
*   output Opsi
*     Opsi.ValueRe points to newly allocated memory
*     Opsi.Determinant, Opsi.Next, Opsi.Keylist point to the same location as psi does
*     This function conservse the values O and psi point to
*/
uint32_t OperatorPsiConserveBasis(const OperatorType O, const WaveFunctionType psi, WaveFunctionType * const Opsi)
{
  if (O.QComplex == 0)
  {
    if (psi.QComplex)
    {
      return OperatorPsiConserveBasisRC(O, psi, Opsi);
    }
    return OperatorPsiConserveBasisRR(O, psi, Opsi);
  }
  if (O.QComplex == 1)
  {
    if (psi.QComplex)
    {
      return OperatorPsiConserveBasisCC(O, psi, Opsi);
    }
    return OperatorPsiConserveBasisCR(O, psi, Opsi);
  }
  if (psi.QComplex)
  {
    return OperatorPsiConserveBasisMC(O, psi, Opsi);
  }
  return OperatorPsiConserveBasisMR(O, psi, Opsi);
}

uint32_t PsiAddOperatorPsiConserveBasis(const OperatorType *O, const WaveFunctionType *psi, WaveFunctionType * const Opsi)
{
  if (O->QComplex == 0)
  {
    if (psi->QComplex)
    {
      return PsiAddOperatorPsiConserveBasisRC(O, psi, Opsi);
    }
    if (Opsi->QComplex)
    {
      return PsiAddOperatorPsiConserveBasisRRC(O, psi, Opsi);
    }
    return PsiAddOperatorPsiConserveBasisRR(O, psi, Opsi);
  }
  if (O->QComplex == 1)
  {
    if (psi->QComplex)
    {
      return PsiAddOperatorPsiConserveBasisCC(O, psi, Opsi);
    }
    return PsiAddOperatorPsiConserveBasisCR(O, psi, Opsi);
  }
  if (psi->QComplex)
  {
    return PsiAddOperatorPsiConserveBasisMC(O, psi, Opsi);
  }
  return PsiAddOperatorPsiConserveBasisMR(O, psi, Opsi);
}

/*
*   Calculate O psi whereby the basis set used is not changed and the memory of Opsi is allready allocated the values of Opsi will be set to zero on entry
*   It is assumed that Opsi.Next==psi.Next and Opsi.Determinant=psi.Determinant either in values pointed to or more efficient in pointers
*   input O psi and Opsi
*   output Opsi
*     Opsi.ValueRe points to the allocated memory on entry
*     Opsi.Determinant, Opsi.Next, Opsi.Keylist point to the same location as psi does
*     This function conservse the values O and psi point to
*/
uint32_t OperatorPsiConserveBasisPreAllocated(const OperatorType *O, const WaveFunctionType *psi,
                                              WaveFunctionType * const Opsi)
{
  if (O->QComplex == 0)
  {
    if (psi->QComplex)
    {
      return OperatorPsiConserveBasisPreAllocatedRC(O, psi, Opsi);
    }
    return OperatorPsiConserveBasisPreAllocatedRR(O, psi, Opsi);
  }
  if (O->QComplex == 1)
  {
    if (psi->QComplex)
    {
      return OperatorPsiConserveBasisPreAllocatedCC(O, psi, Opsi);
    }
    return OperatorPsiConserveBasisPreAllocatedCR(O, psi, Opsi);
  }
  if (psi->QComplex)
  {
    return OperatorPsiConserveBasisPreAllocatedMC(O, psi, Opsi);
  }
  return OperatorPsiConserveBasisPreAllocatedMR(O, psi, Opsi);
}

uint32_t PsiAddOperatorPsiConserveBasisPreAllocated(const OperatorType *O, const WaveFunctionType *psi,
                                                    WaveFunctionType * const Opsi)
{
  if (O->QComplex == 0)
  {
    if (psi->QComplex)
    {
      return PsiAddOperatorPsiConserveBasisPreAllocatedRC(O, psi, Opsi);
    }
    if (Opsi->QComplex)
    {
      return PsiAddOperatorPsiConserveBasisPreAllocatedRRC(O, psi, Opsi);
    }
    return PsiAddOperatorPsiConserveBasisPreAllocatedRR(O, psi, Opsi);
  }
  if (O->QComplex == 1)
  {
    if (psi->QComplex)
    {
      return PsiAddOperatorPsiConserveBasisPreAllocatedCC(O, psi, Opsi);
    }
    return PsiAddOperatorPsiConserveBasisPreAllocatedCR(O, psi, Opsi);
  }
  if (psi->QComplex)
  {
    return PsiAddOperatorPsiConserveBasisPreAllocatedMC(O, psi, Opsi);
  }
  return PsiAddOperatorPsiConserveBasisPreAllocatedMR(O, psi, Opsi);
}

/*
*   Calculate O psi whereby the basis set used is expanded
*   input O and psi
*   output Opsi
*     Opsi points to completely newly assigned memory
*     This function conservse the values O and psi point to
*/
uint32_t OperatorPsi(const OperatorType O, const WaveFunctionType psi, WaveFunctionType * const Opsi)
{
  if (O.QComplex == 0)
  {
    if (psi.QComplex)
    {
      return OperatorPsiRC(O, psi, Opsi);
    }
    return OperatorPsiRR(O, psi, Opsi);
  }
  if (O.QComplex == 1)
  {
    if (psi.QComplex)
    {
      return OperatorPsiCC(O, psi, Opsi);
    }
    return OperatorPsiCR(O, psi, Opsi);
  }
  if (psi.QComplex)
  {
    return OperatorPsiMC(O, psi, Opsi);
  }
  return OperatorPsiMR(O, psi, Opsi);
}

uint32_t PsiAddOperatorPsi(const OperatorType O, const WaveFunctionType psi, WaveFunctionType * const Opsi)
{
  if (O.QComplex == 0)
  {
    if (psi.QComplex)
    {
      return PsiAddOperatorPsiRC(O, psi, Opsi);
    }
    if (Opsi->QComplex)
    {
      return PsiAddOperatorPsiRRC(O, psi, Opsi);
    }
    return PsiAddOperatorPsiRR(O, psi, Opsi);
  }
  if (O.QComplex == 1)
  {
    if (psi.QComplex)
    {
      return PsiAddOperatorPsiCC(O, psi, Opsi);
    }
    return PsiAddOperatorPsiCR(O, psi, Opsi);
  }
  if (psi.QComplex)
  {
    return PsiAddOperatorPsiMC(O, psi, Opsi);
  }
  return PsiAddOperatorPsiMR(O, psi, Opsi);
}

/*  Serial version of previous routine
*
*   Calculate O psi whereby the basis set used is expanded
*   input O and psi
*   output Opsi
*     Opsi points to completely newly assigned memory
*     This function conservse the values O and psi point to
*/
uint32_t OperatorPsiSerial(const OperatorType O, const WaveFunctionType psi, WaveFunctionType * const Opsi)
{
  if (O.QComplex == 0)
  {
    if (psi.QComplex)
    {
      return OperatorPsiSerialRC(O, psi, Opsi);
    }
    return OperatorPsiSerialRR(O, psi, Opsi);
  }
  if (O.QComplex == 1)
  {
    if (psi.QComplex)
    {
      return OperatorPsiSerialCC(O, psi, Opsi);
    }
    return OperatorPsiSerialCR(O, psi, Opsi);
  }
  if (psi.QComplex)
  {
    return OperatorPsiSerialMC(O, psi, Opsi);
  }
  return OperatorPsiSerialMR(O, psi, Opsi);
}

uint32_t PsiAddOperatorPsiSerial(const OperatorType O, const WaveFunctionType psi, WaveFunctionType * const Opsi)
{
  if (O.QComplex == 0)
  {
    if (psi.QComplex)
    {
      return PsiAddOperatorPsiSerialRC(O, psi, Opsi);
    }
    if (Opsi->QComplex)
    {
      return PsiAddOperatorPsiSerialRRC(O, psi, Opsi);
    }
    return PsiAddOperatorPsiSerialRR(O, psi, Opsi);
  }
  if (O.QComplex == 1)
  {
    if (psi.QComplex)
    {
      return PsiAddOperatorPsiSerialCC(O, psi, Opsi);
    }
    return PsiAddOperatorPsiSerialCR(O, psi, Opsi);
  }
  if (psi.QComplex)
  {
    return PsiAddOperatorPsiSerialMC(O, psi, Opsi);
  }
  return PsiAddOperatorPsiSerialMR(O, psi, Opsi);
}

/*  Calcualtes <psi|O|psi>
*
*/
uint32_t PsiOperatorPsi(const OperatorType O, const WaveFunctionType psi, double *VRe, double *VIm)
{
  if (O.QComplex == 0)
  {
    if (psi.QComplex)
    {
      return PsiOperatorPsiRC(O, psi, VRe, VIm);
    }
    return PsiOperatorPsiRR(O, psi, VRe);
  }
  if (O.QComplex == 1)
  {
    if (psi.QComplex)
    {
      return PsiOperatorPsiCC(O, psi, VRe, VIm);
    }
    return PsiOperatorPsiCR(O, psi, VRe, VIm);
  }
  if (psi.QComplex)
  {
    return PsiOperatorPsiMC(O, psi, VRe, VIm);
  }
  return PsiOperatorPsiMR(O, psi, VRe, VIm);
}

/*  Calcualtes <psi|O|psi>
 *  Under the assumption that O^{\dag} = O
 *  O is Hermitian
 */
uint32_t PsiHermitianOperatorPsi(const OperatorType O, const WaveFunctionType psi, double *VRe)
{
  if (O.QComplex == 0)
  {
    if (psi.QComplex)
    {
      return PsiOperatorPsiRHC(O, psi, VRe);
    }
    return PsiOperatorPsiRHR(O, psi, VRe);
  }
  if (O.QComplex == 1)
  {
    if (psi.QComplex)
    {
      return PsiOperatorPsiCHC(O, psi, VRe);
    }
    return PsiOperatorPsiCHR(O, psi, VRe);
  }
  if (psi.QComplex)
  {
    return PsiOperatorPsiMHC(O, psi, VRe);
  }
  return PsiOperatorPsiMHR(O, psi, VRe);
}


/*  Calcualtes <psi|OO|psi>
*
*/
uint32_t PsiOperatorSqrPsi(const OperatorType O, const WaveFunctionType psi, double *VRe, double *VIm)
{
  if (O.QComplex == 0)
  {
    if (psi.QComplex)
    {
      return PsiOperatorSqrPsiRC(O, psi, VRe, VIm);
    }
    return PsiOperatorSqrPsiRR(O, psi, VRe);
  }
  if (O.QComplex == 1)
  {
    if (psi.QComplex)
    {
      return PsiOperatorSqrPsiCC(O, psi, VRe, VIm);
    }
    return PsiOperatorSqrPsiCR(O, psi, VRe, VIm);
  }
  if (psi.QComplex)
  {
    return PsiOperatorSqrPsiMC(O, psi, VRe, VIm);
  }
  return PsiOperatorSqrPsiMR(O, psi, VRe, VIm);
}

/*  Calcualtes <psi|OO|psi>
 *  Under the assumption that O^{\dag} = O
 *  O is Hermitian
 */
uint32_t PsiHermitianOperatorSqrPsi(const OperatorType O, const WaveFunctionType psi, double *VRe)
{
  if (O.QComplex == 0)
  {
    if (psi.QComplex)
    {
      return PsiOperatorSqrPsiRHC(O, psi, VRe);
    }
    return PsiOperatorSqrPsiRHR(O, psi, VRe);
  }
  if (O.QComplex == 1)
  {
    if (psi.QComplex)
    {
      return PsiOperatorSqrPsiCHC(O, psi, VRe);
    }
    return PsiOperatorSqrPsiCHR(O, psi, VRe);
  }
  if (psi.QComplex)
  {
    return PsiOperatorSqrPsiMHC(O, psi, VRe);
  }
  return PsiOperatorSqrPsiMHR(O, psi, VRe);
}


/*  Calcualtes the matrix of dimensions [i,j] <psiPrime_i|O|psi_j>
 *  On entry V.Ni contains the length of psiPrime
 *  On entry V.Nj contains the length of psi
 */
uint32_t PsiPrimeOperatorPsiMatrix(const OperatorType O, const WaveFunctionType * psiPrime,
                                   const WaveFunctionType * psi, MatrixType * const V)
{
  uint32_t i, j;
  WaveFunctionType Opsi;
  double DotRe, DotIm;
  for (j = 0; j < V->Nj; j++)
  {
    if (OperatorPsi(O, psi[j], &Opsi))
    {
      printf("OperatorPsi failed in PsiPrimeOperatorPsiMatrix at level j=%5i", j);
      return 1;
    }
    for (i = 0; i < V->Ni; i++)
    {
      WaveFunctionDotProduct(psiPrime[i], Opsi, &DotRe, &DotIm);
      V->VRe[i][j] = DotRe;
      if (V->QComplex)
        V->VIm[i][j] = DotIm;
    }
    WaveFunctionFree(&Opsi);
  }
  return 0;
}

/*  Calcualtes the matrix of dimensions [i,j] <psi_i|O|psi_j>
 *  On entry V.Ni contains the length of psi
 *  On entry V.Nj contains the length of psi
 */
uint32_t PsiOperatorPsiMatrix(const OperatorType O, const WaveFunctionType * psi, MatrixType * const V)
{
  return PsiPrimeOperatorPsiMatrix(O, psi, psi, V);
}

/*  Calcualtes <psiPrime|O|psi>
 *
 */
uint32_t PsiPrimeOperatorPsi(const OperatorType O, const WaveFunctionType psiPrime, const WaveFunctionType psi,
                             double *VRe, double *VIm)
{
  MatrixType V;
  uint32_t returnVal;
  V.Ni = 1;
  V.Nj = 1;
  V.QComplex = 1;
  InitMatrix(&V);
  V.VRe[0][0] = 0;
  V.VIm[0][0] = 0;

  returnVal = PsiPrimeOperatorPsiMatrix(O, &psiPrime, &psi, &V);

  *VRe = V.VRe[0][0];
  *VIm = V.VIm[0][0];
  return returnVal;
}
