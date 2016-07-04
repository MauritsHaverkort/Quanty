
/*
*   calculate <a|b>
*/
void ComplexRealWaveFunctionDotProduct(const WaveFunctionType *a, const WaveFunctionType *b, double *DotRe, double *DotIm)
{
  uint32_t i, j, NextIndex, PsiAIsPsiB;
  int    iOMP;
  double DotReOMP = 0, DotImOMP = 0;
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
        DotReOMP += a->ValueRe[BlokElement(i)] * b->ValueRe[BlokElement(NextIndex)];
        DotImOMP -= a->ValueIm[BlokElement(i)] * b->ValueRe[BlokElement(NextIndex)];
        NextIndex = 0;
      }
      else
      {
        NextIndex = b->Next[BlokElement(NextIndex)];
      }
    }
  }
  *DotRe = DotReOMP;
  *DotIm = DotImOMP;
}

/*
 *   calculate <a|b>
 *   The first entry is complex conjugated
 */
void ComplexRealWaveFunctionDotProductSameBasis(const WaveFunctionType a, const WaveFunctionType b, double *DotRe,
                                                double *DotIm)
{
  uint32_t iBlok, iElement;
  int    iOMP;
  double DotReOMPA = 0, DotImOMPA = 0, DotReOMPB = 0, DotImOMPB = 0;
  iBlok = a.NBlok - 1;
#pragma omp parallel for private(iElement) reduction(+ : DotReOMPA,DotImOMPA)
  for (iOMP = 0; iOMP <= (int)(a.N - iBlok * BlokSize); ++iOMP)
  {
    iElement = (uint32_t) iOMP;
    DotReOMPA += a.ValueRe[iBlok][iElement] * b.ValueRe[iBlok][iElement];
    DotImOMPA -= a.ValueIm[iBlok][iElement] * b.ValueRe[iBlok][iElement];
  }
#pragma omp parallel for private(iElement,iBlok) reduction(+ : DotReOMPB, DotImOMPB)
  for (iOMP = 0; iOMP < (int)(a.NBlok - 1); iOMP++)
  {
    iBlok = (uint32_t) iOMP;
    for (iElement = 0; iElement < BlokSize; ++iElement)
    {
      DotReOMPB += a.ValueRe[iBlok][iElement] * b.ValueRe[iBlok][iElement];
      DotImOMPB -= a.ValueIm[iBlok][iElement] * b.ValueRe[iBlok][iElement];
    }
  }
  *DotRe = DotReOMPA + DotReOMPB;
  *DotIm = DotImOMPA + DotImOMPB;
}


/*
 *   a->a+alpha b
 *   Look out this function does not alter the values stored in a nor b
 *   As a and b only store pointers
 *   This function does alter the values a points to
 */
void ComplexRealWaveFunctionAddComplexMultipliedSameBasis(WaveFunctionType a, const WaveFunctionType b,
                                                          const double alphaRe, const double alphaIm)
{
  uint32_t iBlok, iElement;
  int    iOMP;
  iBlok = a.NBlok - 1;
#pragma omp parallel for private(iElement)
  for (iOMP = 0; iOMP <= (int)(a.N - iBlok * BlokSize); ++iOMP)
  {
    iElement = (uint32_t) iOMP;
    a.ValueRe[iBlok][iElement] += alphaRe * b.ValueRe[iBlok][iElement] - alphaIm * b.ValueIm[iBlok][iElement];
    a.ValueIm[iBlok][iElement] += alphaIm * b.ValueRe[iBlok][iElement];
  }
#pragma omp parallel for private(iElement,iBlok)
  for (iOMP = 0; iOMP < (int)(a.NBlok - 1); iOMP++)
  {
    iBlok = (uint32_t) iOMP;
    for (iElement = 0; iElement < BlokSize; ++iElement)
    {
      a.ValueRe[iBlok][iElement] += alphaRe * b.ValueRe[iBlok][iElement] - alphaIm * b.ValueIm[iBlok][iElement];
      a.ValueIm[iBlok][iElement] += alphaIm * b.ValueRe[iBlok][iElement];
    }
  }
}


/*
 *   (a) <- a+alpha b
 *   This function leaves b in the memory and adds values to a
 *   if b_i is not in the basis of a it will be ignored
 *   Therefor values of a might change and a is returned from the function
 *   If a is an empty function b is copied into a (use of copy is strongly recommanded)
 */
void ComplexRealWaveFunctionAddComplexMultipliedConserveBasis(WaveFunctionType * const a, const WaveFunctionType b,
                                                              const double alphaRe, const double alphaIm)
{
  uint32_t i;
  int    iOMP;
#pragma omp parallel for private(i)
  for (iOMP = 1; iOMP <= (int)b.N; iOMP++)
  {
    i = (uint32_t) iOMP;
    ComplexWaveFunctionAddElementConserveBasis(*a, &b.Determinant[BlokElementDet(i, b.NWord, 0)],
                                               alphaRe * b.ValueRe[BlokElement(i)],
                                               alphaIm * b.ValueRe[BlokElement(i)]);
  }
}

/*
 *   (a+new allocation) <- a+alpha b
 *   This function leaves b in the memory and adds values to a
 *   When aditional space for a needs to be allocated it will
 *   Therefor values of a might change and a is returned from the function
 *   If a is an empty function b is copied into a (use of copy is strongly recommanded)
 */
uint32_t ComplexRealWaveFunctionAddUnionBasis(WaveFunctionType * const a, const WaveFunctionType b)
{
  uint32_t i, j;
  for (i = b.N; i != 0; --i)
  {
    if (ComplexWaveFunctionAddElement(a, &b.Determinant[BlokElementDet(i, b.NWord, 0)], b.ValueRe[BlokElement(i)], 0) != 0)
    {                           /* Allocation of new memory for a failed Remove all elements added so far and return 1 */
      for (j = b.N; j != i; --j)
      {
        ComplexWaveFunctionAddElementConserveBasis(*a, &b.Determinant[BlokElementDet(i, b.NWord, 0)],-b.ValueRe[BlokElement(i)],0);
      }
      return 1;
    }
  }
  return 0;
}
uint32_t ComplexRealWaveFunctionAddMultipliedUnionBasis(WaveFunctionType * const a, const WaveFunctionType b, const double alphaRe)
{
  uint32_t i, j;
  for (i = b.N; i != 0; --i)
  {
    if (ComplexWaveFunctionAddElement(a, &b.Determinant[BlokElementDet(i, b.NWord, 0)], alphaRe * b.ValueRe[BlokElement(i)],0) != 0)
    {                           /* Allocation of new memory for a failed Remove all elements added so far and return 1 */
      for (j = b.N; j != i; --j)
      {
        ComplexWaveFunctionAddElementConserveBasis(*a, &b.Determinant[BlokElementDet(i, b.NWord, 0)],-alphaRe * b.ValueRe[BlokElement(i)],0);
      }
      return 1;
    }
  }
  return 0;
}
uint32_t ComplexRealWaveFunctionAddComplexMultipliedUnionBasis(WaveFunctionType * const a, const WaveFunctionType b, const double alphaRe, const double alphaIm)
{
  uint32_t i, j;
  for (i = b.N; i != 0; --i)
  {
    if (ComplexWaveFunctionAddElement
        (a, &b.Determinant[BlokElementDet(i, b.NWord, 0)], alphaRe * b.ValueRe[BlokElement(i)],
         alphaIm * b.ValueRe[BlokElement(i)]) != 0)
    {                           /* Allocation of new memory for a failed Remove all elements added so far and return 1 */
      for (j = b.N; j != i; --j)
      {
        ComplexWaveFunctionAddElementConserveBasis(*a, &b.Determinant[BlokElementDet(i, b.NWord, 0)],
                                                   -alphaRe * b.ValueRe[BlokElement(i)],
                                                   -alphaIm * b.ValueRe[BlokElement(i)]);
      }
      return 1;
    }
  }
  return 0;
}
