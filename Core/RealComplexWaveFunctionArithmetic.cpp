/*
*   calculate <a|b>
*/
void RealComplexWaveFunctionDotProduct(const WaveFunctionType *a, const WaveFunctionType *b, double *DotRe, double *DotIm)
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
        DotImOMP += a->ValueRe[BlokElement(i)] * b->ValueIm[BlokElement(NextIndex)];
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
void RealComplexWaveFunctionDotProductSameBasis(const WaveFunctionType a, const WaveFunctionType b, double *DotRe,
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
    DotImOMPA += a.ValueRe[iBlok][iElement] * b.ValueIm[iBlok][iElement];
  }
#pragma omp parallel for private(iElement,iBlok) reduction(+ : DotReOMPB, DotImOMPB)
  for (iOMP = 0; iOMP < (int)(a.NBlok - 1); iOMP++)
  {
    iBlok = (uint32_t) iOMP;
    for (iElement = 0; iElement < BlokSize; ++iElement)
    {
      DotReOMPB += a.ValueRe[iBlok][iElement] * b.ValueRe[iBlok][iElement];
      DotImOMPB += a.ValueRe[iBlok][iElement] * b.ValueIm[iBlok][iElement];
    }
  }
  *DotRe = DotReOMPA + DotReOMPB;
  *DotIm = DotImOMPA + DotImOMPB;
}
