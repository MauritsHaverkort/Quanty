
/*
*   calculate <a|b>
*/
void RealComplexWaveFunctionDotProduct(const WaveFunctionType a, const WaveFunctionType b, double *DotRe, double *DotIm)
{
  uint32_t i,j,NextIndex,PsiAIsPsiB;
  int iOMP;
  double DotReOMP=0,DotImOMP=0;
  #pragma omp parallel for private(i,j,NextIndex,PsiAIsPsiB) reduction(+ : DotReOMP,DotImOMP)
  for(iOMP=1;iOMP<=(int)(a.N);++iOMP)
  {
    i=(uint32_t)iOMP;
    NextIndex=b.KeyList[hash(&a.Determinant[BlokElementDet(i,a.NWord,0)],a.NWord,a.NBitsKey)];
    while(NextIndex)
    {
      PsiAIsPsiB=1;
      for(j=0;j<b.NWord;j++)
      {
        if(b.Determinant[BlokElementDet(NextIndex,b.NWord,j)]!=a.Determinant[BlokElementDet(i,a.NWord,j)]) PsiAIsPsiB=0;
      }
      if(PsiAIsPsiB)
      {/* a=b found add the product to Dot and leave loop  */
        DotReOMP+=a.ValueRe[BlokElement(i)]*b.ValueRe[BlokElement(NextIndex)];
        DotImOMP+=a.ValueRe[BlokElement(i)]*b.ValueIm[BlokElement(NextIndex)];
        NextIndex=0;
      }
      else
      {
        NextIndex=b.Next[BlokElement(NextIndex)];
      }
    }
  }
  *DotRe=DotReOMP;
  *DotIm=DotImOMP;
}


