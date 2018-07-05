void ComplexRealWaveFunctionDotProduct(const WaveFunctionType *a, const WaveFunctionType *b, double *DotRe, double *DotIm);
void ComplexRealWaveFunctionDotProductSameBasis(const WaveFunctionType a, const WaveFunctionType b, double *DotRe,double *DotIm);
void ComplexRealWaveFunctionAddComplexMultipliedSameBasis(WaveFunctionType a, const WaveFunctionType b,const double alphaRe, const double alphaIm);
void ComplexRealWaveFunctionAddComplexMultipliedConserveBasis(WaveFunctionType * const a, const WaveFunctionType b,const double alphaRe, const double alphaIm);
uint32_t ComplexRealWaveFunctionAddUnionBasis(WaveFunctionType * const a, const WaveFunctionType b);
uint32_t ComplexRealWaveFunctionAddMultipliedUnionBasis(WaveFunctionType * const a, const WaveFunctionType b, const double alphaRe);
uint32_t ComplexRealWaveFunctionAddComplexMultipliedUnionBasis(WaveFunctionType * const a, const WaveFunctionType b, const double alphaRe, const double alphaIm);
