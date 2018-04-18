
uint32_t ComplexRealWaveFunctionAddComplexMultipliedUnionBasis(WaveFunctionType *const a,const WaveFunctionType b,const double alphaRe,const double alphaIm);
uint32_t ComplexRealWaveFunctionAddMultipliedUnionBasis(WaveFunctionType *const a,const WaveFunctionType b,const double alphaRe);
uint32_t ComplexWaveFunctionAddElement(WaveFunctionType *const psi,const uint8_t *const det,const double aRe,const double aIm);
uint32_t ComplexRealWaveFunctionAddUnionBasis(WaveFunctionType *const a,const WaveFunctionType b);
void ComplexWaveFunctionAddElementConserveBasis(WaveFunctionType psi,const uint8_t *const det,const double aRe,const double aIm);
void ComplexRealWaveFunctionAddComplexMultipliedConserveBasis(WaveFunctionType *const a,const WaveFunctionType b,const double alphaRe,const double alphaIm);
void ComplexRealWaveFunctionAddComplexMultipliedSameBasis(WaveFunctionType a,const WaveFunctionType b,const double alphaRe,const double alphaIm);
void ComplexRealWaveFunctionDotProductSameBasis(const WaveFunctionType a,const WaveFunctionType b,double *DotRe,double *DotIm);
uint32_t hash(const void *key,size_t length,const uint32_t NBitsKey);
void ComplexRealWaveFunctionDotProduct(const WaveFunctionType *a,const WaveFunctionType *b,double *DotRe,double *DotIm);
