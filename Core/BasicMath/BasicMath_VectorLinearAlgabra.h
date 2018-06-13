void RealVectorCopy(double *const u, double *const v, uint32_t N);
void ComplexVectorCopy(double *const uRe, double *const uIm, double *const vRe, double *const vIm, uint32_t N);
void RealVectorMultiply(double *const u, const uint32_t N, const double alphaRe);
void ComplexVectorMultiply(double *const uRe, double *const uIm, const uint32_t N, const double alphaRe,const double alphaIm);
void RealVectorDotVectorMahan(const double *const u, const double *const v, uint32_t N, double *const a);
void RealVectorDotVector(const double *const u, const double *const v, uint32_t N, double *const a);
void ComplexVectorDotVector(const double *const uRe, const double *const uIm, const double *const vRe,const double *const vIm, const uint32_t N, double *const aRe, double *const aIm);
void RealVectorAdd(double *const u, double *const v, uint32_t N, double alphaRe);
void ComplexVectorAdd(double *const uRe, double *const uIm, double *const vRe, double *const vIm, const uint32_t N,const double alphaRe, const double alphaIm);
