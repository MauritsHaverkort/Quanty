/* This file was automatically generated.  Do not edit! */
int Sign(const int x);
unsigned int *GetK(const int *kappa,int *n);
void Dual(const int n,const double *f,const double *weights,double *fdual);
double *CreateWeights(const unsigned int NGridPoints,const double *x);
void CreateKernel(const unsigned int NGridPoints,const double *x,double *kernel,const unsigned int k);
void vecMul(const int n,const double *f,const double *g,double *fout);
double SlaterIntegral(const int n,const double *f1,const double *f2,const double *f3,const double *f4,const double *Mr1r2);
void CalculateSlaterIntegrals(const int NIntegrals,SSlaterIntegral *integrals,const int dim,double **kernels,double **basisfuns,double **dualfuns);
double **CreateAllKernels(const unsigned int kMax,const unsigned int NGridPoints,const double *x);
SSlaterIntegral *InitializeIntegrals(const int **tupel,unsigned int NIntegrals,unsigned int *kMax);
double **CreateDualFunctions(const unsigned int NBasisFuns,const unsigned int NGridPoints,double **basisfuns,const double *x);
SSlaterIntegral *GetSlaterIntegrals(double **basisfuns,unsigned int NBasisFuns,const double *x,unsigned int NGridPoints,const int **tupel,unsigned int NIntegrals);
//double dsecnd();
//double cblas_ddot();
//void cblas_dsymv();
