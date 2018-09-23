#ifndef ATOMICOPERATORS_SLATERINTEGRALS_H_
#define ATOMICOPERATORS_SLATERINTEGRALS_H_

SSlaterIntegral *GetSlaterIntegrals(const vector<vector<double> > &basisfuns,const vector<double> &x,const vector<int> &indToK, unsigned int &NIntegrals);

SSlaterIntegral *GetSlaterIntegrals(const vector<Interpolation::InterpolatingFunction *> &ipolfuns,const vector<int> &indToK, unsigned int &NIntegrals);

void CalculateSlaterIntegrals(const int NIntegrals, SSlaterIntegral *integrals, const int dim,
                              double **kernels, vector<vector<double> > basisfuns, vector<vector<double> > dualfuns);

void CalculateSlaterIntegrals(const int NIntegrals, SSlaterIntegral *integrals, uint32_t dim,
                              vector<vector<double> > &kernels, vector<vector<double> > &basisfuns, vector<vector<double> > &dualfuns);

void vecMul(const int n, const double *f, const double *g, double *fout);

double SlaterIntegral(const int n, const double *f1, const double *f2, const double *f3,
                      const double *f4, const double *Mr1r2);

double ** CreateAllKernels(const unsigned int kMax, const unsigned int NGridPoints, const double *x);

void CreateKernel(const unsigned int NGridPoints, const double *x, double *kernel, const unsigned int k);

double *CreateWeights(const unsigned int NGridPoints, const double *x);

vector<double> CreateWeights(const unsigned int NGridPoints, const vector<double> &x);

vector<double> Dual(const vector<double> &f, const vector<double> &w);

void CreateDualFunctions(const vector<vector<double> > &basisfuns,const vector<double> &x, vector<vector<double> > &dualfuns);

SSlaterIntegral *InitializeIntegrals(const vector<int> &indToK, unsigned int &NIntegrals, unsigned int *kMax);

unsigned int *GetK(const int *kappa, int *n);

void CreateNonRelKernel(const vector<double> &x, const uint32_t k,vector<double> &kernel);

void CreateAllKernels(const unsigned int kMax,const vector<double> &x, vector<vector<double> > &kernels);

int Sign(const int x);

#endif
