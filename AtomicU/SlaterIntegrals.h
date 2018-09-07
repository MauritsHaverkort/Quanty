#ifndef ATOMICU_SLATERINTEGRALS_H_
#define ATOMICU_SLATERINTEGRALS_H_

SSlaterIntegral *GetSlaterIntegrals(double **basisfuns, unsigned int NBasisFuns, const double *x,
                                    unsigned int NGridPoints,const int **tupel, unsigned int NIntegrals);

SSlaterIntegral *GetSlaterIntegrals(double **basisfuns, unsigned int NBasisFuns, const double *x,
                                    unsigned int NGridPoints,vector<int> indToK, unsigned int &NIntegrals);

SSlaterIntegral *GetSlaterIntegrals(vector<vector<double> > basisfuns, vector<double> x,vector<int> indToK, unsigned int &NIntegrals);

void CalculateSlaterIntegrals(const int NIntegrals, SSlaterIntegral *integrals, const int dim,
                              double **kernels, double **basisfuns, double **dualfuns);

void CalculateSlaterIntegrals(const int NIntegrals, SSlaterIntegral *integrals, const int dim,
                              double **kernels, vector<vector<double> > basisfuns, vector<vector<double> > dualfuns);

void vecMul(const int n, const double *f, const double *g, double *fout);

double SlaterIntegral(const int n, const double *f1, const double *f2, const double *f3,
                      const double *f4, const double *Mr1r2);

double ** CreateAllKernels(const unsigned int kMax, const unsigned int NGridPoints, const double *x);

void CreateKernel(const unsigned int NGridPoints, const double *x, double *kernel, const unsigned int k);

double *CreateWeights(const unsigned int NGridPoints, const double *x);

vector<double> CreateWeights(const unsigned int NGridPoints, vector<double> x);

double **CreateDualFunctions(const unsigned int NBasisFuns, const unsigned int NGridPoints,
                             double **basisfuns, const double *x);

vector<double> Dual(vector<double> f, vector<double> w);

vector<vector<double> > CreateDualFunctions(const unsigned int NBasisFuns, const unsigned int NGridPoints,
		vector<vector<double> > basisfuns, vector<double> x);

void Dual(const int n, const double *f, const double *weights, double *fdual);

SSlaterIntegral *InitializeIntegrals(const int **tupel, unsigned int NIntegrals,
                                     unsigned int *kMax);

SSlaterIntegral *InitializeIntegrals(vector<int> indToK, unsigned int &NIntegrals, unsigned int *kMax);

unsigned int *GetK(const int *kappa, int *n);

int Sign(const int x);

#endif
