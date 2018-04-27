typedef struct
{
    unsigned int n[4];
    int kappa[4];
    
}SQuantumNumbers;

typedef struct
{
    double *val;
    int Nval;
    unsigned int *k;
    SQuantumNumbers tau;
    
}SSlaterIntegral;

#include "SlaterIntegrals.h"

SSlaterIntegral *GetSlaterIntegrals(double **basisfuns, unsigned int NBasisFuns, const double *x,
                                    unsigned int NGridPoints,const int **tupel, unsigned int NIntegrals)
{
        
    SSlaterIntegral *si;
    unsigned int kMax, n;
    double **kernels, **dualfuns;
    
    n = 2* NGridPoints;
    
    
    
    dualfuns = CreateDualFunctions(NBasisFuns, NGridPoints, basisfuns,x);
    
    si = InitializeIntegrals(tupel, NIntegrals,&kMax);

    
    kernels = CreateAllKernels(kMax, NGridPoints, x); 
    
   
    CalculateSlaterIntegrals(NIntegrals, si,n,kernels,basisfuns,dualfuns);
/*
  	int i,j;
    FILE *f = fopen("slaterout4","w");
    for ( i = 0; i < NIntegrals; i++)
    {
    	for ( j = 0; j < si[i].Nval; j++)
    	{
    		fprintf(f,"%.16e\n",si[i].val[j]);
    	}
    }
    fclose(f);
*/
    return si;
}



void CalculateSlaterIntegrals(const int NIntegrals, SSlaterIntegral *integrals, const int dim,
                              double **kernels, double **basisfuns, double **dualfuns)
{
     
    int i,j, k;
    unsigned int *n;
    double Fk;
    
    for( i = 0; i < NIntegrals; i++)
    {
        
        n = integrals[i].tau.n;
        for( j = 0; j < integrals[i].Nval; j++)
        {
            
            k = integrals[i].k[j];
            Fk = SlaterIntegral(dim, dualfuns[n[0]], dualfuns[n[1]],
                                basisfuns[n[2]],basisfuns[n[3]],kernels[k]);
            integrals[i].val[j] = Fk;
            
        }
    }
};

void vecMul(const int n, const double *f, const double *g, double *fout)
{
    int i;

    for ( i = 0; i < n; i++)
    {
        fout[i] = f[i] * g[i];
    }
}
  

double SlaterIntegral(const int n, const double *f1, const double *f2, const double *f3, 
                      const double *f4, const double *Mr1r2)
{
  const char UPLO='L';
  const double alpha=1.0;
  const int INCXY=1;
  const double beta=0;
  double *f13, *f24, *Mf, result;
    
  f13 = (double *) malloc(n*sizeof(double));
  f24 = (double *) malloc(n*sizeof(double));
  Mf = (double *) malloc(n*sizeof(double));
    
  vecMul(n, f1, f3, f13);
  vecMul(n, f2, f4, f24);
    
  dsymv(&UPLO, &n, &alpha, Mr1r2, &n, f24, &INCXY, &beta, Mf, &INCXY);
  result = ddot(&n, f13, &INCXY, Mf, &INCXY);
  
  free(f13);
  free(f24);
  free(Mf);
  return result;
} 

double ** CreateAllKernels(const unsigned int kMax, const unsigned int NGridPoints, const double *x)
{
    unsigned int k, n;
    double ** kernels;
    
    kernels = (double **) malloc((kMax+1)*sizeof(double *));
    n = 2 * NGridPoints;
    
    for ( k = 0; k <= kMax; k++ )
    {
        kernels[k] = (double *) malloc(n*n*sizeof(double));
        CreateKernel(NGridPoints, x, kernels[k], k);
    }
    return kernels;
}

void CreateKernel(const unsigned int NGridPoints, const double *x, double *kernel, const unsigned int k)
{
    int i, j, n;
    double val;
    
    n = NGridPoints;
    
    
    for( i = 0; i < NGridPoints; i++)
    {
        for( j = i; j < NGridPoints; j++)
        {
            if (i == 0)
            {
                val = 0;
            }
            else
            {
                val = pow(x[i],k)/pow(x[j],k+1);
            }
        
            kernel[j+2*n*i] = val;
            kernel[i+2*n*j] = val;
            
            kernel[(j+n)+2*n*i] = val;
            kernel[i+2*n*(j+n)] = val;
            
            kernel[j+2*n*(i+n)] = val;
            kernel[(i+n)+2*n*j] = val;

            kernel[(j+n)+2*n*(i+n)] = val;
            kernel[(i+n)+2*n*(j+n)] = val;            
        }
    }
    
    
    
    
}

double *CreateWeights(const unsigned int NGridPoints, const double *x)
{
    int i;
    double x0 = 0;
    double *weights = (double *) malloc(NGridPoints*2*sizeof(double));
    
    for( i = 0; i < NGridPoints-1; i++)
    {
        weights[i] = 2 * (x[i] - x0);
        weights[i+NGridPoints] = 2 * (x[i] - x0);
        x0 = x[i] + weights[i] / 2;
    }
    weights[NGridPoints-2] += x[i]-x0;
    weights[2*NGridPoints-2] += x[i]-x0;
    return weights;
}

double **CreateDualFunctions(const unsigned int NBasisFuns, const unsigned int NGridPoints,
                             double **basisfuns, const double *x)
{
    int i;
    double **dualfuns, *weights;
    
    weights = CreateWeights(NGridPoints,x);
    dualfuns = (double **) malloc(NBasisFuns*sizeof(double *));
    for ( i = 0; i < NBasisFuns; i++ )
    {
        dualfuns[i] = (double *) malloc(2*NGridPoints*sizeof(double));
        Dual(2*NGridPoints, basisfuns[i], weights, dualfuns[i]);
    }
    free(weights);
    return dualfuns;
}

void Dual(const int n, const double *f, const double *weights, double *fdual)
{
    vecMul(n, f, weights, fdual);
}

SSlaterIntegral *InitializeIntegrals(const int **tupel, unsigned int NIntegrals,
                                     unsigned int *kMax)
{
    unsigned int i, k;
    SSlaterIntegral *integrals;
    
    integrals = (SSlaterIntegral *) malloc((NIntegrals)*sizeof(SSlaterIntegral));
    for( i = 0; i < NIntegrals; i++)
    {
                
    	integrals[i].tau.n[0] = tupel[i][0];
    	integrals[i].tau.n[1] = tupel[i][2];
    	integrals[i].tau.n[2] = tupel[i][4];
    	integrals[i].tau.n[3] = tupel[i][6];

    	integrals[i].tau.kappa[0] = tupel[i][1];
    	integrals[i].tau.kappa[1] = tupel[i][3];
    	integrals[i].tau.kappa[2] = tupel[i][5];
    	integrals[i].tau.kappa[3] = tupel[i][7];

    }
    *kMax = 0;
    for( i = 0; i < NIntegrals; i++ )
    {
        integrals[i].k = GetK(integrals[i].tau.kappa, &integrals[i].Nval);
        k = integrals[i].k[integrals[i].Nval-1];
        if (k > *kMax)
        {
            *kMax = k;
        }
        integrals[i].val = (double *) malloc(integrals[i].Nval*sizeof(double));
    }
    
    return integrals;
}

unsigned int *GetK(const int *kappa, int *n)
{
    unsigned int k0=0, k1=0, i, *k;
    
    
    k0 = fmax( abs( abs(kappa[0])-abs(kappa[2]) ), abs( abs(kappa[1])-abs(kappa[3]) ) );
    k1 = fmin( abs(kappa[0])+abs(kappa[2]), abs(kappa[1])+abs(kappa[3]) ) -1;
    
    if ( (abs(kappa[0])+abs(kappa[2]) + (Sign(kappa[0])+Sign(kappa[2]))/2 + k0) % 2 == 0)
    {
        k0++;
    }
    *n = (k1-k0)/2 + 1;
    if ( k0 > k1 )
    {
      k=NULL;
      printf("Warning: k0(%d) > k1(%d)\n This has to be fixed...",k0,k1);
    }
    else
    {
        k = (unsigned int *) malloc((*n) * sizeof(unsigned int));
        for( i = 0; i < *n; i++)
        {
    
            k[i] = k0 + 2*i;
    
        }
    }
    return k;
}

int Sign(const int x)
{
    return (x > 0) - (x < 0);
}





















