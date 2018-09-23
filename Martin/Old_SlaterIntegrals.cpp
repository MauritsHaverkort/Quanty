/*
 * Old_SlaterIntegrals.cpp
 *
 *  Created on: Sep 19, 2018
 *      Author: brass
 */

SSlaterIntegral *GetSlaterIntegrals(double **basisfuns, unsigned int NBasisFuns, const double *x,
                                    unsigned int NGridPoints,const int **tupel, unsigned int NIntegrals)
{

    SSlaterIntegral *si;
    unsigned int kMax, n, i;
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

    for ( i = 0; i <= kMax; i++)
    {
    	free(kernels[i]);
    }
    free(kernels);
    for ( i = 0; i < NBasisFuns; i++)
    {
    	free(dualfuns[i]);
    }
    free(dualfuns);

    return si;
}

SSlaterIntegral *GetSlaterIntegrals(double **basisfuns, unsigned int NBasisFuns, const double *x,
                                    unsigned int NGridPoints,vector<int> indToK, unsigned int &NIntegrals)
{

    SSlaterIntegral *si;
    unsigned int kMax, n, i;
    double **kernels, **dualfuns;

    n = 2* NGridPoints;

    dualfuns = CreateDualFunctions(NBasisFuns, NGridPoints, basisfuns,x);

    si = InitializeIntegrals(indToK, NIntegrals,&kMax);

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
    for ( i = 0; i <= kMax; i++)
    {
    	free(kernels[i]);
    }
    free(kernels);
    for ( i = 0; i < NBasisFuns; i++)
    {
    	free(dualfuns[i]);
    }
    free(dualfuns);

    return si;
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
            //printf("Fk:\t%e\n",Fk);
            integrals[i].val[j] = Fk;

        }
    }
}

double **CreateDualFunctions(const unsigned int NBasisFuns, const unsigned int NGridPoints,
                             double **basisfuns, const double *x)
{
	uint32_t i;
    double **dualfuns, *weights;

    weights = CreateWeights(NGridPoints,x);
    //printf("%e\t%e\n",weights[2499],weights[4999]);
    dualfuns = (double **) malloc(NBasisFuns*sizeof(double *));
    if (!dualfuns) printf("malloc fail at dualfuns in CreateDualFunctions. size = %u\n",NBasisFuns);
    for ( i = 0; i < NBasisFuns; i++ )
    {
        dualfuns[i] = (double *) malloc(2*NGridPoints*sizeof(double));
        if (!dualfuns[i]) printf("malloc fail at dualfuns[%d] in CreateDualFunctions. size = %u\n",i,2*NGridPoints);
        Dual(2*NGridPoints, basisfuns[i], weights, dualfuns[i]);
    }
    free(weights);
    return dualfuns;
}

void Dual(const int n, const double *f, const double *weights, double *fdual)
{
    vecMul(n, f, weights, fdual);
}
