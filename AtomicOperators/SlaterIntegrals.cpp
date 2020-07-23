#include "../AtomicOperators/SlaterIntegrals.h"
#include "../Martin/Hartree-Fock/Hartree-Fock.h"

SSlaterIntegral *GetSlaterIntegrals(const vector<vector<double> > &basisfuns,const vector<double> &x,const vector<int> &indToK, unsigned int &NIntegrals, bool Qcontinuum)
{

  if(GlobalTimeOn&0x100) BeginTiming("Slater Integrals", &GlobalTime);
    SSlaterIntegral *si;
    unsigned int kMax, n, i;
    double **kernels;
    vector<vector<double> > dualfuns;
    unsigned int NGridPoints = (unsigned int)x.size();


    n = 2* NGridPoints;

    CreateDualFunctions(basisfuns,x,dualfuns);

    si = InitializeIntegrals(indToK, NIntegrals,&kMax,Qcontinuum);

    kernels = CreateAllKernels(kMax, NGridPoints, &x[0]);

    CalculateSlaterIntegrals(NIntegrals, si,n,kernels,basisfuns,dualfuns);

    for ( i = 0; i <= kMax; i++)
    {
    	free(kernels[i]);
    }
    free(kernels);

  if(GlobalTimeOn&0x100) EndTiming("Slater Integrals", &GlobalTime);
    return si;
}

SSlaterIntegral *GetSlaterIntegralsTwoFilled(const vector<vector<double> > &basisfuns,const vector<double> &x,const vector<int> &indToK, unsigned int &NIntegrals)
{

  if(GlobalTimeOn&0x100) BeginTiming("Slater Integrals", &GlobalTime);
    SSlaterIntegral *si;
    unsigned int kMax, n, i;
    double **kernels;
    vector<vector<double> > dualfuns;
    unsigned int NGridPoints = (unsigned int)x.size();


    n = 2* NGridPoints;

    CreateDualFunctions(basisfuns,x,dualfuns);

    si = InitializeIntegralsTwoFilled(indToK, NIntegrals,&kMax);

    kernels = CreateAllKernels(kMax, NGridPoints, &x[0]);

    CalculateSlaterIntegrals(NIntegrals, si,n,kernels,basisfuns,dualfuns);

    for ( i = 0; i <= kMax; i++)
    {
    	free(kernels[i]);
    }
    free(kernels);

  if(GlobalTimeOn&0x100) EndTiming("Slater Integrals", &GlobalTime);
    return si;
}

SSlaterIntegral *GetSlaterIntegralsOneFilled(const vector<vector<double> > &basisfuns,const vector<double> &x,const vector<int> &indToK1,const vector<int> &indToK2, unsigned int &NIntegrals)
{

  if(GlobalTimeOn&0x100) BeginTiming("Slater Integrals", &GlobalTime);
    SSlaterIntegral *si;
    unsigned int kMax, n, i;
    double **kernels;
    vector<vector<double> > dualfuns;
    unsigned int NGridPoints = (unsigned int)x.size();


    n = 2* NGridPoints;

    CreateDualFunctions(basisfuns,x,dualfuns);

    si = InitializeIntegralsOneFilled(indToK1,indToK2, NIntegrals,&kMax);

    kernels = CreateAllKernels(kMax, NGridPoints, &x[0]);

    CalculateSlaterIntegrals(NIntegrals, si,n,kernels,basisfuns,dualfuns);

    for ( i = 0; i <= kMax; i++)
    {
    	free(kernels[i]);
    }
    free(kernels);

  if(GlobalTimeOn&0x100) EndTiming("Slater Integrals", &GlobalTime);
    return si;
}

SSlaterIntegral *GetSlaterIntegrals(const vector<Interpolation::InterpolatingFunction *> &ipolfuns,const vector<int> &indToK, unsigned int &NIntegrals, bool Qcontinuum)
{

	uint32_t quadDegree = 3;
	vector<double> x = ipolfuns[0]->GetX();
	vector<double> knots = Quad::GetKnots(quadDegree,x);
	uint32_t Nknots = (uint32_t)knots.size();
	uint32_t Nfuns = (uint32_t)ipolfuns.size();
    vector<vector<double> > basisfuns = vector<vector<double> >(Nfuns);
    vector<vector<double> > dualfuns  = vector<vector<double> >(Nfuns);
    vector<double> weights = Quad::GetWeights(quadDegree,x);
    for (uint32_t i = 0; i < Nfuns; i++)
    {
    	basisfuns[i] = (*ipolfuns[i])(knots);
    	dualfuns[i] = weights * basisfuns[i];
    }


	SSlaterIntegral *si;
	unsigned int kMax;
	vector<vector<double> > kernels;

    si = InitializeIntegrals(indToK, NIntegrals,&kMax,Qcontinuum);
    //clock_t t = clock();
    CreateAllKernels(kMax, knots,kernels);
    //cout << "Kernels: " << (clock()-t)*1.0/CLOCKS_PER_SEC << endl;
    //t = clock();
    CalculateSlaterIntegrals(NIntegrals, si,Nknots,kernels,basisfuns,dualfuns);
    //cout << "CalcSlater: " << (clock()-t)*1.0/CLOCKS_PER_SEC << endl;

    return si;
}

void GetSlaterIntegrals(const vector<BSpline::TBSpline *> &G,const vector<BSpline::TBSpline *> &F, const vector<int> &indToK,uint32_t Nweights, vector<SSlaterIntegral> &si)
{

	uint32_t quadDegree = 3;
	uint32_t kMax = 0;
	vector<vector<double> > kernels;
	uint32_t Nfuns = (uint32_t) G.size();
	uint32_t n = G[0]->degree;


	kMax = 2 * Abs(MaxAbs(indToK)) -1;
	vector<vector<double> > rk1 = vector<vector<double> >(kMax+1);
	vector<vector<double> > rk2 = vector<vector<double> >(kMax+1);
	vector<vector<double> > rkD = vector<vector<double> >(kMax+1);
	vector<double> BSval;
	BeginTiming("Slater Init", &GlobalTime);
//	uint32_t Nweights = n +10;
	HartreeFock::PutBSplinesOnQuadGrid(G[0]->knots,n,Nweights,BSval);
	for (uint32_t k = 0; k <= kMax; k++)
	{
		HartreeFock::CalcOffDiagSlaterInts(G[0]->knots,BSval,Nweights,rk1[k],rk2[k],k,n);
		HartreeFock::CalcDiagSlaterInts(G[0]->knots,BSval,Nweights,rkD[k],k,n);
	}
	EndTiming("Slater Init", &GlobalTime);
	//printf("init done\n");
	BeginTiming("Slater Calculate", &GlobalTime);
#pragma omp parallel for
	for (uint32_t i1 = 0; i1 < Nfuns; i1++)
	{
		SSlaterIntegral Integral;
		int k_i1 = indToK[i1];
		for (uint32_t j1 = i1; j1 < Nfuns; j1++)
		{
			int k_j1 = indToK[j1];
			for (uint32_t i2 = i1; i2 < Nfuns; i2++)
			{
				int k_i2 = indToK[i2];
				for (uint32_t j2 = j1; j2 < Nfuns; j2++)
				{
					int k_j2 = indToK[j2];
				    uint32_t k0 = (uint32_t) fmax( abs( abs(k_i1)-abs(k_i2) ), abs( abs(k_j1)-abs(k_j2) ) );
				    uint32_t k1 = (uint32_t) fmin( abs(k_i1)+abs(k_i2), abs(k_j1)+abs(k_j2) ) -1;

				    if ( (abs(k_i1)+abs(k_i2) + (Sign(k_i1)+Sign(k_i2))/2 + k0) % 2 == 0)
				    {
				        k0++;
				    }
				    if (k0 <= k1)
				    {
						Integral.tau.kappa[0] = k_i1;
						Integral.tau.kappa[1] = k_j1;
						Integral.tau.kappa[2] = k_i2;
						Integral.tau.kappa[3] = k_j2;
						Integral.tau.n[0] = i1;
						Integral.tau.n[1] = j1;
						Integral.tau.n[2] = i2;
						Integral.tau.n[3] = j2;
						Integral.Nval = (k1-k0)/2 +1;
						Integral.val = (double*) malloc(Integral.Nval * sizeof(double)) ;
						Integral.k = (unsigned int*) malloc(Integral.Nval * sizeof(unsigned int)) ;
						for (uint32_t k = k0; k <= k1; k +=2)
						{
							uint32_t idx = (k - k0)/2;
							Integral.k[idx] = k;
							Integral.val[idx] = HartreeFock::AssembleSlaterInts2(G[i1]->y,G[j1]->y,G[i2]->y,G[j2]->y,rk1[k],rk2[k],rkD[k],n)
											  + HartreeFock::AssembleSlaterInts2(G[i1]->y,F[j1]->y,G[i2]->y,F[j2]->y,rk1[k],rk2[k],rkD[k],n)
											  + HartreeFock::AssembleSlaterInts2(F[i1]->y,G[j1]->y,F[i2]->y,G[j2]->y,rk1[k],rk2[k],rkD[k],n)
											  + HartreeFock::AssembleSlaterInts2(F[i1]->y,F[j1]->y,F[i2]->y,F[j2]->y,rk1[k],rk2[k],rkD[k],n);
						}
#pragma omp critical
						si.push_back(Integral);
				    }
				}
			}
		}
	}
	EndTiming("Slater Calculate", &GlobalTime);
}

void GetSlaterIntegrals(const vector<BSpline::TBSpline *> &G, const vector<int> &indToK,uint32_t Nweights, vector<SSlaterIntegral> &si)
{

	uint32_t quadDegree = 3;
	uint32_t kMax = 0;
	vector<vector<double> > kernels;
	uint32_t Nfuns = (uint32_t) G.size();
	uint32_t n = G[0]->degree;
	SSlaterIntegral Integral;

	kMax = 2 * Abs(MaxAbs(indToK)) -1;
	vector<vector<double> > rk1 = vector<vector<double> >(kMax+1);
	vector<vector<double> > rk2 = vector<vector<double> >(kMax+1);
	vector<vector<double> > rkD = vector<vector<double> >(kMax+1);
	vector<double> BSval;
	BeginTiming("Slater Init", &GlobalTime);
//	uint32_t Nweights = n +10;
	HartreeFock::PutBSplinesOnQuadGrid(G[0]->knots,n,Nweights,BSval);
	for (uint32_t k = 0; k <= kMax; k++)
	{
		HartreeFock::CalcOffDiagSlaterInts(G[0]->knots,BSval,Nweights,rk1[k],rk2[k],k,n);
		HartreeFock::CalcDiagSlaterInts(G[0]->knots,BSval,Nweights,rkD[k],k,n);
	}
	EndTiming("Slater Init", &GlobalTime);
//	printf("init done\n Nfuns = %u\n",Nfuns);
	BeginTiming("Slater Calculate", &GlobalTime);
	for (uint32_t i1 = 0; i1 < Nfuns; i1++)
	{
		int k_i1 = indToK[i1];
		for (uint32_t j1 = i1; j1 < Nfuns; j1++)
		{
			int k_j1 = indToK[j1];
			for (uint32_t i2 = i1; i2 < Nfuns; i2++)
			{
				int k_i2 = indToK[i2];
				for (uint32_t j2 = j1; j2 < Nfuns; j2++)
				{
					int k_j2 = indToK[j2];
				    uint32_t k0 = (uint32_t) fmax( abs( abs(k_i1)-abs(k_i2) ), abs( abs(k_j1)-abs(k_j2) ) );
				    uint32_t k1 = (uint32_t) fmin( abs(k_i1)+abs(k_i2), abs(k_j1)+abs(k_j2) ) -1;

				    if ( (abs(k_i1)+abs(k_i2) + (Sign(k_i1)+Sign(k_i2))/2 + k0) % 2 == 0)
				    {
				        k0++;
				    }
				    if (k0 <= k1)
				    {
						//printf("i1 = %u\t21 = %u\ti2 = %u\tj2 = %u\tk0 %u\tk1 = %u\n",i1,j1,i2,j2,k0,k1);
				    	Integral.tau.kappa[0] = k_i1;
						Integral.tau.kappa[1] = k_j1;
						Integral.tau.kappa[2] = k_i2;
						Integral.tau.kappa[3] = k_j2;
						Integral.tau.n[0] = i1;
						Integral.tau.n[1] = j1;
						Integral.tau.n[2] = i2;
						Integral.tau.n[3] = j2;
						Integral.Nval = (k1-k0)/2 +1;
						Integral.val = (double*) malloc(Integral.Nval * sizeof(double)) ;
						Integral.k = (unsigned int*) malloc(Integral.Nval * sizeof(unsigned int)) ;
						for (uint32_t k = k0; k <= k1; k +=2)
						{
							uint32_t idx = (k - k0)/2;
							Integral.k[idx] = k;
							Integral.val[idx] = HartreeFock::AssembleSlaterInts2(G[i1]->y,G[j1]->y,G[i2]->y,G[j2]->y,rk1[k],rk2[k],rkD[k],n);
						}
						si.push_back(Integral);
				    }
				}
			}
		}
	}
	EndTiming("Slater Calculate", &GlobalTime);
}

void CalculateSlaterIntegrals(const int NIntegrals, SSlaterIntegral *integrals, const int dim,
                              double **kernels,const vector<vector<double> > &basisfuns,const vector<vector<double> > &dualfuns)
{

#pragma omp parallel for
    for(int i = 0; i < NIntegrals; i++)
    {

    	unsigned int *n = integrals[i].tau.n;
        for(int j = 0; j < integrals[i].Nval; j++)
        {

            int k = integrals[i].k[j];
            double Fk = SlaterIntegral(dim, &(dualfuns[n[0]][0]), &(dualfuns[n[1]][0]),
                                &(basisfuns[n[2]][0]),&(basisfuns[n[3]][0]),kernels[k]);
            integrals[i].val[j] = Fk;

        }
    }
}

void CalculateSlaterIntegrals(const int NIntegrals, SSlaterIntegral *integrals, uint32_t dim,
                              vector<vector<double> > &kernels, vector<vector<double> > &basisfuns, vector<vector<double> > &dualfuns)
{

#pragma omp parallel for
    for(int i = 0; i < NIntegrals; i++)
    {

    	unsigned int *n = integrals[i].tau.n;
        for(int j = 0; j < integrals[i].Nval; j++)
        {

        	int k = integrals[i].k[j];
            double Fk = SlaterIntegral(dim, &(dualfuns[n[0]][0]), &(dualfuns[n[1]][0]),
                                &(basisfuns[n[2]][0]),&(basisfuns[n[3]][0]),&(kernels[k][0]));
            integrals[i].val[j] = Fk;

        }
    }
}

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
  const unsigned int INCXY=1;
  const double beta=0;
  double *f13, *f24, *Mf, result;

  f13 = (double *) malloc(n*sizeof(double));
  f24 = (double *) malloc(n*sizeof(double));
  Mf = (double *) malloc(n*sizeof(double));
  if (!Mf || !f24 || !f13) printf("malloc fail in SlaterIntegral. size = %d",n);
    
  vecMul(n, f1, f3, f13);
  vecMul(n, f2, f4, f24);

  dsymv(&UPLO,(unsigned int *) &n, &alpha, Mr1r2,(unsigned int *) &n, f24, &INCXY, &beta, Mf, &INCXY);
  result = ddot((unsigned int *) &n, f13, &INCXY, Mf, &INCXY);
  
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
    if (!kernels) printf("malloc fail at kernels in CreateAllKernels. size = %u\n",(int) kMax+1);
    n = 2 * NGridPoints;
    
    for ( k = 0; k <= kMax; k++ )
    {
    	kernels[k] = (double *) malloc(n*n*sizeof(double));
    	if (!(kernels[k])) printf("malloc fail at kernels[%u] in CreateAllKernels. size = %u\n",k,n*n);
        CreateKernel(NGridPoints, x, kernels[k], k);
    }
    return kernels;
}

void CreateAllKernels(const unsigned int kMax,const vector<double> &x, vector<vector<double> > &kernels)
{
    kernels = vector<vector<double> >(kMax+1);

    for (uint32_t k = 0; k <= kMax; k++ )
    {
    	CreateNonRelKernel(x, k, kernels[k]);
    }
}

void CreateKernel(const unsigned int NGridPoints, const double *x, double *kernel, const unsigned int k)
{
    uint32_t i, j, n;
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

void CreateNonRelKernel(const vector<double> &x, const uint32_t k,vector<double> &kernel)
{
    double val;
    uint32_t n = (uint32_t)x.size();
    kernel = vector<double>(n*n);

    for(uint32_t i = 0; i < n; i++)
    {
        for(uint32_t j = i; j < n; j++)
        {
            if (i == 0)
            {
                val = 0;
            }
            else
            {
                val = pow(x[i],k)/pow(x[j],k+1);
            }

            kernel[j+n*i] = val;
            kernel[i+n*j] = val;
        }
    }

}

double *CreateWeights(const unsigned int NGridPoints, const double *x)
{
	uint32_t i;
	double x0 = 0;

    double *weights = (double *) malloc(NGridPoints*2*sizeof(double));
    if (!weights) printf("malloc fail at weights in CreateWeights. size = %u\n",2*NGridPoints);
    
    for( i = 0; i < NGridPoints-1; i++)
    {
        weights[i] = 2 * (x[i] - x0);
        weights[i+NGridPoints] = 2 * (x[i] - x0);
        //printf("weights[%d]: %f\n",i,weights[i]);
        x0 = x[i] + weights[i] / 2;
        //printf("weights %d\n",NGridPoints*2-(i+NGridPoints)-1);
        //fflush(stdout);
    }
    weights[NGridPoints-2] += x[i]-x0;
    weights[2*NGridPoints-2] += x[i]-x0;
    weights[2*NGridPoints-1] = 0.;
    weights[NGridPoints-1] = 0.;

    return weights;
}

vector<double> CreateWeights(const unsigned int NGridPoints, const vector<double> &x)
{
	uint32_t i;
    double x0 = 0;

    vector<double> weights;

    for( i = 0; i < NGridPoints-1; i++)
    {
        weights.insert(weights.end(), 2 * (x[i] - x0));
        x0 = x[i] + weights[i] / 2;
        //printf("weights %d\n",NGridPoints*2-(i+NGridPoints)-1);
        //fflush(stdout);
    }
    weights.insert(weights.end(), 0.);
    weights[NGridPoints-2] += x[i]-x0;
    for( i = 0; i < NGridPoints-1; i++)
    {
        weights.insert(weights.end(), 2 * (x[i] - x0));
        x0 = x[i] + weights[i] / 2;
        //printf("weights %d\n",NGridPoints*2-(i+NGridPoints)-1);
        //fflush(stdout);
    }
    weights[2*NGridPoints-2] += x[i]-x0;
    weights.insert(weights.end(), 0.);


    return weights;
}

vector<double> CreateWeights(const vector<double> &x)
{
	uint32_t i;
    double x0 = 0;
    uint32_t NGridPoints = (uint32_t) x.size();
    vector<double> weights;

    for( i = 0; i < NGridPoints-1; i++)
    {
        weights.insert(weights.end(), 2 * (x[i] - x0));
        x0 = x[i] + weights[i] / 2;
        //printf("weights %d\n",NGridPoints*2-(i+NGridPoints)-1);
        //fflush(stdout);
    }
    weights.insert(weights.end(), 0.);
    weights[NGridPoints-2] += x[i]-x0;
    return weights;
}

vector<double> Dual(const vector<double> &f, const vector<double> &w)
{
	vector<double> g = vector<double>(f.size());
	if (f.size() != w.size()) printf("error in Dual: %lu = size(f) != size(w) = %lu", f.size(),w.size());
	for (uint32_t i = 0; i < (uint32_t) f.size(); i++)
	{
		g[i] = f[i]*w[i];
	}
	return g;
}

void CreateDualFunctions(const vector<vector<double> > &basisfuns,const vector<double> &x, vector<vector<double> > &dualfuns,bool QRelativistic)
{
	uint32_t i;
    vector<double> weights;
    unsigned int NBasisFuns = (unsigned int )basisfuns.size();
    unsigned int NGridPoints = (unsigned int )x.size();

    if (QRelativistic)
    {
    	weights = CreateWeights(NGridPoints,x);
    }
    else
    {
    	weights = CreateWeights(x);
    }

    dualfuns = vector<vector<double> >(NBasisFuns);
    for ( i = 0; i < NBasisFuns; i++ )
    {
        dualfuns[i] = Dual(basisfuns[i], weights);
    }
}

SSlaterIntegral *InitializeIntegrals(const vector<int> &indToK, unsigned int &NIntegrals, unsigned int *kMax, bool Qcontinuum)
{
    unsigned int i, k;
    SSlaterIntegral *integrals;

    std::map<myarray,myarray> tupel;
    if (Qcontinuum)
    {
    	tupel = CreateContinuumTuple(indToK);
    }
    else
    {
    	tupel = CreateSlaterTuple(indToK);
    }
    NIntegrals = (unsigned int)tupel.size();

    integrals = (SSlaterIntegral *) malloc((NIntegrals)*sizeof(SSlaterIntegral));
    i = 0;
    for( std::map<myarray,myarray>::iterator it = tupel.begin(); it != tupel.end(); ++it)
    {
    	integrals[i].tau.n[0] = it->first.val[0];
    	integrals[i].tau.n[1] = it->first.val[1];
    	integrals[i].tau.n[2] = it->first.val[2];
    	integrals[i].tau.n[3] = it->first.val[3];

    	//printf("%u\t%u\t%u\t%u\n",integrals[i].tau.n[0],integrals[i].tau.n[1],integrals[i].tau.n[2],integrals[i].tau.n[3]);

    	integrals[i].tau.kappa[0] = it->second.val[0];
    	integrals[i].tau.kappa[1] = it->second.val[1];
    	integrals[i].tau.kappa[2] = it->second.val[2];
    	integrals[i].tau.kappa[3] = it->second.val[3];
    													// TODO i could be larger than NIntegrals ???
    	i++;
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

SSlaterIntegral *InitializeIntegralsTwoFilled(const vector<int> &indToK, unsigned int &NIntegrals, unsigned int *kMax)
{
    unsigned int i, k;
    SSlaterIntegral *integrals;

    std::map<myarray,myarray> tupel;
   
    tupel = CreateSlaterTupleTwoFilled(indToK);
    NIntegrals = (unsigned int)tupel.size();

    integrals = (SSlaterIntegral *) malloc((NIntegrals)*sizeof(SSlaterIntegral));
    i = 0;
    for( std::map<myarray,myarray>::iterator it = tupel.begin(); it != tupel.end(); ++it)
    {
    	integrals[i].tau.n[0] = it->first.val[0];
    	integrals[i].tau.n[1] = it->first.val[1];
    	integrals[i].tau.n[2] = it->first.val[2];
    	integrals[i].tau.n[3] = it->first.val[3];


    	integrals[i].tau.kappa[0] = it->second.val[0];
    	integrals[i].tau.kappa[1] = it->second.val[1];
    	integrals[i].tau.kappa[2] = it->second.val[2];
    	integrals[i].tau.kappa[3] = it->second.val[3];

    	i++;
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

SSlaterIntegral *InitializeIntegralsOneFilled(const vector<int> &indToK1,const vector<int> &indToK2, unsigned int &NIntegrals, unsigned int *kMax)
{
    unsigned int i, k;
    SSlaterIntegral *integrals;

    std::map<myarray,myarray> tupel;
   
    tupel = CreateSlaterTupleOneFilled(indToK1,indToK2);
    NIntegrals = (unsigned int)tupel.size();

    integrals = (SSlaterIntegral *) malloc((NIntegrals)*sizeof(SSlaterIntegral));
    i = 0;
    for( std::map<myarray,myarray>::iterator it = tupel.begin(); it != tupel.end(); ++it)
    {
    	integrals[i].tau.n[0] = it->first.val[0];
    	integrals[i].tau.n[1] = it->first.val[1];
    	integrals[i].tau.n[2] = it->first.val[2];
    	integrals[i].tau.n[3] = it->first.val[3];


    	integrals[i].tau.kappa[0] = it->second.val[0];
    	integrals[i].tau.kappa[1] = it->second.val[1];
    	integrals[i].tau.kappa[2] = it->second.val[2];
    	integrals[i].tau.kappa[3] = it->second.val[3];

    	i++;
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
    unsigned int k0=0, k1=0, *k;
    
    
    k0 = (unsigned int) fmax( abs( abs(kappa[0])-abs(kappa[2]) ), abs( abs(kappa[1])-abs(kappa[3]) ) );
    k1 = (unsigned int) fmin( abs(kappa[0])+abs(kappa[2]), abs(kappa[1])+abs(kappa[3]) ) -1;
    
    if ( (abs(kappa[0])+abs(kappa[2]) + (Sign(kappa[0])+Sign(kappa[2]))/2 + k0) % 2 == 0)
    {
        k0++;
    }
    *n = (k1-k0)/2 + 1;
    if ( k0 > k1 )
    {
      k=NULL;
      printf("Warning: k0(%u) > k1(%u)\n This has to be fixed...",k0,k1);
    }
    else
    {
        k = (unsigned int *) malloc((*n) * sizeof(unsigned int));
        for(int i = 0; i < *n; i++)
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





















