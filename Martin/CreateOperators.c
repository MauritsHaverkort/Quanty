#include "../Martin/SlaterIntegrals.c"
#include "../Martin/DiracHamiltonian.c"
#include "../Martin/ContinuumStates.c"
//#include "mkl.h"

void CreateOperatorCoulombRelativisticOneIntegral(uint32_t const NFermion,int const kappa1,int const kappa2, int const kappa3,int const kappa4,uint16_t const *Jindex1,
                                                  uint16_t const *Jindex2,uint16_t const *Jindex3, uint16_t const *Jindex4,double const *Fk, const unsigned int n1,
												  const unsigned int n2, const unsigned int n3, const unsigned int n4, OperatorType *const op)
{
    OperatorType optemp;
    
    CreateOperatorCoulombFourRelativisticShells(NFermion,kappa1,kappa2,kappa3,kappa4, Jindex1,Jindex2,Jindex3,Jindex4,Fk,op);
    
    if (n1 != n3)
    {
        CreateOperatorCoulombFourRelativisticShells(NFermion,kappa3,kappa2,kappa1,kappa4, Jindex3,Jindex2,Jindex1,Jindex4,Fk,&optemp);
        OperatorAdd(op, optemp, 1., 0.);
    }
    
    if (n2 != n4)
    {
        CreateOperatorCoulombFourRelativisticShells(NFermion,kappa1,kappa4,kappa3,kappa2, Jindex1,Jindex4,Jindex3,Jindex2,Fk,&optemp);
        OperatorAdd(op, optemp, 1., 0.);
        
        if (n1 != n3)
        {
            CreateOperatorCoulombFourRelativisticShells(NFermion,kappa3,kappa4,kappa1,kappa2, Jindex3,Jindex4,Jindex1,Jindex2,Fk,&optemp);
            OperatorAdd(op, optemp, 1., 0.);
        }
    }
    
    if ( (n1 != n2) || (n3 != n4) )
    {
        CreateOperatorCoulombFourRelativisticShells(NFermion,kappa2,kappa1,kappa4,kappa3, Jindex2,Jindex1,Jindex4,Jindex3,Fk,&optemp);
        OperatorAdd(op, optemp, 1., 0.);
    
        if (n1 != n3)
        {
            CreateOperatorCoulombFourRelativisticShells(NFermion,kappa2,kappa3,kappa4,kappa1, Jindex2,Jindex3,Jindex4,Jindex1,Fk,&optemp);
            OperatorAdd(op, optemp, 1., 0.);
        }
        
        if (n2 != n4)
        {
            CreateOperatorCoulombFourRelativisticShells(NFermion,kappa4,kappa1,kappa2,kappa3, Jindex4,Jindex1,Jindex2,Jindex3,Fk,&optemp);
            OperatorAdd(op, optemp, 1., 0.);
            if (n1 != n3)
            {
                CreateOperatorCoulombFourRelativisticShells(NFermion,kappa4,kappa3,kappa2,kappa1, Jindex4,Jindex3,Jindex2,Jindex1,Fk,&optemp);
                OperatorAdd(op, optemp, 1., 0.);
            }
        }
    }

}

uint32_t CreateOperatorFastCoulombRelativistic(uint32_t const NFermion,double **basisfuns, uint32_t NbasisFuns, const double *grid, uint32_t NGridpoints, uint32_t **tuple,
                                               uint32_t NIntegrals,uint16_t **indices, OperatorType *op)
{

    OperatorType optemp;
    SSlaterIntegral *integrals;

    integrals = GetSlaterIntegrals(basisfuns,NbasisFuns,grid, NGridpoints, (const int **) tuple, NIntegrals);

    int *k, i;
    unsigned int *n;
    double *Fk;

    k = integrals[0].tau.kappa;
    n = integrals[0].tau.n;
    Fk = integrals[0].val;

    CreateOperatorCoulombRelativisticOneIntegral(NFermion,k[0],k[1],k[2],k[3],indices[n[0]],indices[n[1]],indices[n[2]],indices[n[3]],Fk,n[0],n[1],n[2],n[3],op);

    for( i = 1; i < NIntegrals; i++)
    {
        k = integrals[i].tau.kappa;
        n = integrals[i].tau.n;
        Fk = integrals[i].val;

        CreateOperatorCoulombRelativisticOneIntegral(NFermion,k[0],k[1],k[2],k[3],indices[n[0]],indices[n[1]],indices[n[2]],indices[n[3]],Fk,n[0],n[1],n[2],n[3],&optemp);
        OperatorAdd(op, optemp, 1., 0.);

    }
    for ( i = 0; i < NIntegrals; i++)
    {
    	free(integrals[i].val);
    	free(integrals[i].k);
    }
    free(integrals);

    return 0;
}

uint32_t CreateOperatorFullCoulombRelativistic(uint32_t const NFermion,double **basisfuns, uint32_t NbasisFuns, const double *grid, uint32_t NGridpoints, uint32_t **tuple,
                                               uint32_t NIntegrals,uint16_t **indices, OperatorType *op)
{
    
    OperatorType optemp;
    SSlaterIntegral *integrals;
    //double t;
  
  
    BeginTiming("GetSlaterIntegrals", &GlobalTime);

    integrals = GetSlaterIntegrals(basisfuns,NbasisFuns,grid, NGridpoints, (const int **) tuple, NIntegrals);

    EndTiming("GetSlaterIntegrals", &GlobalTime);
  
    int *k, i;
    unsigned int *n;
    double *Fk;
    

  
    BeginTiming("CreateOperatorCoulomb", &GlobalTime);
    k = integrals[0].tau.kappa;
    n = integrals[0].tau.n;
    Fk = integrals[0].val;

    CreateOperatorCoulombFourRelativisticShells(NFermion,k[0],k[1],k[2],k[3], indices[n[0]],indices[n[1]],indices[n[2]],indices[n[3]],Fk,op);
    
    for( i = 1; i < NIntegrals; i++)
    {
        k = integrals[i].tau.kappa;
        n = integrals[i].tau.n;
        Fk = integrals[i].val;
        CreateOperatorCoulombFourRelativisticShells(NFermion,k[0],k[1],k[2],k[3], indices[n[0]],indices[n[1]],indices[n[2]],indices[n[3]],Fk,&optemp);
        OperatorAdd(op, optemp, 1., 0.);
    
    }
    EndTiming("CreateOperatorCoulomb", &GlobalTime);
    PrintTiming(GlobalTime);
    for ( i = 0; i < NIntegrals; i++)
    {
    	free(integrals[i].val);
    	free(integrals[i].k);
    }
    free(integrals);
    return 0;
}

typedef struct{
	unsigned int *indices;
	uint32_t n;
	int kappa;
}SKappaGroup;

uint32_t CreateOperatorDirac(uint32_t const NFermion,const unsigned int NBasisFuns, double **basisfuns, double **derivfuns, const double *x, uint32_t NGridpoints,
		SKappaGroup *groups,uint32_t Ngroups,uint32_t **indices, int Z, OperatorType *op)
{
	double **diracElements, **dualfuns;

	unsigned int *mj1,*mj2,n,k,j,i;
	uint16_t oadd[2] = { 0, 0 };

	dualfuns = CreateDualFunctions(NBasisFuns, NGridpoints,basisfuns,x);

	InitOperator(op);
	InitOperatorFixedLength(op,2,8);
	op->NFermion = NFermion;
	op->NWord = Min(op->NFermion, 1 + (op->NFermion - 1) / 8);

	for ( n = 0; n < Ngroups; n++)
	{
		diracElements = CalculateDiracMatrixElements(basisfuns,derivfuns,dualfuns,groups[n].indices,groups[n].n,x,NGridpoints,groups[n].kappa,Z);
		//PrintMartix2(diracElements,groups[n].n);

		for ( i = 0; i < groups[n].n; i++)
		{
			for ( j = i; j < groups[n].n; j++)
			{

				mj1 = indices[groups[n].indices[i]];
				mj2 = indices[groups[n].indices[j]];
				for ( k = 0; k < 2*abs(groups[n].kappa); k++)
				{

					oadd[1] = CrAnMask +  mj2[k];
					oadd[0] = mj1[k];
					RealOperatorAddLaderNormalOrder(op, oadd, 2, diracElements[i][j], &QDetNotOrderedFcFaBcBa);

					if ( i != j)
					{
						oadd[1] = CrAnMask +  mj1[k];
						oadd[0] = mj2[k];
						RealOperatorAddLaderNormalOrder(op, oadd, 2, diracElements[i][j], &QDetNotOrderedFcFaBcBa);
					}
				}
			}

			free(diracElements[i]);
		}
		free(diracElements);
	}



	return 0;
}

void PrintMartix(double *M, int n)
{
	int i,j;
	for ( i = 0; i < n; i++)
	{
		for ( j = 0; j < n; j++)
		{
			printf("%f\t",M[i*n + j]);
		}
		printf("\n");
	}
	fflush(stdout);
}

void PrintMartix2(double **M, int n)
{
	int i,j;
	for ( i = 0; i < n; i++)
	{
		for ( j = 0; j < n; j++)
		{
			printf("%f\t",M[i][j]);
		}
		printf("\n");
	}
	fflush(stdout);
}

void PrintArray(double* a, int n)
{
	int i;
	for ( i = 0; i < n; i++)
	{
		printf("%.15e\n", a[i]);
	}
}

double* CreateBaseForContinuum( int NBasisFuns,double **basisfuns,double **derivfuns,double *x,uint32_t NGridpoints,SKappaGroup group,int Ncont,
								double **contfuns, double w)
{
	int l,k,dim,info;

/*	const int LAPACK_ROW_MAJOR = 101;
	const int CblasRowMajor= 101;
	const int blasLeft = 141;
	const int CblasRight = 142;
	const int CblasUpper = 121;
	const int CblasNoTrans = 111;
	const int CblasTrans = 112;
	const int CblasNonUnit = 131;
*/
	/* Fuer Simon*/
	char blasRight = 'R';
	char blasLeft  = 'L';
	char blasUpper = 'U';
	char blasTrans = 'T';
	char blasNoTrans = 'N';
	char blasNonUnit = 'N';

	double *T,*g, **dualfuns, **dualcont, *Tsub, *tau,*d,*e;

	k = group.kappa;
	l = 0;

	dim = Ncont + group.n;


	dualfuns = CreateDualFunctions(NBasisFuns,NGridpoints,basisfuns,x);

	dualcont = CreateDualFunctions(Ncont,NGridpoints,contfuns,x);

	T = KineticEnergyCB(basisfuns,derivfuns, dualfuns,dualcont, Ncont,group.indices, group.n,x,NGridpoints,k,w);

	g = MetricCB(basisfuns,dualfuns,group.indices, group.n,dualcont,Ncont,NGridpoints);

	char cU = 'U';
	char cN = 'N';
	char cV = 'V';
	double EINS = 1;
	//int info;
	dpotrf(&cU,&dim,g,&dim,&info);
	dtrtri(&cU,&cN,&dim,g,&dim,&info);
	dtrmm(&blasLeft,&blasUpper,&blasTrans,&blasNonUnit,&dim,&dim,&EINS,g,&dim,T,&dim);
	dtrmm(&blasRight,&blasUpper,&blasNoTrans,&blasNonUnit,&dim,&dim,&EINS,g,&dim,T,&dim);
/*	LAPACKE_dpotrf(LAPACK_ROW_MAJOR,'U',dim,g,dim);
	LAPACKE_dtrtri(LAPACK_ROW_MAJOR,'U','N',dim,g,dim);
	cblas_dtrmm(CblasRowMajor,CblasLeft,CblasUpper,CblasTrans,CblasNonUnit,dim,dim,1.,g,dim,T,dim);
	cblas_dtrmm(CblasRowMajor,CblasRight,CblasUpper,CblasNoTrans,CblasNonUnit,dim,dim,1.,g,dim,T,dim);
*/
	Tsub = SubMatrix(T,dim,group.n,dim-1,group.n,dim-1);

	d = (double *) malloc(Ncont*sizeof(double));
	if (!d) printf("malloc fail at d in CreateBaseForContinuum. size = %d\n",Ncont);
	e = (double *) malloc((Ncont-1)*sizeof(double));
	if (!e) printf("malloc fail at e in CreateBaseForContinuum. size = %d\n",Ncont-1);
	tau = (double *) malloc(Ncont*sizeof(double));
	if (!tau) printf("malloc fail at tau in CreateBaseForContinuum. size = %d\n",Ncont);

	double *work = (double *) malloc(Ncont*sizeof(double));
	int lwork = Ncont;

	dsytrd(&cU,&Ncont,Tsub,&Ncont,d,e,tau,work,&lwork,&info);
	if ( 0 != info) printf("error: argument %d of LAPACKE_dsytrd",-info);
	dorgtr(&cU,&Ncont,Tsub,&Ncont,tau,work,&lwork,&info);
	if ( 0 != info) printf("error: argument %d of LAPACKE_dorgtr",-info);
	dsteqr(&cV,&Ncont,d,e,Tsub,&Ncont,work,&info);
	if ( 0 != info) printf("error: argument %d of LAPACKE_dsteqr",-info);
/*	info = LAPACKE_dsytrd(LAPACK_ROW_MAJOR, 'U',Ncont,Tsub,Ncont,d,e,tau);
	if ( 0 != info) printf("error: argument %d of LAPACKE_dsytrd",-info);
	info = LAPACKE_dorgtr(LAPACK_ROW_MAJOR,'U',Ncont,Tsub,Ncont,tau);
	if ( 0 != info) printf("error: argument %d of LAPACKE_dorgtr",-info);
	info = LAPACKE_dsteqr(LAPACK_ROW_MAJOR,'V',Ncont,d,e,Tsub,Ncont);
	if ( 0 != info) printf("error: argument %d of LAPACKE_dsteqr",-info);
*/
	CalcContinuumWFs(basisfuns, group.indices, group.n, contfuns, Ncont, NGridpoints, g, Tsub);

	free(T);
	free(g);
	free(work);
	int i;
	for ( i = 0; i < NBasisFuns; i++)
	{
		free(dualfuns[i]);
	}
	free(dualfuns);
	for ( i = 0; i < Ncont; i++)
	{
		free(dualcont[i]);
	}
	free(dualcont);
	free(Tsub);
	free(tau);
	free(e);


	return d;

}

/*
	printf("%f\n",NormDO(n,k,l,w));
	printf("%f\n",C1DO(n,k,l,w));
	printf("%f\n",C2DO(n,k,l,w));
	printf("%f\n",EnergyDO(n,k,l,w));

	m = 1;
	printf("%f\n",EKinDO(n,m,k,l,w));
	m = 2;
	printf("%f\n",EKinDO(n,m,k,l,w));
	m = 3;
	printf("%f\n",EKinDO(n,m,k,l,w));
	n = 2;
	m = 1;
	printf("%f\n",EKinDO(n,m,k,l,w));

	n=1;
	printf("%f\n",Binomial(n + l + 1./2.,n));
	printf("%f\n",Binomial(n + l - 1./2.,n));
	printf("%f\n",Binomial(n-1 + l + 3./2.,n-1));
	*/
