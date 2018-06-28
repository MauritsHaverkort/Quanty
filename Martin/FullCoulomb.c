#include "../Martin/SlaterIntegrals.c"

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
    //double t;

    //t = dsecnd();
    integrals = GetSlaterIntegrals(basisfuns,NbasisFuns,grid, NGridpoints, (const int **) tuple, NIntegrals);
    //printf("%f\n",dsecnd()-t);

    int *k, i;
    unsigned int *n;
    double *Fk;

    //t = dsecnd();

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
    //printf("%f\n",dsecnd()-t);
    return 0;
}

uint32_t CreateOperatorFullCoulombRelativistic(uint32_t const NFermion,double **basisfuns, uint32_t NbasisFuns, const double *grid, uint32_t NGridpoints, uint32_t **tuple,
                                               uint32_t NIntegrals,uint16_t **indices, OperatorType *op)
{
    
    OperatorType optemp;
    SSlaterIntegral *integrals;
    //double t;
  
  
  BeginTiming("GetSlaterIntegrals", &GlobalTime);
    //t = dsecnd();
    integrals = GetSlaterIntegrals(basisfuns,NbasisFuns,grid, NGridpoints, (const int **) tuple, NIntegrals);
    //printf("%f\n",dsecnd()-t);
  EndTiming("GetSlaterIntegrals", &GlobalTime);
  
    int *k, i;
    unsigned int *n;
    double *Fk;
    
    //t = dsecnd();
  
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
    return 0;
}


