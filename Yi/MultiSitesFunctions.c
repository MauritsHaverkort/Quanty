/* Settings of DMFT calculation */

#define NSiteMax 16
typedef struct MultiSitesDMFTSettings
{
  uint32_t NSite;                             /* Number of cluster sites */
  char     G0FileName[NSiteMax][256];         /* list of G0 filenames */
  char     RootFileName[NSiteMax][256];       /* root of GBath, GCluster and Sigma */
  char     SpecRootFileName[256];             /* filename of spectral functions in spectral representation */
  char     WaveFunctionRootFileName[256];     /* filename of wave-function  */
  char     LocalHamiltonianFilename[256];     /* Input file for the local Hamiltonian */
  char     GSymFilename[256];                 /* File specifying the symmetry of G */
  double   Nelectron;                         /* Filling of the impurity */
  double   mu;                                /* Chemical potential of the lattice Greensfunction */
  double   mu0;                               /* Chemical potential of G0 */
  double   Teff;                              /* effective temperature when determening the relation between chemical potential and filling for the bath */
  double   T;                                 /* real temperature in the cluster not jet implemented */
  uint32_t NmuSym;                            /* 0000 0000 mu0 and mu from input file are used as starting values                                 */
                                              /*           mu is optimized during the loops in order to fix the occupation to Nelectron           */
                                              /* 0000 0001 Before start chemical potential mu and mu0 are read from G0 and Sigma input files      */
                                              /*           These override the values in the input file. (Otherwise values in input file are used  */
                                              /* 0000 0010 Nelectron is ignored and the calculation is preformed at fixed mu                      */
                                              /*           otherwise mu is optimized every loop                                                   */
  uint32_t iStart;                            /* start value of loop number */
  uint32_t Ni;                                /* number of loops to do */
  uint32_t Dyson;                             /* 0 no Dyson like perturbation - 1 full perturbation - 2 ... */
  uint32_t SymOpt;                            /* 0000 0000 no aditional symmetrization */
                                              /* 0000 0001 An additional Operator of (U+1) * S^2 is added to the Hamiltonian to stabilize the singlet state */
                                              /* 0000 0010 The Greenfunction is averaged such that G(w) = - G(-w)^*                                         */
                                              /* 0000 0100 The Ground-state wavefunction is taken to have a fully occupied / unoccupied bath                */
                                              /* 0000 1000 mu is forced to be <n> U/2 - this also sets NmuSym to NmuSym & 10                                */
  uint32_t NGCluster;                         /* Number of poles in the resulting cluster green function */
  uint32_t NGBath;                            /* Number of poles in the resulting bath green function */
  double   Mix;                               /* mixing of old into new function */
  double   U;                                 /* valence - valence Coulomb repulsion */
  double   J;                                 /* Hund's coupling */
  uint32_t NBath;                             /* Number of Bath states to include */
  uint32_t Spectra;                           /* 0000 no spectra will be calculated at the end of the DMFT loop                 */
                                              /* 0001 Magnetic susceptibility will be calculated (nup-ndn)                      */
                                              /* 0010 occupation will be calculated (nup+ndn)                                   */
                                              /* 0100 XAS will be calculated with NQ values for core - valence interaction      */
                                              /* 1000 core PES will be calculated with NQ values for core - valence interaction */
  uint32_t NQ;                                /* number of XAS spectra */
  double   *Q;                                /* NQ values of core - valence interaction strengths */
  uint32_t NE;                                /* number of points in the spectrum */
  double   Emin;                              /* minimum energy */
  double   Emax;                              /* maximum energy */
  double   Gamma;                             /* Broadening only important for the Spec file */
  double   zeroLancEnergy;                    /* minimum energy difference for convergence of total energy to be different from zero */
  double   zeroLancWavefunction;              /* minnimum pre-factor in wave-function to be different from zero */
  uint32_t NTriDiag;                          /* number of tri-diagonalization steps in Lanczos */
  uint32_t NBitsKey;                          /* 2^NBitsKey == lenght of hash key */
}
MultiSitesDMFTSettings;


/*
 *   Reads the settings for a DMFTSettings
 */
uint32_t StreamReadMultiSitesDMFTSettings(MultiSitesDMFTSettings * O, FILE *stream)
{
  char comments[1000];
  uint32_t i;

  fscanf(stream,"%i",&(O->NSite));                   fgets(comments,1000,stream);
  for (i=0; i<O->NSite; i++)
  {
    fscanf(stream,"%s",(O->G0FileName[i]));          fgets(comments,1000,stream);
  }
  for (i=0; i<O->NSite; i++)
  {
    fscanf(stream,"%s",(O->RootFileName[i]));        fgets(comments,1000,stream);
  }
  fscanf(stream,"%s",(O->SpecRootFileName));         fgets(comments,1000,stream);
  fscanf(stream,"%s",(O->WaveFunctionRootFileName)); fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->Nelectron));              fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->mu0));                    fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->mu));                     fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->Teff));                   fgets(comments,1000,stream);
  fscanf(stream,"%i",&(O->NmuSym));                  fgets(comments,1000,stream);
  fscanf(stream,"%i",&(O->iStart));                  fgets(comments,1000,stream);
  fscanf(stream,"%i",&(O->Ni));                      fgets(comments,1000,stream);
  fscanf(stream,"%i",&(O->Dyson));                   fgets(comments,1000,stream);
  fscanf(stream,"%i",&(O->SymOpt));                  fgets(comments,1000,stream);
  if((O->SymOpt)&8) O->NmuSym&=2;
  fscanf(stream,"%i",&(O->NGCluster));               fgets(comments,1000,stream);
  fscanf(stream,"%i",&(O->NGBath));                  fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->Mix));                    fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->U));                      fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->J));                      fgets(comments,1000,stream);
  fscanf(stream,"%i",&(O->NBath));                   fgets(comments,1000,stream);
  fscanf(stream,"%i",&(O->Spectra));
  O->NQ=1;
  if((O->Spectra)&240)
  {
    fscanf(stream,"%i",&(O->NQ));
    O->Q=malloc(O->NQ*sizeof(double));
    for(i=0;i<O->NQ;i++)
    {
      fscanf(stream,"%lf",&(O->Q[i]));
    }
    if(O->NQ==0)
    {
      O->Spectra&=!((uint32_t)12);
    }
  }                                                  fgets(comments,1000,stream);
  fscanf(stream,"%i",&(O->NE));                      fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->Emin));                   fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->Emax));                   fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->Gamma));                  fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->zeroLancEnergy));         fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->zeroLancWavefunction));   fgets(comments,1000,stream);
  fscanf(stream,"%i",&(O->NTriDiag));                fgets(comments,1000,stream);
  fscanf(stream,"%i",&(O->NBitsKey));                fgets(comments,1000,stream);

  return 0;
}


uint32_t FileReadMultiSitesDMFTSettings(MultiSitesDMFTSettings * O, const char *const flname)
{
  FILE *stream;
  uint32_t retval;
  if((stream=fopen(flname,"r"))==NULL)
  {
    printf("ERROR\nCan't open %s for reading.\n",flname);
    return 1;
  }
  retval=StreamReadMultiSitesDMFTSettings(O,stream);
  fclose(stream);

  return retval;
}


uint32_t StreamPrintMultiSitesDMFTSettings(MultiSitesDMFTSettings O, FILE *stream)
{
  uint32_t i;

  fprintf(stream,"\nSettings of self-consistent DMFT calculation\n");
  fprintf(stream,"# of sites               : %i\n", O.NSite);
  for (i=0; i<O.NSite; i++)
  {
    fprintf(stream,"G0FileName[%i]            : %s\n", i, O.G0FileName[i]);
  }
  for (i=0; i<O.NSite; i++)
  {
    fprintf(stream,"RootFileName[%i]          : %s\n", i, O.RootFileName[i]);
  }
  fprintf(stream,"SpecRootFileName         : %s\n", O.SpecRootFileName);
  fprintf(stream,"WaveFunctionRootFileName : %s\n", O.WaveFunctionRootFileName);
  fprintf(stream,"number of electrons      : %22.15E\n", O.Nelectron);
  fprintf(stream,"Chemical potential of G0 : %22.15E\n", O.mu0);
  fprintf(stream,"Chemical potential of G  : %22.15E\n", O.mu);
  fprintf(stream,"Effective temperature    : %22.15E\n", O.Teff);
  fprintf(stream,"                         converged results do not depent on Teff, convergence speed and stability do\n");
  fprintf(stream,"Options for Chemical pot : %i\n", O.NmuSym);
  if((O.NmuSym)&1 )fprintf(stream,"                         ( 1) mu and mu0 are read from G0 and Sigma, values in imput ignored;\n");
  if((O.NmuSym)&2 )fprintf(stream,"                         ( 2) Calculations done at fixed mu, Nelectron ignored;\n");
  fprintf(stream,"iStart                   : %i\n", O.iStart);
  fprintf(stream,"Ni                       : %i\n", O.Ni);
  fprintf(stream,"Dyson perturbation       : %i\n", O.Dyson);
  fprintf(stream,"Symmetrizaton options    : %i\n", O.SymOpt);
  if((O.SymOpt)&1) fprintf(stream,"                         ( 1) HAnd_grd += (U+1)S^2;\n");
  if((O.SymOpt)&2) fprintf(stream,"                         ( 2) G(w) +=-G(-w)^*; G(w)/=2;\n");
  if((O.SymOpt)&4) fprintf(stream,"                         ( 4) psi_cluster lives on 2 sites per orbital\n");
  if((O.SymOpt)&8) fprintf(stream,"                         ( 8) mu = U/2\n");
  fprintf(stream,"NGCluster                : %i\n", O.NGCluster);
  fprintf(stream,"NGBath                   : %i\n", O.NGBath);
  fprintf(stream,"Mix                      : %24.17E\n", O.Mix);
  fprintf(stream,"U                        : %24.17E\n", O.U);
  fprintf(stream,"J                        : %24.17E\n", O.J);
  fprintf(stream,"NBath                    : %i\n", O.NBath);
  fprintf(stream,"Spectra                  : %i\n", O.Spectra);
  if((O.Spectra)&1  ) fprintf(stream,"                        (  1) Sz;\n");
  if((O.Spectra)&2  ) fprintf(stream,"                        (  2) N;\n");
  if((O.Spectra)&4  ) fprintf(stream,"                        (  4) N-1;\n");
  if((O.Spectra)&8  ) fprintf(stream,"                        (  8) N_up N_dn;\n");
  if((O.Spectra)&16 ) fprintf(stream,"                        ( 16) core PES;\n");
  if((O.Spectra)&32 ) fprintf(stream,"                        ( 32) XAS;\n");
  if((O.Spectra)&64 ) fprintf(stream,"                        ( 64) XAS - d0;\n");
  if((O.Spectra)&128) fprintf(stream,"                        (128) XAS - d1;\n");
  if((O.Spectra)&240)
  {
    fprintf(stream,"  NQ                     : %i\n", O.NQ);
    for(i=0; i<O.NQ; i++)
    {
      fprintf(stream,"  Q[%3i]                 : %24.17E\n", i, O.Q[i]);
    }
  }
  fprintf(stream,"NE                       : %i\n", O.NE);
  fprintf(stream,"Emin                     : %24.17E\n", O.Emin);
  fprintf(stream,"Emax                     : %24.17E\n", O.Emax);
  fprintf(stream,"Gamma                    : %24.17E\n", O.Gamma);
  fprintf(stream,"zeroLancEnergy           : %24.17E\n", O.zeroLancEnergy);
  fprintf(stream,"zeroLancWaveFunction     : %24.17E\n", O.zeroLancWavefunction);
  fprintf(stream,"NTriDiag                 : %i\n", O.NTriDiag);
  fprintf(stream,"NBitsKey                 : %i\n", O.NBitsKey);

  return 0;
}


uint32_t PrintMultiSitesDMFTSettings(MultiSitesDMFTSettings O)
{
  return StreamPrintMultiSitesDMFTSettings(O, stdout);
}


/* read initial mu values from file */
uint32_t FileReadMu(const char * filename, MultiSitesDMFTSettings O, double * mu)
{
  FILE *stream;
  uint32_t i;

  stream=fopen(filename, "r");
  for (i=0; i<O.NSite; i++)
  {
    fscanf(stream, "%lf", &mu[i]);
  }

  return 0;
}

double read_oneline(const char * filename)
{
  FILE *stream;
  double num=0.0;
  stream=fopen(filename, "r");
  fscanf(stream,"%lf", &num);

  return num;
}


/*  Intput:
 *    Structrue containing the DMFT settings
 *  Output:
 *    G0 in anderson format
 */
uint32_t initializeMultiSitesDMFTReadG0(MultiSitesDMFTSettings *const O,
                                        AndersonMatrixType *const G0And)
{
  uint32_t i, flag;

  flag = 0;

  #pragma omp parallel for
  for (i=0; i<O->NSite; i++)
  {
    #pragma omp flush (flag)
    if (!flag)
    {
      if(O->NmuSym&1)
      {
        /* read chemical potential from G0, overwrite values set in O */
        if (FileReadAndersonMatrixTypeUndeterminedPlusMu(&G0And[i], O->G0FileName[i], &(O->mu0)))
        {
          printf("Fatal Error in main\n Failed to read G0 from file %s\n", O->G0FileName[i]);
          flag = 1;
        }

      }
      else
      {
        if (FileReadAndersonMatrixTypeUndetermined(&G0And[i], O->G0FileName[i]))
        {
          printf("Fatal Error in main\n Failed to read G0 from file %s\n", O->G0FileName[i]);
          flag = 1;
        }
      }
      strcpy(G0And[i].Name, O->G0FileName[i]);

      if((O->SymOpt)&2)
      {
        printf("Symmetrize G0(w) -= G0(-w)^*\n");
        SymmetrizeAndersonMatrix(&G0And[i]);
      }

      printf("\nRead G0 from file \"%s\"\n", O->G0FileName[i]);
      G0And[i].mu=O->mu0;
      PrintAndersonMatrixEssentials(G0And[i]);
    }
  }

  if (flag)
    return 1;

  return 0;
}


/*  Intput:
 *    Structrue containing the DMFT settings
 *  Output:
 *    GBath in anderson and list of poles format
 */
uint32_t initializeMultiSitesDMFTReadGBath(MultiSitesDMFTSettings *const O,
                                           AndersonMatrixType *const GBathupAnd,
                                           AndersonMatrixType *const GBathdnAnd,
                                           ListOfPoles *const GBathupLP,
                                           ListOfPoles *const GBathdnLP)
{
  uint32_t i, flag;
  char FileNameup[MAXStrLength], FileNamedn[MAXStrLength];
  char GBathupFileName[MAXStrLength], GBathdnFileName[MAXStrLength];

  flag = 0;

  #pragma omp parallel for
  for (i=0; i<O->NSite; i++)
  {
    #pragma omp flush (flag)
    if (!flag)
    {
      sprintf(FileNameup, "%sup%03d", O->RootFileName[i], O->iStart);
      sprintf(FileNamedn, "%sdn%03d", O->RootFileName[i], O->iStart);
      if (file_exists(FileNameup) && file_exists(FileNamedn))
      {
        if ((FileReadListOfPolesUndetermined(&GBathupLP[i], FileNameup)) ||
            FileReadListOfPolesUndetermined(&GBathdnLP[i], FileNamedn))
        {
          printf("Fatal Error in main\n Failed to read G0A from file %s (%s)\n",
                 FileNameup, FileNamedn);
          flag = 1;
        }
        sprintf(GBathupFileName, "%sGBath-up", O->RootFileName[i]);
        sprintf(GBathdnFileName, "%sGBath-dn", O->RootFileName[i]);
        strcpy(GBathupLP[i].Name, GBathupFileName);
        strcpy(GBathdnLP[i].Name, GBathdnFileName);

        if (ListOfPolesToAndersonMatrix(GBathupLP[i], &GBathupAnd[i]) ||
            ListOfPolesToAndersonMatrix(GBathdnLP[i], &GBathdnAnd[i]))
        {
          printf("Error in ListOfPolesToAndersonMatrix in initializeMultiSitesDMFTReadGBath\n");
          flag = 1;
        }
        strcpy(GBathupAnd[i].Name, GBathupFileName);
        strcpy(GBathdnAnd[i].Name, GBathdnFileName);
      }
      else
      {
        printf("File for bath Greenfunction %s (%s) does not exist\n",
               FileNameup, FileNamedn);
        printf("Start from the situation where GBath=G0\n");
        initializeMultiSitesDMFTReadG0(O, &(GBathupAnd[i]));
        initializeMultiSitesDMFTReadG0(O, &(GBathdnAnd[i]));
        AndersonMatrixToListOfPoles(GBathupAnd[i], &GBathupLP[i]);
        AndersonMatrixToListOfPoles(GBathdnAnd[i], &GBathdnLP[i]);
        sprintf(GBathupFileName, "%sGBath-up", O->RootFileName[i]);
        sprintf(GBathdnFileName, "%sGBath-dn", O->RootFileName[i]);
        strcpy(GBathupAnd[i].Name, GBathupFileName);
        strcpy(GBathdnAnd[i].Name, GBathdnFileName);
        strcpy(GBathupLP[i].Name, GBathupFileName);
        strcpy(GBathdnLP[i].Name, GBathdnFileName);
      }
    }
  }

  if (flag)
    return 1;

  return 0;
}


/*  Intput:
 *    Structrue containing the DMFT settings
 *  Output:
 *    Sigma as a list of poles
 */
uint32_t initializeMultiSitesDMFTReadSigma(MultiSitesDMFTSettings *const O,
                                           ListOfPoles *const Sigmaup,
                                           ListOfPoles *const Sigmadn)
{
  uint32_t i, flag;
  char FileNameup[MAXStrLength], FileNamedn[MAXStrLength];
  char SigmaupFileName[MAXStrLength], SigmadnFileName[MAXStrLength];

  flag = 0;

  #pragma omp parallel for
  for (i=0; i<O->NSite; i++)
  {
    sprintf(FileNameup, "%sSigmaup%03d", O->RootFileName[i], O->iStart);
    sprintf(FileNamedn, "%sSigmadn%03d", O->RootFileName[i], O->iStart);

    #pragma omp flush (flag)
    if (!flag)
    {
      if (file_exists(FileNameup) && file_exists(FileNamedn))
      {
        if (FileReadListOfPolesUndetermined(&Sigmaup[i],FileNameup) ||
            FileReadListOfPolesUndetermined(&Sigmadn[i],FileNamedn))
        {
          printf("Fatal Error in main\n Failed to read Sigma from file %s (%s)!\n",
                 FileNameup, FileNamedn);
          flag = 1;
        }
        sprintf(SigmaupFileName, "%sSigma-up", O->RootFileName[i]);
        sprintf(SigmadnFileName, "%sSigma-dn", O->RootFileName[i]);
        strcpy(Sigmaup[i].Name, SigmaupFileName);
        strcpy(Sigmadn[i].Name, SigmadnFileName);
        printf("\nRead Sigma from file \"%s (%s)\"\n",
               FileNameup, FileNamedn);
        PrintListOfPolesEssentials(Sigmaup[i]);
        PrintListOfPolesEssentials(Sigmadn[i]);
      }
      else
      {
        printf("File for the selfenergy %s (%s) does not exist\n",
               FileNameup, FileNamedn);
        printf("Start from the situation where Sigma=0\n");
        Sigmaup[i].N=1;
        InitListOfPoles(&Sigmaup[i]);
        Sigmaup[i].a0=0;
        Sigmaup[i].a[0]=0;
        Sigmaup[i].b[0]=0;

        Sigmadn[i].N=1;
        InitListOfPoles(&Sigmadn[i]);
        Sigmadn[i].a0=0;
        Sigmadn[i].a[0]=0;
        Sigmadn[i].b[0]=0;
        sprintf(SigmaupFileName, "%sSigma-up", O->RootFileName[i]);
        sprintf(SigmadnFileName, "%sSigma-dn", O->RootFileName[i]);
        strcpy(Sigmaup[i].Name, SigmaupFileName);
        strcpy(Sigmadn[i].Name, SigmadnFileName);
        Sigmaup[i].mu=0;
        Sigmadn[i].mu=0;
        PrintListOfPolesEssentials(Sigmaup[i]);
        PrintListOfPolesEssentials(Sigmadn[i]);
      }
    }
  }

  if (flag)
    return 1;

  return 0;
}


/* create initial wave-functions
 with impurity sites and BathV fully occupied
 */
uint32_t CreateWaveFunctionMultiSitesDMFT(WaveFunctionType *Psi, uint32_t NSite,
                                          uint32_t nBathV, uint32_t nBathC,
                                          uint32_t NBitsKey)
{
  uint32_t i, j;
  uint8_t *Determinant;
  char Buffer[MAXStrLength];
  char *WaveFunctionName = Buffer;

  Psi->NBitsKey = NBitsKey;
  Psi->QComplex = 0;
  Psi->NBoson = 0;
  Psi->NFermion = 2 * NSite * ( 2 + nBathC + nBathV);
  if(WaveFunctionInit(Psi))
  {
    printf("WaveFunctionInit failed in CreateWaveFunctionMultiSitesDMFT\n");fflush(stdout);
    return 1;
  }

  WaveFunctionName += sprintf(WaveFunctionName, "DMFT %d ", NSite);
  for (i=0; i<NSite; i++)
  {
    WaveFunctionName += sprintf(WaveFunctionName, "%d", nBathV);
  }
  WaveFunctionName += sprintf(WaveFunctionName, " ");
  for (i=0; i<NSite; i++)
  {
    WaveFunctionName += sprintf(WaveFunctionName, "%d", nBathC);
  }
  sprintf(Psi->Name, "%s", Buffer);

  Determinant=calloc(Psi->NWord, sizeof(uint8_t));
  for (i=0; i<2*NSite; i++)
  {
    // fill impurity site up and dn
    OperatorCr(Determinant, i, Psi->NFermion);
  }

  for (i=0; i<NSite; i++)
  {
    for (j=0; j<2*nBathV; j++)
    {
      // BathV up and dn
      OperatorCr(Determinant, j+4*NSite+2*(nBathV+nBathC)*i, Psi->NFermion);
    }
  }
  RealWaveFunctionAddElement(Psi, Determinant, 1.0);
  free(Determinant);

  return 0;
}


/* Read the name of wave function and extract the information:
 *  number of site:                                      NSite
 *  number of valence bath of each site as an array:     nBathV
 *  number of conduction bath of each site as an array:  nBathC
 */

uint32_t ReadMultiSitesWavefunctionName(MultiSitesDMFTSettings *const O,
                                        WaveFunctionType *const Psi,
                                        uint32_t * NSite,
                                        uint32_t * nBathV,
                                        uint32_t * nBathC)
{
  uint32_t i;
  char WaveFunctionSpec[MAXStrLength], BathV[MAXStrLength], BathC[MAXStrLength];

  sscanf(Psi->Name, "%s %i %s %s", WaveFunctionSpec, NSite, BathV, BathC);
  if (*NSite != O->NSite || *NSite!=(int) strlen(BathV) || *NSite!=(int) strlen(BathC))
  {
    printf("ERROR in ReadMultiSitesWavefunctionName: NSite differ.\n");
    return 1;
  }

  printf("Wavefunction has the following specifications:\n");
  printf("The number of sites is            : %i\n", *NSite);
  printf("The numbers of valence bath are   : %s\n", BathV);
  printf("The numbers of conduction bath are: %s\n\n", BathC);
  for (i=0; i<*NSite; i++)
  {
    nBathV[i]=BathV[i]-'0';
    nBathC[i]=BathC[i]-'0';
  }

  return 0;
}


/*  Intput:
 *    Structrue containing the DMFT settings
 *  Output:
 *    Psi, a starting wave-function for the DMFT loops
 */
uint32_t initializeMultiSitesDMFTReadPsi(MultiSitesDMFTSettings *const O,
                                         WaveFunctionType *const Psi)
{
  char FileName[MAXStrLength];
  uint32_t i, PsiNSite;
  uint32_t PsiNBathV[O->NSite], PsiNBathC[O->NSite];
  int TotBathV=0, TotBathC=0;

  sprintf(FileName, "%s%03d", O->WaveFunctionRootFileName, O->iStart);
  if(file_exists(FileName))
  {
    Psi->NBitsKey = O->NBitsKey;
    if (FileReadWaveFunctionPlainBinary(Psi, FileName))
    {
      printf("Fatal Error in main\n Failed to read Psi from file %s !\n", FileName);
      return 1;
    }
    printf("\nStarting wave function read from the file %s\n", FileName);

    if (ReadMultiSitesWavefunctionName(O, Psi, &PsiNSite, PsiNBathV, PsiNBathC))
      return CreateWaveFunctionMultiSitesDMFT(Psi, O->NSite, O->NBath, O->NBath, O->NBitsKey);

    printf("Wave-function \"%s\" found has %i Fermions in the basis\n",
           Psi->Name, Psi->NFermion);
    for (i=0;i<PsiNSite; i++)
    {
      TotBathV += PsiNBathV[i];
      TotBathC += PsiNBathC[i];
    }
    if(4*PsiNSite+2*(TotBathV+TotBathC) != Psi->NFermion)
    {
      printf("The total number of Fermions is incompatible witht the bath sizes\n");
      printf("This wave-function can not be used to start from\nStarted from new function\n");
      return CreateWaveFunctionMultiSitesDMFT(Psi, O->NSite, O->NBath, O->NBath, O->NBitsKey);
    }
    else
      return 0;
  }
  else
  {
    printf("\nWave-function file: '%s' does not exist\nStarted from new function\n", FileName);
    return CreateWaveFunctionMultiSitesDMFT(Psi, O->NSite, O->NBath, O->NBath, O->NBitsKey);
  }
}


/*  Intput:
 *    filename containing the DMFT settings
 *
 *  Output:
 *    DMFT settings obtained from the input file
 *    G0     - returned in Anderson Matrix Format
 *    Sigma  - returned as list of poles
 *    Psi    - read from file or set to starting function
 *
 */
uint32_t initializeMultiSitesDMFT(const char *const flname,
                                  MultiSitesDMFTSettings *const O,
                                  AndersonMatrixType *const G0And,
                                  ListOfPoles *const SigmaupLP,
                                  ListOfPoles *const SigmadnLP,
                                  WaveFunctionType *const PsiStart )
{
  if (FileReadMultiSitesDMFTSettings(O, flname))
  {
    printf("Fatal Error in initializeMultiSitesDMFT\n");
    printf("Failed to read DMFT Settings from file %s\n", flname);
    return 1;
  }
  if (initializeMultiSitesDMFTReadG0(O, G0And))
  {
    printf("Fatal Error in initializeMultiSitesDMFTReadG0\n");
    return 1;
  }
  if (initializeMultiSitesDMFTReadSigma(O, SigmaupLP, SigmadnLP))
  {
    printf("Fatal Error in initializeMultiSitesDMFTReadSigma\n");
    return 1;
  }
  if (initializeMultiSitesDMFTReadPsi(O, PsiStart))
  {
    printf("Fatal Error in initializeDMFTReadPsi\n");
    return 1;
  }

  return PrintMultiSitesDMFTSettings(*O);
}


/* Test if the wave function has the correct dimensions
 specified by DMFT settings */
uint32_t MultiSitesTestPsiDimensions(MultiSitesDMFTSettings *const O,
                                     BADoubleTriDiagonalMatrixType *const GBathupBATrunc,
                                     BADoubleTriDiagonalMatrixType *const GBathdnBATrunc,
                                     WaveFunctionType *const Psi)
{
  char Buffer[MAXStrLength];
  char *WaveFunctionName = Buffer;
  uint32_t PsiNSite, PsiNBathV[O->NSite], PsiNBathC[O->NSite];
  uint32_t i, j, k, StartIndex;
  WaveFunctionType PsiTemp;

  ReadMultiSitesWavefunctionName(O, Psi, &PsiNSite, PsiNBathV, PsiNBathC);

  for (i=O->NSite-1; 0-i<=0; i--)
  {
    if (GBathupBATrunc[i].con.N != GBathdnBATrunc[i].con.N ||
        GBathupBATrunc[i].val.N != GBathdnBATrunc[i].val.N)
    {
      printf("Bath site numbers for spin up and dn are different!\n");
      return 1;
    }

    /* Starting index of the corresponding valence bath */
    StartIndex = 4*O->NSite;
    for (j=0; j<i; j++)
    {
      StartIndex += 2*(PsiNBathV[j] + PsiNBathC[j]);
    }

    /* Comparing conduction baths */
    if (PsiNBathC[i] != GBathupBATrunc[i].con.N-1)
    {
      printf("\nBefore diagonalizing:"
             "\nThe wavefunction has %i conduction bath for site %i"
             "\nThe Anderson Hamiltonian %i\n",
             PsiNBathC[i], i, GBathupBATrunc[i].con.N-1);

      if(PsiNBathC[i] > GBathupBATrunc[i].con.N-1)
      {
        printf("Removing %i occupied fermions from the wave-function\n",
               PsiNBathC[i]-GBathupBATrunc[i].con.N+1);
        for(k=0; k<PsiNBathC[i]-GBathupBATrunc[i].con.N+1; k++)
        {
          WaveFunctionRemoveFermion(*Psi, &PsiTemp,
                                    StartIndex+2*PsiNBathV[i]+2*PsiNBathC[i]-k-1);
          WaveFunctionFree(Psi);
          WaveFunctionRemoveFermion(PsiTemp, Psi,
                                    StartIndex+2*PsiNBathV[i]+  PsiNBathC[i]-k-1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      else
      {
        printf("Adding %i unoccupied fermions to the wave-function\n",
               GBathupBATrunc[i].con.N-1-PsiNBathC[i]);
        for(k=0; k<GBathupBATrunc[i].con.N-1-PsiNBathC[i]; k++)
        {
          WaveFunctionAddFermion(*Psi, &PsiTemp,
                                 StartIndex+2*PsiNBathV[i]+2*PsiNBathC[i]+k, 0);
          WaveFunctionFree(Psi);
          WaveFunctionAddFermion(PsiTemp, Psi,
                                 StartIndex+2*PsiNBathV[i]+  PsiNBathC[i]+k, 0);
          WaveFunctionFree(&PsiTemp);
        }
      }

      PsiNBathC[i] = GBathupBATrunc[i].con.N-1;
    }

    /* Comparing valence baths */
    if (PsiNBathV[i] != GBathupBATrunc[i].val.N-1)
    {
      printf("\nBefore diagonalizing:"
             "\nThe wavefunction has %i valence bath for site %i"
             "\nThe Anderson Hamiltonian %i\n",
             PsiNBathV[i], i, GBathupBATrunc[i].val.N-1);

      if(PsiNBathV[i] > GBathupBATrunc[i].val.N-1)
      {
        printf("Removing %i occupied fermions from the wave-function\n",
               PsiNBathV[i]-GBathupBATrunc[i].val.N+1);
        for(k=0; k<PsiNBathV[i]-GBathupBATrunc[i].val.N+1; k++)
        {
          WaveFunctionRemoveFermion(*Psi, &PsiTemp,
                                    StartIndex+2*PsiNBathV[i]-k-1);
          WaveFunctionFree(Psi);
          WaveFunctionRemoveFermion(PsiTemp, Psi,
                                    StartIndex+  PsiNBathV[i]-k-1);
          WaveFunctionFree(&PsiTemp);
        }
      }
      else
      {
        printf("Adding %i unoccupied fermions to the wave-function\n",
               GBathupBATrunc[i].val.N-1-PsiNBathV[i]);
        for(k=0; k<GBathupBATrunc[i].val.N-1-PsiNBathV[i]; k++)
        {
          WaveFunctionAddFermion(*Psi, &PsiTemp,
                                 StartIndex+2*PsiNBathV[i]+k, 0);
          WaveFunctionFree(Psi);
          WaveFunctionAddFermion(PsiTemp, Psi,
                                 StartIndex+  PsiNBathV[i]+k, 0);
          WaveFunctionFree(&PsiTemp);
        }
      }

      PsiNBathV[i] = GBathupBATrunc[i].val.N-1;
    }
  }

  printf("Corrected basis!\n");
  PrintWaveFunctionReadableTruncated(*Psi, 0.1);

  /* update the wavefunction name */
  WaveFunctionName += sprintf(WaveFunctionName, "DMFT %d ", O->NSite);
  for (i=0; i<O->NSite; i++) {
    WaveFunctionName += sprintf(WaveFunctionName, "%d", PsiNBathV[i]);
  }
  WaveFunctionName += sprintf(WaveFunctionName, " ");
  for (i=0; i<O->NSite; i++) {
    WaveFunctionName += sprintf(WaveFunctionName, "%d", PsiNBathC[i]);
  }
  sprintf(Psi->Name, "%s", Buffer);

  return 0;
}


/*   These two functions creates the most straightforward Hamiltonian
 *            for multiple-site impurity models
 *
 * index:
 *
 *    impurity:
 *      i = 0   for site 0 spin up, i = 1     for site 0 spin dn
 *      i = 2   for site 1 spin up, i = 3     for site 1 spin dn
 *      i = 2*N for site N spin up, i = 2*N+1 for site N spin dn
 *
 *    direct bath:
 *      i = 2*NSite+0   for bath 0 spin up, 2*NSite+1     for bath 0 spin dn
 *      i = 2*NSite+2*N for bath N spin up, 2*NSite+2*N+1 for bath N spin dn
 *    rest of the bath:
 *      ordered as site N   val up, val dn, con up, con dn
 *                 site N+1 val up, val dn, con up, con dn
 *
 *    additional parameters:  J is the Hund coupling.
 */

uint32_t LocalHamiltonian(OperatorType *O,
                          BADoubleTriDiagonalMatrixType *Mup,
                          BADoubleTriDiagonalMatrixType *Mdn,
                          uint32_t NSite,
                          double U, double J, double mu)
{
  uint32_t i;
  uint16_t oadd[4];

  InitNamedOperator(O, "Local Hamiltonian");
  printf("\nLocal Hamiltonian for NSite = %i\n\n", NSite);
  O->NFermion = 0;
  for (i=0; i<NSite; i++) {
    O->NFermion += Mup[i].val.N + Mup[i].con.N + Mdn[i].val.N + Mdn[i].con.N;
  }
  O->NWord = Min(O->NFermion, 1+(O->NFermion-1)/8);
  InitOperatorFixedLength(O, 2, 8);
  InitOperatorFixedLength(O, 4, 8);

  /* Local Hamiltonian */

  /* This has to be changed for different models
   * here shows two settings, for single orbital and
   * two d-orbital with cubic symmetry
   */
  if (NSite==1)
  {
    oadd[3]=CrAnMask+0;
    oadd[2]=CrAnMask+1;
    oadd[1]=1;
    oadd[0]=0;
    RealOperatorAddLaderNormalOrder(O, oadd, 4, U, &QDetNotOrderedFcFaBcBa);
  }
  else if (NSite==2)
  {
    /* U */
    oadd[3] = CrAnMask+0;
    oadd[2] = CrAnMask+1;
    oadd[1] = 1;
    oadd[0] = 0;
    RealOperatorAddLaderNormalOrder(O, oadd, 4, U, &QDetNotOrderedFcFaBcBa);
    oadd[3] = CrAnMask+2;
    oadd[2] = CrAnMask+3;
    oadd[1] = 3;
    oadd[0] = 2;
    RealOperatorAddLaderNormalOrder(O, oadd, 4, U, &QDetNotOrderedFcFaBcBa);
    /* U' */
    oadd[3] = CrAnMask+0;
    oadd[2] = CrAnMask+3;
    oadd[1] = 3;
    oadd[0] = 0;
    RealOperatorAddLaderNormalOrder(O, oadd, 4, U-2*J, &QDetNotOrderedFcFaBcBa);
    oadd[3] = CrAnMask+1;
    oadd[2] = CrAnMask+2;
    oadd[1] = 2;
    oadd[0] = 1;
    RealOperatorAddLaderNormalOrder(O, oadd, 4, U-2*J, &QDetNotOrderedFcFaBcBa);
    /* U'' */
    oadd[3] = CrAnMask+0;
    oadd[2] = CrAnMask+2;
    oadd[1] = 2;
    oadd[0] = 0;
    RealOperatorAddLaderNormalOrder(O, oadd, 4, U-3*J, &QDetNotOrderedFcFaBcBa);
    oadd[3] = CrAnMask+1;
    oadd[2] = CrAnMask+3;
    oadd[1] = 3;
    oadd[0] = 1;
    RealOperatorAddLaderNormalOrder(O, oadd, 4, U-3*J, &QDetNotOrderedFcFaBcBa);
  }

  for (i=0; i<NSite; i++)
  {
    /* onsite mu */
    oadd[1] = CrAnMask+2*i;
    oadd[0] = 2*i;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, -mu, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask+2*i+1;
    oadd[0] = 2*i+1;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, -mu, &QDetNotOrderedFcFaBcBa);
  }
  PrintOperatorReadable(*O);

  return 0;
}

/* Create DCA cluster Hamiltonian for single orbital
 * with different cluster sizes (NSite)
 */
uint32_t DCAHamiltonian(OperatorType *O,
                        BADoubleTriDiagonalMatrixType *Mup,
                        BADoubleTriDiagonalMatrixType *Mdn,
                        uint32_t NSite,
                        double U, double mu)
{
  uint32_t i, j, k, l, NPair;
  uint16_t oadd[4];
  double kxs[NSite], kys[NSite], dkx, dky;

  InitNamedOperator(O, "DCA Hamiltonian");
  printf("\nCluster Hamiltonian for NSite = %i\n\n", NSite);
  O->NFermion = 0;
  for (i=0; i<NSite; i++)
  {
    O->NFermion += Mup[i].val.N + Mup[i].con.N + Mdn[i].val.N + Mdn[i].con.N;
  }
  O->NWord = Min(O->NFermion, 1+(O->NFermion-1)/8);
  InitOperatorFixedLength(O, 2, 8);
  InitOperatorFixedLength(O, 4, 8);

  switch (NSite)
  {
    case 2:
      kxs[0] =   0; kys[0] =   0;
      kxs[1] = 1.0; kys[1] = 1.0;
      break;

    case 4:
      kxs[0] =    0; kys[0] =    0;
      kxs[1] =  1.0; kys[1] =  1.0;
      kxs[2] =  1.0; kys[2] =    0;
      kxs[3] =    0; kys[3] =  1.0;
      break;

    case 8:
      kxs[0] =    0; kys[0] =    0;
      kxs[1] =  1.0; kys[1] =  1.0;
      kxs[2] =  1.0; kys[2] =    0;
      kxs[3] =    0; kys[3] =  1.0;
      kxs[4] =  0.5; kys[4] =  0.5;
      kxs[5] = -0.5; kys[5] =  0.5;
      kxs[6] = -0.5; kys[6] = -0.5;
      kxs[7] =  0.5; kys[7] = -0.5;
      break;

    default:
      printf("Error in DCAHamiltonian: Cluster size not implemented yet!\n");
      return 1;
      break;
  }

  NPair = 0;
  for (i=0; i<NSite; i++)
  {
    for (j=0; j<NSite; j++)
    {
      for (k=0; k<NSite; k++)
      {
        dkx = kxs[i] + kxs[j] - kxs[k];
        dky = kys[i] + kys[j] - kys[k];

        l = 0;
        while ((dkx!=kxs[l]-2 && dkx!=kxs[l] && dkx!=kxs[l]+2) ||
               (dky!=kys[l]-2 && dky!=kys[l] && dky!=kys[l]+2))
          l++;
        if (l<NSite)
        {
          oadd[3] = CrAnMask+2*i;
          oadd[2] = CrAnMask+2*j+1;
          oadd[1] = 2*k;
          oadd[0] = 2*l+1;
          RealOperatorAddLaderNormalOrder(O, oadd, 4, U/NSite, &QDetNotOrderedFcFaBcBa);
          NPair +=1;
        }
        else
        {
          printf("Error in DCAHamiltonian!\n");
          return 1;
        }
      }
    }
  }

  printf("There are %i interaction terms in the Hamiltonian\n", NPair);
    
  for (i=0; i<NSite; i++)
  {
    /* onsite mu */
    oadd[1] = CrAnMask+2*i;
    oadd[0] = 2*i;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, -mu, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask+2*i+1;
    oadd[0] = 2*i+1;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, -mu, &QDetNotOrderedFcFaBcBa);
   }

  return 0;
}


uint32_t MultiSitesBADoubleTriDiagonalMatrixToAndersonImpurityOperator(OperatorType *O,
                                                                       OperatorType *Oloc,
                                                                       BADoubleTriDiagonalMatrixType *Mup,
                                                                       BADoubleTriDiagonalMatrixType *Mdn,
                                                                       uint32_t NSite)
{
  uint32_t i, j, StartIndex;
  uint16_t oadd[4];
  InitNamedOperator(O, "Anderson Impurity Hamiltonian");

  O->NFermion = 0;
  for (i=0; i<NSite; i++) {
    O->NFermion += Mup[i].val.N + Mup[i].con.N + Mdn[i].val.N + Mdn[i].con.N;
  }
  O->NWord = Min(O->NFermion, 1+(O->NFermion-1)/8);
  InitOperatorFixedLength(O, 2, 8);
  InitOperatorFixedLength(O, 4, 8);

  for (i=0; i<NSite; i++)
  {
    /* onsite */
    oadd[1] = CrAnMask+2*i;
    oadd[0] = 2*i;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup[i].acvB, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask+2*i+1;
    oadd[0] = 2*i+1;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn[i].acvB, &QDetNotOrderedFcFaBcBa);

    /* direct bath site */
    oadd[1] = CrAnMask+2*NSite+2*i;
    oadd[0] = 2*NSite+2*i;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup[i].acvA,  &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask+2*i;
    oadd[0] = 2*NSite+2*i;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup[i].acvBA, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask+2*NSite+2*i;
    oadd[0] = 2*i;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup[i].acvBA, &QDetNotOrderedFcFaBcBa);

    oadd[1] = CrAnMask+2*NSite+2*i+1;
    oadd[0] = 2*NSite+2*i+1;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn[i].acvA,  &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask+2*i+1;
    oadd[0] = 2*NSite+2*i+1;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn[i].acvBA, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask+2*NSite+2*i+1;
    oadd[0] = 2*i+1;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn[i].acvBA, &QDetNotOrderedFcFaBcBa);

    /* rest of the bath */

    // val up
    StartIndex = 4*NSite;
    for (j=0; j<i; j++)
    {
      StartIndex += Mup[j].val.N + Mup[j].con.N + Mdn[j].val.N + Mdn[j].con.N - 4;
    }

    oadd[1] = CrAnMask+2*i;
    oadd[0] = StartIndex;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup[i].bvB0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask+StartIndex;
    oadd[0] = 2*i;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup[i].bvB0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask+2*NSite+2*i;
    oadd[0] = StartIndex;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup[i].bvA0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask+StartIndex;
    oadd[0] = 2*NSite+2*i;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup[i].bvA0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask+StartIndex;
    oadd[0] = StartIndex;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup[i].val.a[2], &QDetNotOrderedFcFaBcBa);

    for (j=0; j+2<Mup[i].val.N; j++)
    {
      oadd[1] = CrAnMask+StartIndex+j;
      oadd[0] = StartIndex+j+1;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup[i].val.b[j], &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask+StartIndex+j+1;
      oadd[0] = StartIndex+j;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup[i].val.b[j], &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask+StartIndex+j+1;
      oadd[0] = StartIndex+j+1;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup[i].val.a[j+1], &QDetNotOrderedFcFaBcBa);
    }

    // val dn
    StartIndex += Mup[i].val.N-1;

    oadd[1] = CrAnMask+2*i+1;
    oadd[0] = StartIndex;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn[i].bvB0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask+StartIndex;
    oadd[0] = 2*i+1;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn[i].bvB0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask+2*NSite+2*i+1;
    oadd[0] = StartIndex;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn[i].bvA0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask+StartIndex;
    oadd[0] = 2*NSite+2*i+1;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn[i].bvA0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask+StartIndex;
    oadd[0] = StartIndex;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn[i].val.a[2], &QDetNotOrderedFcFaBcBa);

    for (j=0; j+2<Mdn[i].val.N; j++)
    {
      oadd[1] = CrAnMask+StartIndex+j;
      oadd[0] = StartIndex+j+1;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn[i].val.b[j], &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask+StartIndex+j+1;
      oadd[0] = StartIndex+j;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn[i].val.b[j], &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask+StartIndex+j+1;
      oadd[0] = StartIndex+j+1;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn[i].val.a[j+1], &QDetNotOrderedFcFaBcBa);
    }

    //con up
    StartIndex += Mdn[i].val.N-1;

    oadd[1] = CrAnMask+2*i;
    oadd[0] = StartIndex;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup[i].bcB0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask+StartIndex;
    oadd[0] = 2*i;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup[i].bcB0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask+2*NSite+2*i;
    oadd[0] = StartIndex;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup[i].bcA0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask+StartIndex;
    oadd[0] = 2*NSite+2*i;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup[i].bcA0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask+StartIndex;
    oadd[0] = StartIndex;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup[i].con.a[2], &QDetNotOrderedFcFaBcBa);

    for (j=0; j+2<Mup[i].con.N; j++)
    {
      oadd[1] = CrAnMask+StartIndex+j;
      oadd[0] = StartIndex+j+1;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup[i].con.b[j], &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask+StartIndex+j+1;
      oadd[0] = StartIndex+j;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup[i].con.b[j], &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask+StartIndex+j+1;
      oadd[0] = StartIndex+j+1;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup[i].con.a[j+1], &QDetNotOrderedFcFaBcBa);
    }

    // con dn
    StartIndex += Mup[i].con.N-1;

    oadd[1] = CrAnMask+2*i+1;
    oadd[0] = StartIndex;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn[i].bcB0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask+StartIndex;
    oadd[0] = 2*i+1;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn[i].bcB0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask+2*NSite+2*i+1;
    oadd[0] = StartIndex;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn[i].bcA0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask+StartIndex;
    oadd[0] = 2*NSite+2*i+1;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn[i].bcA0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask+StartIndex;
    oadd[0] = StartIndex;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn[i].con.a[2], &QDetNotOrderedFcFaBcBa);

    for (j=0; j+2<Mdn[i].con.N; j++)
    {
      oadd[1] = CrAnMask+StartIndex+j;
      oadd[0] = StartIndex+j+1;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn[i].con.b[j], &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask+StartIndex+j+1;
      oadd[0] = StartIndex+j;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn[i].con.b[j], &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask+StartIndex+j+1;
      oadd[0] = StartIndex+j+1;
      RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn[i].con.a[j+1], &QDetNotOrderedFcFaBcBa);
    }
  }
  OperatorAdd(O, *Oloc, 1.0, 0.0);

  return 0;
}


uint32_t MultiSitesFindGroundStateFromG0andSigma(MultiSitesDMFTSettings *const O,
                                                 AndersonMatrixType *const G0,
                                                 ListOfPoles *const Sigmaup,
                                                 ListOfPoles *const Sigmadn,
                                                 double *const mu,
                                                 double *const Nelectron,
                                                 WaveFunctionType *const Psi,
                                                 BADoubleTriDiagonalMatrixType *const GBathupBARet,
                                                 BADoubleTriDiagonalMatrixType *const GBathdnBARet,
                                                 AndersonMatrixType *const GBathupAndRet,
                                                 AndersonMatrixType *const GBathdnAndRet,
                                                 TimeType *const Time,
                                                 double *const onsiteshiftstart)
{
  uint32_t i, j, k, StartIndex, converged=0, iiter;
  uint32_t const nitermax=10, nitermem=30;
  AndersonMatrixType GBathupAnd[O->NSite], GBathdnAnd[O->NSite];
  BADoubleTriDiagonalMatrixType GBathupBA[O->NSite], GBathupBATrunc[O->NSite];
  BADoubleTriDiagonalMatrixType GBathdnBA[O->NSite], GBathdnBATrunc[O->NSite];
  double WaveFunctionAccuracy, PsiEnergy, PsiQuality;
  double occupation[O->NSite], occupationup[O->NSite], occupationdn[O->NSite];
  double Rho[4*(O->NSite+O->NBath)];
  double Occu[O->NSite][nitermax+1],  TotalOccu[nitermax+1];
  double ABOccu[O->NSite][nitermax+1], TotalABOccu[nitermax+1];
  double onsiteshift[nitermax+1], mulist[nitermax+1];
  double epsilon_conv=1E-5;
  OperatorType HLoc, HAnd, Sx, Sy, Sz, S2;
  uint16_t *SpinUp,*SpinDn;

  WaveFunctionAccuracy = O->zeroLancWavefunction;
  onsiteshift[0] = *onsiteshiftstart;
  mulist[0] = *mu;
  iiter = 0;

  #pragma omp parallel for
  for (i=0; i<O->NSite; i++)
  {
    GBathAndListAndRepresentation(G0[i], *mu, Sigmaup[i], &GBathupAnd[i]);
    GBathAndListAndRepresentation(G0[i], *mu, Sigmadn[i], &GBathdnAnd[i]);
    if ((O->SymOpt)&2)
    {
      printf("Symmetrize the bath greens function around mu=%22.15E\n", *mu);
      SymmetrizeAndersonMatrix(&GBathupAnd[i]);
      SymmetrizeAndersonMatrix(&GBathdnAnd[i]);
    }
    ReducePolesAndersonMatrixNumberFermi(&GBathupAnd[i], O->NGBath, *mu);
    ReducePolesAndersonMatrixNumberFermi(&GBathdnAnd[i], O->NGBath, *mu);
  }

  while (!converged)
  {
    #pragma omp parallel for
    for (i=0; i<O->NSite; i++)
    {
      printf("\nSpin Up:");
      GBathupAnd[i].mu = mulist[iiter];
      AndersonMatrixToBADoubleTriDiagonalMatrixTemperature(GBathupAnd[i], onsiteshift[iiter], O->Teff, &GBathupBA[i], &occupationup[i]);
      GBathupBATrunc[i] = GBathupBA[i];
      GBathupBATrunc[i].con.N = Min(GBathupBA[i].con.N, O->NBath+1);
      GBathupBATrunc[i].val.N = Min(GBathupBA[i].val.N, O->NBath+1);

      printf("\nSpin Down:");
      GBathdnAnd[i].mu = mulist[iiter];
      AndersonMatrixToBADoubleTriDiagonalMatrixTemperature(GBathdnAnd[i], onsiteshift[iiter], O->Teff, &GBathdnBA[i], &occupationdn[i]);
      GBathdnBATrunc[i] = GBathdnBA[i];
      GBathdnBATrunc[i].con.N = Min(GBathdnBA[i].con.N, O->NBath+1);
      GBathdnBATrunc[i].val.N = Min(GBathdnBA[i].val.N, O->NBath+1);

      occupation[i]=occupationup[i]+occupationdn[i];
      printf("\nTotal Impurity occupation is %22.15E for shift %22.15E"
             " at chemical potential %22.15E\n",
             occupation[i], onsiteshift[iiter], mulist[iiter]);
    }

    MultiSitesTestPsiDimensions(O, GBathupBATrunc, GBathdnBATrunc, Psi);
    //LocalHamiltonian(&HLoc, GBathupBATrunc, GBathdnBATrunc, O->NSite, O->U, O->J, 0);
    DCAHamiltonian(&HLoc, GBathupBATrunc, GBathdnBATrunc, O->NSite, O->U, 0);
    MultiSitesBADoubleTriDiagonalMatrixToAndersonImpurityOperator(&HAnd, &HLoc, GBathupBATrunc, GBathdnBATrunc, O->NSite);

    if ((O->SymOpt)&4)
    {
      printf("Restricted HAnd for the ground-state to the first bath site\n");
      HAnd.l[2].N = 8*O->NSite;
    }

    if ((O->SymOpt)&2)
    {
      printf("Symmetrize the bath greens function around mu=%22.15E\n", *mu);
      #pragma omp parallel for
      for (i=0; i<O->NSite; i++) {
        SymmetrizeBADoubleTriDiagonalMatrix(&GBathupBA[i]);
        SymmetrizeBADoubleTriDiagonalMatrix(&GBathdnBA[i]);
      }
    }

    /* if singlet set is true add the (U+1) * S^2 to the operator HAndGrd */
    if ((O->SymOpt)&1)
    {
      printf("Add %22.15E * S^2 to HAnd for the calculation of the ground-state wave-function\n", 1+O->U);
      printf("This stabilizes the ground-state by enlarging the singlet triplet energy\n");
      SpinUp=calloc(HAnd.NFermion/2, sizeof(uint16_t));
      SpinDn=calloc(HAnd.NFermion/2, sizeof(uint16_t));
      for (i=0; i<2*O->NSite; i++)
      {
        SpinUp[i] = 2*i;
        SpinDn[i] = 2*i+1;
      }

      StartIndex = 4*O->NSite;
      for (i=0; i<O->NSite; i++)
      {
        for (j=0; j+1<GBathupBATrunc[i].val.N; j++)
        {
          SpinUp[StartIndex/2+j] = StartIndex+j;
        }
        for (j=0; j+1<GBathdnBATrunc[i].val.N; j++)
        {
          SpinDn[StartIndex/2+j] = StartIndex+GBathupBATrunc[i].val.N-1+j;
        }
        StartIndex += GBathupBATrunc[i].val.N-1+GBathdnBATrunc[i].val.N-1;

        for (j=0; j+1<GBathupBATrunc[i].con.N; j++)
        {
          SpinUp[StartIndex/2+j] = StartIndex+j;
        }
        for (j=0; j+1<GBathdnBATrunc[i].con.N; j++)
        {
          SpinDn[StartIndex/2+j] = StartIndex+GBathupBATrunc[i].con.N-1+j;
        }
        StartIndex += GBathupBATrunc[i].con.N-1+GBathdnBATrunc[i].con.N-1;
      }

      CreateSpinOperators(HAnd.NFermion, HAnd.NFermion/2, SpinUp, SpinDn, &Sx, &Sy, &Sz, &S2);
      FreeOperator(&Sx);
      FreeOperator(&Sy);
      FreeOperator(&Sz);
      OperatorAddNormalOrder(&HAnd, S2, fabs(O->U)+1, 0, &QDetNotOrderedFcFaBcBa);
      FreeOperator(&S2);
    }

    /* create ground-state wave-function */
    BeginTiming("Lanczos GroundState", Time);
    LanczosGroundState(HAnd, Psi, &WaveFunctionAccuracy, O->zeroLancEnergy, O->U+*mu, Psi, 0);
    EndTiming("Lanczos GroundState", Time);
    /* calculate NelectronNew */
    CalculateDensity(Rho,*Psi);

    TotalOccu[iiter] = 0;
    TotalABOccu[iiter] = 0;
    for (i=0; i<O->NSite; i++)
    {
      Occu[i][iiter] = Rho[2*i]+Rho[2*i+1];
      TotalOccu[iiter] += Occu[i][iiter];
      ABOccu[i][iiter] = Rho[2*i]+Rho[2*i+1]+Rho[2*O->NSite+2*i]+Rho[2*O->NSite+2*i+1];
      TotalABOccu[iiter] += ABOccu[i][iiter];
    }

    /* determine if ground-state calculation has to be restarted */
    if (O->NmuSym&2)
    {
      /* calculations at constant mu, but still need to make sure the
       * bath Greens function is calculated at the same filling as the
       * impurity greens function
       *
       * this is done by optimizing NelectronNew such that the occupation of
       * the impurity and the selected bath site is integer
       */
      if(Abs(TotalABOccu[iiter]-2*O->NSite) > epsilon_conv)
      {
        mulist[iiter+1]=*mu;
        if( iiter == 0)
        {
          onsiteshift[iiter+1] = onsiteshift[iiter]-(TotalABOccu[iiter]-2*O->NSite)*(0.5*O->NSite+Abs(O->U/2));
        }
        else if (iiter == nitermax)
        {
          printf("occupation of bath and impurity green function"
                 "did not converge in %i itterations\n", iiter);
          printf("\noccupation vs. onsite-shift:\n");
          printf("\n==========in total==========\n");
          printf("loop   onsiteshift            total occupation       mu                     impurity occupation\n");
          for (j=0; j<=iiter; j++)
          {
            printf("%3i %22.15E %22.15E %22.15E %22.15E\n",
                   j, onsiteshift[j], TotalABOccu[j], mulist[j], TotalOccu[j]);
          }
          for (k=0; k<O->NSite; k++)
          {
            printf("\n========== site %i ==========\n", k);
            printf("loop   onsiteshift            total occupation       mu                     impurity occupation\n");
            for(j=0; j<=iiter; j++)
            {
              printf("%3i %22.15E %22.15E %22.15E %22.15E\n",
                     j, onsiteshift[j], ABOccu[k][j], mulist[j], Occu[k][j]);
            }
          }
          converged=1;
          *onsiteshiftstart=(Abs(onsiteshift[iiter]) < *mu) ? onsiteshift[iiter] : *mu;
        }
        else
        {
          onsiteshift[iiter+1] = InterpolateRootOfset(&(onsiteshift[Max(iiter, nitermem)-nitermem]),
                                                      &(TotalABOccu[Max(iiter, nitermem)-nitermem]),
                                                      Min(iiter+1, nitermem), 2*O->NSite);
        }
        iiter++;
      }
      else
      {
        converged=1;
        *onsiteshiftstart=(Abs(onsiteshift[iiter]) < *mu) ? onsiteshift[iiter] : *mu;
        printf("\n==========in total==========\n");
        printf("loop   onsiteshift            total occupation       mu                     impurity occupation\n");
        for (j=0; j<=iiter; j++)
        {
          printf("%3i %22.15E %22.15E %22.15E %22.15E\n",
                 j, onsiteshift[j], TotalABOccu[j], mulist[j], TotalOccu[j]);
        }
        for (k=0; k<O->NSite; k++)
        {
          printf("\n========== site %i ==========\n", k);
          printf("loop   onsiteshift            total occupation       mu                     impurity occupation\n");
          for(j=0; j<=iiter; j++)
          {
            printf("%3i %22.15E %22.15E %22.15E %22.15E\n",
                   j, onsiteshift[j], ABOccu[k][j], mulist[j], Occu[k][j]);
          }
        }
      }
    }
    else
    {
      /* calculations at constant occupation, find the correct mu and onsiteshift */
      if((Abs(TotalABOccu[iiter]-2*O->NSite) > epsilon_conv) ||
         (Abs(TotalOccu[iiter]-*Nelectron) > epsilon_conv))
      {
        if( iiter == 0)
        {
          onsiteshift[iiter+1] = onsiteshift[iiter]-(TotalABOccu[iiter]-2*O->NSite)*(0.5*O->NSite+Abs(O->U/2));
          mulist[iiter+1]=mulist[iiter];
        }
        else if (iiter == nitermax)
        {
          printf("occupation of bath and impurity green function"
                 "did not converge in %i itterations\n", iiter);
          printf("loop   onsiteshift            total occupation       mu                     impurity occupation\n");
          for(j=0; j<=iiter; j++)
          {
            printf("%3i %22.15E %22.15E %22.15E %22.15E\n",
                   j, onsiteshift[j], TotalABOccu[j], mulist[j], TotalOccu[j]);
          }
          converged=1;
          *mu=mulist[nitermax];
        }
        else
        {
          InterpolateRoot2DOfset(&(onsiteshift[Max(iiter, nitermem)-nitermem]),
                                 &(mulist[Max(iiter, nitermem)-nitermem]),
                                 &(TotalABOccu[Max(iiter, nitermem)-nitermem]),
                                 &(TotalOccu[Max(iiter, nitermem)-nitermem]),
                                 Min(iiter+1, nitermem), 2*O->NSite, *Nelectron,
                                 &(onsiteshift[iiter+1]), &(mulist[iiter+1]));
        }
        iiter++;
      }
      else
      {
        converged=1;
        *mu=mulist[iiter];
      }
    }
    if (!converged)
      FreeOperator(&HAnd);
  }
  /* print info on groundstate */
  for (i=0; i<O->NSite; i++)
  {
    printf("\nSpin Up:");
    PrintInfoOnGBathAndG0InDMFTFindGroundstate(G0[i], GBathupBA[i]);
    printf("\nSpin Down:");
    PrintInfoOnGBathAndG0InDMFTFindGroundstate(G0[i], GBathdnBA[i]);
  }
  PsiHermitianOperatorPsi(HAnd, *Psi, &PsiEnergy);
  PsiHermitianOperatorSqrPsi(HAnd, *Psi, &PsiQuality);
  PsiQuality/=Sqr(PsiEnergy);
  printf("Energy                  <psi|H^2|psi>/(<psi|H|psi>^2)\n");
  printf("%22.15E  %22.15E\n", PsiEnergy, PsiQuality);
  PrintWaveFunctionDensity(*Psi);
  PrintWaveFunctionPrefactorHistogram(*Psi, 2.0, 20);
  for (i=0; i<O->NSite; i++)
  {
    GBathupBARet[i] = GBathupBA[i];
    GBathupAndRet[i] = GBathupAnd[i];
    GBathdnBARet[i] = GBathdnBA[i];
    GBathdnAndRet[i] = GBathdnAnd[i];
  }

  return 0;
}


/*  Input:
 *    DMFT settings
 *    GBathBA as a double list representing the conduction and valence bath,
 *    inclusive a form where one orbital is taken to be bonding/anitbonding with
 *    the impurity
 *  Output
 *    GBathBATrunc as a copy of GBathBA with the chains the length
 *    included in the groundstate cluster
 */
uint32_t MultiSitesCreateGBathBATrunc(MultiSitesDMFTSettings *const O,
                                      BADoubleTriDiagonalMatrixType *const GBathBA,
                                      BADoubleTriDiagonalMatrixType *const GBathBATrunk)
{
  uint32_t i;

  for (i=0; i<O->NSite; i++)
  {
    GBathBATrunk[i] = GBathBA[i];

    if(GBathBA[i].con.N < O->NBath+1)
    {
      printf("Warning: GBathBA.con.N (%i) < DMFTSettings.NBath + 1 (%i+1)\n",
             GBathBA[i].con.N, O->NBath);
    }
    else
    {
      GBathBATrunk[i].con.N = O->NBath+1;
    }
    if(GBathBA[i].val.N < O->NBath+1)
    {
      printf("Warning: GBathBA.val.N (%i) < DMFTSettings.NBath + 1 (%i+1)\n",
             GBathBA[i].val.N, O->NBath);
    }
    else
    {
      GBathBATrunk[i].val.N = O->NBath+1;
    }
    sprintf(GBathBATrunk[i].Name, "%i - GBath for ground-state", i);
  }

  return 0;
}


/*  Input:
 *    DMFT settings
 *    GBathBA as a double list representing the conduction and valence bath,
 *    inclusive a form where one orbital is taken to be bonding/anitbonding with
 *    the impurity
 *  Output
 *    VPerturbation - operator describing the coupling of the bath included for
 *                    the calculation of the ground-state to the rest of the world
 *    gPerturbation - that part of the bath Greens function not included in
 *                    GBathBATrunc
 */
uint32_t MultiSitesCreateVPerturbation(MultiSitesDMFTSettings *const O,
                                       BADoubleTriDiagonalMatrixType *const GBathupBA,
                                       BADoubleTriDiagonalMatrixType *const GBathdnBA,
                                       OperatorType *const VPerturbation,
                                       TriDiagonalMatrixType *const gPerturbation)
{
  uint32_t i, j, NFermion, StartIndex;
  uint16_t oadd[1];

  NFermion = 4*O->NSite;
  for (i=0; i<O->NSite; i++)
  {
    NFermion += Min(GBathupBA[i].val.N-1, O->NBath)+Min(GBathupBA[i].con.N-1, O->NBath)+
                Min(GBathdnBA[i].val.N-1, O->NBath)+Min(GBathdnBA[i].con.N-1, O->NBath);
  }
  for (i=0; i<4*O->NSite; i++)
  {
    InitOperator(&(VPerturbation[i]));
    VPerturbation[i].NFermion = NFermion;
    VPerturbation[i].NBoson = 0;
    VPerturbation[i].NWord = Min(VPerturbation[i].NFermion, 1+(VPerturbation[i].NFermion-1)/8);
    VPerturbation[i].QComplex = 0;
    InitOperatorFixedLength(&(VPerturbation[i]), 1, 1);
  }

  /* Create VPerturbation */
  #pragma omp parallel for private(StartIndex, j)
  for (i=0; i<O->NSite; i++)
  {
    StartIndex = 4*O->NSite;
    for (j=0; j<i; j++)
    {
      StartIndex += Min(GBathupBA[j].val.N-1, O->NBath)+Min(GBathupBA[j].con.N-1, O->NBath)+
                    Min(GBathdnBA[j].val.N-1, O->NBath)+Min(GBathdnBA[j].con.N-1, O->NBath);
    }

    if (GBathupBA[i].val.N-1 > O->NBath)
    {
      oadd[0] = CrAnMask+StartIndex+O->NBath-1; //spin up
      RealOperatorAddLader(&(VPerturbation[4*i  ]), oadd, 1, GBathupBA[i].val.b[O->NBath+1]);
    }

    StartIndex += Min(GBathupBA[i].val.N-1, O->NBath);
    if (GBathdnBA[i].val.N-1 > O->NBath)
    {
      oadd[0] = CrAnMask+StartIndex+O->NBath-1; //spin dn
      RealOperatorAddLader(&(VPerturbation[4*i+1]), oadd, 1, GBathdnBA[i].val.b[O->NBath+1]);
    }

    StartIndex += Min(GBathdnBA[i].val.N-1, O->NBath);
    if (GBathupBA[i].con.N-1 > O->NBath)
    {
      oadd[0] = StartIndex+O->NBath-1; //spin up
      RealOperatorAddLader(&(VPerturbation[4*i+2]), oadd, 1, GBathupBA[i].con.b[O->NBath+1]);
    }

    StartIndex += Min(GBathupBA[i].con.N-1, O->NBath);
    if (GBathdnBA[i].con.N-1 > O->NBath)
    {
      oadd[0] = StartIndex+O->NBath-1; //spin dn
      RealOperatorAddLader(&(VPerturbation[4*i+3]), oadd, 1, GBathdnBA[i].con.b[O->NBath+1]);
    }
  }

  /* Create gPerturbation */
  #pragma omp parallel for private(j)
  for (i=0; i<O->NSite; i++)
  {
    if(GBathupBA[i].val.N-1 > O->NBath)
    {
      gPerturbation[4*i  ].N = GBathupBA[i].val.N-O->NBath-1;
      InitTriDiagonalMatrix(&(gPerturbation[4*i  ]));
      for (j=1; j<GBathupBA[i].val.N-O->NBath-1; j++)
      {
        gPerturbation[4*i  ].a[j] =-GBathupBA[i].val.a[j+O->NBath+1];
        gPerturbation[4*i  ].b[j] = GBathupBA[i].val.b[j+O->NBath+1];
      }
      gPerturbation[4*i  ].a[GBathupBA[i].val.N-O->NBath-1] = -GBathupBA[i].val.a[GBathupBA[i].val.N];
    }
    else
    {
      gPerturbation[4*i  ].N = 1;
      InitTriDiagonalMatrix(&(gPerturbation[4*i  ]));
      gPerturbation[4*i  ].a[1] = 0;
    }

    if(GBathdnBA[i].val.N-1 > O->NBath)
    {
      gPerturbation[4*i+1].N = GBathdnBA[i].val.N-O->NBath-1;
      InitTriDiagonalMatrix(&(gPerturbation[4*i+1]));
      for (j=1; j<GBathdnBA[i].val.N-O->NBath-1; j++)
      {
        gPerturbation[4*i+1].a[j] =-GBathdnBA[i].val.a[j+O->NBath+1];
        gPerturbation[4*i+1].b[j] = GBathdnBA[i].val.b[j+O->NBath+1];
      }
      gPerturbation[4*i+1].a[GBathdnBA[i].val.N-O->NBath-1] = -GBathdnBA[i].val.a[GBathdnBA[i].val.N];
    }
    else
    {
      gPerturbation[4*i+1].N = 1;
      InitTriDiagonalMatrix(&(gPerturbation[4*i+1]));
      gPerturbation[4*i+1].a[1] = 0;
    }

    if(GBathupBA[i].con.N-1 > O->NBath)
    {
      gPerturbation[4*i+2].N = GBathupBA[i].con.N-O->NBath-1;
      InitTriDiagonalMatrix(&(gPerturbation[4*i+2]));
      for (j=1; j<GBathupBA[i].con.N-O->NBath-1; j++)
      {
        gPerturbation[4*i+2].a[j] = GBathupBA[i].con.a[j+O->NBath+1];
        gPerturbation[4*i+2].b[j] = GBathupBA[i].con.b[j+O->NBath+1];
      }
      gPerturbation[4*i+2].a[GBathupBA[i].con.N-O->NBath-1] = GBathupBA[i].con.a[GBathupBA[i].con.N];
    }
    else
    {
      gPerturbation[4*i+2].N=1;
      InitTriDiagonalMatrix(&(gPerturbation[4*i+2]));
      gPerturbation[4*i+2].a[1]=0;
    }

    if(GBathdnBA[i].con.N-1 > O->NBath)
    {
      gPerturbation[4*i+3].N = GBathdnBA[i].con.N-O->NBath-1;
      InitTriDiagonalMatrix(&(gPerturbation[4*i+3]));
      for (j=1; j<GBathdnBA[i].con.N-O->NBath-1; j++)
      {
        gPerturbation[4*i+3].a[j] = GBathdnBA[i].con.a[j+O->NBath+1];
        gPerturbation[4*i+3].b[j] = GBathdnBA[i].con.b[j+O->NBath+1];
      }
      gPerturbation[4*i+3].a[GBathdnBA[i].con.N-O->NBath-1] = GBathdnBA[i].con.a[GBathdnBA[i].con.N];
    }
    else
    {
      gPerturbation[4*i+3].N=1;
      InitTriDiagonalMatrix(&(gPerturbation[4*i+3]));
      gPerturbation[4*i+3].a[1]=0;
    }
  }

  return 0;
}


uint32_t MultiSitesCreateGBathBATruncandVPerturbation(MultiSitesDMFTSettings *const O,
                                                      BADoubleTriDiagonalMatrixType *const GBathupBA,
                                                      BADoubleTriDiagonalMatrixType *const GBathdnBA,
                                                      BADoubleTriDiagonalMatrixType *const GBathupBATrunc,
                                                      BADoubleTriDiagonalMatrixType *const GBathdnBATrunc,
                                                      OperatorType *const VPerturbation,
                                                      TriDiagonalMatrixType *const gPerturbation)
{
  MultiSitesCreateGBathBATrunc(O, GBathupBA, GBathupBATrunc);
  MultiSitesCreateGBathBATrunc(O, GBathdnBA, GBathdnBATrunc);
  MultiSitesCreateVPerturbation(O, GBathupBA, GBathdnBA, VPerturbation, gPerturbation);
  return 0;
}


uint32_t CreatePESOperator(OperatorType *O, uint32_t NFermion, uint32_t i)
{
  uint16_t oadd;
  InitOperator(O);
  O->NFermion=NFermion;
  O->NWord=Min(O->NFermion, 1+(O->NFermion-1)/8);
  InitOperatorFixedLength(O, 1, 1);
  oadd=i;
  RealOperatorAddLader(O, &oadd, 1, 1.0);
  return 0;
}

uint32_t CreateIPESOperator(OperatorType *O, uint32_t NFermion, uint32_t i)
{
  uint16_t oadd;
  InitOperator(O);
  O->NFermion=NFermion;
  O->NWord=Min(O->NFermion, 1+(O->NFermion-1)/8);
  InitOperatorFixedLength(O, 1, 1);
  oadd=CrAnMask+i;
  RealOperatorAddLader(O, &oadd, 1, 1.0);
  return 0;
}


uint32_t CreateMultiSitesGclDMFT(MultiSitesDMFTSettings *const O,
                                 OperatorType const HAnd,
                                 WaveFunctionType const Psi,
                                 OperatorType const *VPerturbation,
                                 TriDiagonalMatrixType const *gPerturbation,
                                 double const mu,
                                 AndersonMatrixType *const GclupAnd,
                                 AndersonMatrixType *const GcldnAnd,
                                 TimeType *const Time)
{
  uint32_t i, j;
  double Egrd;
  OperatorType OPESup, OPESdn, OIPESup, OIPESdn;
  WaveFunctionType OPsi;
  TriDiagonalMatrixType GclMup, GclMdn, GclPup, GclPdn;
  ListOfPoles GclMupLP, GclMdnLP, GclPupLP, GclPdnLP, GclupLP, GcldnLP, GclavLP;
  TimeType subTime[O->NSite];

  for (i=0; i<O->NSite; i++)
  {
    InitializeTiming(&subTime[i]);
  }

  /* Get ground-state energy needed for shift */
  PsiHermitianOperatorPsi(HAnd, Psi, &Egrd);

  #pragma omp parallel for private(OPESup, OIPESup, OPESdn, OIPESdn, OPsi, \
                                   GclMup, GclMdn, GclPup, GclPdn,\
                                   GclMupLP, GclMdnLP, GclPupLP, \
                                   GclPdnLP, GclupLP, GcldnLP, GclavLP, j)
  for (i=0; i<O->NSite; i++)
  {
    CreatePESOperator(&OPESup, HAnd.NFermion, 2*i);
    CreateIPESOperator(&OIPESup, HAnd.NFermion, 2*i);
    CreatePESOperator(&OPESdn, HAnd.NFermion, 2*i+1);
    CreateIPESOperator(&OIPESdn, HAnd.NFermion, 2*i+1);

    BeginTiming("Lanczos Tri Diagonalization", &subTime[i]);
    OperatorPsi(OPESup, Psi, &OPsi);
    switch (O->Dyson)
    {
      case 1:
        GclMup.N = Max(O->NTriDiag, O->NGCluster);
        InitTriDiagonalMatrix(&GclMup);
        LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, O->zeroLancWavefunction,
                                         O->NTriDiag, VPerturbation,
                                         gPerturbation, 4*O->NSite, &GclMup);
        break;
      default:
        GclMup.N = O->NTriDiag;
        InitTriDiagonalMatrix(&GclMup);
        LanczosTriDiagonalize(HAnd, OPsi, O->zeroLancWavefunction,
                              &GclMup, O->Gamma);
        break;
    }
    EndTiming("Lanczos Tri Diagonalization", &subTime[i]);

    BeginTiming("Lanczos Tri Diagonalization", &subTime[i]);
    OperatorPsi(OPESdn, Psi, &OPsi);
    switch (O->Dyson)
    {
      case 1:
        GclMdn.N = Max(O->NTriDiag, O->NGCluster);
        InitTriDiagonalMatrix(&GclMdn);
        LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, O->zeroLancWavefunction,
                                         O->NTriDiag, VPerturbation,
                                         gPerturbation, 4*O->NSite, &GclMdn);
        break;
      default:
        GclMdn.N = O->NTriDiag;
        InitTriDiagonalMatrix(&GclMdn);
        LanczosTriDiagonalize(HAnd, OPsi, O->zeroLancWavefunction,
                              &GclMdn, O->Gamma);
        break;
    }
    EndTiming("Lanczos Tri Diagonalization", &subTime[i]);

    BeginTiming("Lanczos Tri Diagonalization", &subTime[i]);
    OperatorPsi(OIPESup,Psi,&OPsi);
    switch (O->Dyson)
    {
      case 1:
        GclPup.N = Max(O->NTriDiag, O->NGCluster);
        InitTriDiagonalMatrix(&GclPup);
        LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, O->zeroLancWavefunction,
                                         O->NTriDiag, VPerturbation,
                                         gPerturbation, 4*O->NSite, &GclPup);
        break;
      default:
        GclPup.N = O->NTriDiag;
        InitTriDiagonalMatrix(&GclPup);
        LanczosTriDiagonalize(HAnd, OPsi, O->zeroLancWavefunction,
                              &GclPup, O->Gamma);
        break;
    }
    EndTiming("Lanczos Tri Diagonalization", &subTime[i]);

    BeginTiming("Lanczos Tri Diagonalization", &subTime[i]);
    OperatorPsi(OIPESdn,Psi,&OPsi);
    switch (O->Dyson)
    {
      case 1:
        GclPdn.N = Max(O->NTriDiag, O->NGCluster);
        InitTriDiagonalMatrix(&GclPdn);
        LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, O->zeroLancWavefunction,
                                         O->NTriDiag, VPerturbation,
                                         gPerturbation, 4*O->NSite, &GclPdn);
        break;
      default:
        GclPdn.N = O->NTriDiag;
        InitTriDiagonalMatrix(&GclPdn);
        LanczosTriDiagonalize(HAnd, OPsi, O->zeroLancWavefunction,
                              &GclPdn, O->Gamma);
        break;
    }
    EndTiming("Lanczos Tri Diagonalization", &subTime[i]);

    ShiftTriDiagonalMatrix(&GclMup, -Egrd);
    ShiftTriDiagonalMatrix(&GclMdn, -Egrd);
    ShiftTriDiagonalMatrix(&GclPup, -Egrd);
    ShiftTriDiagonalMatrix(&GclPdn, -Egrd);

    GclMup.mu = mu;
    GclMdn.mu = mu;
    GclPup.mu = mu;
    GclPdn.mu = mu;

    for(j=1; j<=GclMup.N; j++)
    {
      GclMup.a[j] *= -1;
    }
    for(j=1; j<=GclMdn.N; j++)
    {
      GclMdn.a[j] *= -1;
    }

    sprintf(GclMup.Name, "%i-Gcluster-Min-Up", i);
    sprintf(GclMdn.Name, "%i-Gcluster-Min-Dn", i);
    sprintf(GclPup.Name, "%i-Gcluster-Plus-Up", i);
    sprintf(GclPdn.Name, "%i-Gcluster-Plus-Dn", i);

    printf("\nStart SumGreensfunctionsTriDiagonalRepresentation"
           "for site %i\n", i); fflush(stdout);

    TriDiagonalMatrixToListOfPoles(GclMup, &GclMupLP);
    TriDiagonalMatrixToListOfPoles(GclMdn, &GclMdnLP);
    TriDiagonalMatrixToListOfPoles(GclPup, &GclPupLP);
    TriDiagonalMatrixToListOfPoles(GclPdn, &GclPdnLP);

    if (O->NmuSym&2)
    {
      /* calculations at constant mu, move spectra at negative energy to zero in order to excellerate converges */
      ListOfPolesMoveNegativeEnergyToMu(&GclPupLP, mu);
      ListOfPolesMoveNegativeEnergyToMu(&GclPdnLP, mu);
      ListOfPolesMovePositiveEnergyToMu(&GclMupLP, mu);
      ListOfPolesMovePositiveEnergyToMu(&GclMdnLP, mu);
    }

    ListOfPolesCopy(&GclupLP, GclMupLP);
    AddListOfPoles(&GclupLP, GclPupLP);
    ListOfPolesCopy(&GcldnLP, GclMdnLP);
    AddListOfPoles(&GcldnLP, GclPdnLP);

    sprintf(GclupLP.Name, "%i-Gcluster-up", i);
    sprintf(GcldnLP.Name, "%i-Gcluster-dn", i);

    if (O->SymOpt&1)
    {
      printf("Enforcing a paramagnetic (singlet) calculation!\n");
      ListOfPolesCopy(&GclavLP, GclupLP);
      AddListOfPoles(&GclavLP, GcldnLP);
      ScaleListOfPoles(&GclavLP,0.5);
      GclavLP.mu = mu;

      ListOfPolesCopy(&GclupLP, GclavLP);
      ListOfPolesCopy(&GcldnLP, GclavLP);
    }

    if (O->SymOpt&2)
    {
      printf("Symmetryze %s (%s) around mu = %22.15E (%22.15E)\n",
             GclupLP.Name, GcldnLP.Name, GclupLP.mu, GcldnLP.mu);
      SymmetrizeListOfPoles(&GclupLP);
      SymmetrizeListOfPoles(&GcldnLP);
    }

    printf("Cluster Green function created:\n");
    printf("spin up for site %i of size %i ", i, GclupLP.N);
    ReduceListOfPolesNumberFermi(&GclupLP, O->NGCluster, mu);
    printf("reduced to %i\n", GclupLP.N);
    printf("spin dn for site %i of size %i ", i, GcldnLP.N);
    ReduceListOfPolesNumberFermi(&GcldnLP, O->NGCluster, mu);
    printf("reduced to %i\n", GcldnLP.N);
    GclupLP.mu = mu;
    GcldnLP.mu = mu;
    ListOfPolesToAndersonMatrix(GclupLP, &GclupAnd[i]);
    ListOfPolesToAndersonMatrix(GcldnLP, &GcldnAnd[i]);
    sprintf(GclupAnd[i].Name, "%i-Gcluster-up", i);
    sprintf(GcldnAnd[i].Name, "%i-Gcluster-dn", i);
    GclupAnd[i].b[0] = 1.0;
    GcldnAnd[i].b[0] = 1.0;

    /* print occupations from spectra */
    PrintListOfPolesEssentials(GclMdnLP);
    PrintListOfPolesEssentials(GclPdnLP);
    PrintListOfPolesEssentials(GclMupLP);
    PrintListOfPolesEssentials(GclPupLP);
    PrintListOfPolesEssentials(GclupLP);
    PrintListOfPolesEssentials(GcldnLP);
    PrintAndersonMatrixEssentials(GclupAnd[i]);
    PrintAndersonMatrixEssentials(GcldnAnd[i]);

    /* free operators and spectra */
    FreeListOfPoles(&GclupLP);
    FreeListOfPoles(&GcldnLP);
    FreeListOfPoles(&GclMupLP);
    FreeListOfPoles(&GclMdnLP);
    FreeListOfPoles(&GclPupLP);
    FreeListOfPoles(&GclPdnLP);
    FreeTriDiagonalMatrix(&GclMup);
    FreeTriDiagonalMatrix(&GclMdn);
    FreeTriDiagonalMatrix(&GclPup);
    FreeTriDiagonalMatrix(&GclPdn);
    FreeOperator(&OPESup);
    FreeOperator(&OPESdn);
    FreeOperator(&OIPESup);
    FreeOperator(&OIPESdn);
  }

  BeginTiming("Lanczos Tri Diagonalization", Time);
  EndTiming("Lanczos Tri Diagonalization", Time);
  for (i = 0; i < Time->N; i++)
  {
    if ((strcmp(Time->Name[i], "Lanczos Tri Diagonalization") == 0) &&
        (strcmp("Lanczos Tri Diagonalization", Time->Name[i]) == 0))
    {
      for (j = 0; j < O->NSite; j++)
      {
        Time->TotalTime[i] += subTime[j].TotalTime[0];
      }
    }
  }

  return 0;
}
