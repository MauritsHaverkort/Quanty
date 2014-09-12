#define NSiteMax 2
typedef struct CalcData3d
{
  /*Atomic data*/
  uint32_t atomicZ;
  uint32_t num3d;
  
  double iF2ddScale,iF4ddScale,fF2ddScale,fF4ddScale,iG1pdScale,iG3pdScale,fG1pdScale,fG3pdScale,iF2pdScale,fF2pdScale;
  double zeta2pScale,zeta3dScale;

  /*Crystal field data*/
  double iCFb1g,iCFa1g,iCFb2g,iCFeg;	/*crystal field energies*/
  double fCFb1g,fCFa1g,fCFb2g,fCFeg;	/*crystal field energies*/
  
  /*Ligand-hybridization data*/
  double Delta,Udd,Upd;
  double iHyb_b1g,iHyb_a1g,iHyb_b2g,iHyb_eg;
  double fHyb_b1g,fHyb_a1g,fHyb_b2g,fHyb_eg;

  /*Ligand band data*/
  /*Either give keyword for a stored type (circle, rect) or give filenames*/
  /*Filenames should just have Anderson-style onsite and hoppings. The hoppings
	are only used for shapes; the total hopping magnitude is scaled to the parameters above.
	Note if just a shape is given, the width and filling is needed.*/

  uint32_t NSite;                   /* Number of cluster sites and or times spin orbitals */
  
  char   BathFileName[NSiteMax][256];/* filename of Bath */
  uint32_t NBath;                   /* Number of Bath states to include */
  uint32_t NStates;                 /* Number of poles of the bath */
  
  
  /*Spectral Info*/
  int XAS, cPES, PES, RIXS;		/*1 for calculation, 0 for not*/
  double energyWindow;	/*States within energyWindow of ground state are computed. Then individual
						 spectra can be combined according to boltzmann population*/

  
  char   SpecFileName[256];         /* filename of spectral functions in spectral representation */
  char   WaveFunctionFileName[256]; /* filename of wave-function  */
  double Nelectron;                 /* Filling of the impurity */
  double mu;                        /* Chemical potential of the lattice Greensfunction */
  double Teff;                      /* effective temperature when determening the relation between chemical potential and filling for the bath */
  uint32_t Dyson;                   /* 0 no Dyson like perturbation - 1 full perturbation - 2 ... */

  uint32_t NE;                      /* number of points in the spectrum */
  double Emin;                      /* minimum energy */
  double Emax;                      /* maximum energy */
  double Gamma;                     /* Broadening only important for the Spec file */
  
  /*Computational info*/
  double zeroLancEnergy;            /* minimum energy difference for convergence of total energy to be different from zero */
  double zeroLancWavefunction;      /* minnimum pre-factor in wave-function to be different from zero */
  uint32_t NTriDiag;                /* number of tri-diagonalization steps in Lanczos */
  uint32_t NBitsKey;                /* 2^NBitsKey == lenght of hash key */
}
CalcData3d;



typedef struct SpectraSettings
{
  uint32_t NSite;                   /* Number of cluster sites and or times spin orbitals */
  char   BathFileName[NSiteMax][256];/* filename of Bath */
  char   SpecFileName[256];         /* filename of spectral functions in spectral representation */
  char   WaveFunctionFileName[256]; /* filename of wave-function  */
  double Nelectron;                 /* Filling of the impurity */
  double mu;                        /* Chemical potential of the lattice Greensfunction */
  double Teff;                      /* effective temperature when determening the relation between chemical potential and filling for the bath */
  uint32_t Dyson;                   /* 0 no Dyson like perturbation - 1 full perturbation - 2 ... */
  uint32_t NBath;                   /* Number of Bath states to include */
  uint32_t NStates;                 /* Number of poles of the bath */
  double Udd;                       /* valence - valence Coulomb repulsion */
  double Upd;
  double Delta;
  double zeta;
  uint32_t Spectra;                 /* type of spectra */
  uint32_t NE;                      /* number of points in the spectrum */
  double Emin;                      /* minimum energy */
  double Emax;                      /* maximum energy */
  double Gamma;                     /* Broadening only important for the Spec file */
  double zeroLancEnergy;            /* minimum energy difference for convergence of total energy to be different from zero */
  double zeroLancWavefunction;      /* minnimum pre-factor in wave-function to be different from zero */
  uint32_t NTriDiag;                /* number of tri-diagonalization steps in Lanczos */
  uint32_t NBitsKey;                /* 2^NBitsKey == lenght of hash key */
}
SpectraSettings;


uint32_t StreamReadSpectraSettings(SpectraSettings * O, FILE * stream)
{
  char   comments[1000];
  uint32_t i;
  
  fscanf(stream,"%i",&(O->NSite));                   fgets(comments,1000,stream);
  for (i=0; i<O->NSite; i++)
  {
    fscanf(stream,"%s",(O->BathFileName[i]));          fgets(comments,1000,stream);
  }
  fscanf(stream, "%s", (O->SpecFileName));
  fgets(comments, 1000, stream);
  fscanf(stream, "%s", (O->WaveFunctionFileName));
  fgets(comments, 1000, stream);
  fscanf(stream, "%lf", &(O->Nelectron));
  fgets(comments, 1000, stream);
  fscanf(stream, "%lf", &(O->mu));
  fgets(comments, 1000, stream);
  fscanf(stream, "%lf", &(O->Teff));
  fgets(comments, 1000, stream);
  fscanf(stream, "%d", &(O->Dyson));
  fgets(comments, 1000, stream);
  fscanf(stream, "%d", &(O->NBath));
  fgets(comments, 1000, stream);
  fscanf(stream, "%d", &(O->NStates));
  fgets(comments, 1000, stream);
  fscanf(stream, "%lf", &(O->Udd));
  fgets(comments, 1000, stream);
  fscanf(stream, "%lf", &(O->Upd));
  fgets(comments, 1000, stream);
  fscanf(stream, "%lf", &(O->Delta));
  fgets(comments, 1000, stream);
  //  fscanf(stream, "%lf", &(O->zeta));
  //  fgets(comments, 1000, stream);
  fscanf(stream, "%d", &(O->Spectra));
  fgets(comments, 1000, stream);
  fscanf(stream, "%d", &(O->NE));
  fgets(comments, 1000, stream);
  fscanf(stream, "%lf", &(O->Emin));
  fgets(comments, 1000, stream);
  fscanf(stream, "%lf", &(O->Emax));
  fgets(comments, 1000, stream);
  fscanf(stream, "%lf", &(O->Gamma));
  fgets(comments, 1000, stream);
  fscanf(stream, "%lf", &(O->zeroLancEnergy));
  fgets(comments, 1000, stream);
  fscanf(stream, "%lf", &(O->zeroLancWavefunction));
  fgets(comments, 1000, stream);
  fscanf(stream, "%d", &(O->NTriDiag));
  fgets(comments, 1000, stream);
  fscanf(stream, "%d", &(O->NBitsKey));
  fgets(comments, 1000, stream);
  
  return 0;
}

uint32_t FileReadSpectraSettings(SpectraSettings * O, const char *const flname)
{
  FILE  *stream;
  uint32_t retval;
  
  if ((stream = fopen(flname, "r")) == NULL)
  {
    printf("ERROR\nCan't open %s for reading.\n", flname);
    return 1;
  }
  retval = StreamReadSpectraSettings(O, stream);
  fclose(stream);
  return retval;
}

uint32_t StreamPrintSpectraSettings(SpectraSettings O, FILE * stream)
{
  uint32_t i;
  
  fprintf(stream, "\nSettigns of AIM-XAS calculation\n");
  fprintf(stream, "===============================\n");
  fprintf(stream, "number of impurities     :  %i\n", O.NSite);
  for (i=0; i<O.NSite; i++)
  {
    fprintf(stream, "BathFileName[%i]         : %s\n", i, O.BathFileName[i]);
  }
  fprintf(stream, "SpecFileName             :  %s\n", O.SpecFileName);
  fprintf(stream, "WaveFunctionFileName     :  %s\n", O.WaveFunctionFileName);
  fprintf(stream, "number of electrons      : %22.15E\n", O.Nelectron);
  fprintf(stream, "Chemical potential       : %22.15E\n", O.mu);
  fprintf(stream, "Effective temperature    : %22.15E\n", O.Teff);
  fprintf(stream,
          "                         converged results do not depent on Teff, convergence speed and stability do\n");
  fprintf(stream, "Dyson perturbation       :  %i\n", O.Dyson);
  fprintf(stream, "NBath                    :  %i\n", O.NBath);
  fprintf(stream, "NStates                  :  %i\n", O.NStates);
  fprintf(stream, "Udd                      : %24.17E\n", O.Udd);
  fprintf(stream, "Upd                      : %24.17E\n", O.Upd);
  fprintf(stream, "Delta                    : %24.17E\n", O.Delta);
  //  fprintf(stream, "zeta                     : %24.17E\n", O.zeta);
  fprintf(stream, "NE                       :  %i\n", O.NE);
  fprintf(stream, "Spectra                  :  %i\n", O.Spectra);
  fprintf(stream, "Emin                     :  %24.17E\n", O.Emin);
  fprintf(stream, "Emax                     : %24.17E\n", O.Emax);
  fprintf(stream, "Gamma                    : %24.17E\n", O.Gamma);
  fprintf(stream, "zeroLancEnergy           : %24.17E\n", O.zeroLancEnergy);
  fprintf(stream, "zeroLancWaveFunction     : %24.17E\n", O.zeroLancWavefunction);
  fprintf(stream, "NTriDiag                 :  %i\n", O.NTriDiag);
  fprintf(stream, "NBitsKey                 :  %i\n", O.NBitsKey);
  return 0;
}

uint32_t PrintSpectraSettings(SpectraSettings O)
{
  return StreamPrintSpectraSettings(O, stdout);
}

/*******************************
      Some functions
 *******************************/

uint32_t DetermineImpurityShiftFromOccupationAndersonMatrix2(AndersonMatrixType M, double const occupation, double *const shift)
{
  uint32_t MaxIt=30, History=10;
  uint32_t i, imin=0;
  double Shifti[MaxIt+1], occupationi[MaxIt+1], totalb=0;
  double epsilonoccu=1.e-10;
  Shifti[0]=*shift;
  M.a[1]+=Shifti[0];
  CalculateImpurityOccupationAndersonMatrix(M, &occupationi[0]);
  M.a[1]-=Shifti[0];
  for(i=0;i<M.N;i++)
  {
    totalb+=Sqr(M.b[i]);
  }
  totalb=sqrt(totalb);
  i=0;
  while((i<MaxIt)&&(Abs(occupationi[i]-occupation)>epsilonoccu))
  {
    i++;
    if(i==1)
    {
      if(occupationi[0]>occupation)
      {
        Shifti[i]=Shifti[i-1]+totalb;
      }
      else
      {
        Shifti[i]=Shifti[i-1]-totalb;
      }
    }
    else
    {
      Shifti[i]=InterpolateRootOfset(&Shifti[Max(History,i)-History], &occupationi[Max(History,i)-History], Min(i-1,History), occupation);
    }
    M.a[1]+=Shifti[i];
    CalculateImpurityOccupationAndersonMatrix(M, &occupationi[i]);
    M.a[1]-=Shifti[i];
    if (Abs(occupationi[i]-occupation)<Abs(occupationi[imin]-occupation))
    {
      imin = i;
    }
  }
  *shift=Shifti[imin];
  
  MaxIt=i;
  
  
  for(i=0;i<MaxIt+1;i++)
  {
    printf("%2d-th shift %22.15E occupation %22.15E deviation %22.15E\n",
           i, Shifti[i], occupationi[i], occupationi[i]-occupation);
  }
  
  return 0;
}

/*******************************
          Perturbations
 *******************************/


uint32_t XASBathBATrunc(SpectraSettings const O, BADoubleTriDiagonalMatrixType const BathBA,
                        BADoubleTriDiagonalMatrixType * const BathBATrunc)
{
  *BathBATrunc = BathBA;
  if (BathBA.con.N < O.NBath + 1)
  {
    printf("Warning: BathBA.con.N(%i)<XASSettings.NBath+1(%i+1)\n", BathBA.con.N, O.NBath);
  }
  else
  {
    BathBATrunc->con.N = O.NBath + 1;
  }
  if (BathBA.val.N < O.NBath + 1)
  {
    printf("Warning: BathBA.val.N(%i)<XASSettings.NBath+1(%i+1)\n", BathBA.val.N, O.NBath);
  }
  else
  {
    BathBATrunc->val.N = O.NBath + 1;
  }
  strcpy(BathBATrunc->Name, "Bath for ground-state");
  return 0;
}

uint32_t SpectraVPerturbation(SpectraSettings const O,
                              BADoubleTriDiagonalMatrixType const GBathBAa1g,
                              BADoubleTriDiagonalMatrixType const GBathBAb1g,
                              BADoubleTriDiagonalMatrixType const GBathBAb2g,
                              BADoubleTriDiagonalMatrixType const GBathBAeg,
                              OperatorType * const VPerturbation,
                              TriDiagonalMatrixType * const gPerturbation)
{
  uint32_t i, NFermion;
  uint16_t oadd[1];
  
  NFermion = 6+10+10+20*O.NBath;
  for (i = 0; i < 20; i++)
  {
    InitOperator(&(VPerturbation[i]));
    VPerturbation[i].NFermion = NFermion;
    VPerturbation[i].NBoson = 0;
    VPerturbation[i].NWord = Min(VPerturbation[i].NFermion, 1 + (VPerturbation[i].NFermion - 1) / 8);
    VPerturbation[i].QComplex = 0;
    InitOperatorFixedLength(&(VPerturbation[i]), 1, 1);
  }
  
  /* orbitals are ordered as dxy, dyz, dz2, dxz, dx2-y2 */
  // a1g: 3z2-1
  // b1g: x2-y2
  // b2g: xy
  // eg : yz, xz
  
  if (GBathBAa1g.val.N > (O.NBath + 1))
  {
    oadd[0] = CrAnMask + 26 + 10*(O.NBath-1) + 4;
    RealOperatorAddLader(&(VPerturbation[4]), oadd, 1, GBathBAa1g.val.b[O.NBath + 1]);
    oadd[0] = CrAnMask + 26 + 10*(O.NBath-1) + 5;
    RealOperatorAddLader(&(VPerturbation[5]), oadd, 1, GBathBAa1g.val.b[O.NBath + 1]);
  }
  if (GBathBAa1g.con.N > (O.NBath + 1))
  {
    oadd[0] = 26 + 10*(2*O.NBath-1) + 4;
    RealOperatorAddLader(&(VPerturbation[14]), oadd, 1, GBathBAa1g.con.b[O.NBath + 1]);
    oadd[0] = 26 + 10*(2*O.NBath-1) + 5;
    RealOperatorAddLader(&(VPerturbation[15]), oadd, 1, GBathBAa1g.con.b[O.NBath + 1]);
  }
  
  if (GBathBAb1g.val.N > (O.NBath + 1))
  {
    oadd[0] = CrAnMask + 26 + 10*(O.NBath-1) + 8;
    RealOperatorAddLader(&(VPerturbation[8]), oadd, 1, GBathBAb1g.val.b[O.NBath + 1]);
    oadd[0] = CrAnMask + 26 + 10*(O.NBath-1) + 9;
    RealOperatorAddLader(&(VPerturbation[9]), oadd, 1, GBathBAb1g.val.b[O.NBath + 1]);
  }
  if (GBathBAb1g.con.N > (O.NBath + 1))
  {
    oadd[0] = 26 + 10*(2*O.NBath-1) + 8;
    RealOperatorAddLader(&(VPerturbation[18]), oadd, 1, GBathBAb1g.con.b[O.NBath + 1]);
    oadd[0] = 26 + 10*(2*O.NBath-1) + 9;
    RealOperatorAddLader(&(VPerturbation[19]), oadd, 1, GBathBAb1g.con.b[O.NBath + 1]);
  }
  
  if (GBathBAb2g.val.N > (O.NBath + 1))
  {
    oadd[0] = CrAnMask + 26 + 10*(O.NBath-1) + 0;
    RealOperatorAddLader(&(VPerturbation[0]), oadd, 1, GBathBAb2g.val.b[O.NBath + 1]);
    oadd[0] = CrAnMask + 26 + 10*(O.NBath-1) + 1;
    RealOperatorAddLader(&(VPerturbation[1]), oadd, 1, GBathBAb2g.val.b[O.NBath + 1]);
  }
  if (GBathBAb2g.con.N > (O.NBath + 1))
  {
    oadd[0] = 26 + 10*(2*O.NBath-1) + 0;
    RealOperatorAddLader(&(VPerturbation[10]), oadd, 1, GBathBAb2g.con.b[O.NBath + 1]);
    oadd[0] = 26 + 10*(2*O.NBath-1) + 1;
    RealOperatorAddLader(&(VPerturbation[11]), oadd, 1, GBathBAb2g.con.b[O.NBath + 1]);
  }
  
  if (GBathBAeg.val.N > (O.NBath + 1))
  {
    oadd[0] = CrAnMask + 26 + 10*(O.NBath-1) + 2;
    RealOperatorAddLader(&(VPerturbation[2]), oadd, 1, GBathBAeg.val.b[O.NBath + 1]);
    oadd[0] = CrAnMask + 26 + 10*(O.NBath-1) + 3;
    RealOperatorAddLader(&(VPerturbation[3]), oadd, 1, GBathBAeg.val.b[O.NBath + 1]);
    oadd[0] = CrAnMask + 26 + 10*(O.NBath-1) + 6;
    RealOperatorAddLader(&(VPerturbation[6]), oadd, 1, GBathBAeg.val.b[O.NBath + 1]);
    oadd[0] = CrAnMask + 26 + 10*(O.NBath-1) + 7;
    RealOperatorAddLader(&(VPerturbation[7]), oadd, 1, GBathBAeg.val.b[O.NBath + 1]);
  }
  if (GBathBAeg.con.N > (O.NBath + 1))
  {
    oadd[0] = 26 + 10*(2*O.NBath-1) + 2;
    RealOperatorAddLader(&(VPerturbation[12]), oadd, 1, GBathBAeg.con.b[O.NBath + 1]);
    oadd[0] = 26 + 10*(2*O.NBath-1) + 3;
    RealOperatorAddLader(&(VPerturbation[13]), oadd, 1, GBathBAeg.con.b[O.NBath + 1]);
    oadd[0] = 26 + 10*(2*O.NBath-1) + 6;
    RealOperatorAddLader(&(VPerturbation[16]), oadd, 1, GBathBAeg.con.b[O.NBath + 1]);
    oadd[0] = 26 + 10*(2*O.NBath-1) + 7;
    RealOperatorAddLader(&(VPerturbation[17]), oadd, 1, GBathBAeg.con.b[O.NBath + 1]);
  }

  /* Create gPerturbation */
  if (GBathBAa1g.val.N > O.NBath + 1)
  {
    gPerturbation[4].N = GBathBAa1g.val.N - O.NBath - 1;
    gPerturbation[5].N = GBathBAa1g.val.N - O.NBath - 1;
    InitTriDiagonalMatrix(&(gPerturbation[4]));
    InitTriDiagonalMatrix(&(gPerturbation[5]));
    for (i = 1; i < GBathBAa1g.val.N - O.NBath - 1; i++)
    {
      gPerturbation[4].a[i] = -GBathBAa1g.val.a[i + O.NBath + 1];
      gPerturbation[4].b[i] = GBathBAa1g.val.b[i + O.NBath + 1];
      gPerturbation[5].a[i] = -GBathBAa1g.val.a[i + O.NBath + 1];
      gPerturbation[5].b[i] = GBathBAa1g.val.b[i + O.NBath + 1];
    }
    gPerturbation[4].a[GBathBAa1g.val.N - O.NBath - 1] = -GBathBAa1g.val.a[GBathBAa1g.val.N];
    gPerturbation[5].a[GBathBAa1g.val.N - O.NBath - 1] = -GBathBAa1g.val.a[GBathBAa1g.val.N];
  }
  else
  {
    gPerturbation[4].N = 1;
    gPerturbation[5].N = 1;
    InitTriDiagonalMatrix(&(gPerturbation[4]));
    InitTriDiagonalMatrix(&(gPerturbation[5]));
    gPerturbation[4].a[1] = 0;
    gPerturbation[5].a[1] = 0;
  }
  if (GBathBAa1g.con.N > O.NBath + 1)
  {
    gPerturbation[14].N = GBathBAa1g.con.N - O.NBath - 1;
    gPerturbation[15].N = GBathBAa1g.con.N - O.NBath - 1;
    InitTriDiagonalMatrix(&(gPerturbation[14]));
    InitTriDiagonalMatrix(&(gPerturbation[15]));
    for (i = 1; i < GBathBAa1g.con.N - O.NBath - 1; i++)
    {
      gPerturbation[14].a[i] = GBathBAa1g.con.a[i + O.NBath + 1];
      gPerturbation[14].b[i] = GBathBAa1g.con.b[i + O.NBath + 1];
      gPerturbation[15].a[i] = GBathBAa1g.con.a[i + O.NBath + 1];
      gPerturbation[15].b[i] = GBathBAa1g.con.b[i + O.NBath + 1];
    }
    gPerturbation[14].a[GBathBAa1g.con.N - O.NBath - 1] = GBathBAa1g.con.a[GBathBAa1g.con.N];
    gPerturbation[15].a[GBathBAa1g.con.N - O.NBath - 1] = GBathBAa1g.con.a[GBathBAa1g.con.N];
  }
  else
  {
    gPerturbation[14].N = 1;
    gPerturbation[15].N = 1;
    InitTriDiagonalMatrix(&(gPerturbation[14]));
    InitTriDiagonalMatrix(&(gPerturbation[15]));
    gPerturbation[14].a[1] = 0;
    gPerturbation[15].a[1] = 0;
  }
  
  if (GBathBAb1g.val.N > O.NBath + 1)
  {
    gPerturbation[8].N = GBathBAb1g.val.N - O.NBath - 1;
    gPerturbation[9].N = GBathBAb1g.val.N - O.NBath - 1;
    InitTriDiagonalMatrix(&(gPerturbation[8]));
    InitTriDiagonalMatrix(&(gPerturbation[9]));
    for (i = 1; i < GBathBAb1g.val.N - O.NBath - 1; i++)
    {
      gPerturbation[8].a[i] = -GBathBAb1g.val.a[i + O.NBath + 1];
      gPerturbation[8].b[i] = GBathBAb1g.val.b[i + O.NBath + 1];
      gPerturbation[9].a[i] = -GBathBAb1g.val.a[i + O.NBath + 1];
      gPerturbation[9].b[i] = GBathBAb1g.val.b[i + O.NBath + 1];
    }
    gPerturbation[8].a[GBathBAb1g.val.N - O.NBath - 1] = -GBathBAb1g.val.a[GBathBAb1g.val.N];
    gPerturbation[9].a[GBathBAb1g.val.N - O.NBath - 1] = -GBathBAb1g.val.a[GBathBAb1g.val.N];
  }
  else
  {
    gPerturbation[8].N = 1;
    gPerturbation[9].N = 1;
    InitTriDiagonalMatrix(&(gPerturbation[8]));
    InitTriDiagonalMatrix(&(gPerturbation[9]));
    gPerturbation[8].a[1] = 0;
    gPerturbation[9].a[1] = 0;
  }
  if (GBathBAb1g.con.N > O.NBath + 1)
  {
    gPerturbation[18].N = GBathBAb1g.con.N - O.NBath - 1;
    gPerturbation[19].N = GBathBAb1g.con.N - O.NBath - 1;
    InitTriDiagonalMatrix(&(gPerturbation[18]));
    InitTriDiagonalMatrix(&(gPerturbation[19]));
    for (i = 1; i < GBathBAb1g.con.N - O.NBath - 1; i++)
    {
      gPerturbation[18].a[i] = GBathBAb1g.con.a[i + O.NBath + 1];
      gPerturbation[18].b[i] = GBathBAb1g.con.b[i + O.NBath + 1];
      gPerturbation[19].a[i] = GBathBAb1g.con.a[i + O.NBath + 1];
      gPerturbation[19].b[i] = GBathBAb1g.con.b[i + O.NBath + 1];
    }
    gPerturbation[18].a[GBathBAb1g.con.N - O.NBath - 1] = GBathBAb1g.con.a[GBathBAb1g.con.N];
    gPerturbation[19].a[GBathBAb1g.con.N - O.NBath - 1] = GBathBAb1g.con.a[GBathBAb1g.con.N];
  }
  else
  {
    gPerturbation[18].N = 1;
    gPerturbation[19].N = 1;
    InitTriDiagonalMatrix(&(gPerturbation[18]));
    InitTriDiagonalMatrix(&(gPerturbation[19]));
    gPerturbation[18].a[1] = 0;
    gPerturbation[19].a[1] = 0;
  }
  
  if (GBathBAb2g.val.N > O.NBath + 1)
  {
    gPerturbation[0].N = GBathBAb2g.val.N - O.NBath - 1;
    gPerturbation[1].N = GBathBAb2g.val.N - O.NBath - 1;
    InitTriDiagonalMatrix(&(gPerturbation[0]));
    InitTriDiagonalMatrix(&(gPerturbation[1]));
    for (i = 1; i < GBathBAb2g.val.N - O.NBath - 1; i++)
    {
      gPerturbation[0].a[i] = -GBathBAb2g.val.a[i + O.NBath + 1];
      gPerturbation[0].b[i] = GBathBAb2g.val.b[i + O.NBath + 1];
      gPerturbation[1].a[i] = -GBathBAb2g.val.a[i + O.NBath + 1];
      gPerturbation[1].b[i] = GBathBAb2g.val.b[i + O.NBath + 1];
    }
    gPerturbation[0].a[GBathBAb2g.val.N - O.NBath - 1] = -GBathBAb2g.val.a[GBathBAb2g.val.N];
    gPerturbation[1].a[GBathBAb2g.val.N - O.NBath - 1] = -GBathBAb2g.val.a[GBathBAb2g.val.N];
  }
  else
  {
    gPerturbation[0].N = 1;
    gPerturbation[1].N = 1;
    InitTriDiagonalMatrix(&(gPerturbation[0]));
    InitTriDiagonalMatrix(&(gPerturbation[1]));
    gPerturbation[0].a[1] = 0;
    gPerturbation[1].a[1] = 0;
  }
  if (GBathBAb2g.con.N > O.NBath + 1)
  {
    gPerturbation[10].N = GBathBAb2g.con.N - O.NBath - 1;
    gPerturbation[11].N = GBathBAb2g.con.N - O.NBath - 1;
    InitTriDiagonalMatrix(&(gPerturbation[10]));
    InitTriDiagonalMatrix(&(gPerturbation[11]));
    for (i = 1; i < GBathBAb2g.con.N - O.NBath - 1; i++)
    {
      gPerturbation[10].a[i] = GBathBAb2g.con.a[i + O.NBath + 1];
      gPerturbation[10].b[i] = GBathBAb2g.con.b[i + O.NBath + 1];
      gPerturbation[11].a[i] = GBathBAb2g.con.a[i + O.NBath + 1];
      gPerturbation[11].b[i] = GBathBAb2g.con.b[i + O.NBath + 1];
    }
    gPerturbation[10].a[GBathBAb2g.con.N - O.NBath - 1] = GBathBAb2g.con.a[GBathBAb2g.con.N];
    gPerturbation[11].a[GBathBAb2g.con.N - O.NBath - 1] = GBathBAb2g.con.a[GBathBAb2g.con.N];
  }
  else
  {
    gPerturbation[10].N = 1;
    gPerturbation[11].N = 1;
    InitTriDiagonalMatrix(&(gPerturbation[10]));
    InitTriDiagonalMatrix(&(gPerturbation[11]));
    gPerturbation[10].a[1] = 0;
    gPerturbation[11].a[1] = 0;
  }
  
  if (GBathBAeg.val.N > O.NBath + 1)
  {
    gPerturbation[2].N = GBathBAeg.val.N - O.NBath - 1;
    gPerturbation[3].N = GBathBAeg.val.N - O.NBath - 1;
    gPerturbation[6].N = GBathBAeg.val.N - O.NBath - 1;
    gPerturbation[7].N = GBathBAeg.val.N - O.NBath - 1;
    InitTriDiagonalMatrix(&(gPerturbation[2]));
    InitTriDiagonalMatrix(&(gPerturbation[3]));
    InitTriDiagonalMatrix(&(gPerturbation[6]));
    InitTriDiagonalMatrix(&(gPerturbation[7]));
    for (i = 1; i < GBathBAeg.val.N - O.NBath - 1; i++)
    {
      gPerturbation[2].a[i] = -GBathBAeg.val.a[i + O.NBath + 1];
      gPerturbation[2].b[i] = GBathBAeg.val.b[i + O.NBath + 1];
      gPerturbation[3].a[i] = -GBathBAeg.val.a[i + O.NBath + 1];
      gPerturbation[3].b[i] = GBathBAeg.val.b[i + O.NBath + 1];
      gPerturbation[6].a[i] = -GBathBAeg.val.a[i + O.NBath + 1];
      gPerturbation[6].b[i] = GBathBAeg.val.b[i + O.NBath + 1];
      gPerturbation[7].a[i] = -GBathBAeg.val.a[i + O.NBath + 1];
      gPerturbation[7].b[i] = GBathBAeg.val.b[i + O.NBath + 1];
    }
    gPerturbation[2].a[GBathBAeg.val.N - O.NBath - 1] = -GBathBAeg.val.a[GBathBAeg.val.N];
    gPerturbation[3].a[GBathBAeg.val.N - O.NBath - 1] = -GBathBAeg.val.a[GBathBAeg.val.N];
    gPerturbation[6].a[GBathBAeg.val.N - O.NBath - 1] = -GBathBAeg.val.a[GBathBAeg.val.N];
    gPerturbation[7].a[GBathBAeg.val.N - O.NBath - 1] = -GBathBAeg.val.a[GBathBAeg.val.N];
  }
  else
  {
    gPerturbation[2].N = 1;
    gPerturbation[3].N = 1;
    gPerturbation[6].N = 1;
    gPerturbation[7].N = 1;
    InitTriDiagonalMatrix(&(gPerturbation[2]));
    InitTriDiagonalMatrix(&(gPerturbation[3]));
    InitTriDiagonalMatrix(&(gPerturbation[6]));
    InitTriDiagonalMatrix(&(gPerturbation[7]));
    gPerturbation[2].a[1] = 0;
    gPerturbation[3].a[1] = 0;
    gPerturbation[6].a[1] = 0;
    gPerturbation[7].a[1] = 0;
  }
  if (GBathBAeg.con.N > O.NBath + 1)
  {
    gPerturbation[12].N = GBathBAeg.con.N - O.NBath - 1;
    gPerturbation[13].N = GBathBAeg.con.N - O.NBath - 1;
    gPerturbation[16].N = GBathBAeg.con.N - O.NBath - 1;
    gPerturbation[17].N = GBathBAeg.con.N - O.NBath - 1;
    InitTriDiagonalMatrix(&(gPerturbation[12]));
    InitTriDiagonalMatrix(&(gPerturbation[13]));
    InitTriDiagonalMatrix(&(gPerturbation[16]));
    InitTriDiagonalMatrix(&(gPerturbation[17]));
    for (i = 1; i < GBathBAeg.con.N - O.NBath - 1; i++)
    {
      gPerturbation[12].a[i] = GBathBAeg.con.a[i + O.NBath + 1];
      gPerturbation[12].b[i] = GBathBAeg.con.b[i + O.NBath + 1];
      gPerturbation[13].a[i] = GBathBAeg.con.a[i + O.NBath + 1];
      gPerturbation[13].b[i] = GBathBAeg.con.b[i + O.NBath + 1];
      gPerturbation[16].a[i] = GBathBAeg.con.a[i + O.NBath + 1];
      gPerturbation[16].b[i] = GBathBAeg.con.b[i + O.NBath + 1];
      gPerturbation[17].a[i] = GBathBAeg.con.a[i + O.NBath + 1];
      gPerturbation[17].b[i] = GBathBAeg.con.b[i + O.NBath + 1];
    }
    gPerturbation[12].a[GBathBAeg.con.N - O.NBath - 1] = GBathBAeg.con.a[GBathBAeg.con.N];
    gPerturbation[13].a[GBathBAeg.con.N - O.NBath - 1] = GBathBAeg.con.a[GBathBAeg.con.N];
    gPerturbation[16].a[GBathBAeg.con.N - O.NBath - 1] = GBathBAeg.con.a[GBathBAeg.con.N];
    gPerturbation[17].a[GBathBAeg.con.N - O.NBath - 1] = GBathBAeg.con.a[GBathBAeg.con.N];
  }
  else
  {
    gPerturbation[12].N = 1;
    gPerturbation[13].N = 1;
    gPerturbation[16].N = 1;
    gPerturbation[17].N = 1;
    InitTriDiagonalMatrix(&(gPerturbation[12]));
    InitTriDiagonalMatrix(&(gPerturbation[13]));
    InitTriDiagonalMatrix(&(gPerturbation[16]));
    InitTriDiagonalMatrix(&(gPerturbation[17]));
    gPerturbation[12].a[1] = 0;
    gPerturbation[13].a[1] = 0;
    gPerturbation[16].a[1] = 0;
    gPerturbation[17].a[1] = 0;
  }
  
  return 0;
}


uint32_t get2p3dAtomicInfo(uint32_t Z, uint32_t dCount, double *Fdd, double *Fpd, double *Gpd, double *zeta2p, double *zeta3d)
{
    
    
  switch (Z)
  { 
    case 20:
      switch (dCount)
      {
        case 0:
          Fdd[1] = 0.0; Fdd[2] = 0.0; Fpd[1] = 3.794; 
          Gpd[0] = 2.513; Gpd[1] = 1.424; *zeta2p = 2.404; *zeta3d = 0.004; 
          return 0;
        default:
          break;
      }
    case 25:
      switch(dCount)
      {
        case 5:
          Fdd[1] = 10.316; Fdd[2] = 6.414; Fpd[1] = 6.321;
          Gpd[0] = 4.606; Gpd[1] = 2.618; *zeta2p = 6.846; *zeta3d = 0.040;
          return 0;
        default:
          break;
      }
    case 27:
      switch(dCount)
      {
        case 7:
          Fdd[1] = 11.605; Fdd[2] = 7.209; Fpd[1] = 7.260;
          Gpd[0] = 5.397; Gpd[1] = 3.069; *zeta2p = 9.748; *zeta3d = 0.066;
          return 0;
        default:
          break;
      }
    case 28:
      switch(dCount)
      {
        case 8:
          Fdd[1] = 12.234; Fdd[2] = 7.598; Fpd[1] = 7.721;
          Gpd[0] = 5.787; Gpd[1] = 3.291; *zeta2p = 11.507; *zeta3d = 0.083;
          return 0;
        default:
          break;
      }
    default:
      break;
  }
  printf("Configuration with Z = %d, d count of %d not found.",Z,dCount);
  return 1;
}

/*******************************
          XAS models
 *******************************/

int LigandField(uint32_t Z, uint32_t dCount)
{
  TimeType Time;
  InitializeTiming(&Time);
  BeginTiming("Total", &Time);
  printf("\n========================================================");
  printf("\n==========     start ligand field theory     ===========");
  printf("\n========================================================\n");
  uint32_t i, j;
  double Bmag[3], Hex[3], zeta2p, zeta3d;
  double Rydberg, red2P3D, red3D3D, Fdd[3], Fpd[2], Gpd[2];
  double Udd, Upd, Delta, ep, ed, el;
  uint32_t NFermion = 26;
  uint16_t Index2Pdn[3] = {0, 2, 4};
  uint16_t Index2Pup[3] = {1, 3, 5};
  uint16_t Index3Ddn[5] = {6, 8, 10, 12, 14};
  uint16_t Index3Dup[5] = {7, 9, 11, 13, 15};
  uint16_t IndexLDdn[5] = {16, 18, 20, 22, 24};
  uint16_t IndexLDup[5] = {17, 19, 21, 23, 25};
  uint16_t Index2P[6] = {0, 1, 2, 3, 4, 5};
  uint16_t Index3D[10] = {6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
  uint16_t IndexLD[10] = {16, 17, 18, 19, 20, 21, 22, 23, 24, 25};
  uint16_t jlist[10];
  double jRe[10], jIm[10];
  HarmonicsExpansion D4hAkmHop, D4hAkmCF;
  HarmonicsExpansion D4hAkma1g, D4hAkmb1g, D4hAkmb2g, D4hAkmeg;
  double HoppingEnergiesa1g[4]={1., 0., 0., 0.};
  double HoppingEnergiesb1g[4]={0., 1., 0., 0.};
  double HoppingEnergiesb2g[4]={0., 0., 1., 0.};
  double HoppingEnergieseg[4] ={0., 0., 0., 1.};
  double HoppingEnergies[4], CFEnergies[4];
  double Ea1g, Eb1g, Eb2g, Eeg, Hopa1g, Hopb1g, Hopb2g, Hopeg, pds, pdp;
  HarmonicsExpansion Vx, Vy, Vz, Vleft, Vright;
  int Vl[]={1, 1, 1};
  int Vm[]={-1, 0, 1};
  double VxRe[]={sqrt(0.5), 0.0, sqrt(0.5)};
  double VyRe[]={sqrt(0.5), 0.0, -sqrt(0.5)};
  double VzRe[]={0.0, 1.0, 0.0};
  double VlRe[]={1.0, 0.0, 0.0};
  double VrRe[]={0.0, 0.0, 1.0};
  OperatorType OppH, OppSx, OppSy, OppSz, OppSsqr;
  OperatorType OppLx, OppLy, OppLz, OppLsqr;
  OperatorType OppJx, OppJy, OppJz, OppJsqr;
  OperatorType Oppldots2p, Oppldots3d;
  OperatorType OppNthreehalf, OppNfivehalf;
  OperatorType OppUdd, OppUpd, OppNd, OppNL, OppNp;
  OperatorType OppNda1g, OppNdb1g, OppNdb2g, OppNdeg;
  OperatorType Opptemp, OppHop, OppCF;
  OperatorType TPES[10], TIPES[10], TCPES[6], TXASx, TXASy, TXASz, TXASl, TXASr, TCharge, TRIXS[3];
  double epsilon, zero;
  double energy, energysqr, Esoc, ECoulomb, EKin;
  double Sx, Sy, Sz, Ssqr, Lx, Ly, Lz, Lsqr, Jx, Jy, Jz, Jsqr;
  double Nd, Ndsqr, Nthreehalf, Nfivehalf, Nthreehalfsqr, Nfivehalfsqr;
  double Nda1g, Ndb1g, Ndb2g, Ndeg;
  uint32_t NBitsKey;
  WaveFunctionType psi, Opsi;
  uint16_t CreateString[26];
  TriDiagonalMatrixType SpecTri;
  SpectraType Spectra[35];
  MatrixType Density;
  double Hist[30], OccHist[10];
  uint16_t oadd = 0;
  uint32_t doPES, doIPES, doCPES, doXAS, doXiSpin, doXiCharge, doRIXS;
  uint32_t rotate=0;
  uint16_t Nelectrons;
  MatrixType Rotated;


  uint32_t numwRes = 20;
  double gamma, gammaRes, NKrylovXAS, NKrylovRIXS, wRes[numwRes],Gfwhm;
  SpectraType RIXSSpectra[3][numwRes];


  
  printf("\n========================================================");
  printf("\n==========          set parameters           ===========");
  printf("\n========================================================\n");
  
  doPES=0;
  doIPES=0;
  doCPES=0;
  doXAS=1;
  doXiSpin=0;
  doXiCharge=0;
  doRIXS=1;
  
  Rydberg = 13.6057;
  
  Bmag[0]= 0.000000; /* magnetic field in the x direction (units of eV) */
  Bmag[1]= 0.000000; /* magnetic field in the y direction (units of eV) */
  Bmag[2]= 0.000000; /* magnetic field in the z direction (units of eV) */
  Hex[0] = 0.000000; /* exchange field in the x direction (units of eV) */
  Hex[1] = 0.000000; /* exchange field in the y direction (units of eV) */
  Hex[2] = 0.000000; /* exchange field in the z direction (units of eV) */
  
  if(get2p3dAtomicInfo(Z,dCount,Fdd,Fpd,Gpd,&zeta2p,&zeta3d))
  {
    return 1;
  }
  
  Nelectrons = 16+dCount;
  
  red2P3D = 0.8;
  red3D3D = 0.8;
  
  
  Fdd[1] *= red3D3D;    /* F2 */
  Fdd[2] *= red3D3D;     /* F4 */
  
  Fpd[1] *= red2P3D;     /* F2 */
  
  Gpd[0] *= red2P3D;     /* G1 */
  Gpd[1] *= red2P3D;      /* G3 */

  Udd = 7.0;
  Upd = 8.5;
  Delta = -2.0;      /* energy difference between dn and dn+1 ligand hole */
  
  Fdd[0] = Udd + 2.0/63.0 * (Fdd[1] + Fdd[2]);
  Fpd[0] = Upd + 1.0/15.0 * Gpd[0] + 3.0/70.0 * Gpd[1];
  
/* For d8
   First, solve

    6ep +  8ed + 10eL + 28Udd + 48Upd = 0
    6ep +  9ed +  9eL + 36Udd + 54Upd = Delta
    5ep +  9ed + 10eL + 36Udd + 45Upd = 0
    5ep + 10ed +  9eL + 45Udd + 50Upd = Delta + Udd - Upd

   for ep, ed, eL, to get:

    ep = (1/12)( 5Delta + 18Udd - 81Upd)
    ed = (1/12)( 5Delta - 78Udd - 45Upd)
    eL = (1/12)(-7Delta + 18Udd + 27Upd)

  Now shift eL to zero, to get:

 * ed = Delta - 8Udd - 6Upd
 * eL = 0
*/

  
/* For d7
   First, solve

    6ep +  7ed + 10eL + 21Udd + 42Upd = 0
    6ep +  8ed +  9eL + 28Udd + 48Upd = Delta
    5ep +  8ed + 10eL + 28Udd + 40Upd = 0
    5ep +  9ed +  9eL + 36Udd + 45Upd = Delta + Udd - Upd

   for ep, ed, eL, to get:

    ep = (1/23)( 10Delta +  28Udd - 136Upd)
    ed = (1/23)( 10Delta - 133Udd -  90Upd)
    eL = (1/23)(-13Delta +  28Udd +  48Upd)

  Now shift eL to zero, to get:

    ep = 0
    ed = (1/23)(  0Delta - 161Udd + 46Upd)
    eL = (1/23)(-23Delta +  0Udd + 184Upd)
 
 * Now shift eL to zero, to have (note:could I have just shifted eL to zero before?)
 * ed = Delta - 7Udd - 6Upd
 * eL = 0
*/

  if(dCount == 1)
    ed = Delta - 1. * Udd - 6. * Upd;
  else if(dCount == 2)
    ed = Delta - 2. * Udd - 6. * Upd;
  else if(dCount == 3)
    ed = Delta - 3. * Udd - 6. * Upd;
  else if(dCount == 4)
    ed = Delta - 4. * Udd - 6. * Upd;
  else if(dCount == 5)
    ed = Delta - 5. * Udd - 6. * Upd;
  else if(dCount == 7)
    ed = Delta - 6. * Udd - 6. * Upd;
  else if(dCount == 7)
    ed = Delta - 7. * Udd - 6. * Upd;
  else if(dCount == 8)
    ed = Delta - 8. * Udd - 6. * Upd;
  else if(dCount == 9)
    ed = Delta - 9. * Udd - 6. * Upd;
  else if(dCount == 10)
    ed = Delta - 10. * Udd - 6. * Upd;
  else
    ed = 0;
  
  ep = 0;
  el = 0.0;
  
  epsilon = 1.0e-4;
  zero = 1.0e-4;

  
  
  Ea1g = 0.3;
  Eb1g = 0.3;
  Eb2g = -0.2;
  Eeg  = -0.2;
  
  CFEnergies[0] = Ea1g;
  CFEnergies[1] = Eb1g;
  CFEnergies[2] = Eb2g;
  CFEnergies[3] = Eeg;
  
  PotentialExpandedOnY("D4h", 2, CFEnergies, &D4hAkmCF);
  
  /*
  pds = -1.88;
  pdp = pds/2.17;
  Hopa1g = sqrt(3.0) * pds;
  Hopb1g = sqrt(3.0) * pds;
  Hopb2g = 2.0 * pdp;
  Hopeg  = 2.0 * pdp;
  */
  Hopa1g = -2.06;
  Hopb1g = -2.06;
  Hopb2g = 1.21;
  Hopeg  = 1.21;
  
  
  HoppingEnergies[0] = Hopa1g;
  HoppingEnergies[1] = Hopb1g;
  HoppingEnergies[2] = Hopb2g;
  HoppingEnergies[3] = Hopeg;
  
  PotentialExpandedOnY("D4h", 2, HoppingEnergies, &D4hAkmHop);
  PotentialExpandedOnY("D4h", 2, HoppingEnergiesa1g, &D4hAkma1g);
  PotentialExpandedOnY("D4h", 2, HoppingEnergiesb1g, &D4hAkmb1g);
  PotentialExpandedOnY("D4h", 2, HoppingEnergiesb2g, &D4hAkmb2g);
  PotentialExpandedOnY("D4h", 2, HoppingEnergieseg,  &D4hAkmeg );

  /* Polarization dependent transition operators for XAS and RIXS.
     These are created from potentials expanded on spherical harmonics.
     Here the expanded potentials are given.  */  
  Vx.N=3; Vx.QComplex=0; Vx.l=Vl; Vx.m=Vm; Vx.Re=VxRe;
  Vy.N=3; Vy.QComplex=0; Vy.l=Vl; Vy.m=Vm; Vy.Re=VyRe;
  Vz.N=3; Vz.QComplex=0; Vz.l=Vl; Vz.m=Vm; Vz.Re=VzRe;
  Vleft.N=3; Vleft.QComplex=0; Vleft.l=Vl; Vleft.m=Vm; Vleft.Re=VlRe;
  Vright.N=3; Vright.QComplex=0; Vright.l=Vl; Vright.m=Vm; Vright.Re=VrRe;
  
  printf("\n========================================================");
  printf("\n==========         define operators          ===========");
  printf("\n========================================================\n");
  BeginTiming("define operators", &Time);
  
  /* first define some operators */
  
  CreateOperatorSx(NFermion, 5, Index3Dup, Index3Ddn, &OppSx);
  CreateOperatorSy(NFermion, 5, Index3Dup, Index3Ddn, &OppSy);
  CreateOperatorSz(NFermion, 5, Index3Dup, Index3Ddn, &OppSz);
  CreateOperatorSsqr(NFermion, 5, Index3Dup, Index3Ddn, &OppSsqr);
  
  CreateOperatorLx(NFermion, 2, Index3Dup, Index3Ddn, &OppLx);
  CreateOperatorLy(NFermion, 2, Index3Dup, Index3Ddn, &OppLy);
  CreateOperatorLz(NFermion, 2, Index3Dup, Index3Ddn, &OppLz);
  CreateOperatorLsqr(NFermion, 2, Index3Dup, Index3Ddn, &OppLsqr);
  
  CreateOperatorJx(NFermion, 2, Index3Dup, Index3Ddn, &OppJx);
  CreateOperatorJy(NFermion, 2, Index3Dup, Index3Ddn, &OppJy);
  CreateOperatorJz(NFermion, 2, Index3Dup, Index3Ddn, &OppJz);
  CreateOperatorJsqr(NFermion, 2, Index3Dup, Index3Ddn, &OppJsqr);

  CreateOperatorldots(NFermion, 1, Index2Pup, Index2Pdn, &Oppldots2p);
  CreateOperatorldots(NFermion, 2, Index3Dup, Index3Ddn, &Oppldots3d);
  
  CreateOperatorCoulombSingleShell(NFermion, 2, Index3Dup, Index3Ddn, Fdd, &OppUdd);
  CreateOperatorCoulombTwoShellsConstantOccupation(NFermion, 1, 2, Index2Pup, Index2Pdn, Index3Dup, Index3Ddn, Fpd, Gpd, &OppUpd);
  
  CreateOperatorCF(NFermion, 2, 2, Index3Dup, Index3Ddn, Index3Dup, Index3Ddn, D4hAkmCF, &OppCF);
  
  CreateOperatorCF(NFermion, 2, 2, IndexLDup, IndexLDdn, Index3Dup, Index3Ddn, D4hAkmHop, &OppHop);
  CreateOperatorCF(NFermion, 2, 2, Index3Dup, Index3Ddn, IndexLDup, IndexLDdn, D4hAkmHop, &Opptemp);
  OperatorAdd(&OppHop, Opptemp, 1.0, 0.0); FreeOperator(&Opptemp);
  
  
//  CreateOperatorNumberShell(NFermion, 6, Index2P, &OppNp);
  
  CreateOperatorNumberShell(NFermion, 10, Index3D, &OppNd);
  CreateOperatorNumberShell(NFermion, 10, IndexLD, &OppNL);
  CreateOperatorCF(NFermion, 2, 2, Index3Dup, Index3Ddn, Index3Dup, Index3Ddn, D4hAkma1g, &OppNda1g);
  CreateOperatorCF(NFermion, 2, 2, Index3Dup, Index3Ddn, Index3Dup, Index3Ddn, D4hAkmb1g, &OppNdb1g);
  CreateOperatorCF(NFermion, 2, 2, Index3Dup, Index3Ddn, Index3Dup, Index3Ddn, D4hAkmb2g, &OppNdb2g);
  CreateOperatorCF(NFermion, 2, 2, Index3Dup, Index3Ddn, Index3Dup, Index3Ddn, D4hAkmeg,  &OppNdeg );
  OperatorCopy(&OppNthreehalf, Oppldots3d); /* <ls> = -3/2n(j=3/2)+n(j=5/2) */
  OperatorAdd(&OppNthreehalf, OppNd, -1.0, 0.0);
  OperatorScale(&OppNthreehalf, -0.4, 0.0);
  OperatorCopy(&OppNfivehalf, OppNd);
  OperatorAdd(&OppNfivehalf, OppNthreehalf, -1.0, 0.0);
  
  /* sum these operators to create the Hamiltonian */
  
  OperatorCopy(&OppH, OppUdd);
  OperatorAdd(&OppH, OppUpd, 1.0, 0.0);
  OperatorAdd(&OppH, Oppldots2p, zeta2p, 0.0);
  OperatorAdd(&OppH, Oppldots3d, zeta3d, 0.0);
  OperatorAdd(&OppH, OppCF, 1.0, 0.0);
  OperatorAdd(&OppH, OppHop, 1.0, 0.0);

//  OperatorAdd(&OppH, OppNp, ep, 0.0);
  
  OperatorAdd(&OppH, OppNd, ed, 0.0);
  OperatorAdd(&OppH, OppNL, el, 0.0);
  OperatorAdd(&OppH, OppSx, 2*Bmag[0]+Hex[0], 0.0);
  OperatorAdd(&OppH, OppSy, 2*Bmag[1]+Hex[1], 0.0);
  OperatorAdd(&OppH, OppSz, 2*Bmag[2]+Hex[2], 0.0);
  OperatorAdd(&OppH, OppLx, Bmag[0], 0.0);
  OperatorAdd(&OppH, OppLy, Bmag[1], 0.0);
  OperatorAdd(&OppH, OppLz, Bmag[2], 0.0);
  
  OperatorChop(&OppH, DBL_EPSILON*1000);
  
  PrintOperatorReadable(OppH);
  
  /* create the transition operator for XAS */
  
  CreateOperatorCF(NFermion, 2, 1, Index3Dup, Index3Ddn, Index2Pup, Index2Pdn, Vx, &TXASx);
  CreateOperatorCF(NFermion, 2, 1, Index3Dup, Index3Ddn, Index2Pup, Index2Pdn, Vy, &TXASy);
  CreateOperatorCF(NFermion, 2, 1, Index3Dup, Index3Ddn, Index2Pup, Index2Pdn, Vz, &TXASz);
  CreateOperatorCF(NFermion, 2, 1, Index3Dup, Index3Ddn, Index2Pup, Index2Pdn, Vleft, &TXASl);
  CreateOperatorCF(NFermion, 2, 1, Index3Dup, Index3Ddn, Index2Pup, Index2Pdn, Vright, &TXASr);

  /*Create decay operators for RIXS*/
  for(i = 0; i < 3; i++)
  {
    InitOperator(&TRIXS[i]);
  }
  OperatorConjugateTranspose(TXASx,&TRIXS[0]);
  OperatorConjugateTranspose(TXASy,&TRIXS[1]);
  OperatorConjugateTranspose(TXASz,&TRIXS[2]);


  EndTiming("define operators", &Time);
  
  if (rotate)
  {
    printf("\n========================================================");
    printf("\n==== Rotate the operators from lm to real orbital ======");
    printf("\n========================================================\n");
    
    Rotated.QComplex = 0; //1;
    Rotated.Ni = 26;
    Rotated.Nj = 26;
    InitMatrix(&Rotated);
    for (i=0; i<Rotated.Ni; i++)
    {
      for (j=0; j<Rotated.Nj; j++)
      {
        Rotated.VRe[i][j] = 0.0;
//        Rotated.VIm[i][j] = 0.0;
      }
    }



    for (i=0; i<6; i++)
    {
      Rotated.VRe[i][i] = 1.0;
    }
/*
    Rotated.VIm[0][0] =-sqrt(0.5);    Rotated.VRe[0][4] = sqrt(0.5);
    Rotated.VIm[1][1] =-sqrt(0.5);    Rotated.VRe[1][5] = sqrt(0.5);
    
    Rotated.VRe[2][2] = 1.0;
    Rotated.VRe[3][3] = 1.0;
    
    Rotated.VIm[4][0] =-sqrt(0.5);    Rotated.VRe[4][4] =-sqrt(0.5);
    Rotated.VIm[5][1] =-sqrt(0.5);    Rotated.VRe[5][5] =-sqrt(0.5);
  */  
    for (i=0; i<2; i++)
    {
	/*
      Rotated.VIm[ 6+i*10][ 6+i*10] =-sqrt(0.5);   Rotated.VRe[ 6+i*10][14+i*10] = sqrt(0.5);
      Rotated.VIm[ 7+i*10][ 7+i*10] =-sqrt(0.5);   Rotated.VRe[ 7+i*10][15+i*10] = sqrt(0.5);
      
      Rotated.VIm[ 8+i*10][ 8+i*10] =-sqrt(0.5);   Rotated.VRe[ 8+i*10][12+i*10] = sqrt(0.5);
      Rotated.VIm[ 9+i*10][ 9+i*10] =-sqrt(0.5);   Rotated.VRe[ 9+i*10][13+i*10] = sqrt(0.5);
      
      Rotated.VRe[10+i*10][10+i*10] = 1.0;
      Rotated.VRe[11+i*10][11+i*10] = 1.0;
      
      Rotated.VIm[12+i*10][ 8+i*10] =-sqrt(0.5);   Rotated.VRe[12+i*10][12+i*10] =-sqrt(0.5);
      Rotated.VIm[13+i*10][ 9+i*10] =-sqrt(0.5);   Rotated.VRe[13+i*10][13+i*10] =-sqrt(0.5);
      
      Rotated.VIm[14+i*10][ 6+i*10] = sqrt(0.5);   Rotated.VRe[14+i*10][14+i*10] = sqrt(0.5);
      Rotated.VIm[15+i*10][ 7+i*10] = sqrt(0.5);   Rotated.VRe[15+i*10][15+i*10] = sqrt(0.5);
	  */
	Rotated.VRe[ 6+i*10][ 6+i*10] =-sqrt(0.5);   Rotated.VRe[ 6+i*10][14+i*10] = sqrt(0.5);
    Rotated.VRe[ 7+i*10][ 7+i*10] =-sqrt(0.5);   Rotated.VRe[ 7+i*10][15+i*10] = sqrt(0.5);
    
    Rotated.VRe[ 8+i*10][ 8+i*10] =-sqrt(0.5);   Rotated.VRe[ 8+i*10][12+i*10] = sqrt(0.5);
    Rotated.VRe[ 9+i*10][ 9+i*10] =-sqrt(0.5);   Rotated.VRe[ 9+i*10][13+i*10] = sqrt(0.5);
    
    Rotated.VRe[10+i*10][10+i*10] = 1.0;
    Rotated.VRe[11+i*10][11+i*10] = 1.0;
    
    Rotated.VRe[12+i*10][ 8+i*10] =-sqrt(0.5);   Rotated.VRe[12+i*10][12+i*10] =-sqrt(0.5);
    Rotated.VRe[13+i*10][ 9+i*10] =-sqrt(0.5);   Rotated.VRe[13+i*10][13+i*10] =-sqrt(0.5);
    
    Rotated.VRe[14+i*10][ 6+i*10] = sqrt(0.5);   Rotated.VRe[14+i*10][14+i*10] = sqrt(0.5);
    Rotated.VRe[15+i*10][ 7+i*10] = sqrt(0.5);   Rotated.VRe[15+i*10][15+i*10] = sqrt(0.5);
    }
    RotateOperatorNormalOrder(&Opptemp, OppSx, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&OppSx); OperatorCopy(&OppSx, Opptemp); FreeOperator(&Opptemp);
    RotateOperatorNormalOrder(&Opptemp, OppSy, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&OppSy); OperatorCopy(&OppSy, Opptemp); FreeOperator(&Opptemp);
    RotateOperatorNormalOrder(&Opptemp, OppSz, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&OppSz); OperatorCopy(&OppSz, Opptemp); FreeOperator(&Opptemp);
    RotateOperatorNormalOrder(&Opptemp, OppSsqr, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&OppSsqr); OperatorCopy(&OppSsqr, Opptemp); FreeOperator(&Opptemp);
    
    RotateOperatorNormalOrder(&Opptemp, OppLx, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&OppLx); OperatorCopy(&OppLx, Opptemp); FreeOperator(&Opptemp);
    RotateOperatorNormalOrder(&Opptemp, OppLy, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&OppLy); OperatorCopy(&OppLy, Opptemp); FreeOperator(&Opptemp);
    RotateOperatorNormalOrder(&Opptemp, OppLz, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&OppLz); OperatorCopy(&OppLz, Opptemp); FreeOperator(&Opptemp);
    RotateOperatorNormalOrder(&Opptemp, OppLsqr, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&OppLsqr); OperatorCopy(&OppLsqr, Opptemp); FreeOperator(&Opptemp);
    
    RotateOperatorNormalOrder(&Opptemp, OppJx, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&OppJx); OperatorCopy(&OppJx, Opptemp); FreeOperator(&Opptemp);
    RotateOperatorNormalOrder(&Opptemp, OppJy, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&OppJy); OperatorCopy(&OppJy, Opptemp); FreeOperator(&Opptemp);
    RotateOperatorNormalOrder(&Opptemp, OppJz, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&OppJz); OperatorCopy(&OppJz, Opptemp); FreeOperator(&Opptemp);
    RotateOperatorNormalOrder(&Opptemp, OppJsqr, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&OppJsqr); OperatorCopy(&OppJsqr, Opptemp); FreeOperator(&Opptemp);
    
    RotateOperatorNormalOrder(&Opptemp, Oppldots2p, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&Oppldots2p); OperatorCopy(&Oppldots2p, Opptemp); FreeOperator(&Opptemp);
    RotateOperatorNormalOrder(&Opptemp, Oppldots3d, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&Oppldots3d); OperatorCopy(&Oppldots3d, Opptemp); FreeOperator(&Opptemp);
    
    RotateOperatorNormalOrder(&Opptemp, OppUdd, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&OppUdd); OperatorCopy(&OppUdd, Opptemp); FreeOperator(&Opptemp);
    RotateOperatorNormalOrder(&Opptemp, OppUpd, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&OppUpd); OperatorCopy(&OppUpd, Opptemp); FreeOperator(&Opptemp);
    RotateOperatorNormalOrder(&Opptemp, OppCF, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&OppCF); OperatorCopy(&OppCF, Opptemp); FreeOperator(&Opptemp);
    RotateOperatorNormalOrder(&Opptemp, OppHop, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&OppHop); OperatorCopy(&OppHop, Opptemp); FreeOperator(&Opptemp);
    
    RotateOperatorNormalOrder(&Opptemp, OppNd, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&OppNd); OperatorCopy(&OppNd, Opptemp); FreeOperator(&Opptemp);
    RotateOperatorNormalOrder(&Opptemp, OppNda1g, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&OppNda1g); OperatorCopy(&OppNda1g, Opptemp); FreeOperator(&Opptemp);
    RotateOperatorNormalOrder(&Opptemp, OppNdb1g, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&OppNdb1g); OperatorCopy(&OppNdb1g, Opptemp); FreeOperator(&Opptemp);
    RotateOperatorNormalOrder(&Opptemp, OppNdb2g, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&OppNdb2g); OperatorCopy(&OppNdb2g, Opptemp); FreeOperator(&Opptemp);
    RotateOperatorNormalOrder(&Opptemp, OppNdeg, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&OppNdeg); OperatorCopy(&OppNdeg, Opptemp); FreeOperator(&Opptemp);
    RotateOperatorNormalOrder(&Opptemp, OppNthreehalf, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&OppNthreehalf); OperatorCopy(&OppNthreehalf, Opptemp); FreeOperator(&Opptemp);
    RotateOperatorNormalOrder(&Opptemp, OppNfivehalf, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&OppNfivehalf); OperatorCopy(&OppNfivehalf, Opptemp); FreeOperator(&Opptemp);
    RotateOperatorNormalOrder(&Opptemp, OppNL, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&OppNL); OperatorCopy(&OppNL, Opptemp); FreeOperator(&Opptemp);
    RotateOperatorNormalOrder(&Opptemp, OppH, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&OppH); OperatorCopy(&OppH, Opptemp); FreeOperator(&Opptemp);
    
    
    RotateOperatorNormalOrder(&Opptemp, TXASx, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&TXASx); OperatorCopy(&TXASx, Opptemp); FreeOperator(&Opptemp);
    RotateOperatorNormalOrder(&Opptemp, TXASy, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&TXASy); OperatorCopy(&TXASy, Opptemp); FreeOperator(&Opptemp);
    RotateOperatorNormalOrder(&Opptemp, TXASz, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&TXASz); OperatorCopy(&TXASz, Opptemp); FreeOperator(&Opptemp);
    RotateOperatorNormalOrder(&Opptemp, TXASl, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&TXASl); OperatorCopy(&TXASl, Opptemp); FreeOperator(&Opptemp);
    RotateOperatorNormalOrder(&Opptemp, TXASr, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&TXASr); OperatorCopy(&TXASr, Opptemp); FreeOperator(&Opptemp);


    for(j = 0; j < 3; j++)
    {
      RotateOperatorNormalOrder(&Opptemp, TRIXS[j], Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
      FreeOperator(&TRIXS[j]); OperatorCopy(&TRIXS[j], Opptemp); FreeOperator(&Opptemp);
    }
  }
  
  printf("\n========================================================");
  printf("\n==========         find groundstate          ===========");
  printf("\n========================================================\n");
  BeginTiming("find groundstate", &Time);
  
  NBitsKey=16;
  //2p electrons
  for (i=0; i<6; i++)
  {
    CreateString[i] = i;
  }
  //make a highspin 3d config and fill the ligand orbitals
  j = 0;
  for (i=6; i<Nelectrons; i++)
  {
    if((i+j)%2 == 0 && i < (6+dCount) && j < (10-dCount) )
    {
      j++;
      printf("Point 1, i %d j %d\n",i,j);
      CreateString[i] = i+j;
    }
    else if(i < (6+dCount))
    {
      printf("Point 2, i %d j %d\n",i,j);
      CreateString[i] = i+j;
    }
    else
    {
      printf("Point 3, i %d j %d\n",i,j);
      CreateString[i] = i + (10-dCount);
    }
  }
 
  printf("\n\n\n\n");
  printf("Wavefunction string: ");
  for(i = 0; i < Nelectrons; i++)
  {
    printf("%d ", CreateString[i]);
  }
  printf("\n\n\n\n");


  CreateWaveFunctionFromString(&psi, CreateString, Nelectrons, NFermion, 0, NBitsKey);
   
  
  BeginTiming("LanczosGroundState", &Time);
  LanczosGroundState(OppH, &psi, &epsilon, zero, 0, &psi, 0);

  WaveFunctionRandomize(&psi);

  LanczosGroundState(OppH, &psi, &epsilon, zero, 0, &psi, 0);
  
  WaveFunctionRandomize(&psi);
  
  LanczosGroundState(OppH, &psi, &epsilon, zero, 0, &psi, 0);
  LanczosGroundState(OppH, &psi, &epsilon, zero, 0, &psi, 0);
  LanczosGroundState(OppH, &psi, &epsilon, zero, 0, &psi, 0);
  LanczosGroundState(OppH, &psi, &epsilon, zero, 0, &psi, 0);
  LanczosGroundState(OppH, &psi, &epsilon, zero, 0, &psi, 0);
  LanczosGroundState(OppH, &psi, &epsilon, zero, 0, &psi, 0);
  EndTiming("LanczosGroundState", &Time);
  
  CalculatePrefactorHistogram(Hist, psi, 2.0, 30);
  printf("\nPrefactor Histogram of groundstate\n");
  printf("Interval                                 | percentage of states   | alpha^2 * number of states \n");
  for(i=0;i<30;i++)
  {
    printf("2^-%2i (%10.3E) to 2^-%2i (%10.3E) | %22.15E | %22.15E \n",
           i, pow(0.5,(double)i), i+1, pow(0.5,(double)(i+1)),
           Hist[i], Hist[i]*pow(0.5,(double)(2*i)));
  }
  
  CalculateOccupationHistogram(OccHist, psi, Index3D, 10);
  
  PsiHermitianOperatorPsi(OppH, psi, &energy);
  PsiHermitianOperatorSqrPsi(OppH, psi, &energysqr);
  PsiHermitianOperatorPsi(Oppldots3d, psi, &Esoc);
  Esoc *= zeta3d;
  PsiHermitianOperatorPsi(OppUdd, psi, &ECoulomb);
  PsiHermitianOperatorPsi(OppHop, psi, &EKin);
  PsiHermitianOperatorPsi(OppSx, psi, &Sx);
  PsiHermitianOperatorPsi(OppSy, psi, &Sy);
  PsiHermitianOperatorPsi(OppSz, psi, &Sz);
  PsiHermitianOperatorPsi(OppSsqr, psi, &Ssqr);
  PsiHermitianOperatorPsi(OppLx, psi, &Lx);
  PsiHermitianOperatorPsi(OppLy, psi, &Ly);
  PsiHermitianOperatorPsi(OppLz, psi, &Lz);
  PsiHermitianOperatorPsi(OppLsqr, psi, &Lsqr);
  PsiHermitianOperatorPsi(OppJx, psi, &Jx);
  PsiHermitianOperatorPsi(OppJy, psi, &Jy);
  PsiHermitianOperatorPsi(OppJz, psi, &Jz);
  PsiHermitianOperatorPsi(OppJsqr, psi, &Jsqr);
  PsiHermitianOperatorPsi(OppNd, psi, &Nd);
  PsiHermitianOperatorPsi(OppNda1g, psi, &Nda1g);
  PsiHermitianOperatorPsi(OppNdb1g, psi, &Ndb1g);
  PsiHermitianOperatorPsi(OppNdb2g, psi, &Ndb2g);
  PsiHermitianOperatorPsi(OppNdeg,  psi, &Ndeg );
  PsiHermitianOperatorPsi(OppNthreehalf, psi, &Nthreehalf);
  PsiHermitianOperatorPsi(OppNfivehalf, psi, &Nfivehalf);
  PsiHermitianOperatorSqrPsi(OppNd, psi, &Ndsqr);
  PsiHermitianOperatorSqrPsi(OppNthreehalf, psi, &Nthreehalfsqr);
  PsiHermitianOperatorSqrPsi(OppNfivehalf, psi, &Nfivehalfsqr);
  EndTiming("find groundstate", &Time);
  PrintTiming(Time);
  
  printf("\n========================================================");
  printf("\n==========   Groundstte expectation values   ===========");
  printf("\n========================================================\n");
  printf("Energy              %22.15E  (variance %22.15E)\n", energy, energysqr-energy*energy);
  printf("Kinetic energy      %22.15E\n", EKin);
  printf("Relativistic energy %22.15E\n", Esoc);
  printf("Coulomb energy (3d) %22.15E\n", ECoulomb);
  printf("Sx %22.15E Sy %22.15E Sz %22.15E Ssqr %22.15E\n", Sx, Sy, Sz, Ssqr);
  printf("Lx %22.15E Ly %22.15E Lz %22.15E Lsqr %22.15E\n", Lx, Ly, Lz, Lsqr);
  printf("Jx %22.15E Jy %22.15E Jz %22.15E Jsqr %22.15E\n", Jx, Jy, Jz, Jsqr);
  printf("Nd   %22.15E <Nd^2>   %22.15E\n", Nd, Ndsqr);
  printf("Nda1g %22.14E\n", Nda1g);
  printf("Ndb1g %22.14E\n", Ndb1g);
  printf("Ndb2g %22.14E\n", Ndb2g);
  printf("Ndeg  %22.14E\n", Ndeg);
  printf("N3/2 %22.15E <N3/2^2> %22.15E\n", Nthreehalf, Nthreehalfsqr);
  printf("N5/2 %22.15E <N5/2^2> %22.15E\n", Nfivehalf, Nfivehalfsqr);
  printf("\n");
  
  printf("Occupation of d-shell\n");
  for(i=0;i<=10;i++)
  {
    printf("d^%2i | %22.15E\n", i, OccHist[i]);
  }
  printf("\n");
  
  printf("diagonal part of density matrix\n");
  CalculateDensityMatrix(&Density, psi);
  for(i=0;i<Density.Ni;i++)
  {
    printf("%i %22.15E\n",i,Density.VRe[i][i]);
  }
  printf("\n");
  PrintMatrix(Density);
  
  if(doXAS)
  {
    printf("\n========================================================");
    printf("\n==========                XAS                ===========");
    printf("\n========================================================\n");
    BeginTiming("XAS", &Time);

    SpecTri.N = 220;
    InitTriDiagonalMatrix(&SpecTri);
     
    for(i=26;i<31;i++)
    {
      Spectra[i].N=7000;
      Spectra[i].Emin=-85.0;
      Spectra[i].Emax=-30.0;
      Spectra[i].dE = (Spectra[i].Emax - Spectra[i].Emin) / Spectra[i].N;
      Spectra[i].Gamma = 0.4;
      SpectrumAllocate(Spectra[i], &Spectra[i]);
    }
    
    BeginTiming("spectrum", &Time);
    OperatorPsi(TXASx, psi, &Opsi);
    LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[26].Gamma);
    TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[26]);
    FilePrintPlainTriDiagonalMatrix(SpecTri, "LigandField_XASTri_x", "w");
    EndTiming("spectrum", &Time);
    
    BeginTiming("spectrum", &Time);
    OperatorPsi(TXASy, psi, &Opsi);
    LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[27].Gamma);
    TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[27]);
    FilePrintPlainTriDiagonalMatrix(SpecTri, "LigandField_XASTri_y", "w");
    EndTiming("spectrum", &Time);
    
    BeginTiming("spectrum", &Time);
    OperatorPsi(TXASz, psi, &Opsi);
    LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[28].Gamma);
    TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[28]);
    FilePrintPlainTriDiagonalMatrix(SpecTri, "LigandField_XASTri_z", "w");
    EndTiming("spectrum", &Time);
    
    BeginTiming("spectrum", &Time);
    OperatorPsi(TXASl, psi, &Opsi);
    LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[29].Gamma);
    TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[29]);
    EndTiming("spectrum", &Time);
    
    BeginTiming("spectrum", &Time);
    OperatorPsi(TXASr, psi, &Opsi);
    LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[30].Gamma);
    TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[30]);
    EndTiming("spectrum", &Time);
    
    FilePrintSpectra(&Spectra[26], 5, "LigandField_XAS","w");
    EndTiming("XAS", &Time);
    PrintTiming(Time);
  }
  if(doRIXS)
  {
    printf("\n========================================================");
    printf("\n==========               RIXS                ===========");
    printf("\n========================================================\n");
    BeginTiming("RIXS", &Time);
//    for(i = 0; i < 3; i++)
//    {
//      InitOperator(&TRIXS[i]);
//    }
//    OperatorConjugateTranspose(TXASx,&TRIXS[0]);
//    OperatorConjugateTranspose(TXASy,&TRIXS[1]);
//    OperatorConjugateTranspose(TXASz,&TRIXS[2]);

    gamma = 0.1;
    Gfwhm = 0.3;
    gammaRes = 0.4;
    NKrylovXAS = 75;
    NKrylovRIXS = 75;

    for(i = 0; i < numwRes; i++)
      wRes[i] = -69.0 + ((double)i)*0.1;
  
    for(i=0;i<3;i++)
    {
      for(j=0; j < numwRes; j++)
      {
        RIXSSpectra[i][j].N=4601;
        RIXSSpectra[i][j].Emin=-5.0;
        RIXSSpectra[i][j].Emax=18.0;
        RIXSSpectra[i][j].dE = (RIXSSpectra[i][j].Emax - RIXSSpectra[i][j].Emin) / RIXSSpectra[i][j].N;
        RIXSSpectra[i][j].Gamma = gamma;
        SpectrumAllocate(RIXSSpectra[i][j], &RIXSSpectra[i][j]);
      }
    }

    BeginTiming("spectrum", &Time);
//    ResonantSpectra(OppH,TXASx,TRIXS,3,psi,energy,epsilon,&Spectra[26],-77.4,0.4,0.1,70,50);
    ResonantSpectra(OppH,TXASx,TRIXS,3,psi,energy,epsilon,numwRes,RIXSSpectra,wRes,gammaRes,gamma,Gfwhm, NKrylovXAS,NKrylovRIXS);



    EndTiming("spectrum", &Time);
printf("\nPoint1:\n"); fflush(stdout);
    FilePrintRIXSSpectra(numwRes,RIXSSpectra,wRes, 3, "LigandField_RIXS","w");
printf("\nPoint2:\n"); fflush(stdout);
    EndTiming("RIXS", &Time);
    PrintTiming(Time);
  }
  
  EndTiming("Total", &Time);
  PrintTiming(Time);
  return 0;
}

int ImpurityModel3d(uint32_t Z, uint32_t dCount)
{
  TimeType Time;
  InitializeTiming(&Time);
  BeginTiming("Total", &Time);
  printf("\n========================================================");
  printf("\n==========     start impurity model          ===========");
  printf("\n========================================================\n");
  fflush(stdout);

  uint32_t i, j;
  double Bmag[3], Hex[3], zeta2p, zeta3d;
  double red2P3D, red3D3D, Fdd[3], Fpd[2], Gpd[2];
  double Udd, Upd, Delta, ep, ed, el;
  uint32_t NFermion;
  uint16_t Index2Pdn[3] = {0, 2, 4};
  uint16_t Index2Pup[3] = {1, 3, 5};
  uint16_t Index3Ddn[5] = {6, 8, 10, 12, 14};
  uint16_t Index3Dup[5] = {7, 9, 11, 13, 15};
  uint16_t IndexLDvcup[5], IndexLDvcdn[5], IndexLDvcnup[5], IndexLDvcndn[5];
  uint16_t Index3D[10] = {6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
  uint16_t IndexD[10];
  uint16_t jlist[10];
  double jRe[10], jIm[10];
  HarmonicsExpansion D4hAkmHop, D4hAkmCF;
  HarmonicsExpansion D4hAkma1g, D4hAkmb1g, D4hAkmb2g, D4hAkmeg;
  double HoppingEnergiesa1g[4]={1., 0., 0., 0.};
  double HoppingEnergiesb1g[4]={0., 1., 0., 0.};
  double HoppingEnergiesb2g[4]={0., 0., 1., 0.};
  double HoppingEnergieseg[4] ={0., 0., 0., 1.};
  double HoppingEnergies[4], CFEnergies[4];
  double Ea1g, Eb1g, Eb2g, Eeg, Hopa1g, Hopb1g, Hopb2g, Hopeg, pds, pdp;
  HarmonicsExpansion Vx, Vy, Vz, Vleft, Vright;
  int Vl[]={1, 1, 1};
  int Vm[]={-1, 0, 1};
  double VxRe[]={sqrt(0.5), 0.0, sqrt(0.5)};
  double VyRe[]={sqrt(0.5), 0.0, -sqrt(0.5)};
  double VzRe[]={0.0, 1.0, 0.0};
  double VlRe[]={1.0, 0.0, 0.0};
  double VrRe[]={0.0, 0.0, 1.0};
  OperatorType OppH, OppSx, OppSy, OppSz, OppSsqr;
  OperatorType OppLx, OppLy, OppLz, OppLsqr;
  OperatorType OppJx, OppJy, OppJz, OppJsqr;
  OperatorType Oppldots2p, Oppldots3d;
  OperatorType OppNthreehalf, OppNfivehalf;
  OperatorType OppUdd, OppUpd, OppNd, OppNp;
  OperatorType OppNda1g, OppNdb1g, OppNdb2g, OppNdeg;
  OperatorType Opptemp, OppCF;
  OperatorType TXASx, TXASy, TXASz, TXASl, TXASr, TCharge, TRIXS[3];
  OperatorType *OppHybv, *OppNv;
  double epsilon, zero;
  double energy, energysqr, Esoc, ECoulomb, EKin;
  double Sx, Sy, Sz, Ssqr, Lx, Ly, Lz, Lsqr, Jx, Jy, Jz, Jsqr;
  double Nd, Ndsqr, Nthreehalf, Nfivehalf, Nthreehalfsqr, Nfivehalfsqr;
  double Nda1g, Ndb1g, Ndb2g, Ndeg;
  uint32_t NBitsKey;
  WaveFunctionType psi, Opsi;
  uint16_t *CreateString;
  TriDiagonalMatrixType SpecTri;
  SpectraType Spectra[35];
  MatrixType Density;
  double Hist[30], OccHist[10];
  uint16_t oadd = 0;
  uint32_t doXAS, doRIXS;
  uint32_t rotate=0;
  uint16_t Nelectrons;
  MatrixType Rotated;
  uint32_t NBathv;

  uint32_t numwRes = 20;
  double gamma, gammaRes, NKrylovXAS, NKrylovRIXS, wRes[numwRes],Gfwhm;
  SpectraType RIXSSpectra[3][numwRes];
  
  uint16_t star = 1;

  printf("\n========================================================");
  printf("\n==========          set parameters           ===========");
  printf("\n========================================================\n");
  
  doXAS=1; doRIXS=1;
  
  Bmag[0]= 0.000000; /* magnetic field in the x direction (units of eV) */
  Bmag[1]= 0.000000; /* magnetic field in the y direction (units of eV) */
  Bmag[2]= 0.000000; /* magnetic field in the z direction (units of eV) */
  Hex[0] = 0.000000; /* exchange field in the x direction (units of eV) */
  Hex[1] = 0.000000; /* exchange field in the y direction (units of eV) */
  Hex[2] = 0.000000; /* exchange field in the z direction (units of eV) */
  
  if(get2p3dAtomicInfo(Z,dCount,Fdd,Fpd,Gpd,&zeta2p,&zeta3d))
  {
    return 1;
  }
  red2P3D = 0.8; red3D3D = 0.8;
  Fdd[1] *= red3D3D; Fdd[2] *= red3D3D; Fpd[1] *= red2P3D; 
  Gpd[0] *= red2P3D; Gpd[1] *= red2P3D;

  Udd = 7.0;
  Upd = 8.5;
  Delta = 4.0;      /* energy difference between dn and dn+1 ligand hole */
  
  Fdd[0] = Udd + 2.0/63.0 * (Fdd[1] + Fdd[2]);
  Fpd[0] = Upd + 1.0/15.0 * Gpd[0] + 3.0/70.0 * Gpd[1];
  
  /*
   For d8
   
    6ep + 8ed + 100el +28U + 48Q = 0
    6ep + 9ed + 99el + 36U + 54Q = Delta
    5ep + 9ed + 100el +36U + 45Q = 0
    5ep + 10ed + 99el + 45U + 50Q = D + U - Q
   
   Solve to get
   
     ep = (1/57)(50Delta + 18U -486Q)
     ed = (1/57)(50Delta - 438U - 315Q)
     el = (1/57)(-7Delta + 18U + 27Q)
   
   Shift el to zero to get
   * ed = Delta - 8U -6Q
   
   */
  if(dCount == 1)
    ed = Delta - 1. * Udd - 6. * Upd;
  else if(dCount == 2)
    ed = Delta - 2. * Udd - 6. * Upd;
  else if(dCount == 3)
    ed = Delta - 3. * Udd - 6. * Upd;
  else if(dCount == 4)
    ed = Delta - 4. * Udd - 6. * Upd;
  else if(dCount == 5)
    ed = Delta - 5. * Udd - 6. * Upd;
  else if(dCount == 7)
    ed = Delta - 6. * Udd - 6. * Upd;
  else if(dCount == 7)
    ed = Delta - 7. * Udd - 6. * Upd;
  else if(dCount == 8)
    ed = Delta - 8. * Udd - 6. * Upd;
  else if(dCount == 9)
    ed = Delta - 9. * Udd - 6. * Upd;
  else if(dCount == 10)
    ed = Delta - 10. * Udd - 6. * Upd;
  else
    ed = 0;
  
  ep = 0;  el = 0.0;
  
  NBathv = 10;
  Nelectrons = 6 + dCount + 10*NBathv;
  NFermion = 6 + 10 + 10*(NBathv);
  
  epsilon = 1.0e-4; zero = 1.0e-4;

  NBitsKey=16;
  CreateString=calloc(6+10+10*NBathv, sizeof(uint16_t));
  
  /*2p electrons*/
  for (i=0; i<6; i++)
    CreateString[i] = i;
  
  /*make a highspin 3d config and fill the ligand orbitals*/
  j = 0;
  for (i=6; i<Nelectrons; i++)
  {
    if((i+j)%2 == 0 && i < (6+dCount) && j < (10-dCount) )
    {
      j++; CreateString[i] = i+j;
    }
    else if(i < (6+dCount))
    {
      CreateString[i] = i+j;
    }
    else
    {
      CreateString[i] = i + (10-dCount);
    }
  }
 
  
  Ea1g = 0.3; Eb1g = 0.3; Eb2g = -0.2; Eeg  = -0.2;
  CFEnergies[0] = Ea1g; CFEnergies[1] = Eb1g; CFEnergies[2] = Eb2g; CFEnergies[3] = Eeg;
  
  PotentialExpandedOnY("D4h", 2, CFEnergies, &D4hAkmCF);

  Hopa1g = -2.06; Hopb1g = -2.06; Hopb2g = 1.21; Hopeg  = 1.21;
  
  if(star)
  {
    Hopa1g /= sqrt(NBathv);
    Hopb1g /= sqrt(NBathv);
    Hopb2g /= sqrt(NBathv);
    Hopeg /= sqrt(NBathv);
    
  }
  HoppingEnergies[0] = Hopa1g; HoppingEnergies[1] = Hopb1g;
  HoppingEnergies[2] = Hopb2g; HoppingEnergies[3] = Hopeg;
  
  PotentialExpandedOnY("D4h", 2, HoppingEnergies, &D4hAkmHop);
  PotentialExpandedOnY("D4h", 2, HoppingEnergiesa1g, &D4hAkma1g);
  PotentialExpandedOnY("D4h", 2, HoppingEnergiesb1g, &D4hAkmb1g);
  PotentialExpandedOnY("D4h", 2, HoppingEnergiesb2g, &D4hAkmb2g);
  PotentialExpandedOnY("D4h", 2, HoppingEnergieseg,  &D4hAkmeg );

  /* Polarization dependent transition operators for XAS and RIXS.
     These are created from potentials expanded on spherical harmonics.
     Here the expanded potentials are given.  */  
  Vx.N=3; Vx.QComplex=0; Vx.l=Vl; Vx.m=Vm; Vx.Re=VxRe;
  Vy.N=3; Vy.QComplex=0; Vy.l=Vl; Vy.m=Vm; Vy.Re=VyRe;
  Vz.N=3; Vz.QComplex=0; Vz.l=Vl; Vz.m=Vm; Vz.Re=VzRe;
  Vleft.N=3; Vleft.QComplex=0; Vleft.l=Vl; Vleft.m=Vm; Vleft.Re=VlRe;
  Vright.N=3; Vright.QComplex=0; Vright.l=Vl; Vright.m=Vm; Vright.Re=VrRe;
  
  printf("\n========================================================");
  printf("\n==========         define operators          ===========");
  printf("\n========================================================\n");
  BeginTiming("define operators", &Time);
  
  CreateOperatorSx(NFermion, 5, Index3Dup, Index3Ddn, &OppSx);
  CreateOperatorSy(NFermion, 5, Index3Dup, Index3Ddn, &OppSy);
  CreateOperatorSz(NFermion, 5, Index3Dup, Index3Ddn, &OppSz);
  CreateOperatorSsqr(NFermion, 5, Index3Dup, Index3Ddn, &OppSsqr);
  
  CreateOperatorLx(NFermion, 2, Index3Dup, Index3Ddn, &OppLx);
  CreateOperatorLy(NFermion, 2, Index3Dup, Index3Ddn, &OppLy);
  CreateOperatorLz(NFermion, 2, Index3Dup, Index3Ddn, &OppLz);
  CreateOperatorLsqr(NFermion, 2, Index3Dup, Index3Ddn, &OppLsqr);
  
  CreateOperatorJx(NFermion, 2, Index3Dup, Index3Ddn, &OppJx);
  CreateOperatorJy(NFermion, 2, Index3Dup, Index3Ddn, &OppJy);
  CreateOperatorJz(NFermion, 2, Index3Dup, Index3Ddn, &OppJz);
  CreateOperatorJsqr(NFermion, 2, Index3Dup, Index3Ddn, &OppJsqr);

  CreateOperatorldots(NFermion, 1, Index2Pup, Index2Pdn, &Oppldots2p);
  CreateOperatorldots(NFermion, 2, Index3Dup, Index3Ddn, &Oppldots3d);
  
  CreateOperatorCoulombSingleShell(NFermion, 2, Index3Dup, Index3Ddn, Fdd, &OppUdd);
  CreateOperatorCoulombTwoShellsConstantOccupation(NFermion, 1, 2, Index2Pup, Index2Pdn, Index3Dup, Index3Ddn, Fpd, Gpd, &OppUpd);
  
  CreateOperatorCF(NFermion, 2, 2, Index3Dup, Index3Ddn, Index3Dup, Index3Ddn, D4hAkmCF, &OppCF);
  
  //star, easier so try first
  if(star)
  {
    OppNv=malloc(NBathv*sizeof(OperatorType));
    OppHybv=malloc(NBathv*sizeof(OperatorType));

    for(j=0;j<NBathv;j++)
    {
      for(i=0;i<10;i++)
      {
        IndexD[i]=6+10+10*j+i;
      }
      CreateOperatorNumberShell(NFermion, 10, IndexD, &OppNv[j]);
      for(i=0;i<5;i++)
      {
        IndexLDvcdn[i]=6+10+10*(j)+2*i;
        IndexLDvcup[i]=6+10+10*(j)+2*i+1;
      }
      CreateOperatorCF(NFermion, 2, 2, IndexLDvcup, IndexLDvcdn, Index3Dup, Index3Ddn, D4hAkmHop, &OppHybv[j]);
      CreateOperatorCF(NFermion, 2, 2, Index3Dup, Index3Ddn, IndexLDvcup, IndexLDvcdn, D4hAkmHop, &Opptemp);
      OperatorAdd(&OppHybv[j], Opptemp, 1.0, 0.0);
      FreeOperator(&Opptemp);
    } 
  }
  else
  {
    //chain, not done
    /*
    CreateOperatorCF(NFermion, 2, 2, IndexLDup, IndexLDdn, Index3Dup, Index3Ddn, D4hAkmHop, &OppHop);
    CreateOperatorCF(NFermion, 2, 2, Index3Dup, Index3Ddn, IndexLDup, IndexLDdn, D4hAkmHop, &Opptemp);
    OperatorAdd(&OppHop, Opptemp, 1.0, 0.0); FreeOperator(&Opptemp);
    for(j=1;j<NBathv;j++)
    {
      for(i=0;i<10;i++)
      {
        IndexD[i]=6+10+10*j+i;
      }
      CreateOperatorNumberShell(NFermion, 14, IndexD, &OppNv[j]);
      for(i=0;i<5;i++)
      {
        IndexLDvcdn[i]=6+10+10*(j-1)+2*i;
        IndexLDvcup[i]=6+10+10*(j-1)+2*i+1;
        IndexLDvcndn[i]=6+10+10*j+2*i;
        IndexLDvcnup[i]=6+10+10*j+2*i+1;
      }
      CreateOperatorCF(NFermion, 2, 2, IndexLDvcup, IndexLDvcdn, IndexLDvcnup, IndexLDvcndn, SphericalAkm, &OppHybv[j-1]);
      CreateOperatorCF(NFermion, 2, 2, IndexLDvcnup, IndexLDvcndn, IndexLDvcup, IndexLDvcdn, SphericalAkm, &Opptemp);
      OperatorAdd(&OppHybv[j-1], Opptemp, 1.0, 0.0);
      FreeOperator(&Opptemp);
    
    }*/
  }
  
  
  CreateOperatorNumberShell(NFermion, 10, Index3D, &OppNd);
  //CreateOperatorNumberShell(NFermion, 10, IndexLD, &OppNL);
  CreateOperatorCF(NFermion, 2, 2, Index3Dup, Index3Ddn, Index3Dup, Index3Ddn, D4hAkma1g, &OppNda1g);
  CreateOperatorCF(NFermion, 2, 2, Index3Dup, Index3Ddn, Index3Dup, Index3Ddn, D4hAkmb1g, &OppNdb1g);
  CreateOperatorCF(NFermion, 2, 2, Index3Dup, Index3Ddn, Index3Dup, Index3Ddn, D4hAkmb2g, &OppNdb2g);
  CreateOperatorCF(NFermion, 2, 2, Index3Dup, Index3Ddn, Index3Dup, Index3Ddn, D4hAkmeg,  &OppNdeg );
  OperatorCopy(&OppNthreehalf, Oppldots3d); /* <ls> = -3/2n(j=3/2)+n(j=5/2) */
  OperatorAdd(&OppNthreehalf, OppNd, -1.0, 0.0);
  OperatorScale(&OppNthreehalf, -0.4, 0.0);
  OperatorCopy(&OppNfivehalf, OppNd);
  OperatorAdd(&OppNfivehalf, OppNthreehalf, -1.0, 0.0);
  
  /* sum these operators to create the Hamiltonian */
  OperatorCopy(&OppH, OppUdd);
  OperatorAdd(&OppH, OppUpd, 1.0, 0.0);
  OperatorAdd(&OppH, Oppldots2p, zeta2p, 0.0);
  OperatorAdd(&OppH, Oppldots3d, zeta3d, 0.0);
  OperatorAdd(&OppH, OppCF, 1.0, 0.0);
  
  for(j = 0; j < NBathv; j++)
  {
    OperatorAdd(&OppH, OppNv[j],-2.0 + (double)j/(NBathv)*4.0,0.0);
    OperatorAdd(&OppH, OppHybv[j], 1.0, 0.0);
  }
  OperatorAdd(&OppH, OppNd, ed, 0.0);
  OperatorAdd(&OppH, OppSx, 2*Bmag[0]+Hex[0], 0.0);
  OperatorAdd(&OppH, OppSy, 2*Bmag[1]+Hex[1], 0.0);
  OperatorAdd(&OppH, OppSz, 2*Bmag[2]+Hex[2], 0.0);
  OperatorAdd(&OppH, OppLx, Bmag[0], 0.0);
  OperatorAdd(&OppH, OppLy, Bmag[1], 0.0);
  OperatorAdd(&OppH, OppLz, Bmag[2], 0.0);
  
  OperatorChop(&OppH, DBL_EPSILON*1000);
    
  /* create the transition operator for XAS */
  
  CreateOperatorCF(NFermion, 2, 1, Index3Dup, Index3Ddn, Index2Pup, Index2Pdn, Vx, &TXASx);
  CreateOperatorCF(NFermion, 2, 1, Index3Dup, Index3Ddn, Index2Pup, Index2Pdn, Vy, &TXASy);
  CreateOperatorCF(NFermion, 2, 1, Index3Dup, Index3Ddn, Index2Pup, Index2Pdn, Vz, &TXASz);
  CreateOperatorCF(NFermion, 2, 1, Index3Dup, Index3Ddn, Index2Pup, Index2Pdn, Vleft, &TXASl);
  CreateOperatorCF(NFermion, 2, 1, Index3Dup, Index3Ddn, Index2Pup, Index2Pdn, Vright, &TXASr);

  /*Create decay operators for RIXS*/
  for(i = 0; i < 3; i++)
  {
    InitOperator(&TRIXS[i]);
  }
  OperatorConjugateTranspose(TXASx,&TRIXS[0]);
  OperatorConjugateTranspose(TXASy,&TRIXS[1]);
  OperatorConjugateTranspose(TXASz,&TRIXS[2]);

  EndTiming("define operators", &Time);
  
  
  printf("\n========================================================");
  printf("\n==========         find groundstate          ===========");
  printf("\n========================================================\n");
  BeginTiming("find groundstate", &Time);
  
  NBitsKey=16;
 
  
  
  CreateString[6] = 6;
  CreateString[7] = 7;
  CreateString[8] = 8;
  CreateString[9] = 9;
  CreateString[10] = 10;
  CreateString[11] = 11;
  CreateString[12] = 12;
  CreateString[13] = 14;
  CreateString[14] = 15;
  CreateString[15] = 17;
  
  
  printf("\n\n\n\n");
  printf("Wavefunction string: ");
  for(i = 0; i < Nelectrons; i++)
  {
    printf("%d ", CreateString[i]);
  }
  printf("\n\n\n\n");

  

  CreateWaveFunctionFromString(&psi, CreateString, Nelectrons, NFermion, 0, NBitsKey);
  //  WaveFunctionRandomize(&psi);

   epsilon = 1.0e-4; zero = 1.0e-4;
  BeginTiming("LanczosGroundState", &Time);
  LanczosGroundState(OppH, &psi, &epsilon, zero, 0, &psi, 0);

  WaveFunctionRandomize(&psi);

  LanczosGroundState(OppH, &psi, &epsilon, zero, 0, &psi, 0);
  
  WaveFunctionRandomize(&psi);
  
  LanczosGroundState(OppH, &psi, &epsilon, zero, 0, &psi, 0);
  LanczosGroundState(OppH, &psi, &epsilon, zero, 0, &psi, 0);
  LanczosGroundState(OppH, &psi, &epsilon, zero, 0, &psi, 0);
  LanczosGroundState(OppH, &psi, &epsilon, zero, 0, &psi, 0);
  LanczosGroundState(OppH, &psi, &epsilon, zero, 0, &psi, 0);
  LanczosGroundState(OppH, &psi, &epsilon, zero, 0, &psi, 0);
  EndTiming("LanczosGroundState", &Time);
  
  CalculatePrefactorHistogram(Hist, psi, 2.0, 30);
  printf("\nPrefactor Histogram of groundstate\n");
  printf("Interval                                 | percentage of states   | alpha^2 * number of states \n");
  for(i=0;i<30;i++)
  {
    printf("2^-%2i (%10.3E) to 2^-%2i (%10.3E) | %22.15E | %22.15E \n",
           i, pow(0.5,(double)i), i+1, pow(0.5,(double)(i+1)),
           Hist[i], Hist[i]*pow(0.5,(double)(2*i)));
  }
  
  CalculateOccupationHistogram(OccHist, psi, Index3D, 10);
  
  PsiHermitianOperatorPsi(OppH, psi, &energy);
  PsiHermitianOperatorSqrPsi(OppH, psi, &energysqr);
  PsiHermitianOperatorPsi(Oppldots3d, psi, &Esoc);
  Esoc *= zeta3d;
  PsiHermitianOperatorPsi(OppUdd, psi, &ECoulomb);
//  PsiHermitianOperatorPsi(OppHop, psi, &EKin);
  PsiHermitianOperatorPsi(OppSx, psi, &Sx);
  PsiHermitianOperatorPsi(OppSy, psi, &Sy);
  PsiHermitianOperatorPsi(OppSz, psi, &Sz);
  PsiHermitianOperatorPsi(OppSsqr, psi, &Ssqr);
  PsiHermitianOperatorPsi(OppLx, psi, &Lx);
  PsiHermitianOperatorPsi(OppLy, psi, &Ly);
  PsiHermitianOperatorPsi(OppLz, psi, &Lz);
  PsiHermitianOperatorPsi(OppLsqr, psi, &Lsqr);
  PsiHermitianOperatorPsi(OppJx, psi, &Jx);
  PsiHermitianOperatorPsi(OppJy, psi, &Jy);
  PsiHermitianOperatorPsi(OppJz, psi, &Jz);
  PsiHermitianOperatorPsi(OppJsqr, psi, &Jsqr);
  PsiHermitianOperatorPsi(OppNd, psi, &Nd);
  PsiHermitianOperatorPsi(OppNda1g, psi, &Nda1g);
  PsiHermitianOperatorPsi(OppNdb1g, psi, &Ndb1g);
  PsiHermitianOperatorPsi(OppNdb2g, psi, &Ndb2g);
  PsiHermitianOperatorPsi(OppNdeg,  psi, &Ndeg );
  PsiHermitianOperatorPsi(OppNthreehalf, psi, &Nthreehalf);
  PsiHermitianOperatorPsi(OppNfivehalf, psi, &Nfivehalf);
  PsiHermitianOperatorSqrPsi(OppNd, psi, &Ndsqr);
  PsiHermitianOperatorSqrPsi(OppNthreehalf, psi, &Nthreehalfsqr);
  PsiHermitianOperatorSqrPsi(OppNfivehalf, psi, &Nfivehalfsqr);
  EndTiming("find groundstate", &Time);
  PrintTiming(Time);
  
  printf("\n========================================================");
  printf("\n==========   Groundstte expectation values   ===========");
  printf("\n========================================================\n");
  printf("Energy              %22.15E  (variance %22.15E)\n", energy, energysqr-energy*energy);
//  printf("Kinetic energy      %22.15E\n", EKin);
  printf("Relativistic energy %22.15E\n", Esoc);
  printf("Coulomb energy (3d) %22.15E\n", ECoulomb);
  printf("Sx %22.15E Sy %22.15E Sz %22.15E Ssqr %22.15E\n", Sx, Sy, Sz, Ssqr);
  printf("Lx %22.15E Ly %22.15E Lz %22.15E Lsqr %22.15E\n", Lx, Ly, Lz, Lsqr);
  printf("Jx %22.15E Jy %22.15E Jz %22.15E Jsqr %22.15E\n", Jx, Jy, Jz, Jsqr);
  printf("Nd   %22.15E <Nd^2>   %22.15E\n", Nd, Ndsqr);
  printf("Nda1g %22.14E\n", Nda1g);
  printf("Ndb1g %22.14E\n", Ndb1g);
  printf("Ndb2g %22.14E\n", Ndb2g);
  printf("Ndeg  %22.14E\n", Ndeg);
  printf("N3/2 %22.15E <N3/2^2> %22.15E\n", Nthreehalf, Nthreehalfsqr);
  printf("N5/2 %22.15E <N5/2^2> %22.15E\n", Nfivehalf, Nfivehalfsqr);
  printf("\n");
  
  printf("Occupation of d-shell\n");
  for(i=0;i<=10;i++)
  {
    printf("d^%2i | %22.15E\n", i, OccHist[i]);
  }
  printf("\n");
  
  printf("diagonal part of density matrix\n");
  CalculateDensityMatrix(&Density, psi);
  for(i=0;i<Density.Ni;i++)
  {
    printf("%i %22.15E\n",i,Density.VRe[i][i]);
  }
  printf("\n");
  PrintMatrix(Density);
  epsilon = 1.0e-4; zero = 1.0e-4;
  if(doXAS)
  {
    printf("\n========================================================");
    printf("\n==========                XAS                ===========");
    printf("\n========================================================\n");
    BeginTiming("XAS", &Time);

    SpecTri.N = 220;
    InitTriDiagonalMatrix(&SpecTri);
     
    for(i=26;i<31;i++)
    {
      Spectra[i].N=7000;
      Spectra[i].Emin=-85.0;
      Spectra[i].Emax=-30.0;
      Spectra[i].dE = (Spectra[i].Emax - Spectra[i].Emin) / Spectra[i].N;
      Spectra[i].Gamma = 0.4;
      SpectrumAllocate(Spectra[i], &Spectra[i]);
    }
    
    BeginTiming("spectrum", &Time);
    OperatorPsi(TXASx, psi, &Opsi);
    LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[26].Gamma);
    TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[26]);
    FilePrintPlainTriDiagonalMatrix(SpecTri, "LigandField_XASTri_x", "w");
    EndTiming("spectrum", &Time);
    
    BeginTiming("spectrum", &Time);
    OperatorPsi(TXASy, psi, &Opsi);
    LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[27].Gamma);
    TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[27]);
    FilePrintPlainTriDiagonalMatrix(SpecTri, "LigandField_XASTri_y", "w");
    EndTiming("spectrum", &Time);
    
    BeginTiming("spectrum", &Time);
    OperatorPsi(TXASz, psi, &Opsi);
    LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[28].Gamma);
    TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[28]);
    FilePrintPlainTriDiagonalMatrix(SpecTri, "LigandField_XASTri_z", "w");
    EndTiming("spectrum", &Time);
    
    BeginTiming("spectrum", &Time);
    OperatorPsi(TXASl, psi, &Opsi);
    LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[29].Gamma);
    TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[29]);
    EndTiming("spectrum", &Time);
    
    BeginTiming("spectrum", &Time);
    OperatorPsi(TXASr, psi, &Opsi);
    LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[30].Gamma);
    TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[30]);
    EndTiming("spectrum", &Time);
    
    FilePrintSpectra(&Spectra[26], 5, "LigandField_XAS","w");
    EndTiming("XAS", &Time);
    PrintTiming(Time);
  }
  if(doRIXS)
  {
    printf("\n========================================================");
    printf("\n==========               RIXS                ===========");
    printf("\n========================================================\n");
    BeginTiming("RIXS", &Time);
//    for(i = 0; i < 3; i++)
//    {
//      InitOperator(&TRIXS[i]);
//    }
//    OperatorConjugateTranspose(TXASx,&TRIXS[0]);
//    OperatorConjugateTranspose(TXASy,&TRIXS[1]);
//    OperatorConjugateTranspose(TXASz,&TRIXS[2]);

    gamma = 0.1;
    Gfwhm = 0.3;
    gammaRes = 0.4;
    NKrylovXAS = 125;
    NKrylovRIXS = 150;

    for(i = 0; i < numwRes; i++)
      wRes[i] = -80.0 + ((double)i)*0.5;
  
    for(i=0;i<3;i++)
    {
      for(j=0; j < numwRes; j++)
      {
        RIXSSpectra[i][j].N=4601;
        RIXSSpectra[i][j].Emin=-5.0;
        RIXSSpectra[i][j].Emax=18.0;
        RIXSSpectra[i][j].dE = (RIXSSpectra[i][j].Emax - RIXSSpectra[i][j].Emin) / RIXSSpectra[i][j].N;
        RIXSSpectra[i][j].Gamma = gamma;
        SpectrumAllocate(RIXSSpectra[i][j], &RIXSSpectra[i][j]);
      }
    }

    BeginTiming("spectrum", &Time);
//    ResonantSpectra(OppH,TXASx,TRIXS,3,psi,energy,epsilon,&Spectra[26],-77.4,0.4,0.1,70,50);
    ResonantSpectra(OppH,TXASx,TRIXS,3,psi,energy,epsilon,numwRes,RIXSSpectra,wRes,gammaRes,gamma,Gfwhm,NKrylovXAS,NKrylovRIXS);



    EndTiming("spectrum", &Time);
printf("\nPoint1:\n"); fflush(stdout);
    FilePrintRIXSSpectra(numwRes,RIXSSpectra,wRes, 3, "LigandField_RIXS","w");
printf("\nPoint2:\n"); fflush(stdout);
    EndTiming("RIXS", &Time);
    PrintTiming(Time);
  }
  
  EndTiming("Total", &Time);
  PrintTiming(Time);
  return 0;
}