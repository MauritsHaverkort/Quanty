#include "LNOXASFunctions.h"

/*******************************
              IO
 *******************************/

#define NSiteMax 2
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


/*******************************
        Wave-function
 *******************************/

uint32_t InitXASWavefunctionRandom(WaveFunctionType * const psi,
                                   uint32_t const NFermions,
                                   uint32_t const ncore,
                                   uint32_t const nval)
{
  int i, x, y, z, p[NFermions-ncore+2], b[NFermions-ncore], count=0;
  uint8_t *Determinant;
  
  psi->NBitsKey = 16;
  psi->QComplex = 0;
  psi->NBoson = 0;
  psi->NFermion = NFermions;
  if (WaveFunctionInit(psi))
  {
    printf("WaveFunctionInit failed in InitXASWavefunction\n");
    fflush(stdout);
    return 1;
  };
  sprintf(psi->Name, "WaveFunction");
  
  inittwiddle(nval, NFermions-ncore, p);
  
  for(i = 0; i != NFermions-ncore-nval; i++)
  {
    b[i] = 0;
  }
  while(i != NFermions-ncore)
  {
    b[i++] = 1;
  }
  Determinant = (uint8_t*)calloc(psi->NWord, sizeof(uint8_t));
  for (i = 0; i < ncore; i++)
  {
    OperatorCr(Determinant, i, psi->NFermion);
  }
  for (i = ncore; i < NFermions; i++)
  {
    if (b[i-ncore])
    {
      OperatorCr(Determinant, i, psi->NFermion);
    }
  }
  RealWaveFunctionAddElement(psi, Determinant, 1.); free(Determinant);
  count++;
  
  while(!twiddle(&x, &y, &z, p))
  {
    Determinant = (uint8_t*)calloc(psi->NWord, sizeof(uint8_t));
    b[x] = 1;
    b[y] = 0;
    for (i = 0; i < ncore; i++)
    {
      OperatorCr(Determinant, i, psi->NFermion);
    }
    for (i = ncore; i < NFermions; i++)
    {
      if (b[i-ncore])
      {
        OperatorCr(Determinant, i, psi->NFermion);
      }
    }
    RealWaveFunctionAddElement(psi, Determinant, 1.); free(Determinant);
    count++;
  }
  
	printf("The wave function has %i determinants.\n", count);
  WaveFunctionRandomize(psi);
  return 0;
}

uint32_t InitXASWavefunctionRandom2(WaveFunctionType * const psi,
                                    uint32_t const NFermions,
                                    uint32_t const ncore,
                                    uint32_t const nval,
                                    uint32_t const nlig)
{
  int i, x, y, z, p[NFermions-nlig-ncore+2], b[NFermions-nlig-ncore], count=0;
  uint8_t *Determinant;
  
  psi->NBitsKey = 16;
  psi->QComplex = 0;
  psi->NBoson = 0;
  psi->NFermion = NFermions;
  if (WaveFunctionInit(psi))
  {
    printf("WaveFunctionInit failed in InitXASWavefunction\n");
    fflush(stdout);
    return 1;
  };
  sprintf(psi->Name, "WaveFunction");
  
  inittwiddle(nval, NFermions-ncore-nlig, p);
  
  for(i = 0; i != NFermions-ncore-nlig-nval; i++)
  {
    b[i] = 0;
  }
  while(i != NFermions-ncore-nlig)
  {
    b[i++] = 1;
  }
  Determinant = (uint8_t*)calloc(psi->NWord, sizeof(uint8_t));
  for (i = 0; i < ncore; i++)
  {
    OperatorCr(Determinant, i, psi->NFermion);
  }
  for (i = ncore; i < NFermions-nlig; i++)
  {
    if (b[i-ncore])
    {
      OperatorCr(Determinant, i, psi->NFermion);
    }
  }
  for (i = NFermions-1; i >= NFermions-nlig; i--)
  {
      OperatorCr(Determinant, i, psi->NFermion);
  }
  RealWaveFunctionAddElement(psi, Determinant, 1.); free(Determinant);
  count++;
  
  while(!twiddle(&x, &y, &z, p))
  {
    Determinant = (uint8_t*)calloc(psi->NWord, sizeof(uint8_t));
    b[x] = 1;
    b[y] = 0;
    for (i = 0; i < ncore; i++)
    {
      OperatorCr(Determinant, i, psi->NFermion);
    }
    for (i = ncore; i < NFermions-nlig; i++)
    {
      if (b[i-ncore])
      {
        OperatorCr(Determinant, i, psi->NFermion);
      }
    }
    for (i = NFermions-1; i >= NFermions-nlig; i--)
    {
      OperatorCr(Determinant, i, psi->NFermion);
    }
    RealWaveFunctionAddElement(psi, Determinant, 1.); free(Determinant);
    count++;
  }
  
	printf("The wave function has %i determinants.\n", count);
  WaveFunctionRandomize(psi);
  return 0;
}

uint32_t LNOInitXASWavefunctionRandom(WaveFunctionType * const psi,
                                       uint32_t const NFermions,
                                       uint32_t const NBath)
{
  int i, x, y, z, p[12], b[10], count=0;
  uint8_t *Determinant;
  
  psi->NBitsKey = 16;
  psi->QComplex = 0;
  psi->NBoson = 0;
  psi->NFermion = NFermions;
  if (WaveFunctionInit(psi))
  {
    printf("WaveFunctionInit failed in InitXASWavefunction\n");
    fflush(stdout);
    return 1;
  };
  sprintf(psi->Name, "WaveFunction");
  
  inittwiddle(7, 10, p);
  
  for(i = 0; i != 3; i++)
  {
    b[i] = 0;
  }
  while(i != 10)
  {
    b[i++] = 1;
  }
  Determinant = (uint8_t*)calloc(psi->NWord, sizeof(uint8_t));
  for (i = 0; i < 6; i++)
  {
    OperatorCr(Determinant, i, psi->NFermion);
  }
  for (i = 6; i < 16; i++)
  {
    if (b[i-6])
    {
      OperatorCr(Determinant, i, psi->NFermion);
    }
  }
  for (i = NFermions-1-10*NBath; i >= 16; i--)
  {
    OperatorCr(Determinant, i, psi->NFermion);
  }
  RealWaveFunctionAddElement(psi, Determinant, 1.); free(Determinant);
  count++;
  
  while(!twiddle(&x, &y, &z, p))
  {
    Determinant = (uint8_t*)calloc(psi->NWord, sizeof(uint8_t));
    b[x] = 1;
    b[y] = 0;
    for (i = 0; i < 6; i++)
    {
      OperatorCr(Determinant, i, psi->NFermion);
    }
    for (i = 6; i < 16; i++)
    {
      if (b[i-6])
      {
        OperatorCr(Determinant, i, psi->NFermion);
      }
    }
    for (i = NFermions-1-10*NBath; i >= 16; i--)
    {
      OperatorCr(Determinant, i, psi->NFermion);
    }
    RealWaveFunctionAddElement(psi, Determinant, 1.); free(Determinant);
    count++;
  }
  
	printf("The wave function has %i determinants.\n", count);
  WaveFunctionRandomize(psi);
  return 0;
}

/*******************************
          XAS models
 *******************************/

void NiLigandField()
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
  OperatorType OppUdd, OppUpd, OppNd, OppNL;
  OperatorType OppNda1g, OppNdb1g, OppNdb2g, OppNdeg;
  OperatorType Opptemp, OppHop, OppCF;
  OperatorType TPES[10], TIPES[10], TCPES[6], TXASx, TXASy, TXASz, TXASl, TXASr, TCharge;
  double epsilon, zero;
  double energy, energysqr, Esoc, ECoulomb, EKin;
  double Sx, Sy, Sz, Ssqr, Lx, Ly, Lz, Lsqr, Jx, Jy, Jz, Jsqr;
  double Nd, Ndsqr, Nthreehalf, Nfivehalf, Nthreehalfsqr, Nfivehalfsqr;
  double Nda1g, Ndb1g, Ndb2g, Ndeg;
  uint32_t NBitsKey;
  WaveFunctionType psi, Opsi;
  uint16_t CreateString[16];
  TriDiagonalMatrixType SpecTri;
  SpectraType Spectra[35];
  MatrixType Density;
  double Hist[30], OccHist[10];
  uint16_t oadd = 0;
  uint32_t doPES, doIPES, doCPES, doXAS, doXiSpin, doXiCharge;
  uint32_t rotate=1;
  MatrixType Rotated;
  
  printf("\n========================================================");
  printf("\n==========          set parameters           ===========");
  printf("\n========================================================\n");
  
  doPES=0;
  doIPES=0;
  doCPES=0;
  doXAS=1;
  doXiSpin=0;
  doXiCharge=0;
  
  Rydberg = 13.6057;
  
  Bmag[0]= 0.000000; /* magnetic field in the x direction (units of eV) */
  Bmag[1]= 0.000000; /* magnetic field in the y direction (units of eV) */
  Bmag[2]= 0.000000; /* magnetic field in the z direction (units of eV) */
  Hex[0] = 0.000000; /* exchange field in the x direction (units of eV) */
  Hex[1] = 0.000000; /* exchange field in the y direction (units of eV) */
  Hex[2] = 0.000000; /* exchange field in the z direction (units of eV) */
  
  zeta2p = 11.506;
  zeta3d = 0.083;
  
  red2P3D = 0.8;
  red3D3D = 0.8;
  
  Udd = 7.0;
  Upd = 8.5;
  Delta = -0.5;      /* energy differnece between d0 and d1 ligand hole */
  
  Fdd[1] = 12.233 * red3D3D;    /* F2 */
  Fdd[2] = 7.597 * red3D3D;     /* F4 */
  
  Fpd[1] = 7.720 * red2P3D;     /* F2 */
  
  Gpd[0] = 5.783 * red2P3D;     /* G1 */
  Gpd[1] = 3.29 * red2P3D;      /* G3 */
  
  Fdd[0] = Udd + 2.0/63.0 * (Fdd[1] + Fdd[2]);
  Fpd[0] = Upd + 1.0/15.0 * Gpd[0] + 3.0/70.0 * Gpd[1];
  
  /* we define the energy p6d7L10 = 0, p6d8L9 = Dlt */
  /* shift Ep to 0 */
//  ep = 2./23. * (5.*Dlt + 14.*Udd - 68.*Upd);
//  ed = (10.*Dlt - 133.*Udd - 90.*Upd)/23.;
//  el = (-13.*Dlt + 28.*Udd + 48.*Upd)/23.;
  ep = 0.;
  ed = -7. * Udd + 2. * Upd;
  el = - Delta + 8. * Upd;

  
  epsilon = 1.0e-15;
  zero = 1.0e-15;
  SpecTri.N = 50;
  InitTriDiagonalMatrix(&SpecTri);
  
  Ea1g = 1.2;
  Eb1g = 1.2;
  Eb2g = -0.8;
  Eeg  = -0.8;
  
  CFEnergies[0] = Ea1g;
  CFEnergies[1] = Eb1g;
  CFEnergies[2] = Eb2g;
  CFEnergies[3] = Eeg;
  
  PotentialExpandedOnY("D4h", 2, CFEnergies, &D4hAkmCF);
  
  pds = -1.88;
  pdp = pds/2.17;
  Hopa1g = sqrt(3.0) * pds;
  Hopb1g = sqrt(3.0) * pds;
  Hopb2g = 2.0 * pdp;
  Hopeg  = 2.0 * pdp;
  
  HoppingEnergies[0] = Hopa1g;
  HoppingEnergies[1] = Hopb1g;
  HoppingEnergies[2] = Hopb2g;
  HoppingEnergies[3] = Hopeg;
  
  PotentialExpandedOnY("D4h", 2, HoppingEnergies, &D4hAkmHop);
  PotentialExpandedOnY("D4h", 2, HoppingEnergiesa1g, &D4hAkma1g);
  PotentialExpandedOnY("D4h", 2, HoppingEnergiesb1g, &D4hAkmb1g);
  PotentialExpandedOnY("D4h", 2, HoppingEnergiesb2g, &D4hAkmb2g);
  PotentialExpandedOnY("D4h", 2, HoppingEnergieseg,  &D4hAkmeg );
  
  Vx.N=3;             /* for XAS I need polarization dependent transition operators */
  Vx.QComplex=0;      /* these are created from potentials expanded on spherical harmonics */
  Vx.l=Vl;            /* here the expanded potentials are given */
  Vx.m=Vm;
  Vx.Re=VxRe;
  Vy.N=3;
  Vy.QComplex=0;
  Vy.l=Vl;
  Vy.m=Vm;
  Vy.Re=VyRe;
  Vz.N=3;
  Vz.QComplex=0;
  Vz.l=Vl;
  Vz.m=Vm;
  Vz.Re=VzRe;
  Vleft.N=3;
  Vleft.QComplex=0;
  Vleft.l=Vl;
  Vleft.m=Vm;
  Vleft.Re=VlRe;
  Vright.N=3;
  Vright.QComplex=0;
  Vright.l=Vl;
  Vright.m=Vm;
  Vright.Re=VrRe;
  
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
//  OperatorAdd(&OppH, Oppldots3d, zeta3d, 0.0);
  OperatorAdd(&OppH, OppCF, 1.0, 0.0);
  OperatorAdd(&OppH, OppHop, 1.0, 0.0);
  OperatorAdd(&OppH, OppNd, ed, 0.0);
  OperatorAdd(&OppH, OppNL, el, 0.0);
  OperatorAdd(&OppH, OppSx, 2*Bmag[0]+Hex[0], 0.0);
  OperatorAdd(&OppH, OppSy, 2*Bmag[1]+Hex[1], 0.0);
  OperatorAdd(&OppH, OppSz, 2*Bmag[2]+Hex[2], 0.0);
  OperatorAdd(&OppH, OppLx, Bmag[0], 0.0);
  OperatorAdd(&OppH, OppLy, Bmag[1], 0.0);
  OperatorAdd(&OppH, OppLz, Bmag[2], 0.0);
  
  OperatorChop(&OppH, DBL_EPSILON*1000);
  
  /* create the transition operator for pes */
  
  jIm[0]=0;   jIm[1]=0;
  jlist[0]= 7; jlist[1]= 8; jRe[0]=-sqrt(4.0/5.0); jRe[1]= sqrt(1.0/5.0);
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 2, &TPES[0]);
  jlist[0]= 9; jlist[1]=10; jRe[0]=-sqrt(3.0/5.0); jRe[1]= sqrt(2.0/5.0);
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 2, &TPES[1]);
  jlist[0]=11; jlist[1]=12; jRe[0]=-sqrt(2.0/5.0); jRe[1]= sqrt(3.0/5.0);
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 2, &TPES[2]);
  jlist[0]=13; jlist[1]=14; jRe[0]=-sqrt(1.0/5.0); jRe[1]= sqrt(4.0/5.0);
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 2, &TPES[3]);
  jlist[0]=6;             jRe[0]=1.0;
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 1, &TPES[4]);
  jlist[0]= 7; jlist[1]= 8; jRe[0]= sqrt(1.0/5.0); jRe[1]= sqrt(4.0/5.0);
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 2, &TPES[5]);
  jlist[0]= 9; jlist[1]=10; jRe[0]= sqrt(2.0/5.0); jRe[1]= sqrt(3.0/5.0);
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 2, &TPES[6]);
  jlist[0]=11; jlist[1]=12; jRe[0]= sqrt(3.0/5.0); jRe[1]= sqrt(2.0/5.0);
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 2, &TPES[7]);
  jlist[0]=13; jlist[1]=14; jRe[0]= sqrt(4.0/5.0); jRe[1]= sqrt(1.0/5.0);
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 2, &TPES[8]);
  jlist[0]=15;             jRe[0]=1.0;
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 1, &TPES[9]);
  
  /* create the transition operator for ipes */
  
  jIm[0]=0;   jIm[1]=0;
  jlist[0]= 7; jlist[1]= 8; jRe[0]=-sqrt(4.0/5.0); jRe[1]= sqrt(1.0/5.0);
  CreateOperatorCreateList(NFermion, jlist, jRe, jIm, 2, &TIPES[0]);
  jlist[0]= 9; jlist[1]=10; jRe[0]=-sqrt(3.0/5.0); jRe[1]= sqrt(2.0/5.0);
  CreateOperatorCreateList(NFermion, jlist, jRe, jIm, 2, &TIPES[1]);
  jlist[0]=11; jlist[1]=12; jRe[0]=-sqrt(2.0/5.0); jRe[1]= sqrt(3.0/5.0);
  CreateOperatorCreateList(NFermion, jlist, jRe, jIm, 2, &TIPES[2]);
  jlist[0]=13; jlist[1]=14; jRe[0]=-sqrt(1.0/5.0); jRe[1]= sqrt(4.0/5.0);
  CreateOperatorCreateList(NFermion, jlist, jRe, jIm, 2, &TIPES[3]);
  jlist[0]=6;             jRe[0]=1.0;
  CreateOperatorCreateList(NFermion, jlist, jRe, jIm, 1, &TIPES[4]);
  jlist[0]= 7; jlist[1]= 8; jRe[0]= sqrt(1.0/5.0); jRe[1]= sqrt(4.0/5.0);
  CreateOperatorCreateList(NFermion, jlist, jRe, jIm, 2, &TIPES[5]);
  jlist[0]= 9; jlist[1]=10; jRe[0]= sqrt(2.0/5.0); jRe[1]= sqrt(3.0/5.0);
  CreateOperatorCreateList(NFermion, jlist, jRe, jIm, 2, &TIPES[6]);
  jlist[0]=11; jlist[1]=12; jRe[0]= sqrt(3.0/5.0); jRe[1]= sqrt(2.0/5.0);
  CreateOperatorCreateList(NFermion, jlist, jRe, jIm, 2, &TIPES[7]);
  jlist[0]=13; jlist[1]=14; jRe[0]= sqrt(4.0/5.0); jRe[1]= sqrt(1.0/5.0);
  CreateOperatorCreateList(NFermion, jlist, jRe, jIm, 2, &TIPES[8]);
  jlist[0]=15;             jRe[0]=1.0;
  CreateOperatorCreateList(NFermion, jlist, jRe, jIm, 1, &TIPES[9]);
  
  /* create the transition operator for core level pes */
  
  jIm[0]=0;   jIm[1]=0;
  jlist[0]=1; jlist[1]=2; jRe[0]=-sqrt(2.0/3.0); jRe[1]= sqrt(1.0/3.0);
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 2, &TCPES[0]);
  jlist[0]=3; jlist[1]=4; jRe[0]=-sqrt(1.0/3.0); jRe[1]= sqrt(2.0/3.0);
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 2, &TCPES[1]);
  jlist[0]=0;             jRe[0]=1.0;
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 1, &TCPES[2]);
  jlist[0]=1; jlist[1]=2; jRe[0]= sqrt(1.0/3.0); jRe[1]= sqrt(2.0/3.0);
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 2, &TCPES[3]);
  jlist[0]=3; jlist[1]=4; jRe[0]= sqrt(2.0/3.0); jRe[1]= sqrt(1.0/3.0);
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 2, &TCPES[4]);
  jlist[0]=5;             jRe[0]=1.0;
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 1, &TCPES[5]);
  
  /* create the transition operator for XAS */
  
  CreateOperatorCF(NFermion, 2, 1, Index3Dup, Index3Ddn, Index2Pup, Index2Pdn, Vx, &TXASx);
  CreateOperatorCF(NFermion, 2, 1, Index3Dup, Index3Ddn, Index2Pup, Index2Pdn, Vy, &TXASy);
  CreateOperatorCF(NFermion, 2, 1, Index3Dup, Index3Ddn, Index2Pup, Index2Pdn, Vz, &TXASz);
  CreateOperatorCF(NFermion, 2, 1, Index3Dup, Index3Ddn, Index2Pup, Index2Pdn, Vleft, &TXASl);
  CreateOperatorCF(NFermion, 2, 1, Index3Dup, Index3Ddn, Index2Pup, Index2Pdn, Vright, &TXASr);
  EndTiming("define operators", &Time);
  
  if (rotate)
  {
    printf("\n========================================================");
    printf("\n==== Rotate the operators from lm to real orbital ======");
    printf("\n========================================================\n");
    
    Rotated.QComplex = 1;
    Rotated.Ni = 26;
    Rotated.Nj = 26;
    InitMatrix(&Rotated);
    for (i=0; i<Rotated.Ni; i++)
    {
      for (j=0; j<Rotated.Nj; j++)
      {
        Rotated.VRe[i][j] = 0.0;
        Rotated.VIm[i][j] = 0.0;
      }
    }
//    for (i=0; i<6; i++)
//    {
//      Rotated.VRe[i][i] = 1.0;
//    }
    Rotated.VIm[0][0] =-sqrt(0.5);    Rotated.VRe[0][4] = sqrt(0.5);
    Rotated.VIm[1][1] =-sqrt(0.5);    Rotated.VRe[1][5] = sqrt(0.5);
    
    Rotated.VRe[2][2] = 1.0;
    Rotated.VRe[3][3] = 1.0;
    
    Rotated.VIm[4][0] =-sqrt(0.5);    Rotated.VRe[4][4] =-sqrt(0.5);
    Rotated.VIm[5][1] =-sqrt(0.5);    Rotated.VRe[5][5] =-sqrt(0.5);
    
    for (i=0; i<2; i++)
    {
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
    
    for (i=0; i<10; i++)
    {
      RotateOperatorNormalOrder(&Opptemp, TPES[i], Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
      FreeOperator(&TPES[i]); OperatorCopy(&TPES[i], Opptemp); FreeOperator(&Opptemp);
      RotateOperatorNormalOrder(&Opptemp, TIPES[i], Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
      FreeOperator(&TIPES[i]); OperatorCopy(&TIPES[i], Opptemp); FreeOperator(&Opptemp);
    }
    
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
  }
  
  printf("\n========================================================");
  printf("\n==========         find groundstate          ===========");
  printf("\n========================================================\n");
  BeginTiming("find groundstate", &Time);
  
  NBitsKey=16;
  for (i=0; i<6; i++)
  {
    CreateString[i] = i;
  }
  for (i=6; i<13; i++)
  {
    CreateString[i] = i;
  }
  for (i=13; i<23; i++)
  {
    CreateString[i] = i+3;
  }
//  CreateWaveFunctionFromString(&psi, CreateString, 23, NFermion, 0, NBitsKey);
  
//  InitXASWavefunctionRandom2(&psi, 26, 6, 7, 10);
  LNOInitXASWavefunctionRandom(&psi, NFermion, 0);
  
  BeginTiming("LanczosGroundState", &Time);
  LanczosGroundState(OppH, &psi, &epsilon, zero, 0, &psi, 0);
  LanczosGroundState(OppH, &psi, &epsilon, zero, 0, &psi, 0);
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
  
  if(doPES)
  {
    printf("\n========================================================");
    printf("\n==========                PES                ===========");
    printf("\n========================================================\n");
    BeginTiming("PES", &Time);
    
    for(i=0;i<10;i++)
    {
      Spectra[i].N=6000;
      Spectra[i].Emin=-10.0;
      Spectra[i].Emax= 50.0;
      Spectra[i].dE = (Spectra[i].Emax - Spectra[i].Emin) / Spectra[i].N;
      Spectra[i].Gamma = 0.08;
      SpectrumAllocate(Spectra[i], &Spectra[i]);
    }
    
    for(i=0;i<10;i++)
    {
      BeginTiming("spectrum", &Time);
      OperatorPsi(TPES[i], psi, &Opsi);
      LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[i].Gamma);
      TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[i]);
      EndTiming("spectrum", &Time);
    }
    
    FilePrintSpectra(&Spectra[0], 10, "NiLigandField_PES","w");
    EndTiming("PES", &Time);
    PrintTiming(Time);
  }
  
  if(doIPES)
  {
    printf("\n========================================================");
    printf("\n==========               IPES                ===========");
    printf("\n========================================================\n");
    BeginTiming("IPES", &Time);
    
    for(i=10;i<20;i++)
    {
      Spectra[i].N=6000;
      Spectra[i].Emin=-10.0;
      Spectra[i].Emax= 50.0;
      Spectra[i].dE = (Spectra[i].Emax - Spectra[i].Emin) / Spectra[i].N;
      Spectra[i].Gamma = 0.08;
      SpectrumAllocate(Spectra[i], &Spectra[i]);
    }
    
    for(i=0;i<10;i++)
    {
      BeginTiming("spectrum", &Time);
      OperatorPsi(TIPES[i], psi, &Opsi);
      LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[10+i].Gamma);
      TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[10+i]);
      EndTiming("spectrum", &Time);
    }
    
    FilePrintSpectra(&Spectra[10], 10, "NiLigandField_IPES","w");
    EndTiming("IPES", &Time);
    PrintTiming(Time);
  }
  
  if(doCPES)
  {
    printf("\n========================================================");
    printf("\n==========               CPES                ===========");
    printf("\n========================================================\n");
    BeginTiming("CPES", &Time);
    
    for(i=16;i<22;i++)
    {
      Spectra[i].N=7000;
      Spectra[i].Emin=-40.0;
      Spectra[i].Emax= 30.0;
      Spectra[i].dE = (Spectra[i].Emax - Spectra[i].Emin) / Spectra[i].N;
      Spectra[i].Gamma = 0.08;
      SpectrumAllocate(Spectra[i], &Spectra[i]);
    }
    
    for(i=0;i<6;i++)
    {
      BeginTiming("spectrum", &Time);
      OperatorPsi(TCPES[i], psi, &Opsi);
      LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[20+i].Gamma);
      TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[20+i]);
      EndTiming("spectrum", &Time);
    }
    
    FilePrintSpectra(&Spectra[20], 6, "NiLigandField_CPES","w");
    EndTiming("CPES", &Time);
    PrintTiming(Time);
  }
  
  if(doXAS)
  {
    printf("\n========================================================");
    printf("\n==========                XAS                ===========");
    printf("\n========================================================\n");
    BeginTiming("XAS", &Time);
    
    for(i=26;i<31;i++)
    {
      Spectra[i].N=7000;
      Spectra[i].Emin=-40.0;
      Spectra[i].Emax= 30.0;
      Spectra[i].dE = (Spectra[i].Emax - Spectra[i].Emin) / Spectra[i].N;
      Spectra[i].Gamma = 1.0;
      SpectrumAllocate(Spectra[i], &Spectra[i]);
    }
    
    BeginTiming("spectrum", &Time);
    OperatorPsi(TXASx, psi, &Opsi);
    LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[26].Gamma);
    TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[26]);
    FilePrintPlainTriDiagonalMatrix(SpecTri, "NiLigandField_XASTri_x", "w");
    EndTiming("spectrum", &Time);
    
    BeginTiming("spectrum", &Time);
    OperatorPsi(TXASy, psi, &Opsi);
    LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[27].Gamma);
    TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[27]);
    FilePrintPlainTriDiagonalMatrix(SpecTri, "NiLigandField_XASTri_y", "w");
    EndTiming("spectrum", &Time);
    
    BeginTiming("spectrum", &Time);
    OperatorPsi(TXASz, psi, &Opsi);
    LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[28].Gamma);
    TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[28]);
    FilePrintPlainTriDiagonalMatrix(SpecTri, "NiLigandField_XASTri_z", "w");
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
    
    FilePrintSpectra(&Spectra[26], 5, "NiLigandField_XAS","w");
    EndTiming("XAS", &Time);
    PrintTiming(Time);
  }
  
  if(doXiSpin)
  {
    printf("\n========================================================");
    printf("\n==========               Spin                ===========");
    printf("\n========================================================\n");
    BeginTiming("Spin", &Time);
    
    for(i=31;i<34;i++)
    {
      Spectra[i].N=6000;
      Spectra[i].Emin=-1.0;
      Spectra[i].Emax= 11.0;
      Spectra[i].dE = (Spectra[i].Emax - Spectra[i].Emin) / Spectra[i].N;
      Spectra[i].Gamma = 0.01;
      SpectrumAllocate(Spectra[i], &Spectra[i]);
    }
    
    BeginTiming("spectrum", &Time);
    OperatorPsi(OppSx, psi, &Opsi);
    LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[31].Gamma);
    TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[31]);
    EndTiming("spectrum", &Time);
    
    BeginTiming("spectrum", &Time);
    OperatorPsi(OppSy, psi, &Opsi);
    LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[32].Gamma);
    TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[32]);
    EndTiming("spectrum", &Time);
    
    BeginTiming("spectrum", &Time);
    OperatorPsi(OppSz, psi, &Opsi);
    LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[33].Gamma);
    TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[33]);
    EndTiming("spectrum", &Time);
    
    FilePrintSpectra(&Spectra[31], 3, "NiLigandField_Spin","w");
    
    EndTiming("Spin", &Time);
    PrintTiming(Time);
  }
  
  if(doXiCharge)
  {
    printf("\n========================================================");
    printf("\n==========              charge               ===========");
    printf("\n========================================================\n");
    BeginTiming("Charge", &Time);
    
    Spectra[34].N=6000;
    Spectra[34].Emin=-1.0;
    Spectra[34].Emax= 11.0;
    Spectra[34].dE = (Spectra[34].Emax - Spectra[34].Emin) / Spectra[34].N;
    Spectra[34].Gamma = 0.01;
    SpectrumAllocate(Spectra[34], &Spectra[34]);
    
    OperatorCopy(&TCharge, OppNd);
    OperatorAddLader(&TCharge, &oadd, 0, -1, 0);
    
    BeginTiming("spectrum", &Time);
    OperatorPsi(TCharge, psi, &Opsi);
    LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[34].Gamma);
    TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[34]);
    EndTiming("spectrum", &Time);
    
    FilePrintSpectra(&Spectra[34], 1, "NiLigandField_Charge","w");
    
    EndTiming("Charge", &Time);
    PrintTiming(Time);
  }
  
  
  EndTiming("Total", &Time);
  PrintTiming(Time);
}

uint32_t NiImpurity(SpectraSettings SpecSet)
{
  TimeType Time;
  InitializeTiming(&Time);
  BeginTiming("Total", &Time);
  printf("\n========================================================");
  printf("\n==========      start andersen impurity      ===========");
  printf("\n========================================================\n");
  uint32_t i, j, iloop, maxloop=10, igrd, maxgrd=10;
  double egrd[maxgrd];
  double Bmag[3], Hex[3], zeta2p;
  double Rydberg, red2P3D, red3D3D, Fdd[3], Fpd[2], Gpd[2];
  double Udd, Upd, Delta, ed;
  uint16_t Index2Pdn[3] = {0, 2, 4};
  uint16_t Index2Pup[3] = {1, 3, 5};
  uint16_t Index3Ddn[5] = {6, 8, 10, 12, 14};
  uint16_t Index3Dup[5] = {7, 9, 11, 13, 15};
  uint16_t IndexLDdn[5] = {16, 18, 20, 22, 24};
  uint16_t IndexLDup[5] = {17, 19, 21, 23, 25};
  uint16_t Index3D[10] = {6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
  uint16_t IndexLD[10] = {16, 17, 18, 19, 20, 21, 22, 23, 24, 25};
  uint16_t IndexD[10], IndexLDvcup[5], IndexLDvcdn[5], IndexLDvcnup[5], IndexLDvcndn[5];
  uint16_t jlist[10];
  double jRe[10], jIm[10];
  HarmonicsExpansion D4hAkma1g, D4hAkmb1g, D4hAkmb2g, D4hAkmeg, D4hAkmCF;
  double HoppingEnergiesa1g[4]={1., 0., 0., 0.};
  double HoppingEnergiesb1g[4]={0., 1., 0., 0.};
  double HoppingEnergiesb2g[4]={0., 0., 1., 0.};
  double HoppingEnergieseg[4] ={0., 0., 0., 1.};
  double CFEnergies[4]={1.2, 1.2, -0.8, -0.8};
  HarmonicsExpansion Vx, Vy, Vz, Vleft, Vright;
  int Vl[]={1, 1, 1};
  int Vm[]={-1, 0, 1};
  double VxRe[]={sqrt(0.5), 0.0, sqrt(0.5)};
  double VyRe[]={sqrt(0.5), 0.0, -sqrt(0.5)};
  double VzRe[]={0.0, 1.0, 0.0};
  double VlRe[]={1.0, 0.0, 0.0};
  double VrRe[]={0.0, 0.0, 1.0};
  OperatorType OppH, OppHpart, OppSx, OppSy, OppSz, OppSsqr;
  OperatorType OppLx, OppLy, OppLz, OppLsqr;
  OperatorType OppJx, OppJy, OppJz, OppJsqr;
  OperatorType Oppldots2p, Oppldots3d;
  OperatorType OppNthreehalf, OppNfivehalf;
  OperatorType OppUdd, OppUpd, OppCF;
  OperatorType OppNd, OppNda1g, OppNdb1g, OppNdb2g, OppNdeg;
  OperatorType OppNL, OppNLa1g, OppNLb1g, OppNLb2g, OppNLeg;
  OperatorType OppILa1g, OppILb1g, OppILb2g, OppILeg, Opptemp, OppHyb;
  OperatorType OppHybIva1g, OppHybIvb1g, OppHybIvb2g, OppHybIveg;
  OperatorType OppHybIca1g, OppHybIcb1g, OppHybIcb2g, OppHybIceg;
  OperatorType OppHybLva1g, OppHybLvb1g, OppHybLvb2g, OppHybLveg;
  OperatorType OppHybLca1g, OppHybLcb1g, OppHybLcb2g, OppHybLceg;
  OperatorType *OppHybva1g, *OppHybvb1g, *OppHybvb2g, *OppHybveg;
  OperatorType *OppHybca1g, *OppHybcb1g, *OppHybcb2g, *OppHybceg;
  OperatorType *OppNv, *OppNva1g, *OppNvb1g, *OppNvb2g, *OppNveg;
  OperatorType *OppNc, *OppNca1g, *OppNcb1g, *OppNcb2g, *OppNceg;
  OperatorType TPES[10], TIPES[10], TCPES[6], TXASx, TXASy, TXASz, TXASl, TXASr, TCharge;
  OperatorType *VPerturbation;
  TriDiagonalMatrixType *gPerturbation;
  double epsilon, zero, kbTeff, occupation[SpecSet.NSite], onsiteshift[SpecSet.NSite], initoccu, epsilonoccu=1e-10;
  double energy, energysqr, ECoulomb;
  double EKina1g, EKinb1g, EKinb2g, EKineg;
  double Sx, Sy, Sz, Ssqr, Lx, Ly, Lz, Lsqr, Jx, Jy, Jz, Jsqr;
  double Nd, Ndsqr, Nthreehalf, Nfivehalf, Nthreehalfsqr, Nfivehalfsqr, NL;
  double Nda1g, Ndb1g, Ndb2g, Ndeg, NLa1g, NLb1g, NLb2g, NLeg;
  uint32_t NBitsKey;
  uint16_t *CreateString;
  WaveFunctionType psi, psirot, Opsi;
  TriDiagonalMatrixType SpecTri;
  SpectraType Spectra[35];
  MatrixType Density;
  double Hist[30], OccHist[10];
  uint16_t oadd = 0;
  AndersonMatrixType GBath[SpecSet.NSite];
  BADoubleTriDiagonalMatrixType GBathBA[SpecSet.NSite], GBathBAa1g, GBathBAb1g, GBathBAb2g, GBathBAeg;
  uint32_t NBathv, NBathc, NFermion;
  MatrixType Rotated;
  
  printf("\n========================================================");
  printf("\n==========           Read Gbath              ===========");
  printf("\n========================================================\n");
  
  
  for (i=0; i<SpecSet.NSite; i++)
  {
    if (FileReadAndersonMatrixTypeUndeterminedPlusMu(&GBath[i], SpecSet.BathFileName[i], &GBath[i].mu))
      return 1;
  }
  
  
  printf("\n========================================================");
  printf("\n==========          set parameters           ===========");
  printf("\n========================================================\n");
  
  Rydberg = 13.6057;
  
  Bmag[0]= 0.000000; /* magnetic field in the x direction (units of eV) */
  Bmag[1]= 0.000000; /* magnetic field in the y direction (units of eV) */
  Bmag[2]= 0.000000; /* magnetic field in the z direction (units of eV) */
  Hex[0] = 0.000000; /* exchange field in the x direction (units of eV) */
  Hex[1] = 0.000000; /* exchange field in the y direction (units of eV) */
  Hex[2] = 0.000000; /* exchange field in the z direction (units of eV) */
  
  zeta2p = 11.506;
  
  red2P3D = 0.8;
  red3D3D = 0.8;
  
  Udd = SpecSet.Udd;
  Upd = SpecSet.Upd;
  Delta = SpecSet.Delta;
  
  Fdd[1] = 12.233 * red3D3D;    /* F2 */
  Fdd[2] = 7.597 * red3D3D;     /* F4 */
  
  Fpd[1] = 7.720 * red2P3D;     /* F2 */
  
  Gpd[0] = 5.783 * red2P3D;     /* G1 */
  Gpd[1] = 3.29 * red2P3D;      /* G3 */
  
  /*             average interaction energy             */
  /*  Slater, Quantum Theory of Atomic Structure, p.325 */
  Fdd[0] = Udd + 2.0/63.0 * (Fdd[1] + Fdd[2]);
  Fpd[0] = Upd + 1.0/15.0 * Gpd[0] + 3.0/70.0 * Gpd[1];
  
  /*  If one defines p6d7L10 = 0, p6d8L9 = Dlt          */
  /* then:                                              */
  /*  ep = 2./23. * (5.*Dlt + 14.*Udd - 68.*Upd)        */
  /*  ed = (10.*Dlt - 133.*Udd - 90.*Upd)/23.           */
  /*  el = (-13.*Dlt + 28.*Udd + 48.*Upd)/23.           */
  /*                                                    */
  /*  If one defines the energy p6d8L9 - p6d7L10 = Dlt  */
  /*  and ep = 0:                                       */
  /*  ed = -7. * Udd + 2. * Upd                         */
  /*  el = - Delta + 8. * Upd                           */
  /*                                                    */
  /*  If one defines the energy p6d8L9 - p6d7L10 = Dlt  */
  /*  and el = 0:                                       */
  /*  ep = Delta -8. * Upd                              */
  /*  ed = Delta -7. * Udd - 6. * Upd                   */
  ed = Delta -7. * Udd - 6. * Upd;
  
  NBathv = SpecSet.NBath;
  NBathc = SpecSet.NBath;
  
  zero = SpecSet.zeroLancEnergy;
  epsilon = SpecSet.zeroLancWavefunction;
  
  /* some settings one should not change */
  
  NFermion = 6+10+10+10*(NBathv+NBathc);

  if (NBathc != 0)
  {
    OppHybca1g = (OperatorType*)malloc(NBathc*sizeof(OperatorType));
    OppHybcb1g = (OperatorType*)malloc(NBathc*sizeof(OperatorType));
    OppHybcb2g = (OperatorType*)malloc(NBathc*sizeof(OperatorType));
    OppHybceg  = (OperatorType*)malloc(NBathc*sizeof(OperatorType));
    OppNc = (OperatorType*)malloc(NBathc*sizeof(OperatorType));
    OppNca1g = (OperatorType*)malloc(NBathc*sizeof(OperatorType));
    OppNcb1g = (OperatorType*)malloc(NBathc*sizeof(OperatorType));
    OppNcb2g = (OperatorType*)malloc(NBathc*sizeof(OperatorType));
    OppNceg  = (OperatorType*)malloc(NBathc*sizeof(OperatorType));
  }
  if (NBathv != 0)
  {
    OppHybva1g = (OperatorType*)malloc(NBathv*sizeof(OperatorType));
    OppHybvb1g = (OperatorType*)malloc(NBathv*sizeof(OperatorType));
    OppHybvb2g = (OperatorType*)malloc(NBathv*sizeof(OperatorType));
    OppHybveg  = (OperatorType*)malloc(NBathv*sizeof(OperatorType));
    OppNv = (OperatorType*)malloc(NBathv*sizeof(OperatorType));
    OppNva1g = (OperatorType*)malloc(NBathv*sizeof(OperatorType));
    OppNvb1g = (OperatorType*)malloc(NBathv*sizeof(OperatorType));
    OppNvb2g = (OperatorType*)malloc(NBathv*sizeof(OperatorType));
    OppNveg  = (OperatorType*)malloc(NBathv*sizeof(OperatorType));
  }
  
  /* build some potentials */
  
  PotentialExpandedOnY("D4h", 2, HoppingEnergiesa1g, &D4hAkma1g);
  PotentialExpandedOnY("D4h", 2, HoppingEnergiesb1g, &D4hAkmb1g);
  PotentialExpandedOnY("D4h", 2, HoppingEnergiesb2g, &D4hAkmb2g);
  PotentialExpandedOnY("D4h", 2, HoppingEnergieseg,  &D4hAkmeg );
  PotentialExpandedOnY("D4h", 2, CFEnergies,  &D4hAkmCF);
  
  
  /* for XAS I need polarization dependent transition operators */
  /* these are created from potentials expanded on spherical harmonics */
  /* here the expanded potentials are given */
  Vx.N=3;
  Vx.QComplex=0;
  Vx.l=Vl;
  Vx.m=Vm;
  Vx.Re=VxRe;
  Vy.N=3;
  Vy.QComplex=0;
  Vy.l=Vl;
  Vy.m=Vm;
  Vy.Re=VyRe;
  Vz.N=3;
  Vz.QComplex=0;
  Vz.l=Vl;
  Vz.m=Vm;
  Vz.Re=VzRe;
  Vleft.N=3;
  Vleft.QComplex=0;
  Vleft.l=Vl;
  Vleft.m=Vm;
  Vleft.Re=VlRe;
  Vright.N=3;
  Vright.QComplex=0;
  Vright.l=Vl;
  Vright.m=Vm;
  Vright.Re=VrRe;
  
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
  
  CreateOperatorNumberShell(NFermion, 10, Index3D, &OppNd);
  CreateOperatorCF(NFermion, 2, 2, Index3Dup, Index3Ddn, Index3Dup, Index3Ddn, D4hAkma1g, &OppNda1g);
  CreateOperatorCF(NFermion, 2, 2, Index3Dup, Index3Ddn, Index3Dup, Index3Ddn, D4hAkmb1g, &OppNdb1g);
  CreateOperatorCF(NFermion, 2, 2, Index3Dup, Index3Ddn, Index3Dup, Index3Ddn, D4hAkmb2g, &OppNdb2g);
  CreateOperatorCF(NFermion, 2, 2, Index3Dup, Index3Ddn, Index3Dup, Index3Ddn, D4hAkmeg,  &OppNdeg );
  CreateOperatorCF(NFermion, 2, 2, Index3Dup, Index3Ddn, Index3Dup, Index3Ddn, D4hAkmCF,  &OppCF );
  CreateOperatorNumberShell(NFermion, 10, IndexLD, &OppNL);
  CreateOperatorCF(NFermion, 2, 2, IndexLDup, IndexLDdn, IndexLDup, IndexLDdn, D4hAkma1g, &OppNLa1g);
  CreateOperatorCF(NFermion, 2, 2, IndexLDup, IndexLDdn, IndexLDup, IndexLDdn, D4hAkmb1g, &OppNLb1g);
  CreateOperatorCF(NFermion, 2, 2, IndexLDup, IndexLDdn, IndexLDup, IndexLDdn, D4hAkmb2g, &OppNLb2g);
  CreateOperatorCF(NFermion, 2, 2, IndexLDup, IndexLDdn, IndexLDup, IndexLDdn, D4hAkmeg,  &OppNLeg );
  OperatorCopy(&OppNthreehalf, Oppldots3d); /* <ls> = -3/2n(j=3/2)+n(j=5/2) */
  OperatorAdd(&OppNthreehalf, OppNd, -1.0, 0.0);
  OperatorScale(&OppNthreehalf, -0.4, 0.0);
  OperatorCopy(&OppNfivehalf, OppNd);
  OperatorAdd(&OppNfivehalf, OppNthreehalf, -1.0, 0.0);
  
  CreateOperatorCF(NFermion, 2, 2, IndexLDup, IndexLDdn, Index3Dup, Index3Ddn, D4hAkma1g, &OppILa1g);
  CreateOperatorCF(NFermion, 2, 2, Index3Dup, Index3Ddn, IndexLDup, IndexLDdn, D4hAkma1g, &Opptemp);
  OperatorAdd(&OppILa1g, Opptemp, 1.0, 0.0); FreeOperator(&Opptemp);
  CreateOperatorCF(NFermion, 2, 2, IndexLDup, IndexLDdn, Index3Dup, Index3Ddn, D4hAkmb1g, &OppILb1g);
  CreateOperatorCF(NFermion, 2, 2, Index3Dup, Index3Ddn, IndexLDup, IndexLDdn, D4hAkmb1g, &Opptemp);
  OperatorAdd(&OppILb1g, Opptemp, 1.0, 0.0); FreeOperator(&Opptemp);
  CreateOperatorCF(NFermion, 2, 2, IndexLDup, IndexLDdn, Index3Dup, Index3Ddn, D4hAkmb2g, &OppILb2g);
  CreateOperatorCF(NFermion, 2, 2, Index3Dup, Index3Ddn, IndexLDup, IndexLDdn, D4hAkmb2g, &Opptemp);
  OperatorAdd(&OppILb2g, Opptemp, 1.0, 0.0); FreeOperator(&Opptemp);
  CreateOperatorCF(NFermion, 2, 2, IndexLDup, IndexLDdn, Index3Dup, Index3Ddn, D4hAkmeg, &OppILeg);
  CreateOperatorCF(NFermion, 2, 2, Index3Dup, Index3Ddn, IndexLDup, IndexLDdn, D4hAkmeg, &Opptemp);
  OperatorAdd(&OppILeg, Opptemp, 1.0, 0.0); FreeOperator(&Opptemp);
  
  if (NBathv != 0)
  {
    for (i=0; i<10; i++)
    {
      IndexD[i] = 6+10+10+i;
    }
    CreateOperatorNumberShell(NFermion, 10, IndexD, &OppNv[0]);
    for (i=0; i<5; i++)
    {
      IndexLDvcdn[i] = 6+10+10+2*i;
      IndexLDvcup[i] = 6+10+10+2*i+1;
    }
    CreateOperatorCF(NFermion, 2, 2, IndexLDvcup, IndexLDvcdn, IndexLDvcup, IndexLDvcdn, D4hAkma1g, &OppNva1g[0]);
    CreateOperatorCF(NFermion, 2, 2, IndexLDvcup, IndexLDvcdn, IndexLDvcup, IndexLDvcdn, D4hAkmb1g, &OppNvb1g[0]);
    CreateOperatorCF(NFermion, 2, 2, IndexLDvcup, IndexLDvcdn, IndexLDvcup, IndexLDvcdn, D4hAkmb2g, &OppNvb2g[0]);
    CreateOperatorCF(NFermion, 2, 2, IndexLDvcup, IndexLDvcdn, IndexLDvcup, IndexLDvcdn, D4hAkmeg, &OppNveg[0]);
    
    CreateOperatorCF(NFermion, 2, 2, IndexLDvcup, IndexLDvcdn, Index3Dup, Index3Ddn, D4hAkma1g, &OppHybIva1g);
    CreateOperatorCF(NFermion, 2, 2, Index3Dup, Index3Ddn, IndexLDvcup, IndexLDvcdn, D4hAkma1g, &Opptemp);
    OperatorAdd(&OppHybIva1g, Opptemp, 1.0, 0.0); FreeOperator(&Opptemp);
    CreateOperatorCF(NFermion, 2, 2, IndexLDvcup, IndexLDvcdn, Index3Dup, Index3Ddn, D4hAkmb1g, &OppHybIvb1g);
    CreateOperatorCF(NFermion, 2, 2, Index3Dup, Index3Ddn, IndexLDvcup, IndexLDvcdn, D4hAkmb1g, &Opptemp);
    OperatorAdd(&OppHybIvb1g, Opptemp, 1.0, 0.0); FreeOperator(&Opptemp);
    CreateOperatorCF(NFermion, 2, 2, IndexLDvcup, IndexLDvcdn, Index3Dup, Index3Ddn, D4hAkmb2g, &OppHybIvb2g);
    CreateOperatorCF(NFermion, 2, 2, Index3Dup, Index3Ddn, IndexLDvcup, IndexLDvcdn, D4hAkmb2g, &Opptemp);
    OperatorAdd(&OppHybIvb2g, Opptemp, 1.0, 0.0); FreeOperator(&Opptemp);
    CreateOperatorCF(NFermion, 2, 2, IndexLDvcup, IndexLDvcdn, Index3Dup, Index3Ddn, D4hAkmeg, &OppHybIveg);
    CreateOperatorCF(NFermion, 2, 2, Index3Dup, Index3Ddn, IndexLDvcup, IndexLDvcdn, D4hAkmeg, &Opptemp);
    OperatorAdd(&OppHybIveg,  Opptemp, 1.0, 0.0); FreeOperator(&Opptemp);
    
    CreateOperatorCF(NFermion, 2, 2, IndexLDvcup, IndexLDvcdn, IndexLDup, IndexLDdn, D4hAkma1g, &OppHybLva1g);
    CreateOperatorCF(NFermion, 2, 2, IndexLDup, IndexLDdn, IndexLDvcup, IndexLDvcdn, D4hAkma1g, &Opptemp);
    OperatorAdd(&OppHybLva1g, Opptemp, 1.0, 0.0); FreeOperator(&Opptemp);
    CreateOperatorCF(NFermion, 2, 2, IndexLDvcup, IndexLDvcdn, IndexLDup, IndexLDdn, D4hAkmb1g, &OppHybLvb1g);
    CreateOperatorCF(NFermion, 2, 2, IndexLDup, IndexLDdn, IndexLDvcup, IndexLDvcdn, D4hAkmb1g, &Opptemp);
    OperatorAdd(&OppHybLvb1g, Opptemp, 1.0, 0.0); FreeOperator(&Opptemp);
    CreateOperatorCF(NFermion, 2, 2, IndexLDvcup, IndexLDvcdn, IndexLDup, IndexLDdn, D4hAkmb2g, &OppHybLvb2g);
    CreateOperatorCF(NFermion, 2, 2, IndexLDup, IndexLDdn, IndexLDvcup, IndexLDvcdn, D4hAkmb2g, &Opptemp);
    OperatorAdd(&OppHybLvb2g, Opptemp, 1.0, 0.0); FreeOperator(&Opptemp);
    CreateOperatorCF(NFermion, 2, 2, IndexLDvcup, IndexLDvcdn, IndexLDup, IndexLDdn, D4hAkmeg, &OppHybLveg);
    CreateOperatorCF(NFermion, 2, 2, IndexLDup, IndexLDdn, IndexLDvcup, IndexLDvcdn, D4hAkmeg, &Opptemp);
    OperatorAdd(&OppHybLveg,  Opptemp, 1.0, 0.0); FreeOperator(&Opptemp);
    
    for (j=1; j<NBathv; j++)
    {
      for (i=0; i<10; i++)
      {
        IndexD[i]=6+10+10+10*j+i;
      }
      CreateOperatorNumberShell(NFermion, 10, IndexD, &OppNv[j]);
      for (i=0; i<5; i++)
      {
        IndexLDvcdn[i] = 6+10+10+10*(j-1)+2*i;
        IndexLDvcup[i] = 6+10+10+10*(j-1)+2*i+1;
        IndexLDvcndn[i] = 6+10+10+10*j+2*i;
        IndexLDvcnup[i] = 6+10+10+10*j+2*i+1;
      }
      CreateOperatorCF(NFermion, 2, 2, IndexLDvcnup, IndexLDvcndn, IndexLDvcnup, IndexLDvcndn, D4hAkma1g, &OppNva1g[j]);
      CreateOperatorCF(NFermion, 2, 2, IndexLDvcnup, IndexLDvcndn, IndexLDvcnup, IndexLDvcndn, D4hAkmb1g, &OppNvb1g[j]);
      CreateOperatorCF(NFermion, 2, 2, IndexLDvcnup, IndexLDvcndn, IndexLDvcnup, IndexLDvcndn, D4hAkmb2g, &OppNvb2g[j]);
      CreateOperatorCF(NFermion, 2, 2, IndexLDvcnup, IndexLDvcndn, IndexLDvcnup, IndexLDvcndn, D4hAkmeg,  &OppNveg[j] );
      
      CreateOperatorCF(NFermion, 2, 2, IndexLDvcup, IndexLDvcdn, IndexLDvcnup, IndexLDvcndn, D4hAkma1g, &OppHybva1g[j-1]);
      CreateOperatorCF(NFermion, 2, 2, IndexLDvcnup, IndexLDvcndn, IndexLDvcup, IndexLDvcdn, D4hAkma1g, &Opptemp);
      OperatorAdd(&OppHybva1g[j-1], Opptemp, 1.0, 0.0); FreeOperator(&Opptemp);
      CreateOperatorCF(NFermion, 2, 2, IndexLDvcup, IndexLDvcdn, IndexLDvcnup, IndexLDvcndn, D4hAkmb1g, &OppHybvb1g[j-1]);
      CreateOperatorCF(NFermion, 2, 2, IndexLDvcnup, IndexLDvcndn, IndexLDvcup, IndexLDvcdn, D4hAkmb1g, &Opptemp);
      OperatorAdd(&OppHybvb1g[j-1], Opptemp, 1.0, 0.0); FreeOperator(&Opptemp);
      CreateOperatorCF(NFermion, 2, 2, IndexLDvcup, IndexLDvcdn, IndexLDvcnup, IndexLDvcndn, D4hAkmb2g, &OppHybvb2g[j-1]);
      CreateOperatorCF(NFermion, 2, 2, IndexLDvcnup, IndexLDvcndn, IndexLDvcup, IndexLDvcdn, D4hAkmb2g, &Opptemp);
      OperatorAdd(&OppHybvb2g[j-1], Opptemp, 1.0, 0.0); FreeOperator(&Opptemp);
      CreateOperatorCF(NFermion, 2, 2, IndexLDvcup, IndexLDvcdn, IndexLDvcnup, IndexLDvcndn, D4hAkmeg, &OppHybveg[j-1]);
      CreateOperatorCF(NFermion, 2, 2, IndexLDvcnup, IndexLDvcndn, IndexLDvcup, IndexLDvcdn, D4hAkmeg, &Opptemp);
      OperatorAdd(&OppHybveg[j-1], Opptemp, 1.0, 0.0); FreeOperator(&Opptemp);
    }
  }
  
  if (NBathc != 0)
  {
    for (i=0; i<10; i++)
    {
      IndexD[i] = 6+10+10+10*NBathv+i;
    }
    CreateOperatorNumberShell(NFermion, 10, IndexD, &OppNc[0]);
    for (i=0; i<5; i++)
    {
      IndexLDvcdn[i] = 6+10+10+10*NBathv+2*i;
      IndexLDvcup[i] = 6+10+10+10*NBathv+2*i+1;
    }
    CreateOperatorCF(NFermion, 2, 2, IndexLDvcup, IndexLDvcdn, IndexLDvcup, IndexLDvcdn, D4hAkma1g, &OppNca1g[0]);
    CreateOperatorCF(NFermion, 2, 2, IndexLDvcup, IndexLDvcdn, IndexLDvcup, IndexLDvcdn, D4hAkmb1g, &OppNcb1g[0]);
    CreateOperatorCF(NFermion, 2, 2, IndexLDvcup, IndexLDvcdn, IndexLDvcup, IndexLDvcdn, D4hAkmb2g, &OppNcb2g[0]);
    CreateOperatorCF(NFermion, 2, 2, IndexLDvcup, IndexLDvcdn, IndexLDvcup, IndexLDvcdn, D4hAkmeg, &OppNceg[0]);
    
    CreateOperatorCF(NFermion, 2, 2, IndexLDvcup, IndexLDvcdn, Index3Dup, Index3Ddn, D4hAkma1g, &OppHybIca1g);
    CreateOperatorCF(NFermion, 2, 2, Index3Dup, Index3Ddn, IndexLDvcup, IndexLDvcdn, D4hAkma1g, &Opptemp);
    OperatorAdd(&OppHybIca1g, Opptemp, 1.0, 0.0); FreeOperator(&Opptemp);
    CreateOperatorCF(NFermion, 2, 2, IndexLDvcup, IndexLDvcdn, Index3Dup, Index3Ddn, D4hAkmb1g, &OppHybIcb1g);
    CreateOperatorCF(NFermion, 2, 2, Index3Dup, Index3Ddn, IndexLDvcup, IndexLDvcdn, D4hAkmb1g, &Opptemp);
    OperatorAdd(&OppHybIcb1g, Opptemp, 1.0, 0.0); FreeOperator(&Opptemp);
    CreateOperatorCF(NFermion, 2, 2, IndexLDvcup, IndexLDvcdn, Index3Dup, Index3Ddn, D4hAkmb2g, &OppHybIcb2g);
    CreateOperatorCF(NFermion, 2, 2, Index3Dup, Index3Ddn, IndexLDvcup, IndexLDvcdn, D4hAkmb2g, &Opptemp);
    OperatorAdd(&OppHybIcb2g, Opptemp, 1.0, 0.0); FreeOperator(&Opptemp);
    CreateOperatorCF(NFermion, 2, 2, IndexLDvcup, IndexLDvcdn, Index3Dup, Index3Ddn, D4hAkmeg, &OppHybIceg);
    CreateOperatorCF(NFermion, 2, 2, Index3Dup, Index3Ddn, IndexLDvcup, IndexLDvcdn, D4hAkmeg, &Opptemp);
    OperatorAdd(&OppHybIceg,  Opptemp, 1.0, 0.0); FreeOperator(&Opptemp);
    
    CreateOperatorCF(NFermion, 2, 2, IndexLDvcup, IndexLDvcdn, IndexLDup, IndexLDdn, D4hAkma1g, &OppHybLca1g);
    CreateOperatorCF(NFermion, 2, 2, IndexLDup, IndexLDdn, IndexLDvcup, IndexLDvcdn, D4hAkma1g, &Opptemp);
    OperatorAdd(&OppHybLca1g, Opptemp, 1.0, 0.0); FreeOperator(&Opptemp);
    CreateOperatorCF(NFermion, 2, 2, IndexLDvcup, IndexLDvcdn, IndexLDup, IndexLDdn, D4hAkmb1g, &OppHybLcb1g);
    CreateOperatorCF(NFermion, 2, 2, IndexLDup, IndexLDdn, IndexLDvcup, IndexLDvcdn, D4hAkmb1g, &Opptemp);
    OperatorAdd(&OppHybLcb1g, Opptemp, 1.0, 0.0); FreeOperator(&Opptemp);
    CreateOperatorCF(NFermion, 2, 2, IndexLDvcup, IndexLDvcdn, IndexLDup, IndexLDdn, D4hAkmb2g, &OppHybLcb2g);
    CreateOperatorCF(NFermion, 2, 2, IndexLDup, IndexLDdn, IndexLDvcup, IndexLDvcdn, D4hAkmb2g, &Opptemp);
    OperatorAdd(&OppHybLcb2g, Opptemp, 1.0, 0.0); FreeOperator(&Opptemp);
    CreateOperatorCF(NFermion, 2, 2, IndexLDvcup, IndexLDvcdn, IndexLDup, IndexLDdn, D4hAkmeg, &OppHybLceg);
    CreateOperatorCF(NFermion, 2, 2, IndexLDup, IndexLDdn, IndexLDvcup, IndexLDvcdn, D4hAkmeg, &Opptemp);
    OperatorAdd(&OppHybLceg,  Opptemp, 1.0, 0.0); FreeOperator(&Opptemp);
    
    for (j=1; j<NBathc; j++)
    {
      for (i=0; i<10; i++)
      {
        IndexD[i] = 6+10+10+10*NBathv+10*j+i;
      }
      CreateOperatorNumberShell(NFermion, 10, IndexD, &OppNc[j]);
      for (i=0; i<5; i++)
      {
        IndexLDvcdn[i] = 6+10+10+10*NBathv+10*(j-1)+2*i;
        IndexLDvcup[i] = 6+10+10+10*NBathv+10*(j-1)+2*i+1;
        IndexLDvcndn[i] = 6+10+10+10*NBathv+10*j+2*i;
        IndexLDvcnup[i] = 6+10+10+10*NBathv+10*j+2*i+1;
      }
      CreateOperatorCF(NFermion, 2, 2, IndexLDvcnup, IndexLDvcndn, IndexLDvcnup, IndexLDvcndn, D4hAkma1g, &OppNca1g[j]);
      CreateOperatorCF(NFermion, 2, 2, IndexLDvcnup, IndexLDvcndn, IndexLDvcnup, IndexLDvcndn, D4hAkmb1g, &OppNcb1g[j]);
      CreateOperatorCF(NFermion, 2, 2, IndexLDvcnup, IndexLDvcndn, IndexLDvcnup, IndexLDvcndn, D4hAkmb2g, &OppNcb2g[j]);
      CreateOperatorCF(NFermion, 2, 2, IndexLDvcnup, IndexLDvcndn, IndexLDvcnup, IndexLDvcndn, D4hAkmeg,  &OppNceg[j] );
      
      CreateOperatorCF(NFermion, 2, 2, IndexLDvcup, IndexLDvcdn, IndexLDvcnup, IndexLDvcndn, D4hAkma1g, &OppHybca1g[j-1]);
      CreateOperatorCF(NFermion, 2, 2, IndexLDvcnup, IndexLDvcndn, IndexLDvcup, IndexLDvcdn, D4hAkma1g, &Opptemp);
      OperatorAdd(&OppHybca1g[j-1], Opptemp, 1.0, 0.0); FreeOperator(&Opptemp);
      CreateOperatorCF(NFermion, 2, 2, IndexLDvcup, IndexLDvcdn, IndexLDvcnup, IndexLDvcndn, D4hAkmb1g, &OppHybcb1g[j-1]);
      CreateOperatorCF(NFermion, 2, 2, IndexLDvcnup, IndexLDvcndn, IndexLDvcup, IndexLDvcdn, D4hAkmb1g, &Opptemp);
      OperatorAdd(&OppHybcb1g[j-1], Opptemp, 1.0, 0.0); FreeOperator(&Opptemp);
      CreateOperatorCF(NFermion, 2, 2, IndexLDvcup, IndexLDvcdn, IndexLDvcnup, IndexLDvcndn, D4hAkmb2g, &OppHybcb2g[j-1]);
      CreateOperatorCF(NFermion, 2, 2, IndexLDvcnup, IndexLDvcndn, IndexLDvcup, IndexLDvcdn, D4hAkmb2g, &Opptemp);
      OperatorAdd(&OppHybcb2g[j-1], Opptemp, 1.0, 0.0); FreeOperator(&Opptemp);
      CreateOperatorCF(NFermion, 2, 2, IndexLDvcup, IndexLDvcdn, IndexLDvcnup, IndexLDvcndn, D4hAkmeg, &OppHybceg[j-1]);
      CreateOperatorCF(NFermion, 2, 2, IndexLDvcnup, IndexLDvcndn, IndexLDvcup, IndexLDvcdn, D4hAkmeg, &Opptemp);
      OperatorAdd(&OppHybceg[j-1],  Opptemp, 1.0, 0.0); FreeOperator(&Opptemp);
    }
  }
  
  /* create the transition operator for pes */
  
  jIm[0]=0;   jIm[1]=0;
  jlist[0]= 7; jlist[1]= 8; jRe[0]=-sqrt(4.0/5.0); jRe[1]= sqrt(1.0/5.0);
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 2, &TPES[0]);
  jlist[0]= 9; jlist[1]=10; jRe[0]=-sqrt(3.0/5.0); jRe[1]= sqrt(2.0/5.0);
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 2, &TPES[1]);
  jlist[0]=11; jlist[1]=12; jRe[0]=-sqrt(2.0/5.0); jRe[1]= sqrt(3.0/5.0);
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 2, &TPES[2]);
  jlist[0]=13; jlist[1]=14; jRe[0]=-sqrt(1.0/5.0); jRe[1]= sqrt(4.0/5.0);
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 2, &TPES[3]);
  jlist[0]=6;             jRe[0]=1.0;
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 1, &TPES[4]);
  jlist[0]= 7; jlist[1]= 8; jRe[0]= sqrt(1.0/5.0); jRe[1]= sqrt(4.0/5.0);
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 2, &TPES[5]);
  jlist[0]= 9; jlist[1]=10; jRe[0]= sqrt(2.0/5.0); jRe[1]= sqrt(3.0/5.0);
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 2, &TPES[6]);
  jlist[0]=11; jlist[1]=12; jRe[0]= sqrt(3.0/5.0); jRe[1]= sqrt(2.0/5.0);
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 2, &TPES[7]);
  jlist[0]=13; jlist[1]=14; jRe[0]= sqrt(4.0/5.0); jRe[1]= sqrt(1.0/5.0);
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 2, &TPES[8]);
  jlist[0]=15;             jRe[0]=1.0;
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 1, &TPES[9]);
  
  /* create the transition operator for ipes */
  
  jIm[0]=0;   jIm[1]=0;
  jlist[0]= 7; jlist[1]= 8; jRe[0]=-sqrt(4.0/5.0); jRe[1]= sqrt(1.0/5.0);
  CreateOperatorCreateList(NFermion, jlist, jRe, jIm, 2, &TIPES[0]);
  jlist[0]= 9; jlist[1]=10; jRe[0]=-sqrt(3.0/5.0); jRe[1]= sqrt(2.0/5.0);
  CreateOperatorCreateList(NFermion, jlist, jRe, jIm, 2, &TIPES[1]);
  jlist[0]=11; jlist[1]=12; jRe[0]=-sqrt(2.0/5.0); jRe[1]= sqrt(3.0/5.0);
  CreateOperatorCreateList(NFermion, jlist, jRe, jIm, 2, &TIPES[2]);
  jlist[0]=13; jlist[1]=14; jRe[0]=-sqrt(1.0/5.0); jRe[1]= sqrt(4.0/5.0);
  CreateOperatorCreateList(NFermion, jlist, jRe, jIm, 2, &TIPES[3]);
  jlist[0]=6;             jRe[0]=1.0;
  CreateOperatorCreateList(NFermion, jlist, jRe, jIm, 1, &TIPES[4]);
  jlist[0]= 7; jlist[1]= 8; jRe[0]= sqrt(1.0/5.0); jRe[1]= sqrt(4.0/5.0);
  CreateOperatorCreateList(NFermion, jlist, jRe, jIm, 2, &TIPES[5]);
  jlist[0]= 9; jlist[1]=10; jRe[0]= sqrt(2.0/5.0); jRe[1]= sqrt(3.0/5.0);
  CreateOperatorCreateList(NFermion, jlist, jRe, jIm, 2, &TIPES[6]);
  jlist[0]=11; jlist[1]=12; jRe[0]= sqrt(3.0/5.0); jRe[1]= sqrt(2.0/5.0);
  CreateOperatorCreateList(NFermion, jlist, jRe, jIm, 2, &TIPES[7]);
  jlist[0]=13; jlist[1]=14; jRe[0]= sqrt(4.0/5.0); jRe[1]= sqrt(1.0/5.0);
  CreateOperatorCreateList(NFermion, jlist, jRe, jIm, 2, &TIPES[8]);
  jlist[0]=15;             jRe[0]=1.0;
  CreateOperatorCreateList(NFermion, jlist, jRe, jIm, 1, &TIPES[9]);
  
  /* create the transition operator for core level pes */
  
  jIm[0]=0;   jIm[1]=0;
  jlist[0]=1; jlist[1]=2; jRe[0]=-sqrt(2.0/3.0); jRe[1]= sqrt(1.0/3.0);
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 2, &TCPES[0]);
  jlist[0]=3; jlist[1]=4; jRe[0]=-sqrt(1.0/3.0); jRe[1]= sqrt(2.0/3.0);
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 2, &TCPES[1]);
  jlist[0]=0;             jRe[0]=1.0;
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 1, &TCPES[2]);
  jlist[0]=1; jlist[1]=2; jRe[0]= sqrt(1.0/3.0); jRe[1]= sqrt(2.0/3.0);
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 2, &TCPES[3]);
  jlist[0]=3; jlist[1]=4; jRe[0]= sqrt(2.0/3.0); jRe[1]= sqrt(1.0/3.0);
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 2, &TCPES[4]);
  jlist[0]=5;             jRe[0]=1.0;
  CreateOperatorAnnihilateList(NFermion, jlist, jRe, jIm, 1, &TCPES[5]);
  
  /* create the transition operator for XAS */
  
  CreateOperatorCF(NFermion, 2, 1, Index3Dup, Index3Ddn, Index2Pup, Index2Pdn, Vx, &TXASx);
  CreateOperatorCF(NFermion, 2, 1, Index3Dup, Index3Ddn, Index2Pup, Index2Pdn, Vy, &TXASy);
  CreateOperatorCF(NFermion, 2, 1, Index3Dup, Index3Ddn, Index2Pup, Index2Pdn, Vz, &TXASz);
  CreateOperatorCF(NFermion, 2, 1, Index3Dup, Index3Ddn, Index2Pup, Index2Pdn, Vleft, &TXASl);
  CreateOperatorCF(NFermion, 2, 1, Index3Dup, Index3Ddn, Index2Pup, Index2Pdn, Vright, &TXASr);
  EndTiming("define operators", &Time);
  
  
  /* Define rotation matrix */
  
  printf("\n========================================================");
  printf("\n====  Define rotation matrix to 'real' orbitals   ======");
  printf("\n========================================================\n");
  
  Rotated.QComplex = 0;
  Rotated.Ni = NFermion;
  Rotated.Nj = NFermion;
  InitMatrix(&Rotated);
  for (i=0; i<Rotated.Ni; i++)
  {
    for (j=0; j<Rotated.Nj; j++)
    {
      Rotated.VRe[i][j] = 0.0;
    }
  }
    for (i=0; i<6; i++)
    {
      Rotated.VRe[i][i] = 1.0;
    }
//  Rotated.VRe[0][0] =-sqrt(0.5);    Rotated.VRe[0][4] = sqrt(0.5);
//  Rotated.VRe[1][1] =-sqrt(0.5);    Rotated.VRe[1][5] = sqrt(0.5);
//
//  Rotated.VRe[2][2] = 1.0;
//  Rotated.VRe[3][3] = 1.0;
//
//  Rotated.VRe[4][0] =-sqrt(0.5);    Rotated.VRe[4][4] =-sqrt(0.5);
//  Rotated.VRe[5][1] =-sqrt(0.5);    Rotated.VRe[5][5] =-sqrt(0.5);
  
  for (i=0; i<NBathv+NBathc+2; i++)
  {
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
  
  printf("\n========================================================");
  printf("\n==== Rotate the operators from lm to real orbital ======");
  printf("\n========================================================\n");
  
  /* Rotate Operators */
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
  RotateOperatorNormalOrder(&Opptemp, OppNLa1g, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&OppNLa1g); OperatorCopy(&OppNLa1g, Opptemp); FreeOperator(&Opptemp);
  RotateOperatorNormalOrder(&Opptemp, OppNLb1g, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&OppNLb1g); OperatorCopy(&OppNLb1g, Opptemp); FreeOperator(&Opptemp);
  RotateOperatorNormalOrder(&Opptemp, OppNLb2g, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&OppNLb2g); OperatorCopy(&OppNLb2g, Opptemp); FreeOperator(&Opptemp);
  RotateOperatorNormalOrder(&Opptemp, OppNLeg, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&OppNLeg); OperatorCopy(&OppNLeg, Opptemp); FreeOperator(&Opptemp);
  
  RotateOperatorNormalOrder(&Opptemp, OppILa1g, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&OppILa1g); OperatorCopy(&OppILa1g, Opptemp); FreeOperator(&Opptemp);
  RotateOperatorNormalOrder(&Opptemp, OppILb1g, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&OppILb1g); OperatorCopy(&OppILb1g, Opptemp); FreeOperator(&Opptemp);
  RotateOperatorNormalOrder(&Opptemp, OppILb2g, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&OppILb2g); OperatorCopy(&OppILb2g, Opptemp); FreeOperator(&Opptemp);
  RotateOperatorNormalOrder(&Opptemp, OppILeg, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
  FreeOperator(&OppILeg); OperatorCopy(&OppILeg, Opptemp); FreeOperator(&Opptemp);
  
  if (NBathv != 0)
  {
    RotateOperatorNormalOrder(&Opptemp, OppHybIva1g, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&OppHybIva1g); OperatorCopy(&OppHybIva1g, Opptemp); FreeOperator(&Opptemp);
    RotateOperatorNormalOrder(&Opptemp, OppHybIvb1g, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&OppHybIvb1g); OperatorCopy(&OppHybIvb1g, Opptemp); FreeOperator(&Opptemp);
    RotateOperatorNormalOrder(&Opptemp, OppHybIvb2g, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&OppHybIvb2g); OperatorCopy(&OppHybIvb2g, Opptemp); FreeOperator(&Opptemp);
    RotateOperatorNormalOrder(&Opptemp, OppHybIveg, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&OppHybIveg); OperatorCopy(&OppHybIveg, Opptemp); FreeOperator(&Opptemp);
    
    RotateOperatorNormalOrder(&Opptemp, OppHybLva1g, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&OppHybLva1g); OperatorCopy(&OppHybLva1g, Opptemp); FreeOperator(&Opptemp);
    RotateOperatorNormalOrder(&Opptemp, OppHybLvb1g, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&OppHybLvb1g); OperatorCopy(&OppHybLvb1g, Opptemp); FreeOperator(&Opptemp);
    RotateOperatorNormalOrder(&Opptemp, OppHybLvb2g, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&OppHybLvb2g); OperatorCopy(&OppHybLvb2g, Opptemp); FreeOperator(&Opptemp);
    RotateOperatorNormalOrder(&Opptemp, OppHybLveg, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&OppHybLveg); OperatorCopy(&OppHybLveg, Opptemp); FreeOperator(&Opptemp);
    
    for (i=0; i<NBathv; i++)
    {
      RotateOperatorNormalOrder(&Opptemp, OppNv[i], Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
      FreeOperator(&OppNv[i]); OperatorCopy(&OppNv[i], Opptemp); FreeOperator(&Opptemp);
      RotateOperatorNormalOrder(&Opptemp, OppNva1g[i], Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
      FreeOperator(&OppNva1g[i]); OperatorCopy(&OppNva1g[i], Opptemp); FreeOperator(&Opptemp);
      RotateOperatorNormalOrder(&Opptemp, OppNvb1g[i], Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
      FreeOperator(&OppNvb1g[i]); OperatorCopy(&OppNvb1g[i], Opptemp); FreeOperator(&Opptemp);
      RotateOperatorNormalOrder(&Opptemp, OppNvb2g[i], Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
      FreeOperator(&OppNvb2g[i]); OperatorCopy(&OppNvb2g[i], Opptemp); FreeOperator(&Opptemp);
      RotateOperatorNormalOrder(&Opptemp, OppNveg[i], Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
      FreeOperator(&OppNveg[i]); OperatorCopy(&OppNveg[i], Opptemp); FreeOperator(&Opptemp);
    }

    for (i=0; i<NBathv-1; i++)
    {
      RotateOperatorNormalOrder(&Opptemp, OppHybva1g[i], Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
      FreeOperator(&OppHybva1g[i]); OperatorCopy(&OppHybva1g[i], Opptemp); FreeOperator(&Opptemp);
      RotateOperatorNormalOrder(&Opptemp, OppHybvb1g[i], Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
      FreeOperator(&OppHybvb1g[i]); OperatorCopy(&OppHybvb1g[i], Opptemp); FreeOperator(&Opptemp);
      RotateOperatorNormalOrder(&Opptemp, OppHybvb2g[i], Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
      FreeOperator(&OppHybvb2g[i]); OperatorCopy(&OppHybvb2g[i], Opptemp); FreeOperator(&Opptemp);
      RotateOperatorNormalOrder(&Opptemp, OppHybveg[i], Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
      FreeOperator(&OppHybveg[i]); OperatorCopy(&OppHybveg[i], Opptemp); FreeOperator(&Opptemp);
    }
  }
  
  if (NBathc != 0)
  {
    RotateOperatorNormalOrder(&Opptemp, OppHybIca1g, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&OppHybIca1g); OperatorCopy(&OppHybIca1g, Opptemp); FreeOperator(&Opptemp);
    RotateOperatorNormalOrder(&Opptemp, OppHybIcb1g, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&OppHybIcb1g); OperatorCopy(&OppHybIcb1g, Opptemp); FreeOperator(&Opptemp);
    RotateOperatorNormalOrder(&Opptemp, OppHybIcb2g, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&OppHybIcb2g); OperatorCopy(&OppHybIcb2g, Opptemp); FreeOperator(&Opptemp);
    RotateOperatorNormalOrder(&Opptemp, OppHybIceg, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&OppHybIceg); OperatorCopy(&OppHybIceg, Opptemp); FreeOperator(&Opptemp);
    
    RotateOperatorNormalOrder(&Opptemp, OppHybLca1g, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&OppHybLca1g); OperatorCopy(&OppHybLca1g, Opptemp); FreeOperator(&Opptemp);
    RotateOperatorNormalOrder(&Opptemp, OppHybLcb1g, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&OppHybLcb1g); OperatorCopy(&OppHybLcb1g, Opptemp); FreeOperator(&Opptemp);
    RotateOperatorNormalOrder(&Opptemp, OppHybLcb2g, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&OppHybLcb2g); OperatorCopy(&OppHybLcb2g, Opptemp); FreeOperator(&Opptemp);
    RotateOperatorNormalOrder(&Opptemp, OppHybLceg, Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&OppHybLceg); OperatorCopy(&OppHybLceg, Opptemp); FreeOperator(&Opptemp);
    
    for (i=0; i<NBathc; i++)
    {
      RotateOperatorNormalOrder(&Opptemp, OppNc[i], Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
      FreeOperator(&OppNc[i]); OperatorCopy(&OppNc[i], Opptemp); FreeOperator(&Opptemp);
      RotateOperatorNormalOrder(&Opptemp, OppNca1g[i], Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
      FreeOperator(&OppNca1g[i]); OperatorCopy(&OppNca1g[i], Opptemp); FreeOperator(&Opptemp);
      RotateOperatorNormalOrder(&Opptemp, OppNcb1g[i], Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
      FreeOperator(&OppNcb1g[i]); OperatorCopy(&OppNcb1g[i], Opptemp); FreeOperator(&Opptemp);
      RotateOperatorNormalOrder(&Opptemp, OppNcb2g[i], Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
      FreeOperator(&OppNcb2g[i]); OperatorCopy(&OppNcb2g[i], Opptemp); FreeOperator(&Opptemp);
      RotateOperatorNormalOrder(&Opptemp, OppNceg[i], Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
      FreeOperator(&OppNceg[i]); OperatorCopy(&OppNceg[i], Opptemp); FreeOperator(&Opptemp);
    }
    
    for (i=0; i<NBathc-1; i++)
    {
      RotateOperatorNormalOrder(&Opptemp, OppHybca1g[i], Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
      FreeOperator(&OppHybca1g[i]); OperatorCopy(&OppHybca1g[i], Opptemp); FreeOperator(&Opptemp);
      RotateOperatorNormalOrder(&Opptemp, OppHybcb1g[i], Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
      FreeOperator(&OppHybcb1g[i]); OperatorCopy(&OppHybcb1g[i], Opptemp); FreeOperator(&Opptemp);
      RotateOperatorNormalOrder(&Opptemp, OppHybcb2g[i], Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
      FreeOperator(&OppHybcb2g[i]); OperatorCopy(&OppHybcb2g[i], Opptemp); FreeOperator(&Opptemp);
      RotateOperatorNormalOrder(&Opptemp, OppHybceg[i], Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
      FreeOperator(&OppHybceg[i]); OperatorCopy(&OppHybceg[i], Opptemp); FreeOperator(&Opptemp);
    }
  }
  
  /* Rotate spectral operators */
  
  for (i=0; i<6; i++)
  {
    RotateOperatorNormalOrder(&Opptemp, TCPES[i], Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&TCPES[i]); OperatorCopy(&TCPES[i], Opptemp); FreeOperator(&Opptemp);
  }
  
  for (i=0; i<10; i++)
  {
    RotateOperatorNormalOrder(&Opptemp, TPES[i], Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&TPES[i]); OperatorCopy(&TPES[i], Opptemp); FreeOperator(&Opptemp);
    RotateOperatorNormalOrder(&Opptemp, TIPES[i], Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&TIPES[i]); OperatorCopy(&TIPES[i], Opptemp); FreeOperator(&Opptemp);
  }
  
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
  
  /*                                          */
  /* summ operators to create the Hamiltonian */
  /*                                          */
  
  OperatorCopy(&OppHpart, OppUdd);
  OperatorAdd(&OppHpart, OppUpd, 1.0, 0.0);
  OperatorAdd(&OppHpart, OppNd, ed, 0.0);
  OperatorAdd(&OppHpart, Oppldots2p, zeta2p, 0.0);
  
  OperatorChop(&OppHpart, DBL_EPSILON*1000);
  
  printf("\n========================================================");
  printf("\n========== occupation dependent Hamiltonian  ===========");
  printf("\n========================================================\n");
  
  printf("guess the rotation of the bath by assuming an occupation of 8/10.\n");
  
  occupation[0] = 0.5;  //a1g/b1g filling
  occupation[1] = 0.99; //b2g/eg  filling
  
  onsiteshift[0] = -Delta;
  onsiteshift[1] = -Delta;
  
  for (i=0; i<SpecSet.NSite; i++)
  {
    CalculateImpurityOccupationAndersonMatrix(GBath[i], &initoccu);
    if (Abs(initoccu-occupation[i])>epsilonoccu)
    {
      DetermineImpurityShiftFromOccupationAndersonMatrix2(GBath[i], occupation[i], &onsiteshift[i]);
      printf("Onsite shift for site %i is %22.15E\n", i, onsiteshift[i]);
    }
  }
  
  kbTeff = SpecSet.Teff;
  for (i=0; i<SpecSet.NSite; i++)
  {
    AndersonMatrixToBADoubleTriDiagonalMatrixTemperature(GBath[i], onsiteshift[i], kbTeff, &GBathBA[i], &occupation[i]);
    printf("The new impurity Hamiltonian for site %i looks like:\n", i);
    PrintBADoubleTriDiagonalMatrix(GBathBA[i]);
  }
  
  GBathBAa1g.val.N = GBathBA[0].val.N;
  GBathBAa1g.con.N = GBathBA[0].con.N;
  GBathBAb1g.val.N = GBathBA[0].val.N;
  GBathBAb1g.con.N = GBathBA[0].con.N;
  GBathBAb2g.val.N = GBathBA[1].val.N;
  GBathBAb2g.con.N = GBathBA[1].con.N;
  GBathBAeg.val.N = GBathBA[1].val.N;
  GBathBAeg.con.N = GBathBA[1].con.N;
  InitBADoubleTriDiagonalMatrix(&GBathBAa1g);
  InitBADoubleTriDiagonalMatrix(&GBathBAb1g);
  InitBADoubleTriDiagonalMatrix(&GBathBAb2g);
  InitBADoubleTriDiagonalMatrix(&GBathBAeg);
  BADoubleTriDiagonalMatrixCopy(&GBathBAa1g, GBathBA[0]);
  BADoubleTriDiagonalMatrixCopy(&GBathBAb1g, GBathBA[0]);
  BADoubleTriDiagonalMatrixCopy(&GBathBAb2g, GBathBA[1]);
  BADoubleTriDiagonalMatrixCopy(&GBathBAeg,  GBathBA[1]);
  
  OperatorCopy(&OppHyb, OppILa1g);
  OperatorScale(&OppHyb, GBathBAa1g.acvBA, 0.0);
  OperatorAdd(&OppHyb, OppILb1g, GBathBAb1g.acvBA, 0.0);
  OperatorAdd(&OppHyb, OppILb2g, GBathBAb2g.acvBA, 0.0);
  OperatorAdd(&OppHyb, OppILeg,  GBathBAeg.acvBA,  0.0);
  
  OperatorAdd(&OppHyb, OppNLa1g,  GBathBAa1g.acvA,  0.0);
  OperatorAdd(&OppHyb, OppNLb1g,  GBathBAb1g.acvA,  0.0);
  OperatorAdd(&OppHyb, OppNLb2g,  GBathBAb2g.acvA,  0.0);
  OperatorAdd(&OppHyb, OppNLeg,  GBathBAeg.acvA,  0.0);

  OperatorAdd(&OppHyb, OppNda1g,  GBathBAa1g.acvB,  0.0);
  OperatorAdd(&OppHyb, OppNdb1g,  GBathBAb1g.acvB,  0.0);
  OperatorAdd(&OppHyb, OppNdb2g,  GBathBAb2g.acvB,  0.0);
  OperatorAdd(&OppHyb, OppNdeg,  GBathBAeg.acvB,  0.0);
  
  printf("\nThe impurity-ligand hybridizations are:\n");
  printf("a1g: %22.15E\n", GBathBAa1g.acvBA);
  printf("b1g: %22.15E\n", GBathBAb1g.acvBA);
  printf("b2g: %22.15E\n", GBathBAb2g.acvBA);
  printf("eg : %22.15E\n", GBathBAeg.acvBA);
  
  /* Bathv */
  if ((NBathv!=0)&&(GBathBAa1g.val.N>1))
  {
    OperatorAdd(&OppHyb, OppHybLva1g, GBathBAa1g.bvA0, 0.0);
    OperatorAdd(&OppHyb, OppHybIva1g, GBathBAa1g.bvB0, 0.0);
    OperatorAdd(&OppHyb, OppNva1g[0], GBathBAa1g.val.a[2], 0.0);
    
    for (j=1; j<Min(NBathv, GBathBAa1g.val.N-1); j++)
    {
      OperatorAdd(&OppHyb, OppHybva1g[j-1], GBathBAa1g.val.b[j+1], 0.0);
      OperatorAdd(&OppHyb, OppNva1g[j],     GBathBAa1g.val.a[j+2], 0.0);
    }
  }
  if ((NBathv!=0)&&(GBathBAb1g.val.N>1))
  {
    OperatorAdd(&OppHyb, OppHybLvb1g, GBathBAb1g.bvA0, 0.0);
    OperatorAdd(&OppHyb, OppHybIvb1g, GBathBAb1g.bvB0, 0.0);
    OperatorAdd(&OppHyb, OppNvb1g[0], GBathBAb1g.val.a[2], 0.0);
    
    for (j=1; j<Min(NBathv, GBathBAb1g.val.N-1); j++)
    {
      OperatorAdd(&OppHyb, OppHybvb1g[j-1], GBathBAb1g.val.b[j+1], 0.0);
      OperatorAdd(&OppHyb, OppNvb1g[j],     GBathBAb1g.val.a[j+2], 0.0);
    }
  }
  if ((NBathv!=0)&&(GBathBAb2g.val.N>1))
  {
    OperatorAdd(&OppHyb, OppHybLvb2g, GBathBAb2g.bvA0, 0.0);
    OperatorAdd(&OppHyb, OppHybIvb2g, GBathBAb2g.bvB0, 0.0);
    OperatorAdd(&OppHyb, OppNvb2g[0], GBathBAb2g.val.a[2], 0.0);
    
    for (j=1; j<Min(NBathv, GBathBAb2g.val.N-1); j++)
    {
      OperatorAdd(&OppHyb, OppHybvb2g[j-1], GBathBAb2g.val.b[j+1], 0.0);
      OperatorAdd(&OppHyb, OppNvb2g[j],     GBathBAb2g.val.a[j+2], 0.0);
    }
  }
  if ((NBathv!=0)&&(GBathBAeg.val.N>1))
  {
    OperatorAdd(&OppHyb, OppHybLveg, GBathBAeg.bvA0, 0.0);
    OperatorAdd(&OppHyb, OppHybIveg, GBathBAeg.bvB0, 0.0);
    OperatorAdd(&OppHyb, OppNveg[0], GBathBAeg.val.a[2], 0.0);
    
    for (j=1; j<Min(NBathv, GBathBAeg.val.N-1); j++)
    {
      OperatorAdd(&OppHyb, OppHybveg[j-1], GBathBAeg.val.b[j+1], 0.0);
      OperatorAdd(&OppHyb, OppNveg[j],     GBathBAeg.val.a[j+2], 0.0);
    }
  }
  
  /* Bathc */
  if ((NBathc!=0)&&(GBathBAa1g.con.N>1))
  {
    OperatorAdd(&OppHyb, OppHybLca1g, GBathBAa1g.bcA0, 0.0);
    OperatorAdd(&OppHyb, OppHybIca1g, GBathBAa1g.bcB0, 0.0);
    OperatorAdd(&OppHyb, OppNca1g[0], GBathBAa1g.con.a[2], 0.0);
    
    for (j=1; j<Min(NBathc, GBathBAa1g.con.N-1); j++)
    {
      OperatorAdd(&OppHyb, OppHybca1g[j-1], GBathBAa1g.con.b[j+1], 0.0);
      OperatorAdd(&OppHyb, OppNca1g[j],     GBathBAa1g.con.a[j+2], 0.0);
    }
  }
  if ((NBathc!=0)&&(GBathBAb1g.con.N>1))
  {
    OperatorAdd(&OppHyb, OppHybLcb1g, GBathBAb1g.bcA0, 0.0);
    OperatorAdd(&OppHyb, OppHybIcb1g, GBathBAb1g.bcB0, 0.0);
    OperatorAdd(&OppHyb, OppNcb1g[0], GBathBAb1g.con.a[2], 0.0);
    
    for (j=1; j<Min(NBathc, GBathBAb1g.con.N-1); j++)
    {
      OperatorAdd(&OppHyb, OppHybcb1g[j-1], GBathBAb1g.con.b[j+1], 0.0);
      OperatorAdd(&OppHyb, OppNcb1g[j],     GBathBAb1g.con.a[j+2], 0.0);
    }
  }
  if ((NBathc!=0)&&(GBathBAb2g.con.N>1))
  {
    OperatorAdd(&OppHyb, OppHybLcb2g, GBathBAb2g.bcA0, 0.0);
    OperatorAdd(&OppHyb, OppHybIcb2g, GBathBAb2g.bcB0, 0.0);
    OperatorAdd(&OppHyb, OppNcb2g[0], GBathBAb2g.con.a[2], 0.0);
    
    for (j=1; j<Min(NBathc, GBathBAb2g.con.N-1); j++)
    {
      OperatorAdd(&OppHyb, OppHybcb2g[j-1], GBathBAb2g.con.b[j+1], 0.0);
      OperatorAdd(&OppHyb, OppNcb2g[j],     GBathBAb2g.con.a[j+2], 0.0);
    }
  }
  if ((NBathc!=0)&&(GBathBAeg.con.N>1))
  {
    OperatorAdd(&OppHyb, OppHybLceg, GBathBAeg.bcA0, 0.0);
    OperatorAdd(&OppHyb, OppHybIceg, GBathBAeg.bcB0, 0.0);
    OperatorAdd(&OppHyb, OppNceg[0], GBathBAeg.con.a[2], 0.0);
    
    for (j=1; j<Min(NBathc, GBathBAeg.con.N-1); j++)
    {
      OperatorAdd(&OppHyb, OppHybceg[j-1], GBathBAeg.con.b[j+1], 0.0);
      OperatorAdd(&OppHyb, OppNceg[j],     GBathBAeg.con.a[j+2], 0.0);
    }
  }
  
  OperatorCopy(&OppH, OppHpart);
  OperatorAdd(&OppH, OppHyb, 1.0, 0.0);

  printf("\n========================================================");
  printf("\n==========    find groundstate 1-st time     ===========");
  printf("\n========================================================\n");
  BeginTiming("find groundstate", &Time);
  
  NBitsKey = SpecSet.NBitsKey;
  CreateString = (uint16_t*)calloc(6+10+10*NBathv, sizeof(uint16_t));
  for (i=0; i<6; i++)
  {
    CreateString[i] = (uint16_t)i;
  }
  for (i=6; i<(6+10*(NBathv+1)); i++)
  {
    CreateString[i] = (uint16_t)i+10;
  }
  CreateWaveFunctionFromString(&psi, CreateString, 6+10*(1+NBathv), NFermion, 0, NBitsKey);

  BeginTiming("LanczosGroundState", &Time);
  if (OppH.QComplex)
  {
    LanczosGroundState(OppH, &psi, &epsilon, zero, 0, &psi, 0);
    PsiHermitianOperatorPsi(OppH, psi, &egrd[0]);
    LanczosGroundState(OppH, &psi, &epsilon, zero, 0, &psi, 0);
    PsiHermitianOperatorPsi(OppH, psi, &egrd[1]);
    igrd = 1;
    while (igrd < maxgrd && Abs(egrd[igrd]-egrd[igrd-1]) > zero)
    {
      igrd++;
      LanczosGroundState(OppH, &psi, &epsilon, zero, 0, &psi, 0);
      PsiHermitianOperatorPsi(OppH, psi, &egrd[igrd]);
    }
    printf("The ground-state energies are:\n");
    for (i=0; i<MIN(igrd+1, maxgrd); i++)
    {
      printf("%i  %22.15E\n", i, egrd[i]);
    }
  }
  else
  {
    LanczosGroundState(OppH, &psi, &epsilon, zero, 0, &psi, 0);
  }
  EndTiming("LanczosGroundState", &Time);
  
  CalculatePrefactorHistogram(Hist, psi, 2.0, 30);
  printf("\nPrefactor Histogram of groundstate\n");
  printf("Interval                                 | percentage of states   | alpha^2 * number of states \n");
  for (i=0; i<30; i++)
  {
    printf("2^-%2i (%10.3E) to 2^-%2i (%10.3E) | %22.15E | %22.15E \n",
           i, pow(0.5,(double)i), i+1, pow(0.5,(double)(i+1)),
           Hist[i], Hist[i]*pow(0.5,(double)(2*i)));
  }
  
  RotateWaveFunction(&psirot, psi, Rotated);
  CalculateOccupationHistogram(OccHist, psirot, Index3D, 10);
  
  PsiHermitianOperatorPsi(OppH, psi, &energy);
  PsiHermitianOperatorSqrPsi(OppH, psi, &energysqr);
  PsiHermitianOperatorPsi(OppUdd, psi, &ECoulomb);
  PsiHermitianOperatorPsi(OppILa1g, psi, &EKina1g);
  PsiHermitianOperatorPsi(OppILb1g, psi, &EKinb1g);
  PsiHermitianOperatorPsi(OppILb2g, psi, &EKinb2g);
  PsiHermitianOperatorPsi(OppILeg, psi, &EKineg);
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
  PsiHermitianOperatorPsi(OppNL, psi, &NL);
  PsiHermitianOperatorPsi(OppNLa1g, psi, &NLa1g);
  PsiHermitianOperatorPsi(OppNLb1g, psi, &NLb1g);
  PsiHermitianOperatorPsi(OppNLb2g, psi, &NLb2g);
  PsiHermitianOperatorPsi(OppNLeg,  psi, &NLeg );
  PsiHermitianOperatorPsi(OppNthreehalf, psi, &Nthreehalf);
  PsiHermitianOperatorPsi(OppNfivehalf, psi, &Nfivehalf);
  PsiHermitianOperatorSqrPsi(OppNd, psi, &Ndsqr);
  PsiHermitianOperatorSqrPsi(OppNthreehalf, psi, &Nthreehalfsqr);
  PsiHermitianOperatorSqrPsi(OppNfivehalf, psi, &Nfivehalfsqr);
  EndTiming("find groundstate", &Time);
  PrintTiming(Time);
  
  printf("\n========================================================");
  printf("\n==========   Groundstate expectation values   ==========");
  printf("\n========================================================\n");
  printf("Energy              %22.15E  (variance %22.15E)\n", energy, energysqr-energy*energy);
  printf("Kinetic energy\n");
  printf("                    %22.15E\n", EKina1g);
  printf("                    %22.15E\n", EKinb1g);
  printf("                    %22.15E\n", EKinb2g);
  printf("                    %22.15E\n", EKineg);
  printf("Coulomb energy (3d) %22.15E\n", ECoulomb);
  printf("Sx %22.15E Sy %22.15E Sz %22.15E Ssqr %22.15E\n", Sx, Sy, Sz, Ssqr);
  printf("Lx %22.15E Ly %22.15E Lz %22.15E Lsqr %22.15E\n", Lx, Ly, Lz, Lsqr);
  printf("Jx %22.15E Jy %22.15E Jz %22.15E Jsqr %22.15E\n", Jx, Jy, Jz, Jsqr);
  printf("Nd    %22.15E <Nd^2>   %22.15E\n", Nd, Ndsqr);
  printf("Nda1g %22.15E\n", Nda1g);
  printf("Ndb1g %22.15E\n", Ndb1g);
  printf("Ndb2g %22.15E\n", Ndb2g);
  printf("Ndeg  %22.15E\n", Ndeg );
  printf("N3/2 %22.15E <N3/2^2> %22.15E\n", Nthreehalf, Nthreehalfsqr);
  printf("N5/2 %22.15E <N5/2^2> %22.15E\n", Nfivehalf, Nfivehalfsqr);
  printf("NL    %22.15E\n", NL);
  printf("NLa1g %22.15E\n", NLa1g);
  printf("NLb1g %22.15E\n", NLb1g);
  printf("NLb2g %22.15E\n", NLb2g);
  printf("NLeg  %22.15E\n", NLeg );
  printf("Total occupation of impurity and direct bath:\n");
  printf("ABa1g %22.15E (diff 2: %22.15E)\n", Nda1g+NLa1g, Nda1g+NLa1g-2);
  printf("ABb1g %22.15E (diff 2: %22.15E)\n", Ndb1g+NLb1g, Ndb1g+NLb1g-2);
  printf("ABb2g %22.15E (diff 2: %22.15E)\n", Ndb2g+NLb2g, Ndb2g+NLb2g-2);
  printf("ABeg  %22.15E (diff 4: %22.15E)\n", Ndeg+NLeg, Ndeg+NLeg-4);
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
    printf("%i %22.15E\n", i, Density.VRe[i][i]);
  }
  printf("\n");
  PrintMatrix(Density);
  
  for (iloop=1; iloop<maxloop; iloop++)
  {
    printf("\n========================================================");
    printf("\n========== occupation dependent Hamiltonian  ===========");
    printf("\n==========    with previous gs occupation    ===========");
    printf("\n========================================================\n");
    
    printf("use the d occupation of the previous many body calculation\n");
    
    occupation[0] = (Nda1g + Ndb1g)/4.0;
    occupation[1] = (Ndb1g + Ndeg)/6.0;
    
    for (i=0; i<SpecSet.NSite; i++)
    {
      CalculateImpurityOccupationAndersonMatrix(GBath[i], &initoccu);
      if (Abs(initoccu-occupation[i])>epsilonoccu)
      {
        DetermineImpurityShiftFromOccupationAndersonMatrix2(GBath[i], occupation[i], &onsiteshift[i]);
        printf("Onsite shift for site %i is %22.15E\n", i, onsiteshift[i]);
      }
    }
    
    kbTeff = SpecSet.Teff;
    for (i=0; i<SpecSet.NSite; i++)
    {
      FreeBADoubleTriDiagonalMatrix(&GBathBA[i]);
      AndersonMatrixToBADoubleTriDiagonalMatrixTemperature(GBath[i], onsiteshift[i], kbTeff, &GBathBA[i], &occupation[i]);
      printf("The new impurity Hamiltonian for site %i looks like:\n", i);
      PrintBADoubleTriDiagonalMatrix(GBathBA[i]);
    }
    
    FreeBADoubleTriDiagonalMatrix(&GBathBAa1g);
    FreeBADoubleTriDiagonalMatrix(&GBathBAb1g);
    FreeBADoubleTriDiagonalMatrix(&GBathBAb2g);
    FreeBADoubleTriDiagonalMatrix(&GBathBAeg);
    
    GBathBAa1g.val.N = GBathBA[0].val.N;
    GBathBAa1g.con.N = GBathBA[0].con.N;
    GBathBAb1g.val.N = GBathBA[0].val.N;
    GBathBAb1g.con.N = GBathBA[0].con.N;
    GBathBAb2g.val.N = GBathBA[1].val.N;
    GBathBAb2g.con.N = GBathBA[1].con.N;
    GBathBAeg.val.N = GBathBA[1].val.N;
    GBathBAeg.con.N = GBathBA[1].con.N;
    InitBADoubleTriDiagonalMatrix(&GBathBAa1g);
    InitBADoubleTriDiagonalMatrix(&GBathBAb1g);
    InitBADoubleTriDiagonalMatrix(&GBathBAb2g);
    InitBADoubleTriDiagonalMatrix(&GBathBAeg);
    BADoubleTriDiagonalMatrixCopy(&GBathBAa1g, GBathBA[0]);
    BADoubleTriDiagonalMatrixCopy(&GBathBAb1g, GBathBA[0]);
    BADoubleTriDiagonalMatrixCopy(&GBathBAb2g, GBathBA[1]);
    BADoubleTriDiagonalMatrixCopy(&GBathBAeg,  GBathBA[1]);
    
    FreeOperator(&OppHyb);
    OperatorCopy(&OppHyb, OppILa1g);
    OperatorScale(&OppHyb, GBathBAa1g.acvBA, 0.0);
    OperatorAdd(&OppHyb, OppILb1g, GBathBAb1g.acvBA, 0.0);
    OperatorAdd(&OppHyb, OppILb2g, GBathBAb2g.acvBA, 0.0);
    OperatorAdd(&OppHyb, OppILeg,  GBathBAeg.acvBA,  0.0);
    
    OperatorAdd(&OppHyb, OppNLa1g,  GBathBAa1g.acvA,  0.0);
    OperatorAdd(&OppHyb, OppNLb1g,  GBathBAb1g.acvA,  0.0);
    OperatorAdd(&OppHyb, OppNLb2g,  GBathBAb2g.acvA,  0.0);
    OperatorAdd(&OppHyb, OppNLeg,  GBathBAeg.acvA,  0.0);
    
    OperatorAdd(&OppHyb, OppNda1g,  GBathBAa1g.acvB,  0.0);
    OperatorAdd(&OppHyb, OppNdb1g,  GBathBAb1g.acvB,  0.0);
    OperatorAdd(&OppHyb, OppNdb2g,  GBathBAb2g.acvB,  0.0);
    OperatorAdd(&OppHyb, OppNdeg,  GBathBAeg.acvB,  0.0);
    
    printf("\nThe impurity-ligand hybridizations are:\n");
    printf("a1g: %22.15E\n", GBathBAa1g.acvBA);
    printf("b1g: %22.15E\n", GBathBAb1g.acvBA);
    printf("b2g: %22.15E\n", GBathBAb2g.acvBA);
    printf("eg : %22.15E\n", GBathBAeg.acvBA);
    
    /* Bathv */
    if ((NBathv!=0)&&(GBathBAa1g.val.N>1))
    {
      OperatorAdd(&OppHyb, OppHybLva1g, GBathBAa1g.bvA0, 0.0);
      OperatorAdd(&OppHyb, OppHybIva1g, GBathBAa1g.bvB0, 0.0);
      OperatorAdd(&OppHyb, OppNva1g[0], GBathBAa1g.val.a[2], 0.0);
      
      for (j=1; j<Min(NBathv, GBathBAa1g.val.N-1); j++)
      {
        OperatorAdd(&OppHyb, OppHybva1g[j-1], GBathBAa1g.val.b[j+1], 0.0);
        OperatorAdd(&OppHyb, OppNva1g[j],     GBathBAa1g.val.a[j+2], 0.0);
      }
    }
    if ((NBathv!=0)&&(GBathBAb1g.val.N>1))
    {
      OperatorAdd(&OppHyb, OppHybLvb1g, GBathBAb1g.bvA0, 0.0);
      OperatorAdd(&OppHyb, OppHybIvb1g, GBathBAb1g.bvB0, 0.0);
      OperatorAdd(&OppHyb, OppNvb1g[0], GBathBAb1g.val.a[2], 0.0);
      
      for (j=1; j<Min(NBathv, GBathBAb1g.val.N-1); j++)
      {
        OperatorAdd(&OppHyb, OppHybvb1g[j-1], GBathBAb1g.val.b[j+1], 0.0);
        OperatorAdd(&OppHyb, OppNvb1g[j],     GBathBAb1g.val.a[j+2], 0.0);
      }
    }
    if ((NBathv!=0)&&(GBathBAb2g.val.N>1))
    {
      OperatorAdd(&OppHyb, OppHybLvb2g, GBathBAb2g.bvA0, 0.0);
      OperatorAdd(&OppHyb, OppHybIvb2g, GBathBAb2g.bvB0, 0.0);
      OperatorAdd(&OppHyb, OppNvb2g[0], GBathBAb2g.val.a[2], 0.0);
      
      for (j=1; j<Min(NBathv, GBathBAb2g.val.N-1); j++)
      {
        OperatorAdd(&OppHyb, OppHybvb2g[j-1], GBathBAb2g.val.b[j+1], 0.0);
        OperatorAdd(&OppHyb, OppNvb2g[j],     GBathBAb2g.val.a[j+2], 0.0);
      }
    }
    if ((NBathv!=0)&&(GBathBAeg.val.N>1))
    {
      OperatorAdd(&OppHyb, OppHybLveg, GBathBAeg.bvA0, 0.0);
      OperatorAdd(&OppHyb, OppHybIveg, GBathBAeg.bvB0, 0.0);
      OperatorAdd(&OppHyb, OppNveg[0], GBathBAeg.val.a[2], 0.0);
      
      for (j=1; j<Min(NBathv, GBathBAeg.val.N-1); j++)
      {
        OperatorAdd(&OppHyb, OppHybveg[j-1], GBathBAeg.val.b[j+1], 0.0);
        OperatorAdd(&OppHyb, OppNveg[j],     GBathBAeg.val.a[j+2], 0.0);
      }
    }
    
    /* Bathc */
    if ((NBathc!=0)&&(GBathBAa1g.con.N>1))
    {
      OperatorAdd(&OppHyb, OppHybLca1g, GBathBAa1g.bcA0, 0.0);
      OperatorAdd(&OppHyb, OppHybIca1g, GBathBAa1g.bcB0, 0.0);
      OperatorAdd(&OppHyb, OppNca1g[0], GBathBAa1g.con.a[2], 0.0);
      
      for (j=1; j<Min(NBathc, GBathBAa1g.con.N-1); j++)
      {
        OperatorAdd(&OppHyb, OppHybca1g[j-1], GBathBAa1g.con.b[j+1], 0.0);
        OperatorAdd(&OppHyb, OppNca1g[j],     GBathBAa1g.con.a[j+2], 0.0);
      }
    }
    if ((NBathc!=0)&&(GBathBAb1g.con.N>1))
    {
      OperatorAdd(&OppHyb, OppHybLcb1g, GBathBAb1g.bcA0, 0.0);
      OperatorAdd(&OppHyb, OppHybIcb1g, GBathBAb1g.bcB0, 0.0);
      OperatorAdd(&OppHyb, OppNcb1g[0], GBathBAb1g.con.a[2], 0.0);
      
      for (j=1; j<Min(NBathc, GBathBAb1g.con.N-1); j++)
      {
        OperatorAdd(&OppHyb, OppHybcb1g[j-1], GBathBAb1g.con.b[j+1], 0.0);
        OperatorAdd(&OppHyb, OppNcb1g[j],     GBathBAb1g.con.a[j+2], 0.0);
      }
    }
    if ((NBathc!=0)&&(GBathBAb2g.con.N>1))
    {
      OperatorAdd(&OppHyb, OppHybLcb2g, GBathBAb2g.bcA0, 0.0);
      OperatorAdd(&OppHyb, OppHybIcb2g, GBathBAb2g.bcB0, 0.0);
      OperatorAdd(&OppHyb, OppNcb2g[0], GBathBAb2g.con.a[2], 0.0);
      
      for (j=1; j<Min(NBathc, GBathBAb2g.con.N-1); j++)
      {
        OperatorAdd(&OppHyb, OppHybcb2g[j-1], GBathBAb2g.con.b[j+1], 0.0);
        OperatorAdd(&OppHyb, OppNcb2g[j],     GBathBAb2g.con.a[j+2], 0.0);
      }
    }
    if ((NBathc!=0)&&(GBathBAeg.con.N>1))
    {
      OperatorAdd(&OppHyb, OppHybLceg, GBathBAeg.bcA0, 0.0);
      OperatorAdd(&OppHyb, OppHybIceg, GBathBAeg.bcB0, 0.0);
      OperatorAdd(&OppHyb, OppNceg[0], GBathBAeg.con.a[2], 0.0);
      
      for (j=1; j<Min(NBathc, GBathBAeg.con.N-1); j++)
      {
        OperatorAdd(&OppHyb, OppHybceg[j-1], GBathBAeg.con.b[j+1], 0.0);
        OperatorAdd(&OppHyb, OppNceg[j],     GBathBAeg.con.a[j+2], 0.0);
      }
    }
    
    FreeOperator(&OppH);
    OperatorCopy(&OppH, OppHpart);
    OperatorAdd(&OppH, OppHyb, 1.0, 0.0);
    
    printf("\n========================================================");
    printf("\n==========    find groundstate %i-th time     ===========", iloop+1);
    printf("\n========================================================\n");
    BeginTiming("find groundstate", &Time);
    
    BeginTiming("LanczosGroundState", &Time);
    if (OppH.QComplex)
    {
      LanczosGroundState(OppH, &psi, &epsilon, zero, 0, &psi, 0);
      PsiHermitianOperatorPsi(OppH, psi, &egrd[0]);
      LanczosGroundState(OppH, &psi, &epsilon, zero, 0, &psi, 0);
      PsiHermitianOperatorPsi(OppH, psi, &egrd[1]);
      igrd = 1;
      while (igrd < maxgrd && Abs(egrd[igrd]-egrd[igrd-1]) > zero)
      {
        igrd++;
        LanczosGroundState(OppH, &psi, &epsilon, zero, 0, &psi, 0);
        PsiHermitianOperatorPsi(OppH, psi, &egrd[igrd]);
      }
      printf("The ground-state energies are:\n");
      for (i=0; i<MIN(igrd+1, maxgrd); i++)
      {
        printf("%i  %22.15E\n", i, egrd[i]);
      }
    }
    else
    {
      LanczosGroundState(OppH, &psi, &epsilon, zero, 0, &psi, 0);
    }
    EndTiming("LanczosGroundState", &Time);
    
    CalculatePrefactorHistogram(Hist, psi, 2.0, 30);
    printf("\nPrefactor Histogram of groundstate\n");
    printf("Interval                                 | percentage of states   | alpha^2 * number of states \n");
    for (i=0; i<30; i++)
    {
      printf("2^-%2i (%10.3E) to 2^-%2i (%10.3E) | %22.15E | %22.15E \n",
             i, pow(0.5,(double)i), i+1, pow(0.5,(double)(i+1)),
             Hist[i], Hist[i]*pow(0.5,(double)(2*i)));
    }
    
    RotateWaveFunction(&psirot, psi, Rotated);
    CalculateOccupationHistogram(OccHist, psirot, Index3D, 10);
    
    PsiHermitianOperatorPsi(OppH, psi, &energy);
    PsiHermitianOperatorSqrPsi(OppH, psi, &energysqr);
    PsiHermitianOperatorPsi(OppUdd, psi, &ECoulomb);
    PsiHermitianOperatorPsi(OppILa1g, psi, &EKina1g);
    PsiHermitianOperatorPsi(OppILb1g, psi, &EKinb1g);
    PsiHermitianOperatorPsi(OppILb2g, psi, &EKinb2g);
    PsiHermitianOperatorPsi(OppILeg, psi, &EKineg);
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
    PsiHermitianOperatorPsi(OppNL, psi, &NL);
    PsiHermitianOperatorPsi(OppNLa1g, psi, &NLa1g);
    PsiHermitianOperatorPsi(OppNLb1g, psi, &NLb1g);
    PsiHermitianOperatorPsi(OppNLb2g, psi, &NLb2g);
    PsiHermitianOperatorPsi(OppNLeg,  psi, &NLeg );
    PsiHermitianOperatorPsi(OppNthreehalf, psi, &Nthreehalf);
    PsiHermitianOperatorPsi(OppNfivehalf, psi, &Nfivehalf);
    PsiHermitianOperatorSqrPsi(OppNd, psi, &Ndsqr);
    PsiHermitianOperatorSqrPsi(OppNthreehalf, psi, &Nthreehalfsqr);
    PsiHermitianOperatorSqrPsi(OppNfivehalf, psi, &Nfivehalfsqr);
    EndTiming("find groundstate", &Time);
    PrintTiming(Time);
    
    printf("\n========================================================");
    printf("\n==========   Groundstate expectation values   ==========");
    printf("\n========================================================\n");
    printf("Energy              %22.15E  (variance %22.15E)\n", energy, energysqr-energy*energy);
    printf("Kinetic energy\n");
    printf("                    %22.15E\n", EKina1g);
    printf("                    %22.15E\n", EKinb1g);
    printf("                    %22.15E\n", EKinb2g);
    printf("                    %22.15E\n", EKineg);
    printf("Coulomb energy (3d) %22.15E\n", ECoulomb);
    printf("Sx %22.15E Sy %22.15E Sz %22.15E Ssqr %22.15E\n", Sx, Sy, Sz, Ssqr);
    printf("Lx %22.15E Ly %22.15E Lz %22.15E Lsqr %22.15E\n", Lx, Ly, Lz, Lsqr);
    printf("Jx %22.15E Jy %22.15E Jz %22.15E Jsqr %22.15E\n", Jx, Jy, Jz, Jsqr);
    printf("Nd   %22.15E <Nd^2>   %22.15E\n", Nd, Ndsqr);
    printf("Nda1g %22.15E\n", Nda1g);
    printf("Ndb1g %22.15E\n", Ndb1g);
    printf("Ndb2g %22.15E\n", Ndb2g);
    printf("Ndeg  %22.15E\n", Ndeg );
    printf("N3/2 %22.15E <N3/2^2> %22.15E\n", Nthreehalf, Nthreehalfsqr);
    printf("N5/2 %22.15E <N5/2^2> %22.15E\n", Nfivehalf, Nfivehalfsqr);
    printf("NL    %22.15E\n", NL);
    printf("NLa1g %22.15E\n", NLa1g);
    printf("NLb1g %22.15E\n", NLb1g);
    printf("NLb2g %22.15E\n", NLb2g);
    printf("NLeg  %22.15E\n", NLeg );
    printf("Total occupation of impurity and direct bath:\n");
    printf("ABa1g %22.15E (diff 2: %22.15E)\n", Nda1g+NLa1g, Nda1g+NLa1g-2);
    printf("ABb1g %22.15E (diff 2: %22.15E)\n", Ndb1g+NLb1g, Ndb1g+NLb1g-2);
    printf("ABb2g %22.15E (diff 2: %22.15E)\n", Ndb2g+NLb2g, Ndb2g+NLb2g-2);
    printf("ABeg  %22.15E (diff 4: %22.15E)\n", Ndeg+NLeg, Ndeg+NLeg-4);
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
      printf("%i %22.15E\n", i, Density.VRe[i][i]);
    }
    printf("\n");
    PrintMatrix(Density);
  }
  
  if (SpecSet.Dyson)
  {
    VPerturbation= (OperatorType*)malloc(20*sizeof(OperatorType));
    gPerturbation= (TriDiagonalMatrixType*)malloc(20*sizeof(TriDiagonalMatrixType));
    SpectraVPerturbation(SpecSet, GBathBAa1g, GBathBAb1g, GBathBAb2g, GBathBAeg, VPerturbation, gPerturbation);
  }
  
  SpecTri.N = SpecSet.Dyson ? SpecSet.NStates : SpecSet.NTriDiag;
  InitTriDiagonalMatrix(&SpecTri);
  
  if (SpecSet.Spectra & 1)
  {
    printf("\n========================================================");
    printf("\n==========                PES                ===========");
    printf("\n========================================================\n");
    BeginTiming("PES", &Time);
    
    for(i=0;i<10;i++)
    {
      Spectra[i].N=6000;
      Spectra[i].Emin=-10.0;
      Spectra[i].Emax= 50.0;
      Spectra[i].dE = (Spectra[i].Emax - Spectra[i].Emin) / Spectra[i].N;
      Spectra[i].Gamma = 0.08;
      SpectrumAllocate(Spectra[i], &Spectra[i]);
    }
    
    switch (SpecSet.Dyson)
    {
      case 1:
        for (i=0; i<10; i++)
        {
          BeginTiming("spectrum", &Time);
          OperatorPsi(TPES[i], psi, &Opsi);
          LanczosTriDiagonalizeDysonLowMem(OppH, Opsi, SpecSet.zeroLancWavefunction,
                                           SpecSet.NTriDiag, VPerturbation,
                                           gPerturbation, 20, &SpecTri);
          TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[i]);
          EndTiming("spectrum", &Time);
        }
        break;
        
      default:
        for (i=0; i<10; i++)
        {
          BeginTiming("spectrum", &Time);
          OperatorPsi(TPES[i], psi, &Opsi);
          LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[i].Gamma);
          TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[i]);
          EndTiming("spectrum", &Time);
        }
        break;
    }
    
    FilePrintSpectra(&Spectra[0], 10, "NiImpurity_PES","w");
    EndTiming("PES", &Time);
    PrintTiming(Time);
  }
  
  if (SpecSet.Spectra & 2)
  {
    printf("\n========================================================");
    printf("\n==========               IPES                ===========");
    printf("\n========================================================\n");
    BeginTiming("IPES", &Time);
    
    for(i=10;i<20;i++)
    {
      Spectra[i].N=6000;
      Spectra[i].Emin=-10.0;
      Spectra[i].Emax= 50.0;
      Spectra[i].dE = (Spectra[i].Emax - Spectra[i].Emin) / Spectra[i].N;
      Spectra[i].Gamma = 0.08;
      SpectrumAllocate(Spectra[i], &Spectra[i]);
    }
    
    switch (SpecSet.Dyson)
    {
      case 1:
        for (i=0; i<10; i++)
        {
          BeginTiming("spectrum", &Time);
          OperatorPsi(TIPES[i], psi, &Opsi);
          LanczosTriDiagonalizeDysonLowMem(OppH, Opsi, SpecSet.zeroLancWavefunction,
                                           SpecSet.NTriDiag, VPerturbation,
                                           gPerturbation, 20, &SpecTri);
          TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[10+i]);
          EndTiming("spectrum", &Time);
        }
        break;
        
      default:
        for (i=0; i<10; i++)
        {
          BeginTiming("spectrum", &Time);
          OperatorPsi(TIPES[i], psi, &Opsi);
          LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[10+i].Gamma);
          TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[10+i]);
          EndTiming("spectrum", &Time);
        }
        break;
    }
    
    FilePrintSpectra(&Spectra[10], 10, "NiImpurity_IPES","w");
    EndTiming("IPES", &Time);
    PrintTiming(Time);
  }
  
  if (SpecSet.Spectra & 4)
  {
    printf("\n========================================================");
    printf("\n==========               CPES                ===========");
    printf("\n========================================================\n");
    BeginTiming("CPES", &Time);
    
    for(i=16;i<22;i++)
    {
      Spectra[i].N=7000;
      Spectra[i].Emin=-40.0;
      Spectra[i].Emax= 30.0;
      Spectra[i].dE = (Spectra[i].Emax - Spectra[i].Emin) / Spectra[i].N;
      Spectra[i].Gamma = 0.08;
      SpectrumAllocate(Spectra[i], &Spectra[i]);
    }
    
    switch (SpecSet.Dyson)
    {
      case 1:
        for (i=0; i<6; i++)
        {
          BeginTiming("spectrum", &Time);
          OperatorPsi(TCPES[i], psi, &Opsi);
          LanczosTriDiagonalizeDysonLowMem(OppH, Opsi, SpecSet.zeroLancWavefunction,
                                           SpecSet.NTriDiag, VPerturbation,
                                           gPerturbation, 20, &SpecTri);
          TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[20+i]);
          EndTiming("spectrum", &Time);
        }
        break;
        
      default:
        for (i=0; i<6; i++)
        {
          BeginTiming("spectrum", &Time);
          OperatorPsi(TCPES[i], psi, &Opsi);
          LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[20+i].Gamma);
          TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[20+i]);
          EndTiming("spectrum", &Time);
        }
        break;
    }
    
    FilePrintSpectra(&Spectra[20], 6, "NiImpurity_CPES","w");
    EndTiming("CPES", &Time);
    PrintTiming(Time);
  }
  
  if (SpecSet.Spectra & 8)
  {
    printf("\n========================================================");
    printf("\n==========                XAS                ===========");
    printf("\n========================================================\n");
    BeginTiming("XAS", &Time);
    
    for(i=26;i<31;i++)
    {
      Spectra[i].N=7000;
      Spectra[i].Emin=-40.0;
      Spectra[i].Emax= 30.0;
      Spectra[i].dE = (Spectra[i].Emax - Spectra[i].Emin) / Spectra[i].N;
      Spectra[i].Gamma = SpecSet.Gamma;
      SpectrumAllocate(Spectra[i], &Spectra[i]);
    }
    
    switch (SpecSet.Dyson)
    {
      case 1:
        BeginTiming("spectrum", &Time);
        OperatorPsi(TXASx, psi, &Opsi);
        LanczosTriDiagonalizeDysonLowMem(OppH, Opsi, SpecSet.zeroLancWavefunction,
                                           SpecSet.NTriDiag, VPerturbation,
                                           gPerturbation, 20, &SpecTri);
        TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[26]);
        FilePrintPlainTriDiagonalMatrix(SpecTri, "NiImpurity_XASTri_x", "w");
        EndTiming("spectrum", &Time);
        
        BeginTiming("spectrum", &Time);
        OperatorPsi(TXASy, psi, &Opsi);
        LanczosTriDiagonalizeDysonLowMem(OppH, Opsi, SpecSet.zeroLancWavefunction,
                                         SpecSet.NTriDiag, VPerturbation,
                                         gPerturbation, 20, &SpecTri);
        TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[27]);
        FilePrintPlainTriDiagonalMatrix(SpecTri, "NiImpurity_XASTri_y", "w");
        EndTiming("spectrum", &Time);
        
        BeginTiming("spectrum", &Time);
        OperatorPsi(TXASz, psi, &Opsi);
        LanczosTriDiagonalizeDysonLowMem(OppH, Opsi, SpecSet.zeroLancWavefunction,
                                         SpecSet.NTriDiag, VPerturbation,
                                         gPerturbation, 20, &SpecTri);
        TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[28]);
        FilePrintPlainTriDiagonalMatrix(SpecTri, "NiImpurity_XASTri_z", "w");
        EndTiming("spectrum", &Time);
        
        BeginTiming("spectrum", &Time);
        OperatorPsi(TXASl, psi, &Opsi);
        LanczosTriDiagonalizeDysonLowMem(OppH, Opsi, SpecSet.zeroLancWavefunction,
                                         SpecSet.NTriDiag, VPerturbation,
                                         gPerturbation, 20, &SpecTri);
        TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[29]);
        FilePrintPlainTriDiagonalMatrix(SpecTri, "NiImpurity_XASTri_l", "w");
        EndTiming("spectrum", &Time);
        
        BeginTiming("spectrum", &Time);
        OperatorPsi(TXASr, psi, &Opsi);
        LanczosTriDiagonalizeDysonLowMem(OppH, Opsi, SpecSet.zeroLancWavefunction,
                                         SpecSet.NTriDiag, VPerturbation,
                                         gPerturbation, 20, &SpecTri);
        TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[30]);
        FilePrintPlainTriDiagonalMatrix(SpecTri, "NiImpurity_XASTri_r", "w");
        EndTiming("spectrum", &Time);
        
        FilePrintSpectra(&Spectra[26], 5, "NiImpurity_XAS","w");
        EndTiming("XAS", &Time);
        PrintTiming(Time);
        
        break;
        
      default:
        BeginTiming("spectrum", &Time);
        OperatorPsi(TXASx, psi, &Opsi);
        LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[26].Gamma);
        TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[26]);
        FilePrintPlainTriDiagonalMatrix(SpecTri, "NiImpurity_XASTri_x", "w");
        EndTiming("spectrum", &Time);
        
        BeginTiming("spectrum", &Time);
        OperatorPsi(TXASy, psi, &Opsi);
        LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[27].Gamma);
        TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[27]);
        FilePrintPlainTriDiagonalMatrix(SpecTri, "NiImpurity_XASTri_y", "w");
        EndTiming("spectrum", &Time);
        
        BeginTiming("spectrum", &Time);
        OperatorPsi(TXASz, psi, &Opsi);
        LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[28].Gamma);
        TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[28]);
        FilePrintPlainTriDiagonalMatrix(SpecTri, "NiImpurity_XASTri_z", "w");
        EndTiming("spectrum", &Time);
        
        BeginTiming("spectrum", &Time);
        OperatorPsi(TXASl, psi, &Opsi);
        LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[29].Gamma);
        TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[29]);
        FilePrintPlainTriDiagonalMatrix(SpecTri, "NiImpurity_XASTri_l", "w");
        EndTiming("spectrum", &Time);
        
        BeginTiming("spectrum", &Time);
        OperatorPsi(TXASr, psi, &Opsi);
        LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[30].Gamma);
        TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[30]);
        FilePrintPlainTriDiagonalMatrix(SpecTri, "NiImpurity_XASTri_r", "w");
        EndTiming("spectrum", &Time);
        
        break;
    }
    
    FilePrintSpectra(&Spectra[26], 5, "NiImpurity_XAS","w");
    EndTiming("XAS", &Time);
    PrintTiming(Time);
  }
  
  if (SpecSet.Spectra & 16)
  {
    printf("\n========================================================");
    printf("\n==========               Spin                ===========");
    printf("\n========================================================\n");
    BeginTiming("Spin", &Time);
    
    for(i=31;i<34;i++)
    {
      Spectra[i].N=6000;
      Spectra[i].Emin=-1.0;
      Spectra[i].Emax= 11.0;
      Spectra[i].dE = (Spectra[i].Emax - Spectra[i].Emin) / Spectra[i].N;
      Spectra[i].Gamma = 0.01;
      SpectrumAllocate(Spectra[i], &Spectra[i]);
    }
    
    switch (SpecSet.Dyson)
    {
      case 1:
        BeginTiming("spectrum", &Time);
        OperatorPsi(OppSx, psi, &Opsi);
        LanczosTriDiagonalizeDysonLowMem(OppH, Opsi, SpecSet.zeroLancWavefunction,
                                         SpecSet.NTriDiag, VPerturbation,
                                         gPerturbation, 20, &SpecTri);
        TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[31]);
        EndTiming("spectrum", &Time);
        
        BeginTiming("spectrum", &Time);
        OperatorPsi(OppSy, psi, &Opsi);
        LanczosTriDiagonalizeDysonLowMem(OppH, Opsi, SpecSet.zeroLancWavefunction,
                                         SpecSet.NTriDiag, VPerturbation,
                                         gPerturbation, 20, &SpecTri);
        TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[32]);
        EndTiming("spectrum", &Time);
        
        BeginTiming("spectrum", &Time);
        OperatorPsi(OppSz, psi, &Opsi);
        LanczosTriDiagonalizeDysonLowMem(OppH, Opsi, SpecSet.zeroLancWavefunction,
                                         SpecSet.NTriDiag, VPerturbation,
                                         gPerturbation, 20, &SpecTri);
        TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[33]);
        EndTiming("spectrum", &Time);
        
        break;
        
      default:
        BeginTiming("spectrum", &Time);
        OperatorPsi(OppSx, psi, &Opsi);
        LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[31].Gamma);
        TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[31]);
        EndTiming("spectrum", &Time);
        
        BeginTiming("spectrum", &Time);
        OperatorPsi(OppSy, psi, &Opsi);
        LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[32].Gamma);
        TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[32]);
        EndTiming("spectrum", &Time);
        
        BeginTiming("spectrum", &Time);
        OperatorPsi(OppSz, psi, &Opsi);
        LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[33].Gamma);
        TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[33]);
        EndTiming("spectrum", &Time);
        
        break;
    }
    
    FilePrintSpectra(&Spectra[31], 3, "NiImpurity_Spin","w");
    EndTiming("Spin", &Time);
    PrintTiming(Time);
  }
  
  if (SpecSet.Spectra & 32)
  {
    printf("\n========================================================");
    printf("\n==========              charge               ===========");
    printf("\n========================================================\n");
    BeginTiming("Charge", &Time);
    
    Spectra[34].N=6000;
    Spectra[34].Emin=-1.0;
    Spectra[34].Emax= 11.0;
    Spectra[34].dE = (Spectra[34].Emax - Spectra[34].Emin) / Spectra[34].N;
    Spectra[34].Gamma = 0.01;
    SpectrumAllocate(Spectra[34], &Spectra[34]);
    
    OperatorCopy(&TCharge, OppNd);
    OperatorAddLader(&TCharge, &oadd, 0, -1, 0);
    
    switch (SpecSet.Dyson)
    {
      case 1:
        BeginTiming("spectrum", &Time);
        OperatorPsi(TCharge, psi, &Opsi);
        LanczosTriDiagonalizeDysonLowMem(OppH, Opsi, SpecSet.zeroLancWavefunction,
                                         SpecSet.NTriDiag, VPerturbation,
                                         gPerturbation, 20, &SpecTri);
        TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[34]);
        EndTiming("spectrum", &Time);
        
        break;
        
      default:
        BeginTiming("spectrum", &Time);
        OperatorPsi(TCharge, psi, &Opsi);
        LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[34].Gamma);
        TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[34]);
        EndTiming("spectrum", &Time);
        
        break;
    }
    
    FilePrintSpectra(&Spectra[34], 1, "NiImpurity_Charge","w");
    EndTiming("Charge", &Time);
    PrintTiming(Time);
  }
  
  EndTiming("Total", &Time);
  PrintTiming(Time);
  
  return 0;
}
