#include "../Yi/AtomicFunctions.c"

#define MAXNSHELLS 8
#define MAXCOULOMBVV 6
#define MAXCOULOMBCV 6

typedef struct ShellType
{
  int n,l; //n and l quantum numbers
} ShellType;
typedef struct ModelSettings
{
  
  uint32_t Model; /*0 for atomic, 1 for crystal field, 2 for ligand field, 3 for impurity*/
  
  int BathType;  /*0 for using NMTO for bath, 1 for circular model bath, 2 for rect*/
  int BathGeometry; /*0 for star, 1 for chain, 2 for natural*/
  char CTRLFile[100], HAMRFile[100];
  double DOSMin,DOSMax,EFermi;
  uint32_t NPoles; /* Number of poles in the Green's function representation */
  double kBT;
  int NkPoints[3];
  uint32_t NBathV,NBathC;
  double BandWidth;
  
  uint32_t AFM; /*1 for AFM doubling of NMTO unit cell and AFT +U.  */
  int SuperCell[3][3];
  double *UFactors;
  int NUFactors;
  
  int Perturbation; /* Set to 1 for using Dyson spectra calcs */
  int NCellV,NCellC; /* The number of sites to include in many body portion for Perturbative calcs */
  
  uint32_t NFermion, NFermionShift;
  
  uint32_t AtomicZ,Nd;
  uint32_t NValence,NCore; /* number of valence/core shells */
  ShellType CoreShells[MAXNSHELLS],ValenceShells[MAXNSHELLS];
  uint32_t valenceFillings[MAXNSHELLS];
  double Fvv[MAXNSHELLS][MAXCOULOMBVV],Fcv[MAXNSHELLS][MAXNSHELLS][MAXCOULOMBCV],Gcv[MAXNSHELLS][MAXNSHELLS][MAXCOULOMBCV];
  double zetac[MAXNSHELLS],zetav[MAXNSHELLS];
  double Delta, Udd, Upd,Ded,UMeanField;
  double Kdd,Kpd;
  double BMag[3], HExch[3];
  int numIEShifts,numLEShifts,numHopping;
  double *IEShifts,*LEShifts; //impurity and ligand orbital energy shifts (i.e. crystal field and Tpp)
  double *Hopping; //hopping magnitude from each impurity orbital to bath/ligand

  
  double epsilon, zero;
  uint32_t NBitsKey;
  uint32_t DoBands, DoDOS, DoGrd, DoXAS, DoFY, DoRIXS, DoXPS, DoPES, DoIPES, DoSpin, DoCharge;
  uint32_t NKrylovXAS, NKrylovResRIXS, NKrylovRIXS, NKrylovFY;
  double wXASLow, wXASStep, wXASHigh, LXAS, GXAS;
  double wResRIXSLow, wResRIXSStep, wResRIXSHigh, wRIXSLow, wRIXSStep, wRIXSHigh, LResRIXS, LRIXS, GRIXS;
  double wFYLow, wFYStep, wFYHigh, LFY, GFY;
  
} ModelSettings;



uint32_t StreamReadModelSettings(ModelSettings * O, FILE * stream)
{
  char commandStr[100];
  char valueStr[100],valueStr2[9][10];
  char tmpbuff[200];
  char *pch,*token;
  const char delim[2] = " ";
  int i,imax;
  
  while(fgets(tmpbuff,100,stream) != 0)
  {
    sscanf(tmpbuff, "%15s : %15s  #",commandStr,valueStr);
    if(strcmp(commandStr,"Model") == 0)
    {
      O->Model = atoi(valueStr);
    }
    else if(strcmp(commandStr,"BathType") == 0)
    {
      O->BathType = atoi(valueStr);
    }
    else if(strcmp(commandStr,"BathGeometry") == 0)
    {
      O->BathGeometry = atoi(valueStr);
    }
    else if(strcmp(commandStr,"AFM") == 0)
    {
      O->AFM = atoi(valueStr);
    }
    else if(strcmp(commandStr,"SuperCell") == 0)
    {
      sscanf(tmpbuff, "%15s : %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s",commandStr,valueStr,valueStr2[0],valueStr2[1],valueStr2[2],valueStr2[3],valueStr2[4],valueStr2[5],valueStr2[6],valueStr2[7],valueStr2[8]);
      for(i = 0; i < 9; i++)
      {
        O->SuperCell[i/3][i%3] = atoi(valueStr2[i]);
      }
    }
    else if(strcmp(commandStr,"CTRL") == 0)
    {
      strcpy(O->CTRLFile,valueStr);
    }
    else if(strcmp(commandStr,"HAMR") == 0)
    {
      strcpy(O->HAMRFile,valueStr);
    }
    else if(strcmp(commandStr,"DOSMin") == 0)
    {
      O->DOSMin = atof(valueStr);
    }
    else if(strcmp(commandStr,"DOSMax") == 0)
    {
      O->DOSMax = atof(valueStr);
    }
    else if(strcmp(commandStr,"BandWidth") == 0)
    {
      O->BandWidth = atof(valueStr);
    }
    else if(strcmp(commandStr,"EFermi") == 0)
    {
      O->EFermi = atof(valueStr);
    }
    else if(strcmp(commandStr,"kBT") == 0)
    {
      O->kBT = atof(valueStr);
    }
    else if(strcmp(commandStr,"Perturbation") == 0)
    {
      O->Perturbation = atoi(valueStr);
    }
    else if(strcmp(commandStr,"NCellV") == 0)
    {
      O->NCellV = atoi(valueStr);
    }
    else if(strcmp(commandStr,"NCellC") == 0)
    {
      O->NCellC = atoi(valueStr);
    }
    else if(strcmp(commandStr,"NPoles") == 0)
    {
      O->NPoles = atoi(valueStr);
    }
    else if(strcmp(commandStr,"Nk") == 0)
    {      
      sscanf(tmpbuff, "%15s : %d %d %d",commandStr,&O->NkPoints[0],&O->NkPoints[1],&O->NkPoints[2]);
    }
    else if(strcmp(commandStr, "IEShifts") == 0)
    {
      pch = strchr(tmpbuff,':');
      token = strtok(pch+1,delim); //stores the number of eshifts now
      O->numIEShifts = atoi(token);
      O->IEShifts = calloc(O->numIEShifts,sizeof(double));

      for(i=0; i < O->numIEShifts; i++)
      {
        token = strtok(NULL,delim);
        O->IEShifts[i] = atof(token);
      }
    }
    else if(strcmp(commandStr, "LEShifts") == 0)
    {
      pch = strchr(tmpbuff,':');
      token = strtok(pch+1,delim); //stores the number of eshifts now
      O->numLEShifts = atoi(token);
      O->LEShifts = calloc(O->numLEShifts,sizeof(double));

      for(i=0; i < O->numLEShifts; i++)
      {
        token = strtok(NULL,delim);
        O->LEShifts[i] = atof(token);
      }
    }
    else if(strcmp(commandStr, "Hopping") == 0)
    {
      pch = strchr(tmpbuff,':');
      token = strtok(pch+1,delim); //stores the number of eshifts now
      O->numHopping = atoi(token);
      O->Hopping = calloc(O->numHopping,sizeof(double));

      for(i=0; i < O->numHopping; i++)
      {
        token = strtok(NULL,delim);
        O->Hopping[i] = atof(token);
      }
    }
    else if(strcmp(commandStr,"NBathV") == 0)
    {
      O->NBathV = atoi(valueStr);
    }
    else if(strcmp(commandStr,"NBathC") == 0)
    {
      O->NBathC = atoi(valueStr);
    }
    else if(strcmp(commandStr,"AtomicZ") == 0)
    {
      O->AtomicZ = atoi(valueStr);
    }
    else if(strcmp(commandStr,"Nd") == 0)
    {
      O->Nd = atoi(valueStr);
    }
    else if(strcmp(commandStr,"UMeanField") == 0)
    {
      O->UMeanField = atof(valueStr);
    }
    else if(strcmp(commandStr,"Delta") == 0)
    {
      O->Delta = atof(valueStr);
    }
    else if(strcmp(commandStr,"Ded") == 0)
    {
      O->Ded = atof(valueStr);
    }
    else if(strcmp(commandStr,"Udd") == 0)
    {
      O->Udd = atof(valueStr);
    }
    else if(strcmp(commandStr,"Upd") == 0)
    {
      O->Upd = atof(valueStr);
    }
    else if(strcmp(commandStr,"Kdd") == 0)
    {
      O->Kdd = atof(valueStr);
    }
    else if(strcmp(commandStr,"Kpd") == 0)
    {
      O->Kpd = atof(valueStr);
    }
    else if(strcmp(commandStr,"BMag") == 0)
    {
      sscanf(tmpbuff, "%15s : %15s %15s %15s",commandStr,valueStr2[0],valueStr2[1],valueStr2[2]);
      for(i = 0; i < 3; i++)
      {
        O->BMag[i] = atof(valueStr2[i]);
      }
    }
    else if(strcmp(commandStr,"HExch") == 0)
    {
      sscanf(tmpbuff, "%15s : %15s %15s %15s",commandStr,valueStr2[0],valueStr2[1],valueStr2[2]);
      for(i = 0; i < 3; i++)
      {
        O->HExch[i] = atof(valueStr2[i]);
      }
    }
    else if(strcmp(commandStr,"epsilon") == 0)
    {
      O->epsilon = atof(valueStr);
    }
    else if(strcmp(commandStr,"zero") == 0)
    {
      O->zero = atof(valueStr);
    }
    else if(strcmp(commandStr,"NBitsKey") == 0)
    {
      O->NBitsKey = atoi(valueStr);
    }
    else if(strcmp(commandStr,"DoGrd") == 0)
    {
      O->DoGrd = atoi(valueStr);
    }
    else if(strcmp(commandStr,"DoBands") == 0)
    {
      O->DoBands = atoi(valueStr);
    }
    else if(strcmp(commandStr,"DoDOS") == 0)
    {
      O->DoDOS = atoi(valueStr);
    }
    else if(strcmp(commandStr,"DoXAS") == 0)
    {
      O->DoXAS = atoi(valueStr);
    }
    else if(strcmp(commandStr,"DoFY") == 0)
    {
      O->DoFY = atoi(valueStr);
    }
    else if(strcmp(commandStr,"DoRIXS") == 0)
    {
      O->DoRIXS = atoi(valueStr);
    }
    else if(strcmp(commandStr,"DoXPS") == 0)
    {
      O->DoXPS = atoi(valueStr);
    }
    else if(strcmp(commandStr,"DoPES") == 0)
    {
      O->DoPES = atoi(valueStr);
    }
    else if(strcmp(commandStr,"DoIPES") == 0)
    {
      O->DoIPES = atoi(valueStr);
    }
    else if(strcmp(commandStr,"DoSpin") == 0)
    {
      O->DoSpin = atoi(valueStr);
    }
    else if(strcmp(commandStr,"DoCharge") == 0)
    {
      O->DoCharge = atoi(valueStr);
    }
    else if(strcmp(commandStr,"NFermionShift") == 0)
    {
      O->NFermionShift = atoi(valueStr);
    }
    else if(strcmp(commandStr,"NKrylovXAS") == 0)
    {
      O->NKrylovXAS = atoi(valueStr);
    }
    else if(strcmp(commandStr,"wXASLow") == 0)
    {
      O->wXASLow = atof(valueStr);
    }
    else if(strcmp(commandStr,"wXASStep") == 0)
    {
      O->wXASStep = atof(valueStr);
    }
    else if(strcmp(commandStr,"wXASHigh") == 0)
    {
      O->wXASHigh = atof(valueStr);
    }
    else if(strcmp(commandStr,"LXAS") == 0)
    {
      O->LXAS = atof(valueStr);
    }
    else if(strcmp(commandStr,"GXAS") == 0)
    {
      O->GXAS = atof(valueStr);
    }
    else if(strcmp(commandStr,"NKrylovResRIXS") == 0)
    {
      O->NKrylovResRIXS = atoi(valueStr);
    }
    else if(strcmp(commandStr,"NKrylovRIXS") == 0)
    {
      O->NKrylovRIXS = atoi(valueStr);
    }
    else if(strcmp(commandStr,"wResRIXSLow") == 0)
    {
      O->wResRIXSLow = atof(valueStr);
    }
    else if(strcmp(commandStr,"wResRIXSStep") == 0)
    {
      O->wResRIXSStep = atof(valueStr);
    }
    else if(strcmp(commandStr,"wResRIXSHigh") == 0)
    {
      O->wResRIXSHigh = atof(valueStr);
    }
    else if(strcmp(commandStr,"wRIXSLow") == 0)
    {
      O->wRIXSLow = atof(valueStr);
    }
    else if(strcmp(commandStr,"wRIXSStep") == 0)
    {
      O->wRIXSStep = atof(valueStr);
    }
    else if(strcmp(commandStr,"wRIXSHigh") == 0)
    {
      O->wRIXSHigh = atof(valueStr);
    }
    else if(strcmp(commandStr,"LResRIXS") == 0)
    {
      O->LResRIXS = atof(valueStr);
    }
    else if(strcmp(commandStr,"LRIXS") == 0)
    {
      O->LRIXS = atof(valueStr);
    }
    else if(strcmp(commandStr,"GRIXS") == 0)
    {
      O->GRIXS = atof(valueStr);
    }
    else if(strcmp(commandStr,"NKrylovFY") == 0)
    {
      O->NKrylovFY = atoi(valueStr);
    }
    else if(strcmp(commandStr,"wFYLow") == 0)
    {
      O->wFYLow = atof(valueStr);
    }
    else if(strcmp(commandStr,"wFYStep") == 0)
    {
      O->wFYStep = atof(valueStr);
    }
    else if(strcmp(commandStr,"wFYHigh") == 0)
    {
      O->wFYHigh = atof(valueStr);
    }
    else if(strcmp(commandStr,"LFY") == 0)
    {
      O->LFY = atof(valueStr);
    }
    else if(strcmp(commandStr,"GFY") == 0)
    {
      O->GFY = atof(valueStr);
    }
  }
  
  
  return 0;
}
int32_t InitModelSettings(ModelSettings * O)
{
  int i,j;
  O->Model = 2;  //ligand field
  O->BathType = 2; //irrelevant for ligand field
  O->BathGeometry = 0;
  O->AFM = 0;
  for(i = 0; i < 3; i++)
    for(j = 0; j < 3; j++)
      if(i == j)
        O->SuperCell[i][j] = 1;
      else
        O->SuperCell[i][j] = 0;
  
    
  O->BMag[0] = 0.0; O->BMag[1] = 0.0; O->BMag[2] = 0.0;
  O->HExch[0] = 0.0; O->HExch[1] = 0.0; O->HExch[2] = 0.00001;
  
  O->numIEShifts = 0;
  O->numHopping = 0;
  O->numLEShifts = 0;
  
  O->BandWidth = 4.0;
  
  O->Perturbation = 0;
  O->NCellV = 2;
  O->NCellC = 2;
  
  O->DoBands = 0;
  O->DoCharge = 0;
  O->DoDOS = 0;
  O->DoFY = 0;
  O->DoIPES = 0;
  O->DoPES = 0;
  O->DoRIXS = 0;
  O->DoSpin = 0;
  O->DoXPS = 0;
  
  O->EFermi = 0.0;
  O->GFY = 0.3;
  O->GRIXS = 0.3;
  O->GXAS = 0.3;
  O->Kdd = 0.0;
  O->Kpd = 0.0;
  O->kBT = 0.000000001;
  O->LFY = 0.45;
  O->LRIXS = 0.05;
  O->LResRIXS = 0.45;
  O->LXAS = 0.45;
  
  O->NBathC = 1;
  O->NBathV = 1;
  
  O->NBitsKey = 16;
  O->NFermion = 26; //2p3d ligand field
  O->NFermionShift = 6;
  O->NKrylovFY = 100;
  O->NKrylovRIXS = 150;
  O->NKrylovResRIXS = 100;
  O->NKrylovXAS = 100;
  O->NPoles = 1000;
  O->Nd = 8;
  O->NkPoints[0] = 40;
  O->NkPoints[1] = 40;
  O->NkPoints[2] = 40;
  
  O->UMeanField = 0.0;
  O->Udd = 6.0;
  O->Upd = 7.5;
  O->epsilon = 0.001;
  O->wFYHigh = -50.0;
  O->wFYLow = -80.0;
  O->wFYStep = 0.05;
  O->wRIXSHigh = 15.0;
  O->wRIXSLow = -1.0;
  O->wRIXSStep = 0.02;
  O->wResRIXSHigh = -65.0;
  O->wResRIXSLow = -75.0;
  O->wResRIXSStep = 0.5;
  O->wXASHigh = -50.0;
  O->wXASLow = -80.0;
  O->wXASStep = 0.02;
  
  O->Ded = 0.0;
  O->AtomicZ = 28;
  O->Delta = 4.5;
  O->DoGrd = 1;
  O->DoXAS = 1;
  
  return 0;
  
}
uint32_t FileReadModelSettings(ModelSettings * O, const char *const flname)
{
  FILE  *stream;
  uint32_t retval;
  if ((stream = fopen(flname, "r")) == NULL)
  {
    printf("ERROR\nCan't open %s for reading.\n", flname);
    return 1;
  }
  retval = StreamReadModelSettings(O, stream);
  fclose(stream);
  return retval;
}

uint32_t StreamPrintModelSettings(ModelSettings O, FILE * stream)
{
  int i;
  fprintf(stream, "\nSettings of Impurity calculation\n");
  fprintf(stream, "==================================\n");
  fprintf(stream, "Model: %d\n",O.Model);
  fprintf(stream, "BathType: %d\n",O.BathType);
  fprintf(stream, "BathGeometry: %d\n",O.BathGeometry);
  fprintf(stream, "Perturbation: %d\n",O.Perturbation);
  fprintf(stream, "NCellV: %d\n",O.NCellV);
  fprintf(stream, "NCellC: %d\n",O.NCellC);
  fprintf(stream, "AFM: %d\n",O.AFM);
  fprintf(stream, "SuperCell: ");
  for(i = 0; i < 9; i++)
    fprintf(stream, "%d ",O.SuperCell[i/3][i%3]);
  fprintf(stream, "\n");
  fprintf(stream, "CTRL: %s\n",O.CTRLFile);
  fprintf(stream, "HAMR: %s\n",O.HAMRFile);
  fprintf(stream, "DOSMin: %f\n",O.DOSMin);
  fprintf(stream, "DOSMax: %f\n",O.DOSMax);
  fprintf(stream, "BandWidth: %f\n", O.BandWidth);
  fprintf(stream, "EFermi: %f\n",O.EFermi);
  fprintf(stream, "kBT: %f\n",O.kBT);
  fprintf(stream, "NPoles: %d\n",O.NPoles);
  fprintf(stream, "Nk: %d %d %d\n",O.NkPoints[0],O.NkPoints[1],O.NkPoints[2]);
  fprintf(stream, "NBathV: %d\n",O.NBathV);
  fprintf(stream, "NBathC: %d\n",O.NBathC);
  fprintf(stream, "NFermionShift: %d\n",O.NFermionShift);
  
  fprintf(stream, "AtomicZ: %d\n",O.AtomicZ);
  fprintf(stream, "Nd: %d\n",O.Nd);
  fprintf(stream, "UMeanField: %f\n",O.UMeanField);
  fprintf(stream, "Delta: %f\n",O.Delta);
  fprintf(stream, "Ded: %f\n",O.Ded);
  fprintf(stream, "Udd: %f\n",O.Udd);
  fprintf(stream, "Upd: %f\n",O.Upd);
  fprintf(stream, "Kdd: %f\n",O.Kdd);
  fprintf(stream, "Kpd: %f\n",O.Kpd);
  if(O.numIEShifts > 0)
  {
    fprintf(stream, "IEShifts: ");
    for(i=0; i < O.numIEShifts; i++)
    {
      fprintf(stream, "%f ",O.IEShifts[i]);
    }  
    fprintf(stream,"\n");
  }
  if(O.numLEShifts > 0)
  {
    fprintf(stream, "LEShifts: ");
    for(i=0; i < O.numLEShifts; i++)
    {
      fprintf(stream, "%f ",O.LEShifts[i]);
    }  
    fprintf(stream,"\n");
  }
  if(O.numHopping > 0)
  {
    fprintf(stream, "Hopping: ");
    for(i=0; i < O.numHopping; i++)
    {
      fprintf(stream, "%f ",O.Hopping[i]);
    }  
    fprintf(stream,"\n");
  }
  fprintf(stream, "epsilon: %g\n",O.epsilon);
  fprintf(stream, "zero: %g\n",O.zero);
  fprintf(stream, "NBitsKey: %d\n",O.NBitsKey);
  
  fprintf(stream, "DoGrd: %d\n",O.DoGrd);
  fprintf(stream, "DoBands: %d\n",O.DoBands);
  fprintf(stream, "DoDOS: %d\n",O.DoDOS);
  fprintf(stream, "DoXAS: %d\n",O.DoXAS);
  fprintf(stream, "DoFY: %d\n",O.DoFY);
  fprintf(stream, "DoRIXS: %d\n",O.DoRIXS);
  fprintf(stream, "DoXPS: %d\n",O.DoXPS);
  fprintf(stream, "DoPES: %d\n",O.DoPES);
  fprintf(stream, "DoIPES: %d\n",O.DoIPES);
  fprintf(stream, "DoSpin: %d\n",O.DoSpin);
  fprintf(stream, "DoCharge: %d\n",O.DoCharge);
  
  fprintf(stream, "NKrylovXAS: %d\n",O.NKrylovXAS);
  fprintf(stream, "wXAS: %f  %f  %f\n",O.wXASLow, O.wXASStep, O.wXASHigh);
  fprintf(stream, "GXAS: %f\n",O.GXAS);
  fprintf(stream, "LXAS: %f\n",O.LXAS);
  fprintf(stream, "NKrylovResRIXS: %d\n",O.NKrylovResRIXS);
  fprintf(stream, "wResRIXS: %f  %f  %f\n",O.wResRIXSLow, O.wResRIXSStep, O.wResRIXSHigh);
  fprintf(stream, "LResRIXS: %f\n",O.LResRIXS);
  fprintf(stream, "NKrylovRIXS: %d\n",O.NKrylovRIXS);
  fprintf(stream, "wRIXS: %f  %f  %f\n",O.wRIXSLow, O.wRIXSStep, O.wRIXSHigh);
  fprintf(stream, "LRIXS: %f\n",O.LRIXS);
  fprintf(stream, "GRIXS: %f\n",O.GRIXS);
  fprintf(stream, "NKrylovFY: %d\n",O.NKrylovFY);
  fprintf(stream, "wFY: %f  %f  %f\n",O.wFYLow, O.wFYStep, O.wFYHigh);
  fprintf(stream, "GFY: %f\n",O.GFY);
  fprintf(stream, "LFY: %f\n",O.LFY);
  fprintf(stream, "==================================\n");

  return 0;
}

uint32_t PrintModelSettings(ModelSettings O, const char *const flname)
{
  FILE  *stream;
  uint32_t retval;
  
  if(strcmp(flname,"stdout") == 0)
    return StreamPrintModelSettings(O, stdout);
  else
  {
    if ((stream = fopen(flname, "w")) == NULL)
    {
      printf("ERROR\nCan't open %s for writing.\n", flname);
      return 1;
    }
    retval = StreamPrintModelSettings(O, stream);
    fclose(stream);
    return retval;
  }
}
/*Initializes a wavefunction on a natural impurity basis. Assumes fermion ordering is:
 * 
 * | NFermionShift | NImpOrb | NImpOrb | NBathV | NBathC |
 * 
 * which correspond to core orbitals, impurity valence orbitals, ligand site orbitals,
 * valence chain orbitals, and conduction chain orbitals.
 * 
 * The wavefunction is initialized with all possible determinants having minNIm to maxNIm
 * electrons on the impurity site. For a given number N on impurity, NImpOrb - N are put
 * on ligand site. The core orbitals are always full, the valence chain is always full,
 * and the conduction chain always empty */
uint32_t InitNaturalImpurityWavefunctionRandom(WaveFunctionType * const psi, uint32_t NBitsKey, uint32_t const NFermions, uint32_t NFermionShift, uint32_t NImpOrb, uint32_t NBathV,uint32_t minNIm,uint32_t maxNIm)
{
  int i, j, xImp, yImp, zImp, xLig, yLig, zLig, pImp[NImpOrb+2], pLig[NImpOrb+2], bImp[NImpOrb],bLig[NImpOrb], count=0;
  uint8_t *Determinant;
  
  psi->NBitsKey = NBitsKey;
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
  
  for(i = minNIm; i <= maxNIm; i++)
  {
    for(j = 0; j < NImpOrb; j++)
    {
      if(j < i)
      {
        bImp[j] = 1;
        bLig[j] = 0;
      }
      else
      {
        bImp[j] = 0;
        bLig[j] = 1;
      }
    }
    Determinant = calloc(psi->NWord, sizeof(uint8_t));
    for (j = 0; j < NFermionShift; j++)
      OperatorCr(Determinant, j, psi->NFermion);
  
    for (j = NFermionShift; j < (NFermionShift+NImpOrb); j++)
    {
      if (bImp[j-NFermionShift])
        OperatorCr(Determinant, j, psi->NFermion);
      if(bLig[j-NFermionShift])
        OperatorCr(Determinant, j+NImpOrb, psi->NFermion);
    }
    for (j = (NFermionShift+2*NImpOrb); j < (NFermionShift+(2+NBathV)*NImpOrb); j++)
    {
      OperatorCr(Determinant, j, psi->NFermion);
    }
    RealWaveFunctionAddElement(psi, Determinant, 1.); free(Determinant);
    count++;
  
    inittwiddle(i, NImpOrb, pImp);
    while(!twiddle(&xImp, &yImp, &zImp, pImp))
    {
      bImp[xImp] = 1;
      bImp[yImp] = 0; 
      inittwiddle(NImpOrb-i, NImpOrb, pLig);
      while(!twiddle(&xLig, &yLig, &zLig, pLig))
      {
        Determinant = calloc(psi->NWord, sizeof(uint8_t));
        bLig[xLig] = 1;
        bLig[yLig] = 0;
        for (j = 0; j < NFermionShift; j++)
          OperatorCr(Determinant, j, psi->NFermion);
        for (j = NFermionShift; j < (NFermionShift+NImpOrb); j++)
        {
          if (bImp[j-NFermionShift])
            OperatorCr(Determinant, j, psi->NFermion);
          if(bLig[j-NFermionShift])
            OperatorCr(Determinant, j+NImpOrb, psi->NFermion);
        }
        for (j = (NFermionShift+2*NImpOrb); j < (NFermionShift+(2+NBathV)*NImpOrb); j++)
          OperatorCr(Determinant, j, psi->NFermion);
        RealWaveFunctionAddElement(psi, Determinant, 1.); free(Determinant);
        count++;
      }
    }
  }
  printf("The wave function has %i determinants.\n", count);
  WaveFunctionRandomize(psi);
  return 0;
}
uint32_t ImpurityHFromModelBath(OperatorType *const Himp, int geometry, int shape, double width, int NPoles,int NBathSitesAnd, int NBathSitesV, int NBathSitesC, double mu, double kBT, int blocksize, int spin, double *eshiftsImpurity, double *eshifts, double *hopping, int NFermion, int NFermionShift)
{
  /* TO DO: For rectangular bath, I should make it in Anderson geometry, and then transform as need (i.e. to tridiagonal form if desired)
            For circular bath, I should make it in Tridiagonal form, and then transform it to anderson if desired.
   */
  
  /* Shape can be rectangular or circular (DOS in Anderson hybridization) (0 or 1)
   * 
   * 
   * blocksize gives the number of orbitals (i.e. 5 for d, 7 for f). 
   * hopping and eshifts are of size blocksize and contain overall hoppings and band energy shifts for the different orbitals
   * If spin==1, then there are 5/7 hoppings/eshifts given, else 10/14 (for d/f, eg)
   *  
   * mu is chemical potential, width is DOS width
   * 
   * Can be set up as star/Anderson, chain, or natural bath geometry (geometry=0,1,2 respectively
   * 
   * NBathSitesAnd is the number of discretization points when generating the blockandersonmatrix
   * NBathSites is the final number of bath sites (in case one wants a fine mesh when discretizing, and then
   * transform to tridiagonal or natural impurity, and then truncate the sites to a smaller value
   * 
   * !!!!!!!!!!Should return a wavefunction string that gives the filling of the bath orbitals according to mu!!!!!!!!!!
   */
  BlockAndersonMatrixType AndBlockMat;
  BlockTriDiagonalMatrixType TriBlockMat;
  BlockNaturalImpurityOrbitalMatrixType NatBlockMat;
  BlockListOfPoles BlockPoles;
  int i,j,jmax,filling;
  double scaling,dos[NBathSitesAnd];


  
  /* Determine filling of bath according to chemical potential mu */
  /* -1 is empty, 0 is partially full, 1 is full */
  for(i = 0; i < blocksize; i++)
  {
    if(-width/2 + eshifts[i] < mu)
      filling += 1;
  }
  if(shape == 1)
  {
    /* Circular DOS */
    if(geometry == 0)
    {
        /*Make a general Anderson geometry matrix according to desired shape*/
        AndBlockMat.N = NBathSitesAnd+1;
        AndBlockMat.QComplex = 0;
        if(spin)
        {
          AndBlockMat.BlockSize = blocksize;
          jmax = blocksize;
        }
        else
        {
          AndBlockMat.BlockSize = 2*blocksize;
          jmax = 2*blocksize;
        }
        InitBlockAndersonMatrix(&AndBlockMat);
        AndBlockMat.mu = mu;
      /* Find normalization factor so total hopping is correct */
      /* Sum of squares for the circle should be 1, then each  */
      /* can be multiplied by the given hopping*/
    
      /* Since the endpoints of the DOS for a circle would have zero weight (and therefore
       not couple to anything), I instead define the mesh that my endpoints are in from
       the edges of the circle a bit. I divide the width into NBathSites+1 regions, which
       then has NBathSites points within the edges of the circle */
      scaling = 0;
      for(i = 1; i < NBathSitesAnd+1; i++)
      {
        dos[i-1] = sqrt((width/2)*(width/2) - (-width/2 + i*width/(NBathSitesAnd+1))*(-width/2 + i*width/(NBathSitesAnd+1)));
        scaling += dos[i-1]*dos[i-1];
      }
      for(i = 1; i < NBathSitesAnd+1; i++)
      {
        dos[i-1] /= sqrt(scaling);
      }
      /* Impurity energy shifts (crystal field) */
      for(j = 0; j < jmax; j++)
      {
        AndBlockMat.aRe[blocksize*blocksize + j*blocksize + j] = eshiftsImpurity[j];
      }
      for(i = 1; i < NBathSitesAnd+1; i++)
      {
        for(j = 0; j < jmax; j++)
        {
          /* Set diagonal values of off diagonal (hopping) block */
          AndBlockMat.bRe[(i)*blocksize*blocksize + j + j*blocksize] = dos[i-1]*hopping[j];
          /* Set diagonal values of diagonal (onsite energy) block */
          AndBlockMat.aRe[(i+1)*blocksize*blocksize + j + j*blocksize] = -width/2+i*width/(NBathSitesAnd+1)+eshifts[j];
        }
      }
    }
    else
    {
      /* Tridiagonal or natural */
      TriBlockMat.N = NBathSitesV+1;
      TriBlockMat.QComplex = 0;
      if(spin)
      {
        TriBlockMat.BlockSize = blocksize;
        jmax = blocksize;
      }
      else
      {
        TriBlockMat.BlockSize = 2*blocksize;
        jmax = 2*blocksize;
      }
      InitBlockTridiagonalMatrix(&TriBlockMat);
      TriBlockMat.mu = mu;
      
      /* Impurity energy shifts (crystal field) */
      for(j = 0; j < jmax; j++)
      {
        TriBlockMat.aRe[blocksize*blocksize + j*blocksize + j] = eshiftsImpurity[j];
        TriBlockMat.bRe[blocksize*blocksize + j*blocksize + j] = hopping[j];
      }
      for(i = 1; i < NBathSitesV; i++)
      {
        for(j = 0; j < jmax; j++)
        {
          /* Set diagonal values of off diagonal (hopping) block */
          TriBlockMat.bRe[(i+1)*blocksize*blocksize + j + j*blocksize] = width/4;
        }
      }
      for(i = 1; i < NBathSitesV+1; i++)
      {
        for(j = 0; j < jmax; j++)
        {
          /* Set diagonal values of diagonal (onsite energy) block */
          TriBlockMat.aRe[(i+1)*blocksize*blocksize + j + j*blocksize] = eshifts[j];
        }
      }
    }
  }
  else
  {
    /*only other option so far is rectangular*/
    
    /* For rectangular, the DOS endpoints are -width/2 and width/2 */
 
    /* Impurity energy shifts (crystal field) */
    for(j = 0; j < jmax; j++)
    {
      AndBlockMat.aRe[blocksize*blocksize + j*blocksize + j] = eshiftsImpurity[j];
    }
    for(i = 1; i < NBathSitesAnd+1; i++)
    {
      for(j = 0; j < jmax; j++)
      {
        /* Set diagonal values of off diagonal (hopping) block */
        AndBlockMat.bRe[(i)*blocksize*blocksize + j + j*blocksize] = hopping[j]/sqrt(NBathSitesV);
        /* Set diagonal values of diagonal (onsite energy) block */
        if(NBathSitesV > 1)
          AndBlockMat.aRe[(i+1)*blocksize*blocksize + j + j*blocksize] = -width/2+(i-1)*width/(NBathSitesAnd-1)+eshifts[j];
        else
          AndBlockMat.aRe[(i+1)*blocksize*blocksize + j + j*blocksize] = eshifts[j];
      }
    }    
  }
  
  if(geometry == 0)
  {
    BlockAndersonMatrixToImpurityHamiltonian(Himp, &AndBlockMat, NBathSitesAnd+1, spin, NFermionShift,  NFermion);
    OperatorChop(Himp, 1000*DBL_EPSILON);
  }
  else if(geometry == 1)
  {
    BlockTridiagonalMatrixToImpurityHamiltonian(Himp, &TriBlockMat, NBathSitesV+1, spin, NFermionShift, NFermion);
    OperatorChop(Himp, 1000*DBL_EPSILON);
    //BlockAndersonMatrixToDoubleBlockTridiagonalMatrix
    //DoubleBlockTridiagonalMatrixToImpurityHamiltonian
  }
  else if(geometry == 2)
  {
    BlockTriDiagonalMatrixToBlockListOfPoles(&TriBlockMat, &BlockPoles, -width*5/8, width*5/8, NPoles);
    BlockListOfPolesToBlockNaturalImpurityOrbitalMatrix(&BlockPoles, &NatBlockMat, kBT, NBathSitesV, NBathSitesC);
    BlockNaturalImpurityOrbitalMatrixChop(&NatBlockMat, 1000*DBL_EPSILON);
    BlockNaturalImpurityOrbitalMatrixToImpurityHamiltonian(Himp, &NatBlockMat, NBathSitesV, NBathSitesC, spin, NFermionShift, NFermion);  
    OperatorChop(Himp, 1000*DBL_EPSILON);
    
  }
  
  return 0;
}


uint32_t ImpurityHFromModelBathPerturbation(OperatorType *const Himp, OperatorType * const VPert, TriDiagonalMatrixType * const gPert, int NCellV, int geometry, int shape, double width, int NPoles,int NBathSitesAnd, int NBathSitesV, int NBathSitesC, double mu, double kBT, int blocksize, int spin, double *eshiftsImpurity, double *eshifts, double *hopping, int NFermion, int NFermionShift)
{
  
  /* Shape can be rectangular or circular (DOS in Anderson hybridization) (0 or 1)
   * 
   * 
   * blocksize gives the number of orbitals (i.e. 5 for d, 7 for f). 
   * hopping and eshifts are of size blocksize and contain overall hoppings and band energy shifts for the different orbitals
   * If spin==1, then there are 5/7 hoppings/eshifts given, else 10/14 (for d/f, eg)
   *  
   * mu is chemical potential, width is DOS width
   * 
   * Can be set up as star/Anderson, chain, or natural bath geometry (geometry=0,1,2 respectively
   * 
   * NBathSitesAnd is the number of discretization points when generating the blockandersonmatrix
   * NBathSites is the final number of bath sites (in case one wants a fine mesh when discretizing, and then
   * transform to tridiagonal or natural impurity, and then truncate the sites to a smaller value
   * 
   * !!!!!!!!!!Should return a wavefunction string that gives the filling of the bath orbitals according to mu!!!!!!!!!!
   */
  BlockAndersonMatrixType AndBlockMat;
  BlockTriDiagonalMatrixType TriBlockMat;
  BlockNaturalImpurityOrbitalMatrixType NatBlockMat;
  BlockListOfPoles BlockPoles;
  int i,j,jmax,filling;
  double scaling,dos[NBathSitesAnd];


  
  /* Determine filling of bath according to chemical potential mu */
  /* -1 is empty, 0 is partially full, 1 is full */
  for(i = 0; i < blocksize; i++)
  {
    if(-width/2 + eshifts[i] < mu)
      filling += 1;
  }
  if(shape == 1)
  {
    /* Circular DOS */
    if(geometry == 0)
    {
        /*Make a general Anderson geometry matrix according to desired shape*/
        AndBlockMat.N = NBathSitesAnd+1;
        AndBlockMat.QComplex = 0;
        if(spin)
        {
          AndBlockMat.BlockSize = blocksize;
          jmax = blocksize;
        }
        else
        {
          AndBlockMat.BlockSize = 2*blocksize;
          jmax = 2*blocksize;
        }
        InitBlockAndersonMatrix(&AndBlockMat);
        AndBlockMat.mu = mu;
      /* Find normalization factor so total hopping is correct */
      /* Sum of squares for the circle should be 1, then each  */
      /* can be multiplied by the given hopping*/
    
      /* Since the endpoints of the DOS for a circle would have zero weight (and therefore
       not couple to anything), I instead define the mesh that my endpoints are in from
       the edges of the circle a bit. I divide the width into NBathSites+1 regions, which
       then has NBathSites points within the edges of the circle */
      scaling = 0;
      for(i = 1; i < NBathSitesAnd+1; i++)
      {
        dos[i-1] = sqrt((width/2)*(width/2) - (-width/2 + i*width/(NBathSitesAnd+1))*(-width/2 + i*width/(NBathSitesAnd+1)));
        scaling += dos[i-1]*dos[i-1];
      }
      for(i = 1; i < NBathSitesAnd+1; i++)
      {
        dos[i-1] /= sqrt(scaling);
      }
      /* Impurity energy shifts (crystal field) */
      for(j = 0; j < jmax; j++)
      {
        AndBlockMat.aRe[blocksize*blocksize + j*blocksize + j] = eshiftsImpurity[j];
      }
      for(i = 1; i < NBathSitesAnd+1; i++)
      {
        for(j = 0; j < jmax; j++)
        {
          /* Set diagonal values of off diagonal (hopping) block */
          AndBlockMat.bRe[(i)*blocksize*blocksize + j + j*blocksize] = dos[i-1]*hopping[j];
          /* Set diagonal values of diagonal (onsite energy) block */
          AndBlockMat.aRe[(i+1)*blocksize*blocksize + j + j*blocksize] = -width/2+i*width/(NBathSitesAnd+1)+eshifts[j];
        }
      }
    }
    else
    {
      /* Tridiagonal or natural */
      TriBlockMat.N = NBathSitesV+1;
      TriBlockMat.QComplex = 0;
      if(spin)
      {
        TriBlockMat.BlockSize = blocksize;
        jmax = blocksize;
      }
      else
      {
        TriBlockMat.BlockSize = 2*blocksize;
        jmax = 2*blocksize;
      }
      InitBlockTridiagonalMatrix(&TriBlockMat);
      TriBlockMat.mu = mu;
      
      /* Impurity energy shifts (crystal field) */
      for(j = 0; j < jmax; j++)
      {
        TriBlockMat.aRe[blocksize*blocksize + j*blocksize + j] = eshiftsImpurity[j];
        TriBlockMat.bRe[blocksize*blocksize + j*blocksize + j] = hopping[j];
      }
      for(i = 1; i < NBathSitesV; i++)
      {
        for(j = 0; j < jmax; j++)
        {
          /* Set diagonal values of off diagonal (hopping) block */
          TriBlockMat.bRe[(i+1)*blocksize*blocksize + j + j*blocksize] = width/4;
        }
      }
      for(i = 1; i < NBathSitesV+1; i++)
      {
        for(j = 0; j < jmax; j++)
        {
          /* Set diagonal values of diagonal (onsite energy) block */
          TriBlockMat.aRe[(i+1)*blocksize*blocksize + j + j*blocksize] = eshifts[j];
        }
      }
    }
  }
  else
  {
    /*only other option so far is rectangular*/
    
    /* For rectangular, the DOS endpoints are -width/2 and width/2 */
 
    /* Impurity energy shifts (crystal field) */
    for(j = 0; j < jmax; j++)
    {
      AndBlockMat.aRe[blocksize*blocksize + j*blocksize + j] = eshiftsImpurity[j];
    }
    for(i = 1; i < NBathSitesAnd+1; i++)
    {
      for(j = 0; j < jmax; j++)
      {
        /* Set diagonal values of off diagonal (hopping) block */
        AndBlockMat.bRe[(i)*blocksize*blocksize + j + j*blocksize] = hopping[j]/sqrt(NBathSitesV);
        /* Set diagonal values of diagonal (onsite energy) block */
        if(NBathSitesV > 1)
          AndBlockMat.aRe[(i+1)*blocksize*blocksize + j + j*blocksize] = -width/2+(i-1)*width/(NBathSitesAnd-1)+eshifts[j];
        else
          AndBlockMat.aRe[(i+1)*blocksize*blocksize + j + j*blocksize] = eshifts[j];
      }
    }    
  }
  
  if(geometry == 0)
  {
    BlockAndersonMatrixToImpurityHamiltonian(Himp, &AndBlockMat, NBathSitesAnd+1, spin, NFermionShift,  NFermion);
    OperatorChop(Himp, 1000*DBL_EPSILON);
  }
  else if(geometry == 1)
  {
    BlockTridiagonalMatrixToImpurityHamiltonianPerturbation(Himp,VPert,gPert,&TriBlockMat,NCellV,NBathSitesV+1,spin,NFermionShift,NFermion);
    
    
    

    //BlockTridiagonalMatrixToImpurityHamiltonian(Himp, &TriBlockMat, NBathSitesV+1, spin, NFermionShift, NFermion);
    OperatorChop(Himp, 1000*DBL_EPSILON);
    //BlockAndersonMatrixToDoubleBlockTridiagonalMatrix
    //DoubleBlockTridiagonalMatrixToImpurityHamiltonian
  }
  else if(geometry == 2)
  {
    BlockTriDiagonalMatrixToBlockListOfPoles(&TriBlockMat, &BlockPoles, -width*5/8, width*5/8, NPoles);
    BlockListOfPolesToBlockNaturalImpurityOrbitalMatrix(&BlockPoles, &NatBlockMat, kBT, NBathSitesV, NBathSitesC);
    BlockNaturalImpurityOrbitalMatrixChop(&NatBlockMat, 1000*DBL_EPSILON);
    BlockNaturalImpurityOrbitalMatrixToImpurityHamiltonian(Himp, &NatBlockMat, NBathSitesV, NBathSitesC, spin, NFermionShift, NFermion);  
    OperatorChop(Himp, 1000*DBL_EPSILON);
    
  }
  
  return 0;
}

uint32_t ImpurityHFromTB(OperatorType *const Himp, ModelSettings O)
{
  StuttgartCTRL CTRLin;
  StuttgartHAMR HAMRin;
  TightBindingDef TB,TBAF;

  CompactMatrixType AT;
  TimeType Time;

  SpectraType DOS,*GSpec, *GTriSpec, *GTriSpecReduced,*SpecNat;
  BlockListOfPoles G, Sigma;
  BlockTriDiagonalMatrixType GTri, GTriReduced;
  BlockAndersonMatrixType GAndReduced;
  BlockNaturalImpurityOrbitalMatrixType GNat;
  CrystalStructureSuperCellType AFSuperCell;

  int i,j,rawBlockSize,impBlockSize;
  
  InitializeTiming(&Time);
  BeginTiming("Total", &Time);

  /* Read StuttgartNMTO files*/
  if(FileReadStuttgartCTRL(&CTRLin, O.CTRLFile))//"CTRL"))
  {
    printf("Error reading StuttgartNMTO CTRL file.\n");
    return 1;
  }
  if(FileReadStuttgartHAMR(&HAMRin, O.HAMRFile))//"HAMR"))
  {
    printf("Error reading StuttgartNMTO HAMR file.\n");
    return 1;
  }
  if(TightBindingDefFromStuttgartCTRLHAMR(CTRLin, HAMRin, &TB))
  {
    printf("Error creating TB definition from StuttgartNMTO CTRL and HAMR.");
    return 1;
  }
  
  /*Determine how many orbitals in the initial block and allocate spectra arrays*/
  rawBlockSize = 0;
  for(i = 0; i < TB.Cell.NAtoms; i++)
  {
    printf("Atom %d has atomic number %d\n", i,TB.Cell.Atom[i].AtomNumber);
    for(j = 0; j < TB.Cell.Atom[i].Nshell; j++)
    {
      printf("Shell %d of atom %d has degeneracy %d\n",j,i,TB.Cell.Atom[i].shelldegeneracy[j]);
      rawBlockSize += TB.Cell.Atom[i].shelldegeneracy[j];
      if(i == 0)
        impBlockSize = TB.Cell.Atom[i].shelldegeneracy[j]; //assumes impurity atom is first, as does the rest of this function
    }
  }
  
  
  /* Create G0 from tight binding definition */
  DOS.N=O.NPoles; DOS.Emin=O.DOSMin; DOS.Emax=O.DOSMax;
  DOS.dE=(DOS.Emax-DOS.Emin)/DOS.N;
  SpectrumAllocate(DOS, &DOS);
  BeginTiming("TightBindingDOS", &Time);
  TightBindingDOS(TB, O.NkPoints, &DOS);
  EndTiming("TightBindingDOS", &Time);
  FilePrintSpectra(&DOS, 1, "TBDOS", "w");
  
  
  if(O.AFM)
  {
    for(i = 0; i < 9; i++)
      AFSuperCell.rSuperCell[i/3][i%3] = O.SuperCell[i/3][i%3];
    AddAditionalInformationToSuperCell(&AFSuperCell);
    PrintCrystalStructureSuperCell(AFSuperCell);
    CreateTightBindingSuperCell(TB, &TBAF, AFSuperCell);
    
    Sigma.BlockSize=2*rawBlockSize;
    Sigma.QComplex=0;
    Sigma.N=0;
    InitBlockListOfPoles(&Sigma);
    for(i = 0; i < 2*rawBlockSize; i++)
    {
      if(fabs(O.UFactors[i]) < 100*DBL_EPSILON)  
        Sigma.a0Re[i + 2*rawBlockSize*i] = 0.0;
      else
        Sigma.a0Re[i + 2*rawBlockSize*i] = O.UFactors[i] * O.UMeanField + O.Ded;
    }
    
  }
  else
  {
    /* Create the LDA+U Self energy */
    Sigma.BlockSize=rawBlockSize;
    Sigma.QComplex=0;
    Sigma.N=0;
    InitBlockListOfPoles(&Sigma);
    for(i = 0; i < rawBlockSize; i++)
    {
      if(i < impBlockSize)
      {
        Sigma.a0Re[i+rawBlockSize*i] = - O.UMeanField + O.Ded;
      }
      else
      {
        Sigma.a0Re[i+rawBlockSize*i] = 0.0;
      }
    }  
    
    BeginTiming("TightBindingHamiltonianToBlockListOfPoles", &Time);
    TightBindingHamiltonianPlusSelfEnergyToBlockListOfPoles(TB, O.NkPoints, Sigma, &G, O.DOSMin, O.DOSMax, O.NPoles);
    EndTiming("TightBindingHamiltonianToBlockListOfPoles", &Time);   
    
    GSpec = malloc(rawBlockSize*rawBlockSize*sizeof(SpectraType));
    BeginTiming("BlockListOfPolesToSpectra", &Time);
    BlockListOfPolesToSpectra(G, GSpec, O.DOSMin, O.DOSMax, O.NPoles);
    EndTiming("BlockListOfPolesToSpectra", &Time);
    FilePrintSpectra(GSpec, TB.NFermion*TB.NFermion, "TBDOS_BlockListOfPoles", "w");
    free(GSpec);
    
    AT.Ni=impBlockSize; AT.Nj=G.BlockSize; AT.QComplex=0;
    InitCompactIdentityMatrix(&AT);
    
    
    if(O.BathGeometry == 0) /* Star (Anderson) */
    {
      BeginTiming("BlockListOfPolesToBlockAndersonMatrixReduced", &Time);
      BlockListOfPolesToBlockAndersonMatrixReduceDimension(G, &GAndReduced, &AT, O.NBathV+1);
      EndTiming("BlockListOfPolesToBlockAndersonMatrixReduced", &Time);
      
      PrintBlockAndersonMatrix(GAndReduced);
      
      BlockAndersonMatrixToImpurityHamiltonian(Himp, &GAndReduced, O.NBathV+1, 1, O.NFermionShift,  O.NFermion);
      
      OperatorChop(Himp, 1000*DBL_EPSILON);
      PrintOperatorReadable(*Himp);
    }
    else if(O.BathGeometry == 1) /* Chain */
    {
      /*
      BeginTiming("BlockListOfPolesToBlockTridiagonalMatrix", &Time);
      BlockListOfPolesToBlockTridiagonalMatrix(G, &GTri, O.NPoles, DBL_EPSILON);
      EndTiming("BlockListOfPolesToBlockTridiagonalMatrix", &Time);

      GTriSpec = malloc(rawBlockSize*rawBlockSize*sizeof(SpectraType));
      BeginTiming("BlockTriDiagonalMatrixToSpectra", &Time);
      BlockTriDiagonalMatrixToSpectra(&GTri, GTriSpec, O.DOSMin, O.DOSMax, O.NPoles, 0.05);
      EndTiming("BlockTriDiagonalMatrixToSpectra", &Time);
      FilePrintSpectra(GTriSpec, rawBlockSize*rawBlockSize, "TBDOS_TridiagonalMatrix", "w");
      free(GTriSpec);
      */

      BeginTiming("BlockListOfPolesToBlockTridiagonalMatrixReduced", &Time);
      BlockListOfPolesToBlockTridiagonalMatrixReduceDimension(G, &GTriReduced, &AT, O.NBathV+1,DBL_EPSILON);//O.NPoles, DBL_EPSILON);
      EndTiming("BlockListOfPolesToBlockTridiagonalMatrixReduced", &Time);
  
      PrintBlockTriDiagonalMatrix(GTriReduced);
      
      GTriSpecReduced = malloc(impBlockSize*impBlockSize*sizeof(SpectraType));
      BeginTiming("BlockTriDiagonalMatrixToSpectra", &Time);
      BlockTriDiagonalMatrixToSpectra(&GTriReduced, GTriSpecReduced, O.DOSMin, O.DOSMax, O.NPoles, 0.05);
      EndTiming("BlockTriDiagonalMatrixToSpectra", &Time);
      FilePrintSpectra(GTriSpecReduced, impBlockSize*impBlockSize, "TBDOS3reduced", "w");
      free(GTriSpecReduced);
      
      BlockTridiagonalMatrixToImpurityHamiltonian(Himp, &GTriReduced, O.NBathV+1, 1, O.NFermionShift, O.NFermion);
      OperatorChop(Himp, 1000*DBL_EPSILON);
    }
    else if(O.BathGeometry == 2) /* Natural Impurity */
    {
      G.mu = O.EFermi;
     //TightBindingHamiltonianPlusSelfEnergyToBlockListOfPoles
      BlockListOfPolesToBlockNaturalImpurityOrbitalMatrixReduceDimension(&G, &GNat, &AT, O.kBT, O.NBathV, O.NBathC);
      BlockNaturalImpurityOrbitalMatrixChop(&GNat, 1000*DBL_EPSILON);
      BlockNaturalImpurityOrbitalMatrixToImpurityHamiltonian(Himp, &GNat, O.NBathV, O.NBathC, 1, O.NFermionShift, O.NFermion);  
      OperatorChop(Himp, 1000*DBL_EPSILON);
      PrintOperatorReadable(*Himp);
      SpecNat = malloc(rawBlockSize*rawBlockSize*sizeof(SpectraType));
      BlockNaturalImpurityOrbitalToSpectra(&GNat, SpecNat, O.DOSMin, O.DOSMax, O.NPoles, 0.1);
      FilePrintSpectra(SpecNat, impBlockSize*impBlockSize, "SpecNat", "w");
      free(SpecNat);
    
    }
    else
    {
      printf("BathGeometry of type %d not known.\n",O.BathGeometry);
      printf("Options are: \n");
      printf(" 0. Star (Anderson) \n");
      printf(" 1. Chain \n");
      printf(" 2. Natural \n");
      return 1;
    }
    
  }
  EndTiming("Total", &Time);
  printf("\n\n***********************************************************************\n\n\n");
  PrintTiming(Time); 
  return 0;
}

/* 
 * Gives the valence-valence Coulomb interactions.
 * Z: atomic number
 * n: principal quantum number of shell
 * l: azimuthal quantum number of shell
 * vCount: number of electrons in shell
 * Fvv: array storing the slater integrals to be returned
 * 
 * Assumes Fvv is allocated correctly
 */
uint32_t GetIntraShellAtomicInteractions(uint32_t Z, uint32_t n, uint32_t ell, uint32_t vCount, double Fvv[MAXNSHELLS][MAXCOULOMBVV],uint32_t iShell)
{
  /*                Z n l vCount [list of Fvv] */
  char data[] = "Z 20 3 1 0 0.0 0.0 0.0 "
                "Z 25 3 1 5 0.0 10.316 6.414 "
                "Z 27 3 1 7 0.0 11.605 7.209 "
                "Z 28 3 1 7 0.0 13.277 8.295 "
                "Z 28 3 1 8 0.0 12.234 7.598 "
                "Z 29 3 1 9 0.0 13.866 8.670 "
                "X";
  char buff[10];
  int i,j,k;
  
  i = 0; j = 0;
  while(data[i] != 'X')
  {
    if(data[i] == 'Z')
    {
      i = i + 2; //char[i] is now first digit of atomic number]
      j = 1;
      while(data[i+j] != ' ')
        j++;
      memcpy(buff,&data[i],j);
      if(atoi(buff) == Z)
      {
        i = i + j + 1; //char[i] is now the first digit of n
        j = 1;
        while(data[i+j] != ' ')
          j++;
        memcpy(buff,&data[i],j);
        if(atoi(buff) == n)
        {
          i = i + j + 1; //char[i] is now the first digit of ell
          j = 1;
          while(data[i+j] != ' ')
            j++;
          memcpy(buff,&data[i],j);
          if(atoi(buff) == ell)
          {
            i = i + j + 1; //char[i] is now the first digit of vCount
            j = 1;
            while(data[i+j] != ' ')
              j++;
            memcpy(buff,&data[i],j);
            if(atoi(buff) == vCount)
            {
              i = i + j; //char[i] is now the space before first Slater integral
              k = 0;
              while(data[i+1] != 'Z' && data[i+1] != 'X')
              {
                j = 1;
                while(data[i+j] != ' ')
                  j++;
                memcpy(buff,&data[i],j);
                Fvv[iShell][k] = atof(buff);
                k++;
                i = i+j;
              }
              return 0;
            }
          }
        }
      }
    }
    i++;
  }
  printf("Configuration with Z = %d, electron count of %d not found.",Z,vCount);
  return 1;
}
uint32_t GetSpinOrbitEnergy(uint32_t Z, uint32_t n, uint32_t ell, uint32_t vCount, double *zeta)
{
  /* Gives the atomic spin orbit energy (in eV) for a given shell (sets zeta, returns 0 for success) */
  /* vCount is the number of valence electrons (d electrons for 3d TM, f electrons for rare earths 
     One can give a vCount of -1 to get standard value (charge neutral case) */
  /*                Z n l vCount zeta */
  char data[] = "Z 22 2 1 -1 3.775 "
                "Z 22 3 2 -1 0.015 "
                "Z 23 2 1 -1 4.649 "
                "Z 23 3 2 -1 0.021 "
                "Z 24 2 1 -1 5.667 "
                "Z 24 3 2 -1 0.029 "
                "Z 25 2 1 -1 6.845 "
                "Z 25 3 2 -1 0.039 "
                "Z 26 2 1 -1 8.199 "
                "Z 26 3 2 -1 0.051 "
                "Z 27 2 1 -1 9.746 "
                "Z 27 3 2 -1 0.065 "
                "Z 28 2 1 -1 11.506 "
                "Z 28 3 2 -1 0.081 "
                "Z 29 2 1 -1 13.496 "
                "Z 29 3 2 -1 0.100 "
                "Z 57 3 2 -1 6.802 "
                "Z 57 4 2 -1 1.100 "
                "Z 57 4 3 -1 0.059 "
                "Z 58 3 2 -1 7.446 "
                "Z 58 4 2 -1 1.242 "
                "Z 58 4 3 -1 0.087 "
                "Z 59 3 2 -1 8.142 "
                "Z 59 4 2 -1 1.358 "
                "Z 59 4 3 -1 0.089 "
                "Z 60 3 2 -1 8.882 "
                "Z 60 4 2 -1 1.499 "
                "Z 60 4 3 -1 0.105 "
                "Z 61 3 2 -1 9.672 "
                "Z 61 4 2 -1 1.650 "
                "Z 61 4 3 -1 0.122 "
                "Z 62 3 2 -1 10.514 "
                "Z 62 4 2 -1 1.812 "
                "Z 62 4 3 -1 0.140 "
                "Z 63 3 2 -1 11.409 "
                "Z 63 4 2 -1 1.984 "
                "Z 63 4 3 -1 0.160 "
                "Z 64 3 2 -1 12.361 "
                "Z 64 4 2 -1 2.168 "
                "Z 64 4 3 -1 0.181 "
                "Z 65 3 2 -1 13.372 "
                "Z 65 4 2 -1 2.363 "
                "Z 65 4 3 -1 0.204 "
                "Z 66 3 2 -1 14.444 "
                "Z 66 4 2 -1 2.570 "
                "Z 66 4 3 -1 0.229 "
                "Z 67 3 2 -1 15.579 "
                "Z 67 4 2 -1 2.790 "
                "Z 67 4 3 -1 0.255 "
                "Z 68 3 2 -1 16.782 "
                "Z 68 4 3 -1 0.283 "
                "Z 69 3 2 -1 18.053 "
                "Z 69 4 3 -1 0.311 "
                "Z 70 3 2 -1 19.391 "
                "Z 70 4 2 -1 3.557 "
                "Z 70 4 3 -1 0.365 "
                "X";
  char buff[10];
  int i,j,k;
  
  i = 0; j = 0;
  while(data[i] != 'X')
  {
    if(data[i] == 'Z')
    {
      i = i + 2; //char[i] is now first digit of atomic number]
      j = 1;
      while(data[i+j] != ' ')
        j++;
      memcpy(buff,&data[i],j);
      if(atoi(buff) == Z)
      {
        i = i + j + 1; //char[i] is now the first digit of n
        j = 1;
        while(data[i+j] != ' ')
          j++;
        memcpy(buff,&data[i],j);
        if(atoi(buff) == n)
        {
          i = i + j + 1; //char[i] is now the first digit of ell
          j = 1;
          while(data[i+j] != ' ')
            j++;
          memcpy(buff,&data[i],j);
          if(atoi(buff) == ell)
          {
            i = i + j + 1; //char[i] is now the first digit of vCount
            j = 1;
            while(data[i+j] != ' ')
              j++;
            memcpy(buff,&data[i],j);
            if(atoi(buff) == vCount)
            {
              i = i + j; //char[i] is now the space before the spin orbit energy      
              j = 1;
              while(data[i+j] != ' ')
                j++;
              memcpy(buff,&data[i],j);
              *zeta = atof(buff);
              return 0;
            }
          }
        }
      }
    }
    i++;
  }
  printf("Configuration with Z = %d, n = %d, l = %d, and electron count of %d not found.",Z,n,ell,vCount);
  return 1;
}
/* 
 * Gives the core-valence Coulomb interactions.
 * Z: atomic number
 * nc: principal quantum number of core shell
 * lc: azimuthal quantum number of core shell
 * nv: principal quantum number of valence shell
 * lv: azimuthal quantum number of valence shell
 * vCount: number of electrons in valence shell
 * Fcv: array storing the coulomb slater integrals to be returned
 * Assumes Fcv is allocated correctly
 */
uint32_t GetInterShellAtomicCoulombInteractions(uint32_t Z, uint32_t nc, uint32_t lc, uint32_t nv, uint32_t lv, uint32_t vCount, double Fcv[MAXNSHELLS][MAXNSHELLS][MAXCOULOMBCV],uint32_t iShellc,uint32_t iShellv)
{
  /*                Z nc lc nv lv vCount [list of Fcv] */
  char data[] = "Z 20 2 2 3 1 0 0.0 3.794 "
                "Z 25 2 2 3 1 5 0.0 6.321 "
                "Z 27 2 2 3 1 7 0.0 7.260 "
                "Z 28 2 2 3 1 7 0.0 8.350 "
                "Z 28 2 2 3 1 8 0.0 7.721 "
                "Z 29 2 2 3 1 9 0.0 8.797 "
                "X";
  char buff[10];
  int i,j,k;
  
  i = 0; j = 0;
  while(data[i] != 'X')
  {
    if(data[i] == 'Z')
    {
      i = i + 2; //char[i] is now first digit of atomic number]
      j = 1;
      while(data[i+j] != ' ')
        j++;
      memcpy(buff,&data[i],j);
      if(atoi(buff) == Z)
      {
        i = i + j + 1; //char[i] is now the first digit of n
        j = 1;
        while(data[i+j] != ' ')
          j++;
        memcpy(buff,&data[i],j);
        if(atoi(buff) == nc)
        {
          i = i + j + 1; //char[i] is now the first digit of ell
          j = 1;
          while(data[i+j] != ' ')
            j++;
          memcpy(buff,&data[i],j);
          if(atoi(buff) == lc)
          {
            i = i + j + 1; //char[i] is now the first digit of vCount
            j = 1;
            while(data[i+j] != ' ')
              j++;
            memcpy(buff,&data[i],j);
            if(atoi(buff) == vCount)
            {
              i = i + j; //char[i] is now the space before first Slater integral
              k = 0;
              while(data[i+1] != 'Z' && data[i+1] != 'X')
              {
                j = 1;
                while(data[i+j] != ' ')
                  j++;
                memcpy(buff,&data[i],j);
                Fcv[iShellc][iShellv][k] = atof(buff);
                k++;
                i = i+j;
              }
              return 0;
            }
          }
        }
      }
    }
    i++;
  }
  printf("Configuration with Z = %d, electron count of %d not found.",Z,vCount);
  return 1;
}

uint32_t GetAtomicInteractions(ModelSettings *O)
{
  int i, j;
  
    
  
  for(i = 0; i < O->NValence; i++)
  {
    /* For each valence shell, get valence-valence and spin-orbit */
    if(GetIntraShellAtomicInteractions(O->AtomicZ, O->ValenceShells[i].n, O->ValenceShells[i].l, O->valenceFillings[i], O->Fvv,i))
    {
      printf("Failed to get intra shell Coulomb interaction energies for Z = %d, n = %d, l = %d\n",O->AtomicZ, O->ValenceShells[i].n, O->ValenceShells[i].l);
      return 1;
    }
    if(GetSpinOrbitEnergy(O->AtomicZ, O->ValenceShells[i].n, O->ValenceShells[i].l, -1, &O->zetav[i]))
    {
      printf("Failed to get spin orbit energy for Z = %d, n = %d, l = %d\n",O->AtomicZ, O->ValenceShells[i].n, O->ValenceShells[i].l);
      return 1;
    }
    for(j = 0; j < O->NCore; j++)
    {
      /* For each core/valence shell combination, get core-valence (if exists) */
      
    }
  }
  
  //for each core shell, get spin orbit
  for(i = 0; i < O->NCore; i++)
  {
    if(GetSpinOrbitEnergy(O->AtomicZ, O->CoreShells[i].n, O->CoreShells[i].n, -1, &O->zetac[i]))
    {
      printf("Failed to get spin orbit energy for Z = %d, n = %d, l = %d\n",O->AtomicZ, O->CoreShells[i].n, O->CoreShells[i].l);
      return 1;
    }
  }
  
  return 0;
} 
uint32_t ScaleAtomicInteractions(ModelSettings *O)
{
  
  return 0;
}



uint32_t RunModelCalculation(ModelSettings MSet)
{
  TimeType Time;
  InitializeTiming(&Time);
  BeginTiming("Total", &Time);
  printf("\n========================================================");
  printf("\n====         Start Model Calculation         ===========");
  printf("\n========================================================\n");
  fflush(stdout);

  uint32_t i, j;
  double Bmag[3], Hex[3], zeta2p, zeta3d;
  double Fdd[3], Fpd[2], Gpd[2];
  double ep, ed, el;
  uint16_t Index2Pdn[3] = {0, 1, 2};
  uint16_t Index2Pup[3] = {3, 4, 5};
  uint16_t Index3Ddn[5] = {6, 7, 8, 9, 10};
  uint16_t Index3Dup[5] = {11, 12, 13, 14, 15};
  uint16_t Index3D[10] = {6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
  
  uint16_t **IndexCoreUp, **IndexCoreDown, *IndexValenceUp, *IndexValenceDown, *IndexValence;
  
  
  HarmonicsExpansion D4hAkma1g, D4hAkmb1g, D4hAkmb2g, D4hAkmeg;
  double HoppingEnergiesa1g[4]={1., 0., 0., 0.};
  double HoppingEnergiesb1g[4]={0., 1., 0., 0.};
  double HoppingEnergiesb2g[4]={0., 0., 1., 0.};
  double HoppingEnergieseg[4] ={0., 0., 0., 1.};
  HarmonicsExpansion Vx, Vy, Vz, Vleft, Vright;
  int Vl[]={1, 1, 1};
  int Vm[]={-1, 0, 1};
  double VxRe[]={sqrt(0.5), 0.0, sqrt(0.5)};
  double VyRe[]={sqrt(0.5), 0.0, -sqrt(0.5)};
  double VzRe[]={0.0, 1.0, 0.0};
  double VlRe[]={1.0, 0.0, 0.0};
  double VrRe[]={0.0, 0.0, 1.0};
  
  OperatorType OppH, OppHXAS, OppSx, OppSy, OppSz, OppSsqr;
  OperatorType OppLx, OppLy, OppLz, OppLsqr;
  OperatorType OppJx, OppJy, OppJz, OppJsqr;
  OperatorType Oppldots2p, Oppldots3d;
  OperatorType OppNthreehalf, OppNfivehalf;
  OperatorType OppUdd, OppUpd, OppNd, OppNdPert, OppNp;
  OperatorType OppNda1g, OppNdb1g, OppNdb2g, OppNdeg;
  OperatorType Opptemp, Opptemp2;
  OperatorType TXASx, TXASy, TXASz, TXASl, TXASr, TRIXS[3];
  double epsilon, zero;
  double energy, energysqr, Esoc, ECoulomb;
  double Sx, Sy, Sz, Ssqr, Lx, Ly, Lz, Lsqr, Jx, Jy, Jz, Jsqr;
  double Nd, Ndsqr, Nthreehalf, Nfivehalf, Nthreehalfsqr, Nfivehalfsqr;
  double Nda1g, Ndb1g, Ndb2g, Ndeg;
  WaveFunctionType psi, Opsi;
  uint16_t *CreateString;
  TriDiagonalMatrixType SpecTri;
  SpectraType Spectra[35];
  MatrixType Density;
  double Hist[30], OccHist[10];
  uint16_t Nelectrons;
  MatrixType Rotated;
  DeterminantRestrictionType Restrictions;
  OperatorType VPert[10];
  TriDiagonalMatrixType gPert[10];

  int numwRes = (int)((MSet.wResRIXSHigh-MSet.wResRIXSLow)/MSet.wResRIXSStep);
  int numwResFY = (int)((MSet.wFYHigh-MSet.wFYLow)/MSet.wFYStep);
  double wRes[numwRes], wResFY[numwResFY];
  SpectraType RIXSSpectra[3][(int)((MSet.wResRIXSHigh-MSet.wResRIXSLow)/MSet.wResRIXSStep)];



  printf("\n========================================================");
  printf("\n==========          set parameters           ===========");
  printf("\n========================================================\n");
  
  
  /* Get atomic coulomb, exchange, and spin orbit parameters */
  GetAtomicInteractions(&MSet);
  
  /* Scale atomic interactions by given factors (i.e. 0.8 usually) */
  /* Also sets the monopole terms, adjusting by Ucv, Uvv etc., if necessary */
  ScaleAtomicInteractions(&MSet);
  
  if(get2p3dAtomicInfo(MSet.AtomicZ,MSet.Nd,Fdd,Fpd,Gpd,&zeta2p,&zeta3d))
  {
    return 1;
  }
  Fdd[1] *= MSet.Kdd; Fdd[2] *= MSet.Kdd; Fpd[1] *= MSet.Kpd; 
  Gpd[0] *= MSet.Kpd; Gpd[1] *= MSet.Kpd;
  Fdd[0] = MSet.Udd + 2.0/63.0 * (Fdd[1] + Fdd[2]);
  Fpd[0] = MSet.Upd + 1.0/15.0 * Gpd[0] + 3.0/70.0 * Gpd[1];
  
  
  if(MSet.Nd == 0)
    ed = MSet.Delta - 0. * MSet.Udd - 6. * MSet.Upd;
  else if(MSet.Nd == 1)
    ed = MSet.Delta - 1. * MSet.Udd - 6. * MSet.Upd;
  else if(MSet.Nd == 2)
    ed = MSet.Delta - 2. * MSet.Udd - 6. * MSet.Upd;
  else if(MSet.Nd == 3)
    ed = MSet.Delta - 3. * MSet.Udd - 6. * MSet.Upd;
  else if(MSet.Nd == 4)
    ed = MSet.Delta - 4. * MSet.Udd - 6. * MSet.Upd;
  else if(MSet.Nd == 5)
    ed = MSet.Delta - 5. * MSet.Udd - 6. * MSet.Upd;
  else if(MSet.Nd == 6)
    ed = MSet.Delta - 6. * MSet.Udd - 6. * MSet.Upd;
  else if(MSet.Nd == 7)
    ed = MSet.Delta - 7. * MSet.Udd - 6. * MSet.Upd;
  else if(MSet.Nd == 8)
    ed = MSet.Delta - 8. * MSet.Udd - 6. * MSet.Upd;
  else if(MSet.Nd == 9)
    ed = MSet.Delta - 9. * MSet.Udd - 6. * MSet.Upd;
  else if(MSet.Nd == 10)
    ed = MSet.Delta - 10. * MSet.Udd - 6. * MSet.Upd;
  else
    ed = 0;
  
  ep = 0;  el = 0.0;

  
  /* Set the number of electrons and fermions according to number of bath sites */
  if(MSet.BathGeometry == 2)
  {
    Nelectrons = 6 + 10*(MSet.NBathV+1);
    MSet.NFermion = 6 + 10 + 10*(MSet.NBathV) + 10 + 10*MSet.NBathC;
    /* Create a starting wavefunction. */
    CreateString=calloc(6+10+10*MSet.NBathV+10, sizeof(uint16_t));
    /*Core electrons*/
    for (i=0; i<MSet.NFermionShift; i++)
      CreateString[i] = i;
    
        //10 electrons on impurity, empty ligand site, full valence bath, empty conduction bath
    for(i=MSet.NFermionShift; i<(MSet.NFermionShift+10);i++)
      CreateString[i] = i;
    for(i=(MSet.NFermionShift+10); i < (10+MSet.NFermionShift+MSet.NBathV*10); i++)
      CreateString[i] = i+10;
    
    //Nd electrons on impurity, 10 empty orbitals, (10-Nd) electrons on ligand site, full valence bath, empty conduction bath)
    /*for(i=MSet.NFermionShift;i<(MSet.Nd+MSet.NFermionShift); i++)
      CreateString[i] = i;
    
    for(i=(MSet.Nd+MSet.NFermionShift);i<(10+MSet.NFermionShift+MSet.NBathV*10); i++)
      CreateString[i] = i+10;
    */
  }
  else
  {
    if(MSet.Perturbation)
    {
      /*Right now just valence chain implemented.  Future: valence&conduction chain, and star geometry*/
      Nelectrons = 6 + MSet.Nd + 10*MSet.NCellV;
      MSet.NFermion = 6 + 10 + 10*MSet.NCellV;
      /* Create a starting wavefunction. */
      CreateString=calloc(6+10+10*MSet.NCellV, sizeof(uint16_t));
      
      
    }
    else
    {
      /*Right now just valence chain implemented.  Future: valence&conduction chain, and star geometry*/
      Nelectrons = 6 + MSet.Nd + 10*MSet.NBathV;
      MSet.NFermion = 6 + 10 + 10*(MSet.NBathV);
      /* Create a starting wavefunction. */
      CreateString=calloc(6+10+10*MSet.NBathV, sizeof(uint16_t));
    }
      
    /*Core electrons*/
    for (i=0; i<MSet.NFermionShift; i++)
      CreateString[i] = i;
  
    /*make a highspin 3d config and fill the ligand orbitals*/
    j = 0;
    for (i=MSet.NFermionShift; i<Nelectrons; i++)
    {
      if((i+j)%2 == 0 && i < (MSet.NFermionShift+MSet.Nd) && j < (10-MSet.Nd) )
      {
        j++; CreateString[i] = i+j;
      }
      else if(i < (MSet.NFermionShift+MSet.Nd))
      {
        CreateString[i] = i+j;
      }
      else
      {
        CreateString[i] = i + (10-MSet.Nd);
      }
    }
  }
  
  
  
  /* Used for number operators */
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

  if(MSet.BathType == 0)
    ImpurityHFromTB(&OppH,MSet);
  else
  {
    if(MSet.Perturbation)
    {
     // ImpurityHFromModelBath(&OppH, MSet.BathGeometry, MSet.BathType, MSet.BandWidth, MSet.NPoles, MSet.NBathV, MSet.NBathV, MSet.NBathC, MSet.EFermi, MSet.kBT, MSet.numIEShifts, 1, MSet.IEShifts, MSet.LEShifts, MSet.Hopping, MSet.NFermion, MSet.NFermionShift);
      ImpurityHFromModelBathPerturbation(&OppH,VPert,gPert, MSet.NCellV,MSet.BathGeometry, MSet.BathType, MSet.BandWidth, MSet.NPoles, MSet.NBathV, MSet.NBathV, MSet.NBathC, MSet.EFermi, MSet.kBT, MSet.numIEShifts, 1, MSet.IEShifts, MSet.LEShifts, MSet.Hopping, MSet.NFermion, MSet.NFermionShift);
    }
    else
      ImpurityHFromModelBath(&OppH, MSet.BathGeometry, MSet.BathType, MSet.BandWidth, MSet.NPoles, MSet.NBathV, MSet.NBathV, MSet.NBathC, MSet.EFermi, MSet.kBT, MSet.numIEShifts, 1, MSet.IEShifts, MSet.LEShifts, MSet.Hopping, MSet.NFermion, MSet.NFermionShift);
  }
  OperatorChop(&OppH, 1000*DBL_EPSILON);
  

  CreateOperatorSx(MSet.NFermion, 5, Index3Dup, Index3Ddn, &OppSx);
  CreateOperatorSy(MSet.NFermion, 5, Index3Dup, Index3Ddn, &OppSy);
  CreateOperatorSz(MSet.NFermion, 5, Index3Dup, Index3Ddn, &OppSz);
  CreateOperatorSsqr(MSet.NFermion, 5, Index3Dup, Index3Ddn, &OppSsqr);
  
  CreateOperatorLx(MSet.NFermion, 2, Index3Dup, Index3Ddn, &OppLx);
  CreateOperatorLy(MSet.NFermion, 2, Index3Dup, Index3Ddn, &OppLy);
  CreateOperatorLz(MSet.NFermion, 2, Index3Dup, Index3Ddn, &OppLz);
  CreateOperatorLsqr(MSet.NFermion, 2, Index3Dup, Index3Ddn, &OppLsqr);
  
  CreateOperatorJx(MSet.NFermion, 2, Index3Dup, Index3Ddn, &OppJx);
  CreateOperatorJy(MSet.NFermion, 2, Index3Dup, Index3Ddn, &OppJy);
  CreateOperatorJz(MSet.NFermion, 2, Index3Dup, Index3Ddn, &OppJz);
  CreateOperatorJsqr(MSet.NFermion, 2, Index3Dup, Index3Ddn, &OppJsqr);

  
 // if(MSet.Perturbation)
 //   CreateOperatorldots(MSet.NFermion-10*(MSet.NBathV-MSet.NCellV), 1, Index2Pup, Index2Pdn, &Oppldots2p);
 // else
    CreateOperatorldots(MSet.NFermion, 1, Index2Pup, Index2Pdn, &Oppldots2p);
  
  CreateOperatorldots(MSet.NFermion, 2, Index3Dup, Index3Ddn, &Oppldots3d);
  
  CreateOperatorCoulombSingleShell(MSet.NFermion, 2, Index3Dup, Index3Ddn, Fdd, &OppUdd);
  
 // if(MSet.Perturbation)
  //{
  //  CreateOperatorCoulombTwoShellsConstantOccupation(MSet.NFermion-10*(MSet.NBathV-MSet.NCellV), 1, 2, Index2Pup, Index2Pdn, Index3Dup, Index3Ddn, Fpd, Gpd, &OppUpd);
  //  CreateOperatorNumberShell(MSet.NFermion-10*(MSet.NBathV-MSet.NCellV), 10, Index3D, &OppNdPert);
 // }
  //else
 // {
    CreateOperatorCoulombTwoShellsConstantOccupation(MSet.NFermion, 1, 2, Index2Pup, Index2Pdn, Index3Dup, Index3Ddn, Fpd, Gpd, &OppUpd);
  //  CreateOperatorNumberShell(MSet.NFermion, 10, Index3D, &OppNdPert);
 // }

  
  CreateOperatorNumberShell(MSet.NFermion, 10, Index3D, &OppNd);
  
  
  
  ////CreateOperatorNumberShell(NFermion, 10, IndexLD, &OppNL);
  CreateOperatorCF(MSet.NFermion, 2, 2, Index3Dup, Index3Ddn, Index3Dup, Index3Ddn, D4hAkma1g, &OppNda1g);
  CreateOperatorCF(MSet.NFermion, 2, 2, Index3Dup, Index3Ddn, Index3Dup, Index3Ddn, D4hAkmb1g, &OppNdb1g);
  CreateOperatorCF(MSet.NFermion, 2, 2, Index3Dup, Index3Ddn, Index3Dup, Index3Ddn, D4hAkmb2g, &OppNdb2g);
  CreateOperatorCF(MSet.NFermion, 2, 2, Index3Dup, Index3Ddn, Index3Dup, Index3Ddn, D4hAkmeg,  &OppNdeg );
  
  OperatorCopy(&OppNthreehalf, Oppldots3d); /* <ls> = -3/2n(j=3/2)+n(j=5/2) */
  OperatorAdd(&OppNthreehalf, OppNd, -1.0, 0.0);
  OperatorScale(&OppNthreehalf, -0.4, 0.0);
  OperatorCopy(&OppNfivehalf, OppNd);
  OperatorAdd(&OppNfivehalf, OppNthreehalf, -1.0, 0.0);
  
    /* create the transition operator for XAS */
  CreateOperatorCF(MSet.NFermion, 2, 1, Index3Dup, Index3Ddn, Index2Pup, Index2Pdn, Vx, &TXASx);
  CreateOperatorCF(MSet.NFermion, 2, 1, Index3Dup, Index3Ddn, Index2Pup, Index2Pdn, Vy, &TXASy);
  CreateOperatorCF(MSet.NFermion, 2, 1, Index3Dup, Index3Ddn, Index2Pup, Index2Pdn, Vz, &TXASz);
  CreateOperatorCF(MSet.NFermion, 2, 1, Index3Dup, Index3Ddn, Index2Pup, Index2Pdn, Vleft, &TXASl);
  CreateOperatorCF(MSet.NFermion, 2, 1, Index3Dup, Index3Ddn, Index2Pup, Index2Pdn, Vright, &TXASr);

  /*Create transition operators for RIXS */
  CreateOperatorCF(MSet.NFermion, 1, 2, Index2Pup, Index2Pdn, Index3Dup, Index3Ddn, Vx, &TRIXS[0]);
  CreateOperatorCF(MSet.NFermion, 1, 2, Index2Pup, Index2Pdn, Index3Dup, Index3Ddn, Vy, &TRIXS[1]);
  CreateOperatorCF(MSet.NFermion, 1, 2, Index2Pup, Index2Pdn, Index3Dup, Index3Ddn, Vz, &TRIXS[2]);
  
  

    /* Define rotation matrix */
  printf("\n========================================================");
  printf("\n====  Define rotation matrix to 'real' orbitals   ======");
  printf("\n========================================================\n");
  Rotated.QComplex = 0;
  Rotated.Ni = MSet.NFermion;
  Rotated.Nj = MSet.NFermion;
  InitMatrix(&Rotated);
  for (i=0; i<Rotated.Ni; i++)
  {
    for (j=0; j<Rotated.Nj; j++)
    {
      Rotated.VRe[i][j] = 0.0;
    }
  }
  Rotated.VRe[0][0] =-sqrt(0.5);    Rotated.VRe[0][2] = sqrt(0.5);
  Rotated.VRe[3][3] =-sqrt(0.5);    Rotated.VRe[3][5] = sqrt(0.5);

  Rotated.VRe[1][1] = 1.0;
  Rotated.VRe[4][4] = 1.0;

  Rotated.VRe[2][0] =-sqrt(0.5);    Rotated.VRe[2][2] =-sqrt(0.5);
  Rotated.VRe[5][3] =-sqrt(0.5);    Rotated.VRe[5][5] =-sqrt(0.5);
  
  

  if(MSet.BathGeometry == 2)
    j = MSet.NBathV + 1 + 1 + MSet.NBathC; //impurity, ligand, valence chain, conduction chain
  else if(MSet.BathGeometry == 1 || MSet.BathGeometry == 0)
  {
    if(MSet.Perturbation)
      j = MSet.NCellV + 1;
    else
      j = MSet.NBathV + 1; //just valence chain and impurity for now.  conduction chain to be implemented
  }
  else
  {
    printf("BathGeometry %d not supported\n",MSet.BathGeometry);
    exit(1);
  }
  for (i=0; i<j; i++)
  {
    Rotated.VRe[ 6+i*10][ 6+i*10] =-sqrt(0.5);   Rotated.VRe[ 6+i*10][10+i*10] = sqrt(0.5);
    Rotated.VRe[ 11+i*10][ 11+i*10] =-sqrt(0.5);   Rotated.VRe[ 11+i*10][15+i*10] = sqrt(0.5);
    
    Rotated.VRe[ 7+i*10][ 7+i*10] =-sqrt(0.5);   Rotated.VRe[ 7+i*10][9+i*10] = sqrt(0.5);
    Rotated.VRe[ 12+i*10][ 12+i*10] =-sqrt(0.5);   Rotated.VRe[ 12+i*10][14+i*10] = sqrt(0.5);
    
    Rotated.VRe[8+i*10][8+i*10] = 1.0;
    Rotated.VRe[13+i*10][13+i*10] = 1.0;
    
    Rotated.VRe[9+i*10][ 7+i*10] =-sqrt(0.5);   Rotated.VRe[9+i*10][9+i*10] =-sqrt(0.5);
    Rotated.VRe[14+i*10][ 12+i*10] =-sqrt(0.5);   Rotated.VRe[14+i*10][14+i*10] =-sqrt(0.5);
    
    Rotated.VRe[10+i*10][ 6+i*10] = sqrt(0.5);   Rotated.VRe[10+i*10][10+i*10] = sqrt(0.5);
    Rotated.VRe[15+i*10][ 11+i*10] = sqrt(0.5);   Rotated.VRe[15+i*10][15+i*10] = sqrt(0.5);
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
  
  /* Rotate spectral operators */
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
  
  for(i = 0; i < 3; i++)
  {
    RotateOperatorNormalOrder(&Opptemp, TRIXS[i], Rotated, &QDetNotOrderedFcFaBcBa); OperatorChop(&Opptemp, 100*DBL_EPSILON);
    FreeOperator(&TRIXS[i]); OperatorCopy(&TRIXS[i], Opptemp); FreeOperator(&Opptemp);
  }
  
  /*Test that Udd is Hermitian*/
  OperatorConjugateTranspose(OppUdd, &Opptemp);
  OperatorAdd(&Opptemp, OppUdd, -1.0, 0.0);
  InitOperatorCopySettings(&Opptemp2, Opptemp);
  NormalOrderOperator(Opptemp, &Opptemp2, &QDetNotOrderedFcFaBaBc); // different orderings possible, all need to be the same to get annihilation
  OperatorChop(&Opptemp2, 1000*DBL_EPSILON);
  printf("\nPrinting Opptemp = Udd - Udd^+ \n");
  PrintOperatorReadable(Opptemp2);
  FreeOperator(&Opptemp);
  FreeOperator(&Opptemp2);

  /*                                          */
  /* sum operators to create the Hamiltonian */
                                           
  OperatorAdd(&OppH, OppUdd, 1.0, 0.0);
  OperatorAdd(&OppH, Oppldots3d, zeta3d, 0.0);  
  
  OperatorAdd(&OppH, OppSx, 2*MSet.BMag[0]+MSet.HExch[0], 0.0);
  OperatorAdd(&OppH, OppSy, 2*MSet.BMag[1]+MSet.HExch[1], 0.0);
  OperatorAdd(&OppH, OppSz, 2*MSet.BMag[2]+MSet.HExch[2], 0.0);
  OperatorAdd(&OppH, OppLx, MSet.BMag[0], 0.0);
  OperatorAdd(&OppH, OppLy, MSet.BMag[1], 0.0);
  OperatorAdd(&OppH, OppLz, MSet.BMag[2], 0.0); 
  
  OperatorChop(&OppH, DBL_EPSILON*1000);

  //if(!MSet.Perturbation)
    OperatorCopy(&OppHXAS,OppH);
  
  OperatorAdd(&OppH, OppNd, ed+6.*MSet.Upd, 0.0);
  

  OperatorAdd(&OppHXAS, OppNd, ed, 0.0);
  OperatorAdd(&OppHXAS, OppUpd, 1.0, 0.0);
  OperatorAdd(&OppHXAS, Oppldots2p, zeta2p, 0.0);
  

  
  /* Test that H is Hermitian */
  OperatorConjugateTranspose(OppH, &Opptemp);
  OperatorAdd(&Opptemp, OppH, -1.0, 0.0);
  InitOperatorCopySettings(&Opptemp2, Opptemp);
  NormalOrderOperator(Opptemp, &Opptemp2, &QDetNotOrderedFcFaBaBc); // different orderings possible, all need to be the same to get annihilation
  OperatorChop(&Opptemp2, 1000*DBL_EPSILON);
  printf("\nPrinting Opptemp = H - H^+ \n");
  PrintOperatorReadable(Opptemp2);
  FreeOperator(&Opptemp);
  FreeOperator(&Opptemp2);
  
  EndTiming("define operators", &Time);
  
  
  
  
  printf("\n========================================================");
  printf("\n==========        find ground state          ===========");
  printf("\n========================================================\n");
  BeginTiming("find groundstate", &Time);
  
  printf("\n\n\n\n");
  printf("Wavefunction string: ");
  for(i = 0; i < Nelectrons; i++)
    printf("%d ", CreateString[i]);
  printf("\n\n\n");

  CreateWaveFunctionFromString(&psi, CreateString, Nelectrons, MSet.NFermion, 0, MSet.NBitsKey);
//  InitNaturalImpurityWavefunctionRandom(&psi, MSet.NBitsKey,MSet.NFermion,MSet.NFermionShift,10,MSet.NBathV,MSet.Nd,MSet.Nd);

  //PrintWaveFunctionFullReadable(psi);
  //  WaveFunctionRandomize(&psi);

  BeginTiming("LanczosGroundState", &Time);
  if(MSet.BathGeometry == 2)
  {
    // In order to help convergence, first require 10 electrons spread between impurity and separate ligand site.
    Restrictions.N=1;
    Restrictions.Min=malloc(Restrictions.N*sizeof(uint32_t));
    Restrictions.Max=malloc(Restrictions.N*sizeof(uint32_t));
    Restrictions.Determinant=malloc(Restrictions.N*sizeof(uint8_t*));
    for(i=0;i<Restrictions.N;i++)
    {
      Restrictions.Determinant[i]=calloc(psi.NWord,sizeof(uint8_t));
    }
    Restrictions.Min[0]=10;
    Restrictions.Max[0]=10;
    Restrictions.Determinant[0][0] = 0x03;
    Restrictions.Determinant[0][1] = 0xFF;
    Restrictions.Determinant[0][2] = 0xFF;
    Restrictions.Determinant[0][3] = 0xC0;

    
    // start with low accuracy
  //  epsilon=0.001;
   // zero=100*DBL_EPSILON;
    
   // LanczosGroundStateRestricted(OppH, &psi, Restrictions, &epsilon, zero, 0, &psi, 0);
   // PrintWaveFunctionDensity(psi);
    
    // increase accuracy
    epsilon=MSet.epsilon;
    zero=100*DBL_EPSILON;
    
    LanczosGroundStateRestricted(OppH, &psi, Restrictions, &epsilon, zero, 0, &psi, 0);
    PrintWaveFunctionDensity(psi);

    //remove restrictions (shouldn't change things)
    LanczosGroundState(OppH, &psi, &epsilon, zero, 0, &psi, 0);
    PrintWaveFunctionDensity(psi);
  }
  else
  {
    epsilon=MSet.epsilon;
    zero=100*DBL_EPSILON;
    LanczosGroundState(OppH, &psi, &epsilon, zero, 0, &psi, 0);

    WaveFunctionRandomize(&psi);

    LanczosGroundState(OppH, &psi, &epsilon, zero, 0, &psi, 0);
  
    WaveFunctionRandomize(&psi);
  
    LanczosGroundState(OppH, &psi, &epsilon, zero, 0, &psi, 0);
    LanczosGroundState(OppH, &psi, &epsilon, zero, 0, &psi, 0);
    LanczosGroundState(OppH, &psi, &epsilon, zero, 0, &psi, 0);

  }

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
  printf("\n==========   Groundstate expectation values   ===========");
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

  if(MSet.DoXAS)
  {
    printf("\n========================================================");
    printf("\n==========                XAS                ===========");
    printf("\n========================================================\n");
    BeginTiming("XAS", &Time);

    SpecTri.N = MSet.NKrylovXAS;
    InitTriDiagonalMatrix(&SpecTri);
     
    for(i=26;i<31;i++)
    {
      Spectra[i].N=(int)((MSet.wXASHigh-MSet.wXASLow)/MSet.wXASStep);
      Spectra[i].Emin=MSet.wXASLow;
      Spectra[i].Emax=MSet.wXASHigh;
      Spectra[i].dE = (Spectra[i].Emax - Spectra[i].Emin) / Spectra[i].N;
      Spectra[i].Gamma = MSet.LXAS;
      SpectrumAllocate(Spectra[i], &Spectra[i]);
    }
    
    BeginTiming("spectrum", &Time);
    OperatorPsi(TXASx, psi, &Opsi);
    if(MSet.Perturbation)
      LanczosTriDiagonalizeDysonLowMemBLAS(OppHXAS, Opsi, epsilon,MSet.NKrylovXAS, VPert, gPert, 10, &SpecTri);
    else
      LanczosTriDiagonalize(OppHXAS, Opsi, epsilon, &SpecTri,  Spectra[26].Gamma);
    TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[26]);
    FilePrintPlainTriDiagonalMatrix(SpecTri, "XASTri_x", "w");
    EndTiming("spectrum", &Time);
    
    BeginTiming("spectrum", &Time);
    OperatorPsi(TXASy, psi, &Opsi);
    if(MSet.Perturbation)
      LanczosTriDiagonalizeDysonLowMemBLAS(OppHXAS, Opsi, epsilon,MSet.NKrylovXAS, VPert, gPert, 10, &SpecTri);
    else
      LanczosTriDiagonalize(OppHXAS, Opsi, epsilon, &SpecTri,  Spectra[27].Gamma);
    TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[27]);
    FilePrintPlainTriDiagonalMatrix(SpecTri, "XASTri_y", "w");
    EndTiming("spectrum", &Time);
    
    BeginTiming("spectrum", &Time);
    OperatorPsi(TXASz, psi, &Opsi);
    if(MSet.Perturbation)
      LanczosTriDiagonalizeDysonLowMemBLAS(OppHXAS, Opsi, epsilon,MSet.NKrylovXAS, VPert, gPert, 10, &SpecTri);
    else
      LanczosTriDiagonalize(OppHXAS, Opsi, epsilon, &SpecTri,  Spectra[28].Gamma);
    TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[28]);
    FilePrintPlainTriDiagonalMatrix(SpecTri, "XASTri_z", "w");
    EndTiming("spectrum", &Time);
    /*
    BeginTiming("spectrum", &Time);
    OperatorPsi(TXASl, psi, &Opsi);
    LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[29].Gamma);
    TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[29]);
    EndTiming("spectrum", &Time);
    
    BeginTiming("spectrum", &Time);
    OperatorPsi(TXASr, psi, &Opsi);
    LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri,  Spectra[30].Gamma);
    TriDiagonalMatrixToSpectrum(energy, SpecTri, &Spectra[30]);
    EndTiming("spectrum", &Time);*/
    BroadenG(&Spectra[26], MSet.GXAS);
    BroadenG(&Spectra[27], MSet.GXAS);
    BroadenG(&Spectra[28], MSet.GXAS);
    
    
    
    PrintModelSettings(MSet,"XAS_Spectra");
    FilePrintSpectra(&Spectra[26], 3, "XAS_Spectra","a");

    /*Print a flipped, isotropic spectrum*/
    SpectraScalarMultiply(&Spectra[26], -1.0, 0.0);
    AddMultipliedSpectra(&Spectra[26], Spectra[27], -1.0, 0.0);
    AddMultipliedSpectra(&Spectra[26], Spectra[28], -1.0, 0.0);
    FilePrintSpectra(&Spectra[26], 1, "XAS_SpectrumISO","w");

    
    EndTiming("XAS", &Time);
    PrintTiming(Time);
  }
  if(MSet.DoFY)
  {
    printf("\n========================================================");
    printf("\n==========       Fluorescence Yield          ===========");
    printf("\n========================================================\n");
    BeginTiming("FY", &Time);
    
    for(i = 0; i < numwResFY; i++)
      wResFY[i] = MSet.wFYLow + ((double)i)*MSet.wFYStep;
  
    /* Just doing xx for now */
    BeginTiming("spectrum", &Time);
    FluorescenceYield(OppHXAS,TXASx,TRIXS,1,psi,energy,epsilon,numwResFY,&Spectra[0],wResFY,MSet.LFY,MSet.GFY,MSet.NKrylovFY);
    EndTiming("spectrum", &Time);

    FilePrintSpectra(&Spectra[0], 1, "Fluorescence_Yield","w");

    EndTiming("FY", &Time);
    PrintTiming(Time);
  }
  if(MSet.DoRIXS)
  {
    printf("\n========================================================");
    printf("\n==========               RIXS                ===========");
    printf("\n========================================================\n");
    BeginTiming("RIXS", &Time);

    for(i = 0; i < numwRes; i++)
    {
      wRes[i] = MSet.wResRIXSLow + ((double)i)*MSet.wResRIXSStep;
    }  
    for(i=0;i<3;i++)
    {
      for(j=0; j < numwRes; j++)
      {
        RIXSSpectra[i][j].N=(int)((MSet.wRIXSHigh - MSet.wRIXSLow)/MSet.wRIXSStep);
        RIXSSpectra[i][j].Emin=MSet.wRIXSLow;
        RIXSSpectra[i][j].Emax=MSet.wRIXSHigh;
        RIXSSpectra[i][j].dE = (RIXSSpectra[i][j].Emax - RIXSSpectra[i][j].Emin) / RIXSSpectra[i][j].N;
        RIXSSpectra[i][j].Gamma = MSet.LRIXS;
        SpectrumAllocate(RIXSSpectra[i][j], &RIXSSpectra[i][j]);
      }
    }

    /* Just set to do xx for testing.  To do xx, xy, xz, set the parameter to 3 in the call, and
       uncomment the file printing lines
     
       Update: now set to do xy by giving address of TRIXS[1]*/
    BeginTiming("spectrum", &Time);
    ResonantSpectra(OppHXAS,OppH,TXASx,&TRIXS[1],1,psi,energy,epsilon,numwRes,RIXSSpectra,wRes,MSet.LResRIXS,MSet.LRIXS,MSet.GRIXS,MSet.NKrylovResRIXS,MSet.NKrylovRIXS);
    EndTiming("spectrum", &Time);

    FilePrintRIXSSpectrum(numwRes,RIXSSpectra,wRes, 1,0, "RIXS_xy_Real","w");
    FilePrintRIXSSpectrum(numwRes,RIXSSpectra,wRes, 0,0, "RIXS_xy_Imag","w");
    
    
//    FilePrintRIXSSpectrum(numwRes,RIXSSpectra,wRes, 1,0, "RIXS_xx_Real","w");
//    FilePrintRIXSSpectrum(numwRes,RIXSSpectra,wRes, 0,0, "RIXS_xx_Imag","w");
//    FilePrintRIXSSpectrum(numwRes,RIXSSpectra,wRes, 1,1, "RIXS_xy_Real","w");
//    FilePrintRIXSSpectrum(numwRes,RIXSSpectra,wRes, 0,1, "RIXS_xy_Imag","w");
//    FilePrintRIXSSpectrum(numwRes,RIXSSpectra,wRes, 1,2, "RIXS_xz_Real","w");
//    FilePrintRIXSSpectrum(numwRes,RIXSSpectra,wRes, 0,2, "RIXS_xz_Imag","w");

    EndTiming("RIXS", &Time);
    PrintTiming(Time);
  }
  if(MSet.DoPES)
  {
    
  }
  
  EndTiming("Total", &Time);
  PrintTiming(Time);
  return 0;
  
  
}

