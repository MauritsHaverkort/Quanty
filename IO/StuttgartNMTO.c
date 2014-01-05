#define StuttgartCTRLCategories {"HEADER","VERS","IO","DIM","STRUC","OPTIONS","CLASS","SITE","SYMGRP","SCALE","STR","BZ","EWALD","DOS","SYML","START","HARTREE","PLOT","CHARGE","FINDES","SCELL","RHOFIT","NMTO","LOCALCS"}
#define NStuttgartCTRLCategories 24


/*  This function replaces the character "in" by the character "out" in the string "String"
 *  Your string has end with a null terminator, otherwise all bytes of code that match "in" will be replaced by "out"
 */
int StringReplace(char * String, char in, char out)
{
  int i;
  for(i=0;i<strlen(String);i++)
  {
    if(String[i]==in)
    {
      String[i]=out;
    }
  }
  return 0;
}

/*
 *  scans the next value category or token from a CTRL file
 *  the character "=" is treated the same as the space (" ") and enter ("\n") characters
 */
int fscanfCTRL(FILE * stream, char * str, Bool * token)
{
  char strlong[1024];
  int i,imax,retval;
  long irest;
  *token=False;
  retval=fscanf(stream, "%s", strlong);
  if(retval==1)
  {
    imax=(int)strlen(strlong);
    for(i=0;i<imax;i++)
    {
      if(strlong[i]=='=')
      {
        *token=True;
        strlong[i]='\0';
        irest=(long)(i-imax+1);
        fseek(stream, irest, SEEK_CUR);
      }
    }
    strcpy(str, strlong);
  }
  return retval;
}

/*
 *  scans the file char by char until the next character
 *  is not equal to an space, enter or equal sign
 */
int fscanfRemoveEqualSign(FILE * stream)
{
  int c;
  int retval=0;
  c=fgetc(stream);
  while((c==(int)' ')||(c==(int)'=')||(c==(int)'\n'))
  {
    if(c==(int)' ')retval=1;
    c=fgetc(stream);
  }
  fseek(stream, -1, SEEK_CUR);
  return retval;
}

uint32_t StuttgartCategoryToIndex(char *str)
{
  const char * const Categories[] = StuttgartCTRLCategories;
  uint32_t i;
  for(i=0;i<NStuttgartCTRLCategories;i++)
  {
    if(strcmp(str,Categories[i])==0) return i+1;
  }
  return 0;
}

uint32_t QStuttgartCategory(char *str)
{
  return StuttgartCategoryToIndex(str);
}

uint32_t IndexToStuttgartCategory(uint32_t i, char *str)
{
  if((i==0)||(i>NStuttgartCTRLCategories))
  {
    strcpy(str, "UNDEFINED");
    return 1;
  }
  const char * const Categories[] = StuttgartCTRLCategories;
  strcpy(str, Categories[i-1]);
  return 0;
}

/*
 *  set standard values for the structure
 */
uint32_t InitializeStuttgartCTRL(StuttgartCTRL * O)
{
  O->SYMGRP.NGEN=-1;
  O->CLASS.N=-1;
  O->SITE.N=-1;
  O->SYML.N=-1;
  O->LOCALCS.Number_of_rotated_atoms=-1;
  O->NMTO.NENU=-1;
  O->HEADER[0]='\0';
  O->VERS[0]='\0';
  
  return 0;
}

uint32_t StreamReadStuttgartCTRLBool(FILE * stream, uint32_t N, Bool * out)
{
  char inbool;
  uint32_t i;
  for(i=0;i<N;i++)
  {
    fscanfRemoveEqualSign(stream);
    if(fscanf(stream, "%c",&inbool)!=1)
    {
      printf("Error while reading boolean from CTRL file\n");
      return 1;
    }
    if(inbool=='T'){out[i]=True;}
    else if(inbool=='F'){out[i]=False;}
    else
    {
      printf("Error while reading boolean from CTRL file, expected T or F, found %c\n",inbool);
    }
  }
  return 0;
}

uint32_t StreamReadStuttgartCTRLint(FILE * stream, uint32_t N, int * out)
{
  uint32_t i;
  for(i=0;i<N;i++)
  {
    fscanfRemoveEqualSign(stream);
    if(fscanf(stream, "%d",&(out[i]))!=1)
    {
      printf("Error while reading integer from CTRL file\n");
      return 1;
    }
  }
  return 0;
}

uint32_t StreamReadStuttgartCTRLdouble(FILE * stream, uint32_t N, double * out)
{
  uint32_t i;
  for(i=0;i<N;i++)
  {
    fscanfRemoveEqualSign(stream);
    if(fscanf(stream, "%lE",&(out[i]))!=1)
    {
      printf("Error while reading double from CTRL file\n");
      return 1;
    }
  }
  return 0;
}

uint32_t StreamReadStuttgartCTRLchar(FILE * stream, char * out)
{
  fscanfRemoveEqualSign(stream);
  if(fscanf(stream, "%s", out )!=1)
  {
    printf("Error while reading char string from CTRL file\n");
    return 1;
  }
  return 0;
}

/*
 *   Reads a unit cell definition from a CTRL file of the stuttgart NMTO/LMTO code
 */
uint32_t StreamReadStuttgartCTRL(StuttgartCTRL * O, FILE * stream)
{
  const char * const Categories[] = StuttgartCTRLCategories;
  char str[128];
  int Continue=1, ContinueCategory, iCategory, i;
  Bool ForSureAToken;
  int Nlocalcs=-1;
  InitializeStuttgartCTRL(O);
  if(fscanfCTRL(stream, str, &ForSureAToken)!=1) return 0;
  while(Continue)
  {
    if(!QStuttgartCategory(str))
    {
      printf("Expected a Category instead of %s\n",str);
      return 1;
    }
    iCategory=StuttgartCategoryToIndex(str);
    ContinueCategory=1;
    while(ContinueCategory)
    {
      if(fscanfCTRL(stream, str, &ForSureAToken)!=1) return 0;
      if((ForSureAToken==False)&&(QStuttgartCategory(str)))
      {
        ContinueCategory=0;
        break;
      }
      switch (iCategory)
      {
        case 1: /* HEADER */
          strncat(O->HEADER, str, sizeof(O->HEADER)-strlen(O->HEADER)-1);
          strncat(O->HEADER, " ", sizeof(O->HEADER)-strlen(O->HEADER)-1);
        break;
        case 2: /* VERS */
          strncat(O->VERS, str, sizeof(O->HEADER)-strlen(O->HEADER)-1);
          strncat(O->VERS, " ", sizeof(O->HEADER)-strlen(O->HEADER)-1);
        break;
        case 3: /* IO */
          if(strcmp(str, "HELP")==0)
          {
            StreamReadStuttgartCTRLBool(stream, 1, &(O->IO.HELP));
          }
          else if(strcmp(str, "VERBOS")==0)
          {
            StreamReadStuttgartCTRLint(stream, 1, &(O->IO.VERBOS));
          }
          else if(strcmp(str, "WKP")==0)
          {
            StreamReadStuttgartCTRLBool(stream, 1, &(O->IO.WKP));
          }
          else if(strcmp(str, "IACTIV")==0)
          {
            StreamReadStuttgartCTRLBool(stream, 1, &(O->IO.IACTIV));
          }
          else if(strcmp(str, "ERRTOL")==0)
          {
            StreamReadStuttgartCTRLint(stream, 1, &(O->IO.ERRTOL));
          }
          else if(strcmp(str, "OUTPUT")==0)
          {
            StreamReadStuttgartCTRLchar(stream, O->IO.OUTPUT);
          }
          else if(strcmp(str, "ERR")==0)
          {
            StreamReadStuttgartCTRLchar(stream, O->IO.ERR);
          }
          else
          {
            printf("In Category %s undefined token %s\n",Categories[iCategory-1],str);
          }
        break;
        case 4: /* DIM */
          if(strcmp(str, "NBAS")==0)
          {
            StreamReadStuttgartCTRLint(stream, 1, &(O->DIM.NBAS));
          }
          else if(strcmp(str, "NCLASS")==0)
          {
            StreamReadStuttgartCTRLint(stream, 1, &(O->DIM.NCLASS));
          }
          else if(strcmp(str, "NL")==0)
          {
            StreamReadStuttgartCTRLint(stream, 1, &(O->DIM.NL));
          }
          else if(strcmp(str, "LDIM")==0)
          {
            StreamReadStuttgartCTRLint(stream, 1, &(O->DIM.LDIM));
          }
          else if(strcmp(str, "IDIM")==0)
          {
            StreamReadStuttgartCTRLint(stream, 1, &(O->DIM.IDIM));
          }
          else if(strcmp(str, "NSYMOP")==0)
          {
            StreamReadStuttgartCTRLint(stream, 1, &(O->DIM.NSYMOP));
          }
          else if(strcmp(str, "NKP")==0)
          {
            StreamReadStuttgartCTRLint(stream, 1, &(O->DIM.NKP));
          }
          else
          {
            printf("In Category %s undefined token %s\n",Categories[iCategory-1],str);
          }
        break;
        case 5: /* STRUC */
          if(strcmp(str, "ALAT")==0)
          {
            StreamReadStuttgartCTRLdouble(stream, 1, &(O->STRUC.ALAT));
          }
          else if(strcmp(str, "PLAT")==0)
          {
            StreamReadStuttgartCTRLdouble(stream, 9, &(O->STRUC.PLAT[0][0]));
          }
          else if(strcmp(str, "FIXLAT")==0)
          {
            StreamReadStuttgartCTRLBool(stream, 1, &(O->STRUC.FIXLAT));
          }
          else
          {
            printf("In Category %s undefined token %s\n",Categories[iCategory-1],str);
          }
        break;
        case 6: /* OPTIONS */
          if(strcmp(str, "REL")==0)
          {
            StreamReadStuttgartCTRLBool(stream, 1, &(O->OPTIONS.REL));
          }
          else if(strcmp(str, "CCOR")==0)
          {
            StreamReadStuttgartCTRLBool(stream, 1, &(O->OPTIONS.CCOR));
          }
          else if(strcmp(str, "NONLOC")==0)
          {
            StreamReadStuttgartCTRLBool(stream, 1, &(O->OPTIONS.NONLOC));
          }
          else if(strcmp(str, "CORDRD")==0)
          {
            StreamReadStuttgartCTRLBool(stream, 1, &(O->OPTIONS.CORDRD));
          }
          else if(strcmp(str, "CARTESIAN")==0)
          {
            StreamReadStuttgartCTRLBool(stream, 1, &(O->OPTIONS.CARTESIAN));
          }
          else if(strcmp(str, "CHARGE")==0)
          {
            StreamReadStuttgartCTRLBool(stream, 1, &(O->OPTIONS.CHARGE));
          }
          else if(strcmp(str, "FATBAND")==0)
          {
            StreamReadStuttgartCTRLBool(stream, 1, &(O->OPTIONS.FATBAND));
          }
          else if(strcmp(str, "AFM")==0)
          {
            StreamReadStuttgartCTRLBool(stream, 1, &(O->OPTIONS.AFM));
          }
          else if(strcmp(str, "FS")==0)
          {
            StreamReadStuttgartCTRLBool(stream, 1, &(O->OPTIONS.FS));
          }
          else if(strcmp(str, "OVLCOR")==0)
          {
            StreamReadStuttgartCTRLBool(stream, 1, &(O->OPTIONS.OVLCOR));
          }
          else if(strcmp(str, "WRIBAS")==0)
          {
            StreamReadStuttgartCTRLBool(stream, 1, &(O->OPTIONS.WRIBAS));
          }
          else if(strcmp(str, "SEWALD")==0)
          {
            StreamReadStuttgartCTRLBool(stream, 1, &(O->OPTIONS.SEWALD));
          }
          else if(strcmp(str, "LOCALCS")==0)
          {
            StreamReadStuttgartCTRLBool(stream, 1, &(O->OPTIONS.LOCALCS));
          }
          else if(strcmp(str, "NRXC")==0)
          {
            StreamReadStuttgartCTRLint(stream, 1, &(O->OPTIONS.NRXC));
          }
          else if(strcmp(str, "NITATOM")==0)
          {
            StreamReadStuttgartCTRLint(stream, 1, &(O->OPTIONS.NITATOM));
          }
          else if(strcmp(str, "NSPIN")==0)
          {
            StreamReadStuttgartCTRLint(stream, 1, &(O->OPTIONS.NSPIN));
          }
          else if(strcmp(str, "NRMIX")==0)
          {
            StreamReadStuttgartCTRLint(stream, 1, &(O->OPTIONS.NRMIX));
          }
          else if(strcmp(str, "Q")==0)
          {
            StreamReadStuttgartCTRLchar(stream, O->OPTIONS.Q);
          }
          else
          {
            printf("In Category %s undefined token %s\n",Categories[iCategory-1],str);
          }
        break;
        case 7: /* CLASS */
          if(strcmp(str, "ATOM")==0)
          {
            if(O->CLASS.N==-1)
            {
              O->CLASS.NEWNAM=malloc(64*sizeof(Bool));
              O->CLASS.ATOM=malloc(64*sizeof(char*));
              O->CLASS.Z=malloc(64*sizeof(double));
              O->CLASS.R=malloc(64*sizeof(double));
              O->CLASS.LMX=malloc(64*sizeof(int));
              O->CLASS.IDXDN=malloc(64*sizeof(int*));
              O->CLASS.CONF=malloc(64*sizeof(int*));
              O->CLASS.IDMOD=malloc(64*sizeof(int*));
              for(i=0;i<64;i++)
              {
                O->CLASS.LMX[i]=-1;
                O->CLASS.ATOM[i]=malloc(16*sizeof(char));
                O->CLASS.IDXDN[i]=malloc(4*sizeof(int));
                O->CLASS.CONF[i]=malloc(4*sizeof(int));
                O->CLASS.IDMOD[i]=malloc(4*sizeof(int));
              }
              O->CLASS.N=0;
            }
            O->CLASS.N++;
            if((O->CLASS.N%64)==0)
            {
              O->CLASS.NEWNAM=realloc(O->CLASS.NEWNAM, (O->CLASS.N+64)*sizeof(Bool));
              O->CLASS.ATOM=realloc(O->CLASS.ATOM, (O->CLASS.N+64)*sizeof(char*));
              O->CLASS.Z=realloc(O->CLASS.Z, (O->CLASS.N+64)*sizeof(double));
              O->CLASS.R=realloc(O->CLASS.R, (O->CLASS.N+64)*sizeof(double));
              O->CLASS.LMX=realloc(O->CLASS.LMX, (O->CLASS.N+64)*sizeof(int));
              O->CLASS.IDXDN=realloc(O->CLASS.IDXDN, (O->CLASS.N+64)*sizeof(int*));
              O->CLASS.CONF=realloc(O->CLASS.CONF, (O->CLASS.N+64)*sizeof(int*));
              O->CLASS.IDMOD=realloc(O->CLASS.IDMOD, (O->CLASS.N+64)*sizeof(int*));
              for(i=O->CLASS.N;i<O->CLASS.N+64;i++)
              {
                O->CLASS.LMX[i]=-1;
                O->CLASS.ATOM[i]=malloc(16*sizeof(char));
                O->CLASS.IDXDN[i]=malloc(4*sizeof(int));
                O->CLASS.CONF[i]=malloc(4*sizeof(int));
                O->CLASS.IDMOD[i]=malloc(4*sizeof(int));
              }
            }
            StreamReadStuttgartCTRLchar(stream, O->CLASS.ATOM[O->CLASS.N-1]);
          }
          else if(O->CLASS.N==-1)
          {
            printf("In Category CLASS encoutered token %s before token ATOM\n",str);
          }
          else if(strcmp(str, "NEWNAM")==0)
          {
            StreamReadStuttgartCTRLBool(stream, 1, &(O->CLASS.NEWNAM[O->CLASS.N-1]));
          }
          else if(strcmp(str, "Z")==0)
          {
            StreamReadStuttgartCTRLint(stream, 1, &(O->CLASS.Z[O->CLASS.N-1]));
          }
          else if(strcmp(str, "R")==0)
          {
            StreamReadStuttgartCTRLdouble(stream, 1, &(O->CLASS.R[O->CLASS.N-1]));
          }
          else if(strcmp(str, "LMX")==0)
          {
            StreamReadStuttgartCTRLint(stream, 1, &(O->CLASS.LMX[O->CLASS.N-1]));
          }
          else if(strcmp(str, "CONF")==0)
          {
            StreamReadStuttgartCTRLint(stream, 4, O->CLASS.CONF[O->CLASS.N-1]);
          }
          else if(strcmp(str, "IDXDN")==0)
          {
            if(O->CLASS.LMX[O->SYML.N-1]==-1)
            {
              printf("In Category CLASS encoutered token IDXDN for atom %s (%i) before token LMX was found\n",O->CLASS.ATOM[O->CLASS.N-1],O->CLASS.N);
            }
            else
            {
              StreamReadStuttgartCTRLint(stream, O->CLASS.LMX[O->CLASS.N-1]+1, O->CLASS.IDXDN[O->CLASS.N-1]);
            }
          }
          else if(strcmp(str, "IDMOD")==0)
          {
            if(O->CLASS.LMX[O->SYML.N-1]==-1)
            {
              printf("In Category CLASS encoutered token IDMOD for atom %s (%i) before token LMX was found\n",O->CLASS.ATOM[O->CLASS.N-1],O->CLASS.N);
            }
            else
            {
              StreamReadStuttgartCTRLint(stream, O->CLASS.LMX[O->CLASS.N-1]+1, O->CLASS.IDMOD[O->CLASS.N-1]);
            }
          }
          else
          {
            printf("In Category %s undefined token %s\n",Categories[iCategory-1],str);
          }
        break;
        case 8: /* SITE */
          if(strcmp(str, "ATOM")==0)
          {
            if(O->SITE.N==-1)
            {
              O->SITE.ATOM=malloc(64*sizeof(char*));
              O->SITE.POS=malloc(64*sizeof(double*));
              O->SITE.X=malloc(64*sizeof(double*));
              for(i=0;i<64;i++)
              {
                O->SITE.ATOM[i]=malloc(16*sizeof(char));
                O->SITE.POS[i]=malloc(3*sizeof(double));
                O->SITE.X[i]=malloc(3*sizeof(double));
              }
              O->SITE.N=0;
            }
            O->SITE.N++;
            if((O->SITE.N%64)==0)
            {
              O->SITE.ATOM=realloc(O->SITE.ATOM,(O->SITE.N+64)*sizeof(char*));
              O->SITE.POS=realloc(O->SITE.POS,(O->SITE.N+64)*sizeof(double*));
              O->SITE.X=realloc(O->SITE.X,(O->SITE.N+64)*sizeof(double*));
              for(i=O->SITE.N;i<O->SITE.N+64;i++)
              {
                O->SITE.ATOM[i]=malloc(16*sizeof(char));
                O->SITE.POS[i]=malloc(3*sizeof(double));
                O->SITE.X[i]=malloc(3*sizeof(double));
              }
            }
            StreamReadStuttgartCTRLchar(stream, O->SITE.ATOM[O->SITE.N-1]);
          }
          else if(O->SITE.N==-1)
          {
            printf("In Category SITE encoutered token %s before token ATOM\n",str);
          }
          else if(strcmp(str, "POS")==0)
          {
            StreamReadStuttgartCTRLdouble(stream, 3, O->SITE.POS[O->SITE.N-1]);
          }
          else if(strcmp(str, "X")==0)
          {
            StreamReadStuttgartCTRLdouble(stream, 3, O->SITE.X[O->SITE.N-1]);
          }
          else
          {
            printf("In Category %s undefined token %s\n",Categories[iCategory-1],str);
          }
        break;
        case 9: /* SYMGRP */
          if(strcmp(str, "USESYM")==0)
          {
            StreamReadStuttgartCTRLBool(stream, 1, &(O->SYMGRP.USESYM));
          }
          else if(strcmp(str, "NGEN")==0)
          {
            StreamReadStuttgartCTRLint(stream, 1, &(O->SYMGRP.NGEN));
          }
          else if(strcmp(str, "IORIGIN")==0)
          {
            StreamReadStuttgartCTRLint(stream, 1, &(O->SYMGRP.IORIGIN));
          }
          else if(strcmp(str, "SPCGRP")==0)
          {
            StreamReadStuttgartCTRLchar(stream, O->SYMGRP.SPCGRP);
          }
          else if(strcmp(str, "GENGRP")==0)
          {
            if(O->SYMGRP.NGEN==-1)
            {
              printf("in category SYMGRP tried to read token GENGRP before token NGEN\n");
            }
            else
            {
              for(int i=0;i<O->SYMGRP.NGEN;i++)
              {
                StreamReadStuttgartCTRLchar(stream, O->SYMGRP.GENGRP[i]);
              }
            }
          }
          else
          {
            printf("In Category %s undefined token %s\n",Categories[iCategory-1],str);
          }
        break;
        case 10: /* SCALE */
          if(strcmp(str, "SCLWSR")==0)
          {
            StreamReadStuttgartCTRLBool(stream, 1, &(O->SCALE.SCLWSR));
          }
          else if(strcmp(str, "FACVOL")==0)
          {
            StreamReadStuttgartCTRLdouble(stream, 1, &(O->SCALE.FACVOL));
          }
          else if(strcmp(str, "OMMAX1")==0)
          {
            StreamReadStuttgartCTRLdouble(stream, 3, O->SCALE.OMMAX1);
          }
          else if(strcmp(str, "OMMAX2")==0)
          {
            StreamReadStuttgartCTRLdouble(stream, 3, O->SCALE.OMMAX2);
          }
          else if(strcmp(str, "GAMMA")==0)
          {
            StreamReadStuttgartCTRLdouble(stream, 1, &(O->SCALE.GAMMA));
          }
          else
          {
            printf("In Category %s undefined token %s\n",Categories[iCategory-1],str);
          }
        break;
        case 11: /* STR */
          if(strcmp(str, "DOWATS")==0)
          {
            StreamReadStuttgartCTRLBool(stream, 1, &(O->STR.DOWATS));
          }
          else if(strcmp(str, "DONALP")==0)
          {
            StreamReadStuttgartCTRLBool(stream, 1, &(O->STR.DONALP));
          }
          else if(strcmp(str, "NOCALC")==0)
          {
            StreamReadStuttgartCTRLBool(stream, 1, &(O->STR.NOCALC));
          }
          else if(strcmp(str, "IALPHA")==0)
          {
            StreamReadStuttgartCTRLint(stream, 1, &(O->STR.IALPHA));
          }
          else if(strcmp(str, "LMAXW")==0)
          {
            StreamReadStuttgartCTRLint(stream, 1, &(O->STR.LMAXW));
          }
          else if(strcmp(str, "NDIMIN")==0)
          {
            StreamReadStuttgartCTRLint(stream, 1, &(O->STR.NDIMIN));
          }
          else if(strcmp(str, "NTRANS")==0)
          {
            StreamReadStuttgartCTRLint(stream, 1, &(O->STR.NTRANS));
          }
          else if(strcmp(str, "LDN")==0)
          {
            StreamReadStuttgartCTRLint(stream, 1, &(O->STR.LDN));
          }
          else if(strcmp(str, "MDN")==0)
          {
            StreamReadStuttgartCTRLint(stream, 1, &(O->STR.MDN));
          }
          else if(strcmp(str, "JBASDN")==0)
          {
            StreamReadStuttgartCTRLint(stream, 1, &(O->STR.JBASDN));
          }
          else if(strcmp(str, "KAPPA2")==0)
          {
            StreamReadStuttgartCTRLdouble(stream, 1, &(O->STR.KAPPA2));
          }
          else if(strcmp(str, "KW**2")==0)
          {
            StreamReadStuttgartCTRLdouble(stream, 1, &(O->STR.KWxx2));
          }
          else if(strcmp(str, "ADOT")==0)
          {
            StreamReadStuttgartCTRLdouble(stream, 1, &(O->STR.ADOT));
          }
          else if(strcmp(str, "DELTR")==0)
          {
            StreamReadStuttgartCTRLdouble(stream, 1, &(O->STR.DELTR));
          }
          else if(strcmp(str, "RMAXS")==0)
          {
            StreamReadStuttgartCTRLdouble(stream, 1, &(O->STR.RMAXS));
          }
          else if(strcmp(str, "EPS")==0)
          {
            StreamReadStuttgartCTRLdouble(stream, 1, &(O->STR.EPS));
          }
          
          
          
          else
          {
            printf("In Category %s undefined token %s\n",Categories[iCategory-1],str);
          }
        break;
        case 12: /* BZ */
          if(strcmp(str, "TETRA")==0)
          {
            StreamReadStuttgartCTRLBool(stream, 1, &(O->BZ.TETRA));
          }
          else if(strcmp(str, "METAL")==0)
          {
            StreamReadStuttgartCTRLBool(stream, 1, &(O->BZ.METAL));
          }
          else if(strcmp(str, "NKABC")==0)
          {
            StreamReadStuttgartCTRLint(stream, 3, O->BZ.NKABC);
          }
          else if(strcmp(str, "MAXKP")==0)
          {
            StreamReadStuttgartCTRLint(stream, 1, &(O->BZ.MAXKP));
          }
          else if(strcmp(str, "N")==0)
          {
            StreamReadStuttgartCTRLint(stream, 1, &(O->BZ.N));
          }
          else if(strcmp(str, "RANGE")==0)
          {
            StreamReadStuttgartCTRLint(stream, 1, &(O->BZ.RANGE));
          }
          else if(strcmp(str, "NPTS")==0)
          {
            StreamReadStuttgartCTRLint(stream, 1, &(O->BZ.NPTS));
          }
          else if(strcmp(str, "W")==0)
          {
            StreamReadStuttgartCTRLdouble(stream, 1, &(O->BZ.W));
          }
          else if(strcmp(str, "TOL")==0)
          {
            StreamReadStuttgartCTRLdouble(stream, 1, &(O->BZ.TOL));
          }
          else
          {
            printf("In Category %s undefined token %s\n",Categories[iCategory-1],str);
          }
        break;
        case 13: /* EWALD */
          if(strcmp(str, "NKDMX")==0)
          {
            StreamReadStuttgartCTRLint(stream, 1, &(O->EWALD.NKDMX));
          }
          else if(strcmp(str, "AS")==0)
          {
            StreamReadStuttgartCTRLdouble(stream, 1, &(O->EWALD.AS));
          }
          else if(strcmp(str, "TOL")==0)
          {
            StreamReadStuttgartCTRLdouble(stream, 1, &(O->EWALD.TOL));
          }
          else
          {
            printf("In Category %s undefined token %s\n",Categories[iCategory-1],str);
          }
        break;
        case 14: /* DOS */
          if(strcmp(str, "NOPTS")==0)
          {
            StreamReadStuttgartCTRLint(stream, 1, &(O->DOS.NOPTS));
          }
          else if(strcmp(str, "EMIN")==0)
          {
            StreamReadStuttgartCTRLdouble(stream, 1, &(O->DOS.EMIN));
          }
          else if(strcmp(str, "EMAX")==0)
          {
            StreamReadStuttgartCTRLdouble(stream, 1, &(O->DOS.EMAX));
          }
          else
          {
            printf("In Category %s undefined token %s\n",Categories[iCategory-1],str);
          }
        break;
        case 15: /* SYML */
          if(strcmp(str, "NQ")==0)
          {
            if(O->SYML.N==-1)
            {
              O->SYML.NQ=malloc(64*sizeof(double));
              O->SYML.Q1=malloc(64*sizeof(double*));
              O->SYML.Q2=malloc(64*sizeof(double*));
              O->SYML.LAB1=malloc(64*sizeof(char*));
              O->SYML.LAB2=malloc(64*sizeof(char*));
              for(i=0;i<64;i++)
              {
                O->SYML.Q1[i]=malloc(3*sizeof(double));
                O->SYML.Q2[i]=malloc(3*sizeof(double));
                O->SYML.LAB1[i]=malloc(32*sizeof(char));
                O->SYML.LAB2[i]=malloc(32*sizeof(char));
              }
              O->SYML.N=0;
            }
            O->SYML.N++;
            if((O->SYML.N%64)==0)
            {
              O->SYML.NQ=realloc(O->SYML.NQ, (O->SYML.N+64)*sizeof(double));
              O->SYML.Q1=realloc(O->SYML.Q1, (O->SYML.N+64)*sizeof(double*));
              O->SYML.Q2=realloc(O->SYML.Q2, (O->SYML.N+64)*sizeof(double*));
              O->SYML.LAB1=realloc(O->SYML.LAB1, (O->SYML.N+64)*sizeof(char*));
              O->SYML.LAB2=realloc(O->SYML.LAB2, (O->SYML.N+64)*sizeof(char*));
              for(i=O->SYML.N;i<O->SYML.N+64;i++)
              {
                O->SYML.Q1[i]=malloc(3*sizeof(double));
                O->SYML.Q2[i]=malloc(3*sizeof(double));
                O->SYML.LAB1[i]=malloc(32*sizeof(char));
                O->SYML.LAB2[i]=malloc(32*sizeof(char));
              }
            }
            StreamReadStuttgartCTRLint(stream, 1, &(O->SYML.NQ[O->SYML.N-1]));
          }
          else if(strcmp(str, "Q1")==0)
          {
            StreamReadStuttgartCTRLdouble(stream, 3, O->SYML.Q1[O->SYML.N-1]);
          }
          else if(strcmp(str, "Q2")==0)
          {
            StreamReadStuttgartCTRLdouble(stream, 3, O->SYML.Q2[O->SYML.N-1]);
          }
          else if(strcmp(str, "LAB1")==0)
          {
            StreamReadStuttgartCTRLchar(stream, O->SYML.LAB1[O->SYML.N-1]);
          }
          else if(strcmp(str, "LAB2")==0)
          {
            StreamReadStuttgartCTRLchar(stream, O->SYML.LAB2[O->SYML.N-1]);
          }
          else
          {
            printf("In Category %s undefined token %s\n",Categories[iCategory-1],str);
          }
        break;
        case 16: /* START */
          if(strcmp(str, "NIT")==0)
          {
            StreamReadStuttgartCTRLint(stream, 1, &(O->START.NIT));
          }
          
          
          
          
          
          
          
          
          
          
          
          else
          {
            printf("In Category %s undefined token %s\n",Categories[iCategory-1],str);
          }
        break;
        case 17: /* HARTREE */
          if(strcmp(str, "BEGATOM")==0)
          {
            StreamReadStuttgartCTRLBool(stream, 1, &(O->HARTREE.BEGATOM));
          }
          else if(strcmp(str, "LT1")==0)
          {
            StreamReadStuttgartCTRLint(stream, 1, &(O->HARTREE.LT1));
          }
          else if(strcmp(str, "LT2")==0)
          {
            StreamReadStuttgartCTRLint(stream, 1, &(O->HARTREE.LT2));
          }
          else if(strcmp(str, "LT3")==0)
          {
            StreamReadStuttgartCTRLint(stream, 1, &(O->HARTREE.LT3));
          }
          else
          {
            printf("In Category %s undefined token %s\n",Categories[iCategory-1],str);
          }
        break;
        case 18: /* PLOT */
          if(strcmp(str, "NDELR1")==0)
          {
            StreamReadStuttgartCTRLint(stream, 1, &(O->PLOT.NDELR1));
          }
          else if(strcmp(str, "NDELR2")==0)
          {
            StreamReadStuttgartCTRLint(stream, 1, &(O->PLOT.NDELR2));
          }
          else if(strcmp(str, "NDELR3")==0)
          {
            StreamReadStuttgartCTRLint(stream, 1, &(O->PLOT.NDELR3));
          }
          else if(strcmp(str, "FORMAT")==0)
          {
            StreamReadStuttgartCTRLint(stream, 1, &(O->PLOT.FORMAT));
          }
          else if(strcmp(str, "ORIGIN")==0)
          {
            StreamReadStuttgartCTRLdouble(stream, 3, O->PLOT.ORIGIN);
          }
          else if(strcmp(str, "R1")==0)
          {
            StreamReadStuttgartCTRLdouble(stream, 3, O->PLOT.R1);
          }
          else if(strcmp(str, "R2")==0)
          {
            StreamReadStuttgartCTRLdouble(stream, 3, O->PLOT.R2);
          }
          else if(strcmp(str, "R3")==0)
          {
            StreamReadStuttgartCTRLdouble(stream, 3, O->PLOT.R3);
          }
          else if(strcmp(str, "RMIN")==0)
          {
            StreamReadStuttgartCTRLdouble(stream, 1, &(O->PLOT.RMIN));
          }
          else if(strcmp(str, "RMAX")==0)
          {
            StreamReadStuttgartCTRLdouble(stream, 1, &(O->PLOT.RMAX));
          }
          else if(strcmp(str, "TMIN")==0)
          {
            StreamReadStuttgartCTRLdouble(stream, 1, &(O->PLOT.TMIN));
          }
          else if(strcmp(str, "TMAX")==0)
          {
            StreamReadStuttgartCTRLdouble(stream, 1, &(O->PLOT.TMAX));
          }
          else
          {
            printf("In Category %s undefined token %s\n",Categories[iCategory-1],str);
          }
        break;
        case 19: /* CHARGE */
          if(strcmp(str, "LMTODAT")==0)
          {
            StreamReadStuttgartCTRLBool(stream, 1, &(O->CHARGE.LMTODAT));
          }
          else if(strcmp(str, "CHARWIN")==0)
          {
            StreamReadStuttgartCTRLBool(stream, 1, &(O->CHARGE.LMTODAT));
          }
          else if(strcmp(str, "ELF")==0)
          {
            StreamReadStuttgartCTRLBool(stream, 1, &(O->CHARGE.LMTODAT));
          }
          else if(strcmp(str, "ADDCOR")==0)
          {
            StreamReadStuttgartCTRLBool(stream, 1, &(O->CHARGE.LMTODAT));
          }
          else if(strcmp(str, "SPINDENS")==0)
          {
            StreamReadStuttgartCTRLBool(stream, 1, &(O->CHARGE.LMTODAT));
          }
          else if(strcmp(str, "EMIN")==0)
          {
            StreamReadStuttgartCTRLdouble(stream, 1, &(O->CHARGE.EMIN));
          }
          else if(strcmp(str, "EMAX")==0)
          {
            StreamReadStuttgartCTRLdouble(stream, 1, &(O->CHARGE.EMAX));
          }
          else
          {
            printf("In Category %s undefined token %s\n",Categories[iCategory-1],str);
          }
        break;
        case 20: /* FINDES */
          if(strcmp(str, "MAXPT")==0)
          {
            StreamReadStuttgartCTRLint(stream, 1, &(O->FINDES.MAXPT));
          }
          else if(strcmp(str, "NRXYZ")==0)
          {
            StreamReadStuttgartCTRLint(stream, 3, O->FINDES.NRXYZ);
          }
          else if(strcmp(str, "RMINES")==0)
          {
            StreamReadStuttgartCTRLdouble(stream, 1, &(O->FINDES.RMINES));
          }
          else if(strcmp(str, "RMAXES")==0)
          {
            StreamReadStuttgartCTRLdouble(stream, 1, &(O->FINDES.RMAXES));
          }
          else
          {
            printf("In Category %s undefined token %s\n",Categories[iCategory-1],str);
          }
        break;
        case 21: /* SCELL */
          if(strcmp(str, "EQUIV")==0)
          {
            StreamReadStuttgartCTRLBool(stream, 1, &(O->SCELL.EQUIV));
          }
          else if(strcmp(str, "PLAT")==0)
          {
            StreamReadStuttgartCTRLdouble(stream, 9, &(O->SCELL.PLAT[0][0]));
          }
          else
          {
            printf("In Category %s undefined token %s\n",Categories[iCategory-1],str);
          }
        break;
        case 22: /* RHOFIT */
          if(strcmp(str, "FIT")==0)
          {
            StreamReadStuttgartCTRLBool(stream, 1, &(O->RHOFIT.FIT));
          }
          else if(strcmp(str, "KAPPA2")==0)
          {
            StreamReadStuttgartCTRLdouble(stream, 1, &(O->RHOFIT.KAPPA2));
          }
          else if(strcmp(str, "KW**2")==0)
          {
            StreamReadStuttgartCTRLdouble(stream, 1, &(O->RHOFIT.KWxx2));
          }
          else if(strcmp(str, "RMAXS")==0)
          {
            StreamReadStuttgartCTRLdouble(stream, 1, &(O->RHOFIT.RMAXS));
          }
          
          
          
          
          
          
          
          
          else
          {
            printf("In Category %s undefined token %s\n",Categories[iCategory-1],str);
          }
        break;
        case 23: /* NMTO */
          if(strcmp(str, "NENU")==0)
          {
            StreamReadStuttgartCTRLint(stream, 1, &(O->NMTO.NENU));
            O->NMTO.ENU=malloc(O->NMTO.NENU*sizeof(double));
          }
          else if(strcmp(str, "ENU")==0)
          {
            if(O->NMTO.NENU==-1)
            {
              printf("Encountered token ENU in NMTO but have not set NENU\n");
            }
            else
            {
              StreamReadStuttgartCTRLdouble(stream, O->NMTO.NENU, O->NMTO.ENU);
            }
          }
          else
          {
            printf("In Category %s undefined token %s\n",Categories[iCategory-1],str);
          }
        break;
        case 24: /* LOCALCS */
          if(strcmp(str, "Number_of_rotated_atoms")==0)
          {
            StreamReadStuttgartCTRLint(stream, 1, &(O->LOCALCS.Number_of_rotated_atoms));
            O->LOCALCS.LOCAL_Z=malloc(O->LOCALCS.Number_of_rotated_atoms*sizeof(double*));
            O->LOCALCS.LOCAL_X=malloc(O->LOCALCS.Number_of_rotated_atoms*sizeof(double*));
            for(i=0;i<O->LOCALCS.Number_of_rotated_atoms;i++)
            {
              O->LOCALCS.LOCAL_Z[i]=malloc(3*sizeof(double));
              O->LOCALCS.LOCAL_X[i]=malloc(3*sizeof(double));
            }
            O->LOCALCS.ATOM=malloc(O->LOCALCS.Number_of_rotated_atoms*sizeof(int));
            O->LOCALCS.L=malloc(O->LOCALCS.Number_of_rotated_atoms*sizeof(int));
          }
          else if(strcmp(str, "Data_as")==0)
          {
            StreamReadStuttgartCTRLchar(stream, O->LOCALCS.Data_as);
          }
          else if(strcmp(str, "ATOM")==0)
          {
            Nlocalcs++;
            if(Nlocalcs<O->LOCALCS.Number_of_rotated_atoms)
            {
              StreamReadStuttgartCTRLint(stream, 1, &(O->LOCALCS.ATOM[Nlocalcs]));
            }
            else
            {
              printf("Encountered token ATOM in LOCALCS for the %d time\n",Nlocalcs+1);
              printf("But Number_of_rotated_atoms is %d, so token is neglected\n",O->LOCALCS.Number_of_rotated_atoms);
            }
          }
          else if(strcmp(str, "L")==0)
          {
            if(Nlocalcs<O->LOCALCS.Number_of_rotated_atoms)
            {
              StreamReadStuttgartCTRLint(stream, 1, &(O->LOCALCS.L[Nlocalcs]));
            }
          }
          else if(strcmp(str, "LOCAL_Z")==0)
          {
            if(Nlocalcs<O->LOCALCS.Number_of_rotated_atoms)
            {
              StreamReadStuttgartCTRLdouble(stream, 3, O->LOCALCS.LOCAL_Z[Nlocalcs]);
            }
          }
          else if(strcmp(str, "LOCAL_X")==0)
          {
            if(Nlocalcs<O->LOCALCS.Number_of_rotated_atoms)
            {
              StreamReadStuttgartCTRLdouble(stream, 3, O->LOCALCS.LOCAL_X[Nlocalcs]);
            }
          }
          else
          {
            printf("In Category %s undefined token %s\n",Categories[iCategory-1],str);
          }
        break;
        default:
          printf("In undefined Category %s token %s\n",Categories[iCategory-1],str);
        break;
      }
    }
  }
  return 0;
}
uint32_t ReadStuttgartCTRL(StuttgartCTRL * O)
{
  return StreamReadStuttgartCTRL(O, stdin);
}
uint32_t FileReadStuttgartCTRL(StuttgartCTRL * O, const char *const flname)
{
  FILE  *stream;
  uint32_t retval;
  if ((stream = fopen(flname, "r")) == NULL)
  {
    printf("ERROR\nCan't open %s for reading.\n", flname);
    return 1;
  }
  retval = StreamReadStuttgartCTRL(O, stream);
  fclose(stream);
  return retval;
}
