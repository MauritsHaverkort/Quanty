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
  O->STR.N=-1;
  O->START.N=-1;
  O->RHOFIT.N=-1;
  O->SYML.N=-1;
  O->LOCALCS.Number_of_rotated_atoms=-1;
  O->NMTO.NENU=-1;
  O->HEADER[0]='\0';
  O->VERS[0]='\0';
  return 0;
}

/*
 *  Check consistency of CTRL structure
 *  order atoms in STR, START and RHOFIT the same as in CLASS
 *  set the values of ICLASS
 *  set the values of IDXDN_M for all elements of SITE according to the definition in CLASS
 *    if SITE IDXDN_M = 0 it will be set to CLASS IDXDN
 *  set SITE.Nactive and Site.ActiveAtomMap
 */
uint32_t OrderAndCheckStuttgartCTRL(StuttgartCTRL * O)
{
  uint32_t i,j,l,m;
  /* check consistency of number of ATOMS */
  if(O->SITE.N<O->CLASS.N)
  {
    printf("Erorr in CTRL file, number of sites smaller then number of classes\n");
    return 1;
  }
  if(!((O->START.N==-1)||(O->START.N==O->CLASS.N)))
  {
    printf("Erorr in CTRL file, number of ATOMS in START (%i) different from number of ATOMS in CLASS (%i)\n",O->START.N,O->CLASS.N);
    return 1;
  }
  if(!((O->STR.N==-1)||(O->STR.N==O->CLASS.N)))
  {
    printf("Erorr in CTRL file, number of ATOMS in STR (%i) different from number of ATOMS in CLASS (%i)\n",O->STR.N,O->CLASS.N);
    return 1;
  }
  if(!((O->RHOFIT.N==-1)||(O->RHOFIT.N==O->CLASS.N)))
  {
    printf("Erorr in CTRL file, number of ATOMS in RHOFIT (%i) different from number of ATOMS in CLASS (%i)\n",O->RHOFIT.N,O->CLASS.N);
    return 1;
  }
  /* place atoms in order of CLASS in STR, START and RHOFIT */
  
  /* TO DO */
  
  /* set values of ICLASS */
  O->SITE.ICLASS=malloc(O->SITE.N*sizeof(int));
  for(i=0;i<O->SITE.N;i++)
  {
    O->SITE.ICLASS[i]=-1;
    for(j=0;j<O->CLASS.N;j++)
    {
      if(strcmp(O->SITE.ATOM[i], O->CLASS.ATOM[j])==0)
      {
        O->SITE.ICLASS[i]=j;
        break;
      }
    }
    if(O->SITE.ICLASS[i]==-1)
    {
      printf("ERROR in CTRL file ATOM %s has no CLASS\n",O->SITE.ATOM[i]);
    }
  }
  /* set values to IDXDN_M */
  for(i=0;i<O->SITE.N;i++)
  {
    for(l=0;l<O->CLASS.LMX[O->SITE.ICLASS[i]];l++)
    {
      for(m=l*l;m<l*l+2*l+1;m++)
      {
        if(O->SITE.IDXDN_M[i][m]==0)O->SITE.IDXDN_M[i][m]=O->CLASS.IDXDN[O->SITE.ICLASS[i]][l];
      }
    }
  }
  /* set SITE.Nactive and SITE.ActiveAtomMap */
  O->SITE.ActiveAtomMap=malloc(O->SITE.N*sizeof(int));
  O->SITE.Nactive=0;
  for(i=0;i<O->SITE.N;i++)
  {
    O->SITE.ActiveAtomMap[i]=-1;
    for(l=0;l<O->CLASS.LMX[O->SITE.ICLASS[i]];l++)
    {
      for(m=l*l;m<l*l+2*l+1;m++)
      {
        if(O->SITE.IDXDN_M[i][m]==1)
        {
          l++;
          O->SITE.ActiveAtomMap[i]=O->SITE.Nactive;
          O->SITE.Nactive++;
          break;
        }
      }
    }
  }
  
  
  
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
 *   Reads A CTRL file of the stuttgart NMTO/LMTO code into the StuttgartCTRL struct
 */
uint32_t StreamReadStuttgartCTRL(StuttgartCTRL * O, FILE * stream)
{
  const char * const Categories[] = StuttgartCTRLCategories;
  char str[128];
  int Continue=1, ContinueCategory, iCategory, i,LMXnow;
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
      if(fscanfCTRL(stream, str, &ForSureAToken)!=1)
      {
        ContinueCategory=0;
        Continue=0;
        break;
      }
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
                O->CLASS.IDXDN[i]=calloc(4,sizeof(int));
                O->CLASS.CONF[i]=calloc(4,sizeof(int));
                O->CLASS.IDMOD[i]=calloc(4,sizeof(int));
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
                O->CLASS.IDXDN[i]=calloc(4,sizeof(int));
                O->CLASS.CONF[i]=calloc(4,sizeof(int));
                O->CLASS.IDMOD[i]=calloc(4,sizeof(int));
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
            if(O->CLASS.LMX[O->CLASS.N-1]==-1)
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
            if(O->CLASS.LMX[O->CLASS.N-1]==-1)
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
              O->SITE.IDXDN_M=malloc(64*sizeof(int*));
              for(i=0;i<64;i++)
              {
                O->SITE.ATOM[i]=malloc(16*sizeof(char));
                O->SITE.POS[i]=malloc(3*sizeof(double));
                O->SITE.X[i]=malloc(3*sizeof(double));
                O->SITE.IDXDN_M[i]=calloc(16,sizeof(double));
              }
              O->SITE.N=0;
            }
            O->SITE.N++;
            if((O->SITE.N%64)==0)
            {
              O->SITE.ATOM=realloc(O->SITE.ATOM,(O->SITE.N+64)*sizeof(char*));
              O->SITE.POS=realloc(O->SITE.POS,(O->SITE.N+64)*sizeof(double*));
              O->SITE.X=realloc(O->SITE.X,(O->SITE.N+64)*sizeof(double*));
              O->SITE.IDXDN_M=realloc(O->SITE.IDXDN_M,(O->SITE.N+64)*sizeof(int*));
              for(i=O->SITE.N;i<O->SITE.N+64;i++)
              {
                O->SITE.ATOM[i]=malloc(16*sizeof(char));
                O->SITE.POS[i]=malloc(3*sizeof(double));
                O->SITE.X[i]=malloc(3*sizeof(double));
                O->SITE.IDXDN_M[i]=calloc(16,sizeof(double));
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
          else if(strcmp(str, "IDXDN_M")==0)
          {
            LMXnow=-1;
            for(i=0;i<O->CLASS.N;i++)
            {
              if(strcmp(O->CLASS.ATOM[i],O->SITE.ATOM[O->SITE.N-1])==0)
              {
                LMXnow=O->CLASS.LMX[i];
              }
            }
            switch (LMXnow)
            {
              case -1:
                printf("Found token IDXDN_M in Category SITE for ATOM %s, but have not found the token LMX in the CLASS of this ATOM\n",O->SITE.ATOM[O->SITE.N-1]);
                break;
              case  0:
                StreamReadStuttgartCTRLint(stream, 1, O->SITE.IDXDN_M[O->SITE.N-1]);
                break;
              case  1:
                StreamReadStuttgartCTRLint(stream, 4, O->SITE.IDXDN_M[O->SITE.N-1]);
                break;
              case  2:
                StreamReadStuttgartCTRLint(stream, 9, O->SITE.IDXDN_M[O->SITE.N-1]);
                break;
              case  3:
                StreamReadStuttgartCTRLint(stream, 16, O->SITE.IDXDN_M[O->SITE.N-1]);
                break;
              default:
                printf("No clue how I got here LMXnow=%i\n",LMXnow);
                break;
            }
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
          else if(strcmp(str, "ATOM")==0)
          {
            if(O->STR.N==-1)
            {
              O->STR.ATOM=malloc(64*sizeof(char*));
              O->STR.SIGMA=malloc(64*sizeof(double*));
              for(i=0;i<64;i++)
              {
                O->STR.ATOM[i]=malloc(16*sizeof(char));
                O->STR.SIGMA[i]=malloc(4*sizeof(double));
              }
              O->STR.N=0;
            }
            O->STR.N++;
            if((O->STR.N%64)==0)
            {
              O->STR.ATOM=realloc(O->STR.ATOM, (O->STR.N+64)*sizeof(char*));
              O->STR.SIGMA=realloc(O->STR.SIGMA, (O->STR.N+64)*sizeof(double*));
              for(i=O->STR.N;i<O->STR.N+64;i++)
              {
                O->STR.ATOM[i]=malloc(16*sizeof(char));
                O->STR.SIGMA[i]=malloc(4*sizeof(double));
              }
            }
            StreamReadStuttgartCTRLchar(stream, O->STR.ATOM[O->STR.N-1]);
          }
          else if(O->STR.N==-1)
          {
            printf("In Category STR encoutered token %s before token ATOM\n",str);
          }
          else if(strcmp(str, "SIGMA")==0)
          {
            LMXnow=-1;
            for(i=0;i<O->CLASS.N;i++)
            {
              if(strcmp(O->CLASS.ATOM[i],O->STR.ATOM[O->STR.N-1])==0)
              {
                LMXnow=O->CLASS.LMX[i];
              }
            }
            if(LMXnow==-1)
            {
              printf("Found token SIGMA in Category STR for ATOM %s, but have not found the token LMX in the CLASS of this ATOM\n",O->STR.ATOM[O->STR.N-1]);
            }
            else
            {
              StreamReadStuttgartCTRLdouble(stream, LMXnow+1, O->STR.SIGMA[O->STR.N-1]);
            }
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
          if(strcmp(str, "BROY")==0)
          {
            StreamReadStuttgartCTRLBool(stream, 1, &(O->START.BROY));
          }
          else if(strcmp(str, "FREE")==0)
          {
            StreamReadStuttgartCTRLBool(stream, 1, &(O->START.FREE));
          }
          else if(strcmp(str, "BEGMOM")==0)
          {
            StreamReadStuttgartCTRLBool(stream, 1, &(O->START.BEGMOM));
          }
          else if(strcmp(str, "CNTROL")==0)
          {
            StreamReadStuttgartCTRLBool(stream, 1, &(O->START.CNTROL));
          }
          else if(strcmp(str, "NIT")==0)
          {
            StreamReadStuttgartCTRLint(stream, 1, &(O->START.NIT));
          }
          else if(strcmp(str, "NMIX")==0)
          {
            StreamReadStuttgartCTRLint(stream, 1, &(O->START.NMIX));
          }
          else if(strcmp(str, "WC")==0)
          {
            StreamReadStuttgartCTRLdouble(stream, 1, &(O->START.WC));
          }
          else if(strcmp(str, "BETA")==0)
          {
            StreamReadStuttgartCTRLdouble(stream, 1, &(O->START.BETA));
          }
          else if(strcmp(str, "CNVG")==0)
          {
            StreamReadStuttgartCTRLdouble(stream, 1, &(O->START.CNVG));
          }
          else if(strcmp(str, "CNVGET")==0)
          {
            StreamReadStuttgartCTRLdouble(stream, 1, &(O->START.CNVGET));
          }
          else if(strcmp(str, "EFERMI")==0)
          {
            StreamReadStuttgartCTRLdouble(stream, 1, &(O->START.EFERMI));
          }
          else if(strcmp(str, "VMTZ")==0)
          {
            StreamReadStuttgartCTRLdouble(stream, 1, &(O->START.VMTZ));
          }
          
          
          else if(strcmp(str, "ATOM")==0)
          {
            if(O->START.N==-1)
            {
              O->START.ATOM=malloc(64*sizeof(char*));
              O->START.P=malloc(64*sizeof(double*));
              O->START.Q=malloc(64*sizeof(double*));
              O->START.enu=malloc(64*sizeof(double*));
              O->START.c=malloc(64*sizeof(double*));
              O->START.sqrdel=malloc(64*sizeof(double*));
              O->START.p=malloc(64*sizeof(double*));
              O->START.gamma=malloc(64*sizeof(double*));
              for(i=0;i<64;i++)
              {
                O->START.ATOM[i]=malloc(16*sizeof(char));
                O->START.P[i]=malloc(4*sizeof(double));
                O->START.Q[i]=malloc(12*sizeof(double));
                O->START.enu[i]=malloc(4*sizeof(double));
                O->START.c[i]=malloc(4*sizeof(double));
                O->START.sqrdel[i]=malloc(4*sizeof(double));
                O->START.p[i]=malloc(4*sizeof(double));
                O->START.gamma[i]=malloc(4*sizeof(double));
              }
              O->START.N=0;
            }
            O->START.N++;
            if((O->START.N%64)==0)
            {
              O->START.ATOM=realloc(O->START.ATOM, (O->START.N+64)*sizeof(char*));
              O->START.P=realloc(O->START.P, (O->START.N+64)*sizeof(double*));
              O->START.Q=realloc(O->START.Q, (O->START.N+64)*sizeof(double*));
              O->START.enu=realloc(O->START.enu, (O->START.N+64)*sizeof(double*));
              O->START.c=realloc(O->START.c, (O->START.N+64)*sizeof(double*));
              O->START.sqrdel=realloc(O->START.sqrdel, (O->START.N+64)*sizeof(double*));
              O->START.p=realloc(O->START.p, (O->START.N+64)*sizeof(double*));
              O->START.gamma=realloc(O->START.gamma, (O->START.N+64)*sizeof(double*));
              for(i=O->START.N;i<O->START.N+64;i++)
              {
                O->START.ATOM[i]=malloc(16*sizeof(char));
                O->START.P[i]=malloc(4*sizeof(double));
                O->START.Q[i]=malloc(12*sizeof(double));
                O->START.enu[i]=malloc(4*sizeof(double));
                O->START.c[i]=malloc(4*sizeof(double));
                O->START.sqrdel[i]=malloc(4*sizeof(double));
                O->START.p[i]=malloc(4*sizeof(double));
                O->START.gamma[i]=malloc(4*sizeof(double));
              }
            }
            StreamReadStuttgartCTRLchar(stream, O->START.ATOM[O->START.N-1]);
          }
          else if(O->START.N==-1)
          {
            printf("In Category START encoutered token %s before token ATOM\n",str);
          }
          else
          {
            LMXnow=-1;
            for(i=0;i<O->CLASS.N;i++)
            {
              if(strcmp(O->CLASS.ATOM[i],O->START.ATOM[O->START.N-1])==0)
              {
                LMXnow=O->CLASS.LMX[i];
              }
            }
            if(LMXnow==-1)
            {
              printf("In Category %s undefined token %s or\n",Categories[iCategory-1],str);
              printf("Found token %s in Category START for ATOM %s, but have not found the token LMX in the CLASS of this ATOM\n",str,O->START.ATOM[O->START.N-1]);
            }
            else
            {
              if(strcmp(str, "P")==0)
              {
                StreamReadStuttgartCTRLdouble(stream, LMXnow+1, O->START.P[O->START.N-1]);
              }
              else if(strcmp(str, "Q")==0)
              {
                StreamReadStuttgartCTRLdouble(stream, 3*(LMXnow+1), O->START.Q[O->START.N-1]);
              }
              else if(strcmp(str, "enu")==0)
              {
                StreamReadStuttgartCTRLdouble(stream, LMXnow+1, O->START.enu[O->START.N-1]);
              }
              else if(strcmp(str, "c")==0)
              {
                StreamReadStuttgartCTRLdouble(stream, LMXnow+1, O->START.c[O->START.N-1]);
              }
              else if(strcmp(str, "sqrdel")==0)
              {
                StreamReadStuttgartCTRLdouble(stream, LMXnow+1, O->START.sqrdel[O->START.N-1]);
              }
              else if(strcmp(str, "p")==0)
              {
                StreamReadStuttgartCTRLdouble(stream, LMXnow+1, O->START.p[O->START.N-1]);
              }
              else if(strcmp(str, "gamma")==0)
              {
                StreamReadStuttgartCTRLdouble(stream, LMXnow+1, O->START.gamma[O->START.N-1]);
              }
            }
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
          else if(strcmp(str, "OUPTUP")==0)
          {
            StreamReadStuttgartCTRLint(stream, 1, &(O->RHOFIT.OUTPUT));
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
          else if(strcmp(str, "ATOM")==0)
          {
            if(O->RHOFIT.N==-1)
            {
              O->RHOFIT.ATOM=malloc(64*sizeof(char*));
              O->RHOFIT.LMXRHO=malloc(64*sizeof(int));
              O->RHOFIT.SIGMA=malloc(64*sizeof(double*));
              for(i=0;i<64;i++)
              {
                O->RHOFIT.LMXRHO[i]=-1;
                O->RHOFIT.ATOM[i]=malloc(16*sizeof(char));
                O->RHOFIT.SIGMA[i]=malloc(4*sizeof(double));
              }
              O->RHOFIT.N=0;
            }
            O->RHOFIT.N++;
            if((O->RHOFIT.N%64)==0)
            {
              O->RHOFIT.ATOM=realloc(O->RHOFIT.ATOM, (O->RHOFIT.N+64)*sizeof(char*));
              O->RHOFIT.LMXRHO=realloc(O->RHOFIT.LMXRHO, (O->RHOFIT.N+64)*sizeof(int));
              O->RHOFIT.SIGMA=realloc(O->RHOFIT.SIGMA, (O->RHOFIT.N+64)*sizeof(double*));
              for(i=O->RHOFIT.N;i<O->RHOFIT.N+64;i++)
              {
                O->RHOFIT.LMXRHO[i]=-1;
                O->RHOFIT.ATOM[i]=malloc(16*sizeof(char));
                O->RHOFIT.SIGMA[i]=malloc(4*sizeof(double));
              }
            }
            StreamReadStuttgartCTRLchar(stream, O->RHOFIT.ATOM[O->RHOFIT.N-1]);
          }
          else if(O->RHOFIT.N==-1)
          {
            printf("In Category RHOFIT encoutered token %s before token ATOM\n",str);
          }
          else if(strcmp(str, "LMXRHO")==0)
          {
            StreamReadStuttgartCTRLint(stream, 1, &(O->RHOFIT.LMXRHO[O->RHOFIT.N-1]));
          }
          else if(strcmp(str, "SIGMA")==0)
          {
            if(O->RHOFIT.LMXRHO[O->RHOFIT.N-1]==-1)
            {
              printf("In Category RHOFIT encoutered token SIGMA for atom %s (%i) before token LMXRHO was found\n",O->RHOFIT.ATOM[O->RHOFIT.N-1],O->CLASS.N);
            }
            else
            {
              StreamReadStuttgartCTRLdouble(stream, O->RHOFIT.LMXRHO[O->RHOFIT.N-1]+1, O->RHOFIT.SIGMA[O->RHOFIT.N-1]);
            }
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
          else if(strcmp(str, "PLOTATO")==0)
          {
            StreamReadStuttgartCTRLint(stream, 1, &(O->NMTO.PLOTATO));
          }
          else if(strcmp(str, "ORBITAL")==0)
          {
            StreamReadStuttgartCTRLint(stream, 1, &(O->NMTO.ORBITAL));
          }
          else if(strcmp(str, "RMAXHR")==0)
          {
            StreamReadStuttgartCTRLdouble(stream, 1, &(O->NMTO.RMAXHR));
          }
          else if(strcmp(str, "eV")==0)
          {
            StreamReadStuttgartCTRLBool(stream, 1, &(O->NMTO.eV));
          }
          else if(strcmp(str, "ORTHO")==0)
          {
            StreamReadStuttgartCTRLBool(stream, 1, &(O->NMTO.ORTHO));
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
  return OrderAndCheckStuttgartCTRL(O);
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


/*
 *  Reads a Stuttgart NMTO HAMR file into a StuttgartHAMR struct
 */
uint32_t StreamReadStuttgartHAMR(StuttgartHAMR * O, FILE * stream)
{
  char ReadLine[1024];
  char buffer[128];
  int ContinueLine=1;
  int ContinueLabel;
  ReadLine[1023]='\0';
  O->Nelement=-1;
  if(fgets(ReadLine, 1023, stream)==NULL)
  {
    printf("In function StreamReadStuttgartHAMR intial reading of file failed\n");
    return 1;
  }
  O->eV=False;
  while (ContinueLine)
  {
    if(!(strstr(ReadLine, "eV")==NULL))O->eV=True;
    if(strncmp(ReadLine,"INT",3)==0)
    {
      /* start of new element in HAMR make sure enough space allocated */
      if(O->Nelement==-1)
      {
        O->element=malloc(64*sizeof(StuttgartHAMRelement));
        O->Nelement=0;
      }
      O->Nelement++;
      if((O->Nelement%64)==0)
      {
        O->element=realloc(O->element, (O->Nelement+64)*sizeof(StuttgartHAMRelement));
      }
      if(fgets(ReadLine, 1023, stream)==NULL){printf("ERROR 1 in reading HAMR file\n");return 1;}
      if(sscanf(ReadLine, "%s %s %i %s %s %i",buffer,O->element[O->Nelement-1].Atom1,&(O->element[O->Nelement-1].IndexAtom1),buffer,O->element[O->Nelement-1].Atom2,&(O->element[O->Nelement-1].IndexAtom2))!=6)
      {
        printf("ERROR in reading Atom names and index from HAMR file in Line:\n");
        printf("%s\n",ReadLine);
        return 1;
      }
      if(fgets(ReadLine, 1023, stream)==NULL){printf("ERROR 2 in reading HAMR file\n");return 1;}
      if(sscanf(ReadLine, "%s %lE %s %s %s %i %i %i",buffer,&(O->element[O->Nelement-1].Dist),buffer,buffer,buffer,&(O->element[O->Nelement-1].TransVecs[0]),&(O->element[O->Nelement-1].TransVecs[1]),&(O->element[O->Nelement-1].TransVecs[2]))!=8)
      {
        printf("ERROR in reading Dist and TransVecs from HAMR file in Line:\n");
        printf("%s\n",ReadLine);
        return 1;
      }
      if(fgets(ReadLine, 1023, stream)==NULL){printf("ERROR 3 in reading HAMR file\n");return 1;}
      if(sscanf(ReadLine, "%s %s %lE %lE %lE",buffer,buffer, &(O->element[O->Nelement-1].ConnectingVector[0]),&(O->element[O->Nelement-1].ConnectingVector[1]),&(O->element[O->Nelement-1].ConnectingVector[2]))!=5)
      {
        printf("ERROR in reading ConnectingVector from HAMR file in Line:\n");
        printf("%s\n",ReadLine);
        return 1;
      }
      if(fgets(ReadLine, 1023, stream)==NULL){printf("ERROR 4 in reading HAMR file\n");return 1;}
      if(fgets(ReadLine, 1023, stream)==NULL){printf("ERROR 5 in reading HAMR file\n");return 1;}
      O->element[O->Nelement-1].NlabelAtom1=sscanf(ReadLine, "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s",
                                                   O->element[O->Nelement-1].LabelAtom1[0],
                                                   O->element[O->Nelement-1].LabelAtom1[1],
                                                   O->element[O->Nelement-1].LabelAtom1[2],
                                                   O->element[O->Nelement-1].LabelAtom1[3],
                                                   O->element[O->Nelement-1].LabelAtom1[4],
                                                   O->element[O->Nelement-1].LabelAtom1[5],
                                                   O->element[O->Nelement-1].LabelAtom1[6],
                                                   O->element[O->Nelement-1].LabelAtom1[7],
                                                   O->element[O->Nelement-1].LabelAtom1[8],
                                                   O->element[O->Nelement-1].LabelAtom1[9],
                                                   O->element[O->Nelement-1].LabelAtom1[10],
                                                   O->element[O->Nelement-1].LabelAtom1[11],
                                                   O->element[O->Nelement-1].LabelAtom1[12],
                                                   O->element[O->Nelement-1].LabelAtom1[13],
                                                   O->element[O->Nelement-1].LabelAtom1[14],
                                                   O->element[O->Nelement-1].LabelAtom1[15]
                                                   );
      if(O->element[O->Nelement-1].NlabelAtom1==-1)
      {
        printf("ERROR in reading lables for atom1 from HAMR file in Line:\n");
        printf("%s\n",ReadLine);
        return 1;
      }
      ContinueLabel=1;
      O->element[O->Nelement-1].NlabelAtom2=0;
      while(ContinueLabel)
      {
        if(fgets(ReadLine, 1023, stream)==NULL){printf("ERROR 6 in reading HAMR file\n");return 1;}
        if(sscanf(ReadLine, "%s %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE",
                  O->element[O->Nelement-1].LabelAtom2[O->element[O->Nelement-1].NlabelAtom2],
                  &(O->element[O->Nelement-1].H[0][O->element[O->Nelement-1].NlabelAtom2]),
                  &(O->element[O->Nelement-1].H[1][O->element[O->Nelement-1].NlabelAtom2]),
                  &(O->element[O->Nelement-1].H[2][O->element[O->Nelement-1].NlabelAtom2]),
                  &(O->element[O->Nelement-1].H[3][O->element[O->Nelement-1].NlabelAtom2]),
                  &(O->element[O->Nelement-1].H[4][O->element[O->Nelement-1].NlabelAtom2]),
                  &(O->element[O->Nelement-1].H[5][O->element[O->Nelement-1].NlabelAtom2]),
                  &(O->element[O->Nelement-1].H[6][O->element[O->Nelement-1].NlabelAtom2]),
                  &(O->element[O->Nelement-1].H[7][O->element[O->Nelement-1].NlabelAtom2]),
                  &(O->element[O->Nelement-1].H[8][O->element[O->Nelement-1].NlabelAtom2]),
                  &(O->element[O->Nelement-1].H[9][O->element[O->Nelement-1].NlabelAtom2]),
                  &(O->element[O->Nelement-1].H[10][O->element[O->Nelement-1].NlabelAtom2]),
                  &(O->element[O->Nelement-1].H[11][O->element[O->Nelement-1].NlabelAtom2]),
                  &(O->element[O->Nelement-1].H[12][O->element[O->Nelement-1].NlabelAtom2]),
                  &(O->element[O->Nelement-1].H[13][O->element[O->Nelement-1].NlabelAtom2]),
                  &(O->element[O->Nelement-1].H[14][O->element[O->Nelement-1].NlabelAtom2]),
                  &(O->element[O->Nelement-1].H[15][O->element[O->Nelement-1].NlabelAtom2])
                  )
           ==O->element[O->Nelement-1].NlabelAtom1+1)
        {
          O->element[O->Nelement-1].NlabelAtom2++;
          if(O->element[O->Nelement-1].NlabelAtom2==16)ContinueLabel=0;
        }
        else
        {
          strcpy(O->element[O->Nelement-1].LabelAtom2[O->element[O->Nelement-1].NlabelAtom2], "");
          ContinueLabel=0;
        }
      }
      
      
      
    }
    if(fgets(ReadLine, 1023, stream)==NULL)
    {
      ContinueLine=0;
    }
  }
  
  
  return 0;
}
uint32_t ReadStuttgartHAMR(StuttgartHAMR * O)
{
  return StreamReadStuttgartHAMR(O, stdin);
}
uint32_t FileReadStuttgartHAMR(StuttgartHAMR * O, const char *const flname)
{
  FILE  *stream;
  uint32_t retval;
  if ((stream = fopen(flname, "r")) == NULL)
  {
    printf("ERROR\nCan't open %s for reading.\n", flname);
    return 1;
  }
  retval = StreamReadStuttgartHAMR(O, stream);
  fclose(stream);
  return retval;
}
