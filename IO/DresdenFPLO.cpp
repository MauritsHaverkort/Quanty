// Split a string into array of phrases according to the given delimiter.
uint32_t SplitStringToArray(char **array, char * str, const char * const del)
{
  uint32_t len = 0;
  char *p;

  for (p=strtok(str, del); p; p=strtok(NULL, del))
  {
      array[len++] = strdup(p);
  }

  return len;
}

int ShellNameToL(char name)
{
  switch(name)
  {
    case 's':
    case 'S':
      return 0;
    case 'p':
    case 'P':
      return 1;
    case 'd':
    case 'D':
      return 2;
    case 'f':
    case 'F':
      return 3;
    case 'g':
    case 'G':
      return 4;
    case 'h':
    case 'H':
      return 5;
    case 'i':
    case 'I':
      return 6;
    case 'k':
    case 'K':
      return 7;
    default:
      return -1;
  }
}



// pricision defined by the fplo output (coordinates).
double RoundToDigit(double x, int n)
{
  return roundf(x * pow(10,n)) / pow(10,n);
}

void TransVecNorm(const double vec[3], double gUC[3][3], double vecnorm[3])
{
  double nfac = 1.;
  vecnorm[0] = RoundToDigit(Dot(vec, gUC[0])/nfac, 5);
  vecnorm[1] = RoundToDigit(Dot(vec, gUC[1])/nfac, 5);
  vecnorm[2] = RoundToDigit(Dot(vec, gUC[2])/nfac, 5);
}

// Read Wannier orbitals from FPLO wan output.
// Only reads information necessary for creating a TB Hamiltonian.
uint32_t StreamReadDresdenFPLO(DresdenFPLO * O, FILE * stream)
{
  char ReadLine[1024];
  char *phrase[1024] = {NULL}, tmpstr[1024];
  double tmpvec[3];
  uint32_t Nphrase, Nwandef, i, j, ihop=0, iatom;
  int NhopMax=1024;
  O->lDOS.N=0;
  O->Lattice.Molecule=false;
  O->Lattice.Units=0; // Use Angstrom, g.r = 1, and atom positions witin the unit cell are in angstrom units
  ReadLine[1023] = '\0';
  for(i=0;i<3;i++)
  {
    for(j=0;j<3;j++)
    {
      O->Lattice.LattVec[i][j]=0;
      O->Lattice.gLattVec[i][j] = 0;
    }
    O->Lattice.LattVec[i][i]=INFINITY;
    O->Lattice.LattConst[i]=INFINITY;
  }
  while(fgets(ReadLine, 1023, stream) != NULL)
  {
    if(strstr(ReadLine, "Structure type: Molecule") != NULL)
    {
      O->Lattice.Molecule = true;
    }
    if(strstr(ReadLine, "Compound:") != NULL)
    {
      Nphrase = SplitStringToArray(phrase, ReadLine, " :\n");
      strncpy(O->Compound, phrase[1], 15);
      O->Compound[15]=0;
      for (uint p = 0; p < Nphrase; ++p) free(phrase[p]);
    }
    if(strstr(ReadLine, "lattice constants") != NULL)
    {
      Nphrase = SplitStringToArray(phrase, ReadLine, " :\n");
      for (i=0;i<3;i++)
      {
        O->Lattice.LattConst[i] = strtod(phrase[Nphrase-3+i],NULL);
      }
      for (uint p = 0; p < Nphrase; ++p) free(phrase[p]);
    }
    if(strstr(ReadLine, "axis angles") != NULL)
    {
      Nphrase = SplitStringToArray(phrase, ReadLine, " :\n");
      for (i=0;i<3;i++)
      {
        O->Lattice.Angles[i] = strtod(phrase[Nphrase-3+i],NULL);
      }
      for (uint p = 0; p < Nphrase; ++p) free(phrase[p]);
    }
    if( (strstr(ReadLine, "reciprocial lattice vectors") != NULL) || (strstr(ReadLine, "reciprocal lattice vectors") != NULL) )
    {
      for(i=0;i<3;i++)
      {
        fgets(ReadLine, 1023, stream);
        Nphrase = SplitStringToArray(phrase, ReadLine, " :\n");
        for(j=0;j<3;j++)
        {
          O->Lattice.gLattVec[i][j] = strtod(phrase[Nphrase-3+j],NULL);
        }
        for (uint p = 0; p < Nphrase; ++p) free(phrase[p]);
      }
    }
    else if(strstr(ReadLine, "lattice vectors") != NULL)
    {
      for(i=0;i<3;i++)
      {
        fgets(ReadLine, 1023, stream);
        Nphrase = SplitStringToArray(phrase, ReadLine, " :\n");
        for(j=0;j<3;j++)
        {
          O->Lattice.LattVec[i][j] = strtod(phrase[Nphrase-3+j],NULL);
        }
        for (uint p = 0; p < Nphrase; ++p) free(phrase[p]);
      }
    }
    // Read elements in the unit cell.
    if(strstr(ReadLine, "Number of sites :") != NULL)
    {
      Nphrase = SplitStringToArray(phrase, ReadLine, " :\n");
      O->Atoms.Natom = atoi(phrase[3]);
      for (uint p = 0; p < Nphrase; ++p) free(phrase[p]);
    }
    if(strstr(ReadLine, "No.  Element WPS CPA-Block") != NULL)
    {
      O->Atoms.Atom = (char**)malloc(O->Atoms.Natom*sizeof(char*));
      O->Atoms.Pos  = (double**)malloc(O->Atoms.Natom*sizeof(double*));
      for(i=0;i<O->Atoms.Natom;i++)
      {
        fgets(ReadLine, 1023, stream);
        Nphrase = SplitStringToArray(phrase, ReadLine, " \n");
        O->Atoms.Atom[i] = (char*)malloc(16*sizeof(char));
        strcpy(O->Atoms.Atom[i], phrase[1]);
        O->Atoms.Pos[i] = (double*)malloc(3*sizeof(double));
        for(j=0;j<3;j++)
        {
          O->Atoms.Pos[i][j] = strtod(phrase[4+j],NULL);
        }
        for (uint p = 0; p < Nphrase; ++p) free(phrase[p]);
      }
    }
    // Read the constructed Wannier orbitals.
    if(strstr(ReadLine, "Number of Wannier definitions:") != NULL)
    {
      Nphrase = SplitStringToArray(phrase, ReadLine, " :\n");
      Nwandef = atoi(phrase[Nphrase-1]);
      for (uint p = 0; p < Nphrase; ++p) free(phrase[p]);
      O->Wannier.Nwandef   = Nwandef;
      O->Wannier.NActAtom  = 1;
      O->Wannier.Name      = (char**)malloc(Nwandef*sizeof(char*));
      O->Wannier.Site      = (int*)malloc(Nwandef*sizeof(int));
      O->Wannier.Orbital   = (int**)malloc(Nwandef*sizeof(int*));
      O->Wannier.PosCenter = (double**)malloc(Nwandef*sizeof(double*));
      O->Wannier.QRelativistic = (uint32_t)(-1);
      i = 0;
      while((fgets(ReadLine, 1023, stream)!=NULL) && (!strstr(ReadLine, "START calculating WF")))
      {
        if(strstr(ReadLine, "name") != NULL)
        {
          O->Wannier.Name[i] = (char*)malloc(16*sizeof(char));
          Nphrase = SplitStringToArray(phrase, ReadLine, " :\n");
          strcpy(O->Wannier.Name[i], phrase[1]);
          for(j=2;j<Nphrase;j++)
            strcat(strcat(O->Wannier.Name[i], " "), phrase[j]);
          for (uint p = 0; p < Nphrase; ++p) free(phrase[p]);
          while(fgets(ReadLine, 1023, stream) != NULL)
          {
            if(strstr(ReadLine, "center") != NULL)
            {
              O->Wannier.PosCenter[i] = (double*)malloc(3*sizeof(double));
              Nphrase = SplitStringToArray(phrase, ReadLine, " :\n");
              for(j=0;j<3;j++)
              {
                tmpvec[j] = strtod(phrase[Nphrase-3+j],NULL);
              }
              TransVecNorm(tmpvec, O->Lattice.gLattVec, O->Wannier.PosCenter[i]);
              for (uint p = 0; p < Nphrase; ++p) free(phrase[p]);
            }
            if(strstr(ReadLine, "site") != NULL)
            {
              Nphrase = SplitStringToArray(phrase, ReadLine, " :\n");
              O->Wannier.Site[i] = atoi(phrase[1]);
              bool siteCameBefore = false;
							for (j = 0; j < i; ++j)
							{
								if (O->Wannier.Site[j] == O->Wannier.Site[i])
								{
									siteCameBefore = true;
								}
								else if (siteCameBefore)
								{
									std::cout << "Error: In its current structure Quanty expects all Wannier functions belonging to the same site to be clustered together." << endl;
									std::cout << "Currently, a definition on site " << O->Wannier.Site[i] << " is followed by a definition on site " << O->Wannier.Site[j] << ", but a bit later another definition on site "
											<< O->Wannier.Site[i] << " follows." << endl;
									std::cout << "The Quanty team is terribly sorry that you will have to do the downfolding anew." << endl;
									return 1;
								}
							}
							if (i > 0 && O->Wannier.Site[i] < O->Wannier.Site[i - 1])
							{
								std::cout << "Error: In its current structure Quanty expects the sites of all Wannier functions to be in ascending order." << endl;
								std::cout << "Currently, a definition on site " << O->Wannier.Site[i - 1] << " is followed by a definition on site " << O->Wannier.Site[i] << "." << endl;
								std::cout << "The Quanty team is terribly sorry that you will have to do the downfolding anew." << endl;
								return 1;
							}
              j=0;
              while(j<i && O->Wannier.Site[i]!=O->Wannier.Site[j]) j++;
              if(i>=1&&j==i&& O->Wannier.Site[i]!=O->Wannier.Site[i-1]) O->Wannier.NActAtom += 1;
              for (uint p = 0; p < Nphrase; ++p) free(phrase[p]);
            }
            if(strstr(ReadLine, "orbital") != NULL)
            {
              O->Wannier.Orbital[i] = (int*)malloc(3*sizeof(int));
              Nphrase = SplitStringToArray(phrase, ReadLine, " :\n");
							O->Wannier.Orbital[i][0] = atoi(&(phrase[1][0]));
							if (O->Wannier.QRelativistic == (uint32_t)(-1))
							{
								if (strstr(phrase[1], "/2") != NULL)
									O->Wannier.QRelativistic = (uint32_t) true;
								else
									O->Wannier.QRelativistic = (uint32_t) false;
							}
							if ((bool) O->Wannier.QRelativistic)
							{
								int l = ShellNameToL(phrase[1][1]);
								O->Wannier.Orbital[i][1] = l + (atoi(&(phrase[1][2])) - 1) / 2;
								strcpy(tmpstr, phrase[1]);
								O->Wannier.Orbital[i][2] = (atoi(tmpstr + 5) + 1)/2;
							}
							else
							{
								O->Wannier.Orbital[i][1] = ShellNameToL(phrase[1][1]);
								strcpy(tmpstr, phrase[1]);
								O->Wannier.Orbital[i][2] = atoi(tmpstr + 2);
							}
							for (uint p = 0; p < Nphrase; ++p) free(phrase[p]);
							break;
            }
          }

          i = i+1;
        }
      }
      if(i!=Nwandef)
      {
        printf("Error when reading Wannier definitions:\n  The number of Wannier definitions found is %d. It should be %d.\n", i, Nwandef);
        return 1;
      }

      // Read the hopping matrix.
      O->Wannier.QComplex  = (uint32_t)false;
      O->Wannier.Nhop   = (int*)calloc(Nwandef*Nwandef,sizeof(int));
      O->Wannier.HopVec = (double***)malloc(Nwandef*Nwandef*sizeof(double**));
      O->Wannier.HopVal = (double**)malloc(Nwandef*Nwandef*sizeof(double*));
      O->Wannier.HopValIm = (double**)malloc(Nwandef*Nwandef*sizeof(double*));
      i = -1;
      while((fgets(ReadLine, 1023, stream)!=NULL) && (!strstr(ReadLine, "END   calculating WF")))
      {
        if(strncmp(ReadLine, "spin", 4)==0)
        {
          if(++i>Nwandef*Nwandef)
          {
            printf("\nError: The length of hoppings is not compatible with the number of Wannier orbitals.\n");
            return 1;
          }
          O->Wannier.HopVec[i] = (double**)malloc(NhopMax*sizeof(double*));
          O->Wannier.HopVal[i] = (double*)malloc(NhopMax*sizeof(double));
          O->Wannier.HopValIm[i] = (double*)malloc(NhopMax*sizeof(double));
          ihop = 0;
        }
        if(strncmp(ReadLine, "T=", 2) ==0)
        {
          Nphrase = SplitStringToArray(phrase, ReadLine, " \n");
          O->Wannier.HopVec[i][ihop] = (double*)malloc(3*sizeof(double));
          for(j=0;j<3;j++)
          {
            tmpvec[j] = strtod(phrase[1+j], NULL);
          }
          if(O->Lattice.Molecule==true)
          {
            O->Wannier.HopVec[i][ihop][0]=tmpvec[0];
            O->Wannier.HopVec[i][ihop][1]=tmpvec[1];
            O->Wannier.HopVec[i][ihop][2]=tmpvec[2];
          }
          else
          {
            TransVecNorm(tmpvec, O->Lattice.gLattVec, O->Wannier.HopVec[i][ihop]);
          }
          O->Wannier.HopVal[i][ihop] = strtod(phrase[5], NULL);
          O->Wannier.HopValIm[i][ihop] = strtod(phrase[Nphrase-1], NULL);
					if (not O->Wannier.QComplex and Abs(O->Wannier.HopValIm[i][ihop]) > (1000 * DBL_MIN))
					{
						O->Wannier.QComplex = (uint32_t)true;
					}

          ihop = ihop+1;
          if(ihop > NhopMax)
          {
						NhopMax *= 2;
						O->Wannier.HopVec[i] = (double**) realloc(O->Wannier.HopVec[i], NhopMax * sizeof(double*));
						O->Wannier.HopVal[i] = (double*) realloc(O->Wannier.HopVal[i], NhopMax * sizeof(double));
						O->Wannier.HopValIm[i] = (double*) realloc(O->Wannier.HopValIm[i], NhopMax * sizeof(double));
						if (O->Wannier.HopVec[i] == NULL || O->Wannier.HopVal[i] == NULL || O->Wannier.HopValIm[i] == NULL)
						{
							cout << "\nError in ReadDresdenFPLO: Could not allocate memory for hopping to and from site " << (i + 1) << endl;
							return 1;
						}
					}
          else
          {
            O->Wannier.Nhop[i] = ihop;
          }
          for (uint p = 0; p < Nphrase; ++p) free(phrase[p]);
        }
      }
      if(i<Nwandef*Nwandef-1)
      {
        printf("\nError: The length of hoppings is not compatible with the number of Wannier orbitals.\n");
        return 1;
      }
    }
  }

  O->Wannier.ActMap = (int*)malloc(O->Atoms.Natom*sizeof(int));
  iatom = 0;
  for(i=0;i<O->Atoms.Natom;i++)
  {
    O->Wannier.ActMap[i] = -1;
    for(j=0;j<O->Wannier.Nwandef;j++)
    {
      if(i==O->Wannier.Site[j]-1)
      {
        O->Wannier.ActMap[i] = iatom++;
        break;
      }
    }
  }
  if(O->Lattice.Molecule==true)
  {
    for(i=0;i<3;i++)
    {
      for(j=0;j<3;j++)
      {
        O->Lattice.LattVec[i][j]=0;
        O->Lattice.gLattVec[i][j] = 0;
      }
      O->Lattice.LattVec[i][i]=INFINITY;
      O->Lattice.LattConst[i]=INFINITY;
    }
  }
  return 0;
}

uint32_t FileReadDresdenFPLO(DresdenFPLO * O, const char *const flname)
{
  FILE  *stream;
  uint32_t retval;
  if ((stream = fopen(flname, "r")) == NULL)
  {
    printf("ERROR\nCan't open %s for reading.\n", flname);
    return 1;
  }
  retval = StreamReadDresdenFPLO(O, stream);
  fclose(stream);

  return retval;
}

/*  Reads information form the FPLO =.in file
 *  Only selects information needed for the implemented featueres :-)
 *  The code starts general, but the lowest leavel reading is not
 */
uint32_t StreamReadDresdenFPLOisdotin(DresdenFPLO * O, FILE * stream)
{
  const int MAXSECTIONS = 32;
  char *sectionbuffertemp;
  char *sectionbuffer[MAXSECTIONS];
  char *variable;
  int isection;
  // Get the file size
  fseek(stream, 0, SEEK_END);
  long file_size = ftell(stream);
  rewind(stream);

  // Allocate memory for the file content
  char *buffer = (char *)malloc((file_size + 1) * sizeof(char));
  if (buffer==NULL)
  {
    printf("Could not allocate memory for buffer\n");
    return 1;
  }

  // Read the file into the buffer
  fread(buffer, sizeof(char), file_size, stream);
  buffer[file_size] = '\0';
  
  // loop over all sections and split the string according to the sections found
  // The result are MAXSECTIONS strings that point to the sections with the name in
  // sections[isections] and starting at the addres in sectionbuffer[isections]
  const char *sections[]={
    "pretext",
    "header",
    "structure_information",
    "structure_dependend",
    "mesh",
    "brillouin",
    "bandstructure_plot",
    "iteration_control",
    "forces",
    "options",
    "cpa",
    "LSDA_U",
    "OPC",
    "Advanced_output",
    "BasisSettings",
    "dHvA",
    "numerics",
    NULL};
  sectionbuffertemp = buffer;
  for(isection=1;isection<MAXSECTIONS;isection++)
  {
    sectionbuffer[isection]=NULL;
  }
  sectionbuffer[0]=buffer;
  while((sectionbuffertemp=strstr(sectionbuffertemp,"section "))!=NULL)
  {
    isection=0;
    while((sections[isection]!=NULL))
    {
      if(strncmp(sectionbuffertemp+8,sections[isection],strlen(sections[isection]))==0)
      {
        goto FoundSection;
      }
      isection++;
    }
    // Section not found
    {
      printf("Warning: Found a section with unkown name in file \"=.in\"\n");
      printf("%.40s\n",sectionbuffertemp);
      sectionbuffertemp[0]='\0';
      sectionbuffertemp++;
      continue;
    }
    FoundSection:
    // split the string and set the sectionbuffer pointer
    sectionbuffer[isection]=strchr(sectionbuffertemp,'{');
    if((sectionbuffer[isection])==NULL)
    {
      printf("ERROR, this should not happen\n");
      printf("Section %s found, but the contentents in \"=.in\" does not contain a {\n",sections[isection]);
      printf("You will now crash\n");
      sectionbuffertemp[0]='\0';
    }
    else
    {
      sectionbuffertemp[0]='\0';
      sectionbuffertemp++;
    }
  }
  // we now have splitted the file =.in into sections
  // next we need to split the sections into different variables
  isection=0;
  while((sections[isection]!=NULL))
  {
    RemoveCommentLines(sectionbuffer[isection], '#');
    RemoveDoubleSpacesAndNewlines(sectionbuffer[isection]);
    // the rest is coded without being general
    if(strcmp(sections[isection], "options")==0)
    {
      variable=strstr(sectionbuffer[isection]," spin =");
      if(variable!=NULL) sscanf(variable, "%*[^=] ={%d,", &(O->Options.MSpin));
      variable=strstr(sectionbuffer[isection]," xc_potential =");
      if(variable!=NULL) sscanf(variable, "%*[^=] ={%d,", &(O->Options.XCPotential));
      variable=strstr(sectionbuffer[isection]," relativistic =");
      if(variable!=NULL) sscanf(variable, "%*[^=] ={%d,", &(O->Options.Relativistic));
      variable=strstr(sectionbuffer[isection]," nuclear_model =");
      if(variable!=NULL) sscanf(variable, "%*[^=] ={%d,", &(O->Options.NuclearModel));
      variable=strstr(sectionbuffer[isection]," quantization_axis[3]");
      if(variable!=NULL) sscanf(variable, "%*[^=]={%lf,%lf,%lf", &(O->Options.QuantizationAxis[0]),&(O->Options.QuantizationAxis[1]),&(O->Options.QuantizationAxis[2]));
    }
    if(strcmp(sections[isection], "structure_information")==0)
    {
      int tempint=0;
      variable=strstr(sectionbuffer[isection]," structure_type =");
      if(variable!=NULL) sscanf(variable, "%*[^=] ={%d,", &tempint);
      O->Lattice.Molecule=false;
      if(tempint!=1)
      {
        O->Lattice.Molecule=true;
        O->Lattice.LattVec[0][0]=INFINITY; O->Lattice.LattVec[0][1]=0;        O->Lattice.LattVec[0][2]=0;
        O->Lattice.LattVec[1][0]=0;        O->Lattice.LattVec[1][1]=INFINITY; O->Lattice.LattVec[1][2]=0;
        O->Lattice.LattVec[2][0]=0;        O->Lattice.LattVec[2][1]=0;        O->Lattice.LattVec[2][2]=INFINITY;
        
        O->Lattice.gLattVec[0][0]=0;       O->Lattice.gLattVec[0][1]=0;       O->Lattice.gLattVec[0][2]=0;
        O->Lattice.gLattVec[1][0]=0;       O->Lattice.gLattVec[1][1]=0;       O->Lattice.gLattVec[1][2]=0;
        O->Lattice.gLattVec[2][0]=0;       O->Lattice.gLattVec[2][1]=0;       O->Lattice.gLattVec[2][2]=0;
      }
      // lattice vectors are not stored in =.in
      variable=strstr(sectionbuffer[isection]," lengthunit =");
      if(variable!=NULL) sscanf(variable, "%*[^=] ={%d,", &tempint);
      O->Lattice.Units=2;                  /* **** ***x 0 r.g = 1 | 1 r.g = 2 Pi */
      if(tempint==2) O->Lattice.Units=0;   /* **** *xx* 00 Angstrom | 01 a0 | 10 nm | 11 ?? */
                                           /* **** x*** 0 Positions in units described above | 1 Postions in terms of r[0] ,r[1] and r[2] */
      variable=strstr(sectionbuffer[isection]," lattice_constants[3]");
      if(variable!=NULL) sscanf(variable, "%*[^=]={%lf,%lf,%lf", &(O->Lattice.LattConst[0]),&(O->Lattice.LattConst[1]),&(O->Lattice.LattConst[2]));
      variable=strstr(sectionbuffer[isection]," axis_angles[3]");
      if(variable!=NULL) sscanf(variable, "%*[^=]={%lf,%lf,%lf", &(O->Lattice.Angles[0]),&(O->Lattice.Angles[1]),&(O->Lattice.Angles[2]));
    }
    if(strcmp(sections[isection], "bandstructure_plot")==0)
    {
      variable=strstr(sectionbuffer[isection]," bandplot_control =");
      if(variable!=NULL)
      {
        if(sscanf(variable, " bandplot_control ={%*c,%*c,%*d,%lf,%lf,%d,", &(O->lDOS.Emin),&(O->lDOS.Emax),&(O->lDOS.NE))!=3)
        {
          printf("Error in parsing bandplot_control\n");
        }
      }
    }
    isection++;
  }
  free(buffer);
  return 0;
}

/*  Read the structureal data from =.str_template
 *
 */
uint32_t StreamReadDresdenFPLOisdotstrtemplate(DresdenFPLO * O, FILE * stream)
{
  int i;
  const int buffSize = 1024;
  char ReadLine[buffSize];

  // Skip the first two lines
  if((fgets(ReadLine, buffSize, stream)==NULL) ||
     (fgets(ReadLine, buffSize, stream)==NULL)    )
  {
    perror("ERROR in reading =.str_template\nNot able to read the first two lines\n");
    return 1;
  }
  // Read the rUnitCell
  for(i=0;i<3;i++)
  {
    if(fgets(ReadLine, buffSize, stream)==NULL)
    {
      perror("ERROR in reading =.str_template\nNot able to read unit cell\n");
      return 1;
    }
    if(sscanf(ReadLine, "%lf %lf %lf", &(O->Lattice.LattVec[i][0]), &(O->Lattice.LattVec[i][1]), &(O->Lattice.LattVec[i][2]) )!=3)
    {
      perror("ERROR in reading =.str_template\nNot able to parse unit cell\n");
      return 1;
    }
  }
  // Generate the reciprocal Cell
  ReciprocalLatticeVector(O->Lattice.LattVec, O->Lattice.gLattVec, O->Lattice.Units);
  // Read the number of atoms
  if (fgets(ReadLine, buffSize, stream)==NULL)
  {
    perror("ERROR in reading =.str_template\nFailed to read NAtoms");
    return 1;
  }
  if (sscanf(ReadLine, "%d", &(O->Atoms.Natom)) != 1)
  {
    perror("ERROR in reading =.str_template\nFailed to parse NAtoms");
    return 1;
  }
  // allocate memory for atom description
  O->Atoms.Atom = (char**)malloc(O->Atoms.Natom * sizeof(char*));
  O->Atoms.Pos  = (double**)malloc(O->Atoms.Natom * sizeof(double*));
  for(i=0;i<O->Atoms.Natom;i++)
  {
    O->Atoms.Atom[i] = (char*)malloc(4 * sizeof(char));
    O->Atoms.Pos[i]  = (double*)malloc(3 * sizeof(double));
  }
  // skip line 7
  if(fgets(ReadLine, buffSize, stream)==NULL)
  {
    perror("ERROR in reading =.str_template\nNot able to read line 7\n");
    return 1;
  }
  for(i=0;i<O->Atoms.Natom;i++)
  {
    if(fgets(ReadLine, buffSize, stream)==NULL)
    {
      perror("ERROR in reading =.str_template\nNot able to read atomic positions\n");
      return 1;
    }
    if (sscanf(ReadLine, "%2s %*d %lf %lf %lf", O->Atoms.Atom[i],&(O->Atoms.Pos[i][0]),&(O->Atoms.Pos[i][0]),&(O->Atoms.Pos[i][0])) != 4)
    {
      perror("ERROR in reading =.str_template\nFailed to parse atom data\n");
      return 1;
    }
  }
  return 0;
}



uint32_t ReadDresdenFPLOlDOS(DresdenFPLO * O, const char *const dirname)
{
  const int buffSize = 1024;
  char ReadLine[buffSize];
  const int b = (int)strlen(dirname)+128;
  char flname[b];
  FILE  *stream;
  int i,isite,inl,iE;
  int site,n,m,l,spin;
  char lc;
  int MallocSize = 128;
  O->lDOS.sort = (int*)malloc(MallocSize*sizeof(int));
  O->lDOS.n    = (int*)malloc(MallocSize*sizeof(int));
  O->lDOS.l    = (int*)malloc(MallocSize*sizeof(int));
  O->lDOS.nl   = (int*)malloc(MallocSize*sizeof(int));
  O->lDOS.m    = (int*)malloc(MallocSize*sizeof(int));
  O->lDOS.spin = (int*)malloc(MallocSize*sizeof(int));
  O->lDOS.N = 0;
  for(isite=1;isite<=O->Atoms.Natom;isite++)
  {
    inl=1;
    sprintf(flname, "%s/+ldos.site%03d.nl%03d",dirname,isite,inl);
    while ((stream = fopen(flname, "r")) != NULL)
    {
      while(fgets(ReadLine, buffSize, stream) != NULL)
      {
        if(ReadLine[0]=='#')
        {
          if (sscanf(ReadLine, "# site = %d nl = %d%c m =  %d spin = %d", &site,&n,&lc,&m,&spin) != 5)
          {
            perror("ERROR in reading lDOS\n");
            printf("ERROR: %s\n",ReadLine);
            return 1;
          }
          l=ShellNameToL(lc);
          if(site!=isite)
          {
            printf("ERROR in reading lDOS\n");
            printf("%s\n",ReadLine);
            printf("But istite is %d\n",isite);
          }
          O->lDOS.sort[O->lDOS.N] = site;
          O->lDOS.n[O->lDOS.N]=n;
          O->lDOS.l[O->lDOS.N]=l;
          O->lDOS.nl[O->lDOS.N]=inl;
          O->lDOS.m[O->lDOS.N]=m;
          O->lDOS.spin[O->lDOS.N]=spin;
          O->lDOS.N++;
          if(O->lDOS.N==MallocSize)
          {
            MallocSize+=MallocSize;
            O->lDOS.sort = (int*)realloc(O->lDOS.sort, MallocSize*sizeof(int));
            O->lDOS.n    = (int*)realloc(O->lDOS.n,    MallocSize*sizeof(int));
            O->lDOS.l    = (int*)realloc(O->lDOS.l,    MallocSize*sizeof(int));
            O->lDOS.nl   = (int*)realloc(O->lDOS.nl,   MallocSize*sizeof(int));
            O->lDOS.m    = (int*)realloc(O->lDOS.m,    MallocSize*sizeof(int));
            O->lDOS.spin = (int*)realloc(O->lDOS.spin, MallocSize*sizeof(int));
          }
        }
      }
      fclose(stream);
      inl++;
      sprintf(flname, "%s/+ldos.site%03d.nl%03d",dirname,isite,inl);
    }
  }
  if(O->lDOS.N == 0)
  {
    free(O->lDOS.sort);
    free(O->lDOS.n);
    free(O->lDOS.l);
    free(O->lDOS.m);
    free(O->lDOS.spin);
    return 0;
  }
  O->lDOS.V = (double**)malloc(O->lDOS.N*sizeof(double*));
  for(i=0;i<O->lDOS.N;i++)
  {
    O->lDOS.V[i] = (double*)malloc(O->lDOS.NE*sizeof(double));
  }
  i=0;
  for(isite=1;isite<=O->Atoms.Natom;isite++)
  {
    inl=1;
    sprintf(flname, "%s/+ldos.site%03d.nl%03d",dirname,isite,inl);
    while ((stream = fopen(flname, "r")) != NULL)
    {
      while(fgets(ReadLine, buffSize, stream) != NULL)
      {
        if(ReadLine[0]=='#')
        {
          for(iE=0;iE<O->lDOS.NE;iE++)
          {
            fscanf(stream,"%*lf %lf",&(O->lDOS.V[i][iE]));
          }
          i++;
        }
      }
      fclose(stream);
      inl++;
      sprintf(flname, "%s/+ldos.site%03d.nl%03d",dirname,isite,inl);
    }
  }
  return 0;
}

uint32_t ReadFPLO(DresdenFPLO * O, const char *const dirname)
{
  O->Atoms.Natom=0;
  O->Wannier.Nwandef=0;
  O->Wannier.NActAtom=0;
  O->lDOS.N=0;
  const int b = (int)strlen(dirname)+128;
  char flname[b];
  FILE  *stream;
  uint32_t retval=0;
  // read =.in file needs to be present
  sprintf(flname, "%s/=.in", dirname);
  if ((stream = fopen(flname, "r")) == NULL)
  {
    printf("ERROR\nCan't open %s for reading.\n", flname);
    printf("fopen failed with error number %d (%s)\n", errno, strerror(errno));
    return 1;
  }
  if(StreamReadDresdenFPLOisdotin(O, stream)==1)
  {
    printf("ERROR: Reading =.in failed\n");
    return 1;
  }
  fclose(stream);
  
  // read =.str_template file needs to be present
  sprintf(flname, "%s/=.str_template", dirname);
  if ((stream = fopen(flname, "r")) == NULL)
  {
    printf("ERROR\nCan't open %s for reading.\n", flname);
    printf("fopen failed with error number %d (%s)\n", errno, strerror(errno));
    return 1;
  }
  if(StreamReadDresdenFPLOisdotstrtemplate(O, stream)==1)
  {
    printf("ERROR: Reading =.str_template failed\n");
    return 1;
  }
  fclose(stream);
  // read +ldos.sortxxx.nlyyy if present
  ReadDresdenFPLOlDOS(O, dirname);

  
  
  
  // read +fval.xxx.y if present
  
  // read +fcor.xxx.y if present

  // read +hamdata if present
  
  
  

  return retval;
}
