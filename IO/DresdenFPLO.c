// Split a string into array of phrases according to the given delimiter.
uint32_t SplitStringToArray(char **array, char *str, char *del)
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
    default:
      return -1;
  }
}

// pricision defined by the fplo output (coordinates).
double RoundToDigit(double x, int n)
{
  return roundf(x * pow(10,n)) / pow(10,n);
}

void TransVecNorm(const double vec[3], const double gUC[3][3], double vecnorm[3])
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
  uint32_t Nphrase, Nwandef, Nimag=0, i, j, ihop=0, iatom;
  const int NhopMax=256;

  ReadLine[1023] = '\0';

  while(fgets(ReadLine, 1023, stream) != NULL)
  {
    // Read structure information.
    if(strstr(ReadLine, "Compound:") != NULL)
    {
      Nphrase = SplitStringToArray(phrase, ReadLine, " :\n");
      strncpy(O->Compound, phrase[1], 15);
      O->Compound[15]=0;
    }
    if(strstr(ReadLine, "lattice constants") != NULL)
    {
      Nphrase = SplitStringToArray(phrase, ReadLine, " :\n");
      for (i=0;i<3;i++)
      {
        O->Lattice.LattConst[i] = strtod(phrase[Nphrase-3+i],NULL);
      }
    }
    if(strstr(ReadLine, "axis angles") != NULL)
    {
      Nphrase = SplitStringToArray(phrase, ReadLine, " :\n");
      for (i=0;i<3;i++)
      {
        O->Lattice.Angles[i] = strtod(phrase[Nphrase-3+i],NULL);
      }
    }
    if(strstr(ReadLine, "reciprocial lattice vectors") != NULL)
    {
      for(i=0;i<3;i++)
      {
        fgets(ReadLine, 1023, stream);
        Nphrase = SplitStringToArray(phrase, ReadLine, " :\n");
        for(j=0;j<3;j++)
        {
          O->Lattice.gLattVec[i][j] = strtod(phrase[Nphrase-3+j],NULL);
        }
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
      }
    }
    // Read elements in the unit cell.
    if(strstr(ReadLine, "Number of sites :") != NULL)
    {
      Nphrase = SplitStringToArray(phrase, ReadLine, " :\n");
      O->Atoms.Natom = atoi(phrase[3]);
    }
    if(strstr(ReadLine, "No.  Element WPS CPA-Block") != NULL)
    {
      O->Atoms.Atom = malloc(O->Atoms.Natom*sizeof(char*));
      O->Atoms.Pos  = malloc(O->Atoms.Natom*sizeof(double*));
      for(i=0;i<O->Atoms.Natom;i++)
      {
        fgets(ReadLine, 1023, stream);
        Nphrase = SplitStringToArray(phrase, ReadLine, " \n");
        O->Atoms.Atom[i] = malloc(16*sizeof(char));
        strcpy(O->Atoms.Atom[i], phrase[1]);
        O->Atoms.Pos[i] = malloc(3*sizeof(double));
        for(j=0;j<3;j++)
        {
          O->Atoms.Pos[i][j] = strtod(phrase[4+j],NULL);
        }
      }
    }
    // Read the constructed Wannier orbitals.
    if(strstr(ReadLine, "Number of Wannier definitions:") != NULL)
    {
      Nphrase = SplitStringToArray(phrase, ReadLine, " :\n");
      Nwandef = atoi(phrase[Nphrase-1]);
      O->Wannier.Nwandef   = Nwandef;
      O->Wannier.NActAtom  = 1;
      O->Wannier.Name      = malloc(Nwandef*sizeof(char*));
      O->Wannier.Site      = malloc(Nwandef*sizeof(int));
      O->Wannier.Orbital   = malloc(Nwandef*sizeof(int*));
      O->Wannier.PosCenter = malloc(Nwandef*sizeof(double*));
      i = 0;
      while((fgets(ReadLine, 1023, stream)!=NULL) && (!strstr(ReadLine, "START calculating WF")))
      {
        if(strstr(ReadLine, "name") != NULL)
        {
          O->Wannier.Name[i] = malloc(16*sizeof(char));
          Nphrase = SplitStringToArray(phrase, ReadLine, " :\n");
          strcpy(O->Wannier.Name[i], phrase[1]);
          for(j=2;j<Nphrase;j++)
            strcat(strcat(O->Wannier.Name[i], " "), phrase[j]);
          while(fgets(ReadLine, 1023, stream) != NULL)
          {
            if(strstr(ReadLine, "center") != NULL)
            {
              O->Wannier.PosCenter[i] = malloc(3*sizeof(double));
              Nphrase = SplitStringToArray(phrase, ReadLine, " :\n");
              for(j=0;j<3;j++)
              {
                tmpvec[j] = strtod(phrase[Nphrase-3+j],NULL);
              }
              TransVecNorm(tmpvec, O->Lattice.gLattVec, O->Wannier.PosCenter[i]);
            }
            if(strstr(ReadLine, "site") != NULL)
            {
              Nphrase = SplitStringToArray(phrase, ReadLine, " :\n");
              O->Wannier.Site[i] = atoi(phrase[1]);
              j=0;
              while(j<i && O->Wannier.Site[i]!=O->Wannier.Site[j]) j++;
              if(i>=1&&j==i&& O->Wannier.Site[i]!=O->Wannier.Site[i-1]) O->Wannier.NActAtom += 1;
            }
            if(strstr(ReadLine, "orbital") != NULL)
            {
              O->Wannier.Orbital[i] = malloc(3*sizeof(int));
              Nphrase = SplitStringToArray(phrase, ReadLine, " :\n");
              O->Wannier.Orbital[i][0] = atoi(&(phrase[1][0]));
              O->Wannier.Orbital[i][1] = ShellNameToL(phrase[1][1]);
              strcpy(tmpstr, phrase[1]);
              O->Wannier.Orbital[i][2] = atoi(tmpstr+2);
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
      O->Wannier.Nhop   =  malloc(Nwandef*Nwandef*sizeof(int));
      O->Wannier.HopVec =  malloc(Nwandef*Nwandef*sizeof(double**));
      O->Wannier.HopVal =  malloc(Nwandef*Nwandef*sizeof(double*));
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
          O->Wannier.HopVec[i] = malloc(NhopMax*sizeof(double*));
          O->Wannier.HopVal[i] = malloc(NhopMax*sizeof(double));
          ihop = 0;
        }
        if(strncmp(ReadLine, "T=", 2) ==0)
        {
          Nphrase = SplitStringToArray(phrase, ReadLine, " \n");
          O->Wannier.HopVec[i][ihop] = malloc(3*sizeof(double));
          for(j=0;j<3;j++)
          {
            tmpvec[j] = strtod(phrase[1+j], NULL);
          }
          TransVecNorm(tmpvec, O->Lattice.gLattVec, O->Wannier.HopVec[i][ihop]);
          O->Wannier.HopVal[i][ihop] = strtod(phrase[5], NULL);
          if(strtod(phrase[Nphrase-1], NULL)!=0.0) Nimag+=1;

          ihop = ihop+1;
          if(ihop > NhopMax)
          {
            printf("Warning: The length of matrix element is larger than the allocated matrix (%d).\nThe rest of the hoppings are ignored.\n", NhopMax);
          }
          else
          {
            O->Wannier.Nhop[i] = ihop;
          }
        }
      }
      if(i<Nwandef*Nwandef-1)
      {
        printf("\nError: The length of hoppings is not compatible with the number of Wannier orbitals.\n");
        return 1;
      }
      if(Nimag>0)
      {
        printf("Warning: There are %d hopping elements with non-zero imaginary part.\nThey are NOT read in.\n", Nimag);
      }
    }
  }

  O->Wannier.ActMap = malloc(O->Atoms.Natom*sizeof(int));
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
