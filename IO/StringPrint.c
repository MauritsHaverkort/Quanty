/*
 *
 */
int sprintfAdd(char ** A, char *fmt, ...)
{
  char * Atemp;
  int retval;
  unsigned long lenA, lenAlocA;
  va_list args;
  va_start ( args, fmt );
  lenA=strlen(*A);
  lenAlocA = 2024 * ((lenA/2024)+2);
  Atemp = realloc(*A, lenAlocA);
  if(Atemp==NULL)
  {
    printf("ERROR: realloc failed in sprintfAdd\n");
    return 1;
  }
  retval = vsprintf(Atemp+lenA,fmt,args);
  va_end(args);
  *A = Atemp;
  return retval;
}


/*
 *   Prints the values of the operator.
 *   Buffer needs to be a pointer to a list of characters that is not allocated.
 */
uint32_t StringPrintOperatorReadable(OperatorType O, char **Buffer)
{
  int32_t i, j, k;
  *Buffer = malloc(4048);
  if(*Buffer==NULL)
  {
    printf("ERROR: malloc failed in StringPrintOperatorReadable\n");
    return 1;
  }
  *Buffer[0]='\0';
  O.Name[MAXStrLength - 1] = '\0';
  sprintfAdd(Buffer, "\nOperator: %s\n", O.Name);
  sprintfAdd(Buffer, "QComplex         = %10i (Real==0 or Complex==1 or Mixed==2)\n", O.QComplex);
  sprintfAdd(Buffer, "MaxLength        = %10i (largest number of product of lader operators)\n", O.MaxLength);
  sprintfAdd(Buffer, "NFermionic modes = %10i (Number of fermionic modes (site, spin, orbital, ...) in the one particle basis)\n", O.NFermion);
  sprintfAdd(Buffer, "NBosonic modes   = %10i (Number of bosonic modes (phonon modes, ...) in the one particle basis)\n", O.NBoson);
  for (i = 0; i <= O.MaxLength; i++)
  {
    if (O.l[i].N != 0)
    {
      sprintfAdd(Buffer, "\nOperator of Length %3i\n", i);
      sprintfAdd(Buffer, "QComplex      = %10i (Real==0 or Complex==1)\n", O.l[i].QComplex);
      sprintfAdd(Buffer, "N             = %10i (number of operators of length %3i)\n", O.l[i].N, i);
      for (j = 0; j < O.l[i].N; j++)
      {
        for (k = i - 1; k >= 0; k--)
        {
          ((O.l[i].o[j * i + k]) & CrAnMask) ? sprintfAdd(Buffer, "C") : sprintfAdd(Buffer, "A");
          sprintfAdd(Buffer, "%3i ", (O.l[i].o[j * i + k]) & NotCrAnMask);
        }
        (O.l[i].QComplex) ? sprintfAdd(Buffer, "| %22.15E %22.15E\n", O.l[i].V[j].Re, O.l[i].V[j].Im) : sprintfAdd(Buffer, "| %22.15E\n", O.l[i].VRe[j]);
      }
    }
  }
  sprintfAdd(Buffer, "\n");
  return 0;
}

/*
 *   Prints the values of the wave function.
 */
uint32_t StringPrintWaveFunctionReadable(WaveFunctionType psi, char **Buffer)
{
  uint32_t i, j;
  *Buffer = malloc(4048);
  if(*Buffer==NULL)
  {
    printf("ERROR: malloc failed in StringPrintOperatorReadable\n");
    return 1;
  }
  *Buffer[0]='\0';
  psi.Name[MAXStrLength - 1] = '\0';
  sprintfAdd(Buffer, "\nWaveFunction: %s\n", psi.Name);
  sprintfAdd(Buffer, "QComplex         = %10i (Real==0 or Complex==1)\n", psi.QComplex);
  sprintfAdd(Buffer, "N                = %10i (Number of basis functions used to discribe psi)\n", psi.N);
  sprintfAdd(Buffer, "NFermionic modes = %10i (Number of fermions in the one particle basis)\n", psi.NFermion);
  sprintfAdd(Buffer, "NBosonic modes   = %10i (Number of bosons in the one particle basis)\n", psi.NBoson);
  if (psi.QComplex)
  {
    sprintfAdd(Buffer, "\n#      pre-factor             +I  pre-factor         Determinant\n");
    for (i = 1; i <= psi.N; i++)
    {
      sprintfAdd(Buffer, "%4i %20.12E       %20.12E       ", i, psi.ValueRe[BlokElement(i)], psi.ValueIm[BlokElement(i)]);
      for (j = 0; j < psi.NWord; j++)
      {
        sprintfAdd(Buffer, "%2.2X", psi.Determinant[BlokElementDet(i, psi.NWord, j)]);
      }
      sprintfAdd(Buffer, "\n");
    }
    sprintfAdd(Buffer, "\n");
  }
  else
  {
    sprintfAdd(Buffer, "\n#      pre-factor         Determinant\n");
    for (i = 1; i <= psi.N; i++)
    {
      sprintfAdd(Buffer, "%4i %20.12E       ", i, psi.ValueRe[BlokElement(i)]);
      for (j = 0; j < psi.NWord; j++)
      {
        sprintfAdd(Buffer, "%2.2X", psi.Determinant[BlokElementDet(i, psi.NWord, j)]);
      }
      sprintfAdd(Buffer, "\n");
    }
    sprintfAdd(Buffer, "\n");
  }
  return 0;
}

/*
 *  Prints an array of spectra
 */
uint32_t StringPrintSpectra(SpectraType * const Spec, const uint32_t NSpec, char **Buffer)
{
  uint32_t i, j, NE;
  *Buffer = malloc(4048);
  if(*Buffer==NULL)
  {
    printf("ERROR: malloc failed in StringPrintSpectra\n");
    return 1;
  }
  *Buffer[0]='\0';
  sprintfAdd(Buffer, "#Spectra: %i\n", NSpec);
  sprintfAdd(Buffer, "Emin______Emax     ");
  for (i = 0; i < NSpec; i++)
  {
    sprintfAdd(Buffer, " %22.15E %22.15E", Spec[i].Emin, Spec[i].Emax);
  }
  sprintfAdd(Buffer, "\n");
  sprintfAdd(Buffer, "EminPole__EmaxPole ");
  for (i = 0; i < NSpec; i++)
  {
    sprintfAdd(Buffer, " %22.15E %22.15E", Spec[i].EminPole, Spec[i].EmaxPole);
  }
  sprintfAdd(Buffer, "\n");
  sprintfAdd(Buffer, "dE________Gamma    ");
  for (i = 0; i < NSpec; i++)
  {
    sprintfAdd(Buffer, " %22.15E %22.15E", Spec[i].dE, Spec[i].Gamma);
  }
  sprintfAdd(Buffer, "\n");
  sprintfAdd(Buffer, "Phase_____0        ");
  for (i = 0; i < NSpec; i++)
  {
    sprintfAdd(Buffer, " %22.15E %22.15E", Spec[i].phase, 0.0);
  }
  sprintfAdd(Buffer, "\n");
  for (j = 0; j < 5; j++)
  {
    sprintfAdd(Buffer, "a[%1i]______b[%1i]     ", j, j);
    for (i = 0; i < NSpec; i++)
    {
      sprintfAdd(Buffer, " %22.15E %22.15E", Spec[i].a[j], Spec[i].b[j]);
    }
    sprintfAdd(Buffer, "\n");
  }
  sprintfAdd(Buffer, "Energy             ");
  for (i = 0; i < NSpec; i++)
  {
    if (i < 10)
      sprintfAdd(Buffer, "  Re[%1i]                  Im[%1i]                ", i, i);
    else if (i < 100)
      sprintfAdd(Buffer, "  Re[%2i]                 Im[%2i]               ", i, i);
    else if (i < 1000)
      sprintfAdd(Buffer, "  Re[%3i]                Im[%3i]              ", i, i);
    else if (i < 10000)
      sprintfAdd(Buffer, "  Re[%4i]               Im[%4i]             ", i, i);
  }
  sprintfAdd(Buffer, "\n");
  NE = Spec[0].N;
  for (i = 1; i < NSpec; i++)
  {
    NE = Max(NE, Spec[i].N);
  }
  for (j = 0; j <= NE; j++)
  {
    sprintfAdd(Buffer, "%19.12E", (double)(Spec[0].dE * j + Spec[0].Emin));
    for (i = 0; i < NSpec; i++)
    {
      if (j <= Spec[i].N)
        sprintfAdd(Buffer, " %22.15E %22.15E", Spec[i].Re[j], Spec[i].Im[j]);
      else
        sprintfAdd(Buffer, " %22.15E %22.15E", 0.0, 0.0);
    }
    sprintfAdd(Buffer, "\n");
  }
  return 0;
}

uint32_t StringPrintAddMatrix(MatrixType M, char **Buffer)
{
  uint32_t i, j;
  M.Name[MAXStrLength - 1] = '\0';
  sprintfAdd(Buffer, "%s =\n", M.Name);
  sprintfAdd(Buffer, "Real Part of Matrix with dimensions [Ni=%i][Nj=%i] ([Rows][Collums])\n", M.Ni, M.Nj);
  sprintfAdd(Buffer, "        ");
  for (j = 0; j < M.Nj; j++)
  {
    sprintfAdd(Buffer, "  [%12i]", j);
  }
  sprintfAdd(Buffer, "\n");
  for (i = 0; i < M.Ni; i++)
  {
    sprintfAdd(Buffer, "[%6i] ", i);
    for (j = 0; j < M.Nj; j++)
    {
      sprintfAdd(Buffer, "%15.8E ", M.VRe[i][j]);
    }
    sprintfAdd(Buffer, "\n");
  }
  if (M.QComplex)
  {
    sprintfAdd(Buffer, "Imaginary Part\n");
    sprintfAdd(Buffer, "        ");
    for (j = 0; j < M.Nj; j++)
    {
      sprintfAdd(Buffer, "  [%12i]", j);
    }
    sprintfAdd(Buffer, "\n");
    for (i = 0; i < M.Ni; i++)
    {
      sprintfAdd(Buffer, "[%6i] ", i);
      for (j = 0; j < M.Nj; j++)
      {
        sprintfAdd(Buffer, "%15.8E ", M.VIm[i][j]);
      }
      sprintfAdd(Buffer, "\n");
    }
  }
  sprintfAdd(Buffer, "\n");
  return 0;
}


/*  Prints a matrix
 */
uint32_t StringPrintMatrix(MatrixType M, char **Buffer)
{
  *Buffer = malloc(4048);
  if(*Buffer==NULL)
  {
    printf("ERROR: malloc failed in StringPrintTightBindingDef\n");
    return 1;
  }
  *Buffer[0]='\0';
  return StringPrintAddMatrix(M,Buffer);
}


/*
 *   Prints a crystal structure
 */
uint32_t StringPrintAddCrystalStructure(CrystalStructureType O, char **Buffer)
{
  uint32_t i,j,k;
  sprintfAdd(Buffer, "\nprintout of Crystal Structure\n");
  sprintfAdd(Buffer, "Units: ");
  if((O.Units&1)==0) sprintfAdd(Buffer, "NoPi (g.r=1)  ");
  if((O.Units&1)==1) sprintfAdd(Buffer, "2Pi (g.r=2Pi) ");
  if((O.Units&6)==0) sprintfAdd(Buffer, "Angstrom ");
  if((O.Units&6)==2) sprintfAdd(Buffer, "Bohr     ");
  if((O.Units&6)==4) sprintfAdd(Buffer, "nm       ");
  if((O.Units&8)==0) sprintfAdd(Buffer, "Absolute atom positions");
  if((O.Units&8)==8) sprintfAdd(Buffer, "Relative atom positions");
  sprintfAdd(Buffer, "\nUnit cell parameters:\n");
  sprintfAdd(Buffer, "a: %15.7F %15.7F %15.7F\n", O.rUnitCell[0][0], O.rUnitCell[0][1], O.rUnitCell[0][2]);
  sprintfAdd(Buffer, "b: %15.7F %15.7F %15.7F\n", O.rUnitCell[1][0], O.rUnitCell[1][1], O.rUnitCell[1][2]);
  sprintfAdd(Buffer, "c: %15.7F %15.7F %15.7F\n", O.rUnitCell[2][0], O.rUnitCell[2][1], O.rUnitCell[2][2]);
  sprintfAdd(Buffer, "Reciprocal latice:\n");
  sprintfAdd(Buffer, "a: %15.7F %15.7F %15.7F\n", O.gUnitCell[0][0], O.gUnitCell[0][1], O.gUnitCell[0][2]);
  sprintfAdd(Buffer, "b: %15.7F %15.7F %15.7F\n", O.gUnitCell[1][0], O.gUnitCell[1][1], O.gUnitCell[1][2]);
  sprintfAdd(Buffer, "c: %15.7F %15.7F %15.7F\n", O.gUnitCell[2][0], O.gUnitCell[2][1], O.gUnitCell[2][2]);
  sprintfAdd(Buffer, "Number of attoms %i\n", O.NAtoms);
  for (i = 0; i < O.NAtoms; i++)
  {
    sprintfAdd(Buffer, "#%4i | %s ( %i ) at position { %15.7F , %15.7F , %15.7F }\n", i, O.Atom[i].Name, O.Atom[i].AtomNumber, O.Atom[i].Position[0], O.Atom[i].Position[1], O.Atom[i].Position[2]);
    for(j=0;j<O.Atom[i].Nshell;j++)
    {
      sprintfAdd(Buffer, "      | %s shell with %i orbitals {",O.Atom[i].ShellName[j],O.Atom[i].shelldegeneracy[j]);
      for(k=0;k<O.Atom[i].shelldegeneracy[j];k++)
      {
        if(k!=0) sprintfAdd(Buffer, ",");
        sprintfAdd(Buffer, " %s ",O.Atom[i].OrbitalName[j][k]);
      }
      sprintfAdd(Buffer, "}\n");
    }
  }
  return 0;
}

/*
 *   Prints a crystal structure
 */
uint32_t StringPrintCrystalStructure(CrystalStructureType O, char **Buffer)
{
  *Buffer = malloc(4048);
  if(*Buffer==NULL)
  {
    printf("ERROR: malloc failed in StringPrintTightBindingDef\n");
    return 1;
  }
  *Buffer[0]='\0';
  return StringPrintAddCrystalStructure(O,Buffer);
}

/*  Prints the structure TightBindingDef
 *
 */
uint32_t StringPrintTightBindingDef(TightBindingDef O, char **Buffer)
{
  int    i;
  *Buffer = malloc(4048);
  if(*Buffer==NULL)
  {
    printf("ERROR: malloc failed in StringPrintTightBindingDef\n");
    return 1;
  }
  *Buffer[0]='\0';
  O.Name[MAXStrLength - 1] = '\0';
  sprintfAdd(Buffer, "\nSettigns of a tight binding model: %s\n", O.Name);
  StringPrintAddCrystalStructure(O.Cell,Buffer);
  sprintfAdd(Buffer, "Containing a total number of %i orbitals\n",O.NFermion);
  sprintfAdd(Buffer, "Hopping definitions ( %i )\n", O.nH);
  for (i = 0; i < O.nH; i++)
  {
    sprintfAdd(Buffer, "Hopping from %s - %s to %s - %s with translation vector in unit cells: { %i , %i , %i }\n", O.Cell.Atom[O.H[i].Atompos1].Name, O.Cell.Atom[O.H[i].Atompos1].ShellName[O.H[i].Shellpos1], O.Cell.Atom[O.H[i].Atompos2].Name, O.Cell.Atom[O.H[i].Atompos2].ShellName[O.H[i].Shellpos2], O.H[i].TransVecs[0], O.H[i].TransVecs[1], O.H[i].TransVecs[2]);
    StringPrintAddMatrix(O.H[i].Hop, Buffer);
    
  }
  sprintfAdd(Buffer, "\n");
  return 0;
}

