// c++ like implementation
//TODO: Anyone have a better idea where to put this?
template <typename T>
string to_string(const T value)
{
  std::ostringstream os ;
  os << value ;
  return os.str() ;
}

// c implementation
// sprintfAdd takes a string and prints to the end of this string
//
// you find 3 implementation.
//
// if you do not know the length of the string, nor the amount of bytes allocated you need to
// determine the length and reallocate -- 2 times slow
//


/*
 *
 */
int sprintfAdd(char ** A, const char *fmt, ...)
{
  // The function strlen makes this function slow. It can be used once, but not to consequitive build a long string
  char * Atemp;
  int retval;
  unsigned long lenA, lenAlocA;
  va_list args;
  va_start ( args, fmt );
  lenA=strlen(*A);
  lenAlocA = 2024 * ((lenA/2024)+2);
  Atemp = (char*)realloc(*A, lenAlocA);
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

// We have string A->str of length A->length
// We add the string defined by *fmt of length extralen
//    ... It is the users respomnsibility to make sure that this is ok
// If A->length + extralen is smaller than a->allocated we need to realloc the string
int sprintfAdd(DynamicString *A, int extralen, const char *fmt, ...)
{
  char * Atemp;
  int retval;
  va_list args;
  va_start ( args, fmt );
  if(A->allocated < A->length+extralen)
  {
    A->allocated = ( A->allocated + extralen ) * 2;
    Atemp = (char*)realloc(A->str, A->allocated+1);
    if(Atemp==NULL)
    {
      printf("ERROR: realloc failed in sprintfAdd\n");
      return 1;
    }
    A->str = Atemp;
  }
  retval = vsnprintf(A->str+A->length,extralen,fmt,args);
  A->length = A->length + extralen;
  va_end(args);
  return retval;
}


int sprintfAdd(DynamicString *A, const char *fmt, ...)
{
  const int MaxBuffSize = 1024;
  char * Atemp;
  int extralen;
  va_list args;
  va_start ( args, fmt );
  if(A->allocated < A->length+MaxBuffSize)
  {
    A->allocated = Max(A->length+MaxBuffSize,2*(A->allocated));
    Atemp = (char*)realloc(A->str, A->allocated+1);
    if(Atemp==NULL)
    {
      printf("ERROR: realloc failed in sprintfAdd\n");
      return 1;
    }
    A->str = Atemp;
  }
  extralen = vsnprintf(A->str+A->length,MaxBuffSize,fmt,args);
  A->length = A->length + extralen;
  va_end(args);
  return 0;
}




/*
 *   Prints the values of the operator.
 *   StrAdrr needs to be a pointer to a list of characters that is not allocated.
 */
uint32_t StringPrintOperatorReadable(OperatorType O, char **StrAdrr)
{
  int32_t i, j, k;
  DynamicString temp;
  DynamicString *Buffer;
  temp = InitDynamicString();
  Buffer = &temp;
//  *Buffer = (char*)malloc(4048);
//  if(*Buffer==NULL)
//  {
//    printf("ERROR: malloc failed in StringPrintOperatorReadable\n");
//    return 1;
//  }
//  *Buffer[0]='\0';
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
        (O.l[i].QComplex) ? sprintfAdd(Buffer, "| %21.14E %21.14E\n", O.l[i].V[j].Re, O.l[i].V[j].Im) : sprintfAdd(Buffer, "| %21.14E\n", O.l[i].VRe[j]);
      }
    }
  }
  sprintfAdd(Buffer, "\n");
  *StrAdrr = Buffer->str;
  return 0;
}

/*
 *   Prints the values of the wave function.
 */
uint32_t StringPrintWaveFunctionReadable(WaveFunctionType psi, char **Buffer)
{
  uint32_t i, j, StringLength,DecaLengNpsi;
  // Determine the needed size of the Buffer
  // Changing text in the header means one needs to change the allocation size
  i=1;
  for(DecaLengNpsi=1;DecaLengNpsi<=10;DecaLengNpsi++)
  {
    i*=10;
    if(psi.N<i) break;
  }
  StringLength = MAXStrLength+304; // Length of Header
  if(psi.QComplex)
  {
    StringLength+=56+DecaLengNpsi;
    StringLength+=psi.N*(4+40+DecaLengNpsi+psi.NFermion+3*psi.NBoson+1);
  }
  else
  {
    StringLength+=35+DecaLengNpsi;
    StringLength+=psi.N*(3+20+DecaLengNpsi+psi.NFermion+3*psi.NBoson+1);
  }
  *Buffer = (char*)malloc(StringLength);
  if(*Buffer==NULL)
  {
    printf("ERROR: malloc failed in StringPrintOperatorReadable\n");
    return 1;
  }
  *Buffer[0]='\0';
  StringLength=0;
  psi.Name[MAXStrLength - 1] = '\0';
  sprintf(*Buffer+StringLength, "\nWaveFunction: %s\n", psi.Name);
  StringLength+=16+strlen(psi.Name);
  sprintf(*Buffer+StringLength, "QComplex         = %10i (Real==0 or Complex==1)\n", psi.QComplex);
  StringLength+=54;
  sprintf(*Buffer+StringLength, "N                = %10i (Number of basis functions used to discribe psi)\n", psi.N);
  StringLength+=79;
  sprintf(*Buffer+StringLength, "NFermionic modes = %10i (Number of fermions in the one particle basis)\n", psi.NFermion);
  StringLength+=77;
  sprintf(*Buffer+StringLength, "NBosonic modes   = %10i (Number of bosons in the one particle basis)\n", psi.NBoson);
  StringLength+=75;
  
  if (psi.QComplex)
  {
    sprintf(*Buffer+StringLength, "\n");
    StringLength+=1;
    for (i = 1; i<=DecaLengNpsi; i++ )
    {
      sprintf(*Buffer+StringLength, "#");
      StringLength++;
    }
    sprintf(*Buffer+StringLength, "   pre-factor           pre-factor          Determinant\n");
    StringLength+=56;
    
    for (i = 1; i <= psi.N; i++)
    {
      switch(DecaLengNpsi)
      {
        case 1  : sprintf(*Buffer+StringLength,  "%1i %20.12E %20.12E  ", i, psi.ValueRe[BlokElement(i)], psi.ValueIm[BlokElement(i)]); break;
        case 2  : sprintf(*Buffer+StringLength,  "%2i %20.12E %20.12E  ", i, psi.ValueRe[BlokElement(i)], psi.ValueIm[BlokElement(i)]); break;
        case 3  : sprintf(*Buffer+StringLength,  "%3i %20.12E %20.12E  ", i, psi.ValueRe[BlokElement(i)], psi.ValueIm[BlokElement(i)]); break;
        case 4  : sprintf(*Buffer+StringLength,  "%4i %20.12E %20.12E  ", i, psi.ValueRe[BlokElement(i)], psi.ValueIm[BlokElement(i)]); break;
        case 5  : sprintf(*Buffer+StringLength,  "%5i %20.12E %20.12E  ", i, psi.ValueRe[BlokElement(i)], psi.ValueIm[BlokElement(i)]); break;
        case 6  : sprintf(*Buffer+StringLength,  "%6i %20.12E %20.12E  ", i, psi.ValueRe[BlokElement(i)], psi.ValueIm[BlokElement(i)]); break;
        case 7  : sprintf(*Buffer+StringLength,  "%7i %20.12E %20.12E  ", i, psi.ValueRe[BlokElement(i)], psi.ValueIm[BlokElement(i)]); break;
        case 8  : sprintf(*Buffer+StringLength,  "%8i %20.12E %20.12E  ", i, psi.ValueRe[BlokElement(i)], psi.ValueIm[BlokElement(i)]); break;
        case 9  : sprintf(*Buffer+StringLength,  "%9i %20.12E %20.12E  ", i, psi.ValueRe[BlokElement(i)], psi.ValueIm[BlokElement(i)]); break;
        case 10 : sprintf(*Buffer+StringLength, "%10i %20.12E %20.12E  ", i, psi.ValueRe[BlokElement(i)], psi.ValueIm[BlokElement(i)]); break;
      }
      StringLength+=44+DecaLengNpsi;
      //sprintfAdd(Buffer, "%9i %20.12E  %20.12E  ", i, psi.ValueRe[BlokElement(i)], psi.ValueIm[BlokElement(i)]);
      //for (j = 0; j < psi.NWord; j++)
      //{
      //  sprintfAdd(Buffer, "%2.2X", psi.Determinant[BlokElementDet(i, psi.NWord, j)]);
      //}
      
      for(j=0;j<psi.NFermion;j++)
      {
        sprintf(*Buffer+StringLength, "%c",(psi.Determinant[BlokElementDet(i, psi.NWord, (j >> 3))] & (128 >> (j & 7)))?'1':'0');
        StringLength++;
      }
      for(j=psi.NWord-psi.NBoson;j<psi.NWord;j++)
      {
        sprintf(*Buffer+StringLength, " %2.2X", psi.Determinant[BlokElementDet(i, psi.NWord, j)]);
        StringLength+=3;
      }
      sprintf(*Buffer+StringLength, "\n");
      StringLength++;
    }
    sprintf(*Buffer+StringLength, "\n");
    StringLength++;
  }
  else
  {
    sprintf(*Buffer+StringLength, "\n");
    StringLength+=1;
    for (i = 1; i<=DecaLengNpsi; i++ )
    {
      sprintf(*Buffer+StringLength, "#");
      StringLength++;
    }
    sprintf(*Buffer+StringLength, "   pre-factor          Determinant\n");
    StringLength+=35;
    for (i = 1; i <= psi.N; i++)
    {
      switch(DecaLengNpsi)
      {
        case 1  : sprintf(*Buffer+StringLength,  "%1i %20.12E  ", i, psi.ValueRe[BlokElement(i)]); break;
        case 2  : sprintf(*Buffer+StringLength,  "%2i %20.12E  ", i, psi.ValueRe[BlokElement(i)]); break;
        case 3  : sprintf(*Buffer+StringLength,  "%3i %20.12E  ", i, psi.ValueRe[BlokElement(i)]); break;
        case 4  : sprintf(*Buffer+StringLength,  "%4i %20.12E  ", i, psi.ValueRe[BlokElement(i)]); break;
        case 5  : sprintf(*Buffer+StringLength,  "%5i %20.12E  ", i, psi.ValueRe[BlokElement(i)]); break;
        case 6  : sprintf(*Buffer+StringLength,  "%6i %20.12E  ", i, psi.ValueRe[BlokElement(i)]); break;
        case 7  : sprintf(*Buffer+StringLength,  "%7i %20.12E  ", i, psi.ValueRe[BlokElement(i)]); break;
        case 8  : sprintf(*Buffer+StringLength,  "%8i %20.12E  ", i, psi.ValueRe[BlokElement(i)]); break;
        case 9  : sprintf(*Buffer+StringLength,  "%9i %20.12E  ", i, psi.ValueRe[BlokElement(i)]); break;
        case 10 : sprintf(*Buffer+StringLength, "%10i %20.12E  ", i, psi.ValueRe[BlokElement(i)]); break;
      }
      StringLength+=23+DecaLengNpsi;
      
      //for (j = 0; j < psi.NWord; j++)
      //{
      //  sprintfAdd(Buffer, "%2.2X", psi.Determinant[BlokElementDet(i, psi.NWord, j)]);
      //}
      for(j=0;j<psi.NFermion;j++)
      {
        sprintf(*Buffer+StringLength, "%c",(psi.Determinant[BlokElementDet(i, psi.NWord, (j >> 3))] & (128 >> (j & 7)))?'1':'0');
        StringLength++;
      }
      for(j=psi.NWord-psi.NBoson;j<psi.NWord;j++)
      {
        sprintf(*Buffer+StringLength, " %2.2X", psi.Determinant[BlokElementDet(i, psi.NWord, j)]);
        StringLength+=3;
      }
      sprintf(*Buffer+StringLength, "\n");
      StringLength++;
    }
    sprintf(*Buffer+StringLength, "\n");
    StringLength++;
  }
  return 0;
}

/*
 *  Prints an array of spectra
 */
uint32_t StringPrintSpectra(SpectraType * const Spec, const uint32_t NSpec, char **Buffer)
{
  uint32_t i, j, NE;
  *Buffer = (char*)malloc(4048);
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
//  sprintfAdd(Buffer, "Phase_____0        ");
//  for (i = 0; i < NSpec; i++)
//  {
//    sprintfAdd(Buffer, " %22.15E %22.15E", Spec[i].phase, 0.0);
//  }
//  sprintfAdd(Buffer, "\n");
//  for (j = 0; j < 5; j++)
//  {
//    sprintfAdd(Buffer, "a[%1i]______b[%1i]     ", j, j);
//    for (i = 0; i < NSpec; i++)
//    {
//      sprintfAdd(Buffer, " %22.15E %22.15E", Spec[i].a[j], Spec[i].b[j]);
//    }
//    sprintfAdd(Buffer, "\n");
//  }
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
  *Buffer = (char*)malloc(4048);
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
  sprintfAdd(Buffer, "Number of atoms %i\n", O.NAtoms);
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
  *Buffer = (char*)malloc(4048);
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
  *Buffer = (char*)malloc(4048);
  if(*Buffer==NULL)
  {
    printf("ERROR: malloc failed in StringPrintTightBindingDef\n");
    return 1;
  }
  *Buffer[0]='\0';
  O.Name[MAXStrLength - 1] = '\0';
  sprintfAdd(Buffer, "\nSettings of a tight binding model: %s\n", O.Name);
  StringPrintAddCrystalStructure(O.Cell,Buffer);
  sprintfAdd(Buffer, "Containing a total number of %i orbitals\n",O.NFermion);
  sprintfAdd(Buffer, "Hopping definitions ( %i )\n", O.nH);
  for (i = 0; i < O.nH; i++)
  {
    sprintfAdd(Buffer, "Hopping from %i : %s - %s to %i : %s - %s with translation vector in unit cells: { %i , %i , %i } ({%15.8E %15.8E %15.8E })\n", O.H[i].Atompos1, O.Cell.Atom[O.H[i].Atompos1].Name,O.Cell.Atom[O.H[i].Atompos1].ShellName[O.H[i].Shellpos1], O.H[i].Atompos2, O.Cell.Atom[O.H[i].Atompos2].Name, O.Cell.Atom[O.H[i].Atompos2].ShellName[O.H[i].Shellpos2], O.H[i].TransVecs[0], O.H[i].TransVecs[1], O.H[i].TransVecs[2], O.H[i].ConnectingVector[0], O.H[i].ConnectingVector[1], O.H[i].ConnectingVector[2]);
    StringPrintAddMatrix(O.H[i].Hop, Buffer);

  }
  sprintfAdd(Buffer, "\n");
  return 0;
}

