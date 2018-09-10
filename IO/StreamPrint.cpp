/*
 *  Prints an array of spectra
 */
uint32_t StreamPrintSpectra(SpectraType * const Spec, const uint32_t NSpec, FILE * stream)
{
  uint32_t i, j, NE;
  fprintf(stream, "#Spectra: %i\n", NSpec);
  fprintf(stream, "Emin______Emax     ");
  for (i = 0; i < NSpec; i++)
  {
    fprintf(stream, " %22.15E %22.15E", Spec[i].Emin, Spec[i].Emax);
  }
  fprintf(stream, "\n");
  fprintf(stream, "EminPole__EmaxPole ");
  for (i = 0; i < NSpec; i++)
  {
    fprintf(stream, " %22.15E %22.15E", Spec[i].EminPole, Spec[i].EmaxPole);
  }
  fprintf(stream, "\n");
  fprintf(stream, "dE________Gamma    ");
  for (i = 0; i < NSpec; i++)
  {
    fprintf(stream, " %22.15E %22.15E", Spec[i].dE, Spec[i].Gamma);
  }
  fprintf(stream, "\n");
//  fprintf(stream, "Phase_____0        ");
//  for (i = 0; i < NSpec; i++)
//  {
//    fprintf(stream, " %22.15E %22.15E", Spec[i].phase, 0.0);
//  }
//  fprintf(stream, "\n");
//  for (j = 0; j < 5; j++)
//  {
//    fprintf(stream, "a[%1i]______b[%1i]     ", j, j);
//    for (i = 0; i < NSpec; i++)
//    {
//      fprintf(stream, " %22.15E %22.15E", Spec[i].a[j], Spec[i].b[j]);
//    }
//    fprintf(stream, "\n");
//  }
  fprintf(stream, "Energy             ");
  for (i = 0; i < NSpec; i++)
  {
    if (i < 10)
      fprintf(stream, "  Re[%1i]                  Im[%1i]                ", i, i);
    else if (i < 100)
      fprintf(stream, "  Re[%2i]                 Im[%2i]               ", i, i);
    else if (i < 1000)
      fprintf(stream, "  Re[%3i]                Im[%3i]              ", i, i);
    else if (i < 10000)
      fprintf(stream, "  Re[%4i]               Im[%4i]             ", i, i);
  }
  fprintf(stream, "\n");
  NE = Spec[0].N;
  for (i = 1; i < NSpec; i++)
  {
    NE = Max(NE, Spec[i].N);
  }
  for (j = 0; j <= NE; j++)
  {
    fprintf(stream, "%19.12E", (double)(Spec[0].dE * j + Spec[0].Emin));
    for (i = 0; i < NSpec; i++)
    {
      if (j <= Spec[i].N)
        fprintf(stream, " %22.15E %22.15E", Spec[i].Re[j], Spec[i].Im[j]);
      else
        fprintf(stream, " %22.15E %22.15E", 0.0, 0.0);
    }
    fprintf(stream, "\n");
  }
  return 0;
}

/*  by Robert
 *
 */
uint32_t StreamPrintRIXSSpectra(const uint32_t numwRes, SpectraType **Spec, double *wRes, const uint32_t NSpec, FILE * stream)
{
  uint32_t i, j, p, NE;
  fprintf(stream, "Energy           ");
  for(i = 0; i < NSpec; i++)
  {
    /*Header for real part*/
    for(j = 0; j < numwRes; j++)
    {
      fprintf(stream, "%f          ", wRes[j]);
    }
    /*Header for imag part*/
    for(j = 0; j < numwRes; j++)
    {
      fprintf(stream, "%f          ", wRes[j]);
    }

  }
  fprintf(stream, "\n");

  NE = Spec[0][0].N;
  for (i = 1; i < NSpec; i++)
  {
    for(j = 0; j < numwRes; j++)
      NE = Max(NE, Spec[i][j].N);
  }
  for (j = 0; j <= NE; j++)
  {
    fprintf(stream, "%19.12E", (double)(Spec[0][0].dE * j + Spec[0][0].Emin));
    for (i = 0; i < NSpec; i++)
    {

      for(p = 0; p < numwRes; p++)
      {
        if(j <= Spec[i][p].N)
          fprintf(stream, " %22.15E ", Spec[i][p].Re[j]);
        else
          fprintf(stream, " %22.15E ", 0.0);
      }
      for(p = 0; p < numwRes; p++)
      {
        if(j <= Spec[i][p].N)
          fprintf(stream, " %22.15E ", Spec[i][p].Im[j]);
        else
          fprintf(stream, " %22.15E ", 0.0);
      }

    }
    fprintf(stream, "\n");
  }
  return 0;
}

/*  by Robert
 *
 */
uint32_t StreamPrintRIXSSpectrum(const uint32_t numwRes, SpectraType **Spec, double *wRes, int printReal, const uint32_t iSpec, FILE * stream)
{
  //printReal == 1 : print the real part to the file
  //printReal == other : print the imaginary part

  uint32_t j, p, NE;
  fprintf(stream, "Energy           ");
  /*Header */
  for(j = 0; j < numwRes; j++)
  {
    fprintf(stream, "%f          ", wRes[j]);
  }

  fprintf(stream, "\n");

  NE = Spec[iSpec][0].N;

  for (j = 0; j <= NE; j++)
  {
    fprintf(stream, "%19.12E", (double)(Spec[iSpec][0].dE * j + Spec[iSpec][0].Emin));
    if(printReal == 1)
    {
      for(p = 0; p < numwRes; p++)
      {
        fprintf(stream, " %22.15E ", Spec[iSpec][p].Re[j]);
      }
    }
    else
    {
      for(p = 0; p < numwRes; p++)
      {
        fprintf(stream, " %22.15E ", -Spec[iSpec][p].Im[j]);
      }
    }
    fprintf(stream, "\n");
  }
  return 0;
}



/*
*   Prints the values of the wave function.
 */
uint32_t StreamPrintWaveFunctionPlain(WaveFunctionType psi, FILE * stream)
{
  uint32_t i, j;
  psi.Name[MAXStrLength - 1] = '\0';
  fprintf(stream, "%s\n", psi.Name);
  fprintf(stream, "%i\n", psi.QComplex);
  fprintf(stream, "%i\n", psi.N);
  fprintf(stream, "%i\n", psi.NFermion);
  fprintf(stream, "%i\n", psi.NBoson);
  if (psi.QComplex)
  {
    for (i = 1; i <= psi.N; i++)
    {
      fprintf(stream, "%23.15E %23.15E ", psi.ValueRe[BlokElement(i)], psi.ValueIm[BlokElement(i)]);
      for (j = 0; j < psi.NWord; j++)
      {
        fprintf(stream, "%2.2X", psi.Determinant[BlokElementDet(i, psi.NWord, j)]);
      }
      fprintf(stream, "\n");
    }
    fprintf(stream, "\n");
  }
  else
  {
    for (i = 1; i <= psi.N; i++)
    {
      fprintf(stream, "%23.15E ", psi.ValueRe[BlokElement(i)]);
      for (j = 0; j < psi.NWord; j++)
      {
        fprintf(stream, "%2.2X", psi.Determinant[BlokElementDet(i, psi.NWord, j)]);
      }
      fprintf(stream, "\n");
    }
    fprintf(stream, "\n");
  }
  return 0;
}


/*
*   Prints the values of the wave function.
 */
uint32_t StreamPrintWaveFunctionPlainBinary(WaveFunctionType psi, FILE * stream)
{
  uint32_t i;
  psi.Name[MAXStrLength - 1] = '\0';
  fprintf(stream, "%s\n", psi.Name);
  fprintf(stream, "%i\n", psi.QComplex);
  fprintf(stream, "%i\n", psi.N);
  fprintf(stream, "%i\n", psi.NFermion);
  fprintf(stream, "%i\n", psi.NBoson);
  if (psi.QComplex)
  {
    for (i = 1; i <= psi.N; i++)
    {
      fwrite(&(psi.ValueRe[BlokElement(i)]), sizeof(double), 1, stream);
      fwrite(&(psi.ValueIm[BlokElement(i)]), sizeof(double), 1, stream);
      fwrite(&(psi.Determinant[BlokElementDet(i, psi.NWord, 0)]), sizeof(uint8_t), psi.NWord, stream);
    }
  }
  else
  {
    for (i = 1; i <= psi.N; i++)
    {
      fwrite(&(psi.ValueRe[BlokElement(i)]), sizeof(double), 1, stream);
      fwrite(&(psi.Determinant[BlokElementDet(i, psi.NWord, 0)]), sizeof(uint8_t), psi.NWord, stream);
    }
  }
  return 0;
}



/*
 *   Prints the values of the wave function.
 */
uint32_t StreamPrintWaveFunctionReadable(WaveFunctionType psi, FILE * stream)
{
  uint32_t i, j;
  psi.Name[MAXStrLength - 1] = '\0';
  fprintf(stream, "\nWaveFunction: %s\n", psi.Name);
  fprintf(stream, "QComplex         = %10i (Real==0 or Complex==1)\n", psi.QComplex);
  fprintf(stream, "N                = %10i (Number of basis functions used to discribe psi)\n", psi.N);
  fprintf(stream, "NFermionic modes = %10i (Number of fermions in the one particle basis)\n", psi.NFermion);
  fprintf(stream, "NBosonic modes   = %10i (Number of bosons in the one particle basis)\n", psi.NBoson);
  if (psi.QComplex)
  {
    fprintf(stream, "\n#      pre-factor             +I  pre-factor         Determinant\n");
    for (i = 1; i <= psi.N; i++)
    {
      fprintf(stream, "%4i %20.12E       %20.12E       ", i, psi.ValueRe[BlokElement(i)], psi.ValueIm[BlokElement(i)]);
      for (j = 0; j < psi.NWord; j++)
      {
        fprintf(stream, "%2.2X", psi.Determinant[BlokElementDet(i, psi.NWord, j)]);
      }
      fprintf(stream, "\n");
    }
    fprintf(stream, "\n");
  }
  else
  {
    fprintf(stream, "\n#      pre-factor         Determinant\n");
    for (i = 1; i <= psi.N; i++)
    {
      fprintf(stream, "%4i %20.12E       ", i, psi.ValueRe[BlokElement(i)]);
      for (j = 0; j < psi.NWord; j++)
      {
        fprintf(stream, "%2.2X", psi.Determinant[BlokElementDet(i, psi.NWord, j)]);
      }
      fprintf(stream, "\n");
    }
    fprintf(stream, "\n");
  }
  return 0;
}


/*
*   Prints the values of the wave function.
 */
uint32_t StreamPrintWaveFunctionReadableTruncated(WaveFunctionType psi, double epsilon, FILE * stream)
{
  uint32_t i, j;
  double epssqr = Sqr(epsilon);
  psi.Name[MAXStrLength - 1] = '\0';
  fprintf(stream, "\nprintout of WaveFunction: %s\n", psi.Name);
  fprintf(stream, "Only determinants with a norm of the prefactor larger than %22.15E will be printed\n", epsilon);
  fprintf(stream, "\nGenneral settings\n");
  fprintf(stream, "QComplex      = %10i (Real==0 or Complex==1)\n", psi.QComplex);
  fprintf(stream, "N             = %10i (Number of basis functions used to discribe psi)\n", psi.N);
  fprintf(stream, "NFermion      = %10i (Number of fermions in the one particle basis)\n", psi.NFermion);
  fprintf(stream, "NBoson        = %10i (Number of bosons in the one particle basis)\n", psi.NBoson);
  if (psi.QComplex)
  {
    fprintf(stream, "\n  pre-factor             +I  pre-factor         Determinant\n");
    for (i = 1; i <= psi.N; i++)
    {
      if ((Sqr(psi.ValueRe[BlokElement(i)]) + Sqr(psi.ValueIm[BlokElement(i)])) > epssqr)
      {
        fprintf(stream, "%20.12E       %20.12E       ", psi.ValueRe[BlokElement(i)], psi.ValueIm[BlokElement(i)]);
        for (j = 0; j < psi.NWord; j++)
        {
          fprintf(stream, "%2.2X", psi.Determinant[BlokElementDet(i, psi.NWord, j)]);
        }
        fprintf(stream, "\n");
      }
    }
    fprintf(stream, "\n");
  }
  else
  {
    fprintf(stream, "\n  pre-factor         Determinant\n");
    for (i = 1; i <= psi.N; i++)
    {
      if (fabs(psi.ValueRe[BlokElement(i)]) > epsilon)
      {
        fprintf(stream, "%20.12E       ", psi.ValueRe[BlokElement(i)]);
        for (j = 0; j < psi.NWord; j++)
        {
          fprintf(stream, "%2.2X", psi.Determinant[BlokElementDet(i, psi.NWord, j)]);
        }
        fprintf(stream, "\n");
      }
    }
    fprintf(stream, "\n");
  }
  return 0;
}


/*
 *   Prints the values of the wave function, including information on memory allocation.
 */
uint32_t StreamPrintWaveFunctionFullReadable(WaveFunctionType psi, FILE * stream)
{
  uint32_t i, j;
  psi.Name[MAXStrLength - 1] = '\0';
  fprintf(stream, "\nFull printout of WaveFunction: %s\n", psi.Name);
  fprintf(stream, "\nGenneral settings\n");
  fprintf(stream, "QComplex      = %10i (Real==0 or Complex==1)\n", psi.QComplex);
  fprintf(stream, "N             = %10i (Number of basis functions used to discribe psi)\n", psi.N);
  fprintf(stream, "NBlok         = %10i (Number of bloks needed to store psi)\n", psi.NBlok);
  fprintf(stream, "NFermion      = %10i (Number of fermions in the one particle basis)\n", psi.NFermion);
  fprintf(stream, "NBoson        = %10i (Number of bosons in the one particle basis)\n", psi.NBoson);
  fprintf(stream, "NWord         = %10i (NFermion/8)+NBoson\n", psi.NWord);
  fprintf(stream, "NBitsKey      = %10i (Size of the key used for the hash table lookup)\n", psi.NBitsKey);
  fprintf(stream, "\nPointer positions (addresses)\n");
  fprintf(stream, "**ValueRe     = %p\n", psi.ValueRe);
  if (psi.QComplex)
    fprintf(stream, "**ValueIm     = %p\n", psi.ValueIm);
  fprintf(stream, "**Next        = %p\n", psi.Next);
  fprintf(stream, "**Determinant = %p\n", psi.Determinant);
  for (i = 0; i < psi.NBlok; i++)
  {
    fprintf(stream, "BlokNr.       = %i\n", i);
    fprintf(stream, "*ValueRe      = %p\n", psi.ValueRe[i]);
    if (psi.QComplex)
      fprintf(stream, "*ValueIm      = %p\n", psi.ValueIm[i]);
    fprintf(stream, "*Next         = %p\n", psi.Next[i]);
    fprintf(stream, "*Determinant  = %p\n", psi.Determinant[i]);
  }
  fprintf(stream, "\nNon-zero elements of Hash Table\n");
  for (i = 0; i < ((uint32_t) 1 << (psi.NBitsKey)); i++)
  {
    if (psi.KeyList[i] != 0)
    {
      fprintf(stream, "i= %10i , key= %10i\n", i, psi.KeyList[i]);
    }
  }
  if (psi.QComplex)
  {
    fprintf(stream,
            "\n  pre-factor             +I  pre-factor            Pointer Addres  Next_Hash Hash_Value   Determinant\n");
    for (i = 1; i <= psi.N; i++)
    {
      fprintf(stream, "%20.12E       %20.12E       %p %10i %10i   ", psi.ValueRe[BlokElement(i)],
              psi.ValueIm[BlokElement(i)], &psi.ValueRe[BlokElement(i)], psi.Next[BlokElement(i)],
              hash(&psi.Determinant[BlokElementDet(i, psi.NWord, 0)], psi.NWord, psi.NBitsKey));
      for (j = 0; j < psi.NWord; j++)
      {
        fprintf(stream, "%2.2X", psi.Determinant[BlokElementDet(i, psi.NWord, j)]);
      }
      fprintf(stream, "\n");
    }
    fprintf(stream, "\n");
  }
  else
  {
    fprintf(stream, "\n  pre-factor            Pointer Addres  Next_Hash Hash_Value   Determinant\n");
    for (i = 1; i <= psi.N; i++)
    {
      fprintf(stream, "%20.12E       %p %10i %10i   ", psi.ValueRe[BlokElement(i)], &psi.ValueRe[BlokElement(i)],
              psi.Next[BlokElement(i)], hash(&psi.Determinant[BlokElementDet(i, psi.NWord, 0)], psi.NWord,
                                             psi.NBitsKey));
      for (j = 0; j < psi.NWord; j++)
      {
        fprintf(stream, "%2.2X", psi.Determinant[BlokElementDet(i, psi.NWord, j)]);
      }
      fprintf(stream, "\n");
    }
    fprintf(stream, "\n");
  }
  return 0;
}

/*
 *   prints the onsite density of a wave-function.
 *   These results depend on the basis choosen, which can be quite arbetrary
 */
uint32_t StreamPrintWaveFunctionDensity(WaveFunctionType const psi, FILE * stream)
{
  uint32_t i;
  double Rho[psi.NBoson + psi.NFermion];
  CalculateDensity(Rho, psi);
  printf("Onsite Density (-only the diagonal of the density matrix-)\n");
  for (i = 0; i < psi.NBoson + psi.NFermion; i++)
  {
    printf("%5i %22.15E\n", i, Rho[i]);
  }
  printf("\n");
  return 0;
}


/*
 *   prints the correlation tensor of a wavefunction for a given set of fermions
 */
uint32_t StreamPrintWaveFunctionPartialCorrelationTensor(WaveFunctionType const psi, FILE * stream, uint32_t const NFermionIncluded, uint16_t * const FermionIncluded)
{
  uint32_t i,j,k,l,m=0;
  Complex16 * Rho ;
  Rho= (Complex16*)calloc(NFermionIncluded*NFermionIncluded*NFermionIncluded*NFermionIncluded,sizeof(Complex16));
  if (Rho==NULL)
  {
    printf("Faild to allocate memory for Rho in StreamPrintWaveFunctionPartialCorrelationTensor\n");
    fflush(stdout);
    return 1;
  }
  if(CalculatePartialDensityCorrelationTensor(Rho, psi, NFermionIncluded, FermionIncluded))
  {
    printf("CalculatePartialDensityCorrelationTensor failed in StreamPrintWaveFunctionPartialCorrelationTensor\n");
    return 1;
  }
  printf("Correlation function for %i Fermions\n",NFermionIncluded);
  for (l = 0; l < NFermionIncluded; l++)
  {
    for (k = 0; k < NFermionIncluded; k++)
    {
      for (j = 0; j < NFermionIncluded; j++)
      {
        for (i = 0; i < NFermionIncluded; i++)
        {
          printf("< C%3i C%3i A%3i A%3i > = %22.15E + I %22.15E\n",i,j,k,l,Rho[m].Re,Rho[m].Im);
          m++;
        }
      }
    }
  }
  printf("\n");
  free(Rho);
  return 0;
}


/*
 *   prints a Histogram lisitng the prefactors of a wave-funciton
 */
uint32_t StreamPrintWaveFunctionPrefactorHistogram(WaveFunctionType const psi, const double exp, const uint32_t N, FILE * stream)
{
  uint32_t i;
  double Hist[N];
  CalculatePrefactorHistogram(Hist, psi, exp, N);
  printf("\nHistogram of prefactors in Wave-function\n");
  for (i = 0; i < N - 1; i++)
  {
    printf("(%f)^(-%2i) to (%f)^(-%2i) | %22.15E\n", exp, i, exp, i + 1, Hist[i]);
  }
  printf("(%f)^(-%2i) to ( %f)       | %22.15E\n", exp, i, 0.0, Hist[N - 1]);
  printf("\n");
  return 0;
}


/*
 *   Prints the values of the operator.
 */
uint32_t StreamPrintOperatorReadable(OperatorType O, FILE * stream)
{
  int32_t i, j, k;
  O.Name[MAXStrLength - 1] = '\0';
  fprintf(stream, "\nOperator: %s\n", O.Name);
  fprintf(stream, "QComplex         = %10i (Real==0 or Complex==1 or Mixed==2)\n", O.QComplex);
  fprintf(stream, "MaxLength        = %10i (largest number of product of lader operators)\n", O.MaxLength);
  fprintf(stream, "NFermionic modes = %10i (Number of fermions in the one particle basis)\n", O.NFermion);
  fprintf(stream, "NBosonic modes   = %10i (Number of bosons in the one particle basis)\n", O.NBoson);
  for (i = 0; i <= O.MaxLength; i++)
  {
    if (O.l[i].N != 0)
    {
      fprintf(stream, "\nOperator of Length %3i\n", i);
      fprintf(stream, "QComplex      = %10i (Real==0 or Complex==1)\n", O.l[i].QComplex);
      fprintf(stream, "N             = %10i (number of operators of length %3i)\n", O.l[i].N, i);
      for (j = 0; j < O.l[i].N; j++)
      {
        for (k = i - 1; k >= 0; k--)
        {
          ((O.l[i].o[j * i + k]) & CrAnMask) ? fprintf(stream, "C") : fprintf(stream, "A");
          fprintf(stream, "%3i ", (O.l[i].o[j * i + k]) & NotCrAnMask);
        }
        (O.l[i].QComplex) ? fprintf(stream, "| %22.15E %22.15E\n", O.l[i].V[j].Re, O.l[i].V[j].Im) : fprintf(stream, "| %22.15E\n", O.l[i].VRe[j]);
      }
    }
  }
  fprintf(stream, "\n");
  return 0;
}


/*
*   Prints the values of the operator including information on allocation.
 */
uint32_t StreamPrintOperatorFullReadable(OperatorType O, FILE * stream)
{
  int32_t i, j, k;
  O.Name[MAXStrLength - 1] = '\0';
  fprintf(stream, "\nprintout of Operator: %s\n", O.Name);
  fprintf(stream, "\nGenneral settings\n");
  fprintf(stream, "QComplex      = %10i (Real==0 or Complex==1 or Mixed==2)\n", O.QComplex);
  fprintf(stream, "MaxLength     = %10i (largest number of product of lader operators)\n", O.MaxLength);
  fprintf(stream, "NFermion      = %10i (Number of fermions in the one particle basis)\n", O.NFermion);
  fprintf(stream, "NBoson        = %10i (Number of bosons in the one particle basis)\n", O.NBoson);
  fprintf(stream, "NWord         = %10i (Number of fermions/8+ number of bosons in the one particle basis)\n", O.NWord);
  fprintf(stream, "*l            = %p   (Pointer addres of the Fixed Lenght Operator Type structures)\n", O.l);
  for (i = 0; i <= O.MaxLength; i++)
  {
    fprintf(stream, "\nOperator of Length %3i\n", i);
    fprintf(stream, "&l[%3i]       = %p   (Pointer addres of fixed length variable)\n", i, &O.l[i]);
    fprintf(stream, "QComplex      = %10i (Real==0 or Complex==1)\n", O.l[i].QComplex);
    fprintf(stream, "N             = %10i (number of operators of length %3i)\n", O.l[i].N, i);
    fprintf(stream, "NBitsKey      = %10i (Number of bits in the key)\n", O.l[i].NBitsKey);
    fprintf(stream, "NAloc         = %10i (Size of the allocated memory)\n", O.l[i].NAloc);
    fprintf(stream, "*o            = %p   (Pointer addres of the determinant list)\n", O.l[i].o);
    fprintf(stream, "*Next         = %p   (Pointer addres of the list of next operators with the same hash value)\n",
            O.l[i].Next);
    fprintf(stream, "*KeyList      = %p   (Pointer addres of the KeyList)\n", O.l[i].KeyList);
    fprintf(stream, "*VRe          = %p   (Pointer addres of the list with Real values for real operators)\n",
            O.l[i].VRe);
    fprintf(stream, "*V            = %p   (Pointer addres of the list of complex values for complex operators)\n",
            O.l[i].V);

    if ((i != 0) && (O.l[i].NAloc != 0))
    {
      fprintf(stream, "\nNon-zero elements of Hash Table\n");
      for (j = 0; j < ((uint32_t) 1 << (O.l[i].NBitsKey)); j++)
      {
        if (O.l[i].KeyList[j] != 0)
        {
          fprintf(stream, "i= %10i , key= %10i\n", j, O.l[i].KeyList[j]);
        }
      }
    }
    for (j = 0; j < O.l[i].N; j++)
    {
      fprintf(stream, "%4i | ", j + 1);
      for (k = i - 1; k >= 0; k--)
      {
        ((O.l[i].o[j * i + k]) & CrAnMask) ? fprintf(stream, "C") : fprintf(stream, "A");
        fprintf(stream, "%3i ", (O.l[i].o[j * i + k]) & NotCrAnMask);
      }
      (O.l[i].QComplex) ? fprintf(stream, "| %22.15E %22.15E", O.l[i].V[j].Re, O.l[i].V[j].Im) : fprintf(stream,
                                                                                                         "| %22.15E",
                                                                                                         O.l[i].VRe[j]);
      if (i != 0)
        fprintf(stream, " | -> %10i\n", O.l[i].Next[j]);
    }
  }
  return 0;
}


/*
 *   Prints a crystal structure
 */
uint32_t StreamPrintCrystalStructure(CrystalStructureType O, FILE * stream)
{
  uint32_t i,j,k;
  fprintf(stream, "\nprintout of Crystal Structure\n");
  fprintf(stream, "Units: ");
  if((O.Units&1)==0) fprintf(stream, "NoPi (g.r=1)  ");
  if((O.Units&1)==1) fprintf(stream, "2Pi (g.r=2Pi) ");
  if((O.Units&6)==0) fprintf(stream, "Angstrom ");
  if((O.Units&6)==2) fprintf(stream, "Bohr     ");
  if((O.Units&6)==4) fprintf(stream, "nm       ");
  if((O.Units&8)==0) fprintf(stream, "Absolute atom positions");
  if((O.Units&8)==8) fprintf(stream, "Relative atom positions");
  fprintf(stream, "\nUnit cell parameters:\n");
  fprintf(stream, "a: %15.7F %15.7F %15.7F\n", O.rUnitCell[0][0], O.rUnitCell[0][1], O.rUnitCell[0][2]);
  fprintf(stream, "b: %15.7F %15.7F %15.7F\n", O.rUnitCell[1][0], O.rUnitCell[1][1], O.rUnitCell[1][2]);
  fprintf(stream, "c: %15.7F %15.7F %15.7F\n", O.rUnitCell[2][0], O.rUnitCell[2][1], O.rUnitCell[2][2]);
  fprintf(stream, "Reciprocal latice:\n");
  fprintf(stream, "a: %15.7F %15.7F %15.7F\n", O.gUnitCell[0][0], O.gUnitCell[0][1], O.gUnitCell[0][2]);
  fprintf(stream, "b: %15.7F %15.7F %15.7F\n", O.gUnitCell[1][0], O.gUnitCell[1][1], O.gUnitCell[1][2]);
  fprintf(stream, "c: %15.7F %15.7F %15.7F\n", O.gUnitCell[2][0], O.gUnitCell[2][1], O.gUnitCell[2][2]);
  fprintf(stream, "Number of atoms %i\n", O.NAtoms);

  for (i = 0; i < O.NAtoms; i++)
  {
    fprintf(stream, "#%4i | %s ( %i ) at position { %15.7F , %15.7F , %15.7F }\n", i, O.Atom[i].Name, O.Atom[i].AtomNumber, O.Atom[i].Position[0], O.Atom[i].Position[1], O.Atom[i].Position[2]);
    for(j=0;j<O.Atom[i].Nshell;j++)
    {
      fprintf(stream, "      | %s shell with %i orbitals {",O.Atom[i].ShellName[j],O.Atom[i].shelldegeneracy[j]);
      for(k=0;k<O.Atom[i].shelldegeneracy[j];k++)
      {
        if(k!=0) fprintf(stream, ",");
        fprintf(stream, " %s ",O.Atom[i].OrbitalName[j][k]);
      }
      fprintf(stream, "}\n");
    }
  }
  return 0;
}


/*
*   Prints a crystal structure super cell
 */
uint32_t StreamPrintCrystalStructureSuperCell(CrystalStructureSuperCellType O, FILE * stream)
{
  uint32_t i;
  fprintf(stream, "\nprintout of Super Cell\n");
  CrystalStructureType Otemp;
  CreateCrystalStructureFromCrystalStructureSuperCell(O, &Otemp);
  StreamPrintCrystalStructure(Otemp, stream);
  // need to free Otemp

  fprintf(stream, "\nSuper cell created from\n");

  StreamPrintCrystalStructure(O.Crystal, stream);
  fprintf(stream, "Number of Unit Cells in super cell: %4i\n", O.NCell);
  fprintf(stream, "Super Cell in Cell units:\n");
  fprintf(stream, "a: %4i %4i %4i\n", O.rSuperCell[0][0], O.rSuperCell[0][1], O.rSuperCell[0][2]);
  fprintf(stream, "b: %4i %4i %4i\n", O.rSuperCell[1][0], O.rSuperCell[1][1], O.rSuperCell[1][2]);
  fprintf(stream, "c: %4i %4i %4i\n", O.rSuperCell[2][0], O.rSuperCell[2][1], O.rSuperCell[2][2]);
  fprintf(stream, "Position of Unit Cells:\n");
  for (i = 0; i < O.NCell; i++)
  {
    fprintf(stream, "%4i | %4i %4i %4i\n", i, O.CellPos[i][0], O.CellPos[i][1], O.CellPos[i][2]);
  }
  fprintf(stream, "Super Cell Reciprocal lattice vectors in reciprocal Cell units:\n");
  fprintf(stream, "a: %12.7F %12.7F %12.7F\n", O.gSuperCell[0][0], O.gSuperCell[0][1], O.gSuperCell[0][2]);
  fprintf(stream, "b: %12.7F %12.7F %12.7F\n", O.gSuperCell[1][0], O.gSuperCell[1][1], O.gSuperCell[1][2]);
  fprintf(stream, "c: %12.7F %12.7F %12.7F\n", O.gSuperCell[2][0], O.gSuperCell[2][1], O.gSuperCell[2][2]);
  fprintf(stream, "k-vectors included in super Cell:\n");
  for (i = 0; i < O.NCell; i++)
  {
    fprintf(stream, "%4i | %12.7F %12.7F %12.7F\n", i, O.kVector[i][0], O.kVector[i][1], O.kVector[i][2]);
  }
  return 0;
}







/*  Prints a matrix
 */
uint32_t StreamPrintMatrix(MatrixType M, FILE * stream)
{
  uint32_t i, j;
  M.Name[MAXStrLength - 1] = '\0';
  fprintf(stream, "%s =\n", M.Name);
  fprintf(stream, "Real Part of Matrix with dimensions [Ni=%i][Nj=%i] ([Rows][Collums])\n", M.Ni, M.Nj);
  fprintf(stream, "        ");
  for (j = 0; j < M.Nj; j++)
  {
    fprintf(stream, "  [%12i]", j);
  }
  fprintf(stream, "\n");
  for (i = 0; i < M.Ni; i++)
  {
    fprintf(stream, "[%6i] ", i);
    for (j = 0; j < M.Nj; j++)
    {
      fprintf(stream, "%15.8E ", M.VRe[i][j]);
    }
    fprintf(stream, "\n");
  }
  if (M.QComplex)
  {
    fprintf(stream, "Imaginary Part\n");
    fprintf(stream, "        ");
    for (j = 0; j < M.Nj; j++)
    {
      fprintf(stream, "  [%12i]", j);
    }
    fprintf(stream, "\n");
    for (i = 0; i < M.Ni; i++)
    {
      fprintf(stream, "[%6i] ", i);
      for (j = 0; j < M.Nj; j++)
      {
        fprintf(stream, "%15.8E ", M.VIm[i][j]);
      }
      fprintf(stream, "\n");
    }
  }
  fprintf(stream, "\n");
  return 0;
}


/*  Prints a matrix
 */
uint32_t StreamPrintPlainMatrix(MatrixType M, FILE * stream)
{
  uint32_t i, j;
  if (M.QComplex)
  {
    for (i = 0; i < M.Ni; i++)
    {
      for (j = 0; j < M.Nj; j++)
      {
        fprintf(stream, "%19.12E %19.12E ", M.VRe[i][j], M.VIm[i][j]);
      }
      fprintf(stream, "\n");
    }
  }
  else
  {
    for (i = 0; i < M.Ni; i++)
    {
      for (j = 0; j < M.Nj; j++)
      {
        fprintf(stream, "%20.13E ", M.VRe[i][j]);
      }
      fprintf(stream, "\n");
    }
  }
  fprintf(stream, "\n");
  return 0;
}


uint32_t StreamPrintCompactMatrixMinimal(CompactMatrixType const M, FILE * stream)
{
  uint32_t i,j;
  if(M.QComplex)
  {
    for(i=0;i<M.Ni;i++)
    {
      for(j=0;j<M.Nj;j++)
      {
        fprintf(stream,"( %17.10E ",M.V[i*M.Nj+j].Re);
        if(M.V[i*M.Nj+j].Im>=0)
        {
          fprintf(stream,"+");
        }
        else
        {
          fprintf(stream,"-");
        }
        fprintf(stream,"%17.10E I) ",fabs(M.V[i*M.Nj+j].Im));
      }
      fprintf(stream,"\n");
    }
    fprintf(stream,"\n");
  }
  else
  {
    for(i=0;i<M.Ni;i++)
    {
      for(j=0;j<M.Nj;j++)
      {
        fprintf(stream,"%17.10E ",M.VRe[i*M.Nj+j]);
      }
      fprintf(stream,"\n");
    }
    fprintf(stream,"\n");
  }
  return 0;
}


/*  Prints the structure TightBindingDef
 *
 */
uint32_t StreamPrintTightBindingDef(TightBindingDef O, FILE * stream)
{
  int    i;
  O.Name[MAXStrLength - 1] = '\0';
  fprintf(stream, "\nSettings of a tight binding model: %s\n", O.Name);
  StreamPrintCrystalStructure(O.Cell, stream);
  fprintf(stream,"Containing a total number of %i orbitals\n",O.NFermion);
  fprintf(stream, "Hopping definitions ( %i )\n", O.nH);
  for (i = 0; i < O.nH; i++)
  {
    fprintf(stream, "Hopping from %i : %s - %s to %i : %s - %s with translation vector in unit cells: { %i , %i , %i } ({%15.8f %15.8f %15.8f })\n", O.H[i].Atompos1, O.Cell.Atom[O.H[i].Atompos1].Name,O.Cell.Atom[O.H[i].Atompos1].ShellName[O.H[i].Shellpos1], O.H[i].Atompos2, O.Cell.Atom[O.H[i].Atompos2].Name, O.Cell.Atom[O.H[i].Atompos2].ShellName[O.H[i].Shellpos2], O.H[i].TransVecs[0], O.H[i].TransVecs[1], O.H[i].TransVecs[2], O.H[i].ConnectingVector[0], O.H[i].ConnectingVector[1], O.H[i].ConnectingVector[2]);
    StreamPrintMatrix(O.H[i].Hop, stream);
  }
  fprintf(stream, "\n");
  //fprintf(stream, "Fermion index is given as:\n");
  //for(i=0;i<O.Cell.NAtoms;i++)
  //{
  //  for(j=0;j<O.Cell.Atom[i].Nshell;j++)
  //  {
  //    fprintf(stream, "  Hamiltonian for atom %i and shell %i starts at index %i\n",i,j,O.FermionIndex[i][j]);
  //  }
  //}
  //fprintf(stream, "\n");
  return 0;
}



/*  Prints the settings for a Dimerized spin chain model
 *
 */
uint32_t StreamPrintDimerizedSpinChainSettings(DimerizedSpinChainSettings O, FILE * stream)
{
  fprintf(stream, "\nSettings of Heisenberg Model\n");
  fprintf(stream, "Nearest Neighbor exchange constant JA:%22.15F\n", O.JA);
  fprintf(stream, "Nearest Neighbor exchange constant JB:%22.15F\n", O.JB);
  fprintf(stream, "Number of dimers in the chain        :    %i\n", O.Ndimer);
  fprintf(stream, "Additional information for numerical settings\n");
  fprintf(stream, "NPsi                 : %i\n", O.NPsi);
  fprintf(stream, "NE                   : %i\n", O.NE);
  fprintf(stream, "Emin                 : %24.17E\n", O.Emin);
  fprintf(stream, "Emax                 : %24.17E\n", O.Emax);
  fprintf(stream, "Gamma                : %24.17E\n", O.Gamma);
  fprintf(stream, "zeroLancEnergy       : %24.17E\n", O.zeroLancEnergy);
  fprintf(stream, "zeroLancWaveFunction : %24.17E\n", O.zeroLancWavefunction);
  fprintf(stream, "NTriDiag             : %i\n", O.NTriDiag);
  fprintf(stream, "NBitsKey             : %i\n", O.NBitsKey);
  fprintf(stream, "SpecFile             : %s\n", O.SpecFile);
  return 0;
}


/*  Prints the settings for a Heisenberg model
 *
 */
uint32_t StreamPrintHeisenbergModelSettings(HeisenbergModelSettings O, FILE * stream)
{
  fprintf(stream, "\nSettings of Heisenberg Model\n");
  fprintf(stream, "Nearest Neighbor exchange constant J:%22.15F\n", O.J);
  StreamPrintCrystalStructureSuperCell(O.Cluster, stream);
  fprintf(stream, "Additional information for numerical settings\n");
  fprintf(stream, "NPsi                 : %i\n", O.NPsi);
  fprintf(stream, "NE                   : %i\n", O.NE);
  fprintf(stream, "Emin                 : %24.17E\n", O.Emin);
  fprintf(stream, "Emax                 : %24.17E\n", O.Emax);
  fprintf(stream, "Gamma                : %24.17E\n", O.Gamma);
  fprintf(stream, "zeroLancEnergy       : %24.17E\n", O.zeroLancEnergy);
  fprintf(stream, "zeroLancWaveFunction : %24.17E\n", O.zeroLancWavefunction);
  fprintf(stream, "NTriDiag             : %i\n", O.NTriDiag);
  fprintf(stream, "NBitsKey             : %i\n", O.NBitsKey);
  fprintf(stream, "LocalFourierBasis    : %i\n", O.LocalFourierBasis);
  fprintf(stream, "SpecFile             : %s\n", O.SpecFile);
  return 0;
}


/*  Prints the settings for a Heisenberg model
 *
 */
uint32_t StreamPrinttJModelSettings(tJModelSettings O, FILE * stream)
{
  fprintf(stream, "\nSettings of tJ Model\n");
  fprintf(stream, "Nearest Neighbor exchange constant J:%22.15F\n", O.J);
  fprintf(stream, "Nearest Neighbor hopping constant  t:%22.15F\n", O.t);
  fprintf(stream, "Number of electrons                 : %i\n", O.nElectron);
  StreamPrintCrystalStructureSuperCell(O.Cluster, stream);
  fprintf(stream, "Additional information for numerical settings\n");
  fprintf(stream, "NPsi                 : %i\n", O.NPsi);
  fprintf(stream, "NE                   : %i\n", O.NE);
  fprintf(stream, "Emin                 : %24.17E\n", O.Emin);
  fprintf(stream, "Emax                 : %24.17E\n", O.Emax);
  fprintf(stream, "Gamma                : %24.17E\n", O.Gamma);
  fprintf(stream, "zeroLancEnergy       : %24.17E\n", O.zeroLancEnergy);
  fprintf(stream, "zeroLancWaveFunction : %24.17E\n", O.zeroLancWavefunction);
  fprintf(stream, "NTriDiag             : %i\n", O.NTriDiag);
  fprintf(stream, "NBitsKey             : %i\n", O.NBitsKey);
  fprintf(stream, "LocalFourierBasis    : %i\n", O.LocalFourierBasis);
  fprintf(stream, "SpecFile             : %s\n", O.SpecFile);
  return 0;
}


/*  Prints the settings for a Spin Lader model
 *
 */
uint32_t StreamPrintSpinLaderSettings(SpinLaderSettings O, FILE * stream)
{
  fprintf(stream, "\nSettings of Spin Lader Model\n");
  fprintf(stream, "Nearest Neighbor exchange constant Jrung:%22.15F\n", O.Jrung);
  fprintf(stream, "Nearest Neighbor exchange constant Jleg: %22.15F\n", O.Jleg);
  fprintf(stream, "Nearest Neighbor exchange constant Jcyc: %22.15F\n", O.Jcyc);
  fprintf(stream, "Number of rungs                        : %i\n", O.NSite);

  fprintf(stream, "Additional information for numerical settings\n");
  fprintf(stream, "NPsi                 : %i\n", O.NPsi);
  fprintf(stream, "NE                   : %i\n", O.NE);
  fprintf(stream, "Emin                 : %24.17E\n", O.Emin);
  fprintf(stream, "Emax                 : %24.17E\n", O.Emax);
  fprintf(stream, "Gamma                : %24.17E\n", O.Gamma);
  fprintf(stream, "zeroLancEnergy       : %24.17E\n", O.zeroLancEnergy);
  fprintf(stream, "zeroLancWaveFunction : %24.17E\n", O.zeroLancWavefunction);
  fprintf(stream, "NTriDiag             : %i\n", O.NTriDiag);
  fprintf(stream, "NBitsKey             : %i\n", O.NBitsKey);
  fprintf(stream, "SpecFile             : %s\n", O.SpecFile);
  return 0;
}


/*  Prints the settings for a Spin Lader model
 *
 */
uint32_t StreamPrintHubbardModelSettings(HubbardModelSettings O, FILE * stream)
{
  uint32_t iatom, ishell, ifermion;
  fprintf(stream, "\nSettings of Hubbard Model\n");
  fprintf(stream, "Crystal structure and super cell settings\n");
  StreamPrintCrystalStructureSuperCell(O.SuperCell, stream);
  fprintf(stream, "Tigth Binding definitions within unit cell\n");
  StreamPrintTightBindingDef(O.TightBinding, stream);
  fprintf(stream, "Tigth Binding definitions within super cell\n");
  StreamPrintTightBindingDef(O.TightBindingSuperCell, stream);
  fprintf(stream, "Additional Hamiltonian within unit cell\n");
  StreamPrintOperatorReadable(O.U, stream);
  fprintf(stream, "\nAdditional Hamiltonian within Super cell\n");
  StreamPrintOperatorReadable(O.USuperCell, stream);
  fprintf(stream, "\n Filling before diagonalization \n");
  for (iatom = 0; iatom < O.SuperCell.Crystal.NAtoms; iatom++)
  {
    fprintf(stream, "atom %2i ", iatom);
    for (ishell = 0; ishell < O.SuperCell.Crystal.Atom[iatom].Nshell; ishell++)
    {
      fprintf(stream, "shell %2i |", ishell);
      for (ifermion = 0; ifermion < O.SuperCell.Crystal.Atom[iatom].shelldegeneracy[ishell]; ifermion++)
      {
        fprintf(stream, "%2i ", O.NElectron[iatom][ishell][ifermion]);
      }
      fprintf(stream, "\n");
    }
    fprintf(stream, "\n");
  }
  fprintf(stream, "Additional information for numerical settings\n");
  fprintf(stream, "Nq                   : %i %i %i\n", O.Nq[0], O.Nq[1], O.Nq[2]);
  fprintf(stream, "NPsi (max)           : %i\n", O.NPsi);
  fprintf(stream, "DeltaEMax            : %lf\n", O.DeltaEMax);
  fprintf(stream, "NE                   : %i\n", O.NE);
  fprintf(stream, "Emin                 : %24.17E\n", O.Emin);
  fprintf(stream, "Emax                 : %24.17E\n", O.Emax);
  fprintf(stream, "Gamma                : %24.17E\n", O.Gamma);
  fprintf(stream, "zeroLancEnergy       : %24.17E\n", O.zeroLancEnergy);
  fprintf(stream, "zeroLancWaveFunction : %24.17E\n", O.zeroLancWavefunction);
  fprintf(stream, "NTriDiag             : %i\n", O.NTriDiag);
  fprintf(stream, "NBitsKey             : %i\n", O.NBitsKey);
  fprintf(stream, "SpecFile             : %s\n", O.SpecFile);
  return 0;
}


/*  Prints the settings for a XAS spectrum of a single site DMFT model
 *
 */
uint32_t StreamPrintDMFTXASSettings(DMFTXASSettings O, FILE * stream)
{
  fprintf(stream, "\nSettings of XAS on DMFT model\n");
  fprintf(stream, "BathFilename         : %s\n", O.BathFilename);
  fprintf(stream, "U                    : %24.17E\n", O.U);
  fprintf(stream, "Q                    : %24.17E\n", O.Q);
  fprintf(stream, "NBathMax             : %i\n", O.NBathMax);
  fprintf(stream, "NE                   : %i\n", O.NE);
  fprintf(stream, "Emin                 : %24.17E\n", O.Emin);
  fprintf(stream, "Emax                 : %24.17E\n", O.Emax);
  fprintf(stream, "Gamma                : %24.17E\n", O.Gamma);
  fprintf(stream, "zeroLancEnergy       : %24.17E\n", O.zeroLancEnergy);
  fprintf(stream, "zeroLancWaveFunction : %24.17E\n", O.zeroLancWavefunction);
  fprintf(stream, "NTriDiag             : %i\n", O.NTriDiag);
  fprintf(stream, "NBitsKey             : %i\n", O.NBitsKey);
  fprintf(stream, "SpecFile             : %s\n", O.SpecFile);
  return 0;
}

/*
 *
 */
uint32_t StreamPrintDMFTSettings(DMFTSettings O, FILE * stream)
{
  uint32_t i;
  fprintf(stream, "\nSettings of self-consistent DMFT calculation\n");
  fprintf(stream, "G0FileName               : %s\n", O.G0FileName);
  fprintf(stream, "RootFileName             : %s\n", O.RootFileName);
  fprintf(stream, "SpecRootFileName         : %s\n", O.SpecRootFileName);
  fprintf(stream, "WaveFunctionRootFileName : %s\n", O.WaveFunctionRootFileName);
  fprintf(stream, "number of electrons      : %22.15E\n", O.Nelectron);
  fprintf(stream, "Chemical potential of G0 : %22.15E\n", O.mu0);
  fprintf(stream, "Chemical potential of G  : %22.15E\n", O.mu);
  fprintf(stream, "Effective temperature    : %22.15E\n", O.Teff);
  fprintf(stream,
          "                         converged results do not depent on Teff, convergence speed and stability do\n");
  fprintf(stream, "Options for Chemical pot : %i\n", O.NmuSym);
  if ((O.NmuSym) & 1)
    fprintf(stream, "                         ( 1) mu and mu0 are read from G0 and Sigma, values in imput ignored;\n");
  if ((O.NmuSym) & 2)
    fprintf(stream, "                         ( 2) Calculations done at fixed mu, Nelectron ignored;\n");
  fprintf(stream, "# of spin-orbitals-sites : %i\n", O.NSite);
  fprintf(stream, "iStart                   : %i\n", O.iStart);
  fprintf(stream, "Ni                       : %i\n", O.Ni);
  fprintf(stream, "Dyson perturbation       : %i\n", O.Dyson);
  fprintf(stream, "Symmetrizaton options    : %i\n", O.SymOpt);
  if ((O.SymOpt) & 1)
    fprintf(stream, "                         ( 1) HAnd_grd += (U+1)S^2;\n");
  if ((O.SymOpt) & 2)
    fprintf(stream, "                         ( 2) G(w) +=-G(-w)^*; G(w)/=2;\n");
  if ((O.SymOpt) & 4)
    fprintf(stream, "                         ( 4) psi_cluster lives on 2 sites per orbital\n");
  if ((O.SymOpt) & 8)
    fprintf(stream, "                         ( 8) mu = U/2\n");
  fprintf(stream, "NGCluster                : %i\n", O.NGCluster);
  fprintf(stream, "NGBath                   : %i\n", O.NGBath);
  fprintf(stream, "Mix                      : %24.17E\n", O.Mix);
  fprintf(stream, "U                        : %24.17E\n", O.U);
  fprintf(stream, "NBath                    : %i\n", O.NBath);
  fprintf(stream, "Spectra                  : %i\n", O.Spectra);
  if ((O.Spectra) & 1)
    fprintf(stream, "                        (   1) Calculate spectra each loop;\n");
  else
    fprintf(stream, "                        (   0) Calculate spectra in the end;\n");
  if ((O.Spectra) & 2)
    fprintf(stream, "                        (   2) Sz;\n");
  if ((O.Spectra) & 4)
    fprintf(stream, "                        (   4) N;\n");
  if ((O.Spectra) & 8)
    fprintf(stream, "                        (   8) N-1;\n");
  if ((O.Spectra) & 16)
    fprintf(stream, "                        (  16) N_up N_dn;\n");
  if ((O.Spectra) & 32)
    fprintf(stream, "                        (  32) XAS;\n");
  if ((O.Spectra) & 64)
    fprintf(stream, "                        (  64) core PES;\n");
  if ((O.Spectra) & 128)
    fprintf(stream, "                        ( 128) XAS - d0;\n");
  if ((O.Spectra) & 256)
    fprintf(stream, "                        ( 256) XAS - d1;\n");
  if ((O.Spectra) & 512)
    fprintf(stream, "                        ( 512) cPES - d0;\n");
  if ((O.Spectra) & 1024)
    fprintf(stream, "                        (1024) cPES - d1;\n");
  if ((O.Spectra) & 2048)
    fprintf(stream, "                        (2048) cPES - d2;\n");
  if ((O.Spectra) & 4064)
  {
    fprintf(stream, "  NQ                     : %i\n", O.NQ);
    for (i = 0; i < O.NQ; i++)
    {
      fprintf(stream, "  Q[%3i]                 : %24.17E\n", i, O.Q[i]);
    }
  }
  fprintf(stream, "NE                       : %i\n", O.NE);
  fprintf(stream, "Emin                     : %24.17E\n", O.Emin);
  fprintf(stream, "Emax                     : %24.17E\n", O.Emax);
  fprintf(stream, "Gamma                    : %24.17E\n", O.Gamma);
  fprintf(stream, "zeroLancEnergy           : %24.17E\n", O.zeroLancEnergy);
  fprintf(stream, "zeroLancWaveFunction     : %24.17E\n", O.zeroLancWavefunction);
  fprintf(stream, "NTriDiag                 : %i\n", O.NTriDiag);
  fprintf(stream, "NBitsKey                 : %i\n", O.NBitsKey);
  return 0;
}













/*  Prints a matrix
 */
uint32_t StreamPrintSparseMatrix(SparseMatrixType M, FILE * stream)
{
  uint32_t i;
  M.Name[MAXStrLength - 1] = '\0';
  fprintf(stream, "%s =\n", M.Name);
  if (M.QComplex)
  {
    fprintf(stream, "   [i]    [j]          Re(M[i][j])          Im(M[i][j])\n");
    for (i = 0; i < M.N; i++)
    {
      fprintf(stream, "%6i %6i %22.15E %22.15E\n", M.i[i], M.j[i], M.VRe[i], M.VRe[i]);
    }
  }
  else
  {
    fprintf(stream, "   [i]    [j]          Re(M[i][j])\n");
    for (i = 0; i < M.N; i++)
    {
      fprintf(stream, "%6i %6i %22.15E\n", M.i[i], M.j[i], M.VRe[i]);
    }
  }
  fprintf(stream, "\n");
  return 0;
}


uint32_t StreamPrintTriDiagonalMatrix(TriDiagonalMatrixType M, FILE * stream)
{
  uint32_t i;
  M.Name[MAXStrLength - 1] = '\0';
  fprintf(stream, "%s =\n", M.Name);
  fprintf(stream, "TriDiagonalMatrix element a[i-1]=M[i][i] element b[i-1]=M[i,i-1]=M[i-1,i] others are zero.\n");
  fprintf(stream, "Zero elements are prefactors for continued fraction representation\n");
  fprintf(stream, "i    a[i]                   b[i]\n");
  for (i = 0; i < M.N; i++)
  {
    fprintf(stream, "%4i %22.15E %22.15E\n", i, M.a[i], M.b[i]);
  }
  fprintf(stream, "%4i %22.15E                       \n", M.N, M.a[M.N]);
  return 0;
}

uint32_t StreamPrintBlokTriDiagonalMatrix(BlokTriDiagonalMatrixType M, FILE * stream)
{
  uint32_t i;
  M.Name[MAXStrLength - 1] = '\0';
  fprintf(stream, "%s =\n", M.Name);
  fprintf(stream, "BlokTriDiagonalMatrix element a[i-1]=M[i][i] element b[i-1]=M[i,i+1]=M[i+1,i] others are zero.\n");
  fprintf(stream, "Zero elements are prefactors for continued fraction representation\n");
  for (i = 0; i <= M.N; i++)
  {
    fprintf(stream, "a[%i]\n", i);
    StreamPrintMatrix(M.a[i], stream);
  }
  for (i = 0; i < M.N; i++)
  {
    fprintf(stream, "b[%i]\n", i);
    StreamPrintMatrix(M.b[i], stream);
  }
  return 0;
}


uint32_t StreamPrintBlockTriDiagonalMatrix(BlockTriDiagonalMatrixType M, FILE * stream)
{
  uint32_t i;
  CompactMatrixType alpha;
  alpha.Ni=M.BlockSize;
  alpha.Nj=M.BlockSize;
  alpha.QComplex=M.QComplex;
  M.Name[MAXStrLength - 1] = '\0';
  fprintf(stream, "%s :\n", M.Name);
  fprintf(stream, "BlockTriDiagonalMatrix element a[i-1]=M[i][i] element b[i-1]=M[i,i+1]=M[i+1,i] others are zero.\n");
  fprintf(stream, "Zero elements are prefactors for continued fraction representation\n");
  fprintf(stream, "BlockSize %i number of Blocks %i\n",M.BlockSize, M.N);
  for (i = 0; i <= M.N; i++)
  {
    if(M.QComplex)
    {
      alpha.V=&(M.a[i*M.BlockSize*M.BlockSize]);
    }
    else
    {
      alpha.VRe=&(M.aRe[i*M.BlockSize*M.BlockSize]);
    }
    fprintf(stream, "a[%i]\n", i);
    StreamPrintCompactMatrixMinimal(alpha, stream);
  }
  for (i = 0; i < M.N; i++)
  {
    if(M.QComplex)
    {
      alpha.V=&(M.b[i*M.BlockSize*M.BlockSize]);
    }
    else
    {
      alpha.VRe=&(M.bRe[i*M.BlockSize*M.BlockSize]);
    }
    fprintf(stream, "b[%i]\n", i);
    StreamPrintCompactMatrixMinimal(alpha, stream);
  }
  return 0;
}

uint32_t StreamPrintDoubleTriDiagonalMatrix(DoubleTriDiagonalMatrixType M, FILE * stream)
{
  M.Name[MAXStrLength - 1] = '\0';
  fprintf(stream, "%s =\n", M.Name);
  fprintf(stream, "DoubleTriDiagonalMatrix mu = %22.15E Epsilon = %22.15E\n", M.mu, M.Epsilon);
  fprintf(stream, "                        a0 = %22.15E      b0 = %22.15E\n", M.a0, M.b0);
  fprintf(stream, "Valence part (E<Ef):\n");
  StreamPrintTriDiagonalMatrix(M.val, stream);
  fprintf(stream, "Conduction part (E>Ef):\n");
  StreamPrintTriDiagonalMatrix(M.con, stream);
  return 0;
}



uint32_t StreamPrintBADoubleTriDiagonalMatrix(BADoubleTriDiagonalMatrixType M, FILE * stream)
{
  M.Name[MAXStrLength - 1] = '\0';
  fprintf(stream, "%s =\n", M.Name);
  fprintf(stream, "BADoubleTriDiagonalMatrix mu = %22.15E Epsilon = %22.15E\n", M.mu, M.Epsilon);
  fprintf(stream, "                          a0 = %22.15E      b0 = %22.15E\n", M.a0, M.b0);
  fprintf(stream, "Valence part (E<Ef):\n");
  StreamPrintTriDiagonalMatrix(M.val, stream);
  fprintf(stream, "Conduction part (E>Ef):\n");
  StreamPrintTriDiagonalMatrix(M.con, stream);
  fprintf(stream, "Transformation matrix T:\n");
  fprintf(stream, "%22.15E %22.15E\n", M.val.b[0], M.con.b[0]);
  fprintf(stream, "%22.15E %22.15E\n", M.con.b[0], -M.val.b[0]);
  fprintf(stream, "Occupation: val.b[0]^2 %22.15E con.b[0]^2 %22.15E Sum %22.15E b0^2 %22.15E\n", Sqr(M.val.b[0]),
          Sqr(M.con.b[0]), Sqr(M.val.b[0]) + Sqr(M.con.b[0]), Sqr(M.b0));
  fprintf(stream, "Transformed values\n");
  fprintf(stream, "acvB  :%22.15E  (=val.b[0]^2 * val.a[1] + con.b[0]^2 * con.a[1] %22.15E)\n", M.acvB,
          Sqr(M.val.b[0]) * M.val.a[1] + Sqr(M.con.b[0]) * M.con.a[1]);
  fprintf(stream, "acvA  :%22.15E  (=val.b[0]^2 * con.a[1] + con.b[0]^2 * val.a[1] %22.15E)\n", M.acvA,
          Sqr(M.val.b[0]) * M.con.a[1] + Sqr(M.con.b[0]) * M.val.a[1]);
  fprintf(stream, "acvBA :%22.15E  (=val.b[0] * con.b[0] * (val.a[1]-con.a[1])     %22.15E)\n", M.acvBA,
          M.val.b[0] * M.con.b[0] * (M.val.a[1] - M.con.a[1]));
  fprintf(stream, "bvB0  :%22.15E  (=val.b[0] * val.b[1]                           %22.15E)\n", M.bvB0,
          M.val.b[0] * M.val.b[1]);
  fprintf(stream, "bcB0  :%22.15E  (=con.b[0] * con.b[1]                           %22.15E)\n", M.bcB0,
          M.con.b[0] * M.con.b[1]);
  fprintf(stream, "bvA0  :%22.15E  (=con.b[0] * val.b[1]                           %22.15E)\n", M.bvA0,
          M.con.b[0] * M.val.b[1]);
  fprintf(stream, "bcA0  :%22.15E  (=-val.b[0] * con.b[1]                          %22.15E)\n", M.bcA0,
          -M.val.b[0] * M.con.b[1]);
  return 0;
}

uint32_t StreamPrintSecondDoubleTriDiagonalMatrix(SecondDoubleTriDiagonalMatrixType M, FILE * stream)
{
  M.Name[MAXStrLength - 1] = '\0';
  fprintf(stream, "%s =\n", M.Name);
  fprintf(stream, "SecondDoubleTriDiagonalMatrix Ef    = %22.15E Epsilon    = %22.15E\n", M.mu, M.Epsilon);
  fprintf(stream, "                              a0    = %22.15E      b0    = %22.15E\n", M.a0, M.b0);
  fprintf(stream, "                              Dt.a0 = %22.15E      Dt.b0 = %22.15E\n", M.Dt.a0, M.Dt.b0);
  fprintf(stream, "Valence part (E<Ef):\n");
  StreamPrintTriDiagonalMatrix(M.Dt.val, stream);
  fprintf(stream, "Conduction part (E>Ef):\n");
  StreamPrintTriDiagonalMatrix(M.Dt.con, stream);
  return 0;
}

uint32_t StreamPrintBlokDoubleTriDiagonalMatrix(BlokDoubleTriDiagonalMatrixType M, FILE * stream)
{
  M.Name[MAXStrLength - 1] = '\0';
  fprintf(stream, "%s =\n", M.Name);
  fprintf(stream, "BlokDoubleTriDiagonalMatrix Ef = %22.15E Epsilon = %22.15E\n", M.mu, M.Epsilon);
  fprintf(stream, "a[0]\n");
  StreamPrintMatrix(M.a0, stream);
  fprintf(stream, "b[0]\n");
  StreamPrintMatrix(M.b0, stream);
  fprintf(stream, "Valence part (E<Ef):\n");
  StreamPrintBlokTriDiagonalMatrix(M.val, stream);
  fprintf(stream, "Conduction part (E>Ef):\n");
  StreamPrintBlokTriDiagonalMatrix(M.con, stream);
  return 0;
}

uint32_t StreamPrintAndersonMatrix(AndersonMatrixType M, FILE * stream)
{
  uint32_t i;
  M.Name[MAXStrLength - 1] = '\0';
  fprintf(stream, "%s =\n", M.Name);
  fprintf(stream, "AndersonMatrix element a[i-1]=M[i][i] element b[i-1]=M[i,i-1]=M[i-1,i] others are zero.\n");
  fprintf(stream, "Zero elements are prefactors for continued fraction representation\n");
  fprintf(stream, "Chemical potential is %22.15E\n", M.mu);
  fprintf(stream, "i    a[i]                   b[i]\n");
  for (i = 0; i < M.N; i++)
  {
    fprintf(stream, "%4i %22.15E %22.15E\n", i, M.a[i], M.b[i]);
  }
  fprintf(stream, "%4i %22.15E                       \n", M.N, M.a[M.N]);
  return 0;
}

uint32_t StreamPrintBlokAndersonMatrix(BlokAndersonMatrixType M, FILE * stream)
{
  uint32_t i;
  M.Name[MAXStrLength - 1] = '\0';
  fprintf(stream, "%s =\n", M.Name);
  fprintf(stream, "BlokAndersonMatrix element a[i]=M[i][i] element b[i]=M[0,i+1]=M[i+1,0] others are zero.\n");
  fprintf(stream, "Zero elements are prefactors for continued fraction representation\n");
  for (i = 0; i <= M.N; i++)
  {
    fprintf(stream, "a[%i]\n", i);
    StreamPrintMatrix(M.a[i], stream);
  }
  for (i = 0; i < M.N; i++)
  {
    fprintf(stream, "b[%i]\n", i);
    StreamPrintMatrix(M.b[i], stream);
  }
  return 0;
}

uint32_t StreamPrintBlockAndersonMatrix(BlockAndersonMatrixType M, FILE * stream)
{
  uint32_t i;
  CompactMatrixType alpha;
  alpha.Ni=M.BlockSize;
  alpha.Nj=M.BlockSize;
  alpha.QComplex=M.QComplex;
  M.Name[MAXStrLength - 1] = '\0';
  fprintf(stream, "%s =\n", M.Name);
  fprintf(stream, "BlockAndersonMatrix element a[i]=M[i][i] element b[i]=M[0,i+1]=M[i+1,0] others are zero.\n");
  fprintf(stream, "Zero elements are prefactors for continued fraction representation\n");
  for (i = 0; i <= M.N; i++)
  {
    if(M.QComplex)
    {
      alpha.V=&(M.a[i*M.BlockSize*M.BlockSize]);
    }
    else
    {
      alpha.VRe=&(M.aRe[i*M.BlockSize*M.BlockSize]);
    }
    fprintf(stream, "a[%i]\n", i);
    StreamPrintCompactMatrixMinimal(alpha, stream);
  }
  for (i = 0; i < M.N; i++)
  {
    if(M.QComplex)
    {
      alpha.V=&(M.b[i*M.BlockSize*M.BlockSize]);
    }
    else
    {
      alpha.VRe=&(M.bRe[i*M.BlockSize*M.BlockSize]);
    }
    fprintf(stream, "b[%i]\n", i);
    StreamPrintCompactMatrixMinimal(alpha, stream);
  }
  return 0;
}

uint32_t StreamPrintListOfPoles(ListOfPoles M, FILE * stream)
{
  uint32_t i;
  double sum = 0;
  for (i = 0; i < M.N; i++)
    sum += M.b[i];
  M.Name[MAXStrLength - 1] = '\0';
  fprintf(stream, "%s =\n", M.Name);
  fprintf(stream, "List of Poles.\n");
  fprintf(stream, "Chemical potential is %22.15E\n", M.mu);
  fprintf(stream, "a0 = %22.15E, sum_i b[i] = %22.15E\n", M.a0, sum);
  fprintf(stream, "i    a[i]                   b[i]\n");
  for (i = 0; i < M.N; i++)
  {
    fprintf(stream, "%4i %22.15E %22.15E\n", i, M.a[i], M.b[i]);
  }
  return 0;
}

uint32_t StreamPrintListOfPolesEssentials(ListOfPoles M, FILE * stream)
{
  uint32_t i;
  double sum = 0, wsum = 0, wwsum = 0, wwwsum = 0;
  double Osum = 0, Owsum = 0, Owwsum = 0, Owwwsum = 0;
  for (i = 0; i < M.N; i++)
  {
    sum += M.b[i];
    wsum += M.b[i] * M.a[i];
    wwsum += M.b[i] * M.a[i] * M.a[i];
    wwwsum += M.b[i] * M.a[i] * M.a[i] * M.a[i];
    if (M.a[i] < M.mu)
    {
      Osum += M.b[i];
      Owsum += M.b[i] * M.a[i];
      Owwsum += M.b[i] * M.a[i] * M.a[i];
      Owwwsum += M.b[i] * M.a[i] * M.a[i] * M.a[i];
    }
  }
  M.Name[MAXStrLength - 1] = '\0';
  fprintf(stream, "%s : (Moments of List of poles of size (%i))\n", M.Name, M.N);
  fprintf(stream, "Chemical potential is %22.15E\n", M.mu);
  fprintf(stream, "a0  = %22.15E\n", M.a0);
  fprintf(stream, "w0  = %22.15E, w1  = %22.15E, w2  = %22.15E, w3  = %22.15E\n", sum, wsum, wwsum, wwwsum);
  fprintf(stream, "w0+ = %22.15E, w1+ = %22.15E, w2+ = %22.15E, w3+ = %22.15E\n", Osum, Owsum, Owwsum, Owwwsum);
  fprintf(stream, "  Energy-mu               Weight\n");
  for (i = 0; i < Min(M.N, 2); i++)
  {
    fprintf(stream, "%22.15E %22.15E\n", M.a[i] - M.mu, M.b[i]);
  }
  for (i = Max(2, M.N - 2); i < M.N; i++)
  {
    fprintf(stream, "%22.15E %22.15E\n", M.a[i] - M.mu, M.b[i]);
  }
  fprintf(stream, "\n");
  return 0;
}

uint32_t StreamPrintBlockListOfPoles(BlockListOfPoles const M, FILE * stream)
{
  uint32_t i;
  CompactMatrixType b;
  b.QComplex=M.QComplex;
  b.Ni=M.BlockSize;
  b.Nj=M.BlockSize;
  fprintf(stream, "block list of poles: %s\n",M.Name);
  fprintf(stream, "Chemical potential is %22.15E\n",M.mu);
  fprintf(stream, "a0:\n");
  if(M.QComplex)
  {
    b.V = M.a0;
  }
  else
  {
    b.VRe = M.a0Re;
  }
  StreamPrintCompactMatrixMinimal(b, stream);
  for(i=0;i<M.N;i++)
  {
    if(M.QComplex)
    {
      fprintf(stream,"block %i a_i=%22.15E B_i:\n",i,M.V[i*(M.BlockSize*M.BlockSize+1)].Re);
      b.V=&(M.V[i*(M.BlockSize*M.BlockSize+1)+1]);
    }
    else
    {
      fprintf(stream,"block %i a_i=%22.15E B_i:\n",i,M.VRe[i*(M.BlockSize*M.BlockSize+1)]);
      b.VRe=&(M.VRe[i*(M.BlockSize*M.BlockSize+1)+1]);
    }
    StreamPrintCompactMatrixMinimal(b, stream);
  }
  fprintf(stream,"\n");
  return 0;
}

uint32_t StreamPrintAndersonMatrixEssentials(AndersonMatrixType M, FILE * stream)
{
  uint32_t i;
  double tothyb = 0;
  M.Name[MAXStrLength - 1] = '\0';
  fprintf(stream, "%s : (Moments of Anderson Matrix of size (%i))\n", M.Name, M.N);
  fprintf(stream, "Chemical potential is %22.15E\n", M.mu);
  fprintf(stream, "a0  = %22.15E, a1  = %22.15E, b0  = %22.15E\n", M.a[0], M.a[1], M.b[0]);
  for (i = 1; i < M.N; i++)
  {
    tothyb += Sqr(M.b[i]);
  }
  fprintf(stream, "total hybridization (sum b^2) is: %22.15E\n", tothyb);
  fprintf(stream, "\n");
  return 0;
}


/*
 * Prints a matrix in plaim format suitable for rereading
 */
uint32_t StreamPrintPlainTriDiagonalMatrix(TriDiagonalMatrixType M, FILE * stream)
{
  uint32_t i;
  fprintf(stream, "tri %22.15E\n", M.mu);
  fprintf(stream, "%22.15E %22.15E %i\n", M.a[0], M.b[0], M.N);
  for (i = 1; i < M.N; i++)
  {
    fprintf(stream, "%22.15E %22.15E\n", M.a[i], M.b[i]);
  }
  fprintf(stream, "%22.15E\n", M.a[M.N]);
  return 0;
}

uint32_t StreamPrintPlainAndersonMatrix(AndersonMatrixType M, FILE * stream)
{
  uint32_t i;
  fprintf(stream, "and %22.15E\n", M.mu);
  fprintf(stream, "%22.15E %22.15E %i\n", M.a[0], M.b[0], M.N);
  for (i = 1; i < M.N; i++)
  {
    fprintf(stream, "%22.15E %22.15E\n", M.a[i], M.b[i]);
  }
  fprintf(stream, "%22.15E\n", M.a[M.N]);
  return 0;
}

uint32_t StreamPrintPlainListOfPoles(ListOfPoles M, FILE * stream)
{
  uint32_t i;
  fprintf(stream, "lp %22.15E\n", M.mu);
  fprintf(stream, "%22.15E %i\n", M.a0, M.N);
  for (i = 0; i < M.N; i++)
  {
    fprintf(stream, "%22.15E %22.15E\n", M.a[i], M.b[i]);
  }
  return 0;
}


uint32_t StreamPrintBlockNaturalImpurityOrbitalMatrix(BlockNaturalImpurityOrbitalMatrixType M, FILE * stream)
{
  fprintf(stream,"BlockNaturalImpurityOrbitalMatrix: %s\n",M.Name);
  fprintf(stream,"mu %22.15E\n",M.mu);
  fprintf(stream,"Valence part\n");
  StreamPrintBlockTriDiagonalMatrix(M.val, stream);
  fprintf(stream,"Conduction part\n");
  StreamPrintBlockTriDiagonalMatrix(M.con, stream);
  fprintf(stream, "first 4 blocks\n");
  PrintCompactMatrixMinimal(M.TMT);

  return 0;
}








uint32_t StreamPrintHarmonicsExpansion(HarmonicsExpansion M, FILE * stream)
{
  uint32_t i;
  fprintf(stream, "Expansion on spherical harmonics\n");
  if (M.QComplex)
  {
    printf("  i   l    m    Value\n");
    for (i = 0; i < M.N; i++)
    {
      fprintf(stream, "%3i %3i %4i %22.15E + %22.15E I\n", i, M.l[i], M.m[i], M.Re[i], M.Im[i]);
    }
  }
  else
  {
    printf("  i   l    m    Value\n");
    for (i = 0; i < M.N; i++)
    {
      fprintf(stream, "%3i %3i %4i %22.15E\n", i, M.l[i], M.m[i], M.Re[i]);
    }
  }
  fprintf(stream, "\n");
  return 0;
}

uint32_t StreamPrintBandstructurePath(BandstructurePathType M, FILE * stream)
{
  uint32_t i;
  fprintf(stream, "Band-structure path through k-space\n");
  fprintf(stream, "Label N-kpoints Qx          Qy          Qz\n");
  for(i=0;i<M.NQ;i++)
  {
    fprintf(stream, "%c     %4i     %11.8lf %11.8lf %11.8lf\n",M.Q[i].label,M.Q[i].Nq,M.Q[i].q[0],M.Q[i].q[1],M.Q[i].q[2]);
  }
	return 0;
}

string numberToStringMaxWidth(double num, int maxWidth)
{
	std::stringstream sstr;
	std::streamsize old_precision = sstr.precision();
	int power = findPower(num);
	int powerDigits = findPower(power) + 1;
	if (power < 0)
		++powerDigits;
	if (num < 0)
		--maxWidth;
	if (abs(num) < 1)
	{
		if (abs(num) > pow(10,-powerDigits-2) || num == 0)
		{
			if (maxWidth + power - 1 > 0)
				sstr.precision(maxWidth + power - 1);
			else
				sstr.precision(1);
			sstr << num;
		} else
		{
			if (maxWidth - powerDigits - 2 > 0)
				sstr.precision(maxWidth - powerDigits - 2);
			else
				sstr.precision(1);
			sstr << num * pow(10, -power) << "e" << power;
		}
	} else if (abs(num) > pow(10, maxWidth - 1))
	{
		if (maxWidth - powerDigits - 2 > 0)
			sstr.precision(maxWidth - powerDigits - 2);
		else
			sstr.precision(1);
		sstr << num * pow(10, -power) << "e" << power;
	} else
	{
		if (maxWidth - 1 > 0)
			sstr.precision(maxWidth - 1);
		else
			sstr.precision(1);
		sstr << num;
	}
	sstr.precision(old_precision);
	string out;
	sstr >> out;
	return out;
}


void numberOutMaxWidth(double num, int maxWidth)
{
	cout << numberToStringMaxWidth(num,maxWidth);
}

void fillWithBlank(int width, string::size_type occupiedSpace)
{
	if (occupiedSpace >= width)
	{
		cout << " ";
	} else
	{
		for (int i = 0; i < width - occupiedSpace; ++i)
			cout << " ";
	}
}

void fillWithBlank(int width)
{
	fillWithBlank(width,0);
}

void printExpectationValues(const vector < WaveFunctionType >  psiList, OperatorType  hamil, const vector < OperatorType >  opList, bool hasHamiltonian, bool sortByEnergy, bool printSigma, double zero, int printHeaderEvery, int colWidth)
{
	vector<double> EValues(psiList.size());
	vector<double> sigmaValues(psiList.size());
	vector<vector<double> > expValues(psiList.size());
	int hamilAddition = 0, sigmaAddition = 0;
	bool opHermetian = false, opRestricted = false;
	if (hasHamiltonian)
	{
		hamilAddition = 1;
		OperatorChop(&hamil, sqrt(zero));
		opHermetian = isOperatorAlmostHermitian(hamil, zero);
		opRestricted = (hamil.Res.N > 0);
		if (!opHermetian)
		{
			cout << "PrintExpectationValues not yet implemented for non-Hermitian operators, the reults will probably be rubbish." << endl;
		}
		if (opRestricted)
		{
			for (int i = 0; i < psiList.size(); ++i)
			{
				EValues[i] = PsiHermitianOperatorPsi(hamil, psiList[i]); //TODO: es gibt das noch nicht mit Restrictions.
			}
		} else
		{
			for (int i = 0; i < psiList.size(); ++i)
			{
				EValues[i] = PsiHermitianOperatorPsi(hamil, psiList[i]);
			}
		}
		if (printSigma)
		{
			sigmaAddition = 1;
//			TODO: make this possible for non-hermetian operators
			if (opRestricted)
			{
				for (int i = 0; i < psiList.size(); ++i)
					sigmaValues[i] = sqrt(abs(EValues[i] * EValues[i] - PsiHermitianOperatorSqrPsiRestricted(hamil, psiList[i], hamil.Res)));
			} else
			{
				for (int i = 0; i < psiList.size(); ++i)
					sigmaValues[i] = sqrt(abs(EValues[i] * EValues[i] - PsiHermitianOperatorSqrPsi(hamil, psiList[i])));
			}
		}
	}

	for (int i = 0; i < psiList.size(); ++i)
	{

		expValues[i] = vector<double>(opList.size() + hamilAddition + sigmaAddition);
		if (hasHamiltonian)
		{
			expValues[i][0] = EValues[i];
			if (printSigma)
				expValues[i][expValues[i].size() - 1] = sigmaValues[i];
		}
		for (int j = 0; j < opList.size(); ++j)
			expValues[i][j + hamilAddition] = PsiHermitianOperatorPsi(opList[j], psiList[i]); //TODO: If Abfrage auf Restrictions. Leider ändern Restrictions an dieser Stelle auch nichts.
	}

	if (hasHamiltonian && sortByEnergy)
		sort(expValues.begin(), expValues.end(), matrixSorterByFirstValue);

	for (int i = 0; i < psiList.size(); ++i)
	{
		if (i == 0 || (printHeaderEvery > 0 && i % printHeaderEvery == 0))
			printExpectationValuesHeader(opList, hasHamiltonian, printSigma, colWidth);
		int maxISize = findPower(psiList.size()) + 1;
		fillWithBlank(colWidth - 2, maxISize);
		cout.width(maxISize);
		cout << (i + 1);
		cout << "  ";
		for (int j = 0; j < opList.size() + hamilAddition; ++j)
		{
			double out = expValues[i][j];
			if (out <= zero)
				out = 0;
			string outStr = numberToStringMaxWidth(out, colWidth - 1);
			cout << outStr;
			fillWithBlank(colWidth, outStr.length());
		}
		if (hasHamiltonian && printSigma)
		{
			double out = expValues[i][expValues[i].size() - 1];
			string outStr = numberToStringMaxWidth(out, colWidth);
			cout << outStr;
		}
		cout << endl;
	}
}

void printExpectationValuesHeader(const vector < OperatorType > opList, bool hasHamiltonian, bool printSigma, int colWidth)
{
	cout << endl;
	fillWithBlank(colWidth);
	if (hasHamiltonian)
	{
		cout << "E";
		fillWithBlank(colWidth,1);
	}
	for (int i = 0; i < opList.size(); ++i)
	{
		string name(opList[i].Name);
		if (name.length() > colWidth-1)
		{
			name = name.substr(0,colWidth-1);
		}
		cout << name;
		fillWithBlank(colWidth, name.length());
	}
	if (hasHamiltonian && printSigma)
		cout << "std-dev";
	cout << endl;
}

