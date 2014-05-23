
/*
 *  Calculate the density matrix of a wave-function
 *  On input psi contains a wave-function and Rho points to a not allocated matrix
 *  On excite Rho contains an psi.NFermion + psi.NBoson by psi.NFermion + psi.NBoson matrix
 */
uint32_t CalculateDensityMatrix(MatrixType * const Rho, const WaveFunctionType psi)
{
  uint32_t i, j;
  OperatorType O;
  Rho->QComplex = psi.QComplex;
  Rho->Ni = psi.NFermion + psi.NBoson;
  Rho->Nj = psi.NFermion + psi.NBoson;
  if (InitMatrix(Rho))
  {
    printf("InitMatrix failed in CalculateDensityMatrix\n");
    fflush(stdout);
    return 1;
  }
  strcpy(O.Name, "Calc. Den. Mat. Int. Opp. ");
  strcpy(Rho->Name, "Density Matrix.");
  O.QComplex = 0;
  O.NFermion = psi.NFermion;
  O.NBoson = psi.NBoson;
  O.NWord = psi.NWord;
  O.MaxLength = 2;
  O.l = malloc(3 * sizeof(OperatorFixedLengthType));
  if (O.l == NULL)
  {
    printf("alloc of O.l failed in CalculateDensityMatrix\n");
    fflush(stdout);
    return 1;
  }
  O.l[0].N = 0;
  O.l[0].NBitsKey = 0;
  O.l[0].NAloc = 0;
  O.l[0].QComplex = 0;
  O.l[1].N = 0;
  O.l[1].NBitsKey = 0;
  O.l[1].NAloc = 0;
  O.l[1].QComplex = 0;
  O.l[2].N = 1;
  O.l[2].NBitsKey = 0;
  O.l[2].NAloc = 1;
  O.l[2].QComplex = 0;
  O.l[2].VRe = malloc(sizeof(double));
  O.l[2].VRe[0] = 1;
  if (O.l[2].VRe == NULL)
  {
    printf("alloc of O.l[2].VRe failed in CalculateDensityMatrix\n");
    fflush(stdout);
    return 1;
  }
  O.l[2].o = malloc(2 * sizeof(uint16_t));
  if (O.l[2].o == NULL)
  {
    printf("alloc of O.l[2].o failed in CalculateDensityMatrix\n");
    fflush(stdout);
    return 1;
  }
  if (psi.QComplex)
  {
    for (i = 0; i < psi.NFermion + psi.NBoson; i++)
    {
      O.l[2].o[1] = (i | CrAnMask);
      for (j = 0; j < psi.NFermion + psi.NBoson; j++)
      {
        O.l[2].o[0] = j;
        PsiOperatorPsiRC(O, psi, &(Rho->VRe[i][j]), &(Rho->VIm[i][j]));
      }
    }
  }
  else
  {
    for (i = 0; i < psi.NFermion + psi.NBoson; i++)
    {
      O.l[2].o[1] = (i | CrAnMask);
      for (j = 0; j < psi.NFermion + psi.NBoson; j++)
      {
        O.l[2].o[0] = j;
        PsiOperatorPsiRR(O, psi, &(Rho->VRe[i][j]));
      }
    }
  }
  return 0;
}


/*
 *  Calculate the density matrix of a wave-function
 *  On input psi contains a wave-function and Rho points to a not allocated matrix
 *  On excite Rho contains an psi.N by psi.N matrix
 *  a^{\dag}_i a_i a^{\dag}_j a_j
 */
uint32_t CalculateDensityCorrelationMatrix(MatrixType * const Rho, const WaveFunctionType psi)
{
  uint32_t i, j;
  OperatorType O;
  Rho->QComplex = psi.QComplex;
  Rho->Ni = psi.NFermion + psi.NBoson;
  Rho->Nj = psi.NFermion + psi.NBoson;
  if (InitMatrix(Rho))
  {
    printf("InitMatrix failed in CalculateDensityCorrelationMatrix\n");
    fflush(stdout);
    return 1;
  }
  strcpy(O.Name, "Calc. Den. Corr. Mat. Int. Opp. ");
  strcpy(Rho->Name, "Density Correlation Matrix.");
  O.QComplex = 0;
  O.NFermion = psi.NFermion;
  O.NBoson = psi.NBoson;
  O.NWord = psi.NWord;
  O.MaxLength = 4;
  O.l = malloc(5 * sizeof(OperatorFixedLengthType));
  if (O.l == NULL)
  {
    printf("alloc of O.l failed in CalculateDensityMatrix\n");
    fflush(stdout);
    return 1;
  }
  O.l[0].N = 0;
  O.l[0].NBitsKey = 0;
  O.l[0].NAloc = 0;
  O.l[0].QComplex = 0;
  O.l[1].N = 0;
  O.l[1].NBitsKey = 0;
  O.l[1].NAloc = 0;
  O.l[1].QComplex = 0;
  O.l[2].N = 0;
  O.l[2].NBitsKey = 0;
  O.l[2].NAloc = 0;
  O.l[2].QComplex = 0;
  O.l[3].N = 0;
  O.l[3].NBitsKey = 0;
  O.l[3].NAloc = 0;
  O.l[3].QComplex = 0;
  O.l[4].N = 1;
  O.l[4].NBitsKey = 0;
  O.l[4].NAloc = 1;
  O.l[4].QComplex = 0;
  O.l[4].VRe = malloc(sizeof(double));
  O.l[4].VRe[0] = 1;
  if (O.l[4].VRe == NULL)
  {
    printf("alloc of O.l[4].VRe failed in CalculateCorrelationDensityMatrix\n");
    fflush(stdout);
    return 1;
  }
  O.l[4].o = malloc(4 * sizeof(uint16_t));
  if (O.l[4].o == NULL)
  {
    printf("alloc of O.l[4].o failed in CalculateDensityCorrelationMatrix\n");
    fflush(stdout);
    return 1;
  }
  if (psi.QComplex)
  {
    for (i = 0; i < psi.NFermion + psi.NBoson; i++)
    {
      O.l[4].o[3] = (i | CrAnMask);
      O.l[4].o[2] = (i);
      for (j = 0; j < psi.NFermion + psi.NBoson; j++)
      {
        O.l[4].o[1] = (j | CrAnMask);
        O.l[4].o[0] = j;
        PsiOperatorPsiRC(O, psi, &(Rho->VRe[i][j]), &(Rho->VIm[i][j]));
      }
    }
  }
  else
  {
    for (i = 0; i < psi.NFermion + psi.NBoson; i++)
    {
      O.l[4].o[3] = (i | CrAnMask);
      O.l[4].o[2] = (i);
      for (j = 0; j < psi.NFermion + psi.NBoson; j++)
      {
        O.l[4].o[1] = (j | CrAnMask);
        O.l[4].o[0] = j;
        PsiOperatorPsiRR(O, psi, &(Rho->VRe[i][j]));
      }
    }
  }
  return 0;
}

/*
 *  Calculate the density correlation tensor of a wave-function
 *  On input psi contains a wave-function and Rho points to a not allocated matrix
 *  as the density correlation tensor scales to the 4-th power of the number of basis orbitals one can select a subset
 *  NFermionIncluded set the number included in the basis
 *  FermionIncluded is a list of length NFermionIncluded listing the included Fermions
 *  On excite Rho contains an NFerionInlcuded by NFerionInlcuded by NFerionInlcuded by NFerionInlcuded tensor
 *  the index of Rho is [i + NFermionIncluded j + NFermionIncluded^2 k + NFermionIncluded^3 l
 *  a^{\dag}_i a^{\dag}_j a_k a_l
 */
uint32_t CalculatePartialDensityCorrelationTensor(Complex16 * Rho, const WaveFunctionType psi, uint32_t const NFermionIncluded, uint16_t const * FermionIncluded)
{
  uint32_t i, j,k,l,m=0;
  OperatorType O;
  if(NFermionIncluded>psi.NFermion+psi.NBoson)
  {
    printf("Error in CalculatePartialDensityCorrelationTensor, NFermionIncluded=%i but psi.NFermion (psi.NBoson) is only %i (%i)\n\n",NFermionIncluded,psi.NFermion,psi.NBoson);
    fflush(stdout);
    return 1;
  }
  for(i=0;i<NFermionIncluded;i++)
  {
    if(FermionIncluded[i]>psi.NFermion+psi.NBoson)
    {
      printf("Error in CalculatePartialDensityCorrelationTensor, NFermionIncluded=%i but psi.NFermion (psi.NBoson) is only %i (%i)\n\n",FermionIncluded[i],psi.NFermion,psi.NBoson);
      fflush(stdout);
      return 1;
    }
  }
  strcpy(O.Name, "Calc. Den. Corr. Mat. Int. Opp. ");
  O.QComplex = 0;
  O.NFermion = psi.NFermion;
  O.NBoson = psi.NBoson;
  O.NWord = psi.NWord;
  O.MaxLength = 4;
  O.l = malloc(5 * sizeof(OperatorFixedLengthType));
  if (O.l == NULL)
  {
    printf("alloc of O.l failed in CalculatePartialDensityCorrelationTensor\n");
    fflush(stdout);
    return 1;
  }
  O.l[0].N = 0;
  O.l[0].NBitsKey = 0;
  O.l[0].NAloc = 0;
  O.l[0].QComplex = 0;
  O.l[1].N = 0;
  O.l[1].NBitsKey = 0;
  O.l[1].NAloc = 0;
  O.l[1].QComplex = 0;
  O.l[2].N = 0;
  O.l[2].NBitsKey = 0;
  O.l[2].NAloc = 0;
  O.l[2].QComplex = 0;
  O.l[3].N = 0;
  O.l[3].NBitsKey = 0;
  O.l[3].NAloc = 0;
  O.l[3].QComplex = 0;
  O.l[4].N = 1;
  O.l[4].NBitsKey = 0;
  O.l[4].NAloc = 1;
  O.l[4].QComplex = 0;
  O.l[4].VRe = malloc(sizeof(double));
  O.l[4].VRe[0] = 1;
  if (O.l[4].VRe == NULL)
  {
    printf("alloc of O.l[4].VRe failed in CalculatePartialDensityCorrelationTensor\n");
    fflush(stdout);
    return 1;
  }
  O.l[4].o = malloc(4 * sizeof(uint16_t));
  if (O.l[4].o == NULL)
  {
    printf("alloc of O.l[4].o failed in CalculatePartialDensityCorrelationTensor\n");
    fflush(stdout);
    return 1;
  }
  if (psi.QComplex)
  {
    for (l = 0; l < NFermionIncluded; l++)
    {
      O.l[4].o[0] = (FermionIncluded[l]);
      for (k = 0; k < NFermionIncluded; k++)
      {
        O.l[4].o[1] = (FermionIncluded[k]);
        for (j = 0; j < NFermionIncluded; j++)
        {
          O.l[4].o[2] = (FermionIncluded[j] | CrAnMask);
          for (i = 0; i < NFermionIncluded; i++)
          {
            O.l[4].o[3] = (FermionIncluded[i] | CrAnMask);
            PsiOperatorPsiRC(O, psi, &(Rho[m].Re), &(Rho[m].Im));
            m++;
          }
        }
      }
    }
  }
  else
  {
    for (l = 0; l < NFermionIncluded; l++)
    {
      O.l[4].o[0] = (FermionIncluded[l]);
      for (k = 0; k < NFermionIncluded; k++)
      {
        O.l[4].o[1] = (FermionIncluded[k]);
        for (j = 0; j < NFermionIncluded; j++)
        {
          O.l[4].o[2] = (FermionIncluded[j] | CrAnMask);
          for (i = 0; i < NFermionIncluded; i++)
          {
            O.l[4].o[3] = (FermionIncluded[i] | CrAnMask);
            PsiOperatorPsiRR(O, psi, &(Rho[m].Re));
            m++;
          }
        }
      }
    }
  }
  
  FreeOperator(&O);
  return 0;
}

/*
 *  Calculate the Spin correlation matrix
 *  On input psi contains a wave-function and Rho points to a not allocated matrix
 *  On excite Rho contains an psi.N/2 by psi.N/2 matrix
 *  needs testing TODO
 */
uint32_t CalculateSpinCorrelationMatrix(MatrixType * const Rho, const WaveFunctionType psi)
{
  uint16_t i, j;
  OperatorType O;
  Rho->QComplex = psi.QComplex;
  Rho->Ni = psi.NFermion / 2;
  Rho->Nj = psi.NFermion / 2;
  if (InitMatrix(Rho))
  {
    printf("InitMatrix failed in CalculateSpinCorrelationMatrix\n");
    fflush(stdout);
    return 1;
  }
  strcpy(O.Name, "Calc. Spin. Corr. Int. Opp. ");
  strcpy(Rho->Name, "Spin Correlation Matrix.");
  O.QComplex = 0;
  O.NFermion = psi.NFermion;
  O.NBoson = psi.NBoson;
  O.NWord = psi.NWord;
  O.MaxLength = 4;
  O.l = malloc(5 * sizeof(OperatorFixedLengthType));
  if (O.l == NULL)
  {
    printf("alloc of O.l failed in CalculateSpinCorrelationMatrix\n");
    fflush(stdout);
    return 1;
  }
  O.l[0].N = 0;
  O.l[0].NBitsKey = 0;
  O.l[0].NAloc = 0;
  O.l[0].QComplex = 0;
  O.l[1].N = 0;
  O.l[1].NBitsKey = 0;
  O.l[1].NAloc = 0;
  O.l[1].QComplex = 0;
  O.l[2].N = 0;
  O.l[2].NBitsKey = 0;
  O.l[2].NAloc = 0;
  O.l[2].QComplex = 0;
  O.l[3].N = 0;
  O.l[3].NBitsKey = 0;
  O.l[2].NAloc = 0;
  O.l[3].QComplex = 0;
  O.l[4].N = 6;
  O.l[4].NBitsKey = 0;
  O.l[2].NAloc = 6;
  O.l[4].QComplex = 0;
  O.l[4].VRe = malloc(6 * sizeof(double));
  if (O.l[4].VRe == NULL)
  {
    printf("alloc of O.l[4].VRe failed in CalculateSpinCorrelationMatrix\n");
    fflush(stdout);
    return 1;
  }
  O.l[4].VRe[0] = 0.25;
  O.l[4].VRe[1] = -0.25;
  O.l[4].VRe[2] = -0.25;
  O.l[4].VRe[3] = 0.25;
  O.l[4].VRe[4] = 0.5;
  O.l[4].VRe[5] = 0.5;
  O.l[4].o = malloc(6 * 4 * sizeof(uint16_t));
  if (O.l[4].o == NULL)
  {
    printf("alloc of O.l[4].o failed in CalculateDensityMatrix\n");
    fflush(stdout);
    return 1;
  }
  if (psi.QComplex)
  {
    for (i = 0; i < psi.NFermion; i += 2)
    {
      O.l[4].o[0 * 4 + 3] = (i | CrAnMask);
      O.l[4].o[0 * 4 + 2] = (i);
      O.l[4].o[1 * 4 + 3] = (i | CrAnMask);
      O.l[4].o[1 * 4 + 2] = (i);
      O.l[4].o[2 * 4 + 3] = (i + 1 | CrAnMask);
      O.l[4].o[2 * 4 + 2] = (i + 1);
      O.l[4].o[3 * 4 + 3] = (i + 1 | CrAnMask);
      O.l[4].o[3 * 4 + 2] = (i + 1);
      O.l[4].o[4 * 4 + 3] = (i + 1 | CrAnMask);
      O.l[4].o[4 * 4 + 2] = (i);
      O.l[4].o[5 * 4 + 3] = (i | CrAnMask);
      O.l[4].o[5 * 4 + 2] = (i + 1);
      for (j = 0; j < psi.NFermion; j += 2)
      {
        O.l[4].o[0 * 4 + 1] = (j | CrAnMask);
        O.l[4].o[0 * 4 + 0] = (j);
        O.l[4].o[1 * 4 + 1] = (j + 1 | CrAnMask);
        O.l[4].o[1 * 4 + 0] = (j + 1);
        O.l[4].o[2 * 4 + 1] = (j | CrAnMask);
        O.l[4].o[2 * 4 + 0] = (j);
        O.l[4].o[3 * 4 + 1] = (j + 1 | CrAnMask);
        O.l[4].o[3 * 4 + 0] = (j + 1);
        O.l[4].o[4 * 4 + 1] = (j | CrAnMask);
        O.l[4].o[4 * 4 + 0] = (j + 1);
        O.l[4].o[5 * 4 + 1] = (j + 1 | CrAnMask);
        O.l[4].o[5 * 4 + 0] = (j);
        PsiOperatorPsiRC(O, psi, &(Rho->VRe[i / 2][j / 2]), &(Rho->VIm[i / 2][j / 2]));
      }
    }
  }
  else
  {
    for (i = 0; i < psi.NFermion; i += 2)
    {
      O.l[4].o[0 * 4 + 0] = (i | CrAnMask);
      O.l[4].o[0 * 4 + 1] = (i);
      O.l[4].o[1 * 4 + 0] = (i | CrAnMask);
      O.l[4].o[1 * 4 + 1] = (i);
      O.l[4].o[2 * 4 + 0] = (i + 1 | CrAnMask);
      O.l[4].o[2 * 4 + 1] = (i + 1);
      O.l[4].o[3 * 4 + 0] = (i + 1 | CrAnMask);
      O.l[4].o[3 * 4 + 1] = (i + 1);
      O.l[4].o[4 * 4 + 0] = (i + 1 | CrAnMask);
      O.l[4].o[4 * 4 + 1] = (i);
      O.l[4].o[5 * 4 + 0] = (i | CrAnMask);
      O.l[4].o[5 * 4 + 1] = (i + 1);
      for (j = 0; j < psi.NFermion; j += 2)
      {
        O.l[4].o[0 * 4 + 2] = (j | CrAnMask);
        O.l[4].o[0 * 4 + 3] = (j);
        O.l[4].o[1 * 4 + 2] = (j + 1 | CrAnMask);
        O.l[4].o[1 * 4 + 3] = (j + 1);
        O.l[4].o[2 * 4 + 2] = (j | CrAnMask);
        O.l[4].o[2 * 4 + 3] = (j);
        O.l[4].o[3 * 4 + 2] = (j + 1 | CrAnMask);
        O.l[4].o[3 * 4 + 3] = (j + 1);
        O.l[4].o[4 * 4 + 2] = (j | CrAnMask);
        O.l[4].o[4 * 4 + 3] = (j + 1);
        O.l[4].o[5 * 4 + 2] = (j + 1 | CrAnMask);
        O.l[4].o[5 * 4 + 3] = (j);
/*        printf("A\n");
        PrintWaveFunctionReadable(psi);
        printf("B\n");
        fflush(stdout);*/
        PsiOperatorPsiRR(O, psi, &(Rho->VRe[i / 2][j / 2]));
      }
    }
  }
  return 0;
}


/*
 *  Calculate the density matrix of a wave-function for a selected number of particles
 *  The total number of particles is given by Nparticles
 *  The index of these particles is given by the vector j
 *  On input psi contains a wave-function and Rho points to a not allocated matrix
 *  On excite Rho contains an Nparticles by Nparticles matrix
 */
uint32_t CalculatePartialDensityMatrix(MatrixType * const Rho, const WaveFunctionType psi, const uint32_t Nparticles,
                                       const uint32_t * const j)
{
  return 0;
}


/*
 *  Calculate the diagonal of the density matrix of a wave-function
 *  On input psi contains a wave-function and Rho points to an allocated vector of lenght psi.Nboson+psi.NFermion
 *  On excite Rho contains a vector of the diagonal part of the density matrix
 */
uint32_t CalculateDensity(double *const Rho, const WaveFunctionType psi)
{
  uint32_t ipart, idet;
  for (ipart = 0; ipart < psi.NFermion; ipart++)
  {
    Rho[ipart] = 0;
    for (idet = 1; idet <= psi.N; idet++)
    {
      if (psi.Determinant[BlokElementDet(idet, psi.NWord, (ipart >> 3))] & (128 >> (ipart & 7)))
      {
        if (psi.QComplex)
        {
          Rho[ipart] += Sqr(psi.ValueRe[BlokElement(idet)]) + Sqr(psi.ValueIm[BlokElement(idet)]);
        }
        else
        {
          Rho[ipart] += Sqr(psi.ValueRe[BlokElement(idet)]);
        }
      }
    }
  }
  for (ipart = 0; ipart < psi.NBoson; ipart++)
  {
    Rho[ipart + psi.NFermion] = 0;
    for (idet = 1; idet <= psi.N; idet++)
    {
      if (psi.QComplex)
      {
        Rho[ipart + psi.NFermion] +=
          psi.Determinant[BlokElementDet(idet, psi.NWord, ipart + ((psi.NFermion + 7) >> 3))] *
          (Sqr(psi.ValueRe[BlokElement(idet)]) + Sqr(psi.ValueIm[BlokElement(idet)]));
      }
      else
      {
        Rho[ipart + psi.NFermion] +=
          psi.Determinant[BlokElementDet(idet, psi.NWord, ipart + ((psi.NFermion + 7) >> 3))] *
          Sqr(psi.ValueRe[BlokElement(idet)]);
      }
    }
  }
  return 0;
}

/*
 *
 */
uint32_t CalculateOccupationHistogram(double *const Hist, const WaveFunctionType psi, const uint16_t *fermionlist, const uint32_t N)
{
  uint8_t detFilter[psi.NWord], det;
  uint32_t i,j,k,sum;
  for (i = 0; i <= N; i++)
  {
    Hist[i] = 0;
  }
  for(i=0;i<psi.NWord;i++)
  {
    detFilter[i]=0;
  }
  for(i=0;i<N;i++)
  {
    detFilter[fermionlist[i]>>3] +=  (128 >> (fermionlist[i] & 7));
  }  
  if(psi.QComplex==0)
  {
    for(i=1;i<=psi.N;i++)
    {
      sum=0;
      for(j=0;j<psi.NWord;j++)
      {
        det=psi.Determinant[BlokElementDet(i, psi.NWord, j)]&detFilter[j];
        for(k=0;k<=7;k++)
        {
          sum+=(det>>k)&1;
        }
      }
      Hist[sum]+=Sqr(psi.ValueRe[BlokElement(i)]);
    }
  }
  else
  {
    for(i=1;i<=psi.N;i++)
    {
      sum=0;
      for(j=0;j<psi.NWord;j++)
      {
        det=psi.Determinant[BlokElementDet(i, psi.NWord, j)]&detFilter[j];
        for(k=0;k<=7;k++)
        {
          sum+=(det>>i)&1;
        }
      }
      Hist[sum]+=Sqr(psi.ValueRe[BlokElement(i)])+Sqr(psi.ValueIm[BlokElement(i)]);
    }
  }
  return 0;
}


/*
 *
 */
uint32_t CalculatePrefactorHistogram(double *const Hist, const WaveFunctionType psi, const double exp, const uint32_t N)
{
  uint32_t i, iBlok, iElement;
  int    iOMP;
  double *HistPrivate;
  for (i = 0; i < N; i++)
  {
    Hist[i] = 0;
  }
  if (psi.QComplex)
  {
#pragma omp parallel private(HistPrivate, iElement, iBlok, i)
    {
      HistPrivate = malloc(N * sizeof(double));
      if (HistPrivate == NULL)
        AllocFailed();
      for (i = 0; i < N; i++)
        HistPrivate[i] = 0;
      iBlok = psi.NBlok - 1;
#pragma omp for
      for (iOMP = 0; iOMP <= (int)(psi.N - iBlok * BlokSize); iOMP++)
      {
        iElement = (uint32_t) iOMP;
        i =
          (uint32_t) (-
                      (log(DBL_MIN + Sqr(psi.ValueRe[iBlok][iElement]) + Sqr(psi.ValueIm[iBlok][iElement])) /
                       (2 * log(exp))));
        if (i >= N)
          i = N - 1;
        HistPrivate[i]++;
      }
#pragma omp for
      for (iOMP = 0; iOMP < (int)(psi.NBlok - 1); iOMP++)
      {
        iBlok = (uint32_t) iOMP;
        for (iElement = 0; iElement < BlokSize; ++iElement)
        {
          i =
            (uint32_t) (-
                        (log(DBL_MIN + Sqr(psi.ValueRe[iBlok][iElement]) + Sqr(psi.ValueIm[iBlok][iElement])) /
                         (2 * log(exp))));
          if (i >= N)
            i = N - 1;
          HistPrivate[i]++;
        }
      }
      i = (uint32_t) (-(log(DBL_MIN + Sqr(psi.ValueRe[0][0]) + Sqr(psi.ValueIm[0][0])) / (2 * log(exp))));
      if (i >= N)
        i = N - 1;
      HistPrivate[i]--;
#pragma omp critical (AddToHistogram)
      {
        for (i = 0; i < N; i++)
          Hist[i] += HistPrivate[i];
      }
      free(HistPrivate);
    }
  }
  else
  {
#pragma omp parallel private(HistPrivate, iElement, iBlok, i)
    {
      HistPrivate = malloc(N * sizeof(double));
      if (HistPrivate == NULL)
        AllocFailed();
      for (i = 0; i < N; i++)
        HistPrivate[i] = 0;
      iBlok = psi.NBlok - 1;
#pragma omp for
      for (iOMP = 0; iOMP <= (int)(psi.N - iBlok * BlokSize); iOMP++)
      {
        iElement = (uint32_t) iOMP;
        i = (uint32_t) (-(log(DBL_MIN + fabs(psi.ValueRe[iBlok][iElement])) / log(exp)));
        if (i >= N)
          i = N - 1;
        HistPrivate[i]++;
      }
#pragma omp for
      for (iOMP = 0; iOMP < (int)(psi.NBlok - 1); iOMP++)
      {
        iBlok = (uint32_t) iOMP;
        for (iElement = 0; iElement < BlokSize; ++iElement)
        {
          i = (uint32_t) (-(log(DBL_MIN + fabs(psi.ValueRe[iBlok][iElement])) / log(exp)));
          if (i >= N)
            i = N - 1;
          HistPrivate[i]++;
        }
      }
#pragma omp critical (AddToHistogram)
      {
        for (i = 0; i < N; i++)
          Hist[i] += HistPrivate[i];
      }
      free(HistPrivate);
    }
  }
  i = (uint32_t) (-(log(DBL_MIN + fabs(psi.ValueRe[0][0])) / log(exp)));
  if (i >= N)
    i = N - 1;
  Hist[i]--;
  for (i = 0; i < N; i++)
  {
    Hist[i] /= (double)psi.N;
  }
  return 0;
}
