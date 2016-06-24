uint32_t CreateWaveFunctionDMFT(WaveFunctionType * psi, uint32_t nBathV, uint32_t nBathC, uint32_t NBitsKey)
{
  uint32_t i;
  uint8_t *Determinant;
  psi->NBitsKey = NBitsKey;
  psi->QComplex = 0;
  psi->NBoson = 0;
  psi->NFermion = 2 + 2 * nBathC + 2 * nBathV;
  if (WaveFunctionInit(psi))
  {
    printf("WaveFunctionInit failed in CreateWaveFunctionDMFT\n");
    fflush(stdout);
    return 1;
  };
  sprintf(psi->Name, "DMFT %i %i", nBathV, nBathC);
  Determinant = calloc(psi->NWord, sizeof(uint8_t));
  OperatorCr(Determinant, 1, psi->NFermion);
  for (i = 2; i < nBathV + 2; i++)
    OperatorCr(Determinant, i, psi->NFermion);
  for (i = 2 + nBathV + nBathC; i < 2 + 2 * nBathV + nBathC; i++)
    OperatorCr(Determinant, i, psi->NFermion);
  RealWaveFunctionAddElement(psi, Determinant, 1.0);
  free(Determinant);
  return 0;
}

uint32_t CreateWaveFunctionDMFTSecond(WaveFunctionType * psi, uint32_t nBathV, uint32_t nBathC, uint32_t NBitsKey)
{
  uint32_t i;
  uint8_t *Determinant;
  psi->NBitsKey = NBitsKey;
  psi->QComplex = 0;
  psi->NBoson = 0;
  psi->NFermion = 4 + 2 * nBathC + 2 * nBathV;
  if (WaveFunctionInit(psi))
  {
    printf("WaveFunctionInit failed in CreateWaveFunctionDMFT\n");
    fflush(stdout);
    return 1;
  };
  sprintf(psi->Name, "DMFT %i %i", nBathV, nBathC);
  Determinant = calloc(psi->NWord, sizeof(uint8_t));
  OperatorCr(Determinant, 0, psi->NFermion);
  OperatorCr(Determinant, 1, psi->NFermion);
  for (i = 4; i < nBathV + 4; i++)
    OperatorCr(Determinant, i, psi->NFermion);
  for (i = 4 + nBathV + nBathC; i < 4 + 2 * nBathV + nBathC; i++)
    OperatorCr(Determinant, i, psi->NFermion);
  RealWaveFunctionAddElement(psi, Determinant, 1.0);
  free(Determinant);
  return 0;
}

uint32_t CreatePESupOperatorDMFT(OperatorType * O, uint32_t NFermion)
{
  uint16_t oadd;
  InitOperator(O);
  O->NFermion = NFermion;
  O->NWord = Min(O->NFermion, 1 + (O->NFermion - 1) / 8);
  InitOperatorFixedLength(O, 1, 1);
  oadd = 0;
  RealOperatorAddLader(O, &oadd, 1, 1.0);
  return 0;
}

uint32_t CreatePESdnOperatorDMFT(OperatorType * O, uint32_t NFermion)
{
  uint16_t oadd;
  InitOperator(O);
  O->NFermion = NFermion;
  O->NWord = Min(O->NFermion, 1 + (O->NFermion - 1) / 8);
  InitOperatorFixedLength(O, 1, 1);
  oadd = 1;
  RealOperatorAddLader(O, &oadd, 1, 1.0);
  return 0;
}

uint32_t CreateIPESupOperatorDMFT(OperatorType * O, uint32_t NFermion)
{
  uint16_t oadd;
  InitOperator(O);
  O->NFermion = NFermion;
  O->NWord = Min(O->NFermion, 1 + (O->NFermion - 1) / 8);
  InitOperatorFixedLength(O, 1, 1);
  oadd = CrAnMask;
  RealOperatorAddLader(O, &oadd, 1, 1.0);
  return 0;
}

uint32_t CreateIPESdnOperatorDMFT(OperatorType * O, uint32_t NFermion)
{
  uint16_t oadd;
  InitOperator(O);
  O->NFermion = NFermion;
  O->NWord = Min(O->NFermion, 1 + (O->NFermion - 1) / 8);
  InitOperatorFixedLength(O, 1, 1);
  oadd = CrAnMask + 1;
  RealOperatorAddLader(O, &oadd, 1, 1.0);
  return 0;
}


/*  Intput:
 *    Structrue containing the DMFT settings
 *  Output:
 *    G0 in anderson format
 */
uint32_t initializeDMFTReadG0(DMFTSettings * const O, AndersonMatrixType * const G0And)
{
  if (O->NmuSym & 1)
  {
    /* read chemical potential from G0, override values set in O */
    if (FileReadAndersonMatrixTypeUndeterminedPlusMu(G0And, O->G0FileName, &(O->mu0)))
    {
      printf("Fatal Error in main\n Failed to read G0 from file %s\n", O->G0FileName);
      return 1;
    };
  }
  else
  {
    if (FileReadAndersonMatrixTypeUndetermined(G0And, O->G0FileName))
    {
      printf("Fatal Error in main\n Failed to read G0 from file %s\n", O->G0FileName);
      return 1;
    };
    G0And->mu = O->mu0;
  }
  strcpy(G0And->Name, "G0");
  if ((O->SymOpt) & 2)
  {
    printf("Symmetrize G0(w-mu0) -= G0(mu0-w)^*\n");
    SymmetrizeAndersonMatrix(G0And);
  }
  printf("\nRead G0 from file \"%s\"\n", O->G0FileName);
  PrintAndersonMatrixEssentials(*G0And);
  return 0;
}

/*  Intput:
 *    Structrue containing the DMFT settings
 *  Output:
 *    GBath in anderson and list of poles format
 */
uint32_t initializeDMFTReadGBath(DMFTSettings * const O, AndersonMatrixType * const GBathAnd,
                                 ListOfPoles * const GBathLP)
{
  char   FileName[MAXStrLength];
  sprintf(FileName, "%s%03i", O->RootFileName, O->iStart);
  if (file_exists(FileName))
  {
    if (FileReadListOfPolesUndetermined(GBathLP, FileName))
    {
      printf("Fatal Error in main\n Failed to read G0 from file %s\n", FileName);
      return 1;
    };
    strcpy(GBathLP->Name, "GBath");
    if (ListOfPolesToAndersonMatrix(*GBathLP, GBathAnd))
    {
      printf("Error in ListOfPolesToAndersonMatrix in initializeDMFTReadGBath\n");
      return 1;
    };
    strcpy(GBathAnd->Name, "GBath");
  }
  else
  {
    printf("File for bath Greenfunction %s does not excist\n", FileName);
    printf("Start from the situation where GBath=G0\n");
    initializeDMFTReadG0(O, GBathAnd);
    AndersonMatrixToListOfPoles(*GBathAnd, GBathLP);
    strcpy(GBathAnd->Name, "GBath");
    strcpy(GBathLP->Name, "GBath");
  }
  return 0;
}

/*  Intput:
 *    Structrue containing the DMFT settings
 *  Output:
 *    Sigma as a list of poles
 */
uint32_t initializeDMFTReadSigma(DMFTSettings const O, ListOfPoles * const Sigma)
{
  char   FileName[MAXStrLength];
  sprintf(FileName, "%sSigma%03i", O.RootFileName, O.iStart);
  if (file_exists(FileName))
  {
    if (FileReadListOfPolesUndetermined(Sigma, FileName))
    {
      printf("Fatal Error in main\n Failed to read G0 from file %s\n", FileName);
      return 1;
    };
    strcpy(Sigma->Name, "Sigma");
    printf("\nRead Sigma from file \"%s\"\n", FileName);
    PrintListOfPolesEssentials(*Sigma);
  }
  else
  {
    printf("File for the selfenergy %s does not excist\n", FileName);
    printf("Start from the situation where Sigma=0\n");
    Sigma->N = 1;
    InitListOfPoles(Sigma);
    Sigma->a0 = 0;
    Sigma->a[0] = 0;
    Sigma->b[0] = 0;
    strcpy(Sigma->Name, "Sigma");
    Sigma->mu = 0;
    PrintListOfPolesEssentials(*Sigma);
  }
  return 0;
}


/*  Intput:
 *    Structrue containing the DMFT settings
 *  Output:
 *    Psi, a starting wave-function for the DMFT loops
 */
uint32_t initializeDMFTReadPsi(DMFTSettings const O, WaveFunctionType * const Psi)
{
  char   FileName[MAXStrLength], Buffer[MAXStrLength];
  uint32_t PsiNBathV, PsiNBathC;
  sprintf(FileName, "%s%03i", O.WaveFunctionRootFileName, O.iStart);
  if (file_exists(FileName))
  {
    Psi->NBitsKey = O.NBitsKey;
    if (FileReadWaveFunctionPlainBinary(Psi, FileName))
    {
      printf("Fatal Error in main\n Failed to read Psi from file %s !\n", FileName);
      return 1;
    };
    printf("\nRead the starting wave function from the file %s\n", FileName);
    sscanf(Psi->Name, "%s %i %i", Buffer, &PsiNBathC, &PsiNBathV);
    printf("Wave-function \"%s\" found has %i conduction bath sites, %i valence bath sites and %i Fermions in basis\n",
           Psi->Name, PsiNBathC, PsiNBathV, Psi->NFermion);
    if ((2 * (PsiNBathC + PsiNBathV + 2)) != Psi->NFermion)
    {
      printf("The total number of Fermions is incompatible witht the bath sizes\n");
      printf("This wave-function can not be used to start from\n");
      return 1;
    }
  }
  else
  {
    printf("\nWave-function file: '%s' does not excist\n Started from new function\n", FileName);
    CreateWaveFunctionDMFTSecond(Psi, O.NBath, O.NBath, O.NBitsKey);
  }
  return 0;
}



/*  Intput:
 *    filename containing the DMFT settings
 *
 *  Output:
 *    DMFT settings obtained from the input file
 *    G0     - returned in Anderson Matrix Format
 *    Sigma  - returned as list of poles
 *    Psi    - read from file or set to starting function
 *
 *    mu0 and mu calculated and set in DMFT settings
 *
 */
uint32_t initializeDMFT(const char *const flname, DMFTSettings * O, AndersonMatrixType * const G0And,
                        ListOfPoles * const SigmaLP, WaveFunctionType * const PsiStart)
{
  if (FileReadDMFTSettings(O, flname))
  {
    printf("Fatal Error in initializeDMFT\n Failed to read DMFT Settings from file %s\n", flname);
    return 1;
  };
  if (initializeDMFTReadG0(O, G0And))
  {
    printf("Fatal Error in initializeDMFTReadG0\n Failed to read G0 from file %s\n", O->G0FileName);
    return 1;
  };
  if (initializeDMFTReadSigma(*O, SigmaLP))
  {
    printf("Fatal Error in initializeDMFTReadSigma\n Failed to read Sigma from file %s\n", O->RootFileName);
    return 1;
  };
  if (initializeDMFTReadPsi(*O, PsiStart))
  {
    printf("Fatal Error in initializeDMFTReadPsi\n Failed to read PsiStart from file %s\n",
           O->WaveFunctionRootFileName);
    return 1;
  };
  PrintDMFTSettings(*O);
  return 0;
}


/*  Input:
 *    DMFT settings
 *    GBath as a list of poles
 *  Output
 *    GBathBA as a double list representing the conduction and valence bath, inclusive a form where one orbital is taken to be bonding/anitbonding with the impurity
 */
uint32_t createGBathBA(DMFTSettings const O, ListOfPoles const GBathLP, BADoubleTriDiagonalMatrixType * const GBathBA)
{
  double Ef = 0, Epsilon = 0;

  ListOfPolesToBADoubleTridiagonalMatrix(GBathLP, GBathBA, GBathLP.N, GBathLP.N, Ef, Epsilon);
  if ((O.SymOpt) & 2)
  {
    printf("Symmetrize GBath(w) -= GBath(-w)^*\n");
    printf("Size of GAnd went from %i + %i ", GBathBA->val.N, GBathBA->con.N);
    SymmetrizeBADoubleTriDiagonalMatrix(GBathBA);
    printf("to %i + %i ", GBathBA->val.N, GBathBA->con.N);
    GBathBA->val.N = Min(GBathBA->val.N, O.NGBath);
    GBathBA->con.N = Min(GBathBA->con.N, O.NGBath);
    printf("to %i + %i\n", GBathBA->val.N, GBathBA->con.N);
  }
  if ((O.SymOpt) & 16)
  {
    printf("Enforce GBathBA to fullfill the knonw sumrules\n");
    printf("Not yet implemented\n");
  }
  strcpy(GBathBA->Name, "GBath");
  return 0;
}

/*  Input:
 *    DMFT settings
 *    GBathBA as a double list representing the conduction and valence bath, inclusive a form where one orbital is taken to be bonding/anitbonding with the impurity
 *  Output
 *    GBathBATrunc as a copy of GBathBA with the chains the lengt included in the groundstate cluster
 */
uint32_t createGBathBATrunc(DMFTSettings const O, BADoubleTriDiagonalMatrixType const GBathBA,
                            BADoubleTriDiagonalMatrixType * const GBathBATrunk)
{
  *GBathBATrunk = GBathBA;
  if (GBathBA.con.N < O.NBath + 1)
  {
    printf("Warning: GBathBA.con.N(%i)<DMFTSettings.NBath+1(%i+1)\n", GBathBA.con.N, O.NBath);
  }
  else
  {
    GBathBATrunk->con.N = O.NBath + 1;
  }
  if (GBathBA.val.N < O.NBath + 1)
  {
    printf("Warning: GBathBA.val.N(%i)<DMFTSettings.NBath+1(%i+1)\n", GBathBA.val.N, O.NBath);
  }
  else
  {
    GBathBATrunk->val.N = O.NBath + 1;
  }
  strcpy(GBathBATrunk->Name, "GBath for ground-state");
  return 0;
}

/*  Input:
 *    DMFT settings
 *    GBathBA as a double list representing the conduction and valence bath, inclusive a form where one orbital is taken to be bonding/anitbonding with the impurity
 *  Output
 *    VPerturbation - operator describing the coupling of the bath included for the calculation of the ground-state to the rest of the world
 *    gPerturbation - that part of the bath Greens function not included in GBathBATrunc
 */
uint32_t createVPerturbation(DMFTSettings const O, BADoubleTriDiagonalMatrixType const GBathBA,
                             OperatorType * const VPerturbation, TriDiagonalMatrixType * const gPerturbation)
{
  uint32_t i;
  uint16_t oadd[1];
  for (i = 0; i < 2 * O.NSite; i++)
  {
    InitOperator(&(VPerturbation[i]));
    VPerturbation[i].NFermion = 4 + 2 * (Min(-1 + GBathBA.val.N, O.NBath) + Min(-1 + GBathBA.con.N, O.NBath));
    VPerturbation[i].NBoson = 0;
    VPerturbation[i].NWord = Min(VPerturbation[i].NFermion, 1 + (VPerturbation[i].NFermion - 1) / 8);
    VPerturbation[i].QComplex = 0;
    InitOperatorFixedLength(&(VPerturbation[i]), 1, 1);
  }
  if (GBathBA.val.N > (O.NBath + 1))
  {
    oadd[0] = CrAnMask + 3 + O.NBath;
    RealOperatorAddLader(&(VPerturbation[0]), oadd, 1, GBathBA.val.b[O.NBath + 1]);
    oadd[0] = CrAnMask + 3 + Min(-1 + GBathBA.con.N, O.NBath) + 2 * O.NBath;
    RealOperatorAddLader(&(VPerturbation[1]), oadd, 1, GBathBA.val.b[O.NBath + 1]);
  }
  if (GBathBA.con.N > O.NBath + 1)
  {
    oadd[0] = 3 + Min(-1 + GBathBA.val.N, O.NBath) + O.NBath;
    RealOperatorAddLader(&(VPerturbation[2]), oadd, 1, GBathBA.con.b[O.NBath + 1]);
    oadd[0] = 3 + 2 * Min(-1 + GBathBA.val.N, O.NBath) + 2 * O.NBath;
    RealOperatorAddLader(&(VPerturbation[3]), oadd, 1, GBathBA.con.b[O.NBath + 1]);
  }
  /* Create gPerturbation */
  if (GBathBA.val.N > O.NBath + 1)
  {
    gPerturbation[0].N = GBathBA.val.N - O.NBath - 1;
    gPerturbation[1].N = GBathBA.val.N - O.NBath - 1;
    InitTriDiagonalMatrix(&(gPerturbation[0]));
    InitTriDiagonalMatrix(&(gPerturbation[1]));
    for (i = 1; i < GBathBA.val.N - O.NBath - 1; i++)
    {
      gPerturbation[0].a[i] = -GBathBA.val.a[i + O.NBath + 1];
      gPerturbation[0].b[i] = GBathBA.val.b[i + O.NBath + 1];
      gPerturbation[1].a[i] = -GBathBA.val.a[i + O.NBath + 1];
      gPerturbation[1].b[i] = GBathBA.val.b[i + O.NBath + 1];
    }
    gPerturbation[0].a[GBathBA.val.N - O.NBath - 1] = -GBathBA.val.a[GBathBA.val.N];
    gPerturbation[1].a[GBathBA.val.N - O.NBath - 1] = -GBathBA.val.a[GBathBA.val.N];
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
  if (GBathBA.con.N > O.NBath + 1)
  {
    gPerturbation[2].N = GBathBA.con.N - O.NBath - 1;
    gPerturbation[3].N = GBathBA.con.N - O.NBath - 1;
    InitTriDiagonalMatrix(&(gPerturbation[2]));
    InitTriDiagonalMatrix(&(gPerturbation[3]));
    for (i = 1; i < GBathBA.con.N - O.NBath - 1; i++)
    {
      gPerturbation[2].a[i] = GBathBA.con.a[i + O.NBath + 1];
      gPerturbation[2].b[i] = GBathBA.con.b[i + O.NBath + 1];
      gPerturbation[3].a[i] = GBathBA.con.a[i + O.NBath + 1];
      gPerturbation[3].b[i] = GBathBA.con.b[i + O.NBath + 1];
    }
    gPerturbation[2].a[GBathBA.con.N - O.NBath - 1] = GBathBA.con.a[GBathBA.con.N];
    gPerturbation[3].a[GBathBA.con.N - O.NBath - 1] = GBathBA.con.a[GBathBA.con.N];
  }
  else
  {
    gPerturbation[2].N = 1;
    gPerturbation[3].N = 1;
    InitTriDiagonalMatrix(&(gPerturbation[2]));
    InitTriDiagonalMatrix(&(gPerturbation[3]));
    gPerturbation[2].a[1] = 0;
    gPerturbation[3].a[1] = 0;
  }
  return 0;
}


/*  Input:
 *    DMFT settings
 *    GBath as a list of poles
 *  Output
 *    GBathBA        as a double list representing the conduction and valence bath, inclusive a form where one orbital is taken to be bonding/anitbonding with the impurity
 *    GBathBATrunc   as a double list representing the conduction and valence bath, inclusive a form where one orbital is taken to be bonding/anitbonding with the impurity
 *                   truncated to a length included in the ground-state cluster
 *    VPerturbation
 */
uint32_t createGBathBAandVPerturbation(DMFTSettings const O, ListOfPoles const GBathLP,
                                       BADoubleTriDiagonalMatrixType * const GBathBA,
                                       BADoubleTriDiagonalMatrixType * const GBathBATrunc,
                                       OperatorType * const VPerturbation, TriDiagonalMatrixType * const gPerturbation)
{
  createGBathBA(O, GBathLP, GBathBA);
  createGBathBATrunc(O, *GBathBA, GBathBATrunc);
  createVPerturbation(O, *GBathBA, VPerturbation, gPerturbation);
  return 0;
}

uint32_t createGBathBATruncandVPerturbation(DMFTSettings const O, BADoubleTriDiagonalMatrixType const GBathBA,
                                            BADoubleTriDiagonalMatrixType * const GBathBATrunc,
                                            OperatorType * const VPerturbation,
                                            TriDiagonalMatrixType * const gPerturbation)
{
  createGBathBATrunc(O, GBathBA, GBathBATrunc);
  createVPerturbation(O, GBathBA, VPerturbation, gPerturbation);
  return 0;
}

/*  Input
 *    GBathBATrunc contains the size of the really used bath for the groundstate.
 *      This can be smaller than the size in the DMFTSettings, as the bath can be converged with less poles
 *    Psi - which is a function converged or used as an input from file good to start from
 *  Output
 *    Psi with an enlarged or smaler basis
 */
uint32_t TestIfPsiHasCorrectDimensions(BADoubleTriDiagonalMatrixType const GBathBATrunc, WaveFunctionType * const Psi)
{
  char   Buffer[MAXStrLength];
  uint32_t PsiNBathV, PsiNBathC;
  WaveFunctionType PsiTemp;
  uint32_t i;
  sscanf(Psi->Name, "%s %i %i", Buffer, &PsiNBathC, &PsiNBathV);
  if (PsiNBathC != GBathBATrunc.con.N - 1)
  {
    printf("\nBefore diagonalizing:\nThe wavefunction has %i conduction bath sites\nThe Anderson Hamiltonian %i\n",
           PsiNBathC, GBathBATrunc.con.N - 1);
    if (PsiNBathC > GBathBATrunc.con.N - 1)
    {
      printf("Removing %i unoccupied fermions from the wave-function\n", PsiNBathC - GBathBATrunc.con.N + 1);
      for (i = 0; i < PsiNBathC - GBathBATrunc.con.N + 1; i++)
      {
        WaveFunctionRemoveFermion(*Psi, &PsiTemp, Psi->NFermion - 1);
        WaveFunctionFree(Psi);
        WaveFunctionRemoveFermion(PsiTemp, Psi, 4 + PsiNBathV + PsiNBathC - i - 1);
        WaveFunctionFree(&PsiTemp);
      }
    }
    else
    {
      printf("Adding %i unoccupied fermions to the wave-function\n", GBathBATrunc.con.N - 1 - PsiNBathC);
      for (i = 0; i < GBathBATrunc.con.N - 1 - PsiNBathC; i++)
      {
        WaveFunctionAddFermion(*Psi, &PsiTemp, Psi->NFermion, 0);
        WaveFunctionFree(Psi);
        WaveFunctionAddFermion(PsiTemp, Psi, 4 + PsiNBathV + PsiNBathC + i, 0);
        WaveFunctionFree(&PsiTemp);
      }
    }
  }
  if (PsiNBathV != GBathBATrunc.val.N - 1)
  {
    printf("\nBefore diagonalizing:\nThe wavefunction has %i valence bath sites\nThe Anderson Hamiltonian %i\n",
           PsiNBathV, GBathBATrunc.val.N - 1);
    if (PsiNBathV > GBathBATrunc.val.N - 1)
    {
      printf("Removing %i occupied fermions from the wave-function\n", PsiNBathV - GBathBATrunc.val.N + 1);
      for (i = 0; i < PsiNBathV - GBathBATrunc.val.N + 1; i++)
      {
        WaveFunctionRemoveFermion(*Psi, &PsiTemp, Psi->NFermion - GBathBATrunc.con.N);
        WaveFunctionFree(Psi);
        WaveFunctionRemoveFermion(PsiTemp, Psi, 4 + PsiNBathV - i - 1);
        WaveFunctionFree(&PsiTemp);
      }
    }
    else
    {
      printf("Adding %i occupied fermions to the wave-function\n", GBathBATrunc.val.N - 1 - PsiNBathV);
      for (i = 0; i < GBathBATrunc.val.N - 1 - PsiNBathV; i++)
      {
        WaveFunctionAddFermion(*Psi, &PsiTemp, Psi->NFermion - GBathBATrunc.val.N + 1, 1);
        WaveFunctionFree(Psi);
        WaveFunctionAddFermion(PsiTemp, Psi, 4 + PsiNBathV + i, 1);
        WaveFunctionFree(&PsiTemp);
      }
    }
    printf("Corrected basis\n");
    PrintWaveFunctionReadableTruncated(*Psi, 0.1);
  }
  sprintf(Psi->Name, "DMFT %i %i", GBathBATrunc.con.N - 1, GBathBATrunc.val.N - 1);
  return 0;
}

/*  Input
 *    G0LP
 *    GBathBA
 */
uint32_t PrintInfoOnGBathAndG0InDMFTLoop(ListOfPoles const G0LP, BADoubleTriDiagonalMatrixType const GBathBA)
{
  BADoubleTriDiagonalMatrixType G0BA;
  ListOfPolesToBADoubleTridiagonalMatrix(G0LP, &G0BA, G0LP.N, G0LP.N, 0, 0);
  printf
    ("                               G0v                    G0c                    GBathv                 GBathc\n");
  printf("b0^2                        %22.15E %22.15E %22.15E %22.15E\n", Sqr(G0BA.val.b[0]), Sqr(G0BA.con.b[0]),
         Sqr(GBathBA.val.b[0]), Sqr(GBathBA.con.b[0]));
  printf("b0                          %22.15E %22.15E %22.15E %22.15E\n", G0BA.val.b[0], G0BA.con.b[0],
         GBathBA.val.b[0], GBathBA.con.b[0]);
  printf("b1                          %22.15E %22.15E %22.15E %22.15E\n", G0BA.val.b[1], G0BA.con.b[1],
         GBathBA.val.b[1], GBathBA.con.b[1]);
  printf("b2                          %22.15E %22.15E %22.15E %22.15E\n", G0BA.val.b[2], G0BA.con.b[2],
         GBathBA.val.b[2], GBathBA.con.b[2]);
  printf("a1                          %22.15E %22.15E %22.15E %22.15E\n", G0BA.val.a[1], G0BA.con.a[1],
         GBathBA.val.a[1], GBathBA.con.a[1]);
  printf("a2                          %22.15E %22.15E %22.15E %22.15E\n", G0BA.val.a[2], G0BA.con.a[2],
         GBathBA.val.a[2], GBathBA.con.a[2]);
  printf("e_i                         %22.15E                        %22.15E\n", G0BA.acvB, GBathBA.acvB);
  printf("e_i^p                       %22.15E                        %22.15E\n", G0BA.acvA, GBathBA.acvA);
  printf("t_b                         %22.15E                        %22.15E\n", G0BA.acvBA, GBathBA.acvBA);
  printf("t_c                         %22.15E                        %22.15E\n", G0BA.bcB0, GBathBA.bcB0);
  printf("t_c^p                       %22.15E                        %22.15E\n", G0BA.bcA0, GBathBA.bcA0);
  printf("t_v                         %22.15E                        %22.15E\n", G0BA.bvB0, GBathBA.bvB0);
  printf("t_v^p                       %22.15E                        %22.15E\n", G0BA.bvA0, GBathBA.bvA0);
  printf("(t_c)^2+(t_c^p)^2           %22.15E                        %22.15E\n", Sqr(G0BA.bcB0) + Sqr(G0BA.bcA0),
         Sqr(GBathBA.bcB0) + Sqr(GBathBA.bcA0));
  printf("(t_v)^2+(t_v^p)^2           %22.15E                        %22.15E\n", Sqr(G0BA.bvB0) + Sqr(G0BA.bvA0),
         Sqr(GBathBA.bvB0) + Sqr(GBathBA.bvA0));
  printf("(t_c)^2+(t_v)^2+(t_b)^2     %22.15E                        %22.15E\n",
         Sqr(G0BA.bcB0) + Sqr(G0BA.bvB0) + Sqr(G0BA.acvBA), Sqr(GBathBA.bcB0) + Sqr(GBathBA.bvB0) + Sqr(GBathBA.acvBA));
  printf("(t_c^p)^2+(t_v^p)^2+(t_b)^2 %22.15E                        %22.15E\n",
         Sqr(G0BA.bcA0) + Sqr(G0BA.bvA0) + Sqr(G0BA.acvBA), Sqr(GBathBA.bcA0) + Sqr(GBathBA.bvA0) + Sqr(GBathBA.acvBA));
  FreeBADoubleTriDiagonalMatrix(&G0BA);
  return 0;
}

uint32_t PrintInfoOnGBathAndG0InDMFTFindGroundstate(AndersonMatrixType const G0And,
                                                    BADoubleTriDiagonalMatrixType const GBathBA)
{
  double occupation;
  BADoubleTriDiagonalMatrixType G0BA;
  AndersonMatrixToBADoubleTriDiagonalMatrix(G0And, 0, &G0BA, &occupation);
  printf
    ("                               G0v                    G0c                    GBathv                 GBathc\n");
  printf("N           %22i %22i %22i %22i\n", G0BA.val.N, G0BA.con.N, GBathBA.val.N, GBathBA.con.N);
  printf("b0^2                        %22.15E %22.15E %22.15E %22.15E\n", Sqr(G0BA.val.b[0]), Sqr(G0BA.con.b[0]),
         Sqr(GBathBA.val.b[0]), Sqr(GBathBA.con.b[0]));
  printf("b0                          %22.15E %22.15E %22.15E %22.15E\n", G0BA.val.b[0], G0BA.con.b[0],
         GBathBA.val.b[0], GBathBA.con.b[0]);
  printf("b1                          %22.15E %22.15E %22.15E %22.15E\n", G0BA.val.b[1], G0BA.con.b[1],
         GBathBA.val.b[1], GBathBA.con.b[1]);
  if ((GBathBA.con.N > 2) && (GBathBA.val.N > 2))
  {
    printf("b2                          %22.15E %22.15E %22.15E %22.15E\n", G0BA.val.b[2], G0BA.con.b[2],
           GBathBA.val.b[2], GBathBA.con.b[2]);
  }
  printf("a1                          %22.15E %22.15E %22.15E %22.15E\n", G0BA.val.a[1], G0BA.con.a[1],
         GBathBA.val.a[1], GBathBA.con.a[1]);
  if ((GBathBA.con.N > 2) && (GBathBA.val.N > 2))
  {
    printf("a2                          %22.15E %22.15E %22.15E %22.15E\n", G0BA.val.a[2], G0BA.con.a[2],
           GBathBA.val.a[2], GBathBA.con.a[2]);
  }
  printf("e_i                         %22.15E                        %22.15E\n", G0BA.acvB, GBathBA.acvB);
  printf("e_i^p                       %22.15E                        %22.15E\n", G0BA.acvA, GBathBA.acvA);
  printf("t_b                         %22.15E                        %22.15E\n", G0BA.acvBA, GBathBA.acvBA);
  printf("t_c                         %22.15E                        %22.15E\n", G0BA.bcB0, GBathBA.bcB0);
  printf("t_c^p                       %22.15E                        %22.15E\n", G0BA.bcA0, GBathBA.bcA0);
  printf("t_v                         %22.15E                        %22.15E\n", G0BA.bvB0, GBathBA.bvB0);
  printf("t_v^p                       %22.15E                        %22.15E\n", G0BA.bvA0, GBathBA.bvA0);
  printf("(t_c)^2+(t_c^p)^2           %22.15E                        %22.15E\n", Sqr(G0BA.bcB0) + Sqr(G0BA.bcA0),
         Sqr(GBathBA.bcB0) + Sqr(GBathBA.bcA0));
  printf("(t_v)^2+(t_v^p)^2           %22.15E                        %22.15E\n", Sqr(G0BA.bvB0) + Sqr(G0BA.bvA0),
         Sqr(GBathBA.bvB0) + Sqr(GBathBA.bvA0));
  printf("(t_c)^2+(t_v)^2+(t_b)^2     %22.15E                        %22.15E\n",
         Sqr(G0BA.bcB0) + Sqr(G0BA.bvB0) + Sqr(G0BA.acvBA), Sqr(GBathBA.bcB0) + Sqr(GBathBA.bvB0) + Sqr(GBathBA.acvBA));
  printf("(t_c^p)^2+(t_v^p)^2+(t_b)^2 %22.15E                        %22.15E\n",
         Sqr(G0BA.bcA0) + Sqr(G0BA.bvA0) + Sqr(G0BA.acvBA), Sqr(GBathBA.bcA0) + Sqr(GBathBA.bvA0) + Sqr(GBathBA.acvBA));
  FreeBADoubleTriDiagonalMatrix(&G0BA);
  return 0;
}

/*
 *
 */
uint32_t FindGroundStateFromG0andSigma(DMFTSettings const O, AndersonMatrixType const G0, ListOfPoles const Sigma,
                                       double *const mu, double *const Nelectron, WaveFunctionType * const Psi,
                                       BADoubleTriDiagonalMatrixType * const GBathBARet,
                                       AndersonMatrixType * const GBathAndRet, TimeType * const Time,
                                       double *const onsiteshiftstart)
{
  uint32_t i, j, converged = 0, iiter;
  uint32_t const nitermax = 50, nitermem = 60;
  uint16_t FermionIncluded[2];
  AndersonMatrixType GBathAnd;
  BADoubleTriDiagonalMatrixType GBathBA, GBathBATrunc;
  double WaveFunctionAccuracy, PsiEnergy, PsiQuality, onsiteshift[nitermax + 1], occupation, Rho[4 + 4 * O.NBath],
    occlist[nitermax + 1], ABocclist[nitermax + 1], mulist[nitermax + 1];
  OperatorType HAnd, Sx, Sy, Sz, S2;
  uint16_t *SpinUp, *SpinDn;
  WaveFunctionAccuracy = O.zeroLancWavefunction;
  onsiteshift[0] = *onsiteshiftstart;
  mulist[0] = *mu;
  iiter = 0;
  GBathAndListAndRepresentation(G0, *mu, Sigma, &GBathAnd);
  /*
     The bath Greens function is not symmetric around mu
     if((O.SymOpt)&2)
     {
     printf("Symmetrize the bath greens function around mu=%22.15E\n",GBathAnd.mu);
     SymmetrizeAndersonMatrix(&GBathAnd);
     }
   */
  ReducePolesAndersonMatrixNumberFermi(&GBathAnd, O.NGBath, *mu);
  while (!converged)
  {
    GBathAnd.mu = mulist[iiter];
    AndersonMatrixToBADoubleTriDiagonalMatrixTemperature(GBathAnd, onsiteshift[iiter], O.Teff, &GBathBA, &occupation);
    GBathBATrunc = GBathBA;
    GBathBATrunc.con.N = Min(GBathBA.con.N, O.NBath + 1);
    GBathBATrunc.val.N = Min(GBathBA.val.N, O.NBath + 1);
    TestIfPsiHasCorrectDimensions(GBathBATrunc, Psi);
    BADoubleTriDiagonalMatrixToAndersonImpurityOperator(&HAnd, GBathBATrunc, GBathBATrunc, O.U, 0);
    if ((O.SymOpt) & 4)
    {
      printf("Restricted HAnd for the ground-state to the first bath site\n");
      HAnd.l[2].N = 8;
    }
    /*
       The bath Greens function is not symmetric around mu
       if((O.SymOpt)&2)
       {
       printf("Symmetrize the bath greens function around mu=%22.15E\n",GBathBA.mu);
       SymmetrizeBADoubleTriDiagonalMatrix(&GBathBA);
       }
     */
    /* if singlet set is true add the (U+1) * S^2 to the operator HAndGrd */
    if ((O.SymOpt) & 1)
    {
      printf("Add %22.15E * S^2 to HAnd for the calculation of the ground-state wave-function\n", 1 + O.U);
      printf("This stabilizes the ground-state by enlarging the singlet triplet energy\n");
      SpinUp = calloc(HAnd.NFermion / 2, sizeof(uint16_t));
      SpinDn = calloc(HAnd.NFermion / 2, sizeof(uint16_t));
      SpinUp[0] = 0;
      SpinUp[1] = 2;
      SpinDn[0] = 1;
      SpinDn[1] = 3;
      for (i = 4; i < 2 + HAnd.NFermion / 2; i++)
      {
        SpinUp[i - 2] = i;
        SpinDn[i - 2] = i + HAnd.NFermion / 2 - 2;
      }
      CreateSpinOperators(HAnd.NFermion, HAnd.NFermion / 2, SpinUp, SpinDn, &Sx, &Sy, &Sz, &S2);
      FreeOperator(&Sx);
      FreeOperator(&Sy);
      FreeOperator(&Sz);
      OperatorAddNormalOrder(&HAnd, S2, fabs(O.U) + 1, 0, &QDetNotOrderedFcFaBcBa);
      FreeOperator(&S2);
    }
    /* create ground-state wave-function */
    BeginTiming("Lanczos GroundState", Time);
    LanczosGroundState(HAnd, Psi, &WaveFunctionAccuracy, O.zeroLancEnergy, O.U + *mu, Psi, 0);
    EndTiming("Lanczos GroundState", Time);
    /* calculate NelectronNew */
    CalculateDensity(Rho, *Psi);
    ABocclist[iiter] = Rho[0] + Rho[1] + Rho[2] + Rho[3];
    occlist[iiter] = Rho[0] + Rho[1];
    /* determine if ground-state calculation has to be restarted */
    if (O.NmuSym & 2)
    {
      /* calculations at constant mu, but still need to make sure the bath Greens function is calculated at the same filling as the impurity greens function */
      /* this is done by optimizing NelectronNew such that the occupation of the impurity and the selected bath site is integer */
      if (Abs(ABocclist[iiter] - O.NSite) > 100 * DBL_EPSILON)
      {
        mulist[iiter + 1] = *mu;
        if (iiter == 0)
        {
          onsiteshift[iiter + 1] = onsiteshift[iiter] - (ABocclist[iiter] - O.NSite) * (Abs(G0.b[1]) + Abs(O.U / 2));
        }
        else if (iiter == nitermax)
        {
          printf("occupation of bath and impurity green function did not converge in %i itterations\n", iiter);
          printf("loop   onsiteshift            total occupation       mu                     impurity occupation\n");
          for (j = 0; j < iiter; j++)
          {
            printf("%3i %22.15E %22.15E %22.15E %22.15E\n", j, onsiteshift[j], ABocclist[j], mulist[j], occlist[j]);
          }
          converged = 1;
          /* pass the converged onsite shift back to onsiteshift start as the next starting point */
          *onsiteshiftstart = onsiteshift[iiter];
        }
        else
        {
          onsiteshift[iiter + 1] =
            InterpolateRootOfset(&(onsiteshift[Max(iiter, nitermem) - nitermem]),
                                 &(ABocclist[Max(iiter, nitermem) - nitermem]), Min(iiter + 1, nitermem), O.NSite);
        }
        iiter++;
      }
      else
      {
        converged = 1;
        /* pass the converged onsite shift back to onsiteshift start as the next starting point */
        *onsiteshiftstart = onsiteshift[iiter];
      }
    }
    else
    {
      /* calculations at constant occupation, find the correct mu and onsiteshift */
      if ((Abs(ABocclist[iiter] - O.NSite) > 100 * DBL_EPSILON)
          || (Abs(occlist[iiter] - *Nelectron) > 100 * DBL_EPSILON))
      {
        if (iiter == 0)
        {
          onsiteshift[iiter + 1] = onsiteshift[iiter] - (ABocclist[iiter] - O.NSite) * (Abs(G0.b[1]) + Abs(O.U / 2));
          mulist[iiter + 1] = mulist[iiter];
        }
        else if (iiter == nitermax)
        {
          printf("occupation of bath and impurity green function did not converge in %i itterations\n", iiter);
          printf("loop   onsiteshift            total occupation       mu                     impurity occupation\n");
          for (j = 0; j < iiter; j++)
          {
            printf("%3i %22.15E %22.15E %22.15E %22.15E\n", j, onsiteshift[j], ABocclist[j], mulist[j], occlist[j]);
          }
          converged = 1;
          *mu = mulist[nitermax];
        }
        else
        {
          InterpolateRoot2DOfset(&(onsiteshift[Max(iiter, nitermem) - nitermem]),
                                 &(mulist[Max(iiter, nitermem) - nitermem]),
                                 &(ABocclist[Max(iiter, nitermem) - nitermem]),
                                 &(occlist[Max(iiter, nitermem) - nitermem]), Min(iiter + 1, nitermem), O.NSite,
                                 *Nelectron, &(onsiteshift[iiter + 1]), &(mulist[iiter + 1]));
        }
        iiter++;
      }
      else
      {
        converged = 1;
        *mu = mulist[iiter];
      }
    }
    if (!converged)
      FreeOperator(&HAnd);
  }
  /* print info on groundstate */
  PrintInfoOnGBathAndG0InDMFTFindGroundstate(G0, GBathBA);
  PsiHermitianOperatorPsi(HAnd, *Psi, &PsiEnergy);
  PsiHermitianOperatorSqrPsi(HAnd, *Psi, &PsiQuality);
  PsiQuality /= Sqr(PsiEnergy);
  printf("Energy                  <psi|H^2|psi>/(<psi|H|psi>^2)\n");
  printf("%22.15E  %22.15E\n", PsiEnergy, PsiQuality);
  PrintWaveFunctionDensity(*Psi);
  PrintWaveFunctionPrefactorHistogram(*Psi, 2.0, 20);
  FermionIncluded[0]=0;
  FermionIncluded[1]=1;
  PrintWaveFunctionPartialCorrelationTensor(*Psi, 2, FermionIncluded);
  *GBathBARet = GBathBA;
  *GBathAndRet = GBathAnd;
  return 0;
}


/*
 *
 */
uint32_t FindGroundStateFromG0andSigmaDoubleLoop(DMFTSettings const O, AndersonMatrixType const G0,
                                                 ListOfPoles const Sigma, double *const mu, double *const Nelectron,
                                                 WaveFunctionType * const Psi,
                                                 BADoubleTriDiagonalMatrixType * const GBathBARet,
                                                 AndersonMatrixType * const GBathAndRet, TimeType * const Time)
{
  uint32_t i, converged = 0, nitermax = 30, iiter;
  AndersonMatrixType GBathAnd;
  BADoubleTriDiagonalMatrixType GBathBA, GBathBATrunc;
  double WaveFunctionAccuracy, PsiEnergy, PsiQuality, onsiteshift, occupation, NelectronNew, Rho[4 + 4 * O.NBath],
    occlist[nitermax + 1], mulist[nitermax + 1];
  double const *occin = &(occlist[0]), *occout = &(occlist[1]);
  OperatorType HAnd, Sx, Sy, Sz, S2;
  uint16_t *SpinUp, *SpinDn;
  WaveFunctionAccuracy = O.zeroLancWavefunction;

  onsiteshift = *mu;
  occlist[0] = *Nelectron;
  mulist[0] = *mu;
  iiter = 0;
  while (!converged)
  {
    GBathAndListAndRepresentation(G0, *mu, Sigma, &GBathAnd);
    ReducePolesAndersonMatrixNumberFermi(&GBathAnd, O.NGBath, *mu);
    AndersonMatrixFindOnsiteEnergyShift(GBathAnd, *Nelectron / 2.0, &onsiteshift);
    AndersonMatrixToBADoubleTriDiagonalMatrix(GBathAnd, onsiteshift, &GBathBA, &occupation);
    GBathBATrunc = GBathBA;
    GBathBATrunc.con.N = Min(GBathBA.con.N, O.NBath + 1);
    GBathBATrunc.val.N = Min(GBathBA.val.N, O.NBath + 1);
    TestIfPsiHasCorrectDimensions(GBathBATrunc, Psi);
    BADoubleTriDiagonalMatrixToAndersonImpurityOperator(&HAnd, GBathBATrunc, GBathBATrunc, O.U, 0);
    if ((O.SymOpt) & 4)
    {
      printf("Restricted HAnd for the ground-state to the first bath site\n");
      HAnd.l[2].N = 8;
    }
    /* if singlet set is true add the (U+1) * S^2 to the operator HAndGrd */
    if ((O.SymOpt) & 1)
    {
      printf("Add %22.15E * S^2 to HAnd for the calculation of the ground-state wave-function\n", 1 + O.U);
      printf("This stabilizes the ground-state by enlarging the singlet triplet energy\n");
      SpinUp = calloc(HAnd.NFermion / 2, sizeof(uint16_t));
      SpinDn = calloc(HAnd.NFermion / 2, sizeof(uint16_t));
      SpinUp[0] = 0;
      SpinUp[1] = 2;
      SpinDn[0] = 1;
      SpinDn[1] = 3;
      for (i = 4; i < 2 + HAnd.NFermion / 2; i++)
      {
        SpinUp[i - 2] = i;
        SpinDn[i - 2] = i + HAnd.NFermion / 2 - 2;
      }
      CreateSpinOperators(HAnd.NFermion, HAnd.NFermion / 2, SpinUp, SpinDn, &Sx, &Sy, &Sz, &S2);
      FreeOperator(&Sx);
      FreeOperator(&Sy);
      FreeOperator(&Sz);
      OperatorAddNormalOrder(&HAnd, S2, fabs(O.U) + 1, 0, &QDetNotOrderedFcFaBcBa);
      FreeOperator(&S2);
    }
    /* create ground-state wave-function */
    BeginTiming("Lanczos GroundState", Time);
    LanczosGroundState(HAnd, Psi, &WaveFunctionAccuracy, O.zeroLancEnergy, O.U, Psi, 0);
    EndTiming("Lanczos GroundState", Time);
    /* calculate NelectronNew */
    CalculateDensity(Rho, *Psi);
    NelectronNew = Rho[0] + Rho[1];
    /* determine if ground-state calculation has to be restarted */
    if (O.NmuSym & 2)
    {
      /* calculations at constant mu, but still need to make sure the bath Greens function is calculated at the same filling as the impurity greens function */
      if (Abs(NelectronNew - *Nelectron) > 100 * DBL_EPSILON)
      {
        iiter++;
        occlist[iiter] = NelectronNew;
        if (iiter == 1)
        {
          *Nelectron = NelectronNew;
        }
        else if (iiter == nitermax)
        {
          printf("occupation of bath and impurity green function did not converge in %i itterations\n", iiter);
          printf("Nbath = %22.15E and Nimp = %22.15E\n", *Nelectron, NelectronNew);
          converged = 1;
        }
        else
        {
          *Nelectron = FindSelfconsistentnumber(occin, occout, iiter - 1);
        }
      }
      else
      {
        converged = 1;
      }
    }
    else
    {
      /* calculations at constant occupation, find the correct mu */
      printf("Not jet implemented\n");

    }
    if (!converged)
      FreeOperator(&HAnd);
  }

  /* print info on groundstate */

  PrintInfoOnGBathAndG0InDMFTFindGroundstate(G0, GBathBA);

  PsiHermitianOperatorPsi(HAnd, *Psi, &PsiEnergy);
  PsiHermitianOperatorSqrPsi(HAnd, *Psi, &PsiQuality);
  PsiQuality /= Sqr(PsiEnergy);
  printf("Energy                  <psi|H^2|psi>/(<psi|H|psi>^2)\n");
  printf("%22.15E  %22.15E\n", PsiEnergy, PsiQuality);
  PrintWaveFunctionDensity(*Psi);
  PrintWaveFunctionPrefactorHistogram(*Psi, 2.0, 20);
  *GBathBARet = GBathBA;
  *GBathAndRet = GBathAnd;
  return 0;
}

/*
 *
 */
uint32_t FindGroundStateandMu(DMFTSettings const O, BADoubleTriDiagonalMatrixType const GAndTrunk, double *const mu,
                              WaveFunctionType * const Psi, OperatorType * const HAnd, TimeType * const Time)
{
  uint32_t i, imu, imuMax = 20;
  double mulist[imuMax], WaveFunctionAccuracy, PsiEnergy, PsiQuality, *PsiRho, GClOccupation[imuMax], PsiHist[20];
  OperatorType HAndGrd, Sx, Sy, Sz, S2;
  uint16_t *SpinUp, *SpinDn;
  WaveFunctionAccuracy = O.zeroLancWavefunction;
  mulist[0] = *mu;
  imu = 0;
  if ((O.SymOpt) & 8)
  {
    mulist[0] = Sqr(O.Nelectron) * O.U / 2;
  }
  while (imu < imuMax)
  {
    /* Create HAnd */
    BADoubleTriDiagonalMatrixToAndersonImpurityOperator(HAnd, GAndTrunk, GAndTrunk, O.U, mulist[imu]);
    /* Set HAndGrd to HAnd */
    OperatorCopy(&HAndGrd, *HAnd);
    /* if SymOpt (4) is set keep only H on the impurity and first bath site */
    if ((O.SymOpt) & 4)
    {
      printf("Restricted HAnd for the ground-state to the first bath site\n");
      HAndGrd.l[2].N = 8;
    }
    /* if singlet set is true add the (U+1) * S^2 to the operator HAndGrd */
    if ((O.SymOpt) & 1)
    {
      printf("Add %22.15E * S^2 to HAnd for the calculation of the ground-state wave-function\n", 1 + O.U);
      printf("This stabilizes the ground-state by enlarging the singlet triplet energy\n");
      SpinUp = calloc(HAnd->NFermion / 2, sizeof(uint16_t));
      SpinDn = calloc(HAnd->NFermion / 2, sizeof(uint16_t));
      SpinUp[0] = 0;
      SpinUp[1] = 2;
      SpinDn[0] = 1;
      SpinDn[1] = 3;
      for (i = 4; i < 2 + HAnd->NFermion / 2; i++)
      {
        SpinUp[i - 2] = i;
        SpinDn[i - 2] = i + HAnd->NFermion / 2 - 2;
      }
      CreateSpinOperators(HAnd->NFermion, HAnd->NFermion / 2, SpinUp, SpinDn, &Sx, &Sy, &Sz, &S2);
      FreeOperator(&Sx);
      FreeOperator(&Sy);
      FreeOperator(&Sz);
      OperatorAddNormalOrder(&HAndGrd, S2, fabs(O.U) + 1, 0, &QDetNotOrderedFcFaBcBa);
      FreeOperator(&S2);
    }
    /* create ground-state wave-function */
    BeginTiming("Lanczos GroundState", Time);
    LanczosGroundState(HAndGrd, Psi, &WaveFunctionAccuracy, O.zeroLancEnergy, O.U, Psi, 0);
    FreeOperator(&HAndGrd);
    EndTiming("Lanczos GroundState", Time);
    PsiHermitianOperatorPsi(*HAnd, *Psi, &PsiEnergy);
    PsiHermitianOperatorSqrPsi(*HAnd, *Psi, &PsiQuality);
    PsiQuality /= Sqr(PsiEnergy);
    printf("Energy                  <psi|H^2|psi>/(<psi|H|psi>^2)\n");
    printf("%22.15E  %22.15E\n", PsiEnergy, PsiQuality);
    /* Start adjusting mu */
    if ((O.SymOpt) & 8)
    {
      imu = imuMax;
    }
    else
    {
      /* calculate cluster occupation */
      PsiRho = calloc(Psi->NFermion, sizeof(double));
      if (PsiRho == NULL)
      {
        printf("calloc failed in main for PsiRho\n Fatal Error\n");
        fflush(stdout);
        return 1;
      }
      CalculateDensity(PsiRho, *Psi);
      GClOccupation[imu] = PsiRho[0] + PsiRho[1];
      free(PsiRho);
      /* Stop if occupation is close enough to occupation of G0 or if imu==imuMax-1 */
      if ((imu == imuMax - 1) || (Abs(GClOccupation[imu] - O.Nelectron) < (1000 * DBL_EPSILON))
          || (imu == 0 ? 0 : (Abs(GClOccupation[imu] - GClOccupation[imu - 1]) < (1000 * DBL_EPSILON))))
      {                         /* Stop */
        /* Warning if mu not converged when imu reached imuMax */
        if (Abs(GClOccupation[imu] - O.Nelectron) > (1000 * DBL_EPSILON))
        {
          printf("WARNING! stop changing imu without converging\n");
        }
        /* print occupation v.s. mu */
        printf("occupation v.s. mu\n");
        printf("  i    mu[i]                  Occupation             Difference\n");
        for (i = 0; i <= imu; i++)
        {
          printf("%3i %22.15E %22.15E %22.15E\n", i, mulist[i], GClOccupation[i], GClOccupation[i] - O.Nelectron);
        }
        *mu = mulist[imu];
        imu = imuMax;
      }
      else
      {                         /* make another loop */
        FreeOperator(HAnd);
        /* interpolate mu[imu+1] */
        if (imu == 0)
        {
          mulist[1] = mulist[0] - (GClOccupation[0] - O.Nelectron) * O.U;
        }
        else
        {
          if (Abs(GClOccupation[imu] - GClOccupation[imu - 1]) < (100 * DBL_EPSILON))
          {
            printf("Convergence creterium for occupation not reached, two mu's give the same occupation\n");
            mulist[imu + 1] = mulist[imu];
            imu = imuMax;
          }
          else
          {
            mulist[imu + 1] =
              mulist[imu] + (O.Nelectron - GClOccupation[imu]) * (mulist[imu] - mulist[imu - 1]) / (GClOccupation[imu] -
                                                                                                    GClOccupation[imu -
                                                                                                                  1]);
          }
        }
        printf("\nFilling should be        : %22.15E.\nThe current occupation is: %22.15E.\n", O.Nelectron,
               GClOccupation[imu]);
        printf("mu adjusted from %22.15E to %22.15E\n", mulist[imu], mulist[imu + 1]);
        printf("Restarting Lanczos Loop for the %i-th time\n", imu + 1);
        fflush(stdout);
      }
    }
    imu++;
  }


  PsiHermitianOperatorPsi(*HAnd, *Psi, &PsiEnergy);
  PsiHermitianOperatorSqrPsi(*HAnd, *Psi, &PsiQuality);
  PsiQuality /= Sqr(PsiEnergy);
  PrintWaveFunctionReadableTruncated(*Psi, 0.1);
  CalculatePrefactorHistogram(PsiHist, *Psi, 10.0, 20);
  printf("\nHistogram of prefactors in Wave-function\n");
  for (i = 0; i < 19; i++)
  {
    printf("10^(-%2i) to 10^(-%2i) | %22.15E\n", i, i + 1, PsiHist[i]);
  }
  printf("10^(-%2i) to 0        | %22.15E\n", i, PsiHist[19]);
  printf("\n");
  printf("Energy                  <psi|H^2|psi>/(<psi|H|psi>^2)\n");
  printf("%22.15E  %22.15E\n", PsiEnergy, PsiQuality);
  PsiRho = calloc(Psi->NFermion, sizeof(double));
  if (PsiRho == NULL)
  {
    printf("calloc failed in main for PsiRho\n Fatal Error\n");
    fflush(stdout);
    return 1;
  }
  CalculateDensity(PsiRho, *Psi);
  printf("Occupation \n");
  for (i = 0; i < Psi->NFermion; i++)
  {
    printf("%5i %22.15E\n", i, PsiRho[i]);
  }
  free(PsiRho);

  OperatorAddLader(HAnd, (uint16_t *) & i, 0, -PsiEnergy, 0);

  return 0;
}


/*
 *
 */
uint32_t CreateGclDMFT(DMFTSettings const O, OperatorType const HAnd, WaveFunctionType const Psi,
                       OperatorType const *VPerturbation, TriDiagonalMatrixType const *gPerturbation, double const mu,
                       AndersonMatrixType * const GclAnd, TimeType * const Time)
{
  uint32_t i;
  double Egrd;
  OperatorType OPESdn, OPESup, OIPESdn, OIPESup;
  WaveFunctionType OPsi;
  TriDiagonalMatrixType GclMup, GclMdn, GclPup, GclPdn;
  ListOfPoles GclMupLP, GclMdnLP, GclPupLP, GclPdnLP, GclLP;
  /* Get ground-state energy needed for shift */
  PsiHermitianOperatorPsi(HAnd, Psi, &Egrd);
  /* create operators */
  CreatePESdnOperatorDMFT(&OPESdn, HAnd.NFermion);
  CreatePESupOperatorDMFT(&OPESup, HAnd.NFermion);
  CreateIPESdnOperatorDMFT(&OIPESdn, HAnd.NFermion);
  CreateIPESupOperatorDMFT(&OIPESup, HAnd.NFermion);
  /* do the rest - see timing for comments */
  BeginTiming("Lanczos Tri Diagonalization", Time);
  OperatorPsi(OPESup, Psi, &OPsi);
  switch (O.Dyson)
  {
  case 1:
    GclMup.N = Max(O.NTriDiag, O.NGCluster);
    InitTriDiagonalMatrix(&GclMup);
    LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, O.zeroLancWavefunction, O.NTriDiag, VPerturbation, gPerturbation, 4,
                                     &GclMup);
    break;
  default:
    GclMup.N = O.NTriDiag;
    InitTriDiagonalMatrix(&GclMup);
    LanczosTriDiagonalize(HAnd, OPsi, O.zeroLancWavefunction, &GclMup, O.Gamma);
    break;
  }
  EndTiming("Lanczos Tri Diagonalization", Time);
  BeginTiming("Lanczos Tri Diagonalization", Time);
  OperatorPsi(OPESdn, Psi, &OPsi);
  switch (O.Dyson)
  {
  case 1:
    GclMdn.N = Max(O.NTriDiag, O.NGCluster);
    InitTriDiagonalMatrix(&GclMdn);
    LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, O.zeroLancWavefunction, O.NTriDiag, VPerturbation, gPerturbation, 4, &GclMdn);
    break;
  default:
    GclMdn.N = O.NTriDiag;
    InitTriDiagonalMatrix(&GclMdn);
    LanczosTriDiagonalize(HAnd, OPsi, O.zeroLancWavefunction, &GclMdn, O.Gamma);
    break;
  }
  EndTiming("Lanczos Tri Diagonalization", Time);
  BeginTiming("Lanczos Tri Diagonalization", Time);
  OperatorPsi(OIPESup, Psi, &OPsi);
  switch (O.Dyson)
  {
  case 1:
    GclPup.N = Max(O.NTriDiag, O.NGCluster);
    InitTriDiagonalMatrix(&GclPup);
    LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, O.zeroLancWavefunction, O.NTriDiag, VPerturbation, gPerturbation, 4, &GclPup);
    break;
  default:
    GclPup.N = O.NTriDiag;
    InitTriDiagonalMatrix(&GclPup);
    LanczosTriDiagonalize(HAnd, OPsi, O.zeroLancWavefunction, &GclPup, O.Gamma);
    break;
  }
  EndTiming("Lanczos Tri Diagonalization", Time);
  BeginTiming("Lanczos Tri Diagonalization", Time);
  OperatorPsi(OIPESdn, Psi, &OPsi);
  switch (O.Dyson)
  {
  case 1:
    GclPdn.N = Max(O.NTriDiag, O.NGCluster);
    InitTriDiagonalMatrix(&GclPdn);
    LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, O.zeroLancWavefunction, O.NTriDiag, VPerturbation, gPerturbation, 4,
                                     &GclPdn);
    break;
  default:
    GclPdn.N = O.NTriDiag;
    InitTriDiagonalMatrix(&GclPdn);
    LanczosTriDiagonalize(HAnd, OPsi, O.zeroLancWavefunction, &GclPdn, O.Gamma);
    break;
  }
  EndTiming("Lanczos Tri Diagonalization", Time);
  ShiftTriDiagonalMatrix(&GclMdn, -Egrd);
  ShiftTriDiagonalMatrix(&GclMup, -Egrd);
  ShiftTriDiagonalMatrix(&GclPdn, -Egrd);
  ShiftTriDiagonalMatrix(&GclPup, -Egrd);
  GclMdn.mu = mu;
  GclMup.mu = mu;
  GclPdn.mu = mu;
  GclPup.mu = mu;
  for (i = 1; i <= GclMdn.N; i++)
  {
    GclMdn.a[i] *= -1;
  }
  for (i = 1; i <= GclMup.N; i++)
  {
    GclMup.a[i] *= -1;
  }
  strcpy(GclPup.Name, "Gcluster-Plus-Up");
  strcpy(GclPdn.Name, "Gcluster-Plus-Dn");
  strcpy(GclMup.Name, "Gcluster-Min-Up");
  strcpy(GclMdn.Name, "Gcluster-Min-Dn");
  printf("Start SumGreensfunctionsTriDiagonalRepresentation\n");
  fflush(stdout);
  TriDiagonalMatrixToListOfPoles(GclMup, &GclMupLP);
  TriDiagonalMatrixToListOfPoles(GclMdn, &GclMdnLP);
  TriDiagonalMatrixToListOfPoles(GclPup, &GclPupLP);
  TriDiagonalMatrixToListOfPoles(GclPdn, &GclPdnLP);


  if (O.NmuSym & 2)
  {
    /* calculations at constant mu, move spectra at negative energy to zero in order to excellerate converges */
    ListOfPolesMoveNegativeEnergyToMu(&GclPupLP, mu);
    ListOfPolesMoveNegativeEnergyToMu(&GclPdnLP, mu);
    ListOfPolesMovePositiveEnergyToMu(&GclMupLP, mu);
    ListOfPolesMovePositiveEnergyToMu(&GclMdnLP, mu);
  }


  ListOfPolesCopy(&GclLP, GclMupLP);
  AddListOfPoles(&GclLP, GclMdnLP);
  AddListOfPoles(&GclLP, GclPupLP);
  AddListOfPoles(&GclLP, GclPdnLP);
  ScaleListOfPoles(&GclLP, 0.5);
  strcpy(GclLP.Name, "Gcluster");
  GclLP.mu = mu;
  /* symmetrize */
  if ((O.SymOpt) & 2)
  {
    printf("Symmetryze %s around mu = %22.15E\n", GclLP.Name, GclLP.mu);
    SymmetrizeListOfPoles(&GclLP);
  }
  printf("Cluster Green function created of size %i ", GclLP.N);
  ReduceListOfPolesNumberFermi(&GclLP, O.NGCluster, mu);
  printf("reduced to %i\n", GclLP.N);
  GclLP.mu = mu;
  ListOfPolesToAndersonMatrix(GclLP, GclAnd);
  strcpy(GclAnd->Name, "Gcluster");
  GclAnd->b[0] = 1.0;
  /* print occupations from spectra */
  PrintListOfPolesEssentials(GclMdnLP);
  PrintListOfPolesEssentials(GclPdnLP);
  PrintListOfPolesEssentials(GclMupLP);
  PrintListOfPolesEssentials(GclPupLP);
  PrintListOfPolesEssentials(GclLP);
  PrintAndersonMatrixEssentials(*GclAnd);
  /* free operators and spectra */
  FreeListOfPoles(&GclLP);
  FreeListOfPoles(&GclMupLP);
  FreeListOfPoles(&GclMdnLP);
  FreeListOfPoles(&GclPupLP);
  FreeListOfPoles(&GclPdnLP);
  FreeTriDiagonalMatrix(&GclMup);
  FreeTriDiagonalMatrix(&GclMdn);
  FreeTriDiagonalMatrix(&GclPup);
  FreeTriDiagonalMatrix(&GclPdn);
  FreeOperator(&OPESdn);
  FreeOperator(&OPESup);
  FreeOperator(&OIPESdn);
  FreeOperator(&OIPESup);
  return 0;
}
