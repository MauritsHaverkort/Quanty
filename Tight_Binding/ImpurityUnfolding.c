


int TightBindingImpurityUnfolding1Impurity3D(int Nsupercell, double Vpot, char *filename)
{
  TightBindingDef Set, SetSuper;
  int    i, N[3];
  int    Ncharacter;
  char   CharFileName[255];
  double Emin, Emax, ELines[6];
  int    NEpoints, NELines;
  MatrixType BWValue, RGBValue[3], AllValue;
  MatrixType *Character;
  BandstructurePathType QPath;
  MatrixType BANDS;
  CreateOneBandCubicLaticeTightBindingDef(&Set);
//  PrintTightBindingDef(Set);
  N[0] = Nsupercell;
  N[1] = Nsupercell;
  N[2] = Nsupercell;
  CreateTightBindingSquareSuperCell(Set, &SetSuper, N);
  SetSuper.H[0].Hop.VRe[0][0] = Vpot;
//  PrintTightBindingDef(SetSuper);
  QPath.NQ = 5;
  QPath.Q = malloc(QPath.NQ * sizeof(BandstructurePointType));
  QPath.Q[0].Nq = 142;
  QPath.Q[0].label = 'G';
  QPath.Q[0].q[0] = 0.0 * Pi;
  QPath.Q[0].q[1] = 0.0 * Pi;
  QPath.Q[0].q[2] = 0.0 * Pi;
  QPath.Q[1].Nq = 100;
  QPath.Q[1].label = 'M';
  QPath.Q[1].q[0] = 1.0 * Pi;
  QPath.Q[1].q[1] = 1.0 * Pi;
  QPath.Q[1].q[2] = 0.0 * Pi;
  QPath.Q[2].Nq = 100;
  QPath.Q[2].label = 'X';
  QPath.Q[2].q[0] = 1.0 * Pi;
  QPath.Q[2].q[1] = 0.0 * Pi;
  QPath.Q[2].q[2] = 0.0 * Pi;
  QPath.Q[3].Nq = 174;
  QPath.Q[3].label = 'G';
  QPath.Q[3].q[0] = 0.0 * Pi;
  QPath.Q[3].q[1] = 0.0 * Pi;
  QPath.Q[3].q[2] = 0.0 * Pi;
  QPath.Q[4].Nq = 1;
  QPath.Q[4].label = 'R';
  QPath.Q[4].q[0] = 1.0 * Pi;
  QPath.Q[4].q[1] = 1.0 * Pi;
  QPath.Q[4].q[2] = 1.0 * Pi;

  EigenFunctionToCharactersOnsiteEnergyPhaseSetNCharacter(SetSuper, &Ncharacter);
  Character = malloc(Ncharacter * sizeof(MatrixType));
  TightBindingBANDSCharacter(SetSuper, QPath, &BANDS, Ncharacter, Character,
                             &EigenFunctionToCharactersOnsiteEnergyPhase);
  sprintf(CharFileName, "%s.energy", filename);
  FilePrintPlainMatrix(BANDS, CharFileName, "w");
  for (i = 0; i < Ncharacter; i++)
  {
    if (i < 10)
      sprintf(CharFileName, "%s.character_00%1i", filename, i);
    else if (i < 100)
      sprintf(CharFileName, "%s.character_0%2i", filename, i);
    else if (i < 1000)
      sprintf(CharFileName, "%s.character_%3i", filename, i);
    else
      sprintf(CharFileName, "%s.character_%i", filename, i);
    FilePrintPlainMatrix(Character[i], CharFileName, "w");
  }

  Emin = -6.3;
  Emax = 9.3;
  NEpoints = 516;
  NELines = 6;
  ELines[0] = -6;
  ELines[1] = -3;
  ELines[2] = 0;
  ELines[3] = 3;
  ELines[4] = 6;
  ELines[5] = 9;

  TightBindingBANDSToMatrix(Emin, Emax, NEpoints, BANDS, &AllValue);
  TightBindingAddAxesToMatrix(Emin, Emax, NELines, ELines, QPath, &AllValue);
  sprintf(CharFileName, "%s_all.bmp", filename);
  FilePrintRealMatrixToBitmap(AllValue, CharFileName, &RealToRGBColorBWInverse);

  TightBindingBANDSCharacterToMatrix(Emin, Emax, NEpoints, BANDS, Character[2], &BWValue);
  TightBindingAddAxesToMatrix(Emin, Emax, NELines, ELines, QPath, &BWValue);
  sprintf(CharFileName, "%s_BW.bmp", filename);
  FilePrintRealMatrixToBitmap(BWValue, CharFileName, &RealToRGBColorBWInverse);

  TightBindingBANDSCharacterToMatrix(Emin, Emax, NEpoints, BANDS, Character[0], &RGBValue[0]);
  TightBindingBANDSCharacterToMatrix(Emin, Emax, NEpoints, BANDS, Character[2], &RGBValue[1]);
  TightBindingBANDSCharacterToMatrix(Emin, Emax, NEpoints, BANDS, Character[1], &RGBValue[2]);
  TightBindingAddAxesToMatrix(Emin, Emax, NELines, ELines, QPath, &RGBValue[0]);
  TightBindingAddAxesToMatrix(Emin, Emax, NELines, ELines, QPath, &RGBValue[1]);
  TightBindingAddAxesToMatrix(Emin, Emax, NELines, ELines, QPath, &RGBValue[2]);
  sprintf(CharFileName, "%s_RGB.bmp", filename);
  FilePrintReal3MatrixToBitmap(RGBValue, CharFileName, &Real3ToRGBColorSimpleInverse);

  return 0;
}

int TightBindingImpurityUnfolding1Impurity2D(int Nsupercell, double Vpot, char *filename)
{
  TightBindingDef Set, SetSuper;
  int    i, N[3];
  int    Ncharacter;
  char   CharFileName[255];
  double Emin, Emax, ELines[6];
  int    NEpoints, NELines;
  MatrixType BWValue, RGBValue[3], AllValue;
  MatrixType *Character;
  BandstructurePathType QPath;
  MatrixType BANDS;
  CreateOneBandSquareLaticeTightBindingDef(&Set);
  PrintTightBindingDef(Set);
  N[0] = Nsupercell;
  N[1] = Nsupercell;
  N[2] = 1;
  CreateTightBindingSquareSuperCell(Set, &SetSuper, N);
  SetSuper.H[0].Hop.VRe[0][0] = Vpot;
  PrintTightBindingDef(SetSuper);
  QPath.NQ = 4;
  QPath.Q = malloc(QPath.NQ * sizeof(BandstructurePointType));
  QPath.Q[0].Nq = 213;
  QPath.Q[0].label = 'G';
  QPath.Q[0].q[0] = 0.0 * Pi;
  QPath.Q[0].q[1] = 0.0 * Pi;
  QPath.Q[0].q[2] = 0.0 * Pi;
  QPath.Q[1].Nq = 150;
  QPath.Q[1].label = 'M';
  QPath.Q[1].q[0] = 1.0 * Pi;
  QPath.Q[1].q[1] = 1.0 * Pi;
  QPath.Q[1].q[2] = 0.0 * Pi;
  QPath.Q[2].Nq = 150;
  QPath.Q[2].label = 'X';
  QPath.Q[2].q[0] = 1.0 * Pi;
  QPath.Q[2].q[1] = 0.0 * Pi;
  QPath.Q[2].q[2] = 0.0 * Pi;
  QPath.Q[3].Nq = 1;
  QPath.Q[3].label = 'G';
  QPath.Q[3].q[0] = 0.0 * Pi;
  QPath.Q[3].q[1] = 0.0 * Pi;
  QPath.Q[3].q[2] = 0.0 * Pi;

  EigenFunctionToCharactersOnsiteEnergyPhaseSetNCharacter(SetSuper, &Ncharacter);
  Character = malloc(Ncharacter * sizeof(MatrixType));
  TightBindingBANDSCharacter(SetSuper, QPath, &BANDS, Ncharacter, Character,
                             &EigenFunctionToCharactersOnsiteEnergyPhase);
  sprintf(CharFileName, "%s.energy", filename);
  FilePrintPlainMatrix(BANDS, CharFileName, "w");
  for (i = 0; i < Ncharacter; i++)
  {
    if (i < 10)
      sprintf(CharFileName, "%s.character_00%1i", filename, i);
    else if (i < 100)
      sprintf(CharFileName, "%s.character_0%2i", filename, i);
    else if (i < 1000)
      sprintf(CharFileName, "%s.character_%3i", filename, i);
    else
      sprintf(CharFileName, "%s.character_%i", filename, i);
    FilePrintPlainMatrix(Character[i], CharFileName, "w");
  }

  Emin = -4.2;
  Emax = 6.2;
  NEpoints = 513;
  NELines = 6;
  ELines[0] = -4;
  ELines[1] = -2;
  ELines[2] = 0;
  ELines[3] = 2;
  ELines[4] = 4;
  ELines[5] = 6;

  TightBindingBANDSToMatrix(Emin, Emax, NEpoints, BANDS, &AllValue);
  TightBindingAddAxesToMatrix(Emin, Emax, NELines, ELines, QPath, &AllValue);
  sprintf(CharFileName, "%s_all.bmp", filename);
  FilePrintRealMatrixToBitmap(AllValue, CharFileName, &RealToRGBColorBWInverse);

  TightBindingBANDSCharacterToMatrix(Emin, Emax, NEpoints, BANDS, Character[2], &BWValue);
  TightBindingAddAxesToMatrix(Emin, Emax, NELines, ELines, QPath, &BWValue);
  sprintf(CharFileName, "%s_BW.bmp", filename);
  FilePrintRealMatrixToBitmap(BWValue, CharFileName, &RealToRGBColorBWInverse);

  TightBindingBANDSCharacterToMatrix(Emin, Emax, NEpoints, BANDS, Character[0], &RGBValue[0]);
  TightBindingBANDSCharacterToMatrix(Emin, Emax, NEpoints, BANDS, Character[2], &RGBValue[1]);
  TightBindingBANDSCharacterToMatrix(Emin, Emax, NEpoints, BANDS, Character[1], &RGBValue[2]);
  TightBindingAddAxesToMatrix(Emin, Emax, NELines, ELines, QPath, &RGBValue[0]);
  TightBindingAddAxesToMatrix(Emin, Emax, NELines, ELines, QPath, &RGBValue[1]);
  TightBindingAddAxesToMatrix(Emin, Emax, NELines, ELines, QPath, &RGBValue[2]);
  sprintf(CharFileName, "%s_RGB.bmp", filename);
  FilePrintReal3MatrixToBitmap(RGBValue, CharFileName, &Real3ToRGBColorSimpleInverse);

  return 0;
}



int TightBindingImpurityUnfolding1Impurity1D(int Nsupercell, double Vpot, char *filename)
{
  TightBindingDef Set, SetSuper;
  int    i, N[3];
  int    Ncharacter;
  char   CharFileName[255];
  double Emin, Emax, ELines[6];
  int    NEpoints, NELines;
  MatrixType BWValue, RGBValue[3], AllValue;
  MatrixType *Character;
  BandstructurePathType QPath;
  MatrixType BANDS;
  CreateOneBandChainLaticeTightBindingDef(&Set);
  PrintTightBindingDef(Set);
  N[0] = Nsupercell;
  N[1] = 1;
  N[2] = 1;
  CreateTightBindingSquareSuperCell(Set, &SetSuper, N);
  SetSuper.H[0].Hop.VRe[0][0] = Vpot;
  PrintTightBindingDef(SetSuper);
  QPath.NQ = 2;
  QPath.Q = malloc(QPath.NQ * sizeof(BandstructurePointType));
  QPath.Q[0].Nq = 500;
  QPath.Q[0].label = 'G';
  QPath.Q[0].q[0] = 0.0 * Pi;
  QPath.Q[0].q[1] = 0.0 * Pi;
  QPath.Q[0].q[2] = 0.0 * Pi;
  QPath.Q[1].Nq = 1;
  QPath.Q[1].label = 'X';
  QPath.Q[1].q[0] = 1.0 * Pi;
  QPath.Q[1].q[1] = 0.0 * Pi;
  QPath.Q[1].q[2] = 0.0 * Pi;

  EigenFunctionToCharactersOnsiteEnergyPhaseSetNCharacter(SetSuper, &Ncharacter);
  Character = malloc(Ncharacter * sizeof(MatrixType));
  TightBindingBANDSCharacter(SetSuper, QPath, &BANDS, Ncharacter, Character,
                             &EigenFunctionToCharactersOnsiteEnergyPhase);
  sprintf(CharFileName, "%s.energy", filename);
  FilePrintPlainMatrix(BANDS, CharFileName, "w");
  for (i = 0; i < Ncharacter; i++)
  {
    if (i < 10)
      sprintf(CharFileName, "%s.character_00%1i", filename, i);
    else if (i < 100)
      sprintf(CharFileName, "%s.character_0%2i", filename, i);
    else if (i < 1000)
      sprintf(CharFileName, "%s.character_%3i", filename, i);
    else
      sprintf(CharFileName, "%s.character_%i", filename, i);
    FilePrintPlainMatrix(Character[i], CharFileName, "w");
  }

  Emin = -2.1;
  Emax = 3.1;
  NEpoints = 500;
  NELines = 6;
  ELines[0] = -2;
  ELines[1] = -1;
  ELines[2] = 0;
  ELines[3] = 1;
  ELines[4] = 2;
  ELines[5] = 3;

  TightBindingBANDSToMatrix(Emin, Emax, NEpoints, BANDS, &AllValue);
  TightBindingAddAxesToMatrix(Emin, Emax, NELines, ELines, QPath, &AllValue);
  sprintf(CharFileName, "%s_all.bmp", filename);
  FilePrintRealMatrixToBitmap(AllValue, CharFileName, &RealToRGBColorBWInverse);

  TightBindingBANDSCharacterToMatrix(Emin, Emax, NEpoints, BANDS, Character[2], &BWValue);
  TightBindingAddAxesToMatrix(Emin, Emax, NELines, ELines, QPath, &BWValue);
  sprintf(CharFileName, "%s_BW.bmp", filename);
  FilePrintRealMatrixToBitmap(BWValue, CharFileName, &RealToRGBColorBWInverse);

  TightBindingBANDSCharacterToMatrix(Emin, Emax, NEpoints, BANDS, Character[0], &RGBValue[0]);
  TightBindingBANDSCharacterToMatrix(Emin, Emax, NEpoints, BANDS, Character[2], &RGBValue[1]);
  TightBindingBANDSCharacterToMatrix(Emin, Emax, NEpoints, BANDS, Character[1], &RGBValue[2]);
  TightBindingAddAxesToMatrix(Emin, Emax, NELines, ELines, QPath, &RGBValue[0]);
  TightBindingAddAxesToMatrix(Emin, Emax, NELines, ELines, QPath, &RGBValue[1]);
  TightBindingAddAxesToMatrix(Emin, Emax, NELines, ELines, QPath, &RGBValue[2]);
  sprintf(CharFileName, "%s_RGB.bmp", filename);
  FilePrintReal3MatrixToBitmap(RGBValue, CharFileName, &Real3ToRGBColorSimpleInverse);

  return 0;
}


int TightBindingImpurityUnfoldingRandomImpurity3D(int Nsupercell, double Rho, double Vpot, char *filename,
                                                  MatrixType * BWValue, MatrixType * RGBValue, SpectraType DOS[],
                                                  MatrixType Fermi[])
{
  TightBindingDef Set, SetSuper;
  int    i, N[3], NkDOS[3];
  int    Ncharacter;
  int    QImpurity, NImpurity;
  char   CharFileName[255];
  double Emin, Emax, ELines[10];
  double k0[3], k1[3], k2[3];
  int    NEpoints, NELines;
  MatrixType AllValue;
  MatrixType *Character;
  BandstructurePathType QPath;
  MatrixType BANDS;
  CreateOneBandCubicLaticeTightBindingDef(&Set);
//  PrintTightBindingDef(Set);
  N[0] = Nsupercell;
  N[1] = Nsupercell;
  N[2] = Nsupercell;
  CreateTightBindingSquareSuperCell(Set, &SetSuper, N);
  NImpurity = 0;
  for (i = 0; i < Cubed(Nsupercell); i++)
  {
    QImpurity = 0;
    if (Rho > mwhdrand())
      QImpurity = 1;
    SetSuper.H[7 * i].Hop.VRe[0][0] = Vpot * QImpurity;
    printf("%1i", QImpurity);
    NImpurity += QImpurity;
  }
  printf("\n");
  printf("Number of Impurities: %i\n", NImpurity);
  printf("Number of Sites     : %i\n", Cubed(Nsupercell));
  printf("Rho                 : %22.15E\n", ((double)NImpurity) / ((double)Cubed(Nsupercell)));
  printf("Rho average (input) : %22.15E\n", Rho);
//  PrintTightBindingDef(SetSuper);
  QPath.NQ = 5;
  QPath.Q = malloc(QPath.NQ * sizeof(BandstructurePointType));
  QPath.Q[0].Nq = 284;
  QPath.Q[0].label = 'G';
  QPath.Q[0].q[0] = 0.0 * Pi;
  QPath.Q[0].q[1] = 0.0 * Pi;
  QPath.Q[0].q[2] = 0.0 * Pi;
  QPath.Q[1].Nq = 200;
  QPath.Q[1].label = 'M';
  QPath.Q[1].q[0] = 1.0 * Pi;
  QPath.Q[1].q[1] = 1.0 * Pi;
  QPath.Q[1].q[2] = 0.0 * Pi;
  QPath.Q[2].Nq = 200;
  QPath.Q[2].label = 'X';
  QPath.Q[2].q[0] = 1.0 * Pi;
  QPath.Q[2].q[1] = 0.0 * Pi;
  QPath.Q[2].q[2] = 0.0 * Pi;
  QPath.Q[3].Nq = 348;
  QPath.Q[3].label = 'G';
  QPath.Q[3].q[0] = 0.0 * Pi;
  QPath.Q[3].q[1] = 0.0 * Pi;
  QPath.Q[3].q[2] = 0.0 * Pi;
  QPath.Q[4].Nq = 1;
  QPath.Q[4].label = 'R';
  QPath.Q[4].q[0] = 1.0 * Pi;
  QPath.Q[4].q[1] = 1.0 * Pi;
  QPath.Q[4].q[2] = 1.0 * Pi;

  EigenFunctionToCharactersOnsiteEnergyPhaseSetNCharacter(SetSuper, &Ncharacter);
  Character = malloc(Ncharacter * sizeof(MatrixType));
  TightBindingBANDSCharacter(SetSuper, QPath, &BANDS, Ncharacter, Character,
                             &EigenFunctionToCharactersOnsiteEnergyPhase);
  sprintf(CharFileName, "%s.energy", filename);
  FilePrintPlainMatrix(BANDS, CharFileName, "w");
  for (i = 0; i < Ncharacter; i++)
  {
    if (i < 10)
      sprintf(CharFileName, "%s.character_00%1i", filename, i);
    else if (i < 100)
      sprintf(CharFileName, "%s.character_0%2i", filename, i);
    else if (i < 1000)
      sprintf(CharFileName, "%s.character_%3i", filename, i);
    else
      sprintf(CharFileName, "%s.character_%i", filename, i);
    FilePrintPlainMatrix(Character[i], CharFileName, "w");
  }

  Emin = -1.532;
  Emax = 0.532 + 0.250;
  NEpoints = 1032 + 125;
  NELines = 10;
  ELines[0] = -1.5;
  ELines[1] = -1.25;
  ELines[2] = -1.0;
  ELines[3] = -0.75;
  ELines[4] = -0.5;
  ELines[5] = -0.25;
  ELines[6] = 0.0;
  ELines[7] = 0.25;
  ELines[8] = 0.5;
  ELines[9] = 0.75;

  TightBindingBANDSToMatrix(Emin, Emax, NEpoints, BANDS, &AllValue);
  TightBindingAddAxesToMatrix(Emin, Emax, NELines, ELines, QPath, &AllValue);
  sprintf(CharFileName, "%s_all.bmp", filename);
  FilePrintRealMatrixToBitmap(AllValue, CharFileName, &RealToRGBColorBWInverse);

  TightBindingBANDSCharacterToMatrix(Emin, Emax, NEpoints, BANDS, Character[2], BWValue);
//  TightBindingAddAxesToMatrix(Emin,Emax,NELines,ELines,QPath, BWValue);
  sprintf(CharFileName, "%s_BW.bmp", filename);
  FilePrintRealMatrixToBitmap(*BWValue, CharFileName, &RealToRGBColorBWInverse);

  TightBindingBANDSCharacterToMatrix(Emin, Emax, NEpoints, BANDS, Character[0], &RGBValue[0]);
  TightBindingBANDSCharacterToMatrix(Emin, Emax, NEpoints, BANDS, Character[2], &RGBValue[1]);
  TightBindingBANDSCharacterToMatrix(Emin, Emax, NEpoints, BANDS, Character[1], &RGBValue[2]);
  TightBindingAddAxesToMatrix(Emin, Emax, NELines, ELines, QPath, &RGBValue[0]);
  TightBindingAddAxesToMatrix(Emin, Emax, NELines, ELines, QPath, &RGBValue[1]);
  TightBindingAddAxesToMatrix(Emin, Emax, NELines, ELines, QPath, &RGBValue[2]);
  sprintf(CharFileName, "%s_RGB.bmp", filename);
  FilePrintReal3MatrixToBitmap(RGBValue, CharFileName, &Real3ToRGBColorSimpleInverse);

  for (i = 0; i < Ncharacter; i++)
  {
    FreeMatrix(&Character[i]);
  }
  free(Character);

  NkDOS[0] = (30 * mwhdrand() + 200) / Nsupercell;
  NkDOS[1] = NkDOS[0] + 1;
  NkDOS[2] = NkDOS[0] + 2;
  DOS[0].Emin = Emin;
  DOS[0].Emax = Emax;
  DOS[0].N = NEpoints;
  DOS[0].dE = (Emax - Emin) / ((double)(NEpoints));
  SpectrumAllocate(DOS[0], &(DOS[0]));
  SpectrumAllocate(DOS[0], &(DOS[1]));
  SpectrumAllocate(DOS[0], &(DOS[2]));
  EigenFunctionToCharactersOnsiteEnergySetNCharacter(SetSuper, &Ncharacter);
  TightBindingDOSCharacter(SetSuper, NkDOS, DOS, Ncharacter, &EigenFunctionToCharactersOnsiteEnergy);
  sprintf(CharFileName, "%s.DOS", filename);
  FilePrintSpectra(DOS, Ncharacter, CharFileName, "w");
  k0[0] = 0;
  k0[1] = 0;
  k0[2] = 0;
  k1[0] = 2 * Pi / Nsupercell;
  k1[1] = 0;
  k1[2] = 0;
  k2[0] = 0;
  k2[1] = 2 * Pi / Nsupercell;
  k2[2] = 0;
  EigenFunctionToCharactersOnsiteEnergyPhaseRGBSetNCharacter(SetSuper, &Ncharacter);
  TightBindingFermiCharacterUnfolding(SetSuper, Emin, Emax, NEpoints, k0, k1, k2, Nsupercell, Nsupercell, Fermi,
                                      Ncharacter, &EigenFunctionToCharactersOnsiteEnergyPhaseRGBkshifted);
  return 0;
}


int TightBindingImpurityUnfoldingRandomImpurity2D(int Nsupercell, double Rho, double Vpot, char *filename,
                                                  MatrixType * BWValue, MatrixType * RGBValue, SpectraType DOS[],
                                                  MatrixType Fermi[])
{
  TightBindingDef Set, SetSuper;
  int    i, N[3], NkDOS[3];
  int    Ncharacter;
  int    QImpurity, NImpurity;
  char   CharFileName[255];
  double Emin, Emax, ELines[10];
  double k0[3], k1[3], k2[3];
  int    NEpoints, NELines;
  MatrixType AllValue;
  MatrixType *Character;
  BandstructurePathType QPath;
  MatrixType BANDS;
  CreateOneBandSquareLaticeTightBindingDef(&Set);
//  PrintTightBindingDef(Set);
  N[0] = Nsupercell;
  N[1] = Nsupercell;
  N[2] = 1;
  CreateTightBindingSquareSuperCell(Set, &SetSuper, N);
  NImpurity = 0;
  for (i = 0; i < Sqr(Nsupercell); i++)
  {
    QImpurity = 0;
    if (Rho > mwhdrand())
      QImpurity = 1;
    SetSuper.H[5 * i].Hop.VRe[0][0] = Vpot * QImpurity;
    printf("%1i", QImpurity);
    NImpurity += QImpurity;
  }
  printf("\n");
  printf("Number of Impurities: %i\n", NImpurity);
  printf("Number of Sites     : %i\n", Sqr(Nsupercell));
  printf("Rho                 : %22.15E\n", ((double)NImpurity) / ((double)Sqr(Nsupercell)));
  printf("Rho average (input) : %22.15E\n", Rho);
//  PrintTightBindingDef(SetSuper);
  QPath.NQ = 4;
  QPath.Q = malloc(QPath.NQ * sizeof(BandstructurePointType));
  QPath.Q[0].Nq = 284;
  QPath.Q[0].label = 'G';
  QPath.Q[0].q[0] = 0.0 * Pi;
  QPath.Q[0].q[1] = 0.0 * Pi;
  QPath.Q[0].q[2] = 0.0 * Pi;
  QPath.Q[1].Nq = 200;
  QPath.Q[1].label = 'M';
  QPath.Q[1].q[0] = 1.0 * Pi;
  QPath.Q[1].q[1] = 1.0 * Pi;
  QPath.Q[1].q[2] = 0.0 * Pi;
  QPath.Q[2].Nq = 200;
  QPath.Q[2].label = 'X';
  QPath.Q[2].q[0] = 1.0 * Pi;
  QPath.Q[2].q[1] = 0.0 * Pi;
  QPath.Q[2].q[2] = 0.0 * Pi;
  QPath.Q[3].Nq = 1;
  QPath.Q[3].label = 'G';
  QPath.Q[3].q[0] = 0.0 * Pi;
  QPath.Q[3].q[1] = 0.0 * Pi;
  QPath.Q[3].q[2] = 0.0 * Pi;

  EigenFunctionToCharactersOnsiteEnergyPhaseSetNCharacter(SetSuper, &Ncharacter);
  Character = malloc(Ncharacter * sizeof(MatrixType));
  TightBindingBANDSCharacter(SetSuper, QPath, &BANDS, Ncharacter, Character,
                             &EigenFunctionToCharactersOnsiteEnergyPhase);
  sprintf(CharFileName, "%s.energy", filename);
  FilePrintPlainMatrix(BANDS, CharFileName, "w");
  for (i = 0; i < Ncharacter; i++)
  {
    if (i < 10)
      sprintf(CharFileName, "%s.character_00%1i", filename, i);
    else if (i < 100)
      sprintf(CharFileName, "%s.character_0%2i", filename, i);
    else if (i < 1000)
      sprintf(CharFileName, "%s.character_%3i", filename, i);
    else
      sprintf(CharFileName, "%s.character_%i", filename, i);
    FilePrintPlainMatrix(Character[i], CharFileName, "w");
  }

  Emin = -1.532;
  Emax = 0.532 + 0.250;
  NEpoints = 1032 + 125;
  NELines = 10;
  ELines[0] = -1.5;
  ELines[1] = -1.25;
  ELines[2] = -1.0;
  ELines[3] = -0.75;
  ELines[4] = -0.5;
  ELines[5] = -0.25;
  ELines[6] = 0;
  ELines[7] = 0.25;
  ELines[8] = 0.50;
  ELines[9] = 0.75;

  TightBindingBANDSToMatrix(Emin, Emax, NEpoints, BANDS, &AllValue);
  TightBindingAddAxesToMatrix(Emin, Emax, NELines, ELines, QPath, &AllValue);
  sprintf(CharFileName, "%s_all.bmp", filename);
  FilePrintRealMatrixToBitmap(AllValue, CharFileName, &RealToRGBColorBWInverse);

  FreeMatrix(&AllValue);

  TightBindingBANDSCharacterToMatrix(Emin, Emax, NEpoints, BANDS, Character[2], BWValue);
//  TightBindingAddAxesToMatrix(Emin,Emax,NELines,ELines,QPath, BWValue);
  sprintf(CharFileName, "%s_BW.bmp", filename);
  FilePrintRealMatrixToBitmap(*BWValue, CharFileName, &RealToRGBColorBWInverse);

  TightBindingBANDSCharacterToMatrix(Emin, Emax, NEpoints, BANDS, Character[0], &RGBValue[0]);
  TightBindingBANDSCharacterToMatrix(Emin, Emax, NEpoints, BANDS, Character[2], &RGBValue[1]);
  TightBindingBANDSCharacterToMatrix(Emin, Emax, NEpoints, BANDS, Character[1], &RGBValue[2]);
  TightBindingAddAxesToMatrix(Emin, Emax, NELines, ELines, QPath, &RGBValue[0]);
  TightBindingAddAxesToMatrix(Emin, Emax, NELines, ELines, QPath, &RGBValue[1]);
  TightBindingAddAxesToMatrix(Emin, Emax, NELines, ELines, QPath, &RGBValue[2]);
  sprintf(CharFileName, "%s_RGB.bmp", filename);
  FilePrintReal3MatrixToBitmap(RGBValue, CharFileName, &Real3ToRGBColorSimpleInverse);

  for (i = 0; i < Ncharacter; i++)
  {
    FreeMatrix(&Character[i]);
  }
  free(Character);
  FreeMatrix(&BANDS);

  NkDOS[0] = (30 * mwhdrand() + 200) / Nsupercell;
  NkDOS[1] = NkDOS[0] + 1;
  NkDOS[2] = 1;
  DOS[0].Emin = Emin;
  DOS[0].Emax = Emax;
  DOS[0].N = NEpoints;
  DOS[0].dE = (Emax - Emin) / ((double)(NEpoints));
  SpectrumAllocate(DOS[0], &(DOS[0]));
  SpectrumAllocate(DOS[0], &(DOS[1]));
  SpectrumAllocate(DOS[0], &(DOS[2]));
  EigenFunctionToCharactersOnsiteEnergySetNCharacter(SetSuper, &Ncharacter);
  TightBindingDOSCharacter(SetSuper, NkDOS, DOS, Ncharacter, &EigenFunctionToCharactersOnsiteEnergy);
  sprintf(CharFileName, "%s.DOS", filename);
  FilePrintSpectra(DOS, Ncharacter, CharFileName, "w");
  k0[0] = 0;
  k0[1] = 0;
  k0[2] = 0;
  k1[0] = 2 * Pi / Nsupercell;
  k1[1] = 0;
  k1[2] = 0;
  k2[0] = 0;
  k2[1] = 2 * Pi / Nsupercell;
  k2[2] = 0;
  EigenFunctionToCharactersOnsiteEnergyPhaseRGBSetNCharacter(SetSuper, &Ncharacter);
  TightBindingFermiCharacterUnfolding(SetSuper, Emin, Emax, NEpoints, k0, k1, k2, Nsupercell, Nsupercell, Fermi,
                                      Ncharacter, &EigenFunctionToCharactersOnsiteEnergyPhaseRGBkshifted);

//  FreeTightBindingDef(&Set);
//  FreeTightBindingDef(&SetSuper);

  return 0;
}


int TightBindingImpurityUnfoldingRandomImpurity1D(int Nsupercell, double Rho, double Vpot, char *filename,
                                                  MatrixType * BWValue, MatrixType * RGBValue)
{
  TightBindingDef Set, SetSuper;
  int    i, N[3];
  int    Ncharacter;
  int    QImpurity, NImpurity;
  char   CharFileName[255];
  double Emin, Emax, ELines[8];
  int    NEpoints, NELines;
  MatrixType AllValue;
  MatrixType *Character;
  BandstructurePathType QPath;
  MatrixType BANDS;
  CreateOneBandChainLaticeTightBindingDef(&Set);
//  PrintTightBindingDef(Set);
  N[0] = Nsupercell;
  N[1] = 1;
  N[2] = 1;
  CreateTightBindingSquareSuperCell(Set, &SetSuper, N);
  NImpurity = 0;
  for (i = 0; i < (Nsupercell); i++)
  {
    QImpurity = 0;
    if (Rho > mwhdrand())
      QImpurity = 1;
    SetSuper.H[3 * i].Hop.VRe[0][0] = Vpot * QImpurity;
    printf("%1i", QImpurity);
    NImpurity += QImpurity;
  }
  printf("\n");
  printf("Number of Impurities: %i\n", NImpurity);
  printf("Number of Sites     : %i\n", (Nsupercell));
  printf("Rho                 : %22.15E\n", ((double)NImpurity) / ((double)(Nsupercell)));
  printf("Rho average (input) : %22.15E\n", Rho);
//  PrintTightBindingDef(SetSuper);
  QPath.NQ = 2;
  QPath.Q = malloc(QPath.NQ * sizeof(BandstructurePointType));
  QPath.Q[0].Nq = 500;
  QPath.Q[0].label = 'G';
  QPath.Q[0].q[0] = 0.0 * Pi;
  QPath.Q[0].q[1] = 0.0 * Pi;
  QPath.Q[0].q[2] = 0.0 * Pi;
  QPath.Q[1].Nq = 1;
  QPath.Q[1].label = 'X';
  QPath.Q[1].q[0] = 1.0 * Pi;
  QPath.Q[1].q[1] = 0.0 * Pi;
  QPath.Q[1].q[2] = 0.0 * Pi;

  EigenFunctionToCharactersOnsiteEnergyPhaseSetNCharacter(SetSuper, &Ncharacter);
  Character = malloc(Ncharacter * sizeof(MatrixType));
  TightBindingBANDSCharacter(SetSuper, QPath, &BANDS, Ncharacter, Character,
                             &EigenFunctionToCharactersOnsiteEnergyPhase);
  sprintf(CharFileName, "%s.energy", filename);
  FilePrintPlainMatrix(BANDS, CharFileName, "w");
  for (i = 0; i < Ncharacter; i++)
  {
    if (i < 10)
      sprintf(CharFileName, "%s.character_00%1i", filename, i);
    else if (i < 100)
      sprintf(CharFileName, "%s.character_0%2i", filename, i);
    else if (i < 1000)
      sprintf(CharFileName, "%s.character_%3i", filename, i);
    else
      sprintf(CharFileName, "%s.character_%i", filename, i);
    FilePrintPlainMatrix(Character[i], CharFileName, "w");
  }

  Emin = -2.3;
  Emax = 5.3;
  NEpoints = 500;
  NELines = 8;
  ELines[0] = -2;
  ELines[1] = -1;
  ELines[2] = 0;
  ELines[3] = 1;
  ELines[4] = 2;
  ELines[5] = 3;
  ELines[6] = 4;
  ELines[7] = 5;

  TightBindingBANDSToMatrix(Emin, Emax, NEpoints, BANDS, &AllValue);
  TightBindingAddAxesToMatrix(Emin, Emax, NELines, ELines, QPath, &AllValue);
  sprintf(CharFileName, "%s_all.bmp", filename);
  FilePrintRealMatrixToBitmap(AllValue, CharFileName, &RealToRGBColorBWInverse);

  TightBindingBANDSCharacterToMatrix(Emin, Emax, NEpoints, BANDS, Character[2], BWValue);
  TightBindingAddAxesToMatrix(Emin, Emax, NELines, ELines, QPath, BWValue);
  sprintf(CharFileName, "%s_BW.bmp", filename);
  FilePrintRealMatrixToBitmap(*BWValue, CharFileName, &RealToRGBColorBWInverse);

  TightBindingBANDSCharacterToMatrix(Emin, Emax, NEpoints, BANDS, Character[0], &RGBValue[0]);
  TightBindingBANDSCharacterToMatrix(Emin, Emax, NEpoints, BANDS, Character[2], &RGBValue[1]);
  TightBindingBANDSCharacterToMatrix(Emin, Emax, NEpoints, BANDS, Character[1], &RGBValue[2]);
  TightBindingAddAxesToMatrix(Emin, Emax, NELines, ELines, QPath, &RGBValue[0]);
  TightBindingAddAxesToMatrix(Emin, Emax, NELines, ELines, QPath, &RGBValue[1]);
  TightBindingAddAxesToMatrix(Emin, Emax, NELines, ELines, QPath, &RGBValue[2]);
  sprintf(CharFileName, "%s_RGB.bmp", filename);
  FilePrintReal3MatrixToBitmap(RGBValue, CharFileName, &Real3ToRGBColorSimpleInverse);

  return 0;
}


int TightBindingImpurityUnfoldingRandomImpuritySeries(int Dimension, int Nsupercell, double Rho, double Vpot,
                                                      int repetitions, char *filename)
{
  int const NE = 1032 + 250, Ncharacter = 3;
  MatrixType BWValue, RGBValue[3];
  MatrixType SumBWValue, SumRGBValue[3];
  MatrixType Fermi[3 * (NE + 1)];
  SpectraType DOS[3], SumDOS[3];
  char   fullfilename[255];
  int    i, j;
  double intDOS, intDOSA, intDOSB;
  for (i = 0; i < 3 * (NE + 1); i++)
  {
    Fermi[i].QComplex = 0;
    Fermi[i].Ni = 320;
    Fermi[i].Nj = 320;
    if (InitMatrix(&(Fermi[i])))
    {
      printf("InitMatrix failed in TightBindingImpurityUnfoldingRandomImpuritySeries\n");
      fflush(stdout);
      return 1;
    };
  }
  mwhsrand((unsigned int)time(NULL));
  sprintf(fullfilename, "%s_single", filename);
  if (Dimension == 3)
    TightBindingImpurityUnfoldingRandomImpurity3D(Nsupercell, Rho, Vpot, fullfilename, &BWValue, RGBValue, DOS, Fermi);
  if (Dimension == 2)
    TightBindingImpurityUnfoldingRandomImpurity2D(Nsupercell, Rho, Vpot, fullfilename, &BWValue, RGBValue, DOS, Fermi);
  if (Dimension == 1)
    TightBindingImpurityUnfoldingRandomImpurity1D(Nsupercell, Rho, Vpot, fullfilename, &BWValue, RGBValue);
  SpectrumAllocate(DOS[0], &SumDOS[0]);
  SpectrumAllocate(DOS[1], &SumDOS[1]);
  SpectrumAllocate(DOS[2], &SumDOS[2]);
  CopySpectra(DOS[0], &SumDOS[0]);
  CopySpectra(DOS[1], &SumDOS[1]);
  CopySpectra(DOS[2], &SumDOS[2]);
  SpectrumFree(&DOS[0]);
  SpectrumFree(&DOS[1]);
  SpectrumFree(&DOS[2]);
  MatrixCopy(&SumBWValue, BWValue);
  MatrixCopy(&SumRGBValue[0], RGBValue[0]);
  MatrixCopy(&SumRGBValue[1], RGBValue[1]);
  MatrixCopy(&SumRGBValue[2], RGBValue[2]);
  FreeMatrix(&BWValue);
  FreeMatrix(&RGBValue[0]);
  FreeMatrix(&RGBValue[1]);
  FreeMatrix(&RGBValue[2]);
  for (i = 1; i < repetitions; i++)
  {
    fflush(stdout);
    printf("Repetition %i from %i\n", i, repetitions);
    if (Dimension == 3)
      TightBindingImpurityUnfoldingRandomImpurity3D(Nsupercell, Rho, Vpot, fullfilename, &BWValue, RGBValue, DOS,
                                                    Fermi);
    if (Dimension == 2)
      TightBindingImpurityUnfoldingRandomImpurity2D(Nsupercell, Rho, Vpot, fullfilename, &BWValue, RGBValue, DOS,
                                                    Fermi);
    if (Dimension == 1)
      TightBindingImpurityUnfoldingRandomImpurity1D(Nsupercell, Rho, Vpot, fullfilename, &BWValue, RGBValue);
    AddMultipliedSpectra(&SumDOS[0], DOS[0], 1, 0);
    AddMultipliedSpectra(&SumDOS[1], DOS[1], 1, 0);
    AddMultipliedSpectra(&SumDOS[2], DOS[2], 1, 0);
    SpectrumFree(&DOS[0]);
    SpectrumFree(&DOS[1]);
    SpectrumFree(&DOS[2]);
    MatrixAdd(&SumBWValue, BWValue, 1, 0);
    MatrixAdd(&SumRGBValue[0], RGBValue[0], 1, 0);
    MatrixAdd(&SumRGBValue[1], RGBValue[1], 1, 0);
    MatrixAdd(&SumRGBValue[2], RGBValue[2], 1, 0);
    FreeMatrix(&BWValue);
    FreeMatrix(&RGBValue[0]);
    FreeMatrix(&RGBValue[1]);
    FreeMatrix(&RGBValue[2]);
  }
  MatrixScalarMultiply(&SumBWValue, 1.0 / ((double)repetitions), 0);
  MatrixScalarMultiply(&SumRGBValue[0], 1.0 / ((double)repetitions), 0);
  MatrixScalarMultiply(&SumRGBValue[1], 1.0 / ((double)repetitions), 0);
  MatrixScalarMultiply(&SumRGBValue[2], 1.0 / ((double)repetitions), 0);
  SpectraScalarMultiply(&SumDOS[0], 1.0 / ((double)repetitions), 0);
  SpectraScalarMultiply(&SumDOS[1], 1.0 / ((double)repetitions), 0);
  SpectraScalarMultiply(&SumDOS[2], 1.0 / ((double)repetitions), 0);
  sprintf(fullfilename, "%s_All_DOS.dat", filename);
  FilePrintSpectra(SumDOS, 3, fullfilename, "w");
  intDOS = 0;
  intDOSA = 0;
  intDOSB = 0;
  j = 0;
  for (i = 0; i < NE; i++)
  {
    intDOS += SumDOS[0].Re[i];
    intDOSA += SumDOS[1].Re[i];
    intDOSB += SumDOS[2].Re[i];
    if ((intDOS >= 0.5 * Rho) && (!(j & 1)))
    {
      j += 1;
      printf("Fermi energy for s0 doped with s1: %4i %lf filling of impurity: %lf of material: %lf\n", i,
             i * 0.002 - 1.532, intDOSA, intDOSB);
    }
    if ((intDOS >= 0.5 * 1.0) && (!(j & 2)))
    {
      j += 2;
      printf("Fermi energy for s1 doped with s1: %4i %lf filling of impurity: %lf of material: %lf\n", i,
             i * 0.002 - 1.532, intDOSA, intDOSB);
    }
    if ((intDOS >= 0.5 * (1.0 + Rho)) && (!(j & 4)))
    {
      j += 4;
      printf("Fermi energy for s1 doped with s2: %4i %lf filling of impurity: %lf of material: %lf\n", i,
             i * 0.002 - 1.532, intDOSA, intDOSB);
    }
  }
  intDOS = 1.0;
  intDOSA = Rho;
  intDOSB = 1.0 - Rho;
  j = 0;
  for (i = NE; i > 0; i--)
  {
    intDOS -= SumDOS[0].Re[i];
    intDOSA -= SumDOS[1].Re[i];
    intDOSB -= SumDOS[2].Re[i];
    if ((intDOS <= 0.5 * Rho) && (!(j & 1)))
    {
      j += 1;
      printf("Fermi energy for s0 doped with s1: %4i %lf filling of impurity: %lf of material: %lf\n", i,
             i * 0.002 - 1.532, intDOSA, intDOSB);
    }
    if ((intDOS <= 0.5 * 1.0) && (!(j & 2)))
    {
      j += 2;
      printf("Fermi energy for s1 doped with s1: %4i %lf filling of impurity: %lf of material: %lf\n", i,
             i * 0.002 - 1.532, intDOSA, intDOSB);
    }
    if ((intDOS <= 0.5 * (1.0 + Rho)) && (!(j & 4)))
    {
      j += 4;
      printf("Fermi energy for s1 doped with s2: %4i %lf filling of impurity: %lf of material: %lf\n", i,
             i * 0.002 - 1.532, intDOSA, intDOSB);
    }
  }

  sprintf(fullfilename, "%s_All.dat", filename);
  FilePrintPlainMatrix(SumBWValue, fullfilename, "w");
  sprintf(fullfilename, "%s_All_BW_1.bmp", filename);
  FilePrintRealMatrixToBitmap(SumBWValue, fullfilename, &RealToRGBColorBWInverse);
  sprintf(fullfilename, "%s_All_RGB_1.bmp", filename);
  FilePrintReal3MatrixToBitmap(SumRGBValue, fullfilename, &Real3ToRGBColorSimpleInverse);
  MatrixScalarMultiply(&SumBWValue, 2, 0);
  MatrixScalarMultiply(&SumRGBValue[0], 2, 0);
  MatrixScalarMultiply(&SumRGBValue[1], 2, 0);
  MatrixScalarMultiply(&SumRGBValue[2], 2, 0);
  sprintf(fullfilename, "%s_All_BW_2.bmp", filename);
  FilePrintRealMatrixToBitmap(SumBWValue, fullfilename, &RealToRGBColorBWInverse);
  sprintf(fullfilename, "%s_All_RGB_2.bmp", filename);
  FilePrintReal3MatrixToBitmap(SumRGBValue, fullfilename, &Real3ToRGBColorSimpleInverse);
  MatrixScalarMultiply(&SumBWValue, 2, 0);
  MatrixScalarMultiply(&SumRGBValue[0], 2, 0);
  MatrixScalarMultiply(&SumRGBValue[1], 2, 0);
  MatrixScalarMultiply(&SumRGBValue[2], 2, 0);
  sprintf(fullfilename, "%s_All_BW_4.bmp", filename);
  FilePrintRealMatrixToBitmap(SumBWValue, fullfilename, &RealToRGBColorBWInverse);
  sprintf(fullfilename, "%s_All_RGB_4.bmp", filename);
  FilePrintReal3MatrixToBitmap(SumRGBValue, fullfilename, &Real3ToRGBColorSimpleInverse);
  MatrixScalarMultiply(&SumBWValue, 2, 0);
  MatrixScalarMultiply(&SumRGBValue[0], 2, 0);
  MatrixScalarMultiply(&SumRGBValue[1], 2, 0);
  MatrixScalarMultiply(&SumRGBValue[2], 2, 0);
  sprintf(fullfilename, "%s_All_BW_8.bmp", filename);
  FilePrintRealMatrixToBitmap(SumBWValue, fullfilename, &RealToRGBColorBWInverse);
  sprintf(fullfilename, "%s_All_RGB_8.bmp", filename);
  FilePrintReal3MatrixToBitmap(SumRGBValue, fullfilename, &Real3ToRGBColorSimpleInverse);

  for (i = 0; i <= NE; i++)
  {
    for (j = 0; j < Ncharacter; j++)
    {
      MatrixScalarMultiply(&Fermi[Ncharacter * i + j], 1.0 / ((double)repetitions), 0);
      if (i < 10)
      {
        sprintf(fullfilename, "%s_All_Fermi_0000%1i_%1i.dat", filename, i, j);
      }
      else if (i < 100)
      {
        sprintf(fullfilename, "%s_All_Fermi_000%2i_%1i.dat", filename, i, j);
      }
      else if (i < 1000)
      {
        sprintf(fullfilename, "%s_All_Fermi_00%3i_%1i.dat", filename, i, j);
      }
      else if (i < 10000)
      {
        sprintf(fullfilename, "%s_All_Fermi_0%4i_%1i.dat", filename, i, j);
      }
      else
      {
        sprintf(fullfilename, "%s_All_Fermi_%5i_%1i.dat", filename, i, j);
      }
//      FilePrintPlainMatrix(Fermi[Ncharacter*i+j], fullfilename, "w");
      MatrixScalarMultiply(&Fermi[Ncharacter * i + j], 32.0, 0);
    }
    if (i < 10)
    {
      sprintf(fullfilename, "%s_All_Fermi_0000%1i.bmp", filename, i);
    }
    else if (i < 100)
    {
      sprintf(fullfilename, "%s_All_Fermi_000%2i.bmp", filename, i);
    }
    else if (i < 1000)
    {
      sprintf(fullfilename, "%s_All_Fermi_00%3i.bmp", filename, i);
    }
    else if (i < 10000)
    {
      sprintf(fullfilename, "%s_All_Fermi_0%4i.bmp", filename, i);
    }
    else
    {
      sprintf(fullfilename, "%s_All_Fermi_%5i.bmp", filename, i);
    }
    FilePrintReal3MatrixToBitmap(&(Fermi[Ncharacter * i]), fullfilename, &Real3ToRGBColorSimpleInverse);
  }
  return 0;
}
