void TestStruct()
{
/*
*   Just set up some structures, fill them and print them
*/
  printf("Start of TestStruct\n");
  printf("Initiallize several structures and print them\n");
  printf("Complex16\n");
  Complex16 a, b;
  a.Re = 12;
  a.Im = 5;
  b.Re = 5;
  b.Im = 6;
  printf("a  =%22.15E + I %22.15E\nb  =%22.15E + I %22.15E\n", a.Re, a.Im, b.Re, b.Im);
  a.Re += b.Re;
  a.Im += b.Im;
  printf("a+b=%22.15E + I %22.15E\n", a.Re, a.Im);
  printf("WaveFunctionType\n");
  printf("OperatorFixedLengthType\n");
  printf("OperatorType\n");
  printf("MatrixType\n");
  printf("SparseMatrixType\n");
  printf("TriDiagonalMatrixType\n");
  printf("AndersonMatrixType\n");
  printf("BlokTriDiagonalMatrixType\n");
  printf("BlokAndersonMatrixType\n");
  printf("DoubleTriDiagonalMatrixType\n");
  printf("BlokDoubleTriDiagonalMatrixType\n");
  printf("ChebyshevExpansion\n");
  printf("HarmonicsExpansion\n");
  printf("SpectraType\n");
};


void TestHash()
{
  printf("uint32_t Hash(&Psi[0],NWord,NBitsKey)                   \n");
  printf("returns a Key of Psi between 0 and 2^NBitsKey-1         \n");
  printf("NBitsKey can range from 0 to 31                         \n");
  printf("NWord is the length of the array Psi in words (8-bit).  \n");
  printf("Psi is an array of uint8_t                              \n");


  uint32_t NBitsKey;
  uint32_t NorbWord, Norb;
  uint32_t Key;
  NorbWord = 8;
  Norb = 8 * NorbWord;
  uint8_t Psi[NorbWord];
  Psi[0] = 0xFF;
  Psi[1] = 0xFF;
  Psi[2] = 0xFF;
  Psi[3] = 0xFF;
  Psi[4] = 0xF0;
  Psi[5] = 0x00;
  Psi[6] = 0x00;
  Psi[7] = 0x01;
  NBitsKey = 2;
  Key = hash(Psi, NorbWord, NBitsKey);
  printf("Key= %X\n", Key);
  NBitsKey = 4;
  Key = hash(Psi, NorbWord, NBitsKey);
  printf("Key= %X\n", Key);
  NBitsKey = 6;
  Key = hash(Psi, NorbWord, NBitsKey);
  printf("Key= %X\n", Key);
  NBitsKey = 8;
  Key = hash(Psi, NorbWord, NBitsKey);
  printf("Key= %X\n", Key);
  NBitsKey = 12;
  Key = hash(Psi, NorbWord, NBitsKey);
  printf("Key= %X\n", Key);
  NBitsKey = 16;
  Key = hash(Psi, NorbWord, NBitsKey);
  printf("Key= %X\n", Key);
  NBitsKey = 20;
  Key = hash(Psi, NorbWord, NBitsKey);
  printf("Key= %X\n", Key);
  NBitsKey = 24;
  Key = hash(Psi, NorbWord, NBitsKey);
  printf("Key= %X\n", Key);
  NBitsKey = 28;
  Key = hash(Psi, NorbWord, NBitsKey);
  printf("Key= %X\n", Key);
  NBitsKey = 31;
  Key = hash(Psi, NorbWord, NBitsKey);
  printf("Key= %X\n", Key);

};

void TestAddFermion()
{
  WaveFunctionType psi, psiNew;
  uint32_t i;
  uint8_t det[10];
  psi.QComplex = 0;
  psi.NFermion = 16;
  psi.NBoson = 1;
  psi.NBitsKey = 10;
  WaveFunctionInit(&psi);
  for (i = 1; i < 10; i++)
  {
    det[0] = (uint8_t) i;
    det[1] = (uint8_t) 255;
    det[2] = (uint8_t) 16;
    RealWaveFunctionAddElement(&psi, det, 0.5 / ((double)i));
  }
  WaveFunctionAddFermion(psi, &psiNew, 3, 0);
  PrintWaveFunctionFullReadable(psi);
  PrintWaveFunctionFullReadable(psiNew);
}


void TestWaveFunctionArithmetic()
{
  uint32_t i;
  WaveFunctionType psi;
  uint8_t *Determinant;
  double a;
  strcpy(psi.Name, "Test WaveFunction");
  psi.QComplex = 0;
  psi.N = 0;
  psi.NBlok = 1;
  psi.NFermion = 20;
  psi.NBoson = 0;
  psi.NWord = 3;
  psi.NBitsKey = 10;
  psi.KeyList = calloc((uint32_t) 1 << (psi.NBitsKey), sizeof(uint32_t));
  psi.ValueRe = malloc(psi.NBlok * sizeof(double *));
  psi.Next = malloc(psi.NBlok * sizeof(uint32_t *));
  psi.Determinant = malloc(psi.NBlok * sizeof(uint8_t *));
  if ((psi.ValueRe == NULL) || (psi.Next == NULL) || (psi.Determinant == NULL))
    AllocFailed();
  psi.ValueRe[psi.NBlok - 1] = malloc(BlokSize * sizeof(double));
  psi.Next[psi.NBlok - 1] = malloc(BlokSize * sizeof(uint32_t));
  psi.Determinant[psi.NBlok - 1] = malloc(BlokSize * psi.NWord * sizeof(uint8_t));
  if ((psi.ValueRe[psi.NBlok - 1] == NULL) || (psi.Next[psi.NBlok - 1] == NULL)
      || (psi.Determinant[psi.NBlok - 1] == NULL))
    AllocFailed();
  Determinant = calloc(psi.NWord, sizeof(uint8_t));
  Determinant[0] = 255;
  Determinant[1] = 64;
  Determinant[2] = 0;
  a = 0.123;
  RealWaveFunctionAddElement(&psi, Determinant, a);
  Determinant[0] = 255;
  Determinant[1] = 65;
  Determinant[2] = 0;
  a = 1000.123456;
  RealWaveFunctionAddElement(&psi, Determinant, a);
  Determinant[0] = 255;
  Determinant[1] = 66;
  Determinant[2] = 0;
  a = 0.123456789;
  RealWaveFunctionAddElement(&psi, Determinant, a);
  PrintWaveFunctionFullReadable(psi);

  for (i = 1; i < 46; i++)
  {
    Determinant[0] = 255;
    Determinant[1] = i;
    Determinant[2] = 0;
    RealWaveFunctionAddElement(&psi, Determinant, i);
  }

  Determinant[0] = 255;
  Determinant[1] = 30;
  Determinant[2] = 0;
  a = 0.123456789;
  RealWaveFunctionAddElement(&psi, Determinant, a);

  RealWaveFunctionNormalize(&psi);
  printf("Wave Function\n");
  PrintWaveFunctionFullReadable(psi);
  RealWaveFunctionTruncateBasis(&psi, 0.01);
  printf("Wave Function truncated\n");
  PrintWaveFunctionFullReadable(psi);

  RealWaveFunctionFree(&psi);

  strcpy(psi.Name, "Second Test WaveFunction");
  psi.QComplex = 0;
  psi.N = 0;
  psi.NBlok = 1;
  psi.NFermion = 20;
  psi.NBoson = 0;
  psi.NWord = 3;
  psi.NBitsKey = 10;
  psi.KeyList = calloc((uint32_t) 1 << (psi.NBitsKey), sizeof(uint32_t));
  psi.ValueRe = malloc(psi.NBlok * sizeof(double *));
  psi.Next = malloc(psi.NBlok * sizeof(uint32_t *));
  psi.Determinant = malloc(psi.NBlok * sizeof(uint8_t *));
  if ((psi.ValueRe == NULL) || (psi.Next == NULL) || (psi.Determinant == NULL))
    AllocFailed();
  psi.ValueRe[psi.NBlok - 1] = malloc(BlokSize * sizeof(double));
  psi.Next[psi.NBlok - 1] = malloc(BlokSize * sizeof(uint32_t));
  psi.Determinant[psi.NBlok - 1] = malloc(BlokSize * psi.NWord * sizeof(uint8_t));
  if ((psi.ValueRe[psi.NBlok - 1] == NULL) || (psi.Next[psi.NBlok - 1] == NULL)
      || (psi.Determinant[psi.NBlok - 1] == NULL))
    AllocFailed();


  for (i = 1; i < 50; i++)
  {
    Determinant[0] = 255;
    Determinant[1] = i;
    Determinant[2] = 0;
    RealWaveFunctionAddElement(&psi, Determinant, i);
  }

  PrintWaveFunctionFullReadable(psi);

  for (i = 25; i < 75; i++)
  {
    Determinant[0] = 255;
    Determinant[1] = i;
    Determinant[2] = 0;
    RealWaveFunctionAddElementConserveBasis(&psi, Determinant, i);
  }
  PrintWaveFunctionFullReadable(psi);

};

void TestOperatorArithmeticNormalOrder()
{
  uint32_t i, j;
  uint16_t oadd[6];
  double k[3], r[3], NCell = 4;
  OperatorType O;
  InitOperator(&O);
  O.NFermion = 20 * NCell;
  O.NWord = Min(O.NFermion, 1 + (O.NFermion - 1) / 8);
  InitOperatorFixedLength(&O, 2, 8);
  k[0] = 0;
  k[1] = 0;
  k[2] = Pi;
  for (i = 0; i < NCell; i++)
  {
    r[0] = 0;
    r[1] = 0;
    r[2] = i * 1.0 / NCell;
    oadd[1] = CrAnMask + 1 + 2 * i;
    oadd[0] = 1 + 2 * i;
    ComplexOperatorAddLaderNormalOrder(&O, oadd, 2, cos(Dot(r, k)) * sqrt(0.25 / (NCell)),
                                       sin(Dot(r, k)) * sqrt(0.25 / (NCell)), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 2 * i;
    oadd[0] = 2 * i;
    ComplexOperatorAddLaderNormalOrder(&O, oadd, 2, cos(Dot(r, k)) * sqrt(0.25 / (NCell)),
                                       sin(Dot(r, k)) * sqrt(0.25 / (NCell)), &QDetNotOrderedFcFaBcBa);
  }
  PrintOperatorFullReadable(O);
  InitOperatorFixedLength(&O, 4, 8);
  InitOperatorFixedLength(&O, 6, 8);
  for (i = 0; i < NCell; i++)
  {
    j = (i + 1) % (uint32_t) NCell;
    printf("i=%i j=%i\n", i, j);
    oadd[5] = CrAnMask + 2 * i + 1;
    oadd[4] = 2 * i + 1;
    oadd[3] = CrAnMask + 2 * j + 1;
    oadd[2] = CrAnMask + 2 * j;
    oadd[1] = 2 * j;
    oadd[0] = 2 * j + 1;
    RealOperatorAddLaderNormalOrder(&O, oadd, 6, 0.15, &QDetNotOrderedFcFaBcBa);
    oadd[5] = CrAnMask + 2 * i;
    oadd[4] = 2 * i;
    oadd[3] = CrAnMask + 2 * j + 1;
    oadd[2] = CrAnMask + 2 * j;
    oadd[1] = 2 * j;
    oadd[0] = 2 * j + 1;
    RealOperatorAddLaderNormalOrder(&O, oadd, 6, 0.15, &QDetNotOrderedFcFaBcBa);
    oadd[5] = CrAnMask + 2 * j + 1;
    oadd[4] = 2 * j + 1;
    oadd[3] = CrAnMask + 2 * i + 1;
    oadd[2] = CrAnMask + 2 * i;
    oadd[1] = 2 * i;
    oadd[0] = 2 * i + 1;
    RealOperatorAddLaderNormalOrder(&O, oadd, 6, 0.15, &QDetNotOrderedFcFaBcBa);
    oadd[5] = CrAnMask + 2 * j;
    oadd[4] = 2 * j;
    oadd[3] = CrAnMask + 2 * i + 1;
    oadd[2] = CrAnMask + 2 * i;
    oadd[1] = 2 * i;
    oadd[0] = 2 * i + 1;
    RealOperatorAddLaderNormalOrder(&O, oadd, 6, 0.15, &QDetNotOrderedFcFaBcBa);

    oadd[5] = CrAnMask + 2 * i;
    oadd[4] = CrAnMask + 2 * i + 1;
    oadd[3] = 2 * i + 1;
    oadd[2] = CrAnMask + 2 * j + 1;
    oadd[1] = 2 * j + 1;
    oadd[0] = 2 * j;
    RealOperatorAddLaderNormalOrder(&O, oadd, 6, -1.0, &QDetNotOrderedFcFaBcBa);
    oadd[5] = CrAnMask + 2 * i + 1;
    oadd[4] = CrAnMask + 2 * i;
    oadd[3] = 2 * i;
    oadd[2] = CrAnMask + 2 * j;
    oadd[1] = 2 * j;
    oadd[0] = 2 * j + 1;
    RealOperatorAddLaderNormalOrder(&O, oadd, 6, -1.0, &QDetNotOrderedFcFaBcBa);
    oadd[5] = CrAnMask + 2 * j;
    oadd[4] = CrAnMask + 2 * j + 1;
    oadd[3] = 2 * j + 1;
    oadd[2] = CrAnMask + 2 * i + 1;
    oadd[1] = 2 * i + 1;
    oadd[0] = 2 * i;
    RealOperatorAddLaderNormalOrder(&O, oadd, 6, -1.0, &QDetNotOrderedFcFaBcBa);
    oadd[5] = CrAnMask + 2 * j + 1;
    oadd[4] = CrAnMask + 2 * j;
    oadd[3] = 2 * j;
    oadd[2] = CrAnMask + 2 * i;
    oadd[1] = 2 * i;
    oadd[0] = 2 * i + 1;
    RealOperatorAddLaderNormalOrder(&O, oadd, 6, -1.0, &QDetNotOrderedFcFaBcBa);
    j = (uint32_t) (((int)i + 3) % (int)NCell);
    printf("i=%i j=%i\n", i, j);
    oadd[5] = CrAnMask + 2 * i + 1;
    oadd[4] = 2 * i + 1;
    oadd[3] = CrAnMask + 2 * j + 1;
    oadd[2] = CrAnMask + 2 * j;
    oadd[1] = 2 * j;
    oadd[0] = 2 * j + 1;
    RealOperatorAddLaderNormalOrder(&O, oadd, 6, 0.15, &QDetNotOrderedFcFaBcBa);
    oadd[5] = CrAnMask + 2 * i;
    oadd[4] = 2 * i;
    oadd[3] = CrAnMask + 2 * j + 1;
    oadd[2] = CrAnMask + 2 * j;
    oadd[1] = 2 * j;
    oadd[0] = 2 * j + 1;
    RealOperatorAddLaderNormalOrder(&O, oadd, 6, 0.15, &QDetNotOrderedFcFaBcBa);
    oadd[5] = CrAnMask + 2 * j + 1;
    oadd[4] = 2 * j + 1;
    oadd[3] = CrAnMask + 2 * i + 1;
    oadd[2] = CrAnMask + 2 * i;
    oadd[1] = 2 * i;
    oadd[0] = 2 * i + 1;
    RealOperatorAddLaderNormalOrder(&O, oadd, 6, 0.15, &QDetNotOrderedFcFaBcBa);
    oadd[5] = CrAnMask + 2 * j;
    oadd[4] = 2 * j;
    oadd[3] = CrAnMask + 2 * i + 1;
    oadd[2] = CrAnMask + 2 * i;
    oadd[1] = 2 * i;
    oadd[0] = 2 * i + 1;
    RealOperatorAddLaderNormalOrder(&O, oadd, 6, 0.15, &QDetNotOrderedFcFaBcBa);

    oadd[5] = CrAnMask + 2 * i;
    oadd[4] = CrAnMask + 2 * i + 1;
    oadd[3] = 2 * i + 1;
    oadd[2] = CrAnMask + 2 * j + 1;
    oadd[1] = 2 * j + 1;
    oadd[0] = 2 * j;
    RealOperatorAddLaderNormalOrder(&O, oadd, 6, -1.0, &QDetNotOrderedFcFaBcBa);
    oadd[5] = CrAnMask + 2 * i + 1;
    oadd[4] = CrAnMask + 2 * i;
    oadd[3] = 2 * i;
    oadd[2] = CrAnMask + 2 * j;
    oadd[1] = 2 * j;
    oadd[0] = 2 * j + 1;
    RealOperatorAddLaderNormalOrder(&O, oadd, 6, -1.0, &QDetNotOrderedFcFaBcBa);
    oadd[5] = CrAnMask + 2 * j;
    oadd[4] = CrAnMask + 2 * j + 1;
    oadd[3] = 2 * j + 1;
    oadd[2] = CrAnMask + 2 * i + 1;
    oadd[1] = 2 * i + 1;
    oadd[0] = 2 * i;
    RealOperatorAddLaderNormalOrder(&O, oadd, 6, -1.0, &QDetNotOrderedFcFaBcBa);
    oadd[5] = CrAnMask + 2 * j + 1;
    oadd[4] = CrAnMask + 2 * j;
    oadd[3] = 2 * j;
    oadd[2] = CrAnMask + 2 * i;
    oadd[1] = 2 * i;
    oadd[0] = 2 * i + 1;
    RealOperatorAddLaderNormalOrder(&O, oadd, 6, -1.0, &QDetNotOrderedFcFaBcBa);
  }
  PrintOperatorFullReadable(O);
}


void TestOperatorArithmetic()
{
  MatrixType U;
  uint16_t LaderOperator[4], LaderOperatorSec[4];
  OperatorType O, A, B, AB, cO, OO;
  InitOperator(&O);
  PrintOperatorFullReadable(O);
  InitOperatorFixedLength(&O, 2, 8);
  InitOperatorFixedLength(&O, 4, 8);
  PrintOperatorFullReadable(O);
  O.NFermion = 5;
  O.NBoson = 2;
  O.NWord = 3;
  LaderOperator[0] = 1;
  LaderOperator[1] = 4;
  LaderOperator[2] = CrAnMask | 1;
  LaderOperator[3] = CrAnMask | 4;
  LaderOperatorSec[0] = 6;
  LaderOperatorSec[1] = CrAnMask | 6;
  LaderOperatorSec[2] = CrAnMask | 6;
  LaderOperatorSec[3] = 7;
  RealOperatorAddLader(&O, LaderOperator, 2, 1.0);
  RealOperatorAddLader(&O, LaderOperator, 2, 1.0);
  RealOperatorAddLader(&O, LaderOperatorSec, 2, 1.0);
  RealOperatorAddLader(&O, LaderOperator, 2, 1.0);
  PrintOperatorFullReadable(O);
  OperatorRemoveSmallValues(&O, 2);
  PrintOperatorFullReadable(O);
  RealOperatorAddLader(&O, LaderOperator, 1, 1.0);
  RealOperatorAddLader(&O, LaderOperatorSec, 1, 1.0);
  RealOperatorAddLader(&O, LaderOperatorSec, 2, 1.0);
  PrintOperatorFullReadable(O);
  OperatorRemoveSmallValues(&O, 10);
  RealOperatorAddLaderNormalOrder(&O, LaderOperator, 4, 1.0, &QDetNotOrderedFaFcBcBa);
  RealOperatorAddLaderNormalOrder(&O, LaderOperatorSec, 4, 1.0, &QDetNotOrderedFaFcBcBa);
  PrintOperatorFullReadable(O);
  printf("\n\n\n\n\n===============\n\n\n\n\n");
  InitOperator(&A);
  InitOperatorFixedLength(&A, 2, 8);
  InitOperatorCopySettings(&AB, O);
  InitOperatorCopySettings(&B, A);
  RealOperatorAddLaderNormalOrder(&A, LaderOperator, 2, 1.0, &QDetNotOrderedFcFaBcBa);
  RealOperatorAddLaderNormalOrder(&B, LaderOperatorSec, 2, 1.0, &QDetNotOrderedFcFaBcBa);
  OperatorMultiply(&AB, A, B);
  PrintOperatorFullReadable(A);
  PrintOperatorFullReadable(B);
  PrintOperatorFullReadable(AB);
  printf("\n\n\n\n\n===============\n\n\n\n\n");
  FreeOperator(&A);
  FreeOperator(&B);
  InitOperator(&A);
  InitOperatorFixedLength(&A, 2, 8);
  A.NFermion = 5;
  A.NBoson = 2;
  A.NWord = 3;
  InitOperatorCopySettings(&B, A);
  LaderOperator[0] = CrAnMask | 1;
  LaderOperator[1] = 0;
  RealOperatorAddLader(&A, LaderOperator, 2, 1.0);
  LaderOperator[0] = CrAnMask | 0;
  LaderOperator[1] = 1;
  RealOperatorAddLader(&A, LaderOperator, 2, 1.0);
  RealOperatorAddLader(&A, LaderOperator, 1, 3.0);
  RealOperatorAddLader(&A, LaderOperator, 0, -3.0);
  PrintOperatorFullReadable(A);
  PrintOperatorFullReadable(B);
  U.QComplex = 0;
  U.Ni = 7;
  U.Nj = 7;
  InitMatrix(&U);
  U.VRe[0][0] = sqrt(0.5);
  U.VRe[0][1] = sqrt(0.5);
  U.VRe[1][0] = sqrt(0.5);
  U.VRe[1][1] = -sqrt(0.5);
  U.VRe[2][2] = 1;
  U.VRe[3][3] = 1;
  U.VRe[4][4] = 1;
  U.VRe[5][5] = 1;
  U.VRe[6][6] = 1;
  PrintMatrix(U);
  MappedOperatorAdd(&B, A, U);
  PrintOperatorFullReadable(B);
  printf("\n\n ***** Operator copy \n\n");
  OperatorCopy(&cO, O);
  PrintOperatorReadable(O);
  PrintOperatorReadable(cO);

  InitOperator(&OO);
  OO.NFermion = 5;
  OO.NBoson = 2;
  OO.NWord = 3;
  InitOperatorFixedLength(&OO, 2, 8);
  InitOperatorFixedLength(&OO, 4, 8);
  InitOperatorFixedLength(&OO, 6, 8);
  InitOperatorFixedLength(&OO, 8, 8);
  InitOperatorFixedLength(&OO, 10, 8);
  InitOperatorFixedLength(&OO, 12, 8);
  OperatorMultiply3NormalOrder(&OO, O, O, O, 'N', 'N', 'N', &QDetNotOrderedFcFaBcBa);
  OperatorChop(&OO, DBL_EPSILON);
  PrintOperatorReadable(OO);
};

void TestCreate()
{

};

void TestVectorMatrixParallel()
{

  uint32_t i;
  WaveFunctionType psi;
  uint8_t *Determinant;
  double a;
  psi.QComplex = 0;
  psi.N = 0;
  psi.NBlok = 1;
  psi.NFermion = 20;
  psi.NBoson = 0;
  psi.NWord = 3;
  psi.NBitsKey = 10;
  psi.KeyList = calloc((uint32_t) 1 << (psi.NBitsKey), sizeof(uint32_t));
  psi.ValueRe = malloc(psi.NBlok * sizeof(double *));
  psi.Next = malloc(psi.NBlok * sizeof(uint32_t *));
  psi.Determinant = malloc(psi.NBlok * sizeof(uint8_t *));
  if ((psi.ValueRe == NULL) || (psi.Next == NULL) || (psi.Determinant == NULL))
    AllocFailed();
  psi.ValueRe[psi.NBlok - 1] = malloc(BlokSize * sizeof(double));
  psi.Next[psi.NBlok - 1] = malloc(BlokSize * sizeof(uint32_t));
  psi.Determinant[psi.NBlok - 1] = malloc(BlokSize * psi.NWord * sizeof(uint8_t));
  if ((psi.ValueRe[psi.NBlok - 1] == NULL) || (psi.Next[psi.NBlok - 1] == NULL)
      || (psi.Determinant[psi.NBlok - 1] == NULL))
    AllocFailed();
  Determinant = calloc(psi.NWord, sizeof(uint8_t));
  Determinant[0] = 255;
  Determinant[1] = 64;
  Determinant[2] = 0;
  a = 0.123;
  RealWaveFunctionAddElement(&psi, Determinant, a);
  Determinant[0] = 255;
  Determinant[1] = 65;
  Determinant[2] = 0;
  a = 1000.123456;
  RealWaveFunctionAddElement(&psi, Determinant, a);
  Determinant[0] = 255;
  Determinant[1] = 66;
  Determinant[2] = 0;
  a = 0.123456789;
  RealWaveFunctionAddElement(&psi, Determinant, a);
  PrintWaveFunctionFullReadable(psi);
  for (i = 1; i < 46; i++)
  {
    Determinant[0] = 255;
    Determinant[1] = i;
    Determinant[2] = 0;
    RealWaveFunctionAddElement(&psi, Determinant, i);
  }



  Determinant[0] = 255;
  Determinant[1] = 30;
  Determinant[2] = 0;
  a = 0.123456789;
  RealWaveFunctionAddElement(&psi, Determinant, a);


  printf("Wave Function AAA\n");
  PrintWaveFunctionFullReadable(psi);

  for (i = 0; i < 1000; i++)
  {
    RealWaveFunctionMultiply(psi, 2.0);
    RealWaveFunctionNormalize(&psi);
  }

  printf("Wave Function\n");
  PrintWaveFunctionFullReadable(psi);
  RealWaveFunctionTruncateBasis(&psi, 0.01);
  printf("Wave Function truncated\n");
  PrintWaveFunctionFullReadable(psi);

  RealWaveFunctionFree(&psi);
  psi.N = 0;
  psi.NBlok = 1;
  psi.NFermion = 20;
  psi.NBoson = 0;
  psi.NWord = 3;
  psi.NBitsKey = 10;
  psi.KeyList = calloc((uint32_t) 1 << (psi.NBitsKey), sizeof(uint32_t));
  psi.ValueRe = malloc(psi.NBlok * sizeof(double *));
  psi.Next = malloc(psi.NBlok * sizeof(uint32_t *));
  psi.Determinant = malloc(psi.NBlok * sizeof(uint8_t *));
  if ((psi.ValueRe == NULL) || (psi.Next == NULL) || (psi.Determinant == NULL))
    AllocFailed();
  psi.ValueRe[psi.NBlok - 1] = malloc(BlokSize * sizeof(double));
  psi.Next[psi.NBlok - 1] = malloc(BlokSize * sizeof(uint32_t));
  psi.Determinant[psi.NBlok - 1] = malloc(BlokSize * psi.NWord * sizeof(uint8_t));
  if ((psi.ValueRe[psi.NBlok - 1] == NULL) || (psi.Next[psi.NBlok - 1] == NULL)
      || (psi.Determinant[psi.NBlok - 1] == NULL))
    AllocFailed();


  for (i = 1; i < 50; i++)
  {
    Determinant[0] = 255;
    Determinant[1] = i;
    Determinant[2] = 0;
    RealWaveFunctionAddElement(&psi, Determinant, i);
  }

  PrintWaveFunctionFullReadable(psi);

  for (i = 25; i < 75; i++)
  {
    Determinant[0] = 255;
    Determinant[1] = i;
    Determinant[2] = 0;
    RealWaveFunctionAddElementConserveBasis(&psi, Determinant, i);
  }
  PrintWaveFunctionFullReadable(psi);

}



void TestLanczos()
{
  uint32_t i;
  uint16_t oadd[6];
  double NCell = 15;
  OperatorType O;
  double epsilon, zero;
  WaveFunctionType psi;
  uint8_t *Determinant;
  double a;
  psi.QComplex = 0;
  psi.N = 0;
  psi.NBlok = 1;
  psi.NFermion = NCell;
  psi.NBoson = 0;
  psi.NWord = Min(psi.NFermion, 1 + (psi.NFermion - 1) / 8);
  psi.NBitsKey = 10;
  psi.KeyList = calloc((uint32_t) 1 << (psi.NBitsKey), sizeof(uint32_t));
  psi.ValueRe = malloc(psi.NBlok * sizeof(double *));
  psi.Next = malloc(psi.NBlok * sizeof(uint32_t *));
  psi.Determinant = malloc(psi.NBlok * sizeof(uint8_t *));
  if ((psi.ValueRe == NULL) || (psi.Next == NULL) || (psi.Determinant == NULL))
    AllocFailed();
  psi.ValueRe[psi.NBlok - 1] = malloc(BlokSize * sizeof(double));
  psi.Next[psi.NBlok - 1] = malloc(BlokSize * sizeof(uint32_t));
  psi.Determinant[psi.NBlok - 1] = malloc(BlokSize * psi.NWord * sizeof(uint8_t));
  if ((psi.ValueRe[psi.NBlok - 1] == NULL) || (psi.Next[psi.NBlok - 1] == NULL)
      || (psi.Determinant[psi.NBlok - 1] == NULL))
    AllocFailed();
  Determinant = calloc(psi.NWord, sizeof(uint8_t));
  for (i = 0; i < psi.NWord; i++)
  {
    Determinant[i] = 15;
  }
  a = 1;
  Determinant[0] = 3;
  RealWaveFunctionAddElement(&psi, Determinant, a);
  a = 1;
  Determinant[0] = 5;
  RealWaveFunctionAddElement(&psi, Determinant, a);
  RealWaveFunctionNormalize(&psi);
  PrintWaveFunctionFullReadable(psi);

  InitOperator(&O);
  O.NFermion = NCell;
  O.NBoson = 0;
  O.NWord = Min(O.NFermion, 1 + (O.NFermion - 1) / 8);
  O.QComplex = 0;
  InitOperatorFixedLength(&O, 2, 8);
  InitOperatorFixedLength(&O, 4, 8);
  for (i = 1; i < NCell; i++)
  {
    oadd[1] = CrAnMask + i - 1;
    oadd[0] = i;
    OperatorAddLaderNormalOrder(&O, oadd, 2, 1.0, 0.0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i;
    oadd[0] = i - 1;
    OperatorAddLaderNormalOrder(&O, oadd, 2, 1.0, 0.0, &QDetNotOrderedFcFaBcBa);
    oadd[3] = CrAnMask + i;
    oadd[2] = i;
    oadd[1] = CrAnMask + i - 1;
    oadd[0] = i - 1;
    OperatorAddLaderNormalOrder(&O, oadd, 4, 1.0, 0, &QDetNotOrderedFcFaBcBa);
  }
  //OperatorAddLader(&O, oadd, 0, -20.0, 0);
  SetQComplexOperator(&O);
  PrintOperatorFullReadable(O);

  epsilon = 0.0000000001;
  zero = 0.0000000001;
  LanczosGroundState(O, &psi, &epsilon, zero, 0.0, &psi, 0);
  PrintWaveFunctionFullReadable(psi);


}


/*
 *
 */
void TestBlockListOfPoles()
{
  
  BlockListOfPoles variable;
  variable.QComplex=0;
  variable.N=10;
  variable.BlockSize=3;
  
  
  variable.VRe=malloc(variable.BlockSize*variable.BlockSize*variable.N*sizeof(double*));
  
  
  
  
  int iBlock,i,j;
  for(iBlock=0;iBlock<variable.N;iBlock++)
  {
    for(i=0;i<variable.BlockSize;i++)
    {
      for(j=0;j<variable.BlockSize;j++)
      {
        variable.VRe[iBlock * variable.BlockSize * variable.BlockSize + i * variable.BlockSize + j]= i + 100 * j + 10000 * iBlock;
        
        
        
      }
    }
  }
}


void TestMeanField()
{
  OperatorType RH, H,U,Sx,Sy,Sz,Ssqr,Lx,Ly,Lz,Lsqr,Jx,Jy,Jz,Jsqr,ldots;
  double VSx,VSy,VSz,VLx,VLy,VLz,VSsqr,VLsqr;
  OperatorType HMF, RHMF;
  MatrixType MatHMF, MatRHMF, EigenFunMatHMF, EigenFunMatRHMF;
  WaveFunctionType psitemp, psi, Rpsi, Rotpsi;
  double epsilon=0.00000001, zero=100*DBL_EPSILON, RhoDiag[14], EigenvalMatHMF[14], EigenvalMatRHMF[14];
  double energy, MFenergy, Renergy, RMFenergy;
  MatrixType Rho, EigenRho, RRho;
  uint32_t i,NFermion=6,L=1,Nelectron=2;
  uint16_t LindexUp[]={1,3,5,7,9,11,13};
  uint16_t LindexDn[]={0,2,4,6,8,10,12};
  double F[]={0,10.0,6.0,4.0};
  CreateOperatorCoulombSingleShell(NFermion, L, LindexUp, LindexDn, F, &U);
  CreateOperatorSx(NFermion, 2*L+1, LindexUp, LindexDn, &Sx);
  CreateOperatorSy(NFermion, 2*L+1, LindexUp, LindexDn, &Sy);
  CreateOperatorSz(NFermion, 2*L+1, LindexUp, LindexDn, &Sz);
  CreateOperatorSsqr(NFermion, 2*L+1, LindexUp, LindexDn, &Ssqr);
  CreateOperatorLx(NFermion, L, LindexUp, LindexDn, &Lx);
  CreateOperatorLy(NFermion, L, LindexUp, LindexDn, &Ly);
  CreateOperatorLz(NFermion, L, LindexUp, LindexDn, &Lz);
  CreateOperatorLsqr(NFermion, L, LindexUp, LindexDn, &Lsqr);
  CreateOperatorJx(NFermion, L, LindexUp, LindexDn, &Jx);
  CreateOperatorJy(NFermion, L, LindexUp, LindexDn, &Jy);
  CreateOperatorJz(NFermion, L, LindexUp, LindexDn, &Jz);
  CreateOperatorJsqr(NFermion, L, LindexUp, LindexDn, &Jsqr);
  CreateOperatorldots(NFermion, L, LindexUp, LindexDn, &ldots);
  OperatorCopy(&H, U);
  //OperatorAdd(&H, ldots, 0.01, 0.0);
  OperatorAdd(&H, Sz, 0.002, 0.0);
  OperatorAdd(&H, Lz, 0.001, 0.0);
  CreateStartWaveFunctionOnsiteEnergy(H, &psi, Nelectron, 0, 16);
  OperatorPsi(H, psi, &psitemp);
  WaveFunctionRandomize(&psitemp);
  WaveFunctionFree(&psi);
  OperatorPsi(H, psitemp, &psi);
  WaveFunctionRandomize(&psi);
  
  
  
  
  LanczosGroundState(H, &psi, &epsilon, zero, 0, &psi, 0);
  printf("\n================\n Operator \n================\n ");
  PrintOperatorReadable(H);
  printf("\n================\n Groundstate \n================\n ");
  PrintWaveFunctionReadable(psi);
  
  PsiHermitianOperatorPsi(Sx, psi, &VSx);
  PsiHermitianOperatorPsi(Sy, psi, &VSy);
  PsiHermitianOperatorPsi(Sz, psi, &VSz);
  PsiHermitianOperatorPsi(Ssqr, psi, &VSsqr);
  PsiHermitianOperatorPsi(Lx, psi, &VLx);
  PsiHermitianOperatorPsi(Ly, psi, &VLy);
  PsiHermitianOperatorPsi(Lz, psi, &VLz);
  PsiHermitianOperatorPsi(Lsqr, psi, &VLsqr);
  
  printf("Sx %22.15E\n",VSx);
  printf("Sy %22.15E\n",VSy);
  printf("Sz %22.15E\n",VSz);
  printf("Ssqr %22.15E\n",VSsqr);
  printf("Lx %22.15E\n",VLx);
  printf("Ly %22.15E\n",VLy);
  printf("Lz %22.15E\n",VLz);
  printf("Lsqr %22.15E\n",VLsqr);
  
  printf("\n================\n Density Matrix \n================\n ");
  CalculateDensityMatrix(&Rho, psi);
  Eigensystem(Rho, &EigenRho, RhoDiag);
  PrintMatrix(EigenRho);
  if(EigenRho.QComplex==0)
  {
    MatrixTranspose(&EigenRho);
  }
  else
  {
    MatrixConjugateTranspose(&EigenRho);
  }

  PrintMatrix(Rho);
  PrintMatrix(EigenRho);
  
  for(i=0;i<NFermion;i++)
  {
    printf("Density of orbital %i: %lf\n",i,RhoDiag[i]);
  }
  
  printf("\n================\n Mean field operator \n================\n ");
  
  MeanFieldOperator(&H, &HMF, Rho);
  PrintOperatorReadable(HMF);
  
  PsiHermitianOperatorPsi(H, psi, &energy);
  PsiHermitianOperatorPsi(HMF, psi, &MFenergy);
  
  printf("Energy    %22.15E\nMF Energy %22.15E\n",energy,MFenergy);
  
  
  printf("\n================\n Rotate psi \n================\n ");
  RotateWaveFunction(&Rotpsi, psi, EigenRho);
  PrintWaveFunctionReadable(Rotpsi);
  
  
  printf("\n================\n Rotate Operator \n================\n ");
  RotateOperatorNormalOrder(&RH, H, EigenRho, &QDetNotOrderedFcFaBcBa);
  OperatorChop(&RH, 0.000000001);
  PrintOperatorReadable(RH);
  
  
  printf("\n================\n Eigensystem of RH \n================\n ");
  WaveFunctionCopy(Rotpsi, &Rpsi);
  LanczosGroundState(RH, &Rpsi, &epsilon, zero, 0, &Rotpsi, 0);
  
  PrintWaveFunctionReadable(Rpsi);
  
  
  printf("\n================\n Density Matrix Rotated \n================\n ");
  CalculateDensityMatrix(&RRho, Rpsi);
  
  PrintMatrix(RRho);
  
  
  printf("\n================\n Mean field operator \n================\n ");
  
  MeanFieldOperator(&RH, &RHMF, RRho);
  PrintOperatorReadable(RHMF);
  
  PsiHermitianOperatorPsi(RH, Rpsi, &Renergy);
  PsiHermitianOperatorPsi(RHMF, Rpsi, &RMFenergy);
  
  printf("Energy    %22.15E\nMF Energy %22.15E\n",energy,MFenergy);
  printf("Energy    %22.15E\nMF Energy %22.15E\n",Renergy,RMFenergy);
  
  
  printf("\n================\n Mean field operator to one electron operator matrix \n================\n ");
  
  OperatorToOneElectronMatrix(&HMF, &MatHMF);
  PrintMatrix(MatHMF);
  OperatorToOneElectronMatrix(&RHMF, &MatRHMF);
  PrintMatrix(MatRHMF);
  
  
  printf("\n================\n Hartree Fock orbitals \n================\n ");
  
  Eigensystem(MatHMF, &EigenFunMatHMF, EigenvalMatHMF);
  PrintMatrix(EigenFunMatHMF);
  Eigensystem(MatRHMF, &EigenFunMatRHMF, EigenvalMatRHMF);
  PrintMatrix(EigenFunMatRHMF);
  
  for(i=0;i<MatHMF.Ni;i++)
  {
    printf("%22.15E %22.15E\n",EigenvalMatHMF[i],EigenvalMatRHMF[i]);
  }
  
  
  
}

void TestBlockLanczos()
{
  InitializeTiming(&GlobalTime);
  GlobalTimeOn=0xFFFFFFFF;
  BeginTiming("Total", &GlobalTime);
  uint32_t i,j;
  uint32_t Npsi=5;
  uint32_t NFermion=100;
  uint32_t N=NFermion*NFermion;
  OperatorType H;
  WaveFunctionType psi[Npsi];
  uint16_t CrIndex[N];
  uint16_t AnIndex[N];
  uint8_t *det;
  double tRe[N];
  double tIm[N];
  double Energy[Npsi];
  double EnergySqr[Npsi];
  double zero,epsilon;
  DeterminantRestrictionType DetRet;
  
  for(i=0;i<NFermion;i++)
  {
    for(j=0;j<NFermion;j++)
    {
      tIm[i+NFermion*j]=((double)i-(double)j);
      tRe[i+NFermion*j]=sqrt((i+j+1.0)/10.0)-(1.0/(1.+i+j));
      CrIndex[i+NFermion*j]=i;
      AnIndex[i+NFermion*j]=j;
    }
  }
  CreateOperatorCreateAnnihilateList(NFermion, CrIndex, AnIndex, tRe, tIm, N, &H);
  
  //PrintOperatorReadable(H);
  
  psi[0].NFermion=NFermion;
  psi[0].NBoson=0;
  psi[0].QComplex=0;
  psi[0].NBitsKey=16;
  WaveFunctionInit(&(psi[0]));
  det=calloc(psi[0].NWord,sizeof(uint8_t));
  for(i=0;i<NFermion;i++)
  {
    det[i>>3]=128>>(i&7);
    WaveFunctionAddElement(&(psi[0]), det, 1, 0);
    det[i>>3]=0;
  }
  for(i=1;i<Npsi;i++)
  {
    WaveFunctionInitCopyBasis(&(psi[i]), psi[0]);
  }
  
  BeginTiming("BlockLanczosGroundStateConserveBasisDense", &GlobalTime);
  BlockLanczosGroundStateConserveBasisDense(H, psi, Npsi);
  EndTiming("BlockLanczosGroundStateConserveBasisDense", &GlobalTime);
  
  
  BeginTiming("Energy", &GlobalTime);
  for(i=0;i<Npsi;i++)
  {
    PsiHermitianOperatorPsi(H, psi[i], &Energy[i]);
  }
  EndTiming("Energy", &GlobalTime);
  BeginTiming("Energy Sqr", &GlobalTime);
  for(i=0;i<Npsi;i++)
  {
    //PsiHermitianOperatorSqrPsi(H, psi[i], &EnergySqr[i]);
  }
  EndTiming("Energy Sqr", &GlobalTime);
  
  for(i=0;i<Npsi;i++)
  {
    //PrintWaveFunctionReadable(psi[i]);
    printf("Energy %22.15E, Variance %22.15E\n",Energy[i],EnergySqr[i]-Sqr(Energy[i]));
  }
  
  
  
  printf("\n\n\n\n\n================= Band Lanczos ================\n\n\n");
  
  WaveFunctionType *v;
  v=malloc((2*Npsi+100)*sizeof(WaveFunctionType));
  
  
  for(i=0;i<Npsi;i++)
  {
   WaveFunctionRandomize(&psi[i]);
    v[i]=psi[i];
  }
  for(i=Npsi;i<2*Npsi+100;i++)
  {
    WaveFunctionInitCopyBasis(&v[i], v[0]);
  }
  
  BlockLanczosGroundStateConserveBasisKrylov(H, v, Npsi, DBL_EPSILON);
  
  
  printf("\n\n\n\n\n================= BlockGroundState ================\n\n\n");
  
  for(i=0;i<Npsi;i++)
  {
    WaveFunctionRandomize(&v[i]);
    psi[i]=v[i];
  }
  
  
  zero=100*DBL_EPSILON;
  epsilon=DBL_EPSILON;
  DetRet.N=0;
  BlockGroundState(H, Npsi, zero, &epsilon, DetRet, DetRet, psi, 0);
  
  EndTiming("Total", &GlobalTime);
  PrintTiming(GlobalTime);
}