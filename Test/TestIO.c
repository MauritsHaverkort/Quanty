


void TestBMP()
{
  int    i, j;
  MatrixType A;
  char   flname[256];
  A.Ni = 300;
  A.Nj = 500;
  A.QComplex = 0;
  InitMatrix(&A);
  for (i = 0; i < A.Ni; i++)
  {
    for (j = 0; j < A.Nj; j++)
    {
      A.VRe[i][j] = (double)(i * j) / ((double)(A.Ni * A.Nj));
    }
  }
  strcpy(flname, "test.bmp");
  FilePrintRealMatrixToBitmap(A, flname, &RealToRGBColorBW);
}

void TestWaveIO()
{
  WaveFunctionType psi, psiNew;
  uint8_t det[20];
  uint32_t i;
  psi.NBitsKey = 16;
  psi.NFermion = 10;
  WaveFunctionInit(&psi);
  strcpy(psi.Name, "Name of Wave Function");
  for (i = 1; i < 10; i++)
  {
    det[0] = (uint8_t) i;
    det[1] = (uint8_t) (i % 2);
    RealWaveFunctionAddElement(&psi, det, 0.5 / ((double)i));
  }
  FilePrintWaveFunctionPlain(psi, "temp", "w");
  psiNew.NBitsKey = 16;
  FileReadWaveFunctionPlain(&psiNew, "temp");
  PrintWaveFunctionPlain(psi);
  PrintWaveFunctionPlain(psiNew);
}
