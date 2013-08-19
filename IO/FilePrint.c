
uint32_t FilePrintSpectra(SpectraType *const Spec, const uint32_t NSpec, const char *const flname, const char *const mode)
{
  FILE *stream;
  uint32_t retval;
  if((strcmp(mode,"w")&&strcmp(mode,"a"))!=0)
  {
    printf("ERROR in FilePrintSpectra.\nmode must be either \"w\" for overwrite or \"a\" for append\n %s as input is not excepted\n",mode);
    return 1;
  }
  if((stream=fopen(flname,mode))==NULL)
  {
    printf("ERROR\nCan't open %s for writing with mode %s.\n",flname,mode);
    return 1;
  }
  retval=StreamPrintSpectra(Spec,NSpec,stream);
  fclose(stream);
  return retval;
}


uint32_t FilePrintWaveFunctionPlain(WaveFunctionType psi, const char *const flname, const char *const mode)
{
  FILE *stream;
  uint32_t retval;
  if((strcmp(mode,"w")&&strcmp(mode,"a"))!=0)
  {
    printf("ERROR in FilePrintWaveFunctionPlain.\nmode must be either \"w\" for overwrite or \"a\" for append\n %s as input is not excepted\n",mode);
    return 1;
  }
  if((stream=fopen(flname,mode))==NULL)
  {
    printf("ERROR\nCan't open %s for writing with mode %s.\n",flname,mode);
    return 1;
  }
  retval=StreamPrintWaveFunctionPlain(psi,stream);
  fclose(stream);
  return retval;
}


uint32_t FilePrintWaveFunctionPlainBinary(WaveFunctionType psi, const char *const flname, const char *const mode)
{
  FILE *stream;
  uint32_t retval;
  if((strcmp(mode,"w")&&strcmp(mode,"a"))!=0)
  {
    printf("ERROR in FilePrintWaveFunctionPlainBinary.\nmode must be either \"w\" for overwrite or \"a\" for append\n %s as input is not excepted\n",mode);
    return 1;
  }
  if((stream=fopen(flname,mode))==NULL)
  {
    printf("ERROR\nCan't open %s for writing with mode %s.\n",flname,mode);
    return 1;
  }
  retval=StreamPrintWaveFunctionPlainBinary(psi,stream);
  fclose(stream);
  return retval;
}

uint32_t FilePrintWaveFunctionReadable(WaveFunctionType psi, const char *const flname, const char *const mode)
{
  FILE *stream;
  uint32_t retval;
  if((strcmp(mode,"w")&&strcmp(mode,"a"))!=0)
  {
    printf("ERROR in FilePrintWaveFunctionReadable.\nmode must be either \"w\" for overwrite or \"a\" for append\n %s as input is not excepted\n",mode);
    return 1;
  }
  if((stream=fopen(flname,mode))==NULL)
  {
    printf("ERROR\nCan't open %s for writing with mode %s.\n",flname,mode);
    return 1;
  }
  retval=StreamPrintWaveFunctionReadable(psi,stream);
  fclose(stream);
  return retval;
}


uint32_t FilePrintWaveFunctionReadableTruncated(WaveFunctionType psi, double epsilon, const char *const flname, const char *const mode)
{
  FILE *stream;
  uint32_t retval;
  if((strcmp(mode,"w")&&strcmp(mode,"a"))!=0)
  {
    printf("ERROR in FilePrintWaveFunctionReadable.\nmode must be either \"w\" for overwrite or \"a\" for append\n %s as input is not excepted\n",mode);
    return 1;
  }
  if((stream=fopen(flname,mode))==NULL)
  {
    printf("ERROR\nCan't open %s for writing with mode %s.\n",flname,mode);
    return 1;
  }
  retval=StreamPrintWaveFunctionReadableTruncated(psi,epsilon,stream);
  fclose(stream);
  return retval;
}

uint32_t FilePrintWaveFunctionFullReadable(WaveFunctionType psi, const char *const flname, const char *const mode)
{
  FILE *stream;
  uint32_t retval;
  if((strcmp(mode,"w")&&strcmp(mode,"a"))!=0)
  {
    printf("ERROR in FilePrintWaveFunctionFullReadable.\nmode must be either \"w\" for overwrite or \"a\" for append\n %s as input is not excepted\n",mode);
    return 1;
  }
  if((stream=fopen(flname,mode))==NULL)
  {
    printf("ERROR\nCan't open %s for writing with mode %s.\n",flname,mode);
    return 1;
  }
  retval=StreamPrintWaveFunctionFullReadable(psi,stream);
  fclose(stream);
  return retval;
}


uint32_t FilePrintOperatorReadable(OperatorType O, const char *const flname, const char *const mode)
{
  FILE *stream;
  uint32_t retval;
  if((strcmp(mode,"w")&&strcmp(mode,"a"))!=0)
  {
    printf("ERROR in FilePrintOperatorReadable.\nmode must be either \"w\" for overwrite or \"a\" for append\n %s as input is not excepted\n",mode);
    return 1;
  }
  if((stream=fopen(flname,mode))==NULL)
  {
    printf("ERROR\nCan't open %s for writing with mode %s.\n",flname,mode);
    return 1;
  }
  retval=StreamPrintOperatorReadable(O,stream);
  fclose(stream);
  return retval;
}

uint32_t FilePrintOperatorFullReadable(OperatorType O, const char *const flname, const char *const mode)
{
  FILE *stream;
  uint32_t retval;
  if((strcmp(mode,"w")&&strcmp(mode,"a"))!=0)
  {
    printf("ERROR in FilePrintOperatorFullReadable.\nmode must be either \"w\" for overwrite or \"a\" for append\n %s as input is not excepted\n",mode);
    return 1;
  }
  if((stream=fopen(flname,mode))==NULL)
  {
    printf("ERROR\nCan't open %s for writing with mode %s.\n",flname,mode);
    return 1;
  }
  retval=StreamPrintOperatorFullReadable(O,stream);
  fclose(stream);
  return retval;
}

uint32_t FilePrintCrystalStructure(CrystalStructureType O, const char *const flname, const char *const mode)
{
  FILE *stream;
  uint32_t retval;
  if((strcmp(mode,"w")&&strcmp(mode,"a"))!=0)
  {
    printf("ERROR in FilePrintCrystalStructure.\nmode must be either \"w\" for overwrite or \"a\" for append\n %s as input is not excepted\n",mode);
    return 1;
  }
  if((stream=fopen(flname,mode))==NULL)
  {
    printf("ERROR\nCan't open %s for writing with mode %s.\n",flname,mode);
    return 1;
  }
  retval=StreamPrintCrystalStructure(O,stream);
  fclose(stream);
  return retval;
}

uint32_t FilePrintCrystalStructureSuperCell(CrystalStructureSuperCellType O, const char *const flname, const char *const mode)
{
  FILE *stream;
  uint32_t retval;
  if((strcmp(mode,"w")&&strcmp(mode,"a"))!=0)
  {
    printf("ERROR in FilePrintCrystalStructureSuperCell.\nmode must be either \"w\" for overwrite or \"a\" for append\n %s as input is not excepted\n",mode);
    return 1;
  }
  if((stream=fopen(flname,mode))==NULL)
  {
    printf("ERROR\nCan't open %s for writing with mode %s.\n",flname,mode);
    return 1;
  }
  retval=StreamPrintCrystalStructureSuperCell(O,stream);
  fclose(stream);
  return retval;
}

uint32_t FilePrintMatrix(MatrixType M, const char *const flname, const char *const mode)
{
  FILE *stream;
  uint32_t retval;
  if((strcmp(mode,"w")&&strcmp(mode,"a"))!=0)
  {
    printf("ERROR in FilePrintMatrix.\nmode must be either \"w\" for overwrite or \"a\" for append\n %s as input is not excepted\n",mode);
    return 1;
  }
  if((stream=fopen(flname,mode))==NULL)
  {
    printf("ERROR\nCan't open %s for writing with mode %s.\n",flname,mode);
    return 1;
  }
  retval=StreamPrintMatrix(M,stream);
  fclose(stream);
  return retval;
}


uint32_t FilePrintPlainMatrix(MatrixType M, const char *const flname, const char *const mode)
{
  FILE *stream;
  uint32_t retval;
  if((strcmp(mode,"w")&&strcmp(mode,"a"))!=0)
  {
    printf("ERROR in FilePrintPlainMatrix.\nmode must be either \"w\" for overwrite or \"a\" for append\n %s as input is not excepted\n",mode);
    return 1;
  }
  if((stream=fopen(flname,mode))==NULL)
  {
    printf("ERROR\nCan't open %s for writing with mode %s.\n",flname,mode);
    return 1;
  }
  retval=StreamPrintPlainMatrix(M,stream);
  fclose(stream);
  return retval;
}

uint32_t FilePrintPlainTriDiagonalMatrix(TriDiagonalMatrixType M, const char *const flname, const char *const mode)
{
  FILE *stream;
  uint32_t retval;
  if((strcmp(mode,"w")&&strcmp(mode,"a"))!=0)
  {
    printf("ERROR in FilePrintPlainTriDiagonalMatrix.\nmode must be either \"w\" for overwrite or \"a\" for append\n %s as input is not excepted\n",mode);
    return 1;
  }
  if((stream=fopen(flname,mode))==NULL)
  {
    printf("ERROR\nCan't open %s for writing with mode %s.\n",flname,mode);
    return 1;
  }
  retval=StreamPrintPlainTriDiagonalMatrix(M,stream);
  fclose(stream);
  return retval;
}

uint32_t FilePrintTightBindingDef(TightBindingDef O, const char *const flname, const char *const mode)
{
  FILE *stream;
  uint32_t retval;
  if((strcmp(mode,"w")&&strcmp(mode,"a"))!=0)
  {
    printf("ERROR in FilePrintTightBindingDef.\nmode must be either \"w\" for overwrite or \"a\" for append\n %s as input is not excepted\n",mode);
    return 1;
  }
  if((stream=fopen(flname,mode))==NULL)
  {
    printf("ERROR\nCan't open %s for writing with mode %s.\n",flname,mode);
    return 1;
  }
  retval=StreamPrintTightBindingDef(O,stream);
  fclose(stream);
  return retval;
}


uint32_t FilePrintDimerizedSpinChainSettings(DimerizedSpinChainSettings O, const char *const flname, const char *const mode)
{
  FILE *stream;
  uint32_t retval;
  if((strcmp(mode,"w")&&strcmp(mode,"a"))!=0)
  {
    printf("ERROR in FilePrintDimerizedSpinChainSettings.\nmode must be either \"w\" for overwrite or \"a\" for append\n %s as input is not excepted\n",mode);
    return 1;
  }
  if((stream=fopen(flname,mode))==NULL)
  {
    printf("ERROR\nCan't open %s for writing with mode %s.\n",flname,mode);
    return 1;
  }
  retval=StreamPrintDimerizedSpinChainSettings(O,stream);
  fclose(stream);
  return retval;
}



uint32_t FilePrintHeisenbergModelSettings(HeisenbergModelSettings O, const char *const flname, const char *const mode)
{
  FILE *stream;
  uint32_t retval;
  if((strcmp(mode,"w")&&strcmp(mode,"a"))!=0)
  {
    printf("ERROR in FilePrintHeisenbergModelSettings.\nmode must be either \"w\" for overwrite or \"a\" for append\n %s as input is not excepted\n",mode);
    return 1;
  }
  if((stream=fopen(flname,mode))==NULL)
  {
    printf("ERROR\nCan't open %s for writing with mode %s.\n",flname,mode);
    return 1;
  }
  retval=StreamPrintHeisenbergModelSettings(O,stream);
  fclose(stream);
  return retval;
}


uint32_t FilePrinttJModelSettings(tJModelSettings O, const char *const flname, const char *const mode)
{
  FILE *stream;
  uint32_t retval;
  if((strcmp(mode,"w")&&strcmp(mode,"a"))!=0)
  {
    printf("ERROR in FilePrinttJModelSettings.\nmode must be either \"w\" for overwrite or \"a\" for append\n %s as input is not excepted\n",mode);
    return 1;
  }
  if((stream=fopen(flname,mode))==NULL)
  {
    printf("ERROR\nCan't open %s for writing with mode %s.\n",flname,mode);
    return 1;
  }
  retval=StreamPrinttJModelSettings(O,stream);
  fclose(stream);
  return retval;
}

uint32_t FilePrintSpinLaderSettings(SpinLaderSettings O, const char *const flname, const char *const mode)
{
  FILE *stream;
  uint32_t retval;
  if((strcmp(mode,"w")&&strcmp(mode,"a"))!=0)
  {
    printf("ERROR in FilePrintSpinLaderSettings.\nmode must be either \"w\" for overwrite or \"a\" for append\n %s as input is not excepted\n",mode);
    return 1;
  }
  if((stream=fopen(flname,mode))==NULL)
  {
    printf("ERROR\nCan't open %s for writing with mode %s.\n",flname,mode);
    return 1;
  }
  retval=StreamPrintSpinLaderSettings(O,stream);
  fclose(stream);
  return retval;
}

uint32_t FilePrintHubbardModelSettings(HubbardModelSettings O, const char *const flname, const char *const mode)
{
  FILE *stream;
  uint32_t retval;
  if((strcmp(mode,"w")&&strcmp(mode,"a"))!=0)
  {
    printf("ERROR in FilePrintHubbardModelSettings.\nmode must be either \"w\" for overwrite or \"a\" for append\n %s as input is not excepted\n",mode);
    return 1;
  }
  if((stream=fopen(flname,mode))==NULL)
  {
    printf("ERROR\nCan't open %s for writing with mode %s.\n",flname,mode);
    return 1;
  }
  retval=StreamPrintHubbardModelSettings(O,stream);
  fclose(stream);
  return retval;
}

uint32_t FilePrintDMFTXASSettings(DMFTXASSettings O, const char *const flname, const char *const mode)
{
  FILE *stream;
  uint32_t retval;
  if((strcmp(mode,"w")&&strcmp(mode,"a"))!=0)
  {
    printf("ERROR in FilePrintDMFTXASSettings.\nmode must be either \"w\" for overwrite or \"a\" for append\n %s as input is not excepted\n",mode);
    return 1;
  }
  if((stream=fopen(flname,mode))==NULL)
  {
    printf("ERROR\nCan't open %s for writing with mode %s.\n",flname,mode);
    return 1;
  }
  retval=StreamPrintDMFTXASSettings(O,stream);
  fclose(stream);
  return retval;
}

uint32_t FilePrintDMFTSettings(DMFTSettings O, const char *const flname, const char *const mode)
{
  FILE *stream;
  uint32_t retval;
  if((strcmp(mode,"w")&&strcmp(mode,"a"))!=0)
  {
    printf("ERROR in FilePrintDMFTSettings.\nmode must be either \"w\" for overwrite or \"a\" for append\n %s as input is not excepted\n",mode);
    return 1;
  }
  if((stream=fopen(flname,mode))==NULL)
  {
    printf("ERROR\nCan't open %s for writing with mode %s.\n",flname,mode);
    return 1;
  }
  retval=StreamPrintDMFTSettings(O,stream);
  fclose(stream);
  return retval;
}


uint32_t FilePrintSparseMatrix(SparseMatrixType M, const char *const flname, const char *const mode)
{
  FILE *stream;
  uint32_t retval;
  if((strcmp(mode,"w")&&strcmp(mode,"a"))!=0)
  {
    printf("ERROR in FilePrintSparseMatrix.\nmode must be either \"w\" for overwrite or \"a\" for append\n %s as input is not excepted\n",mode);
    return 1;
  }
  if((stream=fopen(flname,mode))==NULL)
  {
    printf("ERROR\nCan't open %s for writing with mode %s.\n",flname,mode);
    return 1;
  }
  retval=StreamPrintSparseMatrix(M,stream);
  fclose(stream);
  return retval;
}


uint32_t FilePrintAndersonMatrix(AndersonMatrixType M, const char *const flname, const char *const mode)
{
  FILE *stream;
  uint32_t retval;
  if((strcmp(mode,"w")&&strcmp(mode,"a"))!=0)
  {
    printf("ERROR in FilePrintAndersonMatrix.\nmode must be either \"w\" for overwrite or \"a\" for append\n %s as input is not excepted\n",mode);
    return 1;
  }
  if((stream=fopen(flname,mode))==NULL)
  {
    printf("ERROR\nCan't open %s for writing with mode %s.\n",flname,mode);
    return 1;
  }
  retval=StreamPrintAndersonMatrix(M,stream);
  fclose(stream);
  return retval;
}