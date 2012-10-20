uint32_t file_exists(const char * filename)
{
  FILE *stream;
  if((stream=fopen(filename, "r")))
  {
    fclose(stream);
    return 1;
  }
  return 0;
}


uint32_t FileReadWaveFunctionPlain(WaveFunctionType * O, const char *const flname)
{
  FILE *stream;
  uint32_t retval;
  if((stream=fopen(flname,"r"))==NULL)
  {
    printf("ERROR\nCan't open %s for reading.\n",flname);
    return 1;
  }
  retval=StreamReadWaveFunctionPlain(O,stream);
  fclose(stream);
  return retval;
}

uint32_t FileReadWaveFunctionPlainBinary(WaveFunctionType * O, const char *const flname)
{
  FILE *stream;
  uint32_t retval;
  if((stream=fopen(flname,"r"))==NULL)
  {
    printf("ERROR\nCan't open %s for reading.\n",flname);
    return 1;
  }
  retval=StreamReadWaveFunctionPlainBinary(O,stream);
  fclose(stream);
  return retval;
}

uint32_t FileReadMatrixType(MatrixType * O, const char *const flname)
{
  FILE *stream;
  uint32_t retval;
  if((stream=fopen(flname,"r"))==NULL)
  {
    printf("ERROR\nCan't open %s for reading.\n",flname);
    return 1;
  }
  retval=StreamReadMatrixType(O,stream);
  fclose(stream);
  return retval;
}

uint32_t FileReadAndersonMatrixType(AndersonMatrixType * O, const char *const flname)
{
  FILE *stream;
  uint32_t retval;
  if((stream=fopen(flname,"r"))==NULL)
  {
    printf("ERROR\nCan't open %s for reading.\n",flname);
    return 1;
  }
  retval=StreamReadAndersonMatrixType(O,stream);
  fclose(stream);
  return retval;
}


uint32_t FileReadPlainTriDiagonalMatrixType(TriDiagonalMatrixType * O, const char *const flname)
{
  FILE *stream;
  uint32_t retval;
  if((stream=fopen(flname,"r"))==NULL)
  {
    printf("ERROR\nCan't open %s for reading.\n",flname);
    return 1;
  }
  retval=StreamReadPlainTriDiagonalMatrixType(O,stream);
  fclose(stream);
  return retval;
}


uint32_t FileReadAtom(AtomType * O, const char *const flname)
{
  FILE *stream;
  uint32_t retval;
  if((stream=fopen(flname,"r"))==NULL)
  {
    printf("ERROR\nCan't open %s for reading.\n",flname);
    return 1;
  }
  retval=StreamReadAtom(O,stream);
  fclose(stream);
  return retval;
}

uint32_t FileReadCrystalStructure(CrystalStructureType * O, const char *const flname)
{
  FILE *stream;
  uint32_t retval;
  if((stream=fopen(flname,"r"))==NULL)
  {
    printf("ERROR\nCan't open %s for reading.\n",flname);
    return 1;
  }
  retval=StreamReadCrystalStructure(O,stream);
  fclose(stream);
  return retval;
}

uint32_t FileReadCrystalStructureSuperCell(CrystalStructureSuperCellType * O, const char *const flname)
{
  FILE *stream;
  uint32_t retval;
  if((stream=fopen(flname,"r"))==NULL)
  {
    printf("ERROR\nCan't open %s for reading.\n",flname);
    return 1;
  }
  retval=StreamReadCrystalStructureSuperCell(O,stream);
  fclose(stream);
  return retval;
}

uint32_t FileReadHeisenbergModelSettings(HeisenbergModelSettings * O, const char *const flname)
{
  FILE *stream;
  uint32_t retval;
  if((stream=fopen(flname,"r"))==NULL)
  {
    printf("ERROR\nCan't open %s for reading.\n",flname);
    return 1;
  }
  retval=StreamReadHeisenbergModelSettings(O,stream);
  fclose(stream);
  return retval;
}

uint32_t FileReadDimerizedSpinChainSettings(DimerizedSpinChainSettings * O, const char *const flname)
{
  FILE *stream;
  uint32_t retval;
  if((stream=fopen(flname,"r"))==NULL)
  {
    printf("ERROR\nCan't open %s for reading.\n",flname);
    return 1;
  }
  retval=StreamReadDimerizedSpinChainSettings(O,stream);
  fclose(stream);
  return retval;
}

uint32_t FileReadtJModelSettings(tJModelSettings * O, const char *const flname)
{
  FILE *stream;
  uint32_t retval;
  if((stream=fopen(flname,"r"))==NULL)
  {
    printf("ERROR\nCan't open %s for reading.\n",flname);
    return 1;
  }
  retval=StreamReadtJModelSettings(O,stream);
  fclose(stream);
  return retval;
}

uint32_t FileReadSpinLaderSettings(SpinLaderSettings * O, const char *const flname)
{
  FILE *stream;
  uint32_t retval;
  if((stream=fopen(flname,"r"))==NULL)
  {
    printf("ERROR\nCan't open %s for reading.\n",flname);
    return 1;
  }
  retval=StreamReadSpinLaderSettings(O,stream);
  fclose(stream);
  return retval;
}

uint32_t FileReadHubbardModelSettings(HubbardModelSettings * O, const char *const flname)
{
  FILE *stream;
  uint32_t retval;
  if((stream=fopen(flname,"r"))==NULL)
  {
    printf("ERROR\nCan't open %s for reading.\n",flname);
    return 1;
  }
  retval=StreamReadHubbardModelSettings(O,stream);
  fclose(stream);
  return retval;
}

uint32_t FileReadDMFTXASSettings(DMFTXASSettings * O, const char *const flname)
{
  FILE *stream;
  uint32_t retval;
  if((stream=fopen(flname,"r"))==NULL)
  {
    printf("ERROR\nCan't open %s for reading.\n",flname);
    return 1;
  }
  retval=StreamReadDMFTXASSettings(O,stream);
  fclose(stream);
  return retval;
}

uint32_t FileReadDMFTSettings(DMFTSettings * O, const char *const flname)
{
  FILE *stream;
  uint32_t retval;
  if((stream=fopen(flname,"r"))==NULL)
  {
    printf("ERROR\nCan't open %s for reading.\n",flname);
    return 1;
  }
  retval=StreamReadDMFTSettings(O,stream);
  fclose(stream);
  return retval;
}


