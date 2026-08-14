/*  Given a BitMask String returns an index array
 *  Bitmasks consist of '0', '1' and other charaters
 *  Other characters are ignored.
 */
int BitMaskToIndex(const char*BitMask, int *NF, int **Index)
{
  int i,iIndex=0;
  int Len = (int) strlen(BitMask);
  *Index = (int*)malloc(Len * sizeof(int));
  if(*Index==NULL) return 1;
  *NF=0;
  for(i=0;i<Len;i++)
  {
    (*Index)[i]=-1;
    switch(BitMask[i])
    {
      case '1':
      {
        (*Index)[iIndex] = *NF;
        *NF=*NF+1;
        iIndex++;
        break;
      }
      case '0':
      {
        *NF=*NF+1;
        break;
      }
    }
  }
  return 0;
}



int IndexToBitMask(int NF, int *Index, int len, char**BitMask)
{
  int i;
  (*BitMask) = (char*)malloc((1+(NF))*sizeof(char));
  (*BitMask)[NF]='\0';
  for(i=0;i<NF;i++)
  {
    (*BitMask)[i]='0';
  }
  for(i=0;i<len;i++)
  {
    (*BitMask)[Index[i]]='1';
  }
  return 0;
}


/* Looks up Name in PhysicalConstants (SI units) and returns its value.
 * Returns 0 and sets *Value if found, 1 (leaving *Value untouched)
 * otherwise. Uses sizeof(PhysicalConstants) rather than a stored
 * count, so this stays correct as entries are added, removed, or
 * reordered. */
int GetPhysicalConstant(const char *Name, double *Value)
{
  size_t N = sizeof(PhysicalConstants)/sizeof(PhysicalConstants[0]);
  size_t i;
  for(i=0;i<N;i++)
    if(strcmp(PhysicalConstants[i].name,Name)==0)
    {
      *Value = PhysicalConstants[i].value;
      return 0;
    }
  return 1;
}


/* Looks up Name in EnergyUnits (matching either namelong or
 * nameshort, e.g. "electron Volt" or "eV") and returns its value in
 * Joule. Returns 0 and sets *Value if found, 1 (leaving *Value
 * untouched) otherwise. Uses sizeof(EnergyUnits) rather than a stored
 * count, so this stays correct as entries are added, removed, or
 * reordered. */
int GetEnergyUnit(const char *Name, double *Value)
{
  size_t N = sizeof(EnergyUnits)/sizeof(EnergyUnits[0]);
  size_t i;
  for(i=0;i<N;i++)
    if(strcmp(EnergyUnits[i].namelong,Name)==0 || strcmp(EnergyUnits[i].nameshort,Name)==0)
    {
      *Value = EnergyUnits[i].value;
      return 0;
    }
  return 1;
}
