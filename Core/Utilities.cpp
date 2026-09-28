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


/*  Combines two BitMask Strings bit by bit into a new BitMask String
 *  Bitmasks consist of '0', '1' and other charaters
 *  Other characters are ignored, i.e. "1100 0011" and "11000011" are the same bitmask
 *  Operation is 'a' (and), 'o' (or) or 'x' (xor)
 *  NF the number of bits of the resulting bitmask
 *  Returns 0 on success, 1 if the allocation failed, 2 if the two bitmasks
 *  contain a different number of bits and 3 for an unknown operation
 */
int BitMaskCombine(const char *BitMask1, const char *BitMask2, char Operation, int *NF, char **BitMask)
{
  int i1=0,i2=0,b1,b2,b=0;
  *NF=0;
  *BitMask = (char*)malloc((1+strlen(BitMask1))*sizeof(char));
  if(*BitMask==NULL) return 1;
  while(1)
  {
    while(BitMask1[i1]!='\0' && BitMask1[i1]!='0' && BitMask1[i1]!='1') i1++;
    while(BitMask2[i2]!='\0' && BitMask2[i2]!='0' && BitMask2[i2]!='1') i2++;
    if(BitMask1[i1]=='\0' || BitMask2[i2]=='\0') break;
    b1 = (BitMask1[i1]=='1');
    b2 = (BitMask2[i2]=='1');
    switch(Operation)
    {
      case 'a': b = (b1 && b2); break;
      case 'o': b = (b1 || b2); break;
      case 'x': b = (b1 != b2); break;
      default:
      {
        free(*BitMask);
        *BitMask=NULL;
        return 3;
      }
    }
    if(b) (*BitMask)[*NF]='1'; else (*BitMask)[*NF]='0';
    *NF=*NF+1;
    i1++;
    i2++;
  }
  (*BitMask)[*NF]='\0';
  /* one of the two strings ended before the other, i.e. they have a different number of bits */
  if(BitMask1[i1]!='\0' || BitMask2[i2]!='\0')
  {
    free(*BitMask);
    *BitMask=NULL;
    return 2;
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
