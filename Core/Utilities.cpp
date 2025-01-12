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
