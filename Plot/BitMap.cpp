int BitMapToBinaryString(BitMapType *BitMap, DynamicString *A)
{
  int i,j,k;
  unsigned char *str;
  A->length = 54 + 4 * BitMap->NpointsX * BitMap->NpointsY;
  A->allocated = A->length+1;
  A->str = (char*) calloc(A->allocated, sizeof(char));
  str = (unsigned char*)A->str;
  str[0]='B';
  str[1]='M';
  str[2] = (unsigned char)(A->length);
  str[3] = (unsigned char)(A->length >> 8);
  str[4] = (unsigned char)(A->length >> 16);
  str[5] = (unsigned char)(A->length >> 24);
  str[10]=54;
  str[14]=40;
  str[18] = (unsigned char)(BitMap->NpointsX);
  str[19] = (unsigned char)(BitMap->NpointsX >> 8);
  str[20] = (unsigned char)(BitMap->NpointsX >> 16);
  str[21] = (unsigned char)(BitMap->NpointsX >> 24);
  str[22] = (unsigned char)(BitMap->NpointsY);
  str[23] = (unsigned char)(BitMap->NpointsY >> 8);
  str[24] = (unsigned char)(BitMap->NpointsY >> 16);
  str[25] = (unsigned char)(BitMap->NpointsY >> 24);
  str[26]=1;
  str[28]=32;
  k=54;
  for(i=BitMap->NpointsY-1;i>=0;i--)
  {
    for(j=0;j<BitMap->NpointsX;j++)
    {
      str[k]=(unsigned char)( BitMap->Pixels[j+i*BitMap->NpointsX] );
      k++;
      str[k]=(unsigned char)( BitMap->Pixels[j+i*BitMap->NpointsX]>>8 );
      k++;
      str[k]=(unsigned char)( BitMap->Pixels[j+i*BitMap->NpointsX]>>16 );
      k++;
      str[k]=(unsigned char)( BitMap->Pixels[j+i*BitMap->NpointsX]>>24 );
      k++;
    }
  }
  return 0;
}

int BitMapToBase64String(BitMapType *BitMap, DynamicString *A64)
{
  DynamicString A;
  BitMapToBinaryString(BitMap, &A);
  EncodeStringBase64(&A, A64);
  FreeDynamicString(&A);
  return 0;
}

int BitMapToFile(BitMapType *BitMap, const char *filename)
{
  FILE *fp;
  DynamicString A;
  BitMapToBinaryString(BitMap, &A);
  fp = fopen( filename , "wb" );
  fwrite(A.str , 1 , A.length , fp );
  fclose(fp);
  FreeDynamicString(&A);
  return 0;
}

int MakeBitMatWhiteTransparent(BitMapType *BitMap)
{
  for(int i=0;i<BitMap->NpointsX*BitMap->NpointsY;i++)
  {
    BitMap->Pixels[i] = 0x00FFFFFF;
  }
  
  return 0;
}
