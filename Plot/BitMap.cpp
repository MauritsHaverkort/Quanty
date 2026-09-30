// writes a little endian 32 bit number at str
static void BitMapPutUint32(unsigned char *str, uint32_t v)
{
  str[0] = (unsigned char)(v);
  str[1] = (unsigned char)(v >> 8);
  str[2] = (unsigned char)(v >> 16);
  str[3] = (unsigned char)(v >> 24);
}

// A 32 bit BMP with a BITMAPV5HEADER and BI_BITFIELDS, so that viewers use the alpha channel.
// With the older 40 byte header Firefox and WebKit ignore it; a BITMAPV4HEADER WebKit does not read.
int BitMapToBinaryString(BitMapType *BitMap, DynamicString *A)
{
  int i,j,k;
  unsigned char *str;
  int const Offset = 14 + 124; // file header and BITMAPV5HEADER
  A->length = Offset + 4 * BitMap->NpointsX * BitMap->NpointsY;
  A->allocated = A->length+1;
  A->str = (char*) calloc(A->allocated, sizeof(char));
  str = (unsigned char*)A->str;
  // file header
  str[0]='B';
  str[1]='M';
  BitMapPutUint32(str+2, (uint32_t)A->length);
  BitMapPutUint32(str+10, (uint32_t)Offset);
  // BITMAPV5HEADER, fields that are not set are 0
  BitMapPutUint32(str+14, 124);                        // size of the header
  BitMapPutUint32(str+18, (uint32_t)BitMap->NpointsX); // width
  BitMapPutUint32(str+22, (uint32_t)BitMap->NpointsY); // height, positive: the bottom row comes first
  str[26]=1;                                           // planes
  str[28]=32;                                          // bits per pixel
  BitMapPutUint32(str+30, 3);                          // BI_BITFIELDS
  BitMapPutUint32(str+34, (uint32_t)(4 * BitMap->NpointsX * BitMap->NpointsY)); // size of the pixel data
  BitMapPutUint32(str+54, 0x00FF0000);                 // red mask
  BitMapPutUint32(str+58, 0x0000FF00);                 // green mask
  BitMapPutUint32(str+62, 0x000000FF);                 // blue mask
  BitMapPutUint32(str+66, 0xFF000000);                 // alpha mask
  BitMapPutUint32(str+70, 0x73524742);                 // colour space 'sRGB'
  BitMapPutUint32(str+122, 4);                         // rendering intent LCS_GM_IMAGES
  k=Offset;
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
      str[k]=(unsigned char)( 255 - ColorTransparency(BitMap->Pixels[j+i*BitMap->NpointsX]) ); // the BMP stores the opacity
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
    BitMap->Pixels[i] = 0xFFFFFFFF; // white, fully transparent
  }
  
  return 0;
}
