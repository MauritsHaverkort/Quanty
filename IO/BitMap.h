/* This file was automatically generated.  Do not edit! */
uint32_t PrintReal3MatrixToBitmapRescale(MatrixType *M,void(*RealToRGBColor)(double,unsigned char *),uint32_t const NW,uint32_t const NH);
int FilePrintReal3MatrixToBitmapRescale(MatrixType *M,const char *const flname,void(*RealToRGBColor)(double,unsigned char *),uint32_t const NW,uint32_t const NH);
int StreamPrintReal3MatrixToBitmapRescale(MatrixType *M,FILE *stream,void(*RealToRGBColor)(double,unsigned char *),uint32_t const NW,uint32_t const NH);
uint32_t PrintRealMatrixToBitmapRescale(MatrixType M,void(*RealToRGBColor)(double,unsigned char *),uint32_t const NW,uint32_t const NH);
int FilePrintRealMatrixToBitmapRescale(MatrixType M,const char *const flname,void(*RealToRGBColor)(double,unsigned char *),uint32_t const NW,uint32_t const NH);
int StreamPrintRealMatrixToBitmapRescale(MatrixType M,FILE *stream,void(*RealToRGBColor)(double,unsigned char *),uint32_t const NW,uint32_t const NH);
uint32_t PrintReal3MatrixToBitmap(MatrixType *M,void(*Real3ToRGBColor)(double *,unsigned char *));
int FilePrintReal3MatrixToBitmap(MatrixType *M,const char *const flname,void(*Real3ToRGBColor)(double *,unsigned char *));
int StreamPrintReal3MatrixToBitmap(MatrixType *M,FILE *stream,void(*Real3ToRGBColor)(double *,unsigned char *));
uint32_t PrintRealMatrixToBitmap(MatrixType M,void(*RealToRGBColor)(double,unsigned char *));
int FilePrintRealMatrixToBitmap(MatrixType M,const char *const flname,void(*RealToRGBColor)(double,unsigned char *));
int StreamPrintRealMatrixToBitmap(MatrixType M,FILE *stream,void(*RealToRGBColor)(double,unsigned char *));
void Real3ToRGBColorSimpleInverse(double V[3],unsigned char *RGB);
void RealToRGBColorBWInverse(double V,unsigned char *RGB);
void RealToRGBColorBW(double V,unsigned char *RGB);
