

/**********************************************************/
/*****************       IO/Print.c       *****************/
/**********************************************************/

uint32_t PrintSpectra();
uint32_t PrintWaveFunctionReadable();
uint32_t PrintWaveFunctionFullReadable();
uint32_t PrintWaveFunctionDensity();
uint32_t PrintWaveFunctionPrefactorHistogram();
uint32_t PrintOperatorReadable();
uint32_t PrintOperatorFullReadable();
uint32_t PrintMatrix();
uint32_t PrintCompactMatrixMinimal();
uint32_t PrintSparseMatrix();
uint32_t PrintTriDiagonalMatrix();
uint32_t PrintBlokTriDiagonalMatrix();
uint32_t PrintBlockTriDiagonalMatrix();
uint32_t PrintDoubleTriDiagonalMatrix();
uint32_t PrintBlokDoubleTriDiagonalMatrix();
uint32_t PrintAndersonMatrix();
uint32_t PrintBlokAndersonMatrix();
uint32_t PrintListOfPoles();
uint32_t PrintOperatorFullReadable(OperatorType O);
uint32_t PrintListOfPolesEssentials();
uint32_t PrintAndersonMatrixEssentials();
uint32_t PrintBandstructurePath();

uint32_t FilePrintSpectra();
uint32_t FilePrintWaveFunctionReadable();
uint32_t FilePrintWaveFunctionFullReadable();
uint32_t FilePrintWaveFunctionDensity();
uint32_t FilePrintWaveFunctionPrefactorHistogram();
uint32_t FilePrintOperatorReadable();
uint32_t FilePrintOperatorFullReadable();
uint32_t FilePrintMatrix();
uint32_t FilePrintCompactMatrixMinimal();
uint32_t FilePrintSparseMatrix();
uint32_t FilePrintPlainTriDiagonalMatrix();
uint32_t FilePrintBlokTriDiagonalMatrix();
uint32_t FilePrintBlockTriDiagonalMatrix();
uint32_t FilePrintDoubleTriDiagonalMatrix();
uint32_t FilePrintBlokDoubleTriDiagonalMatrix();
uint32_t FilePrintAndersonMatrix();
uint32_t FilePrintBlokAndersonMatrix();
uint32_t FilePrintListOfPoles();
uint32_t FilePrintOperatorFullReadable(OperatorType O, const char *const flname, const char *const mode);
uint32_t FilePrintListOfPolesEssentials();
uint32_t FilePrintAndersonMatrixEssentials();
uint32_t FilePrintBandstructurePath();

uint32_t StreamPrintSpectra();
uint32_t StreamPrintWaveFunctionReadable();
uint32_t StreamPrintWaveFunctionFullReadable();
uint32_t StreamPrintWaveFunctionDensity();
uint32_t StreamPrintWaveFunctionPrefactorHistogram();
uint32_t StreamPrintOperatorReadable();
uint32_t StreamPrintOperatorFullReadable();
uint32_t StreamPrintMatrix();
uint32_t StreamPrintCompactMatrixMinimal();
uint32_t StreamPrintSparseMatrix();
uint32_t StreamPrintTriDiagonalMatrix();
uint32_t StreamPrintBlokTriDiagonalMatrix();
uint32_t StreamPrintBlockTriDiagonalMatrix();
uint32_t StreamPrintDoubleTriDiagonalMatrix();
uint32_t StreamPrintBlokDoubleTriDiagonalMatrix();
uint32_t StreamPrintAndersonMatrix();
uint32_t StreamPrintBlokAndersonMatrix();
uint32_t StreamListOfPoles();
uint32_t StreamPrintOperatorFullReadable(OperatorType O, FILE * stream);
uint32_t StreamPrintListOfPolesEssentials();
uint32_t StreamPrintAndersonMatrixEssentials();
uint32_t StreamPrintBandstructurePath();


uint32_t PrintSparseMatrixMathematica();
uint32_t StreamPrintSparseMatrixMathematica();
uint32_t FilePrintSparseMatrixMathematica();
