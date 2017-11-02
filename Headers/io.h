

/**********************************************************/
/*****************       IO/Print.c       *****************/
/**********************************************************/

void     VerbosityPrintf(const uint32_t verbosity, char *fmt, ...);
uint32_t PrintSpectra(SpectraType * const Spec, const uint32_t NSpec);
uint32_t PrintWaveFunctionReadable(WaveFunctionType psi);
uint32_t PrintWaveFunctionFullReadable(WaveFunctionType psi);
uint32_t PrintWaveFunctionDensity(WaveFunctionType psi);
uint32_t PrintWaveFunctionPrefactorHistogram(WaveFunctionType psi, const double exp, const uint32_t N);
uint32_t PrintOperatorReadable(OperatorType O);
uint32_t PrintOperatorFullReadable(OperatorType O);
uint32_t PrintMatrix(MatrixType M);
uint32_t PrintCompactMatrixMinimal(CompactMatrixType M);
uint32_t PrintSparseMatrix(SparseMatrixType M);
uint32_t PrintTriDiagonalMatrix(TriDiagonalMatrixType M);
uint32_t PrintBlokTriDiagonalMatrix(BlokTriDiagonalMatrixType M);
uint32_t PrintBlockTriDiagonalMatrix(BlockTriDiagonalMatrixType M);
uint32_t PrintDoubleTriDiagonalMatrix(DoubleTriDiagonalMatrixType M);
uint32_t PrintBlokDoubleTriDiagonalMatrix(BlokDoubleTriDiagonalMatrixType M);
uint32_t PrintAndersonMatrix(AndersonMatrixType M);
uint32_t PrintBlokAndersonMatrix(BlokAndersonMatrixType M);
uint32_t PrintListOfPoles(ListOfPoles M);
uint32_t PrintOperatorFullReadable(OperatorType O);
uint32_t PrintListOfPolesEssentials(ListOfPoles M);
uint32_t PrintAndersonMatrixEssentials(AndersonMatrixType M);
uint32_t PrintBandstructurePath(BandstructurePathType M);
uint32_t PrintBlockNaturalImpurityOrbitalMatrix(BlockNaturalImpurityOrbitalMatrixType M);
uint32_t PrintBlockListOfPoles(BlockListOfPoles M);

uint32_t FilePrintSpectra(SpectraType * const Spec, const uint32_t NSpec, const char *const flname, const char *const mode);
uint32_t FilePrintWaveFunctionReadable(WaveFunctionType psi, const char *const flname, const char *const mode);
uint32_t FilePrintWaveFunctionFullReadable(WaveFunctionType psi, const char *const flname, const char *const mode);
uint32_t FilePrintWaveFunctionDensity(WaveFunctionType psi, const char *const flname, const char *const mode);
uint32_t FilePrintWaveFunctionPrefactorHistogram(WaveFunctionType psi, const double exp, const uint32_t N, const char *const flname, const char *const mode);
uint32_t FilePrintOperatorReadable(OperatorType O, const char *const flname, const char *const mode);
uint32_t FilePrintOperatorFullReadable(OperatorType O, const char *const flname, const char *const mode);
uint32_t FilePrintMatrix(MatrixType M, const char *const flname, const char *const mode);
uint32_t FilePrintCompactMatrixMinimal(CompactMatrixType M, const char *const flname, const char *const mode);
uint32_t FilePrintSparseMatrix(SparseMatrixType M, const char *const flname, const char *const mode);
uint32_t FilePrintPlainTriDiagonalMatrix(TriDiagonalMatrixType M, const char *const flname, const char *const mode);
//uint32_t FilePrintBlokTriDiagonalMatrix();
//uint32_t FilePrintBlockTriDiagonalMatrix();
//uint32_t FilePrintDoubleTriDiagonalMatrix();
//uint32_t FilePrintBlokDoubleTriDiagonalMatrix();
uint32_t FilePrintAndersonMatrix(AndersonMatrixType M, const char *const flname, const char *const mode);
//uint32_t FilePrintBlokAndersonMatrix();
uint32_t FilePrintListOfPoles(ListOfPoles M, const char *const flname, const char *const mode);
uint32_t FilePrintOperatorFullReadable(OperatorType O, const char *const flname, const char *const mode);
uint32_t FilePrintListOfPolesEssentials(ListOfPoles M, const char *const flname, const char *const mode);
uint32_t FilePrintAndersonMatrixEssentials(AndersonMatrixType M, const char *const flname, const char *const mode);
uint32_t FilePrintBandstructurePath(BandstructurePathType M, const char *const flname, const char *const mode);
uint32_t FilePrintBlockNaturalImpurityOrbitalMatrix(BlockNaturalImpurityOrbitalMatrixType M, const char *const flname, const char *const mode);
uint32_t FilePrintBlockListOfPoles(BlockListOfPoles M, const char *const flname, const char *const mode);

uint32_t StreamPrintSpectra(SpectraType * const Spec, const uint32_t NSpec, FILE * stream);
uint32_t StreamPrintWaveFunctionReadable(WaveFunctionType psi, FILE * stream);
uint32_t StreamPrintWaveFunctionFullReadable(WaveFunctionType psi, FILE * stream);
uint32_t StreamPrintWaveFunctionDensity(WaveFunctionType const psi, FILE * stream);
uint32_t StreamPrintWaveFunctionPrefactorHistogram(WaveFunctionType const psi, const double exp, const uint32_t N, FILE * stream);
uint32_t StreamPrintOperatorReadable(OperatorType O, FILE * stream);
uint32_t StreamPrintOperatorFullReadable(OperatorType O, FILE * stream);
uint32_t StreamPrintMatrix(MatrixType M, FILE * stream);
uint32_t StreamPrintCompactMatrixMinimal(CompactMatrixType const M, FILE * stream);
uint32_t StreamPrintSparseMatrix(SparseMatrixType M, FILE * stream);
uint32_t StreamPrintTriDiagonalMatrix(TriDiagonalMatrixType M, FILE * stream);
uint32_t StreamPrintBlokTriDiagonalMatrix(BlokTriDiagonalMatrixType M, FILE * stream);
uint32_t StreamPrintBlockTriDiagonalMatrix(BlockTriDiagonalMatrixType M, FILE * stream);
uint32_t StreamPrintDoubleTriDiagonalMatrix(DoubleTriDiagonalMatrixType M, FILE * stream);
uint32_t StreamPrintBlokDoubleTriDiagonalMatrix(BlokDoubleTriDiagonalMatrixType M, FILE * stream);
uint32_t StreamPrintAndersonMatrix(AndersonMatrixType M, FILE * stream);
uint32_t StreamPrintBlokAndersonMatrix(BlokAndersonMatrixType M, FILE * stream);
uint32_t StreamPrintListOfPoles(ListOfPoles M, FILE * stream);
uint32_t StreamPrintOperatorFullReadable(OperatorType O, FILE * stream);
uint32_t StreamPrintListOfPolesEssentials(ListOfPoles M, FILE * stream);
uint32_t StreamPrintAndersonMatrixEssentials(AndersonMatrixType M, FILE * stream);
uint32_t StreamPrintBandstructurePath(BandstructurePathType M, FILE * stream);
uint32_t StreamPrintBlockNaturalImpurityOrbitalMatrix(BlockNaturalImpurityOrbitalMatrixType M, FILE * stream);
uint32_t StreamPrintBlockListOfPoles(BlockListOfPoles const M, FILE * stream);

uint32_t PrintSparseMatrixMathematica(SparseMatrixType M);
uint32_t StreamPrintSparseMatrixMathematica(SparseMatrixType M, FILE * stream);
uint32_t FilePrintSparseMatrixMathematica(SparseMatrixType M, const char *const flname, const char *const mode);
