#include "StreamPrint.c"
#include "FilePrint.c"
#include "PrintMathematicaExpression.c"
#include "StringPrint.c"


void VerbosityPrintf(const uint32_t verbosity, char *fmt, ...)
{
  va_list args;
  va_start ( args, fmt );
  if(verbosity&GlobalVerbosity) vprintf(fmt,args);
  va_end(args);
}


/*  Functions come in three variants
*   StreamPrint***
*   FilePrint***
*   Print***
*
*   The functions StreamPrint expect as an input besides the structure to print also a stream name to print to
*   the stream should be opened to be able to print to
*   The functions FilePrint open and close a file to which is printed
*   The functions Print print to the standard output stream
*
*   The functions FilePrint*** and Print*** make calls to the function StreamPrint
*/

uint32_t PrintSpectra(SpectraType * const Spec, const uint32_t NSpec)
{
  return StreamPrintSpectra(Spec, NSpec, stdout);
}

uint32_t PrintRIXSSpectrum(const uint32_t numwRes, SpectraType (*Spec)[numwRes], double *wRes, int printReal, const uint32_t iSpec)
{
  return StreamPrintRIXSSpectrum(numwRes, Spec, wRes, printReal, iSpec, stdout);
}

uint32_t PrintRIXSSpectra(const uint32_t numwRes, SpectraType (*Spec)[numwRes], double *wRes, const uint32_t NSpec)
{
  return StreamPrintRIXSSpectra(numwRes, Spec, wRes, NSpec, stdout);
}

uint32_t PrintWaveFunctionPlain(WaveFunctionType psi)
{
  return StreamPrintWaveFunctionPlain(psi, stdout);
}

uint32_t PrintWaveFunctionPlainBinary(WaveFunctionType psi)
{
  return StreamPrintWaveFunctionPlainBinary(psi, stdout);
}

uint32_t PrintWaveFunctionReadable(WaveFunctionType psi)
{
  return StreamPrintWaveFunctionReadable(psi, stdout);
}

uint32_t PrintWaveFunctionReadableTruncated(WaveFunctionType psi, double epsilon)
{
  return StreamPrintWaveFunctionReadableTruncated(psi, epsilon, stdout);
}

uint32_t PrintWaveFunctionFullReadable(WaveFunctionType psi)
{
  return StreamPrintWaveFunctionFullReadable(psi, stdout);
}

uint32_t PrintWaveFunctionDensity(WaveFunctionType psi)
{
  return StreamPrintWaveFunctionDensity(psi, stdout);
}

uint32_t PrintWaveFunctionPartialCorrelationTensor(WaveFunctionType psi, uint32_t const NFermionIncluded, uint16_t * const FermionIncluded)
{
  return StreamPrintWaveFunctionPartialCorrelationTensor(psi, stdout, NFermionIncluded, FermionIncluded);
}

uint32_t PrintWaveFunctionPrefactorHistogram(WaveFunctionType psi, const double exp, const uint32_t N)
{
  return StreamPrintWaveFunctionPrefactorHistogram(psi, exp, N, stdout);
}

uint32_t PrintOperatorReadable(OperatorType O)
{
  return StreamPrintOperatorReadable(O, stdout);
}

uint32_t PrintOperatorFullReadable(OperatorType O)
{
  return StreamPrintOperatorFullReadable(O, stdout);
}

uint32_t PrintCrystalStructure(CrystalStructureType O)
{
  return StreamPrintCrystalStructure(O, stdout);
}

uint32_t PrintCrystalStructureSuperCell(CrystalStructureSuperCellType O)
{
  return StreamPrintCrystalStructureSuperCell(O, stdout);
}

uint32_t PrintMatrix(MatrixType M)
{
  return StreamPrintMatrix(M, stdout);
}

uint32_t PrintPlainMatrix(MatrixType M)
{
  return StreamPrintPlainMatrix(M, stdout);
}

uint32_t PrintCompactMatrixMinimal(CompactMatrixTypeOld M)
{
  return StreamPrintCompactMatrixMinimal(M, stdout);
}

uint32_t PrintPlainTriDiagonalMatrix(TriDiagonalMatrixType M)
{
  return StreamPrintPlainTriDiagonalMatrix(M, stdout);
}

uint32_t PrintPlainAndersonMatrix(AndersonMatrixType M)
{
  return StreamPrintPlainAndersonMatrix(M, stdout);
}

uint32_t PrintPlainListOfPoles(ListOfPoles M)
{
  return StreamPrintPlainListOfPoles(M, stdout);
}

uint32_t PrintBlockListOfPoles(BlockListOfPoles M)
{
  return StreamPrintBlockListOfPoles(M, stdout);
}

uint32_t PrintTightBindingDef(TightBindingDef O)
{
  return StreamPrintTightBindingDef(O, stdout);
}

uint32_t PrintDimerizedSpinChainSettings(DimerizedSpinChainSettings O)
{
  return StreamPrintDimerizedSpinChainSettings(O, stdout);
}

uint32_t PrintHeisenbergModelSettings(HeisenbergModelSettings O)
{
  return StreamPrintHeisenbergModelSettings(O, stdout);
}

uint32_t PrinttJModelSettings(tJModelSettings O)
{
  return StreamPrinttJModelSettings(O, stdout);
}

uint32_t PrintSpinLaderSettings(SpinLaderSettings O)
{
  return StreamPrintSpinLaderSettings(O, stdout);
}

uint32_t PrintHubbardModelSettings(HubbardModelSettings O)
{
  return StreamPrintHubbardModelSettings(O, stdout);
}

uint32_t PrintDMFTXASSettings(DMFTXASSettings O)
{
  return StreamPrintDMFTXASSettings(O, stdout);
}

uint32_t PrintDMFTSettings(DMFTSettings O)
{
  return StreamPrintDMFTSettings(O, stdout);
}

uint32_t PrintSparseMatrix(SparseMatrixType M)
{
  return StreamPrintSparseMatrix(M, stdout);
}

uint32_t PrintTriDiagonalMatrix(TriDiagonalMatrixType M)
{
  return StreamPrintTriDiagonalMatrix(M, stdout);
}

uint32_t PrintBlokTriDiagonalMatrix(BlokTriDiagonalMatrixType M)
{
  return StreamPrintBlokTriDiagonalMatrix(M, stdout);
}

uint32_t PrintBlockTriDiagonalMatrix(BlockTriDiagonalMatrixType M)
{
  return StreamPrintBlockTriDiagonalMatrix(M, stdout);
}

uint32_t PrintDoubleTriDiagonalMatrix(DoubleTriDiagonalMatrixType M)
{
  return StreamPrintDoubleTriDiagonalMatrix(M, stdout);
}

uint32_t PrintBADoubleTriDiagonalMatrix(BADoubleTriDiagonalMatrixType M)
{
  return StreamPrintBADoubleTriDiagonalMatrix(M, stdout);
}

uint32_t PrintSecondDoubleTriDiagonalMatrix(SecondDoubleTriDiagonalMatrixType M)
{
  return StreamPrintSecondDoubleTriDiagonalMatrix(M, stdout);
}

uint32_t PrintBlokDoubleTriDiagonalMatrix(BlokDoubleTriDiagonalMatrixType M)
{
  return StreamPrintBlokDoubleTriDiagonalMatrix(M, stdout);
}

uint32_t PrintAndersonMatrix(AndersonMatrixType M)
{
  return StreamPrintAndersonMatrix(M, stdout);
}

uint32_t PrintAndersonMatrixEssentials(AndersonMatrixType M)
{
  return StreamPrintAndersonMatrixEssentials(M, stdout);
}

uint32_t PrintBlokAndersonMatrix(BlokAndersonMatrixType M)
{
  return StreamPrintBlokAndersonMatrix(M, stdout);
}

uint32_t PrintBlockAndersonMatrix(BlockAndersonMatrixType M)
{
  return StreamPrintBlockAndersonMatrix(M, stdout);
}

uint32_t PrintListOfPoles(ListOfPoles M)
{
  return StreamPrintListOfPoles(M, stdout);
}

uint32_t PrintListOfPolesEssentials(ListOfPoles M)
{
  return StreamPrintListOfPolesEssentials(M, stdout);
}

uint32_t PrintHarmonicsExpansion(HarmonicsExpansion M)
{
  return StreamPrintHarmonicsExpansion(M, stdout);
}

uint32_t PrintBandstructurePath(BandstructurePathType M)
{
  return StreamPrintBandstructurePath(M, stdout);
}

uint32_t PrintBlockNaturalImpurityOrbitalMatrix(BlockNaturalImpurityOrbitalMatrixType M)
{
  return StreamPrintBlockNaturalImpurityOrbitalMatrix(M, stdout);
}

