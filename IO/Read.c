#include "Read.h"
#include "StreamRead.c"
#include "FileRead.c"

/*  Functions come in three variants
*   StreamRead***
*   FileRead***
*   Read***
*
*   The functions StreamRead expect as an input besides the structure to read also a stream name to read from
*   the stream should be opened to be able to read from
*   The functions FileRead open and close a file from which is read
*   The functions Read read from the standard input stream
*
*   The functions FileRead*** and Read*** make calls to the function StreamRead
*/

uint32_t ReadWaveFunctionPlain(WaveFunctionType * O)
{
  return StreamReadWaveFunctionPlain(O, stdin);
}

uint32_t ReadWaveFunctionPlainBinary(WaveFunctionType * O)
{
  return StreamReadWaveFunctionPlainBinary(O, stdin);
}

uint32_t ReadMatrixType(MatrixType * O)
{
  return StreamReadMatrixType(O, stdin);
}

uint32_t ReadPlainAndersonMatrixType(AndersonMatrixType * O)
{
  return StreamReadPlainAndersonMatrixType(O, stdin);
}

uint32_t ReadPlainTriDiagonalMatrix(TriDiagonalMatrixType * M)
{
  return StreamReadPlainTriDiagonalMatrixType(M, stdout);
}

uint32_t ReadPlainAndersonMatrix(AndersonMatrixType * M)
{
  return StreamReadPlainAndersonMatrixType(M, stdout);
}

uint32_t ReadListOfPoles(ListOfPoles * M)
{
  return StreamReadPlainListOfPoles(M, stdout);
}

uint32_t ReadListOfPolesUndetermined(ListOfPoles * M)
{
  return StreamReadListOfPolesUndetermined(M, stdout);
}

uint32_t ReadAtom(AtomType * O)
{
  return StreamReadAtom(O, stdin);
}

uint32_t ReadCrystalStructure(CrystalStructureType * O)
{
  return StreamReadCrystalStructure(O, stdin);
}

uint32_t ReadCrystalStructureSuperCell(CrystalStructureSuperCellType * O)
{
  return StreamReadCrystalStructureSuperCell(O, stdin);
}

uint32_t ReadHeisenbergModelSettings(HeisenbergModelSettings * O)
{
  return StreamReadHeisenbergModelSettings(O, stdin);
}

uint32_t ReadDimerizedSpinChainSettings(DimerizedSpinChainSettings * O)
{
  return StreamReadDimerizedSpinChainSettings(O, stdin);
}

uint32_t ReadtJModelSettings(tJModelSettings * O)
{
  return StreamReadtJModelSettings(O, stdin);
}

uint32_t ReadSpinLaderSettings(SpinLaderSettings * O)
{
  return StreamReadSpinLaderSettings(O, stdin);
}

uint32_t ReadHubbardModelSettings(HubbardModelSettings * O)
{
  return StreamReadHubbardModelSettings(O, stdin);
}

uint32_t ReadDMFTXASSettings(DMFTXASSettings * O)
{
  return StreamReadDMFTXASSettings(O, stdin);
}

uint32_t ReadDMFTSettings(DMFTSettings * O)
{
  return StreamReadDMFTSettings(O, stdin);
}
