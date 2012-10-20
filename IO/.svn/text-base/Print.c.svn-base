#include "StreamPrint.c"
#include "FilePrint.c"

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

uint32_t PrintSpectra(SpectraType *const Spec, const uint32_t NSpec)
{ return StreamPrintSpectra(Spec,NSpec,stdout);}

uint32_t PrintWaveFunctionPlain(WaveFunctionType psi)
{ return StreamPrintWaveFunctionPlain(psi,stdout);}

uint32_t PrintWaveFunctionPlainBinary(WaveFunctionType psi)
{ return StreamPrintWaveFunctionPlainBinary(psi,stdout);}

uint32_t PrintWaveFunctionReadable(WaveFunctionType psi)
{ return StreamPrintWaveFunctionReadable(psi,stdout);}

uint32_t PrintWaveFunctionReadableTruncated(WaveFunctionType psi, double epsilon)
{ return StreamPrintWaveFunctionReadableTruncated(psi, epsilon, stdout);}

uint32_t PrintWaveFunctionFullReadable(WaveFunctionType psi)
{ return StreamPrintWaveFunctionFullReadable(psi,stdout);}

uint32_t PrintOperatorReadable(OperatorType O)
{ return StreamPrintOperatorReadable(O,stdout);}

uint32_t PrintOperatorFullReadable(OperatorType O)
{ return StreamPrintOperatorFullReadable(O,stdout);}

uint32_t PrintCrystalStructure(CrystalStructureType O)
{ return StreamPrintCrystalStructure(O,stdout);}

uint32_t PrintCrystalStructureSuperCell(CrystalStructureSuperCellType O)
{ return StreamPrintCrystalStructureSuperCell(O,stdout);}

uint32_t PrintMatrix(MatrixType M)
{ return StreamPrintMatrix(M,stdout);}

uint32_t PrintPlainMatrix(MatrixType M)
{ return StreamPrintPlainMatrix(M,stdout);}

uint32_t PrintPlainTriDiagonalMatrix(TriDiagonalMatrixType M)
{ return StreamPrintPlainTriDiagonalMatrix(M,stdout);}

uint32_t PrintTightBindingDef(TightBindingDef O)
{ return StreamPrintTightBindingDef(O,stdout);}

uint32_t PrintDimerizedSpinChainSettings(DimerizedSpinChainSettings O)
{ return StreamPrintDimerizedSpinChainSettings(O,stdout);}

uint32_t PrintHeisenbergModelSettings(HeisenbergModelSettings O)
{ return StreamPrintHeisenbergModelSettings(O,stdout);}

uint32_t PrinttJModelSettings(tJModelSettings O)
{ return StreamPrinttJModelSettings(O,stdout);}

uint32_t PrintSpinLaderSettings(SpinLaderSettings O)
{ return StreamPrintSpinLaderSettings(O,stdout);}

uint32_t PrintHubbardModelSettings(HubbardModelSettings O)
{ return StreamPrintHubbardModelSettings(O,stdout);}

uint32_t PrintDMFTXASSettings(DMFTXASSettings O)
{ return StreamPrintDMFTXASSettings(O,stdout);}

uint32_t PrintDMFTSettings(DMFTSettings O)
{ return StreamPrintDMFTSettings(O,stdout);}







/* old */

uint32_t PrintSparseMatrix(SparseMatrixType M)
{
  uint32_t i;
  M.Name[MAXStrLength-1]='\0';
  printf("%s =\n",M.Name);
  if(M.QComplex)
  {
    printf("   [i]    [j]          Re(M[i][j])          Im(M[i][j])\n");
    for(i=0;i<M.N;i++)
    {
      printf("%6i %6i %22.15E %22.15E\n",M.i[i],M.j[i],M.VRe[i],M.VRe[i]);
    }
  }
  else
  {
    printf("   [i]    [j]          Re(M[i][j])\n");
    for(i=0;i<M.N;i++)
    {
      printf("%6i %6i %22.15E\n",M.i[i],M.j[i],M.VRe[i]);
    }
  }
  printf("\n");
  return 0;
}

uint32_t PrintTriDiagonalMatrix(TriDiagonalMatrixType M)
{
  uint32_t i;
  M.Name[MAXStrLength-1]='\0';
  printf("%s =\n",M.Name);
  printf("TriDiagonalMatrix element a[i-1]=M[i][i] element b[i-1]=M[i,i-1]=M[i-1,i] others are zero.\n");
  printf("Zero elements are prefactors for continued fraction representation\n");
  printf("i    a[i]                   b[i]\n");
  for(i=0;i<M.N;i++)
  {
    printf("%4i %22.15E %22.15E\n",i,M.a[i],M.b[i]);
  }
  printf("%4i %22.15E                       \n",M.N,M.a[M.N]);
  return 0;
}

uint32_t PrintBlokTriDiagonalMatrix(BlokTriDiagonalMatrixType M)
{
  uint32_t i;
  M.Name[MAXStrLength-1]='\0';
  printf("%s =\n",M.Name);
  printf("BlokTriDiagonalMatrix element a[i-1]=M[i][i] element b[i-1]=M[i,i+1]=M[i+1,i] others are zero.\n");
  printf("Zero elements are prefactors for continued fraction representation\n");
  for(i=0;i<=M.N;i++)
  {
    printf("a[%i]\n",i);
    PrintMatrix(M.a[i]);
  }
  for(i=0;i<M.N;i++)
  {
    printf("b[%i]\n",i);
    PrintMatrix(M.b[i]);
  }
  return 0;
}

uint32_t PrintDoubleTriDiagonalMatrix(DoubleTriDiagonalMatrixType M)
{
  M.Name[MAXStrLength-1]='\0';
  printf("%s =\n",M.Name);
  printf("DoubleTriDiagonalMatrix Ef = %22.15E Epsilon = %22.15E\n",M.Ef,M.Epsilon);
  printf("                        a0 = %22.15E      b0 = %22.15E\n",M.a0,M.b0);
  printf("Valence part (E<Ef):\n");
  PrintTriDiagonalMatrix(M.val);
  printf("Conduction part (E>Ef):\n");
  PrintTriDiagonalMatrix(M.con);
  return 0;
}

uint32_t PrintBlokDoubleTriDiagonalMatrix(BlokDoubleTriDiagonalMatrixType M)
{
  M.Name[MAXStrLength-1]='\0';
  printf("%s =\n",M.Name);
  printf("BlokDoubleTriDiagonalMatrix Ef = %22.15E Epsilon = %22.15E\n",M.Ef,M.Epsilon);
  printf("a[0]\n");
  PrintMatrix(M.a0);
  printf("b[0]\n");
  PrintMatrix(M.b0);
  printf("Valence part (E<Ef):\n");
  PrintBlokTriDiagonalMatrix(M.val);
  printf("Conduction part (E>Ef):\n");
  PrintBlokTriDiagonalMatrix(M.con);
  return 0;
}

uint32_t PrintAndersonMatrix(AndersonMatrixType M)
{
  uint32_t i;
  M.Name[MAXStrLength-1]='\0';
  printf("%s =\n",M.Name);
  printf("AndersonMatrix element a[i-1]=M[i][i] element b[i-1]=M[i,i-1]=M[i-1,i] others are zero.\n");
  printf("Zero elements are prefactors for continued fraction representation\n");
  printf("i    a[i]                   b[i]\n");
  for(i=0;i<M.N;i++)
  {
    printf("%4i %22.15E %22.15E\n",i,M.a[i],M.b[i]);
  }
  printf("%4i %22.15E                       \n",M.N,M.a[M.N]);
  return 0;
}

uint32_t PrintBlokAndersonMatrix(BlokAndersonMatrixType M)
{
  uint32_t i;
  M.Name[MAXStrLength-1]='\0';
  printf("%s =\n",M.Name);
  printf("BlokAndersonMatrix element a[i]=M[i][i] element b[i]=M[0,i+1]=M[i+1,0] others are zero.\n");
  printf("Zero elements are prefactors for continued fraction representation\n");
  for(i=0;i<=M.N;i++)
  {
    printf("a[%i]\n",i);
    PrintMatrix(M.a[i]);
  }
  for(i=0;i<M.N;i++)
  {
    printf("b[%i]\n",i);
    PrintMatrix(M.b[i]);
  }
  return 0;
}


uint32_t PrintListOfPoles(ListOfPoles M)
{
  uint32_t i;
  M.Name[MAXStrLength-1]='\0';
  printf("%s =\n",M.Name);
  printf("List of Poles.\n");
  printf("i    a[i]                   b[i]\n");
  for(i=0;i<M.N;i++)
  {
    printf("%4i %22.15E %22.15E\n",i,M.a[i],M.b[i]);
  }
  return 0;
}














