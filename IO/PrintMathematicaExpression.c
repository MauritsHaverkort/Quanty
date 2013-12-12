/*  Prints a matrix
 */
uint32_t StreamPrintSparseMatrixMathematica(SparseMatrixType M, FILE * stream)
{
  uint32_t i;
  fprintf(stream, "SparseArray[{");
  if (M.QComplex)
  {
    for (i = 0; i < M.N; i++)
    {
      fprintf(stream, "{%i, %i} -> %24.15lf + I %24.15lf,", M.i[i] + 1, M.j[i] + 1, M.VRe[i], M.VIm[i]);
    }
  }
  else
  {
    for (i = 0; i < M.N; i++)
    {
      fprintf(stream, "{%i, %i} -> %24.15lf,", M.i[i] + 1, M.j[i] + 1, M.VRe[i]);
    }
  }
  fprintf(stream, "{_,_}->0},{%i,%i}]\n", M.Ni, M.Nj);
  return 0;
}

uint32_t FilePrintSparseMatrixMathematica(SparseMatrixType M, const char *const flname, const char *const mode)
{
  FILE  *stream;
  uint32_t retval;
  if ((strcmp(mode, "w") && strcmp(mode, "a")) != 0)
  {
    printf
      ("ERROR in FilePrintSparseMatrix.\nmode must be either \"w\" for overwrite or \"a\" for append\n %s as input is not excepted\n",
       mode);
    return 1;
  }
  if ((stream = fopen(flname, mode)) == NULL)
  {
    printf("ERROR\nCan't open %s for writing with mode %s.\n", flname, mode);
    return 1;
  }
  retval = StreamPrintSparseMatrixMathematica(M, stream);
  fclose(stream);
  return retval;
}

uint32_t PrintSparseMatrixMathematica(SparseMatrixType M)
{
  return StreamPrintSparseMatrixMathematica(M, stdout);
}
