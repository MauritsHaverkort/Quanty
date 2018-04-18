


void ReadTemperatureSettings(TemperatureSettings * const A, const char *const flname)
{
  FILE  *stream;
  char   comments[6000];
  uint32_t i;
  stream = fopen(flname, "r");
  if (stream == NULL)
  {
    printf("Can't open file %s\n", flname);
    printf("That's bad program terminates\n");
    exit(1);
  }
  fscanf(stream, "%i", &(A->NFolders));
  fgets(comments, 6000, stream);
  fscanf(stream, "%i %i %i", &(A->Nq[0]), &(A->Nq[1]), &(A->Nq[2]));
  fgets(comments, 6000, stream);
  for (i = 0; i < A->NFolders; i++)
  {
    fscanf(stream, "%s", A->FolderName[i]);
    fscanf(stream, "%s", A->FileName[i]);
    fscanf(stream, "%i", &(A->LowestFileNumber[i]));
    fscanf(stream, "%i", &(A->HighestFileNumber[i]));
    fscanf(stream, "%s", A->EnergyFileName[i]);
    fscanf(stream, "%lf", &(A->Degeneracy[i]));
    fgets(comments, 6000, stream);
  }
  fscanf(stream, "%i", &(A->NTemperatures));
  fgets(comments, 6000, stream);
  for (i = 0; i < A->NTemperatures; i++)
  {
    fscanf(stream, "%lf", &(A->Temperatures[i]));
  }
  fclose(stream);

  printf("Read Temperature settings from file %s:\n", flname);
  printf("NFolders %i\n", (A->NFolders));
  printf("Number of q points %i %i %i\n", A->Nq[0], A->Nq[1], A->Nq[2]);
  printf("Folder Names\n");
  for (i = 0; i < A->NFolders; i++)
  {
    printf("%3i %s\n", i, A->FolderName[i]);
  }
  printf("File Names with lowest and highest index number \n");
  for (i = 0; i < A->NFolders; i++)
  {
    printf("%3i %s %4i %4i\n", i, A->FileName[i], A->LowestFileNumber[i], A->HighestFileNumber[i]);
  }
  printf("Energy file name\n");
  for (i = 0; i < A->NFolders; i++)
  {
    printf("%3i %s\n", i, A->EnergyFileName[i]);
  }
  printf("Degeneracy\n");
  for (i = 0; i < A->NFolders; i++)
  {
    printf("%3i %22.15E\n", i, A->Degeneracy[i]);
  }
  printf("Temperature\n");
  for (i = 0; i < A->NTemperatures; i++)
  {
    printf("%3i %22.15E\n", i, A->Temperatures[i]);
  }
}


void TemperatureAvverage(const char *const flname)
{
  TemperatureSettings A;
  ReadTemperatureSettings(&A, flname);
}
