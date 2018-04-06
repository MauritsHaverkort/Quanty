#include "TemperatureAvverage.h"

typedef struct TemperatureSettings
{
  uint32_t NFolders;            /* Number of folders needed to average over */
  char   FolderName[100][256];  /* Name of the folder */
  char   FileName[100][256];    /* General root of the file name within the corresponding folder */
  uint32_t LowestFileNumber[100]; /* Lowest file number */
  uint32_t HighestFileNumber[100];  /* Highest file number */
  char   EnergyFileName[100][256];  /* File name where the energy can be found */
  double Degeneracy[100];       /* Multiplication factor for these spectra */
  uint32_t NTemperatures;       /* number of temperatures to average over */
  double Temperatures[256];     /* Temperatures used */
  uint32_t Nq[3];               /* Number of twisted boundary phases */
}
TemperatureSettings;


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
