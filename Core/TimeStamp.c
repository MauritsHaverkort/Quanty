#ifndef __TIMESTAMP__
# define __TIMESTAMP__ "???"
#endif
#ifndef DBL_MIN
# define DBL_MIN 2.3E-308
#endif
#ifndef DBL_EPSILON
# define DBL_EPSILON 2.3E-16
#endif
#ifndef _POSIX_HOST_NAME_MAX
# ifndef HOST_NAME_MAX
#  define HOST_NAME_MAX 1023
# endif
# define _POSIX_HOST_NAME_MAX HOST_NAME_MAX
#endif
#define TimeStampAll() printf("=============================================================\n");\
                       printf("====    written by Maurits W. Haverkort                  ====\n");\
                       printf("====    with contributions from:                         ====\n");\
                       printf("====    Yi Lu, Robert Green and Sebastian Macke          ====\n");\
                       printf("====    (C) 1995-2015   All rights reserved              ====\n");\
                       printf("====    www.quanty.org                                   ====\n");\
                       printf("====    Beta version, be critical and report errors!!!   ====\n");\
                       printf("=============================================================\n");\
                       printf("====    Version compiled at: %s at %s     ====\n",__DATE__,__TIME__);\
                       printf("=============================================================\n");\
                       printf("====    When used in scientific publications please cite ====\n");\
                       printf("====    Phys. Rev. B 85, 165113 (2012)                   ====\n");\
                       printf("====    Phys. Rev. B 90, 085102 (2014)                   ====\n");\
                       printf("====    Euro Phys. Lett. 108, 57004 (2014)               ====\n");\
                       printf("=============================================================\n");\
                       time_t MWHrawtime;\
                       struct tm * MWHtimeinfo;\
                       time ( &MWHrawtime );\
                       MWHtimeinfo = localtime ( &MWHrawtime );\
                       printf("Program executed on: %s", asctime (MWHtimeinfo) );\
                       char MWHHostName[_POSIX_HOST_NAME_MAX+1];\
                       MWHgethostname(MWHHostName, _POSIX_HOST_NAME_MAX+1);\
                       printf("Running on host    : %s\n",MWHHostName);\
                       printf("number of available processors              : %i\n",omp_get_num_procs());\
                       printf("maximum number of threads in parallel region: %i\n",omp_get_max_threads());\
                       printf("Smallest positive float  : %E \n",DBL_MIN);\
                       printf("Smallest deviation from 1: %E \n",DBL_EPSILON);\
                       printf("\n");


void MWHgethostname(char * MWHHostName, size_t len)
{
#if defined(_WIN32)
  strcpy(MWHHostName, "Windows");
#else
  gethostname(MWHHostName, len);
#endif
}


TimeType GlobalTime;

uint32_t InitializeTiming(TimeType * T)
{
  uint32_t i;
  T->N = 0;
  T->NAloc = 16;
  if ((T->Name = (char**)malloc(T->NAloc * sizeof(char *))) == NULL)
  {
    printf("InitializeTiming failed\n");
    fflush(stdout);
    return 1;
  }
  for (i = 0; i < T->NAloc; i++)
  {
    if ((T->Name[i] = (char*)malloc(MAXStrLength * sizeof(char))) == NULL)
    {
      printf("InitializeTiming failed\n");
      fflush(stdout);
      return 1;
    }
    strcpy(T->Name[i], "\n");
  }
  if ((T->TotalTime = (uint32_t*)malloc(T->NAloc * sizeof(uint32_t))) == NULL)
  {
    printf("InitializeTiming failed\n");
    fflush(stdout);
    return 1;
  }
  if ((T->BeginTime = (time_t*)malloc(T->NAloc * sizeof(time_t))) == NULL)
  {
    printf("InitializeTiming failed\n");
    fflush(stdout);
    return 1;
  }
  if ((T->NumberOfRuns = (uint32_t*)malloc(T->NAloc * sizeof(uint32_t))) == NULL)
  {
    printf("InitializeTiming failed\n");
    fflush(stdout);
    return 1;
  }
  if ((T->NowRunning = (uint32_t*)malloc(T->NAloc * sizeof(uint32_t))) == NULL)
  {
    printf("InitializeTiming failed\n");
    fflush(stdout);
    return 1;
  }
  return 0;
}

uint32_t EnlargeTimingStructure(TimeType * T)
{
  uint32_t i;
  T->NAloc += 16;
  if ((T->Name = (char**)realloc(T->Name, T->NAloc * sizeof(char *))) == NULL)
  {
    printf("EnlargeTimingStructure failed\n");
    fflush(stdout);
    return 1;
  }
  for (i = T->NAloc - 16; i < T->NAloc; i++)
  {
    if ((T->Name[i] = (char*)malloc(MAXStrLength * sizeof(char))) == NULL)
    {
      printf("EnlargeTimingStructure failed\n");
      fflush(stdout);
      return 1;
    }
    strcpy(T->Name[i], "\n");
  }
  if ((T->TotalTime = (uint32_t*)realloc(T->TotalTime, T->NAloc * sizeof(uint32_t))) == NULL)
  {
    printf("EnlargeTimingStructure failed\n");
    fflush(stdout);
    return 1;
  }
  if ((T->BeginTime = (time_t*)realloc(T->BeginTime, T->NAloc * sizeof(time_t))) == NULL)
  {
    printf("EnlargeTimingStructure failed\n");
    fflush(stdout);
    return 1;
  }
  if ((T->NumberOfRuns = (uint32_t*)realloc(T->NumberOfRuns, T->NAloc * sizeof(uint32_t))) == NULL)
  {
    printf("EnlargeTimingStructure failed\n");
    fflush(stdout);
    return 1;
  }
  if ((T->NowRunning = (uint32_t*)realloc(T->NowRunning, T->NAloc * sizeof(uint32_t))) == NULL)
  {
    printf("EnlargeTimingStructure failed\n");
    fflush(stdout);
    return 1;
  }
  return 0;
}

void FreeTiming(TimeType * T)
{
  if (T->NAloc != 0)
  {
    free(T->Name);
    free(T->TotalTime);
    free(T->BeginTime);
    free(T->NumberOfRuns);
    free(T->NowRunning);
  }
  T->N = 0;
  T->NAloc = 0;
}

uint32_t BeginTiming(const char *s, TimeType * T)
{
  uint32_t i;
  for (i = 0; i < T->N; i++)
  {
    if ((strcmp(T->Name[i], s) == 0) && (strcmp(s, T->Name[i]) == 0))
    {
      if (T->NowRunning[i] != 0)
      {
        printf("Trying to start a timing for %s\n", s);
        printf("But timing still runs\n");
        printf("BeginTiming ignored\n");
        return 1;
      }
      time(&(T->BeginTime[i]));
      T->NowRunning[i] = 1;
      return 0;
    }
  }
  if (T->N == T->NAloc)
    if (EnlargeTimingStructure(T))
      return 1;
  strcpy(T->Name[T->N], s);
  T->TotalTime[T->N] = 0;
  time(&(T->BeginTime[T->N]));
  T->NumberOfRuns[T->N] = 0;
  T->NowRunning[T->N] = 1;
  T->N++;
  return 0;
}

uint32_t EndTiming(const char *s, TimeType * T)
{
  uint32_t i;
  time_t now;
  for (i = 0; i < T->N; i++)
  {
    if ((strcmp(T->Name[i], s) == 0) && (strcmp(s, T->Name[i]) == 0))
    {
      if (T->NowRunning[i] == 0)
      {
        printf("Trying to end a timing for %s\n", s);
        printf("But timing already ended\n");
        printf("EndTiming ignored\n");
        return 1;
      }
      else
      {
        time(&now);
        T->TotalTime[i] += (uint32_t) (difftime(now, T->BeginTime[i]));
        T->NumberOfRuns[i]++;
        T->NowRunning[i] = 0;
      }
      return 0;
    }
  }
  printf("Trying to end a timing for %s\n", s);
  printf("But timing never started\n");
  printf("EndTiming ignored\n");
  return 1;
}






/*
*   Prints the values of the operator.
 */
uint32_t StreamPrintTiming(TimeType T, FILE * stream)
{
  uint32_t i;
  fprintf(stream, "Timing results\n");
  fprintf(stream, "   Total_time | NumberOfRuns | Running | Name\n");
  for (i = 0; i < T.N; i++)
  {
    if (T.TotalTime[i] > 86399)
    {
      fprintf(stream, "%4i ", T.TotalTime[i] / 86400);
    }
    else
    {
      fprintf(stream, "     ");
    }
    fprintf(stream, "%2i:", (T.TotalTime[i] / 3600) % 24);
    fprintf(stream, "%02i:", (T.TotalTime[i] / 60) % 60);
    fprintf(stream, "%02i ", (T.TotalTime[i]) % 60);
    fprintf(stream, "| %12i | %7i | %s\n", T.NumberOfRuns[i], T.NowRunning[i], T.Name[i]);
  }
  return 0;
}

uint32_t FilePrintTiming(TimeType T, const char *const flname, const char *const mode)
{
  FILE  *stream;
  uint32_t retval;
  if ((strcmp(mode, "w") && strcmp(mode, "a")) != 0)
  {
    printf
      ("ERROR in FilePrintTiming.\nmode must be either \"w\" for overwrite or \"a\" for append\n %s as input is not excepted\n",
       mode);
    return 1;
  }
  if ((stream = fopen(flname, mode)) == NULL)
  {
    printf("ERROR\nCan't open %s for writing with mode %s.\n", flname, mode);
    return 1;
  }
  retval = StreamPrintTiming(T, stream);
  fclose(stream);
  return retval;
}

uint32_t PrintTiming(TimeType T)
{
  return StreamPrintTiming(T, stdout);
}
