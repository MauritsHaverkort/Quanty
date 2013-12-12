#ifndef __TIMESTAMP__
#define __TIMESTAMP__ "???"
#endif
#ifndef DBL_MIN
#define DBL_MIN 2.3E-308
#endif
#ifndef DBL_EPSILON
#define DBL_EPSILON 2.3E-16
#endif
#define TimeStampAll() time_t MWHrawtime;\
                       struct tm * MWHtimeinfo;\
                       time ( &MWHrawtime );\
                       MWHtimeinfo = localtime ( &MWHrawtime );\
                       printf("Program excicuted on: %s", asctime (MWHtimeinfo) );\
                       char MWHHostName[_POSIX_HOST_NAME_MAX+1];\
                       gethostname(MWHHostName, _POSIX_HOST_NAME_MAX+1);\
                       printf("Running on host: %s\n",MWHHostName);\
                       printf("number of available processors: %i\n",omp_get_num_procs());\
                       printf("maximum number of threads in parallel region: %i\n",omp_get_max_threads());\
                       printf("Compilation of %s\n",__FILE__);\
                       printf("Compiled on %s at %s\n",__DATE__,__TIME__);\
                       printf("File last modified on %s \n",__TIMESTAMP__);\
                       printf("Smallest possitive float:  %E \n",DBL_MIN);\
                       printf("Smallest deviation from 1: %E \n",DBL_EPSILON);


uint32_t InitializeTiming(TimeType * T)
{
  uint32_t i;
  T->N = 0;
  T->NAloc = 16;
  if ((T->Name = malloc(T->NAloc * sizeof(char *))) == NULL)
  {
    printf("InitializeTiming failed\n");
    fflush(stdout);
    return 1;
  }
  for (i = 0; i < T->NAloc; i++)
  {
    if ((T->Name[i] = malloc(MAXStrLength * sizeof(char))) == NULL)
    {
      printf("InitializeTiming failed\n");
      fflush(stdout);
      return 1;
    }
    strcpy(T->Name[i], "\n");
  }
  if ((T->TotalTime = malloc(T->NAloc * sizeof(uint32_t))) == NULL)
  {
    printf("InitializeTiming failed\n");
    fflush(stdout);
    return 1;
  }
  if ((T->BeginTime = malloc(T->NAloc * sizeof(time_t))) == NULL)
  {
    printf("InitializeTiming failed\n");
    fflush(stdout);
    return 1;
  }
  if ((T->NumberOfRuns = malloc(T->NAloc * sizeof(uint32_t))) == NULL)
  {
    printf("InitializeTiming failed\n");
    fflush(stdout);
    return 1;
  }
  if ((T->NowRunning = malloc(T->NAloc * sizeof(uint32_t))) == NULL)
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
  if ((T->Name = realloc(T->Name, T->NAloc * sizeof(char *))) == NULL)
  {
    printf("EnlargeTimingStructure failed\n");
    fflush(stdout);
    return 1;
  }
  for (i = T->NAloc - 16; i < T->NAloc; i++)
  {
    if ((T->Name[i] = malloc(MAXStrLength * sizeof(char))) == NULL)
    {
      printf("EnlargeTimingStructure failed\n");
      fflush(stdout);
      return 1;
    }
    strcpy(T->Name[i], "\n");
  }
  if ((T->TotalTime = realloc(T->TotalTime, T->NAloc * sizeof(uint32_t))) == NULL)
  {
    printf("EnlargeTimingStructure failed\n");
    fflush(stdout);
    return 1;
  }
  if ((T->BeginTime = realloc(T->BeginTime, T->NAloc * sizeof(time_t))) == NULL)
  {
    printf("EnlargeTimingStructure failed\n");
    fflush(stdout);
    return 1;
  }
  if ((T->NumberOfRuns = realloc(T->NumberOfRuns, T->NAloc * sizeof(uint32_t))) == NULL)
  {
    printf("EnlargeTimingStructure failed\n");
    fflush(stdout);
    return 1;
  }
  if ((T->NowRunning = realloc(T->NowRunning, T->NAloc * sizeof(uint32_t))) == NULL)
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
        T->TotalTime[i] += (uint32_t) difftime(now, T->BeginTime[i]);
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
    if (((T.TotalTime[i] / 60) % 60) < 10)
    {
      fprintf(stream, "0%1i:", (T.TotalTime[i] / 60) % 60);
    }
    else
    {
      fprintf(stream, "%2i:", (T.TotalTime[i] / 60) % 60);
    }
    if ((T.TotalTime[i] % 60) < 10)
    {
      fprintf(stream, "0%1i ", T.TotalTime[i] % 60);
    }
    else
    {
      fprintf(stream, "%2i ", T.TotalTime[i] % 60);
    }
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
