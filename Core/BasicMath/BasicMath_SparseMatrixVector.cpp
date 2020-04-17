/*
 *     ( A[0][0] | A[0][1] | A[0][2] )
 *   A=( A[1][0] | A[1][1] | A[1][2] )
 *     ( A[2][0] | A[2][1] | A[2][2] )
 *
 *   (A B)[i][j] = \sum_k A[i][k] B[k][j]
 *
 *   y[i] = \sum_j A[i][j] x[j]
 *
 */


/* Free a sparse matrix the structure M itselve is not freed
 */
void FreeSparseMatrix(SparseMatrixType * const M)
{
  if (M->NAloc)
  {
    free(M->Next);
    free(M->i);
    free(M->j);
    free(M->VRe);
    if (M->QComplex)
      free(M->VIm);
  }
  free(M->KeyList);
  M->Ni = 0;
  M->Nj = 0;
  M->N = 0;
  M->NAloc = 0;
  M->QComplex = 0;
  M->NBitsKey = 0;
}

/* Init sparse matrix
 * M->NBitsKey must be set
 * M->QComplex must be set
 * M->Ni and M->Nj are optional for initialization
 * addvised value 2^NBitsKey ~ number of nonzero elements
 */
uint32_t InitSparseMatrix(SparseMatrixType * const M)
{
  strcpy(M->Name, "Matrix");
  M->N = 0;
  M->NAloc = 1;
  M->Next = (uint32_t*)malloc(M->NAloc * sizeof(uint32_t));
  M->VRe = (double*)malloc(M->NAloc * sizeof(double));
  if (M->QComplex)
    M->VIm = (double*)malloc(M->NAloc * sizeof(double));
  M->i = (uint32_t*)malloc(M->NAloc * sizeof(uint32_t));
  M->j = (uint32_t*)malloc(M->NAloc * sizeof(uint32_t));
  M->KeyList = (uint32_t*)calloc((uint32_t) 1 << (M->NBitsKey), sizeof(uint32_t));
  return 0;
}


/* Init sparse matrix
 * M->NBitsKey must be set
 * M->QComplex must be set
 * M->Ni and M->Nj must be set
 * addvised value 2^NBitsKey ~ number of nonzero elements
 */
uint32_t InitSparseIdentityMatrix(SparseMatrixType * const M)
{
  int    i;
  strcpy(M->Name, "I_Matrix");
  M->N = Min(M->Ni, M->Nj);
  M->NAloc = M->N;
  M->Next = (uint32_t*)malloc(M->NAloc * sizeof(uint32_t));
  M->VRe = (double*)malloc(M->NAloc * sizeof(double));
  if (M->QComplex)
    M->VIm = (double*)calloc(M->NAloc, sizeof(double));
  M->i = (uint32_t*)malloc(M->NAloc * sizeof(uint32_t));
  M->j = (uint32_t*)malloc(M->NAloc * sizeof(uint32_t));
  M->KeyList = (uint32_t*)calloc((uint32_t) 1 << (M->NBitsKey), sizeof(uint32_t));
  if ((M->KeyList == NULL) && (M->j))
  {
    free(M->Next);
    free(M->VRe);
    if (M->QComplex)
      free(M->VIm);
    free(M->i);
    free(M->j);
  }
  for (i = 0; i < M->N; i++)
  {
    M->i[i] = i;
    M->j[i] = i;
    M->VRe[i] = 1;
  }
  return 0;
}

/* turn a real matrix into a complex matrix i.e. allocate memory for the complex part and set this to zero
 */
uint32_t MakeRealSparseMatrixComplex(SparseMatrixType * const M)
{
  if (M->QComplex)
    return 0;
  M->VIm = (double*)calloc(M->NAloc, sizeof(double));
  if (M->VIm == NULL)
    return 1;
  M->QComplex = 1;
  return 0;
}

/*  add a real element to a real sparse matrix 
 *  (No cheks if M is real there will be an error if this is used on a complex matrix)
 *  Function for internal use
 *  if M[i][j] already excist add it to this value
 *  otherwise add a new element to the matrix
 *  enlarge the allocated memory if needed
 */
uint32_t RealSparseMatrixAddElement(SparseMatrixType * const M, uint32_t i, uint32_t j, double VRe)
{
  void  *temp;
  uint32_t NextIndex, ThereIsANext, FirstAlloc = 0, index;
  index = ((i) ^ ((j >> 16) | (j << 16)));
  /* test if the key that belongs to (i,j) already excists */
  if ((M->KeyList[hash(&index, 4, M->NBitsKey)]) == 0)
  {                             /* add a new ellement to O with a unique key */
    /* Test if one has to reallocate M->i, M->j and M->V in order to store the new element */
    if (M->N >= M->NAloc)
    {
      if (M->NAloc == 0)
        FirstAlloc = 1;
      M->NAloc += M->Ni + M->Nj;
      temp = M->i;
      if (FirstAlloc)
      {
        M->i = (uint32_t*)malloc(M->NAloc * sizeof(uint32_t));
      }
      else
      {
        M->i = (uint32_t*)realloc(M->i, M->NAloc * sizeof(uint32_t));
      }
      if (M->i == NULL)
      {                         /* realloc failed */
        printf("alloc failed OperatorAdd 01i\n");
        fflush(stdout);
        /* restore O and return 0 */
        M->i = (uint32_t*)temp;
        M->NAloc -= M->Ni + M->Nj;
        return (1);
      }
      temp = M->j;
      if (FirstAlloc)
      {
        M->j = (uint32_t*)malloc(M->NAloc * sizeof(uint32_t));
      }
      else
      {
        M->j = (uint32_t*)realloc(M->j, M->NAloc * sizeof(uint32_t));
      }
      if (M->j == NULL)
      {                         /* realloc failed */
        printf("alloc failed OperatorAdd 01j\n");
        fflush(stdout);
        /* restore O and return 0 */
        M->j = (uint32_t*)temp;
        M->NAloc -= M->Ni + M->Nj;
        M->i = (uint32_t*)realloc(M->i, M->NAloc * sizeof(uint32_t));
        return (1);
      }
      temp = M->VRe;
      if (FirstAlloc)
      {
        M->VRe = (double*)malloc(M->NAloc * sizeof(double));
      }
      else
      {
        M->VRe = (double*)realloc(M->VRe, M->NAloc * sizeof(double));
      }
      if (M->VRe == NULL)
      {                         /* realloc failed */
        printf("alloc failed OperatorAdd 02\n");
        fflush(stdout);
        /* restore O and return 0 */
        M->VRe = (double*)temp;
        M->NAloc -= M->Ni + M->Nj;
        M->i = (uint32_t*)realloc(M->i, M->NAloc * sizeof(uint32_t));
        M->j = (uint32_t*)realloc(M->j, M->NAloc * sizeof(uint32_t));
        return (1);
      }
      temp = M->Next;
      if (FirstAlloc)
      {
        M->Next = (uint32_t*)malloc(M->NAloc * sizeof(uint32_t));
      }
      else
      {
        M->Next = (uint32_t*)realloc(M->Next, M->NAloc * sizeof(uint32_t));
      }
      if (M->Next == NULL)
      {                         /* realloc failed */
        printf("alloc failed OperatorAdd 03\n");
        fflush(stdout);
        /* restore O and return 0 */
        M->Next = (uint32_t*)temp;
        M->NAloc -= M->Ni + M->Nj;
        M->i = (uint32_t*)realloc(M->i, M->NAloc * sizeof(uint32_t));
        M->j = (uint32_t*)realloc(M->j, M->NAloc * sizeof(uint32_t));
        M->VRe = (double*)realloc(M->VRe, M->NAloc * sizeof(double));
        return (1);
      }
    }
    M->i[M->N] = i;
    M->j[M->N] = j;
    M->VRe[M->N] = VRe;
    M->Next[M->N] = 0;
    M->N++;
    M->KeyList[hash(&index, 4, M->NBitsKey)] = M->N;
  }
  else
  {                             /* the key is already taken, test if this is due to a collision or if i is already in the opperator */
    ThereIsANext = 1;
    NextIndex = M->KeyList[hash(&index, 4, M->NBitsKey)] - 1;
    while (ThereIsANext)
    {
      if ((i == M->i[NextIndex]) && (j == M->j[NextIndex]))
      {                         /* i,j already in the opperator, add and leave loop */
        ThereIsANext = 0;
        M->VRe[NextIndex] += VRe;
      }
      else
      {                         /* i,j not equal to the i,j found which has the same index in the keylist, if there is a next continue, otherwise add the element */
        if (!(M->Next[NextIndex]))
        {                       /* add a new ellemnt to O with already used key */
          /* Test if one has to reallocate M->i and M->V in order to store the new element */
          if (M->N >= M->NAloc)
          {
            M->NAloc += M->Ni + M->Nj;
            temp = M->i;
            M->i = (uint32_t*)realloc(M->i, M->NAloc * sizeof(uint32_t));
            if (M->i == NULL)
            {                   /* realloc failed */
              printf("alloc failed OperatorAdd 11i\n");
              fflush(stdout);
              /* restore O and return 0 */
              M->i = (uint32_t*)temp;
              M->NAloc -= M->Ni + M->Nj;
              return (1);
            }
            temp = M->j;
            M->j = (uint32_t*)realloc(M->j, M->NAloc * sizeof(uint32_t));
            if (M->j == NULL)
            {                   /* realloc failed */
              printf("alloc failed OperatorAdd 11j\n");
              fflush(stdout);
              /* restore O and return 0 */
              M->j = (uint32_t*)temp;
              M->NAloc -= M->Ni + M->Nj;
              M->i = (uint32_t*)realloc(M->i, M->NAloc * sizeof(uint32_t));
              return (1);
            }
            temp = M->VRe;
            M->VRe = (double*)realloc(M->VRe, M->NAloc * sizeof(double));
            if (M->VRe == NULL)
            {                   /* realloc failed */
              printf("alloc failed OperatorAdd 12\n");
              fflush(stdout);
              /* restore O and return 0 */
              M->VRe = (double*)temp;
              M->NAloc -= M->Ni + M->Nj;
              M->i = (uint32_t*)realloc(M->i, M->NAloc * sizeof(uint32_t));
              M->j = (uint32_t*)realloc(M->i, M->NAloc * sizeof(uint32_t));
              return (1);
            }
            temp = M->Next;
            M->Next = (uint32_t*)realloc(M->Next, M->NAloc * sizeof(uint32_t));
            if (M->Next == NULL)
            {                   /* realloc failed */
              printf("alloc failed OperatorAdd 13\n");
              fflush(stdout);
              /* restore O and return 0 */
              M->Next = (uint32_t*)temp;
              M->NAloc -= M->Ni + M->Nj;
              M->i = (uint32_t*)realloc(M->i, M->NAloc * sizeof(uint32_t));
              M->j = (uint32_t*)realloc(M->j, M->NAloc * sizeof(uint32_t));
              M->VRe = (double*)realloc(M->VRe, M->NAloc * sizeof(double));
              return (1);
            }
          }
          M->i[M->N] = i;
          M->j[M->N] = j;
          M->VRe[M->N] = VRe;
          M->Next[M->N] = NextIndex + 1;
          M->N++;
          M->KeyList[hash(&index, 4, M->NBitsKey)] = M->N;
          ThereIsANext = 0;
        }
        else
        {                       /* The key is in use and there are more opperators with this key */
          NextIndex = M->Next[NextIndex] - 1;
        }
      }
    }
  }
  return 0;
}

/*  add a complex element to a complex sparse matrix 
 *  (No cheks if M is complex there will be an error if this is used on a real matrix)
 *  Function for internal use
 *  if M[i][j] already excist add it to this value
 *  otherwise add a new element to the matrix
 *  enlarge the allocated memory if needed
 */
uint32_t ComplexSparseMatrixAddElement(SparseMatrixType * const M, uint32_t i, uint32_t j, double VRe, double VIm)
{
  void  *temp;
  uint32_t NextIndex, ThereIsANext, FirstAlloc = 0, index;
  index = ((i) ^ ((j >> 16) | (j << 16)));
  /* test if the key that belongs to (i,j) already excists */
  if ((M->KeyList[hash(&index, 4, M->NBitsKey)]) == 0)
  {                             /* add a new ellement to O with a unique key */
    /* Test if one has to reallocate M->i, M->j and M->V in order to store the new element */
    if (M->N >= M->NAloc)
    {
      if (M->NAloc == 0)
        FirstAlloc = 1;
      M->NAloc += M->Ni + M->Nj;
      temp = M->i;
      if (FirstAlloc)
      {
        M->i = (uint32_t*)malloc(M->NAloc * sizeof(uint32_t));
      }
      else
      {
        M->i = (uint32_t*)realloc(M->i, M->NAloc * sizeof(uint32_t));
      }
      if (M->i == NULL)
      {                         /* realloc failed */
        printf("alloc failed OperatorAdd 01i\n");
        fflush(stdout);
        /* restore O and return 0 */
        M->i = (uint32_t*)temp;
        M->NAloc -= M->Ni + M->Nj;
        return (1);
      }
      temp = M->j;
      if (FirstAlloc)
      {
        M->j = (uint32_t*)malloc(M->NAloc * sizeof(uint32_t));
      }
      else
      {
        M->j = (uint32_t*)realloc(M->j, M->NAloc * sizeof(uint32_t));
      }
      if (M->j == NULL)
      {                         /* realloc failed */
        printf("alloc failed OperatorAdd 01j\n");
        fflush(stdout);
        /* restore O and return 0 */
        M->j = (uint32_t*)temp;
        M->NAloc -= M->Ni + M->Nj;
        M->i = (uint32_t*)realloc(M->i, M->NAloc * sizeof(uint32_t));
        return (1);
      }
      temp = M->VRe;
      if (FirstAlloc)
      {
        M->VRe = (double*)malloc(M->NAloc * sizeof(double));
      }
      else
      {
        M->VRe = (double*)realloc(M->VRe, M->NAloc * sizeof(double));
      }
      if (M->VRe == NULL)
      {                         /* realloc failed */
        printf("alloc failed OperatorAdd 02\n");
        fflush(stdout);
        /* restore O and return 0 */
        M->VRe = (double*)temp;
        M->NAloc -= M->Ni + M->Nj;
        M->i = (uint32_t*)realloc(M->i, M->NAloc * sizeof(uint32_t));
        M->j = (uint32_t*)realloc(M->j, M->NAloc * sizeof(uint32_t));
        return (1);
      }
      temp = M->VIm;
      if (FirstAlloc)
      {
        M->VIm = (double*)malloc(M->NAloc * sizeof(double));
      }
      else
      {
        M->VIm = (double*)realloc(M->VIm, M->NAloc * sizeof(double));
      }
      if (M->VIm == NULL)
      {                         /* realloc failed */
        printf("alloc failed OperatorAdd 02\n");
        fflush(stdout);
        /* restore O and return 0 */
        M->VIm = (double*)temp;
        M->NAloc -= M->Ni + M->Nj;
        M->i = (uint32_t*)realloc(M->i, M->NAloc * sizeof(uint32_t));
        M->j = (uint32_t*)realloc(M->j, M->NAloc * sizeof(uint32_t));
        M->VRe = (double*)realloc(M->VRe, M->NAloc * sizeof(double));
        return (1);
      }
      temp = M->Next;
      if (FirstAlloc)
      {
        M->Next = (uint32_t*)malloc(M->NAloc * sizeof(uint32_t));
      }
      else
      {
        M->Next = (uint32_t*)realloc(M->Next, M->NAloc * sizeof(uint32_t));
      }
      if (M->Next == NULL)
      {                         /* realloc failed */
        printf("alloc failed OperatorAdd 03\n");
        fflush(stdout);
        /* restore O and return 0 */
        M->Next = (uint32_t*)temp;
        M->NAloc -= M->Ni + M->Nj;
        M->i = (uint32_t*)realloc(M->i, M->NAloc * sizeof(uint32_t));
        M->j = (uint32_t*)realloc(M->j, M->NAloc * sizeof(uint32_t));
        M->VRe = (double*)realloc(M->VRe, M->NAloc * sizeof(double));
        M->VIm = (double*)realloc(M->VIm, M->NAloc * sizeof(double));
        return (1);
      }
    }
    M->i[M->N] = i;
    M->j[M->N] = j;
    M->VRe[M->N] = VRe;
    M->VIm[M->N] = VIm;
    M->Next[M->N] = 0;
    M->N++;
    M->KeyList[hash(&index, 4, M->NBitsKey)] = M->N;
  }
  else
  {                             /* the key is already taken, test if this is due to a collision or if i is already in the opperator */
    ThereIsANext = 1;
    NextIndex = M->KeyList[hash(&index, 4, M->NBitsKey)] - 1;
    while (ThereIsANext)
    {
      if ((i == M->i[NextIndex]) && (j == M->j[NextIndex]))
      {                         /* i,j already in the opperator, add and leave loop */
        ThereIsANext = 0;
        M->VRe[NextIndex] += VRe;
        M->VIm[NextIndex] += VIm;
      }
      else
      {                         /* i,j not equal to the i,j found which has the same index in the keylist, if there is a next continue, otherwise add the element */
        if (!(M->Next[NextIndex]))
        {                       /* add a new ellemnt to O with already used key */
          /* Test if one has to reallocate M->i and M->V in order to store the new element */
          if (M->N >= M->NAloc)
          {
            M->NAloc += M->Ni + M->Nj;
            temp = M->i;
            M->i = (uint32_t*)realloc(M->i, M->NAloc * sizeof(uint32_t));
            if (M->i == NULL)
            {                   /* realloc failed */
              printf("alloc failed OperatorAdd 11i\n");
              fflush(stdout);
              /* restore O and return 0 */
              M->i = (uint32_t*)temp;
              M->NAloc -= M->Ni + M->Nj;
              return (1);
            }
            temp = M->j;
            M->j = (uint32_t*)realloc(M->j, M->NAloc * sizeof(uint32_t));
            if (M->j == NULL)
            {                   /* realloc failed */
              printf("alloc failed OperatorAdd 11j\n");
              fflush(stdout);
              /* restore O and return 0 */
              M->j = (uint32_t*)temp;
              M->NAloc -= M->Ni + M->Nj;
              M->i = (uint32_t*)realloc(M->i, M->NAloc * sizeof(uint32_t));
              return (1);
            }
            temp = M->VRe;
            M->VRe = (double*)realloc(M->VRe, M->NAloc * sizeof(double));
            if (M->VRe == NULL)
            {                   /* realloc failed */
              printf("alloc failed OperatorAdd 12\n");
              fflush(stdout);
              /* restore O and return 0 */
              M->VRe = (double*)temp;
              M->NAloc -= M->Ni + M->Nj;
              M->i = (uint32_t*)realloc(M->i, M->NAloc * sizeof(uint32_t));
              M->j = (uint32_t*)realloc(M->i, M->NAloc * sizeof(uint32_t));
              return (1);
            }
            temp = M->VIm;
            M->VIm = (double*)realloc(M->VIm, M->NAloc * sizeof(double));
            if (M->VIm == NULL)
            {                   /* realloc failed */
              printf("alloc failed OperatorAdd 12\n");
              fflush(stdout);
              /* restore O and return 0 */
              M->VIm = (double*)temp;
              M->NAloc -= M->Ni + M->Nj;
              M->i = (uint32_t*)realloc(M->i, M->NAloc * sizeof(uint32_t));
              M->j = (uint32_t*)realloc(M->i, M->NAloc * sizeof(uint32_t));
              M->VRe = (double*)realloc(M->VRe, M->NAloc * sizeof(double));
              return (1);
            }
            temp = M->Next;
            M->Next = (uint32_t*)realloc(M->Next, M->NAloc * sizeof(uint32_t));
            if (M->Next == NULL)
            {                   /* realloc failed */
              printf("alloc failed OperatorAdd 13\n");
              fflush(stdout);
              /* restore O and return 0 */
              M->Next = (uint32_t*)temp;
              M->NAloc -= M->Ni + M->Nj;
              M->i = (uint32_t*)realloc(M->i, M->NAloc * sizeof(uint32_t));
              M->j = (uint32_t*)realloc(M->j, M->NAloc * sizeof(uint32_t));
              M->VRe = (double*)realloc(M->VRe, M->NAloc * sizeof(double));
              M->VIm = (double*)realloc(M->VIm, M->NAloc * sizeof(double));
              return (1);
            }
          }
          M->i[M->N] = i;
          M->j[M->N] = j;
          M->VRe[M->N] = VRe;
          M->VIm[M->N] = VIm;
          M->Next[M->N] = NextIndex + 1;
          M->N++;
          M->KeyList[hash(&index, 4, M->NBitsKey)] = M->N;
          ThereIsANext = 0;
        }
        else
        {                       /* The key is in use and there are more opperators with this key */
          NextIndex = M->Next[NextIndex] - 1;
        }
      }
    }
  }
  return 0;
}

/*  add an element to a sparse matrix
 *  for real matrices take VIm=0
 *  Function for internal use
 *  if M[i][j] already excist add it to this value
 *  otherwise add a new element to the matrix
 *  enlarge the allocated memory if needed
 */
uint32_t SparseMatrixAddElement(SparseMatrixType * const M, uint32_t i, uint32_t j, double VRe, double VIm)
{
  if ((M->QComplex == 0) && (VIm == 0))
  {
    return RealSparseMatrixAddElement(M, i, j, VRe);
  }
  if (M->QComplex == 0)
  {
    if (MakeRealSparseMatrixComplex(M))
      return 1;
    return ComplexSparseMatrixAddElement(M, i, j, VRe, VIm);
  }
  return ComplexSparseMatrixAddElement(M, i, j, VRe, VIm);
}

/*  add several elements to a sparse matrix starting at position i,j
 *
 */
uint32_t SparseMatrixAddMatrix(SparseMatrixType * const M, uint32_t i, uint32_t j, MatrixType A, double const epsilon)
{
  uint32_t ai, aj;
  const double epsilonsqr = Sqr(epsilon);
  if (i + A.Ni > M->Ni)
  {
    printf("Trid to add elements outside boundary of SparseMatrix in SparseMatrixAddMatrix\n");
    return 1;
  }
  if (j + A.Nj > M->Nj)
  {
    printf("Trid to add elements outside boundary of SparseMatrix in SparseMatrixAddMatrix\n");
    return 1;
  }
  if (M->QComplex)
  {
    if (A.QComplex)
    {
      for (ai = 0; ai < A.Ni; ai++)
      {
        for (aj = 0; aj < A.Nj; aj++)
        {
          if ((Sqr(A.VRe[ai][aj]) + Sqr(A.VIm[ai][aj])) > epsilonsqr)
          {
            if (ComplexSparseMatrixAddElement(M, i + ai, j + aj, A.VRe[ai][aj], A.VIm[ai][aj]))
            {
              printf("ComplexSparseMatrixAddElement failed in SparseMatrixAddMatrix\n");
              return 1;
            }
          }
        }
      }
    }
    else
    {
      for (ai = 0; ai < A.Ni; ai++)
      {
        for (aj = 0; aj < A.Nj; aj++)
        {
          if (Sqr(A.VRe[ai][aj]) > epsilonsqr)
          {
            if (ComplexSparseMatrixAddElement(M, i + ai, j + aj, A.VRe[ai][aj], 0.0))
            {
              printf("ComplexSparseMatrixAddElement failed in SparseMatrixAddMatrix\n");
              return 1;
            }
          }
        }
      }
    }
  }
  else
  {
    if (A.QComplex)
    {
      if (MakeRealSparseMatrixComplex(M))
        return 1;
      for (ai = 0; ai < A.Ni; ai++)
      {
        for (aj = 0; aj < A.Nj; aj++)
        {
          if ((Sqr(A.VRe[ai][aj]) + Sqr(A.VIm[ai][aj])) > epsilonsqr)
          {
            if (ComplexSparseMatrixAddElement(M, i + ai, j + aj, A.VRe[ai][aj], A.VIm[ai][aj]))
            {
              printf("ComplexSparseMatrixAddElement failed in SparseMatrixAddMatrix\n");
              return 1;
            }
          }
        }
      }
    }
    else
    {
      for (ai = 0; ai < A.Ni; ai++)
      {
        for (aj = 0; aj < A.Nj; aj++)
        {
          if (Sqr(A.VRe[ai][aj]) > epsilonsqr)
          {
            if (RealSparseMatrixAddElement(M, i + ai, j + aj, A.VRe[ai][aj]))
            {
              printf("RealSparseMatrixAddElement failed in SparseMatrixAddMatrix\n");
              return 1;
            }
          }
        }
      }
    }
  }
  return 0;
}


/*  add several elements to a sparse matrix starting at position i,j
 *
 */
uint32_t SparseMatrixAddMatrixTransposed(SparseMatrixType * const M, uint32_t i, uint32_t j, MatrixType A,
                                         double const epsilon)
{
  uint32_t ai, aj;
  const double epsilonsqr = Sqr(epsilon);
  if (i + A.Nj > M->Ni)
  {
    printf("Trid to add elements outside boundary of SparseMatrix in SparseMatrixAddMatrixTransposed\n");
    return 1;
  }
  if (j + A.Ni > M->Nj)
  {
    printf("Trid to add elements outside boundary of SparseMatrix in SparseMatrixAddMatrixTransposed\n");
    return 1;
  }
  if (M->QComplex)
  {
    if (A.QComplex)
    {
      for (ai = 0; ai < A.Nj; ai++)
      {
        for (aj = 0; aj < A.Ni; aj++)
        {
          if ((Sqr(A.VRe[aj][ai]) + Sqr(A.VIm[aj][ai])) > epsilonsqr)
          {
            if (ComplexSparseMatrixAddElement(M, i + ai, j + aj, A.VRe[aj][ai], A.VIm[aj][ai]))
            {
              printf("ComplexSparseMatrixAddElement failed in SparseMatrixAddMatrixTransposed\n");
              return 1;
            }
          }
        }
      }
    }
    else
    {
      for (ai = 0; ai < A.Nj; ai++)
      {
        for (aj = 0; aj < A.Ni; aj++)
        {
          if (Sqr(A.VRe[aj][ai]) > epsilonsqr)
          {
            if (ComplexSparseMatrixAddElement(M, i + ai, j + aj, A.VRe[aj][ai], 0.0))
            {
              printf("ComplexSparseMatrixAddElement failed in SparseMatrixAddMatrixTransposed\n");
              return 1;
            }
          }
        }
      }
    }
  }
  else
  {
    if (A.QComplex)
    {
      if (MakeRealSparseMatrixComplex(M))
        return 1;
      for (ai = 0; ai < A.Nj; ai++)
      {
        for (aj = 0; aj < A.Ni; aj++)
        {
          if ((Sqr(A.VRe[aj][ai]) + Sqr(A.VIm[aj][ai])) > epsilonsqr)
          {
            if (ComplexSparseMatrixAddElement(M, i + ai, j + aj, A.VRe[aj][ai], A.VIm[aj][ai]))
            {
              printf("ComplexSparseMatrixAddElement failed in SparseMatrixAddMatrixTransposed\n");
              return 1;
            }
          }
        }
      }
    }
    else
    {
      for (ai = 0; ai < A.Nj; ai++)
      {
        for (aj = 0; aj < A.Ni; aj++)
        {
          if (Sqr(A.VRe[aj][ai]) > epsilonsqr)
          {
            if (RealSparseMatrixAddElement(M, i + ai, j + aj, A.VRe[aj][ai]))
            {
              printf("RealSparseMatrixAddElement failed in SparseMatrixAddMatrixTransposed\n");
              return 1;
            }
          }
        }
      }
    }
  }
  return 0;
}


/*  add several elements to a sparse matrix starting at position i,j
 *
 */
uint32_t SparseMatrixAddTridiagonalMatrix(SparseMatrixType * const M, uint32_t i, uint32_t j, TriDiagonalMatrixType A)
{
  uint32_t ai;
  if (i + A.N > M->Ni)
  {
    printf("Tried to add elements outside boundary of SparseMatrix in SparseMatrixAddTridiagonalMatrix\n");
    return 1;
  }
  if (j + A.N > M->Nj)
  {
    printf("Tried to add elements outside boundary of SparseMatrix in SparseMatrixAddTridiagonalMatrix\n");
    return 1;
  }
  if (M->QComplex)
  {

  }
  else
  {
    for (ai = 0; ai < A.N - 1; ai++)
    {
      if (RealSparseMatrixAddElement(M, i + ai, j + ai, A.a[ai + 1]))
      {
        printf("RealSparseMatrixAddElement failed in SparseMatrixAddTridiagonalMatrix\n");
        return 1;
      }
      if (RealSparseMatrixAddElement(M, i + ai + 1, j + ai, A.b[ai + 1]))
      {
        printf("RealSparseMatrixAddElement failed in SparseMatrixAddTridiagonalMatrix\n");
        return 1;
      }
      if (RealSparseMatrixAddElement(M, i + ai, j + ai + 1, A.b[ai + 1]))
      {
        printf("RealSparseMatrixAddElement failed in SparseMatrixAddTridiagonalMatrix\n");
        return 1;
      }
    }
    if (RealSparseMatrixAddElement(M, i + A.N - 1, j + A.N - 1, A.a[A.N]))
    {
      printf("RealSparseMatrixAddElement failed in SparseMatrixAddTridiagonalMatrix\n");
      return 1;
    }
  }
  return 0;
}

/*  Set the values of VRe and VIm equal to M[i][j]
 */
void SparseMatrixExtractElement(SparseMatrixType const M, uint32_t const i, uint32_t const j, double *VRe, double *VIm)
{
  uint32_t NextIndex, ThereIsANext, index;
  index = ((i) ^ ((j >> 16) | (j << 16)));
  /* test if the key that belongs to (i,j) already excists */
  if ((M.KeyList[hash(&index, 4, M.NBitsKey)]) == 0)
  {                             /* element [i][j] does not excist return 0 */
    *VRe = 0;
    *VIm = 0;
    return;
  }
  else
  {                             /* the key is already taken, test if this is due to a collision or if i is already in the opperator */
    ThereIsANext = 1;
    NextIndex = M.KeyList[hash(&index, 4, M.NBitsKey)] - 1;
    while (ThereIsANext)
    {
      if ((i == M.i[NextIndex]) && (j == M.j[NextIndex]))
      {                         /* i,j found, return the values */
        ThereIsANext = 0;
        *VRe = M.VRe[NextIndex];
        *VIm = 0;
        if (M.QComplex)
          *VIm = M.VIm[NextIndex];
        return;
      }
      else
      {                         /* i,j not equal to the i,j found which has the same index in the keylist, if there is a next continue, otherwise return 0 */
        if (!(M.Next[NextIndex]))
        {                       /* return 0 */
          ThereIsANext = 0;
          *VRe = 0;
          *VIm = 0;
          return;
        }
        else
        {                       /* The key is in use and there are more opperators with this key */
          NextIndex = M.Next[NextIndex] - 1;
        }
      }
    }
  }
}


/*  realloc the memory such that NAloc=N
 */
void MinimalAlocSpareseMatrix(SparseMatrixType * const M)
{
  if (M->NAloc == 0)
    return;
  if (M->NAloc == M->N)
    return;
  M->NAloc = M->N;
  if (M->NAloc == 0)
    M->NAloc = 1;
  M->Next = (uint32_t*)realloc(&M->Next, M->NAloc * sizeof(uint32_t));
  M->VRe = (double*)realloc(&M->VRe, M->NAloc * sizeof(double));
  if (M->QComplex)
    M->VIm = (double*)realloc(&M->VIm, M->NAloc * sizeof(double));
  M->i = (uint32_t*)realloc(&M->i, M->NAloc * sizeof(uint32_t));
  M->j = (uint32_t*)realloc(&M->j, M->NAloc * sizeof(uint32_t));
  return;
}

/*  Remove elements for wich |M->V|<epsilon
 */
void RemoveZeroSpareseMatrix(SparseMatrixType * const M, double epsilon)
{
  uint32_t i, Nold;
  double Sqrepsilon;
  for (i = 0; i < (uint32_t) 1 << (M->NBitsKey); i++)
  {                             /* wipe out entire key-list */
    M->KeyList[i] = 0;
  }
  Nold = M->N;
  M->N = 0;
  if (M->QComplex)
  {
    Sqrepsilon = Sqr(epsilon);
    for (i = 0; i < Nold; i++)
    {                           /* add elements back */
      if ((Sqr(M->VRe[i]) + Sqr(M->VIm[i])) > Sqrepsilon)
      {
        ComplexSparseMatrixAddElement(M, M->i[i], M->j[i], M->VRe[i], M->VIm[i]);
      }
    }
  }
  else
  {
    for (i = 0; i < Nold; i++)
    {                           /* add elements back */
      if (fabs(M->VRe[i]) > epsilon)
      {
        RealSparseMatrixAddElement(M, M->i[i], M->j[i], M->VRe[i]);
      }
    }
  }
  MinimalAlocSpareseMatrix(M);
}

/* set elements with an absolute value smaller then epsilon to zero
 */
void SparseMatrixChop(SparseMatrixType * const M, double const epsilon)
{
  RemoveZeroSpareseMatrix(M, epsilon);
}

/*
 * M=M^{T}
 */
void SparseMatrixTranspose(SparseMatrixType * const M)
{
  uint32_t i, Nold;
  for (i = 0; i < (uint32_t) 1 << (M->NBitsKey); i++)
  {                             /* wipe out entire key-list */
    M->KeyList[i] = 0;
  }
  Nold = M->N;
  M->N = 0;
  if (M->QComplex)
  {
    for (i = 0; i < Nold; i++)
    {                           /* add elements back */
      if (((M->VRe[i]) != 0) & ((M->VIm[i] != 0)))
      {
        ComplexSparseMatrixAddElement(M, M->j[i], M->i[i], M->VRe[i], M->VIm[i]);
      }
    }
  }
  else
  {
    for (i = 0; i < Nold; i++)
    {                           /* add elements back */
      if ((M->VRe[i]) != 0)
      {
        RealSparseMatrixAddElement(M, M->j[i], M->i[i], M->VRe[i]);
      }
    }
  }
  MinimalAlocSpareseMatrix(M);
}

/*
 * M=M^{*}
 */
void SparseMatrixComplexConjugate(SparseMatrixType * const M)
{
  if (M->QComplex)
  {
    uint32_t i;
    for (i = 0; i < M->N; i++)
    {
      M->VIm[i] = -M->VIm[i];
    }
  }
}

/*
 * M=M^{\dag}
 */
void SparseMatrixConjugateTranspose(SparseMatrixType * const M)
{
  uint32_t i, Nold;
  for (i = 0; i < (uint32_t) 1 << (M->NBitsKey); i++)
  {                             /* wipe out entire key-list */
    M->KeyList[i] = 0;
  }
  Nold = M->N;
  M->N = 0;
  if (M->QComplex)
  {
    for (i = 0; i < Nold; i++)
    {                           /* add elements back */
      if (((M->VRe[i]) != 0) & ((M->VIm[i] != 0)))
      {
        ComplexSparseMatrixAddElement(M, M->j[i], M->i[i], M->VRe[i], -M->VIm[i]);
      }
    }
  }
  else
  {
    for (i = 0; i < Nold; i++)
    {                           /* add elements back */
      if ((M->VRe[i]) != 0)
      {
        RealSparseMatrixAddElement(M, M->j[i], M->i[i], M->VRe[i]);
      }
    }
  }
  MinimalAlocSpareseMatrix(M);
}

/*  Inverts a sparse matrix, and places the result in the same matrix
 *  It might be better to convert the sparse matrix before to a full matrix and invert the full matrix
 *  As an inverse is often not sparse
 *  all values for wich |V|<epsilon are not stored (set ot zero)
 */
double SparseMatrixInverse(SparseMatrixType * const M, double epsilon)
{
  MatrixType Mfull;
  double ret, epssqr;
  int    i, j, k;
  Mfull.Ni = M->Ni;
  Mfull.Nj = M->Nj;
  Mfull.QComplex = M->QComplex;
  if (InitMatrix(&Mfull))
    return -1;
  for (k = 0; k < M->N; k++)
  {
    Mfull.VRe[M->i[k]][M->j[k]] = M->VRe[k];
  }
  if (Mfull.QComplex)
  {
    for (k = 0; k < M->N; k++)
    {
      Mfull.VIm[M->i[k]][M->j[k]] = M->VIm[k];
    }
  }
  ret = MatrixInverse(&Mfull);
  if (ret < 0)
  {
    FreeMatrix(&Mfull);
    return ret;
  }
  for (i = 0; i < (uint32_t) 1 << (M->NBitsKey); i++)
  {                             /* wipe out entire key-list */
    M->KeyList[i] = 0;
  }
  M->N = 0;
  if (M->QComplex)
  {                             /* Complex */
    epssqr = Sqr(epsilon);
    for (i = 0; i < M->Ni; i++)
    {
      for (j = 0; j < M->Nj; j++)
      {
        if ((Sqr(Mfull.VRe[i][j]) + Sqr(Mfull.VIm[i][j])) > epssqr)
        {
          if (ComplexSparseMatrixAddElement(M, i, j, Mfull.VRe[i][j], Mfull.VIm[i][j]))
          {
            FreeMatrix(&Mfull);
            FreeSparseMatrix(M);
            return -2;
          }
        }
      }
    }
  }
  else
  {                             /*Real */
    for (i = 0; i < M->Ni; i++)
    {
      for (j = 0; j < M->Nj; j++)
      {
        if (fabs(Mfull.VRe[i][j]) > epsilon)
        {
          if (RealSparseMatrixAddElement(M, i, j, Mfull.VRe[i][j]))
          {
            FreeMatrix(&Mfull);
            FreeSparseMatrix(M);
            return -2;
          }
        }
      }
    }
  }
  return ret;
}

/*
 *  M=M + alpha B
 */
uint32_t SparseMatrixAdd(SparseMatrixType * const M, const SparseMatrixType B, double const alphaRe, double const alphaIm)
{
  if ((M->Ni != B.Ni) || (M->Nj != B.Nj))
  {
    printf("Can only add matrices of the same dimension M+B failed\n");
    printf("M.Ni= %i M.Nj= %i B.Ni= %i B.Nj= %i\n", M->Ni, M->Nj, B.Ni, B.Nj);
    fflush(stdout);
    return 1;
  }
  int    i;
  if ((M->QComplex == 0) && (B.QComplex == 0) && (alphaIm == 0))
  {                             /* all real */
    for (i = 0; i < B.N; i++)
    {
      if (RealSparseMatrixAddElement(M, B.i[i], B.j[i], alphaRe * B.VRe[i]))
        return 1;
    }
    return 0;
  }
  if (M->QComplex == 0)
    if (MakeRealSparseMatrixComplex(M))
      return 1;
  if (B.QComplex == 0)
  {                             /* B real, rest complex */
    for (i = 0; i < B.N; i++)
    {
      if (ComplexSparseMatrixAddElement(M, B.i[i], B.j[i], alphaRe * B.VRe[i], alphaIm * B.VRe[i]))
        return 1;
    }
    return 0;
  }
  /* all complex */
  for (i = 0; i < B.N; i++)
  {
    if (ComplexSparseMatrixAddElement
        (M, B.i[i], B.j[i], alphaRe * B.VRe[i] - alphaIm * B.VIm[i], alphaIm * B.VRe[i] + alphaRe * B.VIm[i]))
      return 1;
  }
  return 0;
}


/* 
 *  M = alpha M
 */
uint32_t SparseMatrixScalarMultiply(SparseMatrixType * const M, double const alphaRe, double const alphaIm)
{
  int    i;
  if ((M->QComplex == 0) && (alphaIm == 0))
  {                             /* real */
    for (i = 0; i < M->N; i++)
    {
      M->VRe[i] *= alphaRe;
    }
    return 0;
  }
  if (M->QComplex == 0)
    if (MakeRealSparseMatrixComplex(M))
      return 1;
  double temp;
  /* all complex */
  for (i = 0; i < M->N; i++)
  {
    temp = M->VRe[i];
    M->VRe[i] = alphaRe * M->VRe[i] - alphaIm * M->VIm[i];
    M->VIm[i] = alphaRe * M->VIm[i] + alphaIm * temp;
  }
  return 0;
}

/* M = T^* A T^{transpose}
 * M[i][j]=\sum_k,l T[i][l]^* A[l][k] T[j][k]
 * (A.Ni==A.Nj==T.Nj)
 * M.Ni=M.Nj=T.Ni;
 * M on entrance not initiallized
 */
uint32_t SparseMatrixTransform(SparseMatrixType * const M, SparseMatrixType const A, SparseMatrixType const T, double epsilon)
{
  if ((A.Ni != A.Nj) || (A.Ni != T.Nj))
  {
    printf("Invallid dimensions for matrix transformation A.Ni should be T.Nj and A should be square\n");
    printf("A.Ni=%i A.Nj=%i T.Ni=%i T.Nj=%i \n", A.Ni, A.Nj, T.Ni, T.Nj);
    return 1;
  }
  uint32_t k, l, m;
  double tempRe, tempIm;
  double const epssqr = Sqr(epsilon);
  M->Ni = T.Ni;
  M->Nj = T.Ni;
  M->QComplex = 1;
  if ((A.QComplex == 0) && (T.QComplex == 0))
    M->QComplex = 0;
  if (InitSparseMatrix(M))
    return 1;
  strcpy(M->Name, "U^* A U^T U=(");
  strncat(M->Name, T.Name, Max(0, MAXStrLength - 24) / 2);
  strcat(M->Name, ") A=(");
  strncat(M->Name, A.Name, Max(0, MAXStrLength - 24) / 2);
  strcat(M->Name, ")");
  if ((A.QComplex) && (T.QComplex))
  {                             /* all complex */
    for (k = 0; k < A.N; k++)
    {
      for (l = 0; l < T.N; l++)
      {
        if (T.j[l] == A.j[k])
        {
          for (m = 0; m < T.N; m++)
          {
            if (T.j[m] == A.i[k])
            {
              tempRe =
                T.VRe[m] * A.VRe[k] * T.VRe[l] + T.VIm[m] * A.VRe[k] * T.VIm[l] + T.VIm[m] * A.VIm[k] * T.VRe[l] -
                T.VRe[m] * A.VIm[k] * T.VIm[l];
              tempIm =
                T.VIm[m] * A.VIm[k] * T.VIm[l] - T.VIm[m] * A.VRe[k] * T.VRe[l] + T.VRe[m] * A.VIm[k] * T.VRe[l] +
                T.VRe[m] * A.VRe[k] * T.VIm[l];
              if ((Sqr(tempRe) + Sqr(tempIm)) > epssqr)
              {
                if (ComplexSparseMatrixAddElement(M, T.i[m], T.i[l], tempRe, tempIm))
                {
                  FreeSparseMatrix(M);
                  return 1;
                }
              }
            }
          }
        }
      }
    }
    return 0;
  }
  if (A.QComplex)
  {                             /* A complex, T real */
    for (k = 0; k < A.N; k++)
    {
      for (l = 0; l < T.N; l++)
      {
        if (T.j[l] == A.j[k])
        {
          for (m = 0; m < T.N; m++)
          {
            if (T.j[m] == A.i[k])
            {
              tempRe = T.VRe[m] * A.VRe[k] * T.VRe[l];
              tempIm = T.VRe[m] * A.VIm[k] * T.VRe[l];
              if ((Sqr(tempRe) + Sqr(tempIm)) > epssqr)
              {
                if (ComplexSparseMatrixAddElement(M, T.i[m], T.i[l], tempRe, tempIm))
                {
                  FreeSparseMatrix(M);
                  return 1;
                }
              }
            }
          }
        }
      }
    }
    return 0;
  }
  if (T.QComplex)
  {                             /* T complex, A real */
    for (k = 0; k < A.N; k++)
    {
      for (l = 0; l < T.N; l++)
      {
        if (T.j[l] == A.j[k])
        {
          for (m = 0; m < T.N; m++)
          {
            if (T.j[m] == A.i[k])
            {
              tempRe = T.VRe[m] * A.VRe[k] * T.VRe[l] + T.VIm[m] * A.VRe[k] * T.VIm[l];
              tempIm = T.VRe[m] * A.VRe[k] * T.VIm[l] - T.VIm[m] * A.VRe[k] * T.VRe[l];
              if ((Sqr(tempRe) + Sqr(tempIm)) > epssqr)
              {
                if (ComplexSparseMatrixAddElement(M, T.i[m], T.i[l], tempRe, tempIm))
                {
                  FreeSparseMatrix(M);
                  return 1;
                }
              }
            }
          }
        }
      }
    }
    return 0;
  }
  /* all real */
  for (k = 0; k < A.N; k++)
  {
    for (l = 0; l < T.N; l++)
    {
      if (T.j[l] == A.j[k])
      {
        for (m = 0; m < T.N; m++)
        {
          if (T.j[m] == A.i[k])
          {
            tempRe = T.VRe[m] * A.VRe[k] * T.VRe[l];
            if (fabs(tempRe) > epsilon)
            {
              if (RealSparseMatrixAddElement(M, T.i[m], T.i[l], tempRe))
              {
                FreeSparseMatrix(M);
                return 1;
              }
            }
          }
        }
      }
    }
  }
  return 0;
}

/* M = T^* A T^{transpose}
 * M[i][j]=\sum_k,l T[i][l]^* A[l][k] T[j][k]
 * (A.Ni==A.Nj==T.Nj)
 * M.Ni=M.Nj=T.Ni;
 * M on entrance not initiallized
 */
uint32_t SparseMatrixDenseTransform(SparseMatrixType * const M, SparseMatrixType const A, MatrixType const T, double epsilon)
{
  /**** To Do ****/
  return 0;
}

/* M = T^* A T^{transpose}
 * M[i][j]=\sum_k,l T[i][l]^* A[l][k] T[j][k]
 * (A.Ni==A.Nj==T.Nj)
 * M.Ni=M.Nj=T.Ni;
 * M on entrance not initiallized
 */
uint32_t DenseMatrixSparseTransform(SparseMatrixType * const M, MatrixType const A, SparseMatrixType const T, double epsilon)
{
  /**** To Do ****/
  return 0;
}

uint32_t SparseMatrixTransformPart(SparseMatrixType * const M, SparseMatrixType const A, SparseMatrixType const T, double epsilon, const uint32_t MNiMin, const uint32_t MNiMax, const uint32_t MNjMin, const uint32_t MNjMax)
{
  /**** To Do ****/
  return 0;
}

uint32_t SparseMatrixDenseTransformPart(SparseMatrixType * const M, SparseMatrixType const A, MatrixType const T, double epsilon, const uint32_t MNiMin, const uint32_t MNiMax, const uint32_t MNjMin, const uint32_t MNjMax)
{
  /**** To Do ****/
  return 0;
}

uint32_t DenseMatrixSparseTransformPart(SparseMatrixType * const M, MatrixType const A, SparseMatrixType const T, double epsilon, const uint32_t MNiMin, const uint32_t MNiMax, const uint32_t MNjMin, const uint32_t MNjMax)
{
  /**** To Do ****/
  return 0;
}

/* 
 * M = A B
 * On entry M is not initiallized
 * M[i][j]=\sum_k A[i][k] B[k][j]
 * if ((A.Nj==B.Ni)) then error
 * (M.Ni==A.Ni)&&(M.Nj==B.Nj)
 */
uint32_t SparseMatrixMultiply(SparseMatrixType * const M, SparseMatrixType const A, SparseMatrixType const B)
{
  uint32_t NextIndex, ThereIsANext, index, iA, i, j, k;
  if (A.Nj != B.Ni)
  {
    printf("Invallid dimensions for matrix multiplication A.Nj should be B.Ni\n");
    printf("A.Ni=%i A.Nj=%i B.Ni=%i B.Nj=%i \n", A.Ni, A.Nj, B.Ni, B.Nj);
    return 1;
  }
  M->NBitsKey = (A.NBitsKey + B.NBitsKey) / 2;
  if ((A.QComplex) || (B.QComplex))
  {
    M->QComplex = 1;
  }
  else
  {
    M->QComplex = 0;
  }
  M->Ni = A.Ni;
  M->Nj = B.Nj;
  if (InitSparseMatrix(M))
    return 1;
  if ((A.QComplex) && (B.QComplex))
  {                             /* all complex */
    for (iA = 0; iA < A.N; iA++)
    {                           /* i=A.i[iA]; k=A.j[iA]; */
      i = A.i[iA];
      k = A.j[iA];
      for (j = 0; j < B.Nj; j++)
      {                         /* if element B[k][j] excist add value A[i][k]*B[k][j] to M[i][j] otherwise do nothing */
        index = ((k) ^ ((j >> 16) | (j << 16)));
        if ((B.KeyList[hash(&index, 4, B.NBitsKey)]) != 0)
        {                       /* the key is already taken, test if this is due to a collision or if i is already in the opperator */
          ThereIsANext = 1;
          NextIndex = B.KeyList[hash(&index, 4, B.NBitsKey)] - 1;
          while (ThereIsANext)
          {
            if ((k == B.i[NextIndex]) && (j == B.j[NextIndex]))
            {                   /* B[k][j] in the opperator, add A[i][k]*B[k][j] to M[i][j] and leave loop */
              ThereIsANext = 0;
              if (ComplexSparseMatrixAddElement
                  (M, i, j, A.VRe[iA] * B.VRe[NextIndex] - A.VIm[iA] * B.VIm[NextIndex],
                   A.VRe[iA] * B.VIm[NextIndex] + A.VIm[iA] * B.VRe[NextIndex]))
              {
                FreeSparseMatrix(M);
                return 1;
              }
            }
            else
            {                   /* i,j not equal to the i,j found which has the same index in the keylist, if there is a next continue, otherwise add the element */
              if (M->Next[NextIndex])
              {                 /* The key is in use and there are more opperators with this key */
                NextIndex = M->Next[NextIndex] - 1;
              }
            }
          }
        }
      }
    }
    return 0;
  }
  if (A.QComplex)
  {                             /* A complex, B real */
    for (iA = 0; iA < A.N; iA++)
    {                           /* i=A.i[iA]; k=A.j[iA]; */
      i = A.i[iA];
      k = A.j[iA];
      for (j = 0; j < B.Nj; j++)
      {                         /* if element B[k][j] excist add value A[i][k]*B[k][j] to M[i][j] otherwise do nothing */
        index = ((k) ^ ((j >> 16) | (j << 16)));
        if ((B.KeyList[hash(&index, 4, B.NBitsKey)]) != 0)
        {                       /* the key is already taken, test if this is due to a collision or if i is already in the opperator */
          ThereIsANext = 1;
          NextIndex = B.KeyList[hash(&index, 4, B.NBitsKey)] - 1;
          while (ThereIsANext)
          {
            if ((k == B.i[NextIndex]) && (j == B.j[NextIndex]))
            {                   /* B[k][j] in the opperator, add A[i][k]*B[k][j] to M[i][j] and leave loop */
              ThereIsANext = 0;
              if (ComplexSparseMatrixAddElement(M, i, j, A.VRe[iA] * B.VRe[NextIndex], A.VIm[iA] * B.VRe[NextIndex]))
              {
                FreeSparseMatrix(M);
                return 1;
              }
            }
            else
            {                   /* i,j not equal to the i,j found which has the same index in the keylist, if there is a next continue, otherwise add the element */
              if (M->Next[NextIndex])
              {                 /* The key is in use and there are more opperators with this key */
                NextIndex = M->Next[NextIndex] - 1;
              }
            }
          }
        }
      }
    }
    return 0;
  }
  if (B.QComplex)
  {                             /* B complex, A real */
    for (iA = 0; iA < A.N; iA++)
    {                           /* i=A.i[iA]; k=A.j[iA]; */
      i = A.i[iA];
      k = A.j[iA];
      for (j = 0; j < B.Nj; j++)
      {                         /* if element B[k][j] excist add value A[i][k]*B[k][j] to M[i][j] otherwise do nothing */
        index = ((k) ^ ((j >> 16) | (j << 16)));
        if ((B.KeyList[hash(&index, 4, B.NBitsKey)]) != 0)
        {                       /* the key is already taken, test if this is due to a collision or if i is already in the opperator */
          ThereIsANext = 1;
          NextIndex = B.KeyList[hash(&index, 4, B.NBitsKey)] - 1;
          while (ThereIsANext)
          {
            if ((k == B.i[NextIndex]) && (j == B.j[NextIndex]))
            {                   /* B[k][j] in the opperator, add A[i][k]*B[k][j] to M[i][j] and leave loop */
              ThereIsANext = 0;
              if (ComplexSparseMatrixAddElement(M, i, j, A.VRe[iA] * B.VRe[NextIndex], A.VRe[iA] * B.VIm[NextIndex]))
              {
                FreeSparseMatrix(M);
                return 1;
              }
            }
            else
            {                   /* i,j not equal to the i,j found which has the same index in the keylist, if there is a next continue, otherwise add the element */
              if (M->Next[NextIndex])
              {                 /* The key is in use and there are more opperators with this key */
                NextIndex = M->Next[NextIndex] - 1;
              }
            }
          }
        }
      }
    }
    return 0;
  }
  /* all real */
  for (iA = 0; iA < A.N; iA++)
  {                             /* i=A.i[iA]; k=A.j[iA]; */
    i = A.i[iA];
    k = A.j[iA];
    for (j = 0; j < B.Nj; j++)
    {                           /* if element B[k][j] excist add value A[i][k]*B[k][j] to M[i][j] otherwise do nothing */
      index = ((k) ^ ((j >> 16) | (j << 16)));
      if ((B.KeyList[hash(&index, 4, B.NBitsKey)]) != 0)
      {                         /* the key is already taken, test if this is due to a collision or if i is already in the opperator */
        ThereIsANext = 1;
        NextIndex = B.KeyList[hash(&index, 4, B.NBitsKey)] - 1;
        while (ThereIsANext)
        {
          if ((k == B.i[NextIndex]) && (j == B.j[NextIndex]))
          {                     /* B[k][j] in the opperator, add A[i][k]*B[k][j] to M[i][j] and leave loop */
            ThereIsANext = 0;
            if (RealSparseMatrixAddElement(M, i, j, A.VRe[iA] * B.VRe[NextIndex]))
            {
              FreeSparseMatrix(M);
              return 1;
            }
          }
          else
          {                     /* i,j not equal to the i,j found which has the same index in the keylist, if there is a next continue, otherwise add the element */
            if (M->Next[NextIndex])
            {                   /* The key is in use and there are more opperators with this key */
              NextIndex = M->Next[NextIndex] - 1;
            }
          }
        }
      }
    }
  }
  return 0;
}


/* 
 * M = A B
 * On entry M is not initiallized
 * M[i][j]=\sum_k A[i][k] B[k][j]
 * if ((A.Nj==B.Ni)) then error
 * (M.Ni==A.Ni)&&(M.Nj==B.Nj)
 */
uint32_t SparseMatrixDenseMultiply(SparseMatrixType * const M, SparseMatrixType const A, MatrixType const B)
{
  /**** To Do ****/
  return 0;
}

/* 
 * M = A B
 * On entry M is not initiallized
 * M[i][j]=\sum_k A[i][k] B[k][j]
 * if ((A.Nj==B.Ni)) then error
 * (M.Ni==A.Ni)&&(M.Nj==B.Nj)
 */
uint32_t DenseMatrixSparseMultiply(SparseMatrixType * const M, MatrixType const A, SparseMatrixType const B)
{
  /**** To Do ****/
  return 0;
}

/* Create a new matrix that has dimensions (1+iEnd-iStart,1+jEnd-jStart) 
 * and is a copy of a part of A
 * M not allocated on call of the subroutine
 */
uint32_t SparseMatrixPart(SparseMatrixType * const M, SparseMatrixType const A, const uint32_t iStart, const uint32_t iEnd, const uint32_t jStart, const uint32_t jEnd)
{
  uint32_t l;
  if ((iEnd < iStart) || (jEnd < jStart) || (iEnd >= A.Ni) || (jEnd >= A.Nj))
    return 1;
  M->Ni = 1 + iEnd - iStart;
  M->Nj = 1 + jEnd - jStart;
  M->QComplex = A.QComplex;
  if (InitSparseMatrix(M))
    return 1;
  strcpy(M->Name, "Part of ");
  strncat(M->Name, A.Name, Max(0, MAXStrLength - 10));
  if (M->QComplex)
  {                             /* Complex */
    for (l = 0; l < A.N; l++)
    {
      if ((A.i[l] >= iStart) && (A.i[l] <= iEnd) && (A.j[l] >= jStart) && (A.j[l] <= jEnd))
      {
        if (ComplexSparseMatrixAddElement(M, A.i[l] - iStart, A.j[l] - jStart, A.VRe[l], A.VIm[l]))
        {
          FreeSparseMatrix(M);
          return 1;
        }
      }
    }
  }
  else
  {                             /* Real */
    for (l = 0; l < A.N; l++)
    {
      if ((A.i[l] >= iStart) && (A.i[l] <= iEnd) && (A.j[l] >= jStart) && (A.j[l] <= jEnd))
      {
        if (RealSparseMatrixAddElement(M, A.i[l] - iStart, A.j[l] - jStart, A.VRe[l]))
        {
          FreeSparseMatrix(M);
          return 1;
        }
      }
    }
  }
  return 0;
}


/*
 * M<-A
 */
uint32_t SparseMatrixCopy(SparseMatrixType * const M, SparseMatrixType const A)
{
  uint32_t i;
  M->Ni = A.Ni;
  M->Nj = A.Nj;
  M->N = A.N;
  M->NAloc = M->N;
  M->QComplex = A.QComplex;
  M->NBitsKey = A.NBitsKey;
  M->Next = (uint32_t*)malloc(M->NAloc * sizeof(uint32_t));
  M->VRe = (double*)malloc(M->NAloc * sizeof(double));
  if (M->QComplex)
    M->VIm = (double*)malloc(M->NAloc * sizeof(double));
  M->i = (uint32_t*)malloc(M->NAloc * sizeof(uint32_t));
  M->j = (uint32_t*)malloc(M->NAloc * sizeof(uint32_t));
  M->KeyList = (uint32_t*)calloc((uint32_t) 1 << (M->NBitsKey), sizeof(uint32_t));
  if (M->KeyList == NULL)
  {
    if (M->j != NULL)
      free(M->j);
    if (M->i != NULL)
      free(M->i);
    if (M->VIm != NULL)
      free(M->VIm);
    if (M->VRe != NULL)
      free(M->VRe);
    if (M->Next != NULL)
      free(M->Next);
    return 1;
  }
  if (M->QComplex)
  {
    for (i = 0; i < M->N; i++)
    {
      M->Next[i] = A.Next[i];
      M->VRe[i] = A.VRe[i];
      M->VIm[i] = A.VIm[i];
      M->i[i] = A.i[i];
      M->j[i] = A.j[i];
    }
  }
  else
  {
    for (i = 0; i < M->N; i++)
    {
      M->Next[i] = A.Next[i];
      M->VRe[i] = A.VRe[i];
      M->i[i] = A.i[i];
      M->j[i] = A.j[i];
    }
  }
  for (i = 0; i < ((uint32_t) 1 << (M->NBitsKey)); i++)
  {
    M->KeyList[i] = A.KeyList[i];
  }
  return 0;
}

/* b <- M a
 * b already allocated
 * b[i]=sum_j M[i][j] a[j]
 */
uint32_t RealSparseMatrixVector(SparseMatrixType const M, const double *const a, double *const b)
{
  int    k;
  if (M.QComplex)
  {
    printf("tried to use the function RealMatrixVector on a complex Matrix\n");
    printf("Name of the matrix %s\n", M.Name);
    return 1;
  }
  for (k = 0; k < M.Ni; k++)
  {
    b[k] = 0;
  }
  for (k = 0; k < M.N; k++)
  {
    b[M.i[k]] += M.VRe[k] * a[M.j[k]];
  }
  return 0;
}

/* b <- M^T a
 * b already allocated
 * b[i]=sum_j M[j][i] a[j]
 */
uint32_t RealTransposeSparseMatrixVector(SparseMatrixType const M, const double *const a, double *const b)
{
  int    k;
  if (M.QComplex)
  {
    printf("tried to use the function RealMatrixVector on a complex Matrix\n");
    printf("Name of the matrix %s\n", M.Name);
    return 1;
  }
  for (k = 0; k < M.Nj; k++)
  {
    b[k] = 0;
  }
  for (k = 0; k < M.N; k++)
  {
    b[M.j[k]] += M.VRe[k] * a[M.i[k]];
  }
  return 0;
}

/* b <- M a
 * b already allocated
 * b[i]=sum_j M[i][j] a[j]
 */
void ComplexSparseMatrixVector(SparseMatrixType const M, const double *const aRe, const double *const aIm, double *const bRe, double *const bIm)
{
  int    k;
  if (M.QComplex)
  {                             /* M Complex */
    for (k = 0; k < M.Ni; k++)
    {
      bRe[k] = 0;
      bIm[k] = 0;
    }
    for (k = 0; k < M.N; k++)
    {
      bRe[M.i[k]] += M.VRe[k] * aRe[M.j[k]] - M.VIm[k] * aIm[M.j[k]];
      bIm[M.i[k]] += M.VRe[k] * aIm[M.j[k]] + M.VIm[k] * aRe[M.j[k]];
    }
  }
  else
  {                             /* M Real */
    for (k = 0; k < M.Ni; k++)
    {
      bRe[k] = 0;
      bIm[k] = 0;
    }
    for (k = 0; k < M.N; k++)
    {
      bRe[M.i[k]] += M.VRe[k] * aRe[M.j[k]];
      bIm[M.i[k]] += M.VRe[k] * aIm[M.j[k]];
    }
  }
}

/* b <- M^+ a
 * b already allocated
 * b[i]=sum_j M[j][i]^* a[j]
 */
void ComplexConjugateTransposeSparseMatrixVector(SparseMatrixType const M, const double *const aRe, const double *const aIm, double *const bRe, double *const bIm)
{
  int    k;
  if (M.QComplex)
  {                             /* M Complex */
    for (k = 0; k < M.Nj; k++)
    {
      bRe[k] = 0;
      bIm[k] = 0;
    }
    for (k = 0; k < M.N; k++)
    {
      bRe[M.j[k]] += M.VRe[k] * aRe[M.i[k]] + M.VIm[k] * aIm[M.i[k]];
      bIm[M.j[k]] += M.VRe[k] * aIm[M.i[k]] - M.VIm[k] * aRe[M.i[k]];
    }
  }
  else
  {                             /* M Real */
    for (k = 0; k < M.Nj; k++)
    {
      bRe[k] = 0;
      bIm[k] = 0;
    }
    for (k = 0; k < M.N; k++)
    {
      bRe[M.j[k]] += M.VRe[k] * aRe[M.i[k]];
      bIm[M.j[k]] += M.VRe[k] * aIm[M.i[k]];
    }
  }
}


/* b <- a M a
 * b already allocated
 * b=sum_ij a[i] M[i][j] a[j]
 */
uint32_t RealVectorSparseMatrixVector(SparseMatrixType const M, const double *const a, double *const b)
{
  volatile double c, y, t;
  int    k;
  if (M.QComplex)
  {
    printf("tried to use the function RealMatrixVector on a complex Matrix\n");
    printf("Name of the matrix %s\n", M.Name);
    return 1;
  }
  *b = 0;
  c = 0.0;
  y = 0.0;
  t = 0.0;
  for (k = 0; k < M.N; k++)
  {
    y = a[M.i[k]] * M.VRe[k] * a[M.j[k]] - c;
    t = *b + y;
    c = t - *b;
    c = c - y;
    *b = t;
    /* *b+=a[M.i[k]]*M.VRe[k]*a[M.j[k]]; */
  }
  return 0;
}


/* b <- a^* M a
 * b already allocated
 * b=sum_ij a[i]^* M[i][j] a[j]
 */
void ComplexVectorSparseMatrixVector(SparseMatrixType const M, const double *const aRe, const double *const aIm, double *const bRe, double *const bIm)
{
  int    k;
  *bRe = 0;
  *bIm = 0;
  if (M.QComplex)
  {                             /* M Complex */
    for (k = 0; k < M.N; k++)
    {
      *bRe +=
        aRe[M.i[k]] * M.VRe[k] * aRe[M.j[k]] + aIm[M.i[k]] * M.VIm[k] * aRe[M.j[k]] +
        aIm[M.i[k]] * M.VRe[k] * aIm[M.j[k]] - aRe[M.i[k]] * M.VIm[k] * aIm[M.j[k]];
      *bIm +=
        aIm[M.i[k]] * M.VIm[k] * aIm[M.j[k]] - aIm[M.i[k]] * M.VRe[k] * aRe[M.j[k]] +
        aRe[M.i[k]] * M.VRe[k] * aIm[M.j[k]] + aRe[M.i[k]] * M.VIm[k] * aRe[M.j[k]];
    }
  }
  else
  {                             /* M Real */
    for (k = 0; k < M.N; k++)
    {
      *bRe += aRe[M.i[k]] * M.VRe[k] * aRe[M.j[k]] + aIm[M.i[k]] * M.VRe[k] * aIm[M.j[k]];
      *bIm += aRe[M.i[k]] * M.VRe[k] * aIm[M.j[k]] - aIm[M.i[k]] * M.VRe[k] * aRe[M.j[k]];
    }
  }
}

/* M transformed such that
 * M[i]^* . M[j] = delta_{i,j}
 * fully dependent vectors are set to zero
 * At the same time the overlap between the new and old M is kept "as large as possible"
 * Whereby M[0] is first optimized, then M[1], etc.
 * if M[1]=1/2 (M[0]+M[2]) then M[1] becomes M[2]
 * and M[2] is set to zero
 */
void SparseMatrixOrthonormalize(MatrixType * const M, double epsilon)
{
  /**** To Do ****/
}
