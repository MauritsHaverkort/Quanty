uint32_t InitListOfPoles(ListOfPoles * M)
{
  strcpy(M->Name, "List of Poles");
  if (M->N == 0)
  {
    printf("Initialized a ListOfPoles with size zero.\n");
    fflush(stdout);
    return 2;
  };
  M->a = (double*)malloc((M->N) * sizeof(double));
  if (M->a == NULL)
  {
    printf("malloc for M->a failed in ListOfPoles\n");
    fflush(stdout);
    return 1;
  };
  M->b = (double*)malloc((M->N) * sizeof(double));
  if (M->b == NULL)
  {
    free(M->a);
    printf("malloc for M->b failed in ListOfPoles\n");
    fflush(stdout);
    return 1;
  }
  M->a0 = 0;
  M->mu = 0;
  return 0;
}

void FreeListOfPoles(ListOfPoles * M)
{
  if (M->N != 0)
  {
    free(M->a);
    free(M->b);
  }
  M->a0 = 0;
}

void ScaleListOfPoles(ListOfPoles * const A, double const scale)
{
  uint32_t i;
  for (i = 0; i < A->N; i++)
  {
    A->b[i] *= scale;
  }
}

double CalcNormListOfPoles(ListOfPoles * const A)
{
  double norm=0;
  for(int i=0;i<A->N;i++)
  {
    norm+=A->b[i];
  }
  return norm;
}

void NormalizeListOfPoles(ListOfPoles * const A)
{
  double norm=0;
  norm = CalcNormListOfPoles(A);
  if(norm!=0) norm=1/norm;
  ScaleListOfPoles(A, norm);
}

void SetNegativeToZeroListOfPoles(ListOfPoles * const A)
{
  for(int i=0;i<A->N;i++)
  {
    if(A->b[i]<0)A->b[i]=0;
  }
}

// This function can leave negative spectral weight at the edges
// otherwise it removes the negative spectral weight weilst keeping locally
// the frist moments correct.
void MergeNegativeToNeighbourListOfPolesSorted(ListOfPoles * const A)
{
  int i,j;
  int muindex;
  double WeightLeft,WeightRight,dE;
  muindex=0;
  while((muindex<A->N)&&(A->b[muindex]<A->mu))
  {
    muindex++;
  }
  muindex--;
  // left of chemical potential
  for(i=muindex-1;i>0;i--)
  {
    if(A->b[i]>=0) continue;
    for(j=i+1;j<=muindex;j++)
    {
      if(A->b[j]==0) continue;
      dE = A->a[j] - A->a[i-1];
      WeightRight = ( A->a[i] - A->a[i-1] ) / dE;
      WeightLeft = ( A->a[j] - A->a[i] ) / dE;
      if(-(WeightRight * A->b[i]) > A->b[j] )
      {
        // WeigthToShiftRigth = A->b[j] / WeightRight;
        A->b[i]   += (WeightLeft+1) * A->b[j] / WeightRight;
        A->b[i-1] -= WeightLeft  * A->b[j] / WeightRight;
        A->b[j] = 0;
        continue;
      }
      else
      {
        A->b[i-1] += WeightLeft  * A->b[i];
        A->b[j]   += WeightRight * A->b[i];
        A->b[i] = 0;
        break;
      }
    }
  }
  // right of chemical potential
  muindex++;
  if(muindex==A->N) return;
  if(A->mu == A->b[muindex]) muindex++;
  if(muindex==A->N) return;
  for(i=muindex+1;i<A->N-1;i++)
  {
    if(A->b[i]>=0) continue;
    for(j=i-1;j>=muindex;j--)
    {
      if(A->b[j]==0) continue;
      dE = A->a[i-1] - A->a[j];
      WeightLeft  = ( A->a[i] - A->a[i-1] ) / dE;
      WeightRight = ( A->a[j] - A->a[i] ) / dE;
      if(-(WeightLeft * A->b[i]) > A->b[j] )
      {
        // WeigthToShiftLeft = A->b[j] / WeightRight;
        A->b[i]   += (WeightRight+1) * A->b[j] / WeightLeft;
        A->b[i+1] -= WeightRight  * A->b[j] / WeightLeft;
        A->b[j] = 0;
        continue;
      }
      else
      {
        A->b[i+1] += WeightRight * A->b[i];
        A->b[j]   += WeightLeft  * A->b[i];
        A->b[i] = 0;
        break;
      }
    }
  }
}

uint32_t AddListOfPoles(ListOfPoles * const A, ListOfPoles B)
{
  double *temp;
  uint32_t i;
  temp = (double*)realloc(A->a, (A->N + B.N) * sizeof(double));
  if (temp == NULL)
  {
    printf("realloc of A->a failed in AddListOfPoles\n");
    fflush(stdout);
    return 1;
  }
  A->a = temp;
  temp = (double*)realloc(A->b, (A->N + B.N) * sizeof(double));
  if (temp == NULL)
  {
    printf("realloc of A->b failed in AddListOfPoles\n");
    fflush(stdout);
    return 1;
  }
  A->b = temp;
  for (i = 0; i < B.N; i++)
  {
    A->a[i + A->N] = B.a[i];
    A->b[i + A->N] = B.b[i];
  }
  A->a0 += B.a0;
  A->N += B.N;
  return 0;
}


/*   A <- alpha A + (1-alpha) B;
 *
 */
uint32_t MixListOfPoles(ListOfPoles * const A, ListOfPoles const B, double const alpha)
{
  double *temp;
  uint32_t i;
  temp = (double*)realloc(A->a, (A->N + B.N) * sizeof(double));
  if (temp == NULL)
  {
    printf("realloc of A->a failed in AddListOfPoles\n");
    fflush(stdout);
    return 1;
  }
  A->a = temp;
  temp = (double*)realloc(A->b, (A->N + B.N) * sizeof(double));
  if (temp == NULL)
  {
    printf("realloc of A->b failed in AddListOfPoles\n");
    fflush(stdout);
    return 1;
  }
  A->b = temp;
  for (i = 0; i < A->N; i++)
  {
    A->b[i] *= alpha;
  }
  for (i = 0; i < B.N; i++)
  {
    A->a[i + A->N] = B.a[i];
    A->b[i + A->N] = B.b[i] * (1 - alpha);
  }
  A->a0 = alpha * A->a0 + (1 - alpha) * B.a0;
  A->mu = alpha * A->mu + (1 - alpha) * B.mu;
  A->N += B.N;
  return 0;
}

uint32_t SymmetrizeListOfPoles(ListOfPoles * const A)
{
  double *temp;
  uint32_t i;
  temp = (double*)realloc(A->a, (2 * A->N) * sizeof(double));
  if (temp == NULL)
  {
    printf("realloc of A->a failed in AddListOfPoles\n");
    fflush(stdout);
    return 1;
  }
  A->a = temp;
  temp = (double*)realloc(A->b, (2 * A->N) * sizeof(double));
  if (temp == NULL)
  {
    printf("realloc of A->b failed in AddListOfPoles\n");
    fflush(stdout);
    return 1;
  }
  A->b = temp;
  for (i = 0; i < A->N; i++)
  {
    A->a[i + A->N] = 2 * A->mu - A->a[i];
    A->b[i + A->N] = A->b[i];
  }
  A->N *= 2;
  ScaleListOfPoles(A, 0.5);
  return 0;
}

int OrderListOfPolesCompare(double *elem1, double *elem2)
{
  if (elem1[0] < elem2[0])
    return -1;
  else if (elem1[0] > elem2[0])
    return 1;
  return 0;
}

void OrderListOfPoles(ListOfPoles * const M)
{
  double *Mab = (double*)malloc(2 * M->N * sizeof(double));
  if (Mab == NULL)
  {
    printf("malloc for Mab failed in ListOfPoles\n");
    fflush(stdout);
    exit(1);
  }  
  uint32_t i;
  for (i = 0; i < M->N; i++)
  {
    Mab[2 * i] = M->a[i];
    Mab[2 * i + 1] = M->b[i];
  }
  qsort((void *)Mab, M->N, 2 * sizeof(double), (compfn) OrderListOfPolesCompare);
  for (i = 0; i < M->N; i++)
  {
    M->a[i] = Mab[2 * i];
    M->b[i] = Mab[2 * i + 1];
  }
  free(Mab);
}

void OrderListOfPolesOld(ListOfPoles * const M)
{
  uint32_t i, j;
  double temp;
  for (i = 1; i < M->N; i++)
  {
    for (j = i; j > 0; j--)
    {
      if (M->a[j] < M->a[j - 1])
      {
        temp = M->a[j];
        M->a[j] = M->a[j - 1];
        M->a[j - 1] = temp;
        temp = M->b[j];
        M->b[j] = M->b[j - 1];
        M->b[j - 1] = temp;
      }
      else
      {
        break;
      }
    }
  }
}


/*  reduces a list of length Ntot to a list of length Ntot-1 by setting the n-th element to zero
 *  keeps the first Ntot-1 moments conserved, does not shift the poles
 *  
 *  Input
 *    a energy of the poles
 *    b weight of the poles
 *
 *  Output
 *    b weight of the poles
 */
void MergeOnePoleWithNothers(double const *a, double *const b, uint32_t Ntot, uint32_t n)
{

}



/*  Create a new set of poles, wich all have a minimal weight of bmin
 *
 */
void ReduceListOfPolesWeigth(ListOfPoles * M, double const bmin, double const Ef)
{
  uint32_t i, j;
  OrderListOfPoles(M);
  if ((M->N == 1) && (M->b[0] < bmin))
  {
    M->b[0] = 0;
    M->a[0] = 0;
    M->N = 0;
  }
  else
  {
    for (i = 0; i < M->N; i++)
    {
      if (M->b[i] < bmin)
      {
        if (i == 0)
        {
          if ((M->a[0] < Ef) && (M->a[1] > Ef))
          {
            printf("only one pole with weight %22.15E below Ef, removed pole and shifted weight true Ef\n", M->a[0]);
          }
          M->a[1] = (M->b[0] * M->a[0] + M->b[1] * M->a[1]) / (M->b[0] + M->b[1] + DBL_MIN);
          M->b[1] += M->b[0];
        }
        else if (i == M->N - 1)
        {
          if ((M->a[i - 1] < Ef) && (M->a[i] > Ef))
          {
            printf("only one pole with weight %22.15E above Ef, removed pole and shifted weight true Ef\n", M->a[i]);
          }
          M->a[i - 1] = (M->b[i - 1] * M->a[i - 1] + M->b[i] * M->a[i]) / (M->b[i] + M->b[i - 1] + DBL_MIN);
          M->b[i - 1] += M->b[i];
        }
        else
        {
          if ((M->a[i - 1] < Ef) && (M->a[i] > Ef))
          {
            M->a[i + 1] = (M->b[i + 1] * M->a[i + 1] + M->b[i] * M->a[i]) / (M->b[i] + M->b[i + 1] + DBL_MIN);
            M->b[i + 1] += M->b[i];
          }
          else if ((M->a[i] < Ef) && (M->a[i + 1] > Ef))
          {
            M->a[i - 1] = (M->b[i - 1] * M->a[i - 1] + M->b[i] * M->a[i]) / (M->b[i] + M->b[i - 1] + DBL_MIN);
            M->b[i - 1] += M->b[i];
          }
          else
          {
            M->b[i - 1] += M->b[i] * (DBL_MIN + M->a[i + 1] - M->a[i]) / (M->a[i + 1] - M->a[i - 1] + 2 * DBL_MIN);
            M->b[i + 1] += M->b[i] * (DBL_MIN + M->a[i] - M->a[i - 1]) / (M->a[i + 1] - M->a[i - 1] + 2 * DBL_MIN);
          }
        }
        for (j = i; j < M->N - 1; j++)
        {
          M->a[j] = M->a[j + 1];
          M->b[j] = M->b[j + 1];
        }
        M->N--;
        if (i != 0)
          i--;
      }
    }
  }
}

/*  Create a new set of poles, which has lenght N, by removing the M->N-N poles with the smallest weight
 *  This algorithm combines the pole with the smallest weight with its nearest two neighbors
 *  and repeats this untill the number of nonzero poles is N
 *
 *  There is no spectral weight transfered through the Fermi energy.
 */
void ReduceListOfPolesNumberFermi(ListOfPoles * M, uint32_t N, double Ef)
{
  if (M->N < N)
    return;
  uint32_t i, imin;
  double biminLeft, biminRight, temp;
  OrderListOfPoles(M);
  while (M->N > N)
  {
    imin = 0;
    for (i = 1; i < M->N; i++)
    {
      if (M->b[i] < M->b[imin])
      {
        imin = i;
      }
    }
    if (M->b[imin] > DBL_MIN)
    {
      if (imin == 0)
      {
        M->a[imin + 1] =
          ((M->b[imin] / M->b[imin + 1]) * M->a[imin] + M->a[imin + 1]) / ((M->b[imin] / M->b[imin + 1]) + 1);
        M->b[imin + 1] = M->b[imin] + M->b[imin + 1];
      }
      else if (imin == M->N - 1)
      {
        M->a[imin - 1] =
          ((M->b[imin] / M->b[imin - 1]) * M->a[imin] + M->a[imin - 1]) / ((M->b[imin] / M->b[imin - 1]) + 1);
        M->b[imin - 1] = M->b[imin] + M->b[imin - 1];
      }
      else if ((M->a[imin] > Ef) && (M->a[imin - 1] < Ef))
      {
        M->a[imin + 1] =
          ((M->b[imin] / M->b[imin + 1]) * M->a[imin] + M->a[imin + 1]) / ((M->b[imin] / M->b[imin + 1]) + 1);
        M->b[imin + 1] = M->b[imin] + M->b[imin + 1];
      }
      else if ((M->a[imin] < Ef) && (M->a[imin + 1] > Ef))
      {
        M->a[imin - 1] =
          ((M->b[imin] / M->b[imin - 1]) * M->a[imin] + M->a[imin - 1]) / ((M->b[imin] / M->b[imin - 1]) + 1);
        M->b[imin - 1] = M->b[imin] + M->b[imin - 1];
      }
      else
      {
        if (fabs(M->a[imin + 1] - M->a[imin - 1]) < fabs(10 * (DBL_MIN + DBL_EPSILON * M->a[imin])) ||
            (fabs(M->a[imin + 1] - M->a[imin]) < fabs(10 * (DBL_MIN + DBL_EPSILON * M->a[imin])) &&
             fabs(M->a[imin] - M->a[imin-1]) < fabs(10 * (DBL_MIN + DBL_EPSILON * M->a[imin]))))
        {
          biminLeft = M->b[imin] * (0.5);
          biminRight = M->b[imin] * (0.5);
        }
        else if (fabs(M->a[imin + 1] - M->a[imin]) < fabs(10 * (DBL_MIN + DBL_EPSILON * M->a[imin])))
        {
          biminLeft = 0;
          biminRight = M->b[imin];
        }
        else if (fabs(M->a[imin] - M->a[imin - 1]) < fabs(10 * (DBL_MIN + DBL_EPSILON * M->a[imin])))
        {
          biminLeft = M->b[imin];
          biminRight = 0;
        }
        else
        {
          biminLeft = M->b[imin] * ((M->a[imin + 1] - M->a[imin]) / (M->a[imin + 1] - M->a[imin - 1]));
          biminRight = M->b[imin] * ((M->a[imin] - M->a[imin - 1]) / (M->a[imin + 1] - M->a[imin - 1]));
        }
        M->a[imin + 1] =
          ((biminRight / M->b[imin + 1]) * M->a[imin] + M->a[imin + 1]) / ((biminRight / M->b[imin + 1]) + 1);
        M->b[imin + 1] = biminRight + M->b[imin + 1];
        M->a[imin - 1] =
          ((biminLeft / M->b[imin - 1]) * M->a[imin] + M->a[imin - 1]) / ((biminLeft / M->b[imin - 1]) + 1);
        M->b[imin - 1] = biminLeft + M->b[imin - 1];
      }
    }
    for (i = imin; i < M->N - 1; i++)
    {
      M->a[i] = M->a[i + 1];
      M->b[i] = M->b[i + 1];
    }
    M->N--;
    i = imin;
    while ((i > 1) && (M->a[i - 1] < M->a[i - 2]))
    {
      temp = M->a[i - 1];
      M->a[i - 1] = M->a[i - 2];
      M->a[i - 2] = temp;
      temp = M->b[i - 1];
      M->b[i - 1] = M->b[i - 2];
      M->b[i - 2] = temp;
      i--;
    }
    i = imin;
    while ((i < M->N - 1) && (M->a[i] > M->a[i + 1]))
    {
      temp = M->a[i];
      M->a[i] = M->a[i + 1];
      M->a[i + 1] = temp;
      temp = M->b[i];
      M->b[i] = M->b[i + 1];
      M->b[i + 1] = temp;
      i++;
    }
  }
}


/*  Create a new set of poles, which has lenght N, by removing the M->N-N poles with the smallest weight
 *  This algorithm combines the pole with the smallest weight with its nearest two neighbors
 *  and repeats this untill the number of nonzero poles is N
 */
void ReduceListOfPolesNumber(ListOfPoles * M, uint32_t N)
{
  if (M->N < N)
    return;
  uint32_t i, imin;
  double biminLeft, biminRight, temp;
  OrderListOfPoles(M);
  while (M->N > N)
  {
    imin = 0;
    for (i = 1; i < M->N; i++)
    {
      if (M->b[i] < M->b[imin])
      {
        imin = i;
      }
    }
    if (M->b[imin] > DBL_MIN)
    {
      if (imin == 0)
      {
        M->a[imin + 1] =
          ((M->b[imin] / M->b[imin + 1]) * M->a[imin] + M->a[imin + 1]) / ((M->b[imin] / M->b[imin + 1]) + 1);
        M->b[imin + 1] = M->b[imin] + M->b[imin + 1];
      }
      else if (imin == M->N - 1)
      {
        M->a[imin - 1] =
          ((M->b[imin] / M->b[imin - 1]) * M->a[imin] + M->a[imin - 1]) / ((M->b[imin] / M->b[imin - 1]) + 1);
        M->b[imin - 1] = M->b[imin] + M->b[imin - 1];
      }

      else
      {
        if (fabs((M->a[imin + 1] - M->a[imin - 1])) < fabs(10 * (DBL_MIN + DBL_EPSILON * M->a[imin])))
        {
          biminLeft = M->b[imin] * (0.5);
          biminRight = M->b[imin] * (0.5);
        }
        else if (fabs((M->a[imin + 1] - M->a[imin])) < fabs(10 * (DBL_MIN + DBL_EPSILON * M->a[imin])))
        {
          biminLeft = M->b[imin];
          biminRight = 0;
        }
        else if (fabs((M->a[imin] - M->a[imin - 1])) < fabs(10 * (DBL_MIN + DBL_EPSILON * M->a[imin])))
        {
          biminLeft = 0;
          biminRight = M->b[imin];
        }
        else
        {
          biminLeft = M->b[imin] * ((M->a[imin + 1] - M->a[imin]) / (M->a[imin + 1] - M->a[imin - 1]));
          biminRight = M->b[imin] * ((M->a[imin] - M->a[imin - 1]) / (M->a[imin + 1] - M->a[imin - 1]));
        }
        M->a[imin + 1] =
          ((biminRight / M->b[imin + 1]) * M->a[imin] + M->a[imin + 1]) / ((biminRight / M->b[imin + 1]) + 1);
        M->b[imin + 1] = biminRight + M->b[imin + 1];
        M->a[imin - 1] =
          ((biminLeft / M->b[imin - 1]) * M->a[imin] + M->a[imin - 1]) / ((biminLeft / M->b[imin - 1]) + 1);
        M->b[imin - 1] = biminLeft + M->b[imin - 1];
      }
    }
    for (i = imin; i < M->N - 1; i++)
    {
      M->a[i] = M->a[i + 1];
      M->b[i] = M->b[i + 1];
    }
    M->N--;
    i = imin;
    while ((i > 1) && (M->a[i - 1] < M->a[i - 2]))
    {
      temp = M->a[i - 1];
      M->a[i - 1] = M->a[i - 2];
      M->a[i - 2] = temp;
      temp = M->b[i - 1];
      M->b[i - 1] = M->b[i - 2];
      M->b[i - 2] = temp;
      i--;
    }
    i = imin;
    while ((i < M->N - 1) && (M->a[i] > M->a[i + 1]))
    {
      temp = M->a[i];
      M->a[i] = M->a[i + 1];
      M->a[i + 1] = temp;
      temp = M->b[i];
      M->b[i] = M->b[i + 1];
      M->b[i + 1] = temp;
      i++;
    }
  }
}


/*  Create a new set of poles, which has lenght N, by removing the M->N-N poles such that for each reduction weight * (delta E neighbor) is minimum
 *  This algorithm combines the pole with the smallest moment to one neighbor with its nearest two neighbors
 *  and repeats this untill the number of nonzero poles is N
 */
void ReduceListOfPolesMoment(ListOfPoles * M, uint32_t N)
{
  uint32_t i, imin;
  double biminLeft, biminRight, temp, MinMoment, TestMoment;
  if (M->N < N)
    return;
  if (N == 0)
  {
    M->N = 0;
    return;
  }
  if (N == 1)
  {
    M->a[0] *= M->b[0];
    for (i = 1; i < M->N; i++)
    {
      M->a[0] += M->b[i] * M->a[i];
    }
    for (i = 1; i < M->N; i++)
    {
      M->b[0] += M->b[i];
    }
    M->a[0] /= (M->b[0]);
    M->N = 1;
    return;
  }
  OrderListOfPoles(M);
  while (M->N > N)
  {
    imin = 0;
    MinMoment = M->b[imin] * (M->a[imin + 1] - M->a[imin]);
    TestMoment = M->b[M->N - 1] * (M->a[M->N - 1] - M->a[M->N - 2]);
    if (TestMoment < MinMoment)
    {
      MinMoment = TestMoment;
      imin = M->N - 1;
    }
    for (i = 1; i < M->N - 1; i++)
    {
      TestMoment = M->b[i] * Min(M->a[i + 1] - M->a[i], M->a[i] - M->a[i - 1]);
      if (TestMoment < MinMoment)
      {
        imin = i;
        MinMoment = TestMoment;
      }
    }
    if (M->b[imin] > DBL_MIN)
    {
      if (imin == 0)
      {
        M->a[imin + 1] =
          ((M->b[imin] / M->b[imin + 1]) * M->a[imin] + M->a[imin + 1]) / ((M->b[imin] / M->b[imin + 1]) + 1);
        M->b[imin + 1] = M->b[imin] + M->b[imin + 1];
      }
      else if (imin == M->N - 1)
      {
        M->a[imin - 1] =
          ((M->b[imin] / M->b[imin - 1]) * M->a[imin] + M->a[imin - 1]) / ((M->b[imin] / M->b[imin - 1]) + 1);
        M->b[imin - 1] = M->b[imin] + M->b[imin - 1];
      }

      else
      {
        if (fabs((M->a[imin + 1] - M->a[imin - 1])) < fabs(10 * (DBL_MIN + DBL_EPSILON * M->a[imin])))
        {
          biminLeft = M->b[imin] * (0.5);
          biminRight = M->b[imin] * (0.5);
        }
        else if (fabs((M->a[imin + 1] - M->a[imin])) < fabs(10 * (DBL_MIN + DBL_EPSILON * M->a[imin])))
        {
          biminLeft = M->b[imin];
          biminRight = 0;
        }
        else if (fabs((M->a[imin] - M->a[imin - 1])) < fabs(10 * (DBL_MIN + DBL_EPSILON * M->a[imin])))
        {
          biminLeft = 0;
          biminRight = M->b[imin];
        }
        else
        {
          biminLeft = M->b[imin] * ((M->a[imin + 1] - M->a[imin]) / (M->a[imin + 1] - M->a[imin - 1]));
          biminRight = M->b[imin] * ((M->a[imin] - M->a[imin - 1]) / (M->a[imin + 1] - M->a[imin - 1]));
        }
        M->a[imin + 1] =
          ((biminRight / M->b[imin + 1]) * M->a[imin] + M->a[imin + 1]) / ((biminRight / M->b[imin + 1]) + 1);
        M->b[imin + 1] = biminRight + M->b[imin + 1];
        M->a[imin - 1] =
          ((biminLeft / M->b[imin - 1]) * M->a[imin] + M->a[imin - 1]) / ((biminLeft / M->b[imin - 1]) + 1);
        M->b[imin - 1] = biminLeft + M->b[imin - 1];
      }
    }
    for (i = imin; i < M->N - 1; i++)
    {
      M->a[i] = M->a[i + 1];
      M->b[i] = M->b[i + 1];
    }
    M->N--;
    i = imin;
    while ((i > 1) && (M->a[i - 1] < M->a[i - 2]))
    {
      temp = M->a[i - 1];
      M->a[i - 1] = M->a[i - 2];
      M->a[i - 2] = temp;
      temp = M->b[i - 1];
      M->b[i - 1] = M->b[i - 2];
      M->b[i - 2] = temp;
      i--;
    }
    i = imin;
    while ((i < M->N - 1) && (M->a[i] > M->a[i + 1]))
    {
      temp = M->a[i];
      M->a[i] = M->a[i + 1];
      M->a[i + 1] = temp;
      temp = M->b[i];
      M->b[i] = M->b[i + 1];
      M->b[i + 1] = temp;
      i++;
    }
  }
}

/*  Create a new set of poles, wich are minimally spaced by da, by combining poles in a fixed region.
 *
 */
void ReduceListOfPolesSpacing(ListOfPoles * M, double da)
{
  uint32_t iRed, i;
  double PoleBorder;
  OrderListOfPoles(M);
  iRed = 0;
  PoleBorder = M->a[iRed] + da / 2.0;
  for (i = 1; i < M->N; i++)
  {
    if (PoleBorder > M->a[i])
    {
      if (fabs(M->b[i] + M->b[iRed]) < 10 * DBL_MIN)
      {
        M->b[iRed] = 0;
        M->a[iRed] = M->a[Min(i + 1, M->N - 1)];
        M->b[i] = 0;
      }
      else
      {
        M->a[iRed] = (M->a[iRed] * M->b[iRed] + M->a[i] * M->b[i]) / (M->b[i] + M->b[iRed]);
        M->b[iRed] += M->b[i];
        M->b[i] = 0;
      }
    }
    else
    {
      iRed++;
      if (iRed != i)
      {
        M->a[iRed] = M->a[i];
        M->b[iRed] = M->b[i];
        M->b[i] = 0;
      }
      PoleBorder += da;
    }
  }
  M->N = iRed + 1;
}

uint32_t ListOfPolesMoveNegativeEnergyToMu(ListOfPoles * const A, double const mu)
{
  uint32_t i, npoleshift = 0;
  double weightpoleshift = 0, delta = 0, momentshift = 0, totw = 0;
  for (i = 0; i < A->N; i++)
  {
    totw += A->b[i];
    if (A->a[i] < mu)
    {
      delta = mu - A->a[i];
      A->a[i] = mu * (1 + 10 * DBL_EPSILON) + 10 * DBL_MIN;
      npoleshift++;
      weightpoleshift += A->b[i];
      momentshift += A->b[i]*delta;
    }
  }
  if (npoleshift)
  {
    printf("******************** WARNING - NOT CONVERGED ********************\n");
    printf("List of poles %s, has %i poles of which %i are incorrectly"
           " at negative energies.\n",
           A->Name, A->N, npoleshift);
    printf("These have a total weight of %22.15E, (outof %22.15E)\n"
           "with moment %22.15E which has been shifted to mu=%22.15E\n",
           weightpoleshift, totw, momentshift, mu);
  }
  else
  {
    printf("List of poles %s, has all poles at positive energy"
           " with weight %22.15E.\n",
           A->Name, totw);
  }
  return 0;
}

uint32_t ListOfPolesMovePositiveEnergyToMu(ListOfPoles * const A, double const mu)
{
  uint32_t i, npoleshift = 0;
  double delta, weightpoleshift = 0, momentshift = 0, totw = 0;
  for (i = 0; i < A->N; i++)
  {
    totw += A->b[i];
    if (A->a[i] > mu)
    {
      delta = A->a[i] - mu;
      A->a[i] = mu * (1 - 10 * DBL_EPSILON) - 10 * DBL_MIN;
      npoleshift++;
      weightpoleshift += A->b[i];
      momentshift += A->b[i]*delta;
    }
  }
  if (npoleshift)
  {
    printf("******************** WARNING - NOT CONVERGED ********************\n");
    printf("List of poles %s, has %i poles of which %i are incorrectly"
           " at positive energies.\n",
           A->Name, A->N, npoleshift);
    printf("These have a total weight of %22.15E, (outof %22.15E)\n"
           " with moment %22.15E which has been shifted to mu=%22.15E\n",
           weightpoleshift, totw, momentshift, mu);
  }
  else
  {
    printf("List of poles %s, has all poles at negative energy"
           " with weight %22.15E.\n",
           A->Name, totw);
  }
  return 0;
}


/*  Given a list of poles, all poles with energy smaller than mu will be removed, the other poles will be rescaled such that the total weight is conserved
 *
 */
uint32_t ListOfPolesRemoveNegativeEnergy(ListOfPoles * const A, double const mu)
{
  uint32_t i, NPos = 0, NNeg = 0;
  double totalweight = 0, totalweightPos = 0, totalweightNeg = 0, scale;
  for (i = 0; i < A->N; i++)
  {
    totalweight += A->b[i];
    if (A->a[i] >= mu)
    {
      totalweightPos += A->b[i];
      NPos++;
    }
    else
    {
      totalweightNeg += A->b[i];
      A->b[i] = 0;
      NNeg++;
    }
  }
  if (totalweight < 10 * DBL_EPSILON)
  {
    /* happy situation where the function is zero or almost zero */
  }
  else if (totalweightPos < 10 * DBL_EPSILON)
  {
    /* The spectrum is there, but most probably has no meaning, i.e. is totally given by numerical noice */
    A->a[0] = 2 * mu;
    A->b[0] += totalweightNeg;
    printf("List of poles %s, has %i poles of which %i are at negative energies.\n", A->Name, A->N, NNeg);
    printf("Total weight %22.15E, Weight at negative energy %22.15E Weigth at positive energy %22.15E\n", totalweight,
           totalweightNeg, totalweightPos);
    printf("Based on this there is most probably no meaning other than numerical noice in this spectrum\n");
  }
  else
  {
    scale = totalweight / totalweightPos;
    for (i = 0; i < A->N; i++)
    {
      A->b[i] *= scale;
    }
    printf("List of poles %s, has %i poles of which %i are incorrectly at negative energies.\n", A->Name, A->N, NNeg);
    printf("Total weight %22.15E, Weight at negative energy %22.15E\n", totalweight, totalweightNeg);
  }
  return 0;
}

/*  Given a list of poles, all poles with energy larger than mu will be removed, the other poles will be rescaled such that the total weight is conserved
 *
 */
uint32_t ListOfPolesRemovePositiveEnergy(ListOfPoles * const A, double const mu)
{
  uint32_t i, NPos = 0, NNeg = 0;
  double totalweight = 0, totalweightPos = 0, totalweightNeg = 0, scale;
  for (i = 0; i < A->N; i++)
  {
    totalweight += A->b[i];
    if (A->a[i] > mu)
    {
      totalweightPos += A->b[i];
      NPos++;
      A->b[i] = 0;
    }
    else
    {
      totalweightNeg += A->b[i];
      NNeg++;
    }
  }
  if (totalweight < 10 * DBL_EPSILON)
  {
    /* happy situation where the function is zero or almost zero */
  }
  else if (totalweightNeg < 10 * DBL_EPSILON)
  {
    /* The spectrum is there, but most probably has no meaning, i.e. is totally given by numerical noice */
    A->a[0] = -2 * mu;
    A->b[0] += totalweightPos;
    printf("List of poles %s, has %i poles of which %i are at negative energies.\n", A->Name, A->N, NNeg);
    printf("Total weight %22.15E, Weight at negative energy %22.15E Weigth at positive energy %22.15E\n", totalweight,
           totalweightNeg, totalweightPos);
    printf("Based on this there is most probably no meaning other than numerical noice in this spectrum\n");
  }
  else
  {
    scale = totalweight / totalweightNeg;
    for (i = 0; i < A->N; i++)
    {
      A->b[i] *= scale;
    }
    printf("List of poles %s, has %i poles of which %i are incorrectly at positive energies.\n", A->Name, A->N, NPos);
    printf("Total weight %22.15E, Weigth at positive energy %22.15E\n", totalweight, totalweightPos);
  }
  return 0;
}


/*   A <- A-B
 *   The energy of the poles in B is adjusted to be on the same "grid" as the poles in A
 *   If A pole in A becomes negative weight it is substracted from the remainder of A
 *   If A has only one pole left the returned weight can be negative
 */
uint32_t SubtractListOfPoles(ListOfPoles * const A, const ListOfPoles B)
{
  uint32_t i, NRes, iMin;
  double VRes[A->N + B.N], WRes[A->N + B.N];
  if (A->N == 0)
  {
    printf("SubtractListOfPoles called with NA=0\n");
    return 1;
  };
  OrderListOfPoles(A);
  A->a0 -= B.a0;
  /* copy B to Res */
  for (i = 0; i < B.N; i++)
  {
    VRes[i] = B.a[i];
    WRes[i] = B.b[i];
  }
  NRes = B.N;
  /* subtract B from A */
  while (NRes != 0)
  {
    /* find neighbors in A for VRes[NRes-1] */
    /* special case *NA=1 */
    if (A->N == 1)
    {
      if (fabs(A->b[0] - WRes[NRes - 1]) > (DBL_MIN + DBL_EPSILON * fabs(WRes[NRes - 1] + A->b[0])))
      {
        A->a[0] += WRes[NRes - 1] * (A->a[0] - VRes[NRes - 1]) / (A->b[0] - WRes[NRes - 1]);
      }
      A->b[0] -= WRes[NRes - 1];
      NRes--;
    }
    else
    {
      if (VRes[NRes - 1] > A->a[A->N - 2])
      {
        iMin = A->N - 2;
      }
      else
      {
        iMin = A->N - 2;
        for (i = 1; i < A->N; i++)
        {
          if (VRes[NRes - 1] < A->a[i])
          {
            iMin = i - 1;
            break;
          }
        }
      }
      if (fabs(A->a[iMin + 1] - A->a[iMin]) < (DBL_MIN + DBL_EPSILON * fabs(A->a[iMin + 1] - VRes[NRes - 1])))
      {
        A->b[iMin] -= 0.5 * WRes[NRes - 1];
        A->b[iMin + 1] -= 0.5 * WRes[NRes - 1];
      }
      else
      {
        A->b[iMin + 1] -= ((VRes[NRes - 1] - A->a[iMin]) / (A->a[iMin + 1] - A->a[iMin])) * WRes[NRes - 1];
        A->b[iMin] -= ((A->a[iMin + 1] - VRes[NRes - 1]) / (A->a[iMin + 1] - A->a[iMin])) * WRes[NRes - 1];
      }
      NRes--;
      if (A->b[iMin + 1] < 0)
      {
        WRes[NRes] = -A->b[iMin + 1];
        VRes[NRes] = A->a[iMin + 1];
        NRes++;
        for (i = iMin + 2; i < A->N; i++)
        {
          A->b[i - 1] = A->b[i];
          A->a[i - 1] = A->a[i];
        }
        A->N--;
      }
      if ((A->b[iMin] < 0) && !(A->N==1))
      {
        WRes[NRes] = -A->b[iMin];
        VRes[NRes] = A->a[iMin];
        NRes++;
        for (i = iMin + 1; i < A->N; i++)
        {
          A->b[i - 1] = A->b[i];
          A->a[i - 1] = A->a[i];
        }
        A->N--;
      }
    }
  }
  if (A->b[0] < 0)
  {
    printf("SubtractListOfPoles returned with A->b[0]=%24.17E and A->N=%i\n", A->b[0], A->N);
    return 2;
  };
  return 0;
}

/*  a list of poles will be Gaussian broadened
 *  returned is a new list of length 2Nmul+1
 *  
 */
uint32_t GaussianBroadenListOfPoles(ListOfPoles * const M, double const G, uint32_t const Nmul)
{
  uint32_t i, j, k;
  double Weight[Nmul + 1], sumWeight, dw = 3 * G / ((double)(2 * Nmul));
  double *temp;
  temp = (double*)realloc(M->a, (2 * Nmul + 1) * M->N * sizeof(double));
  if (temp == NULL)
  {
    printf("realloc failed in GaussianBroadenListOfPoles\n");
    return 1;
  }
  M->a = temp;
  temp = (double*)realloc(M->b, (2 * Nmul + 1) * M->N * sizeof(double));
  if (temp == NULL)
  {
    printf("realloc failed in GaussianBroadenListOfPoles\n");
    return 1;
  }
  M->b = temp;
  sumWeight = 0;
  for (i = 0; i <= Nmul; i++)
  {
    Weight[i] = Gaussian(i * dw, G);
    sumWeight += Weight[i];
  }
  for (i = 1; i <= Nmul; i++)
  {
    sumWeight += Weight[i];
  }
  for (i = 0; i <= Nmul; i++)
  {
    Weight[i] /= sumWeight;
  }
  k = M->N;
  for (i = 0; i < M->N; i++)
  {
    for (j = 1; j <= Nmul; j++)
    {
      M->b[k] = M->b[i] * Weight[j];
      M->a[k] = M->a[i] - j * dw;
      k++;
      M->b[k] = M->b[i] * Weight[j];
      M->a[k] = M->a[i] + j * dw;
      k++;
    }
    M->b[i] *= Weight[0];
  }
  M->N *= (2 * Nmul + 1);
  return 0;
}

uint32_t ListOfPolesCopy(ListOfPoles * M, ListOfPoles B)
{
  uint32_t i;
  M->N = B.N;
  if (InitListOfPoles(M))
  {
    printf("InitListOfPoles failed in ListOfPolesCopy\n");
    return 1;
  }
  M->a0 = B.a0;
  M->mu = B.mu;
  strcpy(M->Name, B.Name);
  for (i = 0; i < B.N; i++)
  {
    M->a[i] = B.a[i];
    M->b[i] = B.b[i];
  }
  return 0;
}


/* copies all poles with energy smaller than Ef-Epsilon
 * Removes all poles with energy larger than Ef+Epsilon
 * Keeps the poles in the interval Ef-Epsilon to Ef+Epsilon, but scales the weight with ( (Ef+Epsilon) - Epole ) / (2 Epsilon)) 
 */
uint32_t ListOfPolesCopyVal(ListOfPoles * M, ListOfPoles B, double Ef, double Epsilon)
{
  uint32_t i, npoles;
  if (ListOfPolesCopy(M, B))
  {
    printf("ListOfPolesCopy failed in ListOfPolesCopyVal\n");
    return 1;
  }
  npoles = B.N;
  for (i = 0; i < B.N; i++)
  {
    if (M->a[i] >= Ef + Epsilon)
    {
      M->b[i] = 0;
      npoles--;
    }
    else if (M->a[i] > Ef - Epsilon)
    {
      M->b[i] *= ((Ef + Epsilon) - M->a[i]) / (2 * Epsilon);
    }
  }
  ReduceListOfPolesNumber(M, npoles);
  return 0;
}


/* copies all poles with energy larger than Ef+Epsilon
 * Removes all poles with energy smaller than Ef-Epsilon
 * Keeps the poles in the interval Ef-Epsilon to Ef+Epsilon, but scales the weight with ( Epole - (Ef-Epsilon) ) / (2 Epsilon))
 */
uint32_t ListOfPolesCopyCon(ListOfPoles * M, ListOfPoles B, double Ef, double Epsilon)
{
  uint32_t i, npoles;
  if (ListOfPolesCopy(M, B))
  {
    printf("ListOfPolesCopy failed in ListOfPolesCopyVal\n");
    return 1;
  }
  npoles = B.N;
  for (i = 0; i < B.N; i++)
  {
    if (M->a[i] <= Ef - Epsilon)
    {
      M->b[i] = 0;
      npoles--;
    }
    else if (M->a[i] < Ef + Epsilon)
    {
      M->b[i] *= (M->a[i] - (Ef - Epsilon)) / (2 * Epsilon);
    }
  }
  ReduceListOfPolesNumber(M, npoles);
  return 0;
}

/*  Keeps the N poles with the smallest energy such that the sum of the
 *  weight of these poles is exxactly equal to occupation.
 *  The weight of the N'th pole with the largest energy is reduced to
 *  exactly match the occupation needed.
 *  Ef is set to the energy of the last pole (partiall) included
 */
uint32_t ListOfPolesCopyValOccupation(ListOfPoles * M, ListOfPoles B, double *Ef, double occupation)
{
  double occupationSum;
  uint32_t i;
  if (ListOfPolesCopy(M, B))
  {
    printf("ListOfPolesCopy failed in ListOfPolesCopyValOccupation\n");
    return 1;
  }
  OrderListOfPoles(M);
  occupationSum = M->b[0];
  i = 0;
  while ((occupationSum < occupation) && (i < M->N - 1))
  {
    i++;
    occupationSum += M->b[i];
  }
  if (i == M->N - 1)
  {
    if (occupationSum < occupation)
    {
      printf
        ("In ListOfPolesCopyValOccupation occupation (%22.15E) is larger than the total integral of the list of poles (%22.15E)\n",
         occupation, occupationSum);
      printf("The entire list is returned\n");
      return 0;
    }
  }
  M->b[i] += occupation - occupationSum;
  *Ef = M->a[i];
  M->N = i + 1;
  return 0;
}

uint32_t ListOfPolesCopyConOccupation(ListOfPoles * M, ListOfPoles B, double *Ef, double occupation)
{
  double occupationSum;
  uint32_t i, j;
  if (ListOfPolesCopy(M, B))
  {
    printf("ListOfPolesCopy failed in ListOfPolesCopyConOccupation\n");
    return 1;
  }
  OrderListOfPoles(M);
  occupationSum = M->b[0];
  i = 0;
  while ((occupationSum < occupation) && (i < M->N - 1))
  {
    i++;
    occupationSum += M->b[i];
  }
  if (i == M->N - 1)
  {
    if (occupationSum < occupation)
    {
      printf
        ("In ListOfPolesCopyConOccupation occupation (%22.15E) is larger than the total integral of the list of poles (%22.15E)\n",
         occupation, occupationSum);
      printf("A zero list is returned\n");
      M->N = 0;
      return 0;
    }
  }
  M->b[0] = occupationSum - occupation;
  M->a[0] = M->a[i];
  *Ef = M->a[i];
  for (j = i + 1; j < M->N; j++)
  {
    M->b[j - i] = M->b[j];
    M->a[j - i] = M->a[j];
  }
  M->N -= i;
  return 0;
}

uint32_t ListOfPolesDetermineOccupation(ListOfPoles M, double const Ef, double *const occupation)
{
  uint32_t i;
  *occupation = 0;
  for (i = 0; i < M.N; i++)
  {
    if (M.a[i] < Ef)
      *occupation += M.b[i];
  }
  return 0;
}

uint32_t ListOfPoleCopyChoped(ListOfPoles M, ListOfPoles * A, double Epsilon)
{
  uint32_t i, shift;
  ListOfPolesCopy(A, M);
  shift = 0;
  for (i = 0; i < A->N; i++)
  {
    if (A->b[i] < Epsilon)
    {
      shift++;
    }
    else
    {
      if (shift <= i)
      {
        A->b[i - shift] = A->b[i];
        A->a[i - shift] = A->a[i];
      }
    }
  }
  A->N -= shift;
  return 0;
}

/*  Shifts the list of poles such that the chemical potential (Ef at T=0) is zero
 *  There is one pole that sits exactly at Ef. This pole is split in two
 *  one part is shifted down one part up with a total splitting of 1000 DBL_EPSILON
 *  The weight is devided such that Ef sits in the middle
 */
uint32_t ShiftListOfPolesToMakeChemPotZeroMetal(ListOfPoles * const M, double const occupation)
{
  uint32_t i;
  double occSum, occSumextra;
  double Ef;
  OrderListOfPoles(M);
  occSum = 0;
  i = 0;
  while ((occSum <= occupation) && (i < M->N))
  {
    occSum += M->b[i];
    i++;
  }
  if (occSum < occupation)
  {
    printf
      ("occupation larger than available in ShiftListOfPolesToMakeChemPotZeroMetal (occupation %22.15E, sum %22.15E)\n",
       occupation, occSum);
    return 1;
  };
  if (i == 0)
  {
    printf("occupation negative in ShiftListOfPolesToMakeEfZero (%22.15E)\n", occupation);
    return 1;
  };
  i--;
  occSumextra = occSum - occupation;
  Ef = M->a[i];
  M->a = (double*)realloc(M->a, (M->N + 1) * sizeof(double));
  M->b = (double*)realloc(M->b, (M->N + 1) * sizeof(double));
  M->a[M->N] = M->a[i] + 1000 * DBL_EPSILON * (M->b[i] - occSumextra) / (M->b[i] + DBL_MIN);
  M->a[i] -= 1000 * DBL_EPSILON * (occSumextra / (M->b[i] + DBL_MIN));
  M->b[M->N] = occSumextra;
  M->b[i] -= occSumextra;
  M->N += 1;
  OrderListOfPoles(M);
  for (i = 0; i < M->N; i++)
  {
    M->a[i] -= Ef;
  }
  return 0;
}

/*  Shifts the list of poles such that the chemical potential (EF at T=0) is zero
 *  Finds the pole at Ef, splits it in two such that the weights are ok below and above Ef
 *  Then looks for the largest distance to the neighbor pool and merges the spectral weight on this side of the pole
 */
uint32_t ShiftListOfPolesToMakeChemPotZeroInsulator(ListOfPoles * const M, double const occupation)
{
  uint32_t i;
  double occSum, occSumextra, occShift;
  double Ef;
  OrderListOfPoles(M);
  occSum = 0;
  i = 0;
  while ((occSum <= occupation) && (i < M->N))
  {
    occSum += M->b[i];
    i++;
  }
  if (occSum < occupation)
  {
    printf
      ("occupation larger than available in ShiftListOfPolesToMakeChemPotZeroInsulator (occupation %22.15E, sum %22.15E)\n",
       occupation, occSum);
    return 1;
  };
  if (i == 0)
  {
    printf("occupation negative in ShiftListOfPolesToMakeEfZero (%22.15E)\n", occupation);
    return 1;
  };
  i--;
  occSumextra = occSum - occupation;
  occShift = M->b[i] - occSumextra;
  M->a[i - 1] = (M->a[i - 1] * M->b[i - 1] + M->a[i] * occShift) / (M->b[i - 1] + occShift + 1000 * DBL_MIN);
  M->b[i - 1] += occShift;
  M->b[i] = occSumextra;
  Ef = (M->a[i] + M->a[i - 1]) / 2;
  for (i = 0; i < M->N; i++)
  {
    M->a[i] -= Ef;
  }
  return 0;
}

/*  Shifts the energy of the list of poles such that the chemical potential (Ef at T=0) is zero
 *  and removes all poles smaler than weight epsilon (set to 1000*DBL_EPSILON or so)
 *  first assumes a mettal and shifts
 *  then removes small poles keeping Chemical potential fixed
 *  thirdly shifts again assuming an insulator (Ef is now between two poles also for a metal)
 */
uint32_t ShiftListOfPolesToMakeChemPotZero(ListOfPoles * const M, double const occupation, double const epsilon)
{
  ShiftListOfPolesToMakeChemPotZeroMetal(M, occupation);
  ReduceListOfPolesWeigth(M, epsilon, 0);
  ShiftListOfPolesToMakeChemPotZeroInsulator(M, occupation);
  return 0;
}

/*
 *
 */
uint32_t ShiftListOfPoles(ListOfPoles * const M, double const munew)
{
  uint32_t i;
  M->mu -= munew;
  for (i = 0; i < M->N; i++)
  {
    M->a[i] += munew;
  }
  return 0;
}


/*  Shift list of poles to make first moment zero
 *
 */
uint32_t ShiftListOfPolesToMakeFirstMomentZero(ListOfPoles * const M)
{
  uint32_t i;
  double moment;
  moment = 0;
  for (i = 0; i < M->N; i++)
  {
    moment += M->b[i] * M->a[i];
  }
  for (i = 0; i < M->N; i++)
  {
    M->a[i] -= moment;
  }
  printf("List of poles %s shifted by the first moment %22.15E\n", M->Name, moment);
  return 0;
}

uint32_t ListOfPolesCalculateChemicalPotential(ListOfPoles const M, double const occupation, double *const mu)
{
  uint32_t i;
  double occsum = 0;
  for (i = 0; i < M.N; i++)
  {
    occsum += M.b[i];
    if (occsum >= occupation)
    {
      *mu = M.a[i];
      return 0;
    }
  }
  *mu = M.a[M.N - 1];
  printf("Total occupation of %s is %22.15E\n", M.Name, occsum);
  printf("Searched for chemical potential of filling of %22.15E\n", occupation);
  printf("that aint working. Returned %22.15E as chemical potential\n", *mu);
  return 1;
}
