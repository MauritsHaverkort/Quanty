
double Factorial(int a)
{
  int    i;
  double result;
  result = 1.0;
  if (a < 0)
  {
    printf("Error Factorial called with a negative number\n");
    return 0;
  }
  for (i = a; i > 0; i--)
  {
    result *= i;
  }
  return result;
}

double Factorial2(int a)
{
  int    i;
  double result;
  result = 1.0;
  if (a < 0)
  {
    printf("Error Factorial2 called with a negative number\n");
    return 0;
  }
  for (i = a; i > 0; i-=2)
  {
    result *= i;
  }
  return result;
}

/*  computes n over k (the number of ways k elements can be placed in n spaces)
 *
 */
uint32_t BinomialCoefficient(uint32_t n, uint32_t k)
{
  if(n<k) return 0;
  uint32_t result = 1, i;
  for (i = 1 + Max(k, n - k); i <= n; i++)
  {
    result *= i;
  }
  for (i = 1; i <= Min(k, n - k); i++)
  {
    result /= i;
  }
  return result;
}
double BinomialCoefficientDouble(uint32_t n, uint32_t k)
{
  if(n<k) return 0;
  double result = 1, i;
  for (i = 1 + Max(k, n - k); i < n+0.5; i++)
  {
    result *= i;
  }
  for (i = 1; i < Min(k, n - k)+0.5; i++)
  {
    result /= i;
  }
  return result;
}

/*
 *  Given NV operators labeled V_1 to V_{NV}
 *  How manny ways are there to create an operator of order "order"
 */
uint32_t PerturbationOptions(uint32_t NV, uint32_t order)
{
  uint32_t result = 1, i;
  if (NV == 0)
    return 0;
  for (i = 1; i < NV; i++)
  {
    result *= i + order;
    result /= i;
  }
  return result;
}

/*
 *  Given NV operators labeled V_{1} to V_{NV}
 *  Given that for a given string I have "order" operators
 *  Given a list of how many of each operators I have (i.e. NV integers, the sum over this list is "order")
 *  retrun a unique index position ranging from 0 to PerturbationOptions(NV,order)-1
 */
void PerturbationOrderListToIndex(uint32_t * const POList, uint32_t * const POIndex, uint32_t const NV,
                                  uint32_t const order)
{
  uint32_t i, POIndexSmaller;
  *POIndex = 0;
  if ((POList[0] == order) || (NV <= 1))
    return;
  for (i = 0; i < order - POList[0]; i++)
  {
    *POIndex += PerturbationOptions(NV - 1, i);
  }
  PerturbationOrderListToIndex(&(POList[1]), &POIndexSmaller, NV - 1, order - POList[0]);
  *POIndex += POIndexSmaller;
}

/*
 *  Given NV operators labeled V_{1} to V_{NV}
 *  Given that for a given string I have order operators
 *  Given a unique index position ranging from 0 to PerturbationOptions(NV,order)-1
 *  return list of how many of each operators I have (i.e. NV integers)
 */
void PerturbationOrderIndexToList(uint32_t const POIndex, uint32_t * const POList, uint32_t const NV, uint32_t const order)
{
  uint32_t i, nl, POIndexnew, deg;
  if (NV == 1)
  {
    POList[0] = order;
    return;
  }
  nl = 0;
  POIndexnew = POIndex;
  for (i = 0; i <= order; i++)
  {
    deg = PerturbationOptions(NV - 1, i);
    nl += deg;
    if (POIndex < nl)
      break;
    POIndexnew -= deg;
  }
  POList[0] = order - i;
  PerturbationOrderIndexToList(POIndexnew, &(POList[1]), NV - 1, i);
}



/*
 *
 */
void AddPerturbationIndex(uint32_t const Vi, uint32_t const IndexOld, uint32_t const NV, uint32_t const order, uint32_t * const IndexNew, int *const phase)
{
  uint32_t i, List[NV];
  PerturbationOrderIndexToList(IndexOld, List, NV, order);
  *phase = 0;
  for (i = 0; i < Vi; i++)
  {
    *phase += List[i];
  }
  *phase = *phase & 1 ? -1 : 1;
  List[Vi] += 1;
  PerturbationOrderListToIndex(List, IndexNew, NV, order + 1);
}

/*
 *
 */
void RemovePerturbationIndex(uint32_t const Vi, uint32_t const IndexOld, uint32_t const NV, uint32_t const order, uint32_t * const IndexNew, int *const phase)
{
  uint32_t i, List[NV];
  PerturbationOrderIndexToList(IndexOld, List, NV, order);
  *phase = 0;
  if (List[Vi] == 0)
    return;
  for (i = 0; i < Vi; i++)
  {
    *phase += List[i];
  }
  *phase = *phase & 1 ? -1 : 1;
  List[Vi] -= 1;
  PerturbationOrderListToIndex(List, IndexNew, NV, order - 1);
}

double ThreeJSymbol(int j1, int m1, int j2, int m2, int j3, int m3)
{
  double ThreeJ;
  int    k;
  ThreeJ = 0;
  if (m1 + m2 + m3 != 0)
    return 0;
  if (j1 + j2 < j3)
    return 0;
  if (j2 + j3 < j1)
    return 0;
  if (j3 + j1 < j2)
    return 0;
  if (abs(m1) > j1)
    return 0;
  if (abs(m2) > j2)
    return 0;
  if (abs(m3) > j3)
    return 0;
  for (k = Max(j1 - j3 + m2, Max(j2 - j3 - m1, 0)); k <= Min(j1 - m1, Min(j2 + m2, j1 + j2 - j3)); k++)
  {
    ThreeJ +=
      (double)(1 - 2 * (k & 1)) / (Factorial(k) * Factorial(j1 + j2 - j3 - k) * Factorial(j1 - m1 - k) * Factorial(j2 + m2 - k) *
                               Factorial(j3 - j2 + m1 + k) * Factorial(j3 - j1 - m2 + k));
  }
  ThreeJ *=
    (double)(1 -
             2 * (Abs((j2 - j1 + m3)) & 1)) * sqrt(Factorial(j1 + j2 - j3) * Factorial(j1 - j2 + j3) * Factorial(-j1 +j2 + j3) *
                                                   Factorial(j1 + m1) * Factorial(j1 - m1) * Factorial(j2 + m2) *
                                                   Factorial(j2 - m2) * Factorial(j3 + m3) * Factorial(j3 - m3) /
                                                   Factorial(j1 + j2 + j3 + 1));
  return ThreeJ;
}

double ThreeJSymbolHalf(int j1, int m1, int j2, int m2, int j3, int m3)
{
  double ThreeJ;
  int    k;
  ThreeJ = 0;
  if (m1 + m2 + m3 != 0)
    return 0;
  if (j1 + j2 < j3)
    return 0;
  if (j2 + j3 < j1)
    return 0;
  if (j3 + j1 < j2)
    return 0;
  if (abs(m1) > j1)
    return 0;
  if (abs(m2) > j2)
    return 0;
  if (abs(m3) > j3)
    return 0;
  for (k = Max(j1 - j3 + m2, Max(j2 - j3 - m1, 0)) / 2; k <= Min(j1 - m1, Min(j2 + m2, j1 + j2 - j3)) / 2; k++)
  {
    ThreeJ +=
      (double)(1 -
               2 * (k & 1)) / (Factorial(k) * Factorial((j1 + j2 - j3) / 2 - k) * Factorial((j1 - m1) / 2 - k) * Factorial((j2 + m2) / 2 - k) *
                               Factorial((j3 - j2 + m1) / 2 + k) * Factorial((j3 - j1 - m2) / 2 + k));
  }
  ThreeJ *=
    (double)(1 - 2 * (((j2 - j1 + m3) / 2) & 1)) * sqrt(Factorial((j1 + j2 - j3) / 2) * Factorial((j1 - j2 + j3) / 2) *
                                                        Factorial((-j1 + j2 + j3) / 2) * Factorial((j1 + m1) / 2) *
                                                        Factorial((j1 - m1) / 2) * Factorial((j2 + m2) / 2) * Factorial((j2 - m2) / 2) *
                                                        Factorial((j3 + m3) / 2) * Factorial((j3 - m3) / 2) / Factorial((j1 + j2 + j3) / 2 + 1));
  return ThreeJ;
}

double SlaterCoefficientC(int k, int l1, int m1, int l2, int m2)
{
  if (Abs(m1 - m2) > k) return 0;
  return (double)(1 - 2 * (m1 & 1)) * sqrt((double)((2 * l1 + 1) * (2 * l2 + 1))) * ThreeJSymbol(l1, 0, k, 0, l2, 0) * ThreeJSymbol(l1, -m1,k,m1 -m2,l2,m2);
}

double SlaterCoefficientCjcoupling(int k, int l1, int twoj1, int twojz1, int l2, int twoj2, int twojz2)
{
  return (1-2*((l1+l2-1+(twojz1+twojz2)/2)&1))*sqrt((double)((twoj1+1)*(twoj2+1)))*
         ( ThreeJSymbolHalf(2*l1, twojz1+1, 1, -1, twoj1,-twojz1)*ThreeJSymbolHalf(2*l2, twojz2+1, 1, -1, twoj2,-twojz2)*SlaterCoefficientC(k, l1, (twojz1+1)/2, l2, (twojz2+1)/2) +
           ThreeJSymbolHalf(2*l1, twojz1-1, 1,  1, twoj1,-twojz1)*ThreeJSymbolHalf(2*l2, twojz2-1, 1,  1, twoj2,-twojz2)*SlaterCoefficientC(k, l1, (twojz1-1)/2, l2, (twojz2-1)/2) );
}

int cmpuint32t (const void * a, const void *b)
{
  return (*(uint32_t*)a<*(uint32_t*)b)?(1):(-1);
}
int cmpdouble (const void * a, const void *b)
{
  return (*(double*)a<*(double*)b)?(1):(-1);
}
