// Solves the (periodic) harmonic oscillator H = -0.5*(d^2/dx^2) + 0.5*x*x
// by using plane waves as basis functions
// The harmonic oscillator 

#define NOSCI 70
#define MOSCI (2*NOSCI+1)
const double a = 10.;
const double dxint = 0.01;

// normalized basis function
double complex inline PsiFun(double x, int i)
{
  double k = 2. * M_PI / a * (double)i;
  return (cos(k * x) + sin(k * x) * I) / sqrt(a);
}

// Evaluate <Psi_i| -0.5*d²/dx² |Psi_j>
double complex IntegrateKinetic(int i, int j)
{
  double k2 = 2. * M_PI / a * (double)j;
  double complex sum = 0.;
  double x;
  for (x = -a / 2.; x <= a / 2.; x += dxint)
  {
    sum += conj(PsiFun(x, i)) * (-k2 * k2) * PsiFun(x, j) * dxint;  // f1*f2*dxint;
  }
  return -sum * 0.5;
}

// Evaluate <Psi_i| 0.5*x² |Psi_j>
double complex IntegratePotential(int i, int j)
{
  double complex sum = 0.;
  double x;
  for (x = -a / 2.; x <= a / 2.; x += dxint)
  {
    double V = 0.5 * x * x;
    sum += conj(PsiFun(x, i)) * V * PsiFun(x, j) * dxint;
  }
  return sum;
}

void SolveBLAS()
{
  int    i, j;

  double matrix[MOSCI * MOSCI];
  double e[MOSCI];
  double wrk[3 * MOSCI];

  for (j = -NOSCI; j <= NOSCI; j++)
    for (i = -NOSCI; i <= NOSCI; i++)
      matrix[(j + NOSCI) * MOSCI + (i + NOSCI)] = IntegrateKinetic(i, j) + IntegratePotential(i, j);

  const char *V = "V";
  const char *U = "U";
  int    npw = MOSCI;
  int    info;
  int    wrksize = 3 * MOSCI;
  printf("Solution based on dsyev (BLAS):\n");
  dsyev(V, U, &npw, matrix, &npw, e, wrk, &wrksize, &info);
  printf("e0 = %f\n", e[0]);
  printf("e1 = %f\n", e[1]);
  printf("e2 = %f\n", e[2]);
}


void SolveLanczos()
{
  int    i, j;
  int    NBitsKey = 8;
  OperatorType H;
  InitOperator(&H);
  H.NFermion = MOSCI;
  InitOperatorFixedLength(&H, 2, NBitsKey);
  for (j = -NOSCI; j <= NOSCI; j++)
  {
    for (i = -NOSCI; i <= NOSCI; i++)
    {
      uint16_t oadd[2];
      oadd[0] = (i + NOSCI);
      oadd[1] = (j + NOSCI) | CrAnMask;
      RealOperatorAddLader(&H, oadd, 2, IntegrateKinetic(i, j) + IntegratePotential(i, j));
    }
  }

  WaveFunctionType psi;
  CreateRealRandomWaveFunction(&psi, MOSCI, 1);
  WaveFunctionRandomize(&psi);
  //PrintWaveFunctionReadable(psi);

  int    n = 3;
  double epsilon = 100 * DBL_EPSILON;
  WaveFunctionType psiEigen[n];
  printf("\nSolution based on Lanczos: \n");

  LanczosEigenSystem(H, psi, psiEigen, &epsilon, 1000 * DBL_EPSILON, n);

  for (i = 0; i < n; i++)
  {
    double ereal, eimag;
    PsiOperatorPsi(H, psiEigen[i], &ereal, &eimag);
    printf("e%i = %f\n", i, ereal);
  }
  //PrintWaveFunctionReadable(psiEigen[0]);
  //PrintOperatorReadable(op);
  //PrintOperatorFullReadable(op);
  FreeOperator(&H);
}


void TestHarmonicOscillator()
{
  SolveBLAS();
  SolveLanczos();
}
