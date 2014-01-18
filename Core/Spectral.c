/*
*   Defines a set of functions that act on the structure spectra
*
*/


/*  Allocate a spectrum according to the settings as saved in Settings.
*
*/
int SpectrumAllocate(const SpectraType Settings, SpectraType * const G)
{
  G->N = Settings.N;
  G->Emin = Settings.Emin;
  G->Emax = Settings.Emax;
  G->dE = Settings.dE;
  G->EminPole = Settings.EminPole;
  G->EmaxPole = Settings.EmaxPole;
  G->Gamma = Settings.Gamma;
  G->ZerothMoment = 0;
  G->FirstMoment = 0;
  G->a[0] = 0;
  G->a[1] = 0;
  G->a[2] = 0;
  G->a[3] = 0;
  G->a[4] = 0;
  G->b[0] = 0;
  G->b[1] = 0;
  G->b[2] = 0;
  G->b[3] = 0;
  G->b[4] = 0;
  G->Re = calloc(G->N + 1, sizeof(double));
  if (G->Re == NULL)
  {
    printf("\nERROR\nSpectrumAllocate failed 01\n");
    return 1;
  }
  G->Im = calloc(G->N + 1, sizeof(double));
  if (G->Im == NULL)
  {
    free(G->Re);
    printf("\nERROR\nSpectrumAllocate falied 02\n");
    return 1;
  }
  return 0;
}

/* Free the memory of an allocated spectra
*
*/
void SpectrumFree(SpectraType * const G)
{
  G->N = 0;
  G->ZerothMoment = 0;
  G->FirstMoment = 0;
  G->a[0] = 0;
  G->a[1] = 0;
  G->a[2] = 0;
  G->a[3] = 0;
  G->a[4] = 0;
  G->b[0] = 0;
  G->b[1] = 0;
  G->b[2] = 0;
  G->b[3] = 0;
  G->b[4] = 0;
  free(G->Re);
  free(G->Im);
}


/* Copy a spectrum
*
*/
void CopySpectra(const SpectraType In, SpectraType * Out)
{
  uint32_t i;
  Out->N = In.N;
  Out->Emin = In.Emin;
  Out->Emax = In.Emax;
  Out->dE = In.dE;
  Out->EminPole = In.EminPole;
  Out->EmaxPole = In.EmaxPole;
  Out->Gamma = In.Gamma;
  Out->ZerothMoment = In.ZerothMoment;
  Out->FirstMoment = In.FirstMoment;
  Out->a[0] = In.a[0];
  Out->a[1] = In.a[1];
  Out->a[2] = In.a[2];
  Out->a[3] = In.a[3];
  Out->a[4] = In.a[4];
  Out->b[0] = In.b[0];
  Out->b[1] = In.b[1];
  Out->b[2] = In.b[2];
  Out->b[3] = In.b[3];
  Out->b[4] = In.b[4];
  for (i = 0; i <= Out->N; i++)
  {
    Out->Re[i] = In.Re[i];
    Out->Im[i] = In.Im[i];
  }
}

int SpectraScalarMultiply(SpectraType * const A, double const aRe, double const aIm)
{
  int    i;
  double temp;
  for (i = 0; i <= A->N; i++)
  {
    temp = A->Re[i];
    A->Re[i] = aRe * temp - aIm * A->Im[i];
    A->Im[i] = aRe * A->Im[i] + aIm * temp;
  }
  return 1;
}

/* A<-(A+B)/2
*  B<-(A+B)/2
*  It is assumed that:
*  A.N=B.N
*/
int AvverageSpectra(SpectraType * const A, SpectraType * const B)
{
  uint32_t i;
  double ABZ, ABF;
  if ((A->N) != (B->N))
  {
    printf("Average Spectra failed A.N=%i B.N=%i\n", A->N, B->N);
    return 1;
  }
  if ((A->Emin) != (B->Emin))
  {
    printf("Warning in Average Spectra A.Emin=%f B.Emin=%f\n", A->Emin, B->Emin);
  }
  if ((A->Emax) != (B->Emax))
  {
    printf("Warning in Average Spectra A.Emax=%f B.Emax=%f\n", A->Emax, B->Emax);
  }
  if ((A->Gamma) != (B->Gamma))
  {
    printf("Warning in Average Spectra A.Gamma=%f B.Gamma=%f\n", A->Gamma, B->Gamma);
  }
  for (i = 0; i <= A->N; i++)
  {
    A->Re[i] += B->Re[i];
    A->Im[i] += B->Im[i];
    A->Re[i] *= 0.5;
    A->Im[i] *= 0.5;
    B->Re[i] = A->Re[i];
    B->Im[i] = A->Im[i];
  }
  ABZ = (A->ZerothMoment + B->ZerothMoment) / 2.0;
  ABF = (A->ZerothMoment * A->FirstMoment + B->ZerothMoment * B->FirstMoment) / (2 * ABZ);
  A->ZerothMoment = ABZ;
  B->ZerothMoment = ABZ;
  A->FirstMoment = ABF;
  B->FirstMoment = ABF;
  A->EminPole = Min(A->EminPole, B->EminPole);
  B->EminPole = Min(A->EminPole, B->EminPole);
  A->EmaxPole = Max(A->EmaxPole, B->EmaxPole);
  B->EmaxPole = Max(A->EmaxPole, B->EmaxPole);
  /* TO DO recalculate a[0] ... a[4] and b[0] ... b[4] */
  return 0;
}

/* A->B
*  B->A
*/
void SwapSpectra(SpectraType * A, SpectraType * B)
{
  SpectraType C;
  C = *A;
  *A = *B;
  *B = C;
}

/*
*   A <- (alpha) A + (1-alpha) B
*  It is assumed that:
*  A.N=B.N
*/
int MixSpectra(SpectraType * const A, const SpectraType B, const double alpha)
{
  uint32_t i;
  double ABZ, ABF;
  if ((A->N) != (B.N))
  {
    printf("Mix Spectra failed A.N=%i B.N=%i\n", A->N, B.N);
    return 1;
  }
  if ((A->Emin) != (B.Emin))
  {
    printf("Warning in Mix Spectra A.Emin=%f B.Emin=%f\n", A->Emin, B.Emin);
  }
  if ((A->Emax) != (B.Emax))
  {
    printf("Warning in Mix Spectra A.Emax=%f B.Emax=%f\n", A->Emax, B.Emax);
  }
  if ((A->Gamma) != (B.Gamma))
  {
    printf("Warning in Mix Spectra A.Gamma=%f B.Gamma=%f\n", A->Gamma, B.Gamma);
  }
  for (i = 0; i <= A->N; i++)
  {
    A->Re[i] *= alpha;
    A->Im[i] *= alpha;
    A->Re[i] += (1 - alpha) * B.Re[i];
    A->Im[i] += (1 - alpha) * B.Im[i];
  }
  ABZ = (alpha * A->ZerothMoment + (1 - alpha) * B.ZerothMoment);
  ABF = (alpha * A->ZerothMoment * A->FirstMoment + (1 - alpha) * B.ZerothMoment * B.FirstMoment) / (ABZ);
  A->ZerothMoment = ABZ;
  A->FirstMoment = ABF;
  A->EminPole = Min(A->EminPole, B.EminPole);
  A->EmaxPole = Max(A->EmaxPole, B.EmaxPole);
  /* TO DO recalculate a[0] ... a[4] and b[0] ... b[4] */
  return 0;
}


/*
*   A <- A + (alpha) B
*  It is assumed that:
*  A.N=B.N
*/
int AddMultipliedSpectra(SpectraType * const A, const SpectraType B, const double alphaRe, const double alphaIm)
{
  uint32_t i;
  if ((A->N) != (B.N))
  {
    printf("Mix Spectra failed A.N=%i B.N=%i\n", A->N, B.N);
    return 1;
  }
  if ((A->Emin) != (B.Emin))
  {
    printf("Warning in Mix Spectra A.Emin=%f B.Emin=%f\n", A->Emin, B.Emin);
  }
  if ((A->Emax) != (B.Emax))
  {
    printf("Warning in Mix Spectra A.Emax=%f B.Emax=%f\n", A->Emax, B.Emax);
  }
  if ((A->Gamma) != (B.Gamma))
  {
    printf("Warning in Mix Spectra A.Gamma=%f B.Gamma=%f\n", A->Gamma, B.Gamma);
  }
  for (i = 0; i <= A->N; i++)
  {
    A->Re[i] += alphaRe * B.Re[i] - alphaIm * B.Im[i];
    A->Im[i] += alphaRe * B.Im[i] + alphaIm * B.Re[i];
  }
  A->EminPole = Min(A->EminPole, B.EminPole);
  A->EmaxPole = Max(A->EmaxPole, B.EmaxPole);
  /* TO DO recalculate a[0] ... a[4] and b[0] ... b[4] phase, ZeroMoment and FirstMoment */
  return 0;
}


/*
*   A <- A - B
*  It is assumed that:
*  A.N=B.N
*/
int SubtractSpectra(SpectraType * const A, const SpectraType B)
{
  double ABZ, ABF;
  uint32_t i;
  if ((A->N) != (B.N))
  {
    printf("Mix Spectra failed A.N=%i B.N=%i\n", A->N, B.N);
    return 1;
  }
  if ((A->Emin) != (B.Emin))
  {
    printf("Warning in Mix Spectra A.Emin=%f B.Emin=%f\n", A->Emin, B.Emin);
  }
  if ((A->Emax) != (B.Emax))
  {
    printf("Warning in Mix Spectra A.Emax=%f B.Emax=%f\n", A->Emax, B.Emax);
  }
  if ((A->Gamma) != (B.Gamma))
  {
    printf("Warning in Mix Spectra A.Gamma=%f B.Gamma=%f\n", A->Gamma, B.Gamma);
  }
  for (i = 0; i <= A->N; i++)
  {
    A->Re[i] -= B.Re[i];
    A->Im[i] -= B.Im[i];
  }
  ABZ = (A->ZerothMoment - B.ZerothMoment);
  ABF = 0;                      /* TO DO */
  A->ZerothMoment = ABZ;
  A->FirstMoment = ABF;
  A->EminPole = Min(A->EminPole, B.EminPole);
  A->EmaxPole = Max(A->EmaxPole, B.EmaxPole);
  /* TO DO recalculate a[0] ... a[4] and b[0] ... b[4] */
  return 0;
}


/*
*   A <- 0
*  It is assumed that:
*  A.N=B.N
*/
int SetSpectraToZero(SpectraType * const A)
{
  uint32_t i;
  for (i = 0; i <= A->N; i++)
  {
    A->Re[i] = 0;
    A->Im[i] = 0;
  }
  A->ZerothMoment = 0;
  A->FirstMoment = 0;
  for (i = 0; i < 5; i++)
  {
    A->a[i] = 0;
    A->b[i] = 0;
  }
  /* TO DO recalculate a[0] ... a[4] and b[0] ... b[4] */
  return 0;
}



double CalcualteError(const SpectraType A, const SpectraType B)
{
  uint32_t i;
  double error = 0;
  for (i = 0; i <= A.N; i++)
  {
    error += Sqr(A.Im[i] - B.Im[i]);
  }
  error *= A.dE;
  error = sqrt(error);
  return error;
}



/*
*   Calculate the product of an matrix in Anderson impurity model form times a vector 
*/
void AndersonModelOperatorPsi(const TriDiagonalMatrixType M, const double *const vin, double *const vout)
{
  uint32_t i;
  vout[0] = M.a[0] * vin[0];
  for (i = 1; i <= M.N; i++)
  {
    vout[i] = M.a[i] * vin[i];
    vout[i] += M.b[i] * vin[0];
    vout[0] += M.b[i] * vin[i];
  }
}

/*
*   Calculate the norm of a vectr of length N+1 (0...N)
*/
double VectorCalculateNorm(const double *const v, const uint32_t N)
{
  double Norm = 0;
  uint32_t i;
  for (i = 0; i <= N; i++)
  {
    Norm += Sqr(v[i]);
  }
  Norm = sqrt(Norm);
  return Norm;
}

/*
*   Calculate the integral of a vectr of length N+1 (0...1)
*/
double VectorCalculateIntegral(const double *const v, const double dE, const uint32_t N)
{
  double Integral = (v[0] + v[N]) / 2.0;  /* Assume that dv/dE =0 at E=0 or E=N; the error is of the order dv/dE[0]+dv/dE[N], which have oposite sign and are both of the same small magnitude */
  uint32_t i;
  for (i = 1; i < N; i++)
  {
    Integral += v[i] + (1.0 / 24.0) * (v[i - 1] - 2 * v[i] + v[i + 1]); /* use a spline interpolation */
  }
  return Integral * dE;
}


/*
*    v<-a*v
*/
void VectorMultiply(double *const v, const uint32_t N, const double a)
{
  uint32_t i;
  for (i = 0; i <= N; i++)
  {
    v[i] *= a;
  }
}

/*
*   v1<-v1+v2
*/
void VectorAdd(double *const v1, const double *const v2, const uint32_t N)
{
  uint32_t i;
  for (i = 0; i <= N; i++)
  {
    v1[i] += v2[i];
  }
}

/*
*   v1<-v1+v2^{T}
*/
void VectorAddReverse(double *const v1, const double *const v2, const uint32_t N)
{
  uint32_t i;
  for (i = 0; i <= N; i++)
  {
    v1[i] += v2[N - i];
  }
}


/*
*   v1<-v1-v2^{T}
*/
void VectorSubtractReverse(double *const v1, const double *const v2, const uint32_t N)
{
  uint32_t i;
  for (i = 0; i <= N; i++)
  {
    v1[i] -= v2[N - i];
  }
}

/*
*   v1<-v1+a*v2
*/
void VectorAddMultiplied(double *const v1, const double *const v2, const uint32_t N, const double a)
{
  uint32_t i;
  for (i = 0; i <= N; i++)
  {
    v1[i] += a * v2[i];
  }
}

/*
*   <-v1.v2
*/
double VectorDotProduct(const double *const v1, const double *const v2, const uint32_t N)
{
  uint32_t i;
  double Dot = 0;
  for (i = 0; i <= N; i++)
  {
    Dot += v1[i] * v2[i];
  }
  return Dot;
}

double SpectraCalculateFermiEnergy(const SpectraType G, const double Nelec)
{
  uint32_t i;
  double sumN, alpha;
  sumN = -G.ZerothMoment * (0.5 - atan(2 * (G.FirstMoment - G.Emin) / G.Gamma) / Pi);
  sumN -= G.Im[0] * G.dE * 0.5; /* The first point only has half the interval */
  i = 0;
  while (sumN < Nelec)
  {
    i++;
    sumN -= G.Im[i] * G.dE;
  }
  alpha = (sumN - Nelec) / (G.Im[i] * G.dE);
  return (G.Emin + (double)(i + 0.5 + alpha) * G.dE);
}


void RemovePositiveImPart(SpectraType * const G)
{
  double Gtemp[G->N + 1];       /* define a temporary array to store G with negative part 'broadenend' */
  double max, min;
  uint32_t i;
  Gtemp[0] = G->Im[0];
  Gtemp[G->N] = G->Im[G->N];
  max = 1;
  min = 0;
  while (-min < ((1E10) * max))
  {
    for (i = 1; i < G->N; i++)
    {
      Gtemp[i] = G->Im[i];
      if (G->Im[i] > 0)
      {
        Gtemp[i - 1] += G->Im[i] / 2.0;
        Gtemp[i + 1] += G->Im[i] / 2.0;
        Gtemp[i] -= G->Im[i];
      }
    }
    i = 0;
    while (Gtemp[i] > 0)
    {
      Gtemp[i + 1] += Gtemp[i];
      Gtemp[i] = 0;
      i++;
    }
    i = G->N;
    while (Gtemp[i] > 0)
    {
      Gtemp[i - 1] += Gtemp[i];
      Gtemp[i] = 0;
      i--;
    }
    /* and the same routine again but now back from Gtemp to ImG */
    G->Im[0] = Gtemp[0];
    G->Im[G->N] = Gtemp[G->N];
    for (i = 1; i < G->N; i++)
    {
      G->Im[i] = Gtemp[i];
      if (Gtemp[i] > 0)
      {
        G->Im[i - 1] += Gtemp[i] / 2.0;
        G->Im[i + 1] += Gtemp[i] / 2.0;
        G->Im[i] -= Gtemp[i];
      }
    }
    i = 0;
    while (G->Im[i] > 0)
    {
      G->Im[i + 1] += G->Im[i];
      G->Im[i] = 0;
      i++;
    }
    i = G->N;
    while (G->Im[i] > 0)
    {
      G->Im[i - 1] += G->Im[i];
      G->Im[i] = 0;
      i--;
    }
    min = G->Im[0];
    max = G->Im[0];
    for (i = 1; i <= G->N; i++)
    {
      min = Min(min, G->Im[i]);
      max = Max(max, G->Im[i]);
    }
  }
}


void RemoveNegativePart(double *const G, const uint32_t N)
{
  double Gtemp[N + 1];          /* define a temporary array to store G with negative part 'broadenend' */
  double max, min;
  uint32_t i;
  Gtemp[0] = G[0];
  Gtemp[N] = G[N];
  max = 0;
  min = -1;
  while (-((1E10) * min) > (max))
  {
    for (i = 1; i < N; i++)
    {
      Gtemp[i] = G[i];
      if (G[i] < 0)
      {
        Gtemp[i - 1] += G[i] / 2.0;
        Gtemp[i + 1] += G[i] / 2.0;
        Gtemp[i] -= G[i];
      }
    }
    i = 0;
    while (Gtemp[i] < 0)
    {
      Gtemp[i + 1] += Gtemp[i];
      Gtemp[i] = 0;
      i++;
    }
    i = N;
    while (Gtemp[i] < 0)
    {
      Gtemp[i - 1] += Gtemp[i];
      Gtemp[i] = 0;
      i--;
    }
    /* and the same routine again but now back from Gtemp to ImG */
    G[0] = Gtemp[0];
    G[N] = Gtemp[N];
    for (i = 1; i < N; i++)
    {
      G[i] = Gtemp[i];
      if (Gtemp[i] < 0)
      {
        G[i - 1] += Gtemp[i] / 2.0;
        G[i + 1] += Gtemp[i] / 2.0;
        G[i] -= Gtemp[i];
      }
    }
    i = 0;
    while (G[i] < 0)
    {
      G[i + 1] += G[i];
      G[i] = 0;
      i++;
    }
    i = N;
    while (G[i] < 0)
    {
      G[i - 1] += G[i];
      G[i] = 0;
      i--;
    }
    min = G[0];
    max = G[0];
    for (i = 1; i < N + 1; i++)
    {
      min = Min(min, G[i]);
      max = Max(max, G[i]);
    }
  }
}



/*
*   All calculations are done away from the real axes by a distance Gamma/2
*   This leads to a Lorenzian Broadening, that is not physical
*   This is not so bad as Gamma can be choosen small
*   However, for small gaps it is sometimes not so nice as a Lorenzian has very wide tails
*   For DMFT it becomes ever very troublesum as each time the Bath Hamiltonian is fitted to the Greens
*   function an additional broadening is added. This leads to incorrect tails of the spectrum and even worse
*   1/G-w-iGamma/2 does not go to zero, but infinite.
*   Most of the problems are solved if the Lorenzian broadening is replaced by a Gaussian broadening
*   tails become short lived and the small deviation (now a Gaussian Broadening) is not so important
*   and can be made as small as desired (not really as one needs more and more accuracy if Gamma becomes smaller and smaller)
*
*   This function works on the Im part of a spectral function only. The real part should (if needed) be reconstructed by a KK transform
*        This is important as acting with this function on both the Re and Im part seperatly does not !!!!!! retrive a KK consistent
*        Responce function, which therefor could have very funny properties.
*
*   input:
*     Emin      : Minimum energy
*     Emax      : Maximum energy
*     N         : Number of energy points - 1 (DE=(Emax-Emin)/N)
*     Gamma     : Lorentzian broadening FWHM
*     f         : spectral function with N+1 elements (Only real or Im part)
*
*   output:
*     spec      : spectral function with N+1 elements (Only real or Im part)
*/
void ReplaceLorentzianByGaussian(const double Emin, const double Emax, const uint32_t N, const double Gamma,
                                const double *const f, double *const spec)
{
  double DE, kmin, kmax, Dk, *FspecRe, *FspecIm, k, w;
  int    ik, Nk, iw;
  DE = (Emax - Emin) / ((double)N);
  kmin = -3.0 / Gamma;
  kmax = 3.0 / Gamma;
  Nk = (uint32_t) (18 * (Emax - Emin) / Gamma);
  Dk = (kmax - kmin) / ((double)Nk);
  FspecRe = malloc((Nk + 1) * sizeof(double));
  FspecIm = malloc((Nk + 1) * sizeof(double));
#pragma omp parallel for private(k,iw,w)
  for (ik = 0; ik <= Nk; ik++)
  {
    k = ((double)ik * Dk) + kmin;
    FspecRe[ik] = 0.0;
    FspecIm[ik] = 0.0;
    for (iw = 0; iw <= N; iw++)
    {
      w = ((double)iw * DE) + Emin;
      FspecRe[ik] -= f[iw] * sin(2 * Pi * k * w);
      FspecIm[ik] += f[iw] * cos(2 * Pi * k * w);
    }
    FspecRe[ik] *= DE * exp(Pi * Gamma * fabs(k) - Sqr(k * Pi * Gamma) / Log16);
    FspecIm[ik] *= DE * exp(Pi * Gamma * fabs(k) - Sqr(k * Pi * Gamma) / Log16);
  }
#pragma omp parallel for private(w,k,ik)
  for (iw = 0; iw <= N; iw++)
  {
    w = ((double)iw * DE) + Emin;
    spec[iw] = 0.0;
    for (ik = 0; ik <= Nk; ik++)
    {
      k = ((double)ik * Dk) + kmin;
      spec[iw] += Dk * (-FspecRe[ik] * sin(2 * Pi * k * w) + FspecIm[ik] * cos(2 * Pi * k * w));
    }
  }
  free(FspecRe);
  free(FspecIm);
}

/*
*   input  : G.N
*          : G.Emin
*          : G.dE
*          : G.Gamma
*          : G.ZerothMoment
*          : G.FirstMoment
*   output : G.Re
*          : G.Im
*/
void Lorentzian(SpectraType * const G)
{
  double w, denom;
  uint32_t i;
  for (i = 0; i <= G->N; i++)
  {
    w = ((double)i * G->dE) + G->Emin - G->FirstMoment;
    denom = Pi / G->ZerothMoment * (Sqr(G->Gamma / 2.0) + Sqr(w));
    G->Re[i] = (w) / denom;
    G->Im[i] = (-G->Gamma / 2.0) / denom;
  }
}

/*
*   input  : G.N
*          : G.dE
*          : G.Re
*          : G.Im
*          : FWHM
*   output : G.Re
*          : G.Im
*/
void BroadenG(SpectraType * const G, const double FWHM)
{
  int    i, Ein, Eout;
  double dataoutRe[G->N + 1];
  double dataoutIm[G->N + 1];
  const double sa = 2.0 * Sqr(FWHM / (2.0 * sqrt(2.0 * log(2.0))));
  const double sb = (FWHM / (2.0 * sqrt(2.0 * log(2.0)))) * sqrt(2.0 * Pi) / G->dE;
  for (i = 0; i <= (int)G->N; i++)
  {
    dataoutRe[i] = 0;
    dataoutIm[i] = 0;
  }
#pragma omp parallel for private(Ein)
  for (Eout = 0; Eout <= (int)G->N; Eout++)
  {
    for (Ein = 0; Ein <= (int)G->N; Ein++)
    {
      dataoutRe[Eout] += (G->Re[Ein] / sb) * exp(-Sqr((Eout - Ein) * G->dE) / sa);
      dataoutIm[Eout] += (G->Im[Ein] / sb) * exp(-Sqr((Eout - Ein) * G->dE) / sa);
    }
  }
  for (i = 0; i <= (int)G->N; i++)
  {
    G->Re[i] = dataoutRe[i];
    G->Im[i] = dataoutIm[i];
  }
}


/*
*   input  : G.N
*          : G.dE
*          : G.Re
*          : G.Im
*          : FWHM
*   output : G.Re
*          : G.Im
*/
void BroadenGvar(SpectraType * const G, const double *const FWHM)
{
  int    i, Ein, Eout;
  double dataoutRe[G->N + 1];
  double dataoutIm[G->N + 1];
  double sa[G->N + 1];
  double sb[G->N + 1];
  for (i = 0; i <= (int)G->N; i++)
  {
    dataoutRe[i] = 0;
    dataoutIm[i] = 0;
    sa[i] = 2.0 * Sqr(FWHM[i] / (2.0 * sqrt(2.0 * log(2.0))));
    sb[i] = (FWHM[i] / (2.0 * sqrt(2.0 * log(2.0)))) * sqrt(2.0 * Pi) / G->dE;
  }
#pragma omp parallel for private(Ein)
  for (Eout = 0; Eout <= (int)G->N; Eout++)
  {
    for (Ein = 0; Ein <= (int)G->N; Ein++)
    {
      dataoutRe[Eout] += (G->Re[Ein] / sb[Ein]) * exp(-Sqr((Eout - Ein) * G->dE) / sa[Ein]);
      dataoutIm[Eout] += (G->Im[Ein] / sb[Ein]) * exp(-Sqr((Eout - Ein) * G->dE) / sa[Ein]);
    }
  }
  for (i = 0; i <= (int)G->N; i++)
  {
    G->Re[i] = dataoutRe[i];
    G->Im[i] = dataoutIm[i];
  }
}

/*
*   input  : G.N
*          : G.Im
*   output : G.Re
*/
void KramersKronigImToRe(SpectraType * const G)
{
  double ReFFTImG[3 * G->N + 3], ImFFTImG[3 * G->N + 3];
  int32_t i, k;
#pragma omp parallel for private(i)
  for (k = 0; k < 3 * G->N + 3; k++)
  {
    ReFFTImG[k] = 0;
    ImFFTImG[k] = 0;
    for (i = 0; i <= G->N; i++)
    {
      ReFFTImG[k] += G->Im[i] * cos((double)(2 * Pi * k * (1 + G->N + i)) / ((double)(3 * G->N + 3)));
      ImFFTImG[k] += G->Im[i] * sin((double)(2 * Pi * k * (1 + G->N + i)) / ((double)(3 * G->N + 3)));
    }
  }
  for (k = (((3 * G->N + 2) / 2) + 1); k < 3 * G->N + 3; k++)
  {
    ReFFTImG[k] = -ReFFTImG[k];
    ImFFTImG[k] = -ImFFTImG[k];
  }
#pragma omp parallel for private(k)
  for (i = 0; i <= G->N; i++)
  {
    G->Re[i] = 0;
    for (k = 0; k < 3 * G->N + 3; k++)
    {
      G->Re[i] -= ReFFTImG[k] * sin((double)(2 * Pi * k * (1 + G->N + i)) / ((double)(3 * G->N + 3)));
      G->Re[i] += ImFFTImG[k] * cos((double)(2 * Pi * k * (1 + G->N + i)) / ((double)(3 * G->N + 3)));
    }
    G->Re[i] /= (double)(3 * G->N + 3);
  }
}

/*
*   input  : G.N
*          : G.dE
*          : G.Im
*   output : G.ZerothMoment
*/
void CalculateZerothMoment(SpectraType * const G)
{
  G->ZerothMoment = VectorCalculateIntegral(G->Im, G->dE, G->N);
}

/*
*   input  : G.N
*          : G.Emin
*          : G.Emax
*          : G.dE
*          : G.Im
*   output : G.FirstMoment
*/
void CalculateFirstMoment(SpectraType * const G)
{
  double Gbr[G->N + 1], maxspline;
  int    i, Ein, Eout, maxi;
  const double FWHM = (G->EmaxPole - G->EminPole) * 10;
  const double sa = 2.0 * Sqr(FWHM / (2.0 * sqrt(2.0 * log(2.0))));
  const double sb = (FWHM / (2.0 * sqrt(2.0 * log(2.0)))) * sqrt(2.0 * Pi) / G->dE;
  for (i = 0; i <= (int)G->N; i++)
  {
    Gbr[i] = 0;
  }
#pragma omp parallel for private(Ein)
  for (Eout = 0; Eout <= (int)G->N; Eout++)
  {
    for (Ein = 0; Ein <= (int)G->N; Ein++)
    {
      Gbr[Eout] += (G->Im[Ein] / sb) * exp(-Sqr((Eout - Ein) * G->dE) / sa);
    }
  }
  maxi = 0;
  for (i = 0; i <= G->N; i++)
  {
    if (Gbr[i] < Gbr[maxi])
    {
      maxi = i;
    }
  }
  maxspline = maxi + (Gbr[maxi - 1] - Gbr[maxi + 1]) / (2 * Gbr[maxi - 1] - 4 * Gbr[maxi] + 2 * Gbr[maxi + 1]);
  G->FirstMoment = G->Emin + G->dE * maxspline;
  if ((G->FirstMoment > G->EmaxPole) || (G->FirstMoment < G->EminPole))
  {
    G->FirstMoment = (G->EmaxPole + G->EminPole) / 2.0;
    printf("WARNING\nFirst moment could not be determined; placed it inbetween EmaxPole and EminPole\n");
  }
}

/*
*   Calculates the first and the second moment of the Im part of a spectral function
*   Input  : G.Im
*            G.N
*            G.Emin
*            G.Emax
*            G.dE
*            G.EminPole
*            G.EmaxPole
*   Output : G.ZerothMoment
*            G.FirstMoment
*
*   In order to calculate the zeroth and first moment it is assumed that the spectral function
*   behaves like a single Lorenzian outside the energy range given by G.Emin to G.Emax
*/
void CalculateMoment(SpectraType * const G)
{
  /* Calculate the first moment by broadening the spectrum and find the maximum */
  double Gbr[G->N + 1], maxspline, En, Enmax;
  int    i, Ein, Eout, maxi;
  const double FWHM = (G->EmaxPole - G->EminPole) * 10;
  const double sa = 2.0 * Sqr(FWHM / (2.0 * sqrt(2.0 * log(2.0))));
  double sb = (FWHM / (2.0 * sqrt(2.0 * log(2.0)))) * sqrt(2.0 * Pi) / G->dE;
  for (i = 0; i <= (int)G->N; i++)
  {
    Gbr[i] = 0;
  }
#pragma omp parallel for private(Ein)
  for (Eout = 0; Eout <= (int)G->N; Eout++) /* broaden */
  {
    for (Ein = 0; Ein <= (int)G->N; Ein++)
    {
      Gbr[Eout] += (G->Im[Ein] / sb) * exp(-Sqr((Eout - Ein) * G->dE) / sa);
    }
  }
  maxi = 0;                     /* find maximum point */
  for (i = 0; i <= G->N; i++)
  {
    if (Gbr[i] < Gbr[maxi])
    {
      maxi = i;
    }
  }
  maxspline = maxi + (Gbr[maxi - 1] - Gbr[maxi + 1]) / (2 * Gbr[maxi - 1] - 4 * Gbr[maxi] + 2 * Gbr[maxi + 1]); /* spline interpolate between the highest point and its nearest neighbors */
  G->FirstMoment = G->Emin + G->dE * maxspline;
  if ((G->FirstMoment > G->EmaxPole) || (G->FirstMoment < G->EminPole))
  {
    printf("WARNING\nFirst moment outside range EmaxPole and EminPole\n");
    printf("EminPole       =%22.15E\n", G->EminPole);
    printf("EmaxPole       =%22.15E\n", G->EmaxPole);
    printf("FirstMoment old=%22.15E\n", G->FirstMoment);
    G->FirstMoment = (G->EmaxPole + G->EminPole) / 2.0;
    printf("FirstMoment new=%22.15E\n", G->FirstMoment);
  }
  G->ZerothMoment = VectorCalculateIntegral(G->Im, G->dE, G->N);  /* calculate the zeroth moment by a spline integration */
  G->ZerothMoment *= Pi / (atan(2 * (G->FirstMoment - G->Emin) / G->Gamma) + atan(2 * (G->Emax - G->FirstMoment) / G->Gamma));  /* add missing intensity from the wings outside the range G.Emin to G.Emax */
  /* if G.FirstMoment and (G.Emin+G.Emax)/2 are not the same an error has been made in establishing the first moment */
  if (G->FirstMoment > (G->Emin + G->Emax) / 2.0)
  {
    En = G->Emax + G->dE;
    Enmax = 2 * G->FirstMoment - G->Emin;
    sb *= 2 * Pi / (G->Gamma * G->ZerothMoment);
    while (En < Enmax)
    {
      for (Eout = 0; Eout <= G->N; Eout++)
      {
        Gbr[Eout] +=
          (1 / (sb * (Sqr(En - G->FirstMoment) + Sqr(G->Gamma) / 4.0))) * exp(-Sqr((Eout) * G->dE + G->Emin - En) / sa);
      }
      En += G->dE;
    }
  }
  if (G->FirstMoment < (G->Emin + G->Emax) / 2.0)
  {
    En = G->Emin - G->dE;
    Enmax = 2 * G->FirstMoment - G->Emax;
    sb *= 2 * Pi / (G->Gamma * G->ZerothMoment);
    while (En > Enmax)
    {
      for (Eout = 0; Eout <= G->N; Eout++)
      {
        Gbr[Eout] +=
          (1 / (sb * (Sqr(En - G->FirstMoment) + Sqr(G->Gamma) / 4.0))) * exp(-Sqr((Eout) * G->dE + G->Emin - En) / sa);
      }
      En -= G->dE;
    }
  }
  maxi = 0;                     /* find maximum point */
  for (i = 0; i <= G->N; i++)
  {
    if (Gbr[i] < Gbr[maxi])
    {
      maxi = i;
    }
  }
  maxspline = maxi + (Gbr[maxi - 1] - Gbr[maxi + 1]) / (2 * Gbr[maxi - 1] - 4 * Gbr[maxi] + 2 * Gbr[maxi + 1]); /* spline interpolate between the highest point and its nearest neighbors */
  G->FirstMoment = G->Emin + G->dE * maxspline;
  if ((G->FirstMoment > G->EmaxPole) || (G->FirstMoment < G->EminPole))
  {
    G->FirstMoment = (G->EmaxPole + G->EminPole) / 2.0;
    printf("WARNING\n");
  }
}

/*
*   input  : G.N
*          : G.Emin
*          : G.Emax
*          : G.dE
*          : G.EminPole
*          : G.EmaxPole
*          : G.ZerothMoment
*          : G.FirstMoment
*          : G.Gamma
*          : G.Im
*   output : G.Im
*/
void ReplaceImWings(SpectraType * const G)
{
  double Energy;
  int    i;
  for (i = 0; i <= G->N; i++)
  {
    Energy = G->Emin + i * G->dE;
    if (Energy < G->EminPole)
    {
      G->Im[i] =
        G->Im[i] * (Energy - G->Emin) / (G->EminPole - G->Emin) +
        ((G->ZerothMoment * G->Gamma / 2) / (Pi * (Sqr(Energy - G->FirstMoment) + Sqr(G->Gamma / 2)))) * (G->EminPole -
                                                                                                          Energy) /
        (G->EminPole - G->Emin);
    }
    else if (Energy > G->EmaxPole)
    {
      G->Im[i] =
        G->Im[i] * (G->Emax - Energy) / (G->Emax - G->EmaxPole) +
        ((G->ZerothMoment * G->Gamma / 2) / (Pi * (Sqr(Energy - G->FirstMoment) + Sqr(G->Gamma / 2)))) * (Energy -
                                                                                                          G->EmaxPole) /
        (G->Emax - G->EmaxPole);
    }
  }
}



/*
*  Make Im negative, replace the wings by a single lorenzian and recalcualte the real part by a Kramers-Kronig relation
*/
void MakeGPhysical(SpectraType * const G)
{
  RemovePositiveImPart(G);
  ReplaceImWings(G);
  KramersKronigImToRe(G);
}

/*
*  Remove positive Im part, recalcualte the real part by a Kramers-Kronig relation and reset the constant Re part back
*/
void MakeSigmaPhysical(SpectraType * const G)
{
  double GReTot = 0;
  uint32_t i;
  RemovePositiveImPart(G);
  ReplaceImWings(G);
  for (i = 0; i <= G->N; i++)
  {
    GReTot += G->Re[i];
  }
  KramersKronigImToRe(G);
  for (i = 0; i <= G->N; i++)
  {
    GReTot -= G->Re[i];
  }
  GReTot /= (G->N + 1);
  for (i = 0; i <= G->N; i++)
  {
    G->Re[i] += GReTot;
  }
  printf("Constant Re ofset of Sigma %22.15lf\n", GReTot);
}


/*
*   Given a tridiagonal matrix
*           a1 b1
*           b1 a2 b2
*              b2 a3 b3
*                 b3 a4 ..
*                    .. ..
*   Calculates the spectrum of psi0
*     the norm of T psi = psi0 is stored in b0 and the spectrum is multiplied by b0^2
*     an additional shift of the real part is stored in a0
*
*   input:
*     Energy    : Energy of the initial state (shift of spectrum)
*     G.Emin    : Minimum energy
*     G.Emax    : Maximum energy
*     G.N       : Number of energy points - 1 (DE=(Emax-Emin)/N)
*     G.Gamma   : Lorentzian broadening FWHM
*     M         : M.n  : number of points in Matrix
*                 M.a  : diagonal values
*                 M.b  : of diagonal values (M.b[0] contains overall multiplication)
*
*   output      : G.Re[0]...G.Re[N] real elements of spectrum 
*               : G.Im[0]...G.Im[N] Im part of spectrum
*               : G.FirstMoment  \int_{-\infty}^{\infty} G.Im
*               : G.SecondMoment \int_{-\infty}^{\infty} w G.Im / G.FirstMoment
*
*    a0 + b0^2/(w-a1-b1^2/(w-a2-b2^2/(...))
*
*/

void TriDiagonalMatrixToSpectrum(const double Energy, const TriDiagonalMatrixType M, SpectraType * const G)
{
  double hnu, valRe, valIm, denom;
  uint32_t i, j;
  if (fabs(G->dE) > (2 * G->Gamma))
    printf("\n\nNeed more then two points per lorenzian FWHM in order to represent a spectrum\n\n");
  for (i = 0; i <= G->N; i++)
  {
    hnu = ((double)i * G->dE) + G->Emin + Energy;
    denom = (Sqr(G->Gamma * 0.5) + Sqr(hnu - M.a[M.N]));
    valRe = (hnu - M.a[M.N]) / denom;
    valIm = (-G->Gamma * 0.5) / denom;
    for (j = Max(1, M.N) - 1; j >= 1; j--)
    {
      denom = Sqr(G->Gamma / 2.0 - valIm * Sqr(M.b[j])) + Sqr(hnu - M.a[j] - valRe * Sqr(M.b[j]));
      valRe = (hnu - M.a[j] - valRe * Sqr(M.b[j])) / denom;
      valIm = (valIm * Sqr(M.b[j]) - G->Gamma / 2.0) / denom;
    }
    G->Re[i] = M.a[0] + valRe * Sqr(M.b[0]);
    G->Im[i] = valIm * Sqr(M.b[0]);
  }
  G->ZerothMoment = -Sqr(M.b[0]);
  G->FirstMoment = M.a[0] - Energy;
  for (i = 0; i < Min(M.N, 5); i++)
  {
    G->a[i] = M.a[i] - Energy;
    G->b[i] = M.b[i];
  }
  for (i = M.N; i < 5; i++)
  {
    G->a[i] = M.a[M.N] - Energy;
    G->b[i] = 0;
  }
}



/*  the spectrum returns the function
 *  ( w + E0 - M - sum_j Mw_j * gw_j(w) + i Gamma/2 )^{-1}_[1,1]
 *  it is the users obbligation to define the matrices such that the inversion is not singular
 *  The functions gw must be defined on the same grid (emin, emax, dE, N) as the function G
 */
uint32_t MatrixToSpectrum(const MatrixType M, const double E0, const double Gamma, const MatrixType * Mw,
                          const SpectraType * Gw, uint32_t NGw, SpectraType * const G)
{
  printf("Start of MatrixToSpectrum\n");

  MatrixType Mtot, Mtest, MtestI;
  uint32_t i, j;
  double w;
  for (i = 0; i <= G->N; i++)
  {
    w = G->Emin + i * G->dE;
    Mtot.QComplex = 1;
    Mtot.Ni = M.Ni;
    Mtot.Nj = M.Nj;
    InitIdentityMatrix(&Mtot);
    MatrixScalarMultiply(&Mtot, w + E0, Gamma / 2.0);
    MatrixAdd(&Mtot, M, -1, 0);
    for (j = 0; j < NGw; j++)
    {
      MatrixAdd(&Mtot, Mw[j], -Gw[j].Re[i], -Gw[j].Im[i]);
    }

    if (i == 2000)
    {
      MatrixCopy(&Mtest, Mtot);
    }

    printf("%lf %lf %lf %lf %lf | ", w, Mtot.VRe[0][0], Mtot.VIm[0][0], Gw[0].Re[i], Gw[0].Im[i]);
    MatrixInverse(&Mtot);
    printf("%lf %lf \n", Mtot.VRe[0][0], Mtot.VIm[0][0]);
    G->Re[i] = Mtot.VRe[0][0];
    G->Im[i] = Mtot.VIm[0][0];

    if (i == 2000)
    {
      PrintMatrix(Mtest);
      PrintMatrix(Mtot);
      MatrixMultiply(&MtestI, Mtest, Mtot);
      PrintMatrix(MtestI);
    }

    FreeMatrix(&Mtot);
  }
  return 0;
}

/*
 *
 */
uint32_t ListOfPolesToSpectrum(const double Energy, const ListOfPoles P, SpectraType * const G)
{
  TriDiagonalMatrixType PTri;
  ListOfPolesToTridiagonalMatrix(P, &PTri, P.N);
  TriDiagonalMatrixToSpectrum(Energy, PTri, G);
  FreeTriDiagonalMatrix(&PTri);
  return 0;
}

/*
 *
 */
uint32_t ListOfPolesToSpectrumGaussian(const double Energy, const ListOfPoles P, SpectraType * const G, double FWHM)
{
  int    iOut, iIn, i;
  const double sa = 2.0 * Sqr(FWHM / (2.0 * sqrt(2.0 * log(2.0))));
  const double sb = ((2.0 * sqrt(Pi * log(2.0))) / FWHM);
  double Eout;
  for (i = 0; i <= G->N; i++)
  {
    G->Im[i] = 0;
    G->Re[i] = 0;
  }
#pragma omp parallel for private(Eout,iIn)
  for (iOut = 0; iOut <= G->N; iOut++)
  {
    Eout = G->Emin + iOut * G->dE;
    for (iIn = 0; iIn < P.N; iIn++)
    {
      G->Im[iOut] -= (P.b[iIn] * sb) * exp(-Sqr(Eout - P.a[iIn]) / sa);
    }
  }
  KramersKronigImToRe(G);
  return 0;
}


uint32_t TriDiagonalMatrixToSpectrumGaussian(const double Energy, const TriDiagonalMatrixType M, SpectraType * const G,
                                            double FWHM)
{
  ListOfPoles P;
  TriDiagonalMatrixToListOfPoles(M, &P);
  ListOfPolesToSpectrumGaussian(Energy, P, G, FWHM);
  FreeListOfPoles(&P);
  return 0;
}

/*
*   Given a tridiagonal matrix for conduction and valence band seperate
*           avc0 bc0                bv0
*           bc0  ac1 bc1
*                bc1 ac2 bc2
*                    bc2 ac3 ...
*                        ... ...
*  
*           bv0                     av1 bv1
*                                   bv1 av2 bv2
*                                       bv2 av3 ...
*                                           ... ...
*   Calculates the spectrum of psi0
*     the norm of T psi = psi0 is stored in b0 and the spectrum is multiplied by b0^2
*     An additional shift of the real part of the spectra is stored in a0
*
*   input:
*     Energy    : Energy of the initial state (shift of spectrum)
*     G.Emin    : Minimum energy
*     G.Emax    : Maximum energy
*     G.N       : Number of energy points - 1 (DE=(Emax-Emin)/N)
*     Gamma     : Lorentzian broadening FWHM
*     V         : V.n  : number of points in Matrix
*                 V.a  : diagonal values
*                 V.b  : of diagonal values (V.b[0] cotains overall multiplication)
*     C         : C.n  : number of points in Matrix
*                 C.a  : diagonal values
*                 C.b  : of diagonal values (C.b[0] cotains overall multiplication)
*
*   output      : G.Im
*               : G.Re
*
*   a0 + b0^2 / (w-avc0-bc0^2/(w-ac1-bc1^2/(w-ac2-bc2^2/(...))-bv0^2/(w-av1-bv1^2/(w-av2-bv2^2/(...)))
*
*/
void DoubleTriDiagonalMatrixToSpectrum(const double Energy, const DoubleTriDiagonalMatrixType M, SpectraType * const G)
{
  SpectraType GV, GC;
  double hnu, denom;
  uint32_t i;
  SpectrumAllocate(*G, &GV);
  SpectrumAllocate(*G, &GC);
  TriDiagonalMatrixToSpectrum(Energy, M.val, &GV);
  TriDiagonalMatrixToSpectrum(Energy, M.con, &GC);
  for (i = 0; i <= G->N; i++)
  {
    hnu = ((double)i * G->dE) + G->Emin + Energy;
    denom = Sqr(G->Gamma * 0.5 - GV.Im[i] - GC.Im[i]) + Sqr(hnu - (M.val.a[0] + M.con.a[0]) / 2 - GV.Re[i] - GC.Re[i]);
    G->Re[i] = M.a0 + Sqr(M.b0) * (hnu - (M.val.a[0] + M.con.a[0]) / 2 - GV.Re[i] - GC.Re[i]) / denom;
    G->Im[i] = Sqr(M.b0) * (-G->Gamma * 0.5 + GV.Im[i] + GC.Im[i]) / denom;
  }
  SpectrumFree(&GV);
  SpectrumFree(&GC);
}

void SecondDoubleTriDiagonalMatrixToSpectrum(const double Energy, const SecondDoubleTriDiagonalMatrixType M,
                                             SpectraType * const G)
{
  SpectraType GSec;
  double hnu, denom;
  uint32_t i;
  SpectrumAllocate(*G, &GSec);
  DoubleTriDiagonalMatrixToSpectrum(Energy, M.Dt, &GSec);
  for (i = 0; i <= G->N; i++)
  {
    hnu = ((double)i * G->dE) + G->Emin + Energy;
    denom = Sqr(G->Gamma * 0.5 - GSec.Im[i]) + Sqr(hnu - GSec.Re[i]);
    G->Re[i] = M.a0 + Sqr(M.b0) * (hnu - GSec.Re[i]) / denom;
    G->Im[i] = Sqr(M.b0) * (-G->Gamma * 0.5 + GSec.Im[i]) / denom;
  }
  SpectrumFree(&GSec);
}

void BADoubleTriDiagonalMatrixToSpectrum(const double Energy, const BADoubleTriDiagonalMatrixType M,
                                         SpectraType * const G)
{
  SpectraType GSec;
  TriDiagonalMatrixToSpectrum(Energy, M.con, G);
  if (SpectrumAllocate(*G, &GSec))
  {
    printf("SpectrumAllocate failed in BADoubleTriDiagonalMatrixToSpectrum;");
    return;
  }
  TriDiagonalMatrixToSpectrum(Energy, M.val, &GSec);
  AddMultipliedSpectra(G, GSec, 1, 0);
  SpectrumFree(&GSec);
}

void DoubleTriDiagonalMatrixToSpectrumOld(const double Energy, const TriDiagonalMatrixType V,
                                          const TriDiagonalMatrixType C, SpectraType * const G)
{
  double hnu, valReC, valImC, valReV, valImV, denom;
  uint32_t i, j;
  if (fabs(G->dE) > (2 * G->Gamma))
    printf("\n\nNeed more then two points per lorenzian FWHM in order to represent a spectrum\n\n");
  for (i = 0; i <= G->N; i++)
  {
    hnu = ((double)i * G->dE) + G->Emin + Energy;
    denom = Sqr(G->Gamma / 2.0) + Sqr(hnu - V.a[V.N]);
    valReV = (hnu - V.a[V.N]) / denom;
    valImV = (-G->Gamma / 2.0) / denom;
    for (j = V.N; j != 1;)
    {
      j--;
      denom = Sqr(G->Gamma / 2.0 - valImV * Sqr(V.b[j + 1])) + Sqr(hnu - V.a[j] - valReV * Sqr(V.b[j + 1]));
      valReV = (hnu - V.a[j] - valReV * Sqr(V.b[j + 1])) / denom;
      valImV = (valImV * Sqr(V.b[j + 1]) - G->Gamma / 2.0) / denom;
    }

    denom = Sqr(G->Gamma / 2.0) + Sqr(hnu - C.a[C.N]);
    valReC = (hnu - C.a[C.N]) / denom;
    valImC = (-G->Gamma / 2.0) / denom;
    for (j = C.N; j != 1;)
    {
      j--;
      denom = Sqr(G->Gamma / 2.0 - valImC * Sqr(C.b[j + 1])) + Sqr(hnu - C.a[j] - valReC * Sqr(C.b[j + 1]));
      valReC = (hnu - C.a[j] - valReC * Sqr(C.b[j + 1])) / denom;
      valImC = (valImC * Sqr(C.b[j + 1]) - G->Gamma / 2.0) / denom;
    }

    denom =
      (Sqr(hnu - V.a[0] - Sqr(C.b[1]) * valReC - Sqr(V.b[1]) * valReV) +
       Sqr(G->Gamma / 2.0 - Sqr(C.b[1]) * valImC - Sqr(V.b[1]) * valImV));
    G->Re[i] = Sqr(V.b[0]) * (hnu - V.a[0] - Sqr(C.b[1]) * valReC - Sqr(V.b[1]) * valReV) / denom;
    G->Im[i] = Sqr(V.b[0]) * (Sqr(C.b[1]) * valImC + Sqr(V.b[1]) * valImV - (G->Gamma / 2.0)) / denom;
  }
  G->ZerothMoment = -(Sqr(V.b[0]) + Sqr(C.b[0])) / 2.0;
  G->FirstMoment = (V.a[0] + C.a[0]) / 2.0;
}


/*
*   Given an Anderson impurity matrix calculate the spectral function
*          a1 b1 b2 b3 ..
*          b1 a2
*          b2    a3
*          b3       a4
*          ..          ..
*
*   Calculates the spectrum of psi0
*     the norm of T psi = psi0 is stored in b0 and the spectrum is multiplied by b0^2
*     a0 contains an additional shift of the real part
*
*   input:
*     Energy    : Energy of the initial state (shift of spectrum)
*     G.Emin    : Minimum energy
*     G.Emax    : Maximum energy
*     G.N       : Number of energy points - 1 (DE=(Emax-Emin)/N)
*     G.Gamma   : Lorentzian broadening FWHM
*     M         : M.n  : number of points in Matrix
*                 M.a  : diagonal values
*                 M.b  : of diagonal values (M.b[0] cotains overall multiplication)
*
*   output      : G.Re
*               : G.Im
*               : G.ZerothMoment
*               : G.FirstMoment
*/
void AndersonImpurityMatrixToSpectrum(const double Energy, const AndersonMatrixType M, SpectraType * const G)
{
  double hnu, valRe, valIm, denom;
  uint32_t i, j;
  if (fabs(G->dE) > (2 * G->Gamma))
    printf("\n\nNeed more then two points per lorenzian FWHM in order to represent a spectrum\n\n");
  for (i = 0; i <= G->N; i++)
  {
    hnu = ((double)i * G->dE) + G->Emin + Energy;
    valRe = 0.0;
    valIm = 0.0;
    for (j = M.N; j >= 2; j--)
    {
      denom = Sqr(G->Gamma / 2.0) + Sqr(hnu - M.a[j]);
      valRe += Sqr(M.b[j - 1]) * (hnu - M.a[j]) / denom;
      valIm += Sqr(M.b[j - 1]) * (-G->Gamma / 2.0) / denom;
    }
    denom = (Sqr(G->Gamma / 2.0 - valIm) + Sqr(hnu - M.a[1] - valRe));
    G->Re[i] = M.a[0] + (hnu - M.a[1] - valRe) * Sqr(M.b[0]) / denom;
    G->Im[i] = (valIm - G->Gamma / 2.0) * Sqr(M.b[0]) / denom;
  }
  G->ZerothMoment = -Sqr(M.b[0]);
  G->FirstMoment = M.a[1];
}

/*
*   input
*     M        : Anderson Impurity matrix
*                M.a[0] contains first moment of spectrum (untouched by routine)
*                M.b[0] contains overall multiplication (untouched by routine)
*     mu       : chemical potential, determines if point belongs to V or C
*     zero     : Energies smaller then mu contribute to V,
*                Energies larger then mu contribute to C,
*                Energies in a range of mu-zero to mu+zero contribute half to V and half to C
*/
void SplitAndersonImpurityHamiltonian(const TriDiagonalMatrixType M, const double mu, const double zero,
                                      TriDiagonalMatrixType * const MV, TriDiagonalMatrixType * const MC)
{
  uint32_t i;
  MV->N = 1;
  MC->N = 1;
  MV->a[0] = M.a[0];
  MC->a[0] = M.a[0];
  MV->b[0] = M.b[0];
  MC->b[0] = M.b[0];
  for (i = 1; i <= M.N; i++)
  {
    if (M.a[i] < mu - zero)
    {
      MV->a[MV->N] = M.a[i];
      MV->b[MV->N] = M.b[i];
      MV->N++;
    }
    else if (M.a[i] > mu + zero)
    {
      MC->a[MC->N] = M.a[i];
      MC->b[MC->N] = M.b[i];
      MC->N++;
    }
    else
    {
      MV->a[MV->N] = M.a[i];
      MV->b[MV->N] = M.b[i] / sqrt(2.0);
      MV->N++;
      MC->a[MC->N] = M.a[i];
      MC->b[MC->N] = M.b[i] / sqrt(2.0);
      MC->N++;
    }
  }
}


/*
*   These functions must suffer from numberloss !!!!!!!!!!!!!!!!!!!
*   They can and probably should be replaced by functions acting on the fitted spectra !!!!!!!!!!!!!!!
*   S <- 1/G0-1/G
*/
void CalculateSelfEnergySpectralRepresentation(const SpectraType G0, const SpectraType G, SpectraType * const S,
                                               const double *const FWHM)
{
  uint32_t i;
  for (i = 0; i <= G0.N; i++)
  {
    S->Re[i] = (G0.Re[i] / (Sqr(G0.Re[i]) + Sqr(G0.Im[i])) - G.Re[i] / (Sqr(G.Re[i]) + Sqr(G.Im[i])));
    S->Im[i] = (-G0.Im[i] / (Sqr(G0.Re[i]) + Sqr(G0.Im[i])) + G.Im[i] / (Sqr(G.Re[i]) + Sqr(G.Im[i])));
  }
  S->N = G0.N;
  S->Emin = G0.Emin;
  S->Emax = G0.Emax;
  S->dE = G0.dE;
  S->EminPole = Min(G0.EminPole, G.EminPole);
  S->EmaxPole = Max(G0.EmaxPole, G.EmaxPole);
  S->Gamma = G0.Gamma;
  CalculateMoment(S);
  BroadenGvar(S, FWHM);
  RemovePositiveImPart(S);
  ReplaceImWings(S);
  ReplaceImWings(S);
  KramersKronigImToRe(S);
  for (i = 0; i <= S->N; i++)
  {
    S->Re[i] += Pi * (G.FirstMoment - G0.FirstMoment);
  }
  printf("Constant Re ofset of Sigma (units of Pi) %22.15lf\n", G.FirstMoment - G0.FirstMoment);
}

/*
*   These functions must suffer from numberloss !!!!!!!!!!!!!!!!!!!
*   They can and probably should be replaced by functions acting on the fitted spectra !!!!!!!!!!!!!!!
*   G0 <- 1/(S+1/G)
*/
void CalculateG0fromSandGSpectralRepresentation(SpectraType * const G0, const SpectraType G, const SpectraType S)
{
  double denom;
  uint32_t i;
  for (i = 0; i <= G.N; i++)
  {
    denom =
      1 - 2 * G.Im[i] * S.Im[i] + 2 * G.Re[i] * S.Re[i] + Sqr(G.Im[i]) * (Sqr(S.Im[i]) + Sqr(S.Re[i])) +
      Sqr(G.Re[i]) * (Sqr(S.Im[i]) + Sqr(S.Re[i]));
    G0->Re[i] = (G.Re[i] + (Sqr(G.Im[i]) + Sqr(G.Re[i])) * S.Re[i]) / denom;
    G0->Im[i] = (G.Im[i] - (Sqr(G.Im[i]) + Sqr(G.Re[i])) * S.Im[i]) / denom;
  }
  G0->N = G.N;
  G0->Emin = G.Emin;
  G0->Emax = G.Emax;
  G0->dE = G.dE;
  G0->EminPole = Min(S.EminPole, G.EminPole);
  G0->EmaxPole = Max(S.EmaxPole, G.EmaxPole);
  G0->Gamma = G.Gamma;
  CalculateZerothMoment(G0);
  CalculateFirstMoment(G0);
}


/* 
* For k resolved spectroscopy one wants to calculate <i|\Sum_r` e^{- i k . r`} T^\dag_r` G \sum_r e^{i k . r} T_r |i>
* However it is simpler to calculate <i|T_(r)+T^\dag_(r+dr) G T_(r)+T^\dag_(r+dr)|i>
* This function assumes T(dr) as input for the list spectra and replaces this by the fourier transform.
* There are Natom atoms and therefor also Natom spectra assumed the first one belongs to dr=0
* There must thus also be Natom positions listed. The first spectrum is taken at the first position.
* The n-th spectrum is taken as the sum of T_atom_0 +T_atom_n
*/
int SiteFourierTransform(uint32_t Natom, int32_t AtomPos[][3], int32_t UnitCell[3][3], double klist[][3],
                         SpectraType * Spectra)
{
  int32_t i, j, k;
  SpectraType SpectraDeltaR[Natom];
  double AtomConnect[3];
  for (j = 0; j < Natom; j++)
  {
    if (SpectrumAllocate(Spectra[0], &SpectraDeltaR[j]) != 0)
    {
      for (i = 0; i < j; i++)
      {
        SpectrumFree(&SpectraDeltaR[i]);
      }
      printf("Could not allocate memory for internal use in SiteFourierTransform\n function returned 1\n");
      return 1;
    }
    CopySpectra(Spectra[j], &SpectraDeltaR[j]);
    if (j != 0)
    {
      SubtractSpectra(&SpectraDeltaR[j], Spectra[0]);
    }
  }
  for (k = 0; k < Natom; k++)
  {
    SetSpectraToZero(&Spectra[k]);
    for (j = 0; j < Natom; j++)
    {
      AtomConnect[0] = AtomPos[0][0] - AtomPos[j][0];
      AtomConnect[1] = AtomPos[0][1] - AtomPos[j][1];
      AtomConnect[2] = AtomPos[0][2] - AtomPos[j][2];
      AddMultipliedSpectra(&Spectra[k], SpectraDeltaR[j], cos(Dot(klist[k], AtomConnect)),
                           sin(Dot(klist[k], AtomConnect)));
    }
  }
  return 0;
}
