double DOSFunctionBethe(double w)
{
  return (4.0 / Pi) * sqrt(Max(0, 1 - 4 * Sqr(w)));
}

double DOSFunctionBlock(double w)
{
  if ((w > -0.5) && (w < 0.5))
    return 1.0;
  if ((w == 0.5) || (w == -0.5))
    return 0.5;
  return 0;
}

double DOSFunctionChain(double w)
{
  if ((w > -0.5) && (w < 0.5))
    return (2.0 / Pi) * sqrt(1 / (1 - 4 * Sqr(w)));
  if ((w == 0.5) || (w == -0.5))
    return 0.5;
  return 0;
}

/*  given a function DOS(w) a tri-diagonal representation is created
 *
 */
uint32_t DOSFunctionToGTriDiagonalRepresentation(TriDiagonalMatrixType * G, uint32_t N, double (*DOS) (double),
                                                 double Emin, double Emax, double DE)
{
  uint32_t i;
  ListOfPoles P;
  P.N = (uint32_t) (1 + (Emax - Emin) / DE);
  if (InitListOfPoles(&P))
  {
    printf("InitListOfPoles failed in DOSFunctionToGTriDiagonalRepresentation\n");
    fflush(stdout);
    return 1;
  };
  for (i = 0; i < P.N; i++)
  {
    P.a[i] = Emin + i * DE;
    P.b[i] = DE * DOS(Emin + i * DE);
  }
  if (ListOfPolesToTridiagonalMatrix(P, G, N))
  {
    printf("ListOfPolesToTridiagonalMatrix failed in DOSFunctionToGTriDiagonalRepresentation\n");
    fflush(stdout);
    FreeListOfPoles(&P);
    return 1;
  };
  FreeListOfPoles(&P);
  return 0;
}

/*  given a function DOS(w) a tri-diagonal representation is created
 *
 */
uint32_t DOSFunctionToGTriDiagonalRepresentationNonOrthogonal(TriDiagonalMatrixType * G, uint32_t N,
                                                              double (*DOS) (double), double Emin, double Emax,
                                                              double DE)
{
  uint32_t i;
  ListOfPoles P;
  P.N = (uint32_t) (1 + (Emax - Emin) / DE);
  if (InitListOfPoles(&P))
  {
    printf("InitListOfPoles failed in DOSFunctionToGTriDiagonalRepresentation\n");
    fflush(stdout);
    return 1;
  };
  for (i = 0; i < P.N; i++)
  {
    P.a[i] = Emin + i * DE;
    P.b[i] = DE * DOS(Emin + i * DE);
  }
  if (ListOfPolesToTridiagonalMatrixNonOrthogonal(P, G, N))
  {
    printf("ListOfPolesToTridiagonalMatrix failed in DOSFunctionToGTriDiagonalRepresentation\n");
    fflush(stdout);
    FreeListOfPoles(&P);
    return 1;
  };
  FreeListOfPoles(&P);
  return 0;
}

uint32_t NearestNeighborHyperCubicTBToGTriDiagonalRepresentation(uint32_t Dim, uint32_t NTriDiag,
                                                                 TriDiagonalMatrixType * const G, uint32_t Nk)
{
  int    i, j, k;
  const double dk = 2 * Pi / ((double)Nk);
  const double t = 0.5 / ((double)Dim);
  ListOfPoles P;
  P.N = 1;
  for (i = 0; i < Dim; i++)
  {
    P.N *= Nk;
  }
  if (InitListOfPoles(&P))
  {
    printf("InitListOfPoles failed in NearestNeighborHyperCubicTBToGTriDiagonalRepresentation\n");
    fflush(stdout);
    return 1;
  };
  P.b[0] = (double)(1) / ((double)P.N);
#pragma omp for private(j,k)
  for (i = 0; i < P.N; i++)
  {
    P.b[i] = P.b[0];
    P.a[i] = 0;
    k = i;
    for (j = 0; j < Dim; j++)
    {
      P.a[i] +=
        0.5 * (cos((k % Nk) * dk + (j * dk) / (2 * (double)Dim)) + cos((k % Nk) * dk - (j * dk) / (2 * (double)Dim)));
      k /= Nk;
    }
    P.a[i] *= t;
  }
  if (P.N > 1000 * NTriDiag)
  {
    printf("Need to reduce the number of poles from %12i\n", P.N);
    ReduceListOfPolesSpacing(&P, 0.001 / NTriDiag);
    printf("                                   to   %12i\n", P.N);
  }
  if (ListOfPolesToTridiagonalMatrix(P, G, NTriDiag))
  {
    printf("ListOfPolesToTridiagonalMatrix failed in NearestNeighborHyperCubicTBToGTriDiagonalRepresentation\n");
    fflush(stdout);
    FreeListOfPoles(&P);
    return 1;
  };
  FreeListOfPoles(&P);
  return 0;
}


uint32_t NearestNeighborHyperCubicTBToGTriDiagonalRepresentationNonOrthogonal(uint32_t Dim, uint32_t NTriDiag,
                                                                              TriDiagonalMatrixType * const G,
                                                                              uint32_t Nk)
{
  int    i, j, k;
  const double dk = 2 * Pi / ((double)Nk);
  const double t = 0.5 / ((double)Dim);
  ListOfPoles P;
  P.N = 1;
  for (i = 0; i < Dim; i++)
  {
    P.N *= Nk;
  }
  if (InitListOfPoles(&P))
  {
    printf("InitListOfPoles failed in NearestNeighborHyperCubicTBToGTriDiagonalRepresentation\n");
    fflush(stdout);
    return 1;
  };
  P.b[0] = (double)(1) / ((double)P.N);
#pragma omp for private(j,k)
  for (i = 0; i < P.N; i++)
  {
    P.b[i] = P.b[0];
    P.a[i] = 0;
    k = i;
    for (j = 0; j < Dim; j++)
    {
      P.a[i] += cos((k % Nk) * dk);
      k /= Nk;
    }
    P.a[i] *= t;
  }
  if (ListOfPolesToTridiagonalMatrixNonOrthogonal(P, G, NTriDiag))
  {
    printf("ListOfPolesToTridiagonalMatrix failed in NearestNeighborHyperCubicTBToGTriDiagonalRepresentation\n");
    fflush(stdout);
    FreeListOfPoles(&P);
    return 1;
  };
  FreeListOfPoles(&P);
  return 0;
}

uint32_t NearestNeighborTriangularTBToGTriDiagonalRepresentation(uint32_t NTriDiag, TriDiagonalMatrixType * const G,
                                                                 uint32_t Nk)
{
  int    i, kx, ky;
  const double dk = Pi / ((double)Nk);
  ListOfPoles P;
  P.N = Nk * Nk;
  if (InitListOfPoles(&P))
  {
    printf("InitListOfPoles failed in NearestNeighborHyperCubicTBToGTriDiagonalRepresentation\n");
    fflush(stdout);
    return 1;
  };
  P.b[0] = (double)(1) / ((double)P.N);
#pragma omp for private(kx,ky)
  for (i = 0; i < P.N; i++)
  {
    P.b[i] = P.b[0];
    kx = i % Nk;
    ky = i / Nk;
    P.a[i] = (1.0 / 6.0) - (2.0 / 9.0) * (cos(kx * dk) + 2 * cos(kx * dk / 2.0) * cos(ky));
  }
  if (P.N > 10 * NTriDiag)
  {
    printf("Need to reduce the number of poles from %12i\n", P.N);
    ReduceListOfPolesSpacing(&P, 0.1 / NTriDiag);
    printf("                                   to   %12i\n", P.N);
  }
  if (ListOfPolesToTridiagonalMatrix(P, G, NTriDiag))
  {
    printf("ListOfPolesToTridiagonalMatrix failed in NearestNeighborHyperCubicTBToGTriDiagonalRepresentation\n");
    fflush(stdout);
    FreeListOfPoles(&P);
    return 1;
  };
  FreeListOfPoles(&P);
  return 0;
}


uint32_t NearestNeighborTriangularTBToGSpectraRepresentation(SpectraType * const G, uint32_t Nk)
{
  int    i, kx, ky;
  const double dk = Pi / ((double)Nk);
  ListOfPoles P;
  P.N = Nk * Nk;
  if (InitListOfPoles(&P))
  {
    printf("InitListOfPoles failed in NearestNeighborHyperCubicTBToGTriDiagonalRepresentation\n");
    fflush(stdout);
    return 1;
  };
  P.b[0] = (double)(1) / ((double)P.N);
#pragma omp for private(kx,ky)
  for (i = 0; i < P.N; i++)
  {
    P.b[i] = P.b[0];
    kx = i % Nk;
    ky = i / Nk;
    P.a[i] = (1.0 / 6.0) - (2.0 / 9.0) * (cos(kx * dk) + 2 * cos(kx * dk / 2.0) * cos(ky));
  }
  ListOfPolesToSpectrumGausian(0, P, G, G->Gamma);

  FreeListOfPoles(&P);
  return 0;
}
