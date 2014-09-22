/*  These functions transform matrices to operators for some standard models
 *
 */





/*  Creates an anderson impurity Hamiltonian
 *
 *  The impurity is at index i=0 for spin up and i=1 for spin down
 *  The spin up bath is from i=2       to i=Mup.N
 *  The spin dn bath is from i=Mup.N+1 to i=Mup.N+Mdn.N-1
 */
uint32_t TriDiagonalMatrixToAndersonImpurityOperator(OperatorType * O, TriDiagonalMatrixType Mup, TriDiagonalMatrixType Mdn, double U, double mu)
{
  uint32_t i;
  uint16_t oadd[4];
  InitNamedOperator(O, "Single Site Anderson Impurity Hamiltonian");
  O->NFermion = Mup.N + Mdn.N;
  O->NWord = Min(O->NFermion, 1 + (O->NFermion - 1) / 8);
  InitOperatorFixedLength(O, 2, 8);
  InitOperatorFixedLength(O, 4, 8);
  oadd[3] = CrAnMask + 0;
  oadd[2] = CrAnMask + 1;
  oadd[1] = 1;
  oadd[0] = 0;
  RealOperatorAddLaderNormalOrder(O, oadd, 4, U, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 0;
  oadd[0] = 0;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, mu, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 1;
  oadd[0] = 1;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, mu, &QDetNotOrderedFcFaBcBa);
  /* spin up */
  oadd[1] = CrAnMask + 0;
  oadd[0] = 0;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.a[1], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 0;
  oadd[0] = 2;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.b[1], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 2;
  oadd[0] = 0;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.b[1], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 2;
  oadd[0] = 2;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.a[2], &QDetNotOrderedFcFaBcBa);
  for (i = 3; i <= Mup.N; i++)
  {
    oadd[1] = CrAnMask + i - 1;
    oadd[0] = i;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.b[i - 1], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i;
    oadd[0] = i - 1;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.b[i - 1], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i;
    oadd[0] = i;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.a[i], &QDetNotOrderedFcFaBcBa);
  }
  /* spin dn */
  oadd[1] = CrAnMask + 1;
  oadd[0] = 1;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.a[1], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 1;
  oadd[0] = Mup.N + 1;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.b[1], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + Mup.N + 1;
  oadd[0] = 1;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.b[1], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + Mup.N + 1;
  oadd[0] = Mup.N + 1;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.a[2], &QDetNotOrderedFcFaBcBa);
  for (i = 2; i < Mdn.N; i++)
  {
    oadd[1] = CrAnMask + Mup.N + i - 1;
    oadd[0] = Mup.N + i;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.b[i], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + Mup.N + i;
    oadd[0] = Mup.N + i - 1;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.b[i], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + Mup.N + i;
    oadd[0] = Mup.N + i;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.a[i + 1], &QDetNotOrderedFcFaBcBa);
  }
  return 0;
}

/*  Creates an anderson impurity Hamiltonian
 *
 *  The impurity is at index i=0 for spin up and i=1 for spin down
 *  The spin up bath is from i=2       to i=Mup.N
 *  The spin dn bath is from i=Mup.N+1 to i=Mup.N+Mdn.N-1
 */
uint32_t AndersonMatrixToAndersonImpurityOperator(OperatorType * O, AndersonMatrixType Mup, AndersonMatrixType Mdn,
                                                  double U, double mu)
{
  uint32_t i;
  uint16_t oadd[4];
  InitNamedOperator(O, "Single Site Anderson Impurity Hamiltonian");
  O->NFermion = Mup.N + Mdn.N;
  O->NWord = Min(O->NFermion, 1 + (O->NFermion - 1) / 8);
  InitOperatorFixedLength(O, 2, 8);
  InitOperatorFixedLength(O, 4, 8);
  oadd[3] = CrAnMask + 0;
  oadd[2] = CrAnMask + 1;
  oadd[1] = 1;
  oadd[0] = 0;
  RealOperatorAddLaderNormalOrder(O, oadd, 4, U, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 0;
  oadd[0] = 0;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, mu, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 1;
  oadd[0] = 1;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, mu, &QDetNotOrderedFcFaBcBa);
  /* spin up */
  oadd[1] = CrAnMask + 0;
  oadd[0] = 0;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.a[1], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 0;
  oadd[0] = 2;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.b[1], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 2;
  oadd[0] = 0;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.b[1], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 2;
  oadd[0] = 2;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.a[2], &QDetNotOrderedFcFaBcBa);
  for (i = 3; i <= Mup.N; i++)
  {
    oadd[1] = CrAnMask + 0;
    oadd[0] = i;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.b[i - 1], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i;
    oadd[0] = 0;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.b[i - 1], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i;
    oadd[0] = i;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.a[i], &QDetNotOrderedFcFaBcBa);
  }
  /* spin dn */
  oadd[1] = CrAnMask + 1;
  oadd[0] = 1;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.a[1], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 1;
  oadd[0] = Mup.N + 1;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.b[1], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + Mup.N + 1;
  oadd[0] = 1;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.b[1], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + Mup.N + 1;
  oadd[0] = Mup.N + 1;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.a[2], &QDetNotOrderedFcFaBcBa);
  for (i = 2; i < Mdn.N; i++)
  {
    oadd[1] = CrAnMask + 1;
    oadd[0] = Mup.N + i;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.b[i], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + Mup.N + i;
    oadd[0] = 1;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.b[i], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + Mup.N + i;
    oadd[0] = Mup.N + i;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.a[i + 1], &QDetNotOrderedFcFaBcBa);
  }
  return 0;
}


/*  Creates an anderson impurity Hamiltonian
 *
 *  The impurity is at index i=0 for spin up and i=1 for spin down
 *  The first bath site is at index i=2 for spin up and i=3 for spin down
 *  The spin up bath is from i=4                     to i=Mup.val.N+Mup.con.N+3
 *  The spin dn bath is from i=Mup.val.N+Mup.con.N+4 to i=Mup.val.N+Mup.con.N+Mdn.val.N+Mdn.con.N+3
 */
uint32_t SecondDoubleTriDiagonalMatrixToAndersonImpurityOperator(OperatorType * O,
                                                                 SecondDoubleTriDiagonalMatrixType Mup,
                                                                 SecondDoubleTriDiagonalMatrixType Mdn, double U,
                                                                 double mu)
{
  uint32_t i;
  uint16_t oadd[4];
  InitNamedOperator(O, "Single Site Anderson Impurity Hamiltonian");
  O->NFermion = Mup.Dt.val.N + Mup.Dt.con.N + Mdn.Dt.val.N + Mdn.Dt.con.N + 4;
  O->NWord = Min(O->NFermion, 1 + (O->NFermion - 1) / 8);
  InitOperatorFixedLength(O, 2, 8);
  InitOperatorFixedLength(O, 4, 8);
  /* impurity */
  oadd[3] = CrAnMask + 0;
  oadd[2] = CrAnMask + 1;
  oadd[1] = 1;
  oadd[0] = 0;
  RealOperatorAddLaderNormalOrder(O, oadd, 4, U, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 0;
  oadd[0] = 0;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, mu, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 1;
  oadd[0] = 1;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, mu, &QDetNotOrderedFcFaBcBa);
  /* first bath site spin up */
  oadd[1] = CrAnMask + 0;
  oadd[0] = 0;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.Dt.a0, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 0;
  oadd[0] = 2;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.Dt.b0, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 2;
  oadd[0] = 0;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.Dt.b0, &QDetNotOrderedFcFaBcBa);
  /* first bath site spin dn */
  oadd[1] = CrAnMask + 1;
  oadd[0] = 1;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.Dt.a0, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 1;
  oadd[0] = 3;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.Dt.b0, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 3;
  oadd[0] = 1;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.Dt.b0, &QDetNotOrderedFcFaBcBa);
  /* spin up valence */
  oadd[1] = CrAnMask + 2;
  oadd[0] = 2;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.Dt.val.a[0] / 2.0, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 2;
  oadd[0] = 4;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.Dt.val.b[0], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 4;
  oadd[0] = 2;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.Dt.val.b[0], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 4;
  oadd[0] = 4;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.Dt.val.a[1], &QDetNotOrderedFcFaBcBa);
  for (i = 2; i <= Mup.Dt.val.N; i++)
  {
    oadd[1] = CrAnMask + i + 2;
    oadd[0] = i + 3;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.Dt.val.b[i - 1], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i + 3;
    oadd[0] = i + 2;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.Dt.val.b[i - 1], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i + 3;
    oadd[0] = i + 3;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.Dt.val.a[i], &QDetNotOrderedFcFaBcBa);
  }
  /* spin up conductance */
  oadd[1] = CrAnMask + 2;
  oadd[0] = 2;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.Dt.con.a[0] / 2.0, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 2;
  oadd[0] = Mup.Dt.val.N + 4;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.Dt.con.b[0], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + Mup.Dt.val.N + 4;
  oadd[0] = 2;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.Dt.con.b[0], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + Mup.Dt.val.N + 4;
  oadd[0] = Mup.Dt.val.N + 4;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.Dt.con.a[1], &QDetNotOrderedFcFaBcBa);
  for (i = 2; i <= Mup.Dt.con.N; i++)
  {
    oadd[1] = CrAnMask + i + 2 + Mup.Dt.val.N;
    oadd[0] = i + 3 + Mup.Dt.val.N;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.Dt.con.b[i - 1], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i + 3 + Mup.Dt.val.N;
    oadd[0] = i + 2 + Mup.Dt.val.N;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.Dt.con.b[i - 1], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i + 3 + Mup.Dt.val.N;
    oadd[0] = i + 3 + Mup.Dt.val.N;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.Dt.con.a[i], &QDetNotOrderedFcFaBcBa);
  }
  /* spin dn valence */
  oadd[1] = CrAnMask + 3;
  oadd[0] = 3;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.Dt.val.a[0] / 2.0, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 3;
  oadd[0] = Mup.Dt.val.N + Mup.Dt.con.N + 4;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.Dt.val.b[0], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + Mup.Dt.val.N + Mup.Dt.con.N + 4;
  oadd[0] = 3;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.Dt.val.b[0], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + Mup.Dt.val.N + Mup.Dt.con.N + 4;
  oadd[0] = Mup.Dt.val.N + Mup.Dt.con.N + 4;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.Dt.val.a[1], &QDetNotOrderedFcFaBcBa);
  for (i = 2; i <= Mdn.Dt.val.N; i++)
  {
    oadd[1] = CrAnMask + i + Mup.Dt.val.N + Mup.Dt.con.N + 2;
    oadd[0] = i + Mup.Dt.val.N + Mup.Dt.con.N + 3;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.Dt.val.b[i - 1], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i + Mup.Dt.val.N + Mup.Dt.con.N + 3;
    oadd[0] = i + Mup.Dt.val.N + Mup.Dt.con.N + 2;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.Dt.val.b[i - 1], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i + Mup.Dt.val.N + Mup.Dt.con.N + 3;
    oadd[0] = i + Mup.Dt.val.N + Mup.Dt.con.N + 3;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.Dt.val.a[i], &QDetNotOrderedFcFaBcBa);
  }
  /* spin dn conductance */
  oadd[1] = CrAnMask + 3;
  oadd[0] = 3;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.Dt.con.a[0] / 2.0, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 3;
  oadd[0] = Mup.Dt.val.N + Mup.Dt.con.N + Mdn.Dt.val.N + 4;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.Dt.con.b[0], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + Mup.Dt.val.N + Mup.Dt.con.N + Mdn.Dt.val.N + 4;
  oadd[0] = 3;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.Dt.con.b[0], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + Mup.Dt.val.N + Mup.Dt.con.N + Mdn.Dt.val.N + 4;
  oadd[0] = Mup.Dt.val.N + Mup.Dt.con.N + Mdn.Dt.val.N + 4;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.Dt.con.a[1], &QDetNotOrderedFcFaBcBa);
  for (i = 2; i <= Mdn.Dt.con.N; i++)
  {
    oadd[1] = CrAnMask + i + Mup.Dt.val.N + Mup.Dt.con.N + Mdn.Dt.val.N + 2;
    oadd[0] = i + Mup.Dt.val.N + Mup.Dt.con.N + Mdn.Dt.val.N + 3;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.Dt.con.b[i - 1], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i + Mup.Dt.val.N + Mup.Dt.con.N + Mdn.Dt.val.N + 3;
    oadd[0] = i + Mup.Dt.val.N + Mup.Dt.con.N + Mdn.Dt.val.N + 2;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.Dt.con.b[i - 1], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i + Mup.Dt.val.N + Mup.Dt.con.N + Mdn.Dt.val.N + 3;
    oadd[0] = i + Mup.Dt.val.N + Mup.Dt.con.N + Mdn.Dt.val.N + 3;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.Dt.con.a[i], &QDetNotOrderedFcFaBcBa);
  }
  return 0;
}

/*  Creates an anderson impurity Hamiltonian
 *
 *  The impurity is at index i=0 for spin up and i=1 for spin down
 *  The spin up bath is from i=2                     to i=Mup.val.N+Mup.con.N+1
 *  The spin dn bath is from i=Mup.val.N+Mup.con.N+2 to i=Mup.val.N+Mup.con.N+Mdn.val.N+Mdn.con.N+1
 */
uint32_t DoubleTriDiagonalMatrixToAndersonImpurityOperator(OperatorType * O, DoubleTriDiagonalMatrixType Mup,
                                                           DoubleTriDiagonalMatrixType Mdn, double U, double mu)
{
  uint32_t i;
  uint16_t oadd[4];
  InitNamedOperator(O, "Single Site Anderson Impurity Hamiltonian");
  O->NFermion = Mup.val.N + Mup.con.N + Mdn.val.N + Mdn.con.N + 2;
  O->NWord = Min(O->NFermion, 1 + (O->NFermion - 1) / 8);
  InitOperatorFixedLength(O, 2, 8);
  InitOperatorFixedLength(O, 4, 8);
  oadd[3] = CrAnMask + 0;
  oadd[2] = CrAnMask + 1;
  oadd[1] = 1;
  oadd[0] = 0;
  RealOperatorAddLaderNormalOrder(O, oadd, 4, U, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 0;
  oadd[0] = 0;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, mu, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 1;
  oadd[0] = 1;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, mu, &QDetNotOrderedFcFaBcBa);
  /* spin up valence */
  oadd[1] = CrAnMask + 0;
  oadd[0] = 0;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.val.a[0] / 2.0, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 0;
  oadd[0] = 2;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.val.b[0], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 2;
  oadd[0] = 0;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.val.b[0], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 2;
  oadd[0] = 2;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.val.a[1], &QDetNotOrderedFcFaBcBa);
  for (i = 2; i <= Mup.val.N; i++)
  {
    oadd[1] = CrAnMask + i;
    oadd[0] = i + 1;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.val.b[i - 1], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i + 1;
    oadd[0] = i;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.val.b[i - 1], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i + 1;
    oadd[0] = i + 1;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.val.a[i], &QDetNotOrderedFcFaBcBa);
  }
  /* spin up conductance */
  oadd[1] = CrAnMask + 0;
  oadd[0] = 0;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.con.a[0] / 2.0, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 0;
  oadd[0] = Mup.val.N + 2;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.con.b[0], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + Mup.val.N + 2;
  oadd[0] = 0;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.con.b[0], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + Mup.val.N + 2;
  oadd[0] = Mup.val.N + 2;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.con.a[1], &QDetNotOrderedFcFaBcBa);
  for (i = 2; i <= Mup.con.N; i++)
  {
    oadd[1] = CrAnMask + i + Mup.val.N;
    oadd[0] = i + 1 + Mup.val.N;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.con.b[i - 1], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i + 1 + Mup.val.N;
    oadd[0] = i + Mup.val.N;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.con.b[i - 1], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i + 1 + Mup.val.N;
    oadd[0] = i + 1 + Mup.val.N;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.con.a[i], &QDetNotOrderedFcFaBcBa);
  }
  /* spin dn valence */
  oadd[1] = CrAnMask + 1;
  oadd[0] = 1;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.val.a[0] / 2.0, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 1;
  oadd[0] = Mup.val.N + Mup.con.N + 2;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.val.b[0], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + Mup.val.N + Mup.con.N + 2;
  oadd[0] = 1;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.val.b[0], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + Mup.val.N + Mup.con.N + 2;
  oadd[0] = Mup.val.N + Mup.con.N + 2;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.val.a[1], &QDetNotOrderedFcFaBcBa);
  for (i = 2; i <= Mdn.val.N; i++)
  {
    oadd[1] = CrAnMask + i + Mup.val.N + Mup.con.N;
    oadd[0] = i + Mup.val.N + Mup.con.N + 1;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.val.b[i - 1], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i + Mup.val.N + Mup.con.N + 1;
    oadd[0] = i + Mup.val.N + Mup.con.N;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.val.b[i - 1], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i + Mup.val.N + Mup.con.N + 1;
    oadd[0] = i + Mup.val.N + Mup.con.N + 1;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.val.a[i], &QDetNotOrderedFcFaBcBa);
  }
  /* spin dn conductance */
  oadd[1] = CrAnMask + 1;
  oadd[0] = 1;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.con.a[0] / 2.0, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 1;
  oadd[0] = Mup.val.N + Mup.con.N + Mdn.val.N + 2;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.con.b[0], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + Mup.val.N + Mup.con.N + Mdn.val.N + 2;
  oadd[0] = 1;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.con.b[0], &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + Mup.val.N + Mup.con.N + Mdn.val.N + 2;
  oadd[0] = Mup.val.N + Mup.con.N + Mdn.val.N + 2;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.con.a[1], &QDetNotOrderedFcFaBcBa);
  for (i = 2; i <= Mdn.con.N; i++)
  {
    oadd[1] = CrAnMask + i + Mup.val.N + Mup.con.N + Mdn.val.N;
    oadd[0] = i + Mup.val.N + Mup.con.N + Mdn.val.N + 1;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.con.b[i - 1], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i + Mup.val.N + Mup.con.N + Mdn.val.N + 1;
    oadd[0] = i + Mup.val.N + Mup.con.N + Mdn.val.N;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.con.b[i - 1], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i + Mup.val.N + Mup.con.N + Mdn.val.N + 1;
    oadd[0] = i + Mup.val.N + Mup.con.N + Mdn.val.N + 1;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.con.a[i], &QDetNotOrderedFcFaBcBa);
  }
  return 0;
}



/*  Creates an anderson impurity Hamiltonian
 *
 *  The impurity is at index i=0 for spin up and i=1 for spin down
 *  The first bath site is at index i=2 for spin up and i=3 for spin down
 *  The spin up bath is from i=4                     to i=Mup.val.N+Mup.con.N+3
 *  The spin dn bath is from i=Mup.val.N+Mup.con.N+4 to i=Mup.val.N+Mup.con.N+Mdn.val.N+Mdn.con.N+3
 */
uint32_t BADoubleTriDiagonalMatrixToAndersonImpurityOperator(OperatorType * O, BADoubleTriDiagonalMatrixType Mup,
                                                             BADoubleTriDiagonalMatrixType Mdn, double U, double mu)
{
  uint32_t i;
  uint16_t oadd[4];
  InitNamedOperator(O, "Single Site Anderson Impurity Hamiltonian");
  O->NFermion = Mup.val.N + Mup.con.N + Mdn.val.N + Mdn.con.N;
  O->NWord = Min(O->NFermion, 1 + (O->NFermion - 1) / 8);
  InitOperatorFixedLength(O, 2, 8);
  InitOperatorFixedLength(O, 4, 8);
  /* impurity */
  oadd[3] = CrAnMask + 0;
  oadd[2] = CrAnMask + 1;
  oadd[1] = 1;
  oadd[0] = 0;
  RealOperatorAddLaderNormalOrder(O, oadd, 4, U, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 0;
  oadd[0] = 0;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.acvB - mu, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 1;
  oadd[0] = 1;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.acvB - mu, &QDetNotOrderedFcFaBcBa);
  /* first bath site spin up */
  oadd[1] = CrAnMask + 2;
  oadd[0] = 2;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.acvA, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 0;
  oadd[0] = 2;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.acvBA, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 2;
  oadd[0] = 0;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.acvBA, &QDetNotOrderedFcFaBcBa);
  /* first bath site spin dn */
  oadd[1] = CrAnMask + 3;
  oadd[0] = 3;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.acvA, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 1;
  oadd[0] = 3;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.acvBA, &QDetNotOrderedFcFaBcBa);
  oadd[1] = CrAnMask + 3;
  oadd[0] = 1;
  RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.acvBA, &QDetNotOrderedFcFaBcBa);
  /* spin up valence */
  if(Mup.val.N>1)
  {
    oadd[1] = CrAnMask + 0;
    oadd[0] = 4;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.bvB0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4;
    oadd[0] = 0;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.bvB0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 2;
    oadd[0] = 4;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.bvA0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4;
    oadd[0] = 2;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.bvA0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4;
    oadd[0] = 4;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.val.a[2], &QDetNotOrderedFcFaBcBa);
  }
  for (i = 2; i < Mup.val.N; i++)
  {
    oadd[1] = CrAnMask + i + 2;
    oadd[0] = i + 3;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.val.b[i], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i + 3;
    oadd[0] = i + 2;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.val.b[i], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i + 3;
    oadd[0] = i + 3;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.val.a[i + 1], &QDetNotOrderedFcFaBcBa);
  }
  /* spin up conductance */
  if(Mup.con.N>1)
  {
    oadd[1] = CrAnMask + 0;
    oadd[0] = 4 + Mup.val.N - 1;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.bcB0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4 + Mup.val.N - 1;
    oadd[0] = 0;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.bcB0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 2;
    oadd[0] = 4 + Mup.val.N - 1;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.bcA0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4 + Mup.val.N - 1;
    oadd[0] = 2;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.bcA0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4 + Mup.val.N - 1;
    oadd[0] = 4 + Mup.val.N - 1;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.con.a[2], &QDetNotOrderedFcFaBcBa);
  }
  for (i = 2; i < Mup.con.N; i++)
  {
    oadd[1] = CrAnMask + i + 2 + Mup.val.N - 1;
    oadd[0] = i + 3 + Mup.val.N - 1;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.con.b[i], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i + 3 + Mup.val.N - 1;
    oadd[0] = i + 2 + Mup.val.N - 1;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.con.b[i], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i + 3 + Mup.val.N - 1;
    oadd[0] = i + 3 + Mup.val.N - 1;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mup.con.a[i + 1], &QDetNotOrderedFcFaBcBa);
  }
  /* spin dn valence */
  if(Mdn.val.N>1)
  {
    oadd[1] = CrAnMask + 1;
    oadd[0] = 4 + Mup.val.N + Mup.con.N - 2;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.bvB0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4 + Mup.val.N + Mup.con.N - 2;
    oadd[0] = 1;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.bvB0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 3;
    oadd[0] = 4 + Mup.val.N + Mup.con.N - 2;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.bvA0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4 + Mup.val.N + Mup.con.N - 2;
    oadd[0] = 3;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.bvA0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4 + Mup.val.N + Mup.con.N - 2;
    oadd[0] = 4 + Mup.val.N + Mup.con.N - 2;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.val.a[2], &QDetNotOrderedFcFaBcBa);
  }
  for (i = 2; i < Mdn.val.N; i++)
  {
    oadd[1] = CrAnMask + i + 2 + Mup.val.N + Mup.con.N - 2;
    oadd[0] = i + 3 + Mup.val.N + Mup.con.N - 2;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.val.b[i], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i + 3 + Mup.val.N + Mup.con.N - 2;
    oadd[0] = i + 2 + Mup.val.N + Mup.con.N - 2;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.val.b[i], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i + 3 + Mup.val.N + Mup.con.N - 2;
    oadd[0] = i + 3 + Mup.val.N + Mup.con.N - 2;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.val.a[i + 1], &QDetNotOrderedFcFaBcBa);
  }
  /* spin dn conductance */
  if(Mdn.con.N>1)
  {
    oadd[1] = CrAnMask + 1;
    oadd[0] = 4 + Mup.val.N + Mup.con.N + Mdn.val.N - 3;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.bcB0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4 + Mup.val.N + Mup.con.N + Mdn.val.N - 3;
    oadd[0] = 1;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.bcB0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 3;
    oadd[0] = 4 + Mup.val.N + Mup.con.N + Mdn.val.N - 3;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.bcA0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4 + Mup.val.N + Mup.con.N + Mdn.val.N - 3;
    oadd[0] = 3;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.bcA0, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4 + Mup.val.N + Mup.con.N + Mdn.val.N - 3;
    oadd[0] = 4 + Mup.val.N + Mup.con.N + Mdn.val.N - 3;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.con.a[2], &QDetNotOrderedFcFaBcBa);
  }
  for (i = 2; i < Mdn.con.N; i++)
  {
    oadd[1] = CrAnMask + i + 2 + Mup.val.N + Mup.con.N + Mdn.val.N - 3;
    oadd[0] = i + 3 + Mup.val.N + Mup.con.N + Mdn.val.N - 3;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.con.b[i], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i + 3 + Mup.val.N + Mup.con.N + Mdn.val.N - 3;
    oadd[0] = i + 2 + Mup.val.N + Mup.con.N + Mdn.val.N - 3;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.con.b[i], &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + i + 3 + Mup.val.N + Mup.con.N + Mdn.val.N - 3;
    oadd[0] = i + 3 + Mup.val.N + Mup.con.N + Mdn.val.N - 3;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, Mdn.con.a[i + 1], &QDetNotOrderedFcFaBcBa);
  }
  return 0;
}
