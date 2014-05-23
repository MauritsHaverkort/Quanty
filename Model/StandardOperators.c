/*  In the following functions the flag OrbitalSymmetry defines the atomic-like basis set used.
 *  Please not that in real calculations one very offten encounters some form of localize Wannier functions 
 *  that are similar to the atomic functions, but naturally not the same. It is those functions with the proper
 *  integrals I have in mind when labeling these.
 *
 *  Here all numbers are binary
 *    0*** 0**0 - No spin degeneracy
 *    0*** 0**1 - With spin degeneracy
 *    0*** 0*10 - only spin up is included
 *    0*** 0*00 - only spin down is included
 *    0*** 0*01 - Aternating down spin up spin per orbital
 *    0*** 0*11 - First a blok of down spin for all orbitals, then a block of up spin for all orbitals
 *
 *    0*00 0*** - Basis of Spherical Harmonics indexed from lz=-l to lz=l
 *    0*10 0*** - Basis of Tesseral Harmonics indexed accourding to ...
 *    0*11 0*** - Basis of Cubic Harmonics indexed according to ...
 *
 *    1*** 1*** - Basises where spin and orbital states are mixed, i.e. the j jz basis states
 */




/* atomic orbitals are entered on several possible basis sets.
 * The flag OrbitalSymmetry can be set accordingly.
 *
 * Bit wise
 * 1)
 * 2) 00 spin degenerate: 01 only down 11: only up
 * 3)
 * 4) 00 complex spherical harmonics 01: tesseral harmonics 10: cubic harmonics 11: undifined
 *
 * Integer wise:
 * OrbitalSymmetry &  1 ==  0 => Spin degenerate basis all functions are alternated dn up dn up dn etc.
 * OrbitalSymmetry &  1 ==  1 => No spin degeneracy, next bit decides if up or dn is included
 * OrbitalSymmetry &  3 ==  3 => No spin degeneracy, only spin up Sz= 1/2 included
 * OribtalSymmetry &  3 ==  1 => No spin degeneracy, only spin dn Sz=-1/2 included
 * OrbitalSymmetry & 12 ==  0 => Basis of complex spherical harmonics ordered from ml=-l to ml=l
 * OrbitalSymmetry & 12 ==  4 => Basis of tesseral harmonics ordered ....
 * OrbitalSymmetry & 12 ==  8 => Basis of cubic harmonics ordered ...
 * OrbitalSymmetry & 12 == 12 => non defined
 */


uint16_t lmsigmaToIndex(int const l, int const m, int const sigma, uint16_t const *LindexUp, uint16_t const *LindexDn)
{
  if (sigma > 0)
    return LindexUp[l + m];
  return LindexDn[l + m];
}

/*  CreateOperatorCoulombSingelShell creates the coulomb operator in a single shell of angular momentum L
 *  Input:
 *    NFermion     - The number of Fermions in the basis
 *    L            - The total angular momentum (2L+1) <= NFermion
 *    *LindexUp    - An array of legnth 2L+1 indexing the ml values from -L to L with spin up
 *    *LindexDn    - An array of legnth 2L+1 indexing the ml values from -L to L with spin down
 *    *F           - An array of length L+1 with the Slater integrals, F[0], F[2], F[4], ...
 *    The operator U is initialized on input
 */
uint32_t CreateOperatorCoulombSingleShell(uint32_t const NFermion, uint32_t const L, uint16_t const *LindexUp,
                                          uint16_t const *LindexDn, double const *F, OperatorType * const U)
{
  uint16_t oadd[4];
  double val;
  int    k, l, m1, m2, m3, m4, s1, s2;
  l = (int)L;
  InitNamedOperator(U, "Coulomb Operator");
  U->NFermion = NFermion;
  U->NWord = Min(U->NFermion, 1 + (U->NFermion - 1) / 8);
  InitOperatorFixedLength(U, 2, 8);
  InitOperatorFixedLength(U, 4, 8);
  for (k = 0; k < 2 * (l + 1); k += 2)
  {
    for (m1 = -l; m1 <= l; m1++)
    {
      for (m2 = -l; m2 <= l; m2++)
      {
        for (m3 = Max(-l, m1 + m2 - l); m3 <= Min(l, m1 + m2 + l); m3++)
        {
          /* -l >= m1+m2-m3 >= l
           *  thus m3-l >= m1 + m2 >= m3+l
           */
          m4 = m1 + m2 - m3;
          val = -0.5 * SlaterCoefficientC(k, l, m1, l, m3) * SlaterCoefficientC(k, l, m4, l, m2);
          if (Abs(val) > DBL_EPSILON)
          {
            for (s1 = 0; s1 <= 1; s1++)
            {
              for (s2 = 0; s2 <= 1; s2++)
              {
                oadd[3] = CrAnMask + (s1 == 0 ? LindexDn[l + m1] : LindexUp[l + m1]);
                oadd[2] = CrAnMask + (s2 == 0 ? LindexDn[l + m2] : LindexUp[l + m2]);
                oadd[1] = (s1 == 0 ? LindexDn[l + m3] : LindexUp[l + m3]);
                oadd[0] = (s2 == 0 ? LindexDn[l + m4] : LindexUp[l + m4]);
                RealOperatorAddLaderNormalOrder(U, oadd, 4, val * F[k / 2], &QDetNotOrderedFcFaBcBa);
              }
            }
          }
        }
      }
    }
  }
  OperatorChop(U, 100 * DBL_EPSILON);
  return 0;
}

uint32_t CreateOperatorCoulombTwoShellsConstantOccupation(uint32_t const NFermion, uint32_t const L1, uint32_t const L2,
                                                          uint16_t const *Lindex1Up, uint16_t const *Lindex1Dn,
                                                          uint16_t const *Lindex2Up, uint16_t const *Lindex2Dn,
                                                          double const *F, double const *G, OperatorType * const U)
{
  uint16_t oadd[4];
  double val;
  int    k, l1, l2, m1, m2, m3, m4, s1, s2;
  l1 = (int)L1;
  l2 = (int)L2;
  InitNamedOperator(U, "Coulomb Operator");
  U->NFermion = NFermion;
  U->NWord = Min(U->NFermion, 1 + (U->NFermion - 1) / 8);
  InitOperatorFixedLength(U, 2, 8);
  InitOperatorFixedLength(U, 4, 8);
  /* case 1 m1 and m3 belong to l1, m2 and m4 belong to l2 */
  for (k = 0; k <= 2 * (Min(l1, l2) + 1); k += 2)
  {
    for (m1 = -l1; m1 <= l1; m1++)
    {
      for (m2 = -l2; m2 <= l2; m2++)
      {
        for (m3 = Max(-l1, m1 + m2 - l2); m3 <= Min(l1, m1 + m2 + l2); m3++)
        {
          /* -l2 <= m4 <= l2
           * -l2 <= m1+m2-m3 <= l2
           *  thus m3-l2 <= m1 + m2 <= m3+l2
           */
          m4 = m1 + m2 - m3;
          val = -SlaterCoefficientC(k, l1, m1, l1, m3) * SlaterCoefficientC(k, l2, m4, l2, m2);
          if (Abs(val) > DBL_EPSILON)
          {
            for (s1 = 0; s1 <= 1; s1++)
            {
              for (s2 = 0; s2 <= 1; s2++)
              {
                oadd[3] = CrAnMask + (s1 == 0 ? Lindex1Dn[l1 + m1] : Lindex1Up[l1 + m1]);
                oadd[2] = CrAnMask + (s2 == 0 ? Lindex2Dn[l2 + m2] : Lindex2Up[l2 + m2]);
                oadd[1] = (s1 == 0 ? Lindex1Dn[l1 + m3] : Lindex1Up[l1 + m3]);
                oadd[0] = (s2 == 0 ? Lindex2Dn[l2 + m4] : Lindex2Up[l2 + m4]);
                RealOperatorAddLaderNormalOrder(U, oadd, 4, val * F[k / 2], &QDetNotOrderedFcFaBcBa);
              }
            }
          }
        }
      }
    }
  }
  /* case 2 m1 and m4 belong to l1, m2 and m3 belong to l2 */
  for (k = Abs(l1 - l2); k <= l1 + l2; k += 2)
  {
    for (m1 = -l1; m1 <= l1; m1++)
    {
      for (m2 = -l2; m2 <= l2; m2++)
      {
        for (m3 = Max(-l2, m1 + m2 - l1); m3 <= Min(l2, m1 + m2 + l1); m3++)
        {
          /* -l1 >= m4 >= l1
           * -l1 >= m1+m2-m3 >= l1
           *  thus m3-l1 >= m1 + m2 >= m3+l1
           */
          m4 = m1 + m2 - m3;
          val = -SlaterCoefficientC(k, l1, m1, l2, m3) * SlaterCoefficientC(k, l1, m4, l2, m2);
          if (Abs(val) > DBL_EPSILON)
          {
            for (s1 = 0; s1 <= 1; s1++)
            {
              for (s2 = 0; s2 <= 1; s2++)
              {
                oadd[3] = CrAnMask + (s1 == 0 ? Lindex1Dn[l1 + m1] : Lindex1Up[l1 + m1]);
                oadd[2] = CrAnMask + (s2 == 0 ? Lindex2Dn[l2 + m2] : Lindex2Up[l2 + m2]);
                oadd[1] = (s1 == 0 ? Lindex2Dn[l2 + m3] : Lindex2Up[l2 + m3]);
                oadd[0] = (s2 == 0 ? Lindex1Dn[l1 + m4] : Lindex1Up[l1 + m4]);
                RealOperatorAddLaderNormalOrder(U, oadd, 4, val * G[(k - (Abs(l1 - l2))) / 2], &QDetNotOrderedFcFaBcBa);
              }
            }
          }
        }
      }
    }
  }
  OperatorChop(U, 100 * DBL_EPSILON);
  return 0;
}

/*  Operators Sx, Sy, Sz, Splus Smin and S^2
 *  Input:
 *    NFermion    - The number of Fermions in the basis
 *    NSpin       - The number of Fermions on which the spin acts. NSpin <= NFermion/2
 *    *SpinUp     - An array of length NSpin with the index of Fermions with spin up
 *    *SpinDn     - An array of legnth NSpin with the index of Fermions with spin dn
 *                - The array SpinUp and SpinDn are related such that the n-th index of each relate to the same particle with either spin up or dn
 *                - Spin quantization is in the z direction
 *    The operators are not initialized on input
 *
 *  Output:
 *    Operators acting on NFermions representing the spin operators
 */
uint32_t CreateOperatorSx(uint32_t NFermion, uint32_t NSpin, uint16_t * SpinUp, uint16_t * SpinDn,
                          OperatorType * const Sx)
{
  uint32_t i;
  uint16_t oadd[2];
  InitNamedOperator(Sx, "Sx");
  Sx->NFermion = NFermion;
  Sx->NWord = Min(Sx->NFermion, 1 + (Sx->NFermion - 1) / 8);
  InitOperatorFixedLength(Sx, 2, 8);
  for (i = 0; i < NSpin; i++)
  {
    oadd[1] = CrAnMask + SpinDn[i];
    oadd[0] = SpinUp[i];
    RealOperatorAddLaderNormalOrder(Sx, oadd, 2, 0.5, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + SpinUp[i];
    oadd[0] = SpinDn[i];
    RealOperatorAddLaderNormalOrder(Sx, oadd, 2, 0.5, &QDetNotOrderedFcFaBcBa);
  }
  OperatorChop(Sx, 0.001);
  return 0;
}

uint32_t CreateOperatorSy(uint32_t NFermion, uint32_t NSpin, uint16_t * SpinUp, uint16_t * SpinDn,
                          OperatorType * const Sy)
{
  uint32_t i;
  uint16_t oadd[2];
  InitNamedOperator(Sy, "Sy");
  Sy->NFermion = NFermion;
  Sy->NWord = Min(Sy->NFermion, 1 + (Sy->NFermion - 1) / 8);
  InitOperatorFixedLength(Sy, 2, 8);
  MakeRealOperatorFixedLengthComplex(Sy, 2);
  Sy->QComplex = 1;
  for (i = 0; i < NSpin; i++)
  {
    oadd[1] = CrAnMask + SpinDn[i];
    oadd[0] = SpinUp[i];
    ComplexOperatorAddLaderNormalOrder(Sy, oadd, 2, 0.0, -0.5, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + SpinUp[i];
    oadd[0] = SpinDn[i];
    ComplexOperatorAddLaderNormalOrder(Sy, oadd, 2, 0.0, 0.5, &QDetNotOrderedFcFaBcBa);
  }
  OperatorChop(Sy, 0.001);
  return 0;
}

uint32_t CreateOperatorSz(uint32_t NFermion, uint32_t NSpin, uint16_t * SpinUp, uint16_t * SpinDn,
                          OperatorType * const Sz)
{
  uint32_t i;
  uint16_t oadd[2];
  InitNamedOperator(Sz, "Sz");
  Sz->NFermion = NFermion;
  Sz->NWord = Min(Sz->NFermion, 1 + (Sz->NFermion - 1) / 8);
  InitOperatorFixedLength(Sz, 2, 8);
  for (i = 0; i < NSpin; i++)
  {
    oadd[1] = CrAnMask + SpinUp[i];
    oadd[0] = SpinUp[i];
    RealOperatorAddLaderNormalOrder(Sz, oadd, 2, 0.5, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + SpinDn[i];
    oadd[0] = SpinDn[i];
    RealOperatorAddLaderNormalOrder(Sz, oadd, 2, -0.5, &QDetNotOrderedFcFaBcBa);
  }
  OperatorChop(Sz, 0.001);
  return 0;
}

uint32_t CreateOperatorSsqr(uint32_t NFermion, uint32_t NSpin, uint16_t * SpinUp, uint16_t * SpinDn,
                            OperatorType * const Ssqr)
{
  uint32_t i, j;
  uint16_t oadd[4];
  InitNamedOperator(Ssqr, "S^2");
  Ssqr->NFermion = NFermion;
  Ssqr->NWord = Min(Ssqr->NFermion, 1 + (Ssqr->NFermion - 1) / 8);
  InitOperatorFixedLength(Ssqr, 2, 8);
  InitOperatorFixedLength(Ssqr, 4, 8);
  for (i = 0; i < NSpin; i++)
  {
    oadd[1] = CrAnMask + SpinUp[i];
    oadd[0] = SpinUp[i];
    RealOperatorAddLaderNormalOrder(Ssqr, oadd, 2, 0.75, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + SpinDn[i];
    oadd[0] = SpinDn[i];
    RealOperatorAddLaderNormalOrder(Ssqr, oadd, 2, 0.75, &QDetNotOrderedFcFaBcBa);
    for (j = 0; j < NSpin; j++)
    {
      oadd[3] = CrAnMask + SpinUp[i];
      oadd[2] = CrAnMask + SpinUp[j];
      oadd[1] = SpinUp[i];
      oadd[0] = SpinUp[j];
      RealOperatorAddLaderNormalOrder(Ssqr, oadd, 4, -0.25, &QDetNotOrderedFcFaBcBa);
      oadd[3] = CrAnMask + SpinDn[i];
      oadd[2] = CrAnMask + SpinDn[j];
      oadd[1] = SpinDn[i];
      oadd[0] = SpinDn[j];
      RealOperatorAddLaderNormalOrder(Ssqr, oadd, 4, -0.25, &QDetNotOrderedFcFaBcBa);
      oadd[3] = CrAnMask + SpinUp[i];
      oadd[2] = CrAnMask + SpinDn[j];
      oadd[1] = SpinUp[i];
      oadd[0] = SpinDn[j];
      RealOperatorAddLaderNormalOrder(Ssqr, oadd, 4, 0.25, &QDetNotOrderedFcFaBcBa);
      oadd[3] = CrAnMask + SpinDn[i];
      oadd[2] = CrAnMask + SpinUp[j];
      oadd[1] = SpinDn[i];
      oadd[0] = SpinUp[j];
      RealOperatorAddLaderNormalOrder(Ssqr, oadd, 4, 0.25, &QDetNotOrderedFcFaBcBa);
      oadd[3] = CrAnMask + SpinUp[i];
      oadd[2] = CrAnMask + SpinDn[j];
      oadd[1] = SpinDn[i];
      oadd[0] = SpinUp[j];
      RealOperatorAddLaderNormalOrder(Ssqr, oadd, 4, -1.00, &QDetNotOrderedFcFaBcBa);
    }
  }
  OperatorChop(Ssqr, 0.001);
  return 0;
}

uint32_t CreateOperatorSplus(uint32_t NFermion, uint32_t NSpin, uint16_t * SpinUp, uint16_t * SpinDn,
                             OperatorType * const Splus)
{
  uint32_t i;
  uint16_t oadd[2];
  InitNamedOperator(Splus, "S^+");
  Splus->NFermion = NFermion;
  Splus->NWord = Min(Splus->NFermion, 1 + (Splus->NFermion - 1) / 8);
  InitOperatorFixedLength(Splus, 2, 8);
  for (i = 0; i < NSpin; i++)
  {
    oadd[1] = CrAnMask + SpinUp[i];
    oadd[0] = SpinDn[i];
    RealOperatorAddLaderNormalOrder(Splus, oadd, 2, 1, &QDetNotOrderedFcFaBcBa);
  }
  OperatorChop(Splus, 0.001);
  return 0;
}

uint32_t CreateOperatorSmin(uint32_t NFermion, uint32_t NSpin, uint16_t * SpinUp, uint16_t * SpinDn,
                            OperatorType * const Smin)
{
  uint32_t i;
  uint16_t oadd[2];
  InitNamedOperator(Smin, "S^-");
  Smin->NFermion = NFermion;
  Smin->NWord = Min(Smin->NFermion, 1 + (Smin->NFermion - 1) / 8);
  InitOperatorFixedLength(Smin, 2, 8);
  for (i = 0; i < NSpin; i++)
  {
    oadd[1] = CrAnMask + SpinDn[i];
    oadd[0] = SpinUp[i];
    RealOperatorAddLaderNormalOrder(Smin, oadd, 2, 1, &QDetNotOrderedFcFaBcBa);
  }
  OperatorChop(Smin, 0.001);
  return 0;
}

uint32_t CreateSpinOperators(uint32_t NFermion, uint32_t NSpin, uint16_t * SpinUp, uint16_t * SpinDn,
                             OperatorType * const Sx, OperatorType * const Sy, OperatorType * const Sz,
                             OperatorType * const S2)
{
  uint32_t i;
  uint16_t oadd[2];
  InitNamedOperator(Sx, "Sx");
  InitNamedOperator(Sy, "Sy");
  InitNamedOperator(Sz, "Sz");
  InitNamedOperator(S2, "S^2");
  Sx->NFermion = NFermion;
  Sy->NFermion = NFermion;
  Sz->NFermion = NFermion;
  S2->NFermion = NFermion;
  Sx->NWord = Min(Sx->NFermion, 1 + (Sx->NFermion - 1) / 8);
  Sy->NWord = Min(Sy->NFermion, 1 + (Sy->NFermion - 1) / 8);
  Sz->NWord = Min(Sz->NFermion, 1 + (Sz->NFermion - 1) / 8);
  S2->NWord = Min(S2->NFermion, 1 + (S2->NFermion - 1) / 8);
  InitOperatorFixedLength(Sx, 2, 8);
  InitOperatorFixedLength(Sy, 2, 8);
  MakeRealOperatorFixedLengthComplex(Sy, 2);
  Sy->QComplex = 1;
  InitOperatorFixedLength(Sz, 2, 8);
  InitOperatorFixedLength(S2, 2, 8);
  InitOperatorFixedLength(S2, 4, 8);
  for (i = 0; i < NSpin; i++)
  {
    oadd[1] = CrAnMask + SpinUp[i];
    oadd[0] = SpinUp[i];
    RealOperatorAddLaderNormalOrder(Sz, oadd, 2, 0.5, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + SpinDn[i];
    oadd[0] = SpinDn[i];
    RealOperatorAddLaderNormalOrder(Sz, oadd, 2, -0.5, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + SpinDn[i];
    oadd[0] = SpinUp[i];
    RealOperatorAddLaderNormalOrder(Sx, oadd, 2, 0.5, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + SpinUp[i];
    oadd[0] = SpinDn[i];
    RealOperatorAddLaderNormalOrder(Sx, oadd, 2, 0.5, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + SpinDn[i];
    oadd[0] = SpinUp[i];
    ComplexOperatorAddLaderNormalOrder(Sy, oadd, 2, 0.0, -0.5, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + SpinUp[i];
    oadd[0] = SpinDn[i];
    ComplexOperatorAddLaderNormalOrder(Sy, oadd, 2, 0.0, 0.5, &QDetNotOrderedFcFaBcBa);
  }
  OperatorMultiplyNormalOrder(S2, *Sy, *Sy, &QDetNotOrderedFcFaBcBa);
  OperatorMultiplyNormalOrder(S2, *Sx, *Sx, &QDetNotOrderedFcFaBcBa);
  OperatorMultiplyNormalOrder(S2, *Sz, *Sz, &QDetNotOrderedFcFaBcBa);
  OperatorChop(Sx, 0.001);
  OperatorChop(Sy, 0.001);
  OperatorChop(Sz, 0.001);
  OperatorChop(S2, 0.001);
  return 0;
}

/*  CreateOperatorL*** creates the operators Lx, Ly, Lz, Lplus, Lmin or L^2
 *  Input:
 *    NFermion    - The number of Fermions in the basis
 *    L           - The total angular momentum (2L+1) <= NFermion
 *    *Lindex     - An array of legnth 2L+1 indexing the ml values from -L to L
 *    The operators L*** are not initialized on input
 */
uint32_t CreateOperatorLxNoSpin(uint32_t NFermion, uint32_t L, uint16_t * Lindex, OperatorType * const Lx)
{
  uint32_t i;
  double m, l;
  uint16_t oadd[2];
  l = (double)L;
  InitNamedOperator(Lx, "Lx");
  Lx->NFermion = NFermion;
  Lx->NWord = Min(Lx->NFermion, 1 + (Lx->NFermion - 1) / 8);
  InitOperatorFixedLength(Lx, 2, 8);
  for (i = 0; i < 2 * L; i++)
  {
    m = (double)i - l;
    oadd[1] = CrAnMask + Lindex[i + 1];
    oadd[0] = Lindex[i];
    RealOperatorAddLaderNormalOrder(Lx, oadd, 2, 0.5 * sqrt(l + m + 1) * sqrt(l - m), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + Lindex[i];
    oadd[0] = Lindex[i + 1];
    RealOperatorAddLaderNormalOrder(Lx, oadd, 2, 0.5 * sqrt(l + m + 1) * sqrt(l - m), &QDetNotOrderedFcFaBcBa);
  }
  OperatorChop(Lx, 0.001);
  return 0;
}

uint32_t CreateOperatorLyNoSpin(uint32_t NFermion, uint32_t L, uint16_t * Lindex, OperatorType * const Ly)
{
  uint32_t i;
  double m, l;
  l = (double)L;
  uint16_t oadd[2];
  InitNamedOperator(Ly, "Ly");
  Ly->NWord = Min(Ly->NFermion, 1 + (Ly->NFermion - 1) / 8);
  InitOperatorFixedLength(Ly, 2, 8);
  MakeRealOperatorFixedLengthComplex(Ly, 2);
  Ly->QComplex = 1;
  for (i = 0; i < 2 * L; i++)
  {
    m = (double)i - l;
    oadd[1] = CrAnMask + Lindex[i + 1];
    oadd[0] = Lindex[i];
    ComplexOperatorAddLaderNormalOrder(Ly, oadd, 2, 0, -0.5 * sqrt(l + m + 1) * sqrt(l - m), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + Lindex[i];
    oadd[0] = Lindex[i + 1];
    ComplexOperatorAddLaderNormalOrder(Ly, oadd, 2, 0, 0.5 * sqrt(l + m + 1) * sqrt(l - m), &QDetNotOrderedFcFaBcBa);
  }
  OperatorChop(Ly, 0.001);
  return 0;
}

uint32_t CreateOperatorLzNoSpin(uint32_t NFermion, uint32_t L, uint16_t * Lindex, OperatorType * const Lz)
{
  uint32_t i;
  uint16_t oadd[2];
  InitNamedOperator(Lz, "Lz");
  Lz->NFermion = NFermion;
  Lz->NWord = Min(Lz->NFermion, 1 + (Lz->NFermion - 1) / 8);
  InitOperatorFixedLength(Lz, 2, 8);
  for (i = 0; i < 2 * L + 1; i++)
  {
    oadd[1] = CrAnMask + Lindex[i];
    oadd[0] = Lindex[i];
    RealOperatorAddLaderNormalOrder(Lz, oadd, 2, (double)i - (double)L, &QDetNotOrderedFcFaBcBa);
  }
  OperatorChop(Lz, 0.001);
  return 0;
}

uint32_t CreateOperatorLsqrNoSpin(uint32_t NFermion, uint32_t L, uint16_t * Lindex, OperatorType * const Lsqr)
{
  uint32_t i, j;
  double mi, mj, l;
  uint16_t oadd[4];
  l = (double)L;
  InitNamedOperator(Lsqr, "L^2");
  Lsqr->NFermion = NFermion;
  Lsqr->NWord = Min(Lsqr->NFermion, 1 + (Lsqr->NFermion - 1) / 8);
  InitOperatorFixedLength(Lsqr, 2, 8);
  InitOperatorFixedLength(Lsqr, 4, 8);
  for (i = 0; i < 2 * L + 1; i++)
  {
    mi = (double)i - l;
    oadd[1] = CrAnMask + Lindex[i];
    oadd[0] = Lindex[i];
    RealOperatorAddLaderNormalOrder(Lsqr, oadd, 2, l * (l + 1), &QDetNotOrderedFcFaBcBa);
    for (j = 0; j < 2 * L + 1; j++)
    {
      mj = (double)j - l;
      oadd[3] = CrAnMask + Lindex[i];
      oadd[2] = CrAnMask + Lindex[j];
      oadd[1] = Lindex[i];
      oadd[0] = Lindex[j];
      RealOperatorAddLaderNormalOrder(Lsqr, oadd, 4, -mi * mj, &QDetNotOrderedFcFaBcBa);
    }
  }
  for (i = 0; i < 2 * L; i++)
  {
    mi = (double)i - l;
    for (j = 0; j < 2 * L; j++)
    {
      mj = (double)j - l;
      oadd[3] = CrAnMask + Lindex[i + 1];
      oadd[2] = CrAnMask + Lindex[j];
      oadd[1] = Lindex[i];
      oadd[0] = Lindex[j + 1];
      RealOperatorAddLaderNormalOrder(Lsqr, oadd, 4, sqrt(l + mi + 1) * sqrt(l - mi) * sqrt(l - mj) * sqrt(l + mj + 1),
                                      &QDetNotOrderedFcFaBcBa);
    }
  }
  OperatorChop(Lsqr, 0.001);
  return 0;
}

uint32_t CreateOperatorLplusNoSpin(uint32_t NFermion, uint32_t L, uint16_t * Lindex, OperatorType * const Lplus)
{
  uint32_t i;
  double m, l;
  l = (double)L;
  uint16_t oadd[2];
  InitNamedOperator(Lplus, "L^+");
  Lplus->NFermion = NFermion;
  Lplus->NWord = Min(Lplus->NFermion, 1 + (Lplus->NFermion - 1) / 8);
  InitOperatorFixedLength(Lplus, 2, 8);
  for (i = 0; i < 2 * L; i++)
  {
    m = (double)i - l;
    oadd[1] = CrAnMask + Lindex[i + 1];
    oadd[0] = Lindex[i];
    RealOperatorAddLaderNormalOrder(Lplus, oadd, 2, sqrt(l + m + 1) * sqrt(l - m), &QDetNotOrderedFcFaBcBa);
  }
  OperatorChop(Lplus, 0.001);
  return 0;
}

uint32_t CreateOperatorLminNoSpin(uint32_t NFermion, uint32_t L, uint16_t * Lindex, OperatorType * const Lmin)
{
  uint32_t i;
  double m, l;
  l = (double)L;
  uint16_t oadd[2];
  InitNamedOperator(Lmin, "L^-");
  Lmin->NFermion = NFermion;
  Lmin->NWord = Min(Lmin->NFermion, 1 + (Lmin->NFermion - 1) / 8);
  InitOperatorFixedLength(Lmin, 2, 8);
  for (i = 1; i < 2 * L + 1; i++)
  {
    m = (double)i - l;
    oadd[1] = CrAnMask + Lindex[i - 1];
    oadd[0] = Lindex[i];
    RealOperatorAddLaderNormalOrder(Lmin, oadd, 2, sqrt(l - m + 1) * sqrt(l + m), &QDetNotOrderedFcFaBcBa);
  }
  OperatorChop(Lmin, 0.001);
  return 0;
}

/*  CreateOperatorL*** creates the operators Lx, Ly, Lz, Lplus, Lmin or L^2
 *  Input:
 *    NFermion    - The number of Fermions in the basis
 *    L           - The total angular momentum (2L+1) <= NFermion
 *    *LindexUp    - An array of legnth 2L+1 indexing the ml values from -L to L with spin up
 *    *LindexDn    - An array of legnth 2L+1 indexing the ml values from -L to L with spin down
 *    The operators L*** are not initialized on input
 */
uint32_t CreateOperatorLx(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn,
                          OperatorType * const Lx)
{
  uint32_t i;
  uint16_t oadd[2];
  double m, l;
  l = (double)L;
  InitNamedOperator(Lx, "Lx");
  Lx->NFermion = NFermion;
  Lx->NWord = Min(Lx->NFermion, 1 + (Lx->NFermion - 1) / 8);
  InitOperatorFixedLength(Lx, 2, 8);
  for (i = 0; i < 2 * L; i++)
  {
    m = (double)i - l;
    oadd[1] = CrAnMask + LindexUp[i + 1];
    oadd[0] = LindexUp[i];
    RealOperatorAddLaderNormalOrder(Lx, oadd, 2, 0.5 * sqrt(l + m + 1) * sqrt(l - m), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + LindexUp[i];
    oadd[0] = LindexUp[i + 1];
    RealOperatorAddLaderNormalOrder(Lx, oadd, 2, 0.5 * sqrt(l + m + 1) * sqrt(l - m), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + LindexDn[i + 1];
    oadd[0] = LindexDn[i];
    RealOperatorAddLaderNormalOrder(Lx, oadd, 2, 0.5 * sqrt(l + m + 1) * sqrt(l - m), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + LindexDn[i];
    oadd[0] = LindexDn[i + 1];
    RealOperatorAddLaderNormalOrder(Lx, oadd, 2, 0.5 * sqrt(l + m + 1) * sqrt(l - m), &QDetNotOrderedFcFaBcBa);
  }
  OperatorChop(Lx, 0.001);
  return 0;
}

uint32_t CreateOperatorLy(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn,
                          OperatorType * const Ly)
{
  uint32_t i;
  uint16_t oadd[2];
  double m, l;
  l = (double)L;
  InitNamedOperator(Ly, "Ly");
  Ly->NFermion = NFermion;
  Ly->NWord = Min(Ly->NFermion, 1 + (Ly->NFermion - 1) / 8);
  InitOperatorFixedLength(Ly, 2, 8);
  MakeRealOperatorFixedLengthComplex(Ly, 2);
  Ly->QComplex = 1;
  for (i = 0; i < 2 * L; i++)
  {
    m = (double)i - l;
    oadd[1] = CrAnMask + LindexUp[i + 1];
    oadd[0] = LindexUp[i];
    ComplexOperatorAddLaderNormalOrder(Ly, oadd, 2, 0, 0.5 * sqrt(l + m + 1) * sqrt(l - m), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + LindexUp[i];
    oadd[0] = LindexUp[i + 1];
    ComplexOperatorAddLaderNormalOrder(Ly, oadd, 2, 0, -0.5 * sqrt(l + m + 1) * sqrt(l - m), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + LindexDn[i + 1];
    oadd[0] = LindexDn[i];
    ComplexOperatorAddLaderNormalOrder(Ly, oadd, 2, 0, 0.5 * sqrt(l + m + 1) * sqrt(l - m), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + LindexDn[i];
    oadd[0] = LindexDn[i + 1];
    ComplexOperatorAddLaderNormalOrder(Ly, oadd, 2, 0, -0.5 * sqrt(l + m + 1) * sqrt(l - m), &QDetNotOrderedFcFaBcBa);
  }
  OperatorChop(Ly, 0.001);
  return 0;
}

uint32_t CreateOperatorLz(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn,
                          OperatorType * const Lz)
{
  uint32_t i;
  uint16_t oadd[2];
  double l;
  l = (double)L;
  InitNamedOperator(Lz, "Lz");
  Lz->NFermion = NFermion;
  Lz->NWord = Min(Lz->NFermion, 1 + (Lz->NFermion - 1) / 8);
  InitOperatorFixedLength(Lz, 2, 8);
  for (i = 0; i < 2 * L + 1; i++)
  {
    oadd[1] = CrAnMask + LindexUp[i];
    oadd[0] = LindexUp[i];
    RealOperatorAddLaderNormalOrder(Lz, oadd, 2, (double)i - l, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + LindexDn[i];
    oadd[0] = LindexDn[i];
    RealOperatorAddLaderNormalOrder(Lz, oadd, 2, (double)i - l, &QDetNotOrderedFcFaBcBa);
  }
  OperatorChop(Lz, 0.001);
  return 0;
}

uint32_t CreateOperatorLsqr(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn,
                            OperatorType * const Lsqr)
{
  uint32_t i, j;
  double mi, mj, l;
  uint16_t oadd[4];
  l = (double)L;
  InitNamedOperator(Lsqr, "L^2");
  Lsqr->NFermion = NFermion;
  Lsqr->NWord = Min(Lsqr->NFermion, 1 + (Lsqr->NFermion - 1) / 8);
  InitOperatorFixedLength(Lsqr, 2, 8);
  InitOperatorFixedLength(Lsqr, 4, 8);
  for (i = 0; i < 2 * L + 1; i++)
  {
    mi = (double)i - l;
    oadd[1] = CrAnMask + LindexUp[i];
    oadd[0] = LindexUp[i];
    RealOperatorAddLaderNormalOrder(Lsqr, oadd, 2, l * (l + 1), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + LindexDn[i];
    oadd[0] = LindexDn[i];
    RealOperatorAddLaderNormalOrder(Lsqr, oadd, 2, l * (l + 1), &QDetNotOrderedFcFaBcBa);
    for (j = 0; j < 2 * L + 1; j++)
    {
      mj = (double)j - l;
      oadd[3] = CrAnMask + LindexUp[i];
      oadd[2] = CrAnMask + LindexUp[j];
      oadd[1] = LindexUp[i];
      oadd[0] = LindexUp[j];
      RealOperatorAddLaderNormalOrder(Lsqr, oadd, 4, -mi * mj, &QDetNotOrderedFcFaBcBa);
      oadd[3] = CrAnMask + LindexDn[i];
      oadd[2] = CrAnMask + LindexDn[j];
      oadd[1] = LindexDn[i];
      oadd[0] = LindexDn[j];
      RealOperatorAddLaderNormalOrder(Lsqr, oadd, 4, -mi * mj, &QDetNotOrderedFcFaBcBa);
      oadd[3] = CrAnMask + LindexUp[i];
      oadd[2] = CrAnMask + LindexDn[j];
      oadd[1] = LindexUp[i];
      oadd[0] = LindexDn[j];
      RealOperatorAddLaderNormalOrder(Lsqr, oadd, 4, -mi * mj, &QDetNotOrderedFcFaBcBa);
      oadd[3] = CrAnMask + LindexDn[i];
      oadd[2] = CrAnMask + LindexUp[j];
      oadd[1] = LindexDn[i];
      oadd[0] = LindexUp[j];
      RealOperatorAddLaderNormalOrder(Lsqr, oadd, 4, -mi * mj, &QDetNotOrderedFcFaBcBa);
    }
  }
  for (i = 0; i < 2 * L; i++)
  {
    mi = (double)i - l;
    for (j = 0; j < 2 * L; j++)
    {
      mj = (double)j - l;
      oadd[3] = CrAnMask + LindexUp[i + 1];
      oadd[2] = CrAnMask + LindexUp[j];
      oadd[1] = LindexUp[i];
      oadd[0] = LindexUp[j + 1];
      RealOperatorAddLaderNormalOrder(Lsqr, oadd, 4, -sqrt(l + mi + 1) * sqrt(l - mi) * sqrt(l - mj) * sqrt(l + mj + 1),
                                      &QDetNotOrderedFcFaBcBa);
      oadd[3] = CrAnMask + LindexDn[i + 1];
      oadd[2] = CrAnMask + LindexDn[j];
      oadd[1] = LindexDn[i];
      oadd[0] = LindexDn[j + 1];
      RealOperatorAddLaderNormalOrder(Lsqr, oadd, 4, -sqrt(l + mi + 1) * sqrt(l - mi) * sqrt(l - mj) * sqrt(l + mj + 1),
                                      &QDetNotOrderedFcFaBcBa);
      oadd[3] = CrAnMask + LindexUp[i + 1];
      oadd[2] = CrAnMask + LindexDn[j];
      oadd[1] = LindexUp[i];
      oadd[0] = LindexDn[j + 1];
      RealOperatorAddLaderNormalOrder(Lsqr, oadd, 4, -sqrt(l + mi + 1) * sqrt(l - mi) * sqrt(l - mj) * sqrt(l + mj + 1),
                                      &QDetNotOrderedFcFaBcBa);
      oadd[3] = CrAnMask + LindexDn[i + 1];
      oadd[2] = CrAnMask + LindexUp[j];
      oadd[1] = LindexDn[i];
      oadd[0] = LindexUp[j + 1];
      RealOperatorAddLaderNormalOrder(Lsqr, oadd, 4, -sqrt(l + mi + 1) * sqrt(l - mi) * sqrt(l - mj) * sqrt(l + mj + 1),
                                      &QDetNotOrderedFcFaBcBa);
    }
  }
  OperatorChop(Lsqr, 0.001);
  return 0;
}

uint32_t CreateOperatorLplus(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn,
                             OperatorType * const Lplus)
{
  uint32_t i;
  double m, l;
  l = (double)L;
  uint16_t oadd[2];
  InitNamedOperator(Lplus, "L^+");
  Lplus->NFermion = NFermion;
  Lplus->NWord = Min(Lplus->NFermion, 1 + (Lplus->NFermion - 1) / 8);
  InitOperatorFixedLength(Lplus, 2, 8);
  for (i = 0; i < 2 * L; i++)
  {
    m = (double)i - l;
    oadd[1] = CrAnMask + LindexUp[i + 1];
    oadd[0] = LindexUp[i];
    RealOperatorAddLaderNormalOrder(Lplus, oadd, 2, sqrt(l + m + 1) * sqrt(l - m), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + LindexDn[i + 1];
    oadd[0] = LindexDn[i];
    RealOperatorAddLaderNormalOrder(Lplus, oadd, 2, sqrt(l + m + 1) * sqrt(l - m), &QDetNotOrderedFcFaBcBa);
  }
  OperatorChop(Lplus, 0.001);
  return 0;
}

uint32_t CreateOperatorLmin(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn,
                            OperatorType * const Lmin)
{
  uint32_t i;
  double m, l;
  l = (double)L;
  uint16_t oadd[2];
  InitNamedOperator(Lmin, "L^-");
  Lmin->NFermion = NFermion;
  Lmin->NWord = Min(Lmin->NFermion, 1 + (Lmin->NFermion - 1) / 8);
  InitOperatorFixedLength(Lmin, 2, 8);
  for (i = 1; i < 2 * L + 1; i++)
  {
    m = (double)i - L;
    oadd[1] = CrAnMask + LindexUp[i - 1];
    oadd[0] = LindexUp[i];
    RealOperatorAddLaderNormalOrder(Lmin, oadd, 2, sqrt(l - m + 1) * sqrt(l + m), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + LindexDn[i - 1];
    oadd[0] = LindexDn[i];
    RealOperatorAddLaderNormalOrder(Lmin, oadd, 2, sqrt(l - m + 1) * sqrt(l + m), &QDetNotOrderedFcFaBcBa);
  }
  OperatorChop(Lmin, 0.001);
  return 0;
}

/*  CreateOperatorJ*** creates the operators Jx, Jy, Jz, Jplus, Jmin or J^2
 *  Input:
 *    NFermion    - The number of Fermions in the basis
 *    L           - The total angular momentum (2L+1) <= NFermion
 *    *LindexUp    - An array of legnth 2L+1 indexing the ml values from -L to L with spin up
 *    *LindexDn    - An array of legnth 2L+1 indexing the ml values from -L to L with spin down
 *    The operators J*** are not initialized on input
 */
uint32_t CreateOperatorJx(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn,
                          OperatorType * const Jx)
{
  uint32_t i;
  double m, l;
  l = (double)L;
  uint16_t oadd[2];
  InitNamedOperator(Jx, "Jx");
  Jx->NFermion = NFermion;
  Jx->NWord = Min(Jx->NFermion, 1 + (Jx->NFermion - 1) / 8);
  InitOperatorFixedLength(Jx, 2, 8);
  for (i = 0; i < 2 * L; i++)
  {
    m = (double)i - l;
    oadd[1] = CrAnMask + LindexUp[i + 1];
    oadd[0] = LindexUp[i];
    RealOperatorAddLaderNormalOrder(Jx, oadd, 2, 0.5 * sqrt(L + m + 1) * sqrt(l - m), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + LindexUp[i];
    oadd[0] = LindexUp[i + 1];
    RealOperatorAddLaderNormalOrder(Jx, oadd, 2, 0.5 * sqrt(L + m + 1) * sqrt(l - m), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + LindexDn[i + 1];
    oadd[0] = LindexDn[i];
    RealOperatorAddLaderNormalOrder(Jx, oadd, 2, 0.5 * sqrt(L + m + 1) * sqrt(l - m), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + LindexDn[i];
    oadd[0] = LindexDn[i + 1];
    RealOperatorAddLaderNormalOrder(Jx, oadd, 2, 0.5 * sqrt(L + m + 1) * sqrt(l - m), &QDetNotOrderedFcFaBcBa);
  }
  for (i = 0; i < 2 * L + 1; i++)
  {
    oadd[1] = CrAnMask + LindexUp[i];
    oadd[0] = LindexDn[i];
    RealOperatorAddLaderNormalOrder(Jx, oadd, 2, 0.5, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + LindexDn[i];
    oadd[0] = LindexUp[i];
    RealOperatorAddLaderNormalOrder(Jx, oadd, 2, 0.5, &QDetNotOrderedFcFaBcBa);
  }
  OperatorChop(Jx, 0.001);
  return 0;
}

uint32_t CreateOperatorJy(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn,
                          OperatorType * const Jy)
{
  uint32_t i;
  double m, l;
  l = (double)L;
  uint16_t oadd[2];
  InitNamedOperator(Jy, "Jy");
  Jy->NFermion = NFermion;
  Jy->NWord = Min(Jy->NFermion, 1 + (Jy->NFermion - 1) / 8);
  InitOperatorFixedLength(Jy, 2, 8);
  MakeRealOperatorFixedLengthComplex(Jy, 2);
  Jy->QComplex = 1;
  for (i = 0; i < 2 * L; i++)
  {
    m = (double)i - l;
    oadd[1] = CrAnMask + LindexUp[i + 1];
    oadd[0] = LindexUp[i];
    ComplexOperatorAddLaderNormalOrder(Jy, oadd, 2, 0, 0.5 * sqrt(l + m + 1) * sqrt(l - m), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + LindexUp[i];
    oadd[0] = LindexUp[i + 1];
    ComplexOperatorAddLaderNormalOrder(Jy, oadd, 2, 0, -0.5 * sqrt(l + m + 1) * sqrt(l - m), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + LindexDn[i + 1];
    oadd[0] = LindexDn[i];
    ComplexOperatorAddLaderNormalOrder(Jy, oadd, 2, 0, 0.5 * sqrt(l + m + 1) * sqrt(l - m), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + LindexDn[i];
    oadd[0] = LindexDn[i + 1];
    ComplexOperatorAddLaderNormalOrder(Jy, oadd, 2, 0, -0.5 * sqrt(l + m + 1) * sqrt(l - m), &QDetNotOrderedFcFaBcBa);
  }
  for (i = 0; i < 2 * L + 1; i++)
  {
    oadd[1] = CrAnMask + LindexDn[i];
    oadd[0] = LindexUp[i];
    ComplexOperatorAddLaderNormalOrder(Jy, oadd, 2, 0.0, -0.5, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + LindexUp[i];
    oadd[0] = LindexDn[i];
    ComplexOperatorAddLaderNormalOrder(Jy, oadd, 2, 0.0, 0.5, &QDetNotOrderedFcFaBcBa);
  }
  OperatorChop(Jy, 0.001);
  return 0;
}

uint32_t CreateOperatorJz(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn,
                          OperatorType * const Jz)
{
  uint32_t i;
  uint16_t oadd[2];
  double l;
  l = (double)L;
  InitNamedOperator(Jz, "Jz");
  Jz->NFermion = NFermion;
  Jz->NWord = Min(Jz->NFermion, 1 + (Jz->NFermion - 1) / 8);
  InitOperatorFixedLength(Jz, 2, 8);
  for (i = 0; i < 2 * L + 1; i++)
  {
    oadd[1] = CrAnMask + LindexUp[i];
    oadd[0] = LindexUp[i];
    RealOperatorAddLaderNormalOrder(Jz, oadd, 2, (double)i - l, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + LindexDn[i];
    oadd[0] = LindexDn[i];
    RealOperatorAddLaderNormalOrder(Jz, oadd, 2, (double)i - l, &QDetNotOrderedFcFaBcBa);
  }
  for (i = 0; i < 2 * L + 1; i++)
  {
    oadd[1] = CrAnMask + LindexUp[i];
    oadd[0] = LindexUp[i];
    RealOperatorAddLaderNormalOrder(Jz, oadd, 2, 0.5, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + LindexDn[i];
    oadd[0] = LindexDn[i];
    RealOperatorAddLaderNormalOrder(Jz, oadd, 2, -0.5, &QDetNotOrderedFcFaBcBa);
  }
  OperatorChop(Jz, 0.001);
  return 0;
}

uint32_t CreateOperatorJsqr(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn,
                            OperatorType * const Jsqr)
{
  uint32_t i, j;
  double mi, mj, l;
  l = (double)L;
  uint16_t oadd[4];
  InitNamedOperator(Jsqr, "J^2");
  Jsqr->NFermion = NFermion;
  Jsqr->NWord = Min(Jsqr->NFermion, 1 + (Jsqr->NFermion - 1) / 8);
  InitOperatorFixedLength(Jsqr, 2, 8);
  InitOperatorFixedLength(Jsqr, 4, 8);
  for (i = 0; i < 2 * L + 1; i++)
  {
    mi = (double)i - l;
    oadd[1] = CrAnMask + LindexUp[i];
    oadd[0] = LindexUp[i];
    RealOperatorAddLaderNormalOrder(Jsqr, oadd, 2, l * (l + 1) + 0.75 + mi, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + LindexDn[i];
    oadd[0] = LindexDn[i];
    RealOperatorAddLaderNormalOrder(Jsqr, oadd, 2, l * (l + 1) + 0.75 - mi, &QDetNotOrderedFcFaBcBa);
    for (j = 0; j < 2 * L + 1; j++)
    {
      mj = (double)j - l;
      oadd[3] = CrAnMask + LindexUp[i];
      oadd[2] = CrAnMask + LindexUp[j];
      oadd[1] = LindexUp[i];
      oadd[0] = LindexUp[j];
      RealOperatorAddLaderNormalOrder(Jsqr, oadd, 4, -mi * mj - 0.25 - mi, &QDetNotOrderedFcFaBcBa);
      oadd[3] = CrAnMask + LindexDn[i];
      oadd[2] = CrAnMask + LindexDn[j];
      oadd[1] = LindexDn[i];
      oadd[0] = LindexDn[j];
      RealOperatorAddLaderNormalOrder(Jsqr, oadd, 4, -mi * mj - 0.25 + mi, &QDetNotOrderedFcFaBcBa);
      oadd[3] = CrAnMask + LindexUp[i];
      oadd[2] = CrAnMask + LindexDn[j];
      oadd[1] = LindexUp[i];
      oadd[0] = LindexDn[j];
      RealOperatorAddLaderNormalOrder(Jsqr, oadd, 4, -mi * mj + 0.25 + mi, &QDetNotOrderedFcFaBcBa);
      oadd[3] = CrAnMask + LindexDn[i];
      oadd[2] = CrAnMask + LindexUp[j];
      oadd[1] = LindexDn[i];
      oadd[0] = LindexUp[j];
      RealOperatorAddLaderNormalOrder(Jsqr, oadd, 4, -mi * mj + 0.25 - mi, &QDetNotOrderedFcFaBcBa);

      oadd[3] = CrAnMask + LindexUp[i];
      oadd[2] = CrAnMask + LindexDn[j];
      oadd[1] = LindexDn[i];
      oadd[0] = LindexUp[j];
      RealOperatorAddLaderNormalOrder(Jsqr, oadd, 4, -1, &QDetNotOrderedFcFaBcBa);
    }
  }
  for (i = 0; i < 2 * L; i++)
  {
    mi = (double)i - l;
    oadd[1] = CrAnMask + LindexDn[i + 1];
    oadd[0] = LindexUp[i];
    RealOperatorAddLaderNormalOrder(Jsqr, oadd, 2, sqrt(l + mi + 1) * sqrt(l - mi), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + LindexUp[i];
    oadd[0] = LindexDn[i + 1];
    RealOperatorAddLaderNormalOrder(Jsqr, oadd, 2, sqrt(l + mi + 1) * sqrt(l - mi), &QDetNotOrderedFcFaBcBa);
    for (j = 0; j < 2 * L; j++)
    {
      mj = (double)j - l;
      oadd[3] = CrAnMask + LindexUp[i + 1];
      oadd[2] = CrAnMask + LindexUp[j];
      oadd[1] = LindexUp[i];
      oadd[0] = LindexUp[j + 1];
      RealOperatorAddLaderNormalOrder(Jsqr, oadd, 4, -sqrt(l + mi + 1) * sqrt(l - mi) * sqrt(l - mj) * sqrt(l + mj + 1),
                                      &QDetNotOrderedFcFaBcBa);
      oadd[3] = CrAnMask + LindexDn[i + 1];
      oadd[2] = CrAnMask + LindexDn[j];
      oadd[1] = LindexDn[i];
      oadd[0] = LindexDn[j + 1];
      RealOperatorAddLaderNormalOrder(Jsqr, oadd, 4, -sqrt(l + mi + 1) * sqrt(l - mi) * sqrt(l - mj) * sqrt(l + mj + 1),
                                      &QDetNotOrderedFcFaBcBa);
      oadd[3] = CrAnMask + LindexUp[i + 1];
      oadd[2] = CrAnMask + LindexDn[j];
      oadd[1] = LindexUp[i];
      oadd[0] = LindexDn[j + 1];
      RealOperatorAddLaderNormalOrder(Jsqr, oadd, 4, -sqrt(l + mi + 1) * sqrt(l - mi) * sqrt(l - mj) * sqrt(l + mj + 1),
                                      &QDetNotOrderedFcFaBcBa);
      oadd[3] = CrAnMask + LindexDn[i + 1];
      oadd[2] = CrAnMask + LindexUp[j];
      oadd[1] = LindexDn[i];
      oadd[0] = LindexUp[j + 1];
      RealOperatorAddLaderNormalOrder(Jsqr, oadd, 4, -sqrt(l + mi + 1) * sqrt(l - mi) * sqrt(l - mj) * sqrt(l + mj + 1),
                                      &QDetNotOrderedFcFaBcBa);
    }
    for (j = 0; j < 2 * L + 1; j++)
    {
      mj = (double)j - l;
      oadd[3] = CrAnMask + LindexUp[i + 1];
      oadd[2] = CrAnMask + LindexDn[j];
      oadd[1] = LindexUp[i];
      oadd[0] = LindexUp[j];
      RealOperatorAddLaderNormalOrder(Jsqr, oadd, 4, -sqrt(l + mi + 1) * sqrt(l - mi), &QDetNotOrderedFcFaBcBa);
      oadd[3] = CrAnMask + LindexDn[i + 1];
      oadd[2] = CrAnMask + LindexDn[j];
      oadd[1] = LindexDn[i];
      oadd[0] = LindexUp[j];
      RealOperatorAddLaderNormalOrder(Jsqr, oadd, 4, -sqrt(l + mi + 1) * sqrt(l - mi), &QDetNotOrderedFcFaBcBa);
      oadd[3] = CrAnMask + LindexUp[i];
      oadd[2] = CrAnMask + LindexUp[j];
      oadd[1] = LindexUp[i + 1];
      oadd[0] = LindexDn[j];
      RealOperatorAddLaderNormalOrder(Jsqr, oadd, 4, -sqrt(l + mi + 1) * sqrt(l - mi), &QDetNotOrderedFcFaBcBa);
      oadd[3] = CrAnMask + LindexDn[i];
      oadd[2] = CrAnMask + LindexUp[j];
      oadd[1] = LindexDn[i + 1];
      oadd[0] = LindexDn[j];
      RealOperatorAddLaderNormalOrder(Jsqr, oadd, 4, -sqrt(l + mi + 1) * sqrt(l - mi), &QDetNotOrderedFcFaBcBa);
    }
  }
  OperatorChop(Jsqr, 0.001);
  return 0;
}

uint32_t CreateOperatorJplus(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn,
                             OperatorType * const Jplus)
{
  uint32_t i;
  double m, l;
  l = (double)L;
  uint16_t oadd[2];
  InitNamedOperator(Jplus, "J^+");
  Jplus->NFermion = NFermion;
  Jplus->NWord = Min(Jplus->NFermion, 1 + (Jplus->NFermion - 1) / 8);
  InitOperatorFixedLength(Jplus, 2, 8);
  for (i = 0; i < 2 * L; i++)
  {
    m = (double)i - l;
    oadd[1] = CrAnMask + LindexUp[i + 1];
    oadd[0] = LindexUp[i];
    RealOperatorAddLaderNormalOrder(Jplus, oadd, 2, sqrt(l + m + 1) * sqrt(l - m), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + LindexDn[i + 1];
    oadd[0] = LindexDn[i];
    RealOperatorAddLaderNormalOrder(Jplus, oadd, 2, sqrt(l + m + 1) * sqrt(l - m), &QDetNotOrderedFcFaBcBa);
  }
  for (i = 0; i < 2 * L + 1; i++)
  {
    oadd[1] = CrAnMask + LindexUp[i];
    oadd[0] = LindexDn[i];
    RealOperatorAddLaderNormalOrder(Jplus, oadd, 2, 1, &QDetNotOrderedFcFaBcBa);
  }
  OperatorChop(Jplus, 0.001);
  return 0;
}

uint32_t CreateOperatorJmin(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn,
                            OperatorType * const Jmin)
{
  uint32_t i;
  double m, l;
  l = (double)L;
  uint16_t oadd[2];
  InitNamedOperator(Jmin, "J^-");
  Jmin->NFermion = NFermion;
  Jmin->NWord = Min(Jmin->NFermion, 1 + (Jmin->NFermion - 1) / 8);
  InitOperatorFixedLength(Jmin, 2, 8);
  for (i = 1; i < 2 * L + 1; i++)
  {
    m = (double)i - l;
    oadd[1] = CrAnMask + LindexUp[i - 1];
    oadd[0] = LindexUp[i];
    RealOperatorAddLaderNormalOrder(Jmin, oadd, 2, sqrt(l - m + 1) * sqrt(l + m), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + LindexDn[i - 1];
    oadd[0] = LindexDn[i];
    RealOperatorAddLaderNormalOrder(Jmin, oadd, 2, sqrt(l - m + 1) * sqrt(l + m), &QDetNotOrderedFcFaBcBa);
  }
  for (i = 0; i < 2 * L + 1; i++)
  {
    oadd[1] = CrAnMask + LindexDn[i];
    oadd[0] = LindexUp[i];
    RealOperatorAddLaderNormalOrder(Jmin, oadd, 2, 1, &QDetNotOrderedFcFaBcBa);
  }
  OperatorChop(Jmin, 0.001);
  return 0;
}

uint32_t CreateOperatorldots(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn,
                             OperatorType * const ldots)
{
  uint32_t i;
  double mi, l;
  l = (double)L;
  uint16_t oadd[4];
  InitNamedOperator(ldots, "l.s");
  ldots->NFermion = NFermion;
  ldots->NWord = Min(ldots->NFermion, 1 + (ldots->NFermion - 1) / 8);
  InitOperatorFixedLength(ldots, 2, 8);
  InitOperatorFixedLength(ldots, 4, 8);
  for (i = 0; i < 2 * L + 1; i++)
  {
    mi = (double)i - l;
    oadd[1] = CrAnMask + LindexUp[i];
    oadd[0] = LindexUp[i];
    RealOperatorAddLaderNormalOrder(ldots, oadd, 2, mi * 0.5, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + LindexDn[i];
    oadd[0] = LindexDn[i];
    RealOperatorAddLaderNormalOrder(ldots, oadd, 2, -mi * 0.5, &QDetNotOrderedFcFaBcBa);
  }
  for (i = 0; i < 2 * L; i++)
  {
    mi = (double)i - l;
    oadd[1] = CrAnMask + LindexDn[i + 1];
    oadd[0] = LindexUp[i];
    RealOperatorAddLaderNormalOrder(ldots, oadd, 2, 0.5 * sqrt(l + mi + 1) * sqrt(l - mi), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + LindexUp[i];
    oadd[0] = LindexDn[i + 1];
    RealOperatorAddLaderNormalOrder(ldots, oadd, 2, 0.5 * sqrt(l + mi + 1) * sqrt(l - mi), &QDetNotOrderedFcFaBcBa);
  }
  OperatorChop(ldots, 0.001);
  return 0;
}

uint32_t CreateOperatorCF(uint32_t const NFermion, uint32_t LCr, uint32_t LAn, uint16_t * LCrindexUp,
                          uint16_t * LCrindexDn, uint16_t * LAnindexUp, uint16_t * LAnindexDn,
                          HarmonicsExpansion const Akm, OperatorType * const CF)
{
  uint32_t i;
  uint16_t oadd[2];
  int    mCr, mAn, lCr, lAn;
  double val;
  InitNamedOperator(CF, "CF");
  CF->NFermion = NFermion;
  CF->NWord = Min(CF->NFermion, 1 + (CF->NFermion - 1) / 8);
  InitOperatorFixedLength(CF, 2, 8);
  if (Akm.QComplex)
  {
    MakeRealOperatorFixedLengthComplex(CF, 2);
    CF->QComplex = 1;
  }
  lCr = (int)LCr;
  lAn = (int)LAn;
  for (i = 0; i < Akm.N; i++)
  {
    for (mCr = -lCr; mCr <= lCr; mCr++)
    {
      mAn = mCr - Akm.m[i];
      val = (double)SlaterCoefficientC(Akm.l[i], lCr, mCr, lAn, mAn);
      if (Abs(val) > DBL_EPSILON)
      {
        oadd[1] = CrAnMask + LCrindexDn[mCr + lCr];
        oadd[0] = LAnindexDn[mAn + lAn];
        if (Akm.QComplex)
        {
          ComplexOperatorAddLaderNormalOrder(CF, oadd, 2, val * Akm.Re[i], val * Akm.Im[i], &QDetNotOrderedFcFaBcBa);
        }
        else
        {
          RealOperatorAddLaderNormalOrder(CF, oadd, 2, val * Akm.Re[i], &QDetNotOrderedFcFaBcBa);
        }
        oadd[1] = CrAnMask + LCrindexUp[mCr + lCr];
        oadd[0] = LAnindexUp[mAn + lAn];
        if (Akm.QComplex)
        {
          ComplexOperatorAddLaderNormalOrder(CF, oadd, 2, val * Akm.Re[i], val * Akm.Im[i], &QDetNotOrderedFcFaBcBa);
        }
        else
        {
          RealOperatorAddLaderNormalOrder(CF, oadd, 2, val * Akm.Re[i], &QDetNotOrderedFcFaBcBa);
        }
      }
    }
  }
  return 0;
}



uint32_t CreateOperatorCreate(uint32_t const NFermion, uint16_t const Crindex, OperatorType * const OppCr)
{
  uint16_t oadd[1];
  InitNamedOperator(OppCr, "Cr");
  sprintf(OppCr->Name, "Creation operator on site %i", Crindex);
  OppCr->NFermion = NFermion;
  OppCr->NWord = Min(OppCr->NFermion, 1 + (OppCr->NFermion - 1) / 8);
  InitOperatorFixedLength(OppCr, 1, 8);
  oadd[0] = CrAnMask + Crindex;
  RealOperatorAddLader(OppCr, oadd, 1, 1);
  return 0;
}

uint32_t CreateOperatorAnnihilate(uint32_t const NFermion, uint16_t const Anindex, OperatorType * const OppAn)
{
  uint16_t oadd[1];
  InitNamedOperator(OppAn, "An");
  sprintf(OppAn->Name, "Annihilation operator on site %i", Anindex);
  OppAn->NFermion = NFermion;
  OppAn->NWord = Min(OppAn->NFermion, 1 + (OppAn->NFermion - 1) / 8);
  InitOperatorFixedLength(OppAn, 1, 8);
  oadd[0] = Anindex;
  RealOperatorAddLader(OppAn, oadd, 1, 1);
  return 0;
}


uint32_t CreateOperatorCreateList(uint32_t const NFermion, uint16_t const * Crindex, double const * VRe, double const * VIm, const uint32_t N, OperatorType * const OppCr)
{
  uint32_t i;
  uint16_t oadd[1];
  InitNamedOperator(OppCr, "Cr");
  OppCr->NFermion = NFermion;
  OppCr->NWord = Min(OppCr->NFermion, 1 + (OppCr->NFermion - 1) / 8);
  InitOperatorFixedLength(OppCr, 1, 8);
  for(i=0;i<N;i++)
  {
    oadd[0] = CrAnMask + Crindex[i];
    OperatorAddLader(OppCr, oadd, 1, VRe[i], VIm[i]);
  }
  return 0;
}

uint32_t CreateOperatorAnnihilateList(uint32_t const NFermion, uint16_t const * Anindex, double const * VRe, double const * VIm, const uint32_t N, OperatorType * const OppAn)
{
  uint32_t i;
  uint16_t oadd[1];
  InitNamedOperator(OppAn, "An");
  OppAn->NFermion = NFermion;
  OppAn->NWord = Min(OppAn->NFermion, 1 + (OppAn->NFermion - 1) / 8);
  InitOperatorFixedLength(OppAn, 1, 8);
  for(i=0;i<N;i++)
  {
    oadd[0] = Anindex[i];
    OperatorAddLader(OppAn, oadd, 1, VRe[i], VIm[i]);
  }
  return 0;
}

uint32_t CreateOperatorNumber(uint32_t const NFermion, uint16_t const Crindex, uint16_t const Anindex,
                              OperatorType * const OppN)
{
  uint16_t oadd[2];
  InitNamedOperator(OppN, "N");
  sprintf(OppN->Name, "Number operator a^{+}_(%i) a_(%i)", Crindex, Anindex);
  OppN->NFermion = NFermion;
  OppN->NWord = Min(OppN->NFermion, 1 + (OppN->NFermion - 1) / 8);
  InitOperatorFixedLength(OppN, 2, 8);
  oadd[1] = CrAnMask + Crindex;
  oadd[0] = Anindex;
  RealOperatorAddLader(OppN, oadd, 2, 1);
  return 0;
}


uint32_t CreateOperatorNumberShell(uint32_t const NFermion, uint32_t const Nelectron, const uint16_t * const index, OperatorType * const OppN)
{
  uint32_t i;
  uint16_t oadd[2];
  InitNamedOperator(OppN, "Number operator");
  OppN->NFermion = NFermion;
  OppN->NWord = Min(OppN->NFermion, 1 + (OppN->NFermion - 1) / 8);
  InitOperatorFixedLength(OppN, 2, 8);
  for(i=0;i<Nelectron;i++)
  {
    oadd[1] = CrAnMask + index[i];
    oadd[0] = index[i];
    RealOperatorAddLader(OppN, oadd, 2, 1);
  }
  return 0;
}

uint32_t CreateOperatorNumberShellSqr(uint32_t const NFermion, uint32_t const Nelectron, const uint16_t * const index, OperatorType * const OppNSqr)
{
  OperatorType temp;
  CreateOperatorNumberShell(NFermion, Nelectron, index, &temp);
  InitNamedOperator(OppNSqr, "Number operator squared");
  OppNSqr->NFermion = NFermion;
  OppNSqr->NWord = Min(OppNSqr->NFermion, 1 + (OppNSqr->NFermion - 1) / 8);
  InitOperatorFixedLength(OppNSqr, 2, 8);
  InitOperatorFixedLength(OppNSqr, 4, 8);
  OperatorMultiplyNormalOrder(OppNSqr, temp, temp,  &QDetNotOrderedFcFaBcBa);
  FreeOperator(&temp);
  return 0;
}

uint32_t CreateRandomTwoParticleOperator(uint32_t const NFermion, uint32_t const QComplex, OperatorType * const Opp)
{
  double ReRand,ImRand;
  uint32_t i,j;
  uint16_t oadd[2];
  InitOperator(Opp);
  Opp->NFermion = NFermion;
  Opp->NWord = Min(Opp->NFermion, 1 + (Opp->NFermion - 1) / 8);
  InitOperatorFixedLength(Opp, 2, 8);
  if(QComplex)
  {
    MakeRealOperatorFixedLengthComplex(Opp, 2);
    Opp->QComplex = 1;
    for(i=0;i<NFermion;i++)
    {
      for(j=0;j<NFermion;j++)
      {
        ReRand=0.5-rand()/(double)RAND_MAX;
        ImRand=0.5-rand()/(double)RAND_MAX;
        oadd[1] = CrAnMask + i;
        oadd[0] = j;
        ComplexOperatorAddLader(Opp, oadd, 2, ReRand, ImRand);
        oadd[1] = CrAnMask + j;
        oadd[0] = i;
        ComplexOperatorAddLader(Opp, oadd, 2, ReRand,-ImRand);
      }
    }
  }
  else
  {
    for(i=0;i<NFermion;i++)
    {
      for(j=0;j<NFermion;j++)
      {
        ReRand=0.5-rand()/(double)RAND_MAX;
        oadd[1] = CrAnMask + i;
        oadd[0] = j;
        RealOperatorAddLader(Opp, oadd, 2, ReRand);
        oadd[1] = CrAnMask + j;
        oadd[0] = i;
        RealOperatorAddLader(Opp, oadd, 2, ReRand);
      }
    }
  }
  RealOperatorAddLader(Opp, oadd, 0, -(double)NFermion);
  return 0;
}

