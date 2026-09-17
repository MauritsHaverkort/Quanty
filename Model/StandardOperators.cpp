uint16_t lmsigmaToIndex(int const l, int const m, int const sigma, uint16_t const *LindexUp, uint16_t const *LindexDn)
{
  if (sigma > 0) return LindexUp[l + m];
  return LindexDn[l + m];
}

/*  The one particle basis of a shell of angular momentum L can be written in several ways. Quanty labels
 *  them by ibasis
 *    0  SphericalHarmonicsY  - Y_{L,m}, m = -L ... L, each with spin dn and up
 *    1  TesseralHarmonicsZ   - the real combinations of Y_{L,m} and Y_{L,-m}, sin first, than cos
 *    2  KubicHarmonicsK      - the cubic harmonics, ordered according to irep a, b, e, t
 *    3  PyramidalHarmonicsP  - not implemented
 *    4  jjz                  - j = L-1/2 (2L states) followed by j = L+1/2 (2L+2 states), jz = -j ... j
 *
 *  All operators are defined in the spherical harmonics basis and rotated from there. The matrix
 *
 *    C[i][m] = <new_i|Y_m>
 *
 *  is what YToZMatrix, YToKMatrix and YToJJzMatrix return. A one particle operator with matrix A in the
 *  spherical harmonics basis has in the new basis the matrix
 *
 *    A' = C^* A C^T
 *
 *  which is what RotateOperatorNormalOrder produces when it is handed C itself. Its convention is
 *
 *    a_m^{dag} (Y basis) = \sum_i C[i][m]^* a_i^{dag} (new basis)
 *    a_m       (Y basis) = \sum_i C[i][m]   a_i       (new basis)
 *
 *  C is complex for the Z and K basis and real for jjz. It is unitary, so C^{-1} = C^{dag}, but C^T is
 *  not C^{dag} and the transpose may never be substituted for the adjoint. Handing
 *  RotateOperatorNormalOrder C^{dag} instead of C rotates the other way, from the new basis back to the
 *  spherical harmonics. That mistake produces a well formed operator and no warning, so it is only
 *  caught by a test: rotating Jz to the jjz basis with C gives a diagonal Jz and with C^{dag} does not.
 */

/*  ShellBasisSlots fills slot with the indices in the NFermion particle basis of the 2(2L+1) orbitals of
 *  one shell, ordered as the rows of the basis rotation matrix C.
 *  Input:
 *    L           - The orbital angular momentum of the shell
 *    ibasis      - The basis the indices are given in
 *    *LindexUp   - for ibasis 0, 1, 2 and 3 an array of length 2L+1 with the spin up indices, m = -L ... L
 *                  for ibasis 4          an array of length 2L   with the j = L-1/2 indices, jz = -j ... j
 *    *LindexDn   - for ibasis 0, 1, 2 and 3 an array of length 2L+1 with the spin dn indices, m = -L ... L
 *                  for ibasis 4          an array of length 2L+2 with the j = L+1/2 indices, jz = -j ... j
 *    *slot       - An array of length 2(2L+1) that is filled
 */
void ShellBasisSlots(uint32_t const L, uint32_t const ibasis, uint16_t const * const LindexUp, uint16_t const * const LindexDn, uint16_t * const slot)
{
  uint32_t m;
  if (ibasis == 4)
  {
    for (m = 0; m < 2 * L; m++)     slot[m]          = LindexUp[m];
    for (m = 0; m < 2 * L + 2; m++) slot[m + 2 * L]  = LindexDn[m];
  }
  else
  {
    for (m = 0; m < 2 * L + 1; m++)
    {
      slot[2 * m]     = LindexDn[m];
      slot[2 * m + 1] = LindexUp[m];
    }
  }
}

/*  ShellBasisYIndices fills, for one shell, both the slot array of ShellBasisSlots and the two index
 *  arrays the operator is to be created on in the spherical harmonics basis. These are the same NFermion
 *  slots, read as dn, up, dn, up, ... from m = -L to m = L.
 *  Input:
 *    L, ibasis, *LindexUp, *LindexDn - As for ShellBasisSlots
 *    *slot       - An array of length 2(2L+1) that is filled
 *    *YindexUp   - An array of length 2L+1 that is filled with the spin up spherical harmonics indices
 *    *YindexDn   - An array of length 2L+1 that is filled with the spin dn spherical harmonics indices
 */
void ShellBasisYIndices(uint32_t const L, uint32_t const ibasis, uint16_t const * const LindexUp, uint16_t const * const LindexDn, uint16_t * const slot, uint16_t * const YindexUp, uint16_t * const YindexDn)
{
  uint32_t m;
  ShellBasisSlots(L, ibasis, LindexUp, LindexDn, slot);
  for (m = 0; m < 2 * L + 1; m++)
  {
    YindexDn[m] = slot[2 * m];
    YindexUp[m] = slot[2 * m + 1];
  }
}

/*  AddShellBasisRotationToMatrix writes into U the block that takes one shell of angular momentum L from
 *  the spherical harmonics basis to the basis ibasis. U is left untouched outside of the shell.
 *  Input:
 *    L           - The orbital angular momentum of the shell
 *    ibasis      - The basis to rotate to, 1 (Z), 2 (K) or 4 (jjz)
 *    *slot       - An array of length 2(2L+1) as filled by ShellBasisSlots
 *    U           - An initialized NFermion x NFermion matrix with QComplex set to 1
 */
uint32_t AddShellBasisRotationToMatrix(uint32_t const L, uint32_t const ibasis, uint16_t const * const slot, MatrixType * const U)
{
  uint32_t i, j;
  MatrixType C;
  switch (ibasis)
  {
    case 1:
      C = YToZMatrix((int)(2 * (2 * L + 1)), true);
      break;
    case 2:
      C = YToKMatrix((int)(2 * (2 * L + 1)), true);
      break;
    case 4:
      C = YToJJzMatrix((int)(2 * (2 * L + 1)), true);
      break;
    default:
      printf("AddShellBasisRotationToMatrix called with ibasis %d, only 1 (TesseralHarmonicsZ), 2 (KubicHarmonicsK) and 4 (jjz) can be reached by a rotation\n", ibasis);
      return 1;
  }
  if (isNullMatrix(C))
  {
    printf("AddShellBasisRotationToMatrix could not create the basis rotation matrix for L = %d and ibasis = %d\n", L, ibasis);
    return 1;
  }
  for (i = 0; i < 2 * (2 * L + 1); i++)
  {
    for (j = 0; j < 2 * (2 * L + 1); j++)
    {
      U->VRe[slot[i]][slot[j]] = C.VRe[i][j];
      U->VIm[slot[i]][slot[j]] = (C.QComplex) ? C.VIm[i][j] : 0.0;
    }
  }
  FreeMatrix(&C);
  return 0;
}

/*  RotateOperatorToBasis rotates the operator Opp, defined on the shells whose slots are listed, from the
 *  spherical harmonics basis to the basis ibasis. Opp is replaced by the rotated operator.
 *  Input:
 *    NFermion    - The number of Fermions in the basis
 *    NShell      - The number of shells the operator acts on
 *    *L          - An array of length NShell with the angular momentum of each shell
 *    ibasis      - The basis to rotate to, 1 (Z), 2 (K) or 4 (jjz)
 *    **slot      - An array of length NShell of arrays as filled by ShellBasisSlots
 *    Opp         - The operator in the spherical harmonics basis, replaced by the rotated operator
 */
/*  RotateOperatorWithMatrix replaces Opp by the operator rotated with the matrix U. */
uint32_t RotateOperatorWithMatrix(MatrixType const U, OperatorType * const Opp)
{
  OperatorType temp;
  if (RotateOperatorNormalOrder(&temp, *Opp, U, &QDetNotOrderedFcFaBcBa))
  {
    printf("RotateOperatorNormalOrder failed in RotateOperatorWithMatrix\n");
    return 1;
  }
  OperatorChop(&temp, 100 * DBL_EPSILON);
  FreeOperator(Opp);
  OperatorCopy(Opp, temp);
  FreeOperator(&temp);
  return 0;
}

/*  CreateBasisRotationMatrix creates the NFermion x NFermion matrix that is the identity outside of the
 *  NShell shells listed and the basis rotation C of AddShellBasisRotationToMatrix on each of them.
 *  Release with FreeMatrix.
 */
uint32_t CreateBasisRotationMatrix(uint32_t const NFermion, uint32_t const NShell, uint32_t const * const L, uint32_t const ibasis, uint16_t const * const * const slot, MatrixType * const U)
{
  uint32_t i;
  U->QComplex = 1;
  U->Ni = NFermion;
  U->Nj = NFermion;
  if (InitIdentityMatrix(U))
  {
    printf("Could not allocate the basis rotation matrix in CreateBasisRotationMatrix\n");
    return 1;
  }
  for (i = 0; i < NShell; i++)
  {
    if (AddShellBasisRotationToMatrix(L[i], ibasis, slot[i], U))
    {
      FreeMatrix(U);
      return 1;
    }
  }
  return 0;
}

uint32_t RotateOperatorToBasis(uint32_t const NFermion, uint32_t const NShell, uint32_t const * const L, uint32_t const ibasis, uint16_t const * const * const slot, OperatorType * const Opp)
{
  uint32_t ret;
  MatrixType U;
  if (CreateBasisRotationMatrix(NFermion, NShell, L, ibasis, slot, &U)) return 1;
  ret = RotateOperatorWithMatrix(U, Opp);
  FreeMatrix(&U);
  return ret;
}

/*  ShellBasisYIndicesList allocates and fills, for NShell shells, the slot arrays of ShellBasisSlots and
 *  the index arrays the operator is to be created on in the spherical harmonics basis. On failure nothing
 *  is left allocated. Release with FreeShellBasisYIndicesList.
 *  Input:
 *    NShell      - The number of shells
 *    *L          - An array of length NShell with the angular momentum of each shell
 *    ibasis      - The basis the indices are given in
 *    **LindexUp  - An array of length NShell of index arrays, as for ShellBasisSlots
 *    **LindexDn  - An array of length NShell of index arrays, as for ShellBasisSlots
 *    **slot      - An array of length NShell that is filled with allocated arrays of length 2(2L+1)
 *    **YindexUp  - An array of length NShell that is filled with allocated arrays of length 2L+1
 *    **YindexDn  - An array of length NShell that is filled with allocated arrays of length 2L+1
 */
uint32_t ShellBasisYIndicesList(uint32_t const NShell, uint32_t const * const L, uint32_t const ibasis, uint16_t * const * const LindexUp, uint16_t * const * const LindexDn, uint16_t ** const slot, uint16_t ** const YindexUp, uint16_t ** const YindexDn)
{
  uint32_t i;
  for (i = 0; i < NShell; i++)
  {
    slot[i]     = (uint16_t *) malloc(2 * (2 * L[i] + 1) * sizeof(uint16_t));
    YindexUp[i] = (uint16_t *) malloc((2 * L[i] + 1) * sizeof(uint16_t));
    YindexDn[i] = (uint16_t *) malloc((2 * L[i] + 1) * sizeof(uint16_t));
    if ((slot[i] == NULL) || (YindexUp[i] == NULL) || (YindexDn[i] == NULL))
    {
      printf("Could not allocate memory in ShellBasisYIndicesList\n");
      FreeShellBasisYIndicesList(i + 1, slot, YindexUp, YindexDn);
      return 1;
    }
    ShellBasisYIndices(L[i], ibasis, LindexUp[i], LindexDn[i], slot[i], YindexUp[i], YindexDn[i]);
  }
  return 0;
}

void FreeShellBasisYIndicesList(uint32_t const NShell, uint16_t ** const slot, uint16_t ** const YindexUp, uint16_t ** const YindexDn)
{
  uint32_t i;
  for (i = 0; i < NShell; i++)
  {
    free(slot[i]);     slot[i]     = NULL;
    free(YindexUp[i]); YindexUp[i] = NULL;
    free(YindexDn[i]); YindexDn[i] = NULL;
  }
}

/*  FinishOperatorInShellBasis rotates an operator that has just been created in the spherical harmonics
 *  basis to the basis ibasis. ret is the return value of the routine that created the operator. If the
 *  basis can not be reached the operator is left in the spherical harmonics basis and a non zero value is
 *  returned.
 */
uint32_t FinishOperatorInShellBasis(uint32_t const NFermion, uint32_t const NShell, uint32_t const * const L, uint32_t const ibasis, uint16_t const * const * const slot, uint32_t const ret, OperatorType * const Opp)
{
  if (ret != 0) return ret;
  if (ibasis == 0) return 0;
  if (ibasis == 3)
  {
    printf("The PyramidalHarmonicsP basis is not implemented, the operator is returned in the SphericalHarmonicsY basis\n");
    return 1;
  }
  return RotateOperatorToBasis(NFermion, NShell, L, ibasis, slot, Opp);
}

/*  CreateOperatorInShellBasis creates a single shell operator in the basis ibasis. The operator is first
 *  created in the spherical harmonics basis by CreateOperatorY and then rotated.
 *  Input:
 *    NFermion        - The number of Fermions in the basis
 *    L               - The orbital angular momentum of the shell
 *    ibasis          - The basis to create the operator in
 *    *LindexUp       - As for ShellBasisSlots
 *    *LindexDn       - As for ShellBasisSlots
 *    CreateOperatorY - The function creating the operator in the spherical harmonics basis
 *    The operator Opp is not initialized on input
 */
uint32_t CreateOperatorInShellBasis(uint32_t const NFermion, uint32_t const L, uint32_t const ibasis, uint16_t * const LindexUp, uint16_t * const LindexDn, uint32_t (*CreateOperatorY) (uint32_t, uint32_t, uint16_t *, uint16_t *, OperatorType * const), OperatorType * const Opp)
{
  uint32_t ret;
  uint16_t *slot[1], *YindexUp[1], *YindexDn[1];
  if (ShellBasisYIndicesList(1, &L, ibasis, &LindexUp, &LindexDn, slot, YindexUp, YindexDn)) return 1;
  ret = CreateOperatorY(NFermion, L, YindexUp[0], YindexDn[0], Opp);
  ret = FinishOperatorInShellBasis(NFermion, 1, &L, ibasis, (uint16_t const * const *) slot, ret, Opp);
  FreeShellBasisYIndicesList(1, slot, YindexUp, YindexDn);
  return ret;
}

/*  CreateOperatorCFInShellBasis creates the crystal-field operator between two shells in the basis ibasis.
 *  Input:
 *    NFermion      - The number of Fermions in the basis
 *    ibasis        - The basis to create the operator in
 *    LCr, LAn      - The orbital angular momentum of the shell that is created and the one annihilated
 *    *LCrindex**   - As for ShellBasisSlots, for the shell that is created
 *    *LAnindex**   - As for ShellBasisSlots, for the shell that is annihilated
 *    Akm           - The expansion of the potential on renormalized spherical harmonics
 *    The operator CF is not initialized on input
 */
uint32_t CreateOperatorCFInShellBasis(uint32_t const NFermion, uint32_t const ibasis, uint32_t const LCr, uint32_t const LAn, uint16_t * const LCrindexUp, uint16_t * const LCrindexDn, uint16_t * const LAnindexUp, uint16_t * const LAnindexDn, HarmonicsExpansion const Akm, OperatorType * const CF)
{
  uint32_t ret;
  uint32_t L[2];
  uint16_t *Lup[2], *Ldn[2], *slot[2], *Yup[2], *Ydn[2];
  L[0] = LCr; Lup[0] = LCrindexUp; Ldn[0] = LCrindexDn;
  L[1] = LAn; Lup[1] = LAnindexUp; Ldn[1] = LAnindexDn;
  if (ShellBasisYIndicesList(2, L, ibasis, Lup, Ldn, slot, Yup, Ydn)) return 1;
  ret = CreateOperatorCF(NFermion, LCr, LAn, Yup[0], Ydn[0], Yup[1], Ydn[1], Akm, CF);
  ret = FinishOperatorInShellBasis(NFermion, 2, L, ibasis, (uint16_t const * const *) slot, ret, CF);
  FreeShellBasisYIndicesList(2, slot, Yup, Ydn);
  return ret;
}

/*  CreateOperatorCoulomb***InBasis create the coulomb operator in the basis ibasis. The operator is first
 *  created in the spherical harmonics basis, where it is parametrized by the Slater integrals F, G and R,
 *  and then rotated. The arguments are those of the routine wrapped, with ibasis inserted after L.
 */
uint32_t CreateOperatorCoulombSingleShellInBasis(uint32_t const NFermion, uint32_t const L, uint32_t const ibasis, uint16_t * const LindexUp, uint16_t * const LindexDn, double const *F, OperatorType * const U)
{
  uint32_t ret;
  uint16_t *slot[1], *Yup[1], *Ydn[1];
  if (ShellBasisYIndicesList(1, &L, ibasis, &LindexUp, &LindexDn, slot, Yup, Ydn)) return 1;
  ret = CreateOperatorCoulombSingleShell(NFermion, L, Yup[0], Ydn[0], F, U);
  ret = FinishOperatorInShellBasis(NFermion, 1, &L, ibasis, (uint16_t const * const *) slot, ret, U);
  FreeShellBasisYIndicesList(1, slot, Yup, Ydn);
  return ret;
}

uint32_t CreateOperatorCoulombTwoShellsConstantOccupationInBasis(uint32_t const NFermion, uint32_t const L1, uint32_t const L2, uint32_t const ibasis, uint16_t * const Lindex1Up, uint16_t * const Lindex1Dn, uint16_t * const Lindex2Up, uint16_t * const Lindex2Dn, double const *F, double const *G, OperatorType * const U)
{
  uint32_t ret;
  uint32_t L[2];
  uint16_t *Lup[2], *Ldn[2], *slot[2], *Yup[2], *Ydn[2];
  L[0] = L1; Lup[0] = Lindex1Up; Ldn[0] = Lindex1Dn;
  L[1] = L2; Lup[1] = Lindex2Up; Ldn[1] = Lindex2Dn;
  if (ShellBasisYIndicesList(2, L, ibasis, Lup, Ldn, slot, Yup, Ydn)) return 1;
  ret = CreateOperatorCoulombTwoShellsConstantOccupation(NFermion, L1, L2, Yup[0], Ydn[0], Yup[1], Ydn[1], F, G, U);
  ret = FinishOperatorInShellBasis(NFermion, 2, L, ibasis, (uint16_t const * const *) slot, ret, U);
  FreeShellBasisYIndicesList(2, slot, Yup, Ydn);
  return ret;
}

uint32_t CreateOperatorCoulombFourShellsInBasis(uint32_t const NFermion, uint32_t const L1, uint32_t const L2, uint32_t const L3, uint32_t const L4, uint32_t const ibasis, uint16_t * const Lindex1Up, uint16_t * const Lindex1Dn, uint16_t * const Lindex2Up, uint16_t * const Lindex2Dn, uint16_t * const Lindex3Up, uint16_t * const Lindex3Dn, uint16_t * const Lindex4Up, uint16_t * const Lindex4Dn, double const *R, OperatorType * const U)
{
  uint32_t ret;
  uint32_t L[4];
  uint16_t *Lup[4], *Ldn[4], *slot[4], *Yup[4], *Ydn[4];
  L[0] = L1; Lup[0] = Lindex1Up; Ldn[0] = Lindex1Dn;
  L[1] = L2; Lup[1] = Lindex2Up; Ldn[1] = Lindex2Dn;
  L[2] = L3; Lup[2] = Lindex3Up; Ldn[2] = Lindex3Dn;
  L[3] = L4; Lup[3] = Lindex4Up; Ldn[3] = Lindex4Dn;
  if (ShellBasisYIndicesList(4, L, ibasis, Lup, Ldn, slot, Yup, Ydn)) return 1;
  ret = CreateOperatorCoulombFourShells(NFermion, L1, L2, L3, L4, Yup[0], Ydn[0], Yup[1], Ydn[1], Yup[2], Ydn[2], Yup[3], Ydn[3], R, U);
  ret = FinishOperatorInShellBasis(NFermion, 4, L, ibasis, (uint16_t const * const *) slot, ret, U);
  FreeShellBasisYIndicesList(4, slot, Yup, Ydn);
  return ret;
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
uint32_t CreateOperatorCoulombSingleShell(uint32_t const NFermion, uint32_t const L, uint16_t const *LindexUp, uint16_t const *LindexDn, double const *F, OperatorType * const U)
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

uint32_t CreateOperatorCoulombTwoShellsConstantOccupation(uint32_t const NFermion, uint32_t const L1, uint32_t const L2, uint16_t const *Lindex1Up, uint16_t const *Lindex1Dn, uint16_t const *Lindex2Up, uint16_t const *Lindex2Dn, double const *F, double const *G, OperatorType * const U)
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

uint32_t CreateOperatorCoulombFourShells(uint32_t const NFermion, uint32_t const L1, uint32_t const L2, uint32_t const L3, uint32_t const L4, uint16_t const *Lindex1Up, uint16_t const *Lindex1Dn, uint16_t const *Lindex2Up, uint16_t const *Lindex2Dn, uint16_t const *Lindex3Up, uint16_t const *Lindex3Dn, uint16_t const *Lindex4Up, uint16_t const *Lindex4Dn, double const *R, OperatorType * const U)
{
  uint16_t oadd[4];
  double val;
  int    ik, k, l1, l2, l3, l4, m1, m2, m3, m4, s1, s2;
  l1 = (int)L1;
  l2 = (int)L2;
  l3 = (int)L3;
  l4 = (int)L4;
  InitNamedOperator(U, "Coulomb Operator");
  U->NFermion = NFermion;
  U->NWord = Min(U->NFermion, 1 + (U->NFermion - 1) / 8);
  InitOperatorFixedLength(U, 2, 8);
  InitOperatorFixedLength(U, 4, 8);
  ik=-1;
  for (k=Max(Abs(l1-l3),Abs(l2-l4)); k<=Min(l1+l3,l2+l4); k+=2)
  {
    ik++;
    for (m1 = -l1; m1 <= l1; m1++)
    {
      for (m2 = -l2; m2 <= l2; m2++)
      {
        for (m3 = Max(-l3, m1 + m2 - l4); m3 <= Min(l3, m1 + m2 + l4); m3++)
        {
          /* -l4 <= m4 <= l4
           * -l4 <= m1+m2-m3 <= l4
           *  thus m3-l4 <= m1 + m2 <= m3+l4
           */
          m4 = m1 + m2 - m3;
          
          val = -0.5*SlaterCoefficientC(k, l1, m1, l3, m3) * SlaterCoefficientC(k, l4, m4, l2, m2);
          if (Abs(val) > DBL_EPSILON)
          {
            for (s1 = 0; s1 <= 1; s1++)
            {
              for (s2 = 0; s2 <= 1; s2++)
              {
                oadd[3] = CrAnMask + (s1 == 0 ? Lindex1Dn[l1 + m1] : Lindex1Up[l1 + m1]);
                oadd[2] = CrAnMask + (s2 == 0 ? Lindex2Dn[l2 + m2] : Lindex2Up[l2 + m2]);
                oadd[1] = (s1 == 0 ? Lindex3Dn[l3 + m3] : Lindex3Up[l3 + m3]);
                oadd[0] = (s2 == 0 ? Lindex4Dn[l4 + m4] : Lindex4Up[l4 + m4]);
                RealOperatorAddLaderNormalOrder(U, oadd, 4, val * R[ik], &QDetNotOrderedFcFaBcBa);
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


uint32_t CreateOperatorCoulombSingleRelativisticShell(uint32_t const NFermion, int const kappa, uint16_t const *Jindex, double const *F, OperatorType * const U)
{
  uint16_t oadd[4];
  double val;
  int const twoj = 2*Abs(kappa)-1;
  int const l=(twoj+Sgn(kappa))/2; // l of the large component, l of the small component is j-Sgn(kappa)/2
  int k,twojz1,twojz2,twojz3,twojz4;
  InitNamedOperator(U, "Coulomb Operator");
  U->NFermion = NFermion;
  U->NWord = Min(U->NFermion, 1 + (U->NFermion - 1) / 8);
  InitOperatorFixedLength(U, 2, 8);
  InitOperatorFixedLength(U, 4, 8);
  for(k=0;k<twoj+1;k+=2)
  {
    for(twojz1 = -twoj; twojz1 <= twoj; twojz1+=2)
    {
      for(twojz2 = -twoj; twojz2 <= twoj; twojz2+=2)
      {
        for (twojz3 = Max(-twoj, twojz1 + twojz2 - twoj); twojz3 <= Min(twoj, twojz1 + twojz2 + twoj); twojz3+=2)
        {
          /* -l >= m1+m2-m3 >= l
           *  thus m3-l >= m1 + m2 >= m3+l
           */
          twojz4 = twojz1 + twojz2 - twojz3;
          val = -0.5 * SlaterCoefficientCjcoupling(k, l, twoj, twojz1, l, twoj, twojz3) * SlaterCoefficientCjcoupling(k, l, twoj, twojz4, l, twoj, twojz2);
          if (Abs(val) > DBL_EPSILON)
          {
            oadd[3] = CrAnMask + Jindex[(twoj+twojz1)/2];
            oadd[2] = CrAnMask + Jindex[(twoj+twojz2)/2];
            oadd[1] =            Jindex[(twoj+twojz3)/2];
            oadd[0] =            Jindex[(twoj+twojz4)/2];
            RealOperatorAddLaderNormalOrder(U, oadd, 4, val * F[k / 2], &QDetNotOrderedFcFaBcBa);
          }
        }
      }
    }
  }
  OperatorChop(U, 100 * DBL_EPSILON);
  return 0;
}

uint32_t CreateOperatorCoulombTwoRelativisticShellsConstantOccupation(uint32_t const NFermion, int const kappa1, int const kappa2, uint16_t const *Jindex1, uint16_t const *Jindex2, double const *F, double const *G, OperatorType * const U)
{
  uint16_t oadd[4];
  double val;
  int const twoj1 = 2*Abs(kappa1)-1;
  int const twoj2 = 2*Abs(kappa2)-1;
  int const l1=(twoj1+Sgn(kappa1))/2; // l of the large component, l of the small component is j-Sgn(kappa)/2
  int const l2=(twoj2+Sgn(kappa2))/2; // l of the large component, l of the small component is j-Sgn(kappa)/2
  int ik,k,twojz1,twojz2,twojz3,twojz4;
  InitNamedOperator(U, "Coulomb Operator");
  U->NFermion = NFermion;
  U->NWord = Min(U->NFermion, 1 + (U->NFermion - 1) / 8);
  InitOperatorFixedLength(U, 2, 8);
  InitOperatorFixedLength(U, 4, 8);
  /* case 1 jz1 and jz3 belong to jz1, jz2 and jz4 belong to j2 */
  for (k = 0; k <= (Min(twoj1, twoj2) + 1); k += 2)
  {
    for(twojz1 = -twoj1; twojz1 <= twoj1; twojz1+=2)
    {
      for(twojz2 = -twoj2; twojz2 <= twoj2; twojz2+=2)
      {
        for (twojz3 = Max(-twoj1, twojz1 + twojz2 - twoj2); twojz3 <= Min(twoj1, twojz1 + twojz2 + twoj2); twojz3+=2)
        {
          /* -j2 <= jz4 <= j2
           * -j2 <= jz1+jz2-jz3 <= j2
           * thus
           * -j2-jz1-jz2 <= -jz3 <= j2-jz1-jz2
           * or
           *  jz1+jz2-j2 <= jz3 <= j2+jz1+jz2
           */
          twojz4 = twojz1 + twojz2 - twojz3;
          val = -SlaterCoefficientCjcoupling(k, l1, twoj1, twojz1, l1, twoj1, twojz3) * SlaterCoefficientCjcoupling(k, l2, twoj2, twojz4, l2, twoj2, twojz2);
          if (Abs(val) > DBL_EPSILON)
          {
            oadd[3] = CrAnMask + Jindex1[(twoj1+twojz1)/2];
            oadd[2] = CrAnMask + Jindex2[(twoj2+twojz2)/2];
            oadd[1] =            Jindex1[(twoj1+twojz3)/2];
            oadd[0] =            Jindex2[(twoj2+twojz4)/2];
            RealOperatorAddLaderNormalOrder(U, oadd, 4, val * F[k / 2], &QDetNotOrderedFcFaBcBa);
          }
        }
      }
    }
  }
  /* case 2 m1 and m4 belong to l1/j1, m2 and m3 belong to l2/j2 */
  ik=-1;
  for (k = Abs(l1 - l2); k <= (twoj1 + twoj2)/2; k += 2)
  {
    if(k<Abs(twoj1 - twoj2)/2) continue;
    ik++;
    for (twojz1 = -twoj1; twojz1 <= twoj1; twojz1+=2)
    {
      for (twojz2 = -twoj2; twojz2 <= twoj2; twojz2+=2)
      {
        for (twojz3 = Max(-twoj2, twojz1 + twojz2 - twoj1); twojz3 <= Min(twoj2, twojz1 + twojz2 + twoj1); twojz3+=2)
        {
         /* -j1 <= jz4 <= j1
          * -j1 <= jz1+jz2-jz3 <= j1
          * thus
          * -j1-jz1-jz2 <= -jz3 <= j1-jz1-jz2
          * or
          *  jz1+jz2-j1 <= jz3 <= j1+jz1+jz2
          */
          twojz4 = twojz1 + twojz2 - twojz3;
          val = -SlaterCoefficientCjcoupling(k, l1, twoj1, twojz1, l2, twoj2, twojz3) * SlaterCoefficientCjcoupling(k, l1, twoj1, twojz4, l2, twoj2, twojz2);
          if (Abs(val) > DBL_EPSILON)
          {
            oadd[3] = CrAnMask + Jindex1[(twoj1+twojz1)/2];
            oadd[2] = CrAnMask + Jindex2[(twoj2+twojz2)/2];
            oadd[1] =            Jindex2[(twoj2+twojz3)/2];
            oadd[0] =            Jindex1[(twoj1+twojz4)/2];
            RealOperatorAddLaderNormalOrder(U, oadd, 4, val * G[ik], &QDetNotOrderedFcFaBcBa);
          }
        }
      }
    }
  }
  OperatorChop(U, 100 * DBL_EPSILON);
  return 0;
}

uint32_t CreateOperatorCoulombFourRelativisticShells(uint32_t const NFermion, int const kappa1, int const kappa2, int const kappa3, int const kappa4, uint16_t const *Jindex1, uint16_t const *Jindex2, uint16_t const *Jindex3, uint16_t const *Jindex4, double const *R, OperatorType * const U)
{
  uint16_t oadd[4];
  double val;
  int const twoj1 = 2*Abs(kappa1)-1;
  int const twoj2 = 2*Abs(kappa2)-1;
  int const twoj3 = 2*Abs(kappa3)-1;
  int const twoj4 = 2*Abs(kappa4)-1;
  int const l1=(twoj1+Sgn(kappa1))/2; // l of the large component, l of the small component is j-Sgn(kappa)/2
  int const l2=(twoj2+Sgn(kappa2))/2; // l of the large component, l of the small component is j-Sgn(kappa)/2
  int const l3=(twoj3+Sgn(kappa3))/2; // l of the large component, l of the small component is j-Sgn(kappa)/2
  int const l4=(twoj4+Sgn(kappa4))/2; // l of the large component, l of the small component is j-Sgn(kappa)/2
  int ik,k,twojz1,twojz2,twojz3,twojz4;
  InitNamedOperator(U, "Coulomb Operator");
  U->NFermion = NFermion;
  U->NWord = Min(U->NFermion, 1 + (U->NFermion - 1) / 8);
  InitOperatorFixedLength(U, 2, 8);
  InitOperatorFixedLength(U, 4, 8);
  ik=-1;
  
  //printf("l=  %d %d %d %d\n",l1,l2,l3,l4);
  //printf("2j= %d %d %d %d\n",twoj1,twoj2,twoj3,twoj4);
  //printf("k range: %d %d\n" ,Max(Abs(l1-l3),Abs(l2-l4)), Min(twoj1+twoj3,twoj2+twoj4));
  //printf("k filter: remove 2k<%d\n",Max(Abs(twoj1-twoj3),Abs(twoj2-twoj4)));
  //printf("k filter: remove 2k>%d\n",Min(twoj1+twoj3,twoj2+twoj4));
  
  for (k=Max(Abs(l1-l3),Abs(l2-l4)); k<=Min(l1+l3,l2+l4); k+=2)
  {
    if ((2*k)<Max(Abs(twoj1-twoj3),Abs(twoj2-twoj4))) continue;
    if ((2*k)>Min(twoj1+twoj3,twoj2+twoj4)) continue;
    ik++;
    for (twojz1 = -twoj1; twojz1 <= twoj1; twojz1+=2)
    {
      for (twojz2 = -twoj2; twojz2 <= twoj2; twojz2+=2)
      {
        for (twojz3 = Max(-twoj3, twojz1 + twojz2 - twoj4); twojz3 <= Min(twoj3, twojz1 + twojz2 + twoj4); twojz3+=2)
        {
          /* -j4 <= jz4 <= j4
           * -j4 <= jz1+jz2-jz3 <= j4
           * thus
           * -j4-jz1-jz2 <= -jz3 <= j4-jz1-jz2
           * or
           *  jz1+jz2-j4 <= jz3 <= j44+jz1+jz2
           */
          twojz4 = twojz1 + twojz2 - twojz3;
          val = -0.5*SlaterCoefficientCjcoupling(k, l1, twoj1, twojz1, l3, twoj3, twojz3) * SlaterCoefficientCjcoupling(k, l4, twoj4, twojz4, l2, twoj2, twojz2);
          if (Abs(val) > DBL_EPSILON)
          {
            oadd[3] = CrAnMask + Jindex1[(twoj1+twojz1)/2];
            oadd[2] = CrAnMask + Jindex2[(twoj2+twojz2)/2];
            oadd[1] =            Jindex3[(twoj3+twojz3)/2];
            oadd[0] =            Jindex4[(twoj4+twojz4)/2];
            RealOperatorAddLaderNormalOrder(U, oadd, 4, val * R[ik], &QDetNotOrderedFcFaBcBa);
          }
        }
      }
    }
  }
  OperatorChop(U, 100 * DBL_EPSILON);
  return 0;
}


uint32_t CreateOperatorCoulombSingleRelativisticShellOld(uint32_t const NFermion, int const kappa, uint16_t const *Jindex, double const *Flarge, double const *Fsmall, double const *Fmix, OperatorType * const U)
{
  uint16_t oadd[4];
  double val;
  int const twoj=abs(2*kappa)-1;
  int const llarge=(twoj+Sgn(kappa))/2;
  int const lsmall=(twoj-Sgn(kappa))/2;
  int k,twojz1,twojz2,twojz3,twojz4;
  InitNamedOperator(U, "Coulomb Operator");
  U->NFermion = NFermion;
  U->NWord = Min(U->NFermion, 1 + (U->NFermion - 1) / 8);
  InitOperatorFixedLength(U, 2, 8);
  InitOperatorFixedLength(U, 4, 8);
  for(k=0;k<2*(llarge+1);k+=2)
  {
    for(twojz1 = -twoj; twojz1 <= twoj; twojz1+=2)
    {
      for(twojz2 = -twoj; twojz2 <= twoj; twojz2+=2)
      {
        for (twojz3 = Max(-twoj, twojz1 + twojz2 - twoj); twojz3 <= Min(twoj, twojz1 + twojz2 + twoj); twojz3+=2)
        {
          /* -l >= m1+m2-m3 >= l
           *  thus m3-l >= m1 + m2 >= m3+l
           */
          
          
          //printf("Ck[k=%i;l=%i,2j=%i,2jz=%2i;l=%i,2j=%i,2jz=%2i]=%7.4f  ",k,llarge,twoj,twojz1,llarge,twoj,twojz3,SlaterCoefficientCjcoupling(k, llarge, twoj, twojz1, llarge, twoj, twojz3));
          //printf("Ck[k=%i;l=%i,2j=%i,2jz=%2i;l=%i,2j=%i,2jz=%2i]=%7.4f\n",k,lsmall,twoj,twojz1,lsmall,twoj,twojz3,SlaterCoefficientCjcoupling(k, lsmall, twoj, twojz1, lsmall, twoj, twojz3));
          
          
          
          twojz4 = twojz1 + twojz2 - twojz3;
          val = -0.5 * SlaterCoefficientCjcoupling(k, llarge, twoj, twojz1, llarge, twoj, twojz3) * SlaterCoefficientCjcoupling(k, llarge, twoj, twojz4, llarge, twoj, twojz2);
          if (Abs(val) > DBL_EPSILON)
          {
            oadd[3] = CrAnMask + Jindex[(twoj+twojz1)/2];
            oadd[2] = CrAnMask + Jindex[(twoj+twojz2)/2];
            oadd[1] =            Jindex[(twoj+twojz3)/2];
            oadd[0] =            Jindex[(twoj+twojz4)/2];
            RealOperatorAddLaderNormalOrder(U, oadd, 4, val * Flarge[k / 2], &QDetNotOrderedFcFaBcBa);
          }
        }
      }
    }
  }
  for(k=0;k<2*(lsmall+1);k+=2)
  {
    for(twojz1 = -twoj; twojz1 <= twoj; twojz1+=2)
    {
      for(twojz2 = -twoj; twojz2 <= twoj; twojz2+=2)
      {
        for (twojz3 = Max(-twoj, twojz1 + twojz2 - twoj); twojz3 <= Min(twoj, twojz1 + twojz2 + twoj); twojz3+=2)
        {
          /* -l >= m1+m2-m3 >= l
           *  thus m3-l >= m1 + m2 >= m3+l
           */
          twojz4 = twojz1 + twojz2 - twojz3;
          val = -0.5 * SlaterCoefficientCjcoupling(k, lsmall, twoj, twojz1, lsmall, twoj, twojz3) * SlaterCoefficientCjcoupling(k, lsmall, twoj, twojz4, lsmall, twoj, twojz2);
          if (Abs(val) > DBL_EPSILON)
          {
            oadd[3] = CrAnMask + Jindex[(twoj+twojz1)/2];
            oadd[2] = CrAnMask + Jindex[(twoj+twojz2)/2];
            oadd[1] =            Jindex[(twoj+twojz3)/2];
            oadd[0] =            Jindex[(twoj+twojz4)/2];
            RealOperatorAddLaderNormalOrder(U, oadd, 4, val * Fsmall[k / 2], &QDetNotOrderedFcFaBcBa);
          }
        }
      }
    }
  }
  for(k=0;k<2*(Min(llarge,lsmall)+1);k+=2)
  {
    for(twojz1 = -twoj; twojz1 <= twoj; twojz1+=2)
    {
      for(twojz2 = -twoj; twojz2 <= twoj; twojz2+=2)
      {
        for (twojz3 = Max(-twoj, twojz1 + twojz2 - twoj); twojz3 <= Min(twoj, twojz1 + twojz2 + twoj); twojz3+=2)
        {
          /* -l >= m1+m2-m3 >= l
           *  thus m3-l >= m1 + m2 >= m3+l
           */
          twojz4 = twojz1 + twojz2 - twojz3;
          val = - SlaterCoefficientCjcoupling(k, llarge, twoj, twojz1, llarge, twoj, twojz3) * SlaterCoefficientCjcoupling(k, lsmall, twoj, twojz4, lsmall, twoj, twojz2);
          if (Abs(val) > DBL_EPSILON)
          {
            oadd[3] = CrAnMask + Jindex[(twoj+twojz1)/2];
            oadd[2] = CrAnMask + Jindex[(twoj+twojz2)/2];
            oadd[1] =            Jindex[(twoj+twojz3)/2];
            oadd[0] =            Jindex[(twoj+twojz4)/2];
            RealOperatorAddLaderNormalOrder(U, oadd, 4, val * Fmix[k / 2], &QDetNotOrderedFcFaBcBa);
          }
        }
      }
    }
  }
  OperatorChop(U, 100 * DBL_EPSILON);
  return 1;
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
uint32_t CreateOperatorSx(uint32_t NFermion, uint32_t NSpin, uint16_t * SpinUp, uint16_t * SpinDn, OperatorType * const Sx)
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
uint32_t CreateOperatorSxY(uint32_t NFermion, uint32_t NSpin, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Sx)
{
  return CreateOperatorSx(NFermion, NSpin, LindexUp, LindexDn, Sx);
}
uint32_t CreateOperatorSxK(uint32_t NFermion, uint32_t NSpin, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Sx)
{
  return CreateOperatorSx(NFermion, NSpin, LindexUp, LindexDn, Sx);
}
uint32_t CreateOperatorSxZ(uint32_t NFermion, uint32_t NSpin, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Sx)
{
  return CreateOperatorSx(NFermion, NSpin, LindexUp, LindexDn, Sx);
}
uint32_t CreateOperatorSxP(uint32_t NFermion, uint32_t NSpin, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Sx)
{
  return CreateOperatorSx(NFermion, NSpin, LindexUp, LindexDn, Sx);
}
uint32_t CreateOperatorSxjjz(uint32_t NFermion, uint32_t L, uint16_t * Indexjmin, uint16_t * Indexjplus, OperatorType * const Sx)
{
  uint32_t i;
  uint16_t oadd[2];
  double m;
  const double l=(double)L;
  InitNamedOperator(Sx, "Sx");
  Sx->NFermion = NFermion;
  Sx->NWord = Min(Sx->NFermion, 1 + (Sx->NFermion - 1) / 8);
  InitOperatorFixedLength(Sx, 2, 8);
  for (i = 0; i < 2*L+1; i++)
  {
    m=(double)i-l;
    oadd[1] = CrAnMask + Indexjplus[i+1];
    oadd[0] = Indexjplus[i];
    RealOperatorAddLaderNormalOrder(Sx, oadd, 2, sqrt(l+m+1)*sqrt(l-m+1)/(4*l+2), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + Indexjplus[i];
    oadd[0] = Indexjplus[i+1];
    RealOperatorAddLaderNormalOrder(Sx, oadd, 2, sqrt(l+m+1)*sqrt(l-m+1)/(4*l+2), &QDetNotOrderedFcFaBcBa);
  }
  for (i = 1; i < 2*L+1; i++)
  {
    m=(double)i-l;
    oadd[1] = CrAnMask + Indexjplus[i+1];
    oadd[0] = Indexjmin[i-1];
    RealOperatorAddLaderNormalOrder(Sx, oadd, 2, sqrt(l+m+1)*sqrt(l+m)/(4*l+2), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + Indexjmin[i-1];
    oadd[0] = Indexjplus[i+1];
    RealOperatorAddLaderNormalOrder(Sx, oadd, 2, sqrt(l+m+1)*sqrt(l+m)/(4*l+2), &QDetNotOrderedFcFaBcBa);
  }
  for (i = 0; i < 2*L; i++)
  {
    m=(double)i-l;
    oadd[1] = CrAnMask + Indexjmin[i];
    oadd[0] = Indexjplus[i];
    RealOperatorAddLaderNormalOrder(Sx, oadd, 2, -sqrt(l-m+1)*sqrt(l-m)/(4*l+2), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + Indexjplus[i];
    oadd[0] = Indexjmin[i];
    RealOperatorAddLaderNormalOrder(Sx, oadd, 2, -sqrt(l-m+1)*sqrt(l-m)/(4*l+2), &QDetNotOrderedFcFaBcBa);
  }
  for (i = 1; i < 2*L; i++)
  {
    m=(double)i-l;
    oadd[1] = CrAnMask + Indexjmin[i];
    oadd[0] = Indexjmin[i-1];
    RealOperatorAddLaderNormalOrder(Sx, oadd, 2, -sqrt(l+m)*sqrt(l-m)/(4*l+2), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + Indexjmin[i-1];
    oadd[0] = Indexjmin[i];
    RealOperatorAddLaderNormalOrder(Sx, oadd, 2, -sqrt(l+m)*sqrt(l-m)/(4*l+2), &QDetNotOrderedFcFaBcBa);
  }
  OperatorChop(Sx, 0.001);
  return 0;
}

uint32_t CreateOperatorSy(uint32_t NFermion, uint32_t NSpin, uint16_t * SpinUp, uint16_t * SpinDn, OperatorType * const Sy)
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
    ComplexOperatorAddLaderNormalOrder(Sy, oadd, 2, 0.0, 0.5, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + SpinUp[i];
    oadd[0] = SpinDn[i];
    ComplexOperatorAddLaderNormalOrder(Sy, oadd, 2, 0.0, -0.5, &QDetNotOrderedFcFaBcBa);
  }
  OperatorChop(Sy, 0.001);
  return 0;
}
uint32_t CreateOperatorSyY(uint32_t NFermion, uint32_t NSpin, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Sy)
{
  return CreateOperatorSy(NFermion, NSpin, LindexUp, LindexDn, Sy);
}
uint32_t CreateOperatorSyK(uint32_t NFermion, uint32_t NSpin, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Sy)
{
  return CreateOperatorSy(NFermion, NSpin, LindexUp, LindexDn, Sy);
}
uint32_t CreateOperatorSyZ(uint32_t NFermion, uint32_t NSpin, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Sy)
{
  return CreateOperatorSy(NFermion, NSpin, LindexUp, LindexDn, Sy);
}
uint32_t CreateOperatorSyP(uint32_t NFermion, uint32_t NSpin, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Sy)
{
  return CreateOperatorSy(NFermion, NSpin, LindexUp, LindexDn, Sy);
}
uint32_t CreateOperatorSyjjz(uint32_t NFermion, uint32_t L, uint16_t * Indexjmin, uint16_t * Indexjplus, OperatorType * const Sy)
{
  uint32_t i;
  uint16_t oadd[2];
  double m;
  const double l=(double)L;
  InitNamedOperator(Sy, "Sy");
  Sy->NFermion = NFermion;
  Sy->NWord = Min(Sy->NFermion, 1 + (Sy->NFermion - 1) / 8);
  InitOperatorFixedLength(Sy, 2, 8);
  MakeRealOperatorFixedLengthComplex(Sy, 2);
  Sy->QComplex = 1;
  for (i = 0; i < 2*L+1; i++)
  {
    m=(double)i-l;
    oadd[1] = CrAnMask + Indexjplus[i+1];
    oadd[0] = Indexjplus[i];
    ComplexOperatorAddLaderNormalOrder(Sy, oadd, 2, 0,-sqrt(l+m+1)*sqrt(l-m+1)/(4*l+2), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + Indexjplus[i];
    oadd[0] = Indexjplus[i+1];
    ComplexOperatorAddLaderNormalOrder(Sy, oadd, 2, 0, sqrt(l+m+1)*sqrt(l-m+1)/(4*l+2), &QDetNotOrderedFcFaBcBa);
  }
  for (i = 1; i < 2*L+1; i++)
  {
    m=(double)i-l;
    oadd[1] = CrAnMask + Indexjplus[i+1];
    oadd[0] = Indexjmin[i-1];
    ComplexOperatorAddLaderNormalOrder(Sy, oadd, 2, 0,-sqrt(l+m+1)*sqrt(l+m)/(4*l+2), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + Indexjmin[i-1];
    oadd[0] = Indexjplus[i+1];
    ComplexOperatorAddLaderNormalOrder(Sy, oadd, 2, 0, sqrt(l+m+1)*sqrt(l+m)/(4*l+2), &QDetNotOrderedFcFaBcBa);
  }
  for (i = 0; i < 2*L; i++)
  {
    m=(double)i-l;
    oadd[1] = CrAnMask + Indexjmin[i];
    oadd[0] = Indexjplus[i];
    ComplexOperatorAddLaderNormalOrder(Sy, oadd, 2, 0,  sqrt(l-m+1)*sqrt(l-m)/(4*l+2), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + Indexjplus[i];
    oadd[0] = Indexjmin[i];
    ComplexOperatorAddLaderNormalOrder(Sy, oadd, 2, 0,-sqrt(l-m+1)*sqrt(l-m)/(4*l+2), &QDetNotOrderedFcFaBcBa);
  }
  for (i = 1; i < 2*L; i++)
  {
    m=(double)i-l;
    oadd[1] = CrAnMask + Indexjmin[i];
    oadd[0] = Indexjmin[i-1];
    ComplexOperatorAddLaderNormalOrder(Sy, oadd, 2, 0, sqrt(l+m)*sqrt(l-m)/(4*l+2), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + Indexjmin[i-1];
    oadd[0] = Indexjmin[i];
    ComplexOperatorAddLaderNormalOrder(Sy, oadd, 2, 0,-sqrt(l+m)*sqrt(l-m)/(4*l+2), &QDetNotOrderedFcFaBcBa);
  }
  OperatorChop(Sy, 0.001);
  return 0;
}

uint32_t CreateOperatorSz(uint32_t NFermion, uint32_t NSpin, uint16_t * SpinUp, uint16_t * SpinDn, OperatorType * const Sz)
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
uint32_t CreateOperatorSzY(uint32_t NFermion, uint32_t NSpin, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Sz)
{
  return CreateOperatorSz(NFermion, NSpin, LindexUp, LindexDn, Sz);
}
uint32_t CreateOperatorSzK(uint32_t NFermion, uint32_t NSpin, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Sz)
{
  return CreateOperatorSz(NFermion, NSpin, LindexUp, LindexDn, Sz);
}
uint32_t CreateOperatorSzZ(uint32_t NFermion, uint32_t NSpin, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Sz)
{
  return CreateOperatorSz(NFermion, NSpin, LindexUp, LindexDn, Sz);
}
uint32_t CreateOperatorSzP(uint32_t NFermion, uint32_t NSpin, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Sz)
{
  return CreateOperatorSz(NFermion, NSpin, LindexUp, LindexDn, Sz);
}
uint32_t CreateOperatorSzjjz(uint32_t NFermion, uint32_t L, uint16_t * Jindexjmin, uint16_t * Jindexjplus, OperatorType * const Sz)
{
  uint32_t i;
  uint16_t oadd[2];
  const double l=(double)L;
  double m,f;
  InitNamedOperator(Sz, "Sz");
  Sz->NFermion = NFermion;
  Sz->NWord = Min(Sz->NFermion, 1 + (Sz->NFermion - 1) / 8);
  InitOperatorFixedLength(Sz, 2, 8);
  for (i = 0; i < 2*L+2; i++)
  {
    m=(double)i-l-1;
    //jz = m+1/2 running from -l-1/2 (i=0) to l+1/2 (i=2l+1)
    oadd[1] = CrAnMask + Jindexjplus[i];
    oadd[0] = Jindexjplus[i];
    RealOperatorAddLaderNormalOrder(Sz, oadd, 2,(m+0.5)/(2*l+1), &QDetNotOrderedFcFaBcBa);
  }
  for (i = 0; i < 2*L; i++)
  {
    m=(double)i-l;
    //jz = m+1/2 running from -l+1/2 (i=0) to l-1/2 (i=2l-1)
    oadd[1] = CrAnMask + Jindexjmin[i];
    oadd[0] = Jindexjmin[i];
    RealOperatorAddLaderNormalOrder(Sz, oadd, 2,-(m+0.5)/(2*l+1), &QDetNotOrderedFcFaBcBa);
    f=(-sqrt(l-m)*sqrt(l+m+1))/(2*l+1);
    oadd[1] = CrAnMask + Jindexjplus[i+1];
    oadd[0] = Jindexjmin[i];
    RealOperatorAddLaderNormalOrder(Sz, oadd, 2,f, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + Jindexjmin[i];
    oadd[0] = Jindexjplus[i+1];
    RealOperatorAddLaderNormalOrder(Sz, oadd, 2,f, &QDetNotOrderedFcFaBcBa);
  }
  OperatorChop(Sz, 0.001);
  return 0;
}

uint32_t CreateOperatorSsqr(uint32_t NFermion, uint32_t NSpin, uint16_t * SpinUp, uint16_t * SpinDn, OperatorType * const Ssqr)
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
uint32_t CreateOperatorSsqrY(uint32_t NFermion, uint32_t NSpin, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Ssqr)
{
  return CreateOperatorSsqr(NFermion,NSpin, LindexUp, LindexDn, Ssqr);
}
uint32_t CreateOperatorSsqrK(uint32_t NFermion, uint32_t NSpin, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Ssqr)
{
  return CreateOperatorSsqr(NFermion,NSpin, LindexUp, LindexDn, Ssqr);
}
uint32_t CreateOperatorSsqrZ(uint32_t NFermion, uint32_t NSpin, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Ssqr)
{
  return CreateOperatorSsqr(NFermion,NSpin, LindexUp, LindexDn, Ssqr);
}
uint32_t CreateOperatorSsqrP(uint32_t NFermion, uint32_t NSpin, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Ssqr)
{
  return CreateOperatorSsqr(NFermion,NSpin, LindexUp, LindexDn, Ssqr);
}
uint32_t CreateOperatorSsqrjjz(uint32_t NFermion, uint32_t L, uint16_t * Indexjmin, uint16_t * Indexjplus, OperatorType * const Ssqr)
{
  //TODO: Implement S^2 directly
  OperatorType Sx,Sy,Sz;
  InitNamedOperator(Ssqr, "S^2");
  Ssqr->NFermion = NFermion;
  Ssqr->NWord = Min(Ssqr->NFermion, 1 + (Ssqr->NFermion - 1) / 8);
  InitOperatorFixedLength(Ssqr, 2, 8);
  InitOperatorFixedLength(Ssqr, 4, 8);
  CreateOperatorSxjjz(NFermion, L, Indexjmin, Indexjplus, &Sx);
  CreateOperatorSyjjz(NFermion, L, Indexjmin, Indexjplus, &Sy);
  CreateOperatorSzjjz(NFermion, L, Indexjmin, Indexjplus, &Sz);
  OperatorMultiplyNormalOrder(Ssqr, Sx, Sx, &QDetNotOrderedFcFaBcBa);
  OperatorMultiplyNormalOrder(Ssqr, Sy, Sy, &QDetNotOrderedFcFaBcBa);
  OperatorMultiplyNormalOrder(Ssqr, Sz, Sz, &QDetNotOrderedFcFaBcBa);
  FreeOperator(&Sx);
  FreeOperator(&Sy);
  FreeOperator(&Sz);
  OperatorChop(Ssqr, 0.001);
  return 0;
}

uint32_t CreateOperatorSplus(uint32_t NFermion, uint32_t NSpin, uint16_t * SpinUp, uint16_t * SpinDn, OperatorType * const Splus)
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
uint32_t CreateOperatorSplusY(uint32_t NFermion, uint32_t NSpin, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Splus)
{
  return CreateOperatorSplus(NFermion, NSpin, LindexUp, LindexDn, Splus);
}
uint32_t CreateOperatorSplusK(uint32_t NFermion, uint32_t NSpin, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Splus)
{
  return CreateOperatorSplus(NFermion, NSpin, LindexUp, LindexDn, Splus);
}
uint32_t CreateOperatorSplusZ(uint32_t NFermion, uint32_t NSpin, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Splus)
{
  return CreateOperatorSplus(NFermion, NSpin, LindexUp, LindexDn, Splus);
}
uint32_t CreateOperatorSplusP(uint32_t NFermion, uint32_t NSpin, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Splus)
{
  return CreateOperatorSplus(NFermion, NSpin, LindexUp, LindexDn, Splus);
}
uint32_t CreateOperatorSplusjjz(uint32_t NFermion, uint32_t L, uint16_t * Indexjmin, uint16_t * Indexjplus, OperatorType * const Splus)
{
  uint32_t i;
  uint16_t oadd[2];
  double m;
  const double l=(double)L;
  InitNamedOperator(Splus, "S^+");
  Splus->NFermion = NFermion;
  Splus->NWord = Min(Splus->NFermion, 1 + (Splus->NFermion - 1) / 8);
  InitOperatorFixedLength(Splus, 2, 8);
  for (i = 0; i < 2*L+1; i++)
  {
    m=(double)i-l;
    oadd[1] = CrAnMask + Indexjplus[i+1];
    oadd[0] = Indexjplus[i];
    RealOperatorAddLaderNormalOrder(Splus, oadd, 2, sqrt(l+m+1)*sqrt(l-m+1)/(2*l+1), &QDetNotOrderedFcFaBcBa);
  }
  for (i = 1; i < 2*L+1; i++)
  {
    m=(double)i-l;
    oadd[1] = CrAnMask + Indexjplus[i+1];
    oadd[0] = Indexjmin[i-1];
    RealOperatorAddLaderNormalOrder(Splus, oadd, 2, sqrt(l+m+1)*sqrt(l+m)/(2*l+1), &QDetNotOrderedFcFaBcBa);
  }
  for (i = 0; i < 2*L; i++)
  {
    m=(double)i-l;
    oadd[1] = CrAnMask + Indexjmin[i];
    oadd[0] = Indexjplus[i];
    RealOperatorAddLaderNormalOrder(Splus, oadd, 2, -sqrt(l-m+1)*sqrt(l-m)/(2*l+1), &QDetNotOrderedFcFaBcBa);
  }
  for (i = 1; i < 2*L; i++)
  {
    m=(double)i-l;
    oadd[1] = CrAnMask + Indexjmin[i];
    oadd[0] = Indexjmin[i-1];
    RealOperatorAddLaderNormalOrder(Splus, oadd, 2, -sqrt(l+m)*sqrt(l-m)/(2*l+1), &QDetNotOrderedFcFaBcBa);
  }
  OperatorChop(Splus, 0.001);
  return 0;
}


uint32_t CreateOperatorSmin(uint32_t NFermion, uint32_t NSpin, uint16_t * SpinUp, uint16_t * SpinDn, OperatorType * const Smin)
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
uint32_t CreateOperatorSminY(uint32_t NFermion, uint32_t NSpin, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Smin)
{
  return CreateOperatorSmin(NFermion, NSpin, LindexUp, LindexDn, Smin);
}
uint32_t CreateOperatorSminK(uint32_t NFermion, uint32_t NSpin, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Smin)
{
  return CreateOperatorSmin(NFermion, NSpin, LindexUp, LindexDn, Smin);
}
uint32_t CreateOperatorSminZ(uint32_t NFermion, uint32_t NSpin, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Smin)
{
  return CreateOperatorSmin(NFermion, NSpin, LindexUp, LindexDn, Smin);
}
uint32_t CreateOperatorSminP(uint32_t NFermion, uint32_t NSpin, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Smin)
{
  return CreateOperatorSmin(NFermion, NSpin, LindexUp, LindexDn, Smin);
}
uint32_t CreateOperatorSminjjz(uint32_t NFermion, uint32_t L, uint16_t * Indexjmin, uint16_t * Indexjplus, OperatorType * const Smin)
{
  uint32_t i;
  uint16_t oadd[2];
  double m;
  const double l=(double)L;
  InitNamedOperator(Smin, "S^-");
  Smin->NFermion = NFermion;
  Smin->NWord = Min(Smin->NFermion, 1 + (Smin->NFermion - 1) / 8);
  InitOperatorFixedLength(Smin, 2, 8);
  for (i = 0; i < 2*L+1; i++)
  {
    m=(double)i-l;
    oadd[1] = CrAnMask + Indexjplus[i];
    oadd[0] = Indexjplus[i+1];
    RealOperatorAddLaderNormalOrder(Smin, oadd, 2, sqrt(l+m+1)*sqrt(l-m+1)/(2*l+1), &QDetNotOrderedFcFaBcBa);
  }
  for (i = 1; i < 2*L+1; i++)
  {
    m=(double)i-l;
    oadd[1] = CrAnMask + Indexjmin[i-1];
    oadd[0] = Indexjplus[i+1];
    RealOperatorAddLaderNormalOrder(Smin, oadd, 2, sqrt(l+m+1)*sqrt(l+m)/(2*l+1), &QDetNotOrderedFcFaBcBa);
  }
  for (i = 0; i < 2*L; i++)
  {
    m=(double)i-l;
    oadd[1] = CrAnMask + Indexjplus[i];
    oadd[0] = Indexjmin[i];
    RealOperatorAddLaderNormalOrder(Smin, oadd, 2, -sqrt(l-m+1)*sqrt(l-m)/(2*l+1), &QDetNotOrderedFcFaBcBa);
  }
  for (i = 1; i < 2*L; i++)
  {
    m=(double)i-l;
    oadd[1] = CrAnMask + Indexjmin[i-1];
    oadd[0] = Indexjmin[i];
    RealOperatorAddLaderNormalOrder(Smin, oadd, 2, -sqrt(l+m)*sqrt(l-m)/(2*l+1), &QDetNotOrderedFcFaBcBa);
  }
  OperatorChop(Smin, 0.001);
  return 0;
}

uint32_t CreateSpinOperators(uint32_t NFermion, uint32_t NSpin, uint16_t * SpinUp, uint16_t * SpinDn, OperatorType * const Sx, OperatorType * const Sy, OperatorType * const Sz, OperatorType * const S2)
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
    ComplexOperatorAddLaderNormalOrder(Ly, oadd, 2, 0.0, -0.5 * sqrt(l + m + 1) * sqrt(l - m), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + Lindex[i];
    oadd[0] = Lindex[i + 1];
    ComplexOperatorAddLaderNormalOrder(Ly, oadd, 2, 0.0, 0.5 * sqrt(l + m + 1) * sqrt(l - m), &QDetNotOrderedFcFaBcBa);
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
uint32_t CreateOperatorLx(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Lx)
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
uint32_t CreateOperatorLxY(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Lx)
{
  return CreateOperatorLx(NFermion, L, LindexUp, LindexDn, Lx);
}
uint32_t CreateOperatorLxK(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Lx)
{
  return CreateOperatorInShellBasis(NFermion, L, 2, LindexUp, LindexDn, &CreateOperatorLx, Lx);
}
uint32_t CreateOperatorLxZ(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Lx)
{
  return CreateOperatorInShellBasis(NFermion, L, 1, LindexUp, LindexDn, &CreateOperatorLx, Lx);
}
uint32_t CreateOperatorLxP(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Lx)
{
  return CreateOperatorInShellBasis(NFermion, L, 3, LindexUp, LindexDn, &CreateOperatorLx, Lx);
}
uint32_t CreateOperatorLxjjz(uint32_t NFermion, uint32_t L, uint16_t * Indexjmin, uint16_t * Indexjplus, OperatorType * const Lx)
{
  uint32_t i;
  uint16_t oadd[2];
  double j,jz,m,f;
  const double l=(double)L;
  InitNamedOperator(Lx, "Lx");
  Lx->NFermion = NFermion;
  Lx->NWord = Min(Lx->NFermion, 1 + (Lx->NFermion - 1) / 8);
  InitOperatorFixedLength(Lx, 2, 8);
  // add terms for j=L+1/2
  j=(double)L+0.5;
  for (i = 0; i < 2*L+1; i++)
  {
    m=(double)i-l;
    jz=(double)i - j;
    f=0.5*sqrt(j+jz+1)*sqrt(j-jz)-sqrt(l+m+1)*sqrt(l-m+1)/(4*l+2);
    oadd[1] = CrAnMask + Indexjplus[i+1];
    oadd[0] = Indexjplus[i];
    RealOperatorAddLaderNormalOrder(Lx, oadd, 2,f, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + Indexjplus[i];
    oadd[0] = Indexjplus[i+1];
    RealOperatorAddLaderNormalOrder(Lx, oadd, 2,f, &QDetNotOrderedFcFaBcBa);
  }
  for (i = 1; i < 2*L+1; i++)
  {
    m=(double)i-l;
    oadd[1] = CrAnMask + Indexjplus[i+1];
    oadd[0] = Indexjmin[i-1];
    RealOperatorAddLaderNormalOrder(Lx, oadd, 2,-sqrt(l+m+1)*sqrt(l+m)/(4*l+2), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + Indexjmin[i-1];
    oadd[0] = Indexjplus[i+1];
    RealOperatorAddLaderNormalOrder(Lx, oadd, 2,-sqrt(l+m+1)*sqrt(l+m)/(4*l+2), &QDetNotOrderedFcFaBcBa);
  }
  for (i = 0; i < 2*L; i++)
  {
    m=(double)i-l;
    oadd[1] = CrAnMask + Indexjmin[i];
    oadd[0] = Indexjplus[i];
    RealOperatorAddLaderNormalOrder(Lx, oadd, 2, sqrt(l-m+1)*sqrt(l-m)/(4*l+2), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + Indexjplus[i];
    oadd[0] = Indexjmin[i];
    RealOperatorAddLaderNormalOrder(Lx, oadd, 2, sqrt(l-m+1)*sqrt(l-m)/(4*l+2), &QDetNotOrderedFcFaBcBa);
  }
  // add terms for j=L-1/2
  j=(double)L-0.5;
  for (i = 1; i < 2 * L ; i++)
  {
    jz=(double)i - j - 1;
    m=(double)i-l;
    f=0.5*sqrt(j+jz+1)*sqrt(j-jz)+sqrt(l+m)*sqrt(l-m)/(4*l+2);
    oadd[1] = CrAnMask + Indexjmin[i];
    oadd[0] = Indexjmin[i-1];
    RealOperatorAddLaderNormalOrder(Lx, oadd, 2, f, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + Indexjmin[i-1];
    oadd[0] = Indexjmin[i];
    RealOperatorAddLaderNormalOrder(Lx, oadd, 2, f, &QDetNotOrderedFcFaBcBa);
  }
  OperatorChop(Lx, 0.001);
  return 0;
}



uint32_t CreateOperatorLy(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Ly)
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
    ComplexOperatorAddLaderNormalOrder(Ly, oadd, 2, 0.0,-0.5 * sqrt(l + m + 1) * sqrt(l - m), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + LindexUp[i];
    oadd[0] = LindexUp[i + 1];
    ComplexOperatorAddLaderNormalOrder(Ly, oadd, 2, 0.0, 0.5 * sqrt(l + m + 1) * sqrt(l - m), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + LindexDn[i + 1];
    oadd[0] = LindexDn[i];
    ComplexOperatorAddLaderNormalOrder(Ly, oadd, 2, 0.0,-0.5 * sqrt(l + m + 1) * sqrt(l - m), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + LindexDn[i];
    oadd[0] = LindexDn[i + 1];
    ComplexOperatorAddLaderNormalOrder(Ly, oadd, 2, 0.0, 0.5 * sqrt(l + m + 1) * sqrt(l - m), &QDetNotOrderedFcFaBcBa);
  }
  OperatorChop(Ly, 0.001);
  return 0;
}
uint32_t CreateOperatorLyY(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Ly)
{
  return CreateOperatorLy(NFermion, L, LindexUp, LindexDn, Ly);
}
uint32_t CreateOperatorLyK(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Ly)
{
  return CreateOperatorInShellBasis(NFermion, L, 2, LindexUp, LindexDn, &CreateOperatorLy, Ly);
}
uint32_t CreateOperatorLyZ(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Ly)
{
  return CreateOperatorInShellBasis(NFermion, L, 1, LindexUp, LindexDn, &CreateOperatorLy, Ly);
}
uint32_t CreateOperatorLyP(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Ly)
{
  return CreateOperatorInShellBasis(NFermion, L, 3, LindexUp, LindexDn, &CreateOperatorLy, Ly);
}
uint32_t CreateOperatorLyjjz(uint32_t NFermion, uint32_t L, uint16_t * Indexjmin, uint16_t * Indexjplus, OperatorType * const Ly)
{
  uint32_t i;
  uint16_t oadd[2];
  double j,jz,m,f;
  const double l=(double)L;
  InitNamedOperator(Ly, "Ly");
  Ly->NFermion = NFermion;
  Ly->NWord = Min(Ly->NFermion, 1 + (Ly->NFermion - 1) / 8);
  InitOperatorFixedLength(Ly, 2, 8);
  MakeRealOperatorFixedLengthComplex(Ly, 2);
  Ly->QComplex = 1;
  // add terms for j=L+1/2
  j=(double)L+0.5;
  for (i = 0; i < 2 * L + 1; i++)
  {
    jz=(double)i - j;
    m=(double)i-l;
    f=0.5*sqrt(j+jz+1)*sqrt(j-jz)-sqrt(l+m+1)*sqrt(l-m+1)/(4*l+2);
    oadd[1] = CrAnMask + Indexjplus[i+1];
    oadd[0] = Indexjplus[i];
    ComplexOperatorAddLaderNormalOrder(Ly, oadd, 2, 0.0,-f, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + Indexjplus[i];
    oadd[0] = Indexjplus[i+1];
    ComplexOperatorAddLaderNormalOrder(Ly, oadd, 2, 0.0, f, &QDetNotOrderedFcFaBcBa);
  }
  for (i = 1; i < 2*L+1; i++)
  {
    m=(double)i-l;
    oadd[1] = CrAnMask + Indexjplus[i+1];
    oadd[0] = Indexjmin[i-1];
    ComplexOperatorAddLaderNormalOrder(Ly, oadd, 2, 0, sqrt(l+m+1)*sqrt(l+m)/(4*l+2), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + Indexjmin[i-1];
    oadd[0] = Indexjplus[i+1];
    ComplexOperatorAddLaderNormalOrder(Ly, oadd, 2, 0,-sqrt(l+m+1)*sqrt(l+m)/(4*l+2), &QDetNotOrderedFcFaBcBa);
  }
  for (i = 0; i < 2*L; i++)
  {
    m=(double)i-l;
    oadd[1] = CrAnMask + Indexjmin[i];
    oadd[0] = Indexjplus[i];
    ComplexOperatorAddLaderNormalOrder(Ly, oadd, 2, 0,-sqrt(l-m+1)*sqrt(l-m)/(4*l+2), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + Indexjplus[i];
    oadd[0] = Indexjmin[i];
    ComplexOperatorAddLaderNormalOrder(Ly, oadd, 2, 0, sqrt(l-m+1)*sqrt(l-m)/(4*l+2), &QDetNotOrderedFcFaBcBa);
  }
  // add terms for j=L-1/2
  j=(double)L-0.5;
  for (i = 1; i < 2 * L ; i++)
  {
    jz=(double)i - j - 1;
    m=(double)i-l;
    f=0.5*sqrt(j+jz+1)*sqrt(j-jz)+sqrt(l+m)*sqrt(l-m)/(4*l+2);
    oadd[1] = CrAnMask + Indexjmin[i];
    oadd[0] = Indexjmin[i-1];
    ComplexOperatorAddLaderNormalOrder(Ly, oadd, 2, 0.0,-f, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + Indexjmin[i-1];
    oadd[0] = Indexjmin[i];
    ComplexOperatorAddLaderNormalOrder(Ly, oadd, 2, 0.0, f, &QDetNotOrderedFcFaBcBa);
  }
  OperatorChop(Ly, 0.001);
  return 0;
}

uint32_t CreateOperatorLz(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Lz)
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
uint32_t CreateOperatorLzY(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Lz)
{
  return CreateOperatorLz(NFermion, L, LindexUp, LindexDn, Lz);
}
uint32_t CreateOperatorLzK(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Lz)
{
  return CreateOperatorInShellBasis(NFermion, L, 2, LindexUp, LindexDn, &CreateOperatorLz, Lz);
}
uint32_t CreateOperatorLzZ(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Lz)
{
  return CreateOperatorInShellBasis(NFermion, L, 1, LindexUp, LindexDn, &CreateOperatorLz, Lz);
}
uint32_t CreateOperatorLzP(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Lz)
{
  return CreateOperatorInShellBasis(NFermion, L, 3, LindexUp, LindexDn, &CreateOperatorLz, Lz);
}
uint32_t CreateOperatorLzjjz(uint32_t NFermion, uint32_t L, uint16_t * Jindexjmin, uint16_t * Jindexjplus, OperatorType * const Lz)
{
  uint32_t i;
  uint16_t oadd[2];
  const double l=(double)L;
  double m,f;
  InitNamedOperator(Lz, "Lz");
  Lz->NFermion = NFermion;
  Lz->NWord = Min(Lz->NFermion, 1 + (Lz->NFermion - 1) / 8);
  InitOperatorFixedLength(Lz, 2, 8);
  for (i = 0; i < 2*L+2; i++)
  {
    m=(double)i-l-1;
    //jz = m+1/2 running from -l-1/2 (i=0) to l+1/2 (i=2l+1)
    oadd[1] = CrAnMask + Jindexjplus[i];
    oadd[0] = Jindexjplus[i];
    RealOperatorAddLaderNormalOrder(Lz, oadd, 2,(2*m+1)*l/(2*l+1), &QDetNotOrderedFcFaBcBa);
  }
  for (i = 0; i < 2*L; i++)
  {
    m=(double)i-l;
    //jz = m+1/2 running from -l+1/2 (i=0) to l-1/2 (i=2l-1)
    oadd[1] = CrAnMask + Jindexjmin[i];
    oadd[0] = Jindexjmin[i];
    RealOperatorAddLaderNormalOrder(Lz, oadd, 2,(2*m+1)*(l+1)/(2*l+1), &QDetNotOrderedFcFaBcBa);
    f=(sqrt(l-m)*sqrt(l+m+1))/(2*l+1);
    oadd[1] = CrAnMask + Jindexjplus[i+1];
    oadd[0] = Jindexjmin[i];
    RealOperatorAddLaderNormalOrder(Lz, oadd, 2,f, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + Jindexjmin[i];
    oadd[0] = Jindexjplus[i+1];
    RealOperatorAddLaderNormalOrder(Lz, oadd, 2,f, &QDetNotOrderedFcFaBcBa);
  }
  OperatorChop(Lz, 0.001);
  return 0;
}

uint32_t CreateOperatorLsqr(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Lsqr)
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
      RealOperatorAddLaderNormalOrder(Lsqr, oadd, 4, -sqrt(l + mi + 1) * sqrt(l - mi) * sqrt(l - mj) * sqrt(l + mj + 1), &QDetNotOrderedFcFaBcBa);
      oadd[3] = CrAnMask + LindexDn[i + 1];
      oadd[2] = CrAnMask + LindexDn[j];
      oadd[1] = LindexDn[i];
      oadd[0] = LindexDn[j + 1];
      RealOperatorAddLaderNormalOrder(Lsqr, oadd, 4, -sqrt(l + mi + 1) * sqrt(l - mi) * sqrt(l - mj) * sqrt(l + mj + 1),  &QDetNotOrderedFcFaBcBa);
      oadd[3] = CrAnMask + LindexUp[i + 1];
      oadd[2] = CrAnMask + LindexDn[j];
      oadd[1] = LindexUp[i];
      oadd[0] = LindexDn[j + 1];
      RealOperatorAddLaderNormalOrder(Lsqr, oadd, 4, -sqrt(l + mi + 1) * sqrt(l - mi) * sqrt(l - mj) * sqrt(l + mj + 1), &QDetNotOrderedFcFaBcBa);
      oadd[3] = CrAnMask + LindexDn[i + 1];
      oadd[2] = CrAnMask + LindexUp[j];
      oadd[1] = LindexDn[i];
      oadd[0] = LindexUp[j + 1];
      RealOperatorAddLaderNormalOrder(Lsqr, oadd, 4, -sqrt(l + mi + 1) * sqrt(l - mi) * sqrt(l - mj) * sqrt(l + mj + 1), &QDetNotOrderedFcFaBcBa);
    }
  }
  OperatorChop(Lsqr, 0.001);
  return 0;
}
uint32_t CreateOperatorLsqrY(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Lsqr)
{
  return CreateOperatorLsqr(NFermion, L, LindexUp, LindexDn, Lsqr);
}
uint32_t CreateOperatorLsqrK(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Lsqr)
{
  return CreateOperatorInShellBasis(NFermion, L, 2, LindexUp, LindexDn, &CreateOperatorLsqr, Lsqr);
}
uint32_t CreateOperatorLsqrZ(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Lsqr)
{
  return CreateOperatorInShellBasis(NFermion, L, 1, LindexUp, LindexDn, &CreateOperatorLsqr, Lsqr);
}
uint32_t CreateOperatorLsqrP(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Lsqr)
{
  return CreateOperatorInShellBasis(NFermion, L, 3, LindexUp, LindexDn, &CreateOperatorLsqr, Lsqr);
}
uint32_t CreateOperatorLsqrjjz(uint32_t NFermion, uint32_t L, uint16_t * Indexjmin, uint16_t * Indexjplus, OperatorType * const Lsqr)
{
  //TODO: Implement L^2 directly
  OperatorType Lx,Ly,Lz;
  InitNamedOperator(Lsqr, "L^2");
  Lsqr->NFermion = NFermion;
  Lsqr->NWord = Min(Lsqr->NFermion, 1 + (Lsqr->NFermion - 1) / 8);
  InitOperatorFixedLength(Lsqr, 2, 8);
  InitOperatorFixedLength(Lsqr, 4, 8);
  CreateOperatorLxjjz(NFermion, L, Indexjmin, Indexjplus, &Lx);
  CreateOperatorLyjjz(NFermion, L, Indexjmin, Indexjplus, &Ly);
  CreateOperatorLzjjz(NFermion, L, Indexjmin, Indexjplus, &Lz);
  OperatorMultiplyNormalOrder(Lsqr, Lx, Lx, &QDetNotOrderedFcFaBcBa);
  OperatorMultiplyNormalOrder(Lsqr, Ly, Ly, &QDetNotOrderedFcFaBcBa);
  OperatorMultiplyNormalOrder(Lsqr, Lz, Lz, &QDetNotOrderedFcFaBcBa);
  FreeOperator(&Lx);
  FreeOperator(&Ly);
  FreeOperator(&Lz);
  OperatorChop(Lsqr, 0.001);
  return 0;
}

uint32_t CreateOperatorLplus(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Lplus)
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
uint32_t CreateOperatorLplusY(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Lplus)
{
  return CreateOperatorLplus(NFermion, L, LindexUp, LindexDn, Lplus);
}
uint32_t CreateOperatorLplusK(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Lplus)
{
  return CreateOperatorInShellBasis(NFermion, L, 2, LindexUp, LindexDn, &CreateOperatorLplus, Lplus);
}
uint32_t CreateOperatorLplusZ(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Lplus)
{
  return CreateOperatorInShellBasis(NFermion, L, 1, LindexUp, LindexDn, &CreateOperatorLplus, Lplus);
}
uint32_t CreateOperatorLplusP(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Lplus)
{
  return CreateOperatorInShellBasis(NFermion, L, 3, LindexUp, LindexDn, &CreateOperatorLplus, Lplus);
}
uint32_t CreateOperatorLplusjjz(uint32_t NFermion, uint32_t L, uint16_t * Indexjmin, uint16_t * Indexjplus, OperatorType * const Lplus)
{
  uint32_t i;
  uint16_t oadd[2];
  double j,jz,m;
  const double l=(double)L;
  InitNamedOperator(Lplus, "L^+");
  Lplus->NFermion = NFermion;
  Lplus->NWord = Min(Lplus->NFermion, 1 + (Lplus->NFermion - 1) / 8);
  InitOperatorFixedLength(Lplus, 2, 8);
  // add terms for j=L+1/2
  j=(double)L+0.5;
  for (i = 0; i < 2 * L + 1; i++)
  {
    m=(double)i-l;
    jz=(double)i - j;
    oadd[1] = CrAnMask + Indexjplus[i+1];
    oadd[0] = Indexjplus[i];
    RealOperatorAddLaderNormalOrder(Lplus, oadd, 2, sqrt(j+jz+1)*sqrt(j-jz)-sqrt(l+m+1)*sqrt(l-m+1)/(2*l+1), &QDetNotOrderedFcFaBcBa);
  }
  for (i = 1; i < 2*L+1; i++)
  {
    m=(double)i-l;
    oadd[1] = CrAnMask + Indexjplus[i+1];
    oadd[0] = Indexjmin[i-1];
    RealOperatorAddLaderNormalOrder(Lplus, oadd, 2,-sqrt(l+m+1)*sqrt(l+m)/(2*l+1), &QDetNotOrderedFcFaBcBa);
  }
  for (i = 0; i < 2*L; i++)
  {
    m=(double)i-l;
    oadd[1] = CrAnMask + Indexjmin[i];
    oadd[0] = Indexjplus[i];
    RealOperatorAddLaderNormalOrder(Lplus, oadd, 2, sqrt(l-m+1)*sqrt(l-m)/(2*l+1), &QDetNotOrderedFcFaBcBa);
  }
  // add terms for j=L-1/2
  j=(double)L-0.5;
  for (i = 1; i < 2 * L ; i++)
  {
    m=(double)i-l;
    jz=(double)i - j -1;
    oadd[1] = CrAnMask + Indexjmin[i];
    oadd[0] = Indexjmin[i-1];
    RealOperatorAddLaderNormalOrder(Lplus, oadd, 2, sqrt(j+jz+1)*sqrt(j-jz)+sqrt(l+m)*sqrt(l-m)/(2*l+1), &QDetNotOrderedFcFaBcBa);
  }
  OperatorChop(Lplus, 0.001);
  return 0;
}

uint32_t CreateOperatorLmin(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Lmin)
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
uint32_t CreateOperatorLminY(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Lmin)
{
  return CreateOperatorLmin(NFermion, L, LindexUp, LindexDn, Lmin);
}
uint32_t CreateOperatorLminK(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Lmin)
{
  return CreateOperatorInShellBasis(NFermion, L, 2, LindexUp, LindexDn, &CreateOperatorLmin, Lmin);
}
uint32_t CreateOperatorLminZ(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Lmin)
{
  return CreateOperatorInShellBasis(NFermion, L, 1, LindexUp, LindexDn, &CreateOperatorLmin, Lmin);
}
uint32_t CreateOperatorLminP(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Lmin)
{
  return CreateOperatorInShellBasis(NFermion, L, 3, LindexUp, LindexDn, &CreateOperatorLmin, Lmin);
}
uint32_t CreateOperatorLminjjz(uint32_t NFermion, uint32_t L, uint16_t * Indexjmin, uint16_t * Indexjplus, OperatorType * const Lmin)
{
  uint32_t i;
  uint16_t oadd[2];
  double j,jz,m;
  const double l=(double)L;
  InitNamedOperator(Lmin, "L^-");
  Lmin->NFermion = NFermion;
  Lmin->NWord = Min(Lmin->NFermion, 1 + (Lmin->NFermion - 1) / 8);
  InitOperatorFixedLength(Lmin, 2, 8);
  // add terms for j=L+1/2
  j=(double)L+0.5;
  for (i = 0; i < 2 * L + 1; i++)
  {
    jz=(double)i - j;
    m=(double)i-l;
    oadd[1] = CrAnMask + Indexjplus[i];
    oadd[0] = Indexjplus[i+1];
    RealOperatorAddLaderNormalOrder(Lmin, oadd, 2, sqrt(j+jz+1)*sqrt(j-jz)-sqrt(l+m+1)*sqrt(l-m+1)/(2*l+1), &QDetNotOrderedFcFaBcBa);
  }
  for (i = 1; i < 2*L+1; i++)
  {
    m=(double)i-l;
    oadd[1] = CrAnMask + Indexjmin[i-1];
    oadd[0] = Indexjplus[i+1];
    RealOperatorAddLaderNormalOrder(Lmin, oadd, 2,-sqrt(l+m+1)*sqrt(l+m)/(2*l+1), &QDetNotOrderedFcFaBcBa);
  }
  for (i = 0; i < 2*L; i++)
  {
    m=(double)i-l;
    oadd[1] = CrAnMask + Indexjplus[i];
    oadd[0] = Indexjmin[i];
    RealOperatorAddLaderNormalOrder(Lmin, oadd, 2, sqrt(l-m+1)*sqrt(l-m)/(2*l+1), &QDetNotOrderedFcFaBcBa);
  }
  // add terms for j=L-1/2
  j=(double)L-0.5;
  for (i = 1; i < 2 * L; i++)
  {
    jz=(double)i - j - 1;
    m=(double)i-l;
    oadd[1] = CrAnMask + Indexjmin[i-1];
    oadd[0] = Indexjmin[i];
    RealOperatorAddLaderNormalOrder(Lmin, oadd, 2, sqrt(j+jz+1)*sqrt(j-jz)+sqrt(l+m)*sqrt(l-m)/(2*l+1), &QDetNotOrderedFcFaBcBa);
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
uint32_t CreateOperatorJx(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Jx)
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
uint32_t CreateOperatorJxY(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Jx)
{
  return CreateOperatorJx(NFermion, L, LindexUp, LindexDn, Jx);
}
uint32_t CreateOperatorJxK(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Jx)
{
  return CreateOperatorInShellBasis(NFermion, L, 2, LindexUp, LindexDn, &CreateOperatorJx, Jx);
}
uint32_t CreateOperatorJxZ(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Jx)
{
  return CreateOperatorInShellBasis(NFermion, L, 1, LindexUp, LindexDn, &CreateOperatorJx, Jx);
}
uint32_t CreateOperatorJxP(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Jx)
{
  return CreateOperatorInShellBasis(NFermion, L, 3, LindexUp, LindexDn, &CreateOperatorJx, Jx);
}
uint32_t CreateOperatorJxjjz(uint32_t NFermion, uint32_t L, uint16_t * Indexjmin, uint16_t * Indexjplus, OperatorType * const Jx)
{
  uint32_t i;
  uint16_t oadd[2];
  double j,jz,f;
  InitNamedOperator(Jx, "Jx");
  Jx->NFermion = NFermion;
  Jx->NWord = Min(Jx->NFermion, 1 + (Jx->NFermion - 1) / 8);
  InitOperatorFixedLength(Jx, 2, 8);
  // add terms for j=L-1/2
  j=(double)L-0.5;
  if(L>0) for (i = 0; i < 2 * L - 1; i++)
  {
    jz=(double)i - j;
    f=0.5*sqrt(j+jz+1)*sqrt(j-jz);
    oadd[1] = CrAnMask + Indexjmin[i+1];
    oadd[0] = Indexjmin[i];
    RealOperatorAddLaderNormalOrder(Jx, oadd, 2, f, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + Indexjmin[i];
    oadd[0] = Indexjmin[i+1];
    RealOperatorAddLaderNormalOrder(Jx, oadd, 2, f, &QDetNotOrderedFcFaBcBa);
  }
  // add terms for j=L+1/2
  j=(double)L+0.5;
  for (i = 0; i < 2 * L + 1; i++)
  {
    jz=(double)i - j;
    f=0.5*sqrt(j+jz+1)*sqrt(j-jz);
    oadd[1] = CrAnMask + Indexjplus[i+1];
    oadd[0] = Indexjplus[i];
    RealOperatorAddLaderNormalOrder(Jx, oadd, 2, f, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + Indexjplus[i];
    oadd[0] = Indexjplus[i+1];
    RealOperatorAddLaderNormalOrder(Jx, oadd, 2, f, &QDetNotOrderedFcFaBcBa);
  }
  OperatorChop(Jx, 0.001);
  return 0;
}

uint32_t CreateOperatorJy(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Jy)
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
    ComplexOperatorAddLaderNormalOrder(Jy, oadd, 2, 0.0,-0.5 * sqrt(l + m + 1) * sqrt(l - m), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + LindexUp[i];
    oadd[0] = LindexUp[i + 1];
    ComplexOperatorAddLaderNormalOrder(Jy, oadd, 2, 0.0, 0.5 * sqrt(l + m + 1) * sqrt(l - m), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + LindexDn[i + 1];
    oadd[0] = LindexDn[i];
    ComplexOperatorAddLaderNormalOrder(Jy, oadd, 2, 0.0,-0.5 * sqrt(l + m + 1) * sqrt(l - m), &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + LindexDn[i];
    oadd[0] = LindexDn[i + 1];
    ComplexOperatorAddLaderNormalOrder(Jy, oadd, 2, 0.0, 0.5 * sqrt(l + m + 1) * sqrt(l - m), &QDetNotOrderedFcFaBcBa);
  }
  for (i = 0; i < 2 * L + 1; i++)
  {
    oadd[1] = CrAnMask + LindexDn[i];
    oadd[0] = LindexUp[i];
    ComplexOperatorAddLaderNormalOrder(Jy, oadd, 2, 0.0, 0.5, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + LindexUp[i];
    oadd[0] = LindexDn[i];
    ComplexOperatorAddLaderNormalOrder(Jy, oadd, 2, 0.0, -0.5, &QDetNotOrderedFcFaBcBa);
  }
  OperatorChop(Jy, 0.001);
  return 0;
}
uint32_t CreateOperatorJyY(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Jy)
{
  return CreateOperatorJy(NFermion, L, LindexUp, LindexDn, Jy);
}
uint32_t CreateOperatorJyK(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Jy)
{
  return CreateOperatorInShellBasis(NFermion, L, 2, LindexUp, LindexDn, &CreateOperatorJy, Jy);
}
uint32_t CreateOperatorJyZ(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Jy)
{
  return CreateOperatorInShellBasis(NFermion, L, 1, LindexUp, LindexDn, &CreateOperatorJy, Jy);
}
uint32_t CreateOperatorJyP(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Jy)
{
  return CreateOperatorInShellBasis(NFermion, L, 3, LindexUp, LindexDn, &CreateOperatorJy, Jy);
}
uint32_t CreateOperatorJyjjz(uint32_t NFermion, uint32_t L, uint16_t * Indexjmin, uint16_t * Indexjplus, OperatorType * const Jy)
{
  uint32_t i;
  uint16_t oadd[2];
  double j,m,f;
  InitNamedOperator(Jy, "Jy");
  Jy->NFermion = NFermion;
  Jy->NWord = Min(Jy->NFermion, 1 + (Jy->NFermion - 1) / 8);
  InitOperatorFixedLength(Jy, 2, 8);
  MakeRealOperatorFixedLengthComplex(Jy, 2);
  Jy->QComplex = 1;
  // add terms for j=L-1/2
  j=(double)L-0.5;
  if (L>0) for (i = 0; i < 2 * L - 1; i++)
  {
    m=(double)i - j;
    f=0.5*sqrt(j+m+1)*sqrt(j-m);
    oadd[1] = CrAnMask + Indexjmin[i+1];
    oadd[0] = Indexjmin[i];
    ComplexOperatorAddLaderNormalOrder(Jy, oadd, 2, 0.0,-f, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + Indexjmin[i];
    oadd[0] = Indexjmin[i+1];
    ComplexOperatorAddLaderNormalOrder(Jy, oadd, 2, 0.0, f, &QDetNotOrderedFcFaBcBa);
  }
  // add terms for j=L+1/2
  j=(double)L+0.5;
  for (i = 0; i < 2 * L + 1; i++)
  {
    m=(double)i - j;
    f=0.5*sqrt(j+m+1)*sqrt(j-m);
    oadd[1] = CrAnMask + Indexjplus[i+1];
    oadd[0] = Indexjplus[i];
    ComplexOperatorAddLaderNormalOrder(Jy, oadd, 2, 0.0,-f, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + Indexjplus[i];
    oadd[0] = Indexjplus[i+1];
    ComplexOperatorAddLaderNormalOrder(Jy, oadd, 2, 0.0, f, &QDetNotOrderedFcFaBcBa);
  }
  OperatorChop(Jy, 0.001);
  return 0;
}

uint32_t CreateOperatorJz(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Jz)
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
uint32_t CreateOperatorJzY(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Jz)
{
  return CreateOperatorJz(NFermion, L, LindexUp, LindexDn, Jz);
}
uint32_t CreateOperatorJzK(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Jz)
{
  return CreateOperatorInShellBasis(NFermion, L, 2, LindexUp, LindexDn, &CreateOperatorJz, Jz);
}
uint32_t CreateOperatorJzZ(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Jz)
{
  return CreateOperatorInShellBasis(NFermion, L, 1, LindexUp, LindexDn, &CreateOperatorJz, Jz);
}
uint32_t CreateOperatorJzP(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Jz)
{
  return CreateOperatorInShellBasis(NFermion, L, 3, LindexUp, LindexDn, &CreateOperatorJz, Jz);
}
uint32_t CreateOperatorJzjjz(uint32_t NFermion, uint32_t L, uint16_t * Indexjmin, uint16_t * Indexjplus, OperatorType * const Jz)
{
  uint32_t i;
  uint16_t oadd[2];
  double j;
  InitNamedOperator(Jz, "Jz");
  Jz->NFermion = NFermion;
  Jz->NWord = Min(Jz->NFermion, 1 + (Jz->NFermion - 1) / 8);
  InitOperatorFixedLength(Jz, 2, 8);
  // add terms for j=L-1/2
  j=(double)L-0.5;
  for (i = 0; i < 2 * L ; i++)
  {
    oadd[1] = CrAnMask + Indexjmin[i];
    oadd[0] = Indexjmin[i];
    RealOperatorAddLaderNormalOrder(Jz, oadd, 2, (double)i - j, &QDetNotOrderedFcFaBcBa);
  }
  // add terms for j=L+1/2
  j=(double)L+0.5;
  for (i = 0; i < 2 * L + 2; i++)
  {
    oadd[1] = CrAnMask + Indexjplus[i];
    oadd[0] = Indexjplus[i];
    RealOperatorAddLaderNormalOrder(Jz, oadd, 2, (double)i - j, &QDetNotOrderedFcFaBcBa);
  }
  OperatorChop(Jz, 0.001);
  return 0;
}

uint32_t CreateOperatorJsqr(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Jsqr)
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
uint32_t CreateOperatorJsqrY(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Jsqr)
{
  return CreateOperatorJsqr(NFermion, L, LindexUp, LindexDn, Jsqr);
}
uint32_t CreateOperatorJsqrK(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Jsqr)
{
  return CreateOperatorInShellBasis(NFermion, L, 2, LindexUp, LindexDn, &CreateOperatorJsqr, Jsqr);
}
uint32_t CreateOperatorJsqrZ(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Jsqr)
{
  return CreateOperatorInShellBasis(NFermion, L, 1, LindexUp, LindexDn, &CreateOperatorJsqr, Jsqr);
}
uint32_t CreateOperatorJsqrP(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Jsqr)
{
  return CreateOperatorInShellBasis(NFermion, L, 3, LindexUp, LindexDn, &CreateOperatorJsqr, Jsqr);
}
uint32_t CreateOperatorJsqrjjz(uint32_t NFermion, uint32_t L, uint16_t * Indexjmin, uint16_t * Indexjplus, OperatorType * const Jsqr)
{
  uint32_t i1,i2;
  uint16_t oadd[4];
  double j1,j2,m1,m2;
  InitNamedOperator(Jsqr, "J^2");
  Jsqr->NFermion = NFermion;
  Jsqr->NWord = Min(Jsqr->NFermion, 1 + (Jsqr->NFermion - 1) / 8);
  InitOperatorFixedLength(Jsqr, 2, 8);
  InitOperatorFixedLength(Jsqr, 4, 8);
  //j1min j2min
  if(L>0)
  {
    j1=(double)L-0.5;
    j2=(double)L-0.5;
    for (i1=0;i1<2*L;i1++)
    {
      m1=(double)i1-j1;
      for (i2=0;i2<2*L;i2++)
      {
        m2=(double)i2-j2;
        oadd[3] = CrAnMask + Indexjmin[i1];
        oadd[2] = Indexjmin[i1];
        oadd[1] = CrAnMask + Indexjmin[i2];
        oadd[0] = Indexjmin[i2];
        RealOperatorAddLaderNormalOrder(Jsqr, oadd, 4, m1*m2, &QDetNotOrderedFcFaBcBa);
      }
    }
    for (i1=0;i1<2*L-1;i1++)
    {
      m1=(double)i1-j1;
      for (i2=0;i2<2*L-1;i2++)
      {
        m2=(double)i2-j2;
        oadd[3] = CrAnMask + Indexjmin[i1+1];
        oadd[2] = Indexjmin[i1];
        oadd[1] = CrAnMask + Indexjmin[i2];
        oadd[0] = Indexjmin[i2+1];
        RealOperatorAddLaderNormalOrder(Jsqr, oadd, 4, 0.5*sqrt(j1+m1+1)*sqrt(j1-m1)*sqrt(j2+m2+1)*sqrt(j2-m2), &QDetNotOrderedFcFaBcBa);
        oadd[3] = CrAnMask + Indexjmin[i1];
        oadd[2] = Indexjmin[i1+1];
        oadd[1] = CrAnMask + Indexjmin[i2+1];
        oadd[0] = Indexjmin[i2];
        RealOperatorAddLaderNormalOrder(Jsqr, oadd, 4, 0.5*sqrt(j1+m1+1)*sqrt(j1-m1)*sqrt(j2+m2+1)*sqrt(j2-m2), &QDetNotOrderedFcFaBcBa);
      }
    }
  }
  //j1min j2plus
  if(L>0)
  {
    j1=(double)L-0.5;
    j2=(double)L+0.5;
    for (i1=0;i1<2*L;i1++)
    {
      m1=(double)i1-j1;
      for (i2=0;i2<2*L+2;i2++)
      {
        m2=(double)i2-j2;
        oadd[3] = CrAnMask + Indexjmin[i1];
        oadd[2] = Indexjmin[i1];
        oadd[1] = CrAnMask + Indexjplus[i2];
        oadd[0] = Indexjplus[i2];
        RealOperatorAddLaderNormalOrder(Jsqr, oadd, 4, m1*m2, &QDetNotOrderedFcFaBcBa);
      }
    }
    for (i1=0;i1<2*L-1;i1++)
    {
      m1=(double)i1-j1;
      for (i2=0;i2<2*L+1;i2++)
      {
        m2=(double)i2-j2;
        oadd[3] = CrAnMask + Indexjmin[i1+1];
        oadd[2] = Indexjmin[i1];
        oadd[1] = CrAnMask + Indexjplus[i2];
        oadd[0] = Indexjplus[i2+1];
        RealOperatorAddLaderNormalOrder(Jsqr, oadd, 4, 0.5*sqrt(j1+m1+1)*sqrt(j1-m1)*sqrt(j2+m2+1)*sqrt(j2-m2), &QDetNotOrderedFcFaBcBa);
        oadd[3] = CrAnMask + Indexjmin[i1];
        oadd[2] = Indexjmin[i1+1];
        oadd[1] = CrAnMask + Indexjplus[i2+1];
        oadd[0] = Indexjplus[i2];
        RealOperatorAddLaderNormalOrder(Jsqr, oadd, 4, 0.5*sqrt(j1+m1+1)*sqrt(j1-m1)*sqrt(j2+m2+1)*sqrt(j2-m2), &QDetNotOrderedFcFaBcBa);
      }
    }
  }
  //j1plus j2min
  if(L>0)
  {
    j1=(double)L+0.5;
    j2=(double)L-0.5;
    for (i1=0;i1<2*L+2;i1++)
    {
      m1=(double)i1-j1;
      for (i2=0;i2<2*L;i2++)
      {
        m2=(double)i2-j2;
        oadd[3] = CrAnMask + Indexjplus[i1];
        oadd[2] = Indexjplus[i1];
        oadd[1] = CrAnMask + Indexjmin[i2];
        oadd[0] = Indexjmin[i2];
        RealOperatorAddLaderNormalOrder(Jsqr, oadd, 4, m1*m2, &QDetNotOrderedFcFaBcBa);
      }
    }
    for (i1=0;i1<2*L+1;i1++)
    {
      m1=(double)i1-j1;
      for (i2=0;i2<2*L-1;i2++)
      {
        m2=(double)i2-j2;
        oadd[3] = CrAnMask + Indexjplus[i1+1];
        oadd[2] = Indexjplus[i1];
        oadd[1] = CrAnMask + Indexjmin[i2];
        oadd[0] = Indexjmin[i2+1];
        RealOperatorAddLaderNormalOrder(Jsqr, oadd, 4, 0.5*sqrt(j1+m1+1)*sqrt(j1-m1)*sqrt(j2+m2+1)*sqrt(j2-m2), &QDetNotOrderedFcFaBcBa);
        oadd[3] = CrAnMask + Indexjplus[i1];
        oadd[2] = Indexjplus[i1+1];
        oadd[1] = CrAnMask + Indexjmin[i2+1];
        oadd[0] = Indexjmin[i2];
        RealOperatorAddLaderNormalOrder(Jsqr, oadd, 4, 0.5*sqrt(j1+m1+1)*sqrt(j1-m1)*sqrt(j2+m2+1)*sqrt(j2-m2), &QDetNotOrderedFcFaBcBa);
      }
    }
  }
  //j1plus j2plus
  {
    j1=(double)L+0.5;
    j2=(double)L+0.5;
    for (i1=0;i1<2*L+2;i1++)
    {
      m1=(double)i1-j1;
      for (i2=0;i2<2*L+2;i2++)
      {
        m2=(double)i2-j2;
        oadd[3] = CrAnMask + Indexjplus[i1];
        oadd[2] = Indexjplus[i1];
        oadd[1] = CrAnMask + Indexjplus[i2];
        oadd[0] = Indexjplus[i2];
        RealOperatorAddLaderNormalOrder(Jsqr, oadd, 4, m1*m2, &QDetNotOrderedFcFaBcBa);
      }
    }
    for (i1=0;i1<2*L+1;i1++)
    {
      m1=(double)i1-j1;
      for (i2=0;i2<2*L+1;i2++)
      {
        m2=(double)i2-j2;
        oadd[3] = CrAnMask + Indexjplus[i1+1];
        oadd[2] = Indexjplus[i1];
        oadd[1] = CrAnMask + Indexjplus[i2];
        oadd[0] = Indexjplus[i2+1];
        RealOperatorAddLaderNormalOrder(Jsqr, oadd, 4, 0.5*sqrt(j1+m1+1)*sqrt(j1-m1)*sqrt(j2+m2+1)*sqrt(j2-m2), &QDetNotOrderedFcFaBcBa);
        oadd[3] = CrAnMask + Indexjplus[i1];
        oadd[2] = Indexjplus[i1+1];
        oadd[1] = CrAnMask + Indexjplus[i2+1];
        oadd[0] = Indexjplus[i2];
        RealOperatorAddLaderNormalOrder(Jsqr, oadd, 4, 0.5*sqrt(j1+m1+1)*sqrt(j1-m1)*sqrt(j2+m2+1)*sqrt(j2-m2), &QDetNotOrderedFcFaBcBa);
      }
    }
  }
  
  
  OperatorChop(Jsqr, 0.001);
  return 0;
}

uint32_t CreateOperatorJplus(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Jplus)
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
uint32_t CreateOperatorJplusY(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Jplus)
{
  return CreateOperatorJplus(NFermion, L, LindexUp, LindexDn, Jplus);
}
uint32_t CreateOperatorJplusK(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Jplus)
{
  return CreateOperatorInShellBasis(NFermion, L, 2, LindexUp, LindexDn, &CreateOperatorJplus, Jplus);
}
uint32_t CreateOperatorJplusZ(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Jplus)
{
  return CreateOperatorInShellBasis(NFermion, L, 1, LindexUp, LindexDn, &CreateOperatorJplus, Jplus);
}
uint32_t CreateOperatorJplusP(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Jplus)
{
  return CreateOperatorInShellBasis(NFermion, L, 3, LindexUp, LindexDn, &CreateOperatorJplus, Jplus);
}
uint32_t CreateOperatorJplusjjz(uint32_t NFermion, uint32_t L, uint16_t * Indexjmin, uint16_t * Indexjplus, OperatorType * const Jplus)
{
  uint32_t i;
  uint16_t oadd[2];
  double j,m;
  InitNamedOperator(Jplus, "J^+");
  Jplus->NFermion = NFermion;
  Jplus->NWord = Min(Jplus->NFermion, 1 + (Jplus->NFermion - 1) / 8);
  InitOperatorFixedLength(Jplus, 2, 8);
  // add terms for j=L-1/2
  j=(double)L-0.5;
  if (L>0) for (i = 0; i < 2 * L - 1; i++)
  {
    m=(double)i - j;
    oadd[1] = CrAnMask + Indexjmin[i+1];
    oadd[0] = Indexjmin[i];
    RealOperatorAddLaderNormalOrder(Jplus, oadd, 2, sqrt(j+m+1)*sqrt(j-m), &QDetNotOrderedFcFaBcBa);
  }
  // add terms for j=L+1/2
  j=(double)L+0.5;
  for (i = 0; i < 2 * L + 1; i++)
  {
    m=(double)i - j;
    oadd[1] = CrAnMask + Indexjplus[i+1];
    oadd[0] = Indexjplus[i];
    RealOperatorAddLaderNormalOrder(Jplus, oadd, 2, sqrt(j+m+1)*sqrt(j-m), &QDetNotOrderedFcFaBcBa);
  }
  OperatorChop(Jplus, 0.001);
  return 0;
}

uint32_t CreateOperatorJmin(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Jmin)
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
uint32_t CreateOperatorJminY(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Jmin)
{
  return CreateOperatorJmin(NFermion, L, LindexUp, LindexDn, Jmin);
}
uint32_t CreateOperatorJminK(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Jmin)
{
  return CreateOperatorInShellBasis(NFermion, L, 2, LindexUp, LindexDn, &CreateOperatorJmin, Jmin);
}
uint32_t CreateOperatorJminZ(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Jmin)
{
  return CreateOperatorInShellBasis(NFermion, L, 1, LindexUp, LindexDn, &CreateOperatorJmin, Jmin);
}
uint32_t CreateOperatorJminP(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Jmin)
{
  return CreateOperatorInShellBasis(NFermion, L, 3, LindexUp, LindexDn, &CreateOperatorJmin, Jmin);
}
uint32_t CreateOperatorJminjjz(uint32_t NFermion, uint32_t L, uint16_t * Indexjmin, uint16_t * Indexjplus, OperatorType * const Jmin)
{
  uint32_t i;
  uint16_t oadd[2];
  double j,m;
  InitNamedOperator(Jmin, "J^-");
  Jmin->NFermion = NFermion;
  Jmin->NWord = Min(Jmin->NFermion, 1 + (Jmin->NFermion - 1) / 8);
  InitOperatorFixedLength(Jmin, 2, 8);
  // add terms for j=L-1/2
  j=(double)L-0.5;
  if (L>0) for (i = 0; i < 2 * L - 1; i++)
  {
    m=(double)i - j;
    oadd[1] = CrAnMask + Indexjmin[i];
    oadd[0] = Indexjmin[i+1];
    RealOperatorAddLaderNormalOrder(Jmin, oadd, 2, sqrt(j+m+1)*sqrt(j-m), &QDetNotOrderedFcFaBcBa);
  }
  // add terms for j=L+1/2
  j=(double)L+0.5;
  for (i = 0; i < 2 * L + 1; i++)
  {
    m=(double)i - j;
    oadd[1] = CrAnMask + Indexjplus[i];
    oadd[0] = Indexjplus[i+1];
    RealOperatorAddLaderNormalOrder(Jmin, oadd, 2, sqrt(j+m+1)*sqrt(j-m), &QDetNotOrderedFcFaBcBa);
  }
  OperatorChop(Jmin, 0.001);
  return 0;
}


/*  CreateOperatorT*** creates the operators Tx, Ty, Tz, Tplus, Tmin or T^2
 *  Input:
 *    NFermion    - The number of Fermions in the basis
 *    L           - The total angular momentum (2L+1) <= NFermion
 *    *LindexUp    - An array of legnth 2L+1 indexing the ml values from -L to L with spin up
 *    *LindexDn    - An array of legnth 2L+1 indexing the ml values from -L to L with spin down
 *    The operators T*** are not initialized on input
 */


uint32_t CreateOperatorTx(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Tx)
{
  uint32_t i,j;
  int    mCr, mAn, lCr, lAn, mk;
  double val;
  lCr = (int)L;
  lAn = (int)L;
  uint16_t oadd[2];
  InitNamedOperator(Tx, "Tx");
  Tx->NFermion = NFermion;
  Tx->NWord = Min(Tx->NFermion, 1 + (Tx->NFermion - 1) / 8);
  InitOperatorFixedLength(Tx, 2, 8);
  for(i=0;i<2*L+1;i++)
  {
    // C1 sz
    mk  = 1;
    mAn = (int)i - L;
    mCr = mAn+mk;
    val = (double)SlaterCoefficientC(2, lCr, mCr, lAn, mAn);
    if (Abs(val) > DBL_EPSILON)
    {
      j = i + mk;
      oadd[1] = CrAnMask + LindexUp[j];
      oadd[0] = LindexUp[i];
      RealOperatorAddLaderNormalOrder(Tx, oadd, 2,  0.5*sqrt(1.5)*val, &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + LindexDn[j];
      oadd[0] = LindexDn[i];
      RealOperatorAddLaderNormalOrder(Tx, oadd, 2, -0.5*sqrt(1.5)*val, &QDetNotOrderedFcFaBcBa);
    }
    // C0 s+ / s-
    mk  = 0;
    mAn = (int)i - L;
    mCr = mAn+mk;
    val = (double)SlaterCoefficientC(2, lCr, mCr, lAn, mAn);
    if (Abs(val) > DBL_EPSILON)
    {
      j = i + mk;
      oadd[1] = CrAnMask + LindexUp[j];
      oadd[0] = LindexDn[i];
      RealOperatorAddLaderNormalOrder(Tx, oadd, 2, 0.5*val, &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + LindexDn[j];
      oadd[0] = LindexUp[i];
      RealOperatorAddLaderNormalOrder(Tx, oadd, 2,  0.5*val, &QDetNotOrderedFcFaBcBa);
    }
    // C2 s-
    mk  = 2;
    mAn = (int)i - L;
    mCr = mAn+mk;
    val = (double)SlaterCoefficientC(2, lCr, mCr, lAn, mAn);
    if (Abs(val) > DBL_EPSILON)
    {
      j = i + mk;
      oadd[1] = CrAnMask + LindexDn[j];
      oadd[0] = LindexUp[i];
      RealOperatorAddLaderNormalOrder(Tx, oadd, 2,  -sqrt(1.5)*val, &QDetNotOrderedFcFaBcBa);
    }
    // C-1 sz
    mk  = -1;
    mAn = (int)i - L;
    mCr = mAn+mk;
    val = (double)SlaterCoefficientC(2, lCr, mCr, lAn, mAn);
    if (Abs(val) > DBL_EPSILON)
    {
      j = i + mk;
      oadd[1] = CrAnMask + LindexUp[j];
      oadd[0] = LindexUp[i];
      RealOperatorAddLaderNormalOrder(Tx, oadd, 2, -0.5*sqrt(1.5)*val, &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + LindexDn[j];
      oadd[0] = LindexDn[i];
      RealOperatorAddLaderNormalOrder(Tx, oadd, 2,  0.5*sqrt(1.5)*val, &QDetNotOrderedFcFaBcBa);
    }
    // C-2 s+
    mk  = -2;
    mAn = (int)i - L;
    mCr = mAn+mk;
    val = (double)SlaterCoefficientC(2, lCr, mCr, lAn, mAn);
    if (Abs(val) > DBL_EPSILON)
    {
      j = i + mk;
      oadd[1] = CrAnMask + LindexUp[j];
      oadd[0] = LindexDn[i];
      RealOperatorAddLaderNormalOrder(Tx, oadd, 2, -sqrt(1.5)*val, &QDetNotOrderedFcFaBcBa);
    }
  }
  return 0;
}
uint32_t CreateOperatorTxY(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Tx)
{
  return CreateOperatorTx(NFermion, L, LindexUp, LindexDn, Tx);
}
uint32_t CreateOperatorTxK(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Tx)
{
  return CreateOperatorInShellBasis(NFermion, L, 2, LindexUp, LindexDn, &CreateOperatorTx, Tx);
}
uint32_t CreateOperatorTxZ(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Tx)
{
  return CreateOperatorInShellBasis(NFermion, L, 1, LindexUp, LindexDn, &CreateOperatorTx, Tx);
}
uint32_t CreateOperatorTxP(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Tx)
{
  return CreateOperatorInShellBasis(NFermion, L, 3, LindexUp, LindexDn, &CreateOperatorTx, Tx);
}
uint32_t CreateOperatorTxjjz(uint32_t NFermion, uint32_t L, uint16_t * Lindexjmin, uint16_t * Lindexjplus, OperatorType * const Tx)
{
  return CreateOperatorInShellBasis(NFermion, L, 4, Lindexjmin, Lindexjplus, &CreateOperatorTx, Tx);
}

uint32_t CreateOperatorTy(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Ty)
{
  uint32_t i,j;
  int    mCr, mAn, lCr, lAn, mk;
  double val;
  lCr = (int)L;
  lAn = (int)L;
  uint16_t oadd[2];
  InitNamedOperator(Ty, "Ty");
  Ty->NFermion = NFermion;
  Ty->NWord = Min(Ty->NFermion, 1 + (Ty->NFermion - 1) / 8);
  InitOperatorFixedLength(Ty, 2, 8);
  MakeRealOperatorFixedLengthComplex(Ty, 2);
  for(i=0;i<2*L+1;i++)
  {
    // C1 sz
    mk  = 1;
    mAn = (int)i - L;
    mCr = mAn+mk;
    val = (double)SlaterCoefficientC(2, lCr, mCr, lAn, mAn);
    if (Abs(val) > DBL_EPSILON)
    {
      j = i + mk;
      oadd[1] = CrAnMask + LindexUp[j];
      oadd[0] = LindexUp[i];
      ComplexOperatorAddLaderNormalOrder(Ty, oadd, 2, 0, -0.5*sqrt(1.5)*val, &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + LindexDn[j];
      oadd[0] = LindexDn[i];
      ComplexOperatorAddLaderNormalOrder(Ty, oadd, 2, 0, 0.5*sqrt(1.5)*val, &QDetNotOrderedFcFaBcBa);
    }
    // C0 s+ / s-
    mk  = 0;
    mAn = (int)i - L;
    mCr = mAn+mk;
    val = (double)SlaterCoefficientC(2, lCr, mCr, lAn, mAn);
    if (Abs(val) > DBL_EPSILON)
    {
      j = i + mk;
      oadd[1] = CrAnMask + LindexUp[j];
      oadd[0] = LindexDn[i];
      ComplexOperatorAddLaderNormalOrder(Ty, oadd, 2, 0, -0.5*val, &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + LindexDn[j];
      oadd[0] = LindexUp[i];
      ComplexOperatorAddLaderNormalOrder(Ty, oadd, 2,  0, 0.5*val, &QDetNotOrderedFcFaBcBa);
    }
    // C2 s-
    mk  = 2;
    mAn = (int)i - L;
    mCr = mAn+mk;
    val = (double)SlaterCoefficientC(2, lCr, mCr, lAn, mAn);
    if (Abs(val) > DBL_EPSILON)
    {
      j = i + mk;
      oadd[1] = CrAnMask + LindexDn[j];
      oadd[0] = LindexUp[i];
      ComplexOperatorAddLaderNormalOrder(Ty, oadd, 2, 0, sqrt(1.5)*val, &QDetNotOrderedFcFaBcBa);
    }
    // C-1 sz
    mk  = -1;
    mAn = (int)i - L;
    mCr = mAn+mk;
    val = (double)SlaterCoefficientC(2, lCr, mCr, lAn, mAn);
    if (Abs(val) > DBL_EPSILON)
    {
      j = i + mk;
      oadd[1] = CrAnMask + LindexUp[j];
      oadd[0] = LindexUp[i];
      ComplexOperatorAddLaderNormalOrder(Ty, oadd, 2, 0, -0.5*sqrt(1.5)*val, &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + LindexDn[j];
      oadd[0] = LindexDn[i];
      ComplexOperatorAddLaderNormalOrder(Ty, oadd, 2, 0,  0.5*sqrt(1.5)*val, &QDetNotOrderedFcFaBcBa);
    }
    // C-2 s+
    mk  = -2;
    mAn = (int)i - L;
    mCr = mAn+mk;
    val = (double)SlaterCoefficientC(2, lCr, mCr, lAn, mAn);
    if (Abs(val) > DBL_EPSILON)
    {
      j = i + mk;
      oadd[1] = CrAnMask + LindexUp[j];
      oadd[0] = LindexDn[i];
      ComplexOperatorAddLaderNormalOrder(Ty, oadd, 2, 0, -sqrt(1.5)*val, &QDetNotOrderedFcFaBcBa);
    }
  }
  return 0;
}
uint32_t CreateOperatorTyY(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Ty)
{
  return CreateOperatorTy(NFermion, L, LindexUp, LindexDn, Ty);
}
uint32_t CreateOperatorTyK(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Ty)
{
  return CreateOperatorInShellBasis(NFermion, L, 2, LindexUp, LindexDn, &CreateOperatorTy, Ty);
}
uint32_t CreateOperatorTyZ(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Ty)
{
  return CreateOperatorInShellBasis(NFermion, L, 1, LindexUp, LindexDn, &CreateOperatorTy, Ty);
}
uint32_t CreateOperatorTyP(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Ty)
{
  return CreateOperatorInShellBasis(NFermion, L, 3, LindexUp, LindexDn, &CreateOperatorTy, Ty);
}
uint32_t CreateOperatorTyjjz(uint32_t NFermion, uint32_t L, uint16_t * Lindexjmin, uint16_t * Lindexjplus, OperatorType * const Ty)
{
  return CreateOperatorInShellBasis(NFermion, L, 4, Lindexjmin, Lindexjplus, &CreateOperatorTy, Ty);
}

uint32_t CreateOperatorTz(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Tz)
{
  uint32_t i,j;
  int    mCr, mAn, lCr, lAn, mk;
  double val;
  lCr = (int)L;
  lAn = (int)L;
  uint16_t oadd[2];
  InitNamedOperator(Tz, "Tz");
  Tz->NFermion = NFermion;
  Tz->NWord = Min(Tz->NFermion, 1 + (Tz->NFermion - 1) / 8);
  InitOperatorFixedLength(Tz, 2, 8);
  for(i=0;i<2*L+1;i++)
  {
    // C0 sz
    mk  = 0;
    mAn = (int)i - L;
    mCr = mAn+mk;
    val = (double)SlaterCoefficientC(2, lCr, mCr, lAn, mAn);
    if (Abs(val) > DBL_EPSILON)
    {
      j = i + mk;
      oadd[1] = CrAnMask + LindexUp[j];
      oadd[0] = LindexUp[i];
      RealOperatorAddLaderNormalOrder(Tz, oadd, 2, -val, &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + LindexDn[j];
      oadd[0] = LindexDn[i];
      RealOperatorAddLaderNormalOrder(Tz, oadd, 2,  val, &QDetNotOrderedFcFaBcBa);
    }
    // C-1 s+
    mk  = -1;
    mAn = (int)i - L;
    mCr = mAn+mk;
    val = (double)SlaterCoefficientC(2, lCr, mCr, lAn, mAn);
    if (Abs(val) > DBL_EPSILON)
    {
      j = i + mk;
      oadd[1] = CrAnMask + LindexUp[j];
      oadd[0] = LindexDn[i];
      RealOperatorAddLaderNormalOrder(Tz, oadd, 2, -0.5*sqrt(6.)*val, &QDetNotOrderedFcFaBcBa);
    }
    // C1 s-
    mk  = 1;
    mAn = (int)i - L;
    mCr = mAn+mk;
    val = (double)SlaterCoefficientC(2, lCr, mCr, lAn, mAn);
    if (Abs(val) > DBL_EPSILON)
    {
      j = i + mk;
      oadd[1] = CrAnMask + LindexDn[j];
      oadd[0] = LindexUp[i];
      RealOperatorAddLaderNormalOrder(Tz, oadd, 2,  0.5*sqrt(6.)*val, &QDetNotOrderedFcFaBcBa);
    }
  }
  return 0;
}
uint32_t CreateOperatorTzY(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Tz)
{
  return CreateOperatorTz(NFermion, L, LindexUp, LindexDn, Tz);
}
uint32_t CreateOperatorTzK(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Tz)
{
  return CreateOperatorInShellBasis(NFermion, L, 2, LindexUp, LindexDn, &CreateOperatorTz, Tz);
}
uint32_t CreateOperatorTzZ(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Tz)
{
  return CreateOperatorInShellBasis(NFermion, L, 1, LindexUp, LindexDn, &CreateOperatorTz, Tz);
}
uint32_t CreateOperatorTzP(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Tz)
{
  return CreateOperatorInShellBasis(NFermion, L, 3, LindexUp, LindexDn, &CreateOperatorTz, Tz);
}
uint32_t CreateOperatorTzjjz(uint32_t NFermion, uint32_t L, uint16_t * Lindexjmin, uint16_t * Lindexjplus, OperatorType * const Tz)
{
  return CreateOperatorInShellBasis(NFermion, L, 4, Lindexjmin, Lindexjplus, &CreateOperatorTz, Tz);
}

uint32_t CreateOperatorTplus(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Tplus)
{
  uint32_t i,j;
  int    mCr, mAn, lCr, lAn, mk;
  double val;
  lCr = (int)L;
  lAn = (int)L;
  uint16_t oadd[2];
  InitNamedOperator(Tplus, "Tplus");
  Tplus->NFermion = NFermion;
  Tplus->NWord = Min(Tplus->NFermion, 1 + (Tplus->NFermion - 1) / 8);
  InitOperatorFixedLength(Tplus, 2, 8);
  for(i=0;i<2*L+1;i++)
  {
    // C1 sz
    mk  = 1;
    mAn = (int)i - L;
    mCr = mAn+mk;
    val = (double)SlaterCoefficientC(2, lCr, mCr, lAn, mAn);
    if (Abs(val) > DBL_EPSILON)
    {
      j = i + mk;
      oadd[1] = CrAnMask + LindexUp[j];
      oadd[0] = LindexUp[i];
      RealOperatorAddLaderNormalOrder(Tplus, oadd, 2,  sqrt(1.5)*val, &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + LindexDn[j];
      oadd[0] = LindexDn[i];
      RealOperatorAddLaderNormalOrder(Tplus, oadd, 2, -sqrt(1.5)*val, &QDetNotOrderedFcFaBcBa);
    }
    // C0 s+
    mk  = 0;
    mAn = (int)i - L;
    mCr = mAn+mk;
    val = (double)SlaterCoefficientC(2, lCr, mCr, lAn, mAn);
    if (Abs(val) > DBL_EPSILON)
    {
      j = i + mk;
      oadd[1] = CrAnMask + LindexUp[j];
      oadd[0] = LindexDn[i];
      RealOperatorAddLaderNormalOrder(Tplus, oadd, 2, val, &QDetNotOrderedFcFaBcBa);
    }
    // C2 s-
    mk  = 2;
    mAn = (int)i - L;
    mCr = mAn+mk;
    val = (double)SlaterCoefficientC(2, lCr, mCr, lAn, mAn);
    if (Abs(val) > DBL_EPSILON)
    {
      j = i + mk;
      oadd[1] = CrAnMask + LindexDn[j];
      oadd[0] = LindexUp[i];
      RealOperatorAddLaderNormalOrder(Tplus, oadd, 2,  -2*sqrt(1.5)*val, &QDetNotOrderedFcFaBcBa);
    }
  }
  return 0;
}
uint32_t CreateOperatorTplusY(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Tplus)
{
  return CreateOperatorTplus(NFermion, L, LindexUp, LindexDn, Tplus);
}
uint32_t CreateOperatorTplusK(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Tplus)
{
  return CreateOperatorInShellBasis(NFermion, L, 2, LindexUp, LindexDn, &CreateOperatorTplus, Tplus);
}
uint32_t CreateOperatorTplusZ(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Tplus)
{
  return CreateOperatorInShellBasis(NFermion, L, 1, LindexUp, LindexDn, &CreateOperatorTplus, Tplus);
}
uint32_t CreateOperatorTplusP(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Tplus)
{
  return CreateOperatorInShellBasis(NFermion, L, 3, LindexUp, LindexDn, &CreateOperatorTplus, Tplus);
}
uint32_t CreateOperatorTplusjjz(uint32_t NFermion, uint32_t L, uint16_t * Lindexjmin, uint16_t * Lindexjplus, OperatorType * const Tplus)
{
  return CreateOperatorInShellBasis(NFermion, L, 4, Lindexjmin, Lindexjplus, &CreateOperatorTplus, Tplus);
}

uint32_t CreateOperatorTmin(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Tmin)
{
  uint32_t i,j;
  int    mCr, mAn, lCr, lAn, mk;
  double val;
  lCr = (int)L;
  lAn = (int)L;
  uint16_t oadd[2];
  InitNamedOperator(Tmin, "Tmin");
  Tmin->NFermion = NFermion;
  Tmin->NWord = Min(Tmin->NFermion, 1 + (Tmin->NFermion - 1) / 8);
  InitOperatorFixedLength(Tmin, 2, 8);
  for(i=0;i<2*L+1;i++)
  {
    // C-1 sz
    mk  = -1;
    mAn = (int)i - L;
    mCr = mAn+mk;
    val = (double)SlaterCoefficientC(2, lCr, mCr, lAn, mAn);
    if (Abs(val) > DBL_EPSILON)
    {
      j = i + mk;
      oadd[1] = CrAnMask + LindexUp[j];
      oadd[0] = LindexUp[i];
      RealOperatorAddLaderNormalOrder(Tmin, oadd, 2, -sqrt(1.5)*val, &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + LindexDn[j];
      oadd[0] = LindexDn[i];
      RealOperatorAddLaderNormalOrder(Tmin, oadd, 2,  sqrt(1.5)*val, &QDetNotOrderedFcFaBcBa);
    }
    // C-2 s+
    mk  = -2;
    mAn = (int)i - L;
    mCr = mAn+mk;
    val = (double)SlaterCoefficientC(2, lCr, mCr, lAn, mAn);
    if (Abs(val) > DBL_EPSILON)
    {
      j = i + mk;
      oadd[1] = CrAnMask + LindexUp[j];
      oadd[0] = LindexDn[i];
      RealOperatorAddLaderNormalOrder(Tmin, oadd, 2, -2*sqrt(1.5)*val, &QDetNotOrderedFcFaBcBa);
    }
    // C0 s-
    mk  = 0;
    mAn = (int)i - L;
    mCr = mAn+mk;
    val = (double)SlaterCoefficientC(2, lCr, mCr, lAn, mAn);
    if (Abs(val) > DBL_EPSILON)
    {
      j = i + mk;
      oadd[1] = CrAnMask + LindexDn[j];
      oadd[0] = LindexUp[i];
      RealOperatorAddLaderNormalOrder(Tmin, oadd, 2,  val, &QDetNotOrderedFcFaBcBa);
    }
  }
  return 0;
}
uint32_t CreateOperatorTminY(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Tmin)
{
  return CreateOperatorTmin(NFermion, L, LindexUp, LindexDn, Tmin);
}
uint32_t CreateOperatorTminK(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Tmin)
{
  return CreateOperatorInShellBasis(NFermion, L, 2, LindexUp, LindexDn, &CreateOperatorTmin, Tmin);
}
uint32_t CreateOperatorTminZ(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Tmin)
{
  return CreateOperatorInShellBasis(NFermion, L, 1, LindexUp, LindexDn, &CreateOperatorTmin, Tmin);
}
uint32_t CreateOperatorTminP(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Tmin)
{
  return CreateOperatorInShellBasis(NFermion, L, 3, LindexUp, LindexDn, &CreateOperatorTmin, Tmin);
}
uint32_t CreateOperatorTminjjz(uint32_t NFermion, uint32_t L, uint16_t * Lindexjmin, uint16_t * Lindexjplus, OperatorType * const Tmin)
{
  return CreateOperatorInShellBasis(NFermion, L, 4, Lindexjmin, Lindexjplus, &CreateOperatorTmin, Tmin);
}

uint32_t CreateOperatorTsqr(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Tsqr)
{
  OperatorType Otemp;
  InitNamedOperator(Tsqr, "T^2");
  Tsqr->NFermion = NFermion;
  Tsqr->NWord = Min(Tsqr->NFermion, 1 + (Tsqr->NFermion - 1) / 8);
  InitOperatorFixedLength(Tsqr, 2, 8);
  InitOperatorFixedLength(Tsqr, 4, 8);
  CreateOperatorTx(NFermion,L,LindexUp,LindexDn,&Otemp);
  OperatorMultiplyNormalOrder(Tsqr, Otemp, Otemp, &QDetNotOrderedFcFaBcBa);
  FreeOperator(&Otemp);
  CreateOperatorTy(NFermion,L,LindexUp,LindexDn,&Otemp);
  OperatorMultiplyNormalOrder(Tsqr, Otemp, Otemp, &QDetNotOrderedFcFaBcBa);
  FreeOperator(&Otemp);
  CreateOperatorTz(NFermion,L,LindexUp,LindexDn,&Otemp);
  OperatorMultiplyNormalOrder(Tsqr, Otemp, Otemp, &QDetNotOrderedFcFaBcBa);
  FreeOperator(&Otemp);
  return 0;
}
uint32_t CreateOperatorTsqrY(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Tsqr)
{
  return CreateOperatorTsqr(NFermion, L, LindexUp, LindexDn, Tsqr);
}
uint32_t CreateOperatorTsqrK(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Tsqr)
{
  return CreateOperatorInShellBasis(NFermion, L, 2, LindexUp, LindexDn, &CreateOperatorTsqr, Tsqr);
}
uint32_t CreateOperatorTsqrZ(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Tsqr)
{
  return CreateOperatorInShellBasis(NFermion, L, 1, LindexUp, LindexDn, &CreateOperatorTsqr, Tsqr);
}
uint32_t CreateOperatorTsqrP(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Tsqr)
{
  return CreateOperatorInShellBasis(NFermion, L, 3, LindexUp, LindexDn, &CreateOperatorTsqr, Tsqr);
}
uint32_t CreateOperatorTsqrjjz(uint32_t NFermion, uint32_t L, uint16_t * Lindexjmin, uint16_t * Lindexjplus, OperatorType * const Tsqr)
{
  return CreateOperatorInShellBasis(NFermion, L, 4, Lindexjmin, Lindexjplus, &CreateOperatorTsqr, Tsqr);
}


uint32_t CreateOperatorldots(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const ldots)
{
  uint32_t i;
  double mi, l;
  l = (double)L;
  uint16_t oadd[2];
  InitNamedOperator(ldots, "l.s");
  ldots->NFermion = NFermion;
  ldots->NWord = Min(ldots->NFermion, 1 + (ldots->NFermion - 1) / 8);
  InitOperatorFixedLength(ldots, 2, 8);
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
uint32_t CreateOperatorldotsY(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const ldots)
{
  return CreateOperatorldots(NFermion, L, LindexUp, LindexDn, ldots);
}
uint32_t CreateOperatorldotsK(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const ldots)
{
  return CreateOperatorInShellBasis(NFermion, L, 2, LindexUp, LindexDn, &CreateOperatorldots, ldots);
}
uint32_t CreateOperatorldotsZ(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const ldots)
{
  return CreateOperatorInShellBasis(NFermion, L, 1, LindexUp, LindexDn, &CreateOperatorldots, ldots);
}
uint32_t CreateOperatorldotsP(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const ldots)
{
  return CreateOperatorInShellBasis(NFermion, L, 3, LindexUp, LindexDn, &CreateOperatorldots, ldots);
}
uint32_t CreateOperatorldotsjjz(uint32_t NFermion, uint32_t L, uint16_t * Indexjmin, uint16_t * Indexjplus, OperatorType * const ldots)
{
  uint32_t i;
  uint16_t oadd[2];
  double fac;
  InitNamedOperator(ldots, "l.s");
  ldots->NFermion = NFermion;
  ldots->NWord = Min(ldots->NFermion, 1 + (ldots->NFermion - 1) / 8);
  InitOperatorFixedLength(ldots, 2, 8);
  // add terms for j=L-1/2
  fac=(-(double)L-1)*0.5;
  for (i = 0; i < 2 * L ; i++)
  {
    oadd[1] = CrAnMask + Indexjmin[i];
    oadd[0] = Indexjmin[i];
    RealOperatorAddLaderNormalOrder(ldots, oadd, 2, fac, &QDetNotOrderedFcFaBcBa);
  }
  // add terms for j=L+1/2
  fac=((double)L)*0.5;
  for (i = 0; i < 2 * L + 2; i++)
  {
    oadd[1] = CrAnMask + Indexjplus[i];
    oadd[0] = Indexjplus[i];
    RealOperatorAddLaderNormalOrder(ldots, oadd, 2, fac, &QDetNotOrderedFcFaBcBa);
  }
  OperatorChop(ldots, 0.001);
  return 0;
}

uint32_t CreateOperatorCF(uint32_t const NFermion, uint32_t LCr, uint32_t LAn, uint16_t * LCrindexUp, uint16_t * LCrindexDn, uint16_t * LAnindexUp, uint16_t * LAnindexDn, HarmonicsExpansion const Akm, OperatorType * const CF)
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
uint32_t CreateOperatorCFY(uint32_t const NFermion, uint32_t LCr, uint32_t LAn, uint16_t * LCrindexUp, uint16_t * LCrindexDn, uint16_t * LAnindexUp, uint16_t * LAnindexDn, HarmonicsExpansion const Akm, OperatorType * const CF)
{
  return CreateOperatorCF(NFermion, LCr, LAn, LCrindexUp,LCrindexDn, LAnindexUp, LAnindexDn, Akm, CF);
}
uint32_t CreateOperatorCFK(uint32_t const NFermion, uint32_t LCr, uint32_t LAn, uint16_t * LCrindexUp, uint16_t * LCrindexDn, uint16_t * LAnindexUp, uint16_t * LAnindexDn, HarmonicsExpansion const Akm, OperatorType * const CF)
{
  return CreateOperatorCFInShellBasis(NFermion, 2, LCr, LAn, LCrindexUp, LCrindexDn, LAnindexUp, LAnindexDn, Akm, CF);
}
uint32_t CreateOperatorCFZ(uint32_t const NFermion, uint32_t LCr, uint32_t LAn, uint16_t * LCrindexUp, uint16_t * LCrindexDn, uint16_t * LAnindexUp, uint16_t * LAnindexDn, HarmonicsExpansion const Akm, OperatorType * const CF)
{
  return CreateOperatorCFInShellBasis(NFermion, 1, LCr, LAn, LCrindexUp, LCrindexDn, LAnindexUp, LAnindexDn, Akm, CF);
}
uint32_t CreateOperatorCFP(uint32_t const NFermion, uint32_t LCr, uint32_t LAn, uint16_t * LCrindexUp, uint16_t * LCrindexDn, uint16_t * LAnindexUp, uint16_t * LAnindexDn, HarmonicsExpansion const Akm, OperatorType * const CF)
{
  return CreateOperatorCFInShellBasis(NFermion, 3, LCr, LAn, LCrindexUp, LCrindexDn, LAnindexUp, LAnindexDn, Akm, CF);
}
uint32_t CreateOperatorCFjjz(uint32_t const NFermion, uint32_t LCr, uint32_t LAn, uint16_t * LCrindexmin, uint16_t * LCrindexplus, uint16_t * LAnindexmin, uint16_t * LAnindexplus, HarmonicsExpansion const Akm, OperatorType * const CF)
{
  return CreateOperatorCFInShellBasis(NFermion, 4, LCr, LAn, LCrindexmin, LCrindexplus, LAnindexmin, LAnindexplus, Akm, CF);
}
uint32_t CreateOperatorCFjjz(uint32_t const NFermion, uint32_t twoJ, uint16_t * JCrindexmin, uint16_t * JAnindexmin, HarmonicsExpansion const Akm, OperatorType * const CF)
{
  printf("not yet implemented\n");
  //TODO:
  return 0;
}


/*  CreateOperatorQ*** creates the operators Qxx, Qyy, Qzz, Qxy, Qxz, Qyz as special cases of the Crystal-Field operator
 *  Qij = 3r_i r_j - r . r delta_{i,j}
 *  Input:
 *    NFermion    - The number of Fermions in the basis
 *    L           - The total angular momentum (2L+1) <= NFermion
 *    *LindexUp    - An array of legnth 2L+1 indexing the ml values from -L to L with spin up
 *    *LindexDn    - An array of legnth 2L+1 indexing the ml values from -L to L with spin down
 *    The operators Q*** are not initialized on input
 */

uint32_t CreateOperatorQxx(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Qxx)
{
  uint32_t i,j;
  int    mCr, mAn, lCr, lAn, mk;
  double val;
  lCr = (int)L;
  lAn = (int)L;
  uint16_t oadd[2];
  InitNamedOperator(Qxx, "Qxx");
  Qxx->NFermion = NFermion;
  Qxx->NWord = Min(Qxx->NFermion, 1 + (Qxx->NFermion - 1) / 8);
  InitOperatorFixedLength(Qxx, 2, 8);
  for(i=0;i<2*L+1;i++)
  {
    // C0
    mk  = 0;
    mAn = (int)i - L;
    mCr = mAn+mk;
    val = (double)SlaterCoefficientC(2, lCr, mCr, lAn, mAn);
    if (Abs(val) > DBL_EPSILON)
    {
      j = i + mk;
      oadd[1] = CrAnMask + LindexUp[j];
      oadd[0] = LindexUp[i];
      RealOperatorAddLaderNormalOrder(Qxx, oadd, 2, -val, &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + LindexDn[j];
      oadd[0] = LindexDn[i];
      RealOperatorAddLaderNormalOrder(Qxx, oadd, 2, -val, &QDetNotOrderedFcFaBcBa);
    }
    // C-2
    mk  = -2;
    mAn = (int)i - L;
    mCr = mAn+mk;
    val = (double)SlaterCoefficientC(2, lCr, mCr, lAn, mAn);
    if (Abs(val) > DBL_EPSILON)
    {
      j = i + mk;
      oadd[1] = CrAnMask + LindexUp[j];
      oadd[0] = LindexUp[i];
      RealOperatorAddLaderNormalOrder(Qxx, oadd, 2,  sqrt(1.5)*val, &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + LindexDn[j];
      oadd[0] = LindexDn[i];
      RealOperatorAddLaderNormalOrder(Qxx, oadd, 2,  sqrt(1.5)*val, &QDetNotOrderedFcFaBcBa);
    }
    // C2
    mk  = 2;
    mAn = (int)i - L;
    mCr = mAn+mk;
    val = (double)SlaterCoefficientC(2, lCr, mCr, lAn, mAn);
    if (Abs(val) > DBL_EPSILON)
    {
      j = i + mk;
      oadd[1] = CrAnMask + LindexUp[j];
      oadd[0] = LindexUp[i];
      RealOperatorAddLaderNormalOrder(Qxx, oadd, 2,  sqrt(1.5)*val, &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + LindexDn[j];
      oadd[0] = LindexDn[i];
      RealOperatorAddLaderNormalOrder(Qxx, oadd, 2,  sqrt(1.5)*val, &QDetNotOrderedFcFaBcBa);
    }
  }
  return 0;
}
uint32_t CreateOperatorQxxY(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Qxx)
{
  return CreateOperatorldots(NFermion, L, LindexUp, LindexDn, Qxx);
}
uint32_t CreateOperatorQxxK(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Qxx)
{
  return CreateOperatorInShellBasis(NFermion, L, 2, LindexUp, LindexDn, &CreateOperatorQxx, Qxx);
}
uint32_t CreateOperatorQxxZ(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Qxx)
{
  return CreateOperatorInShellBasis(NFermion, L, 1, LindexUp, LindexDn, &CreateOperatorQxx, Qxx);
}
uint32_t CreateOperatorQxxP(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Qxx)
{
  return CreateOperatorInShellBasis(NFermion, L, 3, LindexUp, LindexDn, &CreateOperatorQxx, Qxx);
}
uint32_t CreateOperatorQxxjjz(uint32_t NFermion, uint32_t L, uint16_t * Lindexjmin, uint16_t * Lindexjplus, OperatorType * const Qxx)
{
  return CreateOperatorInShellBasis(NFermion, L, 4, Lindexjmin, Lindexjplus, &CreateOperatorQxx, Qxx);
}

uint32_t CreateOperatorQyy(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Qyy)
{
  uint32_t i,j;
  int    mCr, mAn, lCr, lAn, mk;
  double val;
  lCr = (int)L;
  lAn = (int)L;
  uint16_t oadd[2];
  InitNamedOperator(Qyy, "Qyy");
  Qyy->NFermion = NFermion;
  Qyy->NWord = Min(Qyy->NFermion, 1 + (Qyy->NFermion - 1) / 8);
  InitOperatorFixedLength(Qyy, 2, 8);
  for(i=0;i<2*L+1;i++)
  {
    // C0
    mk  = 0;
    mAn = (int)i - L;
    mCr = mAn+mk;
    val = (double)SlaterCoefficientC(2, lCr, mCr, lAn, mAn);
    if (Abs(val) > DBL_EPSILON)
    {
      j = i + mk;
      oadd[1] = CrAnMask + LindexUp[j];
      oadd[0] = LindexUp[i];
      RealOperatorAddLaderNormalOrder(Qyy, oadd, 2, -val, &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + LindexDn[j];
      oadd[0] = LindexDn[i];
      RealOperatorAddLaderNormalOrder(Qyy, oadd, 2, -val, &QDetNotOrderedFcFaBcBa);
    }
    // C-2
    mk  = -2;
    mAn = (int)i - L;
    mCr = mAn+mk;
    val = (double)SlaterCoefficientC(2, lCr, mCr, lAn, mAn);
    if (Abs(val) > DBL_EPSILON)
    {
      j = i + mk;
      oadd[1] = CrAnMask + LindexUp[j];
      oadd[0] = LindexUp[i];
      RealOperatorAddLaderNormalOrder(Qyy, oadd, 2,  -sqrt(1.5)*val, &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + LindexDn[j];
      oadd[0] = LindexDn[i];
      RealOperatorAddLaderNormalOrder(Qyy, oadd, 2,  -sqrt(1.5)*val, &QDetNotOrderedFcFaBcBa);
    }
    // C2
    mk  = 2;
    mAn = (int)i - L;
    mCr = mAn+mk;
    val = (double)SlaterCoefficientC(2, lCr, mCr, lAn, mAn);
    if (Abs(val) > DBL_EPSILON)
    {
      j = i + mk;
      oadd[1] = CrAnMask + LindexUp[j];
      oadd[0] = LindexUp[i];
      RealOperatorAddLaderNormalOrder(Qyy, oadd, 2,  -sqrt(1.5)*val, &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + LindexDn[j];
      oadd[0] = LindexDn[i];
      RealOperatorAddLaderNormalOrder(Qyy, oadd, 2,  -sqrt(1.5)*val, &QDetNotOrderedFcFaBcBa);
    }
  }
  return 0;
}
uint32_t CreateOperatorQyyY(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Qyy)
{
  return CreateOperatorldots(NFermion, L, LindexUp, LindexDn, Qyy);
}
uint32_t CreateOperatorQyyK(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Qyy)
{
  return CreateOperatorInShellBasis(NFermion, L, 2, LindexUp, LindexDn, &CreateOperatorQyy, Qyy);
}
uint32_t CreateOperatorQyyZ(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Qyy)
{
  return CreateOperatorInShellBasis(NFermion, L, 1, LindexUp, LindexDn, &CreateOperatorQyy, Qyy);
}
uint32_t CreateOperatorQyyP(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Qyy)
{
  return CreateOperatorInShellBasis(NFermion, L, 3, LindexUp, LindexDn, &CreateOperatorQyy, Qyy);
}
uint32_t CreateOperatorQyyjjz(uint32_t NFermion, uint32_t L, uint16_t * Lindexjmin, uint16_t * Lindexjplus, OperatorType * const Qyy)
{
  return CreateOperatorInShellBasis(NFermion, L, 4, Lindexjmin, Lindexjplus, &CreateOperatorQyy, Qyy);
}

uint32_t CreateOperatorQzz(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Qzz)
{
  uint32_t i,j;
  int    mCr, mAn, lCr, lAn, mk;
  double val;
  lCr = (int)L;
  lAn = (int)L;
  uint16_t oadd[2];
  InitNamedOperator(Qzz, "Qzz");
  Qzz->NFermion = NFermion;
  Qzz->NWord = Min(Qzz->NFermion, 1 + (Qzz->NFermion - 1) / 8);
  InitOperatorFixedLength(Qzz, 2, 8);
  for(i=0;i<2*L+1;i++)
  {
    // C0
    mk  = 0;
    mAn = (int)i - L;
    mCr = mAn+mk;
    val = (double)SlaterCoefficientC(2, lCr, mCr, lAn, mAn);
    if (Abs(val) > DBL_EPSILON)
    {
      j = i + mk;
      oadd[1] = CrAnMask + LindexUp[j];
      oadd[0] = LindexUp[i];
      RealOperatorAddLaderNormalOrder(Qzz, oadd, 2, 2*val, &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + LindexDn[j];
      oadd[0] = LindexDn[i];
      RealOperatorAddLaderNormalOrder(Qzz, oadd, 2, 2*val, &QDetNotOrderedFcFaBcBa);
    }
  }
  return 0;
}
uint32_t CreateOperatorQzzY(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Qzz)
{
  return CreateOperatorQzz(NFermion, L, LindexUp, LindexDn, Qzz);
}
uint32_t CreateOperatorQzzK(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Qzz)
{
  return CreateOperatorInShellBasis(NFermion, L, 2, LindexUp, LindexDn, &CreateOperatorQzz, Qzz);
}
uint32_t CreateOperatorQzzZ(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Qzz)
{
  return CreateOperatorInShellBasis(NFermion, L, 1, LindexUp, LindexDn, &CreateOperatorQzz, Qzz);
}
uint32_t CreateOperatorQzzP(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Qzz)
{
  return CreateOperatorInShellBasis(NFermion, L, 3, LindexUp, LindexDn, &CreateOperatorQzz, Qzz);
}
uint32_t CreateOperatorQzzjjz(uint32_t NFermion, uint32_t L, uint16_t * Lindexjmin, uint16_t * Lindexjplus, OperatorType * const Qzz)
{
  return CreateOperatorInShellBasis(NFermion, L, 4, Lindexjmin, Lindexjplus, &CreateOperatorQzz, Qzz);
}

uint32_t CreateOperatorQxy(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Qxy)
{
  uint32_t i,j;
  int    mCr, mAn, lCr, lAn, mk;
  double val;
  lCr = (int)L;
  lAn = (int)L;
  uint16_t oadd[2];
  InitNamedOperator(Qxy, "Qxy");
  Qxy->NFermion = NFermion;
  Qxy->NWord = Min(Qxy->NFermion, 1 + (Qxy->NFermion - 1) / 8);
  InitOperatorFixedLength(Qxy, 2, 8);
  MakeRealOperatorFixedLengthComplex(Qxy, 2);
  for(i=0;i<2*L+1;i++)
  {
    // C-2
    mk  = -2;
    mAn = (int)i - L;
    mCr = mAn+mk;
    val = (double)SlaterCoefficientC(2, lCr, mCr, lAn, mAn);
    if (Abs(val) > DBL_EPSILON)
    {
      j = i + mk;
      oadd[1] = CrAnMask + LindexUp[j];
      oadd[0] = LindexUp[i];
      ComplexOperatorAddLaderNormalOrder(Qxy, oadd, 2, 0, sqrt(1.5)*val, &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + LindexDn[j];
      oadd[0] = LindexDn[i];
      ComplexOperatorAddLaderNormalOrder(Qxy, oadd, 2, 0, sqrt(1.5)*val, &QDetNotOrderedFcFaBcBa);
    }
    // C2
    mk  = 2;
    mAn = (int)i - L;
    mCr = mAn+mk;
    val = (double)SlaterCoefficientC(2, lCr, mCr, lAn, mAn);
    if (Abs(val) > DBL_EPSILON)
    {
      j = i + mk;
      oadd[1] = CrAnMask + LindexUp[j];
      oadd[0] = LindexUp[i];
      ComplexOperatorAddLaderNormalOrder(Qxy, oadd, 2, 0, -sqrt(1.5)*val, &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + LindexDn[j];
      oadd[0] = LindexDn[i];
      ComplexOperatorAddLaderNormalOrder(Qxy, oadd, 2, 0, -sqrt(1.5)*val, &QDetNotOrderedFcFaBcBa);
    }
  }
  return 0;
}
uint32_t CreateOperatorQxyY(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Qxy)
{
  return CreateOperatorQxy(NFermion, L, LindexUp, LindexDn, Qxy);
}
uint32_t CreateOperatorQxyK(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Qxy)
{
  return CreateOperatorInShellBasis(NFermion, L, 2, LindexUp, LindexDn, &CreateOperatorQxy, Qxy);
}
uint32_t CreateOperatorQxyZ(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Qxy)
{
  return CreateOperatorInShellBasis(NFermion, L, 1, LindexUp, LindexDn, &CreateOperatorQxy, Qxy);
}
uint32_t CreateOperatorQxyP(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Qxy)
{
  return CreateOperatorInShellBasis(NFermion, L, 3, LindexUp, LindexDn, &CreateOperatorQxy, Qxy);
}
uint32_t CreateOperatorQxyjjz(uint32_t NFermion, uint32_t L, uint16_t * Lindexjmin, uint16_t * Lindexjplus, OperatorType * const Qxy)
{
  return CreateOperatorInShellBasis(NFermion, L, 4, Lindexjmin, Lindexjplus, &CreateOperatorQxy, Qxy);
}

uint32_t CreateOperatorQxz(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Qxz)
{
  uint32_t i,j;
  int    mCr, mAn, lCr, lAn, mk;
  double val;
  lCr = (int)L;
  lAn = (int)L;
  uint16_t oadd[2];
  InitNamedOperator(Qxz, "Qxz");
  Qxz->NFermion = NFermion;
  Qxz->NWord = Min(Qxz->NFermion, 1 + (Qxz->NFermion - 1) / 8);
  InitOperatorFixedLength(Qxz, 2, 8);
  for(i=0;i<2*L+1;i++)
  {
    // C-1
    mk  = -1;
    mAn = (int)i - L;
    mCr = mAn+mk;
    val = (double)SlaterCoefficientC(2, lCr, mCr, lAn, mAn);
    if (Abs(val) > DBL_EPSILON)
    {
      j = i + mk;
      oadd[1] = CrAnMask + LindexUp[j];
      oadd[0] = LindexUp[i];
      RealOperatorAddLaderNormalOrder(Qxz, oadd, 2,  sqrt(1.5)*val, &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + LindexDn[j];
      oadd[0] = LindexDn[i];
      RealOperatorAddLaderNormalOrder(Qxz, oadd, 2,  sqrt(1.5)*val, &QDetNotOrderedFcFaBcBa);
    }
    // C1
    mk  = 1;
    mAn = (int)i - L;
    mCr = mAn+mk;
    val = (double)SlaterCoefficientC(2, lCr, mCr, lAn, mAn);
    if (Abs(val) > DBL_EPSILON)
    {
      j = i + mk;
      oadd[1] = CrAnMask + LindexUp[j];
      oadd[0] = LindexUp[i];
      RealOperatorAddLaderNormalOrder(Qxz, oadd, 2, -sqrt(1.5)*val, &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + LindexDn[j];
      oadd[0] = LindexDn[i];
      RealOperatorAddLaderNormalOrder(Qxz, oadd, 2, -sqrt(1.5)*val, &QDetNotOrderedFcFaBcBa);
    }
  }
  return 0;
}
uint32_t CreateOperatorQxzY(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Qxz)
{
  return CreateOperatorQxz(NFermion, L, LindexUp, LindexDn, Qxz);
}
uint32_t CreateOperatorQxzK(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Qxz)
{
  return CreateOperatorInShellBasis(NFermion, L, 2, LindexUp, LindexDn, &CreateOperatorQxz, Qxz);
}
uint32_t CreateOperatorQxzZ(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Qxz)
{
  return CreateOperatorInShellBasis(NFermion, L, 1, LindexUp, LindexDn, &CreateOperatorQxz, Qxz);
}
uint32_t CreateOperatorQxzP(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Qxz)
{
  return CreateOperatorInShellBasis(NFermion, L, 3, LindexUp, LindexDn, &CreateOperatorQxz, Qxz);
}
uint32_t CreateOperatorQxzjjz(uint32_t NFermion, uint32_t L, uint16_t * Lindexjmin, uint16_t * Lindexjplus, OperatorType * const Qxz)
{
  return CreateOperatorInShellBasis(NFermion, L, 4, Lindexjmin, Lindexjplus, &CreateOperatorQxz, Qxz);
}

uint32_t CreateOperatorQyz(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Qyz)
{
  uint32_t i,j;
  int    mCr, mAn, lCr, lAn, mk;
  double val;
  lCr = (int)L;
  lAn = (int)L;
  uint16_t oadd[2];
  InitNamedOperator(Qyz, "Qyz");
  Qyz->NFermion = NFermion;
  Qyz->NWord = Min(Qyz->NFermion, 1 + (Qyz->NFermion - 1) / 8);
  InitOperatorFixedLength(Qyz, 2, 8);
  MakeRealOperatorFixedLengthComplex(Qyz, 2);
  for(i=0;i<2*L+1;i++)
  {
    // C-1
    mk  = -1;
    mAn = (int)i - L;
    mCr = mAn+mk;
    val = (double)SlaterCoefficientC(2, lCr, mCr, lAn, mAn);
    if (Abs(val) > DBL_EPSILON)
    {
      j = i + mk;
      oadd[1] = CrAnMask + LindexUp[j];
      oadd[0] = LindexUp[i];
      ComplexOperatorAddLaderNormalOrder(Qyz, oadd, 2, 0, sqrt(1.5)*val, &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + LindexDn[j];
      oadd[0] = LindexDn[i];
      ComplexOperatorAddLaderNormalOrder(Qyz, oadd, 2, 0, sqrt(1.5)*val, &QDetNotOrderedFcFaBcBa);
    }
    // C1
    mk  = 1;
    mAn = (int)i - L;
    mCr = mAn+mk;
    val = (double)SlaterCoefficientC(2, lCr, mCr, lAn, mAn);
    if (Abs(val) > DBL_EPSILON)
    {
      j = i + mk;
      oadd[1] = CrAnMask + LindexUp[j];
      oadd[0] = LindexUp[i];
      ComplexOperatorAddLaderNormalOrder(Qyz, oadd, 2, 0, sqrt(1.5)*val, &QDetNotOrderedFcFaBcBa);
      oadd[1] = CrAnMask + LindexDn[j];
      oadd[0] = LindexDn[i];
      ComplexOperatorAddLaderNormalOrder(Qyz, oadd, 2, 0, sqrt(1.5)*val, &QDetNotOrderedFcFaBcBa);
    }
  }
  return 0;
}
uint32_t CreateOperatorQyzY(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Qyz)
{
  return CreateOperatorQyz(NFermion, L, LindexUp, LindexDn, Qyz);
}
uint32_t CreateOperatorQyzK(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Qyz)
{
  return CreateOperatorInShellBasis(NFermion, L, 2, LindexUp, LindexDn, &CreateOperatorQyz, Qyz);
}
uint32_t CreateOperatorQyzZ(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Qyz)
{
  return CreateOperatorInShellBasis(NFermion, L, 1, LindexUp, LindexDn, &CreateOperatorQyz, Qyz);
}
uint32_t CreateOperatorQyzP(uint32_t NFermion, uint32_t L, uint16_t * LindexUp, uint16_t * LindexDn, OperatorType * const Qyz)
{
  return CreateOperatorInShellBasis(NFermion, L, 3, LindexUp, LindexDn, &CreateOperatorQyz, Qyz);
}
uint32_t CreateOperatorQyzjjz(uint32_t NFermion, uint32_t L, uint16_t * Lindexjmin, uint16_t * Lindexjplus, OperatorType * const Qyz)
{
  return CreateOperatorInShellBasis(NFermion, L, 4, Lindexjmin, Lindexjplus, &CreateOperatorQyz, Qyz);
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


uint32_t CreateOperatorCreateAnnihilateList(uint32_t const NFermion, uint16_t const * Crindex, uint16_t const * Anindex, double const * VRe, double const * VIm, const uint32_t N, OperatorType * const OppCrAn)
{
  uint32_t i;
  uint16_t oadd[2];
  InitNamedOperator(OppCrAn, "CrAn");
  OppCrAn->NFermion = NFermion;
  OppCrAn->NWord = Min(OppCrAn->NFermion, 1 + (OppCrAn->NFermion - 1) / 8);
  InitOperatorFixedLength(OppCrAn, 2, 8);
  for(i=0;i<N;i++)
  {
    oadd[1] = CrAnMask +Crindex[i];
    oadd[0] = Anindex[i];
    OperatorAddLader(OppCrAn, oadd, 2, VRe[i], VIm[i]);
  }
  return 0;
}

uint32_t CreateOperatorNumber(uint32_t const NFermion, uint16_t const Crindex, uint16_t const Anindex, OperatorType * const OppN)
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

uint32_t CreateOperatorIdentity(uint32_t const NFermion, OperatorType * const OppN)
{
  InitNamedOperator(OppN, "Identity operator");
  OppN->NFermion = NFermion;
  OppN->NWord = Min(OppN->NFermion, 1 + (OppN->NFermion - 1) / 8);
  RealOperatorAddLader(OppN, NULL, 0, 1);
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

int CreateOperatorHarmonicOscillatorOnStandingWaves(double const a, uint32_t const nmax, OperatorType *const op1)
{
  // generates the potential and kinetic energy terms for the problem of particles in a harmonic-oscillator trap on a basis of standing waves.
  // The box for the plane waves is of size a
  // nmax is the total number of standing waves to be used (here NF = 2*nmax)
  // Order of spins is {dn, up, dn, up, dn, up,...}

  uint16_t i, j;
  uint16_t n1, n2;
  // Variable for the indecies of Cr and An operator
  uint16_t oadd[2];
  uint16_t tmp;
  uint32_t n1Sqr;
  // Variable for the matrix elements of the Hamiltonain (here H is a square matrix since it is a one particle part op)
  double val;
  double const aSqr = Sqr(a);
  // Defining name for op1
  InitNamedOperator(op1, "Harmonic Oscillator Kin+Pot");
  op1->NFermion = 2 * nmax;
  // number of words (meaning the number of 8-bit blocks needed for storing bitstrings)
  op1->NWord = Min(op1->NFermion, 1 + (op1->NFermion - 1) / 8);
  InitOperatorFixedLength(op1, 2, (uint32_t)floor(log(nmax * nmax * nmax) / log(2.0)));

  // first adding the diagonal elements
  for (i = 0; i < nmax; i++)
  {
    n1 = i + 1;
    n1Sqr = n1 * n1;
    val = (n1Sqr * PiSqr / (2 * aSqr) + (aSqr / 24) * (1 - 6 / (n1Sqr * PiSqr)));

    oadd[1] = CrAnMask + 2 * i;
    oadd[0] = 2 * i;
    RealOperatorAddLaderNormalOrder(op1, oadd, 2, val, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 2 * i + 1;
    oadd[0] = 2 * i + 1;
    RealOperatorAddLaderNormalOrder(op1, oadd, 2, val, &QDetNotOrderedFcFaBcBa);
  }

  // now adding the off-diagonal elements (going over only half of the elements since H is hermitian)
  for (i = 0; i < nmax; i++)
  {
    n1 = i + 1;
    n1Sqr = n1 * n1;
    for (j = i + 2; j < nmax; j += 2)
    {
      n2 = j + 1;
      val = ((4 * aSqr * n1 * n2) / (Sqr(n1Sqr - Sqr(n2)) * PiSqr));

      oadd[1] = CrAnMask + 2 * i;
      oadd[0] = 2 * j;
      //       printf("In CreateOperatorHarmonicOscillatorOnStandingWaves we have oadd[1] = %i oadd[0] = %i and CrAnMask = %i\n", oadd[1], oadd[0], CrAnMask);
      RealOperatorAddLaderNormalOrder(op1, oadd, 2, val, &QDetNotOrderedFcFaBcBa);
      // and complex conjugate of the above operator
      tmp = CrAnMask + oadd[1];
      oadd[1] = CrAnMask + oadd[0];
      oadd[0] = tmp;
      RealOperatorAddLaderNormalOrder(op1, oadd, 2, val, &QDetNotOrderedFcFaBcBa);

      oadd[1] = CrAnMask + 2 * i + 1;
      oadd[0] = 2 * j + 1;
      RealOperatorAddLaderNormalOrder(op1, oadd, 2, val, &QDetNotOrderedFcFaBcBa);
      // and complex conjugate of the above operator
      tmp = CrAnMask + oadd[1];
      oadd[1] = CrAnMask + oadd[0];
      oadd[0] = tmp;
      RealOperatorAddLaderNormalOrder(op1, oadd, 2, val, &QDetNotOrderedFcFaBcBa);
    }
  }

  return 0;
}

int CreateOperatorHarmonicOscillatorOnStandingWavesPointInteraction(double const a, uint32_t const nmax, OperatorType *const op2)
{
  // generates the point-contact interaction part of the Hamiltonian for
  // the problem of particles in a harmonic-oscillator trap on a basis of standing waves.
  // Here the two particle terms are of the form:
  // Sum_{n1,n2,n3,n4} V_{n1,n2,n3,n4} a^{+}_{n1,Dn} a^{+}_{n2,Up} a_{n3,Up} a_{n4,Dn}
  // where n1+n2+n3+n4 must be an even number.
  // The box for the stadning waves is of size a
  // nmax is the total number of standing waves to be used (here NF = 2*nmax)
  // Order of spins is {dn, up, dn, up, dn, up,...}

  // Defining name for op2
  InitNamedOperator(op2, "Harmonic Oscillator Point-Contact");
  op2->NFermion = 2 * nmax;
  // number of words (meaning the number of 8-bit blocks needed for storing bitstrings)
  op2->NWord = Min(op2->NFermion, 1 + (op2->NFermion - 1) / 8);
  InitOperatorFixedLength(op2, 4, (uint32_t)floor(log(nmax * nmax * nmax) / log(2.0)));

  uint16_t n2, n3, n4;
  int n1;
  // Variable for the indecies of Cr and An operator
  uint16_t oadd[4];
  // Variable for the coefficients of the terms in the Hamiltonain (i.e. here we have a sum of two-particle operators)
  double ainv = 1 / (2 * a);

  // simple examples illustrating the order of Cr and An ops when defining the op
  /*oadd[0] = CrAnMask + 8;
  oadd[1] = CrAnMask + 9;
  oadd[2] = 1;
  oadd[3] = 2;
  RealOperatorAddLader(op2, oadd, 4, 1.);

  oadd[0] = 6;
  oadd[1] = 5;
  oadd[2] = CrAnMask + 4;
  oadd[3] = CrAnMask + 2;
  RealOperatorAddLader(op2, oadd, 4, -2.);

  oadd[3] = CrAnMask + 1;
  oadd[2] = CrAnMask + 2;
  oadd[1] = 3;
  oadd[0] = 4;
  RealOperatorAddLader(op2, oadd, 4, -3.5);*/

  /*
    Operator of Length   4
    QComplex      =          0 (Real==0 or Complex==1)
    N             =          3 (number of operators of length   4)
    A  2 A  1 C  9 C  8 |  1.00000000000000E+00
    C  2 C  4 A  5 A  6 | -2.00000000000000E+00
    C  1 C  2 A  3 A  4 | -3.50000000000000E+00
   */

  for (n2 = 0; n2 < nmax; n2++)
  {
    for (n3 = 0; n3 < nmax; n3++)
    {
      for (n4 = 0; n4 < nmax; n4++)
      {
        oadd[0] = 2 * n4;
        oadd[1] = 2 * n3 + 1;
        oadd[2] = CrAnMask + 2 * n2 + 1;

        // Below the +3 is because for standing waves ni = 1,2,3,..., so that one must consider ni -> ni+1.
        n1 = n2 + n3 + n4 + 3;
        if (n1 <= nmax)
        {
          oadd[3] = CrAnMask + 2 * (n1 - 1);
          RealOperatorAddLaderNormalOrder(op2, oadd, 4, -ainv, &QDetNotOrderedFcFaBcBa);
        }

        n1 = -n2 + n3 + n4 + 1;
        if (n1 > 0 && n1 <= nmax)
        {
          oadd[3] = CrAnMask + 2 * (n1 - 1);
          RealOperatorAddLaderNormalOrder(op2, oadd, 4, ainv, &QDetNotOrderedFcFaBcBa);
        }

        n1 = n2 - n3 + n4 + 1;
        if (n1 > 0 && n1 <= nmax)
        {
          oadd[3] = CrAnMask + 2 * (n1 - 1);
          RealOperatorAddLaderNormalOrder(op2, oadd, 4, ainv, &QDetNotOrderedFcFaBcBa);
        }

        n1 = n2 + n3 - n4 + 1;
        if (n1 > 0 && n1 <= nmax)
        {
          oadd[3] = CrAnMask + 2 * (n1 - 1);
          RealOperatorAddLaderNormalOrder(op2, oadd, 4, ainv, &QDetNotOrderedFcFaBcBa);
        }

        n1 = -n2 - n3 + n4 - 1;
        if (n1 > 0 && n1 <= nmax)
        {
          oadd[3] = CrAnMask + 2 * (n1 - 1);
          RealOperatorAddLaderNormalOrder(op2, oadd, 4, -ainv, &QDetNotOrderedFcFaBcBa);
        }

        n1 = -n2 + n3 - n4 - 1;
        if (n1 > 0 && n1 <= nmax)
        {
          oadd[3] = CrAnMask + 2 * (n1 - 1);
          RealOperatorAddLaderNormalOrder(op2, oadd, 4, -ainv, &QDetNotOrderedFcFaBcBa);
        }

        n1 = n2 - n3 - n4 - 1;
        if (n1 > 0 && n1 <= nmax)
        {
          oadd[3] = CrAnMask + 2 * (n1 - 1);
          RealOperatorAddLaderNormalOrder(op2, oadd, 4, -ainv, &QDetNotOrderedFcFaBcBa);
        }
      }
    }
  }

  return 0;
}

int CreateOperatorHarmonicOscillatorOnStandingWavesOnBasis(double const a, uint32_t const nmax, CompactMatrixType const *const vec, OperatorType *const op1)
{
  // generates the potential and kinetic energy terms for the problem of particles in a harmonic-oscillator trap on a basis of orbitals given by the rows in the matrix vec in terms of standing waves.
  // The box for the plane waves is of size a
  // nmax is the total number of standing waves to be used (here NF = 2*nmax)
  // Note that the basis of the orbitals must have the same a and nmax, otherwise ends in problematic results and there is no way to check this automatically.
  // Order of spins is {dn, up, dn, up, dn, up,...}

  uint16_t i, j;
  uint16_t n1, n2;
  uint16_t a1, a2;
  // Variable for the indecies of Cr and An operator
  uint16_t oadd[2];
  uint16_t tmp, a1T2, a2T2;
  uint32_t n1Sqr;
  // Variable for the matrix elements of the Hamiltonain (here H is a square matrix since it is a one particle part op)
  double val;
  double const aSqr = Sqr(a);

  // The number of orbitals (without spin)
  uint16_t const NOrbs = (uint16_t)vec->Ni / 2;

  // Defining name for op1
  InitNamedOperator(op1, "Harmonic Oscillator Kin+Pot");
  op1->NFermion = vec->Ni;
  // number of words (meaning the number of 8-bit blocks needed for storing bitstrings)
  op1->NWord = Min(op1->NFermion, 1 + (op1->NFermion - 1) / 8);
  InitOperatorFixedLength(op1, 2, (uint32_t)floor(log(nmax * nmax * nmax) / log(2.0)));

  // Creating the matrix of orbs without spins
  CompactMatrixType *vecWithoutSpin /*= new CompactMatrixType()*/;
  vecWithoutSpin = (CompactMatrixType *)malloc(sizeof(CompactMatrixType));
  vecWithoutSpin->Ni = (uint32_t)NOrbs;
  vecWithoutSpin->Nj = nmax;
  vecWithoutSpin->QComplex = 0;
  //   vecWithoutSpin->VRe = (double*) malloc(sizeof(double)*(vecWithoutSpin->Ni)*(vecWithoutSpin->Nj));
  InitCompactMatrix(vecWithoutSpin);
  for (i = 0; i < NOrbs; i++)
  {
    for (j = 0; j < nmax; j++)
    {
      //         vecWithoutSpin->VRe[i][j] = vec->VRe[2*i-1][2*j-1];
      vecWithoutSpin->VRe[i * vecWithoutSpin->Nj + j] = vec->VRe[(2 * i) * vec->Nj + 2 * j];
    }
  }

  // first adding the diagonal elements H1_{a1,a2} (with a1=a2)
  for (a1 = 0; a1 < NOrbs; a1++)
  {
    val = 0;

    // Performing the sum to find the rotated (diagonal) matrix elements
    // H1_{a1,a2} = Sum_{n1,n2} alpha_{a1,n1} alpha_{a2,n2} H1_{n1,n2} (here a1=a2)

    // Now first only the terms with n1=n2
    for (i = 0; i < nmax; i++)
    {

      n1 = i + 1;
      n1Sqr = n1 * n1;

      val = val + vecWithoutSpin->VRe[a1 * vecWithoutSpin->Nj + i] * vecWithoutSpin->VRe[a1 * vecWithoutSpin->Nj + i] * (n1Sqr * PiSqr / (2 * aSqr) + (aSqr / 24) * (1 - 6 / (n1Sqr * PiSqr)));
    }

    // and now terms with n1!=n2 (splitted into two for-loops in order to avoid terms with n1=n2)
    for (i = 0; i < nmax; i++)
    {
      n1 = i + 1;
      n1Sqr = n1 * n1;
      for (j = i + 2; j < nmax; j += 2)
      {
        n2 = j + 1;
        val = val + vecWithoutSpin->VRe[a1 * vecWithoutSpin->Nj + i] * vecWithoutSpin->VRe[a1 * vecWithoutSpin->Nj + j] * ((4 * aSqr * n1 * n2) / (Sqr(n1Sqr - Sqr(n2)) * PiSqr));
      }
    }

    for (i = 0; i < nmax; i++)
    {
      n1 = i + 1;
      n1Sqr = n1 * n1;
      for (j = i % 2; j < i; j += 2)
      {
        n2 = j + 1;
        val = val + vecWithoutSpin->VRe[a1 * vecWithoutSpin->Nj + i] * vecWithoutSpin->VRe[a1 * vecWithoutSpin->Nj + j] * ((4 * aSqr * n1 * n2) / (Sqr(n1Sqr - Sqr(n2)) * PiSqr));
      }
    }
    // Now adding the terms to the operator
    a1T2 = 2 * a1;

    // spin dn
    oadd[1] = CrAnMask + a1T2;
    oadd[0] = a1T2;
    RealOperatorAddLaderNormalOrder(op1, oadd, 2, val, &QDetNotOrderedFcFaBcBa);
    // spin up
    oadd[1] = CrAnMask + a1T2 + 1;
    oadd[0] = a1T2 + 1;
    RealOperatorAddLaderNormalOrder(op1, oadd, 2, val, &QDetNotOrderedFcFaBcBa);
  }

  // now adding the off-diagonal elements (going over only half of the elements since H is hermitian)
  for (a1 = 0; a1 < NOrbs; a1++)
  {
    for (a2 = a1 + 2; a2 < NOrbs; a2 += 2)
    {
      val = 0;

      // Performing the sum to find the rotated (diagonal) matrix elements
      // H1_{a1,a2} = Sum_{n1,n2} alpha_{a1,n1} alpha_{a2,n2} H1_{n1,n2} (here a1!=a2)

      // Now first only the terms with n1=n2
      for (i = 0; i < nmax; i++)
      {

        n1 = i + 1;
        n1Sqr = n1 * n1;

        val = val + vecWithoutSpin->VRe[a1 * vecWithoutSpin->Nj + i] * vecWithoutSpin->VRe[a2 * vecWithoutSpin->Nj + i] * (n1Sqr * PiSqr / (2 * aSqr) + (aSqr / 24) * (1 - 6 / (n1Sqr * PiSqr)));
      }

      // and now terms with n1!=n2 (splitted into two for-loops)
      for (i = 0; i < nmax; i++)
      {
        n1 = i + 1;
        n1Sqr = n1 * n1;
        for (j = i + 2; j < nmax; j += 2)
        {
          n2 = j + 1;
          val = val + vecWithoutSpin->VRe[a1 * vecWithoutSpin->Nj + i] * vecWithoutSpin->VRe[a2 * vecWithoutSpin->Nj + j] * ((4 * aSqr * n1 * n2) / (Sqr(n1Sqr - Sqr(n2)) * PiSqr));
        }
      }

      for (i = 0; i < nmax; i++)
      {
        n1 = i + 1;
        n1Sqr = n1 * n1;
        for (j = i % 2; j < i; j += 2)
        {
          n2 = j + 1;
          val = val + vecWithoutSpin->VRe[a1 * vecWithoutSpin->Nj + i] * vecWithoutSpin->VRe[a2 * vecWithoutSpin->Nj + j] * ((4 * aSqr * n1 * n2) / (Sqr(n1Sqr - Sqr(n2)) * PiSqr));
        }
      }

      // Now adding the terms to the operator
      a1T2 = 2 * a1;
      a2T2 = 2 * a2;

      // spin dn
      oadd[1] = CrAnMask + a1T2;
      oadd[0] = a2T2;
      RealOperatorAddLaderNormalOrder(op1, oadd, 2, val, &QDetNotOrderedFcFaBcBa);
      // and complex conjugate of the above operator
      tmp = CrAnMask + oadd[1];
      oadd[1] = CrAnMask + oadd[0];
      oadd[0] = tmp;
      RealOperatorAddLaderNormalOrder(op1, oadd, 2, val, &QDetNotOrderedFcFaBcBa);

      // spin up
      oadd[1] = CrAnMask + a1T2 + 1;
      oadd[0] = a2T2 + 1;
      RealOperatorAddLaderNormalOrder(op1, oadd, 2, val, &QDetNotOrderedFcFaBcBa);
      // and complex conjugate of the above operator
      tmp = CrAnMask + oadd[1];
      oadd[1] = CrAnMask + oadd[0];
      oadd[0] = tmp;
      RealOperatorAddLaderNormalOrder(op1, oadd, 2, val, &QDetNotOrderedFcFaBcBa);
    }
  }

  FreeCompactMatrixIncludingPointer(vecWithoutSpin);

  return 0;
}

int CreateOperatorHarmonicOscillatorOnStandingWavesPointInteractionOnBasis(double const a, uint32_t const nmax, CompactMatrixType const *const vec, OperatorType *const op2)
{
  // generates the point-contact interaction part of the Hamiltonian for
  // the problem of particles in a harmonic-oscillator trap on a basis of orbitals given by the rows in the matrix vec in terms of standing waves.
  // Note that the basis of the orbitals must have the same a and nmax, otherwise ends in problematic results and there is no way to check this automatically.
  // Here the two particle terms are of the form:
  // Sum_{n1,n2,n3,n4} V_{n1,n2,n3,n4} a^{+}_{n1,Dn} a^{+}_{n2,Up} a_{n3,Up} a_{n4,Dn}
  // where n1+n2+n3+n4 must be an even number.
  // The box for the stadning waves is of size a
  // nmax is the total number of standing waves to be used (here NF = 2*nmax)
  // Order of spins is {dn, up, dn, up, dn, up,...}

  // The number of orbitals (without spin)
  uint16_t const NOrbs = (uint16_t)vec->Ni / 2;

  // Defining name for op2
  InitNamedOperator(op2, "Harmonic Oscillator Point-Contact");
  op2->NFermion = vec->Ni;
  // number of words (meaning the number of 8-bit blocks needed for storing bitstrings)
  op2->NWord = Min(op2->NFermion, 1 + (op2->NFermion - 1) / 8);
  InitOperatorFixedLength(op2, 4, (uint32_t)floor(log(nmax * nmax * nmax) / log(2.0)));
  int n1, n2, n3, n4;
  uint16_t a1, a2, a3, a4;
  // Variable for the indecies of Cr and An operator
  uint16_t oadd[4];
  // Variable for the coefficients of the terms in the Hamiltonain (i.e. here we have a sum of two-particle operators)
  double ainv = 1 / (2 * a);
  double val;

  // Creating the matrix of orbs without spins
  CompactMatrixType *vecWithoutSpin2;
  vecWithoutSpin2 = (CompactMatrixType *)malloc(sizeof(CompactMatrixType));
  vecWithoutSpin2->Ni = (uint32_t)NOrbs;
  vecWithoutSpin2->Nj = nmax;
  vecWithoutSpin2->QComplex = 0;
  InitCompactMatrix(vecWithoutSpin2);
  for (n1 = 0; n1 < NOrbs; n1++)
  {
    for (n2 = 0; n2 < nmax; n2++)
    {
      vecWithoutSpin2->VRe[n1 * (vecWithoutSpin2->Nj) + n2] = vec->VRe[(2 * n1) * (vec->Nj) + 2 * n2];
    }
  }

  // Calculating all matrix elements of the Hamiltonian on the given orbs with indecies a1, a2, a3, and a4 (i.e. H1_{a1,a2,a3,a4})
  // n1, n2, n3, and n4 play the role of dummy indecies in the sum
  // H1_{a1,a2,a3,a4} = Sum_{n1,n2,n3,n4} alpha_{a1,n1} alpha_{a2,n2} alpha_{a3,n3} alpha_{a4,n4} H1_{n1,n2,n3,n4}
  // where alpha are the elments of the matrix vecWithoutSpin2.
  for (a1 = 0; a1 < NOrbs; a1++)
  {
    for (a2 = 0; a2 < NOrbs; a2++)
    {
      for (a3 = 0; a3 < NOrbs; a3++)
      {
        for (a4 = (a1 + a2 + a3) % 2; a4 < NOrbs; a4++)
        {
          val = 0;
#pragma omp parallel for private(n1,n3,n4) reduction(+ : val)
          for (n2 = 0; n2 < nmax; n2++)
          {
            for (n3 = 0; n3 < nmax; n3++)
            {
              for (n4 = 0; n4 < nmax; n4++)
              {
                n1 = n2 + n3 + n4 + 3;
                if (n1 <= nmax)
                {
                  val = val - (vecWithoutSpin2->VRe[a1 * (vecWithoutSpin2->Nj) + ((uint16_t)n1) - 1]) * (vecWithoutSpin2->VRe[a2 * (vecWithoutSpin2->Nj) + n2]) * (vecWithoutSpin2->VRe[a3 * (vecWithoutSpin2->Nj) + n3]) * (vecWithoutSpin2->VRe[a4 * (vecWithoutSpin2->Nj) + n4]) * ainv;
                }

                n1 = -n2 + n3 + n4 + 1;
                if ((n1 > 0) && (n1 <= nmax))
                {
                  val = val + (vecWithoutSpin2->VRe[a1 * (vecWithoutSpin2->Nj) + ((uint16_t)n1) - 1]) * (vecWithoutSpin2->VRe[a2 * (vecWithoutSpin2->Nj) + n2]) * (vecWithoutSpin2->VRe[a3 * (vecWithoutSpin2->Nj) + n3]) * (vecWithoutSpin2->VRe[a4 * (vecWithoutSpin2->Nj) + n4]) * ainv;
                }

                n1 = n2 - n3 + n4 + 1;
                if ((n1 > 0) && (n1 <= nmax))
                {
                  val = val + (vecWithoutSpin2->VRe[a1 * (vecWithoutSpin2->Nj) + ((uint16_t)n1) - 1]) * (vecWithoutSpin2->VRe[a2 * (vecWithoutSpin2->Nj) + n2]) * (vecWithoutSpin2->VRe[a3 * (vecWithoutSpin2->Nj) + n3]) * (vecWithoutSpin2->VRe[a4 * (vecWithoutSpin2->Nj) + n4]) * ainv;
                }

                n1 = n2 + n3 - n4 + 1;
                if ((n1 > 0) && (n1 <= nmax))
                {
                  val = val + (vecWithoutSpin2->VRe[a1 * (vecWithoutSpin2->Nj) + ((uint16_t)n1) - 1]) * (vecWithoutSpin2->VRe[a2 * (vecWithoutSpin2->Nj) + n2]) * (vecWithoutSpin2->VRe[a3 * (vecWithoutSpin2->Nj) + n3]) * (vecWithoutSpin2->VRe[a4 * (vecWithoutSpin2->Nj) + n4]) * ainv;
                }

                n1 = -n2 - n3 + n4 - 1;
                if ((n1 > 0) && (n1 <= nmax))
                {
                  val = val - (vecWithoutSpin2->VRe[a1 * (vecWithoutSpin2->Nj) + ((uint16_t)n1) - 1]) * (vecWithoutSpin2->VRe[a2 * (vecWithoutSpin2->Nj) + n2]) * (vecWithoutSpin2->VRe[a3 * (vecWithoutSpin2->Nj) + n3]) * (vecWithoutSpin2->VRe[a4 * (vecWithoutSpin2->Nj) + n4]) * ainv;
                }

                n1 = -n2 + n3 - n4 - 1;
                if ((n1 > 0) && (n1 <= nmax))
                {
                  val = val - (vecWithoutSpin2->VRe[a1 * (vecWithoutSpin2->Nj) + ((uint16_t)n1) - 1]) * (vecWithoutSpin2->VRe[a2 * (vecWithoutSpin2->Nj) + n2]) * (vecWithoutSpin2->VRe[a3 * (vecWithoutSpin2->Nj) + n3]) * (vecWithoutSpin2->VRe[a4 * (vecWithoutSpin2->Nj) + n4]) * ainv;
                }

                n1 = n2 - n3 - n4 - 1;
                if ((n1 > 0) && (n1 <= nmax))
                {
                  val = val - (vecWithoutSpin2->VRe[a1 * (vecWithoutSpin2->Nj) + ((uint16_t)n1) - 1]) * (vecWithoutSpin2->VRe[a2 * (vecWithoutSpin2->Nj) + n2]) * (vecWithoutSpin2->VRe[a3 * (vecWithoutSpin2->Nj) + n3]) * (vecWithoutSpin2->VRe[a4 * (vecWithoutSpin2->Nj) + n4]) * ainv;
                }
              }
            }
          }
          oadd[3] = CrAnMask + 2 * a1;
          oadd[2] = CrAnMask + 2 * a2 + 1;
          oadd[1] = 2 * a3 + 1;
          oadd[0] = 2 * a4;
          RealOperatorAddLaderNormalOrder(op2, oadd, 4, val, &QDetNotOrderedFcFaBcBa);
        }
      }
    }
  }

  FreeCompactMatrixIncludingPointer(vecWithoutSpin2);

  return 0;
}





int CreateOperatorHarmonicOscillatorOnGrid(double const a, int const nmax, OperatorType *const op1)
{
  // the grid points are -d to d in steps of d/nmax
  // there are  2nmax+1 grid points
  // the function on the first and the last grid point is zero we  only return 2nmax-1 points
  uint16_t oadd[2];
  int i;
  double const dx = a / (double)nmax;
  double const t = 1 / Sqr(dx);
  double val;
  
  InitNamedOperator(op1, "Harmonic Oscillator Kin+Pot");
  op1->NFermion = 2 * (2*nmax-1);
  // number of words (meaning the number of 8-bit blocks needed for storing bitstrings)
  op1->NWord = Min(op1->NFermion, 1 + (op1->NFermion - 1) / 8);
  InitOperatorFixedLength(op1, 2, (uint32_t)floor(log(nmax * nmax * nmax) / log(2.0)));
  
  // first adding the diagonal elements in position
  for (i = 0; i < 2*nmax-1; i++)
  {
    val = t + 0.5 * Sqr(dx * (1+i-nmax)) ;
    oadd[1] = CrAnMask + 2 * i;
    oadd[0] = 2 * i;
    RealOperatorAddLaderNormalOrder(op1, oadd, 2, val, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 2 * i + 1;
    oadd[0] = 2 * i + 1;
    RealOperatorAddLaderNormalOrder(op1, oadd, 2, val, &QDetNotOrderedFcFaBcBa);
  }
  
  // now adding the elements between site i and i+1
  val = -t/2;
  for (i = 0; i < 2*nmax-2; i++)
  {
    oadd[1] = CrAnMask + 2 * i;
    oadd[0] = 2 * (i+1);
    RealOperatorAddLaderNormalOrder(op1, oadd, 2, val, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 2 * i + 1;
    oadd[0] = 2 * (i+1) + 1;
    RealOperatorAddLaderNormalOrder(op1, oadd, 2, val, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 2 * (i+1);
    oadd[0] = 2 * i;
    RealOperatorAddLaderNormalOrder(op1, oadd, 2, val, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 2 * (i+1) + 1;
    oadd[0] = 2 * i + 1;
    RealOperatorAddLaderNormalOrder(op1, oadd, 2, val, &QDetNotOrderedFcFaBcBa);
  }
  return 0;
}

int CreateOperatorHarmonicOscillatorOnGridInteraction(double const a, uint32_t const nmax, OperatorType *const op2)
{
  // the grid points are -d to d in steps of d/nmax
  // there are 2nmax+1 grid points
  // the function on the first and the last grid point is zero we  only return 2nmax-1 points
  uint16_t i,oadd[4];
  double const dx = a / nmax;
  double const g = 1 / dx;
  
  InitNamedOperator(op2, "Harmonic Oscillator Point Contact Interaction");
  op2->NFermion = 2 * (2*nmax-1);
  // number of words (meaning the number of 8-bit blocks needed for storing bitstrings)
  op2->NWord = Min(op2->NFermion, 1 + (op2->NFermion - 1) / 8);
  InitOperatorFixedLength(op2, 4, (uint32_t)floor(log(nmax * nmax * nmax) / log(2.0)));
  // adding the diagonal elements in position
  for (i = 0; i < 2*nmax-1; i++)
  {
    oadd[3] = CrAnMask + 2 * i;
    oadd[2] = 2 * i;
    oadd[1] = CrAnMask + 2 * i + 1;
    oadd[0] = 2 * i + 1;
    RealOperatorAddLaderNormalOrder(op2, oadd, 4, g, &QDetNotOrderedFcFaBcBa);
  }
  return 0;
}


int CreateOperatorHarmonicOscillatorOnGridOnBasis(double const a, int const nmax, CompactMatrixType const *const vec, OperatorType *const op1)
{
  // generates the potential and kinetic energy terms for the problem of particles in a harmonic-oscillator trap on a basis of orbitals given by the rows in the matrix vec in terms of grid points.
  // The box for the plane waves is of size -a to a
  // 2*nmax + 1 is the total number of grid-points
  // 2*nmax - 1 is the total number of grid-points returned as the functions are zero on the boundary
  // Note that the basis of the orbitals must have the same a and nmax, otherwise ends in problematic results and there is no way to check this automatically.
  // Order of spins is {dn, up, dn, up, dn, up,...}
  
  // The basis vec is assumed to be without spin

  uint16_t i, j;
  int n;
  // Variable for the indecies of Cr and An operator
  uint16_t oadd[2];
  // Variable for the matrix elements of the Hamiltonain (here H is a square matrix since it is a one particle part op)
  double val;
  double const dx = a/nmax;
  double const t  = 1/Sqr(dx);
  uint16_t const NOrbs = (uint16_t)vec->Ni;
  
  // Defining name for op1
  InitNamedOperator(op1, "Harmonic Oscillator Kin+Pot");
  op1->NFermion = 2*NOrbs;
  // number of words (meaning the number of 8-bit blocks needed for storing bitstrings)
  op1->NWord = Min(op1->NFermion, 1 + (op1->NFermion - 1) / 8);
  InitOperatorFixedLength(op1, 2, (uint32_t)floor(log(nmax * nmax * nmax) / log(2.0)));
  if(vec->QComplex)
  {
    printf("Sorry complex basis rotations not yet implemented in CreateOperatorHarmonicOscillatorOnGridOnBasis");
    return 1;
  }
  for (i=0; i<NOrbs; i++)
  {
    for (j=0; j<NOrbs; j++)
    {
      val = 0;
      for (n=0; n<2*nmax-1; n++)
      {
        val += (t + 0.5*Sqr(dx * (1+n-nmax))) * vec->VRe[i * vec->Nj + n] * vec->VRe[j * vec->Nj + n];
      }
      for (n=0; n<2*nmax-1; n++)
      {
        val += (-0.5*t) * ( vec->VRe[i * vec->Nj + n] * vec->VRe[j * vec->Nj + n + 1] + vec->VRe[i * vec->Nj + n + 1] * vec->VRe[j * vec->Nj + n] );
      }
      // spin dn
      oadd[1] = CrAnMask + 2*i;
      oadd[0] = 2*j;
      RealOperatorAddLaderNormalOrder(op1, oadd, 2, val, &QDetNotOrderedFcFaBcBa);
      // spin up
      oadd[1] = CrAnMask + 2*i + 1;
      oadd[0] = 2*j + 1;
      RealOperatorAddLaderNormalOrder(op1, oadd, 2, val, &QDetNotOrderedFcFaBcBa);
    }
  }
  return 0;
}



int CreateOperatorHarmonicOscillatorOnGridInteractionOnBasis(double const a, uint32_t const nmax, CompactMatrixType const *const vec, OperatorType *const op2)
{
  // the grid points are -a to a in steps of a/nmax
  // there are thus 2nmax+1 grid points
  // 2*nmax - 1 is the total number of grid-points returned as the functions are zero on the boundary
  uint16_t i,j,k,l,oadd[4];
  int n;
  double const dx = a / nmax;
  double const g = 1 / dx;
  uint16_t const NOrbs = (uint16_t)vec->Ni;
  double val;
  
  InitNamedOperator(op2, "Harmonic Oscillator Point Contact Interaction");
  op2->NFermion = 2 * NOrbs;
  // number of words (meaning the number of 8-bit blocks needed for storing bitstrings)
  op2->NWord = Min(op2->NFermion, 1 + (op2->NFermion - 1) / 8);
  InitOperatorFixedLength(op2, 4, (uint32_t)floor(log(nmax * nmax * nmax) / log(2.0)));

  if(vec->QComplex)
  {
    printf("Sorry complex basis rotations not yet implemented in CreateOperatorHarmonicOscillatorOnGridOnBasis");
    return 1;
  }

  for (i=0; i<NOrbs; i++)
  {
    for (j=0; j<NOrbs; j++)
    {
      for (k=0; k<NOrbs; k++)
      {
        for (l=0; l<NOrbs; l++)
        {
          val = 0;
          for (n=0; n<2*nmax-1; n++)
          {
            val += g * vec->VRe[i * vec->Nj + n] * vec->VRe[j * vec->Nj + n] * vec->VRe[k * vec->Nj + n] * vec->VRe[l * vec->Nj + n];
          }
          oadd[3] = CrAnMask + 2*l;
          oadd[2] = 2*k;
          oadd[1] = CrAnMask + 2*j+1;
          oadd[0] = 2*i+1;
          RealOperatorAddLaderNormalOrder(op2, oadd, 4, val, &QDetNotOrderedFcFaBcBa);
        }
      }
    }
  }
  return 0;
}




int CreateOverlapMatrixHarmonicOscillatorOnSplineWaveletA(uint32_t const nmax, CompactMatrixType *S)
{
  // the basis functions are
  // g1[x] = (Sqrt[35/13]*Piecewise[{{1 - (3*x^2)/4 - x^3/4, x > -2 && x <= 0}, {1 - (3*x^2)/4 + x^3/4, x > 0 && x < 2}}, 0])/2
  // g2[x] = (Sqrt[105]*Piecewise[{{x*(2 + x)^2, x > -2 && x <= 0}, {(2 - x)^2*x, x > 0 && x < 2}}, 0])/16
  // f1[x,n,d] = sqrt[2/d] g1[2x/d-2n]
  // f2[x,n,d] = sqrt[2/d] g2[2x/d-2n]
  // The basis is ordered from n = -nmax to n = nmax
  // For each n there are 4 states 1down, 1up, 2down, 2up
  int i;
  S->QComplex = 0;
  S->Ni = 4 * (2 * nmax + 1);
  S->Nj = 4 * (2 * nmax + 1);
  InitCompactMatrix(S);
  for (i = 0; i < 2*nmax+1; i++)
  {
    S->VRe[4*i+0 + S->Nj * (4*i+0)] = 1;
    S->VRe[4*i+1 + S->Nj * (4*i+1)] = 1;
    S->VRe[4*i+2 + S->Nj * (4*i+2)] = 1;
    S->VRe[4*i+3 + S->Nj * (4*i+3)] = 1;
  }
  
  for (i = 0; i < 2*nmax; i++)
  {
    S->VRe[4*(i+1)+0 + S->Nj * (4*(i+0)+0)] = 9./52.;
    S->VRe[4*(i+0)+0 + S->Nj * (4*(i+1)+0)] = 9./52.;
    S->VRe[4*(i+1)+1 + S->Nj * (4*(i+0)+1)] = 9./52.;
    S->VRe[4*(i+0)+1 + S->Nj * (4*(i+1)+1)] = 9./52.;
    
    S->VRe[4*(i+1)+2 + S->Nj * (4*(i+0)+2)] =-3./8.;
    S->VRe[4*(i+0)+2 + S->Nj * (4*(i+1)+2)] =-3./8.;
    S->VRe[4*(i+1)+3 + S->Nj * (4*(i+0)+3)] =-3./8.;
    S->VRe[4*(i+0)+3 + S->Nj * (4*(i+1)+3)] =-3./8.;
    
    S->VRe[4*(i+1)+0 + S->Nj * (4*(i+0)+2)] = sqrt(13./3.)/8.;
    S->VRe[4*(i+0)+2 + S->Nj * (4*(i+1)+0)] = sqrt(13./3.)/8.;
    S->VRe[4*(i+1)+2 + S->Nj * (4*(i+0)+0)] =-sqrt(13./3.)/8.;
    S->VRe[4*(i+0)+0 + S->Nj * (4*(i+1)+2)] =-sqrt(13./3.)/8.;
    S->VRe[4*(i+1)+1 + S->Nj * (4*(i+0)+3)] = sqrt(13./3.)/8.;
    S->VRe[4*(i+0)+3 + S->Nj * (4*(i+1)+1)] = sqrt(13./3.)/8.;
    S->VRe[4*(i+1)+3 + S->Nj * (4*(i+0)+1)] =-sqrt(13./3.)/8.;
    S->VRe[4*(i+0)+1 + S->Nj * (4*(i+1)+3)] =-sqrt(13./3.)/8.;
    
  }
  return 0;
}

int CreateOperatorHarmonicOscillatorOnSplineWaveletA(double const d, uint32_t const nmax, OperatorType *const op1)
{
  // the basis functions are
  // g1[x] = (Sqrt[35/13]*Piecewise[{{1 - (3*x^2)/4 - x^3/4, x > -2 && x <= 0}, {1 - (3*x^2)/4 + x^3/4, x > 0 && x < 2}}, 0])/2
  // g2[x] = (Sqrt[105]*Piecewise[{{x*(2 + x)^2, x > -2 && x <= 0}, {(2 - x)^2*x, x > 0 && x < 2}}, 0])/16
  // f1[x,n,d] = sqrt[2/d] g1[2x/d-2n]
  // f2[x,n,d] = sqrt[2/d] g2[2x/d-2n]
  uint16_t i,oadd[2];
  double const Ekin110 = 21/(d*d*13);
  double const Ekin11p =-21/(d*d*26);
  double const Ekin220 = 7/(d*d);
  double const Ekin22p =-7/(d*d*8);
  double const Ekin12p = 7*sqrt(3./13.)/(d*d*8);
  double const Ekin12m =-7*sqrt(3./13.)/(d*d*8);
  double val;
  
  InitNamedOperator(op1, "Harmonic Oscillator Kin+Pot");
  op1->NFermion = 4 * (2*nmax+1);
  // number of words (meaning the number of 8-bit blocks needed for storing bitstrings)
  op1->NWord = Min(op1->NFermion, 1 + (op1->NFermion - 1) / 8);
  InitOperatorFixedLength(op1, 2, (uint32_t)floor(log(nmax * nmax * nmax) / log(2.0)));
  
  // first adding the diagonal elements in position
  for (i = 0; i < 2*nmax+1; i++)
  {
    val = Ekin110 + Sqr(d) * (19/468. + Sqr(i-nmax)/2. );
    oadd[1] = CrAnMask + 4 * i;
    oadd[0] = 4 * i;
    RealOperatorAddLaderNormalOrder(op1, oadd, 2, val, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4 * i + 1;
    oadd[0] = 4 * i + 1;
    RealOperatorAddLaderNormalOrder(op1, oadd, 2, val, &QDetNotOrderedFcFaBcBa);
    
    val = Ekin220 + Sqr(d) * (1/12.) * (1 + 6 * Sqr(i-nmax));
    oadd[1] = CrAnMask + 4 * i + 2;
    oadd[0] = 4 * i + 2;
    RealOperatorAddLaderNormalOrder(op1, oadd, 2, val, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4 * i + 3;
    oadd[0] = 4 * i + 3;
    RealOperatorAddLaderNormalOrder(op1, oadd, 2, val, &QDetNotOrderedFcFaBcBa);
    
    val = Sqr(d) * 7 * (i-nmax) / (4 * sqrt(39.));
    oadd[1] = CrAnMask + 4 * i + 0;
    oadd[0] = 4 * i + 2;
    RealOperatorAddLaderNormalOrder(op1, oadd, 2, val, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4 * i + 1;
    oadd[0] = 4 * i + 3;
    RealOperatorAddLaderNormalOrder(op1, oadd, 2, val, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4 * i + 2;
    oadd[0] = 4 * i + 0;
    RealOperatorAddLaderNormalOrder(op1, oadd, 2, val, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4 * i + 3;
    oadd[0] = 4 * i + 1;
    RealOperatorAddLaderNormalOrder(op1, oadd, 2, val, &QDetNotOrderedFcFaBcBa);
    
  }
  
  // now adding the elements between site i and i+1
  for (i = 0; i < 2*nmax; i++)
  {
    val = Ekin11p + Sqr(d) * (23 + 81 * (i-nmax) * (1 + i - nmax))/936.;
    oadd[1] = CrAnMask + 4 * i;
    oadd[0] = 4 * (i+1);
    RealOperatorAddLaderNormalOrder(op1, oadd, 2, val, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4 * i + 1;
    oadd[0] = 4 * (i+1) + 1;
    RealOperatorAddLaderNormalOrder(op1, oadd, 2, val, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4 * (i+1);
    oadd[0] = 4 * i;
    RealOperatorAddLaderNormalOrder(op1, oadd, 2, val, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4 * (i+1) + 1;
    oadd[0] = 4 * i + 1;
    RealOperatorAddLaderNormalOrder(op1, oadd, 2, val, &QDetNotOrderedFcFaBcBa);
    
    
    val = Ekin22p + Sqr(d) * (-5 - 18 * (i-nmax) * (1 + i - nmax))/96.;
    oadd[1] = CrAnMask + 4 * i;
    oadd[0] = 4 * (i+1);
    RealOperatorAddLaderNormalOrder(op1, oadd, 2, val, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4 * i + 1;
    oadd[0] = 4 * (i+1) + 1;
    RealOperatorAddLaderNormalOrder(op1, oadd, 2, val, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4 * (i+1);
    oadd[0] = 4 * i;
    RealOperatorAddLaderNormalOrder(op1, oadd, 2, val, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4 * (i+1) + 1;
    oadd[0] = 4 * i + 1;
    RealOperatorAddLaderNormalOrder(op1, oadd, 2, val, &QDetNotOrderedFcFaBcBa);
    
    val = Ekin12p - Sqr(d) * (19 + 72 * (i-nmax) + 78 * Sqr(i-nmax) ) / (4 * sqrt(39));
    oadd[1] = CrAnMask + 4 * i;
    oadd[0] = 4 * (i+1) + 2;
    RealOperatorAddLaderNormalOrder(op1, oadd, 2, val, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4 * i + 1;
    oadd[0] = 4 * (i+1) + 3;
    RealOperatorAddLaderNormalOrder(op1, oadd, 2, val, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4 * (i+1) + 2;
    oadd[0] = 4 * i;
    RealOperatorAddLaderNormalOrder(op1, oadd, 2, val, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4 * (i+1) + 3;
    oadd[0] = 4 * i + 1;
    RealOperatorAddLaderNormalOrder(op1, oadd, 2, val, &QDetNotOrderedFcFaBcBa);
    
    val = Ekin12m + Sqr(d) * (25 + 84 * (i-nmax) + 78 * Sqr(i-nmax) ) / (96 * sqrt(39));
    oadd[1] = CrAnMask + 4 * i + 2;
    oadd[0] = 4 * (i+1) ;
    RealOperatorAddLaderNormalOrder(op1, oadd, 2, val, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4 * i + 3;
    oadd[0] = 4 * (i+1) + 1;
    RealOperatorAddLaderNormalOrder(op1, oadd, 2, val, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4 * (i+1);
    oadd[0] = 4 * i + 2;
    RealOperatorAddLaderNormalOrder(op1, oadd, 2, val, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask + 4 * (i+1) + 1;
    oadd[0] = 4 * i + 3;
    RealOperatorAddLaderNormalOrder(op1, oadd, 2, val, &QDetNotOrderedFcFaBcBa);
  }
  return 0;
}


int CreateOperatorHarmonicOscillatorOnSplineWaveletAInteraction(double const d, uint32_t const nmax, OperatorType *const op2)
{
  // the basis functions are
  // g1[x] = (Sqrt[35/13]*Piecewise[{{1 - (3*x^2)/4 - x^3/4, x > -2 && x <= 0}, {1 - (3*x^2)/4 + x^3/4, x > 0 && x < 2}}, 0])/2
  // g2[x] = (Sqrt[105]*Piecewise[{{x*(2 + x)^2, x > -2 && x <= 0}, {(2 - x)^2*x, x > 0 && x < 2}}, 0])/16
  // f1[x,n,d] = sqrt[2/d] g1[2x/d-2n]
  // f2[x,n,d] = sqrt[2/d] g2[2x/d-2n]
  // index for interaction {f1[n], f2[n], f1[n+1], f2[n+1]}
  const double Uijkl[4][4][4][4] = {{{{ 0.96815905987503620639715314271527,   0,                                   0.072503310299168287333967807340588, -0.12596148731355991815742427845580},
                                      { 0,                                    0.65895642818719741796664873587951,  0.090035298959854471542591724341260, -0.15238367401828940290478752017214},
                                      { 0.072503310299168287333967807340588,  0.090035298959854471542591724341260, 0.043990772541068399056564737038110, -0.067141775087606857397300699909875},
                                      {-0.12596148731355991815742427845580,  -0.15238367401828940290478752017214, -0.067141775087606857397300699909875,  0.10511923973462435000896539358078}},
                                     {{ 0,                                    0.65895642818719741796664873587951,  0.090035298959854471542591724341260, -0.15238367401828940290478752017214},
                                      { 0.65895642818719741796664873587951,   0,                                   0.11923973462435000896539358077820,  -0.19718578380499693994277413224233},
                                      { 0.090035298959854471542591724341260,  0.11923973462435000896539358077820,  0.067141775087606857397300699909875, -0.10002017213555675094136632598171},
                                      {-0.15238367401828940290478752017214,  -0.19718578380499693994277413224233, -0.10002017213555675094136632598171,   0.15309455264363116455184327037448}},
                                     {{ 0.072503310299168287333967807340588,  0.090035298959854471542591724341260, 0.043990772541068399056564737038110, -0.067141775087606857397300699909875},
                                      { 0.090035298959854471542591724341260,  0.11923973462435000896539358077820,  0.067141775087606857397300699909875, -0.10002017213555675094136632598171},
                                      { 0.043990772541068399056564737038110,  0.067141775087606857397300699909875, 0.072503310299168287333967807340588, -0.090035298959854471542591724341260},
                                      {-0.067141775087606857397300699909875, -0.10002017213555675094136632598171, -0.090035298959854471542591724341260,  0.11923973462435000896539358077820}},
                                     {{-0.12596148731355991815742427845580,  -0.15238367401828940290478752017214, -0.067141775087606857397300699909875,  0.10511923973462435000896539358078},
                                      {-0.15238367401828940290478752017214,  -0.19718578380499693994277413224233, -0.10002017213555675094136632598171,   0.15309455264363116455184327037448},
                                      {-0.067141775087606857397300699909875, -0.10002017213555675094136632598171, -0.090035298959854471542591724341260,  0.11923973462435000896539358077820},
                                      { 0.10511923973462435000896539358078,   0.15309455264363116455184327037448,  0.11923973462435000896539358077820,  -0.16534211685512165771599073200444}}},
                                    {{{ 0,                                    0.65895642818719741796664873587951,  0.090035298959854471542591724341260, -0.15238367401828940290478752017214},
                                      { 0.65895642818719741796664873587951,   0,                                   0.11923973462435000896539358077820,  -0.19718578380499693994277413224233},
                                      { 0.090035298959854471542591724341260,  0.11923973462435000896539358077820,  0.067141775087606857397300699909875, -0.10002017213555675094136632598171},
                                      {-0.15238367401828940290478752017214,  -0.19718578380499693994277413224233, -0.10002017213555675094136632598171,   0.15309455264363116455184327037448}},
                                     {{ 0.65895642818719741796664873587951,   0,                                   0.11923973462435000896539358077820,  -0.19718578380499693994277413224233},
                                      { 0,                                    0.85664335664335664335664335664336,  0.16534211685512165771599073200444,  -0.26770104895104895104895104895105},
                                      { 0.11923973462435000896539358077820,   0.16534211685512165771599073200444,  0.10511923973462435000896539358078,  -0.15309455264363116455184327037448},
                                      {-0.19718578380499693994277413224233,  -0.26770104895104895104895104895105, -0.15309455264363116455184327037448,   0.22945804195804195804195804195804}},
                                     {{ 0.090035298959854471542591724341260,  0.11923973462435000896539358077820,  0.067141775087606857397300699909875, -0.10002017213555675094136632598171},
                                      { 0.11923973462435000896539358077820,   0.16534211685512165771599073200444,  0.10511923973462435000896539358078,  -0.15309455264363116455184327037448},
                                      { 0.067141775087606857397300699909875,  0.10511923973462435000896539358078,  0.12596148731355991815742427845580,  -0.15238367401828940290478752017214},
                                      {-0.10002017213555675094136632598171,  -0.15309455264363116455184327037448, -0.15238367401828940290478752017214,   0.19718578380499693994277413224233}},
                                     {{-0.15238367401828940290478752017214,  -0.19718578380499693994277413224233, -0.10002017213555675094136632598171,   0.15309455264363116455184327037448},
                                      {-0.19718578380499693994277413224233,  -0.26770104895104895104895104895105, -0.15309455264363116455184327037448,   0.22945804195804195804195804195804},
                                      {-0.10002017213555675094136632598171,  -0.15309455264363116455184327037448, -0.15238367401828940290478752017214,   0.19718578380499693994277413224233},
                                      { 0.15309455264363116455184327037448,   0.22945804195804195804195804195804,  0.19718578380499693994277413224233,  -0.26770104895104895104895104895105}}},
                                    {{{ 0.072503310299168287333967807340588,  0.090035298959854471542591724341260, 0.043990772541068399056564737038110, -0.067141775087606857397300699909875},
                                      { 0.090035298959854471542591724341260,  0.11923973462435000896539358077820,  0.067141775087606857397300699909875, -0.10002017213555675094136632598171},
                                      { 0.043990772541068399056564737038110,  0.067141775087606857397300699909875, 0.072503310299168287333967807340588, -0.090035298959854471542591724341260},
                                      {-0.067141775087606857397300699909875, -0.10002017213555675094136632598171, -0.090035298959854471542591724341260,  0.11923973462435000896539358077820}},
                                     {{ 0.090035298959854471542591724341260,  0.11923973462435000896539358077820,  0.067141775087606857397300699909875, -0.10002017213555675094136632598171},
                                      { 0.11923973462435000896539358077820,   0.16534211685512165771599073200444,  0.10511923973462435000896539358078,  -0.15309455264363116455184327037448},
                                      { 0.067141775087606857397300699909875,  0.10511923973462435000896539358078,  0.12596148731355991815742427845580,  -0.15238367401828940290478752017214},
                                      {-0.10002017213555675094136632598171,  -0.15309455264363116455184327037448, -0.15238367401828940290478752017214,   0.19718578380499693994277413224233}},
                                     {{ 0.043990772541068399056564737038110,  0.067141775087606857397300699909875, 0.072503310299168287333967807340588, -0.090035298959854471542591724341260},
                                      { 0.067141775087606857397300699909875,  0.10511923973462435000896539358078,  0.12596148731355991815742427845580,  -0.15238367401828940290478752017214},
                                      { 0.072503310299168287333967807340588,  0.12596148731355991815742427845580,  0.96815905987503620639715314271527,   0},
                                      {-0.090035298959854471542591724341260, -0.15238367401828940290478752017214,  0,                                    0.65895642818719741796664873587951}},
                                     {{-0.067141775087606857397300699909875, -0.10002017213555675094136632598171, -0.090035298959854471542591724341260,  0.11923973462435000896539358077820},
                                      {-0.10002017213555675094136632598171,  -0.15309455264363116455184327037448, -0.15238367401828940290478752017214,   0.19718578380499693994277413224233},
                                      {-0.090035298959854471542591724341260, -0.15238367401828940290478752017214,  0,                                    0.65895642818719741796664873587951},
                                      { 0.11923973462435000896539358077820,   0.19718578380499693994277413224233,  0.65895642818719741796664873587951,   0}}},
                                    {{{-0.12596148731355991815742427845580,  -0.15238367401828940290478752017214, -0.067141775087606857397300699909875,  0.10511923973462435000896539358078},
                                      {-0.15238367401828940290478752017214,  -0.19718578380499693994277413224233, -0.10002017213555675094136632598171,   0.15309455264363116455184327037448},
                                      {-0.067141775087606857397300699909875, -0.10002017213555675094136632598171, -0.090035298959854471542591724341260,  0.11923973462435000896539358077820},
                                      { 0.10511923973462435000896539358078,   0.15309455264363116455184327037448,  0.11923973462435000896539358077820,  -0.16534211685512165771599073200444}},
                                     {{-0.15238367401828940290478752017214,  -0.19718578380499693994277413224233, -0.10002017213555675094136632598171,   0.15309455264363116455184327037448},
                                      {-0.19718578380499693994277413224233,  -0.26770104895104895104895104895105, -0.15309455264363116455184327037448,   0.22945804195804195804195804195804},
                                      {-0.10002017213555675094136632598171,  -0.15309455264363116455184327037448, -0.15238367401828940290478752017214,   0.19718578380499693994277413224233},
                                      { 0.15309455264363116455184327037448,   0.22945804195804195804195804195804,  0.19718578380499693994277413224233,  -0.26770104895104895104895104895105}},
                                     {{-0.067141775087606857397300699909875, -0.10002017213555675094136632598171, -0.090035298959854471542591724341260,  0.11923973462435000896539358077820},
                                      {-0.10002017213555675094136632598171,  -0.15309455264363116455184327037448, -0.15238367401828940290478752017214,   0.19718578380499693994277413224233},
                                      {-0.090035298959854471542591724341260, -0.15238367401828940290478752017214,  0,                                    0.65895642818719741796664873587951},
                                      { 0.11923973462435000896539358077820,   0.19718578380499693994277413224233,  0.65895642818719741796664873587951,   0}},
                                     {{ 0.10511923973462435000896539358078,   0.15309455264363116455184327037448,  0.11923973462435000896539358077820,  -0.16534211685512165771599073200444},
                                      { 0.15309455264363116455184327037448,   0.22945804195804195804195804195804,  0.19718578380499693994277413224233,  -0.26770104895104895104895104895105},
                                      { 0.11923973462435000896539358077820,   0.19718578380499693994277413224233,  0.65895642818719741796664873587951,   0},
                                      {-0.16534211685512165771599073200444,  -0.26770104895104895104895104895105,  0,                                    0.85664335664335664335664335664336}}}};
  uint16_t i,j,k,l,n,oadd[4];
  InitNamedOperator(op2, "Harmonic Oscillator Point-Contact");
  op2->NFermion = 2 * (4*nmax+2);
  // number of words (meaning the number of 8-bit blocks needed for storing bitstrings)
  op2->NWord = Min(op2->NFermion, 1 + (op2->NFermion - 1) / 8);
  InitOperatorFixedLength(op2, 4, (uint32_t)floor(log(nmax * nmax * nmax) / log(2.0)));
  // interaction on site n
  for (n = 0; n < 2*nmax+1; n++)
  {
    for (i=0; i<2; i++)
    {
      for (j=0; j<2; j++)
      {
        for (k=0; k<2; k++)
        {
          for (l=0; l<2; l++)
          {
            oadd[3] = CrAnMask + 4 * n + 2*i;
            oadd[2] = 4 * n + 2*j;
            oadd[1] = CrAnMask + 4 * n + 2*k + 1;
            oadd[0] = 4 * n + 2*l + 1;
            RealOperatorAddLaderNormalOrder(op2, oadd, 2, Uijkl[i][j][k][l], &QDetNotOrderedFcFaBcBa);
          }
        }
      }
    }
  }
  
  return 0;
}
