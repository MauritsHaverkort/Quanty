/*  A supercell is in principle defined if the unit cell is defined and the number
 *  Of repetitions is given. However it is nice to also have an index of the postions of atoms and g-vectors within the 
 *  supercell. 
 */
uint32_t AddAditionalInformationToSuperCell(CrystalStructureSuperCellType * A)
{
  int32_t ivec[3], icell, i, j, k, l;
  double UnitCellDouble[3][3], a, b, c;
  icell = 0;
  A->NCell = ParallelalPipedVolumeMatrix(A->rSuperCell);
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      UnitCellDouble[i][j] = (double)A->rSuperCell[i][j];
    }
  }
  ResiprocalLatticeVector(UnitCellDouble, A->gSuperCell,A->Crystal.Units);
  l = 0;
  for (i = -(int32_t) A->NCell; i < (int32_t) A->NCell; i++)
  {
    for (j = -(int32_t) A->NCell; j < (int32_t) A->NCell; j++)
    {
      for (k = -(int32_t) A->NCell; k < (int32_t) A->NCell; k++)
      {
        a = A->gSuperCell[0][0] * i + A->gSuperCell[1][0] * j + A->gSuperCell[2][0] * k;
        b = A->gSuperCell[0][1] * i + A->gSuperCell[1][1] * j + A->gSuperCell[2][1] * k;
        c = A->gSuperCell[0][2] * i + A->gSuperCell[1][2] * j + A->gSuperCell[2][2] * k;
        if(A->Crystal.Units&1)
        {
          if ((a >= 0) && (a < 1.999999999 * Pi) && (b >= 0) && (b < 1.999999999 * Pi) && (c >= 0) && (c < 1.999999999 * Pi))
          {
            A->kVector[l][0] = a;
            A->kVector[l][1] = b;
            A->kVector[l][2] = c;
            l++;
          }
        }
        else
        {
          if ((a >= 0) && (a < 0.999999999) && (b >= 0) && (b < 0.999999999) && (c >= 0) && (c < 0.999999999))
          {
            A->kVector[l][0] = a;
            A->kVector[l][1] = b;
            A->kVector[l][2] = c;
            l++;
          }
        }
      }
    }
  }
  if (l != A->NCell)
  {
    printf("NCell=%i NkVector=%i thats not the same and that is a problem !!!\n\n\nSTOP", A->NCell, l);
    PrintCrystalStructureSuperCell(*A);
    exit(1);
    return 1;
  }

  icell = 0;
  for (ivec[0] = -(abs(A->rSuperCell[0][0]) + abs(A->rSuperCell[1][0]) + abs(A->rSuperCell[2][0]));
       ivec[0] < (abs(A->rSuperCell[0][0]) + abs(A->rSuperCell[1][0]) + abs(A->rSuperCell[2][0])); ivec[0]++)
  {
    for (ivec[1] = -(abs(A->rSuperCell[0][1]) + abs(A->rSuperCell[1][1]) + abs(A->rSuperCell[2][1]));
         ivec[1] < (abs(A->rSuperCell[0][1]) + abs(A->rSuperCell[1][1]) + abs(A->rSuperCell[2][1])); ivec[1]++)
    {
      for (ivec[2] = -(abs(A->rSuperCell[0][2]) + abs(A->rSuperCell[1][2]) + abs(A->rSuperCell[2][2]));
           ivec[2] < (abs(A->rSuperCell[0][2]) + abs(A->rSuperCell[1][2]) + abs(A->rSuperCell[2][2])); ivec[2]++)
      {
        a = Dot(ivec, A->gSuperCell[0]);
        b = Dot(ivec, A->gSuperCell[1]);
        c = Dot(ivec, A->gSuperCell[2]);
        if(A->Crystal.Units&1)
        {
          if ((a >= 0) && (a < 1.999999999 * Pi) && (b >= 0) && (b < 1.999999999 * Pi) && (c >= 0) && (c < 1.999999999 * Pi))
          {
            A->CellPos[icell][0] = ivec[0];
            A->CellPos[icell][1] = ivec[1];
            A->CellPos[icell][2] = ivec[2];
            icell++;
          }
        }
        else
        {
          if ((a >= 0) && (a < 0.999999999) && (b >= 0) && (b < 0.999999999) && (c >= 0) && (c < 0.999999999))
          {
            A->CellPos[icell][0] = ivec[0];
            A->CellPos[icell][1] = ivec[1];
            A->CellPos[icell][2] = ivec[2];
            icell++;
          }
        }
      }
    }
  }
  if (icell != A->NCell)
  {
    printf("NCell=%i icell=%i thats not the same and that is a problem !!!\n\n\nSTOP", A->NCell, icell);
    PrintCrystalStructureSuperCell(*A);
    exit(1);
  }
  return 0;
}

uint32_t CopyCrystalStructure(CrystalStructureType A, CrystalStructureType * B)
{
  uint32_t i, j;
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      B->rUnitCell[i][j] = A.rUnitCell[i][j];
      B->gUnitCell[i][j] = A.gUnitCell[i][j];
    }
  }
  B->NAtoms = A.NAtoms;
  if ((B->Atom = malloc(B->NAtoms * sizeof(AtomType))) == NULL)
  {
    printf("malloc failed for B->Atom in CopyCrystalStructure\n");
    fflush(stdout);
    return 1;
  };
  for (i = 0; i < B->NAtoms; i++)
  {
    strcpy(B->Atom[i].Name, A.Atom[i].Name);
    B->Atom[i].AtomNumber = A.Atom[i].AtomNumber;
    B->Atom[i].Position[0] = A.Atom[i].Position[0];
    B->Atom[i].Position[1] = A.Atom[i].Position[1];
    B->Atom[i].Position[2] = A.Atom[i].Position[2];
    B->Atom[i].Nshell = A.Atom[i].Nshell;
    if ((B->Atom[i].shelldegeneracy = malloc(B->Atom[i].Nshell * sizeof(int))) == NULL)
    {
      printf("malloc failed for B->Atom[%i].shelldegeneracy in CopyCrystalStructure\n", i);
      fflush(stdout);
      return 1;
    };
    for (j = 0; j < B->Atom[i].Nshell; j++)
    {
      B->Atom[i].shelldegeneracy[j] = A.Atom[i].shelldegeneracy[j];
    }
  }
  return 0;
}

uint32_t ShortesConnectingVectorPeriodicSuperCell(const int32_t A[], const int32_t B[],
                                                  const CrystalStructureSuperCellType SC, int32_t C[])
{
  int32_t imin, jmin, kmin, i, j, k, Ctest[3];
  double UnitCellDoubleA[3], UnitCellDoubleB[3], UnitCellDoubleC[3];
  double temp[3];
  double denom;
  C[0] = A[0] - B[0];
  C[1] = A[1] - B[1];
  C[2] = A[2] - B[2];
  for (i = 0; i < 3; i++)
  {
    UnitCellDoubleA[i] = (double)SC.rSuperCell[0][i];
    UnitCellDoubleB[i] = (double)SC.rSuperCell[1][i];
    UnitCellDoubleC[i] = (double)SC.rSuperCell[2][i];
  }
  denom = fabs(ParallelalPipedVolume(UnitCellDoubleA, UnitCellDoubleB, UnitCellDoubleC));
  Cross(UnitCellDoubleB, UnitCellDoubleC, temp);
  imin = 1 + (int32_t) (denom / sqrt(Dot(UnitCellDoubleA, UnitCellDoubleA) * Dot(temp, temp)));
  Cross(UnitCellDoubleA, UnitCellDoubleC, temp);
  jmin = 1 + (int32_t) (denom / sqrt(Dot(UnitCellDoubleB, UnitCellDoubleB) * Dot(temp, temp)));
  Cross(UnitCellDoubleA, UnitCellDoubleB, temp);
  kmin = 1 + (int32_t) (denom / sqrt(Dot(UnitCellDoubleC, UnitCellDoubleC) * Dot(temp, temp)));
  for (i = -imin; i <= imin; i++)
  {
    for (j = -jmin; j <= jmin; j++)
    {
      for (k = -kmin; k <= kmin; k++)
      {
        Ctest[0] = A[0] - B[0] + i * SC.rSuperCell[0][0] + j * SC.rSuperCell[1][0] + k * SC.rSuperCell[2][0];
        Ctest[1] = A[1] - B[1] + i * SC.rSuperCell[0][1] + j * SC.rSuperCell[1][1] + k * SC.rSuperCell[2][1];
        Ctest[2] = A[2] - B[2] + i * SC.rSuperCell[0][2] + j * SC.rSuperCell[1][2] + k * SC.rSuperCell[2][2];
        if (Dot(Ctest, Ctest) < Dot(C, C))
        {
          C[0] = Ctest[0];
          C[1] = Ctest[1];
          C[2] = Ctest[2];
        }
      }
    }
  }
  return 0;
}

uint32_t FindCellPosInCrystalStructureSuperCellType(CrystalStructureSuperCellType A, int32_t vec[3], int *pos)
{
  uint32_t i;
  *pos = A.NCell;
  for (i = 0; i < *pos; i++)
  {
    if ((A.CellPos[i][0] == vec[0]) && (A.CellPos[i][1] == vec[1]) && (A.CellPos[i][2] == vec[2])) *pos = i;
  }
  if (*pos == A.NCell)
  {
    printf("vector %i %i %i not found in super cell positions\n", vec[0], vec[1], vec[2]);
    return 1;
  }
  return 0;
}

/* 
* For k resolved spectroscopy one wants to calculate <i|\Sum_r` e^{- i k . r`} T^\dag_r` G \sum_r e^{i k . r} T_r |i>
* However it is simpler to calculate <i|T_(r)+T^\dag_(r+dr) G T_(r)+T^\dag_(r+dr)|i>
* This function assumes T(dr) as input for the list spectra and replaces this by the fourier transform.
* There are Natom atoms and therefor also Natom spectra assumed the first one belongs to dr=0
* There must thus also be Natom positions listed. The first spectrum is taken at the first position.
* The n-th spectrum is taken as the sum of T_atom_0 +T_atom_n
*/
uint32_t SpectraSiteFourierTransform(const CrystalStructureSuperCellType SC, SpectraType * Spectra)
{
  int32_t i, j, k;
  SpectraType SpectraDeltaR[SC.NCell];
  for (j = 0; j < SC.NCell; j++)
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
  for (k = 0; k < SC.NCell; k++)
  {
    SetSpectraToZero(&Spectra[k]);
    AddMultipliedSpectra(&Spectra[k], SpectraDeltaR[j], cos(Dot(SC.kVector[k], SC.CellPos[j])),
                         sin(Dot(SC.kVector[k], SC.CellPos[j])));
  }
  return 0;
}

/*  
 *
 */
uint32_t OperatorSiteFourierTransform(const CrystalStructureSuperCellType SC, const OperatorType Oi, OperatorType * Ok,
                                      uint32_t(*QDetNotOrdered) (uint16_t, uint16_t, uint32_t))
{
  //TODO: test if function servived changes in MappedOperatorAddNormalOrder
  uint32_t i, j, k, FermionDegeneracy, BosonDegeneracy;
  MatrixType U;
  U.Ni = Oi.NFermion + Oi.NBoson;
  U.Nj = Oi.NFermion + Oi.NBoson;
  U.QComplex = 1;
  InitMatrix(&U);
  InitOperatorCopySettings(Ok, Oi);
  FermionDegeneracy = Oi.NFermion / SC.NCell;
  BosonDegeneracy = Oi.NBoson / SC.NCell;
  if ((SC.NCell * FermionDegeneracy) != Oi.NFermion)
  {
    printf("Number of Fermions not a multiple of the number of sites\n");
    printf("Error in OperatorSiteFourierTransform\n");
    return 1;
  }
  if ((SC.NCell * BosonDegeneracy) != Oi.NBoson)
  {
    printf("Number of Bosons not a multiple of the number of sites\n");
    printf("Error in OperatorSiteFourierTransform\n");
    return 1;
  }
  for (i = 0; i < SC.NCell; i++)
  {
    for (j = 0; j < SC.NCell; j++)
    {
      for (k = 0; k < FermionDegeneracy; k++)
      {
        U.VRe[i * FermionDegeneracy + k][j * FermionDegeneracy + k] = cos(Dot(SC.kVector[j], SC.CellPos[i]));
        U.VIm[i * FermionDegeneracy + k][j * FermionDegeneracy + k] = sin(Dot(SC.kVector[j], SC.CellPos[i]));
      }
      for (k = 0; k < BosonDegeneracy; k++)
      {
        U.VRe[i * BosonDegeneracy + k + Oi.NFermion][j * BosonDegeneracy + k + Oi.NFermion] =
          cos(Dot(SC.kVector[j], SC.CellPos[i]));
        U.VIm[i * BosonDegeneracy + k + Oi.NFermion][j * BosonDegeneracy + k + Oi.NFermion] =
          sin(Dot(SC.kVector[j], SC.CellPos[i]));
      }
    }
  }
  return MappedOperatorAddNormalOrder(Ok, Oi, U, QDetNotOrdered);
}
