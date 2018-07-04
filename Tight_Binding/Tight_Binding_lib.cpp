#include "Tight_Binding_lib.h"

uint32_t AtomNameToNumber(const char* Name)
{
  uint32_t i;
  const char *PeriodicTable[] = {"H" ,                                                                                                                                                          "He",
                           "Li","Be",                                                                                                                            "B", "C", "N", "O", "F", "Ne",
                           "Na","Mg",                                                                                                                           "Al","Si", "P", "S","Cl", "Ar",
                            "K","Ca",                                                                        "Sc","Ti", "V","Cr","Mn","Fe","Co","Ni","Cu","Zn", "Ga","Ge","As","Se","Br", "Kr",
                           "Rb","Sr",                                                                         "Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd", "In","Sn","Sb","Te", "I", "Xe",
                           "Cs","Ba", "La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb", "Lu","Hf","Ta", "W","Re","Os","Ir","Pt","Au","Hg", "Tl","Pb","Bi","Po","At", "Rn",
                           "Fr","Ra", "Ac","Th","Pa", "U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No", "Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn"};
  for(i=0;i<36;i++)
  {
    if(strcmp(Name, PeriodicTable[i] )==0) return i+1;
  }
  return 0;
}

int FreeAtomType(AtomType *A)
{
  int i,j;
  for(i=0;i<A->Nshell;i++)
  {
    free(A->ShellName[i]);
    for(j=0;j<A->shelldegeneracy[i];j++)
    {
      free(A->OrbitalName[i][j]);
    }
    free(A->OrbitalName[i]);
  }
  if(A->Nshell!=0)
  {
    free(A->ShellName);
    free(A->shelldegeneracy);
    free(A->OrbitalName);
  }
  A->Nshell=0;
  return 0;
}

int FreeCrystalStructureType(CrystalStructureType *A)
{
  int i;
  for(i=0;i<A->NAtoms;i++)
  {
    FreeAtomType(&(A->Atom[i]));
  }
  if(A->NAtoms!=0) free(A->Atom);
  A->NAtoms=0;
  return 0;
}

int FreeTightBindingDef(TightBindingDef * A)
{
  int    i;
  for(i=0;i<A->nH;i++)
  {
    FreeMatrix(&(A->H[i].Hop));
  }
  if(A->nH!=0) free(A->H);
  A->nH=0;
  if(A->FermionIndex!=NULL)
  {
    for (i = 0; i < A->Cell.NAtoms; i++)
    {
      free(A->FermionIndex[i]);
    }
    if(A->Cell.NAtoms!=0) free(A->FermionIndex);
  }
  A->Cell.NAtoms=0;
  FreeCrystalStructureType(&(A->Cell));
  return 0;
}




/*  The tight-binding structure defines a unit cell (CrystalStructureType) including a basis
 *  This basis is given by atoms (AtomType) at certain positions
 *  An atom has certain shells which have a degeneracy
 *  One often needs to know the total number of fermions in the shell
 *  The list FermionIndex[AtomNumber][ShellNumber] stores the number of the first fermion within this shell.
 */
int SetNFermionIndexTightBindingDef(TightBindingDef * A)
{
  int    i, j;
  A->NFermion = 0;
  A->FermionIndex = (int **)malloc(A->Cell.NAtoms * sizeof(int *));
  for (i = 0; i < A->Cell.NAtoms; i++)
  {
    A->FermionIndex[i] = (int*)malloc(A->Cell.Atom[i].Nshell * sizeof(int));
    for (j = 0; j < A->Cell.Atom[i].Nshell; j++)
    {
      A->FermionIndex[i][j] = A->NFermion;
      A->NFermion += A->Cell.Atom[i].shelldegeneracy[j];
    }
  }
  return 0;
}

int CreateOneBandCubicLaticeTightBindingDef(TightBindingDef * A)
{
  int    i;
  strncpy(A->Name, "Cubic Nearest Neighbour tb", 30);
  A->Cell.rUnitCell[0][0] = 1.00000;
  A->Cell.rUnitCell[0][1] = 0.00000;
  A->Cell.rUnitCell[0][2] = 0.00000;
  A->Cell.rUnitCell[1][0] = 0.00000;
  A->Cell.rUnitCell[1][1] = 1.00000;
  A->Cell.rUnitCell[1][2] = 0.00000;
  A->Cell.rUnitCell[2][0] = 0.00000;
  A->Cell.rUnitCell[2][1] = 0.00000;
  A->Cell.rUnitCell[2][2] = 1.00000;
  A->Cell.Units = 1;
  ReciprocalLatticeVector(A->Cell.rUnitCell, A->Cell.gUnitCell,A->Cell.Units);
  A->Cell.NAtoms = 1;
  A->Cell.Atom = (AtomType*)malloc(1 * sizeof(AtomType));
  strncpy(A->Cell.Atom[0].Name, "s", 3);
  A->Cell.Atom[0].Position[0] = 0.00000;
  A->Cell.Atom[0].Position[1] = 0.00000;
  A->Cell.Atom[0].Position[2] = 0.00000;
  A->Cell.Atom[0].Nshell = 1;
  A->Cell.Atom[0].shelldegeneracy = (int*)malloc(1 * sizeof(int));
  A->Cell.Atom[0].shelldegeneracy[0] = 1;
  A->Cell.Atom[0].ShellName = (char**)malloc(1*sizeof(char*));
  A->Cell.Atom[0].ShellName[0] = (char*)malloc(1*sizeof(char));
  strncpy(A->Cell.Atom[0].ShellName[0], "s", 3);
  A->Cell.Atom[0].OrbitalName= (char***)malloc(1*sizeof(char**));
  A->Cell.Atom[0].OrbitalName[0]= (char**)malloc(1*sizeof(char*));
  A->Cell.Atom[0].OrbitalName[0][0]= (char*)malloc(1*sizeof(char));
  strncpy(A->Cell.Atom[0].OrbitalName[0][0], "s", 3);
  A->nH = 7;
  A->H = (TightBindingDefHopping*)malloc(A->nH * sizeof(TightBindingDefHopping));
  for (i = 0; i < 7; i++)
  {
    A->H[i].Atompos1 = 0;
    A->H[i].Shellpos1 = 0;
    A->H[i].Atompos2 = 0;
    A->H[i].Shellpos2 = 0;
    A->H[i].TransVecs[0] = 0;
    A->H[i].TransVecs[1] = 0;
    A->H[i].TransVecs[2] = 0;
    if (i > 0)
      A->H[i].TransVecs[(i - 1) / 2] = 1 - 2 * (i & 1);
    A->H[i].Hop.Ni = A->Cell.Atom[A->H[i].Atompos1].shelldegeneracy[A->H[i].Shellpos1];
    A->H[i].Hop.Nj = A->Cell.Atom[A->H[i].Atompos2].shelldegeneracy[A->H[i].Shellpos2];
    A->H[i].Hop.QComplex = 0;
    InitMatrix(&(A->H[i].Hop));
    if (i > 0)
      A->H[i].Hop.VRe[0][0] = -1.0 / 12.0;
  }
  SetNFermionIndexTightBindingDef(A);
  return 0;
}


int CreateOneBandSquareLaticeTightBindingDef(TightBindingDef * A)
{
  int    i;
  strncpy(A->Name, "Square Nearest Neighbour tb", 30);
  A->Cell.rUnitCell[0][0] = 1.00000;
  A->Cell.rUnitCell[0][1] = 0.00000;
  A->Cell.rUnitCell[0][2] = 0.00000;
  A->Cell.rUnitCell[1][0] = 0.00000;
  A->Cell.rUnitCell[1][1] = 1.00000;
  A->Cell.rUnitCell[1][2] = 0.00000;
  A->Cell.rUnitCell[2][0] = 0.00000;
  A->Cell.rUnitCell[2][1] = 0.00000;
  A->Cell.rUnitCell[2][2] = 4.00000;
  A->Cell.Units = 1;
  ReciprocalLatticeVector(A->Cell.rUnitCell, A->Cell.gUnitCell,A->Cell.Units);
  A->Cell.NAtoms = 1;
  A->Cell.Atom = (AtomType*)malloc(1 * sizeof(AtomType));
  strncpy(A->Cell.Atom[0].Name, "s", 3);
  A->Cell.Atom[0].Position[0] = 0.00000;
  A->Cell.Atom[0].Position[1] = 0.00000;
  A->Cell.Atom[0].Position[2] = 0.00000;
  A->Cell.Atom[0].Nshell = 1;
  A->Cell.Atom[0].shelldegeneracy = (int*)malloc(1 * sizeof(int));
  A->Cell.Atom[0].shelldegeneracy[0] = 1;
  A->Cell.Atom[0].ShellName = (char**)malloc(1*sizeof(char*));
  A->Cell.Atom[0].ShellName[0] = (char*)malloc(1*sizeof(char));
  strncpy(A->Cell.Atom[0].ShellName[0], "s", 3);
  A->Cell.Atom[0].OrbitalName= (char***)malloc(1*sizeof(char**));
  A->Cell.Atom[0].OrbitalName[0]= (char**)malloc(1*sizeof(char*));
  A->Cell.Atom[0].OrbitalName[0][0]= (char*)malloc(1*sizeof(char));
  strncpy(A->Cell.Atom[0].OrbitalName[0][0], "s", 3);
  A->nH = 5;
  A->H = (TightBindingDefHopping*)malloc(A->nH * sizeof(TightBindingDefHopping));
  for (i = 0; i < 5; i++)
  {
    A->H[i].Atompos1 = 0;
    A->H[i].Shellpos1 = 0;
    A->H[i].Atompos2 = 0;
    A->H[i].Shellpos2 = 0;
    A->H[i].TransVecs[0] = 0;
    A->H[i].TransVecs[1] = 0;
    A->H[i].TransVecs[2] = 0;
    if (i > 0)
      A->H[i].TransVecs[(i - 1) / 2] = 1 - 2 * (i & 1);
    A->H[i].Hop.Ni = A->Cell.Atom[A->H[i].Atompos1].shelldegeneracy[A->H[i].Shellpos1];
    A->H[i].Hop.Nj = A->Cell.Atom[A->H[i].Atompos2].shelldegeneracy[A->H[i].Shellpos2];
    A->H[i].Hop.QComplex = 0;
    InitMatrix(&(A->H[i].Hop));
    if (i > 0)
      A->H[i].Hop.VRe[0][0] = -1.0 / 8.0;
  }
  SetNFermionIndexTightBindingDef(A);
  return 0;
}



int CreateOneBandChainLaticeTightBindingDef(TightBindingDef * A)
{
  int    i;
  strncpy(A->Name, "chain Nearest Neighbour tb", 30);
  A->Cell.rUnitCell[0][0] = 1.00000;
  A->Cell.rUnitCell[0][1] = 0.00000;
  A->Cell.rUnitCell[0][2] = 0.00000;
  A->Cell.rUnitCell[1][0] = 0.00000;
  A->Cell.rUnitCell[1][1] = 4.00000;
  A->Cell.rUnitCell[1][2] = 0.00000;
  A->Cell.rUnitCell[2][0] = 0.00000;
  A->Cell.rUnitCell[2][1] = 0.00000;
  A->Cell.rUnitCell[2][2] = 4.00000;
  A->Cell.Units = 1;
  ReciprocalLatticeVector(A->Cell.rUnitCell, A->Cell.gUnitCell,A->Cell.Units);
  A->Cell.NAtoms = 1;
  A->Cell.Atom = (AtomType*)malloc(1 * sizeof(AtomType));
  strncpy(A->Cell.Atom[0].Name, "s", 3);
  A->Cell.Atom[0].Position[0] = 0.00000;
  A->Cell.Atom[0].Position[1] = 0.00000;
  A->Cell.Atom[0].Position[2] = 0.00000;
  A->Cell.Atom[0].Nshell = 1;
  A->Cell.Atom[0].shelldegeneracy = (int*)malloc(1 * sizeof(int));
  A->Cell.Atom[0].shelldegeneracy[0] = 1;
  A->Cell.Atom[0].ShellName = (char**)malloc(1*sizeof(char*));
  A->Cell.Atom[0].ShellName[0] = (char*)malloc(1*sizeof(char));
  strncpy(A->Cell.Atom[0].ShellName[0], "s", 3);
  A->Cell.Atom[0].OrbitalName= (char***)malloc(1*sizeof(char**));
  A->Cell.Atom[0].OrbitalName[0]= (char**)malloc(1*sizeof(char*));
  A->Cell.Atom[0].OrbitalName[0][0]= (char*)malloc(1*sizeof(char));
  strncpy(A->Cell.Atom[0].OrbitalName[0][0], "s", 3);
  A->nH = 3;
  A->H = (TightBindingDefHopping*)malloc(A->nH * sizeof(TightBindingDefHopping));
  for (i = 0; i < 3; i++)
  {
    A->H[i].Atompos1 = 0;
    A->H[i].Shellpos1 = 0;
    A->H[i].Atompos2 = 0;
    A->H[i].Shellpos2 = 0;
    A->H[i].TransVecs[0] = 0;
    A->H[i].TransVecs[1] = 0;
    A->H[i].TransVecs[2] = 0;
    if (i > 0)
      A->H[i].TransVecs[(i - 1) / 2] = 1 - 2 * (i & 1);
    A->H[i].Hop.Ni = A->Cell.Atom[A->H[i].Atompos1].shelldegeneracy[A->H[i].Shellpos1];
    A->H[i].Hop.Nj = A->Cell.Atom[A->H[i].Atompos2].shelldegeneracy[A->H[i].Shellpos2];
    A->H[i].Hop.QComplex = 0;
    InitMatrix(&(A->H[i].Hop));
    if (i > 0)
      A->H[i].Hop.VRe[0][0] = -1.0 / 4.0;
  }
  SetNFermionIndexTightBindingDef(A);
  return 0;
}


/*  given a Tight binding definition and a k vector creates a matrix
 *  representing the Hamiltonian.
 *    two atoms at different positions in the unit cell hop with a phase e^{i k x_1 - x_2}
 */
int CreateHamiltonianFromTightBindingDef(TightBindingDef A, double k[3], MatrixType * H)
{
  int    iH, i, j;
  double alphaRe, alphaIm;
  double vector[3];
  H->Ni = A.NFermion;
  H->Nj = A.NFermion;
  H->QComplex = 0;
  InitMatrix(H);
  for (iH = 0; iH < A.nH; iH++)
  {
    i = A.FermionIndex[A.H[iH].Atompos1][A.H[iH].Shellpos1];
    j = A.FermionIndex[A.H[iH].Atompos2][A.H[iH].Shellpos2];
    vector[0] = A.Cell.Atom[A.H[iH].Atompos1].Position[0] - A.Cell.Atom[A.H[iH].Atompos2].Position[0];
    vector[1] = A.Cell.Atom[A.H[iH].Atompos1].Position[1] - A.Cell.Atom[A.H[iH].Atompos2].Position[1];
    vector[2] = A.Cell.Atom[A.H[iH].Atompos1].Position[2] - A.Cell.Atom[A.H[iH].Atompos2].Position[2];
    vector[0] +=
      A.Cell.rUnitCell[0][0] * A.H[iH].TransVecs[0] + A.Cell.rUnitCell[1][0] * A.H[iH].TransVecs[1] +
      A.Cell.rUnitCell[2][0] * A.H[iH].TransVecs[2];
    vector[1] +=
      A.Cell.rUnitCell[0][1] * A.H[iH].TransVecs[0] + A.Cell.rUnitCell[1][1] * A.H[iH].TransVecs[1] +
      A.Cell.rUnitCell[2][1] * A.H[iH].TransVecs[2];
    vector[2] +=
      A.Cell.rUnitCell[0][2] * A.H[iH].TransVecs[0] + A.Cell.rUnitCell[1][2] * A.H[iH].TransVecs[1] +
      A.Cell.rUnitCell[2][2] * A.H[iH].TransVecs[2];
    alphaRe = cos(Dot(k, vector));
    alphaIm = sin(Dot(k, vector));
    if (Sqr(alphaIm) < Sqr(alphaRe) * DBL_EPSILON * 100)
      alphaIm = 0;
    MatrixAddShifted(H, A.H[iH].Hop, alphaRe, alphaIm, i, j);
  }
  return 0;
}



/*  given a Tight binding definition and a k vector creates a matrix
 *  representing the Hamiltonian.
 *    two atoms at different positions in the unit cell hop with a constant phase 1
 */
int CreateHamiltonianFromTightBindingDefConstantPhaseInUnitCell(TightBindingDef A, double k[3], MatrixType * H)
{
  int    iH, i, j;
  double alphaRe, alphaIm;
  double vector[3];
  H->Ni = A.NFermion;
  H->Nj = A.NFermion;
  H->QComplex = 0;
  InitMatrix(H);
  for (iH = 0; iH < A.nH; iH++)
  {
    i = A.FermionIndex[A.H[iH].Atompos1][A.H[iH].Shellpos1];
    j = A.FermionIndex[A.H[iH].Atompos2][A.H[iH].Shellpos2];
    //vector[0] = A.Cell.Atom[A.H[iH].Atompos1].Position[0] - A.Cell.Atom[A.H[iH].Atompos2].Position[0];
    //vector[1] = A.Cell.Atom[A.H[iH].Atompos1].Position[1] - A.Cell.Atom[A.H[iH].Atompos2].Position[1];
    //vector[2] = A.Cell.Atom[A.H[iH].Atompos1].Position[2] - A.Cell.Atom[A.H[iH].Atompos2].Position[2];
    vector[0] = A.Cell.rUnitCell[0][0] * A.H[iH].TransVecs[0] + A.Cell.rUnitCell[1][0] * A.H[iH].TransVecs[1] + A.Cell.rUnitCell[2][0] * A.H[iH].TransVecs[2];
    vector[1] = A.Cell.rUnitCell[0][1] * A.H[iH].TransVecs[0] + A.Cell.rUnitCell[1][1] * A.H[iH].TransVecs[1] + A.Cell.rUnitCell[2][1] * A.H[iH].TransVecs[2];
    vector[2] = A.Cell.rUnitCell[0][2] * A.H[iH].TransVecs[0] + A.Cell.rUnitCell[1][2] * A.H[iH].TransVecs[1] + A.Cell.rUnitCell[2][2] * A.H[iH].TransVecs[2];
    alphaRe = cos(Dot(k, vector));
    alphaIm = sin(Dot(k, vector));
    if (Sqr(alphaIm) < Sqr(alphaRe) * DBL_EPSILON * 100) alphaIm = 0;
    MatrixAddShifted(H, A.H[iH].Hop, alphaRe, alphaIm, i, j);
  }
  return 0;
}

/*  given a Tight binding definition and a k vector creates an operator
 *  representing the Hamiltonian.
 */
int CreateHamiltonianOperatorFromTightBindingDef(TightBindingDef A, double k[3], OperatorType * H)
{
  int    iH, i, j, isub, jsub;
  double alphaRe, alphaIm, HVRe, HVIm;
  double vector[3];
  uint16_t oadd[2];
  InitOperator(H);
  H->NFermion = A.NFermion;
  H->NWord = Min(H->NFermion, 1 + (H->NFermion - 1) / 8);
  InitOperatorFixedLength(H, 2, 8);
  for (iH = 0; iH < A.nH; iH++)
  {
    i = A.FermionIndex[A.H[iH].Atompos1][A.H[iH].Shellpos1];
    j = A.FermionIndex[A.H[iH].Atompos2][A.H[iH].Shellpos2];
    vector[0] = A.Cell.Atom[A.H[iH].Atompos1].Position[0] - A.Cell.Atom[A.H[iH].Atompos2].Position[0];
    vector[1] = A.Cell.Atom[A.H[iH].Atompos1].Position[1] - A.Cell.Atom[A.H[iH].Atompos2].Position[1];
    vector[2] = A.Cell.Atom[A.H[iH].Atompos1].Position[2] - A.Cell.Atom[A.H[iH].Atompos2].Position[2];
    vector[0] += A.Cell.rUnitCell[0][0] * A.H[iH].TransVecs[0] + A.Cell.rUnitCell[1][0] * A.H[iH].TransVecs[1] + A.Cell.rUnitCell[2][0] * A.H[iH].TransVecs[2];
    vector[1] += A.Cell.rUnitCell[0][1] * A.H[iH].TransVecs[0] + A.Cell.rUnitCell[1][1] * A.H[iH].TransVecs[1] + A.Cell.rUnitCell[2][1] * A.H[iH].TransVecs[2];
    vector[2] += A.Cell.rUnitCell[0][2] * A.H[iH].TransVecs[0] + A.Cell.rUnitCell[1][2] * A.H[iH].TransVecs[1] + A.Cell.rUnitCell[2][2] * A.H[iH].TransVecs[2];
    alphaRe = cos(Dot(k, vector));
    alphaIm = sin(Dot(k, vector));
    if (Sqr(alphaIm) < Sqr(alphaRe) * DBL_EPSILON * 100) alphaIm = 0;
    if (A.H[iH].Hop.QComplex)
    {                           /* A.H[iH].Hop is Complex */
      for (isub = 0; isub < A.H[iH].Hop.Ni; isub++)
      {
        for (jsub = 0; jsub < A.H[iH].Hop.Nj; jsub++)
        {
          if ((fabs(A.H[iH].Hop.VRe[isub][jsub]) + fabs(A.H[iH].Hop.VIm[isub][jsub])) > DBL_MIN * 10)
          {
            HVRe = A.H[iH].Hop.VRe[isub][jsub] * alphaRe - A.H[iH].Hop.VIm[isub][jsub] * alphaIm;
            HVIm = A.H[iH].Hop.VRe[isub][jsub] * alphaIm + A.H[iH].Hop.VIm[isub][jsub] * alphaRe;
            oadd[1] = CrAnMask + i + isub;
            oadd[0] = j + jsub;
            OperatorAddLaderNormalOrder(H, oadd, 2, HVRe, HVIm, &QDetNotOrderedFcFaBcBa);
          }
        }
      }
    }
    else
    {                           /* A.H[iH].Hop is Real */
      for (isub = 0; isub < A.H[iH].Hop.Ni; isub++)
      {
        for (jsub = 0; jsub < A.H[iH].Hop.Nj; jsub++)
        {
          if (fabs(A.H[iH].Hop.VRe[isub][jsub]) > DBL_MIN * 10)
          {
            HVRe = A.H[iH].Hop.VRe[isub][jsub] * alphaRe;
            HVIm = A.H[iH].Hop.VRe[isub][jsub] * alphaIm;
            oadd[1] = CrAnMask + i + isub;
            oadd[0] = j + jsub;
            OperatorAddLaderNormalOrder(H, oadd, 2, HVRe, HVIm, &QDetNotOrderedFcFaBcBa);
          }
        }
      }
    }
  }
  OperatorChop(H, DBL_MIN*10);
  return 0;
}
int CreateHamiltonianOperatorFromTightBindingDefPlusSpin(TightBindingDef A, double k[3], OperatorType * H)
{
  int    iH, i, j, isub, jsub;
  double alphaRe, alphaIm, HVRe, HVIm;
  double vector[3];
  uint16_t oadd[2];
  InitOperator(H);
  H->NFermion = 2*A.NFermion;
  H->NWord = Min(H->NFermion, 1 + (H->NFermion - 1) / 8);
  InitOperatorFixedLength(H, 2, 10);
  for (iH = 0; iH < A.nH; iH++)
  {
    i = A.FermionIndex[A.H[iH].Atompos1][A.H[iH].Shellpos1];
    j = A.FermionIndex[A.H[iH].Atompos2][A.H[iH].Shellpos2];
    vector[0] = A.Cell.Atom[A.H[iH].Atompos1].Position[0] - A.Cell.Atom[A.H[iH].Atompos2].Position[0];
    vector[1] = A.Cell.Atom[A.H[iH].Atompos1].Position[1] - A.Cell.Atom[A.H[iH].Atompos2].Position[1];
    vector[2] = A.Cell.Atom[A.H[iH].Atompos1].Position[2] - A.Cell.Atom[A.H[iH].Atompos2].Position[2];
    vector[0] += A.Cell.rUnitCell[0][0] * A.H[iH].TransVecs[0] + A.Cell.rUnitCell[1][0] * A.H[iH].TransVecs[1] + A.Cell.rUnitCell[2][0] * A.H[iH].TransVecs[2];
    vector[1] += A.Cell.rUnitCell[0][1] * A.H[iH].TransVecs[0] + A.Cell.rUnitCell[1][1] * A.H[iH].TransVecs[1] + A.Cell.rUnitCell[2][1] * A.H[iH].TransVecs[2];
    vector[2] += A.Cell.rUnitCell[0][2] * A.H[iH].TransVecs[0] + A.Cell.rUnitCell[1][2] * A.H[iH].TransVecs[1] + A.Cell.rUnitCell[2][2] * A.H[iH].TransVecs[2];
    alphaRe = cos(Dot(k, vector));
    alphaIm = sin(Dot(k, vector));
    if (Sqr(alphaIm) < Sqr(alphaRe) * DBL_EPSILON * 100) alphaIm = 0;
    if (A.H[iH].Hop.QComplex)
    {                           /* A.H[iH].Hop is Complex */
      for (isub = 0; isub < A.H[iH].Hop.Ni; isub++)
      {
        for (jsub = 0; jsub < A.H[iH].Hop.Nj; jsub++)
        {
          if ((fabs(A.H[iH].Hop.VRe[isub][jsub]) + fabs(A.H[iH].Hop.VIm[isub][jsub])) > DBL_MIN * 10)
          {
            HVRe = A.H[iH].Hop.VRe[isub][jsub] * alphaRe - A.H[iH].Hop.VIm[isub][jsub] * alphaIm;
            HVIm = A.H[iH].Hop.VRe[isub][jsub] * alphaIm + A.H[iH].Hop.VIm[isub][jsub] * alphaRe;
            oadd[1] = CrAnMask + 2*(i + isub);
            oadd[0] = 2*(j + jsub);
            OperatorAddLaderNormalOrder(H, oadd, 2, HVRe, HVIm, &QDetNotOrderedFcFaBcBa);
            oadd[1] = CrAnMask + 2*(i + isub)+1;
            oadd[0] = 2*(j + jsub)+1;
            OperatorAddLaderNormalOrder(H, oadd, 2, HVRe, HVIm, &QDetNotOrderedFcFaBcBa);
          }
        }
      }
    }
    else
    {                           /* A.H[iH].Hop is Real */
      for (isub = 0; isub < A.H[iH].Hop.Ni; isub++)
      {
        for (jsub = 0; jsub < A.H[iH].Hop.Nj; jsub++)
        {
          if (fabs(A.H[iH].Hop.VRe[isub][jsub]) > DBL_MIN * 10)
          {
            HVRe = A.H[iH].Hop.VRe[isub][jsub] * alphaRe;
            HVIm = A.H[iH].Hop.VRe[isub][jsub] * alphaIm;
            oadd[1] = CrAnMask + 2*(i + isub);
            oadd[0] = 2*(j + jsub);
            OperatorAddLaderNormalOrder(H, oadd, 2, HVRe, HVIm, &QDetNotOrderedFcFaBcBa);
            oadd[1] = CrAnMask + 2*(i + isub)+1;
            oadd[0] = 2*(j + jsub)+1;
            OperatorAddLaderNormalOrder(H, oadd, 2, HVRe, HVIm, &QDetNotOrderedFcFaBcBa);
          }
        }
      }
    }
  }
  OperatorChop(H, DBL_MIN*10);
  return 0;
}

/*
 *
 */
uint32_t CalculateTotalEnergyFromListOfEigenvalues(const double * const val, const uint32_t NVal, const double kbT, const double mu, double *Etot)
{
  uint32_t i;
  *Etot=0;
  for(i=0;i<NVal;i++)
  {
    *Etot+=FermiFunction(val[i], mu, kbT)*val[i];
  }
  return 0;
}


/*
 *
 */
uint32_t CalculateChemicalPotentialFromListOfEigenvalues(const double * const val, const uint32_t NVal, const double kbT, const double NE, double * const mu, double *Etot)
{
  double musmaler, mularger;
  double NEmu,NEmusmaler=0,NEmularger=0;
  Bool Converged=False;
  uint32_t i,iloop=0;
  musmaler=val[0];
  mularger=val[0];
  for(i=0;i<NVal;i++)
  {
    if(val[i]<musmaler)musmaler=val[i];
    if(val[i]>mularger)mularger=val[i];
  }
  for(i=0;i<NVal;i++)
  {
    NEmularger+=FermiFunction(val[i], mularger, kbT);
    NEmusmaler+=FermiFunction(val[i], musmaler, kbT);
  }
  while(!Converged)
  {
    iloop++;
    if(NEmularger<NE) mularger+=NE-NEmularger;
    if(NEmusmaler>NE) musmaler+=NE-NEmusmaler;
    *mu=(musmaler+mularger)/2.0;
    NEmu=0;
    for(i=0;i<NVal;i++)
    {
      NEmu+=FermiFunction(val[i], *mu, kbT);
    }
    if(Abs(NEmu-NE)<100*DBL_EPSILON)
    {
      Converged=True;
    }
    if(iloop==64)
    {
      printf("Warning: CalculateChemicalPotentialFromListOfEigenvalues did not converge in %i steps musmaller=%22.15E mularger=%22.15E\nOccupation is %22.15E wheras it should be %22.15E\n",iloop,musmaler,mularger,NEmu,NE);
      CalculateTotalEnergyFromListOfEigenvalues(val, NVal, kbT, *mu, Etot);
      return 1;
    }
    if(NEmu>NE)
    {
      mularger=*mu;
      NEmularger=NEmu;
    }
    else
    {
      musmaler=*mu;
      NEmusmaler=NEmu;
    }
  }
  CalculateTotalEnergyFromListOfEigenvalues(val, NVal, kbT, *mu, Etot);
  return 0;
}


/*  Calculates the fermi surface
 *  the k points to use are defined by the parallelogram defined by the left bottom point k0, left top point k1 and right bottom point k2
 *  intput:
 *     TightBindingDef A
 *     double Emin          : lowest energy point for which the fermi function is calculated
 *     double Emax          : highest energy point for which the fermi function is calculated
 *     int NE               : number of energy points used at wich the fermi function is calculated
 *     double k0            : left bottom of the parallelogram at wich the fermi function is calcualted
 *     double k1            : left top of the parallelogram at wich the fermi function is calcualted
 *     double k2            : right bottom of the parallelogram at wich the fermi function is calcualted
 *     Fermi[Ei].Ni         : number of k points on the axes k0 to k1
 *     Fermi[Ei].Nj         : number of k points on the axes k0 to k2
 *     Fermi is a matrix with NE+1 components.
 *  output:
 *     Fermi[Ei].Re[ki][kj] is incremented by the DOS at the given k-points and energies.
 */
int TightBindingFermi(TightBindingDef A, double Emin, double Emax, int NE, double k0[3], double k1[3], double k2[3],
                      MatrixType Fermi[])
{
  MatrixType H, T;
  double q[3], val[A.NFermion];
  const double dE = (Emax - Emin) / ((double)NE);
  int    ik1, ik2, ifermion, iE;
  for (ik1 = 0; ik1 < Fermi[0].Ni; ik1++)
  {
    for (ik2 = 0; ik2 < Fermi[0].Nj; ik2++)
    {
      q[0] = k0[0] + ik1 * k1[0] / ((double)(Fermi[0].Ni - 1)) + ik2 * k2[0] / ((double)(Fermi[0].Nj - 1));
      q[1] = k0[1] + ik1 * k1[1] / ((double)(Fermi[0].Ni - 1)) + ik2 * k2[1] / ((double)(Fermi[0].Nj - 1));
      q[2] = k0[2] + ik1 * k1[2] / ((double)(Fermi[0].Ni - 1)) + ik2 * k2[2] / ((double)(Fermi[0].Nj - 1));
      CreateHamiltonianFromTightBindingDef(A, q, &H);
      Eigensystem(H, &T, val);
      FreeMatrix(&H);
      FreeMatrix(&T);
      for (ifermion = 0; ifermion < A.NFermion; ifermion++)
      {
        iE = (int)(floor((val[ifermion] - Emin) / dE));
        if (iE >= NE)
        {
          Fermi[NE].VRe[ik1][ik2] += 1;
        }
        else if (iE < 0)
        {
          Fermi[0].VRe[ik1][ik2] += 1;
        }
        else
        {
          Fermi[iE].VRe[ik1][ik2] += 1 + iE - (val[ifermion] - Emin) / dE;
          Fermi[iE + 1].VRe[ik1][ik2] -= iE - (val[ifermion] - Emin) / dE;
        }
      }
    }
  }
  return 0;
}



int TightBindingDOS(TightBindingDef A, int Nk[3], SpectraType * DOS)
{
  MatrixType H, T;
  double q[3], val[A.NFermion];
  int    iq[3], ifermion, iE, jE;
  double diq[3],rest;
  for (iE = 0; iE <= DOS->N; iE++)
  {
    DOS->Im[iE]=0;
    DOS->Re[iE]=0;
  }
  for (iq[0] = 0; iq[0] < Nk[0]; iq[0]++)
  {
    for (iq[1] = 0; iq[1] < Nk[1]; iq[1]++)
    {
      for (iq[2] = 0; iq[2] < Nk[2]; iq[2]++)
      {
        diq[0] = 1.0 * (iq[0] + 0.0) / ((double)(Nk[0]));
        diq[1] = 1.0 * (iq[1] + 0.0) / ((double)(Nk[1]));
        diq[2] = 1.0 * (iq[2] + 0.0) / ((double)(Nk[2]));
        q[0] = diq[0] * A.Cell.gUnitCell[0][0] + diq[1] * A.Cell.gUnitCell[1][0] + diq[2] * A.Cell.gUnitCell[2][0];
        q[1] = diq[0] * A.Cell.gUnitCell[0][1] + diq[1] * A.Cell.gUnitCell[1][1] + diq[2] * A.Cell.gUnitCell[2][1];
        q[2] = diq[0] * A.Cell.gUnitCell[0][2] + diq[1] * A.Cell.gUnitCell[1][2] + diq[2] * A.Cell.gUnitCell[2][2];
        CreateHamiltonianFromTightBindingDef(A, q, &H);
        Eigensystem(H, &T, val);
        FreeMatrix(&H);
        FreeMatrix(&T);
        for (ifermion = 0; ifermion < A.NFermion; ifermion++)
        {
          rest = (val[ifermion] - DOS->Emin) / DOS->dE;
          iE = (int)(floor(rest));
          if (iE >= DOS->N)
          {
            DOS->Im[DOS->N] -= 1;
          }
          else if (iE < 0)
          {
            DOS->Im[0] -= 1;
          }
          else
          {
            rest -= iE;
            DOS->Im[iE] -= 1 - rest;
            DOS->Im[iE + 1] -= rest;
          }
        }
      }
    }
  }
  for (iE = 0; iE <= DOS->N; iE++)
  {
    DOS->Im[iE] /= DOS->dE*(Nk[0] * Nk[1] * Nk[2]);
  }
  for (iE = 0; iE <= DOS->N; iE++)
  {
    for(jE=0;jE<=DOS->N;jE++)
    {
      if(jE!=iE)
      {
        DOS->Re[iE]+=DOS->Im[jE]/((double)(jE-iE));
      }
    }
    DOS->Re[iE]/=Pi;
  }
  return 0;
}


/*
 *
 */
int TightBindingBANDS(TightBindingDef A, BandstructurePathType kList, MatrixType * BANDS)
{
  int    iq, iqpath, Nqtot, iqtot, ival;
  double q[3], val[A.NFermion];
  MatrixType H, T;
  Nqtot = 0;
  for (iqpath = 0; iqpath < kList.NQ; iqpath++)
  {
    Nqtot += kList.Q[iqpath].Nq;
  }
  BANDS->Ni = Nqtot;
  BANDS->Nj = A.NFermion;
  BANDS->QComplex = 0;
  InitMatrix(BANDS);
  iqtot = 0;
  for (iqpath = 0; iqpath < kList.NQ; iqpath++)
  {
    for (iq = 0; iq < kList.Q[iqpath].Nq; iq++)
    {
      if (iqpath == kList.NQ-1)
      {
        q[0] = kList.Q[iqpath].q[0];
        q[1] = kList.Q[iqpath].q[1];
        q[2] = kList.Q[iqpath].q[2];
      }
      else
      {
        q[0] =
          (kList.Q[iqpath].q[0] * (kList.Q[iqpath].Nq - iq) +
           kList.Q[iqpath + 1].q[0] * (iq)) / ((double)kList.Q[iqpath].Nq);
        q[1] =
          (kList.Q[iqpath].q[1] * (kList.Q[iqpath].Nq - iq) +
           kList.Q[iqpath + 1].q[1] * (iq)) / ((double)kList.Q[iqpath].Nq);
        q[2] =
          (kList.Q[iqpath].q[2] * (kList.Q[iqpath].Nq - iq) +
           kList.Q[iqpath + 1].q[2] * (iq)) / ((double)kList.Q[iqpath].Nq);
      }

      CreateHamiltonianFromTightBindingDef(A, q, &H);
      Eigensystem(H, &T, val);
      for (ival = 0; ival < A.NFermion; ival++)
      {
        BANDS->VRe[iqtot][ival] = val[ival];
      }
      FreeMatrix(&T);
      FreeMatrix(&H);
      iqtot += 1;
    }
  }
  return 0;
}


/*  given an eigen function the character matrix is calculated as char_{i,j} = psi_i * psi_j^*
 *
 */
void RealEigenFunctionToCharacter(double const * EigenFun, MatrixType *Char)
{
  uint32_t i,j;
  for(i=0;i<Char->Ni;i++)
  {
    for(j=0;j<Char->Nj;j++)
    {
      Char->VRe[i][j] = EigenFun[i]*EigenFun[j];
    }
  }
}
void RealEigenFunctionToCompactCharacter(double const * EigenFun, CompactMatrixType *Char)
{
  uint32_t i,j,k=0;
  for(i=0;i<Char->Ni;i++)
  {
    for(j=0;j<Char->Nj;j++)
    {
      Char->VRe[k] = EigenFun[i]*EigenFun[j];
      k++;
    }
  }
}

/*  given an eigen function the character matrix is calculated as char_{i,j} = psi_i * psi_j^*
 *
 */
void ComplexEigenFunctionToCharacter(Complex16 const * EigenFun, MatrixType *Char)
{
  uint32_t i,j;
  for(i=0;i<Char->Ni;i++)
  {
    for(j=0;j<Char->Nj;j++)
    {
      Char->VRe[i][j] = EigenFun[i].Re*EigenFun[j].Re+EigenFun[i].Im*EigenFun[j].Im;
      Char->VIm[i][j] = EigenFun[i].Im*EigenFun[j].Re-EigenFun[i].Re*EigenFun[j].Im;
    }
  }
}
void ComplexEigenFunctionToCompactCharacter(Complex16 const * EigenFun, CompactMatrixType *Char)
{
  uint32_t i,j,k=0;
  for(i=0;i<Char->Ni;i++)
  {
    for(j=0;j<Char->Nj;j++)
    {
      Char->V[k].Re = EigenFun[i].Re*EigenFun[j].Re+EigenFun[i].Im*EigenFun[j].Im;
      Char->V[k].Im = EigenFun[i].Im*EigenFun[j].Re-EigenFun[i].Re*EigenFun[j].Im;
      k++;
    }
  }
}

/*  given an eigen function the character matrix is calculated as char_{i,j} = psi_i * psi_j
 *
 */
void EigenFunctionToCharacter(void const * EigenFun, MatrixType *Char)
{
  if(Char->QComplex)
  {
    ComplexEigenFunctionToCharacter((Complex16*) EigenFun, Char);
  }
  else
  {
    RealEigenFunctionToCharacter((double*) EigenFun, Char);
  }
}
void EigenFunctionToCompactCharacter(void const * EigenFun, CompactMatrixType *Char)
{
  if(Char->QComplex)
  {
    ComplexEigenFunctionToCompactCharacter((Complex16*) EigenFun, Char);
  }
  else
  {
    RealEigenFunctionToCompactCharacter((double*) EigenFun, Char);
  }
}


/*  given an eigen function the character matrix is calculated as char_{i,j} = psi_i * psi_j
 *
 */
void RealEigenFunctionToCharacterMultiSite(double const * EigenFun, MatrixType *Char, uint32_t const Nvec, int const *vec)
{
  uint32_t i,j;
  for(i=0;i<Char->Ni;i++)
  {
    for(j=0;j<Char->Nj;j++)
    {
      Char->VRe[i][j] = EigenFun[i]*EigenFun[j];
    }
  }
}

/*  given an eigen function the character matrix is calculated as char_{i,j} = psi_i * psi_j
 *
 */
void ComplexEigenFunctionToCharacterMultiSite(Complex16 const * EigenFun, MatrixType *Char, uint32_t const Nvec, int const *vec)
{
  uint32_t i,j;
  for(i=0;i<Char->Ni;i++)
  {
    for(j=0;j<Char->Nj;j++)
    {
      Char->VRe[i][j] = EigenFun[i].Re*EigenFun[j].Re+EigenFun[i].Im*EigenFun[j].Im;
      Char->VIm[i][j] = EigenFun[i].Im*EigenFun[j].Re-EigenFun[i].Re*EigenFun[j].Im;
    }
  }
}

/*  given an eigen function the character matrix is calculated as char_{i,j} = psi_i * psi_j
 *
 */
void EigenFunctionToCharacterMultiSite(void const * EigenFun, MatrixType *Char, uint32_t const Nvec, int const *vec)
{
  if(Char->QComplex)
  {
    ComplexEigenFunctionToCharacterMultiSite((Complex16*) EigenFun, Char, Nvec, vec);
  }
  else
  {
    RealEigenFunctionToCharacterMultiSite((double*) EigenFun, Char, Nvec, vec);
  }
}


int EigenFunctionToCharactersStandardSetNCharacter(TightBindingDef A, int *NCharacter)
{
  *NCharacter = A.NFermion + 1;
  return 0;
}

int EigenFunctionToCharactersStandard(TightBindingDef A, MatrixType H, double *extra, MatrixType EigenFun,
                                      MatrixType * Character)
{
  int    i, j;
  Character->Ni = A.NFermion;
  Character->Nj = A.NFermion + 1;
  Character->QComplex = 0;
  InitMatrix(Character);
  if (EigenFun.QComplex)
  {
    double CharImTemp;
    for (i = 0; i < Character->Ni; i++)
    {
      Character->VRe[i][Character->Nj - 1] = 0;
      CharImTemp = 0;
      for (j = 0; j < Character->Nj - 1; j++)
      {
        Character->VRe[i][j] = Sqr(EigenFun.VRe[i][j]) + Sqr(EigenFun.VIm[i][j]);
        Character->VRe[i][Character->Nj - 1] += EigenFun.VRe[i][j];
        CharImTemp += EigenFun.VIm[i][j];
      }
      Character->VRe[i][Character->Nj - 1] *= Character->VRe[i][Character->Nj - 1];
      Character->VRe[i][Character->Nj - 1] += Sqr(CharImTemp);
    }
  }
  else
  {
    for (i = 0; i < Character->Ni; i++)
    {
      Character->VRe[i][Character->Nj - 1] = 0;
      for (j = 0; j < Character->Nj - 1; j++)
      {
        Character->VRe[i][j] = Sqr(EigenFun.VRe[i][j]);
        Character->VRe[i][Character->Nj - 1] += EigenFun.VRe[i][j];
      }
      Character->VRe[i][Character->Nj - 1] *= Character->VRe[i][Character->Nj - 1];
    }
  }
  return 0;
}


int EigenFunctionToCharactersOnsiteEnergySetNCharacter(TightBindingDef A, int *NCharacter)
{
  *NCharacter = 3;
  return 0;
}

int EigenFunctionToCharactersOnsiteEnergy(TightBindingDef A, MatrixType H, double *extra, MatrixType EigenFun,
                                          MatrixType * Character)
{
  int    i, j;
  Character->Ni = A.NFermion;
  Character->Nj = 3;
  Character->QComplex = 0;
  InitMatrix(Character);
  if (EigenFun.QComplex)
  {
    for (i = 0; i < Character->Ni; i++)
    {
      Character->VRe[i][0] = 1.0 / A.NFermion;
      Character->VRe[i][1] = 0;
      Character->VRe[i][2] = 0;
      for (j = 0; j < A.NFermion; j++)
      {
        if (H.VRe[j][j] < 0)
        {
          Character->VRe[i][1] += (Sqr(EigenFun.VRe[i][j]) + Sqr(EigenFun.VIm[i][j])) / A.NFermion;
        }
        else
        {
          Character->VRe[i][2] += (Sqr(EigenFun.VRe[i][j]) + Sqr(EigenFun.VIm[i][j])) / A.NFermion;
        }
      }
    }
  }
  else
  {
    for (i = 0; i < Character->Ni; i++)
    {
      Character->VRe[i][0] = 1.0 / A.NFermion;
      Character->VRe[i][1] = 0;
      Character->VRe[i][2] = 0;
      for (j = 0; j < A.NFermion; j++)
      {
        if (H.VRe[j][j] < 0)
        {
          Character->VRe[i][0] += Sqr(EigenFun.VRe[i][j]) / A.NFermion;
        }
        else
        {
          Character->VRe[i][1] += Sqr(EigenFun.VRe[i][j]) / A.NFermion;
        }
      }
    }
  }
  return 0;
}

int EigenFunctionToCharactersOnsiteEnergyPhaseSetNCharacter(TightBindingDef A, int *NCharacter)
{
  *NCharacter = 3;
  return 0;
}

int EigenFunctionToCharactersOnsiteEnergyPhase(TightBindingDef A, MatrixType H, double *extra, MatrixType EigenFun,
                                               MatrixType * Character)
{
  int    i, j;
  Character->Ni = A.NFermion;
  Character->Nj = 3;
  Character->QComplex = 0;
  InitMatrix(Character);
  if (EigenFun.QComplex)
  {
    double CharImTemp;
    for (i = 0; i < Character->Ni; i++)
    {
      Character->VRe[i][0] = 0;
      Character->VRe[i][1] = 0;
      Character->VRe[i][2] = 0;
      CharImTemp = 0;
      for (j = 0; j < A.NFermion; j++)
      {
        if (H.VRe[j][j] < 0)
        {
          Character->VRe[i][0] += Sqr(EigenFun.VRe[i][j]) + Sqr(EigenFun.VIm[i][j]);
        }
        else
        {
          Character->VRe[i][1] += Sqr(EigenFun.VRe[i][j]) + Sqr(EigenFun.VIm[i][j]);
        }
        Character->VRe[i][2] += EigenFun.VRe[i][j];
        CharImTemp += EigenFun.VIm[i][j];
      }
      Character->VRe[i][2] *= Character->VRe[i][2];
      Character->VRe[i][2] += Sqr(CharImTemp);

      Character->VRe[i][2] /= (double)A.NFermion;
      if (Character->VRe[i][2] < DBL_EPSILON)
        Character->VRe[i][2] = 0;
      Character->VRe[i][0] *= Character->VRe[i][2];
      Character->VRe[i][1] *= Character->VRe[i][2];
    }
  }
  else
  {
    for (i = 0; i < Character->Ni; i++)
    {
      Character->VRe[i][0] = 0;
      Character->VRe[i][1] = 0;
      Character->VRe[i][2] = 0;
      for (j = 0; j < A.NFermion; j++)
      {
        if (H.VRe[j][j] < 0)
        {
          Character->VRe[i][0] += Sqr(EigenFun.VRe[i][j]);
        }
        else
        {
          Character->VRe[i][1] += Sqr(EigenFun.VRe[i][j]);
        }
        Character->VRe[i][2] += EigenFun.VRe[i][j];
      }
      Character->VRe[i][2] *= Character->VRe[i][2];
      Character->VRe[i][2] /= (double)A.NFermion;
      if (Character->VRe[i][2] < DBL_EPSILON)
        Character->VRe[i][2] = 0;
      Character->VRe[i][0] *= Character->VRe[i][2];
      Character->VRe[i][1] *= Character->VRe[i][2];
    }
  }
  return 0;
}


int EigenFunctionToCharactersOnsiteEnergyPhaseRGBSetNCharacter(TightBindingDef A, int *NCharacter)
{
  *NCharacter = 3;
  return 0;
}

int EigenFunctionToCharactersOnsiteEnergyPhaseRGB(TightBindingDef A, MatrixType H, double *extra, MatrixType EigenFun,
                                                  MatrixType * Character)
{
  int    i, j;
  Character->Ni = A.NFermion;
  Character->Nj = 3;
  Character->QComplex = 0;
  InitMatrix(Character);
  if (EigenFun.QComplex)
  {
    double CharImTemp;
    for (i = 0; i < Character->Ni; i++)
    {
      Character->VRe[i][0] = 0;
      Character->VRe[i][1] = 0;
      Character->VRe[i][2] = 0;
      CharImTemp = 0;
      for (j = 0; j < A.NFermion; j++)
      {
        if (H.VRe[j][j] < 0)
        {
          Character->VRe[i][0] += Sqr(EigenFun.VRe[i][j]) + Sqr(EigenFun.VIm[i][j]);
        }
        else
        {
          Character->VRe[i][2] += Sqr(EigenFun.VRe[i][j]) + Sqr(EigenFun.VIm[i][j]);
        }
        Character->VRe[i][1] += EigenFun.VRe[i][j];
        CharImTemp += EigenFun.VIm[i][j];
      }
      Character->VRe[i][1] *= Character->VRe[i][1];
      Character->VRe[i][1] += Sqr(CharImTemp);

      Character->VRe[i][1] /= (double)A.NFermion;
      if (Character->VRe[i][1] < DBL_EPSILON)
        Character->VRe[i][1] = 0;
      Character->VRe[i][0] *= Character->VRe[i][1];
      Character->VRe[i][2] *= Character->VRe[i][1];
    }
  }
  else
  {
    for (i = 0; i < Character->Ni; i++)
    {
      Character->VRe[i][0] = 0;
      Character->VRe[i][1] = 0;
      Character->VRe[i][2] = 0;
      for (j = 0; j < A.NFermion; j++)
      {
        if (H.VRe[j][j] < 0)
        {
          Character->VRe[i][0] += Sqr(EigenFun.VRe[i][j]);
        }
        else
        {
          Character->VRe[i][2] += Sqr(EigenFun.VRe[i][j]);
        }
        Character->VRe[i][1] += EigenFun.VRe[i][j];
      }
      Character->VRe[i][1] *= Character->VRe[i][1];
      Character->VRe[i][1] /= (double)A.NFermion;
      if (Character->VRe[i][1] < DBL_EPSILON)
        Character->VRe[i][1] = 0;
      Character->VRe[i][0] *= Character->VRe[i][1];
      Character->VRe[i][2] *= Character->VRe[i][1];
    }
  }
  return 0;
}


int EigenFunctionToCharactersOnsiteEnergyPhaseRGBkshiftedSetNCharacter(TightBindingDef A, int *NCharacter)
{
  *NCharacter = 3;
  return 0;
}

int EigenFunctionToCharactersOnsiteEnergyPhaseRGBkshifted(TightBindingDef A, MatrixType H, double kshift[3],
                                                          MatrixType EigenFun, MatrixType * Character)
{
  int    i, j;
  double CharImTemp;
  Character->Ni = A.NFermion;
  Character->Nj = 3;
  Character->QComplex = 0;
  InitMatrix(Character);
  if (EigenFun.QComplex)
  {
    for (i = 0; i < Character->Ni; i++)
    {
      Character->VRe[i][0] = 0;
      Character->VRe[i][1] = 0;
      Character->VRe[i][2] = 0;
      CharImTemp = 0;
      for (j = 0; j < A.NFermion; j++)
      {
        if (H.VRe[j][j] < 0)
        {
          Character->VRe[i][0] += Sqr(EigenFun.VRe[i][j]) + Sqr(EigenFun.VIm[i][j]);
        }
        else
        {
          Character->VRe[i][2] += Sqr(EigenFun.VRe[i][j]) + Sqr(EigenFun.VIm[i][j]);
        }
        Character->VRe[i][1] +=
          EigenFun.VRe[i][j] * cos(Dot(A.Cell.Atom[j].Position, kshift)) +
          EigenFun.VIm[i][j] * sin(Dot(A.Cell.Atom[j].Position, kshift));
        CharImTemp +=
          EigenFun.VIm[i][j] * cos(Dot(A.Cell.Atom[j].Position, kshift)) -
          EigenFun.VRe[i][j] * sin(Dot(A.Cell.Atom[j].Position, kshift));
      }
      Character->VRe[i][1] *= Character->VRe[i][1];
      Character->VRe[i][1] += Sqr(CharImTemp);
      Character->VRe[i][1] /= (double)A.NFermion;
      if (Character->VRe[i][1] < DBL_EPSILON)
        Character->VRe[i][1] = 0;
      Character->VRe[i][0] *= Character->VRe[i][1];
      Character->VRe[i][2] *= Character->VRe[i][1];
    }
  }
  else
  {
    for (i = 0; i < Character->Ni; i++)
    {
      Character->VRe[i][0] = 0;
      Character->VRe[i][1] = 0;
      Character->VRe[i][2] = 0;
      CharImTemp = 0;
      for (j = 0; j < A.NFermion; j++)
      {
        if (H.VRe[j][j] < 0)
        {
          Character->VRe[i][0] += Sqr(EigenFun.VRe[i][j]);
        }
        else
        {
          Character->VRe[i][2] += Sqr(EigenFun.VRe[i][j]);
        }
        Character->VRe[i][1] += EigenFun.VRe[i][j] * cos(Dot(A.Cell.Atom[j].Position, kshift));
        CharImTemp -= EigenFun.VRe[i][j] * sin(Dot(A.Cell.Atom[j].Position, kshift));
      }
      Character->VRe[i][1] *= Character->VRe[i][1];
      Character->VRe[i][1] += Sqr(CharImTemp);
      Character->VRe[i][1] /= (double)A.NFermion;
      if (Character->VRe[i][1] < DBL_EPSILON)
        Character->VRe[i][1] = 0;
      Character->VRe[i][0] *= Character->VRe[i][1];
      Character->VRe[i][2] *= Character->VRe[i][1];
    }
  }
  return 0;
}





/*  Calculates the fermi surface
 *  the k points to use are defined by the parallelogram defined by the left bottom point k0, left top point k1 and right bottom point k2
 *  intput:
 *     TightBindingDef A
 *     double Emin          : lowest energy point for which the fermi function is calculated
 *     double Emax          : highest energy point for which the fermi function is calculated
 *     int NE               : number of energy points used at wich the fermi function is calculated
 *     double k0            : left bottom of the parallelogram at wich the fermi function is calcualted
 *     double k1            : left top of the parallelogram at wich the fermi function is calcualted
 *     double k2            : right bottom of the parallelogram at wich the fermi function is calcualted
 *     Fermi[Ei].Ni         : number of k points on the axes k0 to k1
 *     Fermi[Ei].Nj         : number of k points on the axes k0 to k2
 *     Fermi is a matrix with Ncharacter*(NE+1) components.
 *     int Ncharacter       : number of different characters calculated
 *     int *EigenFunctionToCharacters  : a function that specifies how to go from a wave function to a list of characters
 *  output:
 *     Fermi[Ncharacter*Ei+icharacter].Re[ki][kj] is incremented by the DOS at the given k-points and energies.
 */
int TightBindingFermiCharacter(TightBindingDef A, double Emin, double Emax, int NE, double k0[3], double k1[3],
                               double k2[3], MatrixType Fermi[], int Ncharacter,
                               int (*EigenFunctionToCharacters) (TightBindingDef, MatrixType, double *, MatrixType,
                                                                 MatrixType *))
{
  MatrixType H, T, TempChar;
  double q[3], val[A.NFermion];
  const double dE = (Emax - Emin) / ((double)NE);
  int    ik1, ik2, ifermion, iE, icharacter;
  for (ik1 = 0; ik1 < Fermi[0].Ni; ik1++)
  {
    for (ik2 = 0; ik2 < Fermi[0].Nj; ik2++)
    {
      q[0] = k0[0] + ik1 * k1[0] / ((double)(Fermi[0].Ni - 1)) + ik2 * k2[0] / ((double)(Fermi[0].Nj - 1));
      q[1] = k0[1] + ik1 * k1[1] / ((double)(Fermi[0].Ni - 1)) + ik2 * k2[1] / ((double)(Fermi[0].Nj - 1));
      q[2] = k0[2] + ik1 * k1[2] / ((double)(Fermi[0].Ni - 1)) + ik2 * k2[2] / ((double)(Fermi[0].Nj - 1));
      CreateHamiltonianFromTightBindingDef(A, q, &H);
      Eigensystem(H, &T, val);
      EigenFunctionToCharacters(A, H, NULL, T, &TempChar);
      FreeMatrix(&H);
      FreeMatrix(&T);
      for (icharacter = 0; icharacter < Ncharacter; icharacter++)
      {
        for (ifermion = 0; ifermion < A.NFermion; ifermion++)
        {
          iE = (int)(floor((val[ifermion] - Emin) / dE));
          if (iE >= NE - 1)
          {
            Fermi[Ncharacter * NE + icharacter].VRe[ik1][ik2] += 2 * TempChar.VRe[ifermion][icharacter];
          }
          else if (iE <= 0)
          {
            Fermi[icharacter].VRe[ik1][ik2] += 2 * TempChar.VRe[ifermion][icharacter];
          }
          else
          {
            Fermi[Ncharacter * (iE - 1) + icharacter].VRe[ik1][ik2] +=
              0.25 * TempChar.VRe[ifermion][icharacter] * (1 + iE - (val[ifermion] - Emin) / dE);
            Fermi[Ncharacter * (iE + 0) + icharacter].VRe[ik1][ik2] -=
              0.25 * TempChar.VRe[ifermion][icharacter] * (iE - (val[ifermion] - Emin) / dE);

            Fermi[Ncharacter * (iE + 0) + icharacter].VRe[ik1][ik2] +=
              0.50 * TempChar.VRe[ifermion][icharacter] * (1 + iE - (val[ifermion] - Emin) / dE);
            Fermi[Ncharacter * (iE + 1) + icharacter].VRe[ik1][ik2] -=
              0.50 * TempChar.VRe[ifermion][icharacter] * (iE - (val[ifermion] - Emin) / dE);

            Fermi[Ncharacter * (iE + 1) + icharacter].VRe[ik1][ik2] +=
              0.25 * TempChar.VRe[ifermion][icharacter] * (1 + iE - (val[ifermion] - Emin) / dE);
            Fermi[Ncharacter * (iE + 2) + icharacter].VRe[ik1][ik2] -=
              0.25 * TempChar.VRe[ifermion][icharacter] * (iE - (val[ifermion] - Emin) / dE);
          }
        }
      }
      FreeMatrix(&TempChar);
    }
  }
  return 0;
}


/*  Calculates the fermi surface
 *  the k points to use are defined by the parallelogram defined by the left bottom point k0, left top point k1 and right bottom point k2
 *
 *  It is assumed that TightBindingDef describes a supercell structure.
 *
 *  intput:
 *     TightBindingDef A
 *     double Emin          : lowest energy point for which the fermi function is calculated
 *     double Emax          : highest energy point for which the fermi function is calculated
 *     int NE               : number of energy points used at wich the fermi function is calculated
 *     double k0            : left bottom of the parallelogram at wich the fermi function is calcualted
 *     double k1            : left top of the parallelogram at wich the fermi function is calcualted
 *     double k2            : right bottom of the parallelogram at wich the fermi function is calcualted
 *     int Nk1              :
 *     int Nk2              :
 *     Fermi[Ei].Ni         : number of k points on the axes k0 to k1
 *     Fermi[Ei].Nj         : number of k points on the axes k0 to k2
 *     Fermi is a matrix with Ncharacter*(NE+1) components.
 *     int Ncharacter       : number of different characters calculated
 *     int *EigenFunctionToCharacters  : a function that specifies how to go from a wave function to a list of characters
 *  output:
 *     Fermi[Ncharacter*Ei+icharacter].Re[ki][kj] is incremented by the DOS at the given k-points and energies.
 */
int TightBindingFermiCharacterUnfolding(TightBindingDef A, double Emin, double Emax, int NE, double k0[3], double k1[3],
                                        double k2[3], int Nk1, int Nk2, MatrixType Fermi[], int Ncharacter,
                                        int (*EigenFunctionToCharacters) (TightBindingDef, MatrixType, double *,
                                                                          MatrixType, MatrixType *))
{
  MatrixType H, T, TempChar;
  double q[3], val[A.NFermion], kshift[3];
  const double dE = (Emax - Emin) / ((double)NE);
  int    ik1, ik2, ifermion, iE, icharacter, iksub1, iksub2;
  int    ik1max, ik2max;
  ik1max = (int)(floor(Fermi[0].Ni / Nk1));
  ik2max = (int)(floor(Fermi[0].Nj / Nk2));
  for (ik1 = 0; ik1 < ik1max; ik1++)
  {
    for (ik2 = 0; ik2 < ik2max; ik2++)
    {
      q[0] = k0[0] + ik1 * k1[0] / ((double)(ik1max)) + ik2 * k2[0] / ((double)(ik2max));
      q[1] = k0[1] + ik1 * k1[1] / ((double)(ik1max)) + ik2 * k2[1] / ((double)(ik2max));
      q[2] = k0[2] + ik1 * k1[2] / ((double)(ik1max)) + ik2 * k2[2] / ((double)(ik2max));
      CreateHamiltonianFromTightBindingDef(A, q, &H);
      Eigensystem(H, &T, val);
      for (iksub1 = 0; iksub1 < Nk1; iksub1++)
      {
        for (iksub2 = 0; iksub2 < Nk2; iksub2++)
        {
          kshift[0] = iksub1 * (k1[0]) + iksub2 * (k2[0]);
          kshift[1] = iksub1 * (k1[1]) + iksub2 * (k2[1]);
          kshift[2] = iksub1 * (k1[2]) + iksub2 * (k2[2]);
          EigenFunctionToCharacters(A, H, kshift, T, &TempChar);
          for (icharacter = 0; icharacter < Ncharacter; icharacter++)
          {
            for (ifermion = 0; ifermion < A.NFermion; ifermion++)
            {
              iE = (int)(floor((val[ifermion] - Emin) / dE));
              if (iE >= NE - 1)
              {
                Fermi[Ncharacter * NE + icharacter].VRe[ik1 + ik1max * iksub1][ik2 + ik2max * iksub2] +=
                  2 * TempChar.VRe[ifermion][icharacter];
              }
              else if (iE <= 0)
              {
                Fermi[icharacter].VRe[ik1 + ik1max * iksub1][ik2 + ik2max * iksub2] +=
                  2 * TempChar.VRe[ifermion][icharacter];
              }
              else
              {
                Fermi[Ncharacter * (iE - 1) + icharacter].VRe[ik1 + ik1max * iksub1][ik2 + ik2max * iksub2] +=
                  0.25 * TempChar.VRe[ifermion][icharacter] * (1 + iE - (val[ifermion] - Emin) / dE);
                Fermi[Ncharacter * (iE + 0) + icharacter].VRe[ik1 + ik1max * iksub1][ik2 + ik2max * iksub2] -=
                  0.25 * TempChar.VRe[ifermion][icharacter] * (iE - (val[ifermion] - Emin) / dE);

                Fermi[Ncharacter * (iE + 0) + icharacter].VRe[ik1 + ik1max * iksub1][ik2 + ik2max * iksub2] +=
                  0.50 * TempChar.VRe[ifermion][icharacter] * (1 + iE - (val[ifermion] - Emin) / dE);
                Fermi[Ncharacter * (iE + 1) + icharacter].VRe[ik1 + ik1max * iksub1][ik2 + ik2max * iksub2] -=
                  0.50 * TempChar.VRe[ifermion][icharacter] * (iE - (val[ifermion] - Emin) / dE);

                Fermi[Ncharacter * (iE + 1) + icharacter].VRe[ik1 + ik1max * iksub1][ik2 + ik2max * iksub2] +=
                  0.25 * TempChar.VRe[ifermion][icharacter] * (1 + iE - (val[ifermion] - Emin) / dE);
                Fermi[Ncharacter * (iE + 2) + icharacter].VRe[ik1 + ik1max * iksub1][ik2 + ik2max * iksub2] -=
                  0.25 * TempChar.VRe[ifermion][icharacter] * (iE - (val[ifermion] - Emin) / dE);
              }
            }
          }
          FreeMatrix(&TempChar);


        }
      }
      FreeMatrix(&H);
      FreeMatrix(&T);

    }
  }
  return 0;
}




int TightBindingDOSCharacter(TightBindingDef A, int Nk[3], SpectraType * DOS, int Ncharacter, int (*EigenFunctionToCharacters) (TightBindingDef, MatrixType, double *, MatrixType, MatrixType *))
{
  MatrixType H, T, TempChar;
  double q[3], val[A.NFermion];
  int    iq[3], ifermion, iE, ichar;
  double diq[3];
  for (iq[0] = 0; iq[0] < Nk[0]; iq[0]++)
  {
    for (iq[1] = 0; iq[1] < Nk[1]; iq[1]++)
    {
      for (iq[2] = 0; iq[2] < Nk[2]; iq[2]++)
      {
        diq[0] = 0.5 * (iq[0] + 0.5) / ((double)(Nk[0]));
        diq[1] = 0.5 * (iq[1] + 0.5) / ((double)(Nk[1]));
        diq[2] = 0.5 * (iq[2] + 0.5) / ((double)(Nk[2]));
        q[0] = diq[0] * A.Cell.gUnitCell[0][0] + diq[1] * A.Cell.gUnitCell[1][0] + diq[2] * A.Cell.gUnitCell[2][0];
        q[1] = diq[0] * A.Cell.gUnitCell[0][1] + diq[1] * A.Cell.gUnitCell[1][1] + diq[2] * A.Cell.gUnitCell[2][1];
        q[2] = diq[0] * A.Cell.gUnitCell[0][2] + diq[1] * A.Cell.gUnitCell[1][2] + diq[2] * A.Cell.gUnitCell[2][2];
        CreateHamiltonianFromTightBindingDef(A, q, &H);
        Eigensystem(H, &T, val);
        EigenFunctionToCharacters(A, H, NULL, T, &TempChar);
        FreeMatrix(&H);
        FreeMatrix(&T);
        for (ichar = 0; ichar < Ncharacter; ichar++)
        {
          for (ifermion = 0; ifermion < A.NFermion; ifermion++)
          {
            iE = (int)(floor((val[ifermion] - DOS[ichar].Emin) / DOS[ichar].dE));
            if (iE >= DOS[ichar].N)
            {
              DOS[ichar].Re[DOS[ichar].N] += TempChar.VRe[ifermion][ichar];
            }
            else if (iE < 0)
            {
              DOS[ichar].Re[0] += TempChar.VRe[ifermion][ichar];
            }
            else
            {
              DOS[ichar].Re[iE] +=
                TempChar.VRe[ifermion][ichar] * (1 + iE - (val[ifermion] - DOS[ichar].Emin) / DOS[ichar].dE);
              DOS[ichar].Re[iE + 1] -=
                TempChar.VRe[ifermion][ichar] * (iE - (val[ifermion] - DOS[ichar].Emin) / DOS[ichar].dE);
            }
          }
        }
        FreeMatrix(&TempChar);
      }
    }
  }
  for (ichar = 0; ichar < Ncharacter; ichar++)
  {
    for (iE = 0; iE <= DOS[ichar].N; iE++)
    {
      DOS[ichar].Re[iE] /= (Nk[0] * Nk[1] * Nk[2]);
    }
  }
  return 0;
}

int TightBindingBANDSCharacter(TightBindingDef A, BandstructurePathType kList, MatrixType * BANDS, int Ncharacter, MatrixType Character[], int (*EigenFunctionToCharacters) (TightBindingDef, MatrixType, double *, MatrixType, MatrixType *))
{
  int    iq, iqpath, Nqtot, iqtot, ival, ichar;
  double q[3], val[A.NFermion];
  MatrixType H, T, TempChar;
  Nqtot = 0;
  for (iqpath = 0; iqpath < kList.NQ; iqpath++)
  {
    Nqtot += kList.Q[iqpath].Nq;
  }
  BANDS->Ni = Nqtot;
  BANDS->Nj = A.NFermion;
  BANDS->QComplex = 0;
  InitMatrix(BANDS);
  iqtot = 0;
  for (ichar = 0; ichar < Ncharacter; ichar++)
  {
    Character[ichar].Ni = Nqtot;
    Character[ichar].Nj = A.NFermion;
    Character[ichar].QComplex = 0;
    InitMatrix(&Character[ichar]);
  }
  for (iqpath = 0; iqpath < kList.NQ; iqpath++)
  {
    for (iq = 0; iq < kList.Q[iqpath].Nq; iq++)
    {
      if (iqpath == kList.NQ-1)
      {
        q[0] = kList.Q[iqpath].q[0];
        q[1] = kList.Q[iqpath].q[1];
        q[2] = kList.Q[iqpath].q[2];
      }
      else
      {
        q[0] =
          (kList.Q[iqpath].q[0] * (kList.Q[iqpath].Nq - iq) +
           kList.Q[iqpath + 1].q[0] * (iq)) / ((double)kList.Q[iqpath].Nq);
        q[1] =
          (kList.Q[iqpath].q[1] * (kList.Q[iqpath].Nq - iq) +
           kList.Q[iqpath + 1].q[1] * (iq)) / ((double)kList.Q[iqpath].Nq);
        q[2] =
          (kList.Q[iqpath].q[2] * (kList.Q[iqpath].Nq - iq) +
           kList.Q[iqpath + 1].q[2] * (iq)) / ((double)kList.Q[iqpath].Nq);
      }
      CreateHamiltonianFromTightBindingDef(A, q, &H);
      Eigensystem(H, &T, val);
      for (ival = 0; ival < A.NFermion; ival++)
      {
        BANDS->VRe[iqtot][ival] = val[ival];
      }
      EigenFunctionToCharacters(A, H, NULL, T, &TempChar);
      for (ichar = 0; ichar < Ncharacter; ichar++)
      {
        for (ival = 0; ival < A.NFermion; ival++)
        {
          Character[ichar].VRe[iqtot][ival] = TempChar.VRe[ival][ichar];
        }
      }
      FreeMatrix(&T);
      FreeMatrix(&H);
      FreeMatrix(&TempChar);
      iqtot += 1;
    }
  }
  return 0;
}

uint32_t CreateUSuperCell(OperatorType U, OperatorType * USuperCell, CrystalStructureSuperCellType SuperCell)
{
  int    iCell, il, iN, io;
  uint16_t *oadd;
  if ((oadd = (uint16_t*)malloc(U.MaxLength * sizeof(uint16_t))) == NULL)
  {
    printf("malloc for oadd failed in CreateUSuperCell\n");
    fflush(stdout);
    return 1;
  };
  if (InitOperatorCopySettings(USuperCell, U))
  {
    printf("InitOperatorCopySettings failed in CreateUSuperCell\n");
    fflush(stdout);
    return 1;
  };
  USuperCell->NFermion *= SuperCell.NCell;
  USuperCell->NBoson *= SuperCell.NCell;
  USuperCell->NWord = USuperCell->NBoson + Min(USuperCell->NFermion, 1 + (USuperCell->NFermion - 1) / 8);
  for (iCell = 0; iCell < SuperCell.NCell; iCell++)
  {
    for (il = 1; il <= U.MaxLength; il++)
    {
      for (iN = 0; iN < U.l[il].N; iN++)
      {
        for (io = 0; io < il; io++)
        {
          oadd[io] = U.l[il].o[io + il * iN];
          if ((oadd[io] & NotCrAnMask) < U.NFermion)
          {                     /* oadd[io] deals with a fermion */
            oadd[io] += U.NFermion * iCell;
          }
          else
          {                     /* oadd[io] deals with a boson */
            oadd[io] -= U.NFermion;
            oadd[io] += USuperCell->NFermion;
            oadd[io] += U.NBoson * iCell;
          }
        }

        if (U.l[il].QComplex == 0)
        {
          OperatorAddLaderNormalOrder(USuperCell, oadd, il, U.l[il].VRe[iN], 0, QDetNotOrderedFcFaBcBa);
        }
        else
        {
          OperatorAddLaderNormalOrder(USuperCell, oadd, il, U.l[il].V[iN].Re, U.l[il].V[iN].Im, QDetNotOrderedFcFaBcBa);
        }
      }
    }
  }
  return 0;
}



int CreateCrystalStructureFromCrystalStructureSuperCell(CrystalStructureSuperCellType SuperCell, CrystalStructureType * Cell)
{
  int    iCell, iAtom, iAtomsuper, ishell, iorb, ia, ib, ic;
  double transvec[3];
  Cell->Units = SuperCell.Crystal.Units;
  /* A = N[0][0] a + N[0][1] b + N[0][2] c */
  Cell->rUnitCell[0][0] = SuperCell.rSuperCell[0][0] * SuperCell.Crystal.rUnitCell[0][0] + SuperCell.rSuperCell[0][1] * SuperCell.Crystal.rUnitCell[1][0] + SuperCell.rSuperCell[0][2] * SuperCell.Crystal.rUnitCell[2][0];
  Cell->rUnitCell[0][1] = SuperCell.rSuperCell[0][0] * SuperCell.Crystal.rUnitCell[0][1] + SuperCell.rSuperCell[0][1] * SuperCell.Crystal.rUnitCell[1][1] + SuperCell.rSuperCell[0][2] * SuperCell.Crystal.rUnitCell[2][1];
  Cell->rUnitCell[0][2] = SuperCell.rSuperCell[0][0] * SuperCell.Crystal.rUnitCell[0][2] + SuperCell.rSuperCell[0][1] * SuperCell.Crystal.rUnitCell[1][2] + SuperCell.rSuperCell[0][2] * SuperCell.Crystal.rUnitCell[2][2];
  /* A = N[1][0] a + N[1][1] b + N[1][2] c */
  Cell->rUnitCell[1][0] = SuperCell.rSuperCell[1][0] * SuperCell.Crystal.rUnitCell[0][0] + SuperCell.rSuperCell[1][1] * SuperCell.Crystal.rUnitCell[1][0] + SuperCell.rSuperCell[1][2] * SuperCell.Crystal.rUnitCell[2][0];
  Cell->rUnitCell[1][1] = SuperCell.rSuperCell[1][0] * SuperCell.Crystal.rUnitCell[0][1] + SuperCell.rSuperCell[1][1] * SuperCell.Crystal.rUnitCell[1][1] + SuperCell.rSuperCell[1][2] * SuperCell.Crystal.rUnitCell[2][1];
  Cell->rUnitCell[1][2] = SuperCell.rSuperCell[1][0] * SuperCell.Crystal.rUnitCell[0][2] + SuperCell.rSuperCell[1][1] * SuperCell.Crystal.rUnitCell[1][2] + SuperCell.rSuperCell[1][2] * SuperCell.Crystal.rUnitCell[2][2];
  /* A = N[2][0] a + N[2][1] b + N[2][2] c */
  Cell->rUnitCell[2][0] = SuperCell.rSuperCell[2][0] * SuperCell.Crystal.rUnitCell[0][0] + SuperCell.rSuperCell[2][1] * SuperCell.Crystal.rUnitCell[1][0] + SuperCell.rSuperCell[2][2] * SuperCell.Crystal.rUnitCell[2][0];
  Cell->rUnitCell[2][1] = SuperCell.rSuperCell[2][0] * SuperCell.Crystal.rUnitCell[0][1] + SuperCell.rSuperCell[2][1] * SuperCell.Crystal.rUnitCell[1][1] + SuperCell.rSuperCell[2][2] * SuperCell.Crystal.rUnitCell[2][1];
  Cell->rUnitCell[2][2] = SuperCell.rSuperCell[2][0] * SuperCell.Crystal.rUnitCell[0][2] + SuperCell.rSuperCell[2][1] * SuperCell.Crystal.rUnitCell[1][2] + SuperCell.rSuperCell[2][2] * SuperCell.Crystal.rUnitCell[2][2];
  ReciprocalLatticeVector(Cell->rUnitCell, Cell->gUnitCell,Cell->Units);
  Cell->NAtoms = SuperCell.NCell * SuperCell.Crystal.NAtoms;
  Cell->Atom = (AtomType*)malloc(Cell->NAtoms * sizeof(AtomType));
  iAtomsuper = 0;
  for (iCell = 0; iCell < SuperCell.NCell; iCell++)
  {
    ia = SuperCell.CellPos[iCell][0];
    ib = SuperCell.CellPos[iCell][1];
    ic = SuperCell.CellPos[iCell][2];
    transvec[0] = ia * SuperCell.Crystal.rUnitCell[0][0] + ib * SuperCell.Crystal.rUnitCell[1][0] + ic * SuperCell.Crystal.rUnitCell[2][0];
    transvec[1] = ia * SuperCell.Crystal.rUnitCell[0][1] + ib * SuperCell.Crystal.rUnitCell[1][1] + ic * SuperCell.Crystal.rUnitCell[2][1];
    transvec[2] = ia * SuperCell.Crystal.rUnitCell[0][2] + ib * SuperCell.Crystal.rUnitCell[1][2] + ic * SuperCell.Crystal.rUnitCell[2][2];
    for (iAtom = 0; iAtom < SuperCell.Crystal.NAtoms; iAtom++)
    {
      strncpy(Cell->Atom[iAtomsuper].Name, SuperCell.Crystal.Atom[iAtom].Name, 32);
      Cell->Atom[iAtomsuper].AtomNumber = SuperCell.Crystal.Atom[iAtom].AtomNumber;
      Cell->Atom[iAtomsuper].Position[0] = transvec[0] + SuperCell.Crystal.Atom[iAtom].Position[0];
      Cell->Atom[iAtomsuper].Position[1] = transvec[1] + SuperCell.Crystal.Atom[iAtom].Position[1];
      Cell->Atom[iAtomsuper].Position[2] = transvec[2] + SuperCell.Crystal.Atom[iAtom].Position[2];
      Cell->Atom[iAtomsuper].Nshell = SuperCell.Crystal.Atom[iAtom].Nshell;
      Cell->Atom[iAtomsuper].ShellName = (char**)malloc(SuperCell.Crystal.Atom[iAtom].Nshell * sizeof(char*));
      Cell->Atom[iAtomsuper].shelldegeneracy = (int*)malloc(SuperCell.Crystal.Atom[iAtom].Nshell * sizeof(int));
      Cell->Atom[iAtomsuper].OrbitalName = (char***)malloc(SuperCell.Crystal.Atom[iAtom].Nshell * sizeof(char**));
      for (ishell = 0; ishell < SuperCell.Crystal.Atom[iAtom].Nshell; ishell++)
      {
        Cell->Atom[iAtomsuper].ShellName[ishell] = (char*)malloc(32 * sizeof(char));
        strncpy(Cell->Atom[iAtomsuper].ShellName[ishell], SuperCell.Crystal.Atom[iAtom].ShellName[ishell], 32);
        Cell->Atom[iAtomsuper].shelldegeneracy[ishell] = SuperCell.Crystal.Atom[iAtom].shelldegeneracy[ishell];
        Cell->Atom[iAtomsuper].OrbitalName[ishell] = (char**)malloc(Cell->Atom[iAtomsuper].shelldegeneracy[ishell] * sizeof(char*));
        for(iorb = 0; iorb < Cell->Atom[iAtomsuper].shelldegeneracy[ishell]; iorb++)
        {
          Cell->Atom[iAtomsuper].OrbitalName[ishell][iorb] = (char*)malloc(32 * sizeof(char));
          strncpy(Cell->Atom[iAtomsuper].OrbitalName[ishell][iorb], SuperCell.Crystal.Atom[iAtom].OrbitalName[ishell][iorb], 32);
        }
      }
      iAtomsuper++;
    }
  }
  return 0;
}


int CreateTightBindingSuperCell(TightBindingDef A, TightBindingDef * ASuper, CrystalStructureSuperCellType SuperCell)
{
  int    iCell, iAtom, iAtomsuper, ishell, iorb, iH, iHsuper, ia, ib, ic, atomvec[3];
  double transvec[3];
  /* A = N[0][0] a + N[0][1] b + N[0][2] c */
  ASuper->Cell.rUnitCell[0][0] =
    SuperCell.rSuperCell[0][0] * A.Cell.rUnitCell[0][0] + SuperCell.rSuperCell[0][1] * A.Cell.rUnitCell[1][0] +
    SuperCell.rSuperCell[0][2] * A.Cell.rUnitCell[2][0];
  ASuper->Cell.rUnitCell[0][1] =
    SuperCell.rSuperCell[0][0] * A.Cell.rUnitCell[0][1] + SuperCell.rSuperCell[0][1] * A.Cell.rUnitCell[1][1] +
    SuperCell.rSuperCell[0][2] * A.Cell.rUnitCell[2][1];
  ASuper->Cell.rUnitCell[0][2] =
    SuperCell.rSuperCell[0][0] * A.Cell.rUnitCell[0][2] + SuperCell.rSuperCell[0][1] * A.Cell.rUnitCell[1][2] +
    SuperCell.rSuperCell[0][2] * A.Cell.rUnitCell[2][2];
  /* A = N[1][0] a + N[1][1] b + N[1][2] c */
  ASuper->Cell.rUnitCell[1][0] =
    SuperCell.rSuperCell[1][0] * A.Cell.rUnitCell[0][0] + SuperCell.rSuperCell[1][1] * A.Cell.rUnitCell[1][0] +
    SuperCell.rSuperCell[1][2] * A.Cell.rUnitCell[2][0];
  ASuper->Cell.rUnitCell[1][1] =
    SuperCell.rSuperCell[1][0] * A.Cell.rUnitCell[0][1] + SuperCell.rSuperCell[1][1] * A.Cell.rUnitCell[1][1] +
    SuperCell.rSuperCell[1][2] * A.Cell.rUnitCell[2][1];
  ASuper->Cell.rUnitCell[1][2] =
    SuperCell.rSuperCell[1][0] * A.Cell.rUnitCell[0][2] + SuperCell.rSuperCell[1][1] * A.Cell.rUnitCell[1][2] +
    SuperCell.rSuperCell[1][2] * A.Cell.rUnitCell[2][2];
  /* A = N[2][0] a + N[2][1] b + N[2][2] c */
  ASuper->Cell.rUnitCell[2][0] =
    SuperCell.rSuperCell[2][0] * A.Cell.rUnitCell[0][0] + SuperCell.rSuperCell[2][1] * A.Cell.rUnitCell[1][0] +
    SuperCell.rSuperCell[2][2] * A.Cell.rUnitCell[2][0];
  ASuper->Cell.rUnitCell[2][1] =
    SuperCell.rSuperCell[2][0] * A.Cell.rUnitCell[0][1] + SuperCell.rSuperCell[2][1] * A.Cell.rUnitCell[1][1] +
    SuperCell.rSuperCell[2][2] * A.Cell.rUnitCell[2][1];
  ASuper->Cell.rUnitCell[2][2] =
    SuperCell.rSuperCell[2][0] * A.Cell.rUnitCell[0][2] + SuperCell.rSuperCell[2][1] * A.Cell.rUnitCell[1][2] +
    SuperCell.rSuperCell[2][2] * A.Cell.rUnitCell[2][2];
  ReciprocalLatticeVector(ASuper->Cell.rUnitCell, ASuper->Cell.gUnitCell,A.Cell.Units);
  ASuper->Cell.NAtoms = SuperCell.NCell * A.Cell.NAtoms;
  ASuper->Cell.Atom = (AtomType*)malloc(ASuper->Cell.NAtoms * sizeof(AtomType));
  ASuper->nH = SuperCell.NCell * A.nH;
  ASuper->H = (TightBindingDefHopping*)malloc(ASuper->nH * sizeof(TightBindingDefHopping));
  iAtomsuper = 0;
  iHsuper = 0;
  for (iCell = 0; iCell < SuperCell.NCell; iCell++)
  {
    ia = SuperCell.CellPos[iCell][0];
    ib = SuperCell.CellPos[iCell][1];
    ic = SuperCell.CellPos[iCell][2];
    transvec[0] = ia * A.Cell.rUnitCell[0][0] + ib * A.Cell.rUnitCell[1][0] + ic * A.Cell.rUnitCell[2][0];
    transvec[1] = ia * A.Cell.rUnitCell[0][1] + ib * A.Cell.rUnitCell[1][1] + ic * A.Cell.rUnitCell[2][1];
    transvec[2] = ia * A.Cell.rUnitCell[0][2] + ib * A.Cell.rUnitCell[1][2] + ic * A.Cell.rUnitCell[2][2];
    for (iAtom = 0; iAtom < A.Cell.NAtoms; iAtom++)
    {
      strncpy(ASuper->Cell.Atom[iAtomsuper].Name, A.Cell.Atom[iAtom].Name, MAXStrLength);
      ASuper->Cell.Atom[iAtomsuper].AtomNumber = A.Cell.Atom[iAtom].AtomNumber;
      ASuper->Cell.Atom[iAtomsuper].Position[0] = transvec[0] + A.Cell.Atom[iAtom].Position[0];
      ASuper->Cell.Atom[iAtomsuper].Position[1] = transvec[1] + A.Cell.Atom[iAtom].Position[1];
      ASuper->Cell.Atom[iAtomsuper].Position[2] = transvec[2] + A.Cell.Atom[iAtom].Position[2];
      ASuper->Cell.Atom[iAtomsuper].Nshell = A.Cell.Atom[iAtom].Nshell;
      ASuper->Cell.Atom[iAtomsuper].shelldegeneracy = (int*)malloc(ASuper->Cell.Atom[iAtomsuper].Nshell * sizeof(int));
      ASuper->Cell.Atom[iAtomsuper].ShellName= (char**)malloc(ASuper->Cell.Atom[iAtomsuper].Nshell*sizeof(char*));
      ASuper->Cell.Atom[iAtomsuper].OrbitalName= (char***)malloc(ASuper->Cell.Atom[iAtomsuper].Nshell*sizeof(char**));
      for (ishell = 0; ishell < ASuper->Cell.Atom[iAtomsuper].Nshell; ishell++)
      {
        ASuper->Cell.Atom[iAtomsuper].ShellName[ishell]=(char*)malloc(MAXStrLength);
        strncpy(ASuper->Cell.Atom[iAtomsuper].ShellName[ishell], A.Cell.Atom[iAtom].ShellName[ishell], MAXStrLength);
        ASuper->Cell.Atom[iAtomsuper].shelldegeneracy[ishell] = A.Cell.Atom[iAtom].shelldegeneracy[ishell];
        ASuper->Cell.Atom[iAtomsuper].OrbitalName[ishell]= (char**)malloc(A.Cell.Atom[iAtom].shelldegeneracy[ishell]*sizeof(char*));
        for(iorb=0; iorb<A.Cell.Atom[iAtom].shelldegeneracy[ishell]; iorb++)
        {
          ASuper->Cell.Atom[iAtomsuper].OrbitalName[ishell][iorb]=(char*)malloc(MAXStrLength);
          strncpy(ASuper->Cell.Atom[iAtomsuper].OrbitalName[ishell][iorb], A.Cell.Atom[iAtom].OrbitalName[ishell][iorb], MAXStrLength);
        }
      }
      iAtomsuper++;
    }
    for (iH = 0; iH < A.nH; iH++)
    {
      atomvec[0] = ia + A.H[iH].TransVecs[0];
      atomvec[1] = ib + A.H[iH].TransVecs[1];
      atomvec[2] = ic + A.H[iH].TransVecs[2];
      if(SuperCell.Crystal.Units&1)
      {
        if (Dot(atomvec, SuperCell.gSuperCell[0]) >= -0.000000001)
        {
          ASuper->H[iHsuper].TransVecs[0] = (int)(Dot(atomvec, SuperCell.gSuperCell[0]) / (1.999999999 * Pi));
        }
        else
        {
          ASuper->H[iHsuper].TransVecs[0] = (int)(Dot(atomvec, SuperCell.gSuperCell[0]) / (2.000000001 * Pi)) - 1;
        }
        if (Dot(atomvec, SuperCell.gSuperCell[1]) >= -0.000000001)
        {
          ASuper->H[iHsuper].TransVecs[1] = (int)(Dot(atomvec, SuperCell.gSuperCell[1]) / (1.999999999 * Pi));
        }
        else
        {
          ASuper->H[iHsuper].TransVecs[1] = (int)(Dot(atomvec, SuperCell.gSuperCell[1]) / (2.000000001 * Pi)) - 1;
        }
        if (Dot(atomvec, SuperCell.gSuperCell[2]) >= -0.000000001)
        {
          ASuper->H[iHsuper].TransVecs[2] = (int)(Dot(atomvec, SuperCell.gSuperCell[2]) / (1.999999999 * Pi));
        }
        else
        {
          ASuper->H[iHsuper].TransVecs[2] = (int)(Dot(atomvec, SuperCell.gSuperCell[2]) / (2.000000001 * Pi)) - 1;
        }
      }
      else
      {
        if (Dot(atomvec, SuperCell.gSuperCell[0]) >= -0.000000001)
        {
          ASuper->H[iHsuper].TransVecs[0] = (int)(Dot(atomvec, SuperCell.gSuperCell[0]) / (0.999999999));
        }
        else
        {
          ASuper->H[iHsuper].TransVecs[0] = (int)(Dot(atomvec, SuperCell.gSuperCell[0]) / (1.000000001)) - 1;
        }
        if (Dot(atomvec, SuperCell.gSuperCell[1]) >= -0.000000001)
        {
          ASuper->H[iHsuper].TransVecs[1] = (int)(Dot(atomvec, SuperCell.gSuperCell[1]) / (0.999999999));
        }
        else
        {
          ASuper->H[iHsuper].TransVecs[1] = (int)(Dot(atomvec, SuperCell.gSuperCell[1]) / (1.000000001)) - 1;
        }
        if (Dot(atomvec, SuperCell.gSuperCell[2]) >= -0.000000001)
        {
          ASuper->H[iHsuper].TransVecs[2] = (int)(Dot(atomvec, SuperCell.gSuperCell[2]) / (0.999999999));
        }
        else
        {
          ASuper->H[iHsuper].TransVecs[2] = (int)(Dot(atomvec, SuperCell.gSuperCell[2]) / (1.000000001)) - 1;
        }
      }
      atomvec[0] -=
        ASuper->H[iHsuper].TransVecs[0] * SuperCell.rSuperCell[0][0] +
        ASuper->H[iHsuper].TransVecs[1] * SuperCell.rSuperCell[1][0] +
        ASuper->H[iHsuper].TransVecs[2] * SuperCell.rSuperCell[2][0];
      atomvec[1] -=
        ASuper->H[iHsuper].TransVecs[0] * SuperCell.rSuperCell[0][1] +
        ASuper->H[iHsuper].TransVecs[1] * SuperCell.rSuperCell[1][1] +
        ASuper->H[iHsuper].TransVecs[2] * SuperCell.rSuperCell[2][1];
      atomvec[2] -=
        ASuper->H[iHsuper].TransVecs[0] * SuperCell.rSuperCell[0][2] +
        ASuper->H[iHsuper].TransVecs[1] * SuperCell.rSuperCell[1][2] +
        ASuper->H[iHsuper].TransVecs[2] * SuperCell.rSuperCell[2][2];
      FindCellPosInCrystalStructureSuperCellType(SuperCell, atomvec, &(ASuper->H[iHsuper].Atompos2));
      ASuper->H[iHsuper].Atompos2 *= A.Cell.NAtoms;
      ASuper->H[iHsuper].Atompos2 += A.H[iH].Atompos2;
      ASuper->H[iHsuper].Shellpos1 = A.H[iH].Shellpos1;
      ASuper->H[iHsuper].Shellpos2 = A.H[iH].Shellpos2;
      ASuper->H[iHsuper].Atompos1 = A.H[iH].Atompos1 + iCell * A.Cell.NAtoms;
      MatrixCopy(&(ASuper->H[iHsuper].Hop), A.H[iH].Hop);
      iHsuper++;
    }
  }
  SetNFermionIndexTightBindingDef(ASuper);
  strncpy(ASuper->Name, A.Name, MAXStrLength-1);
  return 0;
}

int CreateTightBindingSquareSuperCell(TightBindingDef A, TightBindingDef * ASuper, int N[3])
{
  int    ia, ib, ic, iAtom, iAtomsuper, ishell, iH, iHsuper, ia2, ib2, ic2;
  double transvec[3];
  ASuper->Cell.rUnitCell[0][0] = N[0] * A.Cell.rUnitCell[0][0];
  ASuper->Cell.rUnitCell[0][1] = N[0] * A.Cell.rUnitCell[0][1];
  ASuper->Cell.rUnitCell[0][2] = N[0] * A.Cell.rUnitCell[0][2];
  ASuper->Cell.rUnitCell[1][0] = N[1] * A.Cell.rUnitCell[1][0];
  ASuper->Cell.rUnitCell[1][1] = N[1] * A.Cell.rUnitCell[1][1];
  ASuper->Cell.rUnitCell[1][2] = N[1] * A.Cell.rUnitCell[1][2];
  ASuper->Cell.rUnitCell[2][0] = N[2] * A.Cell.rUnitCell[2][0];
  ASuper->Cell.rUnitCell[2][1] = N[2] * A.Cell.rUnitCell[2][1];
  ASuper->Cell.rUnitCell[2][2] = N[2] * A.Cell.rUnitCell[2][2];
  ReciprocalLatticeVector(ASuper->Cell.rUnitCell, ASuper->Cell.gUnitCell,A.Cell.Units);
  ASuper->Cell.NAtoms = N[0] * N[1] * N[2] * A.Cell.NAtoms;
  ASuper->Cell.Atom = (AtomType*)malloc(ASuper->Cell.NAtoms * sizeof(AtomType));
  ASuper->nH = N[0] * N[1] * N[2] * A.nH;
  ASuper->H = (TightBindingDefHopping*)malloc(ASuper->nH * sizeof(TightBindingDefHopping));
  iAtomsuper = 0;
  iHsuper = 0;
  for (ia = 0; ia < N[0]; ia++)
  {
    for (ib = 0; ib < N[1]; ib++)
    {
      for (ic = 0; ic < N[2]; ic++)
      {
        transvec[0] = ia * A.Cell.rUnitCell[0][0] + ib * A.Cell.rUnitCell[1][0] + ic * A.Cell.rUnitCell[2][0];
        transvec[1] = ia * A.Cell.rUnitCell[0][1] + ib * A.Cell.rUnitCell[1][1] + ic * A.Cell.rUnitCell[2][1];
        transvec[2] = ia * A.Cell.rUnitCell[0][2] + ib * A.Cell.rUnitCell[1][2] + ic * A.Cell.rUnitCell[2][2];
        for (iAtom = 0; iAtom < A.Cell.NAtoms; iAtom++)
        {
          strncpy(ASuper->Cell.Atom[iAtomsuper].Name, A.Cell.Atom[iAtom].Name, 32);
          ASuper->Cell.Atom[iAtomsuper].AtomNumber = A.Cell.Atom[iAtom].AtomNumber;
          ASuper->Cell.Atom[iAtomsuper].Position[0] = transvec[0] + A.Cell.Atom[iAtom].Position[0];
          ASuper->Cell.Atom[iAtomsuper].Position[1] = transvec[1] + A.Cell.Atom[iAtom].Position[1];
          ASuper->Cell.Atom[iAtomsuper].Position[2] = transvec[2] + A.Cell.Atom[iAtom].Position[2];
          ASuper->Cell.Atom[iAtomsuper].Nshell = A.Cell.Atom[iAtom].Nshell;
          ASuper->Cell.Atom[iAtomsuper].shelldegeneracy = (int*)malloc(ASuper->Cell.Atom[iAtomsuper].Nshell * sizeof(int));
          for (ishell = 0; ishell < ASuper->Cell.Atom[iAtomsuper].Nshell; ishell++)
          {
            ASuper->Cell.Atom[iAtomsuper].shelldegeneracy[ishell] = A.Cell.Atom[iAtom].shelldegeneracy[ishell];
          }
          iAtomsuper++;
        }
        for (iH = 0; iH < A.nH; iH++)
        {
          if ((A.H[iH].TransVecs[0] + ia) >= 0)
          {
            ASuper->H[iHsuper].TransVecs[0] = (A.H[iH].TransVecs[0] + ia) / N[0];
            ia2 = (A.H[iH].TransVecs[0] + ia) % N[0];
          }
          else
          {
            ASuper->H[iHsuper].TransVecs[0] = ((A.H[iH].TransVecs[0] + 1 + ia) / N[0]) - 1;
            ia2 = A.H[iH].TransVecs[0] + ia - ASuper->H[iHsuper].TransVecs[0] * N[0];
          }
          if ((A.H[iH].TransVecs[1] + ib) >= 0)
          {
            ASuper->H[iHsuper].TransVecs[1] = (A.H[iH].TransVecs[1] + ib) / N[1];
            ib2 = (A.H[iH].TransVecs[1] + ib) % N[1];
          }
          else
          {
            ASuper->H[iHsuper].TransVecs[1] = ((A.H[iH].TransVecs[1] + 1 + ib) / N[1]) - 1;
            ib2 = A.H[iH].TransVecs[1] + ib - ASuper->H[iHsuper].TransVecs[1] * N[1];
          }
          if ((A.H[iH].TransVecs[2] + ic) >= 0)
          {
            ASuper->H[iHsuper].TransVecs[2] = (A.H[iH].TransVecs[2] + ic) / N[2];
            ic2 = (A.H[iH].TransVecs[2] + ic) % N[2];
          }
          else
          {
            ASuper->H[iHsuper].TransVecs[2] = ((A.H[iH].TransVecs[2] + 1 + ic) / N[2]) - 1;
            ic2 = A.H[iH].TransVecs[2] + ic - ASuper->H[iHsuper].TransVecs[2] * N[2];
          }
          ASuper->H[iHsuper].Shellpos1 = A.H[iH].Shellpos1;
          ASuper->H[iHsuper].Shellpos2 = A.H[iH].Shellpos2;
          ASuper->H[iHsuper].Atompos1 = A.H[iH].Atompos1 + ic + ib * N[2] + ia * N[2] * N[1];
          ASuper->H[iHsuper].Atompos2 = A.H[iH].Atompos2 + ic2 + ib2 * N[2] + ia2 * N[2] * N[1];
          MatrixCopy(&(ASuper->H[iHsuper].Hop), A.H[iH].Hop);
          iHsuper++;
        }
      }
    }
  }
  SetNFermionIndexTightBindingDef(ASuper);
  return 0;
}


int TightBindingBANDSToMatrix(double Emin, double Emax, int NEpoints, MatrixType BANDS, MatrixType * BWValue)
{
  int    iq, jBand, jE;
  double Enrescale, frac;
  BWValue->Ni = BANDS.Ni;
  BWValue->Nj = NEpoints + 1;
  BWValue->QComplex = 0;
  if (InitMatrix(BWValue))
    return 1;
  for (iq = 0; iq < BANDS.Ni; iq++)
  {
    for (jBand = 0; jBand < BANDS.Nj; jBand++)
    {
      Enrescale = (BANDS.VRe[iq][jBand] - Emin) / (Emax - Emin);
      Enrescale *= (NEpoints + 1);
      jE = (int)Enrescale;
      frac = Enrescale - jE;
      if ((jE >= 0) && (jE <= NEpoints))
        BWValue->VRe[iq][jE] += 1;
      if (((jE - 1) >= 0) && ((jE - 1) <= NEpoints))
        BWValue->VRe[iq][jE - 1] += (1 - frac) * 1;
      if (((jE + 1) >= 0) && ((jE + 1) <= NEpoints))
        BWValue->VRe[iq][jE + 1] += (0 + frac) * 1;
    }
  }
  for (iq = 1; iq < BANDS.Ni; iq++)
  {
    for (jBand = 0; jBand < BANDS.Nj; jBand++)
    {
      Enrescale = (BANDS.VRe[iq][jBand] - Emin) / (Emax - Emin);
      Enrescale *= (NEpoints + 1);
      jE = (int)Enrescale;
      frac = Enrescale - jE;
      if ((jE >= 0) && (jE <= NEpoints))
        BWValue->VRe[iq - 1][jE] += 0.5 * 1;
      if (((jE - 1) >= 0) && ((jE - 1) <= NEpoints))
        BWValue->VRe[iq - 1][jE - 1] += 0.5 * (1 - frac) * 1;
      if (((jE + 1) >= 0) && ((jE + 1) <= NEpoints))
        BWValue->VRe[iq - 1][jE + 1] += 0.5 * (0 + frac) * 1;
    }
  }
  for (iq = 0; iq < BANDS.Ni - 1; iq++)
  {
    for (jBand = 0; jBand < BANDS.Nj; jBand++)
    {
      Enrescale = (BANDS.VRe[iq][jBand] - Emin) / (Emax - Emin);
      Enrescale *= (NEpoints + 1);
      jE = (int)Enrescale;
      frac = Enrescale - jE;
      if ((jE >= 0) && (jE <= NEpoints))
        BWValue->VRe[iq + 1][jE] += 0.5 * 1;
      if (((jE - 1) >= 0) && ((jE - 1) <= NEpoints))
        BWValue->VRe[iq + 1][jE - 1] += 0.5 * (1 - frac) * 1;
      if (((jE + 1) >= 0) && ((jE + 1) <= NEpoints))
        BWValue->VRe[iq + 1][jE + 1] += 0.5 * (0 + frac) * 1;
    }
  }
  return 0;
}

int TightBindingBANDSCharacterToMatrix(double Emin, double Emax, int NEpoints, MatrixType BANDS, MatrixType Character,
                                       MatrixType * BWValue)
{
  int    iq, jBand, jE;
  double Enrescale, frac;
  BWValue->Ni = BANDS.Ni;
  BWValue->Nj = NEpoints + 1;
  BWValue->QComplex = 0;
  if (InitMatrix(BWValue))
    return 1;
  for (iq = 0; iq < BANDS.Ni; iq++)
  {
    for (jBand = 0; jBand < BANDS.Nj; jBand++)
    {
      Enrescale = (BANDS.VRe[iq][jBand] - Emin) / (Emax - Emin);
      Enrescale *= (NEpoints + 1);
      jE = (int)floor(Enrescale);
      frac = Enrescale - jE;
      if ((jE >= 0) && (jE <= NEpoints))
        BWValue->VRe[iq][jE] += Character.VRe[iq][jBand];
      if (((jE - 1) >= 0) && ((jE - 1) <= NEpoints))
        BWValue->VRe[iq][jE - 1] += (1 - frac) * Character.VRe[iq][jBand];
      if (((jE + 1) >= 0) && ((jE + 1) <= NEpoints))
        BWValue->VRe[iq][jE + 1] += (0 + frac) * Character.VRe[iq][jBand];
    }
  }
  for (iq = 1; iq < BANDS.Ni; iq++)
  {
    for (jBand = 0; jBand < BANDS.Nj; jBand++)
    {
      Enrescale = (BANDS.VRe[iq][jBand] - Emin) / (Emax - Emin);
      Enrescale *= (NEpoints + 1);
      jE = (int)floor(Enrescale);
      frac = Enrescale - jE;
      if ((jE >= 0) && (jE <= NEpoints))
        BWValue->VRe[iq - 1][jE] += 0.5 * Character.VRe[iq][jBand];
      if (((jE - 1) >= 0) && ((jE - 1) <= NEpoints))
        BWValue->VRe[iq - 1][jE - 1] += 0.5 * (1 - frac) * Character.VRe[iq][jBand];
      if (((jE + 1) >= 0) && ((jE + 1) <= NEpoints))
        BWValue->VRe[iq - 1][jE + 1] += 0.5 * (0 + frac) * Character.VRe[iq][jBand];
    }
  }
  for (iq = 0; iq < BANDS.Ni - 1; iq++)
  {
    for (jBand = 0; jBand < BANDS.Nj; jBand++)
    {
      Enrescale = (BANDS.VRe[iq][jBand] - Emin) / (Emax - Emin);
      Enrescale *= (NEpoints + 1);
      jE = (int)floor(Enrescale);
      frac = Enrescale - jE;
      if ((jE >= 0) && (jE <= NEpoints))
        BWValue->VRe[iq + 1][jE] += 0.5 * Character.VRe[iq][jBand];
      if (((jE - 1) >= 0) && ((jE - 1) <= NEpoints))
        BWValue->VRe[iq + 1][jE - 1] += 0.5 * (1 - frac) * Character.VRe[iq][jBand];
      if (((jE + 1) >= 0) && ((jE + 1) <= NEpoints))
        BWValue->VRe[iq + 1][jE + 1] += 0.5 * (0 + frac) * Character.VRe[iq][jBand];
    }
  }
  return 0;
}

int TightBindingAddAxesToMatrix(double Emin, double Emax, int NELines, double *Elines, BandstructurePathType kList,
                                MatrixType * RGBValue)
{
  int    iQ, i, j, jE;
  i = 0;
  for (iQ = 0; iQ < kList.NQ; iQ++)
  {
    for (j = 0; j < RGBValue->Nj; j++)
    {
      RGBValue->VRe[i][j] += 0.5;
    }
    i += kList.Q[iQ].Nq;
  }
  for (jE = 0; jE < NELines; jE++)
  {
    j = (int)(((double)(RGBValue->Nj) - 1) * (Elines[jE] - Emin) / (Emax - Emin));
    if ((j < 0) || (j >= RGBValue->Nj))
    {
      printf("Value outside range j=%i jE=%i Elines[jE]=%22.15E\n", j, jE, Elines[jE]);
    }
    else
    {
      for (i = 0; i < RGBValue->Ni; i++)
      {
        RGBValue->VRe[i][j] += 0.5;
      }
    }
  }
  return 0;
}
