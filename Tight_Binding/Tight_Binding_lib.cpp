uint32_t AtomNameToNumber(const char* Name)
{
  uint32_t i;
  const char *PeriodicTable[] = {"H" ,                                                                                                                                                          "He",
                                 "Li","Be",                                                                                                                            "B", "C", "N", "O", "F", "Ne",
                                 "Na","Mg",                                                                                                                           "Al","Si", "P", "S","Cl", "Ar",
                                 "K", "Ca",                                                                        "Sc","Ti", "V","Cr","Mn","Fe","Co","Ni","Cu","Zn", "Ga","Ge","As","Se","Br", "Kr",
                                 "Rb","Sr",                                                                         "Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd", "In","Sn","Sb","Te", "I", "Xe",
                                 "Cs","Ba", "La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb", "Lu","Hf","Ta", "W","Re","Os","Ir","Pt","Au","Hg", "Tl","Pb","Bi","Po","At", "Rn",
                                 "Fr","Ra", "Ac","Th","Pa", "U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No", "Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn"};
  for(i=0;i<112;i++)
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

int CopyAtomType(const AtomType *const A, AtomType *B)
{
  int ishell,iorb;
  strcpy(B->Name, A->Name);
  B->AtomNumber=A->AtomNumber;
  B->Position[0]=A->Position[0];
  B->Position[1]=A->Position[1];
  B->Position[2]=A->Position[2];
  B->Nshell=A->Nshell;
  B->ShellName=(char**)malloc(A->Nshell*sizeof(char*));
  B->shelldegeneracy=(int*)malloc(A->Nshell*sizeof(int));
  B->OrbitalName=(char***)malloc(A->Nshell*sizeof(char**));
  for(ishell=0;ishell<B->Nshell;ishell++)
  {
    B->ShellName[ishell]=(char*)malloc((strlen(A->ShellName[ishell])+1)*sizeof(char));
    strcpy(B->ShellName[ishell], A->ShellName[ishell]);
    B->shelldegeneracy[ishell]=A->shelldegeneracy[ishell];
    B->OrbitalName[ishell]=(char**)malloc(B->shelldegeneracy[ishell]*sizeof(char*));
    for(iorb=0;iorb<B->shelldegeneracy[ishell];iorb++)
    {
      B->OrbitalName[ishell][iorb]=(char*)malloc((strlen(A->OrbitalName[ishell][iorb])+1)*sizeof(char));
      strcpy(B->OrbitalName[ishell][iorb], A->OrbitalName[ishell][iorb]);
    }
  }
  return 0;
}

bool QInteger(double x, double epsilon)
{
  if(fabs(x-round(x))>epsilon) return false;
  return true;
}

int FindAtomPosModulusUnitCell(const AtomType *const Atom, const CrystalStructureType *const Cell, double epsilon)
{
  int i;
  double vec[3],relvec[3];
  for(i=0;i<Cell->NAtoms;i++)
  {
    if(strcmp(Atom->Name, Cell->Atom[i].Name)==0)
    { // name is equal next test if position is compatible
      vec[0]=Atom->Position[0]-Cell->Atom[i].Position[0];
      vec[1]=Atom->Position[1]-Cell->Atom[i].Position[1];
      vec[2]=Atom->Position[2]-Cell->Atom[i].Position[2];
      if(Cell->Units&8)
      { // position in units of lattice vectors vec must be integer in all elements
        if(QInteger(vec[0],epsilon)&&QInteger(vec[1],epsilon)&&QInteger(vec[2],epsilon)) return i;
      }
      else
      { // position in absolute units
        if((fabs(vec[0])<epsilon)&&(fabs(vec[1])<epsilon)&&(fabs(vec[2])<epsilon)) return i;
        if(fabs(ParallelalPipedVolume(Cell->gUnitCell[0], Cell->gUnitCell[1], Cell->gUnitCell[2]))>epsilon)
        {
          relvec[0]=Dot(vec,Cell->gUnitCell[0]);
          relvec[1]=Dot(vec,Cell->gUnitCell[1]);
          relvec[2]=Dot(vec,Cell->gUnitCell[2]);
          if(Cell->Units&1)
          {
            relvec[0]/=(2*Pi);
            relvec[1]/=(2*Pi);
            relvec[2]/=(2*Pi);
          }
          if(QInteger(relvec[0],epsilon)&&QInteger(relvec[1],epsilon)&&QInteger(relvec[2],epsilon)) return i;
        }
      }
    }
  }
  // atom not found
  return INT_MAX;
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
    vector[0] += ((A.H[iH].TransVecs[0]==0) ? (0) : (A.Cell.rUnitCell[0][0] * A.H[iH].TransVecs[0])) + ((A.H[iH].TransVecs[1]==0) ? (0) : (A.Cell.rUnitCell[1][0] * A.H[iH].TransVecs[1])) + ((A.H[iH].TransVecs[2]==0) ? (0) : (A.Cell.rUnitCell[2][0] * A.H[iH].TransVecs[2]));
    vector[1] += ((A.H[iH].TransVecs[0]==0) ? (0) : (A.Cell.rUnitCell[0][1] * A.H[iH].TransVecs[0])) + ((A.H[iH].TransVecs[1]==0) ? (0) : (A.Cell.rUnitCell[1][1] * A.H[iH].TransVecs[1])) + ((A.H[iH].TransVecs[2]==0) ? (0) : (A.Cell.rUnitCell[2][1] * A.H[iH].TransVecs[2]));
    vector[2] += ((A.H[iH].TransVecs[0]==0) ? (0) : (A.Cell.rUnitCell[0][2] * A.H[iH].TransVecs[0])) + ((A.H[iH].TransVecs[1]==0) ? (0) : (A.Cell.rUnitCell[1][2] * A.H[iH].TransVecs[1])) + ((A.H[iH].TransVecs[2]==0) ? (0) : (A.Cell.rUnitCell[2][2] * A.H[iH].TransVecs[2]));
    alphaRe = cos(Dot(k, vector));
    alphaIm = sin(Dot(k, vector));
    if (Sqr(alphaIm) < Sqr(alphaRe) * DBL_EPSILON * 100) alphaIm = 0;
    MatrixAddShifted(H, A.H[iH].Hop, alphaRe, alphaIm, i, j);
  }
  return 0;
}



/*  given a Tight binding definition and a k vector creates a matrix
 *  representing the Hamiltonian.
 *    two atoms at different positions in the unit cell hop with a constant phase 1
 */
int CreateHamiltonianFromTightBindingDefConstantPhaseInUnitCell(TightBindingDef A, const double k[3], MatrixType * H)
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
    vector[0] = ((A.H[iH].TransVecs[0]==0) ? (0) : (A.Cell.rUnitCell[0][0] * A.H[iH].TransVecs[0])) + ((A.H[iH].TransVecs[1]==0) ? (0) : (A.Cell.rUnitCell[1][0] * A.H[iH].TransVecs[1])) + ((A.H[iH].TransVecs[2]==0) ? (0) : (A.Cell.rUnitCell[2][0] * A.H[iH].TransVecs[2]));
    vector[1] = ((A.H[iH].TransVecs[0]==0) ? (0) : (A.Cell.rUnitCell[0][1] * A.H[iH].TransVecs[0])) + ((A.H[iH].TransVecs[1]==0) ? (0) : (A.Cell.rUnitCell[1][1] * A.H[iH].TransVecs[1])) + ((A.H[iH].TransVecs[2]==0) ? (0) : (A.Cell.rUnitCell[2][1] * A.H[iH].TransVecs[2]));
    vector[2] = ((A.H[iH].TransVecs[0]==0) ? (0) : (A.Cell.rUnitCell[0][2] * A.H[iH].TransVecs[0])) + ((A.H[iH].TransVecs[1]==0) ? (0) : (A.Cell.rUnitCell[1][2] * A.H[iH].TransVecs[1])) + ((A.H[iH].TransVecs[2]==0) ? (0) : (A.Cell.rUnitCell[2][2] * A.H[iH].TransVecs[2]));
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
    vector[0] += ((A.H[iH].TransVecs[0]==0) ? (0) : (A.Cell.rUnitCell[0][0] * A.H[iH].TransVecs[0])) + ((A.H[iH].TransVecs[1]==0) ? (0) : (A.Cell.rUnitCell[1][0] * A.H[iH].TransVecs[1])) + ((A.H[iH].TransVecs[2]==0) ? (0) : (A.Cell.rUnitCell[2][0] * A.H[iH].TransVecs[2]));
    vector[1] += ((A.H[iH].TransVecs[0]==0) ? (0) : (A.Cell.rUnitCell[0][1] * A.H[iH].TransVecs[0])) + ((A.H[iH].TransVecs[1]==0) ? (0) : (A.Cell.rUnitCell[1][1] * A.H[iH].TransVecs[1])) + ((A.H[iH].TransVecs[2]==0) ? (0) : (A.Cell.rUnitCell[2][1] * A.H[iH].TransVecs[2]));
    vector[2] += ((A.H[iH].TransVecs[0]==0) ? (0) : (A.Cell.rUnitCell[0][2] * A.H[iH].TransVecs[0])) + ((A.H[iH].TransVecs[1]==0) ? (0) : (A.Cell.rUnitCell[1][2] * A.H[iH].TransVecs[1])) + ((A.H[iH].TransVecs[2]==0) ? (0) : (A.Cell.rUnitCell[2][2] * A.H[iH].TransVecs[2]));
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
    vector[0] += ((A.H[iH].TransVecs[0]==0) ? (0) : (A.Cell.rUnitCell[0][0] * A.H[iH].TransVecs[0])) + ((A.H[iH].TransVecs[1]==0) ? (0) : (A.Cell.rUnitCell[1][0] * A.H[iH].TransVecs[1])) + ((A.H[iH].TransVecs[2]==0) ? (0) : (A.Cell.rUnitCell[2][0] * A.H[iH].TransVecs[2]));
    vector[1] += ((A.H[iH].TransVecs[0]==0) ? (0) : (A.Cell.rUnitCell[0][1] * A.H[iH].TransVecs[0])) + ((A.H[iH].TransVecs[1]==0) ? (0) : (A.Cell.rUnitCell[1][1] * A.H[iH].TransVecs[1])) + ((A.H[iH].TransVecs[2]==0) ? (0) : (A.Cell.rUnitCell[2][1] * A.H[iH].TransVecs[2]));
    vector[2] += ((A.H[iH].TransVecs[0]==0) ? (0) : (A.Cell.rUnitCell[0][2] * A.H[iH].TransVecs[0])) + ((A.H[iH].TransVecs[1]==0) ? (0) : (A.Cell.rUnitCell[1][2] * A.H[iH].TransVecs[1])) + ((A.H[iH].TransVecs[2]==0) ? (0) : (A.Cell.rUnitCell[2][2] * A.H[iH].TransVecs[2]));
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
  bool Converged=false;
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
      Converged=true;
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


/*  given an eigen function the character matrix is calculated as char_{i,j} = psi_i^* * psi_j
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

/*  given an eigen function the character matrix is calculated as char_{i,j} = psi_i^* * psi_j
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
      Char->VIm[i][j] = EigenFun[i].Re*EigenFun[j].Im-EigenFun[i].Im*EigenFun[j].Re;
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
      Char->V[k].Im = EigenFun[i].Re*EigenFun[j].Im-EigenFun[i].Im*EigenFun[j].Re;
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

/*  given an eigen function the character matrix is calculated as char_{i,j} = psi_i^* * psi_j
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
      Char->VIm[i][j] = EigenFun[i].Re*EigenFun[j].Im-EigenFun[i].Im*EigenFun[j].Re;
    }
  }
}

/*  given an eigen function the character matrix is calculated as char_{i,j} = psi_i^* * psi_j
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

int EigenFunctionToCharactersStandard(TightBindingDef A, MatrixType H, double *extra, MatrixType EigenFun, MatrixType * Character)
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

int EigenFunctionToCharactersOnsiteEnergy(TightBindingDef A, MatrixType H, double *extra, MatrixType EigenFun, MatrixType * Character)
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

int EigenFunctionToCharactersOnsiteEnergyPhase(TightBindingDef A, MatrixType H, double *extra, MatrixType EigenFun, MatrixType * Character)
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

int EigenFunctionToCharactersOnsiteEnergyPhaseRGB(TightBindingDef A, MatrixType H, double *extra, MatrixType EigenFun, MatrixType * Character)
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

int EigenFunctionToCharactersOnsiteEnergyPhaseRGBkshifted(TightBindingDef A, MatrixType H, double kshift[3], MatrixType EigenFun, MatrixType * Character)
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
int TightBindingFermiCharacter(TightBindingDef A, double Emin, double Emax, int NE, double k0[3], double k1[3], double k2[3], MatrixType Fermi[], int Ncharacter, int (*EigenFunctionToCharacters) (TightBindingDef, MatrixType, double *, MatrixType, MatrixType *))
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
int TightBindingFermiCharacterUnfolding(TightBindingDef A, double Emin, double Emax, int NE, double k0[3], double k1[3], double k2[3], int Nk1, int Nk2, MatrixType Fermi[], int Ncharacter, int (*EigenFunctionToCharacters) (TightBindingDef, MatrixType, double *, MatrixType, MatrixType *))
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
    // set information about atoms in supercell
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
    // set information about hopping
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
      // set Connecting Vector (redundant information, but usefull)
      ASuper->H[iHsuper].ConnectingVector[0] = (ASuper->Cell.Atom[ASuper->H[iHsuper].Atompos2].Position[0]-ASuper->Cell.Atom[ASuper->H[iHsuper].Atompos1].Position[0]) + ASuper->H[iHsuper].TransVecs[0] * ASuper->Cell.rUnitCell[0][0] + ASuper->H[iHsuper].TransVecs[1] * ASuper->Cell.rUnitCell[1][0] + ASuper->H[iHsuper].TransVecs[2] * ASuper->Cell.rUnitCell[2][0] ;
      ASuper->H[iHsuper].ConnectingVector[1] = (ASuper->Cell.Atom[ASuper->H[iHsuper].Atompos2].Position[1]-ASuper->Cell.Atom[ASuper->H[iHsuper].Atompos1].Position[1]) + ASuper->H[iHsuper].TransVecs[0] * ASuper->Cell.rUnitCell[0][1] + ASuper->H[iHsuper].TransVecs[1] * ASuper->Cell.rUnitCell[1][1] + ASuper->H[iHsuper].TransVecs[2] * ASuper->Cell.rUnitCell[2][1] ;
      ASuper->H[iHsuper].ConnectingVector[2] = (ASuper->Cell.Atom[ASuper->H[iHsuper].Atompos2].Position[2]-ASuper->Cell.Atom[ASuper->H[iHsuper].Atompos1].Position[2]) + ASuper->H[iHsuper].TransVecs[0] * ASuper->Cell.rUnitCell[0][2] + ASuper->H[iHsuper].TransVecs[1] * ASuper->Cell.rUnitCell[1][2] + ASuper->H[iHsuper].TransVecs[2] * ASuper->Cell.rUnitCell[2][2] ;
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

int TightBindingBANDSCharacterToMatrix(double Emin, double Emax, int NEpoints, MatrixType BANDS, MatrixType Character,  MatrixType * BWValue)
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

int TightBindingAddAxesToMatrix(double Emin, double Emax, int NELines, double *Elines, BandstructurePathType kList, MatrixType * RGBValue)
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


/*  Given a crystal structure A, a center c and a radius r
 *  Finds all atoms i in A for which
 *  |p_i - c| < r
 *  here i can belong to different cells, i.e.
 *  p_{i,j} = p0_j + alpha_i . A->Cell
 *  where j enumerates the atoms within a unit cell and i enumerates the different cells
 *
 *
 */
int FindAllAtomsInsideSphere(const CrystalStructureType * const A, double *c, double r, CrystalStructureType *B)
{
  int Na,Nb,Nc,Na0,Nb0,Nc0,ia,ib,ic,iatom;
  double V,O[3],rSqr=Sqr(r),CellPos[3],AtomPos[3];
  // determine the number of unit cells one needs to test and the center around one needs to test
  if(Norm(A->rUnitCell[0])==INFINITY)
  {
    Na=0;
    Na0=0;
    if(Norm(A->rUnitCell[1])==INFINITY)
    {
      Nb=0;
      Nb0=0;
      if(Norm(A->rUnitCell[2])==INFINITY)
      { // a, b and c inf
        Nc=0;
        Nc0=0;
      }
      else
      { // a and b inf c finite
        Nc=(int)floor(fabs(r/Norm(A->rUnitCell[2])))+1;
        Nc0=(int)round(Dot(c,A->rUnitCell[2])/Norm(A->rUnitCell[2]));
      }
    }
    else
    { // a inf b finite
      if(Norm(A->rUnitCell[2])==INFINITY)
      { //a and c inf b finite
        Nc=0;
        Nc0=0;
        Nb=(int)floor(fabs(r/Norm(A->rUnitCell[1])))+1;
        Nb0=(int)round(Dot(c,A->rUnitCell[1])/Norm(A->rUnitCell[1]));
      }
      else
      { // a inf b and c finite
        V=Dot(A->rUnitCell[1],A->rUnitCell[2]);
        Nb=(int)(floor(fabs(r*Norm(A->rUnitCell[2])/V))+1);
        Nc=(int)(floor(fabs(r*Norm(A->rUnitCell[1])/V))+1);
        V=A->rUnitCell[1][1]*A->rUnitCell[2][2]-A->rUnitCell[1][2]*A->rUnitCell[2][1];
        Nb0=(int)round((c[1]*A->rUnitCell[2][2]-c[2]*A->rUnitCell[2][1])/V);
        Nc0=(int)round((c[2]*A->rUnitCell[1][1]-c[1]*A->rUnitCell[1][2])/V);
      }
    }
  }
  else
  { // a finite
    if(Norm(A->rUnitCell[1])==INFINITY)
    { // b inf a finite
      Nb=0;
      Nb0=0;
      if(Norm(A->rUnitCell[2])==INFINITY)
      { // b and c inf a finite
        Nc=0;
        Nc0=0;
        Na=(int)floor(fabs(r/Norm(A->rUnitCell[0])))+1;
        Na0=(int)round(Dot(c,A->rUnitCell[0])/Norm(A->rUnitCell[0]));
      }
      else
      { // b inf a and c finite
        V=Dot(A->rUnitCell[0],A->rUnitCell[2]);
        Na=(int)(floor(fabs(r*Norm(A->rUnitCell[2])/V))+1);
        Nc=(int)(floor(fabs(r*Norm(A->rUnitCell[0])/V))+1);
        V=A->rUnitCell[2][2]*A->rUnitCell[0][0]-A->rUnitCell[2][0]*A->rUnitCell[0][2];
        Nc0=(int)round((c[2]*A->rUnitCell[0][0]-c[0]*A->rUnitCell[0][2])/V);
        Na0=(int)round((c[0]*A->rUnitCell[2][2]-c[2]*A->rUnitCell[2][0])/V);
      }
    }
    else
    { // a and b finite
      if(Norm(A->rUnitCell[2])==INFINITY)
      { // c inf a and b finite
        Nc=0;
        Nc0=0;
        V=Dot(A->rUnitCell[0],A->rUnitCell[1]);
        Na=(int)(floor(fabs(r*Norm(A->rUnitCell[1])/V))+1);
        Nb=(int)(floor(fabs(r*Norm(A->rUnitCell[0])/V))+1);
        V=A->rUnitCell[0][0]*A->rUnitCell[1][1]-A->rUnitCell[0][1]*A->rUnitCell[1][0];
        Na0=(int)round((c[0]*A->rUnitCell[1][1]-c[1]*A->rUnitCell[1][0])/V);
        Nb0=(int)round((c[1]*A->rUnitCell[0][0]-c[0]*A->rUnitCell[0][1])/V);
      }
      else
      { // a, b and c finite
        V=ParallelalPipedVolume(A->rUnitCell[0], A->rUnitCell[1], A->rUnitCell[2]);
        Cross(A->rUnitCell[1], A->rUnitCell[2], O);
        Na=(int)floor(fabs(r*Norm(O)/V))+1;
        Na0=(int)round(Dot(c,O)/V);
        Cross(A->rUnitCell[2], A->rUnitCell[0], O);
        Nb=(int)floor(fabs(r*Norm(O)/V))+1;
        Nb0=(int)round(Dot(c,O)/V);
        Cross(A->rUnitCell[0], A->rUnitCell[1], O);
        Nc=(int)floor(fabs(r*Norm(O)/V))+1;
        Nc0=(int)round(Dot(c,O)/V);
      }
    }
  }
  // set unit cells for open cluster
  B->rUnitCell[0][0]=INFINITY; B->rUnitCell[0][1]=0; B->rUnitCell[0][2]=0;
  B->rUnitCell[1][0]=0; B->rUnitCell[1][1]=INFINITY; B->rUnitCell[1][2]=0;
  B->rUnitCell[2][0]=0; B->rUnitCell[2][1]=0; B->rUnitCell[2][2]=INFINITY;
  B->gUnitCell[0][0]=0; B->gUnitCell[0][1]=0; B->gUnitCell[0][2]=0;
  B->gUnitCell[1][0]=0; B->gUnitCell[1][1]=0; B->gUnitCell[1][2]=0;
  B->gUnitCell[2][0]=0; B->gUnitCell[2][1]=0; B->gUnitCell[2][2]=0;
  B->Units=A->Units;
  B->NAtoms=0;
  B->Atom=(AtomType*)malloc(64*sizeof(AtomType));
  for(ia=Na0-Na;ia<=Na0+Na;ia++)
  {
    for(ib=Nb0-Nb;ib<=Nb0+Nb;ib++)
    {
      for(ic=Nc0-Nc;ic<=Nc0+Nc;ic++)
      {
        CellPos[0] = ((ia==0)?0:ia*A->rUnitCell[0][0]) + ((ib==0)?0:ib*A->rUnitCell[1][0]) + ((ic==0)?0:ic*A->rUnitCell[2][0]);
        CellPos[1] = ((ia==0)?0:ia*A->rUnitCell[0][1]) + ((ib==0)?0:ib*A->rUnitCell[1][1]) + ((ic==0)?0:ic*A->rUnitCell[2][1]);
        CellPos[2] = ((ia==0)?0:ia*A->rUnitCell[0][2]) + ((ib==0)?0:ib*A->rUnitCell[1][2]) + ((ic==0)?0:ic*A->rUnitCell[2][2]);
        for(iatom=0;iatom<A->NAtoms;iatom++)
        {
          if((A->Units&8)==0)
          {
            AtomPos[0]=CellPos[0]-c[0]+A->Atom[iatom].Position[0];
            AtomPos[1]=CellPos[1]-c[1]+A->Atom[iatom].Position[1];
            AtomPos[2]=CellPos[2]-c[2]+A->Atom[iatom].Position[2];
          }
          else
          {
            AtomPos[0]=CellPos[0]-c[0] + ((A->Atom[iatom].Position[0]==0)?0:A->Atom[iatom].Position[0]*A->rUnitCell[0][0]) + ((A->Atom[iatom].Position[1]==0)?0:A->Atom[iatom].Position[1]*A->rUnitCell[1][0]) + ((A->Atom[iatom].Position[2]==0)?0:A->Atom[iatom].Position[2]*A->rUnitCell[2][0]);
            AtomPos[1]=CellPos[1]-c[1] + ((A->Atom[iatom].Position[0]==0)?0:A->Atom[iatom].Position[0]*A->rUnitCell[0][1]) + ((A->Atom[iatom].Position[1]==0)?0:A->Atom[iatom].Position[1]*A->rUnitCell[1][1]) + ((A->Atom[iatom].Position[2]==0)?0:A->Atom[iatom].Position[2]*A->rUnitCell[2][1]);
            AtomPos[2]=CellPos[2]-c[2] + ((A->Atom[iatom].Position[0]==0)?0:A->Atom[iatom].Position[0]*A->rUnitCell[0][2]) + ((A->Atom[iatom].Position[1]==0)?0:A->Atom[iatom].Position[1]*A->rUnitCell[1][2]) + ((A->Atom[iatom].Position[2]==0)?0:A->Atom[iatom].Position[2]*A->rUnitCell[2][2]);
          }
          if(Dot(AtomPos,AtomPos)<rSqr)
          {
            B->NAtoms++;
            if((B->NAtoms&63)==0) B->Atom=(AtomType*)realloc(B->Atom,(B->NAtoms+64)*sizeof(AtomType));
            CopyAtomType(&(A->Atom[iatom]), &(B->Atom[B->NAtoms-1]));
            /*
            if((A->Units&8)==0)
            {
              AtomPos[0]=CellPos[0]+A->Atom[iatom].Position[0];
              AtomPos[1]=CellPos[1]+A->Atom[iatom].Position[1];
              AtomPos[2]=CellPos[2]+A->Atom[iatom].Position[2];
            }
            else
            {
              AtomPos[0]=CellPos[0] + ((A->Atom[iatom].Position[0]==0)?0:A->Atom[iatom].Position[0]*A->rUnitCell[0][0]) + ((A->Atom[iatom].Position[1]==0)?0:A->Atom[iatom].Position[1]*A->rUnitCell[1][0]) + ((A->Atom[iatom].Position[2]==0)?0:A->Atom[iatom].Position[2]*A->rUnitCell[2][0]);
              AtomPos[1]=CellPos[1] + ((A->Atom[iatom].Position[0]==0)?0:A->Atom[iatom].Position[0]*A->rUnitCell[0][1]) + ((A->Atom[iatom].Position[1]==0)?0:A->Atom[iatom].Position[1]*A->rUnitCell[1][1]) + ((A->Atom[iatom].Position[2]==0)?0:A->Atom[iatom].Position[2]*A->rUnitCell[2][1]);
              AtomPos[2]=CellPos[2] + ((A->Atom[iatom].Position[0]==0)?0:A->Atom[iatom].Position[0]*A->rUnitCell[0][2]) + ((A->Atom[iatom].Position[1]==0)?0:A->Atom[iatom].Position[1]*A->rUnitCell[1][2]) + ((A->Atom[iatom].Position[2]==0)?0:A->Atom[iatom].Position[2]*A->rUnitCell[2][2]);
            }
             */
            B->Atom[B->NAtoms-1].Position[0]=AtomPos[0]+c[0];
            B->Atom[B->NAtoms-1].Position[1]=AtomPos[1]+c[1];
            B->Atom[B->NAtoms-1].Position[2]=AtomPos[2]+c[2];
          }
        }
      }
    }
  }
  return 0;
}

int artificialOrbitalMultiplicity(const string &orb)
{
	if (orb.at(orb.length() - 1) != 'a')
		return -1;
	int firstNumberChar = 0;
	while (orb.find('_', firstNumberChar) != (int) string::npos)
		firstNumberChar = (int)orb.find('_', firstNumberChar) + 1;
	string multString = orb.substr(firstNumberChar, orb.length() - firstNumberChar - 1);
//	Here's what I would like to do:
//	try
//	{
//		return std::stoi(multString);
//	} catch (const std::invalid_argument&)
//	{
//		cout << "Warning: Multiplicity " << multString << " could not be parsed as a number, will return -1." << endl;
//		return -1;
//	}

//	And here's what I have to do, because stoi is a C++11 function:
	return (int) strtol(multString.c_str(), NULL, 10);
}

int nonRelAngularMomentumZNumbersWithSpin(const string &orb)
{
	if (orb.length() < 1)
		return -1;
	char lastChar = tolower(orb.at(orb.length() - 1));
	switch (lastChar)
	{
	case 's':
		return 2 * (2 * 0 + 1);
	case 'p':
		return 2 * (2 * 1 + 1);
	case 'd':
		return 2 * (2 * 2 + 1);
	case 'f':
		return 2 * (2 * 3 + 1);
	case 'g':
		return 2 * (2 * 4 + 1);
	case 'h':
		return 2 * (2 * 5 + 1);
	default:
		return -1;
	}
}

int relAngularMomentumZNumbers(const string &orb, int * significantCharPos)
{
	char c;
	*significantCharPos = (int) orb.length() - 1;
	c = orb.at(*significantCharPos);
	while (c == ' ' || c == '\t')
	{
		*significantCharPos = *significantCharPos - 1;
		if (*significantCharPos < 0)
			return -1;
		c = orb.at(*significantCharPos);
	}
	if (c != '2')
	{
		return -1;
	} else
	{
		*significantCharPos = *significantCharPos - 1;
		if (*significantCharPos < 0)
			return -1;
		c = orb.at(*significantCharPos);
	}
	while (c == ' ' || c == '\t' || c == '/')
	{
		*significantCharPos = *significantCharPos - 1;
		if (*significantCharPos < 0)
			return -1;
		c = orb.at(*significantCharPos);
	}
	switch (c)
	{
	case '1':
		return 1 + 1;
	case '3':
		return 3 + 1;
	case '5':
		return 5 + 1;
	case '7':
		return 7 + 1;
	case '9':
		return 9 + 1;
	default:
		return -1;
	}
}

int relAngularMomentumZNumbers(const string &orb)
{
	int significantCharPos;
	return relAngularMomentumZNumbers(orb, &significantCharPos);
}

int kappa(const string &orb)
{
	int jzNumbers, lzNumbers, significantCharPos;
	jzNumbers = relAngularMomentumZNumbers(orb, &significantCharPos);
	if (jzNumbers < 0)
	{
		cout << "\nWarning: Key \"" << orb << "\" could not be interpreted as a relativistic orbital, kappa is set to 0.\n";
		return 0;
	}
	string orbStripped = orb.substr(0, significantCharPos);
	lzNumbers = nonRelAngularMomentumZNumbersWithSpin(orbStripped) / 2;
	if (lzNumbers < 0)
	{
		cout << "\nWarning: Remainder of key \"" << orbStripped << "\" could not be interpreted as a non-relativistic orbital, kappa is set to 0.\n";
		return 0;
	}
	if (jzNumbers > lzNumbers)
		return -jzNumbers / 2;
	else
		return jzNumbers / 2;
}

vector<int> kappas(const vector<string> &orbs)
{
	vector<int> kappas(orbs.size());
	for (int i = 0; i < kappas.size(); ++i)
		kappas[i] = kappa(orbs[i]);
	return kappas;
}

int orbToMultiplicity(const string &orb, const bool hasSpin)
{
	int multiplicity = artificialOrbitalMultiplicity(orb);
	if (multiplicity < 0)
	{
		multiplicity = nonRelAngularMomentumZNumbersWithSpin(orb);
		if (!hasSpin)
			multiplicity = (int) multiplicity / 2;
	}
	if (multiplicity < 0)
		multiplicity = relAngularMomentumZNumbers(orb);
	if (multiplicity < 0)
		cout << "\nCould not interpret orbital \"" << orb << "\" as either relativistic or non-relativistic atomic orbital, returning " << multiplicity << endl;
	return multiplicity;
}

int orbToMultiplicity(const string &orb)
{
	return orbToMultiplicity(orb, true);
}



/*  On input TB contains a normal TB object
 *  StripedTB points to a TB object that is uninitialized and has no pointera allocated
 *  On access StripedTB contains a TB object that specifies the same Hamiltonian, albeit with all atoms in TB moved to position {0,0,0} in StripedTB and merged into one big "atom" with one "shell"
 */
uint32_t RemoveCrystalInformationFromTightBinding(TightBindingDef *TB, TightBindingDef *StripedTB)
{
  uint32_t i,j,iHop,jHop;
  StripedTB->NFermion = TB->NFermion;
  strcpy(StripedTB->Name, TB->Name);
  StripedTB->Cell=TB->Cell;
  StripedTB->Cell.NAtoms=1;
  StripedTB->Cell.Atom=(AtomType*)malloc(sizeof(AtomType));
  strcpy(StripedTB->Cell.Atom[0].Name, "All");
  StripedTB->Cell.Atom[0].AtomNumber=0;
  StripedTB->Cell.Atom[0].Position[0]=0;
  StripedTB->Cell.Atom[0].Position[1]=0;
  StripedTB->Cell.Atom[0].Position[2]=0;
  StripedTB->Cell.Atom[0].Nshell=1;
  StripedTB->Cell.Atom[0].ShellName=(char**)malloc(sizeof(char*));
  StripedTB->Cell.Atom[0].ShellName[0]=(char*)malloc(5*sizeof(char));
  strcpy(StripedTB->Cell.Atom[0].ShellName[0], "All");
  StripedTB->Cell.Atom[0].shelldegeneracy=(int*)malloc(sizeof(int));
  StripedTB->Cell.Atom[0].shelldegeneracy[0]=StripedTB->NFermion;
  StripedTB->Cell.Atom[0].OrbitalName=(char***)malloc(sizeof(char**));
  StripedTB->Cell.Atom[0].OrbitalName[0]=(char**)malloc(StripedTB->Cell.Atom[0].shelldegeneracy[0]*sizeof(char*));
  for(i=0;i<StripedTB->Cell.Atom[0].shelldegeneracy[0];i++)
  {
    StripedTB->Cell.Atom[0].OrbitalName[0][i]=(char*)malloc(12*sizeof(char));
    sprintf(StripedTB->Cell.Atom[0].OrbitalName[0][i], "%i",i);
  }
  SetNFermionIndexTightBindingDef(StripedTB);
  StripedTB->nH=0;
  StripedTB->H=NULL;
  for(i=0;i<TB->nH;i++)
  {
    for(j=0;j<StripedTB->nH;j++)
    {
      if((TB->H[i].TransVecs[0]==StripedTB->H[j].TransVecs[0])&&(TB->H[i].TransVecs[1]==StripedTB->H[j].TransVecs[1])&&(TB->H[i].TransVecs[2]==StripedTB->H[j].TransVecs[2]))
      {
        // hopping between two atoms at this distance is already in the new TB operator add here
        for(iHop=0;iHop<TB->H[i].Hop.Ni;iHop++)
        {
          for(jHop=0;jHop<TB->H[i].Hop.Nj;jHop++)
          {
            StripedTB->H[j].Hop.VRe[TB->FermionIndex[TB->H[i].Atompos1][TB->H[i].Shellpos1]+iHop][TB->FermionIndex[TB->H[i].Atompos2][TB->H[i].Shellpos2]+jHop] = TB->H[i].Hop.VRe[iHop][jHop];
            if(TB->H[i].Hop.QComplex)
            {
              StripedTB->H[j].Hop.VIm[TB->FermionIndex[TB->H[i].Atompos1][TB->H[i].Shellpos1]+iHop][TB->FermionIndex[TB->H[i].Atompos2][TB->H[i].Shellpos2]+jHop]=TB->H[i].Hop.VIm[iHop][jHop];
            }
          }
        }
        goto TBElementFound;
      }
    }
    // new element
    if((StripedTB->nH&15)==0)
    {
      StripedTB->H=(TightBindingDefHopping*)realloc(StripedTB->H, (StripedTB->nH+17)*sizeof(TightBindingDefHopping));
      if((StripedTB->H)==NULL)
      {
        printf("Error in RemoveCrystalInformationFromTightBinding\nrealloc failed for StripedTB->H\nI will crash soon\n");
        return 1;
      }
    }
    StripedTB->H[StripedTB->nH].Atompos1=0;
    StripedTB->H[StripedTB->nH].Shellpos1=0;
    StripedTB->H[StripedTB->nH].Atompos2=0;
    StripedTB->H[StripedTB->nH].Shellpos2=0;
    StripedTB->H[StripedTB->nH].TransVecs[0]=TB->H[i].TransVecs[0];
    StripedTB->H[StripedTB->nH].TransVecs[1]=TB->H[i].TransVecs[1];
    StripedTB->H[StripedTB->nH].TransVecs[2]=TB->H[i].TransVecs[2];
    // set ConnectingVector depending on Units
    if((StripedTB->Cell.Units&8)==0)
    {
      // absolute positions
      StripedTB->H[StripedTB->nH].ConnectingVector[0] = ((StripedTB->H[StripedTB->nH].TransVecs[0]==0)?0:(StripedTB->H[StripedTB->nH].TransVecs[0] * StripedTB->Cell.rUnitCell[0][0])) + ((StripedTB->H[StripedTB->nH].TransVecs[1]==0)?0:(StripedTB->H[StripedTB->nH].TransVecs[1] * StripedTB->Cell.rUnitCell[1][0])) + ((StripedTB->H[StripedTB->nH].TransVecs[2]==0)?0:(StripedTB->H[StripedTB->nH].TransVecs[2] * StripedTB->Cell.rUnitCell[2][0]));
      StripedTB->H[StripedTB->nH].ConnectingVector[0] = ((StripedTB->H[StripedTB->nH].TransVecs[0]==0)?0:(StripedTB->H[StripedTB->nH].TransVecs[0] * StripedTB->Cell.rUnitCell[0][1])) + ((StripedTB->H[StripedTB->nH].TransVecs[1]==0)?0:(StripedTB->H[StripedTB->nH].TransVecs[1] * StripedTB->Cell.rUnitCell[1][1])) + ((StripedTB->H[StripedTB->nH].TransVecs[2]==0)?0:(StripedTB->H[StripedTB->nH].TransVecs[2] * StripedTB->Cell.rUnitCell[2][1]));
      StripedTB->H[StripedTB->nH].ConnectingVector[0] = ((StripedTB->H[StripedTB->nH].TransVecs[0]==0)?0:(StripedTB->H[StripedTB->nH].TransVecs[0] * StripedTB->Cell.rUnitCell[0][2])) + ((StripedTB->H[StripedTB->nH].TransVecs[1]==0)?0:(StripedTB->H[StripedTB->nH].TransVecs[1] * StripedTB->Cell.rUnitCell[1][2])) + ((StripedTB->H[StripedTB->nH].TransVecs[2]==0)?0:(StripedTB->H[StripedTB->nH].TransVecs[2] * StripedTB->Cell.rUnitCell[2][2]));
    }
    else
    {
      // positions in terms of unit cells
      StripedTB->H[StripedTB->nH].ConnectingVector[0] = StripedTB->H[StripedTB->nH].TransVecs[0];
      StripedTB->H[StripedTB->nH].ConnectingVector[1] = StripedTB->H[StripedTB->nH].TransVecs[1];
      StripedTB->H[StripedTB->nH].ConnectingVector[2] = StripedTB->H[StripedTB->nH].TransVecs[2];
    }
    StripedTB->H[StripedTB->nH].Hop.Ni=StripedTB->NFermion;
    StripedTB->H[StripedTB->nH].Hop.Nj=StripedTB->NFermion;
    StripedTB->H[StripedTB->nH].Hop.QComplex=1;
    InitMatrix(&(StripedTB->H[StripedTB->nH].Hop));
    // hopping between two atoms this distance has not been added before
    for(iHop=0;iHop<TB->H[i].Hop.Ni;iHop++)
    {
      for(jHop=0;jHop<TB->H[i].Hop.Nj;jHop++)
      {
        StripedTB->H[StripedTB->nH].Hop.VRe[TB->FermionIndex[TB->H[i].Atompos1][TB->H[i].Shellpos1]+iHop][TB->FermionIndex[TB->H[i].Atompos2][TB->H[i].Shellpos2]+jHop] = TB->H[i].Hop.VRe[iHop][jHop];
        if(TB->H[i].Hop.QComplex)
        {
          StripedTB->H[StripedTB->nH].Hop.VIm[TB->FermionIndex[TB->H[i].Atompos1][TB->H[i].Shellpos1]+iHop][TB->FermionIndex[TB->H[i].Atompos2][TB->H[i].Shellpos2]+jHop]=TB->H[i].Hop.VIm[iHop][jHop];
        }
      }
    }
    StripedTB->nH++;
    TBElementFound:;
  }
  // remove complex parts of matrices that are real
  for(i=0;i<StripedTB->nH;i++)
  {
    MatrixChop(&(StripedTB->H[i].Hop), DBL_MIN);
  }
  return 0;
}
/*  Undo that what RemoveCrystalInformationFromTightBinding did with a TB object
 *  On entry StripedTB contains a TB object with a single atom with single shell per unit cell
 *  On entry TB contains a TB object with TB->nH=0 and TB->H=NULL
 *  TB->NFermion must be equal to StripedTB->NFermion
 *
 *  The single matrix for StripedTB is broken into pieces acording to the atomic shells in the Cell definition of TB
 */
uint32_t AddCrystalInformationFromTightBinding(TightBindingDef *StripedTB, TightBindingDef *TB)
{
  uint32_t iH,iAtom,iShell,jAtom,jShell,iHop,jHop;
  
  for(iH=0;iH<StripedTB->nH;iH++)
  {
    for(iAtom=0;iAtom<TB->Cell.NAtoms;iAtom++)
    {
      for(iShell=0;iShell<TB->Cell.Atom[iAtom].Nshell;iShell++)
      {
        for(jAtom=0;jAtom<TB->Cell.NAtoms;jAtom++)
        {
          for(jShell=0;jShell<TB->Cell.Atom[jAtom].Nshell;jShell++)
          {
            if((TB->nH&15)==0)
            {
              TB->H=(TightBindingDefHopping*)realloc(TB->H, (TB->nH+17)*sizeof(TightBindingDefHopping));
              if((TB->H)==NULL)
              {
                printf("Error in AddCrystalInformationFromTightBinding\nrealloc failed for TB->H\nI will crash soon\n");
                return 1;
              }
            }
            TB->H[TB->nH].Atompos1=iAtom;
            TB->H[TB->nH].Atompos2=jAtom;
            TB->H[TB->nH].Shellpos1=iShell;
            TB->H[TB->nH].Shellpos2=jShell;
            TB->H[TB->nH].TransVecs[0]=StripedTB->H[iH].TransVecs[0];
            TB->H[TB->nH].TransVecs[1]=StripedTB->H[iH].TransVecs[1];
            TB->H[TB->nH].TransVecs[2]=StripedTB->H[iH].TransVecs[2];
            if((TB->Cell.Units&8)==0)
            { // position in absolute units
              TB->H[TB->nH].ConnectingVector[0]=TB->Cell.Atom[jAtom].Position[0] - TB->Cell.Atom[iAtom].Position[0] + ((TB->H[TB->nH].TransVecs[0]==0)?0:TB->H[TB->nH].TransVecs[0]*TB->Cell.rUnitCell[0][0]) + ((TB->H[TB->nH].TransVecs[1]==0)?0:TB->H[TB->nH].TransVecs[1]*TB->Cell.rUnitCell[1][0]) + ((TB->H[TB->nH].TransVecs[2]==0)?0:TB->H[TB->nH].TransVecs[2]*TB->Cell.rUnitCell[2][0]) ;
              TB->H[TB->nH].ConnectingVector[1]=TB->Cell.Atom[jAtom].Position[1] - TB->Cell.Atom[iAtom].Position[1] + ((TB->H[TB->nH].TransVecs[0]==0)?0:TB->H[TB->nH].TransVecs[0]*TB->Cell.rUnitCell[0][1]) + ((TB->H[TB->nH].TransVecs[1]==0)?0:TB->H[TB->nH].TransVecs[1]*TB->Cell.rUnitCell[1][1]) + ((TB->H[TB->nH].TransVecs[2]==0)?0:TB->H[TB->nH].TransVecs[2]*TB->Cell.rUnitCell[2][1]) ;
              TB->H[TB->nH].ConnectingVector[2]=TB->Cell.Atom[jAtom].Position[2] - TB->Cell.Atom[iAtom].Position[2] + ((TB->H[TB->nH].TransVecs[0]==0)?0:TB->H[TB->nH].TransVecs[0]*TB->Cell.rUnitCell[0][2]) + ((TB->H[TB->nH].TransVecs[1]==0)?0:TB->H[TB->nH].TransVecs[1]*TB->Cell.rUnitCell[1][2]) + ((TB->H[TB->nH].TransVecs[2]==0)?0:TB->H[TB->nH].TransVecs[2]*TB->Cell.rUnitCell[2][2]) ;
            }
            else
            { // positioin in terms of lattice vectors
              TB->H[TB->nH].ConnectingVector[0]=TB->Cell.Atom[jAtom].Position[0] - TB->Cell.Atom[iAtom].Position[0] + TB->H[TB->nH].TransVecs[0];
              TB->H[TB->nH].ConnectingVector[1]=TB->Cell.Atom[jAtom].Position[1] - TB->Cell.Atom[iAtom].Position[1] + TB->H[TB->nH].TransVecs[1];
              TB->H[TB->nH].ConnectingVector[2]=TB->Cell.Atom[jAtom].Position[2] - TB->Cell.Atom[iAtom].Position[2] + TB->H[TB->nH].TransVecs[2];
            }
            TB->H[TB->nH].Hop.Ni=TB->Cell.Atom[iAtom].shelldegeneracy[iShell];
            TB->H[TB->nH].Hop.Nj=TB->Cell.Atom[jAtom].shelldegeneracy[jShell];
            TB->H[TB->nH].Hop.QComplex=StripedTB->H[iH].Hop.QComplex;
            InitMatrix(&(TB->H[TB->nH].Hop)); //TODO: set name
            for(iHop=0;iHop<TB->H[TB->nH].Hop.Ni;iHop++)
            {
              for(jHop=0;jHop<TB->H[TB->nH].Hop.Nj;jHop++)
              {
                TB->H[TB->nH].Hop.VRe[iHop][jHop]=StripedTB->H[iH].Hop.VRe[TB->FermionIndex[TB->H[TB->nH].Atompos1][TB->H[TB->nH].Shellpos1]+iHop][TB->FermionIndex[TB->H[TB->nH].Atompos2][TB->H[TB->nH].Shellpos2]+jHop];
              }
            }
            if(TB->H[TB->nH].Hop.QComplex)
            {
              for(iHop=0;iHop<TB->H[TB->nH].Hop.Ni;iHop++)
              {
                for(jHop=0;jHop<TB->H[TB->nH].Hop.Nj;jHop++)
                {
                  TB->H[TB->nH].Hop.VIm[iHop][jHop]=StripedTB->H[iH].Hop.VIm[TB->FermionIndex[TB->H[TB->nH].Atompos1][TB->H[TB->nH].Shellpos1]+iHop][TB->FermionIndex[TB->H[TB->nH].Atompos2][TB->H[TB->nH].Shellpos2]+jHop];
                }
              }
            }
            
            TB->nH++;
          }
        }
      }
    }
  }
  SetNFermionIndexTightBindingDef(StripedTB);
  return 0;
}
/*  Rotate the basis of a tight binding Hamiltonian
 *  TB is a fully defined TB opject with only one atom with one shell per unit cell
 *  RTB has RTB->Cell, RTB->NFermion, and RTB->NFermionIndex defined
 *  RTB has only one atom with one shell per unit cell
 *  Rotmat is a matrix of dimension Rotmat.Ni=RTB->NFermion and Rotmat.Nj=TB->NFermion
 *  Returned is a tight binding Operator on a transformed basis
 */
uint32_t RotateStripedTightBinding(TightBindingDef * TB, CompactMatrixType * rotmat, TightBindingDef * RTB)
{
  double AbsAverage;
  uint32_t iH;
  if(TB->Cell.NAtoms!=1)
  {
    printf("Error in RotateStripedTightBinding\nOne one atom per unit cell allowed if one rotates the basis");
    return 1;
  }
  if(TB->Cell.Atom[0].Nshell!=1)
  {
    printf("Error in RotateStripedTightBinding\nOne one shell for the single atom in the unit cell allowed if one rotates the basis");
    return 1;
  }
  if(TB->NFermion!=rotmat->Nj)
  {
    printf("Error in RotateStripedTightBinding\nSize of rotation matrix (number of columns, %i) not equal to the number of Fermions in the original tight binding Hamiltonian (%i)\n",rotmat->Nj,TB->NFermion);
    return 1;
  }
  
  if(RTB->NFermion!=rotmat->Ni)
  {
    printf("Error in RotateStripedTightBinding\nSize of rotation matrix (number of rows, %i) not equal to the number of Fermions in the rotated tight binding Hamiltonian (%i)\n",rotmat->Ni,RTB->NFermion);
    return 1;
  }
  RTB->nH=TB->nH;
  RTB->H=(TightBindingDefHopping*)malloc(RTB->nH*sizeof(TightBindingDefHopping));
  for(iH=0;iH<TB->nH;iH++)
  {
    RTB->H[iH].Atompos1=TB->H[iH].Atompos1;
    RTB->H[iH].Shellpos1=TB->H[iH].Shellpos1;
    RTB->H[iH].Atompos2=TB->H[iH].Atompos2;
    RTB->H[iH].Shellpos2=TB->H[iH].Shellpos2;
    RTB->H[iH].TransVecs[0]=TB->H[iH].TransVecs[0];
    RTB->H[iH].TransVecs[1]=TB->H[iH].TransVecs[1];
    RTB->H[iH].TransVecs[2]=TB->H[iH].TransVecs[2];
    RTB->H[iH].ConnectingVector[0]=TB->H[iH].ConnectingVector[0];
    RTB->H[iH].ConnectingVector[1]=TB->H[iH].ConnectingVector[1];
    RTB->H[iH].ConnectingVector[2]=TB->H[iH].ConnectingVector[2];
    MatrixTransformByCompactMatrix(&(RTB->H[iH].Hop), TB->H[iH].Hop, rotmat);
    MatrixAbsAverage(RTB->H[iH].Hop,&AbsAverage);
    MatrixChop(&(RTB->H[iH].Hop), AbsAverage*DBL_EPSILON*RTB->H[iH].Hop.Ni*10);
    
  }
  SetNFermionIndexTightBindingDef(RTB);
  return 0;
}
/*  Rotate the basis of a tight binding Hamiltonian
 *  TB is a fully defined TB opject
 *  RTB has RTB->Cell, RTB->NFermion, and RTB->NFermionIndex defined
 *  Rotmat is a matrix of dimension Rotmat.Ni=RTB->NFermion and Rotmat.Nj=TB->NFermion
 *  Returned is a tight binding Operator on a transformed basis
 */
uint32_t RotateTightBinding(TightBindingDef * TB, CompactMatrixType * rotmat, TightBindingDef * RTB)
{
  TightBindingDef StripedTB, StripedRTB;
  RemoveCrystalInformationFromTightBinding(TB, &StripedTB);
  RemoveCrystalInformationFromTightBinding(RTB,&StripedRTB);
  RotateStripedTightBinding(&StripedTB,rotmat,&StripedRTB);
  AddCrystalInformationFromTightBinding(&StripedRTB,RTB);
  FreeTightBindingDef(&StripedRTB);
  FreeTightBindingDef(&StripedTB);
  return 0;
}

vector<string> groupedRelOrbitals(const vector<string> &in)
{
	vector < string > out;
	vector<int> blackList;
	for (int i = 0; i < in.size(); ++i)
	{
		if (!isElementOf(i, blackList))
		{
			out.push_back(in[i]);
			if (isRelOrbital(in[i]) && (kappa(in[i]) != -1))
			{
				for (int j = i + 1; j < in.size(); ++j)
				{
					if (QOrbsBelongTogether(in[i], in[j]))
					{
						out.push_back(in[j]);
						blackList.push_back(j);
						if (orbToMultiplicity(out[out.size()-1]) < orbToMultiplicity(out[out.size()-2]) )
							std::swap(out[out.size()-1],out[out.size()-2]);
						break;
					}
				}
			}
		}
	}
	return out;
}

bool isRelOrbital(const string &orb)
{
	return (relAngularMomentumZNumbers(orb) > 0);
}

bool QOrbsBelongTogether(string orb1, string orb2)
{
	if (isRelOrbital(orb1))
		orb1 = cutRelPart(orb1);
	if (isRelOrbital(orb2))
		orb2 = cutRelPart(orb2);
	return orb1 == orb2;
}

string cutRelPart(const string &orb)
{
	int significantCharPos;
	relAngularMomentumZNumbers(orb,&significantCharPos);
	return orb.substr(0,significantCharPos);
}

vector<string> relToNonRelOrbs(const vector<string> &in)
{
	vector < string > out;
	for (int i = 0; i < in.size(); ++i)
	{
		if (isRelOrbital(in[i]))
		{
			if (kappa(in[i]) == -1)
			{
				out.push_back(cutRelPart(in[i]));
			} else if (i + 1 < in.size() && QOrbsBelongTogether(in[i], in[i + 1]))
			{
				out.push_back(cutRelPart(in[i]));
				++i;
			} else
			{
				cout << "Error: Need a pair of relativistic orbitals to convert them to a non-relativistic one, unless l=0." << endl;
				vector < string > errorVec;
				return errorVec;
			}
		} else
		{
			out.push_back(in[i]);
		}
	}
	return out;
}

vector<string> nonRelToRelOrbs(const vector<string> &in)
{
	vector < string > out;
	for (int i = 0; i < in.size(); ++i)
	{
		if (isRelOrbital(in[i]))
		{
			out.push_back(in[i]);
		} else
		{
			int mult = orbToMultiplicity(in[i]);
			if (mult == 2)
			{
				string relOrb = in[i] + "1/2";
				out.push_back(relOrb);
			} else
			{
				string relOrb1 = in[i] + to_string(mult / 2) + "/2";
				string relOrb2 = in[i] + to_string(mult / 2 - 2) + "/2";
				out.push_back(relOrb1);
				out.push_back(relOrb2);
			}
		}
	}
	return out;
}

bool isOrbElementOf(const string &orb, vector<string> set)
{
	if (isRelOrbital(orb))
		set = nonRelToRelOrbs(set);
	else
		set = relToNonRelOrbs(set);
	return isElementOf(orb, set);
}

// input e.g. {"2p32","2s12","2p12"}, {"2s","2p"}
MatrixType OrbitalRotationMatrix(vector<string> orbsOld, vector<string> orbsNew)
{
	if (hasIdenticalEntries(orbsOld))
	{
		cout << "Error: RotationMatrix called with identical entries in the list of old orbitals." << endl;
		return nullMatrix();
	}
	if (hasIdenticalEntries(orbsNew))
	{
		cout << "Error: RotationMatrix called with identical entries in the list of new orbitals." << endl;
		return nullMatrix();
	}
	int NFOld = NumberOfFermionicStates(orbsOld);
	int NFNew = NumberOfFermionicStates(orbsNew);
	int NF = 0;
	if (NFOld < NFNew)
	{
		for (uint i = 0; i < (uint) orbsNew.size(); ++i)
			if (!isOrbElementOf(orbsNew[i], orbsOld))
				orbsOld.push_back(orbsNew[i]);
		NF = NFNew;
	} else if (NFOld > NFNew)
	{
		for (uint i = 0; i < (uint) orbsOld.size(); ++i)
			if (!isOrbElementOf(orbsOld[i], orbsNew))
				orbsNew.push_back(orbsOld[i]);
		NF = NFOld;
	} else
	{
		NF = NFOld;
	}
	if (NumberOfFermionicStates(orbsOld) != NumberOfFermionicStates(orbsNew))
	{
		cout << "Error: RotationMatrix failed to make a connection between old and new orbitals." << endl;
		return nullMatrix();
	}
	vector < string > orbsTmp;
//	{"2p12","2s12","2p32"} -> {"2p12","2p32","2s12"}
	orbsTmp = orbsOld;
	orbsOld = groupedRelOrbitals(orbsOld);
	MatrixType TGroupingRight = orbitalPermutationMatrix(orbsTmp, orbsOld);
	orbsTmp = orbsNew;
	orbsNew = groupedRelOrbitals(orbsNew);
	MatrixType TGroupingLeft = orbitalPermutationMatrix(orbsNew, orbsTmp);
	if (isNullMatrix(TGroupingRight) || isNullMatrix(TGroupingLeft))
	{
		cout << "Error in RotationMatrix: Grouping of orbitals failed." << endl;
		return nullMatrix();
	}
	MatrixType TRelNonRelLeft, TRelNonRelRight;
//	{"2p12","2p32","2s12"} -> {"2p","2s"}
	if (QChangeFromRelToNonRel(orbsOld, orbsNew))
	{
		orbsTmp = orbsOld;
		orbsOld = relToNonRelOrbs(orbsOld);
		TRelNonRelRight = RelToNonRelRotationMatrix(orbsTmp, orbsOld);
		orbsTmp = orbsNew;
		orbsNew = relToNonRelOrbs(orbsNew);
		TRelNonRelLeft = NonRelToRelRotationMatrix(orbsNew, orbsTmp);
		if (isNullMatrix(TRelNonRelLeft) || isNullMatrix(TRelNonRelRight))
		{
			cout << "Error in RotationMatrix: Rotation from relativistic to non-relativistic orbitals or vice versa failed." << endl;
			return nullMatrix();
		}
	} else
	{
		TRelNonRelLeft = identityMatrix(NF);
		TRelNonRelRight = identityMatrix(NF);
	}
//	{"2p","2s"} -> {"2s","2p"}
	MatrixType TPerm = orbitalPermutationMatrix(orbsOld, orbsNew);
	if (isNullMatrix(TPerm))
	{
		cout << "Error in RotationMatrix: Permutation of orbitals failed." << endl;
		return nullMatrix();
	}
	MatrixType rotMat = identityMatrix(NF);
	MatrixMultiply(&rotMat, TGroupingRight, rotMat);
	MatrixMultiply(&rotMat, TRelNonRelRight, rotMat);
	MatrixMultiply(&rotMat, TPerm, rotMat);
	MatrixMultiply(&rotMat, TRelNonRelLeft, rotMat);
	MatrixMultiply(&rotMat, TGroupingLeft, rotMat);
	MatrixChop(&rotMat, 10 * DBL_EPSILON);
	rotMat = Submatrix(rotMat, NFNew, NFOld);
	return rotMat;
}

//TODO: Anyone have a better idea where to put this?
template<typename T>
bool hasIdenticalEntries(const vector<T> &vec)
{
	for (int i = 0; i < (int)vec.size() - 1; ++i)
		for (int j = i + 1; j < (int)vec.size(); ++j)
			if (vec[i] == vec[j])
				return true;
	return false;
}

//TODO: Anyone have a better idea where to put this?
template<typename T>
vector<int> permutationVector(const vector<T> &oldVec, const vector<T> &newVec)
{
	if (oldVec.size() != newVec.size())
	{
		cout << "Error: PermutationVector called with two vectors of different size." << endl;
		vector<int> errorVec(0);
		return errorVec;
	}
	vector<int> permVec(newVec.size());
	for (int i = 0; i < newVec.size(); ++i)
		for (int j = 0; j < oldVec.size(); ++j)
		{
			if (newVec[i] == oldVec[j])
			{
				permVec[i] = j;
				break;
			}
		}
	return permVec;
}

MatrixType YToKMatrix(const int multiplicity, const bool addSpin)
{
	MatrixType out;
	out.QComplex = true;
	out.Ni = multiplicity / 2;
	out.Nj = multiplicity / 2;
	InitMatrix(&out);
	if (multiplicity == 2)
	{
//		s-shell
		out.VRe[0][0] = 1;
	}
  else if (multiplicity == 6)
	{
//		p-shell
		out.VRe[0][0] = 1 / sqrt(2.);
		out.VRe[0][2] = -1 / sqrt(2.);
		out.VIm[1][0] = 1 / sqrt(2.);
		out.VIm[1][2] = 1 / sqrt(2.);
		out.VRe[2][1] = 1;
	}
  else if (multiplicity == 10)
	{
//		d-shell
		out.VRe[0][0] = 1 / sqrt(2.);
		out.VRe[0][4] = 1 / sqrt(2.);
		out.VRe[1][2] = 1;
		out.VIm[2][1] = 1 / sqrt(2.);
		out.VIm[2][3] = 1 / sqrt(2.);
		out.VRe[3][1] = 1 / sqrt(2.);
		out.VRe[3][3] = -1 / sqrt(2.);
		out.VIm[4][0] = 1 / sqrt(2.);
		out.VIm[4][4] = -1 / sqrt(2.);
	}
  else if (multiplicity == 14)
	{
//		f-shell
		out.VIm[0][1] = 1/sqrt(2.);
		out.VIm[0][5] = -1/sqrt(2.);
		out.VRe[1][0] = sqrt(5.)/4;
		out.VRe[1][2] = -sqrt(3.)/4;
		out.VRe[1][4] = sqrt(3.)/4;
		out.VRe[1][6] = -sqrt(5.)/4;
		out.VIm[2][0] = -sqrt(5.)/4;
		out.VIm[2][2] = -sqrt(3.)/4;
		out.VIm[2][4] = -sqrt(3.)/4;
		out.VIm[2][6] = -sqrt(5.)/4;
		out.VRe[3][3] = 1;
		out.VRe[4][0] = -sqrt(3.)/4;
		out.VRe[4][2] = -sqrt(5.)/4;
		out.VRe[4][4] = sqrt(5.)/4;
		out.VRe[4][6] = sqrt(3.)/4;
		out.VIm[5][0] = -sqrt(3.)/4;
		out.VIm[5][2] = sqrt(5.)/4;
		out.VIm[5][4] = sqrt(5.)/4;
		out.VIm[5][6] = -sqrt(3.)/4;
		out.VRe[6][1] = 1/sqrt(2.);
		out.VRe[6][5] = 1/sqrt(2.);
	}
  else if ((multiplicity - 2) % 4 == 0)
	{
		cout << "YToKMatrix not yet implemented for an atomic shell of multiplicity " << multiplicity << ", sorry!" << endl;
		return nullMatrix();
	} else
	{
		cout << "YToKMatrix only defined for multiplicity (2*l+1), with integer l, not for multiplicity " << multiplicity << "!" << endl;
		return nullMatrix();
	}
	if (addSpin)
		out = MatrixEnlarge(out, 2);
	return out;
}

//TODO: should this function not have l as first argument. (Multiplicity is 4l+2 ???)
MatrixType YToZMatrix(const int multiplicity, const bool addSpin)
{
	MatrixType out;
  int l,m;
	out.QComplex = true;
	out.Ni = multiplicity / 2;
	out.Nj = multiplicity / 2;
	InitMatrix(&out);
  l=(multiplicity-2)/4;
/*
	if (multiplicity == 2)
	{
//		s-shell
		out.VRe[0][0] = 1;
	}
  else if (multiplicity == 6)
	{
//		p-shell
		out.VIm[0][0] = 1 / sqrt(2.);
		out.VIm[0][2] = 1 / sqrt(2.);
		out.VRe[1][1] = 1;
		out.VRe[2][0] = 1 / sqrt(2.);
		out.VRe[2][2] = -1 / sqrt(2.);
	}
  else if (multiplicity == 10)
	{
//		d-shell
		out.VIm[0][0] = 1 / sqrt(2.);
		out.VIm[0][4] = -1 / sqrt(2.);
		out.VIm[1][1] = 1 / sqrt(2.);
		out.VIm[1][3] = 1 / sqrt(2.);
		out.VRe[2][2] = 1;
		out.VRe[3][1] = 1 / sqrt(2.);
		out.VRe[3][3] = -1 / sqrt(2.);
		out.VRe[4][0] = 1 / sqrt(2.);
		out.VRe[4][4] = 1 / sqrt(2.);
	}
 */
  if ((multiplicity - 2) % 4 == 0)
	{
    out.VRe[l][l]=1;
    for(m=1;m<=l;m++)
    {
      out.VIm[l-m][l-m] = sqrt(0.5);
      out.VIm[l-m][l+m] =-(1-2*(m%2))*sqrt(0.5);
      out.VRe[l+m][l+m] = (1-2*(m%2))*sqrt(0.5);
      out.VRe[l+m][l-m] = sqrt(0.5);
    }
	}
  else
	{
		cout << "YToZMatrix only defined for multiplicity (2*l+1), with integer l, not for multiplicity " << multiplicity << "!" << endl;
		return nullMatrix();
	}

	if (addSpin)
		out = MatrixEnlarge(out, 2);
	return out;
}

MatrixType YToJJzMatrix(const int size)
{
	return YToJJzMatrix(size, true);
}

//TODO: should size not be replaced by l to make all consistent??
MatrixType YToJJzMatrix(const int size, const bool addSpin)
{
	if (!addSpin)
	{
//		The bool is needed for compability with other code.
		cout << "YToJJzMatrix called with addSpin = false. Since spin is obligatory this request will be ignored." << endl;
	}
  int l,m;
  double j,jz;
	MatrixType out;
	out.QComplex = false;
	out.Ni = size;
	out.Nj = size;
  l=(size-2)/4;
	InitMatrix(&out);
  /*
	if (size == 2)
	{
//		s-shell
		out.VRe[0][0] = 1;
		out.VRe[1][1] = 1;
	}
  else if (size == 6)
	{
//		p-shell
		out.VRe[0][1] = -sqrt(2. / 3.);
		out.VRe[0][2] = 1 / sqrt(3);
		out.VRe[1][3] = -1 / sqrt(3);
		out.VRe[1][4] = sqrt(2. / 3.);
		out.VRe[2][0] = 1;
		out.VRe[3][1] = 1 / sqrt(3);
		out.VRe[3][2] = sqrt(2. / 3.);
		out.VRe[4][3] = sqrt(2. / 3.);
		out.VRe[4][4] = 1 / sqrt(3);
		out.VRe[5][5] = 1;
	}
  else if (size == 10)
	{
//		d-shell
		out.VRe[0][1] = -sqrt(4. / 5.);
		out.VRe[0][2] = 1 / sqrt(5);
		out.VRe[1][3] = -sqrt(3. / 5.);
		out.VRe[1][4] = sqrt(2. / 5.);
		out.VRe[2][5] = -sqrt(2. / 5.);
		out.VRe[2][6] = sqrt(3. / 5.);
		out.VRe[3][7] = -1 / sqrt(5);
		out.VRe[3][8] = sqrt(4. / 5.);
		out.VRe[4][0] = 1;
		out.VRe[5][1] = 1 / sqrt(5);
		out.VRe[5][2] = sqrt(4. / 5.);
		out.VRe[6][3] = sqrt(2. / 5.);
		out.VRe[6][4] = sqrt(3. / 5.);
		out.VRe[7][5] = sqrt(3. / 5.);
		out.VRe[7][6] = sqrt(2. / 5.);
		out.VRe[8][7] = sqrt(4. / 5.);
		out.VRe[8][8] = 1 / sqrt(5);
		out.VRe[9][9] = 1;
	}
   */
  if ((size - 2) % 4 == 0)
	{
    // j=l-1/2
    j=(double)l-0.5;
    for(m=0;m<2*l;m++)
    {
      jz=(double)m-j;
      out.VRe[m][2*m+1] =-sqrt((1-jz/(j+1))*0.5); // spin up
      out.VRe[m][2*m+2] = sqrt((1+jz/(j+1))*0.5); // spin down
    }
    // j=l+1/2
    j=(double)l+0.5;
    out.VRe[2*l][0]=1;
    out.VRe[4*l+1][4*l+1]=1;
    for(m=1;m<2*l+1;m++)
    {
      jz=(double)m-j;
      out.VRe[m+2*l][2*m-1] = sqrt((1+jz/(j))*0.5); // spin up
      out.VRe[m+2*l][2*m]   = sqrt((1-jz/(j))*0.5); // spin down
    }
	}
  else
	{
		cout << "YToJJzMatrix only defined for multiplicity 2*(2*l+1), with integer l, not for multiplicity " << size << "!" << endl;
		return nullMatrix();
	}
	return out;
}

MatrixType orbitalPermutationMatrix(const vector<string> &oldVec, const vector<string> &newVec, bool addSpin)
{
	MatrixType out;
	uint spinFac = addSpin ? 1 : 2;
	vector<int> permVec = permutationVector(oldVec, newVec);
	if (permVec.size() == 0)
		return nullMatrix();
	vector < vector<MatrixType> > matMat = nullMatrixMatrix((uint32_t) newVec.size(), (uint32_t) oldVec.size());
	for (int i = 0; i < permVec.size(); ++i)
	{
		MatrixType identityMat = identityMatrix(orbToMultiplicity(newVec[i]) / spinFac);
		matMat[i][permVec[i]] = identityMat;
	}
	out = MatrixFlatten(matMat);
	FreeMatrixMatrix (matMat);
	return out;
}

MatrixType orbitalPermutationMatrix(const vector<string> &oldVec, const vector<string> &newVec)
{
	return orbitalPermutationMatrix(oldVec, newVec, true);
}


bool QChangeFromRelToNonRel(const vector<string> &orbs1, const vector<string> &orbs2)
{
	for (int i = 0; i < orbs1.size(); ++i)
	{
		for (int j = 0; j < orbs2.size(); ++j)
		{
			if (orbs1[i] == orbs2[j])
			{
				break;
			} else if (QOrbsBelongTogether(orbs1[i], orbs2[j]))
			{
				if (isRelOrbital(orbs1[i]) && !isRelOrbital(orbs2[j]))
				{
					return true;
				} else if (!isRelOrbital(orbs1[i]) && isRelOrbital(orbs2[j]))
				{
					return true;
				}
			}
		}
	}
	return false;
}

int NumberOfFermionicStates(const vector<string> &orbs)
{
	int NF = 0;
	int mult;
	for (int i = 0; i < orbs.size(); ++i)
	{
		mult = orbToMultiplicity(orbs[i]);
		if (mult > 0)
		{
			NF += mult;
		} else
		{
			cout << "Error during calculation of NF: Orbital " << orbs[i] << " could not be interpreted as an atomic orbital.\n";
			return -1;
		}
	}
	return NF;
}

MatrixType NonRelToRelRotationMatrix(const vector<string> &nonRelOrbs, const vector<string> &mixedOrbs)
{
	vector < vector<MatrixType> > blockMat = nullMatrixMatrix((uint32_t) nonRelOrbs.size(), (uint32_t) nonRelOrbs.size());
	int offset = 0;
	for (int i = 0; i < nonRelOrbs.size(); ++i)
	{
		string oldOrb = nonRelOrbs[i];
		string newOrb = mixedOrbs[i + offset];
		if (i + offset >= mixedOrbs.size())
			break;
		if (isRelOrbital(oldOrb))
		{
			cout << "Error in NonRelToRelRotationMatrix: Encountered relativistic orbital \"" << oldOrb << "\" in list of non-relativistic orbitals." << endl;
			return nullMatrix();
		}
		if (QOrbsBelongTogether(oldOrb, newOrb))
		{
			if (isRelOrbital(newOrb))
			{
				if (kappa(newOrb) == -1)
				{
					blockMat[i][i] = YToJJzMatrix(2);
				} else
				{
					if (i + offset + 1 >= mixedOrbs.size())
					{
						cout << "Error in NonRelToRelRotationMatrix: List of new orbitals is missing the second entry belonging to \"" << oldOrb << "\"." << endl;
						return nullMatrix();
					}
					string newOrb2 = mixedOrbs[i + offset + 1];
					if (QOrbsBelongTogether(oldOrb, newOrb2))
					{
						blockMat[i][i] = YToJJzMatrix(orbToMultiplicity(oldOrb));
						++offset;
					} else
					{
						cout << "Error in NonRelToRelRotationMatrix: Could not connect \"" << oldOrb << "\" with \"" << newOrb2 << "\"." << endl;
						return nullMatrix();
					}
				}
			} else
			{
				blockMat[i][i] = identityMatrix(orbToMultiplicity(oldOrb));
			}
		} else
		{
			cout << "Error in NonRelToRelRotationMatrix: Could not connect \"" << oldOrb << "\" with \"" << newOrb << "\"." << endl;
			return nullMatrix();
		}
	}
	return MatrixFlatten(blockMat, true);
}

MatrixType RelToNonRelRotationMatrix(const vector<string> &mixedOrbs, const vector<string> &nonRelOrbs)
{
	MatrixType tmp = NonRelToRelRotationMatrix(nonRelOrbs, mixedOrbs);
	MatrixTranspose(&tmp);
	return tmp;
}

MatrixType LargeRotationMatrix(const vector<int> &multiplicities, MatrixType (*mat)(int, bool), bool addSpin)
{
	vector < vector<MatrixType> > fullM = nullMatrixMatrix((uint32_t) multiplicities.size(), (uint32_t) multiplicities.size());
	for (int i = 0; i < multiplicities.size(); ++i)
		fullM[i][i] = mat(multiplicities[i], addSpin);
	return MatrixFlatten(fullM);
}

MatrixType LargeRotationMatrix(const vector<string> &orbitals, MatrixType (*mat)(int, bool), bool addSpin)
{
	vector<int> multiplicities(orbitals.size());
	for (int i = 0; i < orbitals.size(); ++i)
		multiplicities[i] = orbToMultiplicity(orbitals[i]);
	return LargeRotationMatrix(multiplicities, mat, addSpin);
}

//TODO: Anyone have a better idea where to put this?
template<typename T>
vector<T> dictToList(const map<string, T> &dict, const vector<string> &keys)
{
	vector<T> list;
	typename map<string, T>::const_iterator iter;
	for (int i = 0; i < keys.size(); ++i)
	{
		iter = dict.find(keys[i]);
		if (iter != dict.end())
			list.push_back(iter->second);
		else
			cout << "Warning: No entry \"" << keys[i] << "\" was found in dictionary!" << endl;
	}
	return list;
}


map<string, vector<uint16_t> > createAtmoicIndicesDict(const vector<string> &orbitals)
{
	map<string, vector<uint16_t> > allIndices;
	uint16_t offset = 0;
	for (uint i = 0; i < (uint) orbitals.size(); ++i)
	{
		string orb = orbitals[i];
		int multiplicity = orbToMultiplicity(orb);
		if (multiplicity > 0)
		{
			vector<uint16_t> orbIndices(multiplicity);
			for (uint16_t j = 0; j < (uint16_t) multiplicity; ++j)
				orbIndices[j] = offset + j;
			allIndices.insert( { orb, orbIndices });
			if (!isRelOrbital(orb) && (multiplicity % 2 == 0))
			{
				vector<uint16_t> orbIndicesDn(multiplicity / 2);
				vector<uint16_t> orbIndicesUp(multiplicity / 2);
				for (uint16_t j = 0; j < (uint16_t) multiplicity / 2; ++j)
				{
					orbIndicesDn[j] = offset + 2 * j;
					orbIndicesUp[j] = offset + 2 * j + 1;
				}
				allIndices.insert( { orb + "_Dn", orbIndicesDn });
				allIndices.insert( { orb + "_Up", orbIndicesUp });
			}
			offset += multiplicity;
		} else
		{
			cout << "\nWarning: Key \"" << orb << "\" could not be interpreted as an orbital.\n";
		}
	}
	return allIndices;
}

void groupOrbitals(map<string, vector<uint16_t> > &allIndices, const map<string, vector<string> > &groupings)
{
	map<string, vector<string> >::const_iterator iter;
	for (iter = groupings.begin(); iter != groupings.end(); ++iter)
	{
		string groupKey = iter->first;
		vector < string > groupMembers = iter->second;
		vector<uint16_t> indices;
		for (uint i = 0; i < (uint) groupMembers.size(); ++i)
		{
			string groupMember = groupMembers[i];
			map<string, vector<uint16_t> >::const_iterator groupMemberIter = allIndices.find(groupMember);
			if (groupMemberIter != allIndices.end())
			{
				vector<uint16_t> groupMemberIndices = groupMemberIter->second;
				for (uint j = 0; j < (uint) groupMemberIndices.size(); ++j)
				{
					if (!isElementOf(groupMemberIndices[j], indices))
						indices.push_back(groupMemberIndices[j]);
				}
			} else
			{
				cout << "Warning: Orbital \"" << groupMember << "\" could not be found." << endl;
			}
		}
		allIndices.insert( { groupKey, indices });
	}
}

map<string, vector<uint16_t> > createAtmoicIndicesDict(const vector<string> &orbitals, const map<string, vector<string> > &groupings)
{
	map<string, vector<uint16_t> > allIndices = createAtmoicIndicesDict(orbitals);
	groupOrbitals(allIndices, groupings);
	return allIndices;
}
