void TestCellFunctions()
{
  uint32_t i, j;
  AtomType Atom;
  CrystalStructureType CuLayer;
  CrystalStructureSuperCellType SC;
  strcpy(Atom.Name, "Cu");
  Atom.AtomNumber = 29;
  for (i = 0; i < 3; i++)
    Atom.Position[i] = 0;
  CuLayer.NAtoms = 1;
  CuLayer.Atom = &Atom;
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      CuLayer.rUnitCell[i][j] = 0.0;
  for (i = 0; i < 2; i++)
    CuLayer.rUnitCell[i][i] = 4.0;
  CuLayer.rUnitCell[2][2] = 8.0;
  ResiprocalLatticeVector(CuLayer.rUnitCell, CuLayer.gUnitCell,0);
  SC.Crystal = CuLayer;
  SC.rSuperCell[0][0] = 2;
  SC.rSuperCell[0][1] = 0;
  SC.rSuperCell[0][2] = 0;
  SC.rSuperCell[1][0] = 0;
  SC.rSuperCell[1][1] = 4;
  SC.rSuperCell[1][2] = 0;
  SC.rSuperCell[2][0] = 0;
  SC.rSuperCell[2][1] = 0;
  SC.rSuperCell[2][2] = 1;
  AddAditionalInformationToSuperCell(&SC);
  PrintCrystalStructure(CuLayer);
  printf("=============\n");
  PrintCrystalStructureSuperCell(SC);
}
