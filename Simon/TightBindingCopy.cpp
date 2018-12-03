 

void CrystalStructureCopy(const CrystalStructureType &CSOld, CrystalStructureType &CSNew)
{
	for (uint i = 0; i < 3; ++i)
	{
		for (uint j = 0; j < 3; ++j)
		{
			CSNew.rUnitCell[i][j] = CSOld.rUnitCell[i][j];
			CSNew.gUnitCell[i][j] = CSOld.gUnitCell[i][j];
		}
	}
	CSNew.NAtoms = CSOld.NAtoms;
	CSNew.Atom = (AtomType*) malloc(CSNew.NAtoms * sizeof(AtomType));
	for(uint i=0;i<CSNew.NAtoms;i++)
		CopyAtomType(&(CSOld.Atom[i]),&(CSNew.Atom[i]));
	CSNew.Units = CSOld.Units;
}

void TightBindingDefHoppingCopy(const TightBindingDefHopping &HOld, TightBindingDefHopping &HNew)
{
	HNew.Atompos1 = HOld.Atompos1;
	HNew.Shellpos1 = HOld.Shellpos1;
	HNew.Atompos2 = HOld.Atompos2;
	HNew.Shellpos2 = HOld.Shellpos2;
	for (uint i = 0; i < 3; ++i)
	{
			HNew.TransVecs[i] = HOld.TransVecs[i];
			HNew.ConnectingVector[i] = HOld.ConnectingVector[i];
	}
	MatrixCopy(&HNew.Hop, HOld.Hop);
}

void TightBindingCopy(const TightBindingDef &TBOld, TightBindingDef &TBNew)
{
	strcpy(TBNew.Name, TBOld.Name);
	CrystalStructureCopy(TBOld.Cell, TBNew.Cell);
	TBNew.nH = TBOld.nH;
	TBNew.H = (TightBindingDefHopping*) malloc(TBNew.nH * sizeof(TightBindingDefHopping));
	for (uint i = 0; i < TBNew.nH; ++i)
		TightBindingDefHoppingCopy(TBOld.H[i], TBNew.H[i]);
	TBNew.NFermion = 0;
	SetNFermionIndexTightBindingDef(&TBNew);
}
