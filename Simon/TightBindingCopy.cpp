 

void CrystalStructureCopy(const CrystalStructureType &CSOld, CrystalStructureType &CSNew)
{
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			CSNew.rUnitCell[i][j] = CSOld.rUnitCell[i][j];
			CSNew.gUnitCell[i][j] = CSOld.gUnitCell[i][j];
		}
	}
	CSNew.NAtoms = CSOld.NAtoms;
	CSNew.Atom = (AtomType*) malloc(CSNew.NAtoms * sizeof(AtomType));
	for(int i=0;i<CSNew.NAtoms;i++)
		CopyAtomType(&(CSOld.Atom[i]),&(CSNew.Atom[i]));
	CSNew.Units = CSOld.Units;
}

void TightBindingDefHoppingCopy(const TightBindingDefHopping &HOld, TightBindingDefHopping &HNew)
{
	HNew.Atompos1 = HOld.Atompos1;
	HNew.Shellpos1 = HOld.Shellpos1;
	HNew.Atompos2 = HOld.Atompos2;
	HNew.Shellpos2 = HOld.Shellpos2;
	for (int i = 0; i < 3; ++i)
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
	for (int i = 0; i < TBNew.nH; ++i)
		TightBindingDefHoppingCopy(TBOld.H[i], TBNew.H[i]);
	TBNew.NFermion = 0;
	SetNFermionIndexTightBindingDef(&TBNew);
}

void TightBindingReadyForRotate(const TightBindingDef &TBOld, TightBindingDef &TBNew, const vector<AtomType> &newAtoms)
{
	strcpy(TBNew.Name, TBOld.Name);
	CrystalStructureCopy(TBOld.Cell, TBNew.Cell);
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			TBNew.Cell.rUnitCell[i][j] = TBOld.Cell.rUnitCell[i][j];
			TBNew.Cell.gUnitCell[i][j] = TBOld.Cell.gUnitCell[i][j];
		}
	}
	TBNew.Cell.Units = TBOld.Cell.Units;
	TBNew.Cell.NAtoms = (uint32_t) newAtoms.size();
	TBNew.Cell.Atom = (AtomType*) malloc(TBNew.Cell.NAtoms * sizeof(AtomType));
	for (int i = 0; i < (int) newAtoms.size(); ++i)
		TBNew.Cell.Atom[i] = newAtoms[i];
	TBNew.nH = 0;
	TBNew.H = NULL;
	TBNew.NFermion = 0;
	SetNFermionIndexTightBindingDef(&TBNew);
}
