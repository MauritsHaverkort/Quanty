char LToShellName(int l)
{
  switch(l)
  {
    case 0:
      return 's';
    case 1:
      return 'p';
    case 2:
      return 'd';
    case 3:
      return 'f';
    default:
      return 'x';
  }
}

int LMToOrbitalName(int l, int m, char **name)
{
  *name = (char*)malloc(16*sizeof(char));
  switch (l*l+l+m)
  {
    case 0:
      strcpy(*name, "s");
      break;
    case 1:
      strcpy(*name, "p_y");
      break;
    case 2:
      strcpy(*name, "p_z");
      break;
    case 3:
      strcpy(*name, "p_x");
      break;
    case 4:
      strcpy(*name, "d_{xy}");
      break;
    case 5:
      strcpy(*name, "d_{yz}");
      break;
    case 6:
      strcpy(*name, "d_{3z^2-r^2}");
      break;
    case 7:
      strcpy(*name, "d_{xz}");
      break;
    case 8:
      strcpy(*name, "d_{x^2-y^2}");
      break;
    case 9:
      strcpy(*name, "f_{y^3-3x^2y}");
      break;
    case 10:
      strcpy(*name, "f_{xyz}");
      break;
    case 11:
      strcpy(*name, "f_{5z^2y-yr^2}");
      break;
    case 12:
      strcpy(*name, "f_{5z^3-3zr^2}");
      break;
    case 13:
      strcpy(*name, "f_{5z^2x-xr^2}");
      break;
    case 14:
      strcpy(*name, "f_{x^2z-y^2z}");
      break;
    case 15:
      strcpy(*name, "f_{x^3-3xy^2}");
      break;
    default:
      strcpy(*name, "??");
      break;
  }

  return 0;
}

std::string LPlusJToShellName(int lPlusJ)
{
	switch (lPlusJ)
	{
	case 0:
		return "s1/2";
	case 1:
		return "p1/2";
	case 2:
		return "p3/2";
	case 3:
		return "d3/2";
	case 4:
		return "d5/2";
	case 5:
		return "f5/2";
	case 6:
		return "f7/2";
	case 7:
		return "g7/2";
	case 8:
		return "g9/2";
	default:
		return "??";
	}
}

int LPlusJAndJzToOrbitalName(int lPlusJ, int jz, char **name)
{
	*name = (char*) malloc(16 * sizeof(char));
	std::stringstream ss;
	ss << LPlusJToShellName(lPlusJ);
	if (jz > 0) ss << "+";
	ss << (2 * jz - 1);
	ss << "/2";
	strcpy(*name, ss.str().c_str());
	return 0;
}

int EqualFloat(double a, double b, double epsilon)
{
  if(fabs(a-b)<epsilon) // pricision defined by the fplo output (coordinates).
    return 1;
  else
    return 0;
}

uint32_t CrystalStructureFromDresdenFPLO(DresdenFPLO FPLO, CrystalStructureType *Cell)
{
  int i, j, k, allatom, iatom, inl, is;
  int nlIndex[FPLO.Wannier.NActAtom][FPLO.Wannier.Nwandef], siteIndex[FPLO.Wannier.NActAtom][FPLO.Wannier.Nwandef];
  int ActMap[FPLO.Atoms.Natom];
  memset(nlIndex, -1, FPLO.Wannier.NActAtom*FPLO.Wannier.Nwandef*sizeof(int));
  memset(siteIndex, -1, FPLO.Wannier.NActAtom*FPLO.Wannier.Nwandef*sizeof(int));
  memset(ActMap, -1, FPLO.Atoms.Natom*sizeof(int));

  Cell->Units = 3;
  for(i=0;i<3;i++)
  {
    for(j=0;j<3;j++)
    {
      Cell->rUnitCell[i][j] = FPLO.Lattice.LattVec[i][j];
    }
  }
  ReciprocalLatticeVector(Cell->rUnitCell, Cell->gUnitCell, Cell->Units);

  Cell->NAtoms = FPLO.Wannier.NActAtom;
  Cell->Atom   = (AtomType*)malloc(Cell->NAtoms*sizeof(AtomType));
  for(allatom=0;allatom<FPLO.Atoms.Natom;allatom++)
  {
    iatom = FPLO.Wannier.ActMap[allatom];
    if(iatom!=-1)
    {
      strcpy(Cell->Atom[iatom].Name, FPLO.Atoms.Atom[allatom]);
      Cell->Atom[iatom].AtomNumber = AtomNameToNumber(Cell->Atom[iatom].Name);
      for(i=0;i<3;i++)
      {
        Cell->Atom[iatom].Position[i] = FPLO.Atoms.Pos[allatom][i];
      }
      // Determine number of shells on each site.
      Cell->Atom[iatom].Nshell = 0;
      inl = 0;
      is  = 0;

      for(i=0;i<FPLO.Wannier.Nwandef;i++)
      {
        if(FPLO.Wannier.ActMap[FPLO.Wannier.Site[i]-1]==iatom)
        {
          siteIndex[iatom][is++] = i;
          if(Cell->Atom[iatom].Nshell==0)
          {
            Cell->Atom[iatom].Nshell = 1;
            nlIndex[iatom][inl++] = 10*FPLO.Wannier.Orbital[i][0]+FPLO.Wannier.Orbital[i][1];
          }
          else
          {
            for(j=0;j<inl;j++)
            {
              if( (10*FPLO.Wannier.Orbital[i][0]+FPLO.Wannier.Orbital[i][1])==nlIndex[iatom][j] ) break;
            }
            if(j==inl)
            {
              Cell->Atom[iatom].Nshell += 1;
              nlIndex[iatom][inl++] = 10*FPLO.Wannier.Orbital[i][0]+FPLO.Wannier.Orbital[i][1];
            }
          }
        }
      }

      // Determine properties of each shell.
      if(Cell->Atom[iatom].Nshell!=0)
      {
        Cell->Atom[iatom].ShellName      = (char**)malloc(Cell->Atom[iatom].Nshell*sizeof(char*));
        Cell->Atom[iatom].shelldegeneracy= (int  *)malloc(Cell->Atom[iatom].Nshell*sizeof(int  ));
        Cell->Atom[iatom].OrbitalName    = (char***)malloc(Cell->Atom[iatom].Nshell*sizeof(char**));

        // shell degeneracy
				for (i = 0; i < Cell->Atom[iatom].Nshell; i++)
				{
					Cell->Atom[iatom].ShellName[i] = (char*) malloc(6 * sizeof(char));
					Cell->Atom[iatom].ShellName[i][0] = (char) (nlIndex[iatom][i] / 10 + '0');
					if ((bool) FPLO.Wannier.QRelativistic)
					{
						std::string tmp = LPlusJToShellName(nlIndex[iatom][i] % 10);
						Cell->Atom[iatom].ShellName[i][1] = tmp.c_str()[0];
						Cell->Atom[iatom].ShellName[i][2] = tmp.c_str()[1];
						Cell->Atom[iatom].ShellName[i][3] = '/';
						Cell->Atom[iatom].ShellName[i][4] = '2';
						Cell->Atom[iatom].ShellName[i][5] = 0;
					}
					else
					{
						Cell->Atom[iatom].ShellName[i][1] = LToShellName(nlIndex[iatom][i] % 10);
						Cell->Atom[iatom].ShellName[i][2] = 0;
					}

          Cell->Atom[iatom].shelldegeneracy[i] = 0;
          for(j=0;j<is;j++)
          {
            if( (10*FPLO.Wannier.Orbital[siteIndex[iatom][j]][0]+FPLO.Wannier.Orbital[siteIndex[iatom][j] ][1])==nlIndex[iatom][i] )
            {
              Cell->Atom[iatom].shelldegeneracy[i] += 1;
            }
          }

          // orbital names in each shell
          Cell->Atom[iatom].OrbitalName[i] = (char**)malloc(Cell->Atom[iatom].shelldegeneracy[i]*sizeof(char*));
          k = 0;
          for(j=0;j<is;j++)
          {
            if( (10*FPLO.Wannier.Orbital[siteIndex[iatom][j]][0]+FPLO.Wannier.Orbital[siteIndex[iatom][j] ][1])==nlIndex[iatom][i] )
            {
            	if ((bool)FPLO.Wannier.QRelativistic)
            		LPlusJAndJzToOrbitalName(FPLO.Wannier.Orbital[siteIndex[iatom][j]][1], FPLO.Wannier.Orbital[siteIndex[iatom][j]][2], &(Cell->Atom[iatom].OrbitalName[i][k++]));
            	else
            		LMToOrbitalName(FPLO.Wannier.Orbital[siteIndex[iatom][j]][1], FPLO.Wannier.Orbital[siteIndex[iatom][j]][2], &(Cell->Atom[iatom].OrbitalName[i][k++]));
            }
          }
        }
      }
    }
  }

  return 0;
}

uint32_t TightBindingDefFromDresdenFPLO(DresdenFPLO FPLO, TightBindingDef *HTB)
{
  int i, j, k, l, iatom, ihtb, ish, jsh, indwan, indsh, diffQ;
  int Nwandef=FPLO.Wannier.Nwandef;
  int Nshell, *Norb, *ShellInd, **ShellOfSite, *OrbInd, **NVecShell, **NHopStart;
  char **SiteName, **ShellName;
  double ***HopVec, **SitePos;
  vector < vector<uint> > HopVecJVals; //for freeing the memory in the end.
  const int NhopMax=2048;
  double AtomPositionAccuracy = 0.001;

  strncpy(HTB->Name, FPLO.Compound, 15);
  HTB->Name[15]=0;

  CrystalStructureFromDresdenFPLO(FPLO, &(HTB->Cell));
  SetNFermionIndexTightBindingDef(HTB);

	// Determine the total number of shells.
	Nshell = 0;
	for (iatom = 0; iatom < FPLO.Wannier.NActAtom; iatom++)
	{
		Nshell += HTB->Cell.Atom[iatom].Nshell;
	}

	int lastSite = FPLO.Wannier.Site[0];
	int lastN = FPLO.Wannier.Orbital[0][0];
	int lastShell = FPLO.Wannier.Orbital[0][1];
	iatom = 0;
	ish = 0;
	int siteSpecificShellCounter = 0;
	ShellInd = (int*) malloc(Nwandef * sizeof(int));
	SiteName = (char**) malloc(FPLO.Wannier.NActAtom * sizeof(char*));
	SitePos = (double**) malloc(FPLO.Wannier.NActAtom * sizeof(double*));
	ShellName = (char**) malloc(Nshell * sizeof(char*));
	ShellOfSite = (int**) malloc(Nshell * sizeof(int*));
	Norb = (int*) calloc(Nshell, sizeof(int));
	OrbInd = (int*) calloc(Nwandef, sizeof(int));
	for (i = 0; i < Nwandef; i++)
	{
		bool siteChanged = lastSite != FPLO.Wannier.Site[i];
		bool NChanged = lastN != FPLO.Wannier.Orbital[i][0];
		bool shellChanged = lastShell != FPLO.Wannier.Orbital[i][1];

		if (i == 0 || siteChanged)
		{
			if (siteChanged) ++iatom;

			SiteName[iatom] = (char*) malloc(MAXStrLength * sizeof(char));
			strcpy(SiteName[iatom], FPLO.Atoms.Atom[FPLO.Wannier.Site[i] - 1]);

			SitePos[iatom] = (double*) calloc(3, sizeof(double));
			for (k = 0; k < 3; k++)
				SitePos[iatom][k] = FPLO.Wannier.PosCenter[i][k];
		}

		if (i == 0 || siteChanged || NChanged || shellChanged)
		{
			if (siteChanged || NChanged || shellChanged) ++ish;
			if (siteChanged)
			{
				siteSpecificShellCounter = 0;
			}
			else if (i != 0)
			{
				++siteSpecificShellCounter;
			}
			ShellName[ish] = (char*) malloc(4 * sizeof(char));
			ShellName[ish][0] = (char) (FPLO.Wannier.Orbital[i][0] + '0');
			if ((bool) FPLO.Wannier.QRelativistic)
			{
				std::string tmp = LPlusJToShellName(FPLO.Wannier.Orbital[i][1]);
				ShellName[ish][1] = tmp.c_str()[0];
				ShellName[ish][2] = tmp.c_str()[1];
				ShellName[ish][3] = '/';
				ShellName[ish][4] = '2';
				ShellName[ish][5] = 0;
			}
			else
			{
				ShellName[ish][1] = LToShellName(FPLO.Wannier.Orbital[i][1]);
				ShellName[ish][2] = 0;
			}

			ShellOfSite[ish] = (int*) calloc(2, sizeof(int));
			ShellOfSite[ish][0] = iatom;
			ShellOfSite[ish][1] = siteSpecificShellCounter;
		}
		else
		{
			OrbInd[i] = OrbInd[i - 1] + 1;
		}

		ShellInd[i] = ish;
		++Norb[ish];

		lastSite = FPLO.Wannier.Site[i];
		lastN = FPLO.Wannier.Orbital[i][0];
		lastShell = FPLO.Wannier.Orbital[i][1];
	}

  NVecShell = (int**)malloc(Nshell*sizeof(int*));
  HopVec    = (double***)malloc(Nshell*Nshell*sizeof(double**));
  for(i=0;i<Nshell;i++)
  {
    NVecShell[i] = (int*)calloc(Nshell, sizeof(int));
    for(j=0;j<Nshell;j++)
    {
      HopVec[i*Nshell+j] = (double**)malloc(NhopMax*sizeof(double*));
      HopVecJVals.push_back(vector<uint>());
    }
  }

  // Determine the unique translational vectors for each shell.
  // The vectors will be used to label the HTB list.
  // Hopvec in units of lattice parameters
  // NVecShell[i][j] sets the number of vectors that connect shell i to shell j
  // Hopvec[i*Nj+j][k] contains the k-th vector connecting shell i to shell j
  for(i=0;i<Nwandef;i++)
  {
    for(j=0;j<Nwandef;j++)
    {
      ish    = ShellInd[i];
      jsh    = ShellInd[j];
      indsh  = ish*Nshell+jsh;
      indwan = i*Nwandef+j;
      for(k=0;k<FPLO.Wannier.Nhop[indwan];k++)
      {
        if(NVecShell[ish][jsh]==0)
        {
          HopVec[indsh][0]    = (double*)malloc(3*sizeof(double));
          HopVecJVals[indsh].push_back(0);
          HopVec[indsh][0][0] = FPLO.Wannier.HopVec[indwan][k][0];
          HopVec[indsh][0][1] = FPLO.Wannier.HopVec[indwan][k][1];
          HopVec[indsh][0][2] = FPLO.Wannier.HopVec[indwan][k][2];

          NVecShell[ish][jsh]++;
        }
        else
        {
          l = 0;
          diffQ = 0;
          while( l<NVecShell[ish][jsh] &&
                (diffQ =
                  !(EqualFloat(HopVec[indsh][l][0],FPLO.Wannier.HopVec[indwan][k][0],AtomPositionAccuracy)&&
                    EqualFloat(HopVec[indsh][l][1],FPLO.Wannier.HopVec[indwan][k][1],AtomPositionAccuracy)&&
                    EqualFloat(HopVec[indsh][l][2],FPLO.Wannier.HopVec[indwan][k][2],AtomPositionAccuracy))) )
              l++;

          if(l==NVecShell[ish][jsh] && diffQ)
          {
            HopVec[indsh][l]    = (double*)malloc(3*sizeof(double));
            HopVecJVals[indsh].push_back(l);
            HopVec[indsh][l][0] = FPLO.Wannier.HopVec[indwan][k][0];
            HopVec[indsh][l][1] = FPLO.Wannier.HopVec[indwan][k][1];
            HopVec[indsh][l][2] = FPLO.Wannier.HopVec[indwan][k][2];

            NVecShell[ish][jsh]++;
          }
        }
      }
    }
  }

  // set the values of NHopStart
  // NHopStart[i][j] contains the index in HTB where the hopping from shell i to shell j starts
  NHopStart = (int**)malloc(Nshell*sizeof(int*));
  for(ish=0;ish<Nshell;ish++)
  {
    NHopStart[ish] = (int*)calloc(Nshell, sizeof(int));
  }

  HTB->nH = 0;
  for(ish=0;ish<Nshell;ish++)
  {
    for(jsh=0;jsh<Nshell;jsh++)
    {
      HTB->nH += NVecShell[ish][jsh];
      if(jsh==Nshell-1)
      { // at the end of hopping from i to j (i.e. j=max) we need to set the starting point of the next i
        if(ish!=Nshell-1)
        { // but only if there is a next i
          NHopStart[ish+1][0] = NHopStart[ish][jsh] + NVecShell[ish][jsh];
        }
      }
      else
      { // set the starting point of the next j (which is available due to the test before)
        NHopStart[ish][jsh+1] = NHopStart[ish][jsh] + NVecShell[ish][jsh];
      }
    }
  }
  HTB->H= (TightBindingDefHopping*)malloc(HTB->nH*sizeof(TightBindingDefHopping));

  // Fill in the HTB list.
  // TransVecs distance between the unit cells of the involved atoms in units of unit cells
  // ConnectingVector Distance between the atoms involved in untis of Cell.Units.
  for(ish=0;ish<Nshell;ish++)
  {
    for(jsh=0;jsh<Nshell;jsh++)
    {
      indsh = ish*Nshell+jsh;
      for(ihtb=NHopStart[ish][jsh];ihtb<NHopStart[ish][jsh]+NVecShell[ish][jsh];ihtb++)
      {
        HTB->H[ihtb].Atompos1  = ShellOfSite[ish][0];
        HTB->H[ihtb].Shellpos1 = ShellOfSite[ish][1];
        HTB->H[ihtb].Atompos2  = ShellOfSite[jsh][0];
        HTB->H[ihtb].Shellpos2 = ShellOfSite[jsh][1];
        if(((HTB->Cell.Units&8)==0)&&(FPLO.Lattice.Molecule==false))
        {
          HTB->H[ihtb].ConnectingVector[0] = ((HopVec[indsh][ihtb-NHopStart[ish][jsh]][0]==0)?0:HopVec[indsh][ihtb-NHopStart[ish][jsh]][0]*HTB->Cell.rUnitCell[0][0]) + ((HopVec[indsh][ihtb-NHopStart[ish][jsh]][1]==0)?0:HopVec[indsh][ihtb-NHopStart[ish][jsh]][1]*HTB->Cell.rUnitCell[1][0]) + ((HopVec[indsh][ihtb-NHopStart[ish][jsh]][2]==0)?0:HopVec[indsh][ihtb-NHopStart[ish][jsh]][2]*HTB->Cell.rUnitCell[2][0]);
          HTB->H[ihtb].ConnectingVector[1] = ((HopVec[indsh][ihtb-NHopStart[ish][jsh]][0]==0)?0:HopVec[indsh][ihtb-NHopStart[ish][jsh]][0]*HTB->Cell.rUnitCell[0][1]) + ((HopVec[indsh][ihtb-NHopStart[ish][jsh]][1]==0)?0:HopVec[indsh][ihtb-NHopStart[ish][jsh]][1]*HTB->Cell.rUnitCell[1][1]) + ((HopVec[indsh][ihtb-NHopStart[ish][jsh]][2]==0)?0:HopVec[indsh][ihtb-NHopStart[ish][jsh]][2]*HTB->Cell.rUnitCell[2][1]);
          HTB->H[ihtb].ConnectingVector[2] = ((HopVec[indsh][ihtb-NHopStart[ish][jsh]][0]==0)?0:HopVec[indsh][ihtb-NHopStart[ish][jsh]][0]*HTB->Cell.rUnitCell[0][2]) + ((HopVec[indsh][ihtb-NHopStart[ish][jsh]][1]==0)?0:HopVec[indsh][ihtb-NHopStart[ish][jsh]][1]*HTB->Cell.rUnitCell[1][2]) + ((HopVec[indsh][ihtb-NHopStart[ish][jsh]][2]==0)?0:HopVec[indsh][ihtb-NHopStart[ish][jsh]][2]*HTB->Cell.rUnitCell[2][2]);
        }
        else
        {
          HTB->H[ihtb].ConnectingVector[0] = HopVec[indsh][ihtb-NHopStart[ish][jsh]][0];
          HTB->H[ihtb].ConnectingVector[1] = HopVec[indsh][ihtb-NHopStart[ish][jsh]][1];
          HTB->H[ihtb].ConnectingVector[2] = HopVec[indsh][ihtb-NHopStart[ish][jsh]][2];
        }
        if(FPLO.Lattice.Molecule==false)
        {
          HTB->H[ihtb].TransVecs[0] = HopVec[indsh][ihtb-NHopStart[ish][jsh]][0] - (SitePos[ ShellOfSite[jsh][0] ][0]-SitePos[ ShellOfSite[ish][0] ][0]);
          HTB->H[ihtb].TransVecs[1] = HopVec[indsh][ihtb-NHopStart[ish][jsh]][1] - (SitePos[ ShellOfSite[jsh][0] ][1]-SitePos[ ShellOfSite[ish][0] ][1]);
          HTB->H[ihtb].TransVecs[2] = HopVec[indsh][ihtb-NHopStart[ish][jsh]][2] - (SitePos[ ShellOfSite[jsh][0] ][2]-SitePos[ ShellOfSite[ish][0] ][2]);
        }
        else
        {
          HTB->H[ihtb].TransVecs[0] = 0;
          HTB->H[ihtb].TransVecs[1] = 0;
          HTB->H[ihtb].TransVecs[2] = 0;
        }
        HTB->H[ihtb].Hop.QComplex = FPLO.Wannier.QComplex;
        HTB->H[ihtb].Hop.Ni = Norb[ish];
        HTB->H[ihtb].Hop.Nj = Norb[jsh];
        InitMatrix(&(HTB->H[ihtb].Hop));
        sprintf(HTB->H[ihtb].Hop.Name, "%4s %4s %4i / %4s %4s %4i",
                SiteName[HTB->H[ihtb].Atompos1],
                ShellName[ish],
                HTB->H[ihtb].Atompos1,
                SiteName[HTB->H[ihtb].Atompos2],
                ShellName[jsh],
                HTB->H[ihtb].Atompos2
              );
				if ((bool) FPLO.Wannier.QComplex)
				{
					for (i = 0; i < Nwandef; i++)
					{
						for (j = 0; j < Nwandef; j++)
						{
							// Fill in the matrix elements.
							indwan = i * Nwandef + j;
							if (ShellInd[i] == ish && ShellInd[j] == jsh)
							{
								for (k = 0; k < FPLO.Wannier.Nhop[indwan]; k++)
								{
									if (EqualFloat(FPLO.Wannier.HopVec[indwan][k][0], HopVec[indsh][ihtb - NHopStart[ish][jsh]][0], AtomPositionAccuracy)
											&& EqualFloat(FPLO.Wannier.HopVec[indwan][k][1], HopVec[indsh][ihtb - NHopStart[ish][jsh]][1], AtomPositionAccuracy)
											&& EqualFloat(FPLO.Wannier.HopVec[indwan][k][2], HopVec[indsh][ihtb - NHopStart[ish][jsh]][2], AtomPositionAccuracy))
									{
										HTB->H[ihtb].Hop.VRe[OrbInd[i]][OrbInd[j]] = FPLO.Wannier.HopVal[indwan][k];
										HTB->H[ihtb].Hop.VIm[OrbInd[i]][OrbInd[j]] = FPLO.Wannier.HopValIm[indwan][k];
									}
								}
							}
						}
					}
				}
				else //Hopping real
				{
					for (i = 0; i < Nwandef; i++)
					{
						for (j = 0; j < Nwandef; j++)
						{
							// Fill in the matrix elements.
							indwan = i * Nwandef + j;
							if (ShellInd[i] == ish && ShellInd[j] == jsh)
							{
								for (k = 0; k < FPLO.Wannier.Nhop[indwan]; k++)
								{
									if (EqualFloat(FPLO.Wannier.HopVec[indwan][k][0], HopVec[indsh][ihtb - NHopStart[ish][jsh]][0], AtomPositionAccuracy)
											&& EqualFloat(FPLO.Wannier.HopVec[indwan][k][1], HopVec[indsh][ihtb - NHopStart[ish][jsh]][1], AtomPositionAccuracy)
											&& EqualFloat(FPLO.Wannier.HopVec[indwan][k][2], HopVec[indsh][ihtb - NHopStart[ish][jsh]][2], AtomPositionAccuracy))
									{
										HTB->H[ihtb].Hop.VRe[OrbInd[i]][OrbInd[j]] = FPLO.Wannier.HopVal[indwan][k];
									}
								}
							}
						}
					}
				}
			}
		}
  }
	for (ish = 0; ish < Nshell; ish++)
	{
		free(ShellOfSite[ish]);
		free(NVecShell[ish]);
		free(NHopStart[ish]);
		free(ShellName[ish]);
	}
	for (uint i = 0; i < Nshell * Nshell; ++i)
	{
		for (uint j = 0; j < HopVecJVals[i].size(); ++j)
		{
			free(HopVec[i][HopVecJVals[i][j]]);
		}
		free(HopVec[i]);
	}
	free(HopVec);
	for (i = 0; i < FPLO.Wannier.NActAtom; i++)
	{
		free(SiteName[i]);
		free(SitePos[i]);
	}
	free(Norb);
	free(ShellInd);
	free(ShellOfSite);
	free(OrbInd);
	free(NVecShell);
	free(NHopStart);
	free(SiteName);
	free(ShellName);
	free(SitePos);
	return 0;
}

uint32_t HybridizationFunctionFromFPLO(DresdenFPLO *FPLO, BlockListOfPoles *LP, int *OrbitalList, bool MergeNegativeToNeighbour, bool SetNegativeToZero, bool Normalize, CompactVectorType EnergyWindow)
{
  int i,j,jMin,jMax;
  int const stepsize = 1 + Sqr( LP->BlockSize );
  ListOfPoles lDOSin;
  LP->N=FPLO->lDOS.NE;
  jMin = 0;
  jMax = FPLO->lDOS.NE;
  if(SetNegativeToZero==true)
  {
    LP->QComplex=0;
  }
  else
  {
    LP->QComplex=1;
  }
  InitBlockListOfPoles(LP);
  lDOSin.N=LP->N;
  InitListOfPoles(&lDOSin);
  // set energy
  for(j=0;j<LP->N;j++)
  {
    lDOSin.a[j] = (((double) j ) /  (LP->N-1) ) * FPLO->lDOS.Emax + (((double) (LP->N-1-j) ) /  (LP->N-1) ) * FPLO->lDOS.Emin;
    LP->VRe[j * stepsize] = lDOSin.a[j];
  }
  // determine jMin and jMax from energy window
  if(EnergyWindow.Ni==2)
  {
    for(j=0;j<FPLO->lDOS.NE;j++)
    {
      lDOSin.b[j]=0;
      if(lDOSin.a[j] <  EnergyWindow.VRe[0]) jMin = j;
      if(lDOSin.a[j] <= EnergyWindow.VRe[1]) jMax = j;
    }
    jMin++;
    jMax++;
  }
  //loop over all lDOSes
  for(i=0;i<LP->BlockSize;i++)
  {
    // read lDOS to lDOSin
    for(j=jMin;j<jMax;j++)
    {
      lDOSin.b[j]=FPLO->lDOS.V[OrbitalList[i]][j];
    }
    // reglarize
    if(MergeNegativeToNeighbour==true) MergeNegativeToNeighbourListOfPolesSorted(&lDOSin);
    if(SetNegativeToZero==true) SetNegativeToZeroListOfPoles(&lDOSin);
    if(Normalize==true)
    {
      NormalizeListOfPoles(&lDOSin);
    }
    else
    {
      ScaleListOfPoles(&lDOSin,0.5 * (FPLO->lDOS.Emax-FPLO->lDOS.Emin) / LP->N );
    }
    // set the correct elements of Block List Of Poles
    if(LP->QComplex==0)
    {
      for(j=0;j<LP->N;j++)
      {
        LP->VRe[j * stepsize + 1 + i * ( 1 + LP->BlockSize)] = lDOSin.b[j] ;
      }
    }
    else
    {
      for(j=0;j<LP->N;j++)
      {
        LP->V[j * stepsize + 1 + i * ( 1 + LP->BlockSize)].Im = 0;
        LP->V[j * stepsize + 1 + i * ( 1 + LP->BlockSize)].Re = lDOSin.b[j] ;
      }
    }
  }
  FreeListOfPoles(&lDOSin);
  return 0;
}

int FreeDresdenFPLO(DresdenFPLO *FPLO)
{
  // Atoms
	for (int i = 0; i < FPLO->Atoms.Natom; ++i)
	{
		free(FPLO->Atoms.Atom[i]);
		free(FPLO->Atoms.Pos[i]);
	}
  if(FPLO->Atoms.Natom!=0)
  {
    free(FPLO->Atoms.Atom);
    free(FPLO->Atoms.Pos);
  }
  // Wannier
	for (int i = 0; i < FPLO->Wannier.Nwandef; ++i)
	{
		free(FPLO->Wannier.Name[i]);
		free(FPLO->Wannier.Orbital[i]);
		free(FPLO->Wannier.PosCenter[i]);
	}
	for (int i = 0; i < FPLO->Wannier.Nwandef * FPLO->Wannier.Nwandef; ++i)
	{
		for (int j = 0; j < FPLO->Wannier.Nhop[i]; ++j)
		{
			free(FPLO->Wannier.HopVec[i][j]);
		}
		free(FPLO->Wannier.HopVec[i]);
		free(FPLO->Wannier.HopVal[i]);
		free(FPLO->Wannier.HopValIm[i]);
	}
  if(FPLO->Wannier.Nwandef !=0)
  {
    free(FPLO->Wannier.Name);
    free(FPLO->Wannier.Site);
    free(FPLO->Wannier.Orbital);
    free(FPLO->Wannier.PosCenter);
    free(FPLO->Wannier.Nhop);
    free(FPLO->Wannier.HopVec);
    free(FPLO->Wannier.HopVal);
    free(FPLO->Wannier.HopValIm);
    free(FPLO->Wannier.ActMap);
  }
  // lDOS
  if(FPLO->lDOS.N!=0)
  {
    free(FPLO->lDOS.sort);
    free(FPLO->lDOS.n);
    free(FPLO->lDOS.l);
    free(FPLO->lDOS.nl);
    free(FPLO->lDOS.m);
    free(FPLO->lDOS.spin);
    for(int i=0;i<FPLO->lDOS.N;i++)
    {
      free(FPLO->lDOS.V[i]);
    }
    free(FPLO->lDOS.V);
  }
	return 0;
}

