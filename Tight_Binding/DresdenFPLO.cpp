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

int EqualFloat(double a, double b)
{
  if(fabs(a-b)<0.0001) // pricision defined by the fplo output (coordinates).
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
        for(i=0;i<Cell->Atom[iatom].Nshell;i++)
        {
          Cell->Atom[iatom].ShellName[i] = (char*)malloc(4*sizeof(char));
          Cell->Atom[iatom].ShellName[i][0] = (char)(nlIndex[iatom][i]/10+'0');
          Cell->Atom[iatom].ShellName[i][1] = LToShellName(nlIndex[iatom][i]%10);
          Cell->Atom[iatom].ShellName[i][2] = 0;

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
  int Nshell, *Norb, *NshellOfSite, *ShellLabel, *SiteInd, *ShellInd, **ShellOfSite, *OrbInd, **NVecShell, **NHopStart;
  char **SiteName, **ShellName;
  double ***HopVec, **SitePos;
  const int NhopMax=1024;

  strncpy(HTB->Name, FPLO.Compound, 15);
  HTB->Name[15]=0;

  CrystalStructureFromDresdenFPLO(FPLO, &(HTB->Cell));
  SetNFermionIndexTightBindingDef(HTB);

  // Determine the total number of shells.
  Nshell = 0;
  for(iatom=0;iatom<FPLO.Wannier.NActAtom;iatom++)
  {
    Nshell += HTB->Cell.Atom[iatom].Nshell;
  }

  ShellLabel  = (int*)calloc(Nwandef, sizeof(int));
  for(i=0;i<Nwandef;i++)
  {
    // Compress the shell indicies {site, n, l, l+m} to one number
    ShellLabel[i] = 1000*FPLO.Wannier.Site[i] + 100*FPLO.Wannier.Orbital[i][0] + 11*FPLO.Wannier.Orbital[i][1] + FPLO.Wannier.Orbital[i][2];
  }

  // Index the sites by number [0, 1, 2, ...]
  SiteInd = (int*)calloc(Nwandef, sizeof(int));
  SiteName= (char**)malloc(FPLO.Wannier.NActAtom*sizeof(char*));
  SiteName[0] = (char*)malloc(16*sizeof(char));
  strcpy(SiteName[0], FPLO.Atoms.Atom[FPLO.Wannier.Site[0]-1]); //FPLO prints site starting from 1.
  SitePos = (double**)malloc(FPLO.Wannier.NActAtom*sizeof(double*));
  for(i=0;i<FPLO.Wannier.NActAtom;i++) SitePos[i] = (double*)calloc(3, sizeof(double));
  for(k=0;k<3;k++) SitePos[0][k] = FPLO.Wannier.PosCenter[0][k];

  l = 1;
  for(i=1;i<Nwandef;i++)
  {
    for(j=0;j<i;j++)
    {
      if(FPLO.Wannier.Site[i]==FPLO.Wannier.Site[j])
      {
        SiteInd[i] = SiteInd[j];
        break;
      }
      else if(j+1==i)
      {
        SiteName[l] = (char*)malloc(MAXStrLength*sizeof(char));
        strcpy(SiteName[l], FPLO.Atoms.Atom[FPLO.Wannier.Site[i]-1]);
        for(k=0;k<3;k++) SitePos[l][k] = FPLO.Wannier.PosCenter[i][k];
        SiteInd[i] = l++;
      }
    }
  }

  // Index the l values of Wannier orbitals by numbers [0, 1, ..., Nshell-1]
  ShellInd    = (int*)calloc(Nwandef, sizeof(int));
  ShellName   = (char**)malloc(Nshell*sizeof(char*));
  ShellName[0]= (char*)malloc(4*sizeof(char));
  ShellName[0][0]=(char)(FPLO.Wannier.Orbital[0][0]+'0');
  ShellName[0][1]=LToShellName(FPLO.Wannier.Orbital[0][1]);
  ShellName[0][2]=0;

  l = 1;
  for(i=1;i<Nwandef;i++)
  {
    for(j=0;j<i;j++)
    {
      if(ShellLabel[i]/10==ShellLabel[j]/10)
      {
        ShellInd[i] = ShellInd[j];
        break;
      }
      else if(j+1==i)
      {
        ShellName[l]= (char*)malloc(4*sizeof(char));
        ShellName[l][0]=(char)(FPLO.Wannier.Orbital[i][0]+'0');
        ShellName[l][1]=LToShellName(FPLO.Wannier.Orbital[i][1]);
        ShellName[l][2]=0;
        ShellInd[i] = l++;
      }
    }
  }

  // Shell degeneracy of each site.
  // Map shell ish to the j-th shell of i-th atom.
  NshellOfSite = (int*)calloc(Nshell, sizeof(int));
  ShellOfSite  = (int**)malloc(Nshell*sizeof(int*));
  for(ish=0;ish<Nshell;ish++)
  {
    ShellOfSite[ish] = (int*)calloc(2, sizeof(int));
    ShellOfSite[ish][1] = -1;
  }
  ShellOfSite[ ShellInd[0] ][0] = SiteInd[0]; // i-th atom
  ShellOfSite[ ShellInd[0] ][1] = 0;          // j-th shell of the atom
  NshellOfSite[SiteInd[0]  ] = 1;
  for(i=1;i<Nwandef;i++)
  {
    for(j=0;j<i;j++)
    {
      if(ShellInd[i]==ShellInd[j])
        break;
      else if(j+1==i)
      {
        ShellOfSite[ ShellInd[i] ][0] = SiteInd[i];
        ShellOfSite[ ShellInd[i] ][1] = NshellOfSite[ SiteInd[i] ]++;
      }
    }
  }

  // Orbital degeneracy of each shell.
  Norb  = (int*)calloc(Nshell, sizeof(int));
  Norb[0]  = 1;
  for(i=1;i<Nwandef;i++)
  {
    for(j=0;j<i;j++)
    {
      if(ShellLabel[i]==ShellLabel[j])
        break;
      else if(j+1==i)
        Norb[ ShellInd[i] ]++;
    }
  }

  // Orbital index within each shell.
  OrbInd = (int*)calloc(Nwandef, sizeof(int));
  for(i=1;i<Nwandef;i++)
  {
    for(j=0;j<i;j++)
    {
      if(ShellLabel[i]/10==ShellLabel[j]/10)
        OrbInd[i] = Max(OrbInd[i], OrbInd[j]+1);
    }
  }

  NVecShell = (int**)malloc(Nshell*sizeof(int*));
  HopVec    = (double***)malloc(Nshell*Nshell*sizeof(double**));
  for(i=0;i<Nshell;i++)
  {
    NVecShell[i] = (int*)calloc(Nshell, sizeof(int));
    for(j=0;j<Nshell;j++)
    {
      HopVec[i*Nshell+j] = (double**)malloc(NhopMax*sizeof(double*));
    }
  }

  // Determine the unique translational vectors for each shell.
  // The vectors will be used to label the HTB list.
  // Hopvec in units of lattice parameters
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
                  !(EqualFloat(HopVec[indsh][l][0],FPLO.Wannier.HopVec[indwan][k][0])&&
                    EqualFloat(HopVec[indsh][l][1],FPLO.Wannier.HopVec[indwan][k][1])&&
                    EqualFloat(HopVec[indsh][l][2],FPLO.Wannier.HopVec[indwan][k][2]))) )
              l++;

          if(l==NVecShell[ish][jsh] && diffQ)
          {
            HopVec[indsh][l]    = (double*)malloc(3*sizeof(double));
            HopVec[indsh][l][0] = FPLO.Wannier.HopVec[indwan][k][0];
            HopVec[indsh][l][1] = FPLO.Wannier.HopVec[indwan][k][1];
            HopVec[indsh][l][2] = FPLO.Wannier.HopVec[indwan][k][2];

            NVecShell[ish][jsh]++;
          }
        }
      }
    }
  }

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
      if(ish<Nshell-1 && jsh==Nshell-1)
        NHopStart[ish+1][0] = NHopStart[ish][jsh] + NVecShell[ish][jsh];
      else
        NHopStart[ish][jsh+1] = NHopStart[ish][jsh] + NVecShell[ish][jsh];
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
        if((HTB->Cell.Units&8)==0)
        {
          HTB->H[ihtb].ConnectingVector[0] = HopVec[indsh][ihtb-NHopStart[ish][jsh]][0]*HTB->Cell.rUnitCell[0][0] + HopVec[indsh][ihtb-NHopStart[ish][jsh]][1]*HTB->Cell.rUnitCell[1][0] + HopVec[indsh][ihtb-NHopStart[ish][jsh]][2]*HTB->Cell.rUnitCell[2][0];
          HTB->H[ihtb].ConnectingVector[1] = HopVec[indsh][ihtb-NHopStart[ish][jsh]][0]*HTB->Cell.rUnitCell[0][1] + HopVec[indsh][ihtb-NHopStart[ish][jsh]][1]*HTB->Cell.rUnitCell[1][1] + HopVec[indsh][ihtb-NHopStart[ish][jsh]][2]*HTB->Cell.rUnitCell[2][1];
          HTB->H[ihtb].ConnectingVector[2] = HopVec[indsh][ihtb-NHopStart[ish][jsh]][0]*HTB->Cell.rUnitCell[0][2] + HopVec[indsh][ihtb-NHopStart[ish][jsh]][1]*HTB->Cell.rUnitCell[1][2] + HopVec[indsh][ihtb-NHopStart[ish][jsh]][2]*HTB->Cell.rUnitCell[2][2];
        }
        else
        {
          HTB->H[ihtb].ConnectingVector[0] = HopVec[indsh][ihtb-NHopStart[ish][jsh]][0];
          HTB->H[ihtb].ConnectingVector[1] = HopVec[indsh][ihtb-NHopStart[ish][jsh]][1];
          HTB->H[ihtb].ConnectingVector[2] = HopVec[indsh][ihtb-NHopStart[ish][jsh]][2];
        }
        HTB->H[ihtb].TransVecs[0] = HopVec[indsh][ihtb-NHopStart[ish][jsh]][0] - (SitePos[ ShellOfSite[jsh][0] ][0]-SitePos[ ShellOfSite[ish][0] ][0]);
        HTB->H[ihtb].TransVecs[1] = HopVec[indsh][ihtb-NHopStart[ish][jsh]][1] - (SitePos[ ShellOfSite[jsh][0] ][1]-SitePos[ ShellOfSite[ish][0] ][1]);
        HTB->H[ihtb].TransVecs[2] = HopVec[indsh][ihtb-NHopStart[ish][jsh]][2] - (SitePos[ ShellOfSite[jsh][0] ][2]-SitePos[ ShellOfSite[ish][0] ][2]);

        HTB->H[ihtb].Hop.QComplex = 0;
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
        for(i=0;i<Nwandef;i++)
        {
          for(j=0;j<Nwandef;j++)
          {
            // Fill in the matrix elements.
            indwan = i*Nwandef+j;
            if(ShellInd[i]==ish&&ShellInd[j]==jsh)
            {
              for(k=0;k<FPLO.Wannier.Nhop[indwan];k++)
              {
                if( EqualFloat(FPLO.Wannier.HopVec[indwan][k][0], HopVec[indsh][ihtb-NHopStart[ish][jsh]][0]) &&
                    EqualFloat(FPLO.Wannier.HopVec[indwan][k][1], HopVec[indsh][ihtb-NHopStart[ish][jsh]][1]) &&
                    EqualFloat(FPLO.Wannier.HopVec[indwan][k][2], HopVec[indsh][ihtb-NHopStart[ish][jsh]][2])   )
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
  return 0;
}
