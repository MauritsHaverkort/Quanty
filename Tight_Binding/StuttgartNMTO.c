uint32_t BandstructurePathFromStuttgartCTRL(StuttgartCTRL CTRL, BandstructurePathType *Path)
{
  int QPathBrakes=0;
  int i;
  Path->Units=3;
  for(i=0;i<CTRL.SYML.N-1;i++)
  {
    if(!(strcmp(CTRL.SYML.LAB2[i],CTRL.SYML.LAB1[i+1])==0))QPathBrakes++;
  }
  Path->NQ=1+CTRL.SYML.N+QPathBrakes;
  Path->Q=malloc(Path->NQ*sizeof(BandstructurePointType));
  QPathBrakes=0;
  for(i=0;i<CTRL.SYML.N-1;i++)
  {
    Path->Q[i+QPathBrakes].Nq=CTRL.SYML.NQ[i];
    Path->Q[i+QPathBrakes].q[0]=CTRL.SYML.Q1[i][0]*2*Pi/CTRL.STRUC.ALAT;
    Path->Q[i+QPathBrakes].q[1]=CTRL.SYML.Q1[i][1]*2*Pi/CTRL.STRUC.ALAT;
    Path->Q[i+QPathBrakes].q[2]=CTRL.SYML.Q1[i][2]*2*Pi/CTRL.STRUC.ALAT;
    Path->Q[i+QPathBrakes].label=CTRL.SYML.LAB1[i][0];
    if(!(strcmp(CTRL.SYML.LAB2[i],CTRL.SYML.LAB1[i+1])==0))
    {
      QPathBrakes++;
      Path->Q[i+QPathBrakes].Nq=1;
      Path->Q[i+QPathBrakes].q[0]=CTRL.SYML.Q2[i][0]*2*Pi/CTRL.STRUC.ALAT;
      Path->Q[i+QPathBrakes].q[1]=CTRL.SYML.Q2[i][1]*2*Pi/CTRL.STRUC.ALAT;
      Path->Q[i+QPathBrakes].q[2]=CTRL.SYML.Q2[i][2]*2*Pi/CTRL.STRUC.ALAT;
      Path->Q[i+QPathBrakes].label=CTRL.SYML.LAB2[i][0];
    }
  }
  if(CTRL.SYML.N!=0)
  {
    Path->Q[CTRL.SYML.N-1+QPathBrakes].Nq=CTRL.SYML.NQ[CTRL.SYML.N-1];
    Path->Q[CTRL.SYML.N-1+QPathBrakes].q[0]=CTRL.SYML.Q1[CTRL.SYML.N-1][0]*2*Pi/CTRL.STRUC.ALAT;
    Path->Q[CTRL.SYML.N-1+QPathBrakes].q[1]=CTRL.SYML.Q1[CTRL.SYML.N-1][1]*2*Pi/CTRL.STRUC.ALAT;
    Path->Q[CTRL.SYML.N-1+QPathBrakes].q[2]=CTRL.SYML.Q1[CTRL.SYML.N-1][2]*2*Pi/CTRL.STRUC.ALAT;
    Path->Q[CTRL.SYML.N-1+QPathBrakes].label=CTRL.SYML.LAB1[CTRL.SYML.N-1][0];
    Path->Q[CTRL.SYML.N+QPathBrakes].Nq=1;
    Path->Q[CTRL.SYML.N+QPathBrakes].q[0]=CTRL.SYML.Q2[CTRL.SYML.N-1][0]*2*Pi/CTRL.STRUC.ALAT;
    Path->Q[CTRL.SYML.N+QPathBrakes].q[1]=CTRL.SYML.Q2[CTRL.SYML.N-1][1]*2*Pi/CTRL.STRUC.ALAT;
    Path->Q[CTRL.SYML.N+QPathBrakes].q[2]=CTRL.SYML.Q2[CTRL.SYML.N-1][2]*2*Pi/CTRL.STRUC.ALAT;
    Path->Q[CTRL.SYML.N+QPathBrakes].label=CTRL.SYML.LAB2[CTRL.SYML.N-1][0];
  }
  return 0;
}

uint32_t CrystalStructureFromStuttgartCTRL(StuttgartCTRL CTRL, CrystalStructureType *Cell)
{
  int i,j,k,l,m,iatom;
  Cell->Units=3;
  for(i=0;i<3;i++)
  {
    for(j=0;j<3;j++)
    {
      Cell->rUnitCell[i][j]=CTRL.STRUC.PLAT[i][j]*CTRL.STRUC.ALAT;
    }
  }
  ResiprocalLatticeVector(Cell->rUnitCell, Cell->gUnitCell);
  /* determine number of atoms */
  Cell->NAtoms=CTRL.SITE.Nactive;
  Cell->Atom=malloc(Cell->NAtoms*sizeof(AtomType));
  iatom=0;
  for(i=0;i<CTRL.SITE.N;i++) /**** Rewrite using CTRL.SITE.ActiveAtomMap ****/
  {
    if(iatom==Cell->NAtoms) break; /* all atoms found */
    /* Determine number of active shells */
    Cell->Atom[iatom].Nshell=0;
    for(l=0;l<3;l++)
    {
      for(m=-l;m<=l;m++)
      {
        if(CTRL.SITE.IDXDN_M[i][l*l+l+m]==1)
        {
          Cell->Atom[iatom].Nshell++;
          break;
        }
      }
    }
    if(Cell->Atom[iatom].Nshell==0) continue; /* this atom has no active shells */
    /* set name and position of atom */
    strcpy(Cell->Atom[iatom].Name, CTRL.SITE.ATOM[i]);
    for(j=0;j<3;j++)
    {
      Cell->Atom[iatom].Position[j]=CTRL.SITE.POS[i][j]*CTRL.STRUC.ALAT;
    }
    /* allocate space for shells */
    if(Cell->Atom[iatom].Nshell!=0)
    {
      Cell->Atom[iatom].ShellName=malloc(Cell->Atom[iatom].Nshell*sizeof(char*));
      Cell->Atom[iatom].shelldegeneracy=malloc(Cell->Atom[iatom].Nshell*sizeof(int));
      Cell->Atom[iatom].OrbitalName=malloc(Cell->Atom[iatom].Nshell*sizeof(char*));
      for(j=0;j<Cell->Atom[iatom].Nshell;j++)
      {
        Cell->Atom[iatom].ShellName[j]=malloc(4*sizeof(char));
      }
    }
    /* set name of shells, shell multiplicity, allocate space for orbital names and set orbital names */
    j=0; /* index of active shell in atom[i] */
    for(l=0;l<3;l++)
    {
      k=0;
      for(m=-l;m<=l;m++)
      {
        if(CTRL.SITE.IDXDN_M[i][l*l+l+m]==1)k++;
      }
      if(k!=0)
      {
        /* set name of shell */
        switch (CTRL.CLASS.CONF[CTRL.SITE.ICLASS[i]][l])
        {
          case 1:
            Cell->Atom[iatom].ShellName[j][0]='1';
            break;
          case 2:
            Cell->Atom[iatom].ShellName[j][0]='2';
            break;
          case 3:
            Cell->Atom[iatom].ShellName[j][0]='3';
            break;
          case 4:
            Cell->Atom[iatom].ShellName[j][0]='4';
            break;
          case 5:
            Cell->Atom[i].ShellName[j][0]='5';
            break;
          case 6:
            Cell->Atom[iatom].ShellName[j][0]='6';
            break;
          case 7:
            Cell->Atom[iatom].ShellName[j][0]='7';
            break;
          case 8:
            Cell->Atom[iatom].ShellName[j][0]='8';
            break;
          case 9:
            Cell->Atom[iatom].ShellName[j][0]='9';
            break;
          default:
            Cell->Atom[iatom].ShellName[j][0]='?';
            break;
        }
        switch (l)
        {
          case 0:
            Cell->Atom[iatom].ShellName[j][1]='s';
            break;
          case 1:
            Cell->Atom[iatom].ShellName[j][1]='p';
            break;
          case 2:
            Cell->Atom[iatom].ShellName[j][1]='d';
            break;
          case 3:
            Cell->Atom[iatom].ShellName[j][1]='f';
            break;
          default:
            Cell->Atom[iatom].ShellName[j][1]='?';
            break;
        }
        Cell->Atom[iatom].ShellName[j][2]='\0';
        /* set shell multiplicity */
        Cell->Atom[iatom].shelldegeneracy[j]=k;
        /* allocate space for orbital names */
        Cell->Atom[iatom].OrbitalName[j]=malloc(k*sizeof(char*));
        for(m=0;m<k;m++)
        {
          Cell->Atom[iatom].OrbitalName[j][m]=malloc(16*sizeof(char));
        }
        /* set orbital names */
        k=0;
        for(m=-l;m<=l;m++)
        {
          if(CTRL.SITE.IDXDN_M[i][l*l+l+m]==1)
          {
            switch (l*l+l+m)
            {
              case 0:
                strcpy(Cell->Atom[iatom].OrbitalName[j][k], "s");
                break;
              case 1:
                strcpy(Cell->Atom[iatom].OrbitalName[j][k], "p_y");
                break;
              case 2:
                strcpy(Cell->Atom[iatom].OrbitalName[j][k], "p_z");
                break;
              case 3:
                strcpy(Cell->Atom[iatom].OrbitalName[j][k], "p_x");
                break;
              case 4:
                strcpy(Cell->Atom[iatom].OrbitalName[j][k], "d_{xy}");
                break;
              case 5:
                strcpy(Cell->Atom[iatom].OrbitalName[j][k], "d_{yz}");
                break;
              case 6:
                strcpy(Cell->Atom[iatom].OrbitalName[j][k], "d_{3z^2-r^2}");
                break;
              case 7:
                strcpy(Cell->Atom[iatom].OrbitalName[j][k], "d_{xz}");
                break;
              case 8:
                strcpy(Cell->Atom[iatom].OrbitalName[j][k], "d_{x^2-y^2}");
                break;
              case 9:
                strcpy(Cell->Atom[iatom].OrbitalName[j][k], "f_{y^3-3x^2y}");
                break;
              case 10:
                strcpy(Cell->Atom[iatom].OrbitalName[j][k], "f_{xyz}");
                break;
              case 11:
                strcpy(Cell->Atom[iatom].OrbitalName[j][k], "f_{5z^2y-yr^2}");
                break;
              case 12:
                strcpy(Cell->Atom[iatom].OrbitalName[j][k], "f_{5z^3-3zr^2}");
                break;
              case 13:
                strcpy(Cell->Atom[iatom].OrbitalName[j][k], "f_{5z^2x-xr^2}");
                break;
              case 14:
                strcpy(Cell->Atom[iatom].OrbitalName[j][k], "f_{x^2z-y^2z}");
                break;
              case 15:
                strcpy(Cell->Atom[iatom].OrbitalName[j][k], "f_{x^3-3xy^2}");
                break;
              default:
                strcpy(Cell->Atom[iatom].OrbitalName[j][k], "??");
                break;
            }
            k++;
          }
        }
        j++;
      }
    }
    iatom++;
  }
  return 0;
}


uint32_t TightBindingDefFromStuttgartCTRLHAMR(StuttgartCTRL CTRL, StuttgartHAMR HAMR, TightBindingDef *HTB)
{
  
  int iHAMRelement,iHTBelement,iShell1,iShell2,StartShell1,StartShell2,i1,i2;
  strncpy(HTB->Name, CTRL.HEADER, 31);
  HTB->Name[31]='\0';
  CrystalStructureFromStuttgartCTRL(CTRL, &(HTB->Cell));
  SetNFermionIndexTightBindingDef(HTB);
  /* determine number of elements (sites can have multiple shells) */
  HTB->nH=0;
  for(iHAMRelement=0;iHAMRelement<HAMR.Nelement;iHAMRelement++)
  {
    HTB->nH+=HTB->Cell.Atom[ CTRL.SITE.ActiveAtomMap[ HAMR.element[iHAMRelement].IndexAtom1-1]].Nshell * HTB->Cell.Atom[ CTRL.SITE.ActiveAtomMap[ HAMR.element[iHAMRelement].IndexAtom2-1]].Nshell;
  }
  HTB->H=malloc(HTB->nH*sizeof(TightBindingDefHopping));
  /* set ConnectingVectors, TransVecs and elements of HTB */
  iHTBelement=0;
  for(iHAMRelement=0;iHAMRelement<HAMR.Nelement;iHAMRelement++)
  {
    StartShell1=0;
    for(iShell1=0;iShell1<HTB->Cell.Atom[ CTRL.SITE.ActiveAtomMap[ HAMR.element[iHAMRelement].IndexAtom1-1 ]].Nshell;iShell1++)
    {
      StartShell2=0;
      for(iShell2=0;iShell2<HTB->Cell.Atom[ CTRL.SITE.ActiveAtomMap[ HAMR.element[iHAMRelement].IndexAtom2-1 ]].Nshell;iShell2++)
      {
        HTB->H[iHTBelement].Atompos1=CTRL.SITE.ActiveAtomMap[ HAMR.element[iHAMRelement].IndexAtom1-1];
        HTB->H[iHTBelement].Shellpos1=iShell1;
        HTB->H[iHTBelement].Atompos2=CTRL.SITE.ActiveAtomMap[ HAMR.element[iHAMRelement].IndexAtom2-1];
        HTB->H[iHTBelement].Shellpos2=iShell2;
        HTB->H[iHTBelement].TransVecs[0]=HAMR.element[iHAMRelement].TransVecs[0];
        HTB->H[iHTBelement].TransVecs[1]=HAMR.element[iHAMRelement].TransVecs[1];
        HTB->H[iHTBelement].TransVecs[2]=HAMR.element[iHAMRelement].TransVecs[2];
        HTB->H[iHTBelement].ConnectingVector[0]=HAMR.element[iHAMRelement].ConnectingVector[0];
        HTB->H[iHTBelement].ConnectingVector[1]=HAMR.element[iHAMRelement].ConnectingVector[1];
        HTB->H[iHTBelement].ConnectingVector[2]=HAMR.element[iHAMRelement].ConnectingVector[2];
        HTB->H[iHTBelement].Hop.QComplex=0;
        HTB->H[iHTBelement].Hop.Ni=HTB->Cell.Atom[HTB->H[iHTBelement].Atompos1].shelldegeneracy[iShell1];
        HTB->H[iHTBelement].Hop.Nj=HTB->Cell.Atom[HTB->H[iHTBelement].Atompos2].shelldegeneracy[iShell2];
        InitMatrix(&(HTB->H[iHTBelement].Hop));
        sprintf(HTB->H[iHTBelement].Hop.Name, "%4s %4s %4i / %4s %4s %4i",
                HTB->Cell.Atom[HTB->H[iHTBelement].Atompos1].Name,
                HTB->Cell.Atom[HTB->H[iHTBelement].Atompos1].ShellName[iShell1],
                HTB->H[iHTBelement].Atompos1,
                HTB->Cell.Atom[HTB->H[iHTBelement].Atompos2].Name,
                HTB->Cell.Atom[HTB->H[iHTBelement].Atompos2].ShellName[iShell2],
                HTB->H[iHTBelement].Atompos2
                );
        for(i1=0;i1< HTB->H[iHTBelement].Hop.Ni;i1++)
        {
          for(i2=0;i2< HTB->H[iHTBelement].Hop.Nj;i2++)
          {
            HTB->H[iHTBelement].Hop.VRe[i1][i2]=HAMR.element[iHAMRelement].H[StartShell1+i1][StartShell2+i2];
          }
        }
        StartShell2+=HTB->Cell.Atom[HTB->H[iHTBelement].Atompos2].shelldegeneracy[iShell2];
        iHTBelement++;
      }
      StartShell1+=HTB->Cell.Atom[HTB->H[iHTBelement-1].Atompos1].shelldegeneracy[iShell1];
    }
  }
  return 0;
}

