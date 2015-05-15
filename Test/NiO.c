void TestNiO()
{
  TimeType Time;
  StuttgartCTRL CTRLin;
  StuttgartHAMR HAMRin;
  BandstructurePathType Path;
  CrystalStructureType Cell;
  TightBindingDef HTB, HTBAF;
  MatrixType BANDS, BWValue;
  double Elines[20], EFermi=6.57414, NidShift, UMeanField;
  int Nk[3];
  double kBT=0.0000000001;
  uint32_t NBlockVal,NBlockCon;
  SpectraType DOS, G0Spec[64], G0AFSpec[256], GAFSpec[1024], SpecNat[100],SpectraPES[10],SpectraIPES[10];
  BlockListOfPoles G0, G0AF, GAF, Sigma;
  CrystalStructureSuperCellType AFSuperCell;
  BlockNaturalImpurityOrbitalMatrixType GAFNat;
  CompactMatrixType AT;
  OperatorType HAF;
  uint16_t *CreateString;
  uint32_t NBitsKey=24, NParticle;
  double epsilon,zero,energy;
  uint32_t i;
  WaveFunctionType psi, Opsi;
  DeterminantRestrictionType Restrictions;
  OperatorType TPES[10], TIPES[10];
  TriDiagonalMatrixType SpecTri;

  /*  This routine reads the NMTO DFT output for NiO and creates a cluster Hamiltonian
   *  It allows one to add local correlations and calculate different forms of spectroscopy (PES, IPES, cPES, XAS, RIXS, ...)
   *  It allows one to add a self energy to the neighbours.
   *
   *  On an LDA+U level the self energy is purely real and just a shift of the onsite energy of the orbitals. 
   *    Complex and energy dependent self energies are in principle possible, but the creation of the cluster Hamiltonian (Everywhere a self energy except for the impurity site) 
   *    still needs to be implemented. For real energy independent self energies this is trivial as one can just subtract the change of onsit energy from the Hamiltonian created
   *    using the full Green's function
   */
  
  // set which calculations should be done and which not. Note that not all calculations are needed (DOS and G0 are nice to test, but otherwise just take time
  uint32_t DoBANDS=1, DoDOS=1, DoG0=1, DoG0AF=1, DoGAF=1, DoHAF=1, DoGrd=1, DoPES=1, DoIPES=1;
  
  // read the output from the Stuttgart NMTO code
  if(FileReadStuttgartCTRL(&CTRLin, "NiO/OMTA/CTRL"))printf("ERROR in CTRL file\n");
  if(FileReadStuttgartHAMR(&HAMRin, "NiO/OMTA/HAMR"))printf("ERROR in HAMR file\n");
  
  // create the path in k-space from the CTRL file into the internal format
  BandstructurePathFromStuttgartCTRL(CTRLin, &Path);
  PrintBandstructurePath(Path);
  // create the crystal structure in internal format from the CTRL file
  CrystalStructureFromStuttgartCTRL(CTRLin, &Cell);
  PrintCrystalStructure(Cell);
  // create the TB Hamiltonian from the CTRL and HAMR file in internal format
  TightBindingDefFromStuttgartCTRLHAMR(CTRLin, HAMRin, &HTB);
  FilePrintTightBindingDef(HTB, "NiO/HTB", "w");
  
  // Create the AF supercell
  //   One needs to have a Ni with spin up and a Ni with spin down in the super cell.
  AFSuperCell.rSuperCell[0][0]= 1;
  AFSuperCell.rSuperCell[0][1]= 0;
  AFSuperCell.rSuperCell[0][2]= 1;
  AFSuperCell.rSuperCell[1][0]= 1;
  AFSuperCell.rSuperCell[1][1]=-1;
  AFSuperCell.rSuperCell[1][2]= 0;
  AFSuperCell.rSuperCell[2][0]= 0;
  AFSuperCell.rSuperCell[2][1]=-1;
  AFSuperCell.rSuperCell[2][2]= 1;
  AFSuperCell.Crystal=HTB.Cell;
  AddAditionalInformationToSuperCell(&AFSuperCell);
  PrintCrystalStructureSuperCell(AFSuperCell);
  // Create the corresponding tight-binding Hamiltonian
  CreateTightBindingSuperCell(HTB, &HTBAF, AFSuperCell);
  FilePrintTightBindingDef(HTBAF, "NiO/HTBAF", "w");
  
  // Define the AF self energy
  UMeanField = 4;
  NidShift = 0; // this is related to Delta. Have fun to recrate the relation between Delta, epsilon p, epsilon d, U and NidShift...
  // Create the LDA+U Self energy
  Sigma.BlockSize=32;
  Sigma.QComplex=0;
  Sigma.N=0;
  InitBlockListOfPoles(&Sigma);
  Sigma.a0Re[ 0+32* 0]= -0.5 * UMeanField + NidShift; //t2g down A
  Sigma.a0Re[ 1+32* 1]= -0.5 * UMeanField + NidShift; //t2g down A
  Sigma.a0Re[ 2+32* 2]= -0.5 * UMeanField + NidShift; //eg  down A
  Sigma.a0Re[ 3+32* 3]= -0.5 * UMeanField + NidShift; //t2g down A
  Sigma.a0Re[ 4+32* 4]= -0.5 * UMeanField + NidShift; //eg  down A
  Sigma.a0Re[ 5+32* 5]= 0;                            //p   down A
  Sigma.a0Re[ 6+32* 6]= 0;                            //p   down A
  Sigma.a0Re[ 7+32* 7]= 0;                            //p   down A
  Sigma.a0Re[ 8+32* 8]= -0.5 * UMeanField + NidShift; //t2g down B
  Sigma.a0Re[ 9+32* 9]= -0.5 * UMeanField + NidShift; //t2g down B
  Sigma.a0Re[10+32*10]=  0.5 * UMeanField + NidShift; //eg  down B
  Sigma.a0Re[11+32*11]= -0.5 * UMeanField + NidShift; //t2g down B
  Sigma.a0Re[12+32*12]=  0.5 * UMeanField + NidShift; //eg  down B
  Sigma.a0Re[13+32*13]= 0;                            //p   down B
  Sigma.a0Re[14+32*14]= 0;                            //p   down B
  Sigma.a0Re[15+32*15]= 0;                            //p   down B
  Sigma.a0Re[16+32*16]= -0.5 * UMeanField + NidShift; //t2g up A
  Sigma.a0Re[17+32*17]= -0.5 * UMeanField + NidShift; //t2g up A
  Sigma.a0Re[18+32*18]=  0.5 * UMeanField + NidShift; //eg  up A
  Sigma.a0Re[19+32*19]= -0.5 * UMeanField + NidShift; //t2g up A
  Sigma.a0Re[20+32*20]=  0.5 * UMeanField + NidShift; //eg  up A
  Sigma.a0Re[21+32*21]= 0;                            //p   up A
  Sigma.a0Re[22+32*22]= 0;                            //p   up A
  Sigma.a0Re[23+32*23]= 0;                            //p   up A
  Sigma.a0Re[24+32*24]= -0.5 * UMeanField + NidShift; //t2g up B
  Sigma.a0Re[25+32*25]= -0.5 * UMeanField + NidShift; //t2g up B
  Sigma.a0Re[26+32*26]= -0.5 * UMeanField + NidShift; //eg  up B
  Sigma.a0Re[27+32*27]= -0.5 * UMeanField + NidShift; //t2g up B
  Sigma.a0Re[28+32*28]= -0.5 * UMeanField + NidShift; //eg  up B
  Sigma.a0Re[29+32*29]= 0;                            //p   up B
  Sigma.a0Re[30+32*30]= 0;                            //p   up B
  Sigma.a0Re[31+32*31]= 0;                            //p   up B

  if(DoBANDS)
  {
    // Calculate the band structure -  note that one can also use these functions with character plot
    BeginTiming("BANDS", &Time);
    TightBindingBANDS(HTB, Path, &BANDS);
    // Turn the eigen-values to a matrix one can plot and export the bmp
    TightBindingBANDSToMatrix(EFermi-10, EFermi+5, 1000, BANDS, &BWValue);
    Elines[0]= EFermi-10;
    Elines[1]= EFermi-5;
    Elines[2]= EFermi+0;
    Elines[3]= EFermi+5;
    TightBindingAddAxesToMatrix(EFermi-10, EFermi+5, 4, Elines, Path, &BWValue);
    FilePrintRealMatrixToBitmap(BWValue, "NiO/Bands.bmp", &RealToRGBColorBWInverse);
    EndTiming("BANDS", &Time);
    PrintTiming(Time);
  }
  
  if(DoDOS)
  {
    // Calculate the DOS - note again that one can also use these functions to calculate partial DOSes
    BeginTiming("DOS", &Time);
    // Define k-mesh
    Nk[0]=10;
    Nk[1]=10;
    Nk[2]=10;
    // define DOS range
    DOS.Emin=EFermi-10;
    DOS.Emax=EFermi+5;
    DOS.N=1500;
    DOS.Gamma=0.03;
    DOS.EminPole=DOS.Emin;
    DOS.EmaxPole=DOS.Emax;
    DOS.dE=(DOS.Emax-DOS.Emin)/DOS.N;
    SpectrumAllocate(DOS, &DOS);
    // Calculate eigenvalues
    TightBindingDOS(HTB, Nk, &DOS);
    // print to file
    FilePrintSpectra(&DOS, 1, "NiO/DOS", "w");
    EndTiming("DOS", &Time);
    PrintTiming(Time);
  }
  
  if(DoG0)
  {
    // Create G0 as a 8 by 8 matrix
    BeginTiming("Create G0", &Time);
    // Define k-mesh
    Nk[0]=10;
    Nk[1]=10;
    Nk[2]=10;
    TightBindingHamiltonianToBlockListOfPoles(HTB, Nk, &G0, EFermi-10, EFermi+5, 1500);
    BlockListOfPolesToSpectra(G0, G0Spec, EFermi-10, EFermi+5, 1500);
    FilePrintSpectra(G0Spec, HTB.NFermion*HTB.NFermion, "NiO/G0", "w");
    EndTiming("Create G0", &Time);
    PrintTiming(Time);
  }
  
  if(DoG0AF)
  {
    // Create G0AF as a 16 by 16 matrix
    // It only includes one spin species
    BeginTiming("Create G0AF", &Time);
    // Define k-mesh
    Nk[0]=10;
    Nk[1]=10;
    Nk[2]=10;
    TightBindingHamiltonianToBlockListOfPoles(HTBAF, Nk, &G0AF, EFermi-10, EFermi+5, 1500);
    BlockListOfPolesToSpectra(G0AF, G0AFSpec, EFermi-10, EFermi+5, 1500);
    FilePrintSpectra(G0AFSpec, HTBAF.NFermion*HTBAF.NFermion, "NiO/G0AF", "w");
    EndTiming("Create G0AF", &Time);
    PrintTiming(Time);
  }
  
  if(DoGAF)
  {
    // Create GAF as a 32 by 32 matrix
    // It includes two sites (unit cells) and spin up/dn
    BeginTiming("Create GAF", &Time);
    // Define k-mesh
    Nk[0]=10;
    Nk[1]=10;
    Nk[2]=10;
    TightBindingHamiltonianAddSpinPlusSelfEnergyToBlockListOfPoles(HTBAF, Nk, Sigma, &GAF, EFermi-10, EFermi+5, 250); //change 250 to 1500 for consistence (does not work with debuger on)
    BlockListOfPolesToSpectra(GAF, GAFSpec, EFermi-10, EFermi+5, 1500);
    FilePrintSpectra(GAFSpec, 4*HTBAF.NFermion*HTBAF.NFermion, "NiO/GAF", "w");
    EndTiming("Create GAF", &Time);
    PrintTiming(Time);
  }

  if(DoHAF)
  {
    // Given GAF one can create an impurity Hamiltonian
    //
    // I use a basis of natural impurity orbitals, defined by Sigma
    //
    // If one wants to add locally the full interactions one needs to subtract the potentials added via the self energy and add the Coulomb interactions
    //   Note that the ground-state of HAF has fully occupied valence and fully empty conduction orbitals. The many particle basis is no larger than the many particle
    //   basis of a ligand field theory calculation. This will remain valid if one takes the many body interactions and the self energy such that the local and non-local density matrix is the same
    //   i.e. one needs to (one can) puzzle the self energy such taht the number of determinants needed does not change if one replaces the self energy by the full interactions.
    //   Have fun puzzling what sigma truely should be.
    //
    //  The basis of HAF is
    //     0 --  9                : The d-shell (0--4 spin up 5--9 spin down)
    //    10 -- 19                : a ligand shell such that if the d-occupation is 8.2 the Ligand occupation is 1.8 (total occupation of d and L shell is 10)
    //    20 -- 20 + 10 NBlockVal : NBlockVal valence shells, which are occupied
    //    ?? -- end               : NBlockCon conduction shells which are empty
    NBlockVal=10;  // number of Valence bath orbitals
    NBlockCon=4;   // number of Conduction bath orbitals
    BeginTiming("Create HAF", &Time);
    // convert GAF (BlockListOfPoles) To impurity Hamiltonian
    if(!DoGAF)
    {
      printf("\n\n\nSorry I can not create the Hamiltonian for an AF lattice if I did not calculate the Green's function first\n");
      exit(1);
    }
    // The Green's function has both the d and p orbitals in the basis, for site A and B
    // I want a local Hamiltonian with only the d orbitals of site A in the basis
    // Matrix AT defines the basis transformation
    AT.Ni=10;
    AT.Nj=32;
    AT.QComplex=0;
    InitCompactMatrix(&AT);
    AT.VRe[0 +  0*AT.Nj]=1;
    AT.VRe[1 +  1*AT.Nj]=1;
    AT.VRe[2 +  2*AT.Nj]=1;
    AT.VRe[3 +  3*AT.Nj]=1;
    AT.VRe[4 +  4*AT.Nj]=1;
    AT.VRe[16 + 5*AT.Nj]=1;
    AT.VRe[17 + 6*AT.Nj]=1;
    AT.VRe[18 + 7*AT.Nj]=1;
    AT.VRe[19 + 8*AT.Nj]=1;
    AT.VRe[20 + 9*AT.Nj]=1;
    GAF.mu=EFermi;
    BlockListOfPolesToBlockNaturalImpurityOrbitalMatrixReduceDimension(&GAF, &GAFNat, &AT, kBT, NBlockVal, NBlockCon);
    BlockNaturalImpurityOrbitalMatrixChop(&GAFNat, 1000*DBL_EPSILON);
    //PrintBlockNaturalImpurityOrbitalMatrix(GAFNat);
    BlockNaturalImpurityOrbitalMatrixToImpurityHamiltonian(&HAF, &GAFNat, NBlockVal, NBlockCon, 0, 0, 0);
    OperatorChop(&HAF, 1000*DBL_EPSILON);
    //PrintOperatorReadable(HAF);
    BlockNaturalImpurityOrbitalToSpectra(&GAFNat, SpecNat, EFermi-10, EFermi+5, 1500, 0.01);
    FilePrintSpectra(SpecNat, 10*10, "NiO/GAFd", "w");
    EndTiming("Create HAF", &Time);
  }
  
  if(DoGrd)
  {
    BeginTiming("Create Ground state", &Time);
    // Calculate the ground-state
    if(!DoGAF)
    {
      printf("\n\n\nSorry I can not create the Ground-state for an AF lattice if I did not calculate the Hamiltonian first\n");
      exit(1);
    }
    // We need to define how many electrons there are in the system
    // The current algorithm that finds the ground-state is not so nice and needs a reasonable starting point
    NParticle = (NBlockVal+1)*10;
    CreateString = malloc(NParticle*sizeof(uint16_t));
    for(i=0;i<10;i++)
    {
      CreateString[i]=i;
    }
    for(i=10;i<NParticle;i++)
    {
      CreateString[i]=i+10;
    }
    CreateWaveFunctionFromString(&psi, CreateString, NParticle, HAF.NFermion, HAF.NBoson, NBitsKey);
    
    PrintOperatorReadable(HAF);
    PrintWaveFunctionReadable(psi);
    
    // I know that in the end the conduction and valence orbitals are eiter completely filled or empty.
    // In order to help convergence I first require this to be true.
    Restrictions.N=1;
    Restrictions.Min=malloc(Restrictions.N*sizeof(uint32_t));
    Restrictions.Max=malloc(Restrictions.N*sizeof(uint32_t));
    Restrictions.Determinant=malloc(Restrictions.N*sizeof(uint8_t*));
    for(i=0;i<Restrictions.N;i++)
    {
      Restrictions.Determinant[i]=calloc(psi.NWord,sizeof(uint8_t));
    }
    Restrictions.Min[0]=10;
    Restrictions.Max[0]=10;
    Restrictions.Determinant[0][0] = 0xFF;
    Restrictions.Determinant[0][1] = 0xFF;
    Restrictions.Determinant[0][2] = 0xF0;
    
    // start with low accuracy
    epsilon=0.01;
    zero=100*DBL_EPSILON;
    
    LanczosGroundStateRestricted(HAF, &psi, Restrictions, &epsilon, zero, 0, &psi, 0);
    PrintWaveFunctionDensity(psi);
    //PrintWaveFunctionReadable(psi);
    
    
    // increase accuracy
    epsilon=0.0001;
    zero=100*DBL_EPSILON;
    
    LanczosGroundStateRestricted(HAF, &psi, Restrictions, &epsilon, zero, 0, &psi, 0);
    //PrintWaveFunctionDensity(psi);
    //PrintWaveFunctionReadable(psi);
    
    //LanczosGroundState(HAF, &psi, &epsilon, zero, 0, &psi, 0);
    //PrintWaveFunctionDensity(psi);
    //PrintWaveFunctionReadable(psi);
    
    
    // lifting the restrictions does not change the occupation of the valence and conduction orbitals.
    
    LanczosGroundState(HAF, &psi, &epsilon, zero, 0, &psi, 0);
    PrintWaveFunctionDensity(psi);
    //PrintWaveFunctionReadable(psi);
    
    EndTiming("Create Ground state", &Time);
    
    PsiHermitianOperatorPsi(HAF, psi, &energy);
  }
  
  if(DoPES)
  {
    SpecTri.N=15;
    InitTriDiagonalMatrix(&SpecTri);
    epsilon=0.0001;
    zero=100*DBL_EPSILON;
    
    Restrictions.Min[0]=9;
    Restrictions.Max[0]=10;
    //for(i=0;i<psi.NWord;i++)
    // {
    //  Restrictions.Determinant[0][i]=0;
    //}
    //for(i=20+10*NBlockVal;i<10*NBlockVal+10*NBlockCon+20;i++)
    //{
    //  Restrictions.Determinant[0][i>>3]+=128>>(i&7);
    //}
    for(i=0;i<psi.NWord;i++)
    {
      printf("%x ",Restrictions.Determinant[0][i]);
    }
    printf("\n");
    
    BeginTiming("Calculate PES", &Time);
    if(!DoGAF)
    {
      printf("\n\n\nSorry I can not calculate the PES spectrum if I did not calculate the ground-state first\n");
      exit(1);
    }
    for(i=0;i<10;i++)
    {
      CreateOperatorAnnihilate(HAF.NFermion, (uint16_t)i, &(TPES[i]));
    }
    for(i=0;i<10;i++)
    {
      SpectraPES[i].N=6000;
      SpectraPES[i].Emin=-EFermi-5.0;
      SpectraPES[i].Emax=-EFermi+10.0;
      SpectraPES[i].dE = (SpectraPES[i].Emax - SpectraPES[i].Emin) / SpectraPES[i].N;
      SpectraPES[i].Gamma = 0.01;
      SpectrumAllocate(SpectraPES[i], &SpectraPES[i]);
    }
    
    for(i=0;i<10;i++)
    {
      BeginTiming("spectrum PES", &Time);
      OperatorPsi(TPES[i], psi, &Opsi);
      LanczosTriDiagonalizeRestricted(HAF, Opsi, Restrictions, epsilon, &SpecTri,  SpectraPES[i].Gamma);
      TriDiagonalMatrixToSpectrum(energy, SpecTri, &SpectraPES[i]);
      EndTiming("spectrum PES", &Time);
    }
    
    FilePrintSpectra(&SpectraPES[0], 10, "NiO/SpecPES","w");
    FreeTriDiagonalMatrix(&SpecTri);
    EndTiming("Calculate PES", &Time);
    PrintTiming(Time);
  }
  
  if(DoIPES)
  {
    SpecTri.N=15;
    InitTriDiagonalMatrix(&SpecTri);
    epsilon=0.0001;
    zero=100*DBL_EPSILON;
    
    Restrictions.Min[0]=10;
    Restrictions.Max[0]=11;
    
    BeginTiming("Calculate IPES", &Time);
    if(!DoGAF)
    {
      printf("\n\n\nSorry I can not calculate the IPES spectrum if I did not calculate the ground-state first\n");
      exit(1);
    }
    for(i=0;i<10;i++)
    {
      CreateOperatorCreate(HAF.NFermion, (uint16_t)i, &(TIPES[i]));
    }
    for(i=0;i<10;i++)
    {
      SpectraIPES[i].N=6000;
      SpectraIPES[i].Emin=EFermi-5.0;
      SpectraIPES[i].Emax=EFermi+10.0;
      SpectraIPES[i].dE = (SpectraIPES[i].Emax - SpectraIPES[i].Emin) / SpectraIPES[i].N;
      SpectraIPES[i].Gamma = 0.01;
      SpectrumAllocate(SpectraIPES[i], &SpectraIPES[i]);
    }
    
    for(i=0;i<10;i++)
    {
      BeginTiming("spectrum IPES", &Time);
      OperatorPsi(TIPES[i], psi, &Opsi);
      LanczosTriDiagonalizeRestricted(HAF, Opsi, Restrictions, epsilon, &SpecTri,  SpectraIPES[i].Gamma);
      TriDiagonalMatrixToSpectrum(energy, SpecTri, &SpectraIPES[i]);
      EndTiming("spectrum IPES", &Time);
    }
    
    FilePrintSpectra(&SpectraIPES[0], 10, "NiO/SpecIPES","w");
    FreeTriDiagonalMatrix(&SpecTri);
    EndTiming("Calculate IPES", &Time);
    PrintTiming(Time);
  }
  
   
}