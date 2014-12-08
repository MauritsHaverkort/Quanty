static int LuaHelp(lua_State *L)
{
  uint32_t len = (uint32_t)lua_gettop(L);
  if(len>1) return luaL_error(L, "Error in .Help(arg) expected 0 (general help) or 1 (help on function) arguments\n");
  if(len==0)
  {
    printf("-------------------------------------------------------------------------\n");
    printf("-          solid state physics model library  - lua interface           -\n");
    printf("-          (c) Maurits W. Haverkort                                     -\n");
    printf("-          http://www.solstatphys.org                                   -\n");
    printf("-------------------------------------------------------------------------\n");
    printf("\n\n");
    printf("Available functions\n");
    printf("\n");
    printf("ssp.Help()\n");
    printf("  - shows this file\n");
    printf("  - excepts as input one of the other functions in order to show help on that function\n");
    printf("  - example: ssp.Help(\"ssp.NewOperator\");\n");
    printf("\n");
    printf("ssp.NewOperator()\n");
    printf("  - creates a new operator\n");
    printf("  - expects as input on the first position a string describing the type of operator\n");
    printf("  - other parameters depend on the type of operator defined in the first string\n");
    printf("  - example: ssp.NewOperator(\"Sz\", NFermion, IndexUp, IndexDn);\n");
    printf("\n");
    printf("ssp.PotentialExpandedOnYlm(\"Oh\",2,{0.6,-0.4})\n");
    printf("  - creates a list representing a potential expanded on spherical Harmonics for a given symmetry and eigen energy of the irreducible representations\n");
    printf("  - The potential V(theta, phi) is given as:\n");
    printf("  - V(theta, phi) = sum_k,m A_{k,m} Y_km(theta,phi)\n");
    printf("  - the list saved is {{k1,m1,Ak1m1},{k2,m2,Ak2m2},...}\n");
    printf("  - input is the point-group symmetry, the angular momentum and the eigen energies of the irreducible representations\n");
    printf("  - ssp.PotentialExpandedOnYlm(\"Oh\",2,{0.6,-0.4});\n");
    printf("  - for a potential in cubic symmetry acting on the d-shell with the energy of the eg's as 0.6 and t2g's at -0.4 (10Dq=1.0)\n");
    printf("\n");
    printf("ssp.NewWavefunction()\n");
    printf("  - creates a new wavefunction\n");
    printf("  - expects as input the number of fermions, number of bosons, a table of determinants with the prefactor\n");
    printf("  - example: ssp.NewWavefunction(10,0,{{{0x80,0x00},1},{{0x04,0x00},1}});\n");
    printf("  - example: ssp.NewWavefunction(10,0,{{\"1010101010\",1},{\"0101010101\",1}});\n");
    printf("\n");
    printf("ssp.Copy() or ssp.Clone()\n");
    printf("  - creats a copy of a wavefunction or operator\n");
    printf("  - expects as input an object, returns a copy of that object\n");
    printf("  - example: ssp.Copy(psi);\n");
    printf("\n");
    printf("ssp.Conjugate() or ssp.C()\n");
    printf("  - creats a copy of a wavefunction or operator and conjugates this copy\n");
    printf("  - expects as input an object, returns a copy of that object\n");
    printf("  - example: ssp.Conjugate(psi);\n");
    printf("\n");
    printf("ssp.ConjugateTranspose() or ssp.CT()\n");
    printf("  - creats a copy of a wavefunction or operator and conjugates transposes this copy\n");
    printf("  - expects as input an object, returns a copy of that object\n");
    printf("  - example: ssp.ConjugateTranspose(psi);\n");
    printf("\n");
    printf("ssp.Transpose() or ssp.T()\n");
    printf("  - creats a copy of a wavefunction or operator and transposes this copy\n");
    printf("  - expects as input an object, returns a copy of that object\n");
    printf("  - example: ssp.Transpose(psi);\n");
    printf("ssp.Rotate()\n");
    printf("  - Given an operator on a basis of N fermions + M Bosons\n");
    printf("  - and a matrix of size N+M\n");
    printf("  - creates a new operator on a rotated basis\n");
    printf("  - example: ssp.Rotate(Operator,RotationMatrix);\n");
    printf("\n");
    printf("ssp.EvolveToGroundstate()\n");
    printf("  - calculates eigenfunctions of the Hamiltonian, starting from the specified state\n");
    printf("  - expects as input an operator and a (list of) wavefunction(s), a list of options on the last position is allowed\n");
    printf("  - example ssp.EvolveToGroundstate(OppSy,{psi1,psi2});\n");
    printf("ssp.Eigensystem()\n");
    printf("  - calculates the lowest N eigenfunctions of a given operator\n");
    printf("  - expects as input an operator and either a list of starting functions or a set of restrictions specifying the starting point\n");
    printf("  - example ssp.Eigensystem(OppSy,{psi1,psi2});\n");
    printf("  - example ssp.Eigensystem(OppSy,Restrictions,Nspi);\n");
    printf("\n");
    printf("ssp.CreateSpectra()\n");
    printf("  - calculates spectra of a giving Hamiltonian, starting from the specified state and using a specific transition operator\n");
    printf("  - expects as input an operator a (list of) transition operator(s) and a (list of) wavefunction(s), a list of options on the last position is allowed\n");
    printf("  - example ssp.CreateSpectra(OppH,{T1,T2},{psi1,psi2});\n");
    printf("\n");
    printf("ssp.CreateFluorescenceYield()\n");
    printf("  - calculates the fluorescence yield spectra\n");
    printf("  - calculated is < i | T1^+ G^+ T2 T2^+ G T1 | i>\n");
    printf("  - expects as input an operator a (list of) transition operator(s) T1 a (list of) transition operator(s) T2^+ and a (list of) wavefunction(s), a list of options on the last position is allowed\n");
    printf("  - example ssp.CreateFluorescenceYield(OppH,{T1a,T1b},{T2a,T2b},{psi1,psi2});\n");
    printf("\n");
    
    printf("\n\n");
    printf("Available objects\n");
    printf("\n");
    printf("operator\n");
    printf("  - stores a quantum mechanical operator based on second quantization\n");
    printf("\n");
    printf("wavefunction\n");
    printf("  - stores a wavefunction\n");
    printf("\n");
    printf("complex\n");
    printf("  - stores a complex number\n");
    printf("\n");
    printf("spectra\n");
    printf("  - stores a spectral function\n");
    printf("\n");
    
  }
  if(len==1)
  {
    char const *Help =lua_tostring(L, 1);
    if (strcmp(Help, "ssp.NewOperator") == 0)
    {
      printf("ssp.NewOperator()\n");
      printf("  - creates a new operator\n");
      printf("  - expects as input on the first position a string describing the type of operator\n");
      printf("  - other parameters depend on the type of operator defined in the first string\n");
      printf("  - possible operators and aditional parameters are:\n");
      printf("     - \"Coulomb\" or \"U\"\n");
      printf("       - Lindex Up : array of length 2l+1 of integers defining the position of ml=-l to ml=l with spin up\n");
      printf("       - Lindex Dn : array of length 2l+1 of integers defining the position of ml=-l to ml=l with spin down\n");
      printf("       - F : Slater integrals\n");
      printf("     - \"Coulomb\" or \"U\"\n");
      printf("       - Lindex Up : array of length 2lA+1 of integers defining the position of mlA=-lA to mlA=lA with spin up\n");
      printf("       - Lindex Dn : array of length 2lA+1 of integers defining the position of mlA=-lA to mlA=lA with spin down\n");
      printf("       - Lindex Up : array of length 2lB+1 of integers defining the position of mlB=-lB to mlB=lB with spin up\n");
      printf("       - Lindex Dn : array of length 2lB+1 of integers defining the position of mlB=-lB to mlB=lB with spin down\n");
      printf("       - F : Slater integrals\n");
      printf("       - G : Slater integrals\n");
      printf("     - \"Sx\"\n");
      printf("       - Lindex Up : array of length 2l+1 of integers defining the position of ml=-l to ml=l with spin up\n");
      printf("       - Lindex Dn : array of length 2l+1 of integers defining the position of ml=-l to ml=l with spin down\n");
      printf("     - \"Sy\"\n");
      printf("       - Lindex Up : array of length 2l+1 of integers defining the position of ml=-l to ml=l with spin up\n");
      printf("       - Lindex Dn : array of length 2l+1 of integers defining the position of ml=-l to ml=l with spin down\n");
      printf("     - \"Sz\"\n");
      printf("       - Lindex Up : array of length 2l+1 of integers defining the position of ml=-l to ml=l with spin up\n");
      printf("       - Lindex Dn : array of length 2l+1 of integers defining the position of ml=-l to ml=l with spin down\n");
      printf("     - \"Splus\"\n");
      printf("       - Lindex Up : array of length 2l+1 of integers defining the position of ml=-l to ml=l with spin up\n");
      printf("       - Lindex Dn : array of length 2l+1 of integers defining the position of ml=-l to ml=l with spin down\n");
      printf("     - \"Smin\"\n");
      printf("       - Lindex Up : array of length 2l+1 of integers defining the position of ml=-l to ml=l with spin up\n");
      printf("       - Lindex Dn : array of length 2l+1 of integers defining the position of ml=-l to ml=l with spin down\n");
      printf("     - \"Ssqr\"\n");
      printf("       - Lindex Up : array of length 2l+1 of integers defining the position of ml=-l to ml=l with spin up\n");
      printf("       - Lindex Dn : array of length 2l+1 of integers defining the position of ml=-l to ml=l with spin down\n");
      printf("     - \"Lx\"\n");
      printf("       - Lindex Up : array of length 2l+1 of integers defining the position of ml=-l to ml=l with spin up\n");
      printf("       - Lindex Dn : array of length 2l+1 of integers defining the position of ml=-l to ml=l with spin down\n");
      printf("     - \"Ly\"\n");
      printf("       - Lindex Up : array of length 2l+1 of integers defining the position of ml=-l to ml=l with spin up\n");
      printf("       - Lindex Dn : array of length 2l+1 of integers defining the position of ml=-l to ml=l with spin down\n");
      printf("     - \"Lz\"\n");
      printf("       - Lindex Up : array of length 2l+1 of integers defining the position of ml=-l to ml=l with spin up\n");
      printf("       - Lindex Dn : array of length 2l+1 of integers defining the position of ml=-l to ml=l with spin down\n");
      printf("     - \"Lplus\"\n");
      printf("       - Lindex Up : array of length 2l+1 of integers defining the position of ml=-l to ml=l with spin up\n");
      printf("       - Lindex Dn : array of length 2l+1 of integers defining the position of ml=-l to ml=l with spin down\n");
      printf("     - \"Lmin\"\n");
      printf("       - Lindex Up : array of length 2l+1 of integers defining the position of ml=-l to ml=l with spin up\n");
      printf("       - Lindex Dn : array of length 2l+1 of integers defining the position of ml=-l to ml=l with spin down\n");
      printf("     - \"Lsqr\"\n");
      printf("       - Lindex Up : array of length 2l+1 of integers defining the position of ml=-l to ml=l with spin up\n");
      printf("       - Lindex Dn : array of length 2l+1 of integers defining the position of ml=-l to ml=l with spin down\n");
      printf("     - \"Jx\"\n");
      printf("       - Lindex Up : array of length 2l+1 of integers defining the position of ml=-l to ml=l with spin up\n");
      printf("       - Lindex Dn : array of length 2l+1 of integers defining the position of ml=-l to ml=l with spin down\n");
      printf("     - \"Jy\"\n");
      printf("       - Lindex Up : array of length 2l+1 of integers defining the position of ml=-l to ml=l with spin up\n");
      printf("       - Lindex Dn : array of length 2l+1 of integers defining the position of ml=-l to ml=l with spin down\n");
      printf("     - \"Jz\"\n");
      printf("       - Lindex Up : array of length 2l+1 of integers defining the position of ml=-l to ml=l with spin up\n");
      printf("       - Lindex Dn : array of length 2l+1 of integers defining the position of ml=-l to ml=l with spin down\n");
      printf("     - \"Jplus\"\n");
      printf("       - Lindex Up : array of length 2l+1 of integers defining the position of ml=-l to ml=l with spin up\n");
      printf("       - Lindex Dn : array of length 2l+1 of integers defining the position of ml=-l to ml=l with spin down\n");
      printf("     - \"Jmin\"\n");
      printf("       - Lindex Up : array of length 2l+1 of integers defining the position of ml=-l to ml=l with spin up\n");
      printf("       - Lindex Dn : array of length 2l+1 of integers defining the position of ml=-l to ml=l with spin down\n");
      printf("     - \"Jsqr\"\n");
      printf("       - Lindex Up : array of length 2l+1 of integers defining the position of ml=-l to ml=l with spin up\n");
      printf("       - Lindex Dn : array of length 2l+1 of integers defining the position of ml=-l to ml=l with spin down\n");
      printf("     - \"ldots\"\n");
      printf("       - Lindex Up : array of length 2l+1 of integers defining the position of ml=-l to ml=l with spin up\n");
      printf("       - Lindex Dn : array of length 2l+1 of integers defining the position of ml=-l to ml=l with spin down\n");
      printf("     - \"CF\"\n");
      printf("       - Lindex : array of length 2l+1 of integers defining the position of ml=-l to ml=l\n");
      printf("       - {{k1,k2,Ak1m1}, {k2,m2,Ak2m2}, ...} : crystal field expanded on spherical harmonics\n");
      printf("     - \"CF\"\n");
      printf("       - Lindex Up: array of length 2l+1 of integers defining the position of ml=-l to ml=l with spin up\n");
      printf("       - Lindex Dn: array of length 2l+1 of integers defining the position of ml=-l to ml=l with spin down\n");
      printf("       - {{k1,k2,Ak1m1}, {k2,m2,Ak2m2}, ...} : crystal field expanded on spherical harmonics\n");
      printf("     - \"CF\"\n");
      printf("       - Lindex Up A: array of length 2lA+1 of integers defining the position of mlA=-lA to mlA=lA with spin up\n");
      printf("       - Lindex Dn A: array of length 2lA+1 of integers defining the position of mlA=-lA to mlA=lA with spin down\n");
      printf("       - Lindex Up B: array of length 2lB+1 of integers defining the position of mlB=-lB to mlB=lB with spin up\n");
      printf("       - Lindex Dn B: array of length 2lB+1 of integers defining the position of mlB=-lB to mlB=lB with spin down\n");
      printf("       - {{k1,k2,Ak1m1}, {k2,m2,Ak2m2}, ...} : crystal field expanded on spherical harmonics\n");
      printf("     - \"Create\" or \"Cr\"\n");
      printf("       - NFermion\n");
      printf("       - index (if index smaller than NFermion creates a Fermion, otherwise a boson\n");
      printf("     - \"Annihilate\" or \"An\"\n");
      printf("       - NFermion\n");
      printf("       - index (if index smaller than NFermion annihilates a Fermion, otherwise a boson\n");
      printf("     - \"Number\" or \"CrAn\"\n");
      printf("       - NFermion\n");
      printf("       - index (if index smaller than NFermion creates a Fermion, otherwise a boson\n");
      printf("       - index (if index smaller than NFermion annihilates a Fermion, otherwise a boson\n");
      printf("       - \n");
      printf("\n");
      printf("  - example for a d state:\n");
      printf("  - NFermion=10;\n");
      printf("  - IndexUp={0,1,2,3,4};\n");
      printf("  - IndexDn={5,6,7,8,9};\n");
      printf("  - F0=20; F2=10; F4=6.25;\n");
      printf("  - F={F0,F2,F4}");
      printf("  - OppSz = ssp.NewOperator(\"Sz\", IndexUp, IndexDn);\n");
      printf("  - OppU  = ssp.NewOperator(\"U\", IndexUp, IndexDn, F);\n");
      printf("\n");
      printf("\n");
      printf("\n");
    }
    if (strcmp(Help, "ssp.NewWavefunction") == 0)
    {
      printf("ssp.NewWavefunction()\n");
      printf("  - creates a new wavefunction\n");
      printf("  - expects as input the number of fermions, number of bosons, a table of determinants with the prefactor\n");
      printf("  - an optional 4-th input number is NBitsKey\n");
      printf("    - number of elements in the Hash lookup table. (2^NBitsKey)\n");
      printf("    - This number should be larger than the number of Determinants in the basis for optimal preformance\n");
      printf("    - Standard value is 16 (65 536) which works good for small calculations upto basis sets of 25 000\n");
      printf("    - If more determinants are needed it is recommanded to enlarge\n");
      printf("    - 20 (  1 048 576) for basis sets upto    500 000\n");
      printf("    - 24 ( 16 777 216) for basis sets upto 10 000 000\n");
      printf("    - for even larger basis sets argitecture depenent tests are advised.\n");
      printf("    - note that results are independent of this value, only calculation time and memory requirement changes.\n");
      printf("\n");
      printf("  - example: A ssp.NewWavefunction(10,0,{{\"1010101010\",1}});\n");
      printf("\n");
      printf("\n");
    }
    if ((strcmp(Help, "ssp.Copy") == 0)||(strcmp(Help, "ssp.Clone") == 0))
    {
      printf("ssp.Copy() or ssp.Clone()\n");
      printf("  - creats a copy of a wavefunction or operator\n");
      printf("  - expects as input an object, returns a copy of that object\n");
      printf("  - note that \"=\" often does not what you want.\n");
      printf("  - A=B for operators and wave-functions copies a pointer to the object.\n");
      printf("  - A change in A thus automatically changes B\n");
      printf("\n");
      printf("  - example: A=ssp.Copy(B);\n");
      printf("\n");
      printf("\n");
    }
    if ((strcmp(Help, "ssp.Conjugate") == 0)||(strcmp(Help, "ssp.C") == 0))
    {
      printf("ssp.Conjugate() or ssp.C()\n");
      printf("  - creats a copy of a wavefunction or operator and conjugates this copy\n");
      printf("  - expects as input an object, returns a copy of that object\n");
      printf("\n");
      printf("  - example: A=ssp.Conjugate(B);\n");
      printf("\n");
      printf("\n");
    }
    if ((strcmp(Help, "ssp.ConjugateTranspose") == 0)||(strcmp(Help, "ssp.CT")==0))
    {
      printf("ssp.ConjugateTranspose() or ssp.CT()\n");
      printf("  - creats a copy of a wavefunction or operator and conjugates transposes this copy\n");
      printf("  - expects as input an object, returns a copy of that object\n");
      printf("\n");
      printf("  - example: A=ssp.ConjugateTranspose(B);\n");
      printf("\n");
      printf("\n");
    }
    if ((strcmp(Help, "ssp.Transpose") == 0)||(strcmp(Help, "ssp.T") == 0))
    {
      printf("ssp.Transpose() or ssp.T()\n");
      printf("  - creats a copy of a wavefunction or operator and transposes this copy\n");
      printf("  - expects as input an object, returns a copy of that object\n");
      printf("\n");
      printf("  - example: A=ssp.Transpose(B);\n");
      printf("\n");
      printf("\n");
    }
    if (strcmp(Help, "ssp.Rotate") == 0)
    {
      printf("ssp.Rotate()\n");
      printf("\n");
      printf("\n");
    }
    if (strcmp(Help, "ssp.EvolveToGroundstate") == 0)
    {
      printf("ssp.EvolveToGroundstate()\n");
      printf("  - calculates eigenfunctions of the Hamiltonian, starting from the specified state\n");
      printf("  - expects as input an operator and a (list of) wavefunction(s), a list of options on the last position is allowed\n");
      printf("  - The function is intended to return the eigenstate of the Hamiltonian that has the largest overlap with the starting functions\n");
      printf("  - the returned function will for sure be of the same irreducible representation, but it favoirs low lying states with smaller overlaps to high lying states\n");
      printf("\n");
      printf("  - options are:\n");
      printf("    - \"zero\"\n");
      printf("      - Determins how far and good the states are converged. Standard setting 100*DBL_EPSILON, should not be smaller than DBL_EPSILON\n");
      printf("    - \"epsilon\"\n");
      printf("      - Determins the smallest prefactor of a determinant that is kept in the wave-function. Standard setting 10000*DBL_EPSILON\n");
      printf("    - \"restrictions\"\n");
      printf("      - restricts the basis by excluding certain occupations\n");
      printf("      - Example for a ligand field calculation one can require to have no more than 1 hole in the Ligand shell (between 9 and 10 electrons)\n");
      printf("      - {\"restrictions\",{20,0,{{\"0000000000 1111111111\",9,10}}}}\n");
      printf("\n");
      printf("  - example: ssp.EvolveToGroundstate(H,{psi1,psi2},{{\"epsilon\",0.000001},{\"restrictions\",{10,0,{{\"11111 00000\",5,5}}}});\n");
      printf("\n");
      printf("\n");
    }
    if (strcmp(Help, "ssp.Eigensystem") == 0)
    {
      printf("ssp.Eigensystem()\n");
      printf("  - calculates the lowest eigen functions of a given Hamiltonian. The number of states to calculate is given as an input\n");
      printf("  - expects as input an operator and either\n");
      printf("    - a (list of) wavefunction(s) that specify a starting point. (filling)\n");
      printf("    - note that the final state found might be of different symmetry than the starting point\n");
      printf("    - note also that the lowest state might not be found if the starting point excludes certain symmetries\n");
      printf("    - look at ssp.EvolveToGroundstate if one wants to track a specific state\n");
      printf("  - or\n");
      printf("    - a set of restrictions defining the starting point (number of electrons, filled core shells)\n");
      printf("    - the number of lowest eigen states to calculate");
      printf("  - aditional options can be given\n");
      printf("\n");
      printf("  - options are:\n");
      printf("    - \"zero\"\n");
      printf("      - Determins how far and good the states are converged. Standard setting 100*DBL_EPSILON, should not be smaller than DBL_EPSILON\n");
      printf("    - \"epsilon\"\n");
      printf("      - Determins the smallest prefactor of a determinant that is kept in the wave-function. Standard setting 10000*DBL_EPSILON\n");
      printf("    - \"restrictions\"\n");
      printf("      - restricts the basis by excluding certain occupations\n");
      printf("      - Example for a ligand field calculation one can require to have no more than 1 hole in the Ligand shell (between 9 and 10 electrons)\n");
      printf("      - {\"restrictions\",{20,0,{{\"0000000000 1111111111\",9,10}}}}\n");
      printf("    - \"NBitsKey\"\n");
      printf("      + only when a set of restricitons and Npsi is given as an input, otherwise NBitsKey is taken from the wave function\n");
      printf("      - number of elements in the Hash lookup table. (2^NBitsKey)\n");
      printf("      - This number should be larger than the number of Determinants in the basis for optimal preformance\n");
      printf("      - Standard value is 16 (65 536) which works good for small calculations upto basis sets of 25 000\n");
      printf("      - If more determinants are needed it is recommanded to enlarge\n");
      printf("      - 20 (  1 048 576) for basis sets upto    500 000\n");
      printf("      - 24 ( 16 777 216) for basis sets upto 10 000 000\n");
      printf("      - for even larger basis sets argitecture depenent tests are advised.\n");
      printf("      - note that results are independent of this value, only calculation time and memory requirement changes.\n");
      printf("\n");
      printf("\n");
      printf("\n");
    }
    if (strcmp(Help, "ssp.CreateSpectra") == 0)
    {
      printf("ssp.CreateSpectra()\n");
      printf("  - calculates spectra of a giving Hamiltonian, starting from the specified state and using a specific transition operator\n");
      printf("  - expects as input an operator a (list of) transition operator(s) and a (list of) wavefunction(s), a list of options on the last position is allowed\n");
      printf("\n");
      printf("  - options are:\n");
      printf("    - \"epsilon\"\n");
      printf("      - Determins the smallest prefactor of a determinant that is kept in the wave-function. Standard setting 10000*DBL_EPSILON\n");
      printf("    - \"restrictions\"\n");
      printf("      - restricts the basis by excluding certain occupations\n");
      printf("      - Example for a ligand field calculation one can require to have no more than 1 hole in the Ligand shell (between 9 and 10 electrons)\n");
      printf("      - {\"restrictions\",{20,0,{{\"0000000000 1111111111\",9,10}}}}\n");
      printf("    - \"NE\"\n");
      printf("      - number of energy points in the spectrum\n");
      printf("    - \"Emin\"\n");
      printf("      - smallest energy in the spectrum\n");
      printf("    - \"Emax\"\n");
      printf("      - largest energy in the spectrum\n");
      printf("    - \"Gamma\"\n");
      printf("      - Lorenzian broadening of the spectra (full with half maximum FWHM)\n");
      printf("\n");
      printf("\n");
    }
    if (strcmp(Help, "ssp.CreateFluorescenceYield") == 0)
    {
      printf("ssp.CreateFluorescenceYield()\n");
      printf("  - calculates fluorescence yield spectra of a giving Hamiltonian, starting from the specified state and using a specific transition operator\n");
      printf("  - < i | T1^+ G^+ T2 T2^+ G T1 | i >\n");
      printf("  - expects as input an operator a (list of) transition operator(s) T1 a (list of) transition operator(s) T2^+ and a (list of) wavefunction(s), a list of options on the last position is allowed\n");
      printf("\n");
      printf("  - options are:\n");
      printf("    - \"epsilon\"\n");
      printf("      - Determins the smallest prefactor of a determinant that is kept in the wave-function. Standard setting 10000*DBL_EPSILON\n");
      printf("    - \"restrictions\"\n");
      printf("      - restricts the basis by excluding certain occupations\n");
      printf("      - Example for a ligand field calculation one can require to have no more than 1 hole in the Ligand shell (between 9 and 10 electrons)\n");
      printf("      - {\"restrictions\",{20,0,{{\"0000000000 1111111111\",9,10}}}}\n");
      printf("    - \"NE\"\n");
      printf("      - number of energy points in the spectrum\n");
      printf("    - \"Emin\"\n");
      printf("      - smallest energy in the spectrum\n");
      printf("    - \"Emax\"\n");
      printf("      - largest energy in the spectrum\n");
      printf("    - \"Gamma\"\n");
      printf("      - Lorenzian broadening of the spectra (full with half maximum FWHM)\n");
      printf("\n");
      printf("\n");
    }
    if (strcmp(Help, "operator") == 0)
    {
      printf("operator\n");
      printf("  - stores a quantum mechanical operator based on second quantization\n");
      printf("\n");
      printf("\n");
      printf("\n");
    }
    if (strcmp(Help, "wavefunction") == 0)
    {
      printf("wavefunction\n");
      printf("  - stores a wavefunction\n");
      printf("\n");
      printf("\n");
      printf("\n");
    }
    if (strcmp(Help, "complex") == 0)
    {
      printf("complex\n");
      printf("  - stores a complex number\n");
      printf("\n");
      printf("\n");
      printf("\n");
    }
    if (strcmp(Help, "spectra") == 0)
    {
      printf("spectra\n");
      printf("  - stores a spectral function\n");
      printf("\n");
      printf("\n");
      printf("\n");
    }
    
    
    if (strcmp(Help, "") == 0)
    {
      printf("\n");
      printf("\n");
      printf("\n");
    }
  }
  return 0;
}


static int LuaVerbosity(lua_State * L)
{
  GlobalVerbosity=(uint32_t)luaL_checkunsigned(L, -1);
  return 0;
}



static int LuaNewOperator(lua_State * L)
{
  if(lua_gettop(L)==0)
  {
    return luaL_error(L, "ERROR: LuaNewOperator called with 0 arguments, expect at least one\n");
  }
  if(lua_type(L, 1)!=LUA_TSTRING)
  {
    return luaL_error(L, "ERROR: first argument in LuaNewOperator should be a string (use \"Oppname\" as first argument)\n");
  }
  char const *OT =lua_tostring(L, 1);
  if ((strcmp(OT, "Coulomb") == 0)||(strcmp(OT, "U") == 0))
  {
    lua_remove(L, 1);
    return LuaOperatorNewCoulomb(L);
  }
  if (strcmp(OT, "Sx") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewSx(L);
  }
  if (strcmp(OT, "Sy") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewSy(L);
  }
  if (strcmp(OT, "Sz") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewSz(L);
  }
  if (strcmp(OT, "Splus") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewSplus(L);
  }
  if (strcmp(OT, "Smin") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewSmin(L);
  }
  if ((strcmp(OT, "Ssqr") == 0)||(strcmp(OT, "S^2") == 0)||(strcmp(OT, "S2") == 0))
  {
    lua_remove(L, 1);
    return LuaOperatorNewSsqr(L);
  }
  if (strcmp(OT, "Lx") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewLx(L);
  }
  if (strcmp(OT, "Ly") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewLy(L);
  }
  if (strcmp(OT, "Lz") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewLz(L);
  }
  if (strcmp(OT, "Lplus") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewLplus(L);
  }
  if (strcmp(OT, "Lmin") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewLmin(L);
  }
  if ((strcmp(OT, "Lsqr") == 0)||(strcmp(OT, "L^2") == 0)||(strcmp(OT, "L2") == 0))
  {
    lua_remove(L, 1);
    return LuaOperatorNewLsqr(L);
  }
  if (strcmp(OT, "Jx") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewJx(L);
  }
  if (strcmp(OT, "Jy") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewJy(L);
  }
  if (strcmp(OT, "Jz") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewJz(L);
  }
  if (strcmp(OT, "Jplus") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewJplus(L);
  }
  if (strcmp(OT, "Jmin") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewJmin(L);
  }
  if ((strcmp(OT, "Jsqr") == 0)||(strcmp(OT, "J^2") == 0)||(strcmp(OT, "J2") == 0))
  {
    lua_remove(L, 1);
    return LuaOperatorNewJsqr(L);
  }
  if (strcmp(OT, "ldots") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewldots(L);
  }
  if (strcmp(OT, "CF") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewCF(L);
  }
  if ((strcmp(OT, "Create") == 0)||(strcmp(OT, "Cr") == 0))
  {
    lua_remove(L, 1);
    return LuaOperatorNewCreate(L);
  }
  if ((strcmp(OT, "Annihilate") == 0)||(strcmp(OT, "An") == 0))
  {
    lua_remove(L, 1);
    return LuaOperatorNewAnnihilate(L);
  }
  if ((strcmp(OT, "Number") == 0)||(strcmp(OT, "CrAn") == 0)||(strcmp(OT, "n") == 0))
  {
    lua_remove(L, 1);
    return LuaOperatorNewNumber(L);
  }
  return luaL_error(L, "ERROR: first argument in ssp.NewOperator (%s) is not recognized as a known operatortype\n See ssp.Help(\"ssp.NewOperator\") for more information\n",OT);
}

/*  on stack: string describing the symmetry
 *  integer describing the angular momentum
 *  table with eigen energies of type double
 *
 *  output table of Harmonics Expansion
 */
static int LuaPotentialExpandedOnYlm(lua_State *L)
{
  char *s;
  double * eigenvalues;
  HarmonicsExpansion Akm;
  uint32_t l,len;
  s=(char*)luaL_checkstring(L, 1);
  l=luaL_checkint(L,2);
  eigenvalues=luaL_checkarray_double(L, 3, &len);
  if(PotentialExpandedOnY(s, l, eigenvalues, &Akm))
  {
    luaL_error(L, "Potential of \"%s\" symmetry for angular momentum %d not found\n",s,l);
  }
  LuaPushTable_HarmonicsExpansion(L, &Akm);
  FreeHarmonicsExpansion(&Akm);free(eigenvalues);
  return 1;
}


static int LuaNewWaveFunction(lua_State * L)
{
  char *detchar;
  Complex16 temp;
  uint32_t QComplex=0,NDet,i,detleng;
  uint8_t *det;
	WaveFunctionType *psi;
  uint32_t len = (uint32_t)lua_gettop(L);
	psi = (WaveFunctionType*) lua_newuserdata(L, sizeof(WaveFunctionType));
  luaL_setmetatable(L, "ssp.wavefunction");
  psi->NFermion=(uint32_t)luaL_checkunsigned(L, 1);
  psi->NBoson  =(uint32_t)luaL_checkunsigned(L, 2);
  psi->QComplex=1;
  if(len==3)
  {
    psi->NBitsKey=16;
  }
  else if(len==4)
  {
    psi->NBitsKey  =(uint32_t)luaL_checkunsigned(L, 4);
  }
  else
  {
    return luaL_error(L, "Error in LuaWavefunctionNew, expected 3 or 4 arguments, not %d\n",len);
  }
  WaveFunctionInit(psi);
  luaL_checktype(L, 3, LUA_TTABLE);
  NDet = (uint32_t)lua_objlen(L, 3);
  for(i=1;i<=NDet;i++)
  {
    lua_rawgeti(L, 3, i);
    luaL_checktype(L, -1, LUA_TTABLE);
    if(lua_objlen(L, -1)!=2)
    {
      return luaL_error(L, "Error in LuaWavefunctionNew, determinant number %d is not of the form {Determinant,prefactor}\n",i);
    }
    else
    {
      lua_rawgeti(L, -1, 1);
      if(lua_isstring(L, -1))
      {
        detchar = (char *) luaL_checkstring(L, -1);
        det = StringToDet(psi->NFermion, psi->NBoson, detchar);
        detleng = psi->NBoson + Min(psi->NFermion, 1 + (psi->NFermion - 1) / 8);
      }
      else
      {
        det=luaL_checkarray_unit8_t(L, -1, &detleng);
      }
      if(detleng!=psi->NWord)
      {
        return luaL_error(L, "Error in LuaWavefunctionNew, determinant number %d should be a table of %d words\n",i,psi->NWord);
      }
      lua_pop(L, 1);
      lua_rawgeti(L, -1, 2);
      if(lua_isnumber(L, -1))
      {
        temp.Re=(double)luaL_checknumber(L, -1);
        temp.Im=0;
      }
      else
      {
        QComplex=1;
        temp=*(Complex16 *)luaL_checkudata(L, -1, "Lua.complex");
      }
      lua_pop(L, 1);
      CheckDetPhysical(det, psi->NFermion, psi->NBoson);
      ComplexWaveFunctionAddElement(psi, det, temp.Re, temp.Im);
      free(det);
    }
    lua_pop(L, 1);
  }
  if(QComplex==0) MakeRealWaveFunctionFromComplex(*psi, psi);
  return 1;
}

/*  copy or clone the last element on the lua stack and place it on the stack
 *  works recursively on tables
 *  allocates new memory for user types
 */
static int Luacopy(lua_State * L)
{
  uint32_t len = (uint32_t)lua_gettop(L);
  switch (lua_type(L,len))
  {
    case LUA_TNIL:
      lua_pushnil(L);
      break;
    case LUA_TBOOLEAN:
      lua_pushboolean(L,lua_toboolean(L,len));
      break;
    case LUA_TNUMBER:
      lua_pushnumber(L,lua_tonumber(L,len));
      break;
    case LUA_TSTRING:
      lua_pushlstring(L,lua_tostring(L,len),lua_strlen(L,len));
      break;
    case LUA_TLIGHTUSERDATA:
      lua_pushlightuserdata(L,(void*)lua_touserdata(L,len));
      break;
    case
      LUA_TUSERDATA:
      if (luaL_testudata(L, len, "Lua.complex") != NULL)
      {
        Complex16 *a = (Complex16 *) luaL_checkudata(L, len, "Lua.complex");
        Complex16 *res;
        res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
        res->Re = a->Re;
        res->Im = a->Im;
        luaL_setmetatable(L, "Lua.complex");
        return 1;
      }
      if (luaL_testudata(L, len, "ssp.wavefunction") != NULL)
      {
        WaveFunctionType *psi = (WaveFunctionType*) luaL_checkudata(L, len, "ssp.wavefunction");
        WaveFunctionType *psiNew;
        psiNew = (WaveFunctionType *) lua_newuserdata(L, sizeof(WaveFunctionType));
        WaveFunctionCopy(*psi, psiNew);
        luaL_setmetatable(L, "ssp.wavefunction");
        return 1;
      }
      if (luaL_testudata(L, len, "ssp.operator") != NULL)
      {
        OperatorType *O = (OperatorType*) luaL_checkudata(L, len, "ssp.operator");
        OperatorType *ONew;
        ONew = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
        OperatorCopy(ONew, *O);
        luaL_setmetatable(L, "ssp.operator");
        return 1;
      }
      return luaL_error(L, "Error: copying or cloning an object of unrecognised type\n");
    case LUA_TTABLE:
      lua_newtable(L);
      lua_pushnil(L);
      while(lua_next(L, len) != 0)
      {
        Luacopy(L);
        lua_pushvalue(L, -3);
        lua_insert(L, -2);
        lua_settable(L, -5);
        lua_pop(L, 1);
      }
      break;
    default:
      return luaL_error(L, "Error: copying or cloning an object of unrecognised type\n");
  }
  return 1;
}

/*  copy or clone the last element on the lua stack and place it on the stack
 *  works recursively on tables
 *  allocates new memory for user types
 *  complex conjugates user types that can be conjugated
 */
static int Luaconjugate(lua_State * L)
{
  uint32_t len = (uint32_t)lua_gettop(L);
  switch (lua_type(L,len))
  {
    case LUA_TNIL:
      lua_pushnil(L);
      break;
    case LUA_TBOOLEAN:
      lua_pushboolean(L,lua_toboolean(L,len));
      break;
    case LUA_TNUMBER:
      lua_pushnumber(L,lua_tonumber(L,len));
      break;
    case LUA_TSTRING:
      lua_pushlstring(L,lua_tostring(L,len),lua_strlen(L,len));
      break;
    case LUA_TLIGHTUSERDATA:
      lua_pushlightuserdata(L,(void*)lua_touserdata(L,len));
      break;
    case
    LUA_TUSERDATA:
      if (luaL_testudata(L, len, "Lua.complex") != NULL)
      {
        Complex16 *a = (Complex16 *) luaL_checkudata(L, len, "Lua.complex");
        Complex16 *res;
        res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
        res->Re = a->Re;
        res->Im = -a->Im;
        luaL_setmetatable(L, "Lua.complex");
        return 1;
      }
      if (luaL_testudata(L, len, "ssp.wavefunction") != NULL)
      {
        WaveFunctionType *psi = (WaveFunctionType*) luaL_checkudata(L, len, "ssp.wavefunction");
        WaveFunctionType *psiNew;
        psiNew = (WaveFunctionType *) lua_newuserdata(L, sizeof(WaveFunctionType));
        WaveFunctionCopy(*psi, psiNew);
        WaveFunctionConjugate(*psiNew);
        luaL_setmetatable(L, "ssp.wavefunction");
        return 1;
      }
      if (luaL_testudata(L, len, "ssp.operator") != NULL)
      {
        OperatorType *O = (OperatorType*) luaL_checkudata(L, len, "ssp.operator");
        OperatorType *ONew;
        ONew = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
        OperatorCopy(ONew, *O);
        OperatorComplexConjugate(ONew);
        luaL_setmetatable(L, "ssp.operator");
        return 1;
      }
      return luaL_error(L, "Error: copying or cloning an object of unrecognised type\n");
    case LUA_TTABLE:
      lua_newtable(L);
      lua_pushnil(L);
      while(lua_next(L, len) != 0)
      {
        Luacopy(L);
        lua_pushvalue(L, -3);
        lua_insert(L, -2);
        lua_settable(L, -5);
        lua_pop(L, 1);
      }
      break;
    default:
      return luaL_error(L, "Error: copying or cloning an object of unrecognised type\n");
  }
  return 1;
}

static int Luaconjugatetranspose(lua_State * L)
{
  uint32_t len = (uint32_t)lua_gettop(L);
  switch (lua_type(L,len))
  {
    case LUA_TNIL:
      lua_pushnil(L);
      break;
    case LUA_TBOOLEAN:
      lua_pushboolean(L,lua_toboolean(L,len));
      break;
    case LUA_TNUMBER:
      lua_pushnumber(L,lua_tonumber(L,len));
      break;
    case LUA_TSTRING:
      lua_pushlstring(L,lua_tostring(L,len),lua_strlen(L,len));
      break;
    case LUA_TLIGHTUSERDATA:
      lua_pushlightuserdata(L,(void*)lua_touserdata(L,len));
      break;
    case
    LUA_TUSERDATA:
      if (luaL_testudata(L, len, "Lua.complex") != NULL)
      {
        Complex16 *a = (Complex16 *) luaL_checkudata(L, len, "Lua.complex");
        Complex16 *res;
        res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
        res->Re = a->Re;
        res->Im = -a->Im;
        luaL_setmetatable(L, "Lua.complex");
        return 1;
      }
      if (luaL_testudata(L, len, "ssp.wavefunction") != NULL)
      {
        WaveFunctionType *psi = (WaveFunctionType*) luaL_checkudata(L, len, "ssp.wavefunction");
        WaveFunctionType *psiNew;
        psiNew = (WaveFunctionType *) lua_newuserdata(L, sizeof(WaveFunctionType));
        WaveFunctionCopy(*psi, psiNew);
        WaveFunctionConjugate(*psiNew);
        luaL_setmetatable(L, "ssp.wavefunction");
        return 1;
      }
      if (luaL_testudata(L, len, "ssp.operator") != NULL)
      {
        OperatorType *O = (OperatorType*) luaL_checkudata(L, len, "ssp.operator");
        OperatorType *ONew;
        ONew = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
        OperatorConjugateTranspose(*O, ONew);
        luaL_setmetatable(L, "ssp.operator");
        return 1;
      }
      return luaL_error(L, "Error: copying or cloning an object of unrecognised type\n");
    case LUA_TTABLE:
      lua_newtable(L);
      lua_pushnil(L);
      while(lua_next(L, len) != 0)
      {
        Luacopy(L);
        lua_pushvalue(L, -3);
        lua_insert(L, -2);
        lua_settable(L, -5);
        lua_pop(L, 1);
      }
      break;
    default:
      return luaL_error(L, "Error: copying or cloning an object of unrecognised type\n");
  }
  return 1;
}

static int Luatranspose(lua_State * L)
{
  uint32_t len = (uint32_t)lua_gettop(L);
  switch (lua_type(L,len))
  {
    case LUA_TNIL:
      lua_pushnil(L);
      break;
    case LUA_TBOOLEAN:
      lua_pushboolean(L,lua_toboolean(L,len));
      break;
    case LUA_TNUMBER:
      lua_pushnumber(L,lua_tonumber(L,len));
      break;
    case LUA_TSTRING:
      lua_pushlstring(L,lua_tostring(L,len),lua_strlen(L,len));
      break;
    case LUA_TLIGHTUSERDATA:
      lua_pushlightuserdata(L,(void*)lua_touserdata(L,len));
      break;
    case
    LUA_TUSERDATA:
      if (luaL_testudata(L, len, "Lua.complex") != NULL)
      {
        Complex16 *a = (Complex16 *) luaL_checkudata(L, len, "Lua.complex");
        Complex16 *res;
        res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
        res->Re = a->Re;
        res->Im = a->Im;
        luaL_setmetatable(L, "Lua.complex");
        return 1;
      }
      if (luaL_testudata(L, len, "ssp.wavefunction") != NULL)
      {
        WaveFunctionType *psi = (WaveFunctionType*) luaL_checkudata(L, len, "ssp.wavefunction");
        WaveFunctionType *psiNew;
        psiNew = (WaveFunctionType *) lua_newuserdata(L, sizeof(WaveFunctionType));
        WaveFunctionCopy(*psi, psiNew);
        luaL_setmetatable(L, "ssp.wavefunction");
        return 1;
      }
      if (luaL_testudata(L, len, "ssp.operator") != NULL)
      {
        OperatorType *O = (OperatorType*) luaL_checkudata(L, len, "ssp.operator");
        OperatorType *ONew;
        ONew = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
        OperatorTranspose(*O, ONew);
        luaL_setmetatable(L, "ssp.operator");
        return 1;
      }
      return luaL_error(L, "Error: copying or cloning an object of unrecognised type\n");
    case LUA_TTABLE:
      lua_newtable(L);
      lua_pushnil(L);
      while(lua_next(L, len) != 0)
      {
        Luacopy(L);
        lua_pushvalue(L, -3);
        lua_insert(L, -2);
        lua_settable(L, -5);
        lua_pop(L, 1);
      }
      break;
    default:
      return luaL_error(L, "Error: copying or cloning an object of unrecognised type\n");
      break;
  }
  return 1;
}


static int Luarotate(lua_State *L)
{
  struct
  {
    double zero;
  } Options;
  Options.zero=100*DBL_EPSILON;
  MatrixType *U;
  uint32_t len = (uint32_t)lua_gettop(L);
  if(len<2)
  {
    return luaL_error(L, "ssp.Rotate called with %d arguments, 2 expected\n",len);
  }
  if (luaL_testudata(L, len-1, "ssp.wavefunction") != NULL)
  {
    return luaL_error(L, "ssp.Rotate called with a wave function, not yet implemented\n");
    
    return 0;
  }
  if (luaL_testudata(L, len-1, "ssp.operator") != NULL)
  {
    U = luaL_checkMatrix(L, len);
    OperatorType *O = (OperatorType*) luaL_checkudata(L, len-1, "ssp.operator");
    if(U->Ni!=U->Nj) luaL_error(L, "Error while rotating operator. Rotation matrix should be square (Ni=%d Nj=%d)",U->Ni,U->Nj);
    if(U->Ni!=O->NFermion+O->NBoson) luaL_error(L, "Error while rotating operator. Rotation matrix should have dimension (%d) of number of fermions (%d) + bosons (%d) of operator",U->Ni,O->NFermion,O->NBoson);
    OperatorType *ONew;
    ONew = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    RotateOperatorNormalOrder(ONew, *O, *U,&QDetNotOrderedFcFaBaBc);
    luaL_setmetatable(L, "ssp.operator");
    OperatorChop(ONew, Options.zero);
    return 1;
  }
  if(lua_type(L,len-1)==LUA_TTABLE)
  {
    return luaL_error(L, "ssp.Rotate called with a table, not yet implemented\n");
    
    return 0;
  }
  return luaL_error(L, "ssp.Rotate called unrecognized type. Can only rotate wavefunctions operators and tables of these\n");
}



static int LuaCreateSpectra(lua_State *L)
{
  void *buff;
  struct
  {
    uint32_t NTri;
    double epsilon;
    Bool restrictions;
    DeterminantRestrictionType DeterminantRestrictions;
    double Emin;
    double Emax;
    uint32_t NE;
    double Gamma;
    uint32_t Npsi;
    uint32_t NT;
  } Options;
  Options.NTri=200;
  Options.epsilon=10000*DBL_EPSILON;
  Options.restrictions=False;
  Options.Emin=DBL_MAX;
  Options.Emax=-DBL_MAX;
  Options.NE=1000;
  Options.Gamma=0;
  SpectraTypeArray *Spec;
  TriDiagonalMatrixType *GTri;
  OperatorType *H,*T;
  WaveFunctionType *psi,Tpsi;
  uint32_t i,iT,ipsi,len;
  double Eminpole,Emaxpole, *val, *psiEn, psiEntemp;
  
  /* check number of input parameters */
  len = (uint32_t)lua_gettop(L);
  if(!((len==3)||(len==4)))
  {
    return luaL_error(L, "Error CreateSpectra called with %d arguments, 3 or 4 expected\n",len);
  }
  
  /* read operator and wavefunctions */
  H = (OperatorType *) luaL_checkudata(L, 1, "ssp.operator");
  if(lua_type(L,2)==LUA_TTABLE)
  {
    T=luaL_checkarray_OperatorType(L, 2, &Options.NT);
  }
  else
  {
    Options.NT = 1;
    T = (OperatorType *) luaL_checkudata(L, 2, "ssp.operator");
  }
  if(lua_type(L,3)==LUA_TTABLE)
  {
    psi=luaL_checkarray_WaveFunctionType(L, 3, &Options.Npsi);
  }
  else
  {
    Options.Npsi = 1;
    psi = (WaveFunctionType *) luaL_checkudata(L, 3, "ssp.wavefunction");
  }
  /* read options */
  
  buff=LuaFindOptionsInTable(L, -1, "epsilon", "double");
  if(buff!=NULL)
  {
    Options.epsilon = *(double*)buff;
    free(buff);
  }
  buff=LuaFindOptionsInTable(L, -1, "restrictions", "DeterminantRestrictionType");
  if(buff!=NULL)
  {
    Options.restrictions=True;
    Options.DeterminantRestrictions = *(DeterminantRestrictionType*)buff;
    free(buff);
  }
  buff=LuaFindOptionsInTable(L, -1, "Emin", "double");
  if(buff!=NULL)
  {
    Options.Emin = *(double*)buff;
    free(buff);
  }
  buff=LuaFindOptionsInTable(L, -1, "Emax", "double");
  if(buff!=NULL)
  {
    Options.Emax = *(double*)buff;
    free(buff);
  }
  buff=LuaFindOptionsInTable(L, -1, "NE", "int");
  if(buff!=NULL)
  {
    Options.NE = *(uint32_t*)buff;
    free(buff);
  }
  buff=LuaFindOptionsInTable(L, -1, "Gamma", "double");
  if(buff!=NULL)
  {
    Options.Gamma = *(double*)buff;
    free(buff);
  }
  /* calculate spectra in Tridiagonal form */
  GTri = malloc(Options.NT*Options.Npsi*sizeof(TriDiagonalMatrixType));
  for(i=0;i<Options.NT*Options.Npsi;i++)
  {
    GTri[i].N=Options.NTri;
    InitTriDiagonalMatrix(&GTri[i]);
  }
  /* calculate groundstate energy */
  psiEn=malloc(Options.Npsi*Options.NT*sizeof(double));
  for(ipsi=0;ipsi<Options.Npsi;ipsi++)
  {
    PsiHermitianOperatorPsi(*H, psi[ipsi], &psiEntemp);
    for(iT=0;iT<Options.NT;iT++)
    {
      psiEn[iT*Options.Npsi+ipsi]=psiEntemp;
    }
  }
  /* loop over all Transition operators */
  for(iT=0;iT<Options.NT;iT++)
  {
    /* loop over all wavefunctions */
    for(ipsi=0;ipsi<Options.Npsi;ipsi++)
    {
      if(Options.restrictions)
      {
        OperatorPsiRestricted(T[iT], psi[ipsi], &Tpsi, Options.DeterminantRestrictions);
        LanczosTriDiagonalizeRestricted(*H, Tpsi, Options.DeterminantRestrictions, Options.epsilon, &GTri[iT*Options.Npsi+ipsi], Max(Options.Gamma,Options.epsilon));
      }
      else
      {
        OperatorPsi(T[iT], psi[ipsi], &Tpsi);
        LanczosTriDiagonalize(*H, Tpsi, Options.epsilon, &GTri[iT*Options.Npsi+ipsi], Max(Options.Gamma,Options.epsilon));
      }
    }
  }
  /* determine Eminpole and Emaxpole */
  /* get eigenvalues of all Tri diagonal matrices and set smallest and largest eigenvalue */
  Eminpole=DBL_MAX;
  Emaxpole=-DBL_MAX;
  val=malloc((Options.NTri)*sizeof(double));
  for(i=0;i<Options.NT*Options.Npsi;i++)
  {
    TriDiagonalMatrixToEigenvalues(GTri[i], val);
    if((val[0]-psiEn[i])<Eminpole)Eminpole=val[0]-psiEn[i];
    if((val[Options.NTri-1]-psiEn[i])>Emaxpole)Emaxpole=val[Options.NTri-1]-psiEn[i];
  }
  
  
  /* if Emin=DBL_MAX or Emax=-DBL_MAX or Gamma==0 determine these */

  
  if(Options.Emin== DBL_MAX) Options.Emin = Eminpole - 0.125 * (Emaxpole-Eminpole);
  if(Options.Emax==-DBL_MAX) Options.Emax = Emaxpole + 0.125 * (Emaxpole-Eminpole);
  
  if(Options.Gamma==0) Options.Gamma = 10*(Options.Emax-Options.Emin)/((double)Options.NE);
  /* transform tridiagonal form to spectral form */
  Spec = (SpectraTypeArray *) lua_newuserdata(L, sizeof(SpectraTypeArray));
  luaL_setmetatable(L, "ssp.spectra");
  Spec->NSpectra=Options.NT*Options.Npsi;
  Spec->Spec=malloc(Spec->NSpectra*sizeof(SpectraType));
  Spec->Spec[0].N=Options.NE;
  Spec->Spec[0].Emin=Options.Emin;
  Spec->Spec[0].Emax=Options.Emax;
  Spec->Spec[0].dE=(Options.Emax-Options.Emin)/((double)Options.NE);
  Spec->Spec[0].EminPole=Eminpole;
  Spec->Spec[0].EmaxPole=Emaxpole;
  Spec->Spec[0].Gamma=Options.Gamma;
  strcpy(Spec->Spec[0].Name, "Spec");
  for(i=0;i<Spec->NSpectra;i++)
  {
    SpectrumAllocate(Spec->Spec[0], &(Spec->Spec[i]));
    sprintf(Spec->Spec[i].Name, "Spectrum %i",i);
    TriDiagonalMatrixToSpectrum(psiEn[i], GTri[i], &(Spec->Spec[i]));
  }
  for(i=0;i<Spec->NSpectra;i++)
  {
    FreeTriDiagonalMatrix(&GTri[i]);
  }
  free(GTri);
  free(psiEn);
  return 1;
}



/*  calculates <psi_i | T_1^+ G_1 T_2 G_2 T_2^+ G_1 T_1 | psi_i >
 *  On stack:
 *     H1
 *     H2
 *     T1
 *     T2^+
 *     psi
 *     {E1min, E1max, N1, Gamma1}
 *     {E2min, E2max, N2, Gamma2}
 *     options
 *        NTri1
 *        NTri2
 *        epsilon
 *        restrictions1
 */
static int LuaCreateRIXS(lua_State *L)
{
  void *buff;
  struct
  {
    uint32_t NTri1;
    uint32_t NTri2;
    double epsilon;
    Bool restrictions;
    DeterminantRestrictionType DeterminantRestrictions;
    double Emin1;
    double Emax1;
    uint32_t NE1;
    double Gamma1;
    double Emin2;
    double Emax2;
    uint32_t NE2;
    double Gamma2;
    uint32_t Npsi;
    uint32_t NT1;
    uint32_t NT2;
  } Options;
  Options.NTri1=200;
  Options.NTri2=200;
  Options.epsilon=10000*DBL_EPSILON;
  Options.restrictions=False;
  SpectraTypeArray *Spec;
  OperatorType *H,*T1, *T2;
  WaveFunctionType *psi;
  uint32_t iT,ipsi,len;
  double *psiEn;
  
  /* check number of input parameters */
  len = (uint32_t)lua_gettop(L);
  if(!((len==4)||(len==5)))
  {
    return luaL_error(L, "Error CreateFluorescenceYield called with %d arguments, 4 or 5 expected\n",len);
  }
  
  /* read operator and wavefunctions */
  H = (OperatorType *) luaL_checkudata(L, 1, "ssp.operator");
  if(lua_type(L,2)==LUA_TTABLE)
  {
    T1=luaL_checkarray_OperatorType(L, 2, &Options.NT1);
  }
  else
  {
    Options.NT1 = 1;
    T1 = (OperatorType *) luaL_checkudata(L, 2, "ssp.operator");
  }
  if(lua_type(L,3)==LUA_TTABLE)
  {
    T2=luaL_checkarray_OperatorType(L, 3, &Options.NT2);
  }
  else
  {
    Options.NT2 = 1;
    T2 = (OperatorType *) luaL_checkudata(L, 3, "ssp.operator");
  }
  if(lua_type(L,4)==LUA_TTABLE)
  {
    psi=luaL_checkarray_WaveFunctionType(L, 4, &Options.Npsi);
  }
  else
  {
    Options.Npsi = 1;
    psi = (WaveFunctionType *) luaL_checkudata(L, 4, "ssp.wavefunction");
  }
  /* read options */
  buff=LuaFindOptionsInTable(L, -1, "epsilon", "double");
  if(buff!=NULL)
  {
    Options.epsilon = *(double*)buff;
    free(buff);
  }
  buff=LuaFindOptionsInTable(L, -1, "restrictions", "DeterminantRestrictionType");
  if(buff!=NULL)
  {
    Options.restrictions=True;
    Options.DeterminantRestrictions = *(DeterminantRestrictionType*)buff;
    free(buff);
  }
  /* allocate spectra and and pointer to return matrix on lua stack */
  Spec = (SpectraTypeArray *) lua_newuserdata(L, sizeof(SpectraTypeArray));
  luaL_setmetatable(L, "ssp.spectra");
  Spec->NSpectra=Options.NT1*Options.NT2*Options.Npsi;
  Spec->Spec=malloc(Spec->NSpectra*sizeof(SpectraType));
  /* calculate groundstate energy */
  psiEn=malloc(Options.Npsi*sizeof(double));
  for(ipsi=0;ipsi<Options.Npsi;ipsi++)
  {
    PsiHermitianOperatorPsi(*H, psi[ipsi], &psiEn[ipsi]);
  }
  /* loop over all Transition operators */

  free(psiEn);
  return 1;
}





static int LuaCreateFluorescenceYield(lua_State *L)
{
  void *buff;
  struct
  {
    uint32_t NTri;
    double epsilon;
    Bool restrictions;
    DeterminantRestrictionType DeterminantRestrictions;
    double Emin;
    double Emax;
    uint32_t NE;
    double Gamma;
    uint32_t Npsi;
    uint32_t NT1;
    uint32_t NT2;
  } Options;
  Options.NTri=200;
  Options.epsilon=10000*DBL_EPSILON;
  Options.restrictions=False;
  Options.Emin=DBL_MAX;
  Options.Emax=-DBL_MAX;
  Options.NE=1000;
  Options.Gamma=0;
  SpectraTypeArray *Spec;
  OperatorType *H,*T1, *T2;
  WaveFunctionType *psi;
  uint32_t iT,ipsi,len;
  double *psiEn;
  
  /* check number of input parameters */
  len = (uint32_t)lua_gettop(L);
  if(!((len==4)||(len==5)))
  {
    return luaL_error(L, "Error CreateFluorescenceYield called with %d arguments, 4 or 5 expected\n",len);
  }
  
  /* read operator and wavefunctions */
  H = (OperatorType *) luaL_checkudata(L, 1, "ssp.operator");
  if(lua_type(L,2)==LUA_TTABLE)
  {
    T1=luaL_checkarray_OperatorType(L, 2, &Options.NT1);
  }
  else
  {
    Options.NT1 = 1;
    T1 = (OperatorType *) luaL_checkudata(L, 2, "ssp.operator");
  }
  if(lua_type(L,3)==LUA_TTABLE)
  {
    T2=luaL_checkarray_OperatorType(L, 3, &Options.NT2);
  }
  else
  {
    Options.NT2 = 1;
    T2 = (OperatorType *) luaL_checkudata(L, 3, "ssp.operator");
  }
  if(lua_type(L,4)==LUA_TTABLE)
  {
    psi=luaL_checkarray_WaveFunctionType(L, 4, &Options.Npsi);
  }
  else
  {
    Options.Npsi = 1;
    psi = (WaveFunctionType *) luaL_checkudata(L, 4, "ssp.wavefunction");
  }
  /* read options */
  buff=LuaFindOptionsInTable(L, -1, "epsilon", "double");
  if(buff!=NULL)
  {
    Options.epsilon = *(double*)buff;
    free(buff);
  }
  buff=LuaFindOptionsInTable(L, -1, "restrictions", "DeterminantRestrictionType");
  if(buff!=NULL)
  {
    Options.restrictions=True;
    Options.DeterminantRestrictions = *(DeterminantRestrictionType*)buff;
    free(buff);
  }
  buff=LuaFindOptionsInTable(L, -1, "Emin", "double");
  if(buff!=NULL)
  {
    Options.Emin = *(double*)buff;
    free(buff);
  }
  buff=LuaFindOptionsInTable(L, -1, "Emax", "double");
  if(buff!=NULL)
  {
    Options.Emax = *(double*)buff;
    free(buff);
  }
  buff=LuaFindOptionsInTable(L, -1, "NE", "int");
  if(buff!=NULL)
  {
    Options.NE = *(uint32_t*)buff;
    free(buff);
  }
  buff=LuaFindOptionsInTable(L, -1, "Gamma", "double");
  if(buff!=NULL)
  {
    Options.Gamma = *(double*)buff;
    free(buff);
  }
  /* allocate spectra and and pointer to return matrix on lua stack */
  Spec = (SpectraTypeArray *) lua_newuserdata(L, sizeof(SpectraTypeArray));
  luaL_setmetatable(L, "ssp.spectra");
  Spec->NSpectra=Options.NT1*Options.NT2*Options.Npsi;
  Spec->Spec=malloc(Spec->NSpectra*sizeof(SpectraType));
  /* calculate groundstate energy */
  psiEn=malloc(Options.Npsi*sizeof(double));
  for(ipsi=0;ipsi<Options.Npsi;ipsi++)
  {
    PsiHermitianOperatorPsi(*H, psi[ipsi], &psiEn[ipsi]);
  }
  /* loop over all Transition operators */
  for(iT=0;iT<Options.NT1;iT++)
  {
    /* loop over all wavefunctions */
    for(ipsi=0;ipsi<Options.Npsi;ipsi++)
    {
      Spec->Spec[Options.NT2*(iT*Options.Npsi+ipsi)].N=Options.NE;
      if(ipsi==0)
      {
        Spec->Spec[Options.NT2*(iT*Options.Npsi+ipsi)].Emin=Options.Emin;
        Spec->Spec[Options.NT2*(iT*Options.Npsi+ipsi)].Emax=Options.Emax;
      }
      else
      {
        Spec->Spec[Options.NT2*(iT*Options.Npsi+ipsi)].Emin=Spec->Spec[Options.NT2*(iT*Options.Npsi)].Emin;
        Spec->Spec[Options.NT2*(iT*Options.Npsi+ipsi)].Emax=Spec->Spec[Options.NT2*(iT*Options.Npsi)].Emax;
      }
      Spec->Spec[Options.NT2*(iT*Options.Npsi+ipsi)].Gamma=Options.Gamma;
      if(Options.restrictions)
      {
        FluorescenceYieldRestricted(*H, Options.DeterminantRestrictions, T1[iT], T2, Options.NT2, psi[ipsi], psiEn[ipsi], Options.epsilon, Options.NTri, &(Spec->Spec[Options.NT2*(iT*Options.Npsi+ipsi)]));
      }
      else
      {
        FluorescenceYield(*H, T1[iT], T2, Options.NT2, psi[ipsi], psiEn[ipsi], Options.epsilon, Options.NTri, &(Spec->Spec[Options.NT2*(iT*Options.Npsi+ipsi)]));
      }
    }
  }
  free(psiEn);
  return 1;
}




/* given one operator and a (table of) wavefunctions calculate the ground-state (and excited states) starting from those wavefunctions.
 *
 */
static int LuaEvolveToGroundstate(lua_State * L)
{
  void *buff;
  struct
  {
    uint32_t Npsi;
    double zero;
    double epsilon;
    Bool restrictions;
    DeterminantRestrictionType DeterminantRestrictions;
  } Options;
  Options.zero=100*DBL_EPSILON;
  Options.epsilon=10000*DBL_EPSILON;
  Options.restrictions=False;
  OperatorType *Opp;
  WaveFunctionType *psi,*temp;
  uint32_t i,len;
  /* check number of input parameters */
  len = (uint32_t)lua_gettop(L);
  if(!((len==2)||(len==3)))
  {
    return luaL_error(L, "Error EvolveToGroundstate called with %d arguments, 2 or 3 expected\n",len);
  }
  /* read operator and wavefunctions */
  Opp = (OperatorType *) luaL_checkudata(L, 1, "ssp.operator");
  if(lua_type(L,2)==LUA_TTABLE)
  {
    psi=luaL_checkarray_WaveFunctionType(L, 2, &Options.Npsi);
  }
  else
  {
    Options.Npsi = 1;
    psi = (WaveFunctionType *) luaL_checkudata(L, 2, "ssp.wavefunction");
  }
  /* read options */
  buff=LuaFindOptionsInTable(L, -1, "epsilon", "double");
  if(buff!=NULL)
  {
    Options.epsilon = *(double*)buff;
    free(buff);
  }
  buff=LuaFindOptionsInTable(L, -1, "zero", "double");
  if(buff!=NULL)
  {
    Options.zero = *(double*)buff;
    free(buff);
  }
  buff=LuaFindOptionsInTable(L, -1, "restrictions", "DeterminantRestrictionType");
  if(buff!=NULL)
  {
    Options.restrictions=True;
    Options.DeterminantRestrictions = *(DeterminantRestrictionType*)buff;
    free(buff);
  }
  /* calculate groundstate */
  for (i = 0; i < Options.Npsi; i++)
  {
    printf("\nstart of state #%d/%d",i+1,Options.Npsi);
    if(Options.restrictions)
    {
      LanczosGroundStateRestricted(*Opp, &(psi[i]), Options.DeterminantRestrictions, &Options.epsilon, Options.zero, 0, psi, i);
    }
    else
    {
      LanczosGroundState(*Opp, &(psi[i]), &Options.epsilon, Options.zero, 0, psi, i);
    }
  }
  /* copy changed wavefunctions back to lua stack */
  if(lua_type(L,2)==LUA_TTABLE)
  {
    for(i=0;i<Options.Npsi;i++)
    {
      lua_rawgeti(L, 2, i+1);
      temp = (WaveFunctionType *)luaL_checkudata(L, -1, "ssp.wavefunction");
      *temp = psi[i];
      lua_pop(L, 1);
    }
  }
  return 0;
}

/*  calculate the Npsi lowest eigen-states of operator Opp, using the Block routines
 *
 */
static int LuaEigensystem(lua_State * L)
{
  void *buff;
  struct
  {
    ConvergenceSettingsType Settings;
    Bool restrictions;
    uint32_t NBitsKey;
    DeterminantRestrictionType DeterminantRestrictions;
  } Options;
  Options.Settings=InitConvergenceSettings();
  Options.restrictions=False;
  Options.DeterminantRestrictions.N=0;
  Options.NBitsKey=16;
  OperatorType *Opp;
  DeterminantRestrictionType StartRestrictions;
  StartRestrictions.N=0;
  WaveFunctionType *psi=NULL,*temp;
  uint32_t i,len, Npsi=0,QStartRestrictions=0;
  /* check number of input parameters */
  len = (uint32_t)lua_gettop(L);
  if(!((len==2)||(len==3)||(len==4)))
  {
    return luaL_error(L, "Error Eigensystem called with %d arguments, 2, 3 or 4 expected\n",len);
  }
  /* read operator, (StartRestrictions or StartPsi) and Number of eigen-states */
  Opp = (OperatorType *) luaL_checkudata(L, 1, "ssp.operator");
  if(lua_type(L,2)==LUA_TTABLE)
  {
    lua_rawgeti(L, 2, 1);
    if(lua_type(L,-1)==LUA_TNUMBER) QStartRestrictions=1;
    lua_pop(L, 1);
    if(QStartRestrictions)
    {
      buff=luaL_checkarray_DeterminantRestrictionType(L, 2);
      StartRestrictions = *(DeterminantRestrictionType*)buff;
      free(buff);
    }
    else
    {
      psi=luaL_checkarray_WaveFunctionType(L, 2, &Npsi);
      if(WaveFunctionArrayToSameBasis(psi, Npsi))
      {
        printf("WaveFunctionArrayToSameBasis failed in LuaEigensystem\n");
        return 1;
      }
    }
  }
  else
  {
    Npsi = 1;
    psi = (WaveFunctionType *) luaL_checkudata(L, 2, "ssp.wavefunction");
  }
  if((len==2)&&(QStartRestrictions))
  {
    return luaL_error(L, "Error Eigensystem called with 2 arguments, whereby the first is the operator and the second a set of restrictions\n Need a thrid argument how many eigenstates should be calculated\n");
  }
  if(QStartRestrictions)
  {
    Npsi=(uint32_t)luaL_checkunsigned(L, 3);
    // allocate psi
    if(Npsi==1)
    {
      psi = (WaveFunctionType*) lua_newuserdata(L, sizeof(WaveFunctionType));
      luaL_setmetatable(L, "ssp.wavefunction");
    }
    else
    {
      psi=malloc(Npsi*sizeof(WaveFunctionType));
    }
  }
  /* read options */
  buff=LuaFindOptionsInTable(L, -1, "DenseBoarder", "int");
  if(buff!=NULL)
  {
    Options.Settings.DenseBoarder = *(uint32_t*)buff;
    free(buff);
  }
  buff=LuaFindOptionsInTable(L, -1, "NKrylovStart", "int");
  if(buff!=NULL)
  {
    Options.Settings.NKrylovStart = *(uint32_t*)buff;
    free(buff);
  }
  buff=LuaFindOptionsInTable(L, -1, "NKrylovMax", "int");
  if(buff!=NULL)
  {
    Options.Settings.NKrylovMax = *(uint32_t*)buff;
    free(buff);
  }
  buff=LuaFindOptionsInTable(L, -1, "NKrylovStep", "int");
  if(buff!=NULL)
  {
    Options.Settings.NKrylovStep = *(uint32_t*)buff;
    free(buff);
  }
  buff=LuaFindOptionsInTable(L, -1, "NKrylovMin", "int");
  if(buff!=NULL)
  {
    Options.Settings.NKrylovMin = *(uint32_t*)buff;
    free(buff);
  }
  buff=LuaFindOptionsInTable(L, -1, "NBitsKey", "int");
  if(buff!=NULL)
  {
    Options.NBitsKey = *(uint32_t*)buff;
    free(buff);
  }
  buff=LuaFindOptionsInTable(L, -1, "epsilon", "double");
  if(buff!=NULL)
  {
    Options.Settings.epsilon = *(double*)buff;
    free(buff);
  }
  buff=LuaFindOptionsInTable(L, -1, "zero", "double");
  if(buff!=NULL)
  {
    Options.Settings.zero = *(double*)buff;
    free(buff);
  }
  buff=LuaFindOptionsInTable(L, -1, "restrictions", "DeterminantRestrictionType");
  if(buff!=NULL)
  {
    Options.restrictions=True;
    Options.DeterminantRestrictions = *(DeterminantRestrictionType*)buff;
    free(buff);
  }
  /* calculate groundstate */
  BlockGroundState(*Opp, Npsi, &(Options.Settings) , StartRestrictions, Options.DeterminantRestrictions, psi, Options.NBitsKey);
  /* seperate basis */
  WaveFunctionArraySeperateBasis(psi, Npsi);
  /* copy changed wavefunctions back to lua stack */
  if(QStartRestrictions) // need to return a list of wavefunctions
  {
    if(Npsi==1)
    {
      // return a single wave function
      // already done when allocating
    }
    else
    {
      LuaPushTable_WaveFunctionType(L,psi,Npsi);
      free(psi);
    }
    return 1;
  }
  else
  { // the list of functions is copied into an array, now need to set the pointers in the structure back to what they should be
    // this is not needed for a single value as there the pointer psi just points to the original
    if(lua_type(L,2)==LUA_TTABLE)
    {
      for(i=0;i<Npsi;i++)
      {
        lua_rawgeti(L, 2, i+1);
        temp = (WaveFunctionType *)luaL_checkudata(L, -1, "ssp.wavefunction");
        *temp = psi[i];
        lua_pop(L, 1);
      }
    }
    return 0;
  }
}



static int LuaTimeEnd(lua_State * L)
{
  char *string;
  string=(char*)luaL_checkstring(L, -1);
  EndTiming(string, &GlobalTime);
  return 0;
}


static int LuaTimeStart(lua_State * L)
{
  char *string;
  string=(char*)luaL_checkstring(L, -1);
  BeginTiming(string, &GlobalTime);
  return 0;
}


static int LuaTimePrint(lua_State * L)
{
  PrintTiming(GlobalTime);
  return 0;
}

static int LuaDensityMatrix(lua_State *L)
{
  uint32_t len, Nparticles, *jindex;
  WaveFunctionType *psi;
  MatrixType Rho;
  len = (uint32_t)lua_gettop(L);
  // first argument is a wave-function or a list of functions
  switch (lua_type(L,1))
  {
    case LUA_TTABLE:
    {
      lua_newtable(L);
      lua_pushnil(L);
      while(lua_next(L, 1) != 0)
      {
        lua_insert(L, 1);
        lua_pushvalue(L, -3);
        LuaDensityMatrix(L);
        lua_pushvalue(L, -3);
        lua_replace(L, -3);
        lua_settable(L, -4);
        lua_remove(L, 1);
      }
      return 1;
    }
    case LUA_TUSERDATA:
    {
      if (luaL_testudata(L, 1, "ssp.wavefunction") != NULL)
      {
        psi = (WaveFunctionType *) luaL_checkudata(L, 1, "ssp.wavefunction");
      }
      else
      {
        return luaL_error(L, "Error: DensityMatrix expects a function or a list of functions on the first position\n");
      }
      break;
    }
    default:
    {
      return luaL_error(L, "Error: DensityMatrix expects a function or a list of functions on the first position\n");
    }
  }
  // the argument at position len is a list of integers
  jindex = luaL_checkarray_unit32_t(L,len,&Nparticles);
  // calculate the density matrix
  CalculatePartialDensityMatrix(&Rho, *psi, Nparticles, jindex);
  // return a lua table of tables of dimension Nparticles by Nparticles
  LuaPushTable_MatrixType(L, &Rho);
  FreeMatrix(&Rho);
  return 1;
}


static const struct luaL_Reg ssp_functions[] = {
  {"Help", LuaHelp},
  {"Verbosity",LuaVerbosity},
  {"TimeStart",LuaTimeStart},
  {"TimeEnd",LuaTimeEnd},
  {"TimePrint",LuaTimePrint},
  {"NewOperator", LuaNewOperator},
  {"PotentialExpandedOnYlm", LuaPotentialExpandedOnYlm},
  {"NewWavefunction", LuaNewWaveFunction},
  {"DensityMatrix",LuaDensityMatrix},
  {"Copy", Luacopy},
  {"Clone", Luacopy},
  {"Conjugate", Luaconjugate},
  {"C", Luaconjugate},
  {"ConjugateTranspose", Luaconjugatetranspose},
  {"CT", Luaconjugatetranspose},
  {"Transpose", Luatranspose},
  {"T", Luatranspose},
  {"Rotate", Luarotate},
  {"EvolveToGroundstate", LuaEvolveToGroundstate},
  {"Eigensystem", LuaEigensystem},
  {"CreateSpectra", LuaCreateSpectra},
  {"CreateFluorescenceYield", LuaCreateFluorescenceYield},
  {"BraKet", LuaBraKet},           /* function in file BraKet.c */
  {"Eval", LuaEval},               /* function in file eval.c */
  {NULL, NULL}
};

void luaconstants(lua_State *L)
{
  Complex16 *luaI;
  luaI = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
  luaI->Re = 0.;
  luaI->Im = 1.;
  luaL_setmetatable(L, "Lua.complex");
  lua_setglobal(L, "I");
  lua_pushnumber(L, 3.14159265358979323846264338328);
  lua_setglobal(L, "Pi");
  lua_pushnumber(L, 3.14159265358979323846264338328);
  lua_setglobal(L, "pi");
  lua_pushnumber(L, 2.71828182845904523536028747135);
  lua_setglobal(L, "E");
  lua_pushnumber(L, 2.71828182845904523536028747135);
  lua_setglobal(L, "e");
}

int luaopen_ssp(lua_State * L)
{
  luaL_newmetatable(L, "ssp.operator");
  luaL_setfuncs(L, Operator_meta, 0);
  
  luaL_newmetatable(L, "ssp.wavefunction");
  luaL_setfuncs(L, Wavefunction_meta, 0);
  
  luaL_newmetatable(L, "ssp.spectra");
  luaL_setfuncs(L, Spectra_meta, 0);
    
  luaL_newlib(L, ssp_functions);
  lua_setglobal(L, "ssp");
  
  
  
  return 1;
}