/**
 * Prints the help file to screen
 * @param L
 *        Pointer to a lua state
 * @return 0 
 *         as it leaves the lua state unaltered
 */
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
    printf("Help()\n");
    printf("  - shows this file\n");
    printf("  - excepts as input one of the other functions in order to show help on that function\n");
    printf("  - example: Help(\"NewOperator\");\n");
    printf("\n");
    printf("NewOperator()\n");
    printf("  - creates a new operator\n");
    printf("  - expects as input on the first position a string describing the type of operator\n");
    printf("  - other parameters depend on the type of operator defined in the first string\n");
    printf("  - example: NewOperator(\"Sz\", NFermion, IndexUp, IndexDn);\n");
    printf("\n");
    printf("PotentialExpandedOnClm(\"Oh\",2,{0.6,-0.4})\n");
    printf("  - creates a list representing a potential expanded on spherical Harmonics for a given symmetry and eigen energy of the irreducible representations\n");
    printf("  - The potential V(theta, phi) is given as:\n");
    printf("  - V(theta, phi) = sum_k,m A_{k,m} Y_km(theta,phi)\n");
    printf("  - the list saved is {{k1,m1,Ak1m1},{k2,m2,Ak2m2},...}\n");
    printf("  - input is the point-group symmetry, the angular momentum and the eigen energies of the irreducible representations\n");
    printf("  - PotentialExpandedOnClm(\"Oh\",2,{0.6,-0.4});\n");
    printf("  - for a potential in cubic symmetry acting on the d-shell with the energy of the eg's as 0.6 and t2g's at -0.4 (10Dq=1.0)\n");
    printf("\n");
    printf("NewWavefunction()\n");
    printf("  - creates a new wavefunction\n");
    printf("  - expects as input the number of fermions, number of bosons, a table of determinants with the prefactor\n");
    printf("  - example: NewWavefunction(10,0,{{{0x80,0x00},1},{{0x04,0x00},1}});\n");
    printf("  - example: NewWavefunction(10,0,{{\"1010101010\",1},{\"0101010101\",1}});\n");
    printf("\n");
    printf("Copy() or Clone()\n");
    printf("  - creats a copy of a wavefunction or operator\n");
    printf("  - expects as input an object, returns a copy of that object\n");
    printf("  - example: Copy(psi);\n");
    printf("\n");
    printf("Conjugate()\n");
    printf("  - creats a copy of a wavefunction or operator and conjugates this copy\n");
    printf("  - expects as input an object, returns a copy of that object\n");
    printf("  - example: Conjugate(psi);\n");
    printf("\n");
    printf("ConjugateTranspose()\n");
    printf("  - creats a copy of a wavefunction or operator and conjugates transposes this copy\n");
    printf("  - expects as input an object, returns a copy of that object\n");
    printf("  - example: ConjugateTranspose(psi);\n");
    printf("\n");
    printf("Transpose()\n");
    printf("  - creats a copy of a wavefunction or operator and transposes this copy\n");
    printf("  - expects as input an object, returns a copy of that object\n");
    printf("  - example: Transpose(psi);\n");
    printf("Rotate()\n");
    printf("  - Given an operator on a basis of N fermions + M Bosons\n");
    printf("  - and a matrix of size N+M\n");
    printf("  - creates a new operator on a rotated basis\n");
    printf("  - example: Rotate(Operator,RotationMatrix);\n");
    printf("\n");
    printf("EvolveToGroundstate()\n");
    printf("  - calculates eigenfunctions of the Hamiltonian, starting from the specified state\n");
    printf("  - expects as input an operator and a (list of) wavefunction(s), a list of options on the last position is allowed\n");
    printf("  - example EvolveToGroundstate(OppSy,{psi1,psi2});\n");
    printf("Eigensystem()\n");
    printf("  - calculates the lowest N eigenfunctions of a given operator\n");
    printf("  - expects as input an operator and either a list of starting functions or a set of restrictions specifying the starting point\n");
    printf("  - example Eigensystem(OppSy,{psi1,psi2});\n");
    printf("  - example Eigensystem(OppSy,Restrictions,Nspi);\n");
    printf("\n");
    printf("CreateSpectra()\n");
    printf("  - calculates spectra of a giving Hamiltonian, starting from the specified state and using a specific transition operator\n");
    printf("  - expects as input an operator a (list of) transition operator(s) and a (list of) wavefunction(s), a list of options on the last position is allowed\n");
    printf("  - example CreateSpectra(OppH,{T1,T2},{psi1,psi2});\n");
    printf("\n");
    printf("CreateFluorescenceYield()\n");
    printf("  - calculates the fluorescence yield spectra\n");
    printf("  - calculated is < i | T1^+ G^+ T2 T2^+ G T1 | i>\n");
    printf("  - expects as input an operator a (list of) transition operator(s) T1 a (list of) transition operator(s) T2^+ and a (list of) wavefunction(s), a list of options on the last position is allowed\n");
    printf("  - example CreateFluorescenceYield(OppH,{T1a,T1b},{T2a,T2b},{psi1,psi2});\n");
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
    if (strcmp(Help, "NewOperator") == 0)
    {
      printf("NewOperator()\n");
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
      printf("  - OppSz = NewOperator(\"Sz\", IndexUp, IndexDn);\n");
      printf("  - OppU  = NewOperator(\"U\", IndexUp, IndexDn, F);\n");
      printf("\n");
      printf("\n");
      printf("\n");
    }
    if (strcmp(Help, "NewWavefunction") == 0)
    {
      printf("NewWavefunction()\n");
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
      printf("  - example: A = NewWavefunction(10,0,{{\"1010101010\",1}});\n");
      printf("\n");
      printf("\n");
    }
    if ((strcmp(Help, "Copy") == 0)||(strcmp(Help, "Clone") == 0))
    {
      printf("Copy() or Clone()\n");
      printf("  - creats a copy of a wavefunction or operator\n");
      printf("  - expects as input an object, returns a copy of that object\n");
      printf("  - note that \"=\" often does not what you want.\n");
      printf("  - A=B for operators and wave-functions copies a pointer to the object.\n");
      printf("  - A change in A thus automatically changes B\n");
      printf("\n");
      printf("  - example: A=Copy(B);\n");
      printf("\n");
      printf("\n");
    }
    if ((strcmp(Help, "Conjugate") == 0))
    {
      printf("Conjugate() or C()\n");
      printf("  - creats a copy of a wavefunction or operator and conjugates this copy\n");
      printf("  - expects as input an object, returns a copy of that object\n");
      printf("\n");
      printf("  - example: A=Conjugate(B);\n");
      printf("\n");
      printf("\n");
    }
    if ((strcmp(Help, "ConjugateTranspose") == 0))
    {
      printf("ConjugateTranspose()\n");
      printf("  - creats a copy of a wavefunction or operator and conjugates transposes this copy\n");
      printf("  - expects as input an object, returns a copy of that object\n");
      printf("\n");
      printf("  - example: A=ConjugateTranspose(B);\n");
      printf("\n");
      printf("\n");
    }
    if ((strcmp(Help, "Transpose") == 0))
    {
      printf("Transpose()\n");
      printf("  - creats a copy of a wavefunction or operator and transposes this copy\n");
      printf("  - expects as input an object, returns a copy of that object\n");
      printf("\n");
      printf("  - example: A=Transpose(B);\n");
      printf("\n");
      printf("\n");
    }
    if (strcmp(Help, "Rotate") == 0)
    {
      printf("Rotate()\n");
      printf("\n");
      printf("\n");
    }
    if (strcmp(Help, "EvolveToGroundstate") == 0)
    {
      printf("EvolveToGroundstate()\n");
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
      printf("  - example: EvolveToGroundstate(H,{psi1,psi2},{{\"epsilon\",0.000001},{\"restrictions\",{10,0,{{\"11111 00000\",5,5}}}});\n");
      printf("\n");
      printf("\n");
    }
    if (strcmp(Help, "Eigensystem") == 0)
    {
      printf("Eigensystem()\n");
      printf("  - calculates the lowest eigen functions of a given Hamiltonian. The number of states to calculate is given as an input\n");
      printf("  - expects as input an operator and either\n");
      printf("    - a (list of) wavefunction(s) that specify a starting point. (filling)\n");
      printf("    - note that the final state found might be of different symmetry than the starting point\n");
      printf("    - note also that the lowest state might not be found if the starting point excludes certain symmetries\n");
      printf("    - look at EvolveToGroundstate if one wants to track a specific state\n");
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
    if (strcmp(Help, "CreateSpectra") == 0)
    {
      printf("CreateSpectra()\n");
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
    if (strcmp(Help, "CreateFluorescenceYield") == 0)
    {
      printf("CreateFluorescenceYield()\n");
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
  return luaL_error(L, "ERROR: first argument in NewOperator (%s) is not recognized as a known operatortype\n See Help(\"NewOperator\") for more information\n",OT);
}

/*  on stack: string describing the symmetry
 *  integer describing the angular momentum
 *  table with eigen energies of type double
 *
 *  output table of Harmonics Expansion
 */
static int LuaPotentialExpandedOnClm(lua_State *L)
{
  char *s;
  double * eigenvalues;
  HarmonicsExpansion Akm;
  uint32_t l,len;
  s=(char*)luaL_checkstring(L, 1);
  l=luaL_checkint(L,2);
  eigenvalues=luaL_checkarray_double(L, 3, &len);
  if(PotentialExpandedOnC(s, l, eigenvalues, &Akm))
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
  luaL_setmetatable(L, "Wavefunction_Type");
  psi->NFermion=(uint32_t)luaL_checkunsigned(L, 1);
  psi->NBoson  =(uint32_t)luaL_checkunsigned(L, 2);
  
  if(psi->NBoson!=0)
  {
    return luaL_error(L, "Bosons not included in this version\n Sorry\n");
  }
  
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
        temp=*(Complex16 *)luaL_checkudata(L, -1, "Complex_Type");
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
      if (luaL_testudata(L, len, "Complex_Type") != NULL)
      {
        Complex16 *a = (Complex16 *) luaL_checkudata(L, len, "Complex_Type");
        Complex16 *res;
        res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
        res->Re = a->Re;
        res->Im = a->Im;
        luaL_setmetatable(L, "Complex_Type");
        return 1;
      }
      if (luaL_testudata(L, len, "Wavefunction_Type") != NULL)
      {
        WaveFunctionType *psi = (WaveFunctionType*) luaL_checkudata(L, len, "Wavefunction_Type");
        WaveFunctionType *psiNew;
        psiNew = (WaveFunctionType *) lua_newuserdata(L, sizeof(WaveFunctionType));
        WaveFunctionCopy(*psi, psiNew);
        luaL_setmetatable(L, "Wavefunction_Type");
        return 1;
      }
      if (luaL_testudata(L, len, "Operator_Type") != NULL)
      {
        OperatorType *O = (OperatorType*) luaL_checkudata(L, len, "Operator_Type");
        OperatorType *ONew;
        ONew = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
        OperatorCopy(ONew, *O);
        luaL_setmetatable(L, "Operator_Type");
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


/*  copy or clone the element on the lua stack and place it on the stack
 *  set values smaller then epsilon to zero
 *  works recursively on tables
 *  allocates new memory for user types
 *
 *  if there is one ellement on the stack clone this element using epsilon=2*DBL_EPSILON
 *  if there are two ellements on the stack and the last is a number use the last as epsilon
 *  if there are two ellements on the stack and the last is a list use this list as an option list
 *  if there are three ellements on the stack and the second is a number and the last a list, use the number as epsilon and the list as options
 */
static int LuaChop(lua_State * L)
{
  double epsilon=2*DBL_EPSILON;
  double temp;
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
      temp=lua_tonumber(L,len);
      if(Abs(temp)<epsilon)temp=0;
      lua_pushnumber(L,temp);
      break;
    case LUA_TSTRING:
      lua_pushlstring(L,lua_tostring(L,len),lua_strlen(L,len));
      break;
    case LUA_TLIGHTUSERDATA:
      lua_pushlightuserdata(L,(void*)lua_touserdata(L,len));
      break;
    case
    LUA_TUSERDATA:
      if (luaL_testudata(L, len, "Complex_Type") != NULL)
      {
        Complex16 *a = (Complex16 *) luaL_checkudata(L, len, "Complex_Type");
        Complex16 *res;
        res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
        res->Re = ((Abs(a->Re)>epsilon)?a->Re:0);
        res->Im = ((Abs(a->Im)>epsilon)?a->Im:0);
        luaL_setmetatable(L, "Complex_Type");
        return 1;
      }
      if (luaL_testudata(L, len, "Wavefunction_Type") != NULL)
      {
        WaveFunctionType *psi = (WaveFunctionType*) luaL_checkudata(L, len, "Wavefunction_Type");
        WaveFunctionType *psiNew;
        psiNew = (WaveFunctionType *) lua_newuserdata(L, sizeof(WaveFunctionType));
        WaveFunctionCopy(*psi, psiNew);
        WaveFunctionTruncateBasis(psiNew, epsilon);
        luaL_setmetatable(L, "Wavefunction_Type");
        return 1;
      }
      if (luaL_testudata(L, len, "Operator_Type") != NULL)
      {
        OperatorType *O = (OperatorType*) luaL_checkudata(L, len, "Operator_Type");
        OperatorType *ONew;
        ONew = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
        OperatorCopy(ONew, *O);
        OperatorChop(ONew, epsilon);
        luaL_setmetatable(L, "Operator_Type");
        return 1;
      }
      return luaL_error(L, "Error: Chop an object of unrecognised type\n");
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
      return luaL_error(L, "Error: Chop an object of unrecognised type\n");
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
      if (luaL_testudata(L, len, "Complex_Type") != NULL)
      {
        Complex16 *a = (Complex16 *) luaL_checkudata(L, len, "Complex_Type");
        Complex16 *res;
        res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
        res->Re = a->Re;
        res->Im = -a->Im;
        luaL_setmetatable(L, "Complex_Type");
        return 1;
      }
      if (luaL_testudata(L, len, "Wavefunction_Type") != NULL)
      {
        WaveFunctionType *psi = (WaveFunctionType*) luaL_checkudata(L, len, "Wavefunction_Type");
        WaveFunctionType *psiNew;
        psiNew = (WaveFunctionType *) lua_newuserdata(L, sizeof(WaveFunctionType));
        WaveFunctionCopy(*psi, psiNew);
        WaveFunctionConjugate(*psiNew);
        luaL_setmetatable(L, "Wavefunction_Type");
        return 1;
      }
      if (luaL_testudata(L, len, "Operator_Type") != NULL)
      {
        OperatorType *O = (OperatorType*) luaL_checkudata(L, len, "Operator_Type");
        OperatorType *ONew;
        ONew = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
        OperatorCopy(ONew, *O);
        OperatorComplexConjugate(ONew);
        luaL_setmetatable(L, "Operator_Type");
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
      if (luaL_testudata(L, len, "Complex_Type") != NULL)
      {
        Complex16 *a = (Complex16 *) luaL_checkudata(L, len, "Complex_Type");
        Complex16 *res;
        res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
        res->Re = a->Re;
        res->Im = -a->Im;
        luaL_setmetatable(L, "Complex_Type");
        return 1;
      }
      if (luaL_testudata(L, len, "Wavefunction_Type") != NULL)
      {
        WaveFunctionType *psi = (WaveFunctionType*) luaL_checkudata(L, len, "Wavefunction_Type");
        WaveFunctionType *psiNew;
        psiNew = (WaveFunctionType *) lua_newuserdata(L, sizeof(WaveFunctionType));
        WaveFunctionCopy(*psi, psiNew);
        WaveFunctionConjugate(*psiNew);
        luaL_setmetatable(L, "Wavefunction_Type");
        return 1;
      }
      if (luaL_testudata(L, len, "Operator_Type") != NULL)
      {
        OperatorType *O = (OperatorType*) luaL_checkudata(L, len, "Operator_Type");
        OperatorType *ONew;
        ONew = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
        OperatorConjugateTranspose(*O, ONew);
        luaL_setmetatable(L, "Operator_Type");
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
      if (luaL_testudata(L, len, "Complex_Type") != NULL)
      {
        Complex16 *a = (Complex16 *) luaL_checkudata(L, len, "Complex_Type");
        Complex16 *res;
        res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
        res->Re = a->Re;
        res->Im = a->Im;
        luaL_setmetatable(L, "Complex_Type");
        return 1;
      }
      if (luaL_testudata(L, len, "Wavefunction_Type") != NULL)
      {
        WaveFunctionType *psi = (WaveFunctionType*) luaL_checkudata(L, len, "Wavefunction_Type");
        WaveFunctionType *psiNew;
        psiNew = (WaveFunctionType *) lua_newuserdata(L, sizeof(WaveFunctionType));
        WaveFunctionCopy(*psi, psiNew);
        luaL_setmetatable(L, "Wavefunction_Type");
        return 1;
      }
      if (luaL_testudata(L, len, "Operator_Type") != NULL)
      {
        OperatorType *O = (OperatorType*) luaL_checkudata(L, len, "Operator_Type");
        OperatorType *ONew;
        ONew = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
        OperatorTranspose(*O, ONew);
        luaL_setmetatable(L, "Operator_Type");
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
    return luaL_error(L, "Rotate called with %d arguments, 2 expected\n",len);
  }
  if (luaL_testudata(L, len-1, "Wavefunction_Type") != NULL)
  {
    U = luaL_checkMatrix(L, len);
    WaveFunctionType *psi = (WaveFunctionType*) luaL_checkudata(L, len-1, "Wavefunction_Type");
    if(U->Ni!=U->Nj) luaL_error(L, "Error while rotating wavefunction. Rotation matrix should be square (Ni=%d Nj=%d)",U->Ni,U->Nj);
    if(U->Ni!=psi->NFermion+psi->NBoson) luaL_error(L, "Error while rotating wavefunction. Rotation matrix should have dimension (%d) of number of fermions (%d) + bosons (%d) of wavefunction",U->Ni,psi->NFermion,psi->NBoson);
    WaveFunctionType *psiNew;
    psiNew = (WaveFunctionType *) lua_newuserdata(L, sizeof(WaveFunctionType));
    RotateWaveFunction(psiNew, *psi, *U);
    luaL_setmetatable(L, "Wavefunction_Type");
    WaveFunctionTruncateBasis(psiNew, Options.zero);
    return 1;
  }
  if (luaL_testudata(L, len-1, "Operator_Type") != NULL)
  {
    U = luaL_checkMatrix(L, len);
    OperatorType *O = (OperatorType*) luaL_checkudata(L, len-1, "Operator_Type");
    if(U->Ni!=U->Nj) luaL_error(L, "Error while rotating operator. Rotation matrix should be square (Ni=%d Nj=%d)",U->Ni,U->Nj);
    if(U->Ni!=O->NFermion+O->NBoson) luaL_error(L, "Error while rotating operator. Rotation matrix should have dimension (%d) of number of fermions (%d) + bosons (%d) of operator",U->Ni,O->NFermion,O->NBoson);
    OperatorType *ONew;
    ONew = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    RotateOperatorNormalOrder(ONew, *O, *U,&QDetNotOrderedFcFaBaBc);
    luaL_setmetatable(L, "Operator_Type");
    OperatorChop(ONew, Options.zero);
    return 1;
  }
  if(lua_type(L,len-1)==LUA_TTABLE)
  {
    return luaL_error(L, "Rotate called with a table, not yet implemented\n");
    
    return 0;
  }
  return luaL_error(L, "Rotate called unrecognized type. Can only rotate wavefunctions operators and tables of these\n");
}



static int LuaCreateSpectra(lua_State *L)
{
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
  LuaOptionsType LuaCreateSpectraOptions[]={
    {"NTri","uint32_t",&(Options.NTri)},
    {"epsilon","double",&(Options.epsilon)},
    {"restrictions", "DeterminantRestrictionType", &(Options.DeterminantRestrictions)},
    {"Emin","double",&(Options.Emin)},
    {"Emax","double",&(Options.Emax)},
    {"NE","uint32_t",&(Options.NE)},
    {"Gamma","double",&(Options.Gamma)},
    {NULL, NULL, NULL}};
  /* check number of input parameters */
  len = (uint32_t)lua_gettop(L);
  if(!((len==3)||(len==4)))
  {
    return luaL_error(L, "Error CreateSpectra called with %d arguments, 3 or 4 expected\n",len);
  }
  /* read operator and wavefunctions */
  H = (OperatorType *) luaL_checkudata(L, 1, "Operator_Type");
  if(lua_type(L,2)==LUA_TTABLE)
  {
    T=luaL_checkarray_OperatorType(L, 2, &Options.NT);
  }
  else
  {
    Options.NT = 1;
    T = (OperatorType *) luaL_checkudata(L, 2, "Operator_Type");
  }
  if(lua_type(L,3)==LUA_TTABLE)
  {
    psi=luaL_checkarray_WaveFunctionType(L, 3, &Options.Npsi);
  }
  else
  {
    Options.Npsi = 1;
    psi = (WaveFunctionType *) luaL_checkudata(L, 3, "Wavefunction_Type");
  }
  /* read options */
  if(len==4)
  {
    LuaCheckOptions(L, -1, LuaCreateSpectraOptions);
    if(Options.DeterminantRestrictions.N!=0)Options.restrictions=True;
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
  luaL_setmetatable(L, "Spectra_Type");
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
 *     options
 *        {E1min, E1max, N1, Gamma1}
 *        {E2min, E2max, N2, Gamma2}
 *        NTri1
 *        NTri2
 *        epsilon
 *        restrictions1
 *        restrictions2
 */
static int LuaCreateResonantSpectra(lua_State *L)
{
  struct
  {
    ConvergenceSettingsType Settings;
    uint32_t Npsi;
    uint32_t NT1;
    uint32_t NT2;
  } Options;
  Options.Settings=InitConvergenceSettings();
  ConvergenceSettingsType SettingsNow;
  SpectraTypeArray *Spec;
  SpectraType *SpecT1psi;
  TriDiagonalMatrixType *SpecTri;
  OperatorType *H1,*H2,*T1, *T2;
  WaveFunctionType *psi;
  uint32_t iT,ipsi,len,i;
  LuaOptionsType LuaCreateResonantSpectraOptions[]={
    {"NTri1","uint32_t",&(Options.Settings.NTri1)},
    {"NTri2","uint32_t",&(Options.Settings.NTri2)},
    {"restrictions1", "DeterminantRestrictionType", &(Options.Settings.restrictions1)},
    {"restrictions2", "DeterminantRestrictionType", &(Options.Settings.restrictions2)},
    {"epsilon","double",&(Options.Settings.epsilon)},
    {"Emin1","double",&(Options.Settings.Emin1)},
    {"Emax1","double",&(Options.Settings.Emax1)},
    {"NE1","uint32_t",&(Options.Settings.NE1)},
    {"Gamma1","double",&(Options.Settings.Gamma1)},
    {"Emin2","double",&(Options.Settings.Emin2)},
    {"Emax2","double",&(Options.Settings.Emax2)},
    {"NE2","uint32_t",&(Options.Settings.NE2)},
    {"Gamma2","double",&(Options.Settings.Gamma2)},
    {NULL, NULL, NULL}};
  /* check number of input parameters */
  len = (uint32_t)lua_gettop(L);
  if(!((len==5)||(len==6)))
  {
    return luaL_error(L, "Error LuaCreateResonantSpectra called with %d arguments, 5 or 6 expected (H1, H2, T1, T2^+, psi, {options})\n",len);
  }
  /* read operator and wavefunctions */
  H1 = (OperatorType *) luaL_checkudata(L, 1, "Operator_Type");
  H2 = (OperatorType *) luaL_checkudata(L, 2, "Operator_Type");
  if(lua_type(L,3)==LUA_TTABLE)
  {
    T1=luaL_checkarray_OperatorType(L, 3, &Options.NT1);
  }
  else
  {
    Options.NT1 = 1;
    T1 = (OperatorType *) luaL_checkudata(L, 3, "Operator_Type");
  }
  if(lua_type(L,4)==LUA_TTABLE)
  {
    T2=luaL_checkarray_OperatorType(L, 4, &Options.NT2);
  }
  else
  {
    Options.NT2 = 1;
    T2 = (OperatorType *) luaL_checkudata(L, 4, "Operator_Type");
  }
  if(lua_type(L,5)==LUA_TTABLE)
  {
    psi=luaL_checkarray_WaveFunctionType(L, 5, &Options.Npsi);
  }
  else
  {
    Options.Npsi = 1;
    psi = (WaveFunctionType *) luaL_checkudata(L, 5, "Wavefunction_Type");
  }
  /* read options */
  if (len==6)
  {
    LuaCheckOptions(L, -1, LuaCreateResonantSpectraOptions);
  }
  /* allocate spectra and and pointer to return matrix on lua stack */
  Spec = (SpectraTypeArray *) lua_newuserdata(L, sizeof(SpectraTypeArray));
  luaL_setmetatable(L, "Spectra_Type");
  Spec->NSpectra=0;
  Spec->Spec=NULL;
  for(ipsi=0;ipsi<Options.Npsi;ipsi++)
  {
    for(iT=0;iT<Options.NT1;iT++)
    {
      SettingsNow = Options.Settings;
      ResonantSpectra(*H1, *H2, T1[iT], T2, Options.NT2, psi[ipsi], &SpecTri, &SpecT1psi, &(SettingsNow));
      Spec->NSpectra+=(SettingsNow.NE1+1)*Options.NT2;
      Spec->Spec=realloc(Spec->Spec, Spec->NSpectra*sizeof(SpectraType));
      for(i=0;i<(SettingsNow.NE1+1)*Options.NT2;i++)
      {
        FreeTriDiagonalMatrix(&(SpecTri[i]));
      }
      free(SpecTri);
      for(i=0;i<(SettingsNow.NE1+1)*Options.NT2;i++)
      {
        Spec->Spec[i+Spec->NSpectra-(SettingsNow.NE1+1)*Options.NT2] = SpecT1psi[i];
      }
      free(SpecT1psi);
    }
  }
  return 1;
}





static int LuaCreateFluorescenceYield(lua_State *L)
{
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
  LuaOptionsType LuaCreateSpectraOptions[]={
    {"NTri","uint32_t",&(Options.NTri)},
    {"epsilon","double",&(Options.epsilon)},
    {"restrictions", "DeterminantRestrictionType", &(Options.DeterminantRestrictions)},
    {"Emin","double",&(Options.Emin)},
    {"Emax","double",&(Options.Emax)},
    {"NE","uint32_t",&(Options.NE)},
    {"Gamma","double",&(Options.Gamma)},
    {NULL, NULL, NULL}};

  /* check number of input parameters */
  len = (uint32_t)lua_gettop(L);
  if(!((len==4)||(len==5)))
  {
    return luaL_error(L, "Error CreateFluorescenceYield called with %d arguments, 4 or 5 expected\n",len);
  }
  /* read operator and wavefunctions */
  H = (OperatorType *) luaL_checkudata(L, 1, "Operator_Type");
  if(lua_type(L,2)==LUA_TTABLE)
  {
    T1=luaL_checkarray_OperatorType(L, 2, &Options.NT1);
  }
  else
  {
    Options.NT1 = 1;
    T1 = (OperatorType *) luaL_checkudata(L, 2, "Operator_Type");
  }
  if(lua_type(L,3)==LUA_TTABLE)
  {
    T2=luaL_checkarray_OperatorType(L, 3, &Options.NT2);
  }
  else
  {
    Options.NT2 = 1;
    T2 = (OperatorType *) luaL_checkudata(L, 3, "Operator_Type");
  }
  if(lua_type(L,4)==LUA_TTABLE)
  {
    psi=luaL_checkarray_WaveFunctionType(L, 4, &Options.Npsi);
  }
  else
  {
    Options.Npsi = 1;
    psi = (WaveFunctionType *) luaL_checkudata(L, 4, "Wavefunction_Type");
  }
  /* read options */
  if(len==5)
  {
    LuaCheckOptions(L, -1, LuaCreateSpectraOptions);
    if(Options.DeterminantRestrictions.N!=0)Options.restrictions=True;
  }
  /* allocate spectra and and pointer to return matrix on lua stack */
  Spec = (SpectraTypeArray *) lua_newuserdata(L, sizeof(SpectraTypeArray));
  luaL_setmetatable(L, "Spectra_Type");
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





static int LuaEvolveToGroundstate(lua_State * L)
{
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
  LuaOptionsType LuaEvolveToGroundStateOptions[]={
    {"epsilon","double",&(Options.epsilon)},
    {"zero","double",&(Options.zero)},
    {"restrictions", "DeterminantRestrictionType", &(Options.DeterminantRestrictions)},
    {NULL, NULL, NULL}};
  /* check number of input parameters */
  len = (uint32_t)lua_gettop(L);
  if(!((len==2)||(len==3)))
  {
    return luaL_error(L, "Error EvolveToGroundstate called with %d arguments, 2 or 3 expected\n",len);
  }
  /* read operator and wavefunctions */
  Opp = (OperatorType *) luaL_checkudata(L, 1, "Operator_Type");
  if(lua_type(L,2)==LUA_TTABLE)
  {
    psi=luaL_checkarray_WaveFunctionType(L, 2, &Options.Npsi);
  }
  else
  {
    Options.Npsi = 1;
    psi = (WaveFunctionType *) luaL_checkudata(L, 2, "Wavefunction_Type");
  }
  /* read options */
  if(len==3)
  {
    LuaCheckOptions(L, -1, LuaEvolveToGroundStateOptions);
    if(Options.DeterminantRestrictions.N!=0)Options.restrictions=True;
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
      temp = (WaveFunctionType *)luaL_checkudata(L, -1, "Wavefunction_Type");
      *temp = psi[i];
      lua_pop(L, 1);
    }
  }
  return 0;
}

/**
 * Calculates all eigen values and functions of a matrix stored as a table of tables on the stack
 * @param L
 *        Pointer to a lua state. The last element on the stack needs to be a square table of tables of numbers (real or complex)
 * @return 2
 *         The first element added on the stack is a table of the eigenvalues
 *         The second element added on the stack is a table of the eigen functions
 */
static int LuaEigensystemDenseTable(lua_State *L)
{
  MatrixType *H,fun;
  double *val;
  H=luaL_checkMatrix(L,1);
  val=malloc(H->Ni*sizeof(double));
  Eigensystem(*H,&fun,val);
  FreeMatrix(H);
  LuaPushArray_double(L, val, fun.Ni);
  LuaPushTable_MatrixType(L, &fun);
  return 2;
}

/**
 * Calculates the lowest N eigenstates of a given operator
 * @param L
 *        Pointer to a lua state
 *        The Lua state can have:
 *        - length 1, being a square table of table of real or complex numbers
 *        - length 2, being an operator and a set of starting functions
 *        - length 3, being an operator, a set of starting functions and a set of options
 *        - length 3, being an operator, a set of starting restrictions, and a number defining the number of eigen states wanted
 *        - length 4, being an operator, a set of starting restrictions, a number defining the number of eigen states wanted and a set of options
 * @return 0 or 1
 *         if called with a set of starting wavefunctions returns 0 and alters the functions
 *         if called with a set of restrictions returns 1 and the newly created wavefunctions
 */
static int LuaEigensystem(lua_State * L)
{
  void *buff;
  struct
  {
    ConvergenceSettingsType Settings;
    Bool restrictions;
    Bool ExpandBasis;
    uint32_t NBitsKey;
    DeterminantRestrictionType DeterminantRestrictions;
  } Options;
  Options.Settings=InitConvergenceSettings();
  Options.restrictions=False;
  Options.ExpandBasis=True;
  Options.DeterminantRestrictions.N=0;
  Options.NBitsKey=20;
  OperatorType *Opp;
  DeterminantRestrictionType StartRestrictions;
  StartRestrictions.N=0;
  WaveFunctionType *psi=NULL,*temp;
  uint32_t i,len, Npsi=0,QStartRestrictions=0;
  /* set options */
  LuaOptionsType LuaEigensystemOptions[]={
    {"DenseBoarder", "uint32_t", &(Options.Settings.DenseBoarder)},
    {"NKrylovStart","uint32_t",&(Options.Settings.NKrylovStart)},
    {"NKrylovMax","uint32_t",&(Options.Settings.NKrylovMax)},
    {"NKrylovStep","uint32_t",&(Options.Settings.NKrylovStep)},
    {"NKrylovMin","uint32_t",&(Options.Settings.NKrylovMin)},
    {"NBitsKey","uint32_t",&(Options.Settings.NBitsKey)},
    {"epsilon","double",&(Options.Settings.epsilon)},
    {"zero","double",&(Options.Settings.zero)},
    {"restrictions", "DeterminantRestrictionType", &(Options.DeterminantRestrictions)},
    {"ExpandBasis", "bool",&(Options.ExpandBasis)},
    {NULL, NULL, NULL}};
  /* check number of input parameters */
  len = (uint32_t)lua_gettop(L);
  if(len==1) return LuaEigensystemDenseTable(L);
  if(!((len==2)||(len==3)||(len==4)))
  {
    return luaL_error(L, "Error Eigensystem called with %d arguments, 2, 3 or 4 expected\n",len);
  }
  /* read operator, (StartRestrictions or StartPsi) and Number of eigen-states */
  Opp = (OperatorType *) luaL_checkudata(L, 1, "Operator_Type");
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
    psi = (WaveFunctionType *) luaL_checkudata(L, 2, "Wavefunction_Type");
  }
  if((len==2)&&(QStartRestrictions))
  {
    return luaL_error(L, "Error Eigensystem called with 2 arguments, whereby the first is the operator and the second a set of restrictions\nNeed a thrid argument how many eigenstates should be calculated\n");
  }
  if(QStartRestrictions)
  {
    Npsi=(uint32_t)luaL_checkunsigned(L, 3);
    // allocate psi
    if(Npsi==1)
    {
      psi = (WaveFunctionType*) lua_newuserdata(L, sizeof(WaveFunctionType));
      luaL_setmetatable(L, "Wavefunction_Type");
    }
    else
    {
      psi=malloc(Npsi*sizeof(WaveFunctionType));
    }
  }
  /* read options */
  if((len==4)||((len==3)&&(!(QStartRestrictions))))
  {
    LuaCheckOptions(L, -1, LuaEigensystemOptions);
    if(Options.DeterminantRestrictions.N!=0)Options.restrictions=True;
  }
  /* calculate groundstate */
  if(Options.ExpandBasis==True)
  {
    BlockGroundState(*Opp, Npsi, &(Options.Settings) , StartRestrictions, Options.DeterminantRestrictions, psi, Options.NBitsKey);
  }
  else
  {
    BlockLanczosGroundStateConserveBasis(*Opp, psi, Npsi, &(Options.Settings));
  }
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
        temp = (WaveFunctionType *)luaL_checkudata(L, -1, "Wavefunction_Type");
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
      if (luaL_testudata(L, 1, "Wavefunction_Type") != NULL)
      {
        psi = (WaveFunctionType *) luaL_checkudata(L, 1, "Wavefunction_Type");
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


static const struct luaL_Reg global_functions[] = {
  {"Help", LuaHelp},
  {"Verbosity",LuaVerbosity},
  {"TimeStart",LuaTimeStart},
  {"TimeEnd",LuaTimeEnd},
  {"TimePrint",LuaTimePrint},
  {"PotentialExpandedOnClm", LuaPotentialExpandedOnClm},
  {"NewOperator", LuaNewOperator},
  {"NewWavefunction", LuaNewWaveFunction},
  {"DensityMatrix",LuaDensityMatrix},
  {"Copy", Luacopy},
  {"Clone", Luacopy},
  {"Chop", LuaChop},
  {"Conjugate", Luaconjugate},
  {"ConjugateTranspose", Luaconjugatetranspose},
  {"Transpose", Luatranspose},
  {"Rotate", Luarotate},
  {"EvolveToGroundstate", LuaEvolveToGroundstate},
  {"Eigensystem", LuaEigensystem},
  {"CreateSpectra", LuaCreateSpectra},
  {"CreateFluorescenceYield", LuaCreateFluorescenceYield},
  {"CreateResonantSpectra", LuaCreateResonantSpectra},
  {NULL, NULL}
};

void luaconstants(lua_State *L)
{
  Complex16 *luaI;
  luaI = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
  luaI->Re = 0.;
  luaI->Im = 1.;
  luaL_setmetatable(L, "Complex_Type");
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
  luaL_newmetatable(L, "Operator_Type");
  luaL_setfuncs(L, Operator_meta, 0);
  
  luaL_newmetatable(L, "Wavefunction_Type");
  luaL_setfuncs(L, Wavefunction_meta, 0);
  
  luaL_newmetatable(L, "Spectra_Type");
  luaL_setfuncs(L, Spectra_meta, 0);
    
  luaMWH_RegisterGlobalFunctions(L, global_functions);
  
  return 0;
}