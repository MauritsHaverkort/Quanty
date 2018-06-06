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
    printf("-          Quanty a script language to solve quantum mechanics          -\n");
    printf("-          (c) 1995-2015 All rights reserved                            -\n");
    printf("-          Maurits W. Haverkort                                         -\n");
    printf("-          http://www.quanty.org                                        -\n");
    printf("-------------------------------------------------------------------------\n");
    printf("\n\n");
    printf("The script language is based on Lua (www.lua.org) with several modifications\nThe standard Lua functions and extensions should work.\n");
    printf("\n\n");
    printf("Available functions in Quanty are:\n");
    printf("\n");
    printf("Help() or Help(FunctionName_string)\n");
    printf("  - shows this file\n");
    printf("  - excepts as input one of the other functions in order to show help on that function\n");
    printf("  - example: Help(\"NewOperator\");\n");
    printf("\n");
    printf("Verbosity(_integer)\n");
    printf("  - Set the amount of output the program produces\n");
    printf("  - Verbosity(0) suppresses all output, Verbosity(0xFFFF) shows maximum output\n");
    printf("TimeStart(Name_string)\n");
    printf("  - Starts a timing sequence with the name \"Name\"\n");
    printf("TimeEnd(string_)\n");
    printf("  - Ends a timing sequence with the name \"Name\"\n");
    printf("  - A timing sequence with the same name should have been started otherwise an error is produced\n");
    printf("TimePrint()\n");
    printf("  - Prints an overview of all timings\n");
    printf("PotentialExpandedOnClm(PointGroupName_string,l_integer,Energy_table_real)\n");
    printf("  - creates a list representing a potential expanded on spherical Harmonics for a given symmetry and eigen energy of the irreducible representations\n");
    printf("  - The potential V(theta, phi) is given as:\n");
    printf("  - V(theta, phi) = sum_k,m A_{k,m} Y_km(theta,phi)\n");
    printf("  - the list saved is {{k1,m1,Ak1m1},{k2,m2,Ak2m2},...}\n");
    printf("  - input is the point-group symmetry, the angular momentum and the eigen energies of the irreducible representations\n");
    printf("  - PotentialExpandedOnClm(\"Oh\",2,{0.6,-0.4});\n");
    printf("  - for a potential in cubic symmetry acting on the d-shell with the energy of the eg's as 0.6 and t2g's at -0.4 (10Dq=1.0)\n");
    printf("\n");
    printf("NewOperator(Name_string,OtherSettings___)\n");
    printf("  - creates a new operator\n");
    printf("  - expects as input on the first position a string describing the type of operator\n");
    printf("  - other parameters depend on the type of operator defined in the first string\n");
    printf("  - example: NewOperator(\"Sz\", NFermionicModes, IndexUp, IndexDn);\n");
    printf("\n");
    printf("NewWavefunction(NFermionicModes_integer,NBosonicModes_integer,{{Determinant_string,Prefactor_Complex}, ...})\n");
    printf("  - creates a new wavefunction\n");
    printf("  - expects as input the number of fermionic modes, number of bosonic modes in the one particle basis, and a table of determinants with the prefactor\n");
    printf("  - example: NewWavefunction(10,0,{{{0x80,0x00},0.7},{{0x04,0x00},0.7}});\n");
    printf("  - example: NewWavefunction(10,0,{{\"1010101010\",0.7},{\"0101010101\",0.7}});\n");
    printf("\n");
    printf("DensityMatrix(psi_wavefunction, index_table_integer)\n");
    printf("  - calculated the density matrix including the spin-orbitals at position given by the index table\n");
    printf("Copy(A_freetype) or Clone(A_freetype)\n");
    printf("  - creates a copy of A\n");;
    printf("\n");
    printf("Conjugate(A_freetype)\n");
    printf("  - creats a copy A and conjugates this copy\n");    printf("\n");
    printf("ConjugateTranspose(A_freetype)\n");
    printf("  - creates a copy of A and conjugates transposes this copy\n");
    printf("\n");
    printf("Transpose(A_freetype)\n");
    printf("  - creates a copy of A and transposes this copy\n");
    printf("Rotate(psi_wavefunction, U_table_complex) or Rotate(O_operator, U_table_complex)\n");
    printf("  - Given an operator or a wave-function on a basis of N fermionic modex + M bosonic modes\n");
    printf("  - and a matrix U of size N+M\n");
    printf("  - creates a new operator on a rotated basis\n");
    printf("  - example: Rotate(Operator,RotationMatrix);\n");
    printf("\n");
    printf("Eigensystem(H_operator, StartRestrictions_restrictions, N_integer) or Eigensystem(H_operator, psi_wavefunction_table)\n");
    printf("  - calculates the lowest N eigenfunctions of a given operator\n");
    printf("  - example Eigensystem(OppSy,{psi1,psi2});\n");
    printf("  - example Eigensystem(OppSy,Restrictions,Nspi);\n");
    printf("\n");
    printf("CreateSpectra(H_operator, T_operator(_table), psi_wavefunction(_table))\n");
    printf("  - calculates spectra of a giving Hamiltonian, starting from the specified state and using a specific transition operator\n");
    printf("  - expects as input an operator a (list of) transition operator(s) and a (list of) wavefunction(s)\n");
    printf("  - The program calculates <psi | T^+ 1/(w-H+i G/2) T | psi>)\n");
    printf("\n");
    printf("CreateFluorescenceYield(H_operator, T1_operator(_table), T2_operator(_table), psi_wavefunction(_table)))\n");
    printf("  - calculates the fluorescence yield spectra\n");
    printf("  - calculated is <psi | T1^+ 1/(w-H-i G/2) T2 T2^+ 1/(w-H+i G/2) T1 | psi>\n");
    printf("  - expects as input an operator a (list of) transition operator(s) T1 a (list of) transition operator(s) T2^+ and a (list of) wavefunction(s)\n");
    printf("\n");
    printf("CreateResonantSpectra(H1_operator, H2_operator, T1_operator(_table), T2_operator(_table), psi_wavefunction(_table)))\n");
    printf("  - calculates the fluorescence yield spectra\n");
    printf("  - calculated is <psi | T1^+ 1/(w-H1-i G/2) T2 1/(w-H2+i G/2) T2^+ 1/(w-H1+i G/2) T1 | psi>\n");
    printf("  - expects as input two operators a (list of) transition operator(s) T1 a (list of) transition operator(s) T2^+ and a (list of) wavefunction(s)\n");
    printf("\n");

    printf("\n\n");
    printf("Available objects\n");
    printf("\n");
    printf("Operator\n");
    printf("  - stores a quantum mechanical operator based on second quantization\n");
    printf("  - methods:\n");
    printf("    - O.Name()\n");
    printf("       - set or return the name of the operator\n");
    printf("    - O.Print()\n");
    printf("       - print the operator\n");
    printf("    - O.Chop()\n");
    printf("       - remove small numbers from the operator\n");
    printf("    - O.Conjugate()\n");
    printf("       - conjugate the operator\n");
    printf("    - O.ConjugateTranspose()\n");
    printf("       - conjugate transpose the operator\n");
    printf("    - O.Transpose()\n");
    printf("       - transpose the operator\n");
    printf("\n");
    printf("Wavefunction\n");
    printf("  - stores a wavefunction\n");
    printf("  - methods:\n");
    printf("    - psi.Name()\n");
    printf("       - set or return the name of the operator\n");
    printf("    - psi.Print()\n");
    printf("       - print the operator\n");
    printf("\n");
    printf("Complex\n");
    printf("  - stores a complex number\n");
    printf("  - methods:\n");
    printf("    - z.re\n");
    printf("       - set or return the real part note that z must be of type complex\n");
    printf("    - z.im\n");
    printf("       - set or return the imaginary part note that z must be of type complex\n");
    printf("  - functions:\n");
    printf("    - Complex.Re(z_complex)\n");
    printf("       - return the real part of z\n");
    printf("    - Complex.Im(z_complex)\n");
    printf("       - return the imaginary part of z\n");
    printf("\n");
    printf("spectra\n");
    printf("  - stores a spectral function\n");
    printf("  - methods:\n");
    printf("    - spec.Print()\n");
    printf("       - print the spectra\n");
    printf("    - spec.Element(i_integer_table)\n");
    printf("       - Return new spectra corresponding to the spectra at position i\n");
    printf("    - spec.Sum(alpha_double_table)\n");
    printf("       - Return new spectra corresponding to the sum of the spectra using the pre-factors as given in alpha\n");
    printf("    - spec.Broaden(FWHMGaus_double, FMHMLorentz_double)\n");
    printf("       - Broaden the spectra by a Gaussian and Lorentzian\n");
    printf("    - spec.Shift(x_double)\n");
    printf("       - add x to the energy\n");
    printf("    - spec.MultiplyByEnergy()\n");
    printf("       - multiply the intensity by the energy\n");
    printf("  - functions:\n");
    printf("    - Spectra.Copy() or Spectra.Clone()\n");
    printf("       - copy the spectra\n");
    printf("    - Spectra.Element(i_integer_table)\n");
    printf("       - Return new spectra corresponding to the spectra at position i\n");
    printf("    - Spectra.Sum(alpha_double_table)\n");
    printf("       - Return new spectra corresponding to the sum of the spectra using the pre-factors as given in alpha\n");
    printf("    - Spectra.Broaden(FWHMGaus_double, FMHMLorentz_double)\n");
    printf("       - Return new spectra corresponding to the broadened spectra by a Gaussian and Lorentzian\n");
    printf("    - Spectra.Shift(x_double)\n");
    printf("       - add x to the energy\n");
    printf("    - Spectra.MultiplyByEnergy()\n");
    printf("       - multiply the intensity by the energy\n");
    printf("    - Spectra.ToTable()\n");
    printf("       - return a table representing the spectra\n");
    printf("\n");

  }
  if(len==1)
  {
    char const *Help =luaL_checkstring(L, 1);
    if (strcmp(Help, "NewOperator") == 0)
    {
      printf("NewOperator()\n");
      printf("  - creates a new operator\n");
      printf("  - expects as input on the first position a string describing the type of operator\n");
      printf("  - other parameters depend on the type of operator defined in the first string\n");
      printf("  - possible operators and additional parameters are:\n");
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
      printf("  - F={F0,F2,F4}\n");
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
      printf("    - This number should be larger than the number of Determinants in the basis for optimal performance\n");
      printf("    - Standard value is 16 (65 536) which works good for small calculations upto basis sets of 25 000\n");
      printf("    - If more determinants are needed it is recommended to enlarge\n");
      printf("    - 20 (  1 048 576) for basis sets upto    500 000\n");
      printf("    - 24 ( 16 777 216) for basis sets upto 10 000 000\n");
      printf("    - for even larger basis sets architecture dependent tests are advised.\n");
      printf("    - note that results are independent of this value, only calculation time and memory requirement changes.\n");
      printf("\n");
      printf("  - example: A = NewWavefunction(10,0,{{\"1010101010\",1}});\n");
      printf("\n");
      printf("\n");
    }
    if ((strcmp(Help, "Copy") == 0)||(strcmp(Help, "Clone") == 0))
    {
      printf("Copy() or Clone()\n");
      printf("  - creates a copy of a wavefunction or operator\n");
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
      printf("  - creates a copy of a wavefunction or operator and conjugates this copy\n");
      printf("  - expects as input an object, returns a copy of that object\n");
      printf("\n");
      printf("  - example: A=Conjugate(B);\n");
      printf("\n");
      printf("\n");
    }
    if ((strcmp(Help, "ConjugateTranspose") == 0))
    {
      printf("ConjugateTranspose()\n");
      printf("  - creates a copy of a wavefunction or operator and conjugates transposes this copy\n");
      printf("  - expects as input an object, returns a copy of that object\n");
      printf("\n");
      printf("  - example: A=ConjugateTranspose(B);\n");
      printf("\n");
      printf("\n");
    }
    if ((strcmp(Help, "Transpose") == 0))
    {
      printf("Transpose()\n");
      printf("  - creates a copy of a wavefunction or operator and transposes this copy\n");
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
      printf("    - the number of lowest eigen states to calculate\n");
      printf("  - additional options can be given\n");
      printf("\n");
      printf("  - options are:\n");
      printf("    - \"zero\"\n");
      printf("      - Determines how far and good the states are converged. Standard setting 100*DBL_EPSILON, should not be smaller than DBL_EPSILON\n");
      printf("    - \"epsilon\"\n");
      printf("      - Determines the smallest prefactor of a determinant that is kept in the wave-function. Standard setting 10000*DBL_EPSILON\n");
      printf("    - \"restrictions\"\n");
      printf("      - restricts the basis by excluding certain occupations\n");
      printf("      - Example for a ligand field calculation one can require to have no more than 1 hole in the Ligand shell (between 9 and 10 electrons)\n");
      printf("      - {\"restrictions\",{20,0,{{\"0000000000 1111111111\",9,10}}}}\n");
      printf("    - \"NBitsKey\"\n");
      printf("      + only when a set of restrictions and Npsi is given as an input, otherwise NBitsKey is taken from the wave function\n");
      printf("      - number of elements in the Hash lookup table. (2^NBitsKey)\n");
      printf("      - This number should be larger than the number of Determinants in the basis for optimal performance\n");
      printf("      - Standard value is 16 (65 536) which works good for small calculations upto basis sets of 25 000\n");
      printf("      - If more determinants are needed it is recommended to enlarge\n");
      printf("      - 20 (  1 048 576) for basis sets upto    500 000\n");
      printf("      - 24 ( 16 777 216) for basis sets upto 10 000 000\n");
      printf("      - for even larger basis sets architecture dependent tests are advised.\n");
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
      printf("      - Determines the smallest prefactor of a determinant that is kept in the wave-function.\n");
      printf("    - \"NTri\"\n");
      printf("      - sets the maximum number of states in the Krylov basis\n");
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
    if (strcmp(Help, "CreateResonantSpectra") == 0)
    {
      printf("CreateResonantSpectra()\n");
      printf("  - calculates resonant spectra of a giving Hamiltonian, starting from the specified state and using a specific transition operator\n");
      printf("  - < i | T1^+ G1^+ T2 G2 T2^+ G1 T1 | i >\n");
      printf("  - expects as input two operators a (list of) transition operator(s) T1 a (list of) transition operator(s) T2^+ and a (list of) wavefunction(s), a list of options on the last position is allowed\n");
      printf("\n");
      printf("  - options are:\n");
      printf("    - \"epsilon\"\n");
      printf("      - Determines the smallest prefactor of a determinant that is kept in the wave-function.\n");
      printf("    - \"NTri1\"\n");
      printf("      - sets the maximum number of states in the first Krylov basis\n");
      printf("    - \"NTri2\"\n");
      printf("      - sets the maximum number of states in the second Krylov basis\n");
      printf("    - \"restrictions1\"\n");
      printf("      - restricts the basis for G1 by excluding certain occupations\n");
      printf("    - \"restrictions2\"\n");
      printf("      - restricts the basis for G2 by excluding certain occupations\n");
      printf("    - \"NE1\"\n");
      printf("      - number of energy points in the spectrum\n");
      printf("    - \"Emin1\"\n");
      printf("      - smallest energy in the spectrum\n");
      printf("    - \"Emax1\"\n");
      printf("      - largest energy in the spectrum\n");
      printf("    - \"Gamma1\"\n");
      printf("      - Lorenzian broadening of the spectra (full with half maximum FWHM)\n");
      printf("    - \"NE2\"\n");
      printf("      - number of energy points in the spectrum\n");
      printf("    - \"Emin2\"\n");
      printf("      - smallest energy in the spectrum\n");
      printf("    - \"Emax2\"\n");
      printf("      - largest energy in the spectrum\n");
      printf("    - \"Gamma2\"\n");
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
      printf("      - Determines the smallest prefactor of a determinant that is kept in the wave-function.\n");
      printf("    - \"NTri\"\n");
      printf("      - sets the maximum number of states in the Krylov basis\n");
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
  if(lua_type(L, 1)==LUA_TNUMBER)
  {
    uint32_t i,j,lenopp,l;
    uint16_t *oadd;
    Complex16 V;
    double oaddtemp;
    struct
    {
      char *Name;
      CompactVectorType NBitsKeyDouble;
      DeterminantRestrictionType DeterminantRestrictions;
    } Options;
    OperatorType *O;
    const uint32_t len = (uint32_t)lua_gettop(L);
    LuaOptionsType LuaNewOperatorOptions[]={
      {"restrictions", "DeterminantRestrictionType", &(Options.DeterminantRestrictions)},
      {"Restrictions", "DeterminantRestrictionType", &(Options.DeterminantRestrictions)},
      {"Name" ,        "char*",                      &(Options.Name)},
      {"NBitsKey",     "Vector",                     &(Options.NBitsKeyDouble)},
      {NULL, NULL, NULL}};
    Options.DeterminantRestrictions.N=0;
    Options.Name=NULL;
    Options.NBitsKeyDouble.Ni=0;
    Options.NBitsKeyDouble.QComplex=0;
    if((len!=3)&&(len!=4)) return luaL_error(L, "ERROR: LuaNewOperator called with %u arguments, expectd 3 or 4\n",len);
    if(len==4) LuaCheckOptions(L, -1, LuaNewOperatorOptions);
    O = (OperatorType*) lua_newuserdata(L, sizeof(OperatorType));
    luaL_setmetatable(L, "Operator_Type");
    InitOperator(O);
    O->NFermion=(uint32_t)luaL_checkunsigned(L, 1);
    O->NBoson  =(uint32_t)luaL_checkunsigned(L, 2);
    O->NWord   = O->NBoson + Min(O->NFermion, 1 + (O->NFermion - 1) / 8);
    for(i=0;i<Options.NBitsKeyDouble.Ni;i++)
    {
      if(Abs(Options.NBitsKeyDouble.VRe[i])>0.5)
      {
        InitOperatorFixedLength(O, i+1, (uint32_t)round(Options.NBitsKeyDouble.VRe[i]));
      }
    }
    if(Options.Name!=NULL) strncpy(O->Name, Options.Name, MAXStrLength-1);

    O->Res=Options.DeterminantRestrictions;
    oadd= (uint16_t*)malloc(sizeof(uint16_t)*O->MaxLength);
    luaL_checktype(L, 3, LUA_TTABLE);
    lenopp = (uint32_t)lua_objlen(L, 3);
    for(i=1;i<=lenopp;i++)
    {
      lua_rawgeti(L, 3, i);
      luaL_checktype(L, -1, LUA_TTABLE);
      l = (uint32_t)lua_objlen(L, -1);
      if(l==0) return luaL_error(L, "ERROR: LuaNewOperator called an empty table as operator at position %d\n",i);
      l--;
      for(j=1;j<=l;j++)
      {
        lua_rawgeti(L, -1, j);
        oaddtemp = luaL_checknumber(L, -1);
        if(Abs(oaddtemp)>O->NFermion+O->NBoson-1) return luaL_error(L, "ERROR: Opperator %d at position %d adding particle with index %d but there are only %d (%d Fermions + %d Bosons in the basis)\n)\n",i,j,(int)oaddtemp,O->NFermion+O->NBoson,O->NFermion,O->NBoson);
        if(Abs((Abs(oaddtemp)-(uint16_t)(Abs(oaddtemp))))>10*DBL_EPSILON) return luaL_error(L, "ERROR: Opperator %d at position %d addign particle with index %f should be integer",i,j,oaddtemp);
        // negative numbers annihilate, positive create, but we also have to capture negative and postive zero
        if(signbit(oaddtemp))
        {
          oadd[l-j]=(uint16_t)(-oaddtemp);
        }
        else
        {
          oadd[l-j]=CrAnMask + (uint16_t)oaddtemp;
        }
        lua_pop(L, 1);
      }
      lua_rawgeti(L, -1, l+1);
      if(lua_type(L,-1)==LUA_TNUMBER)
      {
        V.Im=0;
        V.Re=lua_tonumber(L, -1);
      }
      else
      {
        V = *(Complex16 *) luaL_checkudata(L, -1, "Complex_Type");
      }
      lua_pop(L,1);
      if(l>O->MaxLength) InitOperatorFixedLength(O, l, 6);
      OperatorAddLaderNormalOrder(O, oadd, l, V.Re, V.Im, &QDetNotOrderedFcFaBcBa);
      lua_pop(L,1);
    }
    free(oadd);
    return 1;
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
  if (strcmp(OT, "Tx") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewTx(L);
  }
  if (strcmp(OT, "Ty") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewTy(L);
  }
  if (strcmp(OT, "Tz") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewTz(L);
  }
  if (strcmp(OT, "Tplus") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewTplus(L);
  }
  if (strcmp(OT, "Tmin") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewTmin(L);
  }
  if (strcmp(OT, "Tsqr") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewTsqr(L);
  }
  if ((strcmp(OT, "Qx") == 0)||(strcmp(OT, "Qxx") == 0))
  {
    lua_remove(L, 1);
    return LuaOperatorNewQxx(L);
  }
  if ((strcmp(OT, "Qy") == 0)||(strcmp(OT, "Qyy") == 0))
  {
    lua_remove(L, 1);
    return LuaOperatorNewQyy(L);
  }
  if ((strcmp(OT, "Qz") == 0)||(strcmp(OT, "Qzz") == 0))
  {
    lua_remove(L, 1);
    return LuaOperatorNewQzz(L);
  }
  if ((strcmp(OT, "Qxy") == 0)||(strcmp(OT, "Qyx") == 0))
  {
    lua_remove(L, 1);
    return LuaOperatorNewQxy(L);
  }
  if ((strcmp(OT, "Qxz") == 0)||(strcmp(OT, "Qzx") == 0))
  {
    lua_remove(L, 1);
    return LuaOperatorNewQxz(L);
  }
  if ((strcmp(OT, "Qyz") == 0)||(strcmp(OT, "Qzy") == 0))
  {
    lua_remove(L, 1);
    return LuaOperatorNewQyz(L);
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
  if ((strcmp(OT, "Identity") == 0)||(strcmp(OT, "I") == 0)||(strcmp(OT, "1") == 0))
  {
    lua_remove(L, 1);
    return LuaOperatorNewIdentity(L);
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
  const uint32_t len = (uint32_t)lua_gettop(L);
  struct
  {
    char *Name;
    uint32_t NBitsKey;
  } Options;
  LuaOptionsType LuaNewWaveFunctionOptions[]={
    {"Name" ,        "char*",                      &(Options.Name)},
    {"NBitsKey",     "uint32_t",                   &(Options.NBitsKey)},
    {NULL, NULL, NULL}};
  Options.Name=NULL;
  Options.NBitsKey=16;
  psi = (WaveFunctionType*) lua_newuserdata(L, sizeof(WaveFunctionType));
  luaL_setmetatable(L, "Wavefunction_Type");
  psi->NFermion=(uint32_t)luaL_checkunsigned(L, 1);
  psi->NBoson  =(uint32_t)luaL_checkunsigned(L, 2);

  if(psi->NBoson!=0)
  {
    return luaL_error(L, "Bosons not included in this version\n Sorry\n");
  }

  if((len!=3)&&(len!=4)) return luaL_error(L, "ERROR: LuaNewWaveFunction called with %u arguments, expectd 3 or 4\n",len);
  if(len==4)
  {
    if(lua_isnumber(L, 4))
    {
      psi->NBitsKey  =(uint32_t)luaL_checkunsigned(L, 4);
    }
    else
    {
      LuaCheckOptions(L, 4, LuaNewWaveFunctionOptions);
    }
  }
  psi->QComplex=1;
  psi->NBitsKey=Options.NBitsKey;
  if(WaveFunctionInit(psi))
  {
    if(L->l_G->gcrunning)lua_gc (L, LUA_GCEMERGENCYCOLLECT, 0);
    if(WaveFunctionInit(psi))
    {
      return luaL_error(L, "WaveFunctionInit failed, your most probably out of memory");
    }
  }
  if(Options.Name!=NULL)strncpy(psi->Name, Options.Name, MAXStrLength-1);
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
        if(strlen(detchar)<psi->NFermion+2*psi->NBoson)
        {
          return luaL_error(L, "Error in LuaWavefunctionNew, length of string(%d) different from NFermion(%d)+2*NBoson(%d)\n",strlen(detchar),psi->NFermion,psi->NBoson);
        }
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
      if (luaL_testudata(L, len, "Spectra_Type") != NULL)
      {
        SpectraTypeArray *a = (SpectraTypeArray *) luaL_checkudata(L, len, "Spectra_Type");
        SpectraTypeArray *res;
        res = (SpectraTypeArray *) lua_newuserdata(L, sizeof(SpectraTypeArray));
        CopySpectraArray(a, res);
        luaL_setmetatable(L, "Spectra_Type");
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


/*  copy or clone the one but last element on the lua stack and place it on the stack
 *  set values smaller then epsilon to zero
 *  works recursively on tables
 *  allocates new memory for user types
 *
 *  if there is one ellement on the stack use epsilon=10*DBL_EPSILON
 *  if there are two ellements on the stack and the last is a number use the last as epsilon
 */
static int LuaChop(lua_State * L)
{
  double epsilon=10*DBL_EPSILON;
  double temp;
  uint32_t len = (uint32_t)lua_gettop(L);
  if(len!=1)
  {
    epsilon=luaL_checknumber(L,len);
    len=len-1;
  }
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
    case LUA_TUSERDATA:
      if (luaL_testudata(L, len, "Complex_Type") != NULL)
      {
        Complex16 *a = (Complex16 *) luaL_checkudata(L, len, "Complex_Type");
        Complex16 *res;
        double reres = 0;

        if (Abs(a->Im)<=epsilon)
        {
          reres = ((Abs(a->Re)>epsilon)?a->Re:0);
          lua_pushnumber(L, reres);
        }
        else
        {
          res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
          res->Re = ((Abs(a->Re)>epsilon)?a->Re:0);
          res->Im = ((Abs(a->Im)>epsilon)?a->Im:0);
          luaL_setmetatable(L, "Complex_Type");
        }

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
      if (luaL_testudata(L, len, "Matrix_Type") != NULL)
      {
        CompactMatrixType *M = (CompactMatrixType*) luaL_checkudata(L, len, "Matrix_Type");
        CompactMatrixType *Mnew;
        Mnew = (CompactMatrixType*) lua_newuserdata(L, sizeof(CompactMatrixType));
        CompactMatrixCopy(*M, Mnew);
        CompactMatrixChop(Mnew, epsilon);
        luaL_setmetatable(L, "Matrix_Type");
        return 1;
      }
      return luaL_error(L, "Error: Chop an object of unrecognised type\n");
    case LUA_TTABLE:
      lua_newtable(L);
      lua_pushnil(L);
      while(lua_next(L, len) != 0)
      {
        lua_pushnumber(L, epsilon);
        LuaChop(L);
        lua_pushvalue(L, -4);
        lua_insert(L, -2);
        lua_settable(L, -6);
        lua_pop(L, 2);
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
        Luaconjugate(L);
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
  Bool matrix=True,Warning=False;
  uint32_t Ni=0,i,Nj=0,j;
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
        OperatorConjugateTransposeNormalOrder(*O, ONew, &QDetNotOrderedFcFaBaBc);
        luaL_setmetatable(L, "Operator_Type");
        return 1;
      }
      return luaL_error(L, "Error: copying or cloning an object of unrecognised type\n");
    case LUA_TTABLE:
      // check if we have a table of tables of equal length, in that case we treat the table of tables as a matrix and transpose the matrix
      Ni=(uint32_t)lua_objlen(L, len);
      for(i=0;i<Ni;i++)
      {
        lua_rawgeti(L, len, i+1);
        if(lua_type(L,-1)==LUA_TTABLE)
        {
          Warning=True;
          if(i==0)
          {
            Nj=(uint32_t)lua_objlen(L,-1);
          }
          else
          {
            if(Nj!=(uint32_t)lua_objlen(L,-1))
            {
              matrix=False;
            }
          }
        }
        else
        {
          matrix=False;
        }
        lua_pop(L,1);
      }
      if(Warning&(!matrix)) luaL_warning(L, "Transpose of a table containing tables, but not of equal length, transposed the elements, not the table\n");
      if(matrix)
      {
        lua_createtable(L,Nj,0);
        for(j=0;j<Nj;j++)
        {
          lua_pushnumber(L, j+1);
          lua_createtable(L, Ni, 0);
          for(i=0;i<Ni;i++)
          {
            lua_rawgeti(L, len, i+1);
            lua_rawgeti(L,-1,j+1);
            Luaconjugatetranspose(L);
            lua_pushnumber(L, i+1);
            lua_insert(L, -2);
            lua_settable(L, -5);
            lua_pop(L, 2);
          }
          lua_settable(L, -3);
        }
      }
      else
      {
        lua_newtable(L);
        lua_pushnil(L);
        while(lua_next(L, len) != 0)
        {
          Luaconjugatetranspose(L);
          lua_pushvalue(L, -3);
          lua_insert(L, -2);
          lua_settable(L, -5);
          lua_pop(L, 1);
        }
      }
      break;
    default:
      return luaL_error(L, "Error: copying or cloning an object of unrecognised type\n");
  }
  return 1;
}

static int Luatranspose(lua_State * L)
{
  Bool matrix=True,Warning=False;
  uint32_t Ni=0,i,Nj=0,j;
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
        OperatorTransposeNormalOrder(*O, ONew,&QDetNotOrderedFcFaBaBc);
        luaL_setmetatable(L, "Operator_Type");
        return 1;
      }
      return luaL_error(L, "Error: copying or cloning an object of unrecognised type\n");
    case LUA_TTABLE:
      // check if we have a table of tables of equal length, in that case we treat the table of tables as a matrix and transpose the matrix
      Ni=(uint32_t)lua_objlen(L, len);
      for(i=0;i<Ni;i++)
      {
        lua_rawgeti(L, len, i+1);
        if(lua_type(L,-1)==LUA_TTABLE)
        {
          Warning=True;
          if(i==0)
          {
            Nj=(uint32_t)lua_objlen(L,-1);
          }
          else
          {
            if(Nj!=(uint32_t)lua_objlen(L,-1))
            {
              matrix=False;
            }
          }
        }
        else
        {
          matrix=False;
        }
        lua_pop(L,1);
      }
      if(Warning&(!matrix)) luaL_warning(L, "Transpose of a table containing tables, but not of equal length, transposed the elements, not the table\n");
      if(matrix)
      {
        lua_createtable(L,Nj,0);
        for(j=0;j<Nj;j++)
        {
          lua_pushnumber(L, j+1);
          lua_createtable(L, Ni, 0);
          for(i=0;i<Ni;i++)
          {
            lua_rawgeti(L, len, i+1);
            lua_rawgeti(L,-1,j+1);
            Luatranspose(L);
            lua_pushnumber(L, i+1);
            lua_insert(L, -2);
            lua_settable(L, -5);
            lua_pop(L, 2);
          }
          lua_settable(L, -3);
        }
      }
      else
      {
        lua_newtable(L);
        lua_pushnil(L);
        while(lua_next(L, len) != 0)
        {
          Luatranspose(L);
          lua_pushvalue(L, -3);
          lua_insert(L, -2);
          lua_settable(L, -5);
          lua_pop(L, 1);
        }
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
    OperatorType **V;
    SpectraTypeArray **SonGrid;
    ResponsefunctionType **SasMatrix;
    uint32_t NSigma;
  } Sigma;
  struct
  {
    double E0;
    uint32_t Npsi;
    uint32_t NT;
    ConvergenceSettingsType Settings;
    char *stream;
    Bool Tensor;
    Bool Sigma;
    Bool KrylovBasis;
    Bool DenseMethod;
  } Options;
  Options.Settings=InitConvergenceSettings();
  Options.stream=NULL;
  Options.Tensor=False;
  Options.Sigma=False;
  Options.DenseMethod=False;
  Options.KrylovBasis=False;
  Options.E0=DBL_MAX;
  CompactMatrixType HDense, TpsiDense, EigenfunHDense;
  SpectraTypeArray *Spec;
  OperatorType *H,*T;
  WaveFunctionType *Tpsi;
  WaveFunctionType *psi, psiDenseBasis, OpsiDenseBasis;
  WaveFunctionType **psiBasis=NULL;
  uint32_t i,j,k,iT,ipsi,len, *BlockSizeN=NULL;
  double Eminpole,Emaxpole, *val, *psiEn, psiEntemp;
  char filename[MAXStrLength];
  LuaOptionsType LuaCreateSpectraOptions[]={
    {"Tensor","bool",&(Options.Tensor)},
    {"Sigma","bool",&(Options.Sigma)},
    {"KrylovBasis","bool",&(Options.KrylovBasis)},
    {"NTri","uint32_t",&(Options.Settings.NTri)},
    {"epsilon","double",&(Options.Settings.epsilon)},
    {"Epsilon","double",&(Options.Settings.epsilon)},
    {"restrictions", "DeterminantRestrictionType", &(Options.Settings.restrictions)},
    {"Restrictions", "DeterminantRestrictionType", &(Options.Settings.restrictions)},
    {"Emin","double",&(Options.Settings.Emin)},
    {"Emax","double",&(Options.Settings.Emax)},
    {"NE","uint32_t",&(Options.Settings.NE)},
    {"Gamma","double",&(Options.Settings.Gamma)},
    {"E0","double",&(Options.E0)},
    {"DenseBorder","uint32_t",&(Options.Settings.DenseBorder)},
    {"ForceDenseMethod","uint32_t",&(Options.DenseMethod)},
    {"file","char*",&(Options.stream)}, // temporary undoccumented option
    {NULL, NULL, NULL}};
  /* free garbage full cycle*/
  if(L->l_G->gcrunning) lua_gc(L, LUA_GCCOLLECT, 0);
  /* check number of input parameters */
  len = (uint32_t)lua_gettop(L);
  if(!((len==3)||(len==4)||(len=5)))
  {
    return luaL_error(L, "Error CreateSpectra called with %d arguments, 3,4 or 5 expected\n",len);
  }
  /* check if a self energy is given */
  if(len==5) Options.Sigma=True;
  if(len==4)
  {
    if(LuaTestOptionsList(L, 4)==0) Options.Sigma=True;
  }
  /* read operator and wavefunctions */
  H = (OperatorType *) luaL_checkudata(L, 1, "Operator_Type");
  if(lua_type(L,2)==LUA_TTABLE)
  {
    T=luaL_checkarray_OperatorType(L, 2, &Options.NT);
    if(Options.NT==0) return luaL_error(L, "Errof CreateSpectra called with a table of transition operators of length 0\n");
  }
  else
  {
    Options.NT = 1;
    T = (OperatorType *) luaL_checkudata(L, 2, "Operator_Type");
  }
  if(lua_type(L,3)==LUA_TTABLE)
  {
    psi=luaL_checkarray_WaveFunctionType(L, 3, &Options.Npsi);
    if(Options.Npsi==0) return luaL_error(L, "Errof CreateSpectra called with a table of wave functions of length 0\n");
  }
  else
  {
    Options.Npsi = 1;
    psi = (WaveFunctionType *) luaL_checkudata(L, 3, "Wavefunction_Type");
  }
  /* read options */
  if((len==5)||((len==4)&&(Options.Sigma==False)))
  {
    LuaCheckOptions(L, -1, LuaCreateSpectraOptions);
  }
  /* read self energy */
  if(Options.Sigma==True)
  {
    /* read self energy */
    /* the self energy can be */
    // 1 a single object of spectra type
    // 2 a single object of ResponseFunction type
    // 3 A list of length 2 containing an operator or a list of operators and a self energy in format 1 or 2
    // 4 a list of objects of form 1, 2 or 3
    Sigma.NSigma=0;
    if(lua_type(L,4)!=LUA_TTABLE)
    {
      Sigma.NSigma=1;
    }
    else if (lua_objlen(L,4)==2)
    {
      // we still can have one or two self energies
      lua_rawgeti(L,4,1);
      if(luaL_checkudata(L, -1, "Operator_Type")!=NULL)Sigma.NSigma=1;
      else if (lua_istable(L,-1))
      {
        lua_rawgeti(L,-1,(int)lua_objlen(L, -1));
        if(luaL_checkudata(L, -1, "Operator_Type")!=NULL)Sigma.NSigma=1;
        lua_pop(L,1);
      }
      lua_pop(L,1);
    }
    if(Sigma.NSigma==1)
    {
      Sigma.V= (OperatorType**)malloc(Sigma.NSigma*sizeof(OperatorType*));
      Sigma.SonGrid= (SpectraTypeArray**)malloc(Sigma.NSigma*sizeof(SpectraTypeArray*));
      Sigma.SasMatrix= (ResponsefunctionType**)malloc(Sigma.NSigma*sizeof(ResponsefunctionType*));
      LuaReadSingleSelfenergy(L, 4, Sigma.V[0], Sigma.SonGrid[0], Sigma.SasMatrix[0]);
    }
    else
    {
      Sigma.NSigma=(uint32_t)lua_objlen(L, 4);
      Sigma.V= (OperatorType**)malloc(Sigma.NSigma*sizeof(OperatorType*));
      Sigma.SonGrid= (SpectraTypeArray**)malloc(Sigma.NSigma*sizeof(SpectraTypeArray*));
      Sigma.SasMatrix= (ResponsefunctionType**)malloc(Sigma.NSigma*sizeof(ResponsefunctionType*));
      for(i=0;i<Sigma.NSigma;i++)
      {
        lua_rawgeti(L,4,i+1);
        LuaReadSingleSelfenergy(L, -1, Sigma.V[i], Sigma.SonGrid[i], Sigma.SasMatrix[i]);
        lua_pop(L, 1);
      }
    }
  }
  /* test if operator is Hermitian */
  switch(QOperatorHermitian(*H, 0.0001*Options.Settings.zero))
  {
    case 3 :
    {
      printf("WARNING: QOperatorHermitian failed in LuaCreateSpectra, this should not happen\n");
      printf("my guess is that you are out of memory and will crash soon\n");
      break;
    }
    case 2 :
    {
      luaL_warning(L, "CreateSpectra called with an operator that is not fully Hermitian\nStrongly sugested to stop the calculation and repair this\n\n\n\nDo it\n\n\nReally Check this\n\n\nDid you?\n\n\n");
      break;
    }
    case 1 :
    {
      // all fine H is Hermitian upto zero
      break;
    }
    case 0 :
    {
      luaL_error(L, "CreateSpectra called with an operator that is not Hermitian\n");
      break;
    }
  }
  /* calculate groundstate energy */
  psiEn= (double*)malloc(Options.Npsi*Options.NT*sizeof(double));
  if(Options.E0==DBL_MAX)
  {
    for(ipsi=0;ipsi<Options.Npsi;ipsi++)
    {
      PsiHermitianOperatorPsi(*H, psi[ipsi], &psiEntemp);
      for(iT=0;iT<Options.NT;iT++)
      {
        psiEn[iT*Options.Npsi+ipsi]=psiEntemp;
      }
    }
  }
  else
  {
    for(ipsi=0;ipsi<Options.Npsi;ipsi++)
    {
      for(iT=0;iT<Options.NT;iT++)
      {
        psiEn[iT*Options.Npsi+ipsi]=Options.E0;
      }
    }
  }
  /* calculate Tpsi for all functions and transition operators 
     Tpsi is a single vector, the functions are stored at Tpsi[iT*Options.Npsi+ipsi] */
  /* different function calls with and without restrictions */
  Tpsi= (WaveFunctionType*)malloc((Options.NT*Options.Npsi)*sizeof(WaveFunctionType));
  if(Options.Settings.restrictions.N==0)
  {
    /* loop over all Transition operators */
    for(iT=0;iT<Options.NT;iT++)
    {
      /* loop over all wavefunctions */
      for(ipsi=0;ipsi<Options.Npsi;ipsi++)
      {
        if(OperatorPsi(T[iT], psi[ipsi], &Tpsi[iT*Options.Npsi+ipsi]))
        {
          luaL_error(L, "OperatorPsi failed in calculation of spectra\n Out of memory?\nTry to split the call to several calls for different wave-functions or transition operators\n");
        }
      }
    }
  }
  else
  {
    /* loop over all Transition operators */
    for(iT=0;iT<Options.NT;iT++)
    {
      /* loop over all wavefunctions */
      for(ipsi=0;ipsi<Options.Npsi;ipsi++)
      {
        if(OperatorPsiRestricted(T[iT], psi[ipsi], &Tpsi[iT*Options.Npsi+ipsi], Options.Settings.restrictions))
        {
          luaL_error(L, "OperatorPsi with restrictions failed in calculation of spectra\n Out of memory?\nTry to split the call to several calls for different wave-functions or transition operators\n");
        }
      }
    }
  }
  /* set DenseMethod to true if Hilbert space smaller then DensBoarder and calculate Basis at the same time create psiBasis */
  /* create psiBasis = sum_{i,j} T_i psi_j */
  WaveFunctionCopy(Tpsi[0], &psiDenseBasis);
  for(i=1;i<Options.NT*Options.Npsi;i++)
  {
    if((Options.DenseMethod==False)&&(psiDenseBasis.N>Options.Settings.DenseBorder))
    {
      WaveFunctionFree(&psiDenseBasis);
      break;
    }
    if(WaveFunctionAddMultipliedUnionBasis(&psiDenseBasis, Tpsi[i], (double)(1+rand())/(double)RAND_MAX))
    {
      luaL_error(L, "WaveFunctionAddMultipliedUnionBasis failed in calculation of spectra\n Out of memory?\nTry to split the call to several calls for different wave-functions or transition operators\n");
    }
  }
  /* expand psiBasis to represent H^n psiBasis for all n */
  len=0;
  if((Options.DenseMethod==True)||(psiDenseBasis.N<=Options.Settings.DenseBorder))
  {
    while(len!=psiDenseBasis.N)
    {
      len=psiDenseBasis.N;
      if((Options.Settings.restrictions.N==0)?(OperatorPsi(*H, psiDenseBasis, &OpsiDenseBasis)):(OperatorPsiRestricted(*H, psiDenseBasis, &OpsiDenseBasis, Options.Settings.restrictions)))
      {
        luaL_error(L, "OperatorPsi failed in calculation of spectra\n Out of memory?\nTry to split the call to several calls for different wave-functions or transition operators\n");
      }
      if((Options.DenseMethod==False)&&(OpsiDenseBasis.N>Options.Settings.DenseBorder))
      {
        WaveFunctionFree(&OpsiDenseBasis);
        WaveFunctionFree(&psiDenseBasis);
        break;
      }
      WaveFunctionNormalize(&OpsiDenseBasis);
      if(WaveFunctionAddUnionBasis(&psiDenseBasis, OpsiDenseBasis))
      {
        luaL_error(L, "WaveFunctionAddUnionBasis failed in calculation of spectra\n Out of memory?\nTry to split the call to several calls for different wave-functions or transition operators\n");
      }
      WaveFunctionFree(&OpsiDenseBasis);
      if((Options.DenseMethod==False)&&(psiDenseBasis.N>Options.Settings.DenseBorder))
      {
        WaveFunctionFree(&psiDenseBasis);
        break;
      }
      if(len==psiDenseBasis.N) Options.DenseMethod=True;
    }
  }
  /* If DenseMethod is true calculate spectra using dense methods */
  if(Options.DenseMethod==True)
  {
    /* create HDense */
    OperatorToDeterminantMatrix(H, &psiDenseBasis, &HDense);
    /* create TpsiDense */
    TpsiDense.Ni=Options.NT*Options.Npsi;
    TpsiDense.Nj=psiDenseBasis.N;
    TpsiDense.QComplex=psiDenseBasis.QComplex;
    if(InitCompactMatrix(&TpsiDense))
    {
      luaL_error(L, "InitCompactMatrix failed in calculation of spectra\n Out of memory?\nTry to split the call to several calls for different wave-functions or transition operators or use sparse methods ({\"DenseBorder\",0} in options)\n");
    }
    iT=0;
    for(i=0;i<Options.NT*Options.Npsi;i++)
    {
      WaveFunctionMultiply(&psiDenseBasis, 0, 0);
      WaveFunctionAddConserveBasis(&psiDenseBasis, Tpsi[i]);
      WaveFunctionFree(&Tpsi[i]);
      if(TpsiDense.QComplex==0)
      {
        for(j=1;j<=psiDenseBasis.N;j++)
        {
          TpsiDense.VRe[iT]=psiDenseBasis.ValueRe[BlokElement(j)];
          iT++;
        }
      }
      else
      {
        for(j=1;j<=psiDenseBasis.N;j++)
        {
          TpsiDense.V[iT].Re=psiDenseBasis.ValueRe[BlokElement(j)];
          TpsiDense.V[iT].Im=psiDenseBasis.ValueIm[BlokElement(j)];
          iT++;
        }
      }
    }
    free(Tpsi);
    if(Options.KrylovBasis==False)
    {
      WaveFunctionFree(&psiDenseBasis);
    }
    /* create spectra using dense methods */
    if(Options.Tensor==False)
    {
      ListOfPoles *SpecLP;
      SpecLP= (ListOfPoles*)malloc(TpsiDense.Ni*sizeof(ListOfPoles));
      if(Options.KrylovBasis==True)
      {
        HermitianCompactMatrixToListOfPolesWithEigenstates(&HDense, SpecLP, &TpsiDense, &EigenfunHDense);
      }
      else
      {
        HermitianCompactMatrixToListOfPoles(&HDense, SpecLP, &TpsiDense);
      }
      /* determine Eminpole and Emaxpole */
      Eminpole=DBL_MAX;
      Emaxpole=-DBL_MAX;
      for(i=0;i<Options.NT*Options.Npsi;i++)
      {
        if((SpecLP[i].a[0]-psiEn[i])<Eminpole)Eminpole=SpecLP[i].a[0]-psiEn[i];
        if((SpecLP[i].a[SpecLP[i].N-1]-psiEn[i])>Emaxpole)Emaxpole=SpecLP[i].a[SpecLP[i].N-1]-psiEn[i];
      }
      /* if Emin=DBL_MAX or Emax=-DBL_MAX or Gamma==0 determine these */
      if(Eminpole==Emaxpole)
      {
        if(Eminpole==0)
        {
          if(Options.Settings.Emin== DBL_MAX) Options.Settings.Emin =-1;
          if(Options.Settings.Emax==-DBL_MAX) Options.Settings.Emax = 1;
        }
        if(Eminpole<0)
        {
          if(Options.Settings.Emin== DBL_MAX) Options.Settings.Emin = 2.125*Eminpole;
          if(Options.Settings.Emax==-DBL_MAX) Options.Settings.Emax =-0.125*Eminpole;
        }
        if(Eminpole>0)
        {
          if(Options.Settings.Emin== DBL_MAX) Options.Settings.Emin =-0.125*Eminpole;
          if(Options.Settings.Emax==-DBL_MAX) Options.Settings.Emax = 2.125*Eminpole;
        }
      }
      if(Options.Settings.Emin== DBL_MAX) Options.Settings.Emin = Eminpole - 0.125 * (Emaxpole-Eminpole);
      if(Options.Settings.Emax==-DBL_MAX) Options.Settings.Emax = Emaxpole + 0.125 * (Emaxpole-Eminpole);
      if(Options.Settings.Gamma==0) Options.Settings.Gamma = 10*(Options.Settings.Emax-Options.Settings.Emin)/((double)Options.Settings.NE);
      /* transform list of poles form to spectral form and push to stack */
      Spec = (SpectraTypeArray *) lua_newuserdata(L, sizeof(SpectraTypeArray));
      luaL_setmetatable(L, "Spectra_Type");
      Spec->NSpectra=Options.NT*Options.Npsi;
      Spec->Spec= (SpectraType*)malloc(Spec->NSpectra*sizeof(SpectraType));
      Spec->Spec[0].N=Options.Settings.NE;
      Spec->Spec[0].Emin=Options.Settings.Emin;
      Spec->Spec[0].Emax=Options.Settings.Emax;
      Spec->Spec[0].dE=(Options.Settings.Emax-Options.Settings.Emin)/((double)Options.Settings.NE);
      Spec->Spec[0].EminPole=Eminpole;
      Spec->Spec[0].EmaxPole=Emaxpole;
      Spec->Spec[0].Gamma=Options.Settings.Gamma;
      Spec->Spec[0].ResonantEnergy=NAN;
      strcpy(Spec->Spec[0].Name, "Spec");
      for(i=0;i<Spec->NSpectra;i++)
      {
        SpectrumAllocate(Spec->Spec[0], &(Spec->Spec[i]));
        sprintf(Spec->Spec[i].Name, "Spectrum %i",i);
        ListOfPolesToSpectrum(psiEn[i], SpecLP[i], &(Spec->Spec[i]));
      }
      /* push a table to L containing the spectra in list of poles form */
      /* table contains {{a,b,mu=mu,type=listofpoles},{a,b,,mu=mu,type=listofpoles},{a,b,,mu=mu,type=listofpoles},....} */
      lua_createtable(L,Spec->NSpectra,0);
      for(i=0;i<Spec->NSpectra;i++)
      {
        ShiftListOfPoles(&(SpecLP[i]), -psiEn[i]);
        LuaPushTable_ListOfPoles(L, &(SpecLP[i]));
        lua_rawseti(L, -2, i+1);
      }
      /* free list of poles representation in memory */
      for(i=0;i<Options.NT*Options.Npsi;i++)
      {
        FreeListOfPoles(&(SpecLP[i]));
      }
      free(SpecLP);
      free(psiEn);
    }
    else //Dense methods with tensor True
    {
      BlockListOfPoles SpecLP;
      if(Options.KrylovBasis==True)
      {
        HermitianCompactMatrixToBlockListOfPolesWithEigenstates(&HDense, &SpecLP, &TpsiDense, &EigenfunHDense);
      }
      else
      {
        HermitianCompactMatrixToBlockListOfPoles(&HDense, &SpecLP, &TpsiDense);
      }
      /* determine Eminpole and Emaxpole */
      if(SpecLP.QComplex==0)
      {
        Eminpole=SpecLP.VRe[0]-psiEn[0];
        Emaxpole=SpecLP.VRe[(1+SpecLP.BlockSize*SpecLP.BlockSize)*(SpecLP.N-1)]-psiEn[0];
      }
      else
      {
        Eminpole=SpecLP.V[0].Re-psiEn[0];
        Emaxpole=SpecLP.V[(1+SpecLP.BlockSize*SpecLP.BlockSize)*(SpecLP.N-1)].Re-psiEn[0];
      }
      /* if Emin=DBL_MAX or Emax=-DBL_MAX or Gamma==0 determine these */
      if(Eminpole==Emaxpole)
      {
        if(Eminpole==0)
        {
          if(Options.Settings.Emin== DBL_MAX) Options.Settings.Emin =-1;
          if(Options.Settings.Emax==-DBL_MAX) Options.Settings.Emax = 1;
        }
        if(Eminpole<0)
        {
          if(Options.Settings.Emin== DBL_MAX) Options.Settings.Emin = 2.125*Eminpole;
          if(Options.Settings.Emax==-DBL_MAX) Options.Settings.Emax =-0.125*Eminpole;
        }
        if(Eminpole>0)
        {
          if(Options.Settings.Emin== DBL_MAX) Options.Settings.Emin =-0.125*Eminpole;
          if(Options.Settings.Emax==-DBL_MAX) Options.Settings.Emax = 2.125*Eminpole;
        }
      }
      if(Options.Settings.Emin== DBL_MAX) Options.Settings.Emin = Eminpole - 0.125 * (Emaxpole-Eminpole);
      if(Options.Settings.Emax==-DBL_MAX) Options.Settings.Emax = Emaxpole + 0.125 * (Emaxpole-Eminpole);
      if(Options.Settings.Gamma==0) Options.Settings.Gamma = 10*(Options.Settings.Emax-Options.Settings.Emin)/((double)Options.Settings.NE);
      /* transform list of poles form to spectral form and push to stack */
      Spec = (SpectraTypeArray *) lua_newuserdata(L, sizeof(SpectraTypeArray));
      luaL_setmetatable(L, "Spectra_Type");
      Spec->NSpectra=Sqr(Options.NT*Options.Npsi);
      Spec->Spec= (SpectraType*)malloc(Spec->NSpectra*sizeof(SpectraType));
      BlockListOfPolesToSpectraWithBroadening(SpecLP, Spec->Spec, Options.Settings.Emin, Options.Settings.Emax, Options.Settings.NE, Options.Settings.Gamma);
      for(i=0;i<Spec->NSpectra;i++)
      {
        Spec->Spec[i].EminPole=Eminpole;
        Spec->Spec[i].EmaxPole=Emaxpole;
      }
      /* push a table to L containing the spectra in Block list of poles form */
      ShiftBlockListOfPoles(&SpecLP, -psiEn[0]);
      LuaPushTable_BlockListOfPoles(L, &SpecLP);
      /* free list of poles representation in memory */
      FreeBlockListOfPoles(&SpecLP);
      free(psiEn);
    }
    /* create WaveFunction Basis for output if needed and put this on stack and return 2 or 3*/
    if(Options.KrylovBasis==True)
    {
      WaveFunctionType *buff;
      WaveFunctionMultiply(&psiDenseBasis, 0, 0);
      if((psiDenseBasis.QComplex==0)&&(EigenfunHDense.QComplex!=0)) MakeComplexWaveFunctionFromReal(&psiDenseBasis, psiDenseBasis);
      if((psiDenseBasis.QComplex!=0)&&(EigenfunHDense.QComplex==0)) MakeRealWaveFunctionFromComplex(psiDenseBasis, &psiDenseBasis);
      k=0;
      lua_createtable(L, EigenfunHDense.Ni, 0);
      if(EigenfunHDense.QComplex==0)
      {
        for(i=0;i<EigenfunHDense.Ni;i++)
        {
          for(j=1;j<=EigenfunHDense.Nj;j++)
          {
            psiDenseBasis.ValueRe[BlokElement(j)]=EigenfunHDense.VRe[k];
            k++;
          }
          buff = (WaveFunctionType*) lua_newuserdata(L, sizeof(WaveFunctionType));
          luaL_setmetatable(L, "Wavefunction_Type");
          WaveFunctionCopy(psiDenseBasis, buff);
          lua_rawseti(L, -2, i+1);
        }
      }
      else
      {
        for(i=0;i<EigenfunHDense.Ni;i++)
        {
          for(j=1;j<=EigenfunHDense.Nj;j++)
          {
            psiDenseBasis.ValueRe[BlokElement(j)]=EigenfunHDense.V[k].Re;
            psiDenseBasis.ValueIm[BlokElement(j)]=EigenfunHDense.V[k].Im;
            k++;
          }
          buff = (WaveFunctionType*) lua_newuserdata(L, sizeof(WaveFunctionType));
          luaL_setmetatable(L, "Wavefunction_Type");
          WaveFunctionCopy(psiDenseBasis, buff);
          lua_rawseti(L, -2, i+1);
        }
      }
      WaveFunctionFree(&psiDenseBasis);
      FreeCompactMatrix(&EigenfunHDense);
      return 3;
    }
    else
    {
      return 2;
    }
  }
  else
  {
    /* start calculating the spectra using Lanczos */
    if(Options.Tensor==False)
    {
      TriDiagonalMatrixType *GTri;
      /* calculate spectra in Tridiagonal form for each combination of operator and wave-function seperately*/
      GTri = (TriDiagonalMatrixType*)malloc(Options.NT*Options.Npsi*sizeof(TriDiagonalMatrixType));
      for(i=0;i<Options.NT*Options.Npsi;i++)
      {
        GTri[i].N=Options.Settings.NTri;
        InitTriDiagonalMatrix(&GTri[i]);
        GTri[i].mu=0; // we have no knowlede on the chemical potential at this point
      }
      /* alocate memory for the krylov basis */
      if(Options.KrylovBasis==True)
      {
        psiBasis= (WaveFunctionType**)malloc(Options.NT*Options.Npsi*sizeof(WaveFunctionType*));
        for(i=0;i<Options.NT*Options.Npsi;i++)
        {
          psiBasis[i]= (WaveFunctionType*)malloc(Options.Settings.NTri*sizeof(WaveFunctionType));
        }
      }
      /* loop over all Transition operators */
      for(iT=0;iT<Options.NT;iT++)
      {
        /* loop over all wavefunctions */
        for(ipsi=0;ipsi<Options.Npsi;ipsi++)
        {
          if(Options.Settings.restrictions.N==0)
          {
            if(Options.KrylovBasis==True)
            {
              LanczosTriDiagonalizeKrylov(*H, Tpsi[iT*Options.Npsi+ipsi], Options.Settings.epsilon, &GTri[iT*Options.Npsi+ipsi], psiBasis[iT*Options.Npsi+ipsi]);
            }
            else
            {
              LanczosTriDiagonalize(*H, Tpsi[iT*Options.Npsi+ipsi], Options.Settings.epsilon, &GTri[iT*Options.Npsi+ipsi], Max(Options.Settings.Gamma,Options.Settings.epsilon));
            }
          }
          else
          {
            if(Options.KrylovBasis==True)
            {
              LanczosTriDiagonalizeKrylovRestricted(*H, Tpsi[iT*Options.Npsi+ipsi],  Options.Settings.restrictions, Options.Settings.epsilon, &GTri[iT*Options.Npsi+ipsi], psiBasis[iT*Options.Npsi+ipsi]);
            }
            else
            {
              LanczosTriDiagonalizeRestricted(*H, Tpsi[iT*Options.Npsi+ipsi], Options.Settings.restrictions, Options.Settings.epsilon, &GTri[iT*Options.Npsi+ipsi], Max(Options.Settings.Gamma,Options.Settings.epsilon));
            }
          }
          if(Options.stream!=NULL)
          {
            if (Options.NT*Options.Npsi==1)
            {
              snprintf(filename, MAXStrLength, "%s", Options.stream);
            }
            else
            {
              snprintf(filename, MAXStrLength, "%s_%02d", Options.stream, iT*Options.Npsi+ipsi);
            }
            FilePrintPlainTriDiagonalMatrix(GTri[iT*Options.Npsi+ipsi], filename, "w");
          }
        }
      }
      /* determine Eminpole and Emaxpole */
      /* get eigenvalues of all Tri diagonal matrices and set smallest and largest eigenvalue */
      Eminpole=DBL_MAX;
      Emaxpole=-DBL_MAX;
      val= (double*)malloc((Options.Settings.NTri)*sizeof(double));
      for(i=0;i<Options.NT*Options.Npsi;i++)
      {
        TriDiagonalMatrixToEigenvalues(GTri[i], val);
        if((val[0]-psiEn[i])<Eminpole)Eminpole=val[0]-psiEn[i];
        if((val[Options.Settings.NTri-1]-psiEn[i])>Emaxpole)Emaxpole=val[Options.Settings.NTri-1]-psiEn[i];
      }
      /* if Emin=DBL_MAX or Emax=-DBL_MAX or Gamma==0 determine these */
      if(Eminpole==Emaxpole)
      {
        if(Eminpole==0)
        {
          if(Options.Settings.Emin== DBL_MAX) Options.Settings.Emin =-1;
          if(Options.Settings.Emax==-DBL_MAX) Options.Settings.Emax = 1;
        }
        if(Eminpole<0)
        {
          if(Options.Settings.Emin== DBL_MAX) Options.Settings.Emin = 2.125*Eminpole;
          if(Options.Settings.Emax==-DBL_MAX) Options.Settings.Emax =-0.125*Eminpole;
        }
        if(Eminpole>0)
        {
          if(Options.Settings.Emin== DBL_MAX) Options.Settings.Emin =-0.125*Eminpole;
          if(Options.Settings.Emax==-DBL_MAX) Options.Settings.Emax = 2.125*Eminpole;
        }
      }
      if(Options.Settings.Emin== DBL_MAX) Options.Settings.Emin = Eminpole - 0.125 * (Emaxpole-Eminpole);
      if(Options.Settings.Emax==-DBL_MAX) Options.Settings.Emax = Emaxpole + 0.125 * (Emaxpole-Eminpole);
      if(Options.Settings.Gamma==0) Options.Settings.Gamma = 10*(Options.Settings.Emax-Options.Settings.Emin)/((double)Options.Settings.NE);
      // TODO different tridiag to spectral form when self energy is present
      
      /* transform tridiagonal form to spectral form and push to stack */
      Spec = (SpectraTypeArray *) lua_newuserdata(L, sizeof(SpectraTypeArray));
      luaL_setmetatable(L, "Spectra_Type");
      Spec->NSpectra=Options.NT*Options.Npsi;
      Spec->Spec= (SpectraType*)malloc(Spec->NSpectra*sizeof(SpectraType));
      Spec->Spec[0].N=Options.Settings.NE;
      Spec->Spec[0].Emin=Options.Settings.Emin;
      Spec->Spec[0].Emax=Options.Settings.Emax;
      Spec->Spec[0].dE=(Options.Settings.Emax-Options.Settings.Emin)/((double)Options.Settings.NE);
      Spec->Spec[0].EminPole=Eminpole;
      Spec->Spec[0].EmaxPole=Emaxpole;
      Spec->Spec[0].Gamma=Options.Settings.Gamma;
      Spec->Spec[0].ResonantEnergy=NAN;
      strcpy(Spec->Spec[0].Name, "Spec");
      for(i=0;i<Spec->NSpectra;i++)
      {
        SpectrumAllocate(Spec->Spec[0], &(Spec->Spec[i]));
        sprintf(Spec->Spec[i].Name, "Spectrum %i",i);
        TriDiagonalMatrixToSpectrum(psiEn[i], GTri[i], &(Spec->Spec[i]));
      }
      //TODO different spectrum pushed to stack when self energy is present
      
      /* push a table to L containing the spectra in tridiagonal form */
      /* table contains {{a,b,mu=mu,type=tri},{a,b,,mu=mu,type=tri},{a,b,,mu=mu,type=tri},....} */
      lua_createtable(L,Spec->NSpectra,0);
      for(i=0;i<Spec->NSpectra;i++)
      {
        ShiftTriDiagonalMatrix(&(GTri[i]), -psiEn[i]);
        LuaPushTable_TriDiagonalMatrixType(L, &(GTri[i]));
        lua_rawseti(L, -2, i+1);
      }
      /* free tridiagonal representation in memory */
      for(i=0;i<Spec->NSpectra;i++)
      {
        FreeTriDiagonalMatrix(&GTri[i]);
      }
      free(GTri);
      free(psiEn);
      free(Tpsi);
      /* return KrylovBasis if asked for and number of elements pushed on stack */
      if(Options.KrylovBasis==True)
      {
        if(Options.NT*Options.Npsi==1)
        {
          LuaPushTable_WaveFunctionType_RemoveEmpty(L, psiBasis[0], Options.Settings.NTri);
          free(psiBasis[0]);
        }
        else
        {
          lua_createtable(L,Options.NT*Options.Npsi,0);
          for(i=0;i<Options.NT*Options.Npsi;i++)
          {
            LuaPushTable_WaveFunctionType_RemoveEmpty(L, psiBasis[i], Options.Settings.NTri);
            free(psiBasis[i]);
            lua_rawseti(L, -2, i+1);
          }
        }
        free(psiBasis);
        return 3;
      }
      else
      {
        return 2;
      }
    }
    else
    {
      /* calcualte spectra using tensor form, i.e. returned are Sqr(Options.NT * Options.Npsi) sepctra */
      WaveFunctionType TpsiSameBasis[Options.NT*Options.Npsi];
      BlockTriDiagonalMatrixType GTri;
      /* give warning when more than one function defined and E0 not set */
      if((Options.Npsi!=1)&&(Options.E0==DBL_MAX))
      {
        luaL_warning(L, "CreateSpectra called with option Tensor->true and %d wave-functions but no energy-shift set.\nTensors can have only one shift and no different shift for each eigen-state.\nConsider using E0->0.0 in order to define a fixed shift.\nUsed the energy of the first wave-function (%f) in this case.\n(Probably does not make physical sence!)",Options.Npsi,psiEn[0]);
      }
      /* place starting functions on same basis */
      WaveFunctionCopy(Tpsi[0], &TpsiSameBasis[0]);
      WaveFunctionFree(&(Tpsi[0]));
      for(i=1;i<Options.NT*Options.Npsi;i++)
      {
        if(WaveFunctionAddMultipliedUnionBasis(&TpsiSameBasis[0], Tpsi[i], 0)) return luaL_error(L, "malloc failed. Out of memory during initialization in CalculateSpectra\n");
      }
      for(i=1;i<Options.NT*Options.Npsi;i++)
      {
        if(WaveFunctionInitCopyBasis(&(TpsiSameBasis[i]), TpsiSameBasis[0])) return luaL_error(L, "malloc failed. Out of memory during initialization in CalculateSpectra\n");
        WaveFunctionAddConserveBasis(&(TpsiSameBasis[i]), Tpsi[i]);
        WaveFunctionFree(&(Tpsi[i]));
      }
      /* init BlockTriDiagonal Matrix */
      GTri.BlockSize = Options.NT*Options.Npsi;
      GTri.N = Options.Settings.NTri;
      /* create block tridiagonal representation */
      GTri.QComplex=0;
      if(H->QComplex) GTri.QComplex=1;
      for(i=0;i<GTri.BlockSize;i++)
      {
        if(TpsiSameBasis[i].QComplex) GTri.QComplex=1;
      }
      InitBlockTridiagonalMatrix(&GTri);
      GTri.mu=psiEn[0];
      
      if(Options.Settings.restrictions.N==0)
      {
        if(Options.KrylovBasis==True)
        {
          BlockSizeN= (uint32_t*)malloc(GTri.N*sizeof(uint32_t));
          psiBasis= (WaveFunctionType**)malloc(Options.Settings.NTri*sizeof(WaveFunctionType*));
          LanczosBlockTriDiagonalizeKrylov(*H, TpsiSameBasis, Options.Settings.epsilon, &GTri, psiBasis, BlockSizeN);
        }
        else
        {
          LanczosBlockTriDiagonalize(*H, TpsiSameBasis, Options.Settings.epsilon, &GTri, Max(Options.Settings.Gamma,Options.Settings.epsilon));
        }
      }
      else
      {
        if(Options.KrylovBasis==True)
        {
          BlockSizeN= (uint32_t*)malloc(GTri.N*sizeof(uint32_t));
          psiBasis= (WaveFunctionType**)malloc(Options.Settings.NTri*sizeof(WaveFunctionType*));
          LanczosBlockTriDiagonalizeKrylovRestricted(*H, TpsiSameBasis, Options.Settings.restrictions, Options.Settings.epsilon, &GTri, psiBasis, BlockSizeN);
        }
        else
        {
          LanczosBlockTriDiagonalizeRestricted(*H, TpsiSameBasis, Options.Settings.restrictions, Options.Settings.epsilon, &GTri, Max(Options.Settings.Gamma,Options.Settings.epsilon));
        }
      }
      /* determine Eminpole and Emaxpole */
      /* get eigenvalues Block Tridiagonal matrix and set smallest and largest eigenvalue */
      Eminpole = DBL_MAX;
      Emaxpole =-DBL_MAX;
      //set Eminpole and Emaxpole (they are in the output and I need them for auto range determination
      val= (double*)malloc(GTri.N*GTri.BlockSize*sizeof(double));
      BlockTriDiagonalMatrixEigenValues(&GTri,val);
      Eminpole=val[0]-psiEn[0];
      Emaxpole=val[GTri.N*GTri.BlockSize-1]-psiEn[0];
      free(val);
      /* if Emin=DBL_MAX or Emax=-DBL_MAX or Gamma==0 determine these */
      if(Eminpole==Emaxpole)
      {
        if(Eminpole==0)
        {
          if(Options.Settings.Emin== DBL_MAX) Options.Settings.Emin =-1;
          if(Options.Settings.Emax==-DBL_MAX) Options.Settings.Emax = 1;
        }
        if(Eminpole<0)
        {
          if(Options.Settings.Emin== DBL_MAX) Options.Settings.Emin = 2.125*Eminpole;
          if(Options.Settings.Emax==-DBL_MAX) Options.Settings.Emax =-0.125*Eminpole;
        }
        if(Eminpole>0)
        {
          if(Options.Settings.Emin== DBL_MAX) Options.Settings.Emin =-0.125*Eminpole;
          if(Options.Settings.Emax==-DBL_MAX) Options.Settings.Emax = 2.125*Eminpole;
        }
      }
      if(Options.Settings.Emin== DBL_MAX) Options.Settings.Emin = Eminpole - 0.125 * (Emaxpole-Eminpole);
      if(Options.Settings.Emax==-DBL_MAX) Options.Settings.Emax = Emaxpole + 0.125 * (Emaxpole-Eminpole);
      if(Options.Settings.Gamma==0) Options.Settings.Gamma = 10*(Options.Settings.Emax-Options.Settings.Emin)/((double)Options.Settings.NE);
      
      //TODO: different transformation to spectralform when self energy is present
      
      /* transform tridiagonal form to spectral form and push to stack */
      Spec = (SpectraTypeArray *) lua_newuserdata(L, sizeof(SpectraTypeArray));
      luaL_setmetatable(L, "Spectra_Type");
      Spec->NSpectra=Sqr(Options.NT*Options.Npsi);
      Spec->Spec= (SpectraType*)malloc(Spec->NSpectra*sizeof(SpectraType));
      BlockTriDiagonalMatrixToSpectra(psiEn[0], &GTri, Spec->Spec, Options.Settings.Emin, Options.Settings.Emax, Options.Settings.NE, Options.Settings.Gamma);
      for(i=0;i<Spec->NSpectra;i++)
      {
        Spec->Spec[i].EminPole=Eminpole;
        Spec->Spec[i].EmaxPole=Emaxpole;
      }
      /* push a table to L containing the spectra in Block Tridiagonal form */
      ShiftBlockTridiagonalMatrix(&GTri, -psiEn[0]);
      LuaPushTable_BlockTriDiagonalMatrixType(L, &GTri);
      /* push Krylov basis to stack */
      if(Options.KrylovBasis==True)
      {
        /* push a table to L containing the krylov basis */
        lua_createtable(L,GTri.N,0);
        {
          for(i=0;i<GTri.N;i++)
          {
            if(BlockSizeN[i]==0) break;
            WaveFunctionArraySeperateBasis(psiBasis[i], BlockSizeN[i]);
            LuaPushTable_WaveFunctionType(L, psiBasis[i], BlockSizeN[i]);
            lua_rawseti(L, -2, i+1);
          }
        }
        free(BlockSizeN);
      }
      /* free tridiagonal representation in memory */
      FreeBlockTridiagonalMatrix(&GTri);
      
      if(Options.KrylovBasis==True)
      {
        return 3;
      }
      else
      {
        return 2;
      }
    }
  }
}



/*  calculates <psi_i | T_1^+ G_1^+ T_2^+ G_2 T_2 G_1 T_1 | psi_i >
 *  On stack:
 *     H1
 *     H2
 *     T1  - allowed to be a table
 *     T2  - allowed to be a table
 *     psi - allowed to be a table
 *     options
 *        {E1min, E1max, N1, Gamma1}
 *        {E2min, E2max, N2, Gamma2}
 *        NTri1
 *        NTri2
 *        epsilon
 *        restrictions1
 *        restrictions2
 *        E0 E0_1 E0_2
 */
static int LuaCreateResonantSpectra(lua_State *L)
{
  struct
  {
    ConvergenceSettingsType Settings;
    double E0, E01, E02;
    uint32_t Npsi;
    uint32_t NT1;
    uint32_t NT2;
    Bool Tensor;
    Bool DenseMethod;
  } Options;
  Options.Settings=InitConvergenceSettings();
  Options.E0=DBL_MAX;
  Options.E01=DBL_MAX;
  Options.E02=DBL_MAX;
  Options.Tensor=False;
  Options.DenseMethod=False;
  ConvergenceSettingsType SettingsNow;
  SpectraTypeArray *Spec;
  SpectraType *SpecT1psi;
  TriDiagonalMatrixType *SpecTri;
  OperatorType *H1,*H2,*T1, *T2;
  WaveFunctionType *psi, HTpsiDenseBasis,HTHTpsiDenseBasis,Opsi,TpsiSum;
  uint32_t iT,ipsi,len,i;
  LuaOptionsType LuaCreateResonantSpectraOptions[]={
    {"NTri1","uint32_t",&(Options.Settings.NTri1)},
    {"NTri2","uint32_t",&(Options.Settings.NTri2)},
    {"Tensor","bool",&(Options.Tensor)},
    {"restrictions1", "DeterminantRestrictionType", &(Options.Settings.restrictions1)},
    {"restrictions2", "DeterminantRestrictionType", &(Options.Settings.restrictions2)},
    {"Restrictions1", "DeterminantRestrictionType", &(Options.Settings.restrictions1)},
    {"Restrictions2", "DeterminantRestrictionType", &(Options.Settings.restrictions2)},
    {"epsilon","double",&(Options.Settings.epsilon)},
    {"Epsilon","double",&(Options.Settings.epsilon)},
    {"E0","double",&(Options.E0)},
    {"E0_1","double",&(Options.E01)},
    {"E0_2","double",&(Options.E02)},
    {"Emin1","double",&(Options.Settings.Emin1)},
    {"Emax1","double",&(Options.Settings.Emax1)},
    {"NE1","uint32_t",&(Options.Settings.NE1)},
    {"Gamma1","double",&(Options.Settings.Gamma1)},
    {"Emin2","double",&(Options.Settings.Emin2)},
    {"Emax2","double",&(Options.Settings.Emax2)},
    {"NE2","uint32_t",&(Options.Settings.NE2)},
    {"Gamma2","double",&(Options.Settings.Gamma2)},
    {"DenseBorder","uint32_t",&(Options.Settings.DenseBorder)},
    {"ForceDenseMethod","uint32_t",&(Options.DenseMethod)},
    {NULL, NULL, NULL}};
  /* free garbage full cycle*/
  if(L->l_G->gcrunning) lua_gc(L, LUA_GCCOLLECT, 0);
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
    if(Options.NT1==0) return luaL_error(L, "Errof CreateResonantSpectra called with a table of transition operators (1) of length 0\n");
  }
  else
  {
    Options.NT1 = 1;
    T1 = (OperatorType *) luaL_checkudata(L, 3, "Operator_Type");
  }
  if(lua_type(L,4)==LUA_TTABLE)
  {
    T2=luaL_checkarray_OperatorType(L, 4, &Options.NT2);
    if(Options.NT2==0) return luaL_error(L, "Errof CreateResonantSpectra called with a table of transition operators (2) of length 0\n");
  }
  else
  {
    Options.NT2 = 1;
    T2 = (OperatorType *) luaL_checkudata(L, 4, "Operator_Type");
  }
  if(lua_type(L,5)==LUA_TTABLE)
  {
    psi=luaL_checkarray_WaveFunctionType(L, 5, &Options.Npsi);
    if(Options.Npsi==0) return luaL_error(L, "Errof CreateResonantSpectra called with a table of wave functions of length 0\n");
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
  if(Options.E01==DBL_MAX)Options.E01=Options.E0;
  if(Options.E02==DBL_MAX)Options.E02=Options.E0;
  
  /* test if operator is Hermitian */
  switch(QOperatorHermitian(*H1, 0.0001*Options.Settings.zero))
  {
    case 3 :
    {
      printf("WARNING: QOperatorHermitian failed in LuaCreateResonantSpectra, this should not happen\n");
      printf("my guess is that you are out of memory and will crash soon\n");
      break;
    }
    case 2 :
    {
      luaL_warning(L, "CreateResonantSpectra called with an operator (H1) that is not fully Hermitian\nStrongly sugested to stop the calculation and repair this\n\n\n\nDo it\n\n\nReally Check this\n\n\nDid you?\n\n\n");
      break;
    }
    case 1 :
    {
      // all fine H is Hermitian upto zero
      break;
    }
    case 0 :
    {
      luaL_error(L, "CreateResonantSpectra called with an operator (H1) that is not Hermitian\n");
      break;
    }
  }
  switch(QOperatorHermitian(*H2, 0.0001*Options.Settings.zero))
  {
    case 3 :
    {
      printf("WARNING: QOperatorHermitian failed in LuaCreateResonantSpectra, this should not happen\n");
      printf("my guess is that you are out of memory and will crash soon\n");
      break;
    }
    case 2 :
    {
      luaL_warning(L, "CreateResonantSpectra called with an operator (H2) that is not fully Hermitian\nStrongly sugested to stop the calculation and repair this\n\n\n\nDo it\n\n\nReally Check this\n\n\nDid you?\n\n\n");
      break;
    }
    case 1 :
    {
      // all fine H is Hermitian upto zero
      break;
    }
    case 0 :
    {
      luaL_error(L, "CreateResonantSpectra called with an operator (H2) that is not Hermitian\n");
      break;
    }
  }
  /* allocate spectra and a pointer to return matrix on lua stack */
  Spec = (SpectraTypeArray *) lua_newuserdata(L, sizeof(SpectraTypeArray));
  luaL_setmetatable(L, "Spectra_Type");
  Spec->NSpectra=0;
  Spec->Spec=NULL;
  /* test if dense methods can be used */
  // at the same time create the wave-functions that define the basis sets of the dense matrices used in the dense methods
  // we need the vectors:
  // HmDenseBasis       : on which the Matrix for the intermediate state can be represented
  // HfDenseBasis       : on which the Matrix fo rthe final states can be represented
  {
    // psi
    if (Options.DenseMethod==False) if(psi[0].N>Options.Settings.DenseBorder) goto EndOfUseDenseTest;
    WaveFunctionCopy(psi[0], &HTpsiDenseBasis);
    for(i=1;i<Options.Npsi;i++)
    {
      WaveFunctionAddUnionBasis(&HTpsiDenseBasis, psi[i]);
      if (Options.DenseMethod==False) if(HTpsiDenseBasis.N>Options.Settings.DenseBorder)
      {
        WaveFunctionFree(&HTpsiDenseBasis);
        goto EndOfUseDenseTest;
      }
    }
    // T psi
    if((Options.Settings.restrictions1.N==0)?(OperatorPsi(T1[0], HTpsiDenseBasis, &TpsiSum)):(OperatorPsiRestricted(T1[0], HTpsiDenseBasis, &TpsiSum, Options.Settings.restrictions1)))
    {
      WaveFunctionFree(&HTpsiDenseBasis);
      goto EndOfUseDenseTest;
    }
    for(i=1;i<Options.NT1;i++)
    {
      if((Options.Settings.restrictions1.N==0)?(OperatorPsi(T1[i], HTpsiDenseBasis, &Opsi)):(OperatorPsiRestricted(T1[i], HTpsiDenseBasis, &Opsi, Options.Settings.restrictions1)))
      {
        WaveFunctionFree(&TpsiSum);
        WaveFunctionFree(&HTpsiDenseBasis);
        goto EndOfUseDenseTest;
      }
      WaveFunctionAddUnionBasis(&TpsiSum, Opsi);
      WaveFunctionFree(&Opsi);
      if (Options.DenseMethod==False) if(TpsiSum.N>Options.Settings.DenseBorder)
      {
        WaveFunctionFree(&TpsiSum);
        WaveFunctionFree(&HTpsiDenseBasis);
        goto EndOfUseDenseTest;
      }
    }
    WaveFunctionFree(&HTpsiDenseBasis);
    HTpsiDenseBasis=TpsiSum;
    WaveFunctionNormalize(&HTpsiDenseBasis);
    // H T psi
    len=0;
    while(len!=HTpsiDenseBasis.N)
    {
      len=HTpsiDenseBasis.N;
      if((Options.Settings.restrictions1.N==0)?(OperatorPsi(*H1, HTpsiDenseBasis, &Opsi)):(OperatorPsiRestricted(*H1, HTpsiDenseBasis, &Opsi, Options.Settings.restrictions1)))
      {
        WaveFunctionFree(&HTpsiDenseBasis);
        goto EndOfUseDenseTest;
      }
      WaveFunctionFree(&HTpsiDenseBasis);
      if (Options.DenseMethod==False) if(Opsi.N>Options.Settings.DenseBorder)
      {
        WaveFunctionFree(&Opsi);
        goto EndOfUseDenseTest;
      }
      HTpsiDenseBasis=Opsi;
      WaveFunctionNormalize(&HTpsiDenseBasis);
    }
    // T H T psi
    WaveFunctionCopy(HTpsiDenseBasis, &HTHTpsiDenseBasis);
    if((Options.Settings.restrictions2.N==0)?(OperatorPsi(T2[0], HTHTpsiDenseBasis, &TpsiSum)):(OperatorPsiRestricted(T2[0], HTHTpsiDenseBasis, &TpsiSum, Options.Settings.restrictions2)))
    {
      WaveFunctionFree(&HTHTpsiDenseBasis);
      WaveFunctionFree(&HTpsiDenseBasis);
      goto EndOfUseDenseTest;
    }
    for(i=1;i<Options.NT2;i++)
    {
      if((Options.Settings.restrictions2.N==0)?(OperatorPsi(T2[i], HTHTpsiDenseBasis, &Opsi)):(OperatorPsiRestricted(T2[i], HTHTpsiDenseBasis, &Opsi, Options.Settings.restrictions2)))
      {
        WaveFunctionFree(&TpsiSum);
        WaveFunctionFree(&HTHTpsiDenseBasis);
        WaveFunctionFree(&HTpsiDenseBasis);
        goto EndOfUseDenseTest;
      }
      WaveFunctionAddUnionBasis(&TpsiSum, Opsi);
      WaveFunctionFree(&Opsi);
      if (Options.DenseMethod==False) if(TpsiSum.N>Options.Settings.DenseBorder)
      {
        WaveFunctionFree(&TpsiSum);
        WaveFunctionFree(&HTHTpsiDenseBasis);
        WaveFunctionFree(&HTpsiDenseBasis);
        goto EndOfUseDenseTest;
      }
    }
    WaveFunctionFree(&HTHTpsiDenseBasis);
    HTHTpsiDenseBasis=TpsiSum;
    WaveFunctionNormalize(&HTHTpsiDenseBasis);
    // H T H T psi
    len=0;
    while(len!=HTHTpsiDenseBasis.N)
    {
      len=HTHTpsiDenseBasis.N;
      if((Options.Settings.restrictions2.N==0)?(OperatorPsi(*H2, HTHTpsiDenseBasis, &Opsi)):(OperatorPsiRestricted(*H2, HTHTpsiDenseBasis, &Opsi, Options.Settings.restrictions2)))
      {
        WaveFunctionFree(&HTHTpsiDenseBasis);
        WaveFunctionFree(&HTpsiDenseBasis);
        goto EndOfUseDenseTest;
      }
      WaveFunctionFree(&HTHTpsiDenseBasis);
      if (Options.DenseMethod==False) if(Opsi.N>Options.Settings.DenseBorder)
      {
        WaveFunctionFree(&Opsi);
        WaveFunctionFree(&HTpsiDenseBasis);
        goto EndOfUseDenseTest;
      }
      WaveFunctionNormalize(&Opsi);
      HTHTpsiDenseBasis=Opsi;
    }
    // finsihed all tests and psiDenseBasis still fits
    Options.DenseMethod=True;
  EndOfUseDenseTest:;
  }
  if(Options.DenseMethod==True)
  {
    if(Options.Tensor==False)
    {
      ResonantSpectraDense(&HTpsiDenseBasis,&HTHTpsiDenseBasis,*H1, *H2, T1, Options.NT1, T2, Options.NT2, psi, Options.Npsi, &(Spec->Spec), Options.E01, Options.E02, &(Options.Settings));
      Spec->NSpectra=Options.NT1*Options.NT2*Options.Npsi*(Options.Settings.NE1+1);
    }
    else
    {
      return luaL_error(L,"Tensor form for Resonant spectra not yet available");
    }
  }
  else
  {
    // use sparse methods
    if(Options.Tensor==False)
    {
      for(ipsi=0;ipsi<Options.Npsi;ipsi++)
      {
        for(iT=0;iT<Options.NT1;iT++)
        {
          SettingsNow = Options.Settings;
          ResonantSpectra(*H1, *H2, T1[iT], T2, Options.NT2, psi[ipsi], &SpecTri, &SpecT1psi, Options.E01, Options.E02, &(SettingsNow));
          Spec->NSpectra+=(SettingsNow.NE1+1)*Options.NT2;
          Spec->Spec= (SpectraType*)realloc(Spec->Spec, Spec->NSpectra*sizeof(SpectraType));
          for(i=0;i<(SettingsNow.NE1+1)*Options.NT2;i++)
          {
            FreeTriDiagonalMatrix(&(SpecTri[i]));
          }
          free(SpecTri);
          for(i=0;i<(SettingsNow.NE1+1)*Options.NT2;i++)
          {
            sprintf(SpecT1psi[i].Name, "%s T1 %i/%i psi %i/%i",SpecT1psi[i].Name,iT,Options.NT1,ipsi,Options.Npsi);
            Spec->Spec[i+Spec->NSpectra-(SettingsNow.NE1+1)*Options.NT2] = SpecT1psi[i];
          }
          free(SpecT1psi);
        }
      }
    }
    else
    {
      
      //ResonantSpectraTensor(*H1, *H2, T1, Options.NT1, T2, Options.NT2, psi, Options.Npsi, &SpecTri, &(Spec->Spec[0]), Options.E01, Options.E02, &Options.Settings);
      return luaL_error(L,"Tensor form for Resonant spectra not yet available");
    }
  }
  return 1;
}





static int LuaCreateFluorescenceYield(lua_State *L)
{
  struct
  {
    ConvergenceSettingsType Settings;
    uint32_t Npsi;
    uint32_t NT1;
    uint32_t NT2;
  } Options;
  Options.Settings=InitConvergenceSettings();
  SpectraTypeArray *Spec;
  OperatorType *H,*T1, *T2;
  WaveFunctionType *psi;
  uint32_t iT,ipsi,len;
  double *psiEn;
  LuaOptionsType LuaCreateSpectraOptions[]={
    {"NTri","uint32_t",&(Options.Settings.NTri)},
    {"epsilon","double",&(Options.Settings.epsilon)},
    {"Epsilon","double",&(Options.Settings.epsilon)},
    {"restrictions", "DeterminantRestrictionType", &(Options.Settings.restrictions)},
    {"Restrictions", "DeterminantRestrictionType", &(Options.Settings.restrictions)},
    {"Emin","double",&(Options.Settings.Emin)},
    {"Emax","double",&(Options.Settings.Emax)},
    {"NE","uint32_t",&(Options.Settings.NE)},
    {"Gamma","double",&(Options.Settings.Gamma)},
    {NULL, NULL, NULL}};
  /* free garbage full cycle*/
  if(L->l_G->gcrunning) lua_gc(L, LUA_GCCOLLECT, 0);
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
    if(Options.NT1==0) return luaL_error(L, "Errof LuaCreateFluorescenceYield called with a table of transition operators (1) of length 0\n");
  }
  else
  {
    Options.NT1 = 1;
    T1 = (OperatorType *) luaL_checkudata(L, 2, "Operator_Type");
  }
  if(lua_type(L,3)==LUA_TTABLE)
  {
    T2=luaL_checkarray_OperatorType(L, 3, &Options.NT2);
    if(Options.NT2==0) return luaL_error(L, "Errof LuaCreateFluorescenceYield called with a table of transition operators (2) of length 0\n");
  }
  else
  {
    Options.NT2 = 1;
    T2 = (OperatorType *) luaL_checkudata(L, 3, "Operator_Type");
  }
  if(lua_type(L,4)==LUA_TTABLE)
  {
    psi=luaL_checkarray_WaveFunctionType(L, 4, &Options.Npsi);
    if(Options.Npsi==0) return luaL_error(L, "Errof LuaCreateFluorescenceYield called with a table of wave functions of length 0\n");
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
  }
  /* allocate spectra and and pointer to return matrix on lua stack */
  Spec = (SpectraTypeArray *) lua_newuserdata(L, sizeof(SpectraTypeArray));
  luaL_setmetatable(L, "Spectra_Type");
  Spec->NSpectra=Options.NT1*Options.NT2*Options.Npsi;
  Spec->Spec= (SpectraType*)malloc(Spec->NSpectra*sizeof(SpectraType));
  /* test if H is Hermitian */
  switch(QOperatorHermitian(*H, 0.0001*Options.Settings.zero))
  {
    case 3 :
    {
      printf("WARNING: QOperatorHermitian failed in LuaCreateFluorescenceYield, this should not happen\n");
      printf("my guess is that you are out of memory and will crash soon\n");
      break;
    }
    case 2 :
    {
      luaL_warning(L, "CreateFluorescenceYield called with an operator that is not fully Hermitian\nStrongly sugested to stop the calculation and repair this\n\n\n\nDo it\n\n\nReally Check this\n\n\nDid you?\n\n\n");
      break;
    }
    case 1 :
    {
      // all fine H is Hermitian upto zero
      break;
    }
    case 0 :
    {
      luaL_error(L, "CreateFluorescenceYield called with an operator that is not Hermitian\n");
      break;
    }
  }
  /* calculate groundstate energy */
  psiEn= (double*)malloc(Options.Npsi*sizeof(double));
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
      Spec->Spec[Options.NT2*(iT*Options.Npsi+ipsi)].N=Options.Settings.NE;
      if(ipsi==0)
      {
        Spec->Spec[Options.NT2*(iT*Options.Npsi+ipsi)].Emin=Options.Settings.Emin;
        Spec->Spec[Options.NT2*(iT*Options.Npsi+ipsi)].Emax=Options.Settings.Emax;
      }
      else
      {
        Spec->Spec[Options.NT2*(iT*Options.Npsi+ipsi)].Emin=Spec->Spec[Options.NT2*(iT*Options.Npsi)].Emin;
        Spec->Spec[Options.NT2*(iT*Options.Npsi+ipsi)].Emax=Spec->Spec[Options.NT2*(iT*Options.Npsi)].Emax;
      }
      Spec->Spec[Options.NT2*(iT*Options.Npsi+ipsi)].Gamma=Options.Settings.Gamma;
      if(Options.Settings.restrictions.N==0)
      {
        FluorescenceYield(*H, T1[iT], T2, Options.NT2, psi[ipsi], psiEn[ipsi], Options.Settings.epsilon, Options.Settings.NTri, &(Spec->Spec[Options.NT2*(iT*Options.Npsi+ipsi)]));
      }
      else
      {
        FluorescenceYieldRestricted(*H, Options.Settings.restrictions, T1[iT], T2, Options.NT2, psi[ipsi], psiEn[ipsi], Options.Settings.epsilon, Options.Settings.NTri, &(Spec->Spec[Options.NT2*(iT*Options.Npsi+ipsi)]));
      }
    }
  }
  free(psiEn);
  return 1;
}

/*  Calculate the ground-state of an Hamiltonian in mean-field.
 *  On stack are
 *  1) operator
 *  2) density matrix to start from
 *  3) temperature kbT
 *
 *  on output the density matrix will be changed such that it resembles the density matrix of the ground-state
 *  aditionally the total energy is pushed to the stack
 */
static int LuaMeanFieldGroundstate(lua_State *L)
{
  OperatorType *H, *HMF;
  CompactMatrixType *Rho;
  double kbT;
  double Etot;
  H = (OperatorType *) luaL_checkudata(L, 1, "Operator_Type");
  Rho = (CompactMatrixType *) luaL_checkCompactMatrix(L, 2);
  kbT = (double) luaL_checknumber(L, 3);
  if((H->NFermion+H->NBoson)!=Rho->Ni)
  {
    return luaL_error(L, "Error in MeanFieldGroundstate, length of operator (NF %d NB %d) is not equal to length of density matrix (%d)\n",H->NFermion,H->NBoson,Rho->Ni);
  }
  MeanFieldGroundstate(H, Rho, kbT,&Etot);
  LuaPushTable_CompactMatrixType(L, Rho);
  lua_pushnumber(L, Etot);
  HMF = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  MeanFieldOperator(H, HMF, *Rho);
  luaL_setmetatable(L, "Operator_Type");
  FreeCompactMatrix(Rho);
  return 3;
}

static int LuaMeanFieldOperator(lua_State *L)
{
  OperatorType *H, *HMF;
  CompactMatrixType *Rho;
  uint32_t len = (uint32_t)lua_gettop(L);
  if(len<2 || len>2)
  {
    return luaL_error(L, "MeanFieldOperator called with %d arguments. Exactly two expected, Operator and Density Matrix.\n",len);
  }
  H = (OperatorType *) luaL_checkudata(L, 1, "Operator_Type");
  Rho = (CompactMatrixType *) luaL_checkCompactMatrix(L, 2);
  if((H->NFermion+H->NBoson)!=Rho->Ni)
  {
    return luaL_error(L, "Error in MeanFieldGroundstate, length of operator (NF %i NB %i) is not equal to length of density matrix (%i)\n",H->NFermion,H->NBoson,Rho->Ni);
  }
  HMF = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  MeanFieldOperator(H, HMF, *Rho);
  luaL_setmetatable(L, "Operator_Type");
  strncpy(HMF->Name, "", MAXStrLength-1); //An appropriate Name can be given in Lua.
  FreeCompactMatrix(Rho);
  return 1;
}

/*  transform an operator to a table
 *  on stack an operator 
 *  and optional either a wave-function, a list of wave functions or a density matrix
 *  For the case there is only an operator on stack a mean-field operator for 0 density is returned
 *  For the case there is a single function on stack returned is a table of size psi.N representing the Operator on a single determinant basis
 *  For the case there is a list of funcitons on stack returned is a table representing <psi_j O psi_i >
 *  For the case there is a density matrix on stack returned is a table representing the mean-field operator
 */
static int LuaOperatorToMatrix(lua_State *L)
{
  OperatorType *H;
  uint32_t len;
  len = (uint32_t)lua_gettop(L);
  H = (OperatorType *) luaL_checkudata(L, 1, "Operator_Type");
  if(len==1)
  {
  	if (H->MaxLength > 2 || H->l[0].N != 0 || H->l[1].N != 0)
  		VerbosityPrintf(0x0001,"\nWarning: Without a specified density the matrix representation can only be found for the quadratic part of an operator, all other parts will be ignored.\n");
    CompactMatrixType *Rho, HMat;
    Rho= (CompactMatrixType*)malloc(sizeof(CompactMatrixType));
    Rho->Ni=H->NFermion;
    Rho->Nj=H->NFermion;
    Rho->QComplex=0;
    InitCompactMatrix(Rho);
    MeanFieldOperatorOneElectronMatrix(H, &HMat, *Rho);
    LuaPushTable_CompactMatrixType(L, &HMat);
    FreeCompactMatrix(&HMat);
    FreeCompactMatrix(Rho);
    free(Rho);
    return 1;
  }
  else if (len==2)
  {
    switch (lua_type(L,2))
    {
      case LUA_TUSERDATA:
      {
        if (luaL_testudata(L, 2, "Wavefunction_Type") != NULL)
        {
          WaveFunctionType *psi;
          CompactMatrixType HMat;
          psi = (WaveFunctionType *) luaL_checkudata(L, 2, "Wavefunction_Type");
          OperatorToDeterminantMatrix(H, psi, &HMat);
          LuaPushTable_CompactMatrixType(L, &HMat);
          FreeCompactMatrix(&HMat);
          return 1;
        }
        return luaL_error(L, "OperatorToMatrix expects either a wave-function, table of functions or a density matrix on position two\n");
      }
      case LUA_TTABLE:
      {
        if(lua_objlen(L, 2)==0) return luaL_error(L, "Table of length 0 found at position 2 in OperatorToMatrix\n");
        lua_rawgeti(L, 2, 1);
        if (luaL_testudata(L, -1, "Wavefunction_Type") != NULL)
        {
          lua_pop(L,1);
          MatrixType res;
          WaveFunctionType *psi;
          psi = luaL_checkarray_WaveFunctionType(L, 2, &(res.Ni));
          res.Nj=res.Ni;
          res.QComplex=1;
          InitMatrix(&res);
          PsiOperatorPsiMatrix(*H, psi, &res);
          LuaPushTable_MatrixType(L, &res);
          FreeMatrix(&res);
          free(psi);
          return 1;
        }
        else
        {
          lua_pop(L,1);
          CompactMatrixType *Rho, HMat;
          Rho = (CompactMatrixType *) luaL_checkCompactMatrix(L, 2);
          MeanFieldOperatorOneElectronMatrix(H, &HMat, *Rho);
          LuaPushTable_CompactMatrixType(L, &HMat);
          FreeCompactMatrix(Rho);
          FreeCompactMatrix(&HMat);
          free(Rho);
          return 1;
        }
      }
      default:
      {
        return luaL_error(L, "Error in OperatorToMatrix, expected one (Operator) or two (Operator, Wavefunction / List of functions / DensityMatrix) input parameters, not (%d)\n",len);
      }
    }
  }
  else
  {
    return luaL_error(L, "Error in OperatorToMatrix, expected one (Operator) or two (Operator, Wavefunction / List of functions / DensityMatrix) input parameters, not (%d)\n",len);
  }
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
    {"Epsilon","double",&(Options.epsilon)},
    {"zero","double",&(Options.zero)},
    {"Zero","double",&(Options.zero)},
    {"restrictions", "DeterminantRestrictionType", &(Options.DeterminantRestrictions)},
    {"Restrictions", "DeterminantRestrictionType", &(Options.DeterminantRestrictions)},
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
 *        Pointer to a lua state. The last element on the stack needs to be a square table of tables of numbers (real or complex) or a response function
 * @return 2
 *         The first element added on the stack is a table of the eigenvalues
 *         The second element added on the stack is a table of the eigen functions
 */
static int LuaEigensystemDenseTable(lua_State *L)
{
  if (luaL_testudata(L, -1, "ResponseFunction_Type") != NULL) return LuaResponseFunctionEigensystem(L);
  if (luaL_testudata(L, -1, "Matrix_Type") != NULL) return LuaMatrixEigensystem(L);
  // last element on stack is not a user data storing a pointer to one of the types we can diagonalize, i.e. it must be a table
  if(lua_istable(L, -1)==0)
  {
    return luaL_error(L, "Error Eigensystem called with a single argument, which should be either a matrix (square table or one of the (block) response function types)");
  }
  // it is a table, either a square matrix or a reponse function type stored as table
  lua_getfield(L, -1, "type");
  if(lua_isstring(L,-1)==0)
  {
    // table without type definition, assume it is a square matrix
    lua_pop(L, 1);
    return LuaMatrixEigensystem(L);
  }
  lua_pop(L, 1);
  // last element is a table with a type definition, assume it is a response function stored as table
  return LuaResponseFunctionEigensystem(L);
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
  OperatorType *Opp;
  DeterminantRestrictionType StartRestrictions;
  StartRestrictions.N=0;
  WaveFunctionType *psi=NULL,*temp;
  uint32_t i,j,len, Npsi=0,QStartRestrictions=0;
  /* free garbage full cycle*/
  if(L->l_G->gcrunning) lua_gc(L, LUA_GCCOLLECT, 0);
  /* set options */
  LuaOptionsType LuaEigensystemOptions[]={
    {"DenseBorder", "uint32_t", &(Options.Settings.DenseBorder)},
    {"NKrylovStart","uint32_t",&(Options.Settings.NKrylovStart)},
    {"NKrylovMax","uint32_t",&(Options.Settings.NKrylovMax)},
    {"NKrylovStep","uint32_t",&(Options.Settings.NKrylovStep)},
    {"NKrylovMin","uint32_t",&(Options.Settings.NKrylovMin)},
    {"NBitsKey","uint32_t",&(Options.Settings.NBitsKey)},
    {"epsilon","double",&(Options.Settings.epsilon)},
    {"Epsilon","double",&(Options.Settings.epsilon)},
    {"zero","double",&(Options.Settings.zero)},
    {"Zero","double",&(Options.Settings.zero)},
    {"restrictions", "DeterminantRestrictionType", &(Options.DeterminantRestrictions)},
    {"Restrictions", "DeterminantRestrictionType", &(Options.DeterminantRestrictions)},
    {"ExpandBasis", "bool",&(Options.ExpandBasis)},
    {"ShiftScheme", "char*",Options.Settings.ShiftScheme},
    {NULL, NULL, NULL}};
  VerbosityPrintf(0x10000000,"Start of LuaEigensystem\n");
  /* check number of input parameters */
  len = (uint32_t)lua_gettop(L);
  if(len==1) return LuaEigensystemDenseTable(L);
  if(!((len==2)||(len==3)||(len==4)))
  {
    return luaL_error(L, "Error Eigensystem called with %d arguments, 1, 2, 3 or 4 expected\n",len);
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
      // read list of wave-functions to start from
      // The list is copied into a new array in order to have them in "order" in memory
      // Note that this means that when we change the wave-function the original function is not changed (except for the values it points to
      // we need to copy the changed function back into the original at the end.
      psi=luaL_checkarray_WaveFunctionType(L, 2, &Npsi);
      for(i=0;i<Npsi;i++)
      {
        for(j=0;j<i;j++)
        {
          if(psi[i].ValueRe==psi[j].ValueRe)
          {
            return luaL_error(L, "Eigensystem called with a list of functions. Function with name \"%s\" is at position %d and %d\n",psi[i].Name,j+1,i+1);
          }
        }
      }
      if(WaveFunctionArrayToSameBasis(psi, Npsi))
      {
        if(L->l_G->gcrunning)lua_gc (L, LUA_GCEMERGENCYCOLLECT, 0);
        if(WaveFunctionArrayToSameBasis(psi, Npsi))
        {
          return luaL_error(L, "Out of memory in Eigensystem");
        }
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
  /* read options */
  if((len==4)||((len==3)&&(!(QStartRestrictions))))
  {
    LuaCheckOptions(L, -1, LuaEigensystemOptions);
    if(Options.DeterminantRestrictions.N!=0)Options.restrictions=True;
    Options.Settings.DenseBorder = Max(Options.Settings.DenseBorder,Npsi);
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
      psi= (WaveFunctionType*)malloc(Npsi*sizeof(WaveFunctionType));
    }
  }
  /* test if operator is Hermitian */
  switch(QOperatorHermitian(*Opp, 0.0001*Options.Settings.zero))
  {
    case 3 :
    {
      printf("WARNING: QOperatorHermitian failed in LuaEigensystem, this should not happen\n");
      printf("my guess is that you are out of memory and will crash soon\n");
      break;
    }
    case 2 :
    {
      luaL_warning(L, "Eigensystem called with an operator that is not fully Hermitian\nStrongly sugested to stop the calculation and repair this\n\n\n\nDo it\n\n\nReally Check this\n\n\nDid you?\n\n\n");
      break;
    }
    case 1 :
    {
      // all fine H is Hermitian upto zero
      break;
    }
    case 0 :
    {
      luaL_error(L, "Eigensystem called with an operator that is not Hermitian\n");
      break;
    }
  }
  /* calculate groundstate */
  if(Options.ExpandBasis==True)
  {
    BlockGroundState(Opp, Npsi, &(Options.Settings) , StartRestrictions, Options.DeterminantRestrictions, psi, Options.Settings.NBitsKey);
  }
  else
  {
    BlockLanczosGroundStateConserveBasis(Opp, psi, Npsi, &(Options.Settings));
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
      free(psi); //free the array, as this is now stored in temp, but keep the pointers allocated within the structure
    }
    return 0;
  }
}


/*  Given an impurity Green's function, a Bath represented by a hybridization function and a set of transition operators acting on the impurity calculate
 *  < psi_i | T G T | psi_i >
 *  @param L
 *         Pointer to a lua state
 *         1-st  argument on stack should be an Hermitian operator (Hamiltonian) [H] with NI fermions
 *         2-nd  argument on stack should be the bath represented by the bath hybridization function. Given as a block list in either tri-diagonal, anderson or list of pole form matrices should have dimension NI by NI
 *         3-rd  argument is either an integer specifying the number of states to calculate or a list of wave-functions
 *         4-th  argument is a transition operator or list of transition operators
 *         Possible fifth argument can be a list of options
 *         - Tensor          boolean false
 *         - NTri            integer 200
 *         - OrderV          integer 1
 *         - LevelCon        integer 1
 *         - LevelVal        integer 1
 *         - RecalculatePsi  boolean true - can only be given if StartPsi is given
 *
 *  @return 3
 *         in adition changes the lua stack
 *         - First argument the Greens function in blocktridiagonal form
 *         - Second argument list of Npsi eigen-states of the extended impurity
 *         - Third argument the bath hybridization function in NIO form
 *
 */
static int LuaImpurityGreensFunction(lua_State * L)
{
  uint32_t len, Npsi;
  WaveFunctionType *psi;
  OperatorType *H,*T;
  ResponsefunctionType *DeltaBath;
  BlockNaturalImpurityOrbitalMatrixType *DeltaBathNIO=NULL;
  struct
  {
    ConvergenceSettingsType Settings;
    uint32_t NT;
    Bool Tensor;
    uint32_t OrderV;
    uint32_t LevelCon;
    uint32_t LevelVal;
    WaveFunctionType *psi;
    uint32_t Npsi;
    Bool RecalculatePsi;
    Bool RecalculateNIO;
  } Options;
  Options.Settings=InitConvergenceSettings();
  Options.Tensor=False;
  Options.OrderV=1;
  Options.LevelCon=1;
  Options.LevelVal=1;
  Options.psi=NULL;
  Options.Npsi=1;
  Options.RecalculatePsi=True;
  Options.RecalculateNIO=True;
  LuaOptionsType LuaImpurityGreensFunctionOptions[]={
    {"NTri","uint32_t",&(Options.Settings.NTri)},
    {"NBitsKey","uint32_t",&(Options.Settings.NBitsKey)},
    {"epsilon","double",&(Options.Settings.epsilon)},
    {"Epsilon","double",&(Options.Settings.epsilon)},
    {"restrictions", "DeterminantRestrictionType", &(Options.Settings.restrictions)},
    {"Restrictions", "DeterminantRestrictionType", &(Options.Settings.restrictions)},
    {"Tensor","bool",&(Options.Tensor)},
    {"OrderV","uint32_t",&(Options.OrderV)},
    {"LevelCon","uint32_t",&(Options.LevelCon)},
    {"LevelVal","uint32_t",&(Options.LevelVal)},
    {"RecalculatePsi","bool",&(Options.RecalculatePsi)},
    {"RecalculateNIO","bool",&(Options.RecalculateNIO)},
    {"MeanFieldStartingPointTemperature","double",&(Options.Settings.EffT)},
    {NULL, NULL, NULL}};
  /* check number of input parameters */
  len = (uint32_t)lua_gettop(L);
  if(!((len==4)||(len==5)))
  {
    return luaL_error(L, "Error ImpurityGreensFunction called with %d arguments, 4 or 5 expected\n",len);
  }
  /* read Hamiltonian */
  H = (OperatorType *) luaL_checkudata(L, 1, "Operator_Type");
  /* read hybridization function either in Natural Impurity Orbital basis, BlockListOfPoles, Anderson or Tridiagonal form */
  DeltaBath = luaL_checkResponsefunction(L,2);
  /* read psi or Npsi */
  if(lua_istable(L,3))
  {
    psi = luaL_checkarray_WaveFunctionType(L, 3, &Npsi);
    if(Npsi==0) return luaL_error(L, "Error ImpurityGreensFunction called with a table of wave-functions of length 0\n");
  }
  else
  {
    Npsi = luaL_checknumber(L, 3);
    psi = NULL;
  }
  /* read transtion operator */
  if(lua_type(L,4)==LUA_TTABLE)
  {
    T=luaL_checkarray_OperatorType(L, 4, &Options.NT);
    if(Options.NT==0) return luaL_error(L, "Error ImpurityGreensFunction called with a table of transition operators of length 0\n");
  }
  else
  {
    Options.NT = 1;
    T = (OperatorType *) luaL_checkudata(L, 4, "Operator_Type");
  }
  /* read options */
  if(len==5)
  {
    LuaCheckOptions(L, 5, LuaImpurityGreensFunctionOptions);
  }
  if((Options.RecalculatePsi==False)&&(psi==NULL)) return luaL_error(L, "Error ImpurityGreensFunction called option RecalculatePsi = false but without a wave-function\n");
  if((Options.RecalculateNIO==False)&&((DeltaBath->type|32)!='n')) return luaL_error(L, "Error ImpurityGreensFunction called option RecalculateNIO = false but without a hybridization function in natural impurity basis format\n");
  // calculate ground-state and natural impurity orbital representation
  if((Options.RecalculatePsi==True)&&(Options.RecalculateNIO==True))
  {
    if(ImpurityGroundState(H, DeltaBath, Options.Npsi, &Options.Settings, psi, DeltaBathNIO))
    {
      return luaL_error(L, "Error in ImpurityGroundState in ImpurityGreensFunction fatal sorry, should not happen\n");
    }
  }



  // calculate Greens function
  if(Options.Tensor==True)
  {
    // call tensor form


  }
  else
  {
    // call single form



  }
  if(DeltaBath->type&32)
  {
    FreeBlockMatrix(DeltaBath);
    free(DeltaBath);
  }



  return 3;
}

/*  Given a table of vectors or wave-functions returns a new set that is orthonormalized
 * @ param L
 *         Pointer to a lua state
 *         First argument on stack should be a table of complex vectors or wave-functions
 *         second optional argument can be a lsit of options.
 *         possible options are: (first option is standard)
 *         - Method:         string value: "Lowdin", "GramSchimdt"
 *         - Order:          string value: "Row", "Column"
 *         - SingularValue:  double: values equal or smaller than SingularValue are removed from basis
 *         - KeepDimensions: boolean (remove singular vectors or add a zero vector to the basis)
 */
static int LuaOrthonormalize(lua_State *L)
{
  int matorfun;
  uint32_t len;
  struct
  {
    char* Method;
    char* Order;
    double SingularValue;
    Bool ReduceDimension;
  } Options;
  LuaOptionsType LuaOrthonormalizeOptions[]={
    {"Method",          "char*",  &(Options.Method)},
    {"Order",           "char*",  &(Options.Order)},
    {"SingularValue",   "double", &(Options.SingularValue)},
    {"ReduceDimension", "bool",   &(Options.ReduceDimension)},
    {NULL, NULL, NULL}};
  Options.Method="Lowdin";
  Options.Order="Row";
  Options.SingularValue=1E-12;
  Options.ReduceDimension=True;
  /* free garbage full cycle*/
  if(L->l_G->gcrunning) lua_gc(L, LUA_GCCOLLECT, 0);
  // get options when defined
  len = (uint32_t)lua_gettop(L);
  if(len>2)
  {
    return luaL_error(L, "Orthogonalize called with %d arguments. One or Two expected\n",len);
  }
  if(len==2)
  {
    if(LuaTestOptionsList(L, -1))
    {
      LuaCheckOptions(L, -1, LuaOrthonormalizeOptions);
      if(!((strcmp(Options.Method, "Lowdin")==0)||(strcmp(Options.Method, "GramSchmidt")==0)))
      {
        return luaL_error(L, "Error Orthogonalize Option \"Method\" should either be \"Lowdin\" (standard) or \"GramSchmidt\" and not\"%s\"\n",Options.Method);
      }
      if(!((strcmp(Options.Order, "Row")==0)||(strcmp(Options.Order, "Column")==0)))
      {
        return luaL_error(L, "Error Orthogonalize Option \"Order\" should either be \"Row\" (standard) or \"Column\" and not\"%s\"\n",Options.Order);
      }
    }
    else
    {
      return luaL_error(L, "Second argument of Orthogonalize is optional but iff given should be a list of options.\n");
    }
  }
  // determine if we deal with a dense matrix or with a list of wave-functions
  luaL_checktype(L, 1, LUA_TTABLE);
  len = (uint32_t)lua_objlen(L, 1);
  if(len==0)
  {
    return luaL_error(L, "First argument of Orthogonalize should either be a table of wave-functions or a rectangular matrix (table of tables of equal lenght).\n");
  }
  lua_rawgeti(L, 1, 1);
  matorfun = lua_type(L, -1);
  lua_pop(L, 1);
  if(matorfun==LUA_TTABLE)
  {
    CompactMatrixType *M, *MO;
    M = (CompactMatrixType *) luaL_checkCompactMatrix(L, 1);
    MO = (CompactMatrixType*)malloc(sizeof(CompactMatrixType));
    if(CompactMatrixCopy(*M, MO))
    {
      if(L->l_G->gcrunning)lua_gc (L, LUA_GCEMERGENCYCOLLECT, 0);
      if(CompactMatrixCopy(*M, MO))
      {
        return luaL_error(L, "Not able to initialize matrix. Probably ran outoff memory on a bad moment\n");
      }
    }
    if(strcmp(Options.Method, "Lowdin")==0)
    {
      if(strcmp(Options.Order, "Row")==0)
      {
        CompactMatrixSVDOrthogonalizeRow(MO, Options.SingularValue, Options.ReduceDimension);
      }
      else if(strcmp(Options.Order, "Column")==0)
      {
        CompactMatrixSVDOrthogonalizeCol(MO, Options.SingularValue, Options.ReduceDimension);
      }
      else
      {
        return luaL_error(L, "Error Orthogonalize Option \"Order\" should either be \"Row\" (standard) or \"Column\" and not\"%s\"\n",Options.Order);
      }
    }
    else if(strcmp(Options.Method, "GramSchmidt")==0)
    {
      if(strcmp(Options.Order, "Row")==0)
      {
        CompactMatrixGramSchmidtOrthogonalizeRow(MO, Options.SingularValue, Options.ReduceDimension);
      }
      else if(strcmp(Options.Order, "Column")==0)
      {
        CompactMatrixGramSchmidtOrthogonalizeCol(MO, Options.SingularValue, Options.ReduceDimension);
      }
      else
      {
        return luaL_error(L, "Error Orthogonalize Option \"Order\" should either be \"Row\" (standard) or \"Column\" and not\"%s\"\n",Options.Order);
      }
    }
    else
    {
      return luaL_error(L, "Error Orthogonalize Option \"Method\" should either be \"Lowdin\" (standard) or \"GramSchmidt\" and not\"%s\"\n",Options.Method);
    }
    LuaPushTable_CompactMatrixType(L, MO);
    FreeCompactMatrix(MO);
    FreeCompactMatrix(M);
    free(MO);
    free(M);
  }
  else if(matorfun==LUA_TUSERDATA)
  {
    WaveFunctionType *psilist,*psilistNew;
    uint32_t i,Npsi;
    psilist = luaL_checkarray_WaveFunctionType(L, 1, &Npsi);
    psilistNew= (WaveFunctionType*)malloc(Npsi*sizeof(WaveFunctionType));
    if(psilistNew==NULL) return luaL_error(L, "Failure in Orthogonalize, not able to initialize pointers to list of orthogonal functions.\nOut of memory on a bad moment?\n");
    if(strcmp(Options.Method, "Lowdin")==0)
    {
      if(WaveFunctionListSVDOrthogonalizeDifferentBasis(psilist, psilistNew, &Npsi, Options.SingularValue, Options.ReduceDimension)) return luaL_error(L, "Failure in Orthogonalize Sorry\n");
    }
    else if(strcmp(Options.Method, "GramSchmidt")==0)
    {
      for(i=0;i<Npsi;i++)
      {
        if(WaveFunctionCopy(psilist[i], &(psilistNew[i]))) return luaL_error(L, "Failure in Orthogonalize, not able to initialize list of orthogonal functions.\nOut of memory on a bad moment?\n");
      }
      if(WaveFunctionListOrthogonalizeDifferentBasis(psilistNew, &Npsi, Options.SingularValue, Options.ReduceDimension)) return luaL_error(L, "Failure in Orthogonalize Sorry\n");
    }
    else
    {
      return luaL_error(L, "Error Orthogonalize Option \"Method\" should either be \"Lowdin\" (standard) or \"GramSchmidt\" and not\"%s\"\n",Options.Method);
    }
    LuaPushTable_WaveFunctionType(L, psilistNew, Npsi);
    free(psilist);
    free(psilistNew);
  }
  else
  {
    return luaL_error(L, "First argument of Orthogonalize should either be a table of wave-functions or a rectangular matrix (table of tables of equal lenght).\n");
  }
  return 1;
}

/*  Given a square hermitian matrix and a set of starting vectors returns a block band diagonal form
 * @ param L
 *         Pointer to a lua state
 *     Two posible inputs, either on dense matrices or on operators and wave-functions
 *     In both cases returns a BlockTridiagonalMatrix pushed on the lua stack
 *     Case 1
 *         First argument on stack should be a hermitian matrix represented as a table of table of (complex) numbers
 *         Second argument on stack should either be an integer representing the block size or a list of vectors representing the starting states
 *     Case 2
 *         First argument on stack should be a hermitian operator
 *         Second argument should be a single wave-function or a list of wavefunctions
 *     In both cases:
 *         Posible third argument specifies options
 *         - SingularValue:   double:  values equal or smaller than SignularValue are removed from basis (standard 10^-6)
 *         - NTri:            integer: Maximum number of blocks included (standard infinity)
 *         - NOrtho:          integer: Maximum number of reorthogonalizations (standard infinity)
 *         - ReOrthogonalize: boolean: Use aditional Gramm-Schmidt orthogonalization after the Lowdin orthogonalization (standard true)
 */
static int LuaBlockBandDiagonalize(lua_State *L)
{
  uint32_t len;
  BlockTriDiagonalMatrixType B;
  struct
  {
    double SingularValue;
    uint32_t NTri;
    uint32_t NOrtho;
    Bool ReOrthogonalize;
    Bool ReturnT;
  } Options;
  LuaOptionsType LuaBlockBandDiagonalizeOptions[]={
    {"SingularValue",   "char*",    &(Options.SingularValue)},
    {"NTri",            "uint32_t", &(Options.NTri)},
    {"NOrtho",          "uint32_t", &(Options.NOrtho)},
    {"ReOrthogonalize", "bool",     &(Options.ReOrthogonalize)},
    {"ReturnT",         "bool",     &(Options.ReturnT)},
    {NULL, NULL, NULL}};
  Options.SingularValue=1E-6;
  Options.NTri  =200;
  Options.NOrtho=UINT32_MAX;
  Options.ReOrthogonalize=True;
  Options.ReturnT=True;
  // get options when defined
  len = (uint32_t)lua_gettop(L);
  if(len>3)
  {
    return luaL_error(L, "BlockBandDiagonalize called with %d arguments. Two or Three expected\n",len);
  }
  if(len==3)
  {
    if(LuaTestOptionsList(L, -1))
    {
      LuaCheckOptions(L, -1, LuaBlockBandDiagonalizeOptions);
      if(Options.NOrtho<2)
      {
        return luaL_error(L, "Northo (%d) should be an intereger of at least 2\n",Options.NOrtho);
      }
    }
    else
    {
      return luaL_error(L, "Third argument of BlockBandDiagonalize is optional but iff given should be a list of options.\n");
    }
  }
  // case 1 - dense matrices
  if(lua_istable(L,1)==1)
  {
    CompactMatrixType *M,*vec;
    M = (CompactMatrixType *) luaL_checkCompactMatrix(L, 1);
    if(M->Ni != M->Nj)
    {
      return luaL_error(L, "BlockBandDiagonalize expects as a first argument either an operator or a square hermitian matrix, current matrix is not square: dimensions are %d by %d\n",M->Ni,M->Nj);
    }
    if(lua_type(L, 2)==LUA_TNUMBER)
    {
      vec= (CompactMatrixType*)malloc(sizeof(CompactMatrixType));
      vec->QComplex=0;
      vec->Ni=(uint32_t)luaL_checkinteger(L, 2);
      vec->Nj=M->Nj;
      if(vec->Ni>vec->Nj) return luaL_error(L, "Block size (%d) should be smaller than dimension of matrix (%d)\n",vec->Ni,vec->Nj);
      InitCompactIdentityMatrix(vec);
    }
    else
    {
      vec = (CompactMatrixType *) luaL_checkCompactMatrix(L, 2);
      if(vec->Ni>vec->Nj) return luaL_error(L, "Block size (%d) should be smaller than dimension of matrix (%d)",vec->Ni,vec->Nj);
      if(vec->Nj!=M->Nj) return luaL_error(L, "Length of starting vector (%d) should be equivalent to dimension of matrix (%d).\n",vec->Nj,M->Nj);
    }
    Options.NTri=MIN(Options.NTri, M->Ni);
    if(CompactMatrixBlockBandDiagonalize(M, vec, &B, Options.SingularValue, Options.NTri, Options.NOrtho, Options.ReOrthogonalize, Options.ReturnT))
    {
      return luaL_error(L, "CompactMatrixBlockBandDiagonalize failed that should not happen\n");
    }
    LuaPushTable_BlockTriDiagonalMatrixType(L, &B);
    FreeBlockTridiagonalMatrix(&B);
    FreeCompactMatrix(M);
    free(M);
    if(Options.ReturnT==True)
    {
      LuaPushTable_CompactMatrixType(L, vec);
      FreeCompactMatrix(vec);
      free(vec);
      return 2;
    }
    FreeCompactMatrix(vec);
    free(vec);
    return 1;
  }
  // case 2 - operator wave-functions
  else if (luaL_testudata(L, 1, "Operator_Type") != NULL)
  {
    uint32_t i;
    WaveFunctionType **orig=NULL;
    WaveFunctionType *psilist;
    OperatorType *O;
    uint32_t Npsi;
    // read the operator from the lua stack
    O = (OperatorType *) luaL_checkudata(L, 1, "Operator_Type");
    // read a list of pointers to the original functions from the lua stack
    if(lua_type(L,2)==LUA_TTABLE)
    {
      orig = luaL_checkarray_PointerWaveFunctionType(L, 2, &Npsi);    }
    else
    {
      Npsi=1;
      orig = (WaveFunctionType**)malloc(sizeof(WaveFunctionType*));
      orig[0] = (WaveFunctionType *) luaL_checkudata(L, 2, "Wavefunction_Type");
    }
    // copy the original functions that each live on their own basis to a new array where all functions have the same basis
    psilist = (WaveFunctionType*)malloc(Npsi*sizeof(WaveFunctionType));
    if(psilist==NULL)
    {
      return luaL_error(L, "malloc failed. Out of memory during initialization in BlockBandDiagonalize\n");
    }
    if(WaveFunctionCopy(*(orig[0]), &(psilist[0]))) return luaL_error(L, "malloc failed. Out of memory during initialization in BlockBandDiagonalize\n");
    for(i=1;i<Npsi;i++)
    {
      if(WaveFunctionAddMultipliedUnionBasis(&psilist[0], *(orig[i]), 0)) return luaL_error(L, "malloc failed. Out of memory during initialization in BlockBandDiagonalize\n");
    }
    for(i=1;i<Npsi;i++)
    {
      if(WaveFunctionInitCopyBasis(&(psilist[i]), psilist[0])) return luaL_error(L, "malloc failed. Out of memory during initialization in BlockBandDiagonalize\n");
      WaveFunctionAddConserveBasis(&(psilist[i]), *(orig[i]));
    }
    // remove the pointers to the original functions from memory. lua knows where they are
    free(orig);
    // call the c function
    B.N=Options.NTri;
    B.BlockSize=Npsi;
    B.QComplex=0;
    if(O->QComplex) B.QComplex=1;
    for(i=0;i<Npsi;i++)
    {
      if(psilist[i].QComplex) B.QComplex=1;
    }
    InitBlockTridiagonalMatrix(&B);
    if(LanczosBlockTriDiagonalize(*O, psilist, Options.SingularValue, &B, 0)) luaL_error(L, "LanczosBlockTriDiagonalize failed in BlockBandDiagonalize");
    LuaPushTable_BlockTriDiagonalMatrixType(L, &B);
    FreeBlockTridiagonalMatrix(&B);
    return 1;
  }
  // error not case 1 or 2
  else
  {
    return luaL_error(L, "BlockBandDiagonalize expects as a first argument either an operator or a square hermitian matrix (table of tables)\n");
  }
}


static int LuaInverse(lua_State *L)
{
  if(lua_gettop(L)!=1) return luaL_error(L, "Inverse expects a single matrix as input");
  CompactMatrixType *M;
  M = luaL_checkCompactMatrix(L, 1);
  if(CompactMatrixInverse(M)!=0) return luaL_error(L, "Inverse failed");
  LuaPushTable_CompactMatrixType(L, M);
  return 1;
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
  uint32_t i, j, len, Npsi, Nparticles;
  uint16_t *jindex;
  WaveFunctionType *psi,**psiList;
  MatrixType Rho, *RhoMat;
  Bool tensor=False;
  LuaOptionsType LuaDensityMatrixOptions[]={
    {"Tensor", "bool", &(tensor)},
    {NULL, NULL, NULL}};
  len = (uint32_t)lua_gettop(L);
  if(len==1)
  {
    // first argument is a wave-function or a list of functions
    switch (lua_type(L,1))
    {
      case LUA_TTABLE:
      {
        psiList=luaL_checkarray_PointerWaveFunctionType(L, 1, &Npsi);
        lua_createtable(L,Npsi,0);
        for(i=0;i<Npsi;i++)
        {
          CalculateDensityMatrix(&Rho, *psiList[i]);
          LuaPushTable_MatrixType(L, &Rho);
          lua_rawseti(L, -2, i+1);
          FreeMatrix(&Rho);
        }
        free(psiList);
        return 1;
      }
      case LUA_TUSERDATA:
      {
        if (luaL_testudata(L, 1, "Wavefunction_Type") != NULL)
        {
          psi = (WaveFunctionType *) luaL_checkudata(L, 1, "Wavefunction_Type");
          CalculateDensityMatrix(&Rho, *psi);
          LuaPushTable_MatrixType(L, &Rho);
          FreeMatrix(&Rho);
          return 1;
        }
        return luaL_error(L, "Error: DensityMatrix expects a function or a list of functions on the first position\n");
      }
    }
    return luaL_error(L, "Error: DensityMatrix expects a function or a list of functions on the first position\n");
  }
  else if (len==3)
  {
    LuaCheckOptions(L, 3, LuaDensityMatrixOptions);
    jindex = luaL_checkarray_unit16_t(L,2,&Nparticles);
    // first argument is a wave-function or a list of functions
    switch (lua_type(L,1))
    {
      case LUA_TTABLE:
      {
        if(tensor==False)
        {
          psiList=luaL_checkarray_PointerWaveFunctionType(L, 1, &Npsi);
          lua_createtable(L,Npsi,0);
          for(i=0;i<Npsi;i++)
          {
            CalculatePartialDensityMatrix(&Rho, *psiList[i], Nparticles, jindex);
            LuaPushTable_MatrixType(L, &Rho);
            lua_rawseti(L, -2, i+1);
            FreeMatrix(&Rho);
          }
          free(psiList);
          return 1;
        }
        else
        {
          /* calcualte partial density tensor for a list of functions */
          psi=luaL_checkarray_WaveFunctionType(L, 1, &Npsi);
          RhoMat= (MatrixType*)malloc(Npsi*Npsi*sizeof(MatrixType));
          CalculatePartialDensityMatrixWithOffDiagonals(RhoMat, Npsi, psi, Nparticles, jindex);
          lua_createtable(L,Npsi,0);
          for(i=0;i<Npsi;i++)
          {
            lua_createtable(L,Npsi,0);
            for(j=0;j<Npsi;j++)
            {
              LuaPushTable_MatrixType(L, &(RhoMat[i*Npsi+j]));
              lua_rawseti(L, -2, j+1);
              FreeMatrix(&(RhoMat[i*Npsi+j]));
            }
            lua_rawseti(L, -2, i+1);
          }
          free(psi);
          free(RhoMat);
          return 1;
        }
      }
      case LUA_TUSERDATA:
      {
        if(tensor==False)
        {
          if (luaL_testudata(L, 1, "Wavefunction_Type") != NULL)
          {
            psi = (WaveFunctionType *) luaL_checkudata(L, 1, "Wavefunction_Type");
            CalculatePartialDensityMatrix(&Rho, *psi, Nparticles, jindex);
            LuaPushTable_MatrixType(L, &Rho);
            FreeMatrix(&Rho);
            return 1;
          }
          return luaL_error(L, "Error: DensityMatrix expects a function or a list of functions on the first position\n");
        }
        return luaL_error(L, "Error: DensityMatrix expects a list of functions on the first position when the option {\"Tensor\",True} is set\n");
      }
    }
    return luaL_error(L, "Error: DensityMatrix expects a function or a list of functions on the first position\n");
  }
  else if (len==2)
  {
    /* position 1 is either a list of functions or a function */
    /* position 2 is either a list of integers or a list of options */
    luaL_checktype(L, 2, LUA_TTABLE);
    lua_rawgeti(L, 2, 1);
    switch (lua_type(L,-1))
    {
      case LUA_TTABLE:
      {
        lua_pop(L,1);
        LuaCheckOptions(L, 2, LuaDensityMatrixOptions);
        switch (lua_type(L,1))
        {
          case LUA_TTABLE:
          {
            if(tensor==False)
            {
              psiList=luaL_checkarray_PointerWaveFunctionType(L, 1, &Npsi);
              lua_createtable(L,Npsi,0);
              for(i=0;i<Npsi;i++)
              {
                CalculateDensityMatrix(&Rho, *psiList[i]);
                LuaPushTable_MatrixType(L, &Rho);
                lua_rawseti(L, -2, i+1);
                FreeMatrix(&Rho);
              }
              free(psiList);
              return 1;
            }
            else
            {
              /* calcualte density tensor for a list of functions */
              psi=luaL_checkarray_WaveFunctionType(L, 1, &Npsi);
              RhoMat= (MatrixType*)malloc(Npsi*Npsi*sizeof(MatrixType));
              CalculateDensityMatrixWithOffDiagonals(RhoMat, Npsi, psi);
              lua_createtable(L,Npsi,0);
              for(i=0;i<Npsi;i++)
              {
                lua_createtable(L,Npsi,0);
                for(j=0;j<Npsi;j++)
                {
                  LuaPushTable_MatrixType(L, &(RhoMat[i*Npsi+j]));
                  lua_rawseti(L, -2, j+1);
                  FreeMatrix(&(RhoMat[i*Npsi+j]));
                }
                lua_rawseti(L, -2, i+1);
              }
              free(psi);
              free(RhoMat);
              return 1;
            }
          }
          case LUA_TUSERDATA:
          {
            if(tensor==False)
            {
              if (luaL_testudata(L, 1, "Wavefunction_Type") != NULL)
              {
                psi = (WaveFunctionType *) luaL_checkudata(L, 1, "Wavefunction_Type");
                CalculateDensityMatrix(&Rho, *psi);
                LuaPushTable_MatrixType(L, &Rho);
                FreeMatrix(&Rho);
                return 1;
              }
              return luaL_error(L, "Error: DensityMatrix expects a function or a list of functions on the first position\n");
            }
            return luaL_error(L, "Error: DensityMatrix expects a list of functions on the first position when the option {\"Tensor\",True} is set\n");
          }
        }
        return luaL_error(L, "Error: DensityMatrix expects a function or a list of functions on the first position\n");
      }
      case LUA_TNUMBER:
      {
        lua_pop(L,1);
        jindex = luaL_checkarray_unit16_t(L,2,&Nparticles);
        switch (lua_type(L,1))
        {
          case LUA_TTABLE:
          {
            psiList=luaL_checkarray_PointerWaveFunctionType(L, 1, &Npsi);
            lua_createtable(L,Npsi,0);
            for(i=0;i<Npsi;i++)
            {
              CalculatePartialDensityMatrix(&Rho, *psiList[i], Nparticles, jindex);
              LuaPushTable_MatrixType(L, &Rho);
              lua_rawseti(L, -2, i+1);
              FreeMatrix(&Rho);
            }
            free(psiList);
            return 1;
          }
          case LUA_TUSERDATA:
          {
            if(tensor==False)
            {
              if (luaL_testudata(L, 1, "Wavefunction_Type") != NULL)
              {
                psi = (WaveFunctionType *) luaL_checkudata(L, 1, "Wavefunction_Type");
                CalculatePartialDensityMatrix(&Rho, *psi, Nparticles, jindex);
                LuaPushTable_MatrixType(L, &Rho);
                FreeMatrix(&Rho);
                return 1;
              }
              return luaL_error(L, "Error: DensityMatrix expects a function or a list of functions on the first position\n");
            }
            return luaL_error(L, "Error: DensityMatrix expects a list of functions on the first position when the option {\"Tensor\",True} is set\n");
          }
        }
        return luaL_error(L, "Error: DensityMatrix expects a function or a list of functions on the first position\n");
      }
    }
    return luaL_error(L, "Error: DensityMatrix position 2 should either be a list of integers selecting the orbitals to include or a list of options\n");
  }
  return luaL_error(L, "Error: DensityMatrix called with %d arguments, 1, 2 or 3 arguments expected", len);
}

static int LuaStatistics(lua_State *L)
{
  uint32_t len;
  len = (uint32_t)lua_gettop(L);


  if (luaL_testudata(L, 1, "Operator_Type") != NULL)
  {
    OperatorType *O = (OperatorType *) luaL_checkudata(L, 1, "Operator_Type");
    printf("Statistical information on operator %s\n",O->Name);
    //TODO:




    return 0;
  }
  else if (luaL_testudata(L, 1, "Wavefunction_Type") != NULL)
  {
    double *Hist;
    WaveFunctionType *psi = (WaveFunctionType *) luaL_checkudata(L, 1, "Wavefunction_Type");
    double exp = (double) luaL_checknumber(L, 2);
    uint32_t N = (uint32_t) luaL_checkunsigned(L, 3);
    Hist= (double*)malloc(N*sizeof(double));
    CalculatePrefactorHistogram(Hist, *psi, exp, N);
    LuaPushArray_double(L, Hist, N);
    return 1;
  }
  else
  {
    return luaL_error(L, "Error Statistics expects a wave-function or operator as first imput\n");
  }
}


int isIntElementOf(const uint16_t element, const uint16_t * const set, const uint32_t setSize)
{
	int i;
	for (i=0;i<setSize;i++)
	{
		if (element == set[i])
			return 1;
	}
	return 0;
}


static int LuaDeterminantString(lua_State *L)
{
	int NF;
	uint16_t * indices;
	uint32_t indSetSize;
	uint32_t len = (uint32_t)lua_gettop(L);
	if (len < 2 || len > 2)
		return luaL_error(L, "DeterminantString called with %d arguments. List of Indices and Number of Fermions expected.\n",len);
	indices = luaL_checkarray_unit16_t(L,1,&indSetSize);
	NF = luaL_checknumber(L,2);
	int i;
	for (i=0;i<indSetSize;i++)
	{
		if (indices[i]<0 || indices[i]>=NF)
			printf("\nWarning: Index %d is not within [0,NF-1]=[0,%d].\n", indices[i], NF-1);
	}
	char det[NF+1];
	det[NF] = '\0';
	for (i=0;i<NF;i++)
	{
		if (isIntElementOf(i,indices,indSetSize))
		{
			det[i] = '1';
		} else
		{
			det[i] = '0';
		}
	}
	lua_pushstring(L,det);
	return 1;
}

int nonRelAngularMomentumZNumbers(char * orb)
{
	if (strlen(orb)<1)
		return -1;
	char lastChar = orb[strlen(orb)-1];
	switch (lastChar) {
		case 's':
			return 2*0+1;
		case 'p':
			return 2*1+1;
		case 'd':
			return 2*2+1;
		case 'f':
			return 2*3+1;
		case 'g':
			return 2*4+1;
		case 'h':
			return 2*5+1;
		default:
			return -1;
	}
}

int relAngularMomentumZNumbers(char * orb)
{
	if (strlen(orb)<3)
		return -1;
	char lastChar, significantChar;
	lastChar = orb[strlen(orb)-1];
	if (lastChar != '2')
		return -1;
	significantChar = orb[strlen(orb)-2];
	if (significantChar == '/')
		significantChar = orb[strlen(orb)-3];
	switch (significantChar) {
			case '1':
				return 1+1;
			case '3':
				return 3+1;
			case '5':
				return 5+1;
			case '7':
				return 7+1;
			case '9':
				return 9+1;
			default:
				return -1;
		}
}




void lua_groupOrbitals(lua_State *L, int indTablePos, char* groupKey, char** groupMembers, uint32_t groupLen)
{
	int j,k,l;
	k = 1;
	lua_getfield(L,indTablePos,groupKey);
	if (!lua_isnil(L,-1))
	{
		printf("\nWarning: Group with name \"");
		printf("%s",groupKey);
		printf("\" already exists and cannot be included a second time.\n");
		return;
	}
	lua_pop(L, 1);
	lua_pushstring(L, groupKey);
	lua_newtable(L);
	int groupTablePos = (int) lua_gettop(L);
	for (j = 0; j < groupLen; j++)
	{
		lua_getfield(L, indTablePos, groupMembers[j]);
		if (lua_isnil(L, -1))
		{
			printf("\nWarning: Orbital with name \"");
			printf("%s",groupMembers[j]);
			printf("\" was not found and thus cannot be included in any grouping.\n");
			lua_pop(L, 1);
		} else
		{
			int orbTablePos = (int) lua_gettop(L);
			int listSize = (uint32_t) lua_objlen(L, orbTablePos);
			for (l = 1; l < listSize + 1; l++)
			{
				lua_pushnumber(L, k);
				lua_rawgeti(L, orbTablePos, l);
				lua_settable(L, groupTablePos);
				k++;
			}
			lua_pop(L, 1);
		}
	}
	lua_settable(L, indTablePos);
}

void lua_setNonRelAtomicOrbitalIndices(lua_State *L, char* orb, int momZNumbers, int * NF, int indTablePos, int isDict)
{
	int i,top, start;
	start = *NF;
	char * orbDn = (char*) malloc(strlen(orb) + strlen("_Dn"));
	strcpy(orbDn, orb);
	strcat(orbDn, "_Dn");
	char * orbUp = (char*) malloc(strlen(orb) + strlen("_Up"));
	strcpy(orbUp, orb);
	strcat(orbUp, "_Up");
	char * orbs[2] = { orbDn, orbUp };
	if (isDict)
		lua_pushstring(L, orbDn);
	else
		lua_pushnumber(L, (int) lua_objlen(L, indTablePos) + 1);
	lua_newtable(L);
	top = (int) lua_gettop(L);
	for (i = 0; i < momZNumbers; i++)
	{
		lua_pushnumber(L, i + 1);
		lua_pushnumber(L, start + 2 * i);
		lua_settable(L, top);
	}
	lua_settable(L, indTablePos);
	if (isDict)
		lua_pushstring(L, orbUp);
	else
		lua_pushnumber(L, (int) lua_objlen(L, indTablePos) + 1);
	lua_newtable(L);
	top = (int) lua_gettop(L);
	for (i = 0; i < momZNumbers; i++)
	{
		lua_pushnumber(L, i + 1);
		lua_pushnumber(L, start + 2 * i + 1);
		lua_settable(L, top);
	}
	lua_settable(L, indTablePos);
	if (isDict)
		lua_groupOrbitals(L, indTablePos, orb, orbs, 2);
	*NF += 2 * momZNumbers;
}

void lua_setRelAtomicOrbitalIndices(lua_State *L, char* orb, int momZNumbers, int * NF, int indTablePos, int isDict)
{
	int i, top, start;
	start = *NF;
	if (isDict)
		lua_pushstring(L, orb);
	else
		lua_pushnumber(L, (int) lua_objlen(L, indTablePos) + 1);
	lua_newtable(L);
	top = (int) lua_gettop(L);
	for (i = 0; i < momZNumbers; i++)
	{
		lua_pushnumber(L, i + 1);
		lua_pushnumber(L, start +  i);
		lua_settable(L, top);
	}
	lua_settable(L, indTablePos);
	*NF += momZNumbers;
}





void lua_setAtomicOrbitalIndices(lua_State *L, char* orb, int * NF, int indTablePos, int isDict)
{
	luaL_checktype(L, indTablePos, LUA_TTABLE);
	int momZNumbers = nonRelAngularMomentumZNumbers(orb);
	if (momZNumbers >= 0)
	{
		lua_setNonRelAtomicOrbitalIndices(L,orb,momZNumbers,NF,indTablePos, isDict);
		return;
	}
	momZNumbers = relAngularMomentumZNumbers(orb);
	if (momZNumbers >= 0)
	{
		lua_setRelAtomicOrbitalIndices(L,orb,momZNumbers,NF,indTablePos, isDict);
		return;
	}	else
	{
		printf("\nWarning: Key \"");
		printf("%s",orb);
		printf("\" could not be interpreted as an orbital\n");
		return;
	}
}





//Input is a list in the style of {"1s","2s","2p","3s"} or (relativistic) {"1s12","2s12","2p12","2p32"}
//Optional Input is a grouping of orbitals, like { {"unexcitedStates",{"1s","2s"} } , {"excitedStates",{"2p","3s"} } }
//Output 1 is a list {"1s"={0,1},"2s"={2,3},"2p"={4,5,6,7,8,9},"3s"={10,11},"unexcitedStates"={0,1,2,3},"All"={0,1,2,3,4,5,6,7,8,9,10,11},"1s_Dn"={0},"1s_Up"={1}, ...}
//Output 2 is the number of fermionic modes NF.
static int lua_createAtomicIndices(lua_State *L, int isDict)
{
	char** orbitals;
	uint32_t orbListLen;
	int indTablePos;
	int NF = 0;
	uint32_t groupListLen;
	char* groupKey;
	uint32_t groupLen;
	char** groupMembers;
	int i;
	orbitals = luaL_checkarray_string(L,1,&orbListLen);
	uint32_t numberOfInputs = (uint32_t) lua_gettop(L);
	if (numberOfInputs < 1 || numberOfInputs > 2 || (numberOfInputs > 1 && !isDict))
		return luaL_error(L, "CreateAtomicIndices called with %d arguments. List of Orbitals and (in case of a dictionary) optional parameter List of Groups of Orbitals expected.\n", numberOfInputs);
	lua_newtable(L);
	indTablePos  = (int) lua_gettop(L);
	for (i=0;i<orbListLen;i++)
	{
		lua_getfield(L,indTablePos,orbitals[i]);
		if (lua_isnil(L,-1))
		{
			lua_pop(L, 1);
			lua_setAtomicOrbitalIndices(L, orbitals[i], &NF, indTablePos, isDict);
		} else
		{
			lua_pop(L,1);
			printf("\nWarning: Orbital with name \"");
			printf("%s",orbitals[i]);
			printf("\" already exists and cannot be included a second time.\n");
		}
	}
	if (isDict)
	{
		lua_groupOrbitals(L, indTablePos, "All", orbitals, orbListLen);
		lua_groupOrbitals(L, indTablePos, "all", orbitals, orbListLen);
		if (numberOfInputs == 2)
		{
			//grouping of orbitals
			luaL_checktype(L, 2, LUA_TTABLE); //{{"unexcited",{"1s","2s"}},...}
			groupListLen = (uint32_t) lua_objlen(L, 2);
			for (i = 1; i < groupListLen + 1; i++)
			{
				lua_rawgeti(L, 2, i);
				luaL_checktype(L, -1, LUA_TTABLE); //{"unexcited",{"1s","2s"}}
				if ((int) lua_objlen(L, -1) != 2)
				{
					printf(
							"\nWarning: Optional groupings of orbitals need to be 2-element tables in the format {groupName,{listOfOrbitals}}. Table %d has %d elements.\n",
							i, (int) lua_objlen(L, -1));
				}
				if ((int) lua_objlen(L, -1) >= 2)
				{
					lua_rawgeti(L, -1, 1);
					groupKey = (char*) luaL_checkstring(L, -1); //"unexcited"
					lua_pop(L, 1);
					lua_rawgeti(L, -1, 2);
					luaL_checktype(L, -1, LUA_TTABLE); //{"1s","2s"}
					groupMembers = luaL_checkarray_string(L, -1, &groupLen);
					lua_groupOrbitals(L, indTablePos, groupKey, groupMembers, groupLen);
				}
			}
			lua_pop(L, (int) lua_gettop(L) - indTablePos);
		}
	}
	lua_pop(L,(int)lua_gettop(L)-indTablePos);
	lua_pushnumber(L,NF);
	return 2;
}

static int LuaCreateAtomicIndicesDict(lua_State *L)
{
	return lua_createAtomicIndices(L,1);
}


static int LuaCreateAtomicIndicesList(lua_State *L)
{
	return lua_createAtomicIndices(L,0);
}




static const struct luaL_Reg global_functions[] = {
  {"BlockBandDiagonalize", LuaBlockBandDiagonalize},
  {"Braket", LuaBraKet},
  {"CalculateBands", LuaCalculateBands},
  {"CalculateDOS", LuaCalculateDOS},
  {"CalculateFermi", LuaCalculateFermi},
  {"CalculateG", LuaCalculateG},
  {"CalculateRho", LuaCalculateRho},
  {"Chop", LuaChop},
  {"Clone", Luacopy},
  {"Conjugate", Luaconjugate},
  {"ConjugateTranspose", Luaconjugatetranspose},
  {"Copy", Luacopy},
	{"CreateAtomicIndicesDict",LuaCreateAtomicIndicesDict},
	{"CreateAtomicIndicesList",LuaCreateAtomicIndicesList},
  {"CreateSpectra", LuaCreateSpectra},
  {"CreateFluorescenceYield", LuaCreateFluorescenceYield},
  {"CreateResonantSpectra", LuaCreateResonantSpectra},
  {"CreateClusterHamiltonian", LuaCreateClusterHamiltonian},
  {"DensityMatrix",LuaDensityMatrix},
	{"DeterminantString",LuaDeterminantString},
  {"Eigensystem", LuaEigensystem},
  {"EvolveToGroundstate", LuaEvolveToGroundstate},
  {"Expand",LuaWavefunctionExpand},
  {"ExpandToBasis",LuaWaveFunctionToDeterminantBasis},
  {"FileReadStuttgartCTRL", LuaFileReadStuttgartCTRL},
  {"FileReadStuttgartHAMR", LuaFileReadStuttgartHAMR},
  {"FileReadDresdenFPLO", LuaFileReadDresdenFPLO},
//  {"Help", LuaHelp},
  {"ImpurityGreensFunction", LuaImpurityGreensFunction},
  {"Inverse", LuaInverse},
  {"MeanFieldGroundState", LuaMeanFieldGroundstate},
  {"MeanFieldOperator", LuaMeanFieldOperator},
  {"NewOperator", LuaNewOperator},
  {"NewTightBinding", LuaNewTightBinding},
  {"NewWavefunction", LuaNewWaveFunction},
  {"OperatorToMatrix", LuaOperatorToMatrix},
  {"Orthonormalize", LuaOrthonormalize},
  {"PotentialExpandedOnClm", LuaPotentialExpandedOnClm},
  {"PlotBands", LuaPlotBands},
  {"PlotBandsAndDOS", LuaPlotBandsAndDOS},
  {"PlotDOS", LuaPlotDOS},
  {"PlotFermi", LuaPlotFermi},
  {"Rotate", Luarotate},
  {"Statistics",LuaStatistics},
  {"TimeEnd",LuaTimeEnd},
  {"TimePrint",LuaTimePrint},
  {"TimeStart",LuaTimeStart},
  {"TightBindingDefFromStuttgartCTRLHAMR", LuaTightBindingDefFromStuttgartCTRLHAMR},
  {"TightBindingDefFromDresdenFPLO", LuaTightBindingDefFromDresdenFPLO},
  {"Transpose", Luatranspose},
  {"Verbosity",LuaVerbosity},
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
  luaL_newlib(L, Spectra_func);
  lua_setglobal(L, "Spectra");

  luaL_newmetatable(L, "TightBinding_Type");
  luaL_setfuncs(L, TightBinding_meta, 0);

  luaL_newmetatable(L, "ResponseFunction_Type");
  luaL_setfuncs(L, ResponseFunction_meta, 0);
  lua_newtable(L);
  luaL_setfuncs(L, ResponseFunctionTable_meta, 0);
  lua_setglobal(L, "ResponseFunctionMeta");
  luaMWH_RegisterGlobalFunctions(L, ResponseFunction_func);

  luaL_newmetatable(L, "Matrix_Type");
  luaL_setfuncs(L, Matrix_meta, 0);
  luaL_newlib(L, Matrix_func);
  lua_setglobal(L, "Matrix");
  lua_newtable(L);
  luaL_setfuncs(L, MatrixTable_meta, 0);
  lua_setglobal(L, "MatrixMeta");

  luaL_newmetatable(L, "StuttgartCTRL_Type");
  luaL_setfuncs(L, StuttgartCTRL_meta, 0);

  luaL_newmetatable(L, "StuttgartHAMR_Type");
  luaL_setfuncs(L, StuttgartHAMR_meta, 0);

  luaL_newmetatable(L, "DresdenFPLO_Type");
  luaL_setfuncs(L, DresdenFPLO_meta, 0);
  
  
  luaL_newmetatable(L, "Yi_Type");
  luaL_setfuncs(L, Yi_meta, 0);
  luaL_newlib(L, Yi_func);
  lua_setglobal(L, "Yi");
  
  luaL_newmetatable(L, "Martin_Type");
  luaL_setfuncs(L, Martin_meta, 0);
  luaL_newlib(L, Martin_func);
  lua_setglobal(L, "Martin");
  
  luaL_newmetatable(L, "Simon_Type");
  luaL_setfuncs(L, Simon_meta, 0);
  luaL_newlib(L, Simon_func);
  lua_setglobal(L, "Simon");

  luaMWH_RegisterGlobalFunctions(L, global_functions);
  luaMWH_RegisterGlobalFunctions(L, mathlib);

  return 0;
}
