
/* The basis set is given by single slater determinant wave functions
 * these are represented as an array of NFermion bits, whereby 0 means empty and 1 occupied
 * and NBoson bites (8 bits) where 0 means empty, 1 means an occupation of one and 255 means an occupation of 255, the maximum
 * A many particle wave function is stored as an array of prefactors and an array of corresponding
 * Single slater determinants. Only non zero (or large enough) prefactors are stored.
 * In order to find the position of a single slater determinant an hash function is defined
 * Each wave function has a corresponding hash value, wich ranges from 0 to 2^{NBitsKey}-1
 * The KeyList (hash table) contains an array of lenght 2^{NBitsKey} that contains an index
 * number to the position of this Determinant in the many particle wave function
 * The index 0 stands for not available, the wave function is storred at element 1
 * upto element N of the many particle wave function array
 *
 * As hash values might collide, each element also has an value of the position
 * of the next wave function with the same hash value
 *
 * In order to dinamically allocate the array of values, determinants and index of the
 * next element with the same hash value are storred in an allocated space of that can
 * be increased by 2^{NBitsBlok} elements at once.
 */


/* NBitsBlok             defines the size in 2^{NBitsBlok} that is added at each allocation
 * NBitsMiniBlok         defines the size of the memory that each processor has locally inside a parallel part of the program
 * MiniBlokSize          2^NBitsMiniBlok used often and did not want to locally recalculate
 * BlokMask              2^NBitsBlok-1 i.e. a set of 1's of length NBitsBlok, which can be used to moduleate a value to BlokSize
 * BlokSize              2^NBitsBlok used often and did not want to locally recalculate
 * BlokElement(i)        Function to translate a single index to a double index seperted to index within a blok and blok number
 * BlokElementDet(*,*,*) Same as for BlokElement, but now for a determinant, which has an aditional lenght
 * MAXStrLength          The fixed maximum length for strings, i.e. names and variables of matrices, vectors and functions
 *                       Do not go below 32, as some parts of the program do add at least this length to the name
 */
#define NBitsBlok 14
#define NBitsMiniBlok 7

#define MiniBlokSize ((uint32_t)1<<NBitsMiniBlok)

#define BlokMask (((uint32_t)1<<NBitsBlok)-1)
#define BlokSize ((uint32_t)1<<NBitsBlok)

#define BlokElement(i) (i)>>(NBitsBlok)][(i)&(BlokMask)
#define BlokElementDet(i,NWord,j) (i)>>(NBitsBlok)][NWord*((i)&(BlokMask))+(j)

#define MAXStrLength 256

#define CrAnMask ((uint16_t)1<<15)
#define NotCrAnMask (((uint16_t)1<<15)-1)

typedef int (*compfn) (const void *, const void *);

typedef struct Complex16
{
  double Re;
  double Im;
}
Complex16;

typedef enum { False, True } Bool;

/*
 *
 */
const struct { char *name; double value; char * units; char * description;} PhysicalConstants[] =
{
  // Particles
  {"malpha",     6.644657230E-27,     "kg",              "alpha particle mass in kilogram"},
  {"gd",         0.8574382311,        "",                "Deuteron g factor"},
  {"mud",        0.4330735040E-26,    "J T^-1",          "Deuteron magnetic moment in Joule per Tesla"},
  {"md",         3.343583719E-27,     "kg",              "Deuteron mass in kilogram"},
  {"lambdaCe",   2.4263102367E-12,    "m",               "Electron Compton wavelength in units of meter"},
  {"ge",         -2.00231930436182,   "",                "Electron g factor"},
  {"mue",        -928.4764620E-26,    "J T^-1",          "Electron magnetic moment in Joule per Tesla"},
  {"me",         9.10938356E-31,      "kg",              "Electron mass in kilogram"},
  {"lambdaCmu",  11.73444111E-15,     "m",               "Muon Compton wavelength in units of meter"},
  {"gmu",        -2.0023318418,       "",                "Muon g factor"},
  {"mumu",       -4.49044826E-26,     "J T^-1",          "Muon magnetic moment in Joule per Tesla"},
  {"mmu",        1.883531594E-28,     "kg",              "Muon mass in kilogram"},
  {"lambdaCn",   1.319590904E-15,     "m",               "Neutron Compton wavelength in units of meter"},
  {"gn",         -3.82608545,         "",                "Neutron g factor"},
  {"mun",        -0.96623650E-26,     "J T^-1",          "Neutron magnetic moment in Joule per Tesla"},
  {"mn",         1.674927471E-27,     "kg",              "Neutron mass in kilogram"},
  {"lambdaCp",   1.32140985396E-15,   "m",               "Proton Compton wavelength in units of meter"},
  {"gp",         5.585694702,         "",                "Proton g factor"},
  {"mup",        1.4106067873E-26,    "J T^-1",          "Proton magnetic moment in Joule per Tesla"},
  {"mp",         1.672621898E-27,     "kg",              "Proton mass in kilogram"},
  {"lambdaCtau", 0.697787E-15,        "m",               "Tau Compton wavelength in units of meter"},
  {"mtau",       3.167E-27,           "kg",              "Tau mass in kilogram"},
  // Aditional constants
  {"Z0",         376.730313461,       "Ohm",             "Characteristic impedance of vacuum in Ohm"},
  {"G",          6.67408E-11,         "m^3 kg^-1 s^-2",  "Newtonian constant of graviation in units of meter cubed per kilogram per second squared"},
  {"lPlanck",    1.616229E-35,        "m",               "Planck length in meter"},
  {"mPlanck",    2.176470E-8,         "kg",              "Planck mass in kilogram"},
  {"TPlanck",    1.416808E32,         "K",               "Planck temperature in Kelvin"},
  {"tPlanck",    5.39116E-44,         "s",               "Planck time in seconds"},
  {"G0",         7.7480917310E-5,     "S",               "Conductance quantum in units of Siemens"},
  {"R0",         12906.4037278,       "Ohm",             "Resistance quantum in units of Ohm"},
  {"KJ",         483597.8525E9,       "Hz V^_1",         "Josephson constant in Herz per Volt"},
  {"Phi0",       2.067833831E-15,     "Wb",              "Magnetic flux quantum in units of Weber"},
  {"RK",         25812.807455,        "Ohm",             "Klaus von Klitzing constant in units of Ohm"},
  {"re",         2.8179403227E-15,    "m",               "classical electron radius in meter"},
  {"mu",         1.660539040E-27,     "kg",              "Atomic mass constant in kilogram"},
  {"muN",        5.050783699E-27,    "J T^-1",          "Nuclear magneton in Joule per Tesla"},
  {"F",          96485.33289,         "C mol^-1",        "Faraday constant in units of Coulomb per mole"},
  {"R",          8.3144598,           "J mol^-1 K^-1",   "Molar gas constant in Joule per mole per Kelvin"},
  {"sigma",      5.670367E-8,         "W m^-2 K^-4",     "Stefan-Boltzmann constant in Wat per meter square per Kelvin to the forth power"},
  // often used constants
  {"e",          1.6021766208E-19,                     "C",              "Elementary charge in units of Coulomb"},
  {"a0",         0.52917721067E-10,                    "m",              "Bohr radius in meter"},
  {"c",          299792458,                            "m s^-1",         "speed of light in vacuum in meter per second"},
  {"epsilon0",   8.85418781762038985053656303171E-12,  "F m^-1",         "Electric constant (vacuum permittivity) in Farads per meter"},
  {"mu0",        1.25663706143591729538505735331E-6,   "N A^-2",         "Magnetic constant (vacuum permeability) in Newton per Ampare squared"},
  {"h",          6.626070040E-34,                      "J s^-1",         "Planck's constant in Joule per second"},
  {"hbar",       1.054571800E-34,                      "J s^-1",         "Planck's constant divided by 2 pi in Joule per second"},
  {"muB",        927.4009994E-26,                      "J T^-1",         "Bohr magneton in Joule per Tesla"},
  {"kB",         1.38064852E-23,                       "J K^-1",         "Boltzmann constant in units of Joule per Kelvin"},
  {"NA",         6.022140857E23,                       "mol^-1",         "Avogadro constant in units of number of particles per mole"},
  {"alpha",      7.2973525664E-3,                      "",               "fine-structure constant (1/(4 pi epsilon0) e^2/(hbar c)"}
};

/*
 *
 */
const struct { char *namelong; char *nameshort; double value; char * units; char * description;} EnergyUnits[] =
{
  // energy units in Joule
  {"Joule",         "J",     1,                          "J",               "Energy transferred (or work done) to an object when a force of one newton acts on that object in the direction of its motion through a distance of one meter"},
  {"Rydberg",       "Ryd",   2.179872324952945e-18,      "J",               "Rydberg, binding energy of a clasical electron to a point potential of 1 unit of electrical charge"},
  {"Hartree",       "Ha",    4.35974464990589e-18,       "J",               "Hartree energy (2 Rydberg)"},
  {"electron Volt", "eV",    1.6021766207999999e-19,     "J",               "The amount of energy gained (or lost) by the charge of a single electron moved across an electric potential difference of one volt"},
  {"Kelvin",        "K",     1.3806485200000002e-23,     "J",               "kB Boltzmann constant"},
  {"Tesla",         "T",     9.274009993999999e-24,      "J",               "muB Bohr magneton"},
  {"pernm",         "nm^-1", 1.9864458241717586e-16,     "J",               "10^9 h c"},
  {"percm",         "cm^-1", 1.9864458241717587e-23,     "J",               "100 h c"},
  {"TeraHz",        "THz",   6.626070040000001e-22,      "J",               "10^12 h"},
  {"perfs",         "fs^-1", 6.626070040000001e-19,      "J",               "10^15 h"}
};





/* Bosons are distinguished from fermions by the index number
 * Fermions are labled from 0 to NFermion-1
 * Bosons from NFermion to NFermion+NBoson-1
 *
 *  WARNING element [0][0] is not used. But in the calculation of the dot product is is summed over.
 *  It therefor must be set to zero. This is done in the initiallization.
 *  If you do the initiallization by hand one should not forget this.
 */
typedef struct WaveFunctionType
{
  char   Name[MAXStrLength];    /* A pointer to a name string that can be used to lable functions, good for printing etc. */
  uint32_t QComplex;            /* 0 for Real functions 1 for complex. For Real ValueIm is not allocated */
  uint32_t N;                   /* Total number of slater determinants used in this functions basis */
  uint32_t NBlok;               /* Number of Bloks allocated NBlok*(2^NBitsBlok) >= N */
  uint32_t NFermion;            /* Number of Fermions in the basis set */
  uint32_t NBoson;              /* Number of Bosons int he basis set, each boson has a length of 8 bit. */
  uint32_t NWord;               /* NFermionWord+NBoson
                                   For the next 3 variables: Element i can be found at position psi.ValueRe[i/(2^psi.NBitsBlok)][i%(2^psi.NBitsBlok)] */
  double **ValueRe;             /* Real part of Prefactor of corresponding Single Particle Determinant */
  double **ValueIm;             /* Imaginary part of Prefactor of corresponding Single Particle Determinant */
  uint32_t **Next;              /* If more elements have the same key the index number (i) of the next element with the same key */
  uint8_t **Determinant;        /* Single Slater Determinant as a pointer to a list of pointers which point to the begining of each blok of data of lenght (2^NBitsBlok)*(sizeof(uint32_t)) */
  uint32_t NBitsKey;            /* 2^NBitsKey is the number of unique keys used to index the different Slater integrals */
  uint32_t *KeyList;            /* A list of length 2^NBitsKey which has at index hash(psi.Determinant[i][j],NWord,NBitsKey) the index where the prefactor of the Single Slater Determinant can be found */
  /* If two Single Slater Determinants have the same hash function one of the two prefactors will be found at the position .Next of the index of the other.... */
}
WaveFunctionType;




typedef struct DeterminantRestrictionType
{
  uint32_t  N;                  /* number of restrictions */
  uint32_t  *Min;               /* minimum number of particles required */
  uint32_t  *Max;               /* maximum number of particles allowed */
  uint8_t   **Determinant;      /* an array of N vectors of length NWord of the wave-function applied to. Containing a bitfield wich filters the particles included in the sum */
}
DeterminantRestrictionType;

typedef struct OperatorFixedLengthType
{
  uint32_t N;                   /* number of operators of lenght i with i determined by O.l[i] */
  uint32_t NBitsKey;            /* number of bits in the key for this operator with this lenght */
  uint32_t NAloc;               /* Number of allocated opperators */
  uint32_t QComplex;            /* 0 for Real functions 1 for complex */
  uint16_t *o;                  /* array of length NAloc*i operator consistent of single product of lader operators */
  /* (o[j] & CrAnMask) determines creation (1) or annihilation (0) : (o[j] & NotCrAnMask) defines the index that has to be created or annihilated */
  uint32_t *Next;               /* array of length NAloc */
  uint32_t *KeyList;            /* key list */
  double *VRe;                  /* Real values */
  Complex16 *V;                 /* Complex values */
}
OperatorFixedLengthType;

typedef struct OperatorType
{
  char   Name[MAXStrLength];
  uint32_t QComplex;              /* If all operators with non-zero value of O.l[ilength].N are real    O.QComplex=0; * If all operators with non-zero value of O.l[ilength].N are complex O.QComplex=1; * If there is a mixture of real and complex values                   O.QComplex=2; */
  uint32_t NFermion;              /* Number of fermions in the basis */
  uint32_t NBoson;                /* Number of bosons in the basis */
  uint32_t NWord;                 /* NFermionWord+NBoson */
  uint32_t MaxLength;             /* Largest number of terms of lader operators */
  OperatorFixedLengthType *l;     /* array of length MaxLength+1 */
  DeterminantRestrictionType Res; /* set of restrictions applied to Operator psi calculations */
}
OperatorType;


typedef struct ConvergenceSettingsType
{
  uint32_t DenseBoarder;    /* If matrices have dimensions smaller than DenseBoarder dense methods are used - Normally this is determined by psi.N */
  uint32_t NKrylovStart;    /* In tridiagonalization routines the number of Krylov states created in the first itteration */
  uint32_t NKrylovMax;      /* In tridiagonalization routines the maximum number of Krylov states created */
  uint32_t NKrylovStep;     /* the number of Krylov states added when a restart is done */
  uint32_t NKrylovMin;      /* the minumum number of Krylov states included before schwitsching to an algorithm that recalculates the Krylovs states */
  uint32_t NTri;            /* number of Krylov states in single shot tri-diagonalization routines used for spectra */
  uint32_t NTri1;           /* number of Krylov states in single shot tri-diagonalization routines used for resonant spectra */
  uint32_t NTri2;           /* number of Krylov states in single shot tri-diagonalization routines used for resonant spectra */
  double   zero;            /* smallest number different from one in convergence tests */
  double   epsilon;         /* smallest absolute prefactor of a determinant in wave-functions considered different from 0 */
  double   EffT;            /* effective temperature in mean-field convergence loops*/
  uint32_t NBitsKey;        /* number of bits in key */
  DeterminantRestrictionType restrictions;  /* restrictions on possible determinatns */
  DeterminantRestrictionType restrictions1; /* restrictions on possible determinatns in resonant spectra */
  DeterminantRestrictionType restrictions2; /* restrictions on possible determinatns in resonant spectra */
  double Emin;
  double Emax;
  uint32_t NE;
  double Gamma;
  double Emin1;
  double Emax1;
  uint32_t NE1;
  double Gamma1;
  double Emin2;
  double Emax2;
  uint32_t NE2;
  double Gamma2;
  char   *ShiftScheme;
}
ConvergenceSettingsType;

ConvergenceSettingsType InitConvergenceSettings()
{
  char *ShiftScheme = "Variance";
  ConvergenceSettingsType ret=
  {
    1513 /*DenseBoarder*/,
    100 /*NKrylovStart*/,
    500 /*NKrylovMax*/,
     10 /*NKrylovStep*/,
     50 /*NKrylovMin*/,
    200 /*NTri*/,
    100 /*NTri1*/,
    100 /*NTri2*/,
    100*sqrt(DBL_EPSILON) /*zero, i.e. smallest varieance to determine convergence */,
    sqrt(DBL_EPSILON)/100 /*epsilon i.e. smallest prefactor of function to keep */,
    0.00008621738         /* 1K, EffT */,
    16                    /*bits in key*/,
    {0, NULL, NULL, NULL} /*restricitons*/,
    {0, NULL, NULL, NULL} /*restricitons1*/,
    {0, NULL, NULL, NULL} /*restricitons2*/,
    DBL_MAX,-DBL_MAX,1000,0 /*Energy range */,
    DBL_MAX,-DBL_MAX,1000,0 /*Energy range 1*/,
    DBL_MAX,-DBL_MAX,1000,0 /*Energy range 2*/,
    ShiftScheme             /* shift according to variance (options Automatic, Variance, None*/
  };
  return ret;
};

typedef struct TimeType
{
  uint32_t N;
  uint32_t NAloc;
  char **Name;
  uint32_t *TotalTime;
  time_t *BeginTime;
  uint32_t *NumberOfRuns;
  uint32_t *NowRunning;
}
TimeType;

/*  MatrixType defines a 2-dimensional tensor (matrix) with Ni rows and Nj columns
 *  The left top element is A[0][0], the right bottom element is A[Ni-1][Nj-1]
 *  The label QComplex sets a real (0) or complex (1) matrix.
 *
 *     ( A[0][0]    | A[0][1]    | A[0][2]    | ... | A[0][Nj-1]    )
 *   A=( A[1][0]    | A[1][1]    | A[1][2]    | ... | A[1][Nj-1]    )
 *     ( A[2][0]    | A[2][1]    | A[2][2]    | ... | A[2][Nj-1]    )
 *     ( .......    | .......    | .......    | ... | A[2][Nj-1]    )
 *     ( A[Ni-1][0] | A[Ni-1][1] | A[Ni-1][2] | ... | A[Ni-1][Nj-1] )
 *
 *   (A B)[i][j] = \sum_k A[i][k] B[k][j]
 *
 *   y[i] = \sum_j A[i][j] x[j]
 *
 *   A[i][j]=A[0][i*A.Nj+j]
 */
typedef struct MatrixType
{
  char   Name[MAXStrLength];    /* A pointer to a name string that can be used to lable functions, good for printing etc. */
  uint32_t Ni, Nj;              /* dimmensions of the matrix */
  uint32_t QComplex;            /* complex or real */
  double **VRe, **VIm;          /* ellments are found at M.VRe[i][j] */
  double a0;                    /* storing a possible real shift if a continued fraction notation is used */
  double b0;                    /* storing a possilbe over-all scaling if a continued fraction notation is used */
  double mu;
}
MatrixType;

/*  CompactMatrixType defines a 2-dimensional tensor (matrix) with Ni rows and Nj columns
 *  The left top element is A(0,0), the right bottom element is A(Ni-1,Nj-1)
 *  The label QComplex sets a real (0) or complex (1) matrix.
 *
 *     ( A(0,0)    | A(0,1)    | A(0,2)    | ... | A(0,Nj-1)    )
 *   A=( A(1,0)    | A(1,1)    | A(1,2)    | ... | A(1,Nj-1)    )
 *     ( A(2,0)    | A(2,1)    | A(2,2)    | ... | A(2,Nj-1)    )
 *     ( .......   | .......   | .......   | ... | A(2,Nj-1)    )
 *     ( A(Ni-1,0) | A(Ni-1,1) | A(Ni-1,2) | ... | A(Ni-1,Nj-1) )
 *
 *   (A B)(i,j) = \sum_k A(i,k) B(k,j)
 *
 *   y[i] = \sum_j A(i,j) x[j]
 *
 *   element (i,j) is found at position [i*A.Nj+j]
 */
typedef struct CompactMatrixType
{
  char *Name;                   /* A pointer to a name string that can be used to lable functions, good for printing etc. */
  uint32_t Ni, Nj;              /* dimmensions of the matrix */
  uint32_t QComplex;            /* complex or real */
  union
  {
    double *VRe;
    Complex16 *V;
  };
  double a0;                    /* storing a possible real shift if a continued fraction notation is used */
  double b0;                    /* storing a possilbe over-all scaling if a continued fraction notation is used */
  double mu;
}
CompactMatrixType;

typedef struct CompactVectorType
{
  uint32_t Ni;                  /* Length of vector */
  uint32_t QComplex;            /* complex or real */
  union
  {
    double *VRe;
    Complex16 *V;
  };
}
CompactVectorType;

/*  A sparse matrix indexed with the use of a hash function.
 *  KeyList and Next elements that are zero point to nowhere
 *  Other values point to elments 1 to the first element found at [0] etc.
 *
 *  A sparse matrix uses the same multiplication rules as a full matrix, but only stores nonzero elements.
 */
typedef struct SparseMatrixType
{
  char   Name[MAXStrLength];    /* A pointer to a name string that can be used to lable functions, good for printing etc. */
  uint32_t Ni, Nj;              /* dimmention of the matrix */
  uint32_t N;                   /* number of non zero elements */
  uint32_t NAloc;               /* number of allocated memory places for storage */
  uint32_t QComplex;            /* complex or real? */
  double *VRe, *VIm;            /* list of real and complex values */
  uint32_t *i, *j;              /* list of index numbers */
  uint32_t *Next;               /* where to find other elements with the same hash funciton */
  uint32_t *KeyList;            /* hash lookup table */
  uint32_t NBitsKey;            /* number of bits in the hash key */
  double a0;                    /* storing a possible real shift if a continued fraction notation is used */
  double b0;                    /* storing a possilbe over-all scaling if a continued fraction notation is used */
  double mu;
}
SparseMatrixType;

/*  A square real, hermitian matrix, with only the leading diagonal [i][i] and first of diagonal [i][i+1] [i+1][i] non-zero
 *
 *  The element a[0] stores a possilbe real ofset when used to represent a continued fraction
 *  The element b[0] stores a possile weighting when used to represent a continued fraction
 *  The element a[i] stores the elements M[i-1][i-1] of the full matrix for i>=1
 *  The element b[i] stores the elements M[i][i-1]==M[i-1][i] of the full matrix for i>=1
 *  a is a vector of lenght N+1
 *  b is a vector of length N
 *  The matrix represented has a dimension of N*N
 *
 *     ( A[0][0]=a1 | A[0][1]=b1 | 0          | 0          )
 *   A=( A[1][0]=b1 | A[1][1]=a2 | A[1][2]=b2 | 0          )
 *     ( 0          | A[2][1]=b2 | A[2][2]=a3 | A[2][3]=b3 )
 *     ( 0          | 0          | A[3][2]=b3 | A[3][3]=a4 )
 *
 *   (A B)[i][j] = \sum_k A[i][k] B[k][j]
 *
 *    and defines the Green's function
 *    G0 = a0 + b0^2 . ( G1 )
 *       with Gi = 1/(w - ai - bi^2 ( Gi+1 ) )
 *       and Gn = 1/(w - an )
 *       and ai and bi real
 */
typedef struct TriDiagonalMatrixType
{
  char   Name[MAXStrLength];    /* A pointer to a name string that can be used to lable functions, good for printing etc. */
  uint32_t N;                   /* dimension of matrix */
  double *a;                    /* diagonal elements of matrix (length N+1)   element a[i+1]=M[i,i]            of the full matrix */
  double *b;                    /* off diagonal elements of matrix (length N) element b[i+1]=M[i,i+1]=M[i+1,i] of the full matrix */
  double mu;                    /* chemical potential if transformed to Greens function */
}
TriDiagonalMatrixType;

/*  A square real, hermitian matrix, with only the leading diagonal [i][i]  and first row and column [i][0] [0][i] non-zero
 *
 *  The element a[0] stores a possilbe real ofset when used to represent a continued fraction
 *  The element b[0] stores a possile weighting when used to represent a continued fraction
 *  The element a[i] stores the elements M[i-1][i-1] of the full matrix for i>=1
 *  The element b[i] stores the elements M[i][0]==M[0][i] of the full matrix for i>=1
 *  a is a vector of lenght N+1
 *  b is a vector of length N
 *  The matrix represented has a dimension of N*N
 *
 *     ( A[0][0]=a1 | A[0][1]=b1 | A[0][2]=b2 | A[0][3]=b3 )
 *   A=( A[1][0]=b1 | A[1][1]=a2 | 0          | 0          )
 *     ( A[2][0]=b2 | 0          | A[2][2]=a3 | 0          )
 *     ( A[3][0]=b3 | 0          | 0          | A[3][3]=a4 )
 *
 *   (A B)[i][j] = \sum_k A[i][k] B[k][j]
 *
 *    and defines the Green's function
 *    G0 = a0 + b0^2 . ( G1 )
 *       with Gn = 1/(w - an - bn^2 ( Gn+1 ) )
 *       and an and bn real
 */
typedef struct AndersonMatrixType
{
  char   Name[MAXStrLength];    /* A pointer to a name string that can be used to lable functions, good for printing etc. */
  uint32_t N;                   /* dimension of matrix */
  double *a;                    /* diagonal elements of matrix (length N+1)   elements a[i+1]=M[i,i] of the full matrix */
  double *b;                    /* off diagonal elements of matrix (length N) elements b[i]=M[0,i]=M[i,0] of the full matrix */
  double mu;                    /* chemical potential if transformed to Greens function */
}
AndersonMatrixType;

/*  The structure ListOfPoles can represent a function as a sum of Lorenzians
 *  g(w) = a0 + sum_i b_i / (w - a_i)
 */
typedef struct ListOfPoles
{
  char   Name[MAXStrLength];    /* A pointer to a name string that can be used to lable functions, good for printing etc. */
  uint32_t N;                   /* dimension of list */
  double *a;                    /* energy of pole */
  double *b;                    /* weight of pole */
  double a0;                    /* a possible real ofset */
  double mu;                    /* chemical potential if transformed to Greens function */
}
ListOfPoles;



/*
 *     ( A[0][0]=avc0| A[0][1]=bv0 | 0           | 0           | A[0][4]=bc0 | 0           | 0           )
 *   A=( A[1][0]=bv0 | A[1][1]=av1 | A[1][2]=bv1 | 0           | 0           | 0           | 0           )
 *     ( 0           | A[2][1]=bv1 | A[2][2]=av2 | A[2][3]=bv2 | 0           | 0           | 0           )
 *     ( 0           | 0           | A[3][2]=bv2 | A[3][3]=av3 | 0           | 0           | 0           )
 *     ( A[4][0]=bc0 | 0           | 0           | 0           | A[4][4]=ac1 | A[4][5]=bc1 | 0           )
 *     ( 0           | 0           | 0           | 0           | A[5][4]=bc1 | A[5][5]=ac2 | A[5][6]=bc2 )
 *     ( 0           | 0           | 0           | 0           | 0           | A[6][5]=bc2 | A[6][6]=ac3 )
 *
 *   avc0 = (a0val + a0con)/2
 *
 *   (A B)[i][j] = \sum_k A[i][k] B[k][j]
 *
 *   Note that the block A[1;;3][1;;3] and the block A[4;;6][4;;6] are normal tri-diagonal matrices
 */
typedef struct DoubleTriDiagonalMatrixType
{
  char   Name[MAXStrLength];    /* A pointer to a name string that can be used to lable functions, good for printing etc. */
  double mu;                    /* border if a number belongs to val or con */
  double Epsilon;               /* values in a range of [Ef-1/2 Epsilon , Ef+1/2 Epsilon] contribute half to val and half to con */
  double a0;                    /* storing a possible real shift if a continued fraction notation is used */
  double b0;                    /* storing a possilbe over-all scaling if a continued fraction notation is used */
  TriDiagonalMatrixType val;
  TriDiagonalMatrixType con;
}
DoubleTriDiagonalMatrixType;


/*     ( A[0][0]=Dt.a0  | A[0][1]=Dt.b0  | 0              | 0           | 0           | 0           | 0           | 0           )
 *     ( A[1][0]=Dt.b0  | A[1][1]=Dt.avc0| A[1][2]=Dt.bv0 | 0           | 0           | A[1][5]=bc0 | 0           | 0           )
 *   A=( 0              | A[2][1]=Dt.bv0 | A[2][2]=Dt.av1 | A[2][3]=bv1 | 0           | 0           | 0           | 0           )
 *     ( 0              | 0              | A[3][2]=Dt.bv1 | A[3][3]=av2 | A[3][4]=bv2 | 0           | 0           | 0           )
 *     ( 0              | 0              | 0              | A[4][3]=bv2 | A[4][4]=av3 | 0           | 0           | 0           )
 *     ( 0              | A[5][1]=Dt.bc0 | 0              | 0           | 0           | A[5][5]=ac1 | A[5][6]=bc1 | 0           )
 *     ( 0              | 0              | 0              | 0           | 0           | A[6][5]=bc1 | A[6][6]=ac2 | A[6][7]=bc2 )
 *     ( 0              | 0              | 0              | 0           | 0           | 0           | A[7][6]=bc2 | A[7][7]=ac3 )
 *
 *   avc0 = (a0val + a0con)/2
 *
 *   (A B)[i][j] = \sum_k A[i][k] B[k][j]
 *
 *   Note that the block A[1;;3][1;;3] and the block A[4;;6][4;;6] are normal tri-diagonal matrices
 */
typedef struct SecondDoubleTriDiagonalMatrixType
{
  char   Name[MAXStrLength];    /* A pointer to a name string that can be used to lable functions, good for printing etc. */
  double mu;                    /* border if a number belongs to val or con */
  double Epsilon;               /* values in a range of [Ef-1/2 Epsilon , Ef+1/2 Epsilon] contribute half to val and half to con */
  double a0;                    /* storing a possible real shift if a continued fraction notation is used */
  double b0;                    /* storing a possilbe over-all scaling if a continued fraction notation is used */
  DoubleTriDiagonalMatrixType Dt;
}
SecondDoubleTriDiagonalMatrixType;

/*     ( A[0][0]=Dt.a[0] | A[0][1]=Dt.b[0] | 0               | 0               | 0               )
 *     ( A[1][0]=Dt.b[0] | A[1][1]=Dt.a[1] | A[1][2]=Dt.b[1] | A[1][2]=Dt.b[2] | A[1][2]=Dt.b[3] )
 *   A=( 0               | A[2][1]=Dt.b[1] | A[2][2]=Dt.a[2] | 0               | 0               )
 *     ( 0               | A[1][2]=Dt.b[2] | 0               | A[2][2]=Dt.a[3] | 0               )
 *     ( 0               | A[1][2]=Dt.b[3] | 0               | 0               | A[2][2]=Dt.a[4] )
 */
typedef struct SecondAndersonMatrixType
{
  char   Name[MAXStrLength];    /* A pointer to a name string that can be used to lable functions, good for printing etc. */
  double mu;                    /* border if a number belongs to val or con */
  double Epsilon;               /* values in a range of [Ef-1/2 Epsilon , Ef+1/2 Epsilon] contribute half to val and half to con */
  double a0;                    /* storing a possible real shift if a continued fraction notation is used */
  double b0;                    /* storing a possilbe over-all scaling if a continued fraction notation is used */
  AndersonMatrixType Dt;
}
SecondAndersonMatrixType;


/*
 *     ( A[0][0]=acvB  | A[0][0]=acvBA | A[0][1]=bvB0  | 0             | 0           | A[0][4]=bcB0  | 0             | 0           )
 *     ( A[0][0]=acvBA | A[0][0]=acvA  | A[0][1]=bvA0  | 0             | 0           | A[0][4]=bcA0  | 0             | 0           )
 *     ( A[1][0]=bvB0  | A[1][0]=bvA0  | A[1][1]=av2   | A[1][2]=bv2   | 0           | 0             | 0             | 0           )
 *   A=( 0             | 0             | A[2][1]=bv2   | A[2][2]=av3   | A[2][3]=bv3 | 0             | 0             | 0           )
 *     ( 0             | 0             | 0             | A[3][2]=bv3   | A[3][3]=av4 | 0             | 0             | 0           )
 *     ( A[4][0]=bcB0  | A[4][0]=bcA0  | 0             | 0             | 0           | A[4][4]=ac2   | A[4][5]=bc2   | 0           )
 *     ( 0             | 0             | 0             | 0             | 0           | A[5][4]=bc2   | A[5][5]=ac3   | A[5][6]=bc3 )
 *     ( 0             | 0             | 0             | 0             | 0           | 0             | A[6][5]=bc3   | A[6][6]=ac4 )
 *
 *  A = T^* Ap T^T
 *
 *     ( A[0][0]=av1   | 0             | A[0][1]=bv1   | 0             | 0           | 0             | 0             | 0           )
 *     ( 0             | A[0][0]=ac1   | 0             | 0             | 0           | A[0][4]=bc1   | 0             | 0           )
 *     ( A[1][0]=bv1   | 0             | A[1][1]=av2   | A[1][2]=bv2   | 0           | 0             | 0             | 0           )
 *  Ap=( 0             | 0             | A[2][1]=bv2   | A[2][2]=av3   | A[2][3]=bv3 | 0             | 0             | 0           )
 *     ( 0             | 0             | 0             | A[3][2]=bv3   | A[3][3]=av4 | 0             | 0             | 0           )
 *     ( 0             | A[4][0]=bc1   | 0             | 0             | 0           | A[4][4]=ac2   | A[4][5]=bc2   | 0           )
 *     ( 0             | 0             | 0             | 0             | 0           | A[5][4]=bc2   | A[5][5]=ac3   | A[5][6]=bc3 )
 *     ( 0             | 0             | 0             | 0             | 0           | 0             | A[6][5]=bc3   | A[6][6]=ac4 )
 *
 *
 *     ( A[0][0]=bv0   | A[0][0]=bc0   | 0             | 0             | 0           | 0             | 0             | 0           )
 *     ( A[0][0]=bc0   | A[0][0]=-bv0  | 0             | 0             | 0           | 0             | 0             | 0           )
 *     ( 0             | 0             | 1             | 0             | 0           | 0             | 0             | 0           )
 *   T=( 0             | 0             | 0             | 1             | 0           | 0             | 0             | 0           )
 *     ( 0             | 0             | 0             | 0             | 1           | 0             | 0             | 0           )
 *     ( 0             | 0             | 0             | 0             | 0           | 1             | 0             | 0           )
 *     ( 0             | 0             | 0             | 0             | 0           | 0             | 1             | 0           )
 *     ( 0             | 0             | 0             | 0             | 0           | 0             | 0             | 1           )
 *
 *   (A B)[i][j] = \sum_k A[i][k] . B[k][j]
 *
 *   Note that the block A[2;;4][2;;4] and the block A[5;;7][5;;7] are normal tri-diagonal matrices
 *   The states av1 and ac1 in Ap are linear combinations of the impurity site and a special linear combination of the other sites
 *   such that the impurity has an occupation n, the other site an occupation 1-n, av1 an occupation 1 and ac1 an occupaiton of 0.
 */
typedef struct BADoubleTriDiagonalMatrixType
{
  char   Name[MAXStrLength];    /* A pointer to a name string that can be used to lable functions, good for printing etc. */
  double mu;                    /* border if a number belongs to val or con */
  double Epsilon;               /* values in a range of [Ef-1/2 Epsilon , Ef+1/2 Epsilon] contribute half to val and half to con */
  double a0;                    /* storing a possible real shift if a continued fraction notation is used */
  double b0;                    /* storing a possilbe over-all scaling if a continued fraction notation is used */
  double acvB;                  /* val.b[0]^2 * val.a[1] + con.b[0]^2 * con.a[1] */
  double acvA;                  /* val.b[0]^2 * con.a[1] + con.b[0]^2 * val.a[1] */
  double acvBA;                 /* val.b[0]*con.b[0] * (val.a[1]-con.a[1]) */
  double bvB0;                  /* val.b[0] * val.b[1] */
  double bcB0;                  /* con.b[0] * con.b[1] */
  double bvA0;                  /* con.b[0] * val.b[1] */
  double bcA0;                  /*-val.b[0] * con.b[1] */
  TriDiagonalMatrixType val;
  TriDiagonalMatrixType con;
}
BADoubleTriDiagonalMatrixType;



/*  V is an array of length (N+1)*(1+BlockSize*BlockSize)-1
 *  (capitals for matrices, small cap for numbers)
 *  in the array V one can find N blocks of size 1+BlockSize*BlockSize in these blocks
 *  The first number is a
 *    a_i can be found at position (i-1)*(1+BlockSize*BlockSize)
 *  The last BlockSize*BlockSize is the matrix B
 *    B_i can be found at position (i-1)*(1+BlockSize*BlockSize)+1 to i*(1+BlockSize*BlockSize)-1
 *  The last BlockSize*BlockSize elements store the matrix A0
 *    A0 can be found at position N*(1+BlockSize*BlockSize) to (N+1)*(1+BlockSize*BlockSize)-2
 *
 *  The Green's function represented is:
 *    G = A0 + sum_i B'_i 1/(w-a_i I) B'_i^C
 *    with I the identity matrix of size BlockSize
 *  this is equivalent to
 *    G = A0 + sum_i B'_i . B'_i^C * 1/(w-a_i)
 *  use
 *    B_i = B'_i B'_i^C to get to
 *    G = A0 + sum_i B_i * 1/(w-a_i)
 *
 *  a0 is an additional pointer to V[N*(1+BlockSize*BlockSize)]
 */
typedef struct BlockListOfPoles
{
  char Name[MAXStrLength];   /* A pointer to a name string that can be used to lable functions, good for printing etc. */
  uint32_t N;                /* dimension of list */
  uint32_t BlockSize;        /* dimension of square matrices */
  uint32_t QComplex;         /* matrices are complex or real */
  union                      /* list of (N+1)*(1+BlockSize*BlockSize)-1 with the values */
  {
    double *VRe;
    Complex16 *V;
  };
  union                      /* matrix of size BlockSize*BlockSize with possible ofset */
  {
    double *a0Re;
    Complex16 *a0;
  };
  double mu;                 /* chemical potential */
}
BlockListOfPoles;

/*
 *     ( A[0][0]=a1   | A[0][1]=b1   | 0            | 0          )
 *   A=( A[1][0]=b1^+ | A[1][1]=a2   | A[1][2]=b2   | 0          )
 *     ( 0            | A[2][1]=b2^+ | A[2][2]=a3   | A[2][3]=b3 )
 *     ( 0            | 0            | A[3][2]=b3^+ | A[3][3]=a4 )
 *
 *   with ai and bi matrices
 *
 *   (A B)(i,j) = \sum_k A(i,k) B(k,j)
 *
 *   y[i] = \sum_j A(i,j) x[j]
 *
 *   element (i,j) is found at position [i*Nj + j]
 *
 *   matrix a0 and b0 are used as a possible ofset and scaling when transformed to Green's functions
 *
 *   matrix k can be found at position [k*BlockSize*BlockSize]
 *   element [i,j] of matrix k can be found at position [k*BlockSize*BlockSize+i*BlockSize+j]
 *
 *   when represented as a Green's function the function is:
 *   G0 =        a0 + b0^* . ( G1 ) . b0^T
 *   G1 = 1/(w - a1 - b1 . ( G2 ) . b1^+)
 *   G2 = 1/(w - a2 - b2 . ( G3 ) . b2^+)
 *   ...
 *   Gn = 1/(w - an)
 *
 *   note that the later needs to be changed to
 *   G0 =        a0 + b0^* . ( G1 ) . b0^T
 *   G1 = 1/(w - a1 - b1^* . ( G2 ) . b1^T)
 *   G2 = 1/(w - a2 - b2^* . ( G3 ) . b2^T)
 *   ...
 *   Gn = 1/(w - an)
 *
 *   in order to be consistent and easily allow to create H from G^{-1}
 */
typedef struct BlockTriDiagonalMatrixType
{
  char   Name[MAXStrLength];   /* A pointer to a name string that can be used to label functions, good for printing etc. */
  uint32_t N;                  /* number of diagonal blocks */
  uint32_t BlockSize;          /* dimension of square matrices */
  uint32_t QComplex;           /* matrices are complex or real */
  double *aRe;                 /* List of (N+1) * BlockSize*BlockSize values */
  Complex16 *a;
  double *bRe;                 /* List of (N) * BlockSize*BlockSize values */
  Complex16 *b;
  double mu;                   /* chemical potential */
}
BlockTriDiagonalMatrixType;

/*
 *     ( A[0][0]=a1   | A[0][1]=b1   | A[1][2]=b2   | A[2][3]=b3 )
 *   A=( A[1][0]=b1^+ | A[1][1]=a2   | 0            | 0          )
 *     ( A[2][1]=b2^+ | 0            | A[2][2]=a3   | 0          )
 *     ( A[3][2]=b3^+ | 0            | 0            | A[3][3]=a4 )
 *
 *   (A B)(i,j) = \sum_k A(i,k) B(k,j)
 *
 *   y[i] = \sum_j A(i,j) x[j]
 *
 *   element (i,j) is found at position [i*Nj + j]
 *
 *   matrix a0 and b0 are used as a possible ofset and scaling when transformed to Green's functions
 *
 *   matrix k can be found at position [k*BlockSize*BlockSize]
 *   element [i,j] of matrix k can be found at position [k*BlockSize*BlockSize+i*BlockSize+j]
 */
typedef struct BlockAndersonMatrixType
{
  char   Name[MAXStrLength];   /* A pointer to a name string that can be used to lable functions, good for printing etc. */
  uint32_t N;                  /* number of diagonal blocks */
  uint32_t BlockSize;          /* dimension of square matrices */
  uint32_t QComplex;           /* matrices are complex or real */
  union                        /* List of (N+1) * BlockSize*BlockSize values */
  {
    double *aRe;
    Complex16 *a;
  };
  union                        /* List of (N) * BlockSize*BlockSize values */
  {
    double *bRe;
    Complex16 *b;
  };
  double mu;                    /* chemical potential if transformed to Green's function */
}
BlockAndersonMatrixType;



/*
 *     ( A[0][0]=Ai    | A[0][1]=Bib   | A[0][2]=Biv   | 0             | 0           | A[0][5]=Bic   | 0             | 0           )
 *     ( A[1][0]=Bib^+ | A[1][1]=Ab    | A[1][2]=Bbv   | 0             | 0           | A[0][5]=Bbc   | 0             | 0           )
 *     ( A[2][0]=Biv^+ | A[2][1]=Bbv^+ | A[2][2]=av2   | A[2][3]=bv2   | 0           | 0             | 0             | 0           )
 *   A=( 0             | 0             | A[3][2]=bv2^+ | A[3][3]=av3   | A[3][4]=bv3 | 0             | 0             | 0           )
 *     ( 0             | 0             | 0             | A[4][3]=bv3^+ | A[4][4]=av4 | 0             | 0             | 0           )
 *     ( A[5][0]=Bic^+ | A[5][1]=Bbc^+ | 0             | 0             | 0           | A[5][5]=ac2   | A[5][6]=bc2   | 0           )
 *     ( 0             | 0             | 0             | 0             | 0           | A[6][5]=bc2^+ | A[6][6]=ac3   | A[6][7]=bc3 )
 *     ( 0             | 0             | 0             | 0             | 0           | 0             | A[7][6]=bc3^+ | A[7][7]=ac4 )
 *
 *  A = T Ap T^C
 *
 *     ( A[0][0]=av1   | 0             | A[0][2]=bv1   | 0             | 0           | 0             | 0             | 0           )
 *     ( 0             | A[1][1]=ac1   | 0             | 0             | 0           | A[0][5]=bc1   | 0             | 0           )
 *     ( A[2][0]=bv1^+ | 0             | A[2][2]=av2   | A[2][3]=bv2   | 0           | 0             | 0             | 0           )
 *  Ap=( 0             | 0             | A[3][2]=bv2^+ | A[3][3]=av3   | A[3][4]=bv3 | 0             | 0             | 0           )
 *     ( 0             | 0             | 0             | A[4][3]=bv3^+ | A[4][4]=av4 | 0             | 0             | 0           )
 *     ( 0             | A[5][1]=bc1^+ | 0             | 0             | 0           | A[5][5]=ac2   | A[5][6]=bc2   | 0           )
 *     ( 0             | 0             | 0             | 0             | 0           | A[6][5]=bc2^+ | A[6][6]=ac3   | A[6][7]=bc3 )
 *     ( 0             | 0             | 0             | 0             | 0           | 0             | A[7][6]=bc3^+ | A[7][7]=ac4 )
 *
 *
 *     ( A[0][0]=bv0   | A[0][1]=bc0   | 0             | 0             | 0           | 0             | 0             | 0           )
 *     ( A[1][0]=~bc0  | A[1][1]=~-bv0 | 0             | 0             | 0           | 0             | 0             | 0           )
 *     ( 0             | 0             | 1             | 0             | 0           | 0             | 0             | 0           )
 *   T=( 0             | 0             | 0             | 1             | 0           | 0             | 0             | 0           )
 *     ( 0             | 0             | 0             | 0             | 1           | 0             | 0             | 0           )
 *     ( 0             | 0             | 0             | 0             | 0           | 1             | 0             | 0           )
 *     ( 0             | 0             | 0             | 0             | 0           | 0             | 1             | 0           )
 *     ( 0             | 0             | 0             | 0             | 0           | 0             | 0             | 1           )
 *
 *  with,
 *   T T^C = 1
 *
 *   (A B)(i,j) = \sum_k A(i,k) B(k,j)
 *
 *   y[i] = \sum_j A(i,j) x[j]
 *
 *   element (i,j) is found at position [i*Nj + j]
 *
 *   Note that the block A[2;;4][2;;4] and the block A[5;;7][5;;7] are normal tri-diagonal matrices
 *   The states av1 and ac1 in Ap are linear combinations of the impurity site and a special linear combination of the other sites
 *   such that the impurity has an occupation n, the other site an occupation 1-n, av1 an occupation 1 and ac1 an occupaiton of 0.
 */
typedef struct BlockNaturalImpurityOrbitalMatrixType
{
  char   Name[MAXStrLength];   /* A pointer to a name string that can be used to lable functions, good for printing etc. */
  uint32_t BlockSize;          /* dimension of square matrices */
  uint32_t QComplex;           /* matrices are complex or real */
  double   mu;
  double Epsilon;              /* values in a range of [Ef-1/2 Epsilon , Ef+1/2 Epsilon] contribute half to val and half to con */

  BlockTriDiagonalMatrixType val;
  BlockTriDiagonalMatrixType con;
  CompactMatrixType TMT;
  CompactMatrixType a0;
}
BlockNaturalImpurityOrbitalMatrixType;


/*  ResponsefunctionType defines a single type that can hold (Block) matrices of different forms
 *  l or L for (Block)ListOfPoles
 *  a or A for (Block)AndersonMatrix
 *  t or T for (Block)TridiagonalMatrix
 *  n or N for (Block)NaturalMatrix
 *
 *  Boolean Block is true for block matrtices, false for non-block matrices
 *
 *  the difference between small and captical is used for a container pointing to a unique object (capital)
 *  or to a container pointing to a shared object (small)
 *
 *  unique objects should be freed once not used any more
 *  shared objects are controled by other parts of the code and one should realize that changing these objects
 *  will influence other calculations as well.
 *
 *  The capital lower case is used for different representations on the lua stack.
 *  lower case letters are pushed to the lua stack as a pointer to the variable (user type)
 *  upper case letters are pushed to the lua stack as a lua table
 *  reading a user type will read the pointer
 *  reading a lua table will regenerate the ResponsefunctionType
 *
 *  The array SuggestedGrid contains NSuggestedGridPoints in ascending order.
 *  The Suggested Grid is used for transformations of Matrix types and functions acting on them
 *  The actual grid used can be different and will be different in several cases.
 *  i.e. transformation functions will change the BlockMatrix as little as possible, but
 *  if the size of the BlockMatrix would grow (i.e. for transformations to list of poles)
 *  and a Grid is provided the Grid is used
 */

typedef struct ResponsefunctionType
{
  char type;                               // l,a,t,n or L,A,T,N
  Bool Block;
  void * data;
  double *SuggestedGrid;
  uint32_t NSuggestedGridPoints;
}
ResponsefunctionType;



/*
 *
 */
typedef struct LuaOptionsType
{
  char * name;
  char * type;
  void * variable;
}
LuaOptionsType;






typedef struct ChebyshevExpansion
{
  int    N;
  double *mu;
  double Emin;
  double Emax;
}
ChebyshevExpansion;

typedef struct HarmonicsExpansion
{
  uint32_t QComplex;            /* complex or real? */
  uint32_t N;
  int   *l;
  int   *m;
  double *Re;
  double *Im;
}
HarmonicsExpansion;


typedef struct SpectraType
{
  char   Name[MAXStrLength];    /* A pointer to a name string that can be used to lable functions, good for printing etc. */
  uint32_t N;                   /* Number of points in spectrum run from 0 to N */
  double *Re;                   /* Array of lenght N+1 containing the Re values */
  double *Im;                   /* Array of lenght N+1 containing the Im values */
  double Emin;                  /* lowest value of the spectrum */
  double Emax;                  /* highest value of the spectrum */
  double dE;                    /* dE=(Emax-Emin)/N */
  double EminPole;              /* Energy below wich no poles excist and the spectrum is given by a pure laurent series */
  double EmaxPole;              /* Energy below wich no poles excist and the spectrum is given by a pure laurent series */
  double Gamma;                 /* FWHM broadening also the distance by half away from the real axes */
  double phase;                 /* Needed for a correct definition of the moments, a and b below */
  double ZerothMoment;          /* Norm of the Im part of (e^{-i phase}) G */
  double FirstMoment;           /* First moment of the Im part of (e^{-i phase}) G */
  double a[5];                  /* First five values of tri-diagonal representation of the spectrum given by (e^{-i phase}) G */
  double b[5];                  /* First five values of tri-diagonal representation of the spectrum given by (e^{-i phase}) G */
  double ResonantEnergy;        /* NAN in general, value of resonant energy in the case of RIXS */
}
SpectraType;

typedef struct SpectraTypeArray
{
  uint32_t NSpectra;
  SpectraType *Spec;
}
SpectraTypeArray;

/* Orbitals per atom have to be defined later */

typedef struct AtomType
{
  char   Name[MAXStrLength];         /* Atom name */
  uint32_t AtomNumber;               /* atomic number of the atom */
  double Position[3];                /* position of the atom */
  int    Nshell;                     /* number of shells in the atom */
  char  **ShellName;                 /* Name of the shells in the atom */
  int   *shelldegeneracy;            /* degeneracy of the shell */
  char  ***OrbitalName;              /* Name of the orbitals in the different shells */
}
AtomType;

/*  rUnitCell contains the unit cell in Angstrom
 *  gUnitCell contains the reciprocal lattice vectors
 *  NAtoms the number of atoms within the basis
 */
typedef struct CrystalStructureType
{
  double rUnitCell[3][3];       /* Unit cell vectors */
  double gUnitCell[3][3];       /* Unit cell vectors */
  uint32_t NAtoms;              /* Number of atoms wihtin the unit cell */
  AtomType *Atom;               /* Array of length NAtoms specifying the atoms with the unit cell */
  uint32_t Units;               /* **** ***x 0 r.g = 1 | 1 r.g = 2 Pi
                                 * **** *xx* 00 Angstrom | 01 a0 | 10 nm | 11 ??
                                 * **** x*** 0 Positions in units described above | 1 Postions in terms of r[0] ,r[1] and r[2]
                                 */
}
CrystalStructureType;

/*  rSuperCell contains the definition of the supercell in repetitions of the Crystal unit cell
 *  CellPos conatins the (0,0,0) corner of the crystal unit cell in crystal unit cell coordinates
 */
typedef struct CrystalStructureSuperCellType
{
  int32_t rSuperCell[3][3];     /* Number of unit cells within super cell */
  CrystalStructureType Crystal; /* definition of the unit cell */
  /* additional info that can be derived from the above settings */
  uint32_t NCell;               /* total number of cells in the super cell */
  int32_t CellPos[5000][3];      /* position of cells, not read in but calculated from above info */
  double gSuperCell[3][3];      /* reciprocal unit cell */
  double kVector[5000][3];       /* k-points wihtin super cell */
}
CrystalStructureSuperCellType;

typedef struct TightBindingDefHopping
{
  int    Atompos1;
  int    Shellpos1;
  int    Atompos2;
  int    Shellpos2;
  int    TransVecs[3];        /* translation between unit cells */
  double ConnectingVector[3]; /* translation in units according to Units including atom positions insite the cell */
  MatrixType Hop;
}
TightBindingDefHopping;

typedef struct TightBindingDef
{
  char   Name[MAXStrLength];
  CrystalStructureType Cell;
  int    nH;
  TightBindingDefHopping *H;
  int    NFermion;
  int  **FermionIndex;
}
TightBindingDef;

typedef struct BandstructurePointType
{
  int    Nq;
  char   label;
  double q[3];
}
BandstructurePointType;

typedef struct BandstructurePathType
{
  int    NQ;
  BandstructurePointType *Q;
  uint32_t Units;               /* **** ***x 0 r.g = 1 | 1 r.g = 2 Pi
                                 * **** *xx* 00 Angstrom | 01 a0 | 10 nm | 11 ??
                                 * **** x*** 0 Positions in units described above | 1 Postions in terms of r[0] ,r[1] and r[2]
                                 */
}
BandstructurePathType;













typedef struct DMFTSettings
{
  char   G0FileName[256];       /* filename of G0 */
  char   RootFileName[256];     /* root of GBath, GCluster and Sigma */
  char   SpecRootFileName[256]; /* filename of spectral functions in spectral representation */
  char   WaveFunctionRootFileName[256]; /* filename of wave-function  */
  char   LocalHamiltonianFilename[256]; /* Input file for the local Hamiltonian */
  char   GSymFilename[256];     /* File specifying the symmetry of G */
  uint32_t NSite;               /* Number of cluster sites and or times spin orbitals */
  double Nelectron;             /* Filling of the impurity */
  double mu;                    /* Chemical potential of the lattice Greensfunction */
  double mu0;                   /* Chemical potential of G0 */
  double Teff;                  /* effective temperature when determening the relation between chemical potential and filling for the bath */
  double T;                     /* real temperature in the cluster not jet implemented */
  uint32_t NmuSym;              /* 0000 0000 mu0 and mu from input file are used as starting values */
  /*           mu is optimized during the loops in order to fix the occupation to Nelectron */
  /* 0000 0001 Before start chemical potential mu and mu0 are read from G0 and Sigma input files */
  /*           These override the values in the input file. (Otherwise values in input file are used */
  /* 0000 0010 Nelectron is ignored and the calculation is preformed at fixed mu */
  /*           otherwise mu is optimized every loop */
  uint32_t iStart;              /* start value of loop number */
  uint32_t Ni;                  /* number of loops to do */
  uint32_t Dyson;               /* 0 no Dyson like perturbation - 1 full perturbation - 2 ... */
  uint32_t SymOpt;              /* 0000 0000 no aditional symmetrization */
  /* 0000 0001 An additional Operator of (U+1) * S^2 is added to the Hamiltonian to stabilize the singlet state */
  /* 0000 0010 The Greenfunction is averaged such that G(w) = - G(-w)^* */
  /* 0000 0100 The Ground-state wavefunction is taken to have a fully occupied / unoccupied bath */
  /* 0000 1000 mu is forced to be <n> U/2 - this also sets NmuSym to NmuSym & 10 */
  uint32_t NGCluster;           /* Number of poles in the resulting cluster green function */
  uint32_t NGBath;              /* Number of poles in the resulting bath green function */
  double Mix;                   /* mixing of old into new function */
  double U;                     /* valence - valence Coulomb repulsion */
  uint32_t NBath;               /* Number of Bath states to include */
  uint32_t Spectra;             /* 0000 no spectra will be calculated at the end of the DMFT loop */
  /* 0001 Magnetic susceptibility will be calculated (nup-ndn) */
  /* 0010 occupation will be calculated (nup+ndn) */
  /* 0100 XAS will be calculated with NQ values for core - valence interaction */
  /* 1000 core PES will be calculated with NQ values for core - valence interaction */
  uint32_t NQ;                  /* number of XAS spectra */
  double *Q;                    /* NQ values of core - valence interaction strengths */
  uint32_t NE;                  /* number of points in the spectrum */
  double Emin;                  /* minimum energy */
  double Emax;                  /* maximum energy */
  double Gamma;                 /* Broadening only important for the Spec file */
  double zeroLancEnergy;        /* minimum energy difference for convergence of total energy to be different from zero */
  double zeroLancWavefunction;  /* minnimum pre-factor in wave-function to be different from zero */
  uint32_t NTriDiag;            /* number of tri-diagonalization steps in Lanczos */
  uint32_t NBitsKey;            /* 2^NBitsKey == lenght of hash key */
}
DMFTSettings;



typedef struct DMFTXASSettings
{
  double U;                     /* valence - valence Coulomb repulsion */
  double Q;                     /* valence - core    Coulomb repulsion */
  char   BathFilename[255];     /* File name of converged Andersen impurity bath parameters */
  uint32_t NBathMax;            /* Max number of Bath states to include */
  uint32_t NE;                  /* number of points in the spectrum */
  double Emin;                  /* minimum energy */
  double Emax;                  /* maximum energy */
  double Gamma;                 /* Lorentzian lifetime */
  double zeroLancEnergy;        /* minimum energy difference for convergence of total energy to be different from zero */
  double zeroLancWavefunction;  /* minnimum pre-factor in wave-function to be different from zero */
  uint32_t NTriDiag;            /* number of tri-diagonalization steps in Lanczos */
  uint32_t NBitsKey;            /* 2^NBitsKey == lenght of hash key */
  char   SpecFile[256];         /* file name where the spectra are written */
}
DMFTXASSettings;

typedef struct DimerizedSpinChainSettings
{
  double JA;                    /* nearest neighbor exchange integral A */
  double JB;                    /* nearest neighbor exchange integral B */
  uint32_t Ndimer;              /* number of dimers in the chain */
  uint32_t NPsi;                /* number of low energy excited states to include in order to do temperature dependence */
  uint32_t NE;                  /* number of points in the spectrum */
  double Emin;                  /* minimum energy */
  double Emax;                  /* maximum energy */
  double Gamma;                 /* Lorentzian lifetime */
  double zeroLancEnergy;        /* minimum energy difference for convergence of total energy to be different from zero */
  double zeroLancWavefunction;  /* minnimum pre-factor in wave-function to be different from zero */
  uint32_t NTriDiag;            /* number of tri-diagonalization steps in Lanczos */
  uint32_t NBitsKey;            /* 2^NBitsKey == lenght of hash key */
  char   SpecFile[256];         /* file name where the spectra are written */
}
DimerizedSpinChainSettings;

typedef struct HeisenbergModelSettings
{
  double J;                     /* nearest neighbor exchange integral */
  CrystalStructureSuperCellType Cluster;  /* definition of the cluster */
  uint32_t NPsi;                /* number of low energy excited states to include in order to do temperature dependence */
  uint32_t NE;                  /* number of points in the spectrum */
  double Emin;                  /* minimum energy */
  double Emax;                  /* maximum energy */
  double Gamma;                 /* Lorentzian lifetime */
  double zeroLancEnergy;        /* minimum energy difference for convergence of total energy to be different from zero */
  double zeroLancWavefunction;  /* minnimum pre-factor in wave-function to be different from zero */
  uint32_t NTriDiag;            /* number of tri-diagonalization steps in Lanczos */
  uint32_t NBitsKey;            /* 2^NBitsKey == lenght of hash key */
  uint32_t LocalFourierBasis;   /* 1 use local basis for H and T */
  /* 2 use local basis for H and fourier for T */
  /* 3 use fourier basis for H and fourier for T */
  char   SpecFile[256];         /* file name where the spectra are written */
}
HeisenbergModelSettings;


typedef struct tJModelSettings
{
  double J;                     /* nearest neighbor exchange integral */
  double t;                     /* nearest neighbor hopping integral */
  uint32_t nElectron;           /* number of electrons (must be between zero and 2*Cluster.NCell) */
  CrystalStructureSuperCellType Cluster;  /* definition of the cluster */
  uint32_t NPsi;                /* number of low energy excited states to include in order to do temperature dependence */
  uint32_t NE;                  /* number of points in the spectrum */
  double Emin;                  /* minimum energy */
  double Emax;                  /* maximum energy */
  double Gamma;                 /* Lorentzian lifetime */
  double zeroLancEnergy;        /* minimum energy difference for convergence of total energy to be different from zero */
  double zeroLancWavefunction;  /* minnimum pre-factor in wave-function to be different from zero */
  uint32_t NTriDiag;            /* number of tri-diagonalization steps in Lanczos */
  uint32_t NBitsKey;            /* 2^NBitsKey == lenght of hash key */
  uint32_t LocalFourierBasis;   /* 1 use local basis for H and T */
  /* 2 use local basis for H and fourier for T */
  /* 3 use fourier basis for H and fourier for T */
  char   SpecFile[256];         /* file name where the spectra are written */
}
tJModelSettings;

typedef struct SpinLaderSettings
{
  uint32_t NSite;
  double Jleg;
  double Jrung;
  double Jcyc;
  uint32_t NPsi;
  uint32_t NE;
  double Emin;
  double Emax;
  double Gamma;
  double zeroLancEnergy;
  double zeroLancWavefunction;
  uint32_t NTriDiag;
  uint32_t NBitsKey;
  char   SpecFile[256];         /* file name where the spectra are written */
  int32_t AtomPos[100][3];      /* position of atoms, not read in but calculated from above info */
  double gUnitCell[3][3];       /* reciprocal unit cell */
  double kVector[100][3];       /* k-points wihtin super cell */
}
SpinLaderSettings;

typedef struct HubbardModelSettings
{
  TightBindingDef TightBinding; /* setup of hopping cell and sites */
  CrystalStructureSuperCellType SuperCell;  /* definition of super cell */
  TightBindingDef TightBindingSuperCell;  /* setup of hopping cell and sites */
  OperatorType U;               /* correlation operator */
  OperatorType USuperCell;      /* correlation operator within the super cell */
  uint32_t ***NElectron;        /* Nelectron[i][j][k] describes the number of electrons on atom i shell j fermion k */
  uint32_t Nq[3];               /* number of q (k) points to integrate the twisted boundary conditions over for a,b, and c direction */
  uint32_t NPsi;                /* number of low energy wave functions to calculate in order to do T dependence */
  double DeltaEMax;             /* maximum energy difference between ground-state and excited state */
  uint32_t NE;
  double Emin;
  double Emax;
  double Gamma;
  double zeroLancEnergy;
  double zeroLancWavefunction;
  uint32_t NTriDiag;
  uint32_t NBitsKey;            /* 2^NBitsKey == lenght of hash key */
  char   SpecFile[256];         /* file name where the spectra are written */
}
HubbardModelSettings;

/* Defines a structure that can store the Stuttgart CTRL file for LMTO and NMTO calculations
 * The controle file is structured into categories
 *
 * The StuttgartCTRL file defines a struct for each category in the CTRL file.
 *   Categories defined:
 *    HEADER VERS IO DIM STRUC OPTIONS CLASS SITE SYMGRP SCALE STR BZ EWALD DOS SYML START HARTREE PLOT CHARGE FINDES SCELL RHOFIT NMTO LOCALS
 */

typedef struct StuttgartCTRL
{
  char HEADER[1024];
  char VERS[64];
  struct
  {
    Bool HELP;
    Bool WKP;
    Bool IACTIV;
    int VERBOS;
    int ERRTOL;
    char OUTPUT[MAXStrLength];
    char ERR[MAXStrLength];
  }
  IO;
  struct
  {
    int NBAS;
    int NCLASS;
    int NL;
    int LDIM;
    int IDIM;
    int NSYMOP;
    int NKP;
  }
  DIM;
  struct
  {
    double ALAT;
    double PLAT[3][3];
    Bool FIXLAT;
  }
  STRUC;
  struct
  {
    Bool REL;
    Bool CCOR;
    Bool NONLOC;
    Bool CORDRD;
    Bool CARTESIAN;
    Bool CHARGE;
    Bool FATBAND;
    Bool AFM;
    Bool FS;
    Bool OVLCOR;
    Bool WRIBAS;
    Bool SEWALD;
    Bool LOCALCS;
    int NRXC;
    int NITATOM;
    int NSPIN;
    int NRMIX;
    char Q[16];
  }
  OPTIONS;
  struct
  {
    int N;           /* not in control file - number of clases */
    Bool *NEWNAM;
    char **ATOM;
    int *Z;
    double *R;
    int *LMX;
    int **IDXDN;
    int **CONF;
    int **IDMOD;
  }
  CLASS;
  struct
  {
    int N;              /* not in control file - number of sites */
    int Nactive;        /* not in control file - number of active sites */
    int *ActiveAtomMap; /* not in control file - maping of site to active site index */
    int *ICLASS;        /* not in control file - atom belongs to CLASS ICLASS */
    char **ATOM;
    double **POS;
    double **X;
    int **IDXDN_M;
  }
  SITE;
  struct
  {
    Bool USESYM;
    int NGEN;
    int IORIGIN;
    char SPCGRP[64];
    char GENGRP[64][32];
  }
  SYMGRP;
  struct
  {
    Bool SCLWSR;
    double FACVOL;
    double OMMAX1[3];
    double OMMAX2[3];
    double GAMMA;
  }
  SCALE;
  struct
  {
    Bool DOWATS;
    Bool DONALP;
    Bool NOCALC;
    int N;            /* not in control file - number of ATOMS for wich SIGMA is defined */
    int IALPHA;
    int LMAXW;
    int NDIMIN;
    int NTRANS;
    int LDN;
    int MDN;
    int JBASDN;
    double KAPPA2;
    double KWxx2;
    double ADOT;
    double DELTR;
    double RMAXS;
    double EPS;
    double **SIGMA;
    char **ATOM;
  }
  STR;
  struct
  {
    Bool TETRA;
    Bool METAL;
    int NKABC[3];
    int MAXKP;
    int N;
    int RANGE;
    int NPTS;
    double W;
    double TOL;
  }
  BZ;
  struct
  {
    int NKDMX;
    double AS;
    double TOL;
  }
  EWALD;
  struct
  {
    int NOPTS;
    double EMIN;
    double EMAX;
  }
  DOS;
  struct
  {
    int N;         /* not in control file - number of linesegments */
    int *NQ;
    double **Q1;
    double **Q2;
    char **LAB1;
    char **LAB2;
  }
  SYML;
  struct
  {
    Bool BROY;
    Bool FREE;
    Bool BEGMOM;
    Bool CNTROL;
    int N;            /* not in control file - number of ATOMS for settings are defined */
    int NIT;
    int NMIX;
    double WC;
    double BETA;
    double CNVG;
    double CNVGET;
    double EFERMI;
    double VMTZ;
    char **ATOM;
    double **P;
    double **Q;
    double **enu;
    double **c;
    double **sqrdel;
    double **p;
    double **gamma;
  }
  START;
  struct
  {
    Bool BEGATOM;
    int LT1;
    int LT2;
    int LT3;
  }
  HARTREE;
  struct
  {
    int NDELR1;
    int NDELR2;
    int NDELR3;
    int FORMAT;
    double ORIGIN[3];
    double R1[3];
    double R2[3];
    double R3[3];
    double RMIN;
    double RMAX;
    double TMIN;
    double TMAX;
  }
  PLOT;
  struct
  {
    Bool LMTODAT;
    Bool CHARWIN;
    Bool ELF;
    Bool ADDCOR;
    Bool SPINDENS;
    double EMIN;
    double EMAX;
  }
  CHARGE;
  struct
  {
    int MAXPT;
    int NRXYZ[3];
    double RMINES;
    double RMAXES;
  }
  FINDES;
  struct
  {
    Bool EQUIV;
    double PLAT[3][3];
  }
  SCELL;
  struct
  {
    Bool FIT;
    int N;         /* not in control file - number of ATOMS - should be the same as in CLASS */
    int OUTPUT;
    double KAPPA2;
    double KWxx2;
    double RMAXS;
    char **ATOM;
    double **SIGMA;
    int *LMXRHO;
  }
  RHOFIT;
  struct
  {
    int NENU;
    int PLOTATO;
    int ORBITAL;
    double RMAXHR;
    double *ENU;
    Bool eV;
    Bool ORTHO;
  }
  NMTO;
  struct
  {
    char Data_as[32];
    int Number_of_rotated_atoms;
    int *ATOM;
    int *L;
    double **LOCAL_Z;
    double **LOCAL_X;
  }
  LOCALCS;
}
StuttgartCTRL;

typedef struct StuttgartHAMRelement
{
  int IndexAtom2;
  int IndexAtom1;
  int TransVecs[3];
  int NlabelAtom1;
  int NlabelAtom2;
  double Dist;
  double ConnectingVector[3];
  char LabelAtom1[16][16];
  char LabelAtom2[16][16];
  char Atom1[16];
  char Atom2[16];
  double H[16][16];
}
StuttgartHAMRelement;

typedef struct StuttgartHAMR
{
  Bool eV;
  int Nelement;
  StuttgartHAMRelement *element;
}
StuttgartHAMR;

typedef struct DresdenFPLO
{
  char Compound[16];
  struct
  {
    double LattConst[3];
    double Angles[3];
    double LattVec[3][3];
    double gLattVec[3][3];
  }
  Lattice;

  struct
  {
    int Natom;
    char **Atom;
    double **Pos;
  }
  Atoms;

  struct
  {
    int Nwandef;
    char **Name;
    int NActAtom;
    int *ActMap;
    int *Site;
    int **Orbital;     // n,l,m
    double **PosCenter;
    int *Nhop;
    double ***HopVec;  // Translating+Connecting, as direct output of fplo.
    double **HopVal;
  }
  Wannier;

}
DresdenFPLO;



// Obsolete structures that are replaced by new versions, use is discuraged and need to be faded out



/*
 *     ( A[0][0]=a1   | A[0][1]=b1   | 0            | 0          )
 *   A=( A[1][0]=b1^+ | A[1][1]=a2   | A[1][2]=b2   | 0          )
 *     ( 0            | A[2][1]=b2^+ | A[2][2]=a3   | A[2][3]=b3 )
 *     ( 0            | 0            | A[3][2]=b3^+ | A[3][3]=a4 )
 *
 *   with ai and bi matrices
 *
 *   (A B)[i][j] = \sum_k A[i][k] . B[k][j]
 *
 *   !!!! warning, part of the routines include a0 and b0 as ofsets, but not all routines do !!!!!!
 */
typedef struct BlokTriDiagonalMatrixType
{
  char   Name[MAXStrLength];    /* A pointer to a name string that can be used to lable functions, good for printing etc. */
  uint32_t N;                   /* number of diagonal blocks */
  MatrixType *a;                /* diagonal elements of matrix (length N)       elements a[i]=M[i,i]                  of the full matrix */
  MatrixType *b;                /* off diagonal elements of matrix (length N-1) elements b[i]=M[i,i+1]=M[i+1,i]^{dag} of the full matrix */
  double mu;                    /* chemical potential if transformed to Greens function */
}
BlokTriDiagonalMatrixType;

/*
 *     ( A[0][0]=a1   | A[0][1]=b1   | A[1][2]=b2   | A[2][3]=b3 )
 *   A=( A[1][0]=b1^+ | A[1][1]=a2   | 0            | 0          )
 *     ( A[2][1]=b2^+ | 0            | A[2][2]=a3   | 0          )
 *     ( A[3][2]=b3^+ | 0            | 0            | A[3][3]=a4 )
 *
 *   with a and b matrices
 *
 *   (A B)[i][j] = \sum_k A[i][k] . B[k][j]
 *
 *   !!!! warning, part of the routines include a0 and b0 as ofsets, but not all routines do !!!!!!
 */
typedef struct BlokAndersonMatrixType
{
  char   Name[MAXStrLength];    /* A pointer to a name string that can be used to lable functions, good for printing etc. */
  uint32_t N;                   /* number of diagonal blocks */
  MatrixType *a;                /* diagonal elements of matrix (length N)       elements a[i]=M[i,i]                  of the full matrix */
  MatrixType *b;                /* off diagonal elements of matrix (length N-1) elements b[i]=M[0,i+1]=M[i+1,0]^{dag} of the full matrix */
  double mu;                    /* chemical potential if transformed to Greens function */
}
BlokAndersonMatrixType;

/*
 *     ( A[0][0]=acv0  | A[0][1]=bv0   | 0             | 0           | A[0][4]=bc0   | 0             | 0           )
 *   A=( A[1][0]=bv0^+ | A[1][1]=av1   | A[1][2]=bv1   | 0           | 0             | 0             | 0           )
 *     ( 0             | A[2][1]=bv1^+ | A[2][2]=av2   | A[2][3]=bv2 | 0             | 0             | 0           )
 *     ( 0             | 0             | A[3][2]=bv2^+ | A[3][3]=av3 | 0             | 0             | 0           )
 *     ( A[4][0]=bc0^+ | 0             | 0             | 0           | A[4][4]=ac1   | A[4][5]=bc1   | 0           )
 *     ( 0             | 0             | 0             | 0           | A[5][4]=bc1^+ | A[5][5]=ac2   | A[5][6]=bc2 )
 *     ( 0             | 0             | 0             | 0           | 0             | A[6][5]=bc2^+ | A[6][6]=ac3 )
 *
 *   acv0 = (a0val + a0con)/2
 *
 *   (A B)[i][j] = \sum_k A[i][k] . B[k][j]
 *
 *   Note that the block A[1;;3][1;;3] and the block A[4;;6][4;;6] are normal blok-tri-diagonal matrices
 */
typedef struct BlokDoubleTriDiagonalMatrixType
{
  char   Name[MAXStrLength];    /* A pointer to a name string that can be used to lable functions, good for printing etc. */
  double mu;                    /* border if a number belongs to val or con */
  double Epsilon;               /* values in a range of [Ef-1/2 Epsilon , Ef+1/2 Epsilon] contribute half to val and half to con */
  MatrixType a0;                /* storing a possible real shift if a continued fraction notation is used */
  MatrixType b0;                /* storing a possilbe over-all scaling if a continued fraction notation is used */
  BlokTriDiagonalMatrixType val;
  BlokTriDiagonalMatrixType con;
}
BlokDoubleTriDiagonalMatrixType;



/*
 *
 */
typedef struct WaveFunctionPerturbationType
{
  uint32_t NV;
  uint32_t order;
  WaveFunctionType **psiOrder;
}
WaveFunctionPerturbationType;

/*
 *
 */
typedef struct OperatorPerturbationType
{
  uint32_t NV;
  OperatorType O;
  OperatorType *V;
  BlokTriDiagonalMatrixType g;
}
OperatorPerturbationType;

/*
 *
 */
typedef struct TemperatureSettings
{
  uint32_t NFolders;            /* Number of folders needed to average over */
  char   FolderName[100][256];  /* Name of the folder */
  char   FileName[100][256];    /* General root of the file name within the corresponding folder */
  uint32_t LowestFileNumber[100]; /* Lowest file number */
  uint32_t HighestFileNumber[100];  /* Highest file number */
  char   EnergyFileName[100][256];  /* File name where the energy can be found */
  double Degeneracy[100];       /* Multiplication factor for these spectra */
  uint32_t NTemperatures;       /* number of temperatures to average over */
  double Temperatures[256];     /* Temperatures used */
  uint32_t Nq[3];               /* Number of twisted boundary phases */
}
TemperatureSettings;
