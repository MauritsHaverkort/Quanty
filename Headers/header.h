
uint32_t GlobalTimeOn=0;
/*      1111 1111 1111 1111 - 0xFFFF - print all timings
 *      0000 0000 0000 0001 - 0x0001 - main specific
 *      0000 0000 0000 0010 - 0x0002 - main specific
 *      0000 0000 0000 0100 - 0x0004 - main specific
 *      0000 0000 0000 1000 - 0x0008 - main specific
 *      0000 0000 0001 0000 - 0x0010 -
 *      0000 0000 0010 0000 - 0x0020 -
 *      0000 0000 0100 0000 - 0x0040 -
 *      0000 0000 1000 0000 - 0x0080 -
 */

uint32_t GlobalVerbosity=0;
/*      1111 1111 1111 1111 - 0xFFFF - print all statments on what the code does
 *      0000 0000 0000 0001 - 0x0001 - main specific
 *      0000 0000 0000 0010 - 0x0002 - main specific
 *      0000 0000 0000 0100 - 0x0004 - main specific
 *      0000 0000 0000 1000 - 0x0008 - main specific
 *      0000 0000 0001 0000 - 0x0010 - Lanczos function start
 *      0000 0000 0010 0000 - 0x0020 - Lanczos function converge info
 *      0000 0000 0100 0000 - 0x0040 - 
 *      0000 0000 1000 0000 - 0x0080 -
 *
 *      1111 1111 1111 1111 - 0xFFFF0000 - print all debug statements
 *      0000 0000 0000 0001 - 0x00010000 - print start of each subroutine
 *      0000 0000 0000 0010 - 0x00020000 -
 *      0000 0000 0000 0100 - 0x00040000 -
 *      0000 0000 0000 1000 - 0x00080000 -
 *      0000 0000 0001 0000 - 0x00100000 -
 *      0000 0000 0010 0000 - 0x00200000 -
 *      0000 0000 0100 0000 - 0x00400000 -
 *      0000 0000 1000 0000 - 0x00800000 -
 *
 *      0000 0001 0000 0000 - 0x01000000 - used for debuging of functions during development
 */

#define PrintDebugInfo


#include "core.h"
#include "io.h"
#include "core_basic_math.h"


/**********************************************************/
/*****************       blas/lapack      *****************/
/**********************************************************/
#ifdef BLASunderscore
#define daxpb daxpb_
#define dcopy dcopy_
#define ddot ddot_
#defien dgemm dgemm_
#define dgemv dgemv_
#define dgeqrf dgeqrf_
#define dgesvd dgesvd_
#define dgetrf dgetrf_
#define dgetri dgetri_
#define dsyrk dsyrk_
#define dpotrf dpotrf_
#define dstev dstev_
#define dscal dscal_
#define dstevd dstevd_
#define dsyev dsyev_
#define dsyevd dsyevd_
#define dsygv dsygv_
#define dsygvd dsygvd_
#define dtrsm dtrsm_
#define zaxpb zaxpb_
#define zgemm zgemm_
#define zgemv zgemv_
#define zgetrf zgetrf_
#define zgetri zgetri_
#define zherk zherk_
#define zheev zheev_
#define zheevd zheevd_
#define zgesvd zgesvd_
#define zhegv zhegv_
#define zhegvd zhegvd_
#define zdotc zdotc_
#define zscal zscal_
#endif


int    daxpy();
int    dcopy();
double ddot();
int    dgemm();
int    dgemv();
int    dgeqrf();
int    dgesvd();
void   dgetrf();
void   dgetri();
int    dpotrf();
int    dscal();
int    dsyrk();
int    dstev();
int    dstevd();
int    dsyev();
int    dsyevd();
int    dsygv();
int    dsygvd();
int    dtrsm();
int    zscal();
int    zaxpy();
int    zgemm();
int    zgemv();
void   zgetrf();
void   zgetri();
int    zherk();
int    zheev();
int    zheevd();
int    zgesvd();
int    zhegv();
int    zhegvd();
void   zdotc();

/**********************************************************/
/********      Hubbard Model/CellFunctions.c      *********/
/**********************************************************/
uint32_t AddAditionalInformationToSuperCell(CrystalStructureSuperCellType * A);
uint32_t CopyCrystalStructure(CrystalStructureType A, CrystalStructureType * B);


/**********************************************************/
/********    Tight Binding/Tight_Binding_lib.c    *********/
/**********************************************************/
int    SetNFermionIndexTightBindingDef(TightBindingDef * A);
uint32_t CreateUSperCell(OperatorType U, OperatorType * USuperCell, CrystalStructureSuperCellType SuperCell);
int    CreateTightBindingSuperCell();
int CreateHamiltonianFromTightBindingDef();
int CreateHamiltonianFromTightBindingDefConstantPhaseInUnitCell();
int CreateCrystalStructureFromCrystalStructureSuperCell();
void RealEigenFunctionToCharacter();
void ComplexEigenFunctionToCharacter();
void EigenFunctionToCharacter();
void RealEigenFunctionToCharacterMultiSite();
void ComplexEigenFunctionToCharacterMultiSite();
void EigenFunctionToCharacterMultiSite();



/**********************************************************/
/********    DMFT.c    *********/
/**********************************************************/
uint32_t DMFTSpec();
