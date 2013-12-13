#include "core.h"
#include "io.h"
#include "core_basic_math.h"



/**********************************************************/
/*****************       blas/lapack      *****************/
/**********************************************************/
#ifdef BLASunderscore
#define dcopy dcopy_
#define dstev dstev_
#define dstevd dstevd_
#define dsyev dsyev_
#define dsyevd dsyevd_
#define zheev zheev_
#define zheevd zheevd_
#define dgesvd dgesvd_
#define zgesvd zgesvd_
#define dsygv dsygv_
#define zhegv zhegv_
#define dsygvd dsygvd_
#define zhegvd zhegvd_
#define zdotc zdotc_
#endif

int    dstev();
int    dstevd();
int    dsyev();
int    dsyevd();
int    zheev();
int    zheevd();
int    dgesvd();
int    zgesvd();
int    dsygv();
int    zhegv();
int    dsygvd();
int    zhegvd();
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



/**********************************************************/
/********    DMFT.c    *********/
/**********************************************************/
uint32_t DMFTSpec();
