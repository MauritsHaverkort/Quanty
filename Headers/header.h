#include "core.h"
#include "io.h"
#include "core_basic_math.h"



/**********************************************************/
/*****************       blas/lapack      *****************/
/**********************************************************/
int dstev();
int dstevd();
int dsyev();
int dsyevd();
int zheev();
int zheevd();
int dgesvd();
int zgesvd();
int dsygv();
int zhegv();
int dsygvd();
int zhegvd();
/**********************************************************/
/********      Hubbard Model/CellFunctions.c      *********/
/**********************************************************/
uint32_t AddAditionalInformationToSuperCell(CrystalStructureSuperCellType * A);
uint32_t CopyCrystalStructure(CrystalStructureType A, CrystalStructureType *B);


/**********************************************************/
/********    Tight Binding/Tight_Binding_lib.c    *********/
/**********************************************************/
int SetNFermionIndexTightBindingDef(TightBindingDef *A);
uint32_t CreateUSperCell(OperatorType U, OperatorType *USuperCell, CrystalStructureSuperCellType SuperCell);
int CreateTightBindingSuperCell();



/**********************************************************/
/********    DMFT.c    *********/
/**********************************************************/
uint32_t DMFTSpec();









