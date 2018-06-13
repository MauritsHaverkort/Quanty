uint32_t AddAditionalInformationToSuperCell(CrystalStructureSuperCellType * A);
uint32_t CopyCrystalStructure(CrystalStructureType A, CrystalStructureType * B);
uint32_t ShortesConnectingVectorPeriodicSuperCell(const int32_t A[], const int32_t B[],const CrystalStructureSuperCellType SC, int32_t C[]);
uint32_t FindCellPosInCrystalStructureSuperCellType(CrystalStructureSuperCellType A, int32_t vec[3], int *pos);
uint32_t SpectraSiteFourierTransform(const CrystalStructureSuperCellType SC, SpectraType * Spectra);
uint32_t OperatorSiteFourierTransform(const CrystalStructureSuperCellType SC, const OperatorType Oi, OperatorType * Ok,uint32_t(*QDetNotOrdered);
