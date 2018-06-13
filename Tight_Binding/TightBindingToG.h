double DOSFunctionBethe(double w);
double DOSFunctionBlock(double w);
double DOSFunctionChain(double w);
uint32_t DOSFunctionToGTriDiagonalRepresentation(TriDiagonalMatrixType * G, uint32_t N, double (*DOS);
uint32_t DOSFunctionToGTriDiagonalRepresentationNonOrthogonal(TriDiagonalMatrixType * G, uint32_t N,double (*DOS);
uint32_t NearestNeighborHyperCubicTBToGTriDiagonalRepresentation(uint32_t Dim, uint32_t NTriDiag,TriDiagonalMatrixType * const G, uint32_t Nk);
uint32_t NearestNeighborHyperCubicTBToGTriDiagonalRepresentationNonOrthogonal(uint32_t Dim, uint32_t NTriDiag,TriDiagonalMatrixType * const G,uint32_t Nk);
uint32_t NearestNeighborTriangularTBToGTriDiagonalRepresentation(uint32_t NTriDiag, TriDiagonalMatrixType * const G,uint32_t Nk);
uint32_t NearestNeighborTriangularTBToGSpectraRepresentation(SpectraType * const G, uint32_t Nk);
