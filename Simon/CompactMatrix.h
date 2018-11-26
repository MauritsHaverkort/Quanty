CompactMatrixType BlockTriDiagonalMatrixToCompactMatrix(const BlockTriDiagonalMatrixType &M, const double zero);
CompactMatrixType BlockTriDiagonalMatrixToCompactMatrix(const BlockTriDiagonalMatrixType &M);
bool isRowZero(const CompactMatrixType &M, const uint row, const double precision);
bool isRowZero(const CompactMatrixType &M, const uint row);
bool isColZero(const CompactMatrixType &M, const uint col, const double precision);
bool isColZero(const CompactMatrixType &M, const uint col);
bool isZeroMatrix(const CompactMatrixType &M, const double zero);
void FreeCompactMatrixMatrix(vector<vector<CompactMatrixType> > matMat);
