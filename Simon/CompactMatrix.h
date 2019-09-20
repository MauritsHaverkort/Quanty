// mwh these functions were already defined 
//CompactMatrixType BlockTriDiagonalMatrixToCompactMatrix(const BlockTriDiagonalMatrixType &M, const double zero);
//CompactMatrixType BlockTriDiagonalMatrixToCompactMatrix(const BlockTriDiagonalMatrixType &M);
bool isRowZero(const CompactMatrixType &M, const uint row, const double precision);
bool isRowZero(const CompactMatrixType &M, const uint row);
bool isColZero(const CompactMatrixType &M, const uint col, const double precision);
bool isColZero(const CompactMatrixType &M, const uint col);
bool isZeroMatrix(const CompactMatrixType &M, const double zero);
void FreeCompactMatrixMatrix(vector<vector<CompactMatrixType> > &matMat);
void FreeMatrixMatrix(vector<vector<MatrixType> > &matMat);
void CompactMatrixCutDimensions(CompactMatrixType &M);
CompactMatrixType SubspaceRotationMatrix(const uint dim1, const CompactMatrixType &T, const uint dim2);
CompactMatrixType MatrixTypeToCompactMatrixType(const MatrixType &M);
CompactMatrixType removeSpin(const CompactMatrixType &M);
static int LuaRemoveSpin(lua_State *L);
CompactMatrixType Submatrix(const CompactMatrixType &M, const vector<uint> &rows, const vector<uint> &cols);
