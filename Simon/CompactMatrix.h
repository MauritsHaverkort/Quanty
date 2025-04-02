bool isRowZero(const CompactMatrixType &M, const uint32_t row, const double precision);
bool isRowZero(const CompactMatrixType &M, const uint32_t row);
bool isColZero(const CompactMatrixType &M, const uint32_t col, const double precision);
bool isColZero(const CompactMatrixType &M, const uint32_t col);
bool isZeroMatrix(const CompactMatrixType &M, const double zero);
void FreeCompactMatrixMatrix(vector<vector<CompactMatrixType> > &matMat);
void FreeMatrixMatrix(vector<vector<MatrixType> > &matMat);
void CompactMatrixCutDimensions(CompactMatrixType &M);
CompactMatrixType SubspaceRotationMatrix(const uint32_t dim1, const CompactMatrixType &T, const uint32_t dim2);
CompactMatrixType MatrixTypeToCompactMatrixType(const MatrixType &M);
CompactMatrixType removeSpin(const CompactMatrixType &M);
static int LuaRemoveSpin(lua_State *L);
CompactMatrixType Submatrix(const CompactMatrixType &M, const vector<uint> &rows, const vector<uint> &cols);
double RealMatrixTrace(const CompactMatrixType *M);
Complex16 ComplexMatrixTrace(const CompactMatrixType *M);
CompactMatrixType GSOrthogonalizeWithNonOrthoBasis(const CompactMatrixType &in, const CompactMatrixType &Overlap);
static int LuaOrthogonalizeWithNonOrthoBasis(lua_State *L);
static int LuaMatrixTrace(lua_State *L);

