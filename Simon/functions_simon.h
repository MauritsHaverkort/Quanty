template<typename T>void debug(const T debugMessage);
void debug();
template<typename T>void debug(const vector<T> &vec);
template<typename keyType, typename valType>void debug(const map<keyType, valType> &m);
void debug(const CompactMatrixType &M, const string &str);
void debug(const MatrixType &M, const string &str);
void debug(const CompactMatrixType &M);
void debug(const MatrixType &M);
static int LuaSimontest(lua_State * L);
