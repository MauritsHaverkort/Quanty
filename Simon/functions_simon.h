template<typename T>void debug(T debugMessage);
void debug();
template<typename T>void debug(vector<T> vec);
template<typename keyType, typename valType>void debug(const map<keyType, valType> &m);
void debug(CompactMatrixType M, string str);
void debug(MatrixType M, string str);
void debug(CompactMatrixType M);
void debug(MatrixType M);
static int LuaSimontest(lua_State * L);
