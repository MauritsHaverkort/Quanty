
uint32_t LuaASCIIStreamPrint(lua_State *L,int narg,FILE *stream);
char *QuantyVariableToValidFilename(const char *str);
uint32_t LuaPutASCII(lua_State *L, const char * const DirName, const char * const VariableName);
uint32_t LuaPut(lua_State *L, const char * const FileName, const char * const VariableName, const char * const FileType);
uint32_t LuaGetASCII(lua_State *L, const char * const DirName, const char * const VariableName);
uint32_t LuaGet(lua_State *L, const char * const FileName, const char * const VariableName);
uint32_t QStringIsValidLuaVariableName(const char * const c);
