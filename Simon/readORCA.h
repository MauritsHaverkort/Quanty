void lua_push(lua_State * L, const ORCAObject &pushkin);
ORCAObject FileReadORCA(const string &filename);
static int LuaFileReadORCA(lua_State* L);
