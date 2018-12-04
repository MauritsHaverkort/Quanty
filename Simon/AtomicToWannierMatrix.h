string cutLeadingSpaces(string in);
string cutTrailingSpaces(string in);
string cutSpaces(string in);
bool atomicStateSorter(string s1, string s2);
vector<map<string, double> > getWF_coefficients_Filecontent(const string &filename, vector<string> &orbitalNames);
CompactMatrixType AtomicToWannierMatrix(const string &filename, vector<string> &orbitalNames);
bool isFileAcccessible(const string &name);
static int LuaAtomicToWannierMatrix(lua_State * L);
