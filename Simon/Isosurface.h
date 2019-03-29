void lua_push(lua_State * L, const doubleVecValPair &pushkin);
double vecNorm(const vector<double> &v);
double multicenterMultipoleIntegral(Interpolation::InterpolatingFunction &f1, Interpolation::InterpolatingFunction &f2, int k, uint32_t NknotsPerInterval);
vector<double> interpolationRanges(const vector<Interpolation::InterpolatingFunction*> &funs);
vector<Complex16> realspaceBasisValues(const vector<double> &v, const vector<string> &orbs, const vector<vector<double> > &atomPos, const vector<double> &interRanges,const vector<Interpolation::InterpolatingFunction*> &funs, const vector<vector<uint16_t> > indices, const double step);
CompactMatrixType densityMatrixNoSpin(const CompactMatrixType &rho);
vector<doubleVecValPair> realSpaceWavefunction(const CompactMatrixType &rho, const vector<Interpolation::InterpolatingFunction*> &funs, const vector<string> &orbs,const vector<vector<double> > &atomPos, const uint gridpoints);
template<typename keyType, typename valType>void lua_pushMapKeyValPair(lua_State * L, const map<keyType, valType> &pushkin);
static int LuaRealspaceWavefunction(lua_State * L);
