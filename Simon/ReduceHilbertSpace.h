TightBindingDef reduceHilbertSpaceForTightBinding(const TightBindingDef &tb, const vector<AtomType> &consideredAtoms, const vector<AtomType> &unchangedAtoms, const uint filledLigands, const uint emptyLigands);
OperatorType reduceHilbertSpaceForOperator(OperatorType &hamil, const uint unchangedStates, const uint filledLigands, const uint emptyLigands);
CompactMatrixType reducedHilbertSpaceMatrix(const CompactMatrixType &H, const CompactMatrixType &ligandRotMat, const uint unchangedStates, const uint ligandNum);
CompactMatrixType reduceHilbertSpaceForMatrix(const CompactMatrixType &H, const uint unchangedStates, const uint filledLigands, const uint emptyLigands);
static int LuaReduceHilbertSpace(lua_State * L);
