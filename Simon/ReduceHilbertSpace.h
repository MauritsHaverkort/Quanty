TightBindingDef reduceHilbertSpaceForTightBinding(const TightBindingDef &tb, const vector<AtomType> &consideredAtoms, const vector<AtomType> &unchangedAtoms, const uint filledLigands, const uint emptyLigands);
OperatorType reduceHilbertSpaceForOperator(OperatorType &hamil, const uint unchangedStates, const uint filledLigands, const uint emptyLigands);
CompactMatrixType reducedHilbertSpaceMatrix(const CompactMatrixType &H, const CompactMatrixType &ligandRotMat, CompactMatrixType &T, const uint unchangedStates, const uint ligandNum);
CompactMatrixType reducedHilbertSpaceMatrixFilledAndEmpty(const CompactMatrixType &HRedFilled, const CompactMatrixType &HRedEmpty, const uint unchangedStates);
CompactMatrixType reducedHilbertSpaceRotMatFilledAndEmpty(const CompactMatrixType &Tfilled, const CompactMatrixType &Tempty, const uint unchangedStates);
CompactMatrixType reduceHilbertSpaceForMatrix(const CompactMatrixType &H, CompactMatrixType &T, const uint unchangedStates, const uint filledLigands, const uint emptyLigands);
static int LuaReduceHilbertSpace(lua_State * L);
