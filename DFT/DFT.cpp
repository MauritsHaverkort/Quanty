/*
 * DFT.cpp
 *
 * The DFT/HF add-on of Quanty, included from LuaInterface/main.cpp (single compilation
 * unit). See DFT.h. Planned files (design section 11):
 *     DiracRadial.cpp     nucleus models, DKB matrices, radial Dirac eigenproblem
 *     XCFunctional.cpp    exchange-correlation functionals
 *     AtomSCF.cpp         Dirac-Hartree-Fock and Dirac-Kohn-Sham atom
 *     LuaAtom.cpp         Lua bindings
 */

#include "DiracRadial.cpp"
#include "XCFunctional.cpp"
#include "AtomSCF.cpp"
#include "AtomHF.cpp"
#include "AtomBasis.cpp"
#include "NucleusMultipole.cpp"
#include "LuaDiracRadial.cpp"
#include "LuaAtom.cpp"
#include "LuaAtomicOperators.cpp"
#include "LuaAtomicBasis.cpp"
