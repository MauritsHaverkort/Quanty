/*
 * DFT.h
 *
 * The DFT/HF add-on of Quanty: relativistic atoms on B-splines (stage A), later the
 * solid. Design: QuantyDFT_StageA_DiracAtom.md in Doc/Agents.
 *
 * Conventions (design section 0 and 0.1):
 *  - Hartree atomic units; the speed of light c is an argument of every routine that
 *    needs it, read by the Lua binding from PhysicalConstants.cAU.value;
 *  - double precision; DFTReal is the scalar of the LAPACK-free kernels, so that a
 *    long double build of them stays possible;
 *  - every tolerance is an argument, with its default in the Lua binding;
 *  - routines return 0 on success and 1 on error, allocate with malloc, do not print.
 */

#ifndef DFT_DFT_H_
#define DFT_DFT_H_

typedef double DFTReal;

#include "DiracRadial.h"
#include "XCFunctional.h"
#include "AtomSCF.h"
#include "AtomHF.h"
#include "AtomBasis.h"
#include "NucleusMultipole.h"
#include "FourierBessel.h"

#endif /* DFT_DFT_H_ */
