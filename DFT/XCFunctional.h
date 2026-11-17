/*
 * XCFunctional.h
 *
 * Exchange-correlation functionals of the DFT/HF add-on (stage A3), see XCFunctional.cpp.
 */

#ifndef DFT_XCFUNCTIONAL_H_
#define DFT_XCFUNCTIONAL_H_

#define XCLDA     0     // Slater + PW92
#define XCLDAVWN  1     // Slater + VWN5 (NIST LDA)
#define XCRLDAMV  2     // Slater with MacDonald-Vosko + VWN5 (NIST RLDA)
#define XCPBE     3     // PBE
#define XCPBEVWN  4     // PBE with VWN5 under the correlation gradient correction

int XCIsGGA(int method);
void XCVWN(DFTReal rs, DFTReal A, DFTReal x0, DFTReal b, DFTReal c, DFTReal *eps, DFTReal *deps);
void XCPW92(DFTReal rs, DFTReal *eps, DFTReal *deps);
void XCMacDonaldVosko(DFTReal beta, DFTReal *Phi, DFTReal *dPhi);
int XCEvaluate(int method, double c, DFTReal n, DFTReal sigma, DFTReal cutoff, DFTReal *e, DFTReal *en, DFTReal *es);
int XCMethodFromName(const char *name);

#endif /* DFT_XCFUNCTIONAL_H_ */
