/*
 * FourierBessel.h
 *
 * Fourier-Bessel transforms of the radial functions (DFT stage A5d), see FourierBessel.cpp.
 */

#ifndef DFT_FOURIERBESSEL_H_
#define DFT_FOURIERBESSEL_H_

int FourierBesselTransform(const RadialBasisType *b, int relativistic, int kappa, int l, double c, const double *coef,
                           uint32_t nk, const double *k, double *Pk, double *Qk);

#endif /* DFT_FOURIERBESSEL_H_ */
