/* This file was automatically generated.  Do not edit! */
int main(int argc,const char *argv[]);
void PEStest();
uint32_t CreateRealWaveFunctionNew(WaveFunctionType *const psi,uint32_t const NFermions,uint32_t const nelectrons);
uint32_t CreateRealWaveFunctionShellOccupacy(WaveFunctionType *const psi,uint32_t const Nshell[],uint32_t const nelectrons[],uint32_t const arraylength);
uint32_t SpectraAddPlasmons(SpectraType *Spectra,double const Scale,double const EnergyShift,double const LorentzianFWHM,uint8_t maxOrder,double const IntrinsicExtrinsic);
uint32_t BroadenGenergydependent(SpectraType *const G,double const energy[],double const FWHM[],uint32_t const arraylength,double const deltaE);
uint32_t BroadenLenergydependent(SpectraType *const G,double const energy[],double const FWHM[],uint32_t const arraylength,double const deltaE);
uint32_t BroadenLshift(SpectraType *const G,double const FWHM,double const deltaE);
void BroadenL(SpectraType *const G,double const FWHM);
uint32_t BroadenMenergydependent(SpectraType *const G,double const energy[],double const xi0[],double const g[],uint32_t const arraylength);
uint32_t BroadenM(SpectraType *const G,double const xi0,double const g);
double Gamma(double const x);
void test();
