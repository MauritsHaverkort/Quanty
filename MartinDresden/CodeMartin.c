#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <unistd.h>
#include <stdarg.h>


#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#define omp_get_num_procs() 0
#define omp_get_max_threads() 0
#endif

#include "../Core/StructDef.c"
#include "../Headers/header.h"
#include "../Core/TimeStamp.c"
#include "../Core/Hash.c"
#include "../Core/BasicMath.c"
#include "../IO/Print.c"
#include "../IO/Read.c"
#include "../IO/BitMap.c"
#include "../Core/Error.c"
#include "../Core/WaveFunctionArithmetic.c"
#include "../Core/WaveFunctionPerturbationArithmetic.c"
#include "../Core/OperatorArithmetic.c"
#include "../Core/OperatorWaveFunction.c"
#include "../Core/OperatorWaveFunctionPerturbation.c"
#include "../Core/Lanczos.c"
#include "../Core/Spectral.c"
#include "../Core/DensityMatrix.c"
#include "../Core/Random.c"
#include "../Core/MatrixOperator.c"
#include "../Model/StandardOperators.c"
#include "../Hubbard_Model/CellFunctions.c"
#include "../Tight_Binding/Tight_Binding_lib.c"
#include "../Tight_Binding/TightBindingToG.c"
#include "../DMFT/SingleSiteFunctions.c"



void test()
{
  printf("!!!I am here!!!\n");
}

double Gamma(double const x){ // copied from Martin Sundermann and Martin Rotter edit by MS
	int i, gg=7, np=9;
	double const p[] = {0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7};
	double y,a=p[0];
  if(x<0.5){
		for(i=1;i<np;i++) a+=p[i]/(-x+(double)i);
		y=Pi/(sin(Pi*x)*sqrt(2.*Pi)*pow(-x+(double)gg+0.5,-x+0.5)*exp(x-(double)gg-0.5)*a);
  }else{
		for(i=1;i<np;i++) a+=p[i]/(x+(double)i-1.0);
		y=sqrt(2.*Pi)*pow(x+(double)gg-0.5,x-0.5)*exp(-x-(double)gg+0.5)*a;
	}
  return y;
}

uint32_t BroadenM(SpectraType *const G, double const xi0, double const g) { // copied from MS and Martin Rotter edit by MS
	if(xi0<3.0*G->dE) {
    printf("Warning in BroadenM: xi0 is to small -> BroadenM skipped\n");
    return 1;
	}
	if(g<DBL_EPSILON) {
    printf("Warning in BroadenM: g is not positive -> BroadenM skipped\n");
    return 1;
	}
	int i,Ein,Eout;
	double const dE=(double)G->dE;
	double const Gammag=Gamma(g);
	double dataoutIm[G->N+1], sa[G->N+1];
	sa[0]=((pow(dE/(2.*xi0),g))/g-(pow(dE/(2.*xi0),(g+1.)))/(g+1.))/Gammag;
	for(i=1;i<=(int)G->N;i++) {
    sa[i]=dE/((double)i*dE)*pow((double)i*dE/xi0,g)*exp(-(double)i*dE/xi0)/Gammag;
	}
	for(i=0;i<=(int)G->N;i++) {
    dataoutIm[i]=G->Im[i]*sa[0];
	}
#pragma omp parallel for private(Ein)
	for(Eout=0;Eout<=(int)G->N;Eout++){
    for(Ein=0;Ein<Eout;Ein++) {
      dataoutIm[Eout]+=G->Im[Ein]*sa[Eout-Ein];
	  }
  }
	for(i=0;i<=(int)G->N;i++) {
    G->Im[i]=dataoutIm[i];
	}
  return 0;
}

uint32_t BroadenMenergydependent(SpectraType *const G, double const energy[], double const xi0[], double const g[], uint32_t const arraylength) {
	int i,Ein,Eout,j=0;
	double const dE=(double)G->dE;
	double dataoutIm[G->N+1];
	double sa[G->N+1];
	double Enow, xi0now, gnow;
  for(i=1;i<arraylength;i++)if(energy[i-1]>energy[i]){
	  printf("Error in BroadenMenergydependent: energy has to be sorted increasing\n");
    return 1;
	}
  for(i=1;i<arraylength;i++)if(xi0[i]<DBL_EPSILON) {
    printf("Warning in BroadenMenergydependent: xi0 is not positive -> BroadenM skipped\n");
    return 1;
	}
  for(i=1;i<arraylength;i++)if(g[i]<DBL_EPSILON) {
	  printf("Warning in BroadenMenergydependent: g is not positive -> BroadenM skipped\n");
    return 1;
	}
	for(i=0;i<=(int)G->N;i++) dataoutIm[i]=0.0;
	for(Ein=0;Ein<=(int)G->N;Ein++){
    Enow=G->Emin+Ein*dE;
    while(energy[j]<Enow && j<arraylength)j++;
	  if(j==0 || energy[j]<=Enow){
	    xi0now=xi0[j];
      gnow=g[j];
	  }else{
	    xi0now=((Enow-energy[j-1])*xi0[j]+(energy[j]-Enow)*xi0[j-1])/(energy[j]-energy[j-1]);
	    gnow=((Enow-energy[j-1])*g[j]+(energy[j]-Enow)*g[j-1])/(energy[j]-energy[j-1]);
	  }
	  if(xi0now>3.0*dE){																		// maybe not optimal
      sa[0]=((pow(dE/(2.*xi0now),gnow))/gnow-(pow(dE/(2.*xi0now),(gnow+1.)))/(gnow+1.))/Gamma(gnow);
	    for(i=1;i<=(int)G->N-Ein;i++) sa[i]=dE/((double)i*dE)*pow((double)i*dE/xi0now,gnow)*exp(-(double)i*dE/xi0now)/Gamma(gnow);
	    for(Eout=Ein;Eout<=(int)G->N;Eout++) dataoutIm[Eout]+=G->Im[Ein]*sa[Eout-Ein];
	  }else{
	    dataoutIm[Ein]+=G->Im[Ein];
	  }
	}
	for(i=0;i<=(int)G->N;i++) G->Im[i]=dataoutIm[i];
  return 0;
}





void BroadenL(SpectraType *const G, double const FWHM) { // copied from Maurits Haverkort edit by MS
  int i,Ein,Eout;
  double dataoutRe[G->N+1], dataoutIm[G->N+1];
  double const sa=Sqr(FWHM/2.0);
  double const sb=FWHM/(2.0*Pi)*G->dE;
  for(i=0;i<=(int)G->N;i++){
    dataoutRe[i]=0;
    dataoutIm[i]=0;
  }
#pragma omp parallel for private(Ein)
  for(Eout=0;Eout<=(int)G->N;Eout++)
    for(Ein=0;Ein<=(int)G->N;Ein++){
      dataoutRe[Eout]+=(G->Re[Ein]*sb)/(Sqr((Eout-Ein)*G->dE)+sa);
      dataoutIm[Eout]+=(G->Im[Ein]*sb)/(Sqr((Eout-Ein)*G->dE)+sa);
    }
  for(i=0;i<=(int)G->N;i++){
    G->Re[i]=dataoutRe[i];
    G->Im[i]=dataoutIm[i];
  }
}

uint32_t BroadenLshift(SpectraType *const G, double const FWHM, double const deltaE) {
	if(FWHM<3.0*G->dE) {
    printf("Warning in BroadenLshift: FWHM is to small -> BroadenLshift skipped\n");
    return 1;
	}
	int i,Ein,Eout;
	double dataoutRe[G->N+1];
	double dataoutIm[G->N+1];
	double const sa=Sqr(FWHM/2.0);
	double const sb=FWHM/(2.0*Pi)*G->dE;
	for(i=0;i<=(int)G->N;i++){
    dataoutRe[i]=0;
    dataoutIm[i]=0;
	}
#pragma omp parallel for private(Ein)
	for(Eout=0;Eout<=(int)G->N;Eout++)
    for(Ein=0;Ein<=(int)G->N;Ein++){
      dataoutRe[Eout]+=(G->Re[Ein]*sb)/(Sqr((Eout-Ein)*G->dE-deltaE)+sa);
      dataoutIm[Eout]+=(G->Im[Ein]*sb)/(Sqr((Eout-Ein)*G->dE-deltaE)+sa);
    }
	for(i=0;i<=(int)G->N;i++){
    G->Re[i]=dataoutRe[i];
    G->Im[i]=dataoutIm[i];
	}
  return 0;
}

uint32_t BroadenLenergydependent(SpectraType *const G, double const energy[], double const FWHM[], uint32_t const arraylength, double const deltaE) {
	int i,Ein,Eout,j=0;
  for(i=1;i<arraylength;i++)
    if(energy[i-1]>energy[i]){
      printf("Error in BroadenLenergydependent: energy has to be sorted increasing\n");
      return 1;
    }
	double dataoutRe[G->N+1];
	double dataoutIm[G->N+1];
	double sa, sb, Enow, FWHMnow;
	for(i=0;i<=(int)G->N;i++){
    dataoutRe[i]=0;
    dataoutIm[i]=0;
	}
	for(Ein=0;Ein<=(int)G->N;Ein++){
    Enow=G->Emin+Ein*G->dE;
    while(energy[j]<Enow && j<arraylength)j++;
	  if(j==0 || energy[j]<=Enow){
	    FWHMnow=FWHM[j];
	  }else{
	    FWHMnow=((Enow-energy[j-1])*FWHM[j]+(energy[j]-Enow)*FWHM[j-1])/(energy[j]-energy[j-1]);
	  }
	  if(FWHMnow>3.0*G->dE){
      sa=Sqr(FWHMnow/2.0);
	    sb=FWHMnow/(2.0*Pi)*G->dE;
      for(Eout=0;Eout<=(int)G->N;Eout++){
        dataoutRe[Eout]+=(G->Re[Ein]*sb)/(Sqr((Eout-Ein)*G->dE-deltaE)+sa);
        dataoutIm[Eout]+=(G->Im[Ein]*sb)/(Sqr((Eout-Ein)*G->dE-deltaE)+sa);
      }
	  }else if(Ein+(int)(0.5+deltaE/(double)G->dE)<=G->N){
	    dataoutRe[Ein+(int)(0.5+deltaE/(double)G->dE)]+=G->Re[Ein];
      dataoutIm[Ein+(int)(0.5+deltaE/(double)G->dE)]+=G->Im[Ein];
	  }
  }
	for(i=0;i<=(int)G->N;i++){
	  G->Re[i]=dataoutRe[i];
	  G->Im[i]=dataoutIm[i];
	}
  return 0;
}


uint32_t BroadenGenergydependent(SpectraType *const G, double const energy[], double const FWHM[], uint32_t const arraylength, double const deltaE) {
	int i,Ein,Eout,j=0;
  for(i=1;i<arraylength;i++)
    if(energy[i-1]>energy[i]){
      printf("Error in BroadenGenergydependent: energy has to be sorted increasing\n");
      return 1;
    }
	double dataoutRe[G->N+1];
	double dataoutIm[G->N+1];
	double sa, sb, Enow, FWHMnow;
	for(i=0;i<=(int)G->N;i++){
    dataoutRe[i]=0;
    dataoutIm[i]=0;
	}
	for(Ein=0;Ein<=(int)G->N;Ein++){
    Enow=G->Emin+Ein*G->dE;
    while(energy[j]<Enow && j<arraylength)j++;
	  if(j==0 || energy[j]<=Enow){
	    FWHMnow=FWHM[j];
	  }else{
	    FWHMnow=((Enow-energy[j-1])*FWHM[j]+(energy[j]-Enow)*FWHM[j-1])/(energy[j]-energy[j-1]);
	  }
	  if(FWHMnow>3.0*G->dE){
      sa=2.0*Sqr(FWHMnow/(2.0*sqrt(2.0*log(2.0))));
	    sb=(FWHMnow/(2.0*sqrt(2.0*log(2.0))))*sqrt(2.0*Pi)/G->dE;
      for(Eout=0;Eout<=(int)G->N;Eout++){
        dataoutRe[Eout]+=(G->Re[Ein]/sb)*exp(-Sqr((Eout-Ein)*G->dE-deltaE)/sa);
        dataoutIm[Eout]+=(G->Im[Ein]/sb)*exp(-Sqr((Eout-Ein)*G->dE-deltaE)/sa);
      }
	  }else if(Ein+(int)(0.5+deltaE/(double)G->dE)<=G->N){
	    dataoutRe[Ein+(int)(0.5+deltaE/(double)G->dE)]+=G->Re[Ein];
      dataoutIm[Ein+(int)(0.5+deltaE/(double)G->dE)]+=G->Im[Ein];
	  }
	}
	for(i=0;i<=(int)G->N;i++){
    G->Re[i]=dataoutRe[i];
    G->Im[i]=dataoutIm[i];
	}
  return 0;
}


uint32_t SpectraAddPlasmons(SpectraType *Spectra, double const Scale, double const EnergyShift, double const LorentzianFWHM, uint8_t maxOrder, double const IntrinsicExtrinsic){
  if(maxOrder==0){
	  maxOrder=(int)((Spectra->Emax-Spectra->Emin)/EnergyShift)+1;
    printf("Comment to SpectraAddPlasmons: maxOrder was set to %i\n", maxOrder);
  }
  if(Scale>1 || Scale<0) printf("Warning in SpectraAddPlasmons: Scale is not inbetween [0,1]\n"); // 0=Intrinsic ... 1=Extrinsic
  if(EnergyShift<0) printf("Warning in SpectraAddPlasmons: EnergyShift is not positive\n");
  if(LorentzianFWHM<0){
	  printf("Error in SpectraAddPlasmons: Lorentzian has to be positive\n");
    return 1;
	}
  if(IntrinsicExtrinsic>1.+DBL_EPSILON || IntrinsicExtrinsic< -DBL_EPSILON){
	  printf("Error in SpectraAddPlasmons: IntrinsicExtrinsic has to be inbetween [0,1]\n");
    return 1;
	}
  uint32_t i;
  double multi=1.;
  SpectraType nextOrder;
  SpectrumAllocate(*Spectra, &nextOrder);
  CopySpectra(*Spectra, &nextOrder);
  for(i=1;i<maxOrder+1;i++){
    BroadenLshift(&nextOrder, LorentzianFWHM, EnergyShift);
	  if (IntrinsicExtrinsic>1.-DBL_EPSILON){
	    multi*=Scale;
	  }else if(IntrinsicExtrinsic<DBL_EPSILON){
	    multi*=Scale/(double)i;
	  }else{
	    multi*=(1.0-IntrinsicExtrinsic)*Scale/(double)i+IntrinsicExtrinsic*Scale;
	  }
	  AddMultipliedSpectra(Spectra, nextOrder, multi, 0.0);
	}
  return 0;
}


uint32_t CreateRealWaveFunctionShellOccupacy(WaveFunctionType *const psi, uint32_t const Nshell[], uint32_t const nelectrons[], uint32_t const arraylength){
  uint32_t i,NFermions=0;
  for(i=0;i<arraylength;i++){
    NFermions+=Nshell[i]; printf("Nshell=%i ; ",Nshell[i]);
    if(nelectrons[i]>Nshell[i]){
      printf("Error in CreateRealWaveFunctionShellOccupacy: You tried to make a funciton with %i electrons in %i places\n",nelectrons[i],Nshell[i]);
      return 1;
    };
  }
  uint32_t idet,ielec,iferm,iofset,izero=0;
  uint8_t *Determinant;
  psi->NBitsKey=16;
  psi->QComplex=0;
  psi->NBoson=0;
  psi->NFermion=NFermions;
  if(WaveFunctionInit(psi)){
	  printf("WaveFunctionInit failed in CreateRealRandomWaveFunction\n");
	  fflush(stdout);
    return 1;
  };
  sprintf(psi->Name,"WaveFunction");
  for(idet=0;idet<1;idet++){
    Determinant= (uint8_t*)calloc(psi->NWord,sizeof(uint8_t));
	  for(i=0;i<arraylength;i++){
      iofset=0;
      for(ielec=0;ielec<nelectrons[i];ielec++){
        if(rand()&1){
          if(Nshell[i]-nelectrons[i]>iofset)iofset++;
          if(rand()&1){
            if(Nshell[i]-nelectrons[i]>iofset)iofset++;
            if(rand()&1){
              if(Nshell[i]-nelectrons[i]>iofset)iofset++;
              if(rand()&1){
                if(Nshell[i]-nelectrons[i]>iofset)iofset++;
                if(rand()&1){
                  if(Nshell[i]-nelectrons[i]>iofset)iofset++;
                  if(rand()&1){
                    if(Nshell[i]-nelectrons[i]>iofset)iofset++;
                  }
                }
              }
            }
          }
        }
        iferm=izero+ielec+iofset;
        OperatorCr(Determinant,iferm,psi->NFermion);
        printf("iferm=%i ; ",iferm);
      }
      izero+=Nshell[i];
	  }
    RealWaveFunctionAddElement(psi, Determinant, 1.0);
    free(Determinant);
  }
  RealWaveFunctionNormalize(psi);
  return 0;
}


uint32_t CreateRealWaveFunctionNew(WaveFunctionType *const psi, uint32_t const NFermions, uint32_t const nelectrons){
  if(nelectrons>NFermions){
    printf("Error in CreateRealWaveFunctionNew: You tried to make a funciton with %i electrons in %i places\n\n",nelectrons,NFermions);
    return 1;
  };
  uint32_t idet,ielec,jdet;
  uint8_t *Determinant;
  psi->NBitsKey=16;
  psi->QComplex=0;
  psi->NBoson=0;
  psi->NFermion=NFermions;
  if(WaveFunctionInit(psi)){
	  printf("WaveFunctionInit failed in CreateRealRandomWaveFunction\n");fflush(stdout);
    return 1;
  };
  sprintf(psi->Name,"WaveFunction");
  for(idet=0;idet<NFermions-11;idet++)
    for(jdet=idet+1;jdet<NFermions-10;jdet++){
      Determinant= (uint8_t*)calloc(psi->NWord,sizeof(uint8_t));
      for(ielec=0;ielec<nelectrons;ielec++){
        (ielec<10)?OperatorCr(Determinant,ielec,psi->NFermion):OperatorCr(Determinant,10+idet*(1-(ielec%2))+jdet*(ielec%2),psi->NFermion);
      }
      RealWaveFunctionAddElement(psi, Determinant, 1.0);
      free(Determinant);
    }
	RealWaveFunctionNormalize(psi);
  return 0;
}

void PEStest()
{
  int ml, i, j;
  int const ninitials=1, nbands=1, NStates=24+14*nbands;
  double Bext, Rydberg, red3D4F, red4F4F, redZta3D, redZta4F, Zta3D, Zta4F, Delta4FLFbottom, Bandwidth, Delta4FLF[nbands], U3D4F, U4F4F, FK3D4F[3], GK3D4F[3], FK4F4F[4], epsilon, ExpEnergyRe, ExpValRe, ExpValIm, ExpValReSum, ExpValImSum;
  uint16_t Index3Dup[5], Index3Ddn[5], Index4Fup[7], Index4Fdn[7], IndexLFup[nbands][7], IndexLFdn[nbands][7];
  OperatorType OppSz3D, OppSz4F, OppLz3D, OppLz4F, Oppls3D, Oppls4F, OppCF4Fforw, OppCF4Fback, OppCF4F, OppU3D4F, OppU4F4F, OppN3D, Oppct, Oppct2, OppN4F, OppN3D4F, OppN4F4F, OppNLF[nbands], OppH, OppXAS, OppPES[10];
  HarmonicsExpansion Voh, VXAS;
  WaveFunctionType Opsi, psi, psiEigen[ninitials];
  TriDiagonalMatrixType SpecTri;
  SpectraType Spectra, SpectraSum, SpectraPlasmon;
  
  // #pragma omp parallel for private(dynamic)
  
  /*********************/
  /* Define Parameters */
  /*********************/
  
  Bext=0.00000;                                           //later set to "0."
  Rydberg=13.6057;
  
  red3D4F=0.797;
  red4F4F=0.75;
  redZta3D=0.987;
  redZta4F=1.0023;
  
  Zta3D=7.44*redZta3D;                            // Zta3D approx. for 3d9 4f0,1
  Zta4F=0.09*redZta4F;                          // Zta4F valid for 3d10 4f1
  
  U3D4F=11.1;
  
  Bandwidth=0.5;
  Delta4FLFbottom=-2.3;
  if(nbands>1){
    for(i=0;i<nbands;i++) Delta4FLF[i]=Delta4FLFbottom-10*U3D4F-((double)i/(double)(nbands-1))*Bandwidth+Bandwidth/2.;     // compensate U3D4F in initial state
  }else{
    for(i=0;i<nbands;i++) Delta4FLF[i]=Delta4FLFbottom-10*U3D4F;
  }
  GK3D4F[0]=5.1*red3D4F;                     // GK3D4F,FK3D4F valid for 3d10 4f1
  GK3D4F[1]=3.0*red3D4F;
  GK3D4F[2]=2.0*red3D4F;
  FK3D4F[0]=(99.*GK3D4F[0]+44.*GK3D4F[1]+50.*GK3D4F[2])/2310.;
  FK3D4F[1]=7.5*red3D4F;
  FK3D4F[2]=3.4*red3D4F;
  
  U4F4F=8.6;
  FK4F4F[1]=11.4*red4F4F;                                            // FK4F4F valid for 3d10 4f2
  FK4F4F[2]=7.1*red4F4F;
  FK4F4F[3]=5.1*red4F4F;
  FK4F4F[0]=(572.*FK4F4F[1]+390.*FK4F4F[2]+500.*FK4F4F[3])/27885.;
  
  /********************************************************/
  /* Set shellpositions: 3D={0,9}, 4F={10,23}, LF={24,37} */
  /********************************************************/
  for(ml=0;ml<5;ml++){
    Index3Dup[ml]=2*ml;
    Index3Ddn[ml]=2*ml+1;
  }
  for(ml=0;ml<7;ml++){
    Index4Fup[ml]=2*ml+10;
    Index4Fdn[ml]=2*ml+11;
  }
  for(i=0;i<nbands;i++)
    for(ml=0;ml<7;ml++){
      IndexLFup[i][ml]=2*ml+24+14*i;
      IndexLFdn[i][ml]=2*ml+25+14*i;
    }
  
  /**************/
  /* Create Voh */
  /**************/
  int Vohl[]={0};
  int Vohm[]={0};
  double VohRe[]={0.6/sqrt((double)nbands)};
  Voh.QComplex=0;
  Voh.N=1;
  Voh.l=Vohl;
  Voh.m=Vohm;
  Voh.Re=VohRe;
  
  /***************/
  /* Create VXAS */
  /***************/
  int VXASl[]={1,1,1};
  int VXASm[]={-1,0,1};
  double VXASRe[]={1.,1.,1.};
  VXAS.QComplex=0;
  VXAS.N=3;
  VXAS.l=VXASl;
  VXAS.m=VXASm;
  VXAS.Re=VXASRe;
  
  /********************/
  /* Create Operators */
  /********************/
  CreateOperatorSz   (NStates, 5, Index3Dup, Index3Ddn, &OppSz3D);
  CreateOperatorSz   (NStates, 7, Index4Fup, Index4Fdn, &OppSz4F);
  // CreateOperatorSz   (NStates, 7, IndexLFup, IndexLFdn, &OppSzLF);
  CreateOperatorLz   (NStates, 2, Index3Dup, Index3Ddn, &OppLz3D);
  CreateOperatorLz   (NStates, 3, Index4Fup, Index4Fdn, &OppLz4F);
  CreateOperatorldots(NStates, 2, Index3Dup, Index3Ddn, &Oppls3D);
  CreateOperatorldots(NStates, 3, Index4Fup, Index4Fdn, &Oppls4F);
  /* Hopping */
  CreateOperatorCF(NStates, 3, 3, Index4Fup, Index4Fdn, IndexLFup[0], IndexLFdn[0], Voh, &OppCF4Fforw);
  CreateOperatorCF(NStates, 3, 3, IndexLFup[0], IndexLFdn[0], Index4Fup, Index4Fdn, Voh, &OppCF4Fback);
  OperatorCopy(&OppCF4F, OppCF4Fback);
  OperatorAdd (&OppCF4F, OppCF4Fforw, 1 ,0);
  for(i=1;i<nbands;i++){
    CreateOperatorCF(NStates, 3, 3, Index4Fup, Index4Fdn, IndexLFup[i], IndexLFdn[i], Voh, &OppCF4Fforw);
    CreateOperatorCF(NStates, 3, 3, IndexLFup[i], IndexLFdn[i], Index4Fup, Index4Fdn, Voh, &OppCF4Fback);
    OperatorAdd(&OppCF4F, OppCF4Fback, 1 ,0);
    OperatorAdd(&OppCF4F, OppCF4Fforw, 1 ,0);
  }
  /* Onsite Coulomb repulsion */
  CreateOperatorCoulombSingleShell(NStates, 3, Index4Fup, Index4Fdn, FK4F4F, &OppU4F4F);
  InitOperatorCopySettings(&OppN4F4F, OppU4F4F);
  for(i=10;i<23;i++){
    CreateOperatorNumber(NStates, i, i, &Oppct);
    for(j=i+1;j<24;j++){
      CreateOperatorNumber(NStates, j, j, &Oppct2);
      OperatorMultiplyNormalOrder(&OppN4F4F, Oppct, Oppct2, &QDetNotOrderedFcFaBcBa);
      FreeOperator(&Oppct2);
    }
    FreeOperator(&Oppct);
  }
  /* Corehole attraction */
  CreateOperatorCoulombTwoShellsConstantOccupation(NStates, 2, 3, Index3Dup, Index3Ddn, Index4Fup, Index4Fdn, FK3D4F, GK3D4F, &OppU3D4F);
  InitOperatorCopySettings(&OppN3D4F, OppU4F4F);
  for(i=0;i<10;i++){
    CreateOperatorNumber(NStates, i, i, &Oppct);
    for(j=10;j<24;j++){
      CreateOperatorNumber(NStates, j, j, &Oppct2);
      OperatorMultiplyNormalOrder(&OppN3D4F, Oppct, Oppct2, &QDetNotOrderedFcFaBcBa);
      FreeOperator(&Oppct2);
    }
    FreeOperator(&Oppct);
  }
  /* Delta4F */
  CreateOperatorNumber(NStates, 10, 10, &OppN4F);
  for(i=11;i<24;i++){
    CreateOperatorNumber(NStates, i, i, &Oppct);
    OperatorAdd(&OppN4F, Oppct, 1 ,0);
    FreeOperator(&Oppct);
  }
  CreateOperatorCF(NStates, 3, 2, Index4Fup, Index4Fdn, Index3Dup, Index3Ddn, VXAS, &OppXAS);
  /* PES */
  for(i=0;i<10;i++){
    CreateOperatorAnnihilate(NStates, i, &OppPES[i]);
	}
  /* N3D */
  CreateOperatorNumber(NStates, 0, 0, &OppN3D);
  for(i=1;i<10;i++){
    CreateOperatorNumber(NStates, i, i, &Oppct);
    OperatorAdd(&OppN3D, Oppct, 1 ,0);
    FreeOperator(&Oppct);
  }
  /* NLF */
  for(i=0;i<nbands;i++){
    CreateOperatorNumber(NStates, 24+14*i, 24+14*i, &OppNLF[i]);
    for(j=0;j<13;j++){
      CreateOperatorNumber(NStates, 25+14*i+j, 25+14*i+j, &Oppct);
      OperatorAdd(&OppNLF[i], Oppct, 1 ,0);
      FreeOperator(&Oppct);
    }
  }
  OperatorCopy(&OppH, OppCF4F);
  OperatorAdd (&OppH, OppU3D4F, 1.,0);
  OperatorAdd (&OppH, OppU4F4F, 1.,0);
  //OperatorAdd (&OppH, OppSz3D, 2*Bext,0);
  //OperatorAdd (&OppH, OppSz4F, 2*Bext,0);
  //OperatorAdd (&OppH, OppSzLF, 2*Bext,0);
  //OperatorAdd (&OppH, OppLz3D, Bext,0);
  //OperatorAdd (&OppH, OppLz4F, Bext,0);
  OperatorAdd (&OppH, Oppls3D, Zta3D,0);
  OperatorAdd (&OppH, Oppls4F, Zta4F,0);
  OperatorAdd (&OppH, OppN4F, Delta4FLF[0], 0);
  for(i=0;i<nbands;i++) OperatorAdd (&OppH, OppNLF[i], Delta4FLF[0]-Delta4FLF[i], 0);
  OperatorAdd (&OppH, OppN3D4F, U3D4F,0);
  OperatorAdd (&OppH, OppN4F4F, U4F4F,0);
  OperatorAdd (&OppH, OppN3D, 0.0,0);
  
  /********/
  /* test */
  /********/
  if(1){
    printf("\n\n\n\n--------------------Sz---------------------\n");
    PrintOperatorReadable(OppSz4F);
    printf("\n\n\n\n--------------------Lz---------------------\n");
    PrintOperatorReadable(OppLz4F);
    printf("\n\n\n\n--------------------ls---------------------\n");
    PrintOperatorReadable(Oppls4F);
    printf("\n\n\n\n--------------------Ufc--------------------\n");
    PrintOperatorReadable(OppU3D4F);
    printf("\n\n\n\n--------------------Uff--------------------\n");
    PrintOperatorReadable(OppU4F4F);
    printf("\n\n\n\n--------------------CF---------------------\n");
    PrintOperatorReadable(OppCF4F);
    printf("\n\n\n\n--------------------N4F--------------------\n");
    PrintOperatorReadable(OppN4F);
    for(i=0;i>nbands;i++){
      printf("\n\n\n\n------------------NLF[%i]------------------\n",i);
      PrintOperatorReadable(OppNLF[i]);
    }
    printf("\n\n\n\n-------------------N4F4F-------------------\n");
    PrintOperatorReadable(OppN4F4F);
    printf("\n\n\n\n-------------------N3D4F-------------------\n");
    PrintOperatorReadable(OppN3D4F);
    printf("\n\n\n\n---------------------H---------------------\n");
    PrintOperatorReadable(OppH);
    printf("\n\n\n\n--------------------XAS--------------------\n");
    PrintOperatorReadable(OppXAS);
    printf("\n\n\n\n--------------------PES--------------------\n");
    for(i=0;i<10;i++)PrintOperatorReadable(OppPES[i]);
  }
  if(0){ /* XAS part not yet correct mwh - create real random wave function can't be. (holes in d shell...?...)*/
    epsilon=0.01;
    CreateRealRandomWaveFunction(&psi, NStates, 24);
    LanczosEigenSystem(OppH, psi, psiEigen, &epsilon, 1000*DBL_EPSILON, 5);
    OperatorPsi(OppXAS, psiEigen[0], &Opsi);
    SpecTri.N=100;
    InitTriDiagonalMatrix(&SpecTri);
    LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri, 0.6);
	  ExpEnergyRe=0.;
	  Spectra.N=5000;
	  Spectra.Emin=308.;
	  Spectra.Emax=408.;
	  Spectra.dE=((Spectra.Emax-Spectra.Emin)/(double)Spectra.N);
	  Spectra.Gamma=0.16;
	  SpectrumAllocate(Spectra, &Spectra);
    TriDiagonalMatrixToSpectrum(ExpEnergyRe, SpecTri, &Spectra);
    //   sprintf("PEStest.dat","Spectra_l%i_N%i",3,NStates);
    FilePrintSpectra(&Spectra, 1, "PEStest2.dat", "w");
    FreeTriDiagonalMatrix(&SpecTri);
  }
  if(1){
    epsilon=0.0001;																// find optimum here: ~0.0001
    // CreateRealRandomWaveFunction(&psi, NStates, 12);
    uint32_t crshells[]={10,NStates-10};
    uint32_t crelectrons[]={10,2};
    CreateRealWaveFunctionShellOccupacy(&psi, crshells, crelectrons, 2);
    // CreateRealWaveFunctionNew(&psi, NStates, 12);
    PrintWaveFunctionReadable(psi);
    
    if(1){
      PsiOperatorPsi(OppH, psi, &ExpValRe, &ExpValIm);
      printf("E=%f + %f i  |  ",ExpValRe,ExpValIm);
      PsiOperatorPsi(OppN3D, psi, &ExpValRe, &ExpValIm);
      printf("N3D=%f + %f i  |  ",ExpValRe,ExpValIm);
      PsiOperatorPsi(OppN4F, psi, &ExpValRe, &ExpValIm);
      printf("N4F=%f + %f i  |  ",ExpValRe,ExpValIm);
      PsiOperatorPsi(OppN4F4F, psi, &ExpValRe, &ExpValIm);
      printf("N4F2=%f + %f i  |  ",ExpValRe,ExpValIm);
      ExpValReSum=0.; ExpValImSum=0.;
      for(i=0;i<nbands;i++){PsiOperatorPsi(OppNLF[i], psi, &ExpValRe, &ExpValIm);
        ExpValReSum+=ExpValRe; ExpValImSum+=ExpValIm;
      }
      printf("NLF=%f + %f i \n",ExpValReSum,ExpValImSum);
    }
    
    LanczosEigenSystem(OppH, psi, psiEigen, &epsilon, 1000*DBL_EPSILON, ninitials);		// find optimum here:
    
    if(1){
      for(i=0;i<ninitials;i++){
        PsiOperatorPsi(OppH, psiEigen[i], &ExpValRe, &ExpValIm);
        printf("%i:   E=%f + %f i  |  ",i,ExpValRe,ExpValIm);
        PsiOperatorPsi(OppN3D, psiEigen[i], &ExpValRe, &ExpValIm);
        printf("N3D=%f + %f i  |  ",ExpValRe,ExpValIm);
        PsiOperatorPsi(OppN4F, psiEigen[i], &ExpValRe, &ExpValIm);
        printf("N4F=%f + %f i  |  ",ExpValRe,ExpValIm);
        PsiOperatorPsi(OppN4F4F, psiEigen[i], &ExpValRe, &ExpValIm);
        printf("N4F2=%f + %f i  |  ",ExpValRe,ExpValIm);
        ExpValReSum=0.; ExpValImSum=0.;
        for(j=0;j<nbands;j++){PsiOperatorPsi(OppNLF[j], psiEigen[i], &ExpValRe, &ExpValIm);
          ExpValReSum+=ExpValRe; ExpValImSum+=ExpValIm;
        }
        printf("NLF=%f + %f i \n",ExpValReSum,ExpValImSum);
      }
    }
    SpectraSum.N=2000;
    SpectraSum.Emin=-40.;
    SpectraSum.Emax=40.;
    SpectraSum.dE=((SpectraSum.Emax-SpectraSum.Emin)/(double)SpectraSum.N);
    SpectraSum.Gamma=0.4;
    SpectrumAllocate(SpectraSum, &SpectraSum);
    for(i=0;i<10;i++){                                                       // for fast testing set 0
      OperatorPsi(OppPES[i], psiEigen[0], &Opsi);
      SpecTri.N=100;															// find optimum here: ~200
      InitTriDiagonalMatrix(&SpecTri);
      LanczosTriDiagonalize(OppH, Opsi, epsilon, &SpecTri, 0.1);				// find optimum here: ~0.1
      ExpEnergyRe=-2.858;                                               // due to 3P10 in initial state and 2.858 shifted from xtls
      SpectrumAllocate(SpectraSum, &Spectra);
      SpectrumAllocate(SpectraSum, &SpectraPlasmon);
      TriDiagonalMatrixToSpectrum(ExpEnergyRe, SpecTri, &Spectra);
      //   sprintf("PEStest.dat","Spectra_l%i_N%i",3,NStates);
      //   BroadenG(&Spectra, 1.0);												// ok
      //   double Genergy[]={-10.0, 0.0}, GFWHM[]={0.0, 1.0};
      //   BroadenGenergydependent(&Spectra, Genergy, GFWHM, 2, 10.0);				// ok
      //	BroadenM(&Spectra, 10.0, 0.1);											// ok
      //   double Menergy[]={-10.0, 0.0}, Mxi0[]={0.0, 10.0}, Mg[]={0.1, 0.1};
      //   BroadenMenergydependent(&Spectra, Menergy, Mxi0, Mg, 2);						// ok
      //   BroadenL(&Spectra, 0.3);   												// ok
      //   double Lenergy[]={-10.00001,-9.99999}, LFWHM[]={0.0, 0.3};
      //   BroadenLenergydependent(&Spectra, Lenergy, LFWHM, 2, 0.0);   				// ok
      (i==0)?CopySpectra(Spectra, &SpectraSum):AddMultipliedSpectra(&SpectraSum, Spectra, 1.0, 0.0);
      FreeTriDiagonalMatrix(&SpecTri);
    }
    CopySpectra(SpectraSum, &SpectraPlasmon);
    SpectraAddPlasmons(&SpectraPlasmon, 0.3, 7.0, 2.0, 2, 0.9);
    FilePrintSpectra(&SpectraSum, 1, "PEStest2.dat","w");
    FilePrintSpectra(&SpectraPlasmon, 1, "PEStest2Pl.dat","w");
  }
  printf("Finished!!!\n");
}


int main (int argc, const char * argv[])
{
  TimeStampAll();
  argc-=1;
  for (;argc!=0;argc--)
  {
    if(!strncmp("test",argv[argc],strlen(argv[argc]))) test();
    if(!strncmp("pes",argv[argc],strlen(argv[argc]))) PEStest();
  }
  return 0;
}
