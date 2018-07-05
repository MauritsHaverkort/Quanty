#include "TestSpectral.h"

void TestGreenFunctionCreation()
{
  uint32_t i, Dim, Nk, NTriDiag;
  char   filename[255];
  TriDiagonalMatrixType G;
  SpectraType Spec[4];
  Nk = 16;
  NTriDiag = 1000;
  Spec[0].N = 22000;
  Spec[0].Emin = -1.1;
  Spec[0].Emax = 1.1;
  Spec[0].EminPole = -1.1;
  Spec[0].EmaxPole = 1.1;
  Spec[0].dE = (Spec[0].Emax - Spec[0].Emin) / Spec[0].N;
  Spec[0].Gamma = 0.0005;
  Spec[0].phase = 0.0;
  for (i = 0; i < 4; i++)
  {
    SpectrumAllocate(Spec[0], &(Spec[i]));
  }
  NearestNeighborTriangularTBToGSpectraRepresentation(&Spec[0], Nk);
  /*TriDiagonalMatrixToSpectrum(0.0, G, &Spec[0]); */

  for (Dim = 1; Dim <= 4; Dim++)
  {
    printf("%i\n", Dim);
    NearestNeighborHyperCubicTBToGTriDiagonalRepresentationNonOrthogonal(Dim, NTriDiag, &G, Nk);
    NearestNeighborHyperCubicTBToGTriDiagonalRepresentation(Dim, NTriDiag, &G, Nk);
    TriDiagonalMatrixToSpectrum(0.0, G, &Spec[Dim - 1]);
    TriDiagonalMatrixToSpectrumGaussian(0.0, G, &Spec[Dim - 0], Spec[Dim - 0].Gamma);
    sprintf(filename, "G0CubicD%i", Dim);
    PrintPlainTriDiagonalMatrix(G);
    FilePrintPlainTriDiagonalMatrix(G, filename, "w");
    FreeTriDiagonalMatrix(&G);
    Nk /= 4;
  }

  FilePrintSpectra(Spec, 4, "G0_test", "w");
}


void TestSpec()
{
/*  double *specOrig,*specGaus,*specTmp,*spec01,*spec02V,*spec02C,*spec02,*spec03,*spec04,*spec05,*spec06;
  double Emin,Emax,E0,Gamma,DE,mu,zero,hc;
  uint32_t N,i;
  TriDiagonalMatrixType M, MV, MC, MVs, MCs;
  
  mu=0.0;
  zero=1E-10;
  hc=3.0;
  
  N=2000;
  Emin=-20.0;
  Emax=20.0;
  E0=-5.0;
  Gamma=0.2;
  DE=(Emax-Emin)/((double)N);

  specOrig= (double*)calloc(2*(N+1),sizeof(double));
  specGaus= (double*)calloc(2*(N+1),sizeof(double));
  specTmp = (double*)calloc(2*(N+1),sizeof(double));
  spec01  = (double*)calloc(2*(N+1),sizeof(double));
  spec02  = (double*)calloc(2*(N+1),sizeof(double));
  spec02V = (double*)calloc(2*(N+1),sizeof(double));
  spec02C = (double*)calloc(2*(N+1),sizeof(double));
  spec03  = (double*)calloc(2*(N+1),sizeof(double));
  spec04  = (double*)calloc(2*(N+1),sizeof(double));
  spec05  = (double*)calloc(2*(N+1),sizeof(double));
  spec06  = (double*)calloc(2*(N+1),sizeof(double));

  M.a= (double*)calloc(N+2,sizeof(double));
  M.b= (double*)calloc(N+2,sizeof(double));
  MV.a= (double*)calloc(200,sizeof(double));
  MV.b= (double*)calloc(200,sizeof(double));
  MC.a= (double*)calloc(200,sizeof(double));
  MC.b= (double*)calloc(200,sizeof(double));
  MVs.a= (double*)calloc(N+2,sizeof(double));
  MVs.b= (double*)calloc(N+2,sizeof(double));
  MCs.a= (double*)calloc(N+2,sizeof(double));
  MCs.b= (double*)calloc(N+2,sizeof(double));

  Lorentzian(Emin,Emax,N,Gamma,E0+hc,specOrig);
  for(i=1;i<=100;i++)
  {
    Lorentzian(Emin,Emax,N,Gamma,hc+E0-(i*E0/50.0),specTmp);
    VectorAdd(specOrig,specTmp,2*N+2);
  }
  VectorMultiply(specOrig,2*N+2,1.0/101.0);
  ReplaceLorentzianByGaussian(Emin,Emax,N,Gamma,&specOrig[N+1],specGaus);
  
  ProjectGOntoAndersonImpurityModel(Emin, Emax, N, Gamma, specOrig, hc,&M);
  AndersonImpurityMatrixToSpectrum(0.0, Emin, Emax, N, Gamma, M, spec01);
  SplitAndersonImpurityHamiltonian(M, mu, zero, &MVs,&MCs);
  MV.N=200;
  MC.N=200;
  FromAndersonImpurityToTriDiagonalMatrix(MVs, &MV);
  FromAndersonImpurityToTriDiagonalMatrix(MCs, &MC);
  TriDiagonalMatrixToSpectrum(0.0, Emin, Emax, N, Gamma, MV, spec02V);
  TriDiagonalMatrixToSpectrum(0.0, Emin, Emax, N, Gamma, MC, spec02C);
  DoubleTriDiagonalMatrixToSpectrum(0.0, Emin, Emax, N, Gamma, MV, MC, spec02);
  MV.N=100;
  MC.N=100;
  DoubleTriDiagonalMatrixToSpectrum(0.0, Emin, Emax, N, Gamma, MV, MC, spec03);
  MV.N=50;
  MC.N=50;
  DoubleTriDiagonalMatrixToSpectrum(0.0, Emin, Emax, N, Gamma, MV, MC, spec04);
  MV.N=25;
  MC.N=25;
  DoubleTriDiagonalMatrixToSpectrum(0.0, Emin, Emax, N, Gamma, MV, MC, spec05);
  MV.N=10;
  MC.N=10;
  DoubleTriDiagonalMatrixToSpectrum(0.0, Emin, Emax, N, Gamma, MV, MC, spec06);
  
  for(i=0;i<=N;i++)
  {
    printf("%22.15E %22.15E %22.15E %22.15E %22.15E %22.15E %22.15E ",(double)(i*DE+Emin),specOrig[i],specOrig[i+N+1],specGaus[i],M.b[i+1],spec01[i],spec01[i+N+1]);
    printf("%22.15E %22.15E %22.15E %22.15E %22.15E %22.15E ",spec02V[i],spec02V[i+N+1],spec02C[i],spec02C[i+N+1],spec02[i],spec02[i+N+1]);
    printf("%22.15E %22.15E %22.15E %22.15E %22.15E %22.15E %22.15E %22.15E\n",spec03[i],spec03[i+N+1],spec04[i],spec04[i+N+1],spec05[i],spec05[i+N+1],spec06[i],spec06[i+N+1]);
  }
  exit(0);
  */
}

void TestKKTrans()
{
/*
  uint32_t N=2000;
  uint32_t i;
  double G[2*N+2];
  double Gnew[2*N+2];
  double Gamma=10.0,a=500.0,b=700.0;
  for(i=0;i<=N;i++)
  {
    G[i    ]=5*((double)i-a)/(Sqr((double)i-a)+Sqr(Gamma))+((double)i-b)/(Sqr((double)i-b)+Sqr(Gamma*2));
    G[i+N+1]=5*(     -Gamma)/(Sqr((double)i-a)+Sqr(Gamma))+(   -2*Gamma)/(Sqr((double)i-b)+Sqr(Gamma*2));
    Gnew[i]=0;
    Gnew[i+N+1]=G[i+N+1];
  }
  KramersKronigImToRe(Gnew, N);
  for(i=0;i<=N;i++)
  {
    printf("%i %lf %lf %lf %lf\n",i,G[i],G[i+N+1],Gnew[i],Gnew[i+N+1]);
  }
  */
}

void TestSpecSelf()
{
/*
  double *G0,*G,*S,*specTmp;
  double Emin,Emax,E0,Gamma,DE;
  uint32_t N,i;
  
  N=2000;
  Emin=-20.0;
  Emax=20.0;
  E0=-5.0;
  Gamma=0.2;
  DE=(Emax-Emin)/((double)N);
  G0=Lorentzian(Emin,Emax,N,Gamma,E0);
  for(i=1;i<=100;i++)
  {
    specTmp=Lorentzian(Emin,Emax,N,Gamma,E0-(i*E0/50.0));
    VectorAdd(G0,specTmp,2*N+2);
    free(specTmp);
  }
  VectorMultiply(G0,2*N+2,1.0/101.0);
  
  E0=-2.5;
  G=Lorentzian(Emin,Emax,N,Gamma,2*E0-1);
  specTmp=Lorentzian(Emin,Emax,N,Gamma,-2*E0+1);
  VectorAdd(G,specTmp,2*N+2);
  for(i=1;i<=100;i++)
  {
    specTmp=Lorentzian(Emin,Emax,N,Gamma,2*E0-1-(i*E0/50.0));
    VectorAdd(G,specTmp,2*N+2);
    free(specTmp);
    specTmp=Lorentzian(Emin,Emax,N,Gamma,-2*E0+1+(i*E0/50.0));
    VectorAdd(G,specTmp,2*N+2);
    free(specTmp);
  }
  VectorMultiply(G,2*N+2,1.0/202.0);
  S=CalculateSelfEnergySpectralRepresentation(G0, G, N);
  for(i=0;i<=N;i++)
  {
    printf("%22.15E %22.15E %22.15E %22.15E %22.15E %22.15E %22.15E\n",(double)(i*DE+Emin),G0[i],G0[i+N+1],G[i],G[i+N+1],S[i],S[i+N+1]);
  }
*/
}
