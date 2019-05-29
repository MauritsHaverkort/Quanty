/*
 * ResonantSpectrum_Selfenergy.cpp
 *
 *  Created on: May 20, 2019
 *      Author: brass
 */

#include "Modules.h"
#include "ResonantSpectrum_Selfenergy.h"

static int LuaCreateResonantSpectrum(lua_State *L)
{
  struct
  {
    ConvergenceSettingsType Settings;
    double E0, E01, E02;
    uint32_t Npsi;
    uint32_t NT1;
    uint32_t NT2;
    uint32_t NTriSE;
    bool Tensor;
    bool DenseMethod;
  } Options;
  Options.Settings=InitConvergenceSettings();
  Options.E0=DBL_MAX;
  Options.E01=DBL_MAX;
  Options.E02=DBL_MAX;
  Options.Tensor=false;
  Options.DenseMethod=false;
  Options.NTriSE = 5;
  ConvergenceSettingsType SettingsNow;
  SpectraTypeArray *Spec;
  SpectraType *SpecT1psi;
  TriDiagonalMatrixType *SpecTri;
  OperatorType *H1,*H2,*T1, *T2;
  WaveFunctionType *psi, HTpsiDenseBasis,HTHTpsiDenseBasis,Opsi,TpsiSum;
  uint32_t iT,ipsi,len,i;
  LuaOptionsType LuaCreateResonantSpectraOptions[]={
    {"NTri1","uint32_t",&(Options.Settings.NTri1)},
    {"NTri2","uint32_t",&(Options.Settings.NTri2)},
    {"Tensor","bool",&(Options.Tensor)},
    {"restrictions1", "DeterminantRestrictionType", &(Options.Settings.restrictions1)},
    {"restrictions2", "DeterminantRestrictionType", &(Options.Settings.restrictions2)},
    {"Restrictions1", "DeterminantRestrictionType", &(Options.Settings.restrictions1)},
    {"Restrictions2", "DeterminantRestrictionType", &(Options.Settings.restrictions2)},
    {"epsilon","double",&(Options.Settings.epsilon)},
    {"Epsilon","double",&(Options.Settings.epsilon)},
    {"E0","double",&(Options.E0)},
    {"E0_1","double",&(Options.E01)},
    {"E0_2","double",&(Options.E02)},
    {"Emin1","double",&(Options.Settings.Emin1)},
    {"Emax1","double",&(Options.Settings.Emax1)},
    {"NE1","uint32_t",&(Options.Settings.NE1)},
    {"Gamma1","double",&(Options.Settings.Gamma1)},
    {"Emin2","double",&(Options.Settings.Emin2)},
    {"Emax2","double",&(Options.Settings.Emax2)},
    {"NE2","uint32_t",&(Options.Settings.NE2)},
    {"Gamma2","double",&(Options.Settings.Gamma2)},
    {"DenseBorder","uint32_t",&(Options.Settings.DenseBorder)},
    {"ForceDenseMethod","uint32_t",&(Options.DenseMethod)},
	{"NTriSE","uint32_t",&(Options.NTriSE)},
    {NULL, NULL, NULL}};
  printf("Start ResonantSpectrum\n");
  /* free garbage full cycle*/
  if(L->l_G->gcrunning) lua_gc(L, LUA_GCCOLLECT, 0);
  /* check number of input parameters */
  len = (uint32_t)lua_gettop(L);
  if(!((len==9)||(len==10)))
  {
    return luaL_error(L, "Error LuaCreateResonantSpectra called with %d arguments,10 or 9 expected (H1, H2, T1, T2^+, psi,U,EnergyContinuum {options})\n",len);
  }
  /* read operator and wavefunctions */
  H1 = (OperatorType *) luaL_checkudata(L, 1, "Operator_Type");
  H2 = (OperatorType *) luaL_checkudata(L, 2, "Operator_Type");
  if(lua_type(L,3)==LUA_TTABLE)
  {
    T1=luaL_checkarray_OperatorType(L, 3, &Options.NT1);
    if(Options.NT1==0) return luaL_error(L, "Errof CreateResonantSpectra called with a table of transition operators (1) of length 0\n");
  }
  else
  {
    Options.NT1 = 1;
    T1 = (OperatorType *) luaL_checkudata(L, 3, "Operator_Type");
  }
  if(lua_type(L,4)==LUA_TTABLE)
  {
    T2=luaL_checkarray_OperatorType(L, 4, &Options.NT2);
    if(Options.NT2==0) return luaL_error(L, "Errof CreateResonantSpectra called with a table of transition operators (2) of length 0\n");
  }
  else
  {
    Options.NT2 = 1;
    T2 = (OperatorType *) luaL_checkudata(L, 4, "Operator_Type");
  }
  if(lua_type(L,5)==LUA_TTABLE)
  {
    psi=luaL_checkarray_WaveFunctionType(L, 5, &Options.Npsi);
    if(Options.Npsi==0) return luaL_error(L, "Errof CreateResonantSpectra called with a table of wave functions of length 0\n");
  }
  else
  {
    Options.Npsi = 1;
    psi = (WaveFunctionType *) luaL_checkudata(L, 5, "Wavefunction_Type");
  }
  vector<vector<OperatorType> > U;
  vector<vector<vector<double> > > Fk;
  lua_check(L,6,U);
  vector<vector<double> > EnergyContinuum;
  lua_check(L,7,Fk);
  lua_check(L,8,EnergyContinuum);
  double gammaSE;
  lua_check(L,9,gammaSE);
  /* read options */
  if (len==10)
  {
    LuaCheckOptions(L, -1, LuaCreateResonantSpectraOptions);
  }
  if(Options.E01==DBL_MAX)Options.E01=Options.E0;
  if(Options.E02==DBL_MAX)Options.E02=Options.E0;

  /* test if operator is Hermitian */
  switch(QOperatorHermitian(*H1, 0.0001*Options.Settings.zero))
  {
    case 3 :
    {
      printf("WARNING: QOperatorHermitian failed in LuaCreateResonantSpectra, this should not happen\n");
      printf("my guess is that you are out of memory and will crash soon\n");
      break;
    }
    case 2 :
    {
      luaL_warning(L, "CreateResonantSpectra called with an operator (H1) that is not fully Hermitian\nStrongly sugested to stop the calculation and repair this\n\n\n\nDo it\n\n\nReally Check this\n\n\nDid you?\n\n\n");
      break;
    }
    case 1 :
    {
      // all fine H is Hermitian upto zero
      break;
    }
    case 0 :
    {
      luaL_error(L, "CreateResonantSpectra called with an operator (H1) that is not Hermitian\n");
      break;
    }
  }
  switch(QOperatorHermitian(*H2, 0.0001*Options.Settings.zero))
  {
    case 3 :
    {
      printf("WARNING: QOperatorHermitian failed in LuaCreateResonantSpectra, this should not happen\n");
      printf("my guess is that you are out of memory and will crash soon\n");
      break;
    }
    case 2 :
    {
      luaL_warning(L, "CreateResonantSpectra called with an operator (H2) that is not fully Hermitian\nStrongly sugested to stop the calculation and repair this\n\n\n\nDo it\n\n\nReally Check this\n\n\nDid you?\n\n\n");
      break;
    }
    case 1 :
    {
      // all fine H is Hermitian upto zero
      break;
    }
    case 0 :
    {
      luaL_error(L, "CreateResonantSpectra called with an operator (H2) that is not Hermitian\n");
      break;
    }
  }
  /* allocate spectra and a pointer to return matrix on lua stack */
  Spec = (SpectraTypeArray *) lua_newuserdata(L, sizeof(SpectraTypeArray));
  luaL_setmetatable(L, "Spectra_Type");
  Spec->NSpectra=0;
  Spec->Spec=NULL;
  /* test if dense methods can be used */
  // at the same time create the wave-functions that define the basis sets of the dense matrices used in the dense methods
  // we need the vectors:
  // HmDenseBasis       : on which the Matrix for the intermediate state can be represented
  // HfDenseBasis       : on which the Matrix fo rthe final states can be represented
		  {
			// psi
			Options.DenseMethod=false;
			goto EndOfUseDenseTest; // TODO Implement DenseMethod
			if (Options.DenseMethod==false) if(psi[0].N>Options.Settings.DenseBorder) goto EndOfUseDenseTest;
			WaveFunctionCopy(psi[0], &HTpsiDenseBasis);
			for(i=1;i<Options.Npsi;i++)
			{
			  WaveFunctionAddUnionBasis(&HTpsiDenseBasis, psi[i]);
			  if (Options.DenseMethod==false) if(HTpsiDenseBasis.N>Options.Settings.DenseBorder)
			  {
				WaveFunctionFree(&HTpsiDenseBasis);
				goto EndOfUseDenseTest;
			  }
			}
			// T psi
			if((Options.Settings.restrictions1.N==0)?(OperatorPsi(T1[0], HTpsiDenseBasis, &TpsiSum)):(OperatorPsiRestricted(T1[0], HTpsiDenseBasis, &TpsiSum, Options.Settings.restrictions1)))
			{
			  WaveFunctionFree(&HTpsiDenseBasis);
			  goto EndOfUseDenseTest;
			}
			for(i=1;i<Options.NT1;i++)
			{
			  if((Options.Settings.restrictions1.N==0)?(OperatorPsi(T1[i], HTpsiDenseBasis, &Opsi)):(OperatorPsiRestricted(T1[i], HTpsiDenseBasis, &Opsi, Options.Settings.restrictions1)))
			  {
				WaveFunctionFree(&TpsiSum);
				WaveFunctionFree(&HTpsiDenseBasis);
				goto EndOfUseDenseTest;
			  }
			  WaveFunctionAddUnionBasis(&TpsiSum, Opsi);
			  WaveFunctionFree(&Opsi);
			  if (Options.DenseMethod==false) if(TpsiSum.N>Options.Settings.DenseBorder)
			  {
				WaveFunctionFree(&TpsiSum);
				WaveFunctionFree(&HTpsiDenseBasis);
				goto EndOfUseDenseTest;
			  }
			}
			WaveFunctionFree(&HTpsiDenseBasis);
			HTpsiDenseBasis=TpsiSum;
			WaveFunctionNormalize(&HTpsiDenseBasis);
			// H T psi
			len=0;
			while(len!=HTpsiDenseBasis.N)
			{
			  len=HTpsiDenseBasis.N;
			  if((Options.Settings.restrictions1.N==0)?(OperatorPsi(*H1, HTpsiDenseBasis, &Opsi)):(OperatorPsiRestricted(*H1, HTpsiDenseBasis, &Opsi, Options.Settings.restrictions1)))
			  {
				WaveFunctionFree(&HTpsiDenseBasis);
				goto EndOfUseDenseTest;
			  }
			  WaveFunctionFree(&HTpsiDenseBasis);
			  if (Options.DenseMethod==false) if(Opsi.N>Options.Settings.DenseBorder)
			  {
				WaveFunctionFree(&Opsi);
				goto EndOfUseDenseTest;
			  }
			  HTpsiDenseBasis=Opsi;
			  WaveFunctionNormalize(&HTpsiDenseBasis);
			}
			// T H T psi
			WaveFunctionCopy(HTpsiDenseBasis, &HTHTpsiDenseBasis);
			if((Options.Settings.restrictions2.N==0)?(OperatorPsi(T2[0], HTHTpsiDenseBasis, &TpsiSum)):(OperatorPsiRestricted(T2[0], HTHTpsiDenseBasis, &TpsiSum, Options.Settings.restrictions2)))
			{
			  WaveFunctionFree(&HTHTpsiDenseBasis);
			  WaveFunctionFree(&HTpsiDenseBasis);
			  goto EndOfUseDenseTest;
			}
			for(i=1;i<Options.NT2;i++)
			{
			  if((Options.Settings.restrictions2.N==0)?(OperatorPsi(T2[i], HTHTpsiDenseBasis, &Opsi)):(OperatorPsiRestricted(T2[i], HTHTpsiDenseBasis, &Opsi, Options.Settings.restrictions2)))
			  {
				WaveFunctionFree(&TpsiSum);
				WaveFunctionFree(&HTHTpsiDenseBasis);
				WaveFunctionFree(&HTpsiDenseBasis);
				goto EndOfUseDenseTest;
			  }
			  WaveFunctionAddUnionBasis(&TpsiSum, Opsi);
			  WaveFunctionFree(&Opsi);
			  if (Options.DenseMethod==false) if(TpsiSum.N>Options.Settings.DenseBorder)
			  {
				WaveFunctionFree(&TpsiSum);
				WaveFunctionFree(&HTHTpsiDenseBasis);
				WaveFunctionFree(&HTpsiDenseBasis);
				goto EndOfUseDenseTest;
			  }
			}
			WaveFunctionFree(&HTHTpsiDenseBasis);
			HTHTpsiDenseBasis=TpsiSum;
			WaveFunctionNormalize(&HTHTpsiDenseBasis);
			// H T H T psi
			len=0;
			while(len!=HTHTpsiDenseBasis.N)
			{
			  len=HTHTpsiDenseBasis.N;
			  if((Options.Settings.restrictions2.N==0)?(OperatorPsi(*H2, HTHTpsiDenseBasis, &Opsi)):(OperatorPsiRestricted(*H2, HTHTpsiDenseBasis, &Opsi, Options.Settings.restrictions2)))
			  {
				WaveFunctionFree(&HTHTpsiDenseBasis);
				WaveFunctionFree(&HTpsiDenseBasis);
				goto EndOfUseDenseTest;
			  }
			  WaveFunctionFree(&HTHTpsiDenseBasis);
			  if (Options.DenseMethod==false) if(Opsi.N>Options.Settings.DenseBorder)
			  {
				WaveFunctionFree(&Opsi);
				WaveFunctionFree(&HTpsiDenseBasis);
				goto EndOfUseDenseTest;
			  }
			  WaveFunctionNormalize(&Opsi);
			  HTHTpsiDenseBasis=Opsi;
			}
			// finsihed all tests and psiDenseBasis still fits
			Options.DenseMethod=true;
		  EndOfUseDenseTest:;
		  }
  if(Options.DenseMethod==true)
  {
    if(Options.Tensor==false)
    {
      ResonantSpectraDense(&HTpsiDenseBasis,&HTHTpsiDenseBasis,*H1, *H2, T1, Options.NT1, T2, Options.NT2, psi, Options.Npsi, &(Spec->Spec), Options.E01, Options.E02, &(Options.Settings));
      Spec->NSpectra=Options.NT1*Options.NT2*Options.Npsi*(Options.Settings.NE1+1);
    }
    else
    {
      return luaL_error(L,"Tensor form for Resonant spectra not yet available");
    }
  }
  else
  {
    // use sparse methods
    if(Options.Tensor==false)
    {
      for(ipsi=0;ipsi<Options.Npsi;ipsi++)
      {
        for(iT=0;iT<Options.NT1;iT++)
        {
          SettingsNow = Options.Settings;
          ResonantSpectrumWithSelfEnergy(*H1, *H2, T1[iT], T2, U, Fk, Options.NT2, psi[ipsi], &SpecTri, &SpecT1psi, Options.E01, Options.E02, &(SettingsNow),EnergyContinuum,gammaSE);
          Spec->NSpectra+=(SettingsNow.NE1+1)*Options.NT2;
          Spec->Spec= (SpectraType*)realloc(Spec->Spec, Spec->NSpectra*sizeof(SpectraType));
          for(i=0;i<(SettingsNow.NE1+1)*Options.NT2;i++)
          {
            FreeTriDiagonalMatrix(&(SpecTri[i]));
          }
          free(SpecTri);
          for(i=0;i<(SettingsNow.NE1+1)*Options.NT2;i++)
          {
            sprintf(SpecT1psi[i].Name, "%s T1 %i/%i psi %i/%i",SpecT1psi[i].Name,iT,Options.NT1,ipsi,Options.Npsi);
            Spec->Spec[i+Spec->NSpectra-(SettingsNow.NE1+1)*Options.NT2] = SpecT1psi[i];
          }
          free(SpecT1psi);
        }
      }
    }
    else
    {

      //ResonantSpectraTensor(*H1, *H2, T1, Options.NT1, T2, Options.NT2, psi, Options.Npsi, &SpecTri, &(Spec->Spec[0]), Options.E01, Options.E02, &Options.Settings);
      return luaL_error(L,"Tensor form for Resonant spectra not yet available");
    }
  }
  return 1;
}



uint32_t ResonantSpectrumWithSelfEnergy(const OperatorType OppH1, const OperatorType OppH2, const OperatorType T1, const OperatorType *T2,
										const vector<vector<OperatorType> > &U, const vector<vector<vector<double> > > &Fk, uint32_t NT2,
										const WaveFunctionType psi, TriDiagonalMatrixType **SpecTri2, SpectraType **Spec, double energy1,
										double energy2, ConvergenceSettingsType *settings, const vector<vector<double> > &ECont, double gammaSE)
{
  /*Given wavefunction psi (usually ground state, but could be excited)*/

  /*Calculates resonant spectra using <i|T*GTGT*GT|i> in the following steps:

   1. Calculate T|i>
   2. Use T|i> to transform H to tridiagonal form in Krylov basis
   3. Invert (w - H_tri + iGamma) to get G(w)
   4. Calculate G(w)T|i>
   5. Calculate T*G(w)T|i>
   6. Use the vector from step 5 to calculate RIXS spectrum using Lanczos
   */

  /* settings that need to be set:
   * NTri1, NTri2
   * Emin1, Emax1, NE1, Gamma1 - for automatic determination use DBL_MAX, -DBL_MAX, 0, 0
   * Emin2, Emax2, NE2, Gamma2 - for automatic determination use DBL_MAX, -DBL_MAX, 0, 0
   * epsilon
   */
  printf("Debug Mode: ResonantSpectrum\n");
  WaveFunctionType Tpsi,GTpsi,TGTpsi, PsiBasis[settings->NTri1];
  TriDiagonalMatrixType SpecTri1;
  MatrixType EFun;
  double EVal[settings->NTri1], EValInvRe[settings->NTri1], EValInvIm[settings->NTri1],tempDoubleRe,tempDoubleIm, EVal2[settings->NTri2], Norm;
  int i,p,k;
  double wRes;

  VerbosityPrintf(0x0001,"Start of ResonantSpectra using Lanczos\n");

  /* calcualte ground state energy */
  if(energy1==DBL_MAX)
  {
    PsiHermitianOperatorPsi(OppH1, psi, &energy1);
    WaveFunctionCalculateNormSqr(psi, &Norm);
    if(Norm > DBL_EPSILON) energy1/=Norm;
  }
  if(energy2==DBL_MAX)
  {
    PsiHermitianOperatorPsi(OppH2, psi, &energy2);
    WaveFunctionCalculateNormSqr(psi, &Norm);
    if(Norm > DBL_EPSILON) energy2/=Norm;
  }

  /* Calculate T|i> */
  if(settings->restrictions1.N==0)
  {
    OperatorPsi(T1,psi,&Tpsi);
  }
  else
  {
    OperatorPsiRestricted(T1, psi, &Tpsi, settings->restrictions1);
  }

  /* Use T|i> to transform H to tridiagonal form in Krylov basis */
  /* Using the version which gives back the Krylov basis vectors */
  SpecTri1.N = settings->NTri1;
  InitTriDiagonalMatrix(&SpecTri1);
  if(settings->restrictions1.N==0)
  {
    LanczosTriDiagonalizeKrylov(OppH1,Tpsi,settings->epsilon,&SpecTri1,PsiBasis);
  }
  else
  {
    LanczosTriDiagonalizeKrylovRestricted(OppH1, Tpsi, settings->restrictions1, settings->epsilon, &SpecTri1, PsiBasis);
  }

  /*Compute eigenvalues and eigenvectors, then S^-1 = QAQ^T where
   Q are the eigenvectors and A is a matrix with the eigenvalue inverses on the diagonal.
   */
  if(TriDiagonalMatrixToEigensystem(SpecTri1, &EFun, EVal))
  {
    printf("Problem with TridiagonalMatrixToEigenSystem\n");
  }

  /* determine Emin1, Emax1, NE1 and Gamma1 if not set */
  if(settings->Emin1>( 0.001*DBL_MAX) && (settings->Emax1<(-0.001*DBL_MAX)))
  {
    settings->Emin1=1.1*EVal[0]-0.1*EVal[settings->NTri1-1];
    settings->Emax1=1.1*EVal[settings->NTri1-1]-0.1*EVal[0];
    VerbosityPrintf(0x0002, "Automatically set energy range for intermediate states Emin=%22.15E Emax=%22.15E\n",settings->Emin1,settings->Emax1);
  }
  if(settings->Emin1>( 0.001*DBL_MAX))
  {
    settings->Emin1=1.1*EVal[0]-0.1*settings->Emax1;
    VerbosityPrintf(0x0002, "Automatically set energy range for intermediate states Emin=%22.15E Emax=%22.15E\n",settings->Emin1,settings->Emax1);
  }
  if(settings->Emax1<(-0.001*DBL_MAX))
  {
    settings->Emax1=1.1*EVal[settings->NTri1-1]-0.1*settings->Emin1;
    VerbosityPrintf(0x0002, "Automatically set energy range for intermediate states Emin=%22.15E Emax=%22.15E\n",settings->Emin1,settings->Emax1);
  }
  if((settings->Gamma1==0)&&(settings->NE1==0))
  {
    settings->NE1=1000;
    settings->Gamma1=5*(settings->Emax1-settings->Emin1)/settings->NE1;
    VerbosityPrintf(0x0002, "Automatically set number of energy points (NE=%d) and lifetime (Gamma=%22.15E) for intermediate state\n",settings->NE1,settings->Gamma1);
  }
  if(settings->NE1==0)
  {
    settings->NE1    =(uint32_t)(5*(settings->Emax1-settings->Emin1)/settings->Gamma1);
    VerbosityPrintf(0x0002, "Automatically set number of energy points (NE=%d) for intermediate state\n",settings->NE1);
  }
  if(settings->Gamma1==0)
  {
    settings->Gamma1 =(5*(settings->Emax1-settings->Emin1)/settings->NE1);
    VerbosityPrintf(0x0002, "Automatically set lifetime (Gamma=%22.15E) for intermediate state\n",settings->Gamma1);
  }

  /* allocate SpecTri2 */
  *SpecTri2 = (TriDiagonalMatrixType*)malloc((settings->NE1+1)*NT2*sizeof(TriDiagonalMatrixType));
  if(*SpecTri2==NULL)
  {
    printf("malloc of SpecTri2 failed in ResonantSpectra\n");
    return 1;
  }

// Martin start	-------------------------------------------------------------------------
	uint32_t Nw1 = settings->NE1+1;
	uint32_t Nw2 = settings->NE2+1;
	uint32_t Npoles2 = settings->NTri2;
	uint32_t Nkappa = (uint32_t) U.size();
	vector<vector<double> > poles2 = vector<vector<double> >(Nw1 * NT2);
	vector<vector<double> > residues2 = vector<vector<double> >(Nw1 * NT2);
	CompactMatrixType EigenVecs2;
	vector<double> T2G1T1psi = vector<double>(1);

	vector<vector<double> > polesSE = vector<vector<double> >(Nw1 * NT2 * Npoles2 * Nkappa);
	vector<CompactMatrixType> SelfEnergySMatrix = vector<CompactMatrixType>(Nw1 * NT2 * Npoles2 * Nkappa);
// Martin end	-------------------------------------------------------------------------

  for(p = 0; p <= settings->NE1; p++)
  {
    /* New method for GT|i>:
     1. Tridiagonalize H using T|i> (outside wres loop)
     2. Calculate eigenvalues (EVal) and eigenfunctions (EFun) of Htri (outside)
     3. For each wres
     a. EValInv[i] = 1/(wres - EVal[i] + iGammaRes)
     b. GT|i> = sum_j EFun[i][j]*EValIn[j]*K[j], where K is Krylov basis vectors
     */

    /* Calculate 1/(E_0 + wres - EVal[i] + iGammaRes), split into real and imaginary */

    wRes = settings->Emin1 + (settings->Emax1-settings->Emin1)*p/((double)settings->NE1) + energy1;
    for(i = 0; i < settings->NTri1; i++)
    {
      EValInvRe[i] =  (wRes - EVal[i])       / ((wRes - EVal[i])*(wRes - EVal[i])+(settings->Gamma1/2.0)*(settings->Gamma1/2.0));
      EValInvIm[i] = -(settings->Gamma1/2.0) / ((wRes - EVal[i])*(wRes - EVal[i])+(settings->Gamma1/2.0)*(settings->Gamma1/2.0));
    }

    /* We are trying to calculate Q[0,:]A^(-1)Q^(-1).  Because A is diagonal, we were
     able to calculate B = Q[0,:]*I*Q^(-1) outside the loop. Now we calculate diag(A)*B */
    /* Calculate Q*A^(-1)*Q^(-1)[:,0] */
    GTpsi.NFermion = psi.NFermion;
    GTpsi.NBoson   = psi.NBoson;
    GTpsi.QComplex = psi.QComplex;
    GTpsi.NBitsKey = psi.NBitsKey;
    if(WaveFunctionInit(&GTpsi))
    {
      printf("WaveFunctionInit failed in ResonantSpectra 00");
      return 1;
    }
    for(i = 0; i < settings->NTri1; i++)
    {
      tempDoubleRe = 0;
      tempDoubleIm = 0;

      if(EFun.QComplex)
      {
        for(k = 0; k < settings->NTri1; k++)
        {
          tempDoubleRe += EValInvRe[k]*EFun.VRe[k][i]*EFun.VRe[k][0] -
          EValInvIm[k]*EFun.VIm[k][i]*EFun.VRe[k][0] -
          EValInvIm[k]*EFun.VRe[k][i]*EFun.VIm[k][0] -
          EValInvRe[k]*EFun.VIm[k][i]*EFun.VIm[k][0];
          tempDoubleIm += EValInvIm[k]*EFun.VRe[k][i]*EFun.VRe[k][0] +
          EValInvRe[k]*EFun.VIm[k][i]*EFun.VRe[k][0] +
          EValInvRe[k]*EFun.VRe[k][i]*EFun.VIm[k][0] -
          EValInvIm[k]*EFun.VIm[k][i]*EFun.VIm[k][0];
        }
      }
      else
      {
        for(k = 0; k < settings->NTri1; k++)
        {
          tempDoubleRe += EValInvRe[k]*EFun.VRe[k][i]*EFun.VRe[k][0];
          tempDoubleIm += EValInvIm[k]*EFun.VRe[k][i]*EFun.VRe[k][0];
        }
      }
      WaveFunctionAddComplexMultipliedUnionBasis(&GTpsi, PsiBasis[i],SpecTri1.b[0]*tempDoubleRe,SpecTri1.b[0]*tempDoubleIm);
    }
    for(i = 0; i < NT2; i++)
    {
      /* Calculate T*G(w)T|i> */
      if(settings->restrictions2.N==0)
      {
        OperatorPsi(T2[i],GTpsi,&TGTpsi);
      }
      else
      {
        OperatorPsiRestricted(T2[i], GTpsi, &TGTpsi, settings->restrictions2);
      }
      /* Calculate the RIXS spectrum */
      (*SpecTri2)[p + i*(settings->NE1+1)].N=settings->NTri2;
      InitTriDiagonalMatrix(&((*SpecTri2)[p + i*(settings->NE1+1)]));

// Martin start		-----------------------------------------------------------------
      WaveFunctionType PhiBasis2[Npoles2], EigenBasis2[Npoles2];
      if(settings->restrictions2.N==0)
      {
        LanczosTriDiagonalizeKrylov(OppH2, TGTpsi, settings->epsilon, &((*SpecTri2)[p + i*(settings->NE1+1)]), PhiBasis2);
        /*for (uint32_t j = 0; j < Npoles2; j++)
        {
        	printf("PhiBasis2[%u].N = %u\n",j,PhiBasis2[j].N);
        }*/
      }
      else
      {
    	  printf("error: Restrictions for ResonantSpectrumWithSelfenergy not yet implementd \n");
    	  return -1;

    	  LanczosTriDiagonalizeRestricted(OppH2, TGTpsi, settings->restrictions2, settings->epsilon, &((*SpecTri2)[p + i*(settings->NE1+1)]), (sqrt(DBL_EPSILON)+settings->Gamma2));
      }
      poles2[i*Nw1 + p] = vector<double>(Npoles2);
      TriDiagonalMatrixToEigensystemCompact(&((*SpecTri2)[p + i*(settings->NE1+1)]),&(poles2[p + i*Nw1][0]),&EigenVecs2);
      //PrintArray(&poles2[i*Nw1+p][0],poles2[i*Nw1+p].size());
      //PrintMartix(EigenVecs2.V,EigenVecs2.Ni);
      double b0 = (*SpecTri2)[p + i*Nw1].b[0];
      T2G1T1psi[0] = b0;
      //printf("b0 = %.16E\n",b0);
      CalculateResidues(T2G1T1psi,&EigenVecs2,Npoles2,residues2[i*Nw1 + p]);
      RotateWaveFunctions(PhiBasis2,Npoles2,EigenVecs2,EigenBasis2);
      for (uint32_t j = 0; j < Npoles2; j++)
      {
    	  for (uint32_t k = 0; k < Nkappa; k++)
    	  {
    		  CalculateEnergyIndependentSelfEnergy(OppH2,EigenBasis2[j],U[k],SelfEnergySMatrix[(p*NT2*Npoles2 + i*Npoles2 + j)*Nkappa + k],polesSE[(p*NT2*Npoles2 + i*Npoles2 + j)*Nkappa + k]);
    	  }
    	  WaveFunctionFree(&EigenBasis2[j]);
    	  WaveFunctionFree(&PhiBasis2[j]);
      }
// Martin end		-----------------------------------------------------------------

    } //end of loop over i (T2)
    WaveFunctionFree(&GTpsi);
  } // end of loop over p (wres)


  /* create spectra from tridiagonal matrices */
  *Spec = (SpectraType*)malloc(NT2 * (settings->NE1+1) * sizeof(SpectraType));
  for(i=0; i<NT2 * (settings->NE1+1);i++)
  {
    /* determine min and max poles */
    TriDiagonalMatrixToEigenvalues((*SpecTri2)[i], EVal2);
    (*Spec)[i].EminPole = EVal2[0]-energy2;
    (*Spec)[i].EmaxPole = EVal2[settings->NTri2-1]-energy2;
  }
  /* determine the energy range of the spectrum if not give */
  EVal2[0]=(*Spec)[0].EminPole;
  EVal2[1]=(*Spec)[0].EmaxPole;
  for(i=1; i<NT2 * (settings->NE1+1);i++)
  {
    if((*Spec)[i].EminPole<EVal2[0]) EVal2[0]=(*Spec)[i].EminPole;
    if((*Spec)[i].EmaxPole>EVal2[1]) EVal2[1]=(*Spec)[i].EmaxPole;
  }
  if((settings->Emin2>( 0.001*DBL_MAX)) || (settings->Emax2<(-0.001*DBL_MAX)))
  {
    settings->Emin2 = 1.1 * EVal2[0] - 0.1*EVal2[1];
    settings->Emax2 = 1.1 * EVal2[1] - 0.1*EVal2[0];
    VerbosityPrintf(0x0002, "Automatically set energy range for final states Emin=%22.15E Emax=%22.15E\n",settings->Emin2,settings->Emax2);
  }
  if(settings->Emin2>( 0.001*DBL_MAX))
  {
    settings->Emin2=1.1*EVal2[0]-0.1*settings->Emax2;
    VerbosityPrintf(0x0002, "Automatically set energy range for final states Emin=%22.15E Emax=%22.15E\n",settings->Emin2,settings->Emax2);
  }
  if(settings->Emax2<(-0.001*DBL_MAX))
  {
    settings->Emax2=1.1*EVal2[1]-0.1*settings->Emin2;
    VerbosityPrintf(0x0002, "Automatically set energy range for final states Emin=%22.15E Emax=%22.15E\n",settings->Emin2,settings->Emax2);
  }
  if((settings->Gamma2==0)&&(settings->NE2==0))
  {
    settings->NE2=1000;
    settings->Gamma2=5*(settings->Emax2-settings->Emin2)/settings->NE2;
    VerbosityPrintf(0x0002, "Automatically set number of energy points (NE=%d) and lifetime (Gamma=%22.15E) for final state\n",settings->NE2,settings->Gamma2);
  }
  if(settings->NE2==0)
  {
    settings->NE2    =(uint32_t)(5*(settings->Emax2-settings->Emin2)/settings->Gamma2);
    VerbosityPrintf(0x0002, "Automatically set number of energy points (NE=%d) for final state\n",settings->NE2);
  }
  if(settings->Gamma2==0)
  {
    settings->Gamma2 =(5*(settings->Emax2-settings->Emin2)/settings->NE2);
    VerbosityPrintf(0x0002, "Automatically set lifetime (Gamma=%22.15E) for final state\n",settings->Gamma2);
  }
  vector<vector<mypoint> > SelfEnergiesEvaluated = vector<vector<mypoint> >(Npoles2);
  uint32_t Npoints2 = settings->NE2 + 1;
  double dE2 = (settings->Emax2-settings->Emin2)/settings->NE2;
  vector<mypoint> spectrum = vector<mypoint>(Npoints2);
  for (uint32_t widx = 0; widx < Npoints2; widx++)
  {
	  spectrum[widx].w = ((double)widx * dE2) + settings->Emin2; //TODO: + energy2;
  }
  for(uint32_t j = 0; j < Npoles2; j++)
  {
	  SelfEnergiesEvaluated[j] = vector<mypoint>(Npoints2);
  }
  if (!gammaSE) gammaSE = abs(ECont[0][0]-ECont[0][1])*6;
  for (uint32_t p = 0; p < Nw1; p++)
  {
	  for (uint32_t i = 0; i < NT2; i++)
	  {

			for(uint32_t j = 0; j < Npoles2; j++)
			{
			  for (uint32_t widx = 0; widx < Npoints2; widx++)
			  {
				  SelfEnergiesEvaluated[j][widx].w = ((double)widx * dE2) + settings->Emin2 + energy2;  //TODO: Check if that shift is correct + energy2
				  SelfEnergiesEvaluated[j][widx].y.Re = 0.;
				  SelfEnergiesEvaluated[j][widx].y.Im = 0.;
			  }
			  //PrintArray(&poles2[i*Nw1 + p][0],poles2[i*Nw1 + p].size());
			  //printf("w = %.2E\t%.2E\tenergy2 =%.2E\n",SelfEnergiesEvaluated[j][0].w,SelfEnergiesEvaluated[j][Npoints2-1].w,energy2);
			  for (uint32_t k = 0; k < Nkappa; k++)
			  {
				  uint32_t NContStates = (uint32_t) Fk[k].size();
				  for (uint32_t e = 0; e < NContStates; e++)
				  {
					  CalcSpectrumFromPolesAdd(0,gammaSE,ECont[k][e],polesSE[(p*NT2*Npoles2 + i*Npoles2 + j)*Nkappa + k],Fk[k][e],&(SelfEnergySMatrix[(p*NT2*Npoles2 + i*Npoles2 + j)*Nkappa + k]),SelfEnergiesEvaluated[j]);
					  //PrintArray(&polesSE[(p*NT2*Npoles2 + i*Npoles2 + j)*Nkappa + k][0],polesSE[(p*NT2*Npoles2 + i*Npoles2 + j)*Nkappa + k].size());
					  //int foo; std::cin >> foo;
				  }
			  }
			}
//			PrintMySpectrum(SelfEnergiesEvaluated[0]);
			for (uint32_t widx = 0; widx < Npoints2; widx++)
			{
				Complex16 z;
				z.Re = spectrum[widx].w  + energy2;
				z.Im = settings->Gamma2/2;
				spectrum[widx].y = Resolvent(z,poles2[i*Nw1 + p],residues2[i*Nw1 + p],SelfEnergiesEvaluated,widx);
			}

			(*Spec)[p + i*Nw1].N    =settings->NE2;
			(*Spec)[p + i*Nw1].Emin =settings->Emin2;
			(*Spec)[p + i*Nw1].Emax =settings->Emax2;
			(*Spec)[p + i*Nw1].dE   =(settings->Emax2-settings->Emin2)/settings->NE2;
			(*Spec)[p + i*Nw1].Gamma=settings->Gamma2;
			(*Spec)[p + i*Nw1].ResonantEnergy = settings->Emin1 + ((double)(p)/(settings->NE1)) * (settings->Emax1-settings->Emin1);
			//sprintf((*Spec)[i].Name, "Spectrum wres %i/%i T2 %i/%i",(int)(p),settings->NE1,(int)(i),NT2);
			SpectrumAllocate((*Spec)[p + i*Nw1], &(*Spec)[p + i*Nw1]);
			ConvertSpectrum(&(*Spec)[p + i*Nw1],spectrum);
	  }
  }
  return 0;
}



