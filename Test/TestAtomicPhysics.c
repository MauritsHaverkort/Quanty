uint32_t TestCommutationRelations()
{
  return 0;
}

uint32_t TestAtomicPhysics()
{
  uint32_t lshell,nelectron,NFermion,NExpMax=6;
  uint16_t LindexUp[7],LindexDn[7];
  int ml,i;
  TimeType Time;
  double FSlater[4],F0[4],F2[4],F4[4],F6[4],epsilon=1000*DBL_EPSILON;
  double ExpSxRe[NExpMax],ExpSyRe[NExpMax],ExpSzRe[NExpMax],ExpSsqrRe[NExpMax],ExpLxRe[NExpMax],ExpLyRe[NExpMax],ExpLzRe[NExpMax],ExpLsqrRe[NExpMax],ExpJxRe[NExpMax],ExpJyRe[NExpMax],ExpJzRe[NExpMax],ExpJsqrRe[NExpMax];
  double ExpSxIm[NExpMax],ExpSyIm[NExpMax],ExpSzIm[NExpMax],ExpSsqrIm[NExpMax],ExpLxIm[NExpMax],ExpLyIm[NExpMax],ExpLzIm[NExpMax],ExpLsqrIm[NExpMax],ExpJxIm[NExpMax],ExpJyIm[NExpMax],ExpJzIm[NExpMax],ExpJsqrIm[NExpMax];
  double ExpEnergyRe[NExpMax],ExpldotsRe[NExpMax],ExpFRe[NExpMax],ExpF0Re[NExpMax],ExpF2Re[NExpMax],ExpF4Re[NExpMax],ExpF6Re[NExpMax];
  double ExpEnergyIm[NExpMax],ExpldotsIm[NExpMax],ExpFIm[NExpMax],ExpF0Im[NExpMax],ExpF2Im[NExpMax],ExpF4Im[NExpMax],ExpF6Im[NExpMax];
  
  FSlater[0]=0;
  FSlater[1]=10.0;
  FSlater[2]=6.25;
  FSlater[3]=4.0;
  
  F0[0]=1.0; F2[0]=0.0; F4[0]=0.0; F6[0]=0.0;
  F0[1]=0.0; F2[1]=1.0; F4[1]=0.0; F6[1]=0.0;
  F0[2]=0.0; F2[2]=0.0; F4[2]=1.0; F6[2]=0.0;
  F0[3]=0.0; F2[3]=0.0; F4[3]=0.0; F6[3]=1.0;
  
  OperatorType Sx,Sy,Sz,Ssqr,Lx,Ly,Lz,Lsqr,Jx,Jy,Jz,Jsqr,ldots,U,U0,U2,U4,U6,H;
  WaveFunctionType psi, psiEigen[NExpMax];
  
  printf("This subroutine creates the L, S, J, ldots, and coulomb operator for a single atomic shell\n");
  printf("It then tests the eigen-states and the commutation relations for\n");
  printf("this is done succesively for the s, p, d and f shell, whith a filling of 0 to 4l+2\n");
  
  InitializeTiming(&Time);
  BeginTiming("Total",&Time);
  for(lshell=1;lshell<=1;lshell++)
  {
    BeginTiming("Shell",&Time);
    printf("************************************************************\n");
    printf("**********  Start with new angular momentum l=%2i  *********\n",lshell);
    printf("************************************************************\n");
    NFermion=4*lshell+2;
    for(ml=0;ml<NFermion/2;ml++)
    {
      LindexDn[ml]=2*ml;
      LindexUp[ml]=2*ml+1;
    }
    BeginTiming("Create Operators",&Time);
    CreateOperatorSx   (NFermion, NFermion/2, LindexUp, LindexDn, &Sx);
    CreateOperatorSy   (NFermion, NFermion/2, LindexUp, LindexDn, &Sy);
    CreateOperatorSz   (NFermion, NFermion/2, LindexUp, LindexDn, &Sz);
    CreateOperatorSsqr (NFermion, NFermion/2, LindexUp, LindexDn, &Ssqr);
    CreateOperatorLx   (NFermion, lshell,     LindexUp, LindexDn, &Lx);
    CreateOperatorLy   (NFermion, lshell,     LindexUp, LindexDn, &Ly);
    CreateOperatorLz   (NFermion, lshell,     LindexUp, LindexDn, &Lz);
    CreateOperatorLsqr (NFermion, lshell,     LindexUp, LindexDn, &Lsqr);
    CreateOperatorJx   (NFermion, lshell,     LindexUp, LindexDn, &Jx);
    CreateOperatorJy   (NFermion, lshell,     LindexUp, LindexDn, &Jy);
    CreateOperatorJz   (NFermion, lshell,     LindexUp, LindexDn, &Jz);
    CreateOperatorJsqr (NFermion, lshell,     LindexUp, LindexDn, &Jsqr);
    CreateOperatorldots(NFermion, lshell,     LindexUp, LindexDn, &ldots);
    CreateOperatorCoulombSingleShell(NFermion, lshell, LindexUp, LindexDn, FSlater, &U);
    CreateOperatorCoulombSingleShell(NFermion, lshell, LindexUp, LindexDn, F0,     &U0);
    CreateOperatorCoulombSingleShell(NFermion, lshell, LindexUp, LindexDn, F2,     &U2);
    CreateOperatorCoulombSingleShell(NFermion, lshell, LindexUp, LindexDn, F4,     &U4);
    CreateOperatorCoulombSingleShell(NFermion, lshell, LindexUp, LindexDn, F6,     &U6);
    OperatorCopy(&H, U);
    OperatorAdd(&H, Sz, 0.000002, 0);
    OperatorAdd(&H, Lz, 0.000001, 0);
    OperatorAdd(&H, ldots, 0.01, 0);
    EndTiming("Create Operators",&Time);
    
    PrintOperatorReadable(Sx);
    PrintOperatorReadable(Sy);
    PrintOperatorReadable(Sz);
    PrintOperatorReadable(Ssqr);
    PrintOperatorReadable(Lx);
    PrintOperatorReadable(Ly);
    PrintOperatorReadable(Lz);
    PrintOperatorReadable(Lsqr);
    PrintOperatorReadable(Jx);
    PrintOperatorReadable(Jy);
    PrintOperatorReadable(Jz);
    PrintOperatorReadable(Jsqr);
    PrintOperatorReadable(ldots);
    PrintOperatorReadable(U);
    PrintOperatorReadable(U0);
    PrintOperatorReadable(U2);
    PrintOperatorReadable(U4);
    PrintOperatorReadable(U6);
    
    
    
    
    for(nelectron=0;nelectron<=4*lshell+2;nelectron++)
    {
      BeginTiming("nelectrons",&Time);
      printf("************************************************************\n");
      printf("**********    Start with %2i electrons in shell    *********\n",nelectron);
      printf("************************************************************\n");
      
      CreateRealRandomWaveFunction(&psi, NFermion, nelectron);
      
      LanczosEigenSystem(H, psi, psiEigen, &epsilon, 1000*DBL_EPSILON, Min(NExpMax,BinomialCoefficient(NFermion,nelectron)));
      
      for(i=0;i<Min(NExpMax,BinomialCoefficient(NFermion,nelectron));i++)
      {
        PsiOperatorPsi(Sx  , psiEigen[i], &ExpSxRe[i],   &ExpSxIm[i]);
        PsiOperatorPsi(Sy  , psiEigen[i], &ExpSyRe[i],   &ExpSyIm[i]);
        PsiOperatorPsi(Sz  , psiEigen[i], &ExpSzRe[i],   &ExpSzIm[i]);
        PsiOperatorPsi(Ssqr, psiEigen[i], &ExpSsqrRe[i], &ExpSsqrIm[i]);
        PsiOperatorPsi(Lx  , psiEigen[i], &ExpLxRe[i],   &ExpLxIm[i]);
        PsiOperatorPsi(Ly  , psiEigen[i], &ExpLyRe[i],   &ExpLyIm[i]);
        PsiOperatorPsi(Lz  , psiEigen[i], &ExpLzRe[i],   &ExpLzIm[i]);
        PsiOperatorPsi(Lsqr, psiEigen[i], &ExpLsqrRe[i], &ExpLsqrIm[i]);
        PsiOperatorPsi(Jx  , psiEigen[i], &ExpJxRe[i],   &ExpJxIm[i]);
        PsiOperatorPsi(Jy  , psiEigen[i], &ExpJyRe[i],   &ExpJyIm[i]);
        PsiOperatorPsi(Jz  , psiEigen[i], &ExpJzRe[i],   &ExpJzIm[i]);
        PsiOperatorPsi(Jsqr, psiEigen[i], &ExpJsqrRe[i], &ExpJsqrIm[i]);
        PsiOperatorPsi(ldots, psiEigen[i], &ExpldotsRe[i], &ExpldotsIm[i]);
        PsiOperatorPsi(U , psiEigen[i], &ExpFRe[i],  &ExpFIm[i]);
        PsiOperatorPsi(U0, psiEigen[i], &ExpF0Re[i], &ExpF0Im[i]);
        PsiOperatorPsi(U2, psiEigen[i], &ExpF2Re[i], &ExpF2Im[i]);
        PsiOperatorPsi(U4, psiEigen[i], &ExpF4Re[i], &ExpF4Im[i]);
        PsiOperatorPsi(U6, psiEigen[i], &ExpF6Re[i], &ExpF6Im[i]);
        PsiOperatorPsi(H, psiEigen[i], &ExpEnergyRe[i], &ExpEnergyIm[i]);
      }
      printf("   Energy          F0              F2              F4              F6              Sz              Ssqr            Lz              Lsqr            Jz              Jsqr\n");
      for(i=0;i<Min(NExpMax,BinomialCoefficient(NFermion,nelectron));i++)
      {
        printf("%15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E\n",ExpEnergyRe[i],ExpF0Re[i],ExpF2Re[i],ExpF4Re[i],ExpF6Re[i],ExpSzRe[i],ExpSsqrRe[i],ExpLzRe[i],ExpLsqrRe[i],ExpJzRe[i],ExpJsqrRe[i]);
      }
      
      
      
      EndTiming("nelectrons",&Time);
      PrintTiming(Time);
    }
    
    FreeOperator(&Sx);
    FreeOperator(&Sy);
    FreeOperator(&Sz);
    FreeOperator(&Ssqr);
    FreeOperator(&Lx);
    FreeOperator(&Ly);
    FreeOperator(&Lz);
    FreeOperator(&Lsqr);
    FreeOperator(&Jx);
    FreeOperator(&Jy);
    FreeOperator(&Jz);
    FreeOperator(&Jsqr);
    FreeOperator(&ldots);
    FreeOperator(&U);
    
    EndTiming("Shell",&Time);
    PrintTiming(Time);
  }
  EndTiming("Total",&Time);
  PrintTiming(Time);
  
  
  return 0;
}