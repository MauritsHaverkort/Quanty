void BasisOptimizationSmallLowSym();
void BasisOptimizationSmallD4h();

void BasisOptimization()
{
  BasisOptimizationSmallD4h();
  //BasisOptimizationSmallLowSym();
}

double TestPotentialD4h(uint32_t i, uint32_t j)
{
  if((i==j)&&((i==0)||(i==2))) return  0.123;
  if((i==j)&&((i==1)))         return -0.234;
  if((i==j)&&((i==3)||(i==5))) return -1.234;
  if((i==j)&&((i==4)))         return  2.634;
  if((i==j)&&((i==6)||(i==8))) return -0.554;
  if((i==j)&&((i==7)))         return  1.143;
  if(((i==j+3)||(i+3==j))&&((i==0)||(i==2)))         return  0.812;
  if(((i==j+3)||(i+3==j))&&((i==1)))                 return  2.345;
  if(((i==j+6)||(i+6==j))&&((i==0)||(i==2)))         return  0.534;
  if(((i==j+6)||(i+6==j))&&((i==1)))                 return  3.012;
  return 0;
}

void BasisOptimizationSmallD4h()
{
  double delta1=1,delta2=2,U=5,Fpp2=15,Fpp0=U+4*Fpp2/50,tmax=2;
  double Ep3p6p6,Ep4p5p6,Ep4p6p5,Ep5p5p5,Ep5p4p6,Ep5p6p4,Ep1,Ep2,Ep3;
  double FSlater[]={Fpp0,Fpp2};
  MatrixType H1,H1New;
  uint32_t i,j;
  OperatorType Sx, Sy, Sz, Ssqr, Lx, Ly, Lz, Lsqr, Jx, Jy, Jz, Jsqr, OppH1, OppH1New, Coulomb, H, HNew, temp;
  WaveFunctionType psiGrd[12],psiGrdNew[12];
  double EGrd[12],VarEGrd[12],ECoulomb[12],EKinetic[12],accuracy,ValRho[18];
  double EGrdNew[12],VarEGrdNew[12],ECoulombNew[12],EKineticNew[12];
  MatrixType Rho,TRho;
  uint32_t NFermion=18, NElectron=18-3;
  uint16_t LindexUpp1[]={0,2,4}, LindexDnp1[]={1,3,5}; /* LindexUpp2[]={6,8,10}, LindexDnp2[]={7,9,11}, LindexUpp3[]={12,14,16}, LindexDnp3[]={13,15,17}; */
  
  /* create the one electron interactions */
  H1.Ni=18;
  H1.Nj=18;
  H1.QComplex=0;
  InitMatrix(&H1);
  strcpy(H1.Name, "H_{Original one particle Hopping}");
  for(i=0;i<H1.Ni;i+=2)
  {
    for(j=i;j<H1.Nj;j+=2)
    {
      H1.VRe[i][j]=TestPotentialD4h(i/2,j/2);
      H1.VRe[j][i]=H1.VRe[i][j];
      H1.VRe[j+1][i+1]=H1.VRe[i][j];
      H1.VRe[i+1][j+1]=H1.VRe[i][j];
    }
  }
  /* set center of each shell to zero */
  Ep1=0;Ep2=0;Ep3=0;
  for(i=0;i<6;i++)
  {
    Ep1+=H1.VRe[i][i];
  }
  Ep1/=6.0;
  for(i=0;i<6;i++)
  {
    H1.VRe[i][i]-=Ep1;
  }
  for(i=6;i<12;i++)
  {
    Ep2+=H1.VRe[i][i];
  }
  Ep2/=6.0;
  for(i=6;i<12;i++)
  {
    H1.VRe[i][i]-=Ep2;
  }
  for(i=12;i<18;i++)
  {
    Ep3+=H1.VRe[i][i];
  }
  Ep3/=6.0;
  for(i=12;i<18;i++)
  {
    H1.VRe[i][i]-=Ep3;
  }
  /* repair differenece between onsite energy and delta when U is included */
  Ep1=0;Ep2=0;Ep3=0;
  for(i=0;i<6;i++)
  {
    H1.VRe[i][i]+=(1/5.0)*(-13*U+2*delta1+2*delta2);
    Ep1+=H1.VRe[i][i];
  }
  Ep1/=6.0;
  for(i=6;i<12;i++)
  {
    H1.VRe[i][i]+=(1/5.0)*(2*U+2*delta1-3*delta2);
    Ep2+=H1.VRe[i][i];
  }
  Ep2/=6.0;
  for(i=12;i<18;i++)
  {
    H1.VRe[i][i]+=(1/5.0)*(2*U-3*delta1+2*delta2);
    Ep3+=H1.VRe[i][i];
  }
  Ep3/=6.0;
  /* printout of the settings */
  printf("\nCalculation of 3 p-shells interacting with each other by hopping\n");
  printf("The first p-shell has coulomb interactions\n");
  printf("There are 3 holes in the system\n");
  printf("The interactions are:");
  printf("U       : %lf\n",U);
  printf("Fpp0    : %lf\n",Fpp0);
  printf("Fpp2    : %lf\n",Fpp2);
  printf("tmax    : %lf\n",tmax);
  printf("delta1  : %lf\n",delta1);
  printf("delta2  : %lf\n",delta2);
  printf("The resulting one particle Hopping matrix is:\n");
  PrintMatrix(H1);
  printf("Average energy of shell 1: %lf\n",Ep1);
  printf("Average energy of shell 2: %lf\n",Ep2);
  printf("Average energy of shell 3: %lf\n",Ep3);
  Ep3p6p6=3*Ep1+6*Ep2+6*Ep3+(3*U);
  Ep4p5p6=4*Ep1+5*Ep2+6*Ep3+(6*U);
  Ep4p6p5=4*Ep1+6*Ep2+5*Ep3+(6*U);
  Ep5p5p5=5*Ep1+5*Ep2+5*Ep3+(10*U);
  Ep5p6p4=5*Ep1+6*Ep2+4*Ep3+(10*U);
  Ep5p4p6=5*Ep1+4*Ep2+6*Ep3+(10*U);
  printf("Average energy of configuration p3-p6-p6: %lf\n",Ep3p6p6);
  printf("Average energy of configuration p4-p6-p5: %lf\n",Ep4p6p5);
  printf("Average energy of configuration p4-p5-p6: %lf\n",Ep4p5p6);
  printf("Average energy of configuration p5-p5-p5: %lf\n",Ep5p5p5);
  printf("Average energy of configuration p5-p6-p4: %lf\n",Ep5p6p4);
  printf("Average energy of configuration p5-p4-p6: %lf\n",Ep5p4p6);
  /* create operators */
  CreateOperatorSx(NFermion, 3, LindexUpp1, LindexDnp1, &Sx);
  CreateOperatorSy(NFermion, 3, LindexUpp1, LindexDnp1, &Sy);
  CreateOperatorSz(NFermion, 3, LindexUpp1, LindexDnp1, &Sz);
  CreateOperatorSsqr(NFermion, 3, LindexUpp1, LindexDnp1, &Ssqr);
  CreateOperatorLx(NFermion, 1, LindexUpp1, LindexDnp1, &Lx);
  CreateOperatorLy(NFermion, 1, LindexUpp1, LindexDnp1, &Ly);
  CreateOperatorLz(NFermion, 1, LindexUpp1, LindexDnp1, &Lz);
  CreateOperatorLsqr(NFermion, 1, LindexUpp1, LindexDnp1, &Lsqr);
  CreateOperatorJx(NFermion, 1, LindexUpp1, LindexDnp1, &Jx);
  CreateOperatorJy(NFermion, 1, LindexUpp1, LindexDnp1, &Jy);
  CreateOperatorJz(NFermion, 1, LindexUpp1, LindexDnp1, &Jz);
  CreateOperatorJsqr(NFermion, 1, LindexUpp1, LindexDnp1, &Jsqr);
  CreateOperatorCoulombSingleShell(NFermion, 1, LindexUpp1, LindexDnp1, FSlater, &Coulomb);
  InitOperatorCopySettings(&OppH1, Sz);
  strcpy(OppH1.Name, "H_{One electron}");
  for(i=0;i<NFermion;i++)
  {
    for(j=0;j<NFermion;j++)
    {
      CreateOperatorNumber(NFermion, i, j, &temp);
      OperatorAdd(&OppH1,temp,H1.VRe[i][j],0);
      FreeOperator(&temp);
    }
  }
  OperatorChop(&OppH1, DBL_EPSILON);
  OperatorCopy(&H, Coulomb);
  OperatorAdd(&H, Sz, 0.1, 0);
  OperatorAdd(&H, Sx, 0.00001, 0);
  OperatorAdd(&H, OppH1, 1, 0);
  strcpy(H.Name, "Hamiltonian");
  /* Print operators */
  printf("\n\n==================================\nPrint of operators\n==================================\n\n");
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
  PrintOperatorReadable(OppH1);
  PrintOperatorReadable(Coulomb);
  PrintOperatorReadable(H);
  /* loop over accuracies */
  for(i=0;i<12;i++)
  {
    accuracy=1/pow(10, (double)(i+2));
    printf("\nCalculating the ground-state with truncation of basis states smaller than : %22.15E\n",accuracy);
    /* create starting point wavefunction */
    CreateRealRandomWaveFunction(&psiGrd[i], NFermion, NElectron);
    /* Find groundstate */
    LanczosGroundState(H, &(psiGrd[i]), &accuracy, 0.000000000000001, 0, &(psiGrd[i]), 0);
    
    PrintWaveFunctionReadable(psiGrd[i]);
    
    PsiHermitianOperatorPsi(H, psiGrd[i], &EGrd[i]);
    PsiHermitianOperatorPsi(Coulomb, psiGrd[i], &ECoulomb[i]);
    PsiHermitianOperatorPsi(OppH1, psiGrd[i], &EKinetic[i]);
    PsiHermitianOperatorSqrPsi(H, psiGrd[i], &VarEGrd[i]);
    VarEGrd[i]-=Sqr(EGrd[i]);
    printf("\n E                      VarE                   ECoulomb               Ekinetic                   Nstates\n");
    printf("%22.15E %22.15E %22.15E %22.15E %10i\n",EGrd[i],VarEGrd[i],ECoulomb[i],EKinetic[i],psiGrd[i].N);
  }
  /* print an overview of groundstates */
  printf("\n\n==================================\nOriginal Hamiltonian\n==================================\n\n");
  printf("\n E                      VarE                   ECoulomb               Ekinetic                   Nstates       Accuracy\n");
  for(i=0;i<12;i++)
  {
    accuracy=1/pow(10, (double)(i+2));
    printf("%22.15E %22.15E %22.15E %22.15E %10i %22.15E\n",EGrd[i],VarEGrd[i],ECoulomb[i],EKinetic[i],psiGrd[i].N,accuracy);
  }
  /* Density matrix */
  CalculateDensityMatrix(&Rho, psiGrd[4]);
  PrintMatrix(Rho);
  /* create rotation of basis */
  for(i=0;i<18;i++)
  {
    for(j=0;j<6;j++)
    {
      Rho.VRe[i][j]=0;
      Rho.VRe[j][i]=0;
    }
  }
  for(i=0;i<6;i++)
  {
    Rho.VRe[i][i]=-1;
  }
  Eigensystem(Rho, &TRho, ValRho);
  printf("Eigenvalues of density matrix (occupations in optimized basis)\n");
  for(i=0;i<18;i++)
  {
    printf("%3i %22.15E\n",i,ValRho[i]);
  }
  MatrixTransform(&H1New, H1, TRho);
  InitOperatorCopySettings(&OppH1New, Sz);
  strcpy(OppH1New.Name, "H_{One electron new}");
  for(i=0;i<NFermion;i++)
  {
    for(j=0;j<NFermion;j++)
    {
      CreateOperatorNumber(NFermion, i, j, &temp);
      OperatorAdd(&OppH1New,temp,H1New.VRe[i][j],0);
      FreeOperator(&temp);
    }
  }
  OperatorChop(&OppH1New, DBL_EPSILON);
  OperatorCopy(&HNew, Coulomb);
  OperatorAdd(&HNew, Sz, 0.1, 0);
  OperatorAdd(&HNew, Sx, 0.00001, 0);
  OperatorAdd(&HNew, OppH1New, 1, 0);
  strcpy(HNew.Name, "Hamiltonian rotated basis");
  /* print the new operators */
  PrintOperatorReadable(OppH1New);
  PrintOperatorReadable(HNew);
  /* loop over accuracies */
  for(i=0;i<12;i++)
  {
    accuracy=1/pow(10, (double)(i+2));
    printf("\n\n\nCalculating the ground-state with truncation of basis states smaller than : %22.15E\n",accuracy);
    /* create starting point wavefunction */
    CreateRealRandomWaveFunction(&psiGrdNew[i], NFermion, NElectron);
    /* Find groundstate */
    LanczosGroundState(HNew, &(psiGrdNew[i]), &accuracy, 0.000000000000001, 0, &(psiGrdNew[i]), 0);
    WaveFunctionRandomize(&(psiGrdNew[i]));
    LanczosGroundState(HNew, &(psiGrdNew[i]), &accuracy, 0.000000000000001, 0, &(psiGrdNew[i]), 0);
    PsiHermitianOperatorPsi(HNew, psiGrdNew[i], &EGrdNew[i]);
    PsiHermitianOperatorPsi(Coulomb, psiGrdNew[i], &ECoulombNew[i]);
    PsiHermitianOperatorPsi(OppH1New, psiGrdNew[i], &EKineticNew[i]);
    PsiHermitianOperatorSqrPsi(HNew, psiGrdNew[i], &VarEGrdNew[i]);
    VarEGrdNew[i]-=Sqr(EGrdNew[i]);
    printf("\n E                      VarE                   ECoulomb               Ekinetic                   Nstates\n");
    printf("%22.15E %22.15E %22.15E %22.15E %10i\n",EGrdNew[i],VarEGrdNew[i],ECoulombNew[i],EKineticNew[i],psiGrdNew[i].N);
  }
  /* print an overview of groundstates */
  printf("\n\n==================================\nRotated by all configurations included\n==================================\n\n");
  printf("\n E                      VarE                   VarE_old               Enew-Eold                  Nstates      Nstatesold   ccuracy\n");
  for(i=0;i<12;i++)
  {
    accuracy=1/pow(10, (double)(i+2));
    printf("%22.15E %22.15E %22.15E %22.15E %10i %10i %22.15E\n",EGrdNew[i],VarEGrdNew[i],VarEGrd[i],EGrdNew[i]-EGrd[i],psiGrdNew[i].N,psiGrd[i].N,accuracy);
  }
  
  
  /* remove second shell from H */
  for(i=0;i<12;i++)
  {
    WaveFunctionFree(&psiGrdNew[i]);
  }
  FreeOperator(&OppH1New);
  FreeOperator(&HNew);
  for(i=12;i<18;i++)
  {
    for(j=0;j<12;j++)
    {
      H1New.VRe[i][j]=0;
      H1New.VRe[j][i]=0;
    }
  }
  for(i=12;i<18;i++)
  {
    H1New.VRe[i][i]-=2;
  }
  InitOperatorCopySettings(&OppH1New, Sz);
  strcpy(OppH1New.Name, "H_{One electron new}");
  for(i=0;i<NFermion;i++)
  {
    for(j=0;j<NFermion;j++)
    {
      CreateOperatorNumber(NFermion, i, j, &temp);
      OperatorAdd(&OppH1New,temp,H1New.VRe[i][j],0);
      FreeOperator(&temp);
    }
  }
  OperatorChop(&OppH1New, DBL_EPSILON);
  OperatorCopy(&HNew, Coulomb);
  OperatorAdd(&HNew, Sz, 0.1, 0);
  OperatorAdd(&HNew, Sx, 0.00001, 0);
  OperatorAdd(&HNew, OppH1New, 1, 0);
  strcpy(HNew.Name, "Hamiltonian rotated basis");
  /* print the new operators */
  PrintOperatorReadable(OppH1New);
  PrintOperatorReadable(HNew);
  /* loop over accuracies */
  for(i=0;i<12;i++)
  {
    accuracy=1/pow(10, (double)(i+2));
    printf("\n\n\nCalculating the ground-state with truncation of basis states smaller than : %22.15E\n",accuracy);
    /* create starting point wavefunction */
    CreateRealRandomWaveFunction(&psiGrdNew[i], NFermion, NElectron);
    /* Find groundstate */
    LanczosGroundState(HNew, &(psiGrdNew[i]), &accuracy, 0.000000000000001, 0, &(psiGrdNew[i]), 0);
    WaveFunctionRandomize(&(psiGrdNew[i]));
    LanczosGroundState(HNew, &(psiGrdNew[i]), &accuracy, 0.000000000000001, 0, &(psiGrdNew[i]), 0);
    PsiHermitianOperatorPsi(HNew, psiGrdNew[i], &EGrdNew[i]);
    PsiHermitianOperatorPsi(Coulomb, psiGrdNew[i], &ECoulombNew[i]);
    PsiHermitianOperatorPsi(OppH1New, psiGrdNew[i], &EKineticNew[i]);
    PsiHermitianOperatorSqrPsi(HNew, psiGrdNew[i], &VarEGrdNew[i]);
    VarEGrdNew[i]-=Sqr(EGrdNew[i]);
    EGrdNew[i]+=12;
    printf("\n E                      VarE                   ECoulomb               Ekinetic                   Nstates\n");
    printf("%22.15E %22.15E %22.15E %22.15E %10i\n",EGrdNew[i],VarEGrdNew[i],ECoulombNew[i],EKineticNew[i],psiGrdNew[i].N);
  }
  /* print an overview of groundstates */
  printf("\n\n==================================\nRotated one shell less\n==================================\n\n");
  printf("\n E                      VarE                   VarE_old               Enew-Eold                  Nstates      Nstatesold   ccuracy\n");
  for(i=0;i<12;i++)
  {
    accuracy=1/pow(10, (double)(i+2));
    printf("%22.15E %22.15E %22.15E %22.15E %10i %10i %22.15E\n",EGrdNew[i],VarEGrdNew[i],VarEGrd[i],EGrdNew[i]-EGrd[i],psiGrdNew[i].N,psiGrd[i].N,accuracy);
  }
}

void BasisOptimizationSmallLowSym()
{
  double delta1=1,delta2=0.5,U=5,Fpp2=15,Fpp0=U+4*Fpp2/50,tmax=1;
  double Ep3p6p6,Ep4p5p6,Ep4p6p5,Ep5p5p5,Ep5p4p6,Ep5p6p4,Ep1,Ep2,Ep3;
  double FSlater[]={Fpp0,Fpp2};
  MatrixType H1,H1New;
  uint32_t i,j;
  OperatorType Sx, Sy, Sz, Ssqr, Lx, Ly, Lz, Lsqr, Jx, Jy, Jz, Jsqr, OppH1, OppH1New, Coulomb, H, HNew, temp;
  WaveFunctionType psiGrd[12],psiGrdNew[12];
  double EGrd[12],VarEGrd[12],ECoulomb[12],EKinetic[12],accuracy,ValRho[18];
  double EGrdNew[12],VarEGrdNew[12],ECoulombNew[12],EKineticNew[12];
  MatrixType Rho,TRho;
  uint32_t NFermion=18, NElectron=18-3;
  uint16_t LindexUpp1[]={0,2,4}, LindexDnp1[]={1,3,5}; /* LindexUpp2[]={6,8,10}, LindexDnp2[]={7,9,11}, LindexUpp3[]={12,14,16}, LindexDnp3[]={13,15,17}; */
  
  /* create the one electron interactions */
  H1.Ni=18;
  H1.Nj=18;
  H1.QComplex=0;
  InitMatrix(&H1);
  strcpy(H1.Name, "H_{Original one particle Hopping}");
  for(i=0;i<H1.Ni;i+=2)
  {
    for(j=i;j<H1.Nj;j+=2)
    {
      H1.VRe[i][j]=2*tmax*(rand()-RAND_MAX/2)/((double)RAND_MAX);
      H1.VRe[j][i]=H1.VRe[i][j];
      H1.VRe[j+1][i+1]=H1.VRe[i][j];
      H1.VRe[i+1][j+1]=H1.VRe[i][j];
    }
  }
  /* set center of each shell to zero */
  Ep1=0;Ep2=0;Ep3=0;
  for(i=0;i<6;i++)
  {
    Ep1+=H1.VRe[i][i];
  }
  Ep1/=6.0;
  for(i=0;i<6;i++)
  {
    H1.VRe[i][i]-=Ep1;
  }
  for(i=6;i<12;i++)
  {
    Ep2+=H1.VRe[i][i];
  }
  Ep2/=6.0;
  for(i=6;i<12;i++)
  {
    H1.VRe[i][i]-=Ep2;
  }
  for(i=12;i<18;i++)
  {
    Ep3+=H1.VRe[i][i];
  }
  Ep3/=6.0;
  for(i=12;i<18;i++)
  {
    H1.VRe[i][i]-=Ep3;
  }
  /* repair differenece between onsite energy and delta when U is included */
  Ep1=0;Ep2=0;Ep3=0;
  for(i=0;i<6;i++)
  {
    H1.VRe[i][i]+=(1/5.0)*(-13*U+2*delta1+2*delta2);
    Ep1+=H1.VRe[i][i];
  }
  Ep1/=6.0;
  for(i=6;i<12;i++)
  {
    H1.VRe[i][i]+=(1/5.0)*(2*U+2*delta1-3*delta2);
    Ep2+=H1.VRe[i][i];
  }
  Ep2/=6.0;
  for(i=12;i<18;i++)
  {
    H1.VRe[i][i]+=(1/5.0)*(2*U-3*delta1+2*delta2);
    Ep3+=H1.VRe[i][i];
  }
  Ep3/=6.0;
  /* printout of the settings */
  printf("\nCalculation of 3 p-shells interacting with each other by hopping\n");
  printf("The first p-shell has coulomb interactions\n");
  printf("There are 3 holes in the system\n");
  printf("The interactions are:");
  printf("U       : %lf\n",U);
  printf("Fpp2    : %lf\n",Fpp0);
  printf("Fpp2    : %lf\n",Fpp2);
  printf("tmax    : %lf\n",tmax);
  printf("delta1  : %lf\n",delta1);
  printf("delta2  : %lf\n",delta2);
  printf("The resulting one particle Hopping matrix is:\n");
  PrintMatrix(H1);
  printf("Average energy of shell 1: %lf\n",Ep1);
  printf("Average energy of shell 2: %lf\n",Ep2);
  printf("Average energy of shell 3: %lf\n",Ep3);
  Ep3p6p6=3*Ep1+6*Ep2+6*Ep3+(3*U);
  Ep4p5p6=4*Ep1+5*Ep2+6*Ep3+(6*U);
  Ep4p6p5=4*Ep1+6*Ep2+5*Ep3+(6*U);
  Ep5p5p5=5*Ep1+5*Ep2+5*Ep3+(10*U);
  Ep5p6p4=5*Ep1+6*Ep2+4*Ep3+(10*U);
  Ep5p4p6=5*Ep1+4*Ep2+6*Ep3+(10*U);
  printf("Average energy of configuration p3-p6-p6: %lf\n",Ep3p6p6);
  printf("Average energy of configuration p4-p6-p5: %lf\n",Ep4p6p5);
  printf("Average energy of configuration p4-p5-p6: %lf\n",Ep4p5p6);
  printf("Average energy of configuration p5-p5-p5: %lf\n",Ep5p5p5);
  printf("Average energy of configuration p5-p6-p4: %lf\n",Ep5p6p4);
  printf("Average energy of configuration p5-p4-p6: %lf\n",Ep5p4p6);
  /* create operators */
  CreateOperatorSx(NFermion, 3, LindexUpp1, LindexDnp1, &Sx);
  CreateOperatorSy(NFermion, 3, LindexUpp1, LindexDnp1, &Sy);
  CreateOperatorSz(NFermion, 3, LindexUpp1, LindexDnp1, &Sz);
  CreateOperatorSsqr(NFermion, 3, LindexUpp1, LindexDnp1, &Ssqr);
  CreateOperatorLx(NFermion, 1, LindexUpp1, LindexDnp1, &Lx);
  CreateOperatorLy(NFermion, 1, LindexUpp1, LindexDnp1, &Ly);
  CreateOperatorLz(NFermion, 1, LindexUpp1, LindexDnp1, &Lz);
  CreateOperatorLsqr(NFermion, 1, LindexUpp1, LindexDnp1, &Lsqr);
  CreateOperatorJx(NFermion, 1, LindexUpp1, LindexDnp1, &Jx);
  CreateOperatorJy(NFermion, 1, LindexUpp1, LindexDnp1, &Jy);
  CreateOperatorJz(NFermion, 1, LindexUpp1, LindexDnp1, &Jz);
  CreateOperatorJsqr(NFermion, 1, LindexUpp1, LindexDnp1, &Jsqr);
  CreateOperatorCoulombSingleShell(NFermion, 1, LindexUpp1, LindexDnp1, FSlater, &Coulomb);
  InitOperatorCopySettings(&OppH1, Sz);
  strcpy(OppH1.Name, "H_{One electron}");
  for(i=0;i<NFermion;i++)
  {
    for(j=0;j<NFermion;j++)
    {
      CreateOperatorNumber(NFermion, i, j, &temp);
      OperatorAdd(&OppH1,temp,H1.VRe[i][j],0);
      FreeOperator(&temp);
    }
  }
  OperatorChop(&OppH1, DBL_EPSILON);
  OperatorCopy(&H, Coulomb);
  OperatorAdd(&H, Sz, 0.1, 0);
  OperatorAdd(&H, Sx, 0.00001, 0);
  OperatorAdd(&H, OppH1, 1, 0);
  strcpy(H.Name, "Hamiltonian");
  /* Print operators */
  printf("\n\n==================================\nPrint of operators\n==================================\n\n");
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
  PrintOperatorReadable(OppH1);
  PrintOperatorReadable(Coulomb);
  PrintOperatorReadable(H);
  /* loop over accuracies */
  for(i=0;i<12;i++)
  {
    accuracy=1/pow(10, (double)(i+2));
    printf("\nCalculating the ground-state with truncation of basis states smaller than : %22.15E\n",accuracy);
    /* create starting point wavefunction */
    CreateRealRandomWaveFunction(&psiGrd[i], NFermion, NElectron);
    /* Find groundstate */
    LanczosGroundState(H, &(psiGrd[i]), &accuracy, 0.000000000000001, 0, &(psiGrd[i]), 0);
    PsiHermitianOperatorPsi(H, psiGrd[i], &EGrd[i]);
    PsiHermitianOperatorPsi(Coulomb, psiGrd[i], &ECoulomb[i]);
    PsiHermitianOperatorPsi(OppH1, psiGrd[i], &EKinetic[i]);
    PsiHermitianOperatorSqrPsi(H, psiGrd[i], &VarEGrd[i]);
    VarEGrd[i]-=Sqr(EGrd[i]);
    printf("\n E                      VarE                   ECoulomb               Ekinetic                   Nstates\n");
    printf("%22.15E %22.15E %22.15E %22.15E %10i\n",EGrd[i],VarEGrd[i],ECoulomb[i],EKinetic[i],psiGrd[i].N);
  }
  /* print an overview of groundstates */
  printf("\n E                      VarE                   ECoulomb               Ekinetic                   Nstates       Accuracy\n");
  for(i=0;i<12;i++)
  {
    accuracy=1/pow(10, (double)(i+2));
    printf("%22.15E %22.15E %22.15E %22.15E %10i %22.15E\n",EGrd[i],VarEGrd[i],ECoulomb[i],EKinetic[i],psiGrd[i].N,accuracy);
  }
  /* Density matrix */
  CalculateDensityMatrix(&Rho, psiGrd[7]);
  PrintMatrix(Rho);
  /* create rotation of basis */
  for(i=0;i<18;i++)
  {
    for(j=0;j<6;j++)
    {
      Rho.VRe[i][j]=0;
      Rho.VRe[j][i]=0;
    }
  }
  for(i=0;i<6;i++)
  {
    Rho.VRe[i][i]=-1;
  }
  Eigensystem(Rho, &TRho, ValRho);
  printf("Eigenvalues of density matrix (occupations in optimized basis)\n");
  for(i=0;i<18;i++)
  {
    printf("%3i %22.15E\n",i,ValRho[i]);
  }
  MatrixTransform(&H1New, H1, TRho);
  InitOperatorCopySettings(&OppH1New, Sz);
  strcpy(OppH1New.Name, "H_{One electron new}");
  for(i=0;i<NFermion;i++)
  {
    for(j=0;j<NFermion;j++)
    {
      CreateOperatorNumber(NFermion, i, j, &temp);
      OperatorAdd(&OppH1New,temp,H1New.VRe[i][j],0);
      FreeOperator(&temp);
    }
  }
  OperatorChop(&OppH1New, DBL_EPSILON);
  OperatorCopy(&HNew, Coulomb);
  OperatorAdd(&HNew, Sz, 0.1, 0);
  OperatorAdd(&HNew, Sx, 0.00001, 0);
  OperatorAdd(&HNew, OppH1New, 1, 0);
  strcpy(HNew.Name, "Hamiltonian rotated basis");
  /* print the new operators */
  PrintOperatorReadable(OppH1New);
  PrintOperatorReadable(HNew);
  /* loop over accuracies */
  for(i=0;i<12;i++)
  {
    accuracy=1/pow(10, (double)(i+2));
    printf("\n\n\nCalculating the ground-state with truncation of basis states smaller than : %22.15E\n",accuracy);
    /* create starting point wavefunction */
    CreateRealRandomWaveFunction(&psiGrdNew[i], NFermion, NElectron);
    /* Find groundstate */
    LanczosGroundState(HNew, &(psiGrdNew[i]), &accuracy, 0.000000000000001, 0, &(psiGrdNew[i]), 0);
    PsiHermitianOperatorPsi(HNew, psiGrdNew[i], &EGrdNew[i]);
    PsiHermitianOperatorPsi(Coulomb, psiGrdNew[i], &ECoulombNew[i]);
    PsiHermitianOperatorPsi(OppH1New, psiGrdNew[i], &EKineticNew[i]);
    PsiHermitianOperatorSqrPsi(HNew, psiGrdNew[i], &VarEGrdNew[i]);
    VarEGrdNew[i]-=Sqr(EGrdNew[i]);
    printf("\n E                      VarE                   ECoulomb               Ekinetic                   Nstates\n");
    printf("%22.15E %22.15E %22.15E %22.15E %10i\n",EGrdNew[i],VarEGrdNew[i],ECoulombNew[i],EKineticNew[i],psiGrdNew[i].N);
  }
  /* print an overview of groundstates */
  printf("\n E                      VarE                   VarE_old               Enew-Eold                  Nstates      Nstatesold   ccuracy\n");
  for(i=0;i<12;i++)
  {
    accuracy=1/pow(10, (double)(i+2));
    printf("%22.15E %22.15E %22.15E %22.15E %10i %10i %22.15E\n",EGrdNew[i],VarEGrdNew[i],VarEGrd[i],EGrdNew[i]-EGrd[i],psiGrdNew[i].N,psiGrd[i].N,accuracy);
  }
  
  
  /* remove second shell from H */
  for(i=0;i<12;i++)
  {
    WaveFunctionFree(&psiGrdNew[i]);
  }
  FreeOperator(&OppH1New);
  FreeOperator(&HNew);
  for(i=12;i<18;i++)
  {
    for(j=0;j<12;j++)
    {
      H1New.VRe[i][j]=0;
      H1New.VRe[j][i]=0;
    }
  }
  for(i=12;i<18;i++)
  {
    H1New.VRe[i][i]-=2;
  }
  InitOperatorCopySettings(&OppH1New, Sz);
  strcpy(OppH1New.Name, "H_{One electron new}");
  for(i=0;i<NFermion;i++)
  {
    for(j=0;j<NFermion;j++)
    {
      CreateOperatorNumber(NFermion, i, j, &temp);
      OperatorAdd(&OppH1New,temp,H1New.VRe[i][j],0);
      FreeOperator(&temp);
    }
  }
  OperatorChop(&OppH1New, DBL_EPSILON);
  OperatorCopy(&HNew, Coulomb);
  OperatorAdd(&HNew, Sz, 0.1, 0);
  OperatorAdd(&HNew, Sx, 0.00001, 0);
  OperatorAdd(&HNew, OppH1New, 1, 0);
  strcpy(HNew.Name, "Hamiltonian rotated basis");
  /* print the new operators */
  PrintOperatorReadable(OppH1New);
  PrintOperatorReadable(HNew);
  /* loop over accuracies */
  for(i=0;i<12;i++)
  {
    accuracy=1/pow(10, (double)(i+2));
    printf("\n\n\nCalculating the ground-state with truncation of basis states smaller than : %22.15E\n",accuracy);
    /* create starting point wavefunction */
    CreateRealRandomWaveFunction(&psiGrdNew[i], NFermion, NElectron);
    /* Find groundstate */
    LanczosGroundState(HNew, &(psiGrdNew[i]), &accuracy, 0.000000000000001, 0, &(psiGrdNew[i]), 0);
    PsiHermitianOperatorPsi(HNew, psiGrdNew[i], &EGrdNew[i]);
    PsiHermitianOperatorPsi(Coulomb, psiGrdNew[i], &ECoulombNew[i]);
    PsiHermitianOperatorPsi(OppH1New, psiGrdNew[i], &EKineticNew[i]);
    PsiHermitianOperatorSqrPsi(HNew, psiGrdNew[i], &VarEGrdNew[i]);
    VarEGrdNew[i]-=Sqr(EGrdNew[i]);
    EGrdNew[i]+=12;
    printf("\n E                      VarE                   ECoulomb               Ekinetic                   Nstates\n");
    printf("%22.15E %22.15E %22.15E %22.15E %10i\n",EGrdNew[i],VarEGrdNew[i],ECoulombNew[i],EKineticNew[i],psiGrdNew[i].N);
  }
  /* print an overview of groundstates */
  printf("\n E                      VarE                   VarE_old               Enew-Eold                  Nstates      Nstatesold   ccuracy\n");
  for(i=0;i<12;i++)
  {
    accuracy=1/pow(10, (double)(i+2));
    printf("%22.15E %22.15E %22.15E %22.15E %10i %10i %22.15E\n",EGrdNew[i],VarEGrdNew[i],VarEGrd[i],EGrdNew[i]-EGrd[i],psiGrdNew[i].N,psiGrd[i].N,accuracy);
  }
}