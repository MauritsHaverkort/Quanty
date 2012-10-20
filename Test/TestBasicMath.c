void TestMatrixTypeConversion()
{
  MatrixType A,U,UT,MI,ATriFull,AComplex;
  TriDiagonalMatrixType ATri;
  AndersonMatrixType AAnd;
  BlokTriDiagonalMatrixType ATriBlok;
  BlokAndersonMatrixType AAndBlok;
  double val[100];
  uint32_t i,j;
  printf(" Test the conversion of a matrix to different representations \n");
  printf(" Create a matrix A \n");
  A.Ni=12;
  A.Nj=12;
  A.QComplex=0;
  strcpy(A.Name,"A");
  InitMatrix(&A);
  for(i=0;i<A.Ni;i++)
  {
    for(j=0;j<A.Nj;j++)
    {
      A.VRe[i][j]=-30.0+i+j+(i*j)/(i-j+0.5)+(i*j)/(i+j+0.5);
      A.VRe[j][i]=A.VRe[i][j];
    }
  }
  PrintMatrix(A);
  AComplex.Ni=12;
  AComplex.Nj=12;
  AComplex.QComplex=1;
  strcpy(AComplex.Name,"AComplex");
  InitMatrix(&AComplex);
  for(i=0;i<AComplex.Ni;i++)
  {
    for(j=0;j<AComplex.Nj;j++)
    {
      AComplex.VRe[i][j]=-30.0+i+j+(i*j)/(i-j+0.5)+(i*j)/(i+j+0.5);
      AComplex.VRe[j][i]=AComplex.VRe[i][j];
      AComplex.VIm[i][j]=((double)i-(double)j);
    }
  }
  PrintMatrix(AComplex);
  
  printf("\n\n\n test MatrixToTriDiagonalMatrix\n\n\n");
  MatrixToTriDiagonalMatrix(A, &U, &ATri);
  
  MatrixCopy(&UT,U);
  MatrixTranspose(&UT);
  MatrixMultiply(&MI,UT,U);
  MatrixTransform(&ATriFull,A,U);
  MatrixChop(&MI,DBL_EPSILON*100);
  MatrixChop(&ATriFull,DBL_EPSILON*100);
  PrintMatrix(U);
  PrintMatrix(MI);
  PrintMatrix(ATriFull);
  PrintTriDiagonalMatrix(ATri);
  FreeMatrix(&U);
  FreeMatrix(&UT);
  FreeMatrix(&MI);
  FreeMatrix(&ATriFull);
  FreeTriDiagonalMatrix(&ATri);
  
  printf("\n\n\n test MatrixToTriDiagonalMatrix Complex\n\n\n");
  MatrixToTriDiagonalMatrix(AComplex, &U, &ATri);
  
  MatrixCopy(&UT,U);
  MatrixConjugateTranspose(&UT);
  MatrixMultiply(&MI,UT,U);
  MatrixTransform(&ATriFull,AComplex,U);
  MatrixChop(&MI,DBL_EPSILON*100);
  MatrixChop(&ATriFull,DBL_EPSILON*100);
  PrintMatrix(U);
  PrintMatrix(MI);
  PrintMatrix(ATriFull);
  PrintTriDiagonalMatrix(ATri);
  FreeMatrix(&U);
  FreeMatrix(&UT);
  FreeMatrix(&MI);
  FreeMatrix(&ATriFull);
  FreeTriDiagonalMatrix(&ATri);
    
  printf("\n\n\n test MatrixToEigensystem\n\n\n");
  MatrixToEigensystem(A, &U, val);
  
  MatrixCopy(&UT,U);
  MatrixTranspose(&UT);
  MatrixMultiply(&MI,UT,U);
  MatrixTransform(&ATriFull,A,U);
  MatrixChop(&MI,DBL_EPSILON*100);
  MatrixChop(&ATriFull,DBL_EPSILON*100);
  PrintMatrix(U);
  PrintMatrix(MI);
  PrintMatrix(ATriFull);
  printf("Eigenvalues\n");
  for(i=0;i<A.Ni;i++)
  {
    printf("%22.12f \n",val[i]);
  }
  FreeMatrix(&U);
  FreeMatrix(&UT);
  FreeMatrix(&MI);
  FreeMatrix(&ATriFull);
    
  printf("\n\n\n test MatrixToEigensystem Complex\n\n\n");
  MatrixToEigensystem(AComplex, &U, val);
  
  MatrixCopy(&UT,U);
  MatrixConjugateTranspose(&UT);
  MatrixMultiply(&MI,UT,U);
  MatrixTransform(&ATriFull,AComplex,U);
  MatrixChop(&MI,DBL_EPSILON*100);
  MatrixChop(&ATriFull,DBL_EPSILON*100);
  PrintMatrix(U);
  PrintMatrix(MI);
  PrintMatrix(ATriFull);
  printf("Eigenvalues\n");
  for(i=0;i<A.Ni;i++)
  {
    printf("%22.12f \n",val[i]);
  }
  FreeMatrix(&U);
  FreeMatrix(&UT);
  FreeMatrix(&MI);
  FreeMatrix(&ATriFull);
  
  printf("\n\n\n test MatrixToAndersonMatrix\n\n\n");
  MatrixToAndersonMatrix(A, &U, &AAnd);
  
  MatrixCopy(&UT,U);
  MatrixTranspose(&UT);
  MatrixMultiply(&MI,UT,U);
  MatrixTransform(&ATriFull,A,U);
  MatrixChop(&MI,DBL_EPSILON*100);
  MatrixChop(&ATriFull,DBL_EPSILON*100);
  PrintMatrix(U);
  PrintMatrix(MI);
  PrintMatrix(ATriFull);
  PrintAndersonMatrix(AAnd);
  FreeMatrix(&U);
  FreeMatrix(&UT);
  FreeMatrix(&MI);
  FreeMatrix(&ATriFull);
  FreeAndersonMatrix(&AAnd);
  
  printf("\n\n\n test MatrixToAndersonMatrix Complex\n\n\n");
  MatrixToAndersonMatrix(AComplex, &U, &AAnd);
  
  MatrixCopy(&UT,U);
  MatrixConjugateTranspose(&UT);
  MatrixMultiply(&MI,UT,U);
  MatrixTransform(&ATriFull,AComplex,U);
  MatrixChop(&MI,DBL_EPSILON*100);
  MatrixChop(&ATriFull,DBL_EPSILON*100);
  PrintMatrix(U);
  PrintMatrix(MI);
  PrintMatrix(ATriFull);
  PrintAndersonMatrix(AAnd);
  FreeMatrix(&U);
  FreeMatrix(&UT);
  FreeMatrix(&MI);
  FreeMatrix(&ATriFull);
  FreeAndersonMatrix(&AAnd);
  
  printf("\n\n\n test MatrixToBlokTridiagonalMatrix\n\n\n");
  MatrixToBlokTriDiagonalMatrix(A, &U, &ATriBlok,3);
  
  MatrixCopy(&UT,U);
  MatrixTranspose(&UT);
  MatrixMultiply(&MI,UT,U);
  MatrixTransform(&ATriFull,A,U);
  MatrixChop(&MI,DBL_EPSILON*100);
  MatrixChop(&ATriFull,DBL_EPSILON*100);
  PrintMatrix(U);
  PrintMatrix(MI);
  PrintMatrix(ATriFull);
  PrintBlokTriDiagonalMatrix(ATriBlok);
  FreeMatrix(&U);
  FreeMatrix(&UT);
  FreeMatrix(&MI);
  FreeMatrix(&ATriFull);
  FreeBlokTriDiagonalMatrix(&ATriBlok);
  
  printf("\n\n\n test MatrixToBlokTridiagonalMatrix Complex\n\n\n");
  MatrixToBlokTriDiagonalMatrix(AComplex, &U, &ATriBlok,3);
  
  MatrixCopy(&UT,U);
  MatrixConjugateTranspose(&UT);
  MatrixMultiply(&MI,UT,U);
  MatrixTransform(&ATriFull,AComplex,U);
  MatrixChop(&MI,DBL_EPSILON*1000);
  MatrixChop(&ATriFull,DBL_EPSILON*1000);
  PrintMatrix(U);
  PrintMatrix(MI);
  PrintMatrix(ATriFull);
  PrintBlokTriDiagonalMatrix(ATriBlok);
  FreeMatrix(&U);
  FreeMatrix(&UT);
  FreeMatrix(&MI);
  FreeMatrix(&ATriFull);
  FreeBlokTriDiagonalMatrix(&ATriBlok);
  
  
  
  printf("\n\n\n test MatrixToBlokAndersonMatrix\n\n\n");
  MatrixToBlokAndersonMatrix(A, &U, &AAndBlok,3);
  
  MatrixCopy(&UT,U);
  MatrixTranspose(&UT);
  MatrixMultiply(&MI,UT,U);
  MatrixTransform(&ATriFull,A,U);
  MatrixChop(&MI,DBL_EPSILON*1000);
  MatrixChop(&ATriFull,DBL_EPSILON*1000);
  PrintMatrix(U);
  PrintMatrix(MI);
  PrintMatrix(ATriFull);
  PrintBlokAndersonMatrix(AAndBlok);
  FreeMatrix(&U);
  FreeMatrix(&UT);
  FreeMatrix(&MI);
  FreeMatrix(&ATriFull);
  FreeBlokAndersonMatrix(&AAndBlok);
  
  printf("\n\n\n test MatrixToBlokAndersonMatrix Complex\n\n\n");
  MatrixToBlokAndersonMatrix(AComplex, &U, &AAndBlok,3);
  
  MatrixCopy(&UT,U);
  MatrixConjugateTranspose(&UT);
  MatrixMultiply(&MI,UT,U);
  MatrixTransform(&ATriFull,AComplex,U);
  MatrixChop(&MI,DBL_EPSILON*1000);
  MatrixChop(&ATriFull,DBL_EPSILON*1000);
  PrintMatrix(U);
  PrintMatrix(MI);
  PrintMatrix(ATriFull);
  PrintBlokAndersonMatrix(AAndBlok);
  FreeMatrix(&U);
  FreeMatrix(&UT);
  FreeMatrix(&MI);
  FreeMatrix(&ATriFull);
  FreeBlokAndersonMatrix(&AAndBlok);

}

void TestMatrixVector()
{
  MatrixType A,AComplex,IdentM,Acopy,AA,AB,AC,AD,ADCT,IAD,TA,TAATA;
  double AvecRe[12],AvecIm[12],BvecRe[12],BvecIm[12];
  uint32_t i,j;
  printf("Testing MatrixVector functions\n");
  printf(" Create a matrix A \n");
  A.Ni=12;
  A.Nj=12;
  A.QComplex=0;
  InitMatrix(&A);
  strcpy(A.Name,"A");
  for(i=0;i<A.Ni;i++)
  {
    for(j=0;j<A.Nj;j++)
    {
      A.VRe[i][j]=-30.0+i+j+(i*j)/(i-j+0.5)+(i*j)/(i+j+0.5);
      A.VRe[j][i]=A.VRe[i][j];
    }
  }
  PrintMatrix(A);
  AComplex.Ni=12;
  AComplex.Nj=12;
  AComplex.QComplex=1;
  InitMatrix(&AComplex);
  strcpy(AComplex.Name,"AComplex");
  for(i=0;i<AComplex.Ni;i++)
  {
    for(j=0;j<AComplex.Nj;j++)
    {
      AComplex.VRe[i][j]=-30.0+i+j+(i*j)/(i-j+0.5)+(i*j)/(i+j+0.5);
      AComplex.VRe[j][i]=AComplex.VRe[i][j];
      AComplex.VIm[i][j]=((double)i-(double)j);
    }
  }
  PrintMatrix(AComplex);
  IdentM.Ni=12;
  IdentM.Nj=12;
  IdentM.QComplex=0;
  InitIdentityMatrix(&IdentM);
  PrintMatrix(IdentM);
  printf("make A complex\n");
  MakeRealMatrixComplex(&A);
  PrintMatrix(A);
  printf("Chop A with epsilon=10\n");
  MatrixChop(&A,10.0);
  PrintMatrix(A);
  printf("Transpose AComplex\n");
  MatrixTranspose(&AComplex);
  PrintMatrix(AComplex);
  printf("ComplexConjugate AComplex\n");
  MatrixComplexConjugate(&AComplex);
  PrintMatrix(AComplex);
  printf("Conjugate Transpose AComplex\n");
  MatrixConjugateTranspose(&AComplex);
  PrintMatrix(AComplex);
  printf("Make a copy of AComplex\n");
  MatrixCopy(&Acopy,AComplex);
  PrintMatrix(Acopy);
  printf("Invert AComplex\n");
  MatrixInverse(&AComplex);
  PrintMatrix(AComplex);
  printf("Multiply Acopy with AComplex^-1\n");
  MatrixMultiply(&AA,Acopy,AComplex);
  PrintMatrix(AA);
  printf("Add I * AComplex to A\n");
  MatrixAdd(&A,AComplex,0,1);
  PrintMatrix(A);
  printf("Multiply A by 3\n");
  MatrixScalarMultiply(&A,3,0);
  PrintMatrix(A);
  printf("Matrix Transform AB=A^* AA A^T\n");
  MatrixTransform(&AB,AA,A);
  PrintMatrix(AB);
  printf("Matrix Transform part AC=(A^* AB A^T)[2;;5,3;;7]\n");
  MatrixTransformPart(&AC,AA,A,2,5,3,7);
  PrintMatrix(AC);
  printf("AD=AB[2;;5,3;;7]\n");
  MatrixPart(&AD,AB,2,5,3,7);
  PrintMatrix(AD);

  for(i=0;i<12;i++)
  {
    AvecRe[i]=i;
    AvecIm[i]=0;
    BvecRe[i]=i;
    BvecIm[i]=-1.0*i/10.0;
  }
  
  printf("Orthonormalize AD\n");
  MatrixOrthonormalize(&AD);
  printf("Matrix Copy\n");
  MatrixCopy(&ADCT,AD);
  printf("MatrixConjugateTranspose\n");
  MatrixConjugateTranspose(&ADCT);
  printf("MatrixMultiply\n");
  MatrixMultiply(&IAD,AD,ADCT);
  PrintMatrix(AD);
  PrintMatrix(IAD);
  
  FreeMatrix(&IAD);
  FreeMatrix(&ADCT);
  MatrixConjugateTranspose(&AD);
  printf("Orthonormalize AD\n");
  MatrixOrthonormalize(&AD);
  printf("Matrix Copy\n");
  MatrixCopy(&ADCT,AD);
  printf("MatrixConjugateTranspose\n");
  MatrixConjugateTranspose(&ADCT);
  printf("MatrixMultiply\n");
  MatrixMultiply(&IAD,AD,ADCT);
  PrintMatrix(AD);
  PrintMatrix(IAD);
  
  printf(" Create a matrix TA \n");
  TA.Ni=5;
  TA.Nj=12;
  TA.QComplex=0;
  InitMatrix(&TA);
  strcpy(TA.Name,"TA");
  for(i=0;i<TA.Ni;i++)
  {
    for(j=0;j<TA.Nj;j++)
    {
      if(i==j)
      {
        TA.VRe[i][j]=1;
      }
      else
      {
        TA.VRe[i][j]=0;
      }
    }
  }
  PrintMatrix(TA);
  PrintMatrix(A);
  MatrixTransform(&TAATA,A,TA);
  PrintMatrix(TAATA);

  
}


void TestAllMatrix()
{
 printf("Testing the matrix library and the different matrix representations\n");
 MatrixType M,A,B,BT,MI,IdentM,ATri,T;
 TriDiagonalMatrixType MTri;
 BlokTriDiagonalMatrixType MBlokTri;
 uint32_t i,j;
 double singular,val[20];
 A.Ni=6;
 A.Nj=4;
 A.QComplex=0;
 strcpy(A.Name,"A");
 InitMatrix(&A);
 for(i=0;i<A.Ni;i++)
 {
   for(j=0;j<A.Nj;j++)
   {
     A.VRe[i][j]=100+10*i+j;
   }
 }
 
 printf("\nCreated a Real %i by %i matrix with elements value 100+10*i+j\n",A.Ni,A.Nj);
 PrintMatrix(A);
 

 printf("\nAdd A to B with B equal to 0.1*j+0.01*i\n");
 B.Ni=6;
 B.Nj=4;
 B.QComplex=0;
 strcpy(B.Name,"B");
 InitMatrix(&B);
 for(i=0;i<B.Ni;i++)
 {
   for(j=0;j<B.Nj;j++)
   {
     B.VRe[i][j]=0.1*j+0.01*i;
   }
 }
 MatrixAdd(&B,A,1,0);
 PrintMatrix(B);
 printf("\nTranspose the matrix\n");
 MatrixTranspose(&A);
 PrintMatrix(A);
 printf("\nM = A B\n");
 strcpy(M.Name,"M");
 MatrixMultiply(&M,A,B);
 PrintMatrix(M);
 
 M.VRe[3][3]+=0.00100;
 M.VRe[2][2]+=0.00200;
 M.VRe[1][1]+=0.00300;
 M.VRe[0][0]+=0.00400;
 PrintMatrix(M);
 MatrixCopy(&MI,M);
 strcpy(MI.Name,"MI");
 singular=MatrixInverse(&MI);
 printf("Matrix inverted with signular=%22.15E\n",singular);
 PrintMatrix(MI);
 printf("\nIdentM = M MI\n");
 strcpy(IdentM.Name,"IdentM");
 MatrixMultiply(&IdentM,M,MI);
 PrintMatrix(IdentM);
 FreeMatrix(&IdentM);
 printf("\nIdentM = MI M\n");
 MatrixMultiply(&IdentM,MI,M);
 PrintMatrix(IdentM);
 FreeMatrix(&IdentM);
 FreeMatrix(&MI);
 MakeRealMatrixComplex(&M);
 printf("\nSame but complex \n\n");
 M.VIm[0][1]=10.5;
 M.VIm[2][0]=10.5;
 MatrixCopy(&MI,M);
 singular=MatrixInverse(&MI);
 printf("Matrix inverted with signular=%22.15E\n",singular);
 MatrixMultiply(&IdentM,MI,M);
 PrintMatrix(M);
 PrintMatrix(MI);
 PrintMatrix(IdentM);
 singular=MatrixInverse(&MI);
 printf("Matrix inverted with signular=%22.15E\n",singular);
 MatrixAdd(&MI,M,-1,0);
 printf("MI=((M)^{-1})^{-1})-M\n");
 PrintMatrix(MI);
 
 
 printf("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n");
 FreeMatrix(&MI);
 FreeMatrix(&M);
 FreeMatrix(&A);
 FreeMatrix(&B);
 A.Ni=6;
 A.Nj=6;
 A.QComplex=0;
 strcpy(A.Name,"A");
 InitMatrix(&A);
 for(i=0;i<A.Ni;i++)
 {
   for(j=0;j<A.Nj;j++)
   {
     A.VRe[i][j]=1*i+1*j+0.000001*abs(i-j);
     if(i==j) A.VRe[i][j]+=10+0 * (i+1);
   }
 }
 strcpy(B.Name,"B=Tri Diagonalizer of A");
 strcpy(BT.Name,"BT=Transpose of Tri Diagonalizer of A");
 strcpy(MI.Name,"BT . B");
 strcpy(MTri.Name,"MTri=Tri Diagonal form of A");
 strcpy(ATri.Name,"ATri=Tri Diagonal form of A by applying B . A . B^T");
 TriDiagonalize(A, &B, &MTri);
 MatrixCopy(&BT,B);
 MatrixTranspose(&BT);
 MatrixMultiply(&MI,BT,B);
 MatrixTransform(&ATri,A,B);
 PrintMatrix(A);
 PrintMatrix(B);
 PrintMatrix(MI);
 PrintMatrix(ATri);
 PrintTriDiagonalMatrix(MTri);
 
 
 
 printf("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n");
 printf("free MI\n");
 FreeMatrix(&MI);
 printf("free ATri\n");
 FreeMatrix(&ATri);
 printf("free A\n");
 FreeMatrix(&A);
 printf("free B\n");
 FreeMatrix(&B);
 printf("free MTri\n");
 FreeTriDiagonalMatrix(&MTri);
 printf("finished free\n");
 A.Ni=6;
 A.Nj=6;
 A.QComplex=1;
 strcpy(A.Name,"A");
 InitMatrix(&A);
 for(i=0;i<A.Ni;i++)
 {
   for(j=0;j<A.Nj;j++)
   {
     A.VRe[i][j]=1*i+1*j+0.000001*abs(i-j);
     A.VIm[i][j]=1*i+1*j+0.000001*abs(i-j);
     if(i==j)A.VIm[i][j]=0;
     if(i<j)A.VIm[i][j]*=-1;
     if(i==j) A.VRe[i][j]+=10+0 * (i+1);
   }
 }
 strcpy(B.Name,"B=Tri Diagonalizer of A");
 strcpy(BT.Name,"BT=Conugate Transpose of Tri Diagonalizer of A");
 strcpy(MI.Name,"B^+ . B");
 strcpy(MTri.Name,"MTri=Tri Diagonal form of A");
 strcpy(ATri.Name,"ATri=Tri Diagonal form of A by applying B . A . B^+");
 TriDiagonalize(A, &B, &MTri);
 MatrixCopy(&BT,B);
 printf("BT<-B\n");
 PrintMatrix(BT);
 MatrixConjugateTranspose(&BT);
 printf("BT<-BT^+\n");
 PrintMatrix(BT);
 MatrixMultiply(&MI,BT,B);
 PrintMatrix(MI);
 
 printf("\n\n\n\n Tri diagonalization of a complex matrix \n\n\n");
 MatrixTransform(&ATri,A,B);
 PrintMatrix(A);
 PrintMatrix(B);
 PrintMatrix(ATri);
 PrintTriDiagonalMatrix(MTri);
 TriDiagonalMatrixToEigensystem(MTri, &T, val);
 PrintMatrix(T);
 for(i=0;i<6;i++)
 {
   printf("%22.15E \n",val[i]);
 }
 
 printf("\n\n\n\n\n  Start Eigensystem of a matrix  \n\n\n\n\n");
 FreeMatrix(&T);
 FreeMatrix(&A);
 A.Ni=12;
 A.Nj=12;
 A.QComplex=0;
 strcpy(A.Name,"A");
 InitMatrix(&A);
 for(i=0;i<A.Ni;i++)
 {
   for(j=0;j<A.Nj;j++)
   {
     A.VRe[i][j]=i+j+(i*j)/(i-j+0.5)+(i*j)/(i+j+0.5);
     A.VRe[j][i]=A.VRe[i][j];
   }
 }
 

 Eigensystem(A, &T, val);
 PrintMatrix(A);
 PrintMatrix(T);
 FreeMatrix(&T);
 

 
   
 printf("\n\n\n  Start Blok Tri-Diagonalize a matrix \n\n\n");
 
 BlokTriDiagonalize(A, &T, &MBlokTri, 4);

 printf("End of Test All Matrix \n");

}








void TestInterpolate()
{
 double x[11],y[11],ac,bc,cc,a,b;
 int i,N=10;
 for(i=0;i<=10;i++)
 {
   x[i]=0.1*i*i-0.235*i+0.5;
   y[i]=3*Sqr(x[i])+2*x[i]-5;
 }
 LinearInterpolate(x, y, N, &a, &b);
 QuadraticInterpolate(x, y, N, &ac, &bc, &cc);
 printf("%22.15E %22.15E\n",a,b);
 printf("%22.15E %22.15E %22.15E\n",ac,bc,cc);
}

void TestMatrixInverse()
{
  MatrixType M;
  M.QComplex=1;
  M.Ni=4;
  M.Nj=4;
  strcpy(M.Name,"M");
  double **MRe;
  double **MIm;
  int i,j,N;
  double singular;
  N=4;
  MRe=MatrixCalloc(N,N);
  MIm=MatrixCalloc(N,N);
  MRe[0][2]=2.0;
  MRe[1][1]=3.0;
  MRe[2][0]=4.0;
  MRe[3][3]=2.0;
  MIm[3][2]=2.0;
  
  
  MRe[0][0]=2.0;
  MRe[1][1]=3.0;
  MRe[2][3]=4.0;
  MRe[3][2]=2.0;
  MIm[3][2]=2.0;
  
  MRe[0][0]=0.5000001;
  MRe[0][1]=0.6;
  MRe[0][2]=0.7;
  MRe[0][3]=0.8;
  MRe[0][4]=0.9;
  MIm[0][2]=5.0;
  
  MRe[1][0]=0.5;
  MRe[1][1]=0.6;
  MRe[1][2]=0.7;
  MRe[1][3]=0.8;
  MRe[1][4]=0.9;
  MIm[1][2]=5.0;
  
  M.VIm=MIm;
  M.VRe=MRe;
  
  for(i=0;i<N;i++)
  {
    for(j=0;j<N;j++)
    {
      printf("(%12.5f %12.5f)",MRe[i][j],MIm[i][j]);
    }
    printf("\n");
  }
  printf("\n");
  singular=MatrixInverse(&M);
  printf("Singular=%22.15E\n",singular);
  for(i=0;i<N;i++)
  {
    for(j=0;j<N;j++)
    {
      printf("(%12.5E %12.5E)",MRe[i][j],MIm[i][j]);
    }
    printf("\n");
  }
  printf("\n");
  singular=MatrixInverse(&M);
  printf("Singular=%22.15E\n",singular);
  for(i=0;i<N;i++)
  {
    for(j=0;j<N;j++)
    {
      printf("(%12.5f %12.5f)",MRe[i][j],MIm[i][j]);
    }
    printf("\n");
  }
}

void TestBasicMath()
{
  double a,b;
  int l1,l2,l3,m1,m2,m3;
  for(l1=0;l1<50;l1++)
  {
    a=Factorial(l1);
    printf("%4i!= %22.15E\n",l1,a);
  }
  l1=0;l2=0;l3=0;m1=0;m2=0;m3=0;
  a=ThreeJSymbol(l1,m1,l2,m2,l3,m3);
  b=1.0;
  printf("(%4i %4i %4i)\n(%4i %4i %4i) =%22.15E (%22.15E)\n",l1,l2,l3,m1,m2,m3,a,b);
  
  l1=1;l2=0;l3=1;m1=0;m2=0;m3=0;
  a=ThreeJSymbol(l1,m1,l2,m2,l3,m3);
  b=-sqrt(1.0/3.0);
  printf("(%4i %4i %4i)\n(%4i %4i %4i) =%22.15E (%22.15E)\n",l1,l2,l3,m1,m2,m3,a,b);
  
  l1=2;l2=0;l3=2;m1=0;m2=0;m3=0;
  a=ThreeJSymbol(l1,m1,l2,m2,l3,m3);
  b=sqrt(1.0/5.0);
  printf("(%4i %4i %4i)\n(%4i %4i %4i) =%22.15E (%22.15E)\n",l1,l2,l3,m1,m2,m3,a,b);
  
  l1=2;l2=2;l3=2;m1=0;m2=0;m3=0;
  a=ThreeJSymbol(l1,m1,l2,m2,l3,m3);
  b=-sqrt(2.0/35.0);
  printf("(%4i %4i %4i)\n(%4i %4i %4i) =%22.15E (%22.15E)\n",l1,l2,l3,m1,m2,m3,a,b);
  
  l1=2;l2=2;l3=4;m1=1;m2=2;m3=-3;
  a=ThreeJSymbol(l1,m1,l2,m2,l3,m3);
  b=-sqrt(1.0/18.0);
  printf("(%4i %4i %4i)\n(%4i %4i %4i) =%22.15E (%22.15E)\n",l1,l2,l3,m1,m2,m3,a,b);
  
  l1=7;l2=5;l3=4;m1=1;m2=2;m3=-3;
  a=ThreeJSymbol(l1,m1,l2,m2,l3,m3);
  b=7*sqrt(5.0/14586.0);
  printf("(%4i %4i %4i)\n(%4i %4i %4i) =%22.15E (%22.15E)\n",l1,l2,l3,m1,m2,m3,a,b);
  
  l1=14;l2=22;l3=9;m1=13;m2=-5;m3=-8;
  a=ThreeJSymbol(l1,m1,l2,m2,l3,m3);
  b=-(229.0/6.0)*sqrt(1.0/9866335845905.0);
  printf("(%4i %4i %4i)\n(%4i %4i %4i) =%22.15E (%22.15E)\n",l1,l2,l3,m1,m2,m3,a,b);
}

void TestDGESVD3()
{
  const char * JOBU="A";
  const char * JOBVT="A";
  int M=4;
  int N=3;
  double A[N][M];
  int LDA=M;
  double S[M];
  double U[M][M];
  int LDU=M;
  double VT[N][N];
  int LDVT=N;
  int LWORK=N*M*10;
  double WORK[LWORK];
  int INFO;
  
  int i,j;
  for(i=0;i<N;i++)
  {
    for(j=0;j<M;j++)
    {
      A[i][j]=5*i+j;
    }
  }
  printf("A=\n");
  for(i=0;i<N;i++)
  {
    for(j=0;j<M;j++)
    {
      printf("%12.5e ",A[i][j]);;
    }
    printf("\n");
  }
  DGESVD( JOBU, JOBVT, &M, &N, A, &LDA, S, U, &LDU, VT, &LDVT, WORK, &LWORK, &INFO );
  printf("U=\n");
  for(i=0;i<M;i++)
  {
    for(j=0;j<M;j++)
    {
      printf("%12.5e ",U[i][j]);;
    }
    printf("\n");
  }
  printf("VT=\n");
  for(i=0;i<N;i++)
  {
    for(j=0;j<N;j++)
    {
      printf("%12.5e ",VT[i][j]);;
    }
    printf("\n");
  }
  printf("S=\n");
  for(i=0;i<Max(N,M);i++)
  {
    printf("%12.5e ",S[i]);
  }
  printf("\n");
}


void TestDGESVD2()
{
  const char * JOBU="A";
  const char * JOBVT="A";
  int M=3;
  int N=4;
  double A[N][M];
  int LDA=M;
  double S[M];
  double U[M][M];
  int LDU=M;
  double VT[N][N];
  int LDVT=N;
  int LWORK=N*M*10;
  double WORK[LWORK];
  int INFO;
  
  int i,j;
  for(i=0;i<N;i++)
  {
    for(j=0;j<M;j++)
    {
      A[i][j]=5*i+j;
    }
  }
  printf("A=\n");
  for(i=0;i<N;i++)
  {
    for(j=0;j<M;j++)
    {
      printf("%12.5e ",A[i][j]);;
    }
    printf("\n");
  }
  DGESVD( JOBU, JOBVT, &M, &N, A, &LDA, S, U, &LDU, VT, &LDVT, WORK, &LWORK, &INFO );
  printf("U=\n");
  for(i=0;i<M;i++)
  {
    for(j=0;j<M;j++)
    {
      printf("%12.5e ",U[i][j]);;
    }
    printf("\n");
  }
  printf("VT=\n");
  for(i=0;i<N;i++)
  {
    for(j=0;j<N;j++)
    {
      printf("%12.5e ",VT[i][j]);;
    }
    printf("\n");
  }
  printf("S=\n");
  for(i=0;i<Max(N,M);i++)
  {
    printf("%12.5e ",S[i]);
  }
  printf("\n");
  TestDGESVD3();
}


void TestDGESVD()
{
  const char * JOBU="A";
  const char * JOBVT="A";
  int M=2;
  int N=5;
  double A[N][M];
  int LDA=M;
  double S[M];
  double U[M][M];
  int LDU=M;
  double VT[N][N];
  int LDVT=N;
  int LWORK=N*M*10;
  double WORK[LWORK];
  int INFO;
  
  int i,j;
  for(i=0;i<N;i++)
  {
    for(j=0;j<M;j++)
    {
      A[i][j]=120+10*i+j;
    }
  }
  printf("A string\n");
  for(i=0;i<M*N;i++)
  {
    printf("%12.5e ",A[0][i]);
  }
  printf("\n");
  printf("A=\n");
  for(i=0;i<N;i++)
  {
    for(j=0;j<M;j++)
    {
      printf("%12.5e ",A[i][j]);
    }
    printf("\n");
  }
  DGESVD( JOBU, JOBVT, &M, &N, A, &LDA, S, U, &LDU, VT, &LDVT, WORK, &LWORK, &INFO );
  printf("U=\n");
  for(i=0;i<M;i++)
  {
    for(j=0;j<M;j++)
    {
      printf("%12.5e ",U[i][j]);;
    }
    printf("\n");
  }
  printf("VT=\n");
  for(i=0;i<N;i++)
  {
    for(j=0;j<N;j++)
    {
      printf("%12.5e ",VT[i][j]);;
    }
    printf("\n");
  }
  printf("S=\n");
  for(i=0;i<Max(N,M);i++)
  {
    printf("%12.5e ",S[i]);
  }
  printf("\n");
  TestDGESVD2();
}

void TestComplex()
{

int n, inca = 1, incb = 1, i,N=5;
Complex16 a[N], b[N], c;
void zdotc();
n = N;
for( i = 0; i < n; i++ ){
a[i].Re = (double)i; a[i].Im = (double)i * 2.0;
b[i].Re = (double)(n - i); b[i].Im = (double)i * 2.0;
}
zdotc( &c, &n, a, &inca, b, &incb );
printf( "The complex dot product is: ( %6.2f, %6.2f) \n", c.Re, c.Im ); 
  
}

void TestTriDiagonalMatrix()
{
  TriDiagonalMatrixType A,B;
  AndersonMatrixType C;
  DoubleTriDiagonalMatrixType D;
  ListOfPoles E;
  MatrixType T,T2;
  uint32_t i;
  A.N=4;
  InitTriDiagonalMatrix(&A);
  A.a[0]=0;
  A.b[0]=1;
  for(i=1;i<A.N;i++)
  {
    A.a[i]=i;
    A.b[i]=i*10;
  }
  A.a[A.N]=A.N;
  
  PrintTriDiagonalMatrix(A);
  TriDiagonalMatrixPart(&B,A,2,A.N);
  PrintTriDiagonalMatrix(B);
  FreeTriDiagonalMatrix(&B);
  printf("\n\n\n\n\n");
  TriDiagonalMatrixToAndersonMatrix(A,&T,&C);
  
  PrintAndersonMatrix(C);
  
  PrintMatrix(T);
  
  
  printf("\n\n\n\n\n");
  
  AndersonMatrixToDoubleTriDiagonalMatrix(C,&T2,&D,0.0,0.0001);
  PrintDoubleTriDiagonalMatrix(D);
  PrintMatrix(T2);  
  
  
  printf("\n\n\n\n\n");
  A.a[0]=3;
  A.b[0]=5;
  for(i=1;i<A.N;i++)
  {
    A.a[i]=i;
    A.b[i]=i*10;
  }
  A.a[A.N]=A.N;
  TriDiagonalMatrixToListOfPoles(A, &E);
  ListOfPolesToTridiagonalMatrix(E,&B,A.N);
  PrintTriDiagonalMatrix(A);
  PrintTriDiagonalMatrix(B);
  PrintListOfPoles(E);
}