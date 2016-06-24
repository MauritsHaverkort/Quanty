void TestCompactMatrixCallers();
void TestBlockLanczos();

void TestMatrixTypeConversion()
{
  MatrixType A, U, UT, MI, ATriFull, AComplex;
  TriDiagonalMatrixType ATri;
  AndersonMatrixType AAnd;
  BlokTriDiagonalMatrixType ATriBlok;
  BlokAndersonMatrixType AAndBlok;
  double val[100];
  uint32_t i, j;
  printf(" Test the conversion of a matrix to different representations \n");
  printf(" Create a matrix A \n");
  A.Ni = 12;
  A.Nj = 12;
  A.QComplex = 0;
  strcpy(A.Name, "A");
  InitMatrix(&A);
  for (i = 0; i < A.Ni; i++)
  {
    for (j = 0; j < A.Nj; j++)
    {
      A.VRe[i][j] = -30.0 + i + j + (i * j) / (i - j + 0.5) + (i * j) / (i + j + 0.5);
      A.VRe[j][i] = A.VRe[i][j];
    }
  }
  PrintMatrix(A);
  AComplex.Ni = 12;
  AComplex.Nj = 12;
  AComplex.QComplex = 1;
  strcpy(AComplex.Name, "AComplex");
  InitMatrix(&AComplex);
  for (i = 0; i < AComplex.Ni; i++)
  {
    for (j = 0; j < AComplex.Nj; j++)
    {
      AComplex.VRe[i][j] = -30.0 + i + j + (i * j) / (i - j + 0.5) + (i * j) / (i + j + 0.5);
      AComplex.VRe[j][i] = AComplex.VRe[i][j];
      AComplex.VIm[i][j] = ((double)i - (double)j);
    }
  }
  PrintMatrix(AComplex);

  printf("\n\n\n test MatrixToTriDiagonalMatrix\n\n\n");
  MatrixToTriDiagonalMatrix(A, &U, &ATri);

  MatrixCopy(&UT, U);
  MatrixTranspose(&UT);
  MatrixMultiply(&MI, UT, U);
  MatrixTransform(&ATriFull, A, U);
  MatrixChop(&MI, DBL_EPSILON * 100);
  MatrixChop(&ATriFull, DBL_EPSILON * 100);
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

  MatrixCopy(&UT, U);
  MatrixConjugateTranspose(&UT);
  MatrixMultiply(&MI, UT, U);
  MatrixTransform(&ATriFull, AComplex, U);
  MatrixChop(&MI, DBL_EPSILON * 100);
  MatrixChop(&ATriFull, DBL_EPSILON * 100);
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

  MatrixCopy(&UT, U);
  MatrixTranspose(&UT);
  MatrixMultiply(&MI, UT, U);
  MatrixTransform(&ATriFull, A, U);
  MatrixChop(&MI, DBL_EPSILON * 100);
  MatrixChop(&ATriFull, DBL_EPSILON * 100);
  PrintMatrix(U);
  PrintMatrix(MI);
  PrintMatrix(ATriFull);
  printf("Eigenvalues\n");
  for (i = 0; i < A.Ni; i++)
  {
    printf("%22.12f \n", val[i]);
  }
  FreeMatrix(&U);
  FreeMatrix(&UT);
  FreeMatrix(&MI);
  FreeMatrix(&ATriFull);

  printf("\n\n\n test MatrixToEigensystem Complex\n\n\n");
  MatrixToEigensystem(AComplex, &U, val);

  MatrixCopy(&UT, U);
  MatrixConjugateTranspose(&UT);
  MatrixMultiply(&MI, UT, U);
  MatrixTransform(&ATriFull, AComplex, U);
  MatrixChop(&MI, DBL_EPSILON * 100);
  MatrixChop(&ATriFull, DBL_EPSILON * 100);
  PrintMatrix(U);
  PrintMatrix(MI);
  PrintMatrix(ATriFull);
  printf("Eigenvalues\n");
  for (i = 0; i < A.Ni; i++)
  {
    printf("%22.12f \n", val[i]);
  }
  FreeMatrix(&U);
  FreeMatrix(&UT);
  FreeMatrix(&MI);
  FreeMatrix(&ATriFull);

  printf("\n\n\n test MatrixToAndersonMatrix\n\n\n");
  MatrixToAndersonMatrix(A, &U, &AAnd);

  MatrixCopy(&UT, U);
  MatrixTranspose(&UT);
  MatrixMultiply(&MI, UT, U);
  MatrixTransform(&ATriFull, A, U);
  MatrixChop(&MI, DBL_EPSILON * 100);
  MatrixChop(&ATriFull, DBL_EPSILON * 100);
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

  MatrixCopy(&UT, U);
  MatrixConjugateTranspose(&UT);
  MatrixMultiply(&MI, UT, U);
  MatrixTransform(&ATriFull, AComplex, U);
  MatrixChop(&MI, DBL_EPSILON * 100);
  MatrixChop(&ATriFull, DBL_EPSILON * 100);
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
  MatrixToBlokTriDiagonalMatrix(A, &U, &ATriBlok, 3);

  MatrixCopy(&UT, U);
  MatrixTranspose(&UT);
  MatrixMultiply(&MI, UT, U);
  MatrixTransform(&ATriFull, A, U);
  MatrixChop(&MI, DBL_EPSILON * 100);
  MatrixChop(&ATriFull, DBL_EPSILON * 100);
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
  MatrixToBlokTriDiagonalMatrix(AComplex, &U, &ATriBlok, 3);

  MatrixCopy(&UT, U);
  MatrixConjugateTranspose(&UT);
  MatrixMultiply(&MI, UT, U);
  MatrixTransform(&ATriFull, AComplex, U);
  MatrixChop(&MI, DBL_EPSILON * 1000);
  MatrixChop(&ATriFull, DBL_EPSILON * 1000);
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
  MatrixToBlokAndersonMatrix(A, &U, &AAndBlok, 3);

  MatrixCopy(&UT, U);
  MatrixTranspose(&UT);
  MatrixMultiply(&MI, UT, U);
  MatrixTransform(&ATriFull, A, U);
  MatrixChop(&MI, DBL_EPSILON * 1000);
  MatrixChop(&ATriFull, DBL_EPSILON * 1000);
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
  MatrixToBlokAndersonMatrix(AComplex, &U, &AAndBlok, 3);

  MatrixCopy(&UT, U);
  MatrixConjugateTranspose(&UT);
  MatrixMultiply(&MI, UT, U);
  MatrixTransform(&ATriFull, AComplex, U);
  MatrixChop(&MI, DBL_EPSILON * 1000);
  MatrixChop(&ATriFull, DBL_EPSILON * 1000);
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
  MatrixType A, AComplex, IdentM, Acopy, AA, AB, AC, AD, ADCT, IAD, TA, TAATA;
  double AvecRe[12], AvecIm[12], BvecRe[12], BvecIm[12];
  uint32_t i, j;
  printf("Testing MatrixVector functions\n");
  printf(" Create a matrix A \n");
  A.Ni = 12;
  A.Nj = 12;
  A.QComplex = 0;
  InitMatrix(&A);
  strcpy(A.Name, "A");
  for (i = 0; i < A.Ni; i++)
  {
    for (j = 0; j < A.Nj; j++)
    {
      A.VRe[i][j] = -30.0 + i + j + (i * j) / (i - j + 0.5) + (i * j) / (i + j + 0.5);
      A.VRe[j][i] = A.VRe[i][j];
    }
  }
  PrintMatrix(A);
  AComplex.Ni = 12;
  AComplex.Nj = 12;
  AComplex.QComplex = 1;
  InitMatrix(&AComplex);
  strcpy(AComplex.Name, "AComplex");
  for (i = 0; i < AComplex.Ni; i++)
  {
    for (j = 0; j < AComplex.Nj; j++)
    {
      AComplex.VRe[i][j] = -30.0 + i + j + (i * j) / (i - j + 0.5) + (i * j) / (i + j + 0.5);
      AComplex.VRe[j][i] = AComplex.VRe[i][j];
      AComplex.VIm[i][j] = ((double)i - (double)j) * sqrt((double)(i+j));
    }
  }
  PrintMatrix(AComplex);
  IdentM.Ni = 12;
  IdentM.Nj = 12;
  IdentM.QComplex = 0;
  InitIdentityMatrix(&IdentM);
  PrintMatrix(IdentM);
  printf("make A complex\n");
  MakeRealMatrixComplex(&A);
  PrintMatrix(A);
  printf("Chop A with epsilon=10\n");
  MatrixChop(&A, 10.0);
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
  MatrixCopy(&Acopy, AComplex);
  PrintMatrix(Acopy);
  printf("Invert AComplex\n");
  MatrixInverse(&AComplex);
  PrintMatrix(AComplex);
  printf("Multiply Acopy with AComplex^-1\n");
  MatrixMultiply(&AA, Acopy, AComplex);
  PrintMatrix(AA);
  printf("Add I * AComplex to A\n");
  MatrixAdd(&A, AComplex, 0, 1);
  PrintMatrix(A);
  printf("Multiply A by 3\n");
  MatrixScalarMultiply(&A, 3, 0);
  PrintMatrix(A);
  printf("Matrix Transform AB=A^* AA A^T\n");
  MatrixTransform(&AB, AA, A);
  PrintMatrix(AB);
  printf("Matrix Transform part AC=(A^* AB A^T)[2;;5,3;;7]\n");
  MatrixTransformPart(&AC, AA, A, 2, 5, 3, 7);
  PrintMatrix(AC);
  printf("AD=AB[2;;5,3;;7]\n");
  MatrixPart(&AD, AB, 2, 5, 3, 7);
  PrintMatrix(AD);

  for (i = 0; i < 12; i++)
  {
    AvecRe[i] = i;
    AvecIm[i] = 0;
    BvecRe[i] = i;
    BvecIm[i] = -1.0 * i / 10.0;
  }

  printf("Orthonormalize AD\n");
  MatrixOrthonormalize(&AD);
  printf("Matrix Copy\n");
  MatrixCopy(&ADCT, AD);
  printf("MatrixConjugateTranspose\n");
  MatrixConjugateTranspose(&ADCT);
  printf("MatrixMultiply\n");
  MatrixMultiply(&IAD, AD, ADCT);
  PrintMatrix(AD);
  PrintMatrix(IAD);

  FreeMatrix(&IAD);
  FreeMatrix(&ADCT);
  MatrixConjugateTranspose(&AD);
  printf("Orthonormalize AD\n");
  MatrixOrthonormalize(&AD);
  printf("Matrix Copy\n");
  MatrixCopy(&ADCT, AD);
  printf("MatrixConjugateTranspose\n");
  MatrixConjugateTranspose(&ADCT);
  printf("MatrixMultiply\n");
  MatrixMultiply(&IAD, AD, ADCT);
  PrintMatrix(AD);
  PrintMatrix(IAD);

  printf(" Create a matrix TA \n");
  TA.Ni = 5;
  TA.Nj = 12;
  TA.QComplex = 0;
  InitMatrix(&TA);
  strcpy(TA.Name, "TA");
  for (i = 0; i < TA.Ni; i++)
  {
    for (j = 0; j < TA.Nj; j++)
    {
      if (i == j)
      {
        TA.VRe[i][j] = 1;
      }
      else
      {
        TA.VRe[i][j] = 0;
      }
    }
  }
  PrintMatrix(TA);
  PrintMatrix(A);
  MatrixTransform(&TAATA, A, TA);
  PrintMatrix(TAATA);


}


void TestCompactMatrixVector()
{
  CompactVectorType x,y,xcom;
  int i,j,k;
  CompactMatrixType M,Mcom,Mcomfun,Mfun,AM,BM,CM,DM,EM,FM,GM,q,r;
  double Mval[20];
  
  
  

  printf("\n====== test matrix vector =========\n");
  
  x.QComplex=0;
  x.Ni=3;
  x.VRe=malloc(x.Ni*sizeof(double));
  x.VRe[0]=0.1;
  x.VRe[1]=0.2;
  x.VRe[2]=0.3;
  xcom.QComplex=1;
  xcom.Ni=3;
  xcom.V=malloc(xcom.Ni*sizeof(Complex16));
  xcom.V[0].Re=0.1;
  xcom.V[1].Re=0.2;
  xcom.V[2].Re=0.3;
  xcom.V[0].Im=0.01;
  xcom.V[1].Im=0.02;
  xcom.V[2].Im=0.03;
  
  M.Ni=4;
  M.Nj=3;
  M.QComplex=0;
  InitCompactMatrix(&M);
  k=0;
  for(i=0;i<M.Ni;i++)
  {
    for(j=0;j<M.Nj;j++)
    {
      M.VRe[i*M.Nj+j]=k;
      k++;
    }
  }
  Mcom.Ni=4;
  Mcom.Nj=3;
  Mcom.QComplex=1;
  InitCompactMatrix(&Mcom);
  k=0;
  for(i=0;i<M.Ni;i++)
  {
    for(j=0;j<M.Nj;j++)
    {
      Mcom.V[i*M.Nj+j].Re=k;
      Mcom.V[i*M.Nj+j].Im=-10*k+i;
      k++;
    }
  }
  
  printf("\nReal M times Real x\n\n");
  printf("vector x\n");
  for(i=0;i<x.Ni;i++)
  {
    printf("%17.10E \n",x.VRe[i]);
  }
  printf("\nMatrix M\n");
  PrintCompactMatrixMinimal(M);
  CompactMatrixVector(&M, &x, &y);
  printf("vector M x\n");
  for(i=0;i<y.Ni;i++)
  {
    printf("%17.10E ",y.VRe[i]);
  }
  printf("\n\n");
  free(y.VRe);
  
  printf("\nComplex M times Complex x\n\n");
  printf("vector x\n");
  for(i=0;i<xcom.Ni;i++)
  {
    printf("%17.10E %17.10E \n",xcom.V[i].Re,xcom.V[i].Im);
  }
  printf("\nMatrix M\n");
  PrintCompactMatrixMinimal(Mcom);
  CompactMatrixVector(&Mcom, &xcom, &y);
  printf("vector M x\n");
  for(i=0;i<y.Ni;i++)
  {
    printf("%17.10E %17.10E \n",y.V[i].Re,y.V[i].Im);
  }
  printf("\n\n");
  free(y.V);
  
  printf("\nComplex M times Real x\n\n");
  printf("vector x\n");
  for(i=0;i<x.Ni;i++)
  {
    printf("%17.10E\n",x.VRe[i]);
  }
  printf("\nMatrix M\n");
  PrintCompactMatrixMinimal(Mcom);
  CompactMatrixVector(&Mcom, &x, &y);
  printf("vector M x\n");
  for(i=0;i<y.Ni;i++)
  {
    printf("%17.10E %17.10E \n",y.V[i].Re,y.V[i].Im);
  }
  printf("\n\n");
  free(y.V);
  
  printf("\nReal M times Complex x\n\n");
  printf("vector x\n");
  for(i=0;i<xcom.Ni;i++)
  {
    printf("%17.10E %17.10E \n",xcom.V[i].Re,xcom.V[i].Im);
  }
  printf("\nMatrix M\n");
  PrintCompactMatrixMinimal(M);
  CompactMatrixVector(&M, &xcom, &y);
  printf("vector M x\n");
  for(i=0;i<y.Ni;i++)
  {
    printf("%17.10E %17.10E \n",y.V[i].Re,y.V[i].Im);
  }
  printf("\n\n");
  free(y.V);
  
  
  printf("\n====== test matrix eigensystem 1 =========\n");
  FreeCompactMatrix(&M);
  FreeCompactMatrix(&Mcom);
  
  M.Ni=4;
  M.Nj=4;
  M.QComplex=0;
  InitCompactMatrix(&M);
  k=0;
  for(i=0;i<M.Ni;i++)
  {
    for(j=0;j<M.Nj;j++)
    {
      M.VRe[i*M.Nj+j]=sqrt((double)(i+j))/(1+sqrt((double)(10-i-j)));
      k++;
    }
  }
  Mcom.Ni=4;
  Mcom.Nj=4;
  Mcom.QComplex=1;
  InitCompactMatrix(&Mcom);
  k=0;
  for(i=0;i<M.Ni;i++)
  {
    for(j=0;j<M.Nj;j++)
    {
      Mcom.V[i*M.Nj+j].Re=sqrt((double)(i+j))/(1+sqrt((double)(10-i-j)));
      Mcom.V[i*M.Nj+j].Im=((double)i-(double)j)/(1+sqrt((double)(10-i-j)));
      k++;
    }
  }
  
  printf("\n====== test matrix eigensystem 1 real =========\n");
  CompactMatrixEigensystem(&M, &Mfun, Mval);
  printf("\nMatrix M\n");
  PrintCompactMatrixMinimal(M);
  printf("\nMatrix Mfun\n");
  PrintCompactMatrixMinimal(Mfun);
  for(i=0;i<M.Ni;i++)
  {
    printf("%22.15E\n",Mval[i]);
  }
  printf("\n====== test matrix eigensystem 1 complex =========\n");
  CompactMatrixEigensystem(&Mcom, &Mcomfun, Mval);
  printf("\nMatrix M\n");
  PrintCompactMatrixMinimal(Mcom);
  printf("\nMatrix Mcomfun\n");
  PrintCompactMatrixMinimal(Mcomfun);
  for(i=0;i<M.Ni;i++)
  {
    printf("%22.15E\n",Mval[i]);
  }
  
  
  printf("\n====== test matrix matrix =========\n");
  
  BM.Ni=2;
  BM.Nj=3;
  BM.QComplex=0;
  AM.Ni=3;
  AM.Nj=4;
  AM.QComplex=0;
  InitCompactMatrix(&AM);
  InitCompactMatrix(&BM);
  for(i=0;i<AM.Ni*AM.Nj;i++)
  {
    AM.VRe[i]=i;
  }
  for(i=0;i<BM.Ni*BM.Nj;i++)
  {
    BM.VRe[i]=1.0/((double)i+1);
  }
  printf("matrix A (with values j + i*Nj), Ni=%i Nj=%i\n",AM.Ni,AM.Nj);
  PrintCompactMatrixMinimal(AM);
  printf("matrix B (with values 1/(1 + j + i*Nj)), Ni=%i Nj=%i\n",BM.Ni,BM.Nj);
  PrintCompactMatrixMinimal(BM);
  
  
  Mval[0]=1;
  Mval[1]=2;
  Mval[2]=3;
  Mval[3]=4;
  
  MultiplyDiagonalMatrixCompactMatrix(Mval, &AM);
  PrintCompactMatrixMinimal(AM);
  MultiplyCompactMatrixDiagonalMatrix(&AM, Mval);
  PrintCompactMatrixMinimal(AM);
  
  MakeCompactRealMatrixComplex(&AM);
  AM.V[1].Im=1;
  AM.V[2].Im=2;
  AM.V[3].Im=3;
  AM.V[4].Im=4;
  AM.V[5].Im=5;
  printf("AM with Ni=%i Nj=%i\n",AM.Ni,AM.Nj);
  PrintCompactMatrixMinimal(AM);
  
  TransformDiagonalMatrixByCompactMatrixTdagAT(Mval, &AM, &CM);
  printf("AM^dag Mval AM with Ni=%i Nj=%i\n",CM.Ni,CM.Nj);
  PrintCompactMatrixMinimal(CM);
  FreeCompactMatrix(&CM);
  
  TransformDiagonalMatrixByCompactMatrixTATdag(Mval, &AM, &CM);
  printf("AM Mval AM^dag with Ni=%i Nj=%i\n",CM.Ni,CM.Nj);
  PrintCompactMatrixMinimal(CM);
  FreeCompactMatrix(&CM);
  
  TransformDiagonalMatrixByCompactMatrixTtranATcon(Mval, &AM, &CM);
  printf("AM^T Mval AM^* with Ni=%i Nj=%i\n",CM.Ni,CM.Nj);
  PrintCompactMatrixMinimal(CM);
  FreeCompactMatrix(&CM);
  
  TransformDiagonalMatrixByCompactMatrixTconATtran(Mval, &AM, &CM);
  printf("AM^* Mval AM^tran with Ni=%i Nj=%i\n",CM.Ni,CM.Nj);
  PrintCompactMatrixMinimal(CM);
  FreeCompactMatrix(&CM);
  
  
  printf("Multiplication\n");
  PrintCompactMatrixMinimal(AM);
  PrintCompactMatrixMinimal(BM);
  
  CompactMatrixMultiply('N', 'N', &BM, &AM, &CM);
  printf("matrix C = B A, Ni=%i Nj=%i\n",CM.Ni,CM.Nj);
  PrintCompactMatrixMinimal(CM);
  FreeCompactMatrix(&CM);
  
  CompactMatrixMultiply('*', '*', &BM, &AM, &CM);
  printf("matrix C = B^* A^*, Ni=%i Nj=%i\n",CM.Ni,CM.Nj);
  PrintCompactMatrixMinimal(CM);
  FreeCompactMatrix(&CM);
  
  CompactMatrixMultiply('C', 'C', &AM, &BM, &DM);
  printf("matrix D = A^C B^C, Ni=%i Nj=%i\n",DM.Ni,DM.Nj);
  PrintCompactMatrixMinimal(DM);
  FreeCompactMatrix(&DM);
  
  CompactMatrixMultiply('T', 'T', &AM, &BM, &DM);
  printf("matrix D = A^T B^T, Ni=%i Nj=%i\n",DM.Ni,DM.Nj);
  PrintCompactMatrixMinimal(DM);
  FreeCompactMatrix(&DM);
  
   
  
  printf("\n====== test matrix inverse =========\n");
  
  AM.Nj=3;
  AM.VRe[0]+=2;
  
  FreeCompactMatrix(&BM);
  CompactMatrixCopy(AM,&BM);
  
  printf("matrix A, Ni=%i Nj=%i\n",AM.Ni,AM.Nj);
  PrintCompactMatrixMinimal(AM);
  CompactMatrixInverse(&BM);
  printf("matrix B=A^{-1}, Ni=%i Nj=%i\n",BM.Ni,BM.Nj);
  PrintCompactMatrixMinimal(BM);
  
  CompactMatrixMultiply('N', 'N', &BM, &AM, &CM);
  CompactMatrixChop(&CM, 1000*DBL_EPSILON);
  
  printf("matrix C=B A, Ni=%i Nj=%i\n",CM.Ni,CM.Nj);
  PrintCompactMatrixMinimal(CM);
  
  
  printf("\n====== test matrix eigensystem =========\n");
  
  FreeCompactMatrix(&CM);
  FreeCompactMatrix(&BM);
  CompactMatrixCopy(AM, &BM);
  CompactMatrixConjugateTranspose(&BM);
  
  
  CompactMatrixAdd(&AM, &BM);
  double val[3];
  CompactVectorType vec,vecout;
  vec.QComplex=1;
  vec.Ni=3;
  vec.V=malloc(vec.Ni*sizeof(Complex16));
  
  
  CompactMatrixEigensystem(&AM, &BM, val);
  
  printf("matrix A, Ni=%i Nj=%i\n",AM.Ni,AM.Nj);
  PrintCompactMatrixMinimal(AM);
  
  printf("matrix B=Eigenfunctions(A), Ni=%i Nj=%i\n",BM.Ni,BM.Nj);
  PrintCompactMatrixMinimal(BM);
  
  for(i=0;i<3;i++)
  {
    printf("val %i = %22.15E\n",i,val[i]);
  }
  
  for(i=0;i<3;i++)
  {
    vec.V[i]=BM.V[i];
  }
  
  CompactMatrixVector(&AM, &vec, &vecout);
  
  printf("\nmatrix A * vec = val * vec\n");
  printf(" vec               A*vec             A*vec/val\n");
  for(i=0;i<3;i++)
  {
    printf("%17.10E %17.10E %17.10E \n",vec.V[i].Re,vecout.V[i].Re,vecout.V[i].Re/val[0]);
  }
  printf("\n");
  
  
  for(i=0;i<3;i++)
  {
    vec.V[i]=BM.V[i+BM.Nj];
  }
  
  CompactMatrixVector(&AM, &vec, &vecout);
  
  printf("\nmatrix A * vec = val * vec\n");
  printf(" vec               A*vec             A*vec/val\n");
  for(i=0;i<3;i++)
  {
    printf("%17.10E %17.10E %17.10E \n",vec.V[i].Re,vecout.V[i].Re,vecout.V[i].Re/val[1]);
  }
  printf("\n");
  
  CompactMatrixMultiply3('*', 'N', 'T', &BM, &AM, &BM, &CM);
  CompactMatrixChop(&CM, 1000*DBL_EPSILON);
  
  printf("matrix C = vec^* . A . vec^T, Ni=%i Nj=%i\n",CM.Ni,CM.Nj);
  PrintCompactMatrixMinimal(CM);
  
  FreeCompactMatrix(&CM);
  TransformCompactMatrixcT(&AM, &BM, &CM);
  CompactMatrixChop(&CM, 1000*DBL_EPSILON);
  printf("matrix C = vec^* . A . vec^T, Ni=%i Nj=%i\n",CM.Ni,CM.Nj);
  PrintCompactMatrixMinimal(CM);
  

  printf("\n====== test QR decomposition Hermitian =======\n\n");
  
  CompactMatrixQRDecompositionEigenHermitian(&AM, 100*DBL_EPSILON, &q, &r);
  CompactMatrixChop(&q, 1000*DBL_EPSILON);
  CompactMatrixChop(&r, 1000*DBL_EPSILON);
  
  printf("Matrix A\n");
  PrintCompactMatrixMinimal(AM);
  
  printf("Matrix q\n");
  PrintCompactMatrixMinimal(q);
  
  printf("Matrix r\n");
  PrintCompactMatrixMinimal(r);

  FreeCompactMatrix(&BM);
  CompactMatrixMultiply('C', 'N', &q, &r, &BM);
  CompactMatrixChop(&BM, 1000*DBL_EPSILON);
  printf("Matrix q^C r\n");
  PrintCompactMatrixMinimal(BM);
  
  
  FreeCompactMatrix(&CM);
  CompactMatrixMultiply('C', 'N', &q, &q, &CM);
  CompactMatrixChop(&CM, 1000*DBL_EPSILON);
  printf("Matrix q^C q\n");
  PrintCompactMatrixMinimal(CM);
  
  
  printf("\n====== test QR decomposition Hermitian with singular matrix =======\n\n");
  
  
  FreeCompactMatrix(&q);
  FreeCompactMatrix(&r);
  MakeCompactComplexMatrixReal(&AM);
  
  AM.VRe[0]= 1;  AM.VRe[1]= 1;  AM.VRe[2]= 1;
  AM.VRe[3]= 2;  AM.VRe[4]= 2;  AM.VRe[5]= 2;
  AM.VRe[6]= 1;  AM.VRe[7]= 1;  AM.VRe[8]= 2;
  
  
  CompactMatrixQRDecompositionEigenHermitian(&AM, 100*DBL_EPSILON, &q, &r);
  CompactMatrixChop(&q, 1000*DBL_EPSILON);
  CompactMatrixChop(&r, 1000*DBL_EPSILON);
  
  printf("Matrix A\n");
  PrintCompactMatrixMinimal(AM);
  
  printf("Matrix q\n");
  PrintCompactMatrixMinimal(q);
  
  printf("Matrix r\n");
  PrintCompactMatrixMinimal(r);
  
  FreeCompactMatrix(&BM);
  CompactMatrixMultiply('C', 'N', &q, &r, &BM);
  CompactMatrixChop(&BM, 1000*DBL_EPSILON);
  printf("Matrix q^C r\n");
  PrintCompactMatrixMinimal(BM);
  
  
  FreeCompactMatrix(&CM);
  CompactMatrixMultiply('C', 'N', &q, &q, &CM);
  CompactMatrixChop(&CM, 1000*DBL_EPSILON);
  printf("Matrix q^C q\n");
  PrintCompactMatrixMinimal(CM);
  
  printf("\n====== test sqrt of matrix =======\n\n");
  
  FreeCompactMatrix(&BM);
  FreeCompactMatrix(&CM);
  
  AM.VRe[0]= 1;  AM.VRe[1]=  2;  AM.VRe[2]=  5;
  AM.VRe[3]= 2;  AM.VRe[4]=  1;  AM.VRe[5]= -1;
  AM.VRe[6]= 5;  AM.VRe[7]= -1;  AM.VRe[8]=  2;
  
  CompactMatrixMultiply('N', 'C', &AM, &AM, &BM);
  CompactMatrixSqrt(&BM, &CM);
  CompactMatrixMultiply('N', 'C', &CM, &CM, &DM);
  
  printf("Original matrix A\n");
  PrintCompactMatrixMinimal(AM);
  printf("B=A*A^C\n");
  PrintCompactMatrixMinimal(BM);
  printf("C=sqrt(B)\n");
  PrintCompactMatrixMinimal(CM);
  printf("D=C*C^C\n");
  PrintCompactMatrixMinimal(DM);
  
  FreeCompactMatrix(&BM);
  FreeCompactMatrix(&CM);
  FreeCompactMatrix(&DM);
  
  MakeCompactRealMatrixComplex(&AM);
  AM.V[0].Im= 0; AM.V[1].Im= 1; AM.V[2].Im= 2;
  AM.V[3].Im= 3; AM.V[4].Im= 4; AM.V[5].Im= 5;
  AM.V[6].Im= 6; AM.V[7].Im= 7; AM.V[8].Im= 8;
  
  CompactMatrixMultiply('N', 'C', &AM, &AM, &BM);
  CompactMatrixSqrt(&BM, &CM);
  CompactMatrixMultiply('N', 'C', &CM, &CM, &DM);
  
  printf("Original matrix A\n");
  PrintCompactMatrixMinimal(AM);
  printf("B=A*A^C\n");
  PrintCompactMatrixMinimal(BM);
  printf("C=sqrt(B)\n");
  PrintCompactMatrixMinimal(CM);
  printf("D=C*C^C\n");
  PrintCompactMatrixMinimal(DM);
  
  FreeCompactMatrix(&AM);
  FreeCompactMatrix(&BM);
  FreeCompactMatrix(&CM);
  FreeCompactMatrix(&DM);

  
  printf("\n====== test Addition =======\n\n");
  
  AM.Ni=5;
  AM.Nj=6;
  AM.QComplex=0;
  InitCompactIdentityMatrix(&AM);
  AM.VRe[29]=30;
  BM.Ni=3;
  BM.Nj=2;
  BM.QComplex=0;
  InitCompactIdentityMatrix(&BM);
  BM.VRe[5]=6;
  printf("A\n");
  PrintCompactMatrixMinimal(AM);
  printf("B\n");
  PrintCompactMatrixMinimal(BM);
  
  printf("A+B\n");
  CompactMatrixAddShifted(&AM, BM, 1.0, 0.0, 0, 0);
  PrintCompactMatrixMinimal(AM);
  
  printf("A+B + B shifted by 2 in the i direction (down)\n");
  CompactMatrixAddShifted(&AM, BM, 1.0, 0.0, 2, 0);
  PrintCompactMatrixMinimal(AM);
  
  printf("A+B + B shifted by 2 in the i direction (down) + B shifted by (1,3)\n");
  CompactMatrixAddShifted(&AM, BM, 1.0, 0.0, 1, 3);
  PrintCompactMatrixMinimal(AM);

  
  printf("\n====== test orthonormalization ROW =======\n\n");
  
  FreeCompactMatrix(&AM);
  FreeCompactMatrix(&BM);
  
  AM.Ni=4;
  AM.Nj=5;
  AM.QComplex=1;
  InitCompactMatrix(&AM);
  k=0;
  for(i=0;i<AM.Ni;i++)
  {
    for(j=0;j<AM.Nj;j++)
    {
      AM.V[k].Re=i + 20 * j ;
      AM.V[k].Im=i + 20 * j ;
      k++;
    }
  }
  AM.V[0].Re+=12;
  AM.V[6].Re+=12;
  
  CompactMatrixCopy(AM, &EM);
  
  printf("original matrix\n");
  PrintCompactMatrixMinimal(AM);
  CompactMatrixSVDOrthogonalizeRowplusT(&AM, &EM, &BM, &CM, 10000*DBL_EPSILON, True);
  
  
  printf("Orthogonalized\n");
  PrintCompactMatrixMinimal(EM);
  
  printf("T\n");
  PrintCompactMatrixMinimal(BM);
  
  printf("Tinv\n");
  PrintCompactMatrixMinimal(CM);
  
  printf("T * original matrix = orthogonalized matrix\n");
  CompactMatrixMultiply('N', 'N', &BM, &EM, &DM);
  PrintCompactMatrixMinimal(DM);
  
  printf("Tinv * orthogonal matrix = original matrix\n");
  CompactMatrixMultiply('N', 'N', &CM, &AM, &FM);
  PrintCompactMatrixMinimal(FM);
  
  printf("orthogonal matrix * orthogonal matrix^C =identity matrix\n");
  CompactMatrixMultiply('N', 'C', &AM, &AM, &GM);
  CompactMatrixChop(&GM, 1000*DBL_EPSILON);
  PrintCompactMatrixMinimal(GM);
  FreeCompactMatrix(&GM);
  
  printf("original matrix * orthogonal matrix^C = overlap matrix\n");
  CompactMatrixMultiply('N', 'C', &EM, &AM, &GM);
  CompactMatrixChop(&GM, 1000*DBL_EPSILON);
  PrintCompactMatrixMinimal(GM);

  
  printf("\n====== test orthonormalization ROW singular =======\n\n");
  
  FreeCompactMatrix(&AM);
  FreeCompactMatrix(&BM);
  FreeCompactMatrix(&CM);
  FreeCompactMatrix(&DM);
  FreeCompactMatrix(&EM);
  FreeCompactMatrix(&FM);
  FreeCompactMatrix(&GM);
  
  AM.Ni=4;
  AM.Nj=5;
  AM.QComplex=1;
  InitCompactMatrix(&AM);
  k=0;
  for(i=0;i<AM.Ni;i++)
  {
    for(j=0;j<AM.Nj;j++)
    {
      AM.V[k].Re=i + 20 * j ;
      AM.V[k].Im=i + 20 * j ;
      k++;
    }
  }
  AM.V[0].Re+=12;
  
  CompactMatrixCopy(AM, &EM);
  
  printf("original matrix\n");
  PrintCompactMatrixMinimal(AM);
  CompactMatrixSVDOrthogonalizeRowplusT(&AM, &EM, &BM, &CM, 10000*DBL_EPSILON, True);
  
  
  printf("Orthogonalized\n");
  PrintCompactMatrixMinimal(AM);
  
  printf("T\n");
  PrintCompactMatrixMinimal(BM);
  
  printf("Tinv\n");
  PrintCompactMatrixMinimal(CM);
  
  printf("T * original matrix = orthogonalized matrix\n");
  CompactMatrixMultiply('N', 'N', &BM, &EM, &DM);
  PrintCompactMatrixMinimal(DM);
  
  printf("Tinv * orthogonal matrix = original matrix\n");
  CompactMatrixMultiply('N', 'N', &CM, &AM, &FM);
  PrintCompactMatrixMinimal(FM);
  
  printf("orthogonal matrix * orthogonal matrix^C =identity matrix\n");
  CompactMatrixMultiply('N', 'C', &AM, &AM, &GM);
  CompactMatrixChop(&GM, 1000*DBL_EPSILON);
  PrintCompactMatrixMinimal(GM);
  FreeCompactMatrix(&GM);
  
  printf("original matrix * orthogonal matrix^C = overlap matrix\n");
  CompactMatrixMultiply('N', 'C', &EM, &AM, &GM);
  CompactMatrixChop(&GM, 1000*DBL_EPSILON);
  PrintCompactMatrixMinimal(GM);

  
   printf("\n====== test orthonormalization COL =======\n\n");
  
  FreeCompactMatrix(&AM);
  FreeCompactMatrix(&BM);
  FreeCompactMatrix(&CM);
  FreeCompactMatrix(&DM);
  FreeCompactMatrix(&EM);
  FreeCompactMatrix(&FM);
  FreeCompactMatrix(&GM);
  AM.Ni=6;
  AM.Nj=4;
  AM.QComplex=1;
  InitCompactMatrix(&AM);
  k=0;
  for(i=0;i<AM.Ni;i++)
  {
    for(j=0;j<AM.Nj;j++)
    {
      AM.VRe[k] =i + 20 * j ;
      AM.V[k].Re=0.0000000000001 * (i + 20 * j );
      AM.V[k].Im=0.0000000000001 * (i + 20 * j );
      if((i==j)&&(i!=AM.Ni-1)) AM.V[k].Re+=1;
      k++;
    }
  }
  AM.VRe[0]+=12;
  
  CompactMatrixCopy(AM, &EM);
  
  printf("original matrix\n");
  PrintCompactMatrixMinimal(AM);
  CompactMatrixSVDOrthogonalizeRowplusT(&AM, &EM, &BM, &CM, 10000*DBL_EPSILON, True);
  
  
  printf("Orthogonalized\n");
  PrintCompactMatrixMinimal(AM);
  
  printf("T\n");
  PrintCompactMatrixMinimal(BM);
  
  printf("Tinv\n");
  PrintCompactMatrixMinimal(CM);
  
  printf("original matrix * T = orthogonalized matrix\n");
  CompactMatrixMultiply('N', 'N', &EM, &BM, &DM);
  PrintCompactMatrixMinimal(DM);
  
  printf("orthogonal matrix * Tinv = original matrix\n");
  CompactMatrixMultiply('N', 'N', &AM, &CM, &FM);
  PrintCompactMatrixMinimal(FM);
  
  printf("orthogonal matrix^C * orthogonal matrix =identity matrix\n");
  CompactMatrixMultiply('C', 'N', &AM, &AM, &GM);
  CompactMatrixChop(&GM, 1000*DBL_EPSILON);
  PrintCompactMatrixMinimal(GM);
  
   
  printf("\n====== test orthonormalization Gram schmidt ROW =======\n\n");
  
  FreeCompactMatrix(&AM);
  FreeCompactMatrix(&BM);
  FreeCompactMatrix(&CM);
  FreeCompactMatrix(&DM);
  FreeCompactMatrix(&EM);
  FreeCompactMatrix(&FM);
  FreeCompactMatrix(&GM);
  AM.Ni=4;
  AM.Nj=6;
  AM.QComplex=0;
  InitCompactMatrix(&AM);
  k=0;
  for(i=0;i<AM.Ni;i++)
  {
    for(j=0;j<AM.Nj;j++)
    {
      AM.VRe[k] =i + 20 * j ;
      //AM.V[k].Re=0.0000000000001 * (i + 20 * j );
      //AM.V[k].Im=0.0000000000001 * (i + 20 * j );
      //if((i==j)&&(i!=AM.Ni-1)) AM.V[k].Re+=1;
      //AM.VRe[k]=0.0000000000001 * (i + 20 * j );
      //if((i==j)&&(i!=AM.Ni-2)) AM.VRe[k]+=2;
      k++;
    }
  }
  AM.VRe[1]+=0.1;
  AM.VRe[0]+=12;
  
  printf("Original\n");
  PrintCompactMatrixMinimal(AM);
  CompactMatrixGramSchmidtOrthogonalizeRow(&AM, DBL_EPSILON, True);
  printf("Orhtonormal\n");
  PrintCompactMatrixMinimal(AM);
  CompactMatrixMultiply('N', 'C', &AM, &AM, &BM);
  CompactMatrixChop(&BM, 1000*DBL_EPSILON);
  printf("AM * AM^C\n");
  PrintCompactMatrixMinimal(BM);
  
  FreeCompactMatrix(&BM);
  k=0;
  for(i=0;i<AM.Ni;i++)
  {
    for(j=0;j<AM.Nj;j++)
    {
      AM.VRe[k] =i + 20 * j ;
      //AM.V[k].Re=0.0000000000001 * (i + 20 * j );
      //AM.V[k].Im=0.0000000000001 * (i + 20 * j );
      //if((i==j)&&(i!=AM.Ni-1)) AM.V[k].Re+=1;
      //AM.VRe[k]=0.0000000000001 * (i + 20 * j );
      //if((i==j)&&(i!=AM.Ni-2)) AM.VRe[k]+=2;
      k++;
    }
  }
  AM.VRe[1]+=0.1;
  AM.VRe[0]+=12;
  
  CompactMatrixGramSchmidtOrthogonalizeRowplusT(&AM, &BM, &CM, &DM, DBL_EPSILON, True);
  
  printf("Original\n");
  PrintCompactMatrixMinimal(AM);
  
  printf("Orthogonalized\n");
  PrintCompactMatrixMinimal(BM);
  
  printf("T\n");
  PrintCompactMatrixMinimal(CM);
  
  printf("Tinv\n");
  PrintCompactMatrixMinimal(DM);
  
  CompactMatrixMultiply('N', 'N', &CM, &AM, &EM);
  printf("T * Original\n");
  PrintCompactMatrixMinimal(EM);
  
  CompactMatrixMultiply('N', 'N', &DM, &BM, &FM);
  printf("Tinv * Orthogonalized\n");
  PrintCompactMatrixMinimal(FM);
  
  
  CompactMatrixMultiply('N', 'N', &CM, &DM, &GM);
  printf("T * Tinv\n");
  PrintCompactMatrixMinimal(GM);
  
  
  printf("\n====== test orthonormalization Gram schmidt ROW Complex =======\n\n");
  
  FreeCompactMatrix(&AM);
  FreeCompactMatrix(&BM);
  FreeCompactMatrix(&CM);
  FreeCompactMatrix(&DM);
  FreeCompactMatrix(&EM);
  FreeCompactMatrix(&FM);
  FreeCompactMatrix(&GM);
  AM.Ni=4;
  AM.Nj=6;
  AM.QComplex=1;
  InitCompactMatrix(&AM);
  k=0;
  for(i=0;i<AM.Ni;i++)
  {
    for(j=0;j<AM.Nj;j++)
    {
      AM.V[k].Re =i + 20 * j ;
      AM.V[k].Im =i - 20 * j ;
      k++;
    }
  }
  AM.V[1].Re+=0.1;
  AM.V[0].Im+=12;
  
  printf("Original\n");
  PrintCompactMatrixMinimal(AM);
  CompactMatrixGramSchmidtOrthogonalizeRow(&AM, DBL_EPSILON, True);
  printf("Orhtonormal\n");
  PrintCompactMatrixMinimal(AM);
  CompactMatrixMultiply('N', 'C', &AM, &AM, &BM);
  CompactMatrixChop(&BM, 1000*DBL_EPSILON);
  printf("AM * AM^C\n");
  PrintCompactMatrixMinimal(BM);
  
  FreeCompactMatrix(&BM);
  k=0;
  for(i=0;i<AM.Ni;i++)
  {
    for(j=0;j<AM.Nj;j++)
    {
      AM.V[k].Re =i + 20 * j ;
      AM.V[k].Im =i - 20 * j ;
      k++;
    }
  }
  AM.V[1].Re+=0.1;
  AM.V[0].Im+=12;
  
  CompactMatrixGramSchmidtOrthogonalizeRowplusT(&AM, &BM, &CM, &DM, DBL_EPSILON, True);
  
  printf("Original\n");
  PrintCompactMatrixMinimal(AM);
  
  printf("Orthogonalized\n");
  PrintCompactMatrixMinimal(BM);
  
  printf("T\n");
  PrintCompactMatrixMinimal(CM);
  
  printf("Tinv\n");
  PrintCompactMatrixMinimal(DM);
  
  CompactMatrixMultiply('N', 'N', &CM, &AM, &EM);
  printf("T * Original\n");
  PrintCompactMatrixMinimal(EM);
  
  CompactMatrixMultiply('N', 'N', &DM, &BM, &FM);
  printf("Tinv * Orthogonalized\n");
  PrintCompactMatrixMinimal(FM);
  
  
  CompactMatrixMultiply('N', 'N', &CM, &DM, &GM);
  printf("T * Tinv\n");
  PrintCompactMatrixMinimal(GM);
  
  FreeCompactMatrix(&AM);
  FreeCompactMatrix(&BM);
  FreeCompactMatrix(&CM);
  FreeCompactMatrix(&DM);
  FreeCompactMatrix(&EM);
  FreeCompactMatrix(&FM);
  FreeCompactMatrix(&GM);
  
  printf("\n====== test CompactMatrixToSpectra =======\n\n");
  
  AM.Ni=6;
  AM.Nj=6;
  AM.QComplex=1;
  InitCompactMatrix(&AM);
  k=0;
  for(i=0;i<AM.Ni;i++)
  {
    for(j=0;j<AM.Nj;j++)
    {
      AM.V[k].Re =(i+j)*sqrt(i*j)/((1+i+j)*(2+i+j)) ;
      AM.V[k].Im =(i - j) ;
      //if(i==j) AM.V[k].Re=i;
      k++;
    }
  }
  
  PrintCompactMatrixMinimal(AM);
  
  SpectraType G[4];
  CompactMatrixToSpectra(&AM, G, 2, -12, 12, 120);
  PrintSpectra(G, 4);
  
  printf("\n============= uint32_t EigensystemToCompact(MatrixType const M, CompactMatrixType * const T, double *val) ============\n");
  
  MatrixType MM;
  CompactMatrixCopyToMatrix(&MM, AM);
  EigensystemToCompact(MM, &BM, Mval);
  for(i=0;i<MM.Ni;i++)
  {
    printf("%22.15E \n",Mval[i]);
  }
  PrintCompactMatrixMinimal(BM);
  FreeCompactMatrix(&BM);
  
  
  printf("\n============= uint32_t MatrixCopyToCompact(CompactMatrixType * const M, MatrixType const A) ============\n");
  
  PrintCompactMatrixMinimal(AM);
  PrintMatrix(MM);
  MatrixCopyToCompact(&BM, MM);
  PrintCompactMatrixMinimal(BM);
  
  
   /*
   
   printf("\n====== test Block band tridiagonalize =======\n\n");
   
   AM.Ni=20;
   AM.Nj=20;
   AM.QComplex=0;
   InitCompactMatrix(&AM);
   k=0;
   for(i=0;i<AM.Ni;i++)
   {
   for(j=0;j<AM.Nj;j++)
   {
   AM.VRe[k]=(i+j) / (3.0+i+j);
   k++;
   }
   }
   CompactMatrixBlockBandTridiagonalizeHermitianNoReorthogonalization(&AM, &AMTri, DBL_EPSILON, 7, 3);
   
   PrintCompactMatrixMinimal(AM);
   PrintBlockTriDiagonalMatrix(AMTri);
  
    */
  printf("\n====== ======= ======= ======= ======= ======= ======= ======= ======= ======= =======\n");
  printf("====== ======= ======= ======= ======= ======= ======= ======= ======= ======= =======\n");
  printf("====== ======= ======= ======= ======= ======= ======= ======= ======= ======= =======\n");
  printf("====== finished level one testing      ======= ======= ======= ======= ======= =======\n");
  printf("====== continue with testing functions that include the compact matrix lib     =======\n");
  printf("====== ======= ======= ======= ======= ======= ======= ======= ======= ======= =======\n");
  printf("====== ======= ======= ======= ======= ======= ======= ======= ======= ======= =======\n");
  
  TestCompactMatrixCallers();
}

void TestBlockListOfPoles()
{
  
}

/* March 2015 I had the lumious idea to change the compact matrix format from column major to row
 * major. The folowing functions used the DenseMatrix type and thus needed to be changed
 * and checked
 *
 */
void TestBlockLanczosAfterChange()
{
  printf("\n============= uint32_t BlockLanczosGroundStateConserveBasisDense(const OperatorType H, WaveFunctionType * const psi, uint32_t const Npsi) ============\n");
  printf("\n============= uint32_t BlockLanczosGroundStateConserveBasisKrylovReal(const OperatorType H, WaveFunctionType * const v, uint32_t const Npsi, ConvergenceSettingsType * Settings) ============\n");
  printf("\n============= uint32_t BlockLanczosGroundStateConserveBasisKrylovComplex(const OperatorType H, WaveFunctionType * const v, uint32_t const Npsi, ConvergenceSettingsType * Settings) ============\n");
  TestBlockLanczos();
}

void TestCompactMatrixCallers()
{
  uint32_t i;
  BlockListOfPoles RealBLP;
  RealBLP.QComplex=0;
  RealBLP.N=10;
  RealBLP.BlockSize=4;
  InitBlockListOfPoles(&RealBLP);
  for(i=0;i<10;i++)
  {
    RealBLP.VRe[i*(1+RealBLP.BlockSize*RealBLP.BlockSize)] = i+0.1;
    
  }
  
  TestBlockLanczosAfterChange();
  
  //TODO: aditional test the following functions
  
  printf("\n============= uint32_t BlockListOfPolesReduceDim(BlockListOfPoles const * A, BlockListOfPoles * const B, CompactMatrixType const *T) ============\n");
        
  printf("\n============= uint32_t BlockListOfPolesG0plusSigmaToGFixedEGrid(BlockListOfPoles const *G0, BlockListOfPoles const *Sigma, BlockListOfPoles * const G, double const Emin, double const Emax, uint32_t const NE) ============\n");
  
  printf("\n============= uint32_t BlockListOfPolesInversedToCompactMatrix(BlockListOfPoles const * const A, CompactMatrixType * const B) ============\n");
  
  printf("\n============= uint32_t BlockListOfPolesToBlockTridiagonalMatrix(BlockListOfPoles const A, BlockTriDiagonalMatrixType *const B, int const NBlockMax, double const SingularValue) ============\n");
  
  printf("\n============= uint32_t BlockListOfPolesToBlockTridiagonalMatrixReduceDimension(BlockListOfPoles const A, BlockTriDiagonalMatrixType *B, CompactMatrixType const *AT,  int const NBlockMax, double const SingularValue) ============\n");
  
  printf("\n============= uint32_t BlockListOfPolesToBlockTridiagonalMatrixReduceDimension(BlockListOfPoles const A, BlockTriDiagonalMatrixType *B, CompactMatrixType const *AT,  int const NBlockMax, double const SingularValue) ============\n");
  
  printf("\n============= uint32_t BlockListOfPolesToBlockNaturalImpurityOrbitalMatrix(BlockListOfPoles const *A, BlockNaturalImpurityOrbitalMatrixType * const B, double const kbT, uint32_t const NBlockMaxVal, uint32_t const NBlockMaxCon) ============\n");
  
  printf("\n============= uint32_t BlockListOfPolesToBlockNaturalImpurityOrbitalMatrixReduceDimension(BlockListOfPoles const *A, BlockNaturalImpurityOrbitalMatrixType * const B, CompactMatrixType const *AT, double const kbT, uint32_t const NBlockMaxVal, uint32_t const NBlockMaxCon) ============\n");
  

  
  
  printf("\n============= uint32_t TightBindingHamiltonianToBlockListOfPoles(TightBindingDef const A, int Nk[3], BlockListOfPoles * const G, double const Emin, double const Emax, uint32_t const NE) ============\n");
  
  printf("\n============= uint32_t TightBindingHamiltonianMultiSiteToBlockListOfPoles(TightBindingDef const A, int Nk[3], uint32_t const Nvec, int const *vec, BlockListOfPoles * const G, double const Emin, double const Emax, uint32_t const NE) ============\n");
  
  printf("\n============= uint32_t TightBindingHamiltonianPlusSelfEnergyToBlockListOfPoles(TightBindingDef const A, int Nk[3], BlockListOfPoles const Sigma, BlockListOfPoles * const G, double const Emin, double const Emax, uint32_t const NE) ============\n");
  
  printf("\n============= uint32_t TightBindingHamiltonianPlusSelfEnergyMultiSiteToBlockListOfPoles(TightBindingDef const A, int Nk[3], uint32_t const Nvec, int const *vec, BlockListOfPoles const Sigma, BlockListOfPoles * const G, double const Emin, double const Emax, uint32_t const NE) ============\n");
  
  printf("\n============= uint32_t TightBindingHamiltonianAddSpinPlusSelfEnergyToBlockListOfPoles(TightBindingDef const A, int Nk[3], BlockListOfPoles const Sigma, BlockListOfPoles * const G, double const Emin, double const Emax, uint32_t const NE) ============\n");
  
  printf("\n============= uint32_t TightBindingHamiltonianAddSpinPlusSelfEnergyMultiSiteToBlockListOfPoles(TightBindingDef const A, int Nk[3], uint32_t const Nvec, int const *vec, BlockListOfPoles const Sigma, BlockListOfPoles * const G, double const Emin, double const Emax, uint32_t const NE) ============\n");
  
  
  
  printf("\n============= uint32_t DiagonalMatrixBlockTridiagonalize(double const *A, CompactMatrixType const *psi, CompactMatrixType *const T, BlockTriDiagonalMatrixType *const ABlock, double NBlock, double SingularValue) ============\n");

  printf("\n============= uint32_t BlockTriDiagonalMatrixToBlockListOfPoles(BlockTriDiagonalMatrixType const * A, BlockListOfPoles * const G, double const Emin, double const Emax, uint32_t const NE) ============\n");
  
  printf("\n============= uint32_t BlockTriDiagonalMatrixToSpectra(BlockTriDiagonalMatrixType const * A, SpectraType *const G, double const Emin, double const Emax, uint32_t const N, double const Gamma) ============\n");

  printf("\n============= uint32_t BlockWaveFunctionCalculatOverlapMatrixSerial(WaveFunctionType *const psi, uint32_t const Npsi, CompactMatrixType *const S) ============\n");
  
  
  printf("\n============= uint32_t ImpurityHFromTB(OperatorType *const Himp, ModelSettings O) ============\n");

  
  
  
  
  
  
}

void TestAllMatrix()
{
  printf("Testing the matrix library and the different matrix representations\n");
  MatrixType M, A, B, BT, MI, IdentM, ATri, T, Mdiag;
  TriDiagonalMatrixType MTri;
  BlokTriDiagonalMatrixType MBlokTri;
  uint32_t i, j;
  double singular, val[20];
  A.Ni = 6;
  A.Nj = 4;
  A.QComplex = 0;
  InitMatrix(&A);
  strcpy(A.Name, "A");
  for (i = 0; i < A.Ni; i++)
  {
    for (j = 0; j < A.Nj; j++)
    {
      A.VRe[i][j] = 100 + 10 * i + j;
    }
  }

  printf("\nCreated a Real %i by %i matrix with elements value 100+10*i+j\n", A.Ni, A.Nj);
  PrintMatrix(A);

  for(i=0;i<A.Ni;i++)
  {
    for(j=0;j<A.Nj;j++)
    {
      printf("A[i=%3i][j=%3i]=%22.15E = A[0][i*A.Nj+j=%3i]=%22.15E\n",i,j,A.VRe[i][j],i*A.Nj+j,A.VRe[0][i*A.Nj+j]);
    }
  }
  

  printf("\nAdd A to B with B equal to 0.1*j+0.01*i\n");
  B.Ni = 6;
  B.Nj = 4;
  B.QComplex = 0;
  InitMatrix(&B);
  strcpy(B.Name, "B");
  for (i = 0; i < B.Ni; i++)
  {
    for (j = 0; j < B.Nj; j++)
    {
      B.VRe[i][j] = 0.1 * j + 0.01 * i;
    }
  }
  MatrixAdd(&B, A, 1, 0);
  PrintMatrix(B);
  printf("\nTranspose the matrix\n");
  MatrixTranspose(&A);
  PrintMatrix(A);
  printf("\nM = A^T * B\n");
  MatrixMultiply(&M, A, B);
  strcpy(M.Name, "M");
  PrintMatrix(M);
  
  printf("\nAdd some small numbers to the diagonal of M\n");
  M.VRe[3][3] += 0.00100;
  M.VRe[2][2] += 0.00200;
  M.VRe[1][1] += 0.00300;
  M.VRe[0][0] += 0.00400;
  PrintMatrix(M);
  MatrixCopy(&MI, M);
  singular = MatrixInverse(&MI);
  strcpy(MI.Name, "MI");
  printf("Matrix inverted with signular=%22.15E\n", singular);
  PrintMatrix(MI);
  printf("\nIdentM = M * MI\n");
  MatrixMultiply(&IdentM, M, MI);
  strcpy(IdentM.Name, "IdentM");
  PrintMatrix(IdentM);
  FreeMatrix(&IdentM);
  printf("\nIdentM = MI * M\n");
  MatrixMultiply(&IdentM, MI, M);
  strcpy(IdentM.Name, "IdentM");
  PrintMatrix(IdentM);
  FreeMatrix(&IdentM);
  FreeMatrix(&MI);
  MakeRealMatrixComplex(&M);
  
  
  printf("\nMatrix times vector\n");
  double avec[6]={1,2,3,4,5,6};
  double bvec[4];
  RealMatrixVector(A, avec, bvec);
  PrintMatrix(A);
  for(i=0;i<4;i++)
  {
    printf("%22.15E\n",bvec[i]);
  }
  
  printf("\nelement A[2][1] is %22.15E\n",A.VRe[2][1]);
  
  
  printf("\nSame but complex \n\n");
  M.VIm[0][1] = 10.5;
  M.VIm[2][0] = 10.5;
  MatrixCopy(&MI, M);
  singular = MatrixInverse(&MI);
  printf("Matrix inverted with signular=%22.15E\n", singular);
  MatrixMultiply(&IdentM, MI, M);
  PrintMatrix(M);
  PrintMatrix(MI);
  PrintMatrix(IdentM);
  singular = MatrixInverse(&MI);
  printf("Matrix inverted with signular=%22.15E\n", singular);
  MatrixAdd(&MI, M, -1, 0);
  printf("MI=((M)^{-1})^{-1})-M\n");
  PrintMatrix(MI);


  printf("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n");
  FreeMatrix(&MI);
  FreeMatrix(&M);
  FreeMatrix(&A);
  FreeMatrix(&B);
  A.Ni = 6;
  A.Nj = 6;
  A.QComplex = 0;
  strcpy(A.Name, "A");
  InitMatrix(&A);
  for (i = 0; i < A.Ni; i++)
  {
    for (j = 0; j < A.Nj; j++)
    {
      A.VRe[i][j] = 1 * i + 1 * j + 0.000001 * abs(i - j);
      if (i == j)
        A.VRe[i][j] += 10 + 0 * (i + 1);
    }
  }
  strcpy(B.Name, "B=Tri Diagonalizer of A");
  strcpy(BT.Name, "BT=Transpose of Tri Diagonalizer of A");
  strcpy(MI.Name, "BT . B");
  strcpy(MTri.Name, "MTri=Tri Diagonal form of A");
  strcpy(ATri.Name, "ATri=Tri Diagonal form of A by applying B . A . B^T");
  TriDiagonalize(A, &B, &MTri);
  MatrixCopy(&BT, B);
  MatrixTranspose(&BT);
  MatrixMultiply(&MI, BT, B);
  MatrixTransform(&ATri, A, B);
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
  A.Ni = 6;
  A.Nj = 6;
  A.QComplex = 1;
  strcpy(A.Name, "A");
  InitMatrix(&A);
  for (i = 0; i < A.Ni; i++)
  {
    for (j = 0; j < A.Nj; j++)
    {
      A.VRe[i][j] = 1 * i + 1 * j + 0.000001 * abs(i - j);
      A.VIm[i][j] = 1 * i + 1 * j + 0.000001 * abs(i - j);
      if (i == j)
        A.VIm[i][j] = 0;
      if (i < j)
        A.VIm[i][j] *= -1;
      if (i == j)
        A.VRe[i][j] += 10 + 0 * (i + 1);
    }
  }
  strcpy(B.Name, "B=Tri Diagonalizer of A");
  strcpy(BT.Name, "BT=Conugate Transpose of Tri Diagonalizer of A");
  strcpy(MI.Name, "B^+ . B");
  strcpy(MTri.Name, "MTri=Tri Diagonal form of A");
  strcpy(ATri.Name, "ATri=Tri Diagonal form of A by applying B . A . B^+");
  TriDiagonalize(A, &B, &MTri);
  MatrixCopy(&BT, B);
  printf("BT<-B\n");
  PrintMatrix(BT);
  MatrixConjugateTranspose(&BT);
  printf("BT<-BT^+\n");
  PrintMatrix(BT);
  MatrixMultiply(&MI, BT, B);
  PrintMatrix(MI);

  printf("\n\n\n\n Tri diagonalization of a complex matrix \n\n\n");
  MatrixTransform(&ATri, A, B);
  PrintMatrix(A);
  PrintMatrix(B);
  PrintMatrix(ATri);
  PrintTriDiagonalMatrix(MTri);
  TriDiagonalMatrixToEigensystem(MTri, &T, val);
  PrintMatrix(T);
  for (i = 0; i < 6; i++)
  {
    printf("%22.15E \n", val[i]);
  }

  printf("\n\n\n\n\n  Start Eigensystem of a matrix  \n\n\n\n\n");
  FreeMatrix(&T);
  FreeMatrix(&A);
  A.Ni = 12;
  A.Nj = 12;
  A.QComplex = 0;
  InitMatrix(&A);
  strcpy(A.Name, "A");
  for (i = 0; i < A.Ni; i++)
  {
    for (j = 0; j < A.Nj; j++)
    {
      A.VRe[i][j] = i + j + (i * j) / (i - j + 0.5) + (i * j) / (i + j + 0.5);
      A.VRe[j][i] = A.VRe[i][j];
    }
  }


  Eigensystem(A, &T, val);
  PrintMatrix(A);
  PrintMatrix(T);
  
  MatrixTransform(&Mdiag, A, T);
  PrintMatrix(Mdiag);
  FreeMatrix(&Mdiag);

  
  printf("\n\n\n\n\n  Start Eigensystem of a complex matrix  \n\n\n\n\n");
  FreeMatrix(&T);
  FreeMatrix(&A);
  A.Ni = 12;
  A.Nj = 12;
  A.QComplex = 1;
  InitMatrix(&A);
  strcpy(A.Name, "A");
  for (i = 0; i < A.Ni; i++)
  {
    for (j = 0; j < A.Nj; j++)
    {
      A.VRe[i][j] = i + j + (i * j) / (i - j + 0.5) + (i * j) / (i + j + 0.5);
      A.VRe[j][i] = A.VRe[i][j];
      A.VIm[i][j] =((double)i-(double)j)*sqrt((double)(i+j));
    }
  }
  Eigensystem(A, &T, val);
  PrintMatrix(A);
  PrintMatrix(T);
  
  
  
  MatrixTransform(&Mdiag, A, T);
  PrintMatrix(Mdiag);
  FreeMatrix(&Mdiag);
  FreeMatrix(&T);


  printf("\n\n\n  Start Blok Tri-Diagonalize a matrix \n\n\n");

  BlokTriDiagonalize(A, &T, &MBlokTri, 4);
  
  MatrixTransform(&ATri, A, T);
  MatrixChop(&ATri, 1000*DBL_EPSILON);
  
  
  PrintMatrix(A);
  PrintMatrix(T);
  PrintMatrix(ATri);
  PrintBlokTriDiagonalMatrix(MBlokTri);
  
  printf("End of Test All Matrix \n");

}


void TestInterpolate()
{
  double x[11], y[11], ac, bc, cc, a, b;
  int    i, N = 10;
  for (i = 0; i <= 10; i++)
  {
    x[i] = 0.1 * i * i - 0.235 * i + 0.5;
    y[i] = 3 * Sqr(x[i]) + 2 * x[i] - 5;
  }
  LinearInterpolate(x, y, N, &a, &b);
  QuadraticInterpolate(x, y, N, &ac, &bc, &cc);
  printf("%22.15E %22.15E\n", a, b);
  printf("%22.15E %22.15E %22.15E\n", ac, bc, cc);
}

void TestMatrixInverse()
{
  MatrixType M;
  M.QComplex = 1;
  M.Ni = 4;
  M.Nj = 4;
  strcpy(M.Name, "M");
  double **MRe;
  double **MIm;
  int    i, j, N;
  double singular;
  N = 4;
  MRe = MatrixCalloc(N, N);
  MIm = MatrixCalloc(N, N);
  MRe[0][2] = 2.0;
  MRe[1][1] = 3.0;
  MRe[2][0] = 4.0;
  MRe[3][3] = 2.0;
  MIm[3][2] = 2.0;


  MRe[0][0] = 2.0;
  MRe[1][1] = 3.0;
  MRe[2][3] = 4.0;
  MRe[3][2] = 2.0;
  MIm[3][2] = 2.0;

  MRe[0][0] = 0.5000001;
  MRe[0][1] = 0.6;
  MRe[0][2] = 0.7;
  MRe[0][3] = 0.8;
  MRe[0][4] = 0.9;
  MIm[0][2] = 5.0;

  MRe[1][0] = 0.5;
  MRe[1][1] = 0.6;
  MRe[1][2] = 0.7;
  MRe[1][3] = 0.8;
  MRe[1][4] = 0.9;
  MIm[1][2] = 5.0;

  M.VIm = MIm;
  M.VRe = MRe;

  for (i = 0; i < N; i++)
  {
    for (j = 0; j < N; j++)
    {
      printf("(%12.5f %12.5f)", MRe[i][j], MIm[i][j]);
    }
    printf("\n");
  }
  printf("\n");
  singular = MatrixInverse(&M);
  printf("Singular=%22.15E\n", singular);
  for (i = 0; i < N; i++)
  {
    for (j = 0; j < N; j++)
    {
      printf("(%12.5E %12.5E)", MRe[i][j], MIm[i][j]);
    }
    printf("\n");
  }
  printf("\n");
  singular = MatrixInverse(&M);
  printf("Singular=%22.15E\n", singular);
  for (i = 0; i < N; i++)
  {
    for (j = 0; j < N; j++)
    {
      printf("(%12.5f %12.5f)", MRe[i][j], MIm[i][j]);
    }
    printf("\n");
  }
}

void TestBasicMath()
{
  double a, b;
  int    l1, l2, l3, m1, m2, m3;
  for (l1 = 0; l1 < 50; l1++)
  {
    a = Factorial(l1);
    printf("%4i!= %22.15E\n", l1, a);
  }
  l1 = 0;
  l2 = 0;
  l3 = 0;
  m1 = 0;
  m2 = 0;
  m3 = 0;
  a = ThreeJSymbol(l1, m1, l2, m2, l3, m3);
  b = 1.0;
  printf("(%4i %4i %4i)\n(%4i %4i %4i) =%22.15E (%22.15E)\n", l1, l2, l3, m1, m2, m3, a, b);

  l1 = 1;
  l2 = 0;
  l3 = 1;
  m1 = 0;
  m2 = 0;
  m3 = 0;
  a = ThreeJSymbol(l1, m1, l2, m2, l3, m3);
  b = -sqrt(1.0 / 3.0);
  printf("(%4i %4i %4i)\n(%4i %4i %4i) =%22.15E (%22.15E)\n", l1, l2, l3, m1, m2, m3, a, b);

  l1 = 2;
  l2 = 0;
  l3 = 2;
  m1 = 0;
  m2 = 0;
  m3 = 0;
  a = ThreeJSymbol(l1, m1, l2, m2, l3, m3);
  b = sqrt(1.0 / 5.0);
  printf("(%4i %4i %4i)\n(%4i %4i %4i) =%22.15E (%22.15E)\n", l1, l2, l3, m1, m2, m3, a, b);

  l1 = 2;
  l2 = 2;
  l3 = 2;
  m1 = 0;
  m2 = 0;
  m3 = 0;
  a = ThreeJSymbol(l1, m1, l2, m2, l3, m3);
  b = -sqrt(2.0 / 35.0);
  printf("(%4i %4i %4i)\n(%4i %4i %4i) =%22.15E (%22.15E)\n", l1, l2, l3, m1, m2, m3, a, b);

  l1 = 2;
  l2 = 2;
  l3 = 4;
  m1 = 1;
  m2 = 2;
  m3 = -3;
  a = ThreeJSymbol(l1, m1, l2, m2, l3, m3);
  b = -sqrt(1.0 / 18.0);
  printf("(%4i %4i %4i)\n(%4i %4i %4i) =%22.15E (%22.15E)\n", l1, l2, l3, m1, m2, m3, a, b);

  l1 = 7;
  l2 = 5;
  l3 = 4;
  m1 = 1;
  m2 = 2;
  m3 = -3;
  a = ThreeJSymbol(l1, m1, l2, m2, l3, m3);
  b = 7 * sqrt(5.0 / 14586.0);
  printf("(%4i %4i %4i)\n(%4i %4i %4i) =%22.15E (%22.15E)\n", l1, l2, l3, m1, m2, m3, a, b);

  l1 = 14;
  l2 = 22;
  l3 = 9;
  m1 = 13;
  m2 = -5;
  m3 = -8;
  a = ThreeJSymbol(l1, m1, l2, m2, l3, m3);
  b = -(229.0 / 6.0) * sqrt(1.0 / 9866335845905.0);
  printf("(%4i %4i %4i)\n(%4i %4i %4i) =%22.15E (%22.15E)\n", l1, l2, l3, m1, m2, m3, a, b);
}

void TestDGESVD3()
{
  const char *JOBU = "A";
  const char *JOBVT = "A";
  int    M = 4;
  int    N = 3;
  double A[N][M];
  int    LDA = M;
  double S[M];
  double U[M][M];
  int    LDU = M;
  double VT[N][N];
  int    LDVT = N;
  int    LWORK = N * M * 10;
  double WORK[LWORK];
  int    INFO;

  int    i, j;
  for (i = 0; i < N; i++)
  {
    for (j = 0; j < M; j++)
    {
      A[i][j] = 5 * i + j;
    }
  }
  printf("A=\n");
  for (i = 0; i < N; i++)
  {
    for (j = 0; j < M; j++)
    {
      printf("%12.5e ", A[i][j]);;
    }
    printf("\n");
  }
  dgesvd(JOBU, JOBVT, &M, &N, A, &LDA, S, U, &LDU, VT, &LDVT, WORK, &LWORK, &INFO);
  printf("U=\n");
  for (i = 0; i < M; i++)
  {
    for (j = 0; j < M; j++)
    {
      printf("%12.5e ", U[i][j]);;
    }
    printf("\n");
  }
  printf("VT=\n");
  for (i = 0; i < N; i++)
  {
    for (j = 0; j < N; j++)
    {
      printf("%12.5e ", VT[i][j]);;
    }
    printf("\n");
  }
  printf("S=\n");
  for (i = 0; i < Max(N, M); i++)
  {
    printf("%12.5e ", S[i]);
  }
  printf("\n");
}


void TestDGESVD2()
{
  const char *JOBU = "A";
  const char *JOBVT = "A";
  int    M = 3;
  int    N = 4;
  double A[N][M];
  int    LDA = M;
  double S[M];
  double U[M][M];
  int    LDU = M;
  double VT[N][N];
  int    LDVT = N;
  int    LWORK = N * M * 10;
  double WORK[LWORK];
  int    INFO;

  int    i, j;
  for (i = 0; i < N; i++)
  {
    for (j = 0; j < M; j++)
    {
      A[i][j] = 5 * i + j;
    }
  }
  printf("A=\n");
  for (i = 0; i < N; i++)
  {
    for (j = 0; j < M; j++)
    {
      printf("%12.5e ", A[i][j]);;
    }
    printf("\n");
  }
  dgesvd(JOBU, JOBVT, &M, &N, A, &LDA, S, U, &LDU, VT, &LDVT, WORK, &LWORK, &INFO);
  printf("U=\n");
  for (i = 0; i < M; i++)
  {
    for (j = 0; j < M; j++)
    {
      printf("%12.5e ", U[i][j]);;
    }
    printf("\n");
  }
  printf("VT=\n");
  for (i = 0; i < N; i++)
  {
    for (j = 0; j < N; j++)
    {
      printf("%12.5e ", VT[i][j]);;
    }
    printf("\n");
  }
  printf("S=\n");
  for (i = 0; i < Max(N, M); i++)
  {
    printf("%12.5e ", S[i]);
  }
  printf("\n");
  TestDGESVD3();
}


void TestDGESVD()
{
  const char *JOBU = "A";
  const char *JOBVT = "A";
  int    M = 2;
  int    N = 5;
  double A[N][M];
  int    LDA = M;
  double S[M];
  double U[M][M];
  int    LDU = M;
  double VT[N][N];
  int    LDVT = N;
  int    LWORK = N * M * 10;
  double WORK[LWORK];
  int    INFO;

  int    i, j;
  for (i = 0; i < N; i++)
  {
    for (j = 0; j < M; j++)
    {
      A[i][j] = 120 + 10 * i + j;
    }
  }
  printf("A string\n");
  for (i = 0; i < M * N; i++)
  {
    printf("%12.5e ", A[0][i]);
  }
  printf("\n");
  printf("A=\n");
  for (i = 0; i < N; i++)
  {
    for (j = 0; j < M; j++)
    {
      printf("%12.5e ", A[i][j]);
    }
    printf("\n");
  }
  dgesvd(JOBU, JOBVT, &M, &N, A, &LDA, S, U, &LDU, VT, &LDVT, WORK, &LWORK, &INFO);
  printf("U=\n");
  for (i = 0; i < M; i++)
  {
    for (j = 0; j < M; j++)
    {
      printf("%12.5e ", U[i][j]);;
    }
    printf("\n");
  }
  printf("VT=\n");
  for (i = 0; i < N; i++)
  {
    for (j = 0; j < N; j++)
    {
      printf("%12.5e ", VT[i][j]);;
    }
    printf("\n");
  }
  printf("S=\n");
  for (i = 0; i < Max(N, M); i++)
  {
    printf("%12.5e ", S[i]);
  }
  printf("\n");
  TestDGESVD2();
}

void TestComplex()
{

  int    n, inca = 1, incb = 1, i, N = 5;
  Complex16 a[N], b[N], c;
  double aRe[N],aIm[N],bRe[N],bIm[N],cRe,cIm,cReIm,cImRe;
  n = N;
  for (i = 0; i < n; i++)
  {
    a[i].Re = (double)i;
    a[i].Im = (double)i *2.0;
    b[i].Re = (double)(n - i);
    b[i].Im = (double)i *2.0;
    
    aRe[i] = (double)i;
    aIm[i] = (double)i *2.0;
    bRe[i] = (double)(n - i);
    bIm[i] = (double)i *2.0;
  }
  zdotc(&c, &n, a, &inca, b, &incb);
  cRe=ddot(&n, aRe, &inca, bRe, &incb);
  cIm=ddot(&n, aIm, &inca, bIm, &incb);
  cImRe=ddot(&n, aIm, &inca, bRe, &incb);
  cReIm=ddot(&n, aRe, &inca, bIm, &incb);
  printf("The complex dot product is: ( %6.2f, %6.2f) \n", c.Re, c.Im);
  printf("The complex dot product is: ( %6.2f, %6.2f) \n", cRe+cIm, cReIm-cImRe);

}

void TestTriDiagonalMatrix()
{
  TriDiagonalMatrixType A, B;
  AndersonMatrixType C;
  DoubleTriDiagonalMatrixType D;
  ListOfPoles E;
  MatrixType T, T2;
  uint32_t i;
  A.N = 4;
  InitTriDiagonalMatrix(&A);
  A.a[0] = 0;
  A.b[0] = 1;
  for (i = 1; i < A.N; i++)
  {
    A.a[i] = i;
    A.b[i] = i * 10;
  }
  A.a[A.N] = A.N;

  PrintTriDiagonalMatrix(A);
  TriDiagonalMatrixPart(&B, A, 2, A.N);
  PrintTriDiagonalMatrix(B);
  FreeTriDiagonalMatrix(&B);
  printf("\n\n\n\n\n");
  TriDiagonalMatrixToAndersonMatrixPlusT(A, &T, &C);

  PrintAndersonMatrix(C);

  PrintMatrix(T);


  printf("\n\n\n\n\n");

  AndersonMatrixToDoubleTriDiagonalMatrixPlusT(C, &T2, &D, 0.0, 0.0001);
  PrintDoubleTriDiagonalMatrix(D);
  PrintMatrix(T2);


  printf("\n\n\n\n\n");
  A.a[0] = 3;
  A.b[0] = 5;
  for (i = 1; i < A.N; i++)
  {
    A.a[i] = i;
    A.b[i] = i * 10;
  }
  A.a[A.N] = A.N;
  TriDiagonalMatrixToListOfPoles(A, &E);
  ListOfPolesToTridiagonalMatrix(E, &B, A.N);
  PrintTriDiagonalMatrix(A);
  PrintTriDiagonalMatrix(B);
  PrintListOfPoles(E);
}

/*
 *
 */
void TestGeneralizedEigensystem()
{
  uint32_t i, j;
  MatrixType H, H1, H2, H12;
  MatrixType HNO, H1NO, H2NO, H12NO;
  MatrixType S2, B2, B2T, temp1, temp2, temp3;
  MatrixType funH2, funH2NO, S2root;
  double valH2[5], valH2NO[5];
  H.QComplex = 0;
  H.Ni = 7;
  H.Nj = 7;
  H1.QComplex = 0;
  H1.Ni = 2;
  H1.Nj = 2;
  H2.QComplex = 0;
  H2.Ni = 5;
  H2.Nj = 5;
  H12.QComplex = 0;
  H12.Ni = 2;
  H12.Nj = 5;
  HNO.QComplex = 0;
  HNO.Ni = 7;
  HNO.Nj = 7;
  H1NO.QComplex = 0;
  H1NO.Ni = 2;
  H1NO.Nj = 2;
  H2NO.QComplex = 0;
  H2NO.Ni = 5;
  H2NO.Nj = 5;
  H12NO.QComplex = 0;
  H12NO.Ni = 2;
  H12NO.Nj = 5;
  S2.QComplex = 0;
  S2.Ni = 5;
  S2.Nj = 5;
  B2.QComplex = 0;
  B2.Ni = 5;
  B2.Nj = 5;
  InitMatrix(&H);
  strcpy(H.Name, "H    ");
  InitMatrix(&H1);
  strcpy(H1.Name, "H1   ");
  InitMatrix(&H2);
  strcpy(H2.Name, "H2   ");
  InitMatrix(&H12);
  strcpy(H12.Name, "H12  ");
  InitMatrix(&HNO);
  strcpy(HNO.Name, "HNO  ");
  InitMatrix(&H1NO);
  strcpy(H1NO.Name, "H1NO ");
//  InitMatrix(&H2NO);   strcpy(H2NO.Name,"H2NO ");
  InitMatrix(&H12NO);
  strcpy(H12NO.Name, "H12NO");
//  InitMatrix(&S2);     strcpy(S2.Name,"Overlap S2");
  InitMatrix(&B2);
  strcpy(B2.Name, "Basis B2");
  H1.VRe[0][0] = -2;
  H1.VRe[1][1] = -1;
  H1.VRe[0][1] = 1;
  H1.VRe[1][0] = 1;
  H2.VRe[0][0] = 1;
  for (i = 1; i < 5; i++)
  {
    H2.VRe[i - 1][i] = 1;
    H2.VRe[i][i - 1] = 1;
    H2.VRe[i][i] = i + 1;
  }
  H12.VRe[0][0] = 1;
  H12.VRe[0][3] = -1;
  H12.VRe[1][1] = 1;
  H12.VRe[1][3] = -1;
  H12.VRe[1][4] = 1;
  H1NO.VRe[0][0] = -2;
  H1NO.VRe[1][1] = -1;
  H1NO.VRe[0][1] = 1;
  H1NO.VRe[1][0] = 1;
  for (i = 0; i < 5; i++)
  {
    for (j = 0; j <= i; j++)
    {
      B2.VRe[i][j] = 1;
    }
  }
  MatrixCopy(&B2T, B2);
  MatrixTranspose(&B2T);
  MatrixMultiply(&S2, B2, B2T);
  strcpy(S2.Name, "Overlap S2");
  MatrixMultiply(&temp1, H2, B2T);
  MatrixMultiply(&H2NO, B2, temp1);
  FreeMatrix(&temp1);
  strcpy(H2NO.Name, "H2NO ");

  PrintMatrix(B2);
  PrintMatrix(S2);
  PrintMatrix(H2);
  PrintMatrix(H2NO);

  Eigensystem(H2, &funH2, valH2);
  GeneralizedEigensystem(H2NO, S2, &funH2NO, &S2root, valH2NO);

  PrintMatrix(funH2);
  PrintMatrix(funH2NO);
  PrintMatrix(S2root);
  for (i = 0; i < 5; i++)
    printf("i=%i Energy=%lf, %lf\n", i, valH2[i], valH2NO[i]);

  MatrixCopy(&temp1, funH2);
  MatrixTranspose(&temp1);
  MatrixMultiply(&temp2, H2, temp1);
  MatrixMultiply(&temp3, funH2, temp2);
  MatrixChop(&temp3, 10E-15);
  strcpy(temp3.Name, "funH2 . H2 . funH2^T");
  PrintMatrix(temp3);
  FreeMatrix(&temp2);
  FreeMatrix(&temp3);
  MatrixMultiply(&temp2, funH2, temp1);
  MatrixChop(&temp2, 10E-15);
  strcpy(temp2.Name, "funH2 . funH2^T");
  PrintMatrix(temp2);
  FreeMatrix(&temp1);
  FreeMatrix(&temp2);

  MatrixCopy(&temp1, funH2NO);
  MatrixTranspose(&temp1);
  MatrixMultiply(&temp2, H2NO, temp1);
  MatrixMultiply(&temp3, funH2NO, temp2);
  MatrixChop(&temp3, 10E-15);
  strcpy(temp3.Name, "funH2NO . H2NO . funH2NO^T");
  PrintMatrix(temp3);
  FreeMatrix(&temp2);
  FreeMatrix(&temp3);
  MatrixMultiply(&temp2, S2, temp1);
  MatrixMultiply(&temp3, funH2NO, temp2);
  MatrixChop(&temp3, 10E-15);
  strcpy(temp3.Name, "funH2NO . S2 . funH2NO^T");
  PrintMatrix(temp3);
  FreeMatrix(&temp2);
  FreeMatrix(&temp3);
  MatrixMultiply(&temp2, funH2NO, temp1);
  MatrixChop(&temp2, 10E-15);
  strcpy(temp2.Name, "funH2NO . funH2NO^T");
  PrintMatrix(temp2);
  FreeMatrix(&temp1);
  FreeMatrix(&temp2);


}


void TestCombinatorial()
{
  uint32_t i, j, NV, order, POList[50], index;
  order = 3;
  NV = 4;
  for (i = 0; i < PerturbationOptions(NV, order); i++)
  {
    PerturbationOrderIndexToList(i, POList, NV, order);
    PerturbationOrderListToIndex(POList, &index, NV, order);
    printf("%i %i || ", i, index);
    for (j = 0; j < NV; j++)
    {
      printf("%i ", POList[j]);
    }
    printf("\n");
  }

}

void TestThreeJSymbol()
{
  int    l1, m1, l2, m2, l3, m3, lmax = 2;
  double tj;

  for (l1 = 0; l1 <= lmax; l1++)
  {
    for (m1 = -l1; m1 <= l1; m1++)
    {
      for (l2 = 0; l2 <= lmax; l2++)
      {
        for (m2 = -l2; m2 <= l2; m2++)
        {
          for (l3 = 0; l3 <= lmax; l3++)
          {
            for (m3 = -l3; m3 <= l3; m3++)
            {
              tj = ThreeJSymbol(l1, m1, l2, m2, l3, m3);
              printf("%i %2i %i %2i %i %2i %22.15E\n", l1, m1, l2, m2, l3, m3, tj);
            }
          }
        }
      }
    }
  }

}
