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
      AComplex.VIm[i][j] = ((double)i - (double)j);
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
  double vec[10]={0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,1.1};
  double vecout[10]={1,1,1,1,1,1,1,1,1,1};
  int i,j,k;
  CompactMatrixType IM,AM,BM,CM,DM,EM,FM,GM,q,r;
  BlockTriDiagonalMatrixType AMTri;
  
  printf("\n====== test matrix vector =========\n");
  
  IM.Ni=2;
  IM.Nj=3;
  IM.QComplex=0;
  InitCompactIdentityMatrix(&IM);
  printf("vector x\n");
  for(i=0;i<2;i++)
  {
    printf("%17.10E \n",vec[i]);
  }
  printf("\nMatrix I\n");
  PrintCompactMatrixMinimal(IM);
  CompactMatrixVector('N',&IM, vec, vecout);
  printf("vector I x\n");
  for(i=0;i<3;i++)
  {
    printf("%17.10E ",vecout[i]);
  }
  printf("\n\n");
  
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
  
  CompactMatrixMultiply('N', 'N', &AM, &BM, &CM);
  
  printf("matrix C = A B, Ni=%i Nj=%i\n",CM.Ni,CM.Nj);
  PrintCompactMatrixMinimal(CM);
  
  
  
  CompactMatrixMultiply('C', 'C', &BM, &AM, &DM);
  
  printf("matrix D = B^T A^T, Ni=%i Nj=%i\n",DM.Ni,DM.Nj);
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
    vec[i]=BM.VRe[i];
  }
  
  CompactMatrixVector('N',&AM, vec, vecout);
  
  printf("\nmatrix A * vec = val * vec\n");
  printf(" vec               A*vec             A*vec/val\n");
  for(i=0;i<3;i++)
  {
    printf("%17.10E %17.10E %17.10E \n",vec[i],vecout[i],vecout[i]/val[0]);
  }
  printf("\n");
  
  
  for(i=0;i<3;i++)
  {
    vec[i]=BM.VRe[i+BM.Nj];
  }
  
  CompactMatrixVector('N',&AM, vec, vecout);
  
  printf("\nmatrix A * vec = val * vec\n");
  printf(" vec               A*vec             A*vec/val\n");
  for(i=0;i<3;i++)
  {
    printf("%17.10E %17.10E %17.10E \n",vec[i],vecout[i],vecout[i]/val[1]);
  }
  printf("\n");
  
  CompactMatrixMultiply3('T', 'N', 'N', &BM, &AM, &BM, &CM);
  CompactMatrixChop(&CM, 1000*DBL_EPSILON);
  
  printf("matrix C = vec^T . A . vec, Ni=%i Nj=%i\n",CM.Ni,CM.Nj);
  PrintCompactMatrixMinimal(CM);
  
  printf("\n\nQR decomposition\n\n");
  
  CompactMatrixQRDecompositionEigenHermitian(&AM, 100*DBL_EPSILON, &q, &r);
  
  printf("Matrix A\n");
  PrintCompactMatrixMinimal(AM);
  
  printf("Matrix q\n");
  PrintCompactMatrixMinimal(q);
  
  printf("Matrix r\n");
  PrintCompactMatrixMinimal(r);

  FreeCompactMatrix(&BM);
  CompactMatrixMultiply('N', 'N', &q, &r, &BM);
  
  printf("Matrix q r\n");
  PrintCompactMatrixMinimal(BM);
  
  
  FreeCompactMatrix(&CM);
  CompactMatrixMultiply('T', 'N', &q, &q, &CM);
  
  
  printf("Matrix q^T q\n");
  PrintCompactMatrixMinimal(CM);
  
  
  
  
  
  printf("\n\nQR decomposition with singular matrix\n\n");
  
  
  FreeCompactMatrix(&q);
  FreeCompactMatrix(&r);
  
  AM.VRe[0]= 1;  AM.VRe[3]= 1;  AM.VRe[6]= 1;
  AM.VRe[1]= 1;  AM.VRe[4]= 1;  AM.VRe[7]= 1;
  AM.VRe[2]= 1;  AM.VRe[5]= 1;  AM.VRe[8]= 2;
  
  
  CompactMatrixQRDecompositionEigenHermitian(&AM, 100*DBL_EPSILON, &q, &r);
  
  printf("Matrix A\n");
  PrintCompactMatrixMinimal(AM);
  
  printf("Matrix q\n");
  PrintCompactMatrixMinimal(q);
  
  printf("Matrix r\n");
  PrintCompactMatrixMinimal(r);
  
  FreeCompactMatrix(&BM);
  CompactMatrixMultiply('N', 'N', &q, &r, &BM);
  
  printf("Matrix q r\n");
  PrintCompactMatrixMinimal(BM);
  
  
  FreeCompactMatrix(&CM);
  CompactMatrixMultiply('T', 'N', &q, &q, &CM);
  
  
  printf("Matrix q^T q\n");
  PrintCompactMatrixMinimal(CM);
  
  
  printf("\n\nsqrt of a matrix\n\n");
  
  
  AM.VRe[0]= 1;  AM.VRe[3]= 2;  AM.VRe[6]= 5;
  AM.VRe[1]= 2;  AM.VRe[4]= 1;  AM.VRe[7]= -1;
  AM.VRe[2]= 5;  AM.VRe[5]= -1;  AM.VRe[8]= 2;
  
  FreeCompactMatrix(&BM);
  CompactMatrixMultiply('N', 'C', &AM, &AM, &BM);
  
  FreeCompactMatrix(&CM);
  
  CompactMatrixSqrt(&BM, &CM);
  
  CompactMatrixMultiply('N', 'C', &CM, &CM, &DM);
  
  PrintCompactMatrixMinimal(AM);
  PrintCompactMatrixMinimal(BM);
  PrintCompactMatrixMinimal(CM);
  PrintCompactMatrixMinimal(DM);
  
  
  FreeCompactMatrix(&AM);
  FreeCompactMatrix(&BM);
  FreeCompactMatrix(&CM);
  FreeCompactMatrix(&DM);
  
  
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
  
  printf("\n\n Test orthonormalization \n\n");
  
  FreeCompactMatrix(&AM);
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
  AM.VRe[0]+=12;
  
  CompactMatrixCopy(AM, &EM);
  
  printf("original matrix\n");
  PrintCompactMatrixMinimal(AM);
  CompactMatrixSVDOrthogonalizeRow(&AM, 10000*DBL_EPSILON, &BM, &CM);
  
  
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

  
  
}


void TestAllMatrix()
{
  printf("Testing the matrix library and the different matrix representations\n");
  MatrixType M, A, B, BT, MI, IdentM, ATri, T;
  TriDiagonalMatrixType MTri;
  BlokTriDiagonalMatrixType MBlokTri;
  uint32_t i, j;
  double singular, val[20];
  A.Ni = 6;
  A.Nj = 4;
  A.QComplex = 0;
  strcpy(A.Name, "A");
  InitMatrix(&A);
  for (i = 0; i < A.Ni; i++)
  {
    for (j = 0; j < A.Nj; j++)
    {
      A.VRe[i][j] = 100 + 10 * i + j;
    }
  }

  printf("\nCreated a Real %i by %i matrix with elements value 100+10*i+j\n", A.Ni, A.Nj);
  PrintMatrix(A);


  printf("\nAdd A to B with B equal to 0.1*j+0.01*i\n");
  B.Ni = 6;
  B.Nj = 4;
  B.QComplex = 0;
  strcpy(B.Name, "B");
  InitMatrix(&B);
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
  printf("\nM = A B\n");
  strcpy(M.Name, "M");
  MatrixMultiply(&M, A, B);
  PrintMatrix(M);

  M.VRe[3][3] += 0.00100;
  M.VRe[2][2] += 0.00200;
  M.VRe[1][1] += 0.00300;
  M.VRe[0][0] += 0.00400;
  PrintMatrix(M);
  MatrixCopy(&MI, M);
  strcpy(MI.Name, "MI");
  singular = MatrixInverse(&MI);
  printf("Matrix inverted with signular=%22.15E\n", singular);
  PrintMatrix(MI);
  printf("\nIdentM = M MI\n");
  strcpy(IdentM.Name, "IdentM");
  MatrixMultiply(&IdentM, M, MI);
  PrintMatrix(IdentM);
  FreeMatrix(&IdentM);
  printf("\nIdentM = MI M\n");
  MatrixMultiply(&IdentM, MI, M);
  PrintMatrix(IdentM);
  FreeMatrix(&IdentM);
  FreeMatrix(&MI);
  MakeRealMatrixComplex(&M);
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
  strcpy(A.Name, "A");
  InitMatrix(&A);
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
