#include "ForCore.h"








/* BasicMath_MatrixVector.c
 *
 * Issue?:  I've found that I have to calculate T^{trans}AT to get the same results as the non-BLAS version which
 * calculates TAT^{trans}  (I haven't had a complex case to test the conjugation yet)
 * This must be an issue with the conversion between matrix types or something?
 */
/* -- MatrixTransform changed to include BLAS - when it gives a gain. (T^* A T^T with A diagonal is faster without blass and even with A dense is not sure. For Real I can use BLASS, for complex the aditional copies make it for sure slower
uint32_t MatrixTransformBLAS(MatrixType * const M, MatrixType const A, MatrixType const T)
{
  CompactMatrixType AComp, TComp,TComp2,BComp;
  if ((A.Ni != A.Nj) || (A.Ni != T.Nj))
  {
    printf("Invallid dimensions for matrix transformation A.Ni should be T.Nj and A should be square\n");
    printf("A.Ni=%i A.Nj=%i T.Ni=%i T.Nj=%i \n", A.Ni, A.Nj, T.Ni, T.Nj);
    return 1;
  }
  
  printf("MatrixTransormBLAS\nA.Ni=%i A.Nj=%i T.Ni=%i T.Nj=%i \n", A.Ni, A.Nj, T.Ni, T.Nj);fflush(stdout);
  
  if(MatrixCopyToCompact(&AComp,A))
  {
    printf("Matrix conversion to compact failed in MatrixTransformBLAS\n");
    return 1;
  }
  if(MatrixCopyToCompact(&TComp,T))
  {
    printf("Matrix conversion to compact failed in MatrixTransformBLAS\n");
    return 1;
  }
  if(MatrixCopyToCompact(&TComp2,T))
  {
    printf("Matrix conversion to compact failed in MatrixTransformBLAS\n");
    return 1;
  }
  if(CompactMatrixMultiply3('C','N','N',&TComp,&AComp,&TComp2,&BComp))
  {
    printf("CompactMatrixMultiply3 failed in MatrixTransformBLAS\n");
    return 1;
  }
  if(CompactMatrixCopyToMatrix(M,BComp))
  {
    printf("CompactMatrixCopyToMatrix failed in MatrixTransformBLAS\n");
    return 1;
  }
  
  FreeCompactMatrix(&AComp);
  FreeCompactMatrix(&TComp);
  FreeCompactMatrix(&TComp2);
  FreeCompactMatrix(&BComp);
  
  return 0;
}
*/

/*  From BasicMath_MatrixTypeConversion. Only difference is it calls my MatrixTransformBLAS instead of MatrixTransform
 * 
 * 
 * Returns the Generalized Eigensystem of a non-orthnormal matrix.
 *  Given a basis psi_i
 *  Matrix M contains the values psi_i^* . H . psi_j
 *  Matrix S contains the values psi_i^* . psi_j
 *
 *  on exit val defines a diagonal matrix H and
 *  H = T . M . T^{\dag}
 *  I = T . S . T^{\dag}
 *  S = Sroot . Sroot^{\dag}
 *
 */
/* -- MatrixTransform changed to include BLAS - when it gives a gain. (T^* A T^T with A diagonal is faster without blass and even with A dense is not sure. For Real I can use BLASS, for complex the aditional copies make it for sure slower
uint32_t GeneralizedEigensystemSingularBLAS(MatrixType const M, MatrixType const S, MatrixType * const T, MatrixType * const Sroot, double *val)
{
  uint32_t i, j, NewDim, retval = 0;
  double Sval[S.Ni], bound;
  MatrixType EigenSysS, TMT, TST, Ttemp;
  MatrixCopy(&TST, S);
  MatrixScalarMultiply(&TST, -1, 0);
  Eigensystem(TST, &EigenSysS, Sval);
  MatrixScalarMultiply(&TST, 0, 0);
  MatrixTransformBLAS(&TMT, M, EigenSysS);
  NewDim = 1;
  bound = 0;
  for (i = 0; i < M.Ni; i++)
  {
    bound -= Sval[i];
  }
  bound *= DBL_EPSILON * 100;
  bound += DBL_MIN * 100;
  for (i = 0; i < M.Ni; i++)
  {
    if (-Sval[i] > bound)
    {
      TST.VRe[i][i] = -Sval[i];
      NewDim = i + 1;
    }
    else
    {
      for (j = 0; j < M.Nj; j++)
      {
        TMT.VRe[i][j] = 0;
        TMT.VRe[j][i] = 0;
      }
      TST.VRe[i][i] = 1;
    }
  }
  printf("Singularity border =%10.3E : Number of states kept = %i\n", bound, NewDim);
  printf("Norm of basis ( 0...7): ");
  for (i = 0; i < Min(8, S.Ni); i++)
  {
    printf("%10.3E ", -Sval[i]);
  }
  printf("..........\n");
  if (GeneralizedEigenvalues(TMT, TST, val))
  {
    printf("GeneralizedEigenvalues failed in GeneralizedEigensystemSingular\n");
    retval = 1;
  }
  for (i = NewDim; i < M.Ni; i++)
  {
    TMT.VRe[i][i] = val[M.Ni - 1] + 100 * DBL_EPSILON * Abs(val[M.Ni - 1]);
  }
  if (GeneralizedEigensystem(TMT, TST, &Ttemp, Sroot, val))
  {
    printf("GeneralizedEigensystem failed in GeneralizedEigensystemSingular\n");
    retval = 1;
  }
  for (i = NewDim; i < M.Ni; i++)
  {
    for (j = 0; j < M.Ni; j++)
    {
      Ttemp.VRe[j][i] = 0;
    }
  }
  MatrixMultiply(T, Ttemp, EigenSysS);
  FreeMatrix(&EigenSysS);
  FreeMatrix(&TMT);
  FreeMatrix(&TST);
  FreeMatrix(&Ttemp);
  return retval;
}
 */

/* From Lanczos.c  Only difference is it calls my GeneralizedEigensystemSingularBLAS instead of the existing version*/
/* -- MatrixTransform changed to include BLAS - when it gives a gain. (T^* A T^T with A diagonal is faster without blass and even with A dense is not sure. For Real I can use BLASS, for complex the aditional copies make it for sure slower
uint32_t LanczosTriDiagonalizeDysonLowMemBLAS(OperatorType const O, WaveFunctionType psi, double epsilon,  uint32_t const NKrylov, OperatorType const *V, TriDiagonalMatrixType const *g, uint32_t const NV, TriDiagonalMatrixType * Gout)
{
  uint32_t i, j;
  TriDiagonalMatrixType OKrylov;
  MatrixType Stemp, Htemp, HtempEigenFun, StempRoot;
  MatrixType Hval, *VMat;
  WaveFunctionType *PsiKrylov;
  OperatorType Otemp;
  TimeType Time;
  InitializeTiming(&Time);
  BeginTiming("Total", &Time);
  printf("\nStart of LanczosTriDiagonalizeDysonLowMem\n\n");
  // Initiallize matrices
  Htemp.Ni = NKrylov;
  Htemp.Nj = NKrylov;
  Htemp.QComplex = (O.QComplex || psi.QComplex);
  for (i = 0; i < NV; i++)
  {
    Htemp.QComplex = (Htemp.QComplex || V[i].QComplex);
  }
  if (InitMatrix(&Htemp))
  {
    printf("InitMatrix failed in LanczosTriDiagonalizeDyson\n");
    fflush(stdout);
    return 1;
  };
  Stemp.Ni = NKrylov;
  Stemp.Nj = NKrylov;
  Stemp.QComplex = Htemp.QComplex;
  if (InitMatrix(&Stemp))
  {
    printf("InitMatrix failed in LanczosTriDiagonalizeDyson\n");
    fflush(stdout);
    return 1;
  };
  Hval.Ni = NV;
  Hval.Nj = NKrylov;
  Hval.QComplex = 0;
  if (InitMatrix(&Hval))
  {
    printf("InitMatrix failed in LanczosTriDiagonalizeDyson\n");
    fflush(stdout);
    return 1;
  };
  if ((VMat = calloc(NV, sizeof(MatrixType))) == NULL)
  {
    printf("calloc failed in LanczosTriDiagonalizeDyson\n");
    fflush(stdout);
    return 1;
  };
  //Create Krylov basis of O
  OKrylov.N = NKrylov;
  InitTriDiagonalMatrix(&OKrylov);
  if ((PsiKrylov = calloc(NKrylov, sizeof(WaveFunctionType))) == NULL)
  {
    printf("calloc failed in LanczosTriDiagonalizeDyson\n");
    fflush(stdout);
    return 1;
  };
  BeginTiming("Krylov", &Time);
  LanczosTriDiagonalizeKrylov(O, psi, epsilon, &OKrylov, PsiKrylov);
  EndTiming("Krylov", &Time);
  // loop over all V[i]
  BeginTiming("Creation of matrix elements", &Time);
  for (i = 0; i < NV; i++)
  {
    printf("Create perturbation %i from %i\n", i + 1, NV);
    // create extended non orthonormal matrices on original basis
    // initiallize Otemp for V^{\dag} V
    InitOperatorForMult(&Otemp, 2, &V[i], &V[i]);
    // create Stemp
    OperatorMultiply2NormalOrder(&Otemp, V[i], V[i], 'D', 'N', &QDetNotOrderedFcFaBcBa);
    PsiOperatorPsiMatrix(Otemp, PsiKrylov, &Stemp);
    FreeOperator(&Otemp);
    // initiallize Otemp for V^{\dag} O V
    InitOperatorForMult(&Otemp, 3, &V[i], &O, &V[i]);
    // create Htemp
    OperatorMultiply3NormalOrder(&Otemp, V[i], O, V[i], 'D', 'N', 'N', &QDetNotOrderedFcFaBcBa);
    PsiOperatorPsiMatrix(Otemp, PsiKrylov, &Htemp);
    FreeOperator(&Otemp);
    // diagonalize Htemp with repsect to Stemp
    GeneralizedEigensystemSingularBLAS(Htemp, Stemp, &HtempEigenFun, &StempRoot, Hval.VRe[i]);
    MatrixMultiply(&(VMat[i]), HtempEigenFun, Stemp);
    // Free HtempEigenFun
    FreeMatrix(&HtempEigenFun);
  }
  // Free PsiKrylov
  for (j = 0; j < NKrylov; j++)
  {
    WaveFunctionFree(&(PsiKrylov[j]));
  }
  free(PsiKrylov);
  FreeMatrix(&Htemp);
  FreeMatrix(&Stemp);
  EndTiming("Creation of matrix elements", &Time);
  // create Gout
  BeginTiming("Final Tridiagonalize", &Time);
  if (TridiagonalizeImpurityCoupledToBathDyson(OKrylov, VMat, Hval, g, NV, Gout))
  {
    printf("TridiagonalizeImpurityCoupledToBathDyson failed in LanczosTriDiagonalizeDyson\n");
    FreeTriDiagonalMatrix(&OKrylov);
    FreeMatrix(&Hval);
    for (i = 0; i < NV; i++)
    {
      FreeMatrix(&(VMat[i]));
    }
    free(VMat);
    return 1;
  }
  EndTiming("Final Tridiagonalize", &Time);
  // Free things allocated
  FreeTriDiagonalMatrix(&OKrylov);
  FreeMatrix(&Hval);
  for (i = 0; i < NV; i++)
  {
    FreeMatrix(&(VMat[i]));
  }
  free(VMat);
  EndTiming("Total", &Time);
  PrintTiming(Time);
  FreeTiming(&Time);
  return 0;
}
*/