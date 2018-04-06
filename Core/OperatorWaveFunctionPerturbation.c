#include "OperatorWaveFunctionPerturbation.h"

/*
 *  The wave-funciton and operator is saved in parts.
 *  Operator O act on each part and does not connect different parts
 *  Operator V connects the parts of order m to the parts of order m+1.
 *  Operator V^{dag} connects the parts of order m to the parts of oder m-1
 */
uint32_t OperatorPsiPerturbationConserveBasis(const OperatorPerturbationType O, const WaveFunctionPerturbationType psi, WaveFunctionPerturbationType * const Opsi)
{
  uint32_t i, j;
  Opsi->NV = psi.NV;
  Opsi->order = psi.order;
  if (WaveFunctionPerturbationInit(Opsi))
  {
    printf("WaveFunctionPerturbationInit failed in OperatorPsiPerturbationConserveBasis\n");
    return 1;
  };
  for (i = 0; i < psi.order; i++)
  {
    for (j = 0; j < PerturbationOptions(psi.NV, i); j++)
    {
      OperatorPsiConserveBasis(O.O, psi.psiOrder[i][j], &(Opsi->psiOrder[i][j]));
    }
  }



  return 0;
}

uint32_t OperatorPsiPerturbationConserveBasisPreAllocated(const OperatorPerturbationType O, const WaveFunctionPerturbationType psi, WaveFunctionPerturbationType * const Opsi)
{
  printf("to do\n");
  return 0;
}

uint32_t OperatorPsiPerturbation(const OperatorPerturbationType O, const WaveFunctionPerturbationType psi, WaveFunctionPerturbationType * const Opsi)
{
  printf("to do\n");
  return 0;
}

uint32_t PsiOperatorPsiPerturbation(const OperatorPerturbationType O, const WaveFunctionPerturbationType psi, double *VRe, double *VeIm)
{
  printf("to do\n");
  return 0;
}

uint32_t PsiHermitianOperatorPsiPerturbation(const OperatorPerturbationType O, const WaveFunctionPerturbationType psi, double *VRe)
{
  printf("to do\n");
  return 0;
}

uint32_t PsiOperatorSqrPsiPerturbation(const OperatorPerturbationType O, const WaveFunctionPerturbationType psi, double *VRe, double *VeIm)
{
  printf("to do\n");
  return 0;
}

uint32_t PsiHermitianOperatorSqrPsiPerturbation(const OperatorPerturbationType O, const WaveFunctionPerturbationType psi, double *VRe)
{
  printf("to do\n");
  return 0;
}

uint32_t PsiPrimeOperatorPsiPerturbationMatrix(const OperatorPerturbationType O, const WaveFunctionPerturbationType * psiPrime, WaveFunctionPerturbationType * const psi, MatrixType * const V)
{
  printf("to do\n");
  return 0;
}

uint32_t PsiOperatorPsiPerturbationMatrix(const OperatorPerturbationType O, const WaveFunctionPerturbationType * psi, MatrixType * const V)
{
  printf("to do\n");
  return 0;
}

uint32_t PsiPrimeOperatorPsiPerturbation(const OperatorPerturbationType O, const WaveFunctionPerturbationType psiPrime, WaveFunctionPerturbationType const psi, double *VRe, double *VIm)
{
  printf("to do\n");
  return 0;
}
