#include "WaveFunctionPerturbationArithmetic.h"

/*  NV and order must be set.
 *  The operators of order order and perturbation NV are not initiallized
 */
uint32_t WaveFunctionPerturbationInit(WaveFunctionPerturbationType * Psi)
{
  uint32_t i;
  if ((Psi->psiOrder = (WaveFunctionType **)malloc(Psi->order * sizeof(WaveFunctionType *))) == NULL)
  {
    printf("malloc failed in WaveFunctionPerturbationInit\n");
    return 1;
  };
  for (i = 0; i < Psi->order; i++)
  {
    if ((Psi->psiOrder[i] = (WaveFunctionType*)malloc(PerturbationOptions(Psi->NV, i) * sizeof(WaveFunctionType))) == NULL)
    {
      printf("malloc failed in WaveFunctionPerturbationInit\n");
      return 1;
    };
  }
  return 0;
}

uint32_t WaveFunctionPerturbationInitAll(WaveFunctionPerturbationType * Psi, uint32_t NFermion, uint32_t NBoson,
                                         uint32_t QComplex, uint32_t NBitsKey)
{
  uint32_t i, j;
  if ((Psi->psiOrder = (WaveFunctionType **)malloc(Psi->order * sizeof(WaveFunctionType *))) == NULL)
  {
    printf("malloc failed in WaveFunctionPerturbationInit\n");
    return 1;
  };
  for (i = 0; i < Psi->order; i++)
  {
    if ((Psi->psiOrder[i] = (WaveFunctionType*)malloc(PerturbationOptions(Psi->NV, i) * sizeof(WaveFunctionType))) == NULL)
    {
      printf("malloc failed in WaveFunctionPerturbationInit\n");
      return 1;
    };
    for (j = 0; j < PerturbationOptions(Psi->NV, i); j++)
    {
      Psi->psiOrder[i][j].NFermion = NFermion;
      Psi->psiOrder[i][j].NBoson = NBoson;
      Psi->psiOrder[i][j].QComplex = QComplex;
      Psi->psiOrder[i][j].NBitsKey = NBitsKey;
      if (WaveFunctionInit(&(Psi->psiOrder[i][j])))
      {
        printf("WaveFunctionInit failed in WaveFunctionPerturbationInitAll\n");
        return 1;
      };
    }
  }
  return 0;
}

/*
 *
 */
uint32_t QWaveFunctionPerturbationMinimumBloksAlocated(WaveFunctionPerturbationType psi)
{
  int32_t i, j;
  for (i = 0; i < psi.order; i++)
  {
    for (j = 0; j < PerturbationOptions(psi.NV, i); j++)
    {
      if (psi.psiOrder[i][j].NBlok != 1)
        return 0;
    }
  }
  return 1;
}



uint32_t MakeComplexWaveFunctionPerturbationFromReal(WaveFunctionPerturbationType * Cpsi,
                                                     WaveFunctionPerturbationType psi)
{
  uint32_t i, j;
  Cpsi->NV = psi.NV;
  Cpsi->order = psi.order;
  WaveFunctionPerturbationInit(Cpsi);
  for (i = 0; i < psi.order; i++)
  {
    for (j = 0; j < PerturbationOptions(psi.NV, i); j++)
    {
      if (MakeComplexWaveFunctionFromReal(&(Cpsi->psiOrder[i][j]), psi.psiOrder[i][j]))
      {
        printf("MakeComplexWaveFunctionFromReal failed in MakeComplexWaveFunctionPerturbationFromReal\n");
        return 1;
      }
    }
  }
  return 0;
}

uint32_t MakeRealWaveFunctionPerturbationFromComplex(WaveFunctionPerturbationType Cpsi,
                                                     WaveFunctionPerturbationType * psi)
{
  uint32_t i, j;
  psi->NV = Cpsi.NV;
  psi->order = Cpsi.order;
  WaveFunctionPerturbationInit(psi);
  for (i = 0; i < Cpsi.order; i++)
  {
    for (j = 0; j < PerturbationOptions(Cpsi.NV, i); j++)
    {
      if (MakeRealWaveFunctionFromComplex(Cpsi.psiOrder[i][j], &(psi->psiOrder[i][j])))
      {
        printf("MakeRealWaveFunctionFromComplex failed in MakeRealWaveFunctionPerturbationFromComplex\n");
        return 1;
      }
    }
  }
  return 0;
}

void WaveFunctionPerturbationFreeFramework(WaveFunctionPerturbationType * const psi)
{
  uint32_t i;
  for (i = 0; i < psi->order; i++)
  {
    free(psi->psiOrder[i]);
  }
  free(psi->psiOrder);
}

void WaveFunctionPerturbationFree(WaveFunctionPerturbationType * const psi)
{
  uint32_t i, j;
  for (i = 0; i < psi->order; i++)
  {
    for (j = 0; j < PerturbationOptions(psi->NV, i); j++)
    {
      WaveFunctionFree(&(psi->psiOrder[i][j]));
    }
    free(psi->psiOrder[i]);
  }
  free(psi->psiOrder);
}

void WaveFunctionPerturbationFreeValues(WaveFunctionPerturbationType * const psi)
{
  uint32_t i, j;
  for (i = 0; i < psi->order; i++)
  {
    for (j = 0; j < PerturbationOptions(psi->NV, i); j++)
    {
      WaveFunctionFreeValues(&(psi->psiOrder[i][j]));
    }
  }
}

uint32_t WaveFunctionPerturbationAddElement(WaveFunctionPerturbationType * const psi, const uint8_t * const det,
                                            const uint32_t order, const uint32_t iV, const double aRe, const double aIm)
{
  return WaveFunctionAddElement(&(psi->psiOrder[order][iV]), det, aRe, aIm);
}

uint32_t WaveFunctionPerturbationAddElementConserveBasis(WaveFunctionPerturbationType * const psi,
                                                         const uint8_t * const det, const uint32_t order,
                                                         const uint32_t iV, const double aRe, const double aIm)
{
  return WaveFunctionAddElementConserveBasis(&(psi->psiOrder[order][iV]), det, aRe, aIm);
}

void WaveFunctionPerturbationPrefactorDet(const WaveFunctionPerturbationType psi, const uint8_t * const det,
                                          const uint32_t order, const uint32_t iV, double *aRe, double *aIm)
{
  WaveFunctionPrefactorDet(psi.psiOrder[order][iV], det, aRe, aIm);
}

/*
 *  a <- a + alpha b
 */
uint32_t WaveFunctionPerturbationAddSameBasis(WaveFunctionPerturbationType * const a,
                                              const WaveFunctionPerturbationType b)
{
  printf("To do\n");
  return 0;
}

/*
 *  a <- a + alpha b
 */
uint32_t WaveFunctionPerturbationAddMultipliedSameBasis(WaveFunctionPerturbationType * const a,
                                                        const WaveFunctionPerturbationType b, const double alpha)
{
  printf("To do\n");
  return 0;
}

/*
 *  a <- a + alpha b
 */
uint32_t WaveFunctionPerturbationAddComplexMultipliedSameBasis(WaveFunctionPerturbationType * const a,
                                                               const WaveFunctionPerturbationType b,
                                                               const double alphaRe, const double alphaIm)
{
  printf("To do\n");
  return 0;
}


/*
 *  a <- a + alpha b
 */
uint32_t WaveFunctionPerturbationAddConserveBasis(WaveFunctionPerturbationType * const a,
                                                  const WaveFunctionPerturbationType b)
{
  printf("To do\n");
  return 0;
}

/*
 *  a <- a + alpha b
 */
uint32_t WaveFunctionPerturbationAddMultipliedConserveBasis(WaveFunctionPerturbationType * const a,
                                                            const WaveFunctionPerturbationType b, const double alpha)
{
  printf("To do\n");
  return 0;
}

/*
 *  a <- a + alpha b
 */
uint32_t WaveFunctionPerturbationAddComplexMultipliedConserveBasis(WaveFunctionPerturbationType * const a,
                                                                   const WaveFunctionPerturbationType b,
                                                                   const double alphaRe, const double alphaIm)
{
  printf("To do\n");
  return 0;
}

/*
 *  a <- a + alpha b
 */
uint32_t WaveFunctionPerturbationAddUnionBasis(WaveFunctionPerturbationType * const a,
                                               const WaveFunctionPerturbationType b)
{
  printf("To do\n");
  return 0;
}

/*
 *  a <- a + alpha b
 */
uint32_t WaveFunctionPerturbationAddMultipliedUnionBasis(WaveFunctionPerturbationType * const a,
                                                         const WaveFunctionPerturbationType b, const double alpha)
{
  printf("To do\n");
  return 0;
}

/*
 *  a <- a + alpha b
 */
uint32_t WaveFunctionPerturbationAddComplexMultipliedUnionBasis(WaveFunctionPerturbationType * const a,
                                                                const WaveFunctionPerturbationType b,
                                                                const double alphaRe, const double alphaIm)
{
  printf("To do\n");
  return 0;
}

/*
 *  a -> alpha a
 */
uint32_t WaveFunctionPerturbationMultiply(WaveFunctionPerturbationType a, const double bRe, const double bIm)
{
  uint32_t i, j;
  for (i = 0; i < a.order; i++)
  {
    for (j = 0; j < PerturbationOptions(a.NV, i); j++)
    {
      if (WaveFunctionMultiply(&a.psiOrder[i][j], bRe, bIm))
      {
        printf("WaveFunctionMultiply failed in WaveFunctionPerturbationMultiply\n");
        return 1;
      }
    }
  }
  return 0;
}

/*
 *  <a|a>
 */
void WaveFunctionPerturbationCalculateNorm(const WaveFunctionPerturbationType a, double *Norm)
{
  uint32_t i, j;
  double Normpart;
  *Norm = 0;
  for (i = 0; i < a.order; i++)
  {
    for (j = 0; j < PerturbationOptions(a.NV, i); j++)
    {
      WaveFunctionCalculateNorm(a.psiOrder[i][j], &Normpart);
      *Norm += Normpart;
    }
  }
}

/*
 *  a -> a / Sqrt[<a|a>]
 */
void WaveFunctionPerturbationNormalize(WaveFunctionPerturbationType * const a)
{
  double Norm;
  WaveFunctionPerturbationCalculateNorm(*a, &Norm);
  Norm = sqrt(Norm);
  WaveFunctionPerturbationMultiply(*a, Norm, 0);
}



/*
 *   calculate <a|b>
 */
void WaveFunctionPerturbationDotProductSameBasis(const WaveFunctionPerturbationType a,
                                                 const WaveFunctionPerturbationType b, double *DotRe, double *DotIm)
{
  double DotRePart, DotImPart;
  uint32_t i, j;
  *DotRe = 0;
  *DotIm = 0;
  for (i = 0; i < a.order; i++)
  {
    for (j = 0; j < PerturbationOptions(a.NV, i); j++)
    {
      WaveFunctionDotProductSameBasis(a.psiOrder[i][j], b.psiOrder[i][j], &DotRePart, &DotImPart);
      *DotRe += DotRePart;
      *DotIm += DotImPart;
    }
  }
}

/*
 *   calculate <a|b>
 */
void WaveFunctionPerturbationDotProduct(const WaveFunctionPerturbationType a, const WaveFunctionPerturbationType b,
                                        double *DotRe, double *DotIm)
{
  double DotRePart, DotImPart;
  uint32_t i, j;
  *DotRe = 0;
  *DotIm = 0;
  for (i = 0; i < a.order; i++)
  {
    for (j = 0; j < PerturbationOptions(a.NV, i); j++)
    {
      WaveFunctionDotProduct(a.psiOrder[i][j], b.psiOrder[i][j], &DotRePart, &DotImPart);
      *DotRe += DotRePart;
      *DotIm += DotImPart;
    }
  }
}

/*
 *  a -> random vector of norm 1 on basis of a
 */
void WaveFunctionPerturbationRandomize(WaveFunctionPerturbationType * const a)
{
  uint32_t i, j;
  for (i = 0; i < a->order; i++)
  {
    for (j = 0; j < PerturbationOptions(a->NV, i); j++)
    {
      WaveFunctionRandomize(&(a->psiOrder[i][j]));
    }
  }
  WaveFunctionPerturbationNormalize(a);
}

/*
 *  a -> a with all basis functions with an absolute prefactor smaller then epsilon removed
 */
void WaveFunctionPerturbationTruncateBasis(WaveFunctionPerturbationType * const a, const double epsilon)
{
  uint32_t i, j;
  for (i = 0; i < a->order; i++)
  {
    for (j = 0; j < PerturbationOptions(a->NV, i); j++)
    {
      WaveFunctionTruncateBasis(&(a->psiOrder[i][j]), epsilon);
    }
  }
}

/*
 *  Make psi orthogonal to all OtherPsi, but conserve the basis set of psi
 */
void WaveFunctionPerturbationOrthogonalizeConserveBasis(WaveFunctionPerturbationType * const psi,
                                                        const WaveFunctionPerturbationType * const OtherPsi,
                                                        const uint32_t Npsi)
{
  double DotRe, DotIm;
  uint32_t i;
  for (i = 0; i < Npsi; i++)
  {
    WaveFunctionPerturbationDotProduct(OtherPsi[i], *psi, &DotRe, &DotIm);
    WaveFunctionPerturbationAddComplexMultipliedConserveBasis(psi, OtherPsi[i], -DotRe, -DotIm);
  }
}

/*
 *  Make psi orthogonal to all OtherPsi
 */
uint32_t WaveFunctionPerturbationOrthogonalize(WaveFunctionPerturbationType * const psi,
                                               const WaveFunctionPerturbationType * const OtherPsi, const uint32_t Npsi)
{
  double DotRe, DotIm;
  uint32_t i;
  for (i = 0; i < Npsi; i++)
  {
    WaveFunctionPerturbationDotProduct(*OtherPsi, *psi, &DotRe, &DotIm);
    if (WaveFunctionPerturbationAddComplexMultipliedUnionBasis(psi, OtherPsi[i], -DotRe, -DotIm))
    {
      printf
        ("WaveFunctionPerturbationAddComplexMultipliedUnionBasis failed in WaveFunctionPerturbationOrthogonalize\n");
      return 1;
    };
  }
  return 0;
}


uint32_t CreateStartWaveFunctionPerturbationOnsiteEnergy(const OperatorPerturbationType O,
                                                         WaveFunctionPerturbationType * psi, uint32_t order,
                                                         uint32_t nFermion, uint32_t nBoson, uint32_t NBitsKey)
{
  psi->NV = O.NV;
  psi->order = order;
  if (WaveFunctionPerturbationInitAll(psi, nFermion, nBoson, 0, NBitsKey))
  {
    printf("WaveFunctionPerturbationInitAll failed in CreateStartWaveFunctionPerturbationOnsiteEnergy\n");
    return 1;
  };
  WaveFunctionFree(&(psi->psiOrder[0][0]));
  if (CreateStartWaveFunctionOnsiteEnergy(O.O, &(psi->psiOrder[0][0]), nFermion, nBoson, NBitsKey))
  {
    printf("CreateStartWaveFunctionOnsiteEnergy failed in CreateStartWaveFunctionPerturbationOnsiteEnergy\n");
    return 1;
  };
  return 0;
}


/*  increase the number of fermions in psiIn by one, creating psiOut, the Fermion is added at position FermionPos
 *  If occupied is one the Fermion is occupied
 *  If occupied is zero the Fermion is empty
 */
uint32_t WaveFunctionPerturbationAddFermion(WaveFunctionPerturbationType const psiIn,
                                            WaveFunctionPerturbationType * const psiOut, uint32_t FermionPos,
                                            uint32_t occupied)
{
  uint32_t i, j;
  psiOut->NV = psiIn.NV;
  psiOut->order = psiIn.order;
  if (WaveFunctionPerturbationInit(psiOut))
  {
    printf("WaveFunctionPerturbationInit failed in WaveFunctionPerturbationAddFermion\n");
    return 1;
  };
  for (i = 0; i < psiIn.order; i++)
  {
    for (j = 0; j < PerturbationOptions(psiIn.NV, i); j++)
    {
      if (WaveFunctionAddFermion(psiIn.psiOrder[i][j], &(psiOut->psiOrder[i][j]), FermionPos, occupied))
      {
        printf("WaveFunctionAddFermion failed in WaveFunctionPerturbationAddFermion\n");
        return 1;
      };
    }
  }
  return 0;
}

uint32_t WaveFunctionPerturbationRemoveFermion(WaveFunctionPerturbationType const psiIn,
                                               WaveFunctionPerturbationType * const psiOut, uint32_t FermionPos)
{
  uint32_t i, j;
  psiOut->NV = psiIn.NV;
  psiOut->order = psiIn.order;
  if (WaveFunctionPerturbationInit(psiOut))
  {
    printf("WaveFunctionPerturbationInit failed in WaveFunctionPerturbationRemoveFermion\n");
    return 1;
  };
  for (i = 0; i < psiIn.order; i++)
  {
    for (j = 0; j < PerturbationOptions(psiIn.NV, i); j++)
    {
      if (WaveFunctionRemoveFermion(psiIn.psiOrder[i][j], &(psiOut->psiOrder[i][j]), FermionPos))
      {
        printf("WaveFunctionRemoveFermion failed in WaveFunctionPerturbationRemoveFermion\n");
        return 1;
      };
    }
  }
  return 0;
}
