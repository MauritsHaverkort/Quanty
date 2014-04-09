#include "RealWaveFunctionArithmetic.c"
#include "ComplexWaveFunctionArithmetic.c"
#include "RealComplexWaveFunctionArithmetic.c"
#include "ComplexRealWaveFunctionArithmetic.c"


/*  NFermion, NBoson, QComplex, NBitsKey need to be set
 *
 */
uint32_t WaveFunctionInit(WaveFunctionType * Psi)
{
  strcpy(Psi->Name, "Wave Function");
  Psi->N = 0;
  Psi->NBlok = 1;
  Psi->NWord = Psi->NBoson + Min(Psi->NFermion, 1 + (Psi->NFermion - 1) / 8);
  if ((Psi->KeyList = calloc((uint32_t) 1 << (Psi->NBitsKey), sizeof(uint32_t))) == NULL)
  {
    printf("calloc failed in WaveFunctionInit 00\n");
    fflush(stdout);
    Psi->NBlok = 0;
    Psi->NWord = 0;
    return 1;
  };
  if ((Psi->ValueRe = malloc(Psi->NBlok * sizeof(double *))) == NULL)
  {
    printf("malloc failed in WaveFunctionInit 01\n");
    fflush(stdout);
    free(Psi->KeyList);
    Psi->NBlok = 0;
    Psi->NWord = 0;
    return 1;
  };
  if (Psi->QComplex)
    if ((Psi->ValueIm = malloc(Psi->NBlok * sizeof(double *))) == NULL)
    {
      printf("malloc failed in WaveFunctionInit 02\n");
      fflush(stdout);
      free(Psi->KeyList);
      free(Psi->ValueRe);
      Psi->NBlok = 0;
      Psi->NWord = 0;
      return 1;
    };
  if ((Psi->Next = malloc(Psi->NBlok * sizeof(uint32_t *))) == NULL)
  {
    printf("malloc failed in WaveFunctionInit 03\n");
    fflush(stdout);
    free(Psi->KeyList);
    free(Psi->ValueRe);
    if (Psi->QComplex)
      free(Psi->ValueIm);
    Psi->NBlok = 0;
    Psi->NWord = 0;
    return 1;
  };
  if ((Psi->Determinant = malloc(Psi->NBlok * sizeof(uint8_t *))) == NULL)
  {
    printf("malloc failed in WaveFunctionInit 04\n");
    fflush(stdout);
    free(Psi->KeyList);
    free(Psi->ValueRe);
    if (Psi->QComplex)
      free(Psi->ValueIm);
    free(Psi->Next);
    Psi->NBlok = 0;
    Psi->NWord = 0;
    return 1;
  };
  if ((Psi->ValueRe[Psi->NBlok - 1] = malloc(BlokSize * sizeof(double))) == NULL)
  {
    printf("malloc failed in WaveFunctionInit 05\n");
    fflush(stdout);
    free(Psi->KeyList);
    free(Psi->ValueRe);
    if (Psi->QComplex)
      free(Psi->ValueIm);
    free(Psi->Next);
    free(Psi->Determinant);
    Psi->NBlok = 0;
    Psi->NWord = 0;
    return 1;
  };
  Psi->ValueRe[0][0] = 0.0;
  if (Psi->QComplex)
    if ((Psi->ValueIm[Psi->NBlok - 1] = malloc(BlokSize * sizeof(double))) == NULL)
    {
      printf("malloc failed in WaveFunctionInit 06\n");
      fflush(stdout);
      free(Psi->KeyList);
      free(Psi->ValueRe);
      if (Psi->QComplex)
        free(Psi->ValueIm);
      free(Psi->Next);
      free(Psi->Determinant);
      free(Psi->ValueRe[Psi->NBlok - 1]);
      Psi->NBlok = 0;
      Psi->NWord = 0;
      return 1;
    };
  if (Psi->QComplex)
    Psi->ValueIm[0][0] = 0.0;
  if ((Psi->Next[Psi->NBlok - 1] = malloc(BlokSize * sizeof(uint32_t))) == NULL)
  {
    printf("malloc failed in WaveFunctionInit 07\n");
    fflush(stdout);
    free(Psi->KeyList);
    free(Psi->ValueRe);
    if (Psi->QComplex)
      free(Psi->ValueIm);
    free(Psi->Next);
    free(Psi->Determinant);
    free(Psi->ValueRe[Psi->NBlok - 1]);
    if (Psi->QComplex)
      free(Psi->ValueIm[Psi->NBlok - 1]);
    Psi->NBlok = 0;
    Psi->NWord = 0;
    return 1;
  };
  if ((Psi->Determinant[Psi->NBlok - 1] = malloc(BlokSize * Psi->NWord * sizeof(uint8_t))) == NULL)
  {
    printf("malloc failed in WaveFunctionInit 08\n");
    fflush(stdout);
    free(Psi->KeyList);
    free(Psi->ValueRe);
    if (Psi->QComplex)
      free(Psi->ValueIm);
    free(Psi->Next);
    free(Psi->Determinant);
    free(Psi->ValueRe[Psi->NBlok - 1]);
    if (Psi->QComplex)
      free(Psi->ValueIm[Psi->NBlok - 1]);
    free(Psi->Next[Psi->NBlok - 1]);
    Psi->NBlok = 0;
    Psi->NWord = 0;
    return 1;
  };
  return 0;
}

/* Add the complex values (and allocs)
 *  
 */
uint32_t MakeComplexWaveFunctionFromReal(WaveFunctionType * Cpsi, WaveFunctionType psi)
{
  uint32_t i, j;
  Cpsi->N = psi.N;
  Cpsi->NBlok = psi.NBlok;
  Cpsi->NFermion = psi.NFermion;
  Cpsi->NBoson = psi.NBoson;
  Cpsi->NWord = psi.NWord;
  Cpsi->ValueRe = psi.ValueRe;
  Cpsi->Next = psi.Next;
  Cpsi->Determinant = psi.Determinant;
  Cpsi->NBitsKey = psi.NBitsKey;
  Cpsi->KeyList = psi.KeyList;
  Cpsi->ValueIm = malloc(Cpsi->NBlok * sizeof(double *));
  if (Cpsi->ValueIm == NULL)
  {
    printf("could not allocate memory in MakeComplexWaveFunctionFromReal 00\n");
    return 1;
  }
  for (i = 0; i < Cpsi->NBlok; i++)
  {
    Cpsi->ValueIm[i] = calloc(BlokSize, sizeof(double));
    if (Cpsi->ValueIm[i] == NULL)
    {
      for (j = 0; j < i; j++)
      {
        free(Cpsi->ValueIm[j]);
      }
      printf("could not allocate memory in MakeComplexWaveFunctionFromReal %i\n", i);
      return 1;
    }
  }
  Cpsi->QComplex = 1;
  return 0;
}

/* free the complex values, keep the real
 *
 */
uint32_t MakeRealWaveFunctionFromComplex(WaveFunctionType Cpsi, WaveFunctionType * psi)
{
  uint32_t i;
  psi->N = Cpsi.N;
  psi->NBlok = Cpsi.NBlok;
  psi->NFermion = Cpsi.NFermion;
  psi->NBoson = Cpsi.NBoson;
  psi->NWord = Cpsi.NWord;
  psi->ValueRe = Cpsi.ValueRe;
  psi->ValueIm = Cpsi.ValueIm;
  psi->Next = Cpsi.Next;
  psi->Determinant = Cpsi.Determinant;
  psi->NBitsKey = Cpsi.NBitsKey;
  psi->KeyList = Cpsi.KeyList;
  for (i = 0; i < psi->NBlok; i++)
  {
    free(psi->ValueIm[i]);
  }
  free(psi->ValueIm);
  psi->QComplex = 0;
  return 0;
}

void WaveFunctionFree(WaveFunctionType * const psi)
{
  if (psi->QComplex)
  {
    ComplexWaveFunctionFree(psi);
  }
  else
  {
    RealWaveFunctionFree(psi);
  }
}

void WaveFunctionFreeValues(WaveFunctionType * const psi)
{
  if (psi->QComplex)
  {
    ComplexWaveFunctionFreeValues(psi);
  }
  else
  {
    RealWaveFunctionFreeValues(psi);
  }
}

uint32_t WaveFunctionAddElement(WaveFunctionType * const psi, const uint8_t * const det, const double aRe,
                                const double aIm)
{
  if (psi->QComplex == 0)
  {
    if (fabs(aIm) < DBL_MIN * 10)
    {
      return RealWaveFunctionAddElement(psi, det, aRe);
    }
    if (MakeComplexWaveFunctionFromReal(psi, *psi))
    {
      printf("MakeComplexWaveFunctionFromReal failed in WaveFunctionAddElement\n");
      fflush(stdout);
      return 1;
    }
  }
  return ComplexWaveFunctionAddElement(psi, det, aRe, aIm);
}

uint32_t WaveFunctionAddElementConserveBasis(WaveFunctionType * const psi, const uint8_t * const det, const double aRe,
                                             const double aIm)
{
  if (psi->QComplex == 0)
  {
    if (fabs(aIm) < DBL_MIN * 10)
    {
      RealWaveFunctionAddElementConserveBasis(psi, det, aRe);
      return 0;
    }
    if (MakeComplexWaveFunctionFromReal(psi, *psi))
    {
      printf("MakeComplexWaveFunctionFromReal failed in WaveFunctionAddElement\n");
      fflush(stdout);
      return 1;
    };
  }
  ComplexWaveFunctionAddElementConserveBasis(*psi, det, aRe, aIm);
  return 0;
}

void WaveFunctionPrefactorDet(const WaveFunctionType psi, const uint8_t * const det, double *aRe, double *aIm)
{
  if (psi.QComplex == 0)
  {
    *aIm = 0;
    *aRe = RealWaveFunctionPrefactorDet(psi, det);
  }
  else
  {
    ComplexWaveFunctionPrefactorDet(psi, det, aRe, aIm);
  }
}

uint32_t WaveFunctionAddSameBasis(WaveFunctionType a, const WaveFunctionType b)
{
  if ((b.QComplex) && (!a.QComplex))
    if (MakeComplexWaveFunctionFromReal(&a, a))
    {
      printf("MakeComplexWaveFunctionFromReal failed in WaveFunctionAddSameBasis\n");
      return 1;
    }
  if (b.QComplex)
  {
    ComplexWaveFunctionAddSameBasis(a, b);
  }
  else
  {
    RealWaveFunctionAddSameBasis(a, b);
  }
  return 0;
}

uint32_t WaveFunctionAddMultipliedSameBasis(WaveFunctionType a, const WaveFunctionType b, const double alpha)
{
  if ((b.QComplex) && (!a.QComplex))
    if (MakeComplexWaveFunctionFromReal(&a, a))
    {
      printf("MakeComplexWaveFunctionFromReal failed in WaveFunctionAddMultipliedSameBasis\n");
      return 1;
    }
  if (b.QComplex)
  {
    ComplexWaveFunctionAddMultipliedSameBasis(a, b, alpha);
  }
  else
  {
    RealWaveFunctionAddMultipliedSameBasis(a, b, alpha);
  }
  return 0;
}

uint32_t WaveFunctionAddComplexMultipliedSameBasis(WaveFunctionType a, const WaveFunctionType b, const double alphaRe,
                                                   const double alphaIm)
{
  if ((b.QComplex) && (!a.QComplex))
    if (MakeComplexWaveFunctionFromReal(&a, a))
    {
      printf("MakeComplexWaveFunctionFromReal failed in WaveFunctionAddComplexMultipliedSameBasis\n");
      return 1;
    }
  if ((Abs(alphaIm) > (DBL_MIN + Abs(alphaRe) * DBL_EPSILON)) && (!a.QComplex))
    if (MakeComplexWaveFunctionFromReal(&a, a))
    {
      printf("MakeComplexWaveFunctionFromReal failed in WaveFunctionAddComplexMultipliedSameBasis\n");
      return 1;
    }
  if (b.QComplex)
  {
    ComplexWaveFunctionAddComplexMultipliedSameBasis(a, b, alphaRe, alphaIm);
  }
  else if (a.QComplex)
  {
    ComplexRealWaveFunctionAddComplexMultipliedSameBasis(a, b, alphaRe, alphaIm);
  }
  else
  {
    RealWaveFunctionAddMultipliedSameBasis(a, b, alphaRe);
  }
  return 0;
}

uint32_t WaveFunctionAddConserveBasis(WaveFunctionType * const a, const WaveFunctionType b)
{
  if ((b.QComplex) && (!a->QComplex))
    if (MakeComplexWaveFunctionFromReal(a, *a))
    {
      printf("MakeComplexWaveFunctionFromReal failed in WaveFunctionAddConserveBasis\n");
      return 1;
    }
  if (b.QComplex)
  {
    ComplexWaveFunctionAddConserveBasis(a, b);
  }
  else
  {
    RealWaveFunctionAddConserveBasis(a, b);
  }
  return 0;
}

uint32_t WaveFunctionAddMultipliedConserveBasis(WaveFunctionType * const a, const WaveFunctionType b,
                                                const double alpha)
{
  if ((b.QComplex) && (!a->QComplex))
    if (MakeComplexWaveFunctionFromReal(a, *a))
    {
      printf("MakeComplexWaveFunctionFromReal failed in WaveFunctionAddMultipliedConserveBasis\n");
      return 1;
    }
  if (b.QComplex)
  {
    ComplexWaveFunctionAddMultipliedConserveBasis(a, b, alpha);
  }
  else
  {
    RealWaveFunctionAddMultipliedConserveBasis(a, b, alpha);
  }
  return 0;
}

uint32_t WaveFunctionAddComplexMultipliedConserveBasis(WaveFunctionType * const a, const WaveFunctionType b,
                                                       const double alphaRe, const double alphaIm)
{
  if ((b.QComplex) && (!a->QComplex))
    if (MakeComplexWaveFunctionFromReal(a, *a))
    {
      printf("MakeComplexWaveFunctionFromReal failed in WaveFunctionAddComplexMultipliedConserveBasis\n");
      return 1;
    }
  if ((Abs(alphaIm) > (DBL_MIN + Abs(alphaRe) * DBL_EPSILON)) && (!a->QComplex))
    if (MakeComplexWaveFunctionFromReal(a, *a))
    {
      printf("MakeComplexWaveFunctionFromReal failed in WaveFunctionAddComplexMultipliedConserveBasis\n");
      return 1;
    }
  if (b.QComplex)
  {
    ComplexWaveFunctionAddComplexMultipliedConserveBasis(a, b, alphaRe, alphaIm);
  }
  else if (a->QComplex)
  {
    ComplexRealWaveFunctionAddComplexMultipliedConserveBasis(a, b, alphaRe, alphaIm);
  }
  else
  {
    RealWaveFunctionAddMultipliedConserveBasis(a, b, alphaRe);
  }
  return 0;
}

uint32_t WaveFunctionAddUnionBasis(WaveFunctionType * const a, const WaveFunctionType b)
{
  if ((b.QComplex) && (!a->QComplex))
    if (MakeComplexWaveFunctionFromReal(a, *a))
    {
      printf("MakeComplexWaveFunctionFromReal failed in WaveFunctionAddUnionBasis\n");
      return 1;
    }
  if (b.QComplex)
  {
    return ComplexWaveFunctionAddUnionBasis(a, b);
  }
  else
  {
    return RealWaveFunctionAddUnionBasis(a, b);
  }
}

uint32_t WaveFunctionAddMultipliedUnionBasis(WaveFunctionType * const a, const WaveFunctionType b, const double alpha)
{
  if ((b.QComplex) && (!a->QComplex))
    if (MakeComplexWaveFunctionFromReal(a, *a))
    {
      printf("MakeComplexWaveFunctionFromReal failed in WaveFunctionAddMultipliedUnionBasis\n");
      return 1;
    }
  if (b.QComplex)
  {
    return ComplexWaveFunctionAddMultipliedUnionBasis(a, b, alpha);
  }
  else
  {
    return RealWaveFunctionAddMultipliedUnionBasis(a, b, alpha);
  }
}

uint32_t WaveFunctionAddComplexMultipliedUnionBasis(WaveFunctionType * const a, const WaveFunctionType b,
                                                    const double alphaRe, const double alphaIm)
{
  if ((b.QComplex) && (!a->QComplex))
    if (MakeComplexWaveFunctionFromReal(a, *a))
    {
      printf("MakeComplexWaveFunctionFromReal failed in WaveFunctionAddComplexMultipliedUnionBasis\n");
      return 1;
    }
  if ((Abs(alphaIm) > (DBL_MIN + Abs(alphaRe) * DBL_EPSILON)) && (!a->QComplex))
    if (MakeComplexWaveFunctionFromReal(a, *a))
    {
      printf("MakeComplexWaveFunctionFromReal failed in WaveFunctionAddComplexMultipliedUnionBasis\n");
      return 1;
    }
  if (b.QComplex)
  {
    return ComplexWaveFunctionAddComplexMultipliedUnionBasis(a, b, alphaRe, alphaIm);
  }
  else if (a->QComplex)
  {
    return ComplexRealWaveFunctionAddComplexMultipliedUnionBasis(a, b, alphaRe, alphaIm);
  }
  else
  {
    return RealWaveFunctionAddMultipliedUnionBasis(a, b, alphaRe);
  }
}

uint32_t WaveFunctionMultiply(WaveFunctionType a, const double bRe, const double bIm)
{
  if (a.QComplex == 0)
  {
    if (fabs(bIm) < DBL_MIN * 10)
    {
      RealWaveFunctionMultiply(a, bRe);
      return 0;
    }
    if (MakeComplexWaveFunctionFromReal(&a, a))
    {
      printf("MakeComplexWaveFunctionFromReal failed in WaveFunctionMultiply\n");
      fflush(stdout);
      return 1;
    };
  }
  ComplexWaveFunctionComplexMultiply(a, bRe, bIm);
  return 0;
}

void WaveFunctionCalculateNorm(const WaveFunctionType a, double *Norm)
{
  if (a.QComplex == 0)
  {
    *Norm = RealWaveFunctionCalculateNorm(a);
  }
  else
  {
    *Norm = ComplexWaveFunctionCalculateNorm(a);
  }
}

void WaveFunctionNormalize(WaveFunctionType * const a)
{
  if (a->QComplex == 0)
  {
    RealWaveFunctionNormalize(a);
  }
  else
  {
    ComplexWaveFunctionNormalize(a);
  }
}




/*
 *   calculate <a|b>
 */
void WaveFunctionDotProductSameBasis(const WaveFunctionType a, const WaveFunctionType b, double *DotRe, double *DotIm)
{
  if (a.QComplex)
  {
    if (b.QComplex)
    {
      ComplexWaveFunctionDotProductSameBasis(a, b, DotRe, DotIm);
    }
    else
    {
      ComplexRealWaveFunctionDotProductSameBasis(a, b, DotRe, DotIm);
    }
  }
  else
  {
    if (b.QComplex)
    {
      RealComplexWaveFunctionDotProductSameBasis(a, b, DotRe, DotIm);
    }
    else
    {
      *DotIm = 0;
      *DotRe = RealWaveFunctionDotProductSameBasis(a, b);
    }
  }
}

/*
*   calculate <a|b>
*/
void WaveFunctionDotProduct(const WaveFunctionType a, const WaveFunctionType b, double *DotRe, double *DotIm)
{
  if (a.QComplex)
  {
    if (b.QComplex)
    {
      ComplexWaveFunctionDotProduct(a, b, DotRe, DotIm);
    }
    else
    {
      ComplexRealWaveFunctionDotProduct(a, b, DotRe, DotIm);
    }
  }
  else
  {
    if (b.QComplex)
    {
      RealComplexWaveFunctionDotProduct(a, b, DotRe, DotIm);
    }
    else
    {
      *DotIm = 0;
      *DotRe = RealWaveFunctionDotProduct(a, b);
    }
  }
}

void WaveFunctionRandomize(WaveFunctionType * const a)
{
  if (a->QComplex == 0)
  {
    RealWaveFunctionRandomize(a);
  }
  else
  {
    ComplexWaveFunctionRandomize(a);
  }
}

void WaveFunctionTruncateBasis(WaveFunctionType * const a, const double epsilon)
{
  if (a->QComplex == 0)
  {
    RealWaveFunctionTruncateBasis(a, epsilon);
  }
  else
  {
    ComplexWaveFunctionTruncateBasis(a, epsilon);
  }
}

void WaveFunctionOrthogonalizeConserveBasis(WaveFunctionType * const psi, const WaveFunctionType * const OtherPsi,
                                            const uint32_t Npsi)
{
  if (psi->QComplex == 0)
  {
    RealWaveFunctionOrthogonalizeConserveBasis(psi, OtherPsi, Npsi);
  }
  else
  {
    ComplexWaveFunctionOrthogonalizeConserveBasis(psi, OtherPsi, Npsi);
  }
}

uint32_t WaveFunctionOrthogonalize(WaveFunctionType * const psi, const WaveFunctionType * const OtherPsi,
                                   const uint32_t Npsi)
{
  if (psi->QComplex == 0)
  {
    return RealWaveFunctionOrthogonalize(psi, OtherPsi, Npsi);
  }
  else
  {
    return ComplexWaveFunctionOrthogonalize(psi, OtherPsi, Npsi);
  }
}



uint32_t CreateStartWaveFunctionOnsiteEnergy(const OperatorType O, WaveFunctionType * psi, uint32_t nFermion,
                                             uint32_t nBoson, uint32_t NBitsKey)
{
  uint16_t i, j, jpick;
  double MinEnergy, Energy, EnergySqr;
  uint8_t Determinant[O.NWord];
  psi->NFermion = O.NFermion;
  psi->NBoson = O.NBoson;
  psi->QComplex = 0;
  psi->NBitsKey = NBitsKey;
  WaveFunctionInit(psi);
  if (nFermion > O.NFermion)
  {
    printf("Can't create a strating wave-function with %i Fermions, when there are only %i basis fermions\n", nFermion,
           O.NFermion);
    return 1;
  };
  if (nBoson > 256 * O.NBoson)
  {
    printf("Can't create a strating wave-function with %i Bosons, when there are only %i basis bosons\n", nBoson,
           O.NBoson);
    return 1;
  };
  for (i = 0; i < O.NWord; i++)
  {
    Determinant[i] = 0;
  }
  for (i = 0; i < nFermion; i++)
  {
    MinEnergy = DBL_MAX;
    for (j = 0; j < O.NFermion; j++)
    {
      if (OperatorCr(Determinant, j, O.NFermion) != 0)
      {
        WaveFunctionAddElement(psi, Determinant, 1.0, 0);
        PsiHermitianOperatorPsi(O, *psi, &Energy);
        PsiHermitianOperatorSqrPsi(O, *psi, &EnergySqr);
        OperatorAn(Determinant, j, O.NFermion);
        WaveFunctionTruncateBasis(psi, 2.0);
        if (Energy * EnergySqr < MinEnergy)
        {
          jpick = j;
          MinEnergy = Energy * EnergySqr;
        }
      }
    }
    if (OperatorCr(Determinant, jpick, O.NFermion) == 0)
    {
      printf("Error in CreateStartWaveFunctionOnsiteEnergy Fermionic part\n");
      return 1;
    };
  }
  for (i = 0; i < nBoson; i++)
  {
    MinEnergy = DBL_MAX;
    for (j = O.NFermion; j < O.NFermion + O.NBoson; j++)
    {
      if (OperatorCr(Determinant, j, O.NFermion) != 0)
      {
        WaveFunctionAddElement(psi, Determinant, 1.0, 0);
        PsiHermitianOperatorPsi(O, *psi, &Energy);
        PsiHermitianOperatorSqrPsi(O, *psi, &EnergySqr);
        OperatorAn(Determinant, j, O.NFermion);
        WaveFunctionTruncateBasis(psi, 2.0);
        if (Energy * EnergySqr < MinEnergy)
        {
          jpick = j;
          MinEnergy = Energy * EnergySqr;
        }
      }
    }
    if (OperatorCr(Determinant, jpick, O.NFermion) == 0)
    {
      printf("Error in CreateStartWaveFunctionOnsiteEnergy Bosonic part\n");
      return 1;
    };
  }
  WaveFunctionAddElement(psi, Determinant, 1.0, 0);
  return 0;
}

/*  increase the number of fermions in psiIn by one, creating psiOut, the Fermion is added at position FermionPos
 *  If occupied is one the Fermion is occupied
 *  If occupied is zero the Fermion is empty
 */
uint32_t WaveFunctionAddFermion(WaveFunctionType const psiIn, WaveFunctionType * const psiOut, uint32_t FermionPos,
                                uint32_t occupied)
{
  uint8_t det[psiIn.NWord + 1];
  uint32_t i, j;
  psiOut->NFermion = psiIn.NFermion + 1;
  psiOut->NBoson = psiIn.NBoson;
  psiOut->NBitsKey = psiIn.NBitsKey;
  psiOut->QComplex = psiIn.QComplex;
  WaveFunctionInit(psiOut);
  if (psiIn.QComplex)
  {
    for (i = 1; i <= psiIn.N; i++)
    {
      for (j = 0; j < psiIn.NWord + 1; j++)
      {
        det[j] = 0;
      }
      for (j = 1; j <= psiIn.NBoson; j++)
      {
        det[psiOut->NWord - j] = psiIn.Determinant[BlokElementDet(i, psiIn.NWord, psiIn.NWord - j)];
      }
      for (j = 0; j < FermionPos; j++)
      {
        if (psiIn.Determinant[BlokElementDet(i, psiIn.NWord, j >> 3)] & (128 >> (j & 7)))
        {
          det[j >> 3] += (128 >> (j & 7));
        }
      }
      if (occupied)
      {
        det[FermionPos >> 3] += (128 >> (FermionPos & 7));
      }
      for (j = FermionPos; j < psiIn.NFermion; j++)
      {
        if (psiIn.Determinant[BlokElementDet(i, psiIn.NWord, j >> 3)] & (128 >> (j & 7)))
        {
          det[(j + 1) >> 3] += (128 >> ((j + 1) & 7));
        }
      }
      if (ComplexWaveFunctionAddElement(psiOut, det, psiIn.ValueRe[BlokElement(i)], psiIn.ValueIm[BlokElement(i)]))
      {
        printf("RealWaveFunctionAddElement failed in WaveFunctionAddFermion\n");
        fflush(stdout);
        return 1;
      };
    }
  }
  else
  {
    for (i = 1; i <= psiIn.N; i++)
    {
      for (j = 0; j < psiIn.NWord + 1; j++)
      {
        det[j] = 0;
      }
      for (j = 1; j <= psiIn.NBoson; j++)
      {
        det[psiOut->NWord - j] = psiIn.Determinant[BlokElementDet(i, psiIn.NWord, psiIn.NWord - j)];
      }
      for (j = 0; j < FermionPos; j++)
      {
        if (psiIn.Determinant[BlokElementDet(i, psiIn.NWord, j >> 3)] & (128 >> (j & 7)))
        {
          det[j >> 3] += (128 >> (j & 7));
        }
      }
      if (occupied)
      {
        det[FermionPos >> 3] += (128 >> (FermionPos & 7));
      }
      for (j = FermionPos; j < psiIn.NFermion; j++)
      {
        if (psiIn.Determinant[BlokElementDet(i, psiIn.NWord, j >> 3)] & (128 >> (j & 7)))
        {
          det[(j + 1) >> 3] += (128 >> ((j + 1) & 7));
        }
      }
      if (RealWaveFunctionAddElement(psiOut, det, psiIn.ValueRe[BlokElement(i)]))
      {
        printf("RealWaveFunctionAddElement failed in WaveFunctionAddFermion\n");
        fflush(stdout);
        return 1;
      };
    }
  }
  return 0;
}


uint32_t WaveFunctionRemoveFermion(WaveFunctionType const psiIn, WaveFunctionType * const psiOut, uint32_t FermionPos)
{
  uint8_t det[psiIn.NWord];
  uint32_t i, j;
  psiOut->NFermion = psiIn.NFermion - 1;
  psiOut->NBoson = psiIn.NBoson;
  psiOut->NBitsKey = psiIn.NBitsKey;
  psiOut->QComplex = psiIn.QComplex;
  WaveFunctionInit(psiOut);
  if (psiIn.QComplex)
  {
    for (i = 1; i <= psiIn.N; i++)
    {
      for (j = 0; j < psiIn.NWord + 1; j++)
      {
        det[j] = 0;
      }
      for (j = 1; j <= psiIn.NBoson; j++)
      {
        det[psiOut->NWord - j] = psiIn.Determinant[BlokElementDet(i, psiIn.NWord, psiIn.NWord - j)];
      }
      for (j = 0; j < FermionPos; j++)
      {
        if (psiIn.Determinant[BlokElementDet(i, psiIn.NWord, j >> 3)] & (128 >> (j & 7)))
        {
          det[j >> 3] += (128 >> (j & 7));
        }
      }
      for (j = FermionPos + 1; j < psiIn.NFermion; j++)
      {
        if (psiIn.Determinant[BlokElementDet(i, psiIn.NWord, j >> 3)] & (128 >> (j & 7)))
        {
          det[(j - 1) >> 3] += (128 >> ((j - 1) & 7));
        }
      }
      if (ComplexWaveFunctionAddElement(psiOut, det, psiIn.ValueRe[BlokElement(i)], psiIn.ValueIm[BlokElement(i)]))
      {
        printf("RealWaveFunctionAddElement failed in WaveFunctionAddFermion\n");
        fflush(stdout);
        return 1;
      };
    }
  }
  else
  {
    for (i = 1; i <= psiIn.N; i++)
    {
      for (j = 0; j < psiIn.NWord + 1; j++)
      {
        det[j] = 0;
      }
      for (j = 1; j <= psiIn.NBoson; j++)
      {
        det[psiOut->NWord - j] = psiIn.Determinant[BlokElementDet(i, psiIn.NWord, psiIn.NWord - j)];
      }
      for (j = 0; j < FermionPos; j++)
      {
        if (psiIn.Determinant[BlokElementDet(i, psiIn.NWord, j >> 3)] & (128 >> (j & 7)))
        {
          det[j >> 3] += (128 >> (j & 7));
        }
      }
      for (j = FermionPos + 1; j < psiIn.NFermion; j++)
      {
        if (psiIn.Determinant[BlokElementDet(i, psiIn.NWord, j >> 3)] & (128 >> (j & 7)))
        {
          det[(j - 1) >> 3] += (128 >> ((j - 1) & 7));
        }
      }
      if (RealWaveFunctionAddElement(psiOut, det, psiIn.ValueRe[BlokElement(i)]))
      {
        printf("RealWaveFunctionAddElement failed in WaveFunctionAddFermion\n");
        return 1;
      };
    }
  }
  return 0;
}


uint32_t WaveFunctionCopy(WaveFunctionType const A, WaveFunctionType * const B)
{
  B->NFermion = A.NFermion;
  B->NBoson = A.NBoson;
  B->QComplex = A.QComplex;
  B->NBitsKey = A.NBitsKey;
  if (WaveFunctionInit(B))
  {
    printf("WaveFunctionInit failed in WaveFunctionCopy\n");
    return 1;
  };
  if (WaveFunctionAddUnionBasis(B, A))
  {
    printf("WaveFunctionAddUnionBasis failed in WaveFunctionCopy\n");
    return 1;
  };
  return 0;
}
