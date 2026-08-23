#include "BasicMath_BlockAndersonMatrix.cpp"
#include "BasicMath_BlockListOfPoles.cpp"
#include "BasicMath_BlockMatrixToOperator.cpp"
#include "BasicMath_BlockNaturalImpurityOrbital.cpp"
#include "BasicMath_BlockTridiagonalMatrix.cpp"
#include "BasicMath_BlockDoubleTridiagnonalMatrix.cpp"
#include "BasicMath_BlockTypeConversion.cpp"

/* Free the BlockMatrix, independent if the type is capital or small */
/* Block types read from the lua stack need different freeing depending on capital or small cap */
/* use FreeBlockMatrixFromLua instead */
uint32_t FreeBlockMatrix(ResponsefunctionType *M)
{
  if(M->NSuggestedGridPoints!=0)
  {
    free(M->SuggestedGrid);
    M->SuggestedGrid=NULL;
    M->NSuggestedGridPoints=0;
  }
  switch (M->type|32) {
    case 'l':
    {
      FreeBlockListOfPoles((BlockListOfPoles *)M->data);
      free((BlockListOfPoles *)(M->data));
      return 0;
    }
    case 'a':
    {
      FreeBlockAndersonMatrix((BlockAndersonMatrixType *)M->data);
      free((BlockAndersonMatrixType *)(M->data));
      return 0;
    }
    case 't':
    {
      FreeBlockTridiagonalMatrix((BlockTriDiagonalMatrixType *)M->data);
      free((BlockTriDiagonalMatrixType *)(M->data));
      return 0;
    }
    case 'n':
    {
      FreeBlockNaturalImpurityOrbitalMatrix((BlockNaturalImpurityOrbitalMatrixType *)M->data);
      free((BlockNaturalImpurityOrbitalMatrixType *)(M->data));
      return 0;
    }
    default:
    {
      printf("Error in FreeBlockMatrix unrecognized type %c",M->type);
      return 1;
    }
  }
}

/* capital read from table and needs to be freed, small read from user data type */
uint32_t FreeBlockMatrixFromLua(ResponsefunctionType *M)
{
  if((M->type>='A') && (M->type<='Z') && (M->NSuggestedGridPoints!=0))
  {
    free(M->SuggestedGrid);
    M->SuggestedGrid=NULL;
    M->NSuggestedGridPoints=0;
  }
  switch (M->type) {
    case 'l':
    {
      return 0;
    }
    case 'L':
    {
      FreeBlockListOfPoles((BlockListOfPoles *)M->data);
      free((BlockListOfPoles *)(M->data));
      return 0;
    }
    case 'a':
    {
      return 0;
    }
    case 'A':
    {
      FreeBlockAndersonMatrix((BlockAndersonMatrixType *)M->data);
      free((BlockAndersonMatrixType *)(M->data));
      return 0;
    }
    case 't':
    {
      return 0;
    }
    case 'T':
    {
      FreeBlockTridiagonalMatrix((BlockTriDiagonalMatrixType *)M->data);
      free((BlockTriDiagonalMatrixType *)(M->data));
      return 0;
    }
    case 'n':
    {
      return 0;
    }
    case 'N':
    {
      FreeBlockNaturalImpurityOrbitalMatrix((BlockNaturalImpurityOrbitalMatrixType *)M->data);
      free((BlockNaturalImpurityOrbitalMatrixType *)(M->data));
      return 0;
    }
    case 'd':
    {
      return 0;
    }
    case 'D':
    {
      FreeBlockDoubleTriDiagonalMatrix((BlockDoubleTriDiagonalMatrixType *)M->data);
      free((BlockDoubleTriDiagonalMatrixType *)(M->data));
      return 0;
    }
    default:
    {
      printf("Error in FreeBlockMatrixFromLua unrecognized type %c",M->type);
      return 1;
    }
  }
}

/*
 * change the type of a Block response function from M->type to type
 */
int ConvertBlockMatrix(const ResponsefunctionType *const Min, ResponsefunctionType *Mout, const char type, double T, double SingularValue)
{
  if(Min->Block!=true) printf("ERROR: called ConvertBlockMatrix on a ResponsefunctionType that is not in Block Form.\n");
  Mout->NSuggestedGridPoints = 0;
  Mout->type = type;
  Mout->Block = true;
  Mout->NSuggestedGridPoints=Min->NSuggestedGridPoints;
  if(Mout->NSuggestedGridPoints==0)
  {
    Mout->SuggestedGrid=NULL;
  }
  else
  {
    Mout->SuggestedGrid= (double*)malloc(Mout->NSuggestedGridPoints*sizeof(double));
    memcpy(Mout->SuggestedGrid, Min->SuggestedGrid, Mout->NSuggestedGridPoints*sizeof(double));
  }
  switch (Min->type|32)      // original type
  {
    case 'l': // original list of poles
    {
      switch (type|32)     // target type
      {
        case 'l':
        {
          Mout->data = (void*) malloc(sizeof(BlockListOfPoles));
          CopyBlockListOfPoles(*(BlockListOfPoles*)Min->data, (BlockListOfPoles*)Mout->data);
        }
        break;
        case 't':
        {
          Mout->data = (void*) malloc(sizeof(BlockTriDiagonalMatrixType));
          BlockListOfPolesToBlockTridiagonalMatrix(*(BlockListOfPoles*)Min->data, (BlockTriDiagonalMatrixType*)Mout->data, ((BlockListOfPoles*)(Min->data))->N, SingularValue);
        }
        break;
        case 'a':
        {
          Mout->data = (void*) malloc(sizeof(BlockAndersonMatrixType));
          BlockListOfPolesToBlockAndersonMatrix(*(BlockListOfPoles*)Min->data, (BlockAndersonMatrixType*)Mout->data,  ((BlockListOfPoles*)(Min->data))->N+1, SingularValue);
        }
        break;
        case 'n':
        {
          Mout->data = (void*) malloc(sizeof(BlockNaturalImpurityOrbitalMatrixType));
          BlockListOfPolesToBlockNaturalImpurityOrbitalMatrix((BlockListOfPoles*)(Min->data), (BlockNaturalImpurityOrbitalMatrixType*)(Mout->data), T, ((BlockListOfPoles*)(Min->data))->N+1, ((BlockListOfPoles*)(Min->data))->N+1);
        }
        break;
        default:
        {
          printf("Error in ConvertBlockMatrix unrecognized type to convert to: %c", type);
        }
        break;
      }
    }
    break;
    case 't': // original tri-diagonal
    {
      switch (type|32)
      {
        case 'l':
        {
          Mout->data = (void*) malloc(sizeof(BlockListOfPoles));
          if(Min->NSuggestedGridPoints!=0)
          {
            BlockTriDiagonalMatrixToBlockListOfPolesVariableGrid((BlockTriDiagonalMatrixType*)Min->data, (BlockListOfPoles*)Mout->data, Min->SuggestedGrid, Min->NSuggestedGridPoints);
          }
          else
          {
            BlockTriDiagonalMatrixToBlockListOfPoles((BlockTriDiagonalMatrixType*)Min->data, (BlockListOfPoles*)Mout->data);
          }
        }
        break;
        case 't':
        {
          Mout->data = (void*) malloc(sizeof(BlockTriDiagonalMatrixType));
          CopyBlockTridiagonalMatrix(*(BlockTriDiagonalMatrixType*)Min->data, (BlockTriDiagonalMatrixType*)Mout->data);
        }
        break;
        case 'a':
        {
          Mout->data = (void*) malloc(sizeof(BlockAndersonMatrixType));
          BlockTriDiagonalMatrixToBlockAndersonMatrix(*(BlockTriDiagonalMatrixType*)Min->data, (BlockAndersonMatrixType*)Mout->data);
        }
        break;
        case 'n':
        {
          Mout->data = (void*) malloc(sizeof(BlockNaturalImpurityOrbitalMatrixType));
          BlockTriDiagonalMatrixToBlockNaturalImpurityMatrix(*(BlockTriDiagonalMatrixType*)Min->data,(BlockNaturalImpurityOrbitalMatrixType*)Mout->data,T);
        }
        break;
        default:
        {
          printf("Error in ConvertBlockMatrix unrecognized type to convert to: %c",type);
        }
        break;
      }
    }
    break;
    case 'a': // original anderson matrix
    {
      switch (type|32)
      {
        case 'l':
        {
          Mout->data = (void*) malloc(sizeof(BlockListOfPoles));
          if(Min->NSuggestedGridPoints!=0)
          {
            BlockAndersonMatrixToBlockListOfPolesVariableGrid(*(BlockAndersonMatrixType*)Min->data, (BlockListOfPoles*)Mout->data, Min->SuggestedGrid, Min->NSuggestedGridPoints);
          }
          else
          {
            BlockAndersonMatrixToBlockListOfPoles((BlockAndersonMatrixType*)Min->data, (BlockListOfPoles*)Mout->data);
          }
        }
        break;
        case 't':
        {
          Mout->data = (void*) malloc(sizeof(BlockTriDiagonalMatrixType));
          BlockAndersonMatrixToBlockTridiagonalMatrix((BlockAndersonMatrixType*)Min->data, (BlockTriDiagonalMatrixType*)Mout->data, SingularValue);
        }
        break;
        case 'a':
        {
          Mout->data = (void*) malloc(sizeof(BlockAndersonMatrixType));
          CopyBlockAndersonMatrix(*(BlockAndersonMatrixType*)Min->data, (BlockAndersonMatrixType*)Mout->data);
        }
        break;
        case 'n':
        {
          Mout->data = (void*) malloc(sizeof(BlockNaturalImpurityOrbitalMatrixType));
          BlockAndersonMatrixToBlockNaturalImpurityOrbitalMatrix((BlockAndersonMatrixType*)Min->data,(BlockNaturalImpurityOrbitalMatrixType*)Mout->data,T);
        }
        break;
        default:
        {
          printf("Error in ConvertBlockMatrix unrecognized type: %c",type);
        }
        break;
      }
    }
    break;
    case 'n': // original natural impurity orbitals
    {
      switch (type|32)
      {
        case 'l':
        {
          Mout->data = (void*) malloc(sizeof(BlockListOfPoles));
          if(Min->NSuggestedGridPoints!=0)
          {
            BlockNaturalImpurityOrbitalMatrixToBlockListOfPolesVariableGrid(*(BlockNaturalImpurityOrbitalMatrixType*)Min->data, (BlockListOfPoles*)Mout->data, Min->SuggestedGrid, Min->NSuggestedGridPoints);
          }
          else
          {
            BlockNaturalImpurityOrbitalMatrixToBlockListOfPoles((BlockNaturalImpurityOrbitalMatrixType*)Min->data, (BlockListOfPoles*)Mout->data);
          }
        }
        break;
        case 't':
        {
          Mout->data = (void*) malloc(sizeof(BlockTriDiagonalMatrixType));
          BlockNaturalImpurityOrbitalMatrixToBlockTridiagonalMatrix((BlockNaturalImpurityOrbitalMatrixType*)Min->data, (BlockTriDiagonalMatrixType*)Mout->data, SingularValue);
        }
        break;
        case 'a':
        {
          Mout->data = (void*) malloc(sizeof(BlockAndersonMatrixType));
          BlockNaturalImpurityOrbitalMatrixToBlockAndersonMatrix((BlockNaturalImpurityOrbitalMatrixType*)Min->data, (BlockAndersonMatrixType*)Mout->data);

        }
        break;
        case 'n':
        {
          Mout->data = (void*) malloc(sizeof(BlockNaturalImpurityOrbitalMatrixType));
          CopyBlockNaturalImpurityOrbitalMatrix(*(BlockNaturalImpurityOrbitalMatrixType*)Min->data, (BlockNaturalImpurityOrbitalMatrixType*)Mout->data);
        }
        break;
        default:
        {
          printf("Error in ConvertBlockMatrix unrecognized type to convert to: %c",type);
        }
        break;
      }
    }
    break;
    default:
    {
      printf("Error in ConvertBlockMatrix unrecognized type of M: %c",Min->type);
    }
    break;
  }
  return 0;
}

int BlockResponseFunctionRotate(const ResponsefunctionType *const G, const CompactMatrixType *const U, ResponsefunctionType * const GRot)
{
  GRot->type  = G->type;
  GRot->Block = G->Block;
  GRot->NSuggestedGridPoints = G->NSuggestedGridPoints;
  if(GRot->NSuggestedGridPoints==0)
  {
    GRot->SuggestedGrid=NULL;
  }
  else
  {
    GRot->SuggestedGrid= (double*)malloc(GRot->NSuggestedGridPoints*sizeof(double));
    memcpy(GRot->SuggestedGrid, G->SuggestedGrid, GRot->NSuggestedGridPoints*sizeof(double));
  }
  switch (G->type|32)      // original type
  {
    case 'l': // original list of poles
    {
      char Name[MAXStrLength];
      BlockListOfPoles *GRotData, *GData;
      CompactMatrixType Gi, GRoti;
      GData = ((BlockListOfPoles*)(G->data));
      GRotData = (BlockListOfPoles*)malloc(sizeof(BlockListOfPoles));
      GRot->data = GRotData;
      GRotData->BlockSize = U->Ni;
      GRotData->QComplex = GData->QComplex;
      if(U->QComplex!=0) GRotData->QComplex = 1;
      GRotData->N = GData->N;
      InitBlockListOfPoles(GRotData);
      strcpy(GRotData->Name, GData->Name);
      GRotData->mu = GData->mu;
      Gi.QComplex = GData->QComplex;
      Gi.Ni = GData->BlockSize;
      Gi.Nj = GData->BlockSize;
      Gi.Name = Name;
      GRoti.QComplex = GRotData->QComplex;
      GRoti.Ni = GRotData->BlockSize;
      GRoti.Nj = GRotData->BlockSize;
      GRoti.Name = Name;
      // Rotate a0
      if(Gi.QComplex==0)
      {
        Gi.VRe = GData->a0Re;
      }
      else
      {
        Gi.V = GData->a0;
      }
      if(GRoti.QComplex==0)
      {
        GRoti.VRe = GRotData->a0Re;
      }
      else
      {
        GRoti.V = GRotData->a0;
      }
      CompactMatrixMultiply3PreAllocated('*', 'N', 'T', U, &Gi, U, &GRoti);
      // Rotate the residus of the poles
      for(int i=0; i<GData->N; i++)
      {
        if(Gi.QComplex==0)
        {
          Gi.VRe = &(GData->VRe[i*(Sqr(GData->BlockSize)+1)+1]);
        }
        else
        {
          Gi.V = &(GData->V[i*(Sqr(GData->BlockSize)+1)+1]);
          GRotData->V[i*(Sqr(GRotData->BlockSize)+1)]=GData->V[i*(Sqr(GData->BlockSize)+1)];
        }
        if(GRoti.QComplex==0)
        {
          GRoti.VRe = &(GRotData->VRe[i*(Sqr(GRotData->BlockSize)+1)+1]);
          GRotData->VRe[i*(Sqr(GRotData->BlockSize)+1)]=GData->VRe[i*(Sqr(GData->BlockSize)+1)];
        }
        else
        {
          GRoti.V = &(GRotData->V[i*(Sqr(GRotData->BlockSize)+1)+1]);
          if(Gi.QComplex==0)
          {
            GRotData->V[i*(Sqr(GRotData->BlockSize)+1)].Re=GData->VRe[i*(Sqr(GData->BlockSize)+1)];
            GRotData->V[i*(Sqr(GRotData->BlockSize)+1)].Im=0;
          }
        }
        CompactMatrixMultiply3PreAllocated('*', 'N', 'T', U, &Gi, U, &GRoti);
      }
    }
    break;
    case 't': // original tri-diagonal
    {
      ResponsefunctionType Gp;
      ConvertBlockMatrix(G, &Gp, 'L', 0, 0);
      BlockResponseFunctionRotate(&Gp, U, GRot);
      FreeResponsefunction(&Gp);
      return 0;
    }
    break;
    case 'a': // original anderson matrix
    {
      ResponsefunctionType Gp;
      ConvertBlockMatrix(G, &Gp, 'L', 0, 0);
      BlockResponseFunctionRotate(&Gp, U, GRot);
      FreeResponsefunction(&Gp);
      return 0;
    }
    break;
    case 'n': // original natural impurity orbitals
    {
      ResponsefunctionType Gp;
      ConvertBlockMatrix(G, &Gp, 'L', 0, 0);
      BlockResponseFunctionRotate(&Gp, U, GRot);
      FreeResponsefunction(&Gp);
      return 0;
    }
    break;
    default:
    {
      printf("Error in BlockResponseFunctionRotate unrecognized type of G: %c",G->type);
    }
    break;
  }
  return 0;
}

