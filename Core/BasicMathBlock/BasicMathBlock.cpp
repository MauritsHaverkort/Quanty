#include "BasicMath_BlockAndersonMatrix.cpp"
#include "BasicMath_BlockListOfPoles.cpp"
#include "BasicMath_BlockMatrixToOperator.cpp"
#include "BasicMath_BlockNaturalImpurityOrbital.cpp"
#include "BasicMath_BlockTridiagonalMatrix.cpp"
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
ResponsefunctionType *ConvertBlockMatrix(ResponsefunctionType *M, const char type, double T)
{
  ResponsefunctionType *res;
  if(M->Block!=true) printf("ERROR: called ConvertBlockMatrix on a ResponsefunctionType that is not in Block Form.\n");
  res = (ResponsefunctionType*)malloc(sizeof(ResponsefunctionType));
  res->NSuggestedGridPoints = 0;
  res->type = type;
  res->Block = true;
  res->NSuggestedGridPoints=M->NSuggestedGridPoints;
  if(res->NSuggestedGridPoints==0)
  {
    res->SuggestedGrid=NULL;
  }
  else
  {
    res->SuggestedGrid= (double*)malloc(res->NSuggestedGridPoints*sizeof(double));
    memcpy(res->SuggestedGrid, M->SuggestedGrid, res->NSuggestedGridPoints*sizeof(double));
  }
  switch (M->type|32)      // original type
  {
    case 'l': // original list of poles
    {
      switch (type|32)     // target type
      {
        case 'l':
        {
          res->data = (void*) malloc(sizeof(BlockListOfPoles));
          CopyBlockListOfPoles(*(BlockListOfPoles*)M->data, (BlockListOfPoles*)res->data);
        }
        break;
        case 't':
        {
          res->data = (void*) malloc(sizeof(BlockTriDiagonalMatrixType));
          BlockListOfPolesToBlockTridiagonalMatrix(*(BlockListOfPoles*)M->data, (BlockTriDiagonalMatrixType*)res->data, ((BlockListOfPoles*)(M->data))->N+1, 1000*DBL_EPSILON);
        }
        break;
        case 'a':
        {
          res->data = (void*) malloc(sizeof(BlockAndersonMatrixType));
          BlockListOfPolesToBlockAndersonMatrix(*(BlockListOfPoles*)M->data, (BlockAndersonMatrixType*)res->data, ((BlockListOfPoles*)(M->data))->N+1);
        }
        break;
        case 'n':
        {
          res->data = (void*) malloc(sizeof(BlockNaturalImpurityOrbitalMatrixType));
          //TODO: get res->data->N right for block list of poles
          //TODO: allow flexible set of states in tri-diag form with options
         // BlockListOfPolesToBlockNaturalImpurityOrbitalMatrix((BlockListOfPoles*)(M->data), (BlockNaturalImpurityOrbitalMatrixType*)(res->data), T, ((BlockListOfPoles*)(M->data))->N+1, ((BlockListOfPoles*)(M->data))->N+1);
          BlockListOfPolesToBlockNaturalImpurityOrbitalMatrix((BlockListOfPoles*)(M->data), (BlockNaturalImpurityOrbitalMatrixType*)(res->data), T, 10, 10);
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
          res->data = (void*) malloc(sizeof(BlockListOfPoles));
          if(M->NSuggestedGridPoints!=0)
          {
            BlockTriDiagonalMatrixToBlockListOfPolesVariableGrid((BlockTriDiagonalMatrixType*)M->data, (BlockListOfPoles*)res->data, M->SuggestedGrid, M->NSuggestedGridPoints);
          }
          else
          {
            BlockTriDiagonalMatrixToBlockListOfPoles((BlockTriDiagonalMatrixType*)M->data, (BlockListOfPoles*)res->data);
          }
        }
        break;
        case 't':
        {
          res->data = (void*) malloc(sizeof(BlockTriDiagonalMatrixType));
          CopyBlockTridiagonalMatrix(*(BlockTriDiagonalMatrixType*)M->data, (BlockTriDiagonalMatrixType*)res->data);
        }
        break;
        case 'a':
        {
          res->data = (void*) malloc(sizeof(BlockAndersonMatrixType));
          BlockTriDiagonalMatrixToBlockAndersonMatrix(*(BlockTriDiagonalMatrixType*)M->data, (BlockAndersonMatrixType*)res->data);
        }
        break;
        case 'n':
        {
          res->data = (void*) malloc(sizeof(BlockNaturalImpurityOrbitalMatrixType));
          BlockTriDiagonalMatrixToBlockNaturalImpurityMatrix(*(BlockTriDiagonalMatrixType*)M->data,(BlockNaturalImpurityOrbitalMatrixType*)res->data);
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
          res->data = (void*) malloc(sizeof(BlockListOfPoles));
          if(M->NSuggestedGridPoints!=0)
          {
            BlockAndersonMatrixToBlockListOfPolesVariableGrid(*(BlockAndersonMatrixType*)M->data, (BlockListOfPoles*)res->data, M->SuggestedGrid, M->NSuggestedGridPoints);
          }
          else
          {
            BlockAndersonMatrixToBlockListOfPoles(*(BlockAndersonMatrixType*)M->data, (BlockListOfPoles*)res->data);
          }
        }
        break;
        case 't':
        {
          res->data = (void*) malloc(sizeof(BlockTriDiagonalMatrixType));
          BlockAndersonMatrixToBlockTridiagonalMatrix(*(BlockAndersonMatrixType*)M->data, (BlockTriDiagonalMatrixType*)res->data);
        }
        break;
        case 'a':
        {
          CopyBlockAndersonMatrix(*(BlockAndersonMatrixType*)M->data, (BlockAndersonMatrixType*)res->data);
        }
        break;
        case 'n':
        {
          res->data = (void*) malloc(sizeof(BlockNaturalImpurityOrbitalMatrixType));
          BlockAndersonMatrixToBlockNaturalImpurityOrbitalMatrix(*(BlockAndersonMatrixType*)M->data,(BlockNaturalImpurityOrbitalMatrixType*)res->data);
        }
        break;
        default:
        {
          printf("Error in ConvertBlockMatrix unrecognized type of M: %c",M->type);
        }
        break;
      }
    }
    break;
    case 'n': // original natural impurity orbitals
    {
      switch (M->type|32)
      {
        case 'l':
        {
          res->data = (void*) malloc(sizeof(BlockListOfPoles));
          if(M->NSuggestedGridPoints!=0)
          {
            BlockNaturalImpurityOrbitalMatrixToBlockListOfPolesVariableGrid(*(BlockNaturalImpurityOrbitalMatrixType*)M->data, (BlockListOfPoles*)res->data, M->SuggestedGrid, M->NSuggestedGridPoints);
          }
          else
          {
            BlockNaturalImpurityOrbitalMatrixToBlockListOfPoles(*(BlockNaturalImpurityOrbitalMatrixType*)M->data, (BlockListOfPoles*)res->data);
          }
        }
        break;
        case 't':
        {
          res->data = (void*) malloc(sizeof(BlockTriDiagonalMatrixType));
          BlockNaturalImpurityOrbitalMatrixToBlockTridiagonalMatrix(*(BlockNaturalImpurityOrbitalMatrixType*)M->data, (BlockTriDiagonalMatrixType*)res->data);
        }
        break;
        case 'a':
        {
          res->data = (void*) malloc(sizeof(BlockAndersonMatrixType));
          BlockNaturalImpurityOrbitalMatrixToBlockAndersonMatrix(*(BlockNaturalImpurityOrbitalMatrixType*)M->data, (BlockAndersonMatrixType*)res->data);

        }
        break;
        case 'n':
        {
          res->data = (void*) malloc(sizeof(BlockNaturalImpurityOrbitalMatrixType));
          CopyBlockNaturalImpurityOrbitalMatrix(*(BlockNaturalImpurityOrbitalMatrixType*)M->data, (BlockNaturalImpurityOrbitalMatrixType*)res->data);
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
      printf("Error in ConvertBlockMatrix unrecognized type of M: %c",M->type);
    }
    break;
  }
  return res;
}
