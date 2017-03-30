#include "BasicMath_BlockListOfPoles.c"
#include "BasicMath_BlockAndersonMatrix.c"
#include "BasicMath_BlockTridiagonalMatrix.c"
#include "BasicMath_BlockNaturalImpurityOrbital.c"
#include "BasicMath_BlockTypeConversion.c"
#include "BasicMath_BlockMatrixToOperator.c"

/* Free the BlockMatrix, independent if the type is capital or small */
/* Block types read from the lua stack need different freeing depending on capital or small cap */
/* use FreeBlockMatrixFromLua instead */
uint32_t FreeBlockMatrix(BlockMatrixType *M)
{
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
uint32_t FreeBlockMatrixFromLua(BlockMatrixType *M)
{
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
 * change the type of a response funciton from M->type to type
 */
BlockMatrixType *ConvertBlockMatrix(BlockMatrixType *M, const char type, double T)
{
  BlockMatrixType *res;
  res = malloc(sizeof(BlockMatrixType));
  res->type = type;
  switch (type|32)
  {
    case 'l':
    {
      switch (M->type|32)
      {
        case 'l':
        {
          CopyBlockListOfPoles(*(BlockListOfPoles*)M->data, (BlockListOfPoles*)res->data);
        }
        break;
        case 't':
        {
          BlockListOfPolesToBlockTridiagonalMatrix(*(BlockListOfPoles*)M->data, (BlockTriDiagonalMatrixType*)res->data, ((BlockListOfPoles*)(M->data))->N+1, 1000*DBL_EPSILON);
        }
        break;
        case 'a':
        {
          BlockListOfPolesToBlockAndersonMatrix(*(BlockListOfPoles*)M->data, (BlockAndersonMatrixType*)res->data, ((BlockListOfPoles*)(M->data))->N+1);
        }
        break;
        case 'n':
        {
          BlockListOfPolesToBlockNaturalImpurityOrbitalMatrix((BlockListOfPoles*)M->data, (BlockNaturalImpurityOrbitalMatrixType*)res->data, T, ((BlockListOfPoles*)(M->data))->N+1, ((BlockListOfPoles*)(M->data))->N+1);
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
    case 't':
    {
      switch (M->type|32)
      {
        case 'l':
        {
          
        }
        break;
        case 't':
        {
          CopyBlockTridiagonalMatrix(*(BlockTriDiagonalMatrixType*)M->data, (BlockTriDiagonalMatrixType*)res->data);
        }
        break;
        case 'a':
        {
          
        }
        break;
        case 'n':
        {
          
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
    case 'a':
    {
      switch (M->type|32)
      {
        case 'l':
        {
          
        }
        break;
        case 't':
        {
          
        }
        break;
        case 'a':
        {
          CopyBlockAndersonMatrix(*(BlockAndersonMatrixType*)M->data, (BlockAndersonMatrixType*)res->data);
        }
        break;
        case 'n':
        {
          
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
    case 'n':
    {
      switch (M->type|32)
      {
        case 'l':
        {
          
        }
        break;
        case 't':
        {
          
        }
        break;
        case 'a':
        {
          
        }
        break;
        case 'n':
        {
          CopyBlockNaturalImpurityOrbitalMatrix(*(BlockNaturalImpurityOrbitalMatrixType*)M->data, (BlockNaturalImpurityOrbitalMatrixType*)res->data);
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
    default:
    {
      printf("Error in ConvertBlockMatrix unrecognized type to convert to: %c",type);
    }
    break;
  }
  return res;
}
