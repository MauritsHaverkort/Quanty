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
  if(M->NSugestedGridPoints!=0)
  {
    free(M->SugestedGrid);
    M->SugestedGrid=NULL;
    M->NSugestedGridPoints=0;
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
uint32_t FreeBlockMatrixFromLua(BlockMatrixType *M)
{
  if(M->NSugestedGridPoints!=0)
  {
    free(M->SugestedGrid);
    M->SugestedGrid=NULL;
    M->NSugestedGridPoints=0;
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
 * change the type of a response function from M->type to type
 */
BlockMatrixType *ConvertBlockMatrix(BlockMatrixType *M, const char type, double T)
{
  BlockMatrixType *res;
  res = malloc(sizeof(BlockMatrixType));
  res->type = type;
  res->NSugestedGridPoints=M->NSugestedGridPoints;
  if(res->NSugestedGridPoints==0)
  {
    res->SugestedGrid=NULL;
  }
  else
  {
    res->SugestedGrid=malloc(res->NSugestedGridPoints*sizeof(double));
    memcpy(res->SugestedGrid, M->SugestedGrid, res->NSugestedGridPoints*sizeof(double));
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
          BlockListOfPolesToBlockNaturalImpurityOrbitalMatrix((BlockListOfPoles*)M->data, (BlockNaturalImpurityOrbitalMatrixType*)res->data, T, ((BlockListOfPoles*)(M->data))->N+1, ((BlockListOfPoles*)(M->data))->N+1);
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
          
          //TODO: transform tri to list of poles
          //BlockTriDiagonalMatrixToBlockListOfPolesFixedGrid(*(BlockTriDiagonalMatrixType*)M->data, (BlockListOfPoles*)res->data, <#const double Emin#>, <#const double Emax#>, <#const uint32_t NE#>)
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

        }
        break;
        case 'n':
        {

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
    case 'n': // original natural impurity orbitals
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
