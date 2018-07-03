#include "BasicMath_BlockTridiagonalMatrix.h"

/*
 *   Initialize a struct of BlockTridiagonalMatrix type
 *   malloc the space for storing N+1 matrices of dimension BlockSize * BlockSize storing a0 to aN
 *   malloc the space for storing N   matrices of dimension BlockSize * BlockSize storing b0 to bN-1
 */
uint32_t InitBlockTridiagonalMatrix(BlockTriDiagonalMatrixType *const A)
{
  int i;
  strcpy(A->Name, "Block Tridiagonal Matrix");
  if(A->QComplex)
  {
    A->a= (Complex16*)calloc((A->N+1)*(A->BlockSize*A->BlockSize),sizeof(Complex16));
    if(A->a==NULL)
    {
      printf("malloc failed in InitBlockTridiagonalMatrix\n");
      return 1;
    }
    A->b= (Complex16*)calloc((A->N)*(A->BlockSize*A->BlockSize),sizeof(Complex16));
    if(A->b==NULL)
    {
      free(A->a);
      printf("malloc failed in InitBlockTridiagonalMatrix\n");
      return 1;
    }
    for(i=0;i<A->BlockSize;i++)
    {
      A->b[i+i*A->BlockSize].Re=1;
    }
  }
  else
  {
    A->aRe= (double*)calloc((A->N+1)*(A->BlockSize*A->BlockSize),sizeof(double));
    if(A->aRe==NULL)
    {
      printf("malloc failed in InitBlockTridiagonalMatrix\n");
      return 1;
    }
    A->bRe= (double*)calloc((A->N)*(A->BlockSize*A->BlockSize),sizeof(double));
    if(A->bRe==NULL)
    {
      free(A->aRe);
      printf("malloc failed in InitBlockTridiagonalMatrix\n");
      return 1;
    }
    for(i=0;i<A->BlockSize;i++)
    {
      A->bRe[i+i*A->BlockSize]=1;
    }
  }
  A->mu=0;
  return 0;
}

uint32_t FreeBlockTridiagonalMatrix(BlockTriDiagonalMatrixType *A)
{
  if(A->QComplex)
  {
    free(A->a);
    free(A->b);
  }
  else
  {
    free(A->aRe);
    free(A->bRe);
  }
  return 0;
}

uint32_t ComplexBlockTridiagonalMatrixToReal(BlockTriDiagonalMatrixType *const A)
{
  uint32_t i;
  if(A->QComplex==0) return 0;
  A->QComplex=0;
  A->aRe=(double*)(A->a);
  A->bRe=(double*)(A->b);
  for(i=0;i<A->BlockSize*A->BlockSize*(A->N+1);i++)
  {
    A->aRe[i]=A->a[i].Re;
  }
  for(i=0;i<A->BlockSize*A->BlockSize*(A->N);i++)
  {
    A->bRe[i]=A->b[i].Re;
  }
  A->aRe= (double*)realloc(A->aRe, A->BlockSize*A->BlockSize*(1+A->N)*sizeof(double));
  A->bRe= (double*)realloc(A->bRe, A->BlockSize*A->BlockSize*(A->N)*sizeof(double));
  A->QComplex=0;
  return 0;
}

uint32_t BlockTridiagonalMatrixChop(BlockTriDiagonalMatrixType *const A, double const Epsilon)
{
  uint32_t i,KeepComplex=0;
  if(A->QComplex)
  {
    for(i=0;i<A->BlockSize*A->BlockSize*(A->N+1);i++)
    {
      if(Abs(A->a[i].Re)<Epsilon)A->a[i].Re=0;
      if(Abs(A->a[i].Im)<Epsilon)
      {
        A->a[i].Im=0;
      }
      else
      {
        KeepComplex=1;
      }
    }
    for(i=0;i<A->BlockSize*A->BlockSize*(A->N);i++)
    {
      if(Abs(A->b[i].Re)<Epsilon)A->b[i].Re=0;
      if(Abs(A->b[i].Im)<Epsilon)
      {
        A->b[i].Im=0;
      }
      else
      {
        KeepComplex=1;
      }
    }
    if(!KeepComplex) ComplexBlockTridiagonalMatrixToReal(A);
  }
  else
  {
    for(i=0;i<A->BlockSize*A->BlockSize*(A->N+1);i++)
    {
      if(Abs(A->aRe[i])<Epsilon)A->aRe[i]=0;
    }
    for(i=0;i<A->BlockSize*A->BlockSize*(A->N);i++)
    {
      if(Abs(A->bRe[i])<Epsilon)A->bRe[i]=0;
    }
  }
  return 0;
}

void ShiftBlockTridiagonalMatrix(BlockTriDiagonalMatrixType *const A, double const delta)
{
  uint32_t i,k;
  if(A->QComplex)
  {
    for(k=1;k<=A->N;k++)
    {
      for(i=k*A->BlockSize*A->BlockSize;i<(k+1)*A->BlockSize*A->BlockSize;i+=A->BlockSize+1)
      {
        A->a[i].Re+=delta;
      }
    }
  }
  else
  {
    for(k=1;k<=A->N;k++)
    {
      for(i=k*A->BlockSize*A->BlockSize;i<(k+1)*A->BlockSize*A->BlockSize;i+=A->BlockSize+1)
      {
        A->aRe[i]+=delta;
      }
    }
  }
  A->mu+=delta;
}

/*
 *  B <- A
 *  B not initialized on entry
 */
uint32_t CopyBlockTridiagonalMatrix(BlockTriDiagonalMatrixType A, BlockTriDiagonalMatrixType *B)
{
  *B=A;
  if(B->QComplex==0)
  {
    B->aRe = (double*)malloc((B->N+1)*B->BlockSize*B->BlockSize*sizeof(double));
    if(B->aRe==NULL)
    {
      printf("malloc failed in CopyBlockTridiagonalMatrix");
      return 1;
    }
    memcpy(B->aRe, A.aRe, (B->N+1)*B->BlockSize*B->BlockSize*sizeof(double));
    B->bRe = (double*)malloc((B->N)*B->BlockSize*B->BlockSize*sizeof(double));
    if(B->bRe==NULL)
    {
      free(B->aRe);
      printf("malloc failed in CopyBlockTridiagonalMatrix");
      return 1;
    }
    memcpy(B->bRe, A.bRe, (B->N)*B->BlockSize*B->BlockSize*sizeof(double));
  }
  else
  {
    B->a = (Complex16*)malloc((B->N+1)*B->BlockSize*B->BlockSize*sizeof(Complex16));
    if(B->a==NULL)
    {
      printf("malloc failed in CopyBlockTridiagonalMatrix");
      return 1;
    }
    memcpy(B->a, A.a, (B->N+1)*B->BlockSize*B->BlockSize*sizeof(Complex16));
    B->b = (Complex16*)malloc((B->N)*B->BlockSize*B->BlockSize*sizeof(Complex16));
    if(B->b==NULL)
    {
      free(B->a);
      printf("malloc failed in CopyBlockTridiagonalMatrix");
      return 1;
    }
    memcpy(B->b, A.b, (B->N)*B->BlockSize*B->BlockSize*sizeof(Complex16));
  }
  return 0;
}

