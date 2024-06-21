
/* Free a tensor M itselve (i.e. the structure is not freed)
 */
void FreeCompactTensor(CompactTensorType * const M)
{
  if(M->OnLuaStack) return;
  free(M->Name);
  if(M->QComplex)
  {
    free(M->V);
  }
  else
  {
    free(M->VRe);
  }
  free(M->dimensions);
  free(M->variant);
  M->N=0;
}

void FreeCompactTensorIncludingPointer(CompactTensorType * const M)
{
  if(M->OnLuaStack) return;
  FreeCompactTensor(M);
  free(M);
}

/* Initiallize a tensor M
 *  M.QComplex, M.rank, M.dimensions, M.variant must be set
 */
uint32_t InitCompactTensor(CompactTensorType * const M)
{
  int i;
  M->Name=(char*)malloc(MAXStrLength*sizeof(char));
  if(M->Name==NULL){printf("malloc for name failed in InitCompactMatrix\n");return 1;};
  strcpy(M->Name, "Tensor");
  M->N = 1;
  for(i=0;i<M->rank;i++)
  {
    M->N *= M->dimensions[i];
  }
  if(M->QComplex)
  {
    M->V=(Complex16*)calloc(M->N,sizeof(Complex16));
    if(M->V==NULL)
    {
      printf("Malloc failed in InitCompactMatrix for complex allocation\n");
      fflush(stdout);
      return 1;
    }
  }
  else
  {
    M->VRe=(double*)calloc(M->N,sizeof(double));
    if(M->VRe==NULL)
    {
      printf("Malloc failed in InitCompactMatrix for real allocation\n");
      fflush(stdout);
      return 1;
    }
  }
  return 0;
}

int CompactTensorToString(CompactTensorType * const M, DynamicString *Ad)
{
  int i,j;
  int *indices = (int*)calloc(M->rank,sizeof(int));
  int spaces=0;
  bool newline=false;
  *Ad = InitDynamicString();
  for(i=0;i<M->N;i++)
  {
    for(j=M->rank-1;j>=0;j--)
    {
      if(indices[j]!=0)break;
      spaces++;
      sprintfAdd(Ad," { ");
    }
    if(M->QComplex==1)
    {
      sprintfAdd(Ad,"%24.14E",M->V[i].Re);
      if(M->V[i].Im<0)
      {
        sprintfAdd(Ad," - ");
        sprintfAdd(Ad,"%24.14E I",Abs(M->V[i].Im));
      }
      if(M->V[i].Im>0)
      {
        sprintfAdd(Ad," + ");
        sprintfAdd(Ad,"%24.14E I",Abs(M->V[i].Im));
      }
    }
    else
    {
      sprintfAdd(Ad,"%24.14E",M->VRe[i]);
    }
    for(j=M->rank-1;j>=0;j--)
    {
      indices[j]++;
      if(indices[j] < M->dimensions[j]) break;
      indices[j]=0;
    }
    for(j=M->rank-1;j>0;j--)
    {
      if(indices[j]==0)
      {
        sprintfAdd(Ad," } ");
        newline=true;
        spaces--;
      }
      else
      {
        break;
      }
    }
    if(i!=M->N-1)
    {
      sprintfAdd(Ad,",");
      if(newline==true)
      {
        sprintfAdd(Ad,"\n");
        newline=false;
        for(j=0;j<spaces;j++)
        {
          sprintfAdd(Ad,"   ");
        }
      }
    }
  }
  sprintfAdd(Ad,",\n   rank = %d ",M->rank);
  sprintfAdd(Ad,",\n   name = %s ",M->Name);
  sprintfAdd(Ad,",\n   dimension = { %d ",M->dimensions[0]);
  for(i=1;i<M->rank;i++)
  {
    sprintfAdd(Ad,", %d ",M->dimensions[i]);
  }
  sprintfAdd(Ad,"} ");
  sprintfAdd(Ad,",\n   variant = { ");
  if(M->variant[0]==true)
  {
    sprintfAdd(Ad,"true ");
  }
  else
  {
    sprintfAdd(Ad,"false ");
  }
  for(i=1;i<M->rank;i++)
  {
    if(M->variant[i]==true)
    {
      sprintfAdd(Ad,", true ");
    }
    else
    {
      sprintfAdd(Ad,", false ");
    }
  }
  sprintfAdd(Ad,"}");
  sprintfAdd(Ad,"\n }\n");
  return 0;
}

int CompactTensorToString(CompactTensorType * const M, char **A)
{
  DynamicString Ad;
  if(CompactTensorToString(M, &Ad)==1) return 1;
  *A = Ad.str;
  return 0;
}

void MakeCompactComplexTensorReal(CompactTensorType * M)
{
  if(M->QComplex==0) return;
  int i;
  for(i=0;i<M->N;i++)
  {
    M->VRe[i] = M->V[i].Re;
  }
  M->VRe=(double*)realloc(M->V,(M->N*sizeof(double)));
  M->QComplex=0;
}

int CompactTensorCopy(CompactTensorType * A, CompactTensorType * copyA)
{
  copyA->Name=(char*)malloc((strlen(A->Name)+1)*sizeof(char));
  strcpy(copyA->Name, A->Name);
  copyA->rank=A->rank;
  copyA->dimensions=(uint32_t*)malloc(A->rank*sizeof(uint32_t));
  memcpy(copyA->dimensions, A->dimensions, A->rank*sizeof(uint32_t));
  copyA->variant=(bool*)malloc(A->rank*sizeof(bool));
  memcpy(copyA->variant, A->variant, A->rank*sizeof(bool));
  copyA->N = A->N;
  copyA->QComplex=A->QComplex;
  copyA->OnLuaStack = A->OnLuaStack;
  if(A->QComplex==1)
  {
    copyA->V=(Complex16*)malloc(A->N*sizeof(Complex16));
    memcpy(copyA->V, A->V, A->N*sizeof(Complex16));
  }
  else
  {
    copyA->VRe=(double*)malloc(A->N*sizeof(double));
    memcpy(copyA->VRe, A->VRe, A->N*sizeof(double));
  }
  return 0;
}
