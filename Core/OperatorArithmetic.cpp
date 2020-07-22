/* format of a lader operator
 * stored as uint16_t
 *  if(o[i]&CrAnMask) then creation operator
 * !if(o[i]&CrAnMask) then annihilation operator
 *  if((o[i]&NotCrAnMask)<NFermion) then fermion
 *  if((o[i]&NotCrAnMask)>=NFermion) then boson
 *  if((o[i]&NotCrAnMask)>=NFermion+NBoson) then error
 */

/* Commutation relations
*  [A,B]- == AB-BA (for bosons)
*  [A,B]+ == AB+BA (for fermions)
*  bosons:
*  [bi ,bj ]- = 0    : bi  bj  =  bj  bi
*  [bi+,bj+]- = 0    : bi+ bj+ =  bj+ bi+
*  [bi ,bj+]- = dij  : bi  bj+ =  bj+ bi  + dij
*  [bi+,bj ]- =-dij  : bi+ bj  =  bj  bi+ - dij
*  fermions:
*  [fi ,fj ]+ = 0    : fi  fj  =- fj  fi
*  [fi+,fj+]+ = 0    : fi+ fj+ =- fj+ fi+
*  [fi ,fj+]+ = dij  : fi  fj+ =- fj+ fi  + dij
*  [fi+,fj ]+ = dij  : fi+ fj  =- fj  fi+ + dij
*  mixed
*  [f  ,b  ]- = 0
*/






/*  Remove the values of an operator from memory
*   On entry H contains the definitions and allocation of memory
*   On exit all values are freed
*/
uint32_t FreeOperator(OperatorType * const H)
{
  uint32_t i;
  for(i=0;i<H->Res.N;i++)
  {
    free(H->Res.Determinant[i]);
  }
  if(H->Res.N!=0)
  {
    free(H->Res.Min);
    free(H->Res.Max);
    free(H->Res.Determinant);
  }
  if ((H->l!=NULL)&&(H->l[0].N)) (H->l[0].QComplex) ? free(H->l[0].V) : free(H->l[0].VRe);
  for (i = 1; i <= H->MaxLength; i++)
  {
    free(H->l[i].KeyList);
    if (H->l[i].NAloc)
    {
      free(H->l[i].o);
      free(H->l[i].Next);
      (H->l[i].QComplex) ? free(H->l[i].V) : free(H->l[i].VRe);
    }
  }
  free(H->l);
  strcpy(H->Name, "Free Operator");
  H->QComplex = 0;
  H->NFermion = 0;
  H->NBoson = 0;
  H->NWord = 0;
  H->MaxLength = 0;
  H->Res.N=0;
  return 0;
}

/*  Initialize the Operator structure
 *  Allocate space for O.l
 */
uint32_t InitOperator(OperatorType * const O)
{
  strcpy(O->Name, "Operator");
  O->Name[MAXStrLength-1]='\0';
  O->QComplex = 0;
  O->NFermion = 0;
  O->NBoson = 0;
  O->NWord = 0;
  O->MaxLength = 0;
  O->Res.N=0;
  if ((O->l = (OperatorFixedLengthType*)malloc(sizeof(OperatorFixedLengthType))) == NULL)
  {
    printf("Could not allocate OperatorFixedLengthType in InitOperator\n");
    return 1;
  }
  O->l[0].N = 0;
  O->l[0].NBitsKey = 0;
  O->l[0].NAloc = 0;
  O->l[0].QComplex = 0;
  return 0;
}

/*  Initialize the Operator structure
 *  Allocate space for O.l
 */
uint32_t InitNamedOperator(OperatorType * const O, const char * const name)
{
  int ret = InitOperator(O);
  strcpy(O->Name, name);
  return ret;
}

void OperatorSetName(OperatorType * const O, const char * const name)
{
  strncpy(O->Name, name, MAXStrLength-1);
}

/*  Initialize the Operator structure
 *  to apply a Operator product with 
 *  an arbitrary number of operators
 */

void InitOperatorForMult(OperatorType *prod, int count, const OperatorType *l, ...)
{
  va_list ap;
  OperatorType *r;
  int j;    
  InitOperator(prod);
  int MaxLength = l->MaxLength;
  prod->QComplex = l->QComplex;
  prod->NFermion = l->NFermion;
  prod->NBoson = l->NBoson;
  prod->NWord = l->NWord;
  prod->Res.N = 0;
  
  va_start(ap, l);
  for(j=1; j<count; j++)
  {
    r = va_arg(ap, OperatorType*);
    MaxLength += r->MaxLength;
  }
  for (j = 1; j <= MaxLength; j++)
  {
    InitOperatorFixedLength(prod, j, 8);
  }
  va_end(ap);
}

/* each operator of fixed length can be either complex or real
 * this is stored in O.l[ilenght].QComplex
 * If all operators with non-zero value of O.l[ilength].N are real    O.QComplex=0;
 * If all operators with non-zero value of O.l[ilength].N are complex O.QComplex=1;
 * If there is a mixture of real and complex values                   O.QComplex=2;
 */
uint32_t SetQComplexOperator(OperatorType * const O)
{
  uint32_t il, nmax;
  O->QComplex = 0;
  nmax = 0;
  for (il = 0; il <= O->MaxLength; il++)
  {
    if (O->l[il].N != 0)
    {
      nmax++;
      if (O->l[il].QComplex == 1) O->QComplex++;
    }
  }
  if (nmax == 0)
  {
  }
  else if (O->QComplex == nmax)
  {
    O->QComplex = 1;
  }
  else if (O->QComplex != 0)
  {
    O->QComplex = 2;
  }
  return 0;
}

/*  Turn a real operator of a given length into a complex operator
 *  Reallocate the values VRe and allocate V.Im
 *  set V.Im to zero
 */
uint32_t MakeRealOperatorFixedLengthComplex(OperatorType * const O, const uint32_t length)
{
  uint32_t i;
  if (O->MaxLength < length)
  {
    printf
      ("Warning in MakeRealOperatorFixedLengthComplex:\nTried to make a real operator of length %i complex, but Maxlength is %i\n",
       length, O->MaxLength);
    return 2;
  }
  if (O->l[length].QComplex)
    return 0;
  if (O->l[length].NAloc == 0)
  {
    O->l[length].QComplex = 1;
    SetQComplexOperator(O);
    return 0;
  }
  if ((O->l[length].V = (Complex16*)malloc(O->l[length].NAloc * sizeof(Complex16))) == NULL)
  {
    printf("Cloud not allocate memory for O-l[%i].V in MakeRealOperatorFixedLengthComplex\n", length);
    return 1;
  }
  for (i = 0; i < O->l[length].N; i++)
  {
    O->l[length].V[i].Re = O->l[length].VRe[i];
    O->l[length].V[i].Im = 0;
  }
  free(O->l[length].VRe);
  O->l[length].QComplex = 1;
  SetQComplexOperator(O);
  return 0;
}

/*  Initiallize an operator of Fixed length within the full Operator O
 *  Initiallize also all operators of length smaller then l, which were not yet initiallized in O
 *  Set NBitsKey for the operator with the given length
 *  Shorter operators initiallized in addition will have NBitsKey=0
 *  NBitsKey=0 means there is one key entry, which all strings will be mapped to
 *  The indexed list will then become a normal linked list
 */
uint32_t InitOperatorFixedLength(OperatorType * const O, const uint32_t length, const uint32_t NBitsKey)
{
  uint32_t l;
  if (length == 0) return 0;
  if (length > O->MaxLength)
  {                             /* allongate the vector O->l */
    if ((O->l = (OperatorFixedLengthType*)realloc(O->l, (length + 1) * sizeof(OperatorFixedLengthType))) == NULL)
    {
      printf("could not allocate memory for OperatorFixedLengthType of length %i in InitOperatorFixedLength\n", length);
      return 1;
    }
    for (l = O->MaxLength + 1; l < length; l++)
    {
      O->l[l].N = 0;
      O->l[l].NBitsKey = 0;
      O->l[l].NAloc = 0;
      O->l[l].QComplex = 0;
      if ((O->l[l].KeyList = (uint32_t*)calloc(1, sizeof(uint32_t))) == NULL) return 1;
    }
    if ((O->l[length].KeyList = (uint32_t*)calloc((uint32_t) 1 << NBitsKey, sizeof(uint32_t))) == NULL)
    {
      printf("could not allocate memory for O->[length=%i].Keylist in InitOperatorFixedLength\n", length);
      return 1;
    }
    O->l[length].NBitsKey = NBitsKey;
    O->l[length].N = 0;
    O->l[length].NAloc = 0;
    O->l[length].QComplex = 0;
    O->MaxLength = length;
  }
  else
  {                             /* test for double initialization */
    if (O->l[length].NAloc)
    {                           /* Initiallizing an operator with already values stored; free this memory */
      free(O->l[length].o);
      free(O->l[length].Next);
      (O->l[length].QComplex) ? free(O->l[length].V) : free(O->l[length].VRe);
      O->l[length].N = 0;
      O->l[length].NAloc = 0;
      O->l[length].QComplex = 0;
    }
  }
  free(O->l[length].KeyList);
  if ((O->l[length].KeyList = (uint32_t*)calloc((uint32_t) 1 << NBitsKey, sizeof(uint32_t))) == NULL)
  {
    printf("could not allocate memory for O->[length=%i].Keylist in InitOperatorFixedLength\n", length);
    return 1;
  }
  O->l[length].NBitsKey = NBitsKey;
  SetQComplexOperator(O);
  return 0;
}

/*  Initialize operator O with the settings for NBitsKey and Maxlength of Oset
 *  NAloc, N, and QComplex are set to zero
 */
uint32_t InitOperatorCopySettings(OperatorType * const O, const OperatorType Oset)
{
  uint32_t i,length;
  if (InitOperator(O))
  {
    printf("could not Initiallize Operator in InitOperatorCopySettings\n");
    return 1;
  }
  strncpy(O->Name, Oset.Name, MAXStrLength);
  O->Name[MAXStrLength-1]='\0';
  O->NFermion = Oset.NFermion;
  O->NBoson = Oset.NBoson;
  O->NWord = Oset.NWord;
  O->Res.N = Oset.Res.N;
  O->Res.Min = (uint32_t*)malloc(Oset.Res.N*sizeof(uint32_t));
  O->Res.Max = (uint32_t*)malloc(Oset.Res.N*sizeof(uint32_t));
  O->Res.Determinant = (uint8_t**)malloc(Oset.Res.N*sizeof(uint8_t*));
  for(i=0;i<Oset.Res.N;i++)
  {
    O->Res.Min[i]=Oset.Res.Min[i];
    O->Res.Max[i]=Oset.Res.Max[i];
    O->Res.Determinant[i]= (uint8_t*)malloc(Oset.NWord*sizeof(uint8_t));
    memcpy(O->Res.Determinant[i], Oset.Res.Determinant[i], Oset.NWord);
  }
  for (length = 1; length < Oset.MaxLength; length++)
  {
    if (Oset.l[length].NBitsKey)
      if (InitOperatorFixedLength(O, length, Oset.l[length].NBitsKey))
      {
        printf("could not Initiallize InitOperatorFixedLength of length=%i in InitOperatorCopySettings\n", length);
        return 1;
      }
  }
  if (InitOperatorFixedLength(O, Oset.MaxLength, Oset.l[Oset.MaxLength].NBitsKey))
  {
    printf("could not Initiallize InitOperatorFixedLength of length=O.MaxLength=%i in InitOperatorCopySettings\n", Oset.MaxLength);
    return 1;
  }
  return 0;
}

/*  add a single lader operator string to the operator O
 *  O at the given length must be real
 */
uint32_t RealOperatorAddLader(OperatorType * const O, const uint16_t * const oadd, const uint32_t length, const double VRe)
{
  void *temp;
  uint32_t NextIndex, ThereIsANext, FirstAlloc = 0, i, j;
  if (length == 0)
  {
    if (O->l[0].N == 0)
    {
      if (O->l[0].NAloc == 0)
      {
        if ((O->l[0].VRe = (double*)malloc(sizeof(double))) == NULL)
        {
          printf("could not allocate memory for O->l[0].VRe in RealOperatorAddLader\n");
          return 1;
        }
        O->l[0].NAloc = 1;
      }
      O->l[0].N = 1;
      O->l[0].QComplex = 0;
      O->l[0].VRe[0] = VRe;
      SetQComplexOperator(O);
      return 0;
    }
    O->l[0].VRe[0] += VRe;
    return 0;
  }
  /* test if the key that belongs to oadd already excists */
  if ((O->l[length].KeyList[hash(oadd, length * 2, O->l[length].NBitsKey)]) == 0)
  {                             /* add a new ellement to O with a unique key */
    /* Test if one has to reallocate O->l[length].(o,Next,VRe) in order to store the new element */
    if (O->l[length].N >= O->l[length].NAloc)
    {                           /* reallocate O */
      if (O->l[length].NAloc == 0)
        FirstAlloc = 1;
      O->l[length].NAloc += O->NFermion + O->NBoson;
      temp = O->l[length].o;
      if (FirstAlloc)
      {
        O->l[length].o = (uint16_t*)malloc(length * O->l[length].NAloc * sizeof(uint16_t));
        SetQComplexOperator(O);
      }
      else
      {
        O->l[length].o = (uint16_t*)realloc(O->l[length].o, length * O->l[length].NAloc * sizeof(uint16_t));
      }
      if (O->l[length].o == NULL)
      {                         /* realloc failed */
        printf("alloc failed OperatorAddLader 01\n");
        fflush(stdout);
        /* restore O and return 1 */
        O->l[length].o = (uint16_t*)temp;
        O->l[length].NAloc -= O->NFermion + O->NBoson;
        return 1;
      }
      temp = O->l[length].Next;
      if (FirstAlloc)
      {
        O->l[length].Next = (uint32_t*)malloc(O->l[length].NAloc * sizeof(uint32_t));
      }
      else
      {
        O->l[length].Next = (uint32_t*)realloc(O->l[length].Next, O->l[length].NAloc * sizeof(uint32_t));
      }
      if (O->l[length].Next == NULL)
      {                         /* realloc failed */
        printf("alloc failed OperatorAddLader 02\n");
        fflush(stdout);
        /* restore O and return 1 */
        O->l[length].Next = (uint32_t*)temp;
        O->l[length].NAloc -= O->NFermion + O->NBoson;
        O->l[length].o = (uint16_t*)realloc(O->l[length].o, length * O->l[length].NAloc * sizeof(uint16_t));
        return 1;
      }
      temp = O->l[length].VRe;
      if (FirstAlloc)
      {
        O->l[length].VRe = (double*)malloc(O->l[length].NAloc * sizeof(double));
      }
      else
      {
        O->l[length].VRe = (double*)realloc(O->l[length].VRe, O->l[length].NAloc * sizeof(double));
      }
      if (O->l[length].VRe == NULL)
      {                         /* realloc failed */
        printf("alloc failed OperatorAddLader 03\n");
        fflush(stdout);
        /* restore O and return 1 */
        O->l[length].VRe = (double*)temp;
        O->l[length].NAloc -= O->NFermion + O->NBoson;
        O->l[length].o = (uint16_t*)realloc(O->l[length].o, length * O->l[length].NAloc * sizeof(uint16_t));
        O->l[length].Next = (uint32_t*)realloc(O->l[length].Next, O->l[length].NAloc * sizeof(uint32_t));
        return 1;
      }
    }                           /* finished reallocating O if needed */
    for (i = 0; i < length; i++)
    {
      O->l[length].o[O->l[length].N * length + i] = oadd[i];
    }
    O->l[length].VRe[O->l[length].N] = VRe;
    O->l[length].Next[O->l[length].N] = 0;
    O->l[length].N++;
    O->l[length].KeyList[hash(oadd, length * 2, O->l[length].NBitsKey)] = O->l[length].N;
  }
  else
  {                             /* the key is already taken, test if this is due to a collision or if oadd is already in the opperator */
    ThereIsANext = 1;
    NextIndex = O->l[length].KeyList[hash(oadd, length * 2, O->l[length].NBitsKey)] - 1;
    while (ThereIsANext)
    {
      j = 1;
      for (i = 0; i < length; i++)
      {
        if (O->l[length].o[NextIndex * length + i] != oadd[i])
          j = 0;
      }
      if (j == 1)
      {                         /* lader operator string already in the opperator, add and leave loop */
        ThereIsANext = 0;
        O->l[length].VRe[NextIndex] += VRe;
      }
      else
      {                         /* oadd not equal to the O->l[length].o found which has the same index in the keylist, if there is a next continue, otherwise add the element */
        if (!(O->l[length].Next[NextIndex]))
        {                       /* add a new ellemnt to O with already used key */
          /* Test if one has to reallocate O->l[length].(o,Next,VRe) in order to store the new element */
          if (O->l[length].N >= O->l[length].NAloc)
          {
            O->l[length].NAloc += O->NFermion + O->NBoson;
            temp = O->l[length].o;
            O->l[length].o = (uint16_t*)realloc(O->l[length].o, length * O->l[length].NAloc * sizeof(uint16_t));
            if (O->l[length].o == NULL)
            {                   /* realloc failed */
              printf("alloc failed OperatorAddLader 11\n");
              fflush(stdout);
              /* restore O and return 1 */
              O->l[length].o = (uint16_t*)temp;
              O->l[length].NAloc -= O->NFermion + O->NBoson;
              return 1;
            }
            temp = O->l[length].Next;
            O->l[length].Next = (uint32_t*)realloc(O->l[length].Next, O->l[length].NAloc * sizeof(uint32_t));
            if (O->l[length].Next == NULL)
            {                   /* realloc failed */
              printf("alloc failed OperatorAddLader 12\n");
              fflush(stdout);
              /* restore O and return 1 */
              O->l[length].Next = (uint32_t*)temp;
              O->l[length].NAloc -= O->NFermion + O->NBoson;
              O->l[length].o = (uint16_t*)realloc(O->l[length].o, length * O->l[length].NAloc * sizeof(uint16_t));
              return 1;
            }
            temp = O->l[length].VRe;
            O->l[length].VRe = (double*)realloc(O->l[length].VRe, O->l[length].NAloc * sizeof(double));
            if (O->l[length].VRe == NULL)
            {                   /* realloc failed */
              printf("alloc failed OperatorAddLader 13\n");
              fflush(stdout);
              /* restore O and return 1 */
              O->l[length].VRe = (double*)temp;
              O->l[length].NAloc -= O->NFermion + O->NBoson;
              O->l[length].o = (uint16_t*)realloc(O->l[length].o, length * O->l[length].NAloc * sizeof(uint16_t));
              O->l[length].Next = (uint32_t*)realloc(O->l[length].Next, O->l[length].NAloc * sizeof(uint32_t));
              return 1;
            }
          }
          for (i = 0; i < length; i++)
          {
            O->l[length].o[O->l[length].N * length + i] = oadd[i];
          }
          O->l[length].VRe[O->l[length].N] = VRe;
          O->l[length].Next[O->l[length].N] = 0;
          O->l[length].N++;
          O->l[length].Next[NextIndex] = O->l[length].N;
          ThereIsANext = 0;
        }
        else
        {                       /* The key is in use and there are more opperators with this key */
          NextIndex = O->l[length].Next[NextIndex] - 1;
        }
      }
    }
  }
  return 0;
}

/*  add a single lader operator string to the operator O, after normal ordering O
 *  O at the given length must be real
 *  iff due to normal ordering operators of different lengths arize they are added
 */
uint32_t RealOperatorAddLaderNormalOrder(OperatorType * const O, const uint16_t * const oadd, const uint32_t length, const double VRe, uint32_t(*QDetNotOrdered) (uint16_t, uint16_t, uint32_t))
{
  uint16_t oaddordered[length + 1];
  uint16_t oaddnew[length + 1];
  uint16_t temp;
  uint32_t i, j, k, NeedToSwap;
  double phase = 1;
  if(length>0) oaddordered[0] = oadd[0];
  for (i = 1; i < length; i++)
  {
    oaddordered[i] = oadd[i];
    NeedToSwap = 1;
    j = i;
    while (NeedToSwap && j)
    {
      switch (QDetNotOrdered(oaddordered[j], oaddordered[j - 1], O->NFermion))
      {
      case 0:
        /* A is left of B
         */
        NeedToSwap = 0;
        break;
      case 1:
        /*  A should be right of B and A is a boson and B is a fermion
         *  A should be right of B and A is a fermion and B is a boson
         *  A should be right of B and both are bosons and A annihilates and B annihilates
         *  A should be right of B and both are bosons and A creates and B creates
         */
        temp = oaddordered[j - 1];
        oaddordered[j - 1] = oaddordered[j];
        oaddordered[j] = temp;
        break;
      case 2:
        /* A should be right of B and both are fermions and A annihilates and B creates
         * A should be right of B and both are fermions and A creates and B annihilates
         *  [fi ,fj+]+ = dij  : fi  fj+ =- fj+ fi  + dij
         *  [fi+,fj ]+ = dij  : fi+ fj  =- fj  fi+ + dij
         */
        temp = oaddordered[j - 1];
        oaddordered[j - 1] = oaddordered[j];
        oaddordered[j] = temp;
        if ((oaddordered[j - 1] & NotCrAnMask) == (oaddordered[j] & NotCrAnMask))
        {
          for (k = 0; k < j - 1; k++)
            oaddnew[k] = oaddordered[k];
          for (k = j + 1; k <= i; k++)
            oaddnew[k - 2] = oaddordered[k];
          for (k = i + 1; k < length; k++)
            oaddnew[k - 2] = oadd[k];
          if (O->l[length - 2].QComplex)
          {
            if (ComplexOperatorAddLaderNormalOrder(O, oaddnew, length - 2, phase * VRe, 0.0, QDetNotOrdered))
            {
              printf
                ("ComplexOperatorAddLaderNormalOrder with operator of length %i failed in RealOperatorAddLaderNormalOrder case 2\n",
                 length - 2);
              fflush(stdout);
              return 1;
            }
          }
          else
          {
            if (RealOperatorAddLaderNormalOrder(O, oaddnew, length - 2, phase * VRe, QDetNotOrdered))
            {
              printf
                ("RealOperatorAddLaderNormalOrder with operator of length %i failed in RealOperatorAddLaderNormalOrder case 2\n",
                 length - 2);
              fflush(stdout);
              return 1;
            }
          }
        }
        phase *= -1;
        break;
      case 3:
        /* A should be right of B and both are fermions and A annihilates and B annihilates but A!=B
         * A should be right of B and both are fermions and A creates and B creates but A!=B
         */
        temp = oaddordered[j - 1];
        oaddordered[j - 1] = oaddordered[j];
        oaddordered[j] = temp;
        phase *= -1;
        break;
      case 4:
        /* A should be right of B and both are bosons and A annihilates and B creates
         *  [bi ,bj+]- = dij  : bi  bj+ =  bj+ bi  + dij
         */
        temp = oaddordered[j - 1];
        oaddordered[j - 1] = oaddordered[j];
        oaddordered[j] = temp;
        if ((oaddordered[j - 1] & NotCrAnMask) == (oaddordered[j] & NotCrAnMask))
        {
          for (k = 0; k < j - 1; k++)
            oaddnew[k] = oaddordered[k];
          for (k = j + 1; k <= i; k++)
            oaddnew[k - 2] = oaddordered[k];
          for (k = i + 1; k < length; k++)
            oaddnew[k - 2] = oadd[k];
          if (O->l[length - 2].QComplex)
          {
            if (ComplexOperatorAddLaderNormalOrder(O, oaddnew, length - 2, phase * VRe, 0.0, QDetNotOrdered))
            {
              printf
                ("ComplexOperatorAddLaderNormalOrder with operator of length %i failed in RealOperatorAddLaderNormalOrder case 4\n",
                 length - 2);
              fflush(stdout);
              return 1;
            }
          }
          else
          {
            if (RealOperatorAddLaderNormalOrder(O, oaddnew, length - 2, phase * VRe, QDetNotOrdered))
            {
              printf
                ("RealOperatorAddLaderNormalOrder with operator of length %i failed in RealOperatorAddLaderNormalOrder case 4\n",
                 length - 2);
              fflush(stdout);
              return 1;
            }
          }
        }
        break;
      case 5:
        /* A should be right of B and both are bosons and A creates and B annihilates
         *  [bi+,bj ]- =-dij  : bi+ bj  =  bj  bi+ - dij
         */
        temp = oaddordered[j - 1];
        oaddordered[j - 1] = oaddordered[j];
        oaddordered[j] = temp;
        if ((oaddordered[j - 1] & NotCrAnMask) == (oaddordered[j] & NotCrAnMask))
        {
          for (k = 0; k < j - 1; k++)
            oaddnew[k] = oaddordered[k];
          for (k = j + 1; k <= i; k++)
            oaddnew[k - 2] = oaddordered[k];
          for (k = i + 1; k < length; k++)
            oaddnew[k - 2] = oadd[k];
          if (O->l[length - 2].QComplex)
          {
            if (ComplexOperatorAddLaderNormalOrder(O, oaddnew, length - 2, -phase * VRe, 0.0, QDetNotOrdered))
            {
              printf
                ("ComplexOperatorAddLaderNormalOrder with operator of length %i failed in RealOperatorAddLaderNormalOrder case 5\n",
                 length - 2);
              fflush(stdout);
              return 1;
            }
          }
          else
          {
            if (RealOperatorAddLaderNormalOrder(O, oaddnew, length - 2, -phase * VRe, QDetNotOrdered))
            {
              printf
                ("RealOperatorAddLaderNormalOrder with operator of length %i failed in RealOperatorAddLaderNormalOrder case 5\n",
                 length - 2);
              fflush(stdout);
              return 1;
            }
          }
        }
        break;
      case 6:
        /* A and B are fermions and A==B (A and B both create or A and B both annihilate)
         */
        return 0;
        break;
      }
      j--;
    }
  }
  return RealOperatorAddLader(O, oaddordered, length, phase * VRe);
}

/*
 *
 */
double RealOperatorPrefactorLader(const OperatorType * const O, const uint16_t * const oadd, const uint32_t length)
{
  double ret=0;
  uint32_t NextIndex, ThereIsANext,i,j;
  /* test if the key that belongs to oadd already excists */
  if ((O->l[length].KeyList[hash(oadd, length * 2, O->l[length].NBitsKey)]) == 0)
  {
    ret = 0;
  }
  else
  {                             /* the key is already taken, test if this is due to a collision or if oadd is already in the opperator */
    ThereIsANext = 1;
    NextIndex = O->l[length].KeyList[hash(oadd, length * 2, O->l[length].NBitsKey)] - 1;
    while (ThereIsANext)
    {
      j = 1;
      for (i = 0; i < length; i++)
      {
        if (O->l[length].o[NextIndex * length + i] != oadd[i])
          j = 0;
      }
      if (j == 1)
      {                         /* lader operator string already in the opperator, set ret and leave the loop */
        ThereIsANext = 0;
        ret = O->l[length].VRe[NextIndex];
      }
      else
      {                         /* oadd not equal to the O->l[length].o found which has the same index in the keylist, if there is a next continue, otherwise return 0 */
        if (!(O->l[length].Next[NextIndex]))
        {                       /* no more elements with same key and element not foud - return 0 and leave the loop*/
          ret = 0;
          ThereIsANext = 0;
        }
        else
        {                       /* The key is in use and there are more opperators with this key */
          NextIndex = O->l[length].Next[NextIndex] - 1;
        }
      }
    }
  }
  return ret;
}


/*  add a single lader operator string to the operator O
 *  O at the given length must be complex or it must be the first element added.
 */
uint32_t ComplexOperatorAddLader(OperatorType * const O, const uint16_t * const oadd, const uint32_t length, const double VRe, const double VIm)
{
  void *temp;
  uint32_t NextIndex, ThereIsANext, FirstAlloc = 0, i, j;
  if (length == 0)
  {
    if (O->l[0].N == 0)
    {
      if (O->l[0].NAloc == 0)
      {
        if ((O->l[0].V = (Complex16*)malloc(sizeof(Complex16))) == NULL)
          return 1;
        O->l[0].NAloc = 1;
      }
      O->l[0].N = 1;
      O->l[0].QComplex = 1;
      O->l[0].V[0].Re = VRe;
      O->l[0].V[0].Im = VIm;
      return 0;
    }
    O->l[0].V[0].Re += VRe;
    O->l[0].V[0].Im += VIm;
    return 0;
  }
  /* test if the key that belongs to oadd already excists */
  if ((O->l[length].KeyList[hash(oadd, length * 2, O->l[length].NBitsKey)]) == 0)
  {                             /* add a new ellement to O with a unique key */
    /* Test if one has to reallocate O->l[length].(o,Next,VRe) in order to store the new element */
    if (O->l[length].N >= O->l[length].NAloc)
    {                           /* reallocate O */
      if (O->l[length].NAloc == 0)
        FirstAlloc = 1;
      O->l[length].NAloc += O->NFermion + O->NBoson;
      temp = O->l[length].o;
      if (FirstAlloc)
      {
        O->l[length].o = (uint16_t*)malloc(length * O->l[length].NAloc * sizeof(uint16_t));
      }
      else
      {
        O->l[length].o = (uint16_t*)realloc(O->l[length].o, length * O->l[length].NAloc * sizeof(uint16_t));
      }
      if (O->l[length].o == NULL)
      {                         /* realloc failed */
        printf("alloc failed OperatorAddLader 01\n");
        fflush(stdout);
        /* restore O and return 1 */
        O->l[length].o = (uint16_t*)temp;
        O->l[length].NAloc -= O->NFermion + O->NBoson;
        return (1);
      }
      temp = O->l[length].Next;
      if (FirstAlloc)
      {
        O->l[length].Next = (uint32_t*)malloc(O->l[length].NAloc * sizeof(uint32_t));
      }
      else
      {
        O->l[length].Next = (uint32_t*)realloc(O->l[length].Next, O->l[length].NAloc * sizeof(uint32_t));
      }
      if (O->l[length].Next == NULL)
      {                         /* realloc failed */
        printf("alloc failed OperatorAddLader 02\n");
        fflush(stdout);
        /* restore O and return 1 */
        O->l[length].Next = (uint32_t*)temp;
        O->l[length].NAloc -= O->NFermion + O->NBoson;
        O->l[length].o = (uint16_t*)realloc(O->l[length].o, length * O->l[length].NAloc * sizeof(uint16_t));
        return (1);
      }
      temp = O->l[length].V;
      if (FirstAlloc)
      {
        O->l[length].V = (Complex16*)malloc(O->l[length].NAloc * sizeof(Complex16));
        O->l[length].QComplex = 1;
        O->l[length].N = 1;
        SetQComplexOperator(O);
        O->l[length].N = 0;
      }
      else
      {
        O->l[length].V = (Complex16*)realloc(O->l[length].V, O->l[length].NAloc * sizeof(Complex16));
      }
      if (O->l[length].V == NULL)
      {                         /* realloc failed */
        printf("alloc failed OperatorAddLader 03\n");
        fflush(stdout);
        /* restore O and return 1 */
        O->l[length].V = (Complex16*)temp;
        O->l[length].NAloc -= O->NFermion + O->NBoson;
        O->l[length].o = (uint16_t*)realloc(O->l[length].o, length * O->l[length].NAloc * sizeof(uint16_t));
        O->l[length].Next = (uint32_t*)realloc(O->l[length].Next, O->l[length].NAloc * sizeof(uint32_t));
        return (1);
      }
    }                           /* finished reallocating O if needed */
    for (i = 0; i < length; i++)
    {
      O->l[length].o[O->l[length].N * length + i] = oadd[i];
    }
    O->l[length].V[O->l[length].N].Re = VRe;
    O->l[length].V[O->l[length].N].Im = VIm;
    O->l[length].Next[O->l[length].N] = 0;
    O->l[length].N++;
    O->l[length].KeyList[hash(oadd, length * 2, O->l[length].NBitsKey)] = O->l[length].N;
  }
  else
  {                             /* the key is already taken, test if this is due to a collision or if oadd is already in the opperator */
    ThereIsANext = 1;
    NextIndex = O->l[length].KeyList[hash(oadd, length * 2, O->l[length].NBitsKey)] - 1;
    while (ThereIsANext)
    {
      j = 1;
      for (i = 0; i < length; i++)
      {
        if (O->l[length].o[NextIndex * length + i] != oadd[i])
          j = 0;
      }
      if (j == 1)
      {                         /* lader operator string already in the opperator, add and leave loop */
        ThereIsANext = 0;
        O->l[length].V[NextIndex].Re += VRe;
        O->l[length].V[NextIndex].Im += VIm;
      }
      else
      {                         /* oadd not equal to the O->l[length].o found which has the same index in the keylist, if there is a next continue, otherwise add the element */
        if (!(O->l[length].Next[NextIndex]))
        {                       /* add a new ellemnt to O with already used key */
          /* Test if one has to reallocate O->l[length].(o,Next,VRe) in order to store the new element */
          if (O->l[length].N >= O->l[length].NAloc)
          {
            O->l[length].NAloc += O->NFermion + O->NBoson;
            temp = O->l[length].o;
            O->l[length].o = (uint16_t*)realloc(O->l[length].o, length * O->l[length].NAloc * sizeof(uint16_t));
            if (O->l[length].o == NULL)
            {                   /* realloc failed */
              printf("alloc failed OperatorAddLader 11\n");
              fflush(stdout);
              /* restore O and return 1 */
              O->l[length].o = (uint16_t*)temp;
              O->l[length].NAloc -= O->NFermion + O->NBoson;
              return (1);
            }
            temp = O->l[length].Next;
            O->l[length].Next = (uint32_t*)realloc(O->l[length].Next, O->l[length].NAloc * sizeof(uint32_t));
            if (O->l[length].Next == NULL)
            {                   /* realloc failed */
              printf("alloc failed OperatorAddLader 12\n");
              fflush(stdout);
              /* restore O and return 1 */
              O->l[length].Next = (uint32_t*)temp;
              O->l[length].NAloc -= O->NFermion + O->NBoson;
              O->l[length].o = (uint16_t*)realloc(O->l[length].o, length * O->l[length].NAloc * sizeof(uint16_t));
              return (1);
            }
            temp = O->l[length].V;
            O->l[length].V = (Complex16*)realloc(O->l[length].V, O->l[length].NAloc * sizeof(Complex16));
            if (O->l[length].V == NULL)
            {                   /* realloc failed */
              printf("alloc failed OperatorAddLader 13\n");
              fflush(stdout);
              /* restore O and return 1 */
              O->l[length].V = (Complex16*)temp;
              O->l[length].NAloc -= O->NFermion + O->NBoson;
              O->l[length].o = (uint16_t*)realloc(O->l[length].o, length * O->l[length].NAloc * sizeof(uint16_t));
              O->l[length].Next = (uint32_t*)realloc(O->l[length].Next, O->l[length].NAloc * sizeof(uint32_t));
              return (1);
            }
          }
          for (i = 0; i < length; i++)
          {
            O->l[length].o[O->l[length].N * length + i] = oadd[i];
          }
          O->l[length].V[O->l[length].N].Re = VRe;
          O->l[length].V[O->l[length].N].Im = VIm;
          O->l[length].Next[O->l[length].N] = 0;
          O->l[length].N++;
          O->l[length].Next[NextIndex] = O->l[length].N;
          ThereIsANext = 0;
        }
        else
        {                       /* The key is in use and there are more opperators with this key */
          NextIndex = O->l[length].Next[NextIndex] - 1;
        }
      }
    }
  }
  return 0;
}

/*  add a single lader operator string to the operator O, after normal ordering O
 *  O at the given length must be complex
 */
uint32_t ComplexOperatorAddLaderNormalOrder(OperatorType * const O, const uint16_t * const oadd, const uint32_t length, const double VRe, const double VIm, uint32_t(*QDetNotOrdered) (uint16_t,uint16_t,uint32_t))
{
  uint16_t oaddordered[length + 1];
  uint16_t oaddnew[length + 1];
  uint16_t temp;
  uint32_t i, j, k, NeedToSwap;
  double phase = 1;
  if(length>0) oaddordered[0] = oadd[0];
  for (i = 1; i < length; i++)
  {
    oaddordered[i] = oadd[i];
    NeedToSwap = 1;
    j = i;
    while (NeedToSwap && j)
    {
      switch (QDetNotOrdered(oaddordered[j], oaddordered[j - 1], O->NFermion))
      {
      case 0:
        /* A is left of B
         */
        NeedToSwap = 0;
        break;
      case 1:
        /*  A should be right of B and A is a boson and B is a fermion
         *  A should be right of B and A is a fermion and B is a boson
         *  A should be right of B and both are bosons and A annihilates and B annihilates
         *  A should be right of B and both are bosons and A creates and B creates
         */
        temp = oaddordered[j - 1];
        oaddordered[j - 1] = oaddordered[j];
        oaddordered[j] = temp;
        break;
      case 2:
        /* A should be right of B and both are fermions and A annihilates and B creates
         * A should be right of B and both are fermions and A creates and B annihilates
         *  [fi ,fj+]+ = dij  : fi  fj+ =- fj+ fi  + dij
         *  [fi+,fj ]+ = dij  : fi+ fj  =- fj  fi+ + dij
         */
        temp = oaddordered[j - 1];
        oaddordered[j - 1] = oaddordered[j];
        oaddordered[j] = temp;
        if ((oaddordered[j - 1] & NotCrAnMask) == (oaddordered[j] & NotCrAnMask))
        {
          for (k = 0; k < j - 1; k++) oaddnew[k] = oaddordered[k];
          for (k = j + 1; k <= i; k++) oaddnew[k - 2] = oaddordered[k];
          for (k = i + 1; k < length; k++) oaddnew[k - 2] = oadd[k];
          if (O->l[length - 2].QComplex)
          {
            if (ComplexOperatorAddLaderNormalOrder(O, oaddnew, length - 2, phase * VRe, phase * VIm, QDetNotOrdered))
            {
              printf("ComplexOperatorAddLaderNormalOrder with operator of length %i failed in ComplexOperatorAddLaderNormalOrder case 2\n",length - 2);
              fflush(stdout);
              return 1;
            }
          }
          else
          {
            if (fabs(VIm) < DBL_MIN * 10)
            {
              if (RealOperatorAddLaderNormalOrder(O, oaddnew, length - 2, phase * VRe, QDetNotOrdered))
              {
                printf("RealOperatorAddLaderNormalOrder with operator of length %i failed in ComplexOperatorAddLaderNormalOrder case 2\n", length - 2);
                fflush(stdout);
                return 1;
              }
            }
            else
            {
              if (MakeRealOperatorFixedLengthComplex(O, length - 2))
              {
                printf("Failed to make Real operator Complex In ComplexOperatorAddLaderNormalOrder case 2\n");
                fflush(stdout);
                return 1;
              }
              SetQComplexOperator(O);
              if (ComplexOperatorAddLaderNormalOrder(O, oaddnew, length - 2, phase * VRe, phase * VIm, QDetNotOrdered))
              {
                printf("ComplexOperatorAddLaderNormalOrder with operator of length %i failed in ComplexOperatorAddLaderNormalOrder case 2\n", length - 2);
                fflush(stdout);
                return 1;
              }
            }
          }
        }
        phase *= -1;
        break;
      case 3:
        /* A should be right of B and both are fermions and A annihilates and B annihilates but A!=B
         * A should be right of B and both are fermions and A creates and B creates but A!=B
         */
        temp = oaddordered[j - 1];
        oaddordered[j - 1] = oaddordered[j];
        oaddordered[j] = temp;
        phase *= -1;
        break;
      case 4:
        /* A should be right of B and both are bosons and A annihilates and B creates
         *  [bi ,bj+]- = dij  : bi  bj+ =  bj+ bi  + dij
         */
        temp = oaddordered[j - 1];
        oaddordered[j - 1] = oaddordered[j];
        oaddordered[j] = temp;
        if ((oaddordered[j - 1] & NotCrAnMask) == (oaddordered[j] & NotCrAnMask))
        {
          for (k = 0; k < j - 1; k++)
            oaddnew[k] = oaddordered[k];
          for (k = j + 1; k <= i; k++)
            oaddnew[k - 2] = oaddordered[k];
          for (k = i + 1; k < length; k++)
            oaddnew[k - 2] = oadd[k];
          if (O->l[length - 2].QComplex)
          {
            if (ComplexOperatorAddLaderNormalOrder(O, oaddnew, length - 2, phase * VRe, phase * VIm, QDetNotOrdered))
            {
              printf("ComplexOperatorAddLaderNormalOrder with operator of length %i failed in ComplexOperatorAddLaderNormalOrder case 4\n", length - 2);
              fflush(stdout);
              return 1;
            }
          }
          else
          {
            if (fabs(VIm) < DBL_MIN * 10)
            {
              if (RealOperatorAddLaderNormalOrder(O, oaddnew, length - 2, phase * VRe, QDetNotOrdered))
              {
                printf("RealOperatorAddLaderNormalOrder with operator of length %i failed in ComplexOperatorAddLaderNormalOrder case 4\n", length - 2);
                fflush(stdout);
                return 1;
              }
            }
            else
            {
              if (MakeRealOperatorFixedLengthComplex(O, length - 2))
              {
                printf("Failed to make Real operator Complex In ComplexOperatorAddLaderNormalOrder 01\n");
                return 1;
              };
              SetQComplexOperator(O);
              if (ComplexOperatorAddLaderNormalOrder(O, oaddnew, length - 2, phase * VRe, phase * VIm, QDetNotOrdered))
              {
                printf
                  ("ComplexOperatorAddLaderNormalOrder with operator of length %i failed in ComplexOperatorAddLaderNormalOrder case 4\n",length - 2);
                fflush(stdout);
                return 1;
              }
            }
          }
        }
        break;
      case 5:
        /* A should be right of B and both are bosons and A creates and B annihilates
         *  [bi+,bj ]- =-dij  : bi+ bj  =  bj  bi+ - dij
         */
        temp = oaddordered[j - 1];
        oaddordered[j - 1] = oaddordered[j];
        oaddordered[j] = temp;
        if ((oaddordered[j - 1] & NotCrAnMask) == (oaddordered[j] & NotCrAnMask))
        {
          for (k = 0; k < j - 1; k++)  oaddnew[k] = oaddordered[k];
          for (k = j + 1; k <= i; k++) oaddnew[k - 2] = oaddordered[k];
          for (k = i + 1; k < length; k++) oaddnew[k - 2] = oadd[k];
          if (O->l[length - 2].QComplex)
          {
            if (ComplexOperatorAddLaderNormalOrder(O, oaddnew, length - 2, -phase * VRe, -phase * VIm, QDetNotOrdered))
            {
              printf
                ("ComplexOperatorAddLaderNormalOrder with operator of length %i failed in ComplexOperatorAddLaderNormalOrder case 5\n",
                 length - 2);
              fflush(stdout);
              return 1;
            }
          }
          else
          {
            if (fabs(VIm) < DBL_MIN * 10)
            {
              if (RealOperatorAddLaderNormalOrder(O, oaddnew, length - 2, -phase * VRe, QDetNotOrdered))
              {
                printf("RealOperatorAddLaderNormalOrder with operator of length %i failed in ComplexOperatorAddLaderNormalOrder case 5\n", length - 2);
                fflush(stdout);
                return 1;
              }
            }
            else
            {
              if (MakeRealOperatorFixedLengthComplex(O, length - 2))
              {
                printf("Failed to make Real operator Complex In ComplexOperatorAddLaderNormalOrder 01\n");
                return 1;
              }
              SetQComplexOperator(O);
              if (ComplexOperatorAddLaderNormalOrder(O, oaddnew, length - 2, -phase * VRe, -phase * VIm, QDetNotOrdered))
              {
                printf
                  ("ComplexOperatorAddLaderNormalOrder with operator of length %i failed in ComplexOperatorAddLaderNormalOrder case 5\n",length - 2);
                fflush(stdout);
                return 1;
              }
            }
          }
        }
        break;
      case 6:
        /* A and B are fermions and A==B (A and B both create or A and B both annihilate)
         */
        return 0;
        break;
      }
      j--;
    }
  }
  return ComplexOperatorAddLader(O, oaddordered, length, phase * VRe, phase * VIm);
}

/*
 *
 */
Complex16 ComplexOperatorPrefactorLader(const OperatorType * const O, const uint16_t * const oadd, const uint32_t length)
{
  Complex16 ret={0,0};
  uint32_t NextIndex, ThereIsANext,i,j;
  /* test if the key that belongs to oadd already excists */
  if ((O->l[length].KeyList[hash(oadd, length * 2, O->l[length].NBitsKey)]) == 0)
  {
    ret = {0,0};
  }
  else
  {                             /* the key is already taken, test if this is due to a collision or if oadd is already in the opperator */
    ThereIsANext = 1;
    NextIndex = O->l[length].KeyList[hash(oadd, length * 2, O->l[length].NBitsKey)] - 1;
    while (ThereIsANext)
    {
      j = 1;
      for (i = 0; i < length; i++)
      {
        if (O->l[length].o[NextIndex * length + i] != oadd[i]) j = 0;
      }
      if (j == 1)
      {                         /* lader operator string already in the opperator, set ret and leave the loop */
        ThereIsANext = 0;
        ret = O->l[length].V[NextIndex];
      }
      else
      {                         /* oadd not equal to the O->l[length].o found which has the same index in the keylist, if there is a next continue, otherwise return 0 */
        if (!(O->l[length].Next[NextIndex]))
        {                       /* no more elements with same key and element not foud - return 0 and leave the loop*/
          ret = {0,0};
          ThereIsANext = 0;
        }
        else
        {                       /* The key is in use and there are more opperators with this key */
          NextIndex = O->l[length].Next[NextIndex] - 1;
        }
      }
    }
  }
  return ret;
}

/*  add a single lader operator string to the operator O
 *  if VIm is zero and O.l[lenght].QComplex=0 then a real is added
 */
uint32_t OperatorAddLader(OperatorType * const O, const uint16_t * const oadd, const uint32_t length, const double VRe, const double VIm)
{
  if (O->l[length].QComplex)
  {
    return ComplexOperatorAddLader(O, oadd, length, VRe, VIm);
  }
  if (fabs(VIm) < (DBL_MIN * 10)) return RealOperatorAddLader(O, oadd, length, VRe);
  if (MakeRealOperatorFixedLengthComplex(O, length))
  {
    printf("MakeRealOperatorFixedLengthComplex at length=%i failed in OperatorAddLader\n", length);
    fflush(stdout);
    return 1;
  }
  return ComplexOperatorAddLader(O, oadd, length, VRe, VIm);
}

/*  add a single lader operator string to the operator O, after normal ordering O
 */
uint32_t OperatorAddLaderNormalOrder(OperatorType * const O, const uint16_t * const oadd, const uint32_t length,
                                     const double VRe, const double VIm, uint32_t(*QDetNotOrdered) (uint16_t, uint16_t, uint32_t))
{
  if (O->l[length].QComplex)
  {
    return ComplexOperatorAddLaderNormalOrder(O, oadd, length, VRe, VIm, QDetNotOrdered);
  }
  if (fabs(VIm) < (DBL_MIN * 10)) return RealOperatorAddLaderNormalOrder(O, oadd, length, VRe, QDetNotOrdered);
  if (MakeRealOperatorFixedLengthComplex(O, length))
  {
    printf("MakeRealOperatorFixedLengthComplex at length=%i failed in OperatorAddLaderNormalOrder\n", length);
    fflush(stdout);
    return 1;
  }
  return ComplexOperatorAddLaderNormalOrder(O, oadd, length, VRe, VIm, QDetNotOrdered);
}

/*  Ocopy <- O
 */
uint32_t OperatorFixedLengthCopy(OperatorFixedLengthType * const Ocopy, OperatorFixedLengthType const O, uint32_t const length)
{
  Ocopy->N = O.N;
  Ocopy->NBitsKey = O.NBitsKey;
  Ocopy->NAloc = O.NAloc;
  Ocopy->QComplex = O.QComplex;
  if (O.NAloc)
  {
    if (Ocopy->QComplex)
    {
      if ((Ocopy->V = (Complex16*)malloc(Ocopy->NAloc * sizeof(Complex16))) == NULL)
      {
        printf("malloc failed for V in OperatorFixedLengthCopy\n");
        return 1;
      }
      memcpy(Ocopy->V, O.V, Ocopy->NAloc * sizeof(Complex16));
    }
    else
    {
      if ((Ocopy->VRe = (double*)malloc(Ocopy->NAloc * sizeof(double))) == NULL)
      {
        printf("malloc failed for VRe in OperatorFixedLengthCopy\n");
        return 1;
      }
      memcpy(Ocopy->VRe, O.VRe, Ocopy->NAloc * sizeof(double));
    }
    if (length)
    {
      if ((Ocopy->o = (uint16_t*)malloc(length * Ocopy->NAloc * sizeof(uint16_t))) == NULL)
      {
        printf("malloc failed for o in OperatorFixedLengthCopy\n");
        return 1;
      }
      memcpy(Ocopy->o, O.o, length * Ocopy->NAloc * sizeof(uint16_t));
      if ((Ocopy->Next = (uint32_t*)malloc(Ocopy->NAloc * sizeof(uint32_t))) == NULL)
      {
        printf("malloc failed for Next in OperatorFixedLengthCopy\n");
        return 1;
      }
      memcpy(Ocopy->Next, O.Next, Ocopy->NAloc * sizeof(uint32_t));
      memcpy(Ocopy->KeyList, O.KeyList, ((uint32_t) 1 << Ocopy->NBitsKey) * sizeof(uint32_t));
    }
  }
  return 0;
}

/*   Ocopy <- O
 */
uint32_t OperatorCopy(OperatorType * const Ocopy, OperatorType const O)
{
  uint32_t length;
  if (InitOperatorCopySettings(Ocopy, O))
  {
    printf("InitOperatorCopySettings failed in OperatorCopy\n");
    return 1;
  }
  for (length = 0; length <= Ocopy->MaxLength; length++)
  {
    if (OperatorFixedLengthCopy(&(Ocopy->l[length]), O.l[length], length))
    {
      printf("OperatorFixedLengthCopy for length=%i failed in OperatorCopy\n", length);
      return 1;
    }
  }
  SetQComplexOperator(Ocopy);
  return 0;
}


/*  epsilon might be set to DBL_MIN in order to prevent 0's not to be recognized
*   After diagonalization epsilon might be set to DBL_EPSILON*(...) with (...) the
*   maximum value in the operator in order to recognize all numerical unimportant elements.
*/
uint32_t OperatorChop(OperatorType * const O, const double epsilon)
{
  return OperatorRemoveSmallValues(O, epsilon);
}

uint32_t OperatorRemoveSmallValues(OperatorType * const O, const double epsilon)
{
  const double epssqr = Sqr(epsilon);
  int i, j, Nold, MakeComplexReal;
  /* wipe out all KeyLists */
  for (j = 1; j <= O->MaxLength; j++)
  {
    for (i = 0; i < ((uint32_t) 1 << (O->l[j].NBitsKey)); i++)
    {
      O->l[j].KeyList[i] = 0;
    }
  }
  if (O->l[0].QComplex)
  {                             /* Complex length 0 */
    if (O->l[0].N != 0)
    {
      if (fabs(O->l[0].V[0].Im) < epsilon)
      {
        if (fabs(O->l[0].V[0].Re) < epsilon)
        {                       /* both im an re part are small */
          free(O->l[0].V);
          O->l[0].QComplex = 0;
          O->l[0].N = 0;
          O->l[0].NAloc = 0;
        }
        else
        {                       /* Im part can be removed, real part not */
          if ((O->l[0].VRe = (double*)malloc(sizeof(double))) == NULL)
          {
            printf("could not allocate memory for O->l[0].VRe in OperatorRemoveSmallValues\n");
            return 1;
          }
          O->l[0].VRe[0] = O->l[0].V[0].Re;
          free(O->l[0].V);
          O->l[0].QComplex = 0;
          O->l[0].N = 1;
          O->l[0].NAloc = 1;
        }
        SetQComplexOperator(O);
      }
    }
  }
  else
  {                             /* real length 0 */
    if (O->l[0].N != 0)
    {
      if (fabs(O->l[0].VRe[0]) < epsilon)
      {
        free(O->l[0].VRe);
        O->l[0].N = 0;
        O->l[0].NAloc = 0;
        SetQComplexOperator(O);
      }
    }
  }
  for (j = 1; j <= O->MaxLength; j++)
  {
    Nold = O->l[j].N;
    if (Nold != 0)
    {
      O->l[j].N = 0;
      if (O->l[j].QComplex)
      {                         /* Complex */
        MakeComplexReal = 1;
        for (i = 0; i < Nold; i++)
        {
          if (fabs(O->l[j].V[i].Im) > epsilon) MakeComplexReal = 0;
          O->l[j].Next[i] = 0;
          if ((Sqr(O->l[j].V[i].Re) + Sqr(O->l[j].V[i].Im)) > epssqr)
          {
            if (ComplexOperatorAddLader(O, &(O->l[j].o[j * i]), j, O->l[j].V[i].Re, O->l[j].V[i].Im))
            {
              printf
                ("ComplexOperatorAddLader failed in OperatorRemoveSmallValues\n This should not happen\n No clue how to recover\n");
              AllocFailed();
            }
          }
        }
        if (MakeComplexReal)
        {
          if ((O->l[j].VRe = (double*)malloc(O->l[j].N * sizeof(double))) == NULL)
          {
            printf
              ("allocating memory for real part of operator failed. Continue with complex operator with Im entries that are zero.\n");
          }
          else
          {
            for (i = 0; i < O->l[j].N; i++)
            {
              O->l[j].VRe[i] = O->l[j].V[i].Re;
            }
            O->l[j].QComplex = 0;
            SetQComplexOperator(O);
          }
        }
        if ((O->l[j].V = (Complex16*)realloc(O->l[j].V, O->l[j].N * sizeof(Complex16))) == NULL)
          if (O->l[j].N != 0)
          {
            printf
              ("realloc failed for O->l[j=%i].V in OperatorRemoveSmallValues\n This should not happen\n No clue how to recover\n",
               j);
            AllocFailed();
          }
      }
      else
      {                         /* real */
        for (i = 0; i < Nold; i++)
        {
          O->l[j].Next[i] = 0;
          if (fabs(O->l[j].VRe[i]) > epsilon)
            if (RealOperatorAddLader(O, &(O->l[j].o[j * i]), j, O->l[j].VRe[i]))
            {
              printf("RealOperatorAddLader failed in OperatorRemoveSmallValues\n This should not happen\n No clue how to recover\n");
              AllocFailed();
            }
        }
        if ((O->l[j].VRe = (double*)realloc(O->l[j].VRe, O->l[j].N * sizeof(double))) == NULL)
          if (O->l[j].N != 0)
          {
            printf("realloc failed for O->l[j=%i].VRe in OperatorRemoveSmallValues\n This should not happen\n No clue how to recover\n",j);
            AllocFailed();
          }
      }
      if ((O->l[j].Next = (uint32_t*)realloc(O->l[j].Next, O->l[j].N * sizeof(uint32_t))) == NULL)
        if (O->l[j].N != 0)
        {
          printf("realloc failed for O->l[j=%i].Next in OperatorRemoveSmallValues\n This should not happen\n No clue how to recover\n",j);
          AllocFailed();
        }
      if ((O->l[j].o = (uint16_t*)realloc(O->l[j].o, O->l[j].N * j * sizeof(uint16_t))) == NULL)
        if (O->l[j].N != 0)
        {
          printf("realloc failed for O->l[j=%i].o in OperatorRemoveSmallValues\n This should not happen\n No clue how to recover\n",j);
          AllocFailed();
        }
      O->l[j].NAloc = O->l[j].N;
      if (O->l[j].N == 0) SetQComplexOperator(O);
    }
  }
  return 0;
}

uint32_t OperatorSetMaxLength(OperatorType * const O)
{
  while((O->MaxLength!=0)&&(O->l[O->MaxLength].N==0))
  {
    free(O->l[O->MaxLength].KeyList);
    if (O->l[O->MaxLength].NAloc)
    {
      free(O->l[O->MaxLength].o);
      free(O->l[O->MaxLength].Next);
      (O->l[O->MaxLength].QComplex) ? free(O->l[O->MaxLength].V) : free(O->l[O->MaxLength].VRe);
    }
    O->MaxLength--;
  }
  return 0;
}

/* Functions that determine if two particles are ordered or not.
 *  returns 0 if A is left of B
 *  returns 6 if A and B are fermions and A==B
 *  returns 1 if A should be right of B and A is a boson and B is a fermion
 *  returns 1 if A should be right of B and A is a fermion and B is a boson
 *  returns 2 if A should be right of B and both are fermions and A annihilates and B creates
 *  returns 2 if A should be right of B and both are fermions and A creates and B annihilates
 *  returns 3 if A should be right of B and both are fermions and A annihilates and B annihilates but A!=B
 *  returns 3 if A should be right of B and both are fermions and A creates and B creates but A!=B
 *  returns 4 if A should be right of B and both are bosons and A annihilates and B creates
 *  returns 5 if A should be right of B and both are bosons and A creates and B annihilates
 *  returns 1 if A should be right of B and both are bosons and A annihilates and B annihilates
 *  returns 1 if A should be right of B and both are bosons and A creates and B creates
 */
uint32_t QDetNotOrderedFaFcBaBc(const uint16_t A, const uint16_t B, const uint32_t NFermion)
{
  if ((A & NotCrAnMask) < NFermion)
  {                             /* A is a fermion */
    if ((B & NotCrAnMask) < NFermion)
    {                           /* A and B are fermions */
      if (A == B)
        return 6;
      if ((A ^ CrAnMask) > (B ^ CrAnMask))
        return 0;
      if ((A ^ B) & CrAnMask)
        return 2;
      return 3;
    }
    /* A is a fermion B is a boson */
    return 0;
  }
  /* A is a boson */
  if ((B & NotCrAnMask) < NFermion)
  {                             /* A is a boson and B is a fermion */
    return 1;
  }
  /* A and B are bosons */
  if ((A ^ CrAnMask) >= (B ^ CrAnMask))
    return 0;
  if ((A ^ B) & CrAnMask)
    return 5;                   /* A creates B annihilates */
  return 1;
}

uint32_t QDetNotOrderedFcFaBcBa(const uint16_t A, const uint16_t B, const uint32_t NFermion)
{
  if ((A & NotCrAnMask) < NFermion)
  {                             /* A is a fermion */
    if ((B & NotCrAnMask) < NFermion)
    {                           /* A and B are fermions */
      if (A == B)
        return 6;
      if (A > B)
        return 0;
      if ((A ^ B) & CrAnMask)
        return 2;
      return 3;
    }
    /* A is a fermion B is a boson */
    return 0;
  }
  /* A is a boson */
  if ((B & NotCrAnMask) < NFermion)
  {                             /* A is a boson and B is a fermion */
    return 1;
  }
  /* A and B are bosons */
  if (A >= B)
    return 0;
  if ((A ^ B) & CrAnMask)
    return 4;                   /* A annihilates B creates */
  return 1;
}

uint32_t QDetNotOrderedFcFaBaBc(const uint16_t A, const uint16_t B, const uint32_t NFermion)
{
  if ((A & NotCrAnMask) < NFermion)
  {                             /* A is a fermion */
    if ((B & NotCrAnMask) < NFermion)
    {                           /* A and B are fermions */
      if (A == B)
        return 6;
      if (A > B)
        return 0;
      if ((A ^ B) & CrAnMask)
        return 2;
      return 3;
    }
    /* A is a fermion B is a boson */
    return 0;
  }
  /* A is a boson */
  if ((B & NotCrAnMask) < NFermion)
  {                             /* A is a boson and B is a fermion */
    return 1;
  }
  /* A and B are bosons */
  if ((A ^ CrAnMask) >= (B ^ CrAnMask))
    return 0;
  if ((A ^ B) & CrAnMask)
    return 5;                   /* A creates B annihilates */
  return 1;
}

uint32_t QDetNotOrderedFaFcBcBa(const uint16_t A, const uint16_t B, const uint32_t NFermion)
{
  if ((A & NotCrAnMask) < NFermion)
  {                             /* A is a fermion */
    if ((B & NotCrAnMask) < NFermion)
    {                           /* A and B are fermions */
      if (A == B)
        return 6;
      if ((A ^ CrAnMask) > (B ^ CrAnMask))
        return 0;
      if ((A ^ B) & CrAnMask)
        return 2;
      return 3;
    }
    /* A is a fermion B is a boson */
    return 0;
  }
  /* A is a boson */
  if ((B & NotCrAnMask) < NFermion)
  {                             /* A is a boson and B is a fermion */
    return 1;
  }
  /* A and B are bosons */
  if (A >= B)
    return 0;
  if ((A ^ B) & CrAnMask)
    return 4;                   /* A annihilates B creates */
  return 1;
}

/*  Oo must be initiallized
 *  Oo <- Oo + NormalOrdered(O)
 *  It will be an ordered and optimized form of O
 */
uint32_t NormalOrderOperator(const OperatorType O, OperatorType * const Oo, uint32_t(*QDetNotOrdered) (uint16_t, uint16_t, uint32_t))
{
  uint32_t il, iN;
  if(O.l[0].N!=0)
  {
    if(O.l[0].QComplex==0)
    {
      OperatorAddLader(Oo, NULL, 0, O.l[0].VRe[0], 0);
    }
    else
    {
      OperatorAddLader(Oo, NULL, 0, O.l[0].V[0].Re, O.l[0].V[0].Im);
    }
  }
  for (il = 1; il <= O.MaxLength; il++)
  {
    if (O.l[il].QComplex)
    {                           /* O complex at level il */
      if (!(Oo->l[il].QComplex))
      {
        if (MakeRealOperatorFixedLengthComplex(Oo, il))
        {
          printf("MakeRealOperatorFixedLengthComplex at length=%i failed in NormalOrderOperator\n", il);
          fflush(stdout);
          return 1;
        }
      }
      for (iN = 0; iN < O.l[il].N; iN++)
      {
        if (ComplexOperatorAddLaderNormalOrder
            (Oo, &(O.l[il].o[il * iN]), il, O.l[il].V[iN].Re, O.l[il].V[iN].Im, QDetNotOrdered))
        {
          printf("ComplexOperatorAddLaderNormalOrder failed in NormalOrderOperator\n");
          fflush(stdout);
          return 1;
        }
      }
    }
    else
    {                           /* O real at level il */
      if (Oo->l[il].QComplex)
      {
        for (iN = 0; iN < O.l[il].N; iN++)
        {
          if (ComplexOperatorAddLaderNormalOrder(Oo, &(O.l[il].o[il * iN]), il, O.l[il].V[iN].Re, 0.0, QDetNotOrdered))
          {
            printf("ComplexOperatorAddLaderNormalOrder failed in NormalOrderOperator\n");
            fflush(stdout);
            return 1;
          }
        }
      }
      else
      {
        for (iN = 0; iN < O.l[il].N; iN++)
        {
          if (RealOperatorAddLaderNormalOrder(Oo, &(O.l[il].o[il * iN]), il, O.l[il].VRe[iN], QDetNotOrdered))
          {
            printf("RealOperatorAddLaderNormalOrder failed in NormalOrderOperator\n");
            fflush(stdout);
            return 1;
          }
        }
      }
    }
  }
  SetQComplexOperator(Oo);
  return 0;
}

/*
 *
 */
double RealOperatorTraceAverage(const OperatorType O, uint16_t *Orbitals, int NOrbitals)
{
  double trace=0;
  uint32_t i,j;
  if(O.MaxLength<2) return trace;
  for(i=0;i<O.l[2].N;i++)
  {
    if((O.l[2].o[2*i]&NotCrAnMask)==O.l[2].o[2*i+1])
    {
      for(j=0;j<NOrbitals;j++)
      {
        if(O.l[2].o[2*i+1]==Orbitals[j]) trace -= O.l[2].VRe[i];
      }
    }
    else if((O.l[2].o[2*i+1]&NotCrAnMask)==O.l[2].o[2*i])
    {
      for(j=0;j<NOrbitals;j++)
      {
        if(O.l[2].o[2*i]==Orbitals[j]) trace += O.l[2].VRe[i];
      }
    }
  }
  return trace/NOrbitals;
}
Complex16 ComplexOperatorTraceAverage(const OperatorType O, uint16_t *Orbitals, int NOrbitals)
{
  Complex16 trace={0,0};
  uint32_t i,j;
  if(O.MaxLength<2) return trace;
  for(i=0;i<O.l[2].N;i++)
  {
    if((O.l[2].o[2*i]&NotCrAnMask)==O.l[2].o[2*i+1])
    {
      for(j=0;j<NOrbitals;j++)
      {
        if(O.l[2].o[2*i+1]==Orbitals[j]) trace -= O.l[2].V[i];
      }
    }
    else if((O.l[2].o[2*i+1]&NotCrAnMask)==O.l[2].o[2*i])
    {
      for(j=0;j<NOrbitals;j++)
      {
        if(O.l[2].o[2*i]==Orbitals[j]) trace += O.l[2].V[i];
      }
    }
  }
  return trace/NOrbitals;
}
Complex16 OperatorTraceAverage(const OperatorType O, uint16_t *Orbitals, int NOrbitals)
{
  Complex16 traceAverage={0,0};
  if(O.MaxLength<2) return traceAverage;
  if(O.l[2].QComplex==0)
  {
    traceAverage.Re = RealOperatorTraceAverage(O,Orbitals,NOrbitals);
  }
  else
  {
    traceAverage = ComplexOperatorTraceAverage(O,Orbitals,NOrbitals);
  }
  return traceAverage;
}

uint32_t RealOperatorAddDiagonalOffset(OperatorType *O, uint16_t *Orbitals, int NOrbitals, double offset)
{
  uint32_t i;
  uint16_t oadd[2];
  for(i=0;i<NOrbitals;i++)
  {
    oadd[1]=CrAnMask + Orbitals[i];
    oadd[0]=Orbitals[i];
    RealOperatorAddLader(O, oadd, 2, offset);
  }
  return 0;
}

uint32_t ComplexOperatorAddDiagonalOffset(OperatorType *O, uint16_t *Orbitals, int NOrbitals, Complex16 offset)
{
  uint32_t i;
  uint16_t oadd[2];
  for(i=0;i<NOrbitals;i++)
  {
    oadd[1]=CrAnMask + Orbitals[i];
    oadd[0]=Orbitals[i];
    ComplexOperatorAddLader(O, oadd, 2, offset.Re, offset.Im);
  }
  return 0;
}




/*  Depending on Aset and Bset
 *    O <- O + B^'Bset' A^'Aset'
 *    whereby set can be
 *    'N' normal
 *    'C' conjugated (a^dag remains a^dag, prefactors become complex conjugated)
 *    'T' transposed (a^dag becomes a, prefactors remain)
 *    'D' conjugate transposed (a^dag becomes a, prefactors become complex conjugated)
 *
 *  O must be initiallized, i.e O.Maxlength must be A.Maxlength+B.Maxlength
 *  it is advised to set O.l[il].NBitsKey to sensible values
 */
#define HASH_QCOMPLEX(A, B, QC1, QC2, QC3) (((A)<<16) | ((B)<<8) | ((QC1!=0)<<4) | ((QC2!=0)<<2) | ((QC3!=0)<<0)   )


uint32_t OperatorMultiply2(OperatorType * const O, OperatorType const B, OperatorType const A, char const Bset, char const Aset)
{
  uint32_t ilA, ilB, iNA, iNB, k;
  uint16_t det[A.MaxLength + B.MaxLength];

  if (O->MaxLength < A.MaxLength + B.MaxLength)
  {
    printf("Length of output operator does not match in OperatorMultiply2\n");
    fflush(stdout);
    return 1;
  }

  for (ilA = 0; ilA <= A.MaxLength; ilA++)
    for (ilB = 0; ilB <= B.MaxLength; ilB++)
    {
      if (!(O->l[ilA + ilB].QComplex))
        if (A.l[ilA].QComplex || B.l[ilB].QComplex)
          if (MakeRealOperatorFixedLengthComplex(O, ilA + ilB))
          {
            printf("MakeRealOperatorFixedLengthComplex at length=%i failed in OperatorMultiply2\n", ilA + ilB);
            fflush(stdout);
            return 1;
          }

      switch (HASH_QCOMPLEX(Aset, Bset, O->l[ilA + ilB].QComplex, A.l[ilA].QComplex, B.l[ilB].QComplex))
      {
      case HASH_QCOMPLEX('N', 'N', 0, 0, 0):  // All Real
        for (iNA = 0; iNA < A.l[ilA].N; iNA++)
        {
          for (k = 0; k < ilA; k++)
            det[k] = A.l[ilA].o[iNA * ilA + k];
          for (iNB = 0; iNB < B.l[ilB].N; iNB++)
          {
            for (k = 0; k < ilB; k++)
              det[k + ilA] = B.l[ilB].o[iNB * ilB + k];

            if (RealOperatorAddLader(O, det, ilA + ilB, A.l[ilA].VRe[iNA] * B.l[ilB].VRe[iNB]))
            {
              printf("RealOperatorAddLader failed in OperatorMultiply2\n");
              fflush(stdout);
              return 1;
            }
          }
        }
        break;

      case HASH_QCOMPLEX('N', 'N', 0, 1, 1): // A and B complex
      case HASH_QCOMPLEX('N', 'N', 1, 1, 1):
        for (iNA = 0; iNA < A.l[ilA].N; iNA++)
        {
          for (k = 0; k < ilA; k++)
            det[k] = A.l[ilA].o[iNA * ilA + k];
          for (iNB = 0; iNB < B.l[ilB].N; iNB++)
          {
            for (k = 0; k < ilB; k++)
              det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
            if (ComplexOperatorAddLader
                (O, det, ilA + ilB,
                 A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Re - A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Im,
                 A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Im + A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Re))
            {
              printf("ComplexOperatorAddLader failed in OperatorMultiply\n");
              fflush(stdout);
              return 1;
            }
          }
        }
        break;

      case HASH_QCOMPLEX('N', 'N', 0, 1, 0): // A complex
      case HASH_QCOMPLEX('N', 'N', 1, 1, 0):
        for (iNA = 0; iNA < A.l[ilA].N; iNA++)
        {
          for (k = 0; k < ilA; k++)
            det[k] = A.l[ilA].o[iNA * ilA + k];
          for (iNB = 0; iNB < B.l[ilB].N; iNB++)
          {
            for (k = 0; k < ilB; k++)
              det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
            if (ComplexOperatorAddLader
                (O, det, ilA + ilB, A.l[ilA].V[iNA].Re * B.l[ilB].VRe[iNB], A.l[ilA].V[iNA].Im * B.l[ilB].VRe[iNB]))
            {
              printf("ComplexOperatorAddLader failed in OperatorMultiply2\n");
              fflush(stdout);
              return 1;
            }
          }
        }
        break;


      case HASH_QCOMPLEX('N', 'N', 0, 0, 1): // B complex
      case HASH_QCOMPLEX('N', 'N', 1, 0, 1):
        for (iNA = 0; iNA < A.l[ilA].N; iNA++)
        {
          for (k = 0; k < ilA; k++)
            det[k] = A.l[ilA].o[iNA * ilA + k];
          for (iNB = 0; iNB < B.l[ilB].N; iNB++)
          {
            for (k = 0; k < ilB; k++)
              det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
            if (ComplexOperatorAddLader
                (O, det, ilA + ilB, A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Re, A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Im))
            {
              printf("ComplexOperatorAddLader failed in OperatorMultiply2\n");
              fflush(stdout);
              return 1;
            }
          }
        }
        break;

      case HASH_QCOMPLEX('N', 'N', 1, 0, 0): // A and B real
        for (iNA = 0; iNA < A.l[ilA].N; iNA++)
        {
          for (k = 0; k < ilA; k++)
            det[k] = A.l[ilA].o[iNA * ilA + k];
          for (iNB = 0; iNB < B.l[ilB].N; iNB++)
          {
            for (k = 0; k < ilB; k++)
              det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
            if (ComplexOperatorAddLader(O, det, ilA + ilB, A.l[ilA].VRe[iNA] * B.l[ilB].VRe[iNB], 0))
            {
              printf("ComplexOperatorAddLader failed in OperatorMultiply2\n");
              fflush(stdout);
              return 1;
            }
          }
        }
        break;

      case HASH_QCOMPLEX('N', 'D', 0, 0, 0):
        for (iNA = 0; iNA < A.l[ilA].N; iNA++)
        {
          for (k = 0; k < ilA; k++)
            det[k] = A.l[ilA].o[iNA * ilA + k];
          for (iNB = 0; iNB < B.l[ilB].N; iNB++)
          {
            for (k = 0; k < ilB; k++)
              det[ilB - 1 - k + ilA] = B.l[ilB].o[iNB * ilB + k] ^ CrAnMask;
            if (RealOperatorAddLader(O, det, ilA + ilB, A.l[ilA].VRe[iNA] * B.l[ilB].VRe[iNB]))
            {
              printf("RealOperatorAddLader failed in OperatorMultiply2\n");
              fflush(stdout);
              return 1;
            }
          }
        }

        break;
      case HASH_QCOMPLEX('N', 'D', 1, 1, 1):
      case HASH_QCOMPLEX('N', 'D', 0, 1, 1):
        for (iNA = 0; iNA < A.l[ilA].N; iNA++)
        {
          for (k = 0; k < ilA; k++)
            det[k] = A.l[ilA].o[iNA * ilA + k];
          for (iNB = 0; iNB < B.l[ilB].N; iNB++)
          {
            for (k = 0; k < ilB; k++)
              det[ilB - 1 - k + ilA] = B.l[ilB].o[iNB * ilB + k] ^ CrAnMask;
            if (ComplexOperatorAddLader
                (O, det, ilA + ilB,
                 A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Re + A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Im,
                 -A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Im + A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Re))
            {
              printf("ComplexOperatorAddLader failed in OperatorMultiply\n");
              fflush(stdout);
              return 1;
            }
          }
        }
        break;

      case HASH_QCOMPLEX('N', 'D', 0, 1, 0):
      case HASH_QCOMPLEX('N', 'D', 1, 1, 0):
        for (iNA = 0; iNA < A.l[ilA].N; iNA++)
        {
          for (k = 0; k < ilA; k++)
            det[k] = A.l[ilA].o[iNA * ilA + k];
          for (iNB = 0; iNB < B.l[ilB].N; iNB++)
          {
            for (k = 0; k < ilB; k++)
              det[ilB - 1 - k + ilA] = B.l[ilB].o[iNB * ilB + k] ^ CrAnMask;
            if (ComplexOperatorAddLader
                (O, det, ilA + ilB, A.l[ilA].V[iNA].Re * B.l[ilB].VRe[iNB], A.l[ilA].V[iNA].Im * B.l[ilB].VRe[iNB]))
            {
              printf("ComplexOperatorAddLader failed in OperatorMultiply2\n");
              fflush(stdout);
              return 1;
            }
          }
        }
        break;

      case HASH_QCOMPLEX('N', 'D', 0, 0, 1):
      case HASH_QCOMPLEX('N', 'D', 1, 0, 1):
        for (iNA = 0; iNA < A.l[ilA].N; iNA++)
        {
          for (k = 0; k < ilA; k++)
            det[k] = A.l[ilA].o[iNA * ilA + k];
          for (iNB = 0; iNB < B.l[ilB].N; iNB++)
          {
            for (k = 0; k < ilB; k++)
              det[ilB - 1 - k + ilA] = B.l[ilB].o[iNB * ilB + k] ^ CrAnMask;
            if (ComplexOperatorAddLader
                (O, det, ilA + ilB, A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Re, -A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Im))
            {
              printf("ComplexOperatorAddLader failed in OperatorMultiply2\n");
              fflush(stdout);
              return 1;
            }
          }
        }
        break;

      case HASH_QCOMPLEX('N', 'D', 1, 0, 0):
        for (iNA = 0; iNA < A.l[ilA].N; iNA++)
        {
          for (k = 0; k < ilA; k++)
            det[k] = A.l[ilA].o[iNA * ilA + k];
          for (iNB = 0; iNB < B.l[ilB].N; iNB++)
          {
            for (k = 0; k < ilB; k++)
              det[ilB - 1 - k + ilA] = B.l[ilB].o[iNB * ilB + k] ^ CrAnMask;
            if (ComplexOperatorAddLader(O, det, ilA + ilB, A.l[ilA].VRe[iNA] * B.l[ilB].VRe[iNB], 0))
            {
              printf("ComplexOperatorAddLader failed in OperatorMultiply2\n");
              fflush(stdout);
              return 1;
            }
          }
        }
        break;

      case HASH_QCOMPLEX('D', 'N', 0, 0, 0):
        for (iNA = 0; iNA < A.l[ilA].N; iNA++)
        {
          for (k = 0; k < ilA; k++)
            det[ilA - 1 - k] = A.l[ilA].o[iNA * ilA + k] ^ CrAnMask;
          for (iNB = 0; iNB < B.l[ilB].N; iNB++)
          {
            for (k = 0; k < ilB; k++)
              det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
            if (RealOperatorAddLader(O, det, ilA + ilB, A.l[ilA].VRe[iNA] * B.l[ilB].VRe[iNB]))
            {
              printf("RealOperatorAddLader failed in OperatorMultiply2\n");
              fflush(stdout);
              return 1;
            }
          }
        }
        break;

      case HASH_QCOMPLEX('D', 'N', 1, 1, 1):
      case HASH_QCOMPLEX('D', 'N', 0, 1, 1):
        for (iNA = 0; iNA < A.l[ilA].N; iNA++)
        {
          for (k = 0; k < ilA; k++)
            det[ilA - 1 - k] = A.l[ilA].o[iNA * ilA + k] ^ CrAnMask;
          for (iNB = 0; iNB < B.l[ilB].N; iNB++)
          {
            for (k = 0; k < ilB; k++)
              det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
            if (ComplexOperatorAddLader
                (O, det, ilA + ilB,
                 A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Re + A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Im,
                 A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Im - A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Re))
            {
              printf("ComplexOperatorAddLader failed in OperatorMultiply2\n");
              fflush(stdout);
              return 1;
            }
          }
        }
        break;

      case HASH_QCOMPLEX('D', 'N', 1, 1, 0):
      case HASH_QCOMPLEX('D', 'N', 0, 1, 0):
        for (iNA = 0; iNA < A.l[ilA].N; iNA++)
        {
          for (k = 0; k < ilA; k++)
            det[ilA - 1 - k] = A.l[ilA].o[iNA * ilA + k] ^ CrAnMask;
          for (iNB = 0; iNB < B.l[ilB].N; iNB++)
          {
            for (k = 0; k < ilB; k++)
              det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
            if (ComplexOperatorAddLader
                (O, det, ilA + ilB, A.l[ilA].V[iNA].Re * B.l[ilB].VRe[iNB], -A.l[ilA].V[iNA].Im * B.l[ilB].VRe[iNB]))
            {
              printf("ComplexOperatorAddLader failed in OperatorMultiply2\n");
              fflush(stdout);
              return 1;
            }
          }
        }
        break;

      case HASH_QCOMPLEX('D', 'N', 0, 0, 1):
      case HASH_QCOMPLEX('D', 'N', 1, 0, 1):
        for (iNA = 0; iNA < A.l[ilA].N; iNA++)
        {
          for (k = 0; k < ilA; k++)
            det[ilA - 1 - k] = A.l[ilA].o[iNA * ilA + k] ^ CrAnMask;
          for (iNB = 0; iNB < B.l[ilB].N; iNB++)
          {
            for (k = 0; k < ilB; k++)
              det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
            if (ComplexOperatorAddLader
                (O, det, ilA + ilB, A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Re, A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Im))
            {
              printf("ComplexOperatorAddLader failed in OperatorMultiply2\n");
              fflush(stdout);
              return 1;
            }
          }
        }
        break;

      case HASH_QCOMPLEX('D', 'N', 1, 0, 0):
        for (iNA = 0; iNA < A.l[ilA].N; iNA++)
        {
          for (k = 0; k < ilA; k++)
            det[ilA - 1 - k] = A.l[ilA].o[iNA * ilA + k] ^ CrAnMask;
          for (iNB = 0; iNB < B.l[ilB].N; iNB++)
          {
            for (k = 0; k < ilB; k++)
              det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
            if (ComplexOperatorAddLader(O, det, ilA + ilB, A.l[ilA].VRe[iNA] * B.l[ilB].VRe[iNB], 0))
            {
              printf("ComplexOperatorAddLader failed in OperatorMultiply2\n");
              fflush(stdout);
              return 1;
            }
          }
        }
        break;

      case HASH_QCOMPLEX('D', 'D', 0, 0, 0):
        for (iNA = 0; iNA < A.l[ilA].N; iNA++)
        {
          for (k = 0; k < ilA; k++)
            det[ilA - 1 - k] = A.l[ilA].o[iNA * ilA + k] ^ CrAnMask;
          for (iNB = 0; iNB < B.l[ilB].N; iNB++)
          {
            for (k = 0; k < ilB; k++)
              det[ilB - 1 - k + ilA] = B.l[ilB].o[iNB * ilB + k] ^ CrAnMask;
            if (RealOperatorAddLader(O, det, ilA + ilB, A.l[ilA].VRe[iNA] * B.l[ilB].VRe[iNB]))
            {
              printf("RealOperatorAddLader failed in OperatorMultiply2\n");
              fflush(stdout);
              return 1;
            }
          }
        }
        break;

      case HASH_QCOMPLEX('D', 'D', 1, 1, 1):
      case HASH_QCOMPLEX('D', 'D', 0, 1, 1):
        for (iNA = 0; iNA < A.l[ilA].N; iNA++)
        {
          for (k = 0; k < ilA; k++)
            det[ilA - 1 - k] = A.l[ilA].o[iNA * ilA + k] ^ CrAnMask;
          for (iNB = 0; iNB < B.l[ilB].N; iNB++)
          {
            for (k = 0; k < ilB; k++)
              det[ilB - 1 - k + ilA] = B.l[ilB].o[iNB * ilB + k] ^ CrAnMask;
            if (ComplexOperatorAddLader
                (O, det, ilA + ilB,
                 A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Re - A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Im,
                 -A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Im - A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Re))
            {
              printf("ComplexOperatorAddLader failed in OperatorMultiply2\n");
              fflush(stdout);
              return 1;
            }
          }
        }
        break;

      case HASH_QCOMPLEX('D', 'D', 1, 1, 0):
      case HASH_QCOMPLEX('D', 'D', 0, 1, 0):
        for (iNA = 0; iNA < A.l[ilA].N; iNA++)
        {
          for (k = 0; k < ilA; k++)
            det[ilA - 1 - k] = A.l[ilA].o[iNA * ilA + k] ^ CrAnMask;
          for (iNB = 0; iNB < B.l[ilB].N; iNB++)
          {
            for (k = 0; k < ilB; k++)
              det[ilB - 1 - k + ilA] = B.l[ilB].o[iNB * ilB + k] ^ CrAnMask;
            if (ComplexOperatorAddLader
                (O, det, ilA + ilB, A.l[ilA].V[iNA].Re * B.l[ilB].VRe[iNB], -A.l[ilA].V[iNA].Im * B.l[ilB].VRe[iNB]))
            {
              printf("ComplexOperatorAddLader failed in OperatorMultiply2\n");
              fflush(stdout);
              return 1;
            }
          }
        }
        break;

      case HASH_QCOMPLEX('D', 'D', 0, 0, 1):
      case HASH_QCOMPLEX('D', 'D', 1, 0, 1):
        for (iNA = 0; iNA < A.l[ilA].N; iNA++)
        {
          for (k = 0; k < ilA; k++)
            det[ilA - 1 - k] = A.l[ilA].o[iNA * ilA + k] ^ CrAnMask;
          for (iNB = 0; iNB < B.l[ilB].N; iNB++)
          {
            for (k = 0; k < ilB; k++)
              det[ilB - 1 - k + ilA] = B.l[ilB].o[iNB * ilB + k] ^ CrAnMask;
            if (ComplexOperatorAddLader
                (O, det, ilA + ilB, A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Re, -A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Im))
            {
              printf("ComplexOperatorAddLader failed in OperatorMultiply2\n");
              fflush(stdout);
              return 1;
            }
          }
        }
        break;

      case HASH_QCOMPLEX('D', 'D', 1, 0, 0):
        for (iNA = 0; iNA < A.l[ilA].N; iNA++)
        {
          for (k = 0; k < ilA; k++)
            det[ilA - 1 - k] = A.l[ilA].o[iNA * ilA + k] ^ CrAnMask;
          for (iNB = 0; iNB < B.l[ilB].N; iNB++)
          {
            for (k = 0; k < ilB; k++)
              det[ilB - 1 - k + ilA] = B.l[ilB].o[iNB * ilB + k] ^ CrAnMask;
            if (ComplexOperatorAddLader(O, det, ilA + ilB, A.l[ilA].VRe[iNA] * B.l[ilB].VRe[iNB], 0))
            {
              printf("ComplexOperatorAddLader failed in OperatorMultiply2\n");
              fflush(stdout);
              return 1;
            }
          }
        }
        break;

      default:
        printf("in OperatorMultiply2 unrecognized input for Aset and or Bset\n");
        return 1;
      }                         // switch  
    }                          // for loop
  SetQComplexOperator(O);
  return 0;
}
uint32_t OperatorMultiply(OperatorType * const O, OperatorType const B, OperatorType const A)
{
  return OperatorMultiply2(O, B, A, 'N', 'N');
}
/* O <- O*B */
uint32_t OperatorMultiplyInline(OperatorType * const O, OperatorType const B)
{
  OperatorType temp;
  InitOperatorForMult(&temp, 2, O,&B);
  OperatorMultiply(&temp, *O, B);
  strcpy(temp.Name, O->Name);
  FreeOperator(O);
  strcpy(O->Name,temp.Name);
  O->QComplex=temp.QComplex;
  O->NFermion=temp.NFermion;
  O->NBoson=temp.NBoson;
  O->NWord=temp.NWord;
  O->MaxLength=temp.MaxLength;
  O->l=temp.l;
  return 0;
}


/*  Depending on Aset and Bset
 *    O <- O + A^'Aset' B^'Bset'
 *    whereby set can be
 *    'N' normal
 *    'C' conjugated (a^dag remains a^dag, prefactors become complex conjugated)
 *    'T' transposed (a^dag becomes a, prefactors remain)
 *    'D' conjugate transposed (a^dag becomes a, prefactors become complex conjugated)
 *
 *  O must be initiallized, i.e O.Maxlength must be A.Maxlength+B.Maxlength
 *  it is advised to set O.l[il].NBitsKey to sensible values
 */
uint32_t OperatorMultiply2NormalOrder(OperatorType * const O, OperatorType const B, OperatorType const A, char const Bset, char const Aset, uint32_t(*QDetNotOrdered) (uint16_t, uint16_t, uint32_t))
{
  uint32_t ilA, ilB, iNA, iNB, k;
  uint16_t det[A.MaxLength + B.MaxLength];

  if (O->MaxLength != A.MaxLength + B.MaxLength)
  {
    printf("Length of output operator does not match in OperatorMultiply2NormalOrder\n");
    fflush(stdout);
    return 1;
  }

  for (ilA = 0; ilA <= A.MaxLength; ilA++)
    for (ilB = 0; ilB <= B.MaxLength; ilB++)
    {
      if (!(O->l[ilA + ilB].QComplex))
        if (A.l[ilA].QComplex || B.l[ilB].QComplex)
          if (MakeRealOperatorFixedLengthComplex(O, ilA + ilB))
          {
            printf("MakeRealOperatorFixedLengthComplex at length=%i failed in OperatorMultiply2NormalOrder\n", ilA + ilB);
            fflush(stdout);
            return 1;
          }
          SetQComplexOperator(O);
      switch (HASH_QCOMPLEX(Aset, Bset, O->l[ilA + ilB].QComplex, A.l[ilA].QComplex, B.l[ilB].QComplex))
      {
      case HASH_QCOMPLEX('N', 'N', 0, 0, 0):  // All Real
        for (iNA = 0; iNA < A.l[ilA].N; iNA++)
        {
          for (k = 0; k < ilA; k++) det[k] = A.l[ilA].o[iNA * ilA + k];
          for (iNB = 0; iNB < B.l[ilB].N; iNB++)
          {
            for (k = 0; k < ilB; k++) det[k + ilA] = B.l[ilB].o[iNB * ilB + k];

            if (RealOperatorAddLaderNormalOrder(O, det, ilA + ilB, A.l[ilA].VRe[iNA] * B.l[ilB].VRe[iNB], QDetNotOrdered))
            {
              printf("RealOperatorAddLaderNormalOrder failed in OperatorMultiply2NormalOrder\n");
              fflush(stdout);
              return 1;
            }
          }
        }
        break;

      case HASH_QCOMPLEX('N', 'N', 0, 1, 1): // A and B complex
      case HASH_QCOMPLEX('N', 'N', 1, 1, 1):
        for (iNA = 0; iNA < A.l[ilA].N; iNA++)
        {
          for (k = 0; k < ilA; k++)
            det[k] = A.l[ilA].o[iNA * ilA + k];
          for (iNB = 0; iNB < B.l[ilB].N; iNB++)
          {
            for (k = 0; k < ilB; k++) det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
            if (ComplexOperatorAddLaderNormalOrder
                (O, det, ilA + ilB,
                 A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Re - A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Im,
                 A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Im + A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Re, QDetNotOrdered))
            {
              printf("ComplexOperatorAddLaderNormalOrder failed in OperatorMultiplyNormalOrder\n");
              fflush(stdout);
              return 1;
            }
          }
        }
        break;

      case HASH_QCOMPLEX('N', 'N', 0, 1, 0): // A complex
      case HASH_QCOMPLEX('N', 'N', 1, 1, 0):
        for (iNA = 0; iNA < A.l[ilA].N; iNA++)
        {
          for (k = 0; k < ilA; k++) det[k] = A.l[ilA].o[iNA * ilA + k];
          for (iNB = 0; iNB < B.l[ilB].N; iNB++)
          {
            for (k = 0; k < ilB; k++) det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
            if (ComplexOperatorAddLaderNormalOrder
                (O, det, ilA + ilB, A.l[ilA].V[iNA].Re * B.l[ilB].VRe[iNB], A.l[ilA].V[iNA].Im * B.l[ilB].VRe[iNB], QDetNotOrdered))
            {
              printf("ComplexOperatorAddLaderNormalOrder failed in OperatorMultiply2NormalOrder\n");
              fflush(stdout);
              return 1;
            }
          }
        }
        break;


      case HASH_QCOMPLEX('N', 'N', 0, 0, 1): // B complex
      case HASH_QCOMPLEX('N', 'N', 1, 0, 1):
        for (iNA = 0; iNA < A.l[ilA].N; iNA++)
        {
          for (k = 0; k < ilA; k++) det[k] = A.l[ilA].o[iNA * ilA + k];
          for (iNB = 0; iNB < B.l[ilB].N; iNB++)
          {
            for (k = 0; k < ilB; k++) det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
            if (ComplexOperatorAddLaderNormalOrder
                (O, det, ilA + ilB, A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Re, A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Im, QDetNotOrdered))
            {
              printf("ComplexOperatorAddLaderNormalOrder failed in OperatorMultiply2NormalOrder\n");
              fflush(stdout);
              return 1;
            }
          }
        }
        break;

      case HASH_QCOMPLEX('N', 'N', 1, 0, 0): // A and B real
        for (iNA = 0; iNA < A.l[ilA].N; iNA++)
        {
          for (k = 0; k < ilA; k++) det[k] = A.l[ilA].o[iNA * ilA + k];
          for (iNB = 0; iNB < B.l[ilB].N; iNB++)
          {
            for (k = 0; k < ilB; k++) det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
            if (ComplexOperatorAddLaderNormalOrder(O, det, ilA + ilB, A.l[ilA].VRe[iNA] * B.l[ilB].VRe[iNB], 0.0, QDetNotOrdered))
            {
              printf("ComplexOperatorAddLaderNormalOrder failed in OperatorMultiply2NormalOrder\n");
              fflush(stdout);
              return 1;
            }
          }
        }
        break;

      case HASH_QCOMPLEX('N', 'D', 0, 0, 0):
        for (iNA = 0; iNA < A.l[ilA].N; iNA++)
        {
          for (k = 0; k < ilA; k++)  det[k] = A.l[ilA].o[iNA * ilA + k];
          for (iNB = 0; iNB < B.l[ilB].N; iNB++)
          {
            for (k = 0; k < ilB; k++) det[ilB - 1 - k + ilA] = B.l[ilB].o[iNB * ilB + k] ^ CrAnMask;
            if (RealOperatorAddLaderNormalOrder(O, det, ilA + ilB, A.l[ilA].VRe[iNA] * B.l[ilB].VRe[iNB], QDetNotOrdered))
            {
              printf("RealOperatorAddLaderNormalOrder failed in OperatorMultiply2NormalOrder\n");
              fflush(stdout);
              return 1;
            }
          }
        }

        break;
      case HASH_QCOMPLEX('N', 'D', 1, 1, 1):
      case HASH_QCOMPLEX('N', 'D', 0, 1, 1):
        for (iNA = 0; iNA < A.l[ilA].N; iNA++)
        {
          for (k = 0; k < ilA; k++)
            det[k] = A.l[ilA].o[iNA * ilA + k];
          for (iNB = 0; iNB < B.l[ilB].N; iNB++)
          {
            for (k = 0; k < ilB; k++)
              det[ilB - 1 - k + ilA] = B.l[ilB].o[iNB * ilB + k] ^ CrAnMask;
            if (ComplexOperatorAddLaderNormalOrder
                (O, det, ilA + ilB,
                 A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Re + A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Im,
                 -A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Im + A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Re, QDetNotOrdered))
            {
              printf("ComplexOperatorAddLaderNormalOrder failed in OperatorMultiplyNormalOrder\n");
              fflush(stdout);
              return 1;
            }
          }
        }
        break;

      case HASH_QCOMPLEX('N', 'D', 0, 1, 0):
      case HASH_QCOMPLEX('N', 'D', 1, 1, 0):
        for (iNA = 0; iNA < A.l[ilA].N; iNA++)
        {
          for (k = 0; k < ilA; k++)
            det[k] = A.l[ilA].o[iNA * ilA + k];
          for (iNB = 0; iNB < B.l[ilB].N; iNB++)
          {
            for (k = 0; k < ilB; k++)
              det[ilB - 1 - k + ilA] = B.l[ilB].o[iNB * ilB + k] ^ CrAnMask;
            if (ComplexOperatorAddLaderNormalOrder
                (O, det, ilA + ilB, A.l[ilA].V[iNA].Re * B.l[ilB].VRe[iNB], A.l[ilA].V[iNA].Im * B.l[ilB].VRe[iNB], QDetNotOrdered))
            {
              printf("ComplexOperatorAddLaderNormalOrder failed in OperatorMultiply2NormalOrder\n");
              fflush(stdout);
              return 1;
            }
          }
        }
        break;

      case HASH_QCOMPLEX('N', 'D', 0, 0, 1):
      case HASH_QCOMPLEX('N', 'D', 1, 0, 1):
        for (iNA = 0; iNA < A.l[ilA].N; iNA++)
        {
          for (k = 0; k < ilA; k++)
            det[k] = A.l[ilA].o[iNA * ilA + k];
          for (iNB = 0; iNB < B.l[ilB].N; iNB++)
          {
            for (k = 0; k < ilB; k++)
              det[ilB - 1 - k + ilA] = B.l[ilB].o[iNB * ilB + k] ^ CrAnMask;
            if (ComplexOperatorAddLaderNormalOrder
                (O, det, ilA + ilB, A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Re, -A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Im, QDetNotOrdered))
            {
              printf("ComplexOperatorAddLaderNormalOrder failed in OperatorMultiply2NormalOrder\n");
              fflush(stdout);
              return 1;
            }
          }
        }
        break;

      case HASH_QCOMPLEX('N', 'D', 1, 0, 0):
        for (iNA = 0; iNA < A.l[ilA].N; iNA++)
        {
          for (k = 0; k < ilA; k++)
            det[k] = A.l[ilA].o[iNA * ilA + k];
          for (iNB = 0; iNB < B.l[ilB].N; iNB++)
          {
            for (k = 0; k < ilB; k++)
              det[ilB - 1 - k + ilA] = B.l[ilB].o[iNB * ilB + k] ^ CrAnMask;
            if (ComplexOperatorAddLaderNormalOrder(O, det, ilA + ilB, A.l[ilA].VRe[iNA] * B.l[ilB].VRe[iNB], 0.0, QDetNotOrdered))
            {
              printf("ComplexOperatorAddLaderNormalOrder failed in OperatorMultiply2NormalOrder\n");
              fflush(stdout);
              return 1;
            }
          }
        }
        break;

      case HASH_QCOMPLEX('D', 'N', 0, 0, 0):
        for (iNA = 0; iNA < A.l[ilA].N; iNA++)
        {
          for (k = 0; k < ilA; k++)
            det[ilA - 1 - k] = A.l[ilA].o[iNA * ilA + k] ^ CrAnMask;
          for (iNB = 0; iNB < B.l[ilB].N; iNB++)
          {
            for (k = 0; k < ilB; k++)
              det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
            if (RealOperatorAddLaderNormalOrder(O, det, ilA + ilB, A.l[ilA].VRe[iNA] * B.l[ilB].VRe[iNB], QDetNotOrdered))
            {
              printf("RealOperatorAddLaderNormalOrder failed in OperatorMultiply2NormalOrder\n");
              fflush(stdout);
              return 1;
            }
          }
        }
        break;

      case HASH_QCOMPLEX('D', 'N', 1, 1, 1):
      case HASH_QCOMPLEX('D', 'N', 0, 1, 1):
        for (iNA = 0; iNA < A.l[ilA].N; iNA++)
        {
          for (k = 0; k < ilA; k++)
            det[ilA - 1 - k] = A.l[ilA].o[iNA * ilA + k] ^ CrAnMask;
          for (iNB = 0; iNB < B.l[ilB].N; iNB++)
          {
            for (k = 0; k < ilB; k++)
              det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
            if (ComplexOperatorAddLaderNormalOrder
                (O, det, ilA + ilB,
                 A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Re + A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Im,
                 A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Im - A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Re, QDetNotOrdered))
            {
              printf("ComplexOperatorAddLaderNormalOrder failed in OperatorMultiply2NormalOrder\n");
              fflush(stdout);
              return 1;
            }
          }
        }
        break;

      case HASH_QCOMPLEX('D', 'N', 1, 1, 0):
      case HASH_QCOMPLEX('D', 'N', 0, 1, 0):
        for (iNA = 0; iNA < A.l[ilA].N; iNA++)
        {
          for (k = 0; k < ilA; k++)
            det[ilA - 1 - k] = A.l[ilA].o[iNA * ilA + k] ^ CrAnMask;
          for (iNB = 0; iNB < B.l[ilB].N; iNB++)
          {
            for (k = 0; k < ilB; k++)
              det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
            if (ComplexOperatorAddLaderNormalOrder
                (O, det, ilA + ilB, A.l[ilA].V[iNA].Re * B.l[ilB].VRe[iNB], -A.l[ilA].V[iNA].Im * B.l[ilB].VRe[iNB], QDetNotOrdered))
            {
              printf("ComplexOperatorAddLaderNormalOrder failed in OperatorMultiply2NormalOrder\n");
              fflush(stdout);
              return 1;
            }
          }
        }
        break;

      case HASH_QCOMPLEX('D', 'N', 0, 0, 1):
      case HASH_QCOMPLEX('D', 'N', 1, 0, 1):
        for (iNA = 0; iNA < A.l[ilA].N; iNA++)
        {
          for (k = 0; k < ilA; k++)
            det[ilA - 1 - k] = A.l[ilA].o[iNA * ilA + k] ^ CrAnMask;
          for (iNB = 0; iNB < B.l[ilB].N; iNB++)
          {
            for (k = 0; k < ilB; k++)
              det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
            if (ComplexOperatorAddLaderNormalOrder
                (O, det, ilA + ilB, A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Re, A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Im, QDetNotOrdered))
            {
              printf("ComplexOperatorAddLader failedNormalOrder in OperatorMultiply2NormalOrder\n");
              fflush(stdout);
              return 1;
            }
          }
        }
        break;

      case HASH_QCOMPLEX('D', 'N', 1, 0, 0):
        for (iNA = 0; iNA < A.l[ilA].N; iNA++)
        {
          for (k = 0; k < ilA; k++)
            det[ilA - 1 - k] = A.l[ilA].o[iNA * ilA + k] ^ CrAnMask;
          for (iNB = 0; iNB < B.l[ilB].N; iNB++)
          {
            for (k = 0; k < ilB; k++)
              det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
            if (ComplexOperatorAddLaderNormalOrder(O, det, ilA + ilB, A.l[ilA].VRe[iNA] * B.l[ilB].VRe[iNB], 0.0, QDetNotOrdered))
            {
              printf("ComplexOperatorAddLaderNormalOrder failed in OperatorMultiply2NormalOrder\n");
              fflush(stdout);
              return 1;
            }
          }
        }
        break;

      case HASH_QCOMPLEX('D', 'D', 0, 0, 0):
        for (iNA = 0; iNA < A.l[ilA].N; iNA++)
        {
          for (k = 0; k < ilA; k++)
            det[ilA - 1 - k] = A.l[ilA].o[iNA * ilA + k] ^ CrAnMask;
          for (iNB = 0; iNB < B.l[ilB].N; iNB++)
          {
            for (k = 0; k < ilB; k++)
              det[ilB - 1 - k + ilA] = B.l[ilB].o[iNB * ilB + k] ^ CrAnMask;
            if (RealOperatorAddLaderNormalOrder(O, det, ilA + ilB, A.l[ilA].VRe[iNA] * B.l[ilB].VRe[iNB], QDetNotOrdered))
            {
              printf("RealOperatorAddLaderNormalOrder failed in OperatorMultiply2NormalOrder\n");
              fflush(stdout);
              return 1;
            }
          }
        }
        break;

      case HASH_QCOMPLEX('D', 'D', 1, 1, 1):
      case HASH_QCOMPLEX('D', 'D', 0, 1, 1):
        for (iNA = 0; iNA < A.l[ilA].N; iNA++)
        {
          for (k = 0; k < ilA; k++)
            det[ilA - 1 - k] = A.l[ilA].o[iNA * ilA + k] ^ CrAnMask;
          for (iNB = 0; iNB < B.l[ilB].N; iNB++)
          {
            for (k = 0; k < ilB; k++)
              det[ilB - 1 - k + ilA] = B.l[ilB].o[iNB * ilB + k] ^ CrAnMask;
            if (ComplexOperatorAddLaderNormalOrder
                (O, det, ilA + ilB,
                 A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Re - A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Im,
                 -A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Im - A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Re, QDetNotOrdered))
            {
              printf("ComplexOperatorAddLaderNormalOrder failed in OperatorMultiply2NormalOrder\n");
              fflush(stdout);
              return 1;
            }
          }
        }
        break;

      case HASH_QCOMPLEX('D', 'D', 1, 1, 0):
      case HASH_QCOMPLEX('D', 'D', 0, 1, 0):
        for (iNA = 0; iNA < A.l[ilA].N; iNA++)
        {
          for (k = 0; k < ilA; k++)
            det[ilA - 1 - k] = A.l[ilA].o[iNA * ilA + k] ^ CrAnMask;
          for (iNB = 0; iNB < B.l[ilB].N; iNB++)
          {
            for (k = 0; k < ilB; k++)
              det[ilB - 1 - k + ilA] = B.l[ilB].o[iNB * ilB + k] ^ CrAnMask;
            if (ComplexOperatorAddLaderNormalOrder
                (O, det, ilA + ilB, A.l[ilA].V[iNA].Re * B.l[ilB].VRe[iNB], -A.l[ilA].V[iNA].Im * B.l[ilB].VRe[iNB], QDetNotOrdered))
            {
              printf("ComplexOperatorAddLaderNormalOrder failed in OperatorMultiply2NormalOrder\n");
              fflush(stdout);
              return 1;
            }
          }
        }
        break;

      case HASH_QCOMPLEX('D', 'D', 0, 0, 1):
      case HASH_QCOMPLEX('D', 'D', 1, 0, 1):
        for (iNA = 0; iNA < A.l[ilA].N; iNA++)
        {
          for (k = 0; k < ilA; k++)
            det[ilA - 1 - k] = A.l[ilA].o[iNA * ilA + k] ^ CrAnMask;
          for (iNB = 0; iNB < B.l[ilB].N; iNB++)
          {
            for (k = 0; k < ilB; k++)
              det[ilB - 1 - k + ilA] = B.l[ilB].o[iNB * ilB + k] ^ CrAnMask;
            if (ComplexOperatorAddLaderNormalOrder
                (O, det, ilA + ilB, A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Re, -A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Im, QDetNotOrdered))
            {
              printf("ComplexOperatorAddLaderNormalOrder failed in OperatorMultiply2NormalOrder\n");
              fflush(stdout);
              return 1;
            }
          }
        }
        break;

      case HASH_QCOMPLEX('D', 'D', 1, 0, 0):
        for (iNA = 0; iNA < A.l[ilA].N; iNA++)
        {
          for (k = 0; k < ilA; k++)
            det[ilA - 1 - k] = A.l[ilA].o[iNA * ilA + k] ^ CrAnMask;
          for (iNB = 0; iNB < B.l[ilB].N; iNB++)
          {
            for (k = 0; k < ilB; k++)
              det[ilB - 1 - k + ilA] = B.l[ilB].o[iNB * ilB + k] ^ CrAnMask;
            if (ComplexOperatorAddLaderNormalOrder(O, det, ilA + ilB, A.l[ilA].VRe[iNA] * B.l[ilB].VRe[iNB], 0.0, QDetNotOrdered))
            {
              printf("ComplexOperatorAddLaderNormalOrder failed in OperatorMultiply2NormalOrder\n");
              fflush(stdout);
              return 1;
            }
          }
        }
        break;

      default:
        printf("in OperatorMultiply2 unrecognized input for Aset and or Bset\n");
        return 1;
      }                         // switch  
    }                          // for loop
  SetQComplexOperator(O);
  return 0;
}
uint32_t OperatorMultiplyNormalOrder(OperatorType * const O, OperatorType const B, OperatorType const A, uint32_t(*QDetNotOrdered) (uint16_t, uint16_t, uint32_t))
{
  return OperatorMultiply2NormalOrder(O, B, A, 'N', 'N', QDetNotOrdered);
}
/* O <- O*B */
uint32_t OperatorMultiplyNormalOrderInline(OperatorType * const O, OperatorType const B, uint32_t(*QDetNotOrdered) (uint16_t, uint16_t, uint32_t))
{
  OperatorType temp;
  InitOperatorForMult(&temp, 2, O,&B);
  OperatorMultiplyNormalOrder(&temp, *O, B, QDetNotOrdered);
  strcpy(temp.Name, O->Name);
  FreeOperator(O);
  strcpy(O->Name,temp.Name);
  O->QComplex=temp.QComplex;
  O->NFermion=temp.NFermion;
  O->NBoson=temp.NBoson;
  O->NWord=temp.NWord;
  O->MaxLength=temp.MaxLength;
  O->l=temp.l;
  return 0;
}

#undef HASH_QCOMPLEX


/*  Depending on Aset, Bset and Cset
 *    O <- O + C^'Cset' B^'Bset' A^'Aset
 *    whereby set can be
 *    'N' normal
 *    'C' conjugated (a^dag remains a^dag, prefactors become complex conjugated)
 *    'T' transposed (a^dag becomes a, prefactors remain)
 *    'D' conjugate transposed (a^dag becomes a, prefactors become complex conjugated)
 */
uint32_t OperatorMultiply3(OperatorType * const O, OperatorType const C, OperatorType const B, OperatorType const A, char const Cset, char const Bset, char const Aset)
{
  uint32_t ilA, ilB, ilC, iNA, iNB, iNC, k;
  uint16_t det[A.MaxLength + B.MaxLength + C.MaxLength];

  if (O->MaxLength != A.MaxLength + B.MaxLength + C.MaxLength)
  {
    printf("Length of output operator does not match in OperatorMultiply3\n");
    fflush(stdout);
    return 1;
  }

  if ((Aset == 'N') && (Bset == 'N') && (Cset == 'N'))
  {
    for (ilA = 0; ilA <= A.MaxLength; ilA++)
    {
      for (ilB = 0; ilB <= B.MaxLength; ilB++)
      {
        for (ilC = 0; ilC <= C.MaxLength; ilC++)
        {
          if ((O->l[ilA + ilB].QComplex) || (A.l[ilA].QComplex) || (B.l[ilB].QComplex) || (C.l[ilC].QComplex))
          {                     /* at least one complex */
            if (!(O->l[ilA + ilB + ilC].QComplex))
              if (MakeRealOperatorFixedLengthComplex(O, ilA + ilB + ilC))
              {
                printf("MakeRealOperatorFixedLengthComplex at length=%i failed in OperatorMultiply3\n",
                       ilA + ilB + ilC);
                fflush(stdout);
                return 1;
              }
            if (A.l[ilA].QComplex)
            {                   /* A complex */
              if (B.l[ilB].QComplex)
              {                 /* A and B complex */
                if (C.l[ilC].QComplex)
                {               /* A and B and C complex */
                  for (iNA = 0; iNA < A.l[ilA].N; iNA++)
                  {
                    for (k = 0; k < ilA; k++)
                    {
                      det[k] = A.l[ilA].o[iNA * ilA + k];
                    }
                    for (iNB = 0; iNB < B.l[ilB].N; iNB++)
                    {
                      for (k = 0; k < ilB; k++)
                      {
                        det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                      }
                      for (iNC = 0; iNC < C.l[ilC].N; iNC++)
                      {
                        for (k = 0; k < ilC; k++)
                        {
                          det[k + ilA + ilB] = C.l[ilC].o[iNC * ilC + k];
                        }
                        if (ComplexOperatorAddLader
                            (O, det, ilA + ilB + ilC,
                             A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Re * C.l[ilC].V[iNC].Re -
                             A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Im * C.l[ilC].V[iNC].Im -
                             A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Re * C.l[ilC].V[iNC].Im -
                             A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Im * C.l[ilC].V[iNC].Re,
                             A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Re * C.l[ilC].V[iNC].Re +
                             A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Im * C.l[ilC].V[iNC].Re +
                             A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Re * C.l[ilC].V[iNC].Im -
                             A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Im * C.l[ilC].V[iNC].Im))
                        {
                          printf("ComplexOperatorAddLader failed in OperatorMultiply3\n");
                          fflush(stdout);
                          return 1;
                        }
                      }
                    }
                  }
                }
                else
                {               /* A and B complex, C real */
                  for (iNA = 0; iNA < A.l[ilA].N; iNA++)
                  {
                    for (k = 0; k < ilA; k++)
                    {
                      det[k] = A.l[ilA].o[iNA * ilA + k];
                    }
                    for (iNB = 0; iNB < B.l[ilB].N; iNB++)
                    {
                      for (k = 0; k < ilB; k++)
                      {
                        det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                      }
                      for (iNC = 0; iNC < C.l[ilC].N; iNC++)
                      {
                        for (k = 0; k < ilC; k++)
                        {
                          det[k + ilA + ilB] = C.l[ilC].o[iNC * ilC + k];
                        }
                        if (ComplexOperatorAddLader
                            (O, det, ilA + ilB + ilC,
                             A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Re * C.l[ilC].VRe[iNC] -
                             A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Im * C.l[ilC].VRe[iNC],
                             A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Re * C.l[ilC].VRe[iNC] +
                             A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Im * C.l[ilC].VRe[iNC]))
                        {
                          printf("ComplexOperatorAddLader failed in OperatorMultiply3\n");
                          fflush(stdout);
                          return 1;
                        }
                      }
                    }
                  }
                }
              }
              else
              {                 /* A complex B real */
                if (C.l[ilC].QComplex)
                {               /* A complex B real C complex */
                  for (iNA = 0; iNA < A.l[ilA].N; iNA++)
                  {
                    for (k = 0; k < ilA; k++)
                    {
                      det[k] = A.l[ilA].o[iNA * ilA + k];
                    }
                    for (iNB = 0; iNB < B.l[ilB].N; iNB++)
                    {
                      for (k = 0; k < ilB; k++)
                      {
                        det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                      }
                      for (iNC = 0; iNC < C.l[ilC].N; iNC++)
                      {
                        for (k = 0; k < ilC; k++)
                        {
                          det[k + ilA + ilB] = C.l[ilC].o[iNC * ilC + k];
                        }
                        if (ComplexOperatorAddLader
                            (O, det, ilA + ilB + ilC,
                             A.l[ilA].V[iNA].Re * B.l[ilB].VRe[iNB] * C.l[ilC].V[iNC].Re -
                             A.l[ilA].V[iNA].Im * B.l[ilB].VRe[iNB] * C.l[ilC].V[iNC].Im,
                             A.l[ilA].V[iNA].Im * B.l[ilB].VRe[iNB] * C.l[ilC].V[iNC].Re +
                             A.l[ilA].V[iNA].Re * B.l[ilB].VRe[iNB] * C.l[ilC].V[iNC].Im))
                        {
                          printf("ComplexOperatorAddLader failed in OperatorMultiply3\n");
                          fflush(stdout);
                          return 1;
                        }
                      }
                    }
                  }
                }
                else
                {               /* A complex B real C real */
                  for (iNA = 0; iNA < A.l[ilA].N; iNA++)
                  {
                    for (k = 0; k < ilA; k++)
                    {
                      det[k] = A.l[ilA].o[iNA * ilA + k];
                    }
                    for (iNB = 0; iNB < B.l[ilB].N; iNB++)
                    {
                      for (k = 0; k < ilB; k++)
                      {
                        det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                      }
                      for (iNC = 0; iNC < C.l[ilC].N; iNC++)
                      {
                        for (k = 0; k < ilC; k++)
                        {
                          det[k + ilA + ilB] = C.l[ilC].o[iNC * ilC + k];
                        }
                        if (ComplexOperatorAddLader
                            (O, det, ilA + ilB + ilC, A.l[ilA].V[iNA].Re * B.l[ilB].VRe[iNB] * C.l[ilC].VRe[iNC],
                             A.l[ilA].V[iNA].Im * B.l[ilB].VRe[iNB] * C.l[ilC].VRe[iNC]))
                        {
                          printf("ComplexOperatorAddLader failed in OperatorMultiply3\n");
                          fflush(stdout);
                          return 1;
                        }
                      }
                    }
                  }
                }
              }
            }
            else
            {                   /* A real */
              if (B.l[ilB].QComplex)
              {                 /* A real and B complex */
                if (C.l[ilC].QComplex)
                {               /* A real B complex C complex */
                  for (iNA = 0; iNA < A.l[ilA].N; iNA++)
                  {
                    for (k = 0; k < ilA; k++)
                    {
                      det[k] = A.l[ilA].o[iNA * ilA + k];
                    }
                    for (iNB = 0; iNB < B.l[ilB].N; iNB++)
                    {
                      for (k = 0; k < ilB; k++)
                      {
                        det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                      }
                      for (iNC = 0; iNC < C.l[ilC].N; iNC++)
                      {
                        for (k = 0; k < ilC; k++)
                        {
                          det[k + ilA + ilB] = C.l[ilC].o[iNC * ilC + k];
                        }
                        if (ComplexOperatorAddLader
                            (O, det, ilA + ilB + ilC,
                             A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Re * C.l[ilC].V[iNC].Re -
                             A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Im * C.l[ilC].V[iNC].Im,
                             A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Im * C.l[ilC].V[iNC].Re +
                             A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Re * C.l[ilC].V[iNC].Im))
                        {
                          printf("ComplexOperatorAddLader failed in OperatorMultiply3\n");
                          fflush(stdout);
                          return 1;
                        }
                      }
                    }
                  }
                }
                else
                {               /* A real B complex C real */
                  for (iNA = 0; iNA < A.l[ilA].N; iNA++)
                  {
                    for (k = 0; k < ilA; k++)
                    {
                      det[k] = A.l[ilA].o[iNA * ilA + k];
                    }
                    for (iNB = 0; iNB < B.l[ilB].N; iNB++)
                    {
                      for (k = 0; k < ilB; k++)
                      {
                        det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                      }
                      for (iNC = 0; iNC < C.l[ilC].N; iNC++)
                      {
                        for (k = 0; k < ilC; k++)
                        {
                          det[k + ilA + ilB] = C.l[ilC].o[iNC * ilC + k];
                        }
                        if (ComplexOperatorAddLader
                            (O, det, ilA + ilB + ilC, A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Re * C.l[ilC].VRe[iNC],
                             A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Im * C.l[ilC].VRe[iNC]))
                        {
                          printf("ComplexOperatorAddLader failed in OperatorMultiply3\n");
                          fflush(stdout);
                          return 1;
                        }
                      }
                    }
                  }
                }
              }
              else
              {                 /* A and B real */
                if (C.l[ilC].QComplex)
                {               /* A and B real C complex */
                  for (iNA = 0; iNA < A.l[ilA].N; iNA++)
                  {
                    for (k = 0; k < ilA; k++)
                    {
                      det[k] = A.l[ilA].o[iNA * ilA + k];
                    }
                    for (iNB = 0; iNB < B.l[ilB].N; iNB++)
                    {
                      for (k = 0; k < ilB; k++)
                      {
                        det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                      }
                      for (iNC = 0; iNC < C.l[ilC].N; iNC++)
                      {
                        for (k = 0; k < ilC; k++)
                        {
                          det[k + ilA + ilB] = C.l[ilC].o[iNC * ilC + k];
                        }
                        if (ComplexOperatorAddLader
                            (O, det, ilA + ilB + ilC, A.l[ilA].VRe[iNA] * B.l[ilB].VRe[iNB] * C.l[ilC].V[iNC].Re,
                             A.l[ilA].VRe[iNA] * B.l[ilB].VRe[iNB] * C.l[ilC].V[iNC].Im))
                        {
                          printf("ComplexOperatorAddLader failed in OperatorMultiply3\n");
                          fflush(stdout);
                          return 1;
                        }
                      }
                    }
                  }
                }
                else
                {               /* A, B and C real (O complex) */
                  for (iNA = 0; iNA < A.l[ilA].N; iNA++)
                  {
                    for (k = 0; k < ilA; k++)
                    {
                      det[k] = A.l[ilA].o[iNA * ilA + k];
                    }
                    for (iNB = 0; iNB < B.l[ilB].N; iNB++)
                    {
                      for (k = 0; k < ilB; k++)
                      {
                        det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                      }
                      for (iNC = 0; iNC < C.l[ilC].N; iNC++)
                      {
                        for (k = 0; k < ilC; k++)
                        {
                          det[k + ilA + ilB] = C.l[ilC].o[iNC * ilC + k];
                        }
                        if (ComplexOperatorAddLader
                            (O, det, ilA + ilB + ilC, A.l[ilA].VRe[iNA] * B.l[ilB].VRe[iNB] * C.l[ilC].VRe[iNC], 0))
                        {
                          printf("ComplexOperatorAddLader failed in OperatorMultiply3\n");
                          fflush(stdout);
                          return 1;
                        }
                      }
                    }
                  }
                }
              }
            }
          }
          else
          {                     /* All real */
            for (iNA = 0; iNA < A.l[ilA].N; iNA++)
            {
              for (k = 0; k < ilA; k++)
              {
                det[k] = A.l[ilA].o[iNA * ilA + k];
              }
              for (iNB = 0; iNB < B.l[ilB].N; iNB++)
              {
                for (k = 0; k < ilB; k++)
                {
                  det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                }
                for (iNC = 0; iNC < C.l[ilC].N; iNC++)
                {
                  for (k = 0; k < ilC; k++)
                  {
                    det[k + ilA + ilB] = C.l[ilC].o[iNC * ilC + k];
                  }
                  if (RealOperatorAddLader
                      (O, det, ilA + ilB + ilC, A.l[ilA].VRe[iNA] * B.l[ilB].VRe[iNB] * C.l[ilC].VRe[iNC]))
                  {
                    printf("RealOperatorAddLader failed in OperatorMultiply3\n");
                    fflush(stdout);
                    return 1;
                  }
                }
              }
            }
          }
        }
      }
    }
    /* End of NNN */
  }
  else if ((Aset == 'N') && (Bset == 'D') && (Cset == 'N'))
  {
    printf("in OperatorMultiply3 Still to do\n");
    return 1;
  }
  else if ((Aset == 'D') && (Bset == 'N') && (Cset == 'N'))
  {
    for (ilA = 0; ilA <= A.MaxLength; ilA++)
    {
      for (ilB = 0; ilB <= B.MaxLength; ilB++)
      {
        for (ilC = 0; ilC <= C.MaxLength; ilC++)
        {
          if ((O->l[ilA + ilB].QComplex) || (A.l[ilA].QComplex) || (B.l[ilB].QComplex) || (C.l[ilC].QComplex))
          {                     /* at least one complex */
            if (!(O->l[ilA + ilB + ilC].QComplex))
              if (MakeRealOperatorFixedLengthComplex(O, ilA + ilB + ilC))
              {
                printf("MakeRealOperatorFixedLengthComplex at length=%i failed in OperatorMultiply3\n",
                       ilA + ilB + ilC);
                fflush(stdout);
                return 1;
              }
            if (A.l[ilA].QComplex)
            {                   /* A complex */
              if (B.l[ilB].QComplex)
              {                 /* A and B complex */
                if (C.l[ilC].QComplex)
                {               /* A and B and C complex */
                  for (iNA = 0; iNA < A.l[ilA].N; iNA++)
                  {
                    for (k = 0; k < ilA; k++)
                    {
                      det[ilA - 1 - k] = A.l[ilA].o[iNA * ilA + k] ^ CrAnMask;
                    }
                    for (iNB = 0; iNB < B.l[ilB].N; iNB++)
                    {
                      for (k = 0; k < ilB; k++)
                      {
                        det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                      }
                      for (iNC = 0; iNC < C.l[ilC].N; iNC++)
                      {
                        for (k = 0; k < ilC; k++)
                        {
                          det[k + ilA + ilB] = C.l[ilC].o[iNC * ilC + k];
                        }
                        if (ComplexOperatorAddLader
                            (O, det, ilA + ilB + ilC,
                             A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Re * C.l[ilC].V[iNC].Re -
                             A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Im * C.l[ilC].V[iNC].Im +
                             A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Re * C.l[ilC].V[iNC].Im +
                             A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Im * C.l[ilC].V[iNC].Re,
                             -A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Re * C.l[ilC].V[iNC].Re +
                             A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Im * C.l[ilC].V[iNC].Re +
                             A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Re * C.l[ilC].V[iNC].Im +
                             A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Im * C.l[ilC].V[iNC].Im))
                        {
                          printf("ComplexOperatorAddLader failed in OperatorMultiply3\n");
                          fflush(stdout);
                          return 1;
                        }
                      }
                    }
                  }
                }
                else
                {               /* A and B complex, C real */
                  for (iNA = 0; iNA < A.l[ilA].N; iNA++)
                  {
                    for (k = 0; k < ilA; k++)
                    {
                      det[ilA - 1 - k] = A.l[ilA].o[iNA * ilA + k] ^ CrAnMask;
                    }
                    for (iNB = 0; iNB < B.l[ilB].N; iNB++)
                    {
                      for (k = 0; k < ilB; k++)
                      {
                        det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                      }
                      for (iNC = 0; iNC < C.l[ilC].N; iNC++)
                      {
                        for (k = 0; k < ilC; k++)
                        {
                          det[k + ilA + ilB] = C.l[ilC].o[iNC * ilC + k];
                        }
                        if (ComplexOperatorAddLader
                            (O, det, ilA + ilB + ilC,
                             A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Re * C.l[ilC].VRe[iNC] +
                             A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Im * C.l[ilC].VRe[iNC],
                             -A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Re * C.l[ilC].VRe[iNC] +
                             A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Im * C.l[ilC].VRe[iNC]))
                        {
                          printf("ComplexOperatorAddLader failed in OperatorMultiply3\n");
                          fflush(stdout);
                          return 1;
                        }
                      }
                    }
                  }
                }
              }
              else
              {                 /* A complex B real */
                if (C.l[ilC].QComplex)
                {               /* A complex B real C complex */
                  for (iNA = 0; iNA < A.l[ilA].N; iNA++)
                  {
                    for (k = 0; k < ilA; k++)
                    {
                      det[ilA - 1 - k] = A.l[ilA].o[iNA * ilA + k] ^ CrAnMask;
                    }
                    for (iNB = 0; iNB < B.l[ilB].N; iNB++)
                    {
                      for (k = 0; k < ilB; k++)
                      {
                        det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                      }
                      for (iNC = 0; iNC < C.l[ilC].N; iNC++)
                      {
                        for (k = 0; k < ilC; k++)
                        {
                          det[k + ilA + ilB] = C.l[ilC].o[iNC * ilC + k];
                        }
                        if (ComplexOperatorAddLader
                            (O, det, ilA + ilB + ilC,
                             A.l[ilA].V[iNA].Re * B.l[ilB].VRe[iNB] * C.l[ilC].V[iNC].Re +
                             A.l[ilA].V[iNA].Im * B.l[ilB].VRe[iNB] * C.l[ilC].V[iNC].Im,
                             -A.l[ilA].V[iNA].Im * B.l[ilB].VRe[iNB] * C.l[ilC].V[iNC].Re +
                             A.l[ilA].V[iNA].Re * B.l[ilB].VRe[iNB] * C.l[ilC].V[iNC].Im))
                        {
                          printf("ComplexOperatorAddLader failed in OperatorMultiply3\n");
                          fflush(stdout);
                          return 1;
                        }
                      }
                    }
                  }
                }
                else
                {               /* A complex B real C real */
                  for (iNA = 0; iNA < A.l[ilA].N; iNA++)
                  {
                    for (k = 0; k < ilA; k++)
                    {
                      det[ilA - 1 - k] = A.l[ilA].o[iNA * ilA + k] ^ CrAnMask;
                    }
                    for (iNB = 0; iNB < B.l[ilB].N; iNB++)
                    {
                      for (k = 0; k < ilB; k++)
                      {
                        det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                      }
                      for (iNC = 0; iNC < C.l[ilC].N; iNC++)
                      {
                        for (k = 0; k < ilC; k++)
                        {
                          det[k + ilA + ilB] = C.l[ilC].o[iNC * ilC + k];
                        }
                        if (ComplexOperatorAddLader
                            (O, det, ilA + ilB + ilC, A.l[ilA].V[iNA].Re * B.l[ilB].VRe[iNB] * C.l[ilC].VRe[iNC],
                             -A.l[ilA].V[iNA].Im * B.l[ilB].VRe[iNB] * C.l[ilC].VRe[iNC]))
                        {
                          printf("ComplexOperatorAddLader failed in OperatorMultiply3\n");
                          fflush(stdout);
                          return 1;
                        }
                      }
                    }
                  }
                }
              }
            }
            else
            {                   /* A real */
              if (B.l[ilB].QComplex)
              {                 /* A real and B complex */
                if (C.l[ilC].QComplex)
                {               /* A real B complex C complex */
                  for (iNA = 0; iNA < A.l[ilA].N; iNA++)
                  {
                    for (k = 0; k < ilA; k++)
                    {
                      det[ilA - 1 - k] = A.l[ilA].o[iNA * ilA + k] ^ CrAnMask;
                    }
                    for (iNB = 0; iNB < B.l[ilB].N; iNB++)
                    {
                      for (k = 0; k < ilB; k++)
                      {
                        det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                      }
                      for (iNC = 0; iNC < C.l[ilC].N; iNC++)
                      {
                        for (k = 0; k < ilC; k++)
                        {
                          det[k + ilA + ilB] = C.l[ilC].o[iNC * ilC + k];
                        }
                        if (ComplexOperatorAddLader
                            (O, det, ilA + ilB + ilC,
                             A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Re * C.l[ilC].V[iNC].Re -
                             A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Im * C.l[ilC].V[iNC].Im,
                             A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Im * C.l[ilC].V[iNC].Re +
                             A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Re * C.l[ilC].V[iNC].Im))
                        {
                          printf("ComplexOperatorAddLader failed in OperatorMultiply3\n");
                          fflush(stdout);
                          return 1;
                        }
                      }
                    }
                  }
                }
                else
                {               /* A real B complex C real */
                  for (iNA = 0; iNA < A.l[ilA].N; iNA++)
                  {
                    for (k = 0; k < ilA; k++)
                    {
                      det[ilA - 1 - k] = A.l[ilA].o[iNA * ilA + k] ^ CrAnMask;
                    }
                    for (iNB = 0; iNB < B.l[ilB].N; iNB++)
                    {
                      for (k = 0; k < ilB; k++)
                      {
                        det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                      }
                      for (iNC = 0; iNC < C.l[ilC].N; iNC++)
                      {
                        for (k = 0; k < ilC; k++)
                        {
                          det[k + ilA + ilB] = C.l[ilC].o[iNC * ilC + k];
                        }
                        if (ComplexOperatorAddLader
                            (O, det, ilA + ilB + ilC, A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Re * C.l[ilC].VRe[iNC],
                             A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Im * C.l[ilC].VRe[iNC]))
                        {
                          printf("ComplexOperatorAddLader failed in OperatorMultiply3\n");
                          fflush(stdout);
                          return 1;
                        }
                      }
                    }
                  }
                }
              }
              else
              {                 /* A and B real */
                if (C.l[ilC].QComplex)
                {               /* A and B real C complex */
                  for (iNA = 0; iNA < A.l[ilA].N; iNA++)
                  {
                    for (k = 0; k < ilA; k++)
                    {
                      det[ilA - 1 - k] = A.l[ilA].o[iNA * ilA + k] ^ CrAnMask;
                    }
                    for (iNB = 0; iNB < B.l[ilB].N; iNB++)
                    {
                      for (k = 0; k < ilB; k++)
                      {
                        det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                      }
                      for (iNC = 0; iNC < C.l[ilC].N; iNC++)
                      {
                        for (k = 0; k < ilC; k++)
                        {
                          det[k + ilA + ilB] = C.l[ilC].o[iNC * ilC + k];
                        }
                        if (ComplexOperatorAddLader
                            (O, det, ilA + ilB + ilC, A.l[ilA].VRe[iNA] * B.l[ilB].VRe[iNB] * C.l[ilC].V[iNC].Re,
                             A.l[ilA].VRe[iNA] * B.l[ilB].VRe[iNB] * C.l[ilC].V[iNC].Im))
                        {
                          printf("ComplexOperatorAddLader failed in OperatorMultiply3\n");
                          fflush(stdout);
                          return 1;
                        }
                      }
                    }
                  }
                }
                else
                {               /* A, B and C real (O complex) */
                  for (iNA = 0; iNA < A.l[ilA].N; iNA++)
                  {
                    for (k = 0; k < ilA; k++)
                    {
                      det[ilA - 1 - k] = A.l[ilA].o[iNA * ilA + k] ^ CrAnMask;
                    }
                    for (iNB = 0; iNB < B.l[ilB].N; iNB++)
                    {
                      for (k = 0; k < ilB; k++)
                      {
                        det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                      }
                      for (iNC = 0; iNC < C.l[ilC].N; iNC++)
                      {
                        for (k = 0; k < ilC; k++)
                        {
                          det[k + ilA + ilB] = C.l[ilC].o[iNC * ilC + k];
                        }
                        if (ComplexOperatorAddLader
                            (O, det, ilA + ilB + ilC, A.l[ilA].VRe[iNA] * B.l[ilB].VRe[iNB] * C.l[ilC].VRe[iNC], 0))
                        {
                          printf("ComplexOperatorAddLader failed in OperatorMultiply3\n");
                          fflush(stdout);
                          return 1;
                        }
                      }
                    }
                  }
                }
              }
            }
          }
          else
          {                     /* All real */
            for (iNA = 0; iNA < A.l[ilA].N; iNA++)
            {
              for (k = 0; k < ilA; k++)
              {
                det[ilA - 1 - k] = A.l[ilA].o[iNA * ilA + k] ^ CrAnMask;
              }
              for (iNB = 0; iNB < B.l[ilB].N; iNB++)
              {
                for (k = 0; k < ilB; k++)
                {
                  det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                }
                for (iNC = 0; iNC < C.l[ilC].N; iNC++)
                {
                  for (k = 0; k < ilC; k++)
                  {
                    det[k + ilA + ilB] = C.l[ilC].o[iNC * ilC + k];
                  }
                  if (RealOperatorAddLader
                      (O, det, ilA + ilB + ilC, A.l[ilA].VRe[iNA] * B.l[ilB].VRe[iNB] * C.l[ilC].VRe[iNC]))
                  {
                    printf("RealOperatorAddLader failed in OperatorMultiply3\n");
                    fflush(stdout);
                    return 1;
                  }
                }
              }
            }
          }
        }
      }
    }
    /* End of DNN */
  }
  else if ((Aset == 'D') && (Bset == 'D') && (Cset == 'N'))
  {
    printf("in OperatorMultiply3 Still to do\n");
    return 1;
  }
  else if ((Aset == 'N') && (Bset == 'N') && (Cset == 'D'))
  {
    for (ilA = 0; ilA <= A.MaxLength; ilA++)
    {
      for (ilB = 0; ilB <= B.MaxLength; ilB++)
      {
        for (ilC = 0; ilC <= C.MaxLength; ilC++)
        {
          if ((O->l[ilA + ilB].QComplex) || (A.l[ilA].QComplex) || (B.l[ilB].QComplex) || (C.l[ilC].QComplex))
          {                     /* at least one complex */
            if (!(O->l[ilA + ilB + ilC].QComplex))
              if (MakeRealOperatorFixedLengthComplex(O, ilA + ilB + ilC))
              {
                printf("MakeRealOperatorFixedLengthComplex at length=%i failed in OperatorMultiply3\n",
                       ilA + ilB + ilC);
                fflush(stdout);
                return 1;
              }
            if (A.l[ilA].QComplex)
            {                   /* A complex */
              if (B.l[ilB].QComplex)
              {                 /* A and B complex */
                if (C.l[ilC].QComplex)
                {               /* A and B and C complex */
                  for (iNA = 0; iNA < A.l[ilA].N; iNA++)
                  {
                    for (k = 0; k < ilA; k++)
                    {
                      det[k] = A.l[ilA].o[iNA * ilA + k];
                    }
                    for (iNB = 0; iNB < B.l[ilB].N; iNB++)
                    {
                      for (k = 0; k < ilB; k++)
                      {
                        det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                      }
                      for (iNC = 0; iNC < C.l[ilC].N; iNC++)
                      {
                        for (k = 0; k < ilC; k++)
                        {
                          det[ilC - 1 - k + ilA + ilB] = C.l[ilC].o[iNC * ilC + k] ^ CrAnMask;
                        }
                        if (ComplexOperatorAddLader
                            (O, det, ilA + ilB + ilC,
                             A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Re * C.l[ilC].V[iNC].Re +
                             A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Im * C.l[ilC].V[iNC].Im +
                             A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Re * C.l[ilC].V[iNC].Im -
                             A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Im * C.l[ilC].V[iNC].Re,
                             A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Re * C.l[ilC].V[iNC].Re +
                             A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Im * C.l[ilC].V[iNC].Re -
                             A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Re * C.l[ilC].V[iNC].Im +
                             A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Im * C.l[ilC].V[iNC].Im))
                        {
                          printf("ComplexOperatorAddLader failed in OperatorMultiply3\n");
                          fflush(stdout);
                          return 1;
                        }
                      }
                    }
                  }
                }
                else
                {               /* A and B complex, C real */
                  for (iNA = 0; iNA < A.l[ilA].N; iNA++)
                  {
                    for (k = 0; k < ilA; k++)
                    {
                      det[k] = A.l[ilA].o[iNA * ilA + k];
                    }
                    for (iNB = 0; iNB < B.l[ilB].N; iNB++)
                    {
                      for (k = 0; k < ilB; k++)
                      {
                        det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                      }
                      for (iNC = 0; iNC < C.l[ilC].N; iNC++)
                      {
                        for (k = 0; k < ilC; k++)
                        {
                          det[ilC - 1 - k + ilA + ilB] = C.l[ilC].o[iNC * ilC + k] ^ CrAnMask;
                        }
                        if (ComplexOperatorAddLader
                            (O, det, ilA + ilB + ilC,
                             A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Re * C.l[ilC].VRe[iNC] -
                             A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Im * C.l[ilC].VRe[iNC],
                             A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Re * C.l[ilC].VRe[iNC] +
                             A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Im * C.l[ilC].VRe[iNC]))
                        {
                          printf("ComplexOperatorAddLader failed in OperatorMultiply3\n");
                          fflush(stdout);
                          return 1;
                        }
                      }
                    }
                  }
                }
              }
              else
              {                 /* A complex B real */
                if (C.l[ilC].QComplex)
                {               /* A complex B real C complex */
                  for (iNA = 0; iNA < A.l[ilA].N; iNA++)
                  {
                    for (k = 0; k < ilA; k++)
                    {
                      det[k] = A.l[ilA].o[iNA * ilA + k];
                    }
                    for (iNB = 0; iNB < B.l[ilB].N; iNB++)
                    {
                      for (k = 0; k < ilB; k++)
                      {
                        det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                      }
                      for (iNC = 0; iNC < C.l[ilC].N; iNC++)
                      {
                        for (k = 0; k < ilC; k++)
                        {
                          det[ilC - 1 - k + ilA + ilB] = C.l[ilC].o[iNC * ilC + k] ^ CrAnMask;
                        }
                        if (ComplexOperatorAddLader
                            (O, det, ilA + ilB + ilC,
                             A.l[ilA].V[iNA].Re * B.l[ilB].VRe[iNB] * C.l[ilC].V[iNC].Re +
                             A.l[ilA].V[iNA].Im * B.l[ilB].VRe[iNB] * C.l[ilC].V[iNC].Im,
                             A.l[ilA].V[iNA].Im * B.l[ilB].VRe[iNB] * C.l[ilC].V[iNC].Re -
                             A.l[ilA].V[iNA].Re * B.l[ilB].VRe[iNB] * C.l[ilC].V[iNC].Im))
                        {
                          printf("ComplexOperatorAddLader failed in OperatorMultiply3\n");
                          fflush(stdout);
                          return 1;
                        }
                      }
                    }
                  }
                }
                else
                {               /* A complex B real C real */
                  for (iNA = 0; iNA < A.l[ilA].N; iNA++)
                  {
                    for (k = 0; k < ilA; k++)
                    {
                      det[k] = A.l[ilA].o[iNA * ilA + k];
                    }
                    for (iNB = 0; iNB < B.l[ilB].N; iNB++)
                    {
                      for (k = 0; k < ilB; k++)
                      {
                        det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                      }
                      for (iNC = 0; iNC < C.l[ilC].N; iNC++)
                      {
                        for (k = 0; k < ilC; k++)
                        {
                          det[ilC - 1 - k + ilA + ilB] = C.l[ilC].o[iNC * ilC + k] ^ CrAnMask;
                        }
                        if (ComplexOperatorAddLader
                            (O, det, ilA + ilB + ilC, A.l[ilA].V[iNA].Re * B.l[ilB].VRe[iNB] * C.l[ilC].VRe[iNC],
                             A.l[ilA].V[iNA].Im * B.l[ilB].VRe[iNB] * C.l[ilC].VRe[iNC]))
                        {
                          printf("ComplexOperatorAddLader failed in OperatorMultiply3\n");
                          fflush(stdout);
                          return 1;
                        }
                      }
                    }
                  }
                }
              }
            }
            else
            {                   /* A real */
              if (B.l[ilB].QComplex)
              {                 /* A real and B complex */
                if (C.l[ilC].QComplex)
                {               /* A real B complex C complex */
                  for (iNA = 0; iNA < A.l[ilA].N; iNA++)
                  {
                    for (k = 0; k < ilA; k++)
                    {
                      det[k] = A.l[ilA].o[iNA * ilA + k];
                    }
                    for (iNB = 0; iNB < B.l[ilB].N; iNB++)
                    {
                      for (k = 0; k < ilB; k++)
                      {
                        det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                      }
                      for (iNC = 0; iNC < C.l[ilC].N; iNC++)
                      {
                        for (k = 0; k < ilC; k++)
                        {
                          det[ilC - 1 - k + ilA + ilB] = C.l[ilC].o[iNC * ilC + k] ^ CrAnMask;
                        }
                        if (ComplexOperatorAddLader
                            (O, det, ilA + ilB + ilC,
                             A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Re * C.l[ilC].V[iNC].Re +
                             A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Im * C.l[ilC].V[iNC].Im,
                             A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Im * C.l[ilC].V[iNC].Re -
                             A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Re * C.l[ilC].V[iNC].Im))
                        {
                          printf("ComplexOperatorAddLader failed in OperatorMultiply3\n");
                          fflush(stdout);
                          return 1;
                        }
                      }
                    }
                  }
                }
                else
                {               /* A real B complex C real */
                  for (iNA = 0; iNA < A.l[ilA].N; iNA++)
                  {
                    for (k = 0; k < ilA; k++)
                    {
                      det[k] = A.l[ilA].o[iNA * ilA + k];
                    }
                    for (iNB = 0; iNB < B.l[ilB].N; iNB++)
                    {
                      for (k = 0; k < ilB; k++)
                      {
                        det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                      }
                      for (iNC = 0; iNC < C.l[ilC].N; iNC++)
                      {
                        for (k = 0; k < ilC; k++)
                        {
                          det[ilC - 1 - k + ilA + ilB] = C.l[ilC].o[iNC * ilC + k] ^ CrAnMask;
                        }
                        if (ComplexOperatorAddLader
                            (O, det, ilA + ilB + ilC, A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Re * C.l[ilC].VRe[iNC],
                             A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Im * C.l[ilC].VRe[iNC]))
                        {
                          printf("ComplexOperatorAddLader failed in OperatorMultiply3\n");
                          fflush(stdout);
                          return 1;
                        }
                      }
                    }
                  }
                }
              }
              else
              {                 /* A and B real */
                if (C.l[ilC].QComplex)
                {               /* A and B real C complex */
                  for (iNA = 0; iNA < A.l[ilA].N; iNA++)
                  {
                    for (k = 0; k < ilA; k++)
                    {
                      det[k] = A.l[ilA].o[iNA * ilA + k];
                    }
                    for (iNB = 0; iNB < B.l[ilB].N; iNB++)
                    {
                      for (k = 0; k < ilB; k++)
                      {
                        det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                      }
                      for (iNC = 0; iNC < C.l[ilC].N; iNC++)
                      {
                        for (k = 0; k < ilC; k++)
                        {
                          det[ilC - 1 - k + ilA + ilB] = C.l[ilC].o[iNC * ilC + k] ^ CrAnMask;
                        }
                        if (ComplexOperatorAddLader
                            (O, det, ilA + ilB + ilC, A.l[ilA].VRe[iNA] * B.l[ilB].VRe[iNB] * C.l[ilC].V[iNC].Re,
                             -A.l[ilA].VRe[iNA] * B.l[ilB].VRe[iNB] * C.l[ilC].V[iNC].Im))
                        {
                          printf("ComplexOperatorAddLader failed in OperatorMultiply3\n");
                          fflush(stdout);
                          return 1;
                        }
                      }
                    }
                  }
                }
                else
                {               /* A, B and C real (O complex) */
                  for (iNA = 0; iNA < A.l[ilA].N; iNA++)
                  {
                    for (k = 0; k < ilA; k++)
                    {
                      det[k] = A.l[ilA].o[iNA * ilA + k];
                    }
                    for (iNB = 0; iNB < B.l[ilB].N; iNB++)
                    {
                      for (k = 0; k < ilB; k++)
                      {
                        det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                      }
                      for (iNC = 0; iNC < C.l[ilC].N; iNC++)
                      {
                        for (k = 0; k < ilC; k++)
                        {
                          det[ilC - 1 - k + ilA + ilB] = C.l[ilC].o[iNC * ilC + k] ^ CrAnMask;
                        }
                        if (ComplexOperatorAddLader
                            (O, det, ilA + ilB + ilC, A.l[ilA].VRe[iNA] * B.l[ilB].VRe[iNB] * C.l[ilC].VRe[iNC], 0))
                        {
                          printf("ComplexOperatorAddLader failed in OperatorMultiply3\n");
                          fflush(stdout);
                          return 1;
                        }
                      }
                    }
                  }
                }
              }
            }
          }
          else
          {                     /* All real */
            for (iNA = 0; iNA < A.l[ilA].N; iNA++)
            {
              for (k = 0; k < ilA; k++)
              {
                det[k] = A.l[ilA].o[iNA * ilA + k];
              }
              for (iNB = 0; iNB < B.l[ilB].N; iNB++)
              {
                for (k = 0; k < ilB; k++)
                {
                  det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                }
                for (iNC = 0; iNC < C.l[ilC].N; iNC++)
                {
                  for (k = 0; k < ilC; k++)
                  {
                    det[ilC - 1 - k + ilA + ilB] = C.l[ilC].o[iNC * ilC + k] ^ CrAnMask;
                  }
                  if (RealOperatorAddLader
                      (O, det, ilA + ilB + ilC, A.l[ilA].VRe[iNA] * B.l[ilB].VRe[iNB] * C.l[ilC].VRe[iNC]))
                  {
                    printf("RealOperatorAddLader failed in OperatorMultiply3\n");
                    fflush(stdout);
                    return 1;
                  }
                }
              }
            }
          }
        }
      }
    }
    /* End of NND */
  }
  else if ((Aset == 'N') && (Bset == 'D') && (Cset == 'D'))
  {
    printf("in OperatorMultiply3 Still to do\n");
    return 1;
  }
  else if ((Aset == 'D') && (Bset == 'N') && (Cset == 'D'))
  {
    printf("in OperatorMultiply3 Still to do\n");
    return 1;
  }
  else if ((Aset == 'D') && (Bset == 'D') && (Cset == 'D'))
  {
    printf("in OperatorMultiply3 Still to do\n");
    return 1;
  }
  else
  {
    printf("in OperatorMultiply3 unrecognized input for Aset Bset and or Cset\n");
    return 1;
  }
  return 0;
}




/*  Depending on Aset, Bset and Cset
 *    O <- O + A^'Aset' B^'Bset' C^'Cset'
 *    whereby set can be
 *    'N' normal
 *    'C' conjugated (a^dag remains a^dag, prefactors become complex conjugated)
 *    'T' transposed (a^dag becomes a, prefactors remain)
 *    'D' conjugate transposed (a^dag becomes a, prefactors become complex conjugated)
 */
uint32_t OperatorMultiply3NormalOrder(OperatorType * const O, OperatorType const C, OperatorType const B,
                                      OperatorType const A, char const Cset, char const Bset, char const Aset,
                                      uint32_t(*QDetNotOrdered) (uint16_t, uint16_t, uint32_t))
{
  uint32_t ilA, ilB, ilC, iNA, iNB, iNC, k;
  uint16_t det[A.MaxLength + B.MaxLength + C.MaxLength];

  if (O->MaxLength != A.MaxLength + B.MaxLength + C.MaxLength)
  {
    printf("Length of output operator does not match in OperatorMultiply3NormalOrder\n");
    fflush(stdout);
    return 1;
  }

  if ((Aset == 'N') && (Bset == 'N') && (Cset == 'N'))
  {
    for (ilA = 0; ilA <= A.MaxLength; ilA++)
    {
      for (ilB = 0; ilB <= B.MaxLength; ilB++)
      {
        for (ilC = 0; ilC <= C.MaxLength; ilC++)
        {
          if ((O->l[ilA + ilB].QComplex) || (A.l[ilA].QComplex) || (B.l[ilB].QComplex) || (C.l[ilC].QComplex))
          {                     /* at least one complex */
            if (!(O->l[ilA + ilB + ilC].QComplex))
              if (MakeRealOperatorFixedLengthComplex(O, ilA + ilB + ilC))
              {
                printf("MakeRealOperatorFixedLengthComplex at length=%i failed in OperatorMultiply3\n",
                       ilA + ilB + ilC);
                fflush(stdout);
                return 1;
              }
            if (A.l[ilA].QComplex)
            {                   /* A complex */
              if (B.l[ilB].QComplex)
              {                 /* A and B complex */
                if (C.l[ilC].QComplex)
                {               /* A and B and C complex */
                  for (iNA = 0; iNA < A.l[ilA].N; iNA++)
                  {
                    for (k = 0; k < ilA; k++)
                    {
                      det[k] = A.l[ilA].o[iNA * ilA + k];
                    }
                    for (iNB = 0; iNB < B.l[ilB].N; iNB++)
                    {
                      for (k = 0; k < ilB; k++)
                      {
                        det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                      }
                      for (iNC = 0; iNC < C.l[ilC].N; iNC++)
                      {
                        for (k = 0; k < ilC; k++)
                        {
                          det[k + ilA + ilB] = C.l[ilC].o[iNC * ilC + k];
                        }
                        if (ComplexOperatorAddLaderNormalOrder
                            (O, det, ilA + ilB + ilC,
                             A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Re * C.l[ilC].V[iNC].Re -
                             A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Im * C.l[ilC].V[iNC].Im -
                             A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Re * C.l[ilC].V[iNC].Im -
                             A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Im * C.l[ilC].V[iNC].Re,
                             A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Re * C.l[ilC].V[iNC].Re +
                             A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Im * C.l[ilC].V[iNC].Re +
                             A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Re * C.l[ilC].V[iNC].Im -
                             A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Im * C.l[ilC].V[iNC].Im, QDetNotOrdered))
                        {
                          printf("ComplexOperatorAddLader failed in OperatorMultiply3\n");
                          fflush(stdout);
                          return 1;
                        }
                      }
                    }
                  }
                }
                else
                {               /* A and B complex, C real */
                  for (iNA = 0; iNA < A.l[ilA].N; iNA++)
                  {
                    for (k = 0; k < ilA; k++)
                    {
                      det[k] = A.l[ilA].o[iNA * ilA + k];
                    }
                    for (iNB = 0; iNB < B.l[ilB].N; iNB++)
                    {
                      for (k = 0; k < ilB; k++)
                      {
                        det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                      }
                      for (iNC = 0; iNC < C.l[ilC].N; iNC++)
                      {
                        for (k = 0; k < ilC; k++)
                        {
                          det[k + ilA + ilB] = C.l[ilC].o[iNC * ilC + k];
                        }
                        if (ComplexOperatorAddLaderNormalOrder
                            (O, det, ilA + ilB + ilC,
                             A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Re * C.l[ilC].VRe[iNC] -
                             A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Im * C.l[ilC].VRe[iNC],
                             A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Re * C.l[ilC].VRe[iNC] +
                             A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Im * C.l[ilC].VRe[iNC], QDetNotOrdered))
                        {
                          printf("ComplexOperatorAddLader failed in OperatorMultiply3\n");
                          fflush(stdout);
                          return 1;
                        }
                      }
                    }
                  }
                }
              }
              else
              {                 /* A complex B real */
                if (C.l[ilC].QComplex)
                {               /* A complex B real C complex */
                  for (iNA = 0; iNA < A.l[ilA].N; iNA++)
                  {
                    for (k = 0; k < ilA; k++)
                    {
                      det[k] = A.l[ilA].o[iNA * ilA + k];
                    }
                    for (iNB = 0; iNB < B.l[ilB].N; iNB++)
                    {
                      for (k = 0; k < ilB; k++)
                      {
                        det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                      }
                      for (iNC = 0; iNC < C.l[ilC].N; iNC++)
                      {
                        for (k = 0; k < ilC; k++)
                        {
                          det[k + ilA + ilB] = C.l[ilC].o[iNC * ilC + k];
                        }
                        if (ComplexOperatorAddLaderNormalOrder
                            (O, det, ilA + ilB + ilC,
                             A.l[ilA].V[iNA].Re * B.l[ilB].VRe[iNB] * C.l[ilC].V[iNC].Re -
                             A.l[ilA].V[iNA].Im * B.l[ilB].VRe[iNB] * C.l[ilC].V[iNC].Im,
                             A.l[ilA].V[iNA].Im * B.l[ilB].VRe[iNB] * C.l[ilC].V[iNC].Re +
                             A.l[ilA].V[iNA].Re * B.l[ilB].VRe[iNB] * C.l[ilC].V[iNC].Im, QDetNotOrdered))
                        {
                          printf("ComplexOperatorAddLader failed in OperatorMultiply3\n");
                          fflush(stdout);
                          return 1;
                        }
                      }
                    }
                  }
                }
                else
                {               /* A complex B real C real */
                  for (iNA = 0; iNA < A.l[ilA].N; iNA++)
                  {
                    for (k = 0; k < ilA; k++)
                    {
                      det[k] = A.l[ilA].o[iNA * ilA + k];
                    }
                    for (iNB = 0; iNB < B.l[ilB].N; iNB++)
                    {
                      for (k = 0; k < ilB; k++)
                      {
                        det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                      }
                      for (iNC = 0; iNC < C.l[ilC].N; iNC++)
                      {
                        for (k = 0; k < ilC; k++)
                        {
                          det[k + ilA + ilB] = C.l[ilC].o[iNC * ilC + k];
                        }
                        if (ComplexOperatorAddLaderNormalOrder
                            (O, det, ilA + ilB + ilC, A.l[ilA].V[iNA].Re * B.l[ilB].VRe[iNB] * C.l[ilC].VRe[iNC],
                             A.l[ilA].V[iNA].Im * B.l[ilB].VRe[iNB] * C.l[ilC].VRe[iNC], QDetNotOrdered))
                        {
                          printf("ComplexOperatorAddLader failed in OperatorMultiply3\n");
                          fflush(stdout);
                          return 1;
                        }
                      }
                    }
                  }
                }
              }
            }
            else
            {                   /* A real */
              if (B.l[ilB].QComplex)
              {                 /* A real and B complex */
                if (C.l[ilC].QComplex)
                {               /* A real B complex C complex */
                  for (iNA = 0; iNA < A.l[ilA].N; iNA++)
                  {
                    for (k = 0; k < ilA; k++)
                    {
                      det[k] = A.l[ilA].o[iNA * ilA + k];
                    }
                    for (iNB = 0; iNB < B.l[ilB].N; iNB++)
                    {
                      for (k = 0; k < ilB; k++)
                      {
                        det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                      }
                      for (iNC = 0; iNC < C.l[ilC].N; iNC++)
                      {
                        for (k = 0; k < ilC; k++)
                        {
                          det[k + ilA + ilB] = C.l[ilC].o[iNC * ilC + k];
                        }
                        if (ComplexOperatorAddLaderNormalOrder
                            (O, det, ilA + ilB + ilC,
                             A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Re * C.l[ilC].V[iNC].Re -
                             A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Im * C.l[ilC].V[iNC].Im,
                             A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Im * C.l[ilC].V[iNC].Re +
                             A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Re * C.l[ilC].V[iNC].Im, QDetNotOrdered))
                        {
                          printf("ComplexOperatorAddLader failed in OperatorMultiply3\n");
                          fflush(stdout);
                          return 1;
                        }
                      }
                    }
                  }
                }
                else
                {               /* A real B complex C real */
                  for (iNA = 0; iNA < A.l[ilA].N; iNA++)
                  {
                    for (k = 0; k < ilA; k++)
                    {
                      det[k] = A.l[ilA].o[iNA * ilA + k];
                    }
                    for (iNB = 0; iNB < B.l[ilB].N; iNB++)
                    {
                      for (k = 0; k < ilB; k++)
                      {
                        det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                      }
                      for (iNC = 0; iNC < C.l[ilC].N; iNC++)
                      {
                        for (k = 0; k < ilC; k++)
                        {
                          det[k + ilA + ilB] = C.l[ilC].o[iNC * ilC + k];
                        }
                        if (ComplexOperatorAddLaderNormalOrder
                            (O, det, ilA + ilB + ilC, A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Re * C.l[ilC].VRe[iNC],
                             A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Im * C.l[ilC].VRe[iNC], QDetNotOrdered))
                        {
                          printf("ComplexOperatorAddLader failed in OperatorMultiply3\n");
                          fflush(stdout);
                          return 1;
                        }
                      }
                    }
                  }
                }
              }
              else
              {                 /* A and B real */
                if (C.l[ilC].QComplex)
                {               /* A and B real C complex */
                  for (iNA = 0; iNA < A.l[ilA].N; iNA++)
                  {
                    for (k = 0; k < ilA; k++)
                    {
                      det[k] = A.l[ilA].o[iNA * ilA + k];
                    }
                    for (iNB = 0; iNB < B.l[ilB].N; iNB++)
                    {
                      for (k = 0; k < ilB; k++)
                      {
                        det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                      }
                      for (iNC = 0; iNC < C.l[ilC].N; iNC++)
                      {
                        for (k = 0; k < ilC; k++)
                        {
                          det[k + ilA + ilB] = C.l[ilC].o[iNC * ilC + k];
                        }
                        if (ComplexOperatorAddLaderNormalOrder
                            (O, det, ilA + ilB + ilC, A.l[ilA].VRe[iNA] * B.l[ilB].VRe[iNB] * C.l[ilC].V[iNC].Re,
                             A.l[ilA].VRe[iNA] * B.l[ilB].VRe[iNB] * C.l[ilC].V[iNC].Im, QDetNotOrdered))
                        {
                          printf("ComplexOperatorAddLader failed in OperatorMultiply3\n");
                          fflush(stdout);
                          return 1;
                        }
                      }
                    }
                  }
                }
                else
                {               /* A, B and C real (O complex) */
                  for (iNA = 0; iNA < A.l[ilA].N; iNA++)
                  {
                    for (k = 0; k < ilA; k++)
                    {
                      det[k] = A.l[ilA].o[iNA * ilA + k];
                    }
                    for (iNB = 0; iNB < B.l[ilB].N; iNB++)
                    {
                      for (k = 0; k < ilB; k++)
                      {
                        det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                      }
                      for (iNC = 0; iNC < C.l[ilC].N; iNC++)
                      {
                        for (k = 0; k < ilC; k++)
                        {
                          det[k + ilA + ilB] = C.l[ilC].o[iNC * ilC + k];
                        }
                        if (ComplexOperatorAddLaderNormalOrder
                            (O, det, ilA + ilB + ilC, A.l[ilA].VRe[iNA] * B.l[ilB].VRe[iNB] * C.l[ilC].VRe[iNC], 0.0,
                             QDetNotOrdered))
                        {
                          printf("ComplexOperatorAddLader failed in OperatorMultiply3\n");
                          fflush(stdout);
                          return 1;
                        }
                      }
                    }
                  }
                }
              }
            }
          }
          else
          {                     /* All real */
            for (iNA = 0; iNA < A.l[ilA].N; iNA++)
            {
              for (k = 0; k < ilA; k++)
              {
                det[k] = A.l[ilA].o[iNA * ilA + k];
              }
              for (iNB = 0; iNB < B.l[ilB].N; iNB++)
              {
                for (k = 0; k < ilB; k++)
                {
                  det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                }
                for (iNC = 0; iNC < C.l[ilC].N; iNC++)
                {
                  for (k = 0; k < ilC; k++)
                  {
                    det[k + ilA + ilB] = C.l[ilC].o[iNC * ilC + k];
                  }
                  if (RealOperatorAddLaderNormalOrder
                      (O, det, ilA + ilB + ilC, A.l[ilA].VRe[iNA] * B.l[ilB].VRe[iNB] * C.l[ilC].VRe[iNC],
                       QDetNotOrdered))
                  {
                    printf("RealOperatorAddLader failed in OperatorMultiply3\n");
                    fflush(stdout);
                    return 1;
                  }
                }
              }
            }
          }
        }
      }
    }
    /* End of NNN */
  }
  else if ((Aset == 'N') && (Bset == 'D') && (Cset == 'N'))
  {
    printf("in OperatorMultiply3 Still to do\n");
    return 1;
  }
  else if ((Aset == 'D') && (Bset == 'N') && (Cset == 'N'))
  {
    for (ilA = 0; ilA <= A.MaxLength; ilA++)
    {
      for (ilB = 0; ilB <= B.MaxLength; ilB++)
      {
        for (ilC = 0; ilC <= C.MaxLength; ilC++)
        {
          if ((O->l[ilA + ilB].QComplex) || (A.l[ilA].QComplex) || (B.l[ilB].QComplex) || (C.l[ilC].QComplex))
          {                     /* at least one complex */
            if (!(O->l[ilA + ilB + ilC].QComplex))
              if (MakeRealOperatorFixedLengthComplex(O, ilA + ilB + ilC))
              {
                printf("MakeRealOperatorFixedLengthComplex at length=%i failed in OperatorMultiply3\n",
                       ilA + ilB + ilC);
                fflush(stdout);
                return 1;
              }
            if (A.l[ilA].QComplex)
            {                   /* A complex */
              if (B.l[ilB].QComplex)
              {                 /* A and B complex */
                if (C.l[ilC].QComplex)
                {               /* A and B and C complex */
                  for (iNA = 0; iNA < A.l[ilA].N; iNA++)
                  {
                    for (k = 0; k < ilA; k++)
                    {
                      det[ilA - 1 - k] = A.l[ilA].o[iNA * ilA + k] ^ CrAnMask;
                    }
                    for (iNB = 0; iNB < B.l[ilB].N; iNB++)
                    {
                      for (k = 0; k < ilB; k++)
                      {
                        det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                      }
                      for (iNC = 0; iNC < C.l[ilC].N; iNC++)
                      {
                        for (k = 0; k < ilC; k++)
                        {
                          det[k + ilA + ilB] = C.l[ilC].o[iNC * ilC + k];
                        }
                        if (ComplexOperatorAddLaderNormalOrder
                            (O, det, ilA + ilB + ilC,
                             A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Re * C.l[ilC].V[iNC].Re -
                             A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Im * C.l[ilC].V[iNC].Im +
                             A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Re * C.l[ilC].V[iNC].Im +
                             A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Im * C.l[ilC].V[iNC].Re,
                             -A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Re * C.l[ilC].V[iNC].Re +
                             A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Im * C.l[ilC].V[iNC].Re +
                             A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Re * C.l[ilC].V[iNC].Im +
                             A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Im * C.l[ilC].V[iNC].Im, QDetNotOrdered))
                        {
                          printf("ComplexOperatorAddLader failed in OperatorMultiply3\n");
                          fflush(stdout);
                          return 1;
                        }
                      }
                    }
                  }
                }
                else
                {               /* A and B complex, C real */
                  for (iNA = 0; iNA < A.l[ilA].N; iNA++)
                  {
                    for (k = 0; k < ilA; k++)
                    {
                      det[ilA - 1 - k] = A.l[ilA].o[iNA * ilA + k] ^ CrAnMask;
                    }
                    for (iNB = 0; iNB < B.l[ilB].N; iNB++)
                    {
                      for (k = 0; k < ilB; k++)
                      {
                        det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                      }
                      for (iNC = 0; iNC < C.l[ilC].N; iNC++)
                      {
                        for (k = 0; k < ilC; k++)
                        {
                          det[k + ilA + ilB] = C.l[ilC].o[iNC * ilC + k];
                        }
                        if (ComplexOperatorAddLaderNormalOrder
                            (O, det, ilA + ilB + ilC,
                             A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Re * C.l[ilC].VRe[iNC] +
                             A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Im * C.l[ilC].VRe[iNC],
                             -A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Re * C.l[ilC].VRe[iNC] +
                             A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Im * C.l[ilC].VRe[iNC], QDetNotOrdered))
                        {
                          printf("ComplexOperatorAddLader failed in OperatorMultiply3\n");
                          fflush(stdout);
                          return 1;
                        }
                      }
                    }
                  }
                }
              }
              else
              {                 /* A complex B real */
                if (C.l[ilC].QComplex)
                {               /* A complex B real C complex */
                  for (iNA = 0; iNA < A.l[ilA].N; iNA++)
                  {
                    for (k = 0; k < ilA; k++)
                    {
                      det[ilA - 1 - k] = A.l[ilA].o[iNA * ilA + k] ^ CrAnMask;
                    }
                    for (iNB = 0; iNB < B.l[ilB].N; iNB++)
                    {
                      for (k = 0; k < ilB; k++)
                      {
                        det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                      }
                      for (iNC = 0; iNC < C.l[ilC].N; iNC++)
                      {
                        for (k = 0; k < ilC; k++)
                        {
                          det[k + ilA + ilB] = C.l[ilC].o[iNC * ilC + k];
                        }
                        if (ComplexOperatorAddLaderNormalOrder
                            (O, det, ilA + ilB + ilC,
                             A.l[ilA].V[iNA].Re * B.l[ilB].VRe[iNB] * C.l[ilC].V[iNC].Re +
                             A.l[ilA].V[iNA].Im * B.l[ilB].VRe[iNB] * C.l[ilC].V[iNC].Im,
                             -A.l[ilA].V[iNA].Im * B.l[ilB].VRe[iNB] * C.l[ilC].V[iNC].Re +
                             A.l[ilA].V[iNA].Re * B.l[ilB].VRe[iNB] * C.l[ilC].V[iNC].Im, QDetNotOrdered))
                        {
                          printf("ComplexOperatorAddLader failed in OperatorMultiply3\n");
                          fflush(stdout);
                          return 1;
                        }
                      }
                    }
                  }
                }
                else
                {               /* A complex B real C real */
                  for (iNA = 0; iNA < A.l[ilA].N; iNA++)
                  {
                    for (k = 0; k < ilA; k++)
                    {
                      det[ilA - 1 - k] = A.l[ilA].o[iNA * ilA + k] ^ CrAnMask;
                    }
                    for (iNB = 0; iNB < B.l[ilB].N; iNB++)
                    {
                      for (k = 0; k < ilB; k++)
                      {
                        det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                      }
                      for (iNC = 0; iNC < C.l[ilC].N; iNC++)
                      {
                        for (k = 0; k < ilC; k++)
                        {
                          det[k + ilA + ilB] = C.l[ilC].o[iNC * ilC + k];
                        }
                        if (ComplexOperatorAddLaderNormalOrder
                            (O, det, ilA + ilB + ilC, A.l[ilA].V[iNA].Re * B.l[ilB].VRe[iNB] * C.l[ilC].VRe[iNC],
                             -A.l[ilA].V[iNA].Im * B.l[ilB].VRe[iNB] * C.l[ilC].VRe[iNC], QDetNotOrdered))
                        {
                          printf("ComplexOperatorAddLader failed in OperatorMultiply3\n");
                          fflush(stdout);
                          return 1;
                        }
                      }
                    }
                  }
                }
              }
            }
            else
            {                   /* A real */
              if (B.l[ilB].QComplex)
              {                 /* A real and B complex */
                if (C.l[ilC].QComplex)
                {               /* A real B complex C complex */
                  for (iNA = 0; iNA < A.l[ilA].N; iNA++)
                  {
                    for (k = 0; k < ilA; k++)
                    {
                      det[ilA - 1 - k] = A.l[ilA].o[iNA * ilA + k] ^ CrAnMask;
                    }
                    for (iNB = 0; iNB < B.l[ilB].N; iNB++)
                    {
                      for (k = 0; k < ilB; k++)
                      {
                        det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                      }
                      for (iNC = 0; iNC < C.l[ilC].N; iNC++)
                      {
                        for (k = 0; k < ilC; k++)
                        {
                          det[k + ilA + ilB] = C.l[ilC].o[iNC * ilC + k];
                        }
                        if (ComplexOperatorAddLaderNormalOrder
                            (O, det, ilA + ilB + ilC,
                             A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Re * C.l[ilC].V[iNC].Re -
                             A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Im * C.l[ilC].V[iNC].Im,
                             A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Im * C.l[ilC].V[iNC].Re +
                             A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Re * C.l[ilC].V[iNC].Im, QDetNotOrdered))
                        {
                          printf("ComplexOperatorAddLader failed in OperatorMultiply3\n");
                          fflush(stdout);
                          return 1;
                        }
                      }
                    }
                  }
                }
                else
                {               /* A real B complex C real */
                  for (iNA = 0; iNA < A.l[ilA].N; iNA++)
                  {
                    for (k = 0; k < ilA; k++)
                    {
                      det[ilA - 1 - k] = A.l[ilA].o[iNA * ilA + k] ^ CrAnMask;
                    }
                    for (iNB = 0; iNB < B.l[ilB].N; iNB++)
                    {
                      for (k = 0; k < ilB; k++)
                      {
                        det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                      }
                      for (iNC = 0; iNC < C.l[ilC].N; iNC++)
                      {
                        for (k = 0; k < ilC; k++)
                        {
                          det[k + ilA + ilB] = C.l[ilC].o[iNC * ilC + k];
                        }
                        if (ComplexOperatorAddLaderNormalOrder
                            (O, det, ilA + ilB + ilC, A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Re * C.l[ilC].VRe[iNC],
                             A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Im * C.l[ilC].VRe[iNC], QDetNotOrdered))
                        {
                          printf("ComplexOperatorAddLader failed in OperatorMultiply3\n");
                          fflush(stdout);
                          return 1;
                        }
                      }
                    }
                  }
                }
              }
              else
              {                 /* A and B real */
                if (C.l[ilC].QComplex)
                {               /* A and B real C complex */
                  for (iNA = 0; iNA < A.l[ilA].N; iNA++)
                  {
                    for (k = 0; k < ilA; k++)
                    {
                      det[ilA - 1 - k] = A.l[ilA].o[iNA * ilA + k] ^ CrAnMask;
                    }
                    for (iNB = 0; iNB < B.l[ilB].N; iNB++)
                    {
                      for (k = 0; k < ilB; k++)
                      {
                        det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                      }
                      for (iNC = 0; iNC < C.l[ilC].N; iNC++)
                      {
                        for (k = 0; k < ilC; k++)
                        {
                          det[k + ilA + ilB] = C.l[ilC].o[iNC * ilC + k];
                        }
                        if (ComplexOperatorAddLaderNormalOrder
                            (O, det, ilA + ilB + ilC, A.l[ilA].VRe[iNA] * B.l[ilB].VRe[iNB] * C.l[ilC].V[iNC].Re,
                             A.l[ilA].VRe[iNA] * B.l[ilB].VRe[iNB] * C.l[ilC].V[iNC].Im, QDetNotOrdered))
                        {
                          printf("ComplexOperatorAddLader failed in OperatorMultiply3\n");
                          fflush(stdout);
                          return 1;
                        }
                      }
                    }
                  }
                }
                else
                {               /* A, B and C real (O complex) */
                  for (iNA = 0; iNA < A.l[ilA].N; iNA++)
                  {
                    for (k = 0; k < ilA; k++)
                    {
                      det[ilA - 1 - k] = A.l[ilA].o[iNA * ilA + k] ^ CrAnMask;
                    }
                    for (iNB = 0; iNB < B.l[ilB].N; iNB++)
                    {
                      for (k = 0; k < ilB; k++)
                      {
                        det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                      }
                      for (iNC = 0; iNC < C.l[ilC].N; iNC++)
                      {
                        for (k = 0; k < ilC; k++)
                        {
                          det[k + ilA + ilB] = C.l[ilC].o[iNC * ilC + k];
                        }
                        if (ComplexOperatorAddLaderNormalOrder
                            (O, det, ilA + ilB + ilC, A.l[ilA].VRe[iNA] * B.l[ilB].VRe[iNB] * C.l[ilC].VRe[iNC], 0.0,
                             QDetNotOrdered))
                        {
                          printf("ComplexOperatorAddLader failed in OperatorMultiply3\n");
                          fflush(stdout);
                          return 1;
                        }
                      }
                    }
                  }
                }
              }
            }
          }
          else
          {                     /* All real */
            for (iNA = 0; iNA < A.l[ilA].N; iNA++)
            {
              for (k = 0; k < ilA; k++)
              {
                det[ilA - 1 - k] = A.l[ilA].o[iNA * ilA + k] ^ CrAnMask;
              }
              for (iNB = 0; iNB < B.l[ilB].N; iNB++)
              {
                for (k = 0; k < ilB; k++)
                {
                  det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                }
                for (iNC = 0; iNC < C.l[ilC].N; iNC++)
                {
                  for (k = 0; k < ilC; k++)
                  {
                    det[k + ilA + ilB] = C.l[ilC].o[iNC * ilC + k];
                  }
                  if (RealOperatorAddLaderNormalOrder
                      (O, det, ilA + ilB + ilC, A.l[ilA].VRe[iNA] * B.l[ilB].VRe[iNB] * C.l[ilC].VRe[iNC],
                       QDetNotOrdered))
                  {
                    printf("RealOperatorAddLader failed in OperatorMultiply3\n");
                    fflush(stdout);
                    return 1;
                  }
                }
              }
            }
          }
        }
      }
    }
    /* End of DNN */
  }
  else if ((Aset == 'D') && (Bset == 'D') && (Cset == 'N'))
  {
    printf("in OperatorMultiply3 Still to do\n");
    return 1;
  }
  else if ((Aset == 'N') && (Bset == 'N') && (Cset == 'D'))
  {
    for (ilA = 0; ilA <= A.MaxLength; ilA++)
    {
      for (ilB = 0; ilB <= B.MaxLength; ilB++)
      {
        for (ilC = 0; ilC <= C.MaxLength; ilC++)
        {
          if ((O->l[ilA + ilB].QComplex) || (A.l[ilA].QComplex) || (B.l[ilB].QComplex) || (C.l[ilC].QComplex))
          {                     /* at least one complex */
            if (!(O->l[ilA + ilB + ilC].QComplex))
              if (MakeRealOperatorFixedLengthComplex(O, ilA + ilB + ilC))
              {
                printf("MakeRealOperatorFixedLengthComplex at length=%i failed in OperatorMultiply3\n",
                       ilA + ilB + ilC);
                fflush(stdout);
                return 1;
              }
            if (A.l[ilA].QComplex)
            {                   /* A complex */
              if (B.l[ilB].QComplex)
              {                 /* A and B complex */
                if (C.l[ilC].QComplex)
                {               /* A and B and C complex */
                  for (iNA = 0; iNA < A.l[ilA].N; iNA++)
                  {
                    for (k = 0; k < ilA; k++)
                    {
                      det[k] = A.l[ilA].o[iNA * ilA + k];
                    }
                    for (iNB = 0; iNB < B.l[ilB].N; iNB++)
                    {
                      for (k = 0; k < ilB; k++)
                      {
                        det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                      }
                      for (iNC = 0; iNC < C.l[ilC].N; iNC++)
                      {
                        for (k = 0; k < ilC; k++)
                        {
                          det[ilC - 1 - k + ilA + ilB] = C.l[ilC].o[iNC * ilC + k] ^ CrAnMask;
                        }
                        if (ComplexOperatorAddLaderNormalOrder
                            (O, det, ilA + ilB + ilC,
                             A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Re * C.l[ilC].V[iNC].Re +
                             A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Im * C.l[ilC].V[iNC].Im +
                             A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Re * C.l[ilC].V[iNC].Im -
                             A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Im * C.l[ilC].V[iNC].Re,
                             A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Re * C.l[ilC].V[iNC].Re +
                             A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Im * C.l[ilC].V[iNC].Re -
                             A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Re * C.l[ilC].V[iNC].Im +
                             A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Im * C.l[ilC].V[iNC].Im, QDetNotOrdered))
                        {
                          printf("ComplexOperatorAddLader failed in OperatorMultiply3\n");
                          fflush(stdout);
                          return 1;
                        }
                      }
                    }
                  }
                }
                else
                {               /* A and B complex, C real */
                  for (iNA = 0; iNA < A.l[ilA].N; iNA++)
                  {
                    for (k = 0; k < ilA; k++)
                    {
                      det[k] = A.l[ilA].o[iNA * ilA + k];
                    }
                    for (iNB = 0; iNB < B.l[ilB].N; iNB++)
                    {
                      for (k = 0; k < ilB; k++)
                      {
                        det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                      }
                      for (iNC = 0; iNC < C.l[ilC].N; iNC++)
                      {
                        for (k = 0; k < ilC; k++)
                        {
                          det[ilC - 1 - k + ilA + ilB] = C.l[ilC].o[iNC * ilC + k] ^ CrAnMask;
                        }
                        if (ComplexOperatorAddLaderNormalOrder
                            (O, det, ilA + ilB + ilC,
                             A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Re * C.l[ilC].VRe[iNC] -
                             A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Im * C.l[ilC].VRe[iNC],
                             A.l[ilA].V[iNA].Im * B.l[ilB].V[iNB].Re * C.l[ilC].VRe[iNC] +
                             A.l[ilA].V[iNA].Re * B.l[ilB].V[iNB].Im * C.l[ilC].VRe[iNC], QDetNotOrdered))
                        {
                          printf("ComplexOperatorAddLader failed in OperatorMultiply3\n");
                          fflush(stdout);
                          return 1;
                        }
                      }
                    }
                  }
                }
              }
              else
              {                 /* A complex B real */
                if (C.l[ilC].QComplex)
                {               /* A complex B real C complex */
                  for (iNA = 0; iNA < A.l[ilA].N; iNA++)
                  {
                    for (k = 0; k < ilA; k++)
                    {
                      det[k] = A.l[ilA].o[iNA * ilA + k];
                    }
                    for (iNB = 0; iNB < B.l[ilB].N; iNB++)
                    {
                      for (k = 0; k < ilB; k++)
                      {
                        det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                      }
                      for (iNC = 0; iNC < C.l[ilC].N; iNC++)
                      {
                        for (k = 0; k < ilC; k++)
                        {
                          det[ilC - 1 - k + ilA + ilB] = C.l[ilC].o[iNC * ilC + k] ^ CrAnMask;
                        }
                        if (ComplexOperatorAddLaderNormalOrder
                            (O, det, ilA + ilB + ilC,
                             A.l[ilA].V[iNA].Re * B.l[ilB].VRe[iNB] * C.l[ilC].V[iNC].Re +
                             A.l[ilA].V[iNA].Im * B.l[ilB].VRe[iNB] * C.l[ilC].V[iNC].Im,
                             A.l[ilA].V[iNA].Im * B.l[ilB].VRe[iNB] * C.l[ilC].V[iNC].Re -
                             A.l[ilA].V[iNA].Re * B.l[ilB].VRe[iNB] * C.l[ilC].V[iNC].Im, QDetNotOrdered))
                        {
                          printf("ComplexOperatorAddLader failed in OperatorMultiply3\n");
                          fflush(stdout);
                          return 1;
                        }
                      }
                    }
                  }
                }
                else
                {               /* A complex B real C real */
                  for (iNA = 0; iNA < A.l[ilA].N; iNA++)
                  {
                    for (k = 0; k < ilA; k++)
                    {
                      det[k] = A.l[ilA].o[iNA * ilA + k];
                    }
                    for (iNB = 0; iNB < B.l[ilB].N; iNB++)
                    {
                      for (k = 0; k < ilB; k++)
                      {
                        det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                      }
                      for (iNC = 0; iNC < C.l[ilC].N; iNC++)
                      {
                        for (k = 0; k < ilC; k++)
                        {
                          det[ilC - 1 - k + ilA + ilB] = C.l[ilC].o[iNC * ilC + k] ^ CrAnMask;
                        }
                        if (ComplexOperatorAddLaderNormalOrder
                            (O, det, ilA + ilB + ilC, A.l[ilA].V[iNA].Re * B.l[ilB].VRe[iNB] * C.l[ilC].VRe[iNC],
                             A.l[ilA].V[iNA].Im * B.l[ilB].VRe[iNB] * C.l[ilC].VRe[iNC], QDetNotOrdered))
                        {
                          printf("ComplexOperatorAddLader failed in OperatorMultiply3\n");
                          fflush(stdout);
                          return 1;
                        }
                      }
                    }
                  }
                }
              }
            }
            else
            {                   /* A real */
              if (B.l[ilB].QComplex)
              {                 /* A real and B complex */
                if (C.l[ilC].QComplex)
                {               /* A real B complex C complex */
                  for (iNA = 0; iNA < A.l[ilA].N; iNA++)
                  {
                    for (k = 0; k < ilA; k++)
                    {
                      det[k] = A.l[ilA].o[iNA * ilA + k];
                    }
                    for (iNB = 0; iNB < B.l[ilB].N; iNB++)
                    {
                      for (k = 0; k < ilB; k++)
                      {
                        det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                      }
                      for (iNC = 0; iNC < C.l[ilC].N; iNC++)
                      {
                        for (k = 0; k < ilC; k++)
                        {
                          det[ilC - 1 - k + ilA + ilB] = C.l[ilC].o[iNC * ilC + k] ^ CrAnMask;
                        }
                        if (ComplexOperatorAddLaderNormalOrder
                            (O, det, ilA + ilB + ilC,
                             A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Re * C.l[ilC].V[iNC].Re +
                             A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Im * C.l[ilC].V[iNC].Im,
                             A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Im * C.l[ilC].V[iNC].Re -
                             A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Re * C.l[ilC].V[iNC].Im, QDetNotOrdered))
                        {
                          printf("ComplexOperatorAddLader failed in OperatorMultiply3\n");
                          fflush(stdout);
                          return 1;
                        }
                      }
                    }
                  }
                }
                else
                {               /* A real B complex C real */
                  for (iNA = 0; iNA < A.l[ilA].N; iNA++)
                  {
                    for (k = 0; k < ilA; k++)
                    {
                      det[k] = A.l[ilA].o[iNA * ilA + k];
                    }
                    for (iNB = 0; iNB < B.l[ilB].N; iNB++)
                    {
                      for (k = 0; k < ilB; k++)
                      {
                        det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                      }
                      for (iNC = 0; iNC < C.l[ilC].N; iNC++)
                      {
                        for (k = 0; k < ilC; k++)
                        {
                          det[ilC - 1 - k + ilA + ilB] = C.l[ilC].o[iNC * ilC + k] ^ CrAnMask;
                        }
                        if (ComplexOperatorAddLaderNormalOrder
                            (O, det, ilA + ilB + ilC, A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Re * C.l[ilC].VRe[iNC],
                             A.l[ilA].VRe[iNA] * B.l[ilB].V[iNB].Im * C.l[ilC].VRe[iNC], QDetNotOrdered))
                        {
                          printf("ComplexOperatorAddLader failed in OperatorMultiply3\n");
                          fflush(stdout);
                          return 1;
                        }
                      }
                    }
                  }
                }
              }
              else
              {                 /* A and B real */
                if (C.l[ilC].QComplex)
                {               /* A and B real C complex */
                  for (iNA = 0; iNA < A.l[ilA].N; iNA++)
                  {
                    for (k = 0; k < ilA; k++)
                    {
                      det[k] = A.l[ilA].o[iNA * ilA + k];
                    }
                    for (iNB = 0; iNB < B.l[ilB].N; iNB++)
                    {
                      for (k = 0; k < ilB; k++)
                      {
                        det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                      }
                      for (iNC = 0; iNC < C.l[ilC].N; iNC++)
                      {
                        for (k = 0; k < ilC; k++)
                        {
                          det[ilC - 1 - k + ilA + ilB] = C.l[ilC].o[iNC * ilC + k] ^ CrAnMask;
                        }
                        if (ComplexOperatorAddLaderNormalOrder
                            (O, det, ilA + ilB + ilC, A.l[ilA].VRe[iNA] * B.l[ilB].VRe[iNB] * C.l[ilC].V[iNC].Re,
                             -A.l[ilA].VRe[iNA] * B.l[ilB].VRe[iNB] * C.l[ilC].V[iNC].Im, QDetNotOrdered))
                        {
                          printf("ComplexOperatorAddLader failed in OperatorMultiply3\n");
                          fflush(stdout);
                          return 1;
                        }
                      }
                    }
                  }
                }
                else
                {               /* A, B and C real (O complex) */
                  for (iNA = 0; iNA < A.l[ilA].N; iNA++)
                  {
                    for (k = 0; k < ilA; k++)
                    {
                      det[k] = A.l[ilA].o[iNA * ilA + k];
                    }
                    for (iNB = 0; iNB < B.l[ilB].N; iNB++)
                    {
                      for (k = 0; k < ilB; k++)
                      {
                        det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                      }
                      for (iNC = 0; iNC < C.l[ilC].N; iNC++)
                      {
                        for (k = 0; k < ilC; k++)
                        {
                          det[ilC - 1 - k + ilA + ilB] = C.l[ilC].o[iNC * ilC + k] ^ CrAnMask;
                        }
                        if (ComplexOperatorAddLaderNormalOrder
                            (O, det, ilA + ilB + ilC, A.l[ilA].VRe[iNA] * B.l[ilB].VRe[iNB] * C.l[ilC].VRe[iNC], 0.0,
                             QDetNotOrdered))
                        {
                          printf("ComplexOperatorAddLader failed in OperatorMultiply3\n");
                          fflush(stdout);
                          return 1;
                        }
                      }
                    }
                  }
                }
              }
            }
          }
          else
          {                     /* All real */
            for (iNA = 0; iNA < A.l[ilA].N; iNA++)
            {
              for (k = 0; k < ilA; k++)
              {
                det[k] = A.l[ilA].o[iNA * ilA + k];
              }
              for (iNB = 0; iNB < B.l[ilB].N; iNB++)
              {
                for (k = 0; k < ilB; k++)
                {
                  det[k + ilA] = B.l[ilB].o[iNB * ilB + k];
                }
                for (iNC = 0; iNC < C.l[ilC].N; iNC++)
                {
                  for (k = 0; k < ilC; k++)
                  {
                    det[ilC - 1 - k + ilA + ilB] = C.l[ilC].o[iNC * ilC + k] ^ CrAnMask;
                  }
                  if (RealOperatorAddLaderNormalOrder
                      (O, det, ilA + ilB + ilC, A.l[ilA].VRe[iNA] * B.l[ilB].VRe[iNB] * C.l[ilC].VRe[iNC],
                       QDetNotOrdered))
                  {
                    printf("RealOperatorAddLader failed in OperatorMultiply3\n");
                    fflush(stdout);
                    return 1;
                  }
                }
              }
            }
          }
        }
      }
    }
    /* End of NND */
  }
  else if ((Aset == 'N') && (Bset == 'D') && (Cset == 'D'))
  {
    printf("in OperatorMultiply3 Still to do\n");
    return 1;
  }
  else if ((Aset == 'D') && (Bset == 'N') && (Cset == 'D'))
  {
    printf("in OperatorMultiply3 Still to do\n");
    return 1;
  }
  else if ((Aset == 'D') && (Bset == 'D') && (Cset == 'D'))
  {
    printf("in OperatorMultiply3 Still to do\n");
    return 1;
  }
  else
  {
    printf("in OperatorMultiply3 unrecognized input for Aset and or Bset\n");
    return 1;
  }
  return 0;
}



/*  Depending on Aset, Bset, Cset and Dset
 *    O <- O + D^'Dset' C^'Cset' B^'Bset' A^'Aset'
 *    whereby set can be
 *    'N' normal
 *    'C' conjugated (a^dag remains a^dag, prefactors become complex conjugated)
 *    'T' transposed (a^dag becomes a, prefactors remain)
 *    'D' conjugate transposed (a^dag becomes a, prefactors become complex conjugated)
 */
uint32_t OperatorMultiply4(OperatorType * const O, OperatorType const D, OperatorType const C, OperatorType const B,
                           OperatorType const A, char const Dset, char const Cset, char const Bset, char const Aset)
{
  OperatorType BA, DC;
  InitOperatorForMult(&BA, 2, &B, &A);
  InitOperatorForMult(&DC, 2, &D, &C);

  if (OperatorMultiply2(&BA, B, A, Bset, Aset))
  {
    printf("OperatorMultiply2 failed in OperatorMultiply4\n");
    return 1;
  };
  if (OperatorMultiply2(&DC, D, C, Dset, Cset))
  {
    printf("OperatorMultiply2 failed in OperatorMultiply4\n");
    return 1;
  };
  if (OperatorMultiply2(O, DC, BA, 'N', 'N'))
  {
    printf("OperatorMultiply2 failed in OperatorMultiply4\n");
    return 1;
  };
  FreeOperator(&BA);
  FreeOperator(&DC);
  return 0;
}


/*  Depending on Aset, Bset, Cset and Dset
 *    O <- O + D^'Dset' C^'Cset' B^'Bset' A^'Aset'
 *    whereby set can be
 *    'N' normal
 *    'C' conjugated (a^dag remains a^dag, prefactors become complex conjugated)
 *    'T' transposed (a^dag becomes a, prefactors remain)
 *    'D' conjugate transposed (a^dag becomes a, prefactors become complex conjugated)
 */
uint32_t OperatorMultiply4NormalOrder(OperatorType * const O, OperatorType const D, OperatorType const C,
                                      OperatorType const B, OperatorType const A, char const Dset, char const Cset,
                                      char const Bset, char const Aset, uint32_t(*QDetNotOrdered) (uint16_t, uint16_t, uint32_t))
{
  OperatorType BA, DC;
  InitOperatorForMult(&BA, 2, &B, &A);
  InitOperatorForMult(&DC, 2, &D, &C);

  if (OperatorMultiply2NormalOrder(&BA, B, A, Bset, Aset, QDetNotOrdered))
  {
    printf("OperatorMultiply2 failed in OperatorMultiply4\n");
    return 1;
  };
  if (OperatorMultiply2NormalOrder(&DC, D, C, Dset, Cset, QDetNotOrdered))
  {
    printf("OperatorMultiply2 failed in OperatorMultiply4\n");
    return 1;
  };
  if (OperatorMultiply2NormalOrder(O, DC, BA, 'N', 'N', QDetNotOrdered))
  {
    printf("OperatorMultiply2 failed in OperatorMultiply4\n");
    return 1;
  };
  FreeOperator(&BA);
  FreeOperator(&DC);
  return 0;
}

/*  Depending on Aset, Bset, Cset, Dset and Eset
 *    O <- O + E^'Eset' D^'Dset' C^'Cset' B^'Bset' A^'Aset'
 *    whereby set can be
 *    'N' normal
 *    'C' conjugated (a^dag remains a^dag, prefactors become complex conjugated)
 *    'T' transposed (a^dag becomes a, prefactors remain)
 *    'D' conjugate transposed (a^dag becomes a, prefactors become complex conjugated)
 */
uint32_t OperatorMultiply5(OperatorType * const O, OperatorType const E, OperatorType const D, OperatorType const C,
                           OperatorType const B, OperatorType const A, char const Eset, char const Dset,
                           char const Cset, char const Bset, char const Aset)
{
  OperatorType ED, BA;
  InitOperatorForMult(&BA, 2, &B, &A);
  InitOperatorForMult(&ED, 2, &E, &D);

  if (OperatorMultiply2(&BA, B, A, Bset, Aset))
  {
    printf("OperatorMultiply2 failed in OperatorMultiply5\n");
    return 1;
  };
  if (OperatorMultiply2(&ED, E, D, Eset, Dset))
  {
    printf("OperatorMultiply2 failed in OperatorMultiply5\n");
    return 1;
  };
  if (OperatorMultiply3(O, ED, C, BA, 'N', Cset, 'N'))
  {
    printf("OperatorMultiply3 failed in OperatorMultiply5\n");
    return 1;
  };
  FreeOperator(&ED);
  FreeOperator(&BA);
  return 0;
}


/*  Depending on Aset, Bset, Cset, Dset and Eset
 *    O <- O + E^'Eset' D^'Dset' C^'Cset' B^'Bset' A^'Aset'
 *    whereby set can be
 *    'N' normal
 *    'C' conjugated (a^dag remains a^dag, prefactors become complex conjugated)
 *    'T' transposed (a^dag becomes a, prefactors remain)
 *    'D' conjugate transposed (a^dag becomes a, prefactors become complex conjugated)
 */
uint32_t OperatorMultiply5NormalOrder(OperatorType * const O, OperatorType const E, OperatorType const D,
                                      OperatorType const C, OperatorType const B, OperatorType const A, char const Eset,
                                      char const Dset, char const Cset, char const Bset, char const Aset,
                                      uint32_t(*QDetNotOrdered) (uint16_t, uint16_t, uint32_t))
{
  OperatorType ED, BA;
  InitOperatorForMult(&BA, 2, &B, &A);
  InitOperatorForMult(&ED, 2, &E, &D);

  if (OperatorMultiply2NormalOrder(&BA, B, A, Bset, Aset, QDetNotOrdered))
  {
    printf("OperatorMultiply2 failed in OperatorMultiply5\n");
    return 1;
  };
  if (OperatorMultiply2NormalOrder(&ED, E, D, Eset, Dset, QDetNotOrdered))
  {
    printf("OperatorMultiply2 failed in OperatorMultiply5\n");
    return 1;
  };
  if (OperatorMultiply3NormalOrder(O, ED, C, BA, 'N', Cset, 'N', QDetNotOrdered))
  {
    printf("OperatorMultiply3 failed in OperatorMultiply5\n");
    return 1;
  };
  FreeOperator(&BA);
  FreeOperator(&ED);
  return 0;
}

/*  Calculate e^(alpha O)
 *  only works if O has elements for n=2 and/or n=0
 *  whereby the elements of n=2 need to be of the form Cr_tau An_tau'
 *  and O needs to be hermitian 
 *  alpha can be a complex number
 *
 *  used is that e^(alpha O) can be calcualted as:
 *  O = \sum_{tau,tau'} alpha_{tau,tau'} a^dag_tau a_tau'
 *    = \sum_tau beta_tau a^dag_tau a_tau
 *  e^{alpha O} = e^{alpha} e^{beta_0} \Pi_{tau} (1 + (e^{beta_tau}-1) b^dag_tau b_tau)
 *
 */
uint32_t HermitianOperatorExponent(OperatorType const O, Complex16 alpha, OperatorType * const ExpO)
{
  uint16_t oadd[2];
  uint32_t i,j,k;
  CompactMatrixType Omat,fun;
  Complex16 c={0,0};
  double val[O.NFermion],prefac;
  OperatorType Ob;
  /* test if operator fullfills basic conditions */
  if(O.NBoson!=0)
  {
    printf("OperatorExponent failed, tried to take the exponent of an operator with bosons Nboson=%i\n",O.NBoson);
    return 1;
  }
  if(O.MaxLength>2)
  {
    printf("OperatorExponent failed, tried to take the exponent of an operator with max Length %i\n",O.MaxLength);
    return 1;
  }
  if(O.MaxLength<2)
  {
    if(O.MaxLength==1)
    {
      if(O.l[1].N!=0)
      {
        printf("OperatorExponent failed, tried to take the exponent of an operator with single creation or annihilation operators MaxLength=%i and there are %i elements of length 1\n",O.MaxLength,O.l[1].N);
        return 1;
      }
    }
    InitOperatorCopySettings(ExpO, O);
    if(O.l[0].N!=0)
    {
      if(O.l[0].QComplex)
      {
        ComplexOperatorAddLader(ExpO, oadd, 0, exp(O.l[0].V[0].Re)*cos(O.l[0].V[0].Im), exp(O.l[0].V[0].Re)*sin(O.l[0].V[0].Im));
      }
      else
      {
        RealOperatorAddLader(ExpO, oadd, 0, exp(O.l[0].VRe[0]));
      }
    }
    return 0;
  }
  if(O.l[1].N!=0)
  {
    printf("OperatorExponent failed, tried to take the exponent of an operator with single creation or annihilation operators MaxLength=%i and there are %i elements of length 1\n",O.MaxLength,O.l[1].N);
    return 1;
  }
  for(i=0;i<O.l[2].N;i++)
  {
    if(!((O.l[2].o[2*i]^O.l[2].o[2*i+1])&CrAnMask))
    {
      printf("OperatorExponent failed, operator of length 2 at index %i consists of two creation or two annihilation operators\n",i);
      return 1;
    }
  }
  /* create matrix representing O */
  Omat.Ni=O.NFermion;
  Omat.Nj=O.NFermion;
  Omat.QComplex=O.l[2].QComplex;
  InitCompactMatrix(&Omat);
  if(Omat.QComplex)
  {
    for(i=0;i<O.l[2].N;i++)
    {
      if(O.l[2].o[2*i]&CrAnMask)
      {
        Omat.V[(uint32_t)((O.l[2].o[2*i]&NotCrAnMask)*Omat.Nj+(O.l[2].o[2*i+1]&NotCrAnMask))].Re-=O.l[2].V[i].Re;
        Omat.V[(uint32_t)((O.l[2].o[2*i]&NotCrAnMask)*Omat.Nj+(O.l[2].o[2*i+1]&NotCrAnMask))].Im-=O.l[2].V[i].Im;
        if((O.l[2].o[2*i]&NotCrAnMask)==(O.l[2].o[2*i+1]&NotCrAnMask))
        {
          c.Re+=O.l[2].V[i].Re;
          c.Im+=O.l[2].V[i].Im;
        }
      }
      else
      {
        Omat.V[(uint32_t)((O.l[2].o[2*i+1]&NotCrAnMask)*Omat.Nj+(O.l[2].o[2*i]&NotCrAnMask))].Re+=O.l[2].V[i].Re;
        Omat.V[(uint32_t)((O.l[2].o[2*i+1]&NotCrAnMask)*Omat.Nj+(O.l[2].o[2*i]&NotCrAnMask))].Im+=O.l[2].V[i].Im;
      }
    }
  }
  else
  {
    for(i=0;i<O.l[2].N;i++)
    {
      if(O.l[2].o[2*i]&CrAnMask)
      {
        Omat.VRe[(uint32_t)((O.l[2].o[2*i]&NotCrAnMask)*Omat.Nj+(O.l[2].o[2*i+1]&NotCrAnMask))]-=O.l[2].VRe[i];
        if((O.l[2].o[2*i]&NotCrAnMask)==(O.l[2].o[2*i+1]&NotCrAnMask))
        {
          c.Re+=O.l[2].VRe[i];
        }
      }
      else
      {
        Omat.VRe[(uint32_t)((O.l[2].o[2*i+1]&NotCrAnMask)*Omat.Nj+(O.l[2].o[2*i]&NotCrAnMask))]+=O.l[2].VRe[i];
      }
    }
  }
  /* test if Omat is Hermitian */
  if(QCompactMatrixHermitian(Omat, DBL_EPSILON)==0)
  {
    printf("OperatorExponent failed, operator of length 2 should be Hermitian\n");
    PrintCompactMatrixMinimal(Omat);
    return 1;
  }
  /* create eigenbasis */
  CompactMatrixEigensystem(&Omat, &fun, val);
  
  //PrintCompactMatrixMinimal(Omat);
  //printf("%22.15E %22.15E %22.15E %22.15E %22.15E %22.15E \n\n",val[0],val[1],val[2],val[3],val[4],val[5]);
  //PrintCompactMatrixMinimal(fun);
  
  /* init operator ExpO to be the identity operator */
  InitOperator(ExpO);
  ExpO->NFermion=O.NFermion;
  ExpO->NBoson=O.NBoson;
  OperatorAddLader(ExpO, oadd, 0, 1, 0);
  /* run over all eigen-states and create ExpO */
  for(i=0;i<O.NFermion;i++)
  {
    /* init operator Ob to be the identity operator */
    InitOperator(&Ob);
    Ob.NFermion=O.NFermion;
    Ob.NBoson=O.NBoson;
    InitOperatorFixedLength(&Ob, 2, 1);
    OperatorAddLader(&Ob, oadd, 0, 1, 0);
    /* add (e^beta_tau -1) b^dag_tau b^_tau to Ob */
    prefac=exp(val[i])-1.0;
    if(fun.QComplex)
    {
      for(j=0;j<O.NFermion;j++)
      {
        for(k=0;k<O.NFermion;k++)
        {
          oadd[1]=(uint16_t)j+CrAnMask;
          oadd[0]=(uint16_t)k;
          ComplexOperatorAddLader(&Ob, oadd, 2, prefac*(fun.V[i*fun.Nj+j].Re*fun.V[i*fun.Nj+k].Re+fun.V[i*fun.Nj+j].Im*fun.V[i*fun.Nj+k].Im), prefac*(fun.V[i*fun.Nj+j].Im*fun.V[i*fun.Nj+k].Re-fun.V[i*fun.Nj+j].Re*fun.V[i*fun.Nj+k].Im));
        }
      }
    }
    else
    {
      for(j=0;j<O.NFermion;j++)
      {
        for(k=0;k<O.NFermion;k++)
        {
          oadd[1]=(uint16_t)j+CrAnMask;
          oadd[0]=(uint16_t)k;
          RealOperatorAddLader(&Ob, oadd, 2, prefac*fun.VRe[i*fun.Nj+j]*fun.VRe[i*fun.Nj+k]);
        }
      }
    }
    OperatorChop(&Ob, 100*DBL_EPSILON);
    /* multiply Ob into ExpO */
    OperatorMultiplyInline(ExpO, Ob);
    //OperatorMultiplyNormalOrderInline(ExpO, Ob, &QDetNotOrderedFcFaBcBa);
    OperatorChop(ExpO, 100*DBL_EPSILON);
    /* Free operator Ob */
    FreeOperator(&Ob);
  }
  OperatorScale(ExpO, exp(c.Re)*cos(c.Im), exp(c.Re)*sin(c.Im));
  for(i=0;i<O.l[0].N;i++)
  {
    if(O.l[0].QComplex)
    {
      OperatorScale(ExpO, exp(O.l[0].V[i].Re)*cos(O.l[0].V[i].Im), exp(O.l[0].V[i].Re)*sin(O.l[0].V[i].Im));
    }
    else
    {
      OperatorScale(ExpO, exp(O.l[0].VRe[i]), 0);
    }
  }
  return 0;
}
uint32_t HermitianOperatorNormalOrderExponent(OperatorType const O, Complex16 alpha, OperatorType * const ExpO, uint32_t(*QDetNotOrdered) (uint16_t, uint16_t, uint32_t))
{
  uint16_t oadd[2];
  uint32_t i,j,k;
  CompactMatrixType Omat,fun;
  Complex16 c={0,0};
  Complex16 prefac;
  double val[O.NFermion];
  OperatorType Ob;
  /* test if operator fullfills basic conditions */
  if(O.NBoson!=0)
  {
    printf("OperatorExponent failed, tried to take the exponent of an operator with bosons Nboson=%i\n",O.NBoson);
    return 1;
  }
  if(O.MaxLength>2)
  {
    printf("OperatorExponent failed, tried to take the exponent of an operator with max Length %i\n",O.MaxLength);
    return 1;
  }
  if(O.MaxLength<2)
  {
    if(O.MaxLength==1)
    {
      if(O.l[1].N!=0)
      {
        printf("OperatorExponent failed, tried to take the exponent of an operator with single creation or annihilation operators MaxLength=%i and there are %i elements of length 1\n",O.MaxLength,O.l[1].N);
        return 1;
      }
    }
    InitOperatorCopySettings(ExpO, O);
    if(O.l[0].N!=0)
    {
      if(O.l[0].QComplex)
      {
        ComplexOperatorAddLaderNormalOrder(ExpO, oadd, 0, exp(O.l[0].V[0].Re)*cos(O.l[0].V[0].Im), exp(O.l[0].V[0].Re)*sin(O.l[0].V[0].Im),QDetNotOrdered);
      }
      else
      {
        RealOperatorAddLaderNormalOrder(ExpO, oadd, 0, exp(O.l[0].VRe[0]),QDetNotOrdered);
      }
    }
    return 0;
  }
  if(O.l[1].N!=0)
  {
    printf("OperatorExponent failed, tried to take the exponent of an operator with single creation or annihilation operators MaxLength=%i and there are %i elements of length 1\n",O.MaxLength,O.l[1].N);
    return 1;
  }
  for(i=0;i<O.l[2].N;i++)
  {
    if(!((O.l[2].o[2*i]^O.l[2].o[2*i+1])&CrAnMask))
    {
      printf("OperatorExponent failed, operator of length 2 at index %i consists of two creation or two annihilation operators\n",i);
      return 1;
    }
  }
  /* create matrix representing O */
  Omat.Ni=O.NFermion;
  Omat.Nj=O.NFermion;
  Omat.QComplex=O.l[2].QComplex;
  InitCompactMatrix(&Omat);
  if(Omat.QComplex)
  {
    for(i=0;i<O.l[2].N;i++)
    {
      if(O.l[2].o[2*i]&CrAnMask)
      {
        Omat.V[(uint32_t)((O.l[2].o[2*i]&NotCrAnMask)*Omat.Nj+(O.l[2].o[2*i+1]&NotCrAnMask))].Re-=O.l[2].V[i].Re;
        Omat.V[(uint32_t)((O.l[2].o[2*i]&NotCrAnMask)*Omat.Nj+(O.l[2].o[2*i+1]&NotCrAnMask))].Im-=O.l[2].V[i].Im;
        if((O.l[2].o[2*i]&NotCrAnMask)==(O.l[2].o[2*i+1]&NotCrAnMask))
        {
          c.Re+=O.l[2].V[i].Re;
          c.Im+=O.l[2].V[i].Im;
        }
      }
      else
      {
        Omat.V[(uint32_t)((O.l[2].o[2*i+1]&NotCrAnMask)*Omat.Nj+(O.l[2].o[2*i]&NotCrAnMask))].Re+=O.l[2].V[i].Re;
        Omat.V[(uint32_t)((O.l[2].o[2*i+1]&NotCrAnMask)*Omat.Nj+(O.l[2].o[2*i]&NotCrAnMask))].Im+=O.l[2].V[i].Im;
      }
    }
  }
  else
  {
    for(i=0;i<O.l[2].N;i++)
    {
      if(O.l[2].o[2*i]&CrAnMask)
      {
        Omat.VRe[(uint32_t)((O.l[2].o[2*i]&NotCrAnMask)*Omat.Nj+(O.l[2].o[2*i+1]&NotCrAnMask))]-=O.l[2].VRe[i];
        if((O.l[2].o[2*i]&NotCrAnMask)==(O.l[2].o[2*i+1]&NotCrAnMask))
        {
          c.Re+=O.l[2].VRe[i];
        }
      }
      else
      {
        Omat.VRe[(uint32_t)((O.l[2].o[2*i+1]&NotCrAnMask)*Omat.Nj+(O.l[2].o[2*i]&NotCrAnMask))]+=O.l[2].VRe[i];
      }
    }
  }
  /* test if Omat is Hermitian */
  if(QCompactMatrixHermitian(Omat, DBL_EPSILON)==0)
  {
    printf("OperatorExponent failed, operator of length 2 should be Hermitian\n");
    PrintCompactMatrixMinimal(Omat);
    return 1;
  }
  /* create eigenbasis */
  CompactMatrixEigensystem(&Omat, &fun, val);
  /* init operator ExpO to be the identity operator */
  InitOperator(ExpO);
  ExpO->NFermion=O.NFermion;
  ExpO->NBoson=O.NBoson;
  OperatorAddLaderNormalOrder(ExpO, oadd, 0, 1, 0, QDetNotOrdered);
  /* run over all eigen-states and create ExpO */
  for(i=0;i<O.NFermion;i++)
  {
    /* init operator Ob to be the identity operator */
    InitOperator(&Ob);
    Ob.NFermion=O.NFermion;
    Ob.NBoson=O.NBoson;
    InitOperatorFixedLength(&Ob, 2, 1);
    OperatorAddLaderNormalOrder(&Ob, oadd, 0, 1, 0, QDetNotOrdered);
    /* add (e^(alpha*beta_tau) -1) b^dag_tau b^_tau to Ob */
    prefac.Re=exp(val[i]*alpha.Re)*cos(val[i]*alpha.Im)-1.0;
    prefac.Im=exp(val[i]*alpha.Re)*sin(val[i]*alpha.Im);
    if(fun.QComplex)
    {
      for(j=0;j<O.NFermion;j++)
      {
        for(k=0;k<O.NFermion;k++)
        {
          oadd[1]=(uint16_t)j+CrAnMask;
          oadd[0]=(uint16_t)k;
          ComplexOperatorAddLaderNormalOrder(&Ob, oadd, 2, prefac.Re*(fun.V[i*fun.Nj+j].Re*fun.V[i*fun.Nj+k].Re+fun.V[i*fun.Nj+j].Im*fun.V[i*fun.Nj+k].Im)
                                                          -prefac.Im*(fun.V[i*fun.Nj+j].Im*fun.V[i*fun.Nj+k].Re-fun.V[i*fun.Nj+j].Re*fun.V[i*fun.Nj+k].Im),
                                                           prefac.Re*(fun.V[i*fun.Nj+j].Im*fun.V[i*fun.Nj+k].Re-fun.V[i*fun.Nj+j].Re*fun.V[i*fun.Nj+k].Im)
                                                          +prefac.Im*(fun.V[i*fun.Nj+j].Re*fun.V[i*fun.Nj+k].Re+fun.V[i*fun.Nj+j].Im*fun.V[i*fun.Nj+k].Im),QDetNotOrdered);
        }
      }
    }
    else if(fabs(alpha.Im)>100*DBL_MIN)
    {
      for(j=0;j<O.NFermion;j++)
      {
        for(k=0;k<O.NFermion;k++)
        {
          oadd[1]=(uint16_t)j+CrAnMask;
          oadd[0]=(uint16_t)k;
          ComplexOperatorAddLaderNormalOrder(&Ob, oadd, 2, prefac.Re*(fun.VRe[i*fun.Nj+j]*fun.VRe[i*fun.Nj+k]),prefac.Im*(fun.VRe[i*fun.Nj+j]*fun.VRe[i*fun.Nj+k]),QDetNotOrdered);
        }
      }
    }
    else
    {
      for(j=0;j<O.NFermion;j++)
      {
        for(k=0;k<O.NFermion;k++)
        {
          oadd[1]=(uint16_t)j+CrAnMask;
          oadd[0]=(uint16_t)k;
          RealOperatorAddLaderNormalOrder(&Ob, oadd, 2, prefac.Re*fun.VRe[i*fun.Nj+j]*fun.VRe[i*fun.Nj+k],QDetNotOrdered);
        }
      }
    }
    OperatorChop(&Ob, 100*DBL_EPSILON);
    /* multiply Ob into ExpO */
    OperatorMultiplyNormalOrderInline(ExpO, Ob, QDetNotOrdered);
    OperatorChop(ExpO, 100*DBL_EPSILON);
    /* Free operator Ob */
    FreeOperator(&Ob);
  }
  OperatorScale(ExpO, exp(c.Re)*cos(c.Im), exp(c.Re)*sin(c.Im));
  for(i=0;i<O.l[0].N;i++)
  {
    if(O.l[0].QComplex)
    {
      OperatorScale(ExpO, exp(O.l[0].V[i].Re)*cos(O.l[0].V[i].Im), exp(O.l[0].V[i].Re)*sin(O.l[0].V[i].Im));
    }
    else
    {
      OperatorScale(ExpO, exp(O.l[0].VRe[i]), 0);
    }
  }
  return 0;
}



/* Multiply an operator by an scalar
*  O <- alpha O
*/
uint32_t OperatorScale(OperatorType * const O, const double alphaRe, const double alphaIm)
{
  uint32_t i, j;
  double Vtemp;
  for (i = 0; i <= O->MaxLength; i++)
  {
    if (O->l[i].QComplex)
    {                           /* Operator Complex at this length */
      for (j = 0; j < O->l[i].N; j++)
      {
        Vtemp = O->l[i].V[j].Re * alphaRe - O->l[i].V[j].Im * alphaIm;
        O->l[i].V[j].Im = O->l[i].V[j].Im * alphaRe + O->l[i].V[j].Re * alphaIm;
        O->l[i].V[j].Re = Vtemp;
      }
    }
    else
    {                           /* Operator Real at this length */
      if (fabs(alphaIm) < (DBL_MIN * 10))
      {                         /* multiplication also real */
        for (j = 0; j < O->l[i].N; j++)
        {
          O->l[i].VRe[j] *= alphaRe;
        }
      }
      else
      {                         /* O real, but multiplied with a complex number */
        if (MakeRealOperatorFixedLengthComplex(O, i))
        {
          printf("MakeRealOperatorFixedLengthComplex failed in OperatorScale\n");
          fflush(stdout);
          return 1;
        }
        for (j = 0; j < O->l[i].N; j++)
        {
          O->l[i].V[j].Im = O->l[i].V[j].Re * alphaIm;
          O->l[i].V[j].Re = O->l[i].V[j].Re * alphaRe;
        }
      }
    }
  }
  SetQComplexOperator(O);
  return 0;
}



/*  add the scaled normal ordered product of two operators to AB
 *  O <- O + alpha B A
 *  O must be initiallized, i.e O.Maxlength must be A.Maxlength+B.Maxlength
 *  it is advised to set O.l[il].NBitsKey to sensible values
 */
uint32_t OperatorMultiplyScaledNormalOrder(OperatorType * const O, OperatorType const B, OperatorType const A,
                                           const double alphaRe, const double alphaIm,
                                           uint32_t(*QDetNotOrdered) (uint16_t, uint16_t, uint32_t))
{
  if(OperatorMultiply2NormalOrder(O, B, A, 'N', 'N', QDetNotOrdered))
  {
    printf("OperatorMultiply2NormalOrder failed in OperatorMultiplyScaledNormalOrder\n");
    return 1;
  };
  OperatorScale(O, alphaRe, alphaIm);
  return 0;
}




/*  OA <- OA + alpha OB
 *  if OB.Maxlength > OA.Maxlength then OA is extended
 *  the values of NBitsKey are copied
 */
uint32_t OperatorAdd(OperatorType * const OA, const OperatorType OB, const double alphaRe, const double alphaIm)
{
  uint32_t i, j;
  for (i = OA->MaxLength + 1; i <= OB.MaxLength; i++)
  {
    if (OB.l[i].NBitsKey)
      if (InitOperatorFixedLength(OA, i, OB.l[i].NBitsKey))
      {
        printf("InitOperatorFixedLength failed in OperatorAdd");
        fflush(stdout);
        return 1;
      }
  }
  for (i = 0; i <= OB.MaxLength; i++)
  {
    if ((OA->l[i].QComplex) || (OB.l[i].QComplex) || (fabs(alphaIm) > DBL_MIN * 10))
    {
      if (!(OA->l[i].QComplex))
        if (MakeRealOperatorFixedLengthComplex(OA, i))
          return 1;
      if (OB.l[i].QComplex)
      {                         /* all complex */
        for (j = 0; j < OB.l[i].N; j++)
        {
          if (ComplexOperatorAddLader
              (OA, &(OB.l[i].o[i * j]), i, alphaRe * OB.l[i].V[j].Re - alphaIm * OB.l[i].V[j].Im,
               alphaIm * OB.l[i].V[j].Re + alphaRe * OB.l[i].V[j].Im))
          {
            printf("ComplexOperatorAddLader failed in OperatorAdd\n");
            fflush(stdout);
            return 1;
          }
        }
      }
      else
      {                         /* OB real, but otherwise all complex */
        for (j = 0; j < OB.l[i].N; j++)
        {
          if (ComplexOperatorAddLader(OA, &(OB.l[i].o[i * j]), i, alphaRe * OB.l[i].VRe[j], alphaIm * OB.l[i].VRe[j]))
          {
            printf("ComplexOperatorAddLader failed in OperatorAdd\n");
            fflush(stdout);
            return 1;
          }
        }
      }
    }
    else
    {                           /* all real */
      for (j = 0; j < OB.l[i].N; j++)
      {
        if (RealOperatorAddLader(OA, &(OB.l[i].o[i * j]), i, alphaRe * OB.l[i].VRe[j]))
        {
          printf("RealOperatorAddLader failed in OperatorAdd\n");
          fflush(stdout);
          return 1;
        }
      }
    }
  }
  SetQComplexOperator(OA);
  return 0;
}

uint32_t OperatorAdd(OperatorType * const OA, const OperatorType OB)
{
	return OperatorAdd(OA, OB, 1.0, 0.0);
}


/*  OA <- OA + alpha NormalOrder[OB]
 *  if OB.Maxlength > OA.Maxlength then OA is extended
 *  the values of NBitsKey are copied
 */
uint32_t OperatorAddNormalOrder(OperatorType * const OA, const OperatorType OB, const double alphaRe, const double alphaIm, uint32_t(*QDetNotOrdered) (uint16_t, uint16_t, uint32_t))
{
  uint32_t i, j;
  for (i = OA->MaxLength + 1; i <= OB.MaxLength; i++)
  {
    if (OB.l[i].NBitsKey)
      if (InitOperatorFixedLength(OA, i, OB.l[i].NBitsKey))
      {
        printf("InitOperatorFixedLength failed in OperatorAdd");
        fflush(stdout);
        return 1;
      }
  }
  for (i = 0; i <= OB.MaxLength; i++)
  {
    if ((OA->l[i].QComplex) || (OB.l[i].QComplex) || (fabs(alphaIm) > DBL_MIN * 10))
    {
      if (!(OA->l[i].QComplex))
        if (MakeRealOperatorFixedLengthComplex(OA, i))
          return 1;
      if (OB.l[i].QComplex)
      {                         /* all complex */
        for (j = 0; j < OB.l[i].N; j++)
        {
          if (ComplexOperatorAddLaderNormalOrder
              (OA, &(OB.l[i].o[i * j]), i, alphaRe * OB.l[i].V[j].Re - alphaIm * OB.l[i].V[j].Im,
               alphaIm * OB.l[i].V[j].Re + alphaRe * OB.l[i].V[j].Im, QDetNotOrdered))
          {
            printf("ComplexOperatorAddLader failed in OperatorAdd\n");
            fflush(stdout);
            return 1;
          }
        }
      }
      else
      {                         /* OB real, but otherwise all complex */
        for (j = 0; j < OB.l[i].N; j++)
        {
          if (ComplexOperatorAddLaderNormalOrder
              (OA, &(OB.l[i].o[i * j]), i, alphaRe * OB.l[i].VRe[j], alphaIm * OB.l[i].VRe[j], QDetNotOrdered))
          {
            printf("ComplexOperatorAddLader failed in OperatorAdd\n");
            fflush(stdout);
            return 1;
          }
        }
      }
    }
    else
    {                           /* all real */
      for (j = 0; j < OB.l[i].N; j++)
      {
        if (RealOperatorAddLaderNormalOrder(OA, &(OB.l[i].o[i * j]), i, alphaRe * OB.l[i].VRe[j], QDetNotOrdered))
        {
          printf("RealOperatorAddLader failed in OperatorAdd\n");
          fflush(stdout);
          return 1;
        }
      }
    }
  }
  SetQComplexOperator(OA);
  return 0;
}




/* OA <- OA + U OB
 *  U is a (NFermionB+NBosonB)*(NFermionA+NBosonA) demensional matrix that defines the transformation
 *  U.Ni = (NFermionA+NBosonA)
 *  U.Nj = (NFermionB+NBosonB)
 *  a_i        = \sum_j U_{i,j}   a_j
 *  a_i^{\dag} = \sum_j U_{i,j}^* a_j^{dag}
 *
 * an operator is given as sum_j (alpha_j pi_i a_i)
 * this function rotates the creation/annihilation operator at i=level with the knowledge that i runs from 0 to lenght-1
 */
uint32_t MappedOperatorFixedLengthAtLevel(OperatorFixedLengthType * const OA, const OperatorFixedLengthType OB, const MatrixType U, const uint32_t level, const uint32_t length)
{
  uint32_t j, i, k;
  OA->N = 0;
  OA->NAloc = 16;
  if (OB.QComplex)
  {
    if (U.QComplex)
    {                           /* OB and U complex */
      OA->QComplex = 1;
      if ((OA->V = (Complex16*)malloc(OA->NAloc * sizeof(Complex16))) == NULL)
      {
        printf("Could not allocate memory for OA->V in MappedOperatorFixedLengthAtLevel\n");
        return 1;
      }
      if ((OA->o = (uint16_t*)malloc(OA->NAloc * length * sizeof(uint16_t))) == NULL)
      {
        printf("Could not allocate memory for OA->o in MappedOperatorFixedLengthAtLevel\n");
        return 1;
      }
      for (j= 0; j< OB.N; j++)
      {
        for (i = 0; i < U.Ni; i++)
        {
          if (((fabs(U.VRe[i][(OB.o[j * length + level] & NotCrAnMask)]) +
                fabs(U.VIm[i][(OB.o[j * length + level] & NotCrAnMask)])) * (fabs(OB.V[j].Re) + fabs(OB.V[j].Im))) >
              (DBL_MIN * 10))
          {
            if (OA->N == OA->NAloc)
            {
              OA->NAloc *= 2;
              if ((OA->V = (Complex16*)realloc(OA->V, OA->NAloc * sizeof(Complex16))) == NULL)
              {
                printf("Could not reallocate memory for OA->V in MappedOperatorFixedLengthAtLevel\n");
                return 1;
              }
              if ((OA->o = (uint16_t*)realloc(OA->o, OA->NAloc * length * sizeof(uint16_t))) == NULL)
              {
                printf("Could not reallocate memory for OA->o in MappedOperatorFixedLengthAtLevel\n");
                return 1;
              }
            }
            if (OB.o[j* length + level] & CrAnMask)
            {                   /* transformation of a creation operator a_i^{\dag} = \sum_j U_{i,j}^* a_j^{dag} */
              OA->V[OA->N].Re =
                U.VRe[i][(OB.o[j* length + level] & NotCrAnMask)] * OB.V[j].Re +
                U.VIm[i][(OB.o[j* length + level] & NotCrAnMask)] * OB.V[j].Im;
              OA->V[OA->N].Im =
                U.VRe[i][(OB.o[j* length + level] & NotCrAnMask)] * OB.V[j].Im -
                U.VIm[i][(OB.o[j* length + level] & NotCrAnMask)] * OB.V[j].Re;
            }
            else
            {                   /* transformation of an annihilation operator a_i = \sum_j U_{i,j}   a_j */
              OA->V[OA->N].Re =
                U.VRe[i][(OB.o[j * length + level] & NotCrAnMask)] * OB.V[j].Re -
                U.VIm[i][(OB.o[j * length + level] & NotCrAnMask)] * OB.V[j].Im;
              OA->V[OA->N].Im =
                U.VRe[i][(OB.o[j * length + level] & NotCrAnMask)] * OB.V[j].Im +
                U.VIm[i][(OB.o[j * length + level] & NotCrAnMask)] * OB.V[j].Re;
            }
            for (k = 0; k < length; k++)
            {
              OA->o[OA->N * length + k] = OB.o[j * length + k];
            }
            OA->o[OA->N * length + level] &= CrAnMask;
            OA->o[OA->N * length + level] |= (uint16_t)i;
            OA->N++;
          }
        }
      }
    }
    else
    {                           /* OB complex U real */
      OA->QComplex = 1;
      if ((OA->V = (Complex16*)malloc(OA->NAloc * sizeof(Complex16))) == NULL)
      {
        printf("Could not allocate memory for OA->V in MappedOperatorFixedLengthAtLevel\n");
        return 1;
      }
      if ((OA->o = (uint16_t*)malloc(OA->NAloc * length * sizeof(uint16_t))) == NULL)
      {
        printf("Could not allocate memory for OA->o in MappedOperatorFixedLengthAtLevel\n");
        return 1;
      }
      for (j = 0; j < OB.N; j++)
      {
        for (i = 0; i < U.Ni; i++)
        {
          if (((fabs(U.VRe[i][(OB.o[j * length + level] & NotCrAnMask)])) * (fabs(OB.V[j].Re) + fabs(OB.V[j].Im))) >
              (DBL_MIN * 10))
          {
            if (OA->N == OA->NAloc)
            {
              OA->NAloc *= 2;
              if ((OA->V = (Complex16*)realloc(OA->V, OA->NAloc * sizeof(Complex16))) == NULL)
              {
                printf("Could not reallocate memory for OA->V in MappedOperatorFixedLengthAtLevel\n");
                return 1;
              }
              if ((OA->o = (uint16_t*)realloc(OA->o, OA->NAloc * length * sizeof(uint16_t))) == NULL)
              {
                printf("Could not reallocate memory for OA->o in MappedOperatorFixedLengthAtLevel\n");
                return 1;
              }
            }
            OA->V[OA->N].Re = U.VRe[i][(OB.o[j * length + level] & NotCrAnMask)] * OB.V[j].Re;
            OA->V[OA->N].Im = U.VRe[i][(OB.o[j * length + level] & NotCrAnMask)] * OB.V[j].Im;
            for (k = 0; k < length; k++)
            {
              OA->o[OA->N * length + k] = OB.o[j * length + k];
            }
            OA->o[OA->N * length + level] &= CrAnMask;
            OA->o[OA->N * length + level] |= (uint16_t)i;
            OA->N++;
          }
        }
      }
    }
  }
  else
  {
    if (U.QComplex)
    {                           /* OB real U complex */
      OA->QComplex = 1;
      if ((OA->V = (Complex16*)malloc(OA->NAloc * sizeof(Complex16))) == NULL)
      {
        printf("Could not allocate memory for OA->V in MappedOperatorFixedLengthAtLevel\n");
        return 1;
      }
      if ((OA->o = (uint16_t*)malloc(OA->NAloc * length * sizeof(uint16_t))) == NULL)
      {
        printf("Could not allocate memory for OA->o in MappedOperatorFixedLengthAtLevel\n");
        return 1;
      }
      for (j = 0; j < OB.N; j++)
      {
        for (i = 0; i < U.Ni; i++)
        {
          if (((fabs(U.VRe[i][(OB.o[j * length + level] & NotCrAnMask)]) +
                fabs(U.VIm[i][(OB.o[j * length + level] & NotCrAnMask)])) * fabs(OB.VRe[j])) > (DBL_MIN * 10))
          {
            if (OA->N == OA->NAloc)
            {
              OA->NAloc *= 2;
              if ((OA->V = (Complex16*)realloc(OA->V, OA->NAloc * sizeof(Complex16))) == NULL)
              {
                printf("Could not reallocate memory for OA->V in MappedOperatorFixedLengthAtLevel\n");
                return 1;
              }
              if ((OA->o = (uint16_t*)realloc(OA->o, OA->NAloc * length * sizeof(uint16_t))) == NULL)
              {
                printf("Could not allocate memory for OA->o in MappedOperatorFixedLengthAtLevel\n");
                return 1;
              }
            }
            if (OB.o[j * length + level] & CrAnMask)
            {                   /* transformation of a creation operator a_j^{\dag} = \sum_i U_{i,j}^* a_i^{dag} */
              OA->V[OA->N].Re =  U.VRe[i][(OB.o[j * length + level] & NotCrAnMask)] * OB.VRe[j];
              OA->V[OA->N].Im = -U.VIm[i][(OB.o[j * length + level] & NotCrAnMask)] * OB.VRe[j];
            }
            else
            {                   /* transformation of an annihilation operator a_j = \sum_i U_{i,j}   a_i */
              OA->V[OA->N].Re = U.VRe[i][(OB.o[j * length + level] & NotCrAnMask)] * OB.VRe[j];
              OA->V[OA->N].Im = U.VIm[i][(OB.o[j * length + level] & NotCrAnMask)] * OB.VRe[j];
            }
            for (k = 0; k < length; k++)
            {
              OA->o[OA->N * length + k] = OB.o[j * length + k];
            }
            OA->o[OA->N * length + level] &= CrAnMask;
            OA->o[OA->N * length + level] |= (uint16_t)i;
            OA->N++;
          }
        }
      }
    }
    else
    {                           /* OB and U real */
      OA->QComplex = 0;
      if ((OA->VRe = (double*)malloc(OA->NAloc * sizeof(double))) == NULL)
      {
        printf("Could not allocate memory for OA->VRe in MappedOperatorFixedLengthAtLevel\n");
        return 1;
      }
      if ((OA->o = (uint16_t*)malloc(OA->NAloc * length * sizeof(uint16_t))) == NULL)
      {
        printf("Could not allocate memory for OA->o in MappedOperatorFixedLengthAtLevel\n");
        return 1;
      }
      for (j = 0; j < OB.N; j++)
      {
        for (i = 0; i < U.Ni; i++)
        {
          if (fabs(U.VRe[i][(OB.o[j * length + level] & NotCrAnMask)] * OB.VRe[j]) > (DBL_MIN * 10))
          {
            if (OA->N == OA->NAloc)
            {
              OA->NAloc *= 2;
              if ((OA->VRe = (double*)realloc(OA->VRe, OA->NAloc * sizeof(double))) == NULL)
              {
                printf("Could not reallocate memory for OA->VRe in MappedOperatorFixedLengthAtLevel\n");
                return 1;
              }
              if ((OA->o = (uint16_t*)realloc(OA->o, OA->NAloc * length * sizeof(uint16_t))) == NULL)
              {
                printf("Could not reallocate memory for OA->o in MappedOperatorFixedLengthAtLevel\n");
                return 1;
              }
            }
            OA->VRe[OA->N] = U.VRe[i][(OB.o[j * length + level] & NotCrAnMask)] * OB.VRe[j];
            for (k = 0; k < length; k++)
            {
              OA->o[OA->N * length + k] = OB.o[j * length + k];
            }
            OA->o[OA->N * length + level] &= CrAnMask;
            OA->o[OA->N * length + level] |= (uint16_t)i;
            OA->N++;
          }
        }
      }
    }
  }
  return 0;
}


/* OA <- OA + U OB
 *  U is a (NFermionB+NBosonB)*(NFermionA+NBosonA) demensional matrix that defines the transformation
 *  U.Ni = (NFermionA+NBosonA)
 *  U.Nj = (NFermionB+NBosonB)
 *  a_i        = \sum_j U_{i,j}   a_j
 *  a_i^{\dag} = \sum_j U_{i,j}^* a_j^{dag}
 */
uint32_t MappedOperatorAdd(OperatorType * const OA, const OperatorType OB, const MatrixType U)
{
  uint32_t level, length, i, k;
  OperatorFixedLengthType FA, FB;
  for (length = 0; length <= OB.MaxLength; length++)
  {
    for (i = 0; i < OB.l[length].N; i++)
    {
      FA.N = 1;
      FA.NAloc = 1;
      FA.QComplex = OB.l[length].QComplex;
      FA.o = (uint16_t*)malloc((length + 1) * sizeof(uint16_t));
      if (FA.QComplex)
      {
        FA.V = (Complex16*)malloc(sizeof(Complex16));
        FA.V[0].Re = OB.l[length].V[i].Re;
        FA.V[0].Im = OB.l[length].V[i].Im;
      }
      else
      {
        FA.VRe = (double*)malloc(sizeof(double));
        FA.VRe[0] = OB.l[length].VRe[i];
      }
      for (level = 0; level < length; level++)
      {
        FA.o[level] = OB.l[length].o[i * length + level];
      }
      for (level = 0; level < length; level++)
      {
        if (MappedOperatorFixedLengthAtLevel(&FB, FA, U, level, length))
        {
          printf("MappedOperatorFixedLengthAtLevel failed in MappedOperatorAdd\n");
          fflush(stdout);
          return 1;
        }
        free(FA.o);
        (FA.QComplex) ? free(FA.V) : free(FA.VRe);
        FA.NAloc = FB.NAloc;
        FA.N = FB.N;
        FA.QComplex = FB.QComplex;
        FA.o = FB.o;
        FA.V = FB.V;
        FA.VRe = FB.VRe;
        FB.NAloc = 0;
        FB.N = 0;
        FB.QComplex = 0;
      }
      if (FA.QComplex)
      {
        if (OA->l[length].QComplex == 0)
        {
          if (MakeRealOperatorFixedLengthComplex(OA, length))
          {
            printf("MakeRealOperatorFixedLengthComplex failed in MappedOperatorAdd\n");
            fflush(stdout);
            return 1;
          }
        }
        for (k = 0; k < FA.N; k++)
        {
          if (ComplexOperatorAddLader(OA, &(FA.o[k * length]), length, FA.V[k].Re, FA.V[k].Im))
          {
            printf("ComplexOperatorAddLader failed in MappedOperatorAdd\n");
            fflush(stdout);
            return 1;
          }
        }
      }
      else
      {                         /* FA real */
        if (OA->l[length].QComplex)
        {
          for (k = 0; k < FA.N; k++)
          {
            if (ComplexOperatorAddLader(OA, &(FA.o[k * length]), length, FA.VRe[k], 0))
            {
              printf("ComplexOperatorAddLader failed in MappedOperatorAdd\n");
              fflush(stdout);
              return 1;
            }
          }

        }
        else
        {
          for (k = 0; k < FA.N; k++)
          {
            if (RealOperatorAddLader(OA, &(FA.o[k * length]), length, FA.VRe[k]))
            {
              printf("RealOperatorAddLader failed in MappedOperatorAdd\n");
              fflush(stdout);
              return 1;
            }
          }
        }
      }
      free(FA.o);
      (FA.QComplex) ? free(FA.V) : free(FA.VRe);
      FA.N = 0;
      FA.NAloc = 0;
      FA.QComplex = 0;
    }
  }
  return OperatorRemoveSmallValues(OA, DBL_MIN * 10);
}

/* OA <- OA + U OB
 *  U is a (NFermionB+NBosonB)*(NFermionA+NBosonA) dimensional matrix that defines the transformation
 *  U.Ni = (NFermionA+NBosonA)
 *  U.Nj = (NFermionB+NBosonB)
 *  a_i        = \sum_j U_{i,j}   a_j
 *  a_i^{\dag} = \sum_j U_{i,j}^* a_j^{dag}
 */
uint32_t MappedOperatorAddNormalOrder(OperatorType * const OA, const OperatorType OB, const MatrixType U, uint32_t(*QDetNotOrdered) (uint16_t, uint16_t, uint32_t))
{
  uint32_t level, length, i, k;
  OperatorFixedLengthType FA, FB;
  for(i=OA->MaxLength+1;i<=OB.MaxLength;i++)
  {
    if(OB.l[i].N>0)
    {
      InitOperatorFixedLength(OA, i, OB.l[i].NBitsKey);
    }
  }
  for (length = 0; length <= OB.MaxLength; length++)
  {
    for (i = 0; i < OB.l[length].N; i++)
    {
      FA.N = 1;
      FA.NAloc = 1;
      FA.QComplex = OB.l[length].QComplex;
      FA.o = (uint16_t*)malloc((length + 1) * sizeof(uint16_t));
      if (FA.QComplex)
      {
        FA.V = (Complex16*)malloc(sizeof(Complex16));
        FA.V[0].Re = OB.l[length].V[i].Re;
        FA.V[0].Im = OB.l[length].V[i].Im;
      }
      else
      {
        FA.VRe = (double*)malloc(sizeof(double));
        FA.VRe[0] = OB.l[length].VRe[i];
      }
      for (level = 0; level < length; level++)
      {
        FA.o[level] = OB.l[length].o[i * length + level];
      }
      for (level = 0; level < length; level++)
      {
        if (MappedOperatorFixedLengthAtLevel(&FB, FA, U, level, length))
        {
          printf("MappedOperatorFixedLengthAtLevel failed in MappedOperatorAddNormalOrder\n");
          fflush(stdout);
          return 1;
        }
        free(FA.o);
        (FA.QComplex) ? free(FA.V) : free(FA.VRe);
        FA.NAloc = FB.NAloc;
        FA.N = FB.N;
        FA.QComplex = FB.QComplex;
        FA.o = FB.o;
        FA.V = FB.V;
        FA.VRe = FB.VRe;
        FB.NAloc = 0;
        FB.N = 0;
        FB.QComplex = 0;
      }
      if (FA.QComplex)
      {
        if (OA->l[length].QComplex == 0)
        {
          if (MakeRealOperatorFixedLengthComplex(OA, length))
          {
            printf("MakeRealOperatorFixedLengthComplex failed in MappedOperatorAddNormalOrder\n");
            fflush(stdout);
            return 1;
          }
        }
        for (k = 0; k < FA.N; k++)
        {
          if (ComplexOperatorAddLaderNormalOrder
              (OA, &(FA.o[k * length]), length, FA.V[k].Re, FA.V[k].Im, QDetNotOrdered))
          {
            printf("ComplexOperatorAddLaderNormalOrder failed in MappedOperatorAddNormalOrder\n");
            fflush(stdout);
            return 1;
          }
        }
      }
      else
      {                         /* FA real */
        if (OA->l[length].QComplex)
        {
          for (k = 0; k < FA.N; k++)
          {
            if (ComplexOperatorAddLaderNormalOrder(OA, &(FA.o[k * length]), length, FA.VRe[k], 0.0, QDetNotOrdered))
            {
              printf("ComplexOperatorAddLaderNormalOrder failed in MappedOperatorAddNormalOrder\n");
              fflush(stdout);
              return 1;
            }
          }
        }
        else
        {
          for (k = 0; k < FA.N; k++)
          {
            if (RealOperatorAddLaderNormalOrder(OA, &(FA.o[k * length]), length, FA.VRe[k], QDetNotOrdered))
            {
              printf("RealOperatorAddLaderNormalOrder failed in MappedOperatorAddNormalOrder\n");
              fflush(stdout);
              return 1;
            }
          }
        }
      }
      free(FA.o);
      (FA.QComplex) ? free(FA.V) : free(FA.VRe);
      FA.N = 0;
      FA.NAloc = 0;
      FA.QComplex = 0;
    }
  }
  return OperatorRemoveSmallValues(OA, DBL_MIN * 10);
}

/* OA <- U OB
 *  U is a (NFermionB+NBosonB)*(NFermionA+NBosonA) dimensional matrix that defines the transformation
 *  U.Ni = (NFermionA+NBosonA)
 *  U.Nj = (NFermionB+NBosonB)
 *  a_i        = \sum_j U_{i,j}   a_j
 *  a_i^{\dag} = \sum_j U_{i,j}^* a_j^{dag}
 */
uint32_t RotateOperator(OperatorType * const OA, const OperatorType OB, const MatrixType U)
{
  if(InitOperatorCopySettings(OA, OB))
  {
    printf("InitOperatorCopySettings failed in RotateOperator\n");
    return 1;
  }
  OperatorChangeNFandNB(OA, U.Ni, OA->NBoson);
  return MappedOperatorAdd(OA, OB, U);
}
uint32_t RotateOperatorNormalOrder(OperatorType * const OA, const OperatorType OB, const MatrixType U, uint32_t(*QDetNotOrdered) (uint16_t, uint16_t, uint32_t))
{
  if(InitOperatorCopySettings(OA, OB))
  {
    printf("InitOperatorCopySettings failed in RotateOperator\n");
    return 1;
  }
  OperatorChangeNFandNB(OA, U.Ni, OA->NBoson);
  return MappedOperatorAddNormalOrder(OA, OB, U, QDetNotOrdered);
}


/*  OT == O^{T}
 */
uint32_t OperatorTranspose(OperatorType O, OperatorType * OT)
{
  uint32_t il, i, j;
  uint16_t oadd[O.MaxLength];
  if (InitOperatorCopySettings(OT, O))
    return 1;
  for (il = 0; il <= O.MaxLength; il++)
  {
    if (O.l[il].QComplex)
    {
      for (i = 0; i < O.l[il].N; i++)
      {
        for (j = 0; j < il; j++)
        {
          oadd[il - j - 1] = O.l[il].o[j + i * il] ^ CrAnMask;
        }
        if (ComplexOperatorAddLader(OT, oadd, il, O.l[il].V[i].Re, O.l[il].V[i].Im))
          return 1;
      }
    }
    else
    {
      for (i = 0; i < O.l[il].N; i++)
      {
        for (j = 0; j < il; j++)
        {
          oadd[il - j - 1] = O.l[il].o[j + i * il] ^ CrAnMask;
        }
        if (RealOperatorAddLader(OT, oadd, il, O.l[il].VRe[i]))
          return 1;
      }
    }
  }
  SetQComplexOperator(OT);
  return 0;
}

/*  OT == O^{T}
 */
uint32_t OperatorTransposeNormalOrder(OperatorType O, OperatorType * OT, uint32_t(*QDetNotOrdered) (uint16_t, uint16_t, uint32_t))
{
  uint32_t il, i, j;
  uint16_t oadd[O.MaxLength];
  if (InitOperatorCopySettings(OT, O))
    return 1;
  for (il = 0; il <= O.MaxLength; il++)
  {
    if (O.l[il].QComplex)
    {
      for (i = 0; i < O.l[il].N; i++)
      {
        for (j = 0; j < il; j++)
        {
          oadd[il - j - 1] = O.l[il].o[j + i * il] ^ CrAnMask;
        }
        if (ComplexOperatorAddLaderNormalOrder(OT, oadd, il, O.l[il].V[i].Re, O.l[il].V[i].Im,QDetNotOrdered))
          return 1;
      }
    }
    else
    {
      for (i = 0; i < O.l[il].N; i++)
      {
        for (j = 0; j < il; j++)
        {
          oadd[il - j - 1] = O.l[il].o[j + i * il] ^ CrAnMask;
        }
        if (RealOperatorAddLaderNormalOrder(OT, oadd, il, O.l[il].VRe[i],QDetNotOrdered))
          return 1;
      }
    }
  }
  SetQComplexOperator(OT);
  return 0;
}

/*  O -> O^{*}
 */
uint32_t OperatorComplexConjugate(OperatorType * O)
{
  uint32_t il, i;
  for (il = 0; il <= O->MaxLength; il++)
  {
    if (O->l[il].QComplex)
    {
      for (i = 0; i < O->l[il].N; i++)
      {
        O->l[il].V[i].Im = -O->l[il].V[i].Im;
      }
    }
  }
  return 0;
}


/*  OT == O^{+}
 */
uint32_t OperatorConjugateTranspose(OperatorType O, OperatorType * OT)
{
  uint32_t il, i, j;
  uint16_t oadd[O.MaxLength];
  if (InitOperatorCopySettings(OT, O))
    return 1;
  for (il = 0; il <= O.MaxLength; il++)
  {
    if (O.l[il].QComplex)
    {
      for (i = 0; i < O.l[il].N; i++)
      {
        for (j = 0; j < il; j++)
        {
          oadd[il - j - 1] = O.l[il].o[j + i * il] ^ CrAnMask;
        }
        if (ComplexOperatorAddLader(OT, oadd, il, O.l[il].V[i].Re, -O.l[il].V[i].Im))
          return 1;
      }
    }
    else
    {
      for (i = 0; i < O.l[il].N; i++)
      {
        for (j = 0; j < il; j++)
        {
          oadd[il - j - 1] = O.l[il].o[j + i * il] ^ CrAnMask;
        }
        if (RealOperatorAddLader(OT, oadd, il, O.l[il].VRe[i]))
          return 1;
      }
    }
  }
  SetQComplexOperator(OT);
  return 0;
}


/*  OT == O^{+}
 */
uint32_t OperatorConjugateTransposeNormalOrder(OperatorType O, OperatorType * OT, uint32_t(*QDetNotOrdered) (uint16_t, uint16_t, uint32_t))
{
  uint32_t il, i, j;
  uint16_t oadd[O.MaxLength];
  if (InitOperatorCopySettings(OT, O)) return 1;
  
  for (il = 0; il <= O.MaxLength; il++)
  {
    if (O.l[il].QComplex)
    {
      for (i = 0; i < O.l[il].N; i++)
      {
        for (j = 0; j < il; j++)
        {
          oadd[il - j - 1] = O.l[il].o[j + i * il] ^ CrAnMask;
        }
        if (ComplexOperatorAddLaderNormalOrder(OT, oadd, il, O.l[il].V[i].Re, -O.l[il].V[i].Im, QDetNotOrdered))
        {
          FreeOperator(OT);
          return 1;
        }
      }
    }
    else
    {
      for (i = 0; i < O.l[il].N; i++)
      {
        for (j = 0; j < il; j++)
        {
          oadd[il - j - 1] = O.l[il].o[j + i * il] ^ CrAnMask;
        }
        if (RealOperatorAddLaderNormalOrder(OT, oadd, il, O.l[il].VRe[i], QDetNotOrdered))
        {
          FreeOperator(OT);
          return 1;
        }
      }
    }
  }
  SetQComplexOperator(OT);
  return 0;
}

/*  function 
 *  returns 0 if operator is not Hermitian
 *  returns 1 if operator is Hermitian
 *  returns 2 if operator is Hermitian after Chop(H,epsilon), i.e. small non Hermitian parts are allowed
 *  returns 3 if failed
 */
uint32_t QOperatorHermitian(OperatorType O, const double epsilon)
{
  uint32_t il, Hermitian=1;
  OperatorType OT;
  if (OperatorConjugateTransposeNormalOrder(O, &OT, QDetNotOrderedFcFaBaBc))
  {
    printf("OperatorConjugateTransposeNormalOrder failed in QOperatorHermitian\n");
    FreeOperator(&OT);
    fflush(stdout);
    return 3;
  };
  if (OperatorAddNormalOrder(&OT, O, -1, 0, QDetNotOrderedFcFaBaBc))
  {
    printf("OperatorAddNormalOrder failed in QOperatorHermitian\n");
    FreeOperator(&OT);
    fflush(stdout);
    return 3;
  };
  if (OperatorChop(&OT, 100*DBL_MIN))
  {
    printf("OperatorChop failed in QOperatorHermitian\n");
    FreeOperator(&OT);
    fflush(stdout);
    return 3;
  };
  for (il = 0; il <= OT.MaxLength; il++)
  {
    if (OT.l[il].N != 0) Hermitian = 0;
  }
  if (Hermitian==1)
  {
    FreeOperator(&OT);
    return 1;
  }
  Hermitian=2;
  if (OperatorChop(&OT, epsilon))
  {
    printf("OperatorChop failed in QOperatorHermitian\n");
    FreeOperator(&OT);
    fflush(stdout);
    return 3;
  };
  for (il = 0; il <= OT.MaxLength; il++)
  {
    if (OT.l[il].N != 0) Hermitian = 0;
  }
  if (Hermitian==2)
  {
    FreeOperator(&OT);
    return 2;
  }
  FreeOperator(&OT);
  return 0;
}

bool isOperatorAlmostHermitian(const OperatorType &O, const double epsilon)
{
	int q = QOperatorHermitian(O, epsilon);
	return (q == 1 || q == 2);
}

uint32_t QOperatorEqual(const OperatorType *l, const OperatorType *r, const double epsilon)
{
  OperatorType diff;
  uint32_t il;
  if(OperatorCopy(&diff, *l))
  {
    printf("OperatorCopy failed in QOperatorEqual\n");
    FreeOperator(&diff);
    fflush(stdout);
    return 3;
  }
  if (OperatorAddNormalOrder(&diff, *r, -1, 0, QDetNotOrderedFcFaBaBc))
  {
    printf("OperatorAddNormalOrder failed in QOperatorEqual\n");
    FreeOperator(&diff);
    fflush(stdout);
    return 3;
  };
  if (OperatorChop(&diff, epsilon))
  {
    printf("OperatorChop failed in QOperatorEqual\n");
    FreeOperator(&diff);
    fflush(stdout);
    return 3;
  };
  for (il = 0; il <= diff.MaxLength; il++)
  {
    if (diff.l[il].N != 0)
    {
      FreeOperator(&diff);
      return 0;
    }
  }
  FreeOperator(&diff);
  return 1;
}

/*
 *
 */
uint32_t OperatorChangeNFandNB(OperatorType *O, const uint32_t NFnew, uint32_t NBnew)
{
  uint32_t l,j;
  int i,shift;
  if(NFnew+NBnew>65536)
  {
    printf("Error in ChangeNFandNB NFnew (%i) + NBnew (%i) larger than 65536 (%i)\n",NFnew,NBnew,NFnew+NBnew);
    return 1;
  }
  // remove fermions from operator when needed. (leave NFermion to the old value)
  if(NFnew<O->NFermion)
  {
    for(l=1;l<=O->MaxLength;l++)
    {
      if(O->l[l].QComplex)
      {
#pragma omp parallel for private(j)
        for(i=0;i<O->l[l].N;i++)
        {
          for(j=i*l;j<i*l+l;j++)
          {
            if(((O->l[l].o[j])&NotCrAnMask)>=NFnew)
            {
              O->l[l].V[i].Re=0;
              O->l[l].V[i].Im=0;
              break;
            }
          }
        }
      }
      else
      {
#pragma omp parallel for private(j)
        for(i=0;i<O->l[l].N;i++)
        {
          for(j=i*l;j<i*l+l;j++)
          {
            if(((O->l[l].o[j])&NotCrAnMask)>=NFnew)
            {
              O->l[l].VRe[i]=0;
              break;
            }
          }
        }
      }
    }
  }
  // remove bosons from operator when needed. (leave NBoson to the old value)
  if(NBnew<O->NBoson)
  {
    for(l=1;l<=O->MaxLength;l++)
    {
      if(O->l[l].QComplex)
      {
#pragma omp parallel for private(j)
        for(i=0;i<O->l[l].N;i++)
        {
          for(j=i*l;j<i*l+l;j++)
          {
            if(((O->l[l].o[j])&NotCrAnMask)>=NBnew+O->NFermion)
            {
              O->l[l].V[i].Re=0;
              O->l[l].V[i].Im=0;
              break;
            }
          }
        }
      }
      else
      {
#pragma omp parallel for private(j)
        for(i=0;i<O->l[l].N;i++)
        {
          for(j=i*l;j<i*l+l;j++)
          {
            if(((O->l[l].o[j])&NotCrAnMask)>=NBnew+O->NFermion)
            {
              O->l[l].VRe[i]=0;
              break;
            }
          }
        }
      }
    }
  }
  OperatorChop(O, 10*DBL_MIN);
  // shift the bosons when needed
  if(NFnew!=O->NFermion)
  {
    shift=(int)NFnew-(int)O->NFermion;
#pragma omp parallel for private(j)
    for(l=1;l<=O->MaxLength;l++)
    {
      for(i=0;i<O->l[l].N;i++)
      {
        for(j=i*l;j<i*l+l;j++)
        {
          if(((O->l[l].o[j])&NotCrAnMask)>=O->NFermion)
          {
            O->l[l].o[j]+=shift;
          }
        }
      }
    }
  }
  // set the number of bosons and Fermions
  O->NFermion=NFnew;
  O->NBoson=NBnew;
  O->NWord = Min(O->NFermion, 1 + (O->NFermion - 1) / 8) + O->NBoson;
  return 0;
}

/*  
 *
 */
bool isCrAnIncluded(const uint16_t * const CrAn, const uint length, const vector<vector<uint16_t> > &subspace)
{
	for (uint i = 0; i < (uint) subspace.size(); ++i)
	{
		for (uint j = 0; j < length; ++j)
		{
			const uint16_t AnInd = NotCrAnMask & CrAn[j];
			if (isElementOf(AnInd, subspace[i])) return true;
		}
	}
	return false;
}

bool isCrAnExcluded(const uint16_t * const CrAn, const uint length, const vector<vector<uint16_t> > &subspace)
{
	for (uint i = 0; i < (uint) subspace.size(); ++i)
	{
		for (uint j = 0; j < length; ++j)
		{
			const uint16_t AnInd = NotCrAnMask & CrAn[j];
			if (isElementOf(AnInd, subspace[i])) return false;
		}
	}
	return true;
}

bool isCrAnConservesSubspaceOccupation(const uint16_t * const CrAn, const uint length, const vector<vector<uint16_t> > &subspace)
{
	for (uint i = 0; i < (uint) subspace.size(); ++i)
	{
		int additionalSubspaceElectrons = 0;
		for (uint j = 0; j < length; ++j)
		{
			const uint16_t AnInd = NotCrAnMask & CrAn[j];
			if (isElementOf(AnInd, subspace[i]))
			{
				if (AnInd == CrAn[j])
					--additionalSubspaceElectrons;
				else
					++additionalSubspaceElectrons;
			}
		}
		if (additionalSubspaceElectrons != 0) return false;
	}
	return true;
}

void extractPartialOperatorFixedLength(OperatorType * const O, OperatorType * const OPartial, const vector<vector<uint16_t> > &indices, bool (*criterion)(const uint16_t * const, const uint, const vector<vector<uint16_t> >&), const uint length, bool changeInput)
{
	uint16_t *oadd;
	InitOperatorFixedLength(OPartial, length, O->l[length].NBitsKey);
	for (uint i = 0; i < (O->l[length].N); ++i)
	{
    oadd = &(O->l[length].o[i * length]);
		if (criterion(oadd, length, indices))
		{
			if (O->l[length].QComplex)
			{
				ComplexOperatorAddLader(OPartial, oadd, length, O->l[length].V[i].Re, O->l[length].V[i].Im);
				if (changeInput)
				{
					O->l[length].V[i].Re = 0;
					O->l[length].V[i].Im = 0;
				}
			}
			else
			{
				RealOperatorAddLader(OPartial, oadd, length, O->l[length].VRe[i]);
				if (changeInput) O->l[length].VRe[i] = 0;
			}
		}
	}
}

void extractPartialOperatorFixedLength(OperatorType * const O, OperatorType * const OPartial, const vector<vector<uint16_t> > &indices, bool (*criterion)(const uint16_t * const, const uint, const vector<vector<uint16_t> >&), const uint length)
{
	extractPartialOperatorFixedLength(O, OPartial, indices, criterion, length, false);
}


/*  Input:
 *    Operator O
 *    bool changeInput: if true O+Opartial = O on input
 *                      if false O not changed
 *    indices: vector of vector of uint16 { {1,2}, {4,5}, {3,6,7,9}, ...} For each list the function criterion must yield true for the determinant to be included in Opartial
 *    criterion: Pointer to function with input:
 *                 Pointer to single determinant
 *                 length of determinant
 *                 indices
 *               Output true or false depending if the determinant needs to be included in Opartial or not
 *  Output Operator OPartial, not initialized on input
 *         Operator O, if change Input == true change O such that Oout + Opartial = Oin
 */
void extractPartialOperator(OperatorType * const O, OperatorType * const OPartial, const vector<vector<uint16_t> > &indices, bool (*criterion)(const uint16_t * const, const uint, const vector<vector<uint16_t> >&), bool changeInput)
{
	InitOperator(OPartial);
	OPartial->NFermion = O->NFermion;
	OPartial->NBoson = O->NBoson;
  OPartial->NWord = O->NWord;
	if (criterion(NULL, 0, indices) && O->l[0].N == 1)
	{
		if (O->l[0].QComplex) ComplexOperatorAddLader(OPartial, NULL, 0, O->l[0].V[0].Re, O->l[0].V[0].Im);
		else RealOperatorAddLader(OPartial, NULL, 0, O->l[0].VRe[0]);
	}
	for (uint length = 1; length <= O->MaxLength; ++length)
  {
		extractPartialOperatorFixedLength(O, OPartial, indices, criterion, length, changeInput);
  }
	SetQComplexOperator(OPartial);
}

void extractPartialOperator(OperatorType * const O, OperatorType * const OPartial, const vector<vector<uint16_t> > &indices, bool (*criterion)(const uint16_t * const, const uint, const vector<vector<uint16_t> >&))
{
	extractPartialOperator(O, OPartial, indices, criterion, false);
}
