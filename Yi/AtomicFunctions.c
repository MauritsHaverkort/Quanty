#include "AtomicFunctions.h"

/*******************************
            Math
 *******************************/

/* function to generate combinations of C(M, N) */

int twiddle(int *x, int *y, int *z, int *p)
{
  register int i, j, k;
  
  j = 1;
  while (p[j] <= 0)
    j++;
  if(p[j-1] == 0)
  {
    for(i = j-1; i != 1; i--)
      p[i] = -1;
    p[j] = 0;
    *x = *z = 0;
    p[1] = 1;
    *y = j-1;
  }
  else
  {
    if (j > 1)
      p[j-1] = 0;
    do
    {
      j++;
    } while(p[j] > 0);
    k = j-1;
    i = j;
    while (p[i] == 0)
      p[i++] = -1;
    if (p[i] == -1)
    {
      p[i] = p[k];
      *z = p[k]-1;
      *x = i-1;
      *y = k-1;
      p[k] = -1;
    }
    else
    {
      if (i == p[0])
        return 1;
      else
      {
        p[j] = p[i];
        *z = p[i]-1;
        p[i] = 0;
        *x = j-1;
        *y = i-1;
      }
    }
  }
  
  return 0;
}

void inittwiddle(int m, int n, int *p)
{
  int i;
  p[0] = n+1;
  for(i = 1; i != n-m+1; i++)
    p[i] = 0;
  while(i != n+1)
  {
    p[i] = i+m-n;
    i++;
  }
  p[n+1] = -2;
  if(m == 0)
    p[1] = 1;
}


/*******************************
            IO
 *******************************/

typedef struct RadMEList
{
  char   Name[MAXStrLength];    /* A pointer to a name string that can be used to lable functions, good for printing etc. */
  uint32_t N;                   /* dimension of list */
  int *n1;
  int *l1;
  int *n2;
  int *l2;
  int *n3;
  int *l3;
  int *n4;
  int *l4;
  int *ll;
  double *radme;
}
RadMEList;

uint32_t InitRadMEList(RadMEList * M)
{
  strcpy(M->Name, "List of Radial Matrix Element");
  if (M->N == 0)
  {
    printf("Initialized a ListOfPoles with size zero.\n");
    fflush(stdout);
    return 2;
  };
  
  M->n1 = (int*)malloc((M->N) * sizeof(int));
  if (M->n1 == NULL)
  {
    printf("malloc for M->n1 failed in ListOfRadMEList\n");
    fflush(stdout);
    return 1;
  };
  M->l1 = (int*)malloc((M->N) * sizeof(int));
  if (M->l1 == NULL)
  {
    printf("malloc for M->l1 failed in ListOfRadMEList\n");
    fflush(stdout);
    return 1;
  };
  M->n2 = (int*)malloc((M->N) * sizeof(int));
  if (M->n2 == NULL)
  {
    printf("malloc for M->n2 failed in ListOfRadMEList\n");
    fflush(stdout);
    return 1;
  };
  M->l2 = (int*)malloc((M->N) * sizeof(int));
  if (M->l2 == NULL)
  {
    printf("malloc for M->l2 failed in ListOfRadMEList\n");
    fflush(stdout);
    return 1;
  };
  M->n3 = (int*)malloc((M->N) * sizeof(int));
  if (M->n3 == NULL)
  {
    printf("malloc for M->n3 failed in ListOfRadMEList\n");
    fflush(stdout);
    return 1;
  };
  M->l3 = (int*)malloc((M->N) * sizeof(int));
  if (M->l3 == NULL)
  {
    printf("malloc for M->l3 failed in ListOfRadMEList\n");
    fflush(stdout);
    return 1;
  };
  M->n4 = (int*)malloc((M->N) * sizeof(int));
  if (M->n4 == NULL)
  {
    printf("malloc for M->n4 failed in ListOfRadMEList\n");
    fflush(stdout);
    return 1;
  };
  M->l4 = (int*)malloc((M->N) * sizeof(int));
  if (M->l4 == NULL)
  {
    printf("malloc for M->l4 failed in ListOfRadMEList\n");
    fflush(stdout);
    return 1;
  };
  M->ll = (int*)malloc((M->N) * sizeof(int));
  if (M->ll == NULL)
  {
    printf("malloc for M->ll failed in ListOfRadMEList\n");
    fflush(stdout);
    return 1;
  };
  M->radme = (double*)malloc((M->N) * sizeof(double));
  if (M->radme == NULL)
  {
    printf("malloc for M->radme failed in ListOfRadMEList\n");
    fflush(stdout);
    return 1;
  };

  return 0;
}

uint32_t StreamReadRadMEList(RadMEList * O, FILE * stream)
{
  int i;
  
  fscanf(stream, "%i", &(O->N));
  if (ferror(stream))
  {
    printf("Error in StreamReadRadMEList error nr: %i\n", ferror(stream));
    perror("");
    fflush(stdout);
    return 1;
  };
  
  if (InitRadMEList(O))
  {
    printf("InitRadMEList failed in StreamReadRadMEList with O->N=%i\n", O->N);
    fflush(stdout);
    return 1;
  };
  
  for (i = 0; i < O->N; i++)
  {
    if (feof(stream))
    {
      printf("Error in StreamReadRadMEList eof reached at line %i of %i\n", i, O->N);
      return 1;
    };
    
    fscanf(stream, "%i %i %i %i %i %i %i %i %i %lf",
           &(O->n1[i]), &(O->l1[i]),
           &(O->n2[i]), &(O->l2[i]),
           &(O->n3[i]), &(O->l3[i]),
           &(O->n4[i]), &(O->l4[i]),
           &(O->ll[i]), &(O->radme[i]));
    
    if (ferror(stream))
    {
      printf("Error in StreamReadRadMEList error nr: %i\n", ferror(stream));
      perror("");
      return 1;
    };
  }
  
  return 0;
}


uint32_t FileReadRadMEList(RadMEList * O, const char *const flname)
{
  FILE  *stream;
  uint32_t retval;
  if ((stream = fopen(flname, "r")) == NULL)
  {
    printf("ERROR\nCan't open %s for reading.\n", flname);
    return 1;
  }
  retval = StreamReadRadMEList(O, stream);
  fclose(stream);
  
  return retval;
}

/*******************************
  wave-function and hamiltonian
 *******************************/

uint32_t InitWavefunctionRandom(WaveFunctionType * const psi,
                                uint32_t const NFermions,
                                uint32_t const Nelectron)
{
  int i, x, y, z, p[NFermions+2], b[NFermions], count=0;
  uint8_t *Determinant;
  
  psi->NBitsKey = 16;
  psi->QComplex = 0;
  psi->NBoson = 0;
  psi->NFermion = NFermions;
  if (WaveFunctionInit(psi))
  {
    printf("WaveFunctionInit failed in InitXASWavefunction\n");
    fflush(stdout);
    return 1;
  };
  sprintf(psi->Name, "WaveFunction");
  
  inittwiddle(Nelectron, NFermions, p);
  
  for(i = 0; i != NFermions-Nelectron; i++)
  {
    b[i] = 0;
  }
  while(i != NFermions)
  {
    b[i++] = 1;
  }
  Determinant = (uint8_t*)calloc(psi->NWord, sizeof(uint8_t));
  for (i = 0; i < NFermions; i++)
  {
    if (b[i])
    {
      OperatorCr(Determinant, i, psi->NFermion);
    }
  }
  RealWaveFunctionAddElement(psi, Determinant, 1.); free(Determinant);
  count++;
  
  while(!twiddle(&x, &y, &z, p))
  {
    Determinant = (uint8_t*)calloc(psi->NWord, sizeof(uint8_t));
    b[x] = 1;
    b[y] = 0;
    for (i = 0; i < NFermions; i++)
    {
      if (b[i])
      {
        OperatorCr(Determinant, i, psi->NFermion);
      }
    }
    RealWaveFunctionAddElement(psi, Determinant, 1.); free(Determinant);
    count++;
  }
  
	printf("\nThe initial wave function has %i determinants.\n", count);
  WaveFunctionRandomize(psi);
  
  //PrintWaveFunctionFullReadable(*psi);
  
  return 0;
  
}

/* The CI hamiltonian for the Atom */
/* The basis is ordered as :       */
/*
 1s              1
 2s 2p           2   3
 3s 3p 3d        4   5   6
 4s 4p 4d 4f     7   8   9   10
 5s 5p 5d 5f     11  12  13  14
 6s 6p 6d        15  16  17
 7s 7p           18  19
 8s              20
 */

uint32_t AtomicCI(const char *const flname, uint32_t Nelectron, uint32_t NOrb,
                  uint32_t * const BasisLen, OperatorType *O)
{
  int quantN[] = {1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 7, 7, 8};
  int quantL[] = {0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 0, 1, 0};
  uint32_t i, j, k, l, tabi;
  *BasisLen = 0;
  for (i=0; i<NOrb; i++)
  {
    *BasisLen += 2*quantL[i]+1;
  }
  
  uint16_t IndexUp[*BasisLen], IndexDn[*BasisLen], oadd[4];
  int IndexStart;
  int BasisN[*BasisLen], BasisL[*BasisLen], BasisM[*BasisLen];
  double En, coeff1, coeff2, me;
  RadMEList radmelist;
  
  FileReadRadMEList(&radmelist, flname);
  printf("\nRadial matrix element read from file: %s\n", flname);
  printf("with length: %i\n", radmelist.N);
  
  /* Create the basis set */
  
  for (i=0; i<*BasisLen; i++)
  {
    IndexUp[i] = 2*i;
    IndexDn[i] = 2*i + 1;
  }
  
  for (i=0; i<NOrb; i++)
  {
    IndexStart = 0;
    for (j=0; j<i; j++)
    {
      IndexStart += 2*quantL[j]+1;
    }
    for (j=0; j<2*quantL[i]+1; j++)
    {
      BasisN[IndexStart+j] =  quantN[i];
      BasisL[IndexStart+j] =  quantL[i];
      BasisM[IndexStart+j] = -quantL[i] + j;
    }
    
  }
  
  printf("\nBasis created as:\n");
  printf("Quantum number N: ");
  for (i=0; i<*BasisLen; i++)
  {
    printf(" %d ", BasisN[i]);
  }
  printf("\nQuantum number L: ");
  for (i=0; i<*BasisLen; i++)
  {
    printf(" %d ", BasisL[i]);
  }
  printf("\nQuantum number M: ");
  for (i=0; i<*BasisLen; i++)
  {
    if (BasisM[i]<0) {
      printf("%d ", BasisM[i]);
    } else {
      printf(" %d ", BasisM[i]);
    }
  }
  
  /* Create the Hamiltonian */
  
  InitNamedOperator(O, "Atomic Hamiltonian");
  printf("\n\nAtomic Hamiltonian for NOrb = %i, NFermion = %i ", NOrb, *BasisLen*2);
  O->NFermion = *BasisLen * 2;
  O->NWord = Min(O->NFermion, 1+(O->NFermion-1)/8);
  InitOperatorFixedLength(O, 2, 8);
  InitOperatorFixedLength(O, 4, 8);
  
  /*                           */
  /* "onsite" energy -Z^2/2n^2 */
  /*                           */
  
  for (i=0; i<*BasisLen; i++)
  {
    En = -(double)Nelectron*Nelectron/(2*BasisN[i]*BasisN[i]);
    oadd[1] = CrAnMask+2*i;
    oadd[0] = 2*i;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, En, &QDetNotOrderedFcFaBcBa);
    oadd[1] = CrAnMask+2*i+1;
    oadd[0] = 2*i+1;
    RealOperatorAddLaderNormalOrder(O, oadd, 2, En, &QDetNotOrderedFcFaBcBa);
  }
  
  /*               */
  /* Coulomb terms */
  /*               */
  
  /* iterate over all possibilities of <ab|cd> */
  for (i=0; i<*BasisLen; i++)
  {
    for (j=0; j<*BasisLen; j++)
    {
      for (k=0; k<*BasisLen; k++)
      {
        for (l=0; l<*BasisLen; l++)
        {
          me = 0.;
          
          for (tabi=0; tabi<radmelist.N; tabi++)
          {
            if (BasisM[i]+BasisM[j]==BasisM[k]+BasisM[l] &&
                Abs(BasisM[k]-BasisM[i]) <= radmelist.ll[tabi])
            {
              // <ab|cd>, <ba|dc>, <cd|ab>, <dc|ba> equivlanent for radial part
              if ((BasisN[i]==radmelist.n1[tabi] && BasisN[j]==radmelist.n2[tabi] &&
                   BasisN[k]==radmelist.n3[tabi] && BasisN[l]==radmelist.n4[tabi] &&
                   BasisL[i]==radmelist.l1[tabi] && BasisL[j]==radmelist.l2[tabi] &&
                   BasisL[k]==radmelist.l3[tabi] && BasisL[l]==radmelist.l4[tabi])||
                  (BasisN[i]==radmelist.n2[tabi] && BasisN[j]==radmelist.n1[tabi] &&
                   BasisN[k]==radmelist.n4[tabi] && BasisN[l]==radmelist.n3[tabi] &&
                   BasisL[i]==radmelist.l2[tabi] && BasisL[j]==radmelist.l1[tabi] &&
                   BasisL[k]==radmelist.l4[tabi] && BasisL[l]==radmelist.l3[tabi])||
                  (BasisN[i]==radmelist.n3[tabi] && BasisN[j]==radmelist.n4[tabi] &&
                   BasisN[k]==radmelist.n1[tabi] && BasisN[l]==radmelist.n2[tabi] &&
                   BasisL[i]==radmelist.l3[tabi] && BasisL[j]==radmelist.l4[tabi] &&
                   BasisL[k]==radmelist.l1[tabi] && BasisL[l]==radmelist.l2[tabi])||
                  (BasisN[i]==radmelist.n4[tabi] && BasisN[j]==radmelist.n3[tabi] &&
                   BasisN[k]==radmelist.n2[tabi] && BasisN[l]==radmelist.n1[tabi] &&
                   BasisL[i]==radmelist.l4[tabi] && BasisL[j]==radmelist.l3[tabi] &&
                   BasisL[k]==radmelist.l2[tabi] && BasisL[l]==radmelist.l1[tabi]))
              {
                if (radmelist.ll[tabi] <= Min(BasisL[i]+BasisL[k], BasisL[j]+BasisL[l]) &&
                    radmelist.ll[tabi] >= Max(Abs(BasisL[i]-BasisL[k]), Abs(BasisL[j]-BasisL[l])))
                {
                  coeff1 = sqrt((2*BasisL[i]+1)*(2*BasisL[j]+1)*
                                (2*BasisL[k]+1)*(2*BasisL[l]+1));
                  coeff2 = pow(-1, BasisM[j]+BasisM[k]) * coeff1 *
                  ThreeJSymbol(radmelist.ll[tabi], 0, BasisL[i], 0, BasisL[k], 0) *
                  ThreeJSymbol(radmelist.ll[tabi], 0, BasisL[l], 0, BasisL[j], 0) * radmelist.radme[tabi];
                  me +=  coeff2 * ThreeJSymbol(radmelist.ll[tabi], BasisM[k]-BasisM[i], BasisL[i], BasisM[i], BasisL[k], -BasisM[k]) *
                  ThreeJSymbol(radmelist.ll[tabi], BasisM[k]-BasisM[i], BasisL[l], BasisM[l], BasisL[j], -BasisM[j]);
                }
              }
              
            }
          }
          
          me *= 0.5;
          
          if (me != 0.)
          {
            // <up, up | up, up>
            oadd[3] = CrAnMask+2*i;
            oadd[2] = CrAnMask+2*j;
            oadd[1] = 2*l;
            oadd[0] = 2*k;
            RealOperatorAddLaderNormalOrder(O, oadd, 4, me, &QDetNotOrderedFcFaBcBa);
            
            // <dn, dn | dn, dn>
            oadd[3] = CrAnMask+2*i+1;
            oadd[2] = CrAnMask+2*j+1;
            oadd[1] = 2*l+1;
            oadd[0] = 2*k+1;
            RealOperatorAddLaderNormalOrder(O, oadd, 4, me, &QDetNotOrderedFcFaBcBa);
            
            
            // <up, dn | up, dn>
            oadd[3] = CrAnMask+2*i;
            oadd[2] = CrAnMask+2*j+1;
            oadd[1] = 2*l+1;
            oadd[0] = 2*k;
            RealOperatorAddLaderNormalOrder(O, oadd, 4, me, &QDetNotOrderedFcFaBcBa);
            
            // <dn, up | dn, up>
            oadd[3] = CrAnMask+2*i+1;
            oadd[2] = CrAnMask+2*j;
            oadd[1] = 2*l;
            oadd[0] = 2*k+1;
            RealOperatorAddLaderNormalOrder(O, oadd, 4, me, &QDetNotOrderedFcFaBcBa);
            
          }
        }
      }
    }
  }
  
//  PrintOperatorReadable(*O);
  
  return 0;
}
