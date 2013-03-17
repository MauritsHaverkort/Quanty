/*    Within this file the Green function is always represented as the resolvent of some matrix
 *    
 *    1) Tri-diagonal matrix representation
 *
 *     ( A[0][0]=a1 | A[0][1]=b1 | 0          | 0          )
 *   A=( A[1][0]=b1 | A[1][1]=a2 | A[1][2]=b2 | 0          )
 *     ( 0          | A[2][1]=b2 | A[2][2]=a3 | A[2][3]=b3 )
 *     ( 0          | 0          | A[3][2]=b3 | A[3][3]=a4 )
 *
 *   G(w) = a0 + b0^2 / (w - a1 - b1^2 / (w - a2 - b2^2 / (w - a3 - b3^2 / (w ... ))))
 *
 *
 *    2) Anderson matrix representation
 *
 *     ( A[0][0]=a1 | A[0][1]=b1 | A[0][2]=b2 | A[0][3]=b3 )
 *   A=( A[1][0]=b1 | A[1][1]=a2 | 0          | 0          )
 *     ( A[2][0]=b2 | 0          | A[2][2]=a3 | 0          )
 *     ( A[3][0]=b3 | 0          | 0          | A[3][3]=a4 )
 *
 *   G(w) = a0 + b0^2 / ( w - a1 - \sum_i ( b_i^2 / ( w - a_{i+1} ) ) )
 *
 *
 *    3) Double Tri-diagonal matrix representation
 *
 *     ( A[0][0]=avc0| A[0][1]=bv0 | 0           | 0           | A[0][4]=bc0 | 0           | 0           )
 *   A=( A[1][0]=bv0 | A[1][1]=av1 | A[1][2]=bv1 | 0           | 0           | 0           | 0           )
 *     ( 0           | A[2][1]=bv1 | A[2][2]=av2 | A[2][3]=bv2 | 0           | 0           | 0           )
 *     ( 0           | 0           | A[3][2]=bv2 | A[3][3]=av3 | 0           | 0           | 0           )
 *     ( A[4][0]=bc0 | 0           | 0           | 0           | A[4][4]=ac1 | A[4][5]=bc1 | 0           )
 *     ( 0           | 0           | 0           | 0           | A[5][4]=bc1 | A[5][5]=ac2 | A[5][6]=bc2 )
 *     ( 0           | 0           | 0           | 0           | 0           | A[6][5]=bc2 | A[6][6]=ac3 )
 *
 *   G(w) = a0 + b0^2 / ( w - avc0 - Gv(w) - Gc(w) )
 *   a0 and b0 stored seperately
 *   Gv(w) = bv0^2 / (w - av1 - bv1^2 / (w - av2 - bv2^2 / (w - av3 - bv3^2 / (w ... ))))
 *   Gc(w) = bc0^2 / (w - ac1 - bc1^2 / (w - ac2 - bc2^2 / (w - ac3 - bc3^2 / (w ... ))))
 */

/*
 *   G_M = G_A + G_B
 *
 *   In order to calculate G_M we define a matrix M' wich has as its left top matrix A and as its right bottom matrix B
 *   This matrix is tri-diagonalized with respect to a start-vector.
 *   The start-vector is given as {A.b[0],0,0,...,B.b[0],0,0,0,...}
 *
 *   input A and B defined, M not initiallized
 */
uint32_t SumGreensfunctionsTriDiagonalRepresentation(TriDiagonalMatrixType const A, TriDiagonalMatrixType const B, TriDiagonalMatrixType *const M)
{
  double vecA[3][A.N],vecB[3][B.N];
  double Atemp,Btemp;
  uint32_t i,j;
  M->N=A.N+B.N;
  if(InitTriDiagonalMatrix(M)){printf("InitTriDiagonalMatrix failed in SumGreensfunctionsTriDiagonalRepresentation\n");fflush(stdout);return 1;};
  M->a[0]=A.a[0]+B.a[0];
  M->b[0]=sqrt(Sqr(A.b[0])+Sqr(B.b[0]));
  for(i=0;i<3;i++)
  {
    for(j=0;j<A.N;j++)
    {
      vecA[i][j]=0;
    }
    for(j=0;j<B.N;j++)
    {
      vecB[i][j]=0;
    }
  }
  if(M->b[0]>10*DBL_MIN)
  {
    vecA[0][0]=A.b[0]/M->b[0];
    vecB[0][0]=B.b[0]/M->b[0];
  }
  for(i=1;i<M->N;i++)
  {
    RealTriDiagonalMatrixVector(A, vecA[(i-1)%3], vecA[i%3]);       /* vec[i] = A vec[i-1] */
    RealTriDiagonalMatrixVector(B, vecB[(i-1)%3], vecB[i%3]);
    RealVectorAdd(vecA[i%3],vecA[(i+1)%3],A.N,-(M->b[i-1]));        /* orthogonalization of vec[i] to vec[i-2], using that vec[i] at this point is H vec[i-1] */
    RealVectorAdd(vecB[i%3],vecB[(i+1)%3],B.N,-(M->b[i-1]));
    RealVectorDotVector(vecA[(i-1)%3],vecA[i%3],A.N,&Atemp);        /* a[i] = vec[i] . vec[i-1] */
    RealVectorDotVector(vecB[(i-1)%3],vecB[i%3],B.N,&Btemp);
    M->a[i]=Atemp+Btemp;
    RealVectorAdd(vecA[i%3],vecA[(i-1)%3],A.N,-M->a[i]);            /* vec[i] -= a[i] vec[i-1] */
    RealVectorAdd(vecB[i%3],vecB[(i-1)%3],B.N,-M->a[i]);
    RealVectorDotVector(vecA[i%3],vecA[i%3],A.N,&Atemp);            /* b[i] = sqrt( vec[i] . vec[i] ) */
    RealVectorDotVector(vecB[i%3],vecB[i%3],B.N,&Btemp);
    M->b[i]=sqrt(Atemp+Btemp);
    if(M->b[i]>100*DBL_MIN)                                         /* vec[i] *= 1/b[i] */
    {
      RealVectorMultiply(vecA[i%3],A.N,1.0/M->b[i]); 
      RealVectorMultiply(vecB[i%3],B.N,1.0/M->b[i]); 
    }
    else
    {
      M->N=i;
      RealVectorMultiply(vecA[i%3],A.N,0); 
      RealVectorMultiply(vecB[i%3],B.N,0);
    }
  }
  RealVectorTriDiagonalMatrixVector(A, vecA[(M->N-1)%3], &Atemp);   /* a[N] = vec[N-1] A vec[N-1] */
  RealVectorTriDiagonalMatrixVector(B, vecB[(M->N-1)%3], &Btemp);
  M->a[M->N]=Atemp+Btemp;
  return 0;
}

uint32_t SymmetrizeTridiagonalMatrix(TriDiagonalMatrixType *const M)
{
  TriDiagonalMatrixType Mplus,Mmin;
  uint32_t i;
  if(TriDiagonalMatrixCopy(&Mplus,*M)){printf("TriDiagonalMatrixCopy failed in SymmetrizeTridiagonalMatrix\n");fflush(stdout);return 1;};
  if(TriDiagonalMatrixCopy(&Mmin ,*M)){printf("TriDiagonalMatrixCopy failed in SymmetrizeTridiagonalMatrix\n");fflush(stdout);return 1;};
  strcpy(Mplus.Name,M->Name);
  strcpy(Mmin.Name,M->Name);
  FreeTriDiagonalMatrix(M);
  for(i=1;i<=Mmin.N;i++)
  {
    Mmin.a[i]*=-1;
  }
  SumGreensfunctionsTriDiagonalRepresentation(Mplus, Mmin, M);
  strcpy(M->Name,Mmin.Name);
  M->b[0]=Mplus.b[0];
  M->a[0]=Mplus.a[0];
  FreeTriDiagonalMatrix(&Mplus);
  FreeTriDiagonalMatrix(&Mmin);
  return 0;
}

/* M=A-B
 *
 */
uint32_t DifferenceGreensFunctionTriDiagonalRepresentation(TriDiagonalMatrixType const A, TriDiagonalMatrixType const B, TriDiagonalMatrixType *const M)
{
  ListOfPoles Ap,Bp;
  if(TriDiagonalMatrixToListOfPoles(A, &Ap)){printf("TriDiagonalMatrixToListOfPoles failed in DifferenceGreensFunctionTriDiagonalRepresentation\n");fflush(stdout);return 1;};
  if(TriDiagonalMatrixToListOfPoles(B, &Bp)){printf("TriDiagonalMatrixToListOfPoles failed in DifferenceGreensFunctionTriDiagonalRepresentation\n");fflush(stdout);return 1;};
  if(SubtractListOfPoles(&Ap, Bp))
  {
    if((Ap.N==1)&&(Ap.b[0]<0))
    {
      printf("SubtractListOfPoles returned a singel pole with weight %24.17E in DifferenceGreensFunctionTriDiagonalRepresentation\n",Ap.b[0]);
      printf("For \"%s\" - \"%s\"\n",A.Name,B.Name);
      printf("Set weight to zero\n");
      Ap.b[0]=0;
    }
    else
    {
      printf("SubtractListOfPoles failed in DifferenceGreensFunctionTriDiagonalRepresentation\n");
      fflush(stdout);
      return 1;
    }
  }
  if(ListOfPolesToTridiagonalMatrix(Ap, M, Ap.N)){printf("ListOfPolesToTridiagonalMatrix failed in DifferenceGreensFunctionTriDiagonalRepresentation\n");fflush(stdout);return 1;};
  M->a[0]=A.a[0]-B.a[0];
  return 0;
}

/* S = G0^{-1} - G^{-1}
 *
 */
uint32_t SelfEnergyTriDiagonalRepresentation(TriDiagonalMatrixType const G0, TriDiagonalMatrixType const G, TriDiagonalMatrixType *const Sigma)
{
  TriDiagonalMatrixType IG0,IG;
  IG0.N=G0.N-1;
  IG.N=G.N-1;
  IG0.a=&(G0.a[1]);
  IG0.b=&(G0.b[1]);
  IG.a=&(G.a[1]);
  IG.b=&(G.b[1]);
  strcpy(IG0.Name,G0.Name);
  strncat(IG0.Name, "^-1 (-1)", MAX(0,(int)((int)MAXStrLength-(int)strlen(IG0.Name))));
  strcpy(IG.Name,G.Name);
  strncat(IG.Name, "^-1 (-1)", MAX(0,(int)((int)MAXStrLength-(int)strlen(IG.Name))));
  return DifferenceGreensFunctionTriDiagonalRepresentation(IG, IG0, Sigma);
}

/*   GLocal(w) = G0( w - Sigma(w) )
 *             = 1/( w - H0 - Sigma(w) )
 *
 *   GLocal(w)^{-1} = w - H0 - Sigma(w)
 *                  = - ( -G0^{-1} + Sigma(w) )
 *
 */
uint32_t GLocalTriDiagonalRepresentation(TriDiagonalMatrixType const G0, TriDiagonalMatrixType const Sigma, TriDiagonalMatrixType *const G, uint32_t NTriDiag)
{
  TriDiagonalMatrixType G0Inv;
  uint32_t i;
  
  G0Inv.N=G0.N-1;
  G0Inv.a=&(G0.a[1]);
  G0Inv.b=&(G0.b[1]);
  if(SumGreensfunctionsTriDiagonalRepresentation(G0Inv, Sigma, G)){printf("SumGreensfunctionsTriDiagonalRepresentation failed in GLocalTriDiagonalRepresentation\n");return 1;};
  G->a[G->N]=G->a[G->N-1];
  for(i=G->N-1;i>0;i--)
  {
    G->a[i]=G->a[i-1];
    G->b[i]=G->b[i-1];
  }
  G->a[0]=G0.a[0];
  G->b[0]=G0.b[0];
  
  if(NTriDiag>(G->N))
  {
    printf("Warning GLocalTriDiagonalRepresentation called with NTriDiag=%i, but G0.N=%i and Sigma.N=%i\n",NTriDiag,G0.N,Sigma.N);
    printf("Returned Glocal with Glocal.N=%i\n",(G->N));
  }
  else
  {
    G->N=NTriDiag;
  }
  return 0;
}


/*   GLocal(w) = G0( w - Sigma(w) )
 *
 *             = 1/( w - H0 - Sigma(w) )
 *
 *   Sigma(w) = sa0 + (sb0^2)/(w - sa1 - sb1^2 (w-sa2 ...
 *   G0(w)    = 0   +    1   /(w - ga1 - gb1^2 (w-ga2 ...
 *
 */
uint32_t GLocalTriDiagonalRepresentationOld(TriDiagonalMatrixType const G0, TriDiagonalMatrixType const Sigma, TriDiagonalMatrixType *const G, uint32_t NTriDiag)
{
  uint32_t i,j;
  MatrixType T;
  SparseMatrixType GS;
  GS.Ni=G0.N*(Sigma.N+1);
  GS.Nj=G0.N*(Sigma.N+1);
  GS.QComplex=0;
  GS.NBitsKey=8*((uint32_t)(log((double)GS.Ni)/(8*log(2.0)))+1);
  if(InitSparseMatrix(&GS)){printf("InitSparseMatrix failed in GLocalTriDiagonalRepresentation\n");fflush(stdout);return 1;};
  GS.a0=G0.a[0];
  GS.b0=G0.b[0];
  for(i=0;i<G0.N;i++)
  {
    RealSparseMatrixAddElement(&GS,i*(Sigma.N+1),i*(Sigma.N+1),G0.a[i+1]-Sigma.a[0]);
    for(j=1;j<Sigma.N+1;j++)
    {
      RealSparseMatrixAddElement(&GS,i*(Sigma.N+1)+j  ,i*(Sigma.N+1)+j  , Sigma.a[j]  );
      RealSparseMatrixAddElement(&GS,i*(Sigma.N+1)+j-1,i*(Sigma.N+1)+j  , Sigma.b[j-1]);
      RealSparseMatrixAddElement(&GS,i*(Sigma.N+1)+j  ,i*(Sigma.N+1)+j-1, Sigma.b[j-1]);
    }
  }
  for(i=1;i<G0.N;i++)
  {
    RealSparseMatrixAddElement(&GS,(i-1)*(Sigma.N+1),(i  )*(Sigma.N+1),G0.b[i]);
    RealSparseMatrixAddElement(&GS,(i  )*(Sigma.N+1),(i-1)*(Sigma.N+1),G0.b[i]);
  }
  if(SparseMatrixToTriDiagonalMatrix(GS,&T,G,NTriDiag)){printf("SparseMatrixToTriDiagonalMatrix failed in GLocalTriDiagonalRepresentation\n");fflush(stdout);return 1;};
  FreeMatrix(&T);
  FreeSparseMatrix(&GS);
  return 0;
}

/*  GbathNew = (Sigma + Glocal^{-1}) ^{-1}
 *
 *  Assume Glocal.a[0]=0;
 *         Glocal.b[0]=1;
 */
uint32_t GBathNewTriDiagonalRepresentation(TriDiagonalMatrixType const Glocal, TriDiagonalMatrixType const Sigma, TriDiagonalMatrixType * const GBathNew)
{
  uint32_t i;
  TriDiagonalMatrixType GlocInv,GBathNewInv;
  GlocInv.N=Glocal.N-1;
  GlocInv.a=&(Glocal.a[1]);
  GlocInv.b=&(Glocal.b[1]);
  strcpy(GlocInv.Name,Glocal.Name);
  strncat(GlocInv.Name, "^-1 (-1)", MAX(0,(int)((int)MAXStrLength-(int)strlen(Glocal.Name))));
  DifferenceGreensFunctionTriDiagonalRepresentation(GlocInv, Sigma, &GBathNewInv);
  GBathNew->N=GBathNewInv.N+1;
  InitTriDiagonalMatrix(GBathNew);
  GBathNew->a[0]=0;
  GBathNew->b[0]=1;
  for(i=1;i<GBathNew->N;i++)
  {
    GBathNew->a[i]=GBathNewInv.a[i-1];
    GBathNew->b[i]=GBathNewInv.b[i-1];
  }
  GBathNew->a[GBathNew->N]=GBathNewInv.a[GBathNew->N-1];
  FreeTriDiagonalMatrix(&GBathNewInv);
  return 0;
}




/*  This routine returns a tri-diagonal matrix representing a greensfunction as a continued fraction in Gnew
 *  
 */
uint32_t TriDiagonalMatrixDysonScalar(TriDiagonalMatrixType const G, TriDiagonalMatrixType *const Gnew, MatrixType const Rho, uint32_t const NRho, TriDiagonalMatrixType *const Gp)
{
  SparseMatrixType M;
  MatrixType T;
  uint32_t i,j,k,l,NTriDiag,jstuck;
  M.Ni=G.N;
  for(i=0;i<NRho;i++)
  {
    M.Ni+=Gp[i].N*G.N;
  }
  M.Nj=M.Ni;
  M.QComplex=0;
  M.NBitsKey=16;
  InitSparseMatrix(&M);
  M.a0=G.a[0];
  M.b0=G.b[0];
  j=0;
  for(i=0;i<G.N-1;i++)
  {
    RealSparseMatrixAddElement(&M, j, j, G.a[i+1]);
    RealSparseMatrixAddElement(&M, j, j+1, G.b[i+1]);
    RealSparseMatrixAddElement(&M, j+1, j, G.b[i+1]);
    j++;
  }
  RealSparseMatrixAddElement(&M, j, j, G.a[G.N]);
  j++;
  jstuck=j;
  for(k=0;k<NRho;k++)
  {
    for(l=0;l<G.N;l++)
    {
      RealSparseMatrixAddElement(&M, l, j, Gp[k].b[0]*sqrt(fabs(Rho.VRe[l][k])));
      RealSparseMatrixAddElement(&M, j, l, Gp[k].b[0]*sqrt(fabs(Rho.VRe[l][k])));
      for(i=0;i<Gp[k].N-1;i++)
      {
        RealSparseMatrixAddElement(&M, j, j, Gp[k].a[i+1]);
        RealSparseMatrixAddElement(&M, j, j+1, Gp[k].b[i+1]);
        RealSparseMatrixAddElement(&M, j+1, j, Gp[k].b[i+1]);
        j++;
      }
      RealSparseMatrixAddElement(&M, j, j, Gp[k].a[Gp[k].N]);
      j++;
    }
  }
  /*
  printf("\n\n\n\n\n\n\n\n ******* \n\n\n\n\n");
  PrintSparseMatrix(M);
  */
  
  NTriDiag=G.N;
  for(i=0;i<NRho;i++)
  {
    NTriDiag=Max(NTriDiag,Gp[i].N);
  }
  if(SparseMatrixToTriDiagonalMatrix(M, &T, Gnew, NTriDiag)){printf("SparseMatrixToTriDiagonalMatrix failed in TriDiagonalMatrixDysonScalar\n");fflush(stdout);return 1;}
  FreeMatrix(&T);
  FreeSparseMatrix(&M);
  return 0;
}


/*  We start from the following Hamiltonian
 *
 *     ( O       V1      0       0       0         V2      0       0       0         ...       V3      0       0       0   )
 *     ( V1^+    H1+g1a1 g1b1    0       0         0       0       0       0         ...       0       0       0       0   )
 *     ( 0       g1b1    H1+g1a2 g1b2    0         0       0       0       0         ...       0       0       0       0   )
 *     ( 0       0       g1b2    H1+g1a3 ...       0       0       0       0         ...       0       0       0       0   )
 *     ( 0       0       0       ...     ...       0       0       0       0         ...       0       0       0       0   )
 *
 *     ( V2^+    0       0       0       0         H2+g2a1 g2b1    0       0         ...       0       0       0       0   )
 *     ( 0       0       0       0       0         g2b1    H2+g2a2 g2b2    0         ...       0       0       0       0   )
 *     ( 0       0       0       0       0         0       g2b2    H2+g2a3 ...       ...       0       0       0       0   )
 *     ( 0       0       0       0       0         0       0       ...     ...       ...       0       0       0       0   )
 *
 *     ( ...     ...     ...     ...     ...       ...     ...     ...     ...       ...       ...     ...     ...     ... )
 *
 *     ( V3^+    0       0       0       0         0       0       0       0         ...       H3+g3a1 g3b1    0       0   )
 *     ( 0       0       0       0       0         0       0       0       0         ...       g3b1    H3+g3a2 g3b2    0   )
 *     ( 0       0       0       0       0         0       0       0       0         ...       0       g3b2    H3+g3a3 ... )
 *     ( 0       0       0       0       0         0       0       0       0         ...       0       0       ...     ... )
 *
 *  and want to create a tridiagonal matrix on a Krylov basis starting from the first basis function
 *
 *  The matrices included are of special form:
 *  O     - tridiagonal
 *  Vi    - a dense matrix
 *  Hi    - diagonal matrix
 *  gi    - tridiagonal matrix, representing the coupled bath.
 *
 *  gnew  - resulting tridiagonal matrix of length gnew.N (set on entry)
 *
 *  H.Ni counts the number of blocks or the number of matrices V there are
 *  H.Nj is Max[V[i].Nj] for all i or larger
 *
 */
uint32_t TridiagonalizeImpurityCoupledToBathDyson(TriDiagonalMatrixType const O, MatrixType const *V, MatrixType const H, TriDiagonalMatrixType const *g, uint32_t N, TriDiagonalMatrixType * const gnew)
{
  SparseMatrixType Htot;
  MatrixType T;
  uint32_t i,j,k,l;
  Htot.Ni=O.N;
  for(i=0;i<N;i++)
  {
    Htot.Ni+=g[i].N * V[i].Nj;
  }
  Htot.Nj=Htot.Ni;
  if(Htot.Ni<16000)
  {
    Htot.NBitsKey=16;
  }
  else
  {
    Htot.NBitsKey=24;
  }
  Htot.QComplex = 0;
  Htot.a0=O.a[0];
  Htot.b0=O.b[0];
  for(i=0;i<N;i++)
  {
    if(V[i].QComplex)Htot.QComplex=1;
  }
  if(InitSparseMatrix(&Htot)){printf("Could not initiallize sparse matrix in TridiagonalizeImpurityCoupledToBathDyson\n");return 1;}
  SparseMatrixAddTridiagonalMatrix(&Htot, 0, 0, O);
  j=O.N;
  for(i=0;i<N;i++)
  {
    SparseMatrixAddMatrixTransposed(&Htot, 0, j, V[i], DBL_EPSILON);
    SparseMatrixAddMatrix(&Htot, j, 0, V[i], DBL_EPSILON);
    j+=g[i].N * V[i].Nj;    
  }
  l=O.N;
  if(Htot.QComplex)
  {
    for(i=0;i<N;i++)
    {
      for(j=0;j<g[i].N;j++)
      {
        for(k=0;k<V[i].Nj;k++)
        {
          ComplexSparseMatrixAddElement(&Htot, l, l, H.VRe[i][k]+g[i].a[j+1],0);
          if(j!=g[i].N-1)
          {
            ComplexSparseMatrixAddElement(&Htot, l, l+V[i].Nj, g[i].b[j+1],0);
            ComplexSparseMatrixAddElement(&Htot, l+V[i].Nj, l, g[i].b[j+1],0);
          }
          l++;
        }
      }
    }
  }
  else
  {
    for(i=0;i<N;i++)
    {
      for(j=0;j<g[i].N;j++)
      {
        for(k=0;k<V[i].Nj;k++)
        {
          RealSparseMatrixAddElement(&Htot, l, l, H.VRe[i][k]+g[i].a[j+1]);
          if(j!=g[i].N-1)
          {
            RealSparseMatrixAddElement(&Htot, l, l+V[i].Nj, g[i].b[j+1]);
            RealSparseMatrixAddElement(&Htot, l+V[i].Nj, l, g[i].b[j+1]); 
          }
          l++;
        }
      }
    }
  }
  SparseMatrixToTriDiagonalMatrix(Htot, &T, gnew, gnew->N);
  FreeMatrix(&T);
  FreeSparseMatrix(&Htot);
  return 0;
}

/*  We start from the following Hamiltonian
 *
 *     ( H_1     H_12^1         H_12^2          ...   N_12^n        )
 *     ( H_12^1  H_1^2 + H_2^1  0               0     0             )
 *     ( H_12^2  0              H_1^2 + H_2^2   0     0             )
 *     ( ...     0              0               ...   0             )
 *     ( H_12^n  0              0               0     H_1^n + H_2^n )
 *
 *     H_1 is provided on a finite basis in tridiagonal form (G1)
 *     on the same basis the matrix H_12^i H_12^i is given   (M12)
 *     The one particle green function for H2 with the coupling H_12^i is given in tridiagonal form (g2^i)
 *     The expectation value of H_12^i H_1 H_12^i devided by the expectation value of H_12^i H_12^i is given (Emu^i)
 *
 *     Calculated is a tridiaonal matrix representing { 1/(w-H) }_[1,1]
 */
uint32_t TridDiagonalMatrixDysonToMuchSimplified(TriDiagonalMatrixType const G1, MatrixType const *M12, TriDiagonalMatrixType const *g2, MatrixType const *Emu, int32_t const Ni, TriDiagonalMatrixType * const Gnew)
{
  return 0;
}