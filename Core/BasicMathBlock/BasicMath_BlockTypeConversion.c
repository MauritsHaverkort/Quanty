/*  Transform a Block List of Poles representation of G to a spectral representation of G.
 *
 *  On entry:
 *    A contains the BlockLisOfPoles to transform
 *    G is a pointer to SpectraType wich is not allocated
 *    Emin and Emax set the min and max energy in the spectrum
 *    N the number of energies.
 *
 *    G = sum_i lim_{n->0^+} b_i^+ 1/(w - a_i + I n) b_i
 *    G = sum_i lim_{n->0^+} b^2 / (w - a_i + I n)
 *
 *    Implementation first the im part of lim_{n->0^+} 1/(w - a_i + I n) is added as a "delta" function moved to the energies E_i and E_{i+1} such that
 *      E_i < a and E_{i+1} > a
 *
 *      after this is done for all poles the real part of lim_{n->0^+} 1/(w - a_i + I n) is added.
 */
uint32_t BlockListOfPolesToSpectra(BlockListOfPoles const A, SpectraType *G, double Emin, double Emax, uint32_t N)
{
  int i,j,k;
  uint32_t iEmin,iEplus;
  double wmin,wplus;
  SpectraType *Gtemp;
  G[0].N=N;
  G[0].Emin=Emin;
  G[0].Emax=Emax;
  G[0].EminPole=Emin;
  G[0].EmaxPole=Emax;
  G[0].dE=(Emax-Emin)/((double)N);
  G[0].Gamma=0;
  for(i=0;i<A.BlockSize*A.BlockSize;i++)
  {
    SpectrumAllocate(G[0], &(G[i]));
  }
  if(A.QComplex)
  {
    Gtemp=malloc(A.BlockSize*A.BlockSize*sizeof(SpectraType));
    for(i=0;i<A.BlockSize*A.BlockSize;i++)
    {
      SpectrumAllocate(G[0], &(Gtemp[i]));
    }
    /* j loops over all blocks in A */
    for(j=0;j<(A.BlockSize*A.BlockSize+1)*A.N;j+=A.BlockSize*A.BlockSize+1)
    {
      /* the poles have energy A.V[j].Re */
      if(A.V[j].Re<=Emin)
      {
        iEmin=0;
        iEplus=0;
        wmin=1;
        wplus=0;
      }
      else if(A.V[j].Re>=Emax)
      {
        iEmin=N;
        iEplus=N;
        wmin=1;
        wplus=0;
      }
      else
      {
        iEmin=(uint32_t)((A.V[j].Re-Emin)/G[0].dE);
        iEplus=iEmin+1;
        wplus=((A.V[j].Re-Emin)/G[0].dE)-iEmin;
        wmin=1-wplus;
      }
      /* add the im part of the lorenzian to energy iEmin and iEplus */
      for(i=0;i<A.BlockSize*A.BlockSize;i++)
      {
        Gtemp[i].Im[iEmin]-=wmin*A.V[j+i+1].Re;
        Gtemp[i].Re[iEmin]+=wmin*A.V[j+i+1].Im;
        Gtemp[i].Im[iEplus]-=wplus*A.V[j+i+1].Re;
        Gtemp[i].Re[iEplus]+=wplus*A.V[j+i+1].Im;
      }
    }
  }
  else
  {
    Gtemp=G;
    for(j=0;j<(A.BlockSize*A.BlockSize+1)*A.N;j+=A.BlockSize*A.BlockSize+1)
    {
      /* the poles have energy A.V[j].Re */
      if(A.VRe[j]<=Emin)
      {
        iEmin=0;
        iEplus=0;
        wmin=1;
        wplus=0;
      }
      else if(A.VRe[j]>=Emax)
      {
        iEmin=N;
        iEplus=N;
        wmin=1;
        wplus=0;
      }
      else
      {
        iEmin=(uint32_t)((A.VRe[j]-Emin)/G[0].dE);
        iEplus=iEmin+1;
        wplus=((A.VRe[j]-Emin)/G[0].dE)-iEmin;
        wmin=1-wplus;
      }
      /* add the real part of the lorenzian to energy iEmin and iEplus */
      for(i=0;i<A.BlockSize*A.BlockSize;i++)
      {
        Gtemp[i].Im[iEmin]-=wmin*A.VRe[j+i+1];
        Gtemp[i].Im[iEplus]-=wplus*A.VRe[j+i+1];
      }
    }
  }
  /* the im part of the lorenzians is there, we now need to add the real part */
  /* note that in the complex case some lorenzians have a complex pre-factor */
  /* weight w_i at energy E_i contributes w_i / (w - E_i) to all energies w, except at w=E_i */
  #pragma omp parallel for
  for(i=0;i<A.BlockSize*A.BlockSize;i++)
  {
    for(j=0;j<=N;j++)
    {
      for(k=0;k<j;k++)
      {
        G[i].Re[k]-=Gtemp[i].Im[j]/(Pi*(k-j));
      }
      for(k=j+1;k<=N;k++)
      {
        G[i].Re[k]-=Gtemp[i].Im[j]/(Pi*(k-j));
      }
    }
  }
  if(A.QComplex)
  {
    #pragma omp parallel for
    for(i=0;i<A.BlockSize*A.BlockSize;i++)
    {
      for(j=0;j<=N;j++)
      {
        for(k=0;k<j;k++)
        {
          G[i].Im[k]+=Gtemp[i].Re[j]/(Pi*(k-j));
        }
        for(k=j+1;k<=N;k++)
        {
          G[i].Im[k]+=Gtemp[i].Re[j]/(Pi*(k-j));
        }
      }
      for(j=0;j<=N;j++)
      {
        G[i].Re[j]+=Gtemp[i].Re[j];
        G[i].Im[j]+=Gtemp[i].Im[j];
      }
      SpectrumFree(&(Gtemp[i]));
    }
    free(Gtemp);
  }
  #pragma omp parallel for
  for(i=0;i<A.BlockSize*A.BlockSize;i++)
  {
    for(j=0;j<=N;j++)
    {
      G[i].Im[j]*=Pi/G[i].dE;
      G[i].Re[j]*=Pi/G[i].dE;
    }
  }
  return 0;
}


/*  Given an tight binding Hamiltonian and a k-mesh defined by the number of divisions of a, b and c (Nk)
 *  returns the green function matrix of size A.NFermion by A.NFermion
 *          on an energy mesh with N+3 points. N+1 points from Emin to Emax and one point below Emin and one above Emax such that all weight and first moment is ok
 *  Calculates the eigen values and functions for each k
 *  for each k and eigen value i
 *  weight is added to the pole at iE and iE+1 such that Energy[iE] < val < Energy[iE+1]
 *  The matrix h_{jk} is fun[i,j] * fun[i,k]^* is added
 */
uint32_t TightBindingHamiltonianToBlockListOfPoles(TightBindingDef const A, int Nk[3], BlockListOfPoles * const G, double const Emin, double const Emax, uint32_t const NE)
{
  MatrixType H;
  CompactMatrixType T;
  double q[3], val[A.NFermion], rest, totmin=0, totmax=0, totnow;
  int    iq[3],ival,iE,i,j;
  double diq[3],dE=(Emax-Emin)/((double)NE),energy,Nktot;
  G->BlockSize=A.NFermion;
  G->N=NE+3;
  G->QComplex=1;
  strcpy(G->Name, A.Name);
  if(InitBlockListOfPoles(G))
  {
    printf("InitBlockListOfPoles failed in TightBindingHamiltonianToBlockListOfPoles\n");
    return 1;
  }
  energy=Emin-dE;
  for(i=0;i<G->N*(1+G->BlockSize*G->BlockSize);i+=(1+G->BlockSize*G->BlockSize))
  {
    G->V[i].Re=energy;
    energy+=dE;
  }
  for (iq[0] = 0; iq[0] < Nk[0]; iq[0]++)
  {
    for (iq[1] = 0; iq[1] < Nk[1]; iq[1]++)
    {
      for (iq[2] = 0; iq[2] < Nk[2]; iq[2]++)
      {
        diq[0] = iq[0] / ((double)(Nk[0]));
        diq[1] = iq[1] / ((double)(Nk[1]));
        diq[2] = iq[2] / ((double)(Nk[2]));
        q[0] = diq[0] * A.Cell.gUnitCell[0][0] + diq[1] * A.Cell.gUnitCell[1][0] + diq[2] * A.Cell.gUnitCell[2][0];
        q[1] = diq[0] * A.Cell.gUnitCell[0][1] + diq[1] * A.Cell.gUnitCell[1][1] + diq[2] * A.Cell.gUnitCell[2][1];
        q[2] = diq[0] * A.Cell.gUnitCell[0][2] + diq[1] * A.Cell.gUnitCell[1][2] + diq[2] * A.Cell.gUnitCell[2][2];
        CreateHamiltonianFromTightBindingDef(A, q, &H);
        EigensystemToCompact(H, &T, val);
        for(ival=0;ival<A.NFermion;ival++)
        {
          /* calculate matrix B for eigen-value ival */
          if(H.QComplex)
          {
            ComplexEigenFunctionToCharacter(&(T.V[ival*A.NFermion]), &H);
          }
          else
          {
            RealEigenFunctionToCharacter(&(T.VRe[ival*A.NFermion]), &H);
          }
          /* determine the position of the pole in the list of poles */
          rest = (val[ival] - Emin) / dE;
          iE = (int)(floor(rest));
          if (iE >= NE)
          {
            totnow=0;
            for(i=0;i<A.NFermion;i++)
            {
              totnow+=H.VRe[i][i];
            }
            if(totmax+totnow<10*DBL_MIN)
            {
              totmax+=10*DBL_MIN;
              totnow+=10*DBL_MIN;
            }
            G->V[(NE+2)*(1+G->BlockSize*G->BlockSize)].Re = (totmax * G->V[(NE+2)*(1+G->BlockSize*G->BlockSize)].Re + totnow * val[ival]) / (totmax + totnow);
            for(i=0;i<A.NFermion;i++)
            {
              for(j=0;j<A.NFermion;j++)
              {
                G->V[(NE+2)*(1+G->BlockSize*G->BlockSize)+1+j+i*A.NFermion].Re+=H.VRe[i][j];
                if(H.QComplex)
                {
                  G->V[(NE+2)*(1+G->BlockSize*G->BlockSize)+1+j+i*A.NFermion].Im+=H.VIm[i][j];
                }
              }
            }
            totmax+=totnow;
          }
          else if (iE < 0)
          {
            totnow=0;
            for(i=0;i<A.NFermion;i++)
            {
              totnow+=H.VRe[i][i];
            }
            if(totmin+totnow<10*DBL_MIN)
            {
              totmin+=10*DBL_MIN;
              totnow+=10*DBL_MIN;
            }
            G->V[0].Re = (totmin * G->V[0].Re + totnow * val[ival]) / (totmin + totnow);
            for(i=0;i<A.NFermion;i++)
            {
              for(j=0;j<A.NFermion;j++)
              {
                G->V[1+j+i*A.NFermion].Re+=H.VRe[i][j];
                if(H.QComplex)
                {
                  G->V[1+j+i*A.NFermion].Im+=H.VIm[i][j];
                }
              }
            }
            totmin+=totnow;
          }
          else
          {
            rest -= iE;  /* rest is a number between 0 and 1 pole at iE+1 gets 1-rest pole at iE+2 gest rest */
            for(i=0;i<A.NFermion;i++)
            {
              for(j=0;j<A.NFermion;j++)
              {
                G->V[(iE+1)*(1+G->BlockSize*G->BlockSize)+1+j+i*A.NFermion].Re+=(1-rest)*H.VRe[i][j];
                G->V[(iE+2)*(1+G->BlockSize*G->BlockSize)+1+j+i*A.NFermion].Re+=(  rest)*H.VRe[i][j];
                if(H.QComplex)
                {
                  G->V[(iE+1)*(1+G->BlockSize*G->BlockSize)+1+j+i*A.NFermion].Im+=(1-rest)*H.VIm[i][j];
                  G->V[(iE+2)*(1+G->BlockSize*G->BlockSize)+1+j+i*A.NFermion].Im+=(  rest)*H.VIm[i][j];
                }
              }
            }
          }
        }
      }
      FreeMatrix(&H);
      FreeCompactMatrix(&T);
    }
  }
  Nktot=Nk[0]*Nk[1]*Nk[2];
  for(i=0;i<G->N*(1+G->BlockSize*G->BlockSize);i+=(1+G->BlockSize*G->BlockSize))
  {
    for(j=1;j<1+G->BlockSize*G->BlockSize;j++)
    {
      G->V[i+j].Re/=Nktot;
      G->V[i+j].Im/=Nktot;
    }
  }
  ChopBlockListOfPoles(G, DBL_EPSILON);
  return 0;
}

/*  On input 
 *    matrix A is diagonal and real and of dimension psi->Nj. Values are stored as a vector of doubles
 *    matrix psi is of dimension psi->Ni * psi->Nj
 *  on output matrix T is of dimension psi->Ni*NBlock by psi->Nj such that T^{dag} A T is block tridiagonal of dimension Ni*NBlock with psi on element 1 to psi->Ni
 *    ABlock = T^dag A T
 *
 *  if psi is non-orthogonal 
 */
uint32_t DiagonalMatrixBlockTridiagonalize(double const *A, CompactMatrixType const *psi, CompactMatrixType *const T, BlockTriDiagonalMatrixType *ABlock, double NBlock, double SingularValue)
{
  int i,ii,j,k,oldindex,newindex,ipos,iipos,ia,ib,jcol,icol,allzero;
  double norm,dot,dotRe,dotIm;
  CompactMatrixType Tnow;
  ABlock->N=NBlock;
  ABlock->BlockSize=psi->Ni;
  ABlock->QComplex=psi->QComplex;
  InitBlockTridiagonalMatrix(ABlock);
  T->Ni=psi->Ni*NBlock;
  T->Nj=psi->Nj;
  T->QComplex=psi->QComplex;
  if(InitCompactMatrix(T))
  {
    printf("InitCompactMatrix failed in DiagonalMatrixBlockTridiagonalize\n");
    return 1;
  }
  Tnow.Name=T->Name;
  Tnow.QComplex=T->QComplex;
  if(T->QComplex)
  {
    memcpy(T->V, psi->V, psi->Nj*psi->Ni*sizeof(Complex16));
    Tnow.V=T->V;
    Tnow.Ni=psi->Ni;
    Tnow.Nj=psi->Nj;
    oldindex=0;
    newindex=psi->Ni*psi->Nj;
    for(k=1;k<ABlock->N;k++)
    {
      Tnow.Ni+=psi->Ni;
      for(i=0;i<psi->Ni;i++)
      {
        for(j=0;j<psi->Nj;j++)
        {
          Tnow.V[newindex].Re = A[j] * Tnow.V[oldindex].Re;
          Tnow.V[newindex].Im = A[j] * Tnow.V[oldindex].Im;
          oldindex++;
          newindex++;
        }
      }
      for(i=Tnow.Ni-psi->Ni;i<Tnow.Ni;i++)
      {
        ipos=i*Tnow.Nj;
        /* orthogonalize i to ii=0 to i */
        for(ii=0;ii<i;ii++)
        {
          iipos=ii*Tnow.Nj;
          dotRe=0;
          dotIm=0;
          for(j=ipos;j<ipos+Tnow.Nj;j++)
          {
            dotRe+=Tnow.V[iipos].Re*Tnow.V[j].Re + Tnow.V[iipos].Im*Tnow.V[j].Im;
            dotIm+=Tnow.V[iipos].Re*Tnow.V[j].Im - Tnow.V[iipos].Im*Tnow.V[j].Re;
            iipos++;
          }
          iipos=ii*Tnow.Nj;
          for(j=ipos;j<ipos+Tnow.Nj;j++)
          {
            Tnow.V[j].Re -=dotRe*Tnow.V[iipos].Re - dotIm*Tnow.V[iipos].Im;
            Tnow.V[j].Im -=dotRe*Tnow.V[iipos].Im + dotIm*Tnow.V[iipos].Re;
            iipos++;
          }
        }
        /* normalize vector i */
        norm=0;
        for(j=ipos;j<ipos+Tnow.Nj;j++)
        {
          norm+=Sqr(Tnow.V[j].Re)+Sqr(Tnow.V[j].Im);
        }
        if(norm<SingularValue)
        {
          norm=0;
        }
        else
        {
          norm=1.0/sqrt(norm);
        }
        for(j=ipos;j<ipos+Tnow.Nj;j++)
        {
          Tnow.V[j].Re*=norm;
          Tnow.V[j].Im*=norm;
        }
      }
    }
    iipos=ABlock->BlockSize*ABlock->BlockSize; /* position in tridiagonal matrix */
    for(ia=0;ia<ABlock->N;ia++)
    {
      for(i=0;i<ABlock->BlockSize;i++)
      {
        for(j=0;j<ABlock->BlockSize;j++)
        {
          icol=(i+ia*ABlock->BlockSize)*T->Nj;
          jcol=(j+ia*ABlock->BlockSize)*T->Nj;
          ABlock->aRe[iipos]=0;
          for(k=0;k<T->Nj;k++)
          {
            ABlock->a[iipos].Re+=T->V[jcol].Re*A[k]*T->V[icol].Re + T->V[jcol].Im*A[k]*T->V[icol].Im;
            ABlock->a[iipos].Im+=T->V[jcol].Re*A[k]*T->V[icol].Im - T->V[jcol].Im*A[k]*T->V[icol].Re;
            icol++;
            jcol++;
          }
          iipos++;
        }
      }
    }
    iipos=ABlock->BlockSize*ABlock->BlockSize; /* position in tridiagonal matrix */
    for(ib=0;ib<(ABlock->N-1);ib++)
    {
      for(i=0;i<ABlock->BlockSize;i++)
      {
        for(j=0;j<ABlock->BlockSize;j++)
        {
          icol=(i+(ib+1)*ABlock->BlockSize)*T->Nj;
          jcol=(j+ib*ABlock->BlockSize)*T->Nj;
          ABlock->bRe[iipos]=0;
          for(k=0;k<T->Nj;k++)
          {
            ABlock->b[iipos].Re+=T->V[jcol].Re*A[k]*T->V[icol].Re + T->V[jcol].Im*A[k]*T->V[icol].Im;
            ABlock->b[iipos].Im+=T->V[jcol].Re*A[k]*T->V[icol].Im - T->V[jcol].Im*A[k]*T->V[icol].Re;
            icol++;
            jcol++;
          }
          iipos++;
        }
      }
    }
  }
  else
  {
    /* calculate T */
    memcpy(T->VRe, psi->VRe, psi->Nj*psi->Ni*sizeof(double));
    Tnow.VRe=T->VRe;
    Tnow.Ni=psi->Ni;
    Tnow.Nj=psi->Nj;
    oldindex=0;
    newindex=psi->Ni*psi->Nj;
    for(k=1;k<ABlock->N;k++)
    {
      Tnow.Ni+=psi->Ni;
      for(i=0;i<psi->Ni;i++)
      {
        for(j=0;j<psi->Nj;j++)
        {
          Tnow.VRe[newindex] = A[j] * Tnow.VRe[oldindex];
          oldindex++;
          newindex++;
        }
      }
      allzero=1;
      for(i=Tnow.Ni-psi->Ni;i<Tnow.Ni;i++)
      {
        ipos=i*Tnow.Nj;
        /* orthogonalize i to j=0 to i */
        for(ii=0;ii<i;ii++)
        {
          iipos=ii*Tnow.Nj;
          dot=0;
          for(j=ipos;j<ipos+Tnow.Nj;j++)
          {
            dot+=Tnow.VRe[j]*Tnow.VRe[iipos];
            iipos++;
          }
          iipos=ii*Tnow.Nj;
          for(j=ipos;j<ipos+Tnow.Nj;j++)
          {
            Tnow.VRe[j] -=dot*Tnow.VRe[iipos];
            iipos++;
          }
        }
        /* normalize vector i */
        norm=0;
        for(j=ipos;j<ipos+Tnow.Nj;j++)
        {
          norm+=Sqr(Tnow.VRe[j]);
        }
        if(norm<SingularValue)
        {
          norm=0;
        }
        else
        {
          allzero=0;
          norm=1.0/sqrt(norm);
        }
        for(j=ipos;j<ipos+Tnow.Nj;j++)
        {
          Tnow.VRe[j]*=norm;
        }
      }
      if(allzero)
      {
        ABlock->N=k;
      }
    }
    /* calculate T^+ A T */
    iipos=ABlock->BlockSize*ABlock->BlockSize; /* position in tridiagonal matrix */
    for(ia=0;ia<ABlock->N;ia++)
    {
      for(i=0;i<ABlock->BlockSize;i++)
      {
        for(j=0;j<ABlock->BlockSize;j++)
        {
          icol=(i+ia*ABlock->BlockSize)*T->Nj;
          jcol=(j+ia*ABlock->BlockSize)*T->Nj;
          ABlock->aRe[iipos]=0;
          for(k=0;k<T->Nj;k++)
          {
            ABlock->aRe[iipos]+=T->VRe[jcol]*A[k]*T->VRe[icol];
            icol++;
            jcol++;
          }
          iipos++;
        }
      }
    }
    iipos=ABlock->BlockSize*ABlock->BlockSize; /* position in tridiagonal matrix */
    for(ib=0;ib<(ABlock->N-1);ib++)
    {
      for(i=0;i<ABlock->BlockSize;i++)
      {
        for(j=0;j<ABlock->BlockSize;j++)
        {
          icol=(i+(ib+1)*ABlock->BlockSize)*T->Nj;
          jcol=(j+ib*ABlock->BlockSize)*T->Nj;
          ABlock->bRe[iipos]=0;
          for(k=0;k<T->Nj;k++)
          {
            ABlock->bRe[iipos]+=T->VRe[jcol]*A[k]*T->VRe[icol];
            icol++;
            jcol++;
          }
          iipos++;
        }
      }
    }
  }
  return 0;
}

/*  on entry A contains the Green's function in a BlockListOfPoles format.
 *  on excit B contains the Green's function in a BlockTridiagonal format. (note that a0 and b0 are aditional prefactors)
 *  The returned matrix is always of the dimension NBLockMax by A.BlockSize, if couplings are small (SingularValue) they are set to zero.
 */
uint32_t BlockListOfPolesToBlockTridiagonalMatrix(BlockListOfPoles const A, BlockTriDiagonalMatrixType *B, int const NBlockMax, double const SingularValue)
{
  double *alpha;
  int k,i,ialpha,iN;
  char name[]="temp";
  CompactMatrixType betaisqr,betai,beta,T;
  
  betaisqr.QComplex=A.QComplex;
  betaisqr.Ni=A.BlockSize;
  betaisqr.Nj=A.BlockSize;
  betaisqr.Name=name;
  
  beta.QComplex=A.QComplex;
  beta.Ni=A.BlockSize*A.N;
  beta.Nj=A.BlockSize;
  InitCompactMatrix(&beta);
  
  alpha=malloc(A.N*A.BlockSize*sizeof(double));
  ialpha=0;
  iN=0;
  if(A.QComplex)
  {
    for(k=0;k<A.N;k++)
    {
      betaisqr.V = &(A.V[iN+1]);
      CompactMatrixSqrt(&betaisqr, &betai);
      memcpy(&(beta.V[iN-k]), betai.V, sizeof(Complex16)*A.BlockSize*A.BlockSize);
      FreeCompactMatrix(&betai);
      for(i=0;i<A.BlockSize;i++)
      {
        alpha[ialpha]=A.VRe[iN];
        ialpha++;
      }
      iN+=1+A.BlockSize*A.BlockSize;
    }
  }
  else /* A is real */
  {
    for(k=0;k<A.N;k++)
    {
      betaisqr.VRe = &(A.VRe[iN+1]);
      CompactMatrixSqrt(&betaisqr, &betai);
      memcpy(&(beta.VRe[iN-k]), betai.VRe, sizeof(double)*A.BlockSize*A.BlockSize);
      FreeCompactMatrix(&betai);
      for(i=0;i<A.BlockSize;i++)
      {
        alpha[ialpha]=A.VRe[iN];
        ialpha++;
      }
      iN+=1+A.BlockSize*A.BlockSize;
    }
  }
  CompactMatrixTranspose(&beta);
  DiagonalMatrixBlockTridiagonalize(alpha, &beta, &T, B, NBlockMax, SingularValue);
  B->mu=A.mu;
  free(alpha);
  FreeCompactMatrix(&beta);
  return 0;
}

/*  on entry A contains the Green's function in a BlockListOfPoles format.
 *  on excit B contains the Green's function in a BlockTridiagonal format. (note that a0 and b0 are aditional prefactors)
 *  The returned matrix is always of the dimension NBLockMax by A.BlockSize, if couplings are small (SingularValue) they are set to zero.
 *
 *  AT is a matrix of dimension Ni times Nj=A.BlockSize. 
 *  The actual tridiagonalization is done for the starting vectors AT * beta
 *
 */
uint32_t BlockListOfPolesToBlockTridiagonalMatrixReduceDimension(BlockListOfPoles const A, BlockTriDiagonalMatrixType *B, CompactMatrixType const *AT,  int const NBlockMax, double const SingularValue)
{
  double *alpha;
  int k,i,ialpha,iN;
  char name[]="temp";
  CompactMatrixType betaisqr,betai,beta,T,ATbeta;
  
  betaisqr.QComplex=A.QComplex;
  betaisqr.Ni=A.BlockSize;
  betaisqr.Nj=A.BlockSize;
  betaisqr.Name=name;
  
  beta.QComplex=A.QComplex;
  beta.Ni=A.BlockSize*A.N;
  beta.Nj=A.BlockSize;
  InitCompactMatrix(&beta);
  
  alpha=malloc(A.N*A.BlockSize*sizeof(double));
  ialpha=0;
  iN=0;
  if(A.QComplex)
  {
    for(k=0;k<A.N;k++)
    {
      betaisqr.V = &(A.V[iN+1]);
      CompactMatrixSqrt(&betaisqr, &betai);
      memcpy(&(beta.V[iN-k]), betai.V, sizeof(Complex16)*A.BlockSize*A.BlockSize);
      FreeCompactMatrix(&betai);
      for(i=0;i<A.BlockSize;i++)
      {
        alpha[ialpha]=A.VRe[iN];
        ialpha++;
      }
      iN+=1+A.BlockSize*A.BlockSize;
    }
  }
  else /* A is real */
  {
    for(k=0;k<A.N;k++)
    {
      betaisqr.VRe = &(A.VRe[iN+1]);
      CompactMatrixSqrt(&betaisqr, &betai);
      memcpy(&(beta.VRe[iN-k]), betai.VRe, sizeof(double)*A.BlockSize*A.BlockSize);
      FreeCompactMatrix(&betai);
      for(i=0;i<A.BlockSize;i++)
      {
        alpha[ialpha]=A.VRe[iN];
        ialpha++;
      }
      iN+=1+A.BlockSize*A.BlockSize;
    }
  }
  
  CompactMatrixMultiply('N', 'C', AT, &beta, &ATbeta);
  
  FreeCompactMatrix(&beta);
  
  DiagonalMatrixBlockTridiagonalize(alpha, &ATbeta, &T, B, NBlockMax, SingularValue);
  B->mu=A.mu;
  free(alpha);
  FreeCompactMatrix(&ATbeta);
  return 0;
}




/*
 *
 */
uint32_t BlockTriDiagonalMatrixToBlockListOfPoles(BlockTriDiagonalMatrixType const * A, BlockListOfPoles * const G, double const Emin, double const Emax, uint32_t const NE)
{
  int i,j,k,irho,kfun,iE;
  CompactMatrixType Afull, Afun, Rho;
  double energy, dE=(Emax-Emin)/((double)NE), rest, totmin=0, totmax=0, totnow;
  double *Aval;
  if(BlockTriDiagonalMatrixToCompactMatrix(A,&Afull))
  {
    printf("BlockTriDiagonalMatrixToCompactMatrix failed in BlockTriDiagonalMatrixToBlockListOfPoles\n");
    return 1;
  }
  Aval=malloc(A->BlockSize*A->N*sizeof(double));
  CompactMatrixEigensystem(&Afull, &Afun, Aval);
  
  G->BlockSize=A->BlockSize;
  G->N=NE+3;
  G->QComplex=1;
  strcpy(G->Name, A->Name);
  if(InitBlockListOfPoles(G))
  {
    printf("InitBlockListOfPoles failed in BlockTriDiagonalMatrixToBlockListOfPoles\n");
    return 1;
  }
  energy=Emin-dE;
  for(i=0;i<G->N*(1+G->BlockSize*G->BlockSize);i+=(1+G->BlockSize*G->BlockSize))
  {
    G->V[i].Re=energy;
    energy+=dE;
  }
  Rho.QComplex=Afun.QComplex;
  Rho.Ni=A->BlockSize;
  Rho.Nj=A->BlockSize;
  InitCompactMatrix(&Rho);
  kfun=0;
  for(k=0;k<Afun.Ni;k++)
  {
    /* create beta^2 */
    irho=0;
    for(i=0;i<A->BlockSize;i++)
    {
      for(j=0;j<A->BlockSize;j++)
      {
        Rho.VRe[irho] = Afun.VRe[kfun+i]*Afun.VRe[kfun+j];
        irho++;
      }
    }
    kfun+=Afun.Nj;
    /* add beta^2 to the correct position on the list of poles */
    
    
    /* determine the position of the pole in the list of poles */
    rest = (Aval[k] - Emin) / dE;
    iE = (int)(floor(rest));
    if (iE >= NE)
    {
      totnow=0;
      for(i=0;i<A->BlockSize;i++)
      {
        totnow+=Rho.VRe[i+i*A->BlockSize];
      }
      if(totmax+totnow<10*DBL_MIN)
      {
        totmax+=10*DBL_MIN;
        totnow+=10*DBL_MIN;
      }
      G->V[(NE+2)*(1+G->BlockSize*G->BlockSize)].Re = (totmax * G->V[(NE+2)*(1+G->BlockSize*G->BlockSize)].Re + totnow * Aval[k]) / (totmax + totnow);
      irho=0;
      for(i=0;i<A->BlockSize;i++)
      {
        for(j=0;j<A->BlockSize;j++)
        {
          if(Rho.QComplex)
          {
            G->V[(NE+2)*(1+G->BlockSize*G->BlockSize)+1+j+i*A->BlockSize].Re+=Rho.V[irho].Re;
            G->V[(NE+2)*(1+G->BlockSize*G->BlockSize)+1+j+i*A->BlockSize].Im+=Rho.V[irho].Im;
          }
          else
          {
            G->V[(NE+2)*(1+G->BlockSize*G->BlockSize)+1+j+i*A->BlockSize].Re+=Rho.VRe[irho];
          }
          irho++;
        }
      }
      totmax+=totnow;
    }
    else if (iE < 0)
    {
      totnow=0;
      irho=0;
      for(i=0;i<A->BlockSize;i++)
      {
        totnow+=Rho.VRe[irho];
        irho+=A->BlockSize+1;
      }
      if(totmin+totnow<10*DBL_MIN)
      {
        totmin+=10*DBL_MIN;
        totnow+=10*DBL_MIN;
      }
      G->V[0].Re = (totmin * G->V[0].Re + totnow * Aval[k]) / (totmin + totnow);
      irho=0;
      for(i=0;i<A->BlockSize;i++)
      {
        for(j=0;j<A->BlockSize;j++)
        {
          if(Rho.QComplex)
          {
            G->V[1+j+i*A->BlockSize].Re+=Rho.V[irho].Re;
            G->V[1+j+i*A->BlockSize].Im+=Rho.V[irho].Im;
          }
          else
          {
            G->V[1+j+i*A->BlockSize].Re+=Rho.VRe[irho];
          }
        }
      }
      totmin+=totnow;
    }
    else
    {
      rest -= iE;  /* rest is a number between 0 and 1 pole at iE+1 gets 1-rest pole at iE+2 gest rest */
      irho=0;
      for(i=0;i<A->BlockSize;i++)
      {
        for(j=0;j<A->BlockSize;j++)
        {
          if(Rho.QComplex)
          {
            G->V[(iE+1)*(1+G->BlockSize*G->BlockSize)+1+j+i*A->BlockSize].Re+=(1-rest)*Rho.V[irho].Re;
            G->V[(iE+2)*(1+G->BlockSize*G->BlockSize)+1+j+i*A->BlockSize].Re+=(  rest)*Rho.V[irho].Re;
            G->V[(iE+1)*(1+G->BlockSize*G->BlockSize)+1+j+i*A->BlockSize].Im+=(1-rest)*Rho.V[irho].Im;
            G->V[(iE+2)*(1+G->BlockSize*G->BlockSize)+1+j+i*A->BlockSize].Im+=(  rest)*Rho.V[irho].Im;
          }
          else
          {
            G->V[(iE+1)*(1+G->BlockSize*G->BlockSize)+1+j+i*A->BlockSize].Re+=(1-rest)*Rho.VRe[irho];
            G->V[(iE+2)*(1+G->BlockSize*G->BlockSize)+1+j+i*A->BlockSize].Re+=(  rest)*Rho.VRe[irho];
          }
          irho++;
        }
      }
    }
  }
  ChopBlockListOfPoles(G, DBL_EPSILON);
  FreeCompactMatrix(&Rho);
  FreeCompactMatrix(&Afull);
  FreeCompactMatrix(&Afun);
  return 0;
}




/*  Transform a BlockTriDiagonalMatrix representation of G to a spectral representation of G.
 *
 *  On entry:
 *    A contains the BlockTriDiagonal matrix to transform
 *    G is a pointer to SpectraType wich allocats A->BlockSize^2 spectra
 *    Emin and Emax set the min and max energy in the spectrum
 *    N the number of energies.
 *    Gamma sets the lorenzian broadening Gamma > 3*(Emax-Emin)/N
 *
 *    G = a0 + b0^+ (1/( w + i Gamma/2 - a1 - b1^+  ( w + i Gamma/2 - a2 - b2^+ ( ... ) b2 ) b1 )) b0
 *
 */
uint32_t BlockTriDiagonalMatrixToSpectra(BlockTriDiagonalMatrixType const * A, SpectraType *const G, double const Emin, double const Emax, uint32_t const N, double const Gamma)
{
  int i,iw,iblock;
  double omega;
  CompactMatrixType Gw,a,b;
  G[0].N=N;
  G[0].Emin=Emin;
  G[0].Emax=Emax;
  G[0].EminPole=Emin;
  G[0].EmaxPole=Emax;
  G[0].dE=(Emax-Emin)/((double)N);
  G[0].Gamma=Gamma;
  for(i=0;i<A->BlockSize*A->BlockSize;i++)
  {
    SpectrumAllocate(G[0], &(G[i]));
  }
  a.Ni=A->BlockSize;
  a.Nj=A->BlockSize;
  b.Ni=A->BlockSize;
  b.Nj=A->BlockSize;
  a.QComplex=A->QComplex;
  b.QComplex=A->QComplex;
  Gw.QComplex=1;
  Gw.Ni=A->BlockSize;
  Gw.Nj=A->BlockSize;
  InitCompactMatrix(&Gw);
  for(iw=0;iw<=N;iw++)
  {
    omega = Emin + iw * G[0].dE;
    for(i=0;i<A->BlockSize*A->BlockSize;i++)
    {
      Gw.V[i].Re=0;
      Gw.V[i].Im=0;
    }
    for(iblock=A->N;iblock>0;iblock--)
    {
      /* set a and b to point to right matrices */
      if(A->QComplex)
      {
        a.V=&(A->a[iblock    *A->BlockSize*A->BlockSize]);
        b.V=&(A->b[(iblock-1)*A->BlockSize*A->BlockSize]);
      }
      else
      {
        a.VRe=&(A->aRe[iblock    *A->BlockSize*A->BlockSize]);
        b.VRe=&(A->bRe[(iblock-1)*A->BlockSize*A->BlockSize]);
      }
      CompactMatrixAdd(&Gw, &a);
      for(i=0;i<A->BlockSize*A->BlockSize;i++)
      {
        Gw.V[i].Re*=-1;
        Gw.V[i].Im*=-1;
      }
      for(i=0;i<A->BlockSize*A->BlockSize;i+=A->BlockSize+1)
      {
        Gw.V[i].Re+=omega;
        Gw.V[i].Im+=Gamma/2.0;
      }
      CompactMatrixInverse(&Gw);
      TransformCompactMatrixInlineNNC(&Gw, &b); /* Gw <- b^dag Gw b */
    }
    if(A->QComplex)
    {
      a.V=A->a;
    }
    else
    {
      a.VRe=A->aRe;
    }
    CompactMatrixAdd(&Gw, &a);
    for(i=0;i<A->BlockSize*A->BlockSize;i++)
    {
      G[i].Re[iw]=Gw.V[i].Re;
      G[i].Im[iw]=Gw.V[i].Im;
    }
  }
  return 0;
}