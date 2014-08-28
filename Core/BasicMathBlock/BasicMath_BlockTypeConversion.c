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
        printf("%17.10E %17.10E | Ei %17.10E Ei+1 %17.10E A %17.10E\n",wmin,wplus,Emin+G[0].dE*iEmin,Emin+G[0].dE*(iEmin+1),A.VRe[j]);
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



uint32_t BlockListOfPolesToBlockTridiagonalMatrix(BlockListOfPoles const A, BlockTriDiagonalMatrixType *B)
{
  
  /* TO DO */
  
  return 0;
}