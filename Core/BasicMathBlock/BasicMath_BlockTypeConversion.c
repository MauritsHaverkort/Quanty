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
 *    Implementation first the real part of lim_{n->0^+} 1/(w - a_i + I n) is added as a "delta" function moved to the energies E_i and E_{i+1} such that
 *      E_i < a and E_{i+1} > a
 *
 *      after this is done for all poles the imaginary part of lim_{n->0^+} 1/(w - a_i + I n) is added.
 */
uint32_t BlockListOfPolesToSpectra(BlockListOfPoles const A, SpectraType *G, double Emin, double Emax, uint32_t N)
{
  int i,j,k;
  uint32_t iEmin,iEplus;
  double wmin,wplus;
  SpectraType *Gtemp;
  G=malloc(A.BlockSize*A.BlockSize*sizeof(SpectraType));
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
        Gtemp[i].Im[iEmin]-=wmin*A.V[j+i].Re;
        Gtemp[i].Re[iEmin]+=wmin*A.V[j+i].Im;
        Gtemp[i].Im[iEplus]-=wplus*A.V[j+i].Re;
        Gtemp[i].Re[iEplus]+=wplus*A.V[j+i].Im;
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
        Gtemp[i].Im[iEmin]-=wmin*A.VRe[j+i];
        Gtemp[i].Im[iEplus]-=wplus*A.VRe[j+i];
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
        G[i].Re[k]-=Gtemp[i].Im[j]/(G[i].dE*(k-j));
      }
      for(k=j+1;k<=N;k++)
      {
        G[i].Re[k]-=Gtemp[i].Im[j]/(G[i].dE*(k-j));
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
          G[i].Im[k]+=Gtemp[i].Re[j]/(G[i].dE*(k-j));
        }
        for(k=j+1;k<=N;k++)
        {
          G[i].Im[k]+=Gtemp[i].Re[j]/(G[i].dE*(k-j));
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
  return 0;
}

uint32_t TightBindingHamiltonianToBlockListOfPoles(TightBindingDef const A, int Nk[3], BlockListOfPoles const * G)
{
  MatrixType H;
  CompactMatrixType T;
  double q[3], val[A.NFermion];
  int    iq[3];
  double diq[3];
  for (iq[0] = 0; iq[0] < Nk[0]; iq[0]++)
  {
    for (iq[1] = 0; iq[1] < Nk[1]; iq[1]++)
    {
      for (iq[2] = 0; iq[2] < Nk[2]; iq[2]++)
      {
        diq[0] = 1.0 * (iq[0] + 0.0) / ((double)(Nk[0]));
        diq[1] = 1.0 * (iq[1] + 0.0) / ((double)(Nk[1]));
        diq[2] = 1.0 * (iq[2] + 0.0) / ((double)(Nk[2]));
        q[0] = diq[0] * A.Cell.gUnitCell[0][0] + diq[1] * A.Cell.gUnitCell[1][0] + diq[2] * A.Cell.gUnitCell[2][0];
        q[1] = diq[0] * A.Cell.gUnitCell[0][1] + diq[1] * A.Cell.gUnitCell[1][1] + diq[2] * A.Cell.gUnitCell[2][1];
        q[2] = diq[0] * A.Cell.gUnitCell[0][2] + diq[1] * A.Cell.gUnitCell[1][2] + diq[2] * A.Cell.gUnitCell[2][2];
        CreateHamiltonianFromTightBindingDef(A, q, &H);
        EigensystemToCompact(H, &T, val);
        FreeMatrix(&H);
        
        
        /* to do */
        
        
        FreeCompactMatrix(&T);
      }
    }
  }
  return 0;
}





uint32_t BlockListOfPolesToBlockAndersonMatrix(BlockListOfPoles const A)
{
  
  /* TO DO */
  
  return 0;
}




uint32_t BlockListOfPolesToBlockTridiagonalMatrix(BlockListOfPoles const A)
{
  
  /* TO DO */
  
  return 0;
}