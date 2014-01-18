uint32_t BlockListOfPolesToSpectra(BlockListOfPoles const A, SpectraType *G, double Emin, double Emax, uint32_t N)
{
  uint32_t i;
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
  
  
  /* TO DO */
  
  return 0;
}

uint32_t TightBindingHamiltonianToBlockListOfPoles(TightBindingDef const HTB, double *Nk, BlockListOfPoles const * A)
{
  /* TO DO */
  
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