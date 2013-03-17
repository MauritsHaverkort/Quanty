uint32_t CreateWaveFunctionDMFT(WaveFunctionType *psi, uint32_t nBathV, uint32_t nBathC, uint32_t NBitsKey)
{
  uint32_t i;
  uint8_t *Determinant;
  psi->NBitsKey=NBitsKey;
  psi->QComplex=0;
  psi->NBoson=0;
  psi->NFermion=2+2*nBathC+2*nBathV;
  if(WaveFunctionInit(psi)){printf("WaveFunctionInit failed in CreateWaveFunctionDMFT\n");fflush(stdout);return 1;};
  sprintf(psi->Name,"DMFT %i %i",nBathV,nBathC);
  Determinant=calloc(psi->NWord,sizeof(uint8_t));
  OperatorCr(Determinant,1,psi->NFermion);
  for(i=2;i<nBathV+2;i++) OperatorCr(Determinant,i,psi->NFermion);
  for(i=2+nBathV+nBathC;i<2+2*nBathV+nBathC;i++) OperatorCr(Determinant,i,psi->NFermion);
  RealWaveFunctionAddElement(psi, Determinant, 1.0);
  free(Determinant);
  return 0;
}

uint32_t CreateWaveFunctionDMFTSecond(WaveFunctionType *psi, uint32_t nBathV, uint32_t nBathC, uint32_t NBitsKey)
{
  uint32_t i;
  uint8_t *Determinant;
  psi->NBitsKey=NBitsKey;
  psi->QComplex=0;
  psi->NBoson=0;
  psi->NFermion=4+2*nBathC+2*nBathV;
  if(WaveFunctionInit(psi)){printf("WaveFunctionInit failed in CreateWaveFunctionDMFT\n");fflush(stdout);return 1;};
  sprintf(psi->Name,"DMFT %i %i",nBathV,nBathC);
  Determinant=calloc(psi->NWord,sizeof(uint8_t));
  OperatorCr(Determinant,0,psi->NFermion);
  OperatorCr(Determinant,1,psi->NFermion);
  for(i=4;i<nBathV+4;i++) OperatorCr(Determinant,i,psi->NFermion);
  for(i=4+nBathV+nBathC;i<4+2*nBathV+nBathC;i++) OperatorCr(Determinant,i,psi->NFermion);
  RealWaveFunctionAddElement(psi, Determinant, 1.0);
  free(Determinant);
  return 0;
}

uint32_t CreatePESupOperatorDMFT(OperatorType *O, uint32_t NFermion)
{
  uint16_t oadd;
  InitOperator(O);
  O->NFermion=NFermion;
  O->NWord=Min(O->NFermion,1+(O->NFermion-1)/8);
  InitOperatorFixedLength(O, 1, 1);
  oadd=0;
  RealOperatorAddLader(O, &oadd, 1, 1.0);
  return 0;
}
uint32_t CreatePESdnOperatorDMFT(OperatorType *O, uint32_t NFermion)
{
  uint16_t oadd;
  InitOperator(O);
  O->NFermion=NFermion;
  O->NWord=Min(O->NFermion,1+(O->NFermion-1)/8);
  InitOperatorFixedLength(O, 1, 1);
  oadd=1;
  RealOperatorAddLader(O, &oadd, 1, 1.0);
  return 0;
}
uint32_t CreateIPESupOperatorDMFT(OperatorType *O, uint32_t NFermion)
{
  uint16_t oadd;
  InitOperator(O);
  O->NFermion=NFermion;
  O->NWord=Min(O->NFermion,1+(O->NFermion-1)/8);
  InitOperatorFixedLength(O, 1, 1);
  oadd=CrAnMask;
  RealOperatorAddLader(O, &oadd, 1, 1.0);
  return 0;
}
uint32_t CreateIPESdnOperatorDMFT(OperatorType *O, uint32_t NFermion)
{
  uint16_t oadd;
  InitOperator(O);
  O->NFermion=NFermion;
  O->NWord=Min(O->NFermion,1+(O->NFermion-1)/8);
  InitOperatorFixedLength(O, 1, 1);
  oadd=CrAnMask+1;
  RealOperatorAddLader(O, &oadd, 1, 1.0);
  return 0;
}
