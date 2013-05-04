/* data structures */


/*
 *   Reads a atom definition
 */
uint32_t StreamReadAtom(AtomType * O, FILE *stream)
{
  char comments[1000];
  int i;
  fscanf(stream,"%s %i %lf %lf %lf",(O->Name),&(O->AtomNumber),&(O->Position[0]),&(O->Position[1]),&(O->Position[2]));                    fgets(comments,1000,stream);
  fscanf(stream,"%i",&(O->Nshell));
  if((O->shelldegeneracy=malloc(O->Nshell*sizeof(int)))==NULL){printf("malloc failed for O->shelldegeneracy in StreamReadAtom\n");fflush(stdout);return 1;}
  for(i=0;i<O->Nshell;i++)
  {
    fscanf(stream,"%i",&(O->shelldegeneracy[i]));
  }
  fgets(comments,1000,stream);
  return 0;
}

/*
 *   Reads a crystal structure
 */
uint32_t StreamReadCrystalStructure(CrystalStructureType * O, FILE *stream)
{
  char comments[1000];
  uint32_t i;
  fscanf(stream,"%lf %lf %lf",&(O->rUnitCell[0][0]),&(O->rUnitCell[0][1]),&(O->rUnitCell[0][2]));                    fgets(comments,1000,stream);
  fscanf(stream,"%lf %lf %lf",&(O->rUnitCell[1][0]),&(O->rUnitCell[1][1]),&(O->rUnitCell[1][2]));                    fgets(comments,1000,stream);
  fscanf(stream,"%lf %lf %lf",&(O->rUnitCell[2][0]),&(O->rUnitCell[2][1]),&(O->rUnitCell[2][2]));                    fgets(comments,1000,stream);
  fscanf(stream,"%i",&(O->NAtoms));                                                                                  fgets(comments,1000,stream);
  O->Atom=malloc(O->NAtoms*sizeof(AtomType));
  for(i=0;i<O->NAtoms;i++)
  {
    StreamReadAtom(&(O->Atom[i]), stream);
  }
  ResiprocalLatticeVector(O->rUnitCell, O->gUnitCell);
  return 0;
}


/*
 *   Reads a crystal structure super cell
 */
uint32_t StreamReadCrystalStructureSuperCell(CrystalStructureSuperCellType * O, FILE *stream)
{
  char comments[1000];
  StreamReadCrystalStructure(&(O->Crystal),stream);
  fscanf(stream,"%i %i %i",&(O->rSuperCell[0][0]),&(O->rSuperCell[0][1]),&(O->rSuperCell[0][2]));                    fgets(comments,1000,stream);
  fscanf(stream,"%i %i %i",&(O->rSuperCell[1][0]),&(O->rSuperCell[1][1]),&(O->rSuperCell[1][2]));                    fgets(comments,1000,stream);
  fscanf(stream,"%i %i %i",&(O->rSuperCell[2][0]),&(O->rSuperCell[2][1]),&(O->rSuperCell[2][2]));                    fgets(comments,1000,stream);
  AddAditionalInformationToSuperCell(O);
  return 0;
}


/* input for specific programs */


/*
 *   Reads the settings for a Heisenberg model
 */
uint32_t StreamReadHeisenbergModelSettings(HeisenbergModelSettings * O, FILE *stream)
{
  char comments[1000];
  fscanf(stream,"%lf",&(O->J));                    fgets(comments,1000,stream);
  StreamReadCrystalStructureSuperCell(&(O->Cluster),stream);
  fscanf(stream,"%i",&(O->NPsi));                  fgets(comments,1000,stream);
  fscanf(stream,"%i",&(O->NE));                    fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->Emin));                 fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->Emax));                 fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->Gamma));                fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->zeroLancEnergy));       fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->zeroLancWavefunction)); fgets(comments,1000,stream);
  fscanf(stream,"%i",&(O->NTriDiag));              fgets(comments,1000,stream);
  fscanf(stream,"%i",&(O->NBitsKey));              fgets(comments,1000,stream);
  fscanf(stream,"%i",&(O->LocalFourierBasis));     fgets(comments,1000,stream);
  fscanf(stream,"%s",(O->SpecFile));               fgets(comments,1000,stream);
  return 0;
}


/*
 *   Reads the settings for a DimerizedSpinChain
 */
uint32_t StreamReadDimerizedSpinChainSettings(DimerizedSpinChainSettings * O, FILE *stream)
{
  char comments[1000];
  fscanf(stream,"%lf",&(O->JA));                   fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->JB));                   fgets(comments,1000,stream);
  fscanf(stream,"%i",&(O->Ndimer));                fgets(comments,1000,stream);
  fscanf(stream,"%i",&(O->NPsi));                  fgets(comments,1000,stream);
  fscanf(stream,"%i",&(O->NE));                    fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->Emin));                 fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->Emax));                 fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->Gamma));                fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->zeroLancEnergy));       fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->zeroLancWavefunction)); fgets(comments,1000,stream);
  fscanf(stream,"%i",&(O->NTriDiag));              fgets(comments,1000,stream);
  fscanf(stream,"%i",&(O->NBitsKey));              fgets(comments,1000,stream);
  fscanf(stream,"%s",(O->SpecFile));               fgets(comments,1000,stream);
  return 0;
}


/*
 *   Reads the settings for a tJ model
 */
uint32_t StreamReadtJModelSettings(tJModelSettings * O, FILE *stream)
{
  char comments[1000];
  fscanf(stream,"%lf",&(O->J));                    fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->t));                    fgets(comments,1000,stream);
  fscanf(stream,"%i",&(O->nElectron));             fgets(comments,1000,stream);
  StreamReadCrystalStructureSuperCell(&(O->Cluster),stream);
  fscanf(stream,"%i",&(O->NPsi));                  fgets(comments,1000,stream);
  fscanf(stream,"%i",&(O->NE));                    fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->Emin));                 fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->Emax));                 fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->Gamma));                fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->zeroLancEnergy));       fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->zeroLancWavefunction)); fgets(comments,1000,stream);
  fscanf(stream,"%i",&(O->NTriDiag));              fgets(comments,1000,stream);
  fscanf(stream,"%i",&(O->NBitsKey));              fgets(comments,1000,stream);
  fscanf(stream,"%i",&(O->LocalFourierBasis));     fgets(comments,1000,stream);
  fscanf(stream,"%s",(O->SpecFile));               fgets(comments,1000,stream);
  return 0;
}


/*
 *   Reads the settings for a 2 leg spin lader
 */
uint32_t StreamReadSpinLaderSettings(SpinLaderSettings * O, FILE *stream)
{
  char comments[1000];
  fscanf(stream,"%i",&(O->NSite));                 fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->Jleg));                 fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->Jrung));                fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->Jcyc));                 fgets(comments,1000,stream);
  fscanf(stream,"%i",&(O->NPsi));                  fgets(comments,1000,stream);
  fscanf(stream,"%i",&(O->NE));                    fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->Emin));                 fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->Emax));                 fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->Gamma));                fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->zeroLancEnergy));       fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->zeroLancWavefunction)); fgets(comments,1000,stream);
  fscanf(stream,"%i",&(O->NTriDiag));              fgets(comments,1000,stream);
  fscanf(stream,"%i",&(O->NBitsKey));              fgets(comments,1000,stream);
  fscanf(stream,"%s",(O->SpecFile));               fgets(comments,1000,stream);
  return 0;
}

/*
 *  Psi->NBitsKey needs to be set
 */
uint32_t StreamReadWaveFunctionPlain(WaveFunctionType * psi, FILE *stream)
{
  uint32_t i,j,psiN;
  double valRe,valIm;
  uint8_t *det;
  char buffer[1024+MAXStrLength],WaveName[MAXStrLength+1];
  fgets(WaveName,MAXStrLength+1,stream);
  WaveName[strlen(WaveName)-1]='\0';
  fgets(buffer,1023,stream);  sscanf(buffer,"%i",&(psi->QComplex));
  fgets(buffer,1023,stream);  sscanf(buffer,"%i",&(psiN));
  fgets(buffer,1023,stream);  sscanf(buffer,"%i",&(psi->NFermion));
  fgets(buffer,1023,stream);  sscanf(buffer,"%i",&(psi->NBoson));
  if(WaveFunctionInit(psi)){printf("WaveFunctionInit failed in StreamReadWaveFunctionPlain\n");fflush(stdout);return 1;};
  strcpy(psi->Name,WaveName);
  det=malloc(sizeof(uint8_t)*psi->NWord);
  if(det==NULL){WaveFunctionFree(psi);printf("malloc for det failed in StreamReadWaveFunctionPlain\n");fflush(stdout);return 1;};
  if(psi->QComplex)
  {
    for(i=1;i<=psiN;i++)
    {
      fscanf(stream,"%lf",&valRe);
      fscanf(stream,"%lf",&valIm);
      fgets(buffer,2,stream);
      for(j=0;j<psi->NWord;j++)
      {
        fgets(buffer,3,stream);
        sscanf(buffer,"%X",&(det[j]));
      }
      fgets(buffer,1023,stream);
      ComplexWaveFunctionAddElement(psi, det, valRe,valIm);
    }
  }
  else
  {
    for(i=1;i<=psiN;i++)
    {
      fscanf(stream,"%lf",&valRe);
      fgets(buffer,2,stream);
      for(j=0;j<psi->NWord;j++)
      {
        fgets(buffer,3,stream);
        sscanf(buffer,"%X",&(det[j]));
      }
      fgets(buffer,1023,stream);
      RealWaveFunctionAddElement(psi, det, valRe);
    }
  }
  return 0;
}


/*
 *  Psi->NBitsKey needs to be set
 */
uint32_t StreamReadWaveFunctionPlainBinary(WaveFunctionType * psi, FILE *stream)
{
  uint32_t i,psiN;
  double valRe,valIm;
  uint8_t *det;
  char buffer[1024+MAXStrLength],WaveName[MAXStrLength+1];
  fgets(WaveName,MAXStrLength+1,stream);
  WaveName[strlen(WaveName)-1]='\0';
  fgets(buffer,1023,stream);  sscanf(buffer,"%i",&(psi->QComplex));
  fgets(buffer,1023,stream);  sscanf(buffer,"%i",&(psiN));
  fgets(buffer,1023,stream);  sscanf(buffer,"%i",&(psi->NFermion));
  fgets(buffer,1023,stream);  sscanf(buffer,"%i",&(psi->NBoson));
  if(WaveFunctionInit(psi)){printf("WaveFunctionInit failed in StreamReadWaveFunctionPlain\n");fflush(stdout);return 1;};
  strcpy(psi->Name,WaveName);
  det=malloc(sizeof(uint8_t)*psi->NWord);
  if(det==NULL){WaveFunctionFree(psi);printf("malloc for det failed in StreamReadWaveFunctionPlain\n");fflush(stdout);return 1;};
  if(psi->QComplex)
  {
    for(i=1;i<=psiN;i++)
    {
      fread(&valRe, sizeof(double), 1, stream);
      fread(&valIm, sizeof(double), 1, stream);
      fread(det, sizeof(uint8_t), psi->NWord, stream);
      ComplexWaveFunctionAddElement(psi, det, valRe,valIm);
    }
  }
  else
  {
    for(i=1;i<=psiN;i++)
    {
      fread(&valRe, sizeof(double), 1, stream);
      fread(det, sizeof(uint8_t), psi->NWord, stream);
      RealWaveFunctionAddElement(psi, det, valRe);
    }
  }
  return 0;
}


uint32_t StreamReadMatrixType(MatrixType * O, FILE *stream)
{
  char comments[1000];
  int i,j;
  fscanf(stream,"%i %i %i",&(O->Ni),&(O->Nj),&(O->QComplex));    fgets(comments,1000,stream);
  if(InitMatrix(O)){printf("InitMatrix failed in StreamReadMatrixType\n");fflush(stdout);return 1;};
  for(i=0;i<O->Ni;i++)
  {
    for(j=0;j<O->Nj;j++)
    {
      fscanf(stream,"%lf",&(O->VRe[i][j]));
    }
    fgets(comments,1000,stream);
  }
  if(O->QComplex)
  {
    for(i=0;i<O->Ni;i++)
    {
      for(j=0;j<O->Nj;j++)
      {
        fscanf(stream,"%lf",&(O->VIm[i][j]));
      }
      fgets(comments,1000,stream);
    }
  }
  return 0;
}

uint32_t StreamReadAndersonMatrixType(AndersonMatrixType * O, FILE *stream)
{
  char comments[1000];
  int i;
  fscanf(stream,"%i",&(O->N)); fgets(comments,1000,stream);
  if(InitAndersonMatrix(O)){printf("InitAndersonMatrix failed in StreamReadAndersonMatrixType\n");fflush(stdout);return 1;};
  for(i=0;i<O->N;i++)
  {
    fscanf(stream,"%lf %lf",&(O->a[i]),&(O->b[i])); fgets(comments,1000,stream);
  }
  fscanf(stream,"%lf",&(O->a[O->N]));
  return 0;
}

uint32_t StreamReadPlainTriDiagonalMatrixType(TriDiagonalMatrixType * O, FILE *stream)
{
  double atemp,btemp;
  int i;
  fscanf(stream,"%lf %lf %i",&atemp, &btemp, &(O->N));
  if(ferror(stream)){printf("Error in StreamReadPlainTriDiagonalMatrixType error nr: %i\n",ferror(stream));perror("");fflush(stdout);return 1;};
  if(InitTriDiagonalMatrix(O)){printf("InitTriDiagonalMatrix failed in StreamReadPlainTriDiagonalMatrixType with O->N=%i\n",O->N);fflush(stdout);return 1;};
  O->a[0]=atemp;
  O->b[0]=btemp;
  for(i=1;i<O->N;i++)
  {
    fscanf(stream,"%lf %lf", &(O->a[i]), &(O->b[i]));
    if(ferror(stream)){printf("Error in StreamReadPlainTriDiagonalMatrixType error nr: %i\n",ferror(stream));perror("");return 1;};
    if(feof(stream)){printf("Error in StreamReadPlainTriDiagonalMatrixType eof reached at line %i of %i\n",i,O->N);return 1;};
  }
  fscanf(stream,"%lf", &(O->a[O->N]));
  return 0;
}


uint32_t StreamReadTightBindingDefHopping(TightBindingDefHopping * O, FILE *stream)
{
  char comments[1000];
  fscanf(stream,"%i %i %i %i",&(O->Atompos1),&(O->Shellpos1),&(O->Atompos2),&(O->Shellpos2));              fgets(comments,1000,stream);
  fscanf(stream,"%i %i %i",&(O->ConnectingVector[0]),&(O->ConnectingVector[1]),&(O->ConnectingVector[2])); fgets(comments,1000,stream);
  StreamReadMatrixType(&(O->Hop),stream);
  return 0;
}

uint32_t StreamReadTightBindingDef(TightBindingDef * O, FILE *stream)
{
  char comments[1000];
  int i;
  fscanf(stream,"%s",(O->Name));                   fgets(comments,1000,stream);
  StreamReadCrystalStructure(&(O->Cell),stream);
  fscanf(stream,"%i",&(O->nH));                    fgets(comments,1000,stream);
  if((O->H=malloc(O->nH*sizeof(TightBindingDefHopping)))==NULL){printf("malloc failed for TightBindingDefHopping in StreamReadTightBindingDef\n");fflush(stdout);return 1;};
  for(i=0;i<O->nH;i++)
  {
    StreamReadTightBindingDefHopping(&(O->H[i]),stream);
  }
  SetNFermionIndexTightBindingDef(O);
  return 0;
}


uint32_t StreamReadTightBindingDefNoCrystalStructure(TightBindingDef * O, FILE *stream)
{
  char comments[1000];
  int i;
  fscanf(stream,"%s",(O->Name));                   fgets(comments,1000,stream);
  fscanf(stream,"%i",&(O->nH));                    fgets(comments,1000,stream);
  if((O->H=malloc(O->nH*sizeof(TightBindingDefHopping)))==NULL){printf("malloc failed for TightBindingDefHopping in StreamReadTightBindingDef\n");fflush(stdout);return 1;};
  for(i=0;i<O->nH;i++)
  {
    StreamReadTightBindingDefHopping(&(O->H[i]),stream);
  }
  SetNFermionIndexTightBindingDef(O);
  return 0;
}

uint32_t StreamReadOperatorType(OperatorType * O, FILE *stream)
{
  char comments[1000];
  int i,j,NOperators;
  uint32_t length;
  int      oaddread;
  uint16_t oadd[256];
  double VRe,VIm;
  if(InitOperator(O)){printf("InitOperator failed in StreamReadOperatorFixedLengthType");fflush(stdout);return 1;}
  fscanf(stream,"%s",(O->Name));                                                            fgets(comments,1000,stream);
  fscanf(stream,"%i %i %i",&(O->NFermion),&(O->NBoson),&(NOperators));                      fgets(comments,1000,stream);
  O->NWord=O->NBoson+Min(O->NFermion,1+(O->NFermion-1)/8);
  for(i=0;i<NOperators;i++)
  {
    fscanf(stream,"%i",&length);
    if(length>O->MaxLength) InitOperatorFixedLength(O, length,8);
    for(j=length-1;j>=0;j--)
    {
      fscanf(stream,"%i",&oaddread);
      if(oaddread<0) oadd[j]=(uint16_t)(-1-oaddread);
      if(oaddread>0) oadd[j]=CrAnMask+(uint16_t)(oaddread-1);
    }
    fscanf(stream,"%lf %lf",&VRe,&VIm);
    fgets(comments,1000,stream);
    OperatorAddLader(O, oadd, length, VRe, VIm);
  }
  return 0;
}

/*
 *   Reads the settings for a Hubard model
 */
uint32_t StreamReadHubbardModelSettings(HubbardModelSettings * O, FILE *stream)
{
  char comments[1000];
  uint32_t iatom,ishell,ifermion;
  StreamReadCrystalStructureSuperCell(&(O->SuperCell),stream);
  CopyCrystalStructure(O->SuperCell.Crystal, &(O->TightBinding.Cell));
  StreamReadTightBindingDefNoCrystalStructure(&(O->TightBinding),stream);
  CreateTightBindingSuperCell(O->TightBinding, &(O->TightBindingSuperCell), O->SuperCell);
  StreamReadOperatorType(&(O->U),stream);
  CreateUSperCell(O->U, &(O->USuperCell), O->SuperCell);
  O->NElectron=malloc(O->SuperCell.Crystal.NAtoms*sizeof(uint32_t**));
  for(iatom=0;iatom<O->SuperCell.Crystal.NAtoms;iatom++)
  {
    O->NElectron[iatom]=malloc(O->SuperCell.Crystal.Atom[iatom].Nshell*sizeof(uint32_t*));
    for(ishell=0;ishell<O->SuperCell.Crystal.Atom[iatom].Nshell;ishell++)
    {
      O->NElectron[iatom][ishell]=malloc(O->SuperCell.Crystal.Atom[iatom].shelldegeneracy[ishell]*sizeof(uint32_t));
      for(ifermion=0;ifermion<O->SuperCell.Crystal.Atom[iatom].shelldegeneracy[ishell];ifermion++)
      {
        fscanf(stream,"%i",&(O->NElectron[iatom][ishell][ifermion]));
      }
    }
    fgets(comments,1000,stream);
  }
  fscanf(stream,"%i %i %i",&(O->Nq[0]),&(O->Nq[1]),&(O->Nq[2])); fgets(comments,1000,stream);
  fscanf(stream,"%i",&(O->NPsi));                  fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->DeltaEMax));            fgets(comments,1000,stream);
  fscanf(stream,"%i",&(O->NE));                    fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->Emin));                 fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->Emax));                 fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->Gamma));                fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->zeroLancEnergy));       fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->zeroLancWavefunction)); fgets(comments,1000,stream);
  fscanf(stream,"%i",&(O->NTriDiag));              fgets(comments,1000,stream);
  fscanf(stream,"%i",&(O->NBitsKey));              fgets(comments,1000,stream);
  fscanf(stream,"%s",(O->SpecFile));               fgets(comments,1000,stream);
  return 0;
}


/*
 *   Reads the settings for a DMFTXASSettings
 */
uint32_t StreamReadDMFTXASSettings(DMFTXASSettings * O, FILE *stream)
{
  char comments[1000];
  fscanf(stream,"%s",(O->BathFilename));           fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->U));                    fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->Q));                    fgets(comments,1000,stream);
  fscanf(stream,"%i",&(O->NBathMax));              fgets(comments,1000,stream);
  fscanf(stream,"%i",&(O->NE));                    fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->Emin));                 fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->Emax));                 fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->Gamma));                fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->zeroLancEnergy));       fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->zeroLancWavefunction)); fgets(comments,1000,stream);
  fscanf(stream,"%i",&(O->NTriDiag));              fgets(comments,1000,stream);
  fscanf(stream,"%i",&(O->NBitsKey));              fgets(comments,1000,stream);
  fscanf(stream,"%s",(O->SpecFile));               fgets(comments,1000,stream);
  return 0;
}


/*
 *   Reads the settings for a DMFTSettings
 */
uint32_t StreamReadDMFTSettings(DMFTSettings * O, FILE *stream)
{
  char comments[1000];
  uint32_t i;
  fscanf(stream,"%s",(O->G0FileName));               fgets(comments,1000,stream);
  fscanf(stream,"%s",(O->GRootFileName));            fgets(comments,1000,stream);
  fscanf(stream,"%s",(O->SpecRootFileName));         fgets(comments,1000,stream);
  fscanf(stream,"%s",(O->WaveFunctionRootFileName)); fgets(comments,1000,stream);
  fscanf(stream,"%i",&(O->iStart));                  fgets(comments,1000,stream);
  fscanf(stream,"%i",&(O->Ni));                      fgets(comments,1000,stream);
  fscanf(stream,"%i",&(O->Dyson));                   fgets(comments,1000,stream);
  fscanf(stream,"%i",&(O->NGCluster));               fgets(comments,1000,stream);
  fscanf(stream,"%i",&(O->NGBath));                  fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->Mix));                    fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->U));                      fgets(comments,1000,stream);
  fscanf(stream,"%i",&(O->NBath));                   fgets(comments,1000,stream);
  fscanf(stream,"%i",&(O->Spectra));
  O->NQ=1;
  if((O->Spectra)&12)
  {
    fscanf(stream,"%i",&(O->NQ));
    O->Q=malloc(O->NQ*sizeof(double));
    for(i=0;i<O->NQ;i++)
    {
      fscanf(stream,"%lf",&(O->Q[i]));
    }
    if(O->NQ==0)
    {
      O->Spectra&=!((uint32_t)12);
    }
  }                                                  fgets(comments,1000,stream);
  fscanf(stream,"%i",&(O->NE));                      fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->Emin));                   fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->Emax));                   fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->Gamma));                  fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->zeroLancEnergy));         fgets(comments,1000,stream);
  fscanf(stream,"%lf",&(O->zeroLancWavefunction));   fgets(comments,1000,stream);
  fscanf(stream,"%i",&(O->NTriDiag));                fgets(comments,1000,stream);
  fscanf(stream,"%i",&(O->NBitsKey));                fgets(comments,1000,stream);
  return 0;
}

