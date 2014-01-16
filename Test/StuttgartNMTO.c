void TestStuttgartNMTO()
{
  char testString[]="Hello World";
  printf("Length of String %lu\n",sizeof(testString));
  fflush(stdout);
  printf("string: %s\n",testString);
  StringReplace(testString, 'l', 'X');
  
  printf("string: %s\n",testString);
  
  
  StuttgartCTRL CTRLin;
  StuttgartHAMR HAMRin;
  BandstructurePathType Path;
  CrystalStructureType Cell;
  TightBindingDef HTB;
  MatrixType BANDS, BWValue;
  double Elines[20];
  int Nk[3];
  SpectraType DOS;
  
  uint32_t i;
  char str[128];
  for(i=0;i<=NStuttgartCTRLCategories+1;i++)
  {
    IndexToStuttgartCategory(i, str);
    printf("index: %i %i - %s \n",i,StuttgartCategoryToIndex(str),str);
  }
  
  
  if(FileReadStuttgartCTRL(&CTRLin, "CTRL"))printf("EROR somewhere\n");
  
  FileReadStuttgartHAMR(&HAMRin, "HAMR");
  
  BandstructurePathFromStuttgartCTRL(CTRLin, &Path);
  PrintBandstructurePath(Path);
    
  CrystalStructureFromStuttgartCTRL(CTRLin, &Cell);
  
  PrintCrystalStructure(Cell);
  
  TightBindingDefFromStuttgartCTRLHAMR(CTRLin, HAMRin, &HTB);
  
  PrintTightBindingDef(HTB);
  
  TightBindingBANDS(HTB, Path, &BANDS);
  
  TightBindingBANDSToMatrix(-5, 10, 1000, BANDS, &BWValue);
  Elines[0]= -5;
  Elines[1]= 0;
  Elines[2]=5;
  Elines[3]=10;
  TightBindingAddAxesToMatrix(-5, 10, 4, Elines, Path, &BWValue);
  FilePrintRealMatrixToBitmap(BWValue, "Bands.bmp", &RealToRGBColorBWInverse);
  
  
  Nk[0]=80;
  Nk[1]=80;
  Nk[2]=80;
  
  DOS.Emin=-5;
  DOS.Emax=10;
  DOS.N=1500;
  DOS.Gamma=0.03;
  DOS.EminPole=DOS.Emin;
  DOS.EmaxPole=DOS.Emax;
  DOS.dE=(DOS.Emax-DOS.Emin)/DOS.N;
  SpectrumAllocate(DOS, &DOS);
  
  TightBindingDOS(HTB, Nk, &DOS);
  
  FilePrintSpectra(&DOS, 1, "DOS", "w");
  
  /* next
   
     Transform to G0
   
     Transform to Gbath (plus sigma)
   
     Transoform to G (plus sigma)
   
     Transform to bandstructure + sigma
   
     Bandstructure + DOS
   
   */
  
  
  printf("a\n");
}