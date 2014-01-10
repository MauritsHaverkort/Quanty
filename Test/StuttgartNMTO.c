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
  
  /* next
     TightBindingDefFromStuttgartCTRLHAMR
   
     Transform to bandstructure
   
     Transform to G0
   
     Transform to Gbath (plus sigma)
   
     Transoform to G (plus sigma)
   
     Transform to bandstructure + sigma
   
     Bandstructure + DOS
   
   */
  
  
  printf("a\n");
}