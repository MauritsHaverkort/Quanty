void TestStuttgartNMTO()
{
  char testString[]="Hello World";
  printf("Length of String %lu\n",sizeof(testString));
  fflush(stdout);
  printf("string: %s\n",testString);
  StringReplace(testString, 'l', 'X');
  
  printf("string: %s\n",testString);
  
  
  StuttgartCTRL CTRLin;
  
  uint32_t i;
  char str[128];
  for(i=0;i<=NStuttgartCTRLCategories+1;i++)
  {
    IndexToStuttgartCategory(i, str);
    printf("index: %i %i - %s \n",i,StuttgartCategoryToIndex(str),str);
  }
  
  
  FileReadStuttgartCTRL(&CTRLin, "CTRL");
  
  printf("a\n");
}