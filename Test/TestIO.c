static int callback(void *NotUsed, int argc, char **argv, char **azColName){
   int i;
   for(i=0; i<argc; i++){
     printf("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");
   }
    printf("\n");
   return 0;
  }

void TestSQL()
{
  sqlite3 *db;
  char *zErrMsg = 0;
  int  rc;
  char *sql;
  
  /* Open database */
  rc = sqlite3_open("newSQLDatabase", &db);
  if( rc ){
    fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db));
    exit(0);
  }else{
    fprintf(stdout, "Opened database successfully\n");
  }
  /* Create SQL statement */
  sql = "CREATE TABLE COMPANY("  \
  "ID INT PRIMARY KEY     NOT NULL," \
  "NAME           TEXT    NOT NULL," \
  "AGE            INT     NOT NULL," \
  "ADDRESS        CHAR(50)," \
  "SALARY         REAL );";
  
  /* Execute SQL statement */
  rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);
  if( rc != SQLITE_OK ){
    fprintf(stderr, "SQL error: %s\n", zErrMsg);
    sqlite3_free(zErrMsg);
  }else{
    fprintf(stdout, "Table created successfully\n");
  }
  
  /* Create SQL statement */
  sql = "INSERT INTO COMPANY (ID,NAME,AGE,ADDRESS,SALARY) "  \
  "VALUES (1, 'Paul', 32, 'California', 20000.00 ); " \
  "INSERT INTO COMPANY (ID,NAME,AGE,ADDRESS,SALARY) "  \
  "VALUES (2, 'Allen', 25, 'Texas', 15000.00 ); "     \
  "INSERT INTO COMPANY (ID,NAME,AGE,ADDRESS,SALARY)" \
  "VALUES (3, 'Teddy', 23, 'Norway', 20000.00 );" \
  "INSERT INTO COMPANY (ID,NAME,AGE,ADDRESS,SALARY)" \
  "VALUES (4, 'Mark', 25, 'Rich-Mond ', 65000.00 );";
  
  /* Execute SQL statement */
  rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);
  if( rc != SQLITE_OK ){
    fprintf(stderr, "SQL error: %s\n", zErrMsg);
    sqlite3_free(zErrMsg);
  }else{
    fprintf(stdout, "Records created successfully\n");
  }
  
  
  sqlite3_close(db);
  return ;
}


void TestBMP()
{
  int    i, j;
  MatrixType A;
  char   flname[256];
  A.Ni = 300;
  A.Nj = 500;
  A.QComplex = 0;
  InitMatrix(&A);
  for (i = 0; i < A.Ni; i++)
  {
    for (j = 0; j < A.Nj; j++)
    {
      A.VRe[i][j] = (double)(i * j) / ((double)(A.Ni * A.Nj));
    }
  }
  strcpy(flname, "test.bmp");
  FilePrintRealMatrixToBitmap(A, flname, &RealToRGBColorBW);
}

void TestWaveIO()
{
  WaveFunctionType psi, psiNew;
  uint8_t det[20];
  uint32_t i;
  psi.NBitsKey = 16;
  psi.NFermion = 10;
  WaveFunctionInit(&psi);
  strcpy(psi.Name, "Name of Wave Function");
  for (i = 1; i < 10; i++)
  {
    det[0] = (uint8_t) i;
    det[1] = (uint8_t) (i % 2);
    RealWaveFunctionAddElement(&psi, det, 0.5 / ((double)i));
  }
  FilePrintWaveFunctionPlain(psi, "temp", "w");
  psiNew.NBitsKey = 16;
  FileReadWaveFunctionPlain(&psiNew, "temp");
  PrintWaveFunctionPlain(psi);
  PrintWaveFunctionPlain(psiNew);
}
