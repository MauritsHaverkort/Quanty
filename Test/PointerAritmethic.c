void TestPointerArithmetic()
{
  int    a1[3] = { 0, 1, 2 };
  int    a2[2][3] = { {3, 4, 5}, {6, 7, 8} };
  int   *ip;                    /* pointer to int */
  int    (*ap)[3];              /* pointer to array [3] of int */
  int   *(bp[3]);



  ip = a1;
  printf("%d ", *ip);
  ip++;
  printf("%d\n", *ip);

  ap = a2;
  printf("%d %d\n", (*ap)[0], (*ap)[1]);
  printf("%p\n", ap);
  ap++;                         /* steps over entire (sub)array */
  printf("%d %d\n", (*ap)[0], (*ap)[1]);
  printf("%p\n", ap);

  printf("Sizeof (*ap)[3]: %lu\n", (unsigned long)sizeof ap);
  printf("Sizeof *(bp[3]): %lu\n", (unsigned long)sizeof bp);
  printf("Sizeof *ip: %lu\n", (unsigned long)sizeof ip);


}
