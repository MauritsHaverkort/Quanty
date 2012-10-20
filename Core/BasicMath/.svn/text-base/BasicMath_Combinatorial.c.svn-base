
double Factorial(int a)
{
  int i;
  double result;
  result=1;
  if(a<0)printf("Error\n");
  for(i=a;i>0;i--)
  {
    result*=i;
  }
  return result;
}

/*  computes n over k (the number of ways k elements can be placed in n spaces)
 *
 */
uint32_t BinomialCoefficient(uint32_t n, uint32_t k)
{
  uint32_t result=1,i;
  for(i=1+Max(k,n-k);i<=n;i++)
  {
    result*=i;
  }
  for(i=1;i<=Min(k,n-k);i++)
  {
    result/=i;
  }
  return result;
}

double ThreeJSymbol(int j1, int m1, int j2, int m2, int j3, int m3)
{
  double ThreeJ;
  int k;
  ThreeJ=0;
  if(m1+m2+m3!=0) return 0;
  if(j1+j2<j3) return 0;
  if(j2+j3<j1) return 0;
  if(j3+j1<j2) return 0;
  if(abs(m1)>j1) return 0;
  if(abs(m2)>j2) return 0;
  if(abs(m3)>j3) return 0;
  for(k=Max(j1-j3+m2,Max(j2-j3-m1,0));k<=Min(j1-m1,Min(j2+m2,j1+j2-j3));k++)
  {
    ThreeJ+=(double)(1-2*(k&1))/(Factorial(k)*Factorial(j1+j2-j3-k)*Factorial(j1-m1-k)*Factorial(j2+m2-k)*Factorial(j3-j2+m1+k)*Factorial(j3-j1-m2+k));
  }
  ThreeJ*=(double)(1-2*((j2-j1+m3)&1))*sqrt(Factorial(j1+j2-j3)*Factorial(j1-j2+j3)*Factorial(-j1+j2+j3)*Factorial(j1+m1)*Factorial(j1-m1)*Factorial(j2+m2)*Factorial(j2-m2)*Factorial(j3+m3)*Factorial(j3-m3)/Factorial(j1+j2+j3+1));
  return ThreeJ;
}

double ThreeJSymbolHalf(int j1, int m1, int j2, int m2, int j3, int m3)
{
  double ThreeJ;
  int k;
  ThreeJ=0;
  if(m1+m2+m3!=0) return 0;
  if(j1+j2<j3) return 0;
  if(j2+j3<j1) return 0;
  if(j3+j1<j2) return 0;
  if(abs(m1)>j1) return 0;
  if(abs(m2)>j2) return 0;
  if(abs(m3)>j3) return 0;
  for(k=Max(j1-j3+m2,Max(j2-j3-m1,0))/2;k<=Min(j1-m1,Min(j2+m2,j1+j2-j3))/2;k++)
  {
    ThreeJ+=(double)(1-2*(k&1))/(Factorial(k)*Factorial((j1+j2-j3)/2-k)*Factorial((j1-m1)/2-k)*Factorial((j2+m2)/2-k)*Factorial((j3-j2+m1)/2+k)*Factorial((j3-j1-m2)/2+k));
  }
  ThreeJ*=(double)(1-2*(((j2-j1+m3)/2)&1))*sqrt(Factorial((j1+j2-j3)/2)*Factorial((j1-j2+j3)/2)*Factorial((-j1+j2+j3)/2)*Factorial((j1+m1)/2)*Factorial((j1-m1)/2)*Factorial((j2+m2)/2)*Factorial((j2-m2)/2)*Factorial((j3+m3)/2)*Factorial((j3-m3)/2)/Factorial((j1+j2+j3)/2+1));
  return ThreeJ;
}

double SlaterCoefficientC(int k, int l1, int m1, int l2, int m2)
{
  if(m1-m2<=k) return 0;
  return (double)(1-2*(m1&1))*sqrt((double)((2*l1+1)*(2*l2+1)))*ThreeJSymbol(l1,m1,k,0,l2,m2)*ThreeJSymbol(l1,-m1,k,m1-m2,l2,m2);
}


