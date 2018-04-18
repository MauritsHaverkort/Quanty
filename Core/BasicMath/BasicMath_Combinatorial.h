
int cmpdouble(const void *a,const void *b);
int cmpuint32t(const void *a,const void *b);
double SlaterCoefficientCjcoupling(int k,int l1,int twoj1,int twojz1,int l2,int twoj2,int twojz2);
double SlaterCoefficientC(int k,int l1,int m1,int l2,int m2);
double ThreeJSymbolHalf(int j1,int m1,int j2,int m2,int j3,int m3);
double ThreeJSymbol(int j1,int m1,int j2,int m2,int j3,int m3);
void RemovePerturbationIndex(uint32_t const Vi,uint32_t const IndexOld,uint32_t const NV,uint32_t const order,uint32_t *const IndexNew,int *const phase);
void AddPerturbationIndex(uint32_t const Vi,uint32_t const IndexOld,uint32_t const NV,uint32_t const order,uint32_t *const IndexNew,int *const phase);
void PerturbationOrderIndexToList(uint32_t const POIndex,uint32_t *const POList,uint32_t const NV,uint32_t const order);
void PerturbationOrderListToIndex(uint32_t *const POList,uint32_t *const POIndex,uint32_t const NV,uint32_t const order);
uint32_t PerturbationOptions(uint32_t NV,uint32_t order);
double BinomialCoefficientDouble(uint32_t n,uint32_t k);
uint32_t BinomialCoefficient(uint32_t n,uint32_t k);
double Factorial2(int a);
double Factorial(int a);
