int twiddle(int *x, int *y, int *z, int *p);
void inittwiddle(int m, int n, int *p);
uint32_t InitRadMEList(RadMEList * M);
uint32_t StreamReadRadMEList(RadMEList * O, FILE * stream);
uint32_t FileReadRadMEList(RadMEList * O, const char *const flname);
uint32_t InitWavefunctionRandom(WaveFunctionType * const psi,uint32_t const NFermions,uint32_t const Nelectron);
uint32_t AtomicCI(const char *const flname, uint32_t Nelectron, uint32_t NOrb,uint32_t * const BasisLen, OperatorType *O);
