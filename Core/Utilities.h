int BitMaskToIndex(const char*BitMask, int *NF, int **Index);
int IndexToBitMask(int NF, int *Index, int len, char**BitMask);
int BitMaskCombine(const char *BitMask1, const char *BitMask2, char Operation, int *NF, char **BitMask);
int GetPhysicalConstant(const char *Name, double *Value);
int GetEnergyUnit(const char *Name, double *Value);
