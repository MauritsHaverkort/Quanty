uint32_t SplitStringToArray(char **array, char *str, char *del);
int ShellNameToL(char name);
double RoundToDigit(double x, int n);
void TransVecNorm(const double vec[3], double gUC[3][3], double vecnorm[3]);
uint32_t StreamReadDresdenFPLO(DresdenFPLO * O, FILE * stream);
uint32_t FileReadDresdenFPLO(DresdenFPLO * O, const char *const flname);
