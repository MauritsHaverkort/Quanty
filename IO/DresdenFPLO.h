
uint32_t FileReadDresdenFPLO(DresdenFPLO *O,const char *const flname);
uint32_t StreamReadDresdenFPLO(DresdenFPLO *O,FILE *stream);
void TransVecNorm(const double vec[3],const double gUC[3][3],double vecnorm[3]);
double RoundToDigit(double x,int n);
int ShellNameToL(char name);
uint32_t SplitStringToArray(char **array,char *str,char *del);
