char LToShellName(int l);
int LMToOrbitalName(int l, int m, char **name);
std::string LPlusJToShellName(int lPlusJ);
int LPlusJAndJzToOrbitalName(int lPlusJ, int jz, char **name);
int EqualFloat(double a, double b, double epsilon);
uint32_t CrystalStructureFromDresdenFPLO(DresdenFPLO FPLO, CrystalStructureType *Cell);
uint32_t TightBindingDefFromDresdenFPLO(DresdenFPLO FPLO, TightBindingDef *HTB);
