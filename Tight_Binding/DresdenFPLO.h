char LToShellName(int l);
int LMToOrbitalName(int l, int m, char **name);
int EqualFloat(double a, double b);
uint32_t CrystalStructureFromDresdenFPLO(DresdenFPLO FPLO, CrystalStructureType *Cell);
uint32_t TightBindingDefFromDresdenFPLO(DresdenFPLO FPLO, TightBindingDef *HTB);
