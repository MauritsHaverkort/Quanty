int TightBindingImpurityUnfolding1Impurity3D(int Nsupercell, double Vpot, char *filename);
int TightBindingImpurityUnfolding1Impurity2D(int Nsupercell, double Vpot, char *filename);
int TightBindingImpurityUnfolding1Impurity1D(int Nsupercell, double Vpot, char *filename);
int TightBindingImpurityUnfoldingRandomImpurity3D(int Nsupercell, double Rho, double Vpot, char *filename,MatrixType * BWValue, MatrixType * RGBValue, SpectraType DOS[],MatrixType Fermi[]);
int TightBindingImpurityUnfoldingRandomImpurity2D(int Nsupercell, double Rho, double Vpot, char *filename,MatrixType * BWValue, MatrixType * RGBValue, SpectraType DOS[],MatrixType Fermi[]);
int TightBindingImpurityUnfoldingRandomImpurity1D(int Nsupercell, double Rho, double Vpot, char *filename,MatrixType * BWValue, MatrixType * RGBValue);
int TightBindingImpurityUnfoldingRandomImpuritySeries(int Dimension, int Nsupercell, double Rho, double Vpot,int repetitions, char *filename);
