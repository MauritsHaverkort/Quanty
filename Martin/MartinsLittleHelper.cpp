/*
 * MartinsLittleHelper.cpp
 *
 *  Created on: Jul 10, 2018
 *      Author: brass
 */
#include <vector>

using std::vector;
// inclusion of namespace std leads to conflicts with function atomic in lua

void PrintMartix(double *M, int n)
{
	int i,j;
	for ( i = 0; i < n; i++)
	{
		for ( j = 0; j < n; j++)
		{
			printf("%f\t",M[i*n + j]);
		}
		printf("\n");
	}
	fflush(stdout);
}
void PrintMartix(Complex16 *M, int n)
{
	int i,j;
	for ( i = 0; i < n; i++)
	{
		for ( j = 0; j < n; j++)
		{
			printf("%f + I %f\t",M[i*n + j].Re,M[i*n + j].Im);
		}
		printf("\n");
	}
	fflush(stdout);
}

void PrintMartix2(double **M, int n)
{
	int i,j;
	for ( i = 0; i < n; i++)
	{
		for ( j = 0; j < n; j++)
		{
			printf("%f\t",M[i][j]);
		}
		printf("\n");
	}
	fflush(stdout);
}

void PrintArray(double* a, int n)
{
	int i;
	for ( i = 0; i < n; i++)
	{
		printf("%.15e\n", a[i]);
	}
}
void PrintArray(vector<int> a)
{
	int i;
	for ( i = 0; i < a.size(); i++)
	{
		printf("%d\n", a[i]);
	}
}
void PrintArray(Complex16* a, int n)
{
	int i;
	for ( i = 0; i < n; i++)
	{
		printf("%.15e + I %.15e\n", a[i].Re, a[i].Im);
	}
}

vector<double> ConvertArray(double *a, const int n)
{
	vector<double> res = vector<double>(n);
	for (int i = 0; i < n; i++)
	{
		res[i] = a[i];
	}
	return res;
}
vector<int> ConvertArray(uint16_t *a, const int n)
{
	vector<int> res = vector<int>(n);
	for (int i = 0; i < n; i++)
	{
		res[i] = (short) a[i];
	}
	return res;
}

vector<vector<double> > ConvertMatrix(double **m, const int n1, const int n2)
{
	vector<vector<double> > res = vector<vector<double> >(n1);
	for (int i = 0; i < n1; i++)
	{
		res[i] = vector<double>(n2);
		for (int j = 0; j < n2; j++)
		{
			res[i][j] = m[i][j];
		}
	}
	return res;
}

CompactMatrixType CalcSigmaContinuum(Complex16 z, double En, double *d, CompactMatrixType *Un)
{
	CompactMatrixType sigma, zdU;

	zdU.QComplex = 1;
	zdU.Ni = Un->Ni;
	zdU.Nj = Un->Nj;
	zdU.V = (Complex16 *) malloc(zdU.Ni * zdU.Nj * sizeof(Complex16));
	for(int i = 0; i < zdU.Ni; i++)
	{
		for(int j = 0; j < zdU.Nj; j++)
		{
			zdU.V[i*zdU.Nj + j] = Un->VRe[i*Un->Nj + j]/(z - En - d[i]);
		}
	}
	//PrintMartix(zdU.V,zdU.Ni);
	CompactMatrixMultiply('T','N',Un,&zdU,&sigma);
	free(zdU.V);
	return sigma;

}
Complex16 CalcPropagatorContinuum(Complex16 z0, Complex16 z1, double *d0, double *Econt, double *d1, vector<CompactMatrixType *> U, Complex16 *v)
{
	CompactMatrixType Sigma, SigmaAdd;
	//char UPLO = 'U';
	char NTRANS = 'N';
	uint32_t eins = 1;
	Complex16 aleph, beth, *SigV, retval;
	aleph.Im = 0.;
	aleph.Re = 1.;
	beth.Im = 0.;
	beth.Re = 0.;

	Sigma = CalcSigmaContinuum(z1,Econt[0],d1,U[0]);
	for (int n = 1; n < U.size(); n++)
	{
		SigmaAdd = CalcSigmaContinuum(z1,Econt[n],d1,U[n]);
		CompactMatrixAdd(&Sigma,&SigmaAdd);
		free(SigmaAdd.V);
	}
	for (int i = 0; i < Sigma.Ni; i++)
	{
		Sigma.V[i*Sigma.Nj + i] += d0[i] -z0;
	}
	CompactMatrixInverse(&Sigma);
	SigV = (Complex16 *) malloc(Sigma.Ni * sizeof(Complex16));
	zgemv(&NTRANS,&(Sigma.Ni),&(Sigma.Ni),&aleph,Sigma.V,&(Sigma.Ni),v,&eins,&beth,SigV,&eins);
	zdotc(&retval,&(Sigma.Ni),v,&eins,SigV,&eins);
	free(Sigma.V);
	free(SigV);
	return retval;
}

Complex16 CalcPropagatorContinuum(Complex16 z0, Complex16 z1, double *d0, double *Econt, vector<double *> d1, vector<CompactMatrixType *> U, Complex16 *v)
{
	CompactMatrixType Sigma, SigmaAdd;
	//char UPLO = 'U';
	char NTRANS = 'N';
	uint32_t eins = 1;
	Complex16 aleph, beth, *SigV, retval;
	aleph.Im = 0.;
	aleph.Re = 1.;
	beth.Im = 0.;
	beth.Re = 0.;

	Sigma = CalcSigmaContinuum(z1,Econt[0],d1[0],U[0]);
	for (int n = 1; n < U.size(); n++)
	{
		SigmaAdd = CalcSigmaContinuum(z1,Econt[n],d1[n],U[n]);
		CompactMatrixAdd(&Sigma,&SigmaAdd);
		free(SigmaAdd.V);
	}
	for (int i = 0; i < Sigma.Ni; i++)
	{
		Sigma.V[i*Sigma.Nj + i] += d0[i] -z0;
	}
	CompactMatrixInverse(&Sigma);
	SigV = (Complex16 *) malloc(Sigma.Ni * sizeof(Complex16));
	zgemv(&NTRANS,&(Sigma.Ni),&(Sigma.Ni),&aleph,Sigma.V,&(Sigma.Ni),v,&eins,&beth,SigV,&eins);
	zdotc(&retval,&(Sigma.Ni),v,&eins,SigV,&eins);
	free(Sigma.V);
	free(SigV);
	return retval;
}

vector<double> CalcSpectrumContinuum(vector<double> w, double g0, double g1, double Q, double *d0, double *Econt, double *d1, vector<CompactMatrixType *> U, Complex16 *v)
{
	vector<double> spec = vector<double>(w.size());
	Complex16 z0,z1,val;
	z0.Im = g0/2;
	z1.Im = g1/2;


	//#pragma omp parallel for
	for (int i = 0; i < w.size(); i++)
	{
		z0.Re = w[i];
		z1.Re = w[i];
		val = CalcPropagatorContinuum(z0,z1,d0,Econt,d1,U,v);
		spec[i] = val.Im;//*(Q-w[i])*(Q-w[i]); //TODO uncomment this
	}
	//PrintArray(&spec[0],w.size());
	return spec;
}

vector<double> CalcSpectrumContinuum(vector<double> w, double g0, double g1, double Q, double *d0, double *Econt,vector<double *> d1, vector<CompactMatrixType *> U, Complex16 *v)
{
	vector<double> spec = vector<double>(w.size());
	Complex16 z0,z1,val;
	z0.Im = g0/2;
	z1.Im = g1/2;


	//#pragma omp parallel for
	for (int i = 0; i < w.size(); i++)
	{
		z0.Re = w[i];
		z1.Re = w[i];
		val = CalcPropagatorContinuum(z0,z1,d0,Econt,d1,U,v);
		spec[i] = val.Im;//*(Q-w[i])*(Q-w[i]); //TODO uncomment this
	}
	//PrintArray(&spec[0],w.size());
	return spec;
}

void PrintMaxSlaterIntegral(SSlaterIntegral *integrals,const uint32_t NIntegrals)
{
	double m = 0;
	double mold = 0;
	int index = 0;
	for (int i = 0; i < NIntegrals; i++)
	{
		m = Max(m,abs(integrals[i].val[0]));
		if (m != mold)
		{
			index = i;
			mold = m;
		}
		printf("%e\n",integrals[i].val[0]);
	}
	printf("%e\t%d\n",m,index);
}

int GetIndexMaxSlaterIntegral(SSlaterIntegral *integrals,const uint32_t NIntegrals)
{
	double m = 0;
	double mold = 0;
	int index = 0;
	for (int i = 0; i < NIntegrals; i++)
	{
		m = Max(m,abs(integrals[i].val[0]));
		if (m != mold)
		{
			index = i;
			mold = m;
		}

	}
	return index;
}

vector<int> GetMaxSlaterIntegrals(SSlaterIntegral *integrals,const uint32_t NIntegrals, const double cutoff)
{

	vector<int> index;
	for (int i = 0; i < NIntegrals; i++)
	{

		if (abs(integrals[i].val[0]) >= cutoff)
		{
			index.push_back(i);
		}

	}
	return index;
}


















/*
vector<vector<Complex16> > ReadLuaComplexMatrix(lua_State *L, int i)
{
	MatrixType *M = luaL_checkMatrix(L,i);
}

*/
