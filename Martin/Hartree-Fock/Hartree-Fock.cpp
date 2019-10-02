/*
 * Hartree-Fock.cpp
 *
 *  Created on: Sep 26, 2019
 *      Author: brass
 */

namespace HartreeFock
{

double OverlapMatrix(uint32_t i, uint32_t j, uint32_t n, const vector<double> &x)
{
	uint32_t ia = std::max(i,j);
	uint32_t ib = std::min(i,j) + n+1; // n+1 since the B-Splines start with n = 0 so n = polynomial degree
	vector<double> u, w;
	double integral = 0;
	double a,b;
	for(int k = ia; k < ib; k++)
	{
		a = x[k];
		b = x[k+1];
		if (a == b) continue;
		Quad::GaussianQuadratureWeights(n+1,u,w,a,b); //Bi*Bj is polynomial of deg 2n, hence  gauss quad 2(n+1) - 1 > 2n is exact
		for (int l = 0; l < n+1; l++)
		{
			integral += BSpline::Evaluate(i,n,x,u[l])*BSpline::Evaluate(j,n,x,u[l]) * w[l];
		}
	}
	return integral;
}
void CreateOverlapMatrix(const vector<double> &x, const uint32_t n, vector<double> &S)
{
	uint32_t Nspline = (uint32_t) x.size() - n - 1;

	S = vector<double>((Nspline-2)*(Nspline-2));
	for (int j = 0; j < Nspline-2; j++)
	{
		for ( int i = j; i <= std::min(Nspline-2-1,j+n); i++)  // n = #sub-diagonals
		{
			S[(i-j) + j*(Nspline-2)] = OverlapMatrix(i+1,j+1,n,x);
		}
	}
}

double LaplaceMatrix(uint32_t i, uint32_t j, uint32_t n, const vector<double> &x)
{
	double s1 = 0;
	double s2 = 0;
	double s3 = 0;
	double s4 = 0;

	if (x[i] != x[i+n])
	{
		if (x[j] != x[j+n])
		{
			s1 = OverlapMatrix(i,j,n-1,x)/(x[i+n] - x[i])/(x[j+n]-x[j]);
		}
		if (x[j+n+1] != x[j+1])
		{
			s2 = -OverlapMatrix(i,j+1,n-1,x)/(x[i+n] - x[i])/(x[j+n+1] - x[j+1]);
		}
	}
	if (x[i+n+1] != x[i+1])
	{
		if (x[j+n+1] != x[j+1])
		{
			s3 = OverlapMatrix(i+1,j+1,n-1,x)/(x[i+n+1] -x[i+1])/(x[j+n+1]-x[j+1]);
		}
		if (x[j+n] != x[j])
		{
			s4 = -OverlapMatrix(i+1,j,n-1,x)/(x[i+n+1]-x[i+1])/(x[j+n]-x[j]);
		}
	}
	return n*n*(s1 + s2 + s3 + s4);
}

void CreateLaplacian(const vector<double> &x, const uint32_t n, vector<double> &L)
{
	uint32_t Nspline = (uint32_t) x.size() - n - 1;

	L = vector<double>((Nspline-2)*(Nspline-2));
	for (int j = 0; j < Nspline-2; j++)
	{
		for ( int i = j; i <= std::min(Nspline-2-1,j+n); i++)  // n = #sub-diagonals
		{
			L[(i-j) + j*(Nspline-2)] = LaplaceMatrix(i+1,j+1,n,x);
		}
	}
}

double PotentialMatrix(uint32_t i, uint32_t j, uint32_t n, const vector<double> &x, double Z, int l)
{
	uint32_t ia = std::max(i,j);
	uint32_t ib = std::min(i,j) + n+1; // n+1 since the B-Splines start with n = 0 so n = polynomial degree
	vector<double> u, w;
	double integral = 0;
	double a,b;
	uint32_t Nknots = n + 3;
	for(int k = ia; k < ib; k++)
	{
		a = x[k];
		b = x[k+1];
		if (a == b) continue;
		Quad::GaussianQuadratureWeights(Nknots,u,w,a,b); //Bi*Bj is polynomial of deg 2n, hence  gauss quad 2(n+1) - 1 > 2n is exact
		for (int idx = 0; idx < Nknots; idx++)
		{
			integral += BSpline::Evaluate(i,n,x,u[idx])*BSpline::Evaluate(j,n,x,u[idx]) * w[idx] *(l*(l+1)/u[idx]/u[idx]/2 - Z/u[idx]);
		}
	}
	return integral;
}
void CreatePotentialMatrix(const vector<double> &x, const uint32_t n, vector<double> &V, double Z, int l)
{
	uint32_t Nspline = (uint32_t) x.size() - n - 1;

	V = vector<double>((Nspline-2)*(Nspline-2));
	for (int j = 0; j < Nspline-2; j++)
	{
		for ( int i = j; i <= std::min(Nspline-2-1,j+n); i++)  // n = #sub-diagonals
		{
			V[(i-j) + j*(Nspline-2)] = PotentialMatrix(i+1,j+1,n,x,Z,l);
		}
	}
}

void SolveCatOwnerEquation(vector<double> &H, vector<double> &S, uint32_t dim, uint32_t bandsize, vector<double> &eigVal,vector<double> &eigVec)
{
	if ((H.size() != S.size()) || ((uint32_t) H.size() != dim*dim))
		printf("error in SolveCatOwnerEquation: dimensions do not match: #H = %du, #S = %du, dim = %du",H.size(), S.size(), dim);

	vector<double> X = vector<double>(H.size());
	vector<double> e = vector<double>(dim-1);
	eigVal = vector<double>(dim);
	eigVec = vector<double>(H.size());
	char L = 'L';
	char V = 'V';
	int info = 0;
	double work[2*dim];
	uint32_t Lwork = 2*dim;
	dpbstf(&L,&dim,&bandsize,&S[0],&dim,&info);
	if (info) printf("info @ dpbstf = %d\n",info);
	dsbgst(&V,&L,&dim,&bandsize,&bandsize,&H[0],&dim,&S[0],&dim,&X[0],&dim,work,&info);
	if (info) printf("info @ dsbgst = %d\n",info);
	dsbtrd(&V,&L,&dim,&bandsize,&H[0],&dim,&eigVal[0],&e[0],&eigVec[0],&dim,work,&info);
	if (info) printf("info @ dsbtrd = %d\n",info);
	dorgtr(&L,&dim,&H[0],&dim,&eigVec[0],work,&Lwork,&info);
	if (info) printf("info @ dorgtr = %d\n",info);
	dsteqr(&V,&dim,&eigVal[0],&e[0],&eigVec[0],&dim,work,&info);
	if (info) printf("info @ dsteqr = %d\n",info);

}

void CreateGrid(unsigned int Nknots, int n, double R, int powGrid, vector<double> &x)
{
	x = vector<double>(Nknots);
	for (int i = 0; i < n+1; i++)
	{
		x[i] = 0.;
	}
	for (int i = n+1; i < Nknots-n; i++)
	{
		//b.x[i] = (i-n) * R/(Nknots-2*n-1);
		x[i] = std::pow(i-n,powGrid) * R / std::pow(Nknots-2*n-1,powGrid);
		//b.x[i] = R*(std::exp(i-n)-1)/(std::exp(Nknots-2*n-1)-1);
	}
	for (int i = Nknots-n; i < Nknots; i++)
	{
		x[i] = R;
	}
}

void CalcOffDiagSlaterInts(const vector<double> &x, vector<double> &rk1, vector<double> &rk2, int k, uint32_t n)
{
	uint32_t ks = n + 1;
	uint32_t nv = (uint32_t) x.size() - 1;
	uint32_t Nspline = (uint32_t) x.size() - n - 1;
	rk1 = vector<double>(nv *ks*ks);
	rk2 = vector<double>(nv *ks*ks);
	uint32_t Nweights = ks +3;
#pragma omp parallel for
	for (uint32_t iv = n; iv < nv; iv++)
	{
		vector<double> u,w;
		Quad::GaussianQuadratureWeights(Nweights,u,w,x[iv],x[iv+1]);
		uint32_t lowerIdx = iv > n ? iv-n : 0;
		for (uint32_t i = lowerIdx; i <= Min(iv,Nspline-1); i++)
		{
			for (uint32_t j = i; j <= Min(iv,Nspline-1); j++)
			{
				rk1[iv*ks*ks + (i-lowerIdx)*ks +j-lowerIdx] = 0;
				rk2[iv*ks*ks + (i-lowerIdx)*ks +j-lowerIdx] = 0;
				for (uint32_t idx = 0; idx < (uint32_t) u.size(); idx++)
				{
					rk1[iv*ks*ks + (i-lowerIdx)*ks +j-lowerIdx] += std::pow(u[idx],k) * w[idx] * BSpline::Evaluate(i,n,x,u[idx])* BSpline::Evaluate(j,n,x,u[idx]);
					rk2[iv*ks*ks + (i-lowerIdx)*ks +j-lowerIdx] += std::pow(u[idx],-k-1) * w[idx] * BSpline::Evaluate(i,n,x,u[idx])* BSpline::Evaluate(j,n,x,u[idx]);
				}
				rk1[iv*ks*ks + (j-lowerIdx)*ks +i-lowerIdx] = rk1[iv*ks*ks + (i-lowerIdx)*ks +j-lowerIdx];
				rk2[iv*ks*ks + (j-lowerIdx)*ks +i-lowerIdx] = rk2[iv*ks*ks + (i-lowerIdx)*ks +j-lowerIdx];
			}
		}
	}
}

void CalcDiagSlaterInts(const vector<double> &x, vector<double> &rkD, int k, uint32_t n)
{
	uint32_t ks = n + 1;
	uint32_t nv = (uint32_t) x.size() - 1;
	uint32_t Nspline = (uint32_t) x.size() - n - 1;
	uint32_t ks2 = ks*ks;
	uint32_t ks3 = ks2*ks;
	uint32_t ks4 = ks3*ks;
	uint32_t Nweights = ks +3;
	rkD = vector<double>(nv*ks4);
#pragma omp parallel for
	for (uint32_t iv = n; iv < nv; iv++)
	{
		vector<double> u1,w1, u2, w2;
		Quad::GaussianQuadratureWeights(Nweights,u1,w1,x[iv],x[iv+1]);
		uint32_t lowerIdx = iv > n ? iv-n : 0;
		for (uint32_t i1 = lowerIdx; i1 <= Min(iv,Nspline-1); i1++)
		{
			for (uint32_t j1 = lowerIdx; j1 <= Min(iv,Nspline-1); j1++) //TODO: Check symmetry
			{
				for (uint32_t i2 = i1; i2 <= Min(iv,Nspline-1); i2++)
				{
					for (uint32_t j2 = j1; j2 <= Min(iv,Nspline-1); j2++)
					{
						double Itot = 0;
						for (uint32_t idx1 = 0; idx1 < u1.size(); idx1++)
						{
							Quad::GaussianQuadratureWeights(Nweights,u2,w2,x[iv],u1[idx1]);
							double I2 = 0;
							for (uint32_t idx2 = 0; idx2 < u2.size(); idx2++)
							{
								I2 += w2[idx2] * std::pow(u2[idx2],k) * BSpline::Evaluate(j1,n,x,u2[idx2])* BSpline::Evaluate(j2,n,x,u2[idx2]);
							}
							Itot += I2 * w1[idx1] * std::pow(u1[idx1],-k-1) * BSpline::Evaluate(i1,n,x,u1[idx1])* BSpline::Evaluate(i2,n,x,u1[idx1]);
						}
						rkD[iv*ks4 +(i1-lowerIdx)*ks3 + (j1-lowerIdx)*ks2 + (i2-lowerIdx)*ks +(j2-lowerIdx)] = Itot;
						rkD[iv*ks4 +(i2-lowerIdx)*ks3 + (j1-lowerIdx)*ks2 + (i1-lowerIdx)*ks +(j2-lowerIdx)] = Itot;
						rkD[iv*ks4 +(i1-lowerIdx)*ks3 + (j2-lowerIdx)*ks2 + (i2-lowerIdx)*ks +(j1-lowerIdx)] = Itot;
						rkD[iv*ks4 +(i2-lowerIdx)*ks3 + (j2-lowerIdx)*ks2 + (i1-lowerIdx)*ks +(j1-lowerIdx)] = Itot;
//						rkD[iv*ks4 +(j1-lowerIdx)*ks3 + (i1-lowerIdx)*ks2 + (j2-lowerIdx)*ks +(i2-lowerIdx)] = Itot;
//						rkD[iv*ks4 +(j2-lowerIdx)*ks3 + (i1-lowerIdx)*ks2 + (j1-lowerIdx)*ks +(i2-lowerIdx)] = Itot;
//						rkD[iv*ks4 +(j1-lowerIdx)*ks3 + (i2-lowerIdx)*ks2 + (j2-lowerIdx)*ks +(i1-lowerIdx)] = Itot;
//						rkD[iv*ks4 +(j2-lowerIdx)*ks3 + (i2-lowerIdx)*ks2 + (j1-lowerIdx)*ks +(i1-lowerIdx)] = Itot;
					}
				}
			}
		}
	}
}

double AssembleSlaterInts(uint32_t i1, uint32_t j1,uint32_t i2, uint32_t j2, const vector<double> &rk1, const vector<double> &rk2, const vector<double> &rkD, uint32_t n)
{
	double val = 0;
	uint32_t ks = n + 1;
	uint32_t ks2 = ks*ks;
	uint32_t ks3 = ks2*ks;
	uint32_t ks4 = ks3*ks;
	uint32_t nv = (uint32_t) rk1.size() /ks /ks;
	uint32_t iv0i, iv1i,iv0j, iv1j;
	if (i1 < i2)
	{
		iv0i = i2;
		iv1i = i1+n;
	}
	else
	{
		iv0i = i1;
		iv1i = i2+n;
	}
	if (j1 < j2)
	{
		iv0j = j2;
		iv1j = j1+n;
	}
	else
	{
		iv0j = j1;
		iv1j = j2+n;
	}
	for (uint32_t ivi = iv0i; ivi <= Min(iv1i,nv); ivi++)
	{
		uint32_t lowerIdxi = ivi > n ? ivi-n : 0;
		for (uint32_t ivj = iv0j; ivj <= Min(iv1j,nv); ivj++)
		{
			uint32_t lowerIdxj = ivj > n ? ivj-n : 0;
			if (ivi == ivj)
			{
				val += rkD[ivi*ks4 +(i1-lowerIdxi)*ks3 + (j1-lowerIdxj)*ks2 + (i2-lowerIdxi)*ks +(j2-lowerIdxj)]
					+  rkD[ivi*ks4 +(j1-lowerIdxj)*ks3 + (i1-lowerIdxi)*ks2 + (j2-lowerIdxj)*ks +(i2-lowerIdxi)];
			}
			else if (ivi < ivj)
			{
				val += rk1[ivi*ks*ks + (i1-lowerIdxi)*ks +i2-lowerIdxi] * rk2[ivj*ks*ks + (j1-lowerIdxj)*ks +j2-lowerIdxj];
			}
			else
			{
				val += rk2[ivi*ks*ks + (i1-lowerIdxi)*ks +i2-lowerIdxi] * rk1[ivj*ks*ks + (j1-lowerIdxj)*ks +j2-lowerIdxj];
			}

		}
	}
	return val;
}

double AssembleSlaterInts(const vector<double> &Ci1,const vector<double> &Cj1,const vector<double> &Ci2,const vector<double> &Cj2,
						  const vector<double> &rk1, const vector<double> &rk2, const vector<double> &rkD, uint32_t n)
{
	double val = 0;
	uint32_t Nc = (uint32_t) Ci1.size();
	for (uint32_t i1 = 0; i1 < Nc; i1++)
		for (uint32_t j1 = 0; j1 < Nc; j1++)
			for (uint32_t i2 = 0; i2 < Nc; i2++)
				for (uint32_t j2 = 0; j2 < Nc; j2++)
					val += AssembleSlaterInts(i1,j1,i2,j2,rk1,rk2,rkD,n)*Ci1[i1]*Cj1[j1]*Ci2[i2]*Cj2[j2];
	return val;

}

void CreateFockMatrix(const uint32_t la, const vector<vector<double> > &rk1, const vector<vector<double> > &rk2, const vector<vector<double> > &rkD,
					  uint32_t n, const vector<vector<double> > &eigVec, uint32_t dim, const vector<vector<uint32_t> > &Nocc, vector<double> &F, vector<double> &FE)
{
	F = vector<double>(dim*dim);
	FE = vector<double>(dim*dim);
	for (uint32_t nu = 1; nu <= dim; nu++)
	{
		for (uint32_t mu = nu; mu <= dim; mu++)
		{
			uint32_t alpha0 = mu > n ? mu-n : 1;
			uint32_t alpha1 = mu+n < dim ? mu+n : dim;
			uint32_t beta0 = nu > n ? nu-n : 1;
			uint32_t beta1 = nu+n < dim ? nu+n : dim;
			for (uint32_t alpha = alpha0; alpha <= alpha1; alpha++)
			{
				for (uint32_t beta = beta0; beta <= beta1; beta++)
				{
					for (uint32_t lb = 0; lb < (uint32_t) Nocc.size(); lb++)
					{
						uint32_t k0 = Abs(la - lb);
						uint32_t k1 = la + lb;
						if ((k1 +k0) % 2) k0++;
						double Fk = 0;
						for (uint32_t k = k0; k <= k1; k +=2)
						{
							double cred = ThreeJSymbol(la,0,lb,0,k,0);
							Fk += AssembleSlaterInts(nu,alpha,beta,mu,rk1[k],rk2[k],rkD[k],n)*cred*cred;
						}
						for (uint32_t  j = 0; j < (uint32_t) Nocc[lb].size(); j++)
						{
							//F[(nu-1) * dim + (mu-1)]  += Fk* eigVec[lb][j * dim+(alpha-1)] * eigVec[lb][j * dim+(beta-1)]*(2*lb+1)*Nocc[lb][j]/2.;
							F[(nu-1) * dim + (mu-1)]  += Fk* eigVec[lb][j * dim+(alpha-1)] * eigVec[lb][j * dim+(beta-1)]*Nocc[lb][j]/2.;
							FE[(nu-1) * dim + (mu-1)] += Fk* eigVec[lb][j * dim+(alpha-1)] * eigVec[lb][j * dim+(beta-1)]*Nocc[lb][j]/2.;
							//TODO: check facotr  /2 since only spin parallel interaction
							//TODO: Lapack has Col Major and stores eigenvectors collumnwise, hence the j-th ev should be at [j*dim + i]
							// but is that really true or does Lapack do some strange stuff?
						}

					}
				}
			}
			F[(nu-1) + dim * (mu-1)]  = F[(nu-1) * dim + (mu-1)];
			FE[(nu-1) + dim * (mu-1)] = FE[(nu-1) * dim + (mu-1)];
		}
	}
}

void CreateHartreeMatrix(const uint32_t la, const vector<double> &rk1, const vector<double> &rk2, const vector<double> &rkD,
					  uint32_t n, const vector<vector<double> > &eigVec, uint32_t dim, const vector<vector<uint32_t> > &Nocc, vector<double> &H, vector<double> &HE)
{
	H = vector<double>(dim*dim);
	HE = vector<double>(dim*dim);
	for (uint32_t nu = 1; nu <= dim; nu++)
	{
		uint32_t mu1 = nu+n < dim ? nu+n : dim;
		for (uint32_t mu = nu; mu <= mu1; mu++)
		{
			for (uint32_t alpha = 1; alpha <= dim; alpha++)
			{
				uint32_t beta0 = alpha > n ? alpha-n : 1;
				uint32_t beta1 = alpha+n < dim ? alpha+n : dim;
				for (uint32_t beta = beta0; beta <= beta1; beta++) //TODO: alpha <= beta
				{
					for (uint32_t lb = 0; lb < (uint32_t) Nocc.size(); lb++)
					{

						double Fk = AssembleSlaterInts(nu,alpha,mu,beta,rk1,rk2,rkD,n);
						for (uint32_t  j = 0; j < (uint32_t) Nocc[lb].size(); j++)
						{
//							H[(nu-1) * dim + (mu-1)] += Fk* eigVec[lb][j * dim+(alpha-1)] * eigVec[lb][j * dim+(beta-1)]*(2*lb+1)*Nocc[lb][j];
							H[(nu-1) * dim + (mu-1)]  += Fk* eigVec[lb][j * dim+(alpha-1)] * eigVec[lb][j * dim+(beta-1)]*Nocc[lb][j];
							HE[(nu-1) * dim + (mu-1)] += Fk* eigVec[lb][j * dim+(alpha-1)] * eigVec[lb][j * dim+(beta-1)]*Nocc[lb][j];
							//TODO: Lapack has Col Major and stores eigenvectors collumnwise, hence the j-th ev should be at [j*dim + i]
							// but is that really true or does Lapack do some strange stuff?
						}

					}
				}
			}
			H[(nu-1) + dim * (mu-1)] = H[(nu-1) * dim + (mu-1)];
			HE[(nu-1) + dim * (mu-1)] = HE[(nu-1) * dim + (mu-1)];
		}
	}
}

double CalcGSEnergy(const vector<vector<uint32_t> > &Nocc, const vector<vector<double> > &eigVal, const vector<vector<double> > &eigVec, const vector<vector<double> > &H, uint32_t n)
{
	double e0 = 0;
	uint32_t dim = (uint32_t) eigVal[0].size();
	for (uint32_t  l = 0; l < (uint32_t) Nocc.size(); l++)
	{
		for (uint32_t i = 0; i < (uint32_t) Nocc[l].size(); i++)
		{
			e0 += Nocc[l][i] * eigVal[l][i];
			for (uint32_t mu = 0; mu < dim; mu++)
			{
				for (uint32_t nu = 0; nu < dim; nu++)
					e0 -= Nocc[l][i] * H[l][mu*dim + nu]*eigVec[l][i*dim +mu]*eigVec[l][i*dim +nu];
			}
		}
	}
	return e0;
}

void SCScheme(const vector<double> &x, uint32_t n, int Z,int Nelec, uint32_t NSCmax, double eps, vector<vector<double> > &eigVal,vector<vector<double> > &eigVec)
{
	uint32_t itype = 1;
	char jobz = 'V';
	char uplo = 'U';
	char trans ='N';
	char side = 'L';
	char diag = 'N';
	double one = 1;
	uint32_t lwork;
	int info;
	uint32_t lMax = (uint32_t) eigVec.size()-1;
	uint32_t kMax = lMax * 2;
	vector<vector<double> > rk1 = vector<vector<double> >(kMax+1);
	vector<vector<double> > rk2 = vector<vector<double> >(kMax+1);
	vector<vector<double> > rkD = vector<vector<double> >(kMax+1);
	vector<vector<uint32_t> > Nocc = vector<vector<uint32_t> >(lMax+1); //TODO: init
	BeginTiming("SI", &GlobalTime);
	for (uint32_t k = 0; k <= kMax; k++)
	{
		CalcOffDiagSlaterInts(x,rk1[k],rk2[k],k,n);
		CalcDiagSlaterInts(x,rkD[k],k,n);
	}
	EndTiming("SI", &GlobalTime);
	uint32_t Nsplines = (uint32_t) x.size() - n-1;
	uint32_t dim = Nsplines-2;
	double work[dim*dim];
	lwork = dim*dim; //TODO: find something that is not stupid
	vector<double> Lap, S, Sb;
	vector<vector<double> > V = vector<vector<double> >(lMax+1);
	CreateOverlapMatrix(x,n,Sb);
	S = vector<double>(dim*dim);
	vector<vector<double> > H = vector<vector<double> >(lMax+1);
	vector<vector<double> > HE = vector<vector<double> >(lMax+1);
	for (uint32_t  l = 0; l <= lMax; l++)
	{
		H[l] = vector<double>(dim*dim);
		HE[l] = vector<double>(dim*dim);
	}
	BeginTiming("1P", &GlobalTime);
	for (int j = 0; j < dim; j++)
	{
		for ( int i = j; i <= std::min(dim-1,j+n); i++)  // n = #sub-diagonals
		{
			S[i*dim+j] = Sb[(i-j) + j*dim];
			S[j*dim+i] = S[i*dim+j];
		}
	}
	dpotrf( &uplo, &dim, &(S[0]), &dim, &info );
	CreateLaplacian(x,n,Lap);
	for (uint32_t  l = 0; l <= lMax; l++)
		CreatePotentialMatrix(x,n,V[l],Z,l);
	EndTiming("1P", &GlobalTime);
	BeginTiming("SC", &GlobalTime);

	double gsEnergy0 = 0;
	double gsEnergy1 = 0;
	for (uint32_t idxSC = 1; idxSC <= NSCmax; idxSC++)
	{
		vector<double> K, F, KE, FE;
		for (uint32_t  l = 0; l <= lMax; l++)
		{
			for (int j = 0; j < dim; j++)
			{
				for ( int i = j; i <= std::min(dim-1,j+n); i++)  // n = #sub-diagonals
				{
					H[l][i*dim+j] += Lap[(i-j) + j*dim]/2. + V[l][(i-j) + j*dim];
					H[l][j*dim+i] = H[l][i*dim+j];
				}
			}
	//		printf("{ ");
	//		for (uint32_t i = 0; i < dim*dim; i++)
	//		{
	//			printf("%.16f, ",H[l][i]);
	//		}
	//		printf("}\n");
			dsygst( &itype, &uplo, &dim, &(H[l][0]), &dim, &(S[0]), &dim, &info );
			dsyev( &jobz, &uplo, &dim, &(H[l][0]), &dim, &(eigVal[l][0]),work,&lwork,&info);
			//dsygv(&itype,&jobz,&uplo,&dim,&(H[l][0]),&dim,&(S[0]),&dim,&(eigVal[l][0]),work,&lwork,&info); //TODO: be careful with EV update during loop
			if (info) printf("dsygv info = %d\n",info);
			dtrsm( &side, &uplo, &trans, &diag, &dim, &dim, &one, &(S[0]), &dim, &(H[l][0]), &dim );
		}
		if (idxSC == 1)
		{
			for (uint32_t  l = 0; l <= lMax; l++)
			{
				for (uint32_t i = 0; i < dim*dim; i++)
				{
					eigVec[l][i] = H[l][i];
				}
				CreateFockMatrix(l,rk1,rk2,rkD,n,eigVec,dim,Nocc,F,FE);
				CreateHartreeMatrix(l,rk1[0],rk2[0],rkD[0],n,eigVec,dim,Nocc,K,KE);

				for (uint32_t i = 0; i < dim*dim; i++)
				{
					H[l][i] = K[i] - F[i];
					HE[l][i] = (KE[i] - FE[i])/2.;
				}
			}
		}
		else
		{
			for (uint32_t  l = 0; l <= lMax; l++)
			{
				for (uint32_t i = 0; i < dim*dim; i++)
				{
					//eigVec[l][i] = (eigVec[l][i] + H[l][i])/2.;
					eigVec[l][i] = H[l][i];
				}
				CreateFockMatrix(l,rk1,rk2,rkD,n,eigVec,dim,Nocc,F,FE);
				CreateHartreeMatrix(l,rk1[0],rk2[0],rkD[0],n,eigVec,dim,Nocc,K,KE);

				for (uint32_t i = 0; i < dim*dim; i++)
				{
					H[l][i]  = K[i]  - F[i];
					HE[l][i] = (KE[i] - FE[i])/2.;
				}
			}
		}
		vector<uint32_t> curr = vector<uint32_t>(lMax+1);
		Nocc = vector<vector<uint32_t> >(lMax+1);
		int Nrest = Nelec;
		while (Nrest > 0)
		{
			uint32_t idxLow = 0;
			double e0 = eigVal[0][curr[0]];
			for (uint32_t  l = 1; l <= lMax; l++)
			{
				if (curr[l] >= (uint32_t) eigVal[l].size())continue;
				if (eigVal[l][curr[l]] < e0 - Abs(e0)*1E-10)
				{
					idxLow = l;
					e0 = eigVal[l][curr[l]];
				}
			}
			while ((int) Nocc[idxLow].size()-1 < (int) curr[idxLow])
				Nocc[idxLow].push_back(0);
			if (Nrest >= 2*(2*idxLow +1))
			{
				Nocc[idxLow][curr[idxLow]] = 2*(2*idxLow +1);
				Nrest -= 2*(2*idxLow +1);
			}
			else
			{
				Nocc[idxLow][curr[idxLow]] = Nrest;
				Nrest = 0;
			}
			curr[idxLow]++;

		}
		gsEnergy1 = CalcGSEnergy(Nocc,eigVal,eigVec,HE,n);
		printf("GS Energy = %f\n",gsEnergy1);
		if (Abs(1-gsEnergy0/gsEnergy1)<eps)
		{
			printf("Hartree-Fock Converged to precision of %.1E with ground-state energy = %f\n",eps,gsEnergy1);
			break;
		}
		gsEnergy0 = gsEnergy1;


	}
	EndTiming("SC", &GlobalTime);
}
}

static int LuaAtomicHartreeFock(lua_State *L)
{
    struct{
    	uint32_t n = 7;
    	uint32_t NSCmax = 20;
    	uint32_t lMax = 0;
    	double R = 20;
    	int m = 3;
    	double hMax = 0.375;
    	double eps = 1E-6;
    	double scale = 1;
  	} Options;
  	LuaOptionsType LuaCreateResonantSpectraOptions[]={
  	    {"BSpline_degree","uint32_t",&(Options.n)},
		{"NSCmax","uint32_t",&(Options.NSCmax)},
		{"lMax","uint32_t",&(Options.lMax)},
		{"Grid_R","double",&(Options.R)},
		{"Grid_m","int",&(Options.m)},
		{"Grid_scale","double",&(Options.scale)},
		{"Grid_hMax","double",&(Options.hMax)},
		{"Precision","double",&(Options.eps)},
  	    {NULL, NULL, NULL}};
  	int Z ,Nelec;
  	vector<double> x;
  	uint32_t len = (uint32_t)lua_gettop(L);
	lua_check(L,1,Z);
	lua_check(L,2,Nelec);
	if (Nelec <= 4)
	{
		Options.lMax = 0;
	}
	else if (Nelec <= 18)
	{
		Options.lMax = 1;
	}
	else if (Nelec <= 36)
	{
		Options.lMax = 2;
	}
	else
	{
		Options.lMax = 3;
	}
	Options.scale = Z;
	Options.R = 15 * Options.scale;
	if (len == 3) LuaCheckOptions(L, -1, LuaCreateResonantSpectraOptions);
	BSpline::GenerateFisherGrid(Options.m,Options.hMax,Options.R,Options.scale,Options.n,x);
	uint32_t Nsplines = (uint32_t) x.size() - Options.n-1;
	uint32_t dim = Nsplines-2;
	vector<vector<double> > eigVal = vector<vector<double> >(Options.lMax+1);
	vector<vector<double> > eigVec = vector<vector<double> >(Options.lMax+1);
	for (uint32_t l = 0; l <= Options.lMax; l++)
	{
		eigVal[l] = vector<double>(dim);
		eigVec[l] = vector<double>(dim*dim);
	}
	HartreeFock::SCScheme(x,Options.n,Z,Nelec,Options.NSCmax,Options.eps,eigVal,eigVec);
	lua_createtable(L,Options.lMax+1,0);
	vector<vector<BSpline::TBSpline*> > funs = vector<vector<BSpline::TBSpline*> >(Options.lMax+1);
	for (uint32_t l = 0; l <= Options.lMax; l++)
	{
		lua_createtable(L,0,dim);
		funs[l] = vector<BSpline::TBSpline*>(dim);
		for (uint32_t i = 0; i < dim; i++)
		{
			void* memory = lua_newuserdata(L, sizeof(BSpline::TBSpline));
			funs[l][i] = new(memory) BSpline::TBSpline();
			luaL_setmetatable(L, "InterpolatingFunction_Type");
			lua_rawseti(L,-2,i+1);
			funs[l][i]->knots = x;
			funs[l][i]->degree = Options.n;
			funs[l][i]->y.push_back(0);
			for (uint32_t mu = 1; mu <= dim; mu++)
				funs[l][i]->y.push_back(eigVec[l][i*dim + (mu-1)]);
			funs[l][i]->y.push_back(0);

		}
		lua_rawseti(L,-2,l+1);
	}


	return 1;
}


