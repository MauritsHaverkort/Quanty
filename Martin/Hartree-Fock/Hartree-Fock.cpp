/*
 * Hartree-Fock.cpp
 *
 *  Created on: Sep 26, 2019
 *      Author: brass
 */

#include <list>

namespace HartreeFock
{

// The matrix builders below are wrappers around the B-spline core
// (Core/BasicMath/Interpolation/BSplineMatrix), with the layout of before: element (i,j),
// i >= j, of the B-splines 1 .. dim (the first and last are left out), at (i-j) + j*dim in
// a vector of dim*dim, dim = Nspline-2. Martin Brass' original code is kept below as ...Old.
static void CreateRadialMatrix(const vector<double> &x, const uint32_t n, uint32_t nQ, uint32_t a, double Z, int l, bool potential, vector<double> &M)
{
	uint32_t dim = (uint32_t) x.size() - n - 1 - 2;
	double *rq, *wq;
	uint32_t nPoint;
	M = vector<double>(dim*dim);
	BSplineQuadratureSetup(&x[0], (uint32_t) x.size(), nQ, &rq, &wq, &nPoint);
	vector<double> fq = vector<double>(nPoint+1);
	for (uint32_t q = 0; q < nPoint; q++) fq[q] = l*(l+1)/rq[q]/rq[q]/2 - Z/rq[q];
	BSplineMatrixBanded(&x[0], (uint32_t) x.size(), n, rq, wq, nPoint, potential ? &fq[0] : NULL, a, a, 1, dim, &M[0], dim);
	free(rq);
	free(wq);
}

void CreateOverlapMatrix(const vector<double> &x, const uint32_t n, vector<double> &S)
{
	CreateRadialMatrix(x, n, n+1, 0, 0., 0, false, S);
}

// int B_i' B_j'
void CreateLaplacian(const vector<double> &x, const uint32_t n, vector<double> &L)
{
	CreateRadialMatrix(x, n, n, 1, 0., 0, false, L);
}

// int B_i (l(l+1)/(2r^2) - Z/r) B_j, with n+3 Gauss points per knot interval
void CreatePotentialMatrix(const vector<double> &x, const uint32_t n, vector<double> &V, double Z, int l)
{
	CreateRadialMatrix(x, n, n+3, 0, Z, l, true, V);
}

#if 0
double OverlapMatrixOld(uint32_t i, uint32_t j, uint32_t n, const vector<double> &x)
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
			integral += BSplineEvaluate(i,n,&x[0],(uint32_t)x.size(),u[l])*BSplineEvaluate(j,n,&x[0],(uint32_t)x.size(),u[l]) * w[l];
		}
	}
	return integral;
}

void CreateOverlapMatrixOld(const vector<double> &x, const uint32_t n, vector<double> &S)
{
	uint32_t Nspline = (uint32_t) x.size() - n - 1;

	S = vector<double>((Nspline-2)*(Nspline-2));
	for (int j = 0; j < Nspline-2; j++)
	{
		for ( int i = j; i <= std::min(Nspline-2-1,j+n); i++)  // n = #sub-diagonals
		{
			S[(i-j) + j*(Nspline-2)] = OverlapMatrixOld(i+1,j+1,n,x);
		}
	}
}

double LaplaceMatrixOld(uint32_t i, uint32_t j, uint32_t n, const vector<double> &x)
{
	double s1 = 0;
	double s2 = 0;
	double s3 = 0;
	double s4 = 0;

	if (x[i] != x[i+n])
	{
		if (x[j] != x[j+n])
		{
			s1 = OverlapMatrixOld(i,j,n-1,x)/(x[i+n] - x[i])/(x[j+n]-x[j]);
		}
		if (x[j+n+1] != x[j+1])
		{
			s2 = -OverlapMatrixOld(i,j+1,n-1,x)/(x[i+n] - x[i])/(x[j+n+1] - x[j+1]);
		}
	}
	if (x[i+n+1] != x[i+1])
	{
		if (x[j+n+1] != x[j+1])
		{
			s3 = OverlapMatrixOld(i+1,j+1,n-1,x)/(x[i+n+1] -x[i+1])/(x[j+n+1]-x[j+1]);
		}
		if (x[j+n] != x[j])
		{
			s4 = -OverlapMatrixOld(i+1,j,n-1,x)/(x[i+n+1]-x[i+1])/(x[j+n]-x[j]);
		}
	}
	return n*n*(s1 + s2 + s3 + s4);
}

void CreateLaplacianOld(const vector<double> &x, const uint32_t n, vector<double> &L)
{
	uint32_t Nspline = (uint32_t) x.size() - n - 1;

	L = vector<double>((Nspline-2)*(Nspline-2));
	for (int j = 0; j < Nspline-2; j++)
	{
		for ( int i = j; i <= std::min(Nspline-2-1,j+n); i++)  // n = #sub-diagonals
		{
			L[(i-j) + j*(Nspline-2)] = LaplaceMatrixOld(i+1,j+1,n,x);
		}
	}
}

double PotentialMatrixOld(uint32_t i, uint32_t j, uint32_t n, const vector<double> &x, double Z, int l)
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
			integral += BSplineEvaluate(i,n,&x[0],(uint32_t)x.size(),u[idx])*BSplineEvaluate(j,n,&x[0],(uint32_t)x.size(),u[idx]) * w[idx] *(l*(l+1)/u[idx]/u[idx]/2 - Z/u[idx]);
		}
	}
	return integral;
}

void CreatePotentialMatrixOld(const vector<double> &x, const uint32_t n, vector<double> &V, double Z, int l)
{
	uint32_t Nspline = (uint32_t) x.size() - n - 1;

	V = vector<double>((Nspline-2)*(Nspline-2));
	for (int j = 0; j < Nspline-2; j++)
	{
		for ( int i = j; i <= std::min(Nspline-2-1,j+n); i++)  // n = #sub-diagonals
		{
			V[(i-j) + j*(Nspline-2)] = PotentialMatrixOld(i+1,j+1,n,x,Z,l);
		}
	}
}

#endif

void SolveCatOwnerEquation(vector<double> &H, vector<double> &S, uint32_t dim, uint32_t bandsize, vector<double> &eigVal,vector<double> &eigVec)
{
	if ((H.size() != S.size()) || ((uint32_t) H.size() != dim*dim))
		printf("error in SolveCatOwnerEquation: dimensions do not match: #H = %lu, #S = %lu, dim = %du",H.size(), S.size(), dim);

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

// The Slater cell tables, wrappers around BSplineCellValues and BSplineSlaterCellTables
void PutBSplinesOnQuadGrid(const vector<double> &x, uint32_t n,uint32_t Nweights, vector<double> &BSval)
{
	BSval = vector<double>(((uint32_t) x.size() - 1)*(n+1)*Nweights);
	BSplineCellValues(&x[0], (uint32_t) x.size(), n, Nweights, &BSval[0]);
}

void CalcOffDiagSlaterInts(const vector<double> &x, const vector<double> &BSval,uint32_t Nweights, vector<double> &rk1, vector<double> &rk2, int k, uint32_t n)
{
	uint32_t ks = n + 1;
	rk1 = vector<double>(((uint32_t) x.size() - 1)*ks*ks);
	rk2 = vector<double>(((uint32_t) x.size() - 1)*ks*ks);
	BSplineSlaterCellTables(&x[0], (uint32_t) x.size(), n, Nweights, &BSval[0], k, &rk1[0], &rk2[0], NULL);
}

void CalcOffDiagSlaterInts(const vector<double> &x, vector<double> &rk1, vector<double> &rk2, int k, uint32_t n)
{
	vector<double> BSval;
	PutBSplinesOnQuadGrid(x, n, n+1+3, BSval);
	CalcOffDiagSlaterInts(x, BSval, n+1+3, rk1, rk2, k, n);
}

void CalcDiagSlaterInts(const vector<double> &x, const vector<double> &BSval,uint32_t Nweights, vector<double> &rkD, int k, uint32_t n)
{
	uint32_t ks = n + 1;
	rkD = vector<double>(((uint32_t) x.size() - 1)*ks*ks*ks*ks);
	BSplineSlaterCellTables(&x[0], (uint32_t) x.size(), n, Nweights, &BSval[0], k, NULL, NULL, &rkD[0]);
}

void CalcDiagSlaterInts(const vector<double> &x, vector<double> &rkD, int k, uint32_t n)
{
	vector<double> BSval;
	PutBSplinesOnQuadGrid(x, n, n+1+3, BSval);
	CalcDiagSlaterInts(x, BSval, n+1+3, rkD, k, n);
}

#if 0
void PutBSplinesOnQuadGridOld(const vector<double> &x, uint32_t n,uint32_t Nweights, vector<double> &BSval)
{
	uint32_t ks = n + 1;
	uint32_t nv = (uint32_t) x.size() - 1;
	uint32_t Nspline = (uint32_t) x.size() - n - 1;
//	uint32_t Nweights = ks +3;
	BSval = vector<double>(nv *ks*Nweights);
#pragma omp parallel for
	for (uint32_t iv = n; iv < nv; iv++)
	{
		vector<double> u,w;
		Quad::GaussianQuadratureWeights(Nweights,u,w,x[iv],x[iv+1]);
		uint32_t lowerIdx = iv > n ? iv-n : 0;
		for (uint32_t i = lowerIdx; i <= Min(iv,Nspline-1); i++)
		{

			for (uint32_t idx = 0; idx < (uint32_t) u.size(); idx++)
			{
				BSval[iv*ks*Nweights + (i-lowerIdx)*Nweights + idx] = BSplineEvaluate(i,n,&x[0],(uint32_t)x.size(),u[idx]);

			}

		}
	}
}

void CalcOffDiagSlaterIntsOld(const vector<double> &x, vector<double> &rk1, vector<double> &rk2, int k, uint32_t n)
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
					rk1[iv*ks*ks + (i-lowerIdx)*ks +j-lowerIdx] += std::pow(u[idx],k) * w[idx] * BSplineEvaluate(i,n,&x[0],(uint32_t)x.size(),u[idx])* BSplineEvaluate(j,n,&x[0],(uint32_t)x.size(),u[idx]);
					rk2[iv*ks*ks + (i-lowerIdx)*ks +j-lowerIdx] += std::pow(u[idx],-k-1) * w[idx] * BSplineEvaluate(i,n,&x[0],(uint32_t)x.size(),u[idx])* BSplineEvaluate(j,n,&x[0],(uint32_t)x.size(),u[idx]);
				}
				rk1[iv*ks*ks + (j-lowerIdx)*ks +i-lowerIdx] = rk1[iv*ks*ks + (i-lowerIdx)*ks +j-lowerIdx];
				rk2[iv*ks*ks + (j-lowerIdx)*ks +i-lowerIdx] = rk2[iv*ks*ks + (i-lowerIdx)*ks +j-lowerIdx];
			}
		}
	}
}

void CalcOffDiagSlaterIntsOld(const vector<double> &x, const vector<double> &BSval,uint32_t Nweights, vector<double> &rk1, vector<double> &rk2, int k, uint32_t n)
{
	uint32_t ks = n + 1;
	uint32_t nv = (uint32_t) x.size() - 1;
	uint32_t Nspline = (uint32_t) x.size() - n - 1;
	rk1 = vector<double>(nv *ks*ks);
	rk2 = vector<double>(nv *ks*ks);
//	uint32_t Nweights = ks +3;
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
					rk1[iv*ks*ks + (i-lowerIdx)*ks +j-lowerIdx] += std::pow(u[idx],k)    * w[idx] * BSval[iv*ks*Nweights + (i-lowerIdx)*Nweights + idx]* BSval[iv*ks*Nweights + (j-lowerIdx)*Nweights + idx];
					rk2[iv*ks*ks + (i-lowerIdx)*ks +j-lowerIdx] += std::pow(u[idx],-k-1) * w[idx] * BSval[iv*ks*Nweights + (i-lowerIdx)*Nweights + idx]* BSval[iv*ks*Nweights + (j-lowerIdx)*Nweights + idx];
				}
				rk1[iv*ks*ks + (j-lowerIdx)*ks +i-lowerIdx] = rk1[iv*ks*ks + (i-lowerIdx)*ks +j-lowerIdx];
				rk2[iv*ks*ks + (j-lowerIdx)*ks +i-lowerIdx] = rk2[iv*ks*ks + (i-lowerIdx)*ks +j-lowerIdx];
			}
		}
	}
}

void CalcDiagSlaterIntsOld(const vector<double> &x, vector<double> &rkD, int k, uint32_t n)
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
								I2 += w2[idx2] * std::pow(u2[idx2],k) * BSplineEvaluate(j1,n,&x[0],(uint32_t)x.size(),u2[idx2])* BSplineEvaluate(j2,n,&x[0],(uint32_t)x.size(),u2[idx2]);
							}
							Itot += I2 * w1[idx1] * std::pow(u1[idx1],-k-1) * BSplineEvaluate(i1,n,&x[0],(uint32_t)x.size(),u1[idx1])* BSplineEvaluate(i2,n,&x[0],(uint32_t)x.size(),u1[idx1]);
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

void CalcDiagSlaterIntsOld(const vector<double> &x, const vector<double> &BSval,uint32_t Nweights, vector<double> &rkD, int k, uint32_t n)
{
	uint32_t ks = n + 1;
	uint32_t nv = (uint32_t) x.size() - 1;
	uint32_t Nspline = (uint32_t) x.size() - n - 1;
	uint32_t ks2 = ks*ks;
	uint32_t ks3 = ks2*ks;
	uint32_t ks4 = ks3*ks;
//	uint32_t Nweights = ks +3;
	rkD = vector<double>(nv*ks4);
#pragma omp parallel for
	for (uint32_t iv = n; iv < nv; iv++)
	{
		vector<double> u1,w1, u2, w2;
		Quad::GaussianQuadratureWeights(Nweights,u1,w1,x[iv],x[iv+1]);
		uint32_t lowerIdx = iv > n ? iv-n : 0;
		double BSval2[Nweights*ks];
		for (uint32_t idx1 = 0; idx1 < u1.size(); idx1++)
		{
			Quad::GaussianQuadratureWeights(Nweights,u2,w2,x[iv],u1[idx1]);
			for (uint32_t i1 = lowerIdx; i1 <= Min(iv,Nspline-1); i1++)
			{
				for (uint32_t idx2 = 0; idx2 < u2.size(); idx2++)
				{
					BSval2[(i1-lowerIdx)* Nweights + idx2] = BSplineEvaluate(i1,n,&x[0],(uint32_t)x.size(),u2[idx2]);
				}
			}
			for (uint32_t i1 = lowerIdx; i1 <= Min(iv,Nspline-1); i1++)
			{
				for (uint32_t j1 = lowerIdx; j1 <= Min(iv,Nspline-1); j1++) //TODO: Check symmetry
				{
					for (uint32_t i2 = i1; i2 <= Min(iv,Nspline-1); i2++)
					{
						for (uint32_t j2 = j1; j2 <= Min(iv,Nspline-1); j2++)
						{
							double I2 = 0;
							for (uint32_t idx2 = 0; idx2 < u2.size(); idx2++)
							{
								I2 += w2[idx2] * std::pow(u2[idx2],k) * BSval2[(j1-lowerIdx)* Nweights + idx2]* BSval2[(j2-lowerIdx)* Nweights + idx2];
							}
							double Itot = I2 * w1[idx1] * std::pow(u1[idx1],-k-1) * BSval[iv*ks*Nweights + (i1-lowerIdx)*Nweights + idx1]* BSval[iv*ks*Nweights + (i2-lowerIdx)*Nweights + idx1];

							rkD[iv*ks4 +(i1-lowerIdx)*ks3 + (j1-lowerIdx)*ks2 + (i2-lowerIdx)*ks +(j2-lowerIdx)] += Itot;
							if (i1 != i2) rkD[iv*ks4 +(i2-lowerIdx)*ks3 + (j1-lowerIdx)*ks2 + (i1-lowerIdx)*ks +(j2-lowerIdx)] += Itot;
							if (j1 != j2) rkD[iv*ks4 +(i1-lowerIdx)*ks3 + (j2-lowerIdx)*ks2 + (i2-lowerIdx)*ks +(j1-lowerIdx)] += Itot;
							if (i1 != i2 && j1 != j2) rkD[iv*ks4 +(i2-lowerIdx)*ks3 + (j2-lowerIdx)*ks2 + (i1-lowerIdx)*ks +(j1-lowerIdx)] += Itot;
						}
					}
				}
			}
		}
	}
}

#endif

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

// AssembleSlaterInts2 is BSplineSlaterIntegral in the B-spline core
#if 0
double AssembleSlaterInts2Old(const vector<double> &Ci1,const vector<double> &Cj1,const vector<double> &Ci2,const vector<double> &Cj2,
						  const vector<double> &rk1, const vector<double> &rk2, const vector<double> &rkD, uint32_t n)
{
	double val = 0;
	uint32_t ks = n + 1;
	uint32_t ks2 = ks*ks;
	uint32_t ks3 = ks2*ks;
	uint32_t ks4 = ks3*ks;
	uint32_t nv = (uint32_t) rk1.size() /ks /ks;
	uint32_t Nspline = nv - n;
	for (uint32_t iv = 0; iv < nv; iv++)
	{
		uint32_t lowerIdx = iv > n ? iv-n : 0;
		uint32_t iMax = lowerIdx + ks <= Nspline ? ks : Nspline -lowerIdx;
		for (uint32_t i1 = 0; i1 < iMax; i1++)
			for (uint32_t j1 = 0; j1 < iMax; j1++)
				for (uint32_t i2 = 0; i2 < iMax; i2++)
					for (uint32_t j2 = 0; j2 < iMax; j2++)
					{
						val += (rkD[iv*ks4 +(i1)*ks3 + (j1)*ks2 + (i2)*ks +(j2)] + rkD[iv*ks4 +(j1)*ks3 + (i1)*ks2 + (j2)*ks +(i2)])
							* Ci1[i1+lowerIdx] * Ci2[i2+lowerIdx] * Cj1[j1+lowerIdx] * Cj2[j2+lowerIdx] ;
					}
	}

	for (uint32_t iv = 0; iv < nv; iv++)
	{
		uint32_t lowerIdx = iv > n ? iv-n : 0;
		uint32_t iMax = lowerIdx + ks <= Nspline ? ks : Nspline -lowerIdx;
		double I11 = 0;
		double I21 = 0;
		for (uint32_t i1 = 0; i1 < iMax; i1++)
			for (uint32_t i2 = 0; i2 < iMax; i2++)
			{
				I11 += rk1[iv*ks*ks + (i1)*ks +i2] * Ci1[i1+lowerIdx] * Ci2[i2+lowerIdx];
				I21 += rk1[iv*ks*ks + (i1)*ks +i2] * Cj1[i1+lowerIdx] * Cj2[i2+lowerIdx];
			}
		double I12 = 0;
		double I22 = 0;
		for (uint32_t jv = iv+1; jv < nv; jv++)
		{
			uint32_t lowerIdxj = jv > n ? jv-n : 0;
			uint32_t jMax = lowerIdxj + ks <= Nspline ? ks : Nspline -lowerIdxj;
			for (uint32_t j1 = 0; j1 < jMax; j1++)
				for (uint32_t j2 = 0; j2 < jMax; j2++)
				{
					I12 += rk2[jv*ks*ks + (j1)*ks +j2] * Cj1[j1+lowerIdxj] * Cj2[j2+lowerIdxj];
					I22 += rk2[jv*ks*ks + (j1)*ks +j2] * Ci1[j1+lowerIdxj] * Ci2[j2+lowerIdxj];
				}
		}
		val += I11*I12 + I21*I22;

	}
	return val;
}

#endif

double AssembleSlaterInts(const vector<double> &Ci1,const vector<double> &Cj1,const vector<double> &Ci2,const vector<double> &Cj2,
						  const vector<double> &rk1, const vector<double> &rk2, const vector<double> &rkD, uint32_t n)
{
	double val = 0;
	uint32_t ks = n +1;
	uint32_t Nc = (uint32_t) Ci1.size();
	for (uint32_t i1 = 0; i1 < Nc; i1++)
	{
		uint32_t imin = i1 > ks ? i1 -ks : 0;
		uint32_t imax = Min(i1+ks,Nc);
		for (uint32_t j1 = 0; j1 < Nc; j1++)
		{
			uint32_t jmin = j1 > ks ? j1 -ks : 0;
			uint32_t jmax = Min(j1+ks,Nc);
			for (uint32_t i2 = imin; i2 < imax; i2++)
				for (uint32_t j2 = jmin; j2 < jmax; j2++)
					val += AssembleSlaterInts(i1,j1,i2,j2,rk1,rk2,rkD,n)*Ci1[i1]*Cj1[j1]*Ci2[i2]*Cj2[j2];
		}
	}
	return val;

}

void CreateSuperMatrix(const vector<vector<double> > &rk1, const vector<vector<double> > &rk2, const vector<vector<double> > &rkD,
					  uint32_t n, uint32_t dim, uint32_t lMax, double Alpha, double Beta, vector<double> &P, vector<double> &Q)
{
	uint32_t nl = lMax+1;
	uint32_t superdim = dim*dim*nl;
	P = vector<double>(superdim*superdim);
	Q = vector<double>(superdim*superdim);
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
					for (uint32_t la = 0; la < nl; la++)
					{
						for (uint32_t lb = 0; lb < nl; lb++)
						{
							uint32_t k0 = Abs((int)la - (int)lb);
							uint32_t k1 = la + lb;
							if ((k1 +k0) % 2) k0++;
							double Fk = 0;
							for (uint32_t k = k0; k <= k1; k +=2)
							{
								double cred = ThreeJSymbol(la,0,lb,0,k,0);
								Fk += AssembleSlaterInts(nu,alpha,beta,mu,rk1[k],rk2[k],rkD[k],n)*cred*cred;
							}
							P[(nu-1)*superdim*dim*nl + (mu-1)*superdim*nl + la*superdim + (alpha-1)*dim*nl + (beta-1)*nl +lb] = -Fk * (double)(2*lb+1);
							P[(mu-1)*superdim*dim*nl + (nu-1)*superdim*nl + la*superdim + (beta-1)*dim*nl + (alpha-1)*nl +lb] = -Fk * (double)(2*lb+1);

							Q[(nu-1)*superdim*dim*nl + (mu-1)*superdim*nl + la*superdim + (alpha-1)*dim*nl + (beta-1)*nl +lb] = -Fk * Beta * (double)(2*lb+1);
							Q[(mu-1)*superdim*dim*nl + (nu-1)*superdim*nl + la*superdim + (beta-1)*dim*nl + (alpha-1)*nl +lb] = -Fk * Beta * (double)(2*lb+1);

						}
					}
				}
			}
		}
	}
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
					double Fk = AssembleSlaterInts(nu,alpha,mu,beta,rk1[0],rk2[0],rkD[0],n);
					for (uint32_t la = 0; la < nl; la++)
					{
						for (uint32_t lb = 0; lb < nl; lb++)
						{
							P[(nu-1)*superdim*dim*nl + (mu-1)*superdim*nl + la*superdim + (alpha-1)*dim*nl + (beta-1)*nl +lb] += 2*(2*lb+1)*Fk;
							if (mu != nu) P[(mu-1)*superdim*dim*nl + (nu-1)*superdim*nl + la*superdim + (alpha-1)*dim*nl + (beta-1)*nl +lb] += 2*(2*lb+1)*Fk;

							Q[(nu-1)*superdim*dim*nl + (mu-1)*superdim*nl + la*superdim + (alpha-1)*dim*nl + (beta-1)*nl +lb] += 2*(2*lb+1)*Fk * Alpha;
							if (mu != nu) Q[(mu-1)*superdim*dim*nl + (nu-1)*superdim*nl + la*superdim + (alpha-1)*dim*nl + (beta-1)*nl +lb] += 2*(2*lb+1)*Fk * Alpha;
						}
					}
				}
			}
		}
	}
}

void CreateSuperDensity(const vector<vector<double> > &eigVec, const vector<vector<uint32_t> > &Nclosed, const vector<vector<uint32_t> > &Nopen,
						uint32_t lMax, uint32_t dim, vector<double> &Dc, vector<double> &Do, double f)
{
	uint32_t nl = lMax+1;
	uint32_t superdim = dim*dim*nl;
	Dc = vector<double>(superdim);
	for (uint32_t j = 0; j < (uint32_t) Nclosed.size(); j++)
	{
		uint32_t l = Nclosed[j][0];
		uint32_t i = Nclosed[j][1];
		for (uint32_t nu = 0; nu < dim; nu++)
		{
			for (uint32_t mu = 0; mu < dim; mu++)
			{
				Dc[nu*dim*nl + mu*nl +l] += eigVec[l][i*dim + mu] * eigVec[l][i*dim + nu];
			}
		}
	}
	Do = vector<double>(superdim);
	for (uint32_t j = 0; j < (uint32_t) Nopen.size(); j++)
	{
		uint32_t l = Nopen[j][0];
		uint32_t i = Nopen[j][1];
		//printf("Density\tl = %u\t i = %u\t f = %f\n",l,i,f);
		for (uint32_t nu = 0; nu < dim; nu++)
		{
			for (uint32_t mu = 0; mu < dim; mu++)
			{
				Do[nu*dim*nl + mu*nl +l] +=  f * eigVec[l][i*dim + mu] * eigVec[l][i*dim + nu];
			}
		}
	}
}

void CreateQSuperMatrix(uint32_t la, const vector<double> &Q, const vector<double> &Do, uint32_t dim, uint32_t lMax, vector<double> &QSM)
{
	uint32_t nl = lMax + 1;
	uint32_t superdim = dim*dim*nl;
	QSM = vector<double>(dim*dim);
	for (uint32_t mu = 0; mu < dim; mu++)
	{
		for (uint32_t nu = mu; nu < dim; nu++)
		{
			for (uint32_t lb = 0; lb < nl; lb++)
			{
				for (uint32_t alpha = 0; alpha < dim; alpha++)
				{
					for (uint32_t beta = 0; beta < dim; beta++)
					{
						QSM[mu*dim + nu] += Q[nu*superdim*dim*nl + mu*superdim*nl + la*superdim + alpha*dim*nl + beta*nl +lb] * Do[alpha*dim*nl + beta*nl + lb];
					}
				}
			}
			QSM[nu*dim + mu] = QSM[mu*dim + nu];
		}
	}
}

void CreatePQSuperMatrix(uint32_t la, const vector<vector<double> > &rk1, const vector<vector<double> > &rk2, const vector<vector<double> > &rkD,
		  uint32_t n, double Alpha, double Beta, const vector<double> &Dc, const vector<double> &Do, uint32_t dim, uint32_t lMax,
		  vector<double> &PSM, vector<double> &QSM)
{
	uint32_t nl = lMax + 1;
	uint32_t superdim = dim*dim*nl;
	QSM = vector<double>(dim*dim);
	PSM = vector<double>(dim*dim);
#pragma omp parallel for
	for (uint32_t mu = 0; mu < dim; mu++)
	{
		for (uint32_t nu = mu; nu < dim; nu++)
		{
			for (uint32_t lb = 0; lb < nl; lb++)
			{
				for (uint32_t alpha = 0; alpha < dim; alpha++)
				{
					for (uint32_t beta = 0; beta < dim; beta++)
					{
						double Fk = AssembleSlaterInts(nu+1,alpha+1,mu+1,beta+1,rk1[0],rk2[0],rkD[0],n);
						uint32_t k0 = Abs((int)la - (int)lb);
						uint32_t k1 = la + lb;
						if ((k1 +k0) % 2) k0++;
						double Gk = 0;
						for (uint32_t k = k0; k <= k1; k +=2)
						{
							double cred = ThreeJSymbol(la,0,lb,0,k,0);
							Gk += AssembleSlaterInts(nu+1,alpha+1,beta+1,mu+1,rk1[k],rk2[k],rkD[k],n)*cred*cred;
						}
						QSM[mu*dim + nu] += (2*(2*lb+1)*Fk * Alpha - Gk * Beta * (double)(2*lb+1)) * Do[alpha*dim*nl + beta*nl + lb];
						PSM[mu*dim + nu] += (2*(2*lb+1)*Fk  - Gk * (double)(2*lb+1)) * (Do[alpha*dim*nl + beta*nl + lb]+Dc[alpha*dim*nl + beta*nl + lb]);
					}
				}
			}
			QSM[nu*dim + mu] = QSM[mu*dim + nu];
			PSM[nu*dim + mu] = PSM[mu*dim + nu];
		}
	}
}

void CreateRSuperMatrix(uint32_t la, const vector<double> &QSM,const vector<double> &S, const vector<double> &Dc, const vector<double> &Do, uint32_t dim, uint32_t lMax, vector<double> &RSM)
{
	uint32_t nl = lMax + 1;
	uint32_t superdim = dim*dim*nl;
	RSM = vector<double>(dim*dim);
	for (uint32_t mu = 0; mu < dim; mu++)
	{
		for (uint32_t nu = 0; nu < dim; nu++)
		{

			for (uint32_t alpha = 0; alpha < dim; alpha++)
			{
				for (uint32_t beta = 0; beta < dim; beta++)
				{
					RSM[mu*dim + nu] += S[mu*dim + alpha] * (Do[alpha*dim*nl + beta*nl + la]+Dc[alpha*dim*nl + beta*nl + la]) * QSM[beta*dim +nu]
									  + QSM[mu*dim + alpha] * (Do[alpha*dim*nl + beta*nl + la]+Dc[alpha*dim*nl + beta*nl + la]) * S[beta*dim +nu];
				}
			}
		}
	}
//	for (uint32_t mu = 0; mu < dim; mu++)
//	{
//		for (uint32_t nu = 0; nu < dim; nu++)
//		{
//			RSM[mu*dim + nu] += RSM[nu*dim + mu];
//		}
//	}
}

void CreateHartreeFockOperator(uint32_t la, const vector<double> &P,const vector<double> &Q,const vector<double> &S, const vector<double> &H0, const vector<double> &Dc, const vector<double> &Do, uint32_t dim, uint32_t lMax, vector<double> &H)
{
	uint32_t nl = lMax + 1;
	uint32_t superdim = dim*dim*nl;
	H = vector<double>(dim*dim);
	vector<double> QSM, RSM;
	CreateQSuperMatrix(la,Q,Do,dim,lMax,QSM);
	CreateRSuperMatrix(la,QSM,S,Dc,Do,dim,lMax,RSM);
	for (uint32_t mu = 0; mu < dim; mu++)
	{
		for (uint32_t nu = mu; nu < dim; nu++)
		{
			for (uint32_t lb = 0; lb < nl; lb++)
			{
				for (uint32_t alpha = 0; alpha < dim; alpha++)
				{
					for (uint32_t beta = 0; beta < dim; beta++)
					{
						H[mu*dim + nu] += (Dc[alpha*dim*nl + beta*nl + lb] + Do[alpha*dim*nl + beta*nl + lb])
										* P[nu*superdim*dim*nl + mu*superdim*nl + la*superdim + alpha*dim*nl + beta*nl +lb]
										- Do[beta*dim*nl + alpha*nl + lb] * Q[nu*superdim*dim*nl + mu*superdim*nl + la*superdim + alpha*dim*nl + beta*nl +lb];

					}
				}
			}
			H[mu*dim + nu] += RSM[mu*dim + nu] + H0[mu*dim + nu];
			H[nu*dim + mu] = H[mu*dim + nu];
		}
	}

}

void CreateHartreeFockOperator(const vector<double> &P,const vector<double> &Q,const vector<double> &S, const vector<vector<double> > &H0, const vector<double> &Dc, const vector<double> &Do, uint32_t dim, uint32_t lMax, vector<vector<double> > &F)
{
	uint32_t nl = lMax + 1;
	uint32_t superdim = dim*dim*nl;
	for (uint32_t la = 0; la < nl; la++)
	{
		F[la] = vector<double>(dim*dim);
		vector<double> QSM, RSM;
		CreateQSuperMatrix(la,Q,Do,dim,lMax,QSM);
		CreateRSuperMatrix(la,QSM,S,Dc,Do,dim,lMax,RSM);
		for (uint32_t mu = 0; mu < dim; mu++)
		{
			for (uint32_t nu = mu; nu < dim; nu++)
			{
				for (uint32_t lb = 0; lb < nl; lb++)
				{
					for (uint32_t alpha = 0; alpha < dim; alpha++)
					{
						for (uint32_t beta = 0; beta < dim; beta++)
						{
							F[la][mu*dim + nu] += (Dc[alpha*dim*nl + beta*nl + lb] + Do[alpha*dim*nl + beta*nl + lb])
											* P[nu*superdim*dim*nl + mu*superdim*nl + la*superdim + alpha*dim*nl + beta*nl +lb]
											- Do[beta*dim*nl + alpha*nl + lb] * Q[nu*superdim*dim*nl + mu*superdim*nl + la*superdim + alpha*dim*nl + beta*nl +lb];

						}
					}
				}
				F[la][mu*dim + nu] += RSM[mu*dim + nu] + H0[la][mu*dim + nu];
				F[la][nu*dim + mu] = F[la][mu*dim + nu];
			}
		}
	}

}
void CreateHartreeFockOperator(const vector<vector<double> > &H0,const vector<vector<double> > &PSM,const vector<vector<double> > &QSM,
		const vector<vector<double> > &RSM, uint32_t dim, uint32_t lMax, vector<vector<double> > &F)
{
	uint32_t nl = lMax + 1;
	for (uint32_t la = 0; la < nl; la++)
	{
		F[la] = vector<double>(dim*dim);
		for (uint32_t mu = 0; mu < dim; mu++)
		{
			for (uint32_t nu = mu; nu < dim; nu++)
			{
				F[la][mu*dim + nu] = H0[la][mu*dim + nu] + PSM[la][mu*dim + nu] - QSM[la][mu*dim + nu] + RSM[la][mu*dim + nu];
				F[la][nu*dim + mu] = F[la][mu*dim + nu];
			}
		}
	}

}


double CalcGSEnergy(const vector<double> &Dc, const vector<double> &Do,const vector<vector<double> > &H,const vector<vector<double> > &F,
		const vector<vector<double> > &QSM, uint32_t dim, uint32_t lMax, double f,const vector<vector<double> > &RSM)
{
	double e0 = 0;
	uint32_t nl = lMax + 1;
	uint32_t superdim = dim*dim*nl;
	for (uint32_t  l = 0; l < nl; l++)
	{
		for (uint32_t mu = 0; mu < dim; mu++)
		{
			for (uint32_t nu = 0; nu < dim; nu++)
			{
				e0 += (2*l+1) * ( (H[l][mu*dim + nu]+F[l][mu*dim + nu])*(Dc[nu*dim*nl + mu*nl + l] + Do[nu*dim*nl + mu*nl + l])
								- QSM[l][mu*dim + nu] * (Dc[nu*dim*nl + mu*nl + l] + f*Do[nu*dim*nl + mu*nl + l]) );
			}

		}

	}
	return e0;
}

void CreateDIISErrorMatrix(const vector<vector<double> > &F, const vector<double> Dc, const vector<double> Do,
						   const vector<double> S, uint32_t dim, uint32_t lMax, vector<vector<double> > &E)
{
	uint32_t nl = lMax + 1;
	//E = vector<vector<double> >(nl);
	for (uint32_t l = 0; l < nl; l++)
	{
		E[l] = vector<double>(dim*dim);
		for (uint32_t mu = 0; mu < dim; mu++)
		{
			for (uint32_t nu = 0; nu < dim; nu++)
			{
				for (uint32_t alpha = 0; alpha < dim; alpha++)
				{
					for (uint32_t beta = 0; beta < dim; beta++)
					{
						E[l][mu*dim + nu] += F[l][mu*dim + beta] * (Dc[beta*dim*nl + alpha*nl +l] + Do[beta*dim*nl + alpha*nl +l]) * S[alpha*dim + nu]
										   - S[mu*dim + beta] * (Dc[beta*dim*nl + alpha*nl +l] + Do[beta*dim*nl + alpha*nl +l]) * F[l][alpha*dim + nu];
					}
				}
			}
		}
	}
}

double MatrixSkalarProduct(const vector<double> &A, const vector<double> &B, uint32_t dim)
{
	double tr = 0;
	for (uint32_t i = 0; i < dim; i++)
	{
		for (uint32_t j = 0; j < dim; j++)
		{
			tr += A[i +j*dim]*B[j +i*dim];
		}
	}
	return tr;
}

void CreateDIISFockMatrix(std::list<vector<vector<double> > > &F, std::list<vector<vector<double> > > &E, uint32_t dim, uint32_t lMax, vector<vector<double> > &Fout)
{


	uint32_t nl = lMax + 1;
	uint32_t N  = (uint32_t) E.size()+1;
	vector<double> B = vector<double>(N*N);
	vector<double> c = vector<double>(N);
	c[0] = 1;

	char uplo = 'U';
	uint32_t NRHS  = 1;
	uint32_t lwork = N*N;
	int info, ipiv[N];
	double work[N*N];

	Fout = vector<vector<double> >(nl); //TODO: shouldn't be necessary

	for (uint32_t l = 0; l < nl; l++)
	{
		for (uint32_t i = 1; i < N; i++)
			c[i] = 0;
		c[0] = 1;
		B[0] = 0;
		Fout[l] = vector<double>(dim*dim);
		std::list<vector<vector<double> > >::iterator itEi = E.end();
		for (uint32_t i = 1; i < N; i++)
		{
			itEi--;
			std::list<vector<vector<double> > >::iterator itEj = E.end();
			for (uint32_t j = 1; j < N; j++)
			{
				itEj--;
				B[i + N*j] = MatrixSkalarProduct((*itEi)[l],(*itEj)[l],dim); //TODO: store that stuff
			}
			B[i]   = 1;
			B[N*i] = 1;
		}
		//PrintMartix(&B[0],N);
		info = 1;
		uint32_t dN = N+1;
		while (info && dN > 1)
		{
			dN--;
			dsysv(&uplo,&dN,&NRHS,&B[0],&N,ipiv,&c[0],&N,work,&lwork,&info);
		}
		if (info)
		{
			printf("Error in CreateDIISFockMatrix : dsysv info = %d\n",info);
//			printf("dN = %d\nl = %d\nB = \n",dN,l);
//			PrintMartix(&B[0],N);
			std::list<vector<vector<double> > >::iterator itF = F.end();
			Fout = *(--itF);
		}

		double ctot = 0;
//		for (uint32_t i = 1; i < dN; i++)
//			ctot += c[i];
//		printf("ctot = %f\n",ctot);

		for(uint32_t j = 0; j < dim*dim; j++)
		{
			std::list<vector<vector<double> > >::iterator itF = F.end();
			for (uint32_t i = 1; i < dN; i++)
			{
				itF--;
				Fout[l][j] += c[i] * (*itF)[l][j];

			}
		}


	}
}

void SCLoop(const vector<vector<uint32_t> > &Nclosed, const vector<vector<uint32_t> > &Nopen,
			  uint32_t NSCmax, double eps, bool QDIIS, uint32_t dim, const vector<vector<double> > &H0,
			  vector<double> &S_Cholesky, const vector<double> &S,const vector<vector<double> > &rk1, const vector<vector<double> > &rk2,
			  const vector<vector<double> > &rkD, uint32_t n, double Alpha, double Beta, vector<vector<double> > &eigVal,vector<vector<double> > &eigVec,
			  double f, vector<vector<double> > &F, std::list<vector<vector<double> > > &F_list, std::list<vector<vector<double> > > &E_list)
{
	uint32_t itype = 1;
	char jobz = 'V';
	char uplo = 'U';
	char trans ='N';
	char side = 'L';
	char diag = 'N';
	double one = 1;
	uint32_t lwork;
	lwork = dim*dim; //TODO: find something that is not stupid

	uint32_t lMax = (uint32_t) eigVec.size()-1;

	double gsEnergy0 = 0;
	double gsEnergy1 = 0;

	vector<double> Dc, Do, DcOld, DoOld;
	vector<vector<double> > PSM = vector<vector<double> >(lMax+1);
	vector<vector<double> > QSM = vector<vector<double> >(lMax+1);
	vector<vector<double> > RSM = vector<vector<double> >(lMax+1);
	vector<vector<double> > E   = vector<vector<double> >(lMax+1);

	for (uint32_t idxSC = 1; idxSC <= NSCmax; idxSC++)
	{
#pragma omp parallel for
		for (uint32_t  l = 0; l <= lMax; l++)
		{
			int info;
			double work[dim*dim];
			dsygst( &itype, &uplo, &dim, &(F[l][0]), &dim, &(S_Cholesky[0]), &dim, &info );
			if (info) printf("dsygst info = %d\n",info);
			dsyev( &jobz, &uplo, &dim, &(F[l][0]), &dim, &(eigVal[l][0]),work,&lwork,&info);
			if (info) printf("dsygv info = %d\n",info);
			dtrsm( &side, &uplo, &trans, &diag, &dim, &dim, &one, &(S_Cholesky[0]), &dim, &(F[l][0]), &dim ); //these three mimic dsygv
		}
#pragma omp parallel for
		for (uint32_t  l = 0; l <= lMax; l++)
		{
			for (uint32_t i = 0; i < dim*dim; i++)
			{
				eigVec[l][i] = F[l][i];
			}
		}
		CreateSuperDensity(eigVec,Nclosed,Nopen,lMax,dim,Dc,Do,f);
		for (uint32_t  l = 0; l <= lMax; l++)
		{
			CreatePQSuperMatrix(l,rk1,rk2,rkD,n,Alpha,Beta,Dc,Do,dim,lMax,PSM[l],QSM[l]);
			CreateRSuperMatrix(l,QSM[l],S,Dc,Do,dim,lMax,RSM[l]);
		}
		CreateHartreeFockOperator(H0,PSM,QSM,RSM,dim,lMax,F);
		CreateDIISErrorMatrix(F,Dc,Do,S,dim,lMax,E);

		F_list.push_back(F);
		E_list.push_back(E);
		if (QDIIS)
		{
			F_list.pop_front();
			E_list.pop_front();
			CreateDIISFockMatrix(F_list,E_list,dim,lMax,F);
		}

		gsEnergy1 = CalcGSEnergy(Dc,Do,H0,F,QSM,dim,lMax,f,RSM);
		printf("GS Energy = %f\n",gsEnergy1);
		if (Abs(1-gsEnergy0/gsEnergy1)<eps)
		{
			printf("Hartree-Fock Converged to precision of %.1E with ground-state energy = %f\n",eps,gsEnergy1);
			break;
		}
		gsEnergy0 = gsEnergy1;


	}
}

void SCScheme(const vector<double> &x, uint32_t n, int Z,int Nelec, const vector<vector<uint32_t> > &Nclosed, const vector<vector<uint32_t> > &Nopen,
			  uint32_t NSCmax,uint32_t N_DIIS, double eps, vector<vector<double> > &eigVal,vector<vector<double> > &eigVec, double a, double b, double f)
{
	int info;
	char uplo = 'U';
	uint32_t lMax = (uint32_t) eigVec.size()-1;
	vector<vector<double> > F;
	std::list<vector<vector<double> > > F_list, E_list;
	vector<vector<double> > H0 = vector<vector<double> >(lMax+1);
	vector<vector<double> > V = vector<vector<double> >(lMax+1);

	uint32_t kMax = lMax * 2;
	vector<vector<double> > rk1 = vector<vector<double> >(kMax+1);
	vector<vector<double> > rk2 = vector<vector<double> >(kMax+1);
	vector<vector<double> > rkD = vector<vector<double> >(kMax+1);
	vector<vector<double> > SlaterFk = vector<vector<double> >(kMax+1);
	vector<double> BSval, Lap, S_Cholesky, S_band, S;
	uint32_t Nweights = n +4;
	vector<double> P, Q;

	double alpha = (1-a)/(1-f);
	double beta = (1-b)/(1-f);

	double del = 1;
	uint32_t Nsplines, dim;


	PutBSplinesOnQuadGrid(x,n,Nweights,BSval);
	for (uint32_t k = 0; k <= kMax; k++)
	{
		CalcOffDiagSlaterInts(x,BSval,Nweights,rk1[k],rk2[k],k,n);
		CalcDiagSlaterInts(x,BSval,Nweights,rkD[k],k,n);
	}
	Nsplines = (uint32_t) x.size() - n-1;
	dim = Nsplines-2;

	CreateOverlapMatrix(x,n,S_band);
	S = vector<double>(dim*dim);

	for (uint32_t  l = 0; l <= lMax; l++)
	{
		H0[l] = vector<double>(dim*dim);
	}

	for (int j = 0; j < dim; j++)
	{
		for ( int i = j; i <= std::min(dim-1,j+n); i++)  // n = #sub-diagonals
		{
			S[i*dim+j] = S_band[(i-j) + j*dim];
			S[j*dim+i] = S[i*dim+j];
		}
	}
	S_Cholesky = S;
	dpotrf( &uplo, &dim, &(S_Cholesky[0]), &dim, &info );
	if (info) printf("dpotrf info = %d\n",info);
	CreateLaplacian(x,n,Lap);


	for (uint32_t  l = 0; l <= lMax; l++)
		CreatePotentialMatrix(x,n,V[l],Z,l);

	for (uint32_t  l = 0; l <= lMax; l++)
	{
		for (int j = 0; j < dim; j++)
		{
			for ( int i = j; i <= std::min(dim-1,j+n); i++)  // n = #sub-diagonals
			{
				H0[l][i*dim+j] = Lap[(i-j) + j*dim]/2. + V[l][(i-j) + j*dim];
				H0[l][j*dim+i] = H0[l][i*dim+j];

			}
		}
	}
	F = H0;
	SCLoop(Nclosed,Nopen,N_DIIS,0,false,dim,H0,S_Cholesky,S,rk1,rk2,rkD,n,alpha,beta,eigVal,eigVec,f,F,F_list,E_list);
	SCLoop(Nclosed,Nopen,NSCmax,eps,true,dim,H0,S_Cholesky,S,rk1,rk2,rkD,n,alpha,beta,eigVal,eigVec,f,F,F_list,E_list);

//	double h00 = 0;
//	double h01 = 0;
//	double hpp = 0;
//	for (uint32_t mu = 0; mu < dim; mu++)
//		for (uint32_t nu = 0; nu < dim; nu++)
//		{
//			h00 += eigVec[0][0*dim + mu] * H0[0][mu*dim+nu] * eigVec[0][0*dim + nu];
//			h01 += eigVec[0][0*dim + mu] * H0[0][mu*dim+nu] * eigVec[0][1*dim + nu];
//			hpp += eigVec[1][0*dim + mu] * H0[1][mu*dim+nu] * eigVec[1][0*dim + nu];
//		}
//	printf("h00 = %E\th01 = %E\thpp = %E\n",h00,h01,hpp);

}

void InitializeAtom(const string &atom, int  &Z, int &Nelec, vector<vector<uint32_t> > &Nclosed, vector<vector<uint32_t> > &Nopen,double &a, double &b, double &f)
{
	if (atom == "H")
	{
		Z = 1;
		Nelec = 1;
		Nclosed = vector<vector<uint32_t> >(1);
		Nclosed[0].push_back(0);
		Nclosed[0].push_back(0);
	}
	else if (atom == "He")
	{
		Z = 2;
		Nelec = 2;
		Nclosed = vector<vector<uint32_t> >(1);
		Nclosed[0].push_back(0);
		Nclosed[0].push_back(0);
	}
	else if (atom == "Li")
	{
		Z = 3;
		Nelec = 3;
		a = 1.;
		b = 2.;
		f = 1./2.;
		Nclosed = vector<vector<uint32_t> >(1);
		Nclosed[0].push_back(0);
		Nclosed[0].push_back(0);
		Nopen = vector<vector<uint32_t> >(1);
		Nopen[0].push_back(0);
		Nopen[0].push_back(1);

	}
	else if (atom == "Be")
	{
		Z = 4;
		Nelec = 4;
		Nclosed = vector<vector<uint32_t> >(2);
		Nclosed[0].push_back(0);
		Nclosed[0].push_back(0);
		Nclosed[1].push_back(0);
		Nclosed[1].push_back(1);
	}
	else if (atom == "B")
	{
		Z = 5;
		Nelec = 5;
		a = 0;
		b = 0;
		f = 1/6.;
		Nclosed = vector<vector<uint32_t> >(2);
		Nclosed[0].push_back(0);
		Nclosed[0].push_back(0);
		Nclosed[1].push_back(0);
		Nclosed[1].push_back(1);
		Nopen = vector<vector<uint32_t> >(1);
		Nopen[0].push_back(1);
		Nopen[0].push_back(0);
	}
	else if (atom == "C")
	{
		Z = 6;
		Nelec = 6;
		a = 3/4.;
		b = 3/2.;
		f = 1/3.;
		Nclosed = vector<vector<uint32_t> >(2);
		Nclosed[0].push_back(0);
		Nclosed[0].push_back(0);
		Nclosed[1].push_back(0);
		Nclosed[1].push_back(1);
		Nopen = vector<vector<uint32_t> >(1);
		Nopen[0].push_back(1);
		Nopen[0].push_back(0);

	}
	else if (atom == "N")
	{
		Z = 7;
		Nelec = 7;
		a = 1;
		b = 2;
		f = 1/2.;
		Nclosed = vector<vector<uint32_t> >(2);
		Nclosed[0].push_back(0);
		Nclosed[0].push_back(0);
		Nclosed[1].push_back(0);
		Nclosed[1].push_back(1);
		Nopen = vector<vector<uint32_t> >(1);
		Nopen[0].push_back(1);
		Nopen[0].push_back(0);

	}
	else if (atom == "O")
	{
		Z = 8;
		Nelec = 8;
		a = 15/16.;
		b = 9/8.;
		f = 2/3.;
		Nclosed = vector<vector<uint32_t> >(2);
		Nclosed[0].push_back(0);
		Nclosed[0].push_back(0);
		Nclosed[1].push_back(0);
		Nclosed[1].push_back(1);
		Nopen = vector<vector<uint32_t> >(1);
		Nopen[0].push_back(1);
		Nopen[0].push_back(0);

	}
	else if (atom == "F")
	{
		Z = 9;
		Nelec = 9;
		a = 24/25.;
		b = 24/25.;
		f = 5/6.;
		Nclosed = vector<vector<uint32_t> >(2);
		Nclosed[0].push_back(0);
		Nclosed[0].push_back(0);
		Nclosed[1].push_back(0);
		Nclosed[1].push_back(1);
		Nopen = vector<vector<uint32_t> >(1);
		Nopen[0].push_back(1);
		Nopen[0].push_back(0);

	}
	else if (atom == "Ne")
	{
		Z = 10;
		Nelec = 10;
		Nclosed = vector<vector<uint32_t> >(3);
		Nclosed[0].push_back(0);
		Nclosed[0].push_back(0);
		Nclosed[1].push_back(0);
		Nclosed[1].push_back(1);
		Nclosed[2].push_back(1);
		Nclosed[2].push_back(0);
	}
	else if (atom == "Na")
	{
		Z = 11;
		Nelec = 11;
		a = 1;
		b = 2;
		f = 1/2.;
		Nclosed = vector<vector<uint32_t> >(3);
		Nclosed[0].push_back(0);
		Nclosed[0].push_back(0);
		Nclosed[1].push_back(0);
		Nclosed[1].push_back(1);
		Nclosed[2].push_back(1);
		Nclosed[2].push_back(0);

		Nopen = vector<vector<uint32_t> >(1);
		Nopen[0].push_back(0);
		Nopen[0].push_back(2);


	}
	else if (atom == "Mg")
	{
		Z = 12;
		Nelec = 12;
		Nclosed = vector<vector<uint32_t> >(4);
		Nclosed[0].push_back(0);
		Nclosed[0].push_back(0);
		Nclosed[1].push_back(0);
		Nclosed[1].push_back(1);
		Nclosed[2].push_back(0);
		Nclosed[2].push_back(2);
		Nclosed[3].push_back(1);
		Nclosed[3].push_back(0);
	}
	else if (atom == "Al")
	{
		Z = 13;
		Nelec = 13;
		a = 0;
		b = 0;
		f = 1/6.;
		Nclosed = vector<vector<uint32_t> >(4);
		Nclosed[0].push_back(0);
		Nclosed[0].push_back(0);
		Nclosed[1].push_back(0);
		Nclosed[1].push_back(1);
		Nclosed[2].push_back(0);
		Nclosed[2].push_back(2);
		Nclosed[3].push_back(1);
		Nclosed[3].push_back(0);

		Nopen = vector<vector<uint32_t> >(1);
		Nopen[0].push_back(1);
		Nopen[0].push_back(1);

	}
	else if (atom == "Si")
	{
		Z = 14;
		Nelec = 14;
		a = 3/4.;
		b = 3/2.;
		f = 1/3.;
		Nclosed = vector<vector<uint32_t> >(4);
		Nclosed[0].push_back(0);
		Nclosed[0].push_back(0);
		Nclosed[1].push_back(0);
		Nclosed[1].push_back(1);
		Nclosed[2].push_back(0);
		Nclosed[2].push_back(2);
		Nclosed[3].push_back(1);
		Nclosed[3].push_back(0);

		Nopen = vector<vector<uint32_t> >(1);
		Nopen[0].push_back(1);
		Nopen[0].push_back(1);

	}
	else if (atom == "P")
	{
		Z = 15;
		Nelec = 15;
		a = 1;
		b = 2;
		f = 1/2.;
		Nclosed = vector<vector<uint32_t> >(4);
		Nclosed[0].push_back(0);
		Nclosed[0].push_back(0);
		Nclosed[1].push_back(0);
		Nclosed[1].push_back(1);
		Nclosed[2].push_back(0);
		Nclosed[2].push_back(2);
		Nclosed[3].push_back(1);
		Nclosed[3].push_back(0);

		Nopen = vector<vector<uint32_t> >(1);
		Nopen[0].push_back(1);
		Nopen[0].push_back(1);
	}
	else if (atom == "S")
	{
		Z = 16;
		Nelec = 16;
		a = 15/16.;
		b = 9/8.;
		f = 2/3.;
		Nclosed = vector<vector<uint32_t> >(4);
		Nclosed[0].push_back(0);
		Nclosed[0].push_back(0);
		Nclosed[1].push_back(0);
		Nclosed[1].push_back(1);
		Nclosed[2].push_back(0);
		Nclosed[2].push_back(2);
		Nclosed[3].push_back(1);
		Nclosed[3].push_back(0);

		Nopen = vector<vector<uint32_t> >(1);
		Nopen[0].push_back(1);
		Nopen[0].push_back(1);

	}
	else if (atom == "Cl")
	{
		Z = 17;
		Nelec = 17;
		a = 24/25.;
		b = 24/25.;
		f = 5/6.;
		Nclosed = vector<vector<uint32_t> >(4);
		Nclosed[0].push_back(0);
		Nclosed[0].push_back(0);
		Nclosed[1].push_back(0);
		Nclosed[1].push_back(1);
		Nclosed[2].push_back(0);
		Nclosed[2].push_back(2);
		Nclosed[3].push_back(1);
		Nclosed[3].push_back(0);

		Nopen = vector<vector<uint32_t> >(1);
		Nopen[0].push_back(1);
		Nopen[0].push_back(1);

	}
	else if (atom == "Ar")
	{
		Z = 18;
		Nelec = 18;
		Nclosed = vector<vector<uint32_t> >(5);
		Nclosed[0].push_back(0);
		Nclosed[0].push_back(0);
		Nclosed[1].push_back(0);
		Nclosed[1].push_back(1);
		Nclosed[2].push_back(0);
		Nclosed[2].push_back(2);
		Nclosed[3].push_back(1);
		Nclosed[3].push_back(0);
		Nclosed[4].push_back(1);
		Nclosed[4].push_back(1);

	}
	else if (atom == "K")
	{
		Z = 19;
		Nelec = 19;
		a = 1;
		b = 2;
		f = 1/2.;
		Nclosed = vector<vector<uint32_t> >(5);
		Nclosed[0].push_back(0);
		Nclosed[0].push_back(0);
		Nclosed[1].push_back(0);
		Nclosed[1].push_back(1);
		Nclosed[2].push_back(0);
		Nclosed[2].push_back(2);
		Nclosed[3].push_back(1);
		Nclosed[3].push_back(0);
		Nclosed[4].push_back(1);
		Nclosed[4].push_back(1);

		Nopen = vector<vector<uint32_t> >(1);
		Nopen[0].push_back(0);
		Nopen[0].push_back(3);

	}
	else if (atom == "Ca")
	{
		Z = 20;
		Nelec = 20;
		Nclosed = vector<vector<uint32_t> >(6);
		Nclosed[0].push_back(0);
		Nclosed[0].push_back(0);
		Nclosed[1].push_back(0);
		Nclosed[1].push_back(1);
		Nclosed[2].push_back(0);
		Nclosed[2].push_back(2);
		Nclosed[3].push_back(1);
		Nclosed[3].push_back(0);
		Nclosed[4].push_back(1);
		Nclosed[4].push_back(1);
		Nclosed[5].push_back(0);
		Nclosed[5].push_back(3);

	}
	else
	{
		printf("Error in HartreeFock::InitializeAtom : Atom not found\n");
	}
}

#include "Outdated_Functions.cpp"

}

static int LuaAtomicHartreeFock(lua_State *L)
{
    struct{
    	uint32_t n = 7;
    	uint32_t NSCmax = 100;
    	uint32_t N_DIIS = 5;
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
		{"N_DIIS","uint32_t",&(Options.N_DIIS)},
		{"lMax","uint32_t",&(Options.lMax)},
		{"Grid_R","double",&(Options.R)},
		{"Grid_m","int",&(Options.m)},
		{"Grid_scale","double",&(Options.scale)},
		{"Grid_hMax","double",&(Options.hMax)},
		{"Precision","double",&(Options.eps)},
  	    {NULL, NULL, NULL}};
  	int Z ,Nelec;
  	vector<double> x, rho;
  	uint32_t len = (uint32_t)lua_gettop(L);
  	string atom;
  	lua_check(L,1,atom);
  	vector<vector<uint32_t> > Nclosed, Nopen;
  	double a,b,f;
  	HartreeFock::InitializeAtom(atom,Z,Nelec,Nclosed,Nopen,a,b,f);
//	lua_check(L,1,Z);
//	lua_check(L,2,Nelec);
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
	if (len == 2) LuaCheckOptions(L, -1, LuaCreateResonantSpectraOptions);
	{
		double *knot;
		uint32_t nKnot;
		BSplineKnotsFischer(Options.m,Options.hMax,Options.R,Options.scale,Options.n,&knot,&nKnot);
		x.assign(knot, knot + nKnot);
		free(knot);
	}
	uint32_t Nsplines = (uint32_t) x.size() - Options.n-1;
	uint32_t dim = Nsplines-2;
	vector<vector<double> > eigVal = vector<vector<double> >(Options.lMax+1);
	vector<vector<double> > eigVec = vector<vector<double> >(Options.lMax+1);
	for (uint32_t l = 0; l <= Options.lMax; l++)
	{
		eigVal[l] = vector<double>(dim);
		eigVec[l] = vector<double>(dim*dim);
	}
	HartreeFock::SCScheme(x,Options.n,Z,Nelec,Nclosed,Nopen,Options.NSCmax,Options.N_DIIS,Options.eps,eigVal,eigVec,a,b,f);
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
			funs[l][i]->type = "BSpline";
			funs[l][i]->knots = x;
			funs[l][i]->x = x;
			for (int foo = 0; foo < Options.n; foo++)
			{
				funs[l][i]->x.pop_back();
				funs[l][i]->x.erase(funs[l][i]->x.begin());
			}
			funs[l][i]->degree = Options.n;
			funs[l][i]->y.push_back(0);
			for (uint32_t mu = 1; mu <= dim; mu++)
				funs[l][i]->y.push_back(eigVec[l][i*dim + (mu-1)]);
			funs[l][i]->y.push_back(0);

		}
		lua_rawseti(L,-2,l+1);
	}
	lua_push(L,eigVal);

	return 2;
}


