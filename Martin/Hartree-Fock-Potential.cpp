/*
 * Hartree-Fock-Potential.cpp
 *
 *  Created on: Apr 9, 2021
 *      Author: martin
 */

#include "Hartree-Fock-Potential.h"
#include <math.h>
#include <cmath>
#include <algorithm>

// TODO case r0 = 0
double Coulomb_integral(int k, double r0, double r1, double rho0, double rho1)
{
	double m = (rho1-rho0)/(r1-r0);
	double b = rho1 - r1 * m;
	switch (k) {
		case -1:
			return m * (r1-r0) + b * log(r1/r0);
			break;
		case -2:
			return m * log(r1/r0) - b * (1/r1-1/r0);
		default:
			double pr0 = std::pow(r0,k+1);
			double pr1 = std::pow(r1,k+1);
			return m/(k+2.) * (r1*pr1 - r0*pr0) + b/(k+1.) * (pr1 - pr0);
			break;
	}
}

template<typename T>
void InitArray(T A[], int n)
{
	for (int i = 0; i < n; i++)
	{
		A[i] = 0;
	}
}
#ifndef AT_QUANTY
double Factorial(int a)
{
  int    i;
  double result;
  result = 1.0;
  if (a < 0)
  {
    return NAN;
  }
  for (i = a; i > 0; i--)
  {
    result *= i;
  }
  return result;
}

double ThreeJSymbolHalf(int j1, int m1, int j2, int m2, int j3, int m3)
{
  double ThreeJ;
  int    k;
  ThreeJ = 0;
  if (m1 + m2 + m3 != 0)
    return 0;
  if (j1 + j2 < j3)
    return 0;
  if (j2 + j3 < j1)
    return 0;
  if (j3 + j1 < j2)
    return 0;
  if (abs(m1) > j1)
    return 0;
  if (abs(m2) > j2)
    return 0;
  if (abs(m3) > j3)
    return 0;
  for (k = std::max(j1 - j3 + m2, std::max(j2 - j3 - m1, 0)) / 2; k <= std::min(j1 - m1, std::min(j2 + m2, j1 + j2 - j3)) / 2; k++)
  {
    ThreeJ +=
      (double)(1 -
               2 * (k & 1)) / (Factorial(k) * Factorial((j1 + j2 - j3) / 2 - k) * Factorial((j1 - m1) / 2 - k) * Factorial((j2 + m2) / 2 - k) *
                               Factorial((j3 - j2 + m1) / 2 + k) * Factorial((j3 - j1 - m2) / 2 + k));
  }
  ThreeJ *=
    (double)(1 - 2 * (((j2 - j1 + m3) / 2) & 1)) * sqrt(Factorial((j1 + j2 - j3) / 2) * Factorial((j1 - j2 + j3) / 2) *
                                                        Factorial((-j1 + j2 + j3) / 2) * Factorial((j1 + m1) / 2) *
                                                        Factorial((j1 - m1) / 2) * Factorial((j2 + m2) / 2) * Factorial((j2 - m2) / 2) *
                                                        Factorial((j3 + m3) / 2) * Factorial((j3 - m3) / 2) / Factorial((j1 + j2 + j3) / 2 + 1));
  return ThreeJ;
}
#endif
void Calc_total_Density(double G[], double F[], double occ[], int N, double rho_tot[], int n, int nMax[],double eps, int &idx_last_orbital)
{
	InitArray(rho_tot, n);
	rho_tot[0] = 0.;
	idx_last_orbital = 0;
	int curr_Max_n = 0;
	for (int a = 0; a < N; a++)
	{
		int i = 1;
		double psi_sqr = (G[a*n+i]*G[a*n+i] + F[a*n+i]*F[a*n+i]);
		while (i < n && psi_sqr < 10*eps)
		{
			rho_tot[i] += occ[a] * psi_sqr;
			i++;
			psi_sqr = (G[a*n+i]*G[a*n+i] + F[a*n+i]*F[a*n+i]);
		}
		while (i < n && psi_sqr > eps)
		{
			rho_tot[i] += occ[a] * psi_sqr;
			i++;
			psi_sqr = (G[a*n+i]*G[a*n+i] + F[a*n+i]*F[a*n+i]);
		}
		nMax[a] = i;
		if (i > curr_Max_n)
		{
			idx_last_orbital = a;
			curr_Max_n = i;
		}
	}
}

/*	r : radial grid of size n
 * 	G,F : radial wavefunctions of N orbitals on grid of size n => total size = n*N
 * 	occ : occupation numbers of the N orbitals
 * 	twoj: angular momentum of the N orbitals
 * 	V   : output of potential on grid => size = n
 */
void Calc_HF_Potential(double r[], int n, double G[], double F[], double occ[], int kappa[], int N, double V[],double eps)
{
	//TODO: debug 2nd Term, implement 3rd Term

	InitArray(V, n);

	double rho_tot[n];
	int nMax[N];
	int idx_last_orbital;
	Calc_total_Density(G, F, occ, N, rho_tot, n,nMax,eps,idx_last_orbital);

	for (int a = 0; a < N; a++)
	{
		double rho0 = G[a*n]*G[a*n] + F[a*n]*F[a*n];
		double rho1 = G[a*n+1]*G[a*n+1] + F[a*n+1]*F[a*n+1];
		double Y0_1 = Coulomb_integral(0, r[0], r[1], rho0, rho1)/r[1];
		double Y0_2 = 0;
		double fk;

		int la = kappa[a] > 0 ? kappa[a] : - kappa[a] - 1;
		int twoja = 2 * abs(kappa[a]) -1;

		double Yk_1[twoja]; //TODO optimize memory
		double Yk_2[twoja];

		for (int k = 2; k < twoja; k+=2)
		{
			Yk_1[k] = Coulomb_integral(k, r[0], r[1], rho0, rho1)/std::pow(r[1],k+1);
			Yk_2[k] = 0;

		}

		for (int i = 2; i < nMax[a]; i++)
		{
			rho0 = rho1;
			rho1 = G[a*n+i]*G[a*n+i] + F[a*n+i]*F[a*n+i];
			Y0_2 += Coulomb_integral(-1, r[i-1], r[i], rho0, rho1);

			for (int k = 2; k < twoja; k+=2)
			{
				Yk_2[k] += Coulomb_integral(-k-1, r[i-1], r[i], rho0, rho1) * std::pow(r[1],k);
			}
		}

		double Y0 = Y0_1 + Y0_2;
		rho1 = G[a*n+1]*G[a*n+1] + F[a*n+1]*F[a*n+1];
		V[1] += occ[a] * Y0 * (1 - rho1/rho_tot[1]);

		Y0 = 0;
		for (int k = 2; k < twoja; k+=2)
		{
			fk = ThreeJSymbolHalf(twoja, -1, k*2, 0, twoja, 1); //TODO: optimize
			Y0 -= (1.+1./twoja)*fk*fk*(Yk_1[k] + Yk_2[k]);
		}
		V[1] += occ[a] * (occ[a] - 1) * Y0 * rho1/rho_tot[1];

		for (int i = 2; i < nMax[a]; i++)
		{
			rho0 = rho1;
			rho1 = G[a*n+i]*G[a*n+i] + F[a*n+i]*F[a*n+i];
			Y0_1 = r[i-1] / r[i] * Y0_1 + Coulomb_integral(0, r[i-1], r[i], rho0, rho1)/r[i];
			Y0_2 = Y0_2 - Coulomb_integral(-1, r[i-1], r[i], rho0, rho1);

			Y0 = Y0_1 + Y0_2;
			V[i] += occ[a] * Y0 * (1 - rho1/rho_tot[i]);

			Y0 = 0;
			for (int k = 2; k < twoja; k+=2)
			{
				Yk_1[k] = (std::pow(r[i-1],k+1) * Yk_1[k] + Coulomb_integral(k, r[i-1], r[i], rho0, rho1)) / std::pow(r[i],k+1);
				Yk_2[k] = (Yk_2[k]/std::pow(r[i-1],k) - Coulomb_integral(-k-1, r[i-1], r[i], rho0, rho1)) * std::pow(r[i],k);
				fk = ThreeJSymbolHalf(twoja, -1, k*2, 0, twoja, 1); //TODO: optimize
				Y0 -= (1.+1./twoja)*fk*fk*(Yk_1[k] + Yk_2[k]);
			}
			V[i] += occ[a] * (occ[a] - 1) * Y0 * rho1/rho_tot[i];


		}

		if (a != idx_last_orbital)
		{
			for (int i = nMax[a]; i < n; i++)
			{
				Y0_1 = r[i-1] / r[i] * Y0_1 ;
				V[i] += occ[a] * Y0_1 ;
			}
		}
		else
		{
			for (int i = nMax[a]; i < n; i++)
			{
				Y0_1 = r[i-1] / r[i] * Y0_1 ;
				V[i] += Y0_1*(occ[a]-1) ;
				Y0 = 0;
				for (int k = 2; k < twoja; k+=2)
				{
					Yk_1[k] = std::pow(r[i-1],k+1) * Yk_1[k] / std::pow(r[i],k+1);
					fk = ThreeJSymbolHalf(twoja, -1, k*2, 0, twoja, 1); //TODO: optimize
					Y0 -= (1.+1./twoja)*fk*fk*Yk_1[k];
				}
				V[i] += (occ[a] - 1) * Y0;
			}
		}
		for (int b = a+1; b < N; b++)
		{
			int lb = kappa[b] > 0 ? kappa[b] : - kappa[b] - 1;
			int twojb = 2 * abs(kappa[b]) -1;
			int kmin = abs(twojb-twoja)/2;
			if ((la+lb+kmin)%2)
				kmin =  kmin + 1;
			int kmax = (twojb+twoja)/2;
			double Ykab_1[kmax]; //TODO: optimize memory
			double Ykab_2[kmax];
			rho0 = G[a*n]*G[b*n] + F[a*n]*F[b*n];
			rho1 = G[a*n+1]*G[b*n+1] + F[a*n+1]*F[b*n+1];
			for (int k = kmin; k < kmax; k+=2)
			{
				Ykab_1[k] = Coulomb_integral(k, r[0], r[1], rho0, rho1)/std::pow(r[1],k+1);
				Ykab_2[k] = 0;

			}
			int min_nMax = std::min(nMax[a],nMax[b]);
			for (int i = 2; i < min_nMax; i++)
			{
				rho0 = rho1;
				rho1 = G[a*n+i]*G[b*n+i] + F[a*n+i]*F[b*n+i];
				for (int k = kmin; k < kmax; k+=2)
				{
					Ykab_2[k] += Coulomb_integral(-k-1, r[i-1], r[i], rho0, rho1) * std::pow(r[1],k);
				}
			}

			rho1 = G[a*n+1]*G[b*n+1] + F[a*n+1]*F[b*n+1];
			Y0 = 0;
			for (int k = kmin; k < kmax; k+=2)
			{
				fk = ThreeJSymbolHalf(twoja, -1, k*2, 0, twojb, 1); //TODO: optimize
				Y0 -= fk*fk*(Ykab_1[k] + Ykab_2[k]);
			}
			V[1] += 2*occ[a] * occ[b] * Y0 * rho1/rho_tot[1];

			for (int i = 2; i < min_nMax; i++)
			{
				rho0 = rho1;
				rho1 = G[a*n+i]*G[b*n+i] + F[a*n+i]*F[b*n+i];
				Y0 = 0;
				for (int k = kmin; k < kmax; k+=2)
				{
					Ykab_1[k] = (std::pow(r[i-1],k+1) * Ykab_1[k] + Coulomb_integral(k, r[i-1], r[i], rho0, rho1)) / std::pow(r[i],k+1);
					Ykab_2[k] = (Ykab_2[k]/std::pow(r[i-1],k) - Coulomb_integral(-k-1, r[i-1], r[i], rho0, rho1)) * std::pow(r[i],k);
					fk = ThreeJSymbolHalf(twoja, -1, k*2, 0, twojb, 1); //TODO: optimize
					Y0 -= fk*fk*(Ykab_1[k] + Ykab_2[k]);
				}
				V[i] += 2*occ[a] * occ[b] * Y0 * rho1/rho_tot[i];


			}


		}

	}
}


