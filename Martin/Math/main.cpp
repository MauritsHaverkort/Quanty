/*
 * main.cpp
 *
 *  Created on: Mar 25, 2021
 *      Author: martin
 */

#include <algorithm>
#include <math.h>
#include <chrono>
#include <iostream>
#include <fstream>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include "rk4.h"
#include "Quad.h"
#include "Hartree-Fock-Potential.h"
#include <bitset>
#include <bit>

using namespace std;
using namespace std::chrono;

int SIR_GSL(double t, const double y[], double dydt[], void* params)
{

	double *p = (double*) params;
	dydt[0] = -p[0] * y[0] * y[1];             //S' = -cSI
	dydt[1] = p[0] * y[0] * y[1] -p[1]*y[1];   //I' = cIS - wI
	dydt[2] = p[1]*y[1];                       //R' = wI
	return GSL_SUCCESS;
}

int jac (double t, const double y[], double* dfdy, double dfdt[], void* params)
{

	double *p = (double*) params;

	dfdt[0] = 0.;
	dfdt[1] = 0.;
	dfdt[2] = 0.;

	dfdy[0] = -p[0] * y[1];
	dfdy[1] = -p[0] * y[0];
	dfdy[2] = 0.;

	dfdy[3] = p[0] * y[1];
	dfdy[4] = p[0] * y[0] -p[1];
	dfdy[5] = 0.;

	dfdy[2*3 + 0] = 0.;
	dfdy[2*3 + 1] = p[1];
	dfdy[2*3 + 2] = 0.;

	return GSL_SUCCESS;
}

void TestRK4()
{
	int dim = 3;
	double y0[3] = {0.8,0.2,0.};
	double x0 = 0.;
	double h = 0.7/2;
	int N_steps = 2*100;
	double p[2] = {0.38,0.06};
	double y[dim*(N_steps+1)];
	auto start = std::chrono::high_resolution_clock::now();
	rk4_loop(&SIR,y0,dim,x0,h,N_steps,p,y);
	auto end = std::chrono::high_resolution_clock::now();
	auto delta = duration_cast<microseconds>(end - start);
	cout << "Time used : " << delta.count() << endl;
	ofstream file;
	file.open("SIR.dat");
	for (int n = 0; n <= N_steps; n++){
		file << x0 + n*h << "\t" << y[n*dim] << "\t" << y[n*dim+1] << "\t" << y[n*dim+2] << "\n";
	}
	file.close();

	// Start GSL
	double ts[100];
	double ys[3*100];
	gsl_odeiv2_system sys = {SIR_GSL, jac, 3, p};
	gsl_odeiv2_driver*d =gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk4,1e-6, 1e-6, 0.0);
	int i;
	double t = 0.0, t1 = 70.0;

	start = std::chrono::high_resolution_clock::now();
	for(i = 1; i <= 100; i++){
		double ti = i*t1 / 100.0;
		int status = gsl_odeiv2_driver_apply (d, &t, ti, y0);
		if(status != GSL_SUCCESS)
		{
			printf ("error, return value=%d\n", status);
			break;
		}

		ts[i-1] = t;
		ys[3*(i-1)] = y0[0];
		ys[3*(i-1)+1] = y0[1];
		ys[3*(i-1)+2] = y0[2];
	}
	end = std::chrono::high_resolution_clock::now();
	delta = duration_cast<microseconds>(end - start);
	cout << "Time used : " << delta.count() << endl;
	gsl_odeiv2_driver_free (d);

	file.open("SIR_GSL.dat");
	for (int n = 0; n <100; n++){
		file << ts[n] << "\t" << ys[n*dim] << "\t" << ys[n*dim+1] << "\t" << ys[n*dim+2] << "\n";
	}
	file.close();
}

int GM_Product(int A, int B)
{
	if ((A ^ B) == ( A | B))
	{
		int temp = A << 1;
		int pc = __builtin_popcount(temp & B);
		while (temp != 0)
		{
			temp <<= 1;
			pc += __builtin_popcount(temp & B);
		}
		if (pc % 2 ==0)
			return A ^ B;
		else
			return - (A ^ B);
	}
	else
		return 0;
}
int Fac(int k)
{
	int res = 1;
	for (int i = 1; i <= k; i++)
		res *= i;
	return res;
}

void NormalOrderedProduct(int A, int B, int k,double coeff, double rho[], double AB_ord[])
{

	int AB = 0;//GM_Product(A, B);
	//AB_ord[abs(AB)] = AB/abs(AB) * coeff;
	int k_fac = Fac(k);
	int sign;
	if (k%2 == 0)
		sign = ((k/2) % 2 == 0) ? 1 : -1;
	else
		sign = (((k+1)/2) % 2 == 0) ? 1 : -1;
/*	if ((A == 0) || (B == 0))
	{
		AB_ord[A+B] += sign * coeff / (double) k_fac;
	}
	else*/
	{

		int pos_cr = 1 << 7;
		int pos_an = 1 << 3;
		int pc_A1 = 0;
		int pc_B1 = __builtin_popcount(B>>4);
		int pc_A2 = __builtin_popcount(A>>4);
		int pc_B2 = 0;
		int i = 0;
		while (pos_an != 0)
		{

			if ((A & pos_cr) && (B & pos_an) )
			{
				int dA = A & (~pos_cr);
				int dB = B & (~pos_an);
				AB = GM_Product(dA, dB);
				int sign_d = (pc_A1 + pc_B1) % 2 == 0 ? 1 : -1;
				if (AB != 0 ) AB_ord[abs(AB)] += AB/abs(AB) * coeff * rho [i] * sign_d * sign / (double) k_fac;
				NormalOrderedProduct(dA, dB, k+1, coeff * rho [i] * sign_d,rho, AB_ord);
			}
			if ((A & pos_an) && (B & pos_cr) )
			{
				int dA = A & (~pos_an);
				int dB = B & (~pos_cr);
				AB = GM_Product(dA, dB);
				int sign_d = (pc_A2 + pc_B2) % 2 == 0 ? 1 : -1;
				if (AB != 0 ) AB_ord[abs(AB)] += AB/abs(AB) * coeff * (1 - rho [i]) * sign_d * sign / (double) k_fac;
				NormalOrderedProduct(dA, dB, k+1, coeff * (1 - rho [i]) * sign_d,rho, AB_ord);
			}
			if (A & pos_cr) pc_A1++;
			if (A & pos_an) pc_A2++;
			if (B & pos_cr) pc_B2++;
			if (B & pos_an) pc_B1++;
			pos_an >>= 1;
			pos_cr >>= 1;
			i++;
		}
	}
}

void Commutator(double A[], double B[], double rho[], double out[])
{
	double n_zero = 0.0001;
	for (int i = 0; i < 256; i++)
		out[i] = 0.;
	for (int i = 0; i < 256; i++)
	{
		if(fabs(A[i]) > n_zero) // TODO: can I make this exactly zero wlog?
		{
			for (int j = 0; j < 256; j++)
			{
				if(fabs(B[j]) > n_zero)
				{
					double AB_ord[256];
					double BA_ord[256];
					for(int k = 0; k < 256; k++)
					{
						AB_ord[k] = 0.;
						BA_ord[k] = 0.;
					}
					NormalOrderedProduct(i, j, 1, A[i]*B[j], rho, AB_ord);
					NormalOrderedProduct(j, i, 1, A[i]*B[j], rho, BA_ord);
					for(int k = 0; k < 256; k++)
					{
						out[k] += AB_ord[k] - BA_ord[k];
					}
				}
			}
		}
	}
}

void Wegner_Flow(double dH[], double H[], double B, double rho[])
{
	double eta[256];
	double H0[256];
	for (int i = 0; i < 256; i++)
	{
		dH[i] = 0;
		H0[i] = 0;
		eta[i] = 0;
	}
	H0[17] = H[17];
	H0[34] = H[34];
	H0[68] = H[68];
	H0[136] = H[136];
	Commutator(H0, H, rho, eta);
	Commutator(eta, H, rho, dH);

}

void Run_Flow_EQ()
{
	int dim = 256;
	double x0 = 0.;
	double h = 0.001;
	int N_steps = 1000;
	double y[dim*(N_steps+1)];

	double H[256];
	for (int i = 0; i < 256; i++)
	{
		H[i] = 0.;
	}

	double rho[4] = {1.,1.,0.,0.};

	double t = 1;
	double U = 3.;
	H[17] = -t;
	H[34] = -t;
	H[68] = t;
	H[136] = t;

	H[51] = -U/2.;
	H[85] = -U/2.;
	H[170] = -U/2.;
	H[204] = -U/2.;
	H[195] = U/2.;
	H[165] = U/2.;
	H[90] = U/2.;
	H[60] = U/2.;

	auto start = std::chrono::high_resolution_clock::now();
	rk4_loop(&Wegner_Flow,H,dim,x0,h,N_steps,rho,y);
	auto end = std::chrono::high_resolution_clock::now();
	auto delta = duration_cast<microseconds>(end - start);
	cout << "Time used : " << delta.count() << endl;

	for (int i = 0; i < 256; i++)
	{
		//if (fabs(y[(N_steps-190)*dim+i])>0.0001) cout << i<< "\t"<< __builtin_popcount(i) << "\t" << y[(N_steps-190)*dim+i] << endl;
		if (fabs(y[i+dim*N_steps])>0.0001) cout << i<< "\t"<< __builtin_popcount(i) << "\t" << y[i+dim*N_steps] << endl;
	}
	ofstream file;
	file.open("Wegner_Hubbard.dat");
	for (int n = 0; n <= N_steps; n++){
		file << x0 + n*h;
		for (int i = 0; i < dim; i++)
			file << "\t" << y[n*dim+i];
		file  << "\n";
	}
	file.close();
}

long BinomialCoefficient(int n, int k)
{
	if(n<k) return 0;
	long Bnk = 1;
	for (int m = 1; m <= k; m++)
	{
		Bnk *= n + 1 - m;
		Bnk /= m;
	}
	return Bnk;
}

int WeiSquareBrackets0(int m1, int m2, int m3, int m4, int m5, int m6)
{
	int result = 0;
	int start_sum = std::max({m1+m5+m6, m2+m4+m6, m3+m4+m5, m1+m2+m3});
	int end_sum   = std::min({m1+m2+m4+m5, m1+m3+m4+m6, m2+m3+m5+m6});
	int sign = (start_sum % 2) == 0 ? 1 : -1;
	//cout << "{" << start_sum << "," << end_sum << "}\n";
	for (int n = start_sum; n <= end_sum; n++)
	{
		result += sign 	* BinomialCoefficient(n+1, n-m1-m5-m6)
						* BinomialCoefficient(m1+m5-m6, n-m2-m4-m6)
						* BinomialCoefficient(m1-m5+m6, n-m3-m4-m5)
						* BinomialCoefficient(m5+m6-m1, n-m1-m2-m3);
		sign *= -1;
	}
	return result;
}
long WeiSquareBrackets(int two_m1, int two_m2, int two_m3, int two_m4, int two_m5, int two_m6)
{
	long result = 0; //TODO: check below that start end end are always integers this should be guaranteed if WeiDelta !=0 -> check this
	int start_sum = std::max({two_m1+two_m5+two_m6, two_m2+two_m4+two_m6, two_m3+two_m4+two_m5, two_m1+two_m2+two_m3})/2;
	int end_sum   = std::min({two_m1+two_m2+two_m4+two_m5, two_m1+two_m3+two_m4+two_m6, two_m2+two_m3+two_m5+two_m6})/2;
	long sign = (start_sum % 2) == 0 ? 1 : -1;
	for (int n = start_sum; n <= end_sum; n++)
	{
		result += sign 	* BinomialCoefficient(n+1, n-(two_m1+two_m5+two_m6)/2)
						* BinomialCoefficient((two_m1+two_m5-two_m6)/2, n-(two_m2+two_m4+two_m6)/2)
						* BinomialCoefficient((two_m1-two_m5+two_m6)/2, n-(two_m3+two_m4+two_m5)/2)
						* BinomialCoefficient((two_m5+two_m6-two_m1)/2, n-(two_m1+two_m2+two_m3)/2);
		sign *= -1;
	}
	return result;
}

double WeiDelta(int two_a, int two_b, int two_c)
{
	if ((two_c < std::abs(two_a-two_b)) || (two_c > two_a + two_b) || ((two_a + two_b + two_c)%2 != 0)) return 0;

	double result = 1.;
	result *= Factorial((two_a+two_b-two_c)/2) * Factorial((two_a-two_b+two_c)/2)
			* Factorial((two_b+two_c-two_a)/2)/ Factorial((two_a+two_b+two_c)/2+1); //TODO: use cancellation of factors
	return std::sqrt(result);
}

double Wigner9J(int two_j1, int two_j2, int two_j3, int two_j4, int two_j5, int two_j6, int two_j7, int two_j8, int two_j9)
{
	double prefactor = WeiDelta(two_j1, two_j2, two_j3) * WeiDelta(two_j4, two_j5, two_j6) * WeiDelta(two_j7, two_j8, two_j9)
					* WeiDelta(two_j1, two_j4, two_j7) * WeiDelta(two_j2, two_j5, two_j8) * WeiDelta(two_j3, two_j6, two_j9);
	if (prefactor == 0.) return 0.;

	double result = 0.;
	int start_sum = std::max({std::abs(two_j8 - two_j4),std::abs(two_j2-two_j6),std::abs(two_j1-two_j9)});
	int end_sum   = std::min({two_j8 + two_j4, two_j2 + two_j6, two_j1 + two_j9});
	int sign = (start_sum % 2) == 0 ? 1 : -1;

	for (int two_k = start_sum; two_k <= end_sum; two_k += 2)
	{
		result += (((two_k + 1.) 	* WeiSquareBrackets(two_j1, two_j2, two_j3, two_j6, two_j9, two_k))
										* WeiSquareBrackets(two_j6, two_j4, two_j5, two_j8, two_j2, two_k))
										* WeiSquareBrackets(two_j8, two_j9, two_j7, two_j1, two_j4, two_k);

	}
	return sign * prefactor * result;
}

int main(){
/*
	cout << Coulomb_integral(2, 0.3, 0.5, 0.4, 0.7) << endl;
	cout << Coulomb_integral(1, 0.3, 0.5, 0.4, 0.7) << endl;
	cout << Coulomb_integral(0, 0.3, 0.5, 0.4, 0.7) << endl;
	cout << Coulomb_integral(-1, 0.3, 0.5, 0.4, 0.7) << endl;
	cout << Coulomb_integral(-2, 0.3, 0.5, 0.4, 0.7) << endl;
	cout << Coulomb_integral(-3, 0.3, 0.5, 0.4, 0.7) << endl;
*/
	//cout << Coulomb_integral(-3, 1E-9, 8.34105E-9, 4.01455E-11, 4.80171E-10) << endl;
	//cout << sizeof(bool) << endl << sizeof(double) <<endl;
/*	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 8; j++)
			cout << GM_Product(i, j) << "\t";
		cout << endl;
	}
	cout << endl << (1 << 3) << endl;*/
/*	double rho[4] = {0.6,0.6,0.6,0.6};
	double AB_ord[256];
	for (int i = 0; i < 256; i++)
		AB_ord[i] = 0.;
	NormalOrderedProduct(51, 83, 1, 1., rho, AB_ord);
	for (int i = 0; i < 256; i++)
	{
		if (fabs(AB_ord[i])>0.0001) cout << i << "\t" << AB_ord[i] << endl;
		AB_ord[i]=0;
	}

	NormalOrderedProduct(83, 51, 1, 1., rho, AB_ord);
	for (int i = 0; i < 256; i++)
		if (fabs(AB_ord[i])>0.0001) cout << i << "\t" << AB_ord[i] << endl;*/

/*	double H[256];
	double dH[256];
	double H0[256];
	double eta[256];
	for (int i = 0; i < 256; i++)
	{
		H[i] = 0.;
		dH[i] = 0.;
		H0[i] = 0.;
		eta[i] = 0.;
	}

	double rho[4] = {1.,1.,0.,0.};

	double t = 1;
	double U = 3.;
	H[17] = -t;
	H[34] = -t;
	H[68] = t;
	H[136] = t;

	H0[17] = -t;
	H0[34] = -t;
	H0[68] = t;
	H0[136] = t;

	H[51] = -U/2.;
	H[85] = -U/2.;
	H[170] = -U/2.;
	H[204] = -U/2.;
	H[195] = U/2.;
	H[165] = U/2.;
	H[90] = U/2.;
	H[60] = U/2.;

	Commutator(H0, H, rho, eta);
	Commutator(eta, H, rho, dH);
	for (int i = 0; i < 256; i++)
	{
		if (fabs(dH[i])>0.0001) cout << i<< "\t"<< __builtin_popcount(i) << "\t" << dH[i] << endl;
	}*/
	//Run_Flow_EQ();
	/*
	for (int n = 1; n <= 10; n++)
	{
		for (int k = 0; k <= n; k++)
		{
			cout << BinomialCoefficient(n, k) << "\t";
		}
		cout << endl;
	}
	cout << WeiSquareBrackets0(1, 1, 1, 1, 1, 1) << "\t" << WeiSquareBrackets(2*1, 2*1, 2*1, 2*1, 2*1, 2*1) << endl;
	cout << WeiSquareBrackets0(4, 4, 4, 4, 4, 4) << "\t" << WeiSquareBrackets(2*4, 2*4, 2*4, 2*4, 2*4, 2*4) << endl;
	cout << WeiSquareBrackets0(4, 2, 3, 4, 2, 3) << "\t" << WeiSquareBrackets(2*4, 2*2, 2*3, 2*4, 2*2, 2*3) << endl;
	cout << WeiSquareBrackets0(4, 2, 3, 7, 5, 6) << "\t" << WeiSquareBrackets(2*4, 2*2, 2*3, 2*7, 2*5, 2*6) << endl;
	cout << WeiSquareBrackets0(5, 2, 3, 7, 4, 6) << "\t" << WeiSquareBrackets(2*5, 2*2, 2*3, 2*7, 2*4, 2*6) << endl;
	*/

	cout << Wigner9J(1, 1, 2, 1, 2, 1, 0, 1, 1) << endl;
	cout << Wigner9J(1, 1, 2, 1, 2, 1, 2, 1, 1) << endl;

	cout << Wigner9J(1, 1, 2, 1, 2, 1, 2, 1, 3) << endl;

	cout << Wigner9J(1, 9, 8, 7, 4, 3, 6, 11, 5) << endl;
	cout << Wigner9J(1, 9, 8, 7, 4, 3, 8, 11, 5) << endl;
	/*
	cout << "Check 0 \n";

	cout << Wigner9J(5, 1, 2, 1, 2, 1, 0, 1, 1) << endl;
	cout << Wigner9J(6, 1, 2, 1, 2, 1, 2, 1, 1) << endl;

	cout << Wigner9J(7, 1, 2, 1, 2, 1, 2, 1, 3) << endl;

	cout << Wigner9J(8, 9, 8, 7, 4, 3, 6, 11, 5) << endl;*/
	cout << Wigner9J(9, 9, 8, 7, 4, 3, 8, 11, 5) << endl;
	//cout << WeiSquareBrackets(3, 7, 4, 11, 9, 12) << endl;
	//cout << endl << INT64_MAX << "\t" << 8*sizeof(long long) << "\t" << 8*sizeof(unsigned long)<< "\t" << INT32_MAX << endl;
	long i = -1;
	unsigned int u = 10;
	double x = i * u;
	cout << endl << x << endl;

	return 0;
}


