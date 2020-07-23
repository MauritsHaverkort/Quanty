/*
 * Slater_Integrals.cpp
 *
 *  Created on: Jul 17, 2020
 *      Author: brass
 */


double Slater_Integrate(const int K, const vector<double> &x,
								     const vector<double> &G1,
								     const vector<double> &G2,
								     const vector<double> &G3,
								     const vector<double> &G4)
{
	double w1, w2, gg;
	double val = 0.;
	uint32_t N = (uint32_t) x.size();

	for (uint32_t i = 1; i < N-1; i++)
	{
		w1 = (x[i+1]-x[i-1])/2.;
		val += w1 * w1 * G1[i] * G2[i] * G3[i] * G4[i] / x[i];
		w1 = w1 * pow(x[i],K);
		for (uint32_t j = i+1; j < N-1; j++)
		{
			w2 = (x[j+1]-x[j-1])/2.;
			gg = G1[i] * G2[j] * G3[i] * G4[j] + G1[j] * G2[i] * G3[j] * G4[i];
			val += w1 * w2 / pow(x[j],K+1) * gg;
		}
	}
	return val;
}

double Slater_Integrate(const vector<double> &x, const vector<double> &powx, const int K,
										 	 	 const vector<double> &G1,
												 const vector<double> &G2,
												 const vector<double> &G3,
												 const vector<double> &G4,
												 const vector<double> &F1,
												 const vector<double> &F2,
												 const vector<double> &F3,
												 const vector<double> &F4)
{
	double w1, w2, gg, g13, g24;
	double val = 0.;
	uint32_t N = (uint32_t) x.size();
	for (uint32_t i = 1; i < N-1; i++)
	{
		w1 = (x[i+1]-x[i-1])/2.;
		g13 = G1[i]*G3[i] + F1[i]*F3[i];
		g24 = G2[i]*G4[i] + F2[i]*F4[i];
		val += w1 * w1 * g13 * g24 / x[i];
		w1 = w1 * powx[K*N + i];
		for (uint32_t j = i+1; j < N-1; j++)
		{
			w2 = (x[j+1]-x[j-1])/2.;
			gg = g13 * (G2[j]*G4[j] + F2[j]*F4[j]);
			gg +=(G1[j]*G3[j] + F1[j]*F3[j]) * g24;
			val += w1 * w2 * gg / powx[(K+1)*N + j] ;
		}
	}
	return val;
}

void Init_Slater_Kernel(const vector<int> &list_of_kappas, const vector<double> &x,
							  vector<double> &powx)
{
	int kmax = 1;
	for (uint32_t i = 0; i < (uint32_t) list_of_kappas.size(); i++)
	{
		if (kmax < abs(list_of_kappas[i]))
			kmax = abs(list_of_kappas[i]);
	}
	int Lmax = 2 * kmax;
	uint32_t N = (uint32_t) x.size();
	powx.reserve((Lmax+1)*N);
	for (int L = 0; L <= Lmax; L++)
	{
		for (uint32_t i = 0; i < N ; i++)
		{
			powx[L*N + i] = pow(x[i],L);
		}
	}
}




