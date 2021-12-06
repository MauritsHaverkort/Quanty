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

double Slater_Integrate_optimize(const uint32_t N_x, const double *x,
												 const double *powx, const int K,
										 	 	 const double *G13,
												 const double *G24)
{
	double w1, w2, gg, g13, g24;
	double val = 0.;
	for (uint32_t i = 1; i < N_x-1; i++)
	{
		w1 = (x[i+1]-x[i-1])/2.;
		g13 = G13[i]; //G1[i]*G3[i] + F1[i]*F3[i];
		g24 = G24[i]; //G2[i]*G4[i] + F2[i]*F4[i];
		val += w1 * w1 * g13 * g24 / x[i];
		w1 = w1 * powx[K*N_x + i];
		for (uint32_t j = i+1; j < N_x-1; j++)
		{
			w2 = (x[j+1]-x[j-1])/2.;
			gg = g13 * G24[j];
			gg +=G13[j] * g24;
			val += w1 * w2 * gg / powx[(K+1)*N_x + j] ;
		}
	}
	return val;
}

void Init_Slater_Function_Pairs(const uint32_t N_x,
								const vector<Interpolation::InterpolatingFunction*> &G,
								const vector<Interpolation::InterpolatingFunction*> &F, double *GFpairs)
{
	uint32_t N_funs = (uint32_t) G.size();
	for (uint32_t i = 0; i < N_funs; i++)
	{
		for (uint32_t j = i; j < N_funs; j++)
		{
			for (uint32_t k = 0; k < N_x; k++)
			{
				GFpairs[(i*N_funs +j)*N_x + k] = G[i]->y[k]*G[j]->y[k]+F[i]->y[k]*F[j]->y[k];
				GFpairs[(j*N_funs +i)*N_x + k] = GFpairs[(i*N_funs +j)*N_x + k];
			}

		}
	}
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

void Calculate_slater_integrals_with_reduced_prefactor(const uint32_t N_x, const double *x,
		 const double *powx, const double *G13, const double *G24,
		 const int twoj1, const int twoj2, const int twoj3, const int twoj4,
		 const int L_start, const int L_end, double *Fk)
{
	for (int L = L_start; L <= L_end; L +=2)
	{
		double c_reduced_11 = sqrt((twoj1+1)*(twoj3+1))*ThreeJSymbolHalf(twoj1, -1, twoj3, 1, 2*L, 0);
		//Cancelled if (((twoj1+1)/2)%2==1) c_reduced_11 *= -1; //TODO: move this to the outermoost possible loop
		double c_reduced_22 = sqrt((twoj2+1)*(twoj4+1))*ThreeJSymbolHalf(twoj2, -1, twoj4, 1, 2*L, 0);
		//Cancelled if (((twoj2+1)/2)%2==1) c_reduced_22 *= -1; //TODO: move this to the outermoost possible loop
		// This is working but slow: Fk1[(L - L1_start)/2] = Slater_Integrate(G[0]->x,powx,L,G[i1]->y,G[i2]->y,G[n1+N_active_orbitals]->y,G[n2+N_active_orbitals]->y,F[i1]->y,F[i2]->y,F[n1+N_active_orbitals]->y,F[n2+N_active_orbitals]->y);
		Fk[(L - L_start)/2] = Slater_Integrate_optimize(N_x, x, &(powx[0]), L, G13, G24);
		Fk[(L - L_start)/2] *=  c_reduced_11 * c_reduced_22;

	}
}


