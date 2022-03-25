/*
 * CPMC_core.cpp
 *
 *  Created on: Mar 23, 2022
 *      Author: martin
 */

namespace QMC
{

void Inverse_matrix(double *A, unsigned int N)
{
 	int IPIV[N];
	unsigned int LWORK = N*N;
	double WORK[LWORK];
	int info;

	dgetrf(&N,&N,A,&N,IPIV,&info);
	dgetri(&N,A,&N,IPIV,WORK,&LWORK,&info);

	if (info !=0) printf("Error in Inverse_matrix, info = %d\n",info);
}

void Initialize_determinants(double **dets,int n_determinants, int n_single_particle_basis, int n_partilces)
{

}
double Determinant(double *A, int *ipiv, unsigned int N)
{
	double result = A[0];
	for (int i=1; i < N; i++)
	{
		if(i+1==ipiv[i])
			result *= A[i*(1+N)];
		else
			result *= -A[i*(1+N)];
	}
	if (ipiv[0]==1)
		return result;
	else
		return -result;
	//TODO: test this

}

void Propagate_free(double *&det, int n_single_particle_basis, int n_partilces, double &olap, double *&inverse_olap_matrix, double *propagator_free, double *trial_det)
{
	const char Transpose = 'T';
	const char NoTranspose = 'N';
	double one  = 1;
	double zero = 0;

	unsigned int M = n_single_particle_basis;	// number of rows of det
	unsigned int N = n_partilces;	// number of cols of det
	unsigned int K = M;	// number of rows of det

	double *C = (double*) malloc(M*N*sizeof(double));

	dgemm(&NoTranspose,&NoTranspose,&M,&N,&K,&one,propagator_free,&K,det,&K,&zero,C,&M);
	free(det);
	det = C;
	PrintMartix(det, M, N);
	dgemm(&Transpose,&NoTranspose,&N,&M,&M,&one,trial_det,&M,det,&M,&zero,inverse_olap_matrix,&N);
	PrintMartix(inverse_olap_matrix, N);
 	int ipiv[N];
	unsigned int LWORK = N*N;
	double WORK[LWORK];
	int info;

	dgetrf(&N,&N,inverse_olap_matrix,&N,ipiv,&info);
	if (info != 0) printf("Error in Propagate_free, dgetrf info = %d",info);
	double olap_new = Determinant(inverse_olap_matrix, ipiv, N);

	dgetri(&N,inverse_olap_matrix,&N,ipiv,WORK,&LWORK,&info);
	if (info != 0) printf("Error in Propagate_free, dgetri info = %d",info);
	//Inverse_matrix(inverse_olap_matrix, N);
	PrintMartix(inverse_olap_matrix, N);

	//TODO: test this
	cout << "Det = " << olap_new << endl;
}

void Propagate_dets(double **dets, double *weights, double *propagator_free, double *olaps, double **inverse_olap_matrix,
		int n_determinants, int n_single_particle_basis, int n_partilces, double *trial_det)
{
	for (int i = 0; i < n_determinants; i++)
	{
		if (weights[i] > 0)
		{
			Propagate_free(dets[i],n_single_particle_basis, n_partilces, olaps[i],inverse_olap_matrix[i],propagator_free,trial_det);
		}
	}
}
void Orthonormalize_dets(double **dets,int n_determinants, int n_single_particle_basis, int n_partilces)
{

}
void Population_control_dets(double **dets,int n_determinants, int n_single_particle_basis, int n_partilces)
{

}
double Mean(double *X, int n)
{
	double result = 0;
	for (int i = 0; i < n; i++)
	{
		result += X[i];
	}
	return result / n;
}
double STD_dev(double *X, int n, double mean)
{
	double result = 0;
	for (int i = 0; i < n; i++)
	{
		result += (X[i]-mean)*(X[i]-mean);
	}
	return sqrt(result/(n-1.));
}
template<typename T>
void Fill_array(T *A, int n, T value)
{
	for (int i = 0; i < n; i++)
	{
		A[i] = value;
	}
}

void Estimate_ground_state_energy_CPMC(double *propagator_free, double *trial_det,
		        					   int n_determinants, int n_single_particle_basis, int n_partilces, int n_equil_steps,
									   int n_measurements, int n_steps_during_measure,
									   int period_ortho, int period_pop_ctrl, int period_energy_adjust,
									   double &energy, double &err_energy)
{
	double *dets[n_determinants];
	double *inverse_olap_matrix[n_determinants];
	double weights[n_determinants];
	double olaps[n_determinants];
	double energy_meas[n_measurements];
	Fill_array(energy_meas, n_measurements, 0.);
	Initialize_determinants(dets,n_determinants, n_single_particle_basis, n_single_particle_basis);

	// Equilibration phase:
	for (int i = 0; i < n_equil_steps; i++)
	{
		Propagate_dets(dets,weights, propagator_free,olaps, inverse_olap_matrix, n_determinants, n_single_particle_basis, n_single_particle_basis,trial_det);
		if ( i% period_ortho == 0)
		{
			Orthonormalize_dets(dets,n_determinants, n_single_particle_basis, n_single_particle_basis);
		}
		if ( i% period_pop_ctrl == 0)
		{
			Population_control_dets(dets, n_determinants, n_single_particle_basis, n_partilces);
		}
	}

	// Measurement phase:
	for (int i = 0; i < n_measurements; i++)
	{
		for (int j = 0; j < n_steps_during_measure; j++)
		{
			Propagate_dets(dets,weights, propagator_free,olaps, inverse_olap_matrix, n_determinants, n_single_particle_basis, n_single_particle_basis,trial_det);
			if ( i% period_ortho == 0)
			{
				Orthonormalize_dets(dets,n_determinants, n_single_particle_basis, n_single_particle_basis);
			}
			if ( i% period_pop_ctrl == 0)
			{
				Population_control_dets(dets, n_determinants, n_single_particle_basis, n_partilces);
			}
			if ( i% period_energy_adjust == 0)
			{
				//TODO: adjust exponent ground-state energy
			}
		}
	}
	energy = Mean(energy_meas,n_measurements);
	err_energy = STD_dev(energy_meas, n_measurements, energy)/sqrt(n_measurements);
}

}
