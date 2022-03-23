/*
 * CPMC_core.cpp
 *
 *  Created on: Mar 23, 2022
 *      Author: martin
 */

namespace QMC
{

void Initialize_determinants(double *dets,int n_determinants, int n_single_particle_basis, int n_partilces)
{

}

void Propagate_dets(double *dets,int n_determinants, int n_single_particle_basis, int n_partilces)
{

}
void Orthonormalize_dets(double *dets,int n_determinants, int n_single_particle_basis, int n_partilces)
{

}
void Population_control_dets(double *dets,int n_determinants, int n_single_particle_basis, int n_partilces)
{

}
double Mean()

void Estimate_ground_state_energy_CPMC(int n_determinants, int n_single_particle_basis, int n_partilces, int n_equil_steps,
									   int n_measurements, int n_steps_during_measure,
									   int period_ortho, int period_pop_ctrl, int period_energy_adjust)
{
	double dets[n_determinants * n_single_particle_basis * n_single_particle_basis];
	Initialize_determinants(dets,n_determinants, n_single_particle_basis, n_single_particle_basis);

	// Equilibration phase:
	for (int i = 0; i < n_equil_steps; i++)
	{
		Propagate_dets(dets,n_determinants, n_single_particle_basis, n_single_particle_basis);
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
			Propagate_dets(dets,n_determinants, n_single_particle_basis, n_single_particle_basis);
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
}

}
