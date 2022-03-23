#include "ReducedCoulombmatrixelement.cpp"

void Calculate_Xk(const vector<double> &x, const vector<double> &powx, const int k1, const int k2, const int k3, const int k4,
	 	 const vector<double> &G1,
		 const vector<double> &G2,
		 const vector<double> &G3,
		 const vector<double> &G4,
		 const vector<double> &F1,
		 const vector<double> &F2,
		 const vector<double> &F3,
		 const vector<double> &F4, double *Xk)
{
	int twoj1 = 2*std::abs(k1)-1;
	int twoj2 = 2*std::abs(k2)-1;
	int twoj3 = 2*std::abs(k3)-1;
	int twoj4 = 2*std::abs(k4)-1;

	int K_start = std::max( abs(twoj1-twoj3), abs(twoj2-twoj4))/2;
	int K_end = std::min( twoj1+twoj3, twoj2+twoj4)/2;
	if ( (abs(k1)+abs(k3) + (Sign(k1)+Sign(k3))/2 + K_start ) % 2 == 0)
	{
		K_start++;
	}
	int sign = K_start % 2 == 0 ? 1 : -1;
	for (int K = K_start; K <= K_end; K +=2)
	{
		double c_reduced_13 = sqrt((twoj1+1)*(twoj3+1))*ThreeJSymbolHalf(twoj1, -1, twoj3, 1, 2*K, 0);
		if (((twoj1+1)/2)%2==1) c_reduced_13 *= -1;
		double c_reduced_24 = sqrt((twoj2+1)*(twoj4+1))*ThreeJSymbolHalf(twoj2, -1, twoj4, 1, 2*K, 0);
		if (((twoj2+1)/2)%2==1) c_reduced_24 *= -1;
		Xk[(K - K_start)/2] = Slater_Integrate(x,powx,K,G1,G2,G3,G4,F1,F2,F3,F4);
		Xk[(K - K_start)/2] *= sign * c_reduced_13 * c_reduced_24;
	}
}

void Add_Effective_Hamiltonian_Terms(const vector<vector<uint16_t> > &indices_active, int i1, int i2, int i3, int i4, int twoj1, int twoj2, int twoj3, int twoj4, const double *reduced_matrix_ements, OperatorType &U)
{
	int K_start = std::max( abs(twoj1-twoj3), abs(twoj2-twoj4))/2;
	int K_end = std::min( twoj1+twoj3, twoj2+twoj4)/2;
	for (int K = K_start; K <= K_end; K++)
	{
		for (int two_m1 = -twoj1; two_m1 <= twoj1; two_m1 += 2)
		{
			for (int two_m2 = -twoj2; two_m2 <= twoj2; two_m2 += 2)
			{
				for (int two_m3 = -twoj3; two_m3 <= twoj3; two_m3 += 2)
				{
					int two_m4 = two_m1 + two_m2 - two_m3;
					// note the minus sign below, due to ordering of creation/annihilation operators
					double angular_prefactor = -ThreeJSymbolHalf(twoj1, -two_m1, K, -two_m3+two_m1, twoj3, two_m3)*ThreeJSymbolHalf(twoj2, -two_m2, K, two_m3-two_m1, twoj4, two_m4);
					if (((two_m2-two_m4)/2)%2!=0) angular_prefactor *= -1.;
					if (Abs(angular_prefactor) > DBL_EPSILON)
					{
						uint16_t oadd[4];
						oadd[3] = CrAnMask + indices_active[i1][(twoj1+two_m1)/2];
						oadd[2] = CrAnMask + indices_active[i2][(twoj2+two_m2)/2];
						oadd[1] =            indices_active[i3][(twoj3+two_m3)/2];
						oadd[0] =            indices_active[i4][(twoj4+two_m4)/2];
						RealOperatorAddLaderNormalOrder(&U, oadd, 4, angular_prefactor * reduced_matrix_ements[K-K_start], &QDetNotOrderedFcFaBcBa);
					}
				}
			}
		}
	}
}

void Create_Effective_Hamiltonian_One_Orbital(int i1, int twoj1, int k1, double e1,int N_virtual_orbitals, int N_active_orbitals, int N_funs, const vector<double> &energies,const vector<int> &list_of_kappas_active,
		  const vector<int> &list_of_kappas_virtual, double *x, int N_x, const vector<double> &powx, double *GFpairs,const vector<vector<uint16_t> > &indices_active,OperatorType &U)
{
	int K_start = 0;
	int K_end = twoj1;

	int i2 = i1;
	int i3 = i1;
	int i4 = i1;
	int twoj2 = twoj1;
	int twoj3 = twoj1;
	int twoj4 = twoj1;
	int k3 = k1;

	double reduced_matrix_ements[K_end-K_start +1];
	for (int K = K_start; K <= K_end; K++)
	{
		reduced_matrix_ements[K_end-K_start +1] = 0.;
	}
	for (uint32_t n1 = 0; n1 < N_virtual_orbitals; n1++)
	{
		int k1v = list_of_kappas_virtual[n1];
		int twoj1v = 2*std::abs(k1v)-1;
		double e1v = energies[n1+N_active_orbitals];
		for (uint32_t n2 = 0; n2 < N_virtual_orbitals; n2++)
		{
			int k2v = list_of_kappas_virtual[n2];
			int twoj2v = 2*std::abs(k2v)-1;
			double e2v = energies[n2+N_active_orbitals];

			double energy_factor = 2/(e1v + e2v - 2*e1);

			int L1_start = std::max( abs(twoj1-twoj1v), abs(twoj2-twoj2v))/2;
			int L1_end = std::min( twoj1+twoj1v, twoj2+twoj2v)/2;

			int L2_start = std::max(  abs(twoj3-twoj1v), abs(twoj4-twoj2v))/2;
			int L2_end = std::min( twoj3+twoj1v, twoj4+twoj2v )/2;

			if ( (abs(k1)+abs(k1v) + (Sign(k1)+Sign(k1v))/2 + L1_start ) % 2 == 0)
			{
				L1_start++;
			}
			if ( (abs(k3)+abs(k1v) + (Sign(k3)+Sign(k1v))/2 + L2_start ) % 2 == 0)
			{
				L2_start++;
			}
			if ((L1_start <= L1_end) && (L2_start <= L2_end))
			{
				double Fk1[(L1_end-L1_start)/2 +1];
				//int sign = L1_start % 2 == 0 ? 1 : -1; Cancelled in the end
				Calculate_slater_integrals_with_reduced_prefactor(N_x, x, &(powx[0]), &(GFpairs[(i1*N_funs+(n1+N_active_orbitals))*N_x]), &(GFpairs[(i2*N_funs + n2+N_active_orbitals)*N_x]), twoj1, twoj2, twoj1v, twoj2v, L1_start, L1_end, Fk1);
				double Fk2[(L2_end-L2_start)/2 +1];
				//sign = L2_start % 2 == 0 ? 1 : -1; Cancelled  in the end
				Calculate_slater_integrals_with_reduced_prefactor(N_x, x, &(powx[0]), &(GFpairs[((n1+N_active_orbitals)*N_funs+i3)*N_x]), &(GFpairs[((n2+N_active_orbitals)*N_funs + i4)*N_x]), twoj1v, twoj2v, twoj3, twoj4, L2_start, L2_end, Fk2);
				for (int K = K_start; K <= K_end; K++)
				{
					for (int L2 = L2_start; L2 <= L2_end; L2 +=2)
					{
						for (int L1 = L1_start; L1 <= L1_end; L1 +=2)
						{
							double angular_prefactor = (2*K+1)*Wigner9J(twoj2, twoj2v, 2*L1, twoj3, 2*L2, twoj1v, 2*K, twoj4, twoj1)/4.;//  1/4 is the prefactor according to Effective_Hamiltonian_2.pdf
							reduced_matrix_ements[K-K_start] += angular_prefactor * Fk1[(L1 - L1_start)/2] * Fk2[(L2 - L2_start)/2]*energy_factor;
						}
					}
				}
			}
		}
		// three active, one virtual
		for (uint32_t n2 = 0; n2 < N_active_orbitals; n2++)
		{
			int k2v = list_of_kappas_active[n2];
			int twoj2v = 2*std::abs(k2v)-1;
			double e2v = energies[n2];

			double energy_factor = 2/(e1v + e2v - 2*e1);

			int L1_start = std::max( abs(twoj1-twoj1v), abs(twoj2-twoj2v))/2;
			int L1_end = std::min( twoj1+twoj1v, twoj2+twoj2v)/2;

			int L2_start = std::max(  abs(twoj3-twoj1v), abs(twoj4-twoj2v))/2;
			int L2_end = std::min( twoj3+twoj1v, twoj4+twoj2v )/2;

			if ( (abs(k1)+abs(k1v) + (Sign(k1)+Sign(k1v))/2 + L1_start ) % 2 == 0)
			{
				L1_start++;
			}
			if ( (abs(k3)+abs(k1v) + (Sign(k3)+Sign(k1v))/2 + L2_start ) % 2 == 0)
			{
				L2_start++;
			}
			if ((L1_start <= L1_end) && (L2_start <= L2_end))
			{
				double Fk1[(L1_end-L1_start)/2 +1];
				//int sign = L1_start % 2 == 0 ? 1 : -1; Cancelled in the end
				Calculate_slater_integrals_with_reduced_prefactor(N_x, x, &(powx[0]), &(GFpairs[(i1*N_funs+(n1+N_active_orbitals))*N_x]), &(GFpairs[(i2*N_funs + n2)*N_x]), twoj1, twoj2, twoj1v, twoj2v, L1_start, L1_end, Fk1);

				double Fk2[(L2_end-L2_start)/2 +1];
				//sign = L2_start % 2 == 0 ? 1 : -1; Cancelled  in the end
				Calculate_slater_integrals_with_reduced_prefactor(N_x, x, &(powx[0]), &(GFpairs[((n1+N_active_orbitals)*N_funs+i3)*N_x]), &(GFpairs[((n2)*N_funs + i4)*N_x]), twoj1v, twoj2v, twoj3, twoj4, L2_start, L2_end, Fk2);
				for (int K = K_start; K <= K_end; K++)
				{
					for (int L2 = L2_start; L2 <= L2_end; L2 +=2)
					{
						for (int L1 = L1_start; L1 <= L1_end; L1 +=2)
						{
							double angular_prefactor = (2*K+1)*Wigner9J(twoj2, twoj2v, 2*L1, twoj3, 2*L2, twoj1v, 2*K, twoj4, twoj1)/2.; //  1/2 is the prefactor according to Effective_Hamiltonian_2.pdf
							reduced_matrix_ements[K-K_start] += angular_prefactor * Fk1[(L1 - L1_start)/2] * Fk2[(L2 - L2_start)/2]* energy_factor;
						}
					}
				}
			} // if L_start
		} // n2 loop
	}// n1 loop
	Add_Effective_Hamiltonian_Terms(indices_active, i1, i2, i3, i4, twoj1, twoj2, twoj3, twoj4, reduced_matrix_ements, U);

}

void Create_Effective_Hamiltonian_Six_Orbitals(const double Fk1[], int L1_start, int L1_end,const double Fk2[], int L2_start, int L2_end, double energy_factor,
		int i1, int i2, int i3, int i4, int twoj1,int twoj2, int twoj3, int twoj4, int twoj1v, int twoj2v, int K_start, int K_end,
		const vector<vector<uint16_t> > &indices_active, OperatorType &U)
{

	//cout << "K_start = " << K_start << "\t" << K_end <<endl;
	for (int K = K_start; K <= K_end; K++)
	{
		double reduced_matrix_ement = 0;
		for (int L2 = L2_start; L2 <= L2_end; L2 +=2)
		{
			for (int L1 = L1_start; L1 <= L1_end; L1 +=2)
			{
				double angular_prefactor = (2*K+1)*Wigner9J(twoj2, twoj2v, 2*L1, twoj3, 2*L2, twoj1v, 2*K, twoj4, twoj1);//  1/4 is the prefactor according to Effective_Hamiltonian_2.pdf
				reduced_matrix_ement += angular_prefactor * Fk1[(L1 - L1_start)/2] * Fk2[(L2 - L2_start)/2]*energy_factor; // it must be moved to the energy factor, because different for different terms
				//cout << "I'm here\t" << angular_prefactor << endl;
			}
		}
		for (int two_m1 = -twoj1; two_m1 <= twoj1; two_m1 += 2)
		{
			for (int two_m2 = -twoj2; two_m2 <= twoj2; two_m2 += 2)
			{
				for (int two_m3 = -twoj3; two_m3 <= twoj3; two_m3 += 2)
				{
					int two_m4 = two_m1 + two_m2 - two_m3;
					// note the minus sign below, due to ordering of creation/annihilation operators
					//double angular_prefactor = -ThreeJSymbolHalf(twoj1, -two_m1, 2*K, -two_m3+two_m1, twoj3, two_m3)*ThreeJSymbolHalf(twoj2, -two_m2, 2*K, two_m3-two_m1, twoj4, two_m4);
					double angular_prefactor = -ThreeJSymbolHalf(twoj1, -two_m1, 2*K, -two_m4+two_m1, twoj4, two_m4)*ThreeJSymbolHalf(twoj2, -two_m2, 2*K, two_m4-two_m1, twoj3, two_m3);
					//if (((two_m2-two_m4)/2)%2==1) angular_prefactor *= -1.;
					if (((two_m2-two_m4)/2)%2!=0) angular_prefactor *= -1.;
					if (Abs(angular_prefactor) > DBL_EPSILON)
					{
						//cout << "Adding\t" << angular_prefactor << "\t" << reduced_matrix_ement << endl;
						uint16_t oadd[4];
						oadd[3] = CrAnMask + indices_active[i1][(twoj1+two_m1)/2];
						oadd[2] = CrAnMask + indices_active[i2][(twoj2+two_m2)/2];
						oadd[1] =            indices_active[i3][(twoj3+two_m3)/2];
						oadd[0] =            indices_active[i4][(twoj4+two_m4)/2];
						//cout << oadd[3] - CrAnMask << "\t" << oadd[2] - CrAnMask << "\t" << oadd[1] << "\t" << oadd[0] << "\tmod : " << ((two_m2-two_m4)/2)%2 << endl;
						RealOperatorAddLaderNormalOrder(&U, oadd, 4, angular_prefactor * reduced_matrix_ement, &QDetNotOrderedFcFaBcBa);
					}
				}
			}
		}
	}

}
//void Create_Effective_Hamiltonian_Six_Orbitals_rework(const double Fk1[], int L1_start, int L1_end,const double Fk2[], int L2_start, int L2_end, double energy_factor,
//		int i1, int i2, int i3, int i4, int twoj1,int twoj2, int twoj3, int twoj4, int twoj1v, int twoj2v, const vector<vector<uint16_t> > &indices_active, OperatorType &U)
//{
//
//	//cout << "K_start = " << K_start << "\t" << K_end <<endl;
//	for (int K = K_start; K <= K_end; K++)
//	{
//
//
//		for (int two_m1 = -twoj1; two_m1 <= twoj1; two_m1 += 2)
//		{
//			for (int two_m2 = -twoj2; two_m2 <= twoj2; two_m2 += 2)
//			{
//				for (int two_m3 = -twoj3; two_m3 <= twoj3; two_m3 += 2)
//				{
//					int two_m4 = two_m1 + two_m2 - two_m3;
//					double reduced_matrix_ement = 0;
//
//					for (int L2 = L2_start; L2 <= L2_end; L2 +=2)
//					{
//						for (int L1 = L1_start; L1 <= L1_end; L1 +=2)
//						{
//							double angular_prefactor = (2*K+1)*Wigner9J(twoj2, twoj2v, 2*L1, twoj3, 2*L2, twoj1v, 2*K, twoj4, twoj1);//  1/4 is the prefactor according to Effective_Hamiltonian_2.pdf
//							reduced_matrix_ement += angular_prefactor * Fk1[(L1 - L1_start)/2] * Fk2[(L2 - L2_start)/2]*energy_factor; // it must be moved to the energy factor, because different for different terms
//						}
//					}
//
//					double angular_prefactor = -ThreeJSymbolHalf(twoj1, -two_m1, 2*K, -two_m4+two_m1, twoj4, two_m4)*ThreeJSymbolHalf(twoj2, -two_m2, 2*K, two_m4-two_m1, twoj3, two_m3);
//					if (((two_m2-two_m4)/2)%2==1) angular_prefactor *= -1.;
//					if (Abs(angular_prefactor) > DBL_EPSILON)
//					{
//						uint16_t oadd[4];
//						oadd[3] = CrAnMask + indices_active[i1][(twoj1+two_m1)/2];
//						oadd[2] = CrAnMask + indices_active[i2][(twoj2+two_m2)/2];
//						oadd[1] =            indices_active[i3][(twoj3+two_m3)/2];
//						oadd[0] =            indices_active[i4][(twoj4+two_m4)/2];
//						RealOperatorAddLaderNormalOrder(&U, oadd, 4, angular_prefactor * reduced_matrix_ement, &QDetNotOrderedFcFaBcBa);
//					}
//				}
//			}
//		}
//	}
//
//}

double Average_Coulomb_energy_difference_of_excitation(const double *Q_av, const uint32_t N_funs, const vector<double> &occupations, int ind[], int delta_occ[], int N_excitations, const vector<int> &list_of_kappas)
{
	double result = 0;
	int N_occ = (int) occupations.size();
	for (int a = 0; a < N_excitations; a++)
	{
		int i = ind[a];
		double twoj1 = (double) 2*std::abs(list_of_kappas[i])-1;
		int d = delta_occ[a];
		if (i < N_occ)
			result += (2*occupations[i]-1+d)*d * Q_av[i*N_funs + i]*(twoj1+1.)/twoj1/2.;
		else
			result += (d-1)*d * Q_av[i*N_funs + i]*(twoj1+1.)/twoj1/2.;
		for (int b = a+1; b < N_excitations; b++)
		{
			int j = ind[b];
			result += d * delta_occ[b] * Q_av[i*N_funs + j];
		}
		for (int j = i+1; j < N_occ; j++)
		{
			result += occupations[j] * d * Q_av[i*N_funs + j];
		}
		for (int j = 0; j < std::min(i,N_occ); j++) //TODO: wurste diese beiden loops zusammen
		{
			result += occupations[j] * d * Q_av[i*N_funs + j];
		}
	}


	return result;
}

int Determine_excitations(int n1, int n2, int i3, int i4, int ind[], int delta_occ[])
{
	if (n1 == i3)
	{
		ind[1] = n2;
		ind[0] = i4;
		delta_occ[1] = 1;
		delta_occ[0] = -1;
		return 2;
	}
	else if (n1 == i4)
	{
		ind[1] = n2;
		ind[0] = i3;
		delta_occ[1] = 1;
		delta_occ[0] = -1;
		return 2;
	}
	else if (i3 == i4)
	{
		ind[0] = i3;
		delta_occ[0] = -2;
		if (n1 == n2)
		{
			ind[1] = n1;
			delta_occ[1] = 2;
			return 2;
		}
		else
		{
			ind[1] = n1;
			ind[2] = n2;
			delta_occ[1] = 1;
			delta_occ[2] = 1;
			return 3;
		}

	}
	else
	{
		ind[0] = i3;
		delta_occ[0] = -1;
		ind[1] = i4;
		delta_occ[1] = -1;
		if (n1 == n2)
		{
			ind[2] = n1;
			delta_occ[2] = 2;
			return 3;
		}
		else
		{
			ind[2] = n1;
			delta_occ[2] = 1;
			ind[3] = n2;
			delta_occ[3] = 1;
			return 4;
		}

	}
}

//double Coulomb_energy_difference_of_excitation(const double *Q_av, const uint32_t N_funs, const vector<double> &occupations, int ind[], int N_excitations)
//{
//	double result = 0.;
//	int N_occ = (int) occupations.size();
//	for (int a = 0; a < N_excitations; a++)
//	{
//		if (N_excitations == 1)
//		{
//			result -= Q_av[ind[a]*N_funs + ind[1]];
//		}
//		else
//		{
//			for (int b = 0; b < N_excitations; b++)
//			{
//				result += Q_av[ind[a]*N_funs + ind[b]]/2.;
//				result -= Q_av[ind[a]*N_funs + ind[b+N_excitations]]/2.;
//			}
//			for (int b = 0; b < N_excitations; b++)
//			{
//				result -= Q_av[ind[a+N_excitations]*N_funs + ind[b]]/2.;
//				result += Q_av[ind[a+N_excitations]*N_funs + ind[b+N_excitations]]/2.;
//			}
//		}
//		for (int i = 0; i < N_occ; i++)
//		{
//			result += occupations[i] * Q_av[ind[a]*N_funs + i];
//			result -= occupations[i] * Q_av[ind[a+N_excitations]*N_funs + i];
//		}
//
//	}
//
//	return result;
//} dis is bs

void Create_Effective_Hamiltonian(const vector<Interpolation::InterpolatingFunction*> &G,
										  const vector<Interpolation::InterpolatingFunction*> &F,
										  vector<double> &energies,
										  const vector<vector<uint16_t> > &indices_active,
										  const vector<int> &list_of_kappas_active,
										  const vector<int> &list_of_kappas_virtual,
										  const vector<double> &occupations,
										  const uint32_t NF,
										  OperatorType &U)
{
	uint32_t N_active_orbitals  = (uint32_t) indices_active.size();
	uint32_t N_virtual_orbitals = (uint32_t) list_of_kappas_virtual.size();
	uint32_t N_occupied_shells = (uint32_t) occupations.size();
	uint32_t N_funs = N_active_orbitals + N_virtual_orbitals;
	uint32_t N_x = (uint32_t) G[0]->x.size();

	InitNamedOperator(&U, "Effective H");
	U.NFermion = NF;
	U.NWord = Min(U.NFermion, 1 + (U.NFermion - 1) / 8);
	InitOperatorFixedLength(&U, 2, 8);
	InitOperatorFixedLength(&U, 4, 8);
	vector<double> powx;
	vector<int> list_of_kappas;
	list_of_kappas.reserve(N_active_orbitals+N_virtual_orbitals);
	list_of_kappas.insert(list_of_kappas.end(), list_of_kappas_active.begin(), list_of_kappas_active.end());
	list_of_kappas.insert(list_of_kappas.end(), list_of_kappas_virtual.begin(), list_of_kappas_virtual.end());
	Init_Slater_Kernel(list_of_kappas,G[0]->x,powx);

	double GFpairs[N_funs*N_funs*N_x];
	double *x = &(G[0]->x[0]);
	Init_Slater_Function_Pairs(N_x, G, F, GFpairs);

	double Q_av[N_funs*N_funs];
	#pragma omp parallel for
	for (int i = 0; i < N_funs; i++)
	{
		int k1 = list_of_kappas[i];
		int twoj1 = 2*std::abs(k1)-1;
		Q_av[i*N_funs + i] = Average_Coulomb_energy_two_shells(N_x, x, &(powx[0]), &(GFpairs[(i*N_funs+i)*N_x]), &(GFpairs[(i*N_funs+i)*N_x]), &(GFpairs[(i*N_funs+i)*N_x]), twoj1, twoj1, k1, k1);

		for (int j = i+1; j < N_funs; j++)
		{
			int k2 = list_of_kappas[j];
			int twoj2 = 2*std::abs(k2)-1;
			Q_av[i*N_funs + j] = Average_Coulomb_energy_two_shells(N_x, x, &(powx[0]), &(GFpairs[(i*N_funs+i)*N_x]), &(GFpairs[(j*N_funs+j)*N_x]), &(GFpairs[(i*N_funs+j)*N_x]), twoj1, twoj2, k1, k2);
			Q_av[j*N_funs + i] = Q_av[i*N_funs + j];
		}
	}
//	cout << "En: " << energies[0] << "\t" << energies[1] << endl;
//	cout << "En-diff: " <<  2*(energies[1]) - 2*(energies[0]) + 2*Q_av[1*N_funs + 1] - 2*Q_av[0*N_funs + 0] << "\n";
//	cout << "En-diff: " <<  (energies[1]+energies[0]) - 2*(energies[0]) + Q_av[1*N_funs + 0] - 2*Q_av[0*N_funs + 0]<< "\n";
//	cout << "En-diff wrong: " <<  (energies[1]+energies[0]) - 2*(energies[0]) + Q_av[1*N_funs + 0] - 2*Q_av[0*N_funs + 0]+Q_av[0*N_funs + 0]/2.+Q_av[1*N_funs + 1]/2.<< "\n";
	//#pragma omp parallel for
//	for (int i = 0; i < N_funs; i++)
//	{
//		for (int j = 0; j < N_occupied_shells; j++)
//		{
//			energies[i] += occupations[j]*Q_av[i*N_funs+j];
//		}
////		if (i==0)
////			cout << energies[0] -Q_av[0] << endl;
////		else
////			cout << energies[i] << endl;
//	}

	for (uint32_t n1 = 0; n1 < N_funs; n1++)
	{
		int k1v = list_of_kappas[n1];
		int twoj1v = 2*std::abs(k1v)-1;
		int l1v = (twoj1v + Sgn(k1v))/2;
		double e1v = energies[n1];

		for (uint32_t n2 = 0; n2 < N_virtual_orbitals; n2++)
		{
			int k2v = list_of_kappas_virtual[n2];
			int twoj2v = 2*std::abs(k2v)-1;
			int l2v = (twoj2v + Sgn(k2v))/2;
			double e2v = energies[n2+N_active_orbitals];// + 2*Q_av[(n1)*N_funs + n2+N_active_orbitals];

			Reduced_Coulomb_matrixelement X_mat[N_active_orbitals*N_active_orbitals];
			#pragma omp parallel for
			for (uint32_t i1 = 0; i1 < N_active_orbitals; i1++)
			{
				int k1 = list_of_kappas_active[i1];
				int twoj1 = 2*std::abs(k1)-1;
				int l1 = (twoj1 + Sgn(k1))/2;

				for (uint32_t i2 = 0; i2 < N_active_orbitals; i2++)
				{
					int k2 = list_of_kappas_active[i2];
					int twoj2 = 2*std::abs(k2)-1;
					int l2 = (twoj2 + Sgn(k2))/2;

					int L1_start = std::max( abs(twoj1-twoj1v), abs(twoj2-twoj2v))/2;
					int L1_end = std::min( twoj1+twoj1v, twoj2+twoj2v)/2;
					if ( (l1 + l1v + L1_start ) % 2 == 1)
					{
						L1_start++;
					}
					if ( (l2 + l2v + L1_start ) % 2 == 1)
					{
						X_mat[i1*N_active_orbitals + i2].L_start = 0;
						X_mat[i1*N_active_orbitals + i2].L_end = -1;
					}
					else
					{
						X_mat[i1*N_active_orbitals + i2] = Reduced_Coulomb_matrixelement(L1_start,L1_end);

						//int sign = L1_start % 2 == 0 ? 1 : -1; Cancelled in the end
						// TODO: rewrite function below as X_mat contains info about L_start / L__end
						Calculate_slater_integrals_with_reduced_prefactor(N_x, x, &(powx[0]), &(GFpairs[(i1*N_funs+(n1))*N_x]), &(GFpairs[(i2*N_funs + n2+N_active_orbitals)*N_x]), twoj1, twoj2, twoj1v, twoj2v, L1_start, L1_end, X_mat[i1*N_active_orbitals + i2].Fk);
					}
				}
			}
			for (uint32_t i1 = 0; i1 < N_active_orbitals; i1++)
			{
				int k1 = list_of_kappas_active[i1];
				int twoj1 = 2*std::abs(k1)-1;
				double e1 = energies[i1];
				//e1 += Q_av[(n1)*N_funs + i1] + Q_av[(n2+N_active_orbitals)*N_funs + i1]; // plus since e1 gets subtracted, hence opposite sign

				for (uint32_t i2 = 0; i2 < N_active_orbitals; i2++)
				{
					int k2 = list_of_kappas_active[i2];
					int twoj2 = 2*std::abs(k2)-1;
					double e2 = energies[i2];
					//e2 += Q_av[(n1)*N_funs + i2] + Q_av[(n2+N_active_orbitals)*N_funs + i2] - 2*Q_av[i1*N_funs + i2];


					int L1_start = X_mat[i1*N_active_orbitals + i2].L_start;
					int L1_end = X_mat[i1*N_active_orbitals + i2].L_end;

					double *Fk1 = X_mat[i1*N_active_orbitals + i2].Fk;

					for (uint32_t i3 = 0; i3 < N_active_orbitals; i3++)
					{
						int k3 = list_of_kappas_active[i3];
						int twoj3 = 2*std::abs(k3)-1;
						double e3 = energies[i3];
						//e3 += Q_av[(n1)*N_funs + i3] + Q_av[(n2+N_active_orbitals)*N_funs + i3];
						for (uint32_t i4 = 0; i4 < N_active_orbitals; i4++)
						{
							int k4 = list_of_kappas_active[i4];
							int twoj4 = 2*std::abs(k4)-1;
							double e4 = energies[i4];
							//e4 += Q_av[(n1)*N_funs + i4] + Q_av[(n2+N_active_orbitals)*N_funs + i4] - 2*Q_av[i3*N_funs + i4];

							int K_start = std::max( abs(twoj1-twoj3), abs(twoj2-twoj4))/2;
							int K_end = std::min( twoj1+twoj3, twoj2+twoj4)/2;
							if ( K_start<= K_end)
							{

								double denom1, denom2;
								int ind[4] ;
								int delta_occ[4];
								int N_excitations = Determine_excitations(n1, n2+N_active_orbitals, i3, i4, ind, delta_occ);
								denom1 = e1v + e2v - e3 - e4 + Average_Coulomb_energy_difference_of_excitation(Q_av, N_funs, occupations, ind, delta_occ, N_excitations, list_of_kappas);
								N_excitations = Determine_excitations(n1, n2+N_active_orbitals, i1, i2, ind, delta_occ);
								denom2 = e1v + e2v - e1 - e2 + Average_Coulomb_energy_difference_of_excitation(Q_av, N_funs, occupations, ind, delta_occ, N_excitations, list_of_kappas);


								double energy_factor = 1/denom1 + 1/denom2;
								//double energy_factor = 1/(e1v + e2v - (e3 + e4 ));
								//energy_factor = 2.; // TODO: remove that
								//cout << n1 << ", "<< n2+N_active_orbitals << ", " << i1 << ", "<< i2 << ": " << denom2 << endl;
								//cout << n1 << ", "<< n2+N_active_orbitals << ", " << i3 << ", "<< i4 << ": " << denom1 << endl;
								//cout << "-------------------------\n";

								int L2_start = X_mat[i3*N_active_orbitals + i4].L_start;
								int L2_end = X_mat[i3*N_active_orbitals + i4].L_end;

								double *Fk2 = X_mat[i3*N_active_orbitals + i4].Fk;

								if ((L1_start <= L1_end) && (L2_start <= L2_end))
								{
									if (n1 < N_active_orbitals)
										Create_Effective_Hamiltonian_Six_Orbitals(Fk1, L1_start, L1_end, Fk2, L2_start, L2_end, energy_factor/2., i1, i2, i3, i4, twoj1, twoj2, twoj3, twoj4, twoj1v, twoj2v, K_start, K_end, indices_active, U);
									else
										Create_Effective_Hamiltonian_Six_Orbitals(Fk1, L1_start, L1_end, Fk2, L2_start, L2_end, energy_factor/4., i1, i2, i3, i4, twoj1, twoj2, twoj3, twoj4, twoj1v, twoj2v, K_start, K_end, indices_active, U);
									//cout << "Fk1: " << Fk1[0] << "\t" << "Fk2: " << Fk2[0] << endl;
								}
							}
						}
					}
				} // i2
			} // i1
		}// n2
				// three active, one virtual
//		for (uint32_t n2 = 0; n2 < N_active_orbitals; n2++)
//		{
//			int k2v = list_of_kappas_active[n2];
//			int twoj2v = 2*std::abs(k2v)-1;
//			double e2v = energies[n2] + Q_av[(n1+N_active_orbitals)*N_funs + n2];
//			Reduced_Coulomb_matrixelement X_mat[N_active_orbitals*N_active_orbitals];
//			#pragma omp parallel for
//			for (uint32_t i1 = 0; i1 < N_active_orbitals; i1++)
//			{
//				int k1 = list_of_kappas_active[i1];
//				int twoj1 = 2*std::abs(k1)-1;
//
//				for (uint32_t i2 = 0; i2 < N_active_orbitals; i2++)
//				{
//					int k2 = list_of_kappas_active[i2];
//					int twoj2 = 2*std::abs(k2)-1;
//
//					int L1_start = std::max( abs(twoj1-twoj1v), abs(twoj2-twoj2v))/2;
//					int L1_end = std::min( twoj1+twoj1v, twoj2+twoj2v)/2;
//					if ( (abs(k1)+abs(k1v) + (Sign(k1)+Sign(k1v))/2 + L1_start ) % 2 == 0)
//					{
//						L1_start++;
//					}
//					X_mat[i1*N_active_orbitals + i2] = Reduced_Coulomb_matrixelement(L1_start,L1_end);
//					//int sign = L1_start % 2 == 0 ? 1 : -1; Cancelled in the end
//					Calculate_slater_integrals_with_reduced_prefactor(N_x, x, &(powx[0]), &(GFpairs[(i1*N_funs+(n1+N_active_orbitals))*N_x]), &(GFpairs[(i2*N_funs + n2)*N_x]), twoj1, twoj2, twoj1v, twoj2v, L1_start, L1_end, X_mat[i1*N_active_orbitals + i2].Fk);
//				}
//			}
//			for (uint32_t i1 = 0; i1 < N_active_orbitals; i1++)
//			{
//				int k1 = list_of_kappas_active[i1];
//				int twoj1 = 2*std::abs(k1)-1;
//				double e1 = energies[i1];
//				e1 += Q_av[(n1+N_active_orbitals)*N_funs + i1] + Q_av[(n2)*N_funs + i1];
//
//				for (uint32_t i2 = 0; i2 < N_active_orbitals; i2++)
//				{
//					int k2 = list_of_kappas_active[i2];
//					int twoj2 = 2*std::abs(k2)-1;
//					double e2 = energies[i2] + Q_av[(n1+N_active_orbitals)*N_funs + i2] + Q_av[(n2)*N_funs + i2] - Q_av[i1*N_funs + i2];
//
//					int L1_start = X_mat[i1*N_active_orbitals + i2].L_start;
//					int L1_end = X_mat[i1*N_active_orbitals + i2].L_end;
//
//					double *Fk1 = X_mat[i1*N_active_orbitals + i2].Fk;
//
//					for (uint32_t i3 = 0; i3 < N_active_orbitals; i3++)
//					{
//						int k3 = list_of_kappas_active[i3];
//						int twoj3 = 2*std::abs(k3)-1;
//						double e3 = energies[i3];
//						e3 += Q_av[(n1+N_active_orbitals)*N_funs + i3] + Q_av[(n2)*N_funs + i3];
//						for (uint32_t i4 = 0; i4 < N_active_orbitals; i4++)
//						{
//							int k4 = list_of_kappas_active[i4];
//							int twoj4 = 2*std::abs(k4)-1;
//							double e4 = energies[i4] + Q_av[(n1+N_active_orbitals)*N_funs + i4] + Q_av[(n2)*N_funs + i4] - Q_av[i3*N_funs + i4];
//
//							int K_start = std::max( abs(twoj1-twoj3), abs(twoj2-twoj4))/2;
//							int K_end = std::min( twoj1+twoj3, twoj2+twoj4)/2;
//							if ( K_start<= K_end)
//							{
//								double energy_factor = 1/(e1v + e2v - e3 - e4) + 1/(e1v + e2v - e1 - e2);
//								//double energy_factor = 1/(e1v + e2v - (e3 + e4 ));
//								energy_factor = 1.; // TODO: remove that
//
//
//								int L2_start = X_mat[i3*N_active_orbitals + i4].L_start;
//								int L2_end = X_mat[i3*N_active_orbitals + i4].L_end;
//
//								double *Fk2 = X_mat[i3*N_active_orbitals + i4].Fk;
//
//								if ((L1_start <= L1_end) && (L2_start <= L2_end))
//								{
//									Create_Effective_Hamiltonian_Six_Orbitals(Fk1, L1_start, L1_end, Fk2, L2_start, L2_end, energy_factor/2., i1, i2, i3, i4, twoj1, twoj2, twoj3, twoj4, twoj1v, twoj2v, K_start, K_end, indices_active, U);
//									//cout << "Fk1: " << Fk1[0] << "\t" << "Fk2: " << Fk2[0] << endl;
//								} // if L_start
//							} // n2 loop
//						}// if K_start
//					}// i4
//				}// i3 loop
//			}// N1 loop
//		}// i2 loop
	}// i1 loop
}



static int LuaCreateEffectiveHamiltonian(lua_State *L)
{
	vector<Interpolation::InterpolatingFunction*> G;
	vector<Interpolation::InterpolatingFunction*> F;
	vector<double> energies;
	vector<vector<uint16_t> > indices_active;
	vector<int> list_of_kappas_active;
	vector<int> list_of_kappas_virtual;
	vector<double> occupations;
	uint32_t NF;

	lua_check(L,1,NF);
	lua_check(L,2,indices_active);
	lua_getfield(L,2,"kappas");
	if (lua_isnil(L,-1))
	{
		lua_pop(L,1);
		return luaL_error(L,"Error in LuaCreate_Relativistic_Coulomb_Operator: Argument 2: indices.kappas does not exist.\n");
	}
	else
	{
		lua_check(L,-1,list_of_kappas_active);
		lua_pop(L,1);
		if (indices_active.size() != list_of_kappas_active.size()) return luaL_error(L,"Error in LuaCreate_Relativistic_Coulomb_Operator: Argument 2 : indices.kappas and indices do not have same size.\n");
	}
	lua_check(L,3,list_of_kappas_virtual);
	lua_check(L,4,G);
	lua_check(L,5,F);
	lua_check(L,6,energies);
	lua_check(L,7,occupations);

	OperatorType *U = (OperatorType *) lua_newuserdata(L,sizeof(OperatorType));
	luaL_setmetatable(L,"Operator_Type");

	Create_Effective_Hamiltonian(G, F, energies, indices_active, list_of_kappas_active, list_of_kappas_virtual,occupations, NF, *U);
	return 1;
}

static int Test_Effective_H(lua_State *L)
{
	vector<Interpolation::InterpolatingFunction*> G;
		vector<Interpolation::InterpolatingFunction*> F;
		vector<int> list_of_kappas;
		lua_check(L,1,list_of_kappas);
		lua_check(L,2,G);
		lua_check(L,3,F);

		uint32_t N_funs = G.size();
		uint32_t N_x = (uint32_t) G[0]->x.size();

		vector<double> powx;
		Init_Slater_Kernel(list_of_kappas,G[0]->x,powx);

		double GFpairs[N_funs*N_funs*N_x];
		double *x = &(G[0]->x[0]);
		Init_Slater_Function_Pairs(N_x, G, F, GFpairs);

		double Q_av[N_funs*N_funs];
		for (int i = 0; i < N_funs; i++)
		{
			int k1 = list_of_kappas[i];
			int twoj1 = 2*std::abs(k1)-1;
			Q_av[i*N_funs + i] = Average_Coulomb_energy_two_shells(N_x, x, &(powx[0]), &(GFpairs[(i*N_funs+i)*N_x]), &(GFpairs[(i*N_funs+i)*N_x]), &(GFpairs[(i*N_funs+i)*N_x]), twoj1, twoj1, k1, k1);

			for (int j = i+1; j < N_funs; j++)
			{
				int k2 = list_of_kappas[j];
				int twoj2 = 2*std::abs(k2)-1;
				Q_av[i*N_funs + j] = Average_Coulomb_energy_two_shells(N_x, x, &(powx[0]), &(GFpairs[(i*N_funs+i)*N_x]), &(GFpairs[(j*N_funs+j)*N_x]), &(GFpairs[(i*N_funs+j)*N_x]), twoj1, twoj2, k1, k2);
				Q_av[j*N_funs + i] = Q_av[i*N_funs + j];
			}
		}
		vector<double> occupations = {2.,2.,2.,3.};
		int ind[4];
		int d_occ[4];
		int N_ex;

		int ind1[2] = {3,1};
		int d_occ1[2] = {1,-1};
		cout << Average_Coulomb_energy_difference_of_excitation(Q_av, N_funs, occupations, ind1,d_occ1, 2,list_of_kappas) << "\t";
		N_ex = Determine_excitations(0, 3, 1, 0, ind, d_occ);
		cout << Average_Coulomb_energy_difference_of_excitation(Q_av, N_funs, occupations, ind,d_occ, N_ex,list_of_kappas) << "\t";

	//	N_ex = Determine_excitations(3, 0, 1, 4, ind, d_occ);
	//	cout << Average_Coulomb_energy_difference_of_excitation(Q_av, N_funs, occupations, ind,d_occ, N_ex,list_of_kappas) << "\t";
	//	return 0;

		int ind2[2] = {4,1};
		int d_occ2[2] = {2,-2};
		cout << Average_Coulomb_energy_difference_of_excitation(Q_av, N_funs, occupations, ind2,d_occ2, 2,list_of_kappas) << "\t";
		N_ex = Determine_excitations(4, 4, 1, 1, ind, d_occ);
		cout << Average_Coulomb_energy_difference_of_excitation(Q_av, N_funs, occupations, ind,d_occ, N_ex,list_of_kappas) << "\t";
		int ind3[3] = {4,1,2};
		int d_occ3[3] = {2,-1,-1};
		cout << Average_Coulomb_energy_difference_of_excitation(Q_av, N_funs, occupations, ind3,d_occ3, 3,list_of_kappas) << "\t";
		N_ex = Determine_excitations(4, 4, 1, 2, ind, d_occ);
		cout << Average_Coulomb_energy_difference_of_excitation(Q_av, N_funs, occupations, ind,d_occ, N_ex,list_of_kappas) << "\t";
		int ind4[3] = {3,4,1};
		int d_occ4[3] = {1,1,-2};
		cout << Average_Coulomb_energy_difference_of_excitation(Q_av, N_funs, occupations, ind4,d_occ4, 3,list_of_kappas) << "\t";
		N_ex = Determine_excitations(3, 4, 1, 1, ind, d_occ);
		cout << Average_Coulomb_energy_difference_of_excitation(Q_av, N_funs, occupations, ind,d_occ, N_ex,list_of_kappas) << "\t";

		return 0;

		int i = 0;
		int k1 = list_of_kappas[i];
		int twoj1 = 1;
		double U0 = Average_Coulomb_energy_one_shell(N_x, x, &(powx[0]), &(GFpairs[(i*N_funs+i)*N_x]), twoj1) ;
		//cout << << "\n";
		int j = 1;
		int twoj2 = 1;
		int k2 = list_of_kappas[j];
		double U1 = Average_Coulomb_energy_two_shells(N_x, x, &(powx[0]), &(GFpairs[(i*N_funs+i)*N_x]), &(GFpairs[(j*N_funs+j)*N_x]), &(GFpairs[(i*N_funs+j)*N_x]), twoj1, twoj2, k1, k2) ;
		double U2 = Average_Coulomb_energy_one_shell(N_x, x, &(powx[0]), &(GFpairs[(j*N_funs+j)*N_x]), twoj2) ;
		//cout << << "\n";
		j = 3;
		twoj2 = 3;
		k2 = list_of_kappas[j];
		//cout << Average_Coulomb_energy_two_shells(N_x, x, &(powx[0]), &(GFpairs[(i*N_funs+i)*N_x]), &(GFpairs[(j*N_funs+j)*N_x]), &(GFpairs[(i*N_funs+j)*N_x]), twoj1, twoj2, k1, k2) << "\n";
		lua_push(L, U0);
		lua_push(L, U1);
		lua_push(L, U2);
		return 3;
}


