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




void Create_Effective_Hamiltonian(const vector<Interpolation::InterpolatingFunction*> &G,
										  const vector<Interpolation::InterpolatingFunction*> &F,
										  const vector<double> &energies,
										  const vector<vector<uint16_t> > &indices_active,
										  const vector<int> &list_of_kappas_active,
										  const vector<int> &list_of_kappas_virtual,
										  const uint32_t NF,
										  OperatorType &U)
{
	uint32_t N_active_orbitals  = (uint32_t) indices_active.size();
	uint32_t N_virtual_orbitals = (uint32_t) list_of_kappas_virtual.size();
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

//#pragma omp parallel for
	for (uint32_t i1 = 0; i1 < N_active_orbitals; i1++)
	{
		int k1 = list_of_kappas_active[i1];
		int twoj1 = 2*std::abs(k1)-1;
		double e1 = energies[i1];
		for (uint32_t i2 = 0; i2 < N_active_orbitals; i2++)
		{
			int k2 = list_of_kappas_active[i2];
			int twoj2 = 2*std::abs(k2)-1;
			double e2 = energies[i2];
			for (uint32_t i3 = i1; i3 < N_active_orbitals; i3++)
			{
				int k3 = list_of_kappas_active[i3];
				int twoj3 = 2*std::abs(k3)-1;
				double e3 = energies[i3];
				for (uint32_t i4 = i2; i4 < N_active_orbitals; i4++)
				{
					int k4 = list_of_kappas_active[i4];
					int twoj4 = 2*std::abs(k4)-1;
					double e4 = energies[i4];
					double pre_factor = 0;

					int K_start = std::max( abs(twoj1-twoj3), abs(twoj2-twoj4))/2;
					int K_end = std::min( twoj1+twoj3, twoj2+twoj4)/2;
					if ( K_start<= K_end)
					{
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

								double denom = e1v + e2v - e3 - e4;

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
									for (int L = L1_start; L <= L1_end; L +=2)
									{
										double c_reduced_11 = sqrt((twoj1+1)*(twoj1v+1))*ThreeJSymbolHalf(twoj1, -1, twoj1v, 1, 2*L, 0);
										//Cancelled if (((twoj1+1)/2)%2==1) c_reduced_11 *= -1; //TODO: move this to the outermoost possible loop
										double c_reduced_22 = sqrt((twoj2+1)*(twoj2v+1))*ThreeJSymbolHalf(twoj2, -1, twoj2v, 1, 2*L, 0);
										//Cancelled if (((twoj2+1)/2)%2==1) c_reduced_22 *= -1; //TODO: move this to the outermoost possible loop
										// This is working but slow: Fk1[(L - L1_start)/2] = Slater_Integrate(G[0]->x,powx,L,G[i1]->y,G[i2]->y,G[n1+N_active_orbitals]->y,G[n2+N_active_orbitals]->y,F[i1]->y,F[i2]->y,F[n1+N_active_orbitals]->y,F[n2+N_active_orbitals]->y);
										Fk1[(L - L1_start)/2] = Slater_Integrate_optimize(N_x, x, &(powx[0]), L, &(GFpairs[(i1*N_funs+(n1+N_active_orbitals))*N_x]), &(GFpairs[(i2*N_funs + n2+N_active_orbitals)*N_x]));
										Fk1[(L - L1_start)/2] *=  c_reduced_11 * c_reduced_22;

									}
									double Fk2[(L2_end-L2_start)/2 +1];
									//sign = L2_start % 2 == 0 ? 1 : -1; Cancelled  in the end
									for (int L = L2_start; L <= L2_end; L +=2)
									{
										double c_reduced_13 = sqrt((twoj3+1)*(twoj1v+1))*ThreeJSymbolHalf(twoj1v, -1, twoj3, 1, 2*L, 0);
										//Cancelled if (((twoj1v+1)/2)%2==1) c_reduced_13 *= -1; //TODO: move this to the outermoost possible loop
										double c_reduced_24 = sqrt((twoj2v+1)*(twoj4+1))*ThreeJSymbolHalf(twoj2v, -1, twoj4, 1, 2*L, 0);
										//Cancelled if (((twoj2v+1)/2)%2==1) c_reduced_24 *= -1; //TODO: move this to the outermoost possible loop
										// This is working but slow: Fk2[(L - L2_start)/2] = Slater_Integrate(G[0]->x,powx,L,G[n1+N_active_orbitals]->y,G[n2+N_active_orbitals]->y,G[i3]->y,G[i4]->y,F[n1+N_active_orbitals]->y,F[n2+N_active_orbitals]->y,F[i3]->y,F[i4]->y);
										Fk2[(L - L2_start)/2] = Slater_Integrate_optimize(N_x, x, &(powx[0]), L, &(GFpairs[((n1+N_active_orbitals)*N_funs+i3)*N_x]), &(GFpairs[((n2+N_active_orbitals)*N_funs + i4)*N_x]));
										Fk2[(L - L2_start)/2] *= c_reduced_13 * c_reduced_24;
									}
									for (int K = K_start; K <= K_end; K++)
									{
										for (int L2 = L2_start; L2 <= L2_end; L2 +=2)
										{
											for (int L1 = L1_start; L1 <= L1_end; L1 +=2)
											{
												double angular_prefactor = (2*K+1)*Wigner9J(twoj2, twoj2v, 2*L1, twoj3, 2*L2, twoj1v, 2*K, twoj4, twoj1)/4.;//  1/4 is the prefactor according to Effective_Hamiltonian_2.pdf
												reduced_matrix_ements[K-K_start] += angular_prefactor * Fk1[(L1 - L1_start)/2] * Fk2[(L2 - L2_start)/2]/denom;
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

								double denom = e1v + e2v - e3 - e4;

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
									for (int L = L1_start; L <= L1_end; L +=2)
									{
										double c_reduced_11 = sqrt((twoj1+1)*(twoj1v+1))*ThreeJSymbolHalf(twoj1, -1, twoj1v, 1, 2*L, 0);
										//Cancelled if (((twoj1+1)/2)%2==1) c_reduced_11 *= -1; //TODO: move this to the outermoost possible loop
										double c_reduced_22 = sqrt((twoj2+1)*(twoj2v+1))*ThreeJSymbolHalf(twoj2, -1, twoj2v, 1, 2*L, 0);
										//Cancelled if (((twoj2+1)/2)%2==1) c_reduced_22 *= -1; //TODO: move this to the outermoost possible loop
										// This is working but slow: Fk1[(L - L1_start)/2] = Slater_Integrate(G[0]->x,powx,L,G[i1]->y,G[i2]->y,G[n1+N_active_orbitals]->y,G[n2]->y,F[i1]->y,F[i2]->y,F[n1+N_active_orbitals]->y,F[n2]->y);
										Fk1[(L - L1_start)/2] = Slater_Integrate_optimize(N_x, x, &(powx[0]), L, &(GFpairs[(i1*N_funs+(n1+N_active_orbitals))*N_x]), &(GFpairs[(i2*N_funs + n2)*N_x]));
										Fk1[(L - L1_start)/2] *=  c_reduced_11 * c_reduced_22;

									}
									double Fk2[(L2_end-L2_start)/2 +1];
									//sign = L2_start % 2 == 0 ? 1 : -1; Cancelled  in the end
									for (int L = L2_start; L <= L2_end; L +=2)
									{
										double c_reduced_13 = sqrt((twoj3+1)*(twoj1v+1))*ThreeJSymbolHalf(twoj1v, -1, twoj3, 1, 2*L, 0);
										//Cancelled if (((twoj1v+1)/2)%2==1) c_reduced_13 *= -1; //TODO: move this to the outermoost possible loop
										double c_reduced_24 = sqrt((twoj2v+1)*(twoj4+1))*ThreeJSymbolHalf(twoj2v, -1, twoj4, 1, 2*L, 0);
										//Cancelled if (((twoj2v+1)/2)%2==1) c_reduced_24 *= -1; //TODO: move this to the outermoost possible loop
										// This is working but slow: Fk2[(L - L2_start)/2] = Slater_Integrate(G[0]->x,powx,L,G[n1+N_active_orbitals]->y,G[n2]->y,G[i3]->y,G[i4]->y,F[n1+N_active_orbitals]->y,F[n2]->y,F[i3]->y,F[i4]->y);
										Fk2[(L - L2_start)/2] = Slater_Integrate_optimize(N_x, x, &(powx[0]), L, &(GFpairs[((n1+N_active_orbitals)*N_funs+i3)*N_x]), &(GFpairs[((n2)*N_funs + i4)*N_x]));
										Fk2[(L - L2_start)/2] *= c_reduced_13 * c_reduced_24;
									}
									for (int K = K_start; K <= K_end; K++)
									{
										for (int L2 = L2_start; L2 <= L2_end; L2 +=2)
										{
											for (int L1 = L1_start; L1 <= L1_end; L1 +=2)
											{
												double angular_prefactor = (2*K+1)*Wigner9J(twoj2, twoj2v, 2*L1, twoj3, 2*L2, twoj1v, 2*K, twoj4, twoj1)/2.; //  1/2 is the prefactor according to Effective_Hamiltonian_2.pdf
												reduced_matrix_ements[K-K_start] += angular_prefactor * Fk1[(L1 - L1_start)/2] * Fk2[(L2 - L2_start)/2]/denom;
											}
										}
									}
								}
							}
						}
						cout << i1 << "\t"<< i2 << "\t"<< i3 << "\t"<< i4 << "\n";
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
										if (((two_m2-two_m4)/2)%2==1) angular_prefactor *= -1.;
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
				}
			}
		}
	}
}

static int LuaCreateEffectiveHamiltonian(lua_State *L)
{
	vector<Interpolation::InterpolatingFunction*> G;
	vector<Interpolation::InterpolatingFunction*> F;
	vector<double> energies;
	vector<vector<uint16_t> > indices_active;
	vector<int> list_of_kappas_active;
	vector<int> list_of_kappas_virtual;
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

	OperatorType *U = (OperatorType *) lua_newuserdata(L,sizeof(OperatorType));
	luaL_setmetatable(L,"Operator_Type");

	Create_Effective_Hamiltonian(G, F, energies, indices_active, list_of_kappas_active, list_of_kappas_virtual, NF, *U);
	return 1;
}



