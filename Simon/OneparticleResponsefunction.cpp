/* The goal is to calculate a response function
 * G(w)=<psi|T^dag 1/(w-H-Offset+i*Gamma) T|psi>
 * where all operators are one-particle and |psi> is assumed to be |111...000>.
 * Finding R such that R H R^dag = diag(E_j) gives
 * G(w)=<psi|T^dag R^dag 1/diag(w-E_j-Offset+i*Gamma) R T|psi>
 * To have as little matrix multiplications as possible we can rewrite this expectation value as
 * G(w)=Tr(rho T^dag R^dag 1/diag(w-E_j-Offset+i*Gamma) R T)
 * and use the cyclic property
 * G(w)=Tr(R T rho T^dag R^dag 1/diag(w-E_j-Offset+i*Gamma) )
 * Defining
 * M = R T rho T^dag R^dag
 * and
 * D = 1/diag(w-E_j-Offset+i*Gamma)
 * means that
 * G(w)=Sum_{i,j} M_ij D_ji = Sum_i M_ii D_i
 * the calculation of which can be parallelised in w.
 */
vector<Complex16> OneparticleResponsefunctionDiagonal(const OperatorType &H, const OperatorType &T, const uint Nelec, const double wmin, const double wmax, const double dw, const double Offset,
		const double Gamma)
{
	CompactMatrixType Zero, HMat, TMat, R, rhogrd, RT, matprod;
	double *Hevals;
	double E0 = H.l[0].N == 0 ? 0 : H.l[0].VRe[0];
	Zero.Ni = H.NFermion;
	Zero.Nj = H.NFermion;
	Zero.QComplex = 0;
	InitCompactMatrix(&Zero);
	MeanFieldOperatorOneElectronMatrix(&H, &HMat, Zero);
	MeanFieldOperatorOneElectronMatrix(&T, &TMat, Zero);
	Hevals = (double*) malloc(HMat.Ni * sizeof(double));
	CompactMatrixEigensystem(&HMat, &R, Hevals);
	rhogrd.Ni = HMat.Ni;
	rhogrd.Nj = HMat.Nj;
	rhogrd.QComplex = 0;
	InitCompactMatrix(&rhogrd);
	for (uint i = 0; i < Nelec; ++i)
		rhogrd.VRe[i * rhogrd.Nj + i] = 1;
	CompactMatrixMultiply('N', 'N', &R, &TMat, &RT);
	CompactMatrixMultiply3('N', 'N', 'C', &RT, &rhogrd, &RT, &matprod);
	uint iwmax = (uint)((wmax - wmin) / dw) + 1;
	if (wmin > wmax) return vector < Complex16 > (0);
	vector < Complex16 > out(iwmax);
#pragma omp parallel for
	for (uint iw = 0; iw < iwmax; ++iw)
	{
		double w = wmin + (iw) * dw;
		Complex16 val;
		val.Re = 0;
		val.Im = 0;
		for (uint j = 0; j < HMat.Ni; ++j)
		{
			Complex16 a = matprod.V[j * matprod.Nj + j];
			Complex16 b;
			b.Re = w - Offset - Hevals[j] - E0;
			b.Im = Gamma;
			val.Re = val.Re + (a.Re * b.Re + a.Im * b.Im) / (Sqr(b.Im) + Sqr(b.Re));
			val.Im = val.Im + (a.Im * b.Re - a.Re * b.Im) / (Sqr(b.Im) + Sqr(b.Re));
		}
		out[iw] = val;
	}
	FreeCompactMatrix(&Zero);
	FreeCompactMatrix(&HMat);
	FreeCompactMatrix(&TMat);
	FreeCompactMatrix(&rhogrd);
	FreeCompactMatrix(&RT);
	FreeCompactMatrix(&matprod);
	free(Hevals);
	return out;
}

static int LuaOneparticleResponsefunctionDiagonal(lua_State *L)
{
	uint32_t len = (uint32_t) lua_gettop(L);
	if (len != 8)
		return luaL_error(L,
				"Error in OneparticleResponsefunctionDiagonal: Eight elements expected, Hamiltonian, Transition operator, number of electrons, energy minimum, energy maximum, energy step, Offset and Broadening.");
	OperatorType H, T;
	uint Nelec;
	double wmin, wmax, dw, Offset, Gamma;
	if (!luaL_isType(L, 1, H)) return luaL_error(L, "Error in OneparticleResponsefunctionDiagonal: 1. argument should be the Hamiltonian (i.e. an operator)");
	if (!luaL_isType(L, 2, T)) return luaL_error(L, "Error in OneparticleResponsefunctionDiagonal: 2. argument should be the Transition operator (i.e. an operator)");
	if (!luaL_isType(L, 3, Nelec)) return luaL_error(L, "Error in OneparticleResponsefunctionDiagonal: 3. argument should be the number of electrons");
	if (!luaL_isType(L, 4, wmin)) return luaL_error(L, "Error in OneparticleResponsefunctionDiagonal: 4. argument should be the energy minimum");
	if (!luaL_isType(L, 5, wmax)) return luaL_error(L, "Error in OneparticleResponsefunctionDiagonal: 5. argument should be the energy maximum");
	if (!luaL_isType(L, 6, dw)) return luaL_error(L, "Error in OneparticleResponsefunctionDiagonal: 6. argument should be the energy step");
	if (!luaL_isType(L, 7, Offset)) return luaL_error(L, "Error in OneparticleResponsefunctionDiagonal: 7. argument should be the energy offset");
	if (!luaL_isType(L, 8, Gamma)) return luaL_error(L, "Error in OneparticleResponsefunctionDiagonal: 8. argument should be the broadening");
	lua_check(L, 1, H);
	lua_check(L, 2, T);
	lua_check(L, 3, Nelec);
	lua_check(L, 4, wmin);
	lua_check(L, 5, wmax);
	lua_check(L, 6, dw);
	lua_check(L, 7, Offset);
	lua_check(L, 8, Gamma);
	lua_push(L, OneparticleResponsefunctionDiagonal(H, T, Nelec, wmin, wmax, dw, Offset, Gamma));
	return 1;
}
