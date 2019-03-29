 
double vecNorm(const vector<double> &v)
{
	double normSqrd = 0;
	for (uint i = 0;  i < (uint)v.size(); ++i)
		normSqrd += v[i]*v[i];
	return sqrt(normSqrd);
}

double multicenterMultipoleIntegral(Interpolation::InterpolatingFunction &f1, Interpolation::InterpolatingFunction &f2, int k, uint32_t NknotsPerInterval)
{
//	TODO
	return 0.;
}

vector<double> interpolationRanges(const vector<Interpolation::InterpolatingFunction*> &funs)
{
	vector<double> out(funs.size());
	for (uint i = 0; i < (uint)funs.size(); ++i)
		out[i] = funs[i]->GetX()[funs[i]->GetX().size() - 1];
	return out;
}

vector<Complex16> realspaceBasisValues(const vector<double> &v, const vector<string> &orbs, const vector<vector<double> > &atomPos, const vector<double> &interRanges,
		const vector<Interpolation::InterpolatingFunction*> &funs, const vector<vector<uint16_t> > indices, const double step)
{
	vector < Complex16 > psi;
	for (uint orb = 0; orb < (uint) orbs.size(); ++orb)
	{
		vector<double> rVec(3);
		for (uint k = 0; k < 3; ++k)
			rVec[k] = v[k] - atomPos[orb][k];
		double r = vecNorm(rVec);
		Complex16 val;
		if (r <= interRanges[orb])
		{
			if (r > step / 2)
			{
				double theta = atan2(rVec[0], rVec[1]);
				double phi = acos(rVec[2] / r);
				double psiR = (*funs[orb])(r) / r;
				int l = ((int) indices[orb].size() / 2 - 1) / 2;
				for (int m = -l - 1; m <= l + 1; ++m)
				{
					val.Re = psiR * ReSphericalHarmonicY(l, m, theta, phi);
					val.Im = psiR * ImSphericalHarmonicY(l, m, theta, phi);
					psi.push_back(val);
				}
			} else
			{
				double multicellFactor = 1.;
				if (r == step / 2)
				{
					multicellFactor = 2.;
				}
				val.Re = (*funs[orb])(r) / step / multicellFactor;
				val.Im = 0;
				for (uint m = 0; m < (uint) indices[orb].size(); ++m)
					psi.push_back(val);
			}
		} else
		{
			val.Re = 0;
			val.Im = 0;
			for (uint m = 0; m < (uint) indices[orb].size(); ++m)
				psi.push_back(val);
		}
	}
	return psi;
}

CompactMatrixType densityMatrixNoSpin(const CompactMatrixType &rho)
{
	CompactMatrixType rhoNoSpin;
	rhoNoSpin.Ni = rho.Ni / 2;
	rhoNoSpin.Nj = rho.Nj / 2;
	rhoNoSpin.QComplex = (uint) true;
	InitCompactMatrix(&rhoNoSpin);
	for (uint i = 0; i < rhoNoSpin.Ni; ++i)
		for (uint j = 0; j < rhoNoSpin.Nj; ++j)
			if (rho.QComplex)
			{
				rhoNoSpin.V[i * rhoNoSpin.Nj + j].Re = rho.V[2 * i * rho.Nj + 2 * j].Re + rho.V[(2 * i + 1) * rho.Nj + (2 * j + 1)].Re;
				rhoNoSpin.V[i * rhoNoSpin.Nj + j].Im = rho.V[2 * i * rho.Nj + 2 * j].Im + rho.V[(2 * i + 1) * rho.Nj + (2 * j + 1)].Im;
			} else
			{
				rhoNoSpin.V[i * rhoNoSpin.Nj + j].Re = rho.VRe[2 * i * rho.Nj + 2 * j] + rho.VRe[(2 * i + 1) * rho.Nj + (2 * j + 1)];
				rhoNoSpin.V[i * rhoNoSpin.Nj + j].Im = 0;
			}
	return rhoNoSpin;
}

map<vector<double>, double> realspaceWavefunction(const CompactMatrixType &rho, const vector<Interpolation::InterpolatingFunction*> &funs, const vector<string> &orbs,
		const vector<vector<double> > &atomPos, const uint gridpoints, const double domainsize)
{
	CompactMatrixType rhoNoSpin = densityMatrixNoSpin(rho);

	const vector<vector<uint16_t> > indices = dictToList(createAtmoicIndicesDict(orbs), orbs);
	const vector<double> interRanges = interpolationRanges(funs);
	double step = domainsize / (double) (gridpoints - 1);
	map<vector<double>, double> result;
#pragma omp parallel for
	for (uint xi = 0; xi <= gridpoints; ++xi)
	{
		double x = -domainsize / 2 + step * xi;
		for (double y = -domainsize / 2; y <= domainsize / 2; y += step)
		{
			for (double z = -domainsize / 2; z <= domainsize / 2; z += step)
			{
				vector<double> v(3);
				v[0] = x;
				v[1] = y;
				v[2] = z;

				if (vecNorm(v) <= domainsize / 2)
				{
					const vector<Complex16> psi = realspaceBasisValues(v, orbs, atomPos, interRanges, funs, indices, step);
					double WFval = 0;
					for (uint i = 0; i < (uint) rhoNoSpin.Ni; ++i)
					{
						for (uint j = 0; j < (uint) rhoNoSpin.Nj; ++j)
						{
//						Any off-diagonal imaginary parts cancel out because rho is hermetian.
							if (rhoNoSpin.V[i * rhoNoSpin.Nj + j].Re != 0)
							{
								WFval += rhoNoSpin.V[i * rhoNoSpin.Nj + j].Re * (psi[i].Re * psi[j].Re + psi[i].Im * psi[j].Im);
							}
						}
					}
					result.insert( { v, WFval });
				}
			}
		}
	}
	FreeCompactMatrix(&rhoNoSpin);
	return result;
}

map<vector<double>, double> realSpaceWavefunction(const CompactMatrixType &rho, const vector<Interpolation::InterpolatingFunction*> &funs, const vector<string> &orbs,
		const vector<vector<double> > &atomPos, const uint gridpoints)
{
	double maxNorm = 0;
	double maxX = 0;
	for (uint i = 0; i < (uint) atomPos.size(); ++i)
	{
		if (vecNorm(atomPos[i]) > maxNorm)
			maxNorm = vecNorm(atomPos[i]);
		double xMax = funs[i]->GetX()[funs[i]->GetX().size() - 1];
		if (xMax > maxX)
			maxX = xMax;
	}
	double domainsize = 2 * (maxNorm + maxX);
	return realspaceWavefunction(rho, funs, orbs, atomPos, gridpoints, domainsize);
}

template<typename keyType, typename valType>
void lua_pushMapKeyValPair(lua_State * L, const map<keyType, valType> &pushkin)
{
	lua_newtable(L);
	int top = (int) lua_gettop(L);
	typename map<keyType, valType>::const_iterator iter;
	int i = 1;
	for (iter = pushkin.begin(); iter != pushkin.end(); ++iter)
	{
		lua_push(L, i++);

		lua_newtable(L);
		lua_push(L,1);
		lua_push(L, iter->first);
		lua_settable(L, -3);
		lua_push(L,2);
		lua_push(L, iter->second);
		lua_settable(L, -3);

		lua_settable(L, top);
	}
}

static int LuaRealspaceWavefunction(lua_State * L)
{
	CompactMatrixType rho;
	vector<Interpolation::InterpolatingFunction*> funs;
	vector < string > orbs;
	vector < vector<double> > atomPos;
	bool hasAtomPos = false;
	struct
	{
		uint gridpoints;
		double domain;
	} Options;
	LuaOptionsType LuaRealspaceWavefunctionOptions[] =
	{ { "gridpoints", "int", &(Options.gridpoints) },
			{ "domain", "double", &(Options.domain) },
			{ NULL, NULL, NULL } };
	Options.gridpoints = 200;
	double domain = 0;
	bool hasOptions = false;

	uint len = lua_gettop(L);
	if (len < 3 || len > 5)
		return luaL_error(L, "Error in RealspaceWavefunction: Expected 3, 4 or 5 arguments, density matrix, radial functions, orbitals, optional atom positions, and options.");

	lua_check(L, 1, rho);
	lua_check(L, 2, funs);
	lua_check(L, 3, orbs);
	if (len == 4)
	{
		if (luaL_isType(L, 4, atomPos))
		{
			hasAtomPos = true;
		} else
		{
			hasOptions = true;
		}
	} else if (len == 5)
	{
		hasAtomPos = true;
		hasOptions = true;
	}

	if (hasAtomPos)
	{
		lua_check(L, 4, atomPos);
	} else
	{
		vector<double> zeroVec(3);
		zeroVec[0] = 0.;
		zeroVec[1] = 0.;
		zeroVec[2] = 0.;
		for (uint i = 0; i < (uint) orbs.size(); ++i)
			atomPos.push_back(zeroVec);
	}

	if (hasOptions)
	{
		if (LuaTestOptionsList(L, -1))
		{
			LuaCheckOptions(L, -1, LuaRealspaceWavefunctionOptions);
		} else
		{
			return luaL_error(L, "Kind of expected a list of options as last argument of RealspaceWavefunction.\n");
		}
	}

	if (rho.Ni != rho.Nj)
		return luaL_error(L, "Error in RealspaceWavefunction: Rho needs to be a quadratic matrix, but is of size %dx%d.", rho.Ni, rho.Nj);

	if (funs.size() != orbs.size() || orbs.size() != atomPos.size())
		return luaL_error(L, "Error in RealspaceWavefunction: List of functions (%d), list of orbitals (%d) and list of atom positions (%d) need to be of the same size.", (uint) funs.size(),
				(uint) orbs.size(), (uint) atomPos.size());

	uint NF = NumberOfFermionicStates(orbs);
	if (NF != rho.Ni)
		return luaL_error(L, "Error in RealspaceWavefunction: Number of fermionic states obtained from orbitals (%d) is unequal to size of density matrix (%d).", NF, rho.Ni);

//	TODO: Warning if rho is not hermetian

	if (Options.domain == 0)
		lua_pushMapKeyValPair(L, realSpaceWavefunction(rho, funs, orbs, atomPos, Options.gridpoints));
	else
		lua_pushMapKeyValPair(L, realspaceWavefunction(rho, funs, orbs, atomPos, Options.gridpoints, Options.domain));
	FreeCompactMatrix(&rho);
	return 1;
}
