void lua_push(lua_State * L, const doubleVecValPair &pushkin)
{
	lua_newtable(L);
	lua_push(L, 1);
	lua_push(L, pushkin.vec);
	lua_settable(L, -3);
	lua_push(L, 2);
	lua_push(L, pushkin.val);
	lua_settable(L, -3);
}

void lua_check(lua_State * L, int arg, doubleVecValPair &in)
{
	lua_rawgeti(L, arg, 1);
	lua_check(L, -1, in.vec);
	lua_pop(L, 1);
	lua_rawgeti(L, arg, 2);
	lua_check(L, -1, in.val);
	lua_pop(L, 1);
}

double vecNorm(const vector<double> &v)
{
	double normSqrd = 0;
	for (uint i = 0;  i < (uint)v.size(); ++i)
		normSqrd += v[i]*v[i];
	return sqrt(normSqrd);
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

vector<doubleVecValPair> realspaceElectronDensity(const CompactMatrixType &rho, const vector<Interpolation::InterpolatingFunction*> &funs, const vector<string> &orbs, const vector<vector<double> > &atomPos,
		const uint gridpoints, const double domainsize)
{
	CompactMatrixType rhoNoSpin = densityMatrixNoSpin(rho);

	const vector<vector<uint16_t> > indices = dictToList(createAtmoicIndicesDict(orbs), orbs);
	const vector<double> interRanges = interpolationRanges(funs);
	double step = domainsize / (double) (gridpoints - 1);
	vector<doubleVecValPair> result(gridpoints * gridpoints * gridpoints);
#pragma omp parallel for
	for (uint xi = 0; xi < gridpoints; ++xi)
	{
		double x = -domainsize / 2 + step * xi;
		for (uint yi = 0; yi < gridpoints; ++yi)
		{
			double y = -domainsize / 2 + step * yi;
			for (uint zi = 0; zi < gridpoints; ++zi)
			{
				double z = -domainsize / 2 + step * zi;
				uint vecI = xi * gridpoints * gridpoints + yi * gridpoints + zi;
				vector<double> v(3);
				v[0] = x;
				v[1] = y;
				v[2] = z;

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
				doubleVecValPair tmp;
				tmp.vec = v;
				tmp.val = WFval;
				result[vecI] = tmp;
			}
		}
	}
	FreeCompactMatrix(&rhoNoSpin);
	return result;
}

vector<doubleVecValPair> realspaceElectronDensity(const CompactMatrixType &rho, const vector<Interpolation::InterpolatingFunction*> &funs, const vector<string> &orbs,
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
	return realspaceElectronDensity(rho, funs, orbs, atomPos, gridpoints, domainsize);
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

vector<double> vecAverage(const vector<double> &v1, const vector<double> &v2)
{
	vector<double> result(v1.size());
	for (uint i = 0; i < v1.size(); ++i)
		result[i] = (v1[i] + v2[i]) / 2.;
	return result;
}

vector<double> vecAverage(const vector<double> &v1, const vector<double> &v2, double val1, double val2, double valMiddle)
{
	vector<double> result(v1.size());
	for (uint i = 0; i < v1.size(); ++i)
		if (val1 != val2)
			result[i] = v1[i] * (val2 - valMiddle) / (val2 - val1) + v2[i] * (valMiddle - val1) / (val2 - val1);
		else
			result[i] = (v1[i] + v2[i]) / 2.;
	return result;
}

double findIsosurfaceVal(const vector<doubleVecValPair> &density, double portion, double step)
{
	double fullSum = 0;
	for (uint i = 0; i < (uint) density.size(); ++i)
		fullSum += density[i].val;
//	TODO: Verbosity Print for fullInt.
	double threshold = 0;
	double partialSum = 0;
	while (partialSum < fullSum * portion)
	{
		partialSum = 0;
		threshold += step;
		for (uint i = 0; i < (uint) density.size(); ++i)
			if (density[i].val < threshold)
				partialSum += density[i].val;
	}
	return threshold;
}

double findIsosurfaceVal(const vector<doubleVecValPair> &density, double portion)
{
	return findIsosurfaceVal(density, portion, 0.01);
}

vector<vector<double> > isosurface(const vector<doubleVecValPair> &density, double threshold)
{
	uint gridpoints = (uint) pow((double) density.size() + 1, 1 / 3.);
	vector < vector<double> > result;
	for (uint xi = 0; xi < gridpoints - 1; ++xi)
	{
		for (uint yi = 0; yi < gridpoints - 1; ++yi)
		{
			for (uint zi = 0; zi < gridpoints - 1; ++zi)
			{
				uint i0 = xi * gridpoints * gridpoints + yi * gridpoints + zi;
				uint ixpp = (xi + 1) * gridpoints * gridpoints + yi * gridpoints + zi;
				uint iypp = xi * gridpoints * gridpoints + (yi + 1) * gridpoints + zi;
				uint izpp = xi * gridpoints * gridpoints + yi * gridpoints + (zi + 1);
				double val0 = density[i0].val;
				double valxpp = density[ixpp].val;
				double valypp = density[iypp].val;
				double valzpp = density[izpp].val;
				if ((val0 < threshold && valxpp > threshold) || (val0 > threshold && valxpp < threshold))
					result.push_back(vecAverage(density[i0].vec, density[ixpp].vec, density[i0].val, density[ixpp].val, threshold));
				if ((val0 < threshold && valypp > threshold) || (val0 > threshold && valypp < threshold))
					result.push_back(vecAverage(density[i0].vec, density[iypp].vec, density[i0].val, density[iypp].val, threshold));
				if ((val0 < threshold && valzpp > threshold) || (val0 > threshold && valzpp < threshold))
					result.push_back(vecAverage(density[i0].vec, density[izpp].vec, density[i0].val, density[izpp].val, threshold));
			}
		}
	}
	return result;
}

static int LuaRealspaceDensityCheckAndSecurityMeasures(lua_State * L, CompactMatrixType &rho, vector<Interpolation::InterpolatingFunction*> &funs, vector<string> &orbs,
		vector<vector<double> > &atomPos, bool hasAtomPos)
{
	lua_check(L, 1, rho);
	lua_check(L, 2, funs);
	lua_check(L, 3, orbs);
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

	if (rho.Ni != rho.Nj)
		return luaL_error(L, "Error in calculating RealspaceDensity: Rho needs to be a quadratic matrix, but is of size %dx%d.", rho.Ni, rho.Nj);

	if (funs.size() != orbs.size() || orbs.size() != atomPos.size())
		return luaL_error(L, "Error in calculating RealspaceDensity: List of functions (%d), list of orbitals (%d) and list of atom positions (%d) need to be of the same size.", (uint) funs.size(),
				(uint) orbs.size(), (uint) atomPos.size());

	uint NF = NumberOfFermionicStates(orbs);
	if (NF != rho.Ni)
		return luaL_error(L, "Error in calculating RealspaceDensity: Number of fermionic states obtained from orbitals (%d) is unequal to size of density matrix (%d).", NF, rho.Ni);

	if (!QCompactMatrixHermitian(rho, DBL_EPSILON * 10))
		cout << "Warning: density matrix for calculating RealspaceDensity is not fully hermitian. The code will assume it to be." << endl;
	return 0;
}

static int LuaRealspaceDensity(lua_State * L)
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
	LuaOptionsType LuaRealspaceDensityOptions[] =
	{ { "gridpoints", "int", &(Options.gridpoints) },
			{ "domain", "double", &(Options.domain) },
			{ NULL, NULL, NULL } };
	Options.gridpoints = 200;
	Options.domain = 0;
	bool hasOptions = false;

	uint len = lua_gettop(L);
	if (len < 3 || len > 5)
		return luaL_error(L, "Error in RealspaceDensity: Expected 3, 4 or 5 arguments, density matrix, radial functions, orbitals, optional atom positions, and options.");

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

	if (hasOptions)
	{
		if (LuaTestOptionsList(L, -1))
		{
			LuaCheckOptions(L, -1, LuaRealspaceDensityOptions);
		} else
		{
			return luaL_error(L, "Kind of expected a list of options as last argument of RealspaceDensity.\n");
		}
	}

	LuaRealspaceDensityCheckAndSecurityMeasures(L, rho, funs, orbs, atomPos, hasAtomPos);

	if (Options.domain == 0)
		lua_push(L, realspaceElectronDensity(rho, funs, orbs, atomPos, Options.gridpoints));
	else
		lua_push(L, realspaceElectronDensity(rho, funs, orbs, atomPos, Options.gridpoints, Options.domain));
	return 1;
}

static int LuaIsosurface(lua_State * L, bool hasThreshold)
{
	struct
	{
		uint gridpoints;
		double domain;
		double stepsize;
	} Options;
	LuaOptionsType LuaIsosurfaceOptions[] = {
			{ "gridpoints", "int", &(Options.gridpoints) },
			{ "domain", "double", &(Options.domain) },
			{ "stepsize", "double", &(Options.stepsize) },
			{ NULL, NULL, NULL } };
	Options.gridpoints = 200;
	Options.domain = 0;
	Options.stepsize = 0.01;
	bool hasRealspaceDensity;
	bool hasAtomPos;
	bool hasOptions;
	int numberInputPos;

	uint len = lua_gettop(L);
	switch (len)
	{
	case 2:
		hasRealspaceDensity = true;
		hasAtomPos = false;
		hasOptions = false;
		numberInputPos = 2;
		break;
	case 3:
		hasRealspaceDensity = true;
		hasAtomPos = false;
		hasOptions = true;
		numberInputPos = 2;
		break;
	case 4:
		hasRealspaceDensity = false;
		hasAtomPos = false;
		hasOptions = false;
		numberInputPos = 4;
		break;
	case 5:
		hasRealspaceDensity = false;
		if ((bool) lua_isnumber(L, 5))
		{
			hasAtomPos = true;
			hasOptions = false;
			numberInputPos = 5;
		} else
		{
			hasAtomPos = false;
			hasOptions = true;
			numberInputPos = 4;
		}
		break;
	case 6:
		hasRealspaceDensity = false;
		hasAtomPos = true;
		hasOptions = true;
		numberInputPos = 5;
		break;
	default:
		return luaL_error(L,
				"Error in Isosurface: Expected something between 2 and 6 arguments: \n2: Realspace density and isosurface value or integral encapsulated by isosurface. \n3: Realspace density, value or integral, and options. \n4: Density matrix, radial funtions, orbitals, and value or integral. \n5: Density matrix, radial functions, orbitals, and either atom positions and value or integral, or value or integral and options. \n6: Density matrix, radial functions, orbitals, atom positions, value or integral, and options.");
	}

	if (hasOptions)
	{
		if (LuaTestOptionsList(L, -1))
		{
			LuaCheckOptions(L, -1, LuaIsosurfaceOptions);
		} else
		{
			return luaL_error(L, "Kind of expected a list of options as last argument of Isosurface.\n");
		}
	}

	vector<doubleVecValPair> density;
	if (hasRealspaceDensity)
	{
		lua_check(L, 1, density);
	} else
	{
		CompactMatrixType rho;
		vector<Interpolation::InterpolatingFunction*> funs;
		vector < string > orbs;
		vector < vector<double> > atomPos;

		LuaRealspaceDensityCheckAndSecurityMeasures(L, rho, funs, orbs, atomPos, hasAtomPos);

		if (Options.domain == 0)
			density = realspaceElectronDensity(rho, funs, orbs, atomPos, Options.gridpoints);
		else
			density = realspaceElectronDensity(rho, funs, orbs, atomPos, Options.gridpoints, Options.domain);
	}

	if (hasThreshold)
	{
		double threshold;
		lua_check(L, numberInputPos, threshold);
		lua_push(L, isosurface(density, threshold));
		return 1;
	} else
	{
		double portion;
		lua_check(L, numberInputPos, portion);
		double threshold = findIsosurfaceVal(density, portion, Options.stepsize);
		lua_push(L, isosurface(density, threshold));
		lua_push(L, threshold);
		return 2;
	}
}

static int LuaIsosurfaceFixedVal(lua_State * L)
{
	return LuaIsosurface(L, true);
}

static int LuaIsosurfaceFixedInt(lua_State * L)
{
	return LuaIsosurface(L, false);
}
