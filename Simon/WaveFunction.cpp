
WaveFunctionType WaveFunctionTensorProduct(const WaveFunctionType &psi1, const WaveFunctionType &psi2)
{
	WaveFunctionType out;
	out.NFermion = psi1.NFermion + psi2.NFermion;
	out.NBoson = 0;
	if ((bool) psi1.QComplex || (bool) psi2.QComplex)
	{
		out.QComplex = (uint32_t) true;
	}
	else
	{
		out.QComplex = (uint32_t) false;
	}
	out.NBitsKey = std::max(psi1.NBitsKey, psi2.NBitsKey);
	WaveFunctionInit(&out);
	if ((bool) psi1.QComplex && (bool) psi2.QComplex)
	{
		for (uint32_t i = 1; i <= psi1.N; ++i)
		{
			for (uint32_t j = 1; j <= psi2.N; ++j)
			{
				double valRe = psi1.ValueRe[BlokElement(i)] * psi2.ValueRe[BlokElement(j)] - psi1.ValueIm[BlokElement(i)] * psi2.ValueIm[BlokElement(j)];
				double valIm = psi1.ValueRe[BlokElement(i)] * psi2.ValueIm[BlokElement(j)] + psi1.ValueIm[BlokElement(i)] * psi2.ValueRe[BlokElement(j)];
				std::stringstream ss;
				char *det1, *det2;
				uint8_t *detNew;
				det1 = DetToString(psi1.NFermion, psi1.NBoson, &psi1.Determinant[BlokElementDet(i, psi1.NWord, 0)]);
				det2 = DetToString(psi2.NFermion, psi2.NBoson, &psi2.Determinant[BlokElementDet(j, psi2.NWord, 0)]);
				ss << det1 << det2;
				detNew = StringToDet(out.NFermion, out.NBoson, ss.str().c_str());
				ComplexWaveFunctionAddElement(&out, detNew, valRe, valIm);
				free(det1);
				free(det2);
				free(detNew);
			}
		}
	}
	else if ((bool) psi1.QComplex || (bool) psi2.QComplex) //one complex, one real
	{
		WaveFunctionType tmp1, tmp2;
		if ((bool) psi2.QComplex)
		{
			tmp1 = psi2;
			tmp2 = psi1;
		}
		else
		{
			tmp1 = psi1;
			tmp2 = psi2;
		}
		for (uint32_t i = 1; i <= tmp1.N; ++i)
		{
			for (uint32_t j = 1; j <= tmp2.N; ++j)
			{
				double valRe = tmp1.ValueRe[BlokElement(i)] * tmp2.ValueRe[BlokElement(j)];
				double valIm = tmp1.ValueIm[BlokElement(i)] * tmp2.ValueRe[BlokElement(j)];
				std::stringstream ss;
				char *det1, *det2;
				uint8_t *detNew;
				det1 = DetToString(tmp1.NFermion, tmp1.NBoson, &tmp1.Determinant[BlokElementDet(i, tmp1.NWord, 0)]);
				det2 = DetToString(tmp2.NFermion, tmp2.NBoson, &tmp2.Determinant[BlokElementDet(j, tmp2.NWord, 0)]);
				ss << det1 << det2;
				detNew = StringToDet(out.NFermion, out.NBoson, ss.str().c_str());
				ComplexWaveFunctionAddElement(&out, detNew, valRe, valIm);
				free(det1);
				free(det2);
				free(detNew);
			}
		}
	}
	else //all real
	{
		for (uint32_t i = 1; i <= psi1.N; ++i)
		{
			for (uint32_t j = 1; j <= psi2.N; ++j)
			{
				double val = psi1.ValueRe[BlokElement(i)] * psi2.ValueRe[BlokElement(j)];
				std::stringstream ss;
				char *det1, *det2;
				uint8_t *detNew;
				det1 = DetToString(psi1.NFermion, psi1.NBoson, &psi1.Determinant[BlokElementDet(i, psi1.NWord, 0)]);
				det2 = DetToString(psi2.NFermion, psi2.NBoson, &psi2.Determinant[BlokElementDet(j, psi2.NWord, 0)]);
				ss << det1 << det2;
				detNew = StringToDet(out.NFermion, out.NBoson, ss.str().c_str());
				RealWaveFunctionAddElement(&out, detNew, val);
				free(det1);
				free(det2);
				free(detNew);
			}
		}
	}
	return out;
}

static int LuaTensorProduct(lua_State *L)
{
	WaveFunctionType psi1, psi2;
	uint32_t len = (uint32_t) lua_gettop(L);
	if (len != 2) return luaL_error(L, "Error: TensorProduct called with %d arguments, exactly 2 expected.", len);
	if (luaL_isType(L, 1, psi1))
	{
		if (!luaL_isType(L, 2, psi2)) return luaL_error(L, "Error in TensorProduct: First argument is of type Wavefunction, second argument is not.");
		lua_check(L, 1, psi1);
		lua_check(L, 2, psi2);
		if (psi1.NBoson != 0 || psi2.NBoson != 0) return luaL_error(L, "One of the Wavefunctions contains bosonic states, this is not yet implemented.");
		lua_push(L, WaveFunctionTensorProduct(psi1, psi2));
	}
	else
	{
		return luaL_error(L, "Error: TensorProduct called with unknown or unimplemented type.");
	}
	return 1;
}
