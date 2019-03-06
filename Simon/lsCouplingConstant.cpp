double lsdCouplingConstant(const string atom, const uint occupation)
{
	uint atomNum = AtomNameToNumber(atom.c_str());
	switch (atomNum)
	{
	case 19: //K
		switch (occupation)
		{
		case 1:
			return 0.;
		default:
			cout << "Error in CouplingConstant: Occupation of " << occupation << " electrons in d-shell of " << atom << " not supported, returning 0." << endl;
			return 0.;
		}
		break;
	case 20: //Ca
		switch (occupation)
		{
		case 2:
			return 0.001;
		case 1:
			return 0.004;
		default:
			cout << "Error in CouplingConstant: Occupation of " << occupation << " electrons in d-shell of " << atom << " not supported, returning 0." << endl;
			return 0.;
		}
		break;
	case 21: //Sc
		switch (occupation)
		{
		case 3:
			return 0.005;
		case 2:
			return 0.007;
		case 1:
			return 0.010;
		default:
			cout << "Error in CouplingConstant: Occupation of " << occupation << " electrons in d-shell of " << atom << " not supported, returning 0." << endl;
			return 0.;
		}
		break;
	case 22: //Ti
		switch (occupation)
		{
		case 4:
			return 0.009;
		case 3:
			return 0.012;
		case 2:
			return 0.016;
		case 1:
			return 0.019;
		default:
			cout << "Error in CouplingConstant: Occupation of " << occupation << " electrons in d-shell of " << atom << " not supported, returning 0." << endl;
			return 0.;
		}
		break;
	case 23: //V
		switch (occupation)
		{
		case 5:
			return 0.014;
		case 4:
			return 0.018;
		case 3:
			return 0.022;
		case 2:
			return 0.027;
		case 1:
			return 0.031;
		default:
			cout << "Error in CouplingConstant: Occupation of " << occupation << " electrons in d-shell of " << atom << " not supported, returning 0." << endl;
			return 0.;
		}
		break;
	case 24: //Cr
		switch (occupation)
		{
		case 6:
			return 0.021;
		case 5:
			return 0.025;
		case 4:
			return 0.030;
		case 3:
			return 0.035;
		case 2:
			return 0.041;
		default:
			cout << "Error in CouplingConstant: Occupation of " << occupation << " electrons in d-shell of " << atom << " not supported, returning 0." << endl;
			return 0.;
		}
		break;
	case 25: //Mn
		switch (occupation)
		{
		case 7:
			return 0.032;
		case 6:
			return 0.035;
		case 5:
			return 0.040;
		case 4:
			return 0.046;
		case 3:
			return 0.052;
		default:
			cout << "Error in CouplingConstant: Occupation of " << occupation << " electrons in d-shell of " << atom << " not supported, returning 0." << endl;
			return 0.;
		}
		break;
	case 26: //Fe
		switch (occupation)
		{
		case 8:
			return 0.021;
		case 7:
			return 0.045;
		case 6:
			return 0.052;
		case 5:
			return 0.059;
		case 4:
			return 0.066;
		default:
			cout << "Error in CouplingConstant: Occupation of " << occupation << " electrons in d-shell of " << atom << " not supported, returning 0." << endl;
			return 0.;
		}
		break;
	case 27: //Co
		switch (occupation)
		{
		case 9:
			return 0.027;
		case 8:
			return 0.059;
		case 7:
			return 0.066;
		case 6:
			return 0.074;
		case 5:
			return 0.082;
		default:
			cout << "Error in CouplingConstant: Occupation of " << occupation << " electrons in d-shell of " << atom << " not supported, returning 0." << endl;
			return 0.;
		}
		break;
	case 28: //Ni
		switch (occupation)
		{
		case 9:
			return 0.074;
		case 8:
			return 0.083;
		case 7:
			return 0.091;
		case 6:
			return 0.101;
		default:
			cout << "Error in CouplingConstant: Occupation of " << occupation << " electrons in d-shell of " << atom << " not supported, returning 0." << endl;
			return 0.;
		}
		break;
	case 29: //Cu
		switch (occupation)
		{
		case 9:
			return 0.102;
		case 8:
			return 0.112;
		case 7:
			return 0.123;
		default:
			cout << "Error in CouplingConstant: Occupation of " << occupation << " electrons in d-shell of " << atom << " not supported, returning 0." << endl;
			return 0.;
		}
		break;
	case 30: //Zn
		switch (occupation)
		{
		case 9:
			return 0.136;
		case 8:
			return 0.147;
		default:
			cout << "Error in CouplingConstant: Occupation of " << occupation << " electrons in d-shell of " << atom << " not supported, returning 0." << endl;
			return 0.;
		}
		break;
	default:
		cout << "Error in CouplingConstant: Atom with name " << atom << " not included, returning 0." << endl;
		return 0.;
	}
}

double lspCouplingConstant(const string atom, const uint dShellOccupation)
{
	uint atomNum = AtomNameToNumber(atom.c_str());
	switch (atomNum)
	{
	case 19: //K
		switch (dShellOccupation)
		{
		case 2:
			return 1.880;
		case 1:
			return 1.880;
		default:
			cout << "Error in CouplingConstant: Occupation of " << dShellOccupation << " electrons in d-shell of " << atom << " not supported, returning 0." << endl;
			return 0.;
		}
		break;
	case 20: //Ca
		switch (dShellOccupation)
		{
		case 3:
			return 2.405;
		case 2:
			return 2.405;
		case 1:
			return 2.404;
		default:
			cout << "Error in CouplingConstant: Occupation of " << dShellOccupation << " electrons in d-shell of " << atom << " not supported, returning 0." << endl;
			return 0.;
		}
		break;
	case 21: //Sc
		switch (dShellOccupation)
		{
		case 4:
			return 3.034;
		case 3:
			return 3.033;
		case 2:
			return 3.032;
		case 1:
			return 3.032;
		default:
			cout << "Error in CouplingConstant: Occupation of " << dShellOccupation << " electrons in d-shell of " << atom << " not supported, returning 0." << endl;
			return 0.;
		}
		break;
	case 22: //Ti
		switch (dShellOccupation)
		{
		case 5:
			return 3.778;
		case 4:
			return 3.777;
		case 3:
			return 3.776;
		case 2:
			return 3.776;
		case 1:
			return 3.776;
		default:
			cout << "Error in CouplingConstant: Occupation of " << dShellOccupation << " electrons in d-shell of " << atom << " not supported, returning 0." << endl;
			return 0.;
		}
		break;
	case 23: //V
		switch (dShellOccupation)
		{
		case 6:
			return 4.652;
		case 5:
			return 4.651;
		case 4:
			return 4.650;
		case 3:
			return 4.649;
		case 2:
			return 4.650;
		default:
			cout << "Error in CouplingConstant: Occupation of " << dShellOccupation << " electrons in d-shell of " << atom << " not supported, returning 0." << endl;
			return 0.;
		}
		break;
	case 24: //Cr
		switch (dShellOccupation)
		{
		case 7:
			return 5.671;
		case 6:
			return 5.669;
		case 5:
			return 5.668;
		case 4:
			return 5.667;
		case 3:
			return 5.668;
		default:
			cout << "Error in CouplingConstant: Occupation of " << dShellOccupation << " electrons in d-shell of " << atom << " not supported, returning 0." << endl;
			return 0.;
		}
		break;
	case 25: //Mn
		switch (dShellOccupation)
		{
		case 8:
			return 6.849;
		case 7:
			return 6.847;
		case 6:
			return 6.846;
		case 5:
			return 6.845;
		case 4:
			return 6.845;
		default:
			cout << "Error in CouplingConstant: Occupation of " << dShellOccupation << " electrons in d-shell of " << atom << " not supported, returning 0." << endl;
			return 0.;
		}
		break;
	case 26: //Fe
		switch (dShellOccupation)
		{
		case 9:
			return 8.203;
		case 8:
			return 8.202;
		case 7:
			return 8.200;
		case 6:
			return 8.199;
		case 5:
			return 8.199;
		default:
			cout << "Error in CouplingConstant: Occupation of " << dShellOccupation << " electrons in d-shell of " << atom << " not supported, returning 0." << endl;
			return 0.;
		}
		break;
	case 27: //Co
		switch (dShellOccupation)
		{
		case 10:
			return 9.752;
		case 9:
			return 9.750;
		case 8:
			return 9.748;
		case 7:
			return 9.746;
		case 6:
			return 9.746;
		default:
			cout << "Error in CouplingConstant: Occupation of " << dShellOccupation << " electrons in d-shell of " << atom << " not supported, returning 0." << endl;
			return 0.;
		}
		break;
	case 28: //Ni
		switch (dShellOccupation)
		{
		case 10:
			return 11.509;
		case 9:
			return 11.507;
		case 8:
			return 11.506;
		case 7:
			return 11.505;
		default:
			cout << "Error in CouplingConstant: Occupation of " << dShellOccupation << " electrons in d-shell of " << atom << " not supported, returning 0." << endl;
			return 0.;
		}
		break;
	case 29: //Cu
		switch (dShellOccupation)
		{
		case 10:
			return 13.498;
		case 9:
			return 13.496;
		case 8:
			return 13.495;
		default:
			cout << "Error in CouplingConstant: Occupation of " << dShellOccupation << " electrons in d-shell of " << atom << " not supported, returning 0." << endl;
			return 0.;
		}
		break;
	case 30: //Zn
		switch (dShellOccupation)
		{
		case 10:
			return 15.738;
		case 9:
			return 15.737;
		default:
			cout << "Error in CouplingConstant: Occupation of " << dShellOccupation << " electrons in d-shell of " << atom << " not supported, returning 0." << endl;
			return 0.;
		}
		break;
	default:
		cout << "Error in CouplingConstant: Atom with name " << atom << " not included, returning 0." << endl;
		return 0.;
	}
}

static int LuaLSCouplingConstants(lua_State * L)
{
	uint32_t len = (uint32_t) lua_gettop(L);
	if (len != 2)
		return luaL_error(L, "lsCouplingConstants called with %d arguments, exactly 2 expected, Name of atom and intial d-shell occupation.", len);
	string atom;
	uint dOccupation;
	lua_check(L, 1, atom);
	lua_check(L, 2, dOccupation);
	lua_push(L, lspCouplingConstant(atom, dOccupation + 1));
	lua_push(L, lsdCouplingConstant(atom, dOccupation));
	return 2;
}
