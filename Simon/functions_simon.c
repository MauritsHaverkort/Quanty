/* Some makeshift functions */





void debug(char* debugMessage)
{
	printf("\n\nDEBUG: ");
	printf(debugMessage);
	printf("\n\n");
	fflush(stdout);
}








//should go to MeanField.c
uint32_t MeanFieldOperatorSimon(const OperatorType * const O, OperatorType * const OMF, const MatrixType Rho)
{
	uint32_t i, j, k, l, m, n;
	double OVRe, OVIm, OMFVRe, OMFVIm;
	uint16_t oadd[2] = { 0, 0 };
	if ((O->MaxLength > 4) || ((O->MaxLength > 0) && (O->l[1].N != 0))
			|| ((O->MaxLength > 2) && (O->l[3].N != 0)))
	{
		printf(
				"ERROR: MeanFieldOperator only implemented to reduce a 2, 1 and 0 particle operator to a 1 particle operator\n");
		return 1;
	}
	if (Rho.Ni != Rho.Nj)
	{
		printf(
				"ERROR: Density matrix should be a square matrix in MeanFieldOperatorOneElectronMatrix\n");
		return 1;
	}
	if (Rho.Ni != O->NFermion)
	{
		printf(
				"ERROR: Density matrix should have dimensions of NFermion by NFermion of operator in MeanFieldOperatorOneElectronMatrix\n");
		return 1;
	}
	if (O->MaxLength < 3)
	{
		OperatorCopy(OMF,*O);
		return 0;
	}
	/* initiallize */
	InitOperatorCopySettings(OMF, *O);
	OperatorFixedLengthCopy(&(OMF->l[0]),O->l[0],0);
	OperatorFixedLengthCopy(&(OMF->l[2]),O->l[2],2);
	strncpy(OMF->Name, O->Name, MAXStrLength); //TODO: In C++, add "Mean Field" to the name.
	OMF->Name[MAXStrLength - 1] = '\0';
	//TODO: should I check if it is indeed a Creation-Annihilation pair?
	if (O->l[4].QComplex)
	{
		if (Rho.QComplex)
		{
			for (m = 0; m < O->l[4].N; m++)
			{
				i = NotCrAnMask & O->l[4].o[4 * m + 3];
				j = NotCrAnMask & O->l[4].o[4 * m + 2];
				k = O->l[4].o[4 * m + 1];
				l = O->l[4].o[4 * m];
				OVRe = O->l[4].V[m].Re;
				OVIm = O->l[4].V[m].Im;

				oadd[0] = CrAnMask | i;
				oadd[1] = k;
				OMFVRe = -(OVRe * Rho.VRe[j][l] - OVIm * Rho.VIm[j][l]);
				OMFVIm = -(OVIm * Rho.VRe[j][l] + OVRe * Rho.VIm[j][l]);
				OperatorAddLaderNormalOrder(OMF, oadd, 2, OMFVRe, OMFVIm,
						&QDetNotOrderedFcFaBcBa);

				oadd[0] = CrAnMask | i;
				oadd[1] = l;
				OMFVRe = (OVRe * Rho.VRe[j][k] - OVIm * Rho.VIm[j][k]);
				OMFVIm = (OVIm * Rho.VRe[j][k] + OVRe * Rho.VIm[j][k]);
				OperatorAddLaderNormalOrder(OMF, oadd, 2, OMFVRe, OMFVIm,
						&QDetNotOrderedFcFaBcBa);

				oadd[0] = CrAnMask | j;
				oadd[1] = k;
				OMFVRe = (OVRe * Rho.VRe[i][l] - OVIm * Rho.VIm[i][l]);
				OMFVIm = (OVIm * Rho.VRe[i][l] + OVRe * Rho.VIm[i][l]);
				OperatorAddLaderNormalOrder(OMF, oadd, 2, OMFVRe, OMFVIm,
						&QDetNotOrderedFcFaBcBa);

				oadd[0] = CrAnMask | j;
				oadd[1] = l;
				OMFVRe = -(OVRe * Rho.VRe[i][k] - OVIm * Rho.VIm[i][k]);
				OMFVIm = -(OVIm * Rho.VRe[i][k] + OVRe * Rho.VIm[i][k]);
				OperatorAddLaderNormalOrder(OMF, oadd, 2, OMFVRe, OMFVIm,
						&QDetNotOrderedFcFaBcBa);

				OMFVRe = (Rho.VRe[i][k] * Rho.VRe[j][l] - Rho.VRe[i][l] * Rho.VRe[j][k]
						- Rho.VIm[i][k] * Rho.VIm[j][l] + Rho.VIm[i][l] * Rho.VIm[j][k])
						* OVRe
						- (Rho.VIm[i][k] * Rho.VRe[j][l] - Rho.VIm[i][l] * Rho.VRe[j][k]
								+ Rho.VRe[i][k] * Rho.VIm[j][l] - Rho.VRe[i][l] * Rho.VIm[j][k])
								* OVIm;
				OMFVIm = (Rho.VRe[i][k] * Rho.VRe[j][l] - Rho.VRe[i][l] * Rho.VRe[j][k]
						- Rho.VIm[i][k] * Rho.VIm[j][l] + Rho.VIm[i][l] * Rho.VIm[j][k])
						* OVIm
						+ (Rho.VIm[i][k] * Rho.VRe[j][l] - Rho.VIm[i][l] * Rho.VRe[j][k]
								+ Rho.VRe[i][k] * Rho.VIm[j][l] - Rho.VRe[i][l] * Rho.VIm[j][k])
								* OVRe;
				OperatorAddLaderNormalOrder(OMF,NULL,0,OMFVRe,OMFVIm,&QDetNotOrderedFcFaBcBa);
			}
		} else // Rho real O complex OMF complex
		{
			for (m = 0; m < O->l[4].N; m++)
			{
				i = NotCrAnMask & O->l[4].o[4 * m + 3];
				j = NotCrAnMask & O->l[4].o[4 * m + 2];
				k = O->l[4].o[4 * m + 1];
				l = O->l[4].o[4 * m];
				OVRe = O->l[4].V[m].Re;
				OVIm = O->l[4].V[m].Im;

				oadd[0] = CrAnMask || i;
				oadd[1] = k;
				OMFVRe = -(OVRe * Rho.VRe[j][l]);
				OMFVIm = -(OVIm * Rho.VRe[j][l]);
				OperatorAddLaderNormalOrder(OMF, oadd, 2, OMFVRe, OMFVIm,
						&QDetNotOrderedFcFaBcBa);

				oadd[0] = CrAnMask | i;
				oadd[1] = l;
				OMFVRe = (OVRe * Rho.VRe[j][k]);
				OMFVIm = (OVIm * Rho.VRe[j][k]);
				OperatorAddLaderNormalOrder(OMF, oadd, 2, OMFVRe, OMFVIm,
						&QDetNotOrderedFcFaBcBa);

				oadd[0] = CrAnMask | j;
				oadd[1] = k;
				OMFVRe = (OVRe * Rho.VRe[i][l]);
				OMFVIm = (OVIm * Rho.VRe[i][l]);
				OperatorAddLaderNormalOrder(OMF, oadd, 2, OMFVRe, OMFVIm,
						&QDetNotOrderedFcFaBcBa);

				oadd[0] = CrAnMask | j;
				oadd[1] = l;
				OMFVRe = -(OVRe * Rho.VRe[i][k]);
				OMFVIm = -(OVIm * Rho.VRe[i][k]);
				OperatorAddLaderNormalOrder(OMF, oadd, 2, OMFVRe, OMFVIm,
						&QDetNotOrderedFcFaBcBa);

				OMFVRe = (Rho.VRe[i][k] * Rho.VRe[j][l] - Rho.VRe[i][l] * Rho.VRe[j][k])
						* OVRe;
				OMFVIm = (Rho.VRe[i][k] * Rho.VRe[j][l] - Rho.VRe[i][l] * Rho.VRe[j][k])
						* OVIm;
				OperatorAddLaderNormalOrder(OMF, NULL, 0, OMFVRe, OMFVIm,
						&QDetNotOrderedFcFaBcBa);
			}
		}
	} else // O->l[4] real
	{
		if (Rho.QComplex)
		{
			for (m = 0; m < O->l[4].N; m++)
			{
				i = NotCrAnMask & O->l[4].o[4 * m + 3];
				j = NotCrAnMask & O->l[4].o[4 * m + 2];
				k = O->l[4].o[4 * m + 1];
				l = O->l[4].o[4 * m];
				OVRe = O->l[4].VRe[m];

				oadd[0] = CrAnMask | i;
				oadd[1] = k;
				OMFVRe = -(OVRe * Rho.VRe[j][l]);
				OMFVIm = -(OVRe * Rho.VIm[j][l]);
				OperatorAddLaderNormalOrder(OMF, oadd, 2, OMFVRe, OMFVIm,
						&QDetNotOrderedFcFaBcBa);

				oadd[0] = CrAnMask | i;
				oadd[1] = l;
				OMFVRe = (OVRe * Rho.VRe[j][k]);
				OMFVIm = (OVRe * Rho.VIm[j][k]);
				OperatorAddLaderNormalOrder(OMF, oadd, 2, OMFVRe, OMFVIm,
						&QDetNotOrderedFcFaBcBa);

				oadd[0] = CrAnMask | j;
				oadd[1] = k;
				OMFVRe = (OVRe * Rho.VRe[i][l]);
				OMFVIm = (OVRe * Rho.VIm[i][l]);
				OperatorAddLaderNormalOrder(OMF, oadd, 2, OMFVRe, OMFVIm,
						&QDetNotOrderedFcFaBcBa);

				oadd[0] = CrAnMask | j;
				oadd[1] = l;
				OMFVRe = -(OVRe * Rho.VRe[i][k]);
				OMFVIm = -(OVRe * Rho.VIm[i][k]);
				OperatorAddLaderNormalOrder(OMF, oadd, 2, OMFVRe, OMFVIm,
						&QDetNotOrderedFcFaBcBa);

				OMFVRe = (Rho.VRe[i][k] * Rho.VRe[j][l] - Rho.VRe[i][l] * Rho.VRe[j][k]
						- Rho.VIm[i][k] * Rho.VIm[j][l] + Rho.VIm[i][l] * Rho.VIm[j][k])
						* OVRe;
				OMFVIm = (Rho.VIm[i][k] * Rho.VRe[j][l] - Rho.VIm[i][l] * Rho.VRe[j][k]
						+ Rho.VRe[i][k] * Rho.VIm[j][l] - Rho.VRe[i][l] * Rho.VIm[j][k])
						* OVRe;
				OperatorAddLaderNormalOrder(OMF, NULL, 0, OMFVRe, OMFVIm,
						&QDetNotOrderedFcFaBcBa);
			}
		} else // Rho real O->l[4] real
		{
			for (m = 0; m < O->l[4].N; m++)
			{
				i = NotCrAnMask & O->l[4].o[4 * m + 3];
				j = NotCrAnMask & O->l[4].o[4 * m + 2];
				k = O->l[4].o[4 * m + 1];
				l = O->l[4].o[4 * m];
				OVRe = O->l[4].VRe[m];
				OMFVIm = 0;

				oadd[0] = CrAnMask | i;
				oadd[1] = k;
				OMFVRe = -(OVRe * Rho.VRe[j][l]);
				OperatorAddLaderNormalOrder(OMF, oadd, 2, OMFVRe, OMFVIm,
						&QDetNotOrderedFcFaBcBa);

				oadd[0] = CrAnMask | i;
				oadd[1] = l;
				OMFVRe = (OVRe * Rho.VRe[j][k]);
				OperatorAddLaderNormalOrder(OMF, oadd, 2, OMFVRe, OMFVIm,
						&QDetNotOrderedFcFaBcBa);

				oadd[0] = CrAnMask | j;
				oadd[1] = k;
				OMFVRe = (OVRe * Rho.VRe[i][l]);
				OperatorAddLaderNormalOrder(OMF, oadd, 2, OMFVRe, OMFVIm,
						&QDetNotOrderedFcFaBcBa);

				oadd[0] = CrAnMask | j;
				oadd[1] = l;
				OMFVRe = -(OVRe * Rho.VRe[i][k]);
				OperatorAddLaderNormalOrder(OMF, oadd, 2, OMFVRe, OMFVIm,
						&QDetNotOrderedFcFaBcBa);

				OMFVRe = (Rho.VRe[i][k] * Rho.VRe[j][l] - Rho.VRe[i][l] * Rho.VRe[j][k])
						* OVRe;
				OperatorAddLaderNormalOrder(OMF, NULL, 0, OMFVRe, OMFVIm,
						&QDetNotOrderedFcFaBcBa);
			}
		}
	}
	OperatorChop(OMF,0);
	return 0;
}




//should go to MeanField.c
uint32_t PartialMeanFieldOperator(const OperatorType * const O, OperatorType * const OPMF, const MatrixType Rho, const uint16_t * const indices, const uint32_t indSetSize)
{
	OperatorCopy(OPMF,*O);
	if (O->MaxLength < 3 || indSetSize < 1)
	{
		return 0;
	}
	OperatorType MeanFieldPartTmp, MeanFieldPart;
	int length = 4; //If ever we can apply mean field to operators of length unequal to 4, this needs to be changed.
	uint16_t oadd[length];
	int firstIndexOfBlock;
	int i,j;
	InitOperator(&MeanFieldPartTmp);
	MeanFieldPartTmp.NFermion=OPMF->NFermion;
	MeanFieldPartTmp.NBoson=OPMF->NBoson;
	InitOperatorFixedLength(&MeanFieldPartTmp,4,8);
	for (i=0;i<length*(OPMF->l[length].N);i++)
	{
		if ( isIntElementOf(NotCrAnMask&OPMF->l[length].o[i], indices, indSetSize) )
		{
			firstIndexOfBlock = i - (i % length);
			for (j=0;j<length;j++)
			{
				oadd[j] = OPMF->l[length].o[firstIndexOfBlock+j];
			}
			if (OPMF->l[length].QComplex)
			{
				ComplexOperatorAddLaderNormalOrder(&MeanFieldPartTmp, oadd, length, OPMF->l[length].V[firstIndexOfBlock/length].Re, OPMF->l[length].V[firstIndexOfBlock/length].Im, &QDetNotOrderedFcFaBcBa);
				OPMF->l[length].V[firstIndexOfBlock/length].Re = 0;
				OPMF->l[length].V[firstIndexOfBlock/length].Im = 0;
			} else
			{
				RealOperatorAddLaderNormalOrder(&MeanFieldPartTmp, oadd, length, OPMF->l[length].VRe[firstIndexOfBlock/length], &QDetNotOrderedFcFaBcBa);
				OPMF->l[length].VRe[firstIndexOfBlock/length] = 0;
			}
			i = firstIndexOfBlock + length - 1; //I don't want to count that block again in case the index appears several times.
		}
	}
	if (MeanFieldPartTmp.l[length].N != 0)
	{
			MeanFieldOperatorSimon(&MeanFieldPartTmp,&MeanFieldPart,Rho);
			OperatorAdd(OPMF,MeanFieldPart,1.0,0.0); //TODO: In C++, overload the Add Operator so that it accepts arguments (OperatorType *, OperatorType)
			OperatorChop(OPMF,0);
	}
	FreeOperator(&MeanFieldPartTmp);
	FreeOperator(&MeanFieldPart);
	return 0;
}



//should go to ssp.c
static int LuaPartialMeanFieldOperator(lua_State *L)
{
	OperatorType *O,*OMF;
	MatrixType *Rho;
	uint16_t * indices;
	uint32_t indSetSize;
	O = (OperatorType *) luaL_checkudata(L, 1, "Operator_Type");
	Rho = (MatrixType *) luaL_checkMatrix(L, 2);
	indices = luaL_checkarray_unit16_t(L,3,&indSetSize);
	uint32_t len = (uint32_t) lua_gettop(L);
	if (len < 3 || len > 3)
		return luaL_error(L, "MeanFieldOperator called with %d arguments. Operator, Density Matrix and List of Indices expected.\n",len);
	int i;
	for (i=0;i<indSetSize;i++)
	{
		if (indices[i]<0 || indices[i]>= O->NFermion)
			printf("\nWarning: Index %d is not within [0,NF-1]=[0,%d].\n", indices[i], O->NFermion-1);
	}
	OMF = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
	PartialMeanFieldOperator(O, OMF, *Rho, indices, indSetSize);
	luaL_setmetatable(L, "Operator_Type");
	strncpy(OMF->Name, "", MAXStrLength-1); //An appropriate Name can be given in Lua.
	return 1;
}


//TODO: nur zu Testzwecken, sollte bald wieder gelöscht werden.
static int LuaMeanFieldOperatorSimon(lua_State *L)
{
  OperatorType *H, *HMF;
  MatrixType *Rho;
  uint32_t len = (uint32_t)lua_gettop(L);
  if(len<2 || len>2)
  {
    return luaL_error(L, "MeanFieldOperator called with %d arguments. Exactly two expected, Operator and Density Matrix.\n",len);
  }
  H = (OperatorType *) luaL_checkudata(L, 1, "Operator_Type");
  Rho = (MatrixType *) luaL_checkMatrix(L, 2);
  if((H->NFermion+H->NBoson)!=Rho->Ni)
  {
    return luaL_error(L, "Error in MeanFieldGroundstate, length of operator (NF %i NB %i) is not equal to length of density matrix (%i)\n",H->NFermion,H->NBoson,Rho->Ni);
  }
  HMF = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  MeanFieldOperatorSimon(H, HMF, *Rho);
  luaL_setmetatable(L, "Operator_Type");
  strncpy(HMF->Name, "", MAXStrLength-1); //An appropriate Name can be given in Lua.
  FreeMatrix(Rho);
  return 1;
}

void PartialMeanFieldEigensystem()
{
	//TODO
}

void EigensystemSimon()
{
	//TODO: ruft partialMeanFieldEigensystem auf, wenn die Restrictions dies diktieren.
}





static int LuaPartialMeanFieldEigensystem(lua_State *L)
{
	//TODO
  return 0;
}

//initialises and modifies a list of indices that can, without any approximation, be treated on a mean field level.
void mandatoryMeanFieldList(uint16_t **indices, int * const indSetSize, const DeterminantRestrictionType det, int NFermion, int NBoson)//TODO: rewrite this in C++, using vectors
{
	uint32_t NWord=NBoson + Min(NFermion, 1 + (NFermion - 1) / 8);
	*indices = (uint16_t*) malloc(det.N * NWord * sizeof(uint16_t));
	int i;
	char * detStr;
	for (i = 0; i < det.N; ++i)
	{
		if (det.Min[i] == det.Max[i])
		{
			detStr = DetToString(NFermion, NBoson, det.Determinant[i]);
			debug(detStr);
//		TODO: add every index with 1 to a temporary vector.
//		TODO: compare number of entries to number of max and min electrons. If equal, add those indices to the return vector that are not already in it.
		}
	}
}


static int LuaOrbToMultiplicity(lua_State * L)
{
	char * orb;
	int multiplicity;
	uint32_t len = (uint32_t)lua_gettop(L);
  if(len != 1)
  {
    return luaL_error(L, "OrbToMultiplicity called with %d arguments. Only name of the orbital expected.\n",len);
  }
  orb = (char*)luaL_checkstring(L,1);
  lua_pop(L,1);
  multiplicity = nonRelAngularMomentumZNumbers(orb);
  if (multiplicity < 0)
  	multiplicity = relAngularMomentumZNumbers(orb);
  if (multiplicity < 0)
  {
  	printf("\nCould not interpret orbital \"");
  	printf(orb);
  	printf("\" as either relativistic or non-relativistic atomic orbital, returning %d.",multiplicity);
  }
  lua_pushnumber(L,multiplicity);
  return 1;
}

static int LuaSimontest(lua_State * L)
{
	DeterminantRestrictionType testDet = *(DeterminantRestrictionType*)luaL_checkarray_DeterminantRestrictionType(L, 1);
	uint16_t * indList;
	int listLength;

//	char * detStr = DetToString(testDet.Determinant[1]);

	mandatoryMeanFieldList(&indList, &listLength, testDet, 3, 0);
  return 0;
}

// accepts a list of Wave Functions as well as the Hamiltonian and a list of other operators. The Hamiltonian is special here because it can be used to calculate the Deviance. It is optional.
static int LuaPrintExpectationValues(lua_State * L)
{
	WaveFunctionType *psiList;
	OperatorType *hamil, *opList;
	int numberOfWavefunctions, numberOfOperators;
	uint32_t len = (uint32_t) lua_gettop(L);
	if (len < 2 || len > 3)
		return luaL_error(L,
						"PrintExpectationValues called with %d arguments. 2 or 3 expected, List of Wavefunctions, optional Hamiltonian, and list of other operators.\n",
						len);
	if(lua_istable(L,1))
	{
		psiList = luaL_checkarray_WaveFunctionType(L, 1, &numberOfWavefunctions);
	} else
	{
		numberOfWavefunctions = 1;
		psiList = (WaveFunctionType *)luaL_checkudata(L,1,"Wavefunction_Type");
	}
	if (len == 2)
		hamil = (OperatorType*)luaL_checkudata(L,2,"Operator_Type");
	if(lua_istable(L,-1))
	{
		opList = luaL_checkarray_OperatorType(L, -1, &numberOfOperators);
	} else
	{
		numberOfOperators = 1;
		opList = (OperatorType *)luaL_checkudata(L,-1,"Operator_Type");
	}

//TODO: weiter, falls BraKet nicht mehr komplett auf Lua Ebene stattfindet


	return 0;
}



static const struct luaL_Reg Simon_meta[] = {
  {NULL, NULL}
};

static const struct luaL_Reg Simon_func[] = {
  {"test", LuaSimontest},
	{"MeanFieldOperatorSimon",LuaMeanFieldOperatorSimon},
	{"OrbToMultiplicity",LuaOrbToMultiplicity},
	{"PartialMeanFieldEigensystem",LuaPartialMeanFieldEigensystem},
	{"PartialMeanFieldOperator",LuaPartialMeanFieldOperator},
	{"PrintExpectationValues",LuaPrintExpectationValues},
  {NULL, NULL}
};
