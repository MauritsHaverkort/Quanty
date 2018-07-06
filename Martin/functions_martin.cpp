/* Some makeshift functions */

//#include "SlaterTuple.cpp"
#include "CreateOperators.c"

static int LuaOperatorNewFullCoulomb (lua_State *L)
{
    OperatorType *op;
    int i;
    uint32_t **tuple, NIntegrals, NGridpoints, Nindices, temp, len, NF;
    uint16_t **indices;
    MatrixType *wavefunctions;
    double *x;

    len = (uint32_t)lua_gettop(L);

    if ( len != 5 && len != 6)
    {
    	printf("Error in Martin.NewCoulomb: 5 or 6 arguments expected; got %d",len);
    	return -1;
    }

    NF = (uint32_t)luaL_checkunsigned(L, 1);
    wavefunctions = luaL_checkMatrix(L,2);
    x = luaL_checkarray_double(L,3,&NGridpoints);

    if (wavefunctions->Nj != 2*NGridpoints)
    {
    	printf("Error: relativistic wavefunctions do not match grid. Relativistic wavefunctions need as many values as twice the grid points");
    	fflush(stdout);
    	return -1;
    }

    NIntegrals = (uint32_t)lua_objlen(L,4);
    tuple = (uint32_t **) malloc(NIntegrals*sizeof(uint32_t *));
    if(!tuple) printf("malloc fail at tuple in LuaOperatorNewFullCoulomb. size = %d\n",NIntegrals);
    for ( i = 0; i < NIntegrals; i++)
    {
    	lua_rawgeti(L, 4, i+1);
    	tuple[i] = luaL_checkarray_unit32_t(L,-1,&temp);
    	if ( temp != 8)
    	{
    		printf("Error in Martin.NewCoulomb: Length of tuple[%d] (%d) != 8",i,temp);
    		return -1;
    	}
    	lua_pop(L,1);
    }

    luaL_checktype(L, 5, LUA_TTABLE);
    Nindices = (uint32_t)lua_objlen(L,5);

    indices = (uint16_t **) malloc(Nindices * sizeof(uint16_t *));
    if(!indices) printf("malloc fail at indices in LuaOperatorNewFullCoulomb. size = %d\n",Nindices);

    for( i = 0; i < Nindices; i++ )
    {
    	lua_rawgeti(L, 5, i+1);
    	indices[i] = luaL_checkarray_unit16_t(L,-1,&temp);
    	lua_pop(L,1);
    }



    op = (OperatorType *) lua_newuserdata(L,sizeof(OperatorType));


    if ( len == 5)
    {
    	CreateOperatorFullCoulombRelativistic(NF,wavefunctions->VRe,wavefunctions->Ni,x,NGridpoints,tuple, NIntegrals,indices, op);
    }
    else
    {
    	CreateOperatorFastCoulombRelativistic(NF,wavefunctions->VRe,wavefunctions->Ni,x,NGridpoints,tuple, NIntegrals,indices, op);
    }

    luaL_setmetatable(L,"Operator_Type");

    /*printf("Operaotor op: %p\n",op);
    printf("op.l: %p\t %d\n",op->l,op->MaxLength);
    printf("op.l[4]: %p\t %d\n",&(op->l[4]),op->l[4].N);
    printf("WF: %p\n",wavefunctions);
    printf("x: %p\n", x);
    printf("indices: %p \n",indices);
    printf("tuple: %p\n",tuple);
    printf("\n");*/
	FreeMatrix(wavefunctions);
	free(x);

	for( i = 0; i < Nindices; i++ )
	{
		free(indices[i]);
	}
	free(indices);
	for( i = 0; i < NIntegrals; i++ )
	{
		free(tuple[i]);
	}
	free(tuple);
	/**/

    return 1;
}

static int LuaOperatorNewDirac(lua_State *L)
{
	OperatorType *op;
	uint32_t **indices;
	uint32_t NF, i, NGroups, NGridpoints, NIndices, NBasisFuns, temp;
	SKappaGroup *groups;
	MatrixType *basisfuns, *derivfuns;
	int Z;
	double *x;

	NF = (uint32_t)luaL_checkunsigned(L, 1);
	basisfuns = luaL_checkMatrix(L,2);
	derivfuns = luaL_checkMatrix(L,3);
	x = luaL_checkarray_double(L,4,&NGridpoints);
	luaL_checktype(L, 5, LUA_TTABLE);
	NGroups = (uint32_t)lua_objlen(L,5);

	groups = (SKappaGroup *) malloc(NGroups * sizeof(SKappaGroup));
	if (!groups) printf("malloc fail at groups in LuaOperatorNewDirac. size = %d\n",NGroups);

	for( i = 0; i < NGroups; i++ )
	{
		lua_rawgeti(L, 5, i+1);
	   	lua_rawgeti(L, -1, 1);
	   	groups[i].kappa = (int) luaL_checkunsigned(L, -1);
	   	lua_pop(L,1);
	   	lua_rawgeti(L, -1, 2);
		groups[i].indices = luaL_checkarray_unit32_t(L,-1,&(groups[i].n));
	   	lua_pop(L,1);
	   	lua_pop(L,1);
	}
	Z = (int) luaL_checkunsigned(L, 6);
	luaL_checktype(L, 7, LUA_TTABLE);
	NIndices = (uint32_t)lua_objlen(L,7);

	indices = (uint32_t **) malloc(NIndices * sizeof(uint32_t *));
	if(!indices) printf("malloc fail at indices in LuaOperatorNewDirac. size = %d\n",NIndices);
	for( i = 0; i < NIndices; i++ )
	{
	   	lua_rawgeti(L, 7, i+1);
	   	indices[i] = luaL_checkarray_unit32_t(L,-1,&temp);
	   	lua_pop(L,1);
	}

	op = (OperatorType *) lua_newuserdata(L,sizeof(OperatorType));
	CreateOperatorDirac(NF,basisfuns->Ni, basisfuns->VRe,derivfuns->VRe,x,NGridpoints,groups,NGroups,indices,Z,op);

	luaL_setmetatable(L,"Operator_Type");

	FreeMatrix(basisfuns);
	FreeMatrix(derivfuns);
	free(x);
	free(groups);
	for( i = 0; i < NIndices; i++ )
	{
		free(indices[i]);
	}
	free(indices);

	//Todo free Memory !!!!! same for Slaters

	return 1;
}

static int LuaCreateBaseForContinuum(lua_State *L)
{
	uint32_t NGridpoints = 0;
	SKappaGroup group;
	MatrixType *basisfuns, *derivfuns, *contfuns;
	double w;
	double *x;

	printf("w%e\n",w);
	basisfuns = luaL_checkMatrix(L,1);
	derivfuns = luaL_checkMatrix(L,1);
	x = luaL_checkarray_double(L,3,&NGridpoints);
	luaL_checktype(L, 4, LUA_TTABLE);
	/*NGroups = (uint32_t)lua_objlen(L,4);

	groups = (SKappaGroup *) malloc(NGroups * sizeof(SKappaGroup));
	if (!groups) printf("malloc fail at groups in LuaCreateBaseForContinuum. size = %d\n",NGroups);

	for( i = 0; i < NGroups; i++ )
	{
		lua_rawgeti(L, 4, i+1);
	   	lua_rawgeti(L, -1, 1);
	   	groups[i].kappa = (int) luaL_checkunsigned(L, -1);
	   	lua_pop(L,1);
	   	lua_rawgeti(L, -1, 2);
		groups[i].indices = luaL_checkarray_unit32_t(L,-1,&(groups[i].n));
	   	lua_pop(L,1);
	   	lua_pop(L,1);
	}
	*/
	lua_rawgeti(L, 4, 1);
	group.kappa = (int) luaL_checkunsigned(L, -1);
	lua_pop(L,1);
	lua_rawgeti(L, 4, 2);
	group.indices = luaL_checkarray_unit32_t(L,-1,&(group.n));
	lua_pop(L,1);
	w = luaL_checknumber(L, 5);
	luaL_checktype(L, 6, LUA_TTABLE);

	contfuns = luaL_checkMatrix(L,6);


	double *eigenval;
//	printf("CreateCont\n");
//	fflush(stdout);
	eigenval = CreateBaseForContinuum(basisfuns->Ni, basisfuns->VRe,derivfuns->VRe,x,NGridpoints,group,contfuns->Ni,contfuns->VRe,w);
//	printf("Push\n");
//	fflush(stdout);
	LuaPushArray_double(L,&(eigenval[contfuns->Ni/2]), contfuns->Ni/2);
	contfuns->VRe = &(contfuns->VRe[contfuns->Ni/2]);
	contfuns->Ni = contfuns->Ni/2;
	LuaPushTable_MatrixType(L,contfuns);

//	printf("Free\n");
//	fflush(stdout);

	FreeMatrix(basisfuns);
	FreeMatrix(derivfuns);

	free(x);
	//free(groups);

	return 2;
}

static const struct luaL_Reg Martin_meta[] = {
  {NULL, NULL}
};

/*
 * Martin.NewCoulomb(NF, basisfuns,x,tuple,indices,...)
 * 	creates all Coulomb terms specified by tupel and calculates
 * 	the needed Slater integrals for the basisfuns on grid x,
 *
 * 		NF: number of fermionic modes
 * 		x: list of points on which the radial wave functions are defined
 *
 * 		basisfuns: list of radial wave functions on x (let N be the number of points)
 * 		basisfuns = { {G1(x0),...,G1(xN),F1(x0),...,F1(xN)},
 * 					  {G2(x0),...,G2(xN),F2(x0),...,F2(xN)}, ...}
 * 				where G (F) stands for the large (small) spinor component
 *
 * 		indices: indices of the mj orbitals corresponding to the basisfuns
 * 		indices = { {m1(m=-j),...,m1(m=+j)}, {m2(m=-j),...,m2(m=+j)}, ... }
 * 			where indices[i] corresponds to basisfuns[i]
 * 		tuple: list of 4-tuple of quantum numbers for slater integrals to be used
 * 		tuple = { {n1,k1,n2,k2,n3,k3,n4,k4}, ... }
 *
 * 		if a 6th argument "fast" is specified, tuple MUST NOT contain duplicates in
 * 		the sense that two 4-tuple are equivalent if you can map one onto the other by
 * 		using any combination of the following permutations: { (1,3), (2,4), (1,2)(3,4) }
 * 		this assures that Slater integrals with the same value are calculated once only,
 * 		which leads to an increase in speed.
 *
 *
 *
 */

static const struct luaL_Reg Martin_func[] = {
  {"NewCoulomb", LuaOperatorNewFullCoulomb},
  {"NewDirac", LuaOperatorNewDirac},
  {"ContinuumBase",LuaCreateBaseForContinuum},
  {NULL, NULL}
};
