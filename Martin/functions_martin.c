/* Some makeshift functions */

#include "../Martin/FullCoulomb.c"

static int LuaOperatorNewFullCoulomb (lua_State *L)
{
    OperatorType *op;
    int *list, i;
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

    NIntegrals = lua_objlen(L,4);
    tuple = (uint32_t **) malloc(NIntegrals*sizeof(uint32_t *));
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
    Nindices = lua_objlen(L,5);

    indices = (uint16_t **) malloc(Nindices * sizeof(uint16_t *));

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

    return 1;
}

static const struct luaL_Reg Martin_meta[] = {
  {NULL, NULL}
};

/*
 * Martin.NewCoulomb(NF, basisfuns,x,tupel,indices,...)
 * 	creates all Coulomb terms specified by tupel and calculates
 * 	the needed Slater integrals for the basisfuns on grid x,
 *
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
 *
 *
 *
 */

static const struct luaL_Reg Martin_func[] = {
  {"NewCoulomb", LuaOperatorNewFullCoulomb},
  {NULL, NULL}
};
