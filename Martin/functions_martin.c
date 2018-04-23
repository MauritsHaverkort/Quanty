/* Some makeshift functions */

#include "../Martin/FullCoulomb.c"

static int LuaOperatorNewFullCoulomb (lua_State *L)
{
    OperatorType *op;
    int *list, i;
    uint32_t *tuple, Ntuple, NGridpoints, Nindices, temp;
    uint16_t **indices;
    MatrixType *wavefunctions;
    double *x;

    uint32_t NF = (uint32_t)luaL_checkunsigned(L, 1);
    wavefunctions = luaL_checkMatrix(L,2);
    x = luaL_checkarray_double(L,3,&NGridpoints);

    if (wavefunctions->Nj != 2*NGridpoints)
    {
    	printf("Error: relativistic wavefunctions do not match grid. Relativistic wavefunctions need as many values as twice the grid points");
    	fflush(stdout);
    	return -1;
    }

    tuple = luaL_checkarray_unit32_t(L,4, &Ntuple);

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

    CreateOperatorFullCoulombRelativistic(NF,wavefunctions->VRe,wavefunctions->Ni,x,NGridpoints,tuple, Ntuple,indices, op);

    luaL_setmetatable(L,"Operator_Type");

    return 1;
}

static const struct luaL_Reg Martin_meta[] = {
  {NULL, NULL}
};

static const struct luaL_Reg Martin_func[] = {
  {"NewCoulomb", LuaOperatorNewFullCoulomb},
  {NULL, NULL}
};
