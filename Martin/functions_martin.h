/*
 * functions_martin.h
 *
 *  Created on: Aug 31, 2018
 *      Author: brass
 */

#ifndef MARTIN_FUNCTIONS_MARTIN_H_
#define MARTIN_FUNCTIONS_MARTIN_H_

static int LuaOperatorNewFullCoulomb (lua_State *L);

static int LuaOperatorNewDirac(lua_State *L);

static int LuaCreateBaseForContinuum(lua_State *L);

static int LuaPropagator(lua_State *L);

static int LuaPropagator2(lua_State *L);

static int LuaReadFPLOBasis(lua_State *L);

static int LuaMetricWaveFuns(lua_State *L);

#endif /* MARTIN_FUNCTIONS_MARTIN_H_ */
