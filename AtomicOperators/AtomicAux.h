/*
 * AtomicAux.h
 *
 *  Created on: Sep 10, 2018
 *      Author: brass
 */

#ifndef ATOMICOPERATORS_ATOMICAUX_H_
#define ATOMICOPERATORS_ATOMICAUX_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <map>
//#include <time.h>

std::vector<double> ConvertArray(double *a, const int n);

//std::vector<std::vector<double> > PutRelativisticSplineOnGrid(const std::vector<InterpolatingFunctionType *> &funsG, const std::vector<InterpolatingFunctionType *> &funsF, uint32_t NGridPoints);

vector<vector<double> > PutRelativisticSplineOnGrid(const vector<Interpolation::InterpolatingFunction *> &funsG, const vector<Interpolation::InterpolatingFunction *> &funsF, const vector<double> &x,unsigned int k);

std::vector<int> ConvertArray(uint16_t *a, const int n);

std::vector<double> ReadBasisfunction(std::ifstream &file, std::vector<double> &x);

string GetOrbitalName(std::ifstream &file, uint32_t len);

static int LuaReadFPLOBasisFunctions(lua_State *L);

std::vector<double> luaL_checkarray_double(lua_State *L, int narg);

std::vector<int> luaL_checkarray_int(lua_State *L, int narg);

vector<string> luaL_checkarray_string(lua_State *L, int narg);

std::vector<double> operator*(const std::vector<double> &a, const std::vector<double> &b);

static int LuaGetMultipoleBesselIntegral(lua_State *L);

static int LuaGetMultipoleIntegral(lua_State *L);

bool InvolvesOrbital(uint32_t m, uint32_t *orbs, uint32_t N = 4);

#endif /* ATOMICOPERATORS_ATOMICAUX_H_ */
