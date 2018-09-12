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

std::vector<double> ConvertArray(double *a, const int n);

std::vector<std::vector<double> > PutRelativisticSplineOnGrid(std::vector<InterpolatingFunctionType *> funsG, std::vector<InterpolatingFunctionType *> funsF, uint32_t NGridPoints);

vector<vector<double> > PutRelativisticSplineOnGrid(vector<Interpolation::InterpolatingFunction *> funsG, vector<Interpolation::InterpolatingFunction *> funsF, vector<double> x,unsigned int k);

std::vector<int> ConvertArray(uint16_t *a, const int n);

std::vector<double> ReadBasisfunction(std::ifstream &file, std::vector<double> &x);

string GetOrbitalName(std::ifstream &file);

static int LuaReadFPLOBasisFunctions(lua_State *L);

std::vector<double> luaL_checkarray_double(lua_State *L, int narg);

#endif /* ATOMICOPERATORS_ATOMICAUX_H_ */
