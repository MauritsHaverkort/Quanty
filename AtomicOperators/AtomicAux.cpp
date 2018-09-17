/*
 * AtomicAux.cpp
 *
 *  Created on: Sep 10, 2018
 *      Author: brass
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <map>

#include "GaussQuadrature.cpp"
#include "SlaterTuple.cpp"
#include "SlaterIntegrals.cpp"
#include "AtomicU.cpp"
#include "AtomicDirac.cpp"

using std::vector;
using std::string;
using std::ifstream;
using std::map;

vector<double> ConvertArray(double *a, const int n)
{
	vector<double> res = vector<double>(n);
	for (int i = 0; i < n; i++)
	{
		res[i] = a[i];
	}
	return res;
}

vector<vector<double> > PutRelativisticSplineOnGrid(vector<Interpolation::InterpolatingFunction *> funsG, vector<Interpolation::InterpolatingFunction *> funsF, vector<double> x,unsigned int k)
{
	vector<vector<double> > funsCombined = vector<vector<double> >(funsG.size());
	unsigned int NGridPoints = (unsigned int) x.size();
	for (uint32_t i = 0; i < funsG.size(); i++)
	{
		funsCombined[i] = vector<double>(2*NGridPoints);
		for (uint32_t j = 0; j < NGridPoints; j++)
		{
			funsCombined[i][j            ] = (*funsG[i])(x[j],k);
			funsCombined[i][j+NGridPoints] = (*funsF[i])(x[j],k);
		}
	}
	return funsCombined;

}

vector<int> ConvertArray(uint16_t *a, const int n)
{
	vector<int> res = vector<int>(n);
	for (int i = 0; i < n; i++)
	{
		res[i] = (short) a[i];
	}
	return res;
}

vector<double> ReadBasisfunction(ifstream &file, vector<double> &x)
{
	double val;
	vector<double> y;
	string word;

	file >> word;
	while (word != "#")
	{
		val = atof(word.c_str());
		if(x.size() == y.size())
			x.push_back(val);
		file >> word;
		val = atof(word.c_str());
		y.push_back(val);
		file >> word;
	}
	while ((word != "<start") && !(file.eof()))
		file >> word;
	return y;
}

string GetOrbitalName(ifstream &file, uint32_t len)
{
	string word, name;
	while ( file >> word)
	{
		if (word == "function")
		{
			file >> word;
			file >> word;
			file >> word;
			file >> name;
			if (len == 3)
			{
				file >> word;
				name.append(word);
			}
			name.erase(name.length()-1);
			return name;
		}
	}
	return "";
}

static int LuaReadFPLOBasisFunctions(lua_State *L)
{
	ifstream file1,file2;

	uint32_t Norbitals;
	uint32_t len = (uint32_t)lua_gettop(L);
	char** orbitalNames = luaL_checkarray_string(L, 1, &Norbitals); //TODO free
	const char *name1 = luaL_checkstring(L,2);
	const char *name2;

	file1.open(name1);
	if(!file1.is_open())
	{
		cout << "Error Could not open file: " << name1 << endl;
		return 0;
	}
	string line, word, name;
	vector<string> orbitals;
	vector<double> x, y1,y2;
	std::map<string,vector<double> > funsG, funsF;
	vector<Interpolation::CubicSpline*> splinesG, splinesF;

	// Read raw data
	while (!file1.eof())
	{
		name = GetOrbitalName(file1,len);
		y1 = ReadBasisfunction(file1,x);

		funsG[name] = y1;
	}
	file1.close();

	splinesG = vector<Interpolation::CubicSpline*>(Norbitals);

	// Push large part on stack
	lua_createtable(L,0,Norbitals);
	for (uint32_t i = 0; i < Norbitals; i++)
	{
		string charTostring(orbitalNames[i]);
		if (funsG.find(charTostring) != funsG.end())
		{
			void* memory = lua_newuserdata(L, sizeof(Interpolation::CubicSpline));
			splinesG[i] = new(memory) Interpolation::CubicSpline(x,funsG[charTostring]);
			luaL_setmetatable(L, "InterpolatingFunction_Type");
			lua_rawseti(L,-2,i+1);
		}
		else
		{
			printf("LuaReadFPLOBasis: Couldn't find ");
			printf("%s\n", orbitalNames[i]);
			fflush(stdout);
		}

	}
	//Read small part
	if (len == 3)
	{
		name2 = luaL_checkstring(L,3);
		file2.open(name2);
		if(!file2.is_open())
		{
			cout << "Error Could not open file: " << name2 << endl;
			return 0;
		}

		while (!file2.eof())
		{
			name = GetOrbitalName(file2,len);
			y2 = ReadBasisfunction(file2,x);

			funsF[name] = y2;
		}
		file2.close();

		splinesF = vector<Interpolation::CubicSpline*>(Norbitals);

		// Push small part on stack
		lua_createtable(L,0,Norbitals);
		for (uint32_t i = 0; i < Norbitals; i++)
		{
			string charTostring(orbitalNames[i]);
			if (funsF.find(charTostring) != funsF.end())
			{
				void* memory = lua_newuserdata(L, sizeof(Interpolation::CubicSpline));
				splinesF[i] = new(memory) Interpolation::CubicSpline(x,funsF[charTostring]);
				luaL_setmetatable(L, "InterpolatingFunction_Type");
				lua_rawseti(L,-2,i+1);
			}
			else
			{
				printf("LuaReadFPLOBasis: Couldn't find ");
				printf("%s\n", orbitalNames[i]);
				fflush(stdout);
			}

		}
		return 2;
	}
/*	for( int i = 0; i < Norbitals; i++)
	{
		free(orbitalNames[i]);
	}
	free(orbitalNames);*/

	return 1;
}

vector<double> luaL_checkarray_double(lua_State *L, int narg)
{
  luaL_checktype(L, narg, LUA_TTABLE);
  uint32_t len = (uint32_t)lua_objlen(L, narg);

  vector<double> buff = vector<double>(len);
  for(int i=0;i<len;i++)
  {
    lua_rawgeti(L, narg, i+1);
    buff[i]=(double)luaL_checknumber(L, -1);
    lua_pop(L, 1);
  }
  return buff;
}

vector<int> luaL_checkarray_int(lua_State *L, int narg)
{
  luaL_checktype(L, narg, LUA_TTABLE);
  uint32_t len = (uint32_t)lua_objlen(L, narg);

  vector<int> buff = vector<int>(len);
  for(uint32_t i=0;i<len;i++)
  {
    lua_rawgeti(L, narg, i+1);
    buff[i]=(int)luaL_checkunsigned(L, -1);
    lua_pop(L, 1);
  }
  return buff;
}

vector<string> luaL_checkarray_string(lua_State *L, int narg)
{
  int i;
  luaL_checktype(L, narg, LUA_TTABLE);
  uint32_t len = (uint32_t)lua_objlen(L, narg);

  vector<string> buff = vector<string>(len);
  for(i=0;i<len;i++)
  {
    lua_rawgeti(L, narg, i+1);
    char *cstring = (char*)luaL_checkstring(L,-1);
    buff[i] = string(cstring);
    lua_pop(L, 1);
  }
  return buff;
}

vector<double> operator*(vector<double> a, vector<double> b)
{
	unsigned long n = a.size();
	if (n != b.size())
	{
		printf("error in operator*(vector,vector) : vectors of unequal length\n");
		fflush(stdout);
	}
	vector<double> c = vector<double>(n);
	for(unsigned long i = 0; i < n; i++)
	{
		c[i] = a[i]*b[i];
	}
	return c;
}


