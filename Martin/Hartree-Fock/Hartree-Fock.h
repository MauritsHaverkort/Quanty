/*
 * Hartree-Fock.h
 *
 *  Created on: Oct 16, 2019
 *      Author: brass
 */

#ifndef MARTIN_HARTREE_FOCK_HARTREE_FOCK_H_
#define MARTIN_HARTREE_FOCK_HARTREE_FOCK_H_

/*
 * Hartree-Fock.cpp
 *
 *  Created on: Sep 26, 2019
 *      Author: brass
 */

namespace HartreeFock
{

double OverlapMatrix(uint32_t i, uint32_t j, uint32_t n, const vector<double> &x);
void CreateOverlapMatrix(const vector<double> &x, const uint32_t n, vector<double> &S);
double LaplaceMatrix(uint32_t i, uint32_t j, uint32_t n, const vector<double> &x);
void CreateLaplacian(const vector<double> &x, const uint32_t n, vector<double> &L);
double PotentialMatrix(uint32_t i, uint32_t j, uint32_t n, const vector<double> &x, double Z, int l);
void CreatePotentialMatrix(const vector<double> &x, const uint32_t n, vector<double> &V, double Z, int l);
void SolveCatOwnerEquation(vector<double> &H, vector<double> &S, uint32_t dim, uint32_t bandsize, vector<double> &eigVal,vector<double> &eigVec);
void CreateGrid(unsigned int Nknots, int n, double R, int powGrid, vector<double> &x);
void PutBSplinesOnQuadGrid(const vector<double> &x, uint32_t n,uint32_t Nweights, vector<double> &BSval);
void CalcOffDiagSlaterInts(const vector<double> &x, vector<double> &rk1, vector<double> &rk2, int k, uint32_t n);
void CalcOffDiagSlaterInts(const vector<double> &x, const vector<double> &BSval,uint32_t Nweights, vector<double> &rk1, vector<double> &rk2, int k, uint32_t n);
void CalcDiagSlaterInts(const vector<double> &x, vector<double> &rkD, int k, uint32_t n);
void CalcDiagSlaterInts(const vector<double> &x, const vector<double> &BSval,uint32_t Nweights, vector<double> &rkD, int k, uint32_t n);
double AssembleSlaterInts(uint32_t i1, uint32_t j1,uint32_t i2, uint32_t j2, const vector<double> &rk1, const vector<double> &rk2, const vector<double> &rkD, uint32_t n);
double AssembleSlaterInts(const vector<double> &Ci1,const vector<double> &Cj1,const vector<double> &Ci2,const vector<double> &Cj2,
						  const vector<double> &rk1, const vector<double> &rk2, const vector<double> &rkD, uint32_t n);
double AssembleSlaterInts2(const vector<double> &Ci1,const vector<double> &Cj1,const vector<double> &Ci2,const vector<double> &Cj2,
						  const vector<double> &rk1, const vector<double> &rk2, const vector<double> &rkD, uint32_t n);
void CreateFockMatrix(const uint32_t la, const vector<vector<double> > &rk1, const vector<vector<double> > &rk2, const vector<vector<double> > &rkD,
					  uint32_t n, const vector<vector<double> > &eigVec, uint32_t dim, const vector<vector<uint32_t> > &Nocc, vector<double> &F, vector<double> &FE);
void CreateHartreeMatrix(const uint32_t la, const vector<double> &rk1, const vector<double> &rk2, const vector<double> &rkD,
					  uint32_t n, const vector<vector<double> > &eigVec, uint32_t dim, const vector<vector<uint32_t> > &Nocc, vector<double> &H, vector<double> &HE);
double CalcGSEnergy(const vector<vector<uint32_t> > &Nocc, const vector<vector<double> > &eigVal, const vector<vector<double> > &eigVec, const vector<vector<double> > &H, uint32_t n);
double TestGSEnergy(const vector<vector<uint32_t> > &Nocc, const vector<vector<double> > &eigVal, const vector<vector<double> > &eigVec, const vector<vector<double> > &H, uint32_t n);
void SCScheme(const vector<double> &x, uint32_t n, int Z,int Nelec, uint32_t NSCmax, double eps, vector<vector<double> > &eigVal,vector<vector<double> > &eigVec);

}

static int LuaAtomicHartreeFock(lua_State *L);




#endif /* MARTIN_HARTREE_FOCK_HARTREE_FOCK_H_ */
