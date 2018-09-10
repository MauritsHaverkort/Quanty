/*
 * Typedefs_Martin.h
 *
 *  Created on: Aug 31, 2018
 *      Author: brass
 */

#ifndef MARTIN_TYPEDEFS_MARTIN_H_
#define MARTIN_TYPEDEFS_MARTIN_H_

typedef struct SQuantumNumbers
{
  unsigned int n[4];
  int kappa[4];

}SQuantumNumbers;

typedef struct SSlaterIntegral
{
  double *val;
  int Nval;
  unsigned int *k;
  SQuantumNumbers tau;

}SSlaterIntegral;

typedef struct{
	unsigned int *indices;
	uint32_t n;
	int kappa;
}SKappaGroup;



#endif /* MARTIN_TYPEDEFS_MARTIN_H_ */
