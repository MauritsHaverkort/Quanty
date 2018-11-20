/*
 * SparseOperator.cpp
 *
 *  Created on: Nov 19, 2018
 *      Author: brass
 */

uint32_t OperatorToDeterminantMatrix(const OperatorType * const O, const WaveFunctionType * const psi, LinAlg::CSRMatrix<double> &M)
{
  uint32_t i, j, k, l, m, ii, ij, ik;
  double phase;
  uint8_t *Det;
  std::map<uint32_t,double> entries;
  if((O->QComplex!=0))
  {
	  printf("error in OperatorToDeterminatMatrix: Sparse Matrix is real but Operator is complex\n");
	  return 1;
  }
  Det = (uint8_t*)malloc(psi->NWord * sizeof(uint8_t));
  if (Det == NULL) AllocFailed();
  for(i=1;i<=psi->N;i++)
  {
    ii = i >> NBitsBlok;
    ik = i & BlokMask;
    ij = psi->NWord * ik;
    for (l = 0; l <= O->MaxLength; l++)
    {

	  for (j = 0; j < O->l[l].N; j++)
	  {
		for (k = 0; k < psi->NWord; k++) Det[k] = psi->Determinant[ii][ij + k];
		phase = OperatorLader(Det, &(O->l[l].o[l * j]), l, O->NFermion);
		if (phase && O->l[l].VRe[j])
		{
		  m=WaveFunctionIndexDet(psi, Det);
		  if(m!=0)
		  {
			  entries[m-1] += phase * O->l[l].VRe[j];
		  }
		}
	  }

    }
    M.AddRow(entries);
    entries.clear();
  }
  free(Det);
  return 0;
}

void WaveFunctionToVector(const WaveFunctionType * const psi,std::vector<double> &vecOut,const WaveFunctionType * const indexWavefunction)
{
	vecOut = std::vector<double>(indexWavefunction->N);
#pragma omp parallel for
	for (uint32_t i = 1; i <= psi->N; i++)
	{
		uint32_t ii = i >> NBitsBlok;
		uint32_t ik = i & BlokMask;
		uint32_t ij = psi->NWord * ik;

		uint8_t *Det = (uint8_t*)malloc(psi->NWord * sizeof(uint8_t));
		if (Det == NULL) AllocFailed();
		for (uint32_t k = 0; k < psi->NWord; k++)
		{
			Det[k] = psi->Determinant[ii][ij + k];
		}
		uint32_t m = WaveFunctionIndexDet(indexWavefunction, Det);
		if (m!=0)
		{
			vecOut[m-1] = psi->ValueRe[ii][ij]; //TODO which element?
			printf("vecOut[%d] = %.1f\n",m-1,psi->ValueRe[ii][ij]);
		}
		else
		{
			printf("Warning in WaveFunctionToVector: Determinant not found in indexWavefunction; ignoring determinant\n");
		}

	}
}
