/*
 * SparseOperator.cpp
 *
 *  Created on: Nov 19, 2018
 *      Author: brass
 */

uint32_t OperatorToDeterminantMatrix(const OperatorType * const O, const WaveFunctionType * const psi, LinAlg::CSRMatrix<double> &M)
{
  std::vector<std::map<uint32_t,double> > entries = std::vector<std::map<uint32_t,double> >(psi->N);
  if((O->QComplex!=0))
  {
	  printf("error in OperatorToDeterminatMatrix: Sparse Matrix is real but Operator is complex\n");
	  return 1;
  }
#pragma omp parallel
{
  uint8_t *Det = (uint8_t*)malloc(psi->NWord * sizeof(uint8_t));
  if (Det == NULL) AllocFailed();
#pragma omp for
  for(uint32_t i=1;i<=psi->N;i++)
  {

	uint32_t ii = i >> NBitsBlok;
	uint32_t ik = i & BlokMask;
	uint32_t ij = psi->NWord * ik;
	for (uint32_t l = 0; l <= O->MaxLength; l++)
	{

	  for (uint32_t j = 0; j < O->l[l].N; j++)
	  {
		for (uint32_t k = 0; k < psi->NWord; k++) Det[k] = psi->Determinant[ii][ij + k];
		double phase = OperatorLader(Det, &(O->l[l].o[l * j]), l, O->NFermion);
		if (phase && O->l[l].VRe[j])
		{
		  uint32_t m=WaveFunctionIndexDet(psi, Det);
		  if(m!=0)
		  {
			  entries[i-1][m-1] += phase * O->l[l].VRe[j];
		  }
		}
	  }

	}
    //M.AddRow(entries);
    //entries.clear();
  }
  free(Det);
} 				// end of parallel region

  for (uint32_t i = 0; i < entries.size(); i++)
  {
	  M.AddRow(entries[i]);
  }
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
