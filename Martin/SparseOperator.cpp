/*
 * SparseOperator.cpp
 *
 *  Created on: Nov 19, 2018
 *      Author: brass
 */

uint32_t OperatorToSparseMatrix(const OperatorType * const O, const WaveFunctionType * const psi, LinAlg::CSRMatrix<double> &M)
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

template<typename F, typename V>
void OrthogonalizeVecWrtSet(F (&SkalarProd)(const V&,const V&), void (&VecAdd)(const V&, V&, F),V &x, const std::vector<V> &vecSet, uint32_t N)
{
	for(uint32_t i = 0; i < N; i++)
	{
		F olapp = SkalarProd(vecSet[i],x);
		VecAdd(vecSet[i],x,-olapp);
	}
}

template<typename F, typename V>
void Lanczos(LinAlg::LinearMap<F, V, V> &LinMap,F (&SkalarProd)(const V&,const V&), void (&VecAdd)(const V&, V&, F), void (&VecMul)(V&, F),const V &startvec, const uint32_t N, std::vector<double> &a, std::vector<double> &b, std::vector<V> &T)
{
	if (a.size()+b.size()+T.size() != 0)
	{
		printf("Error in Lanczos: a, b or T non-empty\n");
		return;
	}
	a.reserve(N);
	b.reserve(N);
	T = std::vector<V>(N+1);
	T[0] = startvec;

	b.push_back(sqrt(SkalarProd(T[0],T[0])));
	VecMul(T[0],1/b[0]);
	LinMap.Evaluate(T[0],T[1]);
	a.push_back(SkalarProd(T[0],T[1]));
	VecAdd(T[0],T[1],-a[0]);


	for (uint32_t i = 1; i < N; i++)
	{
		b.push_back(sqrt(SkalarProd(T[i],T[i])));
		if(b[i] > 100 * DBL_MIN + DBL_EPSILON * b[1])
		{
			VecMul(T[i],1/b[i]);
		}
		else
		{
			T.resize(i);
			return;
		}
		LinMap.Evaluate(T[i],T[i+1]);
		OrthogonalizeVecWrtSet(SkalarProd,VecAdd,T[i+1],T,i);
		a.push_back(SkalarProd(T[i],T[i+1]));
		VecAdd(T[i],T[i+1],-a[i]);
		//OrthogonalizeVecWrtSet(SkalarProd,VecAdd,T[i+1],T,i+1);

	}
}

double SkalarProdukt(const std::vector<double> &x, const std::vector<double> &y)
{
	if (x.size()!=y.size())
	{
		printf("Error: SkalarProdukt of vectors of unequal lengths %lu and %lu\n",x.size(),y.size());
		return NAN;
	}
	else
	{
		double val = 0;
		for (uint32_t i = 0; i < x.size(); i++)
		{
			val += x[i]*y[i];
		}
		return val;
	}
}

template<typename F>
void VectorMul(std::vector<F> &x, F a)
{
	for (uint32_t i = 0; i < x.size(); i++)
	{
		x[i] *= a;
	}
}

template<typename F>
void VectorAdd(const std::vector<F> & x, std::vector<F> & y, F a)
{
	if (x.size()!=y.size())
	{
		printf("Error: VectorAdd of vectors of unequal lengths %lu and %lu\n",x.size(),y.size());
		return;
	}
	else
	{
		for (uint32_t i = 0; i < x.size(); i++)
		{
			y[i] += a * x[i];
		}
	}
}
