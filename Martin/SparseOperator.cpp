/*
 * SparseOperator.cpp
 *
 *  Created on: Nov 19, 2018
 *      Author: brass
 */

uint32_t OperatorToSparseMatrix(OperatorType * const O, const WaveFunctionType * const psi, LinAlg::CSRMatrix<double> &M)
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
	  entries[i].clear();
  }

  return 0;
}

void WaveFunctionToVector(const WaveFunctionType * const psi,std::vector<double> &vecOut,const WaveFunctionType * const indexWavefunction)
{
	vecOut = std::vector<double>(indexWavefunction->N);
#pragma omp parallel
{
	uint8_t *Det = (uint8_t*)malloc(psi->NWord * sizeof(uint8_t));
	if (!Det)
	{
		printf("malloc failed for Det in WaveFunctionToVector\n");
		fflush(stdout);
	}
	#pragma omp for
	for (uint32_t i = 1; i <= psi->N; i++)
	{
		uint32_t ii = i >> NBitsBlok;
		uint32_t ik = i & BlokMask;
		uint32_t ij = psi->NWord * ik;
		uint32_t im = i % (1 << NBitsBlok);

		//uint8_t *Det = (uint8_t*)malloc(psi->NWord * sizeof(uint8_t));
		if (Det == NULL) AllocFailed();
		for (uint32_t k = 0; k < psi->NWord; k++)
		{
			Det[k] = psi->Determinant[ii][ij + k];
		}
		uint32_t m = WaveFunctionIndexDet(indexWavefunction, Det);
		if (m!=0)
		{
			vecOut[m-1] = psi->ValueRe[ii][im]; //TODO which element?
			//printf("vecOut[%d] = %.4e\n",m-1,psi->ValueRe[ii][im]);
		}
		else
		{
			//printf("Warning in WaveFunctionToVector: Determinant not found in indexWavefunction; ignoring determinant\n");
		}

	}
	free(Det);
} //end of parallel region
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


struct SortPair
{
	double val;
	uint32_t idx;
};

int CompSortPair(const void* a,const void* b)
{
	return ((SortPair*) a)->val > ((SortPair*) b)->val ? -1:1;
}
int CompSortPair2(const void* a,const void* b)
{
	return ((SortPair*) a)->val > ((SortPair*) b)->val ? 1:-1;
}
int CompSortPairIdx(const void* a,const void* b)
{
	return ((SortPair*) a)->idx > ((SortPair*) b)->idx ? -1:1;
}

void LoewdinSVDOrthonormalize(LinAlg::Matrix<double> &A, double epsSV, uint32_t &nSV)
{
	double *U, *S, *VT, *work;
	//int *iwork;
	char const *JOBU ="S";  /* in SVD A = U S V^T return the first min(Ni,Nj) collumns of U in matrix A, which is over written */
	char const *JOBVT="A";  /* in SVD A = U S V^T return full array VT */

	uint32_t Mrow = A.M;
	uint32_t Ncol = A.N;

	int info;
	uint32_t lwork = Max(5,Mrow)*Max(5,Ncol);
	nSV  = Min(Mrow,Ncol);
	U    = (double*) malloc(Mrow*Ncol       * sizeof(double));
	S 	 = (double*) malloc(nSV             * sizeof(double));
	VT 	 = (double*) malloc(Ncol * Ncol     * sizeof(double));
	work = (double*) malloc(lwork 		    * sizeof(double));
	//iwork= (int   *) malloc(8*Min(Mrow,Ncol)* sizeof(int));
	if(!U || !S || !VT || !work)
	{
		printf("malloc failed in LoewdinSVDOrthonormalize\n");
		fflush(stdout);
	}

	dgesvd(JOBU, JOBVT,&Mrow,&Ncol,A.Entries,&Mrow,S,U,&Mrow,VT,&Ncol,work,&lwork,&info);
	//dgesdd(JOBU, &Mrow, &Ncol, A.Entries, &Mrow, S, U, &Mrow, VT, &Ncol, work, &lwork, iwork, &info);
/*	printf("Singular Values:\n");
	for (uint32_t i = 0; i < nSV; i++)
	{
		printf("%.2e\n",S[i]);
	}
	printf("U :\n");
	for(uint32_t i = 0; i < Mrow; i++)
	{
		for(uint32_t j = 0; j < Ncol; j++)
		{
			printf("%.4f\t",U[j*Mrow +i]);
		}
		printf("\n");
	}
	printf("VT :\n");
	for(uint32_t i = 0; i < Ncol; i++)
	{
		for(uint32_t j = 0; j < Ncol; j++)
		{
			printf("%.4f\t",VT[j*Ncol +i]);
		}
		printf("\n");
	}
*/
	if (info)
	{
		printf("error in LoewdinSVDOrthonormalize : dgesdd info = %d\n",info);
	}
	else
	{
		//for (int i = 0; i < nSV; i++) printf("S[%d] = %.4E\n",i,S[i]);
		int idx = nSV-1;
		while((idx >= 0 ) && (S[idx] <= epsSV) )
		{
			nSV--;
			idx--;
		}
		if (nSV == 0)
		{
			printf("LoewdinSVDOrthonormalize : Reduce Dimension from %u to %u\n",Ncol,nSV);
		}
		else
		{

			char const *transN="N";
			double one  = 1;
			double zero = 0;

			dgemm(transN,transN,&Mrow,&Ncol,&nSV,&one,U,&Mrow,VT,&Ncol/*&nSV*/,&zero,A.Entries,&Mrow);
/*			printf("A = \n");
			A.print();
			printf("\n");
*/			if(Ncol-nSV)
			{
				printf("LoewdinSVDOrthonormalize : Reduce Dimension from %u to %u\n",Ncol,nSV);
				SortPair lostNorm[Ncol];
				for (uint32_t i = 0; i < Ncol; i++)
				{
					lostNorm[i].val = 0;
					lostNorm[i].idx = i;
					for (uint32_t j = nSV; j < Ncol; j++)
					{
						lostNorm[i].val += VT[i*Ncol + j]*VT[i*Ncol + j];
					}
				}
				qsort(lostNorm,Ncol,sizeof(SortPair),CompSortPair);
/*				for (uint32_t i = 0; i < Ncol; i++)
				{
					lostNorm[i].val = 0;
					lostNorm[i].idx = i;
					for (uint32_t j = 0; j < Mrow; j++)
					{
						lostNorm[i].val += A.Entries[i*Mrow + j]*A.Entries[i*Mrow + j];
					}
				}
				qsort(lostNorm,Ncol,sizeof(SortPair),CompSortPair2);
*/				qsort(lostNorm,Ncol-nSV,sizeof(SortPair),CompSortPairIdx);
				//printf("Ncol-nSV = %d\n",Ncol-nSV);
				//for (int i = 0; i < Ncol; i++) printf("lostNorm[%d].val = %.16E\t lostNorm[%d].idx = %d\n",i,lostNorm[i].val,i,lostNorm[i].idx);
				for (uint32_t i = 0; i < Ncol-nSV; i++)
				{
					//printf("%.4f\t%u\n",lostNorm[i].val,lostNorm[i].idx);
					A.RemoveCol(lostNorm[i].idx);
				}
/*				printf("A = \n");
				A.print();
				printf("\n");
*/				LoewdinSVDOrthonormalize(A,epsSV,nSV);
			}
		}

	}
	free(U);
	free(S);
	free(VT);
	free(work);
	//free(iwork);
}

inline void BlockSkalarProduct(const LinAlg::Matrix<double> &A,const LinAlg::Matrix<double> &B, LinAlg::Matrix<double> &C)
{
	const char TransA = 'T';
	const char TransB = 'N';
	double one  = 1;
	double zero = 0;

	uint32_t M = A.N;	// number of rows of A^T
	uint32_t N = B.N;	// number of cols of B
	uint32_t K = B.M;	// number of rows of B

	//C = LinAlg::Matrix<double>(M,N);

	dgemm(&TransA,&TransB,&M,&N,&K,&one,A.Entries,&K,B.Entries,&K,&zero,C.Entries,&M);
}
// A = A - B.C where A is m x n, B is m x k and C is k x n matrix
inline void SubtractMultipliedMatrix(const LinAlg::Matrix<double> &A,const LinAlg::Matrix<double> &B,const LinAlg::Matrix<double> &C)
{
	const char TransN = 'N';
	double minusone  = -1;
	double one = 1;

	dgemm(&TransN,&TransN,&B.M,&C.N,&C.M,&minusone,B.Entries,&B.M,C.Entries,&C.M,&one,A.Entries,&B.M);
}

void BlockLanczos(LinAlg::LinearMap<double, std::vector<double>, std::vector<double> > &LinMap,
				  LinAlg::Matrix<double> &startblock, uint32_t &N, double epsSV,
				  std::vector<LinAlg::Matrix<double> > &A,
				  std::vector<LinAlg::Matrix<double> > &B,
				  std::vector<LinAlg::Matrix<double> > &Q)
{
	uint32_t nSV = startblock.N;

	A = std::vector<LinAlg::Matrix<double> >(N  );
	B = std::vector<LinAlg::Matrix<double> >(N  );
	Q = std::vector<LinAlg::Matrix<double> >(N+1);

	//Q[0] = LinAlg::Matrix<double>(startblock);
	Q[0].copy(startblock);
/*	printf("Q[0] = \n");
	Q[0].print();
	printf("\n");
*/	LoewdinSVDOrthonormalize(Q[0],epsSV,nSV); 					// SVD Orthonormalize Q[0]
	B[0].reserve(Q[0].N,startblock.N);
	BlockSkalarProduct(Q[0],startblock,B[0]);					// B[0] = Q[0]^T.startblock
	//BlockSkalarProduct(startblock,Q[0],B[0]);					// B[0] = startblock^T.Q[0] since B[0] will be copied to B[0]^T

	Q[1].reserve(Q[0].M,Q[0].N);
	LinMap.BlockEvaluate(Q[0],Q[1]);							// Q[1] = LinMap(Q[0])
	A[0].reserve(Q[0].N,Q[1].N);
	BlockSkalarProduct(Q[0],Q[1],A[0]);							// A[0] = Q[0]^T.Q[1]
	SubtractMultipliedMatrix(Q[1],Q[0],A[0]);					// Q[1] = Q[1] - Q[0].A[0]
	LoewdinSVDOrthonormalize(Q[1],epsSV,nSV); 					// SVD Orthonormalize Q[1]
	//printf("nSV = %d\n",nSV);
	for (uint32_t i = 1; i < N; i++)
	{
		Q[i+1].reserve(Q[i].M,Q[i].N);
		LinMap.BlockEvaluate(Q[i],Q[i+1]);						// Q[i+1] = LinMap(Q[i])
		A[i].reserve(Q[i].N,Q[i].N);
		BlockSkalarProduct(Q[i],Q[i+1],A[i]);					// A[i]   = Q[i]^T.Q[i+1]
		B[i].reserve(Q[i-1].N,Q[i+1].N);
		BlockSkalarProduct(Q[i-1],Q[i+1],B[i]);					// B[i]   = Q[i-1]^T.Q[i+1]
		SubtractMultipliedMatrix(Q[i+1],Q[i],A[i]);				// Q[i+1] = Q[i+1] - Q[i].A[i] - Q[i-1].B[i]
		SubtractMultipliedMatrix(Q[i+1],Q[i-1],B[i]);
		LoewdinSVDOrthonormalize(Q[i+1],epsSV,nSV); 			// SVD Orthonormalize Q[i+1]
		if (nSV == 0)
		{
			N = i+1;
			printf("Complete deflation after %u iterations\n",N);
			return;
		}

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

/*
void LoewdinSVDOrthonormalize(std::vector<double> &M, uint32_t Mrow, uint32_t Ncol, double epsSV, uint32_t &nSV)
{
	double *U, *S, *VT, *work;
	char const *JOBU ="S";
	char const *JOBVT="A";

	int info;
	uint32_t lwork = Max(5,Mrow)*Max(5,Ncol);
	nSV  = Min(Mrow,Ncol);
	U    = (double*) malloc(Mrow*Ncol      * sizeof(double));
	S 	 = (double*) malloc(nSV            * sizeof(double));
	VT 	 = (double*) malloc(Ncol * Ncol    * sizeof(double));
	work = (double*) malloc(lwork 		   * sizeof(double));

	dgesvd(JOBU, JOBVT,&Mrow,&Ncol,&(M[0]),&Mrow,S,U,&Mrow,VT,&Ncol,work,&lwork,&info);
	printf("Singular Values:\n");
	for (uint32_t i = 0; i < nSV; i++)
	{
		printf("%f\n",S[i]);
	}
	printf("U :\n");
	for(uint32_t i = 0; i < Mrow; i++)
	{
		for(uint32_t j = 0; j < Ncol; j++)
		{
			printf("%.4f\t",U[j*Mrow +i]);
		}
		printf("\n");
	}
	printf("VT :\n");
	for(uint32_t i = 0; i < Ncol; i++)
	{
		for(uint32_t j = 0; j < Ncol; j++)
		{
			printf("%.4f\t",VT[j*Ncol +i]);
		}
		printf("\n");
	}

	if (info)
	{
		printf("error in LoewdinSVDOrthonormalize : dgesvd info = %d\n",info);
	}
	else
	{
		uint32_t idx = nSV-1;
		while(S[idx] <= epsSV )
		{
			nSV--;
			idx--;
		}
		if(Ncol-nSV)
		{
			printf("LoewdinSVDOrthonormalize : Reduce Dimension from %u to %u\n",Ncol,nSV);
			printf("Resize M: %u -> ",M.size());
			M.resize(nSV*Mrow);
			printf("%u\n",M.size());
		}

		char const *transAN="N";
		char const *transBN="N";
		double one  = 1;
		double zero = 0;

		dgemm(transAN,transBN,&Mrow,&nSV,&Ncol,&one,U,&Mrow,VT,&Ncol,&zero,&(M[0]),&Mrow);
		free(U);
		free(S);
		free(VT);
		free(work);
	}


}
*/

/*
inline void BlockSkalarProduct(std::vector<double> &A, std::vector<double> &B, std::vector<double> &AtB, uint32_t m, uint32_t n)
{
	const char TransA = 'T';
	const char TransB = 'N';
	double one  = 1;
	double zero = 0;

	dgemm(&TransA,&TransB,&n,&n,&m,&one,&(A[0]),&m,&(B[0]),&m,&zero,&(AtB[0]),&n);
}
 */

/*
inline void SubtractMultipliedMatrix(std::vector<double> &A,std::vector<double> &B,std::vector<double> &C, uint32_t m, uint32_t n, uint32_t k)
{
	const char TransN = 'N';
	double minusone  = -1;
	double one = 1;

	dgemm(&TransN,&TransN,&m,&n,&k,&minusone,&(B[0]),&m,&(C[0]),&k,&one,&(A[0]),&m);
}
 */

/*
void BlockLanczos(LinAlg::LinearMap<double, std::vector<double>, std::vector<double> > &LinMap,
				  std::vector<double> &startblock, const uint32_t N, const uint32_t blocksize, double epsSV,
				  std::vector<std::vector<double> > &A,
				  std::vector<std::vector<double> > &B,
				  std::vector<std::vector<double> > &Q)
{
	uint32_t Ncol = blocksize;
	uint32_t Mrow = startblock.size() / blocksize;
	uint32_t nSV = Ncol;

	std::vector<double> a = std::vector<double>(Ncol*Ncol);
	std::vector<double> b = std::vector<double>(Ncol*Ncol);

	A.reserve(N  );
	B.reserve(N  );
	Q.reserve(N+1);

	Q.push_back(startblock);
	LoewdinSVDOrthonormalize(Q[0],Mrow,Ncol,epsSV,nSV); 		// SVD Orthonormalize Q[0]
	BlockSkalarProduct(Q[0],startblock,b,Mrow,Ncol);			// B[0] = Q[0]^T.startblock
	B.push_back(b);

	printf("Q[%d]\t size = %d\tMrow=%d\n",0,Q[0].size(),Mrow);
			for (int j = 0; j < Mrow; j++)
			{
				printf("%f\t%f\n",Q[0][j],Q[0][j+(nSV-1)*Mrow]);
			}
			printf("\n");

	Q.push_back(std::vector<double>());
		LinMap.BlockEvaluate(Q[0],Q[1],nSV);					// Q[1] = LinMap(Q[0])
	BlockSkalarProduct(Q[0],Q[1],a,Mrow,nSV);					// A[0] = Q[0]^T.Q[1]
	A.push_back(a);
	SubtractMultipliedMatrix(Q[1],Q[0],a,Mrow,nSV,nSV);			// Q[1] = Q[1] - Q[0].A[0]
	LoewdinSVDOrthonormalize(Q[1],Mrow,nSV,epsSV,nSV); 			// SVD Orthonormalize Q[1]

	uint32_t nOld = nSV;
	for (uint32_t i = 1; i < N; i++)
	{
		Q.push_back(std::vector<double>());
		LinMap.BlockEvaluate(Q[i],Q[i+1],nSV);					// Q[i+1] = LinMap(Q[i])
		BlockSkalarProduct(Q[i],Q[i+1],a,Mrow,nSV);				// A[i]   = Q[i]^T.Q[i+1]
		A.push_back(a);
		BlockSkalarProduct(Q[i-1],Q[i+1],b,Mrow,nSV);			// B[i]   = Q[i-1]^T.Q[i+1]
		B.push_back(b);
		SubtractMultipliedMatrix(Q[i+1],Q[i],a,Mrow,nSV,nSV);	// Q[i+1] = Q[i+1] - Q[i].A[i] - Q[i-1].B[i]
		SubtractMultipliedMatrix(Q[i+1],Q[i-1],b,Mrow,nSV,nSV);
		LoewdinSVDOrthonormalize(Q[i+1],Mrow,nSV,epsSV,nSV); 	// SVD Orthonormalize Q[i+1]
		a.resize(nSV*nSV);
		b.resize(nSV*nOld);
		nOld = nSV;

		printf("Q[%d]\t size = %d\tMrow=%d\n",i,Q[i].size(),Mrow);
		for (int j = 0; j < Mrow; j++)
		{
			printf("%f\t%f\n",Q[i][j],Q[i][j+(nSV-1)*Mrow]);
		}
		printf("\n");

	}


}
 */
