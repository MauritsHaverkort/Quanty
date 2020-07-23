/*
 * Outdated_Functions.cpp
 *
 *  Created on: Nov 7, 2019
 *      Author: brass
 */

double TestGSEnergy(const vector<vector<uint32_t> > &Nocc, const vector<vector<double> > &eigVal, const vector<vector<double> > &eigVec, const vector<vector<double> > &H, uint32_t n)
{
	double e0 = 0;
	uint32_t dim = (uint32_t) eigVal[0].size();
	for (uint32_t  l = 0; l < (uint32_t) Nocc.size(); l++)
	{
		for (uint32_t i = 0; i < (uint32_t) Nocc[l].size(); i++)
		{

			for (uint32_t mu = 0; mu < dim; mu++)
			{
				for (uint32_t nu = 0; nu < dim; nu++)
					e0 += (2*l+1) * H[l][mu*dim + nu]*eigVec[l][i*dim +mu]*eigVec[l][i*dim +nu];
			}
		}
	}
	return e0;
}

void CreateHartreeMatrix(const uint32_t la, const vector<double> &rk1, const vector<double> &rk2, const vector<double> &rkD,
					  uint32_t n, const vector<vector<double> > &eigVec, uint32_t dim, const vector<vector<uint32_t> > &Nocc, vector<double> &H, vector<double> &HE)
{
	H = vector<double>(dim*dim);
	HE = vector<double>(dim*dim);
	for (uint32_t nu = 1; nu <= dim; nu++)
	{
		uint32_t mu1 = nu+n < dim ? nu+n : dim;
		for (uint32_t mu = nu; mu <= mu1; mu++)
		{
			for (uint32_t alpha = 1; alpha <= dim; alpha++)
			{
				uint32_t beta0 = alpha > n ? alpha-n : 1;
				uint32_t beta1 = alpha+n < dim ? alpha+n : dim;
				for (uint32_t beta = beta0; beta <= beta1; beta++) //TODO: alpha <= beta
				{
					for (uint32_t lb = 0; lb < (uint32_t) Nocc.size(); lb++)
					{

						double Fk = AssembleSlaterInts(nu,alpha,mu,beta,rk1,rk2,rkD,n);
						for (uint32_t  j = 0; j < (uint32_t) Nocc[lb].size(); j++)
						{
							H[(nu-1) * dim + (mu-1)]  += Fk* eigVec[lb][j * dim+(alpha-1)] * eigVec[lb][j * dim+(beta-1)]*(2*lb+1);
							//TODO: Lapack has Col Major and stores eigenvectors collumnwise, hence the j-th ev should be at [j*dim + i]
							// but is that really true or does Lapack do some strange stuff?
						}

					}
				}
			}
			H[(nu-1) + dim * (mu-1)] = H[(nu-1) * dim + (mu-1)];
		}
	}
}

void CreateFockMatrix(const uint32_t la, const vector<vector<double> > &rk1, const vector<vector<double> > &rk2, const vector<vector<double> > &rkD,
					  uint32_t n, const vector<vector<double> > &eigVec, uint32_t dim, const vector<vector<uint32_t> > &Nocc, vector<double> &F, vector<double> &FE)
{
	F = vector<double>(dim*dim);
	FE = vector<double>(dim*dim);
	for (uint32_t nu = 1; nu <= dim; nu++)
	{
		for (uint32_t mu = nu; mu <= dim; mu++)
		{
			uint32_t alpha0 = mu > n ? mu-n : 1;
			uint32_t alpha1 = mu+n < dim ? mu+n : dim;
			uint32_t beta0 = nu > n ? nu-n : 1;
			uint32_t beta1 = nu+n < dim ? nu+n : dim;
			for (uint32_t alpha = alpha0; alpha <= alpha1; alpha++)
			{
				for (uint32_t beta = beta0; beta <= beta1; beta++)
				{
					for (uint32_t lb = 0; lb < (uint32_t) Nocc.size(); lb++)
					{
						uint32_t k0 = Abs((int)la - (int)lb);
						uint32_t k1 = la + lb;
						if ((k1 +k0) % 2) k0++;
						double Fk = 0;
						for (uint32_t k = k0; k <= k1; k +=2)
						{
							double cred = ThreeJSymbol(la,0,lb,0,k,0);
							//if (k ==1) printf("ThreeJSymbol(%d,0,%d,0,%d,0)^2 = %f\n",la,lb,k,cred*cred);
							Fk += AssembleSlaterInts(nu,alpha,beta,mu,rk1[k],rk2[k],rkD[k],n)*cred*cred; // *(2*k+1)
						}
						for (uint32_t  j = 0; j < (uint32_t) Nocc[lb].size(); j++)
						{
							F[(nu-1) * dim + (mu-1)]  += Fk* eigVec[lb][j * dim+(alpha-1)] * eigVec[lb][j * dim+(beta-1)]*(2*lb+1);
							//TODO: Lapack has Col Major and stores eigenvectors collumnwise, hence the j-th ev should be at [j*dim + i]
							// but is that really true or does Lapack do some strange stuff?
						}

					}
				}
			}
			F[(nu-1) + dim * (mu-1)]  = F[(nu-1) * dim + (mu-1)];
		}
	}
}

void CreateHartreeFockOperator(const vector<double> &P, const vector<double> &D, uint32_t dim, uint32_t lMax, vector<vector<double> > &H)
{
	uint32_t nl = lMax + 1;
	uint32_t superdim = dim*dim*nl;
	H = vector<vector<double> >(nl);
	for (uint32_t la = 0; la <nl; la++)
	{
		H[la] = vector<double>(dim*dim);
		for (uint32_t mu = 0; mu < dim; mu++)
		{
			for (uint32_t nu = mu; nu < dim; nu++)
			{
				for (uint32_t lb = 0; lb < nl; lb++)
				{
					for (uint32_t alpha = 0; alpha < dim; alpha++)
					{
						for (uint32_t beta = 0; beta < dim; beta++)
						{
							H[la][mu*dim + nu] += D[alpha*dim*nl + beta*nl + lb] * P[nu*superdim*dim*nl + mu*superdim*nl + la*superdim + alpha*dim*nl + beta*nl +lb];

						}
					}
				}
				H[la][nu*dim + mu] = H[la][mu*dim + nu];
			}
		}
	}

}
