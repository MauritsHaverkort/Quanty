/*
 * BreitInteraction.cpp
 *
 *  Created on: Sep 13, 2019
 *      Author: brass
 */

void CreateKernel(const std::vector<double> &r, const std::vector<double> &w, int k, std::vector<double> &ker)
{
	uint32_t n = (uint32_t) r.size();
	ker = std::vector<double>(n*n);
	for (uint32_t i = 1; i < n; i++) // not starting at 0 since r[0] = 0 and the wave-functions vanish
	{
		for (uint32_t j = i; j < n; j++)
		{
			double rr = pow(Min(r[i],r[j]),k-1) / pow(Max(r[i],r[j]),k);
			ker[i*n +j] = w[i] * rr * w[j];
			ker[j*n +i] = w[i] * rr * w[j];
		}
	}
}

double GeneralizedSlaterIntegral(const std::vector<double> &ker, const std::vector<double> &X, const std::vector<double> &Y)
{
	double val = 0;
	uint32_t n = (uint32_t) X.size();
	for (uint32_t i = 0; i < n; i++)
	{
		for (uint32_t j = 0; j < n; j++)
		{
			val += X[i] * ker[i*n+j] * Y[j];
		}
	}
	return val;
}

void CreateIntegrands(const std::vector<double> &g_i,const std::vector<double> &f_i,
		const std::vector<double> &g_k,const std::vector<double> &f_k,uint32_t L, vector<double> &V)
{
	uint32_t d = (uint32_t) g_i.size();
	V = vector<double>(d);
	for (uint32_t ir = 0; ir < d; ir++)
	{
		V[ir] = g_i[ir]*f_k[ir] + f_i[ir]*g_k[ir];
	}
}

void CreateIntegrands(const std::vector<double> &g_i,const std::vector<double> &f_i,
		const std::vector<double> &g_j,const std::vector<double> &f_j,
		const std::vector<double> &g_k,const std::vector<double> &f_k,
		const std::vector<double> &g_l,const std::vector<double> &f_l,
		int k_i, int k_j, int k_k, int k_l, uint32_t L,
		vector<double> &Q1, vector<double> &P1,	vector<double> &Q2, vector<double> &P2)
{
	uint32_t d = (uint32_t) g_i.size();
	int dk1 = k_k - k_i;
	int dk2 = k_l - k_j;
	double U, V;
	Q1 = vector<double>(d);
	P1 = vector<double>(d);
	Q2 = vector<double>(d);
	P2 = vector<double>(d);
	for (uint32_t ir = 0; ir < d; ir++)
	{
		U      = g_i[ir]*f_k[ir] - f_i[ir]*g_k[ir];
		Q1[ir] = -U;
		P1[ir] =  U;

		U      = g_j[ir]*f_l[ir] - f_j[ir]*g_l[ir];
		Q2[ir] = -U;
		P2[ir] =  U;

	}
	if (dk1 != 0)
	{
		for (uint32_t ir = 0; ir < d; ir++)
		{
			V = g_i[ir]*f_k[ir] + f_i[ir]*g_k[ir];
			Q1[ir] += (double)dk1/((double)L+1.) * V;
			P1[ir] += (double)dk1/(double)L * V;
		}
	}
	if (dk2 != 0)
	{
		for (uint32_t ir = 0; ir < d; ir++)
		{
			V = g_j[ir]*f_l[ir] + f_j[ir]*g_l[ir];
			Q2[ir] += (double)dk2/((double)L+1.) * V;
			P2[ir] += (double)dk2/(double)L * V;
		}
	}
}

double CReducedJBasis(int k_a, int k_b, int L)
{
	int l_a = k_a > 0 ? k_a : -k_a -1;
	int l_b = k_b > 0 ? k_b : -k_b -1;
	if ((l_a + l_b + L)%2)
	{
		return 0;
	}
	else
	{
		int j_a2 = k_a > 0 ? 2*l_a -1 : 2*l_a +1;
		int j_b2 = k_b > 0 ? 2*l_b -1 : 2*l_b +1;

		if ((j_a2 +1)/2 %2)
			return -sqrt((j_a2+1)*(j_b2+1)) * ThreeJSymbolHalf(j_a2,-1,j_b2,1,2*L,0);
		else
			return  sqrt((j_a2+1)*(j_b2+1)) * ThreeJSymbolHalf(j_a2,-1,j_b2,1,2*L,0);
	}
}

void CalcBreitInteraction(const std::vector<Interpolation::InterpolatingFunction*> &G,const std::vector<Interpolation::InterpolatingFunction*> &F, uint32_t N)
{
	double M1, M2, N1, O1, O2, O31, O32, O41, O42;
	vector<double> x,w, Q1, P1, Q2, P2, V1,V2;
	vector<vector<double> > g,f, ker;
	uint32_t Nfuns = (uint32_t) G.size();
	BeginTiming("w", &GlobalTime);
	Quad::GaussianQuadratureWeights(N,G[0]->x,x,w);
	EndTiming("w", &GlobalTime);

	uint32_t i,j,k,l;
	int k_k, k_i, k_j, k_l, Lmin,Lmax;
	//printf("%f\t%f\t%f\n",CReducedJBasis(-1,1,1),CReducedJBasis(-1,2,1),CReducedJBasis(-1,2,2));
	i = 0;
	j = 0;
	k = 1;
	l = 1;
	k_i = -1;
	k_j = -1;
	k_k = -2;
	k_l = -2; //TODO: remove these 8 lines

	Lmin = 0;
	Lmax = 3; //TODO: remove these 2 lines
	ker = vector<vector<double> >(Lmax+1);
	BeginTiming("ker", &GlobalTime);
	for (int L = Max(Lmin,1); L <= Lmax; L++)
	{
		CreateKernel(x,w,L,ker[L]);
	}
	EndTiming("ker", &GlobalTime);
	BeginTiming("funs", &GlobalTime);
	g = vector<vector<double> >(Nfuns);
	f = vector<vector<double> >(Nfuns);
	for (uint32_t ifun = 0; ifun < Nfuns; ifun++)
	{
		g[ifun] = (*G[ifun])(x);
		f[ifun] = (*F[ifun])(x);
	}
	EndTiming("funs", &GlobalTime);

	Lmin = 1;
	Lmax = 1; //TODO: remove these 2 lines
	BeginTiming("Int", &GlobalTime);
	for (int L = Lmin; L <= Lmax; L += 2)
	{
		CreateIntegrands(g[i],f[i],g[j],f[j],g[k],f[k],g[l],f[l],k_i,k_j,k_k,k_l,L,Q1,P1,Q2,P2);

		M1 = GeneralizedSlaterIntegral(ker[L+2],Q1,Q2);
		O1 = M1;
		M1 *= ((double)L+1.)/(2.*(double)L-1.); //TODO: check if this is really division by double
		//printf("M1 = %.16f\n",M1);
		O1 *= ((double)L+1.)*((double)L+1.)/(2.*(double)L+1.)/(2.*(double)L+3.);
		if ( L == 0)
		{
			M2  = 0;
			O2  = 0;
			O31 = 0;
			O32 = 0;
			O41 = 0;
			O42 = 0;
		}
		else //TODO: check these; above integrals are tested
		{
			M2 = GeneralizedSlaterIntegral(ker[L],P1,P2);
			O2 = M2;
			M2 *= (double)L/(2.*(double)L-1.);
			O2 *= (double)L*L/(2.*(double)L-1.)/(2.*(double)L+1.);
			O31 = GeneralizedSlaterIntegral(ker[L],Q1,P2);
			O31 *= (double)L*((double)L+1.)/2./(2.*(double)L+1);
			O32 = -GeneralizedSlaterIntegral(ker[L+2],Q1,P2);
			O32 *= (double)L*((double)L+1.)/2./(2.*(double)L+1);
			//printf("O32 = %.16f\n",O32);
			O41 = GeneralizedSlaterIntegral(ker[L],P1,Q2);
			O41 *= (double)L*((double)L+1.)/2./(2.*(double)L+1);
			O42 = -GeneralizedSlaterIntegral(ker[L+2],P1,Q2);
			O42 *= (double)L*((double)L+1.)/2./(2.*(double)L+1);
			//printf("O42 = %.16f\n",O42);
		}
		double MO = M1 + M2 - (O1 + O2 + O31 + O32 + O41 + O42);
		if (L%2)
			MO *= -1;
		MO *= CReducedJBasis(k_i,k_k,L) * CReducedJBasis(k_j,k_l,L);
		printf("M  = %.16f\n",CReducedJBasis(k_i,k_k,L) * CReducedJBasis(k_j,k_l,L)*(M1 + M2)*(L%2?-1:1));
		printf("O  = %.16f\n",-CReducedJBasis(k_i,k_k,L) * CReducedJBasis(k_j,k_l,L)*(O1 + O2 + O31 + O32 + O41 + O42)*(L%2?-1:1));
		printf("MO = %.16f\n",MO);
	}
	Lmin = 2;
	Lmax = 2; //TODO: remove these 2 lines
	for (int L = Lmin; L <= Lmax; L += 2)
	{
		CreateIntegrands(g[i],f[i],g[k],f[k],L,V1);
		CreateIntegrands(g[j],f[j],g[l],f[l],L,V2);
		N1 = GeneralizedSlaterIntegral(ker[L+1],V1,V2); //TODO: immer anderes L als der rest
		N1 *= (double)(k_i+k_k)*(k_j+k_l)/(double)L/((double)L+1.) * CReducedJBasis(-k_i,k_k,L) * CReducedJBasis(-k_j,k_l,L) * (L%2 ? 1 : -1);
		printf("N1 = %.16f\n",N1);
	}
	EndTiming("Int", &GlobalTime);

	//printf("Fk = %.16f\n",GeneralizedSlaterIntegral(ker,Q,P));

}


