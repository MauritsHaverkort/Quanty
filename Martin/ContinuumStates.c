/*
 * ContinuumStates.c
 *
 *  Created on: May 17, 2018
 *      Author: brass
 */

#ifndef speedoflight
#define speedoflight 137.0
#endif

double Gamma(double const x){ // Martin Brass copied this from Martin Sundermann and Martin Rotter edit by MS
	int i, gg=7, np=9;
	double const p[] = {0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7};
	double y,a=p[0];
  if(x<0.5){
		for(i=1;i<np;i++) a+=p[i]/(-x+(double)i);
		y=Pi/(sin(Pi*x)*sqrt(2.*Pi)*pow(-x+(double)gg+0.5,-x+0.5)*exp(x-(double)gg-0.5)*a);
  }else{
		for(i=1;i<np;i++) a+=p[i]/(x+(double)i-1.0);
		y=sqrt(2.*Pi)*pow(x+(double)gg-0.5,x-0.5)*exp(-x-(double)gg+0.5)*a;
	}
  return y;
}

double Binomial(double x, double y)
{
	return Gamma(x+1)/Gamma(y+1)/Gamma(x-y+1);
}

double EnergyDO(int n, int k, int l, double w)
{
	return speedoflight * sqrt(speedoflight*speedoflight + 2 * w *(1 + l + 2*n +k))-speedoflight*speedoflight;
}

double NormDO(int n, int k, int l, double w)
{
	double b1,b2, g1, g2, coeff, res;
	if ( k > 0 )
	{
		b1 = Binomial(n + l + 1./2.,n);
		b2 = Binomial(n + l - 1./2.,n);
		g1 =  Gamma(n + l + 3./2.) ;
		g2 =  Gamma(n + l + 1./2.) ;
		coeff = (1+2*k)/(speedoflight + sqrt(speedoflight*speedoflight + 2*(1 + l + 2*n +k)*w));
		res = sqrt(g1/(2 * sqrt(pow(w,2*l+3)) * Factorial(n) * b1*b1 ) + g2/(2 * sqrt(pow(w,2*l+1)) * Factorial(n) * b2*b2 )*coeff*coeff );
		return res;
	}
	else
	{
		b1 = Binomial(n + l + 1./2.,n);
		b2 = Binomial(n-1 + l + 3./2.,n-1);
		g1 =  Gamma(n + l + 3./2.) ;
		g2 =  g1 ;
		coeff = (-speedoflight + sqrt(speedoflight*speedoflight + 2*(1 + l + 2*n +k)*w)) / (-1+2*k);
		res = sqrt(g1/(2 * sqrt(pow(w,2*l+3)) * Factorial(n) * b1*b1 ) + g2/(2 * sqrt(pow(w,2*l+5)) * Factorial(n-1) * b2*b2)*coeff*coeff );
		return res;

	}
}

double C1DO(int n, int k, int l, double w)
{
	return 1./(Binomial(n+l+1./2.,n) * sqrt(pow(w,l+1)) * NormDO(n,k,l,w));
}

double C2DO(int n, int k, int l, double w)
{
	double bin, coeff;
	if ( k > 0 )
	{
		coeff = (1+2*k)/(speedoflight + sqrt(speedoflight*speedoflight + 2*(1 + l + 2*n +k)*w))/ sqrt(pow(w,l));
		bin = Binomial(n + l - 1./2.,n);
	}
	else
	{
		coeff = (-speedoflight + sqrt(speedoflight*speedoflight + 2*(1 + l + 2*n +k)*w)) / (-1+2*k)/ sqrt(pow(w,l+2));
		bin = Binomial(n-1 + l + 3./2.,n-1);

	}
	return coeff/bin/NormDO(n,k,l,w);
}

double EKinDO(int n, int m, int k, int l, double w)
{
	switch (m-n)
	{
		case 0:
			return EnergyDO(n,k,l,w) + speedoflight * C1DO(n,k,l,w)*C2DO(n,k,l,w)*Gamma(n+l+3./2.)/Factorial(n-1);
			break;
		case 1:
			return -speedoflight/2. * C1DO(n,k,l,w)*C2DO(m,k,l,w)*Gamma(n+l+5./2.)/Factorial(n);
			break;
		case -1:
			return -speedoflight/2. * C1DO(m,k,l,w)*C2DO(n,k,l,w)*Gamma(m+l+5./2.)/Factorial(m);
			break;
		default:
			return 0.;
			break;
	}
}

double **CalculateContinuumBoundDiracElements(double **basisfuns, double **derivativefuns, double **dualcont, const int NcontFuns, const unsigned int *indices,
		const unsigned int Nindices, const double *x, unsigned int NGridPoints, const int kappa, const int Z)
{
	double **result = (double **) malloc(NcontFuns * sizeof(double *));

	unsigned int i, j, m;
	for ( i = 0; i < NcontFuns; i++ )
	{
		result[i] = (double *) malloc(Nindices * sizeof(double));
		for ( j = 0; j < Nindices; j++)
		{
			m = indices[j];
			result[i][j] = DiracH(basisfuns[m],dualcont[i],derivativefuns[m], x, NGridPoints, kappa, Z);
		}
	}


	return result;

}

double *KineticEnergyCB(double **basisfuns, double **derivativefuns, double **dualfuns,double **dualcont, const int NcontFuns, const unsigned int *indices,
		const unsigned int Nindices, const double *x, unsigned int NGridPoints, const int kappa,double w)
{
	double **Tcb, **Tbb, *T;
	int i,j, dim, l;
	const int Z = 0;

	if ( kappa > 0 )
	{
		l = kappa;
	}
	else
	{
		l = -kappa -1;
	}
	Tbb = CalculateDiracMatrixElements(basisfuns,derivativefuns,dualfuns,indices,Nindices,x,NGridPoints,kappa,Z);
	Tcb = CalculateContinuumBoundDiracElements(basisfuns,derivativefuns,dualcont,NcontFuns,indices,Nindices,x,NGridPoints,kappa,Z);

	dim = Nindices + NcontFuns;
	T = (double *) malloc(dim*dim*sizeof(double));

	for ( i = 0; i < Nindices; i++)
	{

		for ( j = i; j < Nindices; j++)
		{
			T[i*dim + j] = Tbb[i][j];
			T[j*dim + i] = Tbb[i][j];
		}
	}
	for ( i = Nindices; i < dim; i++)
	{
		for ( j = 0; j < Nindices; j++)
		{
			T[i*dim + j] = Tcb[i-Nindices][j];
			T[j*dim + i] = Tcb[i-Nindices][j];
		}
		for ( j = Nindices; j < dim; j++)
		{
			T[i*dim + j] = EKinDO(i-Nindices+1,j-Nindices+1,kappa,l,w);
		}
	}
	free(Tbb);
	free(Tcb);
	return T;
}

double *MetricCB(double **basisfuns, double **dualfuns,const int NBasisFuns, double **dualcont, const int NcontFuns, unsigned int NGridPoints)
{
	int dim,i,j,n,eins;
	eins = 1;
	n = 2* NGridPoints;
	dim = NBasisFuns + NcontFuns;
	double *g = (double *) malloc(dim*dim*sizeof(double));
	double y;

	for ( i = 0; i < NBasisFuns; i++)
	{
		for ( j = i; j < NBasisFuns; j++)
		{
			y = i==j ? 1. : 0.;
			g[i*dim +j] = y;
			g[j*dim +i] = y;
		}
	}
	for ( i = NBasisFuns; i < dim; i++)
	{
		for (j = 0; j < NBasisFuns; j++)
		{

			y = ddot(&n,dualcont[i-NBasisFuns],&eins,basisfuns[j],&eins);
			g[i*dim + j] = y;
			g[j*dim + i] = y;
		}
		for ( j = NBasisFuns; j < dim; j++)
		{
			y = i==j ? 1. : 0.;
			g[i*dim +j] = y;
			g[j*dim +i] = y;
		}
	}
	return g;
}

double *SubMatrix(double *M, int dim, int n1, int n2, int m1, int m2 )
{
	double *Msub = (double *) malloc((n2-n1+1)*(m2-m1+1)*sizeof(double));
	int i,j;
	for ( i = 0; i < n2-n1+1; i++)
	{
		for ( j = 0; j < m2-m1+1; j++)
		{
			Msub[i*(m2-m1+1) + j] = M[(i+n1)*dim + (j+m1)];
		}
	}
	return Msub;
}

void vAdd(const double *a,const double aleph, double* b,const int n)
{
	int i;
	for ( i = 0; i < n; i++)
	{
		b[i] += aleph * a[i];
	}
}

void svMul(const double aleph, double *a, const int n)
{
	int i;
	for(i = 0; i < n; i++)
	{
		a[i] *= aleph;
	}
}

void CalcContinuumWFs(double **basisfuns, const int NBasisfuns, double **contfuns, const int Ncont,const int NGridpoints, double *U, double *R)
{
	int i,j,dim,n, iU,jU;

	n = NGridpoints * 2;
	dim = Ncont + NBasisfuns;
	for ( i = Ncont-1; i >= 0; i--)
	{
		iU = i + NBasisfuns;
		svMul(U[iU*dim+ iU],contfuns[i],n);
		for ( j = i-1; j >= 0; j--)
		{
			jU = j + NBasisfuns;
			vAdd(contfuns[j],U[jU*dim + iU],contfuns[i],n);
		}
		for( j = 0; j < NBasisfuns; j++)
		{
			vAdd(basisfuns[j],U[j*dim + iU],contfuns[i],n);
		}
	}
}


























