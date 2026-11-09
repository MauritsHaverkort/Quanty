/*
 * ContinuumStates.c
 *
 *  Created on: May 17, 2018
 *      Author: brass
 */

#include "ContinuumStates.h"

/*
double DiracH2(double *basisfun,double *dualfun,double *derivfun,const double *x, int NGridPoints, const int kappa, const int Z)
{
	double *H11G, *H22F,*H12F, *H21G, result, r;
	//const double speedoflight = 137;
	int i, eins;
	eins = 1;

	H11G = (double *) malloc(NGridPoints * sizeof(double));
	H22F = (double *) malloc(NGridPoints * sizeof(double));
	H12F = (double *) malloc(NGridPoints * sizeof(double));
	H21G = (double *) malloc(NGridPoints * sizeof(double));

	if ( x[0] == 0. )
	{
		H11G[0] = 0.;
		H22F[0] = 0.;
		H12F[0] = 0.;
		H21G[0] = 0.;
		i = 1;
	}
	else
	{
		i = 0;
	}
	while ( i < NGridPoints)
	{
		r = x[i];
		H11G[i] = (-Z)/r * basisfun[i];
		H22F[i] = (-2*speedoflight*speedoflight-Z/r) * basisfun[NGridPoints+i];
		H12F[i] = kappa/r * basisfun[NGridPoints+i] - derivfun[NGridPoints+i];
		H21G[i] = kappa/r * basisfun[i] + derivfun[i];
		i++;
	}
	result = ddot(&NGridPoints,dualfun,&eins,H11G,&eins) + speedoflight * ddot(&NGridPoints,dualfun,&eins,H12F,&eins)
			+ddot(&NGridPoints,&(dualfun[NGridPoints]),&eins,H22F,&eins) + speedoflight * ddot(&NGridPoints,&(dualfun[NGridPoints]),&eins,H21G,&eins);
	free(H11G);
	free(H12F);
	free(H21G);
	free(H22F);

	return result;
}

double **CalculateDiracMatrixElements2(double **basisfuns, double **derivativefuns, double **dualfuns, const unsigned int *indices, const unsigned int Nindices, const double *x,
		unsigned int NGridPoints, const int kappa, const int Z)
{
	double **result = (double **) malloc(Nindices * sizeof(double *));

	unsigned int i, j, n, m;
	for ( i = 0; i < Nindices; i++ )
	{
		result[i] = (double *) malloc(Nindices * sizeof(double));
		n = indices[i];
		for ( j = i; j < Nindices; j++)
		{
			m = indices[j];
			result[i][j] = DiracH(basisfuns[m],dualfuns[n],derivativefuns[m], x, NGridPoints, kappa, Z, speedoflight);
		}
	}


	return result;

}
*/

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

//double Phi(double n, double a, double b);

double Binomial(double x, double y)
{
	double a, b;
	if ( x > 100 && y > 100)
	{
		b = 0;
		a = x-y;
		//TODO: Check a <= 0
		return Phi(y+1,a,b)/Gamma(x-y+1);
	}
	else
	{
		return Gamma(x+1)/Gamma(y+1)/Gamma(x-y+1);
	}
}

double Phi(double n, double a, double b)
{
	double c0, c1, c2;
	c0 = 1.;
	c1 = 0.5 * (a-b) * (a+b-1);
	c2 = 1./12. * Binomial(a-b,2.)*(3*(a+b-1)*(a+b-1) - a + b -1);
	return (c0*n*n  + c1*n + c2)* pow(n,a-b-2);

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
		coeff = (1+2*k)/(speedoflight + sqrt(speedoflight*speedoflight + 2*(1 + l + 2*n +k)*w));
		if ( n > 100 )
		{
			res = sqrt(Phi(n, l + 3./2.,1.)/(2 * sqrt(pow(w,2*l+3)) * b1*b1 ) + Phi(n,l + 1./2.,1.)/(2 * sqrt(pow(w,2*l+1)) * b2*b2 )*coeff*coeff );
		}
		else
		{

			g1 =  Gamma(n + l + 3./2.) ;
			g2 =  Gamma(n + l + 1./2.) ;

			res = sqrt(g1/(2 * sqrt(pow(w,2*l+3)) * Factorial(n) * b1*b1 ) + g2/(2 * sqrt(pow(w,2*l+1)) * Factorial(n) * b2*b2 )*coeff*coeff );
		}
		return res;
	}
	else
	{
		b1 = Binomial(n + l + 1./2.,n);
		b2 = Binomial(n-1 + l + 3./2.,n-1);
		coeff = (-speedoflight + sqrt(speedoflight*speedoflight + 2*(1 + l + 2*n +k)*w)) / (-1+2*k);
		if ( n > 100 )
		{
			res = sqrt(Phi(n,l + 3./2.,1)/(2 * sqrt(pow(w,2*l+3)) * b1*b1 ) + Phi(n,l + 3./2.,0)/(2 * sqrt(pow(w,2*l+5)) * b2*b2)*coeff*coeff );
		}
		else
		{
			g1 =  Gamma(n + l + 3./2.) ;
			g2 =  g1 ;

			res = sqrt(g1/(2 * sqrt(pow(w,2*l+3)) * Factorial(n) * b1*b1 ) + g2/(2 * sqrt(pow(w,2*l+5)) * Factorial(n-1) * b2*b2)*coeff*coeff );
		}
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
	double res;
	switch (m-n)
	{
		case 0:
			if (n > 100)
			{
				res = EnergyDO(n,k,l,w) + speedoflight * C1DO(n,k,l,w)*C2DO(n,k,l,w)*Phi(n,l+3./2.,0);
			}
			else
			{
				res = EnergyDO(n,k,l,w) + speedoflight * C1DO(n,k,l,w)*C2DO(n,k,l,w)*Gamma(n+l+3./2.)/Factorial(n-1);
			}
			return res;
			break;
		case 1:
			if (n > 100)
			{
				res = -speedoflight/2. * C1DO(n,k,l,w)*C2DO(m,k,l,w)*Phi(n,l+5./2.,1);
			}
			else
			{
				res = -speedoflight/2. * C1DO(n,k,l,w)*C2DO(m,k,l,w)*Gamma(n+l+5./2.)/Factorial(n);
			}
			return res;
			break;
		case -1:
			if (n > 100)
			{
				res = -speedoflight/2. * C1DO(m,k,l,w)*C2DO(n,k,l,w)*Phi(m,l+5./2.,1);
			}
			else
			{
				res = -speedoflight/2. * C1DO(m,k,l,w)*C2DO(n,k,l,w)*Gamma(m+l+5./2.)/Factorial(m);
			}
			return res;
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
	if (!result) printf("malloc fail at result in CalculateContinuumBoundDiracElemnts. size = %d\n",NcontFuns);

	unsigned int i, j, m;
	for ( i = 0; i < NcontFuns; i++ )
	{
		result[i] = (double *) malloc(Nindices * sizeof(double));
		if (!result[i]) printf("malloc fail at result[%d] in CalculateContinuumBoundDiracElemnts. size = %d\n",i,NcontFuns);
		for ( j = 0; j < Nindices; j++)
		{
			m = indices[j];
			result[i][j] = DiracH(basisfuns[m],dualcont[i],derivativefuns[m], x, NGridPoints, kappa, Z, speedoflight);
		}
	}


	return result;

}

vector<vector<double> > CalculateContinuumBoundDiracElements(vector<vector<double> > basisfuns, vector<vector<double> > derivativefuns,
		vector<vector<double> > dualcont, const int NcontFuns, const unsigned int *indices,
		const unsigned int Nindices, vector<double> x, unsigned int NGridPoints, const int kappa, const int Z)
{
	vector<vector<double> > result ;
	vector<double> help;

	unsigned int i, j, m;
	for ( i = 0; i < NcontFuns; i++ )
	{
		result.insert(result.end(),help);
		for ( j = 0; j < Nindices; j++)
		{
			m = indices[j];
			result[i].insert(result[i].end(), DiracH(basisfuns[m],dualcont[i],derivativefuns[m], x, NGridPoints, kappa, Z, speedoflight));
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
	Tbb = CalculateDiracMatrixElements(basisfuns,derivativefuns,dualfuns,indices,Nindices,x,NGridPoints,kappa,Z,speedoflight);
	Tcb = CalculateContinuumBoundDiracElements(basisfuns,derivativefuns,dualcont,NcontFuns,indices,Nindices,x,NGridPoints,kappa,Z);

	dim = Nindices + NcontFuns;
	T = (double *) malloc(dim*dim*sizeof(double));
	if(!T) printf("malloc fail at T in KineticEnergyCB. size = %d\n",dim*dim);

	for ( i = 0; i < Nindices; i++)
	{

		for ( j = i; j < Nindices; j++)
		{
			T[i*dim + j] = Tbb[i][j];
			T[j*dim + i] = Tbb[i][j];
	        //printf("T %d\t%d\n",dim*dim-1-(i*dim+j),dim*dim-1-(j*dim+i));
	        //fflush(stdout);
		}
	}
	for ( i = Nindices; i < dim; i++)
	{
		for ( j = 0; j < Nindices; j++)
		{
			T[i*dim + j] = Tcb[i-Nindices][j];
			T[j*dim + i] = Tcb[i-Nindices][j];
	        //printf("T %d\t%d\n",dim*dim-1-(i*dim+j),dim*dim-1-(j*dim+i));
	        //fflush(stdout);
		}
		for ( j = Nindices; j < dim; j++)
		{
			T[i*dim + j] = EKinDO(i-Nindices+1,j-Nindices+1,kappa,l,w);//TODO check transpose
	        //printf("T %d\n",dim*dim-1-i*dim+j);
	        //fflush(stdout);
		}
	}
	for ( i = 0; i < Nindices; i++)
	{
		free(Tbb[i]);
	}
	free(Tbb);
	for ( i = 0; i < NcontFuns; i++)
	{
		free(Tcb[i]);
	}
	free(Tcb);
	return T;
}

vector<double> KineticEnergyCB(vector<vector<double> > basisfuns, vector<vector<double> > derivativefuns, vector<vector<double> > dualfuns,
		vector<vector<double> > dualcont, const int NcontFuns, const unsigned int *indices,
		const unsigned int Nindices, vector<double> x, unsigned int NGridPoints, const int kappa,double w)
{
	vector<vector<double> >  Tbb;
	vector<vector<double> >  Tcb;
	vector<double> T;
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
	Tbb = CalculateDiracMatrixElements(basisfuns,derivativefuns,dualfuns,indices,Nindices,x,NGridPoints,kappa,Z,speedoflight);
	Tcb = CalculateContinuumBoundDiracElements(basisfuns,derivativefuns,dualcont,NcontFuns,indices,Nindices,x,NGridPoints,kappa,Z);

	dim = Nindices + NcontFuns;
	for ( i = 0; i < dim*dim; i++)
	{
		T.insert(T.end(),0.); //TODO that could be better
	}
	for ( i = 0; i < Nindices; i++)
	{

		for ( j = i; j < Nindices; j++)
		{
			T[i*dim + j] = Tbb[i][j];
			T[j*dim + i] = Tbb[i][j];
	        //printf("T %d\t%d\n",dim*dim-1-(i*dim+j),dim*dim-1-(j*dim+i));
	        //fflush(stdout);
		}
	}
	for ( i = Nindices; i < dim; i++)
	{
		for ( j = 0; j < Nindices; j++)
		{
			T[i*dim + j] = Tcb[i-Nindices][j];
			T[j*dim + i] = Tcb[i-Nindices][j];
	        //printf("T %d\t%d\n",dim*dim-1-(i*dim+j),dim*dim-1-(j*dim+i));
	        //fflush(stdout);
		}
		for ( j = Nindices; j < dim; j++)
		{
			T[i*dim + j] = EKinDO(i-Nindices+1,j-Nindices+1,kappa,l,w);//TODO check transpose
	        //printf("T %d\n",dim*dim-1-i*dim+j);
	        //fflush(stdout);
		}
	}

	return T;
}

double *MetricCB(double **basisfuns, double **dualfuns, const unsigned int *indices,
		const unsigned int Nindices, double **dualcont, const int NcontFuns, unsigned int NGridPoints,bool QRelativistic)
{
	int dim,i,j,n;
	unsigned int eins = 1;
	n = QRelativistic ? 2* NGridPoints : NGridPoints;
	dim = Nindices + NcontFuns;
	double *g = (double *) malloc(dim*dim*sizeof(double));
	if (!g) printf("malloc fail at g in MetricCB. size = %d\n",dim*dim);
	double y;

	for ( i = 0; i < Nindices; i++)
	{
		for ( j = i; j < Nindices; j++)
		{
			y = i==j ? 1. : 0.;
			g[i*dim +j] = y;
			g[j*dim +i] = y;
	        //printf("g %d\t%d\n",dim*dim-1-(i*dim+j),dim*dim-1-(j*dim+i));
	        //fflush(stdout);
		}
	}
	for ( i = Nindices; i < dim; i++)
	{
		for (j = 0; j < Nindices; j++)
		{

			y = ddot((unsigned int *) &n,dualcont[i-Nindices],&eins,basisfuns[indices[j]],&eins);
			g[i*dim + j] = y;
			g[j*dim + i] = y;
	        //printf("g %d\t%d\n",dim*dim-1-(i*dim+j),dim*dim-1-(j*dim+i));
	        //fflush(stdout);
		}
		for ( j = Nindices; j < dim; j++)
		{
			y = i==j ? 1. : 0.;
			g[i*dim +j] = y;
			g[j*dim +i] = y;
	        //printf("g %d\t%d\n",dim*dim-1-(i*dim+j),dim*dim-1-(j*dim+i));
	        //fflush(stdout);
		}
	}
	return g;
}

vector<double> MetricCB(vector<vector<double> > basisfuns, vector<vector<double> > dualfuns, const unsigned int *indices,
		const unsigned int Nindices, vector<vector<double> > dualcont, const int NcontFuns, unsigned int NGridPoints,bool QRelativistic)
{
	int dim,i,j,n;
	unsigned int eins = 1;
	n = QRelativistic ? 2* NGridPoints : NGridPoints;
	dim = Nindices + NcontFuns;
	vector<double> g;
	for ( i = 0; i < dim*dim; i++)
	{
		g.insert(g.end(),0.);
	}
	double y;

	for ( i = 0; i < Nindices; i++)
	{
		for ( j = i; j < Nindices; j++)
		{
			y = i==j ? 1. : 0.;
			g[i*dim +j] = y;
			g[j*dim +i] = y;
	        //printf("g %d\t%d\n",dim*dim-1-(i*dim+j),dim*dim-1-(j*dim+i));
	        //fflush(stdout);
		}
	}
	for ( i = Nindices; i < dim; i++)
	{
		for (j = 0; j < Nindices; j++)
		{

			y = ddot((unsigned int *) &n,&(dualcont[i-Nindices][0]),&eins,&(basisfuns[indices[j]][0]),&eins);
			g[i*dim + j] = y;
			g[j*dim + i] = y;
	        //printf("g %d\t%d\n",dim*dim-1-(i*dim+j),dim*dim-1-(j*dim+i));
	        //fflush(stdout);
		}
		for ( j = Nindices; j < dim; j++)
		{
			y = i==j ? 1. : 0.;
			g[i*dim +j] = y;
			g[j*dim +i] = y;
	        //printf("g %d\t%d\n",dim*dim-1-(i*dim+j),dim*dim-1-(j*dim+i));
	        //fflush(stdout);
		}
	}
	return g;
}

double *SubMatrix(double *M, int dim, int n1, int n2, int m1, int m2 )
{
	double *Msub = (double *) malloc((n2-n1+1)*(m2-m1+1)*sizeof(double));
	if (!Msub) printf("malloc fail at Msub in SubMatrix. size = %d\n",(n2-n1+1)*(m2-m1+1));
	int i,j;
	for ( i = 0; i < n2-n1+1; i++)
	{
		for ( j = 0; j < m2-m1+1; j++)
		{
			Msub[i*(m2-m1+1) + j] = M[(i+n1)*dim + (j+m1)];
	        //printf("Msub %d\n",(n2-n1+1)*(m2-m1+1)-1-(i*(m2-m1+1) + j));
	        //fflush(stdout);
		}
	}
	return Msub;
}

vector<double> SubMatrix(vector<double> M, int dim, int n1, int n2, int m1, int m2 )
{
	vector<double> Msub;

	int i,j;

	for ( i = 0; i < n2-n1+1; i++)
	{
		for ( j = 0; j < m2-m1+1; j++)
		{
			Msub.insert(Msub.end(),M[(i+n1)*dim + (j+m1)]);
	        //printf("Msub %d\n",(n2-n1+1)*(m2-m1+1)-1-(i*(m2-m1+1) + j));
	        //fflush(stdout);
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

void vAdd(vector<double> a,const double aleph, vector<double> &b)
{
	int i;
	for ( i = 0; i < b.size(); i++)
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

void svMul(const double aleph, vector<double> &a)
{
	int i;
	for(i = 0; i < a.size(); i++)
	{
		a[i] *= aleph;
	}
}

void CalcContinuumWFs(double **basisfuns, const unsigned int *indices, const unsigned int Nindices, double **contfuns, const int Ncont,const int NGridpoints, double *U, double *R)
{
	int i,j,dim,n, iU,jU;

	double **retfuns;
	retfuns = (double **) malloc(Ncont*sizeof(double *));
	if (!retfuns) printf("malloc fail at retfuns in CalcContinuumWFs. size = %d\n",Ncont);

	n = NGridpoints * 2;
	dim = Ncont + Nindices;

	for ( i = Ncont-1; i >= 0; i--)
	{
		iU = i + Nindices;
		svMul(U[iU*dim+ iU],contfuns[i],n);
		for ( j = i-1; j >= 0; j--)
		{
			jU = j + Nindices;
			vAdd(contfuns[j],U[jU*dim + iU],contfuns[i],n);
		}
		for( j = 0; j < Nindices; j++)
		{
			vAdd(basisfuns[indices[j]],U[j*dim + iU],contfuns[i],n);
		}
	}
	//PrintMartix(R,Ncont);
	//PrintArray(contfuns[0],n);
	//printf("start Rotating\n");
	for ( i = 0; i < Ncont; i++)
	{
		retfuns[i] = (double *) malloc(n*sizeof(double));
		if (!retfuns[i]) printf("malloc fail at retfuns[%d] in CalcContinuumWFs. size = %d\n",i,n);
		for( j = 0; j < n; j++)
		{
			retfuns[i][j] = 0;
		}
		for ( j = 0; j < Ncont; j++)
		{
			vAdd(contfuns[j],R[j*Ncont + i],retfuns[i],n);
		}
	}
	for ( i = 0; i < Ncont; i++)
	{
		//free(contfuns[i]);
		contfuns[i] = retfuns[i];
	}
	//free(contfuns);*/
}

void CalcContinuumWFs(vector<vector<double> > basisfuns, const unsigned int *indices, const unsigned int Nindices,
		vector<vector<double> > &contfuns, const int Ncont,const int NGridpoints, vector<double> U, vector<double> R)
{
	int i,j,dim,n, iU,jU;

	vector<vector<double> > retfuns;
	vector<double> help;

	//retfuns = (double **) malloc(Ncont*sizeof(double *));

	n = NGridpoints * 2;
	dim = Ncont + Nindices;

	for ( i = Ncont-1; i >= 0; i--)
	{
		iU = i + Nindices;
		svMul(U[iU*dim+ iU],contfuns[i]);
		for ( j = i-1; j >= 0; j--)
		{
			jU = j + Nindices;
			vAdd(contfuns[j],U[jU*dim + iU],contfuns[i]);
		}
		for( j = 0; j < Nindices; j++)
		{
			vAdd(basisfuns[indices[j]],U[j*dim + iU],contfuns[i]);
		}
	}
	//PrintMartix(R,Ncont);
	//PrintArray(contfuns[0],n);
	//printf("start Rotating\n");
	for ( i = 0; i < Ncont; i++)
	{
		retfuns.insert(retfuns.end(),help);

		for( j = 0; j < n; j++)
		{
			retfuns[i].insert(retfuns[i].end(),0.);
		}
		for ( j = 0; j < Ncont; j++)
		{
			vAdd(contfuns[j],R[i*Ncont + j],retfuns[i]);
		}
	}

	for ( i = 0; i < Ncont; i++)
	{
		contfuns[i] = retfuns[i];
	}
}

double *CreateBaseForContinuum(vector<vector<double> > basisfuns,vector<vector<double> > derivfuns,vector<double> x,
		SKappaGroup group, vector<vector<double> > &contfuns, vector<vector<double> > h0Cont)
{
	int l,k,dim,info;
    bool QRelativistic = true;
	//unsigned int NBasisFuns = (unsigned int) basisfuns.size();
	unsigned int Ncont = (unsigned int) contfuns.size();
	uint32_t NGridpoints = (uint32_t) x.size();

	char blasRight = 'R';
	char blasLeft  = 'L';
	//char blasUpper = 'U';
	char blasLower = 'L';
	char blasTrans = 'T';
	char blasNoTrans = 'N';
	char blasNonUnit = 'N';

	double *tau,*d,*e;
	vector<double> T, g, Tsub;
	vector<vector<double> > dualfuns, dualcont, Ddualfuns, Ddualcont;

	k = group.kappa;
	l = 0;

	dim = Ncont + group.n;


	CreateDualFunctions(basisfuns,x,dualfuns,QRelativistic);

	CreateDualFunctions(contfuns,x,dualcont,QRelativistic);

    T = CreateEKinCB(basisfuns,derivfuns, dualfuns,dualcont, group.indices, group.n,x,k,h0Cont);


	g = MetricCB(basisfuns,dualfuns,group.indices, group.n,dualcont,Ncont,NGridpoints,QRelativistic);

	//char cU = 'U';
	char cL = 'L';
	char cN = 'N';
	char cV = 'V';
	double EINS = 1;
	info = 1;
	dpotrf(&cL,(unsigned int *) &dim,&(g[0]),(unsigned int *) &dim,&info);
	dtrtri(&cL,&cN,(unsigned int *) &dim,&(g[0]),(unsigned int *) &dim,&info);
	//dtrmm(&blasLeft,&blasUpper,&blasTrans,&blasNonUnit,(unsigned int *) &dim,(unsigned int *) &dim,&EINS,&(g[0]),(unsigned int *) &dim,&(T[0]),(unsigned int *) &dim);
	//dtrmm(&blasRight,&blasUpper,&blasNoTrans,&blasNonUnit,(unsigned int *) &dim, (unsigned int *) &dim,&EINS,&(g[0]),(unsigned int *) &dim,&(T[0]),(unsigned int *) &dim);
	dtrmm(&blasLeft,&blasLower,&blasNoTrans,&blasNonUnit,(unsigned int *) &dim,(unsigned int *) &dim,&EINS,&(g[0]),(unsigned int *) &dim,&(T[0]),(unsigned int *) &dim);
	dtrmm(&blasRight,&blasLower,&blasTrans,&blasNonUnit,(unsigned int *) &dim, (unsigned int *) &dim,&EINS,&(g[0]),(unsigned int *) &dim,&(T[0]),(unsigned int *) &dim);

	Tsub = SubMatrix(T,dim,group.n,dim-1,group.n,dim-1);

	d = (double *) malloc(Ncont*sizeof(double));
	if (!d) printf("malloc fail at d in CreateBaseForContinuum. size = %d\n",Ncont);
	e = (double *) malloc((Ncont-1)*sizeof(double));
	if (!e) printf("malloc fail at e in CreateBaseForContinuum. size = %d\n",Ncont-1);
	tau = (double *) malloc((Ncont-1)*sizeof(double));
	if (!tau) printf("malloc fail at tau in CreateBaseForContinuum. size = %d\n",Ncont-1);

	double *work = (double *) malloc((2*Ncont-2)*sizeof(double));
	if (!work) printf("malloc fail at work in CreateBaseForContinuum. size = %d\n",Ncont);
	unsigned int lwork = (unsigned int) (2*Ncont-2);

	dsytrd(&cL,(unsigned int *) &Ncont,&(Tsub[0]),(unsigned int *) &Ncont,d,e,tau,work, &lwork,&info);
	if ( 0 != info) printf("error: argument %d of LAPACKE_dsytrd",-info);
	dorgtr(&cL,(unsigned int *) &Ncont,&(Tsub[0]),(unsigned int *) &Ncont,tau,work,&lwork,&info);
	if ( 0 != info) printf("error: argument %d of LAPACKE_dorgtr",-info);
	dsteqr(&cV,(unsigned int *) &Ncont,d,e,&(Tsub[0]),(unsigned int *) &Ncont,work,&info);
	if ( 0 != info) printf("error: argument %d of LAPACKE_dsteqr",-info);

	CalcContinuumWFs(basisfuns, group.indices, group.n, contfuns, Ncont, NGridpoints, g, Tsub);

	free(work);
	free(tau);
	free(e);

	return d;

}

double *CreateBaseForContinuum(vector<vector<double> > basisfuns,vector<vector<double> > derivfuns,vector<double> x,
		SKappaGroup group, vector<vector<double> > &contfuns, vector<vector<double> > &Dcontfuns, vector<vector<double> > h0Cont, const int Z = 0)
{
	int k,dim,info;
    bool QRelativistic = false;

	//unsigned int NBasisFuns = (unsigned int) basisfuns.size();
	unsigned int Ncont = (unsigned int) contfuns.size();
	uint32_t NGridpoints = (uint32_t) x.size();

	char blasRight = 'R';
	char blasLeft  = 'L';
	//char blasUpper = 'U';
	char blasLower = 'L';
	char blasTrans = 'T';
	char blasNoTrans = 'N';
	char blasNonUnit = 'N';

	double *tau,*d,*e;
	vector<double> T, g, Tsub;
	vector<vector<double> > dualfuns, dualcont, Ddualfuns, Ddualcont;

	k = group.kappa;
	
    uint32_t l;
    if ( k > 0 )
	{
		l = (uint32_t) k;
	}
	else
	{
		l = (uint32_t)(-k -1);
	}
    

	dim = Ncont + group.n;


	CreateDualFunctions(basisfuns,x,dualfuns,QRelativistic);

	CreateDualFunctions(contfuns,x,dualcont,QRelativistic);
    
    CreateDualFunctions(Dcontfuns,x,Ddualcont,QRelativistic);
    
    CreateDualFunctions(derivfuns,x,Ddualfuns,QRelativistic);

    CreateEKinCBNonRelativistic(dualcont,Ddualcont, basisfuns,derivfuns,dualfuns,Ddualfuns, x,l,h0Cont,T,Z);


	g = MetricCB(basisfuns,dualfuns,group.indices, group.n,dualcont,Ncont,NGridpoints,QRelativistic);

	//char cU = 'U';
	char cL = 'L';
	char cN = 'N';
	char cV = 'V';
	double EINS = 1;
	info = 1;
	dpotrf(&cL,(unsigned int *) &dim,&(g[0]),(unsigned int *) &dim,&info);
	dtrtri(&cL,&cN,(unsigned int *) &dim,&(g[0]),(unsigned int *) &dim,&info);
	//dtrmm(&blasLeft,&blasUpper,&blasTrans,&blasNonUnit,(unsigned int *) &dim,(unsigned int *) &dim,&EINS,&(g[0]),(unsigned int *) &dim,&(T[0]),(unsigned int *) &dim);
	//dtrmm(&blasRight,&blasUpper,&blasNoTrans,&blasNonUnit,(unsigned int *) &dim, (unsigned int *) &dim,&EINS,&(g[0]),(unsigned int *) &dim,&(T[0]),(unsigned int *) &dim);
	dtrmm(&blasLeft,&blasLower,&blasNoTrans,&blasNonUnit,(unsigned int *) &dim,(unsigned int *) &dim,&EINS,&(g[0]),(unsigned int *) &dim,&(T[0]),(unsigned int *) &dim);
	dtrmm(&blasRight,&blasLower,&blasTrans,&blasNonUnit,(unsigned int *) &dim, (unsigned int *) &dim,&EINS,&(g[0]),(unsigned int *) &dim,&(T[0]),(unsigned int *) &dim);

	Tsub = SubMatrix(T,dim,group.n,dim-1,group.n,dim-1);

	d = (double *) malloc(Ncont*sizeof(double));
	if (!d) printf("malloc fail at d in CreateBaseForContinuum. size = %d\n",Ncont);
	e = (double *) malloc((Ncont-1)*sizeof(double));
	if (!e) printf("malloc fail at e in CreateBaseForContinuum. size = %d\n",Ncont-1);
	tau = (double *) malloc((Ncont-1)*sizeof(double));
	if (!tau) printf("malloc fail at tau in CreateBaseForContinuum. size = %d\n",Ncont-1);

	double *work = (double *) malloc((2*Ncont-2)*sizeof(double));
	if (!work) printf("malloc fail at work in CreateBaseForContinuum. size = %d\n",Ncont);
	unsigned int lwork = (unsigned int) (2*Ncont-2);

	dsytrd(&cL,(unsigned int *) &Ncont,&(Tsub[0]),(unsigned int *) &Ncont,d,e,tau,work, &lwork,&info);
	if ( 0 != info) printf("error: argument %d of LAPACKE_dsytrd",-info);
	dorgtr(&cL,(unsigned int *) &Ncont,&(Tsub[0]),(unsigned int *) &Ncont,tau,work,&lwork,&info);
	if ( 0 != info) printf("error: argument %d of LAPACKE_dorgtr",-info);
	dsteqr(&cV,(unsigned int *) &Ncont,d,e,&(Tsub[0]),(unsigned int *) &Ncont,work,&info);
	if ( 0 != info) printf("error: argument %d of LAPACKE_dsteqr",-info);

	CalcContinuumWFs(basisfuns, group.indices, group.n, contfuns, Ncont, NGridpoints, g, Tsub);

	free(work);
	free(tau);
	free(e);

	return d;

}


vector<double> CreateEKinCB(vector<vector<double> > basisfuns, vector<vector<double> > derivativefuns, vector<vector<double> > dualfuns,
		vector<vector<double> > dualcont, const unsigned int *indices, const unsigned int Nindices, vector<double> x, const int kappa,vector<vector<double> > h0cont)
{
	const int NcontFuns = (int) dualcont.size();
	unsigned int NGridPoints = (unsigned int) x.size();

	vector<vector<double> >  Tbb;
	vector<vector<double> >  Tcb;
	vector<double> T;
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
	Tbb = CalculateDiracMatrixElements(basisfuns,derivativefuns,dualfuns,indices,Nindices,x,NGridPoints,kappa,Z,speedoflight);
	Tcb = CalculateContinuumBoundDiracElements(basisfuns,derivativefuns,dualcont,NcontFuns,indices,Nindices,x,NGridPoints,kappa,Z);

	dim = Nindices + NcontFuns;
	T = vector<double>(dim*dim);

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
			T[i*dim + j] = h0cont[i-Nindices][j-Nindices];
		}
	}

	return T;
}

void CalculateH0MatrixElements(const vector<vector<double> > &psi, const vector<vector<double> > &dpsi,const vector<vector<double> > &phi,const vector<vector<double> > &dphi, vector<double> x,uint32_t l, const int Z, vector<vector<double> > &T)
{
    T = vector<vector<double> >(psi.size());
    #pragma omp parallel for
	for (uint32_t i = 0; i < (uint32_t) psi.size(); i++ )
	{
		T[i] = vector<double>(phi.size());
		for (uint32_t j = 0; j < (uint32_t) phi.size(); j++)
		{
			T[i][j] = H0nonRelativistic(psi[i],dpsi[i],phi[j],dphi[j], x, l, Z);
		}
	}
}

double H0nonRelativistic(const vector<double> psi,const vector<double> &dpsi,const vector<double> &phi,const vector<double> &dphi,vector<double> x,const uint32_t l, const int Z)
{
	vector<double> Tphi, Vphi;
    

	uint32_t i, eins, n;
	eins = 1;
    n = (uint32_t) x.size();
    Tphi = vector<double>(x.size());
    Vphi = vector<double>(x.size());
	if ( x[0] == 0. )
	{
		i = 1;
	}
	else
	{
		i = 0;
	}
	while ( i < (uint32_t) x.size())
	{
		double r = x[i];
		Vphi[i] = (l*(l+1)/(2*r)-Z)/r * phi[i];
		Tphi[i] = dphi[i]/2;
		i++;
	}
	return ddot(&n,&(dpsi[0]),&eins,&(Tphi[0]),&eins) + ddot(&n,&(psi[0]),&eins,&(Vphi[0]),&eins);
}

vector<double> CreateEKinCBNonRelativistic(const vector<vector<double> > &dualcont, const vector<vector<double> > &Ddualcont, const vector<vector<double> > &phi,const vector<vector<double> > &dphi
        , const vector<vector<double> > &dualphi,const vector<vector<double> > &Ddualphi, const vector<double> &x, const uint32_t l, const vector<vector<double> > &h0cont, vector<double> &T,const int Z = 0)
{
	const uint32_t NcontFuns = (uint32_t) dualcont.size();
	//uint32_t NGridPoints = (uint32_t) x.size();

	vector<vector<double> >  Tbb;
	vector<vector<double> >  Tcb;
	uint32_t i,j, dim, Nbound;

	CalculateH0MatrixElements(dualphi,Ddualphi,phi,dphi,x,l,Z,Tbb);
	CalculateH0MatrixElements(dualcont,Ddualcont,phi,dphi,x,l,Z,Tcb);
    
    Nbound = (uint32_t) Tbb.size();
    
	dim = Nbound + NcontFuns;
	T = vector<double>(dim*dim);

	for ( i = 0; i < Nbound; i++)
	{

		for ( j = i; j < Nbound; j++)
		{
			T[i*dim + j] = Tbb[i][j];
			T[j*dim + i] = Tbb[i][j];
		}
	}
	for ( i = Nbound; i < dim; i++)
	{
		for ( j = 0; j < (uint32_t) Tcb[i-Nbound].size(); j++)
		{
			T[i*dim + j] = Tcb[i-Nbound][j];
			T[j*dim + i] = Tcb[i-Nbound][j];
		}
		for ( j = Nbound; j < dim; j++)
		{
			T[i*dim + j] = h0cont[i-Nbound][j-Nbound];
		}
	}

	return T;
}

vector<vector<double> > CreateEigenDifferentialEnergies(uint32_t n, double Emin, double dE)
{
	vector<vector<double> > h0;
	double q1,q2;
	h0 = vector<vector<double> >(n);
	for (uint32_t i = 0; i < n; i++)
	{
		h0[i] = vector<double>(n);
		for (uint32_t j = 0; j < n; j++)
		{
			if ( i != j)
			{
				h0[i][j] = 0;
			}
			else
			{
				q1 = sqrt(2*(Emin +  i    * dE));
				q2 = sqrt(2*(Emin + (i+1) * dE));
				h0[i][j] = (q1*q1 + q2*q2 +q1*q2)/6;
			}
		}
	}
	return h0;
}





















