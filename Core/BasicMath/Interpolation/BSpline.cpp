/*
 * BSpline.cpp
 *
 *  Created on: Apr 24, 2019
 *      Author: brass
 */

#include "BSpline.h"
#include "TBSpline.cpp"

namespace BSpline
{

double Evaluate(uint32_t i, uint32_t n, const vector<double> &x, double t)
{
	double s1 = 0.;
	double s2 = 0.;

	if (n==0)
	{
		if ( (x[i] <= t && t < x[i+1]) || (t == x[i+1] && t == x[x.size()-1] ) ) //TODO: Find a better way for this
		{
			return 1.;
		}
		else
		{
			return 0.;
		}
	}
	else
	{
		if (x[i+n] != x[i])
		{
			s1 = (t - x[i])/(x[i+n]-x[i]) * Evaluate(i, n-1,x,t);
		}
		if (x[i+n+1] != x[i+1])
		{
			s2 = (x[i+n+1] - t)/(x[i+n+1]-x[i+1]) * Evaluate(i+1,n-1,x,t);
		}
		return s1 + s2;
	}
}

double deBoor(uint32_t idx, double t, const vector<double> &x, uint32_t n, const vector<double> &c)
{
	if (idx < n)
	{
		printf("Error deBoor: idx < n\n");
		return NAN;
	}
	vector<double> d = vector<double>(n+1);
	for (uint32_t i = 0; i < n+1; i++)
	{
		d[i] = c[i + idx -n];
	}
	for (uint32_t i = 1; i < n+1; i++)
	{
		for (uint32_t j = n; j > i-1; j--)
		{
			double alpha = (t - x[j+idx-n])/(x[j+1+idx-i] - x[j+idx-n]);
			d[j] = (1. - alpha) * d[j-1] + alpha *d[j];
		}
	}

	return d[n];
}

double EvalDerivative(uint32_t i, uint32_t n, const vector<double> &x,  double t, uint32_t k)
{
	if (k > 0)
	{
		double E1 =  EvalDerivative(i,n-1,x,t,k-1);
		double E2 = EvalDerivative(i+1,n-1,x,t,k-1);
		double dx1 = (x[i + n] -x[i]);
		double dx2 = (x[i + n +1] -x[i+1]);
		if (dx1 == 0)
		{
			if (dx2 == 0)
			{
				return 0;
			}
			else
			{
				return - (n / dx2 * E2);
			}
		}
		else
		{
			if (dx2 == 0)
			{
				return n / dx1 * E1;
			}
			else
			{
				return n / dx1 * E1 - (n / dx2 * E2);
			}
		}

	}
	else
	{
		return Evaluate(i,n,x,t);
	}
}

double Derivative(uint32_t idx, double t, const vector<double> &x, uint32_t n, const vector<double> &c, uint32_t k)
{
	double val = 0;
	for (uint32_t i = 0; i <= (uint32_t)c.size()-1; i++) //Max(0,(int)(idx-2*n)) Min(idx +k-n,(uint32_t)c.size()-1)
	{
		val += c[i] * EvalDerivative(i,n,x,t,k);
	}
	return val;
}

void GenerateGrid(const std::vector<double> &t,std::vector<double> &x, uint32_t n)
{
	uint32_t d = (uint32_t) t.size();
	if (n == 0)
	{
		printf("Error BSpline: Interpolation with n = 0 not implemented\n");
		return;
	}
	if (d < (n+1))
	{
		printf("Error BSpline::GenerateGrid: number of points = %u < %u = n+1\n",d,n+1);
		return;
	}
	x = vector<double>(d+n+1);
	for (uint32_t i = 0; i <= n; i++)
	{
		x[i] = t[0];
		x[i + d] = t[d-1];
	}
	for (uint32_t i = 0; i < d - (n+1); i++)
	{
		x[i + n+1] = 0;
		for (uint32_t j = i+1; j < i+1 + n; j++ )
		{
			x[i + n+1] += t[j];
		}
		x[i + n+1] /= n;
	}

}
void GenerateFisherGrid(int m, double hMax,double rMax, double Z, uint32_t n,std::vector<double> &x)
{
	double h = pow(2.,-m);
	uint32_t ks = n+1;
	uint32_t i;
	for (i = 1; i <= ks; i++)
	{
		x.push_back(0.);
	}
	for (i = ks; i <= ks + m; i++)
	{
		x.push_back(x[i-1] + h);
	}
	while (h * x[i-1] < hMax)
	{
		x.push_back(x[i-1]*(1.+h));
		i++;
	}
	while (x[i-1] < rMax)
	{
		x.push_back(x[i-1] + hMax);
		i++;
	}
	double last = x[i-1];
	for (i = 1; i <= ks-1; i++)
	{
		x.push_back(last);
	}
	for (uint32_t j = 0; j < (uint32_t) x.size(); j++)
	{
		x[j] /= (double) Z;
	}

}


void GetKnotMatrix(const std::vector<double> &x, const uint32_t n, const std::vector<double> &t, std::vector<double> &B)
{
	uint32_t d = (uint32_t) t.size();
	uint32_t N = (uint32_t) x.size();
	if (N != d + n + 1)
	{
		printf("Error GetKnotMatrix\n");
		return;
	}
	B = vector<double>(d*d);
	for (uint32_t i = 0; i < d; i++)
	{
		for (uint32_t j = 0; j < d; j++)
		{
			B[i +j*d] = Evaluate(j,n,x,t[i]);
		}
	}
	//if (t[d-1] == x[N-1]) B[d*d-1] = 1;
}

#define TriIdx(i,j) (i + j*(j+1)/2)

void CalculateCellIntegrals(const vector<double> &x, uint32_t ks,uint32_t k, vector<double> &OffDiagInts1, vector<double> &OffDiagInts2, vector<double> &DiagInts)
{
	uint32_t nv = (uint32_t) x.size() - 2 * ks - 1 ;
	uint32_t nk = (ks+1)*((ks+1) + 1) / 2;
	uint32_t nSplines = (uint32_t) x.size() - (ks + 1);

	OffDiagInts1 = vector<double>(nv * nk);
	OffDiagInts2 = vector<double>(nv * nk);
	DiagInts 	 = vector<double>(nv * (ks+1)*(ks+1)*(ks+1)*(ks+1));

	vector<double> BSplinesOnQuadPoints = vector<double>(nSplines * (ks+1)*(ks+1));
	for (uint32_t i = 0; i < nSplines; i++)
	{
		for(uint32_t iv = 0; iv < ks+1; iv++)
		{
			for (uint32_t j = 0; j < ks+1; j++)
			{
				BSplinesOnQuadPoints[i*nSplines*(ks+1) + iv*(ks+1) + j] = 0; // TODO insert stuff here
			}
		}
	}

	for (uint32_t iv = 0; iv < nv; iv++)
	{
		for (uint32_t i = 0; i <= ks; i++ )
		{
			for (uint32_t j = i; j <= ks; j++)
			{
				OffDiagInts1[iv * nk + TriIdx(i,j)] = 0;	// TODO insert stuff here
				OffDiagInts2[iv * nk + TriIdx(i,j)] = 0;	// TODO insert stuff here

			}
		}
	}

}

double Overlap(uint32_t i, uint32_t j, uint32_t ks, const vector<double> &x)
{
	uint32_t N = (uint32_t) x.size();
	uint32_t ia = std::max(ks,std::max(i,j));
	uint32_t ib = std::min(N - ks -1, std::min(i,j) + ks+1);

	double integral = 0;
	for (uint32_t iv = ia; iv < ib; iv++)
	{
		double a = x[iv  ];
		double b = x[iv+1];
		vector<double> u, w;
		Quad::GaussianQuadratureWeights(ks+1,u,w,a,b);
		for (uint32_t m = 0; m < ks+1; m++)
		{
			integral += Evaluate(i,ks,x,u[m])*Evaluate(j,ks,x,u[m])*w[m];
		}
	}
	return integral;
}

void OverlapMatrix(uint32_t ks, vector<double> x, vector<double> &S)
{
	uint32_t Nspline = (uint32_t)x.size() - ks - 1 -2;
	S = vector<double>(Nspline*(ks+1));
	//S = vector<double>(Nspline*Nspline);
/*	for (int k = 0; k < Nspline*Nspline; k++)
	{
		S[k] = 100;
	}*/
	for (uint32_t j = 0; j < Nspline; j++)
	{
		for ( uint32_t i = j; i <= std::min(Nspline-1,j+ks); i++)  // n = #sub-diagonals
		{
			//S[(i-j)*Nspline + j] = Overlap(i+1,j+1,ks,x);
			S[(i-j) + j*(ks+1)] = Overlap(i+1,j+1,ks,x);
		}
	}
}

void Components(const vector<Interpolation::InterpolatingFunction*> &funs, uint32_t ks, const vector<double> &x, vector<double> &S, vector<double> &c)
{
	uint32_t N = (uint32_t) x.size();
	uint32_t Nfuns = (uint32_t) funs.size();
	uint32_t Nspline = N - ks-3;
	c = vector<double> (Nspline*Nfuns);
	for (uint32_t i = 0; i < Nspline; i++)
	{
		for (uint32_t j = 0; j < Nfuns; j++)
		{
			c[i + j*Nspline] = 0;
			for (uint32_t iv = std::max(ks,i+1); iv < std::min(N-ks-1,i+1+ks); iv++)
			{
				double a = x[iv];
				double b = x[iv+1];
				vector<double> u,w;
				Quad::GaussianQuadratureWeights(ks+1,u,w,a,b);
				for (uint32_t m = 0; m < ks+1; m++)
				{
					c[i + j*Nspline] += (*funs[j])(u[m]) * Evaluate(i,ks,x,u[m]) * w[m];
				}
			}
		}
	}
	printf("c = {");
	for (uint32_t i = 0; i < Nspline*Nfuns; i++)
	{
		printf("%.16f, ",c[i]);
	}
	printf("}\n");
  /*
	uint32_t eins = 1;
	double one = 1.;
	double null = 0.;
  */
	const char *l = "L";
	int info = 0;
	uint32_t LDA = ks+1;
/*	dpbsv(l,&Nspline,&ks,&Nfuns,&S[0],&Nspline,&c[0],&Nspline,&info);
	if (info)
	{
		printf("error in BSpline::Components: \n\tLAPACK dpbsv (positive banded symmetric linear equation) info = %d\n",info);
		fflush(stdout);
	}*/
	printf("???\n");
	dpbsv(l,&Nspline,&ks,&Nfuns,&S[0],&LDA,&c[0],&Nspline,&info);
	if (info)
	{
		printf("error in BSpline::Components: \n\tLAPACK dpbsv (positive banded symmetric linear equation) info = %d\n",info);
		fflush(stdout);
	}


}

/*void OverlapMatrix(uint32_t ks, vector<double> x, vector<double> &S)
{
	uint32_t Nspline = (uint32_t) x.size() - ks - 1 -2; 	// -2 since boundary conditions make first and last spline vanish
	uint32_t N = (uint32_t) x.size();
	S = vector<double>(Nspline*(ks+1));
	vector<double> BSplines = vector<double>(Nspline*(ks+1)*(ks+1));
	vector<double> weights = vector<double>((N -1 -2*ks)*(ks+1));
	for (uint32_t iv = ks; iv < N -1 -ks; iv++)
	{
		vector<double> u, w;
		double a = x[iv  ];
		double b = x[iv+1];
		GaussianQuadratureWeights(ks+1,u,w,a,b);
		for (uint32_t j = 0; j < ks+1; j++)
		{
			weights[(iv-ks)*(ks+1) + j] = w[j];
		}
		for (uint32_t j = std::max((uint32_t) 1,iv-ks); j <= iv; j++)
		{
			uint32_t idx = (j-1)*(ks+1)*(ks+1) + (iv - j)*(ks+1);
			for (uint32_t m = 0; m < ks+1; m++)
			{
				BSplines[idx + m] = Evaluate(j,ks,x,u[m]);
			}
		}

	}
	for (uint32_t j = 0; j < Nspline; j++)
	{
		double integral = 0;
		for (uint32_t i = 0; i < ks+1; i++)
		{
			for (uint32_t m = 0; m < ks+1; m++)
			{
				integral += weights[(i+j)*(ks+1) + m];
			}
		}
	}


}*/


} // end of namespace BSpline


