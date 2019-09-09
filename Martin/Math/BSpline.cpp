/*
 * BSpline.cpp
 *
 *  Created on: Apr 24, 2019
 *      Author: brass
 */



namespace BSpline
{

double Evaluate(uint32_t i, uint32_t n, const vector<double> &x, double t)
{
	double s1 = 0.;
	double s2 = 0.;

	if (n==0)
	{
		if (x[i] <= t && t < x[i+1])
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
	uint32_t eins = 1;
	double one = 1.;
	double null = 0.;
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


