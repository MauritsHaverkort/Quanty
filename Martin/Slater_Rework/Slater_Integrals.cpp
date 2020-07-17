/*
 * Slater_Integrals.cpp
 *
 *  Created on: Jul 17, 2020
 *      Author: brass
 */


double Slater_Integrate(const uint32_t K, const vector<double> &x
										, const vector<double> &G1
										, const vector<double> &G2
										, const vector<double> &G3
										, const vector<double> &G4)
{
	double w1, w2, gg;
	double val = 0.;
	uint32_t N = (uint32_t) x.size();

	for (uint32_t i = 1; i < N-1; i++)
	{
		w1 = (x[i+1]-x[i-1])/2.;
		val += w1 * w1 * G1[i] * G2[i] * G3[i] * G4[i] / x[i];
		w1 = w1 * pow(x[i],K);
		for (uint32_t j = i+1; j < N-1; j++)
		{
			w2 = (x[j+1]-x[j-1])/2.;
			gg = G1[i] * G2[j] * G3[i] * G4[j] + G1[j] * G2[i] * G3[j] * G4[i];
			val += w1 * w2 / pow(x[j],K+1) * gg;
		}
	}
	return val;
}


