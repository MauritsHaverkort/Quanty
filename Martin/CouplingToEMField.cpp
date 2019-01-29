/*
 * CouplingToEMField.cpp
 *
 *  Created on: Jan 28, 2019
 *      Author: brass
 */
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
// calculates \int Y^*_{l_1,m_1}Y_{l_2,m_2}Y_{l_3,m_3}d\Omega
double YYYintegral(int l1, int m1,int l2, int m2, int l3, int m3)
{
	if (m1 % 2 == 0)
	{
		return  sqrt(((2*l1+1) + (2*l2+1) + (2*l1+1)) / (4 * M_PI)) * ThreeJSymbol(l1,0,l2,0,l3,0) * ThreeJSymbol(l1,-m1,l2,m2,l3,m3);
	}
	else
	{
		return -sqrt(((2*l1+1) + (2*l2+1) + (2*l1+1)) / (4 * M_PI)) * ThreeJSymbol(l1,0,l2,0,l3,0) * ThreeJSymbol(l1,-m1,l2,m2,l3,m3);
	}
}

double ClebschHalf(int twoj1, int twom1,int twoj2, int twom2, int twoj3, int twom3)
{
	if ((twoj2-twoj1-twom3) % 4 == 0)
	{
		return  sqrt(twoj3 + 1) * ThreeJSymbolHalf(twoj1,twom1,twoj2,twom2,twoj3,-twom3);
	}
	else
	{
		return -sqrt(twoj3 + 1) * ThreeJSymbolHalf(twoj1,twom1,twoj2,twom2,twoj3,-twom3);
	}
}

double DipoleMatrixElement(double deltaE, double grg, double frf, uint32_t mu, int kappa1, int kappa2,int twom1,int twom2)
{
	int twoj1 = 2*abs(kappa1) - 1;
	int twoj2 = 2*abs(kappa2) - 1;
	int l1A, l1B, l2A, l2B;
	if (kappa1 < 0)
	{
		l1A = (twoj1 - 1)/2;
		l1B = (twoj1 + 1)/2;
	}
	else
	{
		l1A = (twoj1 + 1)/2;
		l1B = (twoj1 - 1)/2;
	}
	if (kappa2 < 0)
	{
		l2A = (twoj2 - 1)/2;
		l2B = (twoj2 + 1)/2;
	}
	else
	{
		l2A = (twoj2 + 1)/2;
		l2B = (twoj2 - 1)/2;
	}

	double result, resultG,resultF;
	switch (mu) {
		case 0:

			break;
		case 1:
			resultG  = ClebschHalf(2*l1A,twom1-1,1, 1,twoj1,twom1) * ClebschHalf(2*l2A,twom2-1,1, 1,twoj2,twom2)
					 * (YYYintegral(l1A,(twom1-1)/2,1,-1,l2A,(twom2-1)/2) - YYYintegral(l1A,(twom1-1)/2,1,1,l2A,(twom2-1)/2));
			resultG += ClebschHalf(2*l1A,twom1+1,1,-1,twoj1,twom1) * ClebschHalf(2*l2A,twom2+1,1,-1,twoj2,twom2)
					 * (YYYintegral(l1A,(twom1+1)/2,1,-1,l2A,(twom2+1)/2) - YYYintegral(l1A,(twom1+1)/2,1,1,l2A,(twom2+1)/2));
			resultG *= grg;

			resultF  = ClebschHalf(2*l1B,twom1-1,1, 1,twoj1,twom1) * ClebschHalf(2*l2B,twom2-1,1, 1,twoj2,twom2)
					 * (YYYintegral(l1B,(twom1-1)/2,1,-1,l2B,(twom2-1)/2) - YYYintegral(l1B,(twom1-1)/2,1,1,l2B,(twom2-1)/2));
			resultF += ClebschHalf(2*l1B,twom1+1,1,-1,twoj1,twom1) * ClebschHalf(2*l2B,twom2+1,1,-1,twoj2,twom2)
					 * (YYYintegral(l1B,(twom1+1)/2,1,-1,l2B,(twom2+1)/2) - YYYintegral(l1B,(twom1+1)/2,1,1,l2B,(twom2+1)/2));
			resultF *= frf;

			return sqrt(2 * M_PI / 3) * deltaE * (resultG + resultF);
			break;
		case 2:
			resultG  = ClebschHalf(2*l1A,twom1-1,1, 1,twoj1,twom1) * ClebschHalf(2*l2A,twom2-1,1, 1,twoj2,twom2)
					 * (YYYintegral(l1A,(twom1-1)/2,1,-1,l2A,(twom2-1)/2) + YYYintegral(l1A,(twom1-1)/2,1,1,l2A,(twom2-1)/2));
			resultG += ClebschHalf(2*l1A,twom1+1,1,-1,twoj1,twom1) * ClebschHalf(2*l2A,twom2+1,1,-1,twoj2,twom2)
					 * (YYYintegral(l1A,(twom1+1)/2,1,-1,l2A,(twom2+1)/2) + YYYintegral(l1A,(twom1+1)/2,1,1,l2A,(twom2+1)/2));
			resultG *= grg;

			resultF  = ClebschHalf(2*l1B,twom1-1,1, 1,twoj1,twom1) * ClebschHalf(2*l2B,twom2-1,1, 1,twoj2,twom2)
					 * (YYYintegral(l1B,(twom1-1)/2,1,-1,l2B,(twom2-1)/2) + YYYintegral(l1B,(twom1-1)/2,1,1,l2B,(twom2-1)/2));
			resultF += ClebschHalf(2*l1B,twom1+1,1,-1,twoj1,twom1) * ClebschHalf(2*l2B,twom2+1,1,-1,twoj2,twom2)
					 * (YYYintegral(l1B,(twom1+1)/2,1,-1,l2B,(twom2+1)/2) + YYYintegral(l1B,(twom1+1)/2,1,1,l2B,(twom2+1)/2));
			resultF *= frf;

			return sqrt(2 * M_PI / 3) * deltaE * (resultG + resultF);
			break;
		case 3:
			resultG  = ClebschHalf(2*l1A,twom1-1,1, 1,twoj1,twom1) * ClebschHalf(2*l2A,twom2-1,1, 1,twoj2,twom2)
					 * (YYYintegral(l1A,(twom1-1)/2,1,0,l2A,(twom2-1)/2));
			resultG += ClebschHalf(2*l1A,twom1+1,1,-1,twoj1,twom1) * ClebschHalf(2*l2A,twom2+1,1,-1,twoj2,twom2)
					 * (YYYintegral(l1A,(twom1+1)/2,1,0,l2A,(twom2+1)/2));
			resultG *= grg;

			resultF  = ClebschHalf(2*l1B,twom1-1,1, 1,twoj1,twom1) * ClebschHalf(2*l2B,twom2-1,1, 1,twoj2,twom2)
					 * (YYYintegral(l1B,(twom1-1)/2,1,0,l2B,(twom2-1)/2));
			resultF += ClebschHalf(2*l1B,twom1+1,1,-1,twoj1,twom1) * ClebschHalf(2*l2B,twom2+1,1,-1,twoj2,twom2)
					 * (YYYintegral(l1B,(twom1+1)/2,1,0,l2B,(twom2+1)/2));
			resultF *= frf;

			return sqrt(4 * M_PI / 3) * deltaE * (resultG + resultF);
			break;
		default:
			printf("Are you kidding me?\n mu is a covariant index, so it should take values between zero and three, but you gave me %d.\n What do you expect me to do with that?\n I could invent some extra dimensions...\nActually that would be fun to do.\n",mu);
			return NAN;
			break;
	}

}

double ParseFPLOEnergies(const char *filename, const string &orbital, const string &j)
{
	std::ifstream file;

	file.open(filename);
	if(!file.is_open())
	{
		printf("Error: could not open %s\n",filename);
		fflush(stdout);
		return NAN;
	}
	else
	{
		std::string word;
		double val;
		while(!file.eof() && word != "ENERGIES")
		{
			file >> word;
		}
		while(!file.eof() && word != j)
		{
			while(!file.eof() && word != orbital)
			{
				file >> word;
			}
			file >> word;
		}
		if(file.eof())
		{
			printf("Error while parsing FPLO: could not find %s%s\n",orbital.c_str(),j.c_str());
			file.close();
			return NAN;
		}
		else
		{
			file >> word;
			file.close();
			return atof(word.c_str());
		}
	}
}

void GetFPLOEnergies(const char *filename, const std::vector<std::string> &orbitalnames, std::vector<double> &energies)
{
	std::map<std::string,double> dictEnergies;
	std::ifstream file;

	file.open(filename);
	if(!file.is_open())
	{
		printf("Error: could not open %s\n",filename);
		fflush(stdout);
	}
	else
	{
		std::string word, appendix;
		double val;
		while(!file.eof() && word != "ENERGIES")
		{
			file >> word;
		}
		file >> word;file >> word;file >> word;file >> word;file >> word;file >> word;file >> word;file >> word;file >> word;file >> word;file >> word;file >> word;file >> word;
		file >> word;

		file >> appendix;
		word.append(appendix);
		file >> appendix;
		val = atof(appendix.c_str());
		dictEnergies[word] = val;
		file >> word;file >> word;file >> word;file >> word;file >> word;file >> word;file >> word;file >> word;file >> word;
		cout << val << endl;
	}



}


