/*
 * SlaterTuple.cpp
 *
 *  Created on: Jun 12, 2018
 *      Author: brass
 */

#include <iostream>
#include <cmath>
#include <map>
#include <vector>
//#include "myarray.h"

//using namespace std;
// inclusion of namespace std leads to conflicts with function atomic in lua

class myarray
{
public:

	vector<int> val;

	myarray()
	{
		val = vector<int>(4);
	}
/*	myarray(initializer_list<int> l) : val(l)
	{
		;
	}
*/
	virtual ~myarray()
	{
		;
	}
	void Assign(int n1,int n2,int n3,int n4)
	{
		if (val.size()==4)
		{
			val[0] = n1;
			val[1] = n2;
			val[2] = n3;
			val[3] = n4;
		}
	}

	int &operator[](int i)
	{
		return val[i];
	}

};
bool operator<(const myarray& a, const myarray& b)
{
	for(int i = 0; i < a.val.size(); i++)
	{
		if (a.val[i] != b.val[i])
			return a.val[i] < b.val[i];
	}
	return false;
}

int SSign(const int x)
{
    return (x > 0) - (x < 0);
}

int TupleContributes(myarray k)
{
	int k0,k1,l1,l2,l3,l4, lsum;

	l1 = abs(k[0]) + (SSign(k[0])-1)/2;
	l2 = abs(k[1]) + (SSign(k[1])-1)/2;
	l3 = abs(k[2]) + (SSign(k[2])-1)/2;
	l4 = abs(k[3]) + (SSign(k[3])-1)/2;
	lsum = l1 + l2 + l3 + l4;

	k0 = fmax( abs( abs(k[0])-abs(k[2]) ), abs( abs(k[1])-abs(k[3]) ) );
	k1 = fmin( abs(k[0])+abs(k[2]), abs(k[1])+abs(k[3]) ) -1;
	if ( (abs(k[0])+abs(k[2]) + (SSign(k[0])+SSign(k[2]))/2 + k0) % 2 == 0)
	{
		k0++;
	}
	return (k0 <= k1) && (lsum % 2 == 0);

}

vector<myarray> CreateKappaTuple(vector<int> kappas)
{
	int i,j,m,l;
	myarray k;
	vector<myarray> result;
	for ( i = 0; i < kappas.size(); i++)
		for ( j = 0; j < kappas.size(); j++)
			for ( m = 0; m < kappas.size(); m++)
				for ( l = 0; l < kappas.size(); l++)
				{
					k.Assign(kappas[i],kappas[j],kappas[m],kappas[l]);
					if (TupleContributes(k))
					{
						result.insert(result.end(),k);
					}
				}
	return result;
}

int NotContains(vector<myarray> tau,std::map<myarray,myarray> tuple)
{
	int i;
	for ( i = 0; i < tau.size(); i++)
	{
		if (tuple.find(tau[i]) != tuple.end())
		{
			return 0;
		}
	}
	return 1;
}

std::map<myarray,myarray> CreateSlaterTuple(vector<int> indToK)
{
	int i,j,m,l;
	myarray k;
	vector<myarray> tau = vector<myarray>(8);
  std::map<myarray,myarray> result;

	for ( i = 0; i < indToK.size(); i++)
		for ( j = 0; j < indToK.size(); j++)
			for ( m = 0; m < indToK.size(); m++)
				for ( l = 0; l < indToK.size(); l++)
				{
					k.Assign(indToK[i],indToK[j],indToK[m],indToK[l]);
					//if (i==0 && j == 0 && m ==18 && l==19)
						//cout << indToK[i]<< "\t" <<indToK[j]<< "\t" <<indToK[m]<< "\t" <<indToK[l]<<endl;
					if (TupleContributes(k))
					{

						tau[0].Assign(i, j, m, l);
						tau[1].Assign(m, j, i, l);
						tau[2].Assign(i, l, m, j);
						tau[3].Assign(m, l, i, j);
						tau[4].Assign(j, i, l, m);
						tau[5].Assign(j, m, l, i);
						tau[6].Assign(l, i, j, m);
						tau[7].Assign(l, m, j, i);

						/*if (i==0 && j == 0 && m ==18 && l==19)
						{
							cout << NotContains(tau, result) << "\t" << result[tau[0]][1] << endl;
						}*/

						if (NotContains(tau, result))
						{
							result[tau[0]] = k;
						}
					}
				}


	return result;
}
