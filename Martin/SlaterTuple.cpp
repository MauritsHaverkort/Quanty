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
#include <initializer_list>


using namespace std;

class myarray
{
public:

	vector<int> val;

	myarray()
	{
		;
	}
	myarray(initializer_list<int> l) : val(l)
	{
		;
	}
	virtual ~myarray()
	{
		;
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
					k = {kappas[i],kappas[j],kappas[m],kappas[l]};
					if (TupleContributes(k))
					{
						result.insert(result.end(),k);
					}
				}
	return result;
}

int NotContains(array<myarray,8> tau,map<myarray,myarray> tuple)
{
	int i;
	for ( i = 0; i < 8; i++)
	{
		if (tuple.find(tau[i]) != tuple.end())
		{
			return 0;
		}
	}
	return 1;
}

map<myarray,myarray> CreateSlaterTuple(vector<int> indToK)
{
	int i,j,m,l;
	myarray k;
	array<myarray,8> tau;
	map<myarray,myarray> result;

	for ( i = 0; i < indToK.size(); i++)
		for ( j = 0; j < indToK.size(); j++)
			for ( m = 0; m < indToK.size(); m++)
				for ( l = 0; l < indToK.size(); l++)
				{
					k = {indToK[i],indToK[j],indToK[m],indToK[l]};
					//if (i==0 && j == 0 && m ==18 && l==19)
						//cout << indToK[i]<< "\t" <<indToK[j]<< "\t" <<indToK[m]<< "\t" <<indToK[l]<<endl;
					if (TupleContributes(k))
					{

						tau[0] = {i, j, m, l};
						tau[1] = {m, j, i, l};
						tau[2] = {i, l, m, j};
						tau[3] = {m, l, i, j};
						tau[4] = {j, i, l, m};
						tau[5] = {j, m, l, i};
						tau[6] = {l, i, j, m};
						tau[7] = {l, m, j, i};

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



/*

int main(int argc, char **argv)
{

	vector<myarray> tuple;
	vector<int> kappas = {-1,1,-2,2,-3,3,-4};

	tuple = CreateKappaTuple(kappas);

	map<int, vector<int>> indices;
	indices[-1] = {0,1,2};
	indices[1] = {3,4,5};
	myarray k = {-1,1,3,-4};
	//cout << TupleContributes(k) << endl;

	map<myarray,myarray> st;
	vector<int> indToK = {-1,-1,-1,-1,-1,-1,1,1,1,1,-2,-2,-2,-2,2,2,-3,-3,3,-4};
	st = CreateSlaterTuple(indToK);
	for (auto it : st)
	{
		cout << it.first.val[0] << "\t"<< it.first.val[1] << "\t"<< it.first.val[2] << "\t"<< it.first.val[3] << endl;
	}


}
*/


