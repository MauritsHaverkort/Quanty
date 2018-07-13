/*
 * MartinsLittleHelper.cpp
 *
 *  Created on: Jul 10, 2018
 *      Author: brass
 */
#include <vector>

using namespace std;

vector<double> ConvertArray(double *a, const int n)
{
	vector<double> res;
	for (int i = 0; i < n; i++)
	{
		res.insert(res.end(),a[i]);
	}
	return res;
}

vector<vector<double> > ConvertMatrix(double **m, const int n1, const int n2)
{
	vector<vector<double> > res;
	vector<double> help;
	for (int i = 0; i < n1; i++)
	{
		res.insert(res.end(),help);
		for (int j = 0; j < n2; j++)
		{
			res[i].insert(res[i].end(),m[i][j]);
		}
	}
	return res;
}


