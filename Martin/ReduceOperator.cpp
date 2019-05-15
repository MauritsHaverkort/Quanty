/*
 * ReduceOperator.cpp
 *
 *  Created on: Apr 29, 2019
 *      Author: brass
 */

template<typename T>
bool QContains(const vector<T> &v, T x)
{
	for (uint32_t i = 0; i < (uint32_t) v.size(); i++ )
	{
		if (v[i] == x)
		{
			return true;
		}
	}
	return false;
}

bool OppFixedLengthContributes(const uint16_t* opp, uint32_t length,const vector<uint16_t> &ind)
{
	for (uint32_t i = 0; i < length; i++)
	{
		uint16_t cidx = NotCrAnMask & opp[i];
		if (QContains(ind, cidx))
		{
			int counter = (cidx == opp[i]) ? -1 : 1;
			for (uint32_t j = 0; j < length; j++)
			{
				if (cidx == (NotCrAnMask & opp[j]) && j != i)
				{
					counter = (cidx == opp[j]) ? counter -1 : counter + 1;
				}
			}
			if (counter != 0)
			{
				return false;
			}
		}

	}
	return true;
}

void ReduceOperator(const OperatorType &Opp, OperatorType* ReducedOpp, vector<uint16_t> ind)
{
	InitOperator(ReducedOpp);
	ReducedOpp->NFermion = Opp.NFermion;
	ReducedOpp->NBoson = Opp.NBoson;
	for (uint32_t l = 0; l <= Opp.MaxLength; l++)
	{
		InitOperatorFixedLength(ReducedOpp,l,Opp.l[l].NBitsKey);
		for (uint32_t i = 0; i < Opp.l[l].N; i++)
		{
			if (OppFixedLengthContributes(&(Opp.l[l].o[i*l]),l,ind))
			{
				if (Opp.l[l].QComplex)
				{
					ComplexOperatorAddLaderNormalOrder(ReducedOpp, &(Opp.l[l].o[i*l]), l, Opp.l[l].V[i].Re, Opp.l[l].V[i].Im, &QDetNotOrderedFcFaBcBa);
				}
				else
				{
					/*uint16_t oadd[l];
					for (int j = 0; j < l; j++)
					{
						oadd[j] = Opp.l[l].o[i*l + j];
					}*/
					RealOperatorAddLaderNormalOrder(ReducedOpp, &(Opp.l[l].o[i*l]), l, Opp.l[l].VRe[i], &QDetNotOrderedFcFaBcBa);
					//printf("i = %u\n",i);
					//fflush(stdout);
					//RealOperatorAddLaderNormalOrder(ReducedOpp, oadd, l, Opp.l[l].VRe[i], &QDetNotOrderedFcFaBcBa);
				}
			}
		}
	}
}


