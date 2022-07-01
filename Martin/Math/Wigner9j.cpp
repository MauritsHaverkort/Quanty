long BinomialCoefficientLong(int n, int k)
{
	if(n<k) return 0;
	long Bnk = 1;
	for (int m = 1; m <= k; m++)
	{
		Bnk *= n + 1 - m;
		Bnk /= m;
	}
	return Bnk;
}

long WeiSquareBrackets(int two_m1, int two_m2, int two_m3, int two_m4, int two_m5, int two_m6)
{
	long result = 0; //TODO: check below that start end end are always integers this should be guaranteed if WeiDelta !=0 -> check this
	int start_sum = std::max({two_m1+two_m5+two_m6, two_m2+two_m4+two_m6, two_m3+two_m4+two_m5, two_m1+two_m2+two_m3})/2;
	int end_sum   = std::min({two_m1+two_m2+two_m4+two_m5, two_m1+two_m3+two_m4+two_m6, two_m2+two_m3+two_m5+two_m6})/2;
	long sign = (start_sum % 2) == 0 ? 1 : -1;
	for (int n = start_sum; n <= end_sum; n++)
	{
		result += sign 	* BinomialCoefficientLong(n+1, n-(two_m1+two_m5+two_m6)/2)
						* BinomialCoefficientLong((two_m1+two_m5-two_m6)/2, n-(two_m2+two_m4+two_m6)/2)
						* BinomialCoefficientLong((two_m1-two_m5+two_m6)/2, n-(two_m3+two_m4+two_m5)/2)
						* BinomialCoefficientLong((two_m5+two_m6-two_m1)/2, n-(two_m1+two_m2+two_m3)/2);
		sign *= -1;
	}
	return result;
}

double WeiDelta(int two_a, int two_b, int two_c)
{
	if ((two_c < std::abs(two_a-two_b)) || (two_c > two_a + two_b) || ((two_a + two_b + two_c)%2 != 0)) return 0;

	double result = 1.;
	result *= Factorial((two_a+two_b-two_c)/2) * Factorial((two_a-two_b+two_c)/2)
			* Factorial((two_b+two_c-two_a)/2)/ Factorial((two_a+two_b+two_c)/2+1); //TODO: use cancellation of factors
	return std::sqrt(result);
}

double Wigner9J(int two_j1, int two_j2, int two_j3, int two_j4, int two_j5, int two_j6, int two_j7, int two_j8, int two_j9)
{
	double prefactor = WeiDelta(two_j1, two_j2, two_j3) * WeiDelta(two_j4, two_j5, two_j6) * WeiDelta(two_j7, two_j8, two_j9)
					* WeiDelta(two_j1, two_j4, two_j7) * WeiDelta(two_j2, two_j5, two_j8) * WeiDelta(two_j3, two_j6, two_j9);
	if (prefactor == 0.) return 0.;

	double result = 0.;
	int start_sum = std::max({std::abs(two_j8 - two_j4),std::abs(two_j2-two_j6),std::abs(two_j1-two_j9)});
	int end_sum   = std::min({two_j8 + two_j4, two_j2 + two_j6, two_j1 + two_j9});
	int sign = (start_sum % 2) == 0 ? 1 : -1;

	for (int two_k = start_sum; two_k <= end_sum; two_k += 2)
	{
		result += (((two_k + 1.) 	* WeiSquareBrackets(two_j1, two_j2, two_j3, two_j6, two_j9, two_k))
										* WeiSquareBrackets(two_j6, two_j4, two_j5, two_j8, two_j2, two_k))
										* WeiSquareBrackets(two_j8, two_j9, two_j7, two_j1, two_j4, two_k);

	}
	return sign * prefactor * result;
}
