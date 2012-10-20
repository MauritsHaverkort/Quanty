

/*  Chebyshev polinomials
*   first kind T
*   T0 (x) = 1
*   T-1(x) = x
*   T1 (x) = x
*   Tm+1(x)= 2 x Tm(x) - Tm-1(x) 
*/
double ChebyshevT(const double x, const int n)
{
  return (cos(n*acos(x)));
}

/*  second kind U
*   U0 (x) = 1
*   U-1(x) = 0
*   U1 (x) = 2x
*   Um+1(x)= 2 x Um(x) - Um-1(x)
*/
double ChebyshevU(const double x, const int n)
{
  return ((sin((n+1)*acos(x)))/(sin(acos(x))));
}




/*  TesseralHarmonicZ defines the normalized tesseral harmonic function 
*   The Tesseral Harmonics are defined as
*   m = 0 : Z[k,m,theta,phi] = Y[k,m,theta,phi]
*   m > 0 : Z[k,m,theta,phi] =   sqrt(0.5) * ( Y[k,-m,theta,phi] + (-1)^m Y[k, m,theta,phi] )
*   m < 0 : Z[k,m,theta,phi] = i sqrt(0.5) * ( Y[k, m,theta,phi] - (-1)^m Y[k,-m,theta,phi] )
*
*   for m=0 there is no phi dependence
*   for m>0 there is a cos(m*phi) dependence
*   for m<0 there is a sin(m*phi) dependence
*/
double TesseralHarmonicZ(const int k, const int m, const double t, const double p)
{
  switch (k)
  {
    case 0 : switch (m) 
    {
      case 0 : return 0.28209479177387814347403972578; break;
    }
    case 1 : switch (m) 
    {
      case -1 : return 0.244301255951459960793192311419*(cos(p - t) - cos(p + t)); break;
      case 0 : return 0.488602511902919921586384622838*cos(t); break;
      case 1 : return -0.244301255951459960793192311419*(sin(p - t) - sin(p + t)); break;
    }
    case 2 : switch (m) 
    {
      case -2 : return 0.1365685538240098838179232132253*(2*sin(2*p) - sin(2*p - 2*t) - sin(2*p + 2*t)); break;
      case -1 : return 0.273137107648019767635846426451*(cos(p - 2*t) - cos(p + 2*t)); break;
      case 0 : return 0.1576957826262600030154468451479*(1 + 3*cos(2*t)); break;
      case 1 : return -0.273137107648019767635846426451*(sin(p - 2*t) - sin(p + 2*t)); break;
      case 2 : return 0.1365685538240098838179232132253*(2*cos(2*p) - cos(2*p - 2*t) - cos(2*p + 2*t)); break;
    }
    case 3 : switch (m) 
    {
      case -3 : return -0.0737554487408304387932012846927*(cos(3*p - 3*t) - 3*cos(3*p - t) + 3*cos(3*p + t) - cos(3*p + 3*t)); break;
      case -2 : return -0.1806632151650346284618362596499*(sin(2*p - 3*t) - sin(2*p - t) - sin(2*p + t) + sin(2*p + 3*t)); break;
      case -1 : return 0.0571307249330582170197525871146*(5*cos(p - 3*t) + cos(p - t) - cos(p + t) - 5*cos(p + 3*t)); break;
      case 0 : return 0.0932940831475288478535989782997*(3*cos(t) + 5*cos(3*t)); break;
      case 1 : return -0.0571307249330582170197525871146*(5*sin(p - 3*t) + sin(p - t) - sin(p + t) - 5*sin(p + 3*t)); break;
      case 2 : return -0.1806632151650346284618362596499*(cos(2*p - 3*t) - cos(2*p - t) - cos(2*p + t) + cos(2*p + 3*t)); break;
      case 3 : return 0.0737554487408304387932012846927*(sin(3*p - 3*t) - 3*sin(3*p - t) + 3*sin(3*p + t) - sin(3*p + 3*t)); break;
    }
    case 4 : switch (m) 
    {
      case -4 : return 0.0391147334655735084116650327253*(6*sin(4*p) + sin(4*p - 4*t) - 4*sin(4*p - 2*t) - 4*sin(4*p + 2*t) + sin(4*p + 4*t)); break;
      case -3 : return -0.110633173111245658189801927039*(cos(3*p - 4*t) - 2*cos(3*p - 2*t) + 2*cos(3*p + 2*t) - cos(3*p + 4*t)); break;
      case -2 : return 0.0295679592424237505653962834652*(6*sin(2*p) - 7*sin(2*p - 4*t) + 4*sin(2*p - 2*t) + 4*sin(2*p + 2*t) - 7*sin(2*p + 4*t)); break;
      case -1 : return 0.041815408972330572997007024357*(7*cos(p - 4*t) + 2*cos(p - 2*t) - 2*cos(p + 2*t) - 7*cos(p + 4*t)); break;
      case 0 : return 0.01322319336440053797534561214596*(9 + 20*cos(2*t) + 35*cos(4*t)); break;
      case 1 : return -0.041815408972330572997007024357*(7*sin(p - 4*t) + 2*sin(p - 2*t) - 2*sin(p + 2*t) - 7*sin(p + 4*t)); break;
      case 2 : return 0.0295679592424237505653962834652*(6*cos(2*p) - 7*cos(2*p - 4*t) + 4*cos(2*p - 2*t) + 4*cos(2*p + 2*t) - 7*cos(2*p + 4*t)); break;
      case 3 : return 0.110633173111245658189801927039*(sin(3*p - 4*t) - 2*sin(3*p - 2*t) + 2*sin(3*p + 2*t) - sin(3*p + 4*t)); break;
      case 4 : return 0.0391147334655735084116650327253*(6*cos(4*p) + cos(4*p - 4*t) - 4*cos(4*p - 2*t) - 4*cos(4*p + 2*t) + cos(4*p + 4*t)); break;
    }
    case 5 : switch (m) 
    {
      case -5 : return 0.0205119392762553157129412114493*(cos(5*p - 5*t) - 5*cos(5*p - 3*t) + 10*cos(5*p - t) - 10*cos(5*p + t) + 5*cos(5*p + 3*t) - cos(5*p + 5*t)); break;
      case -4 : return 0.0648644473400325399686187038311*(sin(4*p - 5*t) - 3*sin(4*p - 3*t) + 2*sin(4*p - t) + 2*sin(4*p + t) - 3*sin(4*p + 3*t) + sin(4*p + 5*t)); break;
      case -3 : return -0.01528869685735157461512527237009*(9*cos(3*p - 5*t) - 13*cos(3*p - 3*t) - 6*cos(3*p - t) + 6*cos(3*p + t) + 13*cos(3*p + 3*t) - 9*cos(3*p + 5*t)); break;
      case -2 : return -0.0748990122652081836719047053057*(3*sin(2*p - 5*t) - sin(2*p - 3*t) - 2*sin(2*p - t) - 2*sin(2*p + t) - sin(2*p + 3*t) + 3*sin(2*p + 5*t)); break;
      case -1 : return 0.01415458284986552879057630181929*(21*cos(p - 5*t) + 7*cos(p - 3*t) + 2*cos(p - t) - 2*cos(p + t) - 7*cos(p + 3*t) - 21*cos(p + 5*t)); break;
      case 0 : return 0.00730939515333897477748219950564*(30*cos(t) + 35*cos(3*t) + 63*cos(5*t)); break;
      case 1 : return -0.01415458284986552879057630181929*(21*sin(p - 5*t) + 7*sin(p - 3*t) + 2*sin(p - t) - 2*sin(p + t) - 7*sin(p + 3*t) - 21*sin(p + 5*t)); break;
      case 2 : return -0.0748990122652081836719047053057*(3*cos(2*p - 5*t) - cos(2*p - 3*t) - 2*cos(2*p - t) - 2*cos(2*p + t) - cos(2*p + 3*t) + 3*cos(2*p + 5*t)); break;
      case 3 : return 0.01528869685735157461512527237009*(9*sin(3*p - 5*t) - 13*sin(3*p - 3*t) - 6*sin(3*p - t) + 6*sin(3*p + t) + 13*sin(3*p + 3*t) - 9*sin(3*p + 5*t)); break;
      case 4 : return 0.0648644473400325399686187038311*(cos(4*p - 5*t) - 3*cos(4*p - 3*t) + 2*cos(4*p - t) + 2*cos(4*p + t) - 3*cos(4*p + 3*t) + cos(4*p + 5*t)); break;
      case 5 : return -0.0205119392762553157129412114493*(sin(5*p - 5*t) - 5*sin(5*p - 3*t) + 10*sin(5*p - t) - 10*sin(5*p + t) + 5*sin(5*p + 3*t) - sin(5*p + 5*t)); break;
    }
    case 6 : switch (m) 
    {
      case -6 : return 0.01067475164362366128085636211369*(20*sin(6*p) - sin(6*p - 6*t) + 6*sin(6*p - 4*t) - 15*sin(6*p - 2*t) - 15*sin(6*p + 2*t) + 6*sin(6*p + 4*t) - sin(6*p + 6*t)); break;
      case -5 : return 0.0369784244098711254998084069824*(cos(5*p - 6*t) - 4*cos(5*p - 4*t) + 5*cos(5*p - 2*t) - 5*cos(5*p + 2*t) + 4*cos(5*p + 4*t) - cos(5*p + 6*t)); break;
      case -4 : return 0.00788382657388631498835123218438*(20*sin(4*p) + 11*sin(4*p - 6*t) - 26*sin(4*p - 4*t) + 5*sin(4*p - 2*t) + 5*sin(4*p + 2*t) - 26*sin(4*p + 4*t) + 11*sin(4*p + 6*t)); break;
      case -3 : return -0.014393832179920679674086965566*(11*cos(3*p - 6*t) - 12*cos(3*p - 4*t) - 9*cos(3*p - 2*t) + 9*cos(3*p + 2*t) + 12*cos(3*p + 4*t) - 11*cos(3*p + 6*t)); break;
      case -2 : return 0.007196916089960339837043482783*(20*sin(2*p) - 33*sin(2*p - 6*t) + 6*sin(2*p - 4*t) + 17*sin(2*p - 2*t) + 17*sin(2*p + 2*t) + 6*sin(2*p + 4*t) - 33*sin(2*p + 6*t)); break;
      case -1 : return 0.00910345878935517795054810276546*(33*cos(p - 6*t) + 12*cos(p - 4*t) + 5*cos(p - 2*t) - 5*cos(p + 2*t) - 12*cos(p + 4*t) - 33*cos(p + 6*t)); break;
      case 0 : return 0.001986537570863388310400844720642*(50 + 105*cos(2*t) + 126*cos(4*t) + 231*cos(6*t)); break;
      case 1 : return -0.00910345878935517795054810276546*(33*sin(p - 6*t) + 12*sin(p - 4*t) + 5*sin(p - 2*t) - 5*sin(p + 2*t) - 12*sin(p + 4*t) - 33*sin(p + 6*t)); break;
      case 2 : return 0.007196916089960339837043482783*(20*cos(2*p) - 33*cos(2*p - 6*t) + 6*cos(2*p - 4*t) + 17*cos(2*p - 2*t) + 17*cos(2*p + 2*t) + 6*cos(2*p + 4*t) - 33*cos(2*p + 6*t)); break;
      case 3 : return 0.014393832179920679674086965566*(11*sin(3*p - 6*t) - 12*sin(3*p - 4*t) - 9*sin(3*p - 2*t) + 9*sin(3*p + 2*t) + 12*sin(3*p + 4*t) - 11*sin(3*p + 6*t)); break;
      case 4 : return 0.00788382657388631498835123218438*(20*cos(4*p) + 11*cos(4*p - 6*t) - 26*cos(4*p - 4*t) + 5*cos(4*p - 2*t) + 5*cos(4*p + 2*t) - 26*cos(4*p + 4*t) + 11*cos(4*p + 6*t)); break;
      case 5 : return -0.0369784244098711254998084069824*(sin(5*p - 6*t) - 4*sin(5*p - 4*t) + 5*sin(5*p - 2*t) - 5*sin(5*p + 2*t) + 4*sin(5*p + 4*t) - sin(5*p + 6*t)); break;
      case 6 : return 0.01067475164362366128085636211369*(20*cos(6*p) - cos(6*p - 6*t) + 6*cos(6*p - 4*t) - 15*cos(6*p - 2*t) - 15*cos(6*p + 2*t) + 6*cos(6*p + 4*t) - cos(6*p + 6*t)); break;
    }
    case 7 : switch (m) 
    {
      case -7 : return -0.00552470884784840764239112475841*(cos(7*p - 7*t) - 7*cos(7*p - 5*t) + 21*cos(7*p - 3*t) - 35*cos(7*p - t) + 35*cos(7*p + t) - 21*cos(7*p + 3*t) + 7*cos(7*p + 5*t) - cos(7*p + 7*t)); break;
      case -6 : return -0.0206715676703273454841421913203*(sin(6*p - 7*t) - 5*sin(6*p - 5*t) + 9*sin(6*p - 3*t) - 5*sin(6*p - t) - 5*sin(6*p + t) + 9*sin(6*p + 3*t) - 5*sin(6*p + 5*t) + sin(6*p + 7*t)); break;
      case -5 : return 0.00405402795875203374803866898944*(13*cos(5*p - 7*t) - 43*cos(5*p - 5*t) + 33*cos(5*p - 3*t) + 25*cos(5*p - t) - 25*cos(5*p + t) - 33*cos(5*p + 3*t) + 43*cos(5*p + 5*t) - 13*cos(5*p + 7*t)); break;
      case -4 : return 0.00810805591750406749607733797889*(13*sin(4*p - 7*t) - 25*sin(4*p - 5*t) - 3*sin(4*p - 3*t) + 15*sin(4*p - t) + 15*sin(4*p + t) - 3*sin(4*p + 3*t) - 25*sin(4*p + 5*t) + 13*sin(4*p + 7*t)); break;
      case -3 : return -0.001222335420799172137864532272666*(143*cos(3*p - 7*t) - 121*cos(3*p - 5*t) - 117*cos(3*p - 3*t) - 45*cos(3*p - t) + 45*cos(3*p + t) + 117*cos(3*p + 3*t) + 121*cos(3*p + 5*t) - 143*cos(3*p + 7*t)); break;
      case -2 : return -0.001728643329863213409215655109985*(143*sin(2*p - 7*t) - 11*sin(2*p - 5*t) - 57*sin(2*p - 3*t) - 75*sin(2*p - t) - 75*sin(2*p + t) - 57*sin(2*p + 3*t) - 11*sin(2*p + 5*t) + 143*sin(2*p + 7*t)); break;
      case -1 : return 0.000705715684238416517429047889257*(429*cos(p - 7*t) + 165*cos(p - 5*t) + 81*cos(p - 3*t) + 25*cos(p - t) - 25*cos(p + t) - 81*cos(p + 3*t) - 165*cos(p + 5*t) - 429*cos(p + 7*t)); break;
      case 0 : return 0.001066941826750077217327525103323*(175*cos(t) + 189*cos(3*t) + 231*cos(5*t) + 429*cos(7*t)); break;
      case 1 : return -0.000705715684238416517429047889257*(429*sin(p - 7*t) + 165*sin(p - 5*t) + 81*sin(p - 3*t) + 25*sin(p - t) - 25*sin(p + t) - 81*sin(p + 3*t) - 165*sin(p + 5*t) - 429*sin(p + 7*t)); break;
      case 2 : return -0.001728643329863213409215655109985*(143*cos(2*p - 7*t) - 11*cos(2*p - 5*t) - 57*cos(2*p - 3*t) - 75*cos(2*p - t) - 75*cos(2*p + t) - 57*cos(2*p + 3*t) - 11*cos(2*p + 5*t) + 143*cos(2*p + 7*t)); break;
      case 3 : return 0.001222335420799172137864532272666*(143*sin(3*p - 7*t) - 121*sin(3*p - 5*t) - 117*sin(3*p - 3*t) - 45*sin(3*p - t) + 45*sin(3*p + t) + 117*sin(3*p + 3*t) + 121*sin(3*p + 5*t) - 143*sin(3*p + 7*t)); break;
      case 4 : return 0.00810805591750406749607733797889*(13*cos(4*p - 7*t) - 25*cos(4*p - 5*t) - 3*cos(4*p - 3*t) + 15*cos(4*p - t) + 15*cos(4*p + t) - 3*cos(4*p + 3*t) - 25*cos(4*p + 5*t) + 13*cos(4*p + 7*t)); break;
      case 5 : return -0.00405402795875203374803866898944*(13*sin(5*p - 7*t) - 43*sin(5*p - 5*t) + 33*sin(5*p - 3*t) + 25*sin(5*p - t) - 25*sin(5*p + t) - 33*sin(5*p + 3*t) + 43*sin(5*p + 5*t) - 13*sin(5*p + 7*t)); break;
      case 6 : return -0.0206715676703273454841421913203*(cos(6*p - 7*t) - 5*cos(6*p - 5*t) + 9*cos(6*p - 3*t) - 5*cos(6*p - t) - 5*cos(6*p + t) + 9*cos(6*p + 3*t) - 5*cos(6*p + 5*t) + cos(6*p + 7*t)); break;
      case 7 : return 0.00552470884784840764239112475841*(sin(7*p - 7*t) - 7*sin(7*p - 5*t) + 21*sin(7*p - 3*t) - 35*sin(7*p - t) + 35*sin(7*p + t) - 21*sin(7*p + 3*t) + 7*sin(7*p + 5*t) - sin(7*p + 7*t)); break;
    }
    case 8 : switch (m) 
    {
      case -8 : return 0.00284736976630792917527513627171*(70*sin(8*p) + sin(8*p - 8*t) - 8*sin(8*p - 6*t) + 28*sin(8*p - 4*t) - 56*sin(8*p - 2*t) - 56*sin(8*p + 2*t) + 28*sin(8*p + 4*t) - 8*sin(8*p + 6*t) + sin(8*p + 8*t)); break;
      case -7 : return -0.01138947906523171670110054508683*(cos(7*p - 8*t) - 6*cos(7*p - 6*t) + 14*cos(7*p - 4*t) - 14*cos(7*p - 2*t) + 14*cos(7*p + 2*t) - 14*cos(7*p + 4*t) + 6*cos(7*p + 6*t) - cos(7*p + 8*t)); break;
      case -6 : return 0.00207942486742008818090300760464*(70*sin(6*p) - 15*sin(6*p - 8*t) + 64*sin(6*p - 6*t) - 84*sin(6*p - 4*t) - 84*sin(6*p + 4*t) + 64*sin(6*p + 6*t) - 15*sin(6*p + 8*t)); break;
      case -5 : return 0.01347621336757073444342924115379*(5*cos(5*p - 8*t) - 14*cos(5*p - 6*t) + 6*cos(5*p - 4*t) + 10*cos(5*p - 2*t) - 10*cos(5*p + 2*t) - 6*cos(5*p + 4*t) + 14*cos(5*p + 6*t) - 5*cos(5*p + 8*t)); break;
      case -4 : return 0.001868814549841135536135986971331*(70*sin(4*p) + 65*sin(4*p - 8*t) - 104*sin(4*p - 6*t) - 36*sin(4*p - 4*t) + 40*sin(4*p - 2*t) + 40*sin(4*p + 2*t) - 36*sin(4*p + 4*t) - 104*sin(4*p + 6*t) + 65*sin(4*p + 8*t)); break;
      case -3 : return -0.00482525841912321904672613206427*(39*cos(3*p - 8*t) - 26*cos(3*p - 6*t) - 30*cos(3*p - 4*t) - 18*cos(3*p - 2*t) + 18*cos(3*p + 2*t) + 30*cos(3*p + 4*t) + 26*cos(3*p + 6*t) - 39*cos(3*p + 8*t)); break;
      case -2 : return 0.001781844759511365119420229894973*(70*sin(2*p) - 143*sin(2*p - 8*t) + 44*sin(2*p - 4*t) + 64*sin(2*p - 2*t) + 64*sin(2*p + 2*t) + 44*sin(2*p + 4*t) - 143*sin(2*p + 8*t)); break;
      case -1 : return 0.000425942366792109201273625464506*(715*cos(p - 8*t) + 286*cos(p - 6*t) + 154*cos(p - 4*t) + 70*cos(p - 2*t) - 70*cos(p + 2*t) - 154*cos(p + 4*t) - 286*cos(p + 6*t) - 715*cos(p + 8*t)); break;
      case 0 : return 0.0000709903944653515335456042440843*(1225 + 2520*cos(2*t) + 2772*cos(4*t) + 3432*cos(6*t) + 6435*cos(8*t)); break;
      case 1 : return -0.000425942366792109201273625464506*(715*sin(p - 8*t) + 286*sin(p - 6*t) + 154*sin(p - 4*t) + 70*sin(p - 2*t) - 70*sin(p + 2*t) - 154*sin(p + 4*t) - 286*sin(p + 6*t) - 715*sin(p + 8*t)); break;
      case 2 : return 0.001781844759511365119420229894973*(70*cos(2*p) - 143*cos(2*p - 8*t) + 44*cos(2*p - 4*t) + 64*cos(2*p - 2*t) + 64*cos(2*p + 2*t) + 44*cos(2*p + 4*t) - 143*cos(2*p + 8*t)); break;
      case 3 : return 0.00482525841912321904672613206427*(39*sin(3*p - 8*t) - 26*sin(3*p - 6*t) - 30*sin(3*p - 4*t) - 18*sin(3*p - 2*t) + 18*sin(3*p + 2*t) + 30*sin(3*p + 4*t) + 26*sin(3*p + 6*t) - 39*sin(3*p + 8*t)); break;
      case 4 : return 0.001868814549841135536135986971331*(70*cos(4*p) + 65*cos(4*p - 8*t) - 104*cos(4*p - 6*t) - 36*cos(4*p - 4*t) + 40*cos(4*p - 2*t) + 40*cos(4*p + 2*t) - 36*cos(4*p + 4*t) - 104*cos(4*p + 6*t) + 65*cos(4*p + 8*t)); break;
      case 5 : return -0.01347621336757073444342924115379*(5*sin(5*p - 8*t) - 14*sin(5*p - 6*t) + 6*sin(5*p - 4*t) + 10*sin(5*p - 2*t) - 10*sin(5*p + 2*t) - 6*sin(5*p + 4*t) + 14*sin(5*p + 6*t) - 5*sin(5*p + 8*t)); break;
      case 6 : return 0.00207942486742008818090300760464*(70*cos(6*p) - 15*cos(6*p - 8*t) + 64*cos(6*p - 6*t) - 84*cos(6*p - 4*t) - 84*cos(6*p + 4*t) + 64*cos(6*p + 6*t) - 15*cos(6*p + 8*t)); break;
      case 7 : return 0.01138947906523171670110054508683*(sin(7*p - 8*t) - 6*sin(7*p - 6*t) + 14*sin(7*p - 4*t) - 14*sin(7*p - 2*t) + 14*sin(7*p + 2*t) - 14*sin(7*p + 4*t) + 6*sin(7*p + 6*t) - sin(7*p + 8*t)); break;
      case 8 : return 0.00284736976630792917527513627171*(70*cos(8*p) + cos(8*p - 8*t) - 8*cos(8*p - 6*t) + 28*cos(8*p - 4*t) - 56*cos(8*p - 2*t) - 56*cos(8*p + 2*t) + 28*cos(8*p + 4*t) - 8*cos(8*p + 6*t) + cos(8*p + 8*t)); break;
    }
    case 9 : switch (m) 
    {
      case -9 : return 0.001462697171588258391699761606466*(cos(9*p - 9*t) - 9*cos(9*p - 7*t) + 36*cos(9*p - 5*t) - 84*cos(9*p - 3*t) + 126*cos(9*p - t) - 126*cos(9*p + t) + 84*cos(9*p + 3*t) - 36*cos(9*p + 5*t) + 9*cos(9*p + 7*t) - cos(9*p + 9*t)); break;
      case -8 : return 0.00620569853311464351035377327894*(sin(8*p - 9*t) - 7*sin(8*p - 7*t) + 20*sin(8*p - 5*t) - 28*sin(8*p - 3*t) + 14*sin(8*p - t) + 14*sin(8*p + t) - 28*sin(8*p + 3*t) + 20*sin(8*p + 5*t) - 7*sin(8*p + 7*t) + sin(8*p + 9*t)); break;
      case -7 : return -0.001064268518250103887881676201882*(17*cos(7*p - 9*t) - 89*cos(7*p - 7*t) + 164*cos(7*p - 5*t) - 84*cos(7*p - 3*t) - 98*cos(7*p - t) + 98*cos(7*p + t) + 84*cos(7*p + 3*t) - 164*cos(7*p + 5*t) + 89*cos(7*p + 7*t) - 17*cos(7*p + 9*t)); break;
      case -6 : return -0.00245782286200696648145050602536*(17*sin(6*p - 9*t) - 63*sin(6*p - 7*t) + 60*sin(6*p - 5*t) + 28*sin(6*p - 3*t) - 42*sin(6*p - t) - 42*sin(6*p + t) + 28*sin(6*p + 3*t) + 60*sin(6*p + 5*t) - 63*sin(6*p + 7*t) + 17*sin(6*p + 9*t)); break;
      case -5 : return 0.000951910701248083128099790019333*(85*cos(5*p - 9*t) - 205*cos(5*p - 7*t) + 36*cos(5*p - 5*t) + 140*cos(5*p - 3*t) + 70*cos(5*p - t) - 70*cos(5*p + t) - 140*cos(5*p + 3*t) - 36*cos(5*p + 5*t) + 205*cos(5*p + 7*t) - 85*cos(5*p + 9*t)); break;
      case -4 : return 0.00796425632564291691801487646517*(17*sin(4*p - 9*t) - 23*sin(4*p - 7*t) - 12*sin(4*p - 5*t) + 4*sin(4*p - 3*t) + 14*sin(4*p - t) + 14*sin(4*p + t) + 4*sin(4*p + 3*t) - 12*sin(4*p + 5*t) - 23*sin(4*p + 7*t) + 17*sin(4*p + 9*t)); break;
      case -3 : return -0.00090177445315662977325012536322*(221*cos(3*p - 9*t) - 117*cos(3*p - 7*t) - 156*cos(3*p - 5*t) - 116*cos(3*p - 3*t) - 42*cos(3*p - t) + 42*cos(3*p + t) + 116*cos(3*p + 3*t) + 156*cos(3*p + 5*t) + 117*cos(3*p + 7*t) - 221*cos(3*p + 9*t)); break;
      case -2 : return -0.001180699911819332935171120437115*(221*sin(2*p - 9*t) + 13*sin(2*p - 7*t) - 52*sin(2*p - 5*t) - 84*sin(2*p - 3*t) - 98*sin(2*p - t) - 98*sin(2*p + t) - 84*sin(2*p + 3*t) - 52*sin(2*p + 5*t) + 13*sin(2*p + 7*t) + 221*sin(2*p + 9*t)); break;
      case -1 : return 0.000125863033500443926419532289375*(2431*cos(p - 9*t) + 1001*cos(p - 7*t) + 572*cos(p - 5*t) + 308*cos(p - 3*t) + 98*cos(p - t) - 98*cos(p + t) - 308*cos(p + 3*t) - 572*cos(p + 5*t) - 1001*cos(p + 7*t) - 2431*cos(p + 9*t)); break;
      case 0 : return 0.0000375251065015101235270913259486*(4410*cos(t) + 4620*cos(3*t) + 5148*cos(5*t) + 6435*cos(7*t) + 12155*cos(9*t)); break;
      case 1 : return -0.000125863033500443926419532289375*(2431*sin(p - 9*t) + 1001*sin(p - 7*t) + 572*sin(p - 5*t) + 308*sin(p - 3*t) + 98*sin(p - t) - 98*sin(p + t) - 308*sin(p + 3*t) - 572*sin(p + 5*t) - 1001*sin(p + 7*t) - 2431*sin(p + 9*t)); break;
      case 2 : return -0.001180699911819332935171120437115*(221*cos(2*p - 9*t) + 13*cos(2*p - 7*t) - 52*cos(2*p - 5*t) - 84*cos(2*p - 3*t) - 98*cos(2*p - t) - 98*cos(2*p + t) - 84*cos(2*p + 3*t) - 52*cos(2*p + 5*t) + 13*cos(2*p + 7*t) + 221*cos(2*p + 9*t)); break;
      case 3 : return 0.00090177445315662977325012536322*(221*sin(3*p - 9*t) - 117*sin(3*p - 7*t) - 156*sin(3*p - 5*t) - 116*sin(3*p - 3*t) - 42*sin(3*p - t) + 42*sin(3*p + t) + 116*sin(3*p + 3*t) + 156*sin(3*p + 5*t) + 117*sin(3*p + 7*t) - 221*sin(3*p + 9*t)); break;
      case 4 : return 0.00796425632564291691801487646517*(17*cos(4*p - 9*t) - 23*cos(4*p - 7*t) - 12*cos(4*p - 5*t) + 4*cos(4*p - 3*t) + 14*cos(4*p - t) + 14*cos(4*p + t) + 4*cos(4*p + 3*t) - 12*cos(4*p + 5*t) - 23*cos(4*p + 7*t) + 17*cos(4*p + 9*t)); break;
      case 5 : return -0.000951910701248083128099790019333*(85*sin(5*p - 9*t) - 205*sin(5*p - 7*t) + 36*sin(5*p - 5*t) + 140*sin(5*p - 3*t) + 70*sin(5*p - t) - 70*sin(5*p + t) - 140*sin(5*p + 3*t) - 36*sin(5*p + 5*t) + 205*sin(5*p + 7*t) - 85*sin(5*p + 9*t)); break;
      case 6 : return -0.00245782286200696648145050602536*(17*cos(6*p - 9*t) - 63*cos(6*p - 7*t) + 60*cos(6*p - 5*t) + 28*cos(6*p - 3*t) - 42*cos(6*p - t) - 42*cos(6*p + t) + 28*cos(6*p + 3*t) + 60*cos(6*p + 5*t) - 63*cos(6*p + 7*t) + 17*cos(6*p + 9*t)); break;
      case 7 : return 0.001064268518250103887881676201882*(17*sin(7*p - 9*t) - 89*sin(7*p - 7*t) + 164*sin(7*p - 5*t) - 84*sin(7*p - 3*t) - 98*sin(7*p - t) + 98*sin(7*p + t) + 84*sin(7*p + 3*t) - 164*sin(7*p + 5*t) + 89*sin(7*p + 7*t) - 17*sin(7*p + 9*t)); break;
      case 8 : return 0.00620569853311464351035377327894*(cos(8*p - 9*t) - 7*cos(8*p - 7*t) + 20*cos(8*p - 5*t) - 28*cos(8*p - 3*t) + 14*cos(8*p - t) + 14*cos(8*p + t) - 28*cos(8*p + 3*t) + 20*cos(8*p + 5*t) - 7*cos(8*p + 7*t) + cos(8*p + 9*t)); break;
      case 9 : return -0.001462697171588258391699761606466*(sin(9*p - 9*t) - 9*sin(9*p - 7*t) + 36*sin(9*p - 5*t) - 84*sin(9*p - 3*t) + 126*sin(9*p - t) - 126*sin(9*p + t) + 84*sin(9*p + 3*t) - 36*sin(9*p + 5*t) + 9*sin(9*p + 7*t) - sin(9*p + 9*t)); break;
    }
    case 10 : switch (m) 
    {
      case -10 : return 0.000749409295138662121642181039548*(252*sin(10*p) - sin(10*p - 10*t) + 10*sin(10*p - 8*t) - 45*sin(10*p - 6*t) + 120*sin(10*p - 4*t) - 210*sin(10*p - 2*t) - 210*sin(10*p + 2*t) + 120*sin(10*p + 4*t) - 45*sin(10*p + 6*t) + 10*sin(10*p + 8*t) - sin(10*p + 10*t)); break;
      case -9 : return 0.00335146025380050237786104889053*(cos(9*p - 10*t) - 8*cos(9*p - 8*t) + 27*cos(9*p - 6*t) - 48*cos(9*p - 4*t) + 42*cos(9*p - 2*t) - 42*cos(9*p + 2*t) + 48*cos(9*p + 4*t) - 27*cos(9*p + 6*t) + 8*cos(9*p + 8*t) - cos(9*p + 10*t)); break;
      case -8 : return 0.000543678645234783595594820110941*(252*sin(8*p) + 19*sin(8*p - 10*t) - 118*sin(8*p - 8*t) + 279*sin(8*p - 6*t) - 264*sin(8*p - 4*t) - 42*sin(8*p - 2*t) - 42*sin(8*p + 2*t) - 264*sin(8*p + 4*t) + 279*sin(8*p + 6*t) - 118*sin(8*p + 8*t) + 19*sin(8*p + 10*t)); break;
      case -7 : return -0.001331735264872856806404841968388*(19*cos(7*p - 10*t) - 88*cos(7*p - 8*t) + 129*cos(7*p - 6*t) - 16*cos(7*p - 4*t) - 98*cos(7*p - 2*t) + 98*cos(7*p + 2*t) + 16*cos(7*p + 4*t) - 129*cos(7*p + 6*t) + 88*cos(7*p + 8*t) - 19*cos(7*p + 10*t)); break;
      case -6 : return 0.0001614966224244323865948608696898*(756*sin(6*p) - 323*sin(6*p - 10*t) + 1054*sin(6*p - 8*t) - 711*sin(6*p - 6*t) - 664*sin(6*p - 4*t) + 266*sin(6*p - 2*t) + 266*sin(6*p + 2*t) - 664*sin(6*p + 4*t) - 711*sin(6*p + 6*t) + 1054*sin(6*p + 8*t) - 323*sin(6*p + 10*t)); break;
      case -5 : return 0.00028889394070211816796725317955*(323*cos(5*p - 10*t) - 680*cos(5*p - 8*t) - 15*cos(5*p - 6*t) + 400*cos(5*p - 4*t) + 350*cos(5*p - 2*t) - 350*cos(5*p + 2*t) - 400*cos(5*p + 4*t) + 15*cos(5*p + 6*t) + 680*cos(5*p + 8*t) - 323*cos(5*p + 10*t)); break;
      case -4 : return 0.000456781427420158382998461480527*(252*sin(4*p) + 323*sin(4*p - 10*t) - 374*sin(4*p - 8*t) - 249*sin(4*p - 6*t) - 8*sin(4*p - 4*t) + 182*sin(4*p - 2*t) + 182*sin(4*p + 2*t) - 8*sin(4*p + 4*t) - 249*sin(4*p + 6*t) - 374*sin(4*p + 8*t) + 323*sin(4*p + 10*t)); break;
      case -3 : return -0.000645986489697729546379443478759*(323*cos(3*p - 10*t) - 136*cos(3*p - 8*t) - 207*cos(3*p - 6*t) - 176*cos(3*p - 4*t) - 98*cos(3*p - 2*t) + 98*cos(3*p + 2*t) + 176*cos(3*p + 4*t) + 207*cos(3*p + 6*t) + 136*cos(3*p + 8*t) - 323*cos(3*p + 10*t)); break;
      case -2 : return 0.0000633441868555005273131042029601*(1764*sin(2*p) - 4199*sin(2*p - 10*t) - 442*sin(2*p - 8*t) + 741*sin(2*p - 6*t) + 1352*sin(2*p - 4*t) + 1666*sin(2*p - 2*t) + 1666*sin(2*p + 2*t) + 1352*sin(2*p + 4*t) + 741*sin(2*p + 6*t) - 442*sin(2*p + 8*t) - 4199*sin(2*p + 10*t)); break;
      case -1 : return 0.0000731435666652423668670944387751*(4199*cos(p - 10*t) + 1768*cos(p - 8*t) + 1053*cos(p - 6*t) + 624*cos(p - 4*t) + 294*cos(p - 2*t) - 294*cos(p + 2*t) - 624*cos(p + 4*t) - 1053*cos(p + 6*t) - 1768*cos(p + 8*t) - 4199*cos(p + 10*t)); break;
      case 0 : return 9.8626765171554764711729175821e-6*(7938 + 16170*cos(2*t) + 17160*cos(4*t) + 19305*cos(6*t) + 24310*cos(8*t) + 46189*cos(10*t)); break;
      case 1 : return -0.0000731435666652423668670944387751*(4199*sin(p - 10*t) + 1768*sin(p - 8*t) + 1053*sin(p - 6*t) + 624*sin(p - 4*t) + 294*sin(p - 2*t) - 294*sin(p + 2*t) - 624*sin(p + 4*t) - 1053*sin(p + 6*t) - 1768*sin(p + 8*t) - 4199*sin(p + 10*t)); break;
      case 2 : return 0.0000633441868555005273131042029601*(1764*cos(2*p) - 4199*cos(2*p - 10*t) - 442*cos(2*p - 8*t) + 741*cos(2*p - 6*t) + 1352*cos(2*p - 4*t) + 1666*cos(2*p - 2*t) + 1666*cos(2*p + 2*t) + 1352*cos(2*p + 4*t) + 741*cos(2*p + 6*t) - 442*cos(2*p + 8*t) - 4199*cos(2*p + 10*t)); break;
      case 3 : return 0.000645986489697729546379443478759*(323*sin(3*p - 10*t) - 136*sin(3*p - 8*t) - 207*sin(3*p - 6*t) - 176*sin(3*p - 4*t) - 98*sin(3*p - 2*t) + 98*sin(3*p + 2*t) + 176*sin(3*p + 4*t) + 207*sin(3*p + 6*t) + 136*sin(3*p + 8*t) - 323*sin(3*p + 10*t)); break;
      case 4 : return 0.000456781427420158382998461480527*(252*cos(4*p) + 323*cos(4*p - 10*t) - 374*cos(4*p - 8*t) - 249*cos(4*p - 6*t) - 8*cos(4*p - 4*t) + 182*cos(4*p - 2*t) + 182*cos(4*p + 2*t) - 8*cos(4*p + 4*t) - 249*cos(4*p + 6*t) - 374*cos(4*p + 8*t) + 323*cos(4*p + 10*t)); break;
      case 5 : return -0.00028889394070211816796725317955*(323*sin(5*p - 10*t) - 680*sin(5*p - 8*t) - 15*sin(5*p - 6*t) + 400*sin(5*p - 4*t) + 350*sin(5*p - 2*t) - 350*sin(5*p + 2*t) - 400*sin(5*p + 4*t) + 15*sin(5*p + 6*t) + 680*sin(5*p + 8*t) - 323*sin(5*p + 10*t)); break;
      case 6 : return 0.0001614966224244323865948608696898*(756*cos(6*p) - 323*cos(6*p - 10*t) + 1054*cos(6*p - 8*t) - 711*cos(6*p - 6*t) - 664*cos(6*p - 4*t) + 266*cos(6*p - 2*t) + 266*cos(6*p + 2*t) - 664*cos(6*p + 4*t) - 711*cos(6*p + 6*t) + 1054*cos(6*p + 8*t) - 323*cos(6*p + 10*t)); break;
      case 7 : return 0.001331735264872856806404841968388*(19*sin(7*p - 10*t) - 88*sin(7*p - 8*t) + 129*sin(7*p - 6*t) - 16*sin(7*p - 4*t) - 98*sin(7*p - 2*t) + 98*sin(7*p + 2*t) + 16*sin(7*p + 4*t) - 129*sin(7*p + 6*t) + 88*sin(7*p + 8*t) - 19*sin(7*p + 10*t)); break;
      case 8 : return 0.000543678645234783595594820110941*(252*cos(8*p) + 19*cos(8*p - 10*t) - 118*cos(8*p - 8*t) + 279*cos(8*p - 6*t) - 264*cos(8*p - 4*t) - 42*cos(8*p - 2*t) - 42*cos(8*p + 2*t) - 264*cos(8*p + 4*t) + 279*cos(8*p + 6*t) - 118*cos(8*p + 8*t) + 19*cos(8*p + 10*t)); break;
      case 9 : return -0.00335146025380050237786104889053*(sin(9*p - 10*t) - 8*sin(9*p - 8*t) + 27*sin(9*p - 6*t) - 48*sin(9*p - 4*t) + 42*sin(9*p - 2*t) - 42*sin(9*p + 2*t) + 48*sin(9*p + 4*t) - 27*sin(9*p + 6*t) + 8*sin(9*p + 8*t) - sin(9*p + 10*t)); break;
      case 10 : return 0.000749409295138662121642181039548*(252*cos(10*p) - cos(10*p - 10*t) + 10*cos(10*p - 8*t) - 45*cos(10*p - 6*t) + 120*cos(10*p - 4*t) - 210*cos(10*p - 2*t) - 210*cos(10*p + 2*t) + 120*cos(10*p + 4*t) - 45*cos(10*p + 6*t) + 10*cos(10*p + 8*t) - cos(10*p + 10*t)); break;
    }
    case 11 : switch (m) 
    {
      case -11 : return -0.000383126028216404728392520692202*(cos(11*p - 11*t) - 11*cos(11*p - 9*t) + 55*cos(11*p - 7*t) - 165*cos(11*p - 5*t) + 330*cos(11*p - 3*t) - 462*cos(11*p - t) + 462*cos(11*p + t) - 330*cos(11*p + 3*t) + 165*cos(11*p + 5*t) - 55*cos(11*p + 7*t) + 11*cos(11*p + 9*t) - cos(11*p + 11*t)); break;
      case -10 : return -0.001797020360744780695131126972263*(sin(10*p - 11*t) - 9*sin(10*p - 9*t) + 35*sin(10*p - 7*t) - 75*sin(10*p - 5*t) + 90*sin(10*p - 3*t) - 42*sin(10*p - t) - 42*sin(10*p + t) + 90*sin(10*p + 3*t) - 75*sin(10*p + 5*t) + 35*sin(10*p + 7*t) - 9*sin(10*p + 9*t) + sin(10*p + 11*t)); break;
      case -9 : return 0.000277286261613006547176308033139*(21*cos(9*p - 11*t) - 151*cos(9*p - 9*t) + 435*cos(9*p - 7*t) - 585*cos(9*p - 5*t) + 210*cos(9*p - 3*t) + 378*cos(9*p - t) - 378*cos(9*p + t) - 210*cos(9*p + 3*t) + 585*cos(9*p + 5*t) - 435*cos(9*p + 7*t) + 151*cos(9*p + 9*t) - 21*cos(9*p + 11*t)); break;
      case -8 : return 0.00214785014671857461377915412831*(7*sin(8*p - 11*t) - 39*sin(8*p - 9*t) + 77*sin(8*p - 7*t) - 45*sin(8*p - 5*t) - 42*sin(8*p - 3*t) + 42*sin(8*p - t) + 42*sin(8*p + t) - 42*sin(8*p + 3*t) - 45*sin(8*p + 5*t) + 77*sin(8*p + 7*t) - 39*sin(8*p + 9*t) + 7*sin(8*p + 11*t)); break;
      case -7 : return -0.000246375308826717778183359021972*(133*cos(7*p - 11*t) - 551*cos(7*p - 9*t) + 643*cos(7*p - 7*t) + 135*cos(7*p - 5*t) - 462*cos(7*p - 3*t) - 294*cos(7*p - t) + 294*cos(7*p + t) + 462*cos(7*p + 3*t) - 135*cos(7*p + 5*t) - 643*cos(7*p + 7*t) + 551*cos(7*p + 9*t) - 133*cos(7*p + 11*t)); break;
      case -6 : return -0.0001558214270239629901990670917705*(399*sin(6*p - 11*t) - 1159*sin(6*p - 9*t) + 525*sin(6*p - 7*t) + 795*sin(6*p - 5*t) + 70*sin(6*p - 3*t) - 630*sin(6*p - t) - 630*sin(6*p + t) + 70*sin(6*p + 3*t) + 795*sin(6*p + 5*t) + 525*sin(6*p + 7*t) - 1159*sin(6*p + 9*t) + 399*sin(6*p + 11*t)); break;
      case -5 : return 0.0000462858609332688343303447767626*(2261*cos(5*p - 11*t) - 4199*cos(5*p - 9*t) - 765*cos(5*p - 7*t) + 1975*cos(5*p - 5*t) + 2450*cos(5*p - 3*t) + 1050*cos(5*p - t) - 1050*cos(5*p + t) - 2450*cos(5*p + 3*t) - 1975*cos(5*p + 5*t) + 765*cos(5*p + 7*t) + 4199*cos(5*p + 9*t) - 2261*cos(5*p + 11*t)); break;
      case -4 : return 0.000489843508991797332138028007418*(323*sin(4*p - 11*t) - 323*sin(4*p - 9*t) - 255*sin(4*p - 7*t) - 65*sin(4*p - 5*t) + 110*sin(4*p - 3*t) + 210*sin(4*p - t) + 210*sin(4*p + t) + 110*sin(4*p + 3*t) - 65*sin(4*p + 5*t) - 255*sin(4*p + 7*t) - 323*sin(4*p + 9*t) + 323*sin(4*p + 11*t)); break;
      case -3 : return -0.000223581949601910057019941221859*(969*cos(3*p - 11*t) - 323*cos(3*p - 9*t) - 561*cos(3*p - 7*t) - 525*cos(3*p - 5*t) - 358*cos(3*p - 3*t) - 126*cos(3*p - t) + 126*cos(3*p + t) + 358*cos(3*p + 3*t) + 525*cos(3*p + 5*t) + 561*cos(3*p + 7*t) + 323*cos(3*p + 9*t) - 969*cos(3*p + 11*t)); break;
      case -2 : return -0.0001195095790396151212032517671255*(2261*sin(2*p - 11*t) + 323*sin(2*p - 9*t) - 289*sin(2*p - 7*t) - 615*sin(2*p - 5*t) - 798*sin(2*p - 3*t) - 882*sin(2*p - t) - 882*sin(2*p + t) - 798*sin(2*p + 3*t) - 615*sin(2*p + 5*t) - 289*sin(2*p + 7*t) + 323*sin(2*p + 9*t) + 2261*sin(2*p + 11*t)); break;
      case -1 : return 0.00001048168346807016929290242006096*(29393*cos(p - 11*t) + 12597*cos(p - 9*t) + 7735*cos(p - 7*t) + 4875*cos(p - 5*t) + 2730*cos(p - 3*t) + 882*cos(p - t) - 882*cos(p + t) - 2730*cos(p + 3*t) - 4875*cos(p + 5*t) - 7735*cos(p + 7*t) - 12597*cos(p + 9*t) - 29393*cos(p + 11*t)); break;
      case 0 : return 5.1608241842327211090243086402e-6*(29106*cos(t) + 30030*cos(3*t) + 32175*cos(5*t) + 36465*cos(7*t) + 46189*cos(9*t) + 88179*cos(11*t)); break;
      case 1 : return -0.00001048168346807016929290242006096*(29393*sin(p - 11*t) + 12597*sin(p - 9*t) + 7735*sin(p - 7*t) + 4875*sin(p - 5*t) + 2730*sin(p - 3*t) + 882*sin(p - t) - 882*sin(p + t) - 2730*sin(p + 3*t) - 4875*sin(p + 5*t) - 7735*sin(p + 7*t) - 12597*sin(p + 9*t) - 29393*sin(p + 11*t)); break;
      case 2 : return -0.0001195095790396151212032517671255*(2261*cos(2*p - 11*t) + 323*cos(2*p - 9*t) - 289*cos(2*p - 7*t) - 615*cos(2*p - 5*t) - 798*cos(2*p - 3*t) - 882*cos(2*p - t) - 882*cos(2*p + t) - 798*cos(2*p + 3*t) - 615*cos(2*p + 5*t) - 289*cos(2*p + 7*t) + 323*cos(2*p + 9*t) + 2261*cos(2*p + 11*t)); break;
      case 3 : return 0.000223581949601910057019941221859*(969*sin(3*p - 11*t) - 323*sin(3*p - 9*t) - 561*sin(3*p - 7*t) - 525*sin(3*p - 5*t) - 358*sin(3*p - 3*t) - 126*sin(3*p - t) + 126*sin(3*p + t) + 358*sin(3*p + 3*t) + 525*sin(3*p + 5*t) + 561*sin(3*p + 7*t) + 323*sin(3*p + 9*t) - 969*sin(3*p + 11*t)); break;
      case 4 : return 0.000489843508991797332138028007418*(323*cos(4*p - 11*t) - 323*cos(4*p - 9*t) - 255*cos(4*p - 7*t) - 65*cos(4*p - 5*t) + 110*cos(4*p - 3*t) + 210*cos(4*p - t) + 210*cos(4*p + t) + 110*cos(4*p + 3*t) - 65*cos(4*p + 5*t) - 255*cos(4*p + 7*t) - 323*cos(4*p + 9*t) + 323*cos(4*p + 11*t)); break;
      case 5 : return -0.0000462858609332688343303447767626*(2261*sin(5*p - 11*t) - 4199*sin(5*p - 9*t) - 765*sin(5*p - 7*t) + 1975*sin(5*p - 5*t) + 2450*sin(5*p - 3*t) + 1050*sin(5*p - t) - 1050*sin(5*p + t) - 2450*sin(5*p + 3*t) - 1975*sin(5*p + 5*t) + 765*sin(5*p + 7*t) + 4199*sin(5*p + 9*t) - 2261*sin(5*p + 11*t)); break;
      case 6 : return -0.0001558214270239629901990670917705*(399*cos(6*p - 11*t) - 1159*cos(6*p - 9*t) + 525*cos(6*p - 7*t) + 795*cos(6*p - 5*t) + 70*cos(6*p - 3*t) - 630*cos(6*p - t) - 630*cos(6*p + t) + 70*cos(6*p + 3*t) + 795*cos(6*p + 5*t) + 525*cos(6*p + 7*t) - 1159*cos(6*p + 9*t) + 399*cos(6*p + 11*t)); break;
      case 7 : return 0.000246375308826717778183359021972*(133*sin(7*p - 11*t) - 551*sin(7*p - 9*t) + 643*sin(7*p - 7*t) + 135*sin(7*p - 5*t) - 462*sin(7*p - 3*t) - 294*sin(7*p - t) + 294*sin(7*p + t) + 462*sin(7*p + 3*t) - 135*sin(7*p + 5*t) - 643*sin(7*p + 7*t) + 551*sin(7*p + 9*t) - 133*sin(7*p + 11*t)); break;
      case 8 : return 0.00214785014671857461377915412831*(7*cos(8*p - 11*t) - 39*cos(8*p - 9*t) + 77*cos(8*p - 7*t) - 45*cos(8*p - 5*t) - 42*cos(8*p - 3*t) + 42*cos(8*p - t) + 42*cos(8*p + t) - 42*cos(8*p + 3*t) - 45*cos(8*p + 5*t) + 77*cos(8*p + 7*t) - 39*cos(8*p + 9*t) + 7*cos(8*p + 11*t)); break;
      case 9 : return -0.000277286261613006547176308033139*(21*sin(9*p - 11*t) - 151*sin(9*p - 9*t) + 435*sin(9*p - 7*t) - 585*sin(9*p - 5*t) + 210*sin(9*p - 3*t) + 378*sin(9*p - t) - 378*sin(9*p + t) - 210*sin(9*p + 3*t) + 585*sin(9*p + 5*t) - 435*sin(9*p + 7*t) + 151*sin(9*p + 9*t) - 21*sin(9*p + 11*t)); break;
      case 10 : return -0.001797020360744780695131126972263*(cos(10*p - 11*t) - 9*cos(10*p - 9*t) + 35*cos(10*p - 7*t) - 75*cos(10*p - 5*t) + 90*cos(10*p - 3*t) - 42*cos(10*p - t) - 42*cos(10*p + t) + 90*cos(10*p + 3*t) - 75*cos(10*p + 5*t) + 35*cos(10*p + 7*t) - 9*cos(10*p + 9*t) + cos(10*p + 11*t)); break;
      case 11 : return 0.000383126028216404728392520692202*(sin(11*p - 11*t) - 11*sin(11*p - 9*t) + 55*sin(11*p - 7*t) - 165*sin(11*p - 5*t) + 330*sin(11*p - 3*t) - 462*sin(11*p - t) + 462*sin(11*p + t) - 330*sin(11*p + 3*t) + 165*sin(11*p + 5*t) - 55*sin(11*p + 7*t) + 11*sin(11*p + 9*t) - sin(11*p + 11*t)); break;
    }
    case 12 : switch (m) 
    {
      case -12 : return 0.0001955131825644462192637871742968*(924*sin(12*p) + sin(12*p - 12*t) - 12*sin(12*p - 10*t) + 66*sin(12*p - 8*t) - 220*sin(12*p - 6*t) + 495*sin(12*p - 4*t) - 792*sin(12*p - 2*t) - 792*sin(12*p + 2*t) + 495*sin(12*p + 4*t) - 220*sin(12*p + 6*t) + 66*sin(12*p + 8*t) - 12*sin(12*p + 10*t) + sin(12*p + 12*t)); break;
      case -11 : return -0.000957815070541011820981301730505*(cos(11*p - 12*t) - 10*cos(11*p - 10*t) + 44*cos(11*p - 8*t) - 110*cos(11*p - 6*t) + 165*cos(11*p - 4*t) - 132*cos(11*p - 2*t) + 132*cos(11*p + 2*t) - 165*cos(11*p + 4*t) + 110*cos(11*p + 6*t) - 44*cos(11*p + 8*t) + 10*cos(11*p + 10*t) - cos(11*p + 12*t)); break;
      case -10 : return 0.0001412221276352075699142080544718*(924*sin(10*p) - 23*sin(10*p - 12*t) + 188*sin(10*p - 10*t) - 638*sin(10*p - 8*t) + 1100*sin(10*p - 6*t) - 825*sin(10*p - 4*t) - 264*sin(10*p - 2*t) - 264*sin(10*p + 2*t) - 825*sin(10*p + 4*t) + 1100*sin(10*p + 6*t) - 638*sin(10*p + 8*t) + 188*sin(10*p + 10*t) - 23*sin(10*p + 12*t)); break;
      case -9 : return 0.000382431329497609225230113170196*(23*cos(9*p - 12*t) - 150*cos(9*p - 10*t) + 372*cos(9*p - 8*t) - 370*cos(9*p - 6*t) - 45*cos(9*p - 4*t) + 324*cos(9*p - 2*t) - 324*cos(9*p + 2*t) + 45*cos(9*p + 4*t) + 370*cos(9*p + 6*t) - 372*cos(9*p + 8*t) + 150*cos(9*p + 10*t) - 23*cos(9*p + 12*t)); break;
      case -8 : return 0.0001251800368246708879397256156817*(924*sin(8*p) + 161*sin(8*p - 12*t) - 812*sin(8*p - 10*t) + 1346*sin(8*p - 8*t) - 380*sin(8*p - 6*t) - 945*sin(8*p - 4*t) + 168*sin(8*p - 2*t) + 168*sin(8*p + 2*t) - 945*sin(8*p + 4*t) - 380*sin(8*p + 6*t) + 1346*sin(8*p + 8*t) - 812*sin(8*p + 10*t) + 161*sin(8*p + 12*t)); break;
      case -7 : return -0.000250360073649341775879451231363*(161*cos(7*p - 12*t) - 602*cos(7*p - 10*t) + 556*cos(7*p - 8*t) + 290*cos(7*p - 6*t) - 315*cos(7*p - 4*t) - 420*cos(7*p - 2*t) + 420*cos(7*p + 2*t) + 315*cos(7*p + 4*t) - 290*cos(7*p + 6*t) - 556*cos(7*p + 8*t) + 602*cos(7*p + 10*t) - 161*cos(7*p + 12*t)); break;
      case -6 : return 0.0000234483692766521440912488607244*(4620*sin(6*p) - 3059*sin(6*p - 12*t) + 7980*sin(6*p - 10*t) - 2166*sin(6*p - 8*t) - 5380*sin(6*p - 6*t) - 2205*sin(6*p - 4*t) + 2520*sin(6*p - 2*t) + 2520*sin(6*p + 2*t) - 2205*sin(6*p + 4*t) - 5380*sin(6*p + 6*t) - 2166*sin(6*p + 8*t) + 7980*sin(6*p + 10*t) - 3059*sin(6*p + 12*t)); break;
      case -5 : return 0.000263207292335365907098227366285*(437*cos(5*p - 12*t) - 722*cos(5*p - 10*t) - 228*cos(5*p - 8*t) + 250*cos(5*p - 6*t) + 425*cos(5*p - 4*t) + 300*cos(5*p - 2*t) - 300*cos(5*p + 2*t) - 425*cos(5*p + 4*t) - 250*cos(5*p + 6*t) + 228*cos(5*p + 8*t) + 722*cos(5*p + 10*t) - 437*cos(5*p + 12*t)); break;
      case -4 : return 0.0000225698391173529818724416546524*(4620*sin(4*p) + 7429*sin(4*p - 12*t) - 6460*sin(4*p - 10*t) - 5814*sin(4*p - 8*t) - 2380*sin(4*p - 6*t) + 1195*sin(4*p - 4*t) + 3720*sin(4*p - 2*t) + 3720*sin(4*p + 2*t) + 1195*sin(4*p + 4*t) - 2380*sin(4*p + 6*t) - 5814*sin(4*p + 8*t) - 6460*sin(4*p + 10*t) + 7429*sin(4*p + 12*t)); break;
      case -3 : return -0.0000300931188231373091632555395365*(7429*cos(3*p - 12*t) - 1938*cos(3*p - 10*t) - 3876*cos(3*p - 8*t) - 3910*cos(3*p - 6*t) - 3015*cos(3*p - 4*t) - 1620*cos(3*p - 2*t) + 1620*cos(3*p + 2*t) + 3015*cos(3*p + 4*t) + 3910*cos(3*p + 6*t) + 3876*cos(3*p + 8*t) + 1938*cos(3*p + 10*t) - 7429*cos(3*p + 12*t)); break;
      case -2 : return 0.0000368563929428151113108033083523*(2772*sin(2*p) - 7429*sin(2*p - 12*t) - 1292*sin(2*p - 10*t) + 646*sin(2*p - 8*t) + 1700*sin(2*p - 6*t) + 2325*sin(2*p - 4*t) + 2664*sin(2*p - 2*t) + 2664*sin(2*p + 2*t) + 2325*sin(2*p + 4*t) + 1700*sin(2*p + 6*t) + 646*sin(2*p + 8*t) - 1292*sin(2*p + 10*t) - 7429*sin(2*p + 12*t)); break;
      case -1 : return 5.93994556089267635062966546956e-6*(52003*cos(p - 12*t) + 22610*cos(p - 10*t) + 14212*cos(p - 8*t) + 9350*cos(p - 6*t) + 5775*cos(p - 4*t) + 2772*cos(p - 2*t) - 2772*cos(p + 2*t) - 5775*cos(p + 4*t) - 9350*cos(p + 6*t) - 14212*cos(p + 8*t) - 22610*cos(p + 10*t) - 52003*cos(p + 12*t)); break;
      case 0 : return 6.72566394266791685757731737567e-7*(106722 + 216216*cos(2*t) + 225225*cos(4*t) + 243100*cos(6*t) + 277134*cos(8*t) + 352716*cos(10*t) + 676039*cos(12*t)); break;
      case 1 : return -5.93994556089267635062966546956e-6*(52003*sin(p - 12*t) + 22610*sin(p - 10*t) + 14212*sin(p - 8*t) + 9350*sin(p - 6*t) + 5775*sin(p - 4*t) + 2772*sin(p - 2*t) - 2772*sin(p + 2*t) - 5775*sin(p + 4*t) - 9350*sin(p + 6*t) - 14212*sin(p + 8*t) - 22610*sin(p + 10*t) - 52003*sin(p + 12*t)); break;
      case 2 : return 0.0000368563929428151113108033083523*(2772*cos(2*p) - 7429*cos(2*p - 12*t) - 1292*cos(2*p - 10*t) + 646*cos(2*p - 8*t) + 1700*cos(2*p - 6*t) + 2325*cos(2*p - 4*t) + 2664*cos(2*p - 2*t) + 2664*cos(2*p + 2*t) + 2325*cos(2*p + 4*t) + 1700*cos(2*p + 6*t) + 646*cos(2*p + 8*t) - 1292*cos(2*p + 10*t) - 7429*cos(2*p + 12*t)); break;
      case 3 : return 0.0000300931188231373091632555395365*(7429*sin(3*p - 12*t) - 1938*sin(3*p - 10*t) - 3876*sin(3*p - 8*t) - 3910*sin(3*p - 6*t) - 3015*sin(3*p - 4*t) - 1620*sin(3*p - 2*t) + 1620*sin(3*p + 2*t) + 3015*sin(3*p + 4*t) + 3910*sin(3*p + 6*t) + 3876*sin(3*p + 8*t) + 1938*sin(3*p + 10*t) - 7429*sin(3*p + 12*t)); break;
      case 4 : return 0.0000225698391173529818724416546524*(4620*cos(4*p) + 7429*cos(4*p - 12*t) - 6460*cos(4*p - 10*t) - 5814*cos(4*p - 8*t) - 2380*cos(4*p - 6*t) + 1195*cos(4*p - 4*t) + 3720*cos(4*p - 2*t) + 3720*cos(4*p + 2*t) + 1195*cos(4*p + 4*t) - 2380*cos(4*p + 6*t) - 5814*cos(4*p + 8*t) - 6460*cos(4*p + 10*t) + 7429*cos(4*p + 12*t)); break;
      case 5 : return -0.000263207292335365907098227366285*(437*sin(5*p - 12*t) - 722*sin(5*p - 10*t) - 228*sin(5*p - 8*t) + 250*sin(5*p - 6*t) + 425*sin(5*p - 4*t) + 300*sin(5*p - 2*t) - 300*sin(5*p + 2*t) - 425*sin(5*p + 4*t) - 250*sin(5*p + 6*t) + 228*sin(5*p + 8*t) + 722*sin(5*p + 10*t) - 437*sin(5*p + 12*t)); break;
      case 6 : return 0.0000234483692766521440912488607244*(4620*cos(6*p) - 3059*cos(6*p - 12*t) + 7980*cos(6*p - 10*t) - 2166*cos(6*p - 8*t) - 5380*cos(6*p - 6*t) - 2205*cos(6*p - 4*t) + 2520*cos(6*p - 2*t) + 2520*cos(6*p + 2*t) - 2205*cos(6*p + 4*t) - 5380*cos(6*p + 6*t) - 2166*cos(6*p + 8*t) + 7980*cos(6*p + 10*t) - 3059*cos(6*p + 12*t)); break;
      case 7 : return 0.000250360073649341775879451231363*(161*sin(7*p - 12*t) - 602*sin(7*p - 10*t) + 556*sin(7*p - 8*t) + 290*sin(7*p - 6*t) - 315*sin(7*p - 4*t) - 420*sin(7*p - 2*t) + 420*sin(7*p + 2*t) + 315*sin(7*p + 4*t) - 290*sin(7*p + 6*t) - 556*sin(7*p + 8*t) + 602*sin(7*p + 10*t) - 161*sin(7*p + 12*t)); break;
      case 8 : return 0.0001251800368246708879397256156817*(924*cos(8*p) + 161*cos(8*p - 12*t) - 812*cos(8*p - 10*t) + 1346*cos(8*p - 8*t) - 380*cos(8*p - 6*t) - 945*cos(8*p - 4*t) + 168*cos(8*p - 2*t) + 168*cos(8*p + 2*t) - 945*cos(8*p + 4*t) - 380*cos(8*p + 6*t) + 1346*cos(8*p + 8*t) - 812*cos(8*p + 10*t) + 161*cos(8*p + 12*t)); break;
      case 9 : return -0.000382431329497609225230113170196*(23*sin(9*p - 12*t) - 150*sin(9*p - 10*t) + 372*sin(9*p - 8*t) - 370*sin(9*p - 6*t) - 45*sin(9*p - 4*t) + 324*sin(9*p - 2*t) - 324*sin(9*p + 2*t) + 45*sin(9*p + 4*t) + 370*sin(9*p + 6*t) - 372*sin(9*p + 8*t) + 150*sin(9*p + 10*t) - 23*sin(9*p + 12*t)); break;
      case 10 : return 0.0001412221276352075699142080544718*(924*cos(10*p) - 23*cos(10*p - 12*t) + 188*cos(10*p - 10*t) - 638*cos(10*p - 8*t) + 1100*cos(10*p - 6*t) - 825*cos(10*p - 4*t) - 264*cos(10*p - 2*t) - 264*cos(10*p + 2*t) - 825*cos(10*p + 4*t) + 1100*cos(10*p + 6*t) - 638*cos(10*p + 8*t) + 188*cos(10*p + 10*t) - 23*cos(10*p + 12*t)); break;
      case 11 : return 0.000957815070541011820981301730505*(sin(11*p - 12*t) - 10*sin(11*p - 10*t) + 44*sin(11*p - 8*t) - 110*sin(11*p - 6*t) + 165*sin(11*p - 4*t) - 132*sin(11*p - 2*t) + 132*sin(11*p + 2*t) - 165*sin(11*p + 4*t) + 110*sin(11*p + 6*t) - 44*sin(11*p + 8*t) + 10*sin(11*p + 10*t) - sin(11*p + 12*t)); break;
      case 12 : return 0.0001955131825644462192637871742968*(924*cos(12*p) + cos(12*p - 12*t) - 12*cos(12*p - 10*t) + 66*cos(12*p - 8*t) - 220*cos(12*p - 6*t) + 495*cos(12*p - 4*t) - 792*cos(12*p - 2*t) - 792*cos(12*p + 2*t) + 495*cos(12*p + 4*t) - 220*cos(12*p + 6*t) + 66*cos(12*p + 8*t) - 12*cos(12*p + 10*t) + cos(12*p + 12*t)); break;
    }
  }
  return 0;
}

double ReSphericalHarmonicY(const int k, const int m, const double t, const double p)
{
  switch (k)
  {
    case 0 : switch (m) 
    {
      case 0 : return 0.28209479177387814347403972578; break;
    }
    case 1 : switch (m) 
    {
      case -1 : return -0.1727470747356677396326223230159*(sin(p - t) - sin(p + t)); break;
      case 0 : return 0.488602511902919921586384622838*cos(t); break;
      case 1 : return 0.1727470747356677396326223230159*(sin(p - t) - sin(p + t)); break;
    }
    case 2 : switch (m) 
    {
      case -2 : return 0.0965685505057973950855481838279*(2*cos(2*p) - cos(2*p - 2*t) - cos(2*p + 2*t)); break;
      case -1 : return -0.1931371010115947901710963676557*(sin(p - 2*t) - sin(p + 2*t)); break;
      case 0 : return 0.1576957826262600030154468451479*(1 + 3*cos(2*t)); break;
      case 1 : return 0.1931371010115947901710963676557*(sin(p - 2*t) - sin(p + 2*t)); break;
      case 2 : return 0.0965685505057973950855481838279*(2*cos(2*p) - cos(2*p - 2*t) - cos(2*p + 2*t)); break;
    }
    case 3 : switch (m) 
    {
      case -3 : return 0.0521529779540980112155533769671*(sin(3*p - 3*t) - 3*sin(3*p - t) + 3*sin(3*p + t) - sin(3*p + 3*t)); break;
      case -2 : return -0.1277481845541602954245143114626*(cos(2*p - 3*t) - cos(2*p - t) - cos(2*p + t) + cos(2*p + 3*t)); break;
      case -1 : return -0.0403975230142688316259242704163*(5*sin(p - 3*t) + sin(p - t) - sin(p + t) - 5*sin(p + 3*t)); break;
      case 0 : return 0.0932940831475288478535989782997*(3*cos(t) + 5*cos(3*t)); break;
      case 1 : return 0.0403975230142688316259242704163*(5*sin(p - 3*t) + sin(p - t) - sin(p + t) - 5*sin(p + 3*t)); break;
      case 2 : return -0.1277481845541602954245143114626*(cos(2*p - 3*t) - cos(2*p - t) - cos(2*p + t) + cos(2*p + 3*t)); break;
      case 3 : return -0.0521529779540980112155533769671*(sin(3*p - 3*t) - 3*sin(3*p - t) + 3*sin(3*p + t) - sin(3*p + 3*t)); break;
    }
    case 4 : switch (m) 
    {
      case -4 : return 0.0276582932778114145474504817598*(6*cos(4*p) + cos(4*p - 4*t) - 4*cos(4*p - 2*t) - 4*cos(4*p + 2*t) + cos(4*p + 4*t)); break;
      case -3 : return 0.0782294669311470168233300654506*(sin(3*p - 4*t) - 2*sin(3*p - 2*t) + 2*sin(3*p + 2*t) - sin(3*p + 4*t)); break;
      case -2 : return 0.0209077044861652864985035121785*(6*cos(2*p) - 7*cos(2*p - 4*t) + 4*cos(2*p - 2*t) + 4*cos(2*p + 2*t) - 7*cos(2*p + 4*t)); break;
      case -1 : return -0.0295679592424237505653962834652*(7*sin(p - 4*t) + 2*sin(p - 2*t) - 2*sin(p + 2*t) - 7*sin(p + 4*t)); break;
      case 0 : return 0.01322319336440053797534561214596*(9 + 20*cos(2*t) + 35*cos(4*t)); break;
      case 1 : return 0.0295679592424237505653962834652*(7*sin(p - 4*t) + 2*sin(p - 2*t) - 2*sin(p + 2*t) - 7*sin(p + 4*t)); break;
      case 2 : return 0.0209077044861652864985035121785*(6*cos(2*p) - 7*cos(2*p - 4*t) + 4*cos(2*p - 2*t) + 4*cos(2*p + 2*t) - 7*cos(2*p + 4*t)); break;
      case 3 : return -0.0782294669311470168233300654506*(sin(3*p - 4*t) - 2*sin(3*p - 2*t) + 2*sin(3*p + 2*t) - sin(3*p + 4*t)); break;
      case 4 : return 0.0276582932778114145474504817598*(6*cos(4*p) + cos(4*p - 4*t) - 4*cos(4*p - 2*t) - 4*cos(4*p + 2*t) + cos(4*p + 4*t)); break;
    }
    case 5 : switch (m) 
    {
      case -5 : return -0.01450413135752681752056245642293*(sin(5*p - 5*t) - 5*sin(5*p - 3*t) + 10*sin(5*p - t) - 10*sin(5*p + t) + 5*sin(5*p + 3*t) - sin(5*p + 5*t)); break;
      case -4 : return 0.0458660905720547238453758171103*(cos(4*p - 5*t) - 3*cos(4*p - 3*t) + 2*cos(4*p - t) + 2*cos(4*p + t) - 3*cos(4*p + 3*t) + cos(4*p + 5*t)); break;
      case -3 : return 0.01081074122333875666143645063852*(9*sin(3*p - 5*t) - 13*sin(3*p - 3*t) - 6*sin(3*p - t) + 6*sin(3*p + t) + 13*sin(3*p + 3*t) - 9*sin(3*p + 5*t)); break;
      case -2 : return -0.0529615994769031023692403611097*(3*cos(2*p - 5*t) - cos(2*p - 3*t) - 2*cos(2*p - t) - 2*cos(2*p + t) - cos(2*p + 3*t) + 3*cos(2*p + 5*t)); break;
      case -1 : return -0.01000880151800672273553629477647*(21*sin(p - 5*t) + 7*sin(p - 3*t) + 2*sin(p - t) - 2*sin(p + t) - 7*sin(p + 3*t) - 21*sin(p + 5*t)); break;
      case 0 : return 0.00730939515333897477748219950564*(30*cos(t) + 35*cos(3*t) + 63*cos(5*t)); break;
      case 1 : return 0.01000880151800672273553629477647*(21*sin(p - 5*t) + 7*sin(p - 3*t) + 2*sin(p - t) - 2*sin(p + t) - 7*sin(p + 3*t) - 21*sin(p + 5*t)); break;
      case 2 : return -0.0529615994769031023692403611097*(3*cos(2*p - 5*t) - cos(2*p - 3*t) - 2*cos(2*p - t) - 2*cos(2*p + t) - cos(2*p + 3*t) + 3*cos(2*p + 5*t)); break;
      case 3 : return -0.01081074122333875666143645063852*(9*sin(3*p - 5*t) - 13*sin(3*p - 3*t) - 6*sin(3*p - t) + 6*sin(3*p + t) + 13*sin(3*p + 3*t) - 9*sin(3*p + 5*t)); break;
      case 4 : return 0.0458660905720547238453758171103*(cos(4*p - 5*t) - 3*cos(4*p - 3*t) + 2*cos(4*p - t) + 2*cos(4*p + t) - 3*cos(4*p + 3*t) + cos(4*p + 5*t)); break;
      case 5 : return 0.01450413135752681752056245642293*(sin(5*p - 5*t) - 5*sin(5*p - 3*t) + 10*sin(5*p - t) - 10*sin(5*p + t) + 5*sin(5*p + 3*t) - sin(5*p + 5*t)); break;
    }
    case 6 : switch (m) 
    {
      case -6 : return 0.00754818927468853479645251755055*(20*cos(6*p) - cos(6*p - 6*t) + 6*cos(6*p - 4*t) - 15*cos(6*p - 2*t) - 15*cos(6*p + 2*t) + 6*cos(6*p + 4*t) - cos(6*p + 6*t)); break;
      case -5 : return -0.0261476946578140297069217513434*(sin(5*p - 6*t) - 4*sin(5*p - 4*t) + 5*sin(5*p - 2*t) - 5*sin(5*p + 2*t) + 4*sin(5*p + 4*t) - sin(5*p + 6*t)); break;
      case -4 : return 0.0055747072320937191815291748398*(20*cos(4*p) + 11*cos(4*p - 6*t) - 26*cos(4*p - 4*t) + 5*cos(4*p - 2*t) + 5*cos(4*p + 2*t) - 26*cos(4*p + 4*t) + 11*cos(4*p + 6*t)); break;
      case -3 : return 0.01017797634168305839948117512343*(11*sin(3*p - 6*t) - 12*sin(3*p - 4*t) - 9*sin(3*p - 2*t) + 9*sin(3*p + 2*t) + 12*sin(3*p + 4*t) - 11*sin(3*p + 6*t)); break;
      case -2 : return 0.00508898817084152919974058756172*(20*cos(2*p) - 33*cos(2*p - 6*t) + 6*cos(2*p - 4*t) + 17*cos(2*p - 2*t) + 17*cos(2*p + 2*t) + 6*cos(2*p + 4*t) - 33*cos(2*p + 6*t)); break;
      case -1 : return -0.00643711744220532464693356226676*(33*sin(p - 6*t) + 12*sin(p - 4*t) + 5*sin(p - 2*t) - 5*sin(p + 2*t) - 12*sin(p + 4*t) - 33*sin(p + 6*t)); break;
      case 0 : return 0.001986537570863388310400844720642*(50 + 105*cos(2*t) + 126*cos(4*t) + 231*cos(6*t)); break;
      case 1 : return 0.00643711744220532464693356226676*(33*sin(p - 6*t) + 12*sin(p - 4*t) + 5*sin(p - 2*t) - 5*sin(p + 2*t) - 12*sin(p + 4*t) - 33*sin(p + 6*t)); break;
      case 2 : return 0.00508898817084152919974058756172*(20*cos(2*p) - 33*cos(2*p - 6*t) + 6*cos(2*p - 4*t) + 17*cos(2*p - 2*t) + 17*cos(2*p + 2*t) + 6*cos(2*p + 4*t) - 33*cos(2*p + 6*t)); break;
      case 3 : return -0.01017797634168305839948117512343*(11*sin(3*p - 6*t) - 12*sin(3*p - 4*t) - 9*sin(3*p - 2*t) + 9*sin(3*p + 2*t) + 12*sin(3*p + 4*t) - 11*sin(3*p + 6*t)); break;
      case 4 : return 0.0055747072320937191815291748398*(20*cos(4*p) + 11*cos(4*p - 6*t) - 26*cos(4*p - 4*t) + 5*cos(4*p - 2*t) + 5*cos(4*p + 2*t) - 26*cos(4*p + 4*t) + 11*cos(4*p + 6*t)); break;
      case 5 : return 0.0261476946578140297069217513434*(sin(5*p - 6*t) - 4*sin(5*p - 4*t) + 5*sin(5*p - 2*t) - 5*sin(5*p + 2*t) + 4*sin(5*p + 4*t) - sin(5*p + 6*t)); break;
      case 6 : return 0.00754818927468853479645251755055*(20*cos(6*p) - cos(6*p - 6*t) + 6*cos(6*p - 4*t) - 15*cos(6*p - 2*t) - 15*cos(6*p + 2*t) + 6*cos(6*p + 4*t) - cos(6*p + 6*t)); break;
    }
    case 7 : switch (m) 
    {
      case -7 : return 0.00390655909039492706258854893785*(sin(7*p - 7*t) - 7*sin(7*p - 5*t) + 21*sin(7*p - 3*t) - 35*sin(7*p - t) + 35*sin(7*p + t) - 21*sin(7*p + 3*t) + 7*sin(7*p + 5*t) - sin(7*p + 7*t)); break;
      case -6 : return -0.01461700567744506825595076737686*(cos(6*p - 7*t) - 5*cos(6*p - 5*t) + 9*cos(6*p - 3*t) - 5*cos(6*p - t) - 5*cos(6*p + t) + 9*cos(6*p + 3*t) - 5*cos(6*p + 5*t) + cos(6*p + 7*t)); break;
      case -5 : return -0.00286663066075342024033598856939*(13*sin(5*p - 7*t) - 43*sin(5*p - 5*t) + 33*sin(5*p - 3*t) + 25*sin(5*p - t) - 25*sin(5*p + t) - 33*sin(5*p + 3*t) + 43*sin(5*p + 5*t) - 13*sin(5*p + 7*t)); break;
      case -4 : return 0.00573326132150684048067197713878*(13*cos(4*p - 7*t) - 25*cos(4*p - 5*t) - 3*cos(4*p - 3*t) + 15*cos(4*p - t) + 15*cos(4*p + t) - 3*cos(4*p + 3*t) - 25*cos(4*p + 5*t) + 13*cos(4*p + 7*t)); break;
      case -3 : return 0.000864321664931606704607827554992*(143*sin(3*p - 7*t) - 121*sin(3*p - 5*t) - 117*sin(3*p - 3*t) - 45*sin(3*p - t) + 45*sin(3*p + t) + 117*sin(3*p + 3*t) + 121*sin(3*p + 5*t) - 143*sin(3*p + 7*t)); break;
      case -2 : return -0.001222335420799172137864532272666*(143*cos(2*p - 7*t) - 11*cos(2*p - 5*t) - 57*cos(2*p - 3*t) - 75*cos(2*p - t) - 75*cos(2*p + t) - 57*cos(2*p + 3*t) - 11*cos(2*p + 5*t) + 143*cos(2*p + 7*t)); break;
      case -1 : return -0.000499016345914688654002009029151*(429*sin(p - 7*t) + 165*sin(p - 5*t) + 81*sin(p - 3*t) + 25*sin(p - t) - 25*sin(p + t) - 81*sin(p + 3*t) - 165*sin(p + 5*t) - 429*sin(p + 7*t)); break;
      case 0 : return 0.001066941826750077217327525103323*(175*cos(t) + 189*cos(3*t) + 231*cos(5*t) + 429*cos(7*t)); break;
      case 1 : return 0.000499016345914688654002009029151*(429*sin(p - 7*t) + 165*sin(p - 5*t) + 81*sin(p - 3*t) + 25*sin(p - t) - 25*sin(p + t) - 81*sin(p + 3*t) - 165*sin(p + 5*t) - 429*sin(p + 7*t)); break;
      case 2 : return -0.001222335420799172137864532272666*(143*cos(2*p - 7*t) - 11*cos(2*p - 5*t) - 57*cos(2*p - 3*t) - 75*cos(2*p - t) - 75*cos(2*p + t) - 57*cos(2*p + 3*t) - 11*cos(2*p + 5*t) + 143*cos(2*p + 7*t)); break;
      case 3 : return -0.000864321664931606704607827554992*(143*sin(3*p - 7*t) - 121*sin(3*p - 5*t) - 117*sin(3*p - 3*t) - 45*sin(3*p - t) + 45*sin(3*p + t) + 117*sin(3*p + 3*t) + 121*sin(3*p + 5*t) - 143*sin(3*p + 7*t)); break;
      case 4 : return 0.00573326132150684048067197713878*(13*cos(4*p - 7*t) - 25*cos(4*p - 5*t) - 3*cos(4*p - 3*t) + 15*cos(4*p - t) + 15*cos(4*p + t) - 3*cos(4*p + 3*t) - 25*cos(4*p + 5*t) + 13*cos(4*p + 7*t)); break;
      case 5 : return 0.00286663066075342024033598856939*(13*sin(5*p - 7*t) - 43*sin(5*p - 5*t) + 33*sin(5*p - 3*t) + 25*sin(5*p - t) - 25*sin(5*p + t) - 33*sin(5*p + 3*t) + 43*sin(5*p + 5*t) - 13*sin(5*p + 7*t)); break;
      case 6 : return -0.01461700567744506825595076737686*(cos(6*p - 7*t) - 5*cos(6*p - 5*t) + 9*cos(6*p - 3*t) - 5*cos(6*p - t) - 5*cos(6*p + t) + 9*cos(6*p + 3*t) - 5*cos(6*p + 5*t) + cos(6*p + 7*t)); break;
      case 7 : return -0.00390655909039492706258854893785*(sin(7*p - 7*t) - 7*sin(7*p - 5*t) + 21*sin(7*p - 3*t) - 35*sin(7*p - t) + 35*sin(7*p + t) - 21*sin(7*p + 3*t) + 7*sin(7*p + 5*t) - sin(7*p + 7*t)); break;
    }
    case 8 : switch (m) 
    {
      case -8 : return 0.00201339447030189183486334132529*(70*cos(8*p) + cos(8*p - 8*t) - 8*cos(8*p - 6*t) + 28*cos(8*p - 4*t) - 56*cos(8*p - 2*t) - 56*cos(8*p + 2*t) + 28*cos(8*p + 4*t) - 8*cos(8*p + 6*t) + cos(8*p + 8*t)); break;
      case -7 : return 0.00805357788120756733945336530116*(sin(7*p - 8*t) - 6*sin(7*p - 6*t) + 14*sin(7*p - 4*t) - 14*sin(7*p - 2*t) + 14*sin(7*p + 2*t) - 14*sin(7*p + 4*t) + 6*sin(7*p + 6*t) - sin(7*p + 8*t)); break;
      case -6 : return 0.00147037542472068188952924283307*(70*cos(6*p) - 15*cos(6*p - 8*t) + 64*cos(6*p - 6*t) - 84*cos(6*p - 4*t) - 84*cos(6*p + 4*t) + 64*cos(6*p + 6*t) - 15*cos(6*p + 8*t)); break;
      case -5 : return -0.00952912185692606606411696386973*(5*sin(5*p - 8*t) - 14*sin(5*p - 6*t) + 6*sin(5*p - 4*t) + 10*sin(5*p - 2*t) - 10*sin(5*p + 2*t) - 6*sin(5*p + 4*t) + 14*sin(5*p + 6*t) - 5*sin(5*p + 8*t)); break;
      case -4 : return 0.001321451440972752138178969129294*(70*cos(4*p) + 65*cos(4*p - 8*t) - 104*cos(4*p - 6*t) - 36*cos(4*p - 4*t) + 40*cos(4*p - 2*t) + 40*cos(4*p + 2*t) - 36*cos(4*p + 4*t) - 104*cos(4*p + 6*t) + 65*cos(4*p + 8*t)); break;
      case -3 : return 0.00341197294913950827517017222963*(39*sin(3*p - 8*t) - 26*sin(3*p - 6*t) - 30*sin(3*p - 4*t) - 18*sin(3*p - 2*t) + 18*sin(3*p + 2*t) + 30*sin(3*p + 4*t) + 26*sin(3*p + 6*t) - 39*sin(3*p + 8*t)); break;
      case -2 : return 0.001259954512472199251262382386198*(70*cos(2*p) - 143*cos(2*p - 8*t) + 44*cos(2*p - 4*t) + 64*cos(2*p - 2*t) + 64*cos(2*p + 2*t) + 44*cos(2*p + 4*t) - 143*cos(2*p + 8*t)); break;
      case -1 : return -0.000301186735953348127575578792952*(715*sin(p - 8*t) + 286*sin(p - 6*t) + 154*sin(p - 4*t) + 70*sin(p - 2*t) - 70*sin(p + 2*t) - 154*sin(p + 4*t) - 286*sin(p + 6*t) - 715*sin(p + 8*t)); break;
      case 0 : return 0.0000709903944653515335456042440843*(1225 + 2520*cos(2*t) + 2772*cos(4*t) + 3432*cos(6*t) + 6435*cos(8*t)); break;
      case 1 : return 0.000301186735953348127575578792952*(715*sin(p - 8*t) + 286*sin(p - 6*t) + 154*sin(p - 4*t) + 70*sin(p - 2*t) - 70*sin(p + 2*t) - 154*sin(p + 4*t) - 286*sin(p + 6*t) - 715*sin(p + 8*t)); break;
      case 2 : return 0.001259954512472199251262382386198*(70*cos(2*p) - 143*cos(2*p - 8*t) + 44*cos(2*p - 4*t) + 64*cos(2*p - 2*t) + 64*cos(2*p + 2*t) + 44*cos(2*p + 4*t) - 143*cos(2*p + 8*t)); break;
      case 3 : return -0.00341197294913950827517017222963*(39*sin(3*p - 8*t) - 26*sin(3*p - 6*t) - 30*sin(3*p - 4*t) - 18*sin(3*p - 2*t) + 18*sin(3*p + 2*t) + 30*sin(3*p + 4*t) + 26*sin(3*p + 6*t) - 39*sin(3*p + 8*t)); break;
      case 4 : return 0.001321451440972752138178969129294*(70*cos(4*p) + 65*cos(4*p - 8*t) - 104*cos(4*p - 6*t) - 36*cos(4*p - 4*t) + 40*cos(4*p - 2*t) + 40*cos(4*p + 2*t) - 36*cos(4*p + 4*t) - 104*cos(4*p + 6*t) + 65*cos(4*p + 8*t)); break;
      case 5 : return 0.00952912185692606606411696386973*(5*sin(5*p - 8*t) - 14*sin(5*p - 6*t) + 6*sin(5*p - 4*t) + 10*sin(5*p - 2*t) - 10*sin(5*p + 2*t) - 6*sin(5*p + 4*t) + 14*sin(5*p + 6*t) - 5*sin(5*p + 8*t)); break;
      case 6 : return 0.00147037542472068188952924283307*(70*cos(6*p) - 15*cos(6*p - 8*t) + 64*cos(6*p - 6*t) - 84*cos(6*p - 4*t) - 84*cos(6*p + 4*t) + 64*cos(6*p + 6*t) - 15*cos(6*p + 8*t)); break;
      case 7 : return -0.00805357788120756733945336530116*(sin(7*p - 8*t) - 6*sin(7*p - 6*t) + 14*sin(7*p - 4*t) - 14*sin(7*p - 2*t) + 14*sin(7*p + 2*t) - 14*sin(7*p + 4*t) + 6*sin(7*p + 6*t) - sin(7*p + 8*t)); break;
      case 8 : return 0.00201339447030189183486334132529*(70*cos(8*p) + cos(8*p - 8*t) - 8*cos(8*p - 6*t) + 28*cos(8*p - 4*t) - 56*cos(8*p - 2*t) - 56*cos(8*p + 2*t) + 28*cos(8*p + 4*t) - 8*cos(8*p + 6*t) + cos(8*p + 8*t)); break;
    }
    case 9 : switch (m) 
    {
      case -9 : return -0.001034283088852440585058962213157*(sin(9*p - 9*t) - 9*sin(9*p - 7*t) + 36*sin(9*p - 5*t) - 84*sin(9*p - 3*t) + 126*sin(9*p - t) - 126*sin(9*p + t) + 84*sin(9*p + 3*t) - 36*sin(9*p + 5*t) + 9*sin(9*p + 7*t) - sin(9*p + 9*t)); break;
      case -8 : return 0.0043880915147647751750992848194*(cos(8*p - 9*t) - 7*cos(8*p - 7*t) + 20*cos(8*p - 5*t) - 28*cos(8*p - 3*t) + 14*cos(8*p - t) + 14*cos(8*p + t) - 28*cos(8*p + 3*t) + 20*cos(8*p + 5*t) - 7*cos(8*p + 7*t) + cos(8*p + 9*t)); break;
      case -7 : return 0.000752551486258007370452908696331*(17*sin(7*p - 9*t) - 89*sin(7*p - 7*t) + 164*sin(7*p - 5*t) - 84*sin(7*p - 3*t) - 98*sin(7*p - t) + 98*sin(7*p + t) + 84*sin(7*p + 3*t) - 164*sin(7*p + 5*t) + 89*sin(7*p + 7*t) - 17*sin(7*p + 9*t)); break;
      case -6 : return -0.001737943212680454038596561231567*(17*cos(6*p - 9*t) - 63*cos(6*p - 7*t) + 60*cos(6*p - 5*t) + 28*cos(6*p - 3*t) - 42*cos(6*p - t) - 42*cos(6*p + t) + 28*cos(6*p + 3*t) + 60*cos(6*p + 5*t) - 63*cos(6*p + 7*t) + 17*cos(6*p + 9*t)); break;
      case -5 : return -0.000673102511936561327899554822462*(85*sin(5*p - 9*t) - 205*sin(5*p - 7*t) + 36*sin(5*p - 5*t) + 140*sin(5*p - 3*t) + 70*sin(5*p - t) - 70*sin(5*p + t) - 140*sin(5*p + 3*t) - 36*sin(5*p + 5*t) + 205*sin(5*p + 7*t) - 85*sin(5*p + 9*t)); break;
      case -4 : return 0.00563157965496996303859039131309*(17*cos(4*p - 9*t) - 23*cos(4*p - 7*t) - 12*cos(4*p - 5*t) + 4*cos(4*p - 3*t) + 14*cos(4*p - t) + 14*cos(4*p + t) + 4*cos(4*p + 3*t) - 12*cos(4*p + 5*t) - 23*cos(4*p + 7*t) + 17*cos(4*p + 9*t)); break;
      case -3 : return 0.000637650830927843559575950527255*(221*sin(3*p - 9*t) - 117*sin(3*p - 7*t) - 156*sin(3*p - 5*t) - 116*sin(3*p - 3*t) - 42*sin(3*p - t) + 42*sin(3*p + t) + 116*sin(3*p + 3*t) + 156*sin(3*p + 5*t) + 117*sin(3*p + 7*t) - 221*sin(3*p + 9*t)); break;
      case -2 : return -0.00083488091419380901096632445383*(221*cos(2*p - 9*t) + 13*cos(2*p - 7*t) - 52*cos(2*p - 5*t) - 84*cos(2*p - 3*t) - 98*cos(2*p - t) - 98*cos(2*p + t) - 84*cos(2*p + 3*t) - 52*cos(2*p + 5*t) + 13*cos(2*p + 7*t) + 221*cos(2*p + 9*t)); break;
      case -1 : return -0.000088998604488873504194604687416*(2431*sin(p - 9*t) + 1001*sin(p - 7*t) + 572*sin(p - 5*t) + 308*sin(p - 3*t) + 98*sin(p - t) - 98*sin(p + t) - 308*sin(p + 3*t) - 572*sin(p + 5*t) - 1001*sin(p + 7*t) - 2431*sin(p + 9*t)); break;
      case 0 : return 0.0000375251065015101235270913259486*(4410*cos(t) + 4620*cos(3*t) + 5148*cos(5*t) + 6435*cos(7*t) + 12155*cos(9*t)); break;
      case 1 : return 0.000088998604488873504194604687416*(2431*sin(p - 9*t) + 1001*sin(p - 7*t) + 572*sin(p - 5*t) + 308*sin(p - 3*t) + 98*sin(p - t) - 98*sin(p + t) - 308*sin(p + 3*t) - 572*sin(p + 5*t) - 1001*sin(p + 7*t) - 2431*sin(p + 9*t)); break;
      case 2 : return -0.00083488091419380901096632445383*(221*cos(2*p - 9*t) + 13*cos(2*p - 7*t) - 52*cos(2*p - 5*t) - 84*cos(2*p - 3*t) - 98*cos(2*p - t) - 98*cos(2*p + t) - 84*cos(2*p + 3*t) - 52*cos(2*p + 5*t) + 13*cos(2*p + 7*t) + 221*cos(2*p + 9*t)); break;
      case 3 : return -0.000637650830927843559575950527255*(221*sin(3*p - 9*t) - 117*sin(3*p - 7*t) - 156*sin(3*p - 5*t) - 116*sin(3*p - 3*t) - 42*sin(3*p - t) + 42*sin(3*p + t) + 116*sin(3*p + 3*t) + 156*sin(3*p + 5*t) + 117*sin(3*p + 7*t) - 221*sin(3*p + 9*t)); break;
      case 4 : return 0.00563157965496996303859039131309*(17*cos(4*p - 9*t) - 23*cos(4*p - 7*t) - 12*cos(4*p - 5*t) + 4*cos(4*p - 3*t) + 14*cos(4*p - t) + 14*cos(4*p + t) + 4*cos(4*p + 3*t) - 12*cos(4*p + 5*t) - 23*cos(4*p + 7*t) + 17*cos(4*p + 9*t)); break;
      case 5 : return 0.000673102511936561327899554822462*(85*sin(5*p - 9*t) - 205*sin(5*p - 7*t) + 36*sin(5*p - 5*t) + 140*sin(5*p - 3*t) + 70*sin(5*p - t) - 70*sin(5*p + t) - 140*sin(5*p + 3*t) - 36*sin(5*p + 5*t) + 205*sin(5*p + 7*t) - 85*sin(5*p + 9*t)); break;
      case 6 : return -0.001737943212680454038596561231567*(17*cos(6*p - 9*t) - 63*cos(6*p - 7*t) + 60*cos(6*p - 5*t) + 28*cos(6*p - 3*t) - 42*cos(6*p - t) - 42*cos(6*p + t) + 28*cos(6*p + 3*t) + 60*cos(6*p + 5*t) - 63*cos(6*p + 7*t) + 17*cos(6*p + 9*t)); break;
      case 7 : return -0.000752551486258007370452908696331*(17*sin(7*p - 9*t) - 89*sin(7*p - 7*t) + 164*sin(7*p - 5*t) - 84*sin(7*p - 3*t) - 98*sin(7*p - t) + 98*sin(7*p + t) + 84*sin(7*p + 3*t) - 164*sin(7*p + 5*t) + 89*sin(7*p + 7*t) - 17*sin(7*p + 9*t)); break;
      case 8 : return 0.0043880915147647751750992848194*(cos(8*p - 9*t) - 7*cos(8*p - 7*t) + 20*cos(8*p - 5*t) - 28*cos(8*p - 3*t) + 14*cos(8*p - t) + 14*cos(8*p + t) - 28*cos(8*p + 3*t) + 20*cos(8*p + 5*t) - 7*cos(8*p + 7*t) + cos(8*p + 9*t)); break;
      case 9 : return 0.001034283088852440585058962213157*(sin(9*p - 9*t) - 9*sin(9*p - 7*t) + 36*sin(9*p - 5*t) - 84*sin(9*p - 3*t) + 126*sin(9*p - t) - 126*sin(9*p + t) + 84*sin(9*p + 3*t) - 36*sin(9*p + 5*t) + 9*sin(9*p + 7*t) - sin(9*p + 9*t)); break;
    }
    case 10 : switch (m) 
    {
      case -10 : return 0.000529912394476778770252132215025*(252*cos(10*p) - cos(10*p - 10*t) + 10*cos(10*p - 8*t) - 45*cos(10*p - 6*t) + 120*cos(10*p - 4*t) - 210*cos(10*p - 2*t) - 210*cos(10*p + 2*t) + 120*cos(10*p + 4*t) - 45*cos(10*p + 6*t) + 10*cos(10*p + 8*t) - cos(10*p + 10*t)); break;
      case -9 : return -0.00236984027233952286606648120237*(sin(9*p - 10*t) - 8*sin(9*p - 8*t) + 27*sin(9*p - 6*t) - 48*sin(9*p - 4*t) + 42*sin(9*p - 2*t) - 42*sin(9*p + 2*t) + 48*sin(9*p + 4*t) - 27*sin(9*p + 6*t) + 8*sin(9*p + 8*t) - sin(9*p + 10*t)); break;
      case -8 : return 0.000384438856831830722850812827095*(252*cos(8*p) + 19*cos(8*p - 10*t) - 118*cos(8*p - 8*t) + 279*cos(8*p - 6*t) - 264*cos(8*p - 4*t) - 42*cos(8*p - 2*t) - 42*cos(8*p + 2*t) - 264*cos(8*p + 4*t) + 279*cos(8*p + 6*t) - 118*cos(8*p + 8*t) + 19*cos(8*p + 10*t)); break;
      case -7 : return 0.000941679036536860067470202171862*(19*sin(7*p - 10*t) - 88*sin(7*p - 8*t) + 129*sin(7*p - 6*t) - 16*sin(7*p - 4*t) - 98*sin(7*p - 2*t) + 98*sin(7*p + 2*t) + 16*sin(7*p + 4*t) - 129*sin(7*p + 6*t) + 88*sin(7*p + 8*t) - 19*sin(7*p + 10*t)); break;
      case -6 : return 0.0001141953568550395957496153701317*(756*cos(6*p) - 323*cos(6*p - 10*t) + 1054*cos(6*p - 8*t) - 711*cos(6*p - 6*t) - 664*cos(6*p - 4*t) + 266*cos(6*p - 2*t) + 266*cos(6*p + 2*t) - 664*cos(6*p + 4*t) - 711*cos(6*p + 6*t) + 1054*cos(6*p + 8*t) - 323*cos(6*p + 10*t)); break;
      case -5 : return -0.00020427886451417210708532618237*(323*sin(5*p - 10*t) - 680*sin(5*p - 8*t) - 15*sin(5*p - 6*t) + 400*sin(5*p - 4*t) + 350*sin(5*p - 2*t) - 350*sin(5*p + 2*t) - 400*sin(5*p + 4*t) + 15*sin(5*p + 6*t) + 680*sin(5*p + 8*t) - 323*sin(5*p + 10*t)); break;
      case -4 : return 0.00032299324484886477318972173938*(252*cos(4*p) + 323*cos(4*p - 10*t) - 374*cos(4*p - 8*t) - 249*cos(4*p - 6*t) - 8*cos(4*p - 4*t) + 182*cos(4*p - 2*t) + 182*cos(4*p + 2*t) - 8*cos(4*p + 4*t) - 249*cos(4*p + 6*t) - 374*cos(4*p + 8*t) + 323*cos(4*p + 10*t)); break;
      case -3 : return 0.000456781427420158382998461480527*(323*sin(3*p - 10*t) - 136*sin(3*p - 8*t) - 207*sin(3*p - 6*t) - 176*sin(3*p - 4*t) - 98*sin(3*p - 2*t) + 98*sin(3*p + 2*t) + 176*sin(3*p + 4*t) + 207*sin(3*p + 6*t) + 136*sin(3*p + 8*t) - 323*sin(3*p + 10*t)); break;
      case -2 : return 0.0000447911040742721912552498213968*(1764*cos(2*p) - 4199*cos(2*p - 10*t) - 442*cos(2*p - 8*t) + 741*cos(2*p - 6*t) + 1352*cos(2*p - 4*t) + 1666*cos(2*p - 2*t) + 1666*cos(2*p + 2*t) + 1352*cos(2*p + 4*t) + 741*cos(2*p + 6*t) - 442*cos(2*p + 8*t) - 4199*cos(2*p + 10*t)); break;
      case -1 : return -0.0000517203119891631859074612901462*(4199*sin(p - 10*t) + 1768*sin(p - 8*t) + 1053*sin(p - 6*t) + 624*sin(p - 4*t) + 294*sin(p - 2*t) - 294*sin(p + 2*t) - 624*sin(p + 4*t) - 1053*sin(p + 6*t) - 1768*sin(p + 8*t) - 4199*sin(p + 10*t)); break;
      case 0 : return 9.8626765171554764711729175821e-6*(7938 + 16170*cos(2*t) + 17160*cos(4*t) + 19305*cos(6*t) + 24310*cos(8*t) + 46189*cos(10*t)); break;
      case 1 : return 0.0000517203119891631859074612901462*(4199*sin(p - 10*t) + 1768*sin(p - 8*t) + 1053*sin(p - 6*t) + 624*sin(p - 4*t) + 294*sin(p - 2*t) - 294*sin(p + 2*t) - 624*sin(p + 4*t) - 1053*sin(p + 6*t) - 1768*sin(p + 8*t) - 4199*sin(p + 10*t)); break;
      case 2 : return 0.0000447911040742721912552498213968*(1764*cos(2*p) - 4199*cos(2*p - 10*t) - 442*cos(2*p - 8*t) + 741*cos(2*p - 6*t) + 1352*cos(2*p - 4*t) + 1666*cos(2*p - 2*t) + 1666*cos(2*p + 2*t) + 1352*cos(2*p + 4*t) + 741*cos(2*p + 6*t) - 442*cos(2*p + 8*t) - 4199*cos(2*p + 10*t)); break;
      case 3 : return -0.000456781427420158382998461480527*(323*sin(3*p - 10*t) - 136*sin(3*p - 8*t) - 207*sin(3*p - 6*t) - 176*sin(3*p - 4*t) - 98*sin(3*p - 2*t) + 98*sin(3*p + 2*t) + 176*sin(3*p + 4*t) + 207*sin(3*p + 6*t) + 136*sin(3*p + 8*t) - 323*sin(3*p + 10*t)); break;
      case 4 : return 0.00032299324484886477318972173938*(252*cos(4*p) + 323*cos(4*p - 10*t) - 374*cos(4*p - 8*t) - 249*cos(4*p - 6*t) - 8*cos(4*p - 4*t) + 182*cos(4*p - 2*t) + 182*cos(4*p + 2*t) - 8*cos(4*p + 4*t) - 249*cos(4*p + 6*t) - 374*cos(4*p + 8*t) + 323*cos(4*p + 10*t)); break;
      case 5 : return 0.00020427886451417210708532618237*(323*sin(5*p - 10*t) - 680*sin(5*p - 8*t) - 15*sin(5*p - 6*t) + 400*sin(5*p - 4*t) + 350*sin(5*p - 2*t) - 350*sin(5*p + 2*t) - 400*sin(5*p + 4*t) + 15*sin(5*p + 6*t) + 680*sin(5*p + 8*t) - 323*sin(5*p + 10*t)); break;
      case 6 : return 0.0001141953568550395957496153701317*(756*cos(6*p) - 323*cos(6*p - 10*t) + 1054*cos(6*p - 8*t) - 711*cos(6*p - 6*t) - 664*cos(6*p - 4*t) + 266*cos(6*p - 2*t) + 266*cos(6*p + 2*t) - 664*cos(6*p + 4*t) - 711*cos(6*p + 6*t) + 1054*cos(6*p + 8*t) - 323*cos(6*p + 10*t)); break;
      case 7 : return -0.000941679036536860067470202171862*(19*sin(7*p - 10*t) - 88*sin(7*p - 8*t) + 129*sin(7*p - 6*t) - 16*sin(7*p - 4*t) - 98*sin(7*p - 2*t) + 98*sin(7*p + 2*t) + 16*sin(7*p + 4*t) - 129*sin(7*p + 6*t) + 88*sin(7*p + 8*t) - 19*sin(7*p + 10*t)); break;
      case 8 : return 0.000384438856831830722850812827095*(252*cos(8*p) + 19*cos(8*p - 10*t) - 118*cos(8*p - 8*t) + 279*cos(8*p - 6*t) - 264*cos(8*p - 4*t) - 42*cos(8*p - 2*t) - 42*cos(8*p + 2*t) - 264*cos(8*p + 4*t) + 279*cos(8*p + 6*t) - 118*cos(8*p + 8*t) + 19*cos(8*p + 10*t)); break;
      case 9 : return 0.00236984027233952286606648120237*(sin(9*p - 10*t) - 8*sin(9*p - 8*t) + 27*sin(9*p - 6*t) - 48*sin(9*p - 4*t) + 42*sin(9*p - 2*t) - 42*sin(9*p + 2*t) + 48*sin(9*p + 4*t) - 27*sin(9*p + 6*t) + 8*sin(9*p + 8*t) - sin(9*p + 10*t)); break;
      case 10 : return 0.000529912394476778770252132215025*(252*cos(10*p) - cos(10*p - 10*t) + 10*cos(10*p - 8*t) - 45*cos(10*p - 6*t) + 120*cos(10*p - 4*t) - 210*cos(10*p - 2*t) - 210*cos(10*p + 2*t) + 120*cos(10*p + 4*t) - 45*cos(10*p + 6*t) + 10*cos(10*p + 8*t) - cos(10*p + 10*t)); break;
    }
    case 11 : switch (m) 
    {
      case -11 : return 0.000270911012600888330984112957788*(sin(11*p - 11*t) - 11*sin(11*p - 9*t) + 55*sin(11*p - 7*t) - 165*sin(11*p - 5*t) + 330*sin(11*p - 3*t) - 462*sin(11*p - t) + 462*sin(11*p + t) - 330*sin(11*p + 3*t) + 165*sin(11*p + 5*t) - 55*sin(11*p + 7*t) + 11*sin(11*p + 9*t) - sin(11*p + 11*t)); break;
      case -10 : return -0.001270685283012930339479964681051*(cos(10*p - 11*t) - 9*cos(10*p - 9*t) + 35*cos(10*p - 7*t) - 75*cos(10*p - 5*t) + 90*cos(10*p - 3*t) - 42*cos(10*p - t) - 42*cos(10*p + t) + 90*cos(10*p + 3*t) - 75*cos(10*p + 5*t) + 35*cos(10*p + 7*t) - 9*cos(10*p + 9*t) + cos(10*p + 11*t)); break;
      case -9 : return -0.0001960709959164239929600350650055*(21*sin(9*p - 11*t) - 151*sin(9*p - 9*t) + 435*sin(9*p - 7*t) - 585*sin(9*p - 5*t) + 210*sin(9*p - 3*t) + 378*sin(9*p - t) - 378*sin(9*p + t) - 210*sin(9*p + 3*t) + 585*sin(9*p + 5*t) - 435*sin(9*p + 7*t) + 151*sin(9*p + 9*t) - 21*sin(9*p + 11*t)); break;
      case -8 : return 0.00151875940371722513571203636154*(7*cos(8*p - 11*t) - 39*cos(8*p - 9*t) + 77*cos(8*p - 7*t) - 45*cos(8*p - 5*t) - 42*cos(8*p - 3*t) + 42*cos(8*p - t) + 42*cos(8*p + t) - 42*cos(8*p + 3*t) - 45*cos(8*p + 5*t) + 77*cos(8*p + 7*t) - 39*cos(8*p + 9*t) + 7*cos(8*p + 11*t)); break;
      case -7 : return 0.0001742136515883019988618247258659*(133*sin(7*p - 11*t) - 551*sin(7*p - 9*t) + 643*sin(7*p - 7*t) + 135*sin(7*p - 5*t) - 462*sin(7*p - 3*t) - 294*sin(7*p - t) + 294*sin(7*p + t) + 462*sin(7*p + 3*t) - 135*sin(7*p + 5*t) - 643*sin(7*p + 7*t) + 551*sin(7*p + 9*t) - 133*sin(7*p + 11*t)); break;
      case -6 : return -0.0001101823877028089813227240054718*(399*cos(6*p - 11*t) - 1159*cos(6*p - 9*t) + 525*cos(6*p - 7*t) + 795*cos(6*p - 5*t) + 70*cos(6*p - 3*t) - 630*cos(6*p - t) - 630*cos(6*p + t) + 70*cos(6*p + 3*t) + 795*cos(6*p + 5*t) + 525*cos(6*p + 7*t) - 1159*cos(6*p + 9*t) + 399*cos(6*p + 11*t)); break;
      case -5 : return -0.0000327290461389718940228154423877*(2261*sin(5*p - 11*t) - 4199*sin(5*p - 9*t) - 765*sin(5*p - 7*t) + 1975*sin(5*p - 5*t) + 2450*sin(5*p - 3*t) + 1050*sin(5*p - t) - 1050*sin(5*p + t) - 2450*sin(5*p + 3*t) - 1975*sin(5*p + 5*t) + 765*sin(5*p + 7*t) + 4199*sin(5*p + 9*t) - 2261*sin(5*p + 11*t)); break;
      case -4 : return 0.000346371666928313460878751369526*(323*cos(4*p - 11*t) - 323*cos(4*p - 9*t) - 255*cos(4*p - 7*t) - 65*cos(4*p - 5*t) + 110*cos(4*p - 3*t) + 210*cos(4*p - t) + 210*cos(4*p + t) + 110*cos(4*p + 3*t) - 65*cos(4*p + 5*t) - 255*cos(4*p + 7*t) - 323*cos(4*p + 9*t) + 323*cos(4*p + 11*t)); break;
      case -3 : return 0.000158096312714419511069847760972*(969*sin(3*p - 11*t) - 323*sin(3*p - 9*t) - 561*sin(3*p - 7*t) - 525*sin(3*p - 5*t) - 358*sin(3*p - 3*t) - 126*sin(3*p - t) + 126*sin(3*p + t) + 358*sin(3*p + 3*t) + 525*sin(3*p + 5*t) + 561*sin(3*p + 7*t) + 323*sin(3*p + 9*t) - 969*sin(3*p + 11*t)); break;
      case -2 : return -0.0000845060337556615359449834435337*(2261*cos(2*p - 11*t) + 323*cos(2*p - 9*t) - 289*cos(2*p - 7*t) - 615*cos(2*p - 5*t) - 798*cos(2*p - 3*t) - 882*cos(2*p - t) - 882*cos(2*p + t) - 798*cos(2*p + 3*t) - 615*cos(2*p + 5*t) - 289*cos(2*p + 7*t) + 323*cos(2*p + 9*t) + 2261*cos(2*p + 11*t)); break;
      case -1 : return -7.41166945852334579335102325981e-6*(29393*sin(p - 11*t) + 12597*sin(p - 9*t) + 7735*sin(p - 7*t) + 4875*sin(p - 5*t) + 2730*sin(p - 3*t) + 882*sin(p - t) - 882*sin(p + t) - 2730*sin(p + 3*t) - 4875*sin(p + 5*t) - 7735*sin(p + 7*t) - 12597*sin(p + 9*t) - 29393*sin(p + 11*t)); break;
      case 0 : return 5.1608241842327211090243086402e-6*(29106*cos(t) + 30030*cos(3*t) + 32175*cos(5*t) + 36465*cos(7*t) + 46189*cos(9*t) + 88179*cos(11*t)); break;
      case 1 : return 7.41166945852334579335102325981e-6*(29393*sin(p - 11*t) + 12597*sin(p - 9*t) + 7735*sin(p - 7*t) + 4875*sin(p - 5*t) + 2730*sin(p - 3*t) + 882*sin(p - t) - 882*sin(p + t) - 2730*sin(p + 3*t) - 4875*sin(p + 5*t) - 7735*sin(p + 7*t) - 12597*sin(p + 9*t) - 29393*sin(p + 11*t)); break;
      case 2 : return -0.0000845060337556615359449834435337*(2261*cos(2*p - 11*t) + 323*cos(2*p - 9*t) - 289*cos(2*p - 7*t) - 615*cos(2*p - 5*t) - 798*cos(2*p - 3*t) - 882*cos(2*p - t) - 882*cos(2*p + t) - 798*cos(2*p + 3*t) - 615*cos(2*p + 5*t) - 289*cos(2*p + 7*t) + 323*cos(2*p + 9*t) + 2261*cos(2*p + 11*t)); break;
      case 3 : return -0.000158096312714419511069847760972*(969*sin(3*p - 11*t) - 323*sin(3*p - 9*t) - 561*sin(3*p - 7*t) - 525*sin(3*p - 5*t) - 358*sin(3*p - 3*t) - 126*sin(3*p - t) + 126*sin(3*p + t) + 358*sin(3*p + 3*t) + 525*sin(3*p + 5*t) + 561*sin(3*p + 7*t) + 323*sin(3*p + 9*t) - 969*sin(3*p + 11*t)); break;
      case 4 : return 0.000346371666928313460878751369526*(323*cos(4*p - 11*t) - 323*cos(4*p - 9*t) - 255*cos(4*p - 7*t) - 65*cos(4*p - 5*t) + 110*cos(4*p - 3*t) + 210*cos(4*p - t) + 210*cos(4*p + t) + 110*cos(4*p + 3*t) - 65*cos(4*p + 5*t) - 255*cos(4*p + 7*t) - 323*cos(4*p + 9*t) + 323*cos(4*p + 11*t)); break;
      case 5 : return 0.0000327290461389718940228154423877*(2261*sin(5*p - 11*t) - 4199*sin(5*p - 9*t) - 765*sin(5*p - 7*t) + 1975*sin(5*p - 5*t) + 2450*sin(5*p - 3*t) + 1050*sin(5*p - t) - 1050*sin(5*p + t) - 2450*sin(5*p + 3*t) - 1975*sin(5*p + 5*t) + 765*sin(5*p + 7*t) + 4199*sin(5*p + 9*t) - 2261*sin(5*p + 11*t)); break;
      case 6 : return -0.0001101823877028089813227240054718*(399*cos(6*p - 11*t) - 1159*cos(6*p - 9*t) + 525*cos(6*p - 7*t) + 795*cos(6*p - 5*t) + 70*cos(6*p - 3*t) - 630*cos(6*p - t) - 630*cos(6*p + t) + 70*cos(6*p + 3*t) + 795*cos(6*p + 5*t) + 525*cos(6*p + 7*t) - 1159*cos(6*p + 9*t) + 399*cos(6*p + 11*t)); break;
      case 7 : return -0.0001742136515883019988618247258659*(133*sin(7*p - 11*t) - 551*sin(7*p - 9*t) + 643*sin(7*p - 7*t) + 135*sin(7*p - 5*t) - 462*sin(7*p - 3*t) - 294*sin(7*p - t) + 294*sin(7*p + t) + 462*sin(7*p + 3*t) - 135*sin(7*p + 5*t) - 643*sin(7*p + 7*t) + 551*sin(7*p + 9*t) - 133*sin(7*p + 11*t)); break;
      case 8 : return 0.00151875940371722513571203636154*(7*cos(8*p - 11*t) - 39*cos(8*p - 9*t) + 77*cos(8*p - 7*t) - 45*cos(8*p - 5*t) - 42*cos(8*p - 3*t) + 42*cos(8*p - t) + 42*cos(8*p + t) - 42*cos(8*p + 3*t) - 45*cos(8*p + 5*t) + 77*cos(8*p + 7*t) - 39*cos(8*p + 9*t) + 7*cos(8*p + 11*t)); break;
      case 9 : return 0.0001960709959164239929600350650055*(21*sin(9*p - 11*t) - 151*sin(9*p - 9*t) + 435*sin(9*p - 7*t) - 585*sin(9*p - 5*t) + 210*sin(9*p - 3*t) + 378*sin(9*p - t) - 378*sin(9*p + t) - 210*sin(9*p + 3*t) + 585*sin(9*p + 5*t) - 435*sin(9*p + 7*t) + 151*sin(9*p + 9*t) - 21*sin(9*p + 11*t)); break;
      case 10 : return -0.001270685283012930339479964681051*(cos(10*p - 11*t) - 9*cos(10*p - 9*t) + 35*cos(10*p - 7*t) - 75*cos(10*p - 5*t) + 90*cos(10*p - 3*t) - 42*cos(10*p - t) - 42*cos(10*p + t) + 90*cos(10*p + 3*t) - 75*cos(10*p + 5*t) + 35*cos(10*p + 7*t) - 9*cos(10*p + 9*t) + cos(10*p + 11*t)); break;
      case 11 : return -0.000270911012600888330984112957788*(sin(11*p - 11*t) - 11*sin(11*p - 9*t) + 55*sin(11*p - 7*t) - 165*sin(11*p - 5*t) + 330*sin(11*p - 3*t) - 462*sin(11*p - t) + 462*sin(11*p + t) - 330*sin(11*p + 3*t) + 165*sin(11*p + 5*t) - 55*sin(11*p + 7*t) + 11*sin(11*p + 9*t) - sin(11*p + 11*t)); break;
    }
    case 12 : switch (m) 
    {
      case -12 : return 0.000138248697202683391346364508269*(924*cos(12*p) + cos(12*p - 12*t) - 12*cos(12*p - 10*t) + 66*cos(12*p - 8*t) - 220*cos(12*p - 6*t) + 495*cos(12*p - 4*t) - 792*cos(12*p - 2*t) - 792*cos(12*p + 2*t) + 495*cos(12*p + 4*t) - 220*cos(12*p + 6*t) + 66*cos(12*p + 8*t) - 12*cos(12*p + 10*t) + cos(12*p + 12*t)); break;
      case -11 : return 0.000677277531502220827460282394471*(sin(11*p - 12*t) - 10*sin(11*p - 10*t) + 44*sin(11*p - 8*t) - 110*sin(11*p - 6*t) + 165*sin(11*p - 4*t) - 132*sin(11*p - 2*t) + 132*sin(11*p + 2*t) - 165*sin(11*p + 4*t) + 110*sin(11*p + 6*t) - 44*sin(11*p + 8*t) + 10*sin(11*p + 10*t) - sin(11*p + 12*t)); break;
      case -10 : return 0.0000998591241044474053298361446224*(924*cos(10*p) - 23*cos(10*p - 12*t) + 188*cos(10*p - 10*t) - 638*cos(10*p - 8*t) + 1100*cos(10*p - 6*t) - 825*cos(10*p - 4*t) - 264*cos(10*p - 2*t) - 264*cos(10*p + 2*t) - 825*cos(10*p + 4*t) + 1100*cos(10*p + 6*t) - 638*cos(10*p + 8*t) + 188*cos(10*p + 10*t) - 23*cos(10*p + 12*t)); break;
      case -9 : return -0.000270419786425946424219471915251*(23*sin(9*p - 12*t) - 150*sin(9*p - 10*t) + 372*sin(9*p - 8*t) - 370*sin(9*p - 6*t) - 45*sin(9*p - 4*t) + 324*sin(9*p - 2*t) - 324*sin(9*p + 2*t) + 45*sin(9*p + 4*t) + 370*sin(9*p + 6*t) - 372*sin(9*p + 8*t) + 150*sin(9*p + 10*t) - 23*sin(9*p + 12*t)); break;
      case -8 : return 0.0000885156529079065189295199903972*(924*cos(8*p) + 161*cos(8*p - 12*t) - 812*cos(8*p - 10*t) + 1346*cos(8*p - 8*t) - 380*cos(8*p - 6*t) - 945*cos(8*p - 4*t) + 168*cos(8*p - 2*t) + 168*cos(8*p + 2*t) - 945*cos(8*p + 4*t) - 380*cos(8*p + 6*t) + 1346*cos(8*p + 8*t) - 812*cos(8*p + 10*t) + 161*cos(8*p + 12*t)); break;
      case -7 : return 0.0001770313058158130378590399807944*(161*sin(7*p - 12*t) - 602*sin(7*p - 10*t) + 556*sin(7*p - 8*t) + 290*sin(7*p - 6*t) - 315*sin(7*p - 4*t) - 420*sin(7*p - 2*t) + 420*sin(7*p + 2*t) + 315*sin(7*p + 4*t) - 290*sin(7*p + 6*t) - 556*sin(7*p + 8*t) + 602*sin(7*p + 10*t) - 161*sin(7*p + 12*t)); break;
      case -6 : return 0.00001658050092328703130490742729138*(4620*cos(6*p) - 3059*cos(6*p - 12*t) + 7980*cos(6*p - 10*t) - 2166*cos(6*p - 8*t) - 5380*cos(6*p - 6*t) - 2205*cos(6*p - 4*t) + 2520*cos(6*p - 2*t) + 2520*cos(6*p + 2*t) - 2205*cos(6*p + 4*t) - 5380*cos(6*p + 6*t) - 2166*cos(6*p + 8*t) + 7980*cos(6*p + 10*t) - 3059*cos(6*p + 12*t)); break;
      case -5 : return -0.0001861156612680872278147844505806*(437*sin(5*p - 12*t) - 722*sin(5*p - 10*t) - 228*sin(5*p - 8*t) + 250*sin(5*p - 6*t) + 425*sin(5*p - 4*t) + 300*sin(5*p - 2*t) - 300*sin(5*p + 2*t) - 425*sin(5*p + 4*t) - 250*sin(5*p + 6*t) + 228*sin(5*p + 8*t) + 722*sin(5*p + 10*t) - 437*sin(5*p + 12*t)); break;
      case -4 : return 0.00001595928629016969586604110799767*(4620*cos(4*p) + 7429*cos(4*p - 12*t) - 6460*cos(4*p - 10*t) - 5814*cos(4*p - 8*t) - 2380*cos(4*p - 6*t) + 1195*cos(4*p - 4*t) + 3720*cos(4*p - 2*t) + 3720*cos(4*p + 2*t) + 1195*cos(4*p + 4*t) - 2380*cos(4*p + 6*t) - 5814*cos(4*p + 8*t) - 6460*cos(4*p + 10*t) + 7429*cos(4*p + 12*t)); break;
      case -3 : return 0.0000212790483868929278213881439969*(7429*sin(3*p - 12*t) - 1938*sin(3*p - 10*t) - 3876*sin(3*p - 8*t) - 3910*sin(3*p - 6*t) - 3015*sin(3*p - 4*t) - 1620*sin(3*p - 2*t) + 1620*sin(3*p + 2*t) + 3015*sin(3*p + 4*t) + 3910*sin(3*p + 6*t) + 3876*sin(3*p + 8*t) + 1938*sin(3*p + 10*t) - 7429*sin(3*p + 12*t)); break;
      case -2 : return 0.0000260614053799405792989657278166*(2772*cos(2*p) - 7429*cos(2*p - 12*t) - 1292*cos(2*p - 10*t) + 646*cos(2*p - 8*t) + 1700*cos(2*p - 6*t) + 2325*cos(2*p - 4*t) + 2664*cos(2*p - 2*t) + 2664*cos(2*p + 2*t) + 2325*cos(2*p + 4*t) + 1700*cos(2*p + 6*t) + 646*cos(2*p + 8*t) - 1292*cos(2*p + 10*t) - 7429*cos(2*p + 12*t)); break;
      case -1 : return -4.20017578598614200022888376552e-6*(52003*sin(p - 12*t) + 22610*sin(p - 10*t) + 14212*sin(p - 8*t) + 9350*sin(p - 6*t) + 5775*sin(p - 4*t) + 2772*sin(p - 2*t) - 2772*sin(p + 2*t) - 5775*sin(p + 4*t) - 9350*sin(p + 6*t) - 14212*sin(p + 8*t) - 22610*sin(p + 10*t) - 52003*sin(p + 12*t)); break;
      case 0 : return 6.72566394266791685757731737567e-7*(106722 + 216216*cos(2*t) + 225225*cos(4*t) + 243100*cos(6*t) + 277134*cos(8*t) + 352716*cos(10*t) + 676039*cos(12*t)); break;
      case 1 : return 4.20017578598614200022888376552e-6*(52003*sin(p - 12*t) + 22610*sin(p - 10*t) + 14212*sin(p - 8*t) + 9350*sin(p - 6*t) + 5775*sin(p - 4*t) + 2772*sin(p - 2*t) - 2772*sin(p + 2*t) - 5775*sin(p + 4*t) - 9350*sin(p + 6*t) - 14212*sin(p + 8*t) - 22610*sin(p + 10*t) - 52003*sin(p + 12*t)); break;
      case 2 : return 0.0000260614053799405792989657278166*(2772*cos(2*p) - 7429*cos(2*p - 12*t) - 1292*cos(2*p - 10*t) + 646*cos(2*p - 8*t) + 1700*cos(2*p - 6*t) + 2325*cos(2*p - 4*t) + 2664*cos(2*p - 2*t) + 2664*cos(2*p + 2*t) + 2325*cos(2*p + 4*t) + 1700*cos(2*p + 6*t) + 646*cos(2*p + 8*t) - 1292*cos(2*p + 10*t) - 7429*cos(2*p + 12*t)); break;
      case 3 : return -0.0000212790483868929278213881439969*(7429*sin(3*p - 12*t) - 1938*sin(3*p - 10*t) - 3876*sin(3*p - 8*t) - 3910*sin(3*p - 6*t) - 3015*sin(3*p - 4*t) - 1620*sin(3*p - 2*t) + 1620*sin(3*p + 2*t) + 3015*sin(3*p + 4*t) + 3910*sin(3*p + 6*t) + 3876*sin(3*p + 8*t) + 1938*sin(3*p + 10*t) - 7429*sin(3*p + 12*t)); break;
      case 4 : return 0.00001595928629016969586604110799767*(4620*cos(4*p) + 7429*cos(4*p - 12*t) - 6460*cos(4*p - 10*t) - 5814*cos(4*p - 8*t) - 2380*cos(4*p - 6*t) + 1195*cos(4*p - 4*t) + 3720*cos(4*p - 2*t) + 3720*cos(4*p + 2*t) + 1195*cos(4*p + 4*t) - 2380*cos(4*p + 6*t) - 5814*cos(4*p + 8*t) - 6460*cos(4*p + 10*t) + 7429*cos(4*p + 12*t)); break;
      case 5 : return 0.0001861156612680872278147844505806*(437*sin(5*p - 12*t) - 722*sin(5*p - 10*t) - 228*sin(5*p - 8*t) + 250*sin(5*p - 6*t) + 425*sin(5*p - 4*t) + 300*sin(5*p - 2*t) - 300*sin(5*p + 2*t) - 425*sin(5*p + 4*t) - 250*sin(5*p + 6*t) + 228*sin(5*p + 8*t) + 722*sin(5*p + 10*t) - 437*sin(5*p + 12*t)); break;
      case 6 : return 0.00001658050092328703130490742729138*(4620*cos(6*p) - 3059*cos(6*p - 12*t) + 7980*cos(6*p - 10*t) - 2166*cos(6*p - 8*t) - 5380*cos(6*p - 6*t) - 2205*cos(6*p - 4*t) + 2520*cos(6*p - 2*t) + 2520*cos(6*p + 2*t) - 2205*cos(6*p + 4*t) - 5380*cos(6*p + 6*t) - 2166*cos(6*p + 8*t) + 7980*cos(6*p + 10*t) - 3059*cos(6*p + 12*t)); break;
      case 7 : return -0.0001770313058158130378590399807944*(161*sin(7*p - 12*t) - 602*sin(7*p - 10*t) + 556*sin(7*p - 8*t) + 290*sin(7*p - 6*t) - 315*sin(7*p - 4*t) - 420*sin(7*p - 2*t) + 420*sin(7*p + 2*t) + 315*sin(7*p + 4*t) - 290*sin(7*p + 6*t) - 556*sin(7*p + 8*t) + 602*sin(7*p + 10*t) - 161*sin(7*p + 12*t)); break;
      case 8 : return 0.0000885156529079065189295199903972*(924*cos(8*p) + 161*cos(8*p - 12*t) - 812*cos(8*p - 10*t) + 1346*cos(8*p - 8*t) - 380*cos(8*p - 6*t) - 945*cos(8*p - 4*t) + 168*cos(8*p - 2*t) + 168*cos(8*p + 2*t) - 945*cos(8*p + 4*t) - 380*cos(8*p + 6*t) + 1346*cos(8*p + 8*t) - 812*cos(8*p + 10*t) + 161*cos(8*p + 12*t)); break;
      case 9 : return 0.000270419786425946424219471915251*(23*sin(9*p - 12*t) - 150*sin(9*p - 10*t) + 372*sin(9*p - 8*t) - 370*sin(9*p - 6*t) - 45*sin(9*p - 4*t) + 324*sin(9*p - 2*t) - 324*sin(9*p + 2*t) + 45*sin(9*p + 4*t) + 370*sin(9*p + 6*t) - 372*sin(9*p + 8*t) + 150*sin(9*p + 10*t) - 23*sin(9*p + 12*t)); break;
      case 10 : return 0.0000998591241044474053298361446224*(924*cos(10*p) - 23*cos(10*p - 12*t) + 188*cos(10*p - 10*t) - 638*cos(10*p - 8*t) + 1100*cos(10*p - 6*t) - 825*cos(10*p - 4*t) - 264*cos(10*p - 2*t) - 264*cos(10*p + 2*t) - 825*cos(10*p + 4*t) + 1100*cos(10*p + 6*t) - 638*cos(10*p + 8*t) + 188*cos(10*p + 10*t) - 23*cos(10*p + 12*t)); break;
      case 11 : return -0.000677277531502220827460282394471*(sin(11*p - 12*t) - 10*sin(11*p - 10*t) + 44*sin(11*p - 8*t) - 110*sin(11*p - 6*t) + 165*sin(11*p - 4*t) - 132*sin(11*p - 2*t) + 132*sin(11*p + 2*t) - 165*sin(11*p + 4*t) + 110*sin(11*p + 6*t) - 44*sin(11*p + 8*t) + 10*sin(11*p + 10*t) - sin(11*p + 12*t)); break;
      case 12 : return 0.000138248697202683391346364508269*(924*cos(12*p) + cos(12*p - 12*t) - 12*cos(12*p - 10*t) + 66*cos(12*p - 8*t) - 220*cos(12*p - 6*t) + 495*cos(12*p - 4*t) - 792*cos(12*p - 2*t) - 792*cos(12*p + 2*t) + 495*cos(12*p + 4*t) - 220*cos(12*p + 6*t) + 66*cos(12*p + 8*t) - 12*cos(12*p + 10*t) + cos(12*p + 12*t)); break;
    }
  }
  return 0;
}

double ImSphericalHarmonicY(const int k, const int m, const double t, const double p)
{
  switch (k)
  {
    case 0 : switch (m) 
    {
      case 0 : return 0; break;
    }
    case 1 : switch (m) 
    {
      case -1 : return -0.1727470747356677396326223230159*(cos(p - t) - cos(p + t)); break;
      case 0 : return 0; break;
      case 1 : return -0.1727470747356677396326223230159*(cos(p - t) - cos(p + t)); break;
    }
    case 2 : switch (m) 
    {
      case -2 : return -0.0965685505057973950855481838279*(2*sin(2*p) - sin(2*p - 2*t) - sin(2*p + 2*t)); break;
      case -1 : return -0.1931371010115947901710963676557*(cos(p - 2*t) - cos(p + 2*t)); break;
      case 0 : return 0; break;
      case 1 : return -0.1931371010115947901710963676557*(cos(p - 2*t) - cos(p + 2*t)); break;
      case 2 : return 0.0965685505057973950855481838279*(2*sin(2*p) - sin(2*p - 2*t) - sin(2*p + 2*t)); break;
    }
    case 3 : switch (m) 
    {
      case -3 : return 0.0521529779540980112155533769671*(cos(3*p - 3*t) - 3*cos(3*p - t) + 3*cos(3*p + t) - cos(3*p + 3*t)); break;
      case -2 : return 0.1277481845541602954245143114626*(sin(2*p - 3*t) - sin(2*p - t) - sin(2*p + t) + sin(2*p + 3*t)); break;
      case -1 : return -0.0403975230142688316259242704163*(5*cos(p - 3*t) + cos(p - t) - cos(p + t) - 5*cos(p + 3*t)); break;
      case 0 : return 0; break;
      case 1 : return -0.0403975230142688316259242704163*(5*cos(p - 3*t) + cos(p - t) - cos(p + t) - 5*cos(p + 3*t)); break;
      case 2 : return -0.1277481845541602954245143114626*(sin(2*p - 3*t) - sin(2*p - t) - sin(2*p + t) + sin(2*p + 3*t)); break;
      case 3 : return 0.0521529779540980112155533769671*(cos(3*p - 3*t) - 3*cos(3*p - t) + 3*cos(3*p + t) - cos(3*p + 3*t)); break;
    }
    case 4 : switch (m) 
    {
      case -4 : return -0.0276582932778114145474504817598*(6*sin(4*p) + sin(4*p - 4*t) - 4*sin(4*p - 2*t) - 4*sin(4*p + 2*t) + sin(4*p + 4*t)); break;
      case -3 : return 0.0782294669311470168233300654506*(cos(3*p - 4*t) - 2*cos(3*p - 2*t) + 2*cos(3*p + 2*t) - cos(3*p + 4*t)); break;
      case -2 : return -0.0209077044861652864985035121785*(6*sin(2*p) - 7*sin(2*p - 4*t) + 4*sin(2*p - 2*t) + 4*sin(2*p + 2*t) - 7*sin(2*p + 4*t)); break;
      case -1 : return -0.0295679592424237505653962834652*(7*cos(p - 4*t) + 2*cos(p - 2*t) - 2*cos(p + 2*t) - 7*cos(p + 4*t)); break;
      case 0 : return 0; break;
      case 1 : return -0.0295679592424237505653962834652*(7*cos(p - 4*t) + 2*cos(p - 2*t) - 2*cos(p + 2*t) - 7*cos(p + 4*t)); break;
      case 2 : return 0.0209077044861652864985035121785*(6*sin(2*p) - 7*sin(2*p - 4*t) + 4*sin(2*p - 2*t) + 4*sin(2*p + 2*t) - 7*sin(2*p + 4*t)); break;
      case 3 : return 0.0782294669311470168233300654506*(cos(3*p - 4*t) - 2*cos(3*p - 2*t) + 2*cos(3*p + 2*t) - cos(3*p + 4*t)); break;
      case 4 : return 0.0276582932778114145474504817598*(6*sin(4*p) + sin(4*p - 4*t) - 4*sin(4*p - 2*t) - 4*sin(4*p + 2*t) + sin(4*p + 4*t)); break;
    }
    case 5 : switch (m) 
    {
      case -5 : return -0.01450413135752681752056245642293*(cos(5*p - 5*t) - 5*cos(5*p - 3*t) + 10*cos(5*p - t) - 10*cos(5*p + t) + 5*cos(5*p + 3*t) - cos(5*p + 5*t)); break;
      case -4 : return -0.0458660905720547238453758171103*(sin(4*p - 5*t) - 3*sin(4*p - 3*t) + 2*sin(4*p - t) + 2*sin(4*p + t) - 3*sin(4*p + 3*t) + sin(4*p + 5*t)); break;
      case -3 : return 0.01081074122333875666143645063852*(9*cos(3*p - 5*t) - 13*cos(3*p - 3*t) - 6*cos(3*p - t) + 6*cos(3*p + t) + 13*cos(3*p + 3*t) - 9*cos(3*p + 5*t)); break;
      case -2 : return 0.0529615994769031023692403611097*(3*sin(2*p - 5*t) - sin(2*p - 3*t) - 2*sin(2*p - t) - 2*sin(2*p + t) - sin(2*p + 3*t) + 3*sin(2*p + 5*t)); break;
      case -1 : return -0.01000880151800672273553629477647*(21*cos(p - 5*t) + 7*cos(p - 3*t) + 2*cos(p - t) - 2*cos(p + t) - 7*cos(p + 3*t) - 21*cos(p + 5*t)); break;
      case 0 : return 0; break;
      case 1 : return -0.01000880151800672273553629477647*(21*cos(p - 5*t) + 7*cos(p - 3*t) + 2*cos(p - t) - 2*cos(p + t) - 7*cos(p + 3*t) - 21*cos(p + 5*t)); break;
      case 2 : return -0.0529615994769031023692403611097*(3*sin(2*p - 5*t) - sin(2*p - 3*t) - 2*sin(2*p - t) - 2*sin(2*p + t) - sin(2*p + 3*t) + 3*sin(2*p + 5*t)); break;
      case 3 : return 0.01081074122333875666143645063852*(9*cos(3*p - 5*t) - 13*cos(3*p - 3*t) - 6*cos(3*p - t) + 6*cos(3*p + t) + 13*cos(3*p + 3*t) - 9*cos(3*p + 5*t)); break;
      case 4 : return 0.0458660905720547238453758171103*(sin(4*p - 5*t) - 3*sin(4*p - 3*t) + 2*sin(4*p - t) + 2*sin(4*p + t) - 3*sin(4*p + 3*t) + sin(4*p + 5*t)); break;
      case 5 : return -0.01450413135752681752056245642293*(cos(5*p - 5*t) - 5*cos(5*p - 3*t) + 10*cos(5*p - t) - 10*cos(5*p + t) + 5*cos(5*p + 3*t) - cos(5*p + 5*t)); break;
    }
    case 6 : switch (m) 
    {
      case -6 : return -0.00754818927468853479645251755055*(20*sin(6*p) - sin(6*p - 6*t) + 6*sin(6*p - 4*t) - 15*sin(6*p - 2*t) - 15*sin(6*p + 2*t) + 6*sin(6*p + 4*t) - sin(6*p + 6*t)); break;
      case -5 : return -0.0261476946578140297069217513434*(cos(5*p - 6*t) - 4*cos(5*p - 4*t) + 5*cos(5*p - 2*t) - 5*cos(5*p + 2*t) + 4*cos(5*p + 4*t) - cos(5*p + 6*t)); break;
      case -4 : return -0.0055747072320937191815291748398*(20*sin(4*p) + 11*sin(4*p - 6*t) - 26*sin(4*p - 4*t) + 5*sin(4*p - 2*t) + 5*sin(4*p + 2*t) - 26*sin(4*p + 4*t) + 11*sin(4*p + 6*t)); break;
      case -3 : return 0.01017797634168305839948117512343*(11*cos(3*p - 6*t) - 12*cos(3*p - 4*t) - 9*cos(3*p - 2*t) + 9*cos(3*p + 2*t) + 12*cos(3*p + 4*t) - 11*cos(3*p + 6*t)); break;
      case -2 : return -0.00508898817084152919974058756172*(20*sin(2*p) - 33*sin(2*p - 6*t) + 6*sin(2*p - 4*t) + 17*sin(2*p - 2*t) + 17*sin(2*p + 2*t) + 6*sin(2*p + 4*t) - 33*sin(2*p + 6*t)); break;
      case -1 : return -0.00643711744220532464693356226676*(33*cos(p - 6*t) + 12*cos(p - 4*t) + 5*cos(p - 2*t) - 5*cos(p + 2*t) - 12*cos(p + 4*t) - 33*cos(p + 6*t)); break;
      case 0 : return 0; break;
      case 1 : return -0.00643711744220532464693356226676*(33*cos(p - 6*t) + 12*cos(p - 4*t) + 5*cos(p - 2*t) - 5*cos(p + 2*t) - 12*cos(p + 4*t) - 33*cos(p + 6*t)); break;
      case 2 : return 0.00508898817084152919974058756172*(20*sin(2*p) - 33*sin(2*p - 6*t) + 6*sin(2*p - 4*t) + 17*sin(2*p - 2*t) + 17*sin(2*p + 2*t) + 6*sin(2*p + 4*t) - 33*sin(2*p + 6*t)); break;
      case 3 : return 0.01017797634168305839948117512343*(11*cos(3*p - 6*t) - 12*cos(3*p - 4*t) - 9*cos(3*p - 2*t) + 9*cos(3*p + 2*t) + 12*cos(3*p + 4*t) - 11*cos(3*p + 6*t)); break;
      case 4 : return 0.0055747072320937191815291748398*(20*sin(4*p) + 11*sin(4*p - 6*t) - 26*sin(4*p - 4*t) + 5*sin(4*p - 2*t) + 5*sin(4*p + 2*t) - 26*sin(4*p + 4*t) + 11*sin(4*p + 6*t)); break;
      case 5 : return -0.0261476946578140297069217513434*(cos(5*p - 6*t) - 4*cos(5*p - 4*t) + 5*cos(5*p - 2*t) - 5*cos(5*p + 2*t) + 4*cos(5*p + 4*t) - cos(5*p + 6*t)); break;
      case 6 : return 0.00754818927468853479645251755055*(20*sin(6*p) - sin(6*p - 6*t) + 6*sin(6*p - 4*t) - 15*sin(6*p - 2*t) - 15*sin(6*p + 2*t) + 6*sin(6*p + 4*t) - sin(6*p + 6*t)); break;
    }
    case 7 : switch (m) 
    {
      case -7 : return 0.00390655909039492706258854893785*(cos(7*p - 7*t) - 7*cos(7*p - 5*t) + 21*cos(7*p - 3*t) - 35*cos(7*p - t) + 35*cos(7*p + t) - 21*cos(7*p + 3*t) + 7*cos(7*p + 5*t) - cos(7*p + 7*t)); break;
      case -6 : return 0.01461700567744506825595076737686*(sin(6*p - 7*t) - 5*sin(6*p - 5*t) + 9*sin(6*p - 3*t) - 5*sin(6*p - t) - 5*sin(6*p + t) + 9*sin(6*p + 3*t) - 5*sin(6*p + 5*t) + sin(6*p + 7*t)); break;
      case -5 : return -0.00286663066075342024033598856939*(13*cos(5*p - 7*t) - 43*cos(5*p - 5*t) + 33*cos(5*p - 3*t) + 25*cos(5*p - t) - 25*cos(5*p + t) - 33*cos(5*p + 3*t) + 43*cos(5*p + 5*t) - 13*cos(5*p + 7*t)); break;
      case -4 : return -0.00573326132150684048067197713878*(13*sin(4*p - 7*t) - 25*sin(4*p - 5*t) - 3*sin(4*p - 3*t) + 15*sin(4*p - t) + 15*sin(4*p + t) - 3*sin(4*p + 3*t) - 25*sin(4*p + 5*t) + 13*sin(4*p + 7*t)); break;
      case -3 : return 0.000864321664931606704607827554992*(143*cos(3*p - 7*t) - 121*cos(3*p - 5*t) - 117*cos(3*p - 3*t) - 45*cos(3*p - t) + 45*cos(3*p + t) + 117*cos(3*p + 3*t) + 121*cos(3*p + 5*t) - 143*cos(3*p + 7*t)); break;
      case -2 : return 0.001222335420799172137864532272666*(143*sin(2*p - 7*t) - 11*sin(2*p - 5*t) - 57*sin(2*p - 3*t) - 75*sin(2*p - t) - 75*sin(2*p + t) - 57*sin(2*p + 3*t) - 11*sin(2*p + 5*t) + 143*sin(2*p + 7*t)); break;
      case -1 : return -0.000499016345914688654002009029151*(429*cos(p - 7*t) + 165*cos(p - 5*t) + 81*cos(p - 3*t) + 25*cos(p - t) - 25*cos(p + t) - 81*cos(p + 3*t) - 165*cos(p + 5*t) - 429*cos(p + 7*t)); break;
      case 0 : return 0; break;
      case 1 : return -0.000499016345914688654002009029151*(429*cos(p - 7*t) + 165*cos(p - 5*t) + 81*cos(p - 3*t) + 25*cos(p - t) - 25*cos(p + t) - 81*cos(p + 3*t) - 165*cos(p + 5*t) - 429*cos(p + 7*t)); break;
      case 2 : return -0.001222335420799172137864532272666*(143*sin(2*p - 7*t) - 11*sin(2*p - 5*t) - 57*sin(2*p - 3*t) - 75*sin(2*p - t) - 75*sin(2*p + t) - 57*sin(2*p + 3*t) - 11*sin(2*p + 5*t) + 143*sin(2*p + 7*t)); break;
      case 3 : return 0.000864321664931606704607827554992*(143*cos(3*p - 7*t) - 121*cos(3*p - 5*t) - 117*cos(3*p - 3*t) - 45*cos(3*p - t) + 45*cos(3*p + t) + 117*cos(3*p + 3*t) + 121*cos(3*p + 5*t) - 143*cos(3*p + 7*t)); break;
      case 4 : return 0.00573326132150684048067197713878*(13*sin(4*p - 7*t) - 25*sin(4*p - 5*t) - 3*sin(4*p - 3*t) + 15*sin(4*p - t) + 15*sin(4*p + t) - 3*sin(4*p + 3*t) - 25*sin(4*p + 5*t) + 13*sin(4*p + 7*t)); break;
      case 5 : return -0.00286663066075342024033598856939*(13*cos(5*p - 7*t) - 43*cos(5*p - 5*t) + 33*cos(5*p - 3*t) + 25*cos(5*p - t) - 25*cos(5*p + t) - 33*cos(5*p + 3*t) + 43*cos(5*p + 5*t) - 13*cos(5*p + 7*t)); break;
      case 6 : return -0.01461700567744506825595076737686*(sin(6*p - 7*t) - 5*sin(6*p - 5*t) + 9*sin(6*p - 3*t) - 5*sin(6*p - t) - 5*sin(6*p + t) + 9*sin(6*p + 3*t) - 5*sin(6*p + 5*t) + sin(6*p + 7*t)); break;
      case 7 : return 0.00390655909039492706258854893785*(cos(7*p - 7*t) - 7*cos(7*p - 5*t) + 21*cos(7*p - 3*t) - 35*cos(7*p - t) + 35*cos(7*p + t) - 21*cos(7*p + 3*t) + 7*cos(7*p + 5*t) - cos(7*p + 7*t)); break;
    }
    case 8 : switch (m) 
    {
      case -8 : return -0.00201339447030189183486334132529*(70*sin(8*p) + sin(8*p - 8*t) - 8*sin(8*p - 6*t) + 28*sin(8*p - 4*t) - 56*sin(8*p - 2*t) - 56*sin(8*p + 2*t) + 28*sin(8*p + 4*t) - 8*sin(8*p + 6*t) + sin(8*p + 8*t)); break;
      case -7 : return 0.00805357788120756733945336530116*(cos(7*p - 8*t) - 6*cos(7*p - 6*t) + 14*cos(7*p - 4*t) - 14*cos(7*p - 2*t) + 14*cos(7*p + 2*t) - 14*cos(7*p + 4*t) + 6*cos(7*p + 6*t) - cos(7*p + 8*t)); break;
      case -6 : return -0.00147037542472068188952924283307*(70*sin(6*p) - 15*sin(6*p - 8*t) + 64*sin(6*p - 6*t) - 84*sin(6*p - 4*t) - 84*sin(6*p + 4*t) + 64*sin(6*p + 6*t) - 15*sin(6*p + 8*t)); break;
      case -5 : return -0.00952912185692606606411696386973*(5*cos(5*p - 8*t) - 14*cos(5*p - 6*t) + 6*cos(5*p - 4*t) + 10*cos(5*p - 2*t) - 10*cos(5*p + 2*t) - 6*cos(5*p + 4*t) + 14*cos(5*p + 6*t) - 5*cos(5*p + 8*t)); break;
      case -4 : return -0.001321451440972752138178969129294*(70*sin(4*p) + 65*sin(4*p - 8*t) - 104*sin(4*p - 6*t) - 36*sin(4*p - 4*t) + 40*sin(4*p - 2*t) + 40*sin(4*p + 2*t) - 36*sin(4*p + 4*t) - 104*sin(4*p + 6*t) + 65*sin(4*p + 8*t)); break;
      case -3 : return 0.00341197294913950827517017222963*(39*cos(3*p - 8*t) - 26*cos(3*p - 6*t) - 30*cos(3*p - 4*t) - 18*cos(3*p - 2*t) + 18*cos(3*p + 2*t) + 30*cos(3*p + 4*t) + 26*cos(3*p + 6*t) - 39*cos(3*p + 8*t)); break;
      case -2 : return -0.001259954512472199251262382386198*(70*sin(2*p) - 143*sin(2*p - 8*t) + 44*sin(2*p - 4*t) + 64*sin(2*p - 2*t) + 64*sin(2*p + 2*t) + 44*sin(2*p + 4*t) - 143*sin(2*p + 8*t)); break;
      case -1 : return -0.000301186735953348127575578792952*(715*cos(p - 8*t) + 286*cos(p - 6*t) + 154*cos(p - 4*t) + 70*cos(p - 2*t) - 70*cos(p + 2*t) - 154*cos(p + 4*t) - 286*cos(p + 6*t) - 715*cos(p + 8*t)); break;
      case 0 : return 0; break;
      case 1 : return -0.000301186735953348127575578792952*(715*cos(p - 8*t) + 286*cos(p - 6*t) + 154*cos(p - 4*t) + 70*cos(p - 2*t) - 70*cos(p + 2*t) - 154*cos(p + 4*t) - 286*cos(p + 6*t) - 715*cos(p + 8*t)); break;
      case 2 : return 0.001259954512472199251262382386198*(70*sin(2*p) - 143*sin(2*p - 8*t) + 44*sin(2*p - 4*t) + 64*sin(2*p - 2*t) + 64*sin(2*p + 2*t) + 44*sin(2*p + 4*t) - 143*sin(2*p + 8*t)); break;
      case 3 : return 0.00341197294913950827517017222963*(39*cos(3*p - 8*t) - 26*cos(3*p - 6*t) - 30*cos(3*p - 4*t) - 18*cos(3*p - 2*t) + 18*cos(3*p + 2*t) + 30*cos(3*p + 4*t) + 26*cos(3*p + 6*t) - 39*cos(3*p + 8*t)); break;
      case 4 : return 0.001321451440972752138178969129294*(70*sin(4*p) + 65*sin(4*p - 8*t) - 104*sin(4*p - 6*t) - 36*sin(4*p - 4*t) + 40*sin(4*p - 2*t) + 40*sin(4*p + 2*t) - 36*sin(4*p + 4*t) - 104*sin(4*p + 6*t) + 65*sin(4*p + 8*t)); break;
      case 5 : return -0.00952912185692606606411696386973*(5*cos(5*p - 8*t) - 14*cos(5*p - 6*t) + 6*cos(5*p - 4*t) + 10*cos(5*p - 2*t) - 10*cos(5*p + 2*t) - 6*cos(5*p + 4*t) + 14*cos(5*p + 6*t) - 5*cos(5*p + 8*t)); break;
      case 6 : return 0.00147037542472068188952924283307*(70*sin(6*p) - 15*sin(6*p - 8*t) + 64*sin(6*p - 6*t) - 84*sin(6*p - 4*t) - 84*sin(6*p + 4*t) + 64*sin(6*p + 6*t) - 15*sin(6*p + 8*t)); break;
      case 7 : return 0.00805357788120756733945336530116*(cos(7*p - 8*t) - 6*cos(7*p - 6*t) + 14*cos(7*p - 4*t) - 14*cos(7*p - 2*t) + 14*cos(7*p + 2*t) - 14*cos(7*p + 4*t) + 6*cos(7*p + 6*t) - cos(7*p + 8*t)); break;
      case 8 : return 0.00201339447030189183486334132529*(70*sin(8*p) + sin(8*p - 8*t) - 8*sin(8*p - 6*t) + 28*sin(8*p - 4*t) - 56*sin(8*p - 2*t) - 56*sin(8*p + 2*t) + 28*sin(8*p + 4*t) - 8*sin(8*p + 6*t) + sin(8*p + 8*t)); break;
    }
    case 9 : switch (m) 
    {
      case -9 : return -0.001034283088852440585058962213157*(cos(9*p - 9*t) - 9*cos(9*p - 7*t) + 36*cos(9*p - 5*t) - 84*cos(9*p - 3*t) + 126*cos(9*p - t) - 126*cos(9*p + t) + 84*cos(9*p + 3*t) - 36*cos(9*p + 5*t) + 9*cos(9*p + 7*t) - cos(9*p + 9*t)); break;
      case -8 : return -0.0043880915147647751750992848194*(sin(8*p - 9*t) - 7*sin(8*p - 7*t) + 20*sin(8*p - 5*t) - 28*sin(8*p - 3*t) + 14*sin(8*p - t) + 14*sin(8*p + t) - 28*sin(8*p + 3*t) + 20*sin(8*p + 5*t) - 7*sin(8*p + 7*t) + sin(8*p + 9*t)); break;
      case -7 : return 0.000752551486258007370452908696331*(17*cos(7*p - 9*t) - 89*cos(7*p - 7*t) + 164*cos(7*p - 5*t) - 84*cos(7*p - 3*t) - 98*cos(7*p - t) + 98*cos(7*p + t) + 84*cos(7*p + 3*t) - 164*cos(7*p + 5*t) + 89*cos(7*p + 7*t) - 17*cos(7*p + 9*t)); break;
      case -6 : return 0.001737943212680454038596561231567*(17*sin(6*p - 9*t) - 63*sin(6*p - 7*t) + 60*sin(6*p - 5*t) + 28*sin(6*p - 3*t) - 42*sin(6*p - t) - 42*sin(6*p + t) + 28*sin(6*p + 3*t) + 60*sin(6*p + 5*t) - 63*sin(6*p + 7*t) + 17*sin(6*p + 9*t)); break;
      case -5 : return -0.000673102511936561327899554822462*(85*cos(5*p - 9*t) - 205*cos(5*p - 7*t) + 36*cos(5*p - 5*t) + 140*cos(5*p - 3*t) + 70*cos(5*p - t) - 70*cos(5*p + t) - 140*cos(5*p + 3*t) - 36*cos(5*p + 5*t) + 205*cos(5*p + 7*t) - 85*cos(5*p + 9*t)); break;
      case -4 : return -0.00563157965496996303859039131309*(17*sin(4*p - 9*t) - 23*sin(4*p - 7*t) - 12*sin(4*p - 5*t) + 4*sin(4*p - 3*t) + 14*sin(4*p - t) + 14*sin(4*p + t) + 4*sin(4*p + 3*t) - 12*sin(4*p + 5*t) - 23*sin(4*p + 7*t) + 17*sin(4*p + 9*t)); break;
      case -3 : return 0.000637650830927843559575950527255*(221*cos(3*p - 9*t) - 117*cos(3*p - 7*t) - 156*cos(3*p - 5*t) - 116*cos(3*p - 3*t) - 42*cos(3*p - t) + 42*cos(3*p + t) + 116*cos(3*p + 3*t) + 156*cos(3*p + 5*t) + 117*cos(3*p + 7*t) - 221*cos(3*p + 9*t)); break;
      case -2 : return 0.00083488091419380901096632445383*(221*sin(2*p - 9*t) + 13*sin(2*p - 7*t) - 52*sin(2*p - 5*t) - 84*sin(2*p - 3*t) - 98*sin(2*p - t) - 98*sin(2*p + t) - 84*sin(2*p + 3*t) - 52*sin(2*p + 5*t) + 13*sin(2*p + 7*t) + 221*sin(2*p + 9*t)); break;
      case -1 : return -0.000088998604488873504194604687416*(2431*cos(p - 9*t) + 1001*cos(p - 7*t) + 572*cos(p - 5*t) + 308*cos(p - 3*t) + 98*cos(p - t) - 98*cos(p + t) - 308*cos(p + 3*t) - 572*cos(p + 5*t) - 1001*cos(p + 7*t) - 2431*cos(p + 9*t)); break;
      case 0 : return 0; break;
      case 1 : return -0.000088998604488873504194604687416*(2431*cos(p - 9*t) + 1001*cos(p - 7*t) + 572*cos(p - 5*t) + 308*cos(p - 3*t) + 98*cos(p - t) - 98*cos(p + t) - 308*cos(p + 3*t) - 572*cos(p + 5*t) - 1001*cos(p + 7*t) - 2431*cos(p + 9*t)); break;
      case 2 : return -0.00083488091419380901096632445383*(221*sin(2*p - 9*t) + 13*sin(2*p - 7*t) - 52*sin(2*p - 5*t) - 84*sin(2*p - 3*t) - 98*sin(2*p - t) - 98*sin(2*p + t) - 84*sin(2*p + 3*t) - 52*sin(2*p + 5*t) + 13*sin(2*p + 7*t) + 221*sin(2*p + 9*t)); break;
      case 3 : return 0.000637650830927843559575950527255*(221*cos(3*p - 9*t) - 117*cos(3*p - 7*t) - 156*cos(3*p - 5*t) - 116*cos(3*p - 3*t) - 42*cos(3*p - t) + 42*cos(3*p + t) + 116*cos(3*p + 3*t) + 156*cos(3*p + 5*t) + 117*cos(3*p + 7*t) - 221*cos(3*p + 9*t)); break;
      case 4 : return 0.00563157965496996303859039131309*(17*sin(4*p - 9*t) - 23*sin(4*p - 7*t) - 12*sin(4*p - 5*t) + 4*sin(4*p - 3*t) + 14*sin(4*p - t) + 14*sin(4*p + t) + 4*sin(4*p + 3*t) - 12*sin(4*p + 5*t) - 23*sin(4*p + 7*t) + 17*sin(4*p + 9*t)); break;
      case 5 : return -0.000673102511936561327899554822462*(85*cos(5*p - 9*t) - 205*cos(5*p - 7*t) + 36*cos(5*p - 5*t) + 140*cos(5*p - 3*t) + 70*cos(5*p - t) - 70*cos(5*p + t) - 140*cos(5*p + 3*t) - 36*cos(5*p + 5*t) + 205*cos(5*p + 7*t) - 85*cos(5*p + 9*t)); break;
      case 6 : return -0.001737943212680454038596561231567*(17*sin(6*p - 9*t) - 63*sin(6*p - 7*t) + 60*sin(6*p - 5*t) + 28*sin(6*p - 3*t) - 42*sin(6*p - t) - 42*sin(6*p + t) + 28*sin(6*p + 3*t) + 60*sin(6*p + 5*t) - 63*sin(6*p + 7*t) + 17*sin(6*p + 9*t)); break;
      case 7 : return 0.000752551486258007370452908696331*(17*cos(7*p - 9*t) - 89*cos(7*p - 7*t) + 164*cos(7*p - 5*t) - 84*cos(7*p - 3*t) - 98*cos(7*p - t) + 98*cos(7*p + t) + 84*cos(7*p + 3*t) - 164*cos(7*p + 5*t) + 89*cos(7*p + 7*t) - 17*cos(7*p + 9*t)); break;
      case 8 : return 0.0043880915147647751750992848194*(sin(8*p - 9*t) - 7*sin(8*p - 7*t) + 20*sin(8*p - 5*t) - 28*sin(8*p - 3*t) + 14*sin(8*p - t) + 14*sin(8*p + t) - 28*sin(8*p + 3*t) + 20*sin(8*p + 5*t) - 7*sin(8*p + 7*t) + sin(8*p + 9*t)); break;
      case 9 : return -0.001034283088852440585058962213157*(cos(9*p - 9*t) - 9*cos(9*p - 7*t) + 36*cos(9*p - 5*t) - 84*cos(9*p - 3*t) + 126*cos(9*p - t) - 126*cos(9*p + t) + 84*cos(9*p + 3*t) - 36*cos(9*p + 5*t) + 9*cos(9*p + 7*t) - cos(9*p + 9*t)); break;
    }
    case 10 : switch (m) 
    {
      case -10 : return -0.000529912394476778770252132215025*(252*sin(10*p) - sin(10*p - 10*t) + 10*sin(10*p - 8*t) - 45*sin(10*p - 6*t) + 120*sin(10*p - 4*t) - 210*sin(10*p - 2*t) - 210*sin(10*p + 2*t) + 120*sin(10*p + 4*t) - 45*sin(10*p + 6*t) + 10*sin(10*p + 8*t) - sin(10*p + 10*t)); break;
      case -9 : return -0.00236984027233952286606648120237*(cos(9*p - 10*t) - 8*cos(9*p - 8*t) + 27*cos(9*p - 6*t) - 48*cos(9*p - 4*t) + 42*cos(9*p - 2*t) - 42*cos(9*p + 2*t) + 48*cos(9*p + 4*t) - 27*cos(9*p + 6*t) + 8*cos(9*p + 8*t) - cos(9*p + 10*t)); break;
      case -8 : return -0.000384438856831830722850812827095*(252*sin(8*p) + 19*sin(8*p - 10*t) - 118*sin(8*p - 8*t) + 279*sin(8*p - 6*t) - 264*sin(8*p - 4*t) - 42*sin(8*p - 2*t) - 42*sin(8*p + 2*t) - 264*sin(8*p + 4*t) + 279*sin(8*p + 6*t) - 118*sin(8*p + 8*t) + 19*sin(8*p + 10*t)); break;
      case -7 : return 0.000941679036536860067470202171862*(19*cos(7*p - 10*t) - 88*cos(7*p - 8*t) + 129*cos(7*p - 6*t) - 16*cos(7*p - 4*t) - 98*cos(7*p - 2*t) + 98*cos(7*p + 2*t) + 16*cos(7*p + 4*t) - 129*cos(7*p + 6*t) + 88*cos(7*p + 8*t) - 19*cos(7*p + 10*t)); break;
      case -6 : return -0.0001141953568550395957496153701317*(756*sin(6*p) - 323*sin(6*p - 10*t) + 1054*sin(6*p - 8*t) - 711*sin(6*p - 6*t) - 664*sin(6*p - 4*t) + 266*sin(6*p - 2*t) + 266*sin(6*p + 2*t) - 664*sin(6*p + 4*t) - 711*sin(6*p + 6*t) + 1054*sin(6*p + 8*t) - 323*sin(6*p + 10*t)); break;
      case -5 : return -0.00020427886451417210708532618237*(323*cos(5*p - 10*t) - 680*cos(5*p - 8*t) - 15*cos(5*p - 6*t) + 400*cos(5*p - 4*t) + 350*cos(5*p - 2*t) - 350*cos(5*p + 2*t) - 400*cos(5*p + 4*t) + 15*cos(5*p + 6*t) + 680*cos(5*p + 8*t) - 323*cos(5*p + 10*t)); break;
      case -4 : return -0.00032299324484886477318972173938*(252*sin(4*p) + 323*sin(4*p - 10*t) - 374*sin(4*p - 8*t) - 249*sin(4*p - 6*t) - 8*sin(4*p - 4*t) + 182*sin(4*p - 2*t) + 182*sin(4*p + 2*t) - 8*sin(4*p + 4*t) - 249*sin(4*p + 6*t) - 374*sin(4*p + 8*t) + 323*sin(4*p + 10*t)); break;
      case -3 : return 0.000456781427420158382998461480527*(323*cos(3*p - 10*t) - 136*cos(3*p - 8*t) - 207*cos(3*p - 6*t) - 176*cos(3*p - 4*t) - 98*cos(3*p - 2*t) + 98*cos(3*p + 2*t) + 176*cos(3*p + 4*t) + 207*cos(3*p + 6*t) + 136*cos(3*p + 8*t) - 323*cos(3*p + 10*t)); break;
      case -2 : return -0.0000447911040742721912552498213968*(1764*sin(2*p) - 4199*sin(2*p - 10*t) - 442*sin(2*p - 8*t) + 741*sin(2*p - 6*t) + 1352*sin(2*p - 4*t) + 1666*sin(2*p - 2*t) + 1666*sin(2*p + 2*t) + 1352*sin(2*p + 4*t) + 741*sin(2*p + 6*t) - 442*sin(2*p + 8*t) - 4199*sin(2*p + 10*t)); break;
      case -1 : return -0.0000517203119891631859074612901462*(4199*cos(p - 10*t) + 1768*cos(p - 8*t) + 1053*cos(p - 6*t) + 624*cos(p - 4*t) + 294*cos(p - 2*t) - 294*cos(p + 2*t) - 624*cos(p + 4*t) - 1053*cos(p + 6*t) - 1768*cos(p + 8*t) - 4199*cos(p + 10*t)); break;
      case 0 : return 0; break;
      case 1 : return -0.0000517203119891631859074612901462*(4199*cos(p - 10*t) + 1768*cos(p - 8*t) + 1053*cos(p - 6*t) + 624*cos(p - 4*t) + 294*cos(p - 2*t) - 294*cos(p + 2*t) - 624*cos(p + 4*t) - 1053*cos(p + 6*t) - 1768*cos(p + 8*t) - 4199*cos(p + 10*t)); break;
      case 2 : return 0.0000447911040742721912552498213968*(1764*sin(2*p) - 4199*sin(2*p - 10*t) - 442*sin(2*p - 8*t) + 741*sin(2*p - 6*t) + 1352*sin(2*p - 4*t) + 1666*sin(2*p - 2*t) + 1666*sin(2*p + 2*t) + 1352*sin(2*p + 4*t) + 741*sin(2*p + 6*t) - 442*sin(2*p + 8*t) - 4199*sin(2*p + 10*t)); break;
      case 3 : return 0.000456781427420158382998461480527*(323*cos(3*p - 10*t) - 136*cos(3*p - 8*t) - 207*cos(3*p - 6*t) - 176*cos(3*p - 4*t) - 98*cos(3*p - 2*t) + 98*cos(3*p + 2*t) + 176*cos(3*p + 4*t) + 207*cos(3*p + 6*t) + 136*cos(3*p + 8*t) - 323*cos(3*p + 10*t)); break;
      case 4 : return 0.00032299324484886477318972173938*(252*sin(4*p) + 323*sin(4*p - 10*t) - 374*sin(4*p - 8*t) - 249*sin(4*p - 6*t) - 8*sin(4*p - 4*t) + 182*sin(4*p - 2*t) + 182*sin(4*p + 2*t) - 8*sin(4*p + 4*t) - 249*sin(4*p + 6*t) - 374*sin(4*p + 8*t) + 323*sin(4*p + 10*t)); break;
      case 5 : return -0.00020427886451417210708532618237*(323*cos(5*p - 10*t) - 680*cos(5*p - 8*t) - 15*cos(5*p - 6*t) + 400*cos(5*p - 4*t) + 350*cos(5*p - 2*t) - 350*cos(5*p + 2*t) - 400*cos(5*p + 4*t) + 15*cos(5*p + 6*t) + 680*cos(5*p + 8*t) - 323*cos(5*p + 10*t)); break;
      case 6 : return 0.0001141953568550395957496153701317*(756*sin(6*p) - 323*sin(6*p - 10*t) + 1054*sin(6*p - 8*t) - 711*sin(6*p - 6*t) - 664*sin(6*p - 4*t) + 266*sin(6*p - 2*t) + 266*sin(6*p + 2*t) - 664*sin(6*p + 4*t) - 711*sin(6*p + 6*t) + 1054*sin(6*p + 8*t) - 323*sin(6*p + 10*t)); break;
      case 7 : return 0.000941679036536860067470202171862*(19*cos(7*p - 10*t) - 88*cos(7*p - 8*t) + 129*cos(7*p - 6*t) - 16*cos(7*p - 4*t) - 98*cos(7*p - 2*t) + 98*cos(7*p + 2*t) + 16*cos(7*p + 4*t) - 129*cos(7*p + 6*t) + 88*cos(7*p + 8*t) - 19*cos(7*p + 10*t)); break;
      case 8 : return 0.000384438856831830722850812827095*(252*sin(8*p) + 19*sin(8*p - 10*t) - 118*sin(8*p - 8*t) + 279*sin(8*p - 6*t) - 264*sin(8*p - 4*t) - 42*sin(8*p - 2*t) - 42*sin(8*p + 2*t) - 264*sin(8*p + 4*t) + 279*sin(8*p + 6*t) - 118*sin(8*p + 8*t) + 19*sin(8*p + 10*t)); break;
      case 9 : return -0.00236984027233952286606648120237*(cos(9*p - 10*t) - 8*cos(9*p - 8*t) + 27*cos(9*p - 6*t) - 48*cos(9*p - 4*t) + 42*cos(9*p - 2*t) - 42*cos(9*p + 2*t) + 48*cos(9*p + 4*t) - 27*cos(9*p + 6*t) + 8*cos(9*p + 8*t) - cos(9*p + 10*t)); break;
      case 10 : return 0.000529912394476778770252132215025*(252*sin(10*p) - sin(10*p - 10*t) + 10*sin(10*p - 8*t) - 45*sin(10*p - 6*t) + 120*sin(10*p - 4*t) - 210*sin(10*p - 2*t) - 210*sin(10*p + 2*t) + 120*sin(10*p + 4*t) - 45*sin(10*p + 6*t) + 10*sin(10*p + 8*t) - sin(10*p + 10*t)); break;
    }
    case 11 : switch (m) 
    {
      case -11 : return 0.000270911012600888330984112957788*(cos(11*p - 11*t) - 11*cos(11*p - 9*t) + 55*cos(11*p - 7*t) - 165*cos(11*p - 5*t) + 330*cos(11*p - 3*t) - 462*cos(11*p - t) + 462*cos(11*p + t) - 330*cos(11*p + 3*t) + 165*cos(11*p + 5*t) - 55*cos(11*p + 7*t) + 11*cos(11*p + 9*t) - cos(11*p + 11*t)); break;
      case -10 : return 0.001270685283012930339479964681051*(sin(10*p - 11*t) - 9*sin(10*p - 9*t) + 35*sin(10*p - 7*t) - 75*sin(10*p - 5*t) + 90*sin(10*p - 3*t) - 42*sin(10*p - t) - 42*sin(10*p + t) + 90*sin(10*p + 3*t) - 75*sin(10*p + 5*t) + 35*sin(10*p + 7*t) - 9*sin(10*p + 9*t) + sin(10*p + 11*t)); break;
      case -9 : return -0.0001960709959164239929600350650055*(21*cos(9*p - 11*t) - 151*cos(9*p - 9*t) + 435*cos(9*p - 7*t) - 585*cos(9*p - 5*t) + 210*cos(9*p - 3*t) + 378*cos(9*p - t) - 378*cos(9*p + t) - 210*cos(9*p + 3*t) + 585*cos(9*p + 5*t) - 435*cos(9*p + 7*t) + 151*cos(9*p + 9*t) - 21*cos(9*p + 11*t)); break;
      case -8 : return -0.00151875940371722513571203636154*(7*sin(8*p - 11*t) - 39*sin(8*p - 9*t) + 77*sin(8*p - 7*t) - 45*sin(8*p - 5*t) - 42*sin(8*p - 3*t) + 42*sin(8*p - t) + 42*sin(8*p + t) - 42*sin(8*p + 3*t) - 45*sin(8*p + 5*t) + 77*sin(8*p + 7*t) - 39*sin(8*p + 9*t) + 7*sin(8*p + 11*t)); break;
      case -7 : return 0.0001742136515883019988618247258659*(133*cos(7*p - 11*t) - 551*cos(7*p - 9*t) + 643*cos(7*p - 7*t) + 135*cos(7*p - 5*t) - 462*cos(7*p - 3*t) - 294*cos(7*p - t) + 294*cos(7*p + t) + 462*cos(7*p + 3*t) - 135*cos(7*p + 5*t) - 643*cos(7*p + 7*t) + 551*cos(7*p + 9*t) - 133*cos(7*p + 11*t)); break;
      case -6 : return 0.0001101823877028089813227240054718*(399*sin(6*p - 11*t) - 1159*sin(6*p - 9*t) + 525*sin(6*p - 7*t) + 795*sin(6*p - 5*t) + 70*sin(6*p - 3*t) - 630*sin(6*p - t) - 630*sin(6*p + t) + 70*sin(6*p + 3*t) + 795*sin(6*p + 5*t) + 525*sin(6*p + 7*t) - 1159*sin(6*p + 9*t) + 399*sin(6*p + 11*t)); break;
      case -5 : return -0.0000327290461389718940228154423877*(2261*cos(5*p - 11*t) - 4199*cos(5*p - 9*t) - 765*cos(5*p - 7*t) + 1975*cos(5*p - 5*t) + 2450*cos(5*p - 3*t) + 1050*cos(5*p - t) - 1050*cos(5*p + t) - 2450*cos(5*p + 3*t) - 1975*cos(5*p + 5*t) + 765*cos(5*p + 7*t) + 4199*cos(5*p + 9*t) - 2261*cos(5*p + 11*t)); break;
      case -4 : return -0.000346371666928313460878751369526*(323*sin(4*p - 11*t) - 323*sin(4*p - 9*t) - 255*sin(4*p - 7*t) - 65*sin(4*p - 5*t) + 110*sin(4*p - 3*t) + 210*sin(4*p - t) + 210*sin(4*p + t) + 110*sin(4*p + 3*t) - 65*sin(4*p + 5*t) - 255*sin(4*p + 7*t) - 323*sin(4*p + 9*t) + 323*sin(4*p + 11*t)); break;
      case -3 : return 0.000158096312714419511069847760972*(969*cos(3*p - 11*t) - 323*cos(3*p - 9*t) - 561*cos(3*p - 7*t) - 525*cos(3*p - 5*t) - 358*cos(3*p - 3*t) - 126*cos(3*p - t) + 126*cos(3*p + t) + 358*cos(3*p + 3*t) + 525*cos(3*p + 5*t) + 561*cos(3*p + 7*t) + 323*cos(3*p + 9*t) - 969*cos(3*p + 11*t)); break;
      case -2 : return 0.0000845060337556615359449834435337*(2261*sin(2*p - 11*t) + 323*sin(2*p - 9*t) - 289*sin(2*p - 7*t) - 615*sin(2*p - 5*t) - 798*sin(2*p - 3*t) - 882*sin(2*p - t) - 882*sin(2*p + t) - 798*sin(2*p + 3*t) - 615*sin(2*p + 5*t) - 289*sin(2*p + 7*t) + 323*sin(2*p + 9*t) + 2261*sin(2*p + 11*t)); break;
      case -1 : return -7.41166945852334579335102325981e-6*(29393*cos(p - 11*t) + 12597*cos(p - 9*t) + 7735*cos(p - 7*t) + 4875*cos(p - 5*t) + 2730*cos(p - 3*t) + 882*cos(p - t) - 882*cos(p + t) - 2730*cos(p + 3*t) - 4875*cos(p + 5*t) - 7735*cos(p + 7*t) - 12597*cos(p + 9*t) - 29393*cos(p + 11*t)); break;
      case 0 : return 0; break;
      case 1 : return -7.41166945852334579335102325981e-6*(29393*cos(p - 11*t) + 12597*cos(p - 9*t) + 7735*cos(p - 7*t) + 4875*cos(p - 5*t) + 2730*cos(p - 3*t) + 882*cos(p - t) - 882*cos(p + t) - 2730*cos(p + 3*t) - 4875*cos(p + 5*t) - 7735*cos(p + 7*t) - 12597*cos(p + 9*t) - 29393*cos(p + 11*t)); break;
      case 2 : return -0.0000845060337556615359449834435337*(2261*sin(2*p - 11*t) + 323*sin(2*p - 9*t) - 289*sin(2*p - 7*t) - 615*sin(2*p - 5*t) - 798*sin(2*p - 3*t) - 882*sin(2*p - t) - 882*sin(2*p + t) - 798*sin(2*p + 3*t) - 615*sin(2*p + 5*t) - 289*sin(2*p + 7*t) + 323*sin(2*p + 9*t) + 2261*sin(2*p + 11*t)); break;
      case 3 : return 0.000158096312714419511069847760972*(969*cos(3*p - 11*t) - 323*cos(3*p - 9*t) - 561*cos(3*p - 7*t) - 525*cos(3*p - 5*t) - 358*cos(3*p - 3*t) - 126*cos(3*p - t) + 126*cos(3*p + t) + 358*cos(3*p + 3*t) + 525*cos(3*p + 5*t) + 561*cos(3*p + 7*t) + 323*cos(3*p + 9*t) - 969*cos(3*p + 11*t)); break;
      case 4 : return 0.000346371666928313460878751369526*(323*sin(4*p - 11*t) - 323*sin(4*p - 9*t) - 255*sin(4*p - 7*t) - 65*sin(4*p - 5*t) + 110*sin(4*p - 3*t) + 210*sin(4*p - t) + 210*sin(4*p + t) + 110*sin(4*p + 3*t) - 65*sin(4*p + 5*t) - 255*sin(4*p + 7*t) - 323*sin(4*p + 9*t) + 323*sin(4*p + 11*t)); break;
      case 5 : return -0.0000327290461389718940228154423877*(2261*cos(5*p - 11*t) - 4199*cos(5*p - 9*t) - 765*cos(5*p - 7*t) + 1975*cos(5*p - 5*t) + 2450*cos(5*p - 3*t) + 1050*cos(5*p - t) - 1050*cos(5*p + t) - 2450*cos(5*p + 3*t) - 1975*cos(5*p + 5*t) + 765*cos(5*p + 7*t) + 4199*cos(5*p + 9*t) - 2261*cos(5*p + 11*t)); break;
      case 6 : return -0.0001101823877028089813227240054718*(399*sin(6*p - 11*t) - 1159*sin(6*p - 9*t) + 525*sin(6*p - 7*t) + 795*sin(6*p - 5*t) + 70*sin(6*p - 3*t) - 630*sin(6*p - t) - 630*sin(6*p + t) + 70*sin(6*p + 3*t) + 795*sin(6*p + 5*t) + 525*sin(6*p + 7*t) - 1159*sin(6*p + 9*t) + 399*sin(6*p + 11*t)); break;
      case 7 : return 0.0001742136515883019988618247258659*(133*cos(7*p - 11*t) - 551*cos(7*p - 9*t) + 643*cos(7*p - 7*t) + 135*cos(7*p - 5*t) - 462*cos(7*p - 3*t) - 294*cos(7*p - t) + 294*cos(7*p + t) + 462*cos(7*p + 3*t) - 135*cos(7*p + 5*t) - 643*cos(7*p + 7*t) + 551*cos(7*p + 9*t) - 133*cos(7*p + 11*t)); break;
      case 8 : return 0.00151875940371722513571203636154*(7*sin(8*p - 11*t) - 39*sin(8*p - 9*t) + 77*sin(8*p - 7*t) - 45*sin(8*p - 5*t) - 42*sin(8*p - 3*t) + 42*sin(8*p - t) + 42*sin(8*p + t) - 42*sin(8*p + 3*t) - 45*sin(8*p + 5*t) + 77*sin(8*p + 7*t) - 39*sin(8*p + 9*t) + 7*sin(8*p + 11*t)); break;
      case 9 : return -0.0001960709959164239929600350650055*(21*cos(9*p - 11*t) - 151*cos(9*p - 9*t) + 435*cos(9*p - 7*t) - 585*cos(9*p - 5*t) + 210*cos(9*p - 3*t) + 378*cos(9*p - t) - 378*cos(9*p + t) - 210*cos(9*p + 3*t) + 585*cos(9*p + 5*t) - 435*cos(9*p + 7*t) + 151*cos(9*p + 9*t) - 21*cos(9*p + 11*t)); break;
      case 10 : return -0.001270685283012930339479964681051*(sin(10*p - 11*t) - 9*sin(10*p - 9*t) + 35*sin(10*p - 7*t) - 75*sin(10*p - 5*t) + 90*sin(10*p - 3*t) - 42*sin(10*p - t) - 42*sin(10*p + t) + 90*sin(10*p + 3*t) - 75*sin(10*p + 5*t) + 35*sin(10*p + 7*t) - 9*sin(10*p + 9*t) + sin(10*p + 11*t)); break;
      case 11 : return 0.000270911012600888330984112957788*(cos(11*p - 11*t) - 11*cos(11*p - 9*t) + 55*cos(11*p - 7*t) - 165*cos(11*p - 5*t) + 330*cos(11*p - 3*t) - 462*cos(11*p - t) + 462*cos(11*p + t) - 330*cos(11*p + 3*t) + 165*cos(11*p + 5*t) - 55*cos(11*p + 7*t) + 11*cos(11*p + 9*t) - cos(11*p + 11*t)); break;
    }
    case 12 : switch (m) 
    {
      case -12 : return -0.000138248697202683391346364508269*(924*sin(12*p) + sin(12*p - 12*t) - 12*sin(12*p - 10*t) + 66*sin(12*p - 8*t) - 220*sin(12*p - 6*t) + 495*sin(12*p - 4*t) - 792*sin(12*p - 2*t) - 792*sin(12*p + 2*t) + 495*sin(12*p + 4*t) - 220*sin(12*p + 6*t) + 66*sin(12*p + 8*t) - 12*sin(12*p + 10*t) + sin(12*p + 12*t)); break;
      case -11 : return 0.000677277531502220827460282394471*(cos(11*p - 12*t) - 10*cos(11*p - 10*t) + 44*cos(11*p - 8*t) - 110*cos(11*p - 6*t) + 165*cos(11*p - 4*t) - 132*cos(11*p - 2*t) + 132*cos(11*p + 2*t) - 165*cos(11*p + 4*t) + 110*cos(11*p + 6*t) - 44*cos(11*p + 8*t) + 10*cos(11*p + 10*t) - cos(11*p + 12*t)); break;
      case -10 : return -0.0000998591241044474053298361446224*(924*sin(10*p) - 23*sin(10*p - 12*t) + 188*sin(10*p - 10*t) - 638*sin(10*p - 8*t) + 1100*sin(10*p - 6*t) - 825*sin(10*p - 4*t) - 264*sin(10*p - 2*t) - 264*sin(10*p + 2*t) - 825*sin(10*p + 4*t) + 1100*sin(10*p + 6*t) - 638*sin(10*p + 8*t) + 188*sin(10*p + 10*t) - 23*sin(10*p + 12*t)); break;
      case -9 : return -0.000270419786425946424219471915251*(23*cos(9*p - 12*t) - 150*cos(9*p - 10*t) + 372*cos(9*p - 8*t) - 370*cos(9*p - 6*t) - 45*cos(9*p - 4*t) + 324*cos(9*p - 2*t) - 324*cos(9*p + 2*t) + 45*cos(9*p + 4*t) + 370*cos(9*p + 6*t) - 372*cos(9*p + 8*t) + 150*cos(9*p + 10*t) - 23*cos(9*p + 12*t)); break;
      case -8 : return -0.0000885156529079065189295199903972*(924*sin(8*p) + 161*sin(8*p - 12*t) - 812*sin(8*p - 10*t) + 1346*sin(8*p - 8*t) - 380*sin(8*p - 6*t) - 945*sin(8*p - 4*t) + 168*sin(8*p - 2*t) + 168*sin(8*p + 2*t) - 945*sin(8*p + 4*t) - 380*sin(8*p + 6*t) + 1346*sin(8*p + 8*t) - 812*sin(8*p + 10*t) + 161*sin(8*p + 12*t)); break;
      case -7 : return 0.0001770313058158130378590399807944*(161*cos(7*p - 12*t) - 602*cos(7*p - 10*t) + 556*cos(7*p - 8*t) + 290*cos(7*p - 6*t) - 315*cos(7*p - 4*t) - 420*cos(7*p - 2*t) + 420*cos(7*p + 2*t) + 315*cos(7*p + 4*t) - 290*cos(7*p + 6*t) - 556*cos(7*p + 8*t) + 602*cos(7*p + 10*t) - 161*cos(7*p + 12*t)); break;
      case -6 : return -0.00001658050092328703130490742729138*(4620*sin(6*p) - 3059*sin(6*p - 12*t) + 7980*sin(6*p - 10*t) - 2166*sin(6*p - 8*t) - 5380*sin(6*p - 6*t) - 2205*sin(6*p - 4*t) + 2520*sin(6*p - 2*t) + 2520*sin(6*p + 2*t) - 2205*sin(6*p + 4*t) - 5380*sin(6*p + 6*t) - 2166*sin(6*p + 8*t) + 7980*sin(6*p + 10*t) - 3059*sin(6*p + 12*t)); break;
      case -5 : return -0.0001861156612680872278147844505806*(437*cos(5*p - 12*t) - 722*cos(5*p - 10*t) - 228*cos(5*p - 8*t) + 250*cos(5*p - 6*t) + 425*cos(5*p - 4*t) + 300*cos(5*p - 2*t) - 300*cos(5*p + 2*t) - 425*cos(5*p + 4*t) - 250*cos(5*p + 6*t) + 228*cos(5*p + 8*t) + 722*cos(5*p + 10*t) - 437*cos(5*p + 12*t)); break;
      case -4 : return -0.00001595928629016969586604110799767*(4620*sin(4*p) + 7429*sin(4*p - 12*t) - 6460*sin(4*p - 10*t) - 5814*sin(4*p - 8*t) - 2380*sin(4*p - 6*t) + 1195*sin(4*p - 4*t) + 3720*sin(4*p - 2*t) + 3720*sin(4*p + 2*t) + 1195*sin(4*p + 4*t) - 2380*sin(4*p + 6*t) - 5814*sin(4*p + 8*t) - 6460*sin(4*p + 10*t) + 7429*sin(4*p + 12*t)); break;
      case -3 : return 0.0000212790483868929278213881439969*(7429*cos(3*p - 12*t) - 1938*cos(3*p - 10*t) - 3876*cos(3*p - 8*t) - 3910*cos(3*p - 6*t) - 3015*cos(3*p - 4*t) - 1620*cos(3*p - 2*t) + 1620*cos(3*p + 2*t) + 3015*cos(3*p + 4*t) + 3910*cos(3*p + 6*t) + 3876*cos(3*p + 8*t) + 1938*cos(3*p + 10*t) - 7429*cos(3*p + 12*t)); break;
      case -2 : return -0.0000260614053799405792989657278166*(2772*sin(2*p) - 7429*sin(2*p - 12*t) - 1292*sin(2*p - 10*t) + 646*sin(2*p - 8*t) + 1700*sin(2*p - 6*t) + 2325*sin(2*p - 4*t) + 2664*sin(2*p - 2*t) + 2664*sin(2*p + 2*t) + 2325*sin(2*p + 4*t) + 1700*sin(2*p + 6*t) + 646*sin(2*p + 8*t) - 1292*sin(2*p + 10*t) - 7429*sin(2*p + 12*t)); break;
      case -1 : return -4.20017578598614200022888376552e-6*(52003*cos(p - 12*t) + 22610*cos(p - 10*t) + 14212*cos(p - 8*t) + 9350*cos(p - 6*t) + 5775*cos(p - 4*t) + 2772*cos(p - 2*t) - 2772*cos(p + 2*t) - 5775*cos(p + 4*t) - 9350*cos(p + 6*t) - 14212*cos(p + 8*t) - 22610*cos(p + 10*t) - 52003*cos(p + 12*t)); break;
      case 0 : return 0; break;
      case 1 : return -4.20017578598614200022888376552e-6*(52003*cos(p - 12*t) + 22610*cos(p - 10*t) + 14212*cos(p - 8*t) + 9350*cos(p - 6*t) + 5775*cos(p - 4*t) + 2772*cos(p - 2*t) - 2772*cos(p + 2*t) - 5775*cos(p + 4*t) - 9350*cos(p + 6*t) - 14212*cos(p + 8*t) - 22610*cos(p + 10*t) - 52003*cos(p + 12*t)); break;
      case 2 : return 0.0000260614053799405792989657278166*(2772*sin(2*p) - 7429*sin(2*p - 12*t) - 1292*sin(2*p - 10*t) + 646*sin(2*p - 8*t) + 1700*sin(2*p - 6*t) + 2325*sin(2*p - 4*t) + 2664*sin(2*p - 2*t) + 2664*sin(2*p + 2*t) + 2325*sin(2*p + 4*t) + 1700*sin(2*p + 6*t) + 646*sin(2*p + 8*t) - 1292*sin(2*p + 10*t) - 7429*sin(2*p + 12*t)); break;
      case 3 : return 0.0000212790483868929278213881439969*(7429*cos(3*p - 12*t) - 1938*cos(3*p - 10*t) - 3876*cos(3*p - 8*t) - 3910*cos(3*p - 6*t) - 3015*cos(3*p - 4*t) - 1620*cos(3*p - 2*t) + 1620*cos(3*p + 2*t) + 3015*cos(3*p + 4*t) + 3910*cos(3*p + 6*t) + 3876*cos(3*p + 8*t) + 1938*cos(3*p + 10*t) - 7429*cos(3*p + 12*t)); break;
      case 4 : return 0.00001595928629016969586604110799767*(4620*sin(4*p) + 7429*sin(4*p - 12*t) - 6460*sin(4*p - 10*t) - 5814*sin(4*p - 8*t) - 2380*sin(4*p - 6*t) + 1195*sin(4*p - 4*t) + 3720*sin(4*p - 2*t) + 3720*sin(4*p + 2*t) + 1195*sin(4*p + 4*t) - 2380*sin(4*p + 6*t) - 5814*sin(4*p + 8*t) - 6460*sin(4*p + 10*t) + 7429*sin(4*p + 12*t)); break;
      case 5 : return -0.0001861156612680872278147844505806*(437*cos(5*p - 12*t) - 722*cos(5*p - 10*t) - 228*cos(5*p - 8*t) + 250*cos(5*p - 6*t) + 425*cos(5*p - 4*t) + 300*cos(5*p - 2*t) - 300*cos(5*p + 2*t) - 425*cos(5*p + 4*t) - 250*cos(5*p + 6*t) + 228*cos(5*p + 8*t) + 722*cos(5*p + 10*t) - 437*cos(5*p + 12*t)); break;
      case 6 : return 0.00001658050092328703130490742729138*(4620*sin(6*p) - 3059*sin(6*p - 12*t) + 7980*sin(6*p - 10*t) - 2166*sin(6*p - 8*t) - 5380*sin(6*p - 6*t) - 2205*sin(6*p - 4*t) + 2520*sin(6*p - 2*t) + 2520*sin(6*p + 2*t) - 2205*sin(6*p + 4*t) - 5380*sin(6*p + 6*t) - 2166*sin(6*p + 8*t) + 7980*sin(6*p + 10*t) - 3059*sin(6*p + 12*t)); break;
      case 7 : return 0.0001770313058158130378590399807944*(161*cos(7*p - 12*t) - 602*cos(7*p - 10*t) + 556*cos(7*p - 8*t) + 290*cos(7*p - 6*t) - 315*cos(7*p - 4*t) - 420*cos(7*p - 2*t) + 420*cos(7*p + 2*t) + 315*cos(7*p + 4*t) - 290*cos(7*p + 6*t) - 556*cos(7*p + 8*t) + 602*cos(7*p + 10*t) - 161*cos(7*p + 12*t)); break;
      case 8 : return 0.0000885156529079065189295199903972*(924*sin(8*p) + 161*sin(8*p - 12*t) - 812*sin(8*p - 10*t) + 1346*sin(8*p - 8*t) - 380*sin(8*p - 6*t) - 945*sin(8*p - 4*t) + 168*sin(8*p - 2*t) + 168*sin(8*p + 2*t) - 945*sin(8*p + 4*t) - 380*sin(8*p + 6*t) + 1346*sin(8*p + 8*t) - 812*sin(8*p + 10*t) + 161*sin(8*p + 12*t)); break;
      case 9 : return -0.000270419786425946424219471915251*(23*cos(9*p - 12*t) - 150*cos(9*p - 10*t) + 372*cos(9*p - 8*t) - 370*cos(9*p - 6*t) - 45*cos(9*p - 4*t) + 324*cos(9*p - 2*t) - 324*cos(9*p + 2*t) + 45*cos(9*p + 4*t) + 370*cos(9*p + 6*t) - 372*cos(9*p + 8*t) + 150*cos(9*p + 10*t) - 23*cos(9*p + 12*t)); break;
      case 10 : return 0.0000998591241044474053298361446224*(924*sin(10*p) - 23*sin(10*p - 12*t) + 188*sin(10*p - 10*t) - 638*sin(10*p - 8*t) + 1100*sin(10*p - 6*t) - 825*sin(10*p - 4*t) - 264*sin(10*p - 2*t) - 264*sin(10*p + 2*t) - 825*sin(10*p + 4*t) + 1100*sin(10*p + 6*t) - 638*sin(10*p + 8*t) + 188*sin(10*p + 10*t) - 23*sin(10*p + 12*t)); break;
      case 11 : return 0.000677277531502220827460282394471*(cos(11*p - 12*t) - 10*cos(11*p - 10*t) + 44*cos(11*p - 8*t) - 110*cos(11*p - 6*t) + 165*cos(11*p - 4*t) - 132*cos(11*p - 2*t) + 132*cos(11*p + 2*t) - 165*cos(11*p + 4*t) + 110*cos(11*p + 6*t) - 44*cos(11*p + 8*t) + 10*cos(11*p + 10*t) - cos(11*p + 12*t)); break;
      case 12 : return 0.000138248697202683391346364508269*(924*sin(12*p) + sin(12*p - 12*t) - 12*sin(12*p - 10*t) + 66*sin(12*p - 8*t) - 220*sin(12*p - 6*t) + 495*sin(12*p - 4*t) - 792*sin(12*p - 2*t) - 792*sin(12*p + 2*t) + 495*sin(12*p + 4*t) - 220*sin(12*p + 6*t) + 66*sin(12*p + 8*t) - 12*sin(12*p + 10*t) + sin(12*p + 12*t)); break;
    }
  }
  return 0;
}

/*  The CubicHarmonicK are defined as the normalized cubic harmonics
*
*   For a derivation and definition see
*   Journal of Applied Mathematics and Physics (ZAMP) Vol. 23, 1972
*   Cubic Harmonics as Linear Combinations of Spherical Harmonics
*   Juerg Muggli
*
*   The order of m is realively arbetrary and the choise made here is according to the shonwflies notation of the irrep.
*   i.e. lowest degneracy first, t1 before t2, a before b
*   the order is:
*   k, m   irrep.   pre-factor     wave-function
*   0, 0   a1_g     sqrt(2k+1/4Pi) (1)
*   1,-1   t1_u     sqrt(2k+1/4Pi) (x)
*   1, 0   t1_u     sqrt(2k+1/4Pi) (y)
*   1, 1   t1_u     sqrt(2k+1/4Pi) (z)
*   2,-2   e_g      sqrt(2k+1/4Pi) (sqrt(3/4)(x^2-y^2))
*   2,-1   e_g      sqrt(2k+1/4Pi) (sqrt(1/4)(3z^2-1))
*   2, 0   t2_g     sqrt(2k+1/4Pi) (sqrt(3)yz)
*   2, 1   t2_g     sqrt(2k+1/4Pi) (sqrt(3)xz)
*   2, 2   t2_g     sqrt(2k+1/4Pi) (sqrt(3)xy)
*   3,-3   a2_u     sqrt(2k+1/4Pi) (sqrt(15)xyz)
*   3,-2   t1_u     sqrt(2k+1/4Pi) (sqrt(1/4)x(5x^2-3))
*   3,-1   t1_u     sqrt(2k+1/4Pi) (sqrt(1/4)y(5y^2-3))
*   3, 0   t1_u     sqrt(2k+1/4Pi) (sqrt(1/4)z(5z^2-3))
*   3, 1   t2_u     sqrt(2k+1/4Pi) (sqrt(15/4)x(y^2-z^2))
*   3, 2   t2_u     sqrt(2k+1/4Pi) (sqrt(15/4)y(z^2-x^2))
*   3, 3   t2_u     sqrt(2k+1/4Pi) (sqrt(15/4)z(x^2-y^2))
*   4,-4   a1_g     sqrt(2k+1/4Pi) (sqrt(7/48)(5(x^4+y^4+z^4)-3)
*   4,-3   e_g      sqrt(2k+1/4Pi) (sqrt(5/16)(x^2-y^2)(7z^2-1))
*   4,-2   e_g      sqrt(2k+1/4Pi) (sqrt(5/48)(1-14x^2y^2+4z^2-7z^4))
*   4,-1   t1_g     sqrt(2k+1/4Pi) (sqrt(35/4)yz(y^2-z^2))
*   4, 0   t1_g     sqrt(2k+1/4Pi) (sqrt(35/4)xz(z^2-x^2))
*   4, 1   t1_g     sqrt(2k+1/4Pi) (sqrt(35/4)xy(x^2-y^2))
*   4, 2   t2_g     sqrt(2k+1/4Pi) (sqrt(5/4)yz(7x^2-1))
*   4, 3   t2_g     sqrt(2k+1/4Pi) (sqrt(5/4)xz(7y^2-1))
*   4, 4   t2_g     sqrt(2k+1/4Pi) (sqrt(5/4)xy(7z^2-1))
*   5,-5   e_u      sqrt(2k+1/4Pi) (sqrt(105/4)xyz(3z^2-1))
*   5,-4   e_u      sqrt(2k+1/4Pi) (sqrt(315/4)xyz(x^2-y^2))
*   5,-3   t1_u^A   sqrt(2k+1/4Pi) (sqrt(1/64)x(15-70x^2+63x^4))
*   5,-2   t1_u^A   sqrt(2k+1/4Pi) (sqrt(1/64)y(15-70y^2+63y^4))
*   5,-1   t1_u^A   sqrt(2k+1/4Pi) (sqrt(1/64)z(15-70z^2+63z^4))
*   5, 0   t1_u^B   sqrt(2k+1/4Pi) (sqrt(315/64)x(y^4-6y^2z^2-z^4))
*   5, 1   t1_u^B   sqrt(2k+1/4Pi) (sqrt(315/64)y(z^4-6z^2x^2-x^4))
*   5, 2   t1_u^B   sqrt(2k+1/4Pi) (sqrt(315/64)z(x^4-6x^2y^2-y^4))
*   5, 3   t2_u     sqrt(2k+1/4Pi) (sqrt(105/4)(x^2-y^2)z(3z^2-1)
*   5, 4   t2_u     sqrt(2k+1/4Pi) (sqrt(105/4)(y^2-z^2)x(3x^2-1)
*   5, 5   t2_u     sqrt(2k+1/4Pi) (sqrt(105/4)(z^2-x^2)y(3y^2-1)
*   6,-6   a1_g     sqrt(2k+1/4Pi) (sqrt(1/32)(7(x^6+y^6+z^6)+210x^2y^2z^2-5))
*   6,-5   a2_g     sqrt(2k+1/4Pi) (sqrt(1155/32)(x^2-y^2)(y^2-z^2)(z^2-x^2)
*   6,-4   e_g      sqrt(2k+1/4Pi) (sqrt(21/32)(x^2-y^2)(x^4+y^4-5y^2z^2+5z^4-x^2(9y^2+5z^2)))
*   6,-3   e_g      sqrt(2k+1/4Pi) (sqrt(7/32)(-x^6-y^6+15x^4z^2+15y^4z^2-15x^2z^4-15y^2z^4+2z^6))
*   6,-2   t1_g     sqrt(2k+1/4Pi) (sqrt(63/16)(11x^2-1)yz(y^2-z^2))
*   6,-1   t1_g     sqrt(2k+1/4Pi) (sqrt(63/16)(11y^2-1)xz(z^2-x^2))
*   6, 0   t1_g     sqrt(2k+1/4Pi) (sqrt(63/16)(11z^2-1)xy(x^2-y^2))
*   6, 1   t2_g^A   sqrt(2k+1/4Pi) (sqrt(105/128)yz(1-18x^2+33x^4))
*   6, 2   t2_g^A   sqrt(2k+1/4Pi) (sqrt(105/128)xz(1-18y^2+33y^4))
*   6, 3   t2_g^A   sqrt(2k+1/4Pi) (sqrt(105/128)xy(1-18z^2+33z^4))
*   6, 4   t2_g^B   sqrt(2k+1/4Pi) (sqrt(231/128)yz(3y^4-10y^2z^2+3z^4))
*   6, 5   t2_g^B   sqrt(2k+1/4Pi) (sqrt(231/128)xz(3z^4-10y^2x^2+3x^4))
*   6, 6   t2_g^B   sqrt(2k+1/4Pi) (sqrt(231/128)xy(3x^4-10y^2y^2+3y^4))
*   
*/
double CubicHarmonicK(const int k, const int m, const double t, const double p)
{
  switch (k)
  {
    case 0 : switch (m) 
    {
      case 0 : return 0.28209479177387814347403972578; break;
    }
    case 1 : switch (m) 
    {
      case -1 : return -0.244301255951459960793192311419*(sin(p - t) - sin(p + t)); break;
      case 0 : return 0.244301255951459960793192311419*(cos(p - t) - cos(p + t)); break;
      case 1 : return 0.488602511902919921586384622838*cos(t); break;
    }
    case 2 : switch (m) 
    {
      case -2 : return 0.1365685538240098838179232132253*(2*cos(2*p) - cos(2*p - 2*t) - cos(2*p + 2*t)); break;
      case -1 : return 0.1576957826262600030154468451479*(1 + 3*cos(2*t)); break;
      case 0 : return 0.273137107648019767635846426451*(cos(p - 2*t) - cos(p + 2*t)); break;
      case 1 : return -0.273137107648019767635846426451*(sin(p - 2*t) - sin(p + 2*t)); break;
      case 2 : return 0.1365685538240098838179232132253*(2*sin(2*p) - sin(2*p - 2*t) - sin(2*p + 2*t)); break;
    }
    case 3 : switch (m) 
    {
      case -3 : return -0.1806632151650346284618362596499*(sin(2*p - 3*t) - sin(2*p - t) - sin(2*p + t) + sin(2*p + 3*t)); break;
      case -2 : return 0.01166176039344110598169987228747*(15*sin(p - 3*t) + 5*sin(3*p - 3*t) + 3*sin(p - t) - 15*sin(3*p - t) - 3*sin(p + t) + 15*sin(3*p + t) - 15*sin(p + 3*t) - 5*sin(3*p + 3*t)); break;
      case -1 : return -0.01166176039344110598169987228747*(15*cos(p - 3*t) - 5*cos(3*p - 3*t) + 3*cos(p - t) + 15*cos(3*p - t) - 3*cos(p + t) - 15*cos(3*p + t) - 15*cos(p + 3*t) + 5*cos(3*p + 3*t)); break;
      case 0 : return 0.0932940831475288478535989782997*(3*cos(t) + 5*cos(3*t)); break;
      case 1 : return 0.0451658037912586571154590649125*(5*sin(p - 3*t) - sin(3*p - 3*t) + sin(p - t) + 3*sin(3*p - t) - sin(p + t) - 3*sin(3*p + t) - 5*sin(p + 3*t) + sin(3*p + 3*t)); break;
      case 2 : return 0.0451658037912586571154590649125*(5*cos(p - 3*t) + cos(3*p - 3*t) + cos(p - t) - 3*cos(3*p - t) - cos(p + t) + 3*cos(3*p + t) - 5*cos(p + 3*t) - cos(3*p + 3*t)); break;
      case 3 : return -0.1806632151650346284618362596499*(cos(2*p - 3*t) - cos(2*p - t) - cos(2*p + t) + cos(2*p + 3*t)); break;
    }
    case 4 : switch (m) 
    {
      case -4 : return 0.00504969037678360395324053380204*(18 + 30*cos(4*p) + 5*cos(4*p - 4*t) - 20*cos(4*p - 2*t) + 40*cos(2*t) + 70*cos(4*t) - 20*cos(4*p + 2*t) + 5*cos(4*p + 4*t)); break;
      case -3 : return 0.0295679592424237505653962834652*(6*cos(2*p) - 7*cos(2*p - 4*t) + 4*cos(2*p - 2*t) + 4*cos(2*p + 2*t) - 7*cos(2*p + 4*t)); break;
      case -2 : return 0.00426776730700030886931010041329*(-18 + 42*cos(4*p) + 7*cos(4*p - 4*t) - 28*cos(4*p - 2*t) - 40*cos(2*t) - 70*cos(4*t) - 28*cos(4*p + 2*t) + 7*cos(4*p + 4*t)); break;
      case -1 : return -0.0391147334655735084116650327253*(7*cos(p - 4*t) - cos(3*p - 4*t) + 2*cos(p - 2*t) + 2*cos(3*p - 2*t) - 2*cos(p + 2*t) - 2*cos(3*p + 2*t) - 7*cos(p + 4*t) + cos(3*p + 4*t)); break;
      case 0 : return -0.0391147334655735084116650327253*(7*sin(p - 4*t) + sin(3*p - 4*t) + 2*sin(p - 2*t) - 2*sin(3*p - 2*t) - 2*sin(p + 2*t) + 2*sin(3*p + 2*t) - 7*sin(p + 4*t) - sin(3*p + 4*t)); break;
      case 1 : return 0.0391147334655735084116650327253*(6*sin(4*p) + sin(4*p - 4*t) - 4*sin(4*p - 2*t) - 4*sin(4*p + 2*t) + sin(4*p + 4*t)); break;
      case 2 : return -0.01478397962121187528269814173261*(7*cos(p - 4*t) + 7*cos(3*p - 4*t) + 2*cos(p - 2*t) - 14*cos(3*p - 2*t) - 2*cos(p + 2*t) + 14*cos(3*p + 2*t) - 7*cos(p + 4*t) - 7*cos(3*p + 4*t)); break;
      case 3 : return 0.01478397962121187528269814173261*(7*sin(p - 4*t) - 7*sin(3*p - 4*t) + 2*sin(p - 2*t) + 14*sin(3*p - 2*t) - 2*sin(p + 2*t) - 14*sin(3*p + 2*t) - 7*sin(p + 4*t) + 7*sin(3*p + 4*t)); break;
      case 4 : return 0.0295679592424237505653962834652*(6*sin(2*p) - 7*sin(2*p - 4*t) + 4*sin(2*p - 2*t) + 4*sin(2*p + 2*t) - 7*sin(2*p + 4*t)); break;
    }
    case 5 : switch (m) 
    {
      case -5 : return -0.0748990122652081836719047053057*(3*sin(2*p - 5*t) - sin(2*p - 3*t) - 2*sin(2*p - t) - 2*sin(2*p + t) - sin(2*p + 3*t) + 3*sin(2*p + 5*t)); break;
      case -4 : return 0.0648644473400325399686187038311*(sin(4*p - 5*t) - 3*sin(4*p - 3*t) + 2*sin(4*p - t) + 2*sin(4*p + t) - 3*sin(4*p + 3*t) + sin(4*p + 5*t)); break;
      case -3 : return -0.000228418598541842961796318734551*(630*sin(p - 5*t) + 315*sin(3*p - 5*t) + 63*sin(5*p - 5*t) + 210*sin(p - 3*t) - 455*sin(3*p - 3*t) - 315*sin(5*p - 3*t) + 60*sin(p - t) - 210*sin(3*p - t) + 630*sin(5*p - t) - 60*sin(p + t) + 210*sin(3*p + t) - 630*sin(5*p + t) - 210*sin(p + 3*t) + 455*sin(3*p + 3*t) + 315*sin(5*p + 3*t) - 630*sin(p + 5*t) - 315*sin(3*p + 5*t) - 63*sin(5*p + 5*t)); break;
      case -2 : return 0.000228418598541842961796318734551*(630*cos(p - 5*t) - 315*cos(3*p - 5*t) + 63*cos(5*p - 5*t) + 210*cos(p - 3*t) + 455*cos(3*p - 3*t) - 315*cos(5*p - 3*t) + 60*cos(p - t) + 210*cos(3*p - t) + 630*cos(5*p - t) - 60*cos(p + t) - 210*cos(3*p + t) - 630*cos(5*p + t) - 210*cos(p + 3*t) - 455*cos(3*p + 3*t) + 315*cos(5*p + 3*t) - 630*cos(p + 5*t) + 315*cos(3*p + 5*t) - 63*cos(5*p + 5*t)); break;
      case -1 : return 0.00730939515333897477748219950564*(30*cos(t) + 35*cos(3*t) + 63*cos(5*t)); break;
      case 0 : return -0.00405402795875203374803866898944*(42*sin(p - 5*t) - 27*sin(3*p - 5*t) + sin(5*p - 5*t) + 14*sin(p - 3*t) + 39*sin(3*p - 3*t) - 5*sin(5*p - 3*t) + 4*sin(p - t) + 18*sin(3*p - t) + 10*sin(5*p - t) - 4*sin(p + t) - 18*sin(3*p + t) - 10*sin(5*p + t) - 14*sin(p + 3*t) - 39*sin(3*p + 3*t) + 5*sin(5*p + 3*t) - 42*sin(p + 5*t) + 27*sin(3*p + 5*t) - sin(5*p + 5*t)); break;
      case 1 : return 0.00405402795875203374803866898944*(42*cos(p - 5*t) + 27*cos(3*p - 5*t) + cos(5*p - 5*t) + 14*cos(p - 3*t) - 39*cos(3*p - 3*t) - 5*cos(5*p - 3*t) + 4*cos(p - t) - 18*cos(3*p - t) + 10*cos(5*p - t) - 4*cos(p + t) + 18*cos(3*p + t) - 10*cos(5*p + t) - 14*cos(p + 3*t) + 39*cos(3*p + 3*t) + 5*cos(5*p + 3*t) - 42*cos(p + 5*t) - 27*cos(3*p + 5*t) - cos(5*p + 5*t)); break;
      case 2 : return 0.0648644473400325399686187038311*(cos(4*p - 5*t) - 3*cos(4*p - 3*t) + 2*cos(4*p - t) + 2*cos(4*p + t) - 3*cos(4*p + 3*t) + cos(4*p + 5*t)); break;
      case 3 : return -0.0046811882665755114794940440816*(42*sin(p - 5*t) + 9*sin(3*p - 5*t) - 3*sin(5*p - 5*t) + 14*sin(p - 3*t) - 13*sin(3*p - 3*t) + 15*sin(5*p - 3*t) + 4*sin(p - t) - 6*sin(3*p - t) - 30*sin(5*p - t) - 4*sin(p + t) + 6*sin(3*p + t) + 30*sin(5*p + t) - 14*sin(p + 3*t) + 13*sin(3*p + 3*t) - 15*sin(5*p + 3*t) - 42*sin(p + 5*t) - 9*sin(3*p + 5*t) + 3*sin(5*p + 5*t)); break;
      case 4 : return -0.0046811882665755114794940440816*(42*cos(p - 5*t) - 9*cos(3*p - 5*t) - 3*cos(5*p - 5*t) + 14*cos(p - 3*t) + 13*cos(3*p - 3*t) + 15*cos(5*p - 3*t) + 4*cos(p - t) + 6*cos(3*p - t) - 30*cos(5*p - t) - 4*cos(p + t) - 6*cos(3*p + t) + 30*cos(5*p + t) - 14*cos(p + 3*t) - 13*cos(3*p + 3*t) - 15*cos(5*p + 3*t) - 42*cos(p + 5*t) + 9*cos(3*p + 5*t) + 3*cos(5*p + 5*t)); break;
      case 5 : return -0.0748990122652081836719047053057*(3*cos(2*p - 5*t) - cos(2*p - 3*t) - 2*cos(2*p - t) - 2*cos(2*p + t) - cos(2*p + 3*t) + 3*cos(2*p + 5*t)); break;
    }
    case 6 : switch (m) 
    {
      case -6 : return -0.00035117354685983839122360053676*(-100 + 420*cos(4*p) + 231*cos(4*p - 6*t) - 546*cos(4*p - 4*t) + 105*cos(4*p - 2*t) - 210*cos(2*t) - 252*cos(4*t) - 462*cos(6*t) + 105*cos(4*p + 2*t) - 546*cos(4*p + 4*t) + 231*cos(4*p + 6*t)); break;
      case -5 : return -0.00596736757951752902820692794048*(20*cos(2*p) - 20*cos(6*p) - 33*cos(2*p - 6*t) + cos(6*p - 6*t) + 6*cos(2*p - 4*t) - 6*cos(6*p - 4*t) + 17*cos(2*p - 2*t) + 15*cos(6*p - 2*t) + 17*cos(2*p + 2*t) + 15*cos(6*p + 2*t) + 6*cos(2*p + 4*t) - 6*cos(6*p + 4*t) - 33*cos(2*p + 6*t) + cos(6*p + 6*t)); break;
      case -4 : return 0.000804639680275665580866695283345*(100*cos(2*p) + 220*cos(6*p) - 165*cos(2*p - 6*t) - 11*cos(6*p - 6*t) + 30*cos(2*p - 4*t) + 66*cos(6*p - 4*t) + 85*cos(2*p - 2*t) - 165*cos(6*p - 2*t) + 85*cos(2*p + 2*t) - 165*cos(6*p + 2*t) + 30*cos(2*p + 4*t) + 66*cos(6*p + 4*t) - 165*cos(2*p + 6*t) - 11*cos(6*p + 6*t)); break;
      case -3 : return 0.000929117872015619863588195806633*(100 + 60*cos(4*p) + 33*cos(4*p - 6*t) - 78*cos(4*p - 4*t) + 15*cos(4*p - 2*t) + 210*cos(2*t) + 252*cos(4*t) + 462*cos(6*t) + 15*cos(4*p + 2*t) - 78*cos(4*p + 4*t) + 33*cos(4*p + 6*t)); break;
      case -2 : return 0.001970956643471578747087808046096*(66*cos(p - 6*t) + 55*cos(3*p - 6*t) - 11*cos(5*p - 6*t) + 24*cos(p - 4*t) - 60*cos(3*p - 4*t) + 44*cos(5*p - 4*t) + 10*cos(p - 2*t) - 45*cos(3*p - 2*t) - 55*cos(5*p - 2*t) - 10*cos(p + 2*t) + 45*cos(3*p + 2*t) + 55*cos(5*p + 2*t) - 24*cos(p + 4*t) + 60*cos(3*p + 4*t) - 44*cos(5*p + 4*t) - 66*cos(p + 6*t) - 55*cos(3*p + 6*t) + 11*cos(5*p + 6*t)); break;
      case -1 : return 0.001970956643471578747087808046096*(66*sin(p - 6*t) - 55*sin(3*p - 6*t) - 11*sin(5*p - 6*t) + 24*sin(p - 4*t) + 60*sin(3*p - 4*t) + 44*sin(5*p - 4*t) + 10*sin(p - 2*t) + 45*sin(3*p - 2*t) - 55*sin(5*p - 2*t) - 10*sin(p + 2*t) - 45*sin(3*p + 2*t) + 55*sin(5*p + 2*t) - 24*sin(p + 4*t) - 60*sin(3*p + 4*t) - 44*sin(5*p + 4*t) - 66*sin(p + 6*t) + 55*sin(3*p + 6*t) + 11*sin(5*p + 6*t)); break;
      case 0 : return 0.00788382657388631498835123218438*(20*sin(4*p) + 11*sin(4*p - 6*t) - 26*sin(4*p - 4*t) + 5*sin(4*p - 2*t) + 5*sin(4*p + 2*t) - 26*sin(4*p + 4*t) + 11*sin(4*p + 6*t)); break;
      case 1 : return 0.000899614511245042479630435347875*(66*cos(p - 6*t) + 99*cos(3*p - 6*t) + 33*cos(5*p - 6*t) + 24*cos(p - 4*t) - 108*cos(3*p - 4*t) - 132*cos(5*p - 4*t) + 10*cos(p - 2*t) - 81*cos(3*p - 2*t) + 165*cos(5*p - 2*t) - 10*cos(p + 2*t) + 81*cos(3*p + 2*t) - 165*cos(5*p + 2*t) - 24*cos(p + 4*t) + 108*cos(3*p + 4*t) + 132*cos(5*p + 4*t) - 66*cos(p + 6*t) - 99*cos(3*p + 6*t) - 33*cos(5*p + 6*t)); break;
      case 2 : return -0.000899614511245042479630435347875*(66*sin(p - 6*t) - 99*sin(3*p - 6*t) + 33*sin(5*p - 6*t) + 24*sin(p - 4*t) + 108*sin(3*p - 4*t) - 132*sin(5*p - 4*t) + 10*sin(p - 2*t) + 81*sin(3*p - 2*t) + 165*sin(5*p - 2*t) - 10*sin(p + 2*t) - 81*sin(3*p + 2*t) - 165*sin(5*p + 2*t) - 24*sin(p + 4*t) - 108*sin(3*p + 4*t) + 132*sin(5*p + 4*t) - 66*sin(p + 6*t) + 99*sin(3*p + 6*t) - 33*sin(5*p + 6*t)); break;
      case 3 : return 0.007196916089960339837043482783*(20*sin(2*p) - 33*sin(2*p - 6*t) + 6*sin(2*p - 4*t) + 17*sin(2*p - 2*t) + 17*sin(2*p + 2*t) + 6*sin(2*p + 4*t) - 33*sin(2*p + 6*t)); break;
      case 4 : return 0.001334343955452957660107045264211*(198*cos(p - 6*t) - 55*cos(3*p - 6*t) + 3*cos(5*p - 6*t) + 72*cos(p - 4*t) + 60*cos(3*p - 4*t) - 12*cos(5*p - 4*t) + 30*cos(p - 2*t) + 45*cos(3*p - 2*t) + 15*cos(5*p - 2*t) - 30*cos(p + 2*t) - 45*cos(3*p + 2*t) - 15*cos(5*p + 2*t) - 72*cos(p + 4*t) - 60*cos(3*p + 4*t) + 12*cos(5*p + 4*t) - 198*cos(p + 6*t) + 55*cos(3*p + 6*t) - 3*cos(5*p + 6*t)); break;
      case 5 : return -0.001334343955452957660107045264211*(198*sin(p - 6*t) + 55*sin(3*p - 6*t) + 3*sin(5*p - 6*t) + 72*sin(p - 4*t) - 60*sin(3*p - 4*t) - 12*sin(5*p - 4*t) + 30*sin(p - 2*t) - 45*sin(3*p - 2*t) + 15*sin(5*p - 2*t) - 30*sin(p + 2*t) + 45*sin(3*p + 2*t) - 15*sin(5*p + 2*t) - 72*sin(p + 4*t) + 60*sin(3*p + 4*t) + 12*sin(5*p + 4*t) - 198*sin(p + 6*t) - 55*sin(3*p + 6*t) - 3*sin(5*p + 6*t)); break;
      case 6 : return 0.01067475164362366128085636211369*(20*sin(6*p) - sin(6*p - 6*t) + 6*sin(6*p - 4*t) - 15*sin(6*p - 2*t) - 15*sin(6*p + 2*t) + 6*sin(6*p + 4*t) - sin(6*p + 6*t)); break;
    }
  }
  return 0;
}

