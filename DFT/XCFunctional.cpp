/*
 * XCFunctional.cpp
 *
 * Exchange-correlation functionals for the Dirac-Kohn-Sham atom (DFT stage A3, design
 * section 9), spin unpolarized. For a density n (per volume) and sigma = |grad n|^2 every
 * functional gives the energy per volume e(n, sigma) and its derivatives de/dn and
 * de/dsigma, analytically. Hartree atomic units.
 *
 *   XCLDA       Slater exchange + Perdew-Wang 1992 correlation
 *   XCLDAVWN    Slater exchange + Vosko-Wilk-Nusair correlation (fit to Ceperley-Alder,
 *               "VWN5"), the LDA of the NIST atomic reference data
 *   XCRLDAMV    as XCLDAVWN, with the MacDonald-Vosko relativistic factor on exchange:
 *               the RLDA of the NIST atomic reference data
 *   XCPBE       Perdew-Burke-Ernzerhof (PBE exchange, PBE correlation on PW92)
 *   XCPBEVWN    PBE with VWN5 instead of PW92 under the gradient correction of the
 *               correlation (as in Kraisler et al. 2010 and Lehtola 2020)
 *
 * Parameters: VWN5 from the NIST pages (paramagnetic A = 0.0310907, x0 = -0.10498,
 * b = 3.72744, c = 12.9352). PW92 (paramagnetic, A = 0.031091, alpha1 = 0.21370,
 * beta1..4 = 7.5957, 3.5876, 1.6382, 0.49294) and PBE (kappa = 0.804, beta =
 * 0.06672455060314922, mu = beta pi^2/3, gamma = (1 - ln 2)/pi^2) as in the original papers;
 * these could not be checked against an independent source in 2026-11 (see the design
 * note, section 0.2). PBE with VWN5 is checked against fully numerical atoms.
 *
 * Below the density cutoff (an argument) every functional returns 0.
 */

// 1 for the gradient-corrected functionals (they need sigma), 0 otherwise
int XCIsGGA(int method)
{
  return (method == XCPBE || method == XCPBEVWN);
}

/* VWN correlation energy per particle eps_c(rs) and deps_c/drs, for one parameter set
 * (A, x0, b, c), x = sqrt(rs):
 *   eps = A { ln(x^2/X) + 2b/Q atan(Q/(2x+b))
 *             - b x0/X(x0) [ ln((x-x0)^2/X) + 2(b+2x0)/Q atan(Q/(2x+b)) ] },
 *   X(x) = x^2 + b x + c, Q = sqrt(4c - b^2). */
void XCVWN(DFTReal rs, DFTReal A, DFTReal x0, DFTReal b, DFTReal c, DFTReal *eps, DFTReal *deps)
{
  DFTReal x = sqrt(rs), X = x*x + b*x + c, X0 = x0*x0 + b*x0 + c, Q = sqrt(4.*c - b*b);
  DFTReal at = atan(Q/(2.*x + b)), dX = 2.*x + b, d;
  *eps = A*(log(x*x/X) + 2.*b/Q*at - b*x0/X0*(log((x - x0)*(x - x0)/X) + 2.*(b + 2.*x0)/Q*at));
  // d/dx, with d atan(Q/(2x+b))/dx = -Q/(2X)
  d = A*(2./x - dX/X - b/X - b*x0/X0*(2./(x - x0) - dX/X - (b + 2.*x0)/X));
  *deps = d/(2.*x);
}

/* PW92 correlation energy per particle (paramagnetic) and its derivative:
 *   eps = -2A (1 + alpha1 rs) ln(1 + 1/(2A (beta1 rs^1/2 + beta2 rs + beta3 rs^3/2 + beta4 rs^2))) */
void XCPW92(DFTReal rs, DFTReal *eps, DFTReal *deps)
{
  const DFTReal A = 0.031091, a1 = 0.21370, b1 = 7.5957, b2 = 3.5876, b3 = 1.6382, b4 = 0.49294;
  DFTReal s = sqrt(rs), Q0 = -2.*A*(1. + a1*rs), Q1 = 2.*A*(b1*s + b2*rs + b3*rs*s + b4*rs*rs);
  DFTReal dQ1 = A*(b1/s + 2.*b2 + 3.*b3*s + 4.*b4*rs), L = log1p(1./Q1);   // log1p: 1/Q1 < eps for n < 1e-26
  *eps = Q0*L;
  *deps = -2.*A*a1*L - Q0*dQ1/(Q1*Q1 + Q1);
}

/* MacDonald-Vosko relativistic factor of the exchange energy, with beta = k_F/c:
 *   Phi = 1 - 3/2 g^2,   g = (beta eta - asinh beta)/beta^2,   eta = sqrt(1 + beta^2),
 * and dPhi/dbeta. For beta < 1/2 g and dg/dbeta come from the series of
 *   beta eta - asinh beta = int_0^beta 2t^2/sqrt(1+t^2) dt
 *                         = 2 sum_k binom(-1/2,k) beta^(2k+3)/(2k+3)
 * (40 terms: 0.25^40 < 1e-24), since the closed form loses all digits for small beta. */
void XCMacDonaldVosko(DFTReal beta, DFTReal *Phi, DFTReal *dPhi)
{
  DFTReal g, dg, eta, f, bk, b2 = beta*beta, p;
  int k;
  if (beta < 0.5)
  {
    g = 0.;
    dg = 0.;
    bk = 1.;                       // binom(-1/2, k)
    p = 1.;                        // beta^(2k)
    for (k = 0; k < 40; k++)
    {
      g += 2.*bk*p*beta/(2.*k + 3.);
      dg += 2.*bk*p*(2.*k + 1.)/(2.*k + 3.);
      bk *= -(2.*k + 1.)/(2.*k + 2.);
      p *= b2;
    }
  }
  else
  {
    eta = sqrt(1. + b2);
    f = beta*eta - asinh(beta);
    g = f/b2;
    dg = 2./eta - 2.*f/(b2*beta);
  }
  *Phi = 1. - 1.5*g*g;
  *dPhi = -3.*g*dg;
}

/* e(n, sigma) per volume, de/dn and de/dsigma (0 for the LDAs) of functional method;
 * c the speed of light (only XCRLDAMV uses it). 0 for n <= cutoff. Returns 1 for an
 * unknown method. */
int XCEvaluate(int method, double c, DFTReal n, DFTReal sigma, DFTReal cutoff, DFTReal *e, DFTReal *en, DFTReal *es)
{
  const DFTReal kappa = 0.804, betaPBE = 0.06672455060314922, mu = betaPBE*M_PI*M_PI/3.;
  const DFTReal gam = (1. - log(2.))/(M_PI*M_PI);
  DFTReal rs, kF, ex, exn, ec, dec, ecn, Phi, dPhi, bMV, s2, Fx, dFx, t2, E, A, dA, y, D, R, dRdy, dRdA, H, Hy, HA;
  *e = 0.;
  *en = 0.;
  *es = 0.;
  if (method < XCLDA || method > XCPBEVWN) return 1;
  if (!(n > cutoff)) return 0;
  rs = pow(3./(4.*M_PI*n), 1./3.);
  kF = pow(3.*M_PI*M_PI*n, 1./3.);

  // exchange: Slater, e_x = -3/4 (3/pi)^(1/3) n^(4/3)
  ex = -0.75*pow(3./M_PI, 1./3.)*pow(n, 4./3.);
  exn = -pow(3./M_PI, 1./3.)*pow(n, 1./3.);
  if (method == XCRLDAMV)
  {
    bMV = kF/c;                    // d beta/dn = beta/(3n)
    XCMacDonaldVosko(bMV, &Phi, &dPhi);
    exn = exn*Phi + ex*dPhi*bMV/(3.*n);
    ex *= Phi;
  }
  if (XCIsGGA(method))
  {
    s2 = sigma/(4.*kF*kF*n*n);     // s^2 ~ sigma n^(-8/3)
    Fx = 1. + kappa - kappa/(1. + mu*s2/kappa);
    dFx = mu/((1. + mu*s2/kappa)*(1. + mu*s2/kappa));
    *e += ex*Fx;
    *en += exn*Fx + ex*dFx*(-8./3.)*s2/n;
    *es += ex*dFx/(4.*kF*kF*n*n);
  }
  else
  {
    *e += ex;
    *en += exn;
  }

  // correlation per particle, eps_c(rs), and d eps_c/dn = -rs/(3n) d eps_c/drs
  if (method == XCLDA || method == XCPBE) XCPW92(rs, &ec, &dec);
  else XCVWN(rs, 0.0310907, -0.10498, 3.72744, 12.9352, &ec, &dec);
  ecn = -rs/(3.*n)*dec;
  *e += n*ec;
  *en += ec + n*ecn;
  if (XCIsGGA(method))
  {
    // PBE gradient correction H(eps_c, t^2), phi = 1; t^2 = sigma pi/(16 kF n^2) ~ sigma n^(-7/3)
    t2 = sigma*M_PI/(16.*kF*n*n);
    // expm1: at low density eps_c -> 0 and exp(-eps_c/gamma) - 1 cancels (and was 0 for PW92
    // below n = 1e-26, which made A infinite)
    E = expm1(-ec/gam);                                    // exp(-eps_c/gamma) - 1
    A = betaPBE/gam/E;
    dA = betaPBE/gam*((E + 1.)/gam)/(E*E);                 // dA/d eps_c
    y = t2;
    D = 1. + A*y + A*A*y*y;
    R = betaPBE/gam*y*(1. + A*y)/D;
    dRdy = betaPBE/gam*((1. + 2.*A*y)*D - y*(1. + A*y)*(A + 2.*A*A*y))/(D*D);
    dRdA = betaPBE/gam*y*y*(D - (1. + A*y)*(1. + 2.*A*y))/(D*D);
    H = gam*log(1. + R);
    Hy = gam*dRdy/(1. + R);
    HA = gam*dRdA/(1. + R);
    *e += n*H;
    *en += H + n*(HA*dA*ecn + Hy*(-7./3.)*t2/n);
    *es += n*Hy*M_PI/(16.*kF*n*n);
  }
  return 0;
}

// the method of a name ("LDA", "LDA-VWN", "RLDA-MV", "PBE", "PBE-VWN"), -1 if unknown
int XCMethodFromName(const char *name)
{
  if (strcicmp(name, "LDA") == 0) return XCLDA;
  if (strcicmp(name, "LDA-VWN") == 0) return XCLDAVWN;
  if (strcicmp(name, "RLDA-MV") == 0) return XCRLDAMV;
  if (strcicmp(name, "PBE") == 0) return XCPBE;
  if (strcicmp(name, "PBE-VWN") == 0) return XCPBEVWN;
  return -1;
}
