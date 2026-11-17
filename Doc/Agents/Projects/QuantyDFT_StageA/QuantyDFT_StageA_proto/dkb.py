import numpy as np
from scipy.interpolate import BSpline
from scipy.linalg import eigh
from scipy.special import erf

c = 137.035999177

def knots(k, nint, r1, rmax):
    # k = order (degree k-1); exponential interior knots
    inner = r1 * (rmax / r1) ** (np.arange(nint) / (nint - 1))
    return np.concatenate([np.zeros(k), inner, np.full(k - 1, rmax)])

class Basis:
    def __init__(self, k, nint, r1, rmax, nq=None):
        self.k = k
        self.t = knots(k, nint, r1, rmax)
        self.n = len(self.t) - k
        nq = nq or k + 4
        xg, wg = np.polynomial.legendre.leggauss(nq)
        br = np.unique(self.t)
        a, b = br[:-1], br[1:]
        self.r = (0.5 * (b - a)[:, None] * xg[None, :] + 0.5 * (a + b)[:, None]).ravel()
        self.w = (0.5 * (b - a)[:, None] * wg[None, :]).ravel()
        spl = BSpline(self.t, np.eye(self.n), k - 1, extrapolate=False)
        self.B = np.nan_to_num(spl(self.r)).T            # n x nr
        self.dB = np.nan_to_num(spl.derivative(1)(self.r)).T
        self.d2B = np.nan_to_num(spl.derivative(2)(self.r)).T
        keep = slice(1, self.n - 1)                       # P=Q=0 at 0 and rmax
        self.B, self.dB, self.d2B = self.B[keep], self.dB[keep], self.d2B[keep]

def dkb_functions(bs, kappa):
    r = bs.r; B, dB, d2B = bs.B, bs.dB, bs.d2B
    h = 1 / (2 * c)
    # u1
    P1, dP1 = B, dB
    Q1 = h * (dB + kappa * B / r)
    dQ1 = h * (d2B + kappa * dB / r - kappa * B / r**2)
    # u2
    P2 = h * (dB - kappa * B / r)
    dP2 = h * (d2B - kappa * dB / r + kappa * B / r**2)
    Q2, dQ2 = B, dB
    return (np.vstack([P1, P2]), np.vstack([dP1, dP2]),
            np.vstack([Q1, Q2]), np.vstack([dQ1, dQ2]))

def dirac_matrices(bs, kappa, V, S):
    P, dP, Q, dQ = dkb_functions(bs, kappa)
    w, r = bs.w, bs.r
    def ip(f, g, wt):
        return (f * wt) @ g.T
    H = (ip(P, P, w * (V + S)) + ip(Q, Q, w * (V - S - 2 * c**2))
         + c * (ip(dP, Q, w) + ip(Q, dP, w))
         + c * kappa * (ip(P, Q, w / r) + ip(Q, P, w / r)))
    Sm = ip(P, P, w) + ip(Q, Q, w)
    return 0.5 * (H + H.T), 0.5 * (Sm + Sm.T), P, Q

def exact(Z, n, kappa):
    g = np.sqrt(kappa**2 - (Z / c) ** 2)
    return c**2 * ((1 + (Z / c / (n - abs(kappa) + g)) ** 2) ** -0.5 - 1)

def lval(kappa):
    return kappa if kappa > 0 else -kappa - 1

def electronic(E, n_neg_tol=-c**2):
    return E[E > n_neg_tol]

if __name__ == "__main__":
    print("=== Point nucleus hydrogen-like, DKB, B-spline order 9")
    for Z in (1, 80, 92):
        bs = Basis(k=9, nint=120, r1=1e-6 / Z, rmax=60.0 / Z**0.5 if Z > 1 else 80.0)
        V = -Z / bs.r
        for kappa in (-1, 1, -2, 2, -3):
            H, Sm, *_ = dirac_matrices(bs, kappa, V, 0 * V)
            E = eigh(H, Sm, eigvals_only=True)
            Ee = electronic(E)
            l = lval(kappa)
            ns = range(l + 1, l + 4)
            errs = []
            for i, n in enumerate(ns):
                errs.append(Ee[i] - exact(Z, n, kappa))
            nneg = np.sum(E < -c**2)
            print(f"Z={Z:3d} kappa={kappa:2d}  E1={Ee[0]:.10f}  "
                  f"err(n={ns[0]}..{ns[-1]})=" + " ".join(f"{e:+.1e}" for e in errs)
                  + f"  #neg-cont={nneg} #pos={len(Ee)}")
