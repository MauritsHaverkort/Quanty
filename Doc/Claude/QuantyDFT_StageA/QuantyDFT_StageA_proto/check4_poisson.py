import numpy as np
from scipy.interpolate import BSpline
from scipy.integrate import quad
from dkb import knots
k=9; t=knots(k,int(__import__("sys").argv[1]),1e-5,40.0); n=len(t)-k
xg,wg=np.polynomial.legendre.leggauss(k+4); br=np.unique(t); a,b=br[:-1],br[1:]
r=(0.5*(b-a)[:,None]*xg+0.5*(a+b)[:,None]).ravel(); w=(0.5*(b-a)[:,None]*wg).ravel()
spl=BSpline(t,np.eye(n),k-1,extrapolate=False)
B=np.nan_to_num(spl(r)).T[1:]; dB=np.nan_to_num(spl.derivative()(r)).T[1:]   # keep last spline
R=t[-1]; BR=np.zeros(n-1); BR[-1]=1.0
def Yk(rho,kk):
    A=-(dB*w)@dB.T - kk*(kk+1)*(B*w/r**2)@B.T + np.outer(BR,BR)*(-kk/R)
    rhs=-(2*kk+1)*(B*w/r)@rho
    a=np.linalg.solve(A,rhs); return a@B
rho=lambda x: 4*x**2*np.exp(-2*x)*(1+0.3*x**3)      # some radial density (not normalized)
for kk in (0,1,2,4):
    Y=Yk(rho(r),kk)
    ref=np.array([x*(x**-(kk+1)*quad(lambda s: rho(s)*s**kk,0,x,epsabs=1e-14,epsrel=1e-13)[0]
                    + x**kk*quad(lambda s: rho(s)*s**-(kk+1),x,np.inf,epsabs=1e-14,epsrel=1e-13)[0]) for x in r[::37]])
    print(f"k={kk}: max|Y-Yref| on sample = {np.max(np.abs(Y[::37]-ref)):.2e}  (max|Y|={np.max(np.abs(ref)):.2f})")
# hydrogen 1s analytic check
Y0=Yk(4*r**2*np.exp(-2*r),0); print("H 1s Y0 max err:",np.max(np.abs(Y0-(1-np.exp(-2*r)*(1+r)))))
