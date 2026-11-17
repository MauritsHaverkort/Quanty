# Check AOC (EAL) energy functional gradient vs. Fock-operator formula
import numpy as np
from scipy.interpolate import BSpline
from scipy.linalg import eigh
from sympy import Rational as Rn
from sympy.physics.wigner import wigner_3j
from dkb import Basis, dirac_matrices, c, electronic
Z=28
bs=Basis(k=8,nint=60,r1=1e-5,rmax=30.0)
r,w=bs.r,bs.w
# Poisson basis (keep last spline)
spl=BSpline(bs.t,np.eye(bs.n),bs.k-1,extrapolate=False)
Bp=np.nan_to_num(spl(r)).T[1:]; dBp=np.nan_to_num(spl.derivative()(r)).T[1:]
BR=np.zeros(bs.n-1); BR[-1]=1; Rm=bs.t[-1]
Asolve={}
def Yk(rho,k):
    if k not in Asolve:
        Asolve[k]=np.linalg.inv(-(dBp*w)@dBp.T-k*(k+1)*(Bp*w/r**2)@Bp.T+np.outer(BR,BR)*(-k/Rm))
    return (Asolve[k]@(-(2*k+1)*(Bp*w/r)@rho))@Bp
def lof(k): return k if k>0 else -k-1
def jof(k): return abs(k)-0.5
def Gam(ka,k,kb):
    la,lb,ja,jb=lof(ka),lof(kb),jof(ka),jof(kb)
    if (la+lb+k)%2: return 0.0
    return float(wigner_3j(Rn(int(2*ja),2),k,Rn(int(2*jb),2),Rn(1,2),0,Rn(-1,2))**2)
# subshells: (kappa, shell-id); shells: occupation q, degeneracy G
subsh=[(2,'3d'),(-3,'3d'),(-1,'4s')]; shells={'3d':(8,10),'4s':(2,2)}
V=-Z/r
mats={}; vecs={}
for ka,s in subsh:
    H,Sm,P,Q=dirac_matrices(bs,ka,V,0*V); E,C=eigh(H,Sm)
    idx=np.where(E>-c**2)[0]; n0={2:0,-3:0,-1:3}[ka]
    mats[ka]=(H,Sm,P,Q); vecs[ka]=C[:,idx[n0]]+0.01*np.random.default_rng(1).standard_normal(len(E))
g=lambda ka: 2*jof(ka)+1
def qbar(a): ka,s=subsh[a]; q,G=shells[s]; return q*g(ka)/G
def Npair(a,b):
    (ka,sa),(kb,sb)=subsh[a],subsh[b]; qa,Ga=shells[sa]; qb,Gb=shells[sb]
    if a==b: return qa*(qa-1)/2*g(ka)*(g(ka)-1)/(Ga*(Ga-1))
    if sa==sb: return qa*(qa-1)*g(ka)*g(kb)/(Ga*(Ga-1))
    return qa*qb*g(ka)*g(kb)/(Ga*Gb)
def PQ(a,vv): ka=subsh[a][0]; _,_,P,Q=mats[ka]; return vv[a]@P, vv[a]@Q
def energy(vv):
    nrm=[vv[a]@mats[subsh[a][0]][1]@vv[a] for a in range(3)]
    E=0; dens={}
    for a in range(3):
        ka=subsh[a][0]; H=mats[ka][0]; E+=qbar(a)*(vv[a]@H@vv[a])/nrm[a]
        P,Q=PQ(a,vv); dens[a]=(P,Q,nrm[a])
    def rho(a,b): Pa,Qa,na=dens[a]; Pb,Qb,nb=dens[b]; return (Pa*Pb+Qa*Qb)/np.sqrt(na*nb)
    def Rk(a,b,c_,d,k): return np.sum(w*rho(a,c_)*Yk(rho(b,d),k)/r)
    for a in range(3):
        ka=subsh[a][0]
        if g(ka)>1:
            E+=Npair(a,a)*(Rk(a,a,a,a,0)-g(ka)/(g(ka)-1)*sum(Gam(ka,k,ka)*Rk(a,a,a,a,k) for k in range(2,2*lof(ka)+1,2)))
        for b in range(a+1,3):
            kb=subsh[b][0]; la,lb=lof(ka),lof(kb)
            E+=Npair(a,b)*(Rk(a,b,a,b,0)-sum(Gam(ka,k,kb)*Rk(a,b,b,a,k) for k in range(abs(la-lb),la+lb+2)))
    return E
# analytic gradient wrt coefficient vector of orbital a (orbitals normalized at evaluation point)
def norml(vv): return [vv[a]/np.sqrt(vv[a]@mats[subsh[a][0]][1]@vv[a]) for a in range(3)]
vv=norml([vecs[subsh[a][0]] for a in range(3)])
def grad_analytic(a):
    ka=subsh[a][0]; H,Sm,P,Q=mats[ka]; Pa,Qa=PQ(a,vv)
    pot=np.zeros_like(r); fx_P=np.zeros_like(r); fx_Q=np.zeros_like(r)
    if g(ka)>1:
        ra=Pa**2+Qa**2
        pot+=2*Npair(a,a)/qbar(a)*(Yk(ra,0)-g(ka)/(g(ka)-1)*sum(Gam(ka,k,ka)*Yk(ra,k) for k in range(2,2*lof(ka)+1,2)))/r
    for b in range(3):
        if b==a: continue
        kb=subsh[b][0]; Pb,Qb=PQ(b,vv); la,lb=lof(ka),lof(kb)
        pot+=Npair(a,b)/qbar(a)*Yk(Pb**2+Qb**2,0)/r
        X=sum(Gam(ka,k,kb)*Yk(Pa*Pb+Qa*Qb,k) for k in range(abs(la-lb),la+lb+2))/r
        fx_P-=Npair(a,b)/qbar(a)*X*Pb; fx_Q-=Npair(a,b)/qbar(a)*X*Qb
    fpsi=(H@vv[a]) + (P*w)@(pot*Pa+fx_P) + (Q*w)@(pot*Qa+fx_Q)
    eps=vv[a]@fpsi
    return 2*qbar(a)*(fpsi-eps*(Sm@vv[a]))   # gradient of normalized-energy functional
for a in range(3):
    ga=grad_analytic(a); gn=np.zeros_like(ga); h=1e-6
    for i in range(0,len(ga),7):
        vp=[x.copy() for x in vv]; vm=[x.copy() for x in vv]; vp[a][i]+=h; vm[a][i]-=h
        gn[i]=(energy(vp)-energy(vm))/(2*h)
    sel=slice(0,len(ga),7)
    print(f"subshell kappa={subsh[a][0]:2d}: max|grad_an - grad_fd| = {np.max(np.abs(ga[sel]-gn[sel])):.2e}  (max|grad| = {np.max(np.abs(gn[sel])):.2e})")
