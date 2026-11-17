from dkb import *
Z=26; A=56; Rrms=(0.836*A**(1/3)+0.570)/52917.7210544; zeta=1.5/Rrms**2
bs = Basis(k=9, nint=140, r1=Rrms/50, rmax=40.0)
V=-Z*erf(np.sqrt(zeta)*bs.r)/bs.r
# screened-ish: just bare ion test; look at 3d (kappa=-3) and 4s(kappa=-1) of hydrogen-like Fe
Nexp=8
for kappa,idx,name in ((-3,0,'3d5/2'),(-1,3,'4s')):
  E0=electronic(eigh(*dirac_matrices(bs,kappa,V,0*V)[:2],eigvals_only=True))[idx]
  for r0 in (0.5,1.0,2.0,4.0):
    Sc=(bs.r/r0)**Nexp
    H,Sm,P,Q=dirac_matrices(bs,kappa,V,Sc); E,C=eigh(H,Sm); ii=np.where(E>-c**2)[0][idx]
    v=C[:,ii]; Pv=v@P; Qv=v@Q
    Hv,Smv,*_=dirac_matrices(bs,kappa,V+Sc,0*V); Ev=electronic(eigh(Hv,Smv,eigvals_only=True))[idx]
    tailP=np.sum(bs.w*(bs.r>2*r0)*Pv**2); tailQ=np.sum(bs.w*(bs.r>2*r0)*Qv**2)
    sexp=np.sum(bs.w*Sc*(Pv**2-Qv**2))
    print(f"{name} r0={r0:4.1f} E_scalar-E_free={E[ii]-E0:+.6e}  E_scalar-E_vector={E[ii]-Ev:+.3e}  <beta S>={sexp:.4e} tail P,Q(r>2r0)={tailP:.1e},{tailQ:.1e}")
