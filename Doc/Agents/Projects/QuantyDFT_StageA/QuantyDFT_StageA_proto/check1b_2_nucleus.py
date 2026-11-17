from dkb import *
Z=92
bs = Basis(k=9, nint=120, r1=1e-6/Z, rmax=60/Z**0.5)
H,Sm,P,Q = dirac_matrices(bs,1,-Z/bs.r,0*bs.r)
E,C = eigh(H,Sm); Ee=E[E>-c**2]
print("kappa=1 lowest pos:",Ee[:4]); print("exact 2p1/2..",[exact(Z,n,1) for n in (2,3,4)])
# intruder character: large vs small norm
i=np.where(E>-c**2)[0][0]; v=C[:,i]; Pv=v@P; Qv=v@Q
print("intruder: |P|^2=",np.sum(bs.w*Pv**2)," |Q|^2=",np.sum(bs.w*Qv**2), " r_mean=",np.sum(bs.w*bs.r*(Pv**2+Qv**2)))
# Gaussian nucleus
A=238; Rrms=(0.836*A**(1/3)+0.570)/52917.7210544
zeta=1.5/Rrms**2
print("R_rms bohr=",Rrms)
for nint in (80,120,160):
  bs = Basis(k=9, nint=nint, r1=Rrms/50, rmax=60/Z**0.5)
  V=-Z*erf(np.sqrt(zeta)*bs.r)/bs.r
  out=[]
  for kappa in (-1,1,-2):
    H,Sm,*_=dirac_matrices(bs,kappa,V,0*V); E=eigh(H,Sm,eigvals_only=True); Ee=E[E>-c**2]
    out.append(Ee[0])
  print(nint,["%.9f"%x for x in out], " point-nuc exact:", ["%.6f"%exact(Z,n,k) for n,k in ((1,-1),(2,1),(2,-2))])
