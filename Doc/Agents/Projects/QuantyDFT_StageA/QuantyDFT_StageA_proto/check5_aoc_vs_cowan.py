from sympy import Rational as R, S, simplify, symbols, expand
from sympy.physics.wigner import wigner_3j
half=R(1,2)
def Gam(ja,k,jb,la,lb):
    if (la+lb+k)%2: return 0
    return wigner_3j(ja,k,jb,half,0,-half)**2
def subshells(l):
    return [(l-half,l)] if l==0 else [(l-half,l),(l+half,l)] if l>0 else []
def subs(l): return [(j,l) for j in ([half] if l==0 else [l-half,l+half])]
F=symbols('F0:9'); G=symbols('G0:9')
def rel_intra(l,q):
    Gs=4*l+2; E=0; sb=subs(l)
    for ia,(ja,la) in enumerate(sb):
        ga=2*ja+1
        Naa=R(q*(q-1),2)*ga*(ga-1)/(Gs*(Gs-1))
        if ga>1: E+=Naa*(F[0]-ga/(ga-1)*sum(Gam(ja,k,ja,l,l)*F[k] for k in range(2,2*l+1,2)))
        for jb,lb in sb[ia+1:]:
            gb=2*jb+1; Nab=R(q*(q-1))*ga*gb/(Gs*(Gs-1))
            E+=Nab*(F[0]-sum(Gam(ja,k,jb,l,l)*F[k] for k in range(0,2*l+1,2)))  # G^k(ab)->F^k in nonrel limit
    return expand(E)
def cowan_intra(l,q):
    return expand(R(q*(q-1),2)*(F[0]-R(2*l+1,4*l+1)*sum(wigner_3j(l,k,l,0,0,0)**2*F[k] for k in range(2,2*l+1,2))))
def rel_inter(l1,q1,l2,q2):
    E=0
    for ja,_ in subs(l1):
        for jb,_ in subs(l2):
            ga,gb=2*ja+1,2*jb+1
            E+=R(q1*q2)*ga*gb/((4*l1+2)*(4*l2+2))*(F[0]-sum(Gam(ja,k,jb,l1,l2)*G[k] for k in range(abs(l1-l2),l1+l2+1)))
    return expand(E)
def cowan_inter(l1,q1,l2,q2):
    return expand(q1*q2*(F[0]-R(1,2)*sum(wigner_3j(l1,k,l2,0,0,0)**2*G[k] for k in range(abs(l1-l2),l1+l2+1))))
ok=True
for l in (1,2,3):
    for q in (2,3,5):
        d=simplify(rel_intra(l,q)-cowan_intra(l,q)); ok&=(d==0); print("intra l=%d q=%d diff:"%(l,q),d)
for (l1,l2) in ((0,1),(1,2),(2,3),(0,2)):
    d=simplify(rel_inter(l1,3,l2,2)-cowan_inter(l1,3,l2,2)); ok&=(d==0); print("inter",l1,l2,"diff:",d)
print("ALL OK" if ok else "MISMATCH")
# also show fixed-occupation (q_a = q g_a/G) variant differs:
l,q=2,3; Gs=10
def fixed_intra(l,q):
    Gs=4*l+2; E=0; sb=subs(l)
    for ia,(ja,la) in enumerate(sb):
        ga=2*ja+1; qa=R(q)*ga/Gs
        E+=qa*(qa-1)/2*(F[0]-ga/(ga-1)*sum(Gam(ja,k,ja,l,l)*F[k] for k in range(2,2*l+1,2)))
        for jb,lb in sb[ia+1:]:
            gb=2*jb+1; qb=R(q)*gb/Gs
            E+=qa*qb*(F[0]-sum(Gam(ja,k,jb,l,l)*F[k] for k in range(0,2*l+1,2)))
    return expand(E)
print("fixed-occupation d^3 minus Cowan:", simplify(fixed_intra(2,3)-cowan_intra(2,3)))
