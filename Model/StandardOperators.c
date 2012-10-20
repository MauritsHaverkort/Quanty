double SlaterIntegralFk(int na, int la, int nb, int lb, int nc, int lc, int nd, int ld, int k, double *Values)
{
  /* Values is an array of lenght l+1 */
  /* and contains F^k                 */
  /* at position Values[k/2]          */
  if(!((la==lb)&&(la==lc)&&(la==ld))){printf("Error in Rss not all l equal\n");exit(1);};
  if(!((na==nb)&&(na==nc)&&(na==nd))){printf("Error in Rss not one n equal\n");exit(1);};
  if(((k&1)==1)||(k<0)||(k>2*la)) return 0;
  return Values[k/2];
}

/* atomic orbitals are entered on several possible basis sets.
 * The flag OrbitalSymmetry can be set accordingly.
 *
 * Bit wise
 * 1) 
 * 2) 00 spin degenerate: 01 only down 11: only up
 * 3) 
 * 4) 00 complex spherical harmonics 01: tesseral harmonics 10: cubic harmonics 11: undifined
 *
 * Integer wise:
 * OrbitalSymmetry &  1 ==  0 => Spin degenerate basis all functions are alternated dn up dn up dn etc.
 * OrbitalSymmetry &  1 ==  1 => No spin degeneracy, next bit decides if up or dn is included
 * OrbitalSymmetry &  3 ==  3 => No spin degeneracy, only spin up Sz= 1/2 included
 * OribtalSymmetry &  3 ==  1 => No spin degeneracy, only spin dn Sz=-1/2 included
 * OrbitalSymmetry & 12 ==  0 => Basis of complex spherical harmonics ordered from ml=-l to ml=l
 * OrbitalSymmetry & 12 ==  4 => Basis of tesseral harmonics ordered ....
 * OrbitalSymmetry & 12 ==  8 => Basis of cubic harmonics ordered ...
 * OrbitalSymmetry & 12 == 12 => non defined
 */
