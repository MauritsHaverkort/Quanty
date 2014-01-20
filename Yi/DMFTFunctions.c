#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <unistd.h>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#define omp_get_num_procs() 0
#define omp_get_max_threads() 0
#endif

#include "../Core/StructDef.c"
#include "../header.h"
#include "../Core/TimeStamp.c"
#include "../Core/Hash.c"
#include "../Core/BasicMath.c"
#include "../IO/Print.c"
#include "../IO/Read.c"
#include "../IO/BitMap.c"
#include "../Core/Error.c"
#include "../Core/WaveFunctionArithmetic.c"
#include "../Core/WaveFunctionPerturbationArithmetic.c"
#include "../Core/OperatorArithmetic.c"
#include "../Core/OperatorWaveFunction.c"
#include "../Core/OperatorWaveFunctionPerturbation.c"
#include "../Core/Lanczos.c"
#include "../Core/Spectral.c"
#include "../Core/DensityMatrix.c"
#include "../Core/Random.c"
#include "../Core/MatrixOperator.c"
#include "../Model/StandardOperators.c"
#include "../Hubbard_Model/CellFunctions.c"
#include "../Tight_Binding/Tight_Binding_lib.c"
#include "../Tight_Binding/ImpurityUnfolding.c"
#include "../DMFT/SingleSiteFunctions.c"

/*                 create Hamiltonian for a single site with two orbitals                   */

/* index:   i=0 for orbital A spin up, i=4 for bath site (connected to A) spin up,          */
/*          i=1 for orbital A spin dn, i=5 for bath site (connected to A) spin dn,          */
/*                                                                                          */
/*          i=2 for orbital B spin up, i=6 for bath site (connected to A) spin up,          */
/*          i=3 for orbital B spin dn, i=7 for bath site (connected to A) spin up,          */
/*                                                                                          */
/*      for i>=8, first spin up (ordered as A-val, A-con, B-val, B-con), then spin down.    */

/*   additional parameters: t1 is defined as inter orbital mixing, J is the Hund coupling.  */

uint32_t TwoOrbSecondDoubleTriDiagonalMatrixToAndersonImpurityOperator(OperatorType *O, SecondDoubleTriDiagonalMatrixType MAup, SecondDoubleTriDiagonalMatrixType MAdn, SecondDoubleTriDiagonalMatrixType MBup, SecondDoubleTriDiagonalMatrixType MBdn, double U, double t1, double J, double mu)
{
    uint32_t i;
    uint16_t oadd[4];
    InitNamedOperator(O,"Two-Orbital Anderson Impurity Hamiltonian");
    O->NFermion=MAup.Dt.val.N+MAup.Dt.con.N+MAdn.Dt.val.N+MAdn.Dt.con.N+MBup.Dt.val.N+MBup.Dt.con.N+MBdn.Dt.val.N+MBdn.Dt.con.N+8;
    O->NWord=Min(O->NFermion,1+(O->NFermion-1)/8);
    InitOperatorFixedLength(O, 2, 8);
    InitOperatorFixedLength(O, 4, 8);
    
    RealOperatorAddLaderNormalOrder(O,oadd,0,-50,&QDetNotOrderedFcFaBcBa);
    
    /* U */
    oadd[3]=CrAnMask+0; oadd[2]=CrAnMask+1; oadd[1]=1; oadd[0]=0; RealOperatorAddLaderNormalOrder(O,oadd,4,U,&QDetNotOrderedFcFaBcBa);
    oadd[3]=CrAnMask+2; oadd[2]=CrAnMask+3; oadd[1]=3; oadd[0]=2; RealOperatorAddLaderNormalOrder(O,oadd,4,U,&QDetNotOrderedFcFaBcBa);
    /* U' */
    oadd[3]=CrAnMask+0; oadd[2]=CrAnMask+3; oadd[1]=3; oadd[0]=0; RealOperatorAddLaderNormalOrder(O,oadd,4,U-2*J,&QDetNotOrderedFcFaBcBa);
    oadd[3]=CrAnMask+1; oadd[2]=CrAnMask+2; oadd[1]=2; oadd[0]=1; RealOperatorAddLaderNormalOrder(O,oadd,4,U-2*J,&QDetNotOrderedFcFaBcBa);
    /* U'' */
    oadd[3]=CrAnMask+0; oadd[2]=CrAnMask+2; oadd[1]=2; oadd[0]=0; RealOperatorAddLaderNormalOrder(O,oadd,4,U-3*J,&QDetNotOrderedFcFaBcBa);
    oadd[3]=CrAnMask+1; oadd[2]=CrAnMask+3; oadd[1]=3; oadd[0]=1; RealOperatorAddLaderNormalOrder(O,oadd,4,U-3*J,&QDetNotOrderedFcFaBcBa);
    
    /* onsite mu */
    for (i=0; i<4; i++) {
        oadd[1]=CrAnMask+i; oadd[0]=i; RealOperatorAddLaderNormalOrder(O,oadd,2,mu,&QDetNotOrderedFcFaBcBa);
    }
    
    /* first A bath site spin up */
    oadd[1]=CrAnMask+0; oadd[0]=0; RealOperatorAddLaderNormalOrder(O,oadd,2,MAup.Dt.a0,&QDetNotOrderedFcFaBcBa);
    oadd[1]=CrAnMask+0; oadd[0]=4; RealOperatorAddLaderNormalOrder(O,oadd,2,MAup.Dt.b0,&QDetNotOrderedFcFaBcBa);
    oadd[1]=CrAnMask+4; oadd[0]=0; RealOperatorAddLaderNormalOrder(O,oadd,2,MAup.Dt.b0,&QDetNotOrderedFcFaBcBa);
    /* first A bath site spin dn */
    oadd[1]=CrAnMask+1; oadd[0]=1; RealOperatorAddLaderNormalOrder(O,oadd,2,MAdn.Dt.a0,&QDetNotOrderedFcFaBcBa);
    oadd[1]=CrAnMask+1; oadd[0]=5; RealOperatorAddLaderNormalOrder(O,oadd,2,MAdn.Dt.b0,&QDetNotOrderedFcFaBcBa);
    oadd[1]=CrAnMask+5; oadd[0]=1; RealOperatorAddLaderNormalOrder(O,oadd,2,MAdn.Dt.b0,&QDetNotOrderedFcFaBcBa);
    /* first B bath site spin up */
    oadd[1]=CrAnMask+2; oadd[0]=2; RealOperatorAddLaderNormalOrder(O,oadd,2,MBup.Dt.a0,&QDetNotOrderedFcFaBcBa);
    oadd[1]=CrAnMask+2; oadd[0]=6; RealOperatorAddLaderNormalOrder(O,oadd,2,MBup.Dt.b0,&QDetNotOrderedFcFaBcBa);
    oadd[1]=CrAnMask+6; oadd[0]=2; RealOperatorAddLaderNormalOrder(O,oadd,2,MBup.Dt.b0,&QDetNotOrderedFcFaBcBa);
    /* first B bath site spin dn */
    oadd[1]=CrAnMask+3; oadd[0]=3; RealOperatorAddLaderNormalOrder(O,oadd,2,MBdn.Dt.a0,&QDetNotOrderedFcFaBcBa);
    oadd[1]=CrAnMask+3; oadd[0]=7; RealOperatorAddLaderNormalOrder(O,oadd,2,MBdn.Dt.b0,&QDetNotOrderedFcFaBcBa);
    oadd[1]=CrAnMask+7; oadd[0]=3; RealOperatorAddLaderNormalOrder(O,oadd,2,MBdn.Dt.b0,&QDetNotOrderedFcFaBcBa);
    
    /* A spin up valence */
    oadd[1]=CrAnMask+4; oadd[0]=4; RealOperatorAddLaderNormalOrder(O,oadd,2,MAup.Dt.val.a[0]/2.0,&QDetNotOrderedFcFaBcBa);
    oadd[1]=CrAnMask+4; oadd[0]=8; RealOperatorAddLaderNormalOrder(O,oadd,2,MAup.Dt.val.b[0],&QDetNotOrderedFcFaBcBa);
    oadd[1]=CrAnMask+8; oadd[0]=4; RealOperatorAddLaderNormalOrder(O,oadd,2,MAup.Dt.val.b[0],&QDetNotOrderedFcFaBcBa);
    oadd[1]=CrAnMask+8; oadd[0]=8; RealOperatorAddLaderNormalOrder(O,oadd,2,MAup.Dt.val.a[1],&QDetNotOrderedFcFaBcBa);
    for(i=2;i<=MAup.Dt.val.N;i++)
    {
        oadd[1]=CrAnMask+i+6; oadd[0]=i+7; RealOperatorAddLaderNormalOrder(O,oadd,2,MAup.Dt.val.b[i-1],&QDetNotOrderedFcFaBcBa);
        oadd[1]=CrAnMask+i+7; oadd[0]=i+6; RealOperatorAddLaderNormalOrder(O,oadd,2,MAup.Dt.val.b[i-1],&QDetNotOrderedFcFaBcBa);
        oadd[1]=CrAnMask+i+7; oadd[0]=i+7; RealOperatorAddLaderNormalOrder(O,oadd,2,MAup.Dt.val.a[i]  ,&QDetNotOrderedFcFaBcBa);
    }
    
    /* A spin up conductance */
    oadd[1]=CrAnMask+4;
    oadd[0]=4;
    RealOperatorAddLaderNormalOrder(O,oadd,2,MAup.Dt.con.a[0]/2.0,&QDetNotOrderedFcFaBcBa);
    
    oadd[1]=CrAnMask+4;
    oadd[0]=MAup.Dt.val.N+8;
    RealOperatorAddLaderNormalOrder(O,oadd,2,MAup.Dt.con.b[0],&QDetNotOrderedFcFaBcBa);
    
    oadd[1]=CrAnMask+MAup.Dt.val.N+8;
    oadd[0]=4;
    RealOperatorAddLaderNormalOrder(O,oadd,2,MAup.Dt.con.b[0],&QDetNotOrderedFcFaBcBa);
    
    oadd[1]=CrAnMask+MAup.Dt.val.N+8;
    oadd[0]=MAup.Dt.val.N+8;
    RealOperatorAddLaderNormalOrder(O,oadd,2,MAup.Dt.con.a[1],&QDetNotOrderedFcFaBcBa);
    
    for(i=2;i<=MAup.Dt.con.N;i++)
    {
        oadd[1]=CrAnMask+i+6+MAup.Dt.val.N;
        oadd[0]=i+7+MAup.Dt.val.N;
        RealOperatorAddLaderNormalOrder(O,oadd,2,MAup.Dt.con.b[i-1],&QDetNotOrderedFcFaBcBa);
        
        oadd[1]=CrAnMask+i+7+MAup.Dt.val.N;
        oadd[0]=i+6+MAup.Dt.val.N;
        RealOperatorAddLaderNormalOrder(O,oadd,2,MAup.Dt.con.b[i-1],&QDetNotOrderedFcFaBcBa);
        
        oadd[1]=CrAnMask+i+7+MAup.Dt.val.N;
        oadd[0]=i+7+MAup.Dt.val.N;
        RealOperatorAddLaderNormalOrder(O,oadd,2,MAup.Dt.con.a[i],&QDetNotOrderedFcFaBcBa);
    }
    
    /* B spin up valence */
    oadd[1]=CrAnMask+6;
    oadd[0]=6;
    RealOperatorAddLaderNormalOrder(O,oadd,2,MBup.Dt.val.a[0]/2.0,&QDetNotOrderedFcFaBcBa);
    
    oadd[1]=CrAnMask+6;
    oadd[0]=MAup.Dt.val.N+MAup.Dt.con.N+8;
    RealOperatorAddLaderNormalOrder(O,oadd,2,MBup.Dt.val.b[0],&QDetNotOrderedFcFaBcBa);
    
    oadd[1]=CrAnMask+MAup.Dt.val.N+MAup.Dt.con.N+8;
    oadd[0]=6;
    RealOperatorAddLaderNormalOrder(O,oadd,2,MBup.Dt.val.b[0],&QDetNotOrderedFcFaBcBa);
    
    oadd[1]=CrAnMask+MAup.Dt.val.N+MAup.Dt.con.N+8;
    oadd[0]=MAup.Dt.val.N+MAup.Dt.con.N+8;
    RealOperatorAddLaderNormalOrder(O,oadd,2,MBup.Dt.val.a[1],&QDetNotOrderedFcFaBcBa);
    
    for(i=2;i<=MBup.Dt.val.N;i++)
    {
        oadd[1]=CrAnMask+i+6+MAup.Dt.val.N+MAup.Dt.con.N;
        oadd[0]=i+7+MAup.Dt.val.N+MAup.Dt.con.N;
        RealOperatorAddLaderNormalOrder(O,oadd,2,MBup.Dt.val.b[i-1],&QDetNotOrderedFcFaBcBa);
        
        oadd[1]=CrAnMask+i+7+MAup.Dt.val.N+MAup.Dt.con.N;
        oadd[0]=i+6+MAup.Dt.val.N+MAup.Dt.con.N;
        RealOperatorAddLaderNormalOrder(O,oadd,2,MBup.Dt.val.b[i-1],&QDetNotOrderedFcFaBcBa);
        
        oadd[1]=CrAnMask+i+7+MAup.Dt.val.N+MAup.Dt.con.N;
        oadd[0]=i+7+MAup.Dt.val.N+MAup.Dt.con.N;
        RealOperatorAddLaderNormalOrder(O,oadd,2,MBup.Dt.val.a[i],&QDetNotOrderedFcFaBcBa);
    }
    
    /* B spin up conductance */
    oadd[1]=CrAnMask+6;
    oadd[0]=6;
    RealOperatorAddLaderNormalOrder(O,oadd,2,MBup.Dt.con.a[0]/2.0,&QDetNotOrderedFcFaBcBa);
    
    oadd[1]=CrAnMask+6;
    oadd[0]=MAup.Dt.val.N+MAup.Dt.con.N+MBup.Dt.val.N+8;
    RealOperatorAddLaderNormalOrder(O,oadd,2,MBup.Dt.con.b[0],&QDetNotOrderedFcFaBcBa);
    
    oadd[1]=CrAnMask+MAup.Dt.val.N+MAup.Dt.con.N+MBup.Dt.val.N+8;
    oadd[0]=6;
    RealOperatorAddLaderNormalOrder(O,oadd,2,MBup.Dt.con.b[0],&QDetNotOrderedFcFaBcBa);
    
    oadd[1]=CrAnMask+MAup.Dt.val.N+MAup.Dt.con.N+MBup.Dt.val.N+8;
    oadd[0]=MAup.Dt.val.N+MAup.Dt.con.N+MBup.Dt.val.N+8;
    RealOperatorAddLaderNormalOrder(O,oadd,2,MBup.Dt.con.a[1],&QDetNotOrderedFcFaBcBa);
    
    for(i=2;i<=MBup.Dt.con.N;i++)
    {
        oadd[1]=CrAnMask+i+6+MAup.Dt.val.N+MAup.Dt.con.N+MBup.Dt.val.N;
        oadd[0]=i+7+MAup.Dt.val.N+MAup.Dt.con.N+MBup.Dt.val.N;
        RealOperatorAddLaderNormalOrder(O,oadd,2,MBup.Dt.con.b[i-1],&QDetNotOrderedFcFaBcBa);
        
        oadd[1]=CrAnMask+i+7+MAup.Dt.val.N+MAup.Dt.con.N+MBup.Dt.val.N;
        oadd[0]=i+6+MAup.Dt.val.N+MAup.Dt.con.N+MBup.Dt.val.N;
        RealOperatorAddLaderNormalOrder(O,oadd,2,MBup.Dt.con.b[i-1],&QDetNotOrderedFcFaBcBa);
        
        oadd[1]=CrAnMask+i+7+MAup.Dt.val.N+MAup.Dt.con.N+MBup.Dt.val.N;
        oadd[0]=i+7+MAup.Dt.val.N+MAup.Dt.con.N+MBup.Dt.val.N;
        RealOperatorAddLaderNormalOrder(O,oadd,2,MBup.Dt.con.a[i],&QDetNotOrderedFcFaBcBa);
    }
    
    /* A spin dn valence */
    oadd[1]=CrAnMask+5;
    oadd[0]=5;
    RealOperatorAddLaderNormalOrder(O,oadd,2,MAdn.Dt.val.a[0]/2.0,&QDetNotOrderedFcFaBcBa);
    
    oadd[1]=CrAnMask+5;
    oadd[0]=MAup.Dt.val.N+MAup.Dt.con.N+MBup.Dt.val.N+MBup.Dt.con.N+8;
    RealOperatorAddLaderNormalOrder(O,oadd,2,MAdn.Dt.val.b[0],&QDetNotOrderedFcFaBcBa);
    
    oadd[1]=CrAnMask+MAup.Dt.val.N+MAup.Dt.con.N+MBup.Dt.val.N+MBup.Dt.con.N+8;
    oadd[0]=5;
    RealOperatorAddLaderNormalOrder(O,oadd,2,MAdn.Dt.val.b[0],&QDetNotOrderedFcFaBcBa);
    
    oadd[1]=CrAnMask+MAup.Dt.val.N+MAup.Dt.con.N+MBup.Dt.val.N+MBup.Dt.con.N+8;
    oadd[0]=MAup.Dt.val.N+MAup.Dt.con.N+MBup.Dt.val.N+MBup.Dt.con.N+8;
    RealOperatorAddLaderNormalOrder(O,oadd,2,MAdn.Dt.val.a[1],&QDetNotOrderedFcFaBcBa);
    
    for(i=2;i<=MAdn.Dt.val.N;i++)
    {
        oadd[1]=CrAnMask+i+6+MAup.Dt.val.N+MAup.Dt.con.N+MBup.Dt.val.N+MBup.Dt.con.N;
        oadd[0]=i+7+MAup.Dt.val.N+MAup.Dt.con.N+MBup.Dt.val.N+MBup.Dt.con.N;
        RealOperatorAddLaderNormalOrder(O,oadd,2,MAdn.Dt.val.b[i-1],&QDetNotOrderedFcFaBcBa);
        
        oadd[1]=CrAnMask+i+7+MAup.Dt.val.N+MAup.Dt.con.N+MBup.Dt.val.N+MBup.Dt.con.N;
        oadd[0]=i+6+MAup.Dt.val.N+MAup.Dt.con.N+MBup.Dt.val.N+MBup.Dt.con.N;
        RealOperatorAddLaderNormalOrder(O,oadd,2,MAdn.Dt.val.b[i-1],&QDetNotOrderedFcFaBcBa);
        
        oadd[1]=CrAnMask+i+7+MAup.Dt.val.N+MAup.Dt.con.N+MBup.Dt.val.N+MBup.Dt.con.N;
        oadd[0]=i+7+MAup.Dt.val.N+MAup.Dt.con.N+MBup.Dt.val.N+MBup.Dt.con.N;
        RealOperatorAddLaderNormalOrder(O,oadd,2,MAdn.Dt.val.a[i],&QDetNotOrderedFcFaBcBa);
    }
    
    /* A spin dn conductance */
    oadd[1]=CrAnMask+5;
    oadd[0]=5;
    RealOperatorAddLaderNormalOrder(O,oadd,2,MAdn.Dt.con.a[0]/2.0,&QDetNotOrderedFcFaBcBa);
    
    oadd[1]=CrAnMask+5;
    oadd[0]=MAup.Dt.val.N+MAup.Dt.con.N+MBup.Dt.val.N+MBup.Dt.con.N+MAdn.Dt.val.N+8;
    RealOperatorAddLaderNormalOrder(O,oadd,2,MAdn.Dt.con.b[0],&QDetNotOrderedFcFaBcBa);
    
    oadd[1]=CrAnMask+MAup.Dt.val.N+MAup.Dt.con.N+MBup.Dt.val.N+MBup.Dt.con.N+MAdn.Dt.val.N+8;
    oadd[0]=5;
    RealOperatorAddLaderNormalOrder(O,oadd,2,MAdn.Dt.con.b[0],&QDetNotOrderedFcFaBcBa);
    
    oadd[1]=CrAnMask+MAup.Dt.val.N+MAup.Dt.con.N+MBup.Dt.val.N+MBup.Dt.con.N+MAdn.Dt.val.N+8;
    oadd[0]=MAup.Dt.val.N+MAup.Dt.con.N+MBup.Dt.val.N+MBup.Dt.con.N+MAdn.Dt.val.N+8;
    RealOperatorAddLaderNormalOrder(O,oadd,2,MAdn.Dt.con.a[1],&QDetNotOrderedFcFaBcBa);
    
    for(i=2;i<=MAdn.Dt.con.N;i++)
    {
        oadd[1]=CrAnMask+i+6+MAup.Dt.val.N+MAup.Dt.con.N+MBup.Dt.val.N+MBup.Dt.con.N+MAdn.Dt.val.N;
        oadd[0]=i+7+MAup.Dt.val.N+MAup.Dt.con.N+MBup.Dt.val.N+MBup.Dt.con.N+MAdn.Dt.val.N;
        RealOperatorAddLaderNormalOrder(O,oadd,2,MAdn.Dt.con.b[i-1],&QDetNotOrderedFcFaBcBa);
        
        oadd[1]=CrAnMask+i+7+MAup.Dt.val.N+MAup.Dt.con.N+MBup.Dt.val.N+MBup.Dt.con.N+MAdn.Dt.val.N;
        oadd[0]=i+6+MAup.Dt.val.N+MAup.Dt.con.N+MBup.Dt.val.N+MBup.Dt.con.N+MAdn.Dt.val.N;
        RealOperatorAddLaderNormalOrder(O,oadd,2,MAdn.Dt.con.b[i-1],&QDetNotOrderedFcFaBcBa);
        
        oadd[1]=CrAnMask+i+7+MAup.Dt.val.N+MAup.Dt.con.N+MBup.Dt.val.N+MBup.Dt.con.N+MAdn.Dt.val.N;
        oadd[0]=i+7+MAup.Dt.val.N+MAup.Dt.con.N+MBup.Dt.val.N+MBup.Dt.con.N+MAdn.Dt.val.N;
        RealOperatorAddLaderNormalOrder(O,oadd,2,MAdn.Dt.con.a[i],&QDetNotOrderedFcFaBcBa);
    }
    
    /* B spin dn valence */
    oadd[1]=CrAnMask+7;
    oadd[0]=7;
    RealOperatorAddLaderNormalOrder(O,oadd,2,MBdn.Dt.val.a[0]/2.0,&QDetNotOrderedFcFaBcBa);
    
    oadd[1]=CrAnMask+7;
    oadd[0]=MAup.Dt.val.N+MAup.Dt.con.N+MBup.Dt.val.N+MBup.Dt.con.N+MAdn.Dt.val.N+MAdn.Dt.con.N+8;
    RealOperatorAddLaderNormalOrder(O,oadd,2,MBdn.Dt.val.b[0],&QDetNotOrderedFcFaBcBa);
    
    oadd[1]=CrAnMask+MAup.Dt.val.N+MAup.Dt.con.N+MBup.Dt.val.N+MBup.Dt.con.N+MAdn.Dt.val.N+MAdn.Dt.con.N+8;
    oadd[0]=7;
    RealOperatorAddLaderNormalOrder(O,oadd,2,MBdn.Dt.val.b[0],&QDetNotOrderedFcFaBcBa);
    
    oadd[1]=CrAnMask+MAup.Dt.val.N+MAup.Dt.con.N+MBup.Dt.val.N+MBup.Dt.con.N+MAdn.Dt.val.N+MAdn.Dt.con.N+8;
    oadd[0]=MAup.Dt.val.N+MAup.Dt.con.N+MBup.Dt.val.N+MBup.Dt.con.N+MAdn.Dt.val.N+MAdn.Dt.con.N+8;
    RealOperatorAddLaderNormalOrder(O,oadd,2,MBdn.Dt.val.a[1],&QDetNotOrderedFcFaBcBa);
    
    for(i=2;i<=MBdn.Dt.val.N;i++)
    {
        oadd[1]=CrAnMask+i+6+MAup.Dt.val.N+MAup.Dt.con.N+MBup.Dt.val.N+MBup.Dt.con.N+MAdn.Dt.val.N+MAdn.Dt.con.N;
        oadd[0]=i+7+MAup.Dt.val.N+MAup.Dt.con.N+MBup.Dt.val.N+MBup.Dt.con.N+MAdn.Dt.val.N+MAdn.Dt.con.N;
        RealOperatorAddLaderNormalOrder(O,oadd,2,MBdn.Dt.val.b[i-1],&QDetNotOrderedFcFaBcBa);
        
        oadd[1]=CrAnMask+i+7+MAup.Dt.val.N+MAup.Dt.con.N+MBup.Dt.val.N+MBup.Dt.con.N+MAdn.Dt.val.N+MAdn.Dt.con.N;
        oadd[0]=i+6+MAup.Dt.val.N+MAup.Dt.con.N+MBup.Dt.val.N+MBup.Dt.con.N+MAdn.Dt.val.N+MAdn.Dt.con.N;
        RealOperatorAddLaderNormalOrder(O,oadd,2,MBdn.Dt.val.b[i-1],&QDetNotOrderedFcFaBcBa);
        
        oadd[1]=CrAnMask+i+7+MAup.Dt.val.N+MAup.Dt.con.N+MBup.Dt.val.N+MBup.Dt.con.N+MAdn.Dt.val.N+MAdn.Dt.con.N;
        oadd[0]=i+7+MAup.Dt.val.N+MAup.Dt.con.N+MBup.Dt.val.N+MBup.Dt.con.N+MAdn.Dt.val.N+MAdn.Dt.con.N;
        RealOperatorAddLaderNormalOrder(O,oadd,2,MBdn.Dt.val.a[i],&QDetNotOrderedFcFaBcBa);
    }
    
    /* B spin dn conductance */
    oadd[1]=CrAnMask+7;
    oadd[0]=7;
    RealOperatorAddLaderNormalOrder(O,oadd,2,MBdn.Dt.con.a[0]/2.0,&QDetNotOrderedFcFaBcBa);
    
    oadd[1]=CrAnMask+7;
    oadd[0]=MAup.Dt.val.N+MAup.Dt.con.N+MBup.Dt.val.N+MBup.Dt.con.N+MAdn.Dt.val.N+MAdn.Dt.con.N+MBdn.Dt.val.N+8;
    RealOperatorAddLaderNormalOrder(O,oadd,2,MBdn.Dt.con.b[0],&QDetNotOrderedFcFaBcBa);
    
    oadd[1]=CrAnMask+MAup.Dt.val.N+MAup.Dt.con.N+MBup.Dt.val.N+MBup.Dt.con.N+MAdn.Dt.val.N+MAdn.Dt.con.N+MBdn.Dt.val.N+8;
    oadd[0]=7;
    RealOperatorAddLaderNormalOrder(O,oadd,2,MBdn.Dt.con.b[0],&QDetNotOrderedFcFaBcBa);
    
    oadd[1]=CrAnMask+MAup.Dt.val.N+MAup.Dt.con.N+MBup.Dt.val.N+MBup.Dt.con.N+MAdn.Dt.val.N+MAdn.Dt.con.N+MBdn.Dt.val.N+8;
    oadd[0]=MAup.Dt.val.N+MAup.Dt.con.N+MBup.Dt.val.N+MBup.Dt.con.N+MAdn.Dt.val.N+MAdn.Dt.con.N+MBdn.Dt.val.N+8;
    RealOperatorAddLaderNormalOrder(O,oadd,2,MBdn.Dt.con.a[1],&QDetNotOrderedFcFaBcBa);
    
    for(i=2;i<=MBdn.Dt.con.N;i++)
    {
        oadd[1]=CrAnMask+i+6+MAup.Dt.val.N+MAup.Dt.con.N+MBup.Dt.val.N+MBup.Dt.con.N+MAdn.Dt.val.N+MAdn.Dt.con.N+MBdn.Dt.val.N;
        oadd[0]=i+7+MAup.Dt.val.N+MAup.Dt.con.N+MBup.Dt.val.N+MBup.Dt.con.N+MAdn.Dt.val.N+MAdn.Dt.con.N+MBdn.Dt.val.N;
        RealOperatorAddLaderNormalOrder(O,oadd,2,MBdn.Dt.con.b[i-1],&QDetNotOrderedFcFaBcBa);
        
        oadd[1]=CrAnMask+i+7+MAup.Dt.val.N+MAup.Dt.con.N+MBup.Dt.val.N+MBup.Dt.con.N+MAdn.Dt.val.N+MAdn.Dt.con.N+MBdn.Dt.val.N;
        oadd[0]=i+6+MAup.Dt.val.N+MAup.Dt.con.N+MBup.Dt.val.N+MBup.Dt.con.N+MAdn.Dt.val.N+MAdn.Dt.con.N+MBdn.Dt.val.N;
        RealOperatorAddLaderNormalOrder(O,oadd,2,MBdn.Dt.con.b[i-1],&QDetNotOrderedFcFaBcBa);
        
        oadd[1]=CrAnMask+i+7+MAup.Dt.val.N+MAup.Dt.con.N+MBup.Dt.val.N+MBup.Dt.con.N+MAdn.Dt.val.N+MAdn.Dt.con.N+MBdn.Dt.val.N;
        oadd[0]=i+7+MAup.Dt.val.N+MAup.Dt.con.N+MBup.Dt.val.N+MBup.Dt.con.N+MAdn.Dt.val.N+MAdn.Dt.con.N+MBdn.Dt.val.N;
        RealOperatorAddLaderNormalOrder(O,oadd,2,MBdn.Dt.con.a[i],&QDetNotOrderedFcFaBcBa);
    }
    
    return 0;
}

/* create initial wave-functions */

uint32_t CreateWaveFunctionTwoOrbDMFTSecond(WaveFunctionType *psi, uint32_t nBathVA, uint32_t nBathCA, uint32_t nBathVB, uint32_t nBathCB, uint32_t NBitsKey)
{
    uint32_t i;
    uint8_t *Determinant;
    psi->NBitsKey=NBitsKey;
    psi->QComplex=0;
    psi->NBoson=0;
    psi->NFermion=8+2*nBathCA+2*nBathVA+2*nBathCB+2*nBathVB;
    if(WaveFunctionInit(psi)){printf("WaveFunctionInit failed in CreateWaveFunctionDMFT\n");fflush(stdout);return 1;};
    sprintf(psi->Name,"DMFT %i %i %i %i",nBathVA,nBathCA,nBathVB,nBathCB);
    Determinant=calloc(psi->NWord,sizeof(uint8_t));
    OperatorCr(Determinant,0,psi->NFermion);
    OperatorCr(Determinant,1,psi->NFermion);
    OperatorCr(Determinant,2,psi->NFermion);
    OperatorCr(Determinant,3,psi->NFermion);
    /* fill in A BathV with spin up */
    for(i=8;i<8+nBathVA;i++)
        OperatorCr(Determinant,i,psi->NFermion);
    /* fill in B BathV with spin up */
    for(i=8+nBathVA+nBathCA;i<8+nBathVA+nBathCA+nBathVB;i++)
        OperatorCr(Determinant,i,psi->NFermion);
    /* fill in A BathV with spin dn */
    for(i=8+nBathVA+nBathCA+nBathVB+nBathCB;i<8+2*nBathVA+nBathCA+nBathVB+nBathCB;i++)
        OperatorCr(Determinant,i,psi->NFermion);
    /* fill in B BathV with spin dn */
    for(i=8+2*nBathVA+2*nBathCA+nBathVB+nBathCB;i<8+2*nBathVA+2*nBathCA+2*nBathVB+nBathCB;i++)
        OperatorCr(Determinant,i,psi->NFermion);
    RealWaveFunctionAddElement(psi, Determinant, 1.0);
    free(Determinant);
    return 0;
}

/* define the structure of DMFTSettings for 2 orbitals */

typedef struct TwoOrbDMFTSettings
{
    uint32_t Opt;                              /* option for input, 1 for input as a list of poles */
    char     G0FileNameA[256];                  /* filename of G0 in tridiagonal representation */
    char     G0FileNameB[256];
    char     GRootFileNameA[256];               /* filename of GBath in tridiagonal representation */
    char     GRootFileNameB[256];
    char     SpecRootFileName[256];            /* filename of spectral functions in spectral representation */
    char     WaveFunctionRootFileName[256];    /* filename of wave-function  */
    uint32_t iStart;                           /* start value of loop number */
    uint32_t Ni;                               /* number of loops to do */
    uint32_t Dyson;                            /* 0 no Dyson like perturbation - 1 full perturbation - 2 ... */
    uint32_t NGBath;                           /* Number of poles in the resulting bath green function */
    uint32_t NGCluster;                        /* Number of poles in the resulting cluster green function */
    double   Mix;                              /* mixing of old into new function */
    double   U;                                /* valence - valence Coulomb repulsion */
    double   t1;                               /* hopping between two orbitals */
    double   J;                                /* Hund's coupling */
    double   Nelectron;                        /* number of electrons per spin per site */
    uint32_t NBathA;                           /* Number of Bath states to include */
    uint32_t NBathB;
    double   mu1;                              /* Starting values for chemical potential search */
    double   mu2;
    uint32_t Spectra;                          /* 0000 no spectra will be calculated at the end of the DMFT loop */
                                               /* 0001 Magnetic susceptibility will be calculated (nup-ndn) */
                                               /* 0010 occupation will be calculated (nup+ndn) */
                                               /* 0100 XAS will be calculated with NQ values for core - valence interaction */
    uint32_t NQ;                               /* number of XAS spectra */
    double   *Q;                               /* NQ values of core - valence interaction strengths */
    uint32_t NE;                               /* number of points in the spectrum */
    double   Emin;                             /* minimum energy */
    double   Emax;                             /* maximum energy */
    double   Gamma;                            /* Broadening only important for the Spec file */
    double   zeroLancEnergy;                   /* minimum energy difference for convergence of total energy to be different from zero */
    double   zeroLancWavefunction;             /* minnimum pre-factor in wave-function to be different from zero */
    uint32_t NTriDiag;                         /* number of tri-diagonalization steps in Lanczos */
    uint32_t NBitsKey;                         /* 2^NBitsKey == lenght of hash key */
}
TwoOrbDMFTSettings;

uint32_t StreamReadTwoOrbDMFTSettings(TwoOrbDMFTSettings * O, FILE *stream)
{
    char comments[1000];
    uint32_t i;
    fscanf(stream,"%i",&(O->Opt));                      fgets(comments,1000,stream);
    fscanf(stream,"%s",(O->G0FileNameA));               fgets(comments,1000,stream);
    fscanf(stream,"%s",(O->G0FileNameB));               fgets(comments,1000,stream);
    fscanf(stream,"%s",(O->GRootFileNameA));            fgets(comments,1000,stream);
    fscanf(stream,"%s",(O->GRootFileNameB));            fgets(comments,1000,stream);
    fscanf(stream,"%s",(O->SpecRootFileName));          fgets(comments,1000,stream);
    fscanf(stream,"%s",(O->WaveFunctionRootFileName));  fgets(comments,1000,stream);
    fscanf(stream,"%i",&(O->iStart));                   fgets(comments,1000,stream);
    fscanf(stream,"%i",&(O->Ni));                       fgets(comments,1000,stream);
    fscanf(stream,"%i",&(O->Dyson));                    fgets(comments,1000,stream);
    fscanf(stream,"%i",&(O->NGCluster));                fgets(comments,1000,stream);
    fscanf(stream,"%i",&(O->NGBath));                   fgets(comments,1000,stream);
    fscanf(stream,"%lf",&(O->Mix));                     fgets(comments,1000,stream);
    fscanf(stream,"%lf",&(O->U));                       fgets(comments,1000,stream);
    fscanf(stream,"%lf",&(O->t1));                      fgets(comments,1000,stream);
    fscanf(stream,"%lf",&(O->J));                       fgets(comments,1000,stream);
    fscanf(stream,"%i",&(O->NBathA));                   fgets(comments,1000,stream);
    fscanf(stream,"%i",&(O->NBathB));                   fgets(comments,1000,stream);
    fscanf(stream,"%lf",&(O->mu1));                     fgets(comments,1000,stream);
    fscanf(stream,"%lf",&(O->mu2));                     fgets(comments,1000,stream);
    fscanf(stream,"%i",&(O->Spectra));
    O->NQ=1;
    if((O->Spectra)&12)
    {
        fscanf(stream,"%i",&(O->NQ));
        O->Q=malloc(O->NQ*sizeof(double));
        for(i=0;i<O->NQ;i++)
        {
            fscanf(stream,"%lf",&(O->Q[i]));
        }
        if(O->NQ==0)
        {
            O->Spectra&=!((uint32_t)12);
        }
    }                                                  fgets(comments,1000,stream);
    fscanf(stream,"%i",&(O->NE));                      fgets(comments,1000,stream);
    fscanf(stream,"%lf",&(O->Emin));                   fgets(comments,1000,stream);
    fscanf(stream,"%lf",&(O->Emax));                   fgets(comments,1000,stream);
    fscanf(stream,"%lf",&(O->Gamma));                  fgets(comments,1000,stream);
    fscanf(stream,"%lf",&(O->zeroLancEnergy));         fgets(comments,1000,stream);
    fscanf(stream,"%lf",&(O->zeroLancWavefunction));   fgets(comments,1000,stream);
    fscanf(stream,"%i",&(O->NTriDiag));                fgets(comments,1000,stream);
    fscanf(stream,"%i",&(O->NBitsKey));                fgets(comments,1000,stream);
    return 0;
}

uint32_t FileReadTwoOrbDMFTSettings(TwoOrbDMFTSettings * O, const char *const flname)
{
    FILE *stream;
    uint32_t retval;
    if((stream=fopen(flname,"r"))==NULL)
    {
        printf("ERROR\nCan't open %s for reading.\n",flname);
        return 1;
    }
    retval=StreamReadTwoOrbDMFTSettings(O,stream);
    fclose(stream);
    return retval;
}

uint32_t StreamPrintTwoOrbDMFTSettings(TwoOrbDMFTSettings O, FILE *stream)
{
    uint32_t i;
    fprintf(stream,"\nSettigns of self-consistent DMFT calculation\n");
    fprintf(stream,"G0 Input option           : %i\n",O.Opt);
    fprintf(stream,"G0FileNameA               : %s\n",O.G0FileNameA);
    fprintf(stream,"G0FileNameB               : %s\n",O.G0FileNameB);
    fprintf(stream,"GRootFileNameA            : %s\n",O.GRootFileNameA);
    fprintf(stream,"GRootFileNameB            : %s\n",O.GRootFileNameB);
    fprintf(stream,"SpecRootFileName          : %s\n",O.SpecRootFileName);
    fprintf(stream,"WaveFunctionRootFileName  : %s\n",O.WaveFunctionRootFileName);
    fprintf(stream,"iStart                    : %i\n",O.iStart);
    fprintf(stream,"Ni                        : %i\n",O.Ni);
    fprintf(stream,"Dyson perturbation        : %i\n",O.Dyson);
    fprintf(stream,"NGCluster                 : %i\n",O.NGCluster);
    fprintf(stream,"NGBath                    : %i\n",O.NGBath);
    fprintf(stream,"Mix                       : %24.17E\n",O.Mix);
    fprintf(stream,"U                         : %24.17E\n",O.U);
    fprintf(stream,"t'                        : %24.17E\n",O.t1);
    fprintf(stream,"J                         : %24.17E\n",O.J);
    fprintf(stream,"NBathA                    : %i\n",O.NBathA);
    fprintf(stream,"NBathB                    : %i\n",O.NBathB);
    fprintf(stream,"mu1                       : %24.17E\n",O.mu1);
    fprintf(stream,"mu2                       : %24.17E\n",O.mu2);
    fprintf(stream,"Spectra                   : %i\n",O.Spectra);
    if((O.Spectra)&4)
    {
        fprintf(stream,"  NQ                     : %i\n",O.NQ);
        for(i=0;i<O.NQ;i++)
        {
            fprintf(stream,"  Q[%3i]                 : %24.17E\n",i,O.Q[i]);
        }
    }
    fprintf(stream,"NE                        : %i\n",O.NE);
    fprintf(stream,"Emin                      : %24.17E\n",O.Emin);
    fprintf(stream,"Emax                      : %24.17E\n",O.Emax);
    fprintf(stream,"Gamma                     : %24.17E\n",O.Gamma);
    fprintf(stream,"zeroLancEnergy            : %24.17E\n",O.zeroLancEnergy);
    fprintf(stream,"zeroLancWaveFunction      : %24.17E\n",O.zeroLancWavefunction);
    fprintf(stream,"NTriDiag                  : %i\n",O.NTriDiag);
    fprintf(stream,"NBitsKey                  : %i\n",O.NBitsKey);
    return 0;
}

uint32_t PrintTwoOrbDMFTSettings(TwoOrbDMFTSettings O)
{ return StreamPrintTwoOrbDMFTSettings(O,stdout);}

/* Create PES and IPES operators */

uint32_t CreatePESOperatorDMFT(OperatorType *O, uint32_t NFermion, uint16_t i)
{
    uint16_t oadd;
    InitOperator(O);
    O->NFermion=NFermion;
    O->NWord=Min(O->NFermion,1+(O->NFermion-1)/8);
    InitOperatorFixedLength(O, 1, 1);
    oadd=i;
    RealOperatorAddLader(O, &oadd, 1, 1.0);
    return 0;
}

uint32_t CreateIPESOperatorDMFT(OperatorType *O, uint32_t NFermion, uint16_t i)
{
    uint16_t oadd;
    InitOperator(O);
    O->NFermion=NFermion;
    O->NWord=Min(O->NFermion,1+(O->NFermion-1)/8);
    InitOperatorFixedLength(O, 1, 1);
    oadd=CrAnMask+i;
    RealOperatorAddLader(O, &oadd, 1, 1.0);
    return 0;
}

/* Sum over three Green's functions in tri-diagonal representation */
uint32_t SumThreeGreensfunctionsTriDiagonalRepresentation(TriDiagonalMatrixType const A, TriDiagonalMatrixType const B, TriDiagonalMatrixType const C, TriDiagonalMatrixType *const M)
{
    double vecA[A.N+B.N+C.N][A.N];
    double vecB[A.N+B.N+C.N][B.N];
    double vecC[A.N+B.N+C.N][C.N];
    double Atemp,Btemp,Ctemp;
    uint32_t i,j;
    M->N=A.N+B.N+C.N;
    if(InitTriDiagonalMatrix(M)){
        printf("InitTriDiagonalMatrix failed in SumThreeGreensfunctionsTriDiagonalRepresentation\n");fflush(stdout);return 1;
    };
    M->a[0]=A.a[0]+B.a[0]+C.a[0];
    M->b[0]=sqrt(Sqr(A.b[0])+Sqr(B.b[0])+Sqr(C.b[0]));
    for(i=0;i<3;i++){
        for(j=0;j<A.N;j++){
            vecA[i][j]=0;
        }
        for(j=0;j<B.N;j++){
            vecB[i][j]=0;
        }
        for(j=0;j<C.N;j++){
            vecC[i][j]=0;
        }
    }
    if(M->b[0]>10*DBL_MIN){
        vecA[0][0]=A.b[0]/M->b[0];
        vecB[0][0]=B.b[0]/M->b[0];
        vecC[0][0]=C.b[0]/M->b[0];
    }
    for(i=1;i<M->N;i++){
        RealTriDiagonalMatrixVector(A, vecA[(i-1)], vecA[i]);           /* vec[i] = A vec[i-1] */
        RealTriDiagonalMatrixVector(B, vecB[(i-1)], vecB[i]);
        RealTriDiagonalMatrixVector(C, vecC[(i-1)], vecC[i]);
        for(j=i;j>0;){
            j--;
            RealVectorDotVector(vecA[j], vecA[i], A.N, &Atemp);           /* orthogonalization of vec[i] to vec[j<i-1] */
            RealVectorDotVector(vecB[j], vecB[i], B.N, &Btemp);
            RealVectorDotVector(vecC[j], vecC[i], C.N, &Ctemp);
            RealVectorAdd(vecA[i],vecA[j],A.N,-(Atemp+Btemp+Ctemp));
            RealVectorAdd(vecB[i],vecB[j],B.N,-(Atemp+Btemp+Ctemp));
            RealVectorAdd(vecC[i],vecC[j],C.N,-(Atemp+Btemp+Ctemp));
            if(j==i-1){
                M->a[i]=Atemp+Btemp+Ctemp;
            }
        }
        RealVectorDotVector(vecA[i],vecA[i],A.N,&Atemp);                /* b[i] = sqrt( vec[i] . vec[i] ) */
        RealVectorDotVector(vecB[i],vecB[i],B.N,&Btemp);
        RealVectorDotVector(vecC[i],vecC[i],C.N,&Btemp);
        M->b[i]=sqrt(Atemp+Btemp+Ctemp);
        if(M->b[i]>100*(DBL_MIN+M->b[1]*DBL_EPSILON)){                   /* vec[i] *= 1/b[i] */
            RealVectorMultiply(vecA[i],A.N,1.0/M->b[i]);
            RealVectorMultiply(vecB[i],B.N,1.0/M->b[i]);
            RealVectorMultiply(vecC[i],C.N,1.0/M->b[i]);
        }
        else{
            M->N=i;
            M->b[i]=0;
            RealVectorMultiply(vecA[i],A.N,0);
            RealVectorMultiply(vecB[i],B.N,0);
            RealVectorMultiply(vecB[i],C.N,0);
        }
    }
    RealVectorTriDiagonalMatrixVector(A, vecA[(M->N-1)], &Atemp);     /* a[N] = vec[N-1] A vec[N-1] */
    RealVectorTriDiagonalMatrixVector(B, vecB[(M->N-1)], &Btemp);
    RealVectorTriDiagonalMatrixVector(C, vecC[(M->N-1)], &Ctemp);
    M->a[M->N]=Atemp+Btemp+Ctemp;

    return 0;
}


/* Read Poles list from file into ListOfPoles type */





/* copy a TriDiagnoalMatrixType */

uint32_t CopyTriDiagonalMatrix(TriDiagonalMatrixType const O1, TriDiagonalMatrixType *const O2){
    uint32_t i;
    sprintf(O2->Name, "%s", O1.Name);
    O2->N=O1.N;
    O2->a=malloc((O2->N+1)*sizeof(double));
    if(O2->a==NULL) {printf("malloc for O2->a failed in CopyTriDiagonalMatrix\n");fflush(stdout);return 1;};
    O2->b=malloc((O2->N)*sizeof(double));
    if(O2->b==NULL)
    {
        free(O2->a);
        printf("malloc for O2->b failed in CopyTriDiagonalMatrix\n");
        fflush(stdout);
        return 1;
    }
    for(i=0;i<=O2->N;i++)
    {
        O2->a[i]=O1.a[i];
    }
    for(i=0;i<O2->N;i++)
    {
        O2->b[i]=O1.b[i];
    }
    
    return 0;
}