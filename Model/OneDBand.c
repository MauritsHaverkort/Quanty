#include "OneDBand.h"

void ItinerantOneDBand()
{
  WaveFunctionType psi, Opsi;
  OperatorType H,HRot;
  uint8_t det[10];
  uint16_t oadd[2];
  uint32_t i;
  MatrixType Rho, TRho, TTot, temp;
  double ValRho[80];
  
  psi.NFermion=16;
  psi.NBoson=0;
  psi.QComplex=0;
  psi.NBitsKey=20;
  WaveFunctionInit(&psi);
  det[0]=0xFF;
  det[1]=0x00;
  det[2]=0;
  det[3]=0;
  WaveFunctionAddElement(&psi, det, 1.0, 0.0);
  
  InitOperator(&H);
  H.NBoson=psi.NBoson;
  H.NFermion=psi.NFermion;
  H.NWord = psi.NWord;
  InitOperatorFixedLength(&H, 2, 8);
  for(i=1;i<psi.NFermion;i++)
  {
    oadd[0]=(uint16_t)(i-1);
    oadd[1]=CrAnMask | (uint16_t)(i);
    OperatorAddLaderNormalOrder(&H, oadd, 2, 1, 0,&QDetNotOrderedFcFaBcBa);
    oadd[0]=(uint16_t)(i);
    oadd[1]=CrAnMask | (uint16_t)(i-1);
    OperatorAddLaderNormalOrder(&H, oadd, 2, 1, 0,&QDetNotOrderedFcFaBcBa);
  }
  
  
  PrintWaveFunctionReadable(psi);
  TTot.Ni=psi.NFermion;
  TTot.Nj=psi.NFermion;
  TTot.QComplex=0;
  InitIdentityMatrix(&TTot);
  for(i=0;i<10;i++)
  {
    printf("================= %i ============\n",i+1);
    /* rotate H */
    RotateOperatorNormalOrder(&HRot, H, TTot, &QDetNotOrderedFcFaBcBa);
    OperatorChop(&HRot, 100*DBL_EPSILON);
    /* expand psi and calcualte ground state */
    OperatorPsi(HRot, psi, &Opsi);
    WaveFunctionAddMultipliedUnionBasis(&Opsi, psi,0);
    LanczosGroundStateConserveBasisRR(HRot, &Opsi, 0, 0.00000001, 0, &Opsi, 0);
    /* print psi */
    PrintWaveFunctionReadable(Opsi);
    /* calcualte density matrix and aditional rotation of basis */
    CalculateDensityMatrix(&Rho, Opsi);
    Eigensystem(Rho, &TRho, ValRho);
    /* print eigenvalues of ValRho */
    
    //PrintMatrix(TRho);
    
    /* Rotate psi */
    WaveFunctionFree(&psi);
    MatrixTranspose(&TRho);
    RotateWaveFunction(&psi, Opsi, TRho);
    WaveFunctionFree(&Opsi);
    WaveFunctionTruncateBasis(&psi, 0.01);
    printf("Rotated psi\n");
    PrintWaveFunctionReadable(psi);
    
    /* create total rotation until now */
    MatrixMultiply(&temp, TTot, TRho);
    FreeMatrix(&TTot);
    MatrixCopy(&TTot,temp);
    FreeMatrix(&temp);
  }
  
  
}