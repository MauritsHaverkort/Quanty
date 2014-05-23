

/*  Given a many particle operator O and a density matrix Rho
 *  create the mean-field approximated operator MO
 *  a^dag_1 a^dag_2 a_3 a_4 ->
 *   - a^dag_1 a_3 < a^dag_2 a_4 >
 *   + a^dag_1 a_4 < a^dag_2 a_3 >
 *   + a^dag_2 a_3 < a^dag_1 a_4 >
 *   - a^dag_2 a_4 < a^dag_1 a_3 >
 */
uint32_t MeanFieldOperator(const OperatorType * const O, OperatorType * const MO, const MatrixType Rho)
{
  uint32_t i;
  uint16_t oadd[2];
  InitOperatorCopySettings(MO, *O);
  if(NormalOrderOperator(*O, MO, &QDetNotOrderedFcFaBcBa))
  {
    printf("NormalOrderOperator failed in MeanFieldOperator\n");
    return 1;
  }
  if(MO->MaxLength!=4)
  {
    printf("MeanFieldOperator only implemented to reduce a 2 particle operator to 1 particle operator\n");
    return 1;
  }
  /* all real */
  if((Rho.QComplex==0)&&(MO->l[4].QComplex==0)&&(MO->l[2].QComplex==0))
  {
    for(i=0;i<MO->l[4].N;i++)
    {
      oadd[0]=MO->l[4].o[i*4+0];
      oadd[1]=MO->l[4].o[i*4+2];
      RealOperatorAddLaderNormalOrder(MO, oadd, 2,-Rho.VRe[NotCrAnMask & MO->l[4].o[i*4+3]][MO->l[4].o[i*4+1]]*MO->l[4].VRe[i]/2.0, &QDetNotOrderedFcFaBcBa);
      oadd[0]=MO->l[4].o[i*4+0];
      oadd[1]=MO->l[4].o[i*4+3];
      RealOperatorAddLaderNormalOrder(MO, oadd, 2, Rho.VRe[NotCrAnMask & MO->l[4].o[i*4+2]][MO->l[4].o[i*4+1]]*MO->l[4].VRe[i]/2.0, &QDetNotOrderedFcFaBcBa);
      oadd[0]=MO->l[4].o[i*4+1];
      oadd[1]=MO->l[4].o[i*4+2];
      RealOperatorAddLaderNormalOrder(MO, oadd, 2, Rho.VRe[NotCrAnMask & MO->l[4].o[i*4+3]][MO->l[4].o[i*4+0]]*MO->l[4].VRe[i]/2.0, &QDetNotOrderedFcFaBcBa);
      oadd[0]=MO->l[4].o[i*4+1];
      oadd[1]=MO->l[4].o[i*4+3];
      RealOperatorAddLaderNormalOrder(MO, oadd, 2,-Rho.VRe[NotCrAnMask & MO->l[4].o[i*4+2]][MO->l[4].o[i*4+0]]*MO->l[4].VRe[i]/2.0, &QDetNotOrderedFcFaBcBa);
    }
  }
  else
  {
    if(MO->l[2].QComplex==0) MakeRealOperatorFixedLengthComplex(MO, 2);
    if((Rho.QComplex==0)&&(MO->l[4].QComplex==0))
    {
      for(i=0;i<MO->l[4].N;i++)
      {
        oadd[0]=MO->l[4].o[i*4+0];
        oadd[1]=MO->l[4].o[i*4+2];
        ComplexOperatorAddLaderNormalOrder(MO, oadd, 2,-Rho.VRe[NotCrAnMask & MO->l[4].o[i*4+3]][MO->l[4].o[i*4+1]]*MO->l[4].VRe[i]/2.0, 0, &QDetNotOrderedFcFaBcBa);
        oadd[0]=MO->l[4].o[i*4+0];
        oadd[1]=MO->l[4].o[i*4+3];
        ComplexOperatorAddLaderNormalOrder(MO, oadd, 2, Rho.VRe[NotCrAnMask & MO->l[4].o[i*4+2]][MO->l[4].o[i*4+1]]*MO->l[4].VRe[i]/2.0, 0, &QDetNotOrderedFcFaBcBa);
        oadd[0]=MO->l[4].o[i*4+1];
        oadd[1]=MO->l[4].o[i*4+2];
        ComplexOperatorAddLaderNormalOrder(MO, oadd, 2, Rho.VRe[NotCrAnMask & MO->l[4].o[i*4+3]][MO->l[4].o[i*4+0]]*MO->l[4].VRe[i]/2.0, 0, &QDetNotOrderedFcFaBcBa);
        oadd[0]=MO->l[4].o[i*4+1];
        oadd[1]=MO->l[4].o[i*4+3];
        ComplexOperatorAddLaderNormalOrder(MO, oadd, 2,-Rho.VRe[NotCrAnMask & MO->l[4].o[i*4+2]][MO->l[4].o[i*4+0]]*MO->l[4].VRe[i]/2.0, 0, &QDetNotOrderedFcFaBcBa);
      }
    }
    else if((Rho.QComplex==0)&&(MO->l[4].QComplex==1))
    {
      for(i=0;i<MO->l[4].N;i++)
      {
        oadd[0]=MO->l[4].o[i*4+0];
        oadd[1]=MO->l[4].o[i*4+2];
        ComplexOperatorAddLaderNormalOrder(MO, oadd, 2,-Rho.VRe[NotCrAnMask & MO->l[4].o[i*4+3]][MO->l[4].o[i*4+1]]*MO->l[4].V[i].Re/2.0, -Rho.VRe[NotCrAnMask & MO->l[4].o[i*4+3]][MO->l[4].o[i*4+1]]*MO->l[4].V[i].Im/2.0, &QDetNotOrderedFcFaBcBa);
        oadd[0]=MO->l[4].o[i*4+0];
        oadd[1]=MO->l[4].o[i*4+3];
        ComplexOperatorAddLaderNormalOrder(MO, oadd, 2, Rho.VRe[NotCrAnMask & MO->l[4].o[i*4+2]][MO->l[4].o[i*4+1]]*MO->l[4].V[i].Re/2.0,  Rho.VRe[NotCrAnMask & MO->l[4].o[i*4+2]][MO->l[4].o[i*4+1]]*MO->l[4].V[i].Im/2.0, &QDetNotOrderedFcFaBcBa);
        oadd[0]=MO->l[4].o[i*4+1];
        oadd[1]=MO->l[4].o[i*4+2];
        ComplexOperatorAddLaderNormalOrder(MO, oadd, 2, Rho.VRe[NotCrAnMask & MO->l[4].o[i*4+3]][MO->l[4].o[i*4+0]]*MO->l[4].V[i].Re/2.0,  Rho.VRe[NotCrAnMask & MO->l[4].o[i*4+3]][MO->l[4].o[i*4+0]]*MO->l[4].V[i].Im/2.0, &QDetNotOrderedFcFaBcBa);
        oadd[0]=MO->l[4].o[i*4+1];
        oadd[1]=MO->l[4].o[i*4+3];
        ComplexOperatorAddLaderNormalOrder(MO, oadd, 2,-Rho.VRe[NotCrAnMask & MO->l[4].o[i*4+2]][MO->l[4].o[i*4+0]]*MO->l[4].V[i].Re/2.0, -Rho.VRe[NotCrAnMask & MO->l[4].o[i*4+2]][MO->l[4].o[i*4+0]]*MO->l[4].V[i].Im/2.0, &QDetNotOrderedFcFaBcBa);
      }
    }
    else if((Rho.QComplex==1)&&(MO->l[4].QComplex==0))
    {
      for(i=0;i<MO->l[4].N;i++)
      {
        oadd[0]=MO->l[4].o[i*4+0];
        oadd[1]=MO->l[4].o[i*4+2];
        ComplexOperatorAddLaderNormalOrder(MO, oadd, 2,-Rho.VRe[NotCrAnMask & MO->l[4].o[i*4+3]][MO->l[4].o[i*4+1]]*MO->l[4].VRe[i]/2.0, -Rho.VIm[NotCrAnMask & MO->l[4].o[i*4+3]][MO->l[4].o[i*4+1]]*MO->l[4].VRe[i]/2.0, &QDetNotOrderedFcFaBcBa);
        oadd[0]=MO->l[4].o[i*4+0];
        oadd[1]=MO->l[4].o[i*4+3];
        ComplexOperatorAddLaderNormalOrder(MO, oadd, 2, Rho.VRe[NotCrAnMask & MO->l[4].o[i*4+2]][MO->l[4].o[i*4+1]]*MO->l[4].VRe[i]/2.0,  Rho.VIm[NotCrAnMask & MO->l[4].o[i*4+2]][MO->l[4].o[i*4+1]]*MO->l[4].VRe[i]/2.0, &QDetNotOrderedFcFaBcBa);
        oadd[0]=MO->l[4].o[i*4+1];
        oadd[1]=MO->l[4].o[i*4+2];
        ComplexOperatorAddLaderNormalOrder(MO, oadd, 2, Rho.VRe[NotCrAnMask & MO->l[4].o[i*4+3]][MO->l[4].o[i*4+0]]*MO->l[4].VRe[i]/2.0,  Rho.VIm[NotCrAnMask & MO->l[4].o[i*4+3]][MO->l[4].o[i*4+0]]*MO->l[4].VRe[i]/2.0, &QDetNotOrderedFcFaBcBa);
        oadd[0]=MO->l[4].o[i*4+1];
        oadd[1]=MO->l[4].o[i*4+3];
        ComplexOperatorAddLaderNormalOrder(MO, oadd, 2,-Rho.VRe[NotCrAnMask & MO->l[4].o[i*4+2]][MO->l[4].o[i*4+0]]*MO->l[4].VRe[i]/2.0, -Rho.VIm[NotCrAnMask & MO->l[4].o[i*4+2]][MO->l[4].o[i*4+0]]*MO->l[4].VRe[i]/2.0, &QDetNotOrderedFcFaBcBa);
      }
    }
    else if((Rho.QComplex==1)&&(MO->l[4].QComplex==1))
    {
      for(i=0;i<MO->l[4].N;i++)
      {
        oadd[0]=MO->l[4].o[i*4+0];
        oadd[1]=MO->l[4].o[i*4+2];
        ComplexOperatorAddLaderNormalOrder(MO, oadd, 2,-Rho.VRe[NotCrAnMask & MO->l[4].o[i*4+3]][MO->l[4].o[i*4+1]]*MO->l[4].V[i].Re/2.0 + Rho.VIm[NotCrAnMask & MO->l[4].o[i*4+3]][MO->l[4].o[i*4+1]]*MO->l[4].V[i].Im/2.0, -Rho.VRe[NotCrAnMask & MO->l[4].o[i*4+3]][MO->l[4].o[i*4+1]]*MO->l[4].V[i].Im/2.0 -Rho.VIm[NotCrAnMask & MO->l[4].o[i*4+3]][MO->l[4].o[i*4+1]]*MO->l[4].V[i].Re/2.0, &QDetNotOrderedFcFaBcBa);
        oadd[0]=MO->l[4].o[i*4+0];
        oadd[1]=MO->l[4].o[i*4+3];
        ComplexOperatorAddLaderNormalOrder(MO, oadd, 2, Rho.VRe[NotCrAnMask & MO->l[4].o[i*4+2]][MO->l[4].o[i*4+1]]*MO->l[4].V[i].Re/2.0 - Rho.VIm[NotCrAnMask & MO->l[4].o[i*4+2]][MO->l[4].o[i*4+1]]*MO->l[4].V[i].Im/2.0,  Rho.VRe[NotCrAnMask & MO->l[4].o[i*4+2]][MO->l[4].o[i*4+1]]*MO->l[4].V[i].Im/2.0 +Rho.VIm[NotCrAnMask & MO->l[4].o[i*4+2]][MO->l[4].o[i*4+1]]*MO->l[4].V[i].Re/2.0, &QDetNotOrderedFcFaBcBa);
        oadd[0]=MO->l[4].o[i*4+1];
        oadd[1]=MO->l[4].o[i*4+2];
        ComplexOperatorAddLaderNormalOrder(MO, oadd, 2, Rho.VRe[NotCrAnMask & MO->l[4].o[i*4+3]][MO->l[4].o[i*4+0]]*MO->l[4].V[i].Re/2.0 - Rho.VIm[NotCrAnMask & MO->l[4].o[i*4+3]][MO->l[4].o[i*4+0]]*MO->l[4].V[i].Im/2.0,  Rho.VRe[NotCrAnMask & MO->l[4].o[i*4+3]][MO->l[4].o[i*4+0]]*MO->l[4].V[i].Im/2.0 +Rho.VIm[NotCrAnMask & MO->l[4].o[i*4+3]][MO->l[4].o[i*4+0]]*MO->l[4].V[i].Re/2.0, &QDetNotOrderedFcFaBcBa);
        oadd[0]=MO->l[4].o[i*4+1];
        oadd[1]=MO->l[4].o[i*4+3];
        ComplexOperatorAddLaderNormalOrder(MO, oadd, 2,-Rho.VRe[NotCrAnMask & MO->l[4].o[i*4+2]][MO->l[4].o[i*4+0]]*MO->l[4].V[i].Re/2.0 + Rho.VIm[NotCrAnMask & MO->l[4].o[i*4+2]][MO->l[4].o[i*4+0]]*MO->l[4].V[i].Im/2.0, -Rho.VRe[NotCrAnMask & MO->l[4].o[i*4+2]][MO->l[4].o[i*4+0]]*MO->l[4].V[i].Im/2.0 -Rho.VIm[NotCrAnMask & MO->l[4].o[i*4+2]][MO->l[4].o[i*4+0]]*MO->l[4].V[i].Re/2.0, &QDetNotOrderedFcFaBcBa);
      }
    }
    
  }
  
  
  for (i = 3; i <= MO->MaxLength; i++)
  {
    if (MO->l[i].NAloc)
    {
      free(MO->l[i].o);
      free(MO->l[i].Next);
      free(MO->l[i].KeyList);
      (MO->l[i].QComplex) ? free(MO->l[i].V) : free(MO->l[i].VRe);
    }
  }
  MO->MaxLength=2;
  OperatorChop(MO, DBL_MIN*100);
  SetQComplexOperator(MO);
  return 0;
}


/*
 *
 */
uint32_t OperatorToOneElectronMatrix(const OperatorType * const O, MatrixType * const Oone)
{
  uint32_t i,j;
  Oone->Ni=O->NFermion;
  Oone->Nj=O->NFermion;
  Oone->QComplex=O->l[2].QComplex;
  InitMatrix(Oone);
  if(Oone->QComplex==0)
  {
    for(i=0;i<O->l[0].N;i++)
    {
      for(j=0;j<Oone->Ni;j++)
      {
        Oone->VRe[j][j]+=O->l[0].VRe[i];
      }
    }
    if(O->MaxLength<2) return 0;
    for(i=0;i<O->l[2].N;i++)
    {
      if(((CrAnMask & O->l[2].o[2*i]))&&(!(CrAnMask & O->l[2].o[1+2*i])))
      {
        Oone->VRe[NotCrAnMask & O->l[2].o[2*i]][O->l[2].o[1+2*i]]+=O->l[2].VRe[i];
      }
      if((!(CrAnMask & O->l[2].o[2*i]))&&((CrAnMask & O->l[2].o[1+2*i])))
      {
        Oone->VRe[NotCrAnMask & O->l[2].o[1+2*i]][O->l[2].o[2*i]]+=O->l[2].VRe[i];
      }
    }
  }
  else
  {
    for(i=0;i<O->l[0].N;i++)
    {
      for(j=0;j<Oone->Ni;j++)
      {
        Oone->VRe[j][j]+=O->l[0].VRe[i];
      }
    }
    if(O->MaxLength<2) return 0;
    for(i=0;i<O->l[2].N;i++)
    {
      if(((CrAnMask & O->l[2].o[2*i]))&&(!(CrAnMask & O->l[2].o[1+2*i])))
      {
        Oone->VRe[NotCrAnMask & O->l[2].o[2*i]][O->l[2].o[1+2*i]]+=O->l[2].V[i].Re;
        Oone->VIm[NotCrAnMask & O->l[2].o[2*i]][O->l[2].o[1+2*i]]+=O->l[2].V[i].Im;
      }
      if((!(CrAnMask & O->l[2].o[2*i]))&&((CrAnMask & O->l[2].o[1+2*i])))
      {
        Oone->VRe[NotCrAnMask & O->l[2].o[1+2*i]][O->l[2].o[2*i]]+=O->l[2].V[i].Re;
        Oone->VIm[NotCrAnMask & O->l[2].o[1+2*i]][O->l[2].o[2*i]]+=O->l[2].V[i].Im;
      }
    }
  }
  
  return 0;
}
