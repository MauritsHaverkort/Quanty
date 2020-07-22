/*  This function changes a compact matrix into a one particel operator
 *
 *  state |0> is chaged by the operator Cr1 An0 to state |1>
 *  the matrix
 *    ( 0 0 )
 *    ( 1 0 )
 * acting on the vector
 *   (1 0)
 * becomes the vector
 *   (0 1)
 * if i is the index for the different rows and j for the different columns, Cr1 An0 is thus found at the position i=1 j=0
 * or at position M(i*Nj+j) in compact form
 *   is the position M[i][j] in extended form
 *
 */
uint32_t CompactMatrixToOneParticleOperator(const CompactMatrixType * const Omat, OperatorType * const O)
{
  uint16_t oadd[2];
  if(Omat->Ni!=Omat->Nj)
  {
    printf("Error in CompactMatrixToOneParticleOperator, expected a square matrix got dimensions %i %i\n",Omat->Ni,Omat->Nj);
    return 1;
  }
  uint32_t i,j;
  InitOperator(O);
  O->NFermion=Omat->Ni;
  O->NBoson=0;
  O->NWord = Min(O->NFermion, 1 + (O->NFermion - 1) / 8)+O->NBoson;
  O->QComplex = 0;
  InitOperatorFixedLength(O,2,8);
  if(Omat->QComplex)
  {
    MakeRealOperatorFixedLengthComplex(O, 2);
    for(i=0;i<Omat->Ni;i++)
    {
      oadd[1]=CrAnMask + (uint16_t)i;
      for(j=0;j<Omat->Nj;j++)
      {
        oadd[0]=(uint16_t)j;
        if(Sqr(Omat->V[i*Omat->Nj+j].Re)+Sqr(Omat->V[i*Omat->Nj+j].Im)>10*DBL_MIN) ComplexOperatorAddLader(O, oadd, 2, Omat->V[i*Omat->Nj+j].Re, Omat->V[i*Omat->Nj+j].Im);
      }
    }
  }
  else
  {
    for(i=0;i<Omat->Ni;i++)
    {
      oadd[1]=CrAnMask + (uint16_t)i;
      for(j=0;j<Omat->Nj;j++)
      {
        oadd[0]=(uint16_t)j;
        if(Sqr(Omat->VRe[i*Omat->Nj+j])>10*DBL_MIN) RealOperatorAddLader(O, oadd, 2, Omat->VRe[i*Omat->Nj+j]);
      }
    }
  }
  return 0;
}



/*
 *
 */
uint32_t OperatorToOneElectronMatrix(const OperatorType * const O, MatrixType * const Oone)
{
  //TODO: correct form matrix change and element 0
  exit(9);
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


/*  state |0> is chaged by the operator Cr1 An0 to state |1>
 *  the matrix
 *    ( 0 0 )
 *    ( 1 0 )
 * acting on the vector
 *   (1 0)
 * becomes the vector 
 *   (0 1)
 * if i is the index for the different rows and j for the different columns, Cr1 An0 is thus found at the position i=1 j=0
 * or at position M(i*Nj+j) in compact form 
 *   is the position M[i][j] in extended form
 *
 * This function is "unsafe" as it does not check if the arrays excists it enters
 * O.l[2] has to be ordered in the form that Cr operators are left (pos 1) of annihilation operators (pos 0)
 */
uint32_t OperatorCrAnLength2ToOneElectronCompactMatrix(const OperatorType * const O, CompactMatrixType * const Oone)
{
  uint32_t i,j,k;
  Oone->Ni=O->NFermion;
  Oone->Nj=O->NFermion;
  Oone->QComplex=O->l[2].QComplex;
  if(InitCompactMatrix(Oone))
  {
    printf("InitCompactMatrix failed in OperatorCrAnLength2ToOneElectronCompactMatrix\n");
    return 1;
  }
  if(Oone->QComplex)
  {
    for(k=0;k<O->l[2].N;k++)
    {
      i=(uint32_t)(NotCrAnMask & O->l[2].o[2*k+1]);
      j=(uint32_t)(O->l[2].o[2*k]);
      Oone->V[i*Oone->Nj+j].Re+=O->l[2].V[k].Re;
      Oone->V[i*Oone->Nj+j].Im+=O->l[2].V[k].Im;
    }
  }
  else
  {
    for(k=0;k<O->l[2].N;k++)
    {
      i=(uint32_t)(NotCrAnMask & O->l[2].o[2*k+1]);
      j=(uint32_t)(O->l[2].o[2*k]);
      Oone->VRe[i*Oone->Nj+j]+=O->l[2].VRe[k];
    }
  }
  return 0;
}

/*  Given a many particle operator O and a density matrix Rho
 *  create the mean-field approximated operator MO as a one particle matrix
 *  a^dag_i a^dag_j a_k a_l ->
 *  -  a^dag_i a_k  < a^dag_j a_l >
 *  +  a^dag_i a_l  < a^dag_j a_k >
 *  +  a^dag_j a_k  < a^dag_i a_l >
 *  -  a^dag_j a_l  < a^dag_i a_k >
 *  -< a^dag_i a_l >< a^dag_j a_k >
 *  +< a^dag_i a_k >< a^dag_j a_l >
 */
uint32_t MeanFieldOperatorOneElectronMatrix(const OperatorType * const O, CompactMatrixType * const MO, const CompactMatrixType Rho)
{
  uint32_t i,j,k,l,m,n;
  double HFRe, HFIm;
  Complex16 Nelec={0,0};
  OperatorType Oo;
  if((O->MaxLength>4)||((O->MaxLength>0)&&(O->l[1].N!=0))||((O->MaxLength>2)&&(O->l[3].N!=0)))
  {
    printf("ERROR: MeanFieldOperator only implemented to reduce a 2, 1 and 0 particle operator to a 1 particle operator\n");
    return 1;
  }
  if(Rho.Ni!=Rho.Nj)
  {
    printf("ERROR: Density matrix should be a square matrix in MeanFieldOperatorOneElectronMatrix\n");
    return 1;
  }
  if(Rho.Ni!=O->NFermion)
  {
    printf("ERROR: Density matrix should have dimensions of NFermion by NFermion of operator in MeanFieldOperatorOneElectronMatrix\n");
    return 1;
  }
  /* initiallize */
  InitOperatorCopySettings(&Oo, *O);
  NormalOrderOperator(*O, &Oo, &QDetNotOrderedFcFaBcBa);
  MO->Ni=O->NFermion;
  MO->Nj=O->NFermion;
  MO->QComplex=Rho.QComplex;
  if(O->QComplex) MO->QComplex=1;
  InitCompactMatrix(MO);
  strncpy(MO->Name, O->Name, MAXStrLength);
  MO->Name[MAXStrLength-1]='\0';
  /* calculate number of electrons */
  if(Rho.QComplex)
  {
    for(i=0;i<Rho.Ni;i++)
    {
      Nelec.Re+=Rho.V[i*Rho.Nj+i].Re;
      Nelec.Im+=Rho.V[i*Rho.Nj+i].Im;
    }
  }
  else
  {
    for(i=0;i<Rho.Ni;i++)
    {
      Nelec.Re+=Rho.VRe[i*Rho.Nj+i];
    }
  }
  if(Abs(Nelec.Im)<Abs(DBL_EPSILON*10*Nelec.Re)) Nelec.Im=0;
  /* length 0 */
  if((Sqr(Nelec.Re)+Sqr(Nelec.Im))>DBL_EPSILON)
  {
    if(MO->QComplex)
    {
      for(i=0;i<MO->Ni;i++)
      {
        if(Oo.l[0].QComplex)
        {
          for(j=0;j<Oo.l[0].N;j++)
          {
            MO->V[i*MO->Nj+i].Re+=(Oo.l[0].V[j].Re*Nelec.Re+Oo.l[0].V[j].Im*Nelec.Im)/(Sqr(Nelec.Re)+Sqr(Nelec.Im));
            MO->V[i*MO->Nj+i].Im+=(Oo.l[0].V[j].Im*Nelec.Re-Oo.l[0].V[j].Re*Nelec.Im)/(Sqr(Nelec.Re)+Sqr(Nelec.Im));
          }
        }
        else
        {
          for(j=0;j<Oo.l[0].N;j++)
          {
            MO->V[i*MO->Nj+i].Re+=(Oo.l[0].VRe[j]*Nelec.Re)/(Sqr(Nelec.Re)+Sqr(Nelec.Im));
            MO->V[i*MO->Nj+i].Im-=(Oo.l[0].VRe[j]*Nelec.Im)/(Sqr(Nelec.Re)+Sqr(Nelec.Im));
          }
        }
      }
    }
    else // all real
    {
      for(i=0;i<MO->Ni;i++)
      {
        for(j=0;j<Oo.l[0].N;j++)
        {
          MO->VRe[i*MO->Nj+i]+=(Oo.l[0].VRe[0])/(Nelec.Re);
        }
      }
    }
  }
  if(Oo.MaxLength==0)
  {
    FreeOperator(&Oo);
    return 0;
  }
  /* length 2 */
  if(Oo.l[2].QComplex)
  {
    for(k=0;k<Oo.l[2].N;k++)
    {
      i=(uint32_t)(NotCrAnMask & Oo.l[2].o[2*k+1]);
      j=(uint32_t)(Oo.l[2].o[2*k]);
      MO->V[i*MO->Nj+j].Re+=Oo.l[2].V[k].Re;
      MO->V[i*MO->Nj+j].Im+=Oo.l[2].V[k].Im;
    }
  }
  else
  {
    if(MO->QComplex)
    {
      for(k=0;k<Oo.l[2].N;k++)
      {
        i=(uint32_t)(NotCrAnMask & Oo.l[2].o[2*k+1]);
        j=(uint32_t)(Oo.l[2].o[2*k]);
        MO->V[i*MO->Nj+j].Re+=Oo.l[2].VRe[k];
      }
    }
    else
    {
      for(k=0;k<Oo.l[2].N;k++)
      {
        i=(uint32_t)(NotCrAnMask & Oo.l[2].o[2*k+1]);
        j=(uint32_t)(Oo.l[2].o[2*k]);
        MO->VRe[i*MO->Nj+j]+=Oo.l[2].VRe[k];
      }
    }
  }
  if(Oo.MaxLength==2)
  {
    FreeOperator(&Oo);
    return 0;
  }
  /* length 4 */
  // Cr_i Cr_j An_k An_l ->
  //  -  a^dag_i a_k  < a^dag_j a_l >
  //  +  a^dag_i a_l  < a^dag_j a_k >
  //  +  a^dag_j a_k  < a^dag_i a_l >
  //  -  a^dag_j a_l  < a^dag_i a_k >
  //  -< a^dag_i a_l >< a^dag_j a_k >
  //  +< a^dag_i a_k >< a^dag_j a_l >
  if(Oo.l[4].QComplex)
  {
    if(Rho.QComplex)
    {
      for(m=0;m<Oo.l[4].N;m++)
      {
        i=(uint32_t)(NotCrAnMask & Oo.l[4].o[4*m+3]);
        j=(uint32_t)(NotCrAnMask & Oo.l[4].o[4*m+2]);
        k=(uint32_t)(Oo.l[4].o[4*m+1]);
        l=(uint32_t)(Oo.l[4].o[4*m]);
        MO->V[i*MO->Nj+k].Re-= Oo.l[4].V[m].Re * Rho.V[j*Rho.Nj+l].Re - Oo.l[4].V[m].Im * Rho.V[j*Rho.Nj+l].Im;
        MO->V[i*MO->Nj+l].Re+= Oo.l[4].V[m].Re * Rho.V[j*Rho.Nj+k].Re - Oo.l[4].V[m].Im * Rho.V[j*Rho.Nj+k].Im;
        MO->V[j*MO->Nj+k].Re+= Oo.l[4].V[m].Re * Rho.V[i*Rho.Nj+l].Re - Oo.l[4].V[m].Im * Rho.V[i*Rho.Nj+l].Im;
        MO->V[j*MO->Nj+l].Re-= Oo.l[4].V[m].Re * Rho.V[i*Rho.Nj+k].Re - Oo.l[4].V[m].Im * Rho.V[i*Rho.Nj+k].Im;
        MO->V[i*MO->Nj+k].Im-= Oo.l[4].V[m].Im * Rho.V[j*Rho.Nj+l].Re + Oo.l[4].V[m].Re * Rho.V[j*Rho.Nj+l].Im;
        MO->V[i*MO->Nj+l].Im+= Oo.l[4].V[m].Im * Rho.V[j*Rho.Nj+k].Re + Oo.l[4].V[m].Re * Rho.V[j*Rho.Nj+k].Im;
        MO->V[j*MO->Nj+k].Im+= Oo.l[4].V[m].Im * Rho.V[i*Rho.Nj+l].Re + Oo.l[4].V[m].Re * Rho.V[i*Rho.Nj+l].Im;
        MO->V[j*MO->Nj+l].Im-= Oo.l[4].V[m].Im * Rho.V[i*Rho.Nj+k].Re + Oo.l[4].V[m].Re * Rho.V[i*Rho.Nj+k].Im;
        
        HFRe = ( Rho.V[i*Rho.Nj+k].Re * Rho.V[j*Rho.Nj+l].Re
                -Rho.V[i*Rho.Nj+l].Re * Rho.V[j*Rho.Nj+k].Re
                -Rho.V[i*Rho.Nj+k].Im * Rho.V[j*Rho.Nj+l].Im
                +Rho.V[i*Rho.Nj+l].Im * Rho.V[j*Rho.Nj+k].Im ) * Oo.l[4].V[m].Re
              -( Rho.V[i*Rho.Nj+k].Im * Rho.V[j*Rho.Nj+l].Re
                -Rho.V[i*Rho.Nj+l].Im * Rho.V[j*Rho.Nj+k].Re
                +Rho.V[i*Rho.Nj+k].Re * Rho.V[j*Rho.Nj+l].Im
                -Rho.V[i*Rho.Nj+l].Re * Rho.V[j*Rho.Nj+k].Im ) * Oo.l[4].V[m].Im;
        HFIm = ( Rho.V[i*Rho.Nj+k].Re * Rho.V[j*Rho.Nj+l].Re
                -Rho.V[i*Rho.Nj+l].Re * Rho.V[j*Rho.Nj+k].Re
                -Rho.V[i*Rho.Nj+k].Im * Rho.V[j*Rho.Nj+l].Im
                +Rho.V[i*Rho.Nj+l].Im * Rho.V[j*Rho.Nj+k].Im ) * Oo.l[4].V[m].Im
              +( Rho.V[i*Rho.Nj+k].Im * Rho.V[j*Rho.Nj+l].Re
                -Rho.V[i*Rho.Nj+l].Im * Rho.V[j*Rho.Nj+k].Re
                +Rho.V[i*Rho.Nj+k].Re * Rho.V[j*Rho.Nj+l].Im
                -Rho.V[i*Rho.Nj+l].Re * Rho.V[j*Rho.Nj+k].Im ) * Oo.l[4].V[m].Re;
        
        if((Sqr(Nelec.Re)+Sqr(Nelec.Im))>DBL_EPSILON)
        {
          for (n=0; n<MO->Ni; n++)
          {
            MO->V[n*MO->Nj+n].Re += (HFRe*Nelec.Re+HFIm*Nelec.Im)/(Sqr(Nelec.Re)+Sqr(Nelec.Im));
            MO->V[n*MO->Nj+n].Im += (HFIm*Nelec.Re-HFRe*Nelec.Im)/(Sqr(Nelec.Re)+Sqr(Nelec.Im));
          }
        }
      }
    }
    else // Rho real Oo complex MO complex
    {
      for(m=0;m<Oo.l[4].N;m++)
      {
        i=(uint32_t)(NotCrAnMask & Oo.l[4].o[4*m+3]);
        j=(uint32_t)(NotCrAnMask & Oo.l[4].o[4*m+2]);
        k=(uint32_t)(Oo.l[4].o[4*m+1]);
        l=(uint32_t)(Oo.l[4].o[4*m]);
        MO->V[i*MO->Nj+k].Re-= Oo.l[4].V[m].Re * Rho.VRe[j*Rho.Nj+l];
        MO->V[i*MO->Nj+l].Re+= Oo.l[4].V[m].Re * Rho.VRe[j*Rho.Nj+k];
        MO->V[j*MO->Nj+k].Re+= Oo.l[4].V[m].Re * Rho.VRe[i*Rho.Nj+l];
        MO->V[j*MO->Nj+l].Re-= Oo.l[4].V[m].Re * Rho.VRe[i*Rho.Nj+k];
        MO->V[i*MO->Nj+k].Im-= Oo.l[4].V[m].Im * Rho.VRe[j*Rho.Nj+l];
        MO->V[i*MO->Nj+l].Im+= Oo.l[4].V[m].Im * Rho.VRe[j*Rho.Nj+k];
        MO->V[j*MO->Nj+k].Im+= Oo.l[4].V[m].Im * Rho.VRe[i*Rho.Nj+l];
        MO->V[j*MO->Nj+l].Im-= Oo.l[4].V[m].Im * Rho.VRe[i*Rho.Nj+k];
        HFRe = ( Rho.VRe[i*Rho.Nj+k] * Rho.VRe[j*Rho.Nj+l]
                -Rho.VRe[i*Rho.Nj+l] * Rho.VRe[j*Rho.Nj+k]) * Oo.l[4].V[m].Re;
        HFIm = ( Rho.VRe[i*Rho.Nj+k] * Rho.VRe[j*Rho.Nj+l]
                -Rho.VRe[i*Rho.Nj+l] * Rho.VRe[j*Rho.Nj+k]) * Oo.l[4].V[m].Im;
      
        if((Sqr(Nelec.Re)+Sqr(Nelec.Im))>DBL_EPSILON)
        {
          for (n=0; n<MO->Ni; n++)
          {
            MO->V[n*MO->Nj+n].Re += (HFRe*Nelec.Re+HFIm*Nelec.Im)/(Sqr(Nelec.Re)+Sqr(Nelec.Im));
            MO->V[n*MO->Nj+n].Im += (HFIm*Nelec.Re-HFRe*Nelec.Im)/(Sqr(Nelec.Re)+Sqr(Nelec.Im));
          }
        }
      }
    }
  }
  else // Oo real
  {
    if(MO->QComplex)
    {
      if(Rho.QComplex)
      {
        for(m=0;m<Oo.l[4].N;m++)
        {
          i=(uint32_t)(NotCrAnMask & Oo.l[4].o[4*m+3]);
          j=(uint32_t)(NotCrAnMask & Oo.l[4].o[4*m+2]);
          k=(uint32_t)(Oo.l[4].o[4*m+1]);
          l=(uint32_t)(Oo.l[4].o[4*m]);
          MO->V[i*MO->Nj+k].Re-= Oo.l[4].VRe[m] * Rho.V[j*Rho.Nj+l].Re;
          MO->V[i*MO->Nj+l].Re+= Oo.l[4].VRe[m] * Rho.V[j*Rho.Nj+k].Re;
          MO->V[j*MO->Nj+k].Re+= Oo.l[4].VRe[m] * Rho.V[i*Rho.Nj+l].Re;
          MO->V[j*MO->Nj+l].Re-= Oo.l[4].VRe[m] * Rho.V[i*Rho.Nj+k].Re;
          MO->V[i*MO->Nj+k].Im-= Oo.l[4].VRe[m] * Rho.V[j*Rho.Nj+l].Im;
          MO->V[i*MO->Nj+l].Im+= Oo.l[4].VRe[m] * Rho.V[j*Rho.Nj+k].Im;
          MO->V[j*MO->Nj+k].Im+= Oo.l[4].VRe[m] * Rho.V[i*Rho.Nj+l].Im;
          MO->V[j*MO->Nj+l].Im-= Oo.l[4].VRe[m] * Rho.V[i*Rho.Nj+k].Im;
          HFRe = ( Rho.V[i*Rho.Nj+k].Re * Rho.V[j*Rho.Nj+l].Re
                  -Rho.V[i*Rho.Nj+l].Re * Rho.V[j*Rho.Nj+k].Re
                  -Rho.V[i*Rho.Nj+k].Im * Rho.V[j*Rho.Nj+l].Im
                  +Rho.V[i*Rho.Nj+l].Im * Rho.V[j*Rho.Nj+k].Im ) * Oo.l[4].VRe[m];
          HFIm = ( Rho.V[i*Rho.Nj+k].Im * Rho.V[j*Rho.Nj+l].Re
                  -Rho.V[i*Rho.Nj+l].Im * Rho.V[j*Rho.Nj+k].Re
                  +Rho.V[i*Rho.Nj+k].Re * Rho.V[j*Rho.Nj+l].Im
                  -Rho.V[i*Rho.Nj+l].Re * Rho.V[j*Rho.Nj+k].Im ) * Oo.l[4].VRe[m];
          if((Sqr(Nelec.Re)+Sqr(Nelec.Im))>DBL_EPSILON)
          {
            for (n=0; n<MO->Ni; n++)
            {
              MO->V[n*MO->Nj+n].Re += (HFRe*Nelec.Re+HFIm*Nelec.Im)/(Sqr(Nelec.Re)+Sqr(Nelec.Im));
              MO->V[n*MO->Nj+n].Im += (HFIm*Nelec.Re-HFRe*Nelec.Im)/(Sqr(Nelec.Re)+Sqr(Nelec.Im));
            }
          }
        }
      }
      else // Rho real Oo real MO complex
      {
        for(m=0;m<Oo.l[4].N;m++)
        {
          i=(uint32_t)(NotCrAnMask & Oo.l[4].o[4*m+3]);
          j=(uint32_t)(NotCrAnMask & Oo.l[4].o[4*m+2]);
          k=(uint32_t)(Oo.l[4].o[4*m+1]);
          l=(uint32_t)(Oo.l[4].o[4*m]);
          MO->V[i*MO->Nj+k].Re-= Oo.l[4].VRe[m] * Rho.VRe[j*Rho.Nj+l];
          MO->V[i*MO->Nj+l].Re+= Oo.l[4].VRe[m] * Rho.VRe[j*Rho.Nj+k];
          MO->V[j*MO->Nj+k].Re+= Oo.l[4].VRe[m] * Rho.VRe[i*Rho.Nj+l];
          MO->V[j*MO->Nj+l].Re-= Oo.l[4].VRe[m] * Rho.VRe[i*Rho.Nj+k];
          
          HFRe = ( Rho.VRe[i*Rho.Nj+k] * Rho.VRe[j*Rho.Nj+l]
                  -Rho.VRe[i*Rho.Nj+l] * Rho.VRe[j*Rho.Nj+k]) * Oo.l[4].VRe[m];
          if((Sqr(Nelec.Re)+Sqr(Nelec.Im))>DBL_EPSILON)
          {
            for (n=0; n<MO->Ni; n++)
            {
              MO->V[n*MO->Nj+n].Re += HFRe*Nelec.Re/(Sqr(Nelec.Re)+Sqr(Nelec.Im));
              MO->V[n*MO->Nj+n].Im -= HFRe*Nelec.Im/(Sqr(Nelec.Re)+Sqr(Nelec.Im));
            }
          }
        }
      }
    }
    else // all real
    {
      for(m=0;m<Oo.l[4].N;m++)
      {
        i=(uint32_t)(NotCrAnMask & Oo.l[4].o[4*m+3]);
        j=(uint32_t)(NotCrAnMask & Oo.l[4].o[4*m+2]);
        k=(uint32_t)(Oo.l[4].o[4*m+1]);
        l=(uint32_t)(Oo.l[4].o[4*m]);
        MO->VRe[i*MO->Nj+k]-= Oo.l[4].VRe[m] * Rho.VRe[j*Rho.Nj+l];
        MO->VRe[i*MO->Nj+l]+= Oo.l[4].VRe[m] * Rho.VRe[j*Rho.Nj+k];
        MO->VRe[j*MO->Nj+k]+= Oo.l[4].VRe[m] * Rho.VRe[i*Rho.Nj+l];
        MO->VRe[j*MO->Nj+l]-= Oo.l[4].VRe[m] * Rho.VRe[i*Rho.Nj+k];
        if(Sqr(Nelec.Re)>DBL_EPSILON)
        {
          for (n=0; n<MO->Ni; n++)
          {
            MO->VRe[n*MO->Nj+n] += ( Rho.VRe[i*Rho.Nj+k] * Rho.VRe[j*Rho.Nj+l] -Rho.VRe[i*Rho.Nj+l] * Rho.VRe[j*Rho.Nj+k]) * Oo.l[4].VRe[m]/(Nelec.Re);
          }
        }
      }
    }
  }
  FreeOperator(&Oo);
  return 0;
}
/*  In DFT electrons interact with themselves. We here add this spurious interaction to the mean-field HF equations as well.
 *  i.e. sum_{ijkl} U_{ijkl} a^dag_i a^dag_j a_k a_l ->
 *  sum_{ijkl} U_{ijkl}
 *  ( -  a^dag_i a_k  < a^dag_j a_l >
 *    +  a^dag_i a_l  < a^dag_j a_k >
 *    +  a^dag_j a_k  < a^dag_i a_l >
 *    -  a^dag_j a_l  < a^dag_i a_k >
 *    -< a^dag_i a_l >< a^dag_j a_k >
 *    +< a^dag_i a_k >< a^dag_j a_l > )
 *
 *    + sum_{i} (sum_{mn}^{NN} (U_{mnnm} - U_{mnmn}) )/(N(N-1)) a^dag_i a_i  < a^dag_i a_i >
 */
uint32_t MeanFieldOperatorOneElectronMatrixWithDFTInteraction(const OperatorType * const O, CompactMatrixType * const MO, const CompactMatrixType Rho)
{
  uint32_t i,j,k,l,m,n;
  uint16_t oadd[4];
  double HFRe, HFIm;
  Complex16 Nelec={0,0}, U={0,0};
  OperatorType Oo;
  if((O->MaxLength>4)||((O->MaxLength>0)&&(O->l[1].N!=0))||((O->MaxLength>2)&&(O->l[3].N!=0)))
  {
    printf("ERROR: MeanFieldOperator only implemented to reduce a 2, 1 and 0 particle operator to a 1 particle operator\n");
    return 1;
  }
  if(Rho.Ni!=Rho.Nj)
  {
    printf("ERROR: Density matrix should be a square matrix in MeanFieldOperatorOneElectronMatrixWithDFTInteraction\n");
    return 1;
  }
  if(Rho.Ni!=O->NFermion)
  {
    printf("ERROR: Density matrix should have dimensions of NFermion by NFermion of operator in MeanFieldOperatorOneElectronMatrixWithDFTInteraction\n");
    return 1;
  }
  /* initiallize */
  InitOperatorCopySettings(&Oo, *O);
  NormalOrderOperator(*O, &Oo, &QDetNotOrderedFcFaBcBa);
  MO->Ni=O->NFermion;
  MO->Nj=O->NFermion;
  MO->QComplex=Rho.QComplex;
  if(O->QComplex) MO->QComplex=1;
  InitCompactMatrix(MO);
  strncpy(MO->Name, O->Name, MAXStrLength);
  MO->Name[MAXStrLength-1]='\0';
  /* calculate number of electrons */
  if(Rho.QComplex)
  {
    for(i=0;i<Rho.Ni;i++)
    {
      Nelec.Re+=Rho.V[i*Rho.Nj+i].Re;
      Nelec.Im+=Rho.V[i*Rho.Nj+i].Im;
    }
  }
  else
  {
    for(i=0;i<Rho.Ni;i++)
    {
      Nelec.Re+=Rho.VRe[i*Rho.Nj+i];
    }
  }
  if(Abs(Nelec.Im)<Abs(DBL_EPSILON*10*Nelec.Re)) Nelec.Im=0;
  /* length 0 */
  if((Sqr(Nelec.Re)+Sqr(Nelec.Im))>DBL_EPSILON)
  {
    if(MO->QComplex)
    {
      for(i=0;i<MO->Ni;i++)
      {
        if(Oo.l[0].QComplex)
        {
          for(j=0;j<Oo.l[0].N;j++)
          {
            MO->V[i*MO->Nj+i].Re+=(Oo.l[0].V[j].Re*Nelec.Re+Oo.l[0].V[j].Im*Nelec.Im)/(Sqr(Nelec.Re)+Sqr(Nelec.Im));
            MO->V[i*MO->Nj+i].Im+=(Oo.l[0].V[j].Im*Nelec.Re-Oo.l[0].V[j].Re*Nelec.Im)/(Sqr(Nelec.Re)+Sqr(Nelec.Im));
          }
        }
        else
        {
          for(j=0;j<Oo.l[0].N;j++)
          {
            MO->V[i*MO->Nj+i].Re+=(Oo.l[0].VRe[j]*Nelec.Re)/(Sqr(Nelec.Re)+Sqr(Nelec.Im));
            MO->V[i*MO->Nj+i].Im-=(Oo.l[0].VRe[j]*Nelec.Im)/(Sqr(Nelec.Re)+Sqr(Nelec.Im));
          }
        }
      }
    }
    else // all real
    {
      for(i=0;i<MO->Ni;i++)
      {
        for(j=0;j<Oo.l[0].N;j++)
        {
          MO->VRe[i*MO->Nj+i]+=(Oo.l[0].VRe[0])/(Nelec.Re);
        }
      }
    }
  }
  if(Oo.MaxLength==0)
  {
    FreeOperator(&Oo);
    return 0;
  }
  /* length 2 */
  if(Oo.l[2].QComplex)
  {
    for(k=0;k<Oo.l[2].N;k++)
    {
      i=(uint32_t)(NotCrAnMask & Oo.l[2].o[2*k+1]);
      j=(uint32_t)(Oo.l[2].o[2*k]);
      MO->V[i*MO->Nj+j].Re+=Oo.l[2].V[k].Re;
      MO->V[i*MO->Nj+j].Im+=Oo.l[2].V[k].Im;
    }
  }
  else
  {
    if(MO->QComplex)
    {
      for(k=0;k<Oo.l[2].N;k++)
      {
        i=(uint32_t)(NotCrAnMask & Oo.l[2].o[2*k+1]);
        j=(uint32_t)(Oo.l[2].o[2*k]);
        MO->V[i*MO->Nj+j].Re+=Oo.l[2].VRe[k];
      }
    }
    else
    {
      for(k=0;k<Oo.l[2].N;k++)
      {
        i=(uint32_t)(NotCrAnMask & Oo.l[2].o[2*k+1]);
        j=(uint32_t)(Oo.l[2].o[2*k]);
        MO->VRe[i*MO->Nj+j]+=Oo.l[2].VRe[k];
      }
    }
  }
  if(Oo.MaxLength==2)
  {
    FreeOperator(&Oo);
    return 0;
  }
  /* length 4 */
  // Cr_i Cr_j An_k An_l ->
  //  -  a^dag_i a_k  < a^dag_j a_l >
  //  +  a^dag_i a_l  < a^dag_j a_k >
  //  +  a^dag_j a_k  < a^dag_i a_l >
  //  -  a^dag_j a_l  < a^dag_i a_k >
  //  -< a^dag_i a_l >< a^dag_j a_k >
  //  +< a^dag_i a_k >< a^dag_j a_l >
  //
  //  with the DFT correction being a^dag_i a_j <a^dag_j a_i>
  if(Oo.l[4].QComplex)
  {
    /* calcualte U as the average between all electrons needed for the DFT correction */
    for(m=0;m<Oo.NFermion;m++)
    {
      for(n=0;n<Oo.NFermion;n++)
      {
        oadd[3] = CrAnMask + (uint16_t)m;
        oadd[2] = CrAnMask + (uint16_t)n;
        oadd[1] =            (uint16_t)n;
        oadd[0] =            (uint16_t)m;
        U += ComplexOperatorPrefactorLader(&Oo,oadd,4);
        oadd[3] = CrAnMask + (uint16_t)m;
        oadd[2] = CrAnMask + (uint16_t)n;
        oadd[1] =            (uint16_t)m;
        oadd[0] =            (uint16_t)n;
        U -= ComplexOperatorPrefactorLader(&Oo,oadd,4);
      }
    }
    if(Oo.NFermion!=1) U/=(0.5*Oo.NFermion*(Oo.NFermion-1));
    if(Rho.QComplex)
    {
      for(m=0;m<Oo.l[4].N;m++)
      {
        i=(uint32_t)(NotCrAnMask & Oo.l[4].o[4*m+3]);
        j=(uint32_t)(NotCrAnMask & Oo.l[4].o[4*m+2]);
        k=(uint32_t)(Oo.l[4].o[4*m+1]);
        l=(uint32_t)(Oo.l[4].o[4*m]);
        MO->V[i*MO->Nj+k].Re-= Oo.l[4].V[m].Re * Rho.V[j*Rho.Nj+l].Re - Oo.l[4].V[m].Im * Rho.V[j*Rho.Nj+l].Im;
        MO->V[i*MO->Nj+l].Re+= Oo.l[4].V[m].Re * Rho.V[j*Rho.Nj+k].Re - Oo.l[4].V[m].Im * Rho.V[j*Rho.Nj+k].Im;
        MO->V[j*MO->Nj+k].Re+= Oo.l[4].V[m].Re * Rho.V[i*Rho.Nj+l].Re - Oo.l[4].V[m].Im * Rho.V[i*Rho.Nj+l].Im;
        MO->V[j*MO->Nj+l].Re-= Oo.l[4].V[m].Re * Rho.V[i*Rho.Nj+k].Re - Oo.l[4].V[m].Im * Rho.V[i*Rho.Nj+k].Im;
        MO->V[i*MO->Nj+k].Im-= Oo.l[4].V[m].Im * Rho.V[j*Rho.Nj+l].Re + Oo.l[4].V[m].Re * Rho.V[j*Rho.Nj+l].Im;
        MO->V[i*MO->Nj+l].Im+= Oo.l[4].V[m].Im * Rho.V[j*Rho.Nj+k].Re + Oo.l[4].V[m].Re * Rho.V[j*Rho.Nj+k].Im;
        MO->V[j*MO->Nj+k].Im+= Oo.l[4].V[m].Im * Rho.V[i*Rho.Nj+l].Re + Oo.l[4].V[m].Re * Rho.V[i*Rho.Nj+l].Im;
        MO->V[j*MO->Nj+l].Im-= Oo.l[4].V[m].Im * Rho.V[i*Rho.Nj+k].Re + Oo.l[4].V[m].Re * Rho.V[i*Rho.Nj+k].Im;
        
        HFRe = ( Rho.V[i*Rho.Nj+k].Re * Rho.V[j*Rho.Nj+l].Re
                -Rho.V[i*Rho.Nj+l].Re * Rho.V[j*Rho.Nj+k].Re
                -Rho.V[i*Rho.Nj+k].Im * Rho.V[j*Rho.Nj+l].Im
                +Rho.V[i*Rho.Nj+l].Im * Rho.V[j*Rho.Nj+k].Im ) * Oo.l[4].V[m].Re
        -( Rho.V[i*Rho.Nj+k].Im * Rho.V[j*Rho.Nj+l].Re
          -Rho.V[i*Rho.Nj+l].Im * Rho.V[j*Rho.Nj+k].Re
          +Rho.V[i*Rho.Nj+k].Re * Rho.V[j*Rho.Nj+l].Im
          -Rho.V[i*Rho.Nj+l].Re * Rho.V[j*Rho.Nj+k].Im ) * Oo.l[4].V[m].Im;
        HFIm = ( Rho.V[i*Rho.Nj+k].Re * Rho.V[j*Rho.Nj+l].Re
                -Rho.V[i*Rho.Nj+l].Re * Rho.V[j*Rho.Nj+k].Re
                -Rho.V[i*Rho.Nj+k].Im * Rho.V[j*Rho.Nj+l].Im
                +Rho.V[i*Rho.Nj+l].Im * Rho.V[j*Rho.Nj+k].Im ) * Oo.l[4].V[m].Im
        +( Rho.V[i*Rho.Nj+k].Im * Rho.V[j*Rho.Nj+l].Re
          -Rho.V[i*Rho.Nj+l].Im * Rho.V[j*Rho.Nj+k].Re
          +Rho.V[i*Rho.Nj+k].Re * Rho.V[j*Rho.Nj+l].Im
          -Rho.V[i*Rho.Nj+l].Re * Rho.V[j*Rho.Nj+k].Im ) * Oo.l[4].V[m].Re;
        
        if((Sqr(Nelec.Re)+Sqr(Nelec.Im))>DBL_EPSILON)
        {
          for (n=0; n<MO->Ni; n++)
          {
            MO->V[n*MO->Nj+n].Re += (HFRe*Nelec.Re+HFIm*Nelec.Im)/(Sqr(Nelec.Re)+Sqr(Nelec.Im));
            MO->V[n*MO->Nj+n].Im += (HFIm*Nelec.Re-HFRe*Nelec.Im)/(Sqr(Nelec.Re)+Sqr(Nelec.Im));
          }
        }
      }
      // DFT correction
      for(i=0;i<O->NFermion;i++)
      {
        for(j=0;j<O->NFermion;j++)
        {
          MO->V[j+i*O->NFermion]+=U.Re*Rho.V[i+j*O->NFermion];
        }
      }
    }
    else // Rho real Oo complex MO complex
    {
      for(m=0;m<Oo.l[4].N;m++)
      {
        i=(uint32_t)(NotCrAnMask & Oo.l[4].o[4*m+3]);
        j=(uint32_t)(NotCrAnMask & Oo.l[4].o[4*m+2]);
        k=(uint32_t)(Oo.l[4].o[4*m+1]);
        l=(uint32_t)(Oo.l[4].o[4*m]);
        MO->V[i*MO->Nj+k].Re-= Oo.l[4].V[m].Re * Rho.VRe[j*Rho.Nj+l];
        MO->V[i*MO->Nj+l].Re+= Oo.l[4].V[m].Re * Rho.VRe[j*Rho.Nj+k];
        MO->V[j*MO->Nj+k].Re+= Oo.l[4].V[m].Re * Rho.VRe[i*Rho.Nj+l];
        MO->V[j*MO->Nj+l].Re-= Oo.l[4].V[m].Re * Rho.VRe[i*Rho.Nj+k];
        MO->V[i*MO->Nj+k].Im-= Oo.l[4].V[m].Im * Rho.VRe[j*Rho.Nj+l];
        MO->V[i*MO->Nj+l].Im+= Oo.l[4].V[m].Im * Rho.VRe[j*Rho.Nj+k];
        MO->V[j*MO->Nj+k].Im+= Oo.l[4].V[m].Im * Rho.VRe[i*Rho.Nj+l];
        MO->V[j*MO->Nj+l].Im-= Oo.l[4].V[m].Im * Rho.VRe[i*Rho.Nj+k];
        HFRe = ( Rho.VRe[i*Rho.Nj+k] * Rho.VRe[j*Rho.Nj+l]
                -Rho.VRe[i*Rho.Nj+l] * Rho.VRe[j*Rho.Nj+k]) * Oo.l[4].V[m].Re;
        HFIm = ( Rho.VRe[i*Rho.Nj+k] * Rho.VRe[j*Rho.Nj+l]
                -Rho.VRe[i*Rho.Nj+l] * Rho.VRe[j*Rho.Nj+k]) * Oo.l[4].V[m].Im;
        
        if((Sqr(Nelec.Re)+Sqr(Nelec.Im))>DBL_EPSILON)
        {
          for (n=0; n<MO->Ni; n++)
          {
            MO->V[n*MO->Nj+n].Re += (HFRe*Nelec.Re+HFIm*Nelec.Im)/(Sqr(Nelec.Re)+Sqr(Nelec.Im));
            MO->V[n*MO->Nj+n].Im += (HFIm*Nelec.Re-HFRe*Nelec.Im)/(Sqr(Nelec.Re)+Sqr(Nelec.Im));
          }
        }
      }
      // DFT correction
      for(i=0;i<O->NFermion;i++)
      {
        for(j=0;j<O->NFermion;j++)
        {
          MO->V[j+i*O->NFermion]+=U*Rho.VRe[i+j*O->NFermion];
        }
      }
    }
  }
  else // Oo real
  {
    /* calcualte U as the average between all electrons needed for the DFT correction */
    for(m=0;m<Oo.NFermion;m++)
    {
      for(n=0;n<Oo.NFermion;n++)
      {
        oadd[3] = CrAnMask + (uint16_t)m;
        oadd[2] = CrAnMask + (uint16_t)n;
        oadd[1] =            (uint16_t)n;
        oadd[0] =            (uint16_t)m;
        U.Re += RealOperatorPrefactorLader(&Oo,oadd,4);
        oadd[3] = CrAnMask + (uint16_t)m;
        oadd[2] = CrAnMask + (uint16_t)n;
        oadd[1] =            (uint16_t)m;
        oadd[0] =            (uint16_t)n;
        U.Re -= RealOperatorPrefactorLader(&Oo,oadd,4);
      }
    }
    if(Oo.NFermion!=1) U.Re/=(0.5*Oo.NFermion*(Oo.NFermion-1));
    if(MO->QComplex)
    {
      if(Rho.QComplex)
      {
        for(m=0;m<Oo.l[4].N;m++)
        {
          i=(uint32_t)(NotCrAnMask & Oo.l[4].o[4*m+3]);
          j=(uint32_t)(NotCrAnMask & Oo.l[4].o[4*m+2]);
          k=(uint32_t)(Oo.l[4].o[4*m+1]);
          l=(uint32_t)(Oo.l[4].o[4*m]);
          MO->V[i*MO->Nj+k].Re-= Oo.l[4].VRe[m] * Rho.V[j*Rho.Nj+l].Re;
          MO->V[i*MO->Nj+l].Re+= Oo.l[4].VRe[m] * Rho.V[j*Rho.Nj+k].Re;
          MO->V[j*MO->Nj+k].Re+= Oo.l[4].VRe[m] * Rho.V[i*Rho.Nj+l].Re;
          MO->V[j*MO->Nj+l].Re-= Oo.l[4].VRe[m] * Rho.V[i*Rho.Nj+k].Re;
          MO->V[i*MO->Nj+k].Im-= Oo.l[4].VRe[m] * Rho.V[j*Rho.Nj+l].Im;
          MO->V[i*MO->Nj+l].Im+= Oo.l[4].VRe[m] * Rho.V[j*Rho.Nj+k].Im;
          MO->V[j*MO->Nj+k].Im+= Oo.l[4].VRe[m] * Rho.V[i*Rho.Nj+l].Im;
          MO->V[j*MO->Nj+l].Im-= Oo.l[4].VRe[m] * Rho.V[i*Rho.Nj+k].Im;
          HFRe = ( Rho.V[i*Rho.Nj+k].Re * Rho.V[j*Rho.Nj+l].Re
                  -Rho.V[i*Rho.Nj+l].Re * Rho.V[j*Rho.Nj+k].Re
                  -Rho.V[i*Rho.Nj+k].Im * Rho.V[j*Rho.Nj+l].Im
                  +Rho.V[i*Rho.Nj+l].Im * Rho.V[j*Rho.Nj+k].Im ) * Oo.l[4].VRe[m];
          HFIm = ( Rho.V[i*Rho.Nj+k].Im * Rho.V[j*Rho.Nj+l].Re
                  -Rho.V[i*Rho.Nj+l].Im * Rho.V[j*Rho.Nj+k].Re
                  +Rho.V[i*Rho.Nj+k].Re * Rho.V[j*Rho.Nj+l].Im
                  -Rho.V[i*Rho.Nj+l].Re * Rho.V[j*Rho.Nj+k].Im ) * Oo.l[4].VRe[m];
          if((Sqr(Nelec.Re)+Sqr(Nelec.Im))>DBL_EPSILON)
          {
            for (n=0; n<MO->Ni; n++)
            {
              MO->V[n*MO->Nj+n].Re += (HFRe*Nelec.Re+HFIm*Nelec.Im)/(Sqr(Nelec.Re)+Sqr(Nelec.Im));
              MO->V[n*MO->Nj+n].Im += (HFIm*Nelec.Re-HFRe*Nelec.Im)/(Sqr(Nelec.Re)+Sqr(Nelec.Im));
            }
          }
        }
        // DFT correction
        for(i=0;i<O->NFermion;i++)
        {
          for(j=0;j<O->NFermion;j++)
          {
            MO->V[j+i*O->NFermion]+=U.Re*Rho.V[i+j*O->NFermion];
          }
        }
      }
      else // Rho real Oo real MO complex
      {
        for(m=0;m<Oo.l[4].N;m++)
        {
          i=(uint32_t)(NotCrAnMask & Oo.l[4].o[4*m+3]);
          j=(uint32_t)(NotCrAnMask & Oo.l[4].o[4*m+2]);
          k=(uint32_t)(Oo.l[4].o[4*m+1]);
          l=(uint32_t)(Oo.l[4].o[4*m]);
          MO->V[i*MO->Nj+k].Re-= Oo.l[4].VRe[m] * Rho.VRe[j*Rho.Nj+l];
          MO->V[i*MO->Nj+l].Re+= Oo.l[4].VRe[m] * Rho.VRe[j*Rho.Nj+k];
          MO->V[j*MO->Nj+k].Re+= Oo.l[4].VRe[m] * Rho.VRe[i*Rho.Nj+l];
          MO->V[j*MO->Nj+l].Re-= Oo.l[4].VRe[m] * Rho.VRe[i*Rho.Nj+k];
          
          HFRe = ( Rho.VRe[i*Rho.Nj+k] * Rho.VRe[j*Rho.Nj+l]
                  -Rho.VRe[i*Rho.Nj+l] * Rho.VRe[j*Rho.Nj+k]) * Oo.l[4].VRe[m];
          if((Sqr(Nelec.Re)+Sqr(Nelec.Im))>DBL_EPSILON)
          {
            for (n=0; n<MO->Ni; n++)
            {
              MO->V[n*MO->Nj+n].Re += HFRe*Nelec.Re/(Sqr(Nelec.Re)+Sqr(Nelec.Im));
              MO->V[n*MO->Nj+n].Im -= HFRe*Nelec.Im/(Sqr(Nelec.Re)+Sqr(Nelec.Im));
            }
          }
        }
        // DFT correction
        for(i=0;i<O->NFermion;i++)
        {
          for(j=0;j<O->NFermion;j++)
          {
            MO->V[j+i*O->NFermion]+=U.Re*Rho.VRe[i+j*O->NFermion];
          }
        }
      }
    }
    else // all real
    {
      for(m=0;m<Oo.l[4].N;m++)
      {
        i=(uint32_t)(NotCrAnMask & Oo.l[4].o[4*m+3]);
        j=(uint32_t)(NotCrAnMask & Oo.l[4].o[4*m+2]);
        k=(uint32_t)(Oo.l[4].o[4*m+1]);
        l=(uint32_t)(Oo.l[4].o[4*m]);
        MO->VRe[i*MO->Nj+k]-= Oo.l[4].VRe[m] * Rho.VRe[j*Rho.Nj+l];
        MO->VRe[i*MO->Nj+l]+= Oo.l[4].VRe[m] * Rho.VRe[j*Rho.Nj+k];
        MO->VRe[j*MO->Nj+k]+= Oo.l[4].VRe[m] * Rho.VRe[i*Rho.Nj+l];
        MO->VRe[j*MO->Nj+l]-= Oo.l[4].VRe[m] * Rho.VRe[i*Rho.Nj+k];
        if(Sqr(Nelec.Re)>DBL_EPSILON)
        {
          for (n=0; n<MO->Ni; n++)
          {
            MO->VRe[n*MO->Nj+n] += ( Rho.VRe[i*Rho.Nj+k] * Rho.VRe[j*Rho.Nj+l] -Rho.VRe[i*Rho.Nj+l] * Rho.VRe[j*Rho.Nj+k]) * Oo.l[4].VRe[m]/(Nelec.Re);
          }
        }
      }
      // DFT correction
      for(i=0;i<O->NFermion;i++)
      {
        for(j=0;j<O->NFermion;j++)
        {
          MO->VRe[j+i*O->NFermion]+=U.Re*Rho.VRe[i+j*O->NFermion];
        }
      }
    }
  }
  FreeOperator(&Oo);
  return 0;
}

uint32_t MeanFieldOperator(const OperatorType * const O, OperatorType * const OMF, const MatrixType Rho, bool AddDFTSelfInteraction = false)
{
	uint32_t i, j, k, l, m;
	double OVRe, OVIm, OMFVRe, OMFVIm;
	uint16_t oadd[2] = { 0, 0 };
	if ((O->MaxLength > 4) || ((O->MaxLength > 0) && (O->l[1].N != 0)) || ((O->MaxLength > 2) && (O->l[3].N != 0)))
	{
		printf("ERROR: MeanFieldOperator only implemented to reduce a 2, 1 and 0 particle operator to a 1 particle operator\n");
		return 1;
	}
	if (Rho.Ni != Rho.Nj)
	{
		printf("ERROR: Density matrix should be a square matrix in MeanFieldOperatorOneElectronMatrix\n");
		return 1;
	}
	if (Rho.Ni != O->NFermion)
	{
		printf("ERROR: Density matrix should have dimensions of NFermion by NFermion of operator in MeanFieldOperatorOneElectronMatrix\n");
		return 1;
	}
	if (O->MaxLength < 3)
	{
		OperatorCopy(OMF, *O);
		return 0;
	}
	/* initiallize */
	InitOperatorCopySettings(OMF, *O);
	OperatorFixedLengthCopy(&(OMF->l[0]), O->l[0], 0);
	OperatorFixedLengthCopy(&(OMF->l[2]), O->l[2], 2);
	strncpy(OMF->Name, O->Name, MAXStrLength); //TODO: In C++, add "Mean Field" to the name.
	OMF->Name[MAXStrLength - 1] = '\0';
	//TODO: should I check if it is indeed a Creation-Annihilation pair?
	if (O->l[4].QComplex)
	{
		Complex16 U;
		U.Re = 0;
		U.Im = 0;
		if (AddDFTSelfInteraction)
		{
			uint16_t otmp[4];
			for (uint16_t m = 0; m < O->NFermion; m++)
			{
				for (uint16_t n = 0; n < O->NFermion; n++)
				{
					otmp[3] = CrAnMask + (uint16_t) m;
					otmp[2] = CrAnMask + (uint16_t) n;
					otmp[1] = (uint16_t) n;
					otmp[0] = (uint16_t) m;
					U += ComplexOperatorPrefactorLader(O, otmp, 4);
					otmp[3] = CrAnMask + (uint16_t) m;
					otmp[2] = CrAnMask + (uint16_t) n;
					otmp[1] = (uint16_t) m;
					otmp[0] = (uint16_t) n;
					U -= ComplexOperatorPrefactorLader(O, otmp, 4);
				}
			}
			if (O->NFermion != 1)
				U /= (0.5 * O->NFermion * (O->NFermion - 1));
		}
		if (Rho.QComplex)
		{
			for (m = 0; m < O->l[4].N; m++)
			{
				i = NotCrAnMask & O->l[4].o[4 * m + 3];
				j = NotCrAnMask & O->l[4].o[4 * m + 2];
				k = O->l[4].o[4 * m + 1];
				l = O->l[4].o[4 * m];
				OVRe = O->l[4].V[m].Re;
				OVIm = O->l[4].V[m].Im;

				oadd[0] = CrAnMask | i;
				oadd[1] = k;
				OMFVRe = -(OVRe * Rho.VRe[j][l] - OVIm * Rho.VIm[j][l]);
				OMFVIm = -(OVIm * Rho.VRe[j][l] + OVRe * Rho.VIm[j][l]);
				OperatorAddLaderNormalOrder(OMF, oadd, 2, OMFVRe, OMFVIm, &QDetNotOrderedFcFaBcBa);

				oadd[0] = CrAnMask | i;
				oadd[1] = l;
				OMFVRe = (OVRe * Rho.VRe[j][k] - OVIm * Rho.VIm[j][k]);
				OMFVIm = (OVIm * Rho.VRe[j][k] + OVRe * Rho.VIm[j][k]);
				OperatorAddLaderNormalOrder(OMF, oadd, 2, OMFVRe, OMFVIm, &QDetNotOrderedFcFaBcBa);

				oadd[0] = CrAnMask | j;
				oadd[1] = k;
				OMFVRe = (OVRe * Rho.VRe[i][l] - OVIm * Rho.VIm[i][l]);
				OMFVIm = (OVIm * Rho.VRe[i][l] + OVRe * Rho.VIm[i][l]);
				OperatorAddLaderNormalOrder(OMF, oadd, 2, OMFVRe, OMFVIm, &QDetNotOrderedFcFaBcBa);

				oadd[0] = CrAnMask | j;
				oadd[1] = l;
				OMFVRe = -(OVRe * Rho.VRe[i][k] - OVIm * Rho.VIm[i][k]);
				OMFVIm = -(OVIm * Rho.VRe[i][k] + OVRe * Rho.VIm[i][k]);
				OperatorAddLaderNormalOrder(OMF, oadd, 2, OMFVRe, OMFVIm, &QDetNotOrderedFcFaBcBa);

				OMFVRe = (Rho.VRe[i][k] * Rho.VRe[j][l] - Rho.VRe[i][l] * Rho.VRe[j][k] - Rho.VIm[i][k] * Rho.VIm[j][l] + Rho.VIm[i][l] * Rho.VIm[j][k]) * OVRe
						- (Rho.VIm[i][k] * Rho.VRe[j][l] - Rho.VIm[i][l] * Rho.VRe[j][k] + Rho.VRe[i][k] * Rho.VIm[j][l] - Rho.VRe[i][l] * Rho.VIm[j][k]) * OVIm;
				OMFVIm = (Rho.VRe[i][k] * Rho.VRe[j][l] - Rho.VRe[i][l] * Rho.VRe[j][k] - Rho.VIm[i][k] * Rho.VIm[j][l] + Rho.VIm[i][l] * Rho.VIm[j][k]) * OVIm
						+ (Rho.VIm[i][k] * Rho.VRe[j][l] - Rho.VIm[i][l] * Rho.VRe[j][k] + Rho.VRe[i][k] * Rho.VIm[j][l] - Rho.VRe[i][l] * Rho.VIm[j][k]) * OVRe;
				OperatorAddLaderNormalOrder(OMF, NULL, 0, OMFVRe, OMFVIm, &QDetNotOrderedFcFaBcBa);
			}
			for (m = 0; m < O->NFermion; ++m)
			{
				oadd[1] = CrAnMask | m;
				oadd[0] = NotCrAnMask & m;
				OperatorAddLaderNormalOrder(OMF, oadd, 2, U.Re * Rho.VRe[m][m] - U.Im * Rho.VIm[m][m], U.Re * Rho.VIm[m][m] + U.Im * Rho.VRe[m][m], &QDetNotOrderedFcFaBcBa);
			}
		} else // Rho real O complex OMF complex
		{
			for (m = 0; m < O->l[4].N; m++)
			{
				i = NotCrAnMask & O->l[4].o[4 * m + 3];
				j = NotCrAnMask & O->l[4].o[4 * m + 2];
				k = O->l[4].o[4 * m + 1];
				l = O->l[4].o[4 * m];
				OVRe = O->l[4].V[m].Re;
				OVIm = O->l[4].V[m].Im;

				oadd[0] = CrAnMask || i;
				oadd[1] = k;
				OMFVRe = -(OVRe * Rho.VRe[j][l]);
				OMFVIm = -(OVIm * Rho.VRe[j][l]);
				OperatorAddLaderNormalOrder(OMF, oadd, 2, OMFVRe, OMFVIm, &QDetNotOrderedFcFaBcBa);

				oadd[0] = CrAnMask | i;
				oadd[1] = l;
				OMFVRe = (OVRe * Rho.VRe[j][k]);
				OMFVIm = (OVIm * Rho.VRe[j][k]);
				OperatorAddLaderNormalOrder(OMF, oadd, 2, OMFVRe, OMFVIm, &QDetNotOrderedFcFaBcBa);

				oadd[0] = CrAnMask | j;
				oadd[1] = k;
				OMFVRe = (OVRe * Rho.VRe[i][l]);
				OMFVIm = (OVIm * Rho.VRe[i][l]);
				OperatorAddLaderNormalOrder(OMF, oadd, 2, OMFVRe, OMFVIm, &QDetNotOrderedFcFaBcBa);

				oadd[0] = CrAnMask | j;
				oadd[1] = l;
				OMFVRe = -(OVRe * Rho.VRe[i][k]);
				OMFVIm = -(OVIm * Rho.VRe[i][k]);
				OperatorAddLaderNormalOrder(OMF, oadd, 2, OMFVRe, OMFVIm, &QDetNotOrderedFcFaBcBa);

				OMFVRe = (Rho.VRe[i][k] * Rho.VRe[j][l] - Rho.VRe[i][l] * Rho.VRe[j][k]) * OVRe;
				OMFVIm = (Rho.VRe[i][k] * Rho.VRe[j][l] - Rho.VRe[i][l] * Rho.VRe[j][k]) * OVIm;
				OperatorAddLaderNormalOrder(OMF, NULL, 0, OMFVRe, OMFVIm, &QDetNotOrderedFcFaBcBa);
			}
			for (m = 0; m < O->NFermion; ++m)
			{
				oadd[1] = CrAnMask | m;
				oadd[0] = NotCrAnMask & m;
				OperatorAddLaderNormalOrder(OMF, oadd, 2, U.Re * Rho.VRe[m][m], U.Im * Rho.VRe[m][m], &QDetNotOrderedFcFaBcBa);
			}
		}
	} else // O->l[4] real
	{
		double U = 0;
		if (AddDFTSelfInteraction)
		{
			uint16_t otmp[4];
			for (uint16_t m = 0; m < O->NFermion; m++)
			{
				for (uint16_t n = 0; n < O->NFermion; n++)
				{
					otmp[3] = CrAnMask + (uint16_t) m;
					otmp[2] = CrAnMask + (uint16_t) n;
					otmp[1] = (uint16_t) n;
					otmp[0] = (uint16_t) m;
					U += RealOperatorPrefactorLader(O, otmp, 4);
					otmp[3] = CrAnMask + (uint16_t) m;
					otmp[2] = CrAnMask + (uint16_t) n;
					otmp[1] = (uint16_t) m;
					otmp[0] = (uint16_t) n;
					U -= RealOperatorPrefactorLader(O, otmp, 4);
				}
			}
			if (O->NFermion != 1)
				U /= (0.5 * O->NFermion * (O->NFermion - 1));
		}
		if (Rho.QComplex)
		{
			for (m = 0; m < O->l[4].N; m++)
			{
				i = NotCrAnMask & O->l[4].o[4 * m + 3];
				j = NotCrAnMask & O->l[4].o[4 * m + 2];
				k = O->l[4].o[4 * m + 1];
				l = O->l[4].o[4 * m];
				OVRe = O->l[4].VRe[m];

				oadd[0] = CrAnMask | i;
				oadd[1] = k;
				OMFVRe = -(OVRe * Rho.VRe[j][l]);
				OMFVIm = -(OVRe * Rho.VIm[j][l]);
				OperatorAddLaderNormalOrder(OMF, oadd, 2, OMFVRe, OMFVIm, &QDetNotOrderedFcFaBcBa);

				oadd[0] = CrAnMask | i;
				oadd[1] = l;
				OMFVRe = (OVRe * Rho.VRe[j][k]);
				OMFVIm = (OVRe * Rho.VIm[j][k]);
				OperatorAddLaderNormalOrder(OMF, oadd, 2, OMFVRe, OMFVIm, &QDetNotOrderedFcFaBcBa);

				oadd[0] = CrAnMask | j;
				oadd[1] = k;
				OMFVRe = (OVRe * Rho.VRe[i][l]);
				OMFVIm = (OVRe * Rho.VIm[i][l]);
				OperatorAddLaderNormalOrder(OMF, oadd, 2, OMFVRe, OMFVIm, &QDetNotOrderedFcFaBcBa);

				oadd[0] = CrAnMask | j;
				oadd[1] = l;
				OMFVRe = -(OVRe * Rho.VRe[i][k]);
				OMFVIm = -(OVRe * Rho.VIm[i][k]);
				OperatorAddLaderNormalOrder(OMF, oadd, 2, OMFVRe, OMFVIm, &QDetNotOrderedFcFaBcBa);

				OMFVRe = (Rho.VRe[i][k] * Rho.VRe[j][l] - Rho.VRe[i][l] * Rho.VRe[j][k] - Rho.VIm[i][k] * Rho.VIm[j][l] + Rho.VIm[i][l] * Rho.VIm[j][k]) * OVRe;
				OMFVIm = (Rho.VIm[i][k] * Rho.VRe[j][l] - Rho.VIm[i][l] * Rho.VRe[j][k] + Rho.VRe[i][k] * Rho.VIm[j][l] - Rho.VRe[i][l] * Rho.VIm[j][k]) * OVRe;
				OperatorAddLaderNormalOrder(OMF, NULL, 0, OMFVRe, OMFVIm, &QDetNotOrderedFcFaBcBa);
			}
			for (m = 0; m < O->NFermion; ++m)
			{
				oadd[1] = CrAnMask | m;
				oadd[0] = NotCrAnMask & m;
				OperatorAddLaderNormalOrder(OMF, oadd, 2, U * Rho.VRe[m][m], U * Rho.VIm[m][m], &QDetNotOrderedFcFaBcBa);
			}
		} else // Rho real O->l[4] real
		{
			for (m = 0; m < O->l[4].N; m++)
			{
				i = NotCrAnMask & O->l[4].o[4 * m + 3];
				j = NotCrAnMask & O->l[4].o[4 * m + 2];
				k = O->l[4].o[4 * m + 1];
				l = O->l[4].o[4 * m];
				OVRe = O->l[4].VRe[m];
				OMFVIm = 0;

				oadd[0] = CrAnMask | i;
				oadd[1] = k;
				OMFVRe = -(OVRe * Rho.VRe[j][l]);
				OperatorAddLaderNormalOrder(OMF, oadd, 2, OMFVRe, OMFVIm, &QDetNotOrderedFcFaBcBa);

				oadd[0] = CrAnMask | i;
				oadd[1] = l;
				OMFVRe = (OVRe * Rho.VRe[j][k]);
				OperatorAddLaderNormalOrder(OMF, oadd, 2, OMFVRe, OMFVIm, &QDetNotOrderedFcFaBcBa);

				oadd[0] = CrAnMask | j;
				oadd[1] = k;
				OMFVRe = (OVRe * Rho.VRe[i][l]);
				OperatorAddLaderNormalOrder(OMF, oadd, 2, OMFVRe, OMFVIm, &QDetNotOrderedFcFaBcBa);

				oadd[0] = CrAnMask | j;
				oadd[1] = l;
				OMFVRe = -(OVRe * Rho.VRe[i][k]);
				OperatorAddLaderNormalOrder(OMF, oadd, 2, OMFVRe, OMFVIm, &QDetNotOrderedFcFaBcBa);

				OMFVRe = (Rho.VRe[i][k] * Rho.VRe[j][l] - Rho.VRe[i][l] * Rho.VRe[j][k]) * OVRe;
				OperatorAddLaderNormalOrder(OMF, NULL, 0, OMFVRe, OMFVIm, &QDetNotOrderedFcFaBcBa);
			}
			for (m = 0; m < O->NFermion; ++m)
			{
				oadd[1] = CrAnMask | m;
				oadd[0] = NotCrAnMask & m;
				OperatorAddLaderNormalOrder(OMF, oadd, 2, U * Rho.VRe[m][m], 0, &QDetNotOrderedFcFaBcBa);
			}
		}
	}
	OperatorChop(OMF, 0);
	return 0;
}

uint32_t MeanFieldOperator(const OperatorType * const O, OperatorType * const OMF, const CompactMatrixType Rho, bool AddDFTSelfInteraction = false)
{
	MatrixType RhoMT;
	CompactMatrixCopyToMatrix(&RhoMT,Rho);
	return MeanFieldOperator(O, OMF, RhoMT, AddDFTSelfInteraction);
}

/*
 *
 */
uint32_t MeanFieldOperatorOld(const OperatorType * const O, OperatorType * const MO, const CompactMatrixType Rho)
{
  CompactMatrixType Omat;
  MeanFieldOperatorOneElectronMatrix(O, &Omat, Rho);
  CompactMatrixToOneParticleOperator(&Omat, MO);
  FreeCompactMatrix(&Omat);
  return 0;
}
uint32_t MeanFieldOperatorWithDFTInteractionOld(const OperatorType * const O, OperatorType * const MO, const CompactMatrixType Rho)
{
  CompactMatrixType Omat;
  MeanFieldOperatorOneElectronMatrixWithDFTInteraction(O, &Omat, Rho);
  CompactMatrixToOneParticleOperator(&Omat, MO);
  FreeCompactMatrix(&Omat);
  return 0;
}

/*  given an eigen-system with val[i] the i-th eigenvalue and fun[i*Nj] the i-th eigen-vector and a total occupation of Nelec
 *
 *  first calcualte the chemical potential
 *
 *  for each eigen-state (k) the density matrix is rho[i*Nj+j] = fun[k*Nj+i]^* * fun[k*Nj+j]
 *
 *  the total density matrix is the sum over all density matrices times the occupation of that state
 *
 */
uint32_t CalculateDensityMatrixFromOneElectronEigensystem(const CompactMatrixType * const fun, double *val, CompactMatrixType *Rho, const double kbT, const double Nelec, double *Etot)
{
  uint32_t i,j,k;
  double occi,mu;
  CalculateChemicalPotentialFromListOfEigenvalues(val, fun->Ni, kbT, Nelec, &mu, Etot);
  Rho->Ni=fun->Ni;
  Rho->Nj=fun->Nj;
  Rho->QComplex=fun->QComplex;
  InitCompactMatrix(Rho);
  if(fun->QComplex)
  {
    for(k=0;k<fun->Ni;k++)
    {
      occi=FermiFunction(val[k], mu, kbT);
      if(occi<DBL_EPSILON) continue;
      for(i=0;i<fun->Ni;i++)
      {
        for(j=0;j<fun->Nj;j++)
        {
          Rho->V[i*Rho->Nj+j].Re+=occi*(fun->V[k*fun->Ni+i].Re*fun->V[k*fun->Ni+j].Re+fun->V[k*fun->Ni+i].Im*fun->V[k*fun->Ni+j].Im);
          Rho->V[i*Rho->Nj+j].Im+=occi*(fun->V[k*fun->Ni+i].Re*fun->V[k*fun->Ni+j].Im-fun->V[k*fun->Ni+i].Im*fun->V[k*fun->Ni+j].Re);
        }
      }
    }
  }
  else
  {
    for(k=0;k<fun->Ni;k++)
    {
      occi=FermiFunction(val[k], mu, kbT);
      if(occi<DBL_EPSILON) continue;
      for(i=0;i<fun->Ni;i++)
      {
        for(j=0;j<fun->Nj;j++)
        {
          Rho->VRe[i*Rho->Nj+j]+=occi*fun->VRe[k*fun->Ni+i]*fun->VRe[k*fun->Ni+j];
        }
      }
    }
  }
  return 0;
}

/*  Given a starting density matrix evolve it to a steady state in mean-field theory for the operator O
 *
 */
uint32_t MeanFieldGroundstate(const OperatorType *const O, CompactMatrixType * const Rho, const double kbT, double * Etot)
{
  uint32_t i,iloop=0, Maxit=1000,somemore=0;
  CompactMatrixType HMF, fun, RhoNew;
  double val[O->NFermion],Etotlist[Maxit],mix=0.8;
  double Nelec=0;
  bool Converged=false;
  // calculate the number of electrons - assuming Rho is hermitian (diagonal is real)
  if(Rho->QComplex)
  {
    for(i=0;i<Rho->Ni;i++)
    {
      Nelec+=Rho->V[i*Rho->Ni+i].Re;
    }
  }
  else
  {
    for(i=0;i<Rho->Ni;i++)
    {
      Nelec+=Rho->VRe[i*Rho->Ni+i];
    }
  }
  // loop till convergence
  while(!Converged)
  {
    if(iloop==100) mix=0.5;
    if(iloop==200) mix=1.0;
    if(iloop==300) mix=0.2;
    if(iloop==400) mix=0.8;
    if(iloop==500) mix=0.1;
    if(iloop==600) mix=0.8;
    MeanFieldOperatorOneElectronMatrix(O, &HMF, *Rho);
    
    CompactMatrixEigensystem(&HMF, &fun, val);
    CalculateDensityMatrixFromOneElectronEigensystem(&fun, val, &RhoNew, kbT, Nelec,Etot);
    
    CompactMatrixScale(Rho, 1-mix);
    CompactMatrixScale(&RhoNew, mix);
    CompactMatrixAdd(Rho, &RhoNew);
    FreeCompactMatrix(&RhoNew);
    
    Etotlist[iloop]=*Etot;
    FreeCompactMatrix(&fun);
    if(iloop!=0)
    {
      if(Abs(Etotlist[iloop]-Etotlist[iloop-1])<(100*((DBL_EPSILON)*Abs(Etotlist[iloop]+Etotlist[iloop-1])+sqrt(DBL_MIN))))
      {
        if(somemore==3)
        {
          Converged=true;
        }
        else
        {
          mix=1.0;
          somemore++;
        }
      }
    }
    VerbosityPrintf(0x0100,"In MeanFieldGroundstate loop %4i total energy is %22.15E dE is %22.15E\n",iloop,Etotlist[iloop],Etotlist[iloop]-Etotlist[(iloop==0)?(0):(iloop-1)]);
    iloop++;
    if(iloop==Maxit)
    {
      printf("Max number of loops reached in MeanFieldGroundstate, but energy not converged\n");
      for(i=0;i<Maxit;i++)
      {
        printf("Energy at step %4i is %22.15E\n",i,Etotlist[i]);
      }
      return 1;
    }
  }
  return 0;
}

uint32_t partialMeanFieldOperator(const OperatorType * const O, OperatorType * const OPartialMF, const MatrixType &Rho, const vector<uint16_t> &indices, bool AddDFTSelfInteraction = false)
{
	OperatorCopy(OPartialMF, *O);
	if (O->MaxLength < 3 || indices.size() < 1)
	{
		return 0;
	}
	OperatorType MeanFieldPartTmp, MeanFieldPart;
	int length = 4; //If ever we can apply mean field to operators of length unequal to 4, this needs to be changed.
	InitOperatorCopySettings(&MeanFieldPart,*OPartialMF);
	InitOperatorCopySettings(&MeanFieldPartTmp,*OPartialMF);
	extractPartialOperatorFixedLength(OPartialMF, &MeanFieldPartTmp, vector<vector<uint16_t> >(1,indices), isCrAnIncluded, 4, true);
	if (MeanFieldPartTmp.l[length].N != 0)
	{
		MeanFieldOperator(&MeanFieldPartTmp, &MeanFieldPart, Rho, AddDFTSelfInteraction);
		OperatorAdd(OPartialMF, MeanFieldPart);
		OperatorChop(OPartialMF, 0);
	}
	SetQComplexOperator(OPartialMF);
	FreeOperator(&MeanFieldPartTmp);
	FreeOperator(&MeanFieldPart);
	return 0;
}
