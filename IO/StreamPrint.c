
/*
 *  Prints an array of spectra
 */
uint32_t StreamPrintSpectra(SpectraType *const Spec, const uint32_t NSpec, FILE *stream)
{
  uint32_t i,j,NE;
  fprintf(stream,"#Spectra: %i\n",NSpec);
  fprintf(stream,"Emin______Emax     ");
  for(i=0;i<NSpec;i++)
  {
    fprintf(stream," %22.15E %22.15E",Spec[i].Emin,Spec[i].Emax);
  }
  fprintf(stream,"\n");
  fprintf(stream,"EminPole__EmaxPole ");
  for(i=0;i<NSpec;i++)
  {
    fprintf(stream," %22.15E %22.15E",Spec[i].EminPole,Spec[i].EmaxPole);
  }
  fprintf(stream,"\n");
  fprintf(stream,"dE________Gamma    ");
  for(i=0;i<NSpec;i++)
  {
    fprintf(stream," %22.15E %22.15E",Spec[i].dE,Spec[i].Gamma);
  }
  fprintf(stream,"\n");
  fprintf(stream,"Phase_____0        ");
  for(i=0;i<NSpec;i++)
  {
    fprintf(stream," %22.15E %22.15E",Spec[i].phase,0.0);
  }
  fprintf(stream,"\n");
  for(j=0;j<5;j++)
  {
    fprintf(stream,"a[%1i]______b[%1i]     ",j,j);
    for(i=0;i<NSpec;i++)
    {
      fprintf(stream," %22.15E %22.15E",Spec[i].a[j],Spec[i].b[j]);
    }
    fprintf(stream,"\n");
  }
  fprintf(stream,"Energy             ");
  for(i=0;i<NSpec;i++)
  {
         if(i<10) fprintf(stream,"  Re[%1i]                  Im[%1i]                ",i,i);
    else if(i<100) fprintf(stream,"  Re[%2i]                 Im[%2i]               ",i,i);
    else if(i<1000) fprintf(stream,"  Re[%3i]                Im[%3i]              ",i,i);
    else if(i<10000) fprintf(stream,"  Re[%4i]               Im[%4i]             ",i,i);
  }
  fprintf(stream,"\n");
  NE=Spec[0].N;
  for(i=1;i<NSpec;i++)
  {
    NE=Max(NE,Spec[i].N);
  }
  for(j=0;j<=NE;j++)
  {
    fprintf(stream,"%19.12E",(double)(Spec[0].dE*j+Spec[0].Emin));
    for(i=0;i<NSpec;i++)
    {
      if(j<=Spec[i].N) fprintf(stream," %22.15E %22.15E",Spec[i].Re[j],Spec[i].Im[j]);
      else fprintf(stream," %22.15E %22.15E",0.0,0.0);
    }
    fprintf(stream,"\n");
  }
  return 0;
}


/*
*   Prints the values of the wave function.
 */
uint32_t StreamPrintWaveFunctionPlain(WaveFunctionType psi, FILE *stream)
{
  uint32_t i,j;
  psi.Name[MAXStrLength-1]='\0';
  fprintf(stream,"%s\n",psi.Name);
  fprintf(stream,"%i\n",psi.QComplex);
  fprintf(stream,"%i\n",psi.N);
  fprintf(stream,"%i\n",psi.NFermion);
  fprintf(stream,"%i\n",psi.NBoson);
  if(psi.QComplex)
  {
    for (i=1;i<=psi.N;i++)
    {
      fprintf(stream,"%23.15E %23.15E ",psi.ValueRe[BlokElement(i)],psi.ValueIm[BlokElement(i)]);
      for(j=0;j<psi.NWord;j++)
      {
        fprintf(stream,"%2.2X",psi.Determinant[BlokElementDet(i,psi.NWord,j)]);
      }
      fprintf(stream,"\n");
    }
    fprintf(stream,"\n");  
  }
  else
  {
    for (i=1;i<=psi.N;i++)
    {
      fprintf(stream,"%23.15E ",psi.ValueRe[BlokElement(i)]);
      for(j=0;j<psi.NWord;j++)
      {
        fprintf(stream,"%2.2X",psi.Determinant[BlokElementDet(i,psi.NWord,j)]);
      }
      fprintf(stream,"\n");
    }
    fprintf(stream,"\n");
  }
  return 0;
}


/*
*   Prints the values of the wave function.
 */
uint32_t StreamPrintWaveFunctionPlainBinary(WaveFunctionType psi, FILE *stream)
{
  uint32_t i;
  psi.Name[MAXStrLength-1]='\0';
  fprintf(stream,"%s\n",psi.Name);
  fprintf(stream,"%i\n",psi.QComplex);
  fprintf(stream,"%i\n",psi.N);
  fprintf(stream,"%i\n",psi.NFermion);
  fprintf(stream,"%i\n",psi.NBoson);
  if(psi.QComplex)
  {
    for (i=1;i<=psi.N;i++)
    {
      fwrite(&(psi.ValueRe[BlokElement(i)]), sizeof(double), 1, stream);
      fwrite(&(psi.ValueIm[BlokElement(i)]), sizeof(double), 1, stream);
      fwrite(&(psi.Determinant[BlokElementDet(i,psi.NWord,0)]),sizeof(uint8_t),psi.NWord,stream);
    } 
  }
  else
  {
    for (i=1;i<=psi.N;i++)
    {
      fwrite(&(psi.ValueRe[BlokElement(i)]), sizeof(double), 1, stream);
      fwrite(&(psi.Determinant[BlokElementDet(i,psi.NWord,0)]),sizeof(uint8_t),psi.NWord,stream);
    } 
  }
  return 0;
}



/*
*   Prints the values of the wave function.
 */
uint32_t StreamPrintWaveFunctionReadable(WaveFunctionType psi, FILE *stream)
{
  uint32_t i,j;
  psi.Name[MAXStrLength-1]='\0';
  fprintf(stream,"\nprintout of WaveFunction: %s\n",psi.Name);
  fprintf(stream,"\nGenneral settings\n");
  fprintf(stream,"QComplex      = %10i (Real==0 or Complex==1)\n",psi.QComplex);
  fprintf(stream,"N             = %10i (Number of basis functions used to discribe psi)\n",psi.N);
  fprintf(stream,"NFermion      = %10i (Number of fermions in the one particle basis)\n",psi.NFermion);
  fprintf(stream,"NBoson        = %10i (Number of bosons in the one particle basis)\n",psi.NBoson);
  if(psi.QComplex)
  {
    fprintf(stream,"\n#      pre-factor             +I  pre-factor         Determinant\n");
    for (i=1;i<=psi.N;i++)
    {
      fprintf(stream,"%4i %20.12E       %20.12E       ",i,psi.ValueRe[BlokElement(i)],psi.ValueIm[BlokElement(i)]);
      for(j=0;j<psi.NWord;j++)
      {
        fprintf(stream,"%2.2X",psi.Determinant[BlokElementDet(i,psi.NWord,j)]);
      }
      fprintf(stream,"\n");
    }
    fprintf(stream,"\n");  
  }
  else
  {
    fprintf(stream,"\n#      pre-factor         Determinant\n");
    for (i=1;i<=psi.N;i++)
    {
      fprintf(stream,"%4i %20.12E       ",i,psi.ValueRe[BlokElement(i)]);
      for(j=0;j<psi.NWord;j++)
      {
        fprintf(stream,"%2.2X",psi.Determinant[BlokElementDet(i,psi.NWord,j)]);
      }
      fprintf(stream,"\n");
    }
    fprintf(stream,"\n");
  }
  return 0;
}


/*
*   Prints the values of the wave function.
 */
uint32_t StreamPrintWaveFunctionReadableTruncated(WaveFunctionType psi, double epsilon, FILE *stream)
{
  uint32_t i,j;
  double epssqr=Sqr(epsilon);
  psi.Name[MAXStrLength-1]='\0';
  fprintf(stream,"\nprintout of WaveFunction: %s\n",psi.Name);
  fprintf(stream,"Onle determinants with a norm of the prefactor larger then %22.15E will be printed\n",epsilon); 
  fprintf(stream,"\nGenneral settings\n");
  fprintf(stream,"QComplex      = %10i (Real==0 or Complex==1)\n",psi.QComplex);
  fprintf(stream,"N             = %10i (Number of basis functions used to discribe psi)\n",psi.N);
  fprintf(stream,"NFermion      = %10i (Number of fermions in the one particle basis)\n",psi.NFermion);
  fprintf(stream,"NBoson        = %10i (Number of bosons in the one particle basis)\n",psi.NBoson);
  if(psi.QComplex)
  {
    fprintf(stream,"\n  pre-factor             +I  pre-factor         Determinant\n");
    for (i=1;i<=psi.N;i++)
    {
      if((Sqr(psi.ValueRe[BlokElement(i)])+Sqr(psi.ValueIm[BlokElement(i)]))>epssqr)
      {
        fprintf(stream,"%20.12E       %20.12E       ",psi.ValueRe[BlokElement(i)],psi.ValueIm[BlokElement(i)]);
        for(j=0;j<psi.NWord;j++)
        {
          fprintf(stream,"%2.2X",psi.Determinant[BlokElementDet(i,psi.NWord,j)]);
        }
        fprintf(stream,"\n");
      }
    }
    fprintf(stream,"\n");  
  }
  else
  {
    fprintf(stream,"\n  pre-factor         Determinant\n");
    for (i=1;i<=psi.N;i++)
    {
      if(fabs(psi.ValueRe[BlokElement(i)])>epsilon)
      {
        fprintf(stream,"%20.12E       ",psi.ValueRe[BlokElement(i)]);
        for(j=0;j<psi.NWord;j++)
        {
          fprintf(stream,"%2.2X",psi.Determinant[BlokElementDet(i,psi.NWord,j)]);
        }
        fprintf(stream,"\n");
      }
    }
    fprintf(stream,"\n");
  }
  return 0;
}


/*
 *   Prints the values of the wave function, including information on memory allocation.
 */
uint32_t StreamPrintWaveFunctionFullReadable(WaveFunctionType psi, FILE *stream)
{
  uint32_t i,j;
  psi.Name[MAXStrLength-1]='\0';
  fprintf(stream,"\nFull printout of WaveFunction: %s\n",psi.Name);
  fprintf(stream,"\nGenneral settings\n");
  fprintf(stream,"QComplex      = %10i (Real==0 or Complex==1)\n",psi.QComplex);
  fprintf(stream,"N             = %10i (Number of basis functions used to discribe psi)\n",psi.N);
  fprintf(stream,"NBlok         = %10i (Number of bloks needed to store psi)\n",psi.NBlok);
  fprintf(stream,"NFermion      = %10i (Number of fermions in the one particle basis)\n",psi.NFermion);
  fprintf(stream,"NBoson        = %10i (Number of bosons in the one particle basis)\n",psi.NBoson);
  fprintf(stream,"NWord         = %10i (NFermion/8)+NBoson\n",psi.NWord);
  fprintf(stream,"NBitsKey      = %10i (Size of the key used for the hash table lookup)\n",psi.NBitsKey);
  fprintf(stream,"\nPointer positions (addresses)\n");
  fprintf(stream,"**ValueRe     = %p\n",psi.ValueRe);
  if(psi.QComplex) fprintf(stream,"**ValueIm     = %p\n",psi.ValueIm);
  fprintf(stream,"**Next        = %p\n",psi.Next);
  fprintf(stream,"**Determinant = %p\n",psi.Determinant);
  for (i=0;i<psi.NBlok;i++)
  {
    fprintf(stream,"BlokNr.       = %i\n",i);
    fprintf(stream,"*ValueRe      = %p\n",psi.ValueRe[i]);
    if(psi.QComplex) fprintf(stream,"*ValueIm      = %p\n",psi.ValueIm[i]);
    fprintf(stream,"*Next         = %p\n",psi.Next[i]);
    fprintf(stream,"*Determinant  = %p\n",psi.Determinant[i]);
  }  
  fprintf(stream,"\nNon-zero elements of Hash Table\n");
  for (i=0;i<((uint32_t)1<<(psi.NBitsKey));i++)
  {
    if(psi.KeyList[i]!=0)
    {
      fprintf(stream,"i= %10i , key= %10i\n",i,psi.KeyList[i]);
    }
  }
  if(psi.QComplex)
  {
    fprintf(stream,"\n  pre-factor             +I  pre-factor         Pointer Addres  Next_Hash Hash_Value   Determinant\n");
    for (i=1;i<=psi.N;i++)
    {
      fprintf(stream,"%20.12E       %20.12E       %p %10i %10i   ",psi.ValueRe[BlokElement(i)],psi.ValueIm[BlokElement(i)],&psi.ValueRe[BlokElement(i)],psi.Next[BlokElement(i)],hash(&psi.Determinant[BlokElementDet(i,psi.NWord,0)],psi.NWord,psi.NBitsKey));
      for(j=0;j<psi.NWord;j++)
      {
        fprintf(stream,"%2.2X",psi.Determinant[BlokElementDet(i,psi.NWord,j)]);
      }
      fprintf(stream,"\n");
    }
    fprintf(stream,"\n");  
  }
  else
  {
    fprintf(stream,"\n  pre-factor         Pointer Addres  Next_Hash Hash_Value   Determinant\n");
    for (i=1;i<=psi.N;i++)
    {
      fprintf(stream,"%20.12E       %p %10i %10i   ",psi.ValueRe[BlokElement(i)],&psi.ValueRe[BlokElement(i)],psi.Next[BlokElement(i)],hash(&psi.Determinant[BlokElementDet(i,psi.NWord,0)],psi.NWord,psi.NBitsKey));
      for(j=0;j<psi.NWord;j++)
      {
        fprintf(stream,"%2.2X",psi.Determinant[BlokElementDet(i,psi.NWord,j)]);
      }
      fprintf(stream,"\n");
    }
    fprintf(stream,"\n");
  }
  return 0;
}


/*
*   Prints the values of the operator.
 */
uint32_t StreamPrintOperatorReadable(OperatorType O, FILE *stream)
{
  int32_t i,j,k;
  O.Name[MAXStrLength-1]='\0';
  fprintf(stream,"\nprintout of Operator: %s\n",O.Name);
  fprintf(stream,"\nGenneral settings\n");
  fprintf(stream,"QComplex      = %10i (Real==0 or Complex==1 or Mixed==2)\n",O.QComplex);
  fprintf(stream,"MaxLength     = %10i (largest number of product of lader operators)\n",O.MaxLength);
  fprintf(stream,"NFermion      = %10i (Number of fermions in the one particle basis)\n",O.NFermion);
  fprintf(stream,"NBoson        = %10i (Number of bosons in the one particle basis)\n",O.NBoson);
  for(i=0;i<=O.MaxLength;i++)
  {
    if(O.l[i].N!=0)
    {
      fprintf(stream,"\nOperator of Length %3i\n",i);
      fprintf(stream,"QComplex      = %10i (Real==0 or Complex==1)\n",O.l[i].QComplex);
      fprintf(stream,"N             = %10i (number of operators of length %3i)\n",O.l[i].N,i);
      for(j=0;j<O.l[i].N;j++)
      {
        for(k=i-1;k>=0;k--)
        {
          ((O.l[i].o[j*i+k])&CrAnMask) ? fprintf(stream,"C") : fprintf(stream,"A");
          fprintf(stream,"%3i ",(O.l[i].o[j*i+k])&NotCrAnMask);
        }
        (O.l[i].QComplex) ? fprintf(stream,"| %22.15E %22.15E\n",O.l[i].V[j].Re,O.l[i].V[j].Im) : fprintf(stream,"| %22.15E\n",O.l[i].VRe[j]);
      }
    }
  }
  printf("\n");
  return 0;
}


/*
*   Prints the values of the operator including information on allocation.
 */
uint32_t StreamPrintOperatorFullReadable(OperatorType O, FILE *stream)
{
  int32_t i,j,k;
  O.Name[MAXStrLength-1]='\0';
  fprintf(stream,"\nprintout of Operator: %s\n",O.Name);
  fprintf(stream,"\nGenneral settings\n");
  fprintf(stream,"QComplex      = %10i (Real==0 or Complex==1 or Mixed==2)\n",O.QComplex);
  fprintf(stream,"MaxLength     = %10i (largest number of product of lader operators)\n",O.MaxLength);
  fprintf(stream,"NFermion      = %10i (Number of fermions in the one particle basis)\n",O.NFermion);
  fprintf(stream,"NBoson        = %10i (Number of bosons in the one particle basis)\n",O.NBoson);
  fprintf(stream,"NWord         = %10i (Number of fermions/8+ number of bosons in the one particle basis)\n",O.NWord);
  fprintf(stream,"*l            = %p   (Pointer addres of the Fixed Lenght Operator Type structures)\n",O.l);
  for(i=0;i<=O.MaxLength;i++)
  {
    fprintf(stream,"\nOperator of Length %3i\n",i);
    fprintf(stream,"&l[%3i]       = %p   (Pointer addres of fixed length variable)\n",i,&O.l[i]);
    fprintf(stream,"QComplex      = %10i (Real==0 or Complex==1)\n",O.l[i].QComplex);
    fprintf(stream,"N             = %10i (number of operators of length %3i)\n",O.l[i].N,i);
    fprintf(stream,"NBitsKey      = %10i (Number of bits in the key)\n",O.l[i].NBitsKey);
    fprintf(stream,"NAloc         = %10i (Size of the allocated memory)\n",O.l[i].NAloc);
    fprintf(stream,"*o            = %p   (Pointer addres of the determinant list)\n",O.l[i].o);
    fprintf(stream,"*Next         = %p   (Pointer addres of the list of next operators with the same hash value)\n",O.l[i].Next);
    fprintf(stream,"*KeyList      = %p   (Pointer addres of the KeyList)\n",O.l[i].KeyList);
    fprintf(stream,"*VRe          = %p   (Pointer addres of the list with Real values for real operators)\n",O.l[i].VRe);
    fprintf(stream,"*V            = %p   (Pointer addres of the list of complex values for complex operators)\n",O.l[i].V);
    
    if((i!=0)&&(O.l[i].NAloc!=0))
    {
      fprintf(stream,"\nNon-zero elements of Hash Table\n");
      for (j=0;j<((uint32_t)1<<(O.l[i].NBitsKey));j++)
      {
        if(O.l[i].KeyList[j]!=0)
        {
          fprintf(stream,"i= %10i , key= %10i\n",j,O.l[i].KeyList[j]);
        }
      }
    }
    for(j=0;j<O.l[i].N;j++)
    {
      fprintf(stream,"%4i | ",j+1);
      for(k=i-1;k>=0;k--)
      {
        ((O.l[i].o[j*i+k])&CrAnMask) ? fprintf(stream,"C") : fprintf(stream,"A");
        fprintf(stream,"%3i ",(O.l[i].o[j*i+k])&NotCrAnMask);
      }
      (O.l[i].QComplex) ? fprintf(stream,"| %22.15E %22.15E",O.l[i].V[j].Re,O.l[i].V[j].Im) : fprintf(stream,"| %22.15E",O.l[i].VRe[j]);
      if (i!=0) fprintf(stream," | -> %10i\n",O.l[i].Next[j]);
    }
  }
  return 0;
}


/*
*   Prints a crystal structure
 */
uint32_t StreamPrintCrystalStructure(CrystalStructureType O, FILE *stream)
{
  uint32_t i;
  fprintf(stream,"\nprintout of Crystal Structure\n");
  fprintf(stream,"Unit cell parameters:\n");
  fprintf(stream,"a: %15.7F %15.7F %15.7F\n",O.rUnitCell[0][0],O.rUnitCell[0][1],O.rUnitCell[0][2]);
  fprintf(stream,"b: %15.7F %15.7F %15.7F\n",O.rUnitCell[1][0],O.rUnitCell[1][1],O.rUnitCell[1][2]);
  fprintf(stream,"c: %15.7F %15.7F %15.7F\n",O.rUnitCell[2][0],O.rUnitCell[2][1],O.rUnitCell[2][2]);
  fprintf(stream,"Reciprocal latice:\n");
  fprintf(stream,"a: %15.7F %15.7F %15.7F\n",O.gUnitCell[0][0],O.gUnitCell[0][1],O.gUnitCell[0][2]);
  fprintf(stream,"b: %15.7F %15.7F %15.7F\n",O.gUnitCell[1][0],O.gUnitCell[1][1],O.gUnitCell[1][2]);
  fprintf(stream,"c: %15.7F %15.7F %15.7F\n",O.gUnitCell[2][0],O.gUnitCell[2][1],O.gUnitCell[2][2]);
  fprintf(stream,"Number of attoms %4i\n",O.NAtoms);
  for(i=0;i<O.NAtoms;i++)
  {
    fprintf(stream,"#%4i | %2s | %2i |  %15.7F %15.7F %15.7F\n",i,O.Atom[i].Name,O.Atom[i].AtomNumber,O.Atom[i].Position[0],O.Atom[i].Position[1],O.Atom[i].Position[2]);
  }
  return 0;
}


/*
*   Prints a crystal structure super cell
 */
uint32_t StreamPrintCrystalStructureSuperCell(CrystalStructureSuperCellType O, FILE *stream)
{
  uint32_t i;
  fprintf(stream,"\nprintout of Super Cell\n");
  StreamPrintCrystalStructure(O.Crystal, stream);
  fprintf(stream,"Number of Unit Cells in super cell: %4i\n",O.NCell);
  fprintf(stream,"Super Cell in Cell units:\n");
  fprintf(stream,"a: %4i %4i %4i\n",O.rSuperCell[0][0],O.rSuperCell[0][1],O.rSuperCell[0][2]);
  fprintf(stream,"b: %4i %4i %4i\n",O.rSuperCell[1][0],O.rSuperCell[1][1],O.rSuperCell[1][2]);
  fprintf(stream,"c: %4i %4i %4i\n",O.rSuperCell[2][0],O.rSuperCell[2][1],O.rSuperCell[2][2]);
  fprintf(stream,"Position of Unit Cells:\n");
  for(i=0;i<O.NCell;i++)
  {
    fprintf(stream,"%4i | %4i %4i %4i\n",i,O.CellPos[i][0],O.CellPos[i][1],O.CellPos[i][2]);
  }
  fprintf(stream,"Super Cell Reciprocal lattice vectors in reciprocal Cell units:\n");
  fprintf(stream,"a: %12.7F %12.7F %12.7F\n",O.gSuperCell[0][0],O.gSuperCell[0][1],O.gSuperCell[0][2]);
  fprintf(stream,"b: %12.7F %12.7F %12.7F\n",O.gSuperCell[1][0],O.gSuperCell[1][1],O.gSuperCell[1][2]);
  fprintf(stream,"c: %12.7F %12.7F %12.7F\n",O.gSuperCell[2][0],O.gSuperCell[2][1],O.gSuperCell[2][2]);
  fprintf(stream,"k-vectors included in super Cell:\n");
  for(i=0;i<O.NCell;i++)
  {
    fprintf(stream,"%4i | %12.7F %12.7F %12.7F\n",i,O.kVector[i][0],O.kVector[i][1],O.kVector[i][2]);
  }
  return 0;
}







/*  Prints a matrix
 */
uint32_t StreamPrintMatrix(MatrixType M, FILE *stream)
{
  uint32_t i,j;
  M.Name[MAXStrLength-1]='\0';
  fprintf(stream,"%s =\n",M.Name);
  fprintf(stream,"Real Part of Matrix with dimensions [Ni=%i][Nj=%i] ([Rows][Collums])\n",M.Ni,M.Nj);
  fprintf(stream,"        ");
  for(j=0;j<M.Nj;j++)
  {
     fprintf(stream,"  [%12i]",j);
  }
  fprintf(stream,"\n");
  for(i=0;i<M.Ni;i++)
  {
    fprintf(stream,"[%6i] ",i);
    for(j=0;j<M.Nj;j++)
    {
      fprintf(stream,"%15.8E ",M.VRe[i][j]);
    }
    fprintf(stream,"\n");
  }
  if(M.QComplex)
  {
    fprintf(stream,"Imaginary Part\n");
    fprintf(stream,"        ");
    for(j=0;j<M.Nj;j++)
    {
       fprintf(stream,"  [%12i]",j);
    }
    fprintf(stream,"\n");
    for(i=0;i<M.Ni;i++)
    {
      fprintf(stream,"[%6i] ",i);
      for(j=0;j<M.Nj;j++)
      {
        fprintf(stream,"%15.8E ",M.VIm[i][j]);
      }
      fprintf(stream,"\n");
    }
  }
  fprintf(stream,"\n");
  return 0;
}


/*  Prints a matrix
 */
uint32_t StreamPrintPlainMatrix(MatrixType M, FILE *stream)
{
  uint32_t i,j;
  if(M.QComplex)
  {
    for(i=0;i<M.Ni;i++)
    {
      for(j=0;j<M.Nj;j++)
      {
        fprintf(stream,"%19.12E %19.12E ",M.VRe[i][j],M.VIm[i][j]);
      }
      fprintf(stream,"\n");
    }
  }
  else
  {
    for(i=0;i<M.Ni;i++)
    {
      for(j=0;j<M.Nj;j++)
      {
        fprintf(stream,"%20.13E ",M.VRe[i][j]);
      }
      fprintf(stream,"\n");
    }
  }
  fprintf(stream,"\n");
  return 0;
}


/*  Prints a matrix
 */
uint32_t StreamPrintPlainTriDiagonalMatrix(TriDiagonalMatrixType M, FILE *stream)
{
  uint32_t i;
  fprintf(stream,"%22.15E %22.15E %i\n",M.a[0],M.b[0],M.N);
  for(i=1;i<M.N;i++)
  {
    fprintf(stream,"%22.15E %22.15E\n",M.a[i],M.b[i]);
  }
  fprintf(stream,"%22.15E\n",M.a[M.N]);
  return 0;
}

/*  Prints the structure TightBindingDef
 *
 */
uint32_t StreamPrintTightBindingDef(TightBindingDef O, FILE *stream)
{
  int i,j;
  fprintf(stream,"\nSettigns of a tight binding model: %s\n",O.Name);
  fprintf(stream,"unit cell parameters:\n");
  fprintf(stream,"a: %22.15E %22.15E %22.15E\n",O.Cell.rUnitCell[0][0],O.Cell.rUnitCell[0][1],O.Cell.rUnitCell[0][2]);
  fprintf(stream,"b: %22.15E %22.15E %22.15E\n",O.Cell.rUnitCell[1][0],O.Cell.rUnitCell[1][1],O.Cell.rUnitCell[1][2]);
  fprintf(stream,"c: %22.15E %22.15E %22.15E\n",O.Cell.rUnitCell[2][0],O.Cell.rUnitCell[2][1],O.Cell.rUnitCell[2][2]);
  for(i=0;i<O.Cell.NAtoms;i++)
  {
    fprintf(stream,"Atom  #: %4i %4s   %22.15E %22.15E %22.15E\n",i,O.Cell.Atom[i].Name,O.Cell.Atom[i].Position[0],O.Cell.Atom[i].Position[1],O.Cell.Atom[i].Position[2]);
  }
  for(i=0;i<O.Cell.NAtoms;i++)
  {
    fprintf(stream,"Atom  #: %4i %4s\n",i,O.Cell.Atom[i].Name);
    for(j=0;j<O.Cell.Atom[i].Nshell;j++)
    {
      fprintf(stream,"Shell #: %4i - %4i %4i fold degenerate.\n",i,j,O.Cell.Atom[i].shelldegeneracy[j]);
    }
  }
  fprintf(stream,"Hopping definitions ( %4i )\n",O.nH);
  for(i=0;i<O.nH;i++)
  {
    fprintf(stream,"Hopping from %4i - %4i to %4i - %4i with connecting vector in unit cells: %4i %4i %4i\n",O.H[i].Atompos1,O.H[i].Shellpos1,O.H[i].Atompos2,O.H[i].Shellpos2,O.H[i].ConnectingVector[0],O.H[i].ConnectingVector[1],O.H[i].ConnectingVector[2]);
    StreamPrintMatrix(O.H[i].Hop,stream);
  }
  fprintf(stream,"\n");
  return 0;
}



/*  Prints the settings for a Dimerized spin chain model
 *
 */
uint32_t StreamPrintDimerizedSpinChainSettings(DimerizedSpinChainSettings O, FILE *stream)
{
  fprintf(stream,"\nSettigns of Heisenberg Model\n");
  fprintf(stream,"Nearest Neighbor exchange constant JA:%22.15F\n",O.JA);
  fprintf(stream,"Nearest Neighbor exchange constant JB:%22.15F\n",O.JB);
  fprintf(stream,"Number of dimers in the chain        :    %i\n",O.Ndimer);
  fprintf(stream,"Additional information for numerical settings\n");
  fprintf(stream,"NPsi                 : %i\n",O.NPsi);
  fprintf(stream,"NE                   : %i\n",O.NE);
  fprintf(stream,"Emin                 : %24.17E\n",O.Emin);
  fprintf(stream,"Emax                 : %24.17E\n",O.Emax);
  fprintf(stream,"Gamma                : %24.17E\n",O.Gamma);
  fprintf(stream,"zeroLancEnergy       : %24.17E\n",O.zeroLancEnergy);
  fprintf(stream,"zeroLancWaveFunction : %24.17E\n",O.zeroLancWavefunction);
  fprintf(stream,"NTriDiag             : %i\n",O.NTriDiag);
  fprintf(stream,"NBitsKey             : %i\n",O.NBitsKey);
  fprintf(stream,"SpecFile             : %s\n",O.SpecFile);
  return 0;
}


/*  Prints the settings for a Heisenberg model
 *
 */
uint32_t StreamPrintHeisenbergModelSettings(HeisenbergModelSettings O, FILE *stream)
{
  fprintf(stream,"\nSettigns of Heisenberg Model\n");
  fprintf(stream,"Nearest Neighbor exchange constant J:%22.15F\n",O.J);
  StreamPrintCrystalStructureSuperCell(O.Cluster, stream);
  fprintf(stream,"Additional information for numerical settings\n");
  fprintf(stream,"NPsi                 : %i\n",O.NPsi);
  fprintf(stream,"NE                   : %i\n",O.NE);
  fprintf(stream,"Emin                 : %24.17E\n",O.Emin);
  fprintf(stream,"Emax                 : %24.17E\n",O.Emax);
  fprintf(stream,"Gamma                : %24.17E\n",O.Gamma);
  fprintf(stream,"zeroLancEnergy       : %24.17E\n",O.zeroLancEnergy);
  fprintf(stream,"zeroLancWaveFunction : %24.17E\n",O.zeroLancWavefunction);
  fprintf(stream,"NTriDiag             : %i\n",O.NTriDiag);
  fprintf(stream,"NBitsKey             : %i\n",O.NBitsKey);
  fprintf(stream,"LocalFourierBasis    : %i\n",O.LocalFourierBasis);
  fprintf(stream,"SpecFile             : %s\n",O.SpecFile);
  return 0;
}


/*  Prints the settings for a Heisenberg model
 *
 */
uint32_t StreamPrinttJModelSettings(tJModelSettings O, FILE *stream)
{
  fprintf(stream,"\nSettigns of tJ Model\n");
  fprintf(stream,"Nearest Neighbor exchange constant J:%22.15F\n",O.J);
  fprintf(stream,"Nearest Neighbor hopping constant  t:%22.15F\n",O.t);
  fprintf(stream,"Number of electrons                 : %i\n",O.nElectron);
  StreamPrintCrystalStructureSuperCell(O.Cluster, stream);
  fprintf(stream,"Additional information for numerical settings\n");
  fprintf(stream,"NPsi                 : %i\n",O.NPsi);
  fprintf(stream,"NE                   : %i\n",O.NE);
  fprintf(stream,"Emin                 : %24.17E\n",O.Emin);
  fprintf(stream,"Emax                 : %24.17E\n",O.Emax);
  fprintf(stream,"Gamma                : %24.17E\n",O.Gamma);
  fprintf(stream,"zeroLancEnergy       : %24.17E\n",O.zeroLancEnergy);
  fprintf(stream,"zeroLancWaveFunction : %24.17E\n",O.zeroLancWavefunction);
  fprintf(stream,"NTriDiag             : %i\n",O.NTriDiag);
  fprintf(stream,"NBitsKey             : %i\n",O.NBitsKey);
  fprintf(stream,"LocalFourierBasis    : %i\n",O.LocalFourierBasis);
  fprintf(stream,"SpecFile             : %s\n",O.SpecFile);
  return 0;
}


/*  Prints the settings for a Spin Lader model
 *
 */
uint32_t StreamPrintSpinLaderSettings(SpinLaderSettings O, FILE *stream)
{
  fprintf(stream,"\nSettigns of Spin Lader Model\n");
  fprintf(stream,"Nearest Neighbor exchange constant Jrung:%22.15F\n",O.Jrung);
  fprintf(stream,"Nearest Neighbor exchange constant Jleg: %22.15F\n",O.Jleg);
  fprintf(stream,"Nearest Neighbor exchange constant Jcyc: %22.15F\n",O.Jcyc);
  fprintf(stream,"Number of rungs                        : %i\n",O.NSite);

  fprintf(stream,"Additional information for numerical settings\n");
  fprintf(stream,"NPsi                 : %i\n",O.NPsi);
  fprintf(stream,"NE                   : %i\n",O.NE);
  fprintf(stream,"Emin                 : %24.17E\n",O.Emin);
  fprintf(stream,"Emax                 : %24.17E\n",O.Emax);
  fprintf(stream,"Gamma                : %24.17E\n",O.Gamma);
  fprintf(stream,"zeroLancEnergy       : %24.17E\n",O.zeroLancEnergy);
  fprintf(stream,"zeroLancWaveFunction : %24.17E\n",O.zeroLancWavefunction);
  fprintf(stream,"NTriDiag             : %i\n",O.NTriDiag);
  fprintf(stream,"NBitsKey             : %i\n",O.NBitsKey);
  fprintf(stream,"SpecFile             : %s\n",O.SpecFile);
  return 0;
}


/*  Prints the settings for a Spin Lader model
 *
 */
uint32_t StreamPrintHubbardModelSettings(HubbardModelSettings O, FILE *stream)
{
  uint32_t iatom,ishell,ifermion;
  fprintf(stream,"\nSettigns of Hubbard Model\n");
  fprintf(stream,"Crystal structure and super cell settings\n");
  StreamPrintCrystalStructureSuperCell(O.SuperCell, stream);
  fprintf(stream,"Tigth Binding definitions within unit cell\n");
  StreamPrintTightBindingDef(O.TightBinding, stream);
  fprintf(stream,"Tigth Binding definitions within super cell\n");
  StreamPrintTightBindingDef(O.TightBindingSuperCell, stream);
  fprintf(stream,"Additional Hamiltonian within unit cell\n");
  StreamPrintOperatorReadable(O.U,stream);
  fprintf(stream,"\nAdditional Hamiltonian within Super cell\n");
  StreamPrintOperatorReadable(O.USuperCell,stream);
  fprintf(stream,"\n Filling before diagonalization \n");
  for(iatom=0;iatom<O.SuperCell.Crystal.NAtoms;iatom++)
  {
    fprintf(stream,"atom %2i ",iatom);
    for(ishell=0;ishell<O.SuperCell.Crystal.Atom[iatom].Nshell;ishell++)
    {
      fprintf(stream,"shell %2i |",ishell);
      for(ifermion=0;ifermion<O.SuperCell.Crystal.Atom[iatom].shelldegeneracy[ishell];ifermion++)
      {
        fprintf(stream,"%2i ",O.NElectron[iatom][ishell][ifermion]);
      }
      fprintf(stream,"\n");
    }
    fprintf(stream,"\n");
  }
  fprintf(stream,"Additional information for numerical settings\n");
  fprintf(stream,"Nq                   : %i %i %i\n",O.Nq[0],O.Nq[1],O.Nq[2]);
  fprintf(stream,"NPsi (max)           : %i\n",O.NPsi);
  fprintf(stream,"DeltaEMax            : %lf\n",O.DeltaEMax);
  fprintf(stream,"NE                   : %i\n",O.NE);
  fprintf(stream,"Emin                 : %24.17E\n",O.Emin);
  fprintf(stream,"Emax                 : %24.17E\n",O.Emax);
  fprintf(stream,"Gamma                : %24.17E\n",O.Gamma);
  fprintf(stream,"zeroLancEnergy       : %24.17E\n",O.zeroLancEnergy);
  fprintf(stream,"zeroLancWaveFunction : %24.17E\n",O.zeroLancWavefunction);
  fprintf(stream,"NTriDiag             : %i\n",O.NTriDiag);
  fprintf(stream,"NBitsKey             : %i\n",O.NBitsKey);
  fprintf(stream,"SpecFile             : %s\n",O.SpecFile);
  return 0;
}


/*  Prints the settings for a XAS spectrum of a single site DMFT model
 *
 */
uint32_t StreamPrintDMFTXASSettings(DMFTXASSettings O, FILE *stream)
{
  fprintf(stream,"\nSettigns of XAS on DMFT model\n");
  fprintf(stream,"BathFilename         : %s\n",O.BathFilename);
  fprintf(stream,"U                    : %24.17E\n",O.U);
  fprintf(stream,"Q                    : %24.17E\n",O.Q);
  fprintf(stream,"NBathMax             : %i\n",O.NBathMax);
  fprintf(stream,"NE                   : %i\n",O.NE);
  fprintf(stream,"Emin                 : %24.17E\n",O.Emin);
  fprintf(stream,"Emax                 : %24.17E\n",O.Emax);
  fprintf(stream,"Gamma                : %24.17E\n",O.Gamma);
  fprintf(stream,"zeroLancEnergy       : %24.17E\n",O.zeroLancEnergy);
  fprintf(stream,"zeroLancWaveFunction : %24.17E\n",O.zeroLancWavefunction);
  fprintf(stream,"NTriDiag             : %i\n",O.NTriDiag);
  fprintf(stream,"NBitsKey             : %i\n",O.NBitsKey);
  fprintf(stream,"SpecFile             : %s\n",O.SpecFile);
  return 0;
}

/*
 *
 */
uint32_t StreamPrintDMFTSettings(DMFTSettings O, FILE *stream)
{
  uint32_t i;
  fprintf(stream,"\nSettigns of self-consistent DMFT calculation\n");
  fprintf(stream,"G0FileName               : %s\n",O.G0FileName);  
  fprintf(stream,"GRootFileName            : %s\n",O.GRootFileName);  
  fprintf(stream,"SpecRootFileName         : %s\n",O.SpecRootFileName);  
  fprintf(stream,"WaveFunctionRootFileName : %s\n",O.WaveFunctionRootFileName);  
  fprintf(stream,"iStart                   : %i\n",O.iStart);
  fprintf(stream,"Ni                       : %i\n",O.Ni);
  fprintf(stream,"Dyson perturbation       : %i\n",O.Dyson);
  fprintf(stream,"NGCluster                : %i\n",O.NGCluster);
  fprintf(stream,"NGBath                   : %i\n",O.NGBath);
  fprintf(stream,"Mix                      : %24.17E\n",O.Mix);
  fprintf(stream,"U                        : %24.17E\n",O.U);
  fprintf(stream,"NBath                    : %i\n",O.NBath);
  fprintf(stream,"Spectra                  : %i\n",O.Spectra);
  if((O.Spectra)&4)
  {
    fprintf(stream,"  NQ                     : %i\n",O.NQ);
    for(i=0;i<O.NQ;i++)
    {
      fprintf(stream,"  Q[%3i]                 : %24.17E\n",i,O.Q[i]);
    }
  }
  fprintf(stream,"NE                       : %i\n",O.NE);
  fprintf(stream,"Emin                     : %24.17E\n",O.Emin);
  fprintf(stream,"Emax                     : %24.17E\n",O.Emax);
  fprintf(stream,"Gamma                    : %24.17E\n",O.Gamma);
  fprintf(stream,"zeroLancEnergy           : %24.17E\n",O.zeroLancEnergy);
  fprintf(stream,"zeroLancWaveFunction     : %24.17E\n",O.zeroLancWavefunction);
  fprintf(stream,"NTriDiag                 : %i\n",O.NTriDiag);
  fprintf(stream,"NBitsKey                 : %i\n",O.NBitsKey);
  return 0;
}













/*  Prints a matrix
 */
uint32_t StreamPrintSparseMatrix(SparseMatrixType M, FILE *stream)
{
  uint32_t i;
  M.Name[MAXStrLength-1]='\0';
  fprintf(stream,"%s =\n",M.Name);
  if(M.QComplex)
  {
    fprintf(stream,"   [i]    [j]          Re(M[i][j])          Im(M[i][j])\n");
    for(i=0;i<M.N;i++)
    {
      fprintf(stream,"%6i %6i %22.15E %22.15E\n",M.i[i],M.j[i],M.VRe[i],M.VRe[i]);
    }
  }
  else
  {
    fprintf(stream,"   [i]    [j]          Re(M[i][j])\n");
    for(i=0;i<M.N;i++)
    {
      fprintf(stream,"%6i %6i %22.15E\n",M.i[i],M.j[i],M.VRe[i]);
    }
  }
  fprintf(stream,"\n");
  return 0;
}
