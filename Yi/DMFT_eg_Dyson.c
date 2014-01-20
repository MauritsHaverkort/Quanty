/* Switch "Dyson" enabled */

#include "./DMFTFunctions.c"

int main(int argc, const char * argv[]){
    
    TwoOrbDMFTSettings EgDMFTSet;
    TimeType Time;
    InitializeTiming(&Time);
    BeginTiming("Total", &Time);
    
    argc-=1;
    TimeStampAll();
    
    /* check if setting file exists. */
    if(FileReadTwoOrbDMFTSettings(&EgDMFTSet, argv[argc])){
        printf("TwoOrbDMFTSettings %s not found.\n", argv[argc]);
        return 1;
    }
    PrintTwoOrbDMFTSettings(EgDMFTSet);
    
    MatrixType T;
    OperatorType HAnd, OPESAup, OPESAdn, OIPESAup, OIPESAdn, OPESBup, OPESBdn, OIPESBup, OIPESBdn, VPerturbationA[4], VPerturbationB[4];
    WaveFunctionType Psi, OPsi, PsiTemp;
    TriDiagonalMatrixType G0A, GBathA, GBathAOld, GBathANew, SigmaA, G0B, GBathB, GBathBOld, GBathBNew, SigmaB, GclAMup, GclAMdn, GclAPup, GclAPdn, GclA, GclBMup, GclBMdn, GclBPup, GclBPdn, GclB, gPerturbationA[4], gPerturbationB[4];
    SecondDoubleTriDiagonalMatrixType GAndA,GAndATrunk,GAndB,GAndBTrunk;
    AndersonMatrixType GBathAAnd, GBathBAnd, GclAAnd, GclBAnd, G0AndA, G0AndB, GBathANewAnd, GBathBNewAnd;
    ListOfPoles G0APoles, G0BPoles, SigmaAPoles, SigmaBPoles, GclAList, GclBList, GListTemp;
    SpectraType Spectra[16];
    char GBathFileNameA[255], GBathFileNameB[255], GclFileNameA[255], GclFileNameB[255], GSigmaFileNameA[255], GSigmaFileNameB[255], PsiFileName[255], SpecFileName[255], Buffer[255];
    uint32_t iOuter, i, iCount=0, PsiNBathVA, PsiNBathCA, PsiNBathVB, PsiNBathCB;
    uint16_t oadd[1];
    double Ef=0., epsilonAnd=0.00001, ec=1e-5, G0shift=0.0, Goccu=1.0, mu, muOld, mutmp;
    double PsiEnergy, PsiQuality, PsiHist[20], *PsiRho, Occu, OccuOld;
    
    /* read pole list */
    if(EgDMFTSet.Opt){
        if(FileReadListOfPoles(&G0APoles, EgDMFTSet.G0FileNameA)){
            printf("Pole list A not found.\n");
            return 1;
        };
        ListOfPolesToTridiagonalMatrix(G0APoles, &G0A, 1000);
        
        if(FileReadListOfPoles(&G0BPoles, EgDMFTSet.G0FileNameB)){
            printf("Pole list B not found.\n");
            return 1;
        };
        ListOfPolesToTridiagonalMatrix(G0BPoles, &G0B, 1000);
    }
    
    /* read G0 */
    else{
        if(FileReadPlainTriDiagonalMatrixType(&G0A, EgDMFTSet.G0FileNameA)){
            printf("G0 file A not found.\n");
            return 1;
        }
        if(FileReadPlainTriDiagonalMatrixType(&G0B, EgDMFTSet.G0FileNameB)){
            printf("G0 file B not found.\n");
            return 1;
        }
    }
    strcpy(G0A.Name, "G0A");
    strcpy(G0B.Name, "G0B");
    
    ShiftTriDiagonalMatrix(&G0A, G0shift);
    ShiftTriDiagonalMatrix(&G0B, G0shift);
    
    TriDiagonalMatrixToAndersonMatrix(G0A, &T, &G0AndA);
    TriDiagonalMatrixToAndersonMatrix(G0B, &T, &G0AndB);
    
    /* read G_Bath_iStart and psi_iStart */
    if(EgDMFTSet.iStart==0)
    {
        sprintf(GBathFileNameA, "%s000", EgDMFTSet.GRootFileNameA);
        sprintf(GBathFileNameB, "%s000", EgDMFTSet.GRootFileNameB);
        sprintf(PsiFileName, "%s000", EgDMFTSet.WaveFunctionRootFileName);
    }
    else if(EgDMFTSet.iStart<10)
    {
        sprintf(GBathFileNameA, "%s00%1i", EgDMFTSet.GRootFileNameA, EgDMFTSet.iStart);
        sprintf(GBathFileNameB, "%s00%1i", EgDMFTSet.GRootFileNameB, EgDMFTSet.iStart);
        sprintf(PsiFileName, "%s00%1i", EgDMFTSet.WaveFunctionRootFileName, EgDMFTSet.iStart);
    }
    else if(EgDMFTSet.iStart<100)
    {
        sprintf(GBathFileNameA, "%s0%2i", EgDMFTSet.GRootFileNameA, EgDMFTSet.iStart);
        sprintf(GBathFileNameB, "%s0%2i", EgDMFTSet.GRootFileNameB, EgDMFTSet.iStart);
        sprintf(PsiFileName, "%s0%2i", EgDMFTSet.WaveFunctionRootFileName, EgDMFTSet.iStart);
    }
    else
    {
        sprintf(GBathFileNameA, "%s%i", EgDMFTSet.GRootFileNameA, EgDMFTSet.iStart);
        sprintf(GBathFileNameB, "%s%i", EgDMFTSet.GRootFileNameB, EgDMFTSet.iStart);
        sprintf(PsiFileName, "%s%i", EgDMFTSet.WaveFunctionRootFileName, EgDMFTSet.iStart);
    }
    
    if(FileReadPlainTriDiagonalMatrixType(&GBathA,GBathFileNameA)){
        printf("Failed to read GBath A from file %s !\n", GBathFileNameA);
        return 1;
    }
    if(FileReadPlainTriDiagonalMatrixType(&GBathB,GBathFileNameB)){
        printf("Failed to read GBath B from file %s !\n", GBathFileNameA);
        return 1;
    }
    
    strcpy(GBathA.Name, "GbathA");
    strcpy(GBathB.Name, "GbathB");
    
    /* read or create starting wave function */
    
    if(file_exists(PsiFileName))
    {
        Psi.NBitsKey=EgDMFTSet.NBitsKey;
        if(FileReadWaveFunctionPlainBinary(&Psi, PsiFileName)){printf("Fatal Error in main\n Failed to read Psi from file %s !\n",PsiFileName);return 1;};
        printf("\nRead the starting wave function from the file %s\n",PsiFileName);
        sscanf(Psi.Name,"%s %i %i %i %i",Buffer,&PsiNBathVA,&PsiNBathCA,&PsiNBathVB,&PsiNBathCB);
        printf("Wave-function \"%s\" found has %i valence bath sites and %i conduction bath sites for A, %i valence bath sites and %i conduction bath sites for B, and %i Fermions in basis\n",Psi.Name,PsiNBathVA,PsiNBathCA,PsiNBathVB,PsiNBathCB,Psi.NFermion);
        if(2*(PsiNBathVA+PsiNBathCA+PsiNBathVB+PsiNBathCB+4)!=Psi.NFermion)
        {
            printf("The total number of Fermions is incompatible witht the bath sizes\n");
            printf("This wave-function can not be used to start from\n");
            return 1;
        }
    }
    else
    {
        printf("\nWave-function file: '%s' does not exist\nStarted from new function\n",PsiFileName);
        CreateWaveFunctionTwoOrbDMFTSecond(&Psi, EgDMFTSet.NBathA, EgDMFTSet.NBathA, EgDMFTSet.NBathB, EgDMFTSet.NBathB, EgDMFTSet.NBitsKey);
    }
    
    /* do the DMFT loops */

    for(iOuter=EgDMFTSet.iStart;iOuter<EgDMFTSet.iStart+EgDMFTSet.Ni;iOuter++){
        BeginTiming("DMFTLoop", &Time);
        printf("\n\n***********************************************************************\n");
        printf("Start of loop %i from %i to %i\n",iOuter+1,EgDMFTSet.iStart,EgDMFTSet.iStart+EgDMFTSet.Ni);
        printf("***********************************************************************\n\n\n");
        BeginTiming("Adjust mu", &Time);
        
RestartLanczos:
        
        if (iCount==0){
            mu=EgDMFTSet.mu1;
        }else if (iCount==1){
            mu=EgDMFTSet.mu2;
        }
        
        /*           */
        /* Orbital A */
        /*           */
        
        /* Create splitted anderson matrix */
        TriDiagonalMatrixToSecondDoubleTriDiagonalMatrix(GBathA, &T, &GAndA, Ef, epsilonAnd);
        /* Create GAnd and GAndTrunk as a copy pointing to the same memory */
        GAndATrunk.Ef=GAndA.Ef;
        GAndATrunk.Epsilon=GAndA.Epsilon;
        GAndATrunk.a0=GAndA.a0;
        GAndATrunk.b0=GAndA.b0;
        GAndATrunk.Dt.Ef=GAndA.Dt.Ef;
        GAndATrunk.Dt.Epsilon=GAndA.Dt.Epsilon;
        GAndATrunk.Dt.a0=GAndA.Dt.a0;
        GAndATrunk.Dt.b0=GAndA.Dt.b0;
        GAndATrunk.Dt.val.N=GAndA.Dt.val.N;
        GAndATrunk.Dt.val.a=GAndA.Dt.val.a;
        GAndATrunk.Dt.val.b=GAndA.Dt.val.b;
        GAndATrunk.Dt.con.N=GAndA.Dt.con.N;
        GAndATrunk.Dt.con.a=GAndA.Dt.con.a;
        GAndATrunk.Dt.con.b=GAndA.Dt.con.b;
        
        FreeMatrix(&T);
        if(GAndA.Dt.con.N<EgDMFTSet.NBathA){
            printf("Warning: GAndA.Dt.con.N(%i)<EgDMFTSet.NBathA(%i)\n",GAndA.Dt.con.N,EgDMFTSet.NBathA);
        }
        else{
            GAndATrunk.Dt.con.N=EgDMFTSet.NBathA;
        }
        if(GAndA.Dt.val.N<EgDMFTSet.NBathA){
            printf("Warning: GAndA.Dt.val.N(%i)<EgDMFTSet.NBathA(%i)\n",GAndA.Dt.val.N,EgDMFTSet.NBathA);
        }
        else{
            GAndATrunk.Dt.val.N=EgDMFTSet.NBathA;
        }
        
        /*           */
        /* Orbital B */
        /*           */
        
        /* Create splitted anderson matrix */
        TriDiagonalMatrixToSecondDoubleTriDiagonalMatrix(GBathB, &T, &GAndB, Ef, epsilonAnd);
        /* Create GAnd and GAndTrunk as a copy pointing to the same memory */
        GAndBTrunk.Ef=GAndB.Ef;
        GAndBTrunk.Epsilon=GAndB.Epsilon;
        GAndBTrunk.a0=GAndB.a0;
        GAndBTrunk.b0=GAndB.b0;
        GAndBTrunk.Dt.Ef=GAndB.Dt.Ef;
        GAndBTrunk.Dt.Epsilon=GAndB.Dt.Epsilon;
        GAndBTrunk.Dt.a0=GAndB.Dt.a0;
        GAndBTrunk.Dt.b0=GAndB.Dt.b0;
        GAndBTrunk.Dt.val.N=GAndB.Dt.val.N;
        GAndBTrunk.Dt.val.a=GAndB.Dt.val.a;
        GAndBTrunk.Dt.val.b=GAndB.Dt.val.b;
        GAndBTrunk.Dt.con.N=GAndB.Dt.con.N;
        GAndBTrunk.Dt.con.a=GAndB.Dt.con.a;
        GAndBTrunk.Dt.con.b=GAndB.Dt.con.b;
        
        FreeMatrix(&T);
        if(GAndB.Dt.con.N<EgDMFTSet.NBathB){
            printf("Warning: GAndB.Dt.con.N(%i)<EgDMFTSet.NBathB(%i)\n",GAndB.Dt.con.N,EgDMFTSet.NBathB);
        }
        else{
            GAndBTrunk.Dt.con.N=EgDMFTSet.NBathB;
        }
        if(GAndB.Dt.val.N<EgDMFTSet.NBathB){
            printf("Warning: GAndB.Dt.val.N(%i)<EgDMFTSet.NBathB(%i)\n",GAndB.Dt.val.N,EgDMFTSet.NBathB);
        }
        else{
            GAndBTrunk.Dt.val.N=EgDMFTSet.NBathB;
        }
        
        /* Create HAnd */
        
        TwoOrbSecondDoubleTriDiagonalMatrixToAndersonImpurityOperator(&HAnd, GAndATrunk, GAndATrunk, GAndBTrunk, GAndBTrunk, EgDMFTSet.U, EgDMFTSet.t1, EgDMFTSet.J, mu);
        
        /* Dyson like perturbation */
        
        if(EgDMFTSet.Dyson)
        {
            /* Create VPerturbation */
            for(i=0;i<4;i++)
            {
                InitOperator(&(VPerturbationA[i]));
                VPerturbationA[i].NFermion=4+2*(Min(GAndA.Dt.val.N, EgDMFTSet.NBathA)+Min(GAndA.Dt.con.N, EgDMFTSet.NBathA));
                VPerturbationA[i].NBoson=0;
                VPerturbationA[i].NWord=Min(VPerturbationA[i].NFermion, 1+(VPerturbationA[i].NFermion-1)/8);
                VPerturbationA[i].QComplex=0;
                InitOperatorFixedLength(&(VPerturbationA[i]),1,1);
            }
            
            for(i=0;i<4;i++)
            {
                InitOperator(&(VPerturbationB[i]));
                VPerturbationB[i].NFermion=4+2*(Min(GAndB.Dt.val.N, EgDMFTSet.NBathB)+Min(GAndB.Dt.con.N, EgDMFTSet.NBathB));
                VPerturbationB[i].NBoson=0;
                VPerturbationB[i].NWord=Min(VPerturbationB[i].NFermion, 1+(VPerturbationB[i].NFermion-1)/8);
                VPerturbationB[i].QComplex=0;
                InitOperatorFixedLength(&(VPerturbationB[i]),1,1);
            }
            
            /* V for Bath A val spin up and dn */
            
            if(GAndA.Dt.val.N>EgDMFTSet.NBathA)
            {
                oadd[0]=CrAnMask+7+EgDMFTSet.NBathA;
                RealOperatorAddLader(&(VPerturbationA[0]), oadd, 1, GAndA.Dt.val.b[EgDMFTSet.NBathA]);
                oadd[0]=CrAnMask+7+Min(GAndA.Dt.con.N, EgDMFTSet.NBathA)+Min(GAndB.Dt.val.N, EgDMFTSet.NBathB)+Min(GAndB.Dt.con.N, EgDMFTSet.NBathB)+2*EgDMFTSet.NBathA;
                RealOperatorAddLader(&(VPerturbationA[1]), oadd, 1, GAndA.Dt.val.b[EgDMFTSet.NBathA]);
            }
            
            /* V for Bath A con spin up and dn */
            
            if(GAndA.Dt.con.N>EgDMFTSet.NBathA)
            {
                oadd[0]=7+Min(GAndA.Dt.val.N, EgDMFTSet.NBathA)+EgDMFTSet.NBathA;
                RealOperatorAddLader(&(VPerturbationA[2]), oadd, 1, GAndA.Dt.con.b[EgDMFTSet.NBathA]);
                oadd[0]=7+2*Min(GAndA.Dt.val.N, EgDMFTSet.NBathA)+Min(GAndB.Dt.val.N, EgDMFTSet.NBathB)+Min(GAndB.Dt.con.N, EgDMFTSet.NBathB)+2*EgDMFTSet.NBathA;
                RealOperatorAddLader(&(VPerturbationA[3]), oadd, 1, GAndA.Dt.con.b[EgDMFTSet.NBathA]);
            }
            
            /* V for Bath B val spin up and dn */
            
            if(GAndB.Dt.val.N>EgDMFTSet.NBathB)
            {
                oadd[0]=CrAnMask+7+Min(GAndA.Dt.val.N, EgDMFTSet.NBathA)+Min(GAndA.Dt.con.N, EgDMFTSet.NBathA)+EgDMFTSet.NBathB;
                RealOperatorAddLader(&(VPerturbationB[0]), oadd, 1, GAndB.Dt.val.b[EgDMFTSet.NBathB]);
                oadd[0]=CrAnMask+7+2*Min(GAndA.Dt.val.N, EgDMFTSet.NBathA)+2*Min(GAndA.Dt.con.N, EgDMFTSet.NBathA)+Min(GAndB.Dt.con.N, EgDMFTSet.NBathB)+2*EgDMFTSet.NBathB;
                RealOperatorAddLader(&(VPerturbationB[1]), oadd, 1, GAndB.Dt.val.b[EgDMFTSet.NBathB]);
            }
            
            /* V for Bath B con spin up and dn */
            
            if(GAndB.Dt.con.N>EgDMFTSet.NBathB)
            {
                oadd[0]=7+Min(GAndA.Dt.val.N, EgDMFTSet.NBathA)+Min(GAndA.Dt.con.N, EgDMFTSet.NBathA)+Min(GAndB.Dt.val.N, EgDMFTSet.NBathB)+EgDMFTSet.NBathB;
                RealOperatorAddLader(&(VPerturbationB[2]), oadd, 1, GAndB.Dt.con.b[EgDMFTSet.NBathB]);
                oadd[0]=7+2*Min(GAndA.Dt.val.N, EgDMFTSet.NBathA)+2*Min(GAndA.Dt.con.N, EgDMFTSet.NBathA)+2*Min(GAndB.Dt.val.N, EgDMFTSet.NBathB)+2*EgDMFTSet.NBathB;
                RealOperatorAddLader(&(VPerturbationB[3]), oadd, 1, GAndB.Dt.con.b[EgDMFTSet.NBathB]);
            }
            
            /* Create gPerturbation for A */
            
            if(GAndA.Dt.val.N>EgDMFTSet.NBathA)
            {
                gPerturbationA[0].N=GAndA.Dt.val.N-EgDMFTSet.NBathA;
                gPerturbationA[1].N=GAndA.Dt.val.N-EgDMFTSet.NBathA;
                InitTriDiagonalMatrix(&(gPerturbationA[0]));
                InitTriDiagonalMatrix(&(gPerturbationA[1]));
                for(i=1;i<GAndA.Dt.val.N-EgDMFTSet.NBathA;i++)
                {
                    gPerturbationA[0].a[i]=-GAndA.Dt.val.a[i+EgDMFTSet.NBathA];
                    gPerturbationA[0].b[i]= GAndA.Dt.val.b[i+EgDMFTSet.NBathA];
                    gPerturbationA[1].a[i]=-GAndA.Dt.val.a[i+EgDMFTSet.NBathA];
                    gPerturbationA[1].b[i]= GAndA.Dt.val.b[i+EgDMFTSet.NBathA];
                }
                gPerturbationA[0].a[GAndA.Dt.val.N-EgDMFTSet.NBathA]=-GAndA.Dt.val.a[GAndA.Dt.val.N];
                gPerturbationA[1].a[GAndA.Dt.val.N-EgDMFTSet.NBathA]=-GAndA.Dt.val.a[GAndA.Dt.val.N];
            }
            else
            {
                gPerturbationA[0].N=1;
                gPerturbationA[1].N=1;
                InitTriDiagonalMatrix(&(gPerturbationA[0]));
                InitTriDiagonalMatrix(&(gPerturbationA[1]));
                gPerturbationA[0].a[1]=0;
                gPerturbationA[1].a[1]=0;
            }
            
            if(GAndA.Dt.con.N>EgDMFTSet.NBathA)
            {
                gPerturbationA[2].N=GAndA.Dt.con.N-EgDMFTSet.NBathA;
                gPerturbationA[3].N=GAndA.Dt.con.N-EgDMFTSet.NBathA;
                InitTriDiagonalMatrix(&(gPerturbationA[2]));
                InitTriDiagonalMatrix(&(gPerturbationA[3]));
                for(i=1;i<GAndA.Dt.con.N-EgDMFTSet.NBathA;i++)
                {
                    gPerturbationA[2].a[i]=GAndA.Dt.con.a[i+EgDMFTSet.NBathA];
                    gPerturbationA[2].b[i]=GAndA.Dt.con.b[i+EgDMFTSet.NBathA];
                    gPerturbationA[3].a[i]=GAndA.Dt.con.a[i+EgDMFTSet.NBathA];
                    gPerturbationA[3].b[i]=GAndA.Dt.con.b[i+EgDMFTSet.NBathA];
                }
                gPerturbationA[2].a[GAndA.Dt.con.N-EgDMFTSet.NBathA]=GAndA.Dt.con.a[GAndA.Dt.con.N];
                gPerturbationA[3].a[GAndA.Dt.con.N-EgDMFTSet.NBathA]=GAndA.Dt.con.a[GAndA.Dt.con.N];
            }
            else
            {
                gPerturbationA[2].N=1;
                gPerturbationA[3].N=1;
                InitTriDiagonalMatrix(&(gPerturbationA[2]));
                InitTriDiagonalMatrix(&(gPerturbationA[3]));
                gPerturbationA[2].a[1]=0;
                gPerturbationA[3].a[1]=0;
            }
            
            /* Create gPerturbation for B */
            
            if(GAndB.Dt.val.N>EgDMFTSet.NBathB)
            {
                gPerturbationB[0].N=GAndB.Dt.val.N-EgDMFTSet.NBathB;
                gPerturbationB[1].N=GAndB.Dt.val.N-EgDMFTSet.NBathB;
                InitTriDiagonalMatrix(&(gPerturbationB[0]));
                InitTriDiagonalMatrix(&(gPerturbationB[1]));
                for(i=1;i<GAndB.Dt.val.N-EgDMFTSet.NBathB;i++)
                {
                    gPerturbationB[0].a[i]=-GAndB.Dt.val.a[i+EgDMFTSet.NBathB];
                    gPerturbationB[0].b[i]= GAndB.Dt.val.b[i+EgDMFTSet.NBathB];
                    gPerturbationB[1].a[i]=-GAndB.Dt.val.a[i+EgDMFTSet.NBathB];
                    gPerturbationB[1].b[i]= GAndB.Dt.val.b[i+EgDMFTSet.NBathB];
                }
                gPerturbationB[0].a[GAndB.Dt.val.N-EgDMFTSet.NBathB]=-GAndB.Dt.val.a[GAndB.Dt.val.N];
                gPerturbationB[1].a[GAndB.Dt.val.N-EgDMFTSet.NBathB]=-GAndB.Dt.val.a[GAndB.Dt.val.N];
            }
            else
            {
                gPerturbationB[0].N=1;
                gPerturbationB[1].N=1;
                InitTriDiagonalMatrix(&(gPerturbationB[0]));
                InitTriDiagonalMatrix(&(gPerturbationB[1]));
                gPerturbationB[0].a[1]=0;
                gPerturbationB[1].a[1]=0;
            }
            
            if(GAndB.Dt.con.N>EgDMFTSet.NBathB)
            {
                gPerturbationB[2].N=GAndB.Dt.con.N-EgDMFTSet.NBathB;
                gPerturbationB[3].N=GAndB.Dt.con.N-EgDMFTSet.NBathB;
                InitTriDiagonalMatrix(&(gPerturbationB[2]));
                InitTriDiagonalMatrix(&(gPerturbationB[3]));
                for(i=1;i<GAndB.Dt.con.N-EgDMFTSet.NBathB;i++)
                {
                    gPerturbationB[2].a[i]=GAndB.Dt.con.a[i+EgDMFTSet.NBathB];
                    gPerturbationB[2].b[i]=GAndB.Dt.con.b[i+EgDMFTSet.NBathB];
                    gPerturbationB[3].a[i]=GAndB.Dt.con.a[i+EgDMFTSet.NBathB];
                    gPerturbationB[3].b[i]=GAndB.Dt.con.b[i+EgDMFTSet.NBathB];
                }
                gPerturbationB[2].a[GAndB.Dt.con.N-EgDMFTSet.NBathB]=GAndB.Dt.con.a[GAndB.Dt.con.N];
                gPerturbationB[3].a[GAndB.Dt.con.N-EgDMFTSet.NBathB]=GAndB.Dt.con.a[GAndB.Dt.con.N];
            }
            else
            {
                gPerturbationB[2].N=1;
                gPerturbationB[3].N=1;
                InitTriDiagonalMatrix(&(gPerturbationB[2]));
                InitTriDiagonalMatrix(&(gPerturbationB[3]));
                gPerturbationB[2].a[1]=0;
                gPerturbationB[3].a[1]=0;
            }
        }
        
        /* check the dimension of psi */
        
        sscanf(Psi.Name,"%s %i %i %i %i",Buffer,&PsiNBathVA,&PsiNBathCA, &PsiNBathVB, &PsiNBathCB);
        
        /* comparing conduction bath for B */
        if(PsiNBathCB!=GAndBTrunk.Dt.con.N){
            printf("\nBefore diagonalizing:\nThe wavefunction has %i conduction bath sites for B\nThe Anderson Hamiltonian %i\n",PsiNBathCB,GAndBTrunk.Dt.con.N);
            if(PsiNBathCB>GAndBTrunk.Dt.con.N){
                printf("Removing %i unoccupied fermions from the wave-function\n",PsiNBathCB-GAndBTrunk.Dt.con.N);
                for(i=0;i<PsiNBathCB-GAndBTrunk.Dt.con.N;i++){
                    WaveFunctionRemoveFermion(Psi, &PsiTemp, Psi.NFermion-1);
                    WaveFunctionFree(&Psi);
                    WaveFunctionRemoveFermion(PsiTemp, &Psi, 8+PsiNBathVA+PsiNBathCA+PsiNBathVB+PsiNBathCB-i-1);
                    WaveFunctionFree(&PsiTemp);
                }
            }
            else{
                printf("Adding %i unoccupied fermions to the wave-function\n",GAndBTrunk.Dt.con.N-PsiNBathCB);
                for(i=0;i<GAndBTrunk.Dt.con.N-PsiNBathCB;i++){
                    WaveFunctionAddFermion(Psi, &PsiTemp, Psi.NFermion,0);
                    WaveFunctionFree(&Psi);
                    WaveFunctionAddFermion(PsiTemp, &Psi, 8+PsiNBathVA+PsiNBathCA+PsiNBathVB+PsiNBathCB+i,0);
                    WaveFunctionFree(&PsiTemp);
                }
            }
        }
        
        /* comparing valence bath for B */
        if(PsiNBathVB!=GAndBTrunk.Dt.val.N){
            printf("\nBefore diagonalizing:\nThe wavefunction has %i valence bath sites for B\nThe Anderson Hamiltonian %i\n",PsiNBathVB,GAndBTrunk.Dt.val.N);
            if(PsiNBathVB>GAndBTrunk.Dt.val.N){
                printf("Removing %i unoccupied fermions from the wave-function\n",PsiNBathVB-GAndBTrunk.Dt.val.N);
                for(i=0;i<PsiNBathVB-GAndBTrunk.Dt.val.N;i++){
                    WaveFunctionRemoveFermion(Psi, &PsiTemp, Psi.NFermion-GAndBTrunk.Dt.con.N-1);
                    WaveFunctionFree(&Psi);
                    WaveFunctionRemoveFermion(PsiTemp, &Psi, 8+PsiNBathVA+PsiNBathCA+PsiNBathVB-i-1);
                    WaveFunctionFree(&PsiTemp);
                }
            }
            else{
                printf("Adding %i unoccupied fermions to the wave-function\n",GAndBTrunk.Dt.con.N-PsiNBathCB);
                for(i=0;i<GAndBTrunk.Dt.con.N-PsiNBathVB;i++){
                    WaveFunctionAddFermion(Psi, &PsiTemp, Psi.NFermion-GAndBTrunk.Dt.con.N,1);
                    WaveFunctionFree(&Psi);
                    WaveFunctionAddFermion(PsiTemp, &Psi, 8+PsiNBathVA+PsiNBathCA+PsiNBathVB+i,1);
                    WaveFunctionFree(&PsiTemp);
                }
            }
        }
        
        /* comparing conduction bath for A */
        if(PsiNBathCA!=GAndATrunk.Dt.con.N){
            printf("\nBefore diagonalizing:\nThe wavefunction has %i conduction bath sites for A\nThe Anderson Hamiltonian %i\n",PsiNBathCA,GAndATrunk.Dt.con.N);
            if(PsiNBathCA>GAndATrunk.Dt.con.N){
                printf("Removing %i unoccupied fermions from the wave-function\n",PsiNBathCA-GAndATrunk.Dt.con.N);
                for(i=0;i<PsiNBathCA-GAndATrunk.Dt.con.N;i++){
                    WaveFunctionRemoveFermion(Psi, &PsiTemp, Psi.NFermion-GAndBTrunk.Dt.con.N-GAndBTrunk.Dt.val.N-1);
                    WaveFunctionFree(&Psi);
                    WaveFunctionRemoveFermion(PsiTemp, &Psi, 8+PsiNBathVA+PsiNBathCA-i-1);
                    WaveFunctionFree(&PsiTemp);
                }
            }
            else{
                printf("Adding %i unoccupied fermions to the wave-function\n",GAndATrunk.Dt.con.N-PsiNBathCA);
                for(i=0;i<GAndATrunk.Dt.con.N-PsiNBathCA;i++){
                    WaveFunctionAddFermion(Psi, &PsiTemp, Psi.NFermion-GAndBTrunk.Dt.con.N-GAndBTrunk.Dt.val.N,1);
                    WaveFunctionFree(&Psi);
                    WaveFunctionAddFermion(PsiTemp, &Psi, 8+PsiNBathVA+PsiNBathCA+i,1);
                    WaveFunctionFree(&PsiTemp);
                }
            }
        }
        
        /* comparing valence bath for A */
        if(PsiNBathVA!=GAndATrunk.Dt.val.N){
            printf("\nBefore diagonalizing:\nThe wavefunction has %i valence bath sites for A\nThe Anderson Hamiltonian %i\n",PsiNBathVA,GAndATrunk.Dt.val.N);
            if(PsiNBathVA>GAndATrunk.Dt.val.N){
                printf("Removing %i occupied fermions from the wave-function\n",PsiNBathVA-GAndATrunk.Dt.val.N);
                for(i=0;i<PsiNBathVA-GAndATrunk.Dt.val.N;i++){
                    WaveFunctionRemoveFermion(Psi, &PsiTemp, Psi.NFermion-GAndBTrunk.Dt.con.N-GAndBTrunk.Dt.val.N-GAndA.Dt.con.N-1);
                    WaveFunctionFree(&Psi);
                    WaveFunctionRemoveFermion(PsiTemp, &Psi, 8+PsiNBathVA-i-1);
                    WaveFunctionFree(&PsiTemp);
                }
            }
            else{
                printf("Adding %i occupied fermions to the wave-function\n",GAndATrunk.Dt.val.N-PsiNBathVA);
                for(i=0;i<GAndATrunk.Dt.val.N-PsiNBathVA;i++){
                    WaveFunctionAddFermion(Psi, &PsiTemp, Psi.NFermion-GAndBTrunk.Dt.con.N-GAndBTrunk.Dt.val.N-GAndA.Dt.con.N,1);
                    WaveFunctionFree(&Psi);
                    WaveFunctionAddFermion(PsiTemp, &Psi, 8+PsiNBathVA+i,1);
                    WaveFunctionFree(&PsiTemp);
                }
            }
            printf("Corrected basis\n");
            PrintWaveFunctionReadableTruncated(Psi, 0.1);
        }
        sprintf(Psi.Name,"DMFT %i %i %i %i",GAndATrunk.Dt.val.N,GAndATrunk.Dt.con.N,GAndBTrunk.Dt.val.N,GAndBTrunk.Dt.con.N);
        
        /* Ground State properties */
        
        BeginTiming("Lanczos Groundstate", &Time);
        LanczosGroundState(HAnd,&Psi,&EgDMFTSet.zeroLancWavefunction,EgDMFTSet.zeroLancEnergy,EgDMFTSet.U,&Psi,0);
        EndTiming("Lanczos Groundstate", &Time);
        
        PsiHermitianOperatorPsi(HAnd, Psi, &PsiEnergy);
        PsiHermitianOperatorSqrPsi(HAnd, Psi, &PsiQuality);
        PsiQuality/=Sqr(PsiEnergy);
        PrintWaveFunctionReadableTruncated(Psi, 0.1);
        CalculatePrefactorHistogram(PsiHist, Psi, 10.0, 20);
        printf("\nHistogram of prefactors in Wave-function %i \n",i);
        for(i=0;i<19;i++)
        {
            printf("10^(-%2i) to 10^(-%2i) | %22.15E\n",i,i+1,PsiHist[i]);
        }
        printf("10^(-%2i) to 0        | %22.15E\n",i,PsiHist[19]);
        printf("\n");
        printf("Energy                  <psi|H^2|psi>/(<psi|H|psi>^2)\n");
        printf("%22.15E  %22.15E\n",PsiEnergy,PsiQuality);
        PsiRho=calloc(Psi.NFermion,sizeof(double));
        if(PsiRho==NULL){printf("calloc failed in main for PsiRho\n Fatal Error\n");fflush(stdout);return 1;}
        CalculateDensity(PsiRho,Psi);
        printf("\nOccupation \n");
        for(i=0;i<Psi.NFermion;i++)
        {
            printf("%5i %22.15E\n",i,PsiRho[i]);
        }
        Occu=PsiRho[0]+PsiRho[1]+PsiRho[2]+PsiRho[3];
        free(PsiRho);
        
        if(fabs(Occu-Goccu)<=ec){
            printf("\nOccupation %22.15e is correct.\n", Occu);
            printf("It took %i loops for convergence.\nThe current chemical potential is: %22.15e\n", iCount, mu);
            printf("\nStart creating spectra.\n");fflush(stdout);
            EndTiming("Adjust mu", &Time);
            iCount=0;
        }
        else if (iCount==0){
            muOld=mu;
            OccuOld=Occu;
            
            printf("\nFilling should be %10.5e.\nThe current occupation is: %22.15e.\n", Goccu, Occu);
            printf("mu adjusted from %22.15e (1st initial value) to %22.15e (2nd initial value)\n", EgDMFTSet.mu1, EgDMFTSet.mu2);
            printf("Restarting Lanczos Loop %i\n", ++iCount);fflush(stdout);
            
            goto RestartLanczos;
        }
        else {
            mutmp=mu;
            mu=(mu*(OccuOld-1.0)-muOld*(Occu-1.0))/(OccuOld-Occu);
            OccuOld=Occu;
            muOld=mutmp;
            printf("\nFilling should be %10.5e.\nThe current occupation is: %22.15e.\n", Goccu, Occu);
            printf("mu adjusted from %22.15e to %22.15e\n", muOld, mu);
            printf("Restarting Lanczos Loop %i\n", ++iCount);fflush(stdout);
            
            goto RestartLanczos;
        }

        /*                   */
        /* Calculate spectra */
        
        /* create PES and IPES operator */
        CreatePESOperatorDMFT(&OPESAup, HAnd.NFermion, 0);
        CreatePESOperatorDMFT(&OPESAdn, HAnd.NFermion, 1);
        CreateIPESOperatorDMFT(&OIPESAup, HAnd.NFermion, 0);
        CreateIPESOperatorDMFT(&OIPESAdn, HAnd.NFermion, 1);
        
        CreatePESOperatorDMFT(&OPESBup, HAnd.NFermion, 2);
        CreatePESOperatorDMFT(&OPESBdn, HAnd.NFermion, 3);
        CreateIPESOperatorDMFT(&OIPESBup, HAnd.NFermion, 2);
        CreateIPESOperatorDMFT(&OIPESBdn, HAnd.NFermion, 3);
            
        /* PES for spin up*/
        
        BeginTiming("Lanczos Tri-diagonalization", &Time);
        OperatorPsi(OPESAup, Psi, &OPsi);
        switch (EgDMFTSet.Dyson) {
            case 1:
                GclAMup.N=Max(EgDMFTSet.NTriDiag, EgDMFTSet.NGCluster);
                InitTriDiagonalMatrix(&GclAMup);
                LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, EgDMFTSet.zeroLancWavefunction, EgDMFTSet.NTriDiag, VPerturbationA, gPerturbationA, 4, &GclAMup);
                break;
                
            default:
                GclAMup.N=EgDMFTSet.NTriDiag;
                InitTriDiagonalMatrix(&GclAMup);
                LanczosTriDiagonalize(HAnd, OPsi, EgDMFTSet.zeroLancWavefunction, &GclAMup, EgDMFTSet.Gamma);
                break;
        }
        EndTiming("Lanczos Tri-diagonalization", &Time);
        
        BeginTiming("Lanczos Tri-diagonalization", &Time);
        OperatorPsi(OPESBup, Psi, &OPsi);
        switch (EgDMFTSet.Dyson) {
            case 1:
                GclBMup.N=Max(EgDMFTSet.NTriDiag, EgDMFTSet.NGCluster);
                InitTriDiagonalMatrix(&GclBMup);
                LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, EgDMFTSet.zeroLancWavefunction, EgDMFTSet.NTriDiag, VPerturbationB, gPerturbationB, 4, &GclBMup);
                break;
                
            default:
                GclBMup.N=EgDMFTSet.NTriDiag;
                InitTriDiagonalMatrix(&GclBMup);
                LanczosTriDiagonalize(HAnd, OPsi, EgDMFTSet.zeroLancWavefunction, &GclBMup, EgDMFTSet.Gamma);
                break;
        }
        EndTiming("Lanczos Tri-diagonalization", &Time);
        
        /* PES for spin dn*/
        
        BeginTiming("Lanczos Tri-diagonalization", &Time);
        OperatorPsi(OPESAdn, Psi, &OPsi);
        switch (EgDMFTSet.Dyson) {
            case 1:
                GclAMdn.N=Max(EgDMFTSet.NTriDiag, EgDMFTSet.NGCluster);
                InitTriDiagonalMatrix(&GclAMdn);
                LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, EgDMFTSet.zeroLancWavefunction, EgDMFTSet.NTriDiag, VPerturbationA, gPerturbationA, 4, &GclAMdn);
                break;
                
            default:
                GclAMdn.N=EgDMFTSet.NTriDiag;
                InitTriDiagonalMatrix(&GclAMdn);
                LanczosTriDiagonalize(HAnd, OPsi, EgDMFTSet.zeroLancWavefunction, &GclAMdn, EgDMFTSet.Gamma);
                break;
        }
        EndTiming("Lanczos Tri-diagonalization", &Time);
        
        BeginTiming("Lanczos Tri-diagonalization", &Time);
        OperatorPsi(OPESBdn, Psi, &OPsi);
        switch (EgDMFTSet.Dyson) {
            case 1:
                GclBMdn.N=Max(EgDMFTSet.NTriDiag, EgDMFTSet.NGCluster);
                InitTriDiagonalMatrix(&GclBMdn);
                LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, EgDMFTSet.zeroLancWavefunction, EgDMFTSet.NTriDiag, VPerturbationB, gPerturbationB, 4, &GclBMdn);
                break;
                
            default:
                GclBMdn.N=EgDMFTSet.NTriDiag;
                InitTriDiagonalMatrix(&GclBMdn);
                LanczosTriDiagonalize(HAnd, OPsi, EgDMFTSet.zeroLancWavefunction, &GclBMdn, EgDMFTSet.Gamma);
                break;
        }
        EndTiming("Lanczos Tri-diagonalization", &Time);
        
        /* IPES for spin up*/
        
        BeginTiming("Lanczos Tri-diagonalization", &Time);
        OperatorPsi(OIPESAup, Psi, &OPsi);
        switch (EgDMFTSet.Dyson) {
            case 1:
                GclAPup.N=Max(EgDMFTSet.NTriDiag, EgDMFTSet.NGCluster);
                InitTriDiagonalMatrix(&GclAPup);
                LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, EgDMFTSet.zeroLancWavefunction, EgDMFTSet.NTriDiag, VPerturbationA, gPerturbationA, 4, &GclAPup);
                break;
                
            default:
                GclAPup.N=EgDMFTSet.NTriDiag;
                InitTriDiagonalMatrix(&GclAPup);
                LanczosTriDiagonalize(HAnd, OPsi, EgDMFTSet.zeroLancWavefunction, &GclAPup, EgDMFTSet.Gamma);
                break;
        }
        EndTiming("Lanczos Tri-diagonalization", &Time);
        
        BeginTiming("Lanczos Tri-diagonalization", &Time);
        OperatorPsi(OIPESBup, Psi, &OPsi);
        switch (EgDMFTSet.Dyson) {
            case 1:
                GclBPup.N=Max(EgDMFTSet.NTriDiag, EgDMFTSet.NGCluster);
                InitTriDiagonalMatrix(&GclBPup);
                LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, EgDMFTSet.zeroLancWavefunction, EgDMFTSet.NTriDiag, VPerturbationB, gPerturbationB, 4, &GclBPup);
                break;
                
            default:
                GclBPup.N=EgDMFTSet.NTriDiag;
                InitTriDiagonalMatrix(&GclBPup);
                LanczosTriDiagonalize(HAnd, OPsi, EgDMFTSet.zeroLancWavefunction, &GclBPup, EgDMFTSet.Gamma);
                break;
        }
        EndTiming("Lanczos Tri-diagonalization", &Time);
        
        /* IPES for spin dn*/
        
        BeginTiming("Lanczos Tri-diagonalization", &Time);
        OperatorPsi(OIPESAdn, Psi, &OPsi);
        switch (EgDMFTSet.Dyson) {
            case 1:
                GclAPdn.N=Max(EgDMFTSet.NTriDiag, EgDMFTSet.NGCluster);
                InitTriDiagonalMatrix(&GclAPdn);
                LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, EgDMFTSet.zeroLancWavefunction, EgDMFTSet.NTriDiag, VPerturbationA, gPerturbationA, 4, &GclAPdn);
                break;
                
            default:
                GclAPdn.N=EgDMFTSet.NTriDiag;
                InitTriDiagonalMatrix(&GclAPdn);
                LanczosTriDiagonalize(HAnd, OPsi, EgDMFTSet.zeroLancWavefunction, &GclAPdn, EgDMFTSet.Gamma);
                break;
        }
        EndTiming("Lanczos Tri-diagonalization", &Time);
        
        BeginTiming("Lanczos Tri-diagonalization", &Time);
        OperatorPsi(OIPESBdn, Psi, &OPsi);
        switch (EgDMFTSet.Dyson) {
            case 1:
                GclBPdn.N=Max(EgDMFTSet.NTriDiag, EgDMFTSet.NGCluster);
                InitTriDiagonalMatrix(&GclBPdn);
                LanczosTriDiagonalizeDysonLowMem(HAnd, OPsi, EgDMFTSet.zeroLancWavefunction, EgDMFTSet.NTriDiag, VPerturbationB, gPerturbationB, 4, &GclBPdn);
                break;
                
            default:
                GclBPdn.N=EgDMFTSet.NTriDiag;
                InitTriDiagonalMatrix(&GclBPdn);
                LanczosTriDiagonalize(HAnd, OPsi, EgDMFTSet.zeroLancWavefunction, &GclBPdn, EgDMFTSet.Gamma);
                break;
        }
        EndTiming("Lanczos Tri-diagonalization", &Time);
        
        /* calculate G */
        printf("=== Start G math ===\n");fflush(stdout);
        
        BeginTiming("G math", &Time);
        
        ShiftTriDiagonalMatrix(&GclAMup, -PsiEnergy);
        ShiftTriDiagonalMatrix(&GclAMdn, -PsiEnergy);
        ShiftTriDiagonalMatrix(&GclAPup, -PsiEnergy);
        ShiftTriDiagonalMatrix(&GclAPdn, -PsiEnergy);
        
        ShiftTriDiagonalMatrix(&GclBMup, -PsiEnergy);
        ShiftTriDiagonalMatrix(&GclBMdn, -PsiEnergy);
        ShiftTriDiagonalMatrix(&GclBPup, -PsiEnergy);
        ShiftTriDiagonalMatrix(&GclBPdn, -PsiEnergy);

        for(i=1;i<=GclAMup.N;i++)
        {
            GclAMup.a[i]*=-1;
        }
        
        for(i=1;i<=GclAMdn.N;i++)
        {
            GclAMdn.a[i]*=-1;
        }
        
        for(i=1;i<=GclBMup.N;i++)
        {
            GclBMup.a[i]*=-1;
        }
        
        for(i=1;i<=GclBMdn.N;i++)
        {
            GclBMdn.a[i]*=-1;
        }
        
        printf("Start SumGreensfunctionsTriDiagonalRepresentation for A");fflush(stdout);
        
        TriDiagonalMatrixToListOfPoles(GclAPup, &GclAList);
        TriDiagonalMatrixToListOfPoles(GclAPdn, &GListTemp);
        AddListOfPoles(&GclAList, GListTemp);
        FreeListOfPoles(&GListTemp);
        TriDiagonalMatrixToListOfPoles(GclAMup, &GListTemp);
        AddListOfPoles(&GclAList, GListTemp);
        FreeListOfPoles(&GListTemp);
        TriDiagonalMatrixToListOfPoles(GclAMdn, &GListTemp);
        AddListOfPoles(&GclAList, GListTemp);
        FreeListOfPoles(&GListTemp);
        ScaleListOfPoles(&GclAList,0.5);
        ListOfPolesToTridiagonalMatrix(GclAList, &GclA, GclAList.N);
        sprintf(Buffer, "GclA");
        strcpy(GclA.Name, Buffer);
        GclA.b[0]=1.0;
        
        printf("Start SumGreensfunctionsTriDiagonalRepresentation for B");fflush(stdout);
        
        TriDiagonalMatrixToListOfPoles(GclBPup, &GclBList);
        TriDiagonalMatrixToListOfPoles(GclBPdn, &GListTemp);
        AddListOfPoles(&GclBList, GListTemp);
        FreeListOfPoles(&GListTemp);
        TriDiagonalMatrixToListOfPoles(GclBMup, &GListTemp);
        AddListOfPoles(&GclBList, GListTemp);
        FreeListOfPoles(&GListTemp);
        TriDiagonalMatrixToListOfPoles(GclBMdn, &GListTemp);
        AddListOfPoles(&GclBList, GListTemp);
        FreeListOfPoles(&GListTemp);
        ScaleListOfPoles(&GclBList,0.5);
        ListOfPolesToTridiagonalMatrix(GclBList, &GclB, GclBList.N);
        sprintf(Buffer, "GclB");
        strcpy(GclB.Name, Buffer);
        GclB.b[0]=1.0;
        
        /* calculate Sigma */
        
        BeginTiming("G math Selfenergy", &Time);
        printf("Start SelfEnergyTriDiagonalRepresentation\n");fflush(stdout);
        
        TriDiagonalMatrixToAndersonMatrix(GclA, &T, &GclAAnd);
        FreeMatrix(&T);
        TriDiagonalMatrixToAndersonMatrix(GBathA, &T, &GBathAAnd);
        FreeMatrix(&T);
        SelfEnergyListOfPolesFromAndersonRepresentation(GBathAAnd, GclAAnd, &SigmaAPoles);
        printf("Self energy of A created of size %i\n",SigmaAPoles.N);
        strcpy(SigmaAPoles.Name,"Sigma A");
        ListOfPolesToTridiagonalMatrix(SigmaAPoles, &SigmaA, SigmaAPoles.N);
        
        TriDiagonalMatrixToAndersonMatrix(GclB, &T, &GclBAnd);
        FreeMatrix(&T);
        TriDiagonalMatrixToAndersonMatrix(GBathB, &T, &GBathBAnd);
        FreeMatrix(&T);
        SelfEnergyListOfPolesFromAndersonRepresentation(GBathBAnd, GclBAnd, &SigmaBPoles);
        printf("Self energy of B created of size %i\n",SigmaBPoles.N);
        strcpy(SigmaBPoles.Name,"Sigma B");
        ListOfPolesToTridiagonalMatrix(SigmaBPoles, &SigmaB, SigmaBPoles.N);
        
        EndTiming("G math Selfenergy", &Time);
        
        /* G bath */
        
        BeginTiming("G math GBath", &Time);
        
        printf("Start GBathAndListAndRepresentation\n");fflush(stdout);
        
        GBathAndTriAndRepresentation(G0AndA, SigmaA, &GBathANewAnd);
        printf("New bath Green function of A created of size %i ",GBathANewAnd.N);
        strcpy(GBathANewAnd.Name,"GbathA-New");
        ReducePolesAndersonMatrixNumber(&GBathANewAnd, EgDMFTSet.NGBath);
        printf("reduced to %i\n",GBathANewAnd.N);
        AndersonMatrixToTriDiagonalMatrix(GBathANewAnd, &T, &GBathANew, GBathANewAnd.N);
        FreeMatrix(&T);
        
        GBathAndTriAndRepresentation(G0AndB, SigmaB, &GBathBNewAnd);
        printf("New bath Green function of B created of size %i ",GBathBNewAnd.N);
        strcpy(GBathBNewAnd.Name,"GbathB-New");
        ReducePolesAndersonMatrixNumber(&GBathBNewAnd, EgDMFTSet.NGBath);
        printf("reduced to %i\n",GBathBNewAnd.N);
        AndersonMatrixToTriDiagonalMatrix(GBathBNewAnd, &T, &GBathBNew, GBathBNewAnd.N);
        FreeMatrix(&T);
        
        EndTiming("G math GBath", &Time);
        
        /* mixing */
        
        printf("Start mixing GBath with GBathNew of A.\n");
        TriDiagonalMatrixCopy(&GBathAOld, GBathA);
        FreeTriDiagonalMatrix(&GBathA);
        GBathAOld.N=Min(GBathAOld.N,EgDMFTSet.NGBath);
        strcpy(GBathAOld.Name,"GbathA-Old");
        GBathAOld.b[0]=sqrt(EgDMFTSet.Mix);
        GBathANew.b[0]=sqrt(1-EgDMFTSet.Mix);
        SumGreensfunctionsTriDiagonalRepresentation(GBathANew, GBathAOld, &GBathA);
        strcpy(GBathA.Name,"GbathA");
        GBathAOld.b[0]=1;
        GBathANew.b[0]=1;
        
        printf("Start mixing GBath with GBathNew of B.\n");
        TriDiagonalMatrixCopy(&GBathBOld, GBathB);
        FreeTriDiagonalMatrix(&GBathB);
        GBathBOld.N=Min(GBathBOld.N,EgDMFTSet.NGBath);
        strcpy(GBathBOld.Name,"GbathB-Old");
        GBathBOld.b[0]=sqrt(EgDMFTSet.Mix);
        GBathBNew.b[0]=sqrt(1-EgDMFTSet.Mix);
        SumGreensfunctionsTriDiagonalRepresentation(GBathBNew, GBathBOld, &GBathB);
        strcpy(GBathB.Name,"GbathB");
        GBathBOld.b[0]=1;
        GBathBNew.b[0]=1;
        
        EndTiming("G math", &Time);
        
        /* create spectra */

        Spectra[0].N       =EgDMFTSet.NE;
        Spectra[0].Emin    =EgDMFTSet.Emin;
        Spectra[0].Emax    =EgDMFTSet.Emax;
        Spectra[0].dE      =(EgDMFTSet.Emax-EgDMFTSet.Emin)/EgDMFTSet.NE;
        Spectra[0].EminPole=EgDMFTSet.Emin;
        Spectra[0].EmaxPole=EgDMFTSet.Emax;
        Spectra[0].Gamma   =EgDMFTSet.Gamma;
        
        for(i=0;i<16;i++)
        {
            SpectrumAllocate(Spectra[0], &Spectra[i]);
        }
        
        TriDiagonalMatrixToSpectrum(0,GBathA,&(Spectra[0]));               /*  3 */
        TriDiagonalMatrixToSpectrum(0,GBathB,&(Spectra[1]));               /*  5 */
        TriDiagonalMatrixToSpectrum(0,GBathAOld,&(Spectra[2]));            /*  7 */
        TriDiagonalMatrixToSpectrum(0,GBathBOld,&(Spectra[3]));            /*  9 */
        TriDiagonalMatrixToSpectrum(0,GBathANew,&(Spectra[4]));            /* 11 */
        TriDiagonalMatrixToSpectrum(0,GBathBNew,&(Spectra[5]));            /* 13 */
        SecondDoubleTriDiagonalMatrixToSpectrum(0,GAndA,&(Spectra[6]));    /* 15 */
        SecondDoubleTriDiagonalMatrixToSpectrum(0,GAndB,&(Spectra[7]));    /* 17 */
        TriDiagonalMatrixToSpectrum(0,SigmaA,&(Spectra[8]));               /* 19 */
        TriDiagonalMatrixToSpectrum(0,SigmaB,&(Spectra[9]));               /* 21 */
        AndersonImpurityMatrixToSpectrum(0,GBathANewAnd,&(Spectra[10]));   /* 23 */
        AndersonImpurityMatrixToSpectrum(0,GBathBNewAnd,&(Spectra[11]));   /* 25 */
        TriDiagonalMatrixToSpectrum(0,GclA,&(Spectra[12]));                /* 27 */
        TriDiagonalMatrixToSpectrum(0,GclB,&(Spectra[13]));                /* 29 */
        TriDiagonalMatrixToSpectrum(0,G0A,&(Spectra[14]));                 /* 31 */
        TriDiagonalMatrixToSpectrum(0,G0B,&(Spectra[15]));                 /* 33 */
        
        /* write to disk */
        printf("Start write to dist\n");fflush(stdout);
        if(iOuter==0)
        {
            sprintf(GBathFileNameA,"%s001",EgDMFTSet.GRootFileNameA);
            sprintf(GclFileNameA,"Cl%s001",EgDMFTSet.GRootFileNameA);
            sprintf(GSigmaFileNameA,"Sigma%s001",EgDMFTSet.GRootFileNameA);
            sprintf(GBathFileNameB,"%s001",EgDMFTSet.GRootFileNameB);
            sprintf(GclFileNameB,"Cl%s001",EgDMFTSet.GRootFileNameB);
            sprintf(GSigmaFileNameB,"Sigma%s001",EgDMFTSet.GRootFileNameB);
            sprintf(PsiFileName,"%s001",EgDMFTSet.WaveFunctionRootFileName);
            sprintf(SpecFileName,"%s001",EgDMFTSet.SpecRootFileName);
        }
        else if(iOuter<9)
        {
            sprintf(GBathFileNameA,"%s00%1i",EgDMFTSet.GRootFileNameA,iOuter+1);
            sprintf(GclFileNameA,"Cl%s00%1i",EgDMFTSet.GRootFileNameA,iOuter+1);
            sprintf(GSigmaFileNameA,"Sigma%s00%1i",EgDMFTSet.GRootFileNameA,iOuter+1);
            sprintf(GBathFileNameB,"%s00%1i",EgDMFTSet.GRootFileNameB,iOuter+1);
            sprintf(GclFileNameB,"Cl%s00%1i",EgDMFTSet.GRootFileNameB,iOuter+1);
            sprintf(GSigmaFileNameB,"Sigma%s00%1i",EgDMFTSet.GRootFileNameB,iOuter+1);
            sprintf(PsiFileName,"%s00%1i",EgDMFTSet.WaveFunctionRootFileName,iOuter+1);
            sprintf(SpecFileName,"%s00%1i",EgDMFTSet.SpecRootFileName,iOuter+1);
        }
        else if(iOuter<99)
        {
            sprintf(GBathFileNameA,"%s0%2i",EgDMFTSet.GRootFileNameA,iOuter+1);
            sprintf(GclFileNameA,"Cl%s0%2i",EgDMFTSet.GRootFileNameA,iOuter+1);
            sprintf(GSigmaFileNameA,"Sigma%s0%2i",EgDMFTSet.GRootFileNameA,iOuter+1);
            sprintf(GBathFileNameB,"%s0%2i",EgDMFTSet.GRootFileNameB,iOuter+1);
            sprintf(GclFileNameB,"Cl%s0%2i",EgDMFTSet.GRootFileNameB,iOuter+1);
            sprintf(GSigmaFileNameB,"Sigma%s0%2i",EgDMFTSet.GRootFileNameB,iOuter+1);
            sprintf(PsiFileName,"%s0%2i",EgDMFTSet.WaveFunctionRootFileName,iOuter+1);
            sprintf(SpecFileName,"%s0%2i",EgDMFTSet.SpecRootFileName,iOuter+1);
        }
        else if(iOuter<999)
        {
            sprintf(GBathFileNameA,"%s%3i",EgDMFTSet.GRootFileNameA,iOuter+1);
            sprintf(GclFileNameA,"Cl%s%3i",EgDMFTSet.GRootFileNameA,iOuter+1);
            sprintf(GSigmaFileNameA,"Sigma%s%3i",EgDMFTSet.GRootFileNameA,iOuter+1);
            sprintf(GBathFileNameB,"%s%3i",EgDMFTSet.GRootFileNameB,iOuter+1);
            sprintf(GclFileNameB,"Cl%s%3i",EgDMFTSet.GRootFileNameB,iOuter+1);
            sprintf(GSigmaFileNameB,"Sigma%s%3i",EgDMFTSet.GRootFileNameB,iOuter+1);
            sprintf(PsiFileName,"%s%3i",EgDMFTSet.WaveFunctionRootFileName,iOuter+1);
            sprintf(SpecFileName,"%s%3i",EgDMFTSet.SpecRootFileName,iOuter+1);
        }
        FilePrintSpectra(Spectra, 16, SpecFileName, "w");
        FilePrintWaveFunctionPlainBinary(Psi, PsiFileName, "w");
        FilePrintPlainTriDiagonalMatrix(GBathA,GBathFileNameA,"w");
        FilePrintPlainTriDiagonalMatrix(GBathB,GBathFileNameB,"w");
        FilePrintPlainTriDiagonalMatrix(GclA,GclFileNameA,"w");
        FilePrintPlainTriDiagonalMatrix(GclB,GclFileNameB,"w");
        FilePrintPlainTriDiagonalMatrix(SigmaA,GSigmaFileNameA,"w");
        FilePrintPlainTriDiagonalMatrix(SigmaB,GSigmaFileNameB,"w");
        
        /* free variables */
        printf("Start free variables\n");fflush(stdout);
        FreeTriDiagonalMatrix(&GclAMup);
        FreeTriDiagonalMatrix(&GclAMdn);
        FreeTriDiagonalMatrix(&GclBMup);
        FreeTriDiagonalMatrix(&GclBMdn);
        FreeTriDiagonalMatrix(&GclAPup);
        FreeTriDiagonalMatrix(&GclAPdn);
        FreeTriDiagonalMatrix(&GclBPup);
        FreeTriDiagonalMatrix(&GclBPdn);
        FreeTriDiagonalMatrix(&GclA);
        FreeTriDiagonalMatrix(&GclB);
        FreeListOfPoles(&SigmaAPoles);
        FreeListOfPoles(&SigmaBPoles);
        FreeTriDiagonalMatrix(&SigmaA);
        FreeTriDiagonalMatrix(&GBathANew);
        FreeTriDiagonalMatrix(&GBathAOld);
        FreeTriDiagonalMatrix(&SigmaB);
        FreeTriDiagonalMatrix(&GBathBNew);
        FreeTriDiagonalMatrix(&GBathBOld);
        FreeDoubleTriDiagonalMatrix(&(GAndA.Dt));
        FreeDoubleTriDiagonalMatrix(&(GAndB.Dt));
        FreeAndersonMatrix(&GBathANewAnd);
        FreeAndersonMatrix(&GBathBNewAnd);
        FreeAndersonMatrix(&GBathAAnd);
        FreeAndersonMatrix(&GBathBAnd);
        EndTiming("DMFTLoop", &Time);
        PrintTiming(Time);
        if (EgDMFTSet.Dyson) {
            for (i=0; i<4; i++) {
                FreeOperator(&(VPerturbationA[i]));
                FreeOperator(&(VPerturbationB[i]));
                FreeTriDiagonalMatrix(&(gPerturbationA[i]));
                FreeTriDiagonalMatrix(&(gPerturbationB[i]));
            }
        }
        FreeOperator(&HAnd);
        FreeOperator(&OPESAup);
        FreeOperator(&OPESBup);
        FreeOperator(&OPESAdn);
        FreeOperator(&OPESBdn);
        FreeOperator(&OIPESAup);
        FreeOperator(&OIPESBup);
        FreeOperator(&OIPESAdn);
        FreeOperator(&OIPESBdn);
    }
    EndTiming("Total", &Time);
    printf("\n\n*************************** Finished ***************************\n\n");
    PrintTiming(Time);
    
    return 0;

}