--Verbosity(0xFFFF)
accuracy = 1E-6
--Import the List to combine Raw spectra into cartesian RIXS tensor. this list has been generated with
--the script ComputeRixsTensor.nb in the parent directory.
dofile("./SumTensorE1E1KHarmonics.lua")
dofile("./SumTensorE1E1Cart.lua")


function TestIfSpectrumIsZero(Spec)
	local SpecTable=Spectra.ToTable(Spec)
	for i=1,#SpecTable do
		for j=1,#SpecTable[i] do
			if(not ( Chop(SpecTable[i][j][2])==0) ) then
				print("Test failed at i=",i,"j=",j,"val=",SpecTable[i][j][2])
				return false
			end
		end
	end
    print("Test passed: the spectrum is zero within the accuracy")
	return true
end


--returns the complex vector re + I im normalized to one
function NormalizedPolarizationVector(re, im)
    local norm2 = 0
    for i = 1, 3 do
        norm2 = norm2 + re[i] * re[i] + im[i] * im[i]
    end
    local norm = math.sqrt(norm2)
    return {
        Complex.New(re[1]/norm, im[1]/norm),
        Complex.New(re[2]/norm, im[2]/norm),
        Complex.New(re[3]/norm, im[3]/norm)
    }
end

--fixed, random looking polarization directions (real and complex)
PolarizationVectors = {
    NormalizedPolarizationVector({1, 0, 0}, {0, 0, 0}),
    NormalizedPolarizationVector({1, 2, 3}, {0, 0, 0}),
    NormalizedPolarizationVector({1, 1, 1}, {0, 0, 0}),
    NormalizedPolarizationVector({1, 0, 0}, {0, 1, 0}),
    NormalizedPolarizationVector({0.3, -1.1, 0.5}, {0.7, 0.2, -0.9})
}


---------------------------------------------------------------------------------------------
-------------------------------General Parameters--------------------------------------------

-- in crystal field theory U drops out of the equation
U       =  0.000 
F2dd    = 11.142 
F4dd    =  6.874
F0dd    = U+(F2dd+F4dd)*2/63
-- in crystal field theory U drops out of the equation
Upd     =  0.000 
F2pd    =  6.667
G1pd    =  4.922
G3pd    =  2.796
F0pd    =  Upd + G1pd*1/15 + G3pd*3/70
tenDq   =  1.100
zeta_3d =  0.081
zeta_2p = 11.498

Bz=0.000001
--formal valence 3d shell
n = 8


Neig = 45
if n==0 then
    Neig=1
end

------------------------------Parameters for the spectras-------------------------

--size of Krylov basis (std 200)
NTri1=200
NTri2=200

--Gamma forthe maps 
Gamma1=1--I enhanced this since se generate lot of data
Gamma2=0.2

--RIXS energies
dE1=1
E01_0=-0.2

dE2=4
E02_0=-2

--Non relativistic energies for the core states. Change them for every different An!
E2pCM=365.0
E3dCM=0 --the 3d is a valence state.

--shift depending on number of electrons and edge (here L3, so 2p3/2 hole)
E01=E01_0 + Upd*(5-n) + n*U -1/2*zeta_2p
E02=E02_0 

--the steps in the energy mesh must be small enough.
--Step_1=Gamma1/5--I set to 5 since we generate a lot of data
Step_1=Gamma1/10
Step_2=Gamma2/10

NE1=math.floor(dE1/Step_1)
NE2=math.floor(dE2/Step_2)
print("Map dimension:"..NE1,NE2)






---------------------------------Create the Hamiltonian-------------------------------------------
--define Basis
Index, NF = CreateAtomicIndicesDict({"3d","2p"})

---Create Operators
--Valence-Valence 3d interaction
OppF0_3d = NewOperator("U", NF, Index["3d_Up"], Index["3d_Dn"], {1, 0, 0})
OppF2_3d = NewOperator("U", NF, Index["3d_Up"], Index["3d_Dn"], {0, 1, 0})
OppF4_3d = NewOperator("U", NF, Index["3d_Up"], Index["3d_Dn"], {0, 0, 1})

-- core valence interactions 2p-3d
OppUdfF0 = NewOperator("U", NF, Index["2p_Up"], Index["2p_Dn"], Index["3d_Up"], Index["3d_Dn"], {1, 0}, {0, 0})
OppUdfF2 = NewOperator("U", NF, Index["2p_Up"], Index["2p_Dn"], Index["3d_Up"], Index["3d_Dn"], {0, 1}, {0, 0})

OppUdfG1 = NewOperator("U", NF, Index["2p_Up"], Index["2p_Dn"], Index["3d_Up"], Index["3d_Dn"], {0, 0}, {1, 0})
OppUdfG3 = NewOperator("U", NF, Index["2p_Up"], Index["2p_Dn"], Index["3d_Up"], Index["3d_Dn"], {0, 0}, {0, 1})


--Define the total coulomb operators for the 3 situations------------------
OpU3d = F0dd * OppF0_3d + F2dd * OppF2_3d + F4dd * OppF4_3d 
OpU2p3d = F0pd * OppUdfF0 + F2pd * OppUdfF2 + G1pd * OppUdfG1 + G3pd * OppUdfG3 

OpU_ini = OpU3d
OpU_XAS = OpU3d + OpU2p3d 
OpU_fin = OpU_ini


---------------------SOC-----------------------------------------------------
Oppldots_2p = NewOperator("ldots", NF, Index["2p_Up"],Index["2p_Dn"])
Oppldots_3d = NewOperator("ldots",NF, Index["3d_Up"],Index["3d_Dn"])


--------------------Operators for angular momentum--------------------------
Opp3dSz   = NewOperator("Sz"   ,NF, Index["3d_Up"],Index["3d_Dn"])
Opp3dSsqr = NewOperator("Ssqr" ,NF, Index["3d_Up"],Index["3d_Dn"])

Opp3dLz  = NewOperator("Lz"   ,NF, Index["3d_Up"],Index["3d_Dn"])
Opp3dLsqr = NewOperator("Lsqr" ,NF, Index["3d_Up"],Index["3d_Dn"])

Opp3dJz   = NewOperator("Jz"   ,NF, Index["3d_Up"],Index["3d_Dn"])
Opp3dJsqr = NewOperator("Jsqr" ,NF, Index["3d_Up"],Index["3d_Dn"])


---------------------Some counting operators-------------------
OppNUp_2p = NewOperator("Number", NF, Index["2p_Up"],Index["2p_Up"], {1,1,1})
OppNDn_2p = NewOperator("Number", NF, Index["2p_Dn"],Index["2p_Dn"], {1,1,1})
OppN_2p = OppNUp_2p + OppNDn_2p
OppNUp_3d = NewOperator("Number", NF, Index["3d_Up"],Index["3d_Up"], {1,1,1,1,1})
OppNDn_3d = NewOperator("Number", NF, Index["3d_Dn"],Index["3d_Dn"], {1,1,1,1,1})
OppN_3d = OppNUp_3d + OppNDn_3d


---------------------Crystal Field--------------------------------------------

--Oh
--set the parameters
Eeg  = 3/5*tenDq
Et2g = -2/5*tenDq

AkmOh = {{0, 0, (2/5)*(Eeg) + (3/5)*(Et2g)} , 
{4, 0, (21/10)*(Eeg + (-1)*(Et2g))} , 
{4,-4, (3/2)*((sqrt(7/10))*(Eeg + (-1)*(Et2g)))} , 
{4, 4, (3/2)*((sqrt(7/10))*(Eeg + (-1)*(Et2g)))} }

HCFOh = NewOperator("CF", NF, Index["3d_Up"],Index["3d_Dn"], AkmOh)
HCFOh = 1/2*(HCFOh+ConjugateTranspose(HCFOh))


--D4h
--set the parameters
deltaEegD4h  = 0.72
deltaEt2gD4h = -1.23
Eb1g = 3/5*tenDq-deltaEegD4h
Ea1g = 3/5*tenDq+deltaEegD4h
Eeg  = -2/5*tenDq-1/3*deltaEt2gD4h
Eb2g = -2/5*tenDq+2/3*deltaEt2gD4h

AkmD4h = {{0, 0, (1/5)*(Ea1g + Eb1g + Eb2g + (2)*(Eeg))} , 
{2, 0, Ea1g + (-1)*(Eb1g) + (-1)*(Eb2g) + Eeg} , 
{4, 0, (3/10)*((6)*(Ea1g) + Eb1g + Eb2g + (-8)*(Eeg))} , 
{4,-4, (3/2)*((sqrt(7/10))*(Eb1g + (-1)*(Eb2g)))} , 
{4, 4, (3/2)*((sqrt(7/10))*(Eb1g + (-1)*(Eb2g)))} }

HCFD4h = NewOperator("CF", NF, Index["3d_Up"],Index["3d_Dn"], AkmD4h)
HCFD4h = 1/2*(HCFD4h+ConjugateTranspose(HCFD4h))


--------------------SUM UP THE PARTS-------------------------------------

--SO3
HamiltonianSO3 = zeta_3d * Oppldots_3d + OpU3d + Bz * (2 * Opp3dSz + Opp3dLz)
XASHamiltonianSO3 = zeta_3d * Oppldots_3d + zeta_2p * Oppldots_2p + OpU3d + OpU2p3d 
finHamiltonianSO3 = HamiltonianSO3

--Oh
HamiltonianOh = HamiltonianSO3 + HCFOh
XASHamiltonianOh = XASHamiltonianSO3+HCFOh
finHamiltonianOh = finHamiltonianSO3+HCFOh

--HD4h
HamiltonianD4h = HamiltonianSO3 +HCFD4h
XASHamiltonianD4h = XASHamiltonianSO3+HCFD4h
finHamiltonianD4h = finHamiltonianSO3+HCFD4h


------------------------------------------X-ray---------------------------------------------------------
----------------------------------Create dipole operators----------------------------------------------
t=1/math.sqrt(2)

Akm = {{1,-1,t},{1, 1,-t}}
T2p3dx = NewOperator("CF", NF, Index["3d_Up"],Index["3d_Dn"],Index["2p_Up"],Index["2p_Dn"],  Akm)

Akm = {{1,-1,t*I},{1, 1,t*I}}
T2p3dy = NewOperator("CF", NF, Index["3d_Up"],Index["3d_Dn"],Index["2p_Up"],Index["2p_Dn"], Akm)

Akm = {{1,0,1}}
T2p3dz = NewOperator("CF", NF, Index["3d_Up"],Index["3d_Dn"],Index["2p_Up"],Index["2p_Dn"], Akm)

--construct combination
T2p3dypz  =  sqrt(1/2)*(T2p3dy + T2p3dz)
T2p3dzpx  =  sqrt(1/2)*(T2p3dz + T2p3dx)
T2p3dxpy  =  sqrt(1/2)*(T2p3dx + T2p3dy)
T2p3dypiz =  sqrt(1/2)*(T2p3dy + I * T2p3dz)
T2p3dzpix =  sqrt(1/2)*(T2p3dz + I * T2p3dx)
T2p3dxpiy =  sqrt(1/2)*(T2p3dx + I * T2p3dy)

--for method 1
T2p3dM1Ls={T2p3dx,T2p3dy,T2p3dz}
T3d2pM1Ls={ConjugateTranspose(T2p3dx),ConjugateTranspose(T2p3dy),ConjugateTranspose(T2p3dz)}

--for method 2
T2p3dM2Ls={T2p3dx,T2p3dy,T2p3dz,T2p3dxpy,T2p3dypz,T2p3dzpx,T2p3dxpiy,T2p3dypiz,T2p3dzpix}
T3d2pM2Ls={ConjugateTranspose(T2p3dx),ConjugateTranspose(T2p3dy),ConjugateTranspose(T2p3dz),ConjugateTranspose(T2p3dxpy),ConjugateTranspose(T2p3dypz),ConjugateTranspose(T2p3dzpx),ConjugateTranspose(T2p3dxpiy),ConjugateTranspose(T2p3dypiz),ConjugateTranspose(T2p3dzpix)}




---------------------------------------------RIXS---------------------------------------------------
--explain the inputs and output of the following function
--CalculateRIXSTensor(Hamiltonian,XASHamiltonian,finHamiltonian,dirSymm,n,Neig)
--Hamiltonian: the Hamiltonian of the system
--XASHamiltonian: the Hamiltonian of the system for the XAS
--finHamiltonian: the Hamiltonian of the system for the final state
--dirSymm: the symmetry directory of the system
--n: the number of electrons in the 3d shell
--Neig: the number of eigenstates to compute in Eigensystem
function CalculateRIXSTensor(Hamiltonian,XASHamiltonian,finHamiltonian,dirSymm,n,Neig)
    print("Calculate RIXS...")

    --Observable of Interests.
    oppList={Opp3dJz, Opp3dJsqr, Opp3dLz, Opp3dLsqr, Opp3dSz, Opp3dSsqr,OpU3d}

    --Number of bosonic states in the basis
    NB=0
    --start Restrictions for Eigensystem
    StartRestrictions= {NF, NB, {"1111111111 111111",6+n,6+n},{"1111111111 000000",n,n}}


    -- we now can create the lowest Npsi eigenstates, now with Higher sensibility (The shift of the GS energy to -1 might improve the convergence, but the algorithm should already correct it):
    TimeStart("Eigensystem")
    psiList = Eigensystem(Hamiltonian, StartRestrictions, Neig,{{"Epsilon",10^(-12)},{"Zero",10^(-12)}})
    TimeEnd("Eigensystem")

    --if Neig is 1 make a list
    if Neig==1 then
        psiList={psiList}
    end

    --Ground state energy
    E0 = psiList[1] * Hamiltonian * psiList[1]

    -- Find the number of degenerated ground states
    GS_list = {}
    
    for i, psi in ipairs(psiList) do
        local E = psi * Hamiltonian * psi
        if math.abs(E - E0) < 1e-5 then
            GS_list[#GS_list+1] = psi
        end
    end

    

    print("Number of degenerated ground states: " .. #GS_list)
    print("The indeces for the considered starting states are:")



    --Print The Expectations values
    PrintExpectationValues(psiList,oppList, Hamiltonian,{{"printSigma",true}})
    print("test out tensor option")
    
    --perform RIXS calculation Method1-use the brand new Tensor Option
    RIXSChispectraM1Cart = 0
    --sum over the considered GS's
    for j_GS,psi_GS in pairs(GS_list) do 
        TimeStart("Resonant Tensor Spectra method 1")
        RIXSChispectraM1Cart = 1/(#GS_list)*CreateResonantSpectra(XASHamiltonian,finHamiltonian, T2p3dM1Ls, T3d2pM1Ls , psi_GS,{{"Emin1",E01}, {"Emax1",E01+dE1}, {"NE1",NE1}, {"Gamma1",Gamma1}, {"NTri1",NTri1}, {"Emin2",E02}, {"Emax2",E02+dE2}, {"NE2",NE2}, {"Gamma2",Gamma2},{"NTri2",NTri2}, {"tensor",true},{"denseborder",1000000},{"TensorBasis","cartesian"}}) +RIXSChispectraM1Cart
        TimeEnd("Resonant Tensor Spectra method 1")
        TimePrint()
    end   
    TimePrint()


    RIXSChispectraM1Cart = Chop(RIXSChispectraM1Cart)


    --perform RIXS calculation Method2-Compute all components and use the SumTensorNE1 to combine them into the cartesian basis
    RIXSSpectraRaw = 0
    --sum over the considered GS's
    for j_GS,psi_GS in pairs(GS_list) do  
        TimeStart("Resonant Tensor Spectra method 2")
        RIXSSpectraRaw = RIXSSpectraRaw+1/(#GS_list)*CreateResonantSpectra(XASHamiltonian,finHamiltonian, T2p3dM2Ls, T3d2pM2Ls , psi_GS,{{"Emin1",E01}, {"Emax1",E01+dE1}, {"NE1",NE1}, {"Gamma1",Gamma1}, {"NTri1",NTri1}, {"Emin2",E02}, {"Emax2",E02+dE2}, {"NE2",NE2}, {"Gamma2",Gamma2},{"NTri2",NTri2},{"denseborder",1000000},{"TensorBasis","cartesian"}})
        TimeEnd("Resonant Tensor Spectra method 2")
        TimePrint()
    end  


    --Total Sum Tensor for the Spectra.sum() function. The imported SumTensorE1E1Cart matrix gives us the general rule to combine the computed components
    --into the cartesian ones of the dipole operator.
    SumTensorNE1={}
    for i=1,#SumTensorE1E1Cart do
        for k=1,NE1+1 do
            SumTensorNE1[(i-1)*(NE1+1)+k]={}
            for j=1,#SumTensorE1E1Cart[i] do
                for l=1,NE1+1 do
                    if l==k then
                        --we had to insert the conjugate and is an hint that in the paper the indexes re inverted.
                        SumTensorNE1[(i-1)*(NE1+1)+k][(j-1)*(NE1+1)+l]=Conjugate(SumTensorE1E1Cart[i][j])
                    else
                        SumTensorNE1[(i-1)*(NE1+1)+k][(j-1)*(NE1+1)+l]=0
                    end
                end
            end 
        end  
    end 
    print("Total size spectra=81*(NE1+1)=",81*(NE1+1),#SumTensorNE1,#SumTensorNE1[1])

    TimeStart("Sum Into RIXS Spectra")

    RIXSChispectraM2Cart = Chop(Spectra.Sum(RIXSSpectraRaw, table.unpack(SumTensorNE1)))
    RIXSChispectraM2Cart.NE1 = NE1
    RIXSChispectraM2Cart.TensorBasis = "cartesian"

    --print method 2 cartesian spectra to file
    


    print("Check that the two methods give the same result in cartesian basis")
    TestIfSpectrumIsZero(Chop(RIXSChispectraM1Cart-RIXSChispectraM2Cart,accuracy))


    
    --Total Sum Tensor for the Spectra.sum() function. The imported SumTensorE1E1Cart matrix gives us the general rule to combine the computed components
    --into the cartesian ones of the dipole operator.
    SumTensorNE1={}
    for i=1,#SumTensorE1E1KHarmonics do
        for k=1,NE1+1 do
            SumTensorNE1[(i-1)*(NE1+1)+k]={}
            for j=1,#SumTensorE1E1KHarmonics[i] do
                for l=1,NE1+1 do
                    if l==k then
                        SumTensorNE1[(i-1)*(NE1+1)+k][(j-1)*(NE1+1)+l]=Conjugate(SumTensorE1E1KHarmonics[i][j])
                    else
                        SumTensorNE1[(i-1)*(NE1+1)+k][(j-1)*(NE1+1)+l]=0
                    end
                end
            end 
        end  
    end 
    print("Total size spectra=81*(NE1+1)=",81*(NE1+1),#SumTensorNE1,#SumTensorNE1[1])

    TimeStart("Sum Into RIXS Spectra")
    RIXSChispectraM2SphK = Spectra.Sum(RIXSSpectraRaw, table.unpack(SumTensorNE1))
    TimeEnd("Sum Into RIXS Spectra")
    TimePrint()

    --rotate the RIXS Cartesian tensor to the Spherical Tensor K basis


    print("-------------------------------Start Tests---------------------------------")
    print("---------------------------------------------------------------------------")
    print("-Rotate Spectra Test-------------------------------------------------------")
    print("---------------------------------------------------------------------------")

    UCartToSphKQuanty=CoupledCartesianToSphericalTensorKMatrix(1,1)
    UCartToSphK = {
        {-(1/math.sqrt(3)), 0, 0, 0, -(1/math.sqrt(3)), 0, 0, 0, -(1/math.sqrt(3))},
        {0, 0, 0, 0, 0, I/math.sqrt(2), 0, -I/math.sqrt(2), 0},
        {0, 0, -I/math.sqrt(2), 0, 0, 0, I/math.sqrt(2), 0, 0},
        {0, I/math.sqrt(2), 0, -I/math.sqrt(2), 0, 0, 0, 0, 0},
        {1/math.sqrt(2), 0, 0, 0, -(1/math.sqrt(2)), 0, 0, 0, 0},
        {-(1/math.sqrt(6)), 0, 0, 0, -(1/math.sqrt(6)), 0, 0, 0, math.sqrt(2/3)},
        {0, 0, 0, 0, 0, 1/math.sqrt(2), 0, 1/math.sqrt(2), 0},
        {0, 0, 1/math.sqrt(2), 0, 0, 0, 1/math.sqrt(2), 0, 0},
        {0, 1/math.sqrt(2), 0, 1/math.sqrt(2), 0, 0, 0, 0, 0}
    }
    print("Check that the defined UCartToSphKQuanty matrix matches the manually defined UCartToSphK matrix")
    print(UCartToSphKQuanty-UCartToSphK)
    

    --double rotation
    RIXSChispectraM1Rotated = Rotate(RIXSChispectraM1Cart,UCartToSphKQuanty)
    RIXSChispectraM1Rotated = Rotate(RIXSChispectraM1Rotated,ConjugateTranspose(UCartToSphKQuanty))


    print("Subsequent rotation in opposite directions should leave the matrix unchanged")
    TestIfSpectrumIsZero(Chop(RIXSChispectraM1Cart-RIXSChispectraM1Rotated,accuracy))


    --rotation by identity
    RIXSChispectraM1Id = Rotate(RIXSChispectraM1Cart,Matrix.Identity(9))
    print("Rotation by id() should leave the matrix unchanged")
    TestIfSpectrumIsZero(Chop(RIXSChispectraM1Cart-RIXSChispectraM1Id,accuracy))

    
    --rotate the RIXS Cartesian tensor to the Spherical Tensor K basis
    RIXSChispectraM1SphK = Chop(Rotate(RIXSChispectraM1Cart,Conjugate(UCartToSphK),{{"TensorBasis","SphericalTensorK"}}))


    --check difference method1 and method2 for the SphK basis
    print("The tensors for the two methods in SphK basis should be equivalent")
    TestIfSpectrumIsZero(Chop(RIXSChispectraM1SphK-RIXSChispectraM2SphK,accuracy))
    


    print("---------------------------------------------------------------------------")
    print("--Powder function GetPowderTensor()----------------------------------------")
    print("---------------------------------------------------------------------------")

    --Total Sum Tensor for the Spectra.sum() function. The imported SumTensorE1E1Cart matrix gives us the general rule to combine the computed components
    --into the cartesian ones of the dipole operator.
    SumTensorNE1Pwd={}
    l_vec={{0},{1,2,3},{4,5,6,7,8}}

    for i=1,81 do
        for k=1,NE1+1 do
            SumTensorNE1Pwd[(i-1)*(NE1+1)+k]={}
            for j=1,81 do
                for l=1,NE1+1 do
                    SumTensorNE1Pwd[(i-1)*(NE1+1)+k][(j-1)*(NE1+1)+l]=0
                   
                    row_i = math.floor((i-1)/9)
                    row_j = math.floor((j-1)/9)
                    col_i = (i-1)%9
                    col_j = (j-1)%9
                    if col_i==row_i then
                        if col_j==row_j then
                            if k==l then
                                for lang=1,#l_vec do
                                    for m=1,#l_vec[lang] do
                                        for n=1,#l_vec[lang] do
                                            --it i=j and i in l_vec[lang]
                                            if row_i==l_vec[lang][m] and row_j == l_vec[lang][n] and k==l then
                                                --print("i=",i,"j=",j,"row_i=",row_i,"row_j=",row_j,"col_i=",col_i,"col_j=",col_j)
                                                SumTensorNE1Pwd[(i-1)*(NE1+1)+k][(j-1)*(NE1+1)+l]=1/#l_vec[lang]
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                    
                end
            end 
        end  
    end 


    TimeStart("Sum Into RIXS tensor powder")
    RIXSChispectraM1SphKPwdSum = Spectra.Sum(RIXSChispectraM1SphK, table.unpack(SumTensorNE1Pwd))

    --Manual function
    RIXSChispectraM1SphKPowder=RIXSChispectraM1SphK.GetPowderTensor()
    RIXSChispectraM1CartPowder=RIXSChispectraM1Cart.GetPowderTensor()   
    
    print("Powder tensor from cartesian basis and spherical basis should be equal")
    temp=Rotate(RIXSChispectraM1SphKPowder,ConjugateTranspose(UCartToSphK),{{"TensorBasis","cartesian"}})
    TestIfSpectrumIsZero(Chop(RIXSChispectraM1CartPowder-temp,accuracy))


   
    print("compare powdered by SUM Matrix and powdered by GetPowderTensor()")
    TestIfSpectrumIsZero(Chop(RIXSChispectraM1SphKPwdSum-RIXSChispectraM1SphKPowder,accuracy))

    

    
    print("---------------------------------------------------------------------------")
    print("--ContractSpectraTensorWithPolarization()----------------------------------")
    print("---------------------------------------------------------------------------")
    
  
    --Contract the tensor with fixed polarization vectors; in test k the incoming vector is
    --number k and the outgoing one number k+1, so every vector is used both in and out
    for kPol = 1, #PolarizationVectors do
        Ein = PolarizationVectors[kPol]
        Eout = PolarizationVectors[kPol % #PolarizationVectors + 1]
        print("Polarization test "..kPol.." of "..#PolarizationVectors)
    
        ecoupled ={}
        for i=1,3 do
            for j=1,3 do
                ecoupled[(i-1)*3+j]=Ein[i]*Conjugate(Eout[j])
            end
        end
        print("The polarization vectors are:")
        print("Ein=",Ein)
        print("Eout=",Eout)


        T2p3dIn=Ein[1]*T2p3dx+Ein[2]*T2p3dy+Ein[3]*T2p3dz
        T3d2pOut=Conjugate(Eout[1])*ConjugateTranspose(T2p3dx)+Conjugate(Eout[2])*ConjugateTranspose(T2p3dy)+Conjugate(Eout[3])*ConjugateTranspose(T2p3dz)


        RIXSPolarization = 0
        --sum over the considered GS's
        for j_GS,psi_GS in pairs(GS_list) do 
            TimeStart("Resonant Spectra with polarization vectors")
            RIXSPolarization = 1/(#GS_list)*CreateResonantSpectra(XASHamiltonian,finHamiltonian, T2p3dIn, T3d2pOut , psi_GS,{{"Emin1",E01}, {"Emax1",E01+dE1}, {"NE1",NE1}, {"Gamma1",Gamma1}, {"NTri1",NTri1}, {"Emin2",E02}, {"Emax2",E02+dE2}, {"NE2",NE2}, {"Gamma2",Gamma2},{"NTri2",NTri2},{"denseborder",100000}}) +RIXSPolarization
            TimeEnd("Resonant Spectra with polarization vectors")
        end   
    
        TimePrint()
        print("The coupled polarization vectors are:")
        print(ecoupled)

        print("ContractSpectraTensorWithPolarization() test")
        RIXSPolarizationfromTensorCart = RIXSChispectraM1Cart.ContractSpectraTensorWithPolarization(Ein,Eout)
        print("contraction done with the tensor function")


        print("The coupled polarization vectors in the SphericalTensorK basis are:")
        --rotate ecoupled by multiplyng by UcarttoSphk
        ecoupledSphK = {}
        for i=1,9 do
            ecoupledSphK[i] = 0
            for k=1,9 do
                ecoupledSphK[i] = ecoupledSphK[i] + UCartToSphK[i][k]*ecoupled[k]
            end
            print("ecoupledSphK[",i,"]=",ecoupledSphK[i])
        end


    
        print("ContractSpectraTensorWithPolarization() test")
        RIXSPolarizationfromTensorSphK = RIXSChispectraM1SphK.ContractSpectraTensorWithPolarization(Ein,Eout)
        print("contraction done with the tensor function")

    
        --Total Sum Tensor for the Spectra.sum() function. The imported SumTensorE1E1Cart matrix gives us the general rule to combine the computed components
        --into the cartesian ones of the dipole operator.
        SumTensorNE1Pol={}

        for k=1,NE1+1 do
            SumTensorNE1Pol[k]={}
            for l=1,NE1+1 do
                for i=1,9 do
                    for j=1,9 do
                        idxpol = (i-1)*9+j-1
                        if l==k then 
                            SumTensorNE1Pol[k][idxpol*(NE1+1)+l]=Conjugate(ecoupled[i])*ecoupled[j]
                        else
                            SumTensorNE1Pol[k][idxpol*(NE1+1)+l]=0
                        end
                    end
                end 
            end  
        end 

        print("Total size spectra=81*(NE1+1)=",81*(NE1+1),#SumTensorNE1Pol,#SumTensorNE1Pol[1])

        TimeStart("Sum using crontraction with polarization vectors")
        RIXSPolarizationfromSum = Spectra.Sum(RIXSChispectraM1Cart, table.unpack(SumTensorNE1Pol))

        print("compare Polarization raw and FromTensorCart")
        TestIfSpectrumIsZero(Chop(RIXSPolarization-RIXSPolarizationfromTensorCart,accuracy))
    
        print("compare Polarization raw and Sum")
        TestIfSpectrumIsZero(Chop(RIXSPolarization-RIXSPolarizationfromSum,accuracy))
    
        print("compare Polarization FromTensorCart and Sum")
        TestIfSpectrumIsZero(Chop(RIXSPolarizationfromTensorCart-RIXSPolarizationfromSum,accuracy))
    
        print("compare Polarization From FromTensorCart and FromTensorSphK")
        TestIfSpectrumIsZero(Chop(RIXSPolarizationfromTensorCart-RIXSPolarizationfromTensorSphK,accuracy))
    end



end










--Calculate the RIXS Tensor for the different symmetries
CalculateRIXSTensor(HamiltonianSO3,XASHamiltonianSO3,finHamiltonianSO3,"SO3",n,Neig)
--CalculateRIXSTensor(HamiltonianOh,XASHamiltonianOh,finHamiltonianOh,"Oh",n,Neig)
CalculateRIXSTensor(HamiltonianD4h,XASHamiltonianD4h,finHamiltonianD4h,"D4h",n,Neig)








--give the Informations about the calculation
io.write("\n")
io.write("---------------------------------end of the simulation----------------------------------------")
io.write("The simulation took:")
TimePrint()


io.flush()