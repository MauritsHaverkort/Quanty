--Verbosity(0xFFFF)
--Import the List to combine Raw spectra into cartesian RIXS tensor. this list has been generated with
--the script ComputeRixsTensor.nb in the parent directory.
dofile("./SumTensorE1E1Cart.lua")

--define functions
function slice(tbl, first, last, step)
    local sliced = {}
    for i = first or 1, last or #tbl, step or 1 do
        sliced[#sliced + 1] = tbl[i]
    end
    return sliced
end


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

Bz=0.0000001
--formal valence 3d shell
n = 8


Neig = 100
if n==0 then
    Neig=1
end

------------------------------Parameters for the spectras-------------------------

--size of Krylov basis (std 200)
NTri1=100
NTri2=200

--Gamma forthe maps 
Gamma1=0.5 --I enhanced this since se generate lot of data
Gamma2=0.2

--RIXS energies
dE1=1
E01_0=-0.2

dE2=8
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



---------------------------------Some bureaucracy-------------------------------------------
--Set and create data folders
Dir=""--directory where to data in
Dir_par=string.format("n3d%i",n)--create subfolder for the parameters selected
os.execute("mkdir "..Dir.."data;mkdir "..Dir.."data/"..Dir_par)
dirData=string.format("%sdata/%s", Dir,Dir_par)--path data



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
    
    --perform RIXS calculation Method2
    RIXSChispectraM1 = 0
    --sum over the considered GS's
    for j_GS,psi_GS in pairs(GS_list) do 
        TimeStart("Resonant Tensor Spectra method 1")
        RIXSChispectraM1 = CreateResonantSpectra(XASHamiltonian,finHamiltonian, T2p3dM1Ls, T3d2pM1Ls , psi_GS,{{"Emin1",E01}, {"Emax1",E01+dE1}, {"NE1",NE1}, {"Gamma1",Gamma1}, {"NTri1",NTri1}, {"Emin2",E02}, {"Emax2",E02+dE2}, {"NE2",NE2}, {"Gamma2",Gamma2},{"NTri2",NTri2}, {"tensor",true},{"denseborder",10000},{"TensorBasis","cartesian"}})
        TimeEnd("Resonant Tensor Spectra method 1")
        TimePrint()
    end   
    TimePrint()

    dataChispectraM1=string.format("%s/%s/chiRIXSTensorMethod1.dat", dirData,dirSymm)--path data
    dataChispectraM1Rot=string.format("%s/%s/chiRIXSTensorMethod1Rot.dat", dirData,dirSymm)--path data
    
    os.execute("mkdir -p "..dirData.."/"..dirSymm)
    print("Write RIXS Tensor to file: "..dataChispectraM1)
    RIXSChispectraM1.Print({{"file", dataChispectraM1}})


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

    RIXSChispectraM1Rotated = Rotate(RIXSChispectraM1,UCartToSphK)
    RIXSChispectraM1Rotated.Print({{"file", dataChispectraM1Rot}})


    --perform RIXS calculation Method2
    RIXSSpectraRaw = 0
    --sum over the considered GS's
    for j_GS,psi_GS in pairs(GS_list) do  
        TimeStart("Resonant Tensor Spectra method 2")
        RIXSSpectraRaw = RIXSSpectraRaw+1/(#GS_list)*CreateResonantSpectra(XASHamiltonian,finHamiltonian, T2p3dM2Ls, T3d2pM2Ls , psi_GS,{{"Emin1",E01}, {"Emax1",E01+dE1}, {"NE1",NE1}, {"Gamma1",Gamma1}, {"NTri1",NTri1}, {"Emin2",E02}, {"Emax2",E02+dE2}, {"NE2",NE2}, {"Gamma2",Gamma2},{"NTri2",NTri2},{"denseborder",10000},{"TensorBasis","cartesian"}})
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
                        SumTensorNE1[(i-1)*(NE1+1)+k][(j-1)*(NE1+1)+l]=SumTensorE1E1Cart[i][j]
                    else
                        SumTensorNE1[(i-1)*(NE1+1)+k][(j-1)*(NE1+1)+l]=0
                    end
                end
            end 
        end  
    end 
    print("Total size spectra=81*(NE1+1)=",81*(NE1+1),#SumTensorNE1,#SumTensorNE1[1])

    TimeStart("Sum Into RIXS Spectra")
    RIXSChispectraM2 = Spectra.Sum(RIXSSpectraRaw, table.unpack(SumTensorNE1))

    --Paths where to save the output of the calculations
    dataChispectraM2=string.format("%s/%s/chiRIXSTensorMethod2.dat", dirData,dirSymm)--path data
    RIXSChispectraM2.Print({{"file", dataChispectraM2}})
    
    TimeEnd("Sum Into RIXS Spectra")
    TimePrint()

    
end


--Calculate the RIXS Tensor for the different symmetries
CalculateRIXSTensor(HamiltonianSO3,XASHamiltonianSO3,finHamiltonianSO3,"SO3",n,Neig)
CalculateRIXSTensor(HamiltonianOh,XASHamiltonianOh,finHamiltonianOh,"Oh",n,Neig)
CalculateRIXSTensor(HamiltonianD4h,XASHamiltonianD4h,finHamiltonianD4h,"D4h",n,Neig)


--give the Informations about the calculation
io.write("\n")
io.write("---------------------------------end of the simulation----------------------------------------\n")
io.write("The simulation took:\n")
TimePrint()


pythonInput1 = [[
#!/usr/bin/env python
# coding: utf-8

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import SymLogNorm

# Import APS style
# plt.style.use("../../Styles/APS.py")

# settings nature paper
cm = 1/2.54  # centimeters in inches

# Figsizept
sizeSingleColAPS = 8.5 * cm

#String for the cartesian basis
stringBasisCartesian = ["xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz"]


########################## Definitions ###########################

## Load some Parameters to create the energy grid.

dE1 = ]]..dE1..[[

NE1 = ]]..NE1..[[

E01 = ]]..E01_0..[[

dE2 = ]]..dE2..[[

NE2 = ]]..NE2..[[

E02 = ]]..E02_0..[[

# create for E1 and E2 the energy mesh
E1 = np.linspace(E01, E01 + dE1, NE1 + 1)
E2 = np.linspace(E02, E02 + dE2, NE2 + 1)

# Functions 
UCartToSphK = np.array([
  [-(1/np.sqrt(3)), 0, 0, 0, -(1/np.sqrt(3)), 0, 0, 0, -(1/np.sqrt(3))],
  [0, 0, 0, 0, 0, 1j/np.sqrt(2), 0, -1j/np.sqrt(2), 0],
  [0, 0, -1j/np.sqrt(2), 0, 0, 0, 1j/np.sqrt(2), 0, 0],
  [0, 1j/np.sqrt(2), 0, -1j/np.sqrt(2), 0, 0, 0, 0, 0],
  [1/np.sqrt(2), 0, 0, 0, -(1/np.sqrt(2)), 0, 0, 0, 0],
  [-(1/np.sqrt(6)), 0, 0, 0, -(1/np.sqrt(6)), 0, 0, 0, np.sqrt(2/3)],
  [0, 0, 0, 0, 0, 1/np.sqrt(2), 0, 1/np.sqrt(2), 0],
  [0, 0, 1/np.sqrt(2), 0, 0, 0, 1/np.sqrt(2), 0, 0],
  [0, 1/np.sqrt(2), 0, 1/np.sqrt(2), 0, 0, 0, 0, 0]
])


def Rotate(RIXSTensor,U):
    RIXSTensorOut=[]
    for i in range(len(RIXSTensor)):
        RIXSTensorOut.append([])
        for j in range(len(RIXSTensor)):
            elementij=0
            for m in range(len(RIXSTensor)):
                for n in range(len(RIXSTensor)):
                    elementij+=np.conjugate(U[i][m])*RIXSTensor[m][n]*U[j][n]
            RIXSTensorOut[i].append(elementij)
    return RIXSTensorOut




# Read the Quanty Output
# loop over the stringbasis
def readRIXSData(fileName, stringBasis):


    dataDir = os.path.dirname(fileName) + "/"
    fileDataRe = dataDir + "RIXS_Re.dat"
    fileDataIm = dataDir + "RIXS_Im.dat"

    # Transform the Quanty Output to a Readable ASCII Table.
    os.system("awk '((NR>5)&&(NR<" + str((NE2 + 6) * 1) + ")) {for(i=2;i<=NF;i=i+2){printf \"%s \", $i}{printf \"%s\", RS}}' " + fileName + " >" + fileDataRe)
    os.system("awk '((NR>5)&&(NR<" + str((NE2 + 6) * 1) + ")) {for(i=3;i<=NF;i=i+2){printf \"%s \", $i}{printf \"%s\", RS}}' " + fileName + " >" + fileDataIm)

    RIXS = np.loadtxt(fileDataRe, dtype=complex)
    RIXS.imag = np.loadtxt(fileDataIm)

    print("Loaded RIXS data shape:", RIXS.shape)
    #NE1 = RIXS.shape[1] // (len(stringBasis) * len(stringBasis))
    #NE2 = RIXS.shape[0] 
    RIXSTensor = np.empty((len(stringBasis), len(stringBasis), NE2+1, NE1+1), dtype=complex)


    for keyI, epsOut in enumerate(stringBasis):
        for keyJ, epsIn in enumerate(stringBasis):
                for i in range(NE1+1):
                    index = i + (NE1+1) * (keyJ + len(stringBasis) * keyI)
                    for j in range(NE2):
                        RIXSTensor[keyI][keyJ][j][i]=RIXS[j,index]
    
    print("Loaded RIXS Tensor shape:", RIXSTensor.shape)
    #clean the temporary files
    os.system("rm {}".format(fileDataRe))
    os.system("rm {}".format(fileDataIm))

    return RIXSTensor
########################## Plotting fuinctions the RIXS ###########################

#function for plotting the RIXS Tensor in separated files for each real and imaginary part.
def MakeTensorPlot(RIXSTensor, stringBasis, basisName, plotDir,norm=1,vmax=0.2):

    mapName = "coolwarm"
    Nbasis = len(stringBasis)
    #linear range for the color norm
    linthresh=0.001

    #loop over the real and imaginary part
    for i in range(2):

        ## generate a 9x9 multiplot where each subplot is one element of the RIXSTensor
        fig, axs = plt.subplots(Nbasis, Nbasis, figsize=(2*sizeSingleColAPS, 2*sizeSingleColAPS), constrained_layout=False)

        #outer loop of the polarization back to the Ground state
        for keyI, epsOut in enumerate(stringBasis):
            #inner loop for the polarization to the final state
            for keyJ, epsIn in enumerate(stringBasis):
                if i == 0:
                    cmap = axs[keyI, keyJ].pcolormesh(E1, E2, (1/norm) * RIXSTensor[keyI][keyJ].real, cmap=plt.get_cmap(mapName),
                    linewidth=0, shading="auto",norm=SymLogNorm(linthresh=linthresh, vmin=-vmax, vmax=vmax))
                    fig.suptitle(r"Re($\chi_{RIXS}$) in " + basisName)
                if i == 1:
                    cmap = axs[keyI, keyJ].pcolormesh(E1, E2, (1/norm) * RIXSTensor[keyI][keyJ].imag, cmap=plt.get_cmap(mapName),
                    linewidth=0, shading="auto",norm=SymLogNorm(linthresh=linthresh, vmin=-vmax, vmax=vmax))
                    fig.suptitle(r"Im($\chi_{RIXS}$) in " + basisName)

                # label top row and left column with the basis
                if keyI == 0:
                    axs[keyI, keyJ].xaxis.set_label_position("top")
                    axs[keyI, keyJ].set_xlabel(epsIn)
                if keyJ == len(stringBasis) - 1:
                    axs[keyI, keyJ].yaxis.set_label_position("right")
                    axs[keyI, keyJ].set_ylabel(epsOut)

        # layout space between the subplots and only leave the outermost axes labels
        plt.subplots_adjust(wspace=0., hspace=0.)
        # Adjust the space between the plot and the suptitle
        fig.subplots_adjust(top=0.92)

        # Set a single label for the entire plot
        fig.text(0.51, 0.27, "Relative excitation energy, $\omega_{in}$ (eV)", ha='center', va='center')
        fig.text(0.075, 0.63, "Energy transfer, $\omega$ (eV)", ha='center', va='center', rotation='vertical')

        # drop all tickslabels but the outer ones
        for ax in axs.flat:
            ax.label_outer()
            ax.tick_params(axis='x', direction='in', length=2)  # Customize as needed
            ax.tick_params(axis='y', direction='in', length=2)

        # add cbar to the figure
        cbar = fig.colorbar(cmap, ax=axs, orientation='horizontal', shrink=0.5, pad=0.1)
        cbar.set_label('Intensity (arb. units)', fontsize=10)

        if i == 0:
            fileName = plotDir + "RIXS_Re_" + basisName

        if i == 1:
            fileName = plotDir + "RIXS_Im_" + basisName

        os.makedirs(plotDir, exist_ok=True)

        fig.savefig(fileName.replace(" ", ""))
        plt.close(fig)


############################################################################################################
###################################Actual Plotting#########################################################
############################################################################################################

dirData="]]..dirData..[["

stringBasisSphericalKread=["0s","1Rx","1Ry","1Rz","2x2y2","2z2","2yz","2xz","2xy"]
stringBasisSphericalY=["s","R$_{-1}$","R$_{0}$","R$_{1}$","d$_{-2}$","d$_{-1}$","d$_{0}$","d$_{1}$","d$_{2}$"]

vmaxGlobal=0.1

def readandPlotComparisonMethods(SymmName):

    stringBasisCartesian=["xx","xy","xz","yx","yy","yz","zx","zy","zz"]
    print("Plotting Symmetry: ",SymmName)

    #Read input tensor
    fileNameM1 = dirData+"/"+SymmName+"/chiRIXSTensorMethod1.dat"
    fileNameM2 = dirData+"/"+SymmName+"/chiRIXSTensorMethod2.dat"

    RIXSTensorCartM1=readRIXSData(fileNameM1, stringBasisCartesian)
    RIXSTensorCartM2=readRIXSData(fileNameM2, stringBasisCartesian)

    RIXSTensorSphKM1=Rotate(RIXSTensorCartM1, UCartToSphK)
    RIXSTensorSphKM2=Rotate(RIXSTensorCartM2, UCartToSphK)

    #plot
    plotDir="./plots/n3d]]..n..[[/"+SymmName+"/"

    #single tensor plots
    MakeTensorPlot(RIXSTensorCartM1, stringBasisCartesian,"Cartesian Method1 "+SymmName,plotDir)
    MakeTensorPlot(RIXSTensorCartM2, stringBasisCartesian,"Cartesian Method2 "+SymmName,plotDir)

    MakeTensorPlot(RIXSTensorSphKM1, stringBasisSphericalK,"Spherical K Method1 "+SymmName,plotDir)
    MakeTensorPlot(RIXSTensorSphKM2, stringBasisSphericalK,"Spherical K Method2 "+SymmName,plotDir)

    return RIXSTensorCartM1,RIXSTensorCartM2

#SO3
stringBasisSO3SphericalK=["s","R$_x$","R$_y$","R$_z$","d$_{x^2-y^2}$","d$_{z^2}$","d$_{yz}$","d$_{xz}$","d$_{xy}$"]
stringBasisCartesian=["xx","xy","xz","yx","yy","yz","zx","zy","zz"]
stringBasisSphericalK=stringBasisSO3SphericalK
RIXSTensorSO3CartM1,RIXSTensorSO3CartM2=readandPlotComparisonMethods("SO3")

#Oh
#stringBasisSphericalK=stringBasisSO3SphericalK
RIXSTensorOhCartM1,RIXSTensorOhCartM2=readandPlotComparisonMethods("Oh")

#D4h
#stringBasisSphericalK=stringBasisSO3SphericalK
RIXSTensorD4hCartM1,RIXSTensorD4hCartM2=readandPlotComparisonMethods("D4h")
]]


-- write the gnuplot script to a fil
file = io.open("pythonTensorPlot.py", "w")
file:write(pythonInput1)
file:close()

--excecute, this will take a while since they are many heavy plots
os.execute("python3 pythonTensorPlot.py")
