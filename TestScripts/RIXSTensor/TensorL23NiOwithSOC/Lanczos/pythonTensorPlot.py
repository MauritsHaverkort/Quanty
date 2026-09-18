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

dE1 = 0.2
NE1 = 4
E01 = -0.1
dE2 = 8
NE2 = 400
E02 = -2
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

dirData="data/n3d6"

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
    plotDir="./plots/n3d6/"+SymmName+"/"

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
