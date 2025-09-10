import sys
import os
import glob
import random
import math
import numpy as np
import mathUtility
import getA
def cvrtSpeciesID_LiggghtsToPython(SpeciesID):
    return SpeciesID - 1
def cvrtSpeciesID_PythonToLiggghts(SpeciesID):
    return SpeciesID + 1
        
def computeForceCoefficient(nparticletype, pipeSpecies, textureSpecies, YoungsModulus, PoissonRatio, CoefficientOfRestitution, CharacteristicVelocity, DensityArray, RadiusArray):
    speciesNum = nparticletype
    kn_specified_Matrix = np.zeros((speciesNum, speciesNum))
    kt_specified_Matrix = np.zeros((speciesNum, speciesNum))
    gamman_specified_Matrix= np.zeros((speciesNum, speciesNum))
    gammat_specified_Matrix = np.zeros((speciesNum, speciesNum))
    for i in range(speciesNum):
        for j in range(speciesNum):
            #reduced Young modulus
            Yeff = 1.0/((1.0-PoissonRatio*PoissonRatio)/YoungsModulus + (1.0-PoissonRatio*PoissonRatio)/YoungsModulus)
            #reduced Shear modulus
            Geff = 0.5/((2.0-PoissonRatio)*(1.0+PoissonRatio)/YoungsModulus + (2.0-PoissonRatio)*(1.0+PoissonRatio)/YoungsModulus)
            if(i==pipeSpecies): # infinite mass, infinite radius
                Reff = RadiusArray[j]
                Meff = DensityArray[j]*mathUtility.computeSphereVolume(RadiusArray[j])  
            elif(j==pipeSpecies):
                Reff = RadiusArray[i]
                Meff = DensityArray[i]*mathUtility.computeSphereVolume(RadiusArray[i])
            elif(i==textureSpecies): # infinite mass, finite radius
                Reff = 1.0/(1.0/RadiusArray[i] + 1.0/RadiusArray[j])
                Meff = DensityArray[j]*mathUtility.computeSphereVolume(RadiusArray[j])
            elif(j==textureSpecies):
                Reff = 1.0/(1.0/RadiusArray[i] + 1.0/RadiusArray[j])
                Meff = DensityArray[i]*mathUtility.computeSphereVolume(RadiusArray[i])
            else:
                m = np.zeros((2,1))
                m[0] = DensityArray[i]*mathUtility.computeSphereVolume(RadiusArray[i])
                m[1] = DensityArray[j]*mathUtility.computeSphereVolume(RadiusArray[j])
                Reff = 1.0/(1.0/RadiusArray[i] + 1.0/RadiusArray[j])
                Meff = 1.0/(1.0/m[0] + 1.0/m[1])
            
            # Need verification for the computation of A_n
            density = DensityArray[0]
            An =  getA.getAPhysical(CoefficientOfRestitution,Reff,Reff,density,YoungsModulus,PoissonRatio,CharacteristicVelocity)
            # https://www.cfdem.com/media/DEM/docu/gran_model_hertz_stiffness.html
            kn_specified_Matrix[i][j] = (4/3)* Yeff 
            kt_specified_Matrix[i][j] = 8.0*Geff
            gamman_specified_Matrix[i][j] = 2.0*Yeff*An/Meff # not liggghts definition, but Pade approximation
            gammat_specified_Matrix[i][j] = gamman_specified_Matrix[i][j]
    return kn_specified_Matrix,kt_specified_Matrix,gamman_specified_Matrix,gammat_specified_Matrix