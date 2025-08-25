import sys
import os
import json
import pathlib
import shutil
import subprocess
import glob
import random
import math
import numpy as np
import loadJSONPara
import fileUtility
import getA

#####################################################################################
############# utility function                                #############
######################################################################################
# prime number, adapted from Chatgpt
def is_prime(n):
    """Checks if a number is prime using trial division."""
    if n < 2:
        return False
    if n in (2, 3):
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True

def random_large_prime(bit_length):
    """Generates a random prime number of the given bit length."""
    while True:
        num = random.getrandbits(bit_length) | 1  # Ensure it's odd
        if is_prime(num):
            return num

def generate_prime_list(count, bit_length, minValue):
    """Generates a list of 'count' random large prime numbers."""
    primes = []
    while len(primes) < count:
        prime = random_large_prime(bit_length)
        if prime > minValue:
            if prime not in primes:  # Ensure uniqueness
                primes.append(prime)
    return primes

def computeSphereVolume(sphereRadius):
    return (4/3) * np.pi * (sphereRadius **3 )


def computeForceCoefficient(YoungsModulus, PoissonRatio, CoefficientOfRestitution, DensityArray, RadiusArray):
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
                Meff = DensityArray[j]*computeSphereVolume(RadiusArray[j])  
            elif(j==pipeSpecies):
                Reff = RadiusArray[i]
                Meff = DensityArray[i]*computeSphereVolume(RadiusArray[i])
            elif(i==textureSpecies): # infinite mass, finite radius
                Reff = 1.0/(1.0/RadiusArray[i] + 1.0/RadiusArray[j])
                Meff = DensityArray[j]*computeSphereVolume(RadiusArray[j])
            elif(j==textureSpecies):
                Reff = 1.0/(1.0/RadiusArray[i] + 1.0/RadiusArray[j])
                Meff = DensityArray[i]*computeSphereVolume(RadiusArray[i])
            else:
                m = np.zeros((2,1))
                m[0] = DensityArray[i]*computeSphereVolume(RadiusArray[i])
                m[1] = DensityArray[j]*computeSphereVolume(RadiusArray[j])
                Reff = 1.0/(1.0/RadiusArray[i] + 1.0/RadiusArray[j])
                Meff = 1.0/(1.0/m[0] + 1.0/m[1])
            
            # Need verification for the computation of A_n
            density = DensityArray[0]
            CharacteristicVelocity = 0.1
            An =  getA.getAPhysical(CoefficientOfRestitution,Reff,Reff,density,YoungsModulus,PoissonRatio,CharacteristicVelocity)
            # https://www.cfdem.com/media/DEM/docu/gran_model_hertz_stiffness.html
            kn_specified_Matrix[i][j] = (4/3)* Yeff 
            kt_specified_Matrix[i][j] = 8.0*Geff
            gamman_specified_Matrix[i][j] = 2.0*Yeff*An/Meff # not liggghts definition, but Pade approximation
            gammat_specified_Matrix[i][j] = gamman_specified_Matrix[i][j]
    return kn_specified_Matrix,kt_specified_Matrix,gamman_specified_Matrix,gammat_specified_Matrix

if __name__ == "__main__":
    ###########################################################################
    # Read parameters                             
    ######################################################################################
    # read from json for input parameters
    try:
        json_input_file_path = sys.argv[1]
        print("load json path = {} ".format(json_input_file_path))
    except:
        raise Exception("no json path specify as the first argument")

    for json_file_path in glob.glob(json_input_file_path ):
        print("json path = {} ".format(json_file_path))
        numOfTrial = loadJSONPara.readwithdefault(json_file_path,"numOfTrial",1)
        timeStepSize = loadJSONPara.read(json_file_path,"timeStepSize")
        dumpTimeInteval = loadJSONPara.read(json_file_path,"dumpTimeInteval") 
        insertionTime = loadJSONPara.read(json_file_path,"insertionTime") 
        unrecordedSimulationTime = loadJSONPara.readwithdefault(json_file_path,"unrecordedSimulationTime", 0)
        simulationTime = loadJSONPara.read(json_file_path,"simulationTime")
        pipeRadius = loadJSONPara.read(json_file_path,"pipeRadius") 
        pipeLength = loadJSONPara.read(json_file_path,"pipeLength") 
        particleFrictionCoeffient = loadJSONPara.read(json_file_path,"particleFrictionCoeffient") 
        particleWallfrictionCoeffient = loadJSONPara.read(json_file_path,"particleWallfrictionCoeffient") 
        youngModulus = loadJSONPara.read(json_file_path,"youngModulus") 
        poissionRatio = loadJSONPara.read(json_file_path,"poissionRatio") 
        coefficientOfRestitution = loadJSONPara.read(json_file_path,"coefficientOfRestitution") 
        nparticletype = loadJSONPara.read(json_file_path,"nparticletype") 
        pipeSpecies = loadJSONPara.read(json_file_path,"pipeSpecies") 
        textureSpecies = loadJSONPara.read(json_file_path,"textureSpecies") 
        particleSpecies = loadJSONPara.read(json_file_path,"particleSpecies") 

        particleDensity = loadJSONPara.read(json_file_path,"particleDensity") 
        particleRadius = loadJSONPara.read(json_file_path,"particleRadius") 
        particleSolidFrac = loadJSONPara.read(json_file_path,"particleSolidFrac") 
        textureParticleRadius = loadJSONPara.read(json_file_path,"textureParticleRadius") 
        ###########################################################################
        # Create random prime number base
        ######################################################################################
        numOfSeedsRequiredPerCase = 4
        primeNumberDataBase = generate_prime_list(numOfTrial * numOfSeedsRequiredPerCase, 31, 500000)
        print("primeNumberDataBase = {}".format(primeNumberDataBase))

        ###########################################################################
        # Create folder
        #####################################################################################
        for trialid in range(numOfTrial):
            output_directory = pathlib.Path(json_file_path).stem + "_{}".format(trialid)
            if not os.path.exists(output_directory):
                os.mkdir(output_directory)

            ###########################################################################
            # Create init.in
            #####################################################################################
            target_content = []
            target_content.append( "################### variable definition start #####################\n")
            simulationBoxExtendFactor = 1.2
            simulationBoxEpsilon = 2e-15
            target_content.append("variable xmin equal {}\n".format(-pipeRadius * simulationBoxExtendFactor))
            target_content.append("variable xmax equal {}\n".format(pipeRadius * simulationBoxExtendFactor))
            target_content.append("variable ymin equal {}\n".format(-pipeRadius * simulationBoxExtendFactor))
            target_content.append("variable ymax equal {}\n".format(pipeRadius * simulationBoxExtendFactor))
            target_content.append("variable zmin equal {}\n".format(-simulationBoxEpsilon))
            target_content.append("variable zmax equal {}\n".format(pipeLength+simulationBoxEpsilon))
            target_content.append("variable pipelength equal {}\n".format(pipeLength))
            target_content.append("variable piperadius equal {}\n".format(pipeRadius))
            target_content.append( "################### variable definition finished #####################\n\n\n")
          
            target_path = output_directory + "/init.in"
            with open(target_path, "w") as target_file:
                with open("init_template.in") as templ_file:
                    templ_content = templ_file.readlines()
                    target_content = target_content + templ_content
                    target_file.writelines(target_content)
                    print("write file: {}".format(target_path))
            ###########################################################################
            # Create run.in
            #####################################################################################
            target_content = []
            target_content.append( "################### variable definition start #####################\n")
            target_content.append( "##### Geometry definition #######\n")
            target_content.append("variable pipelength equal {}\n".format(pipeLength))
            target_content.append("variable piperadius equal {}\n".format(pipeRadius))

            target_content.append( "##### particle definition #######\n")
            target_content.append( "variable nparticletype equal {} #{} -> Pipe, #{} -> Helix, #{}..#n ->Particles\n".format(nparticletype,pipeSpecies,textureSpecies,particleSpecies))
   
            DensityArray = particleDensity
            RadiusArray = particleRadius
            kn_specified_Matrix,kt_specified_Matrix,gamman_specified_Matrix,gammat_specified_Matrix = computeForceCoefficient(youngModulus,poissionRatio, coefficientOfRestitution, DensityArray, RadiusArray)
            target_content.append( "##### particle properties kn_specified #######\n")
            for i in range(nparticletype):
                for j in range(nparticletype): 
                    target_content.append("variable kn{}{} equal {}\n".format(i,j,kn_specified_Matrix[i][j]))
            target_content.append( "##### particle properties kt_specified #######\n")
            for i in range(nparticletype):
                for j in range(nparticletype): 
                    target_content.append("variable kt{}{} equal {}\n".format(i,j,kt_specified_Matrix[i][j]))
            target_content.append( "##### particle properties gamman_specified_Matrix #######\n")
            for i in range(nparticletype):
                for j in range(nparticletype): 
                    target_content.append("variable gamman{}{} equal {}\n".format(i,j,gamman_specified_Matrix[i][j]))
            target_content.append( "##### particle properties gammat_specified_Matrix #######\n")
            for i in range(nparticletype):
                for j in range(nparticletype): 
                    target_content.append("variable gammat{}{} equal {}\n".format(i,j,gammat_specified_Matrix[i][j]))
            target_content.append( "##### particle properties sliding friction #######\n")
            for i in range(nparticletype):
                for j in range(nparticletype): 
                    if i == pipeSpecies or j == pipeSpecies or i == textureSpecies or j == textureSpecies:
                        target_content.append("variable sfc{}{} equal {}\n".format(i,j,particleWallfrictionCoeffient))
                    else:
                        target_content.append("variable sfc{}{} equal {}\n".format(i,j,particleFrictionCoeffient))

            target_content.append( "##### particle properties geo #######\n")
            for i in range(nparticletype):
                if i >= particleSpecies:
                    target_content.append("variable particleradii{}{}{} equal {}\n".format(i,i,i,particleRadius[i]))
                    target_content.append("variable particlefraction{}{}{} equal {}\n".format(i,i,i,particleSolidFrac[i]/sum(particleSolidFrac)))
                    target_content.append("variable particledensity{}{} equal {}\n".format(i,i,particleDensity[i]))
            
            target_content.append( "##### particle number #######\n")
            particleNum = 0
            pipeVolume = np.pi * (pipeRadius **2 ) * pipeLength
            for radius,solidFrac in zip(particleRadius,particleSolidFrac):
                particleVolume = computeSphereVolume(radius)
                particleNum += solidFrac * pipeVolume/particleVolume
            particleNum = math.ceil(particleNum)
            print("Total particleNum : {}".format(particleNum))
            target_content.append("variable realparticleinsert equal {}\n".format(particleNum))
            target_content.append("variable maxparticleinsert equal {}\n".format(particleNum*2))
            target_content.append("variable insertionparticlerate equal {}\n".format(5 * math.ceil(particleNum/insertionTime))) #faster insertion

            target_content.append( "##### time #######\n")
            target_content.append("variable dt equal {}\n".format(timeStepSize))
            target_content.append("variable dumpstep equal {}\n".format(math.ceil(dumpTimeInteval/timeStepSize)))
            target_content.append("variable insertionstep equal {}\n".format(math.ceil(insertionTime/timeStepSize)))
            target_content.append("variable unrecordedsimulationstep equal {}\n".format(math.ceil(unrecordedSimulationTime/timeStepSize)))
            target_content.append("variable simulationstep equal {}\n".format(math.ceil(simulationTime/timeStepSize)))

            target_content.append( "##### seed #######\n")
            for i in range(numOfSeedsRequiredPerCase):
                target_content.append("variable theseed{} equal {}\n".format(i,primeNumberDataBase[trialid * numOfSeedsRequiredPerCase + i]))
        
            target_content.append( "################### variable definition finished #####################\n\n\n")
            target_path = output_directory + "/run.in"
            
            with open(target_path, "w") as target_file:
                with open("run_template.in") as templ_file:
                    templ_content = templ_file.readlines()
                    target_content = target_content + templ_content
                    target_file.writelines(target_content)
                    print("write file: {}".format(target_path))
            ###########################################################################
            # Create sbatch
            #####################################################################################
            sbatchName = output_directory
            sbatch_variables = {
                "dummyvariable": 100
            }
            target_path = output_directory + "/" + sbatchName + ".sbatch"
            with open(target_path, "w") as target_file:
                with open("template.sbatch") as templ_file:
                    target_content = templ_file.read()
                    for key, val in sbatch_variables.items():
                        target_content = target_content.replace(f"{{{key}}}", str(val))
                    target_file.write(target_content)
                    print("write file: {}".format(target_path))