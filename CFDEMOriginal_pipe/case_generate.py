import sys
import os
import json
import pathlib
import shutil
import subprocess
import glob
import math
import numpy as np

sys.path.append(os.path.abspath("../"))
import loadJSONPara
import helixWallTextureGenerator
import defaultWallTextureGenerator
import mathUtility
import fileUtility
import liggghtsPipeUtility as lpu

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
        pipeSpecies = lpu.cvrtSpeciesID_LiggghtsToPython(pipeSpecies)
        textureSpecies = loadJSONPara.read(json_file_path,"textureSpecies") 
        textureSpecies = lpu.cvrtSpeciesID_LiggghtsToPython(textureSpecies)
        particleSpecies = loadJSONPara.read(json_file_path,"particleSpecies") 
        particleSpecies = lpu.cvrtSpeciesID_LiggghtsToPython(particleSpecies)

        particleDensity = loadJSONPara.read(json_file_path,"particleDensity") 
        particleRadius = loadJSONPara.read(json_file_path,"particleRadius") 
        particleSolidFrac = loadJSONPara.read(json_file_path,"particleSolidFrac") 

        particleInsertionDiskSize = loadJSONPara.read(json_file_path,"particleInsertionDiskSize") 
        particleInsertionDiskVolumeFraction = loadJSONPara.read(json_file_path,"particleInsertionDiskVolumeFraction") 

        enableHelixWallTexture = loadJSONPara.readwithdefault(json_file_path,"enableHelixWallTexture",False) 
        if enableHelixWallTexture == True:
            helixAmplitude = loadJSONPara.read(json_file_path,"helixAmplitude") 
            helixPeriod = loadJSONPara.read(json_file_path,"helixPeriod") 

        fluidVelocity = loadJSONPara.read(json_file_path,"fluidVelocity")
        kinematicViscosity = loadJSONPara.read(json_file_path,"kinematicViscosity")
        gravity = loadJSONPara.read(json_file_path,"gravity")
        fluidGridSizeXY = loadJSONPara.read(json_file_path,"fluidGridSizeXY")
        fluidGridSizeZ = loadJSONPara.read(json_file_path,"fluidGridSizeZ")

        numOfProcessor = loadJSONPara.readwithdefault(json_file_path,"numOfProcessor",4) 
        ###########################################################################
        # Create random prime number base
        ######################################################################################
        numOfSeedsRequiredPerCase = 4
        primeNumberDataBase = mathUtility.generate_prime_list(numOfTrial * numOfSeedsRequiredPerCase, 31, 500000)
        print("primeNumberDataBase = {}".format(primeNumberDataBase))
        
        ###########################################################################
        # Create folder
        #####################################################################################
        for trialid in range(numOfTrial):
            output_directory = pathlib.Path(json_file_path).stem + "_{}".format(trialid)
            if not os.path.exists(output_directory):
                os.mkdir(output_directory)
            
            ###########################################################################
            # Create CFD
            #####################################################################################
            CFD_dir = output_directory + "/CFD/"
            if not os.path.exists(CFD_dir):
                os.mkdir(CFD_dir)
            ###########################################################################
            # Create CFD/0
            #####################################################################################
            CFD_0_dir = CFD_dir + "/0/"
            if not os.path.exists(CFD_0_dir):
                os.mkdir(CFD_0_dir)

            variables = {
            "fluidVelocity": fluidVelocity
            }
            
            fileUtility.copyfile("CFD_template/0_template/U", CFD_0_dir + "U", variables)
            fileUtility.copyfile("CFD_template/0_template/Us", CFD_0_dir + "Us")
            fileUtility.copyfile("CFD_template/0_template/p", CFD_0_dir + "p")
            fileUtility.copyfile("CFD_template/0_template/voidfraction", CFD_0_dir + "voidfraction")
            fileUtility.copyfile("CFD_template/0_template/Ksl", CFD_0_dir + "Ksl")
            fileUtility.copyfile("CFD_template/0_template/rho", CFD_0_dir + "rho")
            ###########################################################################
            # Create CFD/system
            #####################################################################################
            CFD_system_dir = CFD_dir + "/system/"
            if not os.path.exists(CFD_system_dir):
                os.mkdir(CFD_system_dir)

            gridnumXY = math.ceil(2*pipeRadius/fluidGridSizeXY)
            print("gridnumXY = {}".format(gridnumXY))
            gridnumXYPart1 = math.ceil(gridnumXY/2)
            gridnumXYPart2 = math.ceil(gridnumXYPart1/2)
            gridnumZ = math.ceil(pipeLength/fluidGridSizeZ)
            print("gridnumZ = {}".format(gridnumZ))
             
            variables = {
            "pipeRadius": pipeRadius,
            "pipeLength": pipeLength,
            "halfpipeRadius": pipeRadius * 0.5,
            "diagonalPos": pipeRadius * np.cos(np.pi/4),
            "gridnumXYPart1": gridnumXYPart1,
            "gridnumXYPart2": gridnumXYPart2,
            "gridnumZ": gridnumZ
            }
            fileUtility.copyfile("CFD_template/system_template/blockMeshDict", CFD_system_dir + "blockMeshDict", variables)
            fileUtility.copyfile("CFD_template/system_template/fvSchemes", CFD_system_dir + "fvSchemes")
            fileUtility.copyfile("CFD_template/system_template/fvSolution", CFD_system_dir + "fvSolution")
            variables = {
            "numOfProcessor":numOfProcessor
            }
            fileUtility.copyfile("CFD_template/system_template/decomposeParDict", CFD_system_dir + "decomposeParDict",variables)
            variables = {
            "simulationTime":simulationTime,
            "timeStepSize":timeStepSize,
            "dumpTimeInteval":dumpTimeInteval
            }
            fileUtility.copyfile("CFD_template/system_template/controlDict", CFD_system_dir + "controlDict",variables)

            ###########################################################################
            # Create CFD/constant
            #####################################################################################
            CFD_constant_dir = CFD_dir + "/constant/"
            if not os.path.exists(CFD_constant_dir):
                os.mkdir(CFD_constant_dir)
            fileUtility.copyfile("CFD_template/constant_template/turbulenceProperties", CFD_constant_dir + "turbulenceProperties")
            variables = {
            "kinematicViscosity":kinematicViscosity
            }
            fileUtility.copyfile("CFD_template/constant_template/transportProperties", CFD_constant_dir + "transportProperties",variables)

            variables = {
            "gravity":gravity
            }
            fileUtility.copyfile("CFD_template/constant_template/g", CFD_constant_dir + "g",variables)
            fileUtility.copyfile("CFD_template/constant_template/couplingProperties", CFD_constant_dir + "couplingProperties")
            fileUtility.copyfile("CFD_template/constant_template/liggghtsCommands", CFD_constant_dir + "liggghtsCommands")
            
            ###########################################################################
            # Create DEM 
            #####################################################################################
            DEM_dir = output_directory + "/DEM/"
            if not os.path.exists(DEM_dir):
                os.mkdir(DEM_dir)
            
            ###########################################################################
            # Create wall texture
            ######################################################################################
            simulationBoxExtendFactor = 1.2
            simulationBoxEpsilon = 2e-15
            textureDir = output_directory + "/DEM/"
            if enableHelixWallTexture == True:
                textureRadius = particleRadius[textureSpecies]
                textureDensity = particleDensity[textureSpecies]
                textureAtomtype = lpu.cvrtSpeciesID_PythonToLiggghts(textureSpecies)
                helixWallTextureGenerator.helixWallTextureGenerator(
                                nparticletype, 
                                helixAmplitude,
                                helixPeriod,
                                pipeLength,
                                textureAtomtype,
                                textureRadius,
                                textureDensity,
                                pipeRadius * simulationBoxExtendFactor,
                                pipeRadius * simulationBoxExtendFactor,
                                -simulationBoxEpsilon,
                                pipeLength,
                                0,
                                textureDir)      
            else:
                dummyAtomtype = lpu.cvrtSpeciesID_PythonToLiggghts(particleSpecies)
                dummyRadius = particleRadius[particleSpecies]
                dummyDensity = particleDensity[particleSpecies]
                defaultWallTextureGenerator.defaultWallTextureGenerator( 
                            nparticletype, 
                            pipeLength,
                            dummyAtomtype,
                            dummyRadius,
                            dummyDensity,
                            pipeRadius * simulationBoxExtendFactor,
                            pipeRadius * simulationBoxExtendFactor,
                            -simulationBoxEpsilon,
                            pipeLength,
                            textureDir)


            ###########################################################################
            # Create DEM/init.in
            #####################################################################################
            target_content = []
            target_content.append( "################### variable definition start #####################\n")
            target_content.append("variable xmin equal {}\n".format(-pipeRadius * simulationBoxExtendFactor))
            target_content.append("variable xmax equal {}\n".format(pipeRadius * simulationBoxExtendFactor))
            target_content.append("variable ymin equal {}\n".format(-pipeRadius * simulationBoxExtendFactor))
            target_content.append("variable ymax equal {}\n".format(pipeRadius * simulationBoxExtendFactor))
            target_content.append("variable zmin equal {}\n".format(-simulationBoxEpsilon))
            target_content.append("variable zmax equal {}\n".format(pipeLength+simulationBoxEpsilon))
            target_content.append("variable pipelength equal {}\n".format(pipeLength))
            target_content.append("variable piperadius equal {}\n".format(pipeRadius))
               
            target_content.append( "################### variable definition finished #####################\n\n\n")
            fileUtility.copyfilePrepend("DEM_template/init_template.in" ,DEM_dir + "/init.in",target_content)

            ###########################################################################
            # Create DEM/run.in
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
            kn_specified_Matrix,kt_specified_Matrix,gamman_specified_Matrix,gammat_specified_Matrix = \
                lpu.computeForceCoefficient(nparticletype, pipeSpecies, textureSpecies, youngModulus,poissionRatio, coefficientOfRestitution, DensityArray, RadiusArray)
            target_content.append( "##### particle properties kn_specified #######\n")
            for i in range(nparticletype):
                for j in range(nparticletype): 
                    target_content.append("variable kn{}{} equal {}\n".format(lpu.cvrtSpeciesID_PythonToLiggghts(i),lpu.cvrtSpeciesID_PythonToLiggghts(j),kn_specified_Matrix[i][j]))
            target_content.append( "##### particle properties kt_specified #######\n")
            for i in range(nparticletype):
                for j in range(nparticletype): 
                    target_content.append("variable kt{}{} equal {}\n".format(lpu.cvrtSpeciesID_PythonToLiggghts(i),lpu.cvrtSpeciesID_PythonToLiggghts(j),kt_specified_Matrix[i][j]))
            target_content.append( "##### particle properties gamman_specified_Matrix #######\n")
            for i in range(nparticletype):
                for j in range(nparticletype): 
                    target_content.append("variable gamman{}{} equal {}\n".format(lpu.cvrtSpeciesID_PythonToLiggghts(i),lpu.cvrtSpeciesID_PythonToLiggghts(j),gamman_specified_Matrix[i][j]))
            target_content.append( "##### particle properties gammat_specified_Matrix #######\n")
            for i in range(nparticletype):
                for j in range(nparticletype): 
                    target_content.append("variable gammat{}{} equal {}\n".format(lpu.cvrtSpeciesID_PythonToLiggghts(i),lpu.cvrtSpeciesID_PythonToLiggghts(j),gammat_specified_Matrix[i][j]))
            target_content.append( "##### particle properties sliding friction #######\n")
            for i in range(nparticletype):
                for j in range(nparticletype): 
                    if i == pipeSpecies or j == pipeSpecies or i == textureSpecies or j == textureSpecies:
                        target_content.append("variable sfc{}{} equal {}\n".format(lpu.cvrtSpeciesID_PythonToLiggghts(i),lpu.cvrtSpeciesID_PythonToLiggghts(j),particleWallfrictionCoeffient))
                    else:
                        target_content.append("variable sfc{}{} equal {}\n".format(lpu.cvrtSpeciesID_PythonToLiggghts(i),lpu.cvrtSpeciesID_PythonToLiggghts(j),particleFrictionCoeffient))

            target_content.append( "##### particle properties geo #######\n")
            for i in range(nparticletype):
                if i >= particleSpecies:
                    target_content.append("variable particleradii{} equal {}\n".format(lpu.cvrtSpeciesID_PythonToLiggghts(i),particleRadius[i]))
                    target_content.append("variable particlefraction{} equal {}\n".format(lpu.cvrtSpeciesID_PythonToLiggghts(i),particleSolidFrac[i]/sum(particleSolidFrac)))
                    target_content.append("variable particledensity{} equal {}\n".format(lpu.cvrtSpeciesID_PythonToLiggghts(i),particleDensity[i]))
            
            target_content.append( "##### particle insertion #######\n")
            target_content.append("variable particleInsertionDiskSize equal {}\n".format(particleInsertionDiskSize))
            target_content.append("variable particleInsertionDiskVolumeFraction equal {}\n".format(particleInsertionDiskVolumeFraction))

            target_content.append( "##### time #######\n")
            target_content.append("variable dt equal {}\n".format(timeStepSize))
            target_content.append("variable dumpstep equal {}\n".format(math.ceil(dumpTimeInteval/timeStepSize)))
            target_content.append("variable simulationstep equal {}\n".format(math.ceil(simulationTime/timeStepSize)))

            target_content.append( "##### time #######\n")
            target_content.append("variable definedGravity equal {}\n".format(gravity))
            
            target_content.append( "##### seed #######\n")
            for i in range(numOfSeedsRequiredPerCase):
                target_content.append("variable theseed{} equal {}\n".format(i,primeNumberDataBase[trialid * numOfSeedsRequiredPerCase + i]))

            target_content.append( "variable numOfProcessor equal {}\n".format(numOfProcessor))

            target_content.append( "################### variable definition finished #####################\n\n\n")
            fileUtility.copyfilePrepend("DEM_template/run_template.in" ,DEM_dir + "/run.in",target_content)

            ###########################################################################
            # Create initNrun.sh
            #####################################################################################
            variables = {
            "numOfProcessor": numOfProcessor
            }
            fileUtility.copyfile("initNrun_template.sh", output_directory + "/" + "initNrun.sh", variables)